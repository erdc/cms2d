!=====================================================================    
module met_def
! Meteorological (Air-Sea Interaction) Variable Definitions
! Includes wind, atmospheric pressure, evaporation, and precipitation
! written by Alex Sanchez, USACE-CHL
!=====================================================================       
    use geo_def, only: projection
    use prec_def
    implicit none
    save

    !Files
    character(len=200) :: windpath,windfile,prespath,presfile,windlocfile
    character(len=200) :: wndfl_intpcoef_file

    !Parameters
    integer     :: iwndlagr    !Wind reference frame. 0-Eulerian, 1-Lagrangian
    real(ikind) :: densitair   !Air density [kg/m^3]
    real(ikind) :: rhorat      !Air density divided by water density [-]
    real(ikind) :: cdragwind   !Wind drag coefficient [-]
    real(ikind) :: wndht       !Anemometer height [m]
    real(ikind) :: wndfac      !Scaling factor used to convert from different time averages (e.g. 0.81 from 1 min to 30 min)
    real(ikind) :: wkappa      !Wind drag coefficient scaling factor [-]
    real(ikind) :: ambatmpres  !Ambient atmpheric pressure [Pa]

    !Wind speed and shear stresses
    real(ikind), allocatable :: uwind(:),vwind(:)                   !Wind speed at current time step [m/s]
    real(ikind), allocatable :: uwind1(:),vwind1(:)                 !Wind speed at previous time step [m/s]
    real(ikind), allocatable :: uwind2(:),vwind2(:)                 !Wind speed two time steps ago [m/s]
    real(ikind), allocatable :: tauwindx(:),tauwindy(:)             !Wind stress [N/m^2]
    real(ikind) :: wndx,wndy,tauwx,tauwy

    !Lagrangian reference frame
    real(ikind), allocatable :: cdWndareap(:)            !Wind drag coef * wnd speed * areap, used for Lagrangian ref frame      
    
    !Meteorological Stations
    logical :: windsta,pressta
    integer :: nMetSta
    type metstadriver
      character(len=100) :: name        !Station Name
      character(len=200) :: file        !Station file name *.txt (also used as wind speed file)
      character(len=200) :: path        !Station file name *.txt (also used as wind direction file)
      integer       :: ntimes      !Number of times
      integer       :: inc         !Counter for temporal interpolation
      real(ikind)   :: wndfac      !Temporal factor for each station [-] (.e.g. 2-min avg to 30-min avg)
      real(ikind)   :: wndht       !Anemometer height [m]
      real(ikind)   :: xsta,ysta   !Coordinates in flow coordinate system, [m]
      real(ikind)   :: xlon,ylat   !Coordinates in geographic coordinates, [deg]
      real(ikind),pointer :: times(:)    !Times in hours from simulation starting time, [hr]
      real(ikind),pointer :: wndvalsx(:),wndvalsy(:) !Wind velocities, [m/s] (time)
      real(ikind),pointer :: presvals(:)       !Atmospheric pressure [Pa] (time)
      real(ikind)   :: wndx,wndy    !Wind speeds interpolated to current time step
      real(ikind)   :: pres         !Atmospheric pressure interpolated to current time step
      real(ikind)   :: pow          !Nearghest neighbor interpolation power
      integer       :: ntsi         !Spatial smoothing iterations
      integer       :: ntsw         !Spatial smoothing window width
    endtype metstadriver
    type(metstadriver), allocatable :: metsta(:)
    real(ikind), allocatable :: cntpmet(:,:)   !Interpolation coefficients

    !Wind curve
    logical :: windconst
    integer :: nwnd_inc,nwtimes         
    real(ikind), pointer :: wndvalsx(:),wndvalsy(:),wndtimes(:)    
    real(ikind), allocatable :: wndspeed(:),wnddirection(:)
    !type wndcurv_type
    !  integer :: inc
    !  integer :: ntimes
    !  real(ikind) :: times(:)  
    !  real(ikind) :: wndx(:)
    !  real(ikind) :: wndy(:)
    !endtype wndcurv_type
    !type(prescurv_type) :: wndcurv

    !Atmospheric Pressure curve
    logical :: presconst
    integer :: npres_inc,nptimes         
    real(ikind), pointer :: presvals(:),prestimes(:)
    !type prescurv_type
    !  integer :: inc
    !  integer :: ntimes
    !  real(ikind) :: times(:)  
    !  real(ikind) :: pres(:)      
    !endtype prescurv_type
    !type(prescurv_type) :: prescurv

    !Spatially variable wind and pressure fields, Alex, Mar 15, 2010
    logical :: windvar,wndpresloc
    logical :: write_windpresgrid !Write wind/pressure grid
    integer :: wunit,cunit
    !Input wind format: windformat = 
    !  0-Oceanweather
    !  1-Blended Sea Winds
    !  2-HWIND
    !  3-Speed/Direction/Pressure Single File (ADCIRC NWS=3 format)
    !  6-WindX/WindY/Pressure Single File  (ADCIRC NWS=6 format)
    !  7-Stress/Pressure Single File
    integer :: windformat
    integer :: nwindi,nwindj
    real(ikind) :: windstarthour,windhr1,windhr2  !Referenced to hydro start time         
    real(ikind),allocatable :: cntpwnd(:,:,:)
    integer,    allocatable :: ijntpwnd(:,:)
    real(ikind),allocatable :: xwind(:,:),ywind(:,:)
    real(ikind),allocatable :: wndspdx2(:,:)
    real(ikind),allocatable :: wndspdy2(:,:)
    !ADCIRC NWS=6 format (single file format)
    integer :: nwlat,nwlon
    real(ikind) :: wlatmax,wlonmin,wlatinc,wloninc,wtiminc
    type(projection) :: projwnd
    real(ikind) :: undefvel,lowvel,highvel,defvel
    character(len=200) :: windpresgridfile
      
    !Atmospheric pressure
    logical :: presvar
    integer :: punit
    integer :: nsmoothwind  !Smoothing iterations for wind velocity field [-]
    integer :: nsmoothpres  !Smoothing iterations for pressure field [-], Atm Pressure assumed to vary at scales much larger than the grid
    real(ikind) :: presstarthour,preshr1,preshr2  !Referenced to hydro start time
    real(ikind),allocatable :: atmpres2(:,:)      !Atmospheric pressure on wind grid [Pa]
    real(ikind),allocatable :: pressatm(:),pressatm1(:),pressatm2(:) !Air pressure [Pa]
    real(ikind),allocatable :: pressatmdx(:),pressatmdy(:)  !Air pressure gradients [Pa]    
    real(ikind) :: undefpres,lowpres,highpres,defpres

    !Precipitation and evaporation
    logical :: rain_evap
    character(len=200) :: rf_filename,rainfile,evapfile
    integer :: rf_unit
    real(ikind) :: RF_FRAC_RO !Reduction factor for rain at dry cells [-]
    real(ikind) :: rain_time  ![s], Note: Used only by explicit time stepping
    real(ikind) :: rain,evap  !Implicit: Constant or termporally interpolated [m/s], Explicit: Expected in meters at hourly intervals [m/hr]
    
    !type precip_evap_type
    !  integer :: inc
    !  integer :: ntimes
    !  real(ikind), pointer :: times(:)
    !  real(ikind), pointer :: vals(:)      
    !  real(ikind) :: val
    !endtype precip_evap_type
    !type(precip_evap_type) :: precipevap
    
end module met_def