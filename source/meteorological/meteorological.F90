!===================================================================
! CMS Meteorological routines
!
! Contains the following:
!   General
!     met_default - Sets meteorological variable defaults
!     met_cards   - Reads meteorological cards and overrides defaults
!     met_init    - Initializes meteorological module variables
!     met_print   - Prints the meteorological settings to the 
!                   diagnostic file and the screen
!   Wind time-series curve
!     windcurve_block - Reads a wind time-series curve block
!     windcurve_init  - Initializes a wind time-series curve
!     windcurve_eval  - Interpolates a wind time-series curve
!   Atmospheric pressure time-series curve
!     prescurve_init - Initializes the atmospheric pressure curve
!     prescurve_eval - Interpolations the atmospheric pressure
!   Wind and Atmospheric pressure fields
!     windpresfield_block - Reads a wind and pressure field block 
!                            from the control file
!     windpresfield_eval - Evaluates the wind and pressure fields
!                           at the current time step
!   Oceanweather
!     oceanweather_block - Reads an Oceanweather input block
!   Meteorological stations
!     metsta_block - Reads a meteorological station block
!     metsta_init - Initializes the meteorological station variables
!   Forcing
!     wind_lagrangian_update - Updates the wind forcing variables
!                  when using a Lagriangian reference frames.
!
! written by Alex Sanchez, USACE-CHL
!===================================================================

!*************************************************************
    subroutine met_default
! Sets the default values to the wind variables    
! written by Alex Sanchez, USACE-CHL
!*************************************************************          
    use met_def
    use geo_lib, only: proj_default
    implicit none
    
    !File names and paths
    windfile = ''
    windpath = ''
    prespath = ''
    presfile = ''
    windlocfile = ''
    wndfl_intpcoef_file = ''
    
    !Parameters
    densitair = 1.2       !Air density [kg/m^3]
    rhorat = 0.00120      !Air density / water density [-]
    wkappa = 0.4          !Wind drag coefficient
    wndht = 10.0          !Anemometer height [m]
    ambatmpres = 101325.0 !Ambient atmospheric pressure [Pa] (Average Sea Level Pressure)
    iwndlagr = 0          !Lagrangian wind speed correction
    wndfac = 1.0          !Scaling factor used to convert from different time averages (e.g. 0.81 from 1 min to 30 min)
    
    !Projection
    call proj_default(projwnd)
    
    !Stations
    windsta = .false.  !Wind from weather stations    
    pressta = .false.  !Pressure from weather stations (not operational)
    nMetSta = 0
        
    !Wind Curve
    windconst = .false. !Time series curve of wind speed and direction (met conv)
    nwnd_inc = 1
    nwtimes = 0
    wndx = 0.0
    wndy = 0.0
    tauwx = 0.0
    tauwy = 0.0
    
    !Pressure Curve
    presconst = .false.
    npres_inc = 1
    nptimes = 0
    
    !Spatially variable wind and pressure
    wunit = 445     
    punit = 115
    cunit = 614
    nwindi = 0
    nwindj = 0    
    windvar = .false.   !Spatially variable wind fields
    presvar = .false.   !Spatially variable pressure fields
    wndpresloc = .false. !File with coordiantes for wind specified
    windhr1 = -1.0
    windhr2 = -1.0
    preshr1 = -1.0
    preshr2 = -1.0
    undefvel = -9999.0
    lowvel = -100.0  !Lower limit to wind velocity[m/s]
    highvel = 100.0  ![m/s]
    defvel = 0.0
    undefpres = -9999.0
    lowpres  = 100000.0 !Lower limit to Atmospheric pressure for Quality Control [Pa]. Lowest ever recorded is about 106500.0 Pa
    highpres = 110000.0 !Upper limit to Atmospheric pressure for Quality Control [Pa]. Highest ever recorded is about 108600.0 Pa
    defpres = undefpres
    nsmoothpres = 5      !Smoothing iterations for atmospheric pressure field [-]
    nsmoothwind = 0      !Smoothing iterations for wind velocity field [-]
    write_windpresgrid = .false.
    windpresgridfile = 'wind_grid.xyz'
    
    !Rain and Evaporation
    rain_evap = .false.
    rf_filename = ''
    rainfile = ''
    evapfile = ''
    rain = 0.0
    evap = 0.0
    
    return
    end subroutine met_default
    
!*************************************************************   
    subroutine met_cards(cardname,foundcard)
! Reads wind data from Model Parameters file
! written by Alex Sanchez, USACE-CHL
!*************************************************************     
    use comvarbl, only: mpfile,flowpath
    use met_def
    use out_def, only: outlist
    implicit none
    integer :: ierr
    character(len=37) :: cardname,cdum,aext
    logical :: foundcard
    
    foundcard = .true.
    select case(cardname)            
      case('WIND_OUT_TIMES_LIST')
        backspace(77)
        read(77,*) cardname, outlist(8)%ilist   
      
      case('WIND_REFERENCE_FRAME','WIND_RELATIVE_VELOCITY')
        backspace(77)
        read(77,*) cardname, cdum
        if(cdum(1:3)=='LAG')then !Lagrangian
          iwndlagr = 1
        else  !Eulerian
          iwndlagr = 0
        endif            
      
      !--- Parameters and Constants -----------------------------------------  
      case('WIND_DRAG_COEFFICIENT')
        backspace(77)
        read(77,*) cardname, wkappa
      
      case('WIND_DRAG_SCALE_FACTOR')
        backspace(77)
        read(77,*) cardname, wkappa    
        wkappa=wkappa*0.4
      
      case('WIND_FACTOR','WIND_VELOCITY_FACTOR')  
        backspace(77)
        read(77,*) cardname, wndfac    
        
      case('ANEMOMETER_HEIGHT')
        call card_scalar(77,'m','m',wndht,ierr) 
        
      case('AMBIENT_ATM_PRESSURE','AMBIENT_ATMOSPHERIC_PRESSURE',&
         'ATMOSPHERIC_PRESSURE','AMB_ATM_PRESSURE','AMB_ATM_PRES')
        call card_scalar(77,'mb','Pa',ambatmpres,ierr)  
      
      !--- Wind Curve -----------------------------------------  
      case('WIND_INPUT_CURVE')  
        call card_dataset(77,mpfile,flowpath,windfile,windpath,2)
        windconst = .true.    
      
      case('WIND_SPEED_CURVE')
        backspace(77)
        read(77,*) cardname,windfile
        call fileext(windfile,aext)
        if(aext(1:2)=='h5')then
          backspace(77)
          read(77,*) cardname,windfile,windpath 
        endif
        windconst = .true.
          
      case('WIND_DIRECTION_CURVE')
        backspace(77)
        read(77,*) cardname,windpath
        call fileext(windpath,aext)
        if(aext(1:2)=='h5')then
          backspace(77)
          read(77,*) cardname,windfile,windpath
        endif
        windconst = .true.
        
      case('WIND_CURVE_BEGIN')
        call windcurve_block  
        
      !--- Atmospheric Pressure Curve ---------------------------------  
      case('ATM_PRESSURE_INPUT_CURVE','ATM_PRESSURE_CURVE',&
               'PRESSURE_INPUT_CURVE',    'PRESSURE_CURVE')
        call card_dataset(77,mpfile,flowpath,presfile,prespath,1)
        presconst = .true.        
      
      !--- Spatially Variable Wind and Atmospheric Pressure ------------  
      case('OCEANWEATHER_WIND_FILE') 
        backspace(77)
        read(77,*) cardname,windfile
        windvar = .true.
        windformat = 0        
             
      case('OCEANWEATHER_PRES_FILE','OCEANWEATHER_PRESSURE_FILE')
        backspace(77)
        read(77,*) cardname,presfile
        presvar = .true.
        windformat = 0
          
      case('OCEANWEATHER_XY_FILE')
        backspace(77)
        read(77,*) cardname,windlocfile
        windformat = 0           
        wndpresloc = .true.
        
      case('OCEANWEATHER_BEGIN')
        call oceanweather_block
      
      case('WIND_FLEET_FILE')  !Speed and direction
        backspace(77)
        read(77,*) cardname,windfile      
        windvar = .true.
        presvar = .false.
        windformat = 3
        
      case('WIND_PRESSURE_SINGLE_FILE',&
        'WIND_PRESSURE_FLEET_FILE','WIND_AND_PRESSURE_FILE') !Wx,Wy,Pressure
        backspace(77)
        read(77,*) cardname,windfile      
        windvar = .true.
        presvar = .true.
        windformat = 6        
      
      case('STRESS_PRESSURE_SINGLE_FILE')
        backspace(77)
        read(77,*) cardname,windfile      
        windvar = .true.
        presvar = .true.
        windformat = 7

      case('WIND_PRESSURE_XY_FILE','WIND_PRES_XY_FILE',&
           'WIND_PRESSURE_COORDINATES_FILE') !Old
        backspace(77)
        read(77,*) cardname,windlocfile
        wndpresloc = .true.        
        
      case('WIND_PRESSURE_GRID_PARAM','WIND_GRID_PARAMETERS','STRESS_PRESSURE_GRID_PARAMETERS')
        backspace(77)
        read(77,*) cardname,nwindi,nwindj,wlatmax,wlonmin,wlatinc,wloninc     
!        NWLAT = number of latitude or y-coordinate values in the met file.
!        NWLON = number of longitude or x-coordinate values in met file.
!        WLATMAX = maximum latitude (decimal deg) or y-coordinate of data in met file (< 0 south of the equator).
!        WLONMIN = minimum longitude (decimal deg) or x-coordinateof data in the met file (< 0 west of Greenwich meridian).
!        WLATINC = latitude or y-coordinate increment (decimal deg) of data in the met file (must be > 0).
!        WLONINC = longitude or x-coordinate increment (decimal deg) of data in the met file (must be > 0).
    
      case('WIND_PRESSURE_TIME_INCREMENT','WIND_TIME_INCREMENT','STRESS_PRESSURE_TIME_INCREMENT')
        call card_scalar(77,'sec','hrs',wtiminc,ierr)
!        WTIMINC = meteorological wind time interval (hours)

      case('WIND_PRESSURE_SINGLE_BEGIN','WIND_PRESSURE_BEGIN')
        call windpresfield_block
        
      case('MET_STATION_BEGIN','METEOROLOGICAL_STATION_BEGIN','MET_STA_BEGIN')
          call metsta_block
         
      case('OUTPUT_WIND_PRESSURE_GRID')
        call card_boolean(77,write_windpresgrid,ierr)
          
      case('OUTPUT_WIND_PRESSURE_GRID_FILE')
        backspace(77)
        read(77,*) cardname,windpresgridfile 
        write_windpresgrid = .true.
          
      case('WIND_SPATIAL_SMOOTHING_ITERATIONS','WIND_SPATIAL_SMOOTH_ITER')    
        backspace(77)
        read(77,*) cardname,nsmoothwind 
        nsmoothwind = min(max(nsmoothwind,0),100)
        
      case('ATM_PRESSURE_SPATIAL_SMOOTHING_ITERATIONS','ATM_PRES_SPATIAL_SMOOTHING_ITERATIONS',&
           'ATM_PRES_SPATIAL_SMOOTH_ITER')
        backspace(77)
        read(77,*) cardname,nsmoothpres 
        nsmoothpres = min(max(nsmoothpres,0),100)  
          
      !--- Rain and Evaporation -----------------------------------------    
      case('RAIN_EVAP_INPUT_CURVE','PRECIP_EVAP_INPUT_CURVE',&
          'RAIN_EVAP_CURVE','PRECIP_EVAP_CURVE')
        backspace(77)
        read(77,*) cardname,rf_filename
        rain_evap = .true.
        
      case('RAIN_INPUT_CURVE','PRECIPITATION_INPUT_CURVE',&
          'RAIN_CURVE','PRECIPITATION_CURVE')
        backspace(77)
        read(77,*) cardname,rainfile
        rain_evap = .true.
        
      case('EVAP_INPUT_CURVE','EVAPORATION_INPUT_CURVE',&
          'EVAP_CURVE','EVAPORATION_CURVE')
        backspace(77)
        read(77,*) cardname,evapfile
        rain_evap = .true.
        
      case('RAIN_VALUE','RAIN','RAIN_CONSTANT','PRECIPITATION_CONSTANT','PRECIPITATION_VALUE')
        call card_scalar(77,'m/hr','m/hr',rain,ierr)
        rain_evap = .true.
        
      case('EVAPORATION_VALUE','EVAPORATION','EVAPORATION_CONSTANT')
        call card_scalar(77,'m/hr','m/hr',evap,ierr)
        rain_evap = .true.
        
      case default
        foundcard = .false.
          
    end select
    
    return        
    end subroutine met_cards

!************************************************************    
    subroutine windcurve_block
! Reads a wind time-series curve block from the CMS-Flow
! Model Control File
! written by Alex Sanchez, USACE-CHL
!************************************************************   
    use comvarbl, only: mpfile,flowpath
    use met_def
    use prec_def
    implicit none
    integer :: ii,ierr
    character(len=37) :: cardname,cdum
    logical :: foundcard
    
    windconst = .true.   
    
d1: do ii=1,15
      foundcard = .true.
      read(77,*,iostat=ierr) cardname
      if(ierr/=0) exit
      if(cardname(1:1)=='!' .or. cardname(1:1)=='#') cycle
      select case(cardname)
        case('WIND_FILE')  
          call card_dataset(77,mpfile,flowpath,windfile,windpath,1)
          if(len_trim(windfile)==0)then
            windfile = windpath
            windpath = ''
          endif
          
        case('STATION_NAME')  
          backspace(77)
          read(77,*) cardname, windfile
          
        case('STATION_PATH')  
          backspace(77)
          read(77,*) cardname, windpath  
        
        case('WIND_REFERENCE_FRAME','WIND_RELATIVE_VELOCITY')
          backspace(77)
          read(77,*) cardname, cdum
          if(cdum(1:3)=='LAG')then !Lagrangian
            iwndlagr = 1
          else  !Eulerian
            iwndlagr = 0
          endif            
        
        case('WIND_DRAG_COEFFICIENT')
          backspace(77)
          read(77,*) cardname, wkappa
      
        case('WIND_DRAG_SCALE_FACTOR')
          backspace(77)
          read(77,*) cardname, wkappa    
          wkappa=wkappa*0.4
      
        case('WIND_FACTOR','WIND_VELOCITY_FACTOR')  
          backspace(77)
          read(77,*) cardname, wndfac    
        
        case('ANEMOMETER_HEIGHT')
          call card_scalar(77,'m','m',wndht,ierr)  
        
        case('AMBIENT_ATM_PRESSURE','AMBIENT_ATMOSPHERIC_PRESSURE','ATMOSPHERIC_PRESSURE')
          call card_scalar(77,'mb','Pa',ambatmpres,ierr)  
        
        case('HORIZONTAL_PROJECTION_BEGIN','HORIZ_PROJ_BEGIN')
          call proj_horiz_block(77,projwnd)
              
        case('WIND_CURVE_END','END')
          exit d1
          
        case default
          foundcard = .false.
      end select
    enddo d1

    if(len_trim(windpath)==0)then
      windpath = flowpath    
    endif
    
    return        
    end subroutine windcurve_block
    
!************************************************************
    subroutine oceanweather_block
! Reads an Oceanweather block
! Author: Alex Sanchez, USACE-CHL
!************************************************************
    use geo_def, only: aHorizDatum,&
       aHorizCoordSystem,aHorizUnits
    use met_def
    implicit none
    integer :: ii,ierr
    character(len=37) :: cardname
    logical :: foundcard,projblock
    
    windformat = 0
    projblock = .false.
    
d1: do ii=1,4
      foundcard = .true.
      read(77,*,iostat=ierr) cardname
      if(ierr/=0) exit
      if(cardname(1:1)=='!' .or. cardname(1:1)=='#') cycle
      select case(cardname)
        case('WIND_FILE')
          backspace(77)
          read(77,*) cardname,windfile
          windvar = .true.
          call filepath(windfile,windpath)
          
        case('PRESSURE_FILE','PRES_FILE','PRESS_FILE')
          backspace(77)
          read(77,*) cardname,presfile
          presvar = .true.
          call filepath(presfile,prespath)
        
        case('HORIZONTAL_PROJECTION_BEGIN','HORIZ_PROJ_BEGIN')
          call proj_horiz_block(77,projwnd)
          projblock = .true.
              
        case('OCEANWEATHER_END','END')
          exit d1
          
        case default
          foundcard = .false.
      end select
    enddo d1

    if(.not.projblock)then !Set to geographic, degrees
      projwnd%iHorizDatum = 1         !Horizontal Datum = NAD83
      projwnd%iHorizCoordSystem = 0   !Horizontal Coordinate System = GEOGRAPHIC
      projwnd%iHorizZone = 0          !Horizontal Zone = 0
      projwnd%iHorizUnits = 4         !Horizontal Units = DEGREES
    endif
    
    return
    end subroutine oceanweather_block
    
!************************************************************
    subroutine windpresfield_block
! Reads a wind and pressure field block from the control file
! Author: Alex Sanchez, USACE-CHL
!************************************************************
    use geo_def, only: aHorizDatum,&
       aHorizCoordSystem,aHorizUnits
    use met_def
    implicit none
    integer :: ii,ierr
    character :: cardname*37
    logical :: foundcard
    
    windformat = 6
    windvar = .true.
    presvar = .true.
    
d1: do ii=1,4
      foundcard = .true.
      read(77,*,iostat=ierr) cardname
      if(ierr/=0) exit
      if(cardname(1:1)=='!' .or. cardname(1:1)=='#') cycle
      select case(cardname)
        case('WIND_PRESSURE_SINGLE_END','WIND_PRESSURE_END','END')
          exit d1
          
        case('FILE','WIND_PRESSURE_FILE')
          backspace(77)
          read(77,*) cardname,windfile               
        
        case('GRID_PARAM','GRID_PARAMETERS')
          backspace(77)
          read(77,*) cardname,nwindi,nwindj,wlatmax,wlonmin,wlatinc,wloninc
        
        case('TIME_INCREMENT')
          call card_scalar(77,'sec','hrs',wtiminc,ierr)
                
        case('XY_FILE','COORDINATES_FILE')
          backspace(77)
          read(77,*) cardname,windlocfile
          wndpresloc = .true.
        
        case('UNDEFINED_VELOCITY','VELOCITY_UNDEFINED')
          backspace(77)
          read(77,*) cardname,undefvel
        
        case('DEFAULT_VELOCITY','VELOCITY_DEFAULT')
          backspace(77)
          read(77,*) cardname,defvel
        
        case('UNDEFINED_PRESSURE','PRESSURE_UNDEFINED')
          backspace(77)
          read(77,*) cardname,undefpres
        
        case('DEFAULT_PRESSURE','PRESSURE_DEFAULT')
          backspace(77)
          read(77,*) cardname,defpres 
          
        case('OUTPUT_GRID_FILE')
          backspace(77)
          read(77,*) cardname,windpresgridfile 
          write_windpresgrid = .true.
        
        case('HORIZONTAL_PROJECTION_BEGIN','HORIZ_PROJ_BEGIN')
          call proj_horiz_block(77,projwnd)  
  
        case default
          foundcard = .false.
      end select
    enddo d1

    return
    end subroutine windpresfield_block
    
!*************************************************************   
    subroutine met_init()
! Initializes the meteorological variables
!
! written by Alex Sanchez, USACE-CHL    
!*************************************************************       
    use size_def
    use met_def
    use diag_def, only: debug_mode
    use out_def, only: write_wndvel,write_wndmag
    implicit none       
            
    !Wind speed and shear stresses
    if(windvar .or. windsta .or. windconst)then      
      allocate(uwind(ncellsD),vwind(ncellsD))
      uwind=0.0;  vwind=0.0
      if(iwndlagr==0)then
        allocate(tauwindx(ncellsD),tauwindy(ncellsD))        
        tauwindx=0.0; tauwindy=0.0
      endif
      if(windformat==7) iwndlagr = 0 !Eulerian reference frame
    endif
    
    !Lagrangian reference frame
    if(iwndlagr==1)then
      allocate(cdWndareap(ncellsD))
      cdWndareap=0.001
    endif        
    
    !Meteorological Stations
    if(windsta .or. pressta) call metsta_init
    
    !Spatially Constant Wind and Pressure
    if(windconst) call windcurve_init !Wind curve
    if(presconst) call prescurve_init !Pressure curve

    !Spatially Variable Wind and Pressure
    if(windvar .or. presvar)then
      if(debug_mode)then
        write_windpresgrid = .true.
      endif
      call windpres_init
    endif
    
    !!Rainfall and evaporation
    if(rain_evap)then 
      rf_unit=701
      open(rf_unit,file=rf_filename)
      read(rf_unit,*) !header
      read(rf_unit,*)rain,evap    
      rain_time = 0.0
      !close(rf_unit)
      !precipevap%val      
    endif    
    
    !Wind and Atmospheric Pressure
    if(windconst) call windcurve_eval    !Spatially constant wind field
    if(presconst) call prescurve_eval    !Spatially constant atmospheric pressure    
    if(windvar .or. presvar) call windpresfield_eval  !Spatially variable wind and pressure fields   
    if(windsta .or. pressta) call metsta_eval   !Spatially variable wind and pressure stations
    
    !Output
    if(windformat==7)then !Wind stress input
      write_wndvel = .false.
      write_wndmag = .false.
    endif
    
    return        
    end subroutine met_init

!*************************************************************   
    subroutine windcurve_init
! Reads wind data from Model Parameters file
! written by Alex Sanchez, USACE-CHL
!*************************************************************         
#include "CMS_cpp.h"
    use diag_lib
    use prec_def
    use met_def
    use met_lib,   only: read_stdmetfile,wind_heightcorr
    use geo_def, only: azimuth_fl
    use comvarbl, only: flowpath,tjulday0,tmax
    use const_def, only: deg2rad
    use comvarbl,  only: version
    use out_def,   only: write_ascii_input
    use in_lib,    only: read_tsd,read_xys
#ifdef XMDF_IO
    use xmdf
    use in_xmdf_lib, only: read_dataseth5
#endif

    implicit none
    integer :: i,ndat
    real(ikind) :: ang,fac,wndspd,wnddir,tjulend,tjuldaybegwnd
    real(ikind),pointer :: wdat(:,:)
    character :: aext*10,aname*100,atype*100
    character :: apath*100,awindcurve*100       !Added for new format since WindCurve isn't a folder holding other datasets.
    !logical :: foundfile
    
    !call fileparts(windfile,apath,aname,aext)
    !if(len_trim(apath)==0 .and. len_trim(flowpath)>0)then
    !  windfile = trim(flowpath) // windfile
    !endif       
    !inquire(file=windfile,exist=foundfile)
    !if(.not.foundfile)then
    !  write(*,*) 'ERROR: Oceanweather Pressure File: ',trim(windfile),' not found'
    !  write(*,*) '  Press any key to continue.'
    !  read(*,*)
    !  stop
    !endif    
    
    call fileext(windfile,aext)
    
    select case (aext)                                                         !Switched from IF/ELSE to CASE.  MEB 06/25/2021
    case ('h5')                    !Model parameters file
#ifdef XMDF_IO
      call read_dataseth5(windfile,windpath,'Times',nwtimes,wndtimes)
      call read_dataseth5(windfile,windpath,'Magnitude',nwtimes,wndvalsx)
      call read_dataseth5(windfile,windpath,'Direction',nwtimes,wndvalsy)
#endif

    case ('txt')                   !Standard Met File
      tjulend = tjulday0 + tmax/24.0
      windfile = adjustl(windfile)
      i = index(windfile,'h')
      aname = windfile(1:i-1)
      call read_stdmetfile(aname,windpath,tjulday0,tjulend,&
        nwtimes,wndtimes,wndvalsx,wndvalsy)
      
    case ('tsd')                   !Time Series Delimited file
      call read_tsd(windfile,aname,atype,ndat,nwtimes,tjuldaybegwnd,wndtimes,wdat)
      wndtimes = wndtimes/3600.0 + 24.0*(tjuldaybegwnd - tjulday0) !Note conversion to hours
      allocate(wndvalsx(nwtimes),wndvalsy(nwtimes))
      wndvalsx = wdat(:,1) !Speed [m/s]
      wndvalsy = wdat(:,2) !Direction [deg coming from clockwise from true North]
      deallocate(wdat)
      
    case ('xys')                   !XY Series formatted file
      call read_xys(windfile,nwtimes,wndtimes,wndvalsx)
      call read_xys(windpath,i,wndtimes,wndvalsy)
      if(i/=nwtimes)then
        call diag_print_error('Time series length in wind speed and direction ',&
          '  files do not match')  
      endif    
      
    case default
      call diag_print_error('Unknown wind file format')     !Added MEB 06/25/2021 to stop if the wind file format was unrecognized.
    end select 
    
    !Keep speeds and magnitudes if needed for writing ascii output later
    if(write_ascii_input)then
      allocate(wndspeed(nwtimes),wnddirection(nwtimes))
      wndspeed=wndvalsx
      wnddirection=wndvalsy
    endif
    
    fac=wndfac*wind_heightcorr(wndht) !Temporal and vertical coordinate factors
    ang=(180.0+azimuth_fl)*deg2rad     !Note rotation from meteological convention    
    do i=1,nwtimes      
      wndspd=wndvalsx(i)*fac           !Convert to 10-m reference height     
      wnddir=wndvalsy(i)*deg2rad+ang   !Convert to local coordinate system and in radians
      wndvalsx(i)=wndspd*sin(wnddir)   !Note sin and cos are switched because of coordinate system, 0=top
      wndvalsy(i)=wndspd*cos(wnddir) 
    enddo
    
    return
    end subroutine windcurve_init
    
!*************************************************************   
    subroutine prescurve_init()
! Reads wind data from Model Parameters file
! written by Alex Sanchez, USACE-CHL
!*************************************************************
    use comvarbl, only: flowpath,tjulday0
    use diag_lib
    use flow_def, only: grav
    use in_lib, only: read_tsd,read_xys
    use met_def
    use prec_def
    implicit none
    integer :: ndat
    real(ikind) :: tjuldaybegwnd
    real(ikind), pointer :: pdat(:,:)
    character(len=200) :: apath,aname
    character(len=10) :: aext
    character(len=100) :: atype
    logical :: foundfile
    
    call fileparts(presfile,apath,aname,aext)
    if(len_trim(apath)==0 .and. len_trim(flowpath)>0)then
      presfile = trim(flowpath) // presfile
    endif       
    inquire(file=presfile,exist=foundfile)
    if(.not.foundfile)then
      call diag_print_error('Could not find atmospheric pressure curve file: ',presfile)
    endif
    
    call fileext(presfile,aext)
    if(aext(1:3)=='tsd')then
      call read_tsd(presfile,aname,atype,ndat,nptimes,tjuldaybegwnd,prestimes,pdat)
      allocate(presvals(nptimes))
      presvals(1:nptimes) = pdat(1:nptimes,1)
      deallocate(pdat)
      prestimes = prestimes/3600.0 + 24.0*(tjuldaybegwnd - tjulday0) !Note conversion to hours
    elseif(aext(1:3)=='xys')then
      call read_xys(presfile,nptimes,prestimes,presvals)   
    else
      call diag_print_error('Invalid file format for pressure curve',presfile)        
    endif
    
    !Note: Since the input units cannot be specified, they should be assumed to be
    !in Pa to be consistent with other input parameters
    !!presvals=presvals*100.0 !convert from mbar to Pa ***** IMPORTANT *****
    
    return
    end subroutine prescurve_init
    
!*************************************************************    
    subroutine windpres_init
! Initializes the wind and atmospheric pressure field forcing    
! written by Alex Sanchez, USACE-CHL    
!*************************************************************        
    use size_def
    use comvarbl, only: flowpath
    use met_def
    use diag_lib
    implicit none
    character(len=200) :: apath,aname
    character(len=10) :: aext
    logical :: foundfile
    
    if(presvar .and. .not.windvar)then !Alex, Mar 15, 2010
      call diag_print_error('Cannot have variable pressure field without variable wind field')     
    endif
    
    !Interpolation variables for spatially variable wind and pressure fields
    if(len_trim(windpath)==0)then
      windpath = flowpath    
    endif
    wndfl_intpcoef_file =  trim(windpath) // 'Intpcoef_wndfl.bin'    
    allocate(ijntpwnd(2,ncellsD),cntpwnd(2,2,ncellsD))  
    ijntpwnd=0  !defualt set so non-assigned CMS cells will have zero
    cntpwnd=0.0    
    
    call windpres_grid !Wind/Pressure Grid
    call interp_coef_windpres2fl !Calculate and save spatial interpolation coefficients

    !Allocate and Initialize Variables    
    if(windvar)then
      allocate(wndspdx2(nwindi,nwindj),wndspdy2(nwindi,nwindj))
      wndspdx2=0.0;   wndspdy2=0.0
      allocate(uwind1(ncellsD),vwind1(ncellsD),uwind2(ncellsD),vwind2(ncellsD))
      uwind1 = 0.0; vwind1 = 0.0
      uwind2 = 0.0; vwind2 = 0.0
    endif
    if(presvar)then
      allocate(atmpres2(nwindi,nwindj)) 
      allocate(pressatm(ncellsD),pressatm1(ncellsD),pressatm2(ncellsD))
      pressatm = 0.0; pressatm1 = 0.0; pressatm2 = 0.0
      allocate(pressatmdx(ncellsD),pressatmdy(ncellsD))
      pressatmdx=0.0; pressatmdy=0.0
    endif    
    
    if(windformat==0 .or. windformat==3 .or. windformat==6 .or. windformat==7)then
      call fileparts(windfile,apath,aname,aext)
      if(len_trim(apath)==0 .and. len_trim(flowpath)>0)then
        windfile = trim(flowpath) // windfile
      endif       
      inquire(file=windfile,exist=foundfile)
      if(.not.foundfile)then
        call diag_print_error('Could not find wind/pressure file: ',windfile)
      endif
      if(windformat==0 .and. presvar)then
        call fileparts(presfile,apath,aname,aext)
        if(len_trim(apath)==0 .and. len_trim(flowpath)>0)then
          presfile = trim(flowpath) // presfile
        endif       
        inquire(file=presfile,exist=foundfile)
        if(.not.foundfile)then
          call diag_print_error('Could not find pressure file: ',presfile)
        endif    
      endif
    endif
    
    !Open Files
    if(windformat==0)then !OceanWeather, Inc.        
      open(wunit,file=windfile,status='old')
      read(wunit,*) !Skip first line 
      if(presvar)then
        open(punit,file=presfile,status='old')
        read(punit,*) !Skip first line 
      endif
    elseif(windformat==3 .or. windformat==6 .or. windformat==7)then !ADCIRC NWS=3,6,7 formats
      open(wunit,file=windfile,status='old')
    endif
    
    !Initialize spatially variable wind and pressure fields
    call windpresfield_eval
    
    return        
    end subroutine windpres_init

!**************************************************   
    subroutine met_print()
! Prints the meteorologic settings to the screen
! and diagnositic file
! written by Alex Sanchez, USACE-CHL
!**************************************************   
    use met_def
    use geo_def
    use diag_def, only: dgunit,dgfile
    use diag_lib
    implicit none
    integer :: i,iunit(2),ista
    character(len=200) :: apath
    character(len=100) :: aname,afile
    character(len=10) :: aext
    logical :: ok 

111 format(' ',A,T40,A)    
118 format(' ',A,T40,I0)    
119 format(' ',A,T40,I0,A)
262 format(' ',A,T40,I4)
193 format(' ',A,T40,F0.3,A)
102 format(' ',A,T40,F0.2,A)    
293 format(' ',A,T40,F0.3,2xF0.3,A)
    
    iunit = (/6, dgunit/)    
    open(dgunit,file=dgfile,access='append') 
    do i=1,2    
      write(iunit(i),*)    
      if(windconst .or. windvar .or. windsta)then
        write(iunit(i),111)     'Wind Forcing:','ON'
      else
        write(iunit(i),111)     'Wind Forcing:','OFF'  
      endif    
    
      !Wind
      if(windconst)then !Wind Curve
        call fileext(windfile,aext)
        if(aext(1:1)==' ')then
          write(iunit(i),111)   '  Met Station Name:',trim(windfile)
          if(len_trim(windpath)>0)then
            write(iunit(i),111) '  Met Station Path:',trim(windpath)
          endif
        else 
          call fileext(windfile,aext)
          if(aext(1:3)=='xys')then  
            write(iunit(i),111) '  Wind Speed:',trim(windfile)  
            write(iunit(i),111) '  Wind Direction:',trim(windpath)
          else
            write(iunit(i),111) '  Wind Curve:',trim(windpath)
          endif
        endif
        write(iunit(i),119)     '  Length of Wind Data:',nwtimes, ' records'          
      elseif(windvar)then !Wind Fields
        call fileparts(windfile,apath,aname,aext)
        afile = trim(aname) // '.' // trim(aext)
          write(iunit(i),111)     '  Wind Field:',trim(afile)
        write(iunit(i),111)     '    Horizontal Projection'
        write(iunit(i),111)     '      Coordinate System:',trim(aHorizCoordSystem(projwnd%iHorizCoordSystem))     
        if(projwnd%iHorizCoordSystem/=22)then
          write(iunit(i),111)   '      Datum:',trim(aHorizDatum(projwnd%iHorizDatum))
          if(projwnd%iHorizZone/=0)then
            write(iunit(i),262) '      Zone:',projwnd%iHorizZone
          endif
        endif      
        write(iunit(i),111)     '      Units:',trim(aHorizUnits(projwnd%iHorizUnits))
      endif     
    
      !Wind Parameters
      if(windconst .or. windvar)then
        write(iunit(i),193)     '  Anemometer Height:',wndht,' m'   
        write(iunit(i),193)     '  Drag Coefficient Factor:',wkappa/0.4
        if(iwndlagr==1)then
          write(iunit(i),111)   '  Reference frame:','LAGRANGIAN'    
        else
          write(iunit(i),111)   '  Reference frame:','EULERIAN'
        endif
      endif
      
      !Meteorological Stations
      if(windsta)then
        if(iwndlagr==1)then
          write(iunit(i),111)   '  Reference frame:','LAGRANGIAN'    
        else
          write(iunit(i),111)   '  Reference frame:','EULERIAN'
        endif   
        write(iunit(i),118)     '  Number of Met Stations:',nMetSta
        do ista=1,nMetSta  
          write(iunit(i),111)   '  Met Station:',trim(metsta(ista)%name)          
          write(iunit(i),193)   '    Anemometer Height:',metsta(ista)%wndht,' m'   
          write(iunit(i),193)   '    Wind Scaling Factor:',metsta(ista)%wndfac
          write(iunit(i),293)   '    Coordinates:',metsta(ista)%xsta,metsta(ista)%ysta,' m'
          write(iunit(i),193)   '    Power Parameter:',metsta(ista)%pow
          write(iunit(i),118)   '    Time Series Length:',metsta(ista)%ntimes
          if(len_trim(metsta(ista)%file)>0)then
            write(iunit(i),111) '    File:',trim(metsta(ista)%file)
            write(iunit(i),111) '    Path:',trim(metsta(ista)%path)
          endif
        enddo
      endif
      
      !Atmospheric Pressure Forcing
      if(presvar .and. windformat==0)then
        inquire(file=presfile,exist=ok)
        if(.not.ok)then
          call diag_print_error('Could not find pressure file: ',presfile)
        endif
        write(iunit(i),*)
        write(iunit(i),111)     'Atmospheric Pressure Forcing:','ON'
        call fileparts(presfile,apath,aname,aext)
        afile = trim(aname) // '.' // trim(aext)
          write(iunit(i),111)     '  Pressure Field:',trim(afile)
        write(iunit(i),102)     '  Ambient Atm Pressure:',ambatmpres,' Pa'
        write(iunit(i),111)     '    Horizontal Projection'
        write(iunit(i),111)     '      Coordinate System:',trim(aHorizCoordSystem(projwnd%iHorizCoordSystem))     
        if(projwnd%iHorizCoordSystem/=22)then
          write(iunit(i),111)   '      Datum:',trim(aHorizDatum(projwnd%iHorizDatum))
          if(projwnd%iHorizZone/=0)then
            write(iunit(i),262) '      Zone:',projwnd%iHorizZone
          endif
        endif      
        write(iunit(i),111)     '      Units:',trim(aHorizUnits(projwnd%iHorizUnits))
      elseif(presconst)then
        write(iunit(i),*)
        write(iunit(i),111)     'Atmospheric Pressure Forcing:','ON'
        write(iunit(i),111)     '  Pressure Curve:',trim(presfile)
        write(iunit(i),102)     '  Ambient Atm Pressure:',ambatmpres,' Pa'
      endif
      
      !Precipitation and Evaporation
      if(rain_evap)then
        write(iunit(i),*)
        write(iunit(i),111)     'Precipitation/Evaporation:','ON'  
        if(len_trim(rf_filename)>0)then
          write(iunit(i),111)   '  File:',trim(rf_filename)
        else            
          write(iunit(i),102)   '  Precipitation:',rain,' m/hr'
          write(iunit(i),102)   '  Evaporation:',evap,' m/hr'
        endif
      endif
    enddo
    close(dgunit)
    
    return
    end subroutine met_print
    
!*************************************************************   
    subroutine windcurve_eval()
! Calcualtes wind stress terms
!
! written by Alex, USACE-ERDC-CHL
!*************************************************************             
    use size_def
    use flow_def, only: u,v
    use comvarbl, only: ramp,timehrs
    use met_def, only: wkappa,rhorat,nwnd_inc,nwtimes,&
          wndTimes,wndvalsx,wndvalsy,tauwindx,tauwindy,iwndlagr,&
          wndx,wndy,tauwx,tauwy,wndfac
    use met_lib, only: wind_drag_hsu
    use plagr_lib, only: plagr_fit,plagr_eval
    use prec_def          
    implicit none
    real(ikind) :: wndspd,val,Cd
    integer :: np,nti
    integer, parameter :: nb = 3 !>=np+1
    real(ikind) :: lb(nb)
             
    nti = 1 !Input polynomial order
    call plagr_fit(nwtimes,wndtimes,timehrs,nb,lb,nti,np,nwnd_inc)  
    wndx = plagr_eval(nwtimes,wndvalsx,nb,lb,np,nwnd_inc)
    wndy = plagr_eval(nwtimes,wndvalsy,nb,lb,np,nwnd_inc)
    
    !Convert wind speeds to shear stress
    if(iwndlagr==0)then !Eulerian reference frame
      wndspd = sqrt(wndx*wndx+wndy*wndy) !Eulerian wind speed            
      Cd = wind_drag_hsu(wndspd) !Calculate Wind Drag Coefficient, Hsu (1988)
      val = ramp*Cd*rhorat*wndspd  !Note ramp function
      tauwx = val*wndx
      tauwy = val*wndy
    endif       
    
    !Note: for the Lagrangian wind reference frame the shear stresses
    !  are a function of the current velocities and are therefore updated
    !  during the iterative loop of the implicit scheme 
    
    return
    end subroutine windcurve_eval
    
!*********************************************************    
    subroutine wind_lagrangian_update
! Updates the wind forcing variables when using a 
! Lagriangian reference frames.
! Only used by implicit solver.
! written by Alex Sanchez, USACE-CHL
!*********************************************************    
    use size_def, only: ncells
    use geo_def, only: areap
    use met_def, only: wndx,wndy,uwind,vwind,cdWndareap,&
        windconst,windvar,rhorat
    use met_lib, only: wind_drag_hsu
    use flow_def, only: u,v,us,vs
    use comvarbl, only: ramp
    use prec_def
    implicit none
    integer :: i
    real(ikind):: wndLagrx,wndLagry,Cd,wndspd
    
    if(windconst)then
!$OMP PARALLEL DO PRIVATE(i,wndspd,Cd,wndLagrx,wndLagry)  
      do i=1,ncells        
        wndLagrx = wndx-u(i)+us(i) !Note: should be changed to -(u(i)+u(msl)-us(s))
        wndLagry = wndy-v(i)+vs(i)  
        wndspd = sqrt(wndLagrx*wndLagrx+wndLagry*wndLagry) !Lagrangian wind speed            
        Cd = wind_drag_hsu(wndspd) !Calculate Wind Drag Coefficient, Hsu (1988)
        cdWndareap(i) = ramp*Cd*rhorat*wndspd*areap(i)  !Note ramp function
      enddo
!$OMP END PARALLEL DO 
    elseif(windvar)then
!$OMP PARALLEL DO PRIVATE(i,wndspd,Cd,wndLagrx,wndLagry)  
      do i=1,ncells        
        wndLagrx = uwind(i)-u(i)+us(i)
        wndLagry = vwind(i)-v(i)+vs(i)
        wndspd = sqrt(wndLagrx*wndLagrx+wndLagry*wndLagry) !Lagrangian wind speed   
        Cd = wind_drag_hsu(wndspd) !Calculate Wind Drag Coefficient, Hsu (1988)
        cdWndareap(i) = ramp*Cd*rhorat*wndspd*areap(i)  !Note ramp function
      enddo
!$OMP END PARALLEL DO           
    endif

    return
    end subroutine wind_lagrangian_update
    
!******************************************************************   
    subroutine prescurve_eval
! Interpolates the atmospheric pressure from an input time-series
! Computes a barometric tide and adds it to wse boundary conditions
!
! written by Alex, USACE-ERDC-CHL
!******************************************************************
    use bnd_def
    use comvarbl
    use flow_def
    use met_def, only: npres_inc,nptimes,presTimes,presvals,ambatmpres
    use plagr_lib
    use prec_def 
    use size_def
    implicit none    
    integer :: iwse,j,np,nti
    integer,parameter :: nb = 3 !>=np+1
    real(ikind) :: baro,pres,lb(nb)
    
    nti = 1 !Input polynomial order
    call plagr_fit(nptimes,prestimes,timehrs,nb,lb,nti,np,npres_inc)  
    pres = plagr_eval(nptimes,presvals,nb,lb,np,npres_inc) ![Pa]
    !baro = sum(lb(1:np+1)*Presvals(npres_inc:npres_inc+np))    
    baro = ramp*(ambatmpres - pres)/(grav*rhow) !Barometric tide [m]
    
    !Add barometric pressure to water level boundary conditions    
    !Tidal BC
    do iwse=1,nTHstr
      do j=1,TH_str(iwse)%ncells
        TH_str(iwse)%wsebnd(j) = TH_str(iwse)%wsebnd(j) + baro
      enddo  
    enddo
    
    !Single Water Level BC
    do iwse=1,nHstr  !for each cell string      
      do j=1,H_str(iwse)%ncells
        H_str(iwse)%wsebnd(j) = H_str(iwse)%wsebnd(j)+ baro
      enddo
    enddo ! end of each cell string        
    
    return
    end subroutine prescurve_eval             

!*************************************************************   
    subroutine windpres_grid
! Read or calculate wind/pressure grid coordinates and 
! convert to flow local coordinate system
! written by Alex Sanchez, USACE-CHL  
!*************************************************************   
#include "CMS_cpp.h"
    use met_def
    use geo_def, only: xOrigin,yOrigin,azimuth_fl,projfl
    use geo_lib, only: proj_horiz_conv
    use comvarbl, only: flowpath
    use const_def, only: deg2rad
    use diag_lib, only: diag_print_error, diag_print_message
    use diag_def, only: msg, msg2
    use prec_def
    implicit none
    integer :: i,j,k,nw,ierr
    real(ikind) :: swlat
    real(ikind), allocatable :: xw(:),yw(:)
    character(len=200) :: apath,aname
    character(len=10) :: aext
    logical :: foundfile

!11  format(6x,i4,16x,i4,23x,F6.0,32x,F6.0,44x,F6.0,58x,F6.0)   !This is so wrong - how did this ever work?
11 format(5x,i4, 6x,i4, 3x,F6.4, 3x,F6.4, 6x,F8.0, 6x,F8.0)
!         iLat   iLong     DX       DY     SWLat    SWLon 

   if(windformat==0)then !OceanWeather, Inc.
      open(wunit,file=windfile)
      read(wunit,*) !skip header
      read(wunit,11,iostat=ierr) nwindi, nwindj, wloninc, wlatinc, swlat, wlonmin
      if(ierr/=0)then
        call diag_print_error('Problem reading wind grid from Oceanweather file')
      endif
      close(wunit)
      wlatmax = swlat + (nwindi-1)*wlatinc
      !Note: It is assumed that the pressure grid and wind grid are the same
      projwnd%iHorizDatum = 1         !Horizontal Datum = NAD83
      projwnd%iHorizCoordSystem = 0   !Horizontal Coordinate System = GEOGRAPHIC
      projwnd%iHorizZone = 0          !Horizontal Zone = 0
      projwnd%iHorizUnits = 4         !Horizontal Units = DEGREES
      projwnd%iVertDatum = 9          !Vertical Datum = LOCAL
      projwnd%iVertUnits = 2          !Vertical Units = METERS
      projwnd%VertOffset = 0.0        !Vertical Offset from Datum
    elseif(windformat==1)then !NOAA Blended Sea Winds
      nwindi = 1440 !nlat
      nwindj = 719  !nlon
      wlatmax = 89.75 !deg
      wlonmin = 0.0   !deg
      wlatinc = 0.25  !deg
      wloninc = 0.25  !deg
      projwnd%iHorizDatum = 1         !Horizontal Datum = NAD83
      projwnd%iHorizCoordSystem = 0   !Horizontal Coordinate System = GEOGRAPHIC
      projwnd%iHorizZone = 0          !Horizontal Zone = 0
      projwnd%iHorizUnits = 4         !Horizontal Units = DEGREES
      projwnd%iVertDatum = 9          !Vertical Datum = LOCAL
      projwnd%iVertUnits = 2          !Vertical Units = METERS
      projwnd%VertOffset = 0.0        !Vertical Offset from Datum
    endif  
    
    if(nwindi==0 .or. nwindj==0)then
      call diag_print_error('Must specify wind grid size')
    endif  
    
    !Wind/Pressure Grid Coordinates
    allocate(xwind(nwindi,nwindj),ywind(nwindi,nwindj))
    if(wndpresloc)then !Read grid 
      !Wind/Press Location File
      call fileparts(windlocfile,apath,aname,aext)
      if(len_trim(apath)==0 .and. len_trim(flowpath)>0)then
        windlocfile = trim(flowpath) // windlocfile
      endif       
      inquire(file=windlocfile,exist=foundfile)
      if(.not.foundfile)then  !Compute grid
        msg='WARNING: Could not find wind location/coordinate file: '//trim(windlocfile)
        msg2='Computing Lat/Long from .win file.'
        call diag_print_message('',msg,msg2)
        do i=1,nwindi   !Lat
          do j=1,nwindj !Lon
            xwind(i,j) = wlonmin + wloninc*(j-1) 
            ywind(i,j) = wlatmax - wlatinc*(i-1)
          enddo !j loop
        enddo !k loop 
        
      else
        !xy file in local (global) coordinates
        open(cunit,file=windlocfile,status='old')
        do i=1,nwindi   !Lat
          do j=1,nwindj !Lon
            read(cunit,*) xwind(i,j),ywind(i,j)  
          enddo
        enddo
        close(cunit)   
      endif  
    else !Compute grid
      do i=1,nwindi   !Lat
        do j=1,nwindj !Lon
          xwind(i,j) = wlonmin + wloninc*(j-1) 
          ywind(i,j) = wlatmax - wlatinc*(i-1)
        enddo !j loop
      enddo !k loop 
    endif
    
    !If necessary convert projection    
    if(projwnd%iHorizCoordSystem/=projfl%iHorizCoordSystem)then
      nw = nwindi*nwindj  
      allocate(xw(nw),yw(nw))
      do i=1,nwindi !Lat
        do j=1,nwindj !Lon
          k = j + (i-1)*nwindj
          xw(k) = xwind(i,j)
          yw(k) = ywind(i,j)
        enddo
      enddo
      call proj_horiz_conv(projwnd,projfl,nw,xw,yw)
      do i=1,nwindi !Lat
        do j=1,nwindj !Lon
          k = j + (i-1)*nwindj
          xwind(i,j) = xw(k)
          ywind(i,j) = yw(k)
        enddo
      enddo
      deallocate(xw,yw)
    endif
    
    if(write_windpresgrid)then
      call fileparts(windpresgridfile,apath,aname,aext)
      if(len_trim(apath)==0 .and. len_trim(flowpath)>0)then
        windpresgridfile = trim(flowpath) // windpresgridfile
      endif
      open(454,file=windpresgridfile)
      do i=1,nwindi   !Lat
        do j=1,nwindj !Lon
          write(454,*) xwind(i,j),ywind(i,j),wndht
        enddo !j loop
      enddo !k loop
      close(454)
    endif
    
    return
    end subroutine windpres_grid

!*************************************************************   
    subroutine interp_coef_windpres2fl
! Calculates the interpolation coefficients for the wind and
! atmospheric pressure fields onto the flow grid
!
! written by Chris Reed, URS; Alex Sanchez, USACE-CHL    
!*************************************************************   
    use size_def
    use geo_def, only: xOrigin,yOrigin,azimuth_fl,x,y,xc,yc
    use met_def, only: nwindi,nwindj,xwind,ywind,  & 
        ijntpwnd,cntpwnd,wndfl_intpcoef_file
    use const_def, only: deg2rad
    use interp_lib, only: interp_coef_curv2pts
    use diag_lib
    use prec_def
    implicit none
    integer :: i,j,k,ierr
    integer, parameter :: iwndflver = 4
    integer :: nchk,nDchk,nichk,njchk,iver
    !integer, allocatable :: iintpcells(:)
    !real(ikind) :: xg,yg,cosang,sinang
    logical :: found
    
    inquire(file=wndfl_intpcoef_file,exist=found)
    if(found)then
!      open(93,file='Intpcoef_wndfl.bin',form='unformatted')
      open(93,file=wndfl_intpcoef_file,form='unformatted')
      read(93,iostat=ierr) iver
      if(ierr==0)then
        if(iver==iwndflver)then
          read(93,iostat=ierr) nchk,nDchk,nichk,njchk
          if(ierr==0)then
            if(nchk==ncells .and. nDchk==ncellsD .and. &
              nichk==nwindi .and. njchk==nwindj)then
              call diag_print_message('','Reading wind-to-flow interpolation coefficients','')  
              do i=1,ncellsD
                read(93,iostat=ierr) (ijntpwnd(j,i),j=1,2)
                if(ierr/=0) exit
                read(93,iostat=ierr) ((cntpwnd(j,k,i),j=1,2),k=1,2)
                if(ierr/=0) exit
              enddo
              if(ierr==0)then
                close(93)
                return
              endif
            endif  
          endif
        endif
      endif
      close(93,status='delete')
    endif 
    
    call diag_print_message('','Calculating wind-to-flow interpolation coefficients')  
    
    !open(55,file='wind.xyz') 
    !do i=1,nwindi
    !  do j=1,nwindj
    !    write(55,*) xwind(i,j),ywind(i,j),0.0  
    !  enddo
    !enddo  
    !close(55)
    
    !open(55,file='flow.xyz')     
    if(ncellpoly==0)then
      !do i=1,ncells
      !  write(55,*) x(i),y(i),1.0
      !enddo  
      !allocate(iintpcells(ncells))
      !iintpcells(1:ncells) = (/1:ncells/)
      !call interp_coef__curv2tel(nwindi,nwindj,xwind,ywind, &            
      !     ncellsD,xOrigin,yOrigin,azimuth_fl,x,y,    &
         !  ncells,iintpcells,ijntpwnd,cntpwnd)
      !deallocate(iintpcells)  
      call interp_coef_curv2pts(nwindi,nwindj,xwind,ywind, &            
           ncellsD,ncellsD,xc,yc,ijntpwnd,cntpwnd)
    else      
      !do i=1,ncells
      !  write(55,*) x(i),y(i),1.0
      !enddo
      call interp_coef_curv2pts(nwindi,nwindj,xwind,ywind, &            
           ncellsD,ncellsD,x,y,ijntpwnd,cntpwnd)  
    endif
    !close(55)
    
    call diag_print_message('Saving wind-to-flow interpolation coefficients','')  

    open(93,file=wndfl_intpcoef_file,form='unformatted')
    write(93) iwndflver
    write(93) ncells,ncellsD,nwindi,nwindj 
    do i=1,ncellsD
      write(93) (ijntpwnd(j,i),j=1,2)
      write(93) ((cntpwnd(j,k,i),j=1,2),k=1,2)
    enddo
    close(93)
    
 !   open(93,file='intpcoef_windpres2fl.dat',form='unformatted')
    !write(93,*) iwndflver
    !write(93,*) ncells,ncellsD,nwindi,nwindj 
 !   do i=1,ncellsD
 !     write(93,*) (ijntpwnd(j,i),j=1,2)
 !     write(93,*) ((cntpwnd(j,k,i),j=1,2),k=1,2)
 !   enddo
 !   close(93) 
    
    return
    end subroutine interp_coef_windpres2fl
    
!*************************************************************   
    subroutine windpresfield_eval
! Interpolates the wind field and calculates shear stresses
!
! written by Alex Sanchez, USACE-CHL
!*************************************************************         
    use comvarbl, only: timehrs,ramp
    use der_def, only: nder,goa
    use der_lib, only: der_grad_eval
    use flow_def, only: rhow
    use geo_def, only: azimuth_fl
    use interp_lib, only: interp_scal_curv2tel,interp_vec_curv2tel
    use met_def
    use met_lib
    use prec_def
    use size_def
    implicit none
    integer :: i,ierr
    real(ikind) :: val,wndspd,fac
    !real(ikind) :: scal(ncellsD)

113 format(' timehrs ',F10.3,' windhr ',F10.3)
    
    interface
      subroutine interp_met2hr(t,t1,t2,ramp,nD,nc,a,a1,a2,b,b1,b2)
        use prec_def
        implicit none
        integer,    intent(in) :: nD,nc
        real(ikind),intent(in) :: t,t1,t2,ramp
        real(ikind),intent(in) :: a1(nD),a2(nD)
        real(ikind),intent(out) :: a(nD)
        real(ikind),intent(in),optional :: b1(nD),b2(nD)
        real(ikind),intent(out),optional :: b(nD)
      end subroutine
    endinterface
    
    !Read Files
    fac = wndfac*wind_heightcorr(wndht) !Temporal and vertical coordinate factors
    select case(windformat)
    case(0) !OceanWeather, Inc.
      ierr = 0  
      do while(timehrs+1.0e-5>windhr2 .and. ierr==0)
        windhr1 = windhr2        
        call read_wind_owi(wunit,nwindi,nwindj,windhr2,wndspdx2,wndspdy2,ierr)    
!$OMP WORKSHARE
        wndspdx2 = wndspdx2*fac !Convert from different temporal scales
        wndspdy2 = wndspdy2*fac
        uwind1 = uwind2
        vwind1 = vwind2
!$OMP END WORKSHARE        
        !Spatial interpolation
        call interp_vec_curv2tel(nwindi,nwindj,ncellsD,ncellsD,ijntpwnd,cntpwnd,&
                 undefvel,lowvel,highvel,defvel,wndspdx2,uwind2,wndspdy2,vwind2)
        call rotate_vector(ncellsD,ncells,azimuth_fl,uwind2,vwind2)        
      enddo             
      if(presvar)then
        ierr = 0  
        do while(timehrs+1.0e-5>preshr2 .and. ierr==0)  
          preshr1 = preshr2
!$OMP WORKSHARE
          pressatm1 = pressatm2
!$OMP END WORKSHARE
          call read_pres_owi(punit,nwindi,nwindj,preshr2,atmpres2,ierr)          
          call interp_scal_curv2tel(nwindi,nwindj,ncellsD,ncellsD,ijntpwnd,cntpwnd,&
                 undefpres,lowpres,highpres,defpres,atmpres2,pressatm2)
        enddo              
      endif
      
    case(1) !Blended Winds
      undefvel = -9999.0  !Hard-coded
      do while(timehrs+1.0e-5>windhr2 .and. ierr==0)
        windhr1 = windhr2  ![hrs]
        windhr2 = windhr2 + 6.0 ![hrs]
        call read_blended_winds(wunit,windhr2,wndspdx2,wndspdy2,ierr)        
        preshr1 = windhr1 !Fore wind and pressure times to be the same
!$OMP WORKSHARE
        wndspdx2 = wndspdx2*fac
        wndspdy2 = wndspdy2*fac
        uwind1 = uwind2
        vwind1 = vwind2
!$OMP END WORKSHARE
        preshr2 = windhr2 !Force times to be the same
        call interp_vec_curv2tel(nwindi,nwindj,ncellsD,ncellsD,ijntpwnd,cntpwnd,&
                 undefvel,lowvel,highvel,defvel,wndspdx2,uwind2,wndspdy2,vwind2)        
        call rotate_vector(ncellsD,ncells,azimuth_fl,uwind2,vwind2)
      enddo  
      
    case(3) !Speed/Direction File (ADCIRC NWS=3 format)
      if(windhr1<0.0) windhr2 = -wtiminc !For first time step
      ierr = 0    
      do while(timehrs+1.0e-5>windhr2 .and. ierr==0)
        windhr1 = windhr2          
        call read_met_spddir(wunit,nwindi,nwindj,windhr2,wndspdx2,wndspdy2,ierr)
!$OMP WORKSHARE
        wndspdx2 = wndspdx2*fac
        wndspdy2 = wndspdy2*fac
        uwind1 = uwind2
        vwind1 = vwind2    
!$OMP END WORKSHARE      
        !Spatial Interpolation
        call interp_vec_curv2tel(nwindi,nwindj,ncellsD,ncellsD,ijntpwnd,cntpwnd,&
                 undefvel,lowvel,highvel,defvel,wndspdx2,uwind2,wndspdy2,vwind2)
        call rotate_vector(ncellsD,ncells,azimuth_fl,uwind2,vwind2)        
      enddo
      
    case(6,7) !WindX/WindY/Pressure Single File  (ADCIRC NWS=6 format) or Stressx,Stressy,Pres Single File 
      if(windhr1<0.0) windhr2 = -wtiminc !For first time step
      ierr = 0
      do while(timehrs+1.0e-5>windhr2 .and. ierr==0)
        windhr1 = windhr2  
        preshr1 = windhr1 !Fore wind and pressure times to be the same              
        call read_met_uvp(wunit,nwindi,nwindj,windhr2,&
               wndspdx2,wndspdy2,atmpres2,wtiminc,ierr)
!$OMP WORKSHARE
        wndspdx2 = wndspdx2*fac
        wndspdy2 = wndspdy2*fac
        uwind1 = uwind2
        vwind1 = vwind2
        pressatm1 = pressatm2         
!$OMP END WORKSHARE        
        preshr2 = windhr2 !Fore wind and pressure times to be the same        
        !Spatial Interpolation
        call interp_vec_curv2tel(nwindi,nwindj,ncellsD,ncellsD,ijntpwnd,cntpwnd,&
                 undefvel,lowvel,highvel,defvel,wndspdx2,uwind2,wndspdy2,vwind2)
        call interp_scal_curv2tel(nwindi,nwindj,ncellsD,ncellsD,ijntpwnd,cntpwnd,&
                 undefpres,lowpres,highpres,defpres,atmpres2,pressatm2)
        call rotate_vector(ncellsD,ncells,azimuth_fl,uwind2,vwind2)        
      enddo
      
    end select          
    
    !Wind temporal interpolation
    call interp_met2hr(timehrs,windhr1,windhr2,ramp,ncellsD,ncellsD,&
           uwind,uwind1,uwind2,vwind,vwind1,vwind2)
    
    !Spatial smoothing
    call smooth_flowgrid_vec(uwind,vwind,nsmoothwind)
    
    if(windformat==7)then
!$OMP PARALLEL DO PRIVATE(i,val,wndspd)       
        do i=1,ncells
          tauwindx(i) = uwind(i)
          tauwindy(i) = vwind(i)
        enddo
!$OMP END PARALLEL DO     
    else
      !Convert wind speeds to shear stresses   
      if(iwndlagr==0)then !Eulerian reference frame
!$OMP PARALLEL DO PRIVATE(i,val,wndspd)       
        do i=1,ncells
          wndspd = sqrt(uwind(i)*uwind(i)+vwind(i)*vwind(i))        
          val = rhorat*wind_drag_hsu(wndspd)*wndspd
          tauwindx(i) = val*uwind(i)
          tauwindy(i) = val*vwind(i)
        enddo
!$OMP END PARALLEL DO 
      endif
    endif
    
    if(.not.presvar) return
    
    !Atmospheric spatial gradients temporal interpolation
    call interp_met2hr(timehrs,preshr1,preshr2,ramp,ncellsD,ncellsD,&
             pressatm,pressatm1,pressatm2)
    !call smooth_flowgrid_scal(pressatm,pressatm,nsmoothpres)
    pressatm = pressatm + (1.0-ramp)*ambatmpres !Note: ramp already applied to pressatm above
    call der_grad_eval(goa,0,pressatm,pressatmdx,pressatmdy)
    !scal = pressatm-ambatmpres !To avoid precision errors
    !scal = (pressatm-ambatmpres)/1.0e6 !To avoid precision errors
    !call der_grad_eval(goa,0,scal,pressatmdx,pressatmdy)
    !pressatmdx = pressatmdx*1.0e6
    !pressatmdy = pressatmdy*1.0e6
    !!pressatmdx = sum(pressatmdx(1:ncells))/real(ncells,kind=ikind)  !For testing only
    !!pressatmdy = sum(pressatmdy(1:ncells))/real(ncells,kind=ikind)  !For testing only
    call smooth_flowgrid_vec(pressatmdx,pressatmdy,nsmoothpres) !Atm Pressure assumed to vary at scales much larger than the grid
    
    return
    end subroutine windpresfield_eval

!*************************************************************  
    subroutine interp_met2hr(t,t1,t2,ramp,nD,nc,a,a1,a2,b,b1,b2)
! written by Alex Sanchez, USACE-CHL    
!*************************************************************
    use prec_def        
    implicit none
    !Input/Output
    integer,intent(in) :: nD,nc
    real(ikind),intent(in) :: t,t1,t2,ramp
    real(ikind),intent(in) :: a1(nD),a2(nD)
    real(ikind),intent(out) :: a(nD)
    real(ikind),intent(in),optional :: b1(nD),b2(nD)
    real(ikind),intent(out),optional :: b(nD)
    !Internal Variables
    integer :: i
    real(ikind) :: fac,fac1,fac2

    if(present(b))then
      if(t<=t1)then
!$OMP PARALLEL DO PRIVATE(i)    
        do i=1,nc
          a(i) = ramp*a1(i)
          b(i) = ramp*b1(i)
        enddo
!$OMP END PARALLEL DO  
      elseif(t>t2)then
!$OMP PARALLEL DO PRIVATE(i)    
        do i=1,nc
          a(i) = ramp*a2(i)
          b(i) = ramp*b2(i)
        enddo
!$OMP END PARALLEL DO   
      else
        fac=(t-t1)/(t2-t1)
        fac=max(min(fac,1.0),0.0) !Avoids extrapolation
        fac1=ramp*(1.0-fac)  !Apply ramp here
        fac2=ramp*fac        !Apply ramp here
!$OMP PARALLEL DO PRIVATE(i)          
        do i=1,nc
          a(i) = fac1*a1(i) + fac2*a2(i)
          b(i) = fac1*b1(i) + fac2*b2(i)
        enddo
!$OMP END PARALLEL DO        
      endif    
    else
      if(t<=t1)then
!$OMP PARALLEL DO PRIVATE(i)    
        do i=1,nc
          a(i) = ramp*a1(i)
        enddo
!$OMP END PARALLEL DO  
      elseif(t>t2)then
!$OMP PARALLEL DO PRIVATE(i)    
        do i=1,nc
          a(i) = ramp*a2(i)
        enddo
!$OMP END PARALLEL DO   
      else
        fac=(t-t1)/(t2-t1)
        fac=max(min(fac,1.0),0.0) !Avoids extrapolation
        fac1=ramp*(1.0-fac)  !Apply ramp here
        fac2=ramp*fac        !Apply ramp here
!$OMP PARALLEL DO PRIVATE(i)          
        do i=1,nc
          a(i) = fac1*a1(i) + fac2*a2(i)
        enddo
!$OMP END PARALLEL DO        
      endif
    endif
    
    return
    end subroutine interp_met2hr
    
!***********************************************************************************    
    subroutine metsta_block
! Reads a Meteorological Station Block from the control file    
! written by Alex Sanchez, USACE-CHL
!***********************************************************************************
    use met_def, only: nMetSta,metstadriver,metsta,windsta,pressta,projwnd
    use geo_def, only: azimuth_fl
    use comvarbl, only: tjulhr0,flowpath
    use const_def, only: deg2rad
    use diag_lib
    use prec_def
    implicit none
    integer :: i,ii,ierr
    !real(ikind) :: pa,dum,ang,hrj
    character :: cardname*37,aext*10
    type(metstadriver), allocatable :: metstatemp(:)
    logical :: foundcard
    
    nMetSta=nMetSta+1
    if(nMetSta==1)then
      allocate(metsta(nMetSta))
    else      
      allocate(metstatemp(nMetSta-1))      
      do i=1,nMetSta-1
        metstatemp(i) = metsta(i)
      enddo
      deallocate(metsta)
      allocate(metsta(nMetSta))
      do i=1,nMetSta-1
        metsta(i) = metstatemp(i)
      enddo
      deallocate(metstatemp)
    endif        
    
    !Initialize variables
    metsta(nMetSta)%ntimes = 0
    metsta(nMetSta)%inc = 1
    metsta(nMetSta)%wndfac = 1.0
    metsta(nMetSta)%wndht = 10.0
    metsta(nMetSta)%name = ''
    metsta(nMetSta)%file = ''
    metsta(nMetSta)%path = ''
    metsta(nMetSta)%xsta = -999.0
    metsta(nMetSta)%ysta = -999.0
    metsta(nMetSta)%xlon = -999.0
    metsta(nMetSta)%ylat = -999.0
    metsta(nMetSta)%ntsi = 0
    metsta(nMetSta)%ntsw = 0
    metsta(nMetSta)%pow = 1.0   
    
    !Read Cards
d1: do ii=1,20
      foundcard = .true.
      read(77,*,err=171) cardname
      if(cardname(1:1)=='!' .or. cardname(1:1)=='#') cycle
      select case(cardname)
        case('NAME','STATION_NAME')
          backspace(77)
          read(77,*) cardname,metsta(nMetSta)%name
          
        case('MET_FILE','METEOROLOGICAL_FILE','FILE','WIND_DATASET','WIND_DATA')
          backspace(77)
          read(77,*) cardname,metsta(nMetSta)%file
          windsta = .true.
          
        case('WIND_PATH','PATH')
          backspace(77)
          read(77,*) cardname,metsta(nMetSta)%path
          windsta = .true.
          
        case('WIND_INPUT_CURVE','WIND_CURVE')
          backspace(77)
          read(77,*) cardname,metsta(nMetSta)%file,metsta(nMetSta)%path
          windsta = .true.
          
        case('WIND_SPEED_CURVE')
          backspace(77)
          read(77,*) cardname,metsta(nMetSta)%file
          call fileext(metsta(nMetSta)%file,aext)
          if(aext(1:2)=='h5')then
            backspace(77)
            read(77,*) cardname,metsta(nMetSta)%file,metsta(nMetSta)%path  
          endif
          windsta = .true.
          
        case('WIND_DIRECTION_CURVE')
          backspace(77)
          read(77,*) cardname,metsta(nMetSta)%file        !modified to treat the same as wind_speed curve - 10/03/2017        
          call fileext(metsta(nMetSta)%file,aext)
          if(aext(1:2)=='h5')then
            backspace(77)
            read(77,*) cardname,metsta(nMetSta)%file,metsta(nMetSta)%path  
          endif
          windsta = .true.
          
        case('PRESSURE_INPUT_CURVE')
          backspace(77)
          read(77,*) cardname,metsta(nMetSta)%file        !modified to treat the same as wind_speed and wind_direction curves - 10/03/2017
          call fileext(metsta(nMetSta)%file,aext)
          if(aext(1:2)=='h5')then
            backspace(77)
          read(77,*) cardname,metsta(nMetSta)%file,metsta(nMetSta)%path
          endif
          pressta = .true.
          
        case('COORDINATES')
          backspace(77)
          read(77,*) cardname,metsta(nMetSta)%xsta,metsta(nMetSta)%ysta  
          
        case('LATITUDE_LONGITUDE')
          backspace(77)
          read(77,*) cardname,metsta(nMetSta)%ylat,metsta(nMetSta)%xlon
          if(metsta(nMetSta)%xlon<180.0)then !For the Continential U.S.
            metsta(nMetSta)%xlon = -metsta(nMetSta)%xlon    
          endif          
          projwnd%iHorizCoordSystem = 0 !Geographic
          projwnd%iHorizUnits = 4       !Degrees
      
        case('LATITUDE')
          backspace(77)
          read(77,*) cardname,metsta(nMetSta)%ylat
          projwnd%iHorizCoordSystem = 0 !Geographic
          projwnd%iHorizUnits = 4       !Degrees
      
        case('LONGITUDE')
          backspace(77)
          read(77,*) cardname,metsta(nMetSta)%xlon
          if(metsta(nMetSta)%xlon<180.0)then !For the Continential U.S.
            metsta(nMetSta)%xlon = -metsta(nMetSta)%xlon    
          endif          
          projwnd%iHorizCoordSystem = 0 !Geographic
          projwnd%iHorizUnits = 4       !Degrees
      
        case('WIND_FACTOR','WIND_VELOCITY_FACTOR')
          backspace(77)
          read(77,*) cardname,metsta(nMetSta)%wndfac          
        
        case('ANEMOMETER_HEIGHT')
          call card_scalar(77,'m','m',metsta(nMetSta)%wndht,ierr)
          !metsta(nMetSta)%wndht=max(metsta(nMetSta)%wndht,1.0)
          
        case('TIME_SMOOTH_ITER','SMOOTHING_ITERATIONS')
          backspace(77)
          read(77,*) cardname,metsta(nMetSta)%ntsi
        
        case('TIME_SMOOTH_WIDTH','SMOOTHING_WIDTH')
          backspace(77)
          read(77,*) cardname,metsta(nMetSta)%ntsw
          if(mod(metsta(nMetSta)%ntsw,2)/=0)then !Make sure it is odd
            metsta(nMetSta)%ntsw = metsta(nMetSta)%ntsw + 1
          endif   
        
        case('POWER_PARAMETER')
          backspace(77)
          read(77,*) cardname,metsta(nMetSta)%pow
          metsta(nMetSta)%pow=min(max(metsta(nMetSta)%pow,1.0),5.0) !Limit values
          
        case('MET_STATION_END','METEOROLOGICAL_STATION_END','MET_STA_END','END')
          exit d1
        
        case('HORIZONTAL_PROJECTION_BEGIN','HORIZ_PROJ_BEGIN')
          call proj_horiz_block(77,projwnd)
          
        case default
          foundcard = .false.
      end select
    enddo d1
171 continue    
    
    return
!-------------------------------------------------------------------------
787 close(445)
    call diag_print_error('End-of-File Encountered during read of Meteorological file: ',metsta(nMetSta)%file)

!-------------------------------------------------------------------------
984 close(445)
    call diag_print_error('Error Reading Meteorological file: ',metsta(nMetSta)%file)

    end subroutine metsta_block
    
!***********************************************************************    
    subroutine metsta_init
! Initializes the Meteorological Station variables
! written by Alex Sanchez, USACE-CHL
!***********************************************************************    
#include "CMS_cpp.h"
    use size_def
    use geo_def, only: xOrigin,yOrigin,azimuth_fl,x,y,xc,yc,&
        projection,projfl
    use geo_lib, only: proj_horiz_conv
    use comvarbl, only: flowpath,tjulday0,tmax
    use const_def, only: deg2rad
#ifdef XMDF_IO
    use in_xmdf_lib, only: read_dataseth5
#endif
    use in_lib, only: read_xys,read_tsd
    use met_def
    use met_lib, only: read_stdmetfile,wind_heightcorr
    use diag_lib
    use prec_def
    implicit none
    integer :: i,ista,niter,nwin,ndat
    real(ikind) :: ang,wspd,wdir,tjuldaybegwnd,tjulend,fac
    real(ikind), allocatable :: xmetsta(:),ymetsta(:)
    real(ikind), pointer :: dat(:,:)
    character(len=200) :: apath,aname
    character(len=10) :: aext
    character(len=50) :: atype
    !type(projection) :: projgeo
    logical :: foundfile
            
    if(projwnd%iHorizCoordSystem==0)then !Convert from Geo to Flow Coordinate system  
      projwnd%iHorizDatum = projfl%iHorizDatum
      allocate(xmetsta(nMetSta),ymetsta(nMetSta))
      do ista=1,nMetSta          
        xmetsta(ista)=metsta(ista)%xlon  
        ymetsta(ista)=metsta(ista)%ylat        
      enddo
      call proj_horiz_conv(projwnd,projfl,nMetSta,xmetsta,ymetsta)        
      do ista=1,nMetSta          
        metsta(ista)%xsta=xmetsta(ista)
        metsta(ista)%ysta=ymetsta(ista)
      enddo
      deallocate(xmetsta,ymetsta)
    !else !Convert from Flow Coordinate system to Geographic
    !  projwnd = projfl  
    !  projgeo%iHorizDatum = projfl%iHorizDatum
    !  projgeo%iHorizCoordSystem = 0 !Geographic
    !  projgeo%iHorizUnits = 4       !Degrees    
    !  do ista=1,nMetSta          
    !    xmetsta(ista)=metsta(ista)%xsta  
    !    ymetsta(ista)=metsta(ista)%ysta        
    !  enddo
    !  call proj_horiz_conv(projwnd,projgeo,nMetSta,xmetsta,ymetsta)        
    !  do ista=1,nMetSta          
    !    metsta(ista)%xlon=xmetsta(ista)
    !    metsta(ista)%ylat=ymetsta(ista)
    !  enddo
    endif        
     
    
    !Load in Data
    do ista=1,nMetSta
      !Check if Input Files Exist  
      call fileparts(metsta(ista)%file,apath,aname,aext)
      if(len_trim(aname)>0)then
        if(len_trim(apath)==0 .and. len_trim(flowpath)>0)then
          metsta(ista)%file = trim(flowpath) // metsta(ista)%file
        endif
        if(len_trim(metsta(ista)%file)>0)then
          inquire(file=metsta(ista)%file,exist=foundfile)
          if(.not.foundfile)then
            call diag_print_error('Could not Find Meteorological File: ',metsta(ista)%file)
          endif 
        endif
      endif
      !Read Files
      if(aext(1:2)=='h5')then
#ifdef XMDF_IO
        call read_dataseth5(metsta(ista)%file,metsta(ista)%path,'Times',     metsta(ista)%ntimes,metsta(ista)%times)
        call read_dataseth5(metsta(ista)%file,metsta(ista)%path,'Magnitude', metsta(ista)%ntimes,metsta(ista)%wndvalsx)
        call read_dataseth5(metsta(ista)%file,metsta(ista)%path,'Direction', metsta(ista)%ntimes,metsta(ista)%wndvalsy)
#else
        call diag_print_error('Cannot read wind time series from *.h5 file without XMDF libraries')
#endif
      elseif(aext(1:1)==' ')then !Standard Met File
        tjulend = tjulday0 + tmax/24.0
        metsta(ista)%file = metsta(ista)%name
        if(len_trim(metsta(ista)%path)==0)then
          metsta(ista)%path = windpath
        endif
        call read_stdmetfile(metsta(ista)%file,metsta(ista)%path,tjulday0,tjulend,&
          metsta(ista)%ntimes,metsta(ista)%times,metsta(ista)%wndvalsx,metsta(ista)%wndvalsy)  
      elseif(aext(1:3)=='tsd')then  !TSD file
        ndat = 2 
        call read_tsd(metsta(ista)%file,aname,atype,ndat,metsta(ista)%ntimes,&
           tjuldaybegwnd,metsta(ista)%times,dat)
        metsta(ista)%times = metsta(ista)%times + tjuldaybegwnd - tjulday0
        allocate(metsta(ista)%wndvalsx(metsta(ista)%ntimes))
        allocate(metsta(ista)%wndvalsy(metsta(ista)%ntimes))
        metsta(ista)%wndvalsx = dat(:,1)
        metsta(ista)%wndvalsy = dat(:,2)
        deallocate(dat)
      elseif(aext(1:3)=='xys')then  !XYS file
        call read_xys(metsta(ista)%file,metsta(ista)%ntimes,&
            metsta(ista)%times,metsta(ista)%wndvalsx) !Speed 
        call read_xys(metsta(ista)%path,metsta(ista)%ntimes,&
            metsta(ista)%times,metsta(ista)%wndvalsy) !Direction 
      endif
    enddo    
    
    !Check for NaN Values
    do ista=1,nMetSta
      do i=1,metsta(ista)%ntimes
        if(metsta(ista)%wndvalsx(i)<-90.0 .or. metsta(ista)%wndvalsy(i)<-900.0)then !spd=-99.0,dir=-999.0
           metsta(ista)%wndvalsx(i) = 0.0
        endif
      enddo  
    enddo
    
    !Convert to Components
    ang=(180.0+azimuth_fl)*deg2rad !Note rotation from meteological convention 
    do ista=1,nMetSta
      fac=metsta(ista)%wndfac*wind_heightcorr(metsta(ista)%wndht)
      do i=1,metsta(ista)%ntimes
        wspd=metsta(ista)%wndvalsx(i)*fac
        wdir=metsta(ista)%wndvalsy(i)*deg2rad+ang         !Convert to local coordinate system and in radians
        metsta(ista)%wndvalsx(i)=wspd*sin(wdir)           !Note sin and cos are switched because of coordinate system, 0=top
        metsta(ista)%wndvalsy(i)=wspd*cos(wdir) 
        !metsta(ista)%pa(i)=pa/100.0  !Convert from hPa to Pa (N/m^2)
      enddo  
    enddo  
    
    !Temporal smoothing of Components
    do ista=1,nMetSta    
      niter=metsta(ista)%ntsi
      nwin=metsta(ista)%ntsw
      if(niter>0 .and. nwin>1)then
        call moving_average(niter,nwin,metsta(ista)%ntimes,metsta(ista)%wndvalsx)
        call moving_average(niter,nwin,metsta(ista)%ntimes,metsta(ista)%wndvalsy)
      endif
    enddo
    
    !Inverse Distance (Shepard) Interpolation Coefficients       
    allocate(cntpmet(ncells,nMetSta)) !Interpolation coefficients
    do i=1,ncells
      do ista=1,nMetSta
        cntpmet(i,ista) = sqrt((metsta(ista)%xsta-xc(i))**2 &
               +(metsta(ista)%ysta-yc(i))**2)**(-metsta(ista)%pow)
      enddo      
      cntpmet(i,:)=cntpmet(i,:)/sum(cntpmet(i,:))
    enddo    
    
    return
    end subroutine metsta_init
    
!***********************************************************************    
    subroutine metsta_eval
! Evaluates the wind field by interpolating in the time and space the
! meteorological stations. 
! written by Alex Sanchez, USACE-CHL
!***********************************************************************    
    use size_def
    use comvarbl, only: timehrs,ramp
    use met_def
    use met_lib, only: wind_drag_hsu
    use plagr_lib
    use prec_def
    implicit none
    integer :: i,ista,np,nti
    real(ikind) :: wndspd,val
    integer,parameter :: nb = 3 !>=ni+1
    real(ikind) :: lb(nb)

    !Temporal Interpolation
    nti = 1 !linear (quadratic can cause problems for noisy data)
    do ista=1,nMetSta
      call plagr_fit(metsta(ista)%ntimes,metsta(ista)%times,&
         timehrs,nb,lb,nti,np,metsta(ista)%inc) 
      metsta(ista)%wndx = sum(lb(1:np+1)*metsta(ista)%wndvalsx(metsta(ista)%inc:metsta(ista)%inc+np))
      metsta(ista)%wndy = sum(lb(1:np+1)*metsta(ista)%wndvalsy(metsta(ista)%inc:metsta(ista)%inc+np))     
      metsta(ista)%wndx = ramp*metsta(ista)%wndx
      metsta(ista)%wndy = ramp*metsta(ista)%wndy
    enddo

!$OMP PARALLEL
    !Spatial Interpolation
!$OMP DO PRIVATE(i,ista)
    do i=1,ncells
      uwind(i)=0.0; vwind(i)=0.0
      do ista=1,nMetSta
        uwind(i) = uwind(i) + cntpmet(i,ista)*metsta(ista)%wndx
        vwind(i) = vwind(i) + cntpmet(i,ista)*metsta(ista)%wndy
      enddo
    enddo
!$OMP END DO
    
    !Convert wind speeds to shear stresses   
    if(iwndlagr==0)then    
!$OMP DO PRIVATE(i,val,wndspd)
      do i=1,ncells
        wndspd = sqrt(uwind(i)*uwind(i)+vwind(i)*vwind(i))        
        val = rhorat*wind_drag_hsu(wndspd)*wndspd
        tauwindx(i) = val*uwind(i)
        tauwindy(i) = val*vwind(i)
      enddo
!$OMP END DO
    endif
!$OMP END PARALLEL

    return
    end subroutine metsta_eval
        