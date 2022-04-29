!===================================================================
module hot_def
! Hot start variable definitions
!===================================================================
    use prec_def
    implicit none
    save
    
    !Initial Condition (input)
    logical :: coldstart   !Coldstart for current simulation
    
    !type hot_init_time_type
    !  real(ikind) :: start  !hrs
    !  real(8)     :: reference !Must be double
    !endtype hot_init_time_type
    !
    !type hot_init_var_type
    !  type(hot_init_time_type) :: time
    !  character(len=200) :: filename,datapath
    !  integer :: istat
    !endtype hot_init_var_type
    !
    !type hot_init_type
    !  type(hot_init_var_type) :: wse
    !  type(hot_init_var_type) :: vel
    !  type(hot_init_var_type) :: wet
    !  type(hot_init_var_type) :: depth
    !  type(hot_init_var_type) :: sal
    !  type(hot_init_var_type) :: Ctk
    !  type(hot_init_var_type) :: db
    !  type(hot_init_var_type) :: pbk
    !endtype hot_init_type
    !
    !type(hot_init_type) :: hotinit !hotinit%wse%time
    
    !type init_cond_type
    !  real(ikind) :: start  !hrs
    !  real(8)     :: reference !Must be double
    !  character(len=200) :: filename,datapath
    !  integer :: istat
    !endtype init_cond_type
    !type(init_cond_type) :: wseic
    !type(init_cond_type) :: velic
    !type(init_cond_type) :: wetic
    !type(init_cond_type) :: depthic
    !type(init_cond_type) :: fluxic
    !type(init_cond_type) :: Ctkic
    !type(init_cond_type) :: pbkic
    !type(init_cond_type) :: dbic
    !type(init_cond_type) :: saltic
    
    real(ikind) :: ictime      !Initial condition time
    real(ikind) :: icwsetime   !WSE initial condition time
    real(ikind) :: icveltime   !Velocity initial condition time
    real(ikind) :: icwettime   !Velocity initial condition time
    real(ikind) :: icdepthtime !Depth initial condition time
    real(ikind) :: icsaltime   !Salinity initial condition time
    real(ikind) :: icheatime   !Temperature initial condition time
    character(len=200) :: icfile,icpath           !Initial condition file and path (input hot start file)
    character(len=200) :: icwsefile,icwsepath     !WSE initial condition
    character(len=200) :: icvelfile,icvelpath     !Velocity initial condition  
    character(len=200) :: icwetfile,icwetpath     !Wet/dry initial condition  
    character(len=200) :: icdepthfile,icdepthpath !Depth initial condition 
    character(len=200) :: icsalfile,icsalpath     !Salinity initial condition 
    character(len=200) :: icheatfile,icheatpath   !Temperature initial condition
    character(len=200) :: icpbkfile,icpbkpath     !Bed layer thickness initial condition
    character(len=200) :: iclayfile,iclaypath     !Bed layer thickness initial condition
    logical :: setconc2eq  !Set the concentrations to the equilibrium concentration
    logical :: icpres  !Read in water pressure
    logical :: icwse   !Read in wse
    logical :: icvel   !Read in velocities
    logical :: icwet   !Read in wetting/drying
    logical :: icflux  !Read in cell fluxes
    logical :: icconc  !Read in sediment concentrations
    logical :: icpbk   !Read in sediment fractions
    logical :: iclay   !Read in bed layer thickness
    integer :: icdper  !Percentile diameter input mode, 0-none, 1-d35/d50/d90, 2-d50/sigma
    logical :: icsal   !Read in salinity
    logical :: icheat  !Read in temperature
    logical :: icdepth !Read in depths
    real(8) :: reftimehot
    
    !Hot start (output)    
    logical :: hot_out     !Write a hot start file
    logical :: hot_recur   !Output a hot start file at a recurring interval
    logical :: hot_pres
    logical :: hot_eta
    logical :: hot_vel
    logical :: hot_wet
    logical :: hot_flux
    logical :: hot_timehr  !Output a hot start file at a specific time
    real(ikind) :: hotdt   !Recurring hot start interval
    real(ikind) :: hottime !Time to write hot start file
    real(ikind) :: timeout !Last output time
    real(ikind) :: hstarttime   ! MEB - 05/18/2012
    character(len=15)  :: HotName, AutoHotName      !Hold just the name of the file, not the extension
    character(len=200) :: hotfile,hotfile2,hotpath  !Output hot start file
    character(len=200) :: autohotfile,autohotpath   !Output auto hot start file 2/3/2017 MEB

    logical :: add_duration_HS   !When starting from initial condition, add the specified duration_run to the startup time?  Ex. - Run another month from the end of the IC file.  MEB 04/29/2022

    
end module hot_def