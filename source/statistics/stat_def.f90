!========================================================= 
module stat_def
! Module for Global statistics library
! written by Alex Sanchez, USACE-CHL
!=========================================================
    use prec_def      
    implicit none
    save
    
    !File name
    character(len=200) :: statfile,statpath
    
    !Global stats
    logical :: calc_stats
    integer :: nstatcount
    real(ikind) :: tstat(3)
    character(len=200) :: stat_file
    
    !Flow stats
    logical :: flowstats
    integer :: nflowstatcount
    real(ikind) :: tflowstat(3)  !(start time [hrs], end time [hrs], calc interval [hrs])
    real(ikind), allocatable :: uvmax(:),umax(:),vmax(:)       !Maximum currents [m/s]
    real(ikind), allocatable :: tuvmax(:)                      !Time of maximum currents [m/s]
    real(ikind), allocatable :: uvr(:),ur(:),vr(:)             !Residual currents (net) [m/s]
    real(ikind), allocatable :: uvmaxgrad(:)                   !Maximum current gradients [1/s]
    real(ikind), allocatable :: uvcurv(:)                      !Average curvature of velocity field, 1/Rc [1/m]
    real(ikind), allocatable :: etamax(:)                      !Maximum water level (inundation) [m]
    real(ikind), allocatable :: tetamax(:)                     !Time of maximum water level (inundation) [m]
    real(ikind), allocatable :: etamin(:)                      !Minimum water level [m]
    real(ikind), allocatable :: tetamin(:)                     !Time of minimum water level [m]
    real(ikind), allocatable :: etamean(:)                     !Mean water level [m]
    real(ikind), allocatable :: etamaxgrad(:)                  !Maximum water level gradients (-)
    real(ikind), allocatable :: hydper(:)                      !Average hydroperiod (fraction of time wet) (-)
    real(ikind), allocatable :: rspavg(:),rsuavg(:),rsvavg(:)  !Average Normalized Residuals
    real(ikind), allocatable :: rspmax(:),rsumax(:),rsvmax(:)  !Mean Normalized Residuals
    real(ikind), allocatable :: cflxmax(:),cflymax(:)          !Courant numbers CFLx = (abs(u)+sqrt(grav*h))*dtime/dx
    real(ikind), allocatable :: Rexmax(:),Reymax(:)            !Maximum Reynolds numbers Rex = h*Ux/mu
    real(ikind), allocatable :: Frxmax(:),Frymax(:)            !Maximum Froude numbers Frx = Ux/sqrt(grav*h)    
    real(ikind), allocatable :: umaxgradx(:),umaxgrady(:)      !u current gradients [1/s]
    real(ikind), allocatable :: vmaxgradx(:),vmaxgrady(:)      !v current gradients [1/s]
    real(ikind), allocatable :: umaxcurvx(:),umaxcurvy(:)      !u current curvature [1/s/m]
    real(ikind), allocatable :: vmaxcurvx(:),vmaxcurvy(:)      !v current curvature [1/s/m]
    
    !Grid quality (calculated inide flow stats)
    real(ikind), allocatable :: qhgradxmax(:),qhgradymax(:)    !Grid quality indicators based on depth slopes
    real(ikind), allocatable :: qhcurvxmax(:),qhcurvymax(:)    !Grid quality indicators based on depth curvatures
    real(ikind), allocatable :: qcelmax(:)                     !Grid quality indicator based on wave length
    real(ikind), allocatable :: qcflxmax(:),qcflymax(:)        !Grid quality indicator based on Courant numbers      
    real(ikind), allocatable :: qFrxmax(:),qFrymax(:)          !Grid quality indicator based on Courant numbers      
    real(ikind), allocatable :: qareachg(:)                    !Grid smoothness (static)
    
    !Sediment Stats
    logical :: sedstats
    integer :: nsedstatcount
    real(ikind) :: tsedstat(3)  !(start time [hrs], end time [hrs], calc interval [hrs])
    real(ikind), allocatable :: qtxmax(:),qtymax(:),qtmax(:) !Maximum transpor rate [kg/m/s]
    real(ikind), allocatable :: qtxn(:),qtyn(:),qtn(:)       !Net transport rate [kg/m/s]
    real(ikind), allocatable :: qtg(:)                       !Gross transport rate [kg/m/s]   
    real(ikind), allocatable :: rsCtm(:)                     !Mean of Normalized Residuals
    
    !Salinity Stats
    logical :: salstats
    integer :: nsalstatcount
    real(ikind) :: tsalstat(3)  !(start time [hrs], end time [hrs], calc interval [hrs])
    real(ikind), allocatable :: salmax(:) !Maximum [ppt]
    real(ikind), allocatable :: salmin(:) !Minimum [ppt]
    real(ikind), allocatable :: salavg(:) !Average [ppt]
    real(ikind), allocatable :: rsSalm(:) !Residuals
    
    !Temperature Stats
    logical :: heatstats
    integer :: nheatstatcount
    real(ikind) :: theatstat(3)  !(start time [hrs], end time [hrs], calc interval [hrs])
    real(ikind), allocatable :: heatmax(:) !Maximum [deg C]
    real(ikind), allocatable :: heatmin(:) !Minimum [deg C]
    real(ikind), allocatable :: heatavg(:) !Average [deg C]
    real(ikind), allocatable :: rsheatm(:) !Residuals

    !Wave Stats
    logical :: wavestats
    integer :: nwavestatcount
    real(ikind) :: twavstat(3)  ![start time (hrs)] [end time (hrs)] [calc interval (hrs)]
    real(ikind), allocatable :: whgtmax(:)    !Maximum significant wave height [m]
    real(ikind), allocatable :: whgtmean(:)   !Mean significant wave height [m]
    real(ikind), allocatable :: whgtmeanx(:) !Mean significant wave height x-component [m]
    real(ikind), allocatable :: whgtmeany(:) !Mean significant wave height x-component [m]
    real(ikind), allocatable :: wpermean(:)  !Mean wave period [s]
    real(ikind), allocatable :: ursellmax(:) !Maximum Ursell number Ur = H*L^2/h^3 < 100 for linear waves
    real(ikind), allocatable :: wdissmax(:)   !Maximum wave dissipation 
    
endmodule stat_def