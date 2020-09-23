!===================================================================
module prec_def 
! Specifies th precision for the CMS-Flow floating variables
!===================================================================    
    implicit none
    save
    integer, parameter :: ikind=4  !Alex, 4-single, 8-double  
endmodule prec_def
    
!===================================================================
module cms_def
!===================================================================
    use prec_def
    
    implicit none
    save
    
    logical :: cmsflow, cmswave
    logical :: inlinewave = .false.
    integer :: noptset !1-CMS-Waves only, 2-CMS-Flow only, 3-CMS-Waves and CMS-Flow, 4-CMS-Flow and Wave Input datasets
    integer :: n2Dor3D
    integer :: nsteer
    integer :: nopttime !0-constant interval, 1-time-based
    integer :: noptwse  !0-zero,1-last,2-tidal,3-tidal plus variation
    integer :: noptvel  !0-zero,1-last
    integer :: noptzb   !0-none,1-last depth,2-bed change,3-interp from datasets
    integer :: nradsm   !Radiation stress smoothing iterations
    integer :: nbrksm   !Wave breaking index smoothing iterations
    integer :: ndissm   !Wave breaking dissipation smoothing iterations
    integer :: npersm   !Wave period smoothing iterations 
    integer :: noptxtrpfl,noptxtrpwav !1-User specified, 2-Automatic       
    real(ikind) :: xtrpdistwav,xtrpdistfl !Extrapolation distances (0.0 for none)
    real(ikind) :: dtsteer                !Steering interval in seconds
    integer :: nspectimes
    real(8),allocatable :: spectimes(:)
    real(ikind) :: wavedisstol
    real(ikind) :: waveradfac
    real(8) :: timebegin  !CPU begin time, Note: Keep double precision
    real(8) :: timestart  !Julian start time, Note: Keep double precision
    real(8) :: timenow    !Julian time for speed monitor
    character(len=200) :: wavsimfile,wavepath,wavename  !Alex
    character(len=200) :: radpath,wavpath,perpath,dirpath,disspath
    logical :: ignore_brk_restr = .false.
    
    !For keeping a list of unknown cards and values
    character(len=200) :: aValue
    integer :: nCards = 0
    type card_type
      character(len=200) :: cardname
      character(len=400) :: value      !Rest of the line after the cardname
    endtype card_type
    type(card_type), allocatable :: cardList(:)
    
    !For keeping list of datasets, their paths, and their dimensions together
    integer :: ndsets
    type dset_type !Node string driver
      character(len=200) :: filename
      character(len=200) :: path
      integer :: ndim
    endtype dset_type
    type(dset_type), allocatable :: dsetList(:)
end module cms_def

!===================================================================
module const_def
! Constants
!===================================================================
    use prec_def
    implicit none
    save

    real(ikind),parameter :: pi = acos(-1.0_ikind)
    real(ikind),parameter :: twopi = 2.0_ikind*pi
    real(ikind),parameter :: deg2rad = pi/180.0_ikind
    real(ikind),parameter :: rad2deg = 180.0_ikind/pi
    real(ikind),parameter :: sqrttwo = sqrt(2.0_ikind)
    real(ikind),parameter :: cappa = 0.41_ikind !Von Karman Constant
    real(ikind),parameter :: eps = 1.0e-20  
    real(ikind),parameter :: small = 1.0e-15
    real(ikind),parameter :: great = 1.0e15    
    real(ikind),parameter :: omega = twopi/86400.0_ikind   !Earth's rotation speed, 7.2722e-5 rad/s
    real(ikind),parameter :: earthRadius = 6371000.0_ikind !Earth's radius, m
    
endmodule const_def

!===================================================================
module comvarbl
!===================================================================
    use prec_def
    implicit none
    save
    
    !Code version - Moved assignment to CMS2D_Main.f90
    real      :: Version   !CMS version
    integer   :: Revision  !Revision number    !MEB 09/15/20 - Switched back to integer
    logical   :: release
    integer   :: major_version, minor_version, bugfix
    character(len=10) :: rdate
    character(len=10) :: machine     !06/11/2019 MEB - added to denote type of machine compiled on
    
    real      :: input_ver        !CMS-Flow Input Version. Read from cmcards file    !renamed from 'ver', 6/28/2016
    real      :: SMS_ver   = -1   !SMS Version used to write these files             !added MEB  03/19/2020
   
    !Used for testing only
    integer :: iflow1D !makes flow 1D    
    
    !File names
    character(len=200) :: flowpath             !Path of CMS-Flow files
    character(len=200) :: ctlfile              !Control File
    character(len=200) :: mpfile               !Model parameters file     
    character(len=200) :: advfile              !Advanced card file - added MEB 07/29/2019
    logical            :: read_adv = .false.   !Switched to True if file exists and is read

    !Timing
    !Note: Use double precision for time variables
    integer :: ntime,nramp,nspinup,nprt
    real(ikind) :: dtime,dtime1 !Time steps [sec]
    real(8) :: dtimebeg
    real(ikind) :: stime        !Start time [sec]
    real(ikind) :: ctime        !Current simulation time [sec]
    real(ikind) :: etime        !Elapsed time output interval [sec]
    real(ikind) :: ctime1       !Simulation time of last speed check [sec]
    real(8) :: deltime,timesecs !Use double precision
    real(ikind) :: timehrs      !Current simulation time [hrs]
    !real(8) :: htime            !Hot start time [sec]
    real(ikind) :: stimet       !Total simulation time [sec]
    real(ikind) :: tmax         !Total simulation time [hrs]
    real(ikind) :: ramp         !Ramp function [-]
    real(ikind) :: rampdur      !Ramp duration [hrs]
    real(8) :: reftime          !needs to be double because of XMDF subroutines   
    integer :: iyr,imo,iday,ihr,imin,isec !Gregorian calendar date and time
    integer :: jday             !Julian day
    integer :: ical(6)          !Calendar date [yyyy,mm,dd,HH,MM,SS]
    real(ikind) :: tjulday0,tjulhr0,tjulhryr,dtj     !Julian days and hours
    
    !Skewness correction
    logical:: skewcor
    
    !Solution scheme
    integer:: nfsch  !0-Implicit, 1-Explicit
    
    !Implicit scheme
    logical :: pred_corr !Predictor-corrector scheme (adapative time stepping)
    logical :: dtvar   !Variable time-stepping 
    integer:: ntsch     !Temporal scheme, 1-First Order, 2-Second Order
    integer:: ndsch     !Advection scheme
    character(len=32) :: advsc(0:18) !Advection scheme
    data advsc /'NONE',&                           !0, No advection
      'UPWIND','HYBRID','POWERLAW','EXPONENTIAL',& !1-4; 1st order schemes
      'HLPA',&                                     !5 2nd order NVD sheme
      'GAMMA','CUBISTA','ALVSMART','HOAB',&        !6-9; 2nd-3rd order NVD schemes (gamma family)
      'MINMOD','MUSCL','CHARM','HCUS','HQUICK',&   !10-14; 2nd order TVD schemes
      'OSPRE','KOREN','VAN_ALBADA','SMART'/        !15-18; 2nd order TVD schemes
    integer:: nsolv  !Matrix solver    
    integer:: niter,nprint,maxit,nswp(5)
    integer:: nswp0(5),maxit0,mtime,jtime
    real(ikind):: ctsch,ctsch1,ctsch2,wtsch   !Temporal scheme variables
    real(ikind):: rmommaxp,rmomminp,rmomtargetp,rmomratiop,rmomabschgp,rmomrelchgp !Thresholds
    real(ikind):: rmommaxuv,rmomminuv,rmomtargetuv,rmomratiouv,rmomabschguv,rmomrelchguv !Thresholds    
    real(ikind):: rmom(5),rmom0(5),rmom1(5),rp,ruv,rchp,rchuv
    real(ikind):: relax(10),relaxsor,facpp
    real(ikind):: presmax,velmax !Used for convergence checking
    !real(ikind):: facblend !Blending factor for deferred corrections (0.0<=facblend<=1.0). Placeholder. Not Implemented yet
    
    !Explicit scheme    
    integer:: norder  !Order of solution
    integer:: nriem   !Reimann solver
    !integer:: nbeta   !Slope limiter        
    !character(len=32) :: abeta(0:6)
    character(len=4) :: ariem(3) !Riemann solver names
    data ariem /'ROE',& !1
                'HLL',& !2
                'HLLC'/ !3   
    
    !Parallel processing
    integer:: nthr,nthrmax
    
    !Extremal values. Used for checking stability
    integer:: idu,idv,idp,idcr
    real(ikind):: uxtrm,vxtrm,pxtrm,crmax
    
    !Strings
    character(len=200) :: casename
    
endmodule comvarbl
