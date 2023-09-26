!===================================================================
module flow_def
! Flow module variable definitions
!
! written by Weiming Wu, NCCHE, and Alex Sanchez, USACE-CHL
!===================================================================
  use prec_def, only: ikind
  
  implicit none    
  save  

  !Gravity
  real(ikind) :: grav     !Gravitational constant, m/s^2
  real(ikind) :: gravinv  !Inverse graviational constant [s^2/m]
  real(ikind) :: sqrtgrav !Square-root of graviational constant [m^(1/2)/s]
    
  !Water properties
  integer     :: idensit      !Equation of state for water density, 0-none, 1-Ekart (1958), 2-UNESCO (1981), 3-Fofonoff (1985)
  integer     :: iviscos      !Equation used to estimate the water kinematic viscosity, 0-user specified, 1-temp, 2-temp and salinity
  integer     :: isalt        !Method for estimating salinity, 0-none(default), 1-user specified, 2-initial condition
  integer     :: iheat        !Method for estimating temperature, 0-none(default), 1-user specified, 2-initial condition
  real(ikind) :: rhow         !Water density [kg/m^3]
  real(ikind) :: watertemp    !Water temperature [ºC]
  real(ikind) :: viscos       !Kinematic viscosity [m^2/s]   
  real(ikind) :: watersalt !Water salinity used to estimate density and viscosity [ppt]
        
  !Wetting and Drying
  logical :: ponding,narrowchannels
  real(ikind) :: hmin !Wetting and drying depth [m]
  real(ikind) :: hdry !Wetting and drying criteria [m]
  real(ikind) :: hdry1,hdry2 !Wetting and drying depth for narrow channels 1 and 2 cells wide (only for cartesian grids) [m]
  integer, allocatable :: iwet(:),iwet1(:),iwet2(:) !1-wet, 0-dry
    
  !Types
  integer, allocatable :: iextrap(:),icorner(:)
    
  !Shared Explicit and Implicit Hydro variables
  real(ikind), allocatable :: eta(:),uv(:),flux(:,:),flux1(:,:)
  real(ikind), allocatable :: h(:),h1(:),h2(:),hk(:,:),dhx(:),dhy(:)
  real(ikind), allocatable :: p(:),p1(:),p2(:),pk(:,:),dpx(:),dpy(:)
  real(ikind), allocatable :: u(:),u1(:),u2(:),uk(:,:),dux(:),duy(:)
  real(ikind), allocatable :: v(:),v1(:),v2(:),vk(:,:),dvx(:),dvy(:)
  real(ikind), allocatable :: hustar(:),dhustarx(:),dhustary(:),hustark(:,:)
  real(ikind), allocatable :: hvstar(:),dhvstarx(:),dhvstary(:),hvstark(:,:)
  real(ikind), allocatable :: fluxstar(:,:)
  real(ikind), allocatable :: h0(:)
  real(ikind), allocatable :: maxeta(:)  !Added MEB  9/20/2021
    
  !Wave flux velocities
  logical :: waveflux !Wave mass flux
  real(ikind), allocatable :: us(:),vs(:)  !Stokes velocities
    
  !Implicit solver
  real(ikind), allocatable :: acoef(:,:),ap(:)
  real(ikind), allocatable :: pp(:),ppk(:,:),dppx(:),dppy(:) !pressure correction 
  real(ikind), allocatable :: su(:),sv(:),sp(:),spu(:),spv(:)
  real(ikind), allocatable :: ssu0(:),ssv0(:),spuv0(:),sspp0(:),sppp0(:)
  real(ikind), allocatable :: Hu(:),dHux(:),dHuy(:),dHuxm(:),dHuym(:)
  real(ikind), allocatable :: Hv(:),dHvx(:),dHvy(:),dHvxm(:),dHvym(:)
  real(ikind), allocatable :: apuareap(:),dapuareapx(:),dapuareapy(:),dapuareaxm(:),dapuareaym(:)
  real(ikind), allocatable :: sumu(:) !,dsumux(:),dsumuy(:)
  real(ikind), allocatable :: rsp(:),rsu(:),rsv(:)
    
  !Explicit solver
  real(ikind), allocatable :: Huk(:,:),Hvk(:,:)
  real(ikind), allocatable :: detax(:),detay(:)
  !real(ikind), allocatable :: limdetax(:),limdetay(:),limdhx(:),limdhy(:)
  !real(ikind), allocatable :: limdux(:),limduy(:),limdvx(:),limdvy(:)
    
  !Predictor-corrector scheme (adapative time stepping)
  real(ikind), allocatable :: hpred(:),upred(:),vpred(:)
    
  !Numerical methods for implicit schemes
  type numimp
    integer     :: ntsch     !Temporal scheme
    integer     :: ndsch     !Advection scheme
    integer     :: nrecon    !Reconstruction type, 0-none,1-cell-based,2-face-based
    logical     :: skewcor   !Skewness correction
    logical     :: pred_corr !Predictor-corrector
    real(ikind) :: ctsch,ctsch1,ctsch2,wtsch   !Temporal scheme variables
  endtype numimp
  type(numimp)  :: usch,vsch,ppsch

  logical     :: volcor,flowvolbal
  real(ikind) :: volH2Ocumstg,volH2Ocumbnd,volH2Ocumflux
  real(ikind) :: volH2Ocumrain,volH2Ocumevap,volH2Ocumnorm
    
  !Eddy viscosity
  logical           :: crossdiff
  integer           :: mturbul    !Turbulence model ID
  character(len=32) :: aturb(0:5) !Turbulence model name
  real(ikind)       :: cviscon,cvisbot,cvishor,cviswav,cviswavbrk,cvismax
  real(ikind)       :: cvishor2areaavg2,cvishor2areaavg
  real(ikind), allocatable :: vis(:),viskfl(:,:),visk(:,:)
  real(ikind), allocatable :: diswall(:)
  data aturb /'CONSTANT',&      !0
              'FALCONER',&      !1
              'PARABOLIC',&     !2
              'SUBGRID-WU',&    !3
              'MIXING-LENGTH',& !4
              'SUBGRID'/        !5    
    
  !Coriolis
  integer :: icoriolisplane !1=f-plane,2=beta-plane
  real(ikind) :: fcoriolis,betacoriolis
  real(ikind), allocatable :: fc(:)
    
  !Temporal Wave Interpolation      
  logical :: wave_interp = .true.     !If false, then wave input is constant until next wave update.   MEB  01/22/2014
end module flow_def
