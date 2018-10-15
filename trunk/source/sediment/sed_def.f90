!======================================================================
module sed_def
! Sediment transport variables
!======================================================================
    use prec_def
    implicit none
    save
    
    !Activate sediment transport
    logical :: sedtrans
    
    !Sediment transport model    
    integer :: isedmodel
    character(len=30) :: asedmodel(5)
    data asedmodel /'NET',&      !1, Nonequilibrium total-load    
                    'A-D',&      !2, Advection-diffusion suspended load, equilibrium bed load 
                    'EXNER',&    !3, Equilibrium total load       
                    'WATANABE',& !4, Equilibrium Watanabe
                    'LUND_CIRP'/ !5, Equilibrium Lund-CIRP     
    
    !Scaling factors
    real(ikind) :: scalebed, scalesus, scalemorph
    
    !Avalanching
    logical :: do_aval
    integer :: nmaxaval       !# of iterations
    real(ikind) :: a_repose   !Angle of repose = critical angle
    real(ikind) :: relax_aval !Relaxation coefficient
        
    !Sediment size-class and properties
    integer :: nsed !# of sediment size classes
    integer :: iws  !Default method for calculating sediment fall velocity
    logical :: singlesize
    real(ikind) :: singleD50
    real(ikind), allocatable :: diam(:)       !Characteristic grain size diameter [m]
    real(ikind), allocatable :: diamlim(:)    !Size class limits
    real(ikind), allocatable :: logdiamlim(:) !=log(diamlim)
    real(ikind), allocatable :: dstar(:)      !Nondimensional grain size diameter [-]
    real(ikind), allocatable :: wsfall(:)     !Sediment fall velocity
    real(ikind), allocatable :: thetacr(:)    !Critical Shields parameter
    real(ikind), allocatable :: taucr(:)      !Critical shear stress
    real(ikind), allocatable :: coreyshape(:) !Corey shape factor [-]
    real(ikind) :: sedshape      !Corey shape factor for all
    !integer, allocatable :: itrace(:)   !Size classes to use as tracers
    
    !Sediment properties    
    real(ikind) :: rhosed   !Sediment density [kg/m^3]
    real(ikind) :: poros    !Sediment porosity [-]
    real(ikind) :: solid    !1-poros
    real(ikind) :: s1grav,d2dstar !Internal variables
    real(ikind) :: specgrav !Specific gravity
    
    !Hiding and exposure
    integer :: ihidexpform  !Hiding and exposure formula 1-Egiazaroff,2-Parker,3-Wu
    character(len=30) :: ahidexpform(0:5)
    data ahidexpform /'NONE',&           !0, None
                      'EGIAZAROFF',&     !1, Egiazaroff (1965)
                      'PARKER',&         !2, Parker et al. (1982) and others 
                      'WU',&             !3, Wu et al. (2000)
                      'ASHIDA_MICHIUE',& !4, Ashida and Michiue (1980)
                      'HAYASHI'/         !5, Hayashi et al. (1980)    
    real(ikind) :: mhe      !Hiding and exposure coefficient [-]
    real(ikind), allocatable :: varsigma(:,:)    
    real(ikind) :: varsigmamin,varsigmamax
    
    !Sediment size class
    type sed_size_class
      character(len=100) :: name !Name  
      integer :: idiam          !1-Characteristic diameter, 2-Diameter bounds
      real(ikind) :: diam       !Characteristic diameter
      real(ikind) :: diamlim(2) !Diameter bounds
      integer :: iws            !Fall velocity formula, 0-User-specified, 1-Soulsby, 2-Wu-Wang
      real(ikind) :: wsfall     !Fall velocity [m/s]
      integer :: icr            !Inicipient motion, 0-User-specified thetacr,1-User-specified taucr, 2-Soulsby
      real(ikind) :: thetacr    !Critical shields parameter
      real(ikind) :: taucr      !Critical shear stress [N/m^2]
      real(ikind) :: shape      !Corey shape factor
      !logical :: existcoverage
      !logical, allocatable :: coverage(:) !true - exists, false - does not exist
      !character(len=200) :: covpath,covfile 
    endtype sed_size_class
    type(sed_size_class),allocatable :: sedclass(:)
    
    !Transport Variables
    real(ikind) :: Awatan,Awidg
    integer :: icapac
    character(len=10) :: acapac(6)
    data acapac / 'LUND-CIRP',& !1
                  'VAN_RIJN',&  !2
                  'WATANABE',&  !3
                  'SOULSBY',&   !4
                  'WU',&        !5
                  'CSHORE'/     !6 bdj
    !Mixing    
    logical :: sedconstmix
    real(ikind) :: schmidt,cmixsed,cmixbedload
    
    !Concetnration profile
    integer :: iconcprof
    character(len=20) :: aconcprof(4)
    data aconcprof /'EXPONENTIAL',& !1
                    'ROUSE',&       !2 
                    'LUND-CIRP',&   !3, Same as exponential, used by explicit model
                    'VAN_RIJN'/     !4, Same as Rouse distribution, used by explicit model    
    
    !Concentrations    
    real(ikind), allocatable :: Ctk(:,:)     !Total load sediment concentrations for each size class at time t
    real(ikind), allocatable :: Ctk0(:,:)    !Total load sediment concentrations for each size class at time 0.0
    real(ikind), allocatable :: Ctk1(:,:)    !Total load sediment concentrations for each size class at time t-dt 
    real(ikind), allocatable :: Ctk2(:,:)    !Total load sediment concentrations for each size class at time t-2*dt
    real(ikind), allocatable :: Ctkstar(:,:) !Total-load mass concentration capacities for each size class
    real(ikind), allocatable :: CtstarP(:,:) !Total-load potential mass concentation capacity for each size class
    real(ikind), allocatable :: cak(:,:)     !Concentration at reference height for each size class
    real(ikind) :: Cteqmax                   !Maximum total-load equilibrium mass concentration physical limit is Cteqmax= (1.0-poro)*rhosed
    real(ikind), allocatable :: dCtkx(:,:),dCtky(:,:) !Total-load sediment concentration gradients
    real(ikind), allocatable :: Etkstar(:,:)       !Entrainment rate for each size class              !Wu
    real(ikind), allocatable :: EtstarP(:,:)       !Entrainment potential rate for each size class    !Wu
    
    !Other variables    
    real(ikind), allocatable :: rsk(:,:)           !Fraction of suspended sediments for each size class    
    real(ikind), allocatable :: epsvk(:,:)         !Vertical sediment mixing coefficient
    real(ikind), allocatable :: rsCtk(:)           !Normalized residual for each concentration
    real(ikind), allocatable :: rsCtkmax(:)        !Max Normalized residual for all concentrations
    
    !Bed load variables
    real(ikind), allocatable :: qbk(:,:),qbk1(:,:) !Bed load transport magnitude for each size class
    real(ikind), allocatable :: ubk(:,:),ubk1(:,:) !Bed load transport velocity
    real(ikind), allocatable :: pbkstar(:,:)       !Bed layer exchange fraction
    
    !Cumulative variables    
    real(ikind), allocatable :: Ct(:),rs(:),Ctstar(:)    
    real(ikind), allocatable :: qtx(:),qty(:)      !Total load sediment transports [kg/m/s]
    
    !Load correction factor    
    integer :: ibt  !0-Constant btk, 1-Exponential concentration profile, 2-Rouse concentration profile
    real(ikind) :: betatot  !Constant total-load correction factor 
    real(ikind), allocatable :: btk(:,:),btk1(:,:),btk2(:,:) !Total load correction factors (cell,size_class)
    real(ikind), allocatable :: bsk(:,:)                     !Suspended load correction factor (cell,size_class)

    !Adaptation coefficient
    integer :: iadapttot,iadaptsus,iadaptbed
    real(ikind), allocatable :: alphat(:)          !Total load Adaptation coefficient
    real(ikind) :: Ltot,Lerotot,Ldeptot,Ttot
    real(ikind) :: Lbed,Lerobed,Ldepbed,Tbed,fbed
    real(ikind) :: Lsus,Lerosus,Ldepsus,Tsus,alphasus
    real(ikind) :: ftot,fsus  !added bdj
    !Total-load 
    character(len=32) :: atotm(6)
    data atotm /'CONSTANT_LENGTH',&           !1
                'CONSTANT_TIME',&             !2
                'SPATIAL_LENGTH',&            !3
                'CONSTANT_EROS_DEPO_LENGTH',& !4
                'MAX_BED_SUSP_LENGTH',&       !5 
                'WGHT_AVG_BED_SUSP_LENGTH'/   !6
    character(len=32) :: asusm(8)
    data asusm /'CONSTANT_LENGTH',&           !1
                'CONSTANT_TIME',&             !2
                'CONSTANT_ALPHA',&            !3
                'CONSTANT_EROS_DEPO_LENGTH',& !4
                'ARMANINI_DISILVIO',&         !5
                'LIN',&                       !5 
                'GALLAPPATTI',&               !7
                'XBEACH'/                     !8
    character(len=32) :: abedm(5)
    data abedm /'CONSTANT_LENGTH',&           !1
                'CONSTANT_TIME',&             !2
                'SPATIAL_LENGTH',&            !3, Not implimented yet ******
                'CONSTANT_EROS_DEPO_LENGTH',& !4
                'DEPTH_DEPENDANT'/            !5
    
    !Variable Adaptation Length
    logical :: variableLtot
    real(ikind), allocatable :: vLtot(:),vLbed(:),vLsus(:)
    character(len=200) :: aLtotpath,aLtotfile
    character(len=200) :: aLbedpath,aLbedfile
    character(len=200) :: aLsuspath,aLsusfile        

    !Bed elevation and bed change variables
    logical :: calcmorph !Determines whether to update morphology (bed) or not 
    real(ikind) :: tStartMorph !Starting time for morphology change calculation
    real(ikind), allocatable :: zb1(:)
    real(ikind), allocatable :: dzb(:),dzbk(:,:)
    real(ikind) :: dzbmax
    
    !Bed-material Composition
    logical :: calcbedcomp
    real(ikind) :: tStartBedComp !Starting time for morphology change calculation
    character(len=32) :: apbkinp(0:6)    
    data apbkinp /'NONE',&                 !0
                  'D50_SIGMA',&            !1, Spatially variable, uniform in depth
                  'D16_D50_D84',&          !2, Spatially variable, uniform in depth
                  'D35_D50_D90',&          !3, Spatially variable, uniform in depth *** NOT WORKING *******
                  'SIZE_CLASS_FRACTIONS',& !4, Specifies a constant distribution curve, Spatially uniform, and uniform in depth
                  'FRACTIONS_DATASET',&    !5, Specifies dataset for pbk for cell and layer
                  'PERCENTILES'/           !6, Specifies percentiles, uniform in depth
    real(ikind), allocatable :: pbk(:,:,:),pbk1(:,:) !Fractional bed composition [-]
    real(ikind), allocatable :: d50(:),d90(:)
    character(len=200) :: pbkpath,pbkfile  !Bed composition file and path
    character(len=2) :: aconc(20)  !Size class number string
    
    !Bed Layer Thickness
    integer :: nlay,nlayinp
    character(len=32) :: adbinp(0:3)
    data adbinp /'NONE',&            !0,No information on bed layers, use default values
                 'CONSTANT',&        !1, Applies the same layer thickness to all layers, everywhere
                 'LAYER_THICKNESS',& !2, Specifies the thickness of each layer and is constant spatially
                 'LAYER_DATASET'/    !3, Specifies the thickness of each layer and cell
    logical :: mixlayconst
    real(ikind) :: dmconst !Constant mixing layer thickness [m]
    real(ikind) :: db1min  !Minimum active layer thickness [m]
    real(ikind) :: dbmin   !Minimum layer thickness [m]
    real(ikind) :: dbmax   !Maximum layer thickness for new layers. Existing layers are not affected [m]
    real(ikind), allocatable :: db(:,:),db1(:,:)     !Bed layer thickness (cell,layer) [m]
    
    !Bed-slope effect on transport rate, either through critical shear (Dey), directly on transport rate (Bailard), or the bed shear stress (Wu)
    integer :: ibedslope   !0-none, 1-Dey (2001), 2-Bailard (1981), 3-Wu et al. (2000)
    character(len=20) :: abedslope(0:3)
    data abedslope /'NONE','DEY','BAILARD','WU'/
    real(ikind) :: betaslope !Bed-load slope coefficient for Bailard (1981) method, default = 1.6 ~= 1/tan(a_repose), multiplies transport by slpcor=1.0-betaslope*slp
    real(ikind) :: effslope  !Efficiency factor for Bailard formula
    
    !Bed-slope diffusion term
    logical :: do_bedslope
    real(ikind) :: dcoeff  !Bed-slope diffusion coefficient
    real(ikind), allocatable :: Sb(:,:) !Bed-slope term

    !Hard-bottom
    logical :: hardbottom
    integer :: nhard
    integer, allocatable :: idhard(:)
    real(ikind), allocatable :: hardbed(:)
    real(ikind), allocatable :: hardzb(:)    
    character(len=200) :: hbfile,hbpath   
    
    !Boundary conditions
    integer :: isedinflowbc,nsedflux,nsedsource,nsedbc
    real(ikind) :: Qtotin,facQtotin    
    real(ikind), allocatable :: sedbnd(:,:)
    !Flux and source boundaries
    type sed_driver !Multiple sed rates, one for each cell in string
      integer :: ibndstr
      integer :: ncells
      integer, allocatable :: cells(:)
      integer, allocatable :: faces(:)
      integer :: ntimes,inc
      integer :: ibctype !1-Single transport rate, 2-fractional transport rate, 3-single source/sink, 4-fractional source/sink
      real*4, allocatable :: times(:)    !Times
      real*4, allocatable :: val(:,:)    !Data Values(time,class)
    endtype sed_driver
    character(len=200) :: sedfluxfile,sedsourcefile
    character(len=200), allocatable :: sedfluxpath(:),sedsourcepath(:)
    logical :: sedfluxascii
    
    type(sed_driver), allocatable :: sed_str(:)    
    
    !-- Total Sediment Flux Boundary Condition -----
    type sed_totflux_type
      integer :: ibndstr !Boundary string ID number
      integer :: ntimes
      integer :: inc
      integer :: nti  !Order of temporal interpolation
      integer :: nsw  !Smoothing window width
      integer :: nwi  !Smoothing iterations
      real(ikind),allocatable :: times(:)
      real(ikind),allocatable :: dat(:)  !Data Values(time) [kg/s]
      real(ikind),allocatable :: bnd(:,:)  !Interpolated and fractional flux data at boundary 
    endtype sed_totflux_type
    type(sed_totflux_type),allocatable :: sed_totflux(:)
    
    !-- Fractional Sediment Flux Boundary Condition -----
    type sed_fracflux_type
      integer :: ibndstr !Boundary string ID number
      integer :: ntimes
      integer :: inc
      integer :: nti  !Order of temporal interpolation
      integer :: nsw  !Smoothing window width
      integer :: nwi  !Smoothing iterations
      real(ikind),allocatable :: times(:)
      real(ikind),allocatable :: dat(:,:)  !Data Values(time,class), [kg/s/class]
      real(ikind),allocatable :: bnd(:,:)  !Interpolated and fractional flux data at boundary 
    endtype sed_fracflux_type
    type(sed_fracflux_type),allocatable :: sed_fracflux(:)
    
    !--- Capacity Boundary Condition ------
    type sed_bnd_capac_type
      integer :: ibndstr
      real(ikind) :: facCapac !Over/Under-Loading factor, used for all size classes
    endtype sed_bnd_capac_type
    type(sed_bnd_capac_type),allocatable :: sed_bnd_capac(:)
    
    !--- Total Source/Sink -----------------
    type sed_bnd_srcsnk_type
      integer :: ibndstr !Boundary string ID number
      integer :: ntimes
      integer :: inc
      integer :: nti  !Order of temporal interpolation
      integer :: nsw  !Smoothing window width
      integer :: nwi  !Smoothing iterations
      integer :: ncells !Number of source/sink cells
      integer :: cells  !Cell ID's
      real(ikind),allocatable :: times(:)
      real(ikind),allocatable :: dat(:) !Data Values(time,class) [kg/s]
    endtype sed_bnd_srcsnk_type
    type(sed_bnd_srcsnk_type), allocatable :: sed_bnd_srcsnk(:)
    
    !Constant wave parameters for testing and idealized cases	   !From Chris' Code
    !logical :: constant_waves
    
    !Percentile Datasets
    integer :: nperinp !Number of input percentile datasets
    logical :: variableD50,transD50,calcd50,constd50    
    integer, parameter :: nperdiam = 13
    integer, parameter :: iper(nperdiam) = (/5,10,16,20,30,35,40,50,65,75,84,90,95/) !Percentiles 5,10,etc
    logical :: outperdiam(nperdiam)    
    integer :: ipd(99)      
    type percentile_diameter      
      character(len=200) :: file,path      !Percentile diameter file and path
      real(ikind), pointer :: dper(:) !Percentile diameter dataset
      logical :: inp                  !Input toggle, true if read in correctly
    endtype percentile_diameter
    type(percentile_diameter) :: perdiam(nperdiam)
    
    !Bed layer/composition
    type bed_layer
      !Thickness  
      integer :: idbinp        !Bed layer thickness input mode
      real(ikind) :: dbconst   !Spatially constant initial bed layer thickness [m]
      character(len=200) :: dbfile, dbpath     !Bed layer thickness file and path
      !Composition
      logical :: inppbk        !Input bed layer composition      
      integer :: ipbkinp       !Bed composition input mode
      real(ikind), allocatable :: pbconst(:) !Fractional bed composition
      real(ikind) :: geostddev !Geometric standard deviation for layer (>=1.0)      
      type(percentile_diameter) :: perdiam(nperdiam)
      character(len=200) :: pbkfile, pbkpath
    endtype bed_layer
    type(bed_layer),allocatable :: bedlay(:)
    
    !Global Sediment Balance
    logical :: sedbalance
    type sedvol_type
      real(ikind) :: storage
      real(ikind) :: bedchange
      real(ikind) :: erosion  
      real(ikind) :: deposition
      real(ikind) :: balance
      real(ikind) :: boundary
      real(ikind) :: gross
      real(ikind) :: error
    endtype
    type(sedvol_type),allocatable :: sedvol(:,:)
    type(sedvol_type),allocatable :: sedvolcur(:)
    type(sedvol_type),allocatable :: sedvolcum(:)
    type(sedvol_type),allocatable :: sedvolnet(:)
    !real(ikind) :: volcumbnd
    
    !Wave-induced sediment transport. Includes wave asymetry, roller, and undertow
    logical :: wavesedtrans
    real(ikind) :: scaleonshore,scaleoffshore,scalewaveasym,scaleroller,scaleundertow
    real(ikind) :: fsr,fba,fsm,fbm
    real(ikind),allocatable :: Qws(:,:),QwsP(:,:)

    !Solver variables
    integer :: itersed
    logical :: sedcouple
    integer :: maxitersed,itermaxzb,maxitersed0
    real(ikind) :: tolCtk,errCtk,errCtk0
    real(ikind) :: tolpbk,errpbk,errpbk0
    
    !Erosion of dry cells
    type erosdry_type
      logical :: calc
      real(ikind) :: fac
      real(ikind) :: slopemin
      real(ikind) :: slopemax
    endtype erosdry_type
    type(erosdry_type) :: erosdry
    
    !Output
    character(len=200) :: soutfile,boutfile !Detailed sediment output for mixed sediments
    
    !Cohesive Sediment
    logical :: cohesivesed, consolidation
    integer :: methcoherodcr
    real(ikind) :: wsfallmax, cohk1, cohk2, cohcp, cohn, cohr
    real(ikind) :: cohsalkmax, cohsalcp, cohsaln, cohtaubp, cohturbn1, cohturbn2, cohturbk1
    real(ikind) :: cohdepmax0,cohdepmin0, coherodm0, coherodcr0,coherodn
    real(ikind) :: pcmax,pcmin,rhobedcoh0,arhobed,prhobed,rhobedcoh1yr,betarhobed,rhobednoncoh
    real(ikind), allocatable :: cohdepmax(:),cohdepmin(:), wsfallcohsed(:),coherodm(:), coherodcr(:)  ! Max and min critical shear stresss for deposition
    real(ikind), allocatable :: cohbsxy(:)
    real(ikind), allocatable :: dm(:),dmk(:,:),dbms(:,:),dbms1(:,:),rhobed(:,:),rhobedcoh(:,:)
    real(ikind), allocatable :: tconsolid(:,:),tconsolid0(:)
    real(ikind) :: coherodcratrho0,coherodcrtau,rhobedcohercr0,coherodcrn
    
    endmodule sed_def
    