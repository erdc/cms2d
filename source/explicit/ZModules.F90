!********************************************************************************
MODULE EXP_Global_def
#include "CMS_cpp.h"
#ifdef XMDF_IO
USE XMDF, only: DOUBLE, XF_GUID_STRINGLENGTH
#endif
use prec_def

!CHARACTER*10,PARAMETER :: Rdate    = '8/18/2010'
!REAL, PARAMETER        :: Version  = 3.75
!INTEGER, PARAMETER     :: Revision = 2          !<0 -- Beta,  >=0 -- Release
!CHARACTER*25 :: SIMLABEL = 'Solution'           !Alex, 11/06/09, moved from declare_global_file.fi

INTEGER :: NUM_EXT_N,NUM_EXT_S,NUM_EXT_E,NUM_EXT_W  
integer :: ncn,nce,ncs,ncw,ncne,ncse,ncsw,ncnw
  
!ADDING GUID TO OUTPUT
#ifdef XMDF_IO
CHARACTER(XF_GUID_STRINGLENGTH) :: A_GUID
#endif

CHARACTER*200 :: GUIDPATH
INTEGER GUID_ID
CHARACTER*200 GDC_PATH,GDT_PATH  

!TIMING (WALL CLOCK)
integer :: time_array_0(8), time_array_1(8)
real    :: outinterval = -1.0

!=== Flow parameters ===
logical :: advect, mixing
real(ikind) :: dt,ADVCOEF,FAC_DW,FAC_UW,drydep    !Rhow

logical radstr,waves
real :: wind_fac = 1.0

!=== Sediment transport parameters ===
logical :: sedtransEXP, hardbottom, &
            do_avalanche, &   !
           watanabe,lundcirp,adeq,adneq,totneq
integer :: iripple, netadv, netbedslope, netdiff, netbed, netsus, nettot, netcapac, netaval, netbctype
real(ikind) :: dtsed, dtmorph, A0, dtsalt, dtheat

!for single grain size
real(ikind) :: InflowLoadFac
real(ikind) :: WSFALL,TAUCR,USTCR,THETACR,DSTAR
real(ikind), allocatable :: VWSFALL(:),VTAUCR(:),VDSTAR(:)

!salinity
real(ikind), allocatable :: RHOPrim(:)  !CR - 01/21/2009

INTEGER MAXUNIT,NDUMMY

!basic grid parameters
LOGICAL, ALLOCATABLE :: active(:,:)

!hydrodyanmic variables
real(ikind), allocatable :: etan(:),qx(:),qxn(:),qy(:),qyn(:),Cdx(:),Cdy(:)
real(ikind), allocatable :: uE(:),vE(:)

!advection and mixng variables
real(ikind), allocatable :: Fuu(:),Fuv(:),Gvv(:),Gvu(:),ADVECTX(:),ADVECTY(:) 
integer, allocatable :: FG_A_cells(:,:),FG_M_cells(:)  
integer num_FG_M_cells,num_FG_A_cells          
   
!RAINFALL CARDS
character*120 RF_FILENAME
logical RAINFALL
integer rf_unit
real(ikind) RF_FRAC_RO,rain_time
real(ikind) rain,evap  !expected in meters at hourly intervals

LOGICAL       DIAG  
CHARACTER*200 DIAG_STRING
CHARACTER*6 ADUM6

REAL(IKIND) TIME  
integer, allocatable :: linktodummies(:),linktodummiesTel(:,:)
integer num_linktodummies
      
integer iadv,imix
integer isedform
!integer i,j   !Commented by Alex. Considered dangerous
real(ikind) thetac
         
CONTAINS 

    !****************************************************
    ! double precision rounding function - 12/06/07 meb
    !   two inputs - value, and #digits to round off to
    !   returns - rounded value to specified prec_def
    !****************************************************
    real*8 function rround8 (X,P)    
    use prec_def
    implicit none
    real(ikind) X    !real value passed in
    integer K,P   !digits of prec_def
    real*8 PR,R
    
    PR=10.d0**P      
    K = NINT(PR*(X-AINT(X))) 
    R = AINT(X) + K/PR
    RROUND8 = R
    
    end function
  
    !****************************************************
    ! single prec_def rounding function - 12/06/07 meb
    !   two inputs - value, and #digits to round off to
    !   returns - rounded value to specified prec_def
    !****************************************************
    real function rround (X,P)
    implicit none
    real X    !real value passed in
    integer K,P   !digits of prec_def
    real PR,R
    
    PR=10.0**P      
    K = NINT(PR*(X-AINT(X))) 
    R = AINT(X) + K/PR
    RROUND = R
    
    end function

END MODULE EXP_Global_def 
!********************************************************************************

    Module EXP_Structures_def
    
    TYPE SRM_TYPE
      INTEGER NCELLS
      INTEGER, ALLOCATABLE :: CELLS(:)
      REAL, ALLOCATABLE :: POR(:)
      REAL, ALLOCATABLE :: HGT(:)
      REAL, ALLOCATABLE :: HC(:)
      REAL, ALLOCATABLE :: A(:)
      REAL, ALLOCATABLE :: B(:)                
      REAL BASE
      INTEGER, ALLOCATABLE :: IVAL(:)               
    END TYPE SRM_TYPE
    TYPE (SRM_TYPE) SRM
 
    TYPE SRMU_TYPE
      INTEGER NCELLS
      INTEGER, ALLOCATABLE :: CELLS(:)
      REAL, ALLOCATABLE :: POR(:)
      REAL, ALLOCATABLE :: HGT(:)
      REAL, ALLOCATABLE :: HC(:) 
      REAL, ALLOCATABLE :: L(:) 
      REAL, ALLOCATABLE :: A(:)
      REAL, ALLOCATABLE :: B(:)                        
    END TYPE SRMU_TYPE
    TYPE (SRMU_TYPE) SRMU
      
    TYPE SRMV_TYPE
      INTEGER NCELLS
      INTEGER, ALLOCATABLE :: CELLS(:)
      REAL, ALLOCATABLE :: POR(:)
      REAL, ALLOCATABLE :: HGT(:)
      REAL, ALLOCATABLE :: HC(:)  
      REAL, ALLOCATABLE :: L(:)
      REAL, ALLOCATABLE :: A(:)
      REAL, ALLOCATABLE :: B(:)                       
    END TYPE SRMV_TYPE
    TYPE (SRMV_TYPE) SRMV     
           
    LOGICAL :: STRUCTURES = .false.
    LOGICAL :: SRM_ON = .false.
    
    LOGICAL, ALLOCATABLE :: HASSTRUCT(:)
    REAL, ALLOCATABLE :: HASPOR(:),HASHGT(:),HASHC(:),HASA(:),HASB(:)

    TYPE CUL_TYPE
      INTEGER NUM,NUMUI
      INTEGER, ALLOCATABLE :: CELLS1(:)
      INTEGER, ALLOCATABLE :: CELLS2(:)       
      REAL, ALLOCATABLE :: DIA(:)
      REAL, ALLOCATABLE :: INVERT(:)
      REAL, ALLOCATABLE :: FRIC(:)
      REAL, ALLOCATABLE :: LENGTH(:) 
      REAL, ALLOCATABLE :: FLOW(:) 
      REAL, ALLOCATABLE :: saltF(:)     
      REAL, ALLOCATABLE :: cellUI(:)                              
    END TYPE CUL_TYPE
    TYPE (CUL_TYPE) CUL
      
    LOGICAL :: CUL_ON = .false.        
    
    end Module EXP_Structures_def
    
    
    
!********************************************************************************
MODULE EXP_bndcond_def

TYPE HMOD_TYPE 
  INTEGER INC
  REAL, ALLOCATABLE :: DTIME(:)
  REAL, ALLOCATABLE :: DVALUE(:)
END TYPE HMOD_TYPE
TYPE (HMOD_TYPE) :: HMOD

logical :: WABC = .false.
TYPE WABC_type
  CHARACTER*10 NAME
  real, allocatable :: Q(:),ETA(:),del(:),wdth(:),QN(:)
  integer, allocatable :: SGN(:),CELL(:),RAD(:)
END TYPE WABC_type
TYPE (WABC_type), ALLOCATABLE :: SWABC(:),MWABC(:),TWABC(:)

TYPE QstringEXP_type
integer sgn
logical vface
END TYPE QstringEXP_type
TYPE (QstringEXP_type),  ALLOCATABLE :: QstringEXP(:)

TYPE MVstringEXP_type
logical, allocatable :: cellEW(:)
integer, allocatable :: cell(:)
END TYPE MVstringEXP_type
TYPE (MVstringEXP_type),  ALLOCATABLE :: MVstringEXP(:)

allocatable :: ext_w(:,:),ext_n(:,:),ext_s(:,:),ext_e(:,:)
integer ext_s,ext_n,ext_e,ext_w

logical modify_h

END MODULE
!********************************************************************************

!********************************************************************************
MODULE EXP_transport_def
use prec_def

real*8 SLPFAC,RHOWDIV8,XKS,A0DIVRHOWGRAV,TCR,TWOPI,RHOWDIV2  

real :: rate_avalanche = 0.01

TYPE SALT_type
  REAL diffC
  double precision vol,qx,qy
END TYPE SALT_type	
TYPE (SALT_type), ALLOCATABLE :: SALT(:)	
	
LOGICAL :: saltsimD = .true.              

	
!suspended sediment concentration variables
TYPE ADSS_type
  REAL conc,concn,eros,depo,vol,diffC,qx,qy
END TYPE ADSS_type	
TYPE (ADSS_type), ALLOCATABLE :: ADSS(:)	
double precision voln

TYPE COHES_type
  REAL conc,concn,eros,depo,diffC,ws,T_fluff,C_fluff,tbmax
  real Tcrit_E,Tcrit_D
  double precision vol,qx,qy
END TYPE COHES_type	
TYPE (COHES_type), ALLOCATABLE :: COHES(:)	!mapping arra
	
TYPE CHParms_type
  REAL E,Tcrit_E,Tcrit_D,ws_max,c_max,c_peak,wse_bc,Dfac
  logical :: Tcrit_variable = .false.
  integer numDEPTHS
  real, allocatable :: Depth(:),Tcrit_Ev(:),Tcrit_Dv(:)
END TYPE CHParms_type	
TYPE (CHParms_type), SAVE :: CHparms

logical :: cohesive=.FALSE.,cohesive_read=.FALSE.,cohes_flow_bc=.false.
logical :: cohes_cards_read = .false., SMAG
integer csunit(2),cohes_flow_unit,COHES_UNIT
character*200 csfils(2),cohes_flow_bc_file
real(ikind), allocatable :: cohes_bc(:,:),cohes_bc_time(:),cohes_bc_cells(:,:)

real(ikind), allocatable :: qsx(:),qsy(:),bed(:) !Alex changed to double

REAL    TMORPH_ELAPSE,TSED_ELAPSE,TSALT_ELAPSE


END MODULE 

!********************************************************************************
MODULE NupdateMod
      integer NupdateInt,NupdateCnt
      logical Nupdate
    END MODULE
    
MODULE EXP_TELESCOPING
    use prec_def
    integer, allocatable :: cellfaces(:,:),cellmap(:,:)
    
    integer numxfaces
    integer, allocatable :: xface_CadvF(:,:),xface_cells(:,:),xface_advF(:,:)
    real(ikind), allocatable :: xface_q(:),xface_qn(:),xface_Length(:),xface_gcoef(:,:),xface_vel(:),xface_advdif_I(:,:),xface_advdif_C(:),xface_flux(:)
    logical, allocatable :: xface_wall(:),xface_grad1(:),xface_basic_orientation(:),xface_wet(:)
    
    integer numyfaces
    integer, allocatable :: yface_CadvF(:,:),yface_cells(:,:),yface_advF(:,:)
    real(ikind), allocatable :: yface_q(:),yface_qn(:),yface_Length(:),yface_gcoef(:,:),yface_vel(:),yface_advdif_I(:,:),yface_advdif_C(:),yface_flux(:)
    logical, allocatable :: yface_wall(:),yface_grad1(:),yface_basic_orientation(:),yface_wet(:)   
       
    real(ikind), allocatable :: xTransQ(:),yTransQ(:),xSedTransQ(:),ySedTransQ(:)
    
    integer, allocatable :: normalxfaces(:),normalyfaces(:)
    integer, allocatable :: teleXfaces(:),teleYfaces(:)
    integer Nnormalxfaces,NNormalyfaces,NteleXfaces,NteleYfaces
    
    integer numREGcells,numTBcells,numREgXfaces,numTBXfaces,numREGYfaces,numTBYfaces
    integer, allocatable :: REGcells(:),TBcells(:),REgXfaces(:),TBXfaces(:),REGYfaces(:),TBYfaces(:)
    
    integer, allocatable :: specX(:),specY(:)
    integer numspecX,numspecY
    
    !these arrays are used in Zupte_advection_tel routine (v00 and higher)
    !to deal with issue identified at boundary leading to instabilities, but
    !was also causing smaller errors at all TB faces
    integer, allocatable :: xface_side(:),yface_side(:)
    
    !to expedite the calculation of velocities used in coreolis force calculation
    !they are calcualted in advection term routine to minimze index lookups
    real, allocatable :: QXc(:),QYc(:)

END MODULE