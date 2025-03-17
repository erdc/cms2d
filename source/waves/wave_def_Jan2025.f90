!======================================================================
module global_inline
! CMS-Wave 
!======================================================================
    implicit none
    save
    
    real, parameter    :: A=0.17  !A=0.15 suggested by Li's experiments
    integer :: npf,mpd,mpd2,ipmx,jpmx    !changed for dynamic allocation, Wu, Nov. 2024.
    integer :: igpx,jgpx,mpmx
    integer :: ijpmx,ijgpx
    integer :: komx,nomx
    integer :: koutnest,knest
    
    !added a few other changable parameters to share  MEB  10/19/2021
    real    :: gamma_bj78 = -1.0
    logical :: suppress_obs = .false.  !added 12/06/2021 MEB
    
end module global_inline

    
!======================================================================
module wave_def
! CMS-Wave variables                              !Wu, Nov. 2024
!======================================================================
    use prec_def     
    implicit none
    save
 
    integer       :: igetfile20,nship
    integer       :: ijstruc1,ijstruc2,ijstruc3,ijstruc4,ismall
    integer       :: imod,iprp,island,imd,iprpp,nonln,igrav,isolv,ixmdf
    integer       :: iproc,imud,iwnd
    integer       :: NF,MD,IMAX,JMAX,IGMX,JGMX,JCPB,JCPE,JCB,JCE,NFF,MDD  
    integer       :: ICK3,ICK4,ICHOICE
    integer       :: IBND,KPEAK,IBACK,NBLOCK,IWIND     
    integer       :: kdate,idate
    integer       :: IBK
    integer       :: KRMX
    integer       :: KRF
    integer       :: KOUT,IRS,IBREAK,ICUR,IWET,INST
    integer       :: IWVBK,ibk3
    integer       :: nest
    integer       :: isteer,iidate
    integer       :: iview,iplane,iwave,irunup,iroll
    integer       :: ni,nj                                                   
    integer       :: itms,ibf,iark,iarkr                             
    integer       :: inest1,nestin,nestin1,nestin2
    real(ikind)   :: PAI2,PAI,HPAI,RAD,akap,depmin0
    real(ikind)   :: HS0,wd,ws,ph0
    real(ikind)   :: T13MIN
    real(ikind)   :: TP,PRD,VL,TIDE,WSMAG
    real(ikind)   :: DX,DY,DXX,DMESH,DTH,depmax,depmax0
    real(ikind)   :: DBAR,WL0
    real(ikind)   :: x0,y0,azimuth,sinaz,cosaz
    real(ikind)   :: depmin,g,cflat
    real(ikind)   :: bf,ark,arkr                                     
    real(ikind)   :: hsmin,hs13n,azimnest                            
    character*180 :: OptsFile,  DepFile,  CurrFile,  EngInFile, WaveFile, ObsFile,   EngOutFile, NestFile
    character*180 :: BreakFile, RadsFile, StrucFile, SurgeFile, MudFile,  FricFile,  FrflFile, BrflFile
    character*180 :: SpecFile,  WindFile, XMDFFile,  SetupFile, SeaFile,  SwellFile, ShipFile          
    logical :: getfile4,getfile5,getfile7,getfile8           !Wu/Zhang 8/1/2009  

    real(ikind),allocatable :: dep0(:,:),fsp(:),xc(:),yc(:),wc(:),wcd(:)
    real(ikind),allocatable :: hs13(:),refltx(:,:),reflty(:,:),exx(:,:),eyy(:,:)
    real(ikind),allocatable :: dvarxx(:),dvaryy(:)
    real(ikind),allocatable :: sxx(:,:),sxy(:,:),syy(:,:)
    real(ikind),allocatable :: wxrs(:,:),wyrs(:,:)
    real(ikind),allocatable :: sxxx(:,:),sxyx(:,:),sxyy(:,:),syyy(:,:)
    real(ikind),allocatable :: cosa(:),sina(:),d1(:,:),wk2(:,:,:),cgp(:,:)
    real(ikind),allocatable :: disx(:),disy(:)
    real(ikind),allocatable :: depin(:,:),etain(:,:),uin(:,:),vin(:,:)  
    real(ikind),allocatable :: ShipL(:),ShipB(:),ShipD(:),ShipS(:)
    real(ikind),allocatable :: dstruc1(:),dstruc2(:),dstruc3(:)
    real(ikind),allocatable :: dstruc33(:),dstruc4(:),dstruc44(:)
    real(ikind),allocatable :: sw13(:,:),sa13(:,:),tw13(:,:)
    real(ikind),allocatable :: ta13(:,:),dw13(:,:),da13(:,:)
    real(ikind),allocatable :: wdd(:),aslop(:)
    real(ikind),allocatable :: WCC(:),HSB(:),HSG(:),DCD(:)  
    real(ikind),allocatable :: DVARX(:),DVARY(:),ETA(:,:),HSK(:)
    real(ikind),allocatable :: DEP(:,:),DEPS(:),DBIG(:)
    real(ikind),allocatable :: RK(:),yangl(:)    
    real(ikind),allocatable :: RKR(:),xangl(:)    
    real(ikind),allocatable :: DEPM(:),DMNJ(:),SLF(:),wlmn(:),cmn(:),sigm(:)
    real(ikind),allocatable :: H13(:,:),T13(:,:),DMN(:,:)
    real(ikind),allocatable :: H13R(:,:),T13R(:,:),DMNR(:,:)
    real(ikind),allocatable :: H13S(:,:),DISS(:,:)
    real(ikind),allocatable :: H13F(:,:),DMNF(:,:),T13F(:,:)
    real(ikind),allocatable :: FCN(:),DCM(:),DSFD(:,:)
    real(ikind),allocatable :: DF(:),FFCN(:),DFINP(:),PL1E(:),PL1T(:)
    real(ikind),allocatable :: u(:,:),v(:,:)   
    real(ikind),allocatable :: u1(:,:),v1(:,:)    
    real(ikind),allocatable :: u10(:,:),v10(:,:)
    real(ikind),allocatable :: bfric(:,:),amud(:,:)      
    real(ikind),allocatable :: SCP(:,:),SI(:,:,:)      
    real(ikind),allocatable :: SIMAX(:,:,:),SJJ(:)     
    real(ikind),allocatable :: SOP(:,:,:),SR(:,:,:)     
    real(ikind),allocatable :: SJ(:),SJF(:,:),FJF(:)         
    real(ikind),allocatable :: AA(:,:),B(:),X(:)        
    real(ikind),allocatable :: sgma0(:,:),sgma1(:,:)       
    real(ikind),allocatable :: cwk(:,:,:),cgk(:,:,:)        
    real(ikind),allocatable :: ex(:,:),ey(:,:)        

    integer,allocatable :: IJB(:,:)
    integer,allocatable :: KR(:,:)
    integer,allocatable :: KCR(:,:)
    integer,allocatable :: IBR(:,:)
    integer,allocatable :: istruc1(:),jstruc1(:),istruc2(:),jstruc2(:)
    integer,allocatable :: istruc3(:),jstruc3(:),k3(:)
    integer,allocatable :: istruc4(:),jstruc4(:),kstruc4(:),k4(:)
    integer,allocatable :: IJSP(:,:)
    integer,allocatable :: inest(:),jnest(:),ix1(:),ix2(:)
    integer,allocatable :: IA(:,:)   

end module wave_def
