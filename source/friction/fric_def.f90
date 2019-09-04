!==================================================
module fric_def
! Friction variable definitions module
! written by Alex Sanchez, USACE-CHL
!==================================================    
    use prec_def
    implicit none
    save
    
    !Bottom roughness specification
    integer:: mbedfric                     !Input bottom friction type
    character(len=32) :: abedfric(0:3)     !Bottom friction input type
    data abedfric /'DYNAMIC',&             !0, Dynamic 
                   'DRAG_COEF',&           !1, drag coefficient Cd
                   'MANNING_COEF',&        !2, Manning's n coefficient
                   'ROUGHNESS_HEIGHT'/     !3, roughness height z0
    real(ikind), allocatable :: coefman(:) !Manning's coefficient [s/m^(1/3)]
    real(ikind), allocatable :: cfrict(:)  !Friction coefficient [-]
    real(ikind), allocatable :: z0(:)      !Roughness length [m]
    real(ikind), allocatable :: rksbase(:) !Base roughness height representing rocks and obstacles [m]
    real(ikind) :: fricscale  !Scaling factor for bottom roughness (useful for calibration)
    character(len=200) :: fricfile,fricpath
    
    !Wave Bottom Boundary Streaming (alpha)
    logical :: bbl_stream
    integer:: istreamrough !1-Bed roughness, 2-Apparent Roughness
    
    !Bottom friction
    integer:: mwavcurint        !Formulation for wave-current bottom shear stress
    real(ikind):: cfricwav      !Wave friction coefficient for quadratic formula     
    logical:: constbotfric
    real(ikind):: cbotfric      !Constant bottom friction value, either n, cb, or z0
	character(len=32) :: awavcur(7)  !Wave-current bed shear stress
    data awavcur/'QUAD',&   !1, Quadratic
                 'DATA2',&  !2, Soulsby (1995)
                 'DATA13',& !3, Soulsby (1995)
                 'HT91',&   !4, Huynh-Thanh and Temperville (1991)
                 'F84',&    !5, Fredsoe (1984)
                 'DSK88',&  !6, Davies et al. (1988)
                 'GM79'/    !7, Grant and Madsen (1979)
    !Generalized wave-current bottom shear stress formulation of Soulsby (1995)
    real(ikind) :: logfw(4),logCd(4),ciwc,cjwc
    real(ikind) :: awc(4),cmwc(4),cnwc(4),bwc(4),pwc(4),qwc(4)        
    
    real(ikind), allocatable :: uelwc(:)     !Mean velocity including wave enhacement factor [-]
    real(ikind), allocatable :: bsxy(:)      !Bed shear stress [N/m^2]
    real(ikind), allocatable :: bsvel(:)     !Bed shear velocity [m/s]
    real(ikind), allocatable :: cbcfuwcap(:) !Friction variables for implicit scheme
    real(ikind), allocatable :: cbcfuwc(:)   !Friction variables for explicit scheme     
    
    !Wall friction
    logical :: wallfric
    !integer :: nwallfric,iwallfric(2)
    real(ikind) :: wallfac !,wallman
    
    !Bed slope friction factor
    logical :: fricbedslope
    real(ikind), allocatable :: cmb(:) !Bed-slope coefficient [-]
    
    !Dynamic roughness (mbedfric==0)
    real(ikind), allocatable :: rksr(:),rksmr(:),rksd(:),rksg(:),rkstc(:)
    integer :: mroughripple     !0-none, 1-Soulsby (1997), 2-van Rijn (2007)
    integer :: mroughmegaripple !0-none, 1-Sanchez (2014), 2-van Rijn (2007)
    integer :: mroughdune       !0-none, 1-van Rijn (1984), 2-van Rijn (2007)
    integer :: mripple          !0-None, 1-Soulsby (1997), 2-Raudkivi (1998), 3-Soulsby and Whitehouse (2005)
    integer :: mripplewave    !Method for wave-generated ripples, 0-none, 1-van Rijn (1984), 2-Soulsby and Whitehouse (2005)
    integer :: mripplecurrent !Method for current-generated rippples, 0-none, 1-Soulsby (1997), 2-Soulsby and Whitehouse (2005), 3-Raudkivi (1998)
    integer :: mdune          !Method for dunes (from currents), 0-none, 1-van Rijn
    integer :: mmegaripple    !Method for megaripples (from currents), 0-none, 1-van Rijn
    integer :: mroughtranscur !Method for transport roughness
    real(ikind) :: roughscaleripple        !Roughness scaling facctor for ripples
    real(ikind) :: roughscalemegaripple    !Roughness scaling facctor for megaripples
    real(ikind) :: roughscaledune          !Roughness scaling facctor for current-related rippples
    real(ikind) :: roughscalegrain         !Roughness scaling factor for grains
    real(ikind) :: roughscaletrans         !Roughness scaling factor for sediment transport
    real(ikind), allocatable :: riplen(:)  !Ripple length [m]
    real(ikind), allocatable :: riphgt(:)  !Ripple height [m]
    real(ikind), allocatable :: ripdir(:)  !Ripple direction/orientation [rad] 0-East, pi/2-North
    real(ikind), allocatable :: dunelen(:)  !Dune length [m]
    real(ikind), allocatable :: dunehgt(:)  !Dune height [m]
    real(ikind), allocatable :: megariplen(:) !Megaripple length [m]
    real(ikind), allocatable :: megariphgt(:) !Megaripple height [m]
    real(ikind) :: biodegradhalflife  !Biodegradation half-life
    
endmodule fric_def

