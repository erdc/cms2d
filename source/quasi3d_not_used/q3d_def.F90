!======================================================================
module q3d_def
! Quasi-3D variable definitions
!======================================================================
    use prec_def
    implicit none
    save
    
    !Turns on and off the quasi3d module
    logical:: q3d
    logical:: q3d_to_flow !Turns on the source term in the flow
    logical:: q3d_to_sed  !Turns on the source term in the sediment
    
    !Wave-current interaction term
    logical:: wavcurint
    real(ikind):: facwci
    real(ikind), allocatable :: udsx(:),udsy(:)
    
    !Velocity for wave model
    integer :: ivelwav
    character(len=20) :: avelwav(0:3)
    data avelwav /'NONE',&     !0 - None
                  'MEAN',&     !1 - Mean = depth-averaged
                  'SURFACE',&  !2 - uc(zp=1) = U + uds
                  'WEIGHTED'/  !3 - U + 2*kw/sinh(2*kw*h)*int(ud*cosh(2*kw*h*zp),0,1)
    
    !Streamwise profile due to bottom friction and wind (outside of surf zone)
    integer :: istreampro
    character(len=20) :: astreampro(0:2)
    data astreampro /'NONE',& !0 - None
                     'LOG',&  !1 - Logarithmic
                     'POWER'/ !2 - Powerlaw
    
    !Transverse profile (normal) due to helical flow and Coriolis
    integer :: itranspro
    character(len=20) :: atranspro(0:1)
    data atranspro /'NONE',&  !0 - None
                    'LINEAR'/ !1 - Linear
    
    !Profile in the surf zone
    integer :: isurfpro
    character(len=20) :: asurfpro(0:2)
    data asurfpro /'NONE',& !0 - None
                   'QUAD',& !1 - Quadratic
                   'LOG'/   !2 - Logarithmic
    
    !Vertical eddy viscosity
    real(ikind) :: visvcon   !Base or constant value for vertical eddy viscosity    
    real(ikind) :: cvisvhc   !Coefficient for Svendsen-type vertical viscosity
    real(ikind) :: cvisvv    !Coefficient for friction velocity contributions to vertical eddy viscosity
    real(ikind) :: verteff   !Efficiency factor for surface shear stress based on wave dissipation 
    real(ikind) :: visvslp0  !Base Slope factor for linear vertical eddy viscosity
    real(ikind) :: cvisvslpmu  !Coefficient for slope
    !real(ikind) :: visvm !Mean vertical eddy viscosity (temporary variable)
    integer :: ieddyvert
    character(len=20) :: aeddyvert(1:3)
    data aeddyvert /'CONSTANT',&          !1 - Constant value
                    'HORIZONTAL_SCALED',& !2 - Only for testing
                    'PARABOLIC'/          !3 - Parabolic
    
    !Flow Dispersion terms
    real(ikind),allocatable:: f3dxx(:),f3dxy(:),f3dyy(:)  !Hydrodynamic 3D current-current (dispersion) and current-wave interaction terms
    real(ikind),allocatable:: f3du(:),f3dv(:)

    !Sediment Dispersion terms
    real(ikind),allocatable:: s3dxk(:,:),s3dyk(:,:)       !Sediment 3D dispersion terms
    
    !Output Stations
    integer:: nzpsta !Number of vertical points
    real(ikind),allocatable:: zpsta(:),uzpsta(:),vzpsta(:)  !For stations
    
    !Output Layers
    logical :: q3d_lay
    integer:: nzplay
    real(ikind),allocatable:: zplay(:),uzplay(:,:),vzplay(:,:)  !For stations
    
end module q3d_def