!========================================================    
module der_def
! Spatial Derivative Operator Definitions
! Author: Alex Sanchez, USACE-CHL
!=====================================================
    use prec_def, only: ikind
    implicit none
    save  
    
    !--- Spatial Derivatives ---------------------
    integer :: nder  !Method for spatial derivative calculation
    integer :: npow  !Inverse distance power
    character(len=7) :: ader(0:8)
    data ader /'HYBRID',&  !0, Hybrid method based on local stencil
               'CBGG',&    !1, Cell-based Green-Gauss (1st order for skewed cells)
               'CBGGCR',&  !2, Cell-based Green-Gauss with cell-recontruction (second order)
               'CBGGFR',&  !3, Cell-based Green-Gauss with face-recontruction (second order)
               'CBGGLR',&  !4, Cell-based Green-Gauss with linear interpolation (second order)
               'NBGG',&    !5, Nodal-based Green-Gauss (first order)
               'CBWLSFS',& !6, Cell-based least-squares using face-sharing cells (compact stensil, first order)
               'CBWLSNS',&  !7, Cell-based least-squares using node-sharing cells (extended stensil, first order)
               'CBFD'/     !8, Cell-based Finite-Difference (second order)
    
    !---- Gradient Operator --------------------
    type der_go_type
      integer :: ibc  !Boundary condition
      integer :: ider !Derivative scheme
      integer :: nd   !Maximum size of cells
      integer :: nsc  !Maximum size of stencil for cell centroids
      integer :: nsg  !Maximum size of stencil for cell gradients
      integer,    pointer :: ncx(:),ncy(:)     !Number of neighboring cell centoids in x and y directions (cell)
      integer,    pointer :: icx(:,:),icy(:,:) !Index of neighboring cells centoids in x and y directions (stencil,cell)
      real(ikind),pointer :: wcx(:,:),wcy(:,:) !Weights for cell centroid values (stencil,cell)
      integer,    pointer :: ngx(:),ngy(:)       !Number of neighboring cell gradients used for reconstructions (cell)
      integer,    pointer :: igx(:,:),igy(:,:)   !Index of neighbroing cell gradients used for reconstrctions (stencil,cell)
      real(ikind),pointer :: wgxy(:,:),wgyx(:,:) !Weights for cell gradient values (stencil,cell)
      real(ikind),pointer :: wgxx(:,:),wgyy(:,:) !Weights for cell gradient values (stencil,cell)
    endtype der_go_type
    
    !List not needed for Static Operator since it is computed once for whole grid and doesn't change
    type(der_go_type)  :: goa
    
    !List needed for Dynamic Operator because it changes as a function of wetting and drying
    integer              :: nlistcw
    integer, allocatable :: ilistcw(:)
    integer              :: nlistgw
    integer, allocatable :: ilistgw(:)
    type(der_go_type)  :: gow
    
    !!do i=1,ncells
    !!  dvarx(i)=0.0
    !!  do k=1:derop%ncx(i)
    !!    icx=derop%icx(k,i)
    !!    wck=derop%wcx(k,i)
    !!    dvarx(i)=dvarx(i)+wck*var(icx)
    !!  enddo
    !!enddo
    
    !integer :: nsc !Maximum size of stencil for cell centroids
    !integer :: nsg !Maximum size of stencil for cell gradients
    !
    !!Cell centroid index and weights
    !integer                :: nlistca
    !integer,    allocatable:: ilistca(:)
    !integer,    allocatable:: ncxa(:),  ncya(:)   !Number of neighboring cell centoids in x and y directions
    !integer,    allocatable:: icxa(:,:),icya(:,:) !Index of neighboring cells centoids in x and y directions
    !real(ikind),allocatable:: wcxa(:,:),wcya(:,:) !Weights for cell centroid values
    !
    !!Cell gradient index and weights used for cell-based or face-based reconstructions
    !integer                :: nlistga
    !integer,    allocatable:: ilistga(:)
    !integer,    allocatable:: ngxa(:),  ngya(:)     !Number of neighboring cell gradients used for reconstructions
    !integer,    allocatable:: igxa(:,:),igya(:,:)   !Index of neighbroing cell gradients used for reconstrctions
    !real(ikind),allocatable:: wgxya(:,:),wgyxa(:,:) !Weights for cell gradient values
    !real(ikind),allocatable:: wgxxa(:,:),wgyya(:,:) !Weights for cell gradient values    
    !
    !!--- Dynamic Derivative Operator -------------------
    !!Cell centroid index and weights
    !integer                :: nlistcw
    !integer,    allocatable:: ilistcw(:)
    !integer,    allocatable:: ncxw(:),ncyw(:)     !Number of neighboring cell centoids in x and y directions
    !integer,    allocatable:: icxw(:,:),icyw(:,:) !Index of neighboring cells centoids in x and y directions
    !real(ikind),allocatable:: wcxw(:,:),wcyw(:,:) !Weights for cell centroid values    
    !
    !!Cell gradient index and weights used for cell-based or face-based reconstructions
    !integer                :: nlistgw,nsgw
    !integer,    allocatable:: ilistgw(:)
    !integer,    allocatable:: ngxw(:),ngyw(:)       !Number of neighboring cell gradients used for reconstructions
    !integer,    allocatable:: igxw(:,:),igyw(:,:)   !Index of neighbroing cell gradients used for reconstrctions
    !real(ikind),allocatable:: wgxyw(:,:),wgyxw(:,:) !Weights for cell gradient values
    !real(ikind),allocatable:: wgxxw(:,:),wgyyw(:,:) !Weights for cell gradient values    

    !--- Slope Limiters ------------------------------------------------
    integer :: nlim
    character(len=10) :: alim(0:8)
    data alim /'NONE',& !0 
               'MINMOD',& !1
               'VAN_LEER',& !2
               'VAN_ALBADA',& !3
               'MUSCL',& !4
               'BJ',& !5
               'BJV',& !6
               'LCD',& !7
               'LCDV'/ !8
    real(ikind),allocatable:: rgx(:),rgy(:)   !Grid ratio (expansion) used for slope limiters
    
end module der_def