!==========================================================================
module rbf_lib
! Radial basis function interpolation library
!
! Uses the following equation for interpolation
!   fi(xi) = c0 + sum(c1*xi) + sum(cc*phi(||x-xi||))
! where
!   c0     - constant
!   c1(:)  - linear coefficients for each dimension (ndim)
!   cc(:)  - expansion coefficients (npts)
!   x(:,:) - data coordinates with values (ndim,npts)
!   xi     - interpolation coordinate
!   fi     - interpolated value
!   phi(r) - radial basis function
!   r      - radial distance (r>=0)
!   || ||  - Euclidean norm
!  
! The interpolation conditions lead to a symmetric linear sytem
!                A               *   coef    =    b
!   [ aa(:,:) | 1(:)  | x(:,:) ]   [ cc(:) ]   [ f(:) ] 
!   [  1(:)'  |  0    | 0(:)'  ] * [ c0    ] = [  0   ]
!   [ x(:,:)' | 0(:)' | 0(:,:) ]   [ c1(:) ]   [ 0(:) ]
! where
!   aa - symmetric matrix of radial basis function values 
!      (i.e. aa(i,j) = phi(||x(:,i)-x(:,j)||))
!   f - data values corresponding to x
!
! Contains the following subroutines
!   rbf_create(rbfvar,ndim,npts)
!     Creates an rbf interpolation variable
!     Allocates variables and set defaults
!     Should be called once for each rbf variable.
!
!   rbf_setup(rbfvar,x,fun,constant,smooth)
!     Override defaults, and sets up the RBF matrix (assembly)
!     and does a Lower-Upper decomposition by Crout's algorithm
!     with partial pivoting. The subroutine should be called everytime
!     the coordinates of the function change (e.g. adapting mesh).
!
!   rbf_solve(rbfvar,f)
!     Fits an RBF interpolant to the function f.
!     Solves the RBF matrix and computes the interolation 
!     coefficients by LU subsitution.
!     Should be called every time the function values change.
!
!   rbf_interp(rbfvar,ni,xi,fi)
!     Interpolates points using the RBF variables
!
!   rbf_destroy(rbfvar)
!     Destroys an RBF interpolation variable
!
! Usage:
!   call rbf_create(rbfvar,ndim,npts)
!   call rbf_setup(rbfvar,x,fun='multiquadric',constant=2.0)
!   call rbf_solve(rbfvar,f)
!   call rbf_interp(rbfvar,ni,xi,fi)
!   call rbf_destroy(rbfvar)
!
! Author: Alex Sanchez, USACE-CHL
!==========================================================================
    use rbf_def
    use prec_def
    implicit none
    
    private    
    public :: rbf_create,rbf_setup,rbf_solve,rbf_interp,rbf_destroy
    
  contains
!*********************************************************************
  subroutine rbf_create(rbfvar,ndim,npts)
! Creates an rbf interpolation variable.
! Allocates variables and set defaults.
! Should be called once for each rbf variable.
!
! Arguments:
!   rbfvar - radial basis function variable
!   ndim   - # of dimensions
!   npts   - # of data points
! 
! Author: Alex Sanchez, USACE-CHL    
!*********************************************************************
    implicit none
    !Input/Output
    type(rbf_type),intent(inout) :: rbfvar
    integer,intent(in) :: ndim,npts
    !Internal Variables
    integer :: n
    
    n=npts+ndim+1
    rbfvar%ndim = ndim !# of dimensions
    rbfvar%npts = npts !# of data points
    rbfvar%nsize = n   !Order of system
    rbfvar%fun = 'multiquadric' !Original function proposed by Hardy, R.L. (1968)
    rbfvar%const = -999.0 !Undefined at this point. should be approximate average distance between the points 
    rbfvar%smooth = 0.0   !No smoothing    
    
    !Allocate variables    
    allocate(rbfvar%x(ndim,npts)) !Data coordinates
    allocate(rbfvar%A(n,n))       !LU decomposition of rbf matrix
    allocate(rbfvar%coef(n))      !Coefficients (rbf linear system solution)
    allocate(rbfvar%ind(n))       !Index for LU decomposition
    
    return
  endsubroutine rbf_create
    
!*********************************************************************
  subroutine rbf_setup(rbfvar,x,fun,constant,smooth)
! Description:
!   Override defaults, and sets up the RBF matrix (assembly)
!   and does a Lower-Upper decomposition by Crout's algorithm
!   with partial pivoting. The subroutine should be called everytime
!   the coordinates of the function change (e.g. adapting mesh).
!
! Arguments:
!   rbfvar       - Radial basis function variable
!   x(ndim,npts) - Data coordinates
!   fun          - Radial basis function
!   constant     - Radial basis function constant/parameter
!   smooth       - smoothing coefficient
!
! Author: Alex Sanchez, USACE-CHL    
!*********************************************************************
    use LU_lib
    implicit none
    !Input/Output
    type(rbf_type),intent(inout) :: rbfvar
    real(ikind),intent(in) :: x(rbfvar%ndim,rbfvar%npts)
    character(len=*),intent(in),optional :: fun
    real(ikind),intent(in),optional :: constant,smooth
    !Internal Variables
    integer :: i,ierr
    
    !Override defaults
    if(present(fun))then
      rbfvar%fun = fun
      call lowercase(rbfvar%fun)
    endif
    if(present(constant)) rbfvar%const = constant
    if(present(smooth))   rbfvar%smooth = smooth
    
    !Check constant
    selectcase(rbfvar%fun)
    case('multiquadric','invquadratic','invmultiquadric','gaussian')
      !if not specified then estimate  
      if(rbfvar%const<1.0e-5)then !Automatic
        !approximate average distance between the points
        rbfvar%const = 0.0
        do i=1,rbfvar%ndim
          rbfvar%const = rbfvar%const*(maxval(x(i,:))-minval(x(i,:)))
        enddo
        rbfvar%const = (rbfvar%const/rbfvar%npts)**(1.0/real(rbfvar%ndim))
      endif
    case default;
      rbfvar%const = 0.0 !Not needed for other functions  
    endselect
    
    !Save coordinates
    !if(.not.allocated(rbfvar%x))then
    !  write(*,*) 'ERROR: Must call rbf_create before rbf_setup'
    !  read(*,*)
    !  stop  
    !endif
    rbfvar%x = x
    
    !Coefficient matrix
    selectcase(rbfvar%fun)
    case('linear');          call rbf_assemble(rbfvar,rbf_linear)
    case('cubic');           call rbf_assemble(rbfvar,rbf_cubic)
    case('thinplate');       call rbf_assemble(rbfvar,rbf_thinplate)
    case('thinplatespline'); call rbf_assemble(rbfvar,rbf_thinplatespline)
    case('multiquadric');    call rbf_assemble(rbfvar,rbf_multiquadric)
    case('invmultiquadric'); call rbf_assemble(rbfvar,rbf_invmultiquadric)
    case('invquadratic');    call rbf_assemble(rbfvar,rbf_invquadratic)
    case('gaussian');        call rbf_assemble(rbfvar,rbf_gaussian)
    case default;            call rbf_assemble(rbfvar,rbf_linear)
    endselect
    
    !Decompositions
    call LU_decomp(rbfvar%nsize,rbfvar%A,rbfvar%ind,ierr)
    
    return
  endsubroutine rbf_setup

!*******************************************************************    
  subroutine rbf_assemble(rbfvar,rbf)
! Assembles the coefficient matrix for the rbf interpolation
!
! Builds the matrix x for the following symmetric linear system
!                A               *   coef    =    b
!   [ aa(:,:) | 1(:)  | x(:,:) ]   [ cc(:) ]   [ f(:) ] 
!   [  1(:)'  |  0    | 0(:)'  ] * [ c0    ] = [  0   ]
!   [ x(:,:)' | 0(:)' | 0(:,:) ]   [ c1(:) ]   [ 0(:) ]
! where
!   aa(:,:) - symmetric matrix of radial basis function values 
!             (i.e. aa(i,j) = phi(||x(:,i)-x(:,j)||))
!   r       - radial distance (r>=0)
!   f(:)    - data values corresponding to x
!   c0      - constant
!   c1(:)   - linear coefficients for each dimension (ndim)
!   cc(:)   - expansion coefficients (npts)
!   x(:,:)  - data coordinates with values (ndim,npts)
!
! Author: Alex Sanchez, USACE-CHL
!*******************************************************************
    implicit none
    !Input/Output
    type(rbf_type),intent(inout) :: rbfvar
    !Internal Variables
    integer :: i,j,nptsp1
    real(ikind) :: r
    
    interface
      function rbf(r,c) result(phi)
        use prec_def
        real(ikind),intent(in) :: r,c
        real(ikind) :: phi
      endfunction
    endinterface
    
    nptsp1 = rbfvar%npts+1
!$OMP PARALLEL DO PRIVATE(i,j,r) IF(rbfvar%npts>500)
    do j=1,rbfvar%npts
      do i=1,j
        r = sqrt(sum(abs(rbfvar%x(:,i)-rbfvar%x(:,j))**2)) !Euclidian norm (radial distance)
        rbfvar%A(i,j) = rbf(r,rbfvar%const) !Radial basis functions
        rbfvar%A(j,i) = rbfvar%A(i,j) !Symmetry
      enddo !i
      rbfvar%A(j,j) = rbfvar%A(j,j) - rbfvar%smooth !Add smooth to diagonal
      rbfvar%A(j,nptsp1) = 1.0 !Linear part
      rbfvar%A(nptsp1,j) = 1.0 !Linear part 
      do i=1,rbfvar%ndim    
        rbfvar%A(nptsp1+i,j) = rbfvar%x(i,j) !Constant part
        rbfvar%A(j,nptsp1+i) = rbfvar%x(i,j) !Constant part
      enddo !i
    enddo !j    
!$OMP END PARALLEL DO 
    rbfvar%A(nptsp1:rbfvar%nsize,nptsp1:rbfvar%nsize) = 0.0 !Zeros
      
    return
  endsubroutine rbf_assemble    
    
!**********************************************************************************    
  subroutine rbf_solve(rbfvar,f)
! Fits an RBF interpolant to the function f.
!
! Input:
!   f - Data or funtion to be used for interpolant
!
! Solves the RBF matrix and computes the interolation coefficients by LU subsitution.
! Should be called every time the function values.
!
! Author: Alex Sanchez, USACE-CHL
!**********************************************************************************    
    use LU_lib
    implicit none
    !Input/Output
    type(rbf_type),intent(inout) :: rbfvar
    real(ikind),intent(in) :: f(rbfvar%npts)
    
    !Right-hand-side
    rbfvar%coef(1:rbfvar%npts) = f(1:rbfvar%npts)
    rbfvar%coef(rbfvar%npts+1:rbfvar%nsize) = 0.0
    call LU_subs(rbfvar%nsize,rbfvar%A,rbfvar%ind,rbfvar%coef)
    
    return
  endsubroutine
    
!**********************************************************************************    
  subroutine rbf_interp(rbfvar,ni,xi,fi)
! Augmented Radial Basis Function Interpolation
!
! Uses the following equation for interpolation
!   fi(xi) = c0 + c1*xi + sum(lambda(:)*phi(||x-xi||))
! where
!   c0 - Constant
!   c1 - Linear coefficient(s)
!   lambda - radial basis funciton coefficients
!   x  - data coordinates with values (ndim,npts)
!   xi - interpolation coordinate
!   fi - interpolated value
!   phi - radial basis function
!  
! Author: Alex Sanchez, USACE-CHL
!***********************************************************************************
    implicit none
    !Input/Output
    type(rbf_type),intent(inout) :: rbfvar
    integer,intent(in) :: ni               !# of interpolation points
    real(ikind),intent(in) :: xi(rbfvar%ndim,ni)  !Interpolation coordinates
    real(ikind),intent(out) :: fi(ni)      !Interpolation values 
    
    selectcase(rbfvar%fun)
    case('cubic');           call rbf_interp_fun(rbfvar,rbf_cubic,ni,xi,fi)
    case('thinplate');       call rbf_interp_fun(rbfvar,rbf_thinplate,ni,xi,fi)
    case('thinplatespline'); call rbf_interp_fun(rbfvar,rbf_thinplatespline,ni,xi,fi)
    case('multiquadric');    call rbf_interp_fun(rbfvar,rbf_multiquadric,ni,xi,fi)
    case('invquadratic');    call rbf_interp_fun(rbfvar,rbf_invquadratic,ni,xi,fi)
    case('invmultiquadric'); call rbf_interp_fun(rbfvar,rbf_invmultiquadric,ni,xi,fi)    
    case('gaussian');        call rbf_interp_fun(rbfvar,rbf_gaussian,ni,xi,fi)
    case default;            call rbf_interp_fun(rbfvar,rbf_linear,ni,xi,fi)
    endselect

    return
  endsubroutine rbf_interp
    
!***********************************************************************************    
  subroutine rbf_interp_fun(rbfvar,rbf,ni,xi,fi)
! Performs the rbf interpolation using function rbf
!
! Input:
!   fun - Radial basis function
!
! Author: Alex Sanchez, USACE-CHL
!***********************************************************************************
    implicit none
    !Input/Output
    type(rbf_type),intent(inout) :: rbfvar
    integer,intent(in) :: ni !Dimension and # of Points
    real(ikind),intent(in) :: xi(rbfvar%ndim,ni)  !Points
    real,intent(out) :: fi(ni)
    !Internal Variables
    integer :: i,k
    real(ikind) :: r
    
    interface
      function rbf(r,c) result(phi)
        use prec_def
        real(ikind),intent(in) :: r,c
        real(ikind) :: phi
      endfunction
    endinterface
    
!$OMP PARALLEL DO PRIVATE(i,k,r) IF(ni>500)
    do i=1,ni
      fi(i) = rbfvar%coef(rbfvar%npts+1)  !Constant
      do k=1,rbfvar%ndim
        fi(i) = fi(i) + rbfvar%coef(rbfvar%npts+1+k)*xi(k,i) !linear part
      enddo
      do k=1,rbfvar%npts 
        r = sqrt(sum((xi(:,i)-rbfvar%x(:,k))**2)) !Euclidean distance  
        fi(i) = fi(i) + rbfvar%coef(k)*rbf(r,rbfvar%const) !Radia basis functions
      enddo
    enddo
!$OMP END PARALLEL DO

    return
  endsubroutine rbf_interp_fun    
    
!***********************************************************    
  subroutine rbf_destroy(rbfvar)
! Destroys an rbf variables
! Author: Alex Sanchez, USACE-CHL
!***********************************************************
    implicit none
    type(rbf_type) :: rbfvar
    
    deallocate(rbfvar%x,rbfvar%coef)
    deallocate(rbfvar%A,rbfvar%ind)
    
    rbfvar%smooth = 0.0
    rbfvar%const = 0.0
    rbfvar%ndim = 0
    rbfvar%npts = 0
    rbfvar%nsize = 0
    
    return
  endsubroutine rbf_destroy
    
!***********************************************
  function rbf_linear(r,c) result(u)
! Linear Radial Basis Function
!***********************************************
    implicit none
    real(ikind),intent(in) :: r,c
    real(ikind) :: u
    u = r
  endfunction rbf_linear
    
!***********************************************
  function rbf_cubic(r,c) result(phi)
! Cubic Radial Basis Function
!***********************************************
    implicit none
    real(ikind),intent(in) :: r,c
    real(ikind) :: phi
    phi = r*r*r
  endfunction rbf_cubic
    
!***********************************************
  function rbf_gaussian(r,c) result(phi)
! Gaussian Radial Basis Function
!***********************************************
    implicit none
    real(ikind),intent(in) :: r,c
    real(ikind) :: phi
    phi = exp(-0.5*(r/c)**2)
  endfunction rbf_gaussian
    
!***********************************************
  function rbf_multiquadric(r,c) result(phi)
! multiquadric Radial Basis Function    
!***********************************************
    implicit none
    real(ikind),intent(in) :: r,c
    real(ikind) :: phi
    phi = sqrt(1.0 + (r/c)**2)
  endfunction rbf_multiquadric

!************************************************
  function rbf_invmultiquadric(r,c) result(phi)
! Inverse multiquadric Radial Basis Function    
!************************************************
    implicit none
    real(ikind),intent(in) :: r,c
    real(ikind) :: phi
    phi = (1.0 + (r/c)**2)**(-0.5)
  endfunction rbf_invmultiquadric
    
!***********************************************
  function rbf_invquadratic(r,c) result(phi)
! Inverse Quadratic Radial Basis Function    
!***********************************************
    implicit none
    real(ikind),intent(in) :: r,c
    real(ikind) :: phi
    phi = 1.0/(1.0 + (r/c)**2)
  endfunction rbf_invquadratic
    
!***********************************************
  function rbf_thinplate(r,c) result(phi)
! Thinplate Radial Basis Function    
!***********************************************
    implicit none
    real(ikind),intent(in) :: r,c
    real(ikind) :: phi
    phi = r*r*log(r + 1.0)
  endfunction rbf_thinplate

!************************************************
  function rbf_thinplatespline(r,c) result(phi)
! Thin Plate Spline Radial Basis Function    
!************************************************
    implicit none
    real(ikind),intent(in) :: r,c
    real(ikind) :: phi
    phi = r*r*log(r)
  endfunction rbf_thinplatespline
    
endmodule rbf_lib
