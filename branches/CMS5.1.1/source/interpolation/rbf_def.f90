!============================================
module rbf_def
! Augmented Radial Basis Function Interpolation
! variable defintion module
!
! written by Alex Sanchez, USACE-CHL
!============================================
    use prec_def
    implicit none
    
    type rbf_type
      integer :: ndim   !# of dimensions
      integer :: npts   !# of data points
      integer :: nsize  !Order of system (# of equations)
      real(ikind), allocatable :: x(:,:)  !Data coordinates (ndim,npts)
      real(ikind), allocatable :: xi(:,:) !Interpolation coordinates (ndim,nsize)
      real(ikind), allocatable :: A(:,:)  !LU decomposition (nsize,nsize)
      integer,     allocatable :: ind(:)  !Index for LU decomposition
      real(ikind), allocatable :: coef(:) !Interpolation coefficients (nsize)
      real(ikind) :: const   !radial basis function constant
      real(ikind) :: smooth  !smoothing parameter, 0=none
      character(len=30) :: fun    !radial basis function
    endtype rbf_type
        
endmodule rbf_def