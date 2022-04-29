!=====================================================
module interp_def
!Spatial Interpolation Variable definitions
!=====================================================
    use prec_def    
    implicit none
    save
    
    !Cell-to-face interpolation
    real(ikind), allocatable :: fintp(:,:) !fc2f
    
    !Cell-to-node interolation
    integer :: nnintp !Method 0-None Specified, 1-Inverse Area, 2-Inverse Distance, 3-Least-Squares
    real(ikind), allocatable :: wc2n(:,:)

end module interp_def