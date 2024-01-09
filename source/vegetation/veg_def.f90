!=========================================================
module veg_def
! Vegetation drag for hydrodynamics
!=========================================================
    use prec_def
    implicit none
    save
    
    logical :: veg    
    real(ikind) :: hv,Cdv,Dv,alphac   
    real(ikind), allocatable :: nv(:)   
    character(len=200) :: vegdensfile, vegdenspath  
    
end module veg_def