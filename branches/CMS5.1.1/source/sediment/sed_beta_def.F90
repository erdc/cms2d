!=========================================================
module beta_def    
! Total-load correction factor look-up tables
!=========================================================    
    use prec_def
    implicit none
    save
    
    integer :: nzap,nphip,nr
    real(ikind) :: dzap,logzap
    real(ikind),allocatable :: zap_table(:)        
    
    !Exponential distribution  
    real(ikind) :: dphip,phip0
    real(ikind),allocatable :: phip_table(:)    
    
    !Rouse distribution
    real(ikind):: dr,r0   
    real(ikind),allocatable :: r_table(:)   
    
    !Suspended load correction factor
    real(ikind),allocatable :: betas_table(:,:)

endmodule beta_def