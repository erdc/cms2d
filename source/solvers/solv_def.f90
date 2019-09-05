!=======================================================================
module solv_def
!=======================================================================
    use prec_def
    implicit none
    save
    
    integer :: iconv  !Solver converged
    
    !Sparse matrix iterative solver (with preconditioners)
    type solver_options
      integer :: iconv          !Convergence state 0-Divergent, 1-Continue, 2-Converged, 3-Reduce time step
      integer :: isolv          !Sparse matrix iterative solver
      integer :: lfil
      real(ikind) :: droptol    !Drop tolerance for preconditioner
      real(ikind) :: permtol    !Permutation tolerance for preconditioner
      real(ikind) :: relaxsor   !Succesive-Over-Relaxation coefficient
      integer :: nswp           !Maximum number of inner iterations
      integer :: nswp0          !Maximum number of inner iterations (for adapative mode)
      real(ikind) :: relax      !Implicit relaxation coefficient      
      real(ikind) :: rmom       !Normalized residual
      real(ikind) :: rmom0      !Normalized residual of first iteration
      real(ikind) :: rmommax    !Maximum normalized residual
      real(ikind) :: rmommin    !Minimum normalized residual
      real(ikind) :: rmomtarget !Target normalized residual
      real(ikind) :: rmomratio  !Maximum ratio for normalized residual
      real(ikind) :: rmomabschg !Maximum absolute change for normalized residual
      real(ikind) :: rmomrelchg !Maximum relative change for normalized residual      
    endtype solver_options
    type(solver_options) :: velsolv,ppsolv,salsolv,Ctksolv,heatsolv
    
    character(len=32) :: asolv(0:8) !Sparse Matrix Solvers
    
    !Sparse coefficient matrix
    real(ikind), allocatable :: aa_matrix(:),alu(:)
    integer, allocatable :: ia(:),ja(:),iua(:)
    
    !Preconditioner
    integer, allocatable :: iaa(:,:)
    integer, allocatable :: nposition1(:,:),nposition2(:,:)
    integer, allocatable :: nup(:),nlow(:),jlu(:),ju(:)
    
    integer:: no_zero
    integer:: lfil
    real(ikind):: droptol
    real(ikind):: permtol
    
endmodule solv_def
