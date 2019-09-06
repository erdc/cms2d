!===================================================================
module rol_def
! Surface roller model variable definition module
!===================================================================
    use prec_def
    implicit none
    save
    
    logical :: roller      !Turns on surface roller model, 0-Off, 1-On
    integer :: irolscheme  !Roller scheme, 1-First-order upwind, 2-Lax-Wendroff, 3-Second-order upwind
    logical :: rolflux     !Include roller flux in wave flux velocity, 0-No, 1-Yes    
    real(ikind) :: br,ceff !Roller and efficiency coefficients    
    real(ikind) :: roltol  !Convergence tolerance
    real(ikind) :: rol_courant !Roller courant number
    real(ikind),allocatable :: delx(:),dely(:),d2x(:),d2y(:) !Internal geometric variables
    real(ikind),allocatable :: Sr(:,:) !2x the roller energy per unit area (Sr=2*Esr) *********
    real(ikind),allocatable :: rxrs(:,:),ryrs(:,:)  !Roller stress gradients
    real(ikind),allocatable :: roldiss(:,:) !Roller dissipation (N/m/s)
    
endmodule rol_def
