!===========================================================
module secant_lib
! Secant Method library for solving the following problem
!    0 = var - fvar(c,x)
!  for x, where c is given and var=fvar(c,x)
!
! Contains:
!   secant - Secant Method
!
! written by Alex Sanchez, USACE-CHL
!===========================================================
   implicit none
   
contains
!****************************************************
    subroutine secant(fvar,var,c,x,f,eps,tol,iter)
! Secant Method
! 
! Description:
!   Solves for x in the the following problem
!     0 = var - fvar(c,x)
! where
!   c = constant
!   x(3) = indepedent variable which is being solved for
!          On entry, the first two elements contain the 
!          first two guesses and on output, the third
!          element contains the final value
!   fvar(c,x) = user-defined function
!   var = fvar(c,x(3))
!   eps > abs(x(2)-x(3)) => exit criteria
!   tol > abs(f(2)) => exit criteria
!   iter = On entry equal to the maximum number of iterations
!          On exit equal to the number of iterations
! 
! written by Alejandro Sanchez, USACE-CHL
!****************************************************
    use prec_def
    implicit none
    real(ikind),intent(inout) :: x(3),f(2),eps,tol
    real(ikind),intent(in) :: var,c
    integer,intent(out) :: iter
    !Internal
    integer :: i
	real(ikind) :: err
    
    interface
      function fvar(c,x)
        use prec_def
        implicit none
        real(ikind),intent(in) :: c,x
        real(ikind) :: fvar
      endfunction
    endinterface

	f(1) = var - fvar(c,x(1))
	f(2) = var - fvar(c,x(2))
	do i=1,iter		
      x(3) = x(2) - f(2)*(x(1) - x(2))/(f(1) - f(2))  !Secant formula
      err = abs(x(2)-x(3)) !Change
      if(err<eps) exit
      x(1) = x(2)
      x(2) = x(3)
      f(1) = f(2)
      f(2) = var - fvar(c,x(2))
      if(abs(f(2))<tol) exit !Error
    enddo
	
    iter = i
    eps = err
    tol = f(2)
    
	return
    endsubroutine secant

endmodule secant_lib    