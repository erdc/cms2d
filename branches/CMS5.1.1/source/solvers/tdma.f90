!*************************************************
	subroutine tdma(n,e,f,g,r,x)
! TDMA Algorithm    
!
! Description: Solves a tridiagnol system of equations
! using the TDMA (Thomas) Algorithm
! Tridiagnol Matrix Elements:
!  [e,  f,  g,  |  r]
!
! written by Alex Sanchez, USACE-CHL
!*************************************************    
    implicit none
    !Input/Output
    integer,intent(in) :: n
    real(ikind),intent(inout) :: e(n),f(n),g(n),r(n)
    real(ikind),intent(out) :: x(n)
    !Internal
    integer :: k
    
    !LU Decomposition and Foward substiution
    do k=2,n
	  e(k) = e(k)/f(k-1)
	  f(k) = f(k) - e(k)*g(k-1)
      r(k) = r(k) - e(k)*r(k-1) !Foward substiution
    enddo

    !Back Substitution
    x(n) = r(n)/f(n)
    do k=n-1,1,-1
	  x(k) = (r(k) - g(k)*x(k+1))/f(k)
    enddo

    return
    endsubroutine tdma