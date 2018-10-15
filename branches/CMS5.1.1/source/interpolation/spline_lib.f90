!==============================================================
module spline_lib
! Cubic Spline Interpolation Library
!
! Description:
!   Calculates a spline interpolant of the form
!    s(x) = a*y(k) + b*y(k+1) + c*y2(k) + d*y2(k+1)
!  where 
!    n = size of the arrays x and y
!    x(1:n) = input data abscissas (in strictly increasing order)
!    y(1:n) = input data ordinates
!    y2(1:n) = output second derivative at data abscissas
!    a = (x(k+1)-xi)/h = 1-b
!    b = (xi-x(k))/h = 1-a
!    c = 1/6*(a^3-a)*h^2
!    d = 1/6*(b^3-b)*h^2
!    h = x(k+1)-x(k)
!  and for x(k)<=xi<=x(k+1)
!
!  This approach is preferred here because it only requires
!  saving the second derivates.
!  Uses the "natural" boundary conditions in which the second
!  order derivatives are set to zero at the boundaries 
!  (i.e. the function becomes linear near the boundaries).
! 
! Contains:
!   spline_fit - Calculates the second order derivatives for 
!                spline interpolations
!   spline_eval - Evaluates the spline interpolant
!
! Usage:
!  call spline_fit(n,x,y,y2)
!  yi = spline_eval(n,x,y,y2,xi,inc)
!
! based on the code Numerical Recipies in Fortran 77:
! The Art of Scientific Computing
!==============================================================
    implicit none
    
contains

!**********************************************************************
    subroutine spline_fit(n,x,y,y2)
!  Calculates the second order derivatives for spline interpolations
!
! Description:
!  The spline interpolant is of the form
!    s(x) = a*y(k) + b*y(k+1) + c*y2(k) + d*y2(k+1)
!  where 
!    n = size of the arrays x and y
!    x(1:n) = input data abscissas (in strictly increasing order)
!    y(1:n) = input data ordinates
!    y2(1:n) = output second derivative at data abscissas
!    a = (x(k+1)-xi)/h = 1-b
!    b = (xi-x(k))/h = 1-a
!    c = 1/6*(a^3-a)*h^2
!    d = 1/6*(b^3-b)*h^2
!    h = x(k+1)-x(k)
!  and for x(k)<=xi<=x(k+1)
! 
! Solves the tridiagonal system of equations
! (x(k)-x(k-1))/6*y2 + (x(k+1)-x(k-1))/3*y2 + (x(k+1)-x(k))/6*y2 =
!  (y(k+1)-y(k))/(x(k+1)-x(k)) - (y(k)-y(k-1))/(x(k)-x(k-1))
!  for k=2:n-1
!
! Uses a natural spline, with zero second derivative 
! on both boundaryies
!
! based on the code Numerical Recipies in Fortran 77:
! The Art of Scientific Computing
!**********************************************************************
    use prec_def
    implicit none
    !Input/Output
    integer,    intent(in) :: n
    real(ikind),intent(in) :: x(n),y(n)
    real(ikind),intent(out) :: y2(n)
    !Internal variables
    integer :: i,k
    real(ikind) :: p,qn,sig,un,u(n)
    
    !LU decomposition and forward substitution
    y2(1)=0.0; u(1)=0.0
    do i=2,n-1 
      sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
      p=sig*y2(i-1)+2.0
      y2(i)=(sig-1.0)/p
      u(i)=(6.0*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1)) &
          /(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
    enddo
    
    !Backwards substitution
    qn=0.0; un=0.0
    y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.0)
    do k=n-1,1,-1 
      y2(k)=y2(k)*y2(k+1)+u(k) 
    enddo
    
    return
    endsubroutine spline_fit
    
!*****************************************************************    
    function spline_eval(n,x,y,y2,xi,inc) result(yi)
! Evaluates the spline interpolant
!
! Description:
!   The spline funciton is given by
!     s(x) = a*y(k) + b*y(k+1) + c*y2(k) + d*y2(k+1)
!   where 
!     n = size of the arrays x and y
!     x(1:n) = input data abscissas (in strictly increasing order)
!     y(1:n) = input data ordinates
!     y2(1:n) = output second derivative at data abscissas
!     a = (x(k+1)-xi)/h = 1-b
!     b = (xi-x(k))/h = 1-a
!     c = 1/6*(a^3-a)*h^2
!     d = 1/6*(b^3-b)*h^2
!     h = x(k+1)-x(k)
!     inc = initial guess for k
!   and for x(k)<=xi<=x(k+1)
!
!  Note: If k is not known it can be left out and 
!  the bisection method will be used to find k.
!
! based on the code Numerical Recipies in Fortran 77:
! The Art of Scientific Computing
!*****************************************************************
    use prec_def
    implicit none
    !Input/Output
    integer,    intent(in) :: n
    real(ikind),intent(in) :: xi,x(n),y(n),y2(n)
    integer,    intent(inout),optional :: inc
    real(ikind) :: yi
    !Internal Variables
    integer :: i,j,k,k1
    real(ikind) :: a,b,h
    
    !Find starting location so that x(k)<=xi<=x(k+1)
    if(present(inc))then !Sequential search
      inc=max(inc,1)
      do k=inc,n-1
        if(xi<=x(k+1)) exit  
      enddo
      inc=k !Save for output
    else !Bisection search
      k=1; j=n+1
      do while(j>k+1)
        i=(k+j)/2
        if(xi<x(k))then
          j=i
        else
          k=i
        endif
      enddo
    endif
    
    !Cubic spline polynomial is now evaluated.
    k1=k+1
    h=x(k1)-x(k)
    a=(x(k1)-xi)/h
    b=1.0-a
    yi=a*y(k)+b*y(k1)+((a**3-a)*y2(k)+(b**3-b)*y2(k1))*(h**2)/6.0
    
    return
    endfunction spline_eval
    
endmodule spline_lib