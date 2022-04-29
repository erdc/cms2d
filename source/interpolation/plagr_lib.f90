!=======================================================================
module plagr_lib
!   Piecewise Lagrangian Basis Polynomial Interpolation Library
!
! Contains
!   plagr_fit - Fits a Lagrange basis piecewise polynomial interpolant
!   plagr_eval - Evaluates a Lagrange basis piecewise polynomial interpolant
!
! Usage:
!   call plagr_fit(m,x,xi,nb,lb,ni,np,k)
!   yi = plagr_eval(m,y,nb,lb,np,k)
!
! written by Alex Sanchez, USACE-CHl
!=======================================================================
    implicit none

contains

!************************************************************************    
    subroutine plagr_fit(m,x,xi,nb,lb,ni,np,k)
! Lagrange basis polynomial interpolation
!
! Description:    
!   Calculates Lagrange basis polynomial weights based on  
!   given by the independent coordinates x and
!   interpolation coordinate xi so that the interpolated
!   value can be calcualted as 
!     yi = sum(lb(1:np+1)*y(k:k+np))
!   where w is the computed weights and y are the values 
!   at the independent coordinate x.
!   An advantage of using Lagrange polynomials is that 
!   the interpolation weights are not a function of the 
!   function values at x. 
!   The first order Lagrange polynomial reduces to
!   linear interpolation.
!
! Input:
!   m    - Length of x and y
!   x(m) - Independent coordinate
!   xi   - interpolation coordinate
!   nb   - Size of Lagrangian basis function array (n>=ni+1)
!   k    - Starting location of interpolation points (use 0 if not known)
!   ni   - Input order of Lagrange polynomial fit
!
! Output:
!   lb(nb) - Lagrange basis functions (weights)
!   np    - Output order of polynomial Lagrange polynomial fit limited by 
!             no = min(ni,nw-1,m-1)
!   k     - Starting location of interpolation points
!
! Change Log:
!   Bug fix (12/16/13 by Alex) - Problem for neighest neighbor extrapolation for 
!      first order interpolation. Replaced lb(nb)=1.0 with lb(np+1)=1.0.
!
! written by Alex Sanchez, USACE-CHL
!************************************************************************
    use prec_def
    implicit none
    !Input/Output    
    integer,    intent(in)    :: m     !Length of x coordinates
    real(ikind),intent(in)    :: x(m)  !Input coordinates
    real(ikind),intent(in)    :: xi    !Interpolation coordinate
    integer,    intent(in)    :: nb    !Length of Lagrangian basis functions
    real(ikind),intent(out)   :: lb(nb) !Lagrange polynomial weights   
    integer,    intent(in)    :: ni    !Input order of polynomial
    integer,    intent(out)   :: np    !Output order of polynomial
    integer,    intent(inout) :: k     !Starting point of interpolation points
    !Internal Variables
    integer :: i,j,np1,ik1,jk1,k1
        
    np=min(ni,m-1,nb-1) !Limit order of polynomial
    
    np1=np+1
    !Avoid extrapolation
    if(xi<x(1))then
      lb(1)=1.0
      lb(2:nb)=0.0
      k=1
      return
    elseif(xi>x(m))then
      lb(1:nb)=0.0
      lb(np+1)=1.0 !Alex, Bug fix (11/13/13): Replaced index nb with np+1
      k=m-np
      return   
    endif
    
    !Find starting location so that x(k)<=xi<=x(k+np)
    if(k>=1)then !Sequential search
      i=k
      do k=i,m-np !Forwards search
        if(xi<=x(k+np)) exit  
      enddo      
      if(xi<x(k) .or. xi>x(k+np))then !Check
        do k=i,1,-1 !Backwards  search
          if(x(k)<=xi) exit  
        enddo 
        if(xi<x(k) .or. xi>x(k+np))then !Check
          write(*,*) 'WARNING: Problem calculating starting location in plagr_fit'
          write(*,*) '  Using Bisection search'
          k=0
        endif
      endif
    endif
    if(k<=0)then !Bisection search
      i=1; j=m-np+1
      do while(j>i+1)
        k=(i+j)/2
        if(xi<x(k))then
          j=k
        else
          i=k
        endif
      enddo
      k=i
      if(xi<x(k) .or. xi>x(k+np))then !Check
        write(*,*) 'WARNING: Problem calculating starting location in plagr_fit'
      endif
    endif
    
    k1=k-1
    
    !Calculate coefficients    
    lb(np1+1:nb)=0.0
    do i=1,np1
      lb(i)=1.0
      ik1=i+k1
      do j=1,np1
        if(i/=j)then
          jk1=j+k1
          lb(i)=lb(i)*(xi-x(jk1))/(x(ik1)-x(jk1))
        endif
      enddo
    enddo
    
    lb = lb/sum(lb(1:np+1))
    
    return
    end subroutine plagr_fit
    
!************************************************************************    
    function plagr_eval(m,y,nb,lb,np,k) result(yi)
! Lagrange basis polynomial interpolation
!
! Description:    
!   Calculates Lagrange basis polynomial weights based on  
!   given by the independent coordinates x and
!   interpolation coordinate xi so that the interpolated
!   value can be calcualted as 
!     yi = sum(lb(1:np+1)*y(k:k+np))
!   where w is the computed weights and y are the values 
!   at the independent coordinate x.
!   An advantage of using Lagrange polynomials is that 
!   the interpolation weights are not a function of the 
!   function values at x. 
!   The first order Lagrange polynomial reduces to
!   linear interpolation.
!
! Input:
!   m    - Length of x and y
!   x(m) - Independent coordinate
!   xi   - interpolation coordinate
!   nb   - Size of Lagrangian basis function array (n>=ni+1)
!   k    - Starting location of interpolation points (use <=1 if not known)
!   ni   - Input order of Lagrange polynomial fit
!
! Output:
!   lb(nw) - Lagrange basis functions (weights)
!   no    - Output order of polynomial Lagrange polynomial fit limited by 
!             no = min(ni,nw-1,m-1)
!   k     - Starting location of interpolation points
!
! written by Alex Sanchez, USACE-CHL
!************************************************************************
    use prec_def
    implicit none
    !Input    
    integer,    intent(in) :: m     !Length of x coordinates
    real(ikind),intent(in) :: y(m)  !Input coordinates  
    integer,    intent(in) :: nb    !Length of Lagrangian basis functions
    real(ikind),intent(in) :: lb(nb) !Lagrange polynomial weights           
    integer,    intent(in) :: np    !Output order of polynomial
    integer,    intent(in) :: k     !Starting point of interpolation points
    !Output
    real(ikind) :: yi
    
    yi = sum(lb(1:np+1)*y(k:k+np))
    
    return
    end function plagr_eval
    
end module plagr_lib    