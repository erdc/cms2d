!==================================================================
module math_lib
! Math library
!==================================================================    
    implicit none
    
contains    

!******************************************    
    function lognpdf(x,mu,sigma) result(y)
! Lognormal probability density function
!******************************************   
    use const_def, only: twopi
    use prec_def
    implicit none
    real(ikind) :: x,mu,sigma,y
    
    y = 1.0/(x*sigma*sqrt(twopi))
    y = y*exp(-0.5*((log(x)-mu)/sigma)**2)
    
    return
    endfunction lognpdf
    
!************************************************************************       
    function erfinv(x) result(y)
!Inverse error function
!note that erfc(x) = 1 - erf(x), and erfc(erf^-1(1-x))=x
!The largest error is about 0.00012
! Approximate function from Winitzki, Sergei (6 February 2008) 
! "A handy approximation for the error function and its inverse"
!************************************************************************   
    use const_def, only: pi
    use prec_def
    implicit none
    real(ikind), parameter :: a=0.140012, b=4.546894 !a=8*(pi-3)/(3*pi*(4-pi)),b=2/(pi*a)
    real(ikind) :: c,x,y
    
    c = log(1.0-x**2)
    y = x/abs(x)*sqrt(sqrt((b+c/2.0)**2-c/a)-(b+c/2.0))
    
    return    
    endfunction erfinv
        
!*******************************************       
    function logninv(p,mu,sigma) result(x)
!Inverse lognormal distribution
!*******************************************   
    use const_def, only: sqrttwo
    use prec_def
    implicit none
    real(ikind) :: logx0,p,mu,sigma,x 
    
    logx0 = -sqrttwo*erfinv(1.0-2.0*p)
    x = exp(sigma*logx0 + mu)
    
    return    
    endfunction logninv
    
!********************************************************************     
    function erf(x) result(y)
!Calculates the error function
!used in wave breaking formulation of Alsina and Baldock (2007)
!********************************************************************  
    use const_def, only: pi
    use prec_def    
    implicit none
    integer :: k
    real :: dk,E,S,fac,x,eps,y
    real, parameter :: tol=0.0001
    
    fac = 2.0*x/sqrt(pi)
    eps = tol/fac
    if(x>3.0)then
        y = 1.0
        return    
    endif
    E = 1.0
    S = 1.0
    do k=1,40
        dk = real(k)
        E = -((2.0*dk-1.0)*x**2)/((2.0*dk+1.0)*k)*E
        if(abs(E)<=eps)then
            exit
        endif
        S = S + E
    enddo
    y = fac*S
    if(y>1.0) y = 1.0
    
    return
    endfunction erf
    
!***************************************************    
    function ei(x) result(y)
! Computes the exponential integral Ei(x)
!     Input :  x  --- Argument of Ei(x)
!     Output:  EI --- Ei(x) ( x > 0 )
!***************************************************
    use prec_def
    implicit none
    real(ikind),intent(in) :: x
    integer     :: k
    real(ikind) :: r,y    

    if(x==0.0)then
      y = -1.0e+25
    elseif(x<=40.0)then
      y = 1.0
      r = 1.0
      do k=1,100
        r = r*k*x/(k+1)** 2
        y = y + r
        if(abs(r/y) <= 1.0e-15) exit
      enddo
      y = 0.5772156649015328 + log(x) + x*y
    else
      y = 1.0
      r = 1.0
      do k = 1,20
        r = r*k/x
        y = y + r
      enddo
      y = exp(x)/x*y
    endif
    
    return
    endfunction ei
    
!****************************************************    
    function beta(z,w)
!Returns the value of the beta function B(z,w).
! uses gammln
!****************************************************
    use prec_def
    real(ikind):: beta,w,z
    !real(ikind):: gammln
    
    beta=exp(gammln(z)+gammln(w)-gammln(z+w))
    
    return
    endfunction beta

!****************************************************       
    function gammln(xx)
!Returns the value log[Gamma(xx)] for xx > 0
!****************************************************       
    use prec_def
    integer:: j
    real(ikind):: gammln,xx    
    real(ikind):: ser,tmp,x,y
    real(ikind),parameter :: stp       =  2.5066282746310005_ikind
    real(ikind),parameter :: coef(1:6) = (/76.18009172947146_ikind,&
       -86.50532032941677_ikind,24.01409824083091_ikind,-1.231739572450155_ikind,&
       0.1208650973866179e-2_ikind,-0.5395239384953e-5_ikind/)
    
    x=xx
    y=x
    tmp=x+5.5_ikind
    tmp=(x+0.5)*log(tmp)-tmp
    ser=1.000000000190015_ikind
    do j=1,6
      y=y+1.0_ikind
      ser=ser+coef(j)/y
    enddo
    gammln=tmp+log(stp*ser/x)
    
    return
    endfunction gammln
    
!***************************************************************    
    function e1x(x)
! Purpose: Compute exponential integral E1(x)
! Input :  x   --- Argument of E1(x)
! Output:  CE1 --- E1(x)
!***************************************************************
    use prec_def
    real(ikind) :: x,e1x
    complex(ikind) :: z,ce1
    
    z = x
    call e1z(z,ce1)
    e1x = real(ce1,kind=ikind)
    
    return
    endfunction e1x
    
!***************************************************************    
    subroutine e1z(z, ce1)
! Purpose: Compute complex exponential integral E1(z)
! Input :  z   --- Argument of E1(z)
! Output:  CE1 --- E1(z)
!***************************************************************
    use prec_def
    implicit none
    complex(ikind), intent(in)   :: z
    complex(ikind), intent(out)  :: ce1
    complex(ikind)  :: cr, ct, ct0
    real(ikind)     :: a0, x
    integer       :: k
    real(ikind), parameter :: pi = 3.141592653589793, el = 0.5772156649015328

    x = real(z)
    a0 = abs(z)
    if(a0 == 0.0)then
      ce1 = (1.0D+30,0.0)
    elseif(a0 <= 10.0 .OR. x < 0.0 .AND. a0 < 20.0)then
      ce1 = (1.0,0.0)
      cr = (1.0,0.0)
      do  k = 1, 150
        cr = -cr * k * z / (k+1) ** 2
        ce1 = ce1 + cr
        if(abs(cr) <= abs(ce1)*1.0D-6) exit
      enddo
      ce1 = -el - log(z) + z * ce1
    else
      ct0 = (0.0,0.0)
      do  k = 120, 1, -1
        ct0 = k / (1.0 + k/(z+ct0))
      enddo
      ct = 1.0 / (z+ct0)
      ce1 = exp(-z) * ct
      if(x <= 0.0 .AND. aimag(z) == 0.0) ce1 = ce1 - pi * (0.0,1.0)
    endif
    
    return
    endsubroutine e1z    
    
!***************************************************************  
    subroutine sortup(n,x,ind)
! Sort vector in ascending order
!***************************************************************  
    use prec_def
    implicit none
    !Input
    integer,intent(in) :: n
    !Output
    real(ikind),intent(inout) :: x(1:n)
    integer,optional :: ind(n)
    !Internal
    integer :: i,iswap(1),iswap2,iord(n),itemp
    real(ikind) :: xtemp
    intrinsic minloc
    
    do i=1,n
      iord(i) = i
    enddo
    do i=1,n-1
      iswap=minloc(x(i:n)) !returns the index value for the minimum value
      iswap2=iswap(1)+i-1 !global index
      if(iswap2/=i)then        
        itemp=iord(i)
        iord(i)=iord(iswap2)
        iord(iswap2)=itemp
        xtemp=x(i)
        x(i)=x(iswap2)
        x(iswap2)=xtemp
      endif
    enddo
    
    if(present(ind)) ind = iord
    
    return
    endsubroutine sortup
    
!***************************************************************  
    function avgval(n,x)
! Average value
!***************************************************************  
    use prec_def
    implicit none
    integer,intent(in) :: n
    real(ikind),intent(in) :: x(*)
    real(ikind) :: avgval
    
    avgval = sum(x(1:n))/real(n,kind=ikind)
        
    return
    endfunction avgval
    
endmodule math_lib    