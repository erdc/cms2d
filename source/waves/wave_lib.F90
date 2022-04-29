!=================================================================
module wave_lib
! CMS wave library
! written by Alex Sanchez, USACE-CHL
!=================================================================    
    implicit none
    
contains

!******************************************************************       
    function waveorb_linear(d,H,T,L) result(uw)
! Bottom wave orbital velocity based on linear wave theory
!
! Input:
!   d = Depth [L]
!   H = Wave height [L]
!   T = Wave period [T]
!   L = Wave length [L]
!
! Output:
!   uw = Bottom wave orbital velocity [L/T]
!
! Author: Alex Sanchez, USACE-CHL
!******************************************************************       
    use const_def, only: pi,twopi
    use prec_def
    implicit none
    real(ikind),intent(in) :: d,H,T,L
    real(ikind) :: kh,uw
    
    kh = min(twopi*d/L,10.0)
    uw = pi*min(H,0.78*d)/(max(T,2.0)*sinh(kh))             
!!   uw = pi*H/(max(T,2.0)*sinh(kh))  

    return
    end function waveorb_linear   

!******************************************************************************    
    subroutine waveorbrep_jonswap(Hs,Tp,wd,gam,h,u,v,ubr,Tbr)
! Calculates ubr and Tbr from Hs and Tp based on the JONSWAP spectra.
! The method is more accurate but relatively expensive.
!
! Input:
!   Hs = Significant wave height [L]
!   Tp = Peak period [T]
!   wd = Wave direction (Cartesian, going to) [rad]
!   gam = Peak enhancement factor (~3.3) [-]
!   h = Water depth [L]
!   u,v = Curreng velocities [L/T]
!
! Output:
!   ubr = representative bottom orbital velocity [L/T]
!   Tbr = representative bottom wave period [T]
!
! Reference:
!   Wiberg, P., Sherwood, C. (2008). Calculating wave-generated bottom orbital
!   velocities from surface-wave parameters, Computers & Geosciences, 34, 1243-1262.
!
! Author: Alex Sanchez, USACE-CHL
!******************************************************************************
    use const_def, only: twopi
    use prec_def
    implicit none
    !Input
    real(ikind),intent(in):: Hs,Tp,wd,gam,h,u,v
    !Output
    real(ikind),intent(out):: ubr
    real(ikind),intent(out),optional:: Tbr
    !Internal variables
    integer :: i
    integer, parameter :: nb = 97
    real(ikind),parameter :: dffp = 0.05_ikind
    real(ikind),parameter :: xj = 3.2831_ikind
    real(ikind):: su,sumsu,sumfsu,twosigsq,m0fp4XJ,kh
    real(ikind):: ffp,f,fp,df,t1,t2,sf
        
    fp = 1.0_ikind/Tp
    df = dffp*fp
    m0fp4XJ = (Hs**2/16.0_ikind)*fp**4*XJ
    sumsu = 0.0_ikind
    sumfsu = 0.0_ikind
    do i=1,nb
      ffp = 0.2_ikind + (i-1)*dffp
      f = ffp*fp
      kh = h*wavenumber(f,wd,h,u,v)
      if(ffp>1.0)then        
        twosigsq = 0.0162_ikind !sig = 0.09
      else
        twosigsq = 0.0098_ikind !sig = 0.07
      endif
      t1 = exp(-1.25*ffp**(-4))
      t2 = gam**exp(-((ffp-1.0_ikind)**2)/twosigsq)
      sf = m0fp4XJ/f**5*t1*t2
      su = (twopi*f/sinh(kh))**2*sf  !Eq. 5
      sumsu = sumsu + su
      sumfsu = sumfsu + f*su
    enddo      
    ubr = sqrt(2.0_ikind*sumsu*df) !Eq. 6
    if(present(Tbr))then
      Tbr = sumsu/sumfsu !Eq. 9
    endif
    
    return
    end subroutine waveorbrep_jonswap
    
!****************************************************************************** 
    function waveorbrep_ss87(g,h,Hs,Tp) result(Ubr)
! Calculates the representative bottom orbital velocity amplitude Ubr
! equal to Ubr = sqrt(2)*Urms where Urms is the root-mean-squared bottom
! orbital velocity amplitude Urms = std(ub)
!
! Input:
!   g = Gravity [L/T^2]
!   h = Water depth [L]
!   Hs = Significant wave height [L]
!   Tp = Peak period [T]
!
! Output:
!   ubr = representative bottom orbital velocity [L/T]
!
! Reference:
!   Soulsby, R.L. (1987). Calculating bottom orbital velocities beneath waves,
!      Coastal Engineering, 11, 371-380.
!
! Notes:
!   The equation was obtained by optimizing the three coefficients of the equation
!   bellow using Matlab's fminsearch function and digitized points from the curve
!   presented in Soulsby's paper. 
! 
!   Soulbsy's definition (used here) of Urms is Urms^2 = int(Su*df) 
!   while others use Urms^2 = int(2*Su*df). 
!   This leads to a difference of a factor of sqrt(2).
!
!   Does not consider wave-current interaction.
!
! Authour: Alex Sanchez, USACE-CHL
!******************************************************************************
    use prec_def
    implicit none
    real(ikind),intent(in) :: g,h,Hs,Tp
    real(ikind) :: Tn,Ubr
    
    Tn = sqrt(h/g)
    Ubr = 0.19_ikind*(1.0_ikind+tanh(-7.76_ikind*Tn/Tp+1.34_ikind))*Hs/Tn
    
    return
    end function waveorbrep_ss87

!******************************************************************
    function wavelength(wa,wd,h,u,v,tolinp) result(wlen)
! Calculates the wave length by solving the wave dispersion relation
!    sig^2 = grav*wk*tanh(wk*h)
!
! Usage:
!  wlen = wavelength(wa,wd,h,u,v,tolinp)
!  wlen = wavelength(wa,wd,h,u,v)
!
! Description:
!  Solves the wave dispersion relation 
!   sig^2 = grav*wk*tanh(wk*h)
!  where 
!   sig = wa - k*cos(wd)*u - k*cos(wd)*v = wa - k*uk
!   sig = relative angular frequency [rad/s]
!   wk = wave number [rad/m]
!   wa = absolute angular frequency [rad/s]
!   wd = wave direction (Cartesian, going to) [rad]
!   u = current velocity in x direction [m/s]
!   v = current velocity in y direction [m/s]
!   uk = cos(q)*u + sin(q)*v [m/s]
!
! Uses the Newton-Raphson Method is given by:
!   wk(n+1) = wk(n) - f(wk(n))/fp(wk(n))
! where
!   f = grav*wk*tanh(wk*h) - sig^2
!   fp = grav*(tanh(h*k)-h*k*(tanh(h*k)^2-1)) + 2*uk*sig
!
! Makes initial guess using 
!   Guo, J. (2002) Simple and explicit solution of wave
!      dispersion, Coastal Engineering, 46(2), 71-74.
!
! written by Alex Sanchez, USACE-CHL
!************************************************************
    use const_def, only: twopi
    use prec_def
    implicit none
    !Input
    real(ikind),intent(in) :: wa,wd,h,u,v
    real(ikind),intent(in),optional :: tolinp
    !Output
    real(ikind) :: wlen
    !Internal
    real(ikind) :: tol,wk

    if(present(tolinp))then
      tol = tolinp !user-specified
    else    
      tol = 1.0e-3 !default
    endif
    wk = wavenumber(wa,wd,h,u,v,tol)
    wlen = twopi/wk

    return
    end function wavelength

!**************************************************
    subroutine wave_Hmax(d,s,wk,Hmax)
! Calculates the maximum wave height used for wave
! breaking calculations.
!
! Input:
!   d - total water depth
!   s - bed slope in the direction of waves
!
! Output:
!   Hmax - Maximum wave height
!
! written by Alex Sanchez, USACE-CHL
!**************************************************    
    implicit none
    real :: y !,y2
    real,intent(in) :: d,s,wk
    real,intent(out) :: Hmax   
    
    y=0.76*wk*d+0.29  !Grasmeijer 
    !y=0.5*wk*d+0.64   !Sanchez
    !y=0.7*wk*d+0.4    !new, visser
    !y=0.7*wk*d+0.5    !new, LSTF
    !y=0.76*wk*d+0.35  !new, LIP11D 1B
    y=max(0.64,y)     !Important. Needed for Duck 1990 test case. (Alex 10/09/2014)
    
    !!Ostendorf and Madsen (1979)
    !if(s>0.0)then
    !  y2=0.8+5.0*s
    !  y=min(y,min(y2,1.3)) 
    !endif
    
    !Constant
    !y=0.6   !for irregular waves
    !y=0.64 !for regular waves
    !y=0.75 !for LIP11D 1B
    !Hmax=0.9*d
    !y = 0.9
    Hmax=0.88/max(wk,1.e-6)*tanh(y*wk*d/0.88)
    !Hmax=max(Hmax,0.64*d)
    
    return
    end subroutine wave_Hmax

!*********************************************************************
      subroutine wavebreak_bj78(Hmax,Hrms,fp,Qb,Db,ib,cab)
! Calculates the wave breaking dissipation due to breaking waves
! according to the modified BJ, Baldock (1998)
! written by Alex Sanchez, USACE-CHL
!*********************************************************************
!    use flow_def, only: grav
      implicit none
      integer,intent(out) :: ib
      real,intent(in) :: Hmax,Hrms,fp
      real,intent(out):: Qb,Db,cab        
      real, parameter :: alfabj=1.0
    
      if(Hrms<1.e-4)then !Needed to avoid divid by zero
        ib=0; Db=0.0; cab=0.0
        return
      endif
    
      !Fraction of breaking waves
      if(Hrms>=Hmax)then
        Qb=1.0
      else
      !!val=(Hrms/Hmax)**2
      !!Qb=exp((0.9698*val-1.0)/(0.6574*val)) !Norm of residuals = 0.1384
      !!Qb=exp((0.9698*val-1.0)/(0.6574*val))/0.9551  !Norm of residuals = 0.2343
      Qb=15.7053*exp(-2.754/(Hrms/Hmax)**1.342) !Norm of residuals = 0.0618
      endif
    
      !Dissipation rate, note rho not included and negative, 2.4525=0.25*9.81
      Db=-2.4525*alfabj*fp*Qb*Hmax*Hmax !Note sign
    
      !Dissipation coefficient, 0.8154943935=8.0/9.81
      cab=-0.8154943935*Db/Hrms/Hrms
    
      !Breaking Index
      if(Qb>0.02)then !Good LSTF
      !if(Qb>0.0001)then !Good LIP11D
        ib=1      
      else
        ib=0; Qb = 0.0; Db=0.0; cab=0.0
      endif  
    
      return
    end subroutine wavebreak_bj78
    
!**********************************************************************    
      subroutine wavebreak_jb07(d,Hmax,Hrms,fp,Qb,Db,ib,cab) 
!Calculates the wave breaking dissipation due to breaking waves
!according to Alsina and Baldock (2007) and Janssen and Battjes (2007)
!**********************************************************************
    use math_lib, only: erf
      implicit none
      integer,intent(out) :: ib    
      real,intent(in) :: d,Hmax,Hrms,fp
      real,intent(out) :: Qb,Db,cab
      real, parameter :: B = 1.0
      real :: R
    
      if(Hrms<1.e-4)then !Needed to avoid divid by zero
        ib=0; Db=0.0; cab=0.0
        return
      endif
    
      !Fraction of breaking waves
      R=Hmax/Hrms
      Qb=1.0+0.752252778*(R**3+1.5*R)*exp(-R**2)-erf(R) !0.752252778=4/3/sqrt(pi)  !Bug fix (Alex 02/24/2015)
      Qb=max(min(Qb,1.0),0.0)
    
    !Wave breaking dissipation, 2.601269071=3.0/16.0*sqrt(pi)*9.81, note density not included
      Db=-2.601269071*B*fp/max(d,1.0e-4)*Qb*Hrms**3 !Note sign
    
      !Dissipation coefficient, 0.8154943935=8.0/9.81
    cab=-0.8154943935*Db/Hrms/Hrms
    
    !Breaking Index
      !if(Qb>0.05)then
    if(Qb>0.02)then !Good LSTF
!    if(cab>0.1)then
      ib=1
      else
        ib=0; Qb = 0.0; Db=0.0; cab=0.0
      endif
    
      return
      end subroutine wavebreak_jb07
    
!***********************************************************    
    function wavenumber(wa,wd,h,u,v,tolinp) result(wk)
! Calculates the wave number by solvingthe wave dispersion relation
!    sig^2 = grav*wk*tanh(wk*h)
!
! Usage:
!  wk = wavenumber(wa,wd,h,u,v,tolinp)
!  wk = wavenumber(wa,wd,h,u,v)
!
! Description:
!  Solves the wave dispersion relation 
!   sig^2 = grav*wk*tanh(wk*h)
!  where 
!   sig = wa - k*cos(wd)*u - k*cos(wd)*v = wa - k*uk
!   sig = relative angular frequency [rad/s]
!   wk = wave number [rad/m]
!   wa = absolute angular frequency [rad/s]
!   wd = wave direction (Cartesian, going to) [rad]
!   u = current velocity in x direction [m/s]
!   v = current velocity in y direction [m/s]
!   uk = cos(q)*u + sin(q)*v [m/s]
!
! Uses the Newton-Raphson Method is given by:
!   wk(n+1) = wk(n) - f(wk(n))/fp(wk(n))
! where
!   f = grav*wk*tanh(wk*h) - sig^2
!   fp = grav*(tanh(h*k)-h*k*(tanh(h*k)^2-1)) + 2*uk*sig
!
! Makes initial guess using 
!   Guo, J. (2002) Simple and explicit solution of wave
!      dispersion, Coastal Engineering, 46(2), 71-74.
!
! written by Alex Sanchez, USACE-CHL
!************************************************************
    use const_def, only: twopi
    use diag_def
    use flow_def, only: grav
    use prec_def
    implicit none
    !Input
    real(ikind),intent(in) :: wa,wd,h,u,v
    real(ikind),intent(in),optional :: tolinp
    !Output
    real(ikind) :: wk
    !Internal
    integer :: i,k,kunit(2)
    real(ikind) :: xi,yi,wki,uk,sig,f,fp,tol,err,tanhkh
    logical :: isnankind
    
    !Initial wave number guess based on Guo (2002)
    xi = wa/sqrt(grav/h) !=h*wa/sqrt(g*h) 
    yi = xi*xi/(1.0-exp(-xi**2.4908))**0.4015
    wki = yi/h !Initial guess without current-wave interaction 
    
    !Current velocity in wave direction    
    uk = cos(wd)*u + sin(wd)*v 
    
    !Simple correction for currents (tends to overcorrect, hense the relaxation)
    wki = 0.3*wki + 0.7*(wa-wki*uk)**2/grav/tanh(wki*h)       
    
    !Initalize loop variables
    wk = wki   !save for iteration
    if(present(tolinp))then
      tol = tolinp !user-specified
    else    
      tol = 1.0e-3 !default
    endif
    
    !Newton-Raphson iterations
    do i=1,6
      sig = wa-wk*uk
      tanhkh = tanh(wk*h)  
      f = grav*wk*tanhkh - sig**2
      fp = grav*(tanhkh - h*wk*(tanhkh**2-1.0)) + 2*uk*sig
      wk = wki - f/fp
      if(wk<1.0e-6 .or. wk>1.0e6)then !Did not converge
          !wk = 0.0 !Set to zero
          !wk = yi/h !Use initial guess
          wk = yi/h*2.0 !Use a small wave length (large wave number) to avoid divide by zeros
          return
      elseif(isnankind(wk))then
         wk = yi/h !Use a small wave length (large wave number) to avoid divide by zeros
         kunit = (/6,dgunit/)
         open(dgunit,file=dgfile,access='append') 
         do k=1,2
           write(kunit(k),*) 'ERROR: Calculated NaN value for wk in wavenumber'
           write(kunit(k),*) 'wa=',wa
           write(kunit(k),*) 'wd=',wd
           write(kunit(k),*) 'h=',h
           write(kunit(k),*) 'u=',u
           write(kunit(k),*) 'v=',v
         enddo
         close(dgunit)
         return 
      endif          
      err = abs(wk-wki)
      if(err < tol) return
      wki = wk      
    enddo

    return
    end function wavenumber

!********************************************************************
    subroutine cnoidal_sol(Hgt,h,T,Ur,L,c,B)
! Solves the Cnoidal wave problem and calculates 
! several wave properties.
!
! Input:
!   Hgt - wave heoght [m]
!   h - water depth [m]
!   T - wave period [s]
!
! Output:
!  Ur - Ursell number [-]
!  L - Wave length [m]
!  c - Wave celerity [m/s]
!  B - Wave shape parameter [-]
!
! References:
!   Abramowitz, M., and Stegun, I.A. 1964. Handbook of 
!     mathmematical functions. Dover Publications, New York.
!   Svendsen, I.A. 2006. Introduction to nearshore hydrodynamics. 
!     World Scientific. 722 pp.
!
! written by Alex Sanchez, USACE-CHL
!********************************************************************
    use flow_def, only: grav
    use prec_def
    implicit none
    integer :: i
    !Input/Output
    real(ikind),intent(in) :: Hgt,h,T
    real(ikind),intent(out) :: Ur,L,c,B
    !Internal variables
    real(ikind),parameter:: tol=0.001
    real(ikind),parameter:: a0=1.3862944, a1=0.1119723, a2=0.0725296
    real(ikind),parameter:: b0=0.5, b1=0.1213478, b2=0.0288729
    real(ikind),parameter:: e1=0.4630151, e2=0.1077812
    real(ikind),parameter:: f1=0.2452727, f2=0.0412496
    real(ikind) :: A,Urold,em1,em,f,fp,dKdm1,E,fK
    real(ikind) :: tlogem1,Hgtoh3,htsgh,bm
    
    A = 1.0 !Initial value
    htsgh = h*T*sqrt(grav/h)
    L = htsgh*sqrt(1.0+Hgt/h*A) !wave length
    Hgtoh3 = Hgt/h**3
    Ur = Hgtoh3*L**2 !Ursell number
    em1 = exp((a0-sqrt(0.1875*Ur))/b0)
    em1 = min(em1,1.0)
    tlogem1 = log(em1)
    !write(*,*) 'em1 = ',em1,', Ur = ',Ur
    bm = b0+em1*(b1+b2*em1)
    fK = (a0+a1*em1+a2*em1*em1)-bm*tlogem1
    do i=1,10
      f = 1.0 - em1 - 0.1875*Ur/fK**2 !3/16 = 0.1875
      dKdm1 = a1 + 2.0*a2*em1+(b1+2*b2*em1)*tlogem1 - bm/em1
      fp = 0.375*Ur/fK**3*dKdm1 - 1.0  !3/8 = 0.375
      em1 = em1 - f/fp !Newton-Raphson method
      em1 = min(em1,1.0)
      em = 1.0 - em1
      tlogem1 = log(em1)
      bm = b0 + em1*(b1+b2*em1)
      E = (1.0+em1*(e1+e2*em1)) - em1*(f1+f2*em1)*tlogem1
      fK = (a0+em1*(a1+a2*em1)) - bm*tlogem1
      A = 2.0/em - 1.0 - 3.0/em*E/fK
      L = htsgh*sqrt(1.0+Hgt/h*A) !wave length
      Urold = Ur !Save for convergence test
      Ur = Hgtoh3*L**2 !Ursell number
      write(*,*) 'i = ',i,', em1 = ',em1,', Ur = ',Ur
      if(abs(Urold-Ur)/abs(Ur)<=tol) exit !Check for convergence
    enddo    
    
    !Check solution
    !Ur2 = 16/3*m*K^2;
    !write(*,*) 'Ur = ',Ur,', Ur2 = ',Ur2
    
    c = sqrt(grav*h + Hgt/h*A) !Wave celerity
    B = ((3.0*em**2-5.0*em+2.0+(4.0*em-2.0)*E/fK)/3.0-(1.0-em-E/fK)**2)/em**2 !Shape parameter
    B = min(B,0.125) !Deep water limit

    return
    end subroutine cnoidal_sol
    
end module wave_lib