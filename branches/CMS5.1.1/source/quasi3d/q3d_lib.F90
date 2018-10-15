!==========================================================================================
module q3d_lib
! Quasi-3D module library (under development)
!
! Description:
!   The Quasi-3D module is divided into a variable definitions fortran module 
!   (q3d_def.f90), a library module (q3d_lib.F90), and CMS driver routines
!   (quasi3d.F90). In the quasi-3D method, the vertical structure of the current
!   velocity and sediment concentrations are approximated by assuming analytical
!   expressions, or by solving simple 1DV equations, and then either analytically
!   or numerically integrating the 1DV dispersion expersions. 
!
! Contains:
!  ~ Miscellaneous ~
!   q3d_flow_normal - Calculates the magnitude of the secondary flow velocity components
!                    normal to the stream-wise direction
!   ~ Profile Coefficients ~
!
!   ~ Dispersion ~
!     ~ Flow ~
!       q3d_flow_disp_surfquad - Calculates the dispersion terms in the surf zone zone
!                  based on quadratic velocity profiles
!       q3d_flow_disp_surflog - Calculates the dispersion terms in the surf zone zone 
!                 based on logarithmic velocity profiles
!       q3d_flow_disp_logloglin - Calculates the dispersion terms and wave-current
!               interaction termsin the case of logarithmic bottom friction and wind 
!                forcing profiles and linear secondary current profiles
!       q3d_flow_disp_powpowlin - Calculates the dispersion and wave-current
!               interaction terms for the case of powerlaw bottom friction and wind 
!               forcing profiles and linear secondary flow profiles
!     ~ Sediment ~
!       q3d_sed_surfquad - Calculates the sediment dispersion terms in the surf zone zone
!                 based on a quadratic velocity and exponential concentration profiles
!  ~ Profiles ~
!    ~ Flow ~
!      q3d_velpro_surfquad - Calculates the surf zone quadratic velocity profile
!      q3d_velpro_surflog - Calculates the the surf zone logarithmic velocity profile
!      q3d_velpro_logloglin - Calculates the dispersion terms in the for the  
!               case of bottom friction and wind forcing based on log profiles
!      q3d_velpro_powpowlin - Calculates the dispersion terms in the for the
!               case of bottom friction and wind forcing based on power-law profiles
!   ~ Sediment ~
!
! written by Alex Sanchez, USACE-CHL
!==========================================================================================
    implicit none
    
contains

!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
! Begin Profile Coefficients
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

!***********************************************************************************
    subroutine q3d_flow_coef_surfquad(hc,tausx,tausy,taubx,tauby,visvm,ax,ay,bx,by)
! Fits a quadratic surf zone based on a quadratic velocity profile
!
! Input:
!   hc - Local wave-averaged total water depth [m]
!   tausx,tausy - Surface stresses [Pa]
!   taubx,tauby - Bottom stresses [Pa]
!   visvm - Mean vertical diffusivity [m^2/s]
!
! Output:
!   ax,ay,bx,by - Flow quadratic profile coefficients
!
! written by Alex Sanchez, USACE-CHL
!***********************************************************************************
    use flow_def, only: rhow
    use prec_def
    implicit none
    !Input/Output
    real(ikind),intent(in) :: hc,tausx,tausy,taubx,tauby,visvm
    real(ikind),intent(out):: ax,ay,bx,by
    !Internal Variables
    real(ikind):: fac
    
    !Internal variables
    fac = hc/rhow/visvm
    ax = 0.5*fac*(tausx-taubx)
    ay = 0.5*fac*(tausy-tauby)
    bx = fac*taubx
    by = fac*tauby
    
    return
    endsubroutine q3d_flow_coef_surfquad

!********************************************************************
    subroutine q3d_flow_coef_surflog(hc,tausx,tausy,taubx,tauby,&
      visv0,visvslp,ax,ay,bx,by,dx,dy,c)
! Calculates the dxx dispersion terms in the 
! surf zone zone based on a quadratic profile
!
! Input:
!   hc - Local wave-averaged total water depth [m]
!   tausx,tausy - Surface stresses [Pa]
!   taubx,tauby - Bottom stresses [Pa]
!   visv0 - Bottom vertical diffusivity [m^2/s]
!   visvslp - Vertical diffusivity slope [m/s]
!   usx,usy - Stokes velocities
!
! Output:
!   fxx,fxy,fyy - Flow dispersion terms
!
! written by Alex Sanchez, USACE-CHL
!********************************************************************
    use flow_def, only: rhow
    use prec_def
    implicit none
    !Input/Output
    real(ikind),intent(in) :: hc,tausx,tausy,taubx,tauby,visv0,visvslp
    real(ikind),intent(out):: ax,ay,bx,by,dx,dy,c
    !Internal Variables
    real(ikind):: fac
    
    !Internal variables
    c = visvslp*hc/visv0
    fac = rhow*visvslp
    ax = (tausx-taubx)/fac
    ay = (tausy-tauby)/fac
    bx = taubx/fac-ax/c
    by = tauby/fac-ay/c
    fac = ((log(c+1.0)-1.0)*(c+1.0)+1.0)/c
    dx = ax*0.5+bx*fac
    dy = ay*0.5+by*fac       
    
    return
    endsubroutine q3d_flow_coef_surflog
      
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! End Profile Coefficients
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!  Begin Flow Dispersion
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<      

!*************************************************************************
    subroutine q3d_flow_disp_surfquad(hc,ax,ay,bx,by,fxx,fxy,fyy)
! Calculates the dispersion terms in the surf zone zone
! based on a quadratic velocity profile
!
! Input:
!   hc - Local wave-averaged total water depth [m]
!   ax,ay,bx,by - Profile coefficients
!   usx,usy - Stokes velocities
!
! Output:
!   fxx,fxy,fyy - Flow dispersion terms
!
! written by Alex Sanchez, USACE-CHL
!*************************************************************************
    use prec_def
    implicit none
    !Input/Output
    real(ikind),intent(in) :: hc,ax,ay,bx,by
    real(ikind),intent(out):: fxx,fxy,fyy
    !Internal Variables
    real(ikind), parameter :: aa=0.088888888889 !4/45
    real(ikind), parameter :: cc=0.083333333333 !1/12
    
    fxx = hc*(ax*ax*aa + (2.0*ax*bx+bx*bx)*cc)
    fxy = hc*(ax*ay*aa + (ax*bx+bx*ay+bx*by)*cc)
    fyy = hc*(ay*ay*aa + (2.0*ay*by+by*by)*cc)
    
    return
    endsubroutine q3d_flow_disp_surfquad 
    
!********************************************************************************
    subroutine q3d_flow_disp_surflog(hc,ax,ay,bx,by,dx,dy,c,fxx,fxy,fyy)
! Calculates the dxx dispersion terms in the 
! surf zone zone based on a quadratic profile
!
! Input:
!   hc - Local wave-averaged total water depth [m]
!   ax,ay,bx,by,dx,dy,c - Log profile coefficients
!
! Output:
!   fxx,fxy,fyy - Flow dispersion terms
!
! written by Alex Sanchez, USACE-CHL
!********************************************************************************
    use prec_def
    implicit none
    !Input/Output
    real(ikind),intent(in) :: hc,ax,ay,bx,by,dx,dy,c
    real(ikind),intent(out):: fxx,fxy,fyy
    !Internal Variables
    real(ikind):: ci,lncp1,cab,cbb,cbd
    
    !Temporary variables
    ci = 1.0/c
    lncp1 = log(c+1.0)
    
    !Dispersion terms
    cab = lncp1*(1.0-ci*ci)+ci-0.5
    cbb = lncp1*(lncp1*(1.0-ci)-2.0*ci-2.0)+2.0
    cbd = 1.0-lncp1*(1.0-ci)
    fxx = hc*(ax*ax/3.0 + ax*bx*cab - ax*dx &
        + bx*bx*cbb + 2.0*bx*dx*cbd + dx*dy)
    fxy = hc*(ax*ay/3.0 + 0.5*(ax*by+bx*ay)*cab - 0.5*(ax*dy+dx*ay) + bx*by*cbb &
        + (bx*dy+dx*by)*cbd + dx*dy)
    fyy = hc*(ax*ay/3.0 + ay*by*cab - ay*dy &
        + by*by*cbb + 2.0*by*dy*cbd)
    
    return
    endsubroutine q3d_flow_disp_surflog
    
!************************************************************************
    subroutine q3d_flow_disp_logloglin(hc,ucx,ucy,unx,uny,zap,&
       wx,wy,zwp,fxx,fxy,fyy)
! Calculates the dispersion terms in the for the case of 
! bottom friction and wind forcing based on log and 
! linear secondary flow profiles
!
! Input:
!   hc - Wave-averaged total water depth [m]
!   ucx,ucy - Current velocities (Eulerian) [m/s]
!   unx,uny - Magnitude velocities of secondary flow [m/s]
!   zap - Normalized apparent bottom roughness [-]
!   wx,wy - Lagrangian wind velocities * sqrt(densitair/rhow) [m/s]
!   zwp - Normalized apparent surface roughness [-]
!
! Output:
!   fxx,fxy,fyy - Dispersion terms [m^3/s]
!
! written by Alex Sanchez, USACE-CHL
!************************************************************************
    use prec_def
    implicit none
    !Input/Output
    real(ikind),intent(in) :: wx,wy,zwp,ucx,ucy,unx,uny,zap,hc
    real(ikind),intent(out):: fxx,fxy,fyy
    !Internal Variables
    real(ikind):: lnzap,lnzwp,c,cc,s,ss,bb,r
    real(ikind):: wrx,wry
    real(ikind), parameter:: pisq6m1 = 0.6449340668482     
        
    !Temporary variables
    r = sqrt(1.25/1025.0)
    wrx = r*wx
    wry = r*wy
    !zapc = min(zap,zp(1))
    lnzap = log(zap)
    lnzwp = log(zwp)
    c = lnzap+1.0
    s = lnzwp+1.0
    cc = c*c
    ss = s*s
    bb = pisq6m1/c/s
    
    !Dispersion terms plus wave-current interaction terms
    fxx = hc*(ucx*ucx/cc + wrx*wrx/ss + bb*2.0*ucx*wrx &
        + unx*ucx/c + unx*wrx/s + unx*unx/3.0)
    fxy = hc*(ucx*ucy/cc + wrx*wry/ss + bb*(ucx*wry+wrx*ucy) &
        + ((unx*ucy+ucx*uny)/c + (unx*wry+wrx*uny)/s)/2.0 + unx*uny/3.0)
    fyy = hc*(ucy*ucy/cc + wry*wry/ss + bb*2.0*ucx*wrx &
        + uny*ucy/c + uny*wry/s + uny*uny/3.0)
    
    return
    endsubroutine q3d_flow_disp_logloglin

!**********************************************************************************************
    subroutine q3d_flow_disp_powpowlin(hc,ucx,ucy,unx,uny,cm,wx,wy,cn,fxx,fxy,fyy)
! Calculates the dispersion and wave-current interaction terms for the case of 
! bottom friction and wind forcing based on power-law profiles
!
! Input:
!   hc - Wave-averaged total water depth [m]
!   ucx,ucy - Current velocities (Eulerian) [m/s]
!   unx,uny - Magnitude velocities of secondary flow [m/s]
!   cm - Bottom roughness parameter [-]
!   wx,wy - Lagrangian wind velocities [m/s]
!   cn - Surface roughness parameter [-]
!
! Output:
!   fxx,fxy,fyy - Dispersion terms [m^3/s]
!
! written by Alex Sanchez, USACE-CHL
!**********************************************************************************************
    use const_def, only: pi
    use fric_lib, only: fric_conv_nlength2pow
    use math_lib, only: beta
    use prec_def
    implicit none
    !Input/Output
    real(ikind),intent(in) :: wx,wy,cm,ucx,ucy,unx,uny,cn,hc
    real(ikind),intent(out):: fxx,fxy,fyy
    !Internal Variables
    real(ikind):: cc,ss,bb,ee,dd,cmn,xm,yn,r,wrx,wry
    
    !Temporary variables
    r = sqrt(1.25/1025.0)
    wrx = r*wx
    wry = r*wy   
    cc = cm*(cm+2.0)
    ss = cn*(cn+2.0)
    dd = 2.0*cm+1.0
    ee = 2.0*cn+1.0
    xm = 1.0/cm+1.0
    yn = 1.0/cn+1.0
    cmn = cm*cn
    bb = ((cm+cn+cmn+1.0)*beta(xm,yn)-cmn)/cmn
    
    !Dispersion terms plus wave-current interaction terms
    fxx = hc*(ucx*ucx/cc + wrx*wrx/ss &
        + 2.0*(ucx*wx*bb + unx*ucx/dd + unx*wrx/ee) + unx*unx/3.0)
    fxy = hc*(ucx*ucy/cc + wrx*wry/ss &
        + (ucx*uny+unx*ucy)*bb + (unx*ucy+ucx*uny)/dd &
        + (unx*wry+wrx*uny)/ee + unx*uny/3.0)
    fyy = hc*(ucy*ucy/cc + wy*wy/ss &
        + 2.0*(ucy*wy*bb + uny*ucy/dd + uny*wry/ee) + uny*uny/3.0)
    
    return
    endsubroutine q3d_flow_disp_powpowlin
    
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!  End Flow Dispersion
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>        
    
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!  Begin Sediment Dispersion
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<   

!**********************************************************************
    subroutine q3d_sed_surfquad(hc,ax,ay,bx,by,cak,wsk,visvk,sxk,syk)
! Calculates the sediment dispersion terms in the surf zone zone
! based on a quadratic velocity and exponential concentration profiles
!
! Input:
!   hc - Local wave-averaged total water depth [m]
!   ax,ay,bx,by - Coefficients for hydrodynamic quadratic profile
!   cak - Near-bed sediment concentration [kg/m^3]
!   wsk - Fall velocity for sediment size class [m/s]
!   visvk - Vertical sediment diffusivity [m^2/s]
!
! Output:
!   sxk,syk - Sediment dispersion
!
! written by Alex Sanchez, USACE-CHL
!**********************************************************************
    use prec_def
    implicit none
    !Input/Output
    real(ikind),intent(in) :: hc,ax,ay,bx,by
    real(ikind),intent(in) :: cak,wsk,visvk
    real(ikind),intent(out):: sxk,syk
    !Internal Variables
    real(ikind):: phik,em,ep,f1,f2
    
    !Temporary variables
    phik = hc*wsk/visvk
    ep = exp(phik)
    em = exp(-phik)
    f1 = 2.0*(ep-1.0)/phik**3 - 2.0/phik**2 - (2.0+ep)/(3.0*phik)
    f2 = (ep-1.0)/phik**2 - (1.0+ep)/(2.0*phik)      
    
    !Sediment dispersion
    sxk = cak*em*(ax*f1+bx*f2)
    syk = cak*em*(ay*f1+by*f2)
    
    return
    endsubroutine q3d_sed_surfquad    
    
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!  End Sediment Dispersion
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>        
    
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!  Begin Current Velocity Profiles
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<   

!*********************************************************************
    subroutine q3d_velpro_surfquad(nzp,zp,ucx,ucy,ax,ay,bx,by,uzp,vzp)
! Calculates the dxx dispersion terms in the 
! surf zone zone based on a quadratic profile
!
! Input:
!   nzp - # of vertical coordinates [-]
!   zp - Normalized vertical coordinates (0<zp<1) [-]
!   ucx,ucy - Current velocities (Eulerian) [m/s]
!   ax,ay,bx,by - Quadratic profile coefficients
!
! Output:
!   uzp,vzp - Current velocity profiles [m/s]
!
! written by Alex Sanchez, USACE-CHL
!*********************************************************************
    use prec_def
    implicit none
    !Intput
    integer    ,intent(in) :: nzp
    real(ikind),intent(in) :: zp(nzp)
    real(ikind),intent(in) :: ax,ay,bx,by,ucx,ucy
    !Output
    real(ikind),intent(out):: uzp(nzp),vzp(nzp)
    
    !Current velocity profiles
    uzp = ucx + ax*(zp*zp-0.333333333) + bx*(zp-0.5)
    vzp = ucy + ay*(zp*zp-0.333333333) + by*(zp-0.5)
    
    return
    endsubroutine q3d_velpro_surfquad

!*****************************************************************************
    subroutine q3d_velpro_surflog(nzp,zp,ucx,ucy,ax,ay,bx,by,dx,dy,c,uzp,vzp)
! Calculates the dxx dispersion terms in the 
! surf zone zone based on a quadratic profile
!
! Input:
!   nzp - # of vertical coordinates [-]
!   zp - Normalized vertical coordiantes [-]
!   ucx,ucy - Current velocities (Eulerian) [m/s]
!   ax,ay,bx,by,dx,dy,c - Log profile coefficients
!
! Output:
!   uzp,vzp - Current velocity profiles [m/s]
!
! written by Alex Sanchez, USACE-CHL
!*****************************************************************************
    use prec_def
    implicit none
    !Input
    integer    ,intent(in) :: nzp
    real(ikind),intent(in) :: zp(nzp)
    real(ikind),intent(in) :: ucx,ucy,ax,ay,bx,by,dx,dy,c
    !Output
    real(ikind),intent(out):: uzp(nzp),vzp(nzp)
    !Internal Variables
    real(ikind):: lnczp1(nzp)
    
    !Internal variables
    lnczp1 = log(c*zp+1.0)
    
    !Current velocity profile    
    uzp = ucx + ax*zp + bx*lnczp1 - dx
    vzp = ucy + ay*zp + by*lnczp1 - dy
    
    return
    endsubroutine q3d_velpro_surflog

!****************************************************************************
    subroutine q3d_velpro_logloglin(nzp,zp,ucx,ucy,unx,uny,zap,wx,wy,zwp,uzp,vzp)
! Calculates the dispersion terms in the for the case of 
! logarithmic bottom friction and wind forcing profiles and
! linear secondary flow profiles
!
! Input:
!   ucx,ucy - Current velocities (Eulerian) [m/s]
!   unx,uny - Magnitude velocities of secondary flow [m/s]
!   zap - Normalized apparent bottom roughness [-]
!   wx,wy - Lagrangian wind velocities [m/s]
!   zwp - Normalized apparent surface roughness [-]
!
! Output:
!   uzp,vzp - Current velocity profiles [m/s]
!
! written by Alex Sanchez, USACE-CHL
!****************************************************************************
    use prec_def
    implicit none
    !Input/Output
    integer    ,intent(in) :: nzp
    real(ikind),intent(in) :: zp(nzp)
    real(ikind),intent(in) :: ucx,ucy,unx,uny,zap,wx,wy,zwp
    !Output
    real(ikind),intent(out):: uzp(nzp),vzp(nzp)
    !Internal Variables
    real(ikind):: pro(nzp),prow(nzp),pron(nzp),zapc,r
    
    !Internal variabes
    zapc = min(zap,zp(1))
    pro  = log(zp/zapc)/(log(1.0/zapc)-1.0)-1.0    
    prow = log((1.0-zp)/zwp)/(1.0-log(1.0/zwp))+1.0
    pron = 2.0*zp-1.0
    r = sqrt(1.25/1000.0)
    
    !Profiles
    uzp = ucx + ucx*pro + r*wx*prow + unx*pron
    vzp = ucy + ucy*pro + r*wy*prow + uny*pron
    
    return
    endsubroutine q3d_velpro_logloglin
    
!****************************************************************************
    subroutine q3d_velpro_powpowlin(nzp,zp,ucx,ucy,unx,uny,cm,wx,wy,cn,uzp,vzp)
! Calculates the dispersion and wave-current
! interaction terms for the case of powerlaw bottom friction and wind 
! forcing profiles and linear secondary flow profiles
!
! Input:
!   ucx,ucy - Current velocities (Eulerian) [m/s]
!   unx,uny - Magnitude velocities of secondary flow [m/s]
!   cm - Bottom roughness parameter [-]
!   wx,wy - Lagrangian wind velocities [m/s]
!   cn - Surface roughness parameter [-]
!
! Output:
!   uzp,vzp - Current velocity profiles [m/s]
!
! written by Alex Sanchez, USACE-CHL
!****************************************************************************
    
    use prec_def
    implicit none
    !Input
    integer    ,intent(in) :: nzp
    real(ikind),intent(in) :: zp(nzp)
    real(ikind),intent(in) :: ucx,ucy,unx,uny,cm,wx,wy,cn
    !Output
    real(ikind),intent(out):: uzp(nzp),vzp(nzp)
    !Internal Variables
    real(ikind):: pro(nzp),prow(nzp),pron(nzp),r
    
    !Internal variables
    pro = (cm+1.0)/cm*zp**(1.0/cm)-1.0
    prow = (cn+1.0)/cn*zp**(1.0/cn)-1.0
    pron = 2.0*zp-1.0
    r = sqrt(1.25/1000.0)
    
    !Profiles
    uzp = ucx + ucx*pro + r*wx*prow + unx*pron
    vzp = ucy + unx*pro + r*wy*prow + uny*pron
    
    return
    endsubroutine q3d_velpro_powpowlin
    
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!  End Current Velocity Profiles 
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>        

!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!  Begin Current Velocity Deviation at Surface
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

!*************************************************************************
    subroutine q3d_flow_uds_surfquad(ax,ay,bx,by,udsx,udsy)
! Calculates the current velocity deviation at the surface in the surf zone
! based on a quadratic velocity profile
!
! Input:
!   hc - Local wave-averaged total water depth [m]
!   ax,ay,bx,by - Profile coefficients
!   usx,usy - Stokes velocities [m/s]
!
! Output:
!   udsx,udsy - Current velocity deviation at surface [m/s]
!
! written by Alex Sanchez, USACE-CHL
!*************************************************************************
    use prec_def
    implicit none
    !Input/Output
    real(ikind),intent(in) :: ax,ay,bx,by
    real(ikind),intent(out):: udsx,udsy
    
    !Current deviation at the surface
    udsx = ax*0.6666666667 + bx*0.5
    udsy = ay*0.6666666667 + by*0.5
    
    return
    endsubroutine q3d_flow_uds_surfquad 

!********************************************************************************
    subroutine q3d_flow_uds_surflog(ax,ay,bx,by,dx,dy,c,udsx,udsy)
! Calculates the current velocity deviation at the surface in the 
! surf zone zone based on a log profile
!
! Input:
!   ax,ay,bx,by,dx,dy,c - Log profile coefficients
!   usx,usy - Stokes velocities [m/s]
!
! Output:
!   udsx,udsy - Current velocity deviation at the surface [m/s]
!
! written by Alex Sanchez, USACE-CHL
!********************************************************************************
    use prec_def
    implicit none
    !Input/Output
    real(ikind),intent(in) :: ax,ay,bx,by,dx,dy,c
    real(ikind),intent(out):: udsx,udsy
    !Internal Variables
    real(ikind):: lncp1
    
    !Temporary variables
    lncp1 = log(c+1.0)
    
    !Current deviation at the surface
    udsx = ax + bx*lncp1 - dx
    udsy = ay + by*lncp1 - dy
    
    return
    endsubroutine q3d_flow_uds_surflog
    
!************************************************************************
    subroutine q3d_flow_uds_logloglin(ucx,ucy,unx,uny,zap,wx,wy,zwp,udsx,udsy)
! Calculates the current velocity deviation at the surface for the case of 
! bottom friction and wind forcing based on log and 
! linear secondary flow profiles
!
! Input:
!   ucx,ucy - Current velocities (Eulerian) [m/s]
!   unx,uny - Magnitude velocities of secondary flow [m/s]
!   zap - Normalized apparent bottom roughness [-]
!   wx,wy - Lagrangian wind velocities [m/s]
!   zwp - Normalized apparent surface roughness [-]
!
! Output:
!   udsx,udsy - Current velocity deviation at the surface [m/s]
!
! written by Alex Sanchez, USACE-CHL
!************************************************************************
    use prec_def
    implicit none
    !Input/Output
    real(ikind),intent(in) :: ucx,ucy,unx,uny,zap,wx,wy,zwp
    real(ikind),intent(out):: udsx,udsy
    !Internal Variables
    real(ikind):: lnzap,lnzwp,proc,prow,r
        
    !Temporary variables
    r = sqrt(1.25/1025.0)
    !zapc = min(zap,zp(1))
    lnzap = log(zap)
    lnzwp = log(zwp)
    
    !Current deviation at the surface    
    proc = lnzap/(lnzap+1.0)-1.0
    prow = (log(1.0-0.99999)-lnzwp)/(lnzwp+1.0)+1.0 !Note: wind component evaluated at zp=0.999999
    udsx = ucx*proc + r*wx*prow + unx
    udsy = ucy*proc + r*wy*prow + uny
    
    return
    endsubroutine q3d_flow_uds_logloglin    
       
!**********************************************************************************************
    subroutine q3d_flow_uds_powpowlin(ucx,ucy,unx,uny,cm,wx,wy,udsx,udsy)
! Calculates the current velocity deviation at the surface for the case of 
! bottom friction and wind forcing based on power-law profiles
!
! Input:
!   ucx,ucy - Current velocities (Eulerian) [m/s]
!   unx,uny - Magnitude velocities of secondary flow [m/s]
!   cm - Bottom roughness parameter [-]
!   wx,wy - Lagrangian wind velocities [m/s]
!
! Output:
!   udsx,udsy - Current velocity deviation at the surface [m/s]
!
! written by Alex Sanchez, USACE-CHL
!**********************************************************************************************
    use prec_def
    implicit none
    !Input/Output
    real(ikind),intent(in) :: ucx,ucy,unx,uny,cm,wx,wy
    real(ikind),intent(out):: udsx,udsy
    !Internal Variables
    real(ikind):: r
    
    !Temporary variables
    r = sqrt(1.25/1025.0)
    
    !Current deviation at the surface
    udsx = ucx/cm + r*wx + unx
    udsy = ucy/cm + r*wy + uny
    
    return
    endsubroutine q3d_flow_uds_powpowlin       

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! End Current Velocity Deviation at Surface
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!  Begin Wave Weighted Effective Current Velocity Deviation
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<    

!*************************************************************************
    subroutine q3d_flow_ude_surfquad(hc,ax,ay,bx,by,kw,udex,udey)
! Calculates the current velocity deviation component 
!
! Input:
!   hc - Local wave-averaged total water depth [m]
!   ax,ay,bx,by - Profile coefficients
!   usx,usy - Stokes velocities [m/s]
!
! Output:
!   udex,udey - Wave weighted effective current velocity deviation [m/s]
!
! written by Alex Sanchez, USACE-CHL
!*************************************************************************
    use prec_def
    implicit none
    !Input/Output
    real(ikind),intent(in) :: hc,ax,ay,bx,by,kw
    real(ikind),intent(out):: udex,udey
    !Internal variables
    real(ikind) :: twokh,cosh2kh,twokhsinh2kh
    
    twokh = 2.0*kw*hc
    twokh = min(twokh,20.0) !Limit to avoid NaNs
    cosh2kh = cosh(twokh)
    twokhsinh2kh = twokh*sinh(twokh)
    
    !Current deviation at the surface
    udex = ax*2.66666667 + bx*0.5 - (2.0*ax*cosh2kh+bx*(cosh2kh-1.0))/twokhsinh2kh !8/3=2.66666667
    udey = ay*2.66666667 + by*0.5 - (2.0*ay*cosh2kh+by*(cosh2kh-1.0))/twokhsinh2kh
    
    return
    endsubroutine q3d_flow_ude_surfquad 
    
!********************************************************************************
    subroutine q3d_flow_ude_surflog(hc,ax,ay,bx,by,dx,dy,c,kw,udex,udey)
! Calculates the current velocity deviation at the surface in the 
! surf zone zone based on a log profile
!
! Input:
!   ax,ay,bx,by,dx,dy,c - Log profile coefficients
!   usx,usy - Stokes velocities [m/s]
!
! Output:
!   udsx,udsy - Current velocity deviation at the surface [m/s]
!
! written by Alex Sanchez, USACE-CHL
!********************************************************************************
    use prec_def
    implicit none
    !Input/Output
    real(ikind),intent(in) :: hc,ax,ay,bx,by,dx,dy,c,kw
    real(ikind),intent(out):: udex,udey
    !Internal Variables
    real(ikind):: lnchp1,twokh,cosh2kh,twokhsinh2kh,cc0,cc1,cc2,cc3
        
    !Temporary variables
    lnchp1 = log(c/2.0+1.0)
    twokh = 2.0*kw*hc
    cosh2kh = cosh(twokh)
    twokhsinh2kh = twokh*sinh(twokh)
    
    !Current deviation at the surface
    cc0 = c/(c+2.0)
    cc1 = 2.0*c*c*(twokh**2/4.0+2.0)/(twokh**2*(c+2.0)**2)
    cc2 = (cosh2kh-1.0)/twokhsinh2kh
    cc3 = (cosh2kh+1.0)/twokhsinh2kh
    udex = ax + bx*lnchp1 - dx + bx*cc0 - bx*cc1 &
         - (ax+2.0*bx*cc0)*cc2 + (2.0*bx*c*c+2.0)/(c+2.0)**2*cc3
    udey = ay + by*lnchp1 - dy + by*cc0 - by*cc1 &
         - (ay+2.0*by*cc0)*cc2 + (2.0*by*c*c+2.0)/(c+2.0)**2*cc3
    
    return
    endsubroutine q3d_flow_ude_surflog
    
!********************************************************************************
    subroutine q3d_flow_ude_logloglin(hc,ucx,ucy,unx,uny,zap,&
       wx,wy,zwp,kw,udex,udey)
! Calculates the current velocity deviation at the surface in the 
! surf zone zone based on a log profile
!
! Input:
!   ax,ay,bx,by,dx,dy,c - Log profile coefficients
!   usx,usy - Stokes velocities [m/s]
!
! Output:
!   udsx,udsy - Current velocity deviation at the surface [m/s]
!
! written by Alex Sanchez, USACE-CHL
!********************************************************************************
    use math_lib, only: e1x
    use prec_def
    implicit none
    !Input/Output
    real(ikind),intent(in) :: hc,ucx,ucy,unx,uny,zap,wx,wy,zwp,kw
    real(ikind),intent(out):: udex,udey
    !Internal Variables
    real(ikind):: twokh,cc1,cc2,cc3
    real(ikind):: logzap,kh,r
        
    !Temporary variables
    kh = kw*hc
    kh = min(kh,10.0) !Limit to avoid NaNs
    twokh = 2.0*kh    
    logzap = log(zap)
    
    r = sqrt(1.25/1000.0)
    
    !Current deviation at the surface
    cc1 = (logzap+(e1x(twokh)-e1x(-twokh))/2.0/sinh(twokh))/(logzap+1.0)
    cc2 = r*kh*kh/(-1.0-log(zwp))*(0.5/cosh(kh)+0.0556*kh/sinh(kh))
    cc3 = 1.0-tanh(kh)/kh
    udex = ucx*cc1 + wx*cc2 + unx*cc3
    udey = ucy*cc1 + wy*cc2 + uny*cc3
    
    return
    endsubroutine q3d_flow_ude_logloglin       
    
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! End Wave Weighted Effective Current Velocity Deviation
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!  Begin Wave-current Interaction Term
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
!*************************************************************************
    subroutine q3d_flow_wavcur(hc,udsx,udsy,usx,usy,fxx,fxy,fyy)
! Calculates the wave-current interaction terms
!
! Input:
!   hc - Local wave-averaged total water depth [m]
!   udsx,udsy - Current velocity deviation at the surface [m/s]
!   usx,usy - Stokes velocities [m/s]
!
! Output:
!   fxx,fxy,fyy - Terms with added wave-current interaction
!
! written by Alex Sanchez, USACE-CHL
!*************************************************************************
    use prec_def
    use q3d_def, only: facwci
    implicit none
    !Input/Output
    real(ikind),intent(in) :: hc,udsx,udsy,usx,usy
    real(ikind),intent(out):: fxx,fxy,fyy
    
    !Wave-current interaction terms
    fxx = fxx - facwci*hc*2.0*usx*udsx
    fxy = fxy - facwci*hc*(usy*udsx + usx*udsy)
    fyy = fyy - facwci*hc*2.0*usy*udsy
    
    return
    endsubroutine q3d_flow_wavcur 

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! End Wave-current Interaction Term 
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!  Begin Secondary Flow
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

!*****************************************************************
    subroutine q3d_flow_normal(hc,uc,vc,fcor,zap,Kc,unx,uny)
! Calculates the magnitude of the secondary flow velocity components
! normal to the stream-wise direction
!
! written by Alex, Sanchez, USACE-CHL
!*****************************************************************
    use prec_def
    implicit none
    !Input/Output
    real(ikind),intent(in) :: uc,vc,hc,fcor,zap,Kc
    real(ikind),intent(out) :: unx,uny
    !Internal Variables
    real(ikind) :: bs  
    
    bs = 6.25*(1.5-3.0/(log(1.0/zap)-1.0)) !Approximately equal to 6
    !unx = (vc/max(abs(vc),1.0e-6)*fcor*0.5 + abs(vc)*Kc)*hc*bs
    !uny = (-uc/max(abs(uc),1.0e-6)*fcor*0.5 - abs(uc)*Kc)*hc*bs
    
    unx = (sign(1.0,vc)*fcor*0.5 + abs(vc)*Kc)*hc*bs
    uny = (-sign(1.0,uc)*fcor*0.5 - abs(uc)*Kc)*hc*bs
    
    return
    endsubroutine q3d_flow_normal
    
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!  End Secondary Flow    
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>    
    
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!  Begin Vertical Eddy Viscosity
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

!*************************************************************************
    function q3d_eddyvert_mean(hc,ustarb,ustars) result(visvm)
! Calculates the mean value of the vertical eddy viscosity
!
! Input:
!   hc - Wave-averaged total water depth [m]
!   vish - Horizontal eddy viscosity [m^2/s]
!   ustarb - Bottom shear velocity [m/s]
!   ustars - Surface shear velocity [m/s]
!
! Output:
!   visvm - Mean value of the vertical eddy viscosity [m^2/s]
!
! written by Alex Sanchez, USACE-CHL
!*************************************************************************    
    use q3d_def, only: cvisvv,cvisvhc,visvcon
    use flow_def, only: viscos,grav
    use prec_def
    implicit none
    real(ikind),intent(in) :: hc,ustarb,ustars
    real(ikind) :: visvm
    
    visvm = visvcon + cvisvv*(ustarb+ustars)*hc + cvisvhc*hc*sqrt(grav*hc)
    
    return
    endfunction q3d_eddyvert_mean
    
!*********************************************************************************
    subroutine q3d_eddyvert_bottom_slope(hc,Hsig,visvm,visv0,visvslp)
! Calculates the bottom value and slope of the vertical eddy viscosity
! based on the local water depth, mean vertical viscosity, and swave height
! 
! Input:
!   hc - Wave-averaged total water depth [m]
!   Hsig - Significant wave height [m]
!   visvm - Mean vertical eddy viscosity [m^2/s]
!
! Output:
!   visv0 - Bottom value of the vertical eddy viscosity (visv0>0) [m^2/s]
!   visvslop - Vertical slope of the vertical eddy viscosity (visvslop>0) [m^2/s]
!
! written by Alex Sanchez, USACE-CHL
!*********************************************************************************    
    use q3d_def, only: cvisvv,cvisvhc,visvslp0,cvisvslpmu
    use flow_def, only: viscos,grav
    use prec_def
    implicit none
    real(ikind),intent(in) :: hc,Hsig,visvm
    real(ikind),intent(out) :: visv0,visvslp
    real(ikind) :: gam,visvslpmax
    
    gam = 0.7071*Hsig/hc
    visvslp = visvslp0 + cvisvslpmu*sqrt(grav*hc)/gam !Based on Steetzel
    visvslpmax = 2.0*(visvm-viscos)/hc !Maximum slope value
    visvslp = min(visvslp,visvslpmax)  !Slope (limited to avoid negative values)
    visv0 = visvm-visvslp*hc/2.0   !Bottom viscosity
    visv0 = max(visv0,1.0e-6) !Avoid divide by zero
    
    return
    endsubroutine q3d_eddyvert_bottom_slope

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!  End Vertical Eddy Viscosity
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>    

endmodule q3d_lib