!===================================================================    
module comp_lib
! Computational module
!  
!   Matrix coefficients for various advection-diffusion schemes
!     zerocoef - No advection-diffusion
!     upwindcoef - Upwind scheme
!     hybridcoef - Hybrid scheme
!     exponentialcoef - Exponential scheme
!
!  Deferred corrections
!    gammadefcor - Gamma family schemes
!    cubistadefcor - CUBISTA scheme
!    alvsmartdefcor - ALVSMART scheme
!    hoabdefcor - HOAB scheme
!
!  Flux limiters
!    superbeefluxlim - Superbee scheme
!    minmodfluxlim - Minmod scheme
!    osherfluxlim - Osher scheme
!    musclfluxlim - MUSCL scheme
!    charmfluxlim - CHARM scheme
!    hcusfluxlim - HCUS scheme
!    hquickfluxlim - HQUICK scheme
!    osprefluxlim - OSPRE scheme
!    korenfluxlim - KOREN scheme
!    vanalbadafluxlim - van Albada scheme
!    smartfluxlim - SMART scheme
!
!  Limited slopes
!
! written by Alex Sanchez, USACE-CHL
!===================================================================    
    use prec_def
    implicit none
    
contains

!************************************************************
    function zerocoef(dk,fk) result(acoefik)
! Calculates a matrix coefficient for no advection-diffusion
!************************************************************
    implicit none
    real(ikind),intent(in) :: dk,fk
    real(ikind) :: acoefik,rdum
    
    
    rdum=dk                                  !This is stupid.  We should remove the inbound arguments.
    rdum=fk
    acoefik=0.0
    
    end function zerocoef

!**************************************************
    function upwindcoef(dk,fk) result(acoefik)
! Calculates a matrix coefficient for the 
! Upwind advection-diffusion scheme
!**************************************************
    implicit none
    real(ikind),intent(in) :: dk,fk
    real(ikind) :: acoefik
    
    acoefik=dk+max(0.0,-fk)
    
    end function upwindcoef
   
!**************************************************
    function hybridcoef(dk,fk) result(acoefik)
! Calculates a matrix coefficient for the 
! Hybrid advection-diffusion scheme
!**************************************************
    implicit none
    real(ikind),intent(in) :: dk,fk
    real(ikind) :: peck,acoefik
    
    peck=abs(fk)/max(dk,1.0e-6)
    acoefik=dk*max(0.0,1.0-0.5*peck)+max(0.0,-fk)
    
    end function hybridcoef

!*****************************************************
    function powerlawcoef(dk,fk) result(acoefik)
! Calculates a matrix coefficient for the     
! Power-Law advection-diffusion scheme
!*****************************************************
    implicit none
    real(ikind),intent(in) :: dk,fk
    real(ikind) :: peck,acoefik,absfk
    
    absfk=abs(fk)
    peck=absfk/max(dk,absfk,1.0e-6)  
    acoefik=dk*(1.0-peck**5)+max(0.0,-fk)
    
    end function powerlawcoef
    
!********************************************************
    function exponentialcoef(dk,fk) result(acoefik)
! Calculates a matrix coefficient for the     
! Exponential advection scheme
!********************************************************
    implicit none
    real(ikind),intent(in) :: dk,fk
    real(ikind) :: peck,acoefik,absfk
    
    !!peck=max(min(abs(fk)/max(dk,1.0e-6),30.0),1.0e-6)
    absfk=abs(fk)
    peck=absfk/max(max(dk,1.e-6),absfk)
    acoefik=dk*peck/max(exp(peck)-1.0,1.0e-6)+max(0.0,-fk)
    
    end function exponentialcoef

!*************************************************************
    function gammadefcor(phin) result(gamma)
! Gamma Deferred Correction
! Reference:
!   Jasak H, Weller H, and Gosman A.D. (1999)  
!     High resolution NVD differencing scheme for arbitrarily
!     unstructured meshes. International Journal for
!     Numerical Methods in Fluids, 31, 431–449.
!*************************************************************
    implicit none
    real(ikind),intent(in) :: phin
    real(ikind) :: gamma
    
    if(phin>0.0 .and. phin<=0.5)then
      gamma=2.0*phin
    elseif(phin>0.5 .and. phin<=1.0)then
      gamma=1.0
    else
      gamma=0.0
    endif
    
    end function gammadefcor

!**********************************************************
    function cubistadefcor(phin) result(gamma)
! CUBISTA Deferred Correction
! Reference:
!   Alves M.A., Oliveira P.J., and Pinho F.T. (2003). 
!     A convergent and universally bounded interpolation 
!     scheme for the treatment of advection. International
!     Journal for Numerical Methods in Fluids, 41, 47–75.
!**********************************************************
    implicit none
    real(ikind),intent(in) :: phin
    real(ikind) :: gamma
    
    if(phin>0.0 .and. phin<=0.375)then
      gamma=1.5*phin/(1.0-phin)
    elseif(phin>0.375 .and. phin<=0.75)then
      gamma=(3.0-2.0*phin)/(4.0-4.0*phin)
    elseif(phin>0.75 .and. phin<=1.0)then
      gamma=1.5
    else
      gamma=0.0
    endif
    
    end function cubistadefcor

!**************************************************
    function alvsmartdefcor(phin) result(gamma)
! ALVSMART Deferred Correction
! Reference:
!   Przulj V., and Basara B. (2001) Bounded convection 
!     schemes for unstructured grids. 
!     AIAA Paper No. 2001-2593, AIAA Computational Fluid 
!     Dynamics Conference, Anaheim, California, U.S.A.
!*******************************************************
    implicit none
    real(ikind),intent(in) :: phin
    real(ikind) :: gamma
    
    if(phin>0.0 .and. phin<=0.25)then
      gamma=2.5*phin/(1.0-phin)
    elseif(phin>0.25 .and. phin<=0.75)then
      gamma=(3.0-2.0*phin)/(4.0-4.0*phin)
    elseif(phin>0.75 .and. phin<=1.0)then
      gamma=1.5
    else
      gamma=0.0
    endif
    
    end function alvsmartdefcor

!********************************************************************
    function hoabdefcor(phin) result(gamma)
! HOAB Deferred Correction
! Reference:
!   Wei J.J., Yu B., Tao W.Q., Kawaguchi Y., and Wang H.S. (2003). 
!     A new high-order accurate and bounded scheme for incompressible flow. 
!    Numerical Heat Transfer, Part B, 43, 19–41.
!********************************************************************
    implicit none
    real(ikind),intent(in) :: phin
    real(ikind) :: gamma
    
    if(phin>0.0 .and. phin<=0.1667)then
      gamma=5.0*phin/(1.0-phin)
    elseif(phin>0.1667 .and. phin<=0.5)then
      gamma=1.0
    elseif(phin>0.5 .and. phin<=0.75)then
      gamma=0.5/(1.0-phin)
    elseif(phin>0.75 .and. phin<=1.0)then
      gamma=2.0
    else
      gamma=0.0
    endif
    
    end function hoabdefcor
    
!*****************************************************************
    function superbeefluxlim(r) result(phi)
! Superbee flux limiter
! Reference:
!   Roe, P.L., (1986). Characteristic-based schemes for the Euler equations, 
!     Ann. Rev. Fluid Mech., 18, p337.
!*****************************************************************
    implicit none
    real(ikind),intent(in) :: r
    real(ikind) :: phi
    
    phi=max(0.0,min(1.0,2.0*r),min(2.0,r))
    
    end function superbeefluxlim
    
!*****************************************************************
    function minmodfluxlim(r) result(psi)
! Minmod flux limiter
! Reference:
!   Roe, P L, (1986). Characteristic-based schemes for the Euler equations, 
!     Ann. Rev. Fluid Mech., 18, p337.
!*****************************************************************
    implicit none
    real(ikind),intent(in) :: r
    real(ikind) :: psi
    
    psi=max(0.0,min(1.0,r))
    
    end function minmodfluxlim
    
!*****************************************************************
    function osherfluxlim(r) result(psi)
! Osher flux limiter
! Reference:
!   Chakravarthy, S.R. and Osher, S., (1983). High resolution applications of the 
!     Osher upwind scheme for the Euler equations, AIAA Paper 83-1943, 
!     Proc. AIAA 6th Computational Fluid Dynamics Conference, pp 363–73.
!*****************************************************************
    implicit none
    real(ikind),intent(in) :: r
    real(ikind) :: psi
    
    psi=max(0.0,min(2.0,r))
    
    end function osherfluxlim
    
!*********************************************************************
    function musclfluxlim(r) result(psi)
! MUSCL flux limiter
! Reference:
!   Van Leer, B. (1974). Towards the ultimate conservative 
!     difference scheme II. Monotonicity and conservation 
!     combined in a second order scheme. J. Comp. Phys., 14, p361-70.
!*********************************************************************
    implicit none
    real(ikind),intent(in) :: r
    real(ikind) :: psi,absr
    
    absr=abs(r)
    psi=(r+absr)/(1.0+absr)
    
    end function musclfluxlim
    
!*************************************************************************
    function charmfluxlim(r) result(psi)
! CHARM flux limiter
! Reference:
!   Zhou, G. (1995). Numerical simulations of physical discontinuities 
!     in single and multi-fluid flows for arbitrary Mach numbers, PhD Thesis, 
!     Chalmers Univ. of Tech., Goteborg, Sweden.
!*************************************************************************
    implicit none
    real(ikind),intent(in) :: r
    real(ikind) :: psi
    
    if(r>0.0)then
      psi=r*(3.0*r+1.0)/(r+1.0)**2
    else
      psi=0.0
    endif
    
    end function charmfluxlim

!********************************************************************************
    function hcusfluxlim(r) result(psi)
! HCUS flux limiter
! Waterson, N.P., and Deconinck, H. (1995). A unified approach to the design and
! application of bounded higher-order convection schemes, VKI Preprint 1995-21.
!********************************************************************************
    implicit none
    real(ikind),intent(in) :: r
    real(ikind) :: psi
    
    psi=1.5*(r+abs(r))/(r+2.0)
    
    end function hcusfluxlim        

!********************************************************************************
    function hquickfluxlim(r) result(psi)
! HQUICK flux limiter
! Waterson, N.P., and Deconinck, H. (1995). A unified approach to the design and
! application of bounded higher-order convection schemes, VKI Preprint 1995-21.
!********************************************************************************
    implicit none
    real(ikind),intent(in) :: r
    real(ikind) :: psi
    
    psi=2.0*(r+abs(r))/(r+3.0)
    
    end function hquickfluxlim
    
!********************************************************************************
    function osprefluxlim(r) result(psi)
! OSPRE flux limiter
! Waterson, N.P., and Deconinck, H. (1995). A unified approach to the design and
! application of bounded higher-order convection schemes, VKI Preprint 1995-21.
!********************************************************************************
    implicit none
    real(ikind),intent(in) :: r
    real(ikind) :: psi
    
    psi=1.5*(r*r+r)/(r*r+r+1.0)
    
    end function osprefluxlim
    
!************************************************************************************
    function korenfluxlim(r) result(psi)
! KOREN flux limiter
! Koren, B. (1993). A robust upwind discretisation method for advection, 
! diffusion and source terms, In: Numerical Methods for Advection-Diffusion Problems, 
! Ed. C.B.Vreugdenhil & B.Koren, Vieweg, Braunschweig, p117.
!************************************************************************************
    implicit none
    real(ikind),intent(in) :: r
    real(ikind) :: psi
    
    psi=max(0.0,min(2.0*r,0.6666667+r*0.3333333,2.0))
    
    end function korenfluxlim
    
!************************************************************************************
    function vanalbadafluxlim(r) result(psi)
! van Albada flux limiter
! Van Albada, G.D., Van Leer, B., and Roberts, W.W. (1982). A comparative study
! of computational methods in cosmic gas dynamics, Astron. Astrophysics, 108, p76.
!************************************************************************************
    implicit none
    real(ikind),intent(in) :: r
    real(ikind) :: psi
    
    psi=(r*r+r)/(r*r+1.0)
    
    end function vanalbadafluxlim

!*************************************************************************
    function smartfluxlim(r) result(psi)
! smart flux limiter
! Gaskell, P.H., and Lau, A.K.C. (1988). Curvature-compensated convective 
! transport: SMART, a new boundedness-preserving transport algorithm, 
! Int. J. Num. Meth. Fluids, 8, p617.
!*************************************************************************
    implicit none
    real(ikind),intent(in) :: r
    real(ikind) :: psi
    
    psi=max(0.0,min(2.0*r,0.25+0.75*r,4.0))
    
    end function smartfluxlim
    
!*************************************************************************
    function mimodslope(s1,s2) result(s)
! Minmod limited slope
!*************************************************************************
    implicit none
    real(ikind),intent(in) :: s1,s2
    real(ikind) :: s
    
    if(s1*s2<1.e-15)then
      s=0.0
    else
      s=sign(1.0,s1)*min(abs(s1),abs(s2))
    endif
    
    end function mimodslope    

!!*************************************************************************
!    function superbeeslope(s1,s2) result(s)
!! Superbee limited slope
!!*************************************************************************
!    use comvarbl, only: nbeta
!    implicit none
!    real(ikind),intent(in) :: s1,s2
!    real(ikind) :: s,a,b
!    
!    if(s1*s2<1.e-15)then
!      s=0.0
!    else
!      a=abs(s1);  b=abs(s2)
!      s=sign(1.0,s1)*min(max(a,b),nbeta*min(a,b))
!    endif 
!    
!    end function superbeeslope       
    
!*************************************************************************
    function frommslope(s1,s2) result(s)
! Fromm limited slope
!*************************************************************************
    implicit none
    real(ikind),intent(in) :: s1,s2
    real(ikind) :: s
    
    s=0.5*(s1+s2)
    
    end function frommslope       
    
!*************************************************************************
    function vanleerslope(s1,s2) result(s)
! Van Leer limited slope
!*************************************************************************
    implicit none
    real(ikind),intent(in) :: s1,s2
    real(ikind) :: s
    
    if(s1*s2<1.0e-15)then
      s=0.0
    else
      s=2*s1*s2/(s1+s2)
    endif
    
    end function vanleerslope      
      
!*************************************************************************
    function vanalbadaslope(s1,s2) result(s)
! Van Albada limited slope
!*************************************************************************
    implicit none
    real(ikind),intent(in) :: s1,s2
    real(ikind), parameter :: eps=1.0e-15
    real(ikind) :: s
        
    s=(s1*(s2*s2+eps)+s2*(s1*s1+eps))/(s1*s1+s2*s2+2*eps)
    
    end function vanalbadaslope       

!*************************************************************************
    function dbleminmodslope(s1,s2) result(s)
! Double Minmod limited slope
!*************************************************************************
    implicit none
    real(ikind),intent(in) :: s1,s2
    real(ikind) :: s,a,b,c
        
    if(s1*s2<1.0e-15)then
      s=0.0
    else
      a=abs(s1);  b=abs(s2); c=0.5*(a+b)
      s=sign(1.0,s1)*min(2*a,2*b,c)
    endif
    
    end function dbleminmodslope     
    
!*************************************************************
    subroutine hll(hl,hr,el,er,ul,ur,vl,vr,cn,sn,F)
! Harten, Lax, and van Leer (1983) approximate Riemann solver
!
! References:
!  Harten, A., Lax, P.D., van Leer, B. 1983. On Upstream 
!    Differencing and Godunov-Type Schemes for Hyperbolic 
!    Conservation Laws. SIAM Review, 25(1), 35-61.
!  Ying, X.Y., Wang, S.S.Y. 2008. Improved implementation 
!    of the HLL approximate Riemann solver for one dimensional 
!    open channel flows. Journal of Hydraulic Research. 
!    46(1), 21–34.
!
! Author: Alex Sanchez, USACE-CHL
!*************************************************************
    use flow_def, only: hdry,grav
    implicit none

    !Input
    real(ikind),intent(in) :: hl,hr,el,er,ul,ur,vl,vr,cn,sn

    !Output
    real(ikind),intent(out) :: F(3)

    !Internal
    real(ikind) :: FL(3),FR(3),PL(3),PR(3),hlplushr
    real(ikind) :: sl,sr,hs,cr,cl,qpl,qpr,upl,upr !,hgl,hgr,cs,us,amax

    upl=ul*cn+vl*sn;  upr=ur*cn+vr*sn
    !hl=max(hl,hmin);  hr=max(hr,hmin)
    cl=sqrt(grav*hl); cr=sqrt(grav*hr)
    qpl=upl*hl;       qpr=upr*hr
        
    !!Rusanov (1961) flux for smooth regions
    !!if(abs(hl-hr)/max(hl,hr)<0.001)then
    !  !hgl=0.5_dp*grav*hl*hl; hgr=0.5_dp*grav*hr*hr
    !  FL=(/qpl, qpl*ul, qpl*vl/)
    !  FR=(/qpr, qpr*ur, qpr*vr/)
    !  PL=(/el, qpl*cn, qpl*sn/)
    !  PR=(/er, qpr*cn, qpr*sn/)   
    !  amax=max(abs(upl)+cl,abs(upr)+cr)
    !  F = 0.5*((FL+FR)+0.4*amax*(PL-PR))
    !  return
    !!endif  
    
    hlplushr=hl+hr
    hs=0.5*hlplushr-0.25*(upr-upl)*hlplushr/(cl+cr) !Positive depth
    if(hl>hdry .and. hr>hdry)then !Both sides wet
      if(hs>hl)then !Shock
        sl=upl-cl*sqrt(0.5*hlplushr*hs/(hl*hl))  
      else          !Rarefaction
        sl=upl-cl
      endif
      if(hs>hr)then !Shock
        sr=upr+cr*sqrt(0.5*hlplushr*hs/(hr*hr))
      else          !Rarefaction
        sr=upr+cr
      endif
    elseif(hl>hdry)then !Left is wet
      sl=upl-cl
      sr=upl+2.0*cl
    elseif(hr>hdry)then !Right is wet
      sl=upr-2.0*cr
      sr=upr+cr
    else !Both sides dry
      F=(/0.0, 0.0, 0.0/)
      return
    endif

    !cs=0.5*(cl+cr)+0.25*(upl-upr)
    !us=0.5*(upl+upr)+(cl-cr)
    !if(hl>=hmin .and. hr>=hmin)then !Both sides wet
    !  sl=min(upl-cl,us-cs)
    !  sr=min(upr+cr,us+cs)
    !elseif(hl>=hmin .and. hr<hmin)then !Left is wet, right is dry
    !  sl=upl-cl
    !  sr=upl+2.0*cl
    !elseif(hl<hmin .and. hr>=hmin)then !Left is dry, right is wet
    !  sl=upr-2.0*cr
    !  sr=upr+cr
    !else !if(hl<hmin .and. hr<hmin)then !both sides dry
    !  F=(/0.0, 0.0, 0.0/)
    !  return
    !endif
    
    if(sl>=0.0)then
      F=(/qpl, qpl*ul, qpl*vl/)
    elseif(sl<0.0 .and. sr>1.0e-10)then
      FL=(/qpl, qpl*ul, qpl*vl/)
      FR=(/qpr, qpr*ur, qpr*vr/)
      PL=(/el, qpl*cn, qpl*sn/)
      PR=(/er, qpr*cn, qpr*sn/)   
      F=(sr*FL-sl*FR+sl*sr*(PR-PL))/(sr-sl)
    else
      F=(/qpr, qpr*ur, qpr*vr/)
    endif    

    return
    end subroutine hll
    
!********************************************************
    subroutine hllc(hl,hr,el,er,ul,ur,vl,vr,cn,sn,F)
! Harten, Lax, and van Leer approximate Riemann solver 
! with restored Contact wave (Toro et al. 1994) 
!
! References:
!  Toro, E.F., Spruce, M., Speares, W. 1994, Restoration 
!  of the contact surface in the HLL-Riemann solver. 
!  Shock Waves. 4, 25–34.
!
! Author: Alex Sanchez, USACE-CHL
!********************************************************
  use flow_def, only: hmin,grav
    implicit none
  !Input
    real(ikind),intent(in) :: hl,hr,el,er,ul,ur,vl,vr,cn,sn
  !Output
  real(ikind),intent(out) :: F(3)
  !Internal Variables    
    real(ikind) :: FL(3),FR(3),PL(3),PR(3),PSR(3),PSL(3)
    real(ikind) :: sl,sr,ss,ups             !,amax,as,lambdal,lambdar,hs
    real(ikind) :: cr,cl,upl,upr,qpl,qpr,cs !ql,qr,hgl,hgr,hlplushr

    !hl=max(hl,hmin);  hr=max(hr,hmin)
    upl=ul*cn+vl*sn;  upr=ur*cn+vr*sn
    cl=sqrt(grav*hl); cr=sqrt(grav*hr)
    qpl=upl*hl;       qpr=upr*hr
    
    !hs =0.5*(hl+hr)-0.25*(upr-upl)*(hl+hr)/(cl+cr)
    !ups=0.5*(ul+ur)-0.25*(hr-hl)*(cl+cr)/(hl+hr)

    !hlplushr=hl+hr
    !hs=0.5*hlplushr-0.25*(upr-upl)*hlplushr/(cl+cr) !Positive depth
    !if(hl>=hmin .and. hr>=hmin)then !Both sides wet
    !  if(hs>hl)then !Shock
    !    sl=upl-cl*sqrt(0.5*hlplushr*hs/(hl*hl))  
    !  else          !Rarefaction
    !    sl=upl-cl
    !  endif
    !  if(hs>hr)then !Shock
    !    sr=upr+cr*sqrt(0.5*hlplushr*hs/(hr*hr))
    !  else          !Rarefaction
    !    sr=upr+cr
    !  endif
    !elseif(hl>=hmin .and. hr<hmin)then !Left is wet, right is dry
    !  sl=upl-cl
    !  sr=upl+2*cl
    !elseif(hl<hmin .and. hr>=hmin)then !Left is dry, right is wet
    !  sl=upr-2*cr
    !  sr=upr+cr  
    !else !if(hl<hmin .and. hr<hmin)then !both sides dry
    !  F=(/0.0, 0.0, 0.0/)
    !  return
    !endif
    
    cs=0.5*(cl+cr)+0.25*(upl-upr)
    ups=0.5*(upl+upr)+(cl-cr)
    
    if(hl>=hmin .and. hr>=hmin)then !Both sides wet
      sl=min(upl-cl,ups-cs)
      sr=min(upr+cr,ups+cs)
    elseif(hl>=hmin .and. hr<hmin)then !Left is wet, right is dry
      sl=upl-cl
      sr=upl+2*cl
    elseif(hl<hmin .and. hr>=hmin)then !Left is dry, right is wet
      sl=upr-2*cr
      sr=upr+cr
    else !if(hl<hmin .and. hr<hmin)then !both sides dry
      F=(/0.0, 0.0, 0.0/)
      return
    endif
    
    ss=ups
    if(sl>=0.0)then
      !hgl=0.5*grav*hl*hl
      !F=(/qpl, qpl*ul+hgl*cn, qpl*vl+hgl*sn/)
      F=(/qpl, qpl*ul, qpl*vl/)
    elseif(sl<0.0 .and. ss>=0.0)then
      !hgl=0.5*grav*hl*hl; hgr=0.5*grav*hr*hr
      !FL=(/qpl, qpl*ul+hgl*cn, qpl*vl+hgl*sn/)
      !FR=(/qpr, qpr*ur+hgr*cn, qpr*vr+hgr*sn/)
      FL=(/qpl, qpl*ul, qpl*vl/)
      !FR=(/qpr, qpr*ur, qpr*vr/)
      PL=(/el, qpl*cn, qpl*sn/)          
      PSL=(/el*(sl-upl)/(sl-ss), hl*(sl-upl)/(sl-ss)*ss*cn, hl*(sl-upl)/(sl-ss)*ss*sn /)
      F=FL+sl*(PSL-PL)
    elseif(sr>=0.0 .and. ss<=0.0)then
      !hgl=0.5*grav*hl*hl; hgr=0.5*grav*hr*hr
      !FL=(/qpl, qpl*ul+hgl*cn, qpl*vl+hgl*sn/)
      !FR=(/qpr, qpr*ur+hgr*cn, qpr*vr+hgr*sn/)
      !FL=(/qpl, qpl*ul, qpl*vl/)
      FR=(/qpr, qpr*ur, qpr*vr/)
      PR=(/er, qpr*cn, qpr*sn/)
      PSR=(/er*(sr-upr)/(sl-ss), hr*(sr-upr)/(sr-ss)*ss*cn, hr*(sl-upr)/(sr-ss)*ss*sn /)            
      F=FR+sr*(PSR-PR)
    else
      !hgr=0.5*grav*hr*hr
      !F=(/qpr, qpr*ur+hgr*cn, qpr*vr+hgr*sn/)
      F=(/qpr, qpr*ur, qpr*vr/)
    endif

    return
    end subroutine hllc
    
!****************************************************
    subroutine roe(hl,hr,el,er,ul,ur,vl,vr,cn,sn,F)
! Roe (1981) approximate Riemann solver
!
! References:
!   Roe, P.L. 1981. Approximate Riemann solvers, 
!   parameter vectors and difference schemes. 
!   Journal of Computational Physics. 43(2), 357–372. 
!
! Author: Alex Sanchez, USACE-CHL
!****************************************************
    use flow_def, only: grav                  !hmin  is never used, commented out   MEB  01/26/2022
    implicit none
    !Input
    real(ikind),intent(in) :: hl,hr,el,er,ul,ur,vl,vr,cn,sn
    !Output
    real(ikind),intent(out) :: F(3)
    !Internal
    real(ikind) :: duml,dumr,cl,cr,hhat,uhat,vhat
    real(ikind) :: chat,uperp,de,du,dv,dupar,duperp,uperpl,uperpr
    real(ikind) :: al1,al3,ar1,ar3,da1,da3,a1,a3
    real(ikind) :: dW(3),R(3,3),A(3,3),FL(3),FR(3),B(3,3)

    !hl=max(hl,hmin);    hr=max(hr,hmin)
    duml=sqrt(hl);      dumr=sqrt(hr)
    cl=sqrt(grav*hl);     cr=sqrt(grav*hr)
    uperpl=ul*cn+vl*sn; uperpr=ur*cn+vr*sn
    hhat=duml*dumr    
    uhat=(duml*ul+dumr*ur)/(duml+dumr)
    vhat=(duml*vl+dumr*vr)/(duml+dumr)
    chat=sqrt(0.5*grav*(hl+hr))
    uperp=uhat*cn+vhat*sn
    de=er-el; du=ur-ul; dv=vr-vl
    dupar=-du*sn+dv*cn
    duperp=du*cn+dv*sn
    dW(1)=0.5*(de-hhat*duperp/chat)
    dW(2)=hhat*dupar
    dW(3)=0.5*(de+hhat*duperp/chat)    
    
    R(1,1)=1.0;          R(1,2)=0.0;    R(1,3)=1.0
    R(2,1)=uhat-chat*cn; R(2,2)=-sn;    R(2,3)=uhat+chat*cn
    R(3,1)=vhat-chat*sn; R(3,2)=cn;     R(3,3)=vhat+chat*sn
    
    A=0.0
    A(1,1)=abs(uperp-chat); A(2,2)=abs(uperp); A(3,3)=abs(uperp+chat)

    !Critical flow fix    
    al1=uperpl-cl; al3=uperpl+cl
    ar1=uperpr-cr; ar3=uperpr+cr
    da1=max(0.0,2*(ar1-al1))
    da3=max(0.0,2*(ar3-al3))
    if(A(1,1)<da1) a1=0.5*(A(1,1)*A(1,1)/da1+da1)
    if(A(3,3)<da3) a3=0.5*(A(3,3)*A(3,3)/da3+da3)

    !Compute interface flux
    FL(1)=uperpl*hl
    FL(2)=ul*uperpl*hl
    FL(3)=vl*uperpl*hl
    FR(1)=uperpr*hr
    FR(2)=ur*uperpr*hr
    FR(3)=vr*uperpr*hr
    B=matmul(R,A)
    F=0.5*(FL + FR - matmul(B,dW))

    return
    end subroutine roe    

end module comp_lib