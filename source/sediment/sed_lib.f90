!===========================================================================
module sed_lib
! CMS Sediment Transport library
!
! Contains the following:
! - Grain size distribution - 
!     diam2diamlim - Calculates the sediment size class limits based
!                    on the size class diameters
!     diamlim2diam - Calculates the sediment size class diameters based
!                    on the size class limits
!     gsd2perdiam  - Calculates the sediment percentile diameter based 
!                    on the grain size distribution
! - Fall velocity - 
!     fallvel_soulsby - Soulsby (1997) fall velocity formula
!     fallvel_vanrijn - Van Rijn (1993) sediment fall velocity formula
!     fallvel_wu_wang  - Wu and Wang (2006) fall velocity formula  
!     fallvel_zhang   - Zhang fall velocity formula
!     fallvel_ruby    - Ruby (1933) fall velocity formula 
! - Porosity - 
!     sed_poro_wu_wang - Empirical formula for bed porosity based only on 
!                the median grain size (Wu and Wang 2006)
! - Incipient motion (Shields paramter and critical shear stress) - 
!     shields_soulsby - Shields parameter for incipient motion 
!                    based on Soulsby (1997) 
!     shields_wuwang - Shields parameter for incipient motion 
!                    based on Wu and Wang (1999)
!     critvel_soulsby - Critical (threshold) depth-averaged velocity 
!                    based on Soulsby (1997)
! - Transport Formula - 
!     sedtrans_cur_soulsby - Total load transport based on the Soulsby-Van Rijn (1997) 
!                transport equation for currents only
!     sedtrans_wavcur_soulsby - Total load transport based on the Soulsby-Van Rijn (1997) 
!                transport equation for currents and waves
!     sedtrans_cur_vanrijn - Bed and suspended load based on the Van Rijn (1984)
!                  transport equations for currents only
!     sedtrans_wavcur_vanrijn - Bed and suspended load based on the Van Rijn (1984)
!                  transport equations for currents and waves
!     sedtrans_wavcur_watanabe - Watanabe total load transport formula for waves and currents
!     shearwatanabecw - Computes the maximum bed shear stress for the 
!                       Watanabe transport formula. 
! - Bed-slope effect on incipient motion -
!     critslpcor_dey - Calculates the slope correction for critical shear 
!               stress using the equation by Dey (2001)
!     effshear_wu - Calculates the effective tractive force for bedload 
!                    by adding the stream-wise gravity to the grain shear stress
! - Hiding and exposure corrections - 
!     hidexp_ashida_michiue - Hiding and exposure coefficient 
!                            based on Ashida and Michiue (1972)
!     hidexp_egiazaroff - Hiding and exposure coefficient  
!                         based on Egiazaroff (1965)
!     hidexp_hayashi - Calculates hiding and exposure coefficient 
!                      based on  Hayashi et al. (1980)
!     hidexp_parker - Hiding and exposure coefficient  
!                       based on Parker and others
!     hidexp_wu - Hiding and exposure coefficient 
!                   based on Wu et al. (2000)  
!     hidexp_day - Hiding and exposure coefficient
!                   based on Day (1980)
!     hidexp_proffitt_sutherland - Hiding and exposure correction coefficient
!                   based on Proffitt and Sutherland (1983)
! - Adapataion parameters - 
!     adaptsusp_armanini_disilvio - Suspended load adaptation coefficient 
!                           based on Armanini and di Silvio (1986)
!     adaptsusp_gallappatti - Suspended load adaptation time 
!                     based on Gallappatti (1983)
!     adaptsusp_lin - Suspended load adaptation coefficient based on Lin (1984)
!
! - Bed load characteristics - 
!     bedvel_vanrijn - Bed load velocity using formula of Van Rijn (1984)
!     bedvel_vanrijn_wu - Bed load velocity using formula of Van Rijn (1984)
!              with recalibrated coefficients of Wu et al. (2006)
!     bedvel_fredsoe - Bed-load velocity based on of Fredsoe and Deigaard (1992)
!
! written by Alex Sanchez, USACE-CHL
!            Weiming Wu, NCCHE
!===========================================================================
    use prec_def
    implicit none
    
contains
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
! Begin Miscelaneous
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

!************************************************   
    function diam2dstar(g,mu,s,d) result(dstar)
! Converts a sediment diameter d to a dimension diameter dstar
! 
! Usage:
!  dstar = diam2dstar(g,mu,s,d);
!
! Input:
!   g - Gravitational constant (default 9.81 m/s^2) [L/T^2]
!   mu - Kinematic viscosity [L^2/T]
!   s - Sediment Specific Gravity [-]
!   d - Sediment diameter [L]
!
! Output:
!   dstar - Dimensionless diameter [-]
!
! Author: Alex Sanchez, USACE-CHL
!************************************************
    implicit none
    real(ikind),intent(in)::g,mu,s,d
    real(ikind):: dstar
    
    !The line below is equivalent to
    !  dstar = d*((s-1.0)*g)**(1.0_ikind/3.0_ikind)/mu**(2.0_ikind/3.0_ikind) 
    !but has less precision error 
    !The precision error is on the order of 1.0e-05 for single precision
    dstar = d*((s-1)*g/mu**2)**(1.0_ikind/3.0_ikind) 
    
    return
    endfunction diam2dstar    
    
!************************************************    
    subroutine diam2diamlim(nd,d,dlim)
! Calculates the sediment size class limits based
! on the size class diameters
!************************************************
    implicit none
    integer,intent(in) :: nd
    real(ikind),intent(in) :: d(nd)
    real(ikind),intent(out) :: dlim(nd+1)
        
    dlim(2:nd)=sqrt(d(2:nd)*d(1:nd-1))
    dlim(1)=d(1)*d(1)/dlim(2)
    dlim(nd+1)=d(nd)*d(nd)/dlim(nd)
    
    return
    endsubroutine diam2diamlim

!************************************************    
    subroutine diamlim2diam(nd,dlim,d)
! Calculates the sediment size class diameters based
! on the size class limits
!************************************************
    implicit none
    integer,intent(in) :: nd
    real(ikind),intent(in) :: dlim(nd+1)
    real(ikind),intent(out) :: d(nd)
    
    d(1:nd)=sqrt(dlim(1:nd)*dlim(2:nd+1))
    
    return
    endsubroutine diamlim2diam
    
!********************************************************  
    subroutine sed_dper(n,d,dlim,logdlim,pb,per,dper)    
! Calculates the sediment percentile (dper) based on the
!
! Input:
!   n - Number of grain classes [-]
!   d(1:n) - Grain class diameters (may contain duplicate values) [L]
!   dlim(1:n+1) - Grain class limits [L]
!   logdlim(1:n+1) - log(dlim) (to avoid recomputing)
!   pb(1:n) - Fraction bed composition [-]
!   per - Percentile (0-100)
!
! Output:
!   dper - Percentile diameter [L]
!
! Author: Alex Sanchez, USACE-CHL
!********************************************************
    implicit none
    integer,intent(in) :: n
    !Input/Output
    real(ikind),intent(in):: d(n),dlim(n+1),logdlim(n+1),pb(n),per
    real(ikind),intent(out):: dper
    !Internal
    integer :: ks,kks,ks0,ks2
    real(ikind):: fac,pbcum(n+1)
    
    !dper=dlim(1)
    !pbcum(1)=0.0
    !do ks=1,n
    !  pcum(ks+1)=pcum(ks)+pb(ks)
    !  if(pcum(ks+1)>=per .and. pcum(ks)<per)then
    !    fac=(pcum(ks+1)-per)/(pcum(ks+1)-pcum(ks))
    !    dper=exp((1.0-fac)*log(dlim(ks+1))+fac*log(dlim(ks)))
    !    return
    !  endif      
    !enddo 
    !if(pcum(n+1)>per) dper=d(n)
    
    dper = dlim(1) 
    pbcum(1) = 0.0
    do ks=1,n
      pbcum(ks+1) = pbcum(ks) + pb(ks)
    enddo
    if(pbcum(n+1)<=per)then
      dper = dlim(n+1)
      return
    endif
    do ks=1,n
      if(pbcum(ks+1)>=per .and. pbcum(ks)<per)then
        !Transverse duplicate diameters
        do kks=ks,1,-1
          if(abs(d(kks)-d(ks))<1.0e-6)then
            ks0 = kks
          endif
        enddo
        do kks=ks,n
          if(abs(d(kks)-d(ks))<1.0e-6)then
            ks2 = kks+1
          endif
        enddo
        exit
      endif
    enddo !ks
    fac = (pbcum(ks2)-per)/(pbcum(ks2)-pbcum(ks0))
    dper = exp((1.0-fac)*logdlim(ks2)+fac*logdlim(ks0)) !Interpolation in log-space

    return
    endsubroutine sed_dper    

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! End Miscelaneous
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
! Begin fall velocity formulas
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
!*************************************************************
    function fallvel_soulsby(mu,d,ds,cv) result(ws)
! Soulsby (1997) sediment fall velocity formula
!
! Input:
!   mu - Kinematic viscosity [L^2/T]
!   d - characteristic grain size diameter [L]
!   ds - nondimensional grain size [-]
!   cv - volume concentration [-]
!
! Output:
!   ws - fall velocity [m/s]
!
! Reference: 
!    Soulsby, R.L. (1997). "Dynamics of marine sands, 
!      a manual for practical applications".
!      H.R. Wallingford, UK: Thomas Telford.
!
! Author: Alex Sanchez, USACE-CHL
!*************************************************************
    implicit none
    real(ikind),intent(in) :: mu,ds,d,cv
    real(ikind) :: ws 

    ws=mu/d*(sqrt(107.3296+1.049*(1.0-cv)**4.7*ds**3)-10.36) !107.3296=10.36^2
    
    return
    endfunction fallvel_soulsby

!**************************************************
    function fallvel_vanrijn(g,mu,s,d) result(ws)
! Van Rijn (1993) sediment fall velocity formula
!
! Usage:
!   ws = fallvel_vanrijn(s1grav,mu,d)
!
! Input:
!   g - Gravity [L/T^2]
!   mu - Kinematic viscosity of water [L^2/T]
!   d - Grain size diameter [L]
!   s - Grain specific gravity  [-]
!
! Output:
!   ws - fall velocity [L/T]
!
! Reference: 
!   Van Rijn, L.C. 1993. Principles of Sediment 
!     Transport in Rivers, Estuaries and Coastal Seas. 
!     Aqua Publications, The Netherlands. 
!
! Author: Alex Sanchez, USACE-CHL
!**************************************************s
    implicit none
    real(ikind),intent(in) :: g,mu,s,d
    real(ikind) :: ws 

    if(6.5e-5<d .and. d<=1.0e-4)then
      ws=(s-1)*g*d**2/(18.0*mu)
    elseif(1.0e-4<d .and. d<=1.0e-3)then
      ws=10.0*mu/d*(sqrt(1.0+0.01*(s-1)*g*d**3/mu**2)-1.0)
    else !1.0e-3<d
      ws=1.1*sqrt((s-1)*g*d)
    endif
    
    return
    endfunction fallvel_vanrijn

!*********************************************************
    function fallvel_wu_wang(mu,sp,d,ds) result(ws)
! Wu and Wang (2006)sediment fall velocity formula    
!
! Input:
!   mu - Kinematic viscosity [L/T^2]
!   sp - Corey shape factor [-]
!   d - Characteristic grain size diameter [L]
!   ds - nondimensional grain size [-]
!
! Output:
!   ws - fall velocity [L/T]
!
! Reference: 
!   Wu, W., and Wang. S.S.Y. (2006) Formulas for sediment 
!     porosity and settling velocity, Journal of Hydraulic
!     Engineering, ASCE, 132(8), 858-862.
!
! Author: Alex Sanchez, USACE-CHL
!**********************************************************
    implicit none
    real(ikind),intent(in) :: mu,sp,d,ds
    real(ikind) :: cM,cN,en,a,b,ws
    
    cM=53.5*exp(-0.65*sp)
    cN=5.65*exp(-2.5*sp)    
    en=0.7+0.9*sp    
    a=cM*mu/cN
    b=1.333333*cN/cM**2    
    ws=a/d*(sqrt(0.25+(b*ds**3)**(1.0/en))-0.5)**en
    
    return
    endfunction fallvel_wu_wang
    
!*****************************************************
    function fallvel_zhang(g,mu,s,d) result(ws)
! Zhang (1961) Sediment fall velocity
!
! Input:
!   g - Gravity [L/T^2]
!   mu - Kinematic viscosity of water [L^2/T]
!   d - Grain size diameter [L]
!   s - Grain specific gravity  [-]
!
! Output:
!   ws - fall velocity [L/T]
!
! References:
!   Zhang, R.J. (1961) River Dynamics, Industry Press, 
!     Beijing, China (in Chinese).
!******************************************************
    implicit none
    real(ikind),intent(in) :: g,mu,s,d
    real(ikind) :: ws
    
    ws = sqrt((13.95*mu/d)**2+1.09*(s-1)*g*d)-13.95*mu/d 
    
    return
    endfunction fallvel_zhang
    
!*************************************************************
    function fallvel_ruby(g,mu,s,d) result(ws)
!FALLVEL_RUBY Compute sediment fall velocity based on Ruby (1933)
!
! Usage:
!   ws = fallvel_ruby(g,mu,d,s)
!
! Input:
!   g - Gravity [L/T^2]
!   mu - Kinematic viscosity of water [L^2/T]
!   d - Grain size diameter [L]
!   s - Grain specific gravity  [-]
!
! Output:
!   ws - Sediment fall velocity [L/T]
!
! References:
!   Ruby, W.W. 1933. Settling velocities of gravel, sand and silt 
!     particles. American Journal of Science, 25(148), 325-338. 
!
! Author: Alex Sanchez, USACE-CHL
!*************************************************************
    implicit none
    real(ikind),intent(in) :: g,mu,s,d
    real(ikind) :: a,f,ws
    
    a = (36.0*mu**2)/(g*(s-1.0)*d**3) ![-]
    f = sqrt(0.66667 + a) - sqrt(a) ![-]
    ws = f*sqrt((s-1.0)*g*d)

    return
    endfunction fallvel_ruby
    
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! End Sediment Fall Velocity Formulas
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
! Begin Bed Porosity Formulas
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

!*************************************************************
    function sed_poro_wu_wang(diam50) result(p)
! Calculates the sediment porosity as a function of d50
! using the formula of Wu and Wang (2006)
!
! Input:
!   diam50 - Median grain size, mm
! Output: 
!   p - Sediment porosity
!
! Reference: 
!   Wu, W., and Wang. S.S.Y. (2006) Formulas for sediment 
!     porosity and settling velocity, Journal of Hydraulic
!     Engineering, ASCE, 132(8), 858-862.
!
! written by Alex Sanchez, USACE-CHL
!***************************************************************
    real(ikind),intent(in) :: diam50
    real(ikind) :: p
    
    p = 0.13 + 0.21/((diam50 + 0.002)**0.21)
    
    return
    endfunction sed_poro_wu_wang
    
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>    
! End Bed Porosity Formulas
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
! Begin Incipient Motion Formulas
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

!******************************************************************
    function critslpcor_dey(sl,st,sr) result(val)
! Calculates the slope correction for critical shear stress
! using the equation by Dey (2001) and recalibrated coefficients
! from Walstra et al. (2007)
!
! Input:
!  sr - Repose (critical) slope (always positive) (=tan(phir)) [-]
!  sl - Longitudinal slope (=tan(beta)) [-]
!  st - Transverse slope (=tan(gamma)) [-]
!
! Output:
!  val - Correction factor defined by tcr/tcr0
!        where tcr is corrected critical shear stress and
!        tcr0 is the critical shear stress for a flat bed
!
! References:
!   Dey, S. (2001). "Experimental study on incipient motion of 
!     sediment particles on generalized sloping fluvial beds". 
!     International Journal of Sediment Research, 16(3), 391-398.
!  Walstra, D.J.R., Van Rijn, L.C., Ormondt, M.V., Brière, C., 
!    and Talmon, A.M. (2007). "The effects of bed slope and 
!    wave skewness on sediment transport and morphology", 
!    Proceedings of Coastal Sediments '07, 14 pp.
!
! written by Alex Sanchez, USACE-CHL
!******************************************************************
    real(ikind),intent(in) :: sl,st,sr
    real(ikind) :: val
    
    val = (1.0+sl/sr)**0.75*(1.0-abs(st)/sr)**0.37 !Note: repose angle always positive
    
    return
    endfunction critslpcor_dey

!*********************************************************************
    function effshear_wu(slp,angrep,taubp,tauc) result(taube)
! Calculates the effective tractive force for bedload 
! by adding the stream-wise gravity to the grain shear stress
!
! Input:
!   slp = Stream-wise slope [rad]
!   angrep = Angle of repose [rad]
!   taubp = Grain shear stress [N/m^2]
!   tauc = Critical shear stress for incipient motion [N/m^2]
!
! Output:
!   taube = effective tractive force for bedload transport
!
! References:
!   Wu. W. 2004. Depth-averaged two-dimensional numerical modeling
!     of unsteady flow and nonuniform sediment transport in open channels. 
!     Journal of Hydraulic Engineering, 130(10), 1013-1024.
! 
! written by Alex Sanchez, USACE-CHLs
!*********************************************************************
    implicit none
    real(ikind),intent(in) :: slp,taubp,tauc,angrep
    real(ikind) :: val,lambda0,taube

    val = sin(slp)/sin(angrep)
    if(slp>0.0)then
      lambda0 = 1.0 + 0.22*(taubp/tauc)**0.15*exp(2.0*val) !Notes: slp<=angrep and tauc>small 
      taube = taubp + lambda0*tauc*val
    else
      taube = taubp + tauc*val
    endif
    
    return
    endfunction effshear_wu

!*****************************************************************    
    function shields_soulsby(dstar) result(thcr)
! Calculates the shields parameter based on Soulsby (1997)
!
! Input:
!   dstar - nondimensional grain size 
!      dstar = (g(s-1)/mu^2)**(1/3)*d
!      where
!      g - gravity [L/T^2]
!      s - sediment specific gravity [-]
!      d - grain size diameter [L]
!      mu - kinematic viscosity of water [L^2/T]
!
! Output:
!   thcr - Shields parameter
!
! References:
!   Soulsby, R.L. (1997). "Dynamics of marine sands", 
!      Thomas Telford, 249 p.
!
! written by Alex Sanchez, USACE-CHL
! last modified 12-03-12
!*****************************************************************
    implicit none
    real(ikind),intent(in) :: dstar
    real(ikind) :: thcr

    thcr = 0.3/(1.0+1.2*dstar)+0.055*(1.0-exp(-0.02*dstar))
    
    return
    endfunction shields_soulsby
    
!*******************************************************************************    
    function critvel_soulsby(g,h,s,d,dstar) result(Ucr)
! Critical (threshold) depth-averaged velocity based on Soulsby (1997)
! 
! Input:
!   g - Gravity [L/T^2]
!   h - Water depth [L]
!   d - Grain size [L]
!   dstar - nondimensional grain size [-]
!
! Output:
!   Ucr - Depth-averaged critical (threshold) velocity [L/T]
!
! References:
!   Soulsby, R.L. 1997. Dynamics of marine sands. Thomas Telford, 249 p.
!
! Author: Alex Sanchez, USACE-CHL
!*******************************************************************************
    implicit none
    real(ikind),intent(in) :: g,h,s,d,dstar
    real(ikind) :: Ucr,thcr
    
    thcr = shields_soulsby(dstar)
    Ucr = 7.0*((h/d)**0.1428571429)*sqrt((s-1)*g*d*thcr) !0.1428571429=1/7
    
    return
    endfunction critvel_soulsby
    
!*******************************************************************************    
    function critvel_vanrijn(h,d50,d90) result(Ucrc)
! Critical (threshold) depth-averaged velocity base on Van Rijn (1984)
! 
! Input:
!   d50 - Median grain size [m]
!   d90 - 90th percential grain size [m]
!   h - total water depth [m]
!
! Output:
!   Ucrc - Depth-averaged critical (threshold) velocity [m/s]
!
! References:
!   van Rijn, L.C. 1984. Sediment transport, Part III: 
!     Bed forms and alluvial roughness. Journal of Hydraulic Engineering, 
!     ASCE, 110(12), 1733-1754. 
!
! written by Alex Sanchez, USACE-CHL
!*******************************************************************************
    implicit none
    real(ikind),intent(in) :: h,d50,d90
    real(ikind) :: Ucrc
    
    !Critical shear stresses are calculated using Komar and Miller as in Van Rijn 2007
    if(d50<0.0005)then
      Ucrc = 0.19*(d50**0.1)*log10(4.0*h/d90) !Depth-averaged velocity [m/s]
    else
      Ucrc = 8.5*(d50**0.6)*log10(4.0*h/d90)  !Depth-averaged velocity [m/s]
    endif 
    
    return
    endfunction critvel_vanrijn
    
!*******************************************************************************    
    function critvel_komar_miller(g,Tp,s,d50) result(Ucrw)
! Critical (threshold) depth-averaged velocity base on Komar and Miller 1974
! 
! Input:
!   d50 - Median grain size [m]
!   Tp - Peak wave period [s]
!
! Output:
!   Ucrw - Wave critical (threshold) velocity [m/s]
!
! References:
!   van Rijn, L.C., Nieuwjaar, M., Van der Kaaij, T., Nap, E., and Van
!     Kampen, A. (1993). "Transport of fine sands by currents and waves".
!     Journal of Waterway, Port, Coastal, Ocean Engineering". 1192, 123–143.
!   Komar, P.D., and Miller, M.C. (1975). “On the comparison between the
!     threshold of sediment motion under waves and unidirectional currents
!     with a discussion of the practical evaluation of the threshold”. 
!     Journal of Sedimentation and Petrology, 45, 362–367.
!
! written by Alex Sanchez, USACE-CHL
!*******************************************************************************
    implicit none
    real(ikind),intent(in) :: g,Tp,s,d50
    real(ikind) :: Ucrw
    
    if(d50<0.0005)then
      Ucrw = 0.24*(((s-1.0)*g)**0.66)*(d50**0.33)*Tp**0.33
    else
      Ucrw = 0.95*(((s-1.0)*g)**0.57)*(d50**0.43)*Tp**0.14
    endif
    
    return
    endfunction critvel_komar_miller
    
!*****************************************************************    
    function shields_wuwang(dstar) result(thcr)
! Shields parameter for incipient motion 
! based on Wu and Wang (1999)
!
! Reference:
!   Wu.W., Wang, S.S.Y (1999) "Movable bed roughness 
!     in alluvial rivers", Journal of Hydraulic Engineering,
!     125(12), 1309-1312.
!
! written by Alex Sanchez, USACE-CHL
!*****************************************************************
    implicit none
    real(ikind),intent(in) :: dstar
    real(ikind) :: thcr

    if(dstar<1.5)then
      thcr = 0.126*dstar**(-0.44)
    elseif(1.5<=dstar .and. dstar<10.0)then
      thcr = 0.131*dstar**(-0.55)  
    elseif(10.0<=dstar .and. dstar<20.0)then
      thcr = 0.0685*dstar**(-0.27)  
    elseif(20.0<=dstar .and. dstar<40.0)then
      thcr = 0.0173*dstar**0.19  
    elseif(40.0<=dstar .and. dstar<150.0)then
      thcr = 0.0115*dstar**0.30  
    else !150<=dstar
      thcr = 0.052  
    endif
    
    return
    endfunction shields_wuwang

!!************************************************************    
!    subroutine sedscalefac(sx,sy,ux,uy,um,fcr,fb,fs)
!!************************************************************
!    use sed_def, only: ibedslop,a_repose,scalesus,scalebed
!    use prec_def
!    implicit none
!    !Input/Output
!    real(ikind),intent(in) :: sx,sy,ux,uy,um
!    real(ikind),intent(inout) :: fcr(nsed),fb,fs(nsed)
!    !Internal variables
!    real(ikind) :: sl,st,slp
!    
!    !Streamwise slope, positive uphill
!    if(ibedslope==1)then !Dey (2001)
!      sl = (sx*ux+sy*uy)/um !Longitudinal slope
!      st = (sy*ux-sx*uy)/um !Transverse slope
!      slpcor = critslpcor_dey(sl,st,a_repose)
!      !fb=scalebed and fs=scalesus remain unchanged
!      do ks=1,nsed
!        varsigma(i,ks)=sqrt(slpcor*varsigma(i,ks)) !Both corrections are included in critial shear
!      enddo
!    elseif(ibedslope==2)then !Bailard (1981)
!      !Note that it does not include a correction for the transverse slope
!      slp = (sx*ux+sy*uy)/um !streamwise slope, positive upslope
!      fb = scalebed*(1.0-betaslope*slp) !Bailard (1981) approach
!      fs(:) = scalesus*(1.0-effslope*uv(i)/wsfall(:)*slp)
!      !slpcor=1.0 remains the same
!      do ks=1,nsed
!        fcr(ks)=sqrt(varsigma(i,ks)) !Both corrections are included in critial shear
!      enddo
!    !else !Note: method of Wu (2000) not valid for Soulsby-van Rijn transport formula   
!    endif
!    do ks=1,nsed
!      fcr(ks)=sqrt(varsigma(i,ks)) !Both corrections are included in critial shear
!    enddo
!    
!    return
!    endsubroutine sedscalefac
    
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! End Incipient Motion Formulas
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
! Begin Sediment Transport Formula
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

!*********************************************************************************
    subroutine sedtrans_wavcur_soulsby(g,h,Uc,Cd,Urms,s,d50,d90,dstar,fcr,Qb,Qs)
! Calculates the total load transport based on the Soulsby-Van Rijn (1997) 
! transport equation for both waves and currents
!
! Input:
!  g - Gravity [L/T^2]
!  h - total water depth [L]
!  Uc - Current velocity [L/T]
!  Cd - Drag coefficient [-]
!  Urms - root-mean-squared wave bottom orbital velocity amplitude
!         defined as Urms^2 = var(ubw) = int(Su*df).
!  s - Specific gravity [-]
!  d50 - Median grain size [L]
!  d90 - 90th percential grain size [L]
!  dstar - nondimensional grain size [-]
!  fcr - critical velocity correction factor [-] 
!        (may include hiding and exposure an slope corrections)
!
! Output:
!  Qb - Bed-load volume transport rate [L^2/T]
!  Qs - Suspended-load volume transport rate [L^2/T]
!
! Reference:
!  Soulsby, R.L. 1997. Dynamics of marine sands, Thomas Telford, London.
!
! Notes: 
!  (1) The root-mean-squared orbital velocity Urms is defined by
!      Soulsby (1987,1997); and CMS as
!            Urms^2 = int(Su*df)   
!      where Su is the orbital velocity spectral density
!      and f is the wave frequency.
!      The definition by others such as Madsen (1994) is
!            Urms^2 = 2*int(Su*df)
!      This leades to a difference of a factor of sqrt(2)
!  (2) The code has been tested with Soulsby's example on pages 191-192 
!      of "Dynamics of Marine Sands"
!
! written by Alex Sanchez, USACE-CHL
!*********************************************************************************
    implicit none    
    real(ikind),intent(in) :: g,h,Uc,Cd,Urms,s,d50,d90,dstar,fcr
    real(ikind),intent(out) :: Qb,Qs
    real(ikind) :: Ucr,Me,val,Ue
    
    Ucr = sqrt(fcr)*critvel_vanrijn(h,d50,d90) !Note correction
    !Note: using Urms is more convenient then Uw in van Rijn since it works for both regular and random waves
    Ue = sqrt(Uc**2 + 0.018/Cd*Urms**2) 
    if(Ue<=Ucr)then
      Qb=0.0    
      Qs=0.0
      return !Check for early exit
    endif  
    !Note: Me is exactly the same for van Rijn but with different Ue
    Me = (Ue - Ucr)/sqrt((s-1.0)*g*d50) !Note 2.4/2 = 1.2 exponent above
    val = Uc*Me**2.4
    Qb = 0.005*val*h*(d50/h)**1.2    !Note: Me exponent different from van Rijn
    Qs = 0.012*val*d50*dstar**(-0.6) !Note: eq. the same as for Van Rijn
    
    return
    endsubroutine sedtrans_wavcur_soulsby
    
!***************************************************************************    
    subroutine sedtrans_cur_soulsby(g,h,Uc,s,d50,d90,dstar,fcr,Qb,Qs)
! Calculates the total load transport based on the Soulsby-Van Rijn (1997) 
! transport equation for currents only
!
! Input:
!  g - Gravity [L/T^2]
!  h - total water depth [m]
!  Uc - current velocity [m/s]
!  s - Specific gravity [-]
!  d50 - Median grain size [m]
!  d90 - 90th percential grain size [m]
!  dstar - nondimensional grain size [-]
!  fcr - critical velocity correction factor [-] 
!        (may include hiding and exposure an slope corrections)
!
! Output:
!  Qb = bedload volume transport rate [m^2/s]
!  Qs = suspended load volume transport rate [m^2/s]
!
! Reference:
!  Soulsby, R.L. 1997. Dynamics of marine sands, Thomas Telford, London.
!
! Notes: 
!  (1) Code has been tested with Soulsby's example on pages 191-192 
!      of "Dynamics of Marine Sands"
!
! written by Alex Sanchez, USACE-CHL
!******************************************************************
    implicit none    
    real(ikind),intent(in) :: g,h,Uc,s,d50,d90,dstar,fcr
    real(ikind),intent(out) :: Qb,Qs
    real(ikind) :: Me,val,Ucr
    
    Ucr = sqrt(fcr)*critvel_vanrijn(h,d50,d90) !Note correction
    if(Uc<=Ucr)then
      Qb=0.0
      Qs=0.0
      return !Check for early exit
    endif  
    Me = (Uc - Ucr)/sqrt((s-1.0)*g*d50)
    val = Uc*Me**2.4
    Qb = 0.005*val*h*(d50/h)**1.2
    Qs = 0.012*val*d50*dstar**(-0.6)
    
    return
    endsubroutine sedtrans_cur_soulsby

!*********************************************************************************    
    subroutine sedtrans_wavcur_vanrijn(g,h,Uc,Uw,Tp,s,d50,d90,dstar,fcr,Qb,Qs)
! Calculates the bed and suspended load based on the Van Rijn 1984 equations
! for both waves and currents
!
! Input:
!  g - Gravity [L/T^2]
!  h - total water depth [L]
!  Uc - current velocity [L/T]
!  Uw - bottom wave orbital velocity based on the significant wave height [L/T]
!  Tp - peak wave period [T]
!  d50 - Median grain size [L]
!  d90 - 90th percential grain size [L]
!  dstar - nondimensional grain size [-]
!  fcr - critical velocity correction factor [-] 
!        (may include hiding and exposure an slope corrections)
!
! Output:
!  Qb - Bed-load volume transport rate [L^2/T]
!  Qs - Suspended-load volume transport rate [L^2/T]
!
! References:
!   van Rijn, L.C. (1984). "Sediment transport, Part I: 
!     Bed load transport". Journal of Hydraulic Engineering, 
!     ASCE, 110(10), 1431-1456. 
!   van Rijn, L.C. (1984b). "Sediment transport. Part II: 
!     Suspended load trans-port". Journal of Hydraulic 
!     Engineering, ASCE, 110(11), 1613–1641.
!   van Rijn, L.C. (2007a). "Unified view of sediment 
!     transport by currents and waves. Part I: Initiation of motion, 
!     bed roughness, and bed-load transport". 
!     Journal of Hydraulic Engineering, 133(6), 649-667. 
!   van Rijn, L.C. (2007b). "Unified view of sediment 
!     transport by currents and waves. Part II: Suspended transport". 
!     Journal of Hydraulic Engineering, 133(6), 668-689.
!
! written by Alex Sanchez, USACE-CHL
!**********************************************************************************
    implicit none
    !Input/Output
    real(ikind),intent(in) :: g,h,Uc,Uw,Tp,s,d50,d90,dstar,fcr
    real(ikind),intent(out) :: Qb,Qs
    !Internal Variables
    real(ikind) :: Ucr,Ucrw,Ucrc,beta,Ue,Me
      
    Ucrc = critvel_vanrijn(h,d50,d90)
    Ucrw = critvel_komar_miller(g,Tp,s,d50)
    beta = Uc/max(Uc+Uw,1.0e-6)
    Ucr = sqrt(fcr)*(beta*Ucrc + (1.0-beta)*Ucrw) !Apply correction here
    Ue = Uc + 0.4*Uw !effective velocity 0.4 for irregular waves and 0.8 for regular waves
    if(Ue<=Ucr+1.0e-6)then
      Qb = 0.0
      Qs = 0.0
      return !Check for early exit
    endif  
    !Note: Me is exactly the same for Soulsby-van Rijn but with different Ue
    Me = (Ue - Ucr)/sqrt((s-1.0)*g*d50)
    Qb = 0.015*Uc*Me**1.5*h*(d50/h)**1.2    !L^2/T, Van Rijn (2007) I, Eq. 12  (Note Me exponent different from Soulsby-Van Rijn)
    Qs = 0.012*Uc*Me**2.4*d50*dstar**(-0.6) !L^2/T, Van Rijn (2007) II, Eq. 17 (Note: eq. the same as Soulsby-van Rijn)
    
    return
    endsubroutine sedtrans_wavcur_vanrijn

!****************************************************************************    
    subroutine sedtrans_cur_vanrijn(g,h,Uc,s,d50,d90,dstar,fcr,Qb,Qs)
! Calculates the bed and suspended load based on the Van Rijn 1984 equations
! for currents only
!
! Input:
!  g - Gravity [L/T^2]
!  h - total water depth [L]
!  Uc - current velocity [L/T]
!  s - Specific gravity [-]
!  d50 - Median grain size [L]
!  d90 - 90th percential grain size [L]
!  dstar - nondimensional grain size [-]
!  fcr - critical velocity correction factor [-] 
!        (may include hiding and exposure an slope corrections)
!
! Output:
!  Qb = Bed-load volume transport rate [L^2/T]
!  Qs = Suspended-load volume transport rate [L^2/T]
!
! References:
!   van Rijn, L.C. (1984). "Sediment transport, Part I: 
!     Bed load transport". Journal of Hydraulic Engineering, 
!     ASCE, 110(10), 1431-1456. 
!   van Rijn, L.C. (1984b). "Sediment transport. Part II: 
!     Suspended load trans-port". Journal of Hydraulic 
!     Engineering, ASCE, 110(11), 1613–1641.
!   van Rijn, L.C. (2007a). "Unified view of sediment 
!     transport by currents and waves. Part I: Initiation of motion, 
!     bed roughness, and bed-load transport". 
!     Journal of Hydraulic Engineering, 133(6), 649-667. 
!   van Rijn, L.C. (2007b). "Unified view of sediment 
!     transport by currents and waves. Part II: Suspended transport". 
!     Journal of Hydraulic Engineering, 133(6), 668-689.
!
! Author: Alex Sanchez, USACE-CHL
!**************************************************************************** 
    implicit none
    !Input/Output
    real(ikind),intent(in) :: g,h,Uc,s,d50,d90,dstar,fcr
    real(ikind),intent(out) :: Qb,Qs
    !Internal Variables
    real(ikind) :: Ucrc,Me
    
    Ucrc = sqrt(fcr)*critvel_vanrijn(h,d50,d90) !Apply correction here
    if(Uc<=Ucrc)then
      Qb=0.0; Qs=0.0; 
      return !Check for early exit
    endif
    Me = (Uc - Ucrc)/sqrt((s-1.0)*g*d50) !Note: Me exactly the same as in Soulsby-van Rijn
    Qb = 0.015*Uc*Me**1.5*h*((d50/h)**1.2)    !L^2/T, Van Rijn (2007) I, Eq. 12
    Qs = 0.012*Uc*Me**2.4*d50*(dstar**(-0.6)) !L^2/T, Van Rijn (2007) II, Eq. 17, Eq. same as in Soulsby-van Rijn
    
    return
    endsubroutine sedtrans_cur_vanrijn
    
!******************************************************************  
    subroutine shearwatanabecw(rhow,d50,taub,phi,uw,Tp,taumax)
! Computes the maximum bed shear stress for the Watanabe
! transport formula. 
! written by Alex Sanchez, USACE-CHL
!*******************************************************************
    use const_def, only: twopi
    use fric_lib, only: fric_wavcurmax_soulsby
    use prec_def
    implicit none    
    real(ikind),intent(in) :: rhow,d50,taub,phi,uw,Tp
    real(ikind),intent(out) :: taumax        
    real(ikind) :: Aw,r,fw,tauw
    
    Aw = uw*Tp/twopi
    r = max(Aw/(2.5*d50),1.0)    !Relative roughness
    fw = exp(5.5*r**(-0.2)-6.3)  !Wave friction factor
    tauw = 0.5*rhow*fw*uw*uw   !Wave bed shear stress
    taumax = fric_wavcurmax_soulsby(taub,tauw,phi)
    
    return
    endsubroutine shearwatanabecw
        
!******************************************************************  
    subroutine sedtrans_wavcur_watanabe(Uc,taumax,taucr,fcr,Qtot)
! Computes the Watanabe total load for waves and currents
! written by Alex Sanchez, USACE-CHL
!*******************************************************************    
    use sed_def, only: Awidg
    implicit none    
    real(ikind),intent(in) :: Uc,taumax,taucr,fcr
    real(ikind),intent(out) :: Qtot        
    
    Qtot = Awidg*max(0.0,taumax-fcr*taucr)*Uc

    return
    endsubroutine sedtrans_wavcur_watanabe
    
!******************************************************************  
    function sedtrans_ackers_white(g,rho,h,U,ustar,s,dstar,d,p,etak) result(qtk)
! SEDTRANS_ACKERS_WHITE_DAY Total sediment transport formula based on 
! Akers and White (1973) with the hiding and exposure correction by 
! Day (1980) using the revised version of the coefficients from 
! HR Wallingford (1990)
! 
! Input:
!   g - Gravity (~9.81) [L/T^2]
!   rho - Water density [M/L^3]
!   h - Total water depth [L]
!   U - Depth-averaged current velocity [L/T]
!   ustar - Shear velocity [L/T]
!   s - Specific gravity [-]
!   dstar - Nondimensional grain size [-]
!   d - Characterestic grain size [L]
!   p - Fraction of each grain size [-]
!   etal - Hiding and exposure correction factor [-]
!
! Output:
!   qtk - Fractional total sediment transport [M/L/T]
!  
! References:
!   Ackers, P., and White, W.R. (1973). Sediment transport: A new approach
!     and analysis, Proceedings. ASCE, Journal of the Hydraulic Division, 
!     99, HY11, pp. 2041-2060. 
!   HR Wallingford (1990). Sediment Transport, the Ackers and White 
!     Theory Revised, Report SR237, England.
!
! Author: Alex Sanchez, USACE-CHL
!******************************************************************
    implicit none
    real(ikind),intent(in):: g,rho,h,U,ustar,s,dstar,d,p,etak
    real(ikind) :: en,em,A,C,Fgr,Ggr,X,qtk,dstarlim
    real(ikind),parameter :: alpha = 10.0 !Assumed value in HEC7 and HEC-RAS
    
    if(dstar<60.0)then
      dstarlim = max(dstar,1.0) !Limit *** IMPORTANT ****
      en = 1.0 - 0.56*log10(dstarlim)  !Sediment size-related transition exponent
      A = 0.23*dstar**(-0.5) + 0.14 !Initial motion parameter
      em = 6.83/dstar + 1.67 !Sediment transport function exponent (from HR Wallingford 1990)
      C = 10**(2.79*log10(dstar) - 0.98*log10(dstar)**2 - 3.46) !Sediment transpor function coefficient (from HR Wallingford 1990)
    else
      en = 0.0  !Sediment size-related transition exponent   
      A = 0.17  !Initial motion parameter
      em = 1.78 !Sediment transport function exponent (from HR Wallingford 1990)
      C = 0.025 !Sediment transpor function coefficient (from HR Wallingford 1990)
    endif

    !Sediment mobility number    
    Fgr = etak*ustar**en/sqrt((s-1)*g*d)*(U/(sqrt(32.0)*log10(alpha*h/d)))**(1.0-en)

    !Transport parameter
    Ggr = C*(Fgr/A - 1.0)**em
    Ggr = max(Ggr,0.0)

    !Sediment flux in parts per million by fluid weight
    X = p*Ggr*s*d/(h*(ustar/U)**en)

    !Total sediment transport rate [M/m/s]
    qtk = rho*g*h*U*X

    return
    endfunction sedtrans_ackers_white
    
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! End Sediment Transport Formula
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
! Begin Hiding and Exposure Formula
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

!***********************************************************
    function hidexp_egiazaroff(dk,dm) result(psi)
! Calculates hiding and exposure coefficient based on 
! Egiazaroff (1965)
!
! Input:
!   dk - Characteristic grain size
!   dm - Median grain size for the surface/mixing/active layer
!
! Output:
!   psi - Hiding and exposure correction function
!
! Reference:
!   Egiazaroff, I.V. 1965. Calculation of nonuniform sediment 
!     concentration, Journal of Hydraulics Division, ASCE, 
!     91(H14), 225-247. 
!
! written by Alex Sanchez, USACE-CHL
!***********************************************************
    use prec_def
    implicit none
    real(ikind),intent(in) :: dk,dm
    real(ikind) :: psi
    
    psi = (1.278753/log10(19.0*dk/dm))**2 !log10(19)=1.2788
    
    endfunction hidexp_egiazaroff

!***********************************************************
    function hidexp_ashida_michiue(dk,dm) result(psi)
! Calculates hiding and exposure coefficient based on 
! Ashida and Michiue (1971)
!
! Input:
!   dk - Characteristic grain size
!   dm - Median grain size for the surface/mixing/active layer
!
! Output:
!   psi - Hiding and exposure correction function
!
! Reference:
!   Ashida, K. and Michiue, M. 1971. An investigation of  
!     river bed degradation downstream of a dam, 
!     Proc. 14th Congress of the IAHR,  245-255.
!
! written by Alex Sanchez, USACE-CHL
!***********************************************************
    implicit none
    real(ikind),intent(in) :: dk,dm
    real(ikind) :: psi
    
    psi = dk/dm
    if(psi>=0.4) psi = (0.90309/log10(19.0*psi))**2 !log10(8)=0.90309
    
    endfunction hidexp_ashida_michiue
    
!***********************************************************
    function hidexp_hayashi(dk,dm) result(psi)
! Calculates hiding and exposure coefficient based on 
! Hayashi et al. (1980)
!
! Input:
!   dk - Characteristic grain size
!   dm - Median grain size for the surface/mixing/active layer
!
! Output:
!   psi - Hiding and exposure correction function
!
! Reference:
!   Hayashi, T.S., Ozaki, and Ichibashi, T. 1981. Study on 
!     bed load transport of sediment mixture, Proceedings 
!     of the 24th Japanese Conference on Hydraulics, Japan.
!
! Author: Alex Sanchez, USACE-CHL
!***********************************************************
    implicit none
    real(ikind),intent(in) :: dk,dm
    real(ikind) :: psi
    
    psi = dk/dm
    if(psi>=1.0) psi = (0.90309/log10(8.0*psi))**2 !log10(8)=0.90309
    
    endfunction hidexp_hayashi

!***********************************************************
    function hidexp_parker(dk,d50m,cm) result(psi)
! Calculates hiding and exposure coefficient based on 
! Parker and others
!
! Input:
!   dk   - Characteristic dimater
!   cm   - Hiding and exposure coefficient
!   d50m  - Median grain size of mixing layer    
!
! Output:
!   psi - Hiding and exposure correction
!
! written by Alex Sanchez, USACE-CHL    
!***********************************************************
    implicit none
    real(ikind),intent(in) :: dk,cm,d50m
    real(ikind) :: psi             !Hiding and exposure correction    
    
    psi = (d50m/dk)**cm !Note mhe is positive
    
    endfunction hidexp_parker

!***************************************************************************
    function hidexp_wu(ns,d,k,pbm,cm) result(psi)
! Calculates hiding and exposure coefficient based on Wu et al. (2000)
!
! Input:
!   ns - Number of sediment size classes
!   k - kth size class corresponding to the hiding and exposure function
!   d - Characteristic grain size diameters
!   cm - Hiding and exposure coefficient
!
! Output:
!   psi - Hiding and exposure function for kth size class
!
! Reference:
!   Wu, W., Wang, S.S.Y., and Jia, Y. (2000). "Nonuniform sediment transport
!     in alluvial rivers", Journal of Hydraulic Research, 38(6), 427-434.
!
! written by Alex Sanchez, USACE-CHL
!***************************************************************************
    implicit none
    !Input/Output
    integer,    intent(in) :: ns,k
    real(ikind),intent(in) :: d(ns),pbm(ns),cm
    real(ikind) :: psi
    !Internal variables
    integer :: j
    real(ikind) :: ph,pe,val

    pe = pbm(1)/(d(k)+d(1)) !Exposure probability
    ph = pe*d(1)            !Hiding probabability
    do j=2,ns
      val  = pbm(j)/(d(k)+d(j))
      pe = pe + val
      ph = ph + val*d(j)
    enddo !j 
    pe = pe*d(k)
    psi = (ph/pe)**cm  !Note mhe is positive

    return
    endfunction hidexp_wu
    
!***************************************************************************
    function hidexp_day(d16,d50,d84,d) result(etak)
! Hiding and exposure correction by Day (1980)
!
! Usage:
!   etak = hidexp_day(d16,d50,d84,d);
!
! Input:
!   d16 = 16th percentile diameter [L]
!   d50 = 50th percentile diameter [L]
!   d84 = 84th percentile diameter [L]
!   d = Characterestic grain size [L]
!
! Output:
!   etak = Hiding and exposure correction
!
! References:
!   Day, T.J. (1980). A study of the transport of graded sediments, 
!     HRS Wallingford, Report No. IT 190.
!
! Author: Alex Sanchez, USACE-CHL
!***************************************************************************
    implicit none
    real(ikind),intent(in):: d16,d50,d84,d
    real(ikind) :: dA,etak
    
    dA = d50*1.6*(d84/d16)**(-0.28)
    etak = 1.0/(0.4*(d/dA)**(-0.5)+0.6)

    return
    endfunction hidexp_day
    
!***************************************************************************    
    function hidexp_proffitt_sutherland(theta,d50,d) result(etak)
! Hiding and exposure correction factor by Proffitt and Sutherland (1983)
!
! Usage:
!  etak = hidexp_proffitt_sutherland(theta,d50,d)
! 
! Input:
!   theta = Shields number [-]
!   d50 = Median grain size [L]
!   d = Grain size [L]
!
! Output:
!   etak = Hiding and exposure correction factor [-]
!
! Author: Alex Sanchez, USACE-CHL
!***************************************************************************    
    implicit none
    real(ikind),intent(in):: theta,d50,d
    real(ikind) :: dratio,dadjust,r,etak
    
    if(theta<=0.04)then
      dratio = 1.1
    elseif(theta<=0.045)then
      dratio = 2.3-30*theta
    elseif(theta<=0.095)then
      dratio = 1.4-10*theta
    else 
      dratio = 0.45
    endif
    dadjust = d50*dratio
    r = d/dadjust
          
    if(r>=0.075 .and. r<3.7)then
      etak = 0.53*log10(r)+1.0
    elseif(r>=3.7)then
      etak = 1.3  
    else
        etak = 0.4
    endif

    return
    endfunction hidexp_proffitt_sutherland

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! End Hiding and Exposure Formula
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
! Begin Adaptation Coefficients
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

!******************************************************************************
    function adaptsusp_armanini_disilvio(sqrtgrav,h,d90,ws,ustar) result(alphas)
! Suspended-load adaptation coefficient based on Armanini and di Silvio (1986)
!
! Input:
!   sqrtgrav - sqrt(grav) 
!     where 
!       grav - gravity [L/T^2]
!   h - Water depth [L]
!   d90 - 90th percentile grain size [L]
!   ws - Sediment fall velocity [L/T]
!   ustar - Shear velocity [L/T]
!
! Output:
!    alphas - Suspended-load adapatation coefficient [-]
!
! References:
!  Armanini, A., and di Silvio, G. 1986. Discussion on the paper 
!   ‘A depth-averaged model for suspended sediment transport’,  
!   by Galappattti, G., and Vreugdenhil, C.B. Journal of Hydraulic Research. 
!   24(5), 437-441.
!
! Author: Alex Sanchez, USACE-CHL
!******************************************************************************
    implicit none
    !Input/Output
    real(ikind),intent(in) :: sqrtgrav,h,d90,ws,ustar
    real(ikind) :: alphas
    !Internal variables
    real(ikind) :: Ch,B,val
    
    Ch = 18.0*log(4.0*h/d90) !Chezy coefficient
    B = 33.0*h/exp(1.0+0.4*Ch/sqrtgrav)
    val = B/h
    alphas = 1.0/(val+(1.0-val)*exp(-1.5*val**(-0.1667)*ws/ustar))
    alphas = max(min(alphas,30.0),0.1) !Limit values
    
    return
    endfunction adaptsusp_armanini_disilvio
    
!!******************************************************************************
!    function armanini_disilvio_test(hi,z0i,wsi,ustar) result(alphas)
!! Suspended load adaptation coefficient based on Armanini and di Silvio (1986)
!!******************************************************************************
!    use flow_def, only: sqrtgrav
!    use prec_def
!    implicit none
!    !Input/Output
!    real(ikind),intent(in) :: hi,z0i,wsi,ustar
!    real(ikind) :: alphas
!    !Internal variables
!    real(ikind) :: val
!    
!    val = 33.0*z0i/hi
!    alphas = 1.0/(val+(1.0-val)*exp(-1.5*val**(-0.1667)*wsi/ustar))
!    alphas = max(min(alphas,30.0),0.1)
!    
!    return
!    endfunction armanini_disilvio_test
    
!******************************************************************* 
    function adaptsusp_lin(ws,ustar) result(alphas)
! Suspended-load adaptation coefficient based on Lin (1984)
!
! Input:
!   ws - Sediment fall velocity [L/T]
!   ustar - Shear velocity [L/T]
!
! Output:
!   alphas - Suspended-load adaptation coefficient [-]
!
! References:
!   Lin, B. N. 1984. Current study of unsteady transport of sediment 
!     in China. In Proceedings, Japan-China Bilateral Seminar on 
!     River Hydraulics and Engineering Experience, July, 337–342. 
!     Tokyo-Kyoto-Saporo, Japan.
!
! Author: Alex Sanchez, USACE-CHL
!*******************************************************************
    implicit none
    real(ikind),intent(in) :: ws,ustar
    real(ikind) :: alphas
    
    alphas = 3.25 + 0.55*log(ws/(0.4*ustar)) 
    alphas = max(min(alphas,30.0),0.1)  !Limit values
    
    return
    endfunction adaptsusp_lin 
    
!*******************************************************************
    function adaptsusp_gallappatti(hi,uvi,wsi,ustar) result(Ts)
! Suspended load adaptation time based on Gallappatti (1983)
!*******************************************************************
    implicit none
    !Input/Output
    real(ikind),intent(in) :: hi,uvi,wsi,ustar
    real(ikind) :: Ts
    !Internal variables
    real(ikind) :: wstar,val,ur
    
    ur = ustar/uvi
    wstar = wsi/ustar
    val = (1.547-20.12*ur)*wstar**3+(326.832*ur**2.2047-0.2)*wstar**2+&
          (0.1385*log(ur)-6.4061)*wstar+(0.5467*ur+2.1963)
    Ts = hi*exp(val)/ustar
    Ts = max(Ts,0.2) !Lower limit proposed by Xbeach
    
    return
    endfunction adaptsusp_gallappatti
    
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! End Adaptation Coefficients
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
! Begin Sediment Velocities
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
!*********************************************************************
    function bedvel_vanrijn(g,taubp,taucr,s,d) result(ub)
! Bed load velocity using formula of Van Rijn 1984 
! Note that grain roughness is determined using d50, not d
!
! Input:
!   g - Gravity [L/T^2]
!   taubp - Mean bed shear stress based on grain roughness [M/L/T^2]
!   taucr - Critical shear stress based on Shields diagram [M/L/T^2]
!   s - Sediment specific gravity [-]
!   d - Sediment diameter [L]
!
! Output:
!   ub - Bed-load velocity [L/T]
!
! References:
!   Van Rijn, L.C., 1984a. Sediment transport, Part I: Bed load transport.
!     Journal of Hydraulic Engineerin. 110(10), 1431–1456.
!
! Author: Alex Sanchez, USACE-CHL
!*********************************************************************
    implicit none
    !Input
    real(ikind),intent(in) :: g,taubp,taucr,s,d
    !Output
    real(ikind) :: ub
    !Internal
    real(ikind) :: stage
    
    stage = taubp/taucr-1.0 !Transport stage number [-]
    stage = max(stage,1.0e-15)
    ub = 1.5*sqrt((s-1.0)*g*d)*stage**0.6 !Van Rijn (1984)
    
    return
    endfunction bedvel_vanrijn

!*********************************************************************
    function bedvel_vanrijn_wu(g,taubp,taucr,s,d) result(ub)
! Bed load velocity using formula of Van Rijn 1984 with
! recalibrated coefficients of Wu et al. 2006
! Note that grain roughness is determined using d50, not d
!
! Input:
!   g - Gravity [L/T^2]
!   taubp - Mean bed shear stress based on grain roughness [M/L/T^2]
!   taucr - Critical shear stress based on Shields diagram [M/L/T^2]
!   s - Sediment specific gravity [-]
!   d - Sediment diameter [L]
!
! Output:
!   ub - Bed-load velocity [L/T]
!
! References:
!   Van Rijn, L.C., 1984a. Sediment transport, Part I: Bed load transport.
!     Journal of Hydraulic Engineerin. 110(10), 1431–1456.
!   Wu, W., Altinakar, M., Wang, S.S.Y., 2006. Depth-averaged 
!     analysis of hysteresis between flow and sediment transport 
!     under unsteady conditions. International Journal of Sediment 
!     Research. 21(2), 101–112.
!
! Author: Alex Sanchez, USACE-CHL
!*********************************************************************
    implicit none
    !Input
    real(ikind),intent(in) :: g,taubp,taucr,s,d
    !Output
    real(ikind) :: ub
    !Internal
    real(ikind) :: stage
    
    stage = taubp/taucr-1.0 !Transport stage number [-]
    stage = max(stage,1.0e-15)
    ub = 1.64*sqrt((s-1.0)*g*d*stage)
    
    return
    endfunction bedvel_vanrijn_wu
    
!********************************************************************
    function bedvel_fredsoe(rhow,taubp,taucr,taumaxp) result(ub)
! Bed-load velocity based on of Fredsoe and Deigaard (1992)
!
! Description:
!   Computes the bed-load velocity based on the formula of 
!   of Fredsoe and Deigaard (1992). 
!
! Input:
!   rhow - Water density [M/L^3]
!   taubp - Mean bed shear stress based on the grain roughness [M/L/T^2]
!   taucr - Critical shear stress based on the Shields diagram [M/L/T^2]
!   taumaxp - Maximum combined wave-current shear stress based on 
!         the grain roughness [M/L/T^2]
!
! Output:
!    ub - Bed-load velocity [L/T]
!
! References:
!   Fredsoe, J., and Deigaard, R. 1992. Mechanics of coastal 
!     sediment transport. World Scientific Publishing, Singapore. 
!
! Author: Alex Sanchez, USACE-CHL 
!********************************************************************
    implicit none
    real(ikind),intent(in) :: rhow,taubp,taucr,taumaxp
    real(ikind) :: ustarb,ub
    
    if(taumaxp>taucr)then
      ustarb = sqrt(taubp/rhow)
      ub = 10.0*ustarb*(1.0 - 0.7*sqrt(taubp/taumaxp))
    else
      ub = 0.0
    endif
    
    return
    endfunction bedvel_fredsoe
    
!******************************************************************
    function bedvel_greimann(g,thetap,thetaref,s,d) result(ub)
! Bed load velocity using formula of Greimann et al. (2008)
! Note that grain roughness is determined using d50, not d
!
! Input:
!   g - Gravity [L/T^2]
!   thetap - Mean bed Shields parameter based on grain roughness [-]
!   thetaref - Reference Shields parameter at which there is 
!     a low but measurable reference transport rate, as defined 
!     in Parker (1990) [-]. Griemann et al. (2008) estimated or 
!     calibrated it for each particle type and bed roughness.
!   s - Sediment specific gravity [-]
!   d - Sediment diameter [L]
!
! Output:
!   ub - Bed-load velocity [L/T]
!
! References:
!   Greimann, B., Lai, Y., and Huang, J. 2008. Two-dimensional  
!     total sediment load model equations. Journal of Hydraulics 
!     Engineering, 134, 1142-1146.
!
! Author: Alex Sanchez, USACE-CHL
!******************************************************************
    implicit none
    !Input
    real(ikind),intent(in) :: g,thetap,thetaref,s,d
    !Output
    real(ikind) :: ub
    !Internal
    real(ikind) :: phi
    
    phi = thetap/thetaref
    phi = max(phi,1.0e-15)
    ub = 1.1*sqrt((s-1.0)*g*d)*phi**0.67*(1.0-exp(-5.0*phi))
    
    return
    endfunction bedvel_greimann
    
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! End Sediment Velocities
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

endmodule sed_lib
