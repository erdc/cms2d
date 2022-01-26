!====================================================================================
module fric_lib
! Wall and Bottom Friction Library
! 
! Contains the following:
! ~ Wall friction ~
!     wall_gamma - Calculates the wall gamma coefficient based 
!                  on log-law function using an iterative procedure
!     wall_coef - Calculates the wall friction coefficient
!               based on a simplified log-law function
!
! ~ Wave-current bottom friction ~
!     fric_wavcur_init - Sets the wave-current bottom friction coefficients
!     fric_wavfricfac_intp - Log-scale interpolation of wave friction factor
!     fric_curdragcoef_intp - Log-scale interpolation of current drag coefficient
!     fric_wavcurmax - Calculates the wave-current maximum shear stress 
!                based on Soulsby (1995) 
!     fric_wavcurmax_soulsby - Calculates the wave-current maximum shear stress 
!               based on Soulsby (1997)
!     fric_wavcurmean_fac - Calculates the wave-current bottom friction enhancement factor
!                    based on Soulsby (1995)
!     fric_wavcurmean_data2 - Calculates the wave-current mean bed shear stress
!             based on Soulsby's (1995) two coefficient Data2 method
!     fric_wavcurmean_data13 - Calculates the wave-current mean bed shear stress
!             based on Soulsby's (1995) thirteen coefficient Data13 method
!     fric_wavcurmean_F84 - Calculates the wave-current mean bed shear stress 
!             based on Fredsoe (1984) and parameterized by Soulsby (1995)
!     fric_wavcurmean_HT91 - Calculates the combined wave-current mean bed shear stress 
!            based on Huynh-Thanh and Temperville (1991) and parameterized by Soulsby (1995)
!     fric_wavcurmean_DSK88 - Calculates the wave-current mean bed shear stress 
!            based on Davies et al. (1988) and parameterized by Soulbsy (1995)
!     fric_wavcurmean_GM79 - Calculates the wave-current mean bed shear stress    
!            based on Grant and Madsen (1979) and parameterized by Soulbsy (1995)
!     fric_wavcurmean_quad  - Calculates the wave-current mean bed shear stress 
!            based on a simple quadratic formula
!     fric_curmean_quad  - Calculates the current mean bed shear stress 
!            based on a simple quadratic formula
! ~ Bottom Roughness Conversion ~
!     fric_conv_drag2man - Converts a bottom drag coefficient to a Manning's coefficient
!     fric_conv_man2drag - Converts a Manning's coefficient to a bottom drag coefficient
!     fric_conv_drag2length - Converts a drag coefficient to a roughness height
!     fric_conv_length2drag - Converts a roughness height to a drag coefficient
!     fric_normapprough - Calculates the normalize bed apparent roughness
! ~ Bed forms ~
!   Currents
!     fric_dune_vanrijn - Calculates the bed dune height and length based on van Rijn (1984)
!     fric_ripple_current_soulsby - Calculates the current-related rippple length and height 
!                                   based on Soulsby (1997)
!     fric_ripple_current_soulsbywhitehouse - Calculates the current-related rippple length and height 
!                                   based on Soulsby and Whitehouse (2005)
!     fric_ripple_current_raudkivi - Calculates the current-related ripple height and length 
!                                   based on Raudkivi (1998)
!   Waves
!     fric_ripple_wave_vanrijn - Wave ripple dimensions based on van Rijn (1984,1989)
!     fric_ripple_wave_soulsbywhitehouse - Calculates the wave-related rippple length and height 
!                                   based on Soulsby and Whitehouse (2005)
! ~ Roughness Estimates ~
!     fric_rough_ripple -  Computes the bed roughness height due to ripples
!     fric_rough_dune - Computes the bed roughness height due to dunes based on van Rijn (1984)
!     fric_rough_grain - Estimates the grain-related roughness height
!     fric_rough_trans_nielson - Estimates the sediment transport related roughness length z0t
!                      based on Nielson (1992)
!     fric_rough_trans_wiberg - Estimates the sediment transport related roughness length z0t
!                      based on Wiberg and Rubin (1989)
!     fric_rough_trans_wilsonlog - Sediment Transport roughness length using Wilson formula 
!                                  and a log current velocity profile
!     fric_rough_trans_wilsonswart - Sediment Transport roughness length using Wilson formula 
!                                  and the Swart wave friction factor
!
! written by Alex Sanchez, USACE-CHL
!========================================================================================
    implicit none

contains

!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
! Begin Wall Friction
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

!**************************************************************************
    subroutine wall_gamma(vk,yp,y0,upar,gammawall)
! Calculates the wall gamma coefficient based 
! on log-law function using an iterative procedure
!
! Description:
!   The wall shear stress may be written as
!     tauwall = rhow*gammawall*uwall
!   where 
!     tauwall = wall shear stress [N/m^2]
!     rhow = water density [kg/m^3]
!     gammawall = coefficient
!     uwall = velocity vector parallel to wall [m/s]
!
!   The subroutine is valid for smooth to rough flow 
!   but requires iterations. An initial guess is made 
!   using a simple log-law expression for a rough wall.
!   For further details see Wu (2007)
!   Note that for convenience the definition here for 
!   gamma is slightly different from Wu (2007) in that 
!   it does not include the water density. 
!
! Input:
!   yp - Distance from center of cell P to wall [m]
!   upar - Magnitude of velocity parallel to wall [m/s]
!   ksn - Wall Nikradse roughness [m]
!
! Output:
!   gammawall = wall coefficient
!
! References:
!   Cebeci, T., and Bradshaw, P. (1977). "Momentum transfer 
!     in boundary layers". Hemisphere, Washington D.C.
!   Wu, W. (2007). "Computational River Dynamics". 
!     Taylor & Francis, 494 p.
!
! written by Alex Sanchez, USACE-CHL
!***********************************************************************
    use prec_def
    implicit none    
    !Input/Output
    real(ikind),intent(in) :: vk,yp,y0,upar
    real(ikind),intent(out) :: gammawall
    !Internal variables
    integer :: iterwall
    real(ikind) :: gammawall2,gammawall0,ustarwall,rks,ksplus,yplus
    real(ikind) :: delB,logksplus,Ewall,error
    
    gammawall0 = wall_coef(yp,y0)*upar !Use as initial guess
    gammawall = gammawall0
    rks = y0*30.0 !Wall Nikradse roughness
    do iterwall=1,10   
      ustarwall = sqrt(abs(gammawall*upar)) !Wall shear velocity, m/s
      yplus = ustarwall*yp/vk
      if(yplus>11.63)then
        ksplus = ustarwall*rks/vk
        if(ksplus<2.25)then
          delB = 0.0
        elseif(ksplus>=90.0)then
          delB = -3.3+2.5*log(ksplus) !B0-8.5=-3.3; B0=5.2  
        else
          logksplus = log(ksplus)
          delB = (-3.3+2.5*logksplus)*sin(0.4258*(logksplus-0.811)) !B0-8.5=-3.3; B0=5.2
        endif
        Ewall = exp(0.4*(5.2-delB)) !E=8.432 for smooth wall, B0=5.2
        gammawall2 = ustarwall*0.4/log(Ewall*max(yplus,1.0e-3))
      else
        gammawall2 = vk/yp
      endif
      if(gammawall2<0.0 .or. gammawall2>1.0)then !Did not converge
        gammawall = gammawall0
        exit  
      endif          
      error=abs(gammawall2-gammawall)/gammawall !Relative error change
      if(error<0.001) exit      
      gammawall = gammawall2
    enddo
    
    return
    endsubroutine wall_gamma
    
!**********************************************************
    function wall_coef(yp,y0) result(cwall)
! Calculates the wall friction coefficient using a 
! simplified log-law function
!
! Description:
!  The wall friction coefficient (cwall) is defined as
!    tau_wall = rhow*cwall*|upar|*upar
!  where  
!    rhow = water density [kg/m^3]
!    upar = wall-tangent velocity vector [m/s]
!
!  Valid for a rough flow only.
!
! Input:
!   yp = Distance from velocity point to wall [m]
!   y0 = Wall roughness length [m]
!
! Output:
!    cwall = wall friction coefficient
!
! written by Alex Sanchez, USACE-CHL
!**********************************************************
    use prec_def
    implicit none
    real(ikind),intent(in) :: yp,y0
    real(ikind) :: cwall,y0lim
    
    !Addressing an issue of the variables 'yp' and 'y0lim' having exactly the same value and therefore log(1) = 0 which causes a divide by zero.  MEB  01/10/2022
    y0lim = min(max(y0,1.0e-6),yp*0.9999)      
    cwall = (0.4/log(yp/y0lim))**2
    
    return
    end function wall_coef
    
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! End Wall Friction
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
! Begin Wave-current Bottom Friction
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<    

!*******************************************************************
    subroutine fric_wavcur_init
! Sets the wave-current bottom friction coefficients 
! Coefficient values taken from 
!   Soulsby, R.L. (1997). "Dynamics of marine sands", 
!      Thomas Telford, 249 p.
! 
! written by Alex Sanchez, USACE-CHL
!*******************************************************************
    use prec_def
    use fric_def
    implicit none                                
        
    selectcase(mwavcurint)
      case(3) !DATA13
        !Mean shear stess coefficients
        bwc = (/ 0.47, 0.69,-0.09,-0.08/)
        pwc = (/-0.53, 0.47, 0.07,-0.02/)    
        qwc = (/ 2.34,-2.41, 0.45,-0.61/)    
        cjwc = 8.80
        
      case(4) !HT91         
        !Wave friction factor coefficients
        logfw(1)=log(0.0750) !for A/z0 = 1.0e2
        logfw(2)=log(0.0272) !for A/z0 = 1.0e3
        logfw(3)=log(0.0121) !for A/z0 = 1.0e4
        logfw(4)=log(0.0062) !for A/z0 = 1.0e5
    
        !Current drag coefficients
        logCd(1)=log(0.01231) !for z0/h = 1.0e-2
        logCd(2)=log(0.00482) !for z0/h = 1.0e-3
        logCd(3)=log(0.00237) !for z0/h = 1.0e-4
        logCd(4)=log(0.00141) !for z0/h = 1.0e-5
    
        !Maximum shear stress coefficients
        awc =  (/-0.07, 1.87,-0.34,-0.12/)
        cmwc = (/ 0.72,-0.33, 0.08, 0.34/)
        cnwc = (/ 0.78,-0.23, 0.12,-0.12/)    
        ciwc = 0.82
    
        !Mean shear stess coefficients
        bwc = (/ 0.27, 0.51,-0.10,-0.24/)
        pwc = (/-0.75, 0.13, 0.12, 0.02/)    
        qwc = (/ 0.89, 0.40, 0.50,-0.28/)
        cjwc = 2.70
        
      case(5) !F84
        !Wave friction factor coefficients
        logfw(1)=log(0.0592)
        logfw(2)=log(0.0221)
        logfw(3)=log(0.0102)
        logfw(4)=log(0.0056)    
    
        !Maximum shear stress coefficients
        awc =  (/-0.06, 1.70,-0.29, 0.29/)
        cmwc = (/ 0.67,-0.29, 0.09, 0.42/)
        cnwc = (/ 0.75,-0.27, 0.11,-0.02/)
        ciwc = 0.80
    
        !Mean shear stess coefficients
        bwc = (/ 0.29,0.55,-0.10,-0.14/)
        pwc = (/-0.77,0.10, 0.27, 0.14/)
        qwc = (/ 0.91,0.25, 0.50, 0.45/)
        cjwc = 3.00
        
      case(6) !DSK88
        !Wave friction factor coefficients
        logfw(1)=log(0.0701)
        logfw(2)=log(0.0209)
        logfw(3)=log(0.0120)
        logfw(4)=log(0.00661)    
    
        !Maximum shear stress coefficients
        awc =  (/0.05, 1.62,-0.38, 0.25/)
        cmwc = (/1.05,-0.72,-0.08, 0.59/)
        cnwc = (/0.66,-0.25, 0.19,-0.03/)
        ciwc = 0.82
    
        !Mean shear stess coefficients
        bwc = (/ 0.22, 0.73,-0.05,-0.35/)
        pwc = (/-0.86, 0.26, 0.34,-0.07/)
        qwc = (/-0.89, 2.33, 2.60,-2.50/)
        cjwc = 2.70
      
      case(7) !GM79
        !Wave friction factor coefficients
        logfw(1)=log(0.1057)
        logfw(2)=log(0.0316)
        logfw(3)=log(0.0135)
        logfw(4)=log(0.0069)    
    
        !Maximum shear stress coefficients
        awc =  (/0.11, 1.95,-0.49,-0.28/)
        cmwc = (/0.65,-0.22, 0.15, 0.06/)
        cnwc = (/0.71,-0.19, 0.17,-0.15/)
        ciwc = 0.67
    
        !Mean shear stess coefficients
        bwc = (/ 0.73, 0.40,-0.23,-0.24/)
        pwc = (/-0.68, 0.13, 0.24,-0.07/)
        qwc = (/ 1.04,-0.56, 0.34,-0.27/)
        cjwc = 0.5
    
    endselect
    
    return
    endsubroutine fric_wavcur_init

!****************************************************************
    function fric_wavfricfac_intp(z0,Uw,Tw) result(fw)
! Log-scale interpolation of wave friction factor based on table
!****************************************************************    
    use fric_def, only: logfw
    use const_def, only: twopi
    use prec_def
    implicit none
    integer :: i
    real(ikind),intent(in) :: z0,Uw,Tw
    real(ikind) :: fw,logAz0i,fac
    
    logAz0i = log10(Uw*Tw/(twopi*z0))
    i = int(logAz0i) !truncates towards zero
    i = min(max(i,2),4)
    fac = real(i)-logAz0i+1.0
    fw = exp(fac*logfw(i-1)+(1.0-fac)*logfw(i))       
    fw = min(fw,0.3)
    
    return
    endfunction fric_wavfricfac_intp

!**********************************************************************
    function fric_curdragcoef_intp(h,z0) result(Cd)
! Log-scale interpolation of current drag coefficient based on a table   
!**********************************************************************    
    use fric_def, only: logCd
    use prec_def
    implicit none
    integer :: i
    real(ikind),intent(in) :: h,z0
    real(ikind) :: logz0hi,fac,Cd
    
    logz0hi = log10(z0/max(h,1e-5))
    i = -aint(logz0hi)-1
    i = min(max(i,2),4)
    fac = real(i)+logz0hi+1.0
    Cd = exp(fac*logCd(i-1)+(1.0-fac)*logCd(i))
    
    return
    endfunction fric_curdragcoef_intp

!**************************************************************************
    function fric_wavcurmax_soulsby(taub,tauw,phi) result(taumax)
! Calculates the wave-current maximum shear stress following Soulsby (1997)
!
! Input:
!   taub - Mean wave-current shear stress [M/L/T^2]
!   tauw - Wave shear stress [M/L/T^2]
!   phi - Wave-current angle [rad]
!
! Output:
!  taumax = Maximum combined wave-current shear stress [M/L/T^2]
!
! Author: Alex Sanchez, USACE-CHL
!**************************************************************************
    use prec_def
    implicit none
    real(ikind),intent(in) :: taub,tauw,phi
    real(ikind) :: taumax
    
    taumax = sqrt((taub+tauw*cos(phi))**2+(tauw*sin(phi))**2) !Soulsby 1997
    
    return
    endfunction fric_wavcurmax_soulsby
    
!**************************************************************************
    function fric_wavcurmean_fac(tauc,tauw,Cd,fw,phi) result(fac)
! Calculates the wave-current bottom friction enhancement factor
! based on Soulsby et.al. (1993) and Soulsby 1995. 
! The wave-current bottom friction enhancement factor (fac) is defined as:
!    taucw = fac*tauc
!  where
!    tauc = current bed shear stress
!    taucw = combined wave-current mean bed shear stres
!
! Input:      
!   tauc = current bed shear stress [N/m^2]
!   tauw = wave bed shear stress [N/m^2]
!   Cd = current drag coefficient [-]
!   fw = wave friction factor    [-]
!   phi = wave-current angle [rad]
!
! Output:
!   fac = wave-current bottom friction enhancement factor [-]
!
! References:
!  Soulsby, R.L., Hamm, L., Klopman, G., Myrhaug, D., Simons, R.R., 
!    and Thomas, G.P., 1993. Wave–current interactions within and 
!    outside the bottom boundary layer. Coastal Engineering, 21, 41–69. 
!  Soulsby, R.L. 1995. "Bed shear stresses due to combined waves
!    and currents", in Advances in Coastal Morphodynamics, 
!    ed. M.J.F Stive, and H.J. de Vriend, J. Fredsoe, L. Hamm, 
!    R.L. Soulsby, C. Teisson, and J.C. Winterwerp, 
!    Deflt Hydraulics, Netherlands, pp. 4-20 to 4-23. 
!
! written by Alex Sanchez, USACE-CHL
!**************************************************************************
    use fric_def
    use prec_def
    implicit none
    real(ikind),intent(in) :: tauc,tauw,Cd,fw,phi
    real(ikind) :: b,p,q,xval,cosphij,log10fwCd,fac
    
    if(tauc<1.e-15)then
      fac = 1.0
      return
    endif
    cosphij = abs(cos(phi))**cjwc    
    log10fwCd = log10(fw/Cd)    
    b = bwc(1)+bwc(2)*cosphij+(bwc(3)+bwc(4)*cosphij)*log10fwCd
    p = pwc(1)+pwc(2)*cosphij+(pwc(3)+pwc(4)*cosphij)*log10fwCd
    q = qwc(1)+qwc(2)*cosphij+(qwc(3)+qwc(4)*cosphij)*log10fwCd    
    xval = tauc/(tauc+tauw)            
    fac = 1.0+b*(xval**p)*(1.0-xval)**q !enhancement factor
    
    return
    endfunction fric_wavcurmean_fac
    
!********************************************************************
    function fric_wavcurmax(tauc,tauw,Cd,fw,phi) result(taubmax)
! Calculates the wave-current maximum shear stress 
! based on Soulsby (1995) 
!
! Input:      
!   tauc - current shear stress [N/m^2]
!   tauw - wave shear stress [N/m^2]   
!   Cd - current drag coefficient [-]
!   fw - wave friction factor [-]
!   phi - wave-current angle [rad]
!
! Output:
!   taumax - wave-current maximum shear stress [N/m^2]
!
! References:
!  Soulsby, R.L. 1995. Bed shear stresses due to combined waves
!    and currents. in Advances in Coastal Morphodynamics, 
!    ed. M.J.F Stive, and H.J. de Vriend, J. Fredsoe, L. Hamm, 
!    R.L. Soulsby, C. Teisson, and J.C. Winterwerp, 
!    Deflt Hydraulics, Netherlands, pp. 4-20 to 4-23. 
!
! written by Alex Sanchez, USACE-CHL
!*********************************************************************
    use fric_def
    use prec_def
    implicit none    
    real(ikind),intent(in) :: tauc,tauw,Cd,fw,phi
    real(ikind) :: a,cm,cn,xval,cosphii,log10fwCd,taubmax
    
    cosphii = abs(cos(phi))**ciwc
    log10fwCd = log10(fw/Cd)       
    a  =  awc(1)+ awc(2)*cosphii+( awc(3)+ awc(4)*cosphii)*log10fwCd
    cm = cmwc(1)+cmwc(2)*cosphii+(cmwc(3)+cmwc(4)*cosphii)*log10fwCd
    cn = cnwc(1)+cnwc(2)*cosphii+(cnwc(3)+cnwc(4)*cosphii)*log10fwCd            
    xval = tauc/(tauc+tauw) 
    taubmax = (1.0+a*(xval**cm)*(1.0-xval)**cn)*(tauc+tauw)
    
    return
    endfunction fric_wavcurmax

!****************************************************************************
    subroutine fric_wavcurmean_data2(rhow,Cd,z0,Vx,Vy,us,vs,Uw,Tw,Um,Uwc,taum)
! Calculates the wave-current mean bed shear stress based on
! Soulsby's (1995) two coefficient Data2 Method
!
! Input:
!   rhow = water density [kg/m^3]
!   Cd = bottom drag coefficient [-]
!   z0 = bottom roughness length [m]
!   Vx,Vy = Lagrangian depth-averaged current velocities [m/s]
!   us,vs = Stokes velocities [m/s]
!   Uw = wave bottom orbital velocity amplitude [m/s]
!        (for random waves Uw = sqrt(2)*Urms)
!   Tw = wave period [s]
!        (for random waves Tw is the peak wave period)
!
! Output:
!   Um = Eulerian current velocity magnitude [m/s]
!   Uwc = wave-current velocity [m/s]
!   taum = mean bed shear stress [N/m^2]
!
! References:
!  Soulsby, R.L. (1995). "Bed shear stresses due to combined waves
!    and currents", in Advances in Coastal Morphodynamics, 
!    ed. M.J.F Stive, and H.J. de Vriend, J. Fredsoe, L. Hamm, 
!    R.L. Soulsby, C. Teisson, and J.C. Winterwerp, 
!    Deflt Hydraulics, Netherlands, pp. 4-20 to 4-23. 
!
! written by Alex Sanchez, USACE-CHL
!****************************************************************************
    use const_def, only: pi,twopi
    use prec_def
    implicit none    
    !Input
    real(ikind), intent(in) :: rhow,Cd,Vx,Vy,z0,us,vs,Uw,Tw
    !Output
    real(ikind), intent(out) :: Um,Uwc,taum
    !Internal variables
    real(ikind) :: Um2,fw,tauc,tauw,fac,r,Aw
           
    if(Uw>0.001)then   
      Um2 = (Vx-us)**2+(Vy-vs)**2 !Eulerian velocity squared  
      Um = sqrt(Um2)
      tauc = rhow*Cd*Um2 !Current shear [Pa]
      Aw = Uw*Tw/twopi   !Amplitude of near-bed wave exursion [m]
      r = Aw/(30.0*z0)   !Relative roughness [-]
      fw = fric_wavefac_soulsby(r) !Wave friction factor [-]
      tauw = 0.5*rhow*fw*Uw*Uw    !wave friction [Pa]
      fac = 1.0+1.2*(tauw/(tauc+tauw))**3.2 !Enhancement factor  
      taum = tauc*fac
      Uwc = Um*fac
    else       
      Um2 = Vx*Vx+Vy*Vy   
      Um = sqrt(Um2)     
      taum = rhow*Cd*Um2 !current friction
      Uwc = Um
    endif
    
    return
    endsubroutine fric_wavcurmean_data2
          
!********************************************************************************
    subroutine fric_wavcurmean_data13(rhow,Cd,z0,Vx,Vy,us,vs,Uw,Tw,Dw,Um,Uwc,taum)
! Calculates the wave-current mean bed shear stress based on
! Soulsby's (1995) thirteen coefficient Data13 method
!
! Input:
!   rhow = water density [kg/m^3]
!   Cd = bottom drag coefficient [-]
!   z0 = bottom roughness length [m]
!   Vx,Vy = Lagrangian depth-averaged current velocities [m/s]
!   us,vs = stokes velocities [m/s]
!   Uw = wave bottom orbital velocity [m/s]
!        (for random waves Uw = sqrt(2)*Urms)
!   Tw = wave period [s]
!        (for random waves Tw = Tp)
!   Dw = wave direction [rad]
!
! Output:
!   Um = Eulerian current velocity magnitude [m/s]
!   Uwc = wave-current velocity [m/s]
!   taum = mean bed shear stress [N/m^2]
!
! References:
!  Soulsby, R.L. (1995). "Bed shear stresses due to combined waves
!    and currents", in Advances in Coastal Morphodynamics, 
!    ed. M.J.F Stive, and H.J. de Vriend, J. Fredsoe, L. Hamm, 
!    R.L. Soulsby, C. Teisson, and J.C. Winterwerp, 
!    Deflt Hydraulics, Netherlands, pp. 4-20 to 4-23. 
!
! written by Alex Sanchez, USACE-CHL
!********************************************************************************
    use const_def, only: twopi
    use prec_def
    implicit none    
    !Input
    real(ikind), intent(in) :: rhow,Cd,Vx,Vy,us,vs,Uw,Tw,Dw,z0
    !Output
    real(ikind), intent(out) :: Um,Uwc,taum
    !Internal variables
    real(ikind) :: fw,Um2,phi,tauc,tauw,fac,Aw,r !fric_wavcurmean_fac
           
    if(Uw>0.001)then    
      Um2 = (Vx-us)**2+(Vy-vs)**2 !Eulerian velocity squared  
      Um = sqrt(Um2)
      tauc = rhow*Cd*Um2 !Current stress [Pa]
      Aw = Uw*Tw/twopi   !Amplitude of near-bed wave exursion [m]
      r = Aw/(30.0*z0)   !Relative roughness [-]
      fw = fric_wavefac_soulsby(r) !Wave friction factor [-]
      tauw = 0.5*rhow*fw*Uw*Uw   !Maximum wave stress [Pa]
      phi = abs(Dw-atan2(Vy,Vx)) !Current-wave angle [rad]
      fac = fric_wavcurmean_fac(tauc,tauw,Cd,fw,phi)   
      taum = tauc*fac
      Uwc = Um*fac
    else      
      Um2 = Vx*Vx+Vy*Vy  
      Um = sqrt(Um2)      
      taum = rhow*Cd*Um2 !current friction
      Uwc = Um
    endif
    
    return
    endsubroutine fric_wavcurmean_data13

!*************************************************************************
    subroutine fric_wavcurmean_F84(rhow,Cd,z0,Vx,Vy,us,vs,Uw,Tw,Dw,Um,Uwc,taum)
! Calculates the wave-current mean bed shear stress based on
! Fredsoe (1984) and parameterized by Soulsby (1995)
!
! Input:
!   rhow = water density [kg/m^3]
!   Cd = bottom drag coefficient [-]
!   z0 = bottom roughness length [m]
!   Vx,Vy = Lagrangian depth-averaged current velocities [m/s]
!   us,vs = stokes velocities [m^2/s]
!   Uw = wave bottom orbital velocity [m/s]
!        (for random waves Uw = sqrt(2)*Urms)
!   Tw = wave period [s]
!        (for random waves Tw = Tp)
!   Dw = wave direction [rad]
!
! Output:
!   Um = Eulerian current velocity magnitude [m/s]
!   Uwc = wave-current velocity [m/s]
!   taum = mean bed shear stress [N/m^2]
!
! References:
!  Fredsoe, J. (1984). "Turbulent boundary layer in wave-current motion",
!    Journal of Hydraulic Engineering, 110, 1103-1120.
!  Soulsby, R.L. (1995). "Bed shear stresses due to combined waves
!    and currents", in Advances in Coastal Morphodynamics, 
!    ed. M.J.F Stive, and H.J. de Vriend, J. Fredsoe, L. Hamm, 
!    R.L. Soulsby, C. Teisson, and J.C. Winterwerp, 
!    Deflt Hydraulics, Netherlands, pp. 4-20 to 4-23. 
!
! written by Alex Sanchez, USACE-CHL
!*************************************************************************
    use prec_def
    implicit none
    real(ikind), intent(in) :: rhow,Cd,Vx,Vy,z0,us,vs,Uw,Tw,Dw
    real(ikind), intent(out) :: Um,Uwc,taum
    real(ikind) :: Um2,tauc,tauw,fw,fac,phi
    
    if(Uw>0.001)then      
      Um2 = (Vx-us)**2+(Vy-vs)**2 !Eulerian velocity squared  
      Um = sqrt(Um2)
      tauc = rhow*Cd*Um2       !current shear 
      fw = fric_wavfricfac_intp(z0,Uw,Tw)  
      tauw = 0.5*rhow*fw*Uw*Uw !wave shear   
      phi = abs(Dw-atan2(Vy,Vx))     !Current-wave angle, radians
      fac = fric_wavcurmean_fac(tauc,tauw,Cd,fw,phi)
      taum = tauc*fac
      Uwc = Um*fac
    else            
      Um2 = Vx*Vx+Vy*Vy  
      Um = sqrt(Um2)      
      taum = rhow*Cd*Um2 !current friction
      Uwc = Um
    endif
    
    return
    endsubroutine fric_wavcurmean_F84
    
!******************************************************************************
    subroutine fric_wavcurmean_HT91(rhow,h,z0,Vx,Vy,us,vs,Uw,Tw,Dw,Um,Uwc,taum)
! Calculates the wave-current mean bed shear stress based on
! Huynh-Thanh and Temperville (1991) and parameterized by Soulsby (1995)
!
! Input:
!   rhow = water density [kg/m^3]
!   h = wave-averaged total water depth [m]
!   Cd = bottom drag coefficient [-]
!   z0 = bottom roughness height [m]
!   Vx,Vy = Lagrangian depth-averaged current velocities [m/s]
!   us,vs = stokes velocities [m/s]
!   Uw = wave bottom orbital velocity [m/s]
!       (for random waves Uw = sqrt(2)*Urms)
!   Tw = wave period [s]
!       (for random waves Tw = Tp)
!
! Output:
!   Um = Eulerian current velocity magnitude [m/s]
!   Uwc = wave-current velocity [m/s]
!   taum = mean bed shear stress [N/m^2]
!
! References:
!  Huynh-Thanh, S., and Temperville, A. 1991. A numerical model
!    of the rough turbulent boundary layer in combined wave and
!    current interaction. in Sand Transport in Rivers, Estuaries, 
!    and the Sea, ed. R.L. Soulsby and R. Bettess, 
!    Balkema, Rotterdam. pp. 93-100.
!  Soulsby, R.L. 1995. Bed shear stresses due to combined waves
!    and currents. in Advances in Coastal Morphodynamics, 
!    ed. M.J.F Stive, and H.J. de Vriend, J. Fredsoe, L. Hamm, 
!    R.L. Soulsby, C. Teisson, and J.C. Winterwerp, 
!    Deflt Hydraulics, Netherlands, pp. 4-20 to 4-23. 
!
! written by Alex Sanchez, USACE-CHL
!******************************************************************************
    use prec_def
    implicit none
    real(ikind), intent(in) :: rhow,h,Vx,Vy,z0,us,vs,Uw,Tw,Dw
    real(ikind), intent(out) :: Um,Uwc,taum
    real(ikind) :: Um2,tauc,tauw,phi,fac,fw,Cd
    
    if(Uw>0.001)then
      Cd = fric_curdragcoef_intp(h,z0)
      Um2 = (Vx-us)**2+(Vy-vs)**2 !Eulerian velocity squared  
      Um = sqrt(Um2)
      tauc = rhow*Cd*Um2 !current shear    
      fw = fric_wavfricfac_intp(z0,Uw,Tw)  
      tauw = 0.5*rhow*fw*Uw*Uw    !wave shear
      phi = abs(Dw-atan2(Vy,Vx))     !Current-wave angle, radians
      fac = fric_wavcurmean_fac(tauc,tauw,Cd,fw,phi)
      taum = tauc*fac !combined mean shear
      Uwc = Um*fac
    else            
      Um2 = Vx*Vx+Vy*Vy  
      Um = sqrt(Um2)  
      Cd = fric_curdragcoef_intp(h,z0)         
      taum = rhow*Cd*Um2 !current friction
      Uwc = Um
    endif    
    
    return
    endsubroutine fric_wavcurmean_HT91

!********************************************************************************
    subroutine fric_wavcurmean_DSK88(rhow,Cd,z0,Vx,Vy,us,vs,Uw,Tw,Dw,Um,Uwc,taum)
! Calculates the wave-current mean bed shear stress based on
! Davies et al. (1988) and parameterized by Soulbsy (1995)
!
! Input:
!   rhow = water density [kg/m^3]
!   Cd = bottom drag coefficient [-]
!   z0 = bottom roughness length [m]
!   Vx,Vy = Lagrangian depth-averaged current velocities [m/s]
!   us,vs = stokes velocities [m^2/s]
!   Uw = wave bottom orbital velocity [m/s]
!        (for random waves Uw = sqrt(2)*Urms)
!   Tw = wave period [s]
!        (for random waves Tw = Tp)
!   Dw = wave direction, [rad]
!
! Output:
!   Um = Eulerian current velocity magnitude [m/s]
!   Uwc = wave-current velocity [m/s]
!   taum = mean bed shear stress [N/m^2]
!
! References:
!  Davies, A.G., Soulbsy, R.L., and King, H.L. (1988). "A numerical
!    model of the combined wave and current bottom boundary layer",
!    Journal of Geophysical Research, 93(C1), 491-508.
!  Soulsby, R.L. (1995). "Bed shear stresses due to combined waves
!    and currents", in Advances in Coastal Morphodynamics, 
!    ed. M.J.F Stive, and H.J. de Vriend, J. Fredsoe, L. Hamm, 
!    R.L. Soulsby, C. Teisson, and J.C. Winterwerp, 
!    Deflt Hydraulics, Netherlands, pp. 4-20 to 4-23. 
!
! written by Alex Sanchez, USACE-CHL
!*************************************************************************
    use prec_def
    implicit none
    real(ikind), intent(in) :: rhow,Cd,Vx,Vy,z0,us,vs,Uw,Tw,Dw
    real(ikind), intent(out) :: Um,Uwc,taum
    real(ikind) :: Um2,tauc,tauw,fw,fac,phi
    
    if(Uw>0.001)then      
      Um2 = (Vx-us)**2+(Vy-vs)**2 !Eulerian velocity squared  
      Um = sqrt(Um2)
      tauc = rhow*Cd*Um2 !current shear 
      fw = fric_wavfricfac_intp(z0,Uw,Tw)  
      tauw = 0.5*rhow*fw*Uw*Uw   !wave shear   
      phi = abs(Dw-atan2(Vy,Vx))     !Current-wave angle, radians
      fac = fric_wavcurmean_fac(tauc,tauw,Cd,fw,phi)
      taum = tauc*fac
      Uwc = Um*fac
    else            
      Um2 = Vx*Vx+Vy*Vy  
      Um = sqrt(Um2)      
      taum = rhow*Cd*Um2 !current friction
      Uwc = Um
    endif
    
    return
    endsubroutine fric_wavcurmean_DSK88

!******************************************************************************
    subroutine fric_wavcurmean_GM79(rhow,Cd,z0,Vx,Vy,us,vs,Uw,Tw,Dw,Um,Uwc,taum)
! Calculates the wave-current mean bed shear stress based on
! Grant and Madsen (1979) and parameterized by Soulbsy (1995)
!
! Input:
!   rhow - water density [kg/m^3]
!   Cd - bottom drag coefficient [-]
!   z0 - bottom roughness length [m]
!   Vx,Vy - Lagrangian depth-averaged current velocities [m/s]
!   us,vs - stokes velocities [m^2/s]
!   Uw - wave bottom orbital velocity [m/s]
!        (for random waves Uw = sqrt(2)*Urms)
!   Tw - wave period [s]
!        (for random waves use Tw = Tp)
!   Dw - wave direction [rad]
!
! Output:
!   Um - Eulerian current velocity magnitude [m/s]
!   Uwc - wave-current velocity [m/s]
!   taum - mean bed shear stress [N/m^2]
!
! References:
!  Grant, W.D., and Madsen, O.S. (1979). "Combined wave and current
!    interaction with a rough bottom", Journal of Geophysical Research,
!    84(C4), 1797-1808.
!  Soulsby, R.L. (1995). "Bed shear stresses due to combined waves
!    and currents", in Advances in Coastal Morphodynamics, 
!    ed. M.J.F Stive, and H.J. de Vriend, J. Fredsoe, L. Hamm, 
!    R.L. Soulsby, C. Teisson, and J.C. Winterwerp, 
!    Deflt Hydraulics, Netherlands, pp. 4-20 to 4-23. 
!
! written by Alex Sanchez, USACE-CHL
!******************************************************************************
    use prec_def
    implicit none
    real(ikind), intent(in) :: rhow,Cd,Vx,Vy,z0,us,vs,Uw,Tw,Dw
    real(ikind), intent(out) :: Um,Uwc,taum
    real(ikind) :: Um2,tauc,tauw,fw,fac,phi
    
    if(Uw>0.001)then      
      Um2 = (Vx-us)**2+(Vy-vs)**2 !Eulerian velocity squared  
      Um = sqrt(Um2)
      tauc = rhow*Cd*Um2 !current shear 
      fw = fric_wavfricfac_intp(z0,Uw,Tw)  
      tauw = 0.5*rhow*fw*Uw*Uw   !wave shear   
      phi = abs(Dw-atan2(Vy,Vx))     !Current-wave angle, radians
      fac = fric_wavcurmean_fac(tauc,tauw,Cd,fw,phi)
      taum = tauc*fac
      Uwc = Um*fac
    else            
      Um2 = Vx*Vx+Vy*Vy  
      Um = sqrt(Um2)      
      taum = rhow*Cd*Um2 !current friction
      Uwc = Um
    endif
    
    return
    endsubroutine fric_wavcurmean_GM79

!*********************************************************************
    subroutine fric_wavcurmean_quad(rhow,Cd,Vx,Vy,us,vs,Uw,Um,Uwc,taum,ustar)
! Calculates the wave-current mean bed shear stress based on
! a simple quadratic formula
!
! Input:
!   rhow = water density [kg/m^3]
!   Cd = bottom drag coefficient [-]
!   Vx,Vy = Lagrangian depth-averaged current velocities [m/s]
!   us,vs = stokes velocities [m^2/s]
!   Uw = wave bottom orbital velocity [m/s]
!        (for random waves Uw = sqrt(2)*Urms)
!
! Output:
!   Um = Eulerian current velocity magnitude [m/s]
!   Uwc = wave-current velocity [m/s]
!   taum = mean bed shear stress [N/m^2]
!   ustar = bed shear velocity [m/s]
!
! written by Alex Sanchez, USACE-CHL
!*********************************************************************
    use fric_def, only: cfricwav 
    use prec_def
    implicit none
    real(ikind), intent(in) :: rhow,Cd,Vx,Vy,us,vs,Uw
    real(ikind), intent(out) :: Um,Uwc,taum,ustar    
    real(ikind) :: Um2
    
    Um2 = (Vx-us)**2+(Vy-vs)**2 !Eulerian velocity squared  
    Um = sqrt(Um2)
    Uwc = sqrt(Um2+cfricwav*Uw*Uw)
    taum = rhow*Cd*Um*Uwc
    ustar = sqrt(taum/rhow+1.0e-15)
    
    return
    endsubroutine fric_wavcurmean_quad

!*********************************************************************
    subroutine fric_curmean_quad(rhow,Cd,Ux,Uy,Um,Uwc,taum,ustar)
! Calculates the current mean bed shear stress based on
! a simple quadratic formula
!
! Input:
!   rhow = water density [kg/m^3]
!   Cd = bottom drag coefficient [-]
!   Ux,Uy = Eulerian depth-averaged current velocities [m/s]
!
! Output:
!   Um = Eulerian current velocity magnitude [m/s]
!   Uwc = wave-current velocity [m/s]
!   taum = mean bed shear stress [N/m^2]
!   ustar = bed shear velocity [m/s]
!
! written by Alex Sanchez, USACE-CHL
!*********************************************************************    
    use prec_def
    implicit none
    real(ikind), intent(in) :: rhow,Cd,Ux,Uy
    real(ikind), intent(out) :: Um,Uwc,taum,ustar    
    
    Um = sqrt(Ux*Ux+Uy*Uy)
    Uwc = Um  
    taum = rhow*Cd*Um*Um
    ustar = sqrt(taum/rhow+1.0e-15)
    
    return
    endsubroutine fric_curmean_quad
    
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! End Wave-current Bottom Friction    
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
! Begin Bottom Roughness Conversions    
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<    
        
!*****************************************************************    
    function fric_conv_drag2man(g,Cd,h,hmax) result(cman)
! Converts a bottom drag coefficient to a Manning's coefficient    
!
! Usage:
!   cman = fric_conv_drag2man(g,Cd,h,hmax)
!
! Input:
!   g - Gravity [L/T^2]
!   Cd - Drag coefficient [-]
!   h - Water depth [L]
!   hmax - Maximum water depth [L]
!
! Output:
!   cman - Manning's roughness coefficient [T/L^(1/3)]
!
! Author: Alex Sanchez, USACE-CHL
!*****************************************************************
    use prec_def
    implicit none
    real(ikind),intent(in) :: g,Cd,h,hmax
    real(ikind) :: cman
    
    cman = sqrt(Cd*((min(hmax,h))**0.3333)/g)
    
    return
    endfunction fric_conv_drag2man
    
!***************************************************************
    function fric_conv_man2drag(g,cman,h,hmax) result(Cd)
! Converts a Manning's coefficient to a bottom drag coefficient
!
! Usage:
!   Cd = fric_conv_man2drag(cman,h,hmax)
!
! Input:
!   g - graviational constant [L/T^2]
!   cman - Manning's roughness coefficient [T/L^(1/3)]
!   h - Water depth [L]
!
! Output:
!   z0 - roughness length [m]
!
! Author: Alex Sanchez, USACE-CHL
!***************************************************************
    use prec_def
    implicit none
    real(ikind),intent(in) :: g,cman,h,hmax
    real(ikind) :: Cd,hh
    
    hh = max(min(hmax,h),0.1*cman)
    Cd = g*cman*cman/(hh**0.3333)
    
    return
    endfunction fric_conv_man2drag
    
!*****************************************************    
    function fric_conv_drag2length(Cd,h) result(z0)
! Converts a drag coefficient to a roughness height
!
! Usage:
!   z0 = ric_conv_drag2length(Cd,h)
!
! Input:
!   Cd - Drag coefficient [-]
!   h - Water depth [L]
!
! Output:
!   z0 - roughness length [L]
!
! Author: Alex Sanchez, USACE-CHL
!*****************************************************
    use prec_def
    implicit none
    real(ikind),intent(in) :: Cd,h
    real(ikind) :: z0
        
    z0=h/exp(0.4/sqrt(Cd)+1.0)
    z0=max(z0,1.0e-6) 
    
    return
    endfunction fric_conv_drag2length
    
!*****************************************************    
    function fric_conv_length2drag(z0,h) result(Cd)
! Converts a roughness height to a drag coefficient
!
! Usage:
!   Cd = fric_conv_length2drag(z0,h)
!
! Input:
!   z0 - roughness length [L]
!   h - Water depth [L]
!
! Output:
!   Cd - Drag coefficient [-]
!
! Author: Alex Sanchez, USACE-CHL
!*****************************************************
    use prec_def
    implicit none
    real(ikind),intent(in) :: z0,h
    real(ikind) :: z0lim,Cd
    
    z0lim = max(min(z0,0.5*h),1.0e-5)
    Cd=(0.4/(log(h/z0lim)-1.0))**2
    
    return
    endfunction fric_conv_length2drag
    
!********************************************************    
    function fric_conv_length2man(z0,h) result(cman)
! Converts a roughness height to a Manning's coefficient
!
! Usage:
!   fric_conv_length2man(z0,h)
!
! Input:
!   z0 - roughness length [L]
!   h - Water depth [L]
!
! Output:
!   cman - Manning's roughness coefficient [T/L^(1/3)]
!
! Author: Alex Sanchez, USACE-CHL
!********************************************************
    use prec_def
    implicit none
    real(ikind),intent(in) :: z0,h
    real(ikind) :: z0lim,cman
    
    z0lim = max(min(z0,0.5*h),1.0e-5)
    cman = h**0.166666667/(18.0*log(0.4*h/z0lim)) !0.4=12/30
    
    return
    endfunction fric_conv_length2man    
    
!*****************************************************************
    function fric_normapprough(uwc,umag,cb) result(zap)
! Calculates the normalize apparent roughness length zap = za/h
! Author: Alex Sanchez, USACE-CHL
!*****************************************************************    
    use prec_def
    implicit none
    real(ikind),intent(in) :: uwc,umag,cb
    real(ikind) :: cbwc,zap
    
    cbwc=uwc/max(umag,1.0e-6)*cb
    zap=exp(-1.0-0.4/sqrt(max(cbwc,1.0e-5))) 
    zap=min(zap,0.01) !Limit
    
    return
    endfunction fric_normapprough

!*****************************************************************    
    function fric_conv_nlength2pow(zap) result(cm)
! Converts the normalized roughness length zap to a power-law
! roughness parameter
!
! Description:
!   The power-law velocity distribution is given by
!     u(z) = (m+1)/m*(z/h)^(1/m)*U
!   where
!    z - height above bed [L]
!    m - roughness parameter typically between 6 and 7 [-]
!    u - current velocity at height z [L/T]
!    U - depth-averaged current velocity [L/T]
!
! Author: Alex Sanchez, USACE-CHL
!*****************************************************************        
    use prec_def
    implicit none
    real(ikind),intent(in) :: zap
    real(ikind) :: cm
    
    cm = -2.0*log10(zap)
    
    return
    endfunction fric_conv_nlength2pow    

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! End Bottom Roughness Conversions      
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
! Begin Wave Friction Factors 
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<   

!***********************************************************************    
    function fric_wavefac_soulsby(r) result(fwr)
! Wave friction factor for rough turbulent flow based on Soulsby (1997)
!
! Description:
!   Computes the wave friction factor based on Soulsby (1997)
!   for smooth laminar, smooth turbulent, and rough turbulent flow.
!   If the wave Reynolds number is not specified then the wave 
!   friction factor is for rough turbulent flow
!   For random waves the relative roughness and wave Reynolds number
!   should be estimated by taking Uw = sqrt(2)*Urms and Tw = Tp
!
! Input:
!   r = Aw/rks - Relative wave roughness [-]
!   where
!     Aw = Uw*Tw/(2*pi) - Amplitude of near-bed wave exursion [L]
!     rks = 30.0*z0 - Nikuradse roughness height [L]
!     z0 - roughness length [L]
!     Uw - near-bed wave orbital velocity [L/T]
!        (for random waves Uw = sqrt(2)*Urms)
!     Tw - wave period [T]
!        (for random waves Tw = Tp)
!
! Output:
!   fwr = Wave friction factor for rough turbuulent flow [-]
!
! References:
!   Soulsby, R.L. (1997). "Dynamics of marine sands",
!      Thomas Telford, 249 p.
!
! Author: Alex Sanchez, USACE-CHL
!***********************************************************************
    use const_def, only: twopi
    use prec_def
    implicit none
    !Input
    real(ikind),intent(in) :: r
    !Output
    real(ikind) :: fwr
    
    !Limited to 0.4432
    fwr = 0.237*max(r,0.3)**(-0.52) !rough turbulent flow friction factor = 1.39*(z0/A)**0.52
    
    return
    endfunction fric_wavefac_soulsby
    
!****************************************************************    
    function fric_wavefac_swart(r) result(fwr)
! Wave friction factor for rough turbulent flow  of Swart (1974)
!
! Input:
!   r = Aw/rks - Relative wave roughness [-]
!   where
!     Aw = Uw*Tw/(2*pi) - Amplitude of near-bed wave exursion [L]
!     rks = 30.0*z0 - Nikuradse roughness height [L]
!     z0 - roughness length [L]
!     Uw - near-bed wave orbital velocity [L/T]
!        (for random waves Uw = sqrt(2)*Urms)
!     Tw - wave period [T]
!        (for random waves Tw = Tp)
!
!  Output:
!   fwr = Wave friction factor for rough turbulent flow [-]
!
! References:
!   Swart, D.H. (1974). Offshore sediment transport 
!     and equilibrium beach profiles. Technical report. 
!     Delft, The Netherlands: Delft Hydraulics Laboratory.
!
! Author: Alex Sanchez, USACE-CHL
!****************************************************************
    use prec_def
    implicit none
    real(ikind),intent(in) :: r
    real(ikind) :: fwr
    
    fwr = exp(5.213*r**(-0.194) - 5.977)
    fwr = min(fwr,0.3)

    return
    endfunction fric_wavefac_swart
    
!********************************************************************
    function fric_wavefac_nielson(r) result(fwr)
! Wave friction factor of Nielson (1992)
!
! Description:
!   Calculates the wave frictionf actor fwr for rough turbulent flow
!   based on a recalculated form of Swart's (1974) equation by
!   Nielson (1992).
!   For random waves the relative roughness and wave Reynolds number
!   should be estimated by taking Uw = sqrt(2)*Urms and Tw = Tp
!
! Input:
!   r = Aw/rks - Relative wave roughness [-]
!   where
!     Aw = Uw*Tw/(2*pi) - Amplitude of near-bed wave exursion [L]
!     rks = 30.0*z0 - Nikuradse roughness height [L]
!     z0 - roughness length [L]
!     Uw - near-bed wave orbital velocity [L/T]
!        (for random waves Uw = sqrt(2)*Urms)
!     Tw - wave period [T]
!        (for random waves Tw = Tp)
!
! Output:
!   fwr - Wave friction factor for rough turbulent flow [-]
!
! References:
!   Swart, D.H. (1974). Offshore sediment transport
!     and equilibrium beach profiles. Technical report.
!     Delft, The Netherlands: Delft Hydraulics Laboratory.
!   Nielson, P. (1992). Coastal bottom boundary layers and
!     sediment transport. Advanced Series on Ocean Engineering,
!     Vol. 4. World Scientific. 324 p.
!
! written by Alex Sanchez, USACE-CHL
!***************************************************************
    use prec_def
    implicit none
    real(ikind),intent(in) :: r
    real(ikind) :: fwr
    
    fwr = exp(5.5*r**(-0.2) - 6.3)
    fwr = min(fwr,0.3)

    return
    endfunction fric_wavefac_nielson
    
!********************************************************************
    function fric_wav_fac_huynhthanh_temperville(r) result(fwr)
! Wave friction factor based on Huynh-Thanh and Temperville (1991)
!
! Description:
!   Calculates the wave frictionf actor fwr for rough turbulent flow
!   based on Huynh-Thanh and Temperville (1991)
!   For random waves the relative roughness and wave Reynolds number
!   should be estimated by taking Uw = sqrt(2)*Urms and Tw = Tp
!
! Input:
!   r = Aw/rks - Relative wave roughness [-]
!   where
!     Aw = Uw*Tw/(2*pi) - Amplitude of near-bed wave exursion [L]
!     rks = 30.0*z0 - Nikuradse roughness height [L]
!     z0 - roughness length [L]
!     Uw - near-bed wave orbital velocity [L/T]
!        (for random waves Uw = sqrt(2)*Urms)
!     Tw - wave period [T]
!        (for random waves Tw = Tp)
!
! Output:
!   fwr - Wave friction factor for rough turbulent flow [-]
!
! References:
!   Huynh-Thanh, S., and Temperville, A. 1991. A numerical model
!     of the rough turbulent boundary layer in combined wave and
!     current interaction. in Sand Transport in Rivers, Estuaries, 
!     and the Sea, ed. R.L. Soulsby and R. Bettess, 
!     Balkema, Rotterdam. pp. 93-100.
!
! written by Alex Sanchez, USACE-CHL
!***************************************************************
    use prec_def
    implicit none
    real(ikind),intent(in) :: r
    real(ikind) :: fwr
    
    fwr = 0.00278*exp(4.65*r**(-0.22))
    !fwr = min(fwr,0.3)

    return
    endfunction fric_wav_fac_huynhthanh_temperville
    
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! End Wave Friction Factors 
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
! Begin Beforms 
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
 
!***************************************************************************
    subroutine fric_dune_vanrijn(h,d50,Ts,dhgt,dlen)
! Calculates the bed dune height and length based on van Rijn (1984)
! 
! Usage:
!   call fric_dune_vanrijn(h,d50,taus,taucr,dhgt,dlen)
!
! Input:
!   h - Water depth [m]
!   d50 - Median grain size [m]
!   Ts = thetap/thetacr - 1.0
!     where
!       thetap - Skin bed shields parameter [-]
!       thetacr - Critical shields parameter [-]
!
! Output:
!   dhgt - Dune height [m]
!   dlen - Dune length [m]
!
! References:
!   van Rijn, L.C. 1984. Sediment transport, Part III: Bed forms and alluvial
!      roughness. Journal of Hydraulic Engineering, 110(12), 1733-1754.
!
! Author: Alex Sanchez, USACE-CHL
!**************************************************************************
    use prec_def
    implicit none
    real(ikind),intent(in) :: h,d50,Ts
    real(ikind),intent(out) :: dhgt,dlen

    !Dune length
    dlen = 7.3*h

    !Dune height
    dhgt = 0.11*h*(d50/h)**0.3*(1.0-exp(-0.5*Ts))*(25.0-Ts)
    dhgt = max(dhgt,0.0) !for approximately (Ts<0 .and. Ts>26)

    return
    endsubroutine fric_dune_vanrijn

!************************************************************************************
    subroutine fric_ripple_soulsbywhitehouse(g,s,d50,dstar,thetacr,&
       h,Uc,Uw,Tw,Dw,Tbhalf,riphgt,riplen,ripdir,dt)
! Current and wave ripple dimensions based on Soulbsy and Whitehouse (2005)
!
! Input:
!   g = Gravity [m/s^2]
!   s = specific gravity [m]
!   d50 = median grain size [m] 
!   dstar = nondimensional grain size [-]
!   thetacr = Critical shields parameter
!   h = Total water depth [m]
!   Uc = Eulerian depth-averaged current velocity [m/s]
!   riphgtc = current-related ripple height at previous time step (if dt is specified)
!   riplenc = current-related ripple length at previous time step (if dt is specified)
!   Uw = bottom orbital velocity based on linear wave theory [m/s]
!       (for random waves use H10 ~= 1.27*H13)
!   Tw = wave period [s]
!       (for random waves use Tw = Tp)
!   Dw = Wave direction (going to, clockwise from North) [rad]
!   riphgtw = wave-related ripple height at previous time step (if dt is specified)
!   riplenw = wave-related ripple length at previous time step (if dt is specified)
!   dt = time step (set to 0.0 or large number for steady case) [s]
!
! Output:
!   riphgt = Ripple height at new time step [m]
!   riplen = Ripple length at new time step [m]
!   ripdir = Ripple direction at new time step [rad]
!
! References:    
!   Soulsby, R.L. and Whitehouse, R.J.S. (2005). "Prediction of 
!     Ripple Properties in Shelf Seas; Mark 2 Predictor for 
!     Time Evolution". Report TR 154, HR Wallingford, UK.  
!
! written by Alex Sanchez, USACE-CHL
!***************************************************************************
    use const_def, only: twopi
    use prec_def
    implicit none
    !Input Variables
    real(ikind),intent(in) :: g,s,d50,dstar,thetacr,h,Uc,Uw,Tw,Dw,Tbhalf
    real(ikind),intent(inout) :: riphgt,riplen,ripdir
    real(ikind),intent(in),optional :: dt !If not specified than the equilibrium values are used
    !Internal Variables
    real(ikind) :: betahgt,betalen,betamin
    real(ikind) :: Aw,del,psi,s1gravd50,val,Te,fw,z0g,Cd,r,Tb
    real(ikind) :: thetawp,thetacp,thetawop,thetasfp,c,b
    real(ikind) :: riphgteq,ripleneq,ripdireq
    
    s1gravd50 = (s-1.0)*g*d50 !Internal variable
    
    Aw = Uw*Tw/twopi !amplitude of near-bed wave exursion
    z0g = d50/12.0   !grain-related
    
    !--- Wave skin friction Shields parameter -----
    r = Aw/(30.0*z0g)
    fw = fric_wavefac_swart(r) !wave friction factor of Swart (1974)
    thetawp = 0.5*fw*Uw*Uw/s1gravd50 !skin friction shields parameter due to waves
    
    !--- Current skin friction Shields parameter -----
    Cd = fric_conv_length2drag(z0g,h)
    thetacp = Cd*Uc*Uc/s1gravd50 !skin friction shields parameter due to currents
    
    !--- Check for no change ----------
    if(max(thetawp,thetacp)<=thetacr) return
    
    !Determine dominant bedform
    if(thetawp>thetacp)then !Wave-dominant
      !--- Wave Ripple Dimensions -------------------
      del = Aw/d50  
      !Note: the case where thetacp<=thetacr is treated above    
      if(del<6000.0)then
        ripleneq = Aw/(1.0+1.87e-3*del*(1.0-exp(-(2.0e-4*del)**1.5))) !Maximum wave ripple length 
        riphgteq = ripleneq*0.15*(1.0-exp(-(5000.0/del)**3.5))   !Maximum current ripple height
      else !Wash-out
        ripleneq = Aw/9.2063 !9.2063=(1.0+1.87e-3*6000.0*(1.0-exp(-(2.0e-4*6000.0)**1.5))) !Maximum wave ripple height 
        riphgteq = 0.0   !Maximum current ripple length  
      endif
      !ripdireq = ripdir
    else !Current-dominant
      !--- Current Ripple Dimensions -------------------
      ripleneq = d50*(500.0+1881.0*dstar**(-1.5))
      riphgteq = d50*202.0*dstar**(-0.554)        
      !Apply Current Wash-out critiera
      if(dstar>1.58)then
        val =  dstar**(-1.3)
        thetawop = 1.66*val !Wash-out criteria
        thetasfp = 2.26*val !Sheet flow criteria
      else
        thetawop = 0.916 !Wash-out criteria 
        thetasfp = 1.25  !Sheet flow criteria
      endif
      !Note: The ripple length is unaffected by the wash-out
      !Note: The case where thetacp<=thetacr is treated above
      if(thetacp>thetawop .and. thetacp<=thetasfp)then
        riphgteq = riphgteq*(thetasfp-thetacp)/(thetasfp-thetawop)
      elseif(thetacp>thetasfp)then
        riphgteq = 0.0 !Sheet flow conditions
      !else  
        !do nothing
      endif
    endif
    
    !--- If not dt specified then use equilibrium values ---------------
    if(.not.present(dt))then
      riplen = ripleneq
      riphgt = riphgteq
      ripdir = Dw
      return
    endif

    !--- Calculate New Ripple Dimensions ---------------
    if(thetawp>thetacp)then !Wave-dominant
      Te = Tw !Time scale related to wave period
      psi = Uw**2/s1gravd50 !Wave mobility parameter
      val = psi**1.07
      !Rate of change parameter (same for length and height in the case of waves)
      betalen = 2.996*val/(21700.0+val) 
      betahgt = betalen
    else  !Current-dominant
      Te = riphgt*riplen/sqrt(s1gravd50*d50*d50) !Adaption time for bedload [s]
      val = (thetacp - thetacr)**1.5
      betalen = 12.0*val/(2.5+val) !Rate of change parameter for current ripple length
      betahgt = 20.0*val/(2.5+val) !Rate of change parameter for current ripple height
    endif
    betamin = 1.6*Te/dt !Method only stable for beta/Te*dt>1.6 so limit for simplicity   
    betalen = max(betalen,betamin)
    betahgt = max(betahgt,betamin)
    if(Tbhalf>0.0)then !Include biodegradation
      Tb = 5193.7*Tbhalf !5193.7=3600.0/log(2)
      c = betahgt/Te  !a = c*xii
      b = c + 1.0/Tb
      riphgt = riphgt + (c/b*riphgteq - riphgt)*(1.0-exp(-b*dt))
    else !No biodegradation
      riphgt = riphgt + (riphgteq - riphgt)*(1.0-exp(-betahgt/Te*dt))  
    endif
    riplen = riplen + (ripleneq - riplen)*(1.0-exp(-betalen/Te*dt))
    ripdir = ripdir + (ripdireq - ripdir)*(1.0-exp(-betalen/Te*dt)) !Note: Uses length rate of change parameter
    
    return
    endsubroutine fric_ripple_soulsbywhitehouse
    
!************************************************************************************
    subroutine fric_ripple_current_soulsbywhitehouse(g,s,d50,dstar,thetacr,&   !Sediment
       h,Uc,riphgtc,riplenc,&  !Hydro, current-related ripples
       dt) !Time step 
! Current and wave ripple dimensions based on Soulbsy and Whitehouse (2005)
!
! Usage:
!   call fric_ripple_current_soulsbywhitehouse(s,d50,dstar,thetacr,h,Uc,riphgtc,riplenc) !Steady current-only
!   call fric_ripple_current_soulsbywhitehouse(s,d50,dstar,thetacr,h,Uc,riphgtc,riplenc,dt) !Transient current-only n,eed to use keyword argument
!
! Input:
!   dt = time step (set to 0.0 or large number for steady case) [s]
!   s = specific gravity [m]
!   d50 = median grain size [m] 
!   dstar = nondimensional grain size [-]
!   thetacr = Critical shields parameter
!   h = Total water depth [m]
!   Uc = Eulerian depth-averaged current velocity [m/s]
!   riphgtc = current-related ripple height at previous time step (if dt is specified)
!   riplenc = current-related ripple length at previous time step (if dt is specified)
!   dt = Time step
!
! Output:
!   riphgtc = current-related ripple height at new time step [m]
!   riplenc = current-related ripple length at new time step [m]
!
! References:    
!   Soulsby, R.L. and Whitehouse, R.J.S. (2005). "Prediction of 
!     Ripple Properties in Shelf Seas; Mark 2 Predictor for 
!     Time Evolution". Report TR 154, HR Wallingford, 
!     Wallingford, UK.  
!
! written by Alex Sanchez, USACE-CHL
!***************************************************************************
    use const_def, only: twopi
    use prec_def
    implicit none
    !Input Variables
    real(ikind),intent(in) :: g,s,d50,dstar,thetacr,h,Uc
    real(ikind),intent(inout) :: riphgtc,riplenc
    real(ikind),intent(in),optional :: dt !If not specified than the equilibrium values are used
    !Internal Variables
    real(ikind) :: betahgt,betalen,betahgtc,betalenc,betamin
    real(ikind) :: s1gravd50,val,Te,Tc,z0g,Cd  !del,psi
    real(ikind) :: thetacp,thetasfp,thetawop
    real(ikind) :: riphgtceq,riplenceq
    
    s1gravd50 = (s-1.0)*g*d50 !Internal variable
    
    z0g = d50/12.0   !grain-related
    
    !--- Current skin friction Shields parameter -----
    Cd = fric_conv_length2drag(z0g,h)
    thetacp = Cd*Uc*Uc/s1gravd50 !skin friction shields parameter due to currents
    
    !!--- Check for no change ----------
    !if(thetacp<=thetacr) return
    
    !--- Current Ripple Dimensions -------------------
    riplenceq = d50*(500.0+1881.0*dstar**(-1.5)) !Maximum current ripple length
    riphgtceq = d50*202.0*dstar**(-0.554)        !Maximum current ripple height 
    !Apply Current Wash-out critiera
    if(dstar>1.58)then
      val =  dstar**(-1.3)
      thetawop = 1.66*val !Wash-out criteria
      thetasfp = 2.26*val !Sheet flow criteria
    else
      thetawop = 0.916 !Wash-out criteria 
      thetasfp = 1.25  !Sheet flow criteria
    endif
    !Note: The ripple length is unaffected by the wash-out
    !Note: The case where thetacp<=thetacr is treated above
    if(thetacp>thetawop .and. thetacp<=thetasfp)then
      riphgtceq = riphgtceq*(thetasfp-thetacp)/(thetasfp-thetawop)
    elseif(thetacp>thetasfp)then
      riphgtceq = 0.0 !Sheet flow conditions
    !else  
    !  riphgtceq = riphgtc
    endif  
    
    !--- If not dt specified then use equilibrium values ---------------
    if(.not.present(dt))then
      riplenc = riplenceq
      riphgtc = riphgtceq
      return
    endif

    !--- Calculate New Current Ripple Dimensions ---------------
    Tc = riphgtc*riplenc/sqrt(s1gravd50*d50*d50) !Adaption time for bedload [s]
    val = (thetacp - thetacr)**1.5
    betalenc = 12.0*val/(2.5+val) 
    betahgtc = 20.0*val/(2.5+val) !Rate of change parameter for current ripple height
    betamin = 1.6*Tc/dt !Method only stable for beta/Te*dt>1.6 so limit for simplicity
    betahgt = max(betalenc,betamin)
    betalen = max(betahgtc,betamin)
    riphgtc = riphgtc + (riphgtceq - riphgtc)*(1.0-exp(-betahgt/Te*dt))
    riplenc = riplenc + (riplenceq - riplenc)*(1.0-exp(-betalen/Te*dt))     
      
    return
    endsubroutine fric_ripple_current_soulsbywhitehouse

!************************************************************************************
    subroutine fric_ripple_wave_soulsbywhitehouse(g,s,d50,dstar,thetacr,&   !Sediment
       Uw,Tw,riphgtw,riplenw,& !Waves, waves-related rippples
       dt) !Time step 
! Wave ripple dimensions based on Soulbsy and Whitehouse (2005)
!
! Usage:
!   call fric_ripple_wave_soulsbywhitehouse(s,d50,dstar,thetacr,h,Uc,riphgtc,riplenc) !Steady current-only
!   call fric_ripple_wave_soulsbywhitehouse(s,d50,dstar,thetacr,h,Uc,riphgtc,riplenc,dt=dt) !Transient current-only n,eed to use keyword argument
!   call fric_ripple_wave_soulsbywhitehouse(s,d50,dstar,thetacr,h,Uc,riphgtc,riplenc,Uw,Tw,riphgtw,riplenw) !Steady waves and currents
!   call fric_ripple_wave_soulsbywhitehouse(s,d50,dstar,thetacr,h,Uc,riphgtc,riplenc,Uw,Tw,riphgtw,riplenw,dt) !Steady waves and currents
!
! Input:
!   dt = time step (set to 0.0 or large number for steady case) [s]
!   s = specific gravity [m]
!   d50 = median grain size [m] 
!   dstar = nondimensional grain size [-]
!   thetacr = Critical shields parameter
!   Uw = bottom orbital velocity based on linear wave theory [m/s]
!       (for random waves use H10 ~= 1.27*H13)
!   Tw = wave period [s]
!       (for random waves use Tw = Tp)
!   riphgtw = wave-related ripple height at previous time step (if dt is specified)
!   riplenw = wave-related ripple length at previous time step (if dt is specified)
!
! Output:
!   riphgtw = current-related ripple height at new time step [m]
!   riplenw = current-related ripple length at new time step [m]
!
! References:    
!   Soulsby, R.L. and Whitehouse, R.J.S. (2005). "Prediction of 
!     Ripple Properties in Shelf Seas; Mark 2 Predictor for 
!     Time Evolution". Report TR 154, HR Wallingford, 
!     Wallingford, UK.  
!
! written by Alex Sanchez, USACE-CHL
!***************************************************************************
    use const_def, only: twopi
    use prec_def
    implicit none
    !Input Variables
    real(ikind),intent(in) :: g,s,d50,dstar,thetacr,Uw,Tw
    real(ikind),intent(inout) :: riphgtw,riplenw
    real(ikind),intent(in),optional :: dt !If not specified than the equilibrium values are used
    !Internal Variables
    real(ikind) :: betalen,betaw,betamin              !betahgt
    real(ikind) :: Aw,del,psi,s1gravd50,val,fw,z0g,r  !Te
    real(ikind) :: thetawp                            !Thetawop,Thetasfp
    real(ikind) :: riphgtweq,riplenweq

    real(ikind) :: rdum
    rdum=dstar       !This can probably be removed
    rdum=thetacr   !This can probably be removed
    s1gravd50 = (s-1.0)*g*d50 !Internal variable
    
    Aw = Uw*Tw/twopi !amplitude of near-bed wave exursion
    z0g = d50/12.0   !grain-related
    
    !--- Wave skin friction Shields parameter -----
    r = Aw/(30.0*z0g)
    fw = fric_wavefac_swart(r) !wave friction factor of Swart (1974)
    thetawp = 0.5*fw*Uw*Uw/s1gravd50 !skin friction shields parameter due to waves
    
    !!--- Check for no change ----------
    !if(thetawp<=thetacr) return
    
    !--- Wave Ripple Dimensions -------------------
    del = Aw/d50  
    !Note: the case where thetacp<=thetacr is treated above    
    if(del<6000.0)then
      riplenweq = Aw/(1.0+1.87e-3*del*(1.0-exp(-(2.0e-4*del)**1.5))) !Maximum wave ripple length 
      riphgtweq = riplenweq*0.15*(1.0-exp(-(5000.0/del)**3.5))   !Maximum current ripple height
    else !Wash-out
      riplenweq = Aw/9.2063 !9.2063=(1.0+1.87e-3*6000.0*(1.0-exp(-(2.0e-4*6000.0)**1.5))) !Maximum wave ripple height 
      riphgtweq = 0.0   !Maximum current ripple length  
    endif
        
    !--- If not dt specified then use equilibrium values ---------------
    if(.not.present(dt))then
      riplenw = riplenweq
      riphgtw = riphgtweq
      return
    endif

    !--- Calculate New Wave Ripple Dimensions ---------------
    psi = Uw**2/s1gravd50 !Wave mobility parameter
    val = psi**1.07
    betaw = 2.996*val/(21700.0+val) !Note: same for length and height    
    betamin = 1.6*Tw/dt !Method only stable for beta/Te*dt>1.6 so limit for simplicity
    betaw = max(betaw,betamin)
    betalen = max(betaw,betamin)
    riphgtw = riphgtw + (riphgtweq - riphgtw)*(1.0-exp(-betaw/Tw*dt))
    riplenw = riplenw + (riplenweq - riplenw)*(1.0-exp(-betaw/Tw*dt))   
      
    return
    endsubroutine fric_ripple_wave_soulsbywhitehouse
       
!*****************************************************************************    
    subroutine fric_ripple_wave_vanrijn(g,s,d50,Uw,Tw,riphgtw,riplenw)
!  Wave ripple dimensions based on van Rijn (1984,1989)
!
! Description:
!   Calculates the bed ripple dimensions following the same procedure as in 
!   the Lund-CIRP formulation
!
! Input:
!   g = gravity [m/s^2]
!   s = specific gravity [-]
!   d50 = median grain size [m] 
!   Uw = bottom orbital velocity based on linear wave theory [m/s]
!       (for random waves use H10 H10 ~= 1.27*H13)
!   Tw = wave period [s]
!       (for random waves use Tw = Tp)
!
! Output:
!   riphgt = ripple height at new time step [m]
!   riplen = ripple length at new time step [m]
!
! References:
!  van Rijn, L.C. 1984b. Sediment transport: Part II: Suspended load transport.
!    Journal of Hydraulic Engineering 110(11), 1,613-1,641.
!  van Rijn, L.C. 1984c. Sediment transport: Part III: Bed form sand alluvial
!    roughness. Journal of Hydraulic Engineering 110(12), 1,733-1,754.
! van Rijn, L.C. 1989. Handbook of sediment transport by currents and waves. 
!    Technical Report H461. The Hague: Delft Hydraulics Laboratory.
!
! written by Alex Sanchez, USACE-CHL
!*****************************************************************************    
    use const_def, only: twopi
    use prec_def
    implicit none
    !Input/Output
    real(ikind),intent(in) :: g,s,d50,Uw,Tw
    real(ikind),intent(out) :: riphgtw,riplenw
    !Internal variables
    real(ikind) :: Pw,Aw,s1gravd50

    !-- Wave ripples ----
    s1gravd50 = (s-1.0)*g*d50 !Internal variable
    Aw = Uw*Tw/twopi  !amplitude of near-bed wave exursion
    Pw = Uw**2/s1gravd50 !Wave mobility parameter
    if(Pw<10.0)then
      riphgtw = 0.22*Aw
      riplenw = 1.25*Aw
    elseif(Pw<250.0)then
      riphgtw = 2.8e-13*Aw*(250.0-Pw)**5.0
      riplenw = 1.4e-6*Aw*(250.0-Pw)**2.5
    else !Wash-out
      riphgtw = 0.0
      riplenw = 0.0
    endif

    return
    endsubroutine fric_ripple_wave_vanrijn
    
!*****************************************************************************
    subroutine fric_ripple_current_soulsby(d50,riplenc,riphgtc)
! Calculates the current-related rippple length and height based on Soulsby (1997)
! Author: Alex Sanchez, USACE-CHL
!*****************************************************************************
    use prec_def
    implicit none
    real(ikind),intent(in) :: d50
    real(ikind),intent(out) :: riplenc,riphgtc
    
    riplenc = 1000.0*d50
    riphgtc = riplenc/7.0
    
    return
    endsubroutine fric_ripple_current_soulsby
    
!***************************************************************************
    subroutine fric_ripple_current_raudkivi(d50,riphgtc,riplenc)
! Calculates the current-related ripple height and length based on Raudkivi (1998)
!
! Input:
!   d50 - Median grain size [m]
!
! Output:
!   riphgt - Current-related ripple height [m]
!   riplen - Current-related ripple length [m]
!
! References:
!   Raudkivi, K.A.H. 1998. Loose Boundary Hydraulics. A.A. Balkema,
!      Rotterdam, The Netherlands, 496 p.
!
! Author: Alex Sanchez, USACE-CHL    
!***************************************************************************
    use prec_def
    implicit none
    real(ikind), intent(in) :: d50
    real(ikind), intent(out) :: riphgtc,riplenc
    
    riplenc = 2.7489*d50**0.35 !1000**0.35 = 11.2202 [m]
    riphgtc = 0.0129*riplenc*d50**(-0.253) !1000**(-0.253) [m]
    
    !!Above equivalent to
    !riplenc = 0.245*(d50*1000.0)**0.35  !Raudikivi ripples
    !riphgtc = 0.074*riplenc/(d50*1000.0)**0.253
    
    return
    endsubroutine fric_ripple_current_raudkivi
  
!******************************************************************
    subroutine fric_megaripple(h,Uc,Uw,megariphgt,megariplen)
! Computes the megaripple height and length based on the data collected
! by Miles et al. (2013)
!
! Input:
!   h = Water depth [m]
!   Uc = Current velocity [m/s]
!   Uw = Peak bottom wave orbital velocity [m/s]
!
! Output:
!   megariphgt = Megaripple height [m]
!   megariplen = Megaripple length [m]
!
! References:
!   Miles, J.,  Thorpe, A., Russell, P., and  Masselink, G. 2013. 
!     Field observation of ripple and megaripple dynamics in the nearshore.
!     Coastal Dynamics 2013. 1207-1218. 
!
! Author: Alex Sanchez, USACE-CHL
!******************************************************************   
    use prec_def
    implicit none
    real(ikind),intent(in) :: h,Uc,Uw
    real(ikind) :: megariphgt,megariplen
    
    megariphgt = 0.3*h**0.6*(0.7-exp(-4*Uc))*(0.6-Uc)*(0.5-Uw)**0.1
    megariphgt = max(megariphgt,0.0)
    megariplen = max(min(1.0*h,5.0),0.3)
    
    return
    endsubroutine fric_megaripple
    
!!******************************************************************
!    subroutine ripple_C09(s,d50,Uw,Tw,Uc,riphgt,riplen)
!!  Ripple dimensions based on Camenen (2009)
!!
!! Description:
!!   Calculates the bed ripple dimensions following the same procedure as in 
!!   the Lund-CIRP formulation
!!
!! Input:
!!   s = specific gravity [m]
!!   d50 = median grain size [m] 
!!   Uw = bottom orbital velocity based on linear wave theory [m/s]
!!       (for random waves use H10 H10 ~= 1.27*H13)
!!   Tw = wave period [s]
!!       (for random waves use Tw = Tp)
!!
!! Output:
!!   riphgt = ripple height at new time step [m]
!!   riplen = ripple length at new time step [m]
!!
!! References:
!!   Camenen, B. 2009. Estimation of the wave-related ripple characteristics 
!!     and induced bed shear stress. Estuarine, Coastal, and Shelf Science,
!!     84, 553-564.
!!   Camenen, B., Bayram, A., and Larson, M., 2006. Equivalent roughness 
!!     height for plane bed under steady flow, Journal of Hydraulic 
!!     Engineering, 132, 11, 1146-1158.
!!
!! written by Alex Sanchez, USACE-CHL
!!*****************************************************************************    
!    use const_def, only: twopi,grav
!    use prec_def
!    implicit none
!    !Input/Output
!    real(ikind),intent(in) :: s,d50,Uw,Tw
!    real(ikind),intent(out) :: riphgt,riplen
!    !Internal variables
!    real(ikind) :: Pw,Aw,s1gravd50,thetasf,thetawg,delw,xt
!    real(ikind) :: riphgtw,riplenw,riphgtc,riplenc,fw
!    
!    xt = d50/((s-1.0)*grav)/Tw**2
!    s1gravd50 = (s-1.0)*grav*d50 !Internal variable
!    Aw = Uw*Tw/twopi  !amplitude of near-bed wave exursion
!    Pw = Uw**2/s1gravd50 !Wave mobility parameter
!    
!!    ksg = 2.5*d50
!    delw = sqrt(visc)*Tw/pi
!    fw = 0.5/sqrt(pi)*dstar**(-0.75)*d50/delw
!    thetawg = 0.5*fw*Pw
!    rwc = (1.0-0.6*tanh(50.0*Uc/Uw))*(1.0+3.0*tanh(Uc/Uw))
!    thetasf = 10.0*dstar**(-0.75)*sqrt(d50/delw)*(1.0+rwc)
!
!    fac = -0.025*xt**(-0.2)
!    ripleneq = Aw*1.6*exp(-5.0e-3*xt**(-0.4))*thetawg**fac
!    riphgteq = ripleneq*0.16*exp(-1.0e-7*(Aw/d50)**2)
!    
!    fac = exp(-acr*(thetacr/thetawg)**4)*exp(-acr*(thetawg/thetasf)**4)
!    riplen = ripleneq*fac
!    riphgt = riphgteq*fac
!    
!    return
!    endsubroutine ripple_C09
    
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! End Beforms
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
! Begin Roughness Estimates 
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

!******************************************************************
    function fric_rough_ripple(riphgt,riplen) result(rksr)
! Computes the bed roughness height due to ripples
!
! Input:
!   dhgt = Ripple height [L]
!   dlen = Ripple length [L]
!
! Output:
!   rsksd = Nikuradse roughness height due to ripples L]
!
! References:
!   Soulsby, R.L. (1997). "Dynamics of marine sands", 
!      Thomas Telford, 249 p.
!
! Author: Alex Sanchez, USACE-CHL
!******************************************************************
    use prec_def
    implicit none
    real(ikind),intent(in) :: riphgt,riplen
    real(ikind) :: rksr
    
    rksr = 8.01*riphgt*riphgt/max(riplen,riphgt,1.0e-4) !8.01 = 0.267*30 as recommended by Soulsby (1997)
    
    return
    endfunction fric_rough_ripple
    
!******************************************************************
    function fric_rough_dune(dhgt,dlen) result(rksd)
! Computes the bed roughness height due to dunes baesd on van Rijn (1984)
!
! Input:
!   dhgt = Dune height [L]
!   dlen = Dune length [L]
!
! Output:
!   rsksd = Nikuradse roughness height due to dunes [L]
!
! References:
!  van Rijn, L.C. 1984. Sediment transport: Part III: Bed form sand alluvial
!    roughness. Journal of Hydraulic Engineering 110(12), 1,733-1,754.
!
! Author: Alex Sanchez
!******************************************************************
    use prec_def
    implicit none
    real(ikind),intent(in) :: dhgt,dlen
    real(ikind) :: rksd
    
    rksd = 1.1*dhgt*(1.0-exp(-25.0*dhgt/dlen))
    
    return
    endfunction fric_rough_dune    
    
!***************************************************************************
    function fric_rough_grain(d90) result(rksg)
! Estimates the grain-related roughness height
! Although the routine is simple it is provided for completness
!
! Usage:
!   rksg = fric_rough_grain(d90)
!
! Input: 
!   d90 - 90th percentile sediment grain size [m]
!
! Output:
!   rksg - Grain-related roughness height [m]
!
! References:
!  van Rijn, L.C. 1984. Sediment transport: Part III: Bed form sand alluvial
!    roughness. Journal of Hydraulic Engineering 110(12), 1,733-1,754.
!
! written by Alex Sanchez, USACE-CHL
!****************************************************************************
    use prec_def
    implicit none
    !Input/Output
    real(ikind),intent(in) :: d90
    real(ikind) :: rksg
    
    rksg = 3.0*d90  !3.0 is based on van Rijn (1984)
    
    return
    endfunction fric_rough_grain
    
!*****************************************************************************
    function fric_rough_trans_nielson(d50,thetap,thetacr) result(rkst)
! Estimates the sediment transport related roughness length z0t
! based on Nielson (1992)
!
! Input:
!   d50 = Median grain size [m]
!   thetap = Nondimensional shear stress (Shields parameter) for 
!          based on the grain roughness [-]
!   thetacr = Critical nondiemsnional shear stress (from Sheilds diagram) [-]
!
! Output:
!   z0t = Roughness length due to sediment transport near the bed [m]
!
! Reference:
!  Nielson, P. (1992). Coastal bottom boundary layers and
!    sediment transport. Advanced Series on Ocean Engineering, 
!    Vol. 4. World Scientific. 324 p.
!
! Author: by Alex Sanchez, USACE-CHL
!*****************************************************************************    
    use prec_def
    implicit none
    !Input/Output
    real(ikind),intent(in) :: d50,thetap,thetacr
    real(ikind) :: rkst
    
    if(thetap>thetacr)then
      rkst = 170.0*sqrt(thetap-thetacr)*d50 !Eq. (3.6.10)
    else
      rkst = 0.0
    endif
    
    return
    endfunction fric_rough_trans_nielson    
    
!*********************************************************************************
    function fric_rough_trans_wiberg(d50,Ts) result(rkst)
! Estimates the sediment transport related roughness length z0t
! based on Wiberg and Rubin (1989)
!
! Input:
!   d50 = Median grain size [m]
!   Ts = thetap/thetacr - 1
!     where
!       thetap = Nondimensional shear stress (Shields parameter) for 
!          based on the grain roughness [-]
!       thetacr = Critical nondiemsnional shear stress (from Sheilds diagram) [-]
!
! Output:
!   z0t = Roughness length due to sediment transport near the bed [m]
!
! Reference:
!   Wiberg, P.L. and Rubin, D.M. 1989. Bed roughness produced by 
!     saltating sediment. Journal of Geophysical Research 94(C4), 5011-5016.
!
! Author: Alex Sanchez, USACE-CHL
!**********************************************************************************
    use prec_def
    implicit none
    real(ikind),intent(in) :: d50,Ts
    !real(ikind),parameter :: a = 0.056
    !real(ikind),parameter :: a1 = 0.068
    real(ikind) :: a2
    real(ikind) :: rkst
    
    a2 = 0.0204*log(100.0*D50**2) + 0.0709*log(100.0*d50)
    rkst = 0.0038*d50*(Ts/(1.0+a2*Ts))  !0.0038 = a*a1
    
    return
    endfunction fric_rough_trans_wiberg
    
!***************************************************************************
    function fric_rough_trans_wilsonlog(g,rhow,rhos,Uc,dep) result(rkstc)
! Sediment Transport roughness length using Wilson formula and 
!  a log velocity profile
!
! Description:
!   Calculates the sediment transport roughness length (z0t) by using a 
!   polynomial approximation to the solution of the Wilson formula assuming
!   a log profile for the current velocity over a flat bed. 
!   The polynomials are formulated in non-dimensional terms
!
! Input:
!   g - Gravity [m/s^2]
!   rhow - Water density [kg/m^3]
!   rhos - Sediment density [kg/m^3]
!   Uc - Depth-averaged current velocity [m/s]
!   dep - Water depth [m]
!
! Output:
!   rkstc - Roughness height due to current-related sediment transport [m]
!
! Original code from Magnus Larson
! Translated by Alex Sanchez
!***************************************************************************
    use prec_def
    implicit none
    !Input
    real(ikind),intent(in) :: g,rhow,rhos,Uc,dep
    !Output
    real(ikind) :: rkstc
    !Internal
    real(ikind) :: xi,psi
    
    if(Uc>0.001 .and. dep>=0.001)then
      xi = log(Uc**2/((rhos/rhow-1.0)*g*dep))
      if(xi<-4.0)then
        psi = -4.167248 + 1.269405*xi + 0.0083*xi**2
      elseif(xi>=-4.0 .and. xi<-1.0)then
        psi = -3.958030 + 1.376399*xi + 0.022333*xi**2
      else
        psi = -3.907731 + 1.461098*xi + 0.086438*xi**2 + 0.028461*xi**3
      endif
      rkstc = exp(psi)*dep
    else
      rkstc = 0.0
    endif
    
    return
    endfunction fric_rough_trans_wilsonlog

!***************************************************************************
    function fric_rough_trans_wilsonswart(g,rhow,rhos,Uw,T) result(rkstw)
! Sediment Transport roughness length using Wilson formula and the Swart
! wave friction factor
!
! Description:
!   Calculates the sediment transport roughness length (z0t) by using a 
!   polynomial approximation to the solution of the Wilson formula assuming
!   a log profile for the current velocity over a flat bed. 
!   The polynomials are formulated in non-dimensional terms
!
! Input:
!   g - gravity [m/s^2]
!   rhow - Water density [kg/m^3]
!   rhos - Sediment density [kg/m^3]
!   Uw - Bottom wave orbital velocity [m/s]
!   T - Wave period [s]
!
! Output:
!   rkstw - Roughness height due to wave-related sediment transport [m]
!
! Original code from Magnus Larson
! Translated by Alex Sanchez
!***************************************************************************
    use const_def, only: pi
    use prec_def
    implicit none
    !Input
    real(ikind),intent(in) :: g,rhow,rhos,Uw,T
    !Output
    real(ikind) :: rkstw
    !Internal
    real(ikind) :: xi,psi,Aw

    if(UW>0.0 .and. T>0.0)then
      xi = log(Uw/((rhos/rhow-1.0)*g*T))
      if(xi<-5.0)then
        psi = 0.58184 + 1.62956*xi + 0.03004*xi**2
      elseif(xi>=-5.0 .and. xi<-3.0)then
        psi = 3.23319 + 2.70983*xi + 0.14110*xi**2
      else
        psi = 51.59618 + 56.02861*xi + 19.85308*xi**2  &
            + 2.43875*xi**3
      endif
      Aw = Uw*T/(2.0*pi)
      rkstw = exp(psi)*Aw
    else
      rkstw = 0.0
    endif      
      
    return
    endfunction fric_rough_trans_wilsonswart
    
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! End Roughness Estimates
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
! Begin various
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
!*****************************************************************************
    subroutine fric_skin(tauwc,riphgt,riplen,z0g,z0t,tausfavg,tausfmax)
!*****************************************************************************    
    use prec_def
    implicit none
    !Input/Output
    real(ikind),intent(in) :: tauwc,riphgt,riplen,z0g,z0t
    real(ikind),intent(out) :: tausfavg,tausfmax
    !Internal
    real(ikind) :: fac
    real(ikind),parameter :: CdBF = 0.5
    
    fac = 1.0 + 3.125*CdBF*riphgt/riplen*(log(riphgt/(z0g+z0t))-1)**2   !3.125=0.5/0.4/0.4
    tausfavg = tauwc/fac
    tausfmax = tausfavg*(1.0+8.0*riphgt/riplen)
    
    return
    endsubroutine fric_skin
    
!******************************************************************************
    subroutine bed_current_wave_stress(taub,vx,vy,usx,usy,ux,uy,um,taubx,tauby)
! Separates the current bed shear stress into components
!
! Input:
!   taub - Combined wave current bed shear stress
!   vx,vy - Total flux velocities [m/s]
!   usx,usy - Wave flux velocities [m/s]
!   ux,uy - Depth-averaged current velocities [m/s]
!
! Output:
!   taubx,tauby - Bed shear stresses [Pa]
!
! written by Alex Sanchez, USACE-CHL
!******************************************************************************
    use prec_def
    implicit none
    !Input/Output
    real(ikind),intent(in) :: taub,vx,vy,usx,usy
    real(ikind),intent(out) :: ux,uy,um,taubx,tauby
    
    ux = vx-usx
    uy = vy-usy
    um = sqrt(ux*ux+uy*uy)
    um = max(um,1.0e-8)
    taubx = taub*ux/um
    tauby = taub*uy/um
    
    return
    endsubroutine bed_current_wave_stress        

!**************************************************************
    function fric_streaming_stress(rhow,z0,uw,Tw,wlen) result(taustr)
! Streamming Stress Magnitude
!
! Input:
!   rhow - Water density [M/L^3]
!   z0 - Roughness length [L]
!   uw - Wave bottom orbital velocity [L/T]
!   Tw - Wave period [T]
!   wlen - Wave length [L]
!
! Output:
!   taustr - Bottom streamming stress magnitude [M/L/T^2]
!
! written by Alex Sanchez, USACE-CHL
!**************************************************************
    use const_def, only: twopi
    use prec_def
    implicit none
    !Input/Output
    real(ikind),intent(in) :: rhow,z0,uw,Tw,wlen
    real(ikind) :: taustr
    real(ikind) :: rdum    !This can probably be removed
    !Internal
    real(ikind) :: fw,z0s  !phi,A,d50,wa,wk,re,fe
    
    !z0s = z0
    rdum=z0                   !This can probably be removed
    z0s = 0.15/1000.0/12.0 !=0.2/1000/12 ******** HARD CODED FOR NOW *********************
    
    !Wave friction factor
    !r = Aw/(30.0*z0s)
    !fw = fric_wavefac_swart(r)
    fw = 3.6147*(z0s/max(uw,1.0e-6)/Tw)**0.52   !Soulsby (1997), = 1.39*(Aw/z0)^(-0.52)
    
    !Internal variables
    !wa = twopi/Tw
    !wk = twopi/wlen
    !A = max(uw/wa,1.0e-5)
    
    !!First approach    
    !taustr = rhow*0.176776695*wk*A**3*wa**2*fw !1/4/sqrt(2)=0.176776695    
    !taustr = 0.1768*rhow*uw**3*fw/wlen*Tp
    
    !Alternate approach using wave bottom friction
    !fw = 3.6147*(z0s/max(uw,1.0e-6)/Tw)**0.52   !Soulsby (1997), = 1.39*(Aw/z0)^(-0.52)
    !c = wk/wa !wLen/Tw
    !Df = rhow/(2*sqrt(pi))*fw*uw^3
    !fD = 1.25
    !taustr = fD*Df/c !Used in Delft3D and Unibest
    !taustr = 0.3526*rhow*fw*uw**3*wlen/Tw !fD/(2*sqrt(pi), fD = 1.25 as recommended by Roelvink
    !taustr = 0.2821*rhow*fw*uw**3*wlen/Tw !1/(2*sqrt(pi)
    
    !Nielsen
    !d50 = 0.2/1000.0
    !fw = exp(5.5*(2.5*d50/A)**0.2-6.3)
    !phi = 0.5*fw*(A*wa)**2/(1.65*9.81*d50)
    !re = 170.0*d50*sqrt(max(phi,0.05001)-0.05)
    !fe = exp(5.5*(re/A)**0.2-6.3)
    taustr = 0.0637*rhow*fw*uw**3*wlen/Tw !(2*pi)/(10*pi^2)=0.0637
    
    return
    endfunction fric_streaming_stress
    
!**********************************************************    
    function fric_relroughwave(z0,Uw,Tw) result(r)
! Calculates the relative wave roughness
!
! Input:
!   z0 = roughness length [m]
!   Uw = near-bed wave orbital velocity [m/s]
!      (for random waves Uw = sqrt(2)*Urms)
!   Tw = wave period [s]
!      (for random waves Tw = Tp)
!
! Output:
!   r = relative wave roughness [m]
!
! Author: Alex Sanchez, USACE-CHL
!**********************************************************
    use const_def, only: twopi
    use prec_def
    implicit none
    !Input
    real(ikind),intent(in) :: z0,Uw,Tw
    !Output
    real(ikind) :: r
    !Internal
    real(ikind) :: Aw,rks
    
    Aw = Uw*Tw/twopi !amplitude of near-bed wave exursion
    rks = 30.0*z0    !Nikuradse roughness height
    r = Aw/rks       !relative roughness
    
    return
    endfunction fric_relroughwave
    
!********************************************************
    function fric_wavereynolds(mu,Uw,Tw) result(Rew)
! Computes the wave Reynolds number
!
! Input: 
!   mu = Kinematic viscosity [m^2/s]
!   Uw = near-bed wave orbital velocity [m/s]
!      (for random waves Uw = sqrt(2)*Urms)
!   Tw = wave period [s]
!      (for random waves Tw = Tp)
!
! Output:
!   Rew = Wave Reynolds number [-]
!
! Author: Alex Sanchez, USACE-CHL
!********************************************************
    use const_def, only: twopi
    use prec_def
    implicit none
    real(ikind),intent(in) :: mu,Uw,Tw
    real(ikind) :: Rew,Aw
    
    Aw = Uw*Tw/twopi
    Rew = Uw*Aw/max(mu,1.0e-20) !Wave Reynolds number [-]

    return
    endfunction fric_wavereynolds
    
endmodule fric_lib

