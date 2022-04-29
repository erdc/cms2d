!==========================================================================
module flow_lib
! Flow library
!
! Contains:
!   water_viscosity_kinematic - Calculates kinematic viscosity as a function
!             of temperature
!   water_viscosity_dynamic - Calculates the dynamic viscosity based on the 
!             data from Riley and Skirrow (1965). 
!   water_density_crowley - Crowley (1968) formula for water density as a
!                         function of water temperature and salinity
!   water_density_fofonoff - Fofonoff (1985) formula for water density as a
!                         function of water temperature and salinity
!   water_density_ekart - Ekart (1958) formula for water density as a
!                      function of water temperature and salinity
!   water_density_UNESCO - UNESCO (1981) formula for water density as a
!                      function of water temperature and salinity
!   streamwise_curvature - Computes the streamwise curvature of a
!                          2DH velocity field
!
! written by Alex Sanchez, USACE-CHL
!=========================================================================
    use prec_def
    implicit none
    
contains

!*****************************************************************************    
    function coriolis_fplane(latref) result(fc)
! Calculates the Coriolis parameter using the f-plane approximation
! written by Alex Sanchez, USACE-CHL
!*****************************************************************************    
    use const_def, only: deg2rad,omega
    implicit none
    real(ikind),intent(in) :: latref
    real(ikind) :: fc

    fc = 2.0*omega*sin(latref*deg2rad)
    
    return
    end function coriolis_fplane
    
!*****************************************************************************    
    subroutine coriolis_betaplane(ynor,yref,latref,fc)
! Calculates the Coriolis parameter using the beta-plane approximation
! The reference northing distance yref may be calculated as
!   yref = sum(y(1:ncells))/real(ncells,kind=ikind)
! written by Alex Sanchez, USACE-CHL
!*****************************************************************************    
!    use size_def, only: ncellsD                         !ncellsD was never used, so commented out  MEB  01/26/2022
    use const_def, only: deg2rad,omega,earthRadius
    use prec_def
    implicit none
    !Input/Output
    real(ikind),intent(in) :: latref,yref
    real(ikind),intent(in) :: ynor(:)
    real(ikind) :: fc(:)
    !Internal variables
    real(ikind) :: fref,betaref
    
    fref = coriolis_fplane(latref)
    betaref = 2.0*omega*cos(latref*deg2rad)/earthRadius    
    fc(:) = fref + betaref*(ynor(:)-yref)   !Beta-plane approximation
    
    return
    end subroutine coriolis_betaplane    

!*****************************************************************************    
    function water_viscosity_kinematic(t) result(vk)
! Calculates kinematic viscosity based on the re-evaluation of the viscosity
! data from Korson et al. (1969) by Kelvin et al. (1978). The validity range
! for the equation is from 0 to 40 ºC. Assumes salinity is zero.
! Salinity has the effect of increasing slightly the viscosity.
!
! Input:
!   t - Temperature [ºC]
!   s - Salinity [practical units ~ ppt] 
!
! Output:
!   vk - Kinematic viscosity of water [m^2/s]
!
! References:
!   Kelvin, J., Sokolov, M., Wakeham, W.A. 1978. Viscosity of liquid water in
!     the range -8 ºC to 150 ºC. Journal of Physical Chemistry, 7(3), 941-948.
!   Korson, L., Drost-Hansen, W., Millero, F.J. 1969. Journal of 
!     Physical Chem, 75, 2016.
!
! written by Alex Sanchez, USACE-CHL
!*****************************************************************************
    implicit none
    real(ikind),intent(in) :: t
    real(ikind) :: vk,t2,t3

    t2 = t*t
    t3 = t2*t
    !vk = (1.79 - 0.05873*t + 0.001166*t2 - 1.017e-05*t3)*1.0e-6 !third order fit        
    vk = (1.793 - 0.06149*t + 0.001519*t2 - 2.436e-5*t3 + 1.773e-7*t3*t)*1.0e-6 !forth order fit
    
    return
    end function water_viscosity_kinematic
    
!************************************************************
    function water_viscosity_dynamic(t,s) result(vd)
! Calculates the dynamic viscosity based on the data from
! Riley and Skirrow (1965). A curve fit was presented by 
! Neumeier et al. (2008). 
! The error compared to the data of Riley and Skirrow (1965)
! is less than 0.5% over the range (0–38 salinity and 8–24ºC).
! The error is less than 1.0% over the range 0–28ºC.
!
! Input:
!   t - Temperature [ºC]
!   s - Salinity [practical units ~ ppt] 
!
! Output:
!   vd - Dynamic viscosity of water [kg*s/m]
!
! References:
!   Riley, J.P., and Skirrow, G., 1965. Chemical Oceanography, 
!     Vol 3. Academic Press, London, 564 pp.
!   Neumeier, U., Ferrarin, C., Amos, C.L., Umgiesser, G., and
!     Li, M.Z. 2008. Sedtrans05: An improved sediment-transport
!     model for continental shelves and coastal waters with a 
!     algorithm for cohesive sediments. Computers & Geosciences, 
!     34, 1223-1242. 
!
! written by Alex Sanchez, USACE-CHL
!************************************************************
    implicit none
    real(ikind),intent(in) :: t,s
    real(ikind) :: vd
        
    vd = 1.802863e-3 - 6.1086e-5*t + 1.31419e-6*t**2 &
     - 1.35576e-8*t**3 + 2.15123e-6*s + 3.59406e-11*s**2
    
    return
    end function water_viscosity_dynamic

!*******************************************************************
    function water_density_crowley(t,s) result(rho)
! Crowley (1968) formula for water density as a function
! of water temperature and salinity without the pressure component.
!
! Input:
!   t - Water temperature [ºC]
!   s - Salinity [ppt]
!
! Output:
!   rho - Water density [kg/m^3]
!
! References:
!   Crowley, W.P. 1968. A global numerical ocean model: Part I.
!     Journal of Computational Physics. 3, 243 p.
!
! written by Alex Sanchez, USACE-CHL
!*******************************************************************        
    implicit none
    real(ikind),intent(in) :: t,s
    real(ikind) :: rho
    
    rho = 1028.14 - 0.0735*t - 0.00469*t*t + (0.802 - 0.002*t)*(s - 35.0)
    
    return
    end function water_density_crowley    
    
!*******************************************************************
    function water_density_fofonoff(t,s) result(rho)
! Fofonoff (1985) formula for water density as a function
! of water temperature and salinity without the pressure component.
!
! Input:
!   t - Water temperature [ºC]
!   s - Salinity [practical units ~ ppt]
!
! Output:
!   rho - Water density [kg/m^3]
!
! References:
!   Fofonoff, N.P., 1985. Physical properties of seawater: 
!     A new salinity scale and equation of state for seawater. 
!     Journal of Geophysical Research 90, 3332–3342.
!
! written by Alex Sanchez, USACE-CHL
!*******************************************************************
    implicit none
    real(ikind),intent(in) :: t,s
    real(ikind) :: rho
    !Internal
    real(ikind) :: t2,t3,t4,t5
    
    t2 = t*t
    t3 = t2*t
    t4 = t3*t
    t5 = t4*t
    rho = 999.842594 + 6.793952e-2*t - 9.095290e-3*t2 &
        + 1.00168e-4*t3 - 1.120083e-6*t4 + 6.536332e-9*t5 &
        + (8.24493e-1 - 4.0899e-3*t + 7.6438e-5*t2 &
        - 8.2467e-7*t3 + 5.3875e-9*t4)*s &
        + (-5.72466e-3 + 1.0227e-4*t - 1.6546e-6*t2)*s**1.5 &
        + 4.8314e-4*s**2
    
    return
    end function water_density_fofonoff
    
!*******************************************************
    function water_density_ekart(t,s) result(rho)
! Calculates the water density based on Ekart (1958)
! t - Water temperature in ºC
! s - Salinity in ppt
! Valid for 0<t<40ºC and 0<s<40 ppt
! written by Alex Sanchez, USACE-CHL
!*******************************************************
    implicit none
    !Input/Output
    real(ikind),intent(in) :: t,s
    real(ikind) :: rho
    !Internal Variables
    real(ikind) :: lambda,P0
    real(ikind), parameter :: a0 = 0.6980
    
    lambda = 1779.5 + 11.25*t - 0.0745*t**2 - (3.80 + 0.01*t)*s
    P0 = 5890.0 + 38.0*t - 0.375*t**2 + 3.0*s
    rho = 1000.0*P0/(lambda + a0*P0)
    
    return
    end function water_density_ekart
    
!*********************************************************    
    function water_density_UNESCO(t,s) result(rho)
! Calculates the water density based on UNESCO (1981)
! t - Water temperature in ºC
! s - Salinity in ppt
! Valid for 0<t<40ºC and 0.5<s<43 ppt
! written by Alex Sanchez, USACE-CHL
!*********************************************************
    implicit none
    !Input/Output
    real(ikind),intent(in) :: t,s
    real(ikind) :: rho
    !Internal Variables
    real(ikind) :: A,B,rho0,t2,t3,t4,t5
    real(ikind), parameter :: C = 4.8314e-4
    
    t2 = t*t; t3 = t2*t; t4 = t3*t; t5 = t4*t
    rho0 = 999.842594 + 6.793952e-2*t - 9.095290e-3*t2 &
        + 1.001685e-4*t3 - 1.120083e-6*t4 + 6.536332e-9*t5
    A = 8.24493e-1 - 4.0899e-3*t + 7.6438e-5*t2 &
        - 8.2467e-7*t3 + 5.3875e-9*t4
    B = -5.72466e-3 + 1.0227e-4*t - 1.6546e-6*t2
    rho = rho0 + A*s + B*s**1.5 + C*s**2
    
    return
    end function water_density_UNESCO
    
!********************************************************************************
    function streamwise_curvature(ux,vy,uvm,dudx,dudy,dvdx,dvdy) result(Kc)
! Computes the stream-wise curvature of the 2DH flow field 
! following Theisel (1995)
!
! References:
!   Niell, S.P., Hashemi, M.R., Elliott, A.J. 2007. An enhanced depth-averaged
!      tidal model for morphological studies in the presence of rotary currents. 
!      Continental Shelf Research, 82-102.
!   Theisel, H., 1995. Vector field curvature and applications. Ph.D.
!      Thesis, Rostock University, Germany.
!
! written by Alex Sanchez, USACE-CHL
!********************************************************************************
    implicit none    
    real(ikind),intent(in):: ux,vy,uvm,dudx,dudy,dvdx,dvdy
    real(ikind):: Kc
    
    Kc = (ux*ux*dvdx-vy*vy*dudy+ux*vy*(dvdy-dudx))/max(uvm,1.e-6)**1.5
    
    return
    end function streamwise_curvature
    
!*******************************************************
    integer function number_wet_cells()
! Calculate and return the number of wet cells at the present timestep
! Written by: Mitchell Brown - 10/31/2016
!*******************************************************    
    use size_def, only: ncellsD
    use flow_def, only: iwet
    implicit none
    
    integer :: tmp, i
    
    tmp=0
    do i=1,ncellsD
      if (iwet(i)==1) then
        tmp = tmp+1
      endif
    enddo
    
    number_wet_cells = tmp
    
    end function number_wet_cells
    


end module flow_lib