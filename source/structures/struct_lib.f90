!===============================================================================
module struct_lib
! Structures Library
!
!   Rubble Mounds
!    The Forchheimer equation given by:
!       I = a*u + b*u**2
!      where I is the hydraulic gradient, u is the bulk velocity, 
!      and a and b are the dimensional coefficients.
!
! written by Weiming Wu, Chris Reed, Honghai Li, and Alex Sanchez
!===============================================================================    
    implicit none
    
    private
    public :: rubble_sidiropoulou,rubble_kadlec
    
contains    

!**************************************************************
    subroutine rubble_sidiropoulou(d,p,a,b)
! Description:    
! Calculates the dimensional coefficients a and b for the 
! Forchheimer equation used to simulate porous rubble mound
! structures according to Sidiropoulou et al. (2007). 
!
! Input:
!   d - Rock or riprap diameter [m]
!   p - Structure porosity [0]
!
! Output:
!   a,b - Coefficients
!
! References:
!  Sidiropoulou, M.G., K.N. Moutsopoulos, and 
!    V.A. Tsihrintzis. 2007. Determination of Forchheimer 
!    equation coefficients a and b. Hydrological Processes, 
!    21(4), 534–554.
!***************************************************************
    use prec_def
    implicit none
    !Input/Output
    real(ikind),intent(in) :: d,p
    real(ikind),intent(out) :: a,b
    
    a = 0.00333*d**(-1.5)*p**0.06         
    b = 0.194*d**(-1.265)*p**(-1.14)
                
    return
    end subroutine rubble_sidiropoulou
    
!**************************************************************
    subroutine rubble_kadlec(mu,g,d,p,a,b)
! Description:    
! Calculates the dimensional coefficients a and b for the 
! Forchheimer equation used to simulate porous rubble mound
! structures  according to Kadlec and Knight (1996)
!
! Input:
!   mu - water kinematic viscosity [m^2/s]
!   g - Gravity [m/s^2]
!   d - Rock or riprap diameter [m]
!   p - Structure porosity  [-]
!
! Output:
!   a,b - Coefficients
!
! References:
!  Kadlec, H.R., and L.R. Knight. 1996. Treatment Wetlands. 
!    Lewis Publishers.
!***************************************************************
    use prec_def
    implicit none
    !Input/Output
    real(ikind),intent(in) :: mu,g,d,p
    real(ikind),intent(out) :: a,b
    
    a = 255.0*mu*(1.0-p)/(g*p**3.7*d**2)
    b = 2.0*(1.0-p)/(g*p**3*d)
                
    return
    end subroutine rubble_kadlec
    
end module struct_lib