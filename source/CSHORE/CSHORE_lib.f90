!======================================================================
module CSHORE_lib
! CSHORE Library
! Contains only the modular and portable CSHORE routines
!==================================================================================
#include "CMS_cpp.h"
    implicit none

#ifdef CSHORE
contains
!********************************************************************
    subroutine CSHORE_point(J,Hrms,WT,H,WL,d50,WF,QRAW,Qonsh,Qoffsh)
! Calculates the net, onshore, and offshore sediment transport rate
!
! Input:
!  j - Cell ID
!  Hrms - Root-mean-squared wave height [m]
!  WT - wave period [s]
!  H - water depth [m]
!  WL - Wave length [m]
!  d50 - grain size [m]
!  WF - sediment fall velocity [m/s]
!
! Output:
!  Qraw - Net trnsport [m^2/s]
!  Qonsh - Onshore directed transport [m^2/s]
!  Qoffsh- offshore directed trans [m^2/s]
!
! Author: Brad Johnson, USACE-CHL
!
! Modified by:
!  Chris Reed, URS - Changes for CMS implementation
!  Alex Sanchez, USACE-CHL - Changes for CMS code updates
!
! Issue:
!  This routine should be rewritten to not use module variables as much
!  as possible and to output transport rates in x and y directions.
!
!********************************************************************
    use const_def, only: twopi
    use CSHORE_def
    use flow_def, only: grav
    use prec_def
    implicit none
    real(ikind),intent(in) :: J,Hrms,WT,H,WL,d50,WF
    real(ikind),intent(out) :: QRAW,Qonsh,Qoffsh
    integer:: i
    real(ikind):: KP,WHRMS,CPt,SIGMA,SIGSTA,SIGT,US,IGT,VSIGT,USTD,UMEAN
    real(ikind):: RB,RS,PB,PS,VS,QBX,QSX,ER,FCC,WFSGM1,GSD50S,WKP
    real(ikind):: Per,GBX,GF,DF,STA,USTA,SIGMA2,USIGT,DFSTA,DBSTA
    real(ikind):: ERFCC
    
    WKP = ((twopi/WT)**2)/grav
    WFSGM1 = WF*SGM1
    GSD50S = GSGM1*D50*SHIELD
    WHRMS = Hrms
    CPt = WL/WT
    SIGMA = HRMS/SQR8
    SIGSTA = SIGMA/H
    SIGT = SIGSTA*CPt
    SIGMA2 = SIGMA**2
    USIGT = -SIGSTA*grav*H/CPt/CPt
    USIGT = USIGT*(1.0_ikind+(CPt/grav)*RQ(J)/SIGMA2)
    VSIGT = 0.0_ikind
    call GBXAGF(CTHETA,USIGT,STHETA,VSIGT,GBX,GF)
    DFSTA = FB2*GF*SIGT**3/grav !Wave energy dissipation due to bottom friction
    DBSTA = RBETA(J)*RQ(J)      !Wave energy dissipation due to breaking
    USTD = SIGT*CTHETA
    UMEAN= -USTD*SIGSTA*grav*H/CPt/CPt
    USTA = UMEAN/SIGT
    RB = sqrt(GSD50S/FB2)/USTD
    RS = WF/USTD/FB2**0.3333_ikind
    US = USTA
    !Probabilities of sediment movement
    PB = 0.5_ikind*(ERFCC((RB-US)/SQR2)+ERFCC((RB+US)/SQR2)) !BUG FIX: RBUS not defined at this point. Replaced with Equation (49) (Alex 10/23/14)
    PS = 0.5_ikind*(ERFCC((RS-US)/SQR2)+ERFCC((RS+US)/SQR2)) !BUG FIX: RSUS not defined at this point. Replaced with Equation (50)  (Alex 10/23/14)
    if(PS>PB) PS = PB
    !Suspended sediment volume per unit horizontal bottom area
    VS = PS*(EFFF*DFSTA + EFFB*DBSTA)/WFSGM1
    VS = VS*sqrt(1.0_ikind+BSLOPE(J)*BSLOPE(J))
    QBX = BLD*PB*GSLOPE(J)*USTD**3
    QSX = ASLOPE(J)*UMEAN*VS
    Qonsh = QBX
    Qoffsh = QSX
    QRAW = (QBX + QSX) !/SPORO1
    
    return
    endsubroutine CSHORE_point

!***************************************************************
    function ERFCC(x)
! Complementary error function
!
! Author: Brad Johnson, USACE-CHL
! Modifications:
!  Alex Sanchez, USACE-CHL - Modified for arbritrary precision
!****************************************************************
    use prec_def, only: ikind
    implicit none
    real(ikind) :: X, Z, T, ERFCC
    
    Z = abs(X)
    T = 1.0_ikind/(1.0_ikind+0.5_ikind*Z)
    ERFCC = T*exp(-Z*Z-1.26551223_ikind+T*(1.00002368_ikind+T*(.37409196_ikind+&
      T*(.09678418_ikind+T*(-.18628806_ikind+T*(.27886807_ikind+&
      T*(-1.13520398_ikind+T*(1.48851587_ikind+&
      T*(-.8221522_ikind+T*.17087277_ikind)))))))))
    if(X<0.0_ikind) ERFCC = 2.0_ikind - ERFCC
    
    return
    endfunction ERFCC
    
!**************************************************************
    subroutine GBXAGF(CTHETA,USIGT,STHETA,VSIGT,GBX,GF)
! Computes GBX and GF for specified CTHETA, USIGT,
! STHETA and VSIGT for Gaussian variable R
!
! Author: Brad Johnson, USACE-CHL
! Modifications:
!  Alex Sanchez, USACE-CHL - Modified for arbritrary precision
!***************************************************************
    use prec_def, only: ikind
    use CSHORE_def, only: sqr2,sqrg1
    implicit none
    real(ikind),intent(in) :: CTHETA,USIGT,STHETA,VSIGT
    real(ikind),intent(out) :: GBX,GF
    real(ikind) :: C1,C2,C3,H,HM,B,WKP,WT,ERFCC
    
    ! For normally incident waves, use analytical
    ! expresions involving complementary error function ERFCC below
    C1 = 1.0_ikind-ERFCC(USIGT/SQR2)
    C2 = SQRG1*exp(-USIGT*USIGT/2.0_ikind)
    C3 = 1.0_ikind + USIGT*USIGT
    GBX = C3*C1 + C2*USIGT
    GF = USIGT*(C3 + 2.0_ikind)*C1 + (C3 + 1.0_ikind)*C2
    
    return
    endsubroutine GBXAGF
    
#endif
endmodule CSHORE_lib