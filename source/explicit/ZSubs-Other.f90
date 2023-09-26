!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!
! BELOW FOLLOWS A THREE SUBROUTINES TO COMPUTE SEDIMENT TRANSPORT
! ACCORDING TO THE LUND FORMULA:
!
! 1. TOTAL SHEAR STRESSES UNDER WAVES AND CURRENT (SHEARLUND)
! 2. SUSPENDED LOAD TRANSPORT, INCLUDING REFERENCE CONCENTRATION
!    AND MIXING COEFFICIENT TO BE USED IN THE ADVECTION-DIFFUSION
!    MODEL (SUSPLUND)
! 3. BEAD LOAD TRANSPORT RATE (BEDLUND)
!
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!
! CODED BY M. LARSON
! 2005-11-01
!
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!
!
! *****************************************************************
!
! CALCULATION OF TOTAL SHEAR STRESS TO BE USED WITH THE LUND FORMULA
!
! ROUGHNESS ESTIMATED WITH REGARD TO EFFECTS OF GRAIN SIZE, BED FORMS,
! AND SEDIMENT TRANSPORT
!
! A VALUE ON IRIPPLE DIFFERENT FROM 1 WILL DISABLE RIPPLE ROUGHNESS
! CALCULATION
!
! CODED BY M. LARSON
! 2005-11-01, updated 11/19/2007
! 
! Presently this routine is not used - commenting out MEB 12/11/2013
!******************************************************************
!      SUBROUTINE BSHEARLUND(IRIPPLE,DEP,UC,UW,T,PHI,RHOS,RHOW,D50, &
!                           TAUCT,TAUWT,TAUWMT,TAUCWT,TAUCWMT,      &
!                           FCF,FWF,FCWF)
!      use prec_def
!      use const_def, only: pi, grav
!
!!cwr      INTEGER, PARAMETER :: IKIND=8
!      
!      IMPLICIT REAL(IKIND) (A-H,O-Z)
!
!      REAL(IKIND) RHOS,RHOW,KAPPA, D50
!
!      ! IRIPPLE=0 : DO NOT INCLUDE RIPPLES
!      ! IRIPPLE=1 : INCLUDE RIPPLES
!
!      ! CONSTANTS
!      !PI=4.d0*DATAN(1.D0)
!      !GRAV=9.81D0
!      KAPPA=0.4D0
!
!      ! BOTTOM EXCURSION AMPLITUDE
!      IF(UW.GT.0.D0)THEN 
!         AW=UW*T/(2.d0*PI)
!      ELSE
!         AW=0.d0
!      ENDIF
!
!      IF(IRIPPLE.EQ.1)THEN
!        ! BED FORMS AND ROUGHNESS
!        ! ***********************
!        ! WAVES
!        ! -----
!        ! VAN RIJN'S EQUATION FOR RIPPLE PROPERTIES
!        !   (UNDER IRREGULAR WAVES)
!        ! ESTIMATE VALUE ON NON-DIMENSIONAL PARAMETER PSI
!        IF(UW.GT.0.D0)THEN
!          PW=UW**2.D0/((RHOS/RHOW-1.d0)*GRAV*D50)
!          IF(PW.LE.10.D0)THEN
!            RIPWH=0.22d0*AW
!            RIPWL=1.25d0*AW
!            RHRWL=RIPWH/RIPWL
!          ELSEIF(PW.LE.250.D0)THEN
!            RIPWH=2.8D-13 *(250.d0-PW)**5.D0*AW
!            RIPWL=1.4D-6 *(250.d0-PW)**2.5D0*AW
!            RHRWL=RIPWH/RIPWL
!          ELSE
!            RIPWH=0.D0
!            RIPWL=0.D0
!            RHRWL=0.D0
!          ENDIF
!
!          ! ROUGHNESS FROM WAVE RIPPLES (RKWF)
!          RKWF=30.d0*0.25d0*RIPWH*RHRWL
!        ELSE
!          RKWF=0.d0
!        ENDIF
!
!        ! CURRENT
!        ! -------
!        IF(UC.GT.0.D0)THEN
!          RIPCL=1000.d0*D50    !11/19/2007
!          RIPCH=RIPCL/7.d0     !11/19/2007
!          RHRCL=RIPCH/RIPCL    !11/19/2007
!          RKCF=30.d0*0.25d0*RIPCH*RHRCL ! ROUGHNESS FROM CURRENT RIPPLES (RKCF)
!        ELSE
!          RKCF=0.d0
!        ENDIF
!      ELSE
!        RKCF=0.D0
!        RKWF=0.D0
!      ENDIF
!
!      ! ROUGHNESS RELATED TO GRAIN SIZE AND SEDIMENT TRANSPORT
!      ! ******************************************************
!      ! ROUGHNESS RELATED TO GRAIN FRICTION (RKS)
!      RKS=2.d0*D50
!
!      ! ROUGHNESS RELATED TO SEDIMENT TRANSPORT
!      !   ROUGHNESS DETERMINED BY SOLVING THE WILSON AND SWART FORMULA
!      !   SIMULTANEOUSLY - NUMERICALLY DETERMINED SOLUTION IN NON-DIMENSIONAL
!      !   TERMS APPROXIMATED BY POLYNOMIALS (FLAT BED)
!
!      ! WAVES (RKTW)
!      IF(UW.GT.0.0_ikind .and. t .gt. 0.0_ikind)THEN
!        XI=LOG(UW/(RHOS/RHOW-1.0_ikind)/GRAV/T)
!        IF(XI.LT.-5.0_ikind)THEN
!          PSI=0.58184D0+1.62956D0*XI+0.03004D0*XI**2
!        ELSEIF(XI.GE.-5.0_ikind.AND.XI.LT.-3.0_ikind)THEN
!          PSI=3.23319D0+2.70983D0*XI+0.14110D0*XI**2
!        ELSE
!          PSI=51.59618D0+56.02861D0*XI+19.85308D0*XI**2 &
!               +2.43875D0*XI**3
!        ENDIF
!        RKTW=EXP(PSI)*AW
!      ELSE
!        RKTW=0.0_ikind
!      ENDIF
!
!      ! CURRENT (RKTC)
!      !   ROUGHNESS DETERMINED BY SOLVING THE WILSON FORMULA ASSUMING A LOG
!      !   PROFILE FOR THE CURRENT VELOCITY - NUMERICALLY DETERMINED SOLUTION 
!      !   IN NON-DIMENSIONAL TERMS APPROXIMATED BY POLYNOMIALS (FLAT BED)
!      IF(UC.GT.0.001D0 .and. dep .ge. 0.001d0)THEN
!        XI=LOG(UC**2/(RHOS/RHOW-1.0_ikind)/GRAV/DEP)
!        IF(XI.LT.-4.0_ikind)THEN
!          PSI=-4.167248D0+1.269405D0*XI+0.0083D0*XI**2
!        ELSEIF(XI.GE.-4.D0.AND.XI.LT.-1.D0)THEN
!          PSI=-3.958030D0+1.376399D0*XI+0.022333D0*XI**2
!        ELSE
!          PSI=-3.907731D0+1.461098D0*XI+0.086438D0*XI**2+0.028461D0*XI**3
!        ENDIF
!        RKTC=EXP(PSI)*DEP
!      ELSE
!        RKTC=0.0_ikind
!      ENDIF
!
!      ! ROUGHNESS HEIGHTS
!      Z0S=RKS/30.0_ikind
!      Z0FW=RKWF/30.0_ikind
!      Z0FC=RKCF/30.0_ikind
!      Z0TW=RKTW/30.0_ikind
!      Z0TC=RKTC/30.0_ikind
!      Z0W=Z0S+Z0FW+Z0TW
!      Z0C=Z0S+Z0FC+Z0TC
!
!      ! FRICTION FACTORS
!      ! ****************
!      ! CURRENT FRICTION FACTOR (LOG PROFILE)
!      IF(DEP.GT.0.0_ikind)THEN
!        FCF=2.0_ikind*(KAPPA/(1.0_ikind+LOG(Z0C/DEP)))**2.D0
!        FCF = min(FCF,0.05D0)
!      ELSE
!        FCF=0.0_ikind
!      ENDIF
! 
!      ! WAVE FRICTION FACTOR (SWART)
!      IF(AW.GT.0.0_ikind)THEN
!        RRR=AW/(30.0_ikind*Z0W)
!        FWF=0.3D0
!        IF (RRR.GE.1.57) THEN
!          FWF=DMIN1(DEXP(-5.98D0+5.2D0*RRR**dble(-0.194)),0.3D0)
!        ENDIF
!      ELSE
!        FWF=0.0_ikind
!      ENDIF
!
!      ! FRICTION WEIGHTING FUNCTION
!      IF(ABS(UC).GT.0.0_ikind)THEN
!        X=ABS(UC)/(ABS(UC)+UW)
!      ELSE
!        X=0.0_ikind
!      END IF
!
!      ! WEIGHTED FRICTION VALUE (WAVE AND CURRENT COMBINED)
!      FCWF=X*FCF+(1.0_ikind-X)*FWF
! 
!      ! SHEAR STRESSES
!      ! **************
!      ! CURRENT
!      TAUCT=0.5d0*RHOW*FCF*UC**2
!
!      ! WAVES, MEAN VALUE
!      TAUWT=0.25d0*RHOW*FWF*UW**2
!
!      ! WAVES, MAXIMUM VALUE
!      TAUWMT=0.5d0*RHOW*FWF*UW**2
!
!      ! COMBINED, MEAN VALUE
!      TAUCWT=SQRT(TAUCT**2+TAUWT**2+2.0_ikind*TAUCT*TAUWT*ABS(COS(PHI)))
!
!      ! COMBINED, MAXIMUM VALUE
!      TAUCWMT=SQRT(TAUCT**2+TAUWMT**2+2.0_ikind*TAUCT*TAUWMT*ABS(COS(PHI)))
!
!      RETURN
!      END SUBROUTINE !BSHEARLUND
      
! *****************************************************************
!     LUND-CIRP FORMULA FOR PREDICTING SUSPENDED LOAD INCLUDING 
!     REFERENCE CONCENTRATION AND SEDIMENT MIXING COEFFICIENT
!
!     CODED BY M. LARSON
!     2005-11-01, updated 11/19/2007
!
! Presently this routine is not used - commenting out  MEB  12/11/2013
!******************************************************************
!      SUBROUTINE BSUSPLUND(DEP,UC,UW,RHOS,RHOW,D50,                   &
!                           WS,DST,TAUCT,TAUWT,TAUWMT,TAUCWT,TAUCWMT,  &
!                           FCF,FWF,DB,TAUCR,CRCW,EPSCW,QSS)
!      use prec_def
!      use const_def, only: pi, grav
!!cwr     INTEGER, PARAMETER :: IKIND=8
!
!      IMPLICIT REAL(IKIND) (A-H,O-Z)
!
!      REAL(IKIND) RHOS,RHOW,KAPPA,KC,KW,KB, D50
!
!      ! CONSTANTS
!      !PI=4.0_ikind*ATAN(1.0_ikind)
!      !GRAV=9.81_ikind
!      KAPPA=0.4_ikind
!
!      ! SHEAR VELOCITIES
!      ! ****************
!      ! CURRENT SHEAR VELOCITY
!      USTC=SQRT(0.5_ikind*FCF*UC**2)
!
!      ! WAVE SHEAR VELOCITY
!      USTW=SQRT(0.5_ikind*FWF*UW**2)
!
!      ! SHEAR VELOCITY AT INCIPIENT MOTION
!      USTCR=SQRT(TAUCR/RHOW)
!
!      ! SEDIMENT MIXING COEFFICIENT
!      ! ***************************
!      ! CURRENT ONLY
!      ! IF(USTC.GT.0.0_ikind)THEN   !REPLACED 5/31/07
!      IF(USTC.GT.0.0_ikind .and. WS.GT.0.01_ikind) THEN
!        IF(WS/USTC.LE.1.0_ikind)THEN
!          SIGMAC=0.4_ikind+3.5_ikind*(SIN(PI/2.0_ikind*WS/USTC))**2
!        ELSE
!          SIGMAC=1.0_ikind+2.9_ikind*(SIN(PI/2.0_ikind*USTC/WS))**2
!        ENDIF
!      ELSE
!        SIGMAC=0.0_ikind
!      ENDIF
!      KC=SIGMAC*KAPPA/6.0_ikind
!
!      ! WAVES ONLY
!      ! IF(USTW.GT.0.0_ikind)THEN   !REPLACED 5/31/07
!      IF(USTW.GT.0.0_ikind .and. WS.GT.0.01D0) THEN 
!        IF(WS/USTW.LE.1.0_ikind)THEN
!          SIGMAW=0.15_ikind+1.50_ikind*(SIN(PI/2.0_ikind*WS/USTW))**2
!        ELSE
!          SIGMAW=1.0_ikind+0.65_ikind*(SIN(PI/2.0_ikind*USTW/WS))**2
!        ENDIF
!      ELSE
!        SIGMAW=0.0_ikind
!      ENDIF
!      KW=SIGMAW*KAPPA/9.42_ikind
!
!      ! WAVES AND CURRENT COMBINED
!      ! UC is positive (calculated as sqrt(u^2 + v^2), so no need to take abs(uc)
!      IF(UC .GT. 0.0_ikind .AND. UW .GT. 0.0_ikind)THEN
!         Y=TAUCT/(TAUCT+TAUWMT)
!         SIGMAWC=Y*SIGMAC+(1.0_ikind-Y)*SIGMAW
!         KC=SIGMAWC*KAPPA/6.0_ikind
!         KW=SIGMAWC*KAPPA/9.42_ikind
!      ENDIF
!
!      ! WAVE BREAKING COEFFICIENT
!      KB=0.010_ikind
!
!      ! ENERGY DISSIPATION RATE
!      ! ***********************
!      ! CURRENT AND NON-BREAKING WAVES
!      DC=TAUCT*USTC
!      DW=TAUWMT*USTW
!
!      ! TOTAL WEIGHTED ENERGY DISSIPATION
!      DTOT=KB**3*DB+KC**3*DC+KW**3*DW
!
!      ! SEDIMENT MIXING COEFFICIENT
!      EPSCW=(DTOT/RHOW)**(1.0_ikind/3.0_ikind)*DEP
!
!      ! AVOID ZERO MIXING
!      EPSCW=MAX(EPSCW,1.0D-6)
!
!      ! REFERENCE CONCENTRATION
!      ! ***********************
!      THETACR  = TAUCR/((RHOS-RHOW)*GRAV*D50)   ! CRITICAL SHIELDS NUMBER FOR INCIPIENT MOTION
!      THETAC   = TAUCT/((RHOS-RHOW)*GRAV*D50)   ! CURRENT SHIELDS NUMBER
!      THETAW   = TAUWT/((RHOS-RHOW)*GRAV*D50)   ! WAVE SHIELDS NUMBER (MEAN)
!      THETAWM  = TAUWMT/((RHOS-RHOW)*GRAV*D50)  ! WAVE SHIELDS NUMBER (MAXIMUM)
!      THETACW  = TAUCWT/((RHOS-RHOW)*GRAV*D50)  ! WAVE-CURRENT SHIELDS NUMBER (MEAN)
!      THETACWM = TAUCWMT/((RHOS-RHOW)*GRAV*D50) ! WAVE-CURRENT SHIELDS NUMBER (MAXIMUM)
!      ACRCW    = 0.0035D0*EXP(-0.3d0*DST)      ! COEFFICIENT IN REFERENCE CONCENTRATION EQUATION
!
!      ! REFERENCE CONCENTRATION
!      IF(THETACWM .GT. 0.2_ikind*THETACR)THEN
!        CRCW=ACRCW*THETACW*EXP(-4.5d0*THETACR/THETACWM)
!      ELSE
!        CRCW=0.0_ikind
!      ENDIF
!
!      ! CALCULATE SUSPENDED TRANSPORT RATE
!      if (dep .gt. 0.0_ikind .and. epscw .gt. 0.0_ikind .and. ws .gt. 0.0_ikind) then
!        QSS=UC*CRCW*EPSCW/WS*(1.0_ikind-EXP(-WS*DEP/EPSCW))
!      else
!        QSS = 0.0_ikind
!      end if
!
!      RETURN
!      END SUBROUTINE !BSUSPLUND

! *****************************************************************
!     LUND-CIRP FORMULA FOR PREDICTING BED LOAD TRANSPORT 
!
!     CODED BY M. LARSON
!     2005-11-01
!
! Presently this routine is not used - commenting out  MEB 12/11/2013
!******************************************************************
!      SUBROUTINE BBEDLUND(RHOS,RHOW,TAUCT,TAUCWT,TAUCWMT,TAUCR,QBS)
!      use prec_def
!      use const_def, only: GRAV
!!cwr      INTEGER, PARAMETER :: IKIND=8
!      IMPLICIT REAL(IKIND) (A-H,O-Z)
!      REAL(IKIND) RHOS,RHOW
!!
!! CONSTANTS
!!
!      !GRAV=9.81_ikind
!!
!! TRANSPORT COEFFICIENTS
!!
!      A=12.0_ikind
!      B=4.5_ikind
!!
!! BED LOAD TRANSPORT
!!
!      if (taucwmt .gt. 0.2_ikind*taucr) then
!         QBS=A*SQRT(TAUCT/RHOW)*TAUCWT*EXP(-B*TAUCR/TAUCWMT)
!         QBS=QBS/(RHOS-RHOW)/GRAV
!      else
!         QBS = 0.0_ikind
!      end if
!
!      RETURN
!      END SUBROUTINE !BBEDLUND

!*******************************************************************************
! ROUTINE TO COMPUTE SEDIMENT TRANSPORT
!*******************************************************************************
      SUBROUTINE LUNDCIRP_(H,PERIOD,DB,U0,DEP,D50,PHI_M,VNU,RHOW,RHOS,QBS,QSS,IRIPPLE, &
          BDpart,CRCW,WSFALL,USTC,USTW,TAUCR,TAUCWTB,TAUCWMTB,EPSCW,FCF,FWF)
!*******************************************************************************
! INPUT VARIABLES TO LUND-CIRP SEDIMENT TRANSPORT CALCULATIONS
! ------------------------------------------------------------
! I = COUNTER (?; NOT USED)
! H = WAVE HEIGHT
! PERIOD = WAVE PERIOD
! DB = WAVE ENERGY DISSIPATION DUE TO BREAKING
! U0 = CURRENT SPEED (MAGNITUDE)
! DEP = WATER DEPTH
! D50 = MEDIAN GRAIN SIZE
! PHI_M = ANGLE BETWEEN WAVE AND CURRENT DIRECTION
! VNU = KINEMATIC VISCOSITY
! RHOW = DENSITY OF WATER
! RHOS = DENSITY OF SEDIMENT
! IRIPPLE = SWITCH FOR RIPPLE CALCULATIONS (=1 IMPLIES RIPPLES)
! SCALEBED = CALIBRATION FACTOR FOR BED LOAD (=1.0 IMPLIES ORIGINAL FORMULA)
! SCALESUS = CALIBRATION FACTOR FOR SUSPENDED LOAD (=1.0 IMPLIES ORIGINAL FORMULA)
!
! OUTPUT VARIABLES TO LUND-CIRP SEDIMENT TRANSPORT CALCULATIONS
! ------------------------------------------------------------
! QTOT = TOTAL TRANSPORT
! QBS = BED LOAD TRANSPORT
! QSS = SUSPENDED LOAD TRANSPORT
      USE const_def, only: PI
      USE flow_def,  only: GRAV

      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER IRIPPLE
      REAL*8 H,PERIOD,DB,U0,DEP,D50,L
      REAL*8 CG,CN,L0,PHI_M, VNU !,R
      REAL*8 RHOW, RHOS, QBS, QSS !,S, QTOT
      REAL*8 DVISC, S1, DSTAR, WSFALL !,SCALEBED, SCALESUS
      REAL*8 UC, DEPTH, UW, T, PHI
      REAL*8 TAUCR, THETACR
!ADDED TO DEFINE TYPES - MEB 12/11/2013
      REAL*8 TAUCT,TAUWT,TAUWMT,TAUCWT,TAUCWMT,TAUCTB,TAUCWTB,TAUCWMTB
      REAL*8 FCF,FWF,FCWF,CRCW,EPSCW,BDPART,USTC,USTW

!     INITIALIZATION 
      UW = 0.D0
      L0 = 0.D0
      L = 0.D0
      QBS = 0.d0
      QSS = 0.D0
      UC = U0
      DEPTH = DEP
      T = PERIOD
      PHI = PHI_M

! CONSTANTS
      !GRAV=9.81D0
      !PI=4.D0*DATAN(1.D0)
! HORIZONTAL BOTTOM ORBITAL VELOCITY AMPLITUDE
      IF(PERIOD.GT.0)THEN
         CALL DISPD(PERIOD,DEP,L,CG,CN)
         UW=0.5d0*H*GRAV*T/L/COSH(2.d0*(PI*DEP/L))
      ELSE
         UW=0.D0
      ENDIF
      DVISC = VNU
! SPECIFIC GRAVITY
      S1=(RHOS-RHOW)/RHOW
! DIMENSIONLESS GRAIN SIZE
      DSTAR=D50*(S1*GRAV/DVISC**2.D0)**(1.D0/3.D0)
! FALL SPEED (SOULSBY)
      WSFALL=DVISC/D50*((10.36D0**2.D0+1.049D0*DSTAR**3.D0)**0.5D0-10.36D0)
! CRITICAL SHIELDS PARAMETER FOR INCIPIENT MOTION (SOULSBY AND WHITEHOUSE)
      THETACR=0.30D0/(1.D0+1.2D0*DSTAR)+0.055D0*(1.D0-EXP(-0.02D0*DSTAR))
! CRITICAL SHEAR STRESS
      TAUCR=THETACR*(RHOS-RHOW)*GRAV*D50
      CALL ZSHEARLUND(IRIPPLE,DEPTH,UC,UW,T,PHI,RHOS,RHOW,D50, &
                      TAUCT,TAUWT,TAUWMT,TAUCWT,TAUCWMT,       &
                      TAUCTB,TAUCWTB,TAUCWMTB,                 &
                      FCF,FWF,FCWF)                            
      CALL ZSUSPLUND(DEPTH,UC,UW,RHOS,RHOW,D50,WSFALL,DSTAR,   &
                     TAUCT,TAUWT,TAUWMT,TAUCWT,TAUCWMT,        &
                     FCF,FWF,DB,TAUCR,CRCW,EPSCW,QSS,          &
                     BDpart,USTC,USTW)
      CALL ZBEDLUND(RHOS,RHOW,TAUCTB,TAUCWTB,TAUCWMTB,TAUCR,QBS)

      !QSS = SCALESUS*QSS
      !QBS = SCALEBED*QBS
      !QTOT = QSS + QBS
      RETURN
      END SUBROUTINE !LUNDCIRP_

!**********************************************************
!
! ROUTINE FOR SOLVING THE DISPERSION RELATIONSHIP USING
! A PADE APPROXIMATION
!
!*******************************************************************************
      SUBROUTINE DISPD(T,DD,L,CG,CN)
!*******************************************************************************
      USE const_def, only: PI,TWOPI
      USE flow_def, only: GRAV
      
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 T,DD,L,CG,CN,Y,F
      REAL*8 C  !Added to define explicit type - MEB 12/11/2013
! USE A PADE APPROXIMATION TO CALCULATE THE WAVELENGTH
      !GRAV=9.81D0
      !PI=4.D0*DATAN(1.D0)
      !TWOPI=2.D0*PI
      Y=(TWOPI/T)**2.D0*DD/GRAV

      F=Y+1.D0/(1.D0+Y*(0.66667D0+Y*(0.3555D0+Y*(0.16084D0+Y*(0.0632D0+  &
        Y*(0.02174D0+Y*(0.00654D0+Y*(0.00171D0+Y*(0.00039D0+Y*0.000111D0)))))))))

      L=TWOPI/SQRT(Y*F/DD**2.D0)
      C=4.D0*PI*DD/L
      IF(C.LT.15.D0) THEN
         CN=0.5D0*(1.D0+C/SINH(C))
      ELSE
         CN=0.5D0
      ENDIF
      CG=CN*L/T
      
      RETURN
      END SUBROUTINE !DISPD

! *****************************************************************
!
! CALCULATION OF TOTAL SHEAR STRESS TO BE USED WITH THE LUND FORMULA
!
! ROUGHNESS ESTIMATED WITH REGARD TO EFFECTS OF GRAIN SIZE, BED FORMS,
! AND SEDIMENT TRANSPORT
!
! A VALUE ON IRIPPLE DIFFERENT FROM 1 WILL DISABLE RIPPLE ROUGHNESS
! CALCULATION
!
! CODED BY M. LARSON
! 2005-11-01
!
!******************************************************************
      SUBROUTINE ZSHEARLUND(IRIPPLE,DEP,UC,UW,T,PHI,RHOS,RHOW,D50,  &
                            TAUCT,TAUWT,TAUWMT,TAUCWT,TAUCWMT,     &
                            TAUCTB,TAUCWTB,TAUCWMTB,FCF,FWF,FCWF)
      USE const_def, only: pi
      USE flow_def, only: grav
! VARIABLE DECLARATIONS
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 RHOS,RHOW,KAPPA,D50
! EXPLICIT DEFINITION OF VARIABLES - MEB 12/11/2013
      REAL*8   UW,AW,T,PW,RIPWH,RIPWL,RHRWL,UC,RIPCL,RIPCH,RHRCL,RKCF,RKWF
      REAL*8   RKS,XI,PSI,RKTW,DEP,RKTC,Z0S,Z0FW,Z0TW,Z0TC,Z0W,Z0C,Z0CB,Z0FC
      REAL*8   Z0WB,FCF,FCFB,RRR,RRRB,FWF,FWFB,X,FCWF,FCWFB
      REAL*8   TAUCT,TAUCTB,TAUWT,TAUWTB,TAUWMT,TAUWMTB,TAUCWT,PHI
      REAL*8   TAUCWTB,TAUCWMT,TAUCWMTB
      INTEGER  IRIPPLE

! IRIPPLE=0 : DO NOT INCLUDE RIPPLES
! IRIPPLE=1 : INCLUDE RIPPLES

! CONSTANTS
      !PI=4.d0*ATAN(1.D0)
      !GRAV=9.81D0
      KAPPA=0.4D0

! BOTTOM EXCURSION AMPLITUDE
      IF(UW.GT.0.D0)THEN 
         AW=UW*T/(2.d0*PI)
      ELSE
         AW=0.d0
      ENDIF
      IF(IRIPPLE.EQ.1)THEN

! BED FORMS AND ROUGHNESS
! ***********************
! WAVES
! -----
! VAN RIJN'S EQUATION FOR RIPPLE PROPERTIES
! (UNDER IRREGULAR WAVES)
!
! ESTIMATE VALUE ON NON-DIMENSIONAL PARAMETER PSI
         IF(UW.GT.0.D0)THEN
            PW=UW**2.D0/((RHOS/RHOW-1.d0)*GRAV*D50)
            IF(PW.LE.10.D0)THEN
               RIPWH=0.22d0*AW
               RIPWL=1.25d0*AW
               RHRWL=RIPWH/RIPWL
            ELSEIF(PW.LE.250.D0)THEN
               RIPWH=2.8D-13 *(250.d0-PW)**5.D0*AW
               RIPWL=1.4D-6 *(250.d0-PW)**2.5D0*AW
               RHRWL=RIPWH/RIPWL
            ELSE
               RIPWH=0.D0
               RIPWL=0.D0
               RHRWL=0.D0
            ENDIF

! ROUGHNESS FROM WAVE RIPPLES (RKWF)
            RKWF=30.d0*0.25d0*RIPWH*RHRWL

         ELSE
            RKWF=0.d0
         ENDIF

! CURRENT
! -------
         IF(UC.GT.0.D0)THEN
            RIPCL=1000.d0*D50
            RIPCH=RIPCL/7.d0
            RHRCL=RIPCH/RIPCL

! ROUGHNESS FROM CURRENT RIPPLES (RKCF)
            RKCF=30.d0*0.25d0*RIPCH*RHRCL
         ELSE
            RKCF=0.d0
         ENDIF
      ELSE
         RKCF=0.D0
         RKWF=0.D0
      ENDIF

! ROUGHNESS RELATED TO GRAIN SIZE AND SEDIMENT TRANSPORT
! ******************************************************
! ROUGHNESS RELATED TO GRAIN FRICTION (RKS)
      RKS=2.d0*D50

! ROUGHNESS RELATED TO SEDIMENT TRANSPORT
!
! ROUGHNESS DETERMINED BY SOLVING THE WILSON AND SWART FORMULA
! SIMULTANEOUSLY - NUMERICALLY DETERMINED SOLUTION IN NON-DIMENSIONAL
! TERMS APPROXIMATED BY POLYNOMIALS (FLAT BED)
!

! WAVES (RKTW)
      IF(UW.GT.0.D0 .and. t .gt. 0.d0)THEN
         XI=DLOG(UW/(RHOS/RHOW-1.d0)/GRAV/T)
         IF(XI.LT.-5.D0)THEN
            PSI=0.58184D0+1.62956D0*XI+0.03004D0*XI**2.D0
         ELSEIF(XI.GE.-5.D0.AND.XI.LT.-3.D0)THEN
            PSI=3.23319D0+2.70983D0*XI+0.14110D0*XI**2.D0
         ELSE
            PSI=51.59618D0+56.02861D0*XI+19.85308D0*XI**2.D0  &
                +2.43875D0*XI**3.D0
         ENDIF
         RKTW=DEXP(PSI)*AW
      ELSE
         RKTW=0.d0
      ENDIF

! CURRENT (RKTC)
!
! ROUGHNESS DETERMINED BY SOLVING THE WILSON FORMULA ASSUMING A LOG
! PROFILE FOR THE CURRENT VELOCITY - NUMERICALLY DETERMINED SOLUTION 
! IN NON-DIMENSIONAL TERMS APPROXIMATED BY POLYNOMIALS (FLAT BED)
      IF(UC.GT.0.001D0 .and. dep .ge. 0.001d0)THEN
         XI=DLOG(UC**2.D0/(RHOS/RHOW-1.d0)/GRAV/DEP)
         IF(XI.LT.-4.D0)THEN
            PSI=-4.167248D0+1.269405D0*XI+0.0083D0*XI**2.D0
         ELSEIF(XI.GE.-4.D0.AND.XI.LT.-1.D0)THEN
            PSI=-3.958030D0+1.376399D0*XI+0.022333D0*XI**2.D0
         ELSE
            PSI=-3.907731D0+1.461098D0*XI+0.086438D0*XI**2.D0+0.028461D0*XI**3.D0
         ENDIF
         RKTC=DEXP(PSI)*DEP
      ELSE
         RKTC=0.d0
      ENDIF

! ROUGHNESS HEIGHTS
      Z0S=RKS/30.d0
      Z0FW=RKWF/30.d0
      Z0FC=RKCF/30.d0
      Z0TW=RKTW/30.d0
      Z0TC=RKTC/30.d0
      Z0W=Z0S+Z0FW+Z0TW
      Z0C=Z0S+Z0FC+Z0TC

! ROUGHNESS FOR BED LOAD CALCULATIONS
      Z0CB=Z0S+Z0TC
      Z0WB=Z0S+Z0TW

! FRICTION FACTORS
! ****************
! CURRENT FRICTION FACTOR (LOG PROFILE)
      IF(DEP.GT.0.D0)THEN
         FCF=min(2.D0*(KAPPA/(1.d0+LOG(Z0C/DEP)))**2,0.3D0)  !0.05 cutoff added by C. Reed
         FCFB=min(2.D0*(KAPPA/(1.d0+LOG(Z0CB/DEP)))**2,0.3D0)  !0.05 cutoff added by C. Reed
      ELSE
         FCF=0.d0
         FCFB=0.d0
      ENDIF
! WAVE FRICTION FACTOR (SWART)
      IF(AW.GT.0.d0)THEN
         RRR=max(AW/(30.d0*Z0W),1.0d0)
         RRRB=max(AW/(30.d0*Z0WB),1.0d0)
         FWF=min(EXP(-5.98D0+5.2D0*RRR**(-0.194D0)),0.3D0)
         FWFB=min(EXP(-5.98D0+5.2D0*RRRB**(-0.194D0)),0.3D0)
      ELSE
         FWF=0.D0
         FWFB=0.D0
      ENDIF
! FRICTION WEIGHTING FUNCTION
      IF(ABS(UC).GT.0.D0)THEN
         X=ABS(UC)/(ABS(UC)+UW)
      ELSE
         X=0.D0
      END IF
! WEIGHTED FRICTION VALUE (WAVE AND CURRENT COMBINED)
      FCWF=X*FCF+(1.d0-X)*FWF
      FCWFB=X*FCFB+(1.d0-X)*FWFB

! SHEAR STRESSES
! **************
! CURRENT
      TAUCT=0.5d0*RHOW*FCF*UC**2.D0
      TAUCTB=0.5d0*RHOW*FCFB*UC**2.D0
! WAVES, MEAN VALUE
      TAUWT=0.25d0*RHOW*FWF*UW**2.D0
      TAUWTB=0.25d0*RHOW*FWFB*UW**2.D0
! WAVES, MAXIMUM VALUE
      TAUWMT=0.5d0*RHOW*FWF*UW**2.D0
      TAUWMTB=0.5d0*RHOW*FWFB*UW**2.D0
! COMBINED, MEAN VALUE
      TAUCWT=SQRT(TAUCT**2.D0+TAUWT**2.D0+2.d0*TAUCT*TAUWT*ABS(COS(PHI)))
      TAUCWTB=SQRT(TAUCTB**2.D0+TAUWTB**2.D0+2.d0*TAUCTB*TAUWTB*ABS(COS(PHI)))
! COMBINED, MAXIMUM VALUE
      TAUCWMT=SQRT(TAUCT**2.D0+TAUWMT**2.D0+2.d0*TAUCT*TAUWMT*ABS(COS(PHI)))
      TAUCWMTB=SQRT(TAUCTB**2.D0+TAUWMTB**2.D0+2.d0*TAUCTB*TAUWMTB*ABS(COS(PHI)))

      RETURN
      END SUBROUTINE !SHEARLUND

! *****************************************************************
!
!     LUND-CIRP FORMULA FOR PREDICTING SUSPENDED LOAD INCLUDING 
!     REFERENCE CONCENTRATION AND SEDIMENT MIXING COEFFICIENT
!
!     CODED BY M. LARSON
!     2005-11-01
!
!******************************************************************

      SUBROUTINE ZSUSPLUND(DEP,UC,UW,RHOS,RHOW,D50,                  &
                           WS,DST,TAUCT,TAUWT,TAUWMT,TAUCWT,TAUCWMT, &
                           FCF,FWF,DB,TAUCR,CRCW,EPSCW,QSS,          &
                           BDpart,USTC,USTW)
      USE const_def, only: PI
      USE flow_def, only: GRAV
! VARIABLE DECLARATIONS
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 RHOS,RHOW,KAPPA,KC,KW,KB,D50
! EXPLICIT DEFINITIONS - MEB 12/11/2013
      REAL*8 USTC,FCF,UC,USTW,FWF,UW,USTCR,TAUCR,WS,SIGMAC,SIGMAW,TAUCT,TAUWMT
      REAL*8 Y,SIGMAWC,DC,DW,DTOT,DB,EPSCW,DEP,THETACR,THETAC,THETAW,TAUWT
      REAL*8 THETAWM,THETACW,TAUCWMT,ACRCW,DST,CRCW,BDPART,QSS,TAUCWT,THETACWM

! CONSTANTS
      !PI=4.d0*DATAN(1.D0)
      !GRAV=9.81D0
      KAPPA=0.4D0

! SHEAR VELOCITIES
! ****************
! CURRENT SHEAR VELOCITY
      USTC=SQRT(0.5d0*FCF*UC**2.D0)
! WAVE SHEAR VELOCITY
      USTW=SQRT(0.5d0*FWF*UW**2.D0)
! SHEAR VELOCITY AT INCIPIENT MOTION
      USTCR=SQRT(TAUCR/RHOW)

! SEDIMENT MIXING COEFFICIENT
! ***************************
! CURRENT ONLY
      IF(USTC.GT.0.D0)THEN
         IF(WS/USTC.LE.1.D0)THEN
            SIGMAC=0.4d0+3.5d0*(SIN(PI/2.d0*WS/USTC))**2
         ELSE
            SIGMAC=1.d0+2.9d0*(SIN(PI/2.d0*USTC/WS))**2
         ENDIF
      ELSE
         SIGMAC=0.d0
      ENDIF
      KC=SIGMAC*KAPPA/6.d0
! WAVES ONLY
      IF(USTW.GT.0.D0)THEN
         IF(WS/USTW.LE.1.D0)THEN
            SIGMAW=0.15d0+1.50d0*(SIN(PI/2.d0*WS/USTW))**2
         ELSE
            SIGMAW=1.d0+0.65d0*(SIN(PI/2.d0*USTW/WS))**2
         ENDIF
      ELSE
         SIGMAW=0.d0
      ENDIF
      KW=SIGMAW*KAPPA/9.42d0
! WAVES AND CURRENT COMBINED
!      IF(TAUCT.GT.0.d0 .OR. TAUWMT.GT.0.d0)THEN
      IF(TAUCT.GT.0.d0 .AND. TAUWMT.GT.0.d0)THEN     !reported by Alex on 12/05/2008
         Y=TAUCT/(TAUCT+TAUWMT)
         SIGMAWC=Y*SIGMAC+(1.d0-Y)*SIGMAW
         KC=SIGMAWC*KAPPA/6.d0
         KW=SIGMAWC*KAPPA/9.42d0
      ENDIF
! WAVE BREAKING COEFFICIENT
      KB=0.010d0

! ENERGY DISSIPATION RATE
! ***********************
! CURRENT AND NON-BREAKING WAVES
      DC=TAUCT*USTC
      DW=TAUWMT*USTW
! TOTAL WEIGHTED ENERGY DISSIPATION
      DTOT=KB**3*DB+KC**3*DC+KW**3*DW
! SEDIMENT MIXING COEFFICIENT
      EPSCW=(DTOT/RHOW)**(1.D0/3.D0)*DEP
! AVOID ZERO MIXING
      EPSCW=MAX1(EPSCW,1.0D-6)

! REFERENCE CONCENTRATION
! ***********************
! CRITICAL SHIELDS NUMBER FOR INCIPIENT MOTION
      THETACR=TAUCR/((RHOS-RHOW)*GRAV*D50)
! CURRENT SHIELDS NUMBER
      THETAC=TAUCT/((RHOS-RHOW)*GRAV*D50)
! WAVE SHIELDS NUMBER (MEAN)
      THETAW=TAUWT/((RHOS-RHOW)*GRAV*D50)
! WAVE SHIELDS NUMBER (MAXIMUM)
      THETAWM=TAUWMT/((RHOS-RHOW)*GRAV*D50)
! WAVE-CURRENT SHIELDS NUMBER (MEAN)
      THETACW=TAUCWT/((RHOS-RHOW)*GRAV*D50)
! WAVE-CURRENT SHIELDS NUMBER (MAXIMUM)
      THETACWM=TAUCWMT/((RHOS-RHOW)*GRAV*D50)
! COEFFICIENT IN REFERENCE CONCENTRATION EQUATION
      ACRCW=0.0035D0*EXP(-0.3d0*DST)
! REFERENCE CONCENTRATION
      IF(THETACWM .GT. 0.2D0*THETACR)THEN
         CRCW=ACRCW*THETACW*EXP(-4.5d0*THETACR/THETACWM)
      ELSE
         CRCW=0.D0
      ENDIF
! CALCULATE SUSPENDED TRANSPORT RATE
      if (dep .gt. 0.D0 .and. epscw .gt. 0.D0 .and. ws .gt. 0.D0) then
         BDpart = EPSCW/WS*(1.d0-EXP(-WS*DEP/EPSCW))
         QSS=UC*CRCW*BDpart
      else
       QSS = 0.d0
      end if

      RETURN
      END SUBROUTINE !SUSPLUND

! *****************************************************************
!
!     LUND-CIRP FORMULA FOR PREDICTING BED LOAD TRANSPORT 
!
!     CODED BY M. LARSON
!     2005-11-01
!
!******************************************************************
      SUBROUTINE ZBEDLUND(RHOS,RHOW,TAUCTB,TAUCWTB,TAUCWMTB,TAUCR,QBS)
      USE flow_def, only: GRAV
      
! VARIABLE DECLARATIONS
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 RHOS,RHOW
! EXPLICIT DEFINITIONS - MEB 12/11/2013
      REAL*8 A,B,TAUCWMTB,TAUCR,QBS,TAUCTB,TAUCWTB

! CONSTANTS
      !GRAV=9.81d0

! TRANSPORT COEFFICIENTS
      A=12.d0
      B=4.5d0

! BED LOAD TRANSPORT
      if (taucwmtb .gt. 0.2D0*taucr) then
         QBS=A*SQRT(TAUCTB/RHOW)*TAUCWTB*EXP(-B*TAUCR/TAUCWMTB)
         QBS=QBS/(RHOS-RHOW)/GRAV
      else
       QBS = 0.d0
      end if

      RETURN
      END SUBROUTINE !BEDLUND

!****************************************
    SUBROUTINE ZPRINT_HEADER
    !use EXP_Global_def, ONLY: VERSION,REVISION,RDATE  !synchronize these variables with implicit for the merged code. 
    use comvarbl, ONLY: VERSION,REVISION,RDATE
    use diag_def, ONLY: DGUNIT

    IMPLICIT NONE
    INTEGER I, IUNIT
    
    open(DGUNIT,FILE='CMS_DIAG.TXT',STATUS='unknown') 
    IUNIT = 6
    DO I=6,7
      IF (I.EQ.7) IUNIT=DGUNIT
      WRITE(IUNIT,*)
      WRITE(IUNIT,7009)
      WRITE(IUNIT,7011)
      WRITE(IUNIT,7012)
      WRITE(IUNIT,7013)
      WRITE(IUNIT,7014)
      WRITE(IUNIT,7016) RDATE
      WRITE(IUNIT,7009)
      WRITE(IUNIT,*)
    ENDDO
      
7009 FORMAT(' *******************************************************')
7011 FORMAT('              U.S. Army Corps of Engineers              ')
7012 FORMAT('            Coastal Inlets Research Program             ')
7013 FORMAT('     Explicit Circulation & Morphology Change Model     ')
7014 FORMAT('           CMS-Flow-HS, Version 8.4                     ')
7016 FORMAT('               Last updated - ',A10)

    RETURN
    END SUBROUTINE  

!************************************************************
    SUBROUTINE PRINT_STATUS
    use EXP_Global_def, only: iadv, imix, advcoef, dt, isedform, thetac, a0, dtmorph, dtsed
    use flow_def, only: rhow
    use comvarbl, only: tmax, rampdur
    use sed_def,  only: rhosed,poros,scalesus,scalemorph,scalebed,scalemorph
    use met_def,  only: wkappa
    use hot_def,  only: coldstart
    use diag_def, only: DGUNIT, DGFILE
    
    IMPLICIT NONE    
    INTEGER I, IIDAY, IIHRS, IUNIT !, IHT, ITHR, IPR 
    character*3 cc1,cc2,cc3
    
    !PERFORM THESE IF STATEMENTS ONLY ONCE
    CC1 = 'OFF'; CC2 = 'OFF' ; CC3 = 'OFF'
    IF (IADV  == 1) CC1='ON ' ; IF (IMIX  == 1) CC2='ON ' 
    !IF (WALLFRIC) CC3='ON '
    !IF (DEVELOPER)        ADVANCED=.TRUE.
    !IF (.NOT.DO_AVALANCHE)ADVANCED=.TRUE.
    !IF (NO_BULK_FILES)    ADVANCED=.TRUE.
    !IF (OUTPUT_WAVEINFO)  ADVANCED=.TRUE.
    !IF (DO_SMOOTHDH)      ADVANCED=.TRUE.
    !IF (ICNTAMUTOUT)      ADVANCED=.TRUE.
   ! IF (WRITE_OLD_FORMAT) ADVANCED=.TRUE.
    !IF (RAM_OUTPUT)       ADVANCED=.TRUE.
    !if (WKAPPA.ne.0.40)   ADVANCED=.TRUE.
    IF (TMAX.GT.24) THEN
      IIDAY=INT(TMAX/24)
      IIHRS=MOD(INT(TMAX),24)
    ENDIF     
    
    !IPR = NCPU 
       
    IUNIT = 6

    DO I=6,7
      IF (I.EQ.7) THEN
        IUNIT=DGUNIT
        open(IUNIT,file=dgfile,access='append') 
      ENDIF
      WRITE (IUNIT,*)' '
      WRITE (IUNIT,*)'*********************'  ! 7/13/06
      WRITE (IUNIT,*)'* SETUP INFORMATION *'
      WRITE (IUNIT,*)'*********************'
      WRITE (IUNIT,*)' '
    
!  ! OPEN MP SECTION
!        ITHR = 1
!      IF (NTHR .GE. 1 ) THEN
!        WRITE (IUNIT,8011)NTHR
!!        IHT=int(NTHR/NCPU)                   !07/22/2009 - Brown
!!        IF (IHT.GT.1) WRITE(IUNIT,8017)      !07/22/2009 - Brown
!!        IF (IHT.EQ.1) WRITE(IUNIT,8018)      !07/22/2009 - Brown
!        IF (NTHR .GE. 1) THEN
!          IF (ITHREADS .GT. NTHR) THEN
!            WRITE(IUNIT,8014) ITHREADS,NTHR
!            ITHR = NTHR
!          ELSE
!            WRITE (IUNIT,8012) ITHREADS
!            ITHR = ITHREADS
!          ENDIF
!!          IPR  = ITHR/IHT                    !07/22/2009 - Brown
!        ENDIF
!        IF (IPROCS .GE. 1) THEN
!          IF (IPROCS .GT. NCPU) THEN
!            WRITE(IUNIT,8015) IPROCS,NCPU
!            IPR = NCPU
!          ELSE
!            WRITE (IUNIT,8013) IPROCS
!            IPR = IPROCS
!          ENDIF
!          ITHR = IPR
!        ENDIF
!!$      call omp_set_num_threads(ITHR)                              !set number of threads as defined in CARD file
!        WRITE (IUNIT,8016) ITHR
!        WRITE (IUNIT,*) ' '
!      ENDIF
!
!8011  FORMAT(' OPENMP PARALLELIZATION ENABLED - ',I0,' CPU(s) and ',I0,' threads available')
!8012  FORMAT('   Total number of threads requested - ',I0)
!8013  FORMAT('   Total number of processors requested - ',I0)
!8014  FORMAT('   Total number of threads requested, ',I0,', EXCEEDS maximum of ',I0)
!8015  FORMAT('   Total number of processors requested, ',I0,', EXCEEDS maximum of ',I0)
!8016  FORMAT('   CMS-Flow will use ',I0,' threads')
!!8017  FORMAT('   Hyperthreading (2+ threads / CPU) IS available')
!!8018  FORMAT('   Hyperthreading (2+ threads / CPU) is NOT available')
        
! HYDRODYNAMICS SECTION
      WRITE (IUNIT,*)'HYDRODYNAMICS'
      WRITE (IUNIT,8005) CC3,CC2
      WRITE (IUNIT,8003) CC1,ADVCOEF
     ! WRITE (IUNIT,8004) DRYDEP,ISMOOTH
      IF (TMAX.LE.24) THEN
        WRITE (IUNIT,8006) TMAX, RAMPDUR
      ELSE
        WRITE (IUNIT,8009) IIDAY,IIHRS,RAMPDUR
      ENDIF
      IF (coldstart) THEN
        WRITE (IUNIT,8007) 
      ELSE
       ! WRITE (IUNIT,8008) SMELAPSE
      ENDIF
      !WRITE (IUNIT,8010) EDDYFAC
      WRITE (IUNIT,8001) DT
      WRITE (IUNIT,*) ' '
  
! SEDIMENT TRANSPORT SECTION
      IF (isedform==1) THEN
        WRITE (IUNIT,8021)
        SELECT CASE (ISEDFORM)
        CASE (1)  ! Equilibrium - Total Load - Watanabe
          WRITE(IUNIT,8022)
!          if (D50_READ) WRITE (IUNIT,8030) singleD50*1000.d0,DCOEFF
!          if (VD50_READ) WRITE (IUNIT,8031)DCOEFF
          WRITE (IUNIT,8032) RHOSED,RHOW
          WRITE (IUNIT,8035) THETAC,A0
          WRITE (IUNIT,8027) SCALEMORPH
          WRITE (IUNIT,8034) 1-(1/POROS)
        CASE (2)  ! Equilibrium - Total Load - Lund-CIRP
          WRITE(IUNIT,8023)
!          if (D50_READ)    WRITE (IUNIT,8030) singleD50*1000.d0,DCOEFF
!          if (VD50_READ)   WRITE (IUNIT,8031)DCOEFF
          WRITE (IUNIT,8032) RHOSED,RHOW
          WRITE (IUNIT,8037) SCALEBED,SCALESUS
          WRITE (IUNIT,8027) SCALEMORPH
          !IF (IRIPPLE==1) WRITE (IUNIT,8036)WTEMP,'ON '
          !IF (IRIPPLE==0) WRITE (IUNIT,8036)WTEMP,'OFF'
          WRITE (IUNIT,8034) 1-(1/POROS)
        CASE (3)  ! Equilibrium - Advection-Diffusion
        !  IF (IFUNC==1)   WRITE (IUNIT,8024)  !Exponential concentration profile
       !   IF (IFUNC==2)   WRITE (IUNIT,8025)  !Van Rijn concentration profile
       !   IF (IFUNC==3)   WRITE (IUNIT,8026)  !Lund-CIRP concentration profile
!          if (D50_READ)    WRITE (IUNIT,8030) singleD50*1000.d0,DCOEFF
 !         if (VD50_READ)   WRITE (IUNIT,8031)DCOEFF
          WRITE (IUNIT,8032) RHOSED,RHOW
          WRITE (IUNIT,8037) SCALEBED,SCALESUS
          WRITE (IUNIT,8027) SCALEMORPH
          !IF (IRIPPLE==1) WRITE (IUNIT,8036)WTEMP,'ON '
          !IF (IRIPPLE==0) WRITE (IUNIT,8036)WTEMP,'OFF'
          WRITE (IUNIT,8034) 1-(1/POROS)
        CASE (4)  ! Non-Equilibrium - Advection-Diffusion
!          if (netcapac==1) WRITE (IUNIT,8040)  !Watanabe capacity formula
!          if (netcapac==2) WRITE (IUNIT,8041)  !Lund-CIRP capacity formula
!          if (netcapac==3) WRITE (IUNIT,8042)  !Van Rijn capacity formula
!!          if (D50_READ)    WRITE (IUNIT,8030) singleD50*1000.d0,DCOEFF
!!          if (VD50_READ)   WRITE (IUNIT,8031) DCOEFF
!          WRITE (IUNIT,8032) RHOSED,RHOW
!          WRITE (IUNIT,8037) SCALEBED,SCALESUS
!          WRITE (IUNIT,8027) SCALEMORPH
!          !IF (IRIPPLE==1) WRITE (IUNIT,8036)WTEMP,'ON '
!          !IF (IRIPPLE==0) WRITE (IUNIT,8036)WTEMP,'OFF'
!          WRITE (IUNIT,8034) 1-(1/POROS)
        CASE (5)  ! Cohesive
          WRITE (IUNIT,8028)
 !         if (D50_READ)    WRITE (IUNIT,8030) singleD50*1000.d0,DCOEFF
 !         if (VD50_READ)   WRITE (IUNIT,8031)DCOEFF
          WRITE (IUNIT,8032) RHOSED,RHOW
          WRITE (IUNIT,8037) SCALEBED,SCALESUS
          WRITE (IUNIT,8027) SCALEMORPH
          !IF (IRIPPLE==1) WRITE (IUNIT,8036)WTEMP,'ON '
          !IF (IRIPPLE==0) WRITE (IUNIT,8036)WTEMP,'OFF'
          WRITE (IUNIT,8034) 1-(1/POROS)
        END SELECT
        IF (DTMORPH.LT.900) THEN
          WRITE (IUNIT,8039) DTSED, DTMORPH
        ELSE
          WRITE (IUNIT,8038) DTSED, DTMORPH/3600.d0
        ENDIF
        WRITE (IUNIT,*) ' '
      ENDIF
!      IF (SALTSIM) THEN
!        WRITE (IUNIT,8048)
!        WRITE (IUNIT,8049) SALTparms.IC
!        WRITE (IUNIT,8047) DTSALT
!        WRITE (IUNIT,*) ' '
!      ENDIF
!     
!! FORCING SECTION
!      WRITE(IUNIT,8050)
!      IF (NDRIVER.GT.0)     WRITE(IUNIT,8051) NDRIVER
!      IF (NQDRIVER.GT.0)    WRITE(IUNIT,8052) NQDRIVER 
!      IF (NMDRIVER.GT.0)    WRITE(IUNIT,8053) NMDRIVER
!      IF (NMVDRIVR.GT.0)    WRITE(IUNIT,8054) NMVDRIVR
!      IF (TREAD)            WRITE(IUNIT,8055)
!      IF (WAVFIL.NE.'NONE') WRITE(IUNIT,8056)
!      IF (RADFIL.NE.'NONE') WRITE(IUNIT,8057)
!      IF (WINDREAD)         WRITE(IUNIT,8058)
!      WRITE(IUNIT,*) ' '
!      
!! ADVANCED DIRECTIVES SECTION
!      IF (ADVANCED) THEN
!        WRITE (IUNIT,8090) 
!        IF (DEVELOPER)        WRITE (IUNIT,8091) 
!        IF (.not.DO_AVALANCHE)WRITE (IUNIT,8092) 
!        IF (NO_BULK_FILES)    WRITE (IUNIT,8093) 
!        IF (OUTPUT_WAVEINFO)  WRITE (IUNIT,8094) 
!        IF (DO_SMOOTHDH .and. ISEDTRAN==1) WRITE (IUNIT,8095) 
!        IF (ICNTAMUTOUT)      WRITE (IUNIT,8096) 
!        IF (WRITE_OLD_FORMAT) WRITE (IUNIT,8097)
!        IF (RAM_OUTPUT)       write (IUNIT,8099)
!        if (WKAPPA.ne.0.40)  WRITE (IUNIT,8100) WKAPPA
!        WRITE (IUNIT,*) ' '
!      ENDIF
!
!! ADDITIONAL INFORMATION SECTION
!      WRITE(IUNIT,8070)
!      IF (NFLOW>0 .OR. NTS>0 .OR. NQS>0) THEN
!        WRITE(IUNIT,8074)
!        IF (NFLOW>0) WRITE(IUNIT,8075) NFLOW
!        IF (NTS>0)   WRITE(IUNIT,8076) NTS
!        IF (NQS>0.AND.ISEDTRN==1) WRITE(IUNIT,8077) NQS
!      ENDIF
!      WRITE(IUNIT,8078) NCELLS
    ENDDO
    
    !ITHREADS = ITHR
    !IPROCS = IPR

!Hydrodynamics
8000  FORMAT('  Iterations      - ',I0)
8001  FORMAT('  Solution scheme - EXPLICIT',T40,'Timestep               -  ',F0.3,' secs')
8003  FORMAT('  Advection       - ',A3,T40,'Advection Extrap Coeff -    ',F0.3)
8004  FORMAT('  Drying depth    - ',F0.3,' m ',T40,'WSEL Smoothing Iter    -    ',I0)
8005  FORMAT('  Wall Friction   - ',A3,T40,'Mixing                 - ',A3)
8007  FORMAT('  Run type        - COLDSTART')
8008  FORMAT('  Run type        - HOTSTART',T40,'Starting time          - ',F0.3,' hrs')
8006  FORMAT('  Run Duration    - ',F0.3,' hrs',T40,'Ramp Duration          -   ',F0.3,' hrs')
8009  FORMAT('  Run Duration    - ',I0,' days, ',I0,' hrs',T40,'Ramp Duration   - ',F0.3,' hrs')
8010  FORMAT('  Eddy Visc. Multiplier  - ',F0.2)

!Sediment Transport      
8021  FORMAT(' SEDIMENT TRANSPORT')
8022  FORMAT('  WATANABE formulation')
8023  FORMAT('  LUND-CIRP formulation')
8024  FORMAT('  A-D formulation with EXPONENTIAL ')
8025  FORMAT('  A-D formulation with VAN RIJN')
8026  FORMAT('  A-D formulation with LUND-CIRP')
8027  FORMAT('  Morphology Accel. Factor - ',F0.2)
8028  FORMAT('  Cohesive Sediment Transport')
!8029  Future use
8030  FORMAT('  Using CONSTANT D50    - ',F0.3,' mm ',T40,'Slope Coefficient   -   ',F0.2)
8031  FORMAT('  Using VARIABLE D50',T40,'Slope Coefficient   -   ',F0.2)
8032  FORMAT('  Sediment Density      - ',F0.2,T40,'Water Density       - ',F0.2)
8033  FORMAT('  Slope Coefficient   - ',F0.2)
8034  FORMAT('  Porosity                - ',F0.2)
8035  FORMAT('  Theta Critical          - ',F0.2,T40,'A0 parameter          - ',F0.2)
8036  FORMAT('  Water temperature       - ',F0.2,T40,'Ripple Calculations - ',A3)
8037  FORMAT('  Scaling Factor: Bedload - ',F0.2,T40,'Suspended Load      -   ',F0.2)
8038  FORMAT('  Timestep: Sed Transport - ',F0.3,' s',T40,'Morphology Update   -   ',F0.2,' hr')
8039  FORMAT('  Timestep: Sed Transport - ',F0.3,' s',T40,'Morphology Update   - ',F0.2,' s')
8040  FORMAT('  Non-equilibrium Sediment Transport',T40,'Lund-CIRP capacity formulation')
8041  FORMAT('  Non-equilibrium Sediment Transport',T40,'Watanabe capacity formulation')
8042  FORMAT('  Non-equilibrium Sediment Transport',T40,'Van Rijn capacity formulation')
8048  FORMAT(' SALINITY TRANSPORT')
8049  FORMAT('  Initial Concentration  - ',F0.2,' ppt')
8047  FORMAT('  Salinity Calc Timestep - ',F0.2,' seconds')


!Forcing Information      
8050  FORMAT(' FORCING INFORMATION')
8051  FORMAT('  ',I1,' Water level - spatially CONSTANT boundary condition(s)')
8052  FORMAT('  ',I1,' Normal flow boundary condition(s)')
8053  FORMAT('  ',I1,' Water level - spatially VARIABLE boundary condition(s)')
8054  FORMAT('  ',I1,' Water level and velocity - spatially VARIABLE boundary condition(s)')
8055  FORMAT('  Tidal constituent boundary condition defined')
8056  FORMAT('  Wave height, period, direction and dissipation will be used')
8057  FORMAT('  Radiation stress forcing will be used')
8058  FORMAT('  Wind stress forcing will be used')

!Advanced Directives
8090  FORMAT(' ADVANCED DIRECTIVES')
8091  FORMAT('  CMS Developer, access to restricted features is allowed')
8092  FORMAT('  Avalanching routines will NOT be entered')
8093  FORMAT('  Bulk wave information files will NOT be used')
8094  FORMAT('  Wave information will be output to global solution files')
8095  FORMAT('  Using 9-cell smoothing for DH values during sediment transport')
8096  FORMAT('  Additional output for Friction and Diffusion')
8097  FORMAT('  Will write basic old format to file named "GRID.DAT"')
8098  FORMAT('  Using Eddy Viscosity Factor of ',F0.2)
8099  FORMAT('  Output for RAM to files "VELOCITY.DAT" and "DEPTH.DAT"')
8100  FORMAT('  Using non-default for Wind Drag, Coefficient = ',F0.2)

!Additional Information
8070  FORMAT(' ADDITIONAL INFORMATION')
!8071  FORMAT('  ')
!8072  FORMAT('  ')
!8073  FORMAT('  ')
8074  FORMAT('  Observation Cell output in ASCII format')
8075  FORMAT('  - X/Y flow rate cells        - ',I0)
8076  FORMAT('  - U/V/ETA time series cells  - ',I0)
8077  FORMAT('  - Qx/Qy transport rate cells - ',I0)
8078  FORMAT('  There are ',I0,' Computational cells')

    close(dgunit)
    
    RETURN
    END SUBROUTINE !PRINT_STATUS
    
!   ***************************************************************  
    SUBROUTINE SMOOTH (EL, NCELLS, NTIMES)
    use prec_def, only: ikind

    ! added 07/25/06 by meb
    ! This routine smooths the spatial values of the passed WSEL array
    ! NTIMES is used as a factor to provide for additional smoothing

    INTEGER NCELLS,NTIMES,I,J
    REAL(IKIND) EL(NCELLS), SUM
    REAL(IKIND),ALLOCATABLE :: VALUES (:)

    ALLOCATE (VALUES(NCELLS))
    SUM=0.D0

    DO I=1,NCELLS
      VALUES(I) = EL(I)
    ENDDO
      
    if (NTIMES.LT.0) THEN     !IF NTIMES IS NEGATIVE, FORCE ALL CELLS TO USE THE AVERAGE ELEVATION 
      DO I=1,NCELLS
        SUM=SUM+VALUES(I)
      ENDDO
      DO I=1,NCELLS
        EL(I) = SUM/NCELLS
      ENDDO
    ELSE                      !OTHERWISE, ITERATE THE MOVING 3-CELL AVERAGE "NTIMES" TIMES.
      DO J=1,NTIMES
        DO I=2,NCELLS-1
          EL(I)=(VALUES(I-1)+VALUES(I)+VALUES(I+1))/3
        ENDDO
      ENDDO
    ENDIF
    DEALLOCATE (VALUES)
    
    RETURN
    END SUBROUTINE
    
!*****************************************************************       
    subroutine exp_timing
! Calculates the cumulative time, number of iterations and
! time step
! written by Alex Sanchez, USACE-CHL
! modified by Mitchell Brown, USACE-CHL for explicit operation
!*****************************************************************
#include "CMS_cpp.h"
    use comvarbl,  only: ntime, etime, timesecs, ctime, ctime1, stimet, mtime, deltime, timehrs, ramp, rampdur
    use const_def, only: pi
    use time_lib,  only: time_jul, ramp_func
    use cms_def,   only: timestart, timenow
    use solv_def,  only: iconv
    use prec_def,  only: ikind
    implicit none
    
    integer :: ierr
    real(8) :: dtimetemp,rtime
    real(8) :: timedur,timerem,timelast,speed,err,timeint
    character(len=100) :: str
    
    ntime = ntime + 1  !Time step iteration counter
    
    if (etime > 0 .and. timesecs > 0.0) then
      if(ntime==11 .or. mod(timesecs,etime)==0)then
        timelast = timenow
        timenow = time_jul()
        timeint = timenow - timelast        !Time interval between last speed check [sec]
        timedur = timenow - timestart       !Total simulation clock time [sec]
        speed = dble(ctime-ctime1)/timeint  !Note: computed using last speed check time interval
        ctime1 = ctime
        timerem = dble(stimet-ctime)/speed

       !Modified by Chris Reed - submitted 10/14/16        
        !call time_sec2str(timesecs,str)
        !write(msg,'(A,A)',iostat=ierr)     'Elapsed Simulation Time:  ',trim(str)
        !call diag_print_message(' ',msg)
        !call time_sec2str(timedur,str)
        !write(msg,'(A,A)',iostat=ierr)     'Elapsed Clock Time:       ',trim(str)
        !call diag_print_message(msg)
        !write(msg,'(A,F10.3)',iostat=ierr) 'Computational Speed:      ',speed
        !call diag_print_message(msg)
        !call time_sec2str(timerem,str)
        !write(msg,'(A,A)',iostat=ierr)     'Remaining Clock Time:     ',trim(str)
        !call diag_print_message(msg)
        !call diag_print_message(' ')
      endif
    endif
    
    mtime = mtime + 1           
    timesecs = timesecs + deltime
    ctime = real(timesecs,kind=ikind)
    timehrs = ctime/3600.0_ikind
    ramp = ramp_func(timehrs,rampdur)

    return
    end subroutine exp_timing        