!!**********************************************************
!      program gctp
!!Uses NAD27 Datums      
!!***********************************************************      
!     GCTP  GENERAL CARTOGRAPHIC COORDINATES TRANSFORMATION PACKAGE
!                   ADJLZ0
!                                                  FIPS CODES 6-10-82
! ** U.S.G.S. GENERAL MAP PROJECTION PACKAGE ...... DR. A. A. ELASSAL **
! ** MODULE I                VERSION 1.0.0            NOVEMBER 1,1980 **
! **********************************************************************
      DOUBLE PRECISION FUNCTION ADJLZ0 (LON)
!
! FUNCTION TO ADJUST LONGITUDE ANGLE TO MODULE 180 DEGREES.
!
      IMPLICIT REAL*8 (A-Z)
      DATA TWO,PI /2.0D0,3.14159265358979323846D0/
!
  020 ADJLZ0 = LON
      IF (DABS(LON) .LE. PI) RETURN
      TWOPI = TWO * PI
      LON = LON - DSIGN (TWOPI,LON)
      GO TO 020
!
      END

!                   AL01Z0
! **********************************************************************
! ** U.S.G.S. GENERAL MAP PROJECTION PACKAGE ...... DR. A. A. ELASSAL **
! ** MODULE I                VERSION 1.0.0            NOVEMBER 1,1980 **
! **********************************************************************
      SUBROUTINE AL01Z0 (CRDIN,CRDOUT,FLAG)
!
! SUBROUTINE TO COMPUTE TRANSFORMATION BETWEEN GEOGRAPHIC AND
! ALASKA STATE ZONE NO. 1.
! FLAG = 0, MEANS PLANE TO GEOGRAPHIC.
! FLAG = 1, MEANS GEOGRAPHIC TO PLANE.
!
      IMPLICIT REAL*8 (A-Z)
      INTEGER*4 FLAG
      DIMENSION CRDIN(1),CRDOUT(1)
      DATA B,C,D /1.00029977273D0,0.00447599131D0,6386352.67013D0/
      DATA F,G,E /0.327015517176D0,0.945018968871D0,0.082271854223003D0/
      DATA PI,EPS /3.141592653589793D0,2.718281828459045D0/
      DATA LO /1.771754086D0/
      DATA C1,C2 /0.182880365761D0,0.243840487681D0/
      DATA C3,C4 /7000000.0D0,1000000.0D0/
      DATA ONE,TWO,FOUR /1.0D0,2.0D0,4.0D0/
      DATA MFEET,FALSE /3.280833333333D0,5000000.0D0/
!
      IF (FLAG .EQ. 0) GO TO 020
!
! GEOGRAPHIC TO STATE PLANE TRANSFORMATION.
!
      GEOG1 = DABS(CRDIN(1))
      GEOG2 = CRDIN(2)
      ESINP = E * DSIN(GEOG2)
      TANP = DTAN((PI / FOUR) + (GEOG2 / TWO))
      MU = DLOG(TANP) - (E / TWO) * DLOG((ONE + ESINP) / (ONE - ESINP))
      CON = B * MU + C
      CON1 = EPS**CON
      CON =-CON
      CON2 = EPS**CON
      P = (CON1 - CON2) / TWO
      Q = (CON1 + CON2) / TWO
      CON = B * (GEOG1 - LO)
      CON1 = DSIN(CON)
      CON2 = DCOS(CON)
      U = D * DATAN((G * P + F * CON1) / CON2)
      CON = F * P - G * CON1
      V = (D / TWO) * DLOG((Q + CON) / (Q - CON))
      CRDOUT(1) = MFEET * (-0.6D0 * U + 0.8D0 * V + FALSE)
      CRDOUT(2) = MFEET * ( 0.8D0 * U + 0.6D0 * V - FALSE)
      RETURN
!
! STATE PLANE TO GEOGRAPHIC TRANSFORMATION.
!
  020 U =-C1 * CRDIN(1) + C2 * CRDIN(2) + C3
      V = C2 * CRDIN(1) + C1 * CRDIN(2) - C4
      CON = V / D
      CON1 = EPS**CON
      CON =-CON
      CON2 = EPS**CON
      R = (CON1 - CON2) / TWO
      S = (CON1 + CON2) / TWO
      K1 = DSIN(U / D)
      K2 = DCOS(U / D)
      CON = F * R + G * K1
      MU = (ONE / (TWO * B)) * DLOG((S + CON) / (S - CON)) - (C / B)
      KI = TWO * DATAN(EPS**MU) - PI / TWO
      CON1 = DSIN(KI)
      CON2 = DCOS(KI)
      CRDOUT(2) = KI + (0.006761032571D0 + 0.000053172205D0 * CON2**2 + &
                  0.573027D-6 * CON2**4 + 0.7128D-8 * CON2**6) * CON1 * &
                  CON2
      CRDOUT(1) =-LO - (ONE / B) * DATAN((F * K1 - G * R) / K2)
      RETURN
!
      END

!                   AL29Z0
! **********************************************************************
! ** U.S.G.S. GENERAL MAP PROJECTION PACKAGE ...... DR. A. A. ELASSAL **
! ** MODULE I                VERSION 1.0.0            NOVEMBER 1,1980 **
! **********************************************************************
      SUBROUTINE AL29Z0 (CRDIN,CRDOUT,ZONE,FLAG)
!
! SUBROUTINE TO COMPUTE TRANSFORMATION BETWEEN GEOGRAPHIC AND
! ALASKA STATE ZONES NO. 2 THROUGH 9.
! FLAG = 0, MEANS PLANE TO GEOGRAPHIC.
! FLAG = 1, MEANS GEOGRAPHIC TO PLANE.
!
      IMPLICIT REAL*8 (A-Z)
      INTEGER*4 ZONE,FLAG,IND
      DIMENSION CRDIN(1),CRDOUT(1),CONST(2,8)
      DATA CONST /500000.0D0 , 511200.0D0,  &
                  500000.0D0 , 525600.0D0,  &
                  500000.0D0 , 540000.0D0,  &
                  500000.0D0 , 554400.0D0,  &
                  500000.0D0 , 568800.0D0,  &
                  700000.0D0 , 583200.0D0,  &
                  500000.0D0 , 597600.0D0,  &
                  600000.0D0 , 612000.0D0/
      DATA RADSEC /206264.806247D0/
      DATA ONE,TWO /1.0D0,2.0D0/
!
      IND = ZONE - 1
      C = CONST(1,IND)
      CM = CONST(2,IND)
      IF (FLAG .EQ. 0) GO TO 020
!
! GEOGRAPHIC TO STATE PLANE TRANSFORMATION.
!
      GEOG1 = DABS(CRDIN(1))
      GEOG2 = CRDIN(2)
      C1 = DCOS(GEOG2)
      C2 = C1 * C1
      C3 = C2 * C2
      C4 = C2 * C3
      C5 = (CM - GEOG1 * RADSEC) * 1.0D-4
      C6 = C5 * C5
      C7 = C6 * C6
      C8 = DSQRT(ONE + 0.0068147849 * C2)
      C9 = DSQRT(ONE - C2) * C1
      CRDOUT(1) = C + (1017862.15D0 * C1 / C8) * C5 * (ONE -              &
                  3.91740509D-4 * C6 * (ONE - TWO * C2 - 0.681478D-2 *    &
                  C3) + 4.60382D-8 * C7 * (ONE - 20.0D0 * C2 +            &
                  23.6047D0 * C3 + 0.4907D0 * C4))
      CRDOUT(2) = 101.269278503D0 * (GEOG2 * RADSEC - 193900.05442D0      &
                  - (1052.893943D0 - 4.483386D0 * C2 + 2.3559D-2 * C3) *  &
                  C9) + (24673.6748D0 * C9 * C6 / C8) * (ONE +            &
                  1.958703D-4 * C6 * (-ONE + 6.0D0 * C2 + 6.133306D-2 *   &
                  C3 + 1.8577D-4 * C4) + 1.5346D-8 * C7 * (ONE -          &
                  60.0D0 * C2 + 117.75D0 * C3 + 4.089D0 * C4))
      RETURN
!
! STATE PLANE TO GEOGRAPHIC TRANSFORMATION.
!
  020 OMEGA = 193900.05442D0 + 0.00987466302498D0 * CRDIN(2)
      C1 = DCOS(OMEGA / RADSEC)
      C2 = C1 * C1
      C3 = C2 * C2
      PHI = OMEGA + (1047.546691D0 + 6.193011 * C2 + 5.0699D-2 * C3) * DSQRT(ONE - C2) * C1
      C1 = DCOS(PHI / RADSEC)
      C2 = C1 * C1
      C3 = C2 * C2
      C4 = (CRDIN(1) - C) * 1.0D-6
      C5 = C4 * C4
      C6 = C5 * C5
      C7 = ONE + 0.0068147849 * C2
      C8 = C7 * C7
      CRDOUT(2) = (PHI - 233.973645D0 * C5 * C8 * DSQRT((ONE / C2) -     &
                  ONE) * (ONE - 1.8905604D-4 * C5 * (1.9591113D0 +       &
                  (3.0D0 / C2) + 8.1359D-2 * C2 + 2.79D-4 * C3) +        &
                  1.42969D-8 * C6 * C7 * (15.5D0 + (45.0D0 / C3) -       &
                  (0.307D0 / C2) + 1.53D0 * C2))) / RADSEC
      CRDOUT(1) =-(CM - 9824.513072D0 * DSQRT(C7) * C4 / C1 *            &
                  (ONE - 3.7811208D-4 * C7 * C5 * (-TWO + C7 + (TWO /    &
                  C2)) + 4.2890624D-8 * C7 * C6 * (1.054 + (24.0D0 /     &
                  C3) - (20.0D0 / C2) - 1.36D-2 * C2))) / RADSEC
      RETURN
!
      END

!                   BLOCKD
! **********************************************************************
! ** U.S.G.S. GENERAL MAP PROJECTION PACKAGE ...... DR. A. A. ELASSAL **
! ** MODULE I                VERSION 1.0.0            NOVEMBER 1,1980 **
! **********************************************************************
      BLOCK DATA
!
! INITIALIZATION OF ELLIPSOID TO CLARK'S 1866 PARAMETERS.
!
      IMPLICIT REAL*8 (A-Z)
      INTEGER*4 IPEMSG,IPPARM
!
      COMMON /ELLPZ0/ AZ,EZ,ESZ,E0Z,E1Z,E2Z,E3Z
      COMMON /SPHRZ0/ AZZ
      COMMON /PRINZ0/ IPEMSG,IPPARM
!
      DATA AZ  /0.6378206400000000D+07/
      DATA EZ  /0.8227185422300323D-01/
      DATA ESZ /0.6768657997291094D-02/
      DATA E0Z /0.9983056818784341D+00/
      DATA E1Z /0.2542555507651308D-02/
      DATA E2Z /0.2698084527466011D-05/
      DATA E3Z /0.1003393903560134D+01/
!
      DATA AZZ /0.6370997000000000D+07/
!
      DATA IPEMSG /0/
      DATA IPPARM /0/
!
      END

!                   DMSPZ0
! **********************************************************************
! ** U.S.G.S. GENERAL MAP PROJECTION PACKAGE ...... DR. A. A. ELASSAL **
! ** MODULE I                VERSION 1.0.0            NOVEMBER 1,1980 **
! **********************************************************************
      DOUBLE PRECISION FUNCTION DMSPZ0 (DMS)
!
! SUBROUTINE TO CONVERT UNPACKED DMS TO PACKED DMS ANGLE
!
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*4 REL
      INTEGER*4 DMS(1)
      EQUIVALENCE (REL , INT)
      DATA CON1,CON2 /1000000.0D0,1000.0D0/
      DATA NEG /'-'/
!
      INT = DMS(4)
      CON = DFLOAT (DMS(2)) * CON1 + DFLOAT (DMS(3)) * CON2 + REL
      IF (DMS(1) .EQ. NEG) CON = - CON
      DMSPZ0 = CON
      RETURN
!
      END

!                   E0FNZ0
! **********************************************************************
! ** U.S.G.S. GENERAL MAP PROJECTION PACKAGE ...... DR. A. A. ELASSAL **
! ** MODULE I                VERSION 1.0.0            NOVEMBER 1,1980 **
! **********************************************************************
      DOUBLE PRECISION FUNCTION E0FNZ0 (ECCNTS)
!
! FUNCTION TO COMPUTE CONSTANT (E0).
!
      IMPLICIT REAL*8 (A-Z)
      DATA QUART,ONE,ONEQ,THREE,SIXT /0.25D0,1.0D0,1.25D0,3.0D0,16.0D0/
!
      E0FNZ0 = ONE - QUART * ECCNTS * (ONE + ECCNTS / SIXT * (THREE + ONEQ * ECCNTS))
!
      RETURN
      END

!                   E1FNZ0
! **********************************************************************
! ** U.S.G.S. GENERAL MAP PROJECTION PACKAGE ...... DR. A. A. ELASSAL **
! ** MODULE I                VERSION 1.0.0            NOVEMBER 1,1980 **
! **********************************************************************
      DOUBLE PRECISION FUNCTION E1FNZ0 (ECCNTS)
!
! FUNCTION TO COMPUTE CONSTANT (E1).
!
      IMPLICIT REAL*8 (A-Z)
      DATA CON1,CON2,CON3 /0.375D0,0.25D0,0.46875D0/
      DATA ONE /1.0D0/
!
      E1FNZ0 = CON1 * ECCNTS * (ONE + CON2 * ECCNTS * (ONE + CON3 * ECCNTS))
!
      RETURN
      END

!                   E2FNZ0
! **********************************************************************
! ** U.S.G.S. GENERAL MAP PROJECTION PACKAGE ...... DR. A. A. ELASSAL **
! ** MODULE I                VERSION 1.0.0            NOVEMBER 1,1980 **
! **********************************************************************
      DOUBLE PRECISION FUNCTION E2FNZ0 (ECCNTS)
!
! FUNCTION TO COMPUTE CONSTANT (E2).
!
      IMPLICIT REAL*8 (A-Z)
      DATA CON1,CON2 /0.05859375D0,0.75D0/
      DATA ONE /1.0D0/
!
      E2FNZ0 = CON1 * ECCNTS * ECCNTS * (ONE + CON2 * ECCNTS)
!
      RETURN
      END

!                   E3FNZ0
! **********************************************************************
! ** U.S.G.S. GENERAL MAP PROJECTION PACKAGE ...... DR. A. A. ELASSAL **
! ** MODULE I                VERSION 1.0.0            NOVEMBER 1,1980 **
! **********************************************************************
      DOUBLE PRECISION FUNCTION E3FNZ0 (ECCENT)
!
! FUNCTION TO COMPUTE CONSTANT (E3).
!
      IMPLICIT REAL*8 (A-Z)
      DATA ONE /1.0D0/
!
      CON = ONE + ECCENT
      COM = ONE - ECCENT
      E3FNZ0 = DSQRT ((CON ** CON) * (COM ** COM))
!
      RETURN
      END

! ****                                                             *****
! ** U.S.G.S. GENERAL MAP PROJECTION PACKAGE ...... JOHN F. WAANANEN  **
! ** MODULE I                VERSION 1.0.0              APRIL 6,1981  **
! ****                                                             *****
      SUBROUTINE GTPZ0(CRDIN,INSYS,JNZONE,TPARIN,INUNIT,INSPH,IPR,JPR,  &
                       CRDIO,IOSYS,IOZONE,TPARIO,IOUNIT,IOSPH,IFLG)
!
! GENERAL PROGRAM FOR TRANSFORMATION BETWEEN VARIOUS REFERENCE SYSTEMS
!     MODIFIED VERSION OF GTRNZ0 BY J.F. WAANANEN
!     SUBROUTINE GTPZ0 IS REQUIRED FOR PROGRAMS NO. L176 AND NO. L177
!
! INPUT ****                                                       ****
! CRDIN  : COORDINATES IN INPUT SYSTEM (2 DP WORDS ARRAY).
! INSYS  : CODE NUMBER OF INPUT COORDINATE SYSTEM (INTEGER).
!            =  0 , GEOGRAPHIC
!            =  1 , U T M
!            =  2 , STATE PLANE
!            =  3 , ALBERS CONICAL EQUAL-AREA
!            =  4 , LAMBERT CONFORMAL CONIC
!            =  5 , MERCATOR
!            =  6 , POLAR STEREOGRAPHIC
!            =  7 , POLYCONIC
!            =  8 , EQUIDISTANT CONIC
!            =  9 , TRANSVERSE MERCATOR
!            = 10 , STEREOGRAPHIC
!            = 11 , LAMBERT AZIMUTHAL EQUAL-AREA
!            = 12 , AZIMUTHAL EQUIDISTANT
!            = 13 , GNOMONIC
!            = 14 , ORTHOGRAPHIC
!            = 15 , GENERAL VERTICAL NEAR-SIDE PERSPECTIVE
!            = 16 , SINUSOIDAL
!            = 17 , EQUIRECTANGULAR (PLATE CARREE)
!            = 18 , MILLER CYLINDRICAL
!            = 19 , VAN DER GRINTEN I
!            = 20 , OBLIQUE MERCATOR (HOTINE)
!            = 21 , SPACE OBLIQUE MERCATOR
!
! INZONE : CODE NUMBER OF INPUT COORDINATE ZONE (INTEGER).
! TPARIN : PARAMETERS OF INPUT REFERENCE SYSTEM (15 DP WORDS ARRAY).
! INUNIT : CODE NUMBER OF UNITS OF MEASURE FOR INPUT COORDS (INTEGER).
!            = 0 , RADIANS.
!            = 1 , FEET.
!            = 2 , METERS.
!            = 3 , SECONDS OF ARC.
!            = 4 , DEGREES OF ARC.
! INSPH  : INPUT SPHEROID CODE.  SEE SPHDZ0 FOR PROPER CODES.
! IPR    : PRINTOUT FLAG FOR ERROR MESSAGES. 0=YES, 1=NO
! JPR    : PRINTOUT FLAG FOR PROJECTION PARAMETERS 0=YES, 1=NO
! OUTPUT ***
! IOSYS  : CODE NUMBER OF OUTPUT COORDINATE SYSTEM (INTEGER).
! IOZONE : CODE NUMBER OF OUTPUT COORDINATE ZONE (INTEGER).
! TPARIO : PARAMETERS OF OUTPUT REFERENCE SYSTEM (15 DP WORDS ARRAY).
! IOUNIT : CODE NUMBER OF UNITS OF MEASURE FOR OUTPUT COORDS (INTEGER).
! IOSPH  : OUTPUT SPHEROID CODE.  SEE SPHDZ0 FOR PROPER CODES.
! CRDIO  : COORDINATES IN OUTPUT REFERENCE SYSTEM (2 DP WORDS ARRAY).
! IFLG   : RETURN FLAG (INTEGER).
!            =    0 , SUCCESSFUL TRANSFORMATION.
!            = 1001 , ILLEGAL INPUT SYSTEM CODE.
!            = 1002 , ILLEGAL OUTPUT SYSTEM CODE.
!            = 1003 , ILLEGAL INPUT UNIT CODE.
!            = 1004 , ILLEGAL OUTPUT UNIT CODE.
!            = 1005 , INCONSISTANT UNIT AND SYSTEM CODES FOR INPUT.
!            = 1006 , INCONSISTANT UNIT AND SYSTEM CODES FOR OUTPUT.
!            = 1007 , ILLEGAL INPUT ZONE CODE.
!            = 1008 , ILLEGAL OUTPUT ZONE CODE.
!      OTHERWISE , ERROR CODE FROM PROJECTION COMPUTATIONAL MODULE.
!
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER*4 SYSUNT(22)
      DIMENSION CRDIN(2),CRDIO(2),TPARIN(15),TPARIO(15),COORD(2)
      DIMENSION DUMMY(2)
      COMMON/ERRMZ0/ IERROR
      COMMON/PRINZ0/ IP1,IP2
      COMMON/ELLPZ0/ AZ,EZ,ESZ,E0Z,E1Z,E2Z,E3Z
      COMMON/PROJZ0/IPROJ
!
      DATA SYSUNT / 0 , 2 , 1 , 19*2 /
      DATA MAXUNT,MAXSYS / 4 , 21 /
      DATA JFLAG/0/
!
!     SETUP
      IP1 = IPR
      IP2 = JPR
      IPROJ=INSYS
      IF (JFLAG.NE.0) GO TO 10
      EZ = 0.0D0
      ESZ = 0.0D0
      CALL SPHDZ0(0,DUMMY)
      JFLAG = 1
   10 INZONE = IABS(JNZONE)
      IOSPH  = IABS(IOSPH)
      INSPH  = IABS(INSPH)
!
! CHECK VALIDITY OF CODES FOR UNITS OF MEASURE AND REFERENCE SYSTEMS.
      IF (INSYS.GE.0 .AND. INSYS.LE.MAXSYS) GO TO 020
      IFLG = 1001
      RETURN
  020 IF (IOSYS.GE.0 .AND. IOSYS.LE.MAXSYS) GO TO 040
      IFLG = 1002
      RETURN
  040 IF (INUNIT.GE.0 .AND. INUNIT.LE.MAXUNT) GO TO 060
      IFLG = 1003
      RETURN
  060 IF (IOUNIT.GE.0 .AND. IOUNIT.LE.MAXUNT) GO TO 080
      IFLG = 1004
      RETURN
!
! CHECK CONSISTANCY BETEEN UNITS OF MEASURE AND REFERENCE SYSTEM.
  080 IUNIT = SYSUNT(INSYS + 1)
      CALL UNTFZ0 (INUNIT,IUNIT,FACTOR,IFLG)
      IF (IFLG .EQ. 0) GO TO 100
      IFLG = 1005
      RETURN
  100 COORD(1) = FACTOR * CRDIN(1)
      COORD(2) = FACTOR * CRDIN(2)
      IUNIT = SYSUNT(IOSYS + 1)
      CALL UNTFZ0 (IUNIT,IOUNIT,FACTOR,IFLG)
      IF (IFLG .EQ. 0) GO TO 120
      IFLG = 1006
      RETURN
  120 IF (INSYS.NE.IOSYS.OR.INZONE.NE.IOZONE.OR.JNZONE.LE.0) GO TO 140
      IF (IOSPH.NE.INSPH) GO TO 140
      CRDIO(1) = FACTOR * COORD(1)
      CRDIO(2) = FACTOR * COORD(2)
      RETURN
!
! COMPUTE TRANSFORMED COORDINATES AND ADJUST THEIR UNITS.
  140 IF (INSYS .EQ. 0) GO TO 520
      IF (INZONE.GT.60 .OR. INSYS.EQ.1) GO TO 200
      IFLG = 1007
      RETURN
!
! INVERSE TRANSFORMATION.
  200 IPROJ=INSYS
      IF (INSYS.GT.2) CALL SPHDZ0(INSPH,TPARIN)
      GO TO (210,220,230,240,250,260,270,280,290,300,310,320,330,340,350,360,370,380,390,400,410),INSYS
  210 IF (INZONE.EQ.0.AND.TPARIN(1).NE.0.0D0) GO TO 211
      TPARIN(1) = 1.0D6*DFLOAT(6*INZONE-183)
      TPARIN(2) = DSIGN(4.0D7,DFLOAT(JNZONE))
!  211 CALL SPHDZ0(INSPH,DUMMY)
!      TPARIN(14) = DUMMY(1)
!      TPARIN(15) = DUMMY(2)
  211 CALL IS1AZ0(JNZONE,TPARIN)
      IF (IERROR .NE. 0) GO TO 500
      CALL PI1AZ0 (COORD,CRDIO)
      GO TO 500
  220 CALL IS02Z0 (INZONE,TPARIN)
      IF (IERROR .NE. 0) GO TO 500
      CALL PI02Z0 (COORD,CRDIO)
      GO TO 500
  230 CALL IS03Z0 (INZONE,TPARIN)
      IF (IERROR .NE. 0) GO TO 500
      CALL PI03Z0 (COORD,CRDIO)
      GO TO 500
  240 CALL IS04Z0 (INZONE,TPARIN)
      IF (IERROR .NE. 0) GO TO 500
      CALL PI04Z0 (COORD,CRDIO)
      GO TO 500
  250 CALL IS05Z0 (INZONE,TPARIN)
      IF (IERROR .NE. 0) GO TO 500
      CALL PI05Z0 (COORD,CRDIO)
      GO TO 500
  260 CALL IS06Z0 (INZONE,TPARIN)
      IF (IERROR .NE. 0) GO TO 500
      CALL PI06Z0 (COORD,CRDIO)
      GO TO 500
  270 CALL IS07Z0 (INZONE,TPARIN)
      IF (IERROR .NE. 0) GO TO 500
      CALL PI07Z0 (COORD,CRDIO)
      GO TO 500
  280 CALL IS08Z0 (INZONE,TPARIN)
      IF (IERROR .NE. 0) GO TO 500
      CALL PI08Z0 (COORD,CRDIO)
      GO TO 500
  290 CALL IS09Z0 (INZONE,TPARIN)
      IF (IERROR .NE. 0) GO TO 500
      CALL PI09Z0 (COORD,CRDIO)
      GO TO 500
  300 CALL IS10Z0 (INZONE,TPARIN)
      IF (IERROR .NE. 0) GO TO 500
      CALL PI10Z0 (COORD,CRDIO)
      GO TO 500
  310 CALL IS11Z0 (INZONE,TPARIN)
      IF (IERROR .NE. 0) GO TO 500
      CALL PI11Z0 (COORD,CRDIO)
      GO TO 500
  320 CALL IS12Z0 (INZONE,TPARIN)
      IF (IERROR .NE. 0) GO TO 500
      CALL PI12Z0 (COORD,CRDIO)
      GO TO 500
  330 CALL IS13Z0 (INZONE,TPARIN)
      IF (IERROR .NE. 0) GO TO 500
      CALL PI13Z0 (COORD,CRDIO)
      GO TO 500
  340 CALL IS14Z0 (INZONE,TPARIN)
      IF (IERROR .NE. 0) GO TO 500
      CALL PI14Z0 (COORD,CRDIO)
      GO TO 500
  350 CALL IS15Z0 (INZONE,TPARIN)
      IF (IERROR .NE. 0) GO TO 500
      CALL PI15Z0 (COORD,CRDIO)
      GO TO 500
  360 CALL IS16Z0 (INZONE,TPARIN)
      IF (IERROR .NE. 0) GO TO 500
      CALL PI16Z0 (COORD,CRDIO)
      GO TO 500
  370 CALL IS17Z0 (INZONE,TPARIN)
      IF (IERROR .NE. 0) GO TO 500
      CALL PI17Z0 (COORD,CRDIO)
      GO TO 500
  380 CALL IS18Z0 (INZONE,TPARIN)
      IF (IERROR .NE. 0) GO TO 500
      CALL PI18Z0 (COORD,CRDIO)
      GO TO 500
  390 CALL IS19Z0 (INZONE,TPARIN)
      IF (IERROR .NE. 0) GO TO 500
      CALL PI19Z0 (COORD,CRDIO)
      GO TO 500
  400 CALL IS20Z0 (INZONE,TPARIN)
      IF (IERROR .NE. 0) GO TO 500
      CALL PI20Z0 (COORD,CRDIO)
      GO TO 500
  410 CALL IS21Z0 (INZONE,TPARIN)
      IF (IERROR .NE. 0) GO TO 500
      CALL PI21Z0 (COORD,CRDIO)
  500 IFLG = IERROR
      IF (IFLG .NE. 0) RETURN
      IF (IOSYS .EQ. 0) GO TO 920
      COORD(1) = CRDIO(1)
      COORD(2) = CRDIO(2)
  520 IF (IOZONE.GT.60 .OR. IOSYS.EQ.1) GO TO 540
      IFLG = 1008
      RETURN
!
! FORWARD TRANSFORMATION.
  540 IPROJ=IOSYS
      IF (IOSYS.GT.2) CALL SPHDZ0(IOSPH,TPARIO)
      GO TO (610,620,630,640,650,660,670,680,690,700,710,720,730,740,750,760,770,780,790,800,810),IOSYS
  610 TPARIO(1) = COORD(1)
      TPARIO(2) = COORD(2)
!      CALL SPHDZ0(IOSPH,DUMMY)
!      TPARIO(14) = DUMMY(1)
!      TPARIO(15) = DUMMY(2)
      CALL IS1AZ0(IOZONE,TPARIO)
      IF (IERROR .NE. 0) GO TO 900
      CALL PF1AZ0 (COORD,CRDIO)
      GO TO 900
  620 CALL IS02Z0 (IOZONE,TPARIO)
      IF (IERROR .NE. 0) GO TO 900
      CALL PF02Z0 (COORD,CRDIO)
      GO TO 900
  630 CALL IS03Z0 (IOZONE,TPARIO)
      IF (IERROR .NE. 0) GO TO 900
      CALL PF03Z0 (COORD,CRDIO)
      GO TO 900
  640 CALL IS04Z0 (IOZONE,TPARIO)
      IF (IERROR .NE. 0) GO TO 900
      CALL PF04Z0 (COORD,CRDIO)
      GO TO 900
  650 CALL IS05Z0 (IOZONE,TPARIO)
      IF (IERROR .NE. 0) GO TO 900
      CALL PF05Z0 (COORD,CRDIO)
      GO TO 900
  660 CALL IS06Z0 (IOZONE,TPARIO)
      IF (IERROR .NE. 0) GO TO 900
      CALL PF06Z0 (COORD,CRDIO)
      GO TO 900
  670 CALL IS07Z0 (IOZONE,TPARIO)
      IF (IERROR .NE. 0) GO TO 900
      CALL PF07Z0 (COORD,CRDIO)
      GO TO 900
  680 CALL IS08Z0 (IOZONE,TPARIO)
      IF (IERROR .NE. 0) GO TO 900
      CALL PF08Z0 (COORD,CRDIO)
      GO TO 900
  690 CALL IS09Z0 (IOZONE,TPARIO)
      IF (IERROR .NE. 0) GO TO 900
      CALL PF09Z0 (COORD,CRDIO)
      GO TO 900
  700 CALL IS10Z0 (IOZONE,TPARIO)
      IF (IERROR .NE. 0) GO TO 900
      CALL PF10Z0 (COORD,CRDIO)
      GO TO 900
  710 CALL IS11Z0 (IOZONE,TPARIO)
      IF (IERROR .NE. 0) GO TO 900
      CALL PF11Z0 (COORD,CRDIO)
      GO TO 900
  720 CALL IS12Z0 (IOZONE,TPARIO)
      IF (IERROR .NE. 0) GO TO 900
      CALL PF12Z0 (COORD,CRDIO)
      GO TO 900
  730 CALL IS13Z0 (IOZONE,TPARIO)
      IF (IERROR .NE. 0) GO TO 900
      CALL PF13Z0 (COORD,CRDIO)
      GO TO 900
  740 CALL IS14Z0 (IOZONE,TPARIO)
      IF (IERROR .NE. 0) GO TO 900
      CALL PF14Z0 (COORD,CRDIO)
      GO TO 900
  750 CALL IS15Z0 (IOZONE,TPARIO)
      IF (IERROR .NE. 0) GO TO 900
      CALL PF15Z0 (COORD,CRDIO)
      GO TO 900
  760 CALL IS16Z0 (IOZONE,TPARIO)
      IF (IERROR .NE. 0) GO TO 900
      CALL PF16Z0 (COORD,CRDIO)
      GO TO 900
  770 CALL IS17Z0 (IOZONE,TPARIO)
      IF (IERROR .NE. 0) GO TO 900
      CALL PF17Z0 (COORD,CRDIO)
      GO TO 900
  780 CALL IS18Z0 (IOZONE,TPARIO)
      IF (IERROR .NE. 0) GO TO 900
      CALL PF18Z0 (COORD,CRDIO)
      GO TO 900
  790 CALL IS19Z0 (IOZONE,TPARIO)
      IF (IERROR .NE. 0) GO TO 900
      CALL PF19Z0 (COORD,CRDIO)
      GO TO 900
  800 CALL IS20Z0 (IOZONE,TPARIO)
      IF (IERROR .NE. 0) GO TO 900
      CALL PF20Z0 (COORD,CRDIO)
      GO TO 900
  810 CALL IS21Z0 (IOZONE,TPARIO)
      IF (IERROR .NE. 0) GO TO 900
      CALL PF21Z0 (COORD,CRDIO)
  900 IFLG = IERROR
  920 CONTINUE
      CRDIO(1) = FACTOR * CRDIO(1)
      CRDIO(2) = FACTOR * CRDIO(2)
      RETURN
      END

!                   GTRNZ0
! **********************************************************************
! ** U.S.G.S. GENERAL MAP PROJECTION PACKAGE ...... DR. A. A. ELASSAL **
! ** MODULE I                VERSION 1.0.0            NOVEMBER 1,1980 **
! **********************************************************************
      SUBROUTINE GTRNZ0 (CRDIN,INSYS,INZONE,TPARIN,INUNIT,CRDIO,IOSYS,IOZONE,TPARIO,IOUNIT,IFLG)
!
! GENERAL PROGRAM FOR TRANSFORMATION BETWEEN VARIOUS REFERENCE SYSTEMS
!
! INPUT ***************************************************************
! CRDIN  : COORDINATES IN INPUT SYSTEM (2 DP WORDS ARRAY).
! INSYS  : CODE NUMBER OF INPUT COORDINATE SYSTEM (INTEGER).
!            =  0 , GEOGRAPHIC
!            =  1 , U T M
!            =  2 , STATE PLANE
!            =  3 , ALBERS CONICAL EQUAL-AREA
!            =  4 , LAMBERT CONFORMAL CONIC
!            =  5 , MERCATOR
!            =  6 , POLAR STEREOGRAPHIC
!            =  7 , POLYCONIC
!            =  8 , EQUIDISTANT CONIC
!            =  9 , TRANSVERSE MERCATOR
!            = 10 , STEREOGRAPHIC
!            = 11 , LAMBERT AZIMUTHAL EQUAL-AREA
!            = 12 , AZIMUTHAL EQUIDISTANT
!            = 13 , GNOMONIC
!            = 14 , ORTHOGRAPHIC
!            = 15 , GENERAL VERTICAL NEAR-SIDE PERSPECTIVE
!            = 16 , SINUSOIDAL
!            = 17 , EQUIRECTANGULAR (PLATE CARREE)
!            = 18 , MILLER CYLINDRICAL
!            = 19 , VAN DER GRINTEN I
!            = 20 , OBLIQUE MERCATOR (HOTINE)
!            = 21 , SPACE OBLIQUE MERCATOR
!
! INZONE : CODE NUMBER OF INPUT COORDINATE ZONE (INTEGER).
! TPARIN : PARAMETERS OF INPUT REFERENCE SYSTEM (15 DP WORDS ARRAY).
! INUNIT : CODE NUMBER OF UNITS OF MEASURE FOR INPUT COORDS (INTEGER).
!            = 0 , RADIANS.
!            = 1 , FEET.
!            = 2 , METERS.
!            = 3 , SECONDS OF ARC.
!            = 4 , DEGREES OF ARC.
! IOSYS  : CODE NUMBER OF OUTPUT COORDINATE SYSTEM (INTEGER).
! IOZONE : CODE NUMBER OF OUTPUT COORDINATE ZONE (INTEGER).
! TPARIO : PARAMETERS OF OUTPUT REFERENCE SYSTEM (15 DP WORDS ARRAY).
! IOUNIT : CODE NUMBER OF UNITS OF MEASURE FOR OUTPUT COORDS (INTEGER).
!
! OUTPUT **************************************************************
! CRDIO  : COORDINATES IN OUTPUT REFERENCE SYSTEM (2 DP WORDS ARRAY).
! IFLG   : RETURN FLAG (INTEGER).
!            = 0 , SUCCESSFUL TRANSFORMATION.
!            = 1 , ILLEGAL INPUT SYSTEM CODE.
!            = 2 , ILLEGAL OUTPUT SYSTEM CODE.
!            = 3 , ILLEGAL INPUT UNIT CODE.
!            = 4 , ILLEGAL OUTPUT UNIT CODE.
!            = 5 , INCONSISTANT UNIT AND SYSTEM CODES FOR INPUT.
!            = 6 , INCONSISTANT UNIT AND SYSTEM CODES FOR OUTPUT.
!            = 7 , ILLEGAL INPUT ZONE CODE.
!            = 8 , ILLEGAL OUTPUT ZONE CODE.
!      OTHERWISE , ERROR CODE FROM PROJECTION COMPUTATIONAL MODULE.
!
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER*4 SYSUNT(22)
      COMMON /ERRMZ0/ IERROR
      DIMENSION CRDIN(1),CRDIO(1),TPARIN(1),TPARIO(1),COORD(2)
      DATA SYSUNT / 0 , 2 , 1 , 19*2 /
      DATA MAXUNT,MAXSYS / 4 , 21 /
!
! CHECK VALIDITY OF CODES FOR UNITS OF MEASURE AND REFERENCE SYSTEMS.
!
      IF (INSYS.GE.0 .AND. INSYS.LE.MAXSYS) GO TO 020
      IFLG = 1
      RETURN
  020 IF (IOSYS.GE.0 .AND. IOSYS.LE.MAXSYS) GO TO 040
      IFLG = 2
      RETURN
  040 IF (INUNIT.GE.0 .AND. INUNIT.LE.MAXUNT) GO TO 060
      IFLG = 3
      RETURN
  060 IF (IOUNIT.GE.0 .AND. IOUNIT.LE.MAXUNT) GO TO 080
      IFLG = 4
      RETURN
!
! CHECK CONSISTANCY BETEEN UNITS OF MEASURE AND REFERENCE SYSTEM.
!
  080 IUNIT = SYSUNT(INSYS + 1)
      CALL UNTFZ0 (INUNIT,IUNIT,FACTOR,IFLG)
      IF (IFLG .EQ. 0) GO TO 100
      IFLG = 5
      RETURN
  100 COORD(1) = FACTOR * CRDIN(1)
      COORD(2) = FACTOR * CRDIN(2)
      IUNIT = SYSUNT(IOSYS + 1)
      CALL UNTFZ0 (IUNIT,IOUNIT,FACTOR,IFLG)
      IF (IFLG .EQ. 0) GO TO 120
      IFLG = 6
      RETURN
  120 IF (INSYS.NE.IOSYS.OR.INZONE.NE.IOZONE.OR.INZONE.LE.0) GO TO 140
      CRDIO(1) = FACTOR * COORD(1)
      CRDIO(2) = FACTOR * COORD(2)
      RETURN
!
! COMPUTE TRANSFORMED COORDINATES AND ADJUST THEIR UNITS.
!
  140 IF (INSYS .EQ. 0) GO TO 520
      IF (INZONE.GT.60 .OR. INSYS.EQ.1) GO TO 200
      IFLG = 7
      RETURN
!
! INVERSE TRANSFORMATION.
!
  200 GO TO (210,220,230,240,250,260,270,280,290,300,310,320,330,340,350,360,370,380,390,400,410),INSYS
  210 CALL IS01Z0 (INZONE,TPARIN)
      IF (IERROR .NE. 0) GO TO 500
      CALL PI01Z0 (COORD,CRDIO)
      GO TO 500
  220 CALL IS02Z0 (INZONE,TPARIN)
      IF (IERROR .NE. 0) GO TO 500
      CALL PI02Z0 (COORD,CRDIO)
      GO TO 500
  230 CALL IS03Z0 (INZONE,TPARIN)
      IF (IERROR .NE. 0) GO TO 500
      CALL PI03Z0 (COORD,CRDIO)
      GO TO 500
  240 CALL IS04Z0 (INZONE,TPARIN)
      IF (IERROR .NE. 0) GO TO 500
      CALL PI04Z0 (COORD,CRDIO)
      GO TO 500
  250 CALL IS05Z0 (INZONE,TPARIN)
      IF (IERROR .NE. 0) GO TO 500
      CALL PI05Z0 (COORD,CRDIO)
      GO TO 500
  260 CALL IS06Z0 (INZONE,TPARIN)
      IF (IERROR .NE. 0) GO TO 500
      CALL PI06Z0 (COORD,CRDIO)
      GO TO 500
  270 CALL IS07Z0 (INZONE,TPARIN)
      IF (IERROR .NE. 0) GO TO 500
      CALL PI07Z0 (COORD,CRDIO)
      GO TO 500
  280 CALL IS08Z0 (INZONE,TPARIN)
      IF (IERROR .NE. 0) GO TO 500
      CALL PI08Z0 (COORD,CRDIO)
      GO TO 500
  290 continue !CALL IS09Z0 (INZONE,TPARIN) !Alex
      IF (IERROR .NE. 0) GO TO 500
      CALL PI09Z0 (COORD,CRDIO)
      GO TO 500
  300 CALL IS10Z0 (INZONE,TPARIN)
      IF (IERROR .NE. 0) GO TO 500
      CALL PI10Z0 (COORD,CRDIO)
      GO TO 500
  310 CALL IS11Z0 (INZONE,TPARIN)
      IF (IERROR .NE. 0) GO TO 500
      CALL PI11Z0 (COORD,CRDIO)
      GO TO 500
  320 CALL IS12Z0 (INZONE,TPARIN)
      IF (IERROR .NE. 0) GO TO 500
      CALL PI12Z0 (COORD,CRDIO)
      GO TO 500
  330 CALL IS13Z0 (INZONE,TPARIN)
      IF (IERROR .NE. 0) GO TO 500
      CALL PI13Z0 (COORD,CRDIO)
      GO TO 500
  340 CALL IS14Z0 (INZONE,TPARIN)
      IF (IERROR .NE. 0) GO TO 500
      CALL PI14Z0 (COORD,CRDIO)
      GO TO 500
  350 CALL IS15Z0 (INZONE,TPARIN)
      IF (IERROR .NE. 0) GO TO 500
      CALL PI15Z0 (COORD,CRDIO)
      GO TO 500
  360 CALL IS16Z0 (INZONE,TPARIN)
      IF (IERROR .NE. 0) GO TO 500
      CALL PI16Z0 (COORD,CRDIO)
      GO TO 500
  370 CALL IS17Z0 (INZONE,TPARIN)
      IF (IERROR .NE. 0) GO TO 500
      CALL PI17Z0 (COORD,CRDIO)
      GO TO 500
  380 CALL IS18Z0 (INZONE,TPARIN)
      IF (IERROR .NE. 0) GO TO 500
      CALL PI18Z0 (COORD,CRDIO)
      GO TO 500
  390 CALL IS19Z0 (INZONE,TPARIN)
      IF (IERROR .NE. 0) GO TO 500
      CALL PI19Z0 (COORD,CRDIO)
      GO TO 500
  400 CALL IS20Z0 (INZONE,TPARIN)
      IF (IERROR .NE. 0) GO TO 500
      CALL PI20Z0 (COORD,CRDIO)
      GO TO 500
  410 CALL IS21Z0 (INZONE,TPARIN)
      IF (IERROR .NE. 0) GO TO 500
      CALL PI21Z0 (COORD,CRDIO)
  500 IFLG = IERROR
      IF (IFLG .NE. 0) RETURN
      IF (IOSYS .EQ. 0) GO TO 920
      COORD(1) = CRDIO(1)
      COORD(2) = CRDIO(2)
  520 IF (IOZONE.GT.60 .OR. IOSYS.EQ.1) GO TO 540
      IFLG = 8
      RETURN
!
! FORWARD TRANSFORMATION.
!
  540 GO TO (610,620,630,640,650,660,670,680,690,700,710,720,730,740,750,760,770,780,790,800,810),IOSYS
  610 CALL IS01Z0 (IOZONE,TPARIO)
      IF (IERROR .NE. 0) GO TO 900
      CALL PF01Z0 (COORD,CRDIO)
      GO TO 900
  620 CALL IS02Z0 (IOZONE,TPARIO)
      IF (IERROR .NE. 0) GO TO 900
      CALL PF02Z0 (COORD,CRDIO)
      GO TO 900
  630 CALL IS03Z0 (IOZONE,TPARIO)
      IF (IERROR .NE. 0) GO TO 900
      CALL PF03Z0 (COORD,CRDIO)
      GO TO 900
  640 CALL IS04Z0 (IOZONE,TPARIO)
      IF (IERROR .NE. 0) GO TO 900
      CALL PF04Z0 (COORD,CRDIO)
      GO TO 900
  650 CALL IS05Z0 (IOZONE,TPARIO)
      IF (IERROR .NE. 0) GO TO 900
      CALL PF05Z0 (COORD,CRDIO)
      GO TO 900
  660 CALL IS06Z0 (IOZONE,TPARIO)
      IF (IERROR .NE. 0) GO TO 900
      CALL PF06Z0 (COORD,CRDIO)
      GO TO 900
  670 CALL IS07Z0 (IOZONE,TPARIO)
      IF (IERROR .NE. 0) GO TO 900
      CALL PF07Z0 (COORD,CRDIO)
      GO TO 900
  680 CALL IS08Z0 (IOZONE,TPARIO)
      IF (IERROR .NE. 0) GO TO 900
      CALL PF08Z0 (COORD,CRDIO)
      GO TO 900
  690 continue !CALL IS09Z0 (IOZONE,TPARIO) !Alex
      IF (IERROR .NE. 0) GO TO 900
      CALL PF09Z0 (COORD,CRDIO)
      GO TO 900
  700 CALL IS10Z0 (IOZONE,TPARIO)
      IF (IERROR .NE. 0) GO TO 900
      CALL PF10Z0 (COORD,CRDIO)
      GO TO 900
  710 CALL IS11Z0 (IOZONE,TPARIO)
      IF (IERROR .NE. 0) GO TO 900
      CALL PF11Z0 (COORD,CRDIO)
      GO TO 900
  720 CALL IS12Z0 (IOZONE,TPARIO)
      IF (IERROR .NE. 0) GO TO 900
      CALL PF12Z0 (COORD,CRDIO)
      GO TO 900
  730 CALL IS13Z0 (IOZONE,TPARIO)
      IF (IERROR .NE. 0) GO TO 900
      CALL PF13Z0 (COORD,CRDIO)
      GO TO 900
  740 CALL IS14Z0 (IOZONE,TPARIO)
      IF (IERROR .NE. 0) GO TO 900
      CALL PF14Z0 (COORD,CRDIO)
      GO TO 900
  750 CALL IS15Z0 (IOZONE,TPARIO)
      IF (IERROR .NE. 0) GO TO 900
      CALL PF15Z0 (COORD,CRDIO)
      GO TO 900
  760 CALL IS16Z0 (IOZONE,TPARIO)
      IF (IERROR .NE. 0) GO TO 900
      CALL PF16Z0 (COORD,CRDIO)
      GO TO 900
  770 CALL IS17Z0 (IOZONE,TPARIO)
      IF (IERROR .NE. 0) GO TO 900
      CALL PF17Z0 (COORD,CRDIO)
      GO TO 900
  780 CALL IS18Z0 (IOZONE,TPARIO)
      IF (IERROR .NE. 0) GO TO 900
      CALL PF18Z0 (COORD,CRDIO)
      GO TO 900
  790 CALL IS19Z0 (IOZONE,TPARIO)
      IF (IERROR .NE. 0) GO TO 900
      CALL PF19Z0 (COORD,CRDIO)
      GO TO 900
  800 CALL IS20Z0 (IOZONE,TPARIO)
      IF (IERROR .NE. 0) GO TO 900
      CALL PF20Z0 (COORD,CRDIO)
      GO TO 900
  810 CALL IS21Z0 (IOZONE,TPARIO)
      IF (IERROR .NE. 0) GO TO 900
      CALL PF21Z0 (COORD,CRDIO)
  900 IFLG = IERROR
  920 CONTINUE
      CRDIO(1) = FACTOR * CRDIO(1)
      CRDIO(2) = FACTOR * CRDIO(2)
      RETURN
!
      END

!                   MLFNZ0
! **********************************************************************
! ** U.S.G.S. GENERAL MAP PROJECTION PACKAGE ...... DR. A. A. ELASSAL **
! ** MODULE I                VERSION 1.0.0            NOVEMBER 1,1980 **
! **********************************************************************
      DOUBLE PRECISION FUNCTION MLFNZ0 (E0,E1,E2,PHI)
!
! FUNCTION TO COMPUTE CONSTANT (M).
!
      IMPLICIT REAL*8 (A-Z)
      DATA TWO,FOUR /2.0D0,4.0D0/
!
      MLFNZ0 = E0 * PHI - E1 * DSIN (TWO * PHI) + E2 * DSIN (FOUR * PHI)
!
      RETURN
      END

!                   MSFNZ0
! **********************************************************************
! ** U.S.G.S. GENERAL MAP PROJECTION PACKAGE ...... DR. A. A. ELASSAL **
! ** MODULE I                VERSION 1.0.0            NOVEMBER 1,1980 **
! **********************************************************************
      DOUBLE PRECISION FUNCTION MSFNZ0 (ECCENT,SINPHI,COSPHI)
!
! FUNCTION TO COMPUTE CONSTANT (SMALL M).
!
      IMPLICIT REAL*8 (A-Z)
      DATA ONE /1.0D0/
!
      CON = ECCENT * SINPHI
      MSFNZ0 = COSPHI / DSQRT (ONE - CON * CON)
!
      RETURN
      END

!                   PAKDZ0
! **********************************************************************
! ** U.S.G.S. GENERAL MAP PROJECTION PACKAGE ...... DR. A. A. ELASSAL **
! ** MODULE I                VERSION 1.0.0            NOVEMBER 1,1980 **
! **********************************************************************
      SUBROUTINE PAKDZ0 (PAK,DMS)
!
! SUBROUTINE TO CONVERT PACKED DMS TO UNPACKED DMS ANGLE.
!
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*4 REL
      INTEGER*4 DMS(1)
      EQUIVALENCE (REL , INT)
      DATA ZERO,CON1,CON2 /0.0D0,1000000.0D0,1000.0D0/
      DATA IBLNK,NEG /' ','-'/
!
      DMS(1) = IBLNK
      IF (PAK .LT. ZERO) DMS(1) = NEG
      CON = DABS (PAK)
      DMS(2) = CON / CON1
      CON = DMOD (CON , CON1)
      DMS(3) = CON / CON2
      REL = DMOD (CON , CON2)
      DMS(4) = INT
      RETURN
!
      END

!                   PAKRZ0
! **********************************************************************
! ** U.S.G.S. GENERAL MAP PROJECTION PACKAGE ...... DR. A. A. ELASSAL **
! ** MODULE I                VERSION 1.0.0            NOVEMBER 1,1980 **
! **********************************************************************
      DOUBLE PRECISION FUNCTION PAKRZ0 (ANG)
!
! FUNCTION TO CONVERT DMS PACKED ANGLE INTO RADIANS.
!
      IMPLICIT REAL*8 (A-H,O-Z)
      DATA SECRAD /0.4848136811095359D-5/
!
! CONVERT ANGLE TO SECONDS OF ARC
!
      SEC = PAKSZ0 (ANG)
!
! CONVERT ANGLE TO RADIANS.
!
      PAKRZ0 = SEC * SECRAD
!
      RETURN
      END

!                   PAKSZ0
! **********************************************************************
! ** U.S.G.S. GENERAL MAP PROJECTION PACKAGE ...... DR. A. A. ELASSAL **
! ** MODULE I                VERSION 1.0.0            NOVEMBER 1,1980 **
! **********************************************************************
      DOUBLE PRECISION FUNCTION PAKSZ0 (ANG)
!
! FUNCTION TO CONVERT DMS PACKED ANGLE INTO SECONDS OF ARC.
!
      IMPLICIT REAL*8 (A-H,M-Z)
      DIMENSION CODE(2)
      DATA CODE /1000000.0D0,1000.0D0/
      DATA ZERO,ONE /0.0D0,1.0D0/
      DATA C1,C2 /3600.0D0,60.0D0/
!
! SEPERATE DEGREE FIELD.
!
      FACTOR = ONE
      IF (ANG .LT. ZERO) FACTOR = - ONE
      SEC = DABS(ANG)
      TMP = CODE(1)
      I = SEC / TMP
      IF (I .GT. 360) GO TO 020
      DEG = I
!
! SEPERATE MINUTES FIELD.
!
      SEC = SEC - DEG * TMP
      TMP = CODE(2)
      I = SEC / TMP
      IF (I .GT. 60) GO TO 020
      MIN = I
!
! SEPERATE SECONDS FIELD.
!
      SEC = SEC - MIN * TMP
      IF (SEC .GT. C2) GO TO 020
      SEC = FACTOR * (DEG * C1 + MIN * C2 + SEC)
      GO TO 040
!
! ERROR DETECTED IN DMS FORM.
!
  20  continue
!
      paksz0 = -999
      return
  040 PAKSZ0 = SEC
!
      RETURN
      END

!                   PHI1Z0
! **********************************************************************
! ** U.S.G.S. GENERAL MAP PROJECTION PACKAGE ...... DR. A. A. ELASSAL **
! ** MODULE I                VERSION 1.0.0            NOVEMBER 1,1980 **
! **********************************************************************
      DOUBLE PRECISION FUNCTION PHI1Z0 (ECCENT,QS)
!
! FUNCTION TO COMPUTE LATITUDE ANGLE (PHI-1).
!
      IMPLICIT REAL*8 (A-Z)
      INTEGER*4 IERROR,IPEMSG,IPPARM
      INTEGER*4 II,NIT
      COMMON /ERRMZ0/ IERROR
      COMMON /PRINZ0/ IPEMSG,IPPARM
      DATA HALF,ONE /0.5D0,1.0D0/
      DATA EPSLN,TOL,NIT /1.0D-7,1.0D-10,15/
!
      PHI1Z0 =  DASIN (HALF * QS)
      IF (ECCENT .LT. EPSLN) RETURN
!
      ECCNTS = ECCENT * ECCENT
      PHI = PHI1Z0
      DO 020 II = 1,NIT
      SINPI = DSIN (PHI)
      COSPI = DCOS (PHI)
      CON = ECCENT * SINPI
      COM = ONE - CON * CON
      DPHI = HALF * COM * COM / COSPI * (QS / (ONE - ECCNTS) -  &
             SINPI / COM + HALF / ECCENT * DLOG ((ONE - CON) /  &
             (ONE + CON)))
      PHI = PHI + DPHI
      IF (DABS(DPHI) .GT. TOL) GO TO 020
      PHI1Z0 = PHI
      RETURN
  020 CONTINUE
!
      IF (IPEMSG .EQ. 0) PRINT 2000, NIT,ECCENT,QS
 2000 FORMAT (' ERROR PHI1Z0' /                                       &
              ' LATITUDE FAILED TO CONVERGE AFTER',I3,' ITERATIONS'/  &
              ' ECCENTRICITY =',D25.16,'   QS =',D25.16)
      IERROR = 001
      RETURN
!
      END

!                   PHI2Z0
! **********************************************************************
! ** U.S.G.S. GENERAL MAP PROJECTION PACKAGE ...... DR. A. A. ELASSAL **
! ** MODULE I                VERSION 1.0.0            NOVEMBER 1,1980 **
! **********************************************************************
      DOUBLE PRECISION FUNCTION PHI2Z0 (ECCENT,TS)
!
! FUNCTION TO COMPUTE LATITUDE ANGLE (PHI-2).
!
      IMPLICIT REAL*8 (A-Z)
      INTEGER*4 IERROR,IPEMSG,IPPARM
      INTEGER*4 II,NIT
      COMMON /ERRMZ0/ IERROR
      COMMON /PRINZ0/ IPEMSG,IPPARM
      DATA HALF,ONE,TWO /0.5D0,1.0D0,2.0D0/
      DATA TOL,NIT /1.0D-10,15/
      DATA HALFPI /1.5707963267948966D0/
!
      ECCNTH = HALF * ECCENT
      PHI = HALFPI - TWO * DATAN (TS)
      DO 020 II = 1,NIT
      SINPI = DSIN (PHI)
      CON = ECCENT * SINPI
      DPHI = HALFPI - TWO * DATAN (TS * ((ONE - CON) / (ONE + CON)) ** ECCNTH) - PHI
      PHI = PHI + DPHI
      IF (DABS(DPHI) .GT. TOL) GO TO 020
      PHI2Z0 = PHI
      RETURN
  020 CONTINUE
!
      IF (IPEMSG .EQ. 0) PRINT 2000, NIT,ECCENT,TS
 2000 FORMAT (' ERROR PHI2Z0' /                                       &
              ' LATITUDE FAILED TO CONVERGE AFTER',I3,' ITERATIONS'/  &
              ' ECCENTRICITY =',D25.16,'   TS =',D25.16)
      IERROR = 002
      RETURN
!
      END

!                   PHI3Z0
! **********************************************************************
! ** U.S.G.S. GENERAL MAP PROJECTION PACKAGE ...... DR. A. A. ELASSAL **
! ** MODULE I                VERSION 1.0.0            NOVEMBER 1,1980 **
! **********************************************************************
      DOUBLE PRECISION FUNCTION PHI3Z0 (ML,E0,E1,E2)
!
! FUNCTION TO COMPUTE LATITUDE ANGLE (PHI-3).
!
      IMPLICIT REAL*8 (A-Z)
      INTEGER*4 IERROR,IPEMSG,IPPARM
      INTEGER*4 II,NIT
      COMMON /ERRMZ0/ IERROR
      COMMON /PRINZ0/ IPEMSG,IPPARM
      DATA TWO,FOUR /2.0D0,4.0D0/
      DATA TOL,NIT /1.0D-10,15/
!
      PHI = ML
      DO 020 II = 1,NIT
      DPHI = (ML + E1 * DSIN (TWO * PHI) - E2 * DSIN (FOUR * PHI)) / E0 - PHI
      PHI = PHI + DPHI
      IF (DABS(DPHI) .GT. TOL) GO TO 020
      PHI3Z0 = PHI
      RETURN
  020 CONTINUE
!
      IF (IPEMSG .EQ. 0) PRINT 2000, NIT,ML,E0,E1,E2
 2000 FORMAT (' ERROR PHI3Z0' /                                       &
              ' LATITUDE FAILED TO CONVERGE AFTER',I3,' ITERATIONS'/  &
              ' ML =',D25.16,'   E0 =',D25.16/                        &
              ' E1 =',D25.16,'   E2 =',D25.16)
      IERROR = 003
      RETURN
!
      END

!                   PHI4Z0
! **********************************************************************
! ** U.S.G.S. GENERAL MAP PROJECTION PACKAGE ...... DR. A. A. ELASSAL **
! ** MODULE I                VERSION 1.0.0            NOVEMBER 1,1980 **
! **********************************************************************
      DOUBLE PRECISION FUNCTION PHI4Z0 (ECCNTS,E0,E1,E2,A,B,C)
!
! FUNCTION TO COMPUTE LATITUDE ANGLE (PHI-4).
!
      IMPLICIT REAL*8 (A-Z)
      INTEGER*4 IERROR,IPEMSG,IPPARM
      INTEGER*4 II,NIT
      COMMON /ERRMZ0/ IERROR
      COMMON /PRINZ0/ IPEMSG,IPPARM
      DATA ONE,TWO,FOUR /1.0D0,2.0D0,4.0D0/
      DATA TOL,NIT /1.0D-10,15/
!
      PHI = A
      DO 020 II = 1,NIT
      SINPHI = DSIN (PHI)
      TANPHI = DTAN (PHI)
      C = TANPHI * DSQRT (ONE - ECCNTS * SINPHI * SINPHI)
      SIN2PH = DSIN (TWO * PHI)
      ML = E0 * PHI - E1 * SIN2PH + E2 * DSIN (FOUR * PHI)
      MLP = E0 - TWO * E1 * DCOS (TWO * PHI) + FOUR * E2 * DCOS (FOUR * PHI)
      CON1 = TWO * ML + C * (ML * ML + B) - TWO * A * (C * ML + ONE)
      CON2 = ECCNTS * SIN2PH * (ML * ML + B - TWO * A * ML) / (TWO * C)
      CON3 = TWO * (A - ML) * (C * MLP - TWO / SIN2PH) - TWO * MLP
      DPHI = CON1 / (CON2 + CON3)
      PHI = PHI + DPHI
      IF (DABS(DPHI) .GT. TOL) GO TO 020
      PHI4Z0 = PHI
      RETURN
  020 CONTINUE
!
      IF (IPEMSG .EQ. 0) PRINT 2000, NIT,E0,E1,E2,A,B,C,ECCNTS
 2000 FORMAT (' ERROR PHI4Z0' /                                       &
              ' LATITUDE FAILED TO CONVERGE AFTER',I3,' ITERATIONS'/  &
              ' E0 =',D25.16,'   E1 =',D25.16/                        &
              ' E2 =',D25.16,'   A  =',D25.16/                        &
              ' B  =',D25.16,'   C  =',D25.16/                        &
              ' ECCENTRICITY SQUARE =',D25.16)
      IERROR = 004
      RETURN
!
      END

! ****                                                             *****
! ** U.S.G.S. GENERAL MAP PROJECTION PACKAGE ...... DR. A. A. ELASSAL **
! ** MODULE I                VERSION 1.0.0            NOVEMBER 1,1980 **
! ****                                                             *****
!                              *  U T M  *
! ****                                                             *****
!
      SUBROUTINE PJ1AZ0
!
!     MODIFIED APRIL 1981  FROM PJ01Z0  BY JOHN F. WAANANEN
!     SUBROUTINE PJ1AZ0 IS REQUIRED FOR PROGRAMS NO. L176 AND NO. L177
!
      IMPLICIT REAL*8 (A-Z)
      INTEGER*4 IERROR,IPEMSG,IPPARM
      INTEGER*4 SWITCH,IND,ZONE,INFILE,ANG
      integer*4 zzone
      DIMENSION BUFFL(15),DATA(1),GEOG(1),PROJ(1)
      COMMON /ERRMZ0/ IERROR
      COMMON /PRINZ0/ IPEMSG,IPPARM
!     COMMON /WORKZ0/ BUFF(15),ANG(4)
      COMMON /WORKZ0/ BUFF(15)
      COMMON /WK1AZ0/ ANG(4)
      COMMON /ELLPZ0/ AZ,EZ,ESZ,E0Z,E1Z,E2Z,E3Z
      real*4 rangs
      equivalence (rangs,ang(4))
      DATA HALFPI /1.5707963267948966D0/
      DATA ZERO /0.0D0/
      DATA SWITCH /0/
! ....                                                             .....
!       .  INITIALIZATION OF PROJECTION PARAMETERS (FILE INPUT)  .
! ....                                                             .....
      ENTRY IF1AZ0 (INFILE)
      IERROR = 0
      READ (INFILE,END=120) ZONE,BUFF
      IF (SWITCH.NE.0 .AND. SWITCH.EQ.ZONE) RETURN
  020 IF (ZONE .EQ. ZERO) GO TO 040
      IF (ZONE.GE.1.AND.ZONE.LE.60) GO TO 100
      IF (IPEMSG .EQ. 0) PRINT 2000,ZONE
 2000 FORMAT (' ILLEGAL ZONE NO : ',I10)
      IERROR = 011
      RETURN
  040 ZONE = PAKRZ0 (BUFF(1)) * 15.0D0 / HALFPI
      IND = 1
      IF (ZONE .LT. 0) IND = 0
      ZONE = MOD ((ZONE + 30) , 60) + IND
      IF (ZONE) 060,080,100
  060 ZONE = ZONE + 60
      GO TO 100
  080 ZONE = 1
  100 BUFFL(1) = BUFF(14)
      BUFFL(2) = BUFF(15)
      BUFFL(3) = 0.9996D0
      BUFFL(4) = ZERO
      LON0 = DFLOAT (6 * ZONE - 183) * HALFPI / 90.0D0
      CALL RADDZ0 (LON0,ANG)
      BUFFL(5) = DMSPZ0 (ANG)
      BUFFL(6) = ZERO
      BUFFL(7) = 500000.0D0
      BUFFL(8) = ZERO
      IF (BUFF(2) .LT. ZERO .OR. ZONE.LT.0) BUFFL(8) = 10000000.0D0
      IND = IPPARM
      IPPARM = 1
      SWITCH = 0
      CALL IS09Z0 (ZONE,BUFFL)
      IPPARM = IND
      IF (IERROR .NE. 0) GO TO 110
!
! LIST RESULTS OF PARAMETER INITIALIZATION.
      IF (IPPARM.EQ.0) then
      PRINT 12341
      PRINT 12342,zone
      PRINT 12343,buffl(1)
      PRINT 12344,buffl(2)
      PRINT 12345,buffl(3)
      PRINT 12346,ang(1),ang(2),ang(3),rangs
      PRINT 12347,buffl(7)
      PRINT 12348,buffl(8)
12341 format(' INITIALIZATION PARAMETERS (U T M PROJECTION)')
12342 format(' ZONE = ',I3)
12343 format(' SEMI-MAJOR AXIS OF ELLIPSOID = ',F12.2,' METERS')
12344 format(' ECCENTRICITY SQUARED         = ',f18.15)
12345 format(' SCALE FACTOR AT C. MERIDIAN  = ',F9.6)
12346 format(' LONGITUDE OF CENTRAL MERIDIAN= ',A1,2I3,F7.3)
12347 format(' FALSE EASTING                = ',f12.2,' METERS')
12348 format(' FALSE NORTHING               = ',f12.2,' METERS')
      end if
      SWITCH = ZONE
  110 RETURN
  120 IF (IPEMSG .EQ. 0) PRINT 2020
 2020 FORMAT (' MISSING PROJECTION PARAMETERS')
      IERROR = 012
      RETURN
! ....                                                             .....
!      .  INITIALIZATION OF PROJECTION PARAMETERS (ENTRY INPUT)  .
! ....                                                             .....
      ENTRY IS1AZ0 (ZZONE,DATA)
      zone = zzone
      IERROR = 0
      IF (SWITCH.NE.0.AND.SWITCH.EQ.ZONE.AND.DATA(14).EQ.BUFF(1))RETURN
      ZONE = IABS(ZONE)
      SWITCH = 0
      BUFF(1) = DATA(1)
      BUFF(2) = DATA(2)
      BUFF(14) = DATA(14)
      BUFF(15) = DATA(15)
      GO TO 020
! ....                                                             .....
!                      .  FORWARD TRANSFORMATION  .
! ....                                                             .....
      ENTRY PF1AZ0 (GEOG,PROJ)
      IERROR = 0
      IF (SWITCH .NE. 0) GO TO 140
      IF (IPEMSG .EQ. 0) PRINT 2020
      IERROR = 013
      RETURN
  140 CALL PF09Z0 (GEOG,PROJ)
      RETURN
! ....                                                             .....
!                      .  INVERSE TRANSFORMATION  .
! ....                                                             .....
      ENTRY PI1AZ0 (PROJ,GEOG)
      IERROR = 0
      IF (SWITCH .NE. 0) GO TO 160
      IF (IPEMSG .EQ. 0) PRINT 2020
      IERROR = 014
      RETURN
  160 CALL PI09Z0 (PROJ,GEOG)
      RETURN
      END

!                   PJ01Z0
! **********************************************************************
! ** U.S.G.S. GENERAL MAP PROJECTION PACKAGE ...... DR. A. A. ELASSAL **
! ** MODULE I                VERSION 1.0.0            NOVEMBER 1,1980 **
! **********************************************************************
!                              *  U T M  *
! **********************************************************************
!
      SUBROUTINE PJ01Z0
!
      IMPLICIT REAL*8 (A-Z)
      integer*4 zzone
      INTEGER*4 IERROR,IPEMSG,IPPARM
      INTEGER*4 SWITCH,IND,ZONE,INFILE,ANG
      COMMON /ERRMZ0/ IERROR
      COMMON /PRINZ0/ IPEMSG,IPPARM
!     COMMON /WORKZ0/ BUFF(15),BUFFL(15),ANG(4)
      COMMON /WORKZ0/ BUFF(15)
      COMMON /WK01Z0/ BUFFL(15),ANG(4)
      DIMENSION DATA(1),GEOG(1),PROJ(1)
      DATA HALFPI /1.5707963267948966D0/
      DATA ZERO /0.0D0/
      DATA SWITCH /0/
!
! ......................................................................
!       .  INITIALIZATION OF PROJECTION PARAMETERS (FILE INPUT)  .
! ......................................................................
!
      ENTRY IF01Z0 (INFILE,data)
!
      IERROR = 0
      READ (INFILE,END=120) ZONE,BUFF
      IF (SWITCH.NE.0 .AND. SWITCH.EQ.ZONE) RETURN
  020 IF (ZONE .EQ. ZERO) GO TO 040
      IF (ZONE.GE.1 .AND. ZONE.LE.60) GO TO 100
      IF (IPEMSG .EQ. 0) PRINT 2000, ZONE
 2000 FORMAT (' ERROR PJ01Z0'/           &
              ' ILLEGAL ZONE NO : ',I10)
      IERROR = 011
      RETURN
  040 ZONE = PAKRZ0 (BUFF(1)) * 15.0D0 / HALFPI
      IND = 1
      IF (ZONE .LT. 0) IND = 0
      ZONE = MOD ((ZONE + 30) , 60) + IND
      IF (ZONE) 060,080,100
  060 ZONE = ZONE + 60
      GO TO 100
  080 ZONE = 1
  100 BUFFL(1) = ZERO
      BUFFL(2) = ZERO
      BUFFL(3) = 0.9996D0
      BUFFL(4) = ZERO
      LON0 = DFLOAT (6 * ZONE - 183) * HALFPI / 90.0D0
      CALL RADDZ0 (LON0,ANG)
      BUFFL(5) = DMSPZ0 (ANG)
      BUFFL(6) = ZERO
      BUFFL(7) = 500000.0D0
      BUFFL(8) = ZERO
      IF (BUFF(2) .LT. ZERO) BUFFL(8) = 10000000.0D0
      IND = IPPARM
      IPPARM = 1
      SWITCH = 0
      CALL IS09Z0 (ZONE,BUFFL)
      IPPARM = IND
      IF (IERROR .NE. 0) GO TO 110
!
! LIST RESULTS OF PARAMETER INITIALIZATION.
!
      IF (IPPARM .EQ. 0) PRINT 2010, ZONE
 2010 FORMAT (' INITIALIZATION PARAMETERS (U T M PROJECTION)'/  &
              ' ZONE = ',I2)
      SWITCH = ZONE
  110 RETURN
  120 IF (IPEMSG .EQ. 0) PRINT 2020
 2020 FORMAT (' ERROR PJ01Z0'/                      &
              ' MISSING PROJECTION PARAMETERS')
      IERROR = 012
      RETURN
!
! ......................................................................
!      .  INITIALIZATION OF PROJECTION PARAMETERS (ENTRY INPUT)  .
! ......................................................................
!
      ENTRY IS01Z0 (ZZONE,DATA)
      zone = zzone
!
      IERROR = 0
      IF (SWITCH.NE.0 .AND. SWITCH.EQ.ZONE) RETURN
      BUFF(1) = DATA(1)
      BUFF(2) = DATA(2)
      GO TO 020
!
! ......................................................................
!                      .  FORWARD TRANSFORMATION  .
! ......................................................................
!
      ENTRY PF01Z0 (GEOG,PROJ)
!
      IERROR = 0
      IF (SWITCH .NE. 0) GO TO 140
      IF (IPEMSG .EQ. 0) PRINT 2020
      IERROR = 013
      RETURN
  140 CALL PF09Z0 (GEOG,PROJ)
      RETURN
!
! ......................................................................
!                      .  INVERSE TRANSFORMATION  .
! ......................................................................
!
      ENTRY PI01Z0 (PROJ,GEOG)
!
      IERROR = 0
      IF (SWITCH .NE. 0) GO TO 160
      IF (IPEMSG .EQ. 0) PRINT 2020
      IERROR = 014
      RETURN
  160 CALL PI09Z0 (PROJ,GEOG)
      RETURN
!
      END

!                   PJ02Z0
! **********************************************************************
! ** U.S.G.S. GENERAL MAP PROJECTION PACKAGE ...... DR. A. A. ELASSAL **
! ** MODULE I                VERSION 1.0.1               JUNE 10,1982 **
! **********************************************************************
!                           *  STATE PLANE  *
! **********************************************************************
!
      SUBROUTINE PJ02Z0
!                       CODE NUMBERS MODIFIED 6-10-82 BY J.F.WAANANEN
      IMPLICIT REAL*8(A-H,O-Z)
      integer*4 zzone
      INTEGER*4 SWITCH,ZONE
      COMMON /ERRMZ0/ IERROR
      COMMON /PRINZ0/ IPEMSG,IPPARM
      COMMON /WORKZ0/ BUFF(15)
      DIMENSION GEOG(1),PROJ(1),DATA(1)
      DIMENSION ITEM(131),ID(9,131),TABLE(11,122)
      DATA ITEM /0101,0102,5010, -1 ,0201,0202,0203,0301,0302,0401,0402,  &
                 0403,0404,0405,0406,0407,0501,0502,0503,0600,0700,0901,  &
                 0902,0903,1001,1002,5101,5102,5103,5104,5105,1101,1102,  &
                 1103,1201,1202,1301,1302,1401,1402,1501,1502,1601,1602,  &
                 1701,1702,1703,1801,1802,1900,2001,2002,2101,2102,2103,  &
                 2111,2112,2113,2201,2202,2203,2301,2302,2401,2402,2403,  &
                 2501,2502,2503,2601,2602,2701,2702,2703,2800,2900,3001,  &
                 3002,3003,3101,3102,3103,3104,3200,3301,3302,3401,3402,  &
                 3501,3502,3601,3602,3701,3702,3800,3901,3902,4001,4002,  &
                 4100,4201,4202,4203,4204,4205,4301,4302,4303,4400,4501,  &
                 4502,4601,4602,4701,4702,4801,4802,4803,4901,4902,4903,  &
                 4904,5001,5002,5003,5004,5005,5006,5007,5008,5009/
! ....................................................................
! ALABAMA                    EAST                                      T
       DATA ID(1,  1),ID(2,  1),ID(3,  1) /4HALAB,4HAMA ,4H    /
       DATA ID(4,  1),ID(5,  1),ID(6,  1) /4H    ,4HEAST,4H    /
       DATA ID(7,  1),ID(8,  1),ID(9,  1) /4H    ,4H    ,0/
! ALABAMA                    WEST                                      T
       DATA ID(1,  2),ID(2,  2),ID(3,  2) /4HALAB,4HAMA ,4H    /
       DATA ID(4,  2),ID(5,  2),ID(6,  2) /4H    ,4HWEST,4H    /
       DATA ID(7,  2),ID(8,  2),ID(9,  2) /4H    ,4H    ,0/
! ALASKA                     ZONE NO. 10
       DATA ID(1,  3),ID(2,  3),ID(3,  3) /4HALAS,4HKA  ,4H    /
       DATA ID(4,  3),ID(5,  3),ID(6,  3) /4H    ,4HZONE,4H NO./
       DATA ID(7,  3),ID(8,  3),ID(9,  3) /4H 10 ,4H    ,1/
! ALASKA                                                               T
       DATA ID(1,  4),ID(2,  4),ID(3,  4) /4HALAS,4HKA  ,4H    /
       DATA ID(4,  4),ID(5,  4),ID(6,  4) /4H    ,4H    ,4H    /
       DATA ID(7,  4),ID(8,  4),ID(9,  4) /4H    ,4H    ,0/
! ARIZONA                    EAST                                      T
       DATA ID(1,  5),ID(2,  5),ID(3,  5) /4HARIZ,4HONA ,4H    /
       DATA ID(4,  5),ID(5,  5),ID(6,  5) /4H    ,4HEAST,4H    /
       DATA ID(7,  5),ID(8,  5),ID(9,  5) /4H    ,4H    ,0/
! ARIZONA                    CENTRAL                                   T
       DATA ID(1,  6),ID(2,  6),ID(3,  6) /4HARIZ,4HONA ,4H    /
       DATA ID(4,  6),ID(5,  6),ID(6,  6) /4H    ,4HCENT,4HRAL /
       DATA ID(7,  6),ID(8,  6),ID(9,  6) /4H    ,4H    ,0/
! ARIZONA                    WEST                                      T
       DATA ID(1,  7),ID(2,  7),ID(3,  7) /4HARIZ,4HONA ,4H    /
       DATA ID(4,  7),ID(5,  7),ID(6,  7) /4H    ,4HWEST,4H    /
       DATA ID(7,  7),ID(8,  7),ID(9,  7) /4H    ,4H    ,0/
! ARKANSAS                   NORTH                                     L
       DATA ID(1,  8),ID(2,  8),ID(3,  8) /4HARKA,4HNSAS,4H    /
       DATA ID(4,  8),ID(5,  8),ID(6,  8) /4H    ,4HNORT,4HH   /
       DATA ID(7,  8),ID(8,  8),ID(9,  8) /4H    ,4H    ,1/
! ARKANSAS                   SOUTH                                     L
       DATA ID(1,  9),ID(2,  9),ID(3,  9) /4HARKA,4HNSAS,4H    /
       DATA ID(4,  9),ID(5,  9),ID(6,  9) /4H    ,4HSOUT,4HH   /
       DATA ID(7,  9),ID(8,  9),ID(9,  9) /4H    ,4H    ,1/
! CALIFORNIA                 I                                         L
       DATA ID(1, 10),ID(2, 10),ID(3, 10) /4HCALI,4HFORN,4HIA  /
       DATA ID(4, 10),ID(5, 10),ID(6, 10) /4H    ,4HI   ,4H    /
       DATA ID(7, 10),ID(8, 10),ID(9, 10) /4H    ,4H    ,1/
! CALIFORNIA                 II                                        L
       DATA ID(1, 11),ID(2, 11),ID(3, 11) /4HCALI,4HFORN,4HIA  /
       DATA ID(4, 11),ID(5, 11),ID(6, 11) /4H    ,4HII  ,4H    /
       DATA ID(7, 11),ID(8, 11),ID(9, 11) /4H    ,4H    ,1/
! CALIFORNIA                 III                                       L
       DATA ID(1, 12),ID(2, 12),ID(3, 12) /4HCALI,4HFORN,4HIA  /
       DATA ID(4, 12),ID(5, 12),ID(6, 12) /4H    ,4HIII ,4H    /
       DATA ID(7, 12),ID(8, 12),ID(9, 12) /4H    ,4H    ,1/
! CALIFORNIA                 IV                                        L
       DATA ID(1, 13),ID(2, 13),ID(3, 13) /4HCALI,4HFORN,4HIA  /
       DATA ID(4, 13),ID(5, 13),ID(6, 13) /4H    ,4HIV  ,4H    /
       DATA ID(7, 13),ID(8, 13),ID(9, 13) /4H    ,4H    ,1/
! CALIFORNIA                 V                                         L
       DATA ID(1, 14),ID(2, 14),ID(3, 14) /4HCALI,4HFORN,4HIA  /
       DATA ID(4, 14),ID(5, 14),ID(6, 14) /4H    ,4HV   ,4H    /
       DATA ID(7, 14),ID(8, 14),ID(9, 14) /4H    ,4H    ,1/
! CALIFORNIA                 VI                                        L
       DATA ID(1, 15),ID(2, 15),ID(3, 15) /4HCALI,4HFORN,4HIA  /
       DATA ID(4, 15),ID(5, 15),ID(6, 15) /4H    ,4HVI  ,4H    /
       DATA ID(7, 15),ID(8, 15),ID(9, 15) /4H    ,4H    ,1/
! CALIFORNIA                 VII                                       L
       DATA ID(1, 16),ID(2, 16),ID(3, 16) /4HCALI,4HFORN,4HIA  /
       DATA ID(4, 16),ID(5, 16),ID(6, 16) /4H    ,4HVII ,4H    /
       DATA ID(7, 16),ID(8, 16),ID(9, 16) /4H    ,4H    ,1/
! COLORADO                   NORTH                                     L
       DATA ID(1, 17),ID(2, 17),ID(3, 17) /4HCOLO,4HRADO,4H    /
       DATA ID(4, 17),ID(5, 17),ID(6, 17) /4H    ,4HNORT,4HH   /
       DATA ID(7, 17),ID(8, 17),ID(9, 17) /4H    ,4H    ,1/
! COLORADO                   CENTRAL                                   L
       DATA ID(1, 18),ID(2, 18),ID(3, 18) /4HCOLO,4HRADO,4H    /
       DATA ID(4, 18),ID(5, 18),ID(6, 18) /4H    ,4HCENT,4HRAL /
       DATA ID(7, 18),ID(8, 18),ID(9, 18) /4H    ,4H    ,1/
! COLORADO                   SOUTH                                     L
       DATA ID(1, 19),ID(2, 19),ID(3, 19) /4HCOLO,4HRADO,4H    /
       DATA ID(4, 19),ID(5, 19),ID(6, 19) /4H    ,4HSOUT,4HH   /
       DATA ID(7, 19),ID(8, 19),ID(9, 19) /4H    ,4H    ,1/
! CONNECTICUT                ---                                       L
       DATA ID(1, 20),ID(2, 20),ID(3, 20) /4HCONN,4HECTI,4HCUT /
       DATA ID(4, 20),ID(5, 20),ID(6, 20) /4H    ,4H--- ,4H    /
       DATA ID(7, 20),ID(8, 20),ID(9, 20) /4H    ,4H    ,1/
! DELAWARE                   ---                                       T
       DATA ID(1, 21),ID(2, 21),ID(3, 21) /4HDELA,4HWARE,4H    /
       DATA ID(4, 21),ID(5, 21),ID(6, 21) /4H    ,4H--- ,4H    /
       DATA ID(7, 21),ID(8, 21),ID(9, 21) /4H    ,4H    ,0/
! FLORIDA                    EAST                                      T
       DATA ID(1, 22),ID(2, 22),ID(3, 22) /4HFLOR,4HIDA ,4H    /
       DATA ID(4, 22),ID(5, 22),ID(6, 22) /4H    ,4HEAST,4H    /
       DATA ID(7, 22),ID(8, 22),ID(9, 22) /4H    ,4H    ,0/
! FLORIDA                    WEST                                      T
       DATA ID(1, 23),ID(2, 23),ID(3, 23) /4HFLOR,4HIDA ,4H    /
       DATA ID(4, 23),ID(5, 23),ID(6, 23) /4H    ,4HWEST,4H    /
       DATA ID(7, 23),ID(8, 23),ID(9, 23) /4H    ,4H    ,0/
! FLORIDA                    NORTH                                     L
       DATA ID(1, 24),ID(2, 24),ID(3, 24) /4HFLOR,4HIDA ,4H    /
       DATA ID(4, 24),ID(5, 24),ID(6, 24) /4H    ,4HNORT,4HH   /
       DATA ID(7, 24),ID(8, 24),ID(9, 24) /4H    ,4H    ,1/
! GEORGIA                    EAST                                      T
       DATA ID(1, 25),ID(2, 25),ID(3, 25) /4HGEOR,4HGIA ,4H    /
       DATA ID(4, 25),ID(5, 25),ID(6, 25) /4H    ,4HEAST,4H    /
       DATA ID(7, 25),ID(8, 25),ID(9, 25) /4H    ,4H    ,0/
! GEORGIA                    WEST                                      T
       DATA ID(1, 26),ID(2, 26),ID(3, 26) /4HGEOR,4HGIA ,4H    /
       DATA ID(4, 26),ID(5, 26),ID(6, 26) /4H    ,4HWEST,4H    /
       DATA ID(7, 26),ID(8, 26),ID(9, 26) /4H    ,4H    ,0/
! HAWAII                     1                                         T
       DATA ID(1, 27),ID(2, 27),ID(3, 27) /4HHAWA,4HII  ,4H    /
       DATA ID(4, 27),ID(5, 27),ID(6, 27) /4H    ,4H1   ,4H    /
       DATA ID(7, 27),ID(8, 27),ID(9, 27) /4H    ,4H    ,0/
! HAWAII                     2                                         T
       DATA ID(1, 28),ID(2, 28),ID(3, 28) /4HHAWA,4HII  ,4H    /
       DATA ID(4, 28),ID(5, 28),ID(6, 28) /4H    ,4H2   ,4H    /
       DATA ID(7, 28),ID(8, 28),ID(9, 28) /4H    ,4H    ,0/
! HAWAII                     3                                         T
       DATA ID(1, 29),ID(2, 29),ID(3, 29) /4HHAWA,4HII  ,4H    /
       DATA ID(4, 29),ID(5, 29),ID(6, 29) /4H    ,4H3   ,4H    /
       DATA ID(7, 29),ID(8, 29),ID(9, 29) /4H    ,4H    ,0/
! HAWAII                     4                                         T
       DATA ID(1, 30),ID(2, 30),ID(3, 30) /4HHAWA,4HII  ,4H    /
       DATA ID(4, 30),ID(5, 30),ID(6, 30) /4H    ,4H4   ,4H    /
       DATA ID(7, 30),ID(8, 30),ID(9, 30) /4H    ,4H    ,0/
! HAWAII                     5                                         T
       DATA ID(1, 31),ID(2, 31),ID(3, 31) /4HHAWA,4HII  ,4H    /
       DATA ID(4, 31),ID(5, 31),ID(6, 31) /4H    ,4H5   ,4H    /
       DATA ID(7, 31),ID(8, 31),ID(9, 31) /4H    ,4H    ,0/
! IDAHO                      EAST                                      T
       DATA ID(1, 32),ID(2, 32),ID(3, 32) /4HIDAH,4HO   ,4H    /
       DATA ID(4, 32),ID(5, 32),ID(6, 32) /4H    ,4HEAST,4H    /
       DATA ID(7, 32),ID(8, 32),ID(9, 32) /4H    ,4H    ,0/
! IDAHO                      CENTRAL                                   T
       DATA ID(1, 33),ID(2, 33),ID(3, 33) /4HIDAH,4HO   ,4H    /
       DATA ID(4, 33),ID(5, 33),ID(6, 33) /4H    ,4HCENT,4HRAL /
       DATA ID(7, 33),ID(8, 33),ID(9, 33) /4H    ,4H    ,0/
! IDAHO                      WEST                                      T
       DATA ID(1, 34),ID(2, 34),ID(3, 34) /4HIDAH,4HO   ,4H    /
       DATA ID(4, 34),ID(5, 34),ID(6, 34) /4H    ,4HWEST,4H    /
       DATA ID(7, 34),ID(8, 34),ID(9, 34) /4H    ,4H    ,0/
! ILLINOIS                   EAST                                      T
       DATA ID(1, 35),ID(2, 35),ID(3, 35) /4HILLI,4HNOIS,4H    /
       DATA ID(4, 35),ID(5, 35),ID(6, 35) /4H    ,4HEAST,4H    /
       DATA ID(7, 35),ID(8, 35),ID(9, 35) /4H    ,4H    ,0/
! ILLINOIS                   WEST                                      T
       DATA ID(1, 36),ID(2, 36),ID(3, 36) /4HILLI,4HNOIS,4H    /
       DATA ID(4, 36),ID(5, 36),ID(6, 36) /4H    ,4HWEST,4H    /
       DATA ID(7, 36),ID(8, 36),ID(9, 36) /4H    ,4H    ,0/
! INDIANA                    EAST                                      T
       DATA ID(1, 37),ID(2, 37),ID(3, 37) /4HINDI,4HANA ,4H    /
       DATA ID(4, 37),ID(5, 37),ID(6, 37) /4H    ,4HEAST,4H    /
       DATA ID(7, 37),ID(8, 37),ID(9, 37) /4H    ,4H    ,0/
! INDIANA                    WEST                                      T
       DATA ID(1, 38),ID(2, 38),ID(3, 38) /4HINDI,4HANA ,4H    /
       DATA ID(4, 38),ID(5, 38),ID(6, 38) /4H    ,4HWEST,4H    /
       DATA ID(7, 38),ID(8, 38),ID(9, 38) /4H    ,4H    ,0/
! IOWA                       NORTH                                     L
       DATA ID(1, 39),ID(2, 39),ID(3, 39) /4HIOWA,4H    ,4H    /
       DATA ID(4, 39),ID(5, 39),ID(6, 39) /4H    ,4HNORT,4HH   /
       DATA ID(7, 39),ID(8, 39),ID(9, 39) /4H    ,4H    ,1/
! IOWA                       SOUTH                                     L
       DATA ID(1, 40),ID(2, 40),ID(3, 40) /4HIOWA,4H    ,4H    /
       DATA ID(4, 40),ID(5, 40),ID(6, 40) /4H    ,4HSOUT,4HH   /
       DATA ID(7, 40),ID(8, 40),ID(9, 40) /4H    ,4H    ,1/
! KANSAS                     NORTH                                     L
       DATA ID(1, 41),ID(2, 41),ID(3, 41) /4HKANS,4HAS  ,4H    /
       DATA ID(4, 41),ID(5, 41),ID(6, 41) /4H    ,4HNORT,4HH   /
       DATA ID(7, 41),ID(8, 41),ID(9, 41) /4H    ,4H    ,1/
! KANSAS                     SOUTH                                     L
       DATA ID(1, 42),ID(2, 42),ID(3, 42) /4HKANS,4HAS  ,4H    /
       DATA ID(4, 42),ID(5, 42),ID(6, 42) /4H    ,4HSOUT,4HH   /
       DATA ID(7, 42),ID(8, 42),ID(9, 42) /4H    ,4H    ,1/
! KENTUCKY                   NORTH                                     L
       DATA ID(1, 43),ID(2, 43),ID(3, 43) /4HKENT,4HUCKY,4H    /
       DATA ID(4, 43),ID(5, 43),ID(6, 43) /4H    ,4HNORT,4HH   /
       DATA ID(7, 43),ID(8, 43),ID(9, 43) /4H    ,4H    ,1/
! KENTUCKY                   SOUTH                                     L
       DATA ID(1, 44),ID(2, 44),ID(3, 44) /4HKENT,4HUCKY,4H    /
       DATA ID(4, 44),ID(5, 44),ID(6, 44) /4H    ,4HSOUT,4HH   /
       DATA ID(7, 44),ID(8, 44),ID(9, 44) /4H    ,4H    ,1/
! LOUISIANA                  NORTH                                     L
       DATA ID(1, 45),ID(2, 45),ID(3, 45) /4HLOUI,4HSIAN,4HA   /
       DATA ID(4, 45),ID(5, 45),ID(6, 45) /4H    ,4HNORT,4HH   /
       DATA ID(7, 45),ID(8, 45),ID(9, 45) /4H    ,4H    ,1/
! LOUISIANA                  SOUTH                                     L
       DATA ID(1, 46),ID(2, 46),ID(3, 46) /4HLOUI,4HSIAN,4HA   /
       DATA ID(4, 46),ID(5, 46),ID(6, 46) /4H    ,4HSOUT,4HH   /
       DATA ID(7, 46),ID(8, 46),ID(9, 46) /4H    ,4H    ,1/
! LOUISIANA                  OFFSHORE                                  L
       DATA ID(1, 47),ID(2, 47),ID(3, 47) /4HLOUI,4HSIAN,4HA   /
       DATA ID(4, 47),ID(5, 47),ID(6, 47) /4H    ,4HOFFS,4HHORE/
       DATA ID(7, 47),ID(8, 47),ID(9, 47) /4H    ,4H    ,1/
! MAINE                      EAST                                      T
       DATA ID(1, 48),ID(2, 48),ID(3, 48) /4HMAIN,4HE   ,4H    /
       DATA ID(4, 48),ID(5, 48),ID(6, 48) /4H    ,4HEAST,4H    /
       DATA ID(7, 48),ID(8, 48),ID(9, 48) /4H    ,4H    ,0/
! MAINE                      WEST                                      T
       DATA ID(1, 49),ID(2, 49),ID(3, 49) /4HMAIN,4HE   ,4H    /
       DATA ID(4, 49),ID(5, 49),ID(6, 49) /4H    ,4HWEST,4H    /
       DATA ID(7, 49),ID(8, 49),ID(9, 49) /4H    ,4H    ,0/
! MARYLAND                   ---                                       L
       DATA ID(1, 50),ID(2, 50),ID(3, 50) /4HMARY,4HLAND,4H    /
       DATA ID(4, 50),ID(5, 50),ID(6, 50) /4H    ,4H--- ,4H    /
       DATA ID(7, 50),ID(8, 50),ID(9, 50) /4H    ,4H    ,1/
! MASSACHUSETTS              MAINLAND                                  L
       DATA ID(1, 51),ID(2, 51),ID(3, 51) /4HMASS,4HACHU,4HSETT/
       DATA ID(4, 51),ID(5, 51),ID(6, 51) /4HS   ,4HMAIN,4HLAND/
       DATA ID(7, 51),ID(8, 51),ID(9, 51) /4H    ,4H    ,1/
! MASSACHUSETTS              ISLAND                                    L
       DATA ID(1, 52),ID(2, 52),ID(3, 52) /4HMASS,4HACHU,4HSETT/
       DATA ID(4, 52),ID(5, 52),ID(6, 52) /4HS   ,4HISLA,4HND  /
       DATA ID(7, 52),ID(8, 52),ID(9, 52) /4H    ,4H    ,1/
! MICHIGAN                   EAST                                      T
       DATA ID(1, 53),ID(2, 53),ID(3, 53) /4HMICH,4HIGAN,4H    /
       DATA ID(4, 53),ID(5, 53),ID(6, 53) /4H    ,4HEAST,4H    /
       DATA ID(7, 53),ID(8, 53),ID(9, 53) /4H    ,4H    ,0/
! MICHIGAN                   CENTRAL                                   T
       DATA ID(1, 54),ID(2, 54),ID(3, 54) /4HMICH,4HIGAN,4H    /
       DATA ID(4, 54),ID(5, 54),ID(6, 54) /4H    ,4HCENT,4HRAL /
       DATA ID(7, 54),ID(8, 54),ID(9, 54) /4H    ,4H    ,0/
! MICHIGAN                   WEST                                      T
       DATA ID(1, 55),ID(2, 55),ID(3, 55) /4HMICH,4HIGAN,4H    /
       DATA ID(4, 55),ID(5, 55),ID(6, 55) /4H    ,4HWEST,4H    /
       DATA ID(7, 55),ID(8, 55),ID(9, 55) /4H    ,4H    ,0/
! MICHIGAN                   NORTH                                     L
       DATA ID(1, 56),ID(2, 56),ID(3, 56) /4HMICH,4HIGAN,4H    /
       DATA ID(4, 56),ID(5, 56),ID(6, 56) /4H    ,4HNORT,4HH   /
       DATA ID(7, 56),ID(8, 56),ID(9, 56) /4H    ,4H    ,1/
! MICHIGAN                   CENTRAL                                   L
       DATA ID(1, 57),ID(2, 57),ID(3, 57) /4HMICH,4HIGAN,4H    /
       DATA ID(4, 57),ID(5, 57),ID(6, 57) /4H    ,4HCENT,4HRAL /
       DATA ID(7, 57),ID(8, 57),ID(9, 57) /4H    ,4H    ,1/
! MICHIGAN                   SOUTH                                     L
       DATA ID(1, 58),ID(2, 58),ID(3, 58) /4HMICH,4HIGAN,4H    /
       DATA ID(4, 58),ID(5, 58),ID(6, 58) /4H    ,4HSOUT,4HH   /
       DATA ID(7, 58),ID(8, 58),ID(9, 58) /4H    ,4H    ,1/
! MINNESOTA                  NORTH                                     L
       DATA ID(1, 59),ID(2, 59),ID(3, 59) /4HMINN,4HESOT,4HA   /
       DATA ID(4, 59),ID(5, 59),ID(6, 59) /4H    ,4HNORT,4HH   /
       DATA ID(7, 59),ID(8, 59),ID(9, 59) /4H    ,4H    ,1/
! MINNESOTA                  CENTRAL                                   L
       DATA ID(1, 60),ID(2, 60),ID(3, 60) /4HMINN,4HESOT,4HA   /
       DATA ID(4, 60),ID(5, 60),ID(6, 60) /4H    ,4HCENT,4HRAL /
       DATA ID(7, 60),ID(8, 60),ID(9, 60) /4H    ,4H    ,1/
! MINNESOTA                  SOUTH                                     L
       DATA ID(1, 61),ID(2, 61),ID(3, 61) /4HMINN,4HESOT,4HA   /
       DATA ID(4, 61),ID(5, 61),ID(6, 61) /4H    ,4HSOUT,4HH   /
       DATA ID(7, 61),ID(8, 61),ID(9, 61) /4H    ,4H    ,1/
! MISSISSIPPI                EAST                                      T
       DATA ID(1, 62),ID(2, 62),ID(3, 62) /4HMISS,4HISSI,4HPPI /
       DATA ID(4, 62),ID(5, 62),ID(6, 62) /4H    ,4HEAST,4H    /
       DATA ID(7, 62),ID(8, 62),ID(9, 62) /4H    ,4H    ,0/
! MISSISSIPPI                WEST                                      T
       DATA ID(1, 63),ID(2, 63),ID(3, 63) /4HMISS,4HISSI,4HPPI /
       DATA ID(4, 63),ID(5, 63),ID(6, 63) /4H    ,4HWEST,4H    /
       DATA ID(7, 63),ID(8, 63),ID(9, 63) /4H    ,4H    ,0/
! MISSOURI                   EAST                                      T
       DATA ID(1, 64),ID(2, 64),ID(3, 64) /4HMISS,4HOURI,4H    /
       DATA ID(4, 64),ID(5, 64),ID(6, 64) /4H    ,4HEAST,4H    /
       DATA ID(7, 64),ID(8, 64),ID(9, 64) /4H    ,4H    ,0/
! MISSOURI                   CENTRAL                                   T
       DATA ID(1, 65),ID(2, 65),ID(3, 65) /4HMISS,4HOURI,4H    /
       DATA ID(4, 65),ID(5, 65),ID(6, 65) /4H    ,4HCENT,4HRAL /
       DATA ID(7, 65),ID(8, 65),ID(9, 65) /4H    ,4H    ,0/
! MISSOURI                   WEST                                      T
       DATA ID(1, 66),ID(2, 66),ID(3, 66) /4HMISS,4HOURI,4H    /
       DATA ID(4, 66),ID(5, 66),ID(6, 66) /4H    ,4HWEST,4H    /
       DATA ID(7, 66),ID(8, 66),ID(9, 66) /4H    ,4H    ,0/
! MONTANA                    NORTH                                     L
       DATA ID(1, 67),ID(2, 67),ID(3, 67) /4HMONT,4HANA ,4H    /
       DATA ID(4, 67),ID(5, 67),ID(6, 67) /4H    ,4HNORT,4HH   /
       DATA ID(7, 67),ID(8, 67),ID(9, 67) /4H    ,4H    ,1/
! MONTANA                    CENTRAL                                   L
       DATA ID(1, 68),ID(2, 68),ID(3, 68) /4HMONT,4HANA ,4H    /
       DATA ID(4, 68),ID(5, 68),ID(6, 68) /4H    ,4HCENT,4HRAL /
       DATA ID(7, 68),ID(8, 68),ID(9, 68) /4H    ,4H    ,1/
! MONTANA                    SOUTH                                     L
       DATA ID(1, 69),ID(2, 69),ID(3, 69) /4HMONT,4HANA ,4H    /
       DATA ID(4, 69),ID(5, 69),ID(6, 69) /4H    ,4HSOUT,4HH   /
       DATA ID(7, 69),ID(8, 69),ID(9, 69) /4H    ,4H    ,1/
! NEBRASKA                   NORTH                                     L
       DATA ID(1, 70),ID(2, 70),ID(3, 70) /4HNEBR,4HASKA,4H    /
       DATA ID(4, 70),ID(5, 70),ID(6, 70) /4H    ,4HNORT,4HH   /
       DATA ID(7, 70),ID(8, 70),ID(9, 70) /4H    ,4H    ,1/
! NEBRASKA                   SOUTH                                     L
       DATA ID(1, 71),ID(2, 71),ID(3, 71) /4HNEBR,4HASKA,4H    /
       DATA ID(4, 71),ID(5, 71),ID(6, 71) /4H    ,4HSOUT,4HH   /
       DATA ID(7, 71),ID(8, 71),ID(9, 71) /4H    ,4H    ,1/
! NEVADA                     EAST                                      T
       DATA ID(1, 72),ID(2, 72),ID(3, 72) /4HNEVA,4HDA  ,4H    /
       DATA ID(4, 72),ID(5, 72),ID(6, 72) /4H    ,4HEAST,4H    /
       DATA ID(7, 72),ID(8, 72),ID(9, 72) /4H    ,4H    ,0/
! NEVADA                     CENTRAL                                   T
       DATA ID(1, 73),ID(2, 73),ID(3, 73) /4HNEVA,4HDA  ,4H    /
       DATA ID(4, 73),ID(5, 73),ID(6, 73) /4H    ,4HCENT,4HRAL /
       DATA ID(7, 73),ID(8, 73),ID(9, 73) /4H    ,4H    ,0/
! NEVADA                     WEST                                      T
       DATA ID(1, 74),ID(2, 74),ID(3, 74) /4HNEVA,4HDA  ,4H    /
       DATA ID(4, 74),ID(5, 74),ID(6, 74) /4H    ,4HWEST,4H    /
       DATA ID(7, 74),ID(8, 74),ID(9, 74) /4H    ,4H    ,0/
! NEW HAMPSHIRE              ---                                       T
       DATA ID(1, 75),ID(2, 75),ID(3, 75) /4HNEW ,4HHAMP,4HSHIR/
       DATA ID(4, 75),ID(5, 75),ID(6, 75) /4HE   ,4H--- ,4H    /
       DATA ID(7, 75),ID(8, 75),ID(9, 75) /4H    ,4H    ,0/
! NEW JERSEY                 ---                                       T
       DATA ID(1, 76),ID(2, 76),ID(3, 76) /4HNEW ,4HJERS,4HEY  /
       DATA ID(4, 76),ID(5, 76),ID(6, 76) /4H    ,4H--- ,4H    /
       DATA ID(7, 76),ID(8, 76),ID(9, 76) /4H    ,4H    ,0/
! NEW MEXICO                 EAST                                      T
       DATA ID(1, 77),ID(2, 77),ID(3, 77) /4HNEW ,4HMEXI,4HCO  /
       DATA ID(4, 77),ID(5, 77),ID(6, 77) /4H    ,4HEAST,4H    /
       DATA ID(7, 77),ID(8, 77),ID(9, 77) /4H    ,4H    ,0/
! NEW MEXICO                 CENTRAL                                   T
       DATA ID(1, 78),ID(2, 78),ID(3, 78) /4HNEW ,4HMEXI,4HCO  /
       DATA ID(4, 78),ID(5, 78),ID(6, 78) /4H    ,4HCENT,4HRAL /
       DATA ID(7, 78),ID(8, 78),ID(9, 78) /4H    ,4H    ,0/
! NEW MEXICO                 WEST                                      T
       DATA ID(1, 79),ID(2, 79),ID(3, 79) /4HNEW ,4HMEXI,4HCO  /
       DATA ID(4, 79),ID(5, 79),ID(6, 79) /4H    ,4HWEST,4H    /
       DATA ID(7, 79),ID(8, 79),ID(9, 79) /4H    ,4H    ,0/
! NEW YORK                   EAST                                      T
       DATA ID(1, 80),ID(2, 80),ID(3, 80) /4HNEW ,4HYORK,4H    /
       DATA ID(4, 80),ID(5, 80),ID(6, 80) /4H    ,4HEAST,4H    /
       DATA ID(7, 80),ID(8, 80),ID(9, 80) /4H    ,4H    ,0/
! NEW YORK                   CENTRAL                                   T
       DATA ID(1, 81),ID(2, 81),ID(3, 81) /4HNEW ,4HYORK,4H    /
       DATA ID(4, 81),ID(5, 81),ID(6, 81) /4H    ,4HCENT,4HRAL /
       DATA ID(7, 81),ID(8, 81),ID(9, 81) /4H    ,4H    ,0/
! NEW YORK                   WEST                                      T
       DATA ID(1, 82),ID(2, 82),ID(3, 82) /4HNEW ,4HYORK,4H    /
       DATA ID(4, 82),ID(5, 82),ID(6, 82) /4H    ,4HWEST,4H    /
       DATA ID(7, 82),ID(8, 82),ID(9, 82) /4H    ,4H    ,0/
! NEW YORK                   LONG ISLAND                               L
       DATA ID(1, 83),ID(2, 83),ID(3, 83) /4HNEW ,4HYORK,4H    /
       DATA ID(4, 83),ID(5, 83),ID(6, 83) /4H    ,4HLONG,4H ISL/
       DATA ID(7, 83),ID(8, 83),ID(9, 83) /4HAND ,4H    ,1/
! NORTH CAROLINA             ---                                       L
       DATA ID(1, 84),ID(2, 84),ID(3, 84) /4HNORT,4HH CA,4HROLI/
       DATA ID(4, 84),ID(5, 84),ID(6, 84) /4HNA  ,4H--- ,4H    /
       DATA ID(7, 84),ID(8, 84),ID(9, 84) /4H    ,4H    ,1/
! NORTH DAKOTA               NORTH                                     L
       DATA ID(1, 85),ID(2, 85),ID(3, 85) /4HNORT,4HH DA,4HKOTA/
       DATA ID(4, 85),ID(5, 85),ID(6, 85) /4H    ,4HNORT,4HH   /
       DATA ID(7, 85),ID(8, 85),ID(9, 85) /4H    ,4H    ,1/
! NORTH DAKOTA               SOUTH                                     L
       DATA ID(1, 86),ID(2, 86),ID(3, 86) /4HNORT,4HH DA,4HKOTA/
       DATA ID(4, 86),ID(5, 86),ID(6, 86) /4H    ,4HSOUT,4HH   /
       DATA ID(7, 86),ID(8, 86),ID(9, 86) /4H    ,4H    ,1/
! OHIO                       NORTH                                     L
       DATA ID(1, 87),ID(2, 87),ID(3, 87) /4HOHIO,4H    ,4H    /
       DATA ID(4, 87),ID(5, 87),ID(6, 87) /4H    ,4HNORT,4HH   /
       DATA ID(7, 87),ID(8, 87),ID(9, 87) /4H    ,4H    ,1/
! OHIO                       SOUTH                                     L
       DATA ID(1, 88),ID(2, 88),ID(3, 88) /4HOHIO,4H    ,4H    /
       DATA ID(4, 88),ID(5, 88),ID(6, 88) /4H    ,4HSOUT,4HH   /
       DATA ID(7, 88),ID(8, 88),ID(9, 88) /4H    ,4H    ,1/
! OKLAHOMA                   NORTH                                     L
       DATA ID(1, 89),ID(2, 89),ID(3, 89) /4HOKLA,4HHOMA,4H    /
       DATA ID(4, 89),ID(5, 89),ID(6, 89) /4H    ,4HNORT,4HH   /
       DATA ID(7, 89),ID(8, 89),ID(9, 89) /4H    ,4H    ,1/
! OKLAHOMA                   SOUTH                                     L
       DATA ID(1, 90),ID(2, 90),ID(3, 90) /4HOKLA,4HHOMA,4H    /
       DATA ID(4, 90),ID(5, 90),ID(6, 90) /4H    ,4HSOUT,4HH   /
       DATA ID(7, 90),ID(8, 90),ID(9, 90) /4H    ,4H    ,1/
! OREGON                     NORTH                                     L
       DATA ID(1, 91),ID(2, 91),ID(3, 91) /4HOREG,4HON  ,4H    /
       DATA ID(4, 91),ID(5, 91),ID(6, 91) /4H    ,4HNORT,4HH   /
       DATA ID(7, 91),ID(8, 91),ID(9, 91) /4H    ,4H    ,1/
! OREGON                     SOUTH                                     L
       DATA ID(1, 92),ID(2, 92),ID(3, 92) /4HOREG,4HON  ,4H    /
       DATA ID(4, 92),ID(5, 92),ID(6, 92) /4H    ,4HSOUT,4HH   /
       DATA ID(7, 92),ID(8, 92),ID(9, 92) /4H    ,4H    ,1/
! PENNSYLVANIA               NORTH                                     L
       DATA ID(1, 93),ID(2, 93),ID(3, 93) /4HPENN,4HSYLV,4HANIA/
       DATA ID(4, 93),ID(5, 93),ID(6, 93) /4H    ,4HNORT,4HH   /
       DATA ID(7, 93),ID(8, 93),ID(9, 93) /4H    ,4H    ,1/
! PENNSYLVANIA               SOUTH                                     L
       DATA ID(1, 94),ID(2, 94),ID(3, 94) /4HPENN,4HSYLV,4HANIA/
       DATA ID(4, 94),ID(5, 94),ID(6, 94) /4H    ,4HSOUT,4HH   /
       DATA ID(7, 94),ID(8, 94),ID(9, 94) /4H    ,4H    ,1/
! RHODE ISLAND               ---                                       T
       DATA ID(1, 95),ID(2, 95),ID(3, 95) /4HRHOD,4HE IS,4HLAND/
       DATA ID(4, 95),ID(5, 95),ID(6, 95) /4H    ,4H--- ,4H    /
       DATA ID(7, 95),ID(8, 95),ID(9, 95) /4H    ,4H    ,0/
! SOUTH CAROLINA             NORTH                                     L
       DATA ID(1, 96),ID(2, 96),ID(3, 96) /4HSOUT,4HH CA,4HROLI/
       DATA ID(4, 96),ID(5, 96),ID(6, 96) /4HNA  ,4HNORT,4HH   /
       DATA ID(7, 96),ID(8, 96),ID(9, 96) /4H    ,4H    ,1/
! SOUTH CAROLINA             SOUTH                                     L
       DATA ID(1, 97),ID(2, 97),ID(3, 97) /4HSOUT,4HH CA,4HROLI/
       DATA ID(4, 97),ID(5, 97),ID(6, 97) /4HNA  ,4HSOUT,4HH   /
       DATA ID(7, 97),ID(8, 97),ID(9, 97) /4H    ,4H    ,1/
! SOUTH DAKOTA               NORTH                                     L
       DATA ID(1, 98),ID(2, 98),ID(3, 98) /4HSOUT,4HH DA,4HKOTA/
       DATA ID(4, 98),ID(5, 98),ID(6, 98) /4H    ,4HNORT,4HH   /
       DATA ID(7, 98),ID(8, 98),ID(9, 98) /4H    ,4H    ,1/
! SOUTH DAKOTA               SOUTH                                     L
       DATA ID(1, 99),ID(2, 99),ID(3, 99) /4HSOUT,4HH DA,4HKOTA/
       DATA ID(4, 99),ID(5, 99),ID(6, 99) /4H    ,4HSOUT,4HH   /
       DATA ID(7, 99),ID(8, 99),ID(9, 99) /4H    ,4H    ,1/
! TENNESSEE                  ---                                       L
       DATA ID(1,100),ID(2,100),ID(3,100) /4HTENN,4HESSE,4HE   /
       DATA ID(4,100),ID(5,100),ID(6,100) /4H    ,4H--- ,4H    /
       DATA ID(7,100),ID(8,100),ID(9,100) /4H    ,4H    ,1/
! TEXAS                      NORTH                                     L
       DATA ID(1,101),ID(2,101),ID(3,101) /4HTEXA,4HS   ,4H    /
       DATA ID(4,101),ID(5,101),ID(6,101) /4H    ,4HNORT,4HH   /
       DATA ID(7,101),ID(8,101),ID(9,101) /4H    ,4H    ,1/
! TEXAS                      NORTH CENTRAL                             L
       DATA ID(1,102),ID(2,102),ID(3,102) /4HTEXA,4HS   ,4H    /
       DATA ID(4,102),ID(5,102),ID(6,102) /4H    ,4HNORT,4HH CE/
       DATA ID(7,102),ID(8,102),ID(9,102) /4HNTRA,4HL   ,1/
! TEXAS                      CENTRAL                                   L
       DATA ID(1,103),ID(2,103),ID(3,103) /4HTEXA,4HS   ,4H    /
       DATA ID(4,103),ID(5,103),ID(6,103) /4H    ,4HCENT,4HRAL /
       DATA ID(7,103),ID(8,103),ID(9,103) /4H    ,4H    ,1/
! TEXAS                      SOUTH CENTRAL                             L
       DATA ID(1,104),ID(2,104),ID(3,104) /4HTEXA,4HS   ,4H    /
       DATA ID(4,104),ID(5,104),ID(6,104) /4H    ,4HSOUT,4HH CE/
       DATA ID(7,104),ID(8,104),ID(9,104) /4HNTRA,4HL   ,1/
! TEXAS                      SOUTH                                     L
       DATA ID(1,105),ID(2,105),ID(3,105) /4HTEXA,4HS   ,4H    /
       DATA ID(4,105),ID(5,105),ID(6,105) /4H    ,4HSOUT,4HH   /
       DATA ID(7,105),ID(8,105),ID(9,105) /4H    ,4H    ,1/
! UTAH                       NORTH                                     L
       DATA ID(1,106),ID(2,106),ID(3,106) /4HUTAH,4H    ,4H    /
       DATA ID(4,106),ID(5,106),ID(6,106) /4H    ,4HNORT,4HH   /
       DATA ID(7,106),ID(8,106),ID(9,106) /4H    ,4H    ,1/
! UTAH                       CENTRAL                                   L
       DATA ID(1,107),ID(2,107),ID(3,107) /4HUTAH,4H    ,4H    /
       DATA ID(4,107),ID(5,107),ID(6,107) /4H    ,4HCENT,4HRAL /
       DATA ID(7,107),ID(8,107),ID(9,107) /4H    ,4H    ,1/
! UTAH                       SOUTH                                     L
       DATA ID(1,108),ID(2,108),ID(3,108) /4HUTAH,4H    ,4H    /
       DATA ID(4,108),ID(5,108),ID(6,108) /4H    ,4HSOUT,4HH   /
       DATA ID(7,108),ID(8,108),ID(9,108) /4H    ,4H    ,1/
! VERMONT                    ---                                       T
       DATA ID(1,109),ID(2,109),ID(3,109) /4HVERM,4HONT ,4H    /
       DATA ID(4,109),ID(5,109),ID(6,109) /4H    ,4H--- ,4H    /
       DATA ID(7,109),ID(8,109),ID(9,109) /4H    ,4H    ,0/
! VIRGINIA                   NORTH                                     L
       DATA ID(1,110),ID(2,110),ID(3,110) /4HVIRG,4HINIA,4H    /
       DATA ID(4,110),ID(5,110),ID(6,110) /4H    ,4HNORT,4HH   /
       DATA ID(7,110),ID(8,110),ID(9,110) /4H    ,4H    ,1/
! VIRGINIA                   SOUTH                                     L
       DATA ID(1,111),ID(2,111),ID(3,111) /4HVIRG,4HINIA,4H    /
       DATA ID(4,111),ID(5,111),ID(6,111) /4H    ,4HSOUT,4HH   /
       DATA ID(7,111),ID(8,111),ID(9,111) /4H    ,4H    ,1/
! WASHINGTON                 NORTH                                     L
       DATA ID(1,112),ID(2,112),ID(3,112) /4HWASH,4HINGT,4HON  /
       DATA ID(4,112),ID(5,112),ID(6,112) /4H    ,4HNORT,4HH   /
       DATA ID(7,112),ID(8,112),ID(9,112) /4H    ,4H    ,1/
! WASHINGTON                 SOUTH                                     L
       DATA ID(1,113),ID(2,113),ID(3,113) /4HWASH,4HINGT,4HON  /
       DATA ID(4,113),ID(5,113),ID(6,113) /4H    ,4HSOUT,4HH   /
       DATA ID(7,113),ID(8,113),ID(9,113) /4H    ,4H    ,1/
! WEST VIRGINIA              NORTH                                     L
       DATA ID(1,114),ID(2,114),ID(3,114) /4HWEST,4H VIR,4HGINI/
       DATA ID(4,114),ID(5,114),ID(6,114) /4HA   ,4HNORT,4HH   /
       DATA ID(7,114),ID(8,114),ID(9,114) /4H    ,4H    ,1/
! WEST VIRGINIA              SOUTH                                     L
       DATA ID(1,115),ID(2,115),ID(3,115) /4HWEST,4H VIR,4HGINI/
       DATA ID(4,115),ID(5,115),ID(6,115) /4HA   ,4HSOUT,4HH   /
       DATA ID(7,115),ID(8,115),ID(9,115) /4H    ,4H    ,1/
! WISCONSIN                  NORTH                                     L
       DATA ID(1,116),ID(2,116),ID(3,116) /4HWISC,4HONSI,4HN   /
       DATA ID(4,116),ID(5,116),ID(6,116) /4H    ,4HNORT,4HH   /
       DATA ID(7,116),ID(8,116),ID(9,116) /4H    ,4H    ,1/
! WISCONSIN                  CENTRAL                                   L
       DATA ID(1,117),ID(2,117),ID(3,117) /4HWISC,4HONSI,4HN   /
       DATA ID(4,117),ID(5,117),ID(6,117) /4H    ,4HCENT,4HRAL /
       DATA ID(7,117),ID(8,117),ID(9,117) /4H    ,4H    ,1/
! WISCONSIN                  SOUTH                                     L
       DATA ID(1,118),ID(2,118),ID(3,118) /4HWISC,4HONSI,4HN   /
       DATA ID(4,118),ID(5,118),ID(6,118) /4H    ,4HSOUT,4HH   /
       DATA ID(7,118),ID(8,118),ID(9,118) /4H    ,4H    ,1/
! WYOMING                    EAST                                      T
       DATA ID(1,119),ID(2,119),ID(3,119) /4HWYOM,4HING ,4H    /
       DATA ID(4,119),ID(5,119),ID(6,119) /4H    ,4HEAST,4H    /
       DATA ID(7,119),ID(8,119),ID(9,119) /4H    ,4H    ,0/
! WYOMING                    EAST CENTRAL                              T
       DATA ID(1,120),ID(2,120),ID(3,120) /4HWYOM,4HING ,4H    /
       DATA ID(4,120),ID(5,120),ID(6,120) /4H    ,4HEAST,4H CEN/
       DATA ID(7,120),ID(8,120),ID(9,120) /4HTRAL,4H    ,0/
! WYOMING                    WEST CENTRAL                              T
       DATA ID(1,121),ID(2,121),ID(3,121) /4HWYOM,4HING ,4H    /
       DATA ID(4,121),ID(5,121),ID(6,121) /4H    ,4HWEST,4H CEN/
       DATA ID(7,121),ID(8,121),ID(9,121) /4HTRAL,4H    ,0/
! WYOMING                    WEST                                      T
       DATA ID(1,122),ID(2,122),ID(3,122) /4HWYOM,4HING ,4H    /
       DATA ID(4,122),ID(5,122),ID(6,122) /4H    ,4HWEST,4H    /
       DATA ID(7,122),ID(8,122),ID(9,122) /4H    ,4H    ,0/
! ALASKA                     ZONE NO. 1
       DATA ID(1,123),ID(2,123),ID(3,123) /4HALAS,4HKA  ,4H    /
       DATA ID(4,123),ID(5,123),ID(6,123) /4H    ,4HZONE,4H NO./
       DATA ID(7,123),ID(8,123),ID(9,123) /4H 1  ,4H    ,2/
! ALASKA                     ZONE NO. 2
       DATA ID(1,124),ID(2,124),ID(3,124) /4HALAS,4HKA  ,4H    /
       DATA ID(4,124),ID(5,124),ID(6,124) /4H    ,4HZONE,4H NO./
       DATA ID(7,124),ID(8,124),ID(9,124) /4H 2  ,4H    ,2/
! ALASKA                     ZONE NO. 3
       DATA ID(1,125),ID(2,125),ID(3,125) /4HALAS,4HKA  ,4H    /
       DATA ID(4,125),ID(5,125),ID(6,125) /4H    ,4HZONE,4H NO./
       DATA ID(7,125),ID(8,125),ID(9,125) /4H 3  ,4H    ,2/
! ALASKA                     ZONE NO. 4
       DATA ID(1,126),ID(2,126),ID(3,126) /4HALAS,4HKA  ,4H    /
       DATA ID(4,126),ID(5,126),ID(6,126) /4H    ,4HZONE,4H NO./
       DATA ID(7,126),ID(8,126),ID(9,126) /4H 4  ,4H    ,2/
! ALASKA                     ZONE NO. 5
       DATA ID(1,127),ID(2,127),ID(3,127) /4HALAS,4HKA  ,4H    /
       DATA ID(4,127),ID(5,127),ID(6,127) /4H    ,4HZONE,4H NO./
       DATA ID(7,127),ID(8,127),ID(9,127) /4H 5  ,4H    ,2/
! ALASKA                     ZONE NO. 6
       DATA ID(1,128),ID(2,128),ID(3,128) /4HALAS,4HKA  ,4H    /
       DATA ID(4,128),ID(5,128),ID(6,128) /4H    ,4HZONE,4H NO./
       DATA ID(7,128),ID(8,128),ID(9,128) /4H 6  ,4H    ,2/
! ALASKA                     ZONE NO. 7
       DATA ID(1,129),ID(2,129),ID(3,129) /4HALAS,4HKA  ,4H    /
       DATA ID(4,129),ID(5,129),ID(6,129) /4H    ,4HZONE,4H NO./
       DATA ID(7,129),ID(8,129),ID(9,129) /4H 7  ,4H    ,2/
! ALASKA                     ZONE NO. 8
       DATA ID(1,130),ID(2,130),ID(3,130) /4HALAS,4HKA  ,4H    /
       DATA ID(4,130),ID(5,130),ID(6,130) /4H    ,4HZONE,4H NO./
       DATA ID(7,130),ID(8,130),ID(9,130) /4H 8  ,4H    ,2/
! ALASKA                     ZONE NO. 9
       DATA ID(1,131),ID(2,131),ID(3,131) /4HALAS,4HKA  ,4H    /
       DATA ID(4,131),ID(5,131),ID(6,131) /4H    ,4HZONE,4H NO./
       DATA ID(7,131),ID(8,131),ID(9,131) /4H 9  ,4H    ,2/
! ....................................................................
       DATA TABLE(1 ,1  ) /500000.00     /
       DATA TABLE(2 ,1  ) /309000.00     /
       DATA TABLE(3 ,1  ) /1822.0        /
       DATA TABLE(4 ,1  ) /21.00903      /
       DATA TABLE(5 ,1  ) /.9999600000   /
       DATA TABLE(6 ,1  ) /.3817065      /
       DATA TABLE(1 ,2  ) /500000.00     /
       DATA TABLE(2 ,2  ) /315000.00     /
       DATA TABLE(3 ,2  ) /1792.0        /
       DATA TABLE(4 ,2  ) /25.53386      /
       DATA TABLE(5 ,2  ) /.9999333333   /
       DATA TABLE(6 ,2  ) /.3817477      /
       DATA TABLE(1 ,5  ) /500000.00     /
       DATA TABLE(2 ,5  ) /396600.00     /
       DATA TABLE(3 ,5  ) /1852.0        /
       DATA TABLE(4 ,5  ) /16.62358      /
       DATA TABLE(5 ,5  ) /.9999000000   /
       DATA TABLE(6 ,5  ) /.3816485      /
       DATA TABLE(1 ,6  ) /500000.00     /
       DATA TABLE(2 ,6  ) /402900.00     /
       DATA TABLE(3 ,6  ) /1852.0        /
       DATA TABLE(4 ,6  ) /16.62358      /
       DATA TABLE(5 ,6  ) /.9999000000   /
       DATA TABLE(6 ,6  ) /.3816485      /
       DATA TABLE(1 ,7  ) /500000.00     /
       DATA TABLE(2 ,7  ) /409500.00     /
       DATA TABLE(3 ,7  ) /1852.0        /
       DATA TABLE(4 ,7  ) /16.62358      /
       DATA TABLE(5 ,7  ) /.9999333333   /
       DATA TABLE(6 ,7  ) /.3815948      /
       DATA TABLE(1 ,21 ) /500000.00     /
       DATA TABLE(2 ,21 ) /271500.00     /
       DATA TABLE(3 ,21 ) /2271.0        /
       DATA TABLE(4 ,21 ) /30.53702      /
       DATA TABLE(5 ,21 ) /.9999950281   /
       DATA TABLE(6 ,21 ) /.3811454      /
       DATA TABLE(1 ,22 ) /500000.00     /
       DATA TABLE(2 ,22 ) /291600.00     /
       DATA TABLE(3 ,22 ) /1453.0        /
       DATA TABLE(4 ,22 ) /26.09287      /
       DATA TABLE(5 ,22 ) /.9999411765   /
       DATA TABLE(6 ,22 ) /.3821090      /
       DATA TABLE(1 ,23 ) /500000.00     /
       DATA TABLE(2 ,23 ) /295200.00     /
       DATA TABLE(3 ,23 ) /1453.0        /
       DATA TABLE(4 ,23 ) /26.09287      /
       DATA TABLE(5 ,23 ) /.9999411765   /
       DATA TABLE(6 ,23 ) /.3821090      /
       DATA TABLE(1 ,25 ) /500000.00     /
       DATA TABLE(2 ,25 ) /295800.00     /
       DATA TABLE(3 ,25 ) /1792.0        /
       DATA TABLE(4 ,25 ) /25.53386      /
       DATA TABLE(5 ,25 ) /.9999000000   /
       DATA TABLE(6 ,25 ) /.3817593      /
       DATA TABLE(1 ,26 ) /500000.00     /
       DATA TABLE(2 ,26 ) /303000.00     /
       DATA TABLE(3 ,26 ) /1792.0        /
       DATA TABLE(4 ,26 ) /25.53386      /
       DATA TABLE(5 ,26 ) /.9999000000   /
       DATA TABLE(6 ,26 ) /.3817593      /
       DATA TABLE(1 ,27 ) /500000.00     /
       DATA TABLE(2 ,27 ) /559800.00     /
       DATA TABLE(3 ,27 ) /1124.0        /
       DATA TABLE(4 ,27 ) /39.52714      /
       DATA TABLE(5 ,27 ) /.9999666667   /
       DATA TABLE(6 ,27 ) /.3826496      /
       DATA TABLE(1 ,28 ) /500000.00     /
       DATA TABLE(2 ,28 ) /564000.00     /
       DATA TABLE(3 ,28 ) /1214.0        /
       DATA TABLE(4 ,28 ) /18.21554      /
       DATA TABLE(5 ,28 ) /.9999666667   /
       DATA TABLE(6 ,28 ) /.3825762      /
       DATA TABLE(1 ,29 ) /500000.00     /
       DATA TABLE(2 ,29 ) /568800.00     /
       DATA TABLE(3 ,29 ) /1264.0        /
       DATA TABLE(4 ,29 ) /6.77497       /
       DATA TABLE(5 ,29 ) /.9999900000   /
       DATA TABLE(6 ,29 ) /.3825176      /
       DATA TABLE(1 ,30 ) /500000.00     /
       DATA TABLE(2 ,30 ) /574200.00     /
       DATA TABLE(3 ,30 ) /1303.0        /
       DATA TABLE(4 ,30 ) /57.83623      /
       DATA TABLE(5 ,30 ) /.9999900000   /
       DATA TABLE(6 ,30 ) /.3824812      /
       DATA TABLE(1 ,31 ) /500000.00     /
       DATA TABLE(2 ,31 ) /576600.00     /
       DATA TABLE(3 ,31 ) /1294.0        /
       DATA TABLE(4 ,31 ) /0.05280       /
       DATA TABLE(5 ,31 ) /.9999999999   /
       DATA TABLE(6 ,31 ) /.3824867      /
       DATA TABLE(1 ,32 ) /500000.00     /
       DATA TABLE(2 ,32 ) /403800.00     /
       DATA TABLE(3 ,32 ) /2491.0        /
       DATA TABLE(4 ,32 ) /18.35156      /
       DATA TABLE(5 ,32 ) /.9999473684   /
       DATA TABLE(6 ,32 ) /.3807624      /
       DATA TABLE(1 ,33 ) /500000.00     /
       DATA TABLE(2 ,33 ) /410400.00     /
       DATA TABLE(3 ,33 ) /2491.0        /
       DATA TABLE(4 ,33 ) /18.35156      /
       DATA TABLE(5 ,33 ) /.9999473684   /
       DATA TABLE(6 ,33 ) /.3807624      /
       DATA TABLE(1 ,34 ) /500000.00     /
       DATA TABLE(2 ,34 ) /416700.00     /
       DATA TABLE(3 ,34 ) /2491.0        /
       DATA TABLE(4 ,34 ) /18.35156      /
       DATA TABLE(5 ,34 ) /.9999333333   /
       DATA TABLE(6 ,34 ) /.3806227      /
       DATA TABLE(1 ,35 ) /500000.00     /
       DATA TABLE(2 ,35 ) /318000.00     /
       DATA TABLE(3 ,35 ) /2191.0        /
       DATA TABLE(4 ,35 ) /37.04639      /
       DATA TABLE(5 ,35 ) /.9999750000   /
       DATA TABLE(6 ,35 ) /.3811074      /
       DATA TABLE(1 ,36 ) /500000.00     /
       DATA TABLE(2 ,36 ) /324600.00     /
       DATA TABLE(3 ,36 ) /2191.0        /
       DATA TABLE(4 ,36 ) /37.04639      /
       DATA TABLE(5 ,36 ) /.9999411765   /
       DATA TABLE(6 ,36 ) /.3811332      /
       DATA TABLE(1 ,37 ) /500000.00     /
       DATA TABLE(2 ,37 ) /308400.00     /
       DATA TABLE(3 ,37 ) /2241.0        /
       DATA TABLE(4 ,37 ) /32.84965      /
       DATA TABLE(5 ,37 ) /.9999666667   /
       DATA TABLE(6 ,37 ) /.3811064      /
       DATA TABLE(1 ,38 ) /500000.00     /
       DATA TABLE(2 ,38 ) /313500.00     /
       DATA TABLE(3 ,38 ) /2241.0        /
       DATA TABLE(4 ,38 ) /32.84965      /
       DATA TABLE(5 ,38 ) /.9999666667   /
       DATA TABLE(6 ,38 ) /.3811064      /
       DATA TABLE(1 ,48 ) /500000.00     /
       DATA TABLE(2 ,48 ) /246600.00     /
       DATA TABLE(3 ,48 ) /2621.0        /
       DATA TABLE(4 ,48 ) /15.15187      /
       DATA TABLE(5 ,48 ) /.9999000000   /
       DATA TABLE(6 ,48 ) /.3806180      /
       DATA TABLE(1 ,49 ) /500000.00     /
       DATA TABLE(2 ,49 ) /252600.00     /
       DATA TABLE(3 ,49 ) /2561.0        /
       DATA TABLE(4 ,49 ) /16.25668      /
       DATA TABLE(5 ,49 ) /.9999666667   /
       DATA TABLE(6 ,49 ) /.3806575      /
       DATA TABLE(1 ,53 ) /500000.00     /
       DATA TABLE(2 ,53 ) /301200.00     /
       DATA TABLE(3 ,53 ) /2481.0        /
       DATA TABLE(4 ,53 ) /18.72150      /
       DATA TABLE(5 ,53 ) /.9999428571   /
       DATA TABLE(6 ,53 ) /.3807283      /
       DATA TABLE(1 ,54 ) /500000.00     /
       DATA TABLE(2 ,54 ) /308700.00     /
       DATA TABLE(3 ,54 ) /2481.0        /
       DATA TABLE(4 ,54 ) /18.72150      /
       DATA TABLE(5 ,54 ) /.9999090909   /
       DATA TABLE(6 ,54 ) /.3807541      /
       DATA TABLE(1 ,55 ) /500000.00     /
       DATA TABLE(2 ,55 ) /319500.00     /
       DATA TABLE(3 ,55 ) /2481.0        /
       DATA TABLE(4 ,55 ) /18.72150      /
       DATA TABLE(5 ,55 ) /.9999090909   /
       DATA TABLE(6 ,55 ) /.3805361      /
       DATA TABLE(1 ,62 ) /500000.00     /
       DATA TABLE(2 ,62 ) /319800.00     /
       DATA TABLE(3 ,62 ) /1772.0        /
       DATA TABLE(4 ,62 ) /28.62716      /
       DATA TABLE(5 ,62 ) /.9999600000   /
       DATA TABLE(6 ,62 ) /.3817257      /
       DATA TABLE(1 ,63 ) /500000.00     /
       DATA TABLE(2 ,63 ) /325200.00     /
       DATA TABLE(3 ,63 ) /1822.0        /
       DATA TABLE(4 ,63 ) /21.00903      /
       DATA TABLE(5 ,63 ) /.9999411765   /
       DATA TABLE(6 ,63 ) /.3816986      /
       DATA TABLE(1 ,64 ) /500000.00     /
       DATA TABLE(2 ,64 ) /325800.00     /
       DATA TABLE(3 ,64 ) /2141.0        /
       DATA TABLE(4 ,64 ) /41.66790      /
       DATA TABLE(5 ,64 ) /.9999333333   /
       DATA TABLE(6 ,64 ) /.3812643      /
       DATA TABLE(1 ,65 ) /500000.00     /
       DATA TABLE(2 ,65 ) /333000.00     /
       DATA TABLE(3 ,65 ) /2141.0        /
       DATA TABLE(4 ,65 ) /41.66790      /
       DATA TABLE(5 ,65 ) /.9999333333   /
       DATA TABLE(6 ,65 ) /.3812422      /
       DATA TABLE(1 ,66 ) /500000.00     /
       DATA TABLE(2 ,66 ) /340200.00     /
       DATA TABLE(3 ,66 ) /2161.0        /
       DATA TABLE(4 ,66 ) /39.76857      /
       DATA TABLE(5 ,66 ) /.9999411765   /
       DATA TABLE(6 ,66 ) /.3812362      /
       DATA TABLE(1 ,72 ) /500000.00     /
       DATA TABLE(2 ,72 ) /416100.00     /
       DATA TABLE(3 ,72 ) /2076.0        /
       DATA TABLE(4 ,72 ) /48.30429      /
       DATA TABLE(5 ,72 ) /.9999000000   /
       DATA TABLE(6 ,72 ) /.3812311      /
       DATA TABLE(1 ,73 ) /500000.00     /
       DATA TABLE(2 ,73 ) /420000.00     /
       DATA TABLE(3 ,73 ) /2076.0        /
       DATA TABLE(4 ,73 ) /48.30429      /
       DATA TABLE(5 ,73 ) /.9999000000   /
       DATA TABLE(6 ,73 ) /.3812311      /
       DATA TABLE(1 ,74 ) /500000.00     /
       DATA TABLE(2 ,74 ) /426900.00     /
       DATA TABLE(3 ,74 ) /2076.0        /
       DATA TABLE(4 ,74 ) /48.30429      /
       DATA TABLE(5 ,74 ) /.9999000000   /
       DATA TABLE(6 ,74 ) /.3812311      /
       DATA TABLE(1 ,75 ) /500000.00     /
       DATA TABLE(2 ,75 ) /258000.00     /
       DATA TABLE(3 ,75 ) /2541.0        /
       DATA TABLE(4 ,75 ) /16.76677      /
       DATA TABLE(5 ,75 ) /.9999666667   /
       DATA TABLE(6 ,75 ) /.3807327      /
       DATA TABLE(1 ,76 ) /2000000.00    /
       DATA TABLE(2 ,76 ) /268800.00     /
       DATA TABLE(3 ,76 ) /2321.0        /
       DATA TABLE(4 ,76 ) /27.02745      /
       DATA TABLE(5 ,76 ) /.9999750295   /
       DATA TABLE(6 ,76 ) /.3810845      /
       DATA TABLE(1 ,77 ) /500000.00     /
       DATA TABLE(2 ,77 ) /375600.00     /
       DATA TABLE(3 ,77 ) /1852.0        /
       DATA TABLE(4 ,77 ) /16.62358      /
       DATA TABLE(5 ,77 ) /.9999090909   /
       DATA TABLE(6 ,77 ) /.3816135      /
       DATA TABLE(1 ,78 ) /500000.00     /
       DATA TABLE(2 ,78 ) /382500.00     /
       DATA TABLE(3 ,78 ) /1852.0        /
       DATA TABLE(4 ,78 ) /16.62358      /
       DATA TABLE(5 ,78 ) /.9999000000   /
       DATA TABLE(6 ,78 ) /.3816204      /
       DATA TABLE(1 ,79 ) /500000.00     /
       DATA TABLE(2 ,79 ) /388200.00     /
       DATA TABLE(3 ,79 ) /1852.0        /
       DATA TABLE(4 ,79 ) /16.62358      /
       DATA TABLE(5 ,79 ) /.9999166667   /
       DATA TABLE(6 ,79 ) /.3816288      /
       DATA TABLE(1 ,80 ) /500000.00     /
       DATA TABLE(2 ,80 ) /267600.00     /
       DATA TABLE(3 ,80 ) /2391.0        /
       DATA TABLE(4 ,80 ) /22.84247      /
       DATA TABLE(5 ,80 ) /.9999666667   /
       DATA TABLE(6 ,80 ) /.3808377      /
       DATA TABLE(1 ,81 ) /500000.00     /
       DATA TABLE(2 ,81 ) /275700.00     /
       DATA TABLE(3 ,81 ) /2391.0        /
       DATA TABLE(4 ,81 ) /22.84247      /
       DATA TABLE(5 ,81 ) /.9999375000   /
       DATA TABLE(6 ,81 ) /.3808450      /
       DATA TABLE(1 ,82 ) /500000.00     /
       DATA TABLE(2 ,82 ) /282900.00     /
       DATA TABLE(3 ,82 ) /2391.0        /
       DATA TABLE(4 ,82 ) /22.84247      /
       DATA TABLE(5 ,82 ) /.9999375000   /
       DATA TABLE(6 ,82 ) /.3808750      /
       DATA TABLE(1 ,95 ) /500000.00     /
       DATA TABLE(2 ,95 ) /257400.00     /
       DATA TABLE(3 ,95 ) /2456.0        /
       DATA TABLE(4 ,95 ) /19.72344      /
       DATA TABLE(5 ,95 ) /.9999937500   /
       DATA TABLE(6 ,95 ) /.3809220      /
       DATA TABLE(1 ,109) /500000.00     /
       DATA TABLE(2 ,109) /261000.00     /
       DATA TABLE(3 ,109) /2541.0        /
       DATA TABLE(4 ,109) /16.76677      /
       DATA TABLE(5 ,109) /.9999642857   /
       DATA TABLE(6 ,109) /.3807420      /
       DATA TABLE(1 ,119) /500000.00     /
       DATA TABLE(2 ,119) /378600.00     /
       DATA TABLE(3 ,119) /2431.0        /
       DATA TABLE(4 ,119) /20.83533      /
       DATA TABLE(5 ,119) /.9999411765   /
       DATA TABLE(6 ,119) /.3808422      /
       DATA TABLE(1 ,120) /500000.00     /
       DATA TABLE(2 ,120) /386400.00     /
       DATA TABLE(3 ,120) /2431.0        /
       DATA TABLE(4 ,120) /20.83533      /
       DATA TABLE(5 ,120) /.9999411765   /
       DATA TABLE(6 ,120) /.3808422      /
       DATA TABLE(1 ,121) /500000.00     /
       DATA TABLE(2 ,121) /391500.00     /
       DATA TABLE(3 ,121) /2431.0        /
       DATA TABLE(4 ,121) /20.83533      /
       DATA TABLE(5 ,121) /.9999411765   /
       DATA TABLE(6 ,121) /.3808422      /
       DATA TABLE(1 ,122) /500000.00     /
       DATA TABLE(2 ,122) /396300.00     /
       DATA TABLE(3 ,122) /2431.0        /
       DATA TABLE(4 ,122) /20.83533      /
       DATA TABLE(5 ,122) /.9999411765   /
       DATA TABLE(6 ,122) /.3808422      /
       DATA TABLE(1 ,3  ) /3000000.00    /
       DATA TABLE(2 ,3  ) /633600.00     /
       DATA TABLE(3 ,3  ) /15893950.36   /
       DATA TABLE(4 ,3  ) /16564628.77   /
       DATA TABLE(5 ,3  ) /.9998480641   /
       DATA TABLE(6 ,3  ) /.7969223940   /
       DATA TABLE(7 ,3  ) /3161.0        /
       DATA TABLE(8 ,3  ) /47.87068      /
       DATA TABLE(9 ,3  ) /3.79919       /
       DATA TABLE(10,3  ) /5.91550       /
       DATA TABLE(11,3  ) /44.0          /
       DATA TABLE(1 ,8  ) /2000000.00    /
       DATA TABLE(2 ,8  ) /331200.00     /
       DATA TABLE(3 ,8  ) /29277593.61   /
       DATA TABLE(4 ,8  ) /29732882.87   /
       DATA TABLE(5 ,8  ) /.9999359370   /
       DATA TABLE(6 ,8  ) /.5818991407   /
       DATA TABLE(7 ,8  ) /2126.0        /
       DATA TABLE(8 ,8  ) /46.35656      /
       DATA TABLE(9 ,8  ) /3.81452       /
       DATA TABLE(10,8  ) /3.26432       /
       DATA TABLE(11,8  ) /0.0           /
       DATA TABLE(1 ,9  ) /2000000.00    /
       DATA TABLE(2 ,9  ) /331200.00     /
       DATA TABLE(3 ,9  ) /31014039.23   /
       DATA TABLE(4 ,9  ) /31511724.20   /
       DATA TABLE(5 ,9  ) /.9999184698   /
       DATA TABLE(6 ,9  ) /.5596906871   /
       DATA TABLE(7 ,9  ) /2033.0        /
       DATA TABLE(8 ,9  ) /56.94711      /
       DATA TABLE(9 ,9  ) /3.81550       /
       DATA TABLE(10,9  ) /3.08256       /
       DATA TABLE(11,9  ) /0.0           /
       DATA TABLE(1 ,10 ) /2000000.00    /
       DATA TABLE(2 ,10 ) /439200.00     /
       DATA TABLE(3 ,10 ) /24245358.05   /
       DATA TABLE(4 ,10 ) /24792436.23   /
       DATA TABLE(5 ,10 ) /.9998946358   /
       DATA TABLE(6 ,10 ) /.6538843192   /
       DATA TABLE(7 ,10 ) /2441.0        /
       DATA TABLE(8 ,10 ) /26.75847      /
       DATA TABLE(9 ,10 ) /3.80992       /
       DATA TABLE(10,10 ) /3.93575       /
       DATA TABLE(11,10 ) /0.0           /
       DATA TABLE(1 ,11 ) /2000000.00    /
       DATA TABLE(2 ,11 ) /439200.00     /
       DATA TABLE(3 ,11 ) /25795850.31   /
       DATA TABLE(4 ,11 ) /26312257.65   /
       DATA TABLE(5 ,11 ) /.9999146793   /
       DATA TABLE(6 ,11 ) /.6304679732   /
       DATA TABLE(7 ,11 ) /2336.0        /
       DATA TABLE(8 ,11 ) /30.81964      /
       DATA TABLE(9 ,11 ) /3.81147       /
       DATA TABLE(10,11 ) /3.70114       /
       DATA TABLE(11,11 ) /0.0           /
       DATA TABLE(1 ,12 ) /2000000.00    /
       DATA TABLE(2 ,12 ) /433800.00     /
       DATA TABLE(3 ,12 ) /27057475.85   /
       DATA TABLE(4 ,12 ) /27512992.04   /
       DATA TABLE(5 ,12 ) /.9999291792   /
       DATA TABLE(6 ,12 ) /.6122320427   /
       DATA TABLE(7 ,12 ) /2256.0        /
       DATA TABLE(8 ,12 ) /35.52018      /
       DATA TABLE(9 ,12 ) /3.81265       /
       DATA TABLE(10,12 ) /3.52998       /
       DATA TABLE(11,12 ) /0.0           /
       DATA TABLE(1 ,13 ) /2000000.00    /
       DATA TABLE(2 ,13 ) /428400.00     /
       DATA TABLE(3 ,13 ) /28182405.33   /
       DATA TABLE(4 ,13 ) /28652931.96   /
       DATA TABLE(5 ,13 ) /.9999407628   /
       DATA TABLE(6 ,13 ) /.5965871443   /
       DATA TABLE(7 ,13 ) /2189.0        /
       DATA TABLE(8 ,13 ) /10.35494      /
       DATA TABLE(9 ,13 ) /3.81362       /
       DATA TABLE(10,13 ) /3.39020       /
       DATA TABLE(11,13 ) /0.0           /
       DATA TABLE(1 ,14 ) /2000000.00    /
       DATA TABLE(2 ,14 ) /424800.00     /
       DATA TABLE(3 ,14 ) /30194145.54   /
       DATA TABLE(4 ,14 ) /30649424.27   /
       DATA TABLE(5 ,14 ) /.9999221277   /
       DATA TABLE(6 ,14 ) /.5700119219   /
       DATA TABLE(7 ,14 ) /2076.0        /
       DATA TABLE(8 ,14 ) /52.10305      /
       DATA TABLE(9 ,14 ) /3.81523       /
       DATA TABLE(10,14 ) /3.16593       /
       DATA TABLE(11,14 ) /0.0           /
       DATA TABLE(1 ,15 ) /2000000.00    /
       DATA TABLE(2 ,15 ) /418500.00     /
       DATA TABLE(3 ,15 ) /31846570.92   /
       DATA TABLE(4 ,15 ) /32271267.72   /
       DATA TABLE(5 ,15 ) /.9999541438   /
       DATA TABLE(6 ,15 ) /.5495175982   /
       DATA TABLE(7 ,15 ) /1992.0        /
       DATA TABLE(8 ,15 ) /00.16335      /
       DATA TABLE(9 ,15 ) /3.81642       /
       DATA TABLE(10,15 ) /3.00292       /
       DATA TABLE(11,15 ) /0.0           /
       DATA TABLE(1 ,16 ) /4186692.58    /
       DATA TABLE(2 ,16 ) /426000.00     /
       DATA TABLE(3 ,16 ) /30891382.10   /
       DATA TABLE(4 ,16 ) /35055396.31   /
       DATA TABLE(5 ,16 ) /.9999885350   /
       DATA TABLE(6 ,16 ) /.5612432071   /
       DATA TABLE(7 ,16 ) /2040.0        /
       DATA TABLE(8 ,16 ) /22.88096      /
       DATA TABLE(9 ,16 ) /3.81572       /
       DATA TABLE(10,16 ) /3.09520       /
       DATA TABLE(11,16 ) /0.0           /
       DATA TABLE(1 ,17 ) /2000000.00    /
       DATA TABLE(2 ,17 ) /379800.00     /
       DATA TABLE(3 ,17 ) /24751897.68   /
       DATA TABLE(4 ,17 ) /25086068.20   /
       DATA TABLE(5 ,17 ) /.9999568475   /
       DATA TABLE(6 ,17 ) /.6461334829   /
       DATA TABLE(7 ,17 ) /2406.0        /
       DATA TABLE(8 ,17 ) /24.62308      /
       DATA TABLE(9 ,17 ) /3.81044       /
       DATA TABLE(10,17 ) /3.85610       /
       DATA TABLE(11,17 ) /0.0           /
       DATA TABLE(1 ,18 ) /2000000.00    /
       DATA TABLE(2 ,18 ) /379800.00     /
       DATA TABLE(3 ,18 ) /25781376.91   /
       DATA TABLE(4 ,18 ) /26243052.74   /
       DATA TABLE(5 ,18 ) /.9999359117   /
       DATA TABLE(6 ,18 ) /.6306895773   /
       DATA TABLE(7 ,18 ) /2337.0        /
       DATA TABLE(8 ,18 ) /29.65162      /
       DATA TABLE(9 ,18 ) /3.81146       /
       DATA TABLE(10,18 ) /3.70326       /
       DATA TABLE(11,18 ) /0.0           /
       DATA TABLE(1 ,19 ) /2000000.00    /
       DATA TABLE(2 ,19 ) /379800.00     /
       DATA TABLE(3 ,19 ) /26977133.89   /
       DATA TABLE(4 ,19 ) /27402231.82   /
       DATA TABLE(5 ,19 ) /.9999453995   /
       DATA TABLE(6 ,19 ) /.6133780528   /
       DATA TABLE(7 ,19 ) /2261.0        /
       DATA TABLE(8 ,19 ) /34.26662      /
       DATA TABLE(9 ,19 ) /3.81257       /
       DATA TABLE(10,19 ) /3.54046       /
       DATA TABLE(11,19 ) /0.0           /
       DATA TABLE(1 ,20 ) /600000.00     /
       DATA TABLE(2 ,20 ) /261900.00     /
       DATA TABLE(3 ,20 ) /23659233.56   /
       DATA TABLE(4 ,20 ) /23914389.02   /
       DATA TABLE(5 ,20 ) /.9999831405   /
       DATA TABLE(6 ,20 ) /.6630594147   /
       DATA TABLE(7 ,20 ) /2483.0        /
       DATA TABLE(8 ,20 ) /19.67980      /
       DATA TABLE(9 ,20 ) /3.80929       /
       DATA TABLE(10,20 ) /4.03278       /
       DATA TABLE(11,20 ) /0.0           /
       DATA TABLE(1 ,24 ) /2000000.00    /
       DATA TABLE(2 ,24 ) /304200.00     /
       DATA TABLE(3 ,24 ) /36030443.05   /
       DATA TABLE(4 ,24 ) /36454924.53   /
       DATA TABLE(5 ,24 ) /.9999484343   /
       DATA TABLE(6 ,24 ) /.5025259000   /
       DATA TABLE(7 ,24 ) /1802.0        /
       DATA TABLE(8 ,24 ) /26.11701      /
       DATA TABLE(9 ,24 ) /3.81898       /
       DATA TABLE(10,24 ) /2.65643       /
       DATA TABLE(11,24 ) /0.0           /
       DATA TABLE(1 ,39 ) /2000000.00    /
       DATA TABLE(2 ,39 ) /336600.00     /
       DATA TABLE(3 ,39 ) /22736950.34   /
       DATA TABLE(4 ,39 ) /23162461.59   /
       DATA TABLE(5 ,39 ) /.9999453686   /
       DATA TABLE(6 ,39 ) /.6777445518   /
       DATA TABLE(7 ,39 ) /2551.0        /
       DATA TABLE(8 ,39 ) /20.02265      /
       DATA TABLE(9 ,39 ) /3.80827       /
       DATA TABLE(10,39 ) /4.19479       /
       DATA TABLE(11,39 ) /0.0           /
       DATA TABLE(1 ,40 ) /2000000.00    /
       DATA TABLE(2 ,40 ) /336600.00     /
       DATA TABLE(3 ,40 ) /23936585.11   /
       DATA TABLE(4 ,40 ) /24374096.67   /
       DATA TABLE(5 ,40 ) /.9999483705   /
       DATA TABLE(6 ,40 ) /.6587010213   /
       DATA TABLE(7 ,40 ) /2463.0        /
       DATA TABLE(8 ,40 ) /22.59905      /
       DATA TABLE(9 ,40 ) /3.80959       /
       DATA TABLE(10,40 ) /3.98630       /
       DATA TABLE(11,40 ) /0.0           /
       DATA TABLE(1 ,41 ) /2000000.00    /
       DATA TABLE(2 ,41 ) /352800.00     /
       DATA TABLE(3 ,41 ) /25644959.12   /
       DATA TABLE(4 ,41 ) /25979068.57   /
       DATA TABLE(5 ,41 ) /.9999568556   /
       DATA TABLE(6 ,41 ) /.6327148646   /
       DATA TABLE(7 ,41 ) /2346.0        /
       DATA TABLE(8 ,41 ) /27.97215      /
       DATA TABLE(9 ,41 ) /3.81133       /
       DATA TABLE(10,41 ) /3.72376       /
       DATA TABLE(11,41 ) /0.0           /
       DATA TABLE(1 ,42 ) /2000000.00    /
       DATA TABLE(2 ,42 ) /354600.00     /
       DATA TABLE(3 ,42 ) /26896024.48   /
       DATA TABLE(4 ,42 ) /27351521.50   /
       DATA TABLE(5 ,42 ) /.9999359200   /
       DATA TABLE(6 ,42 ) /.6145281068   /
       DATA TABLE(7 ,42 ) /2266.0        /
       DATA TABLE(8 ,42 ) /34.41020      /
       DATA TABLE(9 ,42 ) /3.81250       /
       DATA TABLE(10,42 ) /3.55102       /
       DATA TABLE(11,42 ) /0.0           /
       DATA TABLE(1 ,43 ) /2000000.00    /
       DATA TABLE(2 ,43 ) /303300.00     /
       DATA TABLE(3 ,43 ) /26371820.68   /
       DATA TABLE(4 ,43 ) /26724051.82   /
       DATA TABLE(5 ,43 ) /.9999620817   /
       DATA TABLE(6 ,43 ) /.6220672671   /
       DATA TABLE(7 ,43 ) /2299.0        /
       DATA TABLE(8 ,43 ) /30.63364      /
       DATA TABLE(9 ,43 ) /3.81202       /
       DATA TABLE(10,43 ) /3.62113       /
       DATA TABLE(11,43 ) /0.0           /
       DATA TABLE(1 ,44 ) /2000000.00    /
       DATA TABLE(2 ,44 ) /308700.00     /
       DATA TABLE(3 ,44 ) /27467860.75   /
       DATA TABLE(4 ,44 ) /27832235.64   /
       DATA TABLE(5 ,44 ) /.9999453808   /
       DATA TABLE(6 ,44 ) /.6064623718   /
       DATA TABLE(7 ,44 ) /2231.0        /
       DATA TABLE(8 ,44 ) /36.57874      /
       DATA TABLE(9 ,44 ) /3.81301       /
       DATA TABLE(10,44 ) /3.47771       /
       DATA TABLE(11,44 ) /0.0           /
       DATA TABLE(1 ,45 ) /2000000.00    /
       DATA TABLE(2 ,45 ) /333000.00     /
       DATA TABLE(3 ,45 ) /33624568.36   /
       DATA TABLE(4 ,45 ) /34079629.33   /
       DATA TABLE(5 ,45 ) /.9999147417   /
       DATA TABLE(6 ,45 ) /.5287006734   /
       DATA TABLE(7 ,45 ) /1907.0        /
       DATA TABLE(8 ,45 ) /12.68515      /
       DATA TABLE(9 ,45 ) /3.81758       /
       DATA TABLE(10,45 ) /2.84511       /
       DATA TABLE(11,45 ) /0.0           /
       DATA TABLE(1 ,46 ) /2000000.00    /
       DATA TABLE(2 ,46 ) /328800.00     /
       DATA TABLE(3 ,46 ) /36271389.35   /
       DATA TABLE(4 ,46 ) /36756553.45   /
       DATA TABLE(5 ,46 ) /.9999257458   /
       DATA TABLE(6 ,46 ) /.5000126971   /
       DATA TABLE(7 ,46 ) /1792.0        /
       DATA TABLE(8 ,46 ) /28.55026      /
       DATA TABLE(9 ,46 ) /3.81911       /
       DATA TABLE(10,46 ) /2.63885       /
       DATA TABLE(11,46 ) /0.0           /
       DATA TABLE(1 ,47 ) /2000000.00    /
       DATA TABLE(2 ,47 ) /328800.00     /
       DATA TABLE(3 ,47 ) /41091749.54   /
       DATA TABLE(4 ,47 ) /41576762.39   /
       DATA TABLE(5 ,47 ) /.9998947956   /
       DATA TABLE(6 ,47 ) /.4540068519   /
       DATA TABLE(7 ,47 ) /1612.0        /
       DATA TABLE(8 ,47 ) /59.30342      /
       DATA TABLE(9 ,47 ) /3.82138       /
       DATA TABLE(10,47 ) /2.27436       /
       DATA TABLE(11,47 ) /25.0          /
       DATA TABLE(1 ,50 ) /800000.00     /
       DATA TABLE(2 ,50 ) /277200.00     /
       DATA TABLE(4 ,50 ) /26369112.76   /
       DATA TABLE(5 ,50 ) /.9999498485   /
       DATA TABLE(6 ,50 ) /.6276341196   /
       DATA TABLE(3 ,50 ) /25989474.99   /
       DATA TABLE(7 ,50 ) /2323.0        /
       DATA TABLE(8 ,50 ) /59.69369      /
       DATA TABLE(9 ,50 ) /3.81166       /
       DATA TABLE(10,50 ) /3.67392       /
       DATA TABLE(11,50 ) /0.0           /
       DATA TABLE(1 ,51 ) /600000.00     /
       DATA TABLE(2 ,51 ) /257400.00     /
       DATA TABLE(3 ,51 ) /23111975.14   /
       DATA TABLE(4 ,51 ) /23549477.32   /
       DATA TABLE(5 ,51 ) /.9999645506   /
       DATA TABLE(6 ,51 ) /.6717286561   /
       DATA TABLE(7 ,51 ) /2523.0        /
       DATA TABLE(8 ,51 ) /19.53138      /
       DATA TABLE(9 ,51 ) /3.80870       /
       DATA TABLE(10,51 ) /4.12738       /
       DATA TABLE(11,51 ) /0.0           /
       DATA TABLE(1 ,52 ) /200000.00     /
       DATA TABLE(2 ,52 ) /253800.00     /
       DATA TABLE(3 ,52 ) /23784678.44   /
       DATA TABLE(4 ,52 ) /23924398.02   /
       DATA TABLE(5 ,52 ) /.9999984844   /
       DATA TABLE(6 ,52 ) /.6610953994   /
       DATA TABLE(7 ,52 ) /2474.0        /
       DATA TABLE(8 ,52 ) /19.47463      /
       DATA TABLE(9 ,52 ) /3.80943       /
       DATA TABLE(10,52 ) /4.01174       /
       DATA TABLE(11,52 ) /0.0           /
       DATA TABLE(1 ,56 ) /2000000.00    /
       DATA TABLE(2 ,56 ) /313200.00     /
       DATA TABLE(3 ,56 ) /20041716.18   /
       DATA TABLE(4 ,56 ) /20589420.09   /
       DATA TABLE(5 ,56 ) /.9999410344   /
       DATA TABLE(6 ,56 ) /.7227899381   /
       DATA TABLE(7 ,56 ) /2768.0        /
       DATA TABLE(8 ,56 ) /22.25085      /
       DATA TABLE(9 ,56 ) /3.80501       /
       DATA TABLE(10,56 ) /4.68430       /
       DATA TABLE(11,56 ) /36.0          /
       DATA TABLE(1 ,57 ) /2000000.00    /
       DATA TABLE(2 ,57 ) /303600.00     /
       DATA TABLE(3 ,57 ) /21001715.22   /
       DATA TABLE(4 ,57 ) /21594768.40   /
       DATA TABLE(5 ,57 ) /.9999509058   /
       DATA TABLE(6 ,57 ) /.7064074100   /
       DATA TABLE(7 ,57 ) /2687.0        /
       DATA TABLE(8 ,57 ) /50.76661      /
       DATA TABLE(9 ,57 ) /3.80622       /
       DATA TABLE(10,57 ) /4.46875       /
       DATA TABLE(11,57 ) /35.0          /
       DATA TABLE(1 ,58 ) /2000000.00    /
       DATA TABLE(2 ,58 ) /303600.00     /
       DATA TABLE(3 ,58 ) /22564848.51   /
       DATA TABLE(4 ,58 ) /23069597.22   /
       DATA TABLE(5 ,58 ) /.9999450783   /
       DATA TABLE(6 ,58 ) /.6805292633   /
       DATA TABLE(7 ,58 ) /2564.0        /
       DATA TABLE(8 ,58 ) /22.23938      /
       DATA TABLE(9 ,58 ) /3.80808       /
       DATA TABLE(10,58 ) /4.15706       /
       DATA TABLE(11,58 ) /33.0          /
       DATA TABLE(1 ,59 ) /2000000.00    /
       DATA TABLE(2 ,59 ) /335160.00     /
       DATA TABLE(3 ,59 ) /18984319.62   /
       DATA TABLE(4 ,59 ) /19471398.75   /
       DATA TABLE(5 ,59 ) /.9999028166   /
       DATA TABLE(6 ,59 ) /.7412196637   /
       DATA TABLE(7 ,59 ) /2861.0        /
       DATA TABLE(8 ,59 ) /24.63011      /
       DATA TABLE(9 ,59 ) /3.80362       /
       DATA TABLE(10,59 ) /5.01609       /
       DATA TABLE(11,59 ) /0.0           /
       DATA TABLE(1 ,60 ) /2000000.00    /
       DATA TABLE(2 ,60 ) /339300.00     /
       DATA TABLE(3 ,60 ) /20006679.72   /
       DATA TABLE(4 ,60 ) /20493457.15   /
       DATA TABLE(5 ,60 ) /.9999220223   /
       DATA TABLE(6 ,60 ) /.7233880702   /
       DATA TABLE(7 ,60 ) /2771.0        /
       DATA TABLE(8 ,60 ) /20.89747      /
       DATA TABLE(9 ,60 ) /3.80497       /
       DATA TABLE(10,60 ) /4.76197       /
       DATA TABLE(11,60 ) /0.0           /
       DATA TABLE(1 ,61 ) /2000000.00    /
       DATA TABLE(2 ,61 ) /338400.00     /
       DATA TABLE(3 ,61 ) /21327006.06   /
       DATA TABLE(4 ,61 ) /21874349.14   /
       DATA TABLE(5 ,61 ) /.9999220448   /
       DATA TABLE(6 ,61 ) /.7009277824   /
       DATA TABLE(7 ,61 ) /2661.0        /
       DATA TABLE(8 ,61 ) /20.12517      /
       DATA TABLE(9 ,61 ) /3.80662       /
       DATA TABLE(10,61 ) /4.46959       /
       DATA TABLE(11,61 ) /0.0           /
       DATA TABLE(1 ,67 ) /2000000.00    /
       DATA TABLE(2 ,67 ) /394200.00     /
       DATA TABLE(3 ,67 ) /18689498.40   /
       DATA TABLE(4 ,67 ) /19157874.26   /
       DATA TABLE(5 ,67 ) /.9999714855   /
       DATA TABLE(6 ,67 ) /.7464518080   /
       DATA TABLE(7 ,67 ) /2888.0        /
       DATA TABLE(8 ,67 ) /20.21285      /
       DATA TABLE(9 ,67 ) /3.80322       /
       DATA TABLE(10,67 ) /5.09490       /
       DATA TABLE(11,67 ) /0.0           /
       DATA TABLE(1 ,68 ) /2000000.00    /
       DATA TABLE(2 ,68 ) /394200.00     /
       DATA TABLE(3 ,68 ) /19432939.76   /
       DATA TABLE(4 ,68 ) /19919806.36   /
       DATA TABLE(5 ,68 ) /.9999220151   /
       DATA TABLE(6 ,68 ) /.7333538278   /
       DATA TABLE(7 ,68 ) /2821.0        /
       DATA TABLE(8 ,68 ) /21.96779      /
       DATA TABLE(9 ,68 ) /3.80422       /
       DATA TABLE(10,68 ) /4.90135       /
       DATA TABLE(11,68 ) /0.0           /
       DATA TABLE(1 ,69 ) /2000000.00    /
       DATA TABLE(2 ,69 ) /394200.00     /
       DATA TABLE(3 ,69 ) /20500650.51   /
       DATA TABLE(4 ,69 ) /21096820.93   /
       DATA TABLE(5 ,69 ) /.9999107701   /
       DATA TABLE(6 ,69 ) /.7149012442   /
       DATA TABLE(7 ,69 ) /2729.0        /
       DATA TABLE(8 ,69 ) /21.15820      /
       DATA TABLE(9 ,69 ) /3.80560       /
       DATA TABLE(10,69 ) /4.64814       /
       DATA TABLE(11,69 ) /0.0           /
       DATA TABLE(1 ,70 ) /2000000.00    /
       DATA TABLE(2 ,70 ) /360000.00     /
       DATA TABLE(3 ,70 ) /23004346.29   /
       DATA TABLE(4 ,70 ) /23368977.46   /
       DATA TABLE(5 ,70 ) /.9999645501   /
       DATA TABLE(6 ,70 ) /.6734507906   /
       DATA TABLE(7 ,70 ) /2531.0        /
       DATA TABLE(8 ,70 ) /19.30504      /
       DATA TABLE(9 ,70 ) /3.80858       /
       DATA TABLE(10,70 ) /4.14653       /
       DATA TABLE(11,70 ) /0.0           /
       DATA TABLE(1 ,71 ) /2000000.00    /
       DATA TABLE(2 ,71 ) /358200.00     /
       DATA TABLE(3 ,71 ) /24104561.06   /
       DATA TABLE(4 ,71 ) /24590781.86   /
       DATA TABLE(5 ,71 ) /.9999220725   /
       DATA TABLE(6 ,71 ) /.6560764003   /
       DATA TABLE(7 ,71 ) /2451.0        /
       DATA TABLE(8 ,71 ) /24.68139      /
       DATA TABLE(9 ,71 ) /3.80977       /
       DATA TABLE(10,71 ) /3.95865       /
       DATA TABLE(11,71 ) /0.0           /
       DATA TABLE(1 ,83 ) /2000000.00    /
       DATA TABLE(2 ,83 ) /266400.00     /
       DATA TABLE(3 ,83 ) /24235000.80   /
       DATA TABLE(4 ,83 ) /24462545.30   /
       DATA TABLE(5 ,83 ) /.9999949000   /
       DATA TABLE(6 ,83 ) /.6540820950   /
       DATA TABLE(7 ,83 ) /2442.0        /
       DATA TABLE(8 ,83 ) /20.64240      /
       DATA TABLE(9 ,83 ) /3.80990       /
       DATA TABLE(10,83 ) /3.93780       /
       DATA TABLE(11,83 ) /0.0           /
       DATA TABLE(1 ,84 ) /2000000.00    /
       DATA TABLE(2 ,84 ) /284400.00     /
       DATA TABLE(3 ,84 ) /29637059.47   /
       DATA TABLE(4 ,84 ) /30183611.25   /
       DATA TABLE(5 ,84 ) /.9998725510   /
       DATA TABLE(6 ,84 ) /.5771707700   /
       DATA TABLE(7 ,84 ) /2106.0        /
       DATA TABLE(8 ,84 ) /51.60353      /
       DATA TABLE(9 ,84 ) /3.81480       /
       DATA TABLE(10,84 ) /3.22483       /
       DATA TABLE(11,84 ) /0.0           /
       DATA TABLE(1 ,85 ) /2000000.00    /
       DATA TABLE(2 ,85 ) /361800.00     /
       DATA TABLE(3 ,85 ) /18819849.05   /
       DATA TABLE(4 ,85 ) /19215516.01   /
       DATA TABLE(5 ,85 ) /.9999358426   /
       DATA TABLE(6 ,85 ) /.7441333961   /
       DATA TABLE(7 ,85 ) /2876.0        /
       DATA TABLE(8 ,85 ) /22.57950      /
       DATA TABLE(9 ,85 ) /3.80339       /
       DATA TABLE(10,85 ) /5.05972       /
       DATA TABLE(11,85 ) /0.0           /
       DATA TABLE(1 ,86 ) /2000000.00    /
       DATA TABLE(2 ,86 ) /361800.00     /
       DATA TABLE(3 ,86 ) /19661027.79   /
       DATA TABLE(4 ,86 ) /20086977.18   /
       DATA TABLE(5 ,86 ) /.9999358523   /
       DATA TABLE(6 ,86 ) /.7293826040   /
       DATA TABLE(7 ,86 ) /2801.0        /
       DATA TABLE(8 ,86 ) /20.45445      /
       DATA TABLE(9 ,86 ) /3.80452       /
       DATA TABLE(10,86 ) /4.84504       /
       DATA TABLE(11,86 ) /0.0           /
       DATA TABLE(1 ,87 ) /2000000.00    /
       DATA TABLE(2 ,87 ) /297000.00     /
       DATA TABLE(3 ,87 ) /24048738.51   /
       DATA TABLE(4 ,87 ) /24559158.47   /
       DATA TABLE(5 ,87 ) /.9999391411   /
       DATA TABLE(6 ,87 ) /.6569503193   /
       DATA TABLE(7 ,87 ) /2455.0        /
       DATA TABLE(8 ,87 ) /23.48125      /
       DATA TABLE(9 ,87 ) /3.80971       /
       DATA TABLE(10,87 ) /3.96783       /
       DATA TABLE(11,87 ) /0.0           /
       DATA TABLE(1 ,88 ) /2000000.00    /
       DATA TABLE(2 ,88 ) /297000.00     /
       DATA TABLE(3 ,88 ) /25522875.81   /
       DATA TABLE(4 ,88 ) /26027071.12   /
       DATA TABLE(5 ,88 ) /.9999359346   /
       DATA TABLE(6 ,88 ) /.6345195439   /
       DATA TABLE(7 ,88 ) /2354.0        /
       DATA TABLE(8 ,88 ) /28.63705      /
       DATA TABLE(9 ,88 ) /3.81121       /
       DATA TABLE(10,88 ) /3.74048       /
       DATA TABLE(11,88 ) /0.0           /
       DATA TABLE(1 ,89 ) /2000000.00    /
       DATA TABLE(2 ,89 ) /352800.00     /
       DATA TABLE(3 ,89 ) /28657871.66   /
       DATA TABLE(4 ,89 ) /29082831.70   /
       DATA TABLE(5 ,89 ) /.9999454101   /
       DATA TABLE(6 ,89 ) /.5901470744   /
       DATA TABLE(7 ,89 ) /2161.0        /
       DATA TABLE(8 ,89 ) /42.56887      /
       DATA TABLE(9 ,89 ) /3.81402       /
       DATA TABLE(10,89 ) /3.33440       /
       DATA TABLE(11,89 ) /0.0           /
       DATA TABLE(1 ,90 ) /2000000.00    /
       DATA TABLE(2 ,90 ) /352800.00     /
       DATA TABLE(3 ,90 ) /30382831.06   /
       DATA TABLE(4 ,90 ) /30838032.96   /
       DATA TABLE(5 ,90 ) /.9999359432   /
       DATA TABLE(6 ,90 ) /.5676166827   /
       DATA TABLE(7 ,90 ) /2066.0        /
       DATA TABLE(8 ,90 ) /52.48935      /
       DATA TABLE(9 ,90 ) /3.81537       /
       DATA TABLE(10,90 ) /3.14645       /
       DATA TABLE(11,90 ) /0.0           /
       DATA TABLE(1 ,91 ) /2000000.00    /
       DATA TABLE(2 ,91 ) /433800.00     /
       DATA TABLE(3 ,91 ) /20836250.94   /
       DATA TABLE(4 ,91 ) /21383852.48   /
       DATA TABLE(5 ,91 ) /.9998945810   /
       DATA TABLE(6 ,91 ) /.7091860222   /
       DATA TABLE(7 ,91 ) /2701.0        /
       DATA TABLE(8 ,91 ) /22.08858      /
       DATA TABLE(9 ,91 ) /3.80602       /
       DATA TABLE(10,91 ) /4.57382       /
       DATA TABLE(11,91 ) /0.0           /
       DATA TABLE(1 ,92 ) /2000000.00    /
       DATA TABLE(2 ,92 ) /433800.00     /
       DATA TABLE(3 ,92 ) /22341309.43   /
       DATA TABLE(4 ,92 ) /22888667.15   /
       DATA TABLE(5 ,92 ) /.9998946058   /
       DATA TABLE(6 ,92 ) /.6841473833   /
       DATA TABLE(7 ,92 ) /2581.0        /
       DATA TABLE(8 ,92 ) /22.74104      /
       DATA TABLE(9 ,92 ) /3.80782       /
       DATA TABLE(10,92 ) /4.26823       /
       DATA TABLE(11,92 ) /0.0           /
       DATA TABLE(1 ,93 ) /2000000.00    /
       DATA TABLE(2 ,93 ) /279900.00     /
       DATA TABLE(3 ,93 ) /23755351.27   /
       DATA TABLE(4 ,93 ) /24211050.37   /
       DATA TABLE(5 ,93 ) /.9999568410   /
       DATA TABLE(6 ,93 ) /.6615397363   /
       DATA TABLE(7 ,93 ) /2476.0        /
       DATA TABLE(8 ,93 ) /21.57953      /
       DATA TABLE(9 ,93 ) /3.80940       /
       DATA TABLE(10,93 ) /4.01753       /
       DATA TABLE(11,93 ) /0.0           /
       DATA TABLE(1 ,94 ) /2000000.00    /
       DATA TABLE(2 ,94 ) /279900.00     /
       DATA TABLE(3 ,94 ) /24577800.67   /
       DATA TABLE(4 ,94 ) /24984826.43   /
       DATA TABLE(5 ,94 ) /.9999595012   /
       DATA TABLE(6 ,94 ) /.6487931668   /
       DATA TABLE(7 ,94 ) /2418.0        /
       DATA TABLE(8 ,94 ) /23.87979      /
       DATA TABLE(9 ,94 ) /3.81026       /
       DATA TABLE(10,94 ) /3.88319       /
       DATA TABLE(11,94 ) /0.0           /
       DATA TABLE(1 ,96 ) /2000000.00    /
       DATA TABLE(2 ,96 ) /291600.00     /
       DATA TABLE(3 ,96 ) /30630125.53   /
       DATA TABLE(4 ,96 ) /31127724.75   /
       DATA TABLE(5 ,96 ) /.9999454207   /
       DATA TABLE(6 ,96 ) /.5644973800   /
       DATA TABLE(7 ,96 ) /2053.0        /
       DATA TABLE(8 ,96 ) /53.44099      /
       DATA TABLE(9 ,96 ) /3.81555       /
       DATA TABLE(10,96 ) /3.12127       /
       DATA TABLE(11,96 ) /0.0           /
       DATA TABLE(1 ,97 ) /2000000.00    /
       DATA TABLE(2 ,97 ) /291600.00     /
       DATA TABLE(3 ,97 ) /32252126.30   /
       DATA TABLE(4 ,97 ) /32676887.65   /
       DATA TABLE(5 ,97 ) /.9999326284   /
       DATA TABLE(6 ,97 ) /.5446515700   /
       DATA TABLE(7 ,97 ) /1972.0        /
       DATA TABLE(8 ,97 ) /3.57839       /
       DATA TABLE(9 ,97 ) /3.81669       /
       DATA TABLE(10,97 ) /2.94381       /
       DATA TABLE(11,97 ) /0.0           /
       DATA TABLE(1 ,98 ) /2000000.00    /
       DATA TABLE(2 ,98 ) /360000.00     /
       DATA TABLE(3 ,98 ) /20922704.09   /
       DATA TABLE(4 ,98 ) /21366697.03   /
       DATA TABLE(5 ,98 ) /.9999391116   /
       DATA TABLE(6 ,98 ) /.7077381841   /
       DATA TABLE(7 ,98 ) /2694.0        /
       DATA TABLE(8 ,98 ) /18.93392      /
       DATA TABLE(9 ,98 ) /3.80912       /
       DATA TABLE(10,98 ) /4.55529       /
       DATA TABLE(11,98 ) /0.0           /
       DATA TABLE(1 ,99 ) /2000000.00    /
       DATA TABLE(2 ,99 ) /361200.00     /
       DATA TABLE(3 ,99 ) /21993575.61   /
       DATA TABLE(4 ,99 ) /22461937.05   /
       DATA TABLE(5 ,99 ) /.9999068931   /
       DATA TABLE(6 ,99 ) /.6898519579   /
       DATA TABLE(7 ,99 ) /2608.0        /
       DATA TABLE(8 ,99 ) /21.54370      /
       DATA TABLE(9 ,99 ) /3.80742       /
       DATA TABLE(10,99 ) /4.33519       /
       DATA TABLE(11,99 ) /0.0           /
       DATA TABLE(1 ,100) /2000000.00    /
       DATA TABLE(2 ,100) /309600.00     /
       DATA TABLE(3 ,100) /29010231.09   /
       DATA TABLE(4 ,100) /29535149.91   /
       DATA TABLE(5 ,100) /.9999484030   /
       DATA TABLE(6 ,100) /.5854397296   /
       DATA TABLE(7 ,100) /2141.0        /
       DATA TABLE(8 ,100) /44.28313      /
       DATA TABLE(9 ,100) /3.81431       /
       DATA TABLE(10,100) /3.29422       /
       DATA TABLE(11,100) /0.0           /
       DATA TABLE(1 ,101) /2000000.00    /
       DATA TABLE(2 ,101) /365400.00     /
       DATA TABLE(3 ,101) /29456907.29   /
       DATA TABLE(4 ,101) /29972959.94   /
       DATA TABLE(5 ,101) /.9999108771   /
       DATA TABLE(6 ,101) /.5795358654   /
       DATA TABLE(7 ,101) /2116.0        /
       DATA TABLE(8 ,101) /48.58548      /
       DATA TABLE(9 ,101) /3.81466       /
       DATA TABLE(10,101) /3.24452       /
       DATA TABLE(11,101) /0.0           /
       DATA TABLE(1 ,102) /2000000.00    /
       DATA TABLE(2 ,102) /351000.00     /
       DATA TABLE(3 ,102) /32187809.58   /
       DATA TABLE(4 ,102) /32691654.54   /
       DATA TABLE(5 ,102) /.9998726224   /
       DATA TABLE(6 ,102) /.5453944146   /
       DATA TABLE(7 ,102) /1975.0        /
       DATA TABLE(8 ,102) /5.95074       /
       DATA TABLE(9 ,102) /3.81665       /
       DATA TABLE(10,102) /2.97107       /
       DATA TABLE(11,102) /0.0           /
       DATA TABLE(1 ,103) /2000000.00    /
       DATA TABLE(2 ,103) /361200.00     /
       DATA TABLE(3 ,103) /34851703.46   /
       DATA TABLE(4 ,103) /35337121.23   /
       DATA TABLE(5 ,103) /.9998817443   /
       DATA TABLE(6 ,103) /.5150588857   /
       DATA TABLE(7 ,103) /1852.0        /
       DATA TABLE(8 ,103) /21.62181      /
       DATA TABLE(9 ,103) /3.81832       /
       DATA TABLE(10,103) /2.74550       /
       DATA TABLE(11,103) /0.0           /
       DATA TABLE(1 ,104) /2000000.00    /
       DATA TABLE(2 ,104) /356400.00     /
       DATA TABLE(3 ,104) /37261509.20   /
       DATA TABLE(4 ,104) /37807440.38   /
       DATA TABLE(5 ,104) /.9998632433   /
       DATA TABLE(6 ,104) /.4899126408   /
       DATA TABLE(7 ,104) /1752.0        /
       DATA TABLE(8 ,104) /37.19059      /
       DATA TABLE(9 ,104) /3.81962       /
       DATA TABLE(10,104) /2.56899       /
       DATA TABLE(11,104) /0.0           /
       DATA TABLE(1 ,105) /2000000.00    /
       DATA TABLE(2 ,105) /354600.00     /
       DATA TABLE(3 ,105) /41091749.54   /
       DATA TABLE(4 ,105) /41576762.39   /
       DATA TABLE(5 ,105) /.9998947956   /
       DATA TABLE(6 ,105) /.4540068519   /
       DATA TABLE(7 ,105) /1612.0        /
       DATA TABLE(8 ,105) /59.30342      /
       DATA TABLE(9 ,105) /3.82138       /
       DATA TABLE(10,105) /2.33094       /
       DATA TABLE(11,105) /0.0           /
       DATA TABLE(1 ,106) /2000000.00    /
       DATA TABLE(2 ,106) /401400.00     /
       DATA TABLE(3 ,106) /23894872.45   /
       DATA TABLE(4 ,106) /24229110.29   /
       DATA TABLE(5 ,106) /.9999568422   /
       DATA TABLE(6 ,106) /.6593554910   /
       DATA TABLE(7 ,106) /2466.0        /
       DATA TABLE(8 ,106) /21.96231      /
       DATA TABLE(9 ,106) /3.80955       /
       DATA TABLE(10,106) /3.99323       /
       DATA TABLE(11,106) /0.0           /
       DATA TABLE(1 ,107) /2000000.00    /
       DATA TABLE(2 ,107) /401400.00     /
       DATA TABLE(3 ,107) /25117176.75   /
       DATA TABLE(4 ,107) /25664114.42   /
       DATA TABLE(5 ,107) /.9998988207   /
       DATA TABLE(6 ,107) /.6405785926   /
       DATA TABLE(7 ,107) /2381.0        /
       DATA TABLE(8 ,107) /29.30066      /
       DATA TABLE(9 ,107) /3.81081       /
       DATA TABLE(10,107) /3.80024       /
       DATA TABLE(11,107) /0.0           /
       DATA TABLE(1 ,108) /2000000.00    /
       DATA TABLE(2 ,108) /401400.00     /
       DATA TABLE(3 ,108) /27025955.35   /
       DATA TABLE(4 ,108) /27432812.88   /
       DATA TABLE(5 ,108) /.9999512939   /
       DATA TABLE(6 ,108) /.6126873424   /
       DATA TABLE(7 ,108) /2258.0        /
       DATA TABLE(8 ,108) /34.16878      /
       DATA TABLE(9 ,108) /3.81262       /
       DATA TABLE(10,108) /3.53414       /
       DATA TABLE(11,108) /0.0           /
       DATA TABLE(1 ,110) /2000000.00    /
       DATA TABLE(2 ,110) /282600.00     /
       DATA TABLE(3 ,110) /26230200.09   /
       DATA TABLE(4 ,110) /26576444.45   /
       DATA TABLE(5 ,110) /.9999483859   /
       DATA TABLE(6 ,110) /.6241178597   /
       DATA TABLE(7 ,110) /2308.0        /
       DATA TABLE(8 ,110) /30.78682      /
       DATA TABLE(9 ,110) /3.81189       /
       DATA TABLE(10,110) /3.64047       /
       DATA TABLE(11,110) /0.0           /
       DATA TABLE(1 ,111) /2000000.00    /
       DATA TABLE(2 ,111) /282600.00     /
       DATA TABLE(3 ,111) /27434800.06   /
       DATA TABLE(4 ,111) /27811312.71   /
       DATA TABLE(5 ,111) /.9999454027   /
       DATA TABLE(6 ,111) /.6069248249   /
       DATA TABLE(7 ,111) /2233.0        /
       DATA TABLE(8 ,111) /36.41072      /
       DATA TABLE(9 ,111) /3.81298       /
       DATA TABLE(10,111) /3.48187       /
       DATA TABLE(11,111) /0.0           /
       DATA TABLE(1 ,112) /2000000.00    /
       DATA TABLE(2 ,112) /435000.00     /
       DATA TABLE(3 ,112) /18798081.67   /
       DATA TABLE(4 ,112) /19205863.43   /
       DATA TABLE(5 ,112) /.9999422551   /
       DATA TABLE(6 ,112) /.7445203390   /
       DATA TABLE(7 ,112) /2878.0        /
       DATA TABLE(8 ,112) /22.15711      /
       DATA TABLE(9 ,112) /3.80336       /
       DATA TABLE(10,112) /5.06556       /
       DATA TABLE(11,112) /0.0           /
       DATA TABLE(1 ,113) /2000000.00    /
       DATA TABLE(2 ,113) /433800.00     /
       DATA TABLE(3 ,113) /19832653.52   /
       DATA TABLE(4 ,113) /20289119.60   /
       DATA TABLE(5 ,113) /.9999145875   /
       DATA TABLE(6 ,113) /.7263957947   /
       DATA TABLE(7 ,113) /2786.0        /
       DATA TABLE(8 ,113) /21.72121      /
       DATA TABLE(9 ,113) /3.80474       /
       DATA TABLE(10,113) /4.80336       /
       DATA TABLE(11,113) /0.0           /
       DATA TABLE(1 ,114) /2000000.00    /
       DATA TABLE(2 ,114) /286200.00     /
       DATA TABLE(3 ,114) /25305029.12   /
       DATA TABLE(4 ,114) /25715126.55   /
       DATA TABLE(5 ,114) /.9999407460   /
       DATA TABLE(6 ,114) /.6377729696   /
       DATA TABLE(7 ,114) /2368.0        /
       DATA TABLE(8 ,114) /57.52979      /
       DATA TABLE(9 ,114) /3.81099       /
       DATA TABLE(10,114) /3.77244       /
       DATA TABLE(11,114) /0.0           /
       DATA TABLE(1 ,115) /2000000.00    /
       DATA TABLE(2 ,115) /291600.00     /
       DATA TABLE(3 ,115) /26639323.45   /
       DATA TABLE(4 ,115) /27070620.78   /
       DATA TABLE(5 ,115) /.9999256928   /
       DATA TABLE(6 ,115) /.6181953936   /
       DATA TABLE(7 ,115) /2282.0        /
       DATA TABLE(8 ,115) /33.82207      /
       DATA TABLE(9 ,115) /3.81227       /
       DATA TABLE(10,115) /3.58491       /
       DATA TABLE(11,115) /0.0           /
       DATA TABLE(1 ,116) /2000000.00    /
       DATA TABLE(2 ,116) /324000.00     /
       DATA TABLE(3 ,116) /20124133.05   /
       DATA TABLE(4 ,116) /20489179.67   /
       DATA TABLE(5 ,116) /.9999453461   /
       DATA TABLE(6 ,116) /.7213707913   /
       DATA TABLE(7 ,116) /2761.0        /
       DATA TABLE(8 ,116) /19.04034      /
       DATA TABLE(9 ,116) /3.80511       /
       DATA TABLE(10,116) /4.73451       /
       DATA TABLE(11,116) /0.0           /
       DATA TABLE(1 ,117) /2000000.00    /
       DATA TABLE(2 ,117) /324000.00     /
       DATA TABLE(3 ,117) /21050746.99   /
       DATA TABLE(4 ,117) /21430913.91   /
       DATA TABLE(5 ,117) /.9999407059   /
       DATA TABLE(6 ,117) /.7055766312   /
       DATA TABLE(7 ,117) /2683.0        /
       DATA TABLE(8 ,117) /48.81363      /
       DATA TABLE(9 ,117) /3.80628       /
       DATA TABLE(10,117) /4.52782       /
       DATA TABLE(11,117) /0.0           /
       DATA TABLE(1 ,118) /2000000.00    /
       DATA TABLE(2 ,118) /324000.00     /
       DATA TABLE(3 ,118) /22161432.25   /
       DATA TABLE(4 ,118) /22672134.66   /
       DATA TABLE(5 ,118) /.9999325474   /
       DATA TABLE(6 ,118) /.6871032423   /
       DATA TABLE(7 ,118) /2595.0        /
       DATA TABLE(8 ,118) /20.01691      /
       DATA TABLE(9 ,118) /3.80761       /
       DATA TABLE(10,118) /4.30274       /
       DATA TABLE(11,118) /0.0           /
! ....................................................................
      DATA C1,C2,C3,C4 /101.2794065D0,60.0D0,1052.893882D0,4.483344D0/
      DATA C5,C6,C7 /0.023520D0,100000000.0D0,0.009873675553D0/
      DATA C8,C9,C10 /1047.54671D0,6.19276D0,0.050912D0/
! ....................................................................
      DATA D1,D2,D3,D4 /3.9174D0,30.92241724D0,4.0831D0,3.280833333D0/
      DATA D5,D6,D7 /25.52381D0,0.3048006099D0,4.0831D0/
      DATA D8,D9,D10 /100000.0D0,0.0000000001D0,10000.0D0/
! ....................................................................
      DATA RADSEC /206264.806247D0/
! ....................................................................
      DATA ZERO,ONE,TWO /0.0D0,1.0D0,2.0D0/
      DATA ESQ /0.006768658D0/
      DATA EPSLN,EPSLN1,EPSLN2 /0.0001D0,0.1D0,0.01D0/
      DATA NIT /5/
      DATA SWITCH /0/
!
! ......................................................................
!       .  INITIALIZATION OF PROJECTION PARAMETERS (FILE INPUT)  .
! ......................................................................
!
      ENTRY IF02Z0 (INFILE,data)
!
      IERROR = 0
      READ (INFILE,END=120) ZONE,BUFF
      IF (SWITCH.NE.0 .AND. SWITCH.EQ.ZONE) RETURN
  020 IF (ZONE .LE. 0) GO TO 050
      DO 040 IND = 1,131
      IF (ZONE .EQ. ITEM(IND)) GO TO 060
  040 CONTINUE
  050 IF (IPEMSG .EQ. 0) PRINT 2000, ZONE
 2000 FORMAT (' ERROR PJ02Z0'/                     &
              ' ILLEGAL ZONE NO : ',I10)
      IERROR = 021
      RETURN
  060 ITYPE = ID(9,IND) + 1
      GO TO (080,080,100) , ITYPE
  080 T1 = TABLE(1,IND)
      T2 = TABLE(2,IND)
      T3 = TABLE(3,IND)
      T4 = TABLE(4,IND)
      T5 = TABLE(5,IND)
      T6 = TABLE(6,IND)
      T7 = TABLE(7,IND)
      T8 = TABLE(8,IND)
      T9 = TABLE(9,IND)
      T10= TABLE(10,IND)
      T11= TABLE(11,IND)
!
! LIST RESULTS OF PARAMETER INITIALIZATION.
!
  100 IF (IPPARM .EQ. 0) PRINT 2010, (ID(I,IND),I=1,8)
 2010 FORMAT (' INITIALIZATION PARAMETERS (STATE PLANE PROJECTION)'/  &
              ' ZONE = ',8A4)
      SWITCH = ZONE
      RETURN
  120 IF (IPEMSG .EQ. 0) PRINT 2020
 2020 FORMAT (' ERROR PJ02Z0'/                     &
              ' MISSING PROJECTION PARAMETERS')
      IERROR = 022
      RETURN
!
! ......................................................................
!      .  INITIALIZATION OF PROJECTION PARAMETERS (ENTRY INPUT)  .
! ......................................................................
!
      ENTRY IS02Z0 (ZZONE,DATA)
      zone = zzone
!
      IERROR = 0
      IF (SWITCH.NE.0 .AND. SWITCH.EQ.ZONE) RETURN
      BUFF(1) = DATA(1)
      GO TO 020
!
! ......................................................................
!                      .  FORWARD TRANSFORMATION  .
! ......................................................................
!
      ENTRY PF02Z0 (GEOG,PROJ)
!
      IERROR = 0
      IF (SWITCH .NE. 0) GO TO 140
      IF (IPEMSG .EQ. 0) PRINT 2020
      IERROR = 023
      RETURN
  140 GO TO (160,220,240) , ITYPE
!
! MERCATOR PROJECTION.
!
  160 FLONG = DABS(GEOG(1) * RADSEC)
      FLAT = GEOG(2) * RADSEC
      DL = T2 - FLONG
      DL1 = DL - D1 * (DL / D10)**3
      SINPH = DSIN(GEOG(2))
      COSPH = DCOS(GEOG(2))
      S1 = D2 * COSPH * DL1 / DSQRT(ONE - ESQ * SINPH * SINPH)
      SM = S1 + D3 * (S1 / D8)**3
      SG = D4 * T5 * SM
      PROJ(1) = T1 + SG + T6 * (SG / D8)**3
      DLPI = ZERO
      DO 180 I = 1,NIT
      PHI = GEOG(2) + DLPI / RADSEC
      SINPH = DSIN(PHI)
      COSPH = DCOS(PHI)
      B = D5 * D9 * SINPH / COSPH * (ONE - ESQ * SINPH * SINPH)**2
      DLPIN = B * SM**2
      IF (DABS(DLPI - DLPIN) .LE. EPSLN) GO TO 200
      DLPI = DLPIN
  180 CONTINUE
      IF (IPEMSG .EQ. 0) PRINT 2030, NIT
 2030 FORMAT (' ERROR PJ02Z0'/                                   &
              ' FAILED TO COVERGE AFTER ',I2,' ITERATIONS')
      IERROR = 024
      RETURN
  200 FLATR = GEOG(2) + DLPIN / RADSEC
      SINPH = DSIN(FLATR)
      COSPH = DCOS(FLATR)
      PROJ(2) = C1 * T5 * (FLAT + DLPIN - C2 * T3 - T4 - SINPH * COSPH *  &
               (C3 - COSPH * COSPH * (C4 - C5 * COSPH * COSPH)))
      RETURN
!
! LAMBERT PROJECTION.
!
  220 FLONG = DABS(GEOG(1) * RADSEC)
      FLAT = GEOG(2) * RADSEC
      SINPH = DSIN(GEOG(2))
      COSPH = DCOS(GEOG(2))
      TEMP = SINPH * COSPH * (C3 - COSPH * COSPH * (C4 - C5 * COSPH * COSPH))
      S = C1 * (C2 * T7 + T8 - FLAT + TEMP)
      TEMP = S / C6
      TEMP1 = T9 - T10 * TEMP + T11 * TEMP * TEMP
      R = T3 + S * T5 * (ONE + TEMP1 * TEMP * TEMP)
      THETA = (T6 * (T2 - FLONG)) / RADSEC
      SINTH = DSIN(THETA)
      SINHTH = DSIN(THETA / TWO)
      PROJ(1) = R * SINTH + T1
      PROJ(2) = T4 - R + TWO * R * SINHTH * SINHTH
      RETURN
!
! SPECIAL ALASKA PROJECTIONS.
!
  240 IZ = IND - 122
      IFLG = 1
      IF (IZ .NE. 1) GO TO 260
      CALL AL01Z0 (GEOG,PROJ,IFLG)
      RETURN
  260 CALL AL29Z0 (GEOG,PROJ,IZ,IFLG)
      RETURN
!
! ......................................................................
!                      .  INVERSE TRANSFORMATION  .
! ......................................................................
!
      ENTRY PI02Z0 (PROJ,GEOG)
!
      IERROR = 0
      IF (SWITCH .NE. 0) GO TO 300
      IF (IPEMSG .EQ. 0) PRINT 2020
      IERROR = 025
      RETURN
  300 GO TO (320,340,420) , ITYPE
!
! MERCATOR PROJECTION.
!
  320 XP = PROJ(1) - T1
      SG1 = XP - T6 * (XP / D8)**3
      SG = XP - T6 * (SG1 / D8)**3
      SM = D6 * SG / T5
      WSEC = C2 * T3 + T4 + C7 * PROJ(2) / T5
      W = WSEC / RADSEC
      SINW = DSIN(W)
      COSW = DCOS(W)
      PHI = (WSEC + SINW * COSW * (C8 + COSW * COSW * (C9 + C10 * COSW * COSW))) / RADSEC
      SINPH = DSIN(PHI)
      COSPH = DCOS(PHI)
      GEOG(2) = PHI - (D5 * D9 * SM * SM * SINPH / COSPH * (ONE - ESQ * SINPH * SINPH)**2) / RADSEC
      SINPH = DSIN(GEOG(2))
      COSPH = DCOS(GEOG(2))
      SA = SM - D7 * (SM / D8)**3
      S1 = SM - D7 * (SA / D8)**3
      DLAM1 = S1 * DSQRT(ONE - ESQ * SINPH * SINPH) / (D2 * COSPH)
      DLAMA = DLAM1 + D1 * (DLAM1 / D10)**3
      GEOG(1) =-(T2 - DLAM1 - D1 * (DLAMA / D10)**3) / RADSEC
      RETURN
!
! LAMBERT PROJECTION.
!
  340 TEMP = (PROJ(1) - T1) / (T4 - PROJ(2))
      THETA = DATAN(TEMP)
      COSTH = DCOS(THETA)
      SINHTH = DSIN(THETA / TWO)
      GEOG(1) =-(T2 / RADSEC) + (THETA / T6)
      R = (T4 - PROJ(2)) / COSTH
      YPRIME = PROJ(2) - TWO * R * SINHTH * SINHTH
      SAPP = (T4 - T3 - YPRIME) / T5
      SAPPP = SAPP
      DO 380 I = 1,NIT
      TEMP = SAPP / C6
      TEMP = T9 * TEMP**2 - T10 * TEMP**3 + T11 * TEMP**4
      S = SAPPP / (ONE + TEMP)
      IF (DABS(S - SAPP) .GT. EPSLN1) GO TO 360
      S = SAPPP * (ONE - TEMP * (ONE - TEMP))
      IF (DABS(S - SAPP) .LE. EPSLN2) GO TO 400
  360 SAPP = S
  380 CONTINUE
      IF (IPEMSG .EQ. 0) PRINT 2030, NIT
      IERROR = 026
      RETURN
  400 WSEC = C2 * T7 + T8 - C7 * S
      W = WSEC / RADSEC
      SINW = DSIN(W)
      COSW = DCOS(W)
      GEOG(2) = (WSEC + SINW * COSW * (C8 + COSW * COSW * (C9 + C10 * COSW * COSW))) / RADSEC
      RETURN
!
! SPECIAL ALASKA PROJECTIONS.
!
  420 IZ = IND - 122
      IFLG = 0
      IF (IZ .NE. 1) GO TO 440
      CALL AL01Z0 (PROJ,GEOG,IFLG)
      RETURN
  440 CALL AL29Z0 (PROJ,GEOG,IZ,IFLG)
      RETURN
!
      END

!                   PJ03Z0
! **********************************************************************
! ** U.S.G.S. GENERAL MAP PROJECTION PACKAGE ...... DR. A. A. ELASSAL **
! **          MATHEMATICAL ANALYSIS BY JOHN SNYDER                    **
! ** MODULE I                VERSION 1.0.0            NOVEMBER 1,1980 **
! **********************************************************************
!                    *  ALBERS CONICAL EQUAL AREA  *
! **********************************************************************
!
      SUBROUTINE PJ03Z0
!
      IMPLICIT REAL*8 (A-Z)
      integer*4 zzone
      INTEGER*4 IERROR,IPEMSG,IPPARM
      INTEGER*4 SWITCH,I,ZONE,ANGS,INFILE
      COMMON /ELLPZ0/ AZ,EZ,ESZ,E0Z,E1Z,E2Z,E3Z
! **** PARAMETERS **** A,E,ES,LAT1,LAT2,LON0,LAT0,X0,Y0,NS,C,RH0 *******
      COMMON /ERRMZ0/ IERROR
      COMMON /PRINZ0/ IPEMSG,IPPARM
!     COMMON /WORKZ0/ BUFF(15),ANGS(4,4)
      COMMON /WORKZ0/ BUFF(15)
      COMMON /WK03Z0/ ANGS(4,4)
      real*4 rangs1,rangs2,rangs3,rangs4
      equivalence (rangs1,angs(4,1))
      equivalence (rangs2,angs(4,2))
      equivalence (rangs3,angs(4,3))
      equivalence (rangs4,angs(4,4))
      DIMENSION DATA(1),GEOG(1),PROJ(1)
      DATA TOL,EPSLN /1.0D-7,1.0D-10/
      DATA HALFPI /1.5707963267948966D0/
      DATA ZERO,HALF,ONE /0.0D0,0.5D0,1.0D0/
      DATA SWITCH /0/
!
! ......................................................................
!       .  INITIALIZATION OF PROJECTION PARAMETERS (FILE INPUT)  .
! ......................................................................
!
      ENTRY IF03Z0 (INFILE,DATA)
!
      IERROR = 0
      READ (INFILE,END=160) ZONE,BUFF
      IF (SWITCH.NE.0 .AND. SWITCH.EQ.ZONE) RETURN
  020 IF (BUFF(1) .LE. ZERO) GO TO 100
      A = BUFF(1)
      B = BUFF(2)
      IF (B .GT. ZERO) GO TO 040
      E = ZERO
      ES = ZERO
      GO TO 120
  040 IF (B .GT. ONE) GO TO 060
      E = DSQRT (B)
      ES = B
      GO TO 120
  060 ES = ONE - (B / A) ** 2
      E = DSQRT (ES)
      GO TO 120
  100 A = AZ
      E = EZ
      ES = ESZ
  120 LAT1 = PAKRZ0 (BUFF(3))
      LAT2 = PAKRZ0 (BUFF(4))
      IF (DABS(LAT1+LAT2) .GE. EPSLN) GO TO 130
      IF (IPEMSG .EQ. 0) PRINT 2000
 2000 FORMAT (' ERROR PJ03Z0'/                                     &
              ' EQUAL LATITUDES FOR ST. PARALLELS ON OPPOSITE',    &
              ' SIDES OF EQUATOR')
      IERROR = 031
      RETURN
  130 LON0 = PAKRZ0 (BUFF(5))
      LAT0 = PAKRZ0 (BUFF(6))
      X0 = BUFF(7)
      Y0 = BUFF(8)
      SINPHI = DSIN (LAT1)
      CON = SINPHI
      COSPHI = DCOS (LAT1)
      MS1 = MSFNZ0 (E,SINPHI,COSPHI)
      QS1 = QSFNZ0 (E,SINPHI,COSPHI)
      SINPHI = DSIN (LAT2)
      COSPHI = DCOS (LAT2)
      MS2 = MSFNZ0 (E,SINPHI,COSPHI)
      QS2 = QSFNZ0 (E,SINPHI,COSPHI)
      SINPHI = DSIN (LAT0)
      COSPHI = DCOS (LAT0)
      QS0 = QSFNZ0 (E,SINPHI,COSPHI)
      IF (DABS(LAT1-LAT2) .GE. EPSLN) GO TO 140
      NS = CON
      GO TO 150
  140 NS = (MS1 * MS1 - MS2 * MS2) / (QS2 - QS1)
  150 C = MS1 * MS1 + NS * QS1
      RH0 = A * DSQRT (C - NS * QS0) / NS
!
! LIST RESULTS OF PARAMETER INITIALIZATION.
!
      CALL RADDZ0 (LAT1,ANGS(1,1))
      CALL RADDZ0 (LAT2,ANGS(1,2))
      CALL RADDZ0 (LON0,ANGS(1,3))
      CALL RADDZ0 (LAT0,ANGS(1,4))
!      IF (IPPARM .EQ. 0) PRINT 2010, A,ES,ANGS,X0,Y0
      IF (IPPARM .EQ. 0) PRINT 2010, A,ES,angs(1,1),angs(2,1),   &
         angs(3,1),rangs1,angs(1,2),angs(2,2),angs(3,2),rangs2,  &
         angs(1,3),angs(2,3),angs(3,3),rangs3,angs(1,4),         &
         angs(2,4),angs(3,4),rangs4,X0,Y0
 2010 FORMAT (' INITIALIZATION PARAMETERS (ALBERS CONICAL EQUAL-AREA PROJECTION)'/  &
              ' SEMI-MAJOR AXIS OF ELLIPSOID =',F12.2,' METERS'/                    &
              ' ECCENTRICITY SQUARED         =',F12.9/                              &
              ' LATITUDE OF 1ST ST. PARALLEL = ',A1,2I3,F7.3/                       &
              ' LATITUDE OF 2ND ST. PARALLEL = ',A1,2I3,F7.3/                       &
              ' LONGITUDE OF ORIGIN          = ',A1,2I3,F7.3/                       &
              ' LATITUDE OF ORIGIN           = ',A1,2I3,F7.3/                       &
              ' FALSE EASTING                =',F12.2,' METERS'/                    & 
              ' FALSE NORTHING               =',F12.2,' METERS')
      DATA(1) = A
      DATA(2) = ES
      SWITCH = ZONE
      RETURN
  160 IF (IPEMSG .EQ. 0) PRINT 2020
 2020 FORMAT (' ERROR PJ03Z0'/                     &
              ' MISSING PROJECTION PARAMETERS')
      IERROR = 032
      RETURN
!
! ......................................................................
!      .  INITIALIZATION OF PROJECTION PARAMETERS (ENTRY INPUT)  .
! ......................................................................
!
      ENTRY IS03Z0 (ZZONE,DATA)
      zone = zzone
!
      IERROR = 0
      IF (SWITCH.NE.0 .AND. SWITCH.EQ.ZONE) RETURN
      DO 180 I = 1,8
      BUFF(I) = DATA(I)
  180 CONTINUE
      GO TO 020
!
! ......................................................................
!                      .  FORWARD TRANSFORMATION  .
! ......................................................................
!
      ENTRY PF03Z0 (GEOG,PROJ)
!
      IERROR = 0
      IF (SWITCH .NE. 0) GO TO 220
      IF (IPEMSG .EQ. 0) PRINT 2020
      IERROR = 033
      RETURN
  220 SINPHI = DSIN (GEOG(2))
      COSPHI = DCOS (GEOG(2))
      QS = QSFNZ0 (E,SINPHI,COSPHI)
      RH = A * DSQRT (C - NS * QS) / NS
      THETA = NS * ADJLZ0 (GEOG(1) - LON0)
      PROJ(1) = X0 + RH * DSIN (THETA)
      PROJ(2) = Y0 + RH0 - RH * DCOS (THETA)
      RETURN
!
! ......................................................................
!                      .  INVERSE TRANSFORMATION  .
! ......................................................................
!
      ENTRY PI03Z0 (PROJ,GEOG)
!
      IERROR = 0
      IF (SWITCH .NE. 0) GO TO 240
      IF (IPEMSG .EQ. 0) PRINT 2020
      IERROR = 034
      RETURN
  240 X = PROJ(1) - X0
      Y = RH0 - PROJ(2) + Y0
      RH = DSIGN (DSQRT (X * X + Y * Y) , NS)
      THETA = ZERO
      CON = DSIGN (ONE , NS)
      IF (RH .NE. ZERO) THETA = DATAN2 (CON * X , CON * Y)
      CON = RH * NS / A
      QS = (C - CON * CON) / NS
      IF (E .LT. TOL) GO TO 260
      CON = ONE - HALF * (ONE - ES) * DLOG ((ONE - E) / (ONE + E)) / E
      IF ((DABS(CON) - DABS(QS)) .GT. TOL) GO TO 260
      GEOG(2) = DSIGN (HALFPI , QS)
      GO TO 280
  260 GEOG(2) = PHI1Z0 (E,QS)
      IF (IERROR .EQ. 0) GO TO 280
      IERROR = 035
      RETURN
  280 GEOG(1) = ADJLZ0 (THETA / NS + LON0)
      RETURN
!
      END

!                   PJ04Z0
! **********************************************************************
! ** U.S.G.S. GENERAL MAP PROJECTION PACKAGE ...... DR. A. A. ELASSAL **
! **          MATHEMATICAL ANALYSIS BY JOHN SNYDER                    **
! ** MODULE I                VERSION 1.0.0            NOVEMBER 1,1980 **
! **********************************************************************
!                     *  LAMBERT CONFORMAL CONIC  *
! **********************************************************************
!
      SUBROUTINE PJ04Z0
!
      IMPLICIT REAL*8 (A-Z)
      integer*4 zzone
      INTEGER*4 IERROR,IPEMSG,IPPARM
      INTEGER*4 SWITCH,I,ZONE,ANGS,INFILE
      COMMON /ELLPZ0/ AZ,EZ,ESZ,E0Z,E1Z,E2Z,E3Z
! **** PARAMETERS **** A,E,ES,LAT1,LAT2,LON0,LAT0,X0,Y0,NS,F,RH0 *******
      COMMON /ERRMZ0/ IERROR
      COMMON /PRINZ0/ IPEMSG,IPPARM
!     COMMON /WORKZ0/ BUFF(15),ANGS(4,4)
      COMMON /WORKZ0/ BUFF(15)
      COMMON /WK04Z0/ ANGS(4,4)
      real*4 rangs1,rangs2,rangs3,rangs4
      equivalence (rangs1,angs(4,1))
      equivalence (rangs2,angs(4,2))
      equivalence (rangs3,angs(4,3))
      equivalence (rangs4,angs(4,4))
      DIMENSION DATA(1),GEOG(1),PROJ(1)
      DATA HALFPI /1.5707963267948966D0/
      DATA EPSLN /1.0D-10/
      DATA ZERO,ONE /0.0D0,1.0D0/
      DATA SWITCH /0/
!
! ......................................................................
!       .  INITIALIZATION OF PROJECTION PARAMETERS (FILE INPUT)  .
! ......................................................................
!
      ENTRY IF04Z0 (INFILE,DATA)
!
      IERROR = 0
      READ (INFILE,END=160) ZONE,BUFF
      IF (SWITCH.NE.0 .AND. SWITCH.EQ.ZONE) RETURN
  020 IF (BUFF(1) .LE. ZERO) GO TO 100
      A = BUFF(1)
      B = BUFF(2)
      IF (B .GT. ZERO) GO TO 040
      E = ZERO
      ES = ZERO
      GO TO 120
  040 IF (B .GT. ONE) GO TO 060
      E = DSQRT (B)
      ES = B
      GO TO 120
  060 ES = ONE - (B / A) ** 2
      E = DSQRT (ES)
      GO TO 120
  100 A = AZ
      E = EZ
      ES = ESZ
  120 LAT1 = PAKRZ0 (BUFF(3))
      LAT2 = PAKRZ0 (BUFF(4))
      IF (DABS(LAT1+LAT2) .GE. EPSLN) GO TO 130
      IF (IPEMSG .EQ. 0) PRINT 2000
 2000 FORMAT (' ERROR PJ04Z0'/                                    &
              ' EQUAL LATITUDES FOR ST. PARALLELS ON OPPOSITE',   &
              ' SIDES OF EQUATOR')
      IERROR = 041
      RETURN
  130 LON0 = PAKRZ0 (BUFF(5))
      LAT0 = PAKRZ0 (BUFF(6))
      X0 = BUFF(7)
      Y0 = BUFF(8)
      SINPHI = DSIN (LAT1)
      CON = SINPHI
      COSPHI = DCOS (LAT1)
      MS1 = MSFNZ0 (E,SINPHI,COSPHI)
      TS1 = TSFNZ0 (E,LAT1,SINPHI)
      SINPHI = DSIN (LAT2)
      COSPHI = DCOS (LAT2)
      MS2 = MSFNZ0 (E,SINPHI,COSPHI)
      TS2 = TSFNZ0 (E,LAT2,SINPHI)
      SINPHI = DSIN (LAT0)
      TS0 = TSFNZ0 (E,LAT0,SINPHI)
      IF (DABS(LAT1-LAT2) .GE. EPSLN) GO TO 140
      NS = CON
      GO TO 150
  140 NS = DLOG (MS1 / MS2) / DLOG (TS1 / TS2)
  150 F = MS1 / (NS * TS1 ** NS)
      RH0 = A * F * TS0 ** NS
!
! LIST RESULTS OF PARAMETER INITIALIZATION.
!
      CALL RADDZ0 (LAT1,ANGS(1,1))
      CALL RADDZ0 (LAT2,ANGS(1,2))
      CALL RADDZ0 (LON0,ANGS(1,3))
      CALL RADDZ0 (LAT0,ANGS(1,4))
!     IF (IPPARM .EQ. 0) PRINT 2010, A,ES,ANGS,X0,Y0
      IF (IPPARM .EQ. 0) PRINT 2010, A,ES,angs(1,1),angs(2,1),  &
         angs(3,1),rangs1,angs(1,2),angs(2,2),angs(3,2),rangs2, &
         angs(1,3),angs(2,3),angs(3,3),rangs3,angs(1,4),        & 
         angs(2,4),angs(3,4),rangs4,X0,Y0
 2010 FORMAT (' INITIALIZATION PARAMETERS (LAMBERT CONFORMAL CONIC',  &
              ' PROJECTION)'/                                         &
              ' SEMI-MAJOR AXIS OF ELLIPSOID =',F12.2,' METERS'/      & 
              ' ECCENTRICITY SQUARED         =',F12.9/                &
              ' LATITUDE OF 1ST ST. PARALLEL = ',A1,2I3,F7.3/         &
              ' LATITUDE OF 2ND ST. PARALLEL = ',A1,2I3,F7.3/         & 
              ' LONGITUDE OF ORIGIN          = ',A1,2I3,F7.3/         &
              ' LATITUDE OF ORIGIN           = ',A1,2I3,F7.3/         &
              ' FALSE EASTING                =',F12.2,' METERS'/      &
              ' FALSE NORTHING               =',F12.2,' METERS')
      DATA(1) = A
      DATA(2) = ES
      SWITCH = ZONE
      RETURN
  160 IF (IPEMSG .EQ. 0) PRINT 2020
 2020 FORMAT (' ERROR PJ04Z0'/                    &
              ' MISSING PROJECTION PARAMETERS')
      IERROR = 042
      RETURN
!
! ......................................................................
!      .  INITIALIZATION OF PROJECTION PARAMETERS (ENTRY INPUT)  .
! ......................................................................
!
      ENTRY IS04Z0 (ZZONE,DATA)
      zone = zzone
!
      IERROR = 0
      IF (SWITCH.NE.0 .AND. SWITCH.EQ.ZONE) RETURN
      DO 180 I = 1,8
      BUFF(I) = DATA(I)
  180 CONTINUE
      GO TO 020
!
! ......................................................................
!                      .  FORWARD TRANSFORMATION  .
! ......................................................................
!
      ENTRY PF04Z0 (GEOG,PROJ)
!
      IERROR = 0
      IF (SWITCH .NE. 0) GO TO 200
      IF (IPEMSG .EQ. 0) PRINT 2020
      IERROR = 043
      RETURN
  200 CON = DABS (DABS (GEOG(2)) - HALFPI)
      IF (CON .GT. EPSLN) GO TO 220
      CON = GEOG(2) * NS
      IF (CON .GT. ZERO) GO TO 210
      if (ipemsg .eq. 0) PRINT 2030
 2030 FORMAT (' ERROR PJ04Z0'/              &
              ' POINT CANNOT BE PROJECTED')
      IERROR = 044
      RETURN
  210 RH = ZERO
      GO TO 230
  220 SINPHI = DSIN (GEOG(2))
      TS = TSFNZ0 (E,GEOG(2),SINPHI)
      RH = A * F * TS ** NS
  230 THETA = NS * ADJLZ0 (GEOG(1) - LON0)
      PROJ(1) = X0 + RH * DSIN (THETA)
      PROJ(2) = Y0 + RH0 - RH * DCOS (THETA)
      RETURN
!
! ......................................................................
!                      .  INVERSE TRANSFORMATION  .
! ......................................................................
!
      ENTRY PI04Z0 (PROJ,GEOG)
!
      IERROR = 0
      IF (SWITCH .NE. 0) GO TO 240
      IF (IPEMSG .EQ. 0) PRINT 2020
      IERROR = 045
      RETURN
  240 X = PROJ(1) - X0
      Y = RH0 - PROJ(2) + Y0
      RH = DSIGN (DSQRT (X*X + Y*Y) , NS)
      THETA = ZERO
      CON = DSIGN (ONE , NS)
      IF (RH .NE. ZERO) THETA = DATAN2 (CON * X , CON * Y)
      IF (RH.NE.ZERO .OR. NS.GT.ZERO) GO TO 250
      GEOG(2) = - HALFPI
      GO TO 260
  250 CON = ONE / NS
      TS = (RH / (A * F)) ** CON
      GEOG(2) = PHI2Z0 (E,TS)
      IF (IERROR .EQ. 0) GO TO 260
      IERROR = 046
      RETURN
  260 GEOG(1) = ADJLZ0 (THETA / NS + LON0)
      RETURN
!
      END

!                   PJ05Z0
! **********************************************************************
! ** U.S.G.S. GENERAL MAP PROJECTION PACKAGE ...... DR. A. A. ELASSAL **
! **          MATHEMATICAL ANALYSIS BY JOHN SNYDER                    **
! ** MODULE I                VERSION 1.0.2       SEPTEMBER 23, 1983  ***
! **********************************************************************
!                            *  MERCATOR  *
! **********************************************************************
!
      SUBROUTINE PJ05Z0
!
      IMPLICIT REAL*8 (A-Z)
      integer*4 zzone
      INTEGER*4 IERROR,IPEMSG,IPPARM
      INTEGER*4 SWITCH,I,ZONE,ANGS,INFILE
      COMMON /ELLPZ0/ AZ,EZ,ESZ,E0Z,E1Z,E2Z,E3Z
! **** PARAMETERS **** A,E,ES,LON0,X0,Y0,NS,F,RH0,LAT1,M1 **************
      COMMON /ERRMZ0/ IERROR
      COMMON /PRINZ0/ IPEMSG,IPPARM
!     COMMON /WORKZ0/ BUFF(15),ANGS(4,2)
      COMMON /WORKZ0/ BUFF(15)
      COMMON /WK05Z0/ ANGS(4,2)
      real*4 rangs1,rangs2
      equivalence (rangs1,angs(4,1))
      equivalence (rangs2,angs(4,2))
      DIMENSION DATA(1),GEOG(1),PROJ(1)
      DATA HALFPI /1.5707963267948966D0/
      DATA EPSLN /1.0D-10/
      DATA ZERO,ONE /0.0D0,1.0D0/
      DATA SWITCH /0/
!
! ......................................................................
!       .  INITIALIZATION OF PROJECTION PARAMETERS (FILE INPUT)  .
! ......................................................................
!
      ENTRY IF05Z0 (INFILE,DATA)
!
      IERROR = 0
      READ (INFILE,END=160) ZONE,BUFF
      IF (SWITCH.NE.0 .AND. SWITCH.EQ.ZONE) RETURN
  020 IF (BUFF(1) .LE. ZERO) GO TO 100
      A = BUFF(1)
      B = BUFF(2)
      IF (B .GT. ZERO) GO TO 040
      E = ZERO
      ES = ZERO
      GO TO 120
  040 IF (B .GT. ONE) GO TO 060
      E = DSQRT (B)
      ES = B
      GO TO 120
  060 ES = ONE - (B / A) ** 2
      E = DSQRT (ES)
      GO TO 120
  100 A = AZ
      E = EZ
      ES = ESZ
  120 LON0 = PAKRZ0 (BUFF(5))
      LAT1 = PAKRZ0 (BUFF(6))
      M1 = DCOS(LAT1) / (DSQRT( ONE - ES * DSIN(LAT1) **2))
      X0 = BUFF(7)
      Y0 = BUFF(8)
!
! LIST RESULTS OF PARAMETER INITIALIZATION.
!
      CALL RADDZ0 (LAT1,ANGS(1,1))
      CALL RADDZ0 (LON0,ANGS(1,2))
!     IF (IPPARM .EQ. 0) PRINT 2000, A,ES,ANGS,X0,Y0
      IF (IPPARM .EQ. 0) PRINT 2000, A,ES,angs(1,1),angs(2,1),  &
         angs(3,1),rangs1,angs(1,2),angs(2,2),angs(3,2),rangs2, &
         X0,Y0
 2000 FORMAT (' INITIALIZATION PARAMETERS (MERCATOR PROJECTION)'/ &
              ' SEMI-MAJOR AXIS OF ELLIPSOID =',F12.2,' METERS'/  &
              ' ECCENTRICITY SQUARED         =',F12.9/            &
              ' LATITUDE OF TRUE SCALE       = ',A1,2I3,F7.3/     &
              ' CENTRAL LONGITUDE            = ',A1,2I3,F7.3/     &
              ' FALSE EASTING                =',F12.2,' METERS'/  &
              ' FALSE NORTHING               =',F12.2,' METERS')
      DATA(1) = A
      DATA(2) = ES
      SWITCH = ZONE
      RETURN
  160 IF (IPEMSG .EQ. 0) PRINT 2010
 2010 FORMAT (' ERROR PJ05Z0'/                  &
              ' MISSING PROJECTION PARAMETERS')
      IERROR = 051
      RETURN
!
! ......................................................................
!      .  INITIALIZATION OF PROJECTION PARAMETERS (ENTRY INPUT)  .
! ......................................................................
!
      ENTRY IS05Z0 (ZZONE,DATA)
      zone = zzone
!
      IERROR = 0
      IF (SWITCH.NE.0 .AND. SWITCH.EQ.ZONE) RETURN
      DO 180 I = 1,8
      BUFF(I) = DATA(I)
  180 CONTINUE
      GO TO 020
!
! ......................................................................
!                      .  FORWARD TRANSFORMATION  .
! ......................................................................
!
      ENTRY PF05Z0 (GEOG,PROJ)
!
      IERROR = 0
      IF (SWITCH .NE. 0) GO TO 220
      IF (IPEMSG .EQ. 0) PRINT 2010
      IERROR = 052
      RETURN
  220 IF (DABS(DABS(GEOG(2)) - HALFPI) .GT. EPSLN) GO TO 240
      IF (IPEMSG .EQ. 0) PRINT 2020
 2020 FORMAT (' ERROR PJ05Z0'/                                   &
              ' TRANSFORMATION CANNOT BE COMPUTED AT THE POLES')
      IERROR = 053
      RETURN
  240 SINPHI = DSIN (GEOG(2))
      TS = TSFNZ0 (E,GEOG(2),SINPHI)
      PROJ(1) = X0 + A * M1 * ADJLZ0 (GEOG(1) - LON0)
      PROJ(2) = Y0 - A * M1 * DLOG (TS)
      RETURN
!
! ......................................................................
!                      .  INVERSE TRANSFORMATION  .
! ......................................................................
!
      ENTRY PI05Z0 (PROJ,GEOG)
!
      IERROR = 0
      IF (SWITCH .NE. 0) GO TO 260
      IF (IPEMSG .EQ. 0) PRINT 2010
      IERROR = 054
      RETURN
  260 X = PROJ(1) - X0
      Y = PROJ(2) - Y0
      TS = DEXP (- Y / (A * M1))
      GEOG(2) = PHI2Z0 (E,TS)
      IF (IERROR .EQ. 0) GO TO 280
      IERROR = 055
      RETURN
  280 GEOG(1) = ADJLZ0 (LON0 + X / (A * M1))
      RETURN
!
      END

!                   PJ06Z0
! **********************************************************************
! ** U.S.G.S. GENERAL MAP PROJECTION PACKAGE ...... DR. A. A. ELASSAL **
! **          MATHEMATICAL ANALYSIS BY JOHN SNYDER                    **
! ** MODULE I                VERSION 1.0.0            NOVEMBER 1,1980 **
! **********************************************************************
!                       *  POLAR STEREOGRAPHIC  *
! **********************************************************************
!
      SUBROUTINE PJ06Z0
!
      IMPLICIT REAL*8 (A-Z)
      integer*4 zzone
      INTEGER*4 IERROR,IPEMSG,IPPARM
      INTEGER*4 SWITCH,IND,I,ZONE,ANGS,INFILE
      COMMON /ELLPZ0/ AZ,EZ,ESZ,E0Z,E1Z,E2Z,E3Z
! **** PARAMETERS **** A,E,ES,LON0,LATC,X0,Y0,E3,MCS,TCS,FAC,IND *******
      COMMON /ERRMZ0/ IERROR
      COMMON /PRINZ0/ IPEMSG,IPPARM
!     COMMON /WORKZ0/ BUFF(15),ANGS(4,2)
      COMMON /WORKZ0/ BUFF(15)
      COMMON /WK06Z0/ ANGS(4,2)
      real*4 rangs1,rangs2
      equivalence (rangs1,angs(4,1))
      equivalence (rangs2,angs(4,2))
      DIMENSION DATA(1),GEOG(1),PROJ(1)
      DATA NINTYD /90000000.0D0/
      DATA ZERO,ONE,TWO /0.0D0,1.0D0,2.0D0/
      DATA SWITCH /0/
!
! ......................................................................
!       .  INITIALIZATION OF PROJECTION PARAMETERS (FILE INPUT)  .
! ......................................................................
!
      ENTRY IF06Z0 (INFILE,DATA)
!
      IERROR = 0
      READ (INFILE,END=160) ZONE,BUFF
      IF (SWITCH.NE.0 .AND. SWITCH.EQ.ZONE) RETURN
  020 IF (BUFF(1) .LE. ZERO) GO TO 100
      A = BUFF(1)
      B = BUFF(2)
      IF (B .GT. ZERO) GO TO 040
      E = ZERO
      ES = ZERO
      E3 = ONE
      GO TO 120
  040 IF (B .GT. ONE) GO TO 060
      E = DSQRT (B)
      ES = B
      GO TO 080
  060 ES = ONE - (B / A) ** 2
      E = DSQRT (ES)
  080 E3 = E3FNZ0 (E)
      GO TO 120
  100 A = AZ
      E = EZ
      ES = ESZ
      E3 = E3Z
  120 LON0 = PAKRZ0 (BUFF(5))
      SAVE = BUFF(6)
      LATC = PAKRZ0 (SAVE)
      X0 = BUFF(7)
      Y0 = BUFF(8)
      FAC = ONE
      IF (SAVE .LT. ZERO) FAC =-ONE
      IND = 0
      IF (DABS(SAVE) .EQ. NINTYD) GO TO 130
      IND = 1
      CON1 = FAC * LATC
      SINPHI = DSIN (CON1)
      COSPHI = DCOS (CON1)
      MCS = MSFNZ0 (E,SINPHI,COSPHI)
      TCS = TSFNZ0 (E,CON1,SINPHI)
!
! LIST RESULTS OF PARAMETER INITIALIZATION.
!
  130 CALL RADDZ0 (LON0,ANGS(1,1))
      CALL RADDZ0 (LATC,ANGS(1,2))
!     IF (IPPARM .EQ. 0) PRINT 2000, A,ES,ANGS,X0,Y0
      IF (IPPARM .EQ. 0) PRINT 2000, A,ES,angs(1,1),angs(2,1),     &
         angs(3,1),rangs1,angs(1,2),angs(2,2),angs(3,2),rangs2,    &
         X0,Y0
 2000 FORMAT (' INITIALIZATION PARAMETERS (POLAR STEREOGRAPHIC',   &
              ' PROJECTION)'/                                      &
              ' SEMI-MAJOR AXIS OF ELLIPSOID =',F12.2,' METERS'/   &
              ' ECCENTRICITY SQUARED         =',F12.9/             &
              ' LONGITUDE OF Y-AXIS          = ',A1,2I3,F7.3/      &
              ' LATITUDE OF TRUE SCALE       = ',A1,2I3,F7.3/      &
              ' FALSE EASTING                =',F12.2,' METERS'/   &
              ' FALSE NORTHING               =',F12.2,' METERS')
      DATA(1) = A
      DATA(2) = ES
      SWITCH = ZONE
      RETURN
  160 IF (IPEMSG .EQ. 0) PRINT 2010
 2010 FORMAT (' ERROR PJ06Z0'/                   &
              ' MISSING PROJECTION PARAMETERS')
      IERROR = 061
      RETURN
!
! ......................................................................
!      .  INITIALIZATION OF PROJECTION PARAMETERS (ENTRY INPUT)  .
! ......................................................................
!
      ENTRY IS06Z0 (ZZONE,DATA)
      zone = zzone
!
      IERROR = 0
      IF (SWITCH.NE.0 .AND. SWITCH.EQ.ZONE) RETURN
      DO 180 I = 1,8
      BUFF(I) = DATA(I)
  180 CONTINUE
      GO TO 020
!
! ......................................................................
!                      .  FORWARD TRANSFORMATION  .
! ......................................................................
!
      ENTRY PF06Z0 (GEOG,PROJ)
!
      IERROR = 0
      IF (SWITCH .NE. 0) GO TO 220
      IF (IPEMSG .EQ. 0) PRINT 2010
      IERROR = 062
      RETURN
  220 CON1 = FAC * ADJLZ0 (GEOG(1) - LON0)
      CON2 = FAC * GEOG(2)
      SINPHI = DSIN (CON2)
      TS = TSFNZ0 (E,CON2,SINPHI)
      IF (IND .EQ. 0) GO TO 240
      RH = A * MCS * TS / TCS
      GO TO 260
  240 RH = TWO * A * TS / E3
  260 PROJ(1) = X0 + FAC * RH * DSIN (CON1)
      PROJ(2) = Y0 - FAC * RH * DCOS (CON1)
      RETURN
!
! ......................................................................
!                      .  INVERSE TRANSFORMATION  .
! ......................................................................
!
      ENTRY PI06Z0 (PROJ,GEOG)
!
      IERROR = 0
      IF (SWITCH .NE. 0) GO TO 320
      IF (IPEMSG .EQ. 0) PRINT 2010
      IERROR = 063
      RETURN
  320 X = FAC * (PROJ(1) - X0)
      Y = FAC * (PROJ(2) - Y0)
      RH = DSQRT (X * X + Y * Y)
      IF (IND .EQ. 0) GO TO 340
      TS = RH * TCS / (A * MCS)
      GO TO 360
  340 TS = RH * E3 / (TWO * A)
  360 GEOG(2) = FAC * PHI2Z0 (E,TS)
      IF (IERROR .EQ. 0) GO TO 380
      IERROR = 064
      RETURN
  380 IF (RH .NE. ZERO) GO TO 400
      GEOG(1) = FAC * LON0
      RETURN
  400 GEOG(1) = ADJLZ0 (FAC * DATAN2 (X , -Y) + LON0)
      RETURN
!
      END

!                   PJ07Z0
! **********************************************************************
! ** U.S.G.S. GENERAL MAP PROJECTION PACKAGE ...... DR. A. A. ELASSAL **
! **          MATHEMATICAL ANALYSIS BY JOHN SNYDER                    **
! ** MODULE I                VERSION 1.0.0            NOVEMBER 1,1980 **
! **********************************************************************
!                            *  POLYCONIC  *
! **********************************************************************
!
      SUBROUTINE PJ07Z0
!
      IMPLICIT REAL*8 (A-Z)
      integer*4 zzone
      INTEGER*4 IERROR,IPEMSG,IPPARM
      INTEGER*4 SWITCH,I,ZONE,ANGS,INFILE
      COMMON /ELLPZ0/ AZ,EZ,ESZ,E0Z,E1Z,E2Z,E3Z
! **** PARAMETERS **** A,E,ES,LON0,LAT0,X0,Y0,E0,E1,E2,ML0 *************
      COMMON /ERRMZ0/ IERROR
      COMMON /PRINZ0/ IPEMSG,IPPARM
!     COMMON /WORKZ0/ BUFF(15),ANGS(4,2)
      COMMON /WORKZ0/ BUFF(15)
      COMMON /WK07Z0/ ANGS(4,2)
      real*4 rangs1,rangs2
      equivalence (rangs1,angs(4,1))
      equivalence (rangs2,angs(4,2))
      DIMENSION DATA(1),GEOG(1),PROJ(1)
      DATA TOL /1.0D-7/
      DATA ZERO,ONE /0.0D0,1.0D0/
      DATA SWITCH /0/
!
! ......................................................................
!       .  INITIALIZATION OF PROJECTION PARAMETERS (FILE INPUT)  .
! ......................................................................
!
      ENTRY IF07Z0 (INFILE,DATA)
!
      IERROR = 0
      READ (INFILE,END=160) ZONE,BUFF
      IF (SWITCH.NE.0 .AND. SWITCH.EQ.ZONE) RETURN
  020 IF (BUFF(1) .LE. ZERO) GO TO 100
      A = BUFF(1)
      B = BUFF(2)
      IF (B .GT. ZERO) GO TO 040
      E = ZERO
      ES = ZERO
      E0 = ONE
      E1 = ZERO
      E2 = ZERO
      GO TO 120
  040 IF (B .GT. ONE) GO TO 060
      E = DSQRT (B)
      ES = B
      GO TO 080
  060 ES = ONE - (B / A) ** 2
      E = DSQRT (ES)
  080 E0 = E0FNZ0 (ES)
      E1 = E1FNZ0 (ES)
      E2 = E2FNZ0 (ES)
      GO TO 120
  100 A = AZ
      E = EZ
      ES = ESZ
      E0 = E0Z
      E1 = E1Z
      E2 = E2Z
  120 LON0 = PAKRZ0 (BUFF(5))
      LAT0 = PAKRZ0 (BUFF(6))
      X0 = BUFF(7)
      Y0 = BUFF(8)
      ML0 = MLFNZ0 (E0,E1,E2,LAT0)
!
! LIST RESULTS OF PARAMETER INITIALIZATION.
!
      CALL RADDZ0 (LON0,ANGS(1,1))
      CALL RADDZ0 (LAT0,ANGS(1,2))
!     IF (IPPARM .EQ. 0) PRINT 2000, A,ES,ANGS,X0,Y0
      IF (IPPARM .EQ. 0) PRINT 2000, A,ES,angs(1,1),angs(2,1),   &
         angs(3,1),rangs1,angs(1,2),angs(2,2),angs(3,2),rangs2,  &
         X0,Y0
 2000 FORMAT (' INITIALIZATION PARAMETERS (POLYCONIC',            &
              ' PROJECTION)'/                                     &
              ' SEMI-MAJOR AXIS OF ELLIPSOID =',F12.2,' METERS'/  &
              ' ECCENTRICITY SQUARED         =',F12.9/            &
              ' LONGITUDE OF ORIGIN          = ',A1,2I3,F7.3/     &
              ' LATITUDE OF ORIGIN           = ',A1,2I3,F7.3/     &
              ' FALSE EASTING                =',F12.2,' METERS'/  &
              ' FALSE NORTHING               =',F12.2,' METERS')
      DATA(1) = A
      DATA(2) = ES
      SWITCH = ZONE
      RETURN
  160 IF (IPEMSG .EQ. 0) PRINT 2010
 2010 FORMAT (' ERROR PJ07Z0'/                   &
              ' MISSING PROJECTION PARAMETERS')
      IERROR = 071
      RETURN
!
! ......................................................................
!      .  INITIALIZATION OF PROJECTION PARAMETERS (ENTRY INPUT)  .
! ......................................................................
!
      ENTRY IS07Z0 (ZZONE,DATA)
      zone = zzone
!
      IERROR = 0
      IF (SWITCH.NE.0 .AND. SWITCH.EQ.ZONE) RETURN
      DO 180 I = 1,8
      BUFF(I) = DATA(I)
  180 CONTINUE
      GO TO 020
!
! ......................................................................
!                      .  FORWARD TRANSFORMATION  .
! ......................................................................
!
      ENTRY PF07Z0 (GEOG,PROJ)
!
      IERROR = 0
      IF (SWITCH .NE. 0) GO TO 220
      IF (IPEMSG .EQ. 0) PRINT 2010
      IERROR = 072
      RETURN
  220 CON = ADJLZ0 (GEOG(1) - LON0)
      IF (DABS(GEOG(2)) .GT. TOL) GO TO 240
      PROJ(1) = X0 + A * CON
      PROJ(2) = Y0 - A * ML0
      RETURN
  240 SINPHI = DSIN (GEOG(2))
      COSPHI = DCOS (GEOG(2))
      ML = MLFNZ0 (E0,E1,E2,GEOG(2))
      MS = MSFNZ0 (E,SINPHI,COSPHI)
      CON = CON * SINPHI
      PROJ(1) = X0 + A * MS * DSIN (CON) / SINPHI
      PROJ(2) = Y0 + A * (ML - ML0 + MS * (ONE - DCOS (CON)) / SINPHI)
      RETURN
!
! ......................................................................
!                      .  INVERSE TRANSFORMATION  .
! ......................................................................
!
      ENTRY PI07Z0 (PROJ,GEOG)
!
      IERROR = 0
      IF (SWITCH .NE. 0) GO TO 320
      IF (IPEMSG .EQ. 0) PRINT 2010
      IERROR = 073
      RETURN
  320 X = PROJ(1) - X0
      Y = PROJ(2) - Y0
      AL = ML0 + Y / A
      IF (DABS (AL) .GT. TOL) GO TO 340
      GEOG(1) = X / A + LON0
      GEOG(2) = ZERO
      RETURN
  340 B = AL * AL + (X / A) ** 2
      GEOG(2) = PHI4Z0 (ES,E0,E1,E2,AL,B,C)
      IF (IERROR .EQ. 0) GO TO 360
      IERROR = 074
      RETURN
  360 GEOG(1) = ADJLZ0 ( DASIN (X * C / A) / DSIN (GEOG(2)) + LON0)
      RETURN
!
      END

!                   PJ08Z0
! **********************************************************************
! ** U.S.G.S. GENERAL MAP PROJECTION PACKAGE ...... DR. A. A. ELASSAL **
! **          MATHEMATICAL ANALYSIS BY JOHN SNYDER                    **
! ** MODULE I                VERSION 1.0.0            NOVEMBER 1,1980 **
! **********************************************************************
!                        *  EQUIDISTANT CONIC  *
! **********************************************************************
!
      SUBROUTINE PJ08Z0
!
      IMPLICIT REAL*8 (A-Z)
      integer*4 zzone
      INTEGER*4 IERROR,IPEMSG,IPPARM
      INTEGER*4 SWITCH,IND,I,ZONE,ANGS,ANG1,ANG2,INFILE
      COMMON /ELLPZ0/ AZ,EZ,ESZ,E0Z,E1Z,E2Z,E3Z
! ** PARAMETERS ** A,E,ES,LAT1,LAT2,LON0,LAT0,X0,Y0,E0,E1,E2,NS,GL,RH0 *
      COMMON /ERRMZ0/ IERROR
      COMMON /PRINZ0/ IPEMSG,IPPARM
!     COMMON /WORKZ0/ BUFF(15),ANGS(4,4)
      COMMON /WORKZ0/ BUFF(15)
      COMMON /WK08Z0/ ANGS(4,4)
      real*4 rangs1,rangs2,rangs3,rangs4
      equivalence (rangs1,angs(4,1))
      equivalence (rangs2,angs(4,2))
      equivalence (rangs3,angs(4,3))
      equivalence (rangs4,angs(4,4))
      DIMENSION DATA(1),GEOG(1),PROJ(1),ANG1(4),ANG2(4,2)
      EQUIVALENCE (ANG1(1),ANGS(1,1)) , (ANG2(1,1),ANGS(1,3))
      DATA ZERO,ONE /0.0D0,1.0D0/
      DATA EPSLN /1.0D-10/
      DATA SWITCH /0/
!
! ......................................................................
!       .  INITIALIZATION OF PROJECTION PARAMETERS (FILE INPUT)  .
! ......................................................................
!
      ENTRY IF08Z0 (INFILE,DATA)
!
      IERROR = 0
      READ (INFILE,END=240) ZONE,BUFF
      IF (SWITCH.NE.0 .AND. SWITCH.EQ.ZONE) RETURN
  020 IF (BUFF(1) .LE. ZERO) GO TO 100
      A = BUFF(1)
      B = BUFF(2)
      IF (B .GT. ZERO) GO TO 040
      E = ZERO
      ES = ZERO
      E0 = ONE
      E1 = ZERO
      E2 = ZERO
      GO TO 120
  040 IF (B .GT. ONE) GO TO 060
      E = DSQRT (B)
      ES = B
      GO TO 080
  060 ES = ONE - (B / A) ** 2
      E = DSQRT (ES)
  080 E0 = E0FNZ0 (ES)
      E1 = E1FNZ0 (ES)
      E2 = E2FNZ0 (ES)
      GO TO 120
  100 A = AZ
      E = EZ
      ES = ESZ
      E0 = E0Z
      E1 = E1Z
      E2 = E2Z
  120 LAT1 = PAKRZ0 (BUFF(3))
      LAT2 = PAKRZ0 (BUFF(4))
      IF (DABS(LAT1+LAT2) .GE. EPSLN) GO TO 130
      IF (IPEMSG .EQ. 0) PRINT 2000
 2000 FORMAT (' ERROR PJ08Z0'/                                   &
              ' EQUAL LATITUDES FOR ST. PARALLELS ON OPPOSITE',  &
              ' SIDES OF EQUATOR')
      IERROR = 081
      RETURN
  130 LON0 = PAKRZ0 (BUFF(5))
      LAT0 = PAKRZ0 (BUFF(6))
      X0 = BUFF(7)
      Y0 = BUFF(8)
      SINPHI = DSIN (LAT1)
      COSPHI = DCOS (LAT1)
      MS1 = MSFNZ0 (E,SINPHI,COSPHI)
      ML1 = MLFNZ0 (E0,E1,E2,LAT1)
      IND = 0
      IF (BUFF(9) .NE. ZERO) GO TO 140
      NS = SINPHI
      GO TO 160
  140 IND = 1
      SINPHI = DSIN (LAT2)
      COSPHI = DCOS (LAT2)
      MS2 = MSFNZ0 (E,SINPHI,COSPHI)
      ML2 = MLFNZ0 (E0,E1,E2,LAT2)
      IF (DABS(LAT1-LAT2) .GE. EPSLN) GO TO 150
      NS = SINPHI
      GO TO 160
  150 NS = (MS1 - MS2) / (ML2 - ML1)
  160 GL = ML1 + MS1 / NS
      ML0 = MLFNZ0 (E0,E1,E2,LAT0)
      RH0 = A * (GL - ML0)
!
! LIST RESULTS OF PARAMETER INITIALIZATION.
!
      CALL RADDZ0 (LAT1,ANGS(1,1))
      CALL RADDZ0 (LAT2,ANGS(1,2))
      CALL RADDZ0 (LON0,ANGS(1,3))
      CALL RADDZ0 (LAT0,ANGS(1,4))
      IF (IND .EQ. 0) GO TO 200
!     IF (IPPARM .EQ. 0) PRINT 2010, A,ES,ANGS,X0,Y0
      IF (IPPARM .EQ. 0) PRINT 2010, A,ES,angs(1,1),angs(2,1),   &
         angs(3,1),rangs1,angs(1,2),angs(2,2),angs(3,2),rangs2,  &
         angs(1,3),angs(2,3),angs(3,3),rangs3,angs(1,4),         &
         angs(2,4),angs(3,4),rangs4,X0,Y0
 2010 FORMAT (' INITIALIZATION PARAMETERS (EQUIDISTANT CONIC',    &
              ' PROJECTION)'/                                     &
              ' SEMI-MAJOR AXIS OF ELLIPSOID =',F12.2,' METERS'/  &
              ' ECCENTRICITY SQUARED         =',F12.9/            & 
              ' LATITUDE OF 1ST ST. PARALLEL = ',A1,2I3,F7.3/     &
              ' LATITUDE OF 2ND ST. PARALLEL = ',A1,2I3,F7.3/     &
              ' LONGITUDE OF ORIGIN          = ',A1,2I3,F7.3/     &
              ' LATITUDE OF ORIGIN           = ',A1,2I3,F7.3/     &
              ' FALSE EASTING                =',F12.2,' METERS'/  &
              ' FALSE NORTHING               =',F12.2,' METERS')
      GO TO 220
! 200 IF (IPPARM .EQ. 0) PRINT 2020, A,ES,ANG1,ANG2,X0,Y0
  200 IF (IPPARM .EQ. 0) PRINT 2020, A,ES,angs(1,1),angs(2,1),angs(3,1), &
          rangs1,angs(1,3),angs(2,3),angs(3,3),rangs3,angs(1,4),         &
          angs(2,4),angs(3,4),rangs4,X0,Y0
 2020 FORMAT (' INITIALIZATION PARAMETERS (EQUIDISTANT CONIC',    &
              ' PROJECTION)'/                                     &
              ' SEMI-MAJOR AXIS OF ELLIPSOID =',F12.2,' METERS'/  &
              ' ECCENTRICITY SQUARED         =',F12.9/            &
              ' LATITUDE OF ST. PARALLEL     = ',A1,2I3,F7.3/     & 
              ' LONGITUDE OF ORIGIN          = ',A1,2I3,F7.3/     & 
              ' LATITUDE OF ORIGIN           = ',A1,2I3,F7.3/     &
              ' FALSE EASTING                =',F12.2,' METERS'/  &
              ' FALSE NORTHING               =',F12.2,' METERS')
  220 DATA(1) = A
      DATA(2) = ES
      SWITCH = ZONE
      RETURN
  240 IF (IPEMSG .EQ. 0) PRINT 2030
 2030 FORMAT (' ERROR PJ08Z0'/                    &
              ' MISSING PROJECTION PARAMETERS')
      IERROR = 082
      RETURN
!
! ......................................................................
!      .  INITIALIZATION OF PROJECTION PARAMETERS (ENTRY INPUT)  .
! ......................................................................
!
      ENTRY IS08Z0 (ZZONE,DATA)
      zone = zzone
!
      IERROR = 0
      IF (SWITCH.NE.0 .AND. SWITCH.EQ.ZONE) RETURN
      DO 260 I = 1,9
      BUFF(I) = DATA(I)
  260 CONTINUE
      GO TO 020
!
! ......................................................................
!                      .  FORWARD TRANSFORMATION  .
! ......................................................................
!
      ENTRY PF08Z0 (GEOG,PROJ)
!
      IERROR = 0
      IF (SWITCH .NE. 0) GO TO 300
      IF (IPEMSG .EQ. 0) PRINT 2030
      IERROR = 083
      RETURN
  300 ML = MLFNZ0 (E0,E1,E2,GEOG(2))
      RH = A * (GL - ML)
      THETA = NS * ADJLZ0 (GEOG(1) - LON0)
      PROJ(1) = X0 + RH * DSIN (THETA)
      PROJ(2) = Y0 + RH0 - RH * DCOS (THETA)
      RETURN
!
! ......................................................................
!                      .  INVERSE TRANSFORMATION  .
! ......................................................................
!
      ENTRY PI08Z0 (PROJ,GEOG)
!
      IERROR = 0
      IF (SWITCH .NE. 0) GO TO 320
      IF (IPEMSG .EQ. 0) PRINT 2030
      IERROR = 084
      RETURN
  320 X = PROJ(1) - X0
      Y = RH0 - PROJ(2) + Y0
      RH = DSIGN (DSQRT (X * X + Y * Y) , NS)
      THETA = ZERO
      CON = DSIGN (ONE , NS)
      IF (RH .NE. ZERO) THETA = DATAN2 (CON * X , CON * Y)
      ML = GL - RH / A
      GEOG(2) = PHI3Z0 (ML,E0,E1,E2)
      IF (IERROR .EQ. 0) GO TO 340
      IERROR = 085
      RETURN
  340 GEOG(1) = ADJLZ0 (LON0 + THETA / NS)
      RETURN
!
      END

!                   PJ09Z0
! **********************************************************************
! ** U.S.G.S. GENERAL MAP PROJECTION PACKAGE ...... DR. A. A. ELASSAL **
! **          MATHEMATICAL ANALYSIS BY JOHN SNYDER                    **
! ** MODULE I                VERSION 1.0.2                MAY 14,1981 **
! **********************************************************************
!                       *  TRANSVERSE MARCATOR  *
! **********************************************************************
!
      SUBROUTINE PJ09Z0
!
      IMPLICIT REAL*8 (A-Z)
      integer*4 zzone
      INTEGER*4 IERROR,IPEMSG,IPPARM
      INTEGER*4 SWITCH,I,ZONE,ANGS,INFILE,IND,NIT
      COMMON /ELLPZ0/ AZ,EZ,ESZ,E0Z,E1Z,E2Z,E3Z
! **** PARAMETERS **** A,E,ES,KS0,LON0,LAT0,X0,Y0,E0,E1,E2,ESP,ML0,IND *
      COMMON /ERRMZ0/ IERROR
      COMMON /PRINZ0/ IPEMSG,IPPARM
!     COMMON /WORKZ0/ BUFF(15),ANGS(4,2)
      COMMON /WORKZ0/ BUFF(15)
      COMMON /WK09Z0/ ANGS(4,2)
      real*4 rangs1,rangs2
      equivalence (rangs1,angs(4,1))
      equivalence (rangs2,angs(4,2))
      DIMENSION DATA(15),GEOG(1),PROJ(1)
      DATA ZERO,HALF,ONE,TWO,THREE /0.0D0,0.5D0,1.0D0,2.0D0,3.0D0/
      DATA FOUR,FIVE,SIX,EIGHT,NINE /4.0D0,5.0D0,6.0D0,8.0D0,9.0D0/
      DATA HALFPI /1.5707963267948966D0/
      DATA TEN /10.0D0/
      DATA TOL,EPSLN,NIT /1.0D-5,1.0D-10,6/
      DATA SWITCH /0/
!
! ......................................................................
!       .  INITIALIZATION OF PROJECTION PARAMETERS (FILE INPUT)  .
! ......................................................................
!
      ENTRY IF09Z0 (INFILE,data)
!
      IERROR = 0
      READ (INFILE,END=160) ZONE,BUFF
      IF (SWITCH.NE.0 .AND. SWITCH.EQ.ZONE) RETURN
  020 IF (BUFF(1) .LE. ZERO) GO TO 100
      A = BUFF(1)
      B = BUFF(2)
      IF (B .GT. ZERO) GO TO 040
      E = ZERO
      ES = ZERO
      E0 = ONE
      E1 = ZERO
      E2 = ZERO
      GO TO 120
  040 IF (B .GT. ONE) GO TO 060
      E = DSQRT (B)
      ES = B
      GO TO 080
  060 ES = ONE - (B / A) ** 2
      E = DSQRT (ES)
  080 E0 = E0FNZ0 (ES)
      E1 = E1FNZ0 (ES)
      E2 = E2FNZ0 (ES)
      GO TO 120
  100 A = AZ
      E = EZ
      ES = ESZ
      E0 = E0Z
      E1 = E1Z
      E2 = E2Z
  120 KS0 = BUFF(3)
      LON0 = PAKRZ0 (BUFF(5))
      LAT0 = PAKRZ0 (BUFF(6))
      X0 = BUFF(7)
      Y0 = BUFF(8)
      ML0 = A * MLFNZ0 (E0,E1,E2,LAT0)
      IND = 1
      IF (E .LT. TOL) GO TO 130
      IND = 0
      ESP = ES / (ONE - ES)
!
! LIST RESULTS OF PARAMETER INITIALIZATION.
!
  130 CALL RADDZ0 (LON0,ANGS(1,1))
      CALL RADDZ0 (LAT0,ANGS(1,2))
!     IF (IPPARM .EQ. 0) PRINT 2000, A,ES,KS0,ANGS,X0,Y0
      IF (IPPARM .EQ. 0) PRINT 2000, A,ES,ks0,angs(1,1),angs(2,1),angs(3,1),  &
          rangs1,angs(1,2),angs(2,2),angs(3,2),rangs2, X0,Y0
 2000 FORMAT (' INITIALIZATION PARAMETERS (TRANSVERSE MERCATOR',   &
              ' PROJECTION)'/                                      &
              ' SEMI-MAJOR AXIS OF ELLIPSOID =',F12.2,' METERS'/   &
              ' ECCENTRICITY SQUARED         =',F12.9/             &
              ' SCALE FACTOR AT C. MERIDIAN  =',F9.6/              &
              ' LONGITUDE OF C. MERIDIAN     = ',A1,2I3,F7.3/      &  
              ' LATITUDE OF ORIGIN           = ',A1,2I3,F7.3/      &
              ' FALSE EASTING                =',F12.2,' METERS'/   &
              ' FALSE NORTHING               =',F12.2,' METERS')
      DATA(1) = A
      DATA(2) = ES
      SWITCH = ZONE
      RETURN
  160 IF (IPEMSG .EQ. 0) PRINT 2010
 2010 FORMAT (' ERROR PJ09Z0'/                      &
              ' MISSING PROJECTION PARAMETERS')
      IERROR = 091
      RETURN
!
! ......................................................................
!      .  INITIALIZATION OF PROJECTION PARAMETERS (ENTRY INPUT)  .
! ......................................................................
!
      ENTRY IS09Z0 (ZZONE,DATA)
      zone = zzone
!
      IERROR = 0
!$$$$$$$$$$$$$$$$$$$$$ ADDITIONS BY JFWAANANEN 5/1/81 $$$$$$$$$$$
      IF (DATA(1).NE.0.0D0.AND.DATA(1).NE.BUFF(1)) SWITCH=0
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      IF (SWITCH.NE.0 .AND. SWITCH.EQ.ZONE) RETURN
      DO 180 I = 1,8
      BUFF(I) = DATA(I)
  180 CONTINUE
      GO TO 020
!
! ......................................................................
!                      .  FORWARD TRANSFORMATION  .
! ......................................................................
!
      ENTRY PF09Z0 (GEOG,PROJ)
!
      IERROR = 0
      IF (SWITCH .NE. 0) GO TO 220
      IF (IPEMSG .EQ. 0) PRINT 2010
      IERROR = 092
      RETURN
  220 DLON = ADJLZ0 (GEOG(1) - LON0)
      LAT = GEOG(2)
      IF (IND .EQ. 0) GO TO 240
      COSPHI = DCOS (LAT)
      B = COSPHI * DSIN (DLON)
      IF (DABS(DABS(B) - ONE) .GT. EPSLN) GO TO 230
      IF (IPEMSG .EQ. 0) PRINT 2020
 2020 FORMAT (' ERROR PJ09Z0'/                    &
              ' POINT PROJECTS INTO INFINITY')
      IERROR = 093
      RETURN
  230 PROJ(1) = HALF * A * KS0 * DLOG ((ONE + B) / (ONE - B))
      CON =  DACOS (COSPHI * DCOS (DLON) / DSQRT (ONE - B * B))
      IF (LAT .LT. ZERO) CON =-CON
      PROJ(2) = A * KS0 * (CON - LAT0)
      RETURN
!
  240 SINPHI = DSIN (LAT)
      COSPHI = DCOS (LAT)
      AL = COSPHI * DLON
      ALS = AL * AL
      C = ESP * COSPHI * COSPHI
      TQ = DTAN (LAT)
      T = TQ * TQ
      N = A / DSQRT (ONE - ES * SINPHI * SINPHI)
      ML = A * MLFNZ0 (E0,E1,E2,LAT)
      PROJ(1) = KS0 * N * AL * (ONE + ALS / SIX * (ONE - T + C +          &
                ALS / 20.0D0 * (FIVE - 18.0D0 * T + T * T + 72.0D0 *      &
                C - 58.0D0 * ESP))) + X0
      PROJ(2) = KS0 * (ML - ML0 + N * TQ * (ALS * (HALF + ALS / 24.0D0 *  &
                (FIVE - T + NINE * C + FOUR * C * C + ALS / 30.0D0 *      &
                (61.0D0 - 58.0D0 * T + T * T + 600.0D0 * C -              &
                330.0D0 * ESP))))) + Y0
      RETURN
!
! ......................................................................
!                      .  INVERSE TRANSFORMATION  .
! ......................................................................
!
      ENTRY PI09Z0 (PROJ,GEOG)
!
      IERROR = 0
      IF (SWITCH .NE. 0) GO TO 320
      IF (IPEMSG .EQ. 0) PRINT 2010
      IERROR = 094
      RETURN
  320 X = PROJ(1) - X0
      Y = PROJ(2) - Y0
      IF (IND .EQ. 0) GO TO 340
      F = DEXP (X / (A * KS0))
      G = HALF * (F - ONE / F)
      H = DCOS (LAT0 + Y / (A * KS0))
      CON = DSQRT ((ONE - H * H) / (ONE + G * G))
      GEOG(2) =  DASIN (CON)
      IF (Y .LT. ZERO) GEOG(2) =-GEOG(2)
      IF (G.NE.ZERO .OR. H.NE.ZERO) GO TO 330
      GEOG(1) = LON0
      RETURN
  330 GEOG(1) = ADJLZ0 (DATAN2 (G,H) + LON0)
      RETURN
!
  340 CON = (ML0 + Y / KS0) / A
      PHI = CON
      DO 360 I = 1,NIT
      DPHI = ((CON + E1 * DSIN (TWO * PHI) - E2 * DSIN (FOUR * PHI)) / E0) - PHI
      PHI = PHI + DPHI
      IF (DABS(DPHI) .LE. EPSLN) GO TO 380
  360 CONTINUE
      IF (IPEMSG .EQ. 0) PRINT 2030, NIT
 2030 FORMAT (' ERROR PI09Z0' /                                       &
              ' LATITUDE FAILED TO CONVERGE AFTER',I3,' ITERATIONS')
      IERROR = 095
      RETURN
  380 IF (DABS(PHI) .LT. HALFPI) GO TO 400
      GEOG(2) = DSIGN (HALFPI , Y)
      GEOG(1) = LON0
      RETURN
  400 SINPHI = DSIN (PHI)
      COSPHI = DCOS (PHI)
      TANPHI = DTAN (PHI)
      C = ESP * COSPHI * COSPHI
      CS = C * C
      T = TANPHI * TANPHI
      TS = T * T
      CON = ONE - ES * SINPHI * SINPHI
      N = A / DSQRT (CON)
      R = N * (ONE - ES) / CON
      D = X / (N * KS0)
      DS = D * D
      GEOG(2) = PHI - (N * TANPHI * DS / R) * (HALF - DS / 24.0D0 *     &
                (FIVE + THREE * T + TEN * C - FOUR * CS - NINE * ESP -  &
                DS / 30.0D0 * (61.0D0 + 90.0D0 * T + 298.0D0 * C +      &
                45.0D0 * TS - 252.0D0 * ESP - THREE * CS)))
      GEOG(1) = ADJLZ0 (LON0 + (D * (ONE - DS / SIX * (ONE + TWO *      &
                T + C - DS / 20.0D0 * (FIVE - TWO * C + 28.0D0 * T -    &
                THREE * CS + EIGHT * ESP + 24.0D0 * TS))) / COSPHI))
      RETURN
!
      END

!                   PJ10Z0
! **********************************************************************
! ** U.S.G.S. GENERAL MAP PROJECTION PACKAGE ...... DR. A. A. ELASSAL **
! **          MATHEMATICAL ANALYSIS BY JOHN SNYDER                    **
! ** MODULE I                VERSION 1.0.0            NOVEMBER 1,1980 **
! **********************************************************************
!                          *  STEREOGRAPHIC  *
! **********************************************************************
!
      SUBROUTINE PJ10Z0
!
      IMPLICIT REAL*8 (A-Z)
      integer*4 zzone
      INTEGER*4 IERROR,IPEMSG,IPPARM
      INTEGER*4 SWITCH,I,ZONE,ANGS,INFILE
      COMMON /SPHRZ0/ AZZ
! **** PARAMETERS **** A,LON0,LAT0,X0,Y0,SINPH0,COSPH0 *****************
      COMMON /ERRMZ0/ IERROR
      COMMON /PRINZ0/ IPEMSG,IPPARM
!     COMMON /WORKZ0/ BUFF(15),ANGS(4,2)
      COMMON /WORKZ0/ BUFF(15)
      COMMON /WK10Z0/ ANGS(4,2)
      real*4 rangs1,rangs2
      equivalence (rangs1,angs(4,1))
      equivalence (rangs2,angs(4,2))
      DIMENSION DATA(1),GEOG(1),PROJ(1)
      DATA HALFPI /1.5707963267948966D0/
      DATA EPSLN /1.0D-10/
      DATA ZERO,ONE,TWO /0.0D0,1.0D0,2.0D0/
      DATA SWITCH /0/
!
! ......................................................................
!       .  INITIALIZATION OF PROJECTION PARAMETERS (FILE INPUT)  .
! ......................................................................
!
      ENTRY IF10Z0 (INFILE,data)
!
      IERROR = 0
      READ (INFILE,END=060) ZONE,BUFF
      IF (SWITCH.NE.0 .AND. SWITCH.EQ.ZONE) RETURN
  020 A = BUFF(1)
      IF (A .LE. ZERO) A = AZZ
      LON0 = PAKRZ0 (BUFF(5))
      LAT0 = PAKRZ0 (BUFF(6))
      X0 = BUFF(7)
      Y0 = BUFF(8)
      SINPH0 = DSIN (LAT0)
      COSPH0 = DCOS (LAT0)
!
! LIST RESULTS OF PARAMETER INITIALIZATION.
!
      CALL RADDZ0 (LON0,ANGS(1,1))
      CALL RADDZ0 (LAT0,ANGS(1,2))
!     IF (IPPARM .EQ. 0) PRINT 2000, A,ANGS,X0,Y0
      IF (IPPARM .EQ. 0) PRINT 2000, A,angs(1,1),angs(2,1),angs(3,1), &
          rangs1,angs(1,2),angs(2,2),angs(3,2),rangs2,X0,Y0
 2000 FORMAT (' INITIALIZATION PARAMETERS (STEREOGRAPHIC',         &
              ' PROJECTION)'/                                      &
              ' RADIUS OF SPHERE             =',F12.2,' METERS'/   &
              ' LONGITUDE OF CENTER          = ',A1,2I3,F7.3/      &
              ' LATITUDE  OF CENTER          = ',A1,2I3,F7.3/      &
              ' FALSE EASTING                =',F12.2,' METERS'/   &
              ' FALSE NORTHING               =',F12.2,' METERS')
      DATA(1) = A
      SWITCH = ZONE
      RETURN
  060 IF (IPEMSG .EQ. 0) PRINT 2010
 2010 FORMAT (' ERROR PJ10Z0'/                   &
              ' MISSING PROJECTION PARAMETERS')
      IERROR = 101
      RETURN
!
! ......................................................................
!      .  INITIALIZATION OF PROJECTION PARAMETERS (ENTRY INPUT)  .
! ......................................................................
!
      ENTRY IS10Z0 (ZZONE,DATA)
      zone = zzone
!
      IERROR = 0
      IF (SWITCH.NE.0 .AND. SWITCH.EQ.ZONE) RETURN
      DO 180 I = 1,8
      BUFF(I) = DATA(I)
  180 CONTINUE
      GO TO 020
!
! ......................................................................
!                      .  FORWARD TRANSFORMATION  .
! ......................................................................
!
      ENTRY PF10Z0 (GEOG,PROJ)
!
      IERROR = 0
      IF (SWITCH .NE. 0) GO TO 120
      IF (IPEMSG .EQ. 0) PRINT 2010
      IERROR = 102
      RETURN
  120 LON = ADJLZ0 (GEOG(1) - LON0)
      SINPHI = DSIN (GEOG(2))
      COSPHI = DCOS (GEOG(2))
      COSLON = DCOS (LON)
      G = SINPH0 * SINPHI + COSPH0 * COSPHI * COSLON
      IF (DABS(G + ONE) .GT. EPSLN) GO TO 140
      if (ipemsg .eq. 0) PRINT 2020
 2020 FORMAT (' ERROR PJ10Z0'/                  &
              ' POINT PROJECTS INTO INFINITY')
      IERROR = 103
      RETURN
  140 KSP = TWO / (ONE + G)
      PROJ(1) = X0 + A * KSP * COSPHI * DSIN (LON)
      PROJ(2) = Y0 + A * KSP * (COSPH0 * SINPHI - SINPH0 * COSPHI * COSLON)
      RETURN
!
! ......................................................................
!                      .  INVERSE TRANSFORMATION  .
! ......................................................................
!
      ENTRY PI10Z0 (PROJ,GEOG)
!
      IERROR = 0
      IF (SWITCH .NE. 0) GO TO 220
      IF (IPEMSG .EQ. 0) PRINT 2010
      IERROR = 104
      RETURN
  220 X = PROJ(1) - X0
      Y = PROJ(2) - Y0
      RH = DSQRT (X * X + Y * Y)
      Z = TWO * DATAN (RH / (TWO * A))
      SINZ = DSIN (Z)
      COSZ = DCOS (Z)
      GEOG(1) = LON0
      IF (DABS(RH) .GT. EPSLN) GO TO 240
      GEOG(2) = LAT0
      RETURN
  240 GEOG(2) =  DASIN (COSZ * SINPH0 + Y * SINZ * COSPH0 / RH)
      CON = DABS (LAT0) - HALFPI
      IF (DABS (CON) .GT. EPSLN) GO TO 260
      IF (LAT0 .LT. ZERO) GO TO 250
      GEOG(1) = ADJLZ0 (LON0 + DATAN2 (X , -Y))
      RETURN
  250 GEOG(1) = ADJLZ0 (LON0 - DATAN2 (-X , Y))
      RETURN
  260 CON = COSZ - SINPH0 * DSIN (GEOG(2))
      IF (CON.EQ.ZERO .AND. X.EQ.ZERO) RETURN
      GEOG(1) = ADJLZ0 (LON0 + DATAN2 ((X*SINZ*COSPH0) , (CON*RH)))
      RETURN
!
      END

!                   PJ11Z0
! **********************************************************************
! ** U.S.G.S. GENERAL MAP PROJECTION PACKAGE ...... DR. A. A. ELASSAL **
! **          MATHEMATICAL ANALYSIS BY JOHN SNYDER                    **
! ** MODULE I                VERSION 1.0.0            NOVEMBER 1,1980 **
! **********************************************************************
!                   *  LAMBERT AZIMUTHAL EQUAL-AREA  *
! **********************************************************************
!
      SUBROUTINE PJ11Z0
!
      IMPLICIT REAL*8 (A-Z)
      integer*4 zzone
      INTEGER*4 IERROR,IPEMSG,IPPARM
      INTEGER*4 SWITCH,I,ZONE,ANGS,INFILE
      COMMON /SPHRZ0/ AZZ
! **** PARAMETERS **** A,LON0,LAT0,X0,Y0,SINPH0,COSPH0 *****************
      COMMON /ERRMZ0/ IERROR
      COMMON /PRINZ0/ IPEMSG,IPPARM
!     COMMON /WORKZ0/ BUFF(15),ANGS(4,2)
      COMMON /WORKZ0/ BUFF(15)
      COMMON /WK11Z0/ ANGS(4,2)
      real*4 rangs1,rangs2
      equivalence (rangs1,angs(4,1))
      equivalence (rangs2,angs(4,2))
      DIMENSION DATA(1),GEOG(1),PROJ(1)
      DATA HALFPI /1.5707963267948966D0/
      DATA EPSLN /1.0D-10/
      DATA ZERO,ONE,TWO /0.0D0,1.0D0,2.0D0/
      DATA SWITCH /0/
!
! ......................................................................
!       .  INITIALIZATION OF PROJECTION PARAMETERS (FILE INPUT)  .
! ......................................................................
!
      ENTRY IF11Z0 (INFILE,data)
!
      IERROR = 0
      READ (INFILE,END=060) ZONE,BUFF
      IF (SWITCH.NE.0 .AND. SWITCH.EQ.ZONE) RETURN
  020 A = BUFF(1)
      IF (A .LE. ZERO) A = AZZ
      LON0 = PAKRZ0 (BUFF(5))
      LAT0 = PAKRZ0 (BUFF(6))
      X0 = BUFF(7)
      Y0 = BUFF(8)
      SINPH0 = DSIN (LAT0)
      COSPH0 = DCOS (LAT0)
!
! LIST RESULTS OF PARAMETER INITIALIZATION.
!
      CALL RADDZ0 (LON0,ANGS(1,1))
      CALL RADDZ0 (LAT0,ANGS(1,2))
!     IF (IPPARM .EQ. 0) PRINT 2000, A,ANGS,X0,Y0
      IF (IPPARM .EQ. 0) PRINT 2000, A,angs(1,1),angs(2,1),angs(3,1),  &
          rangs1,angs(1,2),angs(2,2),angs(3,2),rangs2,X0,Y0
 2000 FORMAT (' INITIALIZATION PARAMETERS (LAMBERT AZIMUTHAL EQUAL-AREA'  &
             ,' PROJECTION)'/                                             &
              ' RADIUS OF SPHERE             =',F12.2,' METERS'/          & 
              ' LONGITUDE OF CENTER          = ',A1,2I3,F7.3/             &
              ' LATITUDE  OF CENTER          = ',A1,2I3,F7.3/             & 
              ' FALSE EASTING                =',F12.2,' METERS'/          &
              ' FALSE NORTHING               =',F12.2,' METERS')
      DATA(1) = A
      SWITCH = ZONE
      RETURN
  060 IF (IPEMSG .EQ. 0) PRINT 2010
 2010 FORMAT (' ERROR PJ11Z0'/                    &
              ' MISSING PROJECTION PARAMETERS')
      IERROR = 111
      RETURN
!
! ......................................................................
!      .  INITIALIZATION OF PROJECTION PARAMETERS (ENTRY INPUT)  .
! ......................................................................
!
      ENTRY IS11Z0 (ZZONE,DATA)
      zone = zzone
!
      IERROR = 0
      IF (SWITCH.NE.0 .AND. SWITCH.EQ.ZONE) RETURN
      DO 180 I = 1,8
      BUFF(I) = DATA(I)
  180 CONTINUE
      GO TO 020
!
! ......................................................................
!                      .  FORWARD TRANSFORMATION  .
! ......................................................................
!
      ENTRY PF11Z0 (GEOG,PROJ)
!
      IERROR = 0
      IF (SWITCH .NE. 0) GO TO 120
      IF (IPEMSG .EQ. 0) PRINT 2010
      IERROR = 112
      RETURN
  120 LON = ADJLZ0 (GEOG(1) - LON0)
      SINPHI = DSIN (GEOG(2))
      COSPHI = DCOS (GEOG(2))
      COSLON = DCOS (LON)
      G = SINPH0 * SINPHI + COSPH0 * COSPHI * COSLON
      IF (G .NE. -ONE) GO TO 140
      CON = TWO * A
      IF (IPEMSG .EQ. 0) PRINT 2020, CON
 2020 FORMAT (' POINT PROJECTS INTO A CIRCLE OF RADIUS =',F12.2,' METERS')
      IERROR = 113
      RETURN
  140 KSP = DSQRT (TWO / (ONE + G))
      PROJ(1) = X0 + A * KSP * COSPHI * DSIN (LON)
      PROJ(2) = Y0 + A * KSP * (COSPH0 * SINPHI - SINPH0 * COSPHI * COSLON)
      RETURN
!
! ......................................................................
!                      .  INVERSE TRANSFORMATION  .
! ......................................................................
!
      ENTRY PI11Z0 (PROJ,GEOG)
!
      IERROR = 0
      IF (SWITCH .NE. 0) GO TO 220
      IF (IPEMSG .EQ. 0) PRINT 2010
      IERROR = 114
      RETURN
  220 X = PROJ(1) - X0
      Y = PROJ(2) - Y0
      RH = DSQRT (X * X + Y * Y)
      CON = RH / (TWO * A)
      IF (CON .LE. ONE) GO TO 230
      if (ipemsg .eq. 0) PRINT 2030
 2030 FORMAT (' ERROR PJ11Z0'/      &
              ' INPUT DATA ERROR')
      IERROR = 115
      RETURN
  230 Z = TWO *  DASIN (CON)
      SINZ = DSIN (Z)
      COSZ = DCOS (Z)
      GEOG(1) = LON0
      IF (DABS(RH) .GT. EPSLN) GO TO 240
      GEOG(2) = LAT0
      RETURN
  240 GEOG(2) =  DASIN (COSZ * SINPH0 + Y * SINZ * COSPH0 / RH)
      CON = DABS (LAT0) - HALFPI
      IF (DABS (CON) .GT. EPSLN) GO TO 260
      IF (LAT0 .LT. ZERO) GO TO 250
      GEOG(1) = ADJLZ0 (LON0 + DATAN2 (X , -Y))
      RETURN
  250 GEOG(1) = ADJLZ0 (LON0 - DATAN2 (-X , Y))
      RETURN
  260 CON = COSZ - SINPH0 * DSIN (GEOG(2))
      IF (CON .EQ. ZERO) RETURN
      GEOG(1) = ADJLZ0 (LON0 + DATAN2 ((X*SINZ*COSPH0) , (CON*RH)))
      RETURN
!
      END

!                   PJ12Z0
! **********************************************************************
! ** U.S.G.S. GENERAL MAP PROJECTION PACKAGE ...... DR. A. A. ELASSAL **
! **          MATHEMATICAL ANALYSIS BY JOHN SNYDER                    **
! ** MODULE I                VERSION 1.0.0            NOVEMBER 1,1980 **
! **********************************************************************
!                      *  AZIMUTHAL EQUIDISTANT  *
! **********************************************************************
!
      SUBROUTINE PJ12Z0
!
      IMPLICIT REAL*8 (A-Z)
      integer*4 zzone
      INTEGER*4 IERROR,IPEMSG,IPPARM
      INTEGER*4 SWITCH,I,ZONE,ANGS,INFILE
      COMMON /SPHRZ0/ AZZ
! **** PARAMETERS **** A,LON0,LAT0,X0,Y0,SINPH0,COSPH0 *****************
      COMMON /ERRMZ0/ IERROR
      COMMON /PRINZ0/ IPEMSG,IPPARM
!     COMMON /WORKZ0/ BUFF(15),ANGS(4,2)
      COMMON /WORKZ0/ BUFF(15)
      COMMON /WK12Z0/ ANGS(4,2)
      real*4 rangs1,rangs2
      equivalence (rangs1,angs(4,1))
      equivalence (rangs2,angs(4,2))
      DIMENSION DATA(1),GEOG(1),PROJ(1)
      DATA HALFPI /1.5707963267948966D0/
      DATA EPSLN /1.0D-10/
      DATA ZERO,ONE,TWO /0.0D0,1.0D0,2.0D0/
      DATA SWITCH /0/
!
! ......................................................................
!       .  INITIALIZATION OF PROJECTION PARAMETERS (FILE INPUT)  .
! ......................................................................
!
      ENTRY IF12Z0 (INFILE,data)
!
      IERROR = 0
      READ (INFILE,END=060) ZONE,BUFF
      IF (SWITCH.NE.0 .AND. SWITCH.EQ.ZONE) RETURN
  020 A = BUFF(1)
      IF (A .LE. ZERO) A = AZZ
      LON0 = PAKRZ0 (BUFF(5))
      LAT0 = PAKRZ0 (BUFF(6))
      X0 = BUFF(7)
      Y0 = BUFF(8)
      SINPH0 = DSIN (LAT0)
      COSPH0 = DCOS (LAT0)
!
! LIST RESULTS OF PARAMETER INITIALIZATION.
!
      CALL RADDZ0 (LON0,ANGS(1,1))
      CALL RADDZ0 (LAT0,ANGS(1,2))
!     IF (IPPARM .EQ. 0) PRINT 2000, A,ANGS,X0,Y0
      IF (IPPARM .EQ. 0) PRINT 2000, A,angs(1,1),angs(2,1),angs(3,1),  &
          rangs1,angs(1,2),angs(2,2),angs(3,2),rangs2,X0,Y0
 2000 FORMAT (' INITIALIZATION PARAMETERS (AZIMUTHAL EQUIDISTANT',     &
              ' PROJECTION)'/                                          &  
              ' RADIUS OF SPHERE             =',F12.2,' METERS'/       &
              ' LONGITUDE OF CENTER          = ',A1,2I3,F7.3/          &
              ' LATITUDE  OF CENTER          = ',A1,2I3,F7.3/          &
              ' FALSE EASTING                =',F12.2,' METERS'/       &  
              ' FALSE NORTHING               =',F12.2,' METERS')
      DATA(1) = A
      SWITCH = ZONE
      RETURN
  060 IF (IPEMSG .EQ. 0) PRINT 2010
 2010 FORMAT (' ERROR PJ12Z0'/                    &
              ' MISSING PROJECTION PARAMETERS')
      IERROR = 121
      RETURN
!
! ......................................................................
!      .  INITIALIZATION OF PROJECTION PARAMETERS (ENTRY INPUT)  .
! ......................................................................
!
      ENTRY IS12Z0 (ZZONE,DATA)
      zone = zzone
!
      IERROR = 0
      IF (SWITCH.NE.0 .AND. SWITCH.EQ.ZONE) RETURN
      DO 180 I = 1,8
      BUFF(I) = DATA(I)
  180 CONTINUE
      GO TO 020
!
! ......................................................................
!                      .  FORWARD TRANSFORMATION  .
! ......................................................................
!
      ENTRY PF12Z0 (GEOG,PROJ)
!
      IERROR = 0
      IF (SWITCH .NE. 0) GO TO 120
      IF (IPEMSG .EQ. 0) PRINT 2010
      IERROR = 122
      RETURN
  120 LON = ADJLZ0 (GEOG(1) - LON0)
      SINPHI = DSIN (GEOG(2))
      COSPHI = DCOS (GEOG(2))
      COSLON = DCOS (LON)
      G = SINPH0 * SINPHI + COSPH0 * COSPHI * COSLON
      IF (DABS(DABS(G) - ONE) .GE. EPSLN) GO TO 140
      KSP = ONE
      IF (G .GE. ZERO) GO TO 160
      CON = TWO * HALFPI * A
      IF (IPEMSG .EQ. 0) PRINT 2020, CON
 2020 FORMAT (' POINT PROJECTS INTO CIRCLE OF RADIUS =',F12.2,' METERS')
      IERROR = 123
      RETURN
  140 Z =  DACOS (G)
      KSP = Z / DSIN (Z)
  160 PROJ(1) = X0 + A * KSP * COSPHI * DSIN (LON)
      PROJ(2) = Y0 + A * KSP * (COSPH0 * SINPHI - SINPH0 * COSPHI * COSLON)
      RETURN
!
! ......................................................................
!                      .  INVERSE TRANSFORMATION  .
! ......................................................................
!
      ENTRY PI12Z0 (PROJ,GEOG)
!
      IERROR = 0
      IF (SWITCH .NE. 0) GO TO 220
      IF (IPEMSG .EQ. 0) PRINT 2010
      IERROR = 124
      RETURN
  220 X = PROJ(1) - X0
      Y = PROJ(2) - Y0
      RH = DSQRT (X * X + Y * Y)
      IF (RH .LE. (TWO * HALFPI * A)) GO TO 230
      IF (IPEMSG .EQ. 0) PRINT 2030
 2030 FORMAT (' ERROR PJ12Z0'/       &
              ' INPUT DATA ERROR')
      IERROR = 125
      RETURN
  230 Z = RH / A
      SINZ = DSIN (Z)
      COSZ = DCOS (Z)
      GEOG(1) = LON0
      IF (DABS(RH) .GT. EPSLN) GO TO 240
      GEOG(2) = LAT0
      RETURN
  240 GEOG(2) =  DASIN (COSZ * SINPH0 + Y * SINZ * COSPH0 / RH)
      CON = DABS (LAT0) - HALFPI
      IF (DABS (CON) .GT. EPSLN) GO TO 260
      IF (LAT0 .LT. ZERO) GO TO 250
      GEOG(1) = ADJLZ0 (LON0 + DATAN2 (X , -Y))
      RETURN
  250 GEOG(1) = ADJLZ0 (LON0 - DATAN2 (-X , Y))
      RETURN
  260 CON = COSZ - SINPH0 * DSIN (GEOG(2))
      IF (CON .EQ. ZERO) RETURN
      GEOG(1) = ADJLZ0 (LON0 + DATAN2 ((X*SINZ*COSPH0) , (CON*RH)))
      RETURN
!
      END

!                   PJ13Z0
! **********************************************************************
! ** U.S.G.S. GENERAL MAP PROJECTION PACKAGE ...... DR. A. A. ELASSAL **
! **          MATHEMATICAL ANALYSIS BY JOHN SNYDER                    **
! ** MODULE I                VERSION 1.0.0            NOVEMBER 1,1980 **
! **********************************************************************
!                            *  GNOMONIC  *
! **********************************************************************
!
      SUBROUTINE PJ13Z0
!
      IMPLICIT REAL*8 (A-Z)
      integer*4 zzone
      INTEGER*4 IERROR,IPEMSG,IPPARM
      INTEGER*4 SWITCH,I,ZONE,ANGS,INFILE
      COMMON /SPHRZ0/ AZZ
! **** PARAMETERS **** A,LON0,LAT0,X0,Y0,SINPH0,COSPH0 *****************
      COMMON /ERRMZ0/ IERROR
      COMMON /PRINZ0/ IPEMSG,IPPARM
!     COMMON /WORKZ0/ BUFF(15),ANGS(4,2)
      COMMON /WORKZ0/ BUFF(15)
      COMMON /WK13Z0/ ANGS(4,2)
      real*4 rangs1,rangs2
      equivalence (rangs1,angs(4,1))
      equivalence (rangs2,angs(4,2))
      DIMENSION DATA(1),GEOG(1),PROJ(1)
      DATA HALFPI /1.5707963267948966D0/
      DATA EPSLN /1.0D-10/
      DATA ZERO,ONE /0.0D0,1.0D0/
      DATA SWITCH /0/
!
! ......................................................................
!       .  INITIALIZATION OF PROJECTION PARAMETERS (FILE INPUT)  .
! ......................................................................
!
      ENTRY IF13Z0 (INFILE,data)
!
      IERROR = 0
      READ (INFILE,END=060) ZONE,BUFF
      IF (SWITCH.NE.0 .AND. SWITCH.EQ.ZONE) RETURN
  020 A = BUFF(1)
      IF (A .LE. ZERO) A = AZZ
      LON0 = PAKRZ0 (BUFF(5))
      LAT0 = PAKRZ0 (BUFF(6))
      X0 = BUFF(7)
      Y0 = BUFF(8)
      SINPH0 = DSIN (LAT0)
      COSPH0 = DCOS (LAT0)
!
! LIST RESULTS OF PARAMETER INITIALIZATION.
!
      CALL RADDZ0 (LON0,ANGS(1,1))
      CALL RADDZ0 (LAT0,ANGS(1,2))
!     IF (IPPARM .EQ. 0) PRINT 2000, A,ANGS,X0,Y0
      IF (IPPARM .EQ. 0) PRINT 2000, A,angs(1,1),angs(2,1),angs(3,1),  &
          rangs1,angs(1,2),angs(2,2),angs(3,2),rangs2,X0,Y0
 2000 FORMAT (' INITIALIZATION PARAMETERS (GNOMONIC PROJECTION)'/  &
              ' RADIUS OF SPHERE             =',F12.2,' METERS'/   &
              ' LONGITUDE OF CENTER          = ',A1,2I3,F7.3/      &
              ' LATITUDE  OF CENTER          = ',A1,2I3,F7.3/      &
              ' FALSE EASTING                =',F12.2,' METERS'/   &
              ' FALSE NORTHING               =',F12.2,' METERS')
      DATA(1) = A
      SWITCH = ZONE
      RETURN
  060 IF (IPEMSG .EQ. 0) PRINT 2010
 2010 FORMAT (' ERROR PJ13Z0'/                  &
              ' MISSING PROJECTION PARAMETERS')
      IERROR = 131
      RETURN
!
! ......................................................................
!      .  INITIALIZATION OF PROJECTION PARAMETERS (ENTRY INPUT)  .
! ......................................................................
!
      ENTRY IS13Z0 (ZZONE,DATA)
      zone = zzone
!
      IERROR = 0
      IF (SWITCH.NE.0 .AND. SWITCH.EQ.ZONE) RETURN
      DO 180 I = 1,8
      BUFF(I) = DATA(I)
  180 CONTINUE
      GO TO 020
!
! ......................................................................
!                      .  FORWARD TRANSFORMATION  .
! ......................................................................
!
      ENTRY PF13Z0 (GEOG,PROJ)
!
      IERROR = 0
      IF (SWITCH .NE. 0) GO TO 120
      IF (IPEMSG .EQ. 0) PRINT 2010
      IERROR = 132
      RETURN
  120 LON = ADJLZ0 (GEOG(1) - LON0)
      SINPHI = DSIN (GEOG(2))
      COSPHI = DCOS (GEOG(2))
      COSLON = DCOS (LON)
      G = SINPH0 * SINPHI + COSPH0 * COSPHI * COSLON
      IF (G .GT. ZERO) GO TO 140
      IF (IPEMSG .EQ. 0) PRINT 2020
 2020 FORMAT (' POINT MAPS INTO INFINITY')
      IERROR = 133
      RETURN
  140 KSP = ONE / G
      PROJ(1) = X0 + A * KSP * COSPHI * DSIN (LON)
      PROJ(2) = Y0 + A * KSP * (COSPH0 * SINPHI - SINPH0 * COSPHI * COSLON)
      RETURN
!
! ......................................................................
!                      .  INVERSE TRANSFORMATION  .
! ......................................................................
!
      ENTRY PI13Z0 (PROJ,GEOG)
!
      IERROR = 0
      IF (SWITCH .NE. 0) GO TO 220
      IF (IPEMSG .EQ. 0) PRINT 2010
      IERROR = 134
      RETURN
  220 X = PROJ(1) - X0
      Y = PROJ(2) - Y0
      RH = DSQRT (X * X + Y * Y)
      Z = DATAN (RH / A)
      SINZ = DSIN (Z)
      COSZ = DCOS (Z)
      GEOG(1) = LON0
      IF (DABS(RH) .GT. EPSLN) GO TO 240
      GEOG(2) = LAT0
      RETURN
  240 GEOG(2) =  DASIN (COSZ * SINPH0 + Y * SINZ * COSPH0 / RH)
      CON = DABS (LAT0) - HALFPI
      IF (DABS (CON) .GT. EPSLN) GO TO 260
      IF (LAT0 .LT. ZERO) GO TO 250
      GEOG(1) = ADJLZ0 (LON0 + DATAN2 (X , -Y))
      RETURN
  250 GEOG(1) = ADJLZ0 (LON0 - DATAN2 (-X , Y))
      RETURN
  260 CON = COSZ - SINPH0 * DSIN (GEOG(2))
      IF (CON .EQ. ZERO) RETURN
      GEOG(1) = ADJLZ0 (LON0 + DATAN2 ((X*SINZ*COSPH0) , (CON*RH)))
      RETURN
!
      END

!                   PJ14Z0
! **********************************************************************
! ** U.S.G.S. GENERAL MAP PROJECTION PACKAGE ...... DR. A. A. ELASSAL **
! **          MATHEMATICAL ANALYSIS BY JOHN SNYDER                    **
! ** MODULE I                VERSION 1.0.0            NOVEMBER 1,1980 **
! **********************************************************************
!                          *  ORTHOGRAPHIC  *
! **********************************************************************
!
      SUBROUTINE PJ14Z0
!
      IMPLICIT REAL*8 (A-Z)
      integer*4 zzone
      INTEGER*4 IERROR,IPEMSG,IPPARM
      INTEGER*4 SWITCH,I,ZONE,ANGS,INFILE
      COMMON /SPHRZ0/ AZZ
! **** PARAMETERS **** A,LON0,LAT0,X0,Y0,SINPH0,COSPH0 *****************
      COMMON /ERRMZ0/ IERROR
      COMMON /PRINZ0/ IPEMSG,IPPARM
!     COMMON /WORKZ0/ BUFF(15),ANGS(4,2)
      COMMON /WORKZ0/ BUFF(15)
      COMMON /WK14Z0/ ANGS(4,2)
      real*4 rangs1,rangs2
      equivalence (rangs1,angs(4,1))
      equivalence (rangs2,angs(4,2))
      DIMENSION DATA(1),GEOG(1),PROJ(1)
      DATA HALFPI /1.5707963267948966D0/
      DATA EPSLN /1.0D-10/
      DATA ZERO,ONE /0.0D0,1.0D0/
      DATA SWITCH /0/
!
! ......................................................................
!       .  INITIALIZATION OF PROJECTION PARAMETERS (FILE INPUT)  .
! ......................................................................
!
      ENTRY IF14Z0 (INFILE,data)
!
      IERROR = 0
      READ (INFILE,END=060) ZONE,BUFF
      IF (SWITCH.NE.0 .AND. SWITCH.EQ.ZONE) RETURN
  020 A = BUFF(1)
      IF (A .LE. ZERO) A = AZZ
      LON0 = PAKRZ0 (BUFF(5))
      LAT0 = PAKRZ0 (BUFF(6))
      X0 = BUFF(7)
      Y0 = BUFF(8)
      SINPH0 = DSIN (LAT0)
      COSPH0 = DCOS (LAT0)
!
! LIST RESULTS OF PARAMETER INITIALIZATION.
!
      CALL RADDZ0 (LON0,ANGS(1,1))
      CALL RADDZ0 (LAT0,ANGS(1,2))
!     IF (IPPARM .EQ. 0) PRINT 2000, A,ANGS,X0,Y0
      IF (IPPARM .EQ. 0) PRINT 2000, A,angs(1,1),angs(2,1),angs(3,1),  &
          rangs1,angs(1,2),angs(2,2),angs(3,2),rangs2,X0,Y0
 2000 FORMAT (' INITIALIZATION PARAMETERS (ORTHOGRAPHIC',          &
              ' PROJECTION)'/                                      &
              ' RADIUS OF SPHERE             =',F12.2,' METERS'/   &
              ' LONGITUDE OF CENTER          = ',A1,2I3,F7.3/      & 
              ' LATITUDE  OF CENTER          = ',A1,2I3,F7.3/      &
              ' FALSE EASTING                =',F12.2,' METERS'/   &
              ' FALSE NORTHING               =',F12.2,' METERS')
      DATA(1) = A
      SWITCH = ZONE
      RETURN
  060 IF (IPEMSG .EQ. 0) PRINT 2010
 2010 FORMAT (' ERROR PJ14Z0'/                   &
              ' MISSING PROJECTION PARAMETERS')
      IERROR = 141
      RETURN
!
! ......................................................................
!      .  INITIALIZATION OF PROJECTION PARAMETERS (ENTRY INPUT)  .
! ......................................................................
!
      ENTRY IS14Z0 (ZZONE,DATA)
      zone = zzone
!
      IERROR = 0
      IF (SWITCH.NE.0 .AND. SWITCH.EQ.ZONE) RETURN
      DO 180 I = 1,8
      BUFF(I) = DATA(I)
  180 CONTINUE
      GO TO 020
!
! ......................................................................
!                      .  FORWARD TRANSFORMATION  .
! ......................................................................
!
      ENTRY PF14Z0 (GEOG,PROJ)
!
      IERROR = 0
      IF (SWITCH .NE. 0) GO TO 120
      IF (IPEMSG .EQ. 0) PRINT 2010
      IERROR = 142
      RETURN
  120 LON = ADJLZ0 (GEOG(1) - LON0)
      SINPHI = DSIN (GEOG(2))
      COSPHI = DCOS (GEOG(2))
      COSLON = DCOS (LON)
      G = SINPH0 * SINPHI + COSPH0 * COSPHI * COSLON
      KSP = ONE
      IF (G.GT.ZERO .OR. DABS(G).LE.EPSLN) GO TO 140
      IF (IPEMSG .EQ. 0) PRINT 2020
 2020 FORMAT (' POINT CANNOT BE PROJECTED')
      IERROR = 143
      RETURN
  140 PROJ(1) = X0 + A * KSP * COSPHI * DSIN (LON)
      PROJ(2) = Y0 + A * KSP * (COSPH0 * SINPHI - SINPH0 * COSPHI * COSLON)
      RETURN
!
! ......................................................................
!                      .  INVERSE TRANSFORMATION  .
! ......................................................................
!
      ENTRY PI14Z0 (PROJ,GEOG)
!
      IERROR = 0
      IF (SWITCH .NE. 0) GO TO 220
      IF (IPEMSG .EQ. 0) PRINT 2010
      IERROR = 144
      RETURN
  220 X = PROJ(1) - X0
      Y = PROJ(2) - Y0
      RH = DSQRT (X * X + Y * Y)
      IF (RH .LE. A) GO TO 230
      IF (IPEMSG .EQ. 0) PRINT 2030
 2030 FORMAT (' ERROR PJ14Z0'/     &
              ' INPUT DATA ERROR')
      IERROR = 145
      RETURN
  230 Z =  DASIN (RH / A)
      SINZ = DSIN (Z)
      COSZ = DCOS (Z)
      GEOG(1) = LON0
      IF (DABS(RH) .GT. EPSLN) GO TO 240
      GEOG(2) = LAT0
      RETURN
  240 GEOG(2) =  DASIN (COSZ * SINPH0 + Y * SINZ * COSPH0 / RH)
      CON = DABS (LAT0) - HALFPI
      IF (DABS (CON) .GT. EPSLN) GO TO 260
      IF (LAT0 .LT. ZERO) GO TO 250
      GEOG(1) = ADJLZ0 (LON0 + DATAN2 (X , -Y))
      RETURN
  250 GEOG(1) = ADJLZ0 (LON0 - DATAN2 (-X , Y))
      RETURN
  260 CON = COSZ - SINPH0 * DSIN (GEOG(2))
      IF (CON .EQ. ZERO) RETURN
      GEOG(1) = ADJLZ0 (LON0 + DATAN2 ((X*SINZ*COSPH0) , (CON*RH)))
      RETURN
!
      END

!                   PJ15Z0
! **********************************************************************
! ** U.S.G.S. GENERAL MAP PROJECTION PACKAGE ...... DR. A. A. ELASSAL **
! **          MATHEMATICAL ANALYSIS BY JOHN SNYDER                    **
! ** MODULE I                VERSION 1.0.0            NOVEMBER 1,1980 **
! **********************************************************************
!              *  GENERAL VERTICAL NEAR-SIDE PERSPECTIVE  *
! **********************************************************************
!
      SUBROUTINE PJ15Z0
!
      IMPLICIT REAL*8 (A-Z)
      integer*4 zzone
      INTEGER*4 IERROR,IPEMSG,IPPARM
      INTEGER*4 SWITCH,I,ZONE,ANGS,INFILE
      COMMON /SPHRZ0/ AZZ
! **** PARAMETERS **** A,P,LON0,LAT0,X0,Y0,SINPH0,COSPH0 ***************
      COMMON /ERRMZ0/ IERROR
      COMMON /PRINZ0/ IPEMSG,IPPARM
!     COMMON /WORKZ0/ BUFF(15),ANGS(4,2)
      COMMON /WORKZ0/ BUFF(15)
      COMMON /WK15Z0/ ANGS(4,2)
      real*4 rangs1,rangs2
      equivalence (rangs1,angs(4,1))
      equivalence (rangs2,angs(4,2))
      DIMENSION DATA(1),GEOG(1),PROJ(1)
      DATA HALFPI /1.5707963267948966D0/
      DATA EPSLN /1.0D-10/
      DATA ZERO,ONE /0.0D0,1.0D0/
      DATA SWITCH /0/
!
! ......................................................................
!       .  INITIALIZATION OF PROJECTION PARAMETERS (FILE INPUT)  .
! ......................................................................
!
      ENTRY IF15Z0 (INFILE,data)
!
      IERROR = 0
      READ (INFILE,END=060) ZONE,BUFF
      IF (SWITCH.NE.0 .AND. SWITCH.EQ.ZONE) RETURN
  020 A = BUFF(1)
      IF (A .LE. ZERO) A = AZZ
      P = ONE + BUFF(3) / A
      LON0 = PAKRZ0 (BUFF(5))
      LAT0 = PAKRZ0 (BUFF(6))
      X0 = BUFF(7)
      Y0 = BUFF(8)
      SINPH0 = DSIN (LAT0)
      COSPH0 = DCOS (LAT0)
!
! LIST RESULTS OF PARAMETER INITIALIZATION.
!
      CALL RADDZ0 (LON0,ANGS(1,1))
      CALL RADDZ0 (LAT0,ANGS(1,2))
!     IF (IPPARM .EQ. 0) PRINT 2000, A,BUFF(3),ANGS,X0,Y0
      IF (IPPARM .EQ. 0) PRINT 2000, A,buff(3),angs(1,1),angs(2,1),angs(3,1),  &
          rangs1,angs(1,2),angs(2,2),angs(3,2),rangs2,X0,Y0
 2000 FORMAT (' INITIALIZATION PARAMETERS (GENERAL VERTICAL NEAR-SIDE',  &
              ' PERSPECTIVE PROJECTION)'/                                &
              ' RADIUS OF SPHERE             =',F12.2,' METERS'/         &
              ' HEIGHT OF PERSPECTIVE POINT'/                            &
              ' ABOVE SPHERE                 =',F12.2,' METERS'/         &
              ' LONGITUDE OF CENTER          = ',A1,2I3,F7.3/            &
              ' LATITUDE  OF CENTER          = ',A1,2I3,F7.3/            &
              ' FALSE EASTING                =',F12.2,' METERS'/         &
              ' FALSE NORTHING               =',F12.2,' METERS')
      DATA(1) = A
      SWITCH = ZONE
      RETURN
  060 IF (IPEMSG .EQ. 0) PRINT 2010
 2010 FORMAT (' ERROR PJ15Z0'/                   &
              ' MISSING PROJECTION PARAMETERS')
      IERROR = 151
      RETURN
!
! ......................................................................
!      .  INITIALIZATION OF PROJECTION PARAMETERS (ENTRY INPUT)  .
! ......................................................................
!
      ENTRY IS15Z0 (ZZONE,DATA)
      zone = zzone
!
      IERROR = 0
      IF (SWITCH.NE.0 .AND. SWITCH.EQ.ZONE) RETURN
      DO 180 I = 1,8
      BUFF(I) = DATA(I)
  180 CONTINUE
      GO TO 020
!
! ......................................................................
!                      .  FORWARD TRANSFORMATION  .
! ......................................................................
!
      ENTRY PF15Z0 (GEOG,PROJ)
!
      IERROR = 0
      IF (SWITCH .NE. 0) GO TO 120
      IF (IPEMSG .EQ. 0) PRINT 2010
      IERROR = 152
      RETURN
  120 LON = ADJLZ0 (GEOG(1) - LON0)
      SINPHI = DSIN (GEOG(2))
      COSPHI = DCOS (GEOG(2))
      COSLON = DCOS (LON)
      G = SINPH0 * SINPHI + COSPH0 * COSPHI * COSLON
      IF (G .GE. (ONE / P)) GO TO 140
      IF (IPEMSG .EQ. 0) PRINT 2020
 2020 FORMAT (' POINT CANNOT BE PROJECTED')
      IERROR = 153
      RETURN
  140 KSP = (P - ONE) / (P - G)
      PROJ(1) = X0 + A * KSP * COSPHI * DSIN (LON)
      PROJ(2) = Y0 + A * KSP * (COSPH0 * SINPHI - SINPH0 * COSPHI * COSLON)
      RETURN
!
! ......................................................................
!                      .  INVERSE TRANSFORMATION  .
! ......................................................................
!
      ENTRY PI15Z0 (PROJ,GEOG)
!
      IERROR = 0
      IF (SWITCH .NE. 0) GO TO 220
      IF (IPEMSG .EQ. 0) PRINT 2010
      IERROR = 154
      RETURN
  220 X = PROJ(1) - X0
      Y = PROJ(2) - Y0
      RH = DSQRT (X * X + Y * Y)
      R = RH / A
      CON = P - ONE
      COM = P + ONE
      IF (R .LE. DSQRT (CON / COM)) GO TO 230
      IF (IPEMSG .EQ. 0) PRINT 2030
 2030 FORMAT (' ERROR PJ15Z0'/       &
              ' INPUT DATA ERROR')
      IERROR = 155
      RETURN
  230 SINZ = (P - DSQRT (ONE - R * R * COM / CON)) / (CON / R + R / CON)
      Z =  DASIN (SINZ)
      SINZ = DSIN (Z)
      COSZ = DCOS (Z)
      GEOG(1) = LON0
      IF (DABS(RH) .GT. EPSLN) GO TO 240
      GEOG(2) = LAT0
      RETURN
  240 GEOG(2) =  DASIN (COSZ * SINPH0 + Y * SINZ * COSPH0 / RH)
      CON = DABS (LAT0) - HALFPI
      IF (DABS (CON) .GT. EPSLN) GO TO 260
      IF (LAT0 .LT. ZERO) GO TO 250
      GEOG(1) = ADJLZ0 (LON0 + DATAN2 (X , -Y))
      RETURN
  250 GEOG(1) = ADJLZ0 (LON0 - DATAN2 (-X , Y))
      RETURN
  260 CON = COSZ - SINPH0 * DSIN (GEOG(2))
      IF (CON .EQ. ZERO) RETURN
      GEOG(1) = ADJLZ0 (LON0 + DATAN2 ((X*SINZ*COSPH0) , (CON*RH)))
      RETURN
!
      END

!                   PJ16Z0
! **********************************************************************
! ** U.S.G.S. GENERAL MAP PROJECTION PACKAGE ...... DR. A. A. ELASSAL **
! **          MATHEMATICAL ANALYSIS BY JOHN SNYDER                    **
! ** MODULE I                VERSION 1.0.0            NOVEMBER 1,1980 **
! **********************************************************************
!                           *  SINUSOIDAL  *
! **********************************************************************
!
      SUBROUTINE PJ16Z0
!
      IMPLICIT REAL*8 (A-Z)
      integer*4 zzone
      INTEGER*4 IERROR,IPEMSG,IPPARM
      INTEGER*4 SWITCH,I,ZONE,ANGS,INFILE
      COMMON /SPHRZ0/ AZZ
! **** PARAMETERS **** A,LON0,X0,Y0 ************************************
      COMMON /ERRMZ0/ IERROR
      COMMON /PRINZ0/ IPEMSG,IPPARM
!     COMMON /WORKZ0/ BUFF(15),ANGS(4)
      COMMON /WORKZ0/ BUFF(15)
      COMMON /WK16Z0/ ANGS(4)
      real*4 rangs
      equivalence (rangs,angs(4))
      DIMENSION DATA(1),GEOG(1),PROJ(1)
      DATA HALFPI /1.5707963267948966D0/
      DATA EPSLN /1.0D-10/
      DATA ZERO /0.0D0/
      DATA SWITCH /0/
!
! ......................................................................
!       .  INITIALIZATION OF PROJECTION PARAMETERS (FILE INPUT)  .
! ......................................................................
!
      ENTRY IF16Z0 (INFILE,data)
!
      IERROR = 0
      READ (INFILE,END=060) ZONE,BUFF
      IF (SWITCH.NE.0 .AND. SWITCH.EQ.ZONE) RETURN
  020 A = BUFF(1)
      IF (A .LE. ZERO) A = AZZ
      LON0 = PAKRZ0 (BUFF(5))
      X0 = BUFF(7)
      Y0 = BUFF(8)
!
! LIST RESULTS OF PARAMETER INITIALIZATION.
!
      CALL RADDZ0 (LON0,ANGS)
!     IF (IPPARM .EQ. 0) PRINT 2000, A,ANGS,X0,Y0
      IF (IPPARM .EQ. 0) PRINT 2000, A,angs(1),angs(2),angs(3),rangs,X0,Y0
 2000 FORMAT (' INITIALIZATION PARAMETERS (SINUSOIDAL',             &
              ' PROJECTION)'/                                       & 
              ' RADIUS OF SPHERE             =',F12.2,' METERS'/    & 
              ' LONGITUDE OF C. MERIDIAN     = ',A1,2I3,F7.3/       &
              ' FALSE EASTING                =',F12.2,' METERS'/    &
              ' FALSE NORTHING               =',F12.2,' METERS')
      DATA(1) = A
      SWITCH = ZONE
      RETURN
  060 IF (IPEMSG .EQ. 0) PRINT 2010
 2010 FORMAT (' ERROR PJ16Z0'/                  &
              ' MISSING PROJECTION PARAMETERS')
      IERROR = 161
      RETURN
!
! ......................................................................
!      .  INITIALIZATION OF PROJECTION PARAMETERS (ENTRY INPUT)  .
! ......................................................................
!
      ENTRY IS16Z0 (ZZONE,DATA)
      zone = zzone
!
      IERROR = 0
      IF (SWITCH.NE.0 .AND. SWITCH.EQ.ZONE) RETURN
      DO 080 I = 1,8
      BUFF(I) = DATA(I)
  080 CONTINUE
      GO TO 020
!
! ......................................................................
!                      .  FORWARD TRANSFORMATION  .
! ......................................................................
!
      ENTRY PF16Z0 (GEOG,PROJ)
!
      IERROR = 0
      IF (SWITCH .NE. 0) GO TO 120
      IF (IPEMSG .EQ. 0) PRINT 2010
      IERROR = 162
      RETURN
  120 LON = ADJLZ0 (GEOG(1) - LON0)
      PROJ(1) = X0 + A * LON * DCOS (GEOG(2))
      PROJ(2) = Y0 + A * GEOG(2)
      RETURN
!
! ......................................................................
!                      .  INVERSE TRANSFORMATION  .
! ......................................................................
!
      ENTRY PI16Z0 (PROJ,GEOG)
!
      IERROR = 0
      IF (SWITCH .NE. 0) GO TO 220
      IF (IPEMSG .EQ. 0) PRINT 2010
      IERROR = 163
      RETURN
  220 X = PROJ(1) - X0
      Y = PROJ(2) - Y0
      GEOG(2) = Y / A
      IF (DABS(GEOG(2)) .LE. HALFPI) GO TO 230
      IF (IPEMSG .EQ. 0) PRINT 2020
 2020 FORMAT (' ERROR PJ16Z0'/       &
              ' INPUT DATA ERROR')
      IERROR = 164
      RETURN
  230 CON = DABS (GEOG(2)) - HALFPI
      IF (DABS (CON) .GT. EPSLN) GO TO 240
      GEOG(1) = LON0
      RETURN
  240 GEOG(1) = ADJLZ0 (LON0 + X / (A * DCOS (GEOG(2))))
      RETURN
!
      END

!                   PJ17Z0
! **********************************************************************
! ** U.S.G.S. GENERAL MAP PROJECTION PACKAGE ...... DR. A. A. ELASSAL **
! **          MATHEMATICAL ANALYSIS BY JOHN SNYDER                    **
! ** MODULE I                VERSION 1.0.1            MAY 1 ,1981 ******
! **********************************************************************
!                  *  EQUIRECTANGULAR   *
! **********************************************************************
!
      SUBROUTINE PJ17Z0
!
      IMPLICIT REAL*8 (A-Z)
      integer*4 zzone
      INTEGER*4 IERROR,IPEMSG,IPPARM
      INTEGER*4 SWITCH,I,ZONE,ANGS,INFILE
      COMMON /SPHRZ0/ AZZ
! **** PARAMETERS **** A,LON0,X0,Y0,LAT1 *******************************
      COMMON /ERRMZ0/ IERROR
      COMMON /PRINZ0/ IPEMSG,IPPARM
!     COMMON /WORKZ0/ BUFF(15),ANGS(4,2)
      COMMON /WORKZ0/ BUFF(15)
      COMMON /WK17Z0/ ANGS(4,2)
      real*4 rangs1,rangs2
      equivalence (rangs1,angs(4,1))
      equivalence (rangs2,angs(4,2))
      DIMENSION DATA(1),GEOG(1),PROJ(1)
      DATA HALFPI /1.5707963267948966D0/
      DATA ZERO /0.0D0/
      DATA SWITCH /0/
!
! ......................................................................
!       .  INITIALIZATION OF PROJECTION PARAMETERS (FILE INPUT)  .
! ......................................................................
!
      ENTRY IF17Z0 (INFILE,data)
!
      IERROR = 0
      READ (INFILE,END=060) ZONE,BUFF
      IF (SWITCH.NE.0 .AND. SWITCH.EQ.ZONE) RETURN
  020 A = BUFF(1)
      IF (A .LE. ZERO) A = AZZ
      LAT1 = PAKRZ0 (BUFF(6))
      LON0 = PAKRZ0 (BUFF(5))
      X0 = BUFF(7)
      Y0 = BUFF(8)
!
! LIST RESULTS OF PARAMETER INITIALIZATION.
!
      CALL RADDZ0 (LAT1,ANGS(1,1))
      CALL RADDZ0 (LON0,ANGS(1,2))
!     IF (IPPARM .EQ. 0) PRINT 2000, A,ANGS,X0,Y0
      IF (IPPARM .EQ. 0) PRINT 2000, A,angs(1,1),angs(2,1),angs(3,1),  &
          rangs1,angs(1,2),angs(2,2),angs(3,2),rangs2,X0,Y0
 2000 FORMAT (' INITIALIZATION PARAMETERS (EQUIRECTANGULAR PROJECTION)'/  &
              ' RADIUS OF SPHERE             =',F12.2,' METERS'/          &
              ' LATITUDE OF TRUE SCALE       = ',A1,2I2,F7.3/             & 
              ' LONGITUDE OF C. MERIDIAN     = ',A1,2I3,F7.3/             &
              ' FALSE EASTING                =',F12.2,' METERS'/          &
              ' FALSE NORTHING               =',F12.2,' METERS')
      DATA(1) = A
      SWITCH = ZONE
      RETURN
  060 IF (IPEMSG .EQ. 0) PRINT 2010
 2010 FORMAT (' ERROR PJ17Z0'/                   &
              ' MISSING PROJECTION PARAMETERS')
      IERROR = 171
      RETURN
!
! ......................................................................
!      .  INITIALIZATION OF PROJECTION PARAMETERS (ENTRY INPUT)  .
! ......................................................................
!
      ENTRY IS17Z0 (ZZONE,DATA)
      zone = zzone
!
      IERROR = 0
      IF (SWITCH.NE.0 .AND. SWITCH.EQ.ZONE) RETURN
      DO 080 I = 1,8
      BUFF(I) = DATA(I)
  080 CONTINUE
      GO TO 020
!
! ......................................................................
!                      .  FORWARD TRANSFORMATION  .
! ......................................................................
!
      ENTRY PF17Z0 (GEOG,PROJ)
!
      IERROR = 0
      IF (SWITCH .NE. 0) GO TO 120
      IF (IPEMSG .EQ. 0) PRINT 2010
      IERROR = 172
      RETURN
  120 LON = ADJLZ0 (GEOG(1) - LON0)
      PROJ(1) = X0 + A * LON * DCOS(LAT1)
      PROJ(2) = Y0 + A * GEOG(2)
      RETURN
!
! ......................................................................
!                      .  INVERSE TRANSFORMATION  .
! ......................................................................
!
      ENTRY PI17Z0 (PROJ,GEOG)
!
      IERROR = 0
      IF (SWITCH .NE. 0) GO TO 220
      IF (IPEMSG .EQ. 0) PRINT 2010
      IERROR = 173
      RETURN
  220 X = PROJ(1) - X0
      Y = PROJ(2) - Y0
      GEOG(2) = Y / A
      IF (DABS(GEOG(2)) .LE. HALFPI) GO TO 240
      IF (IPEMSG .EQ. 0) PRINT 2020
 2020 FORMAT (' ERROR PJ17Z0'/       &
              ' INPUT DATA ERROR')
      IERROR = 174
      RETURN
  240 GEOG(1) = ADJLZ0 (LON0 + X / (A * DCOS(LAT1) ))
      RETURN
!
      END

!                   PJ18Z0
! **********************************************************************
! ** U.S.G.S. GENERAL MAP PROJECTION PACKAGE ...... DR. A. A. ELASSAL **
! **          MATHEMATICAL ANALYSIS BY JOHN SNYDER                    **
! ** MODULE I                VERSION 1.0.0            NOVEMBER 1,1980 **
! **********************************************************************
!                       *  MILLER CYLINDRICAL  *
! **********************************************************************
!
      SUBROUTINE PJ18Z0
!
      IMPLICIT REAL*8 (A-Z)
      integer*4 zzone
      INTEGER*4 IERROR,IPEMSG,IPPARM
      INTEGER*4 SWITCH,I,ZONE,ANGS,INFILE
      COMMON /SPHRZ0/ AZZ
! **** PARAMETERS **** A,LON0,X0,Y0 ************************************
      COMMON /ERRMZ0/ IERROR
      COMMON /PRINZ0/ IPEMSG,IPPARM
!     COMMON /WORKZ0/ BUFF(15),ANGS(4)
      COMMON /WORKZ0/ BUFF(15)
      COMMON /WK18Z0/ ANGS(4)
      real*4 rangs
      equivalence (rangs,angs(4))
      DIMENSION DATA(1),GEOG(1),PROJ(1)
      DATA FORTPI /0.78539816339744833D0/
      DATA ZERO,ONEQ,TWOH /0.0D0,1.25D0,2.5D0/
      DATA SWITCH /0/
!
! ......................................................................
!       .  INITIALIZATION OF PROJECTION PARAMETERS (FILE INPUT)  .
! ......................................................................
!
      ENTRY IF18Z0 (INFILE,data)
!
      IERROR = 0
      READ (INFILE,END=060) ZONE,BUFF
      IF (SWITCH.NE.0 .AND. SWITCH.EQ.ZONE) RETURN
  020 A = BUFF(1)
      IF (A .LE. ZERO) A = AZZ
      LON0 = PAKRZ0 (BUFF(5))
      X0 = BUFF(7)
      Y0 = BUFF(8)
!
! LIST RESULTS OF PARAMETER INITIALIZATION.
!
      CALL RADDZ0 (LON0,ANGS)
!     IF (IPPARM .EQ. 0) PRINT 2000, A,ANGS,X0,Y0
      IF (IPPARM .EQ. 0) PRINT 2000, A,angs(1),angs(2),angs(3),rangs,X0,Y0
 2000 FORMAT (' INITIALIZATION PARAMETERS (MILLER CYLINDRICAL',  &
              ' PROJECTION)'/                                    &
              ' RADIUS OF SPHERE             =',F12.2,' METERS'/ &
              ' LONGITUDE OF C. MERIDIAN     = ',A1,2I3,F7.3/    &
              ' FALSE EASTING                =',F12.2,' METERS'/ &
              ' FALSE NORTHING               =',F12.2,' METERS')
      DATA(1) = A
      SWITCH = ZONE
      RETURN
  060 IF (IPEMSG .EQ. 0) PRINT 2010
 2010 FORMAT (' ERROR PJ18Z0'/                    &
              ' MISSING PROJECTION PARAMETERS')
      IERROR = 181
      RETURN
!
! ......................................................................
!      .  INITIALIZATION OF PROJECTION PARAMETERS (ENTRY INPUT)  .
! ......................................................................
!
      ENTRY IS18Z0 (ZZONE,DATA)
      zone = zzone
!
      IERROR = 0
      IF (SWITCH.NE.0 .AND. SWITCH.EQ.ZONE) RETURN
      DO 080 I = 1,8
      BUFF(I) = DATA(I)
  080 CONTINUE
      GO TO 020
!
! ......................................................................
!                      .  FORWARD TRANSFORMATION  .
! ......................................................................
!
      ENTRY PF18Z0 (GEOG,PROJ)
!
      IERROR = 0
      IF (SWITCH .NE. 0) GO TO 120
      IF (IPEMSG .EQ. 0) PRINT 2010
      IERROR = 182
      RETURN
  120 LON = ADJLZ0 (GEOG(1) - LON0)
      PROJ(1) = X0 + A * LON
      PROJ(2) = Y0 + A * DLOG (DTAN (FORTPI + GEOG(2) / TWOH)) * ONEQ
      RETURN
!
! ......................................................................
!                      .  INVERSE TRANSFORMATION  .
! ......................................................................
!
      ENTRY PI18Z0 (PROJ,GEOG)
!
      IERROR = 0
      IF (SWITCH .NE. 0) GO TO 220
      IF (IPEMSG .EQ. 0) PRINT 2010
      IERROR = 183
      RETURN
  220 X = PROJ(1) - X0
      Y = PROJ(2) - Y0
      GEOG(1) = ADJLZ0 (LON0 + X / A)
      GEOG(2) = TWOH * DATAN (DEXP (Y / A / ONEQ)) - FORTPI * TWOH
      RETURN
!
      END

!                   PJ19Z0
! **********************************************************************
! ** U.S.G.S. GENERAL MAP PROJECTION PACKAGE ...... DR. A. A. ELASSAL **
! **          MATHEMATICAL ANALYSIS BY JOHN SNYDER                    **
! ** MODULE I                VERSION 1.0.1            MAY 1, 1981 ******
! **********************************************************************
!                        *  VAN DER GRINTEN I  *
! **********************************************************************
!
      SUBROUTINE PJ19Z0
!
      IMPLICIT REAL*8 (A-Z)
      integer*4 zzone
      INTEGER*4 IERROR,IPEMSG,IPPARM
      INTEGER*4 SWITCH,I,ZONE,ANGS,NIT,INFILE
      COMMON /SPHRZ0/ AZZ
! **** PARAMETERS **** A,LON0,X0,Y0 ************************************
      COMMON /ERRMZ0/ IERROR
      COMMON /PRINZ0/ IPEMSG,IPPARM
!     COMMON /WORKZ0/ BUFF(15),ANGS(4)
      COMMON /WORKZ0/ BUFF(15)
      COMMON /WK19Z0/ ANGS(4)
      real*4 rangs
      equivalence (rangs,angs(4))
      DIMENSION DATA(1),GEOG(1),PROJ(1)
      DATA PI /3.14159265358979323846D0/
      DATA HALFPI /1.5707963267948966D0/
      DATA EPSLN,TOL,NIT /1.0D-10,0.7D0,35/
      DATA ZERO,HALF,ONE,TWO,FOUR /0.0D0,0.5D0,1.0D0,2.0D0,4.0D0/
      DATA SWITCH /0/
!
! ......................................................................
!       .  INITIALIZATION OF PROJECTION PARAMETERS (FILE INPUT)  .
! ......................................................................
!
      ENTRY IF19Z0 (INFILE,data)
!
      IERROR = 0
      READ (INFILE,END=060) ZONE,BUFF
      IF (SWITCH.NE.0 .AND. SWITCH.EQ.ZONE) RETURN
  020 A = BUFF(1)
      IF (A .LE. ZERO) A = AZZ
      LON0 = PAKRZ0 (BUFF(5))
      X0 = BUFF(7)
      Y0 = BUFF(8)
!
! LIST RESULTS OF PARAMETER INITIALIZATION.
!
      CALL RADDZ0 (LON0,ANGS)
!     IF (IPPARM .EQ. 0) PRINT 2000, A,ANGS,X0,Y0
      IF (IPPARM .EQ. 0) PRINT 2000, A,angs(1),angs(2),angs(3),rangs,X0,Y0
 2000 FORMAT (' INITIALIZATION PARAMETERS (VAN DER GRINTEN I',     &
              ' PROJECTION)'/                                      &
              ' RADIUS OF SPHERE             =',F12.2,' METERS'/   &
              ' LONGITUDE OF C. MERIDIAN     = ',A1,2I3,F7.3/      & 
              ' FALSE EASTING                =',F12.2,' METERS'/   &
              ' FALSE NORTHING               =',F12.2,' METERS')
      DATA(1) = A
      SWITCH = ZONE
      RETURN
  060 IF (IPEMSG .EQ. 0) PRINT 2010
 2010 FORMAT (' ERROR PJ19Z0'/                    &
              ' MISSING PROJECTION PARAMETERS')
      IERROR = 191
      RETURN
!
! ......................................................................
!      .  INITIALIZATION OF PROJECTION PARAMETERS (ENTRY INPUT)  .
! ......................................................................
!
      ENTRY IS19Z0 (ZZONE,DATA)
      zone = zzone
!
      IERROR = 0
      IF (SWITCH.NE.0 .AND. SWITCH.EQ.ZONE) RETURN
      DO 080 I = 1,8
      BUFF(I) = DATA(I)
  080 CONTINUE
      GO TO 020
!
! ......................................................................
!                      .  FORWARD TRANSFORMATION  .
! ......................................................................
!
      ENTRY PF19Z0 (GEOG,PROJ)
!
      IERROR = 0
      IF (SWITCH .NE. 0) GO TO 120
      IF (IPEMSG .EQ. 0) PRINT 2010
      IERROR = 192
      RETURN
  120 LON = ADJLZ0 (GEOG(1) - LON0)
      LAT = GEOG(2)
      IF (DABS(LAT) .GT. EPSLN) GO TO 140
      PROJ(1) = X0 + A * LON
      PROJ(2) = Y0
      RETURN
  140 THETA =  DASIN (DABS (LAT /HALFPI))
      IF (DABS(LON) .GT. EPSLN) GO TO 160
      PROJ(1) = X0
      PROJ(2) = Y0 + PI * A * DSIGN( DTAN (HALF * THETA), LAT)
      RETURN
  160 AL = HALF * DABS (PI / LON - LON / PI)
      ASQ = AL * AL
      SINTHT = DSIN (THETA)
      COSTHT = DCOS (THETA)
      G = COSTHT / (SINTHT + COSTHT - ONE)
      GSQ = G * G
      M = G * (TWO / SINTHT - ONE)
      MSQ = M * M
      CON = PI * A * (AL * (G - MSQ) + DSQRT (ASQ * (G - MSQ)**2 - (MSQ + ASQ) * (GSQ - MSQ))) / (MSQ + ASQ)
      CON = DSIGN (CON , LON)
      PROJ(1) = X0 + CON
      CON = DABS (CON / (PI * A))
      PROJ(2) = Y0 + DSIGN (PI * A * DSQRT (ONE - CON * CON - TWO * AL * CON) , LAT)
      RETURN
!
! ......................................................................
!                      .  INVERSE TRANSFORMATION  .
! ......................................................................
!
      ENTRY PI19Z0 (PROJ,GEOG)
!
      IERROR = 0
      IF (SWITCH .NE. 0) GO TO 220
      IF (IPEMSG .EQ. 0) PRINT 2010
      IERROR = 193
      RETURN
  220 X = PROJ(1) - X0
      Y = PROJ(2) - Y0
      CON = DABS (Y / (PI * A))
      THETA = TWO * DATAN (CON)
      IF (DABS(X) .GT. EPSLN) GO TO 240
      GEOG(1) = LON0
      GEOG(2) = HALFPI * DSIGN( DSIN (THETA), Y)
      RETURN
  240 IF (DABS(Y) .GT. EPSLN) GO TO 260
      GEOG(1) = ADJLZ0 (LON0 + X / A)
      GEOG(2) = ZERO
      RETURN
  260 IF (DSQRT(X*X+Y*Y) .LE. PI*A) GO TO 270
      IF (IPEMSG .EQ. 0) PRINT 2020
 2020 FORMAT (' ERROR PI19Z0'/       &
              ' INPUT DATA ERROR')
      IERROR = 194
      RETURN
  270 CNN = CON * CON
      COM = DABS (X / (PI * A))
      CMM = COM * COM
      AL = (ONE - CMM - CNN) / (TWO * COM)
      GEOG(1) = ADJLZ0 (LON0 + DSIGN (PI*(-AL + DSQRT (AL*AL+ONE)) , X))
      PHI = THETA
      IF (CON .GT. TOL) GO TO 320
!
! LOW LATITUDE CASE
!
      DO 280 I = 1,NIT
      THETA =  DASIN (PHI / HALFPI)
      SINTHT = DSIN (THETA)
      COSTHT = DCOS (THETA)
      G = COSTHT / (SINTHT + COSTHT - ONE)
      D = CON / SINTHT - ONE / (ONE + COSTHT)
      H = TWO - SINTHT
      J = DTAN (HALF * THETA)
      DPHI = (CMM + CNN - TWO * D * G * H - J * J) * PI * COSTHT /          &
             (FOUR * (G * H * (CON * COSTHT / (ONE - COSTHT) + J) /         &
             (ONE + COSTHT) + D * G * ((ONE + TWO * COSTHT * COSTHT) /      &
             COSTHT + H * (COSTHT - SINTHT) / (SINTHT + COSTHT - ONE)) -    &
             J * (J * J + ONE)))
      PHI = PHI - DPHI
      IF (DABS(DPHI) .LT. EPSLN) GO TO 400
  280 CONTINUE
  300 IF (IPEMSG .EQ. 0) PRINT 2030, NIT
 2030 FORMAT (' ERROR PI19Z0'/                                         &
              ' LATITUDE FAILED TO CONVERGE AFTER',I3,' ITERATIONS')
      IERROR = 195
      RETURN
!
! HIGH LATITUDE CASE.
!
  320 LON = ADJLZ0 (GEOG(1) - LON0)
      DO 380 I = 1,NIT
      IF (DABS(PHI) .GT. EPSLN) GO TO 330
      Y1 = ZERO
      GO TO 360
  330 THETA =  DASIN (DABS (PHI /HALFPI))
      IF (DABS(LON) .GT. EPSLN) GO TO 340
      Y1 = PI * A * DTAN (HALF * THETA)
      GO TO 360
  340 AL = HALF * DABS (PI / LON - LON / PI)
      ASQ = AL * AL
      SINTHT = DSIN (THETA)
      COSTHT = DCOS (THETA)
      G = COSTHT / (SINTHT + COSTHT - ONE)
      GSQ = G * G
      M = G * (TWO / SINTHT - ONE)
      MSQ = M * M
      CON = DABS ((AL * (G - MSQ) + DSQRT (ASQ * (G - MSQ)**2 - (MSQ + ASQ) * (GSQ - MSQ))) / (MSQ + ASQ))
      Y1 = DSIGN (PI * A * DSQRT (ONE - CON * CON - TWO * AL * CON) , PHI)
  360 DPHI = ((DABS(Y) - Y1) / (PI * A - Y1)) * (HALFPI - PHI)
      PHI = PHI + DPHI
      IF (DABS(DPHI) .LT. EPSLN) GO TO 400
  380 CONTINUE
      GO TO 300
  400 GEOG(2) = DSIGN (PHI , Y)
      RETURN
!
      END

!                   PJ20Z0
! **********************************************************************
! ** U.S.G.S. GENERAL MAP PROJECTION PACKAGE ...... DR. A. A. ELASSAL **
! **          MATHEMATICAL ANALYSIS BY JOHN SNYDER                    **
! ** MODULE I                VERSION 1.0.1            MAY 1, 1981 ******
! **********************************************************************
!                    *  OBLIQUE MERCATOR (HOTINE)  *
! **********************************************************************
!
      SUBROUTINE PJ20Z0
!
      IMPLICIT REAL*8 (A-Z)
      integer*4 zzone
      INTEGER*4 IERROR,IPEMSG,IPPARM
      INTEGER*4 SWITCH,I,ZONE,ANGS1,ANGS2,MODE,INFILE
      COMMON /ELLPZ0/ AZ,EZ,ESZ,E0Z,E1Z,E2Z,E3Z
! **** PARAMETERS **** A,E,ES,KS0,ALPHA,LONC,LON1,LAT1,LON2,LAT2,LAT0 **
! ********************** X0,Y0,GAMMA,LON0,AL,BL,EL *********************
      COMMON /ERRMZ0/ IERROR
      COMMON /PRINZ0/ IPEMSG,IPPARM
!     COMMON /WORKZ0/ BUFF(15),ANGS1(4,5),ANGS2(4,3)
      COMMON /WORKZ0/ BUFF(15)
      COMMON /WK20Z0/ ANGS1(4,5),ANGS2(4,3)
      real*4 rangs1,rangs2,rangs3,rangs11,rangs12,rangs13,rangs14,rangs15
      equivalence (rangs1,angs2(4,1))
      equivalence (rangs2,angs2(4,2))
      equivalence (rangs3,angs2(4,3))
      equivalence (rangs11,angs1(4,1))
      equivalence (rangs12,angs1(4,2))
      equivalence (rangs13,angs1(4,3))
      equivalence (rangs14,angs1(4,4))
      equivalence (rangs15,angs1(4,5))
      DIMENSION DATA(1),GEOG(1),PROJ(1)
      DATA PI /3.14159265358979323846D0/
      DATA HALFPI /1.5707963267948966D0/
      DATA TOL,EPSLN /1.0D-7,1.0D-10/
      DATA ZERO,HALF,ONE /0.0D0,0.5D0,1.0D0/
      DATA SWITCH /0/
!
! ......................................................................
!       .  INITIALIZATION OF PROJECTION PARAMETERS (FILE INPUT)  .
! ......................................................................
!
      ENTRY IF20Z0 (INFILE,DATA)
!
      IERROR = 0
      READ (INFILE,END=180) ZONE,BUFF
      IF (SWITCH.NE.0 .AND. SWITCH.EQ.ZONE) RETURN
  020 MODE = 0
      IF (BUFF(13) .NE. ZERO) MODE = 1
      IF (BUFF(1) .LE. ZERO) GO TO 100
      A = BUFF(1)
      B = BUFF(2)
      IF (B .GT. ZERO) GO TO 040
      E = ZERO
      ES = ZERO
      GO TO 120
  040 IF (B .GT. ONE) GO TO 060
      E = DSQRT (B)
      ES = B
      GO TO 120
  060 ES = ONE - (B / A) ** 2
      E = DSQRT (ES)
      GO TO 120
  100 A = AZ
      E = EZ
      ES = ESZ
  120 KS0 = BUFF(3)
      LAT0 = PAKRZ0 (BUFF(6))
      X0 = BUFF(7)
      Y0 = BUFF(8)
      SINPH0 = DSIN (LAT0)
      COSPH0 = DCOS (LAT0)
      CON = ONE - ES * SINPH0 * SINPH0
      COM = DSQRT (ONE - ES)
      BL = DSQRT (ONE + ES * COSPH0 ** 4 / (ONE - ES))
      AL = A * BL * KS0 * COM / CON
      TS0 = TSFNZ0 (E,LAT0,SINPH0)
      CON = DSQRT (CON)
      D = BL * COM / (COSPH0 * CON)
      F = D + DSIGN (DSQRT (DMAX1 ((D * D - ONE), 0.0D0)) , LAT0)
      EL = F * TS0 ** BL
      IF (IPPARM .EQ. 0) PRINT 2000, A,ES,KS0
 2000 FORMAT (' INITIALIZATION PARAMETERS (OBLIQUE MERCATOR ''HOTINE''', &
              ' PROJECTION)'/                                            & 
              ' SEMI-MAJOR AXIS OF ELLIPSOID =',F12.2,' METERS'/         &
              ' ECCENTRICITY SQUARED         =',F12.9/                   &
              ' SCALE AT CENTER              =',F12.9)
      IF (MODE .EQ. 0) GO TO 140
      ALPHA = PAKRZ0 (BUFF(4))
      LONC = PAKRZ0 (BUFF(5))
      G = HALF * (F - ONE / F)
      GAMMA =  DASIN (DSIN (ALPHA) / D)
      LON0 = LONC -  DASIN (G * DTAN (GAMMA)) / BL
!
! LIST INITIALIZATION PARAMETERS (CASE B).
!
      CALL RADDZ0 (ALPHA,ANGS2(1,1))
      CALL RADDZ0 (LONC,ANGS2(1,2))
      CALL RADDZ0 (LAT0,ANGS2(1,3))
!     IF (IPPARM .EQ. 0) PRINT 2010, ANGS2
      IF (IPPARM .EQ. 0) PRINT 2010, angs2(1,1),angs2(2,1),angs2(3,1),  &
          rangs1,angs2(1,2),angs2(2,2),angs2(3,2),rangs2,angs2(1,3),    &
          angs2(2,3),angs2(3,3),rangs3,X0,Y0
 2010 FORMAT (' AZIMUTH OF CENTRAL LINE      = ',A1,2I3,F7.3/  &
              ' LONGITUDE OF ORIGIN          = ',A1,2I3,F7.3/  &
              ' LATITUDE OF ORIGIN           = ',A1,2I3,F7.3)
      CON = DABS (LAT0)
      IF (CON.GT.EPSLN .AND. DABS(CON - HALFPI).GT.EPSLN) GO TO 160
      IF (IPEMSG .EQ. 0) PRINT 2020
 2020 FORMAT (' ERROR PJ20Z0'/     &
              ' INPUT DATA ERROR')
      IERROR = 201
      RETURN
  140 LON1 = PAKRZ0 (BUFF(9))
      LAT1 = PAKRZ0 (BUFF(10))
      LON2 = PAKRZ0 (BUFF(11))
      LAT2 = PAKRZ0 (BUFF(12))
      SINPHI = DSIN (LAT1)
      TS1 = TSFNZ0 (E,LAT1,SINPHI)
      SINPHI = DSIN (LAT2)
      TS2 = TSFNZ0 (E,LAT2,SINPHI)
      H = TS1 ** BL
      L = TS2 ** BL
      F = EL / H
      G = HALF * (F - ONE / F)
      J = (EL * EL - L * H) / (EL * EL + L * H)
      P = (L - H) / (L + H)
      CALL RADDZ0 (LON2,ANGS1(1,3))
      DLON = LON1 - LON2
      IF (DLON .LT. -PI) LON2 = LON2 - 2.D0 * PI
      IF (DLON .GT.  PI) LON2 = LON2 + 2.D0 * PI
      DLON = LON1 - LON2
      LON0 = HALF * (LON1 + LON2) - DATAN (J * DTAN (HALF * BL * DLON) / P) / BL
      DLON = ADJLZ0 (LON1 - LON0)
      GAMMA = DATAN (DSIN (BL * DLON) / G)
      ALPHA =  DASIN (D * DSIN (GAMMA))
      CALL RADDZ0 (LON1,ANGS1(1,1))
      CALL RADDZ0 (LAT1,ANGS1(1,2))
!     CALL RADDZ0 (LON2,ANGS1(1,3))
      CALL RADDZ0 (LAT2,ANGS1(1,4))
      CALL RADDZ0 (LAT0,ANGS1(1,5))
!     IF (IPPARM .EQ. 0) PRINT 2030, ANGS1
      IF (IPPARM .EQ. 0) PRINT 2030,                  &
         angs1(1,1),angs1(2,1),angs1(3,1),rangs11,    &
         angs1(1,2),angs1(2,2),angs1(3,2),rangs12,    & 
         angs1(1,3),angs1(2,3),angs1(3,3),rangs13,    &
         angs1(1,4),angs1(2,4),angs1(3,4),rangs14,    &
         angs1(1,5),angs1(2,5),angs1(3,5),rangs15
 2030 FORMAT (' LONGITUDE OF 1ST POINT       = ',A1,2I3,F7.3/  &
              ' LATITUDE OF 1ST POINT        = ',A1,2I3,F7.3/  &
              ' LONGITUDE OF 2ND POINT       = ',A1,2I3,F7.3/  &
              ' LATITUDE OF 2ND POINT        = ',A1,2I3,F7.3/  &
              ' LATITUDE OF ORIGIN           = ',A1,2I3,F7.3)
      IF (DABS(LAT1 - LAT2) .LE. EPSLN) GO TO 150
      CON = DABS (LAT1)
      IF (CON.LE.EPSLN .OR. DABS(CON - HALFPI).LE.EPSLN) GO TO 150
      IF (DABS(DABS(LAT0) - HALFPI) .GT. EPSLN) GO TO 160
  150 IF (IPEMSG .EQ. 0) PRINT 2020
      IERROR = 202
      RETURN
  160 SINGAM = DSIN (GAMMA)
      COSGAM = DCOS (GAMMA)
      SINALF = DSIN (ALPHA)
      COSALF = DCOS (ALPHA)
      IF (IPEMSG .EQ. 0) PRINT 2040, X0,Y0
 2040 FORMAT (' FALSE EASTING                =',F12.2,' METERS'/  &
              ' FALSE NORTHING               =',F12.2,' METERS')
      DATA(1) = A
      DATA(2) = ES
      SWITCH = ZONE
      RETURN
  180 IF (IPEMSG .EQ. 0) PRINT 2050
 2050 FORMAT (' ERROR PJ20Z0'/                      &
              ' MISSING PROJECTION PARAMETERS')
      IERROR = 203
      RETURN
!
! ......................................................................
!      .  INITIALIZATION OF PROJECTION PARAMETERS (ENTRY INPUT)  .
! ......................................................................
!
      ENTRY IS20Z0 (ZZONE,DATA)
      zone = zzone
!
      IERROR = 0
      IF (SWITCH.NE.0 .AND. SWITCH.EQ.ZONE) RETURN
      DO 200 I = 1,13
      BUFF(I) = DATA(I)
  200 CONTINUE
      GO TO 020
!
! ......................................................................
!                      .  FORWARD TRANSFORMATION  .
! ......................................................................
!
      ENTRY PF20Z0 (GEOG,PROJ)
!
      IERROR = 0
      IF (SWITCH .NE. 0) GO TO 220
      IF (IPEMSG .EQ. 0) PRINT 2050
      IERROR = 204
      RETURN
  220 SINPHI = DSIN (GEOG(2))
      DLON = ADJLZ0 (GEOG(1) - LON0)
      VL = DSIN (BL * DLON)
      IF (DABS(DABS(GEOG(2)) - HALFPI) .GT. EPSLN) GO TO 230
      UL = SINGAM * DSIGN (ONE , GEOG(2))
      US = AL * GEOG(2) / BL
      GO TO 250
  230 TS = TSFNZ0 (E,GEOG(2),SINPHI)
      Q = EL / TS ** BL
      S = HALF * (Q - ONE / Q)
      T = HALF * (Q + ONE / Q)
      UL = (S * SINGAM - VL * COSGAM) / T
      CON = DCOS (BL * DLON)
      IF (DABS(CON) .LT. TOL) GO TO 240
      US = AL * DATAN ((S * COSGAM + VL * SINGAM) / CON) / BL
      IF (CON .LT. ZERO) US = US + PI * AL / BL
      GO TO 250
  240 US = AL * BL * DLON
  250 IF (DABS(DABS(UL) - ONE) .GT. EPSLN) GO TO 255
      IF (IPEMSG .EQ. 0) PRINT 2060
 2060 FORMAT (' ERROR PJ20Z0'/                  &
              ' POINT PROJECTS INTO INFINITY')
      IERROR = 205
      RETURN
  255 VS = HALF * AL * DLOG ((ONE - UL) / (ONE + UL)) / BL
  260 PROJ(1) = X0 + VS * COSALF + US * SINALF
      PROJ(2) = Y0 + US * COSALF - VS * SINALF
      RETURN
!
! ......................................................................
!                      .  INVERSE TRANSFORMATION  .
! ......................................................................
!
      ENTRY PI20Z0 (PROJ,GEOG)
!
      IERROR = 0
      IF (SWITCH .NE. 0) GO TO 280
      IF (IPEMSG .EQ. 0) PRINT 2050
      IERROR = 206
      RETURN
  280 X = PROJ(1) - X0
      Y = PROJ(2) - Y0
      VS = X * COSALF - Y * SINALF
      US = Y * COSALF + X * SINALF
      Q = DEXP (- BL * VS / AL)
      S = HALF * (Q - ONE / Q)
      T = HALF * (Q + ONE / Q)
      VL = DSIN (BL * US / AL)
      UL = (VL * COSGAM + S * SINGAM) / T
      IF (DABS (DABS (UL) - ONE) .GE. EPSLN) GO TO 300
      GEOG(1) = LON0
      GEOG(2) = DSIGN (HALFPI , UL)
      RETURN
  300 CON = ONE / BL
      TS = (EL / DSQRT ((ONE + UL) / (ONE - UL))) ** CON
      GEOG(2) = PHI2Z0 (E,TS)
      CON = DCOS (BL * US / AL)
      LON = LON0 - DATAN2 ((S * COSGAM - VL * SINGAM) , CON) / BL
      GEOG(1) = ADJLZ0 (LON)
      RETURN
!
      END

! **********************************************************************
! ** U.S.G.S. GENERAL MAP PROJECTION PACKAGE ...... EROS DATA CENTER **
! **          MATHEMATICAL ANALYSIS BY JOHN SNYDER                   **
! **                                                  JUNE,1982      **
! *********************************************************************
!                    *  SPACE OBLIQUE MERCATOR  *
! *********************************************************************
!
      SUBROUTINE PJ21Z0
!
      IMPLICIT REAL*8 (A-Z)
      REAL RANGS1,RANGS2
      INTEGER*4 IERROR,IPEMSG,IPPARM,ANGS2
      INTEGER*4 SWITCH,I,ZONE,INFILE,N,L
      integer *4 zzone
      COMMON/PJ21/LON0,A,B,A2,A4,C1,C3,Q,T,U,W,XJ,P21,SA,CA,ES,S,START
      COMMON /ELLPZ0/ AZ,EZ,ESZ,E0Z,E1Z,E2Z,E3Z
! ********************** X0,Y0,GAMMA,LON0,AL,BL,EL *********************
      COMMON /ERRMZ0/ IERROR
      COMMON /PRINZ0/ IPEMSG,IPPARM
!      COMMON /WORKZ0/ BUFF(15),ANGS2(4,2)
      COMMON /WORKZ0/ BUFF(15)
      COMMON /WK21Z0/ ANGS2(4,2)
      DIMENSION DATA(1),GEOG(1),PROJ(1)
      EQUIVALENCE (RANGS1,ANGS2(4,1))
      EQUIVALENCE (RANGS2,ANGS2(4,2))
      DATA PI /3.14159265358979323846D0/
      DATA HALFPI /1.5707963267948966D0/
      DATA TOL,EPSLN /1.0D-7,1.0D-10/
      DATA ZERO,HALF,ONE /0.0D0,0.5D0,1.0D0/
      DATA SWITCH /0/
!
! ......................................................................
!       .  INITIALIZATION OF PROJECTION PARAMETERS (FILE INPUT)  .
! ......................................................................
!
      ENTRY IF21Z0(INFILE,DATA)
!
      IERROR = 0
      READ (INFILE,END=180) ZONE,BUFF
      IF (SWITCH.NE.0 .AND. SWITCH.EQ.ZONE) RETURN
  020 IF (BUFF(1) .LE. ZERO) GO TO 100
      A = BUFF(1)
      B = BUFF(2)
      IF (B .GT. ZERO) GO TO 040
      ES = ZERO
      GO TO 120
  040 IF (B .GT. ONE) GO TO 060
      ES = B
      GO TO 120
  060 ES = ONE - (B / A) ** 2
      GO TO 120
  100 A = AZ
      ES = ESZ
  120 CONTINUE
      IF (IPPARM .EQ. 0) PRINT 2000, A,ES
 2000 FORMAT (' INITIALIZATION PARAMETERS (SPACE OBLIQUE MERCATOR',  &
              ' PROJECTION)'/                                        &
              ' SEMI-MAJOR AXIS OF ELLIPSOID =',F12.2,' METERS'/     &
              ' ECCENTRICITY SQUARED         =',F12.9)   
      ALF = PAKRZ0 (BUFF(4))
      LON0 = PAKRZ0 (BUFF(5))
      P21=BUFF(9)/1440.D0
      START=BUFF(11)
      CALL RADDZ0 (ALF,ANGS2(1,1))
      CALL RADDZ0 (LON0,ANGS2(1,2))
!     IF (IPPARM .EQ. 0) PRINT 2010, ANGS2,BUFF(9),BUFF(10),BUFF(11)
      IF (IPPARM .EQ. 0) PRINT 2010,                  &
         ANGS2(1,1),ANGS2(2,1),ANGS2(3,1),RANGS1,     &
         ANGS2(1,2),ANGS2(2,2),ANGS2(3,2),RANGS2,     &
         BUFF(9),BUFF(10),BUFF(11)
 2010 FORMAT (' INCLINATION OF ORBIT           = ',A1,2I3,F7.3/  &
              ' LONGITUDE OF ASCENDING ORBIT   = ',A1,2I3,F7.3/  &
              ' PERIOD OF SATELLITE REVOLUTION = ',F15.10/       &
              ' LANDSAT RATIO                  = ',F15.10/       &
              ' LANDSAT END OF PATH FLAG       = ',F15.10)
      CA=DCOS(ALF)
      IF (DABS(CA).LT.1.D-9) CA=1.D-9
      SA=DSIN(ALF)
      E2C=ES*CA*CA
      E2S=ES*SA*SA
      W=((1.D0-E2C)/(1.D0-ES))**2-1.D0
      Q = E2S / (1.D0-ES)
      T = (E2S*(2.D0-ES)) / (1.D0-ES)**2
      U= E2C / (1.D0-ES)
      XJ = (1.D0-ES)**3
      DATA(1) = A
      DATA(2) = ES
      DLAM=0.0D0
      CALL SE21Z0(FB,FA2,FA4,FC1,FC3,DLAM)
      SUMA2=FA2
      SUMA4=FA4
      SUMB=FB
      SUMC1=FC1
      SUMC3=FC3
      DO 206 I=9,81,18
      DLAM=I
      CALL SE21Z0(FB,FA2,FA4,FC1,FC3,DLAM)
      SUMA2=SUMA2+4.D0*FA2
      SUMA4=SUMA4+4.D0*FA4
      SUMB=SUMB+4.D0*FB
      SUMC1=SUMC1+4.D0*FC1
      SUMC3=SUMC3+4.D0*FC3
  206 CONTINUE
      DO 207 I=18,72,18
      DLAM=I
      CALL SE21Z0(FB,FA2,FA4,FC1,FC3,DLAM)
      SUMA2=SUMA2+2.D0*FA2
      SUMA4=SUMA4+2.D0*FA4
      SUMB=SUMB+2.D0*FB
      SUMC1=SUMC1+2.D0*FC1
      SUMC3=SUMC3+2.D0*FC3
  207 CONTINUE
      DLAM=90.D0
      CALL SE21Z0(FB,FA2,FA4,FC1,FC3,DLAM)
      SUMA2=SUMA2+FA2
      SUMA4=SUMA4+FA4
      SUMB=SUMB+FB
      SUMC1=SUMC1+FC1
      SUMC3=SUMC3+FC3
      A2=SUMA2/30.D0
      A4=SUMA4/60.D0
      B=SUMB/30.D0
      C1=SUMC1/15.D0
      C3=SUMC3/45.D0
 1002 FORMAT(1X,'A FOURIER TERMS (RAD): ',2F15.10)
 1004 FORMAT(1X,'B FOURIER TERM (RAD): ',F15.10)
 1003 FORMAT(1X,'C FOURIER TERMS (RAD): ',2F15.10)
      SWITCH = ZONE
      RETURN
  180 IF (IPEMSG .EQ. 0) PRINT 2050
 2050 FORMAT (' ERROR PJ21Z0'/                     &
              ' MISSING PROJECTION PARAMETERS')
      IERROR = 211
      RETURN
!
! ......................................................................
!      .  INITIALIZATION OF PROJECTION PARAMETERS (ENTRY INPUT)  .
! ......................................................................
!
      ENTRY IS21Z0(ZZONE,DATA)
      zone = zzone
!
      IERROR = 0
      IF (SWITCH.NE.0 .AND. SWITCH.EQ.ZONE) RETURN
      DO 200 I = 1,13
      BUFF(I) = DATA(I)
  200 CONTINUE
      GO TO 020
!
! ......................................................................
!                      .  FORWARD TRANSFORMATION  .
! ......................................................................
!
      ENTRY PF21Z0 (GEOG,PROJ)
!
      IERROR = 0
      IF (SWITCH .NE. 0) GO TO 220
      IF (IPEMSG .EQ. 0) PRINT 2050
      IERROR = 212
      RETURN
  220 CONTINUE
      CONV=1.D-7
      DLAT=GEOG(2)
      DLON=GEOG(1)-LON0
!     TEST FOR LATITUDE AND LONGITUDE APPROACHING 90 DEGREES
!
      IF (DLAT.GT.1.570796D0) DLAT=1.570796D0
!
      IF (DLAT.LT.-1.570796D0) DLAT=-1.570796D0
      RADLT=DLAT
      RADLN=DLON
      IF(DLAT.GE.0.0D0) TLAMP=PI/2.0D0
      IF(START .NE. 0.0D0)TLAMP=2.5D0*PI
      IF(DLAT.LT.0.0D0) TLAMP=1.5D0*PI
      N=0
  230 SAV=TLAMP
      L=0
      XLAMP=RADLN+P21*TLAMP
      AB1=DCOS(XLAMP)
      IF(DABS(AB1).LT.CONV) XLAMP=XLAMP-1.D-7
      IF(AB1.GE.0.D0) SCL=1.D0
      IF(AB1.LT.0.D0) SCL=-1.D0
      AB2=TLAMP-(SCL)*DSIN(TLAMP)*HALFPI
  240 XLAMT=RADLN+P21*SAV
      C=DCOS(XLAMT)
      IF (DABS(C).LT.1.D-7) XLAMT=XLAMT-1.D-7
      XLAM=(((1.D0-ES)*DTAN(RADLT)*SA)+DSIN(XLAMT)*CA)/C
      TLAM=DATAN(XLAM)
      TLAM=TLAM+AB2
      ABS=DABS(SAV)-DABS(TLAM)
      IF(DABS(ABS).LT.CONV) GO TO 250
      L=L+1
      IF (L .GT. 50) GO TO 260
      SAV=TLAM
      GO TO 240
!
!     ADJUST FOR CONFUSION AT BEGINNING AND END OF LANDSAT ORBITS
!
  250 RLM=PI*BUFF(10)
      RLM2=RLM+2.0D0*PI
      N=N+1
      IF(N.GE.3) GO TO 300
      IF(TLAM.GT.RLM.AND.TLAM.LT.RLM2) GO TO 300
      IF(TLAM.LE.RLM)TLAMP=2.5D0*PI
      IF(TLAM.GE.RLM2) TLAMP=HALFPI
      GO TO 230
  260 IF (IPEMSG .EQ. 0) PRINT 50, SAV,TLAM,DLAT,DLON
   50 FORMAT(1X,'  50 ITERATIONS WITHOUT CONV; SAV/TLAM ',5X,2F15.10/  &
             2F25.10)
      IERROR = 214
      RETURN
  300 CONTINUE
!
!     TLAM COMPUTED - NOW COMPUTE TPHI
!
      DS=DSIN(TLAM)
      DD=DS*DS
      DP=DSIN(RADLT)
      TPHI=DASIN(((1.D0-ES)*CA*DP-SA*DCOS(RADLT)*DSIN(XLAMT))/DSQRT(1.D0-ES*DP*DP))
!
!     COMPUTE X AND Y
!
      XTAN = (PI/4.0D0) + (TPHI/2.0D0)
      TANLG = DLOG(DTAN(XTAN))
      SD=DSIN(TLAM)
      SDSQ=SD*SD
      S=P21*SA*DCOS(TLAM)*DSQRT((1.D0+T*SDSQ)/((1.D0+W*SDSQ)*(1.D0+Q*SDSQ)))
      D=DSQRT(XJ*XJ+S*S)
      X=B*TLAM+A2*DSIN(2.D0*TLAM)+A4*DSIN(4.D0*TLAM)-TANLG*S/D
      PROJ(1)=A*X
      Y=C1*SD+C3*DSIN(3.D0*TLAM)+TANLG*XJ/D
      PROJ(2)=A*Y
      RETURN
!
! ......................................................................
!                      .  INVERSE TRANSFORMATION  .
! ......................................................................
!
      ENTRY PI21Z0 (PROJ,GEOG)
!
      IERROR = 0
      IF (SWITCH .NE. 0) GO TO 280
      IF (IPEMSG .EQ. 0) PRINT 2050
      IERROR = 213
      RETURN
  280 CONTINUE
!
!     COMPUTES TRANSFORMED LAT/LON AND GEODETIC
!     LAT/LON GIVEN X-Y
!
!
!     BEGIN INVERSE COMPUTATION WITH APPROXIMATION FOR TLON. SOLVE
!     FOR TRANSFORMED LONG.
!
      X=PROJ(1)
      Y=PROJ(2)
      TLON= X/(A*B)
      CONV=1.D-9
      INUMB=0
  813 SAV=TLON
      SD=DSIN(TLON)
      SDSQ=SD*SD
      S=P21*SA*DCOS(TLON)*DSQRT((1.D0+T*SDSQ)/((1.D0+W*SDSQ)*(1.D0+Q*SDSQ)))
      BLON=(X/A)+(Y/A)*S/XJ-A2*DSIN(2.D0*TLON)-A4*DSIN(4.D0*TLON)-(S/XJ) *  &
           (C1*DSIN(TLON)+C3*DSIN(3.D0*TLON))
      TLON=BLON/B
      DIF=TLON-SAV
      IF(DABS(DIF).LT.CONV) GO TO 814
      INUMB=INUMB+1
      IF(INUMB.LT.50) GO TO 813
      IERROR=215
      RETURN
!
!     COMPUTE TRANSFORMED LAT.
!
  814 CONTINUE
      ST=DSIN(TLON)
      DEFAC=DEXP(DSQRT(1.D0+S*S/XJ/XJ)*(Y/A-C1*ST-C3*DSIN(3.D0*TLON)))
      ACTAN=DATAN(DEFAC)
      TLAT=2.0D0*(ACTAN-(PI/4.0D0))
!
!     COMPUTE GEODETIC LONGITUDE
!
      DD=ST*ST
      IF(DABS(DCOS(TLON)).LT.1.D-7) TLON=TLON-1.D-7
      BIGK=DSIN(TLAT)
      BIGK2=BIGK*BIGK
      XLAMT=DATAN(((1.D0-BIGK2/(1.D0-ES))*DTAN(TLON)*CA-BIGK*SA*DSQRT((1.D0+Q*DD)* &
           (1.D0-BIGK2)-BIGK2*U)/DCOS(TLON))/(1.D0-BIGK2*(1.D0+U)))
!
!     CORRECT INVERSE  QUADRANT
!
      IF(XLAMT.GE.0.D0) SL=1.D0
      IF(XLAMT.LT.0.D0) SL=-1.D0
      IF(DCOS(TLON).GE.0.D0) SCL=1.D0
      IF(DCOS(TLON).LT.0.D0) SCL=-1.0D0
      XLAMT=XLAMT-((PI/2.D0)*(1.D0-SCL)*SL)
      DLON=XLAMT-P21*TLON
!
!     COMPUTE GEODETIC LATITUDE
!
      IF(DABS(SA) .LT.1.D-7)DLAT=DASIN(BIGK/DSQRT((1.D0-ES)*(1.D0-ES)+ES*BIGK2))
      IF(DABS(SA) .LT.1.D-7)GO TO 1
      DLAT=DATAN((DTAN(TLON)*DCOS(XLAMT)-CA*DSIN(XLAMT))/((1.D0-ES)*SA))
    1 CONTINUE
      GEOG(1)=ADJLZ0(DLON+LON0)
      GEOG(2)=DLAT
      RETURN
      END

!                   QSFNZ0
! **********************************************************************
! ** U.S.G.S. GENERAL MAP PROJECTION PACKAGE ...... DR. A. A. ELASSAL **
! ** MODULE I                VERSION 1.0.0            NOVEMBER 1,1980 **
! **********************************************************************
      DOUBLE PRECISION FUNCTION QSFNZ0 (ECCENT,SINPHI,COSPHI)
!
! FUNCTION TO COMPUTE CONSTANT (SMALL Q).
!
      IMPLICIT REAL*8 (A-Z)
      DATA HALF,ONE,TWO /0.5D0,1.0D0,2.0D0/
      DATA EPSLN /1.0D-7/
!
      IF (ECCENT .LT. EPSLN) GO TO 020
      CON = ECCENT * SINPHI
      QSFNZ0 = (ONE - ECCENT * ECCENT) * (SINPHI / (ONE - CON * CON) - &
               (HALF / ECCENT) * DLOG ((ONE - CON) / (ONE + CON)))
      RETURN
!
  020 QSFNZ0 = TWO * SINPHI
      RETURN
      END

!                   RADDZ0
! **********************************************************************
! ** U.S.G.S. GENERAL MAP PROJECTION PACKAGE ...... DR. A. A. ELASSAL **
! ** MODULE I                VERSION 1.0.0            NOVEMBER 1,1980 **
! **********************************************************************
      SUBROUTINE RADDZ0 (RAD,IDMS)
!
! SUBROUTINE TO CONVERT ANGLE FROM RADIANS TO DMS.
! SECONDS (IDMS) .
!
      REAL*8 RAD,CON,RADSEC
      DIMENSION IDMS(1)
      EQUIVALENCE (FLOT,INTG)
      DATA RADSEC,IBLANK,NEG /206264.80625D0,' ','-'/
      DATA ZERO /0.0D0/
!
! DETERMINE THE SIGN OF THE ANGLE.
!
      IDMS(1) = IBLANK
      IF (RAD .LT. ZERO) IDMS(1) = NEG
!
! CONVERT THE ANGLE TO THOUSANDTH OF SECONDS.
!
      CON = DABS(RAD) * RADSEC
      ISEC = CON
!
! COMPUTE DEGREES PART OF THE ANGLE.
!
      INTG = ISEC / 3600
      IDMS(2) = INTG
      ISEC = INTG * 3600
      CON = CON - DFLOAT(ISEC)
      ISEC = CON
!
! COMPUTE MINUTES PART OF THE ANGLE.
!
      IDMS(3) = ISEC / 60
      ISEC = IDMS(3) * 60
      CON = CON - DFLOAT(ISEC)
!
! COMPUTE SECONDS PART OF THE ANGLE.
!
      FLOT = CON
      IDMS(4) = INTG
!
      RETURN
      END

! *********************************************************************
! ** U.S.G.S  GENERAL MAP PROJECTION PACKAGE ******* EROS DATA CENTER *
! ** MATHEMATICAL ANALYSIS BY JOHN SNYDER            JUNE,1982 ********
! *********************************************************************
!
!    SERIES TO CALCULATE A,B,C COEFFICIENTS TO CONVERT FROM
!    TRANSFORM LATITUDE,LONGITUDE TO SPACE OBLIQUE MERCATOR(SOM)
!    RECTANGULAR COORDINATES
!
! *********************************************************************
      SUBROUTINE SE21Z0(FB,FA2,FA4,FC1,FC3,DLAM)
      IMPLICIT REAL*8(A-Z)
      COMMON/PJ21/LON0,A,B,A2,A4,C1,C3,Q,T,U,W,XJ,P21,SA,CA,ES,S,START
!    CONVERT DLAM TO RADIANS
      DLAM=DLAM*0.0174532925D0
      SD=DSIN(DLAM)
      SDSQ=SD*SD
      S=P21*SA*DCOS(DLAM)*DSQRT((1.D0+T*SDSQ)/((1.D0+W*SDSQ)*(1.D0+Q*SDSQ)))
      H=DSQRT((1.D0+Q*SDSQ)/(1.D0+W*SDSQ))*(((1.D0+W*SDSQ)/((1.D0+Q*SDSQ) * (1.D0+Q*SDSQ)))-P21*CA)
      SQ=DSQRT(XJ*XJ+S*S)
      FB=(H*XJ-S*S)/SQ
      FA2=FB*DCOS(2.D0*DLAM)
      FA4=FB*DCOS(4.D0*DLAM)
      FC=S*(H+XJ)/SQ
      FC1=FC*DCOS(DLAM)
      FC3=FC*DCOS(3.D0*DLAM)
      RETURN
      END

!****                                                              *****
! ** U.S.G.S. GENERAL MAP PROJECTION PACKAGE ...... JOHN F. WAANANEN  **
! ** MODULE I                VERSION 1.0.0             APRIL 2, 1981  **
! ****                                                             *****
!
      SUBROUTINE SPHDZ0(ISPH,PARM)
!     SUBROUTINE TO COMPUTE SPHEROID PARAMETERS
!     SUBROUTINE SPHDZ0 IS REQUIRED FOR PROGRAMS NO. L176 AND NO. L177
!
!     ISPH IS THE SPHEROID CODE FROM THE FOLLOWING LIST:
!     0 = CLARKE 1866           1 = CLARKE 1880
!     2 = BESSEL                3 = NEW INTERNATIONAL 1967
!     4 = INTERNATIONAL 1909    5 = WGS 72
!     6 = EVEREST               7 = WGS 66
!     8 = GRS 1980              9 = AIRY
!    10 = MODIFIED EVEREST     11 = MODIFIED AIRY
!    12 = WALBECK              13 = SOUTHEAST ASIA
!    14 = AUSTRALIAN NATIONAL  15 = KRASSOVSKY
!    16 = HOUGH                17 = MERCURY 1960
!    18 = MODIFIED MERC 1968   19 = SPHERE OF RADIUS 6370997 M
!
!    PARM IS ARRAY OF PROJECTION PARAMETERS:
!       PARM(1) IS THE SEMI-MAJOR AXIS
!       PARM(2) IS THE ECCENTRICITY SQUARED
!
!     IF ISPH IS NEGATIVE, THE DEFAULT IS RESET FROM CLARKE 1866
!     TO THE POSITIVE OF "ISPH" SPHEROID.
!
!     IF ISPH = 0 , THE DEFAULT IS RESET TO CLARKE 1866
!
! ****                                                             *****
!
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION PARM(2),AXIS(20),BXIS(20) !Alex, parm(15) to parm(2)   S
!
      COMMON/ELLPZ0/ AZ,EZ,ESZ,E0Z,E1Z,E2Z,E3Z
      COMMON/SPHRZ0/ AZZ
      COMMON/ERRMZ0/ IERROR
      COMMON/PRINZ0/ IPEMSG,IPPARM
      COMMON/PROJZ0/ IPROJ
!
      DATA AXIS/6378206.4D0,6378249.145D0,6377397.155D0,6378157.5D0,      &
       6378388.0D0,6378135.0D0,6377276.3452D0,6378145.0D0,6378137.0D0,    &
       6377563.396D0,6377304.063D0,6377341.89D0,6376896.0D0,6378155.0D0,  &
       6378160.0D0,6378245.0D0,6378270.0D0,6378166.0D0,6378150.0D0,       &
       6370997.0D0/
!
      DATA BXIS/6356583.8D0,6356514.86955D0,6356078.96284D0,              &
       6356772.2D0,6356911.94613D0,6356750.519915D0,6356075.4133D0,       &
       6356759.769356D0,6356752.31414D0,6356256.91D0,6356103.039D0,       &
       6356036.143D0,6355834.8467D0,6356773.3205D0,6356774.719D0,         &
       6356863.0188D0,6356794.343479D0,6356784.283666D0,6356768.337303D0, & 
       6370997.0D0/
!
      IF (ISPH.LT.0) GO TO 5
      IF (PARM(1).NE.0.0D0.AND.IPROJ.NE.1) RETURN
    5 JSPH = IABS(ISPH) + 1
      IF (JSPH.LE.20) GO TO 10
      IERROR = 1211
      IF (IPEMSG.EQ.0) PRINT 1
    1 FORMAT(' Spheroid code of ',I5,' reset to 0')
      ISPH = 0
       JSPH = 1
   10 A = AXIS(JSPH)
      B = BXIS(JSPH)
      ES = (A*A-B*B)/(A*A)
!     IF ISPH LE 0 THEN RESET DEFAULT
      IF (ISPH.GT.0) GO TO 20
      AZZ = 6370997.0D0
      EZ  = DSQRT(ES)
      E0Z = E0FNZ0(ES)
      E1Z = E1FNZ0(ES)
      E2Z = E2FNZ0(ES)
      E3Z = E3FNZ0(EZ)
      AZ  = A
      ESZ = ES
      IF (ES.EQ.0.0D0) AZZ=A
      IF (IPPARM.NE.0) GO TO 20
      PRINT 12342,isph
      PRINT 12343,az
      PRINT 12344,ez
      PRINT 12345,esz
      PRINT 12346,e0z
      PRINT 12347,e1z
      PRINT 12348,e2z
      PRINT 12349,e3z
      PRINT 12350,azz
12342 format(' SPHEROID DEFAULT SET TO NUMBER ',i2)
12343 format(' SEMI-MAJOR AXIS  = ',D24.15)
12344 format(' ECCENTRICITY     = ',d24.15)
12345 format(' ECCENT. SQUARED  = ',d24.15)
12346 format(' E0               = ',D24.15)
12347 format(' E1               = ',D24.15)
12348 format(' E2               = ',D24.15)
12349 format(' E3               = ',D24.15)
12350 format(' REFERENCE SPHERE = ',d24.15)
!
   20 PARM(1) = A
      PARM(2) = ES
      RETURN
      END

!                   TSFNZ0
! **********************************************************************
! ** U.S.G.S. GENERAL MAP PROJECTION PACKAGE ...... DR. A. A. ELASSAL **
! ** MODULE I                VERSION 1.0.0            NOVEMBER 1,1980 **
! **********************************************************************
      DOUBLE PRECISION FUNCTION TSFNZ0 (ECCENT,PHI,SINPHI)
!
! FUNCTION TO COMPUTE CONSTANT (SMALL T).
!
      IMPLICIT REAL*8 (A-Z)
      DATA HALF,ONE /0.5D0,1.0D0/
      DATA HALFPI /1.5707963267948966D0/
!
      CON = ECCENT * SINPHI
      COM = HALF * ECCENT
      CON = ((ONE - CON) / (ONE + CON)) ** COM
      TSFNZ0 = DTAN (HALF * (HALFPI - PHI)) / CON
!
      RETURN
      END

!                   UNTFZ0
! **********************************************************************
! ** U.S.G.S. GENERAL MAP PROJECTION PACKAGE ...... DR. A. A. ELASSAL **
! ** MODULE I                VERSION 1.0.0            NOVEMBER 1,1980 **
! **********************************************************************
      SUBROUTINE UNTFZ0 (INUNIT,IOUNIT,FACTOR,IFLG)
!
! SUBROUTINE TO DETERMINE CONVERGENCE FACTOR BETWEEN TWO LINEAL UNITS
!
! * INPUT ........
! * INUNIT * UNIT CODE OF SOURCE.
! * IOUNIT * UNIT CODE OF TARGET.
!
! * OUTPUT .......
! * FACTOR * CONVERGENCE FACTOR FROM SOURCE TO TARGET.
! * IFLG   * RETURN FLAG = 0 , NORMAL RETURN.
!                        = 1 , ABNORMAL RETURN.
!
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION FACTRS(5,5)
      DATA ZERO,MAXUNT /0.0D0,5/
      DATA FACTRS /0.1000000000000000D01 , 0.0000000000000000D00 ,  &
                   0.0000000000000000D00 , 0.2062648062470963D06 ,  &
                   0.5729577951308231D02 ,                          &
                   0.0000000000000000D00 , 0.1000000000000000D01 ,  &
                   0.3048006096012192D00 , 0.0000000000000000D00 ,  &
                   0.0000000000000000D00 ,                          &
                   0.0000000000000000D00 , 0.3280833333333333D01 ,  &
                   0.1000000000000000D01 , 0.0000000000000000D00 ,  &
                   0.0000000000000000D00 ,                          &
                   0.4848136811095360D-5 , 0.0000000000000000D00 ,  &
                   0.0000000000000000D00 , 0.1000000000000000D01 ,  &
                   0.2777777777777778D-3 ,                          &
                   0.1745329251994330D-1 , 0.0000000000000000D00 ,  &
                   0.0000000000000000D00 , 0.3600000000000000D04 ,  &
                   0.1000000000000000D01 /
!
      IF (INUNIT.LT.0 .OR. INUNIT.GE.MAXUNT) GO TO 020
      IF (IOUNIT.GE.0 .AND. IOUNIT.LT.MAXUNT) GO TO 040
  020 IFLG = 1
      RETURN
  040 FACTOR = FACTRS(IOUNIT+1 , INUNIT+1)
      IF (FACTOR .EQ. ZERO) GO TO 020
      IFLG = 0
      RETURN
!
      END
