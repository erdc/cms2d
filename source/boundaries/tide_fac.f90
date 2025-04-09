      subroutine tide_fac(IYR,IMO,IDAY,BHR,XDAYS,NCONST,NAME,SPD,FACNOD,EQARG)
! PROGRAM TO COMPUTE NODAL FACTORS AND EQUILIBRIUM ARGUEMENTS
      use diag_lib, only: diag_print_error
       
      !Input
      integer :: IDAY,IMO,IYR
      real :: XDAYS,BHR
      !Output
      CHARACTER NAME(NCONST)*8
      real :: SPD(NCONST),FACNOD(NCONST),EQARG(NCONST)
      
      !Internal
      INTEGER NCNST
      PARAMETER(NCNST=37)
      CHARACTER CNAME(NCNST)*8      
      COMMON /CNSNAM/ CNAME
      REAL NODFAC,MONTH,GRTERM,SPEED,P
      DIMENSION NCON(NCNST)
      COMMON /CNST/ NODFAC(NCNST),GRTERM(NCNST),SPEED(NCNST),P(NCNST)
      REAL YR,DAY,HRM,RHRS,DAYJ
      
      !OPEN(UNIT=11,FILE='tide_fac.out',STATUS='UNKNOWN')

      !WRITE(*,*) 'ENTER LENGTH OF RUN TIME (DAYS)'
      !READ(*,*) XDAYS
      RHRS=XDAYS*24.

      !WRITE(*,*)' ENTER START TIME - BHR,IDAY,IMO,IYR (IYR e.g. 1992)'
      !READ(*,*) BHR,IDAY,IMO,IYR
      YR=IYR
      MONTH=IMO
      DAY=IDAY
      HRM=BHR+RHRS/2.
      !WRITE(11,10) BHR,IDAY,IMO,IYR
      WRITE(*,10) BHR,IDAY,IMO,IYR
  10  FORMAT(' TIDAL FACTORS STARTING: ', ' HR-',F5.2,',  DAY-',I3,',  MONTH-',I3,'  YEAR-',I5,/)
      WRITE(*,11) XDAYS
      !WRITE(11,11) XDAYS
  11  FORMAT(' FOR A RUN LASTING ',F8.2,' DAYS',//)

!   DETERMINE THE JULIAN TIME AT BEGINNING AND MIDDLE OF RECORD
      DAYJ=DAYJUL(YR,MONTH,DAY)

!   DETERMINE NODE FACTORS AT MIDDLE OF RECORD
      CALL NFACS(YR,DAYJ,HRM)

!   DETERMINE GREENWICH EQUIL. TERMS AT BEGINNING OF RECORD
      CALL GTERMS(YR,DAYJ,BHR,DAYJ,HRM)

!   GET MATCHING CONSTITUENTS       
      DO IC=1,NCONST
        DO JC=1,NCNST
          IF(NAME(IC).EQ.CNAME(JC))EXIT
        ENDDO
        IF(JC.GT.NCNST) call diag_print_error('Constituent not found')
        NAME(IC)=CNAME(JC)
        SPD(IC)=SPEED(JC)
        FACNOD(IC)=NODFAC(JC)
        EQARG(IC)=GRTERM(JC)
      ENDDO
          
      RETURN
      NUMCON=8
      NCON(1)=4
      NCON(2)=6
      NCON(3)=30
      NCON(4)=26
      NCON(5)=3
      NCON(6)=1
      NCON(7)=2
      NCON(8)=35

      WRITE(11,*) 'CONST   NODE     EQ ARG (ref GM)'
      WRITE(11,1300)
 1300 FORMAT(' NAME   FACTOR    (DEG) ',//)

      DO 20 NC=1,NUMCON
        IC=NCON(NC)

!  EQUILIBRIUM ARGUMENT IS REFERENCED TO THE GRENWICH MERIDIAN

        WRITE(11,2001) CNAME(IC),NODFAC(IC),GRTERM(IC)
 2001   FORMAT(1X,A4,2x,F10.5,4x,F10.2,2x,F10.4)
   20   CONTINUE

      RETURN
      end subroutine

!*********************************************************************
      SUBROUTINE NFACS(YR,DAYJ,HR)

!   CALCULATES NODE FACTORS FOR CONSTITUENT TIDAL SIGNAL

!   THE EQUATIONS USED IN THIS ROUTINE COME FROM:
!         "MANUAL OF HARMONIC ANALYSIS AND PREDICTION OF TIDES"
!         BY PAUL SCHUREMAN, SPECIAL PUBLICATION #98, US COAST
!         AND GEODETIC SURVEY, DEPARTMENT OF COMMERCE (1958).

!   IF DAYM AND HRM CORRESPOND TO MIDYEAR, THEN THIS ROUTINE
!   RETURNS THE SAME VALUES AS FOUND IN TABLE 14 OF SCHUREMAN.
!*********************************************************************

      CHARACTER*8   CST(37)
      REAL          I,N,NU,YR,HR

      COMMON/ORBITF/DS,DP,DH,DP1,DN,DI,DNU,DXI,DNUP,DNUP2,DPC
      COMMON/ CNST /FNDCST(37),EQCST(37),ACST(37),PCST(37)
      COMMON/CNSNAM/CST

!   CONSTITUENT NAMES:
      DATA CST     /'M2      ','S2      ','N2      ','K1      ',  &
                    'M4      ','O1      ','M6      ','MK3     ',  &
                    'S4      ','MN4     ','NU2     ','S6      ',  &
                    'MU2     ','2N2     ','OO1     ','LAMBDA2 ',  &
                    'S1      ','M1      ','J1      ','MM      ',  &
                    'SSA     ','SA      ','MSF     ','MF      ',  &
                    'RHO1    ','Q1      ','T2      ','R2      ',  &
                    '2Q1     ','P1      ','2SM2    ','M3      ',  &
                    'L2      ','2MK3    ','K2      ','M8      ',  &
                    'MS4     '/

!   ORBITAL SPEEDS (DEGREES/HOUR):
      DATA ACST/28.9841042,30.0,28.4397295,15.0410686,57.9682084,         &
      13.9430356,86.9523127,44.0251729,60.0,57.4238337,28.5125831,90.0,   &
      27.9682084,27.8953548,16.1391017,29.4556253,15.0,14.4966939,        &
      15.5854433,0.5443747,0.0821373,0.0410686,1.0158958,1.0980331,       &
      13.4715145,13.3986609,29.9589333,30.0410667,12.8542862,14.9589314,  &
      31.0158958,43.4761563,29.5284789,42.9271398,30.0821373,             &
      15.9364169,58.9841042/

!   NUMBER OF TIDE CYCLES PER DAY PER CONSTITUENT:
      DATA PCST/2.,2.,2.,1.,4.,1.,6.,3.,4.,4.,2.,6.,2.,2.,1.,2.,1.,1.,   &
      1.,0.,0.,0.,0.,0.,1.,1.,2.,2.,1.,1.,2.,3.,2.,3.,2.,8.,4./

      PI180=3.14159265/180.
      CALL ORBIT(YR,DAYJ,HR)
      N=DN*PI180
      I=DI*PI180
      NU=DNU*PI180
      XI=DXI*PI180
      P=DP*PI180
      PC=DPC*PI180
      SINI=SIN(I)
      SINI2=SIN(I/2.)
      SIN2I=SIN(2.*I)
      COSI2=COS(I/2.)
      TANI2=TAN(I/2.)
!   EQUATION 197, SCHUREMAN
      QAINV=SQRT(2.310+1.435*COS(2.*PC))
!   EQUATION 213, SCHUREMAN
      RAINV=SQRT(1.-12.*TANI2**2*COS(2.*PC)+36.*TANI2**4)
!   VARIABLE NAMES REFER TO EQUATION NUMBERS IN SCHUREMAN
      EQ73=(2./3.-SINI**2)/.5021
      EQ74=SINI**2/.1578
      EQ75=SINI*COSI2**2/.37988
      EQ76=SIN(2*I)/.7214
      EQ77=SINI*SINI2**2/.0164
      EQ78=(COSI2**4)/.91544
      EQ149=COSI2**6/.8758
      EQ207=EQ75*QAINV
      EQ215=EQ78*RAINV
      EQ227=SQRT(.8965*SIN2I**2+.6001*SIN2I*COS(NU)+.1006)
      EQ235=.001+SQRT(19.0444*SINI**4+2.7702*SINI**2*COS(2.*NU)+.0981)
!   NODE FACTORS FOR 37 CONSTITUENTS:
      FNDCST(1)=EQ78
      FNDCST(2)=1.0
      FNDCST(3)=EQ78
      FNDCST(4)=EQ227
      FNDCST(5)=FNDCST(1)**2
      FNDCST(6)=EQ75
      FNDCST(7)=FNDCST(1)**3
      FNDCST(8)=FNDCST(1)*FNDCST(4)
      FNDCST(9)=1.0
      FNDCST(10)=FNDCST(1)**2
      FNDCST(11)=EQ78
      FNDCST(12)=1.0
      FNDCST(13)=EQ78
      FNDCST(14)=EQ78
      FNDCST(15)=EQ77
      FNDCST(16)=EQ78
      FNDCST(17)=1.0
!** EQUATION 207 NOT PRODUCING CORRECT ANSWER FOR M1
!**SET NODE FACTOR FOR M1 = 0 UNTIL CAN FURTHER RESEARCH
      FNDCST(18)=0.
!     FNDCST(18)=EQ207
      FNDCST(19)=EQ76
      FNDCST(20)=EQ73
      FNDCST(21)=1.0
      FNDCST(22)=1.0
      FNDCST(23)=EQ78
      FNDCST(24)=EQ74
      FNDCST(25)=EQ75
      FNDCST(26)=EQ75
      FNDCST(27)=1.0
      FNDCST(28)=1.0
      FNDCST(29)=EQ75
      FNDCST(30)=1.0
      FNDCST(31)=EQ78
      FNDCST(32)=EQ149
!** EQUATION 215 NOT PRODUCING CORRECT ANSWER FOR L2
!** SET NODE FACTOR FOR L2 = 0 UNTIL CAN FURTHER RESEARCH
      FNDCST(33)=0.
!     FNDCST(33)=EQ215
      FNDCST(34)=FNDCST(1)**2*FNDCST(4)
      FNDCST(35)=EQ235
      FNDCST(36)=FNDCST(1)**4
      FNDCST(37)=EQ78
      END
      
!*********************************************************************
      SUBROUTINE GTERMS(YR,DAYJ,HR,DAYM,HRM)
!   CALCULATES EQUILIBRIUM ARGUMENTS V0+U FOR CONSTITUENT TIDE

!   THE EQUATIONS USED IN THIS ROUTINE COME FROM:
!         "MANUAL OF HARMONIC ANALYSIS AND PREDICTION OF TIDES"
!         BY PAUL SCHUREMAN, SPECIAL PUBLICATION #98, US COAST
!         AND GEODETIC SURVEY, DEPARTMENT OF COMMERCE (1958).

!   IF DAYM AND HRM CORRESPOND TO MIDYEAR, THEN THIS ROUTINE
!   RETURNS THE SAME VALUES AS FOUND IN TABLE 15 OF SCHUREMAN.
!*********************************************************************
      REAL NU,NUP,NUP2,I
      COMMON /ORBITF/DS,DP,DH,DP1,DN,DI,DNU,DXI,DNUP,DNUP2,DPC
      COMMON /CNST/ FNDCST(37),EQCST(37),ACST(37),PCST(37)
      
      PI180=3.14159265/180.
!* OBTAINING ORBITAL VALUES AT BEGINNING OF SERIES FOR V0
      CALL ORBIT(YR,DAYJ,HR)
      S=DS
      P=DP
      H=DH
      P1=DP1
      T=ANGLE(180.+HR*(360./24.))
!** OBTAINING ORBITAL VALUES AT MIDDLE OF SERIES FOR U
      CALL ORBIT(YR,DAYM,HRM)
      NU=DNU
      XI=DXI
      NUP=DNUP
      NUP2=DNUP2
!* SUMMING TERMS TO OBTAIN EQUILIBRIUM ARGUMENTS
      EQCST(1)=2.*(T-S+H)+2.*(XI-NU)
      EQCST(2)=2.*T
      EQCST(3)=2.*(T+H)-3.*S+P+2.*(XI-NU)
      EQCST(4)=T+H-90.-NUP
      EQCST(5)=4.*(T-S+H)+4.*(XI-NU)
      EQCST(6)=T-2.*S+H+90.+2.*XI-NU
      EQCST(7)=6.*(T-S+H)+6.*(XI-NU)
      EQCST(8)=3.*(T+H)-2.*S-90.+2.*(XI-NU)-NUP
      EQCST(9)=4.*T
      EQCST(10)=4.*(T+H)-5.*S+P+4.*(XI-NU)
      EQCST(11)=2.*T-3.*S+4.*H-P+2.*(XI-NU)
      EQCST(12)=6.*T
      EQCST(13)=2.*(T+2.*(H-S))+2.*(XI-NU)
      EQCST(14)=2.*(T-2.*S+H+P)+2.*(XI-NU)
      EQCST(15)=T+2.*S+H-90.-2.*XI-NU
      EQCST(16)=2.*T-S+P+180.+2.*(XI-NU)
      EQCST(17)=T
      I=DI*PI180
      PC=DPC*PI180
      TOP=(5.*COS(I)-1.)*SIN(PC)
      BOTTOM=(7.*COS(I)+1.)*COS(PC)
      Q=ARCTAN(TOP,BOTTOM,1)
      EQCST(18)=T-S+H-90.+XI-NU+Q
      EQCST(19)=T+S+H-P-90.-NU
      EQCST(20)=S-P
      EQCST(21)=2.*H
      EQCST(22)=H
      EQCST(23)=2.*(S-H)
      EQCST(24)=2.*S-2.*XI
      EQCST(25)=T+3.*(H-S)-P+90.+2.*XI-NU
      EQCST(26)=T-3.*S+H+P+90.+2.*XI-NU
      EQCST(27)=2.*T-H+P1
      EQCST(28)=2.*T+H-P1+180.
      EQCST(29)=T-4.*S+H+2.*P+90.+2.*XI-NU
      EQCST(30)=T-H+90.
      EQCST(31)=2.*(T+S-H)+2.*(NU-XI)
      EQCST(32)=3.*(T-S+H)+3.*(XI-NU)
      R=SIN(2.*PC)/((1./6.)*(1./TAN(.5*I))**2-COS(2.*PC))
      R=ATAN(R)/PI180
      EQCST(33)=2.*(T+H)-S-P+180.+2.*(XI-NU)-R
      EQCST(34)=3.*(T+H)-4.*S+90.+4.*(XI-NU)+NUP
      EQCST(35)=2.*(T+H)-2.*NUP2
      EQCST(36)=8.*(T-S+H)+8.*(XI-NU)
      EQCST(37)=2.*(2.*T-S+H)+2.*(XI-NU)
      DO IH=1,37
        EQCST(IH)=ANGLE(EQCST(IH))
      ENDDO
      
      END

!*********************************************************************
      SUBROUTINE ORBIT(YR,DAYJ,HR)

!   DETERMINATION OF PRIMARY AND SECONDARY ORBITAL FUNCTIONS

!   THE EQUATIONS PROGRAMMED HERE ARE NOT REPRESENTED BY EQUATIONS IN
!   SCHUREMAN.  THE CODING IN THIS ROUTINE DERIVES FROM A PROGRAM BY
!   THE NATIONAL OCEANIC AND ATMOSPHERIC ADMINISTRATION (NOAA).
!   HOWEVER, TABULAR VALUES OF THE ORBITAL FUNCTIONS CAN BE FOUND IN
!   TABLE 1 OF SCHUREMAN.
!*********************************************************************
      REAL I,N,NU,NUP,NUP2
      COMMON /ORBITF/DS,DP,DH,DP1,DN,DI,DNU,DXI,DNUP,DNUP2,DPC

      PI180=3.14159265/180.
      X=AINT((YR-1901.)/4.)
      DYR=YR-1900.
      DDAY=DAYJ+X-1.
!   DN IS THE MOON'S NODE (CAPITAL N, TABLE 1, SCHUREMAN)
      DN=259.1560564-19.328185764*DYR-.0529539336*DDAY-.0022064139*HR
      DN=ANGLE(DN)
      N=DN*PI180
!   DP IS THE LUNAR PERIGEE (SMALL P, TABLE 1)
      DP=334.3837214+40.66246584*DYR+.111404016*DDAY+.004641834*HR
      DP=ANGLE(DP)
      P=DP*PI180
      I=ACOS(.9136949-.0356926*COS(N))
      DI=ANGLE(I/PI180)
      NU=ASIN(.0897056*SIN(N)/SIN(I))
      DNU=NU/PI180
      XI=N-2.*ATAN(.64412*TAN(N/2.))-NU
      DXI=XI/PI180
      DPC=ANGLE(DP-DXI)
!   DH IS THE MEAN LONGITUDE OF THE SUN (SMALL H, TABLE 1)
      DH=280.1895014-.238724988*DYR+.9856473288*DDAY+.0410686387*HR
      DH=ANGLE(DH)
!   DP1 IS THE SOLAR PERIGEE (SMALL P1, TABLE 1)
      DP1=281.2208569+.01717836*DYR+.000047064*DDAY+.000001961*HR
      DP1=ANGLE(DP1)
!   DS IS THE MEAN LONGITUDE OF THE MOON (SMALL S, TABLE 1)
      DS=277.0256206+129.38482032*DYR+13.176396768*DDAY+.549016532*HR
      DS=ANGLE(DS)
      NUP=ATAN(SIN(NU)/(COS(NU)+.334766/SIN(2.*I)))
      DNUP=NUP/PI180
      NUP2=ATAN(SIN(2.*NU)/(COS(2.*NU)+.0726184/SIN(I)**2))/2.
      DNUP2=NUP2/PI180
      END

!*********************************************************************      
      FUNCTION ANGLE(ARG)
!
!*** THIS ROUTINE PLACES AN ANGLE IN 0-360 (+) FORMAT
!
!*********************************************************************
      INTEGER M
      REAL ARG,ANGLE      
      M=-IFIX(ARG/360.)
      ANGLE=ARG+FLOAT(M)*360.
      IF(ANGLE .LT. 0.) ANGLE=ANGLE+360.
      END
      
!*********************************************************************      
      FUNCTION ARCTAN(TOP,BOTTOM,KEY)
!** DETERMINE ARCTANGENT AND PLACE IN CORRECT QUADRANT
!   IF KEY EQ 0  NO QUADRANT SELECTION MADE
!   IF KEY .NE. 0 PROPER QUADRANT IS SELECTED
!*********************************************************************
      INTEGER K
      REAL BOTTOM,TOP,ARCTAN
      
   !   IF(BOTTOM .NE. 0.0) GO TO 4
   !   IF(TOP) 2,9,3
   ! 2 ARCTAN=270.
   !   RETURN
   ! 3 ARCTAN=90.
   !   RETURN
   ! 4 ARCTAN=ATAN(TOP/BOTTOM)*57.2957795
   !   IF(KEY.EQ.0) RETURN
   !   IF(TOP) 5,5,7
   ! 5 IF(BOTTOM) 6,9,8
   ! 6 ARCTAN=ARCTAN+180.
   !   RETURN
   ! 7 IF(BOTTOM) 6,3,10
   ! 8 ARCTAN=ARCTAN+360.
   !   RETURN
   ! 9 ARCTAN=0.
   !10 RETURN
      
      !MEB change for Gnu fortran issue
      if (BOTTOM .eq. 0.0) then
        if (TOP .lt. 0) ARCTAN = 270.                !2
        if (TOP .eq. 0) ARCTAN = 0.                  !9
        if (TOP .gt. 0) ARCTAN = 90.                 !3
        RETURN
      else
        ARCTAN = ATAN(TOP/BOTTOM)*57.2957795
        if (KEY .eq. 0) RETURN
        if (TOP .le. 0) then
          if (BOTTOM .lt. 0) ARCTAN = ARCTAN + 180.  !6
          if (BOTTOM .eq. 0) ARCTAN = 0.             !9
          if (BOTTOM .gt. 0) ARCTAN = ARCTAN + 360.  !8
        else 
          if (BOTTOM .lt. 0) ARCTAN = ARCTAN + 180.  !6
          if (BOTTOM .eq. 0) ARCTAN = 90.            !3
          RETURN
        endif
      endif
      
      END

!*********************************************************************
      FUNCTION DAYJUL(YR,XMONTH,DAY)
!
!*** THIS ROUTINE COMPUTES THE JULIAN DAY (AS A REAL VARIABLE)
!
!*********************************************************************
      REAL :: YR,XMONTH,DAY
      DIMENSION DAYT(12),DAYS(12)
      DATA DAYT/0.,31.,59.,90.,120.,151.,181.,212.,243.,273.,304.,334./
      DATA DAYS(1),DAYS(2) /0.,31./
      
      DINC=0.
      YRLP=MOD((YR-1900.),4.)
      IF(YRLP .EQ. 0.) DINC=1.
      DO I=3,12             !1   MEB change for Gnu fortran issue
        DAYS(I)=DAYT(I)+DINC
      ENDDO                 !1
      DAYJUL=DAYS(IFIX(XMONTH))+DAY
      
      END
      