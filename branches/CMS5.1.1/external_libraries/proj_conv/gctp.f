!!**********************************************************
!      program gctp
!!Uses NAD27 Datums      
!!***********************************************************      
!      IMPLICIT REAL*8 (A-Z)
!      integer :: INSYS,JNZONE,INUNIT,INSPH,IPR,JPR
!      integer :: IOSYS,IOZONE,IOUNIT,IOSPH,IFLG
!      DIMENSION CRDIN(2),CRDIO(2),TPARIN(15),TPARIO(15)
!      
!      !Input Coordinates
!      CRDIN(1)=208390.0
!      CRDIN(2)=190031.0
!      INSYS = 2 !State Plane
!      JNZONE = 4602 !Zone
!      INUNIT = 2 !meters
!      INSPH = 0 !Sphere 
!      IPR = 0 !Print error messages 0-yes, 1-no
!      JPR = 0 !Print parameters 0-yes, 1-no
!      
!      !Output Coordinates
!      IOSYS = 0
!      IOZONE = 0
!      IOUNIT = 4 !Degrees
!      IOSPH = 0
!      call GTPZ0(CRDIN,INSYS,JNZONE,TPARIN,INUNIT,INSPH,IPR,JPR,  
!     .                 CRDIO,IOSYS,IOZONE,TPARIO,IOUNIT,IOSPH,IFLG)    
!
!      write(*,*) CRDIO
!
!      stop      
!      endprogram

C     GCTP  GENERAL CARTOGRAPHIC COORDINATES TRANSFORMATION PACKAGE     0000000
C                   ADJLZ0                                              0000001
C                                                  FIPS CODES 6-10-82   0000002
C ** U.S.G.S. GENERAL MAP PROJECTION PACKAGE ...... DR. A. A. ELASSAL **0000003
C ** MODULE I                VERSION 1.0.0            NOVEMBER 1,1980 **0000004
C **********************************************************************0000005
      DOUBLE PRECISION FUNCTION ADJLZ0 (LON)                            0000006
C                                                                       0000007
C FUNCTION TO ADJUST LONGITUDE ANGLE TO MODULE 180 DEGREES.             0000008
C                                                                       0000009
      IMPLICIT REAL*8 (A-Z)                                             0000010
      DATA TWO,PI /2.0D0,3.14159265358979323846D0/                      0000011
C                                                                       0000012
  020 ADJLZ0 = LON                                                      0000013
      IF (DABS(LON) .LE. PI) RETURN                                     0000014
      TWOPI = TWO * PI                                                  0000015
      LON = LON - DSIGN (TWOPI,LON)                                     0000016
      GO TO 020                                                         0000017
C                                                                       0000018
      END                                                               0000019

C                   AL01Z0                                                     
C **********************************************************************0000021
C ** U.S.G.S. GENERAL MAP PROJECTION PACKAGE ...... DR. A. A. ELASSAL **0000022
C ** MODULE I                VERSION 1.0.0            NOVEMBER 1,1980 **0000023
C **********************************************************************0000024
      SUBROUTINE AL01Z0 (CRDIN,CRDOUT,FLAG)                             0000025
C                                                                       0000026
C SUBROUTINE TO COMPUTE TRANSFORMATION BETWEEN GEOGRAPHIC AND           0000027
C ALASKA STATE ZONE NO. 1.                                              0000028
C FLAG = 0, MEANS PLANE TO GEOGRAPHIC.                                  0000029
C FLAG = 1, MEANS GEOGRAPHIC TO PLANE.                                  0000030
C                                                                       0000031
      IMPLICIT REAL*8 (A-Z)                                             0000032
      INTEGER*4 FLAG                                                    0000033
      DIMENSION CRDIN(1),CRDOUT(1)                                      0000034
      DATA B,C,D /1.00029977273D0,0.00447599131D0,6386352.67013D0/      0000035
      DATA F,G,E /0.327015517176D0,0.945018968871D0,0.082271854223003D0/0000036
      DATA PI,EPS /3.141592653589793D0,2.718281828459045D0/             0000037
      DATA LO /1.771754086D0/                                           0000038
      DATA C1,C2 /0.182880365761D0,0.243840487681D0/                    0000039
      DATA C3,C4 /7000000.0D0,1000000.0D0/                              0000040
      DATA ONE,TWO,FOUR /1.0D0,2.0D0,4.0D0/                             0000041
      DATA MFEET,FALSE /3.280833333333D0,5000000.0D0/                   0000042
C                                                                       0000043
      IF (FLAG .EQ. 0) GO TO 020                                        0000044
C                                                                       0000045
C GEOGRAPHIC TO STATE PLANE TRANSFORMATION.                             0000046
C                                                                       0000047
      GEOG1 = DABS(CRDIN(1))                                            0000048
      GEOG2 = CRDIN(2)                                                  0000049
      ESINP = E * DSIN(GEOG2)                                           0000050
      TANP = DTAN((PI / FOUR) + (GEOG2 / TWO))                          0000051
      MU = DLOG(TANP) - (E / TWO) * DLOG((ONE + ESINP) / (ONE - ESINP)) 0000052
      CON = B * MU + C                                                  0000053
      CON1 = EPS**CON                                                   0000054
      CON =-CON                                                         0000055
      CON2 = EPS**CON                                                   0000056
      P = (CON1 - CON2) / TWO                                           0000057
      Q = (CON1 + CON2) / TWO                                           0000058
      CON = B * (GEOG1 - LO)                                            0000059
      CON1 = DSIN(CON)                                                  0000060
      CON2 = DCOS(CON)                                                  0000061
      U = D * DATAN((G * P + F * CON1) / CON2)                          0000062
      CON = F * P - G * CON1                                            0000063
      V = (D / TWO) * DLOG((Q + CON) / (Q - CON))                       0000064
      CRDOUT(1) = MFEET * (-0.6D0 * U + 0.8D0 * V + FALSE)              0000065
      CRDOUT(2) = MFEET * ( 0.8D0 * U + 0.6D0 * V - FALSE)              0000066
      RETURN                                                            0000067
C                                                                       0000068
C STATE PLANE TO GEOGRAPHIC TRANSFORMATION.                             0000069
C                                                                       0000070
  020 U =-C1 * CRDIN(1) + C2 * CRDIN(2) + C3                            0000071
      V = C2 * CRDIN(1) + C1 * CRDIN(2) - C4                            0000072
      CON = V / D                                                       0000073
      CON1 = EPS**CON                                                   0000074
      CON =-CON                                                         0000075
      CON2 = EPS**CON                                                   0000076
      R = (CON1 - CON2) / TWO                                           0000077
      S = (CON1 + CON2) / TWO                                           0000078
      K1 = DSIN(U / D)                                                  0000079
      K2 = DCOS(U / D)                                                  0000080
      CON = F * R + G * K1                                              0000081
      MU = (ONE / (TWO * B)) * DLOG((S + CON) / (S - CON)) - (C / B)    0000082
      KI = TWO * DATAN(EPS**MU) - PI / TWO                              0000083
      CON1 = DSIN(KI)                                                   0000084
      CON2 = DCOS(KI)                                                   0000085
      CRDOUT(2) = KI + (0.006761032571D0 + 0.000053172205D0 * CON2**2 + 0000086
     .            0.573027D-6 * CON2**4 + 0.7128D-8 * CON2**6) * CON1 * 0000087
     .            CON2                                                  0000088
      CRDOUT(1) =-LO - (ONE / B) * DATAN((F * K1 - G * R) / K2)         0000089
      RETURN                                                            0000090
C                                                                       0000091
      END                                                               0000092

C                   AL29Z0                                                     
C **********************************************************************0000094
C ** U.S.G.S. GENERAL MAP PROJECTION PACKAGE ...... DR. A. A. ELASSAL **0000095
C ** MODULE I                VERSION 1.0.0            NOVEMBER 1,1980 **0000096
C **********************************************************************0000097
      SUBROUTINE AL29Z0 (CRDIN,CRDOUT,ZONE,FLAG)                        0000098
C                                                                       0000099
C SUBROUTINE TO COMPUTE TRANSFORMATION BETWEEN GEOGRAPHIC AND           0000100
C ALASKA STATE ZONES NO. 2 THROUGH 9.                                   0000101
C FLAG = 0, MEANS PLANE TO GEOGRAPHIC.                                  0000102
C FLAG = 1, MEANS GEOGRAPHIC TO PLANE.                                  0000103
C                                                                       0000104
      IMPLICIT REAL*8 (A-Z)                                             0000105
      INTEGER*4 ZONE,FLAG,IND                                           0000106
      DIMENSION CRDIN(1),CRDOUT(1),CONST(2,8)                           0000107
      DATA CONST /500000.0D0 , 511200.0D0,                              0000108
     .            500000.0D0 , 525600.0D0,                              0000109
     .            500000.0D0 , 540000.0D0,                              0000110
     .            500000.0D0 , 554400.0D0,                              0000111
     .            500000.0D0 , 568800.0D0,                              0000112
     .            700000.0D0 , 583200.0D0,                              0000113
     .            500000.0D0 , 597600.0D0,                              0000114
     .            600000.0D0 , 612000.0D0/                              0000115
      DATA RADSEC /206264.806247D0/                                     0000116
      DATA ONE,TWO /1.0D0,2.0D0/                                        0000117
C                                                                       0000118
      IND = ZONE - 1                                                    0000119
      C = CONST(1,IND)                                                  0000120
      CM = CONST(2,IND)                                                 0000121
      IF (FLAG .EQ. 0) GO TO 020                                        0000122
C                                                                       0000123
C GEOGRAPHIC TO STATE PLANE TRANSFORMATION.                             0000124
C                                                                       0000125
      GEOG1 = DABS(CRDIN(1))                                            0000126
      GEOG2 = CRDIN(2)                                                  0000127
      C1 = DCOS(GEOG2)                                                  0000128
      C2 = C1 * C1                                                      0000129
      C3 = C2 * C2                                                      0000130
      C4 = C2 * C3                                                      0000131
      C5 = (CM - GEOG1 * RADSEC) * 1.0D-4                               0000132
      C6 = C5 * C5                                                      0000133
      C7 = C6 * C6                                                      0000134
      C8 = DSQRT(ONE + 0.0068147849 * C2)                               0000135
      C9 = DSQRT(ONE - C2) * C1                                         0000136
      CRDOUT(1) = C + (1017862.15D0 * C1 / C8) * C5 * (ONE -            0000137
     .            3.91740509D-4 * C6 * (ONE - TWO * C2 - 0.681478D-2 *  0000138
     .            C3) + 4.60382D-8 * C7 * (ONE - 20.0D0 * C2 +          0000139
     .            23.6047D0 * C3 + 0.4907D0 * C4))                      0000140
      CRDOUT(2) = 101.269278503D0 * (GEOG2 * RADSEC - 193900.05442D0    0000141
     .            - (1052.893943D0 - 4.483386D0 * C2 + 2.3559D-2 * C3) *0000142
     .            C9) + (24673.6748D0 * C9 * C6 / C8) * (ONE +          0000143
     .            1.958703D-4 * C6 * (-ONE + 6.0D0 * C2 + 6.133306D-2 * 0000144
     .            C3 + 1.8577D-4 * C4) + 1.5346D-8 * C7 * (ONE -        0000145
     .            60.0D0 * C2 + 117.75D0 * C3 + 4.089D0 * C4))          0000146
      RETURN                                                            0000147
C                                                                       0000148
C STATE PLANE TO GEOGRAPHIC TRANSFORMATION.                             0000149
C                                                                       0000150
  020 OMEGA = 193900.05442D0 + 0.00987466302498D0 * CRDIN(2)            0000151
      C1 = DCOS(OMEGA / RADSEC)                                         0000152
      C2 = C1 * C1                                                      0000153
      C3 = C2 * C2                                                      0000154
      PHI = OMEGA + (1047.546691D0 + 6.193011 * C2 + 5.0699D-2 * C3) *  0000155
     .      DSQRT(ONE - C2) * C1                                        0000156
      C1 = DCOS(PHI / RADSEC)                                           0000157
      C2 = C1 * C1                                                      0000158
      C3 = C2 * C2                                                      0000159
      C4 = (CRDIN(1) - C) * 1.0D-6                                      0000160
      C5 = C4 * C4                                                      0000161
      C6 = C5 * C5                                                      0000162
      C7 = ONE + 0.0068147849 * C2                                      0000163
      C8 = C7 * C7                                                      0000164
      CRDOUT(2) = (PHI - 233.973645D0 * C5 * C8 * DSQRT((ONE / C2) -    0000165
     .            ONE) * (ONE - 1.8905604D-4 * C5 * (1.9591113D0 +      0000166
     .            (3.0D0 / C2) + 8.1359D-2 * C2 + 2.79D-4 * C3) +       0000167
     .            1.42969D-8 * C6 * C7 * (15.5D0 + (45.0D0 / C3) -      0000168
     .            (0.307D0 / C2) + 1.53D0 * C2))) / RADSEC              0000169
      CRDOUT(1) =-(CM - 9824.513072D0 * DSQRT(C7) * C4 / C1 *           0000170
     .            (ONE - 3.7811208D-4 * C7 * C5 * (-TWO + C7 + (TWO /   0000171
     .            C2)) + 4.2890624D-8 * C7 * C6 * (1.054 + (24.0D0 /    0000172
     .            C3) - (20.0D0 / C2) - 1.36D-2 * C2))) / RADSEC        0000173
      RETURN                                                            0000174
C                                                                       0000175
      END                                                               0000176

C                   BLOCKD                                                     
C **********************************************************************0000178
C ** U.S.G.S. GENERAL MAP PROJECTION PACKAGE ...... DR. A. A. ELASSAL **0000179
C ** MODULE I                VERSION 1.0.0            NOVEMBER 1,1980 **0000180
C **********************************************************************0000181
      BLOCK DATA                                                        0000182
C                                                                       0000183
C INITIALIZATION OF ELLIPSOID TO CLARK'S 1866 PARAMETERS.               0000184
C                                                                       0000185
      IMPLICIT REAL*8 (A-Z)                                             0000186
      INTEGER*4 IPEMSG,IPPARM                                           0000187
C                                                                       0000188
      COMMON /ELLPZ0/ AZ,EZ,ESZ,E0Z,E1Z,E2Z,E3Z                         0000189
      COMMON /SPHRZ0/ AZZ                                               0000190
      COMMON /PRINZ0/ IPEMSG,IPPARM                                     0000191
C                                                                       0000192
      DATA AZ  /0.6378206400000000D+07/                                 0000193
      DATA EZ  /0.8227185422300323D-01/                                 0000194
      DATA ESZ /0.6768657997291094D-02/                                 0000195
      DATA E0Z /0.9983056818784341D+00/                                 0000196
      DATA E1Z /0.2542555507651308D-02/                                 0000197
      DATA E2Z /0.2698084527466011D-05/                                 0000198
      DATA E3Z /0.1003393903560134D+01/                                 0000199
C                                                                       0000200
      DATA AZZ /0.6370997000000000D+07/                                 0000201
C                                                                       0000202
      DATA IPEMSG /0/                                                   0000203
      DATA IPPARM /0/                                                   0000204
C                                                                       0000205
      END                                                               0000206

C                   DMSPZ0                                                     
C **********************************************************************0000208
C ** U.S.G.S. GENERAL MAP PROJECTION PACKAGE ...... DR. A. A. ELASSAL **0000209
C ** MODULE I                VERSION 1.0.0            NOVEMBER 1,1980 **0000210
C **********************************************************************0000211
      DOUBLE PRECISION FUNCTION DMSPZ0 (DMS)                            0000212
C                                                                       0000213
C SUBROUTINE TO CONVERT UNPACKED DMS TO PACKED DMS ANGLE                0000214
C                                                                       0000215
      IMPLICIT REAL*8 (A-H,O-Z)                                         0000216
      REAL*4 REL                                                        0000217
      INTEGER*4 DMS(1)                                                  0000218
      EQUIVALENCE (REL , INT)                                           0000219
      DATA CON1,CON2 /1000000.0D0,1000.0D0/                             0000220
      DATA NEG /'-'/                                                    0000221
C                                                                       0000222
      INT = DMS(4)                                                      0000223
      CON = DFLOAT (DMS(2)) * CON1 + DFLOAT (DMS(3)) * CON2 + REL       0000224
      IF (DMS(1) .EQ. NEG) CON = - CON                                  0000225
      DMSPZ0 = CON                                                      0000226
      RETURN                                                            0000227
C                                                                       0000228
      END                                                               0000229

C                   E0FNZ0                                                     
C **********************************************************************0000231
C ** U.S.G.S. GENERAL MAP PROJECTION PACKAGE ...... DR. A. A. ELASSAL **0000232
C ** MODULE I                VERSION 1.0.0            NOVEMBER 1,1980 **0000233
C **********************************************************************0000234
      DOUBLE PRECISION FUNCTION E0FNZ0 (ECCNTS)                         0000235
C                                                                       0000236
C FUNCTION TO COMPUTE CONSTANT (E0).                                    0000237
C                                                                       0000238
      IMPLICIT REAL*8 (A-Z)                                             0000239
      DATA QUART,ONE,ONEQ,THREE,SIXT /0.25D0,1.0D0,1.25D0,3.0D0,16.0D0/ 0000240
C                                                                       0000241
      E0FNZ0 = ONE - QUART * ECCNTS * (ONE + ECCNTS / SIXT *            0000242
     .         (THREE + ONEQ * ECCNTS))                                 0000243
C                                                                       0000244
      RETURN                                                            0000245
      END                                                               0000246

C                   E1FNZ0                                                     
C **********************************************************************0000248
C ** U.S.G.S. GENERAL MAP PROJECTION PACKAGE ...... DR. A. A. ELASSAL **0000249
C ** MODULE I                VERSION 1.0.0            NOVEMBER 1,1980 **0000250
C **********************************************************************0000251
      DOUBLE PRECISION FUNCTION E1FNZ0 (ECCNTS)                         0000252
C                                                                       0000253
C FUNCTION TO COMPUTE CONSTANT (E1).                                    0000254
C                                                                       0000255
      IMPLICIT REAL*8 (A-Z)                                             0000256
      DATA CON1,CON2,CON3 /0.375D0,0.25D0,0.46875D0/                    0000257
      DATA ONE /1.0D0/                                                  0000258
C                                                                       0000259
      E1FNZ0 = CON1 * ECCNTS * (ONE + CON2 * ECCNTS *                   0000260
     .         (ONE + CON3 * ECCNTS))                                   0000261
C                                                                       0000262
      RETURN                                                            0000263
      END                                                               0000264

C                   E2FNZ0                                                     
C **********************************************************************0000266
C ** U.S.G.S. GENERAL MAP PROJECTION PACKAGE ...... DR. A. A. ELASSAL **0000267
C ** MODULE I                VERSION 1.0.0            NOVEMBER 1,1980 **0000268
C **********************************************************************0000269
      DOUBLE PRECISION FUNCTION E2FNZ0 (ECCNTS)                         0000270
C                                                                       0000271
C FUNCTION TO COMPUTE CONSTANT (E2).                                    0000272
C                                                                       0000273
      IMPLICIT REAL*8 (A-Z)                                             0000274
      DATA CON1,CON2 /0.05859375D0,0.75D0/                              0000275
      DATA ONE /1.0D0/                                                  0000276
C                                                                       0000277
      E2FNZ0 = CON1 * ECCNTS * ECCNTS * (ONE + CON2 * ECCNTS)           0000278
C                                                                       0000279
      RETURN                                                            0000280
      END                                                               0000281

C                   E3FNZ0                                                     
C **********************************************************************0000283
C ** U.S.G.S. GENERAL MAP PROJECTION PACKAGE ...... DR. A. A. ELASSAL **0000284
C ** MODULE I                VERSION 1.0.0            NOVEMBER 1,1980 **0000285
C **********************************************************************0000286
      DOUBLE PRECISION FUNCTION E3FNZ0 (ECCENT)                         0000287
C                                                                       0000288
C FUNCTION TO COMPUTE CONSTANT (E3).                                    0000289
C                                                                       0000290
      IMPLICIT REAL*8 (A-Z)                                             0000291
      DATA ONE /1.0D0/                                                  0000292
C                                                                       0000293
      CON = ONE + ECCENT                                                0000294
      COM = ONE - ECCENT                                                0000295
      E3FNZ0 = DSQRT ((CON ** CON) * (COM ** COM))                      0000296
C                                                                       0000297
      RETURN                                                            0000298
      END                                                               0000299

C ****                                                             *****
C ** U.S.G.S. GENERAL MAP PROJECTION PACKAGE ...... JOHN F. WAANANEN  **
C ** MODULE I                VERSION 1.0.0              APRIL 6,1981  **
C ****                                                             *****
      SUBROUTINE GTPZ0(CRDIN,INSYS,JNZONE,TPARIN,INUNIT,INSPH,IPR,JPR,  
     .                 CRDIO,IOSYS,IOZONE,TPARIO,IOUNIT,IOSPH,IFLG)     
C                                                                       
C GENERAL PROGRAM FOR TRANSFORMATION BETWEEN VARIOUS REFERENCE SYSTEMS  
C     MODIFIED VERSION OF GTRNZ0 BY J.F. WAANANEN                       
C     SUBROUTINE GTPZ0 IS REQUIRED FOR PROGRAMS NO. L176 AND NO. L177   
C                                                                       
C INPUT ****                                                       **** 
C CRDIN  : COORDINATES IN INPUT SYSTEM (2 DP WORDS ARRAY).              
C INSYS  : CODE NUMBER OF INPUT COORDINATE SYSTEM (INTEGER).            
C            =  0 , GEOGRAPHIC                                          
C            =  1 , U T M                                               
C            =  2 , STATE PLANE                                         
C            =  3 , ALBERS CONICAL EQUAL-AREA                           
C            =  4 , LAMBERT CONFORMAL CONIC                             
C            =  5 , MERCATOR                                            
C            =  6 , POLAR STEREOGRAPHIC                                 
C            =  7 , POLYCONIC                                           
C            =  8 , EQUIDISTANT CONIC                                   
C            =  9 , TRANSVERSE MERCATOR                                 
C            = 10 , STEREOGRAPHIC                                       
C            = 11 , LAMBERT AZIMUTHAL EQUAL-AREA                        
C            = 12 , AZIMUTHAL EQUIDISTANT                              
C            = 13 , GNOMONIC                                          
C            = 14 , ORTHOGRAPHIC                                     
C            = 15 , GENERAL VERTICAL NEAR-SIDE PERSPECTIVE          
C            = 16 , SINUSOIDAL                                     
C            = 17 , EQUIRECTANGULAR (PLATE CARREE)                
C            = 18 , MILLER CYLINDRICAL                           
C            = 19 , VAN DER GRINTEN I                           
C            = 20 , OBLIQUE MERCATOR (HOTINE)                  
C            = 21 , SPACE OBLIQUE MERCATOR
C
C INZONE : CODE NUMBER OF INPUT COORDINATE ZONE (INTEGER).    
C TPARIN : PARAMETERS OF INPUT REFERENCE SYSTEM (15 DP WORDS ARRAY).    
C INUNIT : CODE NUMBER OF UNITS OF MEASURE FOR INPUT COORDS (INTEGER). 
C            = 0 , RADIANS.                                           
C            = 1 , FEET.                                             
C            = 2 , METERS.                                          
C            = 3 , SECONDS OF ARC.                                 
C            = 4 , DEGREES OF ARC.                                
C INSPH  : INPUT SPHEROID CODE.  SEE SPHDZ0 FOR PROPER CODES.    
C IPR    : PRINTOUT FLAG FOR ERROR MESSAGES. 0=YES, 1=NO        
C JPR    : PRINTOUT FLAG FOR PROJECTION PARAMETERS 0=YES, 1=NO 
C OUTPUT ***                                                  
C IOSYS  : CODE NUMBER OF OUTPUT COORDINATE SYSTEM (INTEGER).
C IOZONE : CODE NUMBER OF OUTPUT COORDINATE ZONE (INTEGER). 
C TPARIO : PARAMETERS OF OUTPUT REFERENCE SYSTEM (15 DP WORDS ARRAY).   
C IOUNIT : CODE NUMBER OF UNITS OF MEASURE FOR OUTPUT COORDS (INTEGER).
C IOSPH  : OUTPUT SPHEROID CODE.  SEE SPHDZ0 FOR PROPER CODES.        
C CRDIO  : COORDINATES IN OUTPUT REFERENCE SYSTEM (2 DP WORDS ARRAY).
C IFLG   : RETURN FLAG (INTEGER).                                   
C            =    0 , SUCCESSFUL TRANSFORMATION.                      
C            = 1001 , ILLEGAL INPUT SYSTEM CODE.                     
C            = 1002 , ILLEGAL OUTPUT SYSTEM CODE.                   
C            = 1003 , ILLEGAL INPUT UNIT CODE.                     
C            = 1004 , ILLEGAL OUTPUT UNIT CODE.                   
C            = 1005 , INCONSISTANT UNIT AND SYSTEM CODES FOR INPUT.        
C            = 1006 , INCONSISTANT UNIT AND SYSTEM CODES FOR OUTPUT.      
C            = 1007 , ILLEGAL INPUT ZONE CODE.                           
C            = 1008 , ILLEGAL OUTPUT ZONE CODE.                         
C      OTHERWISE , ERROR CODE FROM PROJECTION COMPUTATIONAL MODULE. 
C                                                                  
      IMPLICIT REAL*8 (A-H,O-Z)                                   
      INTEGER*4 SYSUNT(22)                                       
      DIMENSION CRDIN(2),CRDIO(2),TPARIN(15),TPARIO(15),COORD(2)
      DIMENSION DUMMY(2)                                       
      COMMON/ERRMZ0/ IERROR                                   
      COMMON/PRINZ0/ IP1,IP2                                 
      COMMON/ELLPZ0/ AZ,EZ,ESZ,E0Z,E1Z,E2Z,E3Z              
      COMMON/PROJZ0/IPROJ                                  
C                                                         
      DATA SYSUNT / 0 , 2 , 1 , 19*2 /                   
      DATA MAXUNT,MAXSYS / 4 , 21 /                     
      DATA JFLAG/0/                                    
C                                                     
C     SETUP                                          
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
C                                        
C CHECK VALIDITY OF CODES FOR UNITS OF MEASURE AND REFERENCE SYSTEMS.   
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
C                                                                  
C CHECK CONSISTANCY BETEEN UNITS OF MEASURE AND REFERENCE SYSTEM. 
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
C                                                                  
C COMPUTE TRANSFORMED COORDINATES AND ADJUST THEIR UNITS.         
  140 IF (INSYS .EQ. 0) GO TO 520                                
      IF (INZONE.GT.60 .OR. INSYS.EQ.1) GO TO 200               
      IFLG = 1007                                                 
      RETURN                                                  
C                                                            
C INVERSE TRANSFORMATION.                                   
  200 IPROJ=INSYS                                          
      IF (INSYS.GT.2) CALL SPHDZ0(INSPH,TPARIN)           
      GO TO (210,220,230,240,250,260,270,280,290,300,    
     .       310,320,330,340,350,360,370,380,390,400,410),INSYS
  210 IF (INZONE.EQ.0.AND.TPARIN(1).NE.0.0D0) GO TO 211                 
      TPARIN(1) = 1.0D6*DFLOAT(6*INZONE-183)                           
      TPARIN(2) = DSIGN(4.0D7,DFLOAT(JNZONE))                         
c  211 CALL SPHDZ0(INSPH,DUMMY)                                       
c      TPARIN(14) = DUMMY(1)                                         
c      TPARIN(15) = DUMMY(2)                                        
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
C                                                               
C FORWARD TRANSFORMATION.                                      
  540 IPROJ=IOSYS                                             
      IF (IOSYS.GT.2) CALL SPHDZ0(IOSPH,TPARIO)              
      GO TO (610,620,630,640,650,660,670,680,690,700,       
     .       710,720,730,740,750,760,770,780,790,800,810),IOSYS
  610 TPARIO(1) = COORD(1)                                              
      TPARIO(2) = COORD(2)                                             
c      CALL SPHDZ0(IOSPH,DUMMY)                                        
c      TPARIO(14) = DUMMY(1)                                          
c      TPARIO(15) = DUMMY(2)                                         
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

C                   GTRNZ0                                                     
C **********************************************************************
C ** U.S.G.S. GENERAL MAP PROJECTION PACKAGE ...... DR. A. A. ELASSAL **
C ** MODULE I                VERSION 1.0.0            NOVEMBER 1,1980 **
C **********************************************************************
      SUBROUTINE GTRNZ0 (CRDIN,INSYS,INZONE,TPARIN,INUNIT,
     .                   CRDIO,IOSYS,IOZONE,TPARIO,IOUNIT,IFLG)
C
C GENERAL PROGRAM FOR TRANSFORMATION BETWEEN VARIOUS REFERENCE SYSTEMS
C
C INPUT ***************************************************************
C CRDIN  : COORDINATES IN INPUT SYSTEM (2 DP WORDS ARRAY).
C INSYS  : CODE NUMBER OF INPUT COORDINATE SYSTEM (INTEGER).
C            =  0 , GEOGRAPHIC
C            =  1 , U T M                                               
C            =  2 , STATE PLANE                                         
C            =  3 , ALBERS CONICAL EQUAL-AREA                           
C            =  4 , LAMBERT CONFORMAL CONIC                             
C            =  5 , MERCATOR                                            
C            =  6 , POLAR STEREOGRAPHIC                                 
C            =  7 , POLYCONIC                                           
C            =  8 , EQUIDISTANT CONIC                                   
C            =  9 , TRANSVERSE MERCATOR                                 
C            = 10 , STEREOGRAPHIC                                       
C            = 11 , LAMBERT AZIMUTHAL EQUAL-AREA                        
C            = 12 , AZIMUTHAL EQUIDISTANT                               
C            = 13 , GNOMONIC                                            
C            = 14 , ORTHOGRAPHIC                                        
C            = 15 , GENERAL VERTICAL NEAR-SIDE PERSPECTIVE              
C            = 16 , SINUSOIDAL                                          
C            = 17 , EQUIRECTANGULAR (PLATE CARREE)                      
C            = 18 , MILLER CYLINDRICAL                                  
C            = 19 , VAN DER GRINTEN I                                   
C            = 20 , OBLIQUE MERCATOR (HOTINE)                           
C            = 21 , SPACE OBLIQUE MERCATOR
C
C INZONE : CODE NUMBER OF INPUT COORDINATE ZONE (INTEGER).              
C TPARIN : PARAMETERS OF INPUT REFERENCE SYSTEM (15 DP WORDS ARRAY).    
C INUNIT : CODE NUMBER OF UNITS OF MEASURE FOR INPUT COORDS (INTEGER).  
C            = 0 , RADIANS.                                             
C            = 1 , FEET.                                                
C            = 2 , METERS.                                              
C            = 3 , SECONDS OF ARC.                                      
C            = 4 , DEGREES OF ARC.                                      
C IOSYS  : CODE NUMBER OF OUTPUT COORDINATE SYSTEM (INTEGER).           
C IOZONE : CODE NUMBER OF OUTPUT COORDINATE ZONE (INTEGER).             
C TPARIO : PARAMETERS OF OUTPUT REFERENCE SYSTEM (15 DP WORDS ARRAY).   
C IOUNIT : CODE NUMBER OF UNITS OF MEASURE FOR OUTPUT COORDS (INTEGER). 
C                                                                       
C OUTPUT ************************************************************** 
C CRDIO  : COORDINATES IN OUTPUT REFERENCE SYSTEM (2 DP WORDS ARRAY).   
C IFLG   : RETURN FLAG (INTEGER).                                       
C            = 0 , SUCCESSFUL TRANSFORMATION.                           
C            = 1 , ILLEGAL INPUT SYSTEM CODE.                          
C            = 2 , ILLEGAL OUTPUT SYSTEM CODE.                         
C            = 3 , ILLEGAL INPUT UNIT CODE.                           
C            = 4 , ILLEGAL OUTPUT UNIT CODE.                         
C            = 5 , INCONSISTANT UNIT AND SYSTEM CODES FOR INPUT.    
C            = 6 , INCONSISTANT UNIT AND SYSTEM CODES FOR OUTPUT.  
C            = 7 , ILLEGAL INPUT ZONE CODE.                       
C            = 8 , ILLEGAL OUTPUT ZONE CODE.                     
C      OTHERWISE , ERROR CODE FROM PROJECTION COMPUTATIONAL MODULE.
C                                                                 
      IMPLICIT REAL*8 (A-H,O-Z)                                  
      INTEGER*4 SYSUNT(22)                                      
      COMMON /ERRMZ0/ IERROR                                   
      DIMENSION CRDIN(1),CRDIO(1),TPARIN(1),TPARIO(1),COORD(2)
      DATA SYSUNT / 0 , 2 , 1 , 19*2 /                       
      DATA MAXUNT,MAXSYS / 4 , 21 /                         
C                                                          
C CHECK VALIDITY OF CODES FOR UNITS OF MEASURE AND REFERENCE SYSTEMS.   
C                                                                      
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
C                                                         
C CHECK CONSISTANCY BETEEN UNITS OF MEASURE AND REFERENCE SYSTEM.       
C                                                                      
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
C                                                                   
C COMPUTE TRANSFORMED COORDINATES AND ADJUST THEIR UNITS.          
C                                                                 
  140 IF (INSYS .EQ. 0) GO TO 520                                
      IF (INZONE.GT.60 .OR. INSYS.EQ.1) GO TO 200               
      IFLG = 7                                                 
      RETURN                                                  
C                                                            
C INVERSE TRANSFORMATION.                                   
C                                                          
  200 GO TO (210,220,230,240,250,260,270,280,290,300,     
     .       310,320,330,340,350,360,370,380,390,400,410),INSYS
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
C                                               
C FORWARD TRANSFORMATION.                      
C                                             
  540 GO TO (610,620,630,640,650,660,670,680,690,700,                   
     .       710,720,730,740,750,760,770,780,790,800,810),IOSYS          
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
C                                                               
      END                                                      

C                   MLFNZ0                                                     
C **********************************************************************0000591
C ** U.S.G.S. GENERAL MAP PROJECTION PACKAGE ...... DR. A. A. ELASSAL **0000592
C ** MODULE I                VERSION 1.0.0            NOVEMBER 1,1980 **0000593
C **********************************************************************0000594
      DOUBLE PRECISION FUNCTION MLFNZ0 (E0,E1,E2,PHI)                   0000595
C                                                                       0000596
C FUNCTION TO COMPUTE CONSTANT (M).                                     0000597
C                                                                       0000598
      IMPLICIT REAL*8 (A-Z)                                             0000599
      DATA TWO,FOUR /2.0D0,4.0D0/                                       0000600
C                                                                       0000601
      MLFNZ0 = E0 * PHI - E1 * DSIN (TWO * PHI) + E2 * DSIN (FOUR * PHI)0000602
C                                                                       0000603
      RETURN                                                            0000604
      END                                                               0000605

C                   MSFNZ0                                                     
C **********************************************************************0000607
C ** U.S.G.S. GENERAL MAP PROJECTION PACKAGE ...... DR. A. A. ELASSAL **0000608
C ** MODULE I                VERSION 1.0.0            NOVEMBER 1,1980 **0000609
C **********************************************************************0000610
      DOUBLE PRECISION FUNCTION MSFNZ0 (ECCENT,SINPHI,COSPHI)           0000611
C                                                                       0000612
C FUNCTION TO COMPUTE CONSTANT (SMALL M).                               0000613
C                                                                       0000614
      IMPLICIT REAL*8 (A-Z)                                             0000615
      DATA ONE /1.0D0/                                                  0000616
C                                                                       0000617
      CON = ECCENT * SINPHI                                             0000618
      MSFNZ0 = COSPHI / DSQRT (ONE - CON * CON)                         0000619
C                                                                       0000620
      RETURN                                                            0000621
      END                                                               0000622

C                   PAKDZ0                                                     
C **********************************************************************0000624
C ** U.S.G.S. GENERAL MAP PROJECTION PACKAGE ...... DR. A. A. ELASSAL **0000625
C ** MODULE I                VERSION 1.0.0            NOVEMBER 1,1980 **0000626
C **********************************************************************0000627
      SUBROUTINE PAKDZ0 (PAK,DMS)                                       0000628
C                                                                       0000629
C SUBROUTINE TO CONVERT PACKED DMS TO UNPACKED DMS ANGLE.               0000630
C                                                                       0000631
      IMPLICIT REAL*8 (A-H,O-Z)                                         0000632
      REAL*4 REL                                                        0000633
      INTEGER*4 DMS(1)                                                  0000634
      EQUIVALENCE (REL , INT)                                           0000635
      DATA ZERO,CON1,CON2 /0.0D0,1000000.0D0,1000.0D0/                  0000636
      DATA IBLNK,NEG /' ','-'/                                          0000637
C                                                                       0000638
      DMS(1) = IBLNK                                                    0000639
      IF (PAK .LT. ZERO) DMS(1) = NEG                                   0000640
      CON = DABS (PAK)                                                  0000641
      DMS(2) = CON / CON1                                               0000642
      CON = DMOD (CON , CON1)                                           0000643
      DMS(3) = CON / CON2                                               0000644
      REL = DMOD (CON , CON2)                                           0000645
      DMS(4) = INT                                                      0000646
      RETURN                                                            0000647
C                                                                       0000648
      END                                                               0000649

C                   PAKRZ0                                                     
C **********************************************************************0000651
C ** U.S.G.S. GENERAL MAP PROJECTION PACKAGE ...... DR. A. A. ELASSAL **0000652
C ** MODULE I                VERSION 1.0.0            NOVEMBER 1,1980 **0000653
C **********************************************************************0000654
      DOUBLE PRECISION FUNCTION PAKRZ0 (ANG)                            0000655
C                                                                       0000656
C FUNCTION TO CONVERT DMS PACKED ANGLE INTO RADIANS.                    0000657
C                                                                       0000658
      IMPLICIT REAL*8 (A-H,O-Z)                                         0000659
      DATA SECRAD /0.4848136811095359D-5/                               0000660
C                                                                       0000661
C CONVERT ANGLE TO SECONDS OF ARC                                       0000662
C                                                                       0000663
      SEC = PAKSZ0 (ANG)                                                0000664
C                                                                       0000665
C CONVERT ANGLE TO RADIANS.                                             0000666
C                                                                       0000667
      PAKRZ0 = SEC * SECRAD                                             0000668
C                                                                       0000669
      RETURN                                                            0000670
      END                                                               0000671

C                   PAKSZ0                                                     
C **********************************************************************0000673
C ** U.S.G.S. GENERAL MAP PROJECTION PACKAGE ...... DR. A. A. ELASSAL **0000674
C ** MODULE I                VERSION 1.0.0            NOVEMBER 1,1980 **0000675
C **********************************************************************0000676
      DOUBLE PRECISION FUNCTION PAKSZ0 (ANG)                            0000677
C                                                                       0000678
C FUNCTION TO CONVERT DMS PACKED ANGLE INTO SECONDS OF ARC.             0000679
C                                                                       0000680
      IMPLICIT REAL*8 (A-H,M-Z)                                         0000681
      DIMENSION CODE(2)                                                 0000682
      DATA CODE /1000000.0D0,1000.0D0/                                  0000683
      DATA ZERO,ONE /0.0D0,1.0D0/                                       0000684
      DATA C1,C2 /3600.0D0,60.0D0/                                      0000685
C                                                                       0000686
C SEPERATE DEGREE FIELD.                                                0000687
C                                                                       0000688
      FACTOR = ONE                                                      0000689
      IF (ANG .LT. ZERO) FACTOR = - ONE                                 0000690
      SEC = DABS(ANG)                                                   0000691
      TMP = CODE(1)                                                     0000692
      I = SEC / TMP                                                     0000693
      IF (I .GT. 360) GO TO 020                                         0000694
      DEG = I                                                           0000695
C                                                                       0000696
C SEPERATE MINUTES FIELD.                                               0000697
C                                                                       0000698
      SEC = SEC - DEG * TMP                                             0000699
      TMP = CODE(2)                                                     0000700
      I = SEC / TMP                                                     0000701
      IF (I .GT. 60) GO TO 020                                          0000702
      MIN = I                                                           0000703
C                                                                       0000704
C SEPERATE SECONDS FIELD.                                               0000705
C                                                                       0000706
      SEC = SEC - MIN * TMP                                             0000707
      IF (SEC .GT. C2) GO TO 020                                        0000708
      SEC = FACTOR * (DEG * C1 + MIN * C2 + SEC)                        0000709
      GO TO 040                                                         0000710
C                                                                       0000711
C ERROR DETECTED IN DMS FORM.                                           0000712
C                                                                       0000713
  20  continue
c 020 PRINT 2000, ANG                                                   0000714
c2000 FORMAT (' ERROR PAKSZ0'/                                          0000715
c    .        ' ILLEGAL DMS FIELD =',F15.3)                             0000716
c     STOP 16                                                           0000717
C                                                                       0000718
      paksz0 = -999
      return
  040 PAKSZ0 = SEC                                                      0000719
C                                                                       0000720
      RETURN                                                            0000721
      END                                                               0000722

C                   PHI1Z0                                                     
C **********************************************************************0000724
C ** U.S.G.S. GENERAL MAP PROJECTION PACKAGE ...... DR. A. A. ELASSAL **0000725
C ** MODULE I                VERSION 1.0.0            NOVEMBER 1,1980 **0000726
C **********************************************************************0000727
      DOUBLE PRECISION FUNCTION PHI1Z0 (ECCENT,QS)                      0000728
C                                                                       0000729
C FUNCTION TO COMPUTE LATITUDE ANGLE (PHI-1).                           0000730
C                                                                       0000731
      IMPLICIT REAL*8 (A-Z)                                             0000732
      INTEGER*4 IERROR,IPEMSG,IPPARM                                    0000733
      INTEGER*4 II,NIT                                                  0000734
      COMMON /ERRMZ0/ IERROR                                            0000735
      COMMON /PRINZ0/ IPEMSG,IPPARM                                     0000736
      DATA HALF,ONE /0.5D0,1.0D0/                                       0000737
      DATA EPSLN,TOL,NIT /1.0D-7,1.0D-10,15/                            0000738
C                                                                       0000739
      PHI1Z0 =  DASIN (HALF * QS)                                       0000740
      IF (ECCENT .LT. EPSLN) RETURN                                     0000741
C                                                                       0000742
      ECCNTS = ECCENT * ECCENT                                          0000743
      PHI = PHI1Z0                                                      0000744
      DO 020 II = 1,NIT                                                 0000745
      SINPI = DSIN (PHI)                                                0000746
      COSPI = DCOS (PHI)                                                0000747
      CON = ECCENT * SINPI                                              0000748
      COM = ONE - CON * CON                                             0000749
      DPHI = HALF * COM * COM / COSPI * (QS / (ONE - ECCNTS) -          0000750
     .       SINPI / COM + HALF / ECCENT * DLOG ((ONE - CON) /          0000751
     .       (ONE + CON)))                                              0000752
      PHI = PHI + DPHI                                                  0000753
      IF (DABS(DPHI) .GT. TOL) GO TO 020                                0000754
      PHI1Z0 = PHI                                                      0000755
      RETURN                                                            0000756
  020 CONTINUE                                                          0000757
C                                                                       0000758
      IF (IPEMSG .EQ. 0) PRINT 2000, NIT,ECCENT,QS                      0000759
 2000 FORMAT (' ERROR PHI1Z0' /                                         0000760
     .        ' LATITUDE FAILED TO CONVERGE AFTER',I3,' ITERATIONS'/    0000761
     .        ' ECCENTRICITY =',D25.16,'   QS =',D25.16)                0000762
      IERROR = 001                                                      0000763
      RETURN                                                            0000764
C                                                                       0000765
      END                                                               0000766

C                   PHI2Z0                                                     
C **********************************************************************0000768
C ** U.S.G.S. GENERAL MAP PROJECTION PACKAGE ...... DR. A. A. ELASSAL **0000769
C ** MODULE I                VERSION 1.0.0            NOVEMBER 1,1980 **0000770
C **********************************************************************0000771
      DOUBLE PRECISION FUNCTION PHI2Z0 (ECCENT,TS)                      0000772
C                                                                       0000773
C FUNCTION TO COMPUTE LATITUDE ANGLE (PHI-2).                           0000774
C                                                                       0000775
      IMPLICIT REAL*8 (A-Z)                                             0000776
      INTEGER*4 IERROR,IPEMSG,IPPARM                                    0000777
      INTEGER*4 II,NIT                                                  0000778
      COMMON /ERRMZ0/ IERROR                                            0000779
      COMMON /PRINZ0/ IPEMSG,IPPARM                                     0000780
      DATA HALF,ONE,TWO /0.5D0,1.0D0,2.0D0/                             0000781
      DATA TOL,NIT /1.0D-10,15/                                         0000782
      DATA HALFPI /1.5707963267948966D0/                                0000783
C                                                                       0000784
      ECCNTH = HALF * ECCENT                                            0000785
      PHI = HALFPI - TWO * DATAN (TS)                                   0000786
      DO 020 II = 1,NIT                                                 0000787
      SINPI = DSIN (PHI)                                                0000788
      CON = ECCENT * SINPI                                              0000789
      DPHI = HALFPI - TWO * DATAN (TS * ((ONE - CON) /                  0000790
     .       (ONE + CON)) ** ECCNTH) - PHI                              0000791
      PHI = PHI + DPHI                                                  0000792
      IF (DABS(DPHI) .GT. TOL) GO TO 020                                0000793
      PHI2Z0 = PHI                                                      0000794
      RETURN                                                            0000795
  020 CONTINUE                                                          0000796
C                                                                       0000797
      IF (IPEMSG .EQ. 0) PRINT 2000, NIT,ECCENT,TS                      0000798
 2000 FORMAT (' ERROR PHI2Z0' /                                         0000799
     .        ' LATITUDE FAILED TO CONVERGE AFTER',I3,' ITERATIONS'/    0000800
     .        ' ECCENTRICITY =',D25.16,'   TS =',D25.16)                0000801
      IERROR = 002                                                      0000802
      RETURN                                                            0000803
C                                                                       0000804
      END                                                               0000805

C                   PHI3Z0                                                     
C **********************************************************************0000807
C ** U.S.G.S. GENERAL MAP PROJECTION PACKAGE ...... DR. A. A. ELASSAL **0000808
C ** MODULE I                VERSION 1.0.0            NOVEMBER 1,1980 **0000809
C **********************************************************************0000810
      DOUBLE PRECISION FUNCTION PHI3Z0 (ML,E0,E1,E2)                    0000811
C                                                                       0000812
C FUNCTION TO COMPUTE LATITUDE ANGLE (PHI-3).                           0000813
C                                                                       0000814
      IMPLICIT REAL*8 (A-Z)                                             0000815
      INTEGER*4 IERROR,IPEMSG,IPPARM                                    0000816
      INTEGER*4 II,NIT                                                  0000817
      COMMON /ERRMZ0/ IERROR                                            0000818
      COMMON /PRINZ0/ IPEMSG,IPPARM                                     0000819
      DATA TWO,FOUR /2.0D0,4.0D0/                                       0000820
      DATA TOL,NIT /1.0D-10,15/                                         0000821
C                                                                       0000822
      PHI = ML                                                          0000823
      DO 020 II = 1,NIT                                                 0000824
      DPHI = (ML + E1 * DSIN (TWO * PHI) - E2 * DSIN (FOUR * PHI)) /    0000825
     .       E0 - PHI                                                   0000826
      PHI = PHI + DPHI                                                  0000827
      IF (DABS(DPHI) .GT. TOL) GO TO 020                                0000828
      PHI3Z0 = PHI                                                      0000829
      RETURN                                                            0000830
  020 CONTINUE                                                          0000831
C                                                                       0000832
      IF (IPEMSG .EQ. 0) PRINT 2000, NIT,ML,E0,E1,E2                    0000833
 2000 FORMAT (' ERROR PHI3Z0' /                                         0000834
     .        ' LATITUDE FAILED TO CONVERGE AFTER',I3,' ITERATIONS'/    0000835
     .        ' ML =',D25.16,'   E0 =',D25.16/                          0000836
     .        ' E1 =',D25.16,'   E2 =',D25.16)                          0000837
      IERROR = 003                                                      0000838
      RETURN                                                            0000839
C                                                                       0000840
      END                                                               0000841

C                   PHI4Z0                                                     
C **********************************************************************0000843
C ** U.S.G.S. GENERAL MAP PROJECTION PACKAGE ...... DR. A. A. ELASSAL **0000844
C ** MODULE I                VERSION 1.0.0            NOVEMBER 1,1980 **0000845
C **********************************************************************0000846
      DOUBLE PRECISION FUNCTION PHI4Z0 (ECCNTS,E0,E1,E2,A,B,C)          0000847
C                                                                       0000848
C FUNCTION TO COMPUTE LATITUDE ANGLE (PHI-4).                           0000849
C                                                                       0000850
      IMPLICIT REAL*8 (A-Z)                                             0000851
      INTEGER*4 IERROR,IPEMSG,IPPARM                                    0000852
      INTEGER*4 II,NIT                                                  0000853
      COMMON /ERRMZ0/ IERROR                                            0000854
      COMMON /PRINZ0/ IPEMSG,IPPARM                                     0000855
      DATA ONE,TWO,FOUR /1.0D0,2.0D0,4.0D0/                             0000856
      DATA TOL,NIT /1.0D-10,15/                                         0000857
C                                                                       0000858
      PHI = A                                                           0000859
      DO 020 II = 1,NIT                                                 0000860
      SINPHI = DSIN (PHI)                                               0000861
      TANPHI = DTAN (PHI)                                               0000862
      C = TANPHI * DSQRT (ONE - ECCNTS * SINPHI * SINPHI)               0000863
      SIN2PH = DSIN (TWO * PHI)                                         0000864
      ML = E0 * PHI - E1 * SIN2PH + E2 * DSIN (FOUR * PHI)              0000865
      MLP = E0 - TWO * E1 * DCOS (TWO * PHI) + FOUR * E2 *              0000866
     .      DCOS (FOUR * PHI)                                           0000867
      CON1 = TWO * ML + C * (ML * ML + B) - TWO * A *                   0000868
     .       (C * ML + ONE)                                             0000869
      CON2 = ECCNTS * SIN2PH * (ML * ML + B - TWO * A * ML) / (TWO * C) 0000870
      CON3 = TWO * (A - ML) * (C * MLP - TWO / SIN2PH) - TWO * MLP      0000871
      DPHI = CON1 / (CON2 + CON3)                                       0000872
      PHI = PHI + DPHI                                                  0000873
      IF (DABS(DPHI) .GT. TOL) GO TO 020                                0000874
      PHI4Z0 = PHI                                                      0000875
      RETURN                                                            0000876
  020 CONTINUE                                                          0000877
C                                                                       0000878
      IF (IPEMSG .EQ. 0) PRINT 2000, NIT,E0,E1,E2,A,B,C,ECCNTS          0000879
 2000 FORMAT (' ERROR PHI4Z0' /                                         0000880
     .        ' LATITUDE FAILED TO CONVERGE AFTER',I3,' ITERATIONS'/    0000881
     .        ' E0 =',D25.16,'   E1 =',D25.16/                          0000882
     .        ' E2 =',D25.16,'   A  =',D25.16/                          0000883
     .        ' B  =',D25.16,'   C  =',D25.16/                          0000884
     .        ' ECCENTRICITY SQUARE =',D25.16)                          0000885
      IERROR = 004                                                      0000886
      RETURN                                                            0000887
C                                                                       0000888
      END                                                               0000889

C ****                                                             *****PJ1A   
C ** U.S.G.S. GENERAL MAP PROJECTION PACKAGE ...... DR. A. A. ELASSAL **PJ1A   
C ** MODULE I                VERSION 1.0.0            NOVEMBER 1,1980 **PJ1A   
C ****                                                             *****PJ1A   
C                              *  U T M  *                              PJ1A   
C ****                                                             *****PJ1A   
C                                                                       PJ1A   
      SUBROUTINE PJ1AZ0                                                 PJ1A   
C                                                                       PJ1A   
C     MODIFIED APRIL 1981  FROM PJ01Z0  BY JOHN F. WAANANEN             PJ1A  1
C     SUBROUTINE PJ1AZ0 IS REQUIRED FOR PROGRAMS NO. L176 AND NO. L177  PJ1A  1
C                                                                       PJ1A  1
      IMPLICIT REAL*8 (A-Z)                                             PJ1A  1
      INTEGER*4 IERROR,IPEMSG,IPPARM                                    PJ1A  1
      INTEGER*4 SWITCH,IND,ZONE,INFILE,ANG                              PJ1A  1
	integer*4 zzone
      DIMENSION BUFFL(15),DATA(1),GEOG(1),PROJ(1)                       PJ1A  1
      COMMON /ERRMZ0/ IERROR                                            PJ1A  1
      COMMON /PRINZ0/ IPEMSG,IPPARM                                     PJ1A  1
c     COMMON /WORKZ0/ BUFF(15),ANG(4)                                   PJ1A  1
	COMMON /WORKZ0/ BUFF(15)
	COMMON /WK1AZ0/ ANG(4)
      COMMON /ELLPZ0/ AZ,EZ,ESZ,E0Z,E1Z,E2Z,E3Z                         PJ1A  2
      real*4 rangs
      equivalence (rangs,ang(4))
      DATA HALFPI /1.5707963267948966D0/                                PJ1A  2
      DATA ZERO /0.0D0/                                                 PJ1A  2
      DATA SWITCH /0/                                                   PJ1A  2
C ....                                                             .....PJ1A  2
C       .  INITIALIZATION OF PROJECTION PARAMETERS (FILE INPUT)  .      PJ1A  2
C ....                                                             .....PJ1A  2
      ENTRY IF1AZ0 (INFILE)                                             PJ1A  2
      IERROR = 0                                                        PJ1A  2
      READ (INFILE,END=120) ZONE,BUFF                                   PJ1A  2
      IF (SWITCH.NE.0 .AND. SWITCH.EQ.ZONE) RETURN                      PJ1A  3
  020 IF (ZONE .EQ. ZERO) GO TO 040                                     PJ1A  3
      IF (ZONE.GE.1.AND.ZONE.LE.60) GO TO 100                           PJ1A  3
      IF (IPEMSG .EQ. 0) PRINT 2000,ZONE
 2000 FORMAT (' ILLEGAL ZONE NO : ',I10)
      IERROR = 011                                                      PJ1A  3
      RETURN                                                            PJ1A  3
  040 ZONE = PAKRZ0 (BUFF(1)) * 15.0D0 / HALFPI                         PJ1A  3
      IND = 1                                                           PJ1A  3
      IF (ZONE .LT. 0) IND = 0                                          PJ1A  4
      ZONE = MOD ((ZONE + 30) , 60) + IND                               PJ1A  4
      IF (ZONE) 060,080,100                                             PJ1A  4
  060 ZONE = ZONE + 60                                                  PJ1A  4
      GO TO 100                                                         PJ1A  4
  080 ZONE = 1                                                          PJ1A  4
  100 BUFFL(1) = BUFF(14)                                               PJ1A  4
      BUFFL(2) = BUFF(15)                                               PJ1A  4
      BUFFL(3) = 0.9996D0                                               PJ1A  4
      BUFFL(4) = ZERO                                                   PJ1A  4
      LON0 = DFLOAT (6 * ZONE - 183) * HALFPI / 90.0D0                  PJ1A  5
      CALL RADDZ0 (LON0,ANG)                                            PJ1A  5
      BUFFL(5) = DMSPZ0 (ANG)                                           PJ1A  5
      BUFFL(6) = ZERO                                                   PJ1A  5
      BUFFL(7) = 500000.0D0                                             PJ1A  5
      BUFFL(8) = ZERO                                                   PJ1A  5
      IF (BUFF(2) .LT. ZERO .OR. ZONE.LT.0) BUFFL(8) = 10000000.0D0     PJ1A  5
      IND = IPPARM                                                      PJ1A  5
      IPPARM = 1                                                        PJ1A  5
      SWITCH = 0                                                        PJ1A  5
      CALL IS09Z0 (ZONE,BUFFL)                                          PJ1A  6
      IPPARM = IND                                                      PJ1A  6
      IF (IERROR .NE. 0) GO TO 110                                      PJ1A  6
C                                                                       PJ1A  6
C LIST RESULTS OF PARAMETER INITIALIZATION.                             PJ1A  6
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
      SWITCH = ZONE                                                     PJ1A  7
  110 RETURN                                                            PJ1A  7
  120 IF (IPEMSG .EQ. 0) PRINT 2020
 2020 FORMAT (' MISSING PROJECTION PARAMETERS')
      IERROR = 012                                                      PJ1A  8
      RETURN                                                            PJ1A  8
C ....                                                             .....PJ1A  8
C      .  INITIALIZATION OF PROJECTION PARAMETERS (ENTRY INPUT)  .      PJ1A  8
C ....                                                             .....PJ1A  8
      ENTRY IS1AZ0 (ZZONE,DATA)                                         PJ1A  8
	zone = zzone
      IERROR = 0                                                        PJ1A  8
      IF (SWITCH.NE.0.AND.SWITCH.EQ.ZONE.AND.DATA(14).EQ.BUFF(1))RETURN PJ1A  8
      ZONE = IABS(ZONE)                                                 PJ1A  8
      SWITCH = 0                                                        PJ1A  8
      BUFF(1) = DATA(1)                                                 PJ1A  9
      BUFF(2) = DATA(2)                                                 PJ1A  9
      BUFF(14) = DATA(14)                                               PJ1A  9
      BUFF(15) = DATA(15)                                               PJ1A  9
      GO TO 020                                                         PJ1A  9
C ....                                                             .....PJ1A  9
C                      .  FORWARD TRANSFORMATION  .                     PJ1A  9
C ....                                                             .....PJ1A  9
      ENTRY PF1AZ0 (GEOG,PROJ)                                          PJ1A  9
      IERROR = 0                                                        PJ1A  9
      IF (SWITCH .NE. 0) GO TO 140                                      PJ1A 10
      IF (IPEMSG .EQ. 0) PRINT 2020
      IERROR = 013                                                      PJ1A 10
      RETURN                                                            PJ1A 10
  140 CALL PF09Z0 (GEOG,PROJ)                                           PJ1A 10
      RETURN                                                            PJ1A 10
C ....                                                             .....PJ1A 10
C                      .  INVERSE TRANSFORMATION  .                     PJ1A 10
C ....                                                             .....PJ1A 10
      ENTRY PI1AZ0 (PROJ,GEOG)                                          PJ1A 10
      IERROR = 0                                                        PJ1A 11
      IF (SWITCH .NE. 0) GO TO 160                                      PJ1A 11
      IF (IPEMSG .EQ. 0) PRINT 2020
      IERROR = 014                                                      PJ1A 11
      RETURN                                                            PJ1A 11
  160 CALL PI09Z0 (PROJ,GEOG)                                           PJ1A 11
      RETURN                                                            PJ1A 11
      END                                                               PJ1A 11

C                   PJ01Z0                                                     
C **********************************************************************0000891
C ** U.S.G.S. GENERAL MAP PROJECTION PACKAGE ...... DR. A. A. ELASSAL **0000892
C ** MODULE I                VERSION 1.0.0            NOVEMBER 1,1980 **0000893
C **********************************************************************0000894
C                              *  U T M  *                              0000895
C **********************************************************************0000896
C                                                                       0000897
      SUBROUTINE PJ01Z0                                                 0000898
C                                                                       0000899
      IMPLICIT REAL*8 (A-Z)                                             0000900
	integer*4 zzone
      INTEGER*4 IERROR,IPEMSG,IPPARM                                    0000901
      INTEGER*4 SWITCH,IND,ZONE,INFILE,ANG                              0000902
      COMMON /ERRMZ0/ IERROR                                            0000903
      COMMON /PRINZ0/ IPEMSG,IPPARM                                     0000904
c     COMMON /WORKZ0/ BUFF(15),BUFFL(15),ANG(4)                         0000905
	COMMON /WORKZ0/ BUFF(15)
	COMMON /WK01Z0/ BUFFL(15),ANG(4)
      DIMENSION DATA(1),GEOG(1),PROJ(1)                                 0000906
      DATA HALFPI /1.5707963267948966D0/                                0000907
      DATA ZERO /0.0D0/                                                 0000908
      DATA SWITCH /0/                                                   0000909
C                                                                       0000910
C ......................................................................0000911
C       .  INITIALIZATION OF PROJECTION PARAMETERS (FILE INPUT)  .      0000912
C ......................................................................0000913
C                                                                       0000914
      ENTRY IF01Z0 (INFILE,data)                                        0000915
C                                                                       0000916
      IERROR = 0                                                        0000917
      READ (INFILE,END=120) ZONE,BUFF                                   0000918
      IF (SWITCH.NE.0 .AND. SWITCH.EQ.ZONE) RETURN                      0000919
  020 IF (ZONE .EQ. ZERO) GO TO 040                                     0000920
      IF (ZONE.GE.1 .AND. ZONE.LE.60) GO TO 100                         0000921
      IF (IPEMSG .EQ. 0) PRINT 2000, ZONE                               0000922
 2000 FORMAT (' ERROR PJ01Z0'/                                          0000923
     .        ' ILLEGAL ZONE NO : ',I10)                                0000924
      IERROR = 011                                                      0000925
      RETURN                                                            0000926
  040 ZONE = PAKRZ0 (BUFF(1)) * 15.0D0 / HALFPI                         0000927
      IND = 1                                                           0000928
      IF (ZONE .LT. 0) IND = 0                                          0000929
      ZONE = MOD ((ZONE + 30) , 60) + IND                               0000930
      IF (ZONE) 060,080,100                                             0000931
  060 ZONE = ZONE + 60                                                  0000932
      GO TO 100                                                         0000933
  080 ZONE = 1                                                          0000934
  100 BUFFL(1) = ZERO                                                   0000935
      BUFFL(2) = ZERO                                                   0000936
      BUFFL(3) = 0.9996D0                                               0000937
      BUFFL(4) = ZERO                                                   0000938
      LON0 = DFLOAT (6 * ZONE - 183) * HALFPI / 90.0D0                  0000939
      CALL RADDZ0 (LON0,ANG)                                            0000940
      BUFFL(5) = DMSPZ0 (ANG)                                           0000941
      BUFFL(6) = ZERO                                                   0000942
      BUFFL(7) = 500000.0D0                                             0000943
      BUFFL(8) = ZERO                                                   0000944
      IF (BUFF(2) .LT. ZERO) BUFFL(8) = 10000000.0D0                    0000945
      IND = IPPARM                                                      0000946
      IPPARM = 1                                                        0000947
      SWITCH = 0                                                        0000948
      CALL IS09Z0 (ZONE,BUFFL)                                          0000949
      IPPARM = IND                                                      0000950
      IF (IERROR .NE. 0) GO TO 110                                      0000951
C                                                                       0000952
C LIST RESULTS OF PARAMETER INITIALIZATION.                             0000953
C                                                                       0000954
      IF (IPPARM .EQ. 0) PRINT 2010, ZONE                               0000955
 2010 FORMAT (' INITIALIZATION PARAMETERS (U T M PROJECTION)'/          0000956
     .        ' ZONE = ',I2)                                            0000957
      SWITCH = ZONE                                                     0000958
  110 RETURN                                                            0000959
  120 IF (IPEMSG .EQ. 0) PRINT 2020                                     0000960
 2020 FORMAT (' ERROR PJ01Z0'/                                          0000961
     .        ' MISSING PROJECTION PARAMETERS')                         0000962
      IERROR = 012                                                      0000963
      RETURN                                                            0000964
C                                                                       0000965
C ......................................................................0000966
C      .  INITIALIZATION OF PROJECTION PARAMETERS (ENTRY INPUT)  .      0000967
C ......................................................................0000968
C                                                                       0000969
      ENTRY IS01Z0 (ZZONE,DATA)                                         0000970
	zone = zzone
C                                                                       0000971
      IERROR = 0                                                        0000972
      IF (SWITCH.NE.0 .AND. SWITCH.EQ.ZONE) RETURN                      0000973
      BUFF(1) = DATA(1)                                                 0000974
      BUFF(2) = DATA(2)                                                 0000975
      GO TO 020                                                         0000976
C                                                                       0000977
C ......................................................................0000978
C                      .  FORWARD TRANSFORMATION  .                     0000979
C ......................................................................0000980
C                                                                       0000981
      ENTRY PF01Z0 (GEOG,PROJ)                                          0000982
C                                                                       0000983
      IERROR = 0                                                        0000984
      IF (SWITCH .NE. 0) GO TO 140                                      0000985
      IF (IPEMSG .EQ. 0) PRINT 2020                                     0000986
      IERROR = 013                                                      0000987
      RETURN                                                            0000988
  140 CALL PF09Z0 (GEOG,PROJ)                                           0000989
      RETURN                                                            0000990
C                                                                       0000991
C ......................................................................0000992
C                      .  INVERSE TRANSFORMATION  .                     0000993
C ......................................................................0000994
C                                                                       0000995
      ENTRY PI01Z0 (PROJ,GEOG)                                          0000996
C                                                                       0000997
      IERROR = 0                                                        0000998
      IF (SWITCH .NE. 0) GO TO 160                                      0000999
      IF (IPEMSG .EQ. 0) PRINT 2020                                     0001000
      IERROR = 014                                                      0001001
      RETURN                                                            0001002
  160 CALL PI09Z0 (PROJ,GEOG)                                           0001003
      RETURN                                                            0001004
C                                                                       0001005
      END                                                               0001006

C                   PJ02Z0                                                     
C **********************************************************************0001008
C ** U.S.G.S. GENERAL MAP PROJECTION PACKAGE ...... DR. A. A. ELASSAL **0001009
C ** MODULE I                VERSION 1.0.1               JUNE 10,1982 **0001010
C **********************************************************************0001011
C                           *  STATE PLANE  *                           0001012
C **********************************************************************0001013
C                                                                       0001014
      SUBROUTINE PJ02Z0                                                 0001015
C                       CODE NUMBERS MODIFIED 6-10-82 BY J.F.WAANANEN   0001016
      IMPLICIT REAL*8(A-H,O-Z)                                          0001017
	integer*4 zzone
      INTEGER*4 SWITCH,ZONE                                             0001018
      COMMON /ERRMZ0/ IERROR                                            0001019
      COMMON /PRINZ0/ IPEMSG,IPPARM                                     0001020
      COMMON /WORKZ0/ BUFF(15)                                          0001021
      DIMENSION GEOG(1),PROJ(1),DATA(1)                                 0001022
      DIMENSION ITEM(131),ID(9,131),TABLE(11,122)                       0001023
      DATA ITEM /0101,0102,5010, -1 ,0201,0202,0203,0301,0302,0401,0402,0001024
     .           0403,0404,0405,0406,0407,0501,0502,0503,0600,0700,0901,0001025
     .           0902,0903,1001,1002,5101,5102,5103,5104,5105,1101,1102,0001026
     .           1103,1201,1202,1301,1302,1401,1402,1501,1502,1601,1602,0001027
     .           1701,1702,1703,1801,1802,1900,2001,2002,2101,2102,2103,0001028
     .           2111,2112,2113,2201,2202,2203,2301,2302,2401,2402,2403,0001029
     .           2501,2502,2503,2601,2602,2701,2702,2703,2800,2900,3001,0001030
     .           3002,3003,3101,3102,3103,3104,3200,3301,3302,3401,3402,0001031
     .           3501,3502,3601,3602,3701,3702,3800,3901,3902,4001,4002,0001032
     .           4100,4201,4202,4203,4204,4205,4301,4302,4303,4400,4501,0001033
     .           4502,4601,4602,4701,4702,4801,4802,4803,4901,4902,4903,0001034
     .           4904,5001,5002,5003,5004,5005,5006,5007,5008,5009/     0001035
C ....................................................................  0001036
C ALABAMA                    EAST                                      T0001037
       DATA ID(1,  1),ID(2,  1),ID(3,  1) /4HALAB,4HAMA ,4H    /        0001038
       DATA ID(4,  1),ID(5,  1),ID(6,  1) /4H    ,4HEAST,4H    /        0001039
       DATA ID(7,  1),ID(8,  1),ID(9,  1) /4H    ,4H    ,0/             0001040
C ALABAMA                    WEST                                      T0001041
       DATA ID(1,  2),ID(2,  2),ID(3,  2) /4HALAB,4HAMA ,4H    /        0001042
       DATA ID(4,  2),ID(5,  2),ID(6,  2) /4H    ,4HWEST,4H    /        0001043
       DATA ID(7,  2),ID(8,  2),ID(9,  2) /4H    ,4H    ,0/             0001044
C ALASKA                     ZONE NO. 10                                0001045
       DATA ID(1,  3),ID(2,  3),ID(3,  3) /4HALAS,4HKA  ,4H    /        0001046
       DATA ID(4,  3),ID(5,  3),ID(6,  3) /4H    ,4HZONE,4H NO./        0001047
       DATA ID(7,  3),ID(8,  3),ID(9,  3) /4H 10 ,4H    ,1/             0001048
C ALASKA                                                               T0001049
       DATA ID(1,  4),ID(2,  4),ID(3,  4) /4HALAS,4HKA  ,4H    /        0001050
       DATA ID(4,  4),ID(5,  4),ID(6,  4) /4H    ,4H    ,4H    /        0001051
       DATA ID(7,  4),ID(8,  4),ID(9,  4) /4H    ,4H    ,0/             0001052
C ARIZONA                    EAST                                      T0001053
       DATA ID(1,  5),ID(2,  5),ID(3,  5) /4HARIZ,4HONA ,4H    /        0001054
       DATA ID(4,  5),ID(5,  5),ID(6,  5) /4H    ,4HEAST,4H    /        0001055
       DATA ID(7,  5),ID(8,  5),ID(9,  5) /4H    ,4H    ,0/             0001056
C ARIZONA                    CENTRAL                                   T0001057
       DATA ID(1,  6),ID(2,  6),ID(3,  6) /4HARIZ,4HONA ,4H    /        0001058
       DATA ID(4,  6),ID(5,  6),ID(6,  6) /4H    ,4HCENT,4HRAL /        0001059
       DATA ID(7,  6),ID(8,  6),ID(9,  6) /4H    ,4H    ,0/             0001060
C ARIZONA                    WEST                                      T0001061
       DATA ID(1,  7),ID(2,  7),ID(3,  7) /4HARIZ,4HONA ,4H    /        0001062
       DATA ID(4,  7),ID(5,  7),ID(6,  7) /4H    ,4HWEST,4H    /        0001063
       DATA ID(7,  7),ID(8,  7),ID(9,  7) /4H    ,4H    ,0/             0001064
C ARKANSAS                   NORTH                                     L0001065
       DATA ID(1,  8),ID(2,  8),ID(3,  8) /4HARKA,4HNSAS,4H    /        0001066
       DATA ID(4,  8),ID(5,  8),ID(6,  8) /4H    ,4HNORT,4HH   /        0001067
       DATA ID(7,  8),ID(8,  8),ID(9,  8) /4H    ,4H    ,1/             0001068
C ARKANSAS                   SOUTH                                     L0001069
       DATA ID(1,  9),ID(2,  9),ID(3,  9) /4HARKA,4HNSAS,4H    /        0001070
       DATA ID(4,  9),ID(5,  9),ID(6,  9) /4H    ,4HSOUT,4HH   /        0001071
       DATA ID(7,  9),ID(8,  9),ID(9,  9) /4H    ,4H    ,1/             0001072
C CALIFORNIA                 I                                         L0001073
       DATA ID(1, 10),ID(2, 10),ID(3, 10) /4HCALI,4HFORN,4HIA  /        0001074
       DATA ID(4, 10),ID(5, 10),ID(6, 10) /4H    ,4HI   ,4H    /        0001075
       DATA ID(7, 10),ID(8, 10),ID(9, 10) /4H    ,4H    ,1/             0001076
C CALIFORNIA                 II                                        L0001077
       DATA ID(1, 11),ID(2, 11),ID(3, 11) /4HCALI,4HFORN,4HIA  /        0001078
       DATA ID(4, 11),ID(5, 11),ID(6, 11) /4H    ,4HII  ,4H    /        0001079
       DATA ID(7, 11),ID(8, 11),ID(9, 11) /4H    ,4H    ,1/             0001080
C CALIFORNIA                 III                                       L0001081
       DATA ID(1, 12),ID(2, 12),ID(3, 12) /4HCALI,4HFORN,4HIA  /        0001082
       DATA ID(4, 12),ID(5, 12),ID(6, 12) /4H    ,4HIII ,4H    /        0001083
       DATA ID(7, 12),ID(8, 12),ID(9, 12) /4H    ,4H    ,1/             0001084
C CALIFORNIA                 IV                                        L0001085
       DATA ID(1, 13),ID(2, 13),ID(3, 13) /4HCALI,4HFORN,4HIA  /        0001086
       DATA ID(4, 13),ID(5, 13),ID(6, 13) /4H    ,4HIV  ,4H    /        0001087
       DATA ID(7, 13),ID(8, 13),ID(9, 13) /4H    ,4H    ,1/             0001088
C CALIFORNIA                 V                                         L0001089
       DATA ID(1, 14),ID(2, 14),ID(3, 14) /4HCALI,4HFORN,4HIA  /        0001090
       DATA ID(4, 14),ID(5, 14),ID(6, 14) /4H    ,4HV   ,4H    /        0001091
       DATA ID(7, 14),ID(8, 14),ID(9, 14) /4H    ,4H    ,1/             0001092
C CALIFORNIA                 VI                                        L0001093
       DATA ID(1, 15),ID(2, 15),ID(3, 15) /4HCALI,4HFORN,4HIA  /        0001094
       DATA ID(4, 15),ID(5, 15),ID(6, 15) /4H    ,4HVI  ,4H    /        0001095
       DATA ID(7, 15),ID(8, 15),ID(9, 15) /4H    ,4H    ,1/             0001096
C CALIFORNIA                 VII                                       L0001097
       DATA ID(1, 16),ID(2, 16),ID(3, 16) /4HCALI,4HFORN,4HIA  /        0001098
       DATA ID(4, 16),ID(5, 16),ID(6, 16) /4H    ,4HVII ,4H    /        0001099
       DATA ID(7, 16),ID(8, 16),ID(9, 16) /4H    ,4H    ,1/             0001100
C COLORADO                   NORTH                                     L0001101
       DATA ID(1, 17),ID(2, 17),ID(3, 17) /4HCOLO,4HRADO,4H    /        0001102
       DATA ID(4, 17),ID(5, 17),ID(6, 17) /4H    ,4HNORT,4HH   /        0001103
       DATA ID(7, 17),ID(8, 17),ID(9, 17) /4H    ,4H    ,1/             0001104
C COLORADO                   CENTRAL                                   L0001105
       DATA ID(1, 18),ID(2, 18),ID(3, 18) /4HCOLO,4HRADO,4H    /        0001106
       DATA ID(4, 18),ID(5, 18),ID(6, 18) /4H    ,4HCENT,4HRAL /        0001107
       DATA ID(7, 18),ID(8, 18),ID(9, 18) /4H    ,4H    ,1/             0001108
C COLORADO                   SOUTH                                     L0001109
       DATA ID(1, 19),ID(2, 19),ID(3, 19) /4HCOLO,4HRADO,4H    /        0001110
       DATA ID(4, 19),ID(5, 19),ID(6, 19) /4H    ,4HSOUT,4HH   /        0001111
       DATA ID(7, 19),ID(8, 19),ID(9, 19) /4H    ,4H    ,1/             0001112
C CONNECTICUT                ---                                       L0001113
       DATA ID(1, 20),ID(2, 20),ID(3, 20) /4HCONN,4HECTI,4HCUT /        0001114
       DATA ID(4, 20),ID(5, 20),ID(6, 20) /4H    ,4H--- ,4H    /        0001115
       DATA ID(7, 20),ID(8, 20),ID(9, 20) /4H    ,4H    ,1/             0001116
C DELAWARE                   ---                                       T0001117
       DATA ID(1, 21),ID(2, 21),ID(3, 21) /4HDELA,4HWARE,4H    /        0001118
       DATA ID(4, 21),ID(5, 21),ID(6, 21) /4H    ,4H--- ,4H    /        0001119
       DATA ID(7, 21),ID(8, 21),ID(9, 21) /4H    ,4H    ,0/             0001120
C FLORIDA                    EAST                                      T0001121
       DATA ID(1, 22),ID(2, 22),ID(3, 22) /4HFLOR,4HIDA ,4H    /        0001122
       DATA ID(4, 22),ID(5, 22),ID(6, 22) /4H    ,4HEAST,4H    /        0001123
       DATA ID(7, 22),ID(8, 22),ID(9, 22) /4H    ,4H    ,0/             0001124
C FLORIDA                    WEST                                      T0001125
       DATA ID(1, 23),ID(2, 23),ID(3, 23) /4HFLOR,4HIDA ,4H    /        0001126
       DATA ID(4, 23),ID(5, 23),ID(6, 23) /4H    ,4HWEST,4H    /        0001127
       DATA ID(7, 23),ID(8, 23),ID(9, 23) /4H    ,4H    ,0/             0001128
C FLORIDA                    NORTH                                     L0001129
       DATA ID(1, 24),ID(2, 24),ID(3, 24) /4HFLOR,4HIDA ,4H    /        0001130
       DATA ID(4, 24),ID(5, 24),ID(6, 24) /4H    ,4HNORT,4HH   /        0001131
       DATA ID(7, 24),ID(8, 24),ID(9, 24) /4H    ,4H    ,1/             0001132
C GEORGIA                    EAST                                      T0001133
       DATA ID(1, 25),ID(2, 25),ID(3, 25) /4HGEOR,4HGIA ,4H    /        0001134
       DATA ID(4, 25),ID(5, 25),ID(6, 25) /4H    ,4HEAST,4H    /        0001135
       DATA ID(7, 25),ID(8, 25),ID(9, 25) /4H    ,4H    ,0/             0001136
C GEORGIA                    WEST                                      T0001137
       DATA ID(1, 26),ID(2, 26),ID(3, 26) /4HGEOR,4HGIA ,4H    /        0001138
       DATA ID(4, 26),ID(5, 26),ID(6, 26) /4H    ,4HWEST,4H    /        0001139
       DATA ID(7, 26),ID(8, 26),ID(9, 26) /4H    ,4H    ,0/             0001140
C HAWAII                     1                                         T0001141
       DATA ID(1, 27),ID(2, 27),ID(3, 27) /4HHAWA,4HII  ,4H    /        0001142
       DATA ID(4, 27),ID(5, 27),ID(6, 27) /4H    ,4H1   ,4H    /        0001143
       DATA ID(7, 27),ID(8, 27),ID(9, 27) /4H    ,4H    ,0/             0001144
C HAWAII                     2                                         T0001145
       DATA ID(1, 28),ID(2, 28),ID(3, 28) /4HHAWA,4HII  ,4H    /        0001146
       DATA ID(4, 28),ID(5, 28),ID(6, 28) /4H    ,4H2   ,4H    /        0001147
       DATA ID(7, 28),ID(8, 28),ID(9, 28) /4H    ,4H    ,0/             0001148
C HAWAII                     3                                         T0001149
       DATA ID(1, 29),ID(2, 29),ID(3, 29) /4HHAWA,4HII  ,4H    /        0001150
       DATA ID(4, 29),ID(5, 29),ID(6, 29) /4H    ,4H3   ,4H    /        0001151
       DATA ID(7, 29),ID(8, 29),ID(9, 29) /4H    ,4H    ,0/             0001152
C HAWAII                     4                                         T0001153
       DATA ID(1, 30),ID(2, 30),ID(3, 30) /4HHAWA,4HII  ,4H    /        0001154
       DATA ID(4, 30),ID(5, 30),ID(6, 30) /4H    ,4H4   ,4H    /        0001155
       DATA ID(7, 30),ID(8, 30),ID(9, 30) /4H    ,4H    ,0/             0001156
C HAWAII                     5                                         T0001157
       DATA ID(1, 31),ID(2, 31),ID(3, 31) /4HHAWA,4HII  ,4H    /        0001158
       DATA ID(4, 31),ID(5, 31),ID(6, 31) /4H    ,4H5   ,4H    /        0001159
       DATA ID(7, 31),ID(8, 31),ID(9, 31) /4H    ,4H    ,0/             0001160
C IDAHO                      EAST                                      T0001161
       DATA ID(1, 32),ID(2, 32),ID(3, 32) /4HIDAH,4HO   ,4H    /        0001162
       DATA ID(4, 32),ID(5, 32),ID(6, 32) /4H    ,4HEAST,4H    /        0001163
       DATA ID(7, 32),ID(8, 32),ID(9, 32) /4H    ,4H    ,0/             0001164
C IDAHO                      CENTRAL                                   T0001165
       DATA ID(1, 33),ID(2, 33),ID(3, 33) /4HIDAH,4HO   ,4H    /        0001166
       DATA ID(4, 33),ID(5, 33),ID(6, 33) /4H    ,4HCENT,4HRAL /        0001167
       DATA ID(7, 33),ID(8, 33),ID(9, 33) /4H    ,4H    ,0/             0001168
C IDAHO                      WEST                                      T0001169
       DATA ID(1, 34),ID(2, 34),ID(3, 34) /4HIDAH,4HO   ,4H    /        0001170
       DATA ID(4, 34),ID(5, 34),ID(6, 34) /4H    ,4HWEST,4H    /        0001171
       DATA ID(7, 34),ID(8, 34),ID(9, 34) /4H    ,4H    ,0/             0001172
C ILLINOIS                   EAST                                      T0001173
       DATA ID(1, 35),ID(2, 35),ID(3, 35) /4HILLI,4HNOIS,4H    /        0001174
       DATA ID(4, 35),ID(5, 35),ID(6, 35) /4H    ,4HEAST,4H    /        0001175
       DATA ID(7, 35),ID(8, 35),ID(9, 35) /4H    ,4H    ,0/             0001176
C ILLINOIS                   WEST                                      T0001177
       DATA ID(1, 36),ID(2, 36),ID(3, 36) /4HILLI,4HNOIS,4H    /        0001178
       DATA ID(4, 36),ID(5, 36),ID(6, 36) /4H    ,4HWEST,4H    /        0001179
       DATA ID(7, 36),ID(8, 36),ID(9, 36) /4H    ,4H    ,0/             0001180
C INDIANA                    EAST                                      T0001181
       DATA ID(1, 37),ID(2, 37),ID(3, 37) /4HINDI,4HANA ,4H    /        0001182
       DATA ID(4, 37),ID(5, 37),ID(6, 37) /4H    ,4HEAST,4H    /        0001183
       DATA ID(7, 37),ID(8, 37),ID(9, 37) /4H    ,4H    ,0/             0001184
C INDIANA                    WEST                                      T0001185
       DATA ID(1, 38),ID(2, 38),ID(3, 38) /4HINDI,4HANA ,4H    /        0001186
       DATA ID(4, 38),ID(5, 38),ID(6, 38) /4H    ,4HWEST,4H    /        0001187
       DATA ID(7, 38),ID(8, 38),ID(9, 38) /4H    ,4H    ,0/             0001188
C IOWA                       NORTH                                     L0001189
       DATA ID(1, 39),ID(2, 39),ID(3, 39) /4HIOWA,4H    ,4H    /        0001190
       DATA ID(4, 39),ID(5, 39),ID(6, 39) /4H    ,4HNORT,4HH   /        0001191
       DATA ID(7, 39),ID(8, 39),ID(9, 39) /4H    ,4H    ,1/             0001192
C IOWA                       SOUTH                                     L0001193
       DATA ID(1, 40),ID(2, 40),ID(3, 40) /4HIOWA,4H    ,4H    /        0001194
       DATA ID(4, 40),ID(5, 40),ID(6, 40) /4H    ,4HSOUT,4HH   /        0001195
       DATA ID(7, 40),ID(8, 40),ID(9, 40) /4H    ,4H    ,1/             0001196
C KANSAS                     NORTH                                     L0001197
       DATA ID(1, 41),ID(2, 41),ID(3, 41) /4HKANS,4HAS  ,4H    /        0001198
       DATA ID(4, 41),ID(5, 41),ID(6, 41) /4H    ,4HNORT,4HH   /        0001199
       DATA ID(7, 41),ID(8, 41),ID(9, 41) /4H    ,4H    ,1/             0001200
C KANSAS                     SOUTH                                     L0001201
       DATA ID(1, 42),ID(2, 42),ID(3, 42) /4HKANS,4HAS  ,4H    /        0001202
       DATA ID(4, 42),ID(5, 42),ID(6, 42) /4H    ,4HSOUT,4HH   /        0001203
       DATA ID(7, 42),ID(8, 42),ID(9, 42) /4H    ,4H    ,1/             0001204
C KENTUCKY                   NORTH                                     L0001205
       DATA ID(1, 43),ID(2, 43),ID(3, 43) /4HKENT,4HUCKY,4H    /        0001206
       DATA ID(4, 43),ID(5, 43),ID(6, 43) /4H    ,4HNORT,4HH   /        0001207
       DATA ID(7, 43),ID(8, 43),ID(9, 43) /4H    ,4H    ,1/             0001208
C KENTUCKY                   SOUTH                                     L0001209
       DATA ID(1, 44),ID(2, 44),ID(3, 44) /4HKENT,4HUCKY,4H    /        0001210
       DATA ID(4, 44),ID(5, 44),ID(6, 44) /4H    ,4HSOUT,4HH   /        0001211
       DATA ID(7, 44),ID(8, 44),ID(9, 44) /4H    ,4H    ,1/             0001212
C LOUISIANA                  NORTH                                     L0001213
       DATA ID(1, 45),ID(2, 45),ID(3, 45) /4HLOUI,4HSIAN,4HA   /        0001214
       DATA ID(4, 45),ID(5, 45),ID(6, 45) /4H    ,4HNORT,4HH   /        0001215
       DATA ID(7, 45),ID(8, 45),ID(9, 45) /4H    ,4H    ,1/             0001216
C LOUISIANA                  SOUTH                                     L0001217
       DATA ID(1, 46),ID(2, 46),ID(3, 46) /4HLOUI,4HSIAN,4HA   /        0001218
       DATA ID(4, 46),ID(5, 46),ID(6, 46) /4H    ,4HSOUT,4HH   /        0001219
       DATA ID(7, 46),ID(8, 46),ID(9, 46) /4H    ,4H    ,1/             0001220
C LOUISIANA                  OFFSHORE                                  L0001221
       DATA ID(1, 47),ID(2, 47),ID(3, 47) /4HLOUI,4HSIAN,4HA   /        0001222
       DATA ID(4, 47),ID(5, 47),ID(6, 47) /4H    ,4HOFFS,4HHORE/        0001223
       DATA ID(7, 47),ID(8, 47),ID(9, 47) /4H    ,4H    ,1/             0001224
C MAINE                      EAST                                      T0001225
       DATA ID(1, 48),ID(2, 48),ID(3, 48) /4HMAIN,4HE   ,4H    /        0001226
       DATA ID(4, 48),ID(5, 48),ID(6, 48) /4H    ,4HEAST,4H    /        0001227
       DATA ID(7, 48),ID(8, 48),ID(9, 48) /4H    ,4H    ,0/             0001228
C MAINE                      WEST                                      T0001229
       DATA ID(1, 49),ID(2, 49),ID(3, 49) /4HMAIN,4HE   ,4H    /        0001230
       DATA ID(4, 49),ID(5, 49),ID(6, 49) /4H    ,4HWEST,4H    /        0001231
       DATA ID(7, 49),ID(8, 49),ID(9, 49) /4H    ,4H    ,0/             0001232
C MARYLAND                   ---                                       L0001233
       DATA ID(1, 50),ID(2, 50),ID(3, 50) /4HMARY,4HLAND,4H    /        0001234
       DATA ID(4, 50),ID(5, 50),ID(6, 50) /4H    ,4H--- ,4H    /        0001235
       DATA ID(7, 50),ID(8, 50),ID(9, 50) /4H    ,4H    ,1/             0001236
C MASSACHUSETTS              MAINLAND                                  L0001237
       DATA ID(1, 51),ID(2, 51),ID(3, 51) /4HMASS,4HACHU,4HSETT/        0001238
       DATA ID(4, 51),ID(5, 51),ID(6, 51) /4HS   ,4HMAIN,4HLAND/        0001239
       DATA ID(7, 51),ID(8, 51),ID(9, 51) /4H    ,4H    ,1/             0001240
C MASSACHUSETTS              ISLAND                                    L0001241
       DATA ID(1, 52),ID(2, 52),ID(3, 52) /4HMASS,4HACHU,4HSETT/        0001242
       DATA ID(4, 52),ID(5, 52),ID(6, 52) /4HS   ,4HISLA,4HND  /        0001243
       DATA ID(7, 52),ID(8, 52),ID(9, 52) /4H    ,4H    ,1/             0001244
C MICHIGAN                   EAST                                      T0001245
       DATA ID(1, 53),ID(2, 53),ID(3, 53) /4HMICH,4HIGAN,4H    /        0001246
       DATA ID(4, 53),ID(5, 53),ID(6, 53) /4H    ,4HEAST,4H    /        0001247
       DATA ID(7, 53),ID(8, 53),ID(9, 53) /4H    ,4H    ,0/             0001248
C MICHIGAN                   CENTRAL                                   T0001249
       DATA ID(1, 54),ID(2, 54),ID(3, 54) /4HMICH,4HIGAN,4H    /        0001250
       DATA ID(4, 54),ID(5, 54),ID(6, 54) /4H    ,4HCENT,4HRAL /        0001251
       DATA ID(7, 54),ID(8, 54),ID(9, 54) /4H    ,4H    ,0/             0001252
C MICHIGAN                   WEST                                      T0001253
       DATA ID(1, 55),ID(2, 55),ID(3, 55) /4HMICH,4HIGAN,4H    /        0001254
       DATA ID(4, 55),ID(5, 55),ID(6, 55) /4H    ,4HWEST,4H    /        0001255
       DATA ID(7, 55),ID(8, 55),ID(9, 55) /4H    ,4H    ,0/             0001256
C MICHIGAN                   NORTH                                     L0001257
       DATA ID(1, 56),ID(2, 56),ID(3, 56) /4HMICH,4HIGAN,4H    /        0001258
       DATA ID(4, 56),ID(5, 56),ID(6, 56) /4H    ,4HNORT,4HH   /        0001259
       DATA ID(7, 56),ID(8, 56),ID(9, 56) /4H    ,4H    ,1/             0001260
C MICHIGAN                   CENTRAL                                   L0001261
       DATA ID(1, 57),ID(2, 57),ID(3, 57) /4HMICH,4HIGAN,4H    /        0001262
       DATA ID(4, 57),ID(5, 57),ID(6, 57) /4H    ,4HCENT,4HRAL /        0001263
       DATA ID(7, 57),ID(8, 57),ID(9, 57) /4H    ,4H    ,1/             0001264
C MICHIGAN                   SOUTH                                     L0001265
       DATA ID(1, 58),ID(2, 58),ID(3, 58) /4HMICH,4HIGAN,4H    /        0001266
       DATA ID(4, 58),ID(5, 58),ID(6, 58) /4H    ,4HSOUT,4HH   /        0001267
       DATA ID(7, 58),ID(8, 58),ID(9, 58) /4H    ,4H    ,1/             0001268
C MINNESOTA                  NORTH                                     L0001269
       DATA ID(1, 59),ID(2, 59),ID(3, 59) /4HMINN,4HESOT,4HA   /        0001270
       DATA ID(4, 59),ID(5, 59),ID(6, 59) /4H    ,4HNORT,4HH   /        0001271
       DATA ID(7, 59),ID(8, 59),ID(9, 59) /4H    ,4H    ,1/             0001272
C MINNESOTA                  CENTRAL                                   L0001273
       DATA ID(1, 60),ID(2, 60),ID(3, 60) /4HMINN,4HESOT,4HA   /        0001274
       DATA ID(4, 60),ID(5, 60),ID(6, 60) /4H    ,4HCENT,4HRAL /        0001275
       DATA ID(7, 60),ID(8, 60),ID(9, 60) /4H    ,4H    ,1/             0001276
C MINNESOTA                  SOUTH                                     L0001277
       DATA ID(1, 61),ID(2, 61),ID(3, 61) /4HMINN,4HESOT,4HA   /        0001278
       DATA ID(4, 61),ID(5, 61),ID(6, 61) /4H    ,4HSOUT,4HH   /        0001279
       DATA ID(7, 61),ID(8, 61),ID(9, 61) /4H    ,4H    ,1/             0001280
C MISSISSIPPI                EAST                                      T0001281
       DATA ID(1, 62),ID(2, 62),ID(3, 62) /4HMISS,4HISSI,4HPPI /        0001282
       DATA ID(4, 62),ID(5, 62),ID(6, 62) /4H    ,4HEAST,4H    /        0001283
       DATA ID(7, 62),ID(8, 62),ID(9, 62) /4H    ,4H    ,0/             0001284
C MISSISSIPPI                WEST                                      T0001285
       DATA ID(1, 63),ID(2, 63),ID(3, 63) /4HMISS,4HISSI,4HPPI /        0001286
       DATA ID(4, 63),ID(5, 63),ID(6, 63) /4H    ,4HWEST,4H    /        0001287
       DATA ID(7, 63),ID(8, 63),ID(9, 63) /4H    ,4H    ,0/             0001288
C MISSOURI                   EAST                                      T0001289
       DATA ID(1, 64),ID(2, 64),ID(3, 64) /4HMISS,4HOURI,4H    /        0001290
       DATA ID(4, 64),ID(5, 64),ID(6, 64) /4H    ,4HEAST,4H    /        0001291
       DATA ID(7, 64),ID(8, 64),ID(9, 64) /4H    ,4H    ,0/             0001292
C MISSOURI                   CENTRAL                                   T0001293
       DATA ID(1, 65),ID(2, 65),ID(3, 65) /4HMISS,4HOURI,4H    /        0001294
       DATA ID(4, 65),ID(5, 65),ID(6, 65) /4H    ,4HCENT,4HRAL /        0001295
       DATA ID(7, 65),ID(8, 65),ID(9, 65) /4H    ,4H    ,0/             0001296
C MISSOURI                   WEST                                      T0001297
       DATA ID(1, 66),ID(2, 66),ID(3, 66) /4HMISS,4HOURI,4H    /        0001298
       DATA ID(4, 66),ID(5, 66),ID(6, 66) /4H    ,4HWEST,4H    /        0001299
       DATA ID(7, 66),ID(8, 66),ID(9, 66) /4H    ,4H    ,0/             0001300
C MONTANA                    NORTH                                     L0001301
       DATA ID(1, 67),ID(2, 67),ID(3, 67) /4HMONT,4HANA ,4H    /        0001302
       DATA ID(4, 67),ID(5, 67),ID(6, 67) /4H    ,4HNORT,4HH   /        0001303
       DATA ID(7, 67),ID(8, 67),ID(9, 67) /4H    ,4H    ,1/             0001304
C MONTANA                    CENTRAL                                   L0001305
       DATA ID(1, 68),ID(2, 68),ID(3, 68) /4HMONT,4HANA ,4H    /        0001306
       DATA ID(4, 68),ID(5, 68),ID(6, 68) /4H    ,4HCENT,4HRAL /        0001307
       DATA ID(7, 68),ID(8, 68),ID(9, 68) /4H    ,4H    ,1/             0001308
C MONTANA                    SOUTH                                     L0001309
       DATA ID(1, 69),ID(2, 69),ID(3, 69) /4HMONT,4HANA ,4H    /        0001310
       DATA ID(4, 69),ID(5, 69),ID(6, 69) /4H    ,4HSOUT,4HH   /        0001311
       DATA ID(7, 69),ID(8, 69),ID(9, 69) /4H    ,4H    ,1/             0001312
C NEBRASKA                   NORTH                                     L0001313
       DATA ID(1, 70),ID(2, 70),ID(3, 70) /4HNEBR,4HASKA,4H    /        0001314
       DATA ID(4, 70),ID(5, 70),ID(6, 70) /4H    ,4HNORT,4HH   /        0001315
       DATA ID(7, 70),ID(8, 70),ID(9, 70) /4H    ,4H    ,1/             0001316
C NEBRASKA                   SOUTH                                     L0001317
       DATA ID(1, 71),ID(2, 71),ID(3, 71) /4HNEBR,4HASKA,4H    /        0001318
       DATA ID(4, 71),ID(5, 71),ID(6, 71) /4H    ,4HSOUT,4HH   /        0001319
       DATA ID(7, 71),ID(8, 71),ID(9, 71) /4H    ,4H    ,1/             0001320
C NEVADA                     EAST                                      T0001321
       DATA ID(1, 72),ID(2, 72),ID(3, 72) /4HNEVA,4HDA  ,4H    /        0001322
       DATA ID(4, 72),ID(5, 72),ID(6, 72) /4H    ,4HEAST,4H    /        0001323
       DATA ID(7, 72),ID(8, 72),ID(9, 72) /4H    ,4H    ,0/             0001324
C NEVADA                     CENTRAL                                   T0001325
       DATA ID(1, 73),ID(2, 73),ID(3, 73) /4HNEVA,4HDA  ,4H    /        0001326
       DATA ID(4, 73),ID(5, 73),ID(6, 73) /4H    ,4HCENT,4HRAL /        0001327
       DATA ID(7, 73),ID(8, 73),ID(9, 73) /4H    ,4H    ,0/             0001328
C NEVADA                     WEST                                      T0001329
       DATA ID(1, 74),ID(2, 74),ID(3, 74) /4HNEVA,4HDA  ,4H    /        0001330
       DATA ID(4, 74),ID(5, 74),ID(6, 74) /4H    ,4HWEST,4H    /        0001331
       DATA ID(7, 74),ID(8, 74),ID(9, 74) /4H    ,4H    ,0/             0001332
C NEW HAMPSHIRE              ---                                       T0001333
       DATA ID(1, 75),ID(2, 75),ID(3, 75) /4HNEW ,4HHAMP,4HSHIR/        0001334
       DATA ID(4, 75),ID(5, 75),ID(6, 75) /4HE   ,4H--- ,4H    /        0001335
       DATA ID(7, 75),ID(8, 75),ID(9, 75) /4H    ,4H    ,0/             0001336
C NEW JERSEY                 ---                                       T0001337
       DATA ID(1, 76),ID(2, 76),ID(3, 76) /4HNEW ,4HJERS,4HEY  /        0001338
       DATA ID(4, 76),ID(5, 76),ID(6, 76) /4H    ,4H--- ,4H    /        0001339
       DATA ID(7, 76),ID(8, 76),ID(9, 76) /4H    ,4H    ,0/             0001340
C NEW MEXICO                 EAST                                      T0001341
       DATA ID(1, 77),ID(2, 77),ID(3, 77) /4HNEW ,4HMEXI,4HCO  /        0001342
       DATA ID(4, 77),ID(5, 77),ID(6, 77) /4H    ,4HEAST,4H    /        0001343
       DATA ID(7, 77),ID(8, 77),ID(9, 77) /4H    ,4H    ,0/             0001344
C NEW MEXICO                 CENTRAL                                   T0001345
       DATA ID(1, 78),ID(2, 78),ID(3, 78) /4HNEW ,4HMEXI,4HCO  /        0001346
       DATA ID(4, 78),ID(5, 78),ID(6, 78) /4H    ,4HCENT,4HRAL /        0001347
       DATA ID(7, 78),ID(8, 78),ID(9, 78) /4H    ,4H    ,0/             0001348
C NEW MEXICO                 WEST                                      T0001349
       DATA ID(1, 79),ID(2, 79),ID(3, 79) /4HNEW ,4HMEXI,4HCO  /        0001350
       DATA ID(4, 79),ID(5, 79),ID(6, 79) /4H    ,4HWEST,4H    /        0001351
       DATA ID(7, 79),ID(8, 79),ID(9, 79) /4H    ,4H    ,0/             0001352
C NEW YORK                   EAST                                      T0001353
       DATA ID(1, 80),ID(2, 80),ID(3, 80) /4HNEW ,4HYORK,4H    /        0001354
       DATA ID(4, 80),ID(5, 80),ID(6, 80) /4H    ,4HEAST,4H    /        0001355
       DATA ID(7, 80),ID(8, 80),ID(9, 80) /4H    ,4H    ,0/             0001356
C NEW YORK                   CENTRAL                                   T0001357
       DATA ID(1, 81),ID(2, 81),ID(3, 81) /4HNEW ,4HYORK,4H    /        0001358
       DATA ID(4, 81),ID(5, 81),ID(6, 81) /4H    ,4HCENT,4HRAL /        0001359
       DATA ID(7, 81),ID(8, 81),ID(9, 81) /4H    ,4H    ,0/             0001360
C NEW YORK                   WEST                                      T0001361
       DATA ID(1, 82),ID(2, 82),ID(3, 82) /4HNEW ,4HYORK,4H    /        0001362
       DATA ID(4, 82),ID(5, 82),ID(6, 82) /4H    ,4HWEST,4H    /        0001363
       DATA ID(7, 82),ID(8, 82),ID(9, 82) /4H    ,4H    ,0/             0001364
C NEW YORK                   LONG ISLAND                               L0001365
       DATA ID(1, 83),ID(2, 83),ID(3, 83) /4HNEW ,4HYORK,4H    /        0001366
       DATA ID(4, 83),ID(5, 83),ID(6, 83) /4H    ,4HLONG,4H ISL/        0001367
       DATA ID(7, 83),ID(8, 83),ID(9, 83) /4HAND ,4H    ,1/             0001368
C NORTH CAROLINA             ---                                       L0001369
       DATA ID(1, 84),ID(2, 84),ID(3, 84) /4HNORT,4HH CA,4HROLI/        0001370
       DATA ID(4, 84),ID(5, 84),ID(6, 84) /4HNA  ,4H--- ,4H    /        0001371
       DATA ID(7, 84),ID(8, 84),ID(9, 84) /4H    ,4H    ,1/             0001372
C NORTH DAKOTA               NORTH                                     L0001373
       DATA ID(1, 85),ID(2, 85),ID(3, 85) /4HNORT,4HH DA,4HKOTA/        0001374
       DATA ID(4, 85),ID(5, 85),ID(6, 85) /4H    ,4HNORT,4HH   /        0001375
       DATA ID(7, 85),ID(8, 85),ID(9, 85) /4H    ,4H    ,1/             0001376
C NORTH DAKOTA               SOUTH                                     L0001377
       DATA ID(1, 86),ID(2, 86),ID(3, 86) /4HNORT,4HH DA,4HKOTA/        0001378
       DATA ID(4, 86),ID(5, 86),ID(6, 86) /4H    ,4HSOUT,4HH   /        0001379
       DATA ID(7, 86),ID(8, 86),ID(9, 86) /4H    ,4H    ,1/             0001380
C OHIO                       NORTH                                     L0001381
       DATA ID(1, 87),ID(2, 87),ID(3, 87) /4HOHIO,4H    ,4H    /        0001382
       DATA ID(4, 87),ID(5, 87),ID(6, 87) /4H    ,4HNORT,4HH   /        0001383
       DATA ID(7, 87),ID(8, 87),ID(9, 87) /4H    ,4H    ,1/             0001384
C OHIO                       SOUTH                                     L0001385
       DATA ID(1, 88),ID(2, 88),ID(3, 88) /4HOHIO,4H    ,4H    /        0001386
       DATA ID(4, 88),ID(5, 88),ID(6, 88) /4H    ,4HSOUT,4HH   /        0001387
       DATA ID(7, 88),ID(8, 88),ID(9, 88) /4H    ,4H    ,1/             0001388
C OKLAHOMA                   NORTH                                     L0001389
       DATA ID(1, 89),ID(2, 89),ID(3, 89) /4HOKLA,4HHOMA,4H    /        0001390
       DATA ID(4, 89),ID(5, 89),ID(6, 89) /4H    ,4HNORT,4HH   /        0001391
       DATA ID(7, 89),ID(8, 89),ID(9, 89) /4H    ,4H    ,1/             0001392
C OKLAHOMA                   SOUTH                                     L0001393
       DATA ID(1, 90),ID(2, 90),ID(3, 90) /4HOKLA,4HHOMA,4H    /        0001394
       DATA ID(4, 90),ID(5, 90),ID(6, 90) /4H    ,4HSOUT,4HH   /        0001395
       DATA ID(7, 90),ID(8, 90),ID(9, 90) /4H    ,4H    ,1/             0001396
C OREGON                     NORTH                                     L0001397
       DATA ID(1, 91),ID(2, 91),ID(3, 91) /4HOREG,4HON  ,4H    /        0001398
       DATA ID(4, 91),ID(5, 91),ID(6, 91) /4H    ,4HNORT,4HH   /        0001399
       DATA ID(7, 91),ID(8, 91),ID(9, 91) /4H    ,4H    ,1/             0001400
C OREGON                     SOUTH                                     L0001401
       DATA ID(1, 92),ID(2, 92),ID(3, 92) /4HOREG,4HON  ,4H    /        0001402
       DATA ID(4, 92),ID(5, 92),ID(6, 92) /4H    ,4HSOUT,4HH   /        0001403
       DATA ID(7, 92),ID(8, 92),ID(9, 92) /4H    ,4H    ,1/             0001404
C PENNSYLVANIA               NORTH                                     L0001405
       DATA ID(1, 93),ID(2, 93),ID(3, 93) /4HPENN,4HSYLV,4HANIA/        0001406
       DATA ID(4, 93),ID(5, 93),ID(6, 93) /4H    ,4HNORT,4HH   /        0001407
       DATA ID(7, 93),ID(8, 93),ID(9, 93) /4H    ,4H    ,1/             0001408
C PENNSYLVANIA               SOUTH                                     L0001409
       DATA ID(1, 94),ID(2, 94),ID(3, 94) /4HPENN,4HSYLV,4HANIA/        0001410
       DATA ID(4, 94),ID(5, 94),ID(6, 94) /4H    ,4HSOUT,4HH   /        0001411
       DATA ID(7, 94),ID(8, 94),ID(9, 94) /4H    ,4H    ,1/             0001412
C RHODE ISLAND               ---                                       T0001413
       DATA ID(1, 95),ID(2, 95),ID(3, 95) /4HRHOD,4HE IS,4HLAND/        0001414
       DATA ID(4, 95),ID(5, 95),ID(6, 95) /4H    ,4H--- ,4H    /        0001415
       DATA ID(7, 95),ID(8, 95),ID(9, 95) /4H    ,4H    ,0/             0001416
C SOUTH CAROLINA             NORTH                                     L0001417
       DATA ID(1, 96),ID(2, 96),ID(3, 96) /4HSOUT,4HH CA,4HROLI/        0001418
       DATA ID(4, 96),ID(5, 96),ID(6, 96) /4HNA  ,4HNORT,4HH   /        0001419
       DATA ID(7, 96),ID(8, 96),ID(9, 96) /4H    ,4H    ,1/             0001420
C SOUTH CAROLINA             SOUTH                                     L0001421
       DATA ID(1, 97),ID(2, 97),ID(3, 97) /4HSOUT,4HH CA,4HROLI/        0001422
       DATA ID(4, 97),ID(5, 97),ID(6, 97) /4HNA  ,4HSOUT,4HH   /        0001423
       DATA ID(7, 97),ID(8, 97),ID(9, 97) /4H    ,4H    ,1/             0001424
C SOUTH DAKOTA               NORTH                                     L0001425
       DATA ID(1, 98),ID(2, 98),ID(3, 98) /4HSOUT,4HH DA,4HKOTA/        0001426
       DATA ID(4, 98),ID(5, 98),ID(6, 98) /4H    ,4HNORT,4HH   /        0001427
       DATA ID(7, 98),ID(8, 98),ID(9, 98) /4H    ,4H    ,1/             0001428
C SOUTH DAKOTA               SOUTH                                     L0001429
       DATA ID(1, 99),ID(2, 99),ID(3, 99) /4HSOUT,4HH DA,4HKOTA/        0001430
       DATA ID(4, 99),ID(5, 99),ID(6, 99) /4H    ,4HSOUT,4HH   /        0001431
       DATA ID(7, 99),ID(8, 99),ID(9, 99) /4H    ,4H    ,1/             0001432
C TENNESSEE                  ---                                       L0001433
       DATA ID(1,100),ID(2,100),ID(3,100) /4HTENN,4HESSE,4HE   /        0001434
       DATA ID(4,100),ID(5,100),ID(6,100) /4H    ,4H--- ,4H    /        0001435
       DATA ID(7,100),ID(8,100),ID(9,100) /4H    ,4H    ,1/             0001436
C TEXAS                      NORTH                                     L0001437
       DATA ID(1,101),ID(2,101),ID(3,101) /4HTEXA,4HS   ,4H    /        0001438
       DATA ID(4,101),ID(5,101),ID(6,101) /4H    ,4HNORT,4HH   /        0001439
       DATA ID(7,101),ID(8,101),ID(9,101) /4H    ,4H    ,1/             0001440
C TEXAS                      NORTH CENTRAL                             L0001441
       DATA ID(1,102),ID(2,102),ID(3,102) /4HTEXA,4HS   ,4H    /        0001442
       DATA ID(4,102),ID(5,102),ID(6,102) /4H    ,4HNORT,4HH CE/        0001443
       DATA ID(7,102),ID(8,102),ID(9,102) /4HNTRA,4HL   ,1/             0001444
C TEXAS                      CENTRAL                                   L0001445
       DATA ID(1,103),ID(2,103),ID(3,103) /4HTEXA,4HS   ,4H    /        0001446
       DATA ID(4,103),ID(5,103),ID(6,103) /4H    ,4HCENT,4HRAL /        0001447
       DATA ID(7,103),ID(8,103),ID(9,103) /4H    ,4H    ,1/             0001448
C TEXAS                      SOUTH CENTRAL                             L0001449
       DATA ID(1,104),ID(2,104),ID(3,104) /4HTEXA,4HS   ,4H    /        0001450
       DATA ID(4,104),ID(5,104),ID(6,104) /4H    ,4HSOUT,4HH CE/        0001451
       DATA ID(7,104),ID(8,104),ID(9,104) /4HNTRA,4HL   ,1/             0001452
C TEXAS                      SOUTH                                     L0001453
       DATA ID(1,105),ID(2,105),ID(3,105) /4HTEXA,4HS   ,4H    /        0001454
       DATA ID(4,105),ID(5,105),ID(6,105) /4H    ,4HSOUT,4HH   /        0001455
       DATA ID(7,105),ID(8,105),ID(9,105) /4H    ,4H    ,1/             0001456
C UTAH                       NORTH                                     L0001457
       DATA ID(1,106),ID(2,106),ID(3,106) /4HUTAH,4H    ,4H    /        0001458
       DATA ID(4,106),ID(5,106),ID(6,106) /4H    ,4HNORT,4HH   /        0001459
       DATA ID(7,106),ID(8,106),ID(9,106) /4H    ,4H    ,1/             0001460
C UTAH                       CENTRAL                                   L0001461
       DATA ID(1,107),ID(2,107),ID(3,107) /4HUTAH,4H    ,4H    /        0001462
       DATA ID(4,107),ID(5,107),ID(6,107) /4H    ,4HCENT,4HRAL /        0001463
       DATA ID(7,107),ID(8,107),ID(9,107) /4H    ,4H    ,1/             0001464
C UTAH                       SOUTH                                     L0001465
       DATA ID(1,108),ID(2,108),ID(3,108) /4HUTAH,4H    ,4H    /        0001466
       DATA ID(4,108),ID(5,108),ID(6,108) /4H    ,4HSOUT,4HH   /        0001467
       DATA ID(7,108),ID(8,108),ID(9,108) /4H    ,4H    ,1/             0001468
C VERMONT                    ---                                       T0001469
       DATA ID(1,109),ID(2,109),ID(3,109) /4HVERM,4HONT ,4H    /        0001470
       DATA ID(4,109),ID(5,109),ID(6,109) /4H    ,4H--- ,4H    /        0001471
       DATA ID(7,109),ID(8,109),ID(9,109) /4H    ,4H    ,0/             0001472
C VIRGINIA                   NORTH                                     L0001473
       DATA ID(1,110),ID(2,110),ID(3,110) /4HVIRG,4HINIA,4H    /        0001474
       DATA ID(4,110),ID(5,110),ID(6,110) /4H    ,4HNORT,4HH   /        0001475
       DATA ID(7,110),ID(8,110),ID(9,110) /4H    ,4H    ,1/             0001476
C VIRGINIA                   SOUTH                                     L0001477
       DATA ID(1,111),ID(2,111),ID(3,111) /4HVIRG,4HINIA,4H    /        0001478
       DATA ID(4,111),ID(5,111),ID(6,111) /4H    ,4HSOUT,4HH   /        0001479
       DATA ID(7,111),ID(8,111),ID(9,111) /4H    ,4H    ,1/             0001480
C WASHINGTON                 NORTH                                     L0001481
       DATA ID(1,112),ID(2,112),ID(3,112) /4HWASH,4HINGT,4HON  /        0001482
       DATA ID(4,112),ID(5,112),ID(6,112) /4H    ,4HNORT,4HH   /        0001483
       DATA ID(7,112),ID(8,112),ID(9,112) /4H    ,4H    ,1/             0001484
C WASHINGTON                 SOUTH                                     L0001485
       DATA ID(1,113),ID(2,113),ID(3,113) /4HWASH,4HINGT,4HON  /        0001486
       DATA ID(4,113),ID(5,113),ID(6,113) /4H    ,4HSOUT,4HH   /        0001487
       DATA ID(7,113),ID(8,113),ID(9,113) /4H    ,4H    ,1/             0001488
C WEST VIRGINIA              NORTH                                     L0001489
       DATA ID(1,114),ID(2,114),ID(3,114) /4HWEST,4H VIR,4HGINI/        0001490
       DATA ID(4,114),ID(5,114),ID(6,114) /4HA   ,4HNORT,4HH   /        0001491
       DATA ID(7,114),ID(8,114),ID(9,114) /4H    ,4H    ,1/             0001492
C WEST VIRGINIA              SOUTH                                     L0001493
       DATA ID(1,115),ID(2,115),ID(3,115) /4HWEST,4H VIR,4HGINI/        0001494
       DATA ID(4,115),ID(5,115),ID(6,115) /4HA   ,4HSOUT,4HH   /        0001495
       DATA ID(7,115),ID(8,115),ID(9,115) /4H    ,4H    ,1/             0001496
C WISCONSIN                  NORTH                                     L0001497
       DATA ID(1,116),ID(2,116),ID(3,116) /4HWISC,4HONSI,4HN   /        0001498
       DATA ID(4,116),ID(5,116),ID(6,116) /4H    ,4HNORT,4HH   /        0001499
       DATA ID(7,116),ID(8,116),ID(9,116) /4H    ,4H    ,1/             0001500
C WISCONSIN                  CENTRAL                                   L0001501
       DATA ID(1,117),ID(2,117),ID(3,117) /4HWISC,4HONSI,4HN   /        0001502
       DATA ID(4,117),ID(5,117),ID(6,117) /4H    ,4HCENT,4HRAL /        0001503
       DATA ID(7,117),ID(8,117),ID(9,117) /4H    ,4H    ,1/             0001504
C WISCONSIN                  SOUTH                                     L0001505
       DATA ID(1,118),ID(2,118),ID(3,118) /4HWISC,4HONSI,4HN   /        0001506
       DATA ID(4,118),ID(5,118),ID(6,118) /4H    ,4HSOUT,4HH   /        0001507
       DATA ID(7,118),ID(8,118),ID(9,118) /4H    ,4H    ,1/             0001508
C WYOMING                    EAST                                      T0001509
       DATA ID(1,119),ID(2,119),ID(3,119) /4HWYOM,4HING ,4H    /        0001510
       DATA ID(4,119),ID(5,119),ID(6,119) /4H    ,4HEAST,4H    /        0001511
       DATA ID(7,119),ID(8,119),ID(9,119) /4H    ,4H    ,0/             0001512
C WYOMING                    EAST CENTRAL                              T0001513
       DATA ID(1,120),ID(2,120),ID(3,120) /4HWYOM,4HING ,4H    /        0001514
       DATA ID(4,120),ID(5,120),ID(6,120) /4H    ,4HEAST,4H CEN/        0001515
       DATA ID(7,120),ID(8,120),ID(9,120) /4HTRAL,4H    ,0/             0001516
C WYOMING                    WEST CENTRAL                              T0001517
       DATA ID(1,121),ID(2,121),ID(3,121) /4HWYOM,4HING ,4H    /        0001518
       DATA ID(4,121),ID(5,121),ID(6,121) /4H    ,4HWEST,4H CEN/        0001519
       DATA ID(7,121),ID(8,121),ID(9,121) /4HTRAL,4H    ,0/             0001520
C WYOMING                    WEST                                      T0001521
       DATA ID(1,122),ID(2,122),ID(3,122) /4HWYOM,4HING ,4H    /        0001522
       DATA ID(4,122),ID(5,122),ID(6,122) /4H    ,4HWEST,4H    /        0001523
       DATA ID(7,122),ID(8,122),ID(9,122) /4H    ,4H    ,0/             0001524
C ALASKA                     ZONE NO. 1                                 0001525
       DATA ID(1,123),ID(2,123),ID(3,123) /4HALAS,4HKA  ,4H    /        0001526
       DATA ID(4,123),ID(5,123),ID(6,123) /4H    ,4HZONE,4H NO./        0001527
       DATA ID(7,123),ID(8,123),ID(9,123) /4H 1  ,4H    ,2/             0001528
C ALASKA                     ZONE NO. 2                                 0001529
       DATA ID(1,124),ID(2,124),ID(3,124) /4HALAS,4HKA  ,4H    /        0001530
       DATA ID(4,124),ID(5,124),ID(6,124) /4H    ,4HZONE,4H NO./        0001531
       DATA ID(7,124),ID(8,124),ID(9,124) /4H 2  ,4H    ,2/             0001532
C ALASKA                     ZONE NO. 3                                 0001533
       DATA ID(1,125),ID(2,125),ID(3,125) /4HALAS,4HKA  ,4H    /        0001534
       DATA ID(4,125),ID(5,125),ID(6,125) /4H    ,4HZONE,4H NO./        0001535
       DATA ID(7,125),ID(8,125),ID(9,125) /4H 3  ,4H    ,2/             0001536
C ALASKA                     ZONE NO. 4                                 0001537
       DATA ID(1,126),ID(2,126),ID(3,126) /4HALAS,4HKA  ,4H    /        0001538
       DATA ID(4,126),ID(5,126),ID(6,126) /4H    ,4HZONE,4H NO./        0001539
       DATA ID(7,126),ID(8,126),ID(9,126) /4H 4  ,4H    ,2/             0001540
C ALASKA                     ZONE NO. 5                                 0001541
       DATA ID(1,127),ID(2,127),ID(3,127) /4HALAS,4HKA  ,4H    /        0001542
       DATA ID(4,127),ID(5,127),ID(6,127) /4H    ,4HZONE,4H NO./        0001543
       DATA ID(7,127),ID(8,127),ID(9,127) /4H 5  ,4H    ,2/             0001544
C ALASKA                     ZONE NO. 6                                 0001545
       DATA ID(1,128),ID(2,128),ID(3,128) /4HALAS,4HKA  ,4H    /        0001546
       DATA ID(4,128),ID(5,128),ID(6,128) /4H    ,4HZONE,4H NO./        0001547
       DATA ID(7,128),ID(8,128),ID(9,128) /4H 6  ,4H    ,2/             0001548
C ALASKA                     ZONE NO. 7                                 0001549
       DATA ID(1,129),ID(2,129),ID(3,129) /4HALAS,4HKA  ,4H    /        0001550
       DATA ID(4,129),ID(5,129),ID(6,129) /4H    ,4HZONE,4H NO./        0001551
       DATA ID(7,129),ID(8,129),ID(9,129) /4H 7  ,4H    ,2/             0001552
C ALASKA                     ZONE NO. 8                                 0001553
       DATA ID(1,130),ID(2,130),ID(3,130) /4HALAS,4HKA  ,4H    /        0001554
       DATA ID(4,130),ID(5,130),ID(6,130) /4H    ,4HZONE,4H NO./        0001555
       DATA ID(7,130),ID(8,130),ID(9,130) /4H 8  ,4H    ,2/             0001556
C ALASKA                     ZONE NO. 9                                 0001557
       DATA ID(1,131),ID(2,131),ID(3,131) /4HALAS,4HKA  ,4H    /        0001558
       DATA ID(4,131),ID(5,131),ID(6,131) /4H    ,4HZONE,4H NO./        0001559
       DATA ID(7,131),ID(8,131),ID(9,131) /4H 9  ,4H    ,2/             0001560
C ....................................................................  0001561
       DATA TABLE(1 ,1  ) /500000.00   D0/                              0001562
       DATA TABLE(2 ,1  ) /309000.00   D0/                              0001563
       DATA TABLE(3 ,1  ) /1822.0      D0/                              0001564
       DATA TABLE(4 ,1  ) /21.00903    D0/                              0001565
       DATA TABLE(5 ,1  ) /.9999600000 D0/                              0001566
       DATA TABLE(6 ,1  ) /.3817065    D0/                              0001567
       DATA TABLE(1 ,2  ) /500000.00   D0/                              0001568
       DATA TABLE(2 ,2  ) /315000.00   D0/                              0001569
       DATA TABLE(3 ,2  ) /1792.0      D0/                              0001570
       DATA TABLE(4 ,2  ) /25.53386    D0/                              0001571
       DATA TABLE(5 ,2  ) /.9999333333 D0/                              0001572
       DATA TABLE(6 ,2  ) /.3817477    D0/                              0001573
       DATA TABLE(1 ,5  ) /500000.00   D0/                              0001574
       DATA TABLE(2 ,5  ) /396600.00   D0/                              0001575
       DATA TABLE(3 ,5  ) /1852.0      D0/                              0001576
       DATA TABLE(4 ,5  ) /16.62358    D0/                              0001577
       DATA TABLE(5 ,5  ) /.9999000000 D0/                              0001578
       DATA TABLE(6 ,5  ) /.3816485    D0/                              0001579
       DATA TABLE(1 ,6  ) /500000.00   D0/                              0001580
       DATA TABLE(2 ,6  ) /402900.00   D0/                              0001581
       DATA TABLE(3 ,6  ) /1852.0      D0/                              0001582
       DATA TABLE(4 ,6  ) /16.62358    D0/                              0001583
       DATA TABLE(5 ,6  ) /.9999000000 D0/                              0001584
       DATA TABLE(6 ,6  ) /.3816485    D0/                              0001585
       DATA TABLE(1 ,7  ) /500000.00   D0/                              0001586
       DATA TABLE(2 ,7  ) /409500.00   D0/                              0001587
       DATA TABLE(3 ,7  ) /1852.0      D0/                              0001588
       DATA TABLE(4 ,7  ) /16.62358    D0/                              0001589
       DATA TABLE(5 ,7  ) /.9999333333 D0/                              0001590
       DATA TABLE(6 ,7  ) /.3815948    D0/                              0001591
       DATA TABLE(1 ,21 ) /500000.00   D0/                              0001592
       DATA TABLE(2 ,21 ) /271500.00   D0/                              0001593
       DATA TABLE(3 ,21 ) /2271.0      D0/                              0001594
       DATA TABLE(4 ,21 ) /30.53702    D0/                              0001595
       DATA TABLE(5 ,21 ) /.9999950281 D0/                              0001596
       DATA TABLE(6 ,21 ) /.3811454    D0/                              0001597
       DATA TABLE(1 ,22 ) /500000.00   D0/                              0001598
       DATA TABLE(2 ,22 ) /291600.00   D0/                              0001599
       DATA TABLE(3 ,22 ) /1453.0      D0/                              0001600
       DATA TABLE(4 ,22 ) /26.09287    D0/                              0001601
       DATA TABLE(5 ,22 ) /.9999411765 D0/                              0001602
       DATA TABLE(6 ,22 ) /.3821090    D0/                              0001603
       DATA TABLE(1 ,23 ) /500000.00   D0/                              0001604
       DATA TABLE(2 ,23 ) /295200.00   D0/                              0001605
       DATA TABLE(3 ,23 ) /1453.0      D0/                              0001606
       DATA TABLE(4 ,23 ) /26.09287    D0/                              0001607
       DATA TABLE(5 ,23 ) /.9999411765 D0/                              0001608
       DATA TABLE(6 ,23 ) /.3821090    D0/                              0001609
       DATA TABLE(1 ,25 ) /500000.00   D0/                              0001610
       DATA TABLE(2 ,25 ) /295800.00   D0/                              0001611
       DATA TABLE(3 ,25 ) /1792.0      D0/                              0001612
       DATA TABLE(4 ,25 ) /25.53386    D0/                              0001613
       DATA TABLE(5 ,25 ) /.9999000000 D0/                              0001614
       DATA TABLE(6 ,25 ) /.3817593    D0/                              0001615
       DATA TABLE(1 ,26 ) /500000.00   D0/                              0001616
       DATA TABLE(2 ,26 ) /303000.00   D0/                              0001617
       DATA TABLE(3 ,26 ) /1792.0      D0/                              0001618
       DATA TABLE(4 ,26 ) /25.53386    D0/                              0001619
       DATA TABLE(5 ,26 ) /.9999000000 D0/                              0001620
       DATA TABLE(6 ,26 ) /.3817593    D0/                              0001621
       DATA TABLE(1 ,27 ) /500000.00   D0/                              0001622
       DATA TABLE(2 ,27 ) /559800.00   D0/                              0001623
       DATA TABLE(3 ,27 ) /1124.0      D0/                              0001624
       DATA TABLE(4 ,27 ) /39.52714    D0/                              0001625
       DATA TABLE(5 ,27 ) /.9999666667 D0/                              0001626
       DATA TABLE(6 ,27 ) /.3826496    D0/                              0001627
       DATA TABLE(1 ,28 ) /500000.00   D0/                              0001628
       DATA TABLE(2 ,28 ) /564000.00   D0/                              0001629
       DATA TABLE(3 ,28 ) /1214.0      D0/                              0001630
       DATA TABLE(4 ,28 ) /18.21554    D0/                              0001631
       DATA TABLE(5 ,28 ) /.9999666667 D0/                              0001632
       DATA TABLE(6 ,28 ) /.3825762    D0/                              0001633
       DATA TABLE(1 ,29 ) /500000.00   D0/                              0001634
       DATA TABLE(2 ,29 ) /568800.00   D0/                              0001635
       DATA TABLE(3 ,29 ) /1264.0      D0/                              0001636
       DATA TABLE(4 ,29 ) /6.77497     D0/                              0001637
       DATA TABLE(5 ,29 ) /.9999900000 D0/                              0001638
       DATA TABLE(6 ,29 ) /.3825176    D0/                              0001639
       DATA TABLE(1 ,30 ) /500000.00   D0/                              0001640
       DATA TABLE(2 ,30 ) /574200.00   D0/                              0001641
       DATA TABLE(3 ,30 ) /1303.0      D0/                              0001642
       DATA TABLE(4 ,30 ) /57.83623    D0/                              0001643
       DATA TABLE(5 ,30 ) /.9999900000 D0/                              0001644
       DATA TABLE(6 ,30 ) /.3824812    D0/                              0001645
       DATA TABLE(1 ,31 ) /500000.00   D0/                              0001646
       DATA TABLE(2 ,31 ) /576600.00   D0/                              0001647
       DATA TABLE(3 ,31 ) /1294.0      D0/                              0001648
       DATA TABLE(4 ,31 ) /0.05280     D0/                              0001649
       DATA TABLE(5 ,31 ) /.9999999999 D0/                              0001650
       DATA TABLE(6 ,31 ) /.3824867    D0/                              0001651
       DATA TABLE(1 ,32 ) /500000.00   D0/                              0001652
       DATA TABLE(2 ,32 ) /403800.00   D0/                              0001653
       DATA TABLE(3 ,32 ) /2491.0      D0/                              0001654
       DATA TABLE(4 ,32 ) /18.35156    D0/                              0001655
       DATA TABLE(5 ,32 ) /.9999473684 D0/                              0001656
       DATA TABLE(6 ,32 ) /.3807624    D0/                              0001657
       DATA TABLE(1 ,33 ) /500000.00   D0/                              0001658
       DATA TABLE(2 ,33 ) /410400.00   D0/                              0001659
       DATA TABLE(3 ,33 ) /2491.0      D0/                              0001660
       DATA TABLE(4 ,33 ) /18.35156    D0/                              0001661
       DATA TABLE(5 ,33 ) /.9999473684 D0/                              0001662
       DATA TABLE(6 ,33 ) /.3807624    D0/                              0001663
       DATA TABLE(1 ,34 ) /500000.00   D0/                              0001664
       DATA TABLE(2 ,34 ) /416700.00   D0/                              0001665
       DATA TABLE(3 ,34 ) /2491.0      D0/                              0001666
       DATA TABLE(4 ,34 ) /18.35156    D0/                              0001667
       DATA TABLE(5 ,34 ) /.9999333333 D0/                              0001668
       DATA TABLE(6 ,34 ) /.3806227    D0/                              0001669
       DATA TABLE(1 ,35 ) /500000.00   D0/                              0001670
       DATA TABLE(2 ,35 ) /318000.00   D0/                              0001671
       DATA TABLE(3 ,35 ) /2191.0      D0/                              0001672
       DATA TABLE(4 ,35 ) /37.04639    D0/                              0001673
       DATA TABLE(5 ,35 ) /.9999750000 D0/                              0001674
       DATA TABLE(6 ,35 ) /.3811074    D0/                              0001675
       DATA TABLE(1 ,36 ) /500000.00   D0/                              0001676
       DATA TABLE(2 ,36 ) /324600.00   D0/                              0001677
       DATA TABLE(3 ,36 ) /2191.0      D0/                              0001678
       DATA TABLE(4 ,36 ) /37.04639    D0/                              0001679
       DATA TABLE(5 ,36 ) /.9999411765 D0/                              0001680
       DATA TABLE(6 ,36 ) /.3811332    D0/                              0001681
       DATA TABLE(1 ,37 ) /500000.00   D0/                              0001682
       DATA TABLE(2 ,37 ) /308400.00   D0/                              0001683
       DATA TABLE(3 ,37 ) /2241.0      D0/                              0001684
       DATA TABLE(4 ,37 ) /32.84965    D0/                              0001685
       DATA TABLE(5 ,37 ) /.9999666667 D0/                              0001686
       DATA TABLE(6 ,37 ) /.3811064    D0/                              0001687
       DATA TABLE(1 ,38 ) /500000.00   D0/                              0001688
       DATA TABLE(2 ,38 ) /313500.00   D0/                              0001689
       DATA TABLE(3 ,38 ) /2241.0      D0/                              0001690
       DATA TABLE(4 ,38 ) /32.84965    D0/                              0001691
       DATA TABLE(5 ,38 ) /.9999666667 D0/                              0001692
       DATA TABLE(6 ,38 ) /.3811064    D0/                              0001693
       DATA TABLE(1 ,48 ) /500000.00   D0/                              0001694
       DATA TABLE(2 ,48 ) /246600.00   D0/                              0001695
       DATA TABLE(3 ,48 ) /2621.0      D0/                              0001696
       DATA TABLE(4 ,48 ) /15.15187    D0/                              0001697
       DATA TABLE(5 ,48 ) /.9999000000 D0/                              0001698
       DATA TABLE(6 ,48 ) /.3806180    D0/                              0001699
       DATA TABLE(1 ,49 ) /500000.00   D0/                              0001700
       DATA TABLE(2 ,49 ) /252600.00   D0/                              0001701
       DATA TABLE(3 ,49 ) /2561.0      D0/                              0001702
       DATA TABLE(4 ,49 ) /16.25668    D0/                              0001703
       DATA TABLE(5 ,49 ) /.9999666667 D0/                              0001704
       DATA TABLE(6 ,49 ) /.3806575    D0/                              0001705
       DATA TABLE(1 ,53 ) /500000.00   D0/                              0001706
       DATA TABLE(2 ,53 ) /301200.00   D0/                              0001707
       DATA TABLE(3 ,53 ) /2481.0      D0/                              0001708
       DATA TABLE(4 ,53 ) /18.72150    D0/                              0001709
       DATA TABLE(5 ,53 ) /.9999428571 D0/                              0001710
       DATA TABLE(6 ,53 ) /.3807283    D0/                              0001711
       DATA TABLE(1 ,54 ) /500000.00   D0/                              0001712
       DATA TABLE(2 ,54 ) /308700.00   D0/                              0001713
       DATA TABLE(3 ,54 ) /2481.0      D0/                              0001714
       DATA TABLE(4 ,54 ) /18.72150    D0/                              0001715
       DATA TABLE(5 ,54 ) /.9999090909 D0/                              0001716
       DATA TABLE(6 ,54 ) /.3807541    D0/                              0001717
       DATA TABLE(1 ,55 ) /500000.00   D0/                              0001718
       DATA TABLE(2 ,55 ) /319500.00   D0/                              0001719
       DATA TABLE(3 ,55 ) /2481.0      D0/                              0001720
       DATA TABLE(4 ,55 ) /18.72150    D0/                              0001721
       DATA TABLE(5 ,55 ) /.9999090909 D0/                              0001722
       DATA TABLE(6 ,55 ) /.3805361    D0/                              0001723
       DATA TABLE(1 ,62 ) /500000.00   D0/                              0001724
       DATA TABLE(2 ,62 ) /319800.00   D0/                              0001725
       DATA TABLE(3 ,62 ) /1772.0      D0/                              0001726
       DATA TABLE(4 ,62 ) /28.62716    D0/                              0001727
       DATA TABLE(5 ,62 ) /.9999600000 D0/                              0001728
       DATA TABLE(6 ,62 ) /.3817257    D0/                              0001729
       DATA TABLE(1 ,63 ) /500000.00   D0/                              0001730
       DATA TABLE(2 ,63 ) /325200.00   D0/                              0001731
       DATA TABLE(3 ,63 ) /1822.0      D0/                              0001732
       DATA TABLE(4 ,63 ) /21.00903    D0/                              0001733
       DATA TABLE(5 ,63 ) /.9999411765 D0/                              0001734
       DATA TABLE(6 ,63 ) /.3816986    D0/                              0001735
       DATA TABLE(1 ,64 ) /500000.00   D0/                              0001736
       DATA TABLE(2 ,64 ) /325800.00   D0/                              0001737
       DATA TABLE(3 ,64 ) /2141.0      D0/                              0001738
       DATA TABLE(4 ,64 ) /41.66790    D0/                              0001739
       DATA TABLE(5 ,64 ) /.9999333333 D0/                              0001740
       DATA TABLE(6 ,64 ) /.3812643    D0/                              0001741
       DATA TABLE(1 ,65 ) /500000.00   D0/                              0001742
       DATA TABLE(2 ,65 ) /333000.00   D0/                              0001743
       DATA TABLE(3 ,65 ) /2141.0      D0/                              0001744
       DATA TABLE(4 ,65 ) /41.66790    D0/                              0001745
       DATA TABLE(5 ,65 ) /.9999333333 D0/                              0001746
       DATA TABLE(6 ,65 ) /.3812422    D0/                              0001747
       DATA TABLE(1 ,66 ) /500000.00   D0/                              0001748
       DATA TABLE(2 ,66 ) /340200.00   D0/                              0001749
       DATA TABLE(3 ,66 ) /2161.0      D0/                              0001750
       DATA TABLE(4 ,66 ) /39.76857    D0/                              0001751
       DATA TABLE(5 ,66 ) /.9999411765 D0/                              0001752
       DATA TABLE(6 ,66 ) /.3812362    D0/                              0001753
       DATA TABLE(1 ,72 ) /500000.00   D0/                              0001754
       DATA TABLE(2 ,72 ) /416100.00   D0/                              0001755
       DATA TABLE(3 ,72 ) /2076.0      D0/                              0001756
       DATA TABLE(4 ,72 ) /48.30429    D0/                              0001757
       DATA TABLE(5 ,72 ) /.9999000000 D0/                              0001758
       DATA TABLE(6 ,72 ) /.3812311    D0/                              0001759
       DATA TABLE(1 ,73 ) /500000.00   D0/                              0001760
       DATA TABLE(2 ,73 ) /420000.00   D0/                              0001761
       DATA TABLE(3 ,73 ) /2076.0      D0/                              0001762
       DATA TABLE(4 ,73 ) /48.30429    D0/                              0001763
       DATA TABLE(5 ,73 ) /.9999000000 D0/                              0001764
       DATA TABLE(6 ,73 ) /.3812311    D0/                              0001765
       DATA TABLE(1 ,74 ) /500000.00   D0/                              0001766
       DATA TABLE(2 ,74 ) /426900.00   D0/                              0001767
       DATA TABLE(3 ,74 ) /2076.0      D0/                              0001768
       DATA TABLE(4 ,74 ) /48.30429    D0/                              0001769
       DATA TABLE(5 ,74 ) /.9999000000 D0/                              0001770
       DATA TABLE(6 ,74 ) /.3812311    D0/                              0001771
       DATA TABLE(1 ,75 ) /500000.00   D0/                              0001772
       DATA TABLE(2 ,75 ) /258000.00   D0/                              0001773
       DATA TABLE(3 ,75 ) /2541.0      D0/                              0001774
       DATA TABLE(4 ,75 ) /16.76677    D0/                              0001775
       DATA TABLE(5 ,75 ) /.9999666667 D0/                              0001776
       DATA TABLE(6 ,75 ) /.3807327    D0/                              0001777
       DATA TABLE(1 ,76 ) /2000000.00  D0/                              0001778
       DATA TABLE(2 ,76 ) /268800.00   D0/                              0001779
       DATA TABLE(3 ,76 ) /2321.0      D0/                              0001780
       DATA TABLE(4 ,76 ) /27.02745    D0/                              0001781
       DATA TABLE(5 ,76 ) /.9999750295 D0/                              0001782
       DATA TABLE(6 ,76 ) /.3810845    D0/                              0001783
       DATA TABLE(1 ,77 ) /500000.00   D0/                              0001784
       DATA TABLE(2 ,77 ) /375600.00   D0/                              0001785
       DATA TABLE(3 ,77 ) /1852.0      D0/                              0001786
       DATA TABLE(4 ,77 ) /16.62358    D0/                              0001787
       DATA TABLE(5 ,77 ) /.9999090909 D0/                              0001788
       DATA TABLE(6 ,77 ) /.3816135    D0/                              0001789
       DATA TABLE(1 ,78 ) /500000.00   D0/                              0001790
       DATA TABLE(2 ,78 ) /382500.00   D0/                              0001791
       DATA TABLE(3 ,78 ) /1852.0      D0/                              0001792
       DATA TABLE(4 ,78 ) /16.62358    D0/                              0001793
       DATA TABLE(5 ,78 ) /.9999000000 D0/                              0001794
       DATA TABLE(6 ,78 ) /.3816204    D0/                              0001795
       DATA TABLE(1 ,79 ) /500000.00   D0/                              0001796
       DATA TABLE(2 ,79 ) /388200.00   D0/                              0001797
       DATA TABLE(3 ,79 ) /1852.0      D0/                              0001798
       DATA TABLE(4 ,79 ) /16.62358    D0/                              0001799
       DATA TABLE(5 ,79 ) /.9999166667 D0/                              0001800
       DATA TABLE(6 ,79 ) /.3816288    D0/                              0001801
       DATA TABLE(1 ,80 ) /500000.00   D0/                              0001802
       DATA TABLE(2 ,80 ) /267600.00   D0/                              0001803
       DATA TABLE(3 ,80 ) /2391.0      D0/                              0001804
       DATA TABLE(4 ,80 ) /22.84247    D0/                              0001805
       DATA TABLE(5 ,80 ) /.9999666667 D0/                              0001806
       DATA TABLE(6 ,80 ) /.3808377    D0/                              0001807
       DATA TABLE(1 ,81 ) /500000.00   D0/                              0001808
       DATA TABLE(2 ,81 ) /275700.00   D0/                              0001809
       DATA TABLE(3 ,81 ) /2391.0      D0/                              0001810
       DATA TABLE(4 ,81 ) /22.84247    D0/                              0001811
       DATA TABLE(5 ,81 ) /.9999375000 D0/                              0001812
       DATA TABLE(6 ,81 ) /.3808450    D0/                              0001813
       DATA TABLE(1 ,82 ) /500000.00   D0/                              0001814
       DATA TABLE(2 ,82 ) /282900.00   D0/                              0001815
       DATA TABLE(3 ,82 ) /2391.0      D0/                              0001816
       DATA TABLE(4 ,82 ) /22.84247    D0/                              0001817
       DATA TABLE(5 ,82 ) /.9999375000 D0/                              0001818
       DATA TABLE(6 ,82 ) /.3808750    D0/                              0001819
       DATA TABLE(1 ,95 ) /500000.00   D0/                              0001820
       DATA TABLE(2 ,95 ) /257400.00   D0/                              0001821
       DATA TABLE(3 ,95 ) /2456.0      D0/                              0001822
       DATA TABLE(4 ,95 ) /19.72344    D0/                              0001823
       DATA TABLE(5 ,95 ) /.9999937500 D0/                              0001824
       DATA TABLE(6 ,95 ) /.3809220    D0/                              0001825
       DATA TABLE(1 ,109) /500000.00   D0/                              0001826
       DATA TABLE(2 ,109) /261000.00   D0/                              0001827
       DATA TABLE(3 ,109) /2541.0      D0/                              0001828
       DATA TABLE(4 ,109) /16.76677    D0/                              0001829
       DATA TABLE(5 ,109) /.9999642857 D0/                              0001830
       DATA TABLE(6 ,109) /.3807420    D0/                              0001831
       DATA TABLE(1 ,119) /500000.00   D0/                              0001832
       DATA TABLE(2 ,119) /378600.00   D0/                              0001833
       DATA TABLE(3 ,119) /2431.0      D0/                              0001834
       DATA TABLE(4 ,119) /20.83533    D0/                              0001835
       DATA TABLE(5 ,119) /.9999411765 D0/                              0001836
       DATA TABLE(6 ,119) /.3808422    D0/                              0001837
       DATA TABLE(1 ,120) /500000.00   D0/                              0001838
       DATA TABLE(2 ,120) /386400.00   D0/                              0001839
       DATA TABLE(3 ,120) /2431.0      D0/                              0001840
       DATA TABLE(4 ,120) /20.83533    D0/                              0001841
       DATA TABLE(5 ,120) /.9999411765 D0/                              0001842
       DATA TABLE(6 ,120) /.3808422    D0/                              0001843
       DATA TABLE(1 ,121) /500000.00   D0/                              0001844
       DATA TABLE(2 ,121) /391500.00   D0/                              0001845
       DATA TABLE(3 ,121) /2431.0      D0/                              0001846
       DATA TABLE(4 ,121) /20.83533    D0/                              0001847
       DATA TABLE(5 ,121) /.9999411765 D0/                              0001848
       DATA TABLE(6 ,121) /.3808422    D0/                              0001849
       DATA TABLE(1 ,122) /500000.00   D0/                              0001850
       DATA TABLE(2 ,122) /396300.00   D0/                              0001851
       DATA TABLE(3 ,122) /2431.0      D0/                              0001852
       DATA TABLE(4 ,122) /20.83533    D0/                              0001853
       DATA TABLE(5 ,122) /.9999411765 D0/                              0001854
       DATA TABLE(6 ,122) /.3808422    D0/                              0001855
       DATA TABLE(1 ,3  ) /3000000.00  D0/                              0001856
       DATA TABLE(2 ,3  ) /633600.00   D0/                              0001857
       DATA TABLE(3 ,3  ) /15893950.36 D0/                              0001858
       DATA TABLE(4 ,3  ) /16564628.77 D0/                              0001859
       DATA TABLE(5 ,3  ) /.9998480641 D0/                              0001860
       DATA TABLE(6 ,3  ) /.7969223940 D0/                              0001861
       DATA TABLE(7 ,3  ) /3161.0      D0/                              0001862
       DATA TABLE(8 ,3  ) /47.87068    D0/                              0001863
       DATA TABLE(9 ,3  ) /3.79919     D0/                              0001864
       DATA TABLE(10,3  ) /5.91550     D0/                              0001865
       DATA TABLE(11,3  ) /44.0        D0/                              0001866
       DATA TABLE(1 ,8  ) /2000000.00  D0/                              0001867
       DATA TABLE(2 ,8  ) /331200.00   D0/                              0001868
       DATA TABLE(3 ,8  ) /29277593.61 D0/                              0001869
       DATA TABLE(4 ,8  ) /29732882.87 D0/                              0001870
       DATA TABLE(5 ,8  ) /.9999359370 D0/                              0001871
       DATA TABLE(6 ,8  ) /.5818991407 D0/                              0001872
       DATA TABLE(7 ,8  ) /2126.0      D0/                              0001873
       DATA TABLE(8 ,8  ) /46.35656    D0/                              0001874
       DATA TABLE(9 ,8  ) /3.81452     D0/                              0001875
       DATA TABLE(10,8  ) /3.26432     D0/                              0001876
       DATA TABLE(11,8  ) /0.0         D0/                              0001877
       DATA TABLE(1 ,9  ) /2000000.00  D0/                              0001878
       DATA TABLE(2 ,9  ) /331200.00   D0/                              0001879
       DATA TABLE(3 ,9  ) /31014039.23 D0/                              0001880
       DATA TABLE(4 ,9  ) /31511724.20 D0/                              0001881
       DATA TABLE(5 ,9  ) /.9999184698 D0/                              0001882
       DATA TABLE(6 ,9  ) /.5596906871 D0/                              0001883
       DATA TABLE(7 ,9  ) /2033.0      D0/                              0001884
       DATA TABLE(8 ,9  ) /56.94711    D0/                              0001885
       DATA TABLE(9 ,9  ) /3.81550     D0/                              0001886
       DATA TABLE(10,9  ) /3.08256     D0/                              0001887
       DATA TABLE(11,9  ) /0.0         D0/                              0001888
       DATA TABLE(1 ,10 ) /2000000.00  D0/                              0001889
       DATA TABLE(2 ,10 ) /439200.00   D0/                              0001890
       DATA TABLE(3 ,10 ) /24245358.05 D0/                              0001891
       DATA TABLE(4 ,10 ) /24792436.23 D0/                              0001892
       DATA TABLE(5 ,10 ) /.9998946358 D0/                              0001893
       DATA TABLE(6 ,10 ) /.6538843192 D0/                              0001894
       DATA TABLE(7 ,10 ) /2441.0      D0/                              0001895
       DATA TABLE(8 ,10 ) /26.75847    D0/                              0001896
       DATA TABLE(9 ,10 ) /3.80992     D0/                              0001897
       DATA TABLE(10,10 ) /3.93575     D0/                              0001898
       DATA TABLE(11,10 ) /0.0         D0/                              0001899
       DATA TABLE(1 ,11 ) /2000000.00  D0/                              0001900
       DATA TABLE(2 ,11 ) /439200.00   D0/                              0001901
       DATA TABLE(3 ,11 ) /25795850.31 D0/                              0001902
       DATA TABLE(4 ,11 ) /26312257.65 D0/                              0001903
       DATA TABLE(5 ,11 ) /.9999146793 D0/                              0001904
       DATA TABLE(6 ,11 ) /.6304679732 D0/                              0001905
       DATA TABLE(7 ,11 ) /2336.0      D0/                              0001906
       DATA TABLE(8 ,11 ) /30.81964    D0/                              0001907
       DATA TABLE(9 ,11 ) /3.81147     D0/                              0001908
       DATA TABLE(10,11 ) /3.70114     D0/                              0001909
       DATA TABLE(11,11 ) /0.0         D0/                              0001910
       DATA TABLE(1 ,12 ) /2000000.00  D0/                              0001911
       DATA TABLE(2 ,12 ) /433800.00   D0/                              0001912
       DATA TABLE(3 ,12 ) /27057475.85 D0/                              0001913
       DATA TABLE(4 ,12 ) /27512992.04 D0/                              0001914
       DATA TABLE(5 ,12 ) /.9999291792 D0/                              0001915
       DATA TABLE(6 ,12 ) /.6122320427 D0/                              0001916
       DATA TABLE(7 ,12 ) /2256.0      D0/                              0001917
       DATA TABLE(8 ,12 ) /35.52018    D0/                              0001918
       DATA TABLE(9 ,12 ) /3.81265     D0/                              0001919
       DATA TABLE(10,12 ) /3.52998     D0/                              0001920
       DATA TABLE(11,12 ) /0.0         D0/                              0001921
       DATA TABLE(1 ,13 ) /2000000.00  D0/                              0001922
       DATA TABLE(2 ,13 ) /428400.00   D0/                              0001923
       DATA TABLE(3 ,13 ) /28182405.33 D0/                              0001924
       DATA TABLE(4 ,13 ) /28652931.96 D0/                              0001925
       DATA TABLE(5 ,13 ) /.9999407628 D0/                              0001926
       DATA TABLE(6 ,13 ) /.5965871443 D0/                              0001927
       DATA TABLE(7 ,13 ) /2189.0      D0/                              0001928
       DATA TABLE(8 ,13 ) /10.35494    D0/                              0001929
       DATA TABLE(9 ,13 ) /3.81362     D0/                              0001930
       DATA TABLE(10,13 ) /3.39020     D0/                              0001931
       DATA TABLE(11,13 ) /0.0         D0/                              0001932
       DATA TABLE(1 ,14 ) /2000000.00  D0/                              0001933
       DATA TABLE(2 ,14 ) /424800.00   D0/                              0001934
       DATA TABLE(3 ,14 ) /30194145.54 D0/                              0001935
       DATA TABLE(4 ,14 ) /30649424.27 D0/                              0001936
       DATA TABLE(5 ,14 ) /.9999221277 D0/                              0001937
       DATA TABLE(6 ,14 ) /.5700119219 D0/                              0001938
       DATA TABLE(7 ,14 ) /2076.0      D0/                              0001939
       DATA TABLE(8 ,14 ) /52.10305    D0/                              0001940
       DATA TABLE(9 ,14 ) /3.81523     D0/                              0001941
       DATA TABLE(10,14 ) /3.16593     D0/                              0001942
       DATA TABLE(11,14 ) /0.0         D0/                              0001943
       DATA TABLE(1 ,15 ) /2000000.00  D0/                              0001944
       DATA TABLE(2 ,15 ) /418500.00   D0/                              0001945
       DATA TABLE(3 ,15 ) /31846570.92 D0/                              0001946
       DATA TABLE(4 ,15 ) /32271267.72 D0/                              0001947
       DATA TABLE(5 ,15 ) /.9999541438 D0/                              0001948
       DATA TABLE(6 ,15 ) /.5495175982 D0/                              0001949
       DATA TABLE(7 ,15 ) /1992.0      D0/                              0001950
       DATA TABLE(8 ,15 ) /00.16335    D0/                              0001951
       DATA TABLE(9 ,15 ) /3.81642     D0/                              0001952
       DATA TABLE(10,15 ) /3.00292     D0/                              0001953
       DATA TABLE(11,15 ) /0.0         D0/                              0001954
       DATA TABLE(1 ,16 ) /4186692.58  D0/                              0001955
       DATA TABLE(2 ,16 ) /426000.00   D0/                              0001956
       DATA TABLE(3 ,16 ) /30891382.10 D0/                              0001957
       DATA TABLE(4 ,16 ) /35055396.31 D0/                              0001958
       DATA TABLE(5 ,16 ) /.9999885350 D0/                              0001959
       DATA TABLE(6 ,16 ) /.5612432071 D0/                              0001960
       DATA TABLE(7 ,16 ) /2040.0      D0/                              0001961
       DATA TABLE(8 ,16 ) /22.88096    D0/                              0001962
       DATA TABLE(9 ,16 ) /3.81572     D0/                              0001963
       DATA TABLE(10,16 ) /3.09520     D0/                              0001964
       DATA TABLE(11,16 ) /0.0         D0/                              0001965
       DATA TABLE(1 ,17 ) /2000000.00  D0/                              0001966
       DATA TABLE(2 ,17 ) /379800.00   D0/                              0001967
       DATA TABLE(3 ,17 ) /24751897.68 D0/                              0001968
       DATA TABLE(4 ,17 ) /25086068.20 D0/                              0001969
       DATA TABLE(5 ,17 ) /.9999568475 D0/                              0001970
       DATA TABLE(6 ,17 ) /.6461334829 D0/                              0001971
       DATA TABLE(7 ,17 ) /2406.0      D0/                              0001972
       DATA TABLE(8 ,17 ) /24.62308    D0/                              0001973
       DATA TABLE(9 ,17 ) /3.81044     D0/                              0001974
       DATA TABLE(10,17 ) /3.85610     D0/                              0001975
       DATA TABLE(11,17 ) /0.0         D0/                              0001976
       DATA TABLE(1 ,18 ) /2000000.00  D0/                              0001977
       DATA TABLE(2 ,18 ) /379800.00   D0/                              0001978
       DATA TABLE(3 ,18 ) /25781376.91 D0/                              0001979
       DATA TABLE(4 ,18 ) /26243052.74 D0/                              0001980
       DATA TABLE(5 ,18 ) /.9999359117 D0/                              0001981
       DATA TABLE(6 ,18 ) /.6306895773 D0/                              0001982
       DATA TABLE(7 ,18 ) /2337.0      D0/                              0001983
       DATA TABLE(8 ,18 ) /29.65162    D0/                              0001984
       DATA TABLE(9 ,18 ) /3.81146     D0/                              0001985
       DATA TABLE(10,18 ) /3.70326     D0/                              0001986
       DATA TABLE(11,18 ) /0.0         D0/                              0001987
       DATA TABLE(1 ,19 ) /2000000.00  D0/                              0001988
       DATA TABLE(2 ,19 ) /379800.00   D0/                              0001989
       DATA TABLE(3 ,19 ) /26977133.89 D0/                              0001990
       DATA TABLE(4 ,19 ) /27402231.82 D0/                              0001991
       DATA TABLE(5 ,19 ) /.9999453995 D0/                              0001992
       DATA TABLE(6 ,19 ) /.6133780528 D0/                              0001993
       DATA TABLE(7 ,19 ) /2261.0      D0/                              0001994
       DATA TABLE(8 ,19 ) /34.26662    D0/                              0001995
       DATA TABLE(9 ,19 ) /3.81257     D0/                              0001996
       DATA TABLE(10,19 ) /3.54046     D0/                              0001997
       DATA TABLE(11,19 ) /0.0         D0/                              0001998
       DATA TABLE(1 ,20 ) /600000.00   D0/                              0001999
       DATA TABLE(2 ,20 ) /261900.00   D0/                              0002000
       DATA TABLE(3 ,20 ) /23659233.56 D0/                              0002001
       DATA TABLE(4 ,20 ) /23914389.02 D0/                              0002002
       DATA TABLE(5 ,20 ) /.9999831405 D0/                              0002003
       DATA TABLE(6 ,20 ) /.6630594147 D0/                              0002004
       DATA TABLE(7 ,20 ) /2483.0      D0/                              0002005
       DATA TABLE(8 ,20 ) /19.67980    D0/                              0002006
       DATA TABLE(9 ,20 ) /3.80929     D0/                              0002007
       DATA TABLE(10,20 ) /4.03278     D0/                              0002008
       DATA TABLE(11,20 ) /0.0         D0/                              0002009
       DATA TABLE(1 ,24 ) /2000000.00  D0/                              0002010
       DATA TABLE(2 ,24 ) /304200.00   D0/                              0002011
       DATA TABLE(3 ,24 ) /36030443.05 D0/                              0002012
       DATA TABLE(4 ,24 ) /36454924.53 D0/                              0002013
       DATA TABLE(5 ,24 ) /.9999484343 D0/                              0002014
       DATA TABLE(6 ,24 ) /.5025259000 D0/                              0002015
       DATA TABLE(7 ,24 ) /1802.0      D0/                              0002016
       DATA TABLE(8 ,24 ) /26.11701    D0/                              0002017
       DATA TABLE(9 ,24 ) /3.81898     D0/                              0002018
       DATA TABLE(10,24 ) /2.65643     D0/                              0002019
       DATA TABLE(11,24 ) /0.0         D0/                              0002020
       DATA TABLE(1 ,39 ) /2000000.00  D0/                              0002021
       DATA TABLE(2 ,39 ) /336600.00   D0/                              0002022
       DATA TABLE(3 ,39 ) /22736950.34 D0/                              0002023
       DATA TABLE(4 ,39 ) /23162461.59 D0/                              0002024
       DATA TABLE(5 ,39 ) /.9999453686 D0/                              0002025
       DATA TABLE(6 ,39 ) /.6777445518 D0/                              0002026
       DATA TABLE(7 ,39 ) /2551.0      D0/                              0002027
       DATA TABLE(8 ,39 ) /20.02265    D0/                              0002028
       DATA TABLE(9 ,39 ) /3.80827     D0/                              0002029
       DATA TABLE(10,39 ) /4.19479     D0/                              0002030
       DATA TABLE(11,39 ) /0.0         D0/                              0002031
       DATA TABLE(1 ,40 ) /2000000.00  D0/                              0002032
       DATA TABLE(2 ,40 ) /336600.00   D0/                              0002033
       DATA TABLE(3 ,40 ) /23936585.11 D0/                              0002034
       DATA TABLE(4 ,40 ) /24374096.67 D0/                              0002035
       DATA TABLE(5 ,40 ) /.9999483705 D0/                              0002036
       DATA TABLE(6 ,40 ) /.6587010213 D0/                              0002037
       DATA TABLE(7 ,40 ) /2463.0      D0/                              0002038
       DATA TABLE(8 ,40 ) /22.59905    D0/                              0002039
       DATA TABLE(9 ,40 ) /3.80959     D0/                              0002040
       DATA TABLE(10,40 ) /3.98630     D0/                              0002041
       DATA TABLE(11,40 ) /0.0         D0/                              0002042
       DATA TABLE(1 ,41 ) /2000000.00  D0/                              0002043
       DATA TABLE(2 ,41 ) /352800.00   D0/                              0002044
       DATA TABLE(3 ,41 ) /25644959.12 D0/                              0002045
       DATA TABLE(4 ,41 ) /25979068.57 D0/                              0002046
       DATA TABLE(5 ,41 ) /.9999568556 D0/                              0002047
       DATA TABLE(6 ,41 ) /.6327148646 D0/                              0002048
       DATA TABLE(7 ,41 ) /2346.0      D0/                              0002049
       DATA TABLE(8 ,41 ) /27.97215    D0/                              0002050
       DATA TABLE(9 ,41 ) /3.81133     D0/                              0002051
       DATA TABLE(10,41 ) /3.72376     D0/                              0002052
       DATA TABLE(11,41 ) /0.0         D0/                              0002053
       DATA TABLE(1 ,42 ) /2000000.00  D0/                              0002054
       DATA TABLE(2 ,42 ) /354600.00   D0/                              0002055
       DATA TABLE(3 ,42 ) /26896024.48 D0/                              0002056
       DATA TABLE(4 ,42 ) /27351521.50 D0/                              0002057
       DATA TABLE(5 ,42 ) /.9999359200 D0/                              0002058
       DATA TABLE(6 ,42 ) /.6145281068 D0/                              0002059
       DATA TABLE(7 ,42 ) /2266.0      D0/                              0002060
       DATA TABLE(8 ,42 ) /34.41020    D0/                              0002061
       DATA TABLE(9 ,42 ) /3.81250     D0/                              0002062
       DATA TABLE(10,42 ) /3.55102     D0/                              0002063
       DATA TABLE(11,42 ) /0.0         D0/                              0002064
       DATA TABLE(1 ,43 ) /2000000.00  D0/                              0002065
       DATA TABLE(2 ,43 ) /303300.00   D0/                              0002066
       DATA TABLE(3 ,43 ) /26371820.68 D0/                              0002067
       DATA TABLE(4 ,43 ) /26724051.82 D0/                              0002068
       DATA TABLE(5 ,43 ) /.9999620817 D0/                              0002069
       DATA TABLE(6 ,43 ) /.6220672671 D0/                              0002070
       DATA TABLE(7 ,43 ) /2299.0      D0/                              0002071
       DATA TABLE(8 ,43 ) /30.63364    D0/                              0002072
       DATA TABLE(9 ,43 ) /3.81202     D0/                              0002073
       DATA TABLE(10,43 ) /3.62113     D0/                              0002074
       DATA TABLE(11,43 ) /0.0         D0/                              0002075
       DATA TABLE(1 ,44 ) /2000000.00  D0/                              0002076
       DATA TABLE(2 ,44 ) /308700.00   D0/                              0002077
       DATA TABLE(3 ,44 ) /27467860.75 D0/                              0002078
       DATA TABLE(4 ,44 ) /27832235.64 D0/                              0002079
       DATA TABLE(5 ,44 ) /.9999453808 D0/                              0002080
       DATA TABLE(6 ,44 ) /.6064623718 D0/                              0002081
       DATA TABLE(7 ,44 ) /2231.0      D0/                              0002082
       DATA TABLE(8 ,44 ) /36.57874    D0/                              0002083
       DATA TABLE(9 ,44 ) /3.81301     D0/                              0002084
       DATA TABLE(10,44 ) /3.47771     D0/                              0002085
       DATA TABLE(11,44 ) /0.0         D0/                              0002086
       DATA TABLE(1 ,45 ) /2000000.00  D0/                              0002087
       DATA TABLE(2 ,45 ) /333000.00   D0/                              0002088
       DATA TABLE(3 ,45 ) /33624568.36 D0/                              0002089
       DATA TABLE(4 ,45 ) /34079629.33 D0/                              0002090
       DATA TABLE(5 ,45 ) /.9999147417 D0/                              0002091
       DATA TABLE(6 ,45 ) /.5287006734 D0/                              0002092
       DATA TABLE(7 ,45 ) /1907.0      D0/                              0002093
       DATA TABLE(8 ,45 ) /12.68515    D0/                              0002094
       DATA TABLE(9 ,45 ) /3.81758     D0/                              0002095
       DATA TABLE(10,45 ) /2.84511     D0/                              0002096
       DATA TABLE(11,45 ) /0.0         D0/                              0002097
       DATA TABLE(1 ,46 ) /2000000.00  D0/                              0002098
       DATA TABLE(2 ,46 ) /328800.00   D0/                              0002099
       DATA TABLE(3 ,46 ) /36271389.35 D0/                              0002100
       DATA TABLE(4 ,46 ) /36756553.45 D0/                              0002101
       DATA TABLE(5 ,46 ) /.9999257458 D0/                              0002102
       DATA TABLE(6 ,46 ) /.5000126971 D0/                              0002103
       DATA TABLE(7 ,46 ) /1792.0      D0/                              0002104
       DATA TABLE(8 ,46 ) /28.55026    D0/                              0002105
       DATA TABLE(9 ,46 ) /3.81911     D0/                              0002106
       DATA TABLE(10,46 ) /2.63885     D0/                              0002107
       DATA TABLE(11,46 ) /0.0         D0/                              0002108
       DATA TABLE(1 ,47 ) /2000000.00  D0/                              0002109
       DATA TABLE(2 ,47 ) /328800.00   D0/                              0002110
       DATA TABLE(3 ,47 ) /41091749.54 D0/                              0002111
       DATA TABLE(4 ,47 ) /41576762.39 D0/                              0002112
       DATA TABLE(5 ,47 ) /.9998947956 D0/                              0002113
       DATA TABLE(6 ,47 ) /.4540068519 D0/                              0002114
       DATA TABLE(7 ,47 ) /1612.0      D0/                              0002115
       DATA TABLE(8 ,47 ) /59.30342    D0/                              0002116
       DATA TABLE(9 ,47 ) /3.82138     D0/                              0002117
       DATA TABLE(10,47 ) /2.27436     D0/                              0002118
       DATA TABLE(11,47 ) /25.0        D0/                              0002119
       DATA TABLE(1 ,50 ) /800000.00   D0/                              0002120
       DATA TABLE(2 ,50 ) /277200.00   D0/                              0002121
       DATA TABLE(4 ,50 ) /26369112.76 D0/                              0002122
       DATA TABLE(5 ,50 ) /.9999498485 D0/                              0002123
       DATA TABLE(6 ,50 ) /.6276341196 D0/                              0002124
       DATA TABLE(3 ,50 ) /25989474.99 D0/                              0002125
       DATA TABLE(7 ,50 ) /2323.0      D0/                              0002126
       DATA TABLE(8 ,50 ) /59.69369    D0/                              0002127
       DATA TABLE(9 ,50 ) /3.81166     D0/                              0002128
       DATA TABLE(10,50 ) /3.67392     D0/                              0002129
       DATA TABLE(11,50 ) /0.0         D0/                              0002130
       DATA TABLE(1 ,51 ) /600000.00   D0/                              0002131
       DATA TABLE(2 ,51 ) /257400.00   D0/                              0002132
       DATA TABLE(3 ,51 ) /23111975.14 D0/                              0002133
       DATA TABLE(4 ,51 ) /23549477.32 D0/                              0002134
       DATA TABLE(5 ,51 ) /.9999645506 D0/                              0002135
       DATA TABLE(6 ,51 ) /.6717286561 D0/                              0002136
       DATA TABLE(7 ,51 ) /2523.0      D0/                              0002137
       DATA TABLE(8 ,51 ) /19.53138    D0/                              0002138
       DATA TABLE(9 ,51 ) /3.80870     D0/                              0002139
       DATA TABLE(10,51 ) /4.12738     D0/                              0002140
       DATA TABLE(11,51 ) /0.0         D0/                              0002141
       DATA TABLE(1 ,52 ) /200000.00   D0/                              0002142
       DATA TABLE(2 ,52 ) /253800.00   D0/                              0002143
       DATA TABLE(3 ,52 ) /23784678.44 D0/                              0002144
       DATA TABLE(4 ,52 ) /23924398.02 D0/                              0002145
       DATA TABLE(5 ,52 ) /.9999984844 D0/                              0002146
       DATA TABLE(6 ,52 ) /.6610953994 D0/                              0002147
       DATA TABLE(7 ,52 ) /2474.0      D0/                              0002148
       DATA TABLE(8 ,52 ) /19.47463    D0/                              0002149
       DATA TABLE(9 ,52 ) /3.80943     D0/                              0002150
       DATA TABLE(10,52 ) /4.01174     D0/                              0002151
       DATA TABLE(11,52 ) /0.0         D0/                              0002152
       DATA TABLE(1 ,56 ) /2000000.00  D0/                              0002153
       DATA TABLE(2 ,56 ) /313200.00   D0/                              0002154
       DATA TABLE(3 ,56 ) /20041716.18 D0/                              0002155
       DATA TABLE(4 ,56 ) /20589420.09 D0/                              0002156
       DATA TABLE(5 ,56 ) /.9999410344 D0/                              0002157
       DATA TABLE(6 ,56 ) /.7227899381 D0/                              0002158
       DATA TABLE(7 ,56 ) /2768.0      D0/                              0002159
       DATA TABLE(8 ,56 ) /22.25085    D0/                              0002160
       DATA TABLE(9 ,56 ) /3.80501     D0/                              0002161
       DATA TABLE(10,56 ) /4.68430     D0/                              0002162
       DATA TABLE(11,56 ) /36.0        D0/                              0002163
       DATA TABLE(1 ,57 ) /2000000.00  D0/                              0002164
       DATA TABLE(2 ,57 ) /303600.00   D0/                              0002165
       DATA TABLE(3 ,57 ) /21001715.22 D0/                              0002166
       DATA TABLE(4 ,57 ) /21594768.40 D0/                              0002167
       DATA TABLE(5 ,57 ) /.9999509058 D0/                              0002168
       DATA TABLE(6 ,57 ) /.7064074100 D0/                              0002169
       DATA TABLE(7 ,57 ) /2687.0      D0/                              0002170
       DATA TABLE(8 ,57 ) /50.76661    D0/                              0002171
       DATA TABLE(9 ,57 ) /3.80622     D0/                              0002172
       DATA TABLE(10,57 ) /4.46875     D0/                              0002173
       DATA TABLE(11,57 ) /35.0        D0/                              0002174
       DATA TABLE(1 ,58 ) /2000000.00  D0/                              0002175
       DATA TABLE(2 ,58 ) /303600.00   D0/                              0002176
       DATA TABLE(3 ,58 ) /22564848.51 D0/                              0002177
       DATA TABLE(4 ,58 ) /23069597.22 D0/                              0002178
       DATA TABLE(5 ,58 ) /.9999450783 D0/                              0002179
       DATA TABLE(6 ,58 ) /.6805292633 D0/                              0002180
       DATA TABLE(7 ,58 ) /2564.0      D0/                              0002181
       DATA TABLE(8 ,58 ) /22.23938    D0/                              0002182
       DATA TABLE(9 ,58 ) /3.80808     D0/                              0002183
       DATA TABLE(10,58 ) /4.15706     D0/                              0002184
       DATA TABLE(11,58 ) /33.0        D0/                              0002185
       DATA TABLE(1 ,59 ) /2000000.00  D0/                              0002186
       DATA TABLE(2 ,59 ) /335160.00   D0/                              0002187
       DATA TABLE(3 ,59 ) /18984319.62 D0/                              0002188
       DATA TABLE(4 ,59 ) /19471398.75 D0/                              0002189
       DATA TABLE(5 ,59 ) /.9999028166 D0/                              0002190
       DATA TABLE(6 ,59 ) /.7412196637 D0/                              0002191
       DATA TABLE(7 ,59 ) /2861.0      D0/                              0002192
       DATA TABLE(8 ,59 ) /24.63011    D0/                              0002193
       DATA TABLE(9 ,59 ) /3.80362     D0/                              0002194
       DATA TABLE(10,59 ) /5.01609     D0/                              0002195
       DATA TABLE(11,59 ) /0.0         D0/                              0002196
       DATA TABLE(1 ,60 ) /2000000.00  D0/                              0002197
       DATA TABLE(2 ,60 ) /339300.00   D0/                              0002198
       DATA TABLE(3 ,60 ) /20006679.72 D0/                              0002199
       DATA TABLE(4 ,60 ) /20493457.15 D0/                              0002200
       DATA TABLE(5 ,60 ) /.9999220223 D0/                              0002201
       DATA TABLE(6 ,60 ) /.7233880702 D0/                              0002202
       DATA TABLE(7 ,60 ) /2771.0      D0/                              0002203
       DATA TABLE(8 ,60 ) /20.89747    D0/                              0002204
       DATA TABLE(9 ,60 ) /3.80497     D0/                              0002205
       DATA TABLE(10,60 ) /4.76197     D0/                              0002206
       DATA TABLE(11,60 ) /0.0         D0/                              0002207
       DATA TABLE(1 ,61 ) /2000000.00  D0/                              0002208
       DATA TABLE(2 ,61 ) /338400.00   D0/                              0002209
       DATA TABLE(3 ,61 ) /21327006.06 D0/                              0002210
       DATA TABLE(4 ,61 ) /21874349.14 D0/                              0002211
       DATA TABLE(5 ,61 ) /.9999220448 D0/                              0002212
       DATA TABLE(6 ,61 ) /.7009277824 D0/                              0002213
       DATA TABLE(7 ,61 ) /2661.0      D0/                              0002214
       DATA TABLE(8 ,61 ) /20.12517    D0/                              0002215
       DATA TABLE(9 ,61 ) /3.80662     D0/                              0002216
       DATA TABLE(10,61 ) /4.46959     D0/                              0002217
       DATA TABLE(11,61 ) /0.0         D0/                              0002218
       DATA TABLE(1 ,67 ) /2000000.00  D0/                              0002219
       DATA TABLE(2 ,67 ) /394200.00   D0/                              0002220
       DATA TABLE(3 ,67 ) /18689498.40 D0/                              0002221
       DATA TABLE(4 ,67 ) /19157874.26 D0/                              0002222
       DATA TABLE(5 ,67 ) /.9999714855 D0/                              0002223
       DATA TABLE(6 ,67 ) /.7464518080 D0/                              0002224
       DATA TABLE(7 ,67 ) /2888.0      D0/                              0002225
       DATA TABLE(8 ,67 ) /20.21285    D0/                              0002226
       DATA TABLE(9 ,67 ) /3.80322     D0/                              0002227
       DATA TABLE(10,67 ) /5.09490     D0/                              0002228
       DATA TABLE(11,67 ) /0.0         D0/                              0002229
       DATA TABLE(1 ,68 ) /2000000.00  D0/                              0002230
       DATA TABLE(2 ,68 ) /394200.00   D0/                              0002231
       DATA TABLE(3 ,68 ) /19432939.76 D0/                              0002232
       DATA TABLE(4 ,68 ) /19919806.36 D0/                              0002233
       DATA TABLE(5 ,68 ) /.9999220151 D0/                              0002234
       DATA TABLE(6 ,68 ) /.7333538278 D0/                              0002235
       DATA TABLE(7 ,68 ) /2821.0      D0/                              0002236
       DATA TABLE(8 ,68 ) /21.96779    D0/                              0002237
       DATA TABLE(9 ,68 ) /3.80422     D0/                              0002238
       DATA TABLE(10,68 ) /4.90135     D0/                              0002239
       DATA TABLE(11,68 ) /0.0         D0/                              0002240
       DATA TABLE(1 ,69 ) /2000000.00  D0/                              0002241
       DATA TABLE(2 ,69 ) /394200.00   D0/                              0002242
       DATA TABLE(3 ,69 ) /20500650.51 D0/                              0002243
       DATA TABLE(4 ,69 ) /21096820.93 D0/                              0002244
       DATA TABLE(5 ,69 ) /.9999107701 D0/                              0002245
       DATA TABLE(6 ,69 ) /.7149012442 D0/                              0002246
       DATA TABLE(7 ,69 ) /2729.0      D0/                              0002247
       DATA TABLE(8 ,69 ) /21.15820    D0/                              0002248
       DATA TABLE(9 ,69 ) /3.80560     D0/                              0002249
       DATA TABLE(10,69 ) /4.64814     D0/                              0002250
       DATA TABLE(11,69 ) /0.0         D0/                              0002251
       DATA TABLE(1 ,70 ) /2000000.00  D0/                              0002252
       DATA TABLE(2 ,70 ) /360000.00   D0/                              0002253
       DATA TABLE(3 ,70 ) /23004346.29 D0/                              0002254
       DATA TABLE(4 ,70 ) /23368977.46 D0/                              0002255
       DATA TABLE(5 ,70 ) /.9999645501 D0/                              0002256
       DATA TABLE(6 ,70 ) /.6734507906 D0/                              0002257
       DATA TABLE(7 ,70 ) /2531.0      D0/                              0002258
       DATA TABLE(8 ,70 ) /19.30504    D0/                              0002259
       DATA TABLE(9 ,70 ) /3.80858     D0/                              0002260
       DATA TABLE(10,70 ) /4.14653     D0/                              0002261
       DATA TABLE(11,70 ) /0.0         D0/                              0002262
       DATA TABLE(1 ,71 ) /2000000.00  D0/                              0002263
       DATA TABLE(2 ,71 ) /358200.00   D0/                              0002264
       DATA TABLE(3 ,71 ) /24104561.06 D0/                              0002265
       DATA TABLE(4 ,71 ) /24590781.86 D0/                              0002266
       DATA TABLE(5 ,71 ) /.9999220725 D0/                              0002267
       DATA TABLE(6 ,71 ) /.6560764003 D0/                              0002268
       DATA TABLE(7 ,71 ) /2451.0      D0/                              0002269
       DATA TABLE(8 ,71 ) /24.68139    D0/                              0002270
       DATA TABLE(9 ,71 ) /3.80977     D0/                              0002271
       DATA TABLE(10,71 ) /3.95865     D0/                              0002272
       DATA TABLE(11,71 ) /0.0         D0/                              0002273
       DATA TABLE(1 ,83 ) /2000000.00  D0/                              0002274
       DATA TABLE(2 ,83 ) /266400.00   D0/                              0002275
       DATA TABLE(3 ,83 ) /24235000.80 D0/                              0002276
       DATA TABLE(4 ,83 ) /24462545.30 D0/                              0002277
       DATA TABLE(5 ,83 ) /.9999949000 D0/                              0002278
       DATA TABLE(6 ,83 ) /.6540820950 D0/                              0002279
       DATA TABLE(7 ,83 ) /2442.0      D0/                              0002280
       DATA TABLE(8 ,83 ) /20.64240    D0/                              0002281
       DATA TABLE(9 ,83 ) /3.80990     D0/                              0002282
       DATA TABLE(10,83 ) /3.93780     D0/                              0002283
       DATA TABLE(11,83 ) /0.0         D0/                              0002284
       DATA TABLE(1 ,84 ) /2000000.00  D0/                              0002285
       DATA TABLE(2 ,84 ) /284400.00   D0/                              0002286
       DATA TABLE(3 ,84 ) /29637059.47 D0/                              0002287
       DATA TABLE(4 ,84 ) /30183611.25 D0/                              0002288
       DATA TABLE(5 ,84 ) /.9998725510 D0/                              0002289
       DATA TABLE(6 ,84 ) /.5771707700 D0/                              0002290
       DATA TABLE(7 ,84 ) /2106.0      D0/                              0002291
       DATA TABLE(8 ,84 ) /51.60353    D0/                              0002292
       DATA TABLE(9 ,84 ) /3.81480     D0/                              0002293
       DATA TABLE(10,84 ) /3.22483     D0/                              0002294
       DATA TABLE(11,84 ) /0.0         D0/                              0002295
       DATA TABLE(1 ,85 ) /2000000.00  D0/                              0002296
       DATA TABLE(2 ,85 ) /361800.00   D0/                              0002297
       DATA TABLE(3 ,85 ) /18819849.05 D0/                              0002298
       DATA TABLE(4 ,85 ) /19215516.01 D0/                              0002299
       DATA TABLE(5 ,85 ) /.9999358426 D0/                              0002300
       DATA TABLE(6 ,85 ) /.7441333961 D0/                              0002301
       DATA TABLE(7 ,85 ) /2876.0      D0/                              0002302
       DATA TABLE(8 ,85 ) /22.57950    D0/                              0002303
       DATA TABLE(9 ,85 ) /3.80339     D0/                              0002304
       DATA TABLE(10,85 ) /5.05972     D0/                              0002305
       DATA TABLE(11,85 ) /0.0         D0/                              0002306
       DATA TABLE(1 ,86 ) /2000000.00  D0/                              0002307
       DATA TABLE(2 ,86 ) /361800.00   D0/                              0002308
       DATA TABLE(3 ,86 ) /19661027.79 D0/                              0002309
       DATA TABLE(4 ,86 ) /20086977.18 D0/                              0002310
       DATA TABLE(5 ,86 ) /.9999358523 D0/                              0002311
       DATA TABLE(6 ,86 ) /.7293826040 D0/                              0002312
       DATA TABLE(7 ,86 ) /2801.0      D0/                              0002313
       DATA TABLE(8 ,86 ) /20.45445    D0/                              0002314
       DATA TABLE(9 ,86 ) /3.80452     D0/                              0002315
       DATA TABLE(10,86 ) /4.84504     D0/                              0002316
       DATA TABLE(11,86 ) /0.0         D0/                              0002317
       DATA TABLE(1 ,87 ) /2000000.00  D0/                              0002318
       DATA TABLE(2 ,87 ) /297000.00   D0/                              0002319
       DATA TABLE(3 ,87 ) /24048738.51 D0/                              0002320
       DATA TABLE(4 ,87 ) /24559158.47 D0/                              0002321
       DATA TABLE(5 ,87 ) /.9999391411 D0/                              0002322
       DATA TABLE(6 ,87 ) /.6569503193 D0/                              0002323
       DATA TABLE(7 ,87 ) /2455.0      D0/                              0002324
       DATA TABLE(8 ,87 ) /23.48125    D0/                              0002325
       DATA TABLE(9 ,87 ) /3.80971     D0/                              0002326
       DATA TABLE(10,87 ) /3.96783     D0/                              0002327
       DATA TABLE(11,87 ) /0.0         D0/                              0002328
       DATA TABLE(1 ,88 ) /2000000.00  D0/                              0002329
       DATA TABLE(2 ,88 ) /297000.00   D0/                              0002330
       DATA TABLE(3 ,88 ) /25522875.81 D0/                              0002331
       DATA TABLE(4 ,88 ) /26027071.12 D0/                              0002332
       DATA TABLE(5 ,88 ) /.9999359346 D0/                              0002333
       DATA TABLE(6 ,88 ) /.6345195439 D0/                              0002334
       DATA TABLE(7 ,88 ) /2354.0      D0/                              0002335
       DATA TABLE(8 ,88 ) /28.63705    D0/                              0002336
       DATA TABLE(9 ,88 ) /3.81121     D0/                              0002337
       DATA TABLE(10,88 ) /3.74048     D0/                              0002338
       DATA TABLE(11,88 ) /0.0         D0/                              0002339
       DATA TABLE(1 ,89 ) /2000000.00  D0/                              0002340
       DATA TABLE(2 ,89 ) /352800.00   D0/                              0002341
       DATA TABLE(3 ,89 ) /28657871.66 D0/                              0002342
       DATA TABLE(4 ,89 ) /29082831.70 D0/                              0002343
       DATA TABLE(5 ,89 ) /.9999454101 D0/                              0002344
       DATA TABLE(6 ,89 ) /.5901470744 D0/                              0002345
       DATA TABLE(7 ,89 ) /2161.0      D0/                              0002346
       DATA TABLE(8 ,89 ) /42.56887    D0/                              0002347
       DATA TABLE(9 ,89 ) /3.81402     D0/                              0002348
       DATA TABLE(10,89 ) /3.33440     D0/                              0002349
       DATA TABLE(11,89 ) /0.0         D0/                              0002350
       DATA TABLE(1 ,90 ) /2000000.00  D0/                              0002351
       DATA TABLE(2 ,90 ) /352800.00   D0/                              0002352
       DATA TABLE(3 ,90 ) /30382831.06 D0/                              0002353
       DATA TABLE(4 ,90 ) /30838032.96 D0/                              0002354
       DATA TABLE(5 ,90 ) /.9999359432 D0/                              0002355
       DATA TABLE(6 ,90 ) /.5676166827 D0/                              0002356
       DATA TABLE(7 ,90 ) /2066.0      D0/                              0002357
       DATA TABLE(8 ,90 ) /52.48935    D0/                              0002358
       DATA TABLE(9 ,90 ) /3.81537     D0/                              0002359
       DATA TABLE(10,90 ) /3.14645     D0/                              0002360
       DATA TABLE(11,90 ) /0.0         D0/                              0002361
       DATA TABLE(1 ,91 ) /2000000.00  D0/                              0002362
       DATA TABLE(2 ,91 ) /433800.00   D0/                              0002363
       DATA TABLE(3 ,91 ) /20836250.94 D0/                              0002364
       DATA TABLE(4 ,91 ) /21383852.48 D0/                              0002365
       DATA TABLE(5 ,91 ) /.9998945810 D0/                              0002366
       DATA TABLE(6 ,91 ) /.7091860222 D0/                              0002367
       DATA TABLE(7 ,91 ) /2701.0      D0/                              0002368
       DATA TABLE(8 ,91 ) /22.08858    D0/                              0002369
       DATA TABLE(9 ,91 ) /3.80602     D0/                              0002370
       DATA TABLE(10,91 ) /4.57382     D0/                              0002371
       DATA TABLE(11,91 ) /0.0         D0/                              0002372
       DATA TABLE(1 ,92 ) /2000000.00  D0/                              0002373
       DATA TABLE(2 ,92 ) /433800.00   D0/                              0002374
       DATA TABLE(3 ,92 ) /22341309.43 D0/                              0002375
       DATA TABLE(4 ,92 ) /22888667.15 D0/                              0002376
       DATA TABLE(5 ,92 ) /.9998946058 D0/                              0002377
       DATA TABLE(6 ,92 ) /.6841473833 D0/                              0002378
       DATA TABLE(7 ,92 ) /2581.0      D0/                              0002379
       DATA TABLE(8 ,92 ) /22.74104    D0/                              0002380
       DATA TABLE(9 ,92 ) /3.80782     D0/                              0002381
       DATA TABLE(10,92 ) /4.26823     D0/                              0002382
       DATA TABLE(11,92 ) /0.0         D0/                              0002383
       DATA TABLE(1 ,93 ) /2000000.00  D0/                              0002384
       DATA TABLE(2 ,93 ) /279900.00   D0/                              0002385
       DATA TABLE(3 ,93 ) /23755351.27 D0/                              0002386
       DATA TABLE(4 ,93 ) /24211050.37 D0/                              0002387
       DATA TABLE(5 ,93 ) /.9999568410 D0/                              0002388
       DATA TABLE(6 ,93 ) /.6615397363 D0/                              0002389
       DATA TABLE(7 ,93 ) /2476.0      D0/                              0002390
       DATA TABLE(8 ,93 ) /21.57953    D0/                              0002391
       DATA TABLE(9 ,93 ) /3.80940     D0/                              0002392
       DATA TABLE(10,93 ) /4.01753     D0/                              0002393
       DATA TABLE(11,93 ) /0.0         D0/                              0002394
       DATA TABLE(1 ,94 ) /2000000.00  D0/                              0002395
       DATA TABLE(2 ,94 ) /279900.00   D0/                              0002396
       DATA TABLE(3 ,94 ) /24577800.67 D0/                              0002397
       DATA TABLE(4 ,94 ) /24984826.43 D0/                              0002398
       DATA TABLE(5 ,94 ) /.9999595012 D0/                              0002399
       DATA TABLE(6 ,94 ) /.6487931668 D0/                              0002400
       DATA TABLE(7 ,94 ) /2418.0      D0/                              0002401
       DATA TABLE(8 ,94 ) /23.87979    D0/                              0002402
       DATA TABLE(9 ,94 ) /3.81026     D0/                              0002403
       DATA TABLE(10,94 ) /3.88319     D0/                              0002404
       DATA TABLE(11,94 ) /0.0         D0/                              0002405
       DATA TABLE(1 ,96 ) /2000000.00  D0/                              0002406
       DATA TABLE(2 ,96 ) /291600.00   D0/                              0002407
       DATA TABLE(3 ,96 ) /30630125.53 D0/                              0002408
       DATA TABLE(4 ,96 ) /31127724.75 D0/                              0002409
       DATA TABLE(5 ,96 ) /.9999454207 D0/                              0002410
       DATA TABLE(6 ,96 ) /.5644973800 D0/                              0002411
       DATA TABLE(7 ,96 ) /2053.0      D0/                              0002412
       DATA TABLE(8 ,96 ) /53.44099    D0/                              0002413
       DATA TABLE(9 ,96 ) /3.81555     D0/                              0002414
       DATA TABLE(10,96 ) /3.12127     D0/                              0002415
       DATA TABLE(11,96 ) /0.0         D0/                              0002416
       DATA TABLE(1 ,97 ) /2000000.00  D0/                              0002417
       DATA TABLE(2 ,97 ) /291600.00   D0/                              0002418
       DATA TABLE(3 ,97 ) /32252126.30 D0/                              0002419
       DATA TABLE(4 ,97 ) /32676887.65 D0/                              0002420
       DATA TABLE(5 ,97 ) /.9999326284 D0/                              0002421
       DATA TABLE(6 ,97 ) /.5446515700 D0/                              0002422
       DATA TABLE(7 ,97 ) /1972.0      D0/                              0002423
       DATA TABLE(8 ,97 ) /3.57839     D0/                              0002424
       DATA TABLE(9 ,97 ) /3.81669     D0/                              0002425
       DATA TABLE(10,97 ) /2.94381     D0/                              0002426
       DATA TABLE(11,97 ) /0.0         D0/                              0002427
       DATA TABLE(1 ,98 ) /2000000.00  D0/                              0002428
       DATA TABLE(2 ,98 ) /360000.00   D0/                              0002429
       DATA TABLE(3 ,98 ) /20922704.09 D0/                              0002430
       DATA TABLE(4 ,98 ) /21366697.03 D0/                              0002431
       DATA TABLE(5 ,98 ) /.9999391116 D0/                              0002432
       DATA TABLE(6 ,98 ) /.7077381841 D0/                              0002433
       DATA TABLE(7 ,98 ) /2694.0      D0/                              0002434
       DATA TABLE(8 ,98 ) /18.93392    D0/                              0002435
       DATA TABLE(9 ,98 ) /3.80912     D0/                              0002436
       DATA TABLE(10,98 ) /4.55529     D0/                              0002437
       DATA TABLE(11,98 ) /0.0         D0/                              0002438
       DATA TABLE(1 ,99 ) /2000000.00  D0/                              0002439
       DATA TABLE(2 ,99 ) /361200.00   D0/                              0002440
       DATA TABLE(3 ,99 ) /21993575.61 D0/                              0002441
       DATA TABLE(4 ,99 ) /22461937.05 D0/                              0002442
       DATA TABLE(5 ,99 ) /.9999068931 D0/                              0002443
       DATA TABLE(6 ,99 ) /.6898519579 D0/                              0002444
       DATA TABLE(7 ,99 ) /2608.0      D0/                              0002445
       DATA TABLE(8 ,99 ) /21.54370    D0/                              0002446
       DATA TABLE(9 ,99 ) /3.80742     D0/                              0002447
       DATA TABLE(10,99 ) /4.33519     D0/                              0002448
       DATA TABLE(11,99 ) /0.0         D0/                              0002449
       DATA TABLE(1 ,100) /2000000.00  D0/                              0002450
       DATA TABLE(2 ,100) /309600.00   D0/                              0002451
       DATA TABLE(3 ,100) /29010231.09 D0/                              0002452
       DATA TABLE(4 ,100) /29535149.91 D0/                              0002453
       DATA TABLE(5 ,100) /.9999484030 D0/                              0002454
       DATA TABLE(6 ,100) /.5854397296 D0/                              0002455
       DATA TABLE(7 ,100) /2141.0      D0/                              0002456
       DATA TABLE(8 ,100) /44.28313    D0/                              0002457
       DATA TABLE(9 ,100) /3.81431     D0/                              0002458
       DATA TABLE(10,100) /3.29422     D0/                              0002459
       DATA TABLE(11,100) /0.0         D0/                              0002460
       DATA TABLE(1 ,101) /2000000.00  D0/                              0002461
       DATA TABLE(2 ,101) /365400.00   D0/                              0002462
       DATA TABLE(3 ,101) /29456907.29 D0/                              0002463
       DATA TABLE(4 ,101) /29972959.94 D0/                              0002464
       DATA TABLE(5 ,101) /.9999108771 D0/                              0002465
       DATA TABLE(6 ,101) /.5795358654 D0/                              0002466
       DATA TABLE(7 ,101) /2116.0      D0/                              0002467
       DATA TABLE(8 ,101) /48.58548    D0/                              0002468
       DATA TABLE(9 ,101) /3.81466     D0/                              0002469
       DATA TABLE(10,101) /3.24452     D0/                              0002470
       DATA TABLE(11,101) /0.0         D0/                              0002471
       DATA TABLE(1 ,102) /2000000.00  D0/                              0002472
       DATA TABLE(2 ,102) /351000.00   D0/                              0002473
       DATA TABLE(3 ,102) /32187809.58 D0/                              0002474
       DATA TABLE(4 ,102) /32691654.54 D0/                              0002475
       DATA TABLE(5 ,102) /.9998726224 D0/                              0002476
       DATA TABLE(6 ,102) /.5453944146 D0/                              0002477
       DATA TABLE(7 ,102) /1975.0      D0/                              0002478
       DATA TABLE(8 ,102) /5.95074     D0/                              0002479
       DATA TABLE(9 ,102) /3.81665     D0/                              0002480
       DATA TABLE(10,102) /2.97107     D0/                              0002481
       DATA TABLE(11,102) /0.0         D0/                              0002482
       DATA TABLE(1 ,103) /2000000.00  D0/                              0002483
       DATA TABLE(2 ,103) /361200.00   D0/                              0002484
       DATA TABLE(3 ,103) /34851703.46 D0/                              0002485
       DATA TABLE(4 ,103) /35337121.23 D0/                              0002486
       DATA TABLE(5 ,103) /.9998817443 D0/                              0002487
       DATA TABLE(6 ,103) /.5150588857 D0/                              0002488
       DATA TABLE(7 ,103) /1852.0      D0/                              0002489
       DATA TABLE(8 ,103) /21.62181    D0/                              0002490
       DATA TABLE(9 ,103) /3.81832     D0/                              0002491
       DATA TABLE(10,103) /2.74550     D0/                              0002492
       DATA TABLE(11,103) /0.0         D0/                              0002493
       DATA TABLE(1 ,104) /2000000.00  D0/                              0002494
       DATA TABLE(2 ,104) /356400.00   D0/                              0002495
       DATA TABLE(3 ,104) /37261509.20 D0/                              0002496
       DATA TABLE(4 ,104) /37807440.38 D0/                              0002497
       DATA TABLE(5 ,104) /.9998632433 D0/                              0002498
       DATA TABLE(6 ,104) /.4899126408 D0/                              0002499
       DATA TABLE(7 ,104) /1752.0      D0/                              0002500
       DATA TABLE(8 ,104) /37.19059    D0/                              0002501
       DATA TABLE(9 ,104) /3.81962     D0/                              0002502
       DATA TABLE(10,104) /2.56899     D0/                              0002503
       DATA TABLE(11,104) /0.0         D0/                              0002504
       DATA TABLE(1 ,105) /2000000.00  D0/                              0002505
       DATA TABLE(2 ,105) /354600.00   D0/                              0002506
       DATA TABLE(3 ,105) /41091749.54 D0/                              0002507
       DATA TABLE(4 ,105) /41576762.39 D0/                              0002508
       DATA TABLE(5 ,105) /.9998947956 D0/                              0002509
       DATA TABLE(6 ,105) /.4540068519 D0/                              0002510
       DATA TABLE(7 ,105) /1612.0      D0/                              0002511
       DATA TABLE(8 ,105) /59.30342    D0/                              0002512
       DATA TABLE(9 ,105) /3.82138     D0/                              0002513
       DATA TABLE(10,105) /2.33094     D0/                              0002514
       DATA TABLE(11,105) /0.0         D0/                              0002515
       DATA TABLE(1 ,106) /2000000.00  D0/                              0002516
       DATA TABLE(2 ,106) /401400.00   D0/                              0002517
       DATA TABLE(3 ,106) /23894872.45 D0/                              0002518
       DATA TABLE(4 ,106) /24229110.29 D0/                              0002519
       DATA TABLE(5 ,106) /.9999568422 D0/                              0002520
       DATA TABLE(6 ,106) /.6593554910 D0/                              0002521
       DATA TABLE(7 ,106) /2466.0      D0/                              0002522
       DATA TABLE(8 ,106) /21.96231    D0/                              0002523
       DATA TABLE(9 ,106) /3.80955     D0/                              0002524
       DATA TABLE(10,106) /3.99323     D0/                              0002525
       DATA TABLE(11,106) /0.0         D0/                              0002526
       DATA TABLE(1 ,107) /2000000.00  D0/                              0002527
       DATA TABLE(2 ,107) /401400.00   D0/                              0002528
       DATA TABLE(3 ,107) /25117176.75 D0/                              0002529
       DATA TABLE(4 ,107) /25664114.42 D0/                              0002530
       DATA TABLE(5 ,107) /.9998988207 D0/                              0002531
       DATA TABLE(6 ,107) /.6405785926 D0/                              0002532
       DATA TABLE(7 ,107) /2381.0      D0/                              0002533
       DATA TABLE(8 ,107) /29.30066    D0/                              0002534
       DATA TABLE(9 ,107) /3.81081     D0/                              0002535
       DATA TABLE(10,107) /3.80024     D0/                              0002536
       DATA TABLE(11,107) /0.0         D0/                              0002537
       DATA TABLE(1 ,108) /2000000.00  D0/                              0002538
       DATA TABLE(2 ,108) /401400.00   D0/                              0002539
       DATA TABLE(3 ,108) /27025955.35 D0/                              0002540
       DATA TABLE(4 ,108) /27432812.88 D0/                              0002541
       DATA TABLE(5 ,108) /.9999512939 D0/                              0002542
       DATA TABLE(6 ,108) /.6126873424 D0/                              0002543
       DATA TABLE(7 ,108) /2258.0      D0/                              0002544
       DATA TABLE(8 ,108) /34.16878    D0/                              0002545
       DATA TABLE(9 ,108) /3.81262     D0/                              0002546
       DATA TABLE(10,108) /3.53414     D0/                              0002547
       DATA TABLE(11,108) /0.0         D0/                              0002548
       DATA TABLE(1 ,110) /2000000.00  D0/                              0002549
       DATA TABLE(2 ,110) /282600.00   D0/                              0002550
       DATA TABLE(3 ,110) /26230200.09 D0/                              0002551
       DATA TABLE(4 ,110) /26576444.45 D0/                              0002552
       DATA TABLE(5 ,110) /.9999483859 D0/                              0002553
       DATA TABLE(6 ,110) /.6241178597 D0/                              0002554
       DATA TABLE(7 ,110) /2308.0      D0/                              0002555
       DATA TABLE(8 ,110) /30.78682    D0/                              0002556
       DATA TABLE(9 ,110) /3.81189     D0/                              0002557
       DATA TABLE(10,110) /3.64047     D0/                              0002558
       DATA TABLE(11,110) /0.0         D0/                              0002559
       DATA TABLE(1 ,111) /2000000.00  D0/                              0002560
       DATA TABLE(2 ,111) /282600.00   D0/                              0002561
       DATA TABLE(3 ,111) /27434800.06 D0/                              0002562
       DATA TABLE(4 ,111) /27811312.71 D0/                              0002563
       DATA TABLE(5 ,111) /.9999454027 D0/                              0002564
       DATA TABLE(6 ,111) /.6069248249 D0/                              0002565
       DATA TABLE(7 ,111) /2233.0      D0/                              0002566
       DATA TABLE(8 ,111) /36.41072    D0/                              0002567
       DATA TABLE(9 ,111) /3.81298     D0/                              0002568
       DATA TABLE(10,111) /3.48187     D0/                              0002569
       DATA TABLE(11,111) /0.0         D0/                              0002570
       DATA TABLE(1 ,112) /2000000.00  D0/                              0002571
       DATA TABLE(2 ,112) /435000.00   D0/                              0002572
       DATA TABLE(3 ,112) /18798081.67 D0/                              0002573
       DATA TABLE(4 ,112) /19205863.43 D0/                              0002574
       DATA TABLE(5 ,112) /.9999422551 D0/                              0002575
       DATA TABLE(6 ,112) /.7445203390 D0/                              0002576
       DATA TABLE(7 ,112) /2878.0      D0/                              0002577
       DATA TABLE(8 ,112) /22.15711    D0/                              0002578
       DATA TABLE(9 ,112) /3.80336     D0/                              0002579
       DATA TABLE(10,112) /5.06556     D0/                              0002580
       DATA TABLE(11,112) /0.0         D0/                              0002581
       DATA TABLE(1 ,113) /2000000.00  D0/                              0002582
       DATA TABLE(2 ,113) /433800.00   D0/                              0002583
       DATA TABLE(3 ,113) /19832653.52 D0/                              0002584
       DATA TABLE(4 ,113) /20289119.60 D0/                              0002585
       DATA TABLE(5 ,113) /.9999145875 D0/                              0002586
       DATA TABLE(6 ,113) /.7263957947 D0/                              0002587
       DATA TABLE(7 ,113) /2786.0      D0/                              0002588
       DATA TABLE(8 ,113) /21.72121    D0/                              0002589
       DATA TABLE(9 ,113) /3.80474     D0/                              0002590
       DATA TABLE(10,113) /4.80336     D0/                              0002591
       DATA TABLE(11,113) /0.0         D0/                              0002592
       DATA TABLE(1 ,114) /2000000.00  D0/                              0002593
       DATA TABLE(2 ,114) /286200.00   D0/                              0002594
       DATA TABLE(3 ,114) /25305029.12 D0/                              0002595
       DATA TABLE(4 ,114) /25715126.55 D0/                              0002596
       DATA TABLE(5 ,114) /.9999407460 D0/                              0002597
       DATA TABLE(6 ,114) /.6377729696 D0/                              0002598
       DATA TABLE(7 ,114) /2368.0      D0/                              0002599
       DATA TABLE(8 ,114) /57.52979    D0/                              0002600
       DATA TABLE(9 ,114) /3.81099     D0/                              0002601
       DATA TABLE(10,114) /3.77244     D0/                              0002602
       DATA TABLE(11,114) /0.0         D0/                              0002603
       DATA TABLE(1 ,115) /2000000.00  D0/                              0002604
       DATA TABLE(2 ,115) /291600.00   D0/                              0002605
       DATA TABLE(3 ,115) /26639323.45 D0/                              0002606
       DATA TABLE(4 ,115) /27070620.78 D0/                              0002607
       DATA TABLE(5 ,115) /.9999256928 D0/                              0002608
       DATA TABLE(6 ,115) /.6181953936 D0/                              0002609
       DATA TABLE(7 ,115) /2282.0      D0/                              0002610
       DATA TABLE(8 ,115) /33.82207    D0/                              0002611
       DATA TABLE(9 ,115) /3.81227     D0/                              0002612
       DATA TABLE(10,115) /3.58491     D0/                              0002613
       DATA TABLE(11,115) /0.0         D0/                              0002614
       DATA TABLE(1 ,116) /2000000.00  D0/                              0002615
       DATA TABLE(2 ,116) /324000.00   D0/                              0002616
       DATA TABLE(3 ,116) /20124133.05 D0/                              0002617
       DATA TABLE(4 ,116) /20489179.67 D0/                              0002618
       DATA TABLE(5 ,116) /.9999453461 D0/                              0002619
       DATA TABLE(6 ,116) /.7213707913 D0/                              0002620
       DATA TABLE(7 ,116) /2761.0      D0/                              0002621
       DATA TABLE(8 ,116) /19.04034    D0/                              0002622
       DATA TABLE(9 ,116) /3.80511     D0/                              0002623
       DATA TABLE(10,116) /4.73451     D0/                              0002624
       DATA TABLE(11,116) /0.0         D0/                              0002625
       DATA TABLE(1 ,117) /2000000.00  D0/                              0002626
       DATA TABLE(2 ,117) /324000.00   D0/                              0002627
       DATA TABLE(3 ,117) /21050746.99 D0/                              0002628
       DATA TABLE(4 ,117) /21430913.91 D0/                              0002629
       DATA TABLE(5 ,117) /.9999407059 D0/                              0002630
       DATA TABLE(6 ,117) /.7055766312 D0/                              0002631
       DATA TABLE(7 ,117) /2683.0      D0/                              0002632
       DATA TABLE(8 ,117) /48.81363    D0/                              0002633
       DATA TABLE(9 ,117) /3.80628     D0/                              0002634
       DATA TABLE(10,117) /4.52782     D0/                              0002635
       DATA TABLE(11,117) /0.0         D0/                              0002636
       DATA TABLE(1 ,118) /2000000.00  D0/                              0002637
       DATA TABLE(2 ,118) /324000.00   D0/                              0002638
       DATA TABLE(3 ,118) /22161432.25 D0/                              0002639
       DATA TABLE(4 ,118) /22672134.66 D0/                              0002640
       DATA TABLE(5 ,118) /.9999325474 D0/                              0002641
       DATA TABLE(6 ,118) /.6871032423 D0/                              0002642
       DATA TABLE(7 ,118) /2595.0      D0/                              0002643
       DATA TABLE(8 ,118) /20.01691    D0/                              0002644
       DATA TABLE(9 ,118) /3.80761     D0/                              0002645
       DATA TABLE(10,118) /4.30274     D0/                              0002646
       DATA TABLE(11,118) /0.0         D0/                              0002647
C ....................................................................  0002648
      DATA C1,C2,C3,C4 /101.2794065D0,60.0D0,1052.893882D0,4.483344D0/  0002649
      DATA C5,C6,C7 /0.023520D0,100000000.0D0,0.009873675553D0/         0002650
      DATA C8,C9,C10 /1047.54671D0,6.19276D0,0.050912D0/                0002651
C ....................................................................  0002652
      DATA D1,D2,D3,D4 /3.9174D0,30.92241724D0,4.0831D0,3.280833333D0/  0002653
      DATA D5,D6,D7 /25.52381D0,0.3048006099D0,4.0831D0/                0002654
      DATA D8,D9,D10 /100000.0D0,0.0000000001D0,10000.0D0/              0002655
C ....................................................................  0002656
      DATA RADSEC /206264.806247D0/                                     0002657
C ....................................................................  0002658
      DATA ZERO,ONE,TWO /0.0D0,1.0D0,2.0D0/                             0002659
      DATA ESQ /0.006768658D0/                                          0002660
      DATA EPSLN,EPSLN1,EPSLN2 /0.0001D0,0.1D0,0.01D0/                  0002661
      DATA NIT /5/                                                      0002662
      DATA SWITCH /0/                                                   0002663
C                                                                       0002664
C ......................................................................0002665
C       .  INITIALIZATION OF PROJECTION PARAMETERS (FILE INPUT)  .      0002666
C ......................................................................0002667
C                                                                       0002668
      ENTRY IF02Z0 (INFILE,data)                                        0002669
C                                                                       0002670
      IERROR = 0                                                        0002671
      READ (INFILE,END=120) ZONE,BUFF                                   0002672
      IF (SWITCH.NE.0 .AND. SWITCH.EQ.ZONE) RETURN                      0002673
  020 IF (ZONE .LE. 0) GO TO 050                                        0002674
      DO 040 IND = 1,131                                                0002675
      IF (ZONE .EQ. ITEM(IND)) GO TO 060                                0002676
  040 CONTINUE                                                          0002677
  050 IF (IPEMSG .EQ. 0) PRINT 2000, ZONE                               0002678
 2000 FORMAT (' ERROR PJ02Z0'/                                          0002679
     .        ' ILLEGAL ZONE NO : ',I10)                                0002680
      IERROR = 021                                                      0002681
      RETURN                                                            0002682
  060 ITYPE = ID(9,IND) + 1                                             0002683
      GO TO (080,080,100) , ITYPE                                       0002684
  080 T1 = TABLE(1,IND)                                                 0002685
      T2 = TABLE(2,IND)                                                 0002686
      T3 = TABLE(3,IND)                                                 0002687
      T4 = TABLE(4,IND)                                                 0002688
      T5 = TABLE(5,IND)                                                 0002689
      T6 = TABLE(6,IND)                                                 0002690
      T7 = TABLE(7,IND)                                                 0002691
      T8 = TABLE(8,IND)                                                 0002692
      T9 = TABLE(9,IND)                                                 0002693
      T10= TABLE(10,IND)                                                0002694
      T11= TABLE(11,IND)                                                0002695
C                                                                       0002696
C LIST RESULTS OF PARAMETER INITIALIZATION.                             0002697
C                                                                       0002698
  100 IF (IPPARM .EQ. 0) PRINT 2010, (ID(I,IND),I=1,8)                  0002699
 2010 FORMAT (' INITIALIZATION PARAMETERS (STATE PLANE PROJECTION)'/    0002700
     .        ' ZONE = ',8A4)                                           0002701
      SWITCH = ZONE                                                     0002702
      RETURN                                                            0002703
  120 IF (IPEMSG .EQ. 0) PRINT 2020                                     0002704
 2020 FORMAT (' ERROR PJ02Z0'/                                          0002705
     .        ' MISSING PROJECTION PARAMETERS')                         0002706
      IERROR = 022                                                      0002707
      RETURN                                                            0002708
C                                                                       0002709
C ......................................................................0002710
C      .  INITIALIZATION OF PROJECTION PARAMETERS (ENTRY INPUT)  .      0002711
C ......................................................................0002712
C                                                                       0002713
      ENTRY IS02Z0 (ZZONE,DATA)                                         0002714
	zone = zzone
C                                                                       0002715
      IERROR = 0                                                        0002716
      IF (SWITCH.NE.0 .AND. SWITCH.EQ.ZONE) RETURN                      0002717
      BUFF(1) = DATA(1)                                                 0002718
      GO TO 020                                                         0002719
C                                                                       0002720
C ......................................................................0002721
C                      .  FORWARD TRANSFORMATION  .                     0002722
C ......................................................................0002723
C                                                                       0002724
      ENTRY PF02Z0 (GEOG,PROJ)                                          0002725
C                                                                       0002726
      IERROR = 0                                                        0002727
      IF (SWITCH .NE. 0) GO TO 140                                      0002728
      IF (IPEMSG .EQ. 0) PRINT 2020                                     0002729
      IERROR = 023                                                      0002730
      RETURN                                                            0002731
  140 GO TO (160,220,240) , ITYPE                                       0002732
C                                                                       0002733
C MERCATOR PROJECTION.                                                  0002734
C                                                                       0002735
  160 FLONG = DABS(GEOG(1) * RADSEC)                                    0002736
      FLAT = GEOG(2) * RADSEC                                           0002737
      DL = T2 - FLONG                                                   0002738
      DL1 = DL - D1 * (DL / D10)**3                                     0002739
      SINPH = DSIN(GEOG(2))                                             0002740
      COSPH = DCOS(GEOG(2))                                             0002741
      S1 = D2 * COSPH * DL1 / DSQRT(ONE - ESQ * SINPH * SINPH)          0002742
      SM = S1 + D3 * (S1 / D8)**3                                       0002743
      SG = D4 * T5 * SM                                                 0002744
      PROJ(1) = T1 + SG + T6 * (SG / D8)**3                             0002745
      DLPI = ZERO                                                       0002746
      DO 180 I = 1,NIT                                                  0002747
      PHI = GEOG(2) + DLPI / RADSEC                                     0002748
      SINPH = DSIN(PHI)                                                 0002749
      COSPH = DCOS(PHI)                                                 0002750
      B = D5 * D9 * SINPH / COSPH * (ONE - ESQ * SINPH * SINPH)**2      0002751
      DLPIN = B * SM**2                                                 0002752
      IF (DABS(DLPI - DLPIN) .LE. EPSLN) GO TO 200                      0002753
      DLPI = DLPIN                                                      0002754
  180 CONTINUE                                                          0002755
      IF (IPEMSG .EQ. 0) PRINT 2030, NIT                                0002756
 2030 FORMAT (' ERROR PJ02Z0'/                                          0002757
     .        ' FAILED TO COVERGE AFTER ',I2,' ITERATIONS')             0002758
      IERROR = 024                                                      0002759
      RETURN                                                            0002760
  200 FLATR = GEOG(2) + DLPIN / RADSEC                                  0002761
      SINPH = DSIN(FLATR)                                               0002762
      COSPH = DCOS(FLATR)                                               0002763
      PROJ(2) = C1 * T5 * (FLAT + DLPIN - C2 * T3 - T4 - SINPH * COSPH *0002764
     .           (C3 - COSPH * COSPH * (C4 - C5 * COSPH * COSPH)))      0002765
      RETURN                                                            0002766
C                                                                       0002767
C LAMBERT PROJECTION.                                                   0002768
C                                                                       0002769
  220 FLONG = DABS(GEOG(1) * RADSEC)                                    0002770
      FLAT = GEOG(2) * RADSEC                                           0002771
      SINPH = DSIN(GEOG(2))                                             0002772
      COSPH = DCOS(GEOG(2))                                             0002773
      TEMP = SINPH * COSPH * (C3 - COSPH * COSPH * (C4 - C5 * COSPH *   0002774
     .       COSPH))                                                    0002775
      S = C1 * (C2 * T7 + T8 - FLAT + TEMP)                             0002776
      TEMP = S / C6                                                     0002777
      TEMP1 = T9 - T10 * TEMP + T11 * TEMP * TEMP                       0002778
      R = T3 + S * T5 * (ONE + TEMP1 * TEMP * TEMP)                     0002779
      THETA = (T6 * (T2 - FLONG)) / RADSEC                              0002780
      SINTH = DSIN(THETA)                                               0002781
      SINHTH = DSIN(THETA / TWO)                                        0002782
      PROJ(1) = R * SINTH + T1                                          0002783
      PROJ(2) = T4 - R + TWO * R * SINHTH * SINHTH                      0002784
      RETURN                                                            0002785
C                                                                       0002786
C SPECIAL ALASKA PROJECTIONS.                                           0002787
C                                                                       0002788
  240 IZ = IND - 122                                                    0002789
      IFLG = 1                                                          0002790
      IF (IZ .NE. 1) GO TO 260                                          0002791
      CALL AL01Z0 (GEOG,PROJ,IFLG)                                      0002792
      RETURN                                                            0002793
  260 CALL AL29Z0 (GEOG,PROJ,IZ,IFLG)                                   0002794
      RETURN                                                            0002795
C                                                                       0002796
C ......................................................................0002797
C                      .  INVERSE TRANSFORMATION  .                     0002798
C ......................................................................0002799
C                                                                       0002800
      ENTRY PI02Z0 (PROJ,GEOG)                                          0002801
C                                                                       0002802
      IERROR = 0                                                        0002803
      IF (SWITCH .NE. 0) GO TO 300                                      0002804
      IF (IPEMSG .EQ. 0) PRINT 2020                                     0002805
      IERROR = 025                                                      0002806
      RETURN                                                            0002807
  300 GO TO (320,340,420) , ITYPE                                       0002808
C                                                                       0002809
C MERCATOR PROJECTION.                                                  0002810
C                                                                       0002811
  320 XP = PROJ(1) - T1                                                 0002812
      SG1 = XP - T6 * (XP / D8)**3                                      0002813
      SG = XP - T6 * (SG1 / D8)**3                                      0002814
      SM = D6 * SG / T5                                                 0002815
      WSEC = C2 * T3 + T4 + C7 * PROJ(2) / T5                           0002816
      W = WSEC / RADSEC                                                 0002817
      SINW = DSIN(W)                                                    0002818
      COSW = DCOS(W)                                                    0002819
      PHI = (WSEC + SINW * COSW * (C8 + COSW * COSW * (C9 + C10 *       0002820
     .      COSW * COSW))) / RADSEC                                     0002821
      SINPH = DSIN(PHI)                                                 0002822
      COSPH = DCOS(PHI)                                                 0002823
      GEOG(2) = PHI - (D5 * D9 * SM * SM * SINPH / COSPH *              0002824
     .           (ONE - ESQ * SINPH * SINPH)**2) / RADSEC               0002825
      SINPH = DSIN(GEOG(2))                                             0002826
      COSPH = DCOS(GEOG(2))                                             0002827
      SA = SM - D7 * (SM / D8)**3                                       0002828
      S1 = SM - D7 * (SA / D8)**3                                       0002829
      DLAM1 = S1 * DSQRT(ONE - ESQ * SINPH * SINPH) / (D2 * COSPH)      0002830
      DLAMA = DLAM1 + D1 * (DLAM1 / D10)**3                             0002831
      GEOG(1) =-(T2 - DLAM1 - D1 * (DLAMA / D10)**3) / RADSEC           0002832
      RETURN                                                            0002833
C                                                                       0002834
C LAMBERT PROJECTION.                                                   0002835
C                                                                       0002836
  340 TEMP = (PROJ(1) - T1) / (T4 - PROJ(2))                            0002837
      THETA = DATAN(TEMP)                                               0002838
      COSTH = DCOS(THETA)                                               0002839
      SINHTH = DSIN(THETA / TWO)                                        0002840
      GEOG(1) =-(T2 / RADSEC) + (THETA / T6)                            0002841
      R = (T4 - PROJ(2)) / COSTH                                        0002842
      YPRIME = PROJ(2) - TWO * R * SINHTH * SINHTH                      0002843
      SAPP = (T4 - T3 - YPRIME) / T5                                    0002844
      SAPPP = SAPP                                                      0002845
      DO 380 I = 1,NIT                                                  0002846
      TEMP = SAPP / C6                                                  0002847
      TEMP = T9 * TEMP**2 - T10 * TEMP**3 + T11 * TEMP**4               0002848
      S = SAPPP / (ONE + TEMP)                                          0002849
      IF (DABS(S - SAPP) .GT. EPSLN1) GO TO 360                         0002850
      S = SAPPP * (ONE - TEMP * (ONE - TEMP))                           0002851
      IF (DABS(S - SAPP) .LE. EPSLN2) GO TO 400                         0002852
  360 SAPP = S                                                          0002853
  380 CONTINUE                                                          0002854
      IF (IPEMSG .EQ. 0) PRINT 2030, NIT                                0002855
      IERROR = 026                                                      0002856
      RETURN                                                            0002857
  400 WSEC = C2 * T7 + T8 - C7 * S                                      0002858
      W = WSEC / RADSEC                                                 0002859
      SINW = DSIN(W)                                                    0002860
      COSW = DCOS(W)                                                    0002861
      GEOG(2) = (WSEC + SINW * COSW * (C8 + COSW * COSW * (C9 + C10 *   0002862
     .            COSW * COSW))) / RADSEC                               0002863
      RETURN                                                            0002864
C                                                                       0002865
C SPECIAL ALASKA PROJECTIONS.                                           0002866
C                                                                       0002867
  420 IZ = IND - 122                                                    0002868
      IFLG = 0                                                          0002869
      IF (IZ .NE. 1) GO TO 440                                          0002870
      CALL AL01Z0 (PROJ,GEOG,IFLG)                                      0002871
      RETURN                                                            0002872
  440 CALL AL29Z0 (PROJ,GEOG,IZ,IFLG)                                   0002873
      RETURN                                                            0002874
C                                                                       0002875
      END                                                               0002876

C                   PJ03Z0                                                     
C **********************************************************************0002878
C ** U.S.G.S. GENERAL MAP PROJECTION PACKAGE ...... DR. A. A. ELASSAL **0002879
C **          MATHEMATICAL ANALYSIS BY JOHN SNYDER                    **0002880
C ** MODULE I                VERSION 1.0.0            NOVEMBER 1,1980 **0002881
C **********************************************************************0002882
C                    *  ALBERS CONICAL EQUAL AREA  *                    0002883
C **********************************************************************0002884
C                                                                       0002885
      SUBROUTINE PJ03Z0                                                 0002886
C                                                                       0002887
      IMPLICIT REAL*8 (A-Z)                                             0002888
	integer*4 zzone
      INTEGER*4 IERROR,IPEMSG,IPPARM                                    0002889
      INTEGER*4 SWITCH,I,ZONE,ANGS,INFILE                               0002890
      COMMON /ELLPZ0/ AZ,EZ,ESZ,E0Z,E1Z,E2Z,E3Z                         0002891
C **** PARAMETERS **** A,E,ES,LAT1,LAT2,LON0,LAT0,X0,Y0,NS,C,RH0 *******0002892
      COMMON /ERRMZ0/ IERROR                                            0002893
      COMMON /PRINZ0/ IPEMSG,IPPARM                                     0002894
c     COMMON /WORKZ0/ BUFF(15),ANGS(4,4)                                0002895
	COMMON /WORKZ0/ BUFF(15)
	COMMON /WK03Z0/ ANGS(4,4)
      real*4 rangs1,rangs2,rangs3,rangs4
      equivalence (rangs1,angs(4,1))
      equivalence (rangs2,angs(4,2))
      equivalence (rangs3,angs(4,3))
      equivalence (rangs4,angs(4,4))
      DIMENSION DATA(1),GEOG(1),PROJ(1)                                 0002896
      DATA TOL,EPSLN /1.0D-7,1.0D-10/                                   0002897
      DATA HALFPI /1.5707963267948966D0/                                0002898
      DATA ZERO,HALF,ONE /0.0D0,0.5D0,1.0D0/                            0002899
      DATA SWITCH /0/                                                   0002900
C                                                                       0002901
C ......................................................................0002902
C       .  INITIALIZATION OF PROJECTION PARAMETERS (FILE INPUT)  .      0002903
C ......................................................................0002904
C                                                                       0002905
      ENTRY IF03Z0 (INFILE,DATA)                                        0002906
C                                                                       0002907
      IERROR = 0                                                        0002908
      READ (INFILE,END=160) ZONE,BUFF                                   0002909
      IF (SWITCH.NE.0 .AND. SWITCH.EQ.ZONE) RETURN                      0002910
  020 IF (BUFF(1) .LE. ZERO) GO TO 100                                  0002911
      A = BUFF(1)                                                       0002912
      B = BUFF(2)                                                       0002913
      IF (B .GT. ZERO) GO TO 040                                        0002914
      E = ZERO                                                          0002915
      ES = ZERO                                                         0002916
      GO TO 120                                                         0002917
  040 IF (B .GT. ONE) GO TO 060                                         0002918
      E = DSQRT (B)                                                     0002919
      ES = B                                                            0002920
      GO TO 120                                                         0002921
  060 ES = ONE - (B / A) ** 2                                           0002922
      E = DSQRT (ES)                                                    0002923
      GO TO 120                                                         0002924
  100 A = AZ                                                            0002925
      E = EZ                                                            0002926
      ES = ESZ                                                          0002927
  120 LAT1 = PAKRZ0 (BUFF(3))                                           0002928
      LAT2 = PAKRZ0 (BUFF(4))                                           0002929
      IF (DABS(LAT1+LAT2) .GE. EPSLN) GO TO 130                         0002930
      IF (IPEMSG .EQ. 0) PRINT 2000                                     0002931
 2000 FORMAT (' ERROR PJ03Z0'/                                          0002932
     .        ' EQUAL LATITUDES FOR ST. PARALLELS ON OPPOSITE',         0002933
     .        ' SIDES OF EQUATOR')                                      0002934
      IERROR = 031                                                      0002935
      RETURN                                                            0002936
  130 LON0 = PAKRZ0 (BUFF(5))                                           0002937
      LAT0 = PAKRZ0 (BUFF(6))                                           0002938
      X0 = BUFF(7)                                                      0002939
      Y0 = BUFF(8)                                                      0002940
      SINPHI = DSIN (LAT1)                                              0002941
      CON = SINPHI                                                      0002942
      COSPHI = DCOS (LAT1)                                              0002943
      MS1 = MSFNZ0 (E,SINPHI,COSPHI)                                    0002944
      QS1 = QSFNZ0 (E,SINPHI,COSPHI)                                    0002945
      SINPHI = DSIN (LAT2)                                              0002946
      COSPHI = DCOS (LAT2)                                              0002947
      MS2 = MSFNZ0 (E,SINPHI,COSPHI)                                    0002948
      QS2 = QSFNZ0 (E,SINPHI,COSPHI)                                    0002949
      SINPHI = DSIN (LAT0)                                              0002950
      COSPHI = DCOS (LAT0)                                              0002951
      QS0 = QSFNZ0 (E,SINPHI,COSPHI)                                    0002952
      IF (DABS(LAT1-LAT2) .GE. EPSLN) GO TO 140                         0002953
      NS = CON                                                          0002954
      GO TO 150                                                         0002955
  140 NS = (MS1 * MS1 - MS2 * MS2) / (QS2 - QS1)                        0002956
  150 C = MS1 * MS1 + NS * QS1                                          0002957
      RH0 = A * DSQRT (C - NS * QS0) / NS                               0002958
C                                                                       0002959
C LIST RESULTS OF PARAMETER INITIALIZATION.                             0002960
C                                                                       0002961
      CALL RADDZ0 (LAT1,ANGS(1,1))                                      0002962
      CALL RADDZ0 (LAT2,ANGS(1,2))                                      0002963
      CALL RADDZ0 (LON0,ANGS(1,3))                                      0002964
      CALL RADDZ0 (LAT0,ANGS(1,4))                                      0002965
c      IF (IPPARM .EQ. 0) PRINT 2010, A,ES,ANGS,X0,Y0                   0002966
      IF (IPPARM .EQ. 0) PRINT 2010, A,ES,angs(1,1),angs(2,1),
     .   angs(3,1),rangs1,angs(1,2),angs(2,2),angs(3,2),rangs2,
     .   angs(1,3),angs(2,3),angs(3,3),rangs3,angs(1,4),
     .   angs(2,4),angs(3,4),rangs4,X0,Y0
 2010 FORMAT (' INITIALIZATION PARAMETERS (ALBERS CONICAL EQUAL-AREA',  0002967
     .        ' PROJECTION)'/                                           0002968
     .        ' SEMI-MAJOR AXIS OF ELLIPSOID =',F12.2,' METERS'/        0002969
     .        ' ECCENTRICITY SQUARED         =',F12.9/                  0002970
     .        ' LATITUDE OF 1ST ST. PARALLEL = ',A1,2I3,F7.3/           0002971
     .        ' LATITUDE OF 2ND ST. PARALLEL = ',A1,2I3,F7.3/           0002972
     .        ' LONGITUDE OF ORIGIN          = ',A1,2I3,F7.3/           0002973
     .        ' LATITUDE OF ORIGIN           = ',A1,2I3,F7.3/           0002974
     .        ' FALSE EASTING                =',F12.2,' METERS'/        0002975
     .        ' FALSE NORTHING               =',F12.2,' METERS')        0002976
      DATA(1) = A                                                       0002977
      DATA(2) = ES                                                      0002978
      SWITCH = ZONE                                                     0002979
      RETURN                                                            0002980
  160 IF (IPEMSG .EQ. 0) PRINT 2020                                     0002981
 2020 FORMAT (' ERROR PJ03Z0'/                                          0002982
     .        ' MISSING PROJECTION PARAMETERS')                         0002983
      IERROR = 032                                                      0002984
      RETURN                                                            0002985
C                                                                       0002986
C ......................................................................0002987
C      .  INITIALIZATION OF PROJECTION PARAMETERS (ENTRY INPUT)  .      0002988
C ......................................................................0002989
C                                                                       0002990
      ENTRY IS03Z0 (ZZONE,DATA)                                         0002991
	zone = zzone
C                                                                       0002992
      IERROR = 0                                                        0002993
      IF (SWITCH.NE.0 .AND. SWITCH.EQ.ZONE) RETURN                      0002994
      DO 180 I = 1,8                                                    0002995
      BUFF(I) = DATA(I)                                                 0002996
  180 CONTINUE                                                          0002997
      GO TO 020                                                         0002998
C                                                                       0002999
C ......................................................................0003000
C                      .  FORWARD TRANSFORMATION  .                     0003001
C ......................................................................0003002
C                                                                       0003003
      ENTRY PF03Z0 (GEOG,PROJ)                                          0003004
C                                                                       0003005
      IERROR = 0                                                        0003006
      IF (SWITCH .NE. 0) GO TO 220                                      0003007
      IF (IPEMSG .EQ. 0) PRINT 2020                                     0003008
      IERROR = 033                                                      0003009
      RETURN                                                            0003010
  220 SINPHI = DSIN (GEOG(2))                                           0003011
      COSPHI = DCOS (GEOG(2))                                           0003012
      QS = QSFNZ0 (E,SINPHI,COSPHI)                                     0003013
      RH = A * DSQRT (C - NS * QS) / NS                                 0003014
      THETA = NS * ADJLZ0 (GEOG(1) - LON0)                              0003015
      PROJ(1) = X0 + RH * DSIN (THETA)                                  0003016
      PROJ(2) = Y0 + RH0 - RH * DCOS (THETA)                            0003017
      RETURN                                                            0003018
C                                                                       0003019
C ......................................................................0003020
C                      .  INVERSE TRANSFORMATION  .                     0003021
C ......................................................................0003022
C                                                                       0003023
      ENTRY PI03Z0 (PROJ,GEOG)                                          0003024
C                                                                       0003025
      IERROR = 0                                                        0003026
      IF (SWITCH .NE. 0) GO TO 240                                      0003027
      IF (IPEMSG .EQ. 0) PRINT 2020                                     0003028
      IERROR = 034                                                      0003029
      RETURN                                                            0003030
  240 X = PROJ(1) - X0                                                  0003031
      Y = RH0 - PROJ(2) + Y0                                            0003032
      RH = DSIGN (DSQRT (X * X + Y * Y) , NS)                           0003033
      THETA = ZERO                                                      0003034
      CON = DSIGN (ONE , NS)                                            0003035
      IF (RH .NE. ZERO) THETA = DATAN2 (CON * X , CON * Y)              0003036
      CON = RH * NS / A                                                 0003037
      QS = (C - CON * CON) / NS                                         0003038
      IF (E .LT. TOL) GO TO 260                                         0003039
      CON = ONE - HALF * (ONE - ES) * DLOG ((ONE - E) /                 0003040
     .      (ONE + E)) / E                                              0003041
      IF ((DABS(CON) - DABS(QS)) .GT. TOL) GO TO 260                    0003042
      GEOG(2) = DSIGN (HALFPI , QS)                                     0003043
      GO TO 280                                                         0003044
  260 GEOG(2) = PHI1Z0 (E,QS)                                           0003045
      IF (IERROR .EQ. 0) GO TO 280                                      0003046
      IERROR = 035                                                      0003047
      RETURN                                                            0003048
  280 GEOG(1) = ADJLZ0 (THETA / NS + LON0)                              0003049
      RETURN                                                            0003050
C                                                                       0003051
      END                                                               0003052

C                   PJ04Z0                                                     
C **********************************************************************0003054
C ** U.S.G.S. GENERAL MAP PROJECTION PACKAGE ...... DR. A. A. ELASSAL **0003055
C **          MATHEMATICAL ANALYSIS BY JOHN SNYDER                    **0003056
C ** MODULE I                VERSION 1.0.0            NOVEMBER 1,1980 **0003057
C **********************************************************************0003058
C                     *  LAMBERT CONFORMAL CONIC  *                     0003059
C **********************************************************************0003060
C                                                                       0003061
      SUBROUTINE PJ04Z0                                                 0003062
C                                                                       0003063
      IMPLICIT REAL*8 (A-Z)                                             0003064
	integer*4 zzone
      INTEGER*4 IERROR,IPEMSG,IPPARM                                    0003065
      INTEGER*4 SWITCH,I,ZONE,ANGS,INFILE                               0003066
      COMMON /ELLPZ0/ AZ,EZ,ESZ,E0Z,E1Z,E2Z,E3Z                         0003067
C **** PARAMETERS **** A,E,ES,LAT1,LAT2,LON0,LAT0,X0,Y0,NS,F,RH0 *******0003068
      COMMON /ERRMZ0/ IERROR                                            0003069
      COMMON /PRINZ0/ IPEMSG,IPPARM                                     0003070
c     COMMON /WORKZ0/ BUFF(15),ANGS(4,4)                                0003071
	COMMON /WORKZ0/ BUFF(15)
	COMMON /WK04Z0/ ANGS(4,4)
      real*4 rangs1,rangs2,rangs3,rangs4
      equivalence (rangs1,angs(4,1))
      equivalence (rangs2,angs(4,2))
      equivalence (rangs3,angs(4,3))
      equivalence (rangs4,angs(4,4))
      DIMENSION DATA(1),GEOG(1),PROJ(1)                                 0003072
      DATA HALFPI /1.5707963267948966D0/                                0003073
      DATA EPSLN /1.0D-10/                                              0003074
      DATA ZERO,ONE /0.0D0,1.0D0/                                       0003075
      DATA SWITCH /0/                                                   0003076
C                                                                       0003077
C ......................................................................0003078
C       .  INITIALIZATION OF PROJECTION PARAMETERS (FILE INPUT)  .      0003079
C ......................................................................0003080
C                                                                       0003081
      ENTRY IF04Z0 (INFILE,DATA)                                        0003082
C                                                                       0003083
      IERROR = 0                                                        0003084
      READ (INFILE,END=160) ZONE,BUFF                                   0003085
      IF (SWITCH.NE.0 .AND. SWITCH.EQ.ZONE) RETURN                      0003086
  020 IF (BUFF(1) .LE. ZERO) GO TO 100                                  0003087
      A = BUFF(1)                                                       0003088
      B = BUFF(2)                                                       0003089
      IF (B .GT. ZERO) GO TO 040                                        0003090
      E = ZERO                                                          0003091
      ES = ZERO                                                         0003092
      GO TO 120                                                         0003093
  040 IF (B .GT. ONE) GO TO 060                                         0003094
      E = DSQRT (B)                                                     0003095
      ES = B                                                            0003096
      GO TO 120                                                         0003097
  060 ES = ONE - (B / A) ** 2                                           0003098
      E = DSQRT (ES)                                                    0003099
      GO TO 120                                                         0003100
  100 A = AZ                                                            0003101
      E = EZ                                                            0003102
      ES = ESZ                                                          0003103
  120 LAT1 = PAKRZ0 (BUFF(3))                                           0003104
      LAT2 = PAKRZ0 (BUFF(4))                                           0003105
      IF (DABS(LAT1+LAT2) .GE. EPSLN) GO TO 130                         0003106
      IF (IPEMSG .EQ. 0) PRINT 2000                                     0003107
 2000 FORMAT (' ERROR PJ04Z0'/                                          0003108
     .        ' EQUAL LATITUDES FOR ST. PARALLELS ON OPPOSITE',         0003109
     .        ' SIDES OF EQUATOR')                                      0003110
      IERROR = 041                                                      0003111
      RETURN                                                            0003112
  130 LON0 = PAKRZ0 (BUFF(5))                                           0003113
      LAT0 = PAKRZ0 (BUFF(6))                                           0003114
      X0 = BUFF(7)                                                      0003115
      Y0 = BUFF(8)                                                      0003116
      SINPHI = DSIN (LAT1)                                              0003117
      CON = SINPHI                                                      0003118
      COSPHI = DCOS (LAT1)                                              0003119
      MS1 = MSFNZ0 (E,SINPHI,COSPHI)                                    0003120
      TS1 = TSFNZ0 (E,LAT1,SINPHI)                                      0003121
      SINPHI = DSIN (LAT2)                                              0003122
      COSPHI = DCOS (LAT2)                                              0003123
      MS2 = MSFNZ0 (E,SINPHI,COSPHI)                                    0003124
      TS2 = TSFNZ0 (E,LAT2,SINPHI)                                      0003125
      SINPHI = DSIN (LAT0)                                              0003126
      TS0 = TSFNZ0 (E,LAT0,SINPHI)                                      0003127
      IF (DABS(LAT1-LAT2) .GE. EPSLN) GO TO 140                         0003128
      NS = CON                                                          0003129
      GO TO 150                                                         0003130
  140 NS = DLOG (MS1 / MS2) / DLOG (TS1 / TS2)                          0003131
  150 F = MS1 / (NS * TS1 ** NS)                                        0003132
      RH0 = A * F * TS0 ** NS                                           0003133
C                                                                       0003134
C LIST RESULTS OF PARAMETER INITIALIZATION.                             0003135
C                                                                       0003136
      CALL RADDZ0 (LAT1,ANGS(1,1))                                      0003137
      CALL RADDZ0 (LAT2,ANGS(1,2))                                      0003138
      CALL RADDZ0 (LON0,ANGS(1,3))                                      0003139
      CALL RADDZ0 (LAT0,ANGS(1,4))                                      0003140
c     IF (IPPARM .EQ. 0) PRINT 2010, A,ES,ANGS,X0,Y0                    0003141
      IF (IPPARM .EQ. 0) PRINT 2010, A,ES,angs(1,1),angs(2,1),
     .   angs(3,1),rangs1,angs(1,2),angs(2,2),angs(3,2),rangs2,
     .   angs(1,3),angs(2,3),angs(3,3),rangs3,angs(1,4),
     .   angs(2,4),angs(3,4),rangs4,X0,Y0
 2010 FORMAT (' INITIALIZATION PARAMETERS (LAMBERT CONFORMAL CONIC',    0003142
     .        ' PROJECTION)'/                                           0003143
     .        ' SEMI-MAJOR AXIS OF ELLIPSOID =',F12.2,' METERS'/        0003144
     .        ' ECCENTRICITY SQUARED         =',F12.9/                  0003145
     .        ' LATITUDE OF 1ST ST. PARALLEL = ',A1,2I3,F7.3/           0003146
     .        ' LATITUDE OF 2ND ST. PARALLEL = ',A1,2I3,F7.3/           0003147
     .        ' LONGITUDE OF ORIGIN          = ',A1,2I3,F7.3/           0003148
     .        ' LATITUDE OF ORIGIN           = ',A1,2I3,F7.3/           0003149
     .        ' FALSE EASTING                =',F12.2,' METERS'/        0003150
     .        ' FALSE NORTHING               =',F12.2,' METERS')        0003151
      DATA(1) = A                                                       0003152
      DATA(2) = ES                                                      0003153
      SWITCH = ZONE                                                     0003154
      RETURN                                                            0003155
  160 IF (IPEMSG .EQ. 0) PRINT 2020                                     0003156
 2020 FORMAT (' ERROR PJ04Z0'/                                          0003157
     .        ' MISSING PROJECTION PARAMETERS')                         0003158
      IERROR = 042                                                      0003159
      RETURN                                                            0003160
C                                                                       0003161
C ......................................................................0003162
C      .  INITIALIZATION OF PROJECTION PARAMETERS (ENTRY INPUT)  .      0003163
C ......................................................................0003164
C                                                                       0003165
      ENTRY IS04Z0 (ZZONE,DATA)                                         0003166
	zone = zzone
C                                                                       0003167
      IERROR = 0                                                        0003168
      IF (SWITCH.NE.0 .AND. SWITCH.EQ.ZONE) RETURN                      0003169
      DO 180 I = 1,8                                                    0003170
      BUFF(I) = DATA(I)                                                 0003171
  180 CONTINUE                                                          0003172
      GO TO 020                                                         0003173
C                                                                       0003174
C ......................................................................0003175
C                      .  FORWARD TRANSFORMATION  .                     0003176
C ......................................................................0003177
C                                                                       0003178
      ENTRY PF04Z0 (GEOG,PROJ)                                          0003179
C                                                                       0003180
      IERROR = 0                                                        0003181
      IF (SWITCH .NE. 0) GO TO 200                                      0003182
      IF (IPEMSG .EQ. 0) PRINT 2020                                     0003183
      IERROR = 043                                                      0003184
      RETURN                                                            0003185
  200 CON = DABS (DABS (GEOG(2)) - HALFPI)                              0003186
      IF (CON .GT. EPSLN) GO TO 220                                     0003187
      CON = GEOG(2) * NS                                                0003188
      IF (CON .GT. ZERO) GO TO 210                                      0003189
      if (ipemsg .eq. 0) PRINT 2030
 2030 FORMAT (' ERROR PJ04Z0'/                                          0003191
     .        ' POINT CANNOT BE PROJECTED')                             0003192
      IERROR = 044                                                      0003193
      RETURN                                                            0003194
  210 RH = ZERO                                                         0003195
      GO TO 230                                                         0003196
  220 SINPHI = DSIN (GEOG(2))                                           0003197
      TS = TSFNZ0 (E,GEOG(2),SINPHI)                                    0003198
      RH = A * F * TS ** NS                                             0003199
  230 THETA = NS * ADJLZ0 (GEOG(1) - LON0)                              0003200
      PROJ(1) = X0 + RH * DSIN (THETA)                                  0003201
      PROJ(2) = Y0 + RH0 - RH * DCOS (THETA)                            0003202
      RETURN                                                            0003203
C                                                                       0003204
C ......................................................................0003205
C                      .  INVERSE TRANSFORMATION  .                     0003206
C ......................................................................0003207
C                                                                       0003208
      ENTRY PI04Z0 (PROJ,GEOG)                                          0003209
C                                                                       0003210
      IERROR = 0                                                        0003211
      IF (SWITCH .NE. 0) GO TO 240                                      0003212
      IF (IPEMSG .EQ. 0) PRINT 2020                                     0003213
      IERROR = 045                                                      0003214
      RETURN                                                            0003215
  240 X = PROJ(1) - X0                                                  0003216
      Y = RH0 - PROJ(2) + Y0                                            0003217
      RH = DSIGN (DSQRT (X*X + Y*Y) , NS)                               0003218
      THETA = ZERO                                                      0003219
      CON = DSIGN (ONE , NS)                                            0003220
      IF (RH .NE. ZERO) THETA = DATAN2 (CON * X , CON * Y)              0003221
      IF (RH.NE.ZERO .OR. NS.GT.ZERO) GO TO 250                         0003222
      GEOG(2) = - HALFPI                                                0003223
      GO TO 260                                                         0003224
  250 CON = ONE / NS                                                    0003225
      TS = (RH / (A * F)) ** CON                                        0003226
      GEOG(2) = PHI2Z0 (E,TS)                                           0003227
      IF (IERROR .EQ. 0) GO TO 260                                      0003228
      IERROR = 046                                                      0003229
      RETURN                                                            0003230
  260 GEOG(1) = ADJLZ0 (THETA / NS + LON0)                              0003231
      RETURN                                                            0003232
C                                                                       0003233
      END                                                               0003234

C                   PJ05Z0                                                     
C **********************************************************************0003236
C ** U.S.G.S. GENERAL MAP PROJECTION PACKAGE ...... DR. A. A. ELASSAL **0003237
C **          MATHEMATICAL ANALYSIS BY JOHN SNYDER                    **0003238
C ** MODULE I                VERSION 1.0.2       SEPTEMBER 23, 1983  ***0003239
C **********************************************************************0003240
C                            *  MERCATOR  *                             0003241
C **********************************************************************0003242
C                                                                       0003243
      SUBROUTINE PJ05Z0                                                 0003244
C                                                                       0003245
      IMPLICIT REAL*8 (A-Z)                                             0003246
	integer*4 zzone
      INTEGER*4 IERROR,IPEMSG,IPPARM                                    0003247
      INTEGER*4 SWITCH,I,ZONE,ANGS,INFILE                               0003248
      COMMON /ELLPZ0/ AZ,EZ,ESZ,E0Z,E1Z,E2Z,E3Z                         0003249
C **** PARAMETERS **** A,E,ES,LON0,X0,Y0,NS,F,RH0,LAT1,M1 **************0003250
      COMMON /ERRMZ0/ IERROR                                            0003251
      COMMON /PRINZ0/ IPEMSG,IPPARM                                     0003252
c     COMMON /WORKZ0/ BUFF(15),ANGS(4,2)                                0003253
	COMMON /WORKZ0/ BUFF(15)
	COMMON /WK05Z0/ ANGS(4,2)
      real*4 rangs1,rangs2
      equivalence (rangs1,angs(4,1))
      equivalence (rangs2,angs(4,2))
      DIMENSION DATA(1),GEOG(1),PROJ(1)                                 0003254
      DATA HALFPI /1.5707963267948966D0/                                0003255
      DATA EPSLN /1.0D-10/                                              0003256
      DATA ZERO,ONE /0.0D0,1.0D0/                                       0003257
      DATA SWITCH /0/                                                   0003258
C                                                                       0003259
C ......................................................................0003260
C       .  INITIALIZATION OF PROJECTION PARAMETERS (FILE INPUT)  .      0003261
C ......................................................................0003262
C                                                                       0003263
      ENTRY IF05Z0 (INFILE,DATA)                                        0003264
C                                                                       0003265
      IERROR = 0                                                        0003266
      READ (INFILE,END=160) ZONE,BUFF                                   0003267
      IF (SWITCH.NE.0 .AND. SWITCH.EQ.ZONE) RETURN                      0003268
  020 IF (BUFF(1) .LE. ZERO) GO TO 100                                  0003269
      A = BUFF(1)                                                       0003270
      B = BUFF(2)                                                       0003271
      IF (B .GT. ZERO) GO TO 040                                        0003272
      E = ZERO                                                          0003273
      ES = ZERO                                                         0003274
      GO TO 120                                                         0003275
  040 IF (B .GT. ONE) GO TO 060                                         0003276
      E = DSQRT (B)                                                     0003277
      ES = B                                                            0003278
      GO TO 120                                                         0003279
  060 ES = ONE - (B / A) ** 2                                           0003280
      E = DSQRT (ES)                                                    0003281
      GO TO 120                                                         0003282
  100 A = AZ                                                            0003283
      E = EZ                                                            0003284
      ES = ESZ                                                          0003285
  120 LON0 = PAKRZ0 (BUFF(5))                                           0003286
      LAT1 = PAKRZ0 (BUFF(6))                                           0003286
      M1 = DCOS(LAT1) / (DSQRT( ONE - ES * DSIN(LAT1) **2))             0003286
      X0 = BUFF(7)                                                      0003287
      Y0 = BUFF(8)                                                      0003288
C                                                                       0003289
C LIST RESULTS OF PARAMETER INITIALIZATION.                             0003290
C                                                                       0003291
      CALL RADDZ0 (LAT1,ANGS(1,1))                                      0003291
      CALL RADDZ0 (LON0,ANGS(1,2))                                      0003292
c     IF (IPPARM .EQ. 0) PRINT 2000, A,ES,ANGS,X0,Y0                    0003293
      IF (IPPARM .EQ. 0) PRINT 2000, A,ES,angs(1,1),angs(2,1),
     .   angs(3,1),rangs1,angs(1,2),angs(2,2),angs(3,2),rangs2,
     .   X0,Y0
 2000 FORMAT (' INITIALIZATION PARAMETERS (MERCATOR',                   0003294
     .        ' PROJECTION)'/                                           0003295
     .        ' SEMI-MAJOR AXIS OF ELLIPSOID =',F12.2,' METERS'/        0003296
     .        ' ECCENTRICITY SQUARED         =',F12.9/                  0003297
     .        ' LATITUDE OF TRUE SCALE       = ',A1,2I3,F7.3/           0003297
     .        ' CENTRAL LONGITUDE            = ',A1,2I3,F7.3/           0003298
     .        ' FALSE EASTING                =',F12.2,' METERS'/        0003299
     .        ' FALSE NORTHING               =',F12.2,' METERS')        0003300
      DATA(1) = A                                                       0003301
      DATA(2) = ES                                                      0003302
      SWITCH = ZONE                                                     0003303
      RETURN                                                            0003304
  160 IF (IPEMSG .EQ. 0) PRINT 2010                                     0003305
 2010 FORMAT (' ERROR PJ05Z0'/                                          0003306
     .        ' MISSING PROJECTION PARAMETERS')                         0003307
      IERROR = 051                                                      0003308
      RETURN                                                            0003309
C                                                                       0003310
C ......................................................................0003311
C      .  INITIALIZATION OF PROJECTION PARAMETERS (ENTRY INPUT)  .      0003312
C ......................................................................0003313
C                                                                       0003314
      ENTRY IS05Z0 (ZZONE,DATA)                                         0003315
	zone = zzone
C                                                                       0003316
      IERROR = 0                                                        0003317
      IF (SWITCH.NE.0 .AND. SWITCH.EQ.ZONE) RETURN                      0003318
      DO 180 I = 1,8                                                    0003319
      BUFF(I) = DATA(I)                                                 0003320
  180 CONTINUE                                                          0003321
      GO TO 020                                                         0003322
C                                                                       0003323
C ......................................................................0003324
C                      .  FORWARD TRANSFORMATION  .                     0003325
C ......................................................................0003326
C                                                                       0003327
      ENTRY PF05Z0 (GEOG,PROJ)                                          0003328
C                                                                       0003329
      IERROR = 0                                                        0003330
      IF (SWITCH .NE. 0) GO TO 220                                      0003331
      IF (IPEMSG .EQ. 0) PRINT 2010                                     0003332
      IERROR = 052                                                      0003333
      RETURN                                                            0003334
  220 IF (DABS(DABS(GEOG(2)) - HALFPI) .GT. EPSLN) GO TO 240            0003335
      IF (IPEMSG .EQ. 0) PRINT 2020                                     0003336
 2020 FORMAT (' ERROR PJ05Z0'/                                          0003337
     .        ' TRANSFORMATION CANNOT BE COMPUTED AT THE POLES')        0003338
      IERROR = 053                                                      0003339
      RETURN                                                            0003340
  240 SINPHI = DSIN (GEOG(2))                                           0003341
      TS = TSFNZ0 (E,GEOG(2),SINPHI)                                    0003342
      PROJ(1) = X0 + A * M1 * ADJLZ0 (GEOG(1) - LON0)                   0003343
      PROJ(2) = Y0 - A * M1 * DLOG (TS)                                 0003344
      RETURN                                                            0003345
C                                                                       0003346
C ......................................................................0003347
C                      .  INVERSE TRANSFORMATION  .                     0003348
C ......................................................................0003349
C                                                                       0003350
      ENTRY PI05Z0 (PROJ,GEOG)                                          0003351
C                                                                       0003352
      IERROR = 0                                                        0003353
      IF (SWITCH .NE. 0) GO TO 260                                      0003354
      IF (IPEMSG .EQ. 0) PRINT 2010                                     0003355
      IERROR = 054                                                      0003356
      RETURN                                                            0003357
  260 X = PROJ(1) - X0                                                  0003358
      Y = PROJ(2) - Y0                                                  0003359
      TS = DEXP (- Y / (A * M1))                                        0003360
      GEOG(2) = PHI2Z0 (E,TS)                                           0003361
      IF (IERROR .EQ. 0) GO TO 280                                      0003362
      IERROR = 055                                                      0003363
      RETURN                                                            0003364
  280 GEOG(1) = ADJLZ0 (LON0 + X / (A * M1))                            0003365
      RETURN                                                            0003366
C                                                                       0003367
      END                                                               0003368

C                   PJ06Z0                                                     
C **********************************************************************0003370
C ** U.S.G.S. GENERAL MAP PROJECTION PACKAGE ...... DR. A. A. ELASSAL **0003371
C **          MATHEMATICAL ANALYSIS BY JOHN SNYDER                    **0003372
C ** MODULE I                VERSION 1.0.0            NOVEMBER 1,1980 **0003373
C **********************************************************************0003374
C                       *  POLAR STEREOGRAPHIC  *                       0003375
C **********************************************************************0003376
C                                                                       0003377
      SUBROUTINE PJ06Z0                                                 0003378
C                                                                       0003379
      IMPLICIT REAL*8 (A-Z)                                             0003380
	integer*4 zzone
      INTEGER*4 IERROR,IPEMSG,IPPARM                                    0003381
      INTEGER*4 SWITCH,IND,I,ZONE,ANGS,INFILE                           0003382
      COMMON /ELLPZ0/ AZ,EZ,ESZ,E0Z,E1Z,E2Z,E3Z                         0003383
C **** PARAMETERS **** A,E,ES,LON0,LATC,X0,Y0,E3,MCS,TCS,FAC,IND *******0003384
      COMMON /ERRMZ0/ IERROR                                            0003385
      COMMON /PRINZ0/ IPEMSG,IPPARM                                     0003386
c     COMMON /WORKZ0/ BUFF(15),ANGS(4,2)                                0003387
	COMMON /WORKZ0/ BUFF(15)
	COMMON /WK06Z0/ ANGS(4,2)
      real*4 rangs1,rangs2
      equivalence (rangs1,angs(4,1))
      equivalence (rangs2,angs(4,2))
      DIMENSION DATA(1),GEOG(1),PROJ(1)                                 0003388
      DATA NINTYD /90000000.0D0/                                        0003389
      DATA ZERO,ONE,TWO /0.0D0,1.0D0,2.0D0/                             0003390
      DATA SWITCH /0/                                                   0003391
C                                                                       0003392
C ......................................................................0003393
C       .  INITIALIZATION OF PROJECTION PARAMETERS (FILE INPUT)  .      0003394
C ......................................................................0003395
C                                                                       0003396
      ENTRY IF06Z0 (INFILE,DATA)                                        0003397
C                                                                       0003398
      IERROR = 0                                                        0003399
      READ (INFILE,END=160) ZONE,BUFF                                   0003400
      IF (SWITCH.NE.0 .AND. SWITCH.EQ.ZONE) RETURN                      0003401
  020 IF (BUFF(1) .LE. ZERO) GO TO 100                                  0003402
      A = BUFF(1)                                                       0003403
      B = BUFF(2)                                                       0003404
      IF (B .GT. ZERO) GO TO 040                                        0003405
      E = ZERO                                                          0003406
      ES = ZERO                                                         0003407
      E3 = ONE                                                          0003408
      GO TO 120                                                         0003409
  040 IF (B .GT. ONE) GO TO 060                                         0003410
      E = DSQRT (B)                                                     0003411
      ES = B                                                            0003412
      GO TO 080                                                         0003413
  060 ES = ONE - (B / A) ** 2                                           0003414
      E = DSQRT (ES)                                                    0003415
  080 E3 = E3FNZ0 (E)                                                   0003416
      GO TO 120                                                         0003417
  100 A = AZ                                                            0003418
      E = EZ                                                            0003419
      ES = ESZ                                                          0003420
      E3 = E3Z                                                          0003421
  120 LON0 = PAKRZ0 (BUFF(5))                                           0003422
      SAVE = BUFF(6)                                                    0003423
      LATC = PAKRZ0 (SAVE)                                              0003424
      X0 = BUFF(7)                                                      0003425
      Y0 = BUFF(8)                                                      0003426
      FAC = ONE                                                         0003427
      IF (SAVE .LT. ZERO) FAC =-ONE                                     0003428
      IND = 0                                                           0003429
      IF (DABS(SAVE) .EQ. NINTYD) GO TO 130                             0003430
      IND = 1                                                           0003431
      CON1 = FAC * LATC                                                 0003432
      SINPHI = DSIN (CON1)                                              0003433
      COSPHI = DCOS (CON1)                                              0003434
      MCS = MSFNZ0 (E,SINPHI,COSPHI)                                    0003435
      TCS = TSFNZ0 (E,CON1,SINPHI)                                      0003436
C                                                                       0003437
C LIST RESULTS OF PARAMETER INITIALIZATION.                             0003438
C                                                                       0003439
  130 CALL RADDZ0 (LON0,ANGS(1,1))                                      0003440
      CALL RADDZ0 (LATC,ANGS(1,2))                                      0003441
c     IF (IPPARM .EQ. 0) PRINT 2000, A,ES,ANGS,X0,Y0                    0003442
      IF (IPPARM .EQ. 0) PRINT 2000, A,ES,angs(1,1),angs(2,1),
     .   angs(3,1),rangs1,angs(1,2),angs(2,2),angs(3,2),rangs2,
     .   X0,Y0
 2000 FORMAT (' INITIALIZATION PARAMETERS (POLAR STEREOGRAPHIC',        0003443
     .        ' PROJECTION)'/                                           0003444
     .        ' SEMI-MAJOR AXIS OF ELLIPSOID =',F12.2,' METERS'/        0003445
     .        ' ECCENTRICITY SQUARED         =',F12.9/                  0003446
     .        ' LONGITUDE OF Y-AXIS          = ',A1,2I3,F7.3/           0003447
     .        ' LATITUDE OF TRUE SCALE       = ',A1,2I3,F7.3/           0003448
     .        ' FALSE EASTING                =',F12.2,' METERS'/        0003449
     .        ' FALSE NORTHING               =',F12.2,' METERS')        0003450
      DATA(1) = A                                                       0003451
      DATA(2) = ES                                                      0003452
      SWITCH = ZONE                                                     0003453
      RETURN                                                            0003454
  160 IF (IPEMSG .EQ. 0) PRINT 2010                                     0003455
 2010 FORMAT (' ERROR PJ06Z0'/                                          0003456
     .        ' MISSING PROJECTION PARAMETERS')                         0003457
      IERROR = 061                                                      0003458
      RETURN                                                            0003459
C                                                                       0003460
C ......................................................................0003461
C      .  INITIALIZATION OF PROJECTION PARAMETERS (ENTRY INPUT)  .      0003462
C ......................................................................0003463
C                                                                       0003464
      ENTRY IS06Z0 (ZZONE,DATA)                                         0003465
	zone = zzone
C                                                                       0003466
      IERROR = 0                                                        0003467
      IF (SWITCH.NE.0 .AND. SWITCH.EQ.ZONE) RETURN                      0003468
      DO 180 I = 1,8                                                    0003469
      BUFF(I) = DATA(I)                                                 0003470
  180 CONTINUE                                                          0003471
      GO TO 020                                                         0003472
C                                                                       0003473
C ......................................................................0003474
C                      .  FORWARD TRANSFORMATION  .                     0003475
C ......................................................................0003476
C                                                                       0003477
      ENTRY PF06Z0 (GEOG,PROJ)                                          0003478
C                                                                       0003479
      IERROR = 0                                                        0003480
      IF (SWITCH .NE. 0) GO TO 220                                      0003481
      IF (IPEMSG .EQ. 0) PRINT 2010                                     0003482
      IERROR = 062                                                      0003483
      RETURN                                                            0003484
  220 CON1 = FAC * ADJLZ0 (GEOG(1) - LON0)                              0003485
      CON2 = FAC * GEOG(2)                                              0003486
      SINPHI = DSIN (CON2)                                              0003487
      TS = TSFNZ0 (E,CON2,SINPHI)                                       0003488
      IF (IND .EQ. 0) GO TO 240                                         0003489
      RH = A * MCS * TS / TCS                                           0003490
      GO TO 260                                                         0003491
  240 RH = TWO * A * TS / E3                                            0003492
  260 PROJ(1) = X0 + FAC * RH * DSIN (CON1)                             0003493
      PROJ(2) = Y0 - FAC * RH * DCOS (CON1)                             0003494
      RETURN                                                            0003495
C                                                                       0003496
C ......................................................................0003497
C                      .  INVERSE TRANSFORMATION  .                     0003498
C ......................................................................0003499
C                                                                       0003500
      ENTRY PI06Z0 (PROJ,GEOG)                                          0003501
C                                                                       0003502
      IERROR = 0                                                        0003503
      IF (SWITCH .NE. 0) GO TO 320                                      0003504
      IF (IPEMSG .EQ. 0) PRINT 2010                                     0003505
      IERROR = 063                                                      0003506
      RETURN                                                            0003507
  320 X = FAC * (PROJ(1) - X0)                                          0003508
      Y = FAC * (PROJ(2) - Y0)                                          0003509
      RH = DSQRT (X * X + Y * Y)                                        0003510
      IF (IND .EQ. 0) GO TO 340                                         0003511
      TS = RH * TCS / (A * MCS)                                         0003512
      GO TO 360                                                         0003513
  340 TS = RH * E3 / (TWO * A)                                          0003514
  360 GEOG(2) = FAC * PHI2Z0 (E,TS)                                     0003515
      IF (IERROR .EQ. 0) GO TO 380                                      0003516
      IERROR = 064                                                      0003517
      RETURN                                                            0003518
  380 IF (RH .NE. ZERO) GO TO 400                                       0003519
      GEOG(1) = FAC * LON0                                              0003520
      RETURN                                                            0003521
  400 GEOG(1) = ADJLZ0 (FAC * DATAN2 (X , -Y) + LON0)                   0003522
      RETURN                                                            0003523
C                                                                       0003524
      END                                                               0003525

C                   PJ07Z0                                                     
C **********************************************************************0003527
C ** U.S.G.S. GENERAL MAP PROJECTION PACKAGE ...... DR. A. A. ELASSAL **0003528
C **          MATHEMATICAL ANALYSIS BY JOHN SNYDER                    **0003529
C ** MODULE I                VERSION 1.0.0            NOVEMBER 1,1980 **0003530
C **********************************************************************0003531
C                            *  POLYCONIC  *                            0003532
C **********************************************************************0003533
C                                                                       0003534
      SUBROUTINE PJ07Z0                                                 0003535
C                                                                       0003536
      IMPLICIT REAL*8 (A-Z)                                             0003537
	integer*4 zzone
      INTEGER*4 IERROR,IPEMSG,IPPARM                                    0003538
      INTEGER*4 SWITCH,I,ZONE,ANGS,INFILE                               0003539
      COMMON /ELLPZ0/ AZ,EZ,ESZ,E0Z,E1Z,E2Z,E3Z                         0003540
C **** PARAMETERS **** A,E,ES,LON0,LAT0,X0,Y0,E0,E1,E2,ML0 *************0003541
      COMMON /ERRMZ0/ IERROR                                            0003542
      COMMON /PRINZ0/ IPEMSG,IPPARM                                     0003543
c     COMMON /WORKZ0/ BUFF(15),ANGS(4,2)                                0003544
	COMMON /WORKZ0/ BUFF(15)
	COMMON /WK07Z0/ ANGS(4,2)
      real*4 rangs1,rangs2
      equivalence (rangs1,angs(4,1))
      equivalence (rangs2,angs(4,2))
      DIMENSION DATA(1),GEOG(1),PROJ(1)                                 0003545
      DATA TOL /1.0D-7/                                                 0003546
      DATA ZERO,ONE /0.0D0,1.0D0/                                       0003547
      DATA SWITCH /0/                                                   0003548
C                                                                       0003549
C ......................................................................0003550
C       .  INITIALIZATION OF PROJECTION PARAMETERS (FILE INPUT)  .      0003551
C ......................................................................0003552
C                                                                       0003553
      ENTRY IF07Z0 (INFILE,DATA)                                        0003554
C                                                                       0003555
      IERROR = 0                                                        0003556
      READ (INFILE,END=160) ZONE,BUFF                                   0003557
      IF (SWITCH.NE.0 .AND. SWITCH.EQ.ZONE) RETURN                      0003558
  020 IF (BUFF(1) .LE. ZERO) GO TO 100                                  0003559
      A = BUFF(1)                                                       0003560
      B = BUFF(2)                                                       0003561
      IF (B .GT. ZERO) GO TO 040                                        0003562
      E = ZERO                                                          0003563
      ES = ZERO                                                         0003564
      E0 = ONE                                                          0003565
      E1 = ZERO                                                         0003566
      E2 = ZERO                                                         0003567
      GO TO 120                                                         0003568
  040 IF (B .GT. ONE) GO TO 060                                         0003569
      E = DSQRT (B)                                                     0003570
      ES = B                                                            0003571
      GO TO 080                                                         0003572
  060 ES = ONE - (B / A) ** 2                                           0003573
      E = DSQRT (ES)                                                    0003574
  080 E0 = E0FNZ0 (ES)                                                  0003575
      E1 = E1FNZ0 (ES)                                                  0003576
      E2 = E2FNZ0 (ES)                                                  0003577
      GO TO 120                                                         0003578
  100 A = AZ                                                            0003579
      E = EZ                                                            0003580
      ES = ESZ                                                          0003581
      E0 = E0Z                                                          0003582
      E1 = E1Z                                                          0003583
      E2 = E2Z                                                          0003584
  120 LON0 = PAKRZ0 (BUFF(5))                                           0003585
      LAT0 = PAKRZ0 (BUFF(6))                                           0003586
      X0 = BUFF(7)                                                      0003587
      Y0 = BUFF(8)                                                      0003588
      ML0 = MLFNZ0 (E0,E1,E2,LAT0)                                      0003589
C                                                                       0003590
C LIST RESULTS OF PARAMETER INITIALIZATION.                             0003591
C                                                                       0003592
      CALL RADDZ0 (LON0,ANGS(1,1))                                      0003593
      CALL RADDZ0 (LAT0,ANGS(1,2))                                      0003594
c     IF (IPPARM .EQ. 0) PRINT 2000, A,ES,ANGS,X0,Y0                    0003595
      IF (IPPARM .EQ. 0) PRINT 2000, A,ES,angs(1,1),angs(2,1),
     .   angs(3,1),rangs1,angs(1,2),angs(2,2),angs(3,2),rangs2,
     .   X0,Y0
 2000 FORMAT (' INITIALIZATION PARAMETERS (POLYCONIC',                  0003596
     .        ' PROJECTION)'/                                           0003597
     .        ' SEMI-MAJOR AXIS OF ELLIPSOID =',F12.2,' METERS'/        0003598
     .        ' ECCENTRICITY SQUARED         =',F12.9/                  0003599
     .        ' LONGITUDE OF ORIGIN          = ',A1,2I3,F7.3/           0003600
     .        ' LATITUDE OF ORIGIN           = ',A1,2I3,F7.3/           0003601
     .        ' FALSE EASTING                =',F12.2,' METERS'/        0003602
     .        ' FALSE NORTHING               =',F12.2,' METERS')        0003603
      DATA(1) = A                                                       0003604
      DATA(2) = ES                                                      0003605
      SWITCH = ZONE                                                     0003606
      RETURN                                                            0003607
  160 IF (IPEMSG .EQ. 0) PRINT 2010                                     0003608
 2010 FORMAT (' ERROR PJ07Z0'/                                          0003609
     .        ' MISSING PROJECTION PARAMETERS')                         0003610
      IERROR = 071                                                      0003611
      RETURN                                                            0003612
C                                                                       0003613
C ......................................................................0003614
C      .  INITIALIZATION OF PROJECTION PARAMETERS (ENTRY INPUT)  .      0003615
C ......................................................................0003616
C                                                                       0003617
      ENTRY IS07Z0 (ZZONE,DATA)                                         0003618
	zone = zzone
C                                                                       0003619
      IERROR = 0                                                        0003620
      IF (SWITCH.NE.0 .AND. SWITCH.EQ.ZONE) RETURN                      0003621
      DO 180 I = 1,8                                                    0003622
      BUFF(I) = DATA(I)                                                 0003623
  180 CONTINUE                                                          0003624
      GO TO 020                                                         0003625
C                                                                       0003626
C ......................................................................0003627
C                      .  FORWARD TRANSFORMATION  .                     0003628
C ......................................................................0003629
C                                                                       0003630
      ENTRY PF07Z0 (GEOG,PROJ)                                          0003631
C                                                                       0003632
      IERROR = 0                                                        0003633
      IF (SWITCH .NE. 0) GO TO 220                                      0003634
      IF (IPEMSG .EQ. 0) PRINT 2010                                     0003635
      IERROR = 072                                                      0003636
      RETURN                                                            0003637
  220 CON = ADJLZ0 (GEOG(1) - LON0)                                     0003638
      IF (DABS(GEOG(2)) .GT. TOL) GO TO 240                             0003639
      PROJ(1) = X0 + A * CON                                            0003640
      PROJ(2) = Y0 - A * ML0                                            0003641
      RETURN                                                            0003642
  240 SINPHI = DSIN (GEOG(2))                                           0003643
      COSPHI = DCOS (GEOG(2))                                           0003644
      ML = MLFNZ0 (E0,E1,E2,GEOG(2))                                    0003645
      MS = MSFNZ0 (E,SINPHI,COSPHI)                                     0003646
      CON = CON * SINPHI                                                0003647
      PROJ(1) = X0 + A * MS * DSIN (CON) / SINPHI                       0003648
      PROJ(2) = Y0 + A * (ML - ML0 + MS * (ONE - DCOS (CON)) / SINPHI)  0003649
      RETURN                                                            0003650
C                                                                       0003651
C ......................................................................0003652
C                      .  INVERSE TRANSFORMATION  .                     0003653
C ......................................................................0003654
C                                                                       0003655
      ENTRY PI07Z0 (PROJ,GEOG)                                          0003656
C                                                                       0003657
      IERROR = 0                                                        0003658
      IF (SWITCH .NE. 0) GO TO 320                                      0003659
      IF (IPEMSG .EQ. 0) PRINT 2010                                     0003660
      IERROR = 073                                                      0003661
      RETURN                                                            0003662
  320 X = PROJ(1) - X0                                                  0003663
      Y = PROJ(2) - Y0                                                  0003664
      AL = ML0 + Y / A                                                  0003665
      IF (DABS (AL) .GT. TOL) GO TO 340                                 0003666
      GEOG(1) = X / A + LON0                                            0003667
      GEOG(2) = ZERO                                                    0003668
      RETURN                                                            0003669
  340 B = AL * AL + (X / A) ** 2                                        0003670
      GEOG(2) = PHI4Z0 (ES,E0,E1,E2,AL,B,C)                             0003671
      IF (IERROR .EQ. 0) GO TO 360                                      0003672
      IERROR = 074                                                      0003673
      RETURN                                                            0003674
  360 GEOG(1) = ADJLZ0 ( DASIN (X * C / A) / DSIN (GEOG(2)) + LON0)     0003675
      RETURN                                                            0003676
C                                                                       0003677
      END                                                               0003678

C                   PJ08Z0                                                     
C **********************************************************************0003680
C ** U.S.G.S. GENERAL MAP PROJECTION PACKAGE ...... DR. A. A. ELASSAL **0003681
C **          MATHEMATICAL ANALYSIS BY JOHN SNYDER                    **0003682
C ** MODULE I                VERSION 1.0.0            NOVEMBER 1,1980 **0003683
C **********************************************************************0003684
C                        *  EQUIDISTANT CONIC  *                        0003685
C **********************************************************************0003686
C                                                                       0003687
      SUBROUTINE PJ08Z0                                                 0003688
C                                                                       0003689
      IMPLICIT REAL*8 (A-Z)                                             0003690
	integer*4 zzone
      INTEGER*4 IERROR,IPEMSG,IPPARM                                    0003691
      INTEGER*4 SWITCH,IND,I,ZONE,ANGS,ANG1,ANG2,INFILE                 0003692
      COMMON /ELLPZ0/ AZ,EZ,ESZ,E0Z,E1Z,E2Z,E3Z                         0003693
C ** PARAMETERS ** A,E,ES,LAT1,LAT2,LON0,LAT0,X0,Y0,E0,E1,E2,NS,GL,RH0 *0003694
      COMMON /ERRMZ0/ IERROR                                            0003695
      COMMON /PRINZ0/ IPEMSG,IPPARM                                     0003696
c     COMMON /WORKZ0/ BUFF(15),ANGS(4,4)                                0003697
	COMMON /WORKZ0/ BUFF(15)
	COMMON /WK08Z0/ ANGS(4,4)
      real*4 rangs1,rangs2,rangs3,rangs4
      equivalence (rangs1,angs(4,1))
      equivalence (rangs2,angs(4,2))
      equivalence (rangs3,angs(4,3))
      equivalence (rangs4,angs(4,4))
      DIMENSION DATA(1),GEOG(1),PROJ(1),ANG1(4),ANG2(4,2)               0003698
      EQUIVALENCE (ANG1(1),ANGS(1,1)) , (ANG2(1,1),ANGS(1,3))           0003699
      DATA ZERO,ONE /0.0D0,1.0D0/                                       0003700
      DATA EPSLN /1.0D-10/                                              0003701
      DATA SWITCH /0/                                                   0003702
C                                                                       0003703
C ......................................................................0003704
C       .  INITIALIZATION OF PROJECTION PARAMETERS (FILE INPUT)  .      0003705
C ......................................................................0003706
C                                                                       0003707
      ENTRY IF08Z0 (INFILE,DATA)                                        0003708
C                                                                       0003709
      IERROR = 0                                                        0003710
      READ (INFILE,END=240) ZONE,BUFF                                   0003711
      IF (SWITCH.NE.0 .AND. SWITCH.EQ.ZONE) RETURN                      0003712
  020 IF (BUFF(1) .LE. ZERO) GO TO 100                                  0003713
      A = BUFF(1)                                                       0003714
      B = BUFF(2)                                                       0003715
      IF (B .GT. ZERO) GO TO 040                                        0003716
      E = ZERO                                                          0003717
      ES = ZERO                                                         0003718
      E0 = ONE                                                          0003719
      E1 = ZERO                                                         0003720
      E2 = ZERO                                                         0003721
      GO TO 120                                                         0003722
  040 IF (B .GT. ONE) GO TO 060                                         0003723
      E = DSQRT (B)                                                     0003724
      ES = B                                                            0003725
      GO TO 080                                                         0003726
  060 ES = ONE - (B / A) ** 2                                           0003727
      E = DSQRT (ES)                                                    0003728
  080 E0 = E0FNZ0 (ES)                                                  0003729
      E1 = E1FNZ0 (ES)                                                  0003730
      E2 = E2FNZ0 (ES)                                                  0003731
      GO TO 120                                                         0003732
  100 A = AZ                                                            0003733
      E = EZ                                                            0003734
      ES = ESZ                                                          0003735
      E0 = E0Z                                                          0003736
      E1 = E1Z                                                          0003737
      E2 = E2Z                                                          0003738
  120 LAT1 = PAKRZ0 (BUFF(3))                                           0003739
      LAT2 = PAKRZ0 (BUFF(4))                                           0003740
      IF (DABS(LAT1+LAT2) .GE. EPSLN) GO TO 130                         0003741
      IF (IPEMSG .EQ. 0) PRINT 2000                                     0003742
 2000 FORMAT (' ERROR PJ08Z0'/                                          0003743
     .        ' EQUAL LATITUDES FOR ST. PARALLELS ON OPPOSITE',         0003744
     .        ' SIDES OF EQUATOR')                                      0003745
      IERROR = 081                                                      0003746
      RETURN                                                            0003747
  130 LON0 = PAKRZ0 (BUFF(5))                                           0003748
      LAT0 = PAKRZ0 (BUFF(6))                                           0003749
      X0 = BUFF(7)                                                      0003750
      Y0 = BUFF(8)                                                      0003751
      SINPHI = DSIN (LAT1)                                              0003752
      COSPHI = DCOS (LAT1)                                              0003753
      MS1 = MSFNZ0 (E,SINPHI,COSPHI)                                    0003754
      ML1 = MLFNZ0 (E0,E1,E2,LAT1)                                      0003755
      IND = 0                                                           0003756
      IF (BUFF(9) .NE. ZERO) GO TO 140                                  0003757
      NS = SINPHI                                                       0003758
      GO TO 160                                                         0003759
  140 IND = 1                                                           0003760
      SINPHI = DSIN (LAT2)                                              0003761
      COSPHI = DCOS (LAT2)                                              0003762
      MS2 = MSFNZ0 (E,SINPHI,COSPHI)                                    0003763
      ML2 = MLFNZ0 (E0,E1,E2,LAT2)                                      0003764
      IF (DABS(LAT1-LAT2) .GE. EPSLN) GO TO 150                         0003765
      NS = SINPHI                                                       0003766
      GO TO 160                                                         0003767
  150 NS = (MS1 - MS2) / (ML2 - ML1)                                    0003768
  160 GL = ML1 + MS1 / NS                                               0003769
      ML0 = MLFNZ0 (E0,E1,E2,LAT0)                                      0003770
      RH0 = A * (GL - ML0)                                              0003771
C                                                                       0003772
C LIST RESULTS OF PARAMETER INITIALIZATION.                             0003773
C                                                                       0003774
      CALL RADDZ0 (LAT1,ANGS(1,1))                                      0003775
      CALL RADDZ0 (LAT2,ANGS(1,2))                                      0003776
      CALL RADDZ0 (LON0,ANGS(1,3))                                      0003777
      CALL RADDZ0 (LAT0,ANGS(1,4))                                      0003778
      IF (IND .EQ. 0) GO TO 200                                         0003779
c     IF (IPPARM .EQ. 0) PRINT 2010, A,ES,ANGS,X0,Y0                    0003780
      IF (IPPARM .EQ. 0) PRINT 2010, A,ES,angs(1,1),angs(2,1),
     .   angs(3,1),rangs1,angs(1,2),angs(2,2),angs(3,2),rangs2,
     .   angs(1,3),angs(2,3),angs(3,3),rangs3,angs(1,4),
     .   angs(2,4),angs(3,4),rangs4,X0,Y0
 2010 FORMAT (' INITIALIZATION PARAMETERS (EQUIDISTANT CONIC',          0003781
     .        ' PROJECTION)'/                                           0003782
     .        ' SEMI-MAJOR AXIS OF ELLIPSOID =',F12.2,' METERS'/        0003783
     .        ' ECCENTRICITY SQUARED         =',F12.9/                  0003784
     .        ' LATITUDE OF 1ST ST. PARALLEL = ',A1,2I3,F7.3/           0003785
     .        ' LATITUDE OF 2ND ST. PARALLEL = ',A1,2I3,F7.3/           0003786
     .        ' LONGITUDE OF ORIGIN          = ',A1,2I3,F7.3/           0003787
     .        ' LATITUDE OF ORIGIN           = ',A1,2I3,F7.3/           0003788
     .        ' FALSE EASTING                =',F12.2,' METERS'/        0003789
     .        ' FALSE NORTHING               =',F12.2,' METERS')        0003790
      GO TO 220                                                         0003791
c 200 IF (IPPARM .EQ. 0) PRINT 2020, A,ES,ANG1,ANG2,X0,Y0               0003792
  200 IF (IPPARM .EQ. 0) PRINT 2020, A,ES,angs(1,1),angs(2,1),
     .   angs(3,1),rangs1,
     .   angs(1,3),angs(2,3),angs(3,3),rangs3,angs(1,4),
     .   angs(2,4),angs(3,4),rangs4,X0,Y0
 2020 FORMAT (' INITIALIZATION PARAMETERS (EQUIDISTANT CONIC',          0003793
     .        ' PROJECTION)'/                                           0003794
     .        ' SEMI-MAJOR AXIS OF ELLIPSOID =',F12.2,' METERS'/        0003795
     .        ' ECCENTRICITY SQUARED         =',F12.9/                  0003796
     .        ' LATITUDE OF ST. PARALLEL     = ',A1,2I3,F7.3/           0003797
     .        ' LONGITUDE OF ORIGIN          = ',A1,2I3,F7.3/           0003798
     .        ' LATITUDE OF ORIGIN           = ',A1,2I3,F7.3/           0003799
     .        ' FALSE EASTING                =',F12.2,' METERS'/        0003800
     .        ' FALSE NORTHING               =',F12.2,' METERS')        0003801
  220 DATA(1) = A                                                       0003802
      DATA(2) = ES                                                      0003803
      SWITCH = ZONE                                                     0003804
      RETURN                                                            0003805
  240 IF (IPEMSG .EQ. 0) PRINT 2030                                     0003806
 2030 FORMAT (' ERROR PJ08Z0'/                                          0003807
     .        ' MISSING PROJECTION PARAMETERS')                         0003808
      IERROR = 082                                                      0003809
      RETURN                                                            0003810
C                                                                       0003811
C ......................................................................0003812
C      .  INITIALIZATION OF PROJECTION PARAMETERS (ENTRY INPUT)  .      0003813
C ......................................................................0003814
C                                                                       0003815
      ENTRY IS08Z0 (ZZONE,DATA)                                         0003816
	zone = zzone
C                                                                       0003817
      IERROR = 0                                                        0003818
      IF (SWITCH.NE.0 .AND. SWITCH.EQ.ZONE) RETURN                      0003819
      DO 260 I = 1,9                                                    0003820
      BUFF(I) = DATA(I)                                                 0003821
  260 CONTINUE                                                          0003822
      GO TO 020                                                         0003823
C                                                                       0003824
C ......................................................................0003825
C                      .  FORWARD TRANSFORMATION  .                     0003826
C ......................................................................0003827
C                                                                       0003828
      ENTRY PF08Z0 (GEOG,PROJ)                                          0003829
C                                                                       0003830
      IERROR = 0                                                        0003831
      IF (SWITCH .NE. 0) GO TO 300                                      0003832
      IF (IPEMSG .EQ. 0) PRINT 2030                                     0003833
      IERROR = 083                                                      0003834
      RETURN                                                            0003835
  300 ML = MLFNZ0 (E0,E1,E2,GEOG(2))                                    0003836
      RH = A * (GL - ML)                                                0003837
      THETA = NS * ADJLZ0 (GEOG(1) - LON0)                              0003838
      PROJ(1) = X0 + RH * DSIN (THETA)                                  0003839
      PROJ(2) = Y0 + RH0 - RH * DCOS (THETA)                            0003840
      RETURN                                                            0003841
C                                                                       0003842
C ......................................................................0003843
C                      .  INVERSE TRANSFORMATION  .                     0003844
C ......................................................................0003845
C                                                                       0003846
      ENTRY PI08Z0 (PROJ,GEOG)                                          0003847
C                                                                       0003848
      IERROR = 0                                                        0003849
      IF (SWITCH .NE. 0) GO TO 320                                      0003850
      IF (IPEMSG .EQ. 0) PRINT 2030                                     0003851
      IERROR = 084                                                      0003852
      RETURN                                                            0003853
  320 X = PROJ(1) - X0                                                  0003854
      Y = RH0 - PROJ(2) + Y0                                            0003855
      RH = DSIGN (DSQRT (X * X + Y * Y) , NS)                           0003856
      THETA = ZERO                                                      0003857
      CON = DSIGN (ONE , NS)                                            0003858
      IF (RH .NE. ZERO) THETA = DATAN2 (CON * X , CON * Y)              0003859
      ML = GL - RH / A                                                  0003860
      GEOG(2) = PHI3Z0 (ML,E0,E1,E2)                                    0003861
      IF (IERROR .EQ. 0) GO TO 340                                      0003862
      IERROR = 085                                                      0003863
      RETURN                                                            0003864
  340 GEOG(1) = ADJLZ0 (LON0 + THETA / NS)                              0003865
      RETURN                                                            0003866
C                                                                       0003867
      END                                                               0003868

C                   PJ09Z0                                                     
C **********************************************************************0003870
C ** U.S.G.S. GENERAL MAP PROJECTION PACKAGE ...... DR. A. A. ELASSAL **0003871
C **          MATHEMATICAL ANALYSIS BY JOHN SNYDER                    **0003872
C ** MODULE I                VERSION 1.0.2                MAY 14,1981 **0003873
C **********************************************************************0003874
C                       *  TRANSVERSE MARCATOR  *                       0003875
C **********************************************************************0003876
C                                                                       0003877
      SUBROUTINE PJ09Z0                                                 0003878
C                                                                       0003879
      IMPLICIT REAL*8 (A-Z)                                             0003880
	integer*4 zzone
      INTEGER*4 IERROR,IPEMSG,IPPARM                                    0003881
      INTEGER*4 SWITCH,I,ZONE,ANGS,INFILE,IND,NIT                       0003882
      COMMON /ELLPZ0/ AZ,EZ,ESZ,E0Z,E1Z,E2Z,E3Z                         0003883
C **** PARAMETERS **** A,E,ES,KS0,LON0,LAT0,X0,Y0,E0,E1,E2,ESP,ML0,IND *0003884
      COMMON /ERRMZ0/ IERROR                                            0003885
      COMMON /PRINZ0/ IPEMSG,IPPARM                                     0003886
c     COMMON /WORKZ0/ BUFF(15),ANGS(4,2)                                0003887
	COMMON /WORKZ0/ BUFF(15)
	COMMON /WK09Z0/ ANGS(4,2)
      real*4 rangs1,rangs2
      equivalence (rangs1,angs(4,1))
      equivalence (rangs2,angs(4,2))
      DIMENSION DATA(15),GEOG(1),PROJ(1)                                0003888
      DATA ZERO,HALF,ONE,TWO,THREE /0.0D0,0.5D0,1.0D0,2.0D0,3.0D0/      0003889
      DATA FOUR,FIVE,SIX,EIGHT,NINE /4.0D0,5.0D0,6.0D0,8.0D0,9.0D0/     0003890
      DATA HALFPI /1.5707963267948966D0/                                0003891
      DATA TEN /10.0D0/                                                 0003892
      DATA TOL,EPSLN,NIT /1.0D-5,1.0D-10,6/                             0003893
      DATA SWITCH /0/                                                   0003894
C                                                                       0003895
C ......................................................................0003896
C       .  INITIALIZATION OF PROJECTION PARAMETERS (FILE INPUT)  .      0003897
C ......................................................................0003898
C                                                                       0003899
      ENTRY IF09Z0 (INFILE,data)                                        0003900
C                                                                       0003901
      IERROR = 0                                                        0003902
      READ (INFILE,END=160) ZONE,BUFF                                   0003903
      IF (SWITCH.NE.0 .AND. SWITCH.EQ.ZONE) RETURN                      0003903
  020 IF (BUFF(1) .LE. ZERO) GO TO 100                                  0003904
      A = BUFF(1)                                                       0003905
      B = BUFF(2)                                                       0003906
      IF (B .GT. ZERO) GO TO 040                                        0003907
      E = ZERO                                                          0003908
      ES = ZERO                                                         0003909
      E0 = ONE                                                          0003910
      E1 = ZERO                                                         0003911
      E2 = ZERO                                                         0003912
      GO TO 120                                                         0003913
  040 IF (B .GT. ONE) GO TO 060                                         0003914
      E = DSQRT (B)                                                     0003915
      ES = B                                                            0003916
      GO TO 080                                                         0003917
  060 ES = ONE - (B / A) ** 2                                           0003918
      E = DSQRT (ES)                                                    0003919
  080 E0 = E0FNZ0 (ES)                                                  0003920
      E1 = E1FNZ0 (ES)                                                  0003921
      E2 = E2FNZ0 (ES)                                                  0003922
      GO TO 120                                                         0003923
  100 A = AZ                                                            0003924
      E = EZ                                                            0003925
      ES = ESZ                                                          0003926
      E0 = E0Z                                                          0003927
      E1 = E1Z                                                          0003928
      E2 = E2Z                                                          0003929
  120 KS0 = BUFF(3)                                                     0003930
      LON0 = PAKRZ0 (BUFF(5))                                           0003931
      LAT0 = PAKRZ0 (BUFF(6))                                           0003932
      X0 = BUFF(7)                                                      0003933
      Y0 = BUFF(8)                                                      0003934
      ML0 = A * MLFNZ0 (E0,E1,E2,LAT0)                                  0003935
      IND = 1                                                           0003936
      IF (E .LT. TOL) GO TO 130                                         0003937
      IND = 0                                                           0003938
      ESP = ES / (ONE - ES)                                             0003939
C                                                                       0003940
C LIST RESULTS OF PARAMETER INITIALIZATION.                             0003941
C                                                                       0003942
  130 CALL RADDZ0 (LON0,ANGS(1,1))                                      0003943
      CALL RADDZ0 (LAT0,ANGS(1,2))                                      0003944
c     IF (IPPARM .EQ. 0) PRINT 2000, A,ES,KS0,ANGS,X0,Y0                0003945
      IF (IPPARM .EQ. 0) PRINT 2000, A,ES,ks0,angs(1,1),angs(2,1),
     .   angs(3,1),rangs1,angs(1,2),angs(2,2),angs(3,2),rangs2,
     .   X0,Y0
 2000 FORMAT (' INITIALIZATION PARAMETERS (TRANSVERSE MERCATOR',        0003946
     .        ' PROJECTION)'/                                           0003947
     .        ' SEMI-MAJOR AXIS OF ELLIPSOID =',F12.2,' METERS'/        0003948
     .        ' ECCENTRICITY SQUARED         =',F12.9/                  0003949
     .        ' SCALE FACTOR AT C. MERIDIAN  =',F9.6/                   0003950
     .        ' LONGITUDE OF C. MERIDIAN     = ',A1,2I3,F7.3/           0003951
     .        ' LATITUDE OF ORIGIN           = ',A1,2I3,F7.3/           0003952
     .        ' FALSE EASTING                =',F12.2,' METERS'/        0003953
     .        ' FALSE NORTHING               =',F12.2,' METERS')        0003954
      DATA(1) = A                                                       0003955
      DATA(2) = ES                                                      0003956
      SWITCH = ZONE                                                     0003957
      RETURN                                                            0003958
  160 IF (IPEMSG .EQ. 0) PRINT 2010                                     0003959
 2010 FORMAT (' ERROR PJ09Z0'/                                          0003960
     .        ' MISSING PROJECTION PARAMETERS')                         0003961
      IERROR = 091                                                      0003962
      RETURN                                                            0003963
C                                                                       0003964
C ......................................................................0003965
C      .  INITIALIZATION OF PROJECTION PARAMETERS (ENTRY INPUT)  .      0003966
C ......................................................................0003967
C                                                                       0003968
      ENTRY IS09Z0 (ZZONE,DATA)                                         0003969
	zone = zzone
C                                                                       0003970
      IERROR = 0                                                        0003971
C$$$$$$$$$$$$$$$$$$$$$ ADDITIONS BY JFWAANANEN 5/1/81 $$$$$$$$$$$       0003971
      IF (DATA(1).NE.0.0D0.AND.DATA(1).NE.BUFF(1)) SWITCH=0             0003971
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$     0003971
      IF (SWITCH.NE.0 .AND. SWITCH.EQ.ZONE) RETURN                      0003972
      DO 180 I = 1,8                                                    0003973
      BUFF(I) = DATA(I)                                                 0003974
  180 CONTINUE                                                          0003975
      GO TO 020                                                         0003976
C                                                                       0003977
C ......................................................................0003978
C                      .  FORWARD TRANSFORMATION  .                     0003979
C ......................................................................0003980
C                                                                       0003981
      ENTRY PF09Z0 (GEOG,PROJ)                                          0003982
C                                                                       0003983
      IERROR = 0                                                        0003984
      IF (SWITCH .NE. 0) GO TO 220                                      0003985
      IF (IPEMSG .EQ. 0) PRINT 2010                                     0003986
      IERROR = 092                                                      0003987
      RETURN                                                            0003988
  220 DLON = ADJLZ0 (GEOG(1) - LON0)                                    0003989
      LAT = GEOG(2)                                                     0003990
      IF (IND .EQ. 0) GO TO 240                                         0003991
      COSPHI = DCOS (LAT)                                               0003992
      B = COSPHI * DSIN (DLON)                                          0003993
      IF (DABS(DABS(B) - ONE) .GT. EPSLN) GO TO 230                     0003994
      IF (IPEMSG .EQ. 0) PRINT 2020                                     0003995
 2020 FORMAT (' ERROR PJ09Z0'/                                          0003996
     .        ' POINT PROJECTS INTO INFINITY')                          0003997
      IERROR = 093                                                      0003998
      RETURN                                                            0003999
  230 PROJ(1) = HALF * A * KS0 * DLOG ((ONE + B) / (ONE - B))           0004000
      CON =  DACOS (COSPHI * DCOS (DLON) / DSQRT (ONE - B * B))         0004001
      IF (LAT .LT. ZERO) CON =-CON                                      0004002
      PROJ(2) = A * KS0 * (CON - LAT0)                                  0004003
      RETURN                                                            0004004
C                                                                       0004005
  240 SINPHI = DSIN (LAT)                                               0004006
      COSPHI = DCOS (LAT)                                               0004007
      AL = COSPHI * DLON                                                0004008
      ALS = AL * AL                                                     0004009
      C = ESP * COSPHI * COSPHI                                         0004010
      TQ = DTAN (LAT)                                                   0004011
      T = TQ * TQ                                                       0004012
      N = A / DSQRT (ONE - ES * SINPHI * SINPHI)                        0004013
      ML = A * MLFNZ0 (E0,E1,E2,LAT)                                    0004014
      PROJ(1) = KS0 * N * AL * (ONE + ALS / SIX * (ONE - T + C +        0004015
     .          ALS / 20.0D0 * (FIVE - 18.0D0 * T + T * T + 72.0D0 *    0004016
     .          C - 58.0D0 * ESP))) + X0                                0004017
      PROJ(2) = KS0 * (ML - ML0 + N * TQ * (ALS * (HALF + ALS / 24.0D0 *0004018
     .          (FIVE - T + NINE * C + FOUR * C * C + ALS / 30.0D0 *    0004019
     .          (61.0D0 - 58.0D0 * T + T * T + 600.0D0 * C -            0004020
     .          330.0D0 * ESP))))) + Y0                                 0004021
      RETURN                                                            0004022
C                                                                       0004023
C ......................................................................0004024
C                      .  INVERSE TRANSFORMATION  .                     0004025
C ......................................................................0004026
C                                                                       0004027
      ENTRY PI09Z0 (PROJ,GEOG)                                          0004028
C                                                                       0004029
      IERROR = 0                                                        0004030
      IF (SWITCH .NE. 0) GO TO 320                                      0004031
      IF (IPEMSG .EQ. 0) PRINT 2010                                     0004032
      IERROR = 094                                                      0004033
      RETURN                                                            0004034
  320 X = PROJ(1) - X0                                                  0004035
      Y = PROJ(2) - Y0                                                  0004036
      IF (IND .EQ. 0) GO TO 340                                         0004037
      F = DEXP (X / (A * KS0))                                          0004038
      G = HALF * (F - ONE / F)                                          0004039
      H = DCOS (LAT0 + Y / (A * KS0))                                   0004040
      CON = DSQRT ((ONE - H * H) / (ONE + G * G))                       0004041
      GEOG(2) =  DASIN (CON)                                            0004042
      IF (Y .LT. ZERO) GEOG(2) =-GEOG(2)                                0004043
      IF (G.NE.ZERO .OR. H.NE.ZERO) GO TO 330                           0004044
      GEOG(1) = LON0                                                    0004045
      RETURN                                                            0004046
  330 GEOG(1) = ADJLZ0 (DATAN2 (G,H) + LON0)                            0004047
      RETURN                                                            0004048
C                                                                       0004049
  340 CON = (ML0 + Y / KS0) / A                                         0004050
      PHI = CON                                                         0004051
      DO 360 I = 1,NIT                                                  0004052
      DPHI = ((CON + E1 * DSIN (TWO * PHI) - E2 * DSIN (FOUR * PHI)) /  0004053
     .       E0) - PHI                                                  0004054
      PHI = PHI + DPHI                                                  0004055
      IF (DABS(DPHI) .LE. EPSLN) GO TO 380                              0004056
  360 CONTINUE                                                          0004057
      IF (IPEMSG .EQ. 0) PRINT 2030, NIT                                0004058
 2030 FORMAT (' ERROR PI09Z0' /                                         0004059
     .        ' LATITUDE FAILED TO CONVERGE AFTER',I3,' ITERATIONS')    0004060
      IERROR = 095                                                      0004061
      RETURN                                                            0004062
  380 IF (DABS(PHI) .LT. HALFPI) GO TO 400                              0004063
      GEOG(2) = DSIGN (HALFPI , Y)                                      0004064
      GEOG(1) = LON0                                                    0004065
      RETURN                                                            0004066
  400 SINPHI = DSIN (PHI)                                               0004067
      COSPHI = DCOS (PHI)                                               0004068
      TANPHI = DTAN (PHI)                                               0004069
      C = ESP * COSPHI * COSPHI                                         0004070
      CS = C * C                                                        0004071
      T = TANPHI * TANPHI                                               0004072
      TS = T * T                                                        0004073
      CON = ONE - ES * SINPHI * SINPHI                                  0004074
      N = A / DSQRT (CON)                                               0004075
      R = N * (ONE - ES) / CON                                          0004076
      D = X / (N * KS0)                                                 0004077
      DS = D * D                                                        0004078
      GEOG(2) = PHI - (N * TANPHI * DS / R) * (HALF - DS / 24.0D0 *     0004079
     .          (FIVE + THREE * T + TEN * C - FOUR * CS - NINE * ESP -  0004080
     .          DS / 30.0D0 * (61.0D0 + 90.0D0 * T + 298.0D0 * C +      0004081
     .          45.0D0 * TS - 252.0D0 * ESP - THREE * CS)))             0004082
      GEOG(1) = ADJLZ0 (LON0 + (D * (ONE - DS / SIX * (ONE + TWO *      0004083
     .          T + C - DS / 20.0D0 * (FIVE - TWO * C + 28.0D0 * T -    0004084
     .          THREE * CS + EIGHT * ESP + 24.0D0 * TS))) / COSPHI))    0004085
      RETURN                                                            0004086
C                                                                       0004087
      END                                                               0004088

C                   PJ10Z0                                                     
C **********************************************************************0004090
C ** U.S.G.S. GENERAL MAP PROJECTION PACKAGE ...... DR. A. A. ELASSAL **0004091
C **          MATHEMATICAL ANALYSIS BY JOHN SNYDER                    **0004092
C ** MODULE I                VERSION 1.0.0            NOVEMBER 1,1980 **0004093
C **********************************************************************0004094
C                          *  STEREOGRAPHIC  *                          0004095
C **********************************************************************0004096
C                                                                       0004097
      SUBROUTINE PJ10Z0                                                 0004098
C                                                                       0004099
      IMPLICIT REAL*8 (A-Z)                                             0004100
	integer*4 zzone
      INTEGER*4 IERROR,IPEMSG,IPPARM                                    0004101
      INTEGER*4 SWITCH,I,ZONE,ANGS,INFILE                               0004102
      COMMON /SPHRZ0/ AZZ                                               0004103
C **** PARAMETERS **** A,LON0,LAT0,X0,Y0,SINPH0,COSPH0 *****************0004104
      COMMON /ERRMZ0/ IERROR                                            0004105
      COMMON /PRINZ0/ IPEMSG,IPPARM                                     0004106
c     COMMON /WORKZ0/ BUFF(15),ANGS(4,2)                                0004107
	COMMON /WORKZ0/ BUFF(15)
	COMMON /WK10Z0/ ANGS(4,2)
      real*4 rangs1,rangs2
      equivalence (rangs1,angs(4,1))
      equivalence (rangs2,angs(4,2))
      DIMENSION DATA(1),GEOG(1),PROJ(1)                                 0004108
      DATA HALFPI /1.5707963267948966D0/                                0004109
      DATA EPSLN /1.0D-10/                                              0004110
      DATA ZERO,ONE,TWO /0.0D0,1.0D0,2.0D0/                             0004111
      DATA SWITCH /0/                                                   0004112
C                                                                       0004113
C ......................................................................0004114
C       .  INITIALIZATION OF PROJECTION PARAMETERS (FILE INPUT)  .      0004115
C ......................................................................0004116
C                                                                       0004117
      ENTRY IF10Z0 (INFILE,data)                                        0004118
C                                                                       0004119
      IERROR = 0                                                        0004120
      READ (INFILE,END=060) ZONE,BUFF                                   0004121
      IF (SWITCH.NE.0 .AND. SWITCH.EQ.ZONE) RETURN                      0004122
  020 A = BUFF(1)                                                       0004123
      IF (A .LE. ZERO) A = AZZ                                          0004124
      LON0 = PAKRZ0 (BUFF(5))                                           0004125
      LAT0 = PAKRZ0 (BUFF(6))                                           0004126
      X0 = BUFF(7)                                                      0004127
      Y0 = BUFF(8)                                                      0004128
      SINPH0 = DSIN (LAT0)                                              0004129
      COSPH0 = DCOS (LAT0)                                              0004130
C                                                                       0004131
C LIST RESULTS OF PARAMETER INITIALIZATION.                             0004132
C                                                                       0004133
      CALL RADDZ0 (LON0,ANGS(1,1))                                      0004134
      CALL RADDZ0 (LAT0,ANGS(1,2))                                      0004135
c     IF (IPPARM .EQ. 0) PRINT 2000, A,ANGS,X0,Y0                       0004136
      IF (IPPARM .EQ. 0) PRINT 2000, A,angs(1,1),angs(2,1),
     .   angs(3,1),rangs1,angs(1,2),angs(2,2),angs(3,2),rangs2,
     .   X0,Y0
 2000 FORMAT (' INITIALIZATION PARAMETERS (STEREOGRAPHIC',              0004137
     .        ' PROJECTION)'/                                           0004138
     .        ' RADIUS OF SPHERE             =',F12.2,' METERS'/        0004139
     .        ' LONGITUDE OF CENTER          = ',A1,2I3,F7.3/           0004140
     .        ' LATITUDE  OF CENTER          = ',A1,2I3,F7.3/           0004141
     .        ' FALSE EASTING                =',F12.2,' METERS'/        0004142
     .        ' FALSE NORTHING               =',F12.2,' METERS')        0004143
      DATA(1) = A                                                       0004144
      SWITCH = ZONE                                                     0004145
      RETURN                                                            0004146
  060 IF (IPEMSG .EQ. 0) PRINT 2010                                     0004147
 2010 FORMAT (' ERROR PJ10Z0'/                                          0004148
     .        ' MISSING PROJECTION PARAMETERS')                         0004149
      IERROR = 101                                                      0004150
      RETURN                                                            0004151
C                                                                       0004152
C ......................................................................0004153
C      .  INITIALIZATION OF PROJECTION PARAMETERS (ENTRY INPUT)  .      0004154
C ......................................................................0004155
C                                                                       0004156
      ENTRY IS10Z0 (ZZONE,DATA)                                         0004157
	zone = zzone
C                                                                       0004158
      IERROR = 0                                                        0004159
      IF (SWITCH.NE.0 .AND. SWITCH.EQ.ZONE) RETURN                      0004160
      DO 180 I = 1,8                                                    0004161
      BUFF(I) = DATA(I)                                                 0004162
  180 CONTINUE                                                          0004163
      GO TO 020                                                         0004164
C                                                                       0004165
C ......................................................................0004166
C                      .  FORWARD TRANSFORMATION  .                     0004167
C ......................................................................0004168
C                                                                       0004169
      ENTRY PF10Z0 (GEOG,PROJ)                                          0004170
C                                                                       0004171
      IERROR = 0                                                        0004172
      IF (SWITCH .NE. 0) GO TO 120                                      0004173
      IF (IPEMSG .EQ. 0) PRINT 2010                                     0004174
      IERROR = 102                                                      0004175
      RETURN                                                            0004176
  120 LON = ADJLZ0 (GEOG(1) - LON0)                                     0004177
      SINPHI = DSIN (GEOG(2))                                           0004178
      COSPHI = DCOS (GEOG(2))                                           0004179
      COSLON = DCOS (LON)                                               0004180
      G = SINPH0 * SINPHI + COSPH0 * COSPHI * COSLON                    0004181
      IF (DABS(G + ONE) .GT. EPSLN) GO TO 140                           0004182
      if (ipemsg .eq. 0) PRINT 2020
 2020 FORMAT (' ERROR PJ10Z0'/                                          0004184
     .        ' POINT PROJECTS INTO INFINITY')                          0004185
      IERROR = 103                                                      0004186
      RETURN                                                            0004187
  140 KSP = TWO / (ONE + G)                                             0004188
      PROJ(1) = X0 + A * KSP * COSPHI * DSIN (LON)                      0004189
      PROJ(2) = Y0 + A * KSP * (COSPH0 * SINPHI - SINPH0 * COSPHI *     0004190
     .          COSLON)                                                 0004191
      RETURN                                                            0004192
C                                                                       0004193
C ......................................................................0004194
C                      .  INVERSE TRANSFORMATION  .                     0004195
C ......................................................................0004196
C                                                                       0004197
      ENTRY PI10Z0 (PROJ,GEOG)                                          0004198
C                                                                       0004199
      IERROR = 0                                                        0004200
      IF (SWITCH .NE. 0) GO TO 220                                      0004201
      IF (IPEMSG .EQ. 0) PRINT 2010                                     0004202
      IERROR = 104                                                      0004203
      RETURN                                                            0004204
  220 X = PROJ(1) - X0                                                  0004205
      Y = PROJ(2) - Y0                                                  0004206
      RH = DSQRT (X * X + Y * Y)                                        0004207
      Z = TWO * DATAN (RH / (TWO * A))                                  0004208
      SINZ = DSIN (Z)                                                   0004209
      COSZ = DCOS (Z)                                                   0004210
      GEOG(1) = LON0                                                    0004211
      IF (DABS(RH) .GT. EPSLN) GO TO 240                                0004212
      GEOG(2) = LAT0                                                    0004213
      RETURN                                                            0004214
  240 GEOG(2) =  DASIN (COSZ * SINPH0 + Y * SINZ * COSPH0 / RH)         0004215
      CON = DABS (LAT0) - HALFPI                                        0004216
      IF (DABS (CON) .GT. EPSLN) GO TO 260                              0004217
      IF (LAT0 .LT. ZERO) GO TO 250                                     0004218
      GEOG(1) = ADJLZ0 (LON0 + DATAN2 (X , -Y))                         0004219
      RETURN                                                            0004220
  250 GEOG(1) = ADJLZ0 (LON0 - DATAN2 (-X , Y))                         0004221
      RETURN                                                            0004222
  260 CON = COSZ - SINPH0 * DSIN (GEOG(2))                              0004223
      IF (CON.EQ.ZERO .AND. X.EQ.ZERO) RETURN                           0004224
      GEOG(1) = ADJLZ0 (LON0 + DATAN2 ((X*SINZ*COSPH0) , (CON*RH)))     0004225
      RETURN                                                            0004226
C                                                                       0004227
      END                                                               0004228

C                   PJ11Z0                                                     
C **********************************************************************0004230
C ** U.S.G.S. GENERAL MAP PROJECTION PACKAGE ...... DR. A. A. ELASSAL **0004231
C **          MATHEMATICAL ANALYSIS BY JOHN SNYDER                    **0004232
C ** MODULE I                VERSION 1.0.0            NOVEMBER 1,1980 **0004233
C **********************************************************************0004234
C                   *  LAMBERT AZIMUTHAL EQUAL-AREA  *                  0004235
C **********************************************************************0004236
C                                                                       0004237
      SUBROUTINE PJ11Z0                                                 0004238
C                                                                       0004239
      IMPLICIT REAL*8 (A-Z)                                             0004240
	integer*4 zzone
      INTEGER*4 IERROR,IPEMSG,IPPARM                                    0004241
      INTEGER*4 SWITCH,I,ZONE,ANGS,INFILE                               0004242
      COMMON /SPHRZ0/ AZZ                                               0004243
C **** PARAMETERS **** A,LON0,LAT0,X0,Y0,SINPH0,COSPH0 *****************0004244
      COMMON /ERRMZ0/ IERROR                                            0004245
      COMMON /PRINZ0/ IPEMSG,IPPARM                                     0004246
c     COMMON /WORKZ0/ BUFF(15),ANGS(4,2)                                0004247
	COMMON /WORKZ0/ BUFF(15)
	COMMON /WK11Z0/ ANGS(4,2)
      real*4 rangs1,rangs2
      equivalence (rangs1,angs(4,1))
      equivalence (rangs2,angs(4,2))
      DIMENSION DATA(1),GEOG(1),PROJ(1)                                 0004248
      DATA HALFPI /1.5707963267948966D0/                                0004249
      DATA EPSLN /1.0D-10/                                              0004250
      DATA ZERO,ONE,TWO /0.0D0,1.0D0,2.0D0/                             0004251
      DATA SWITCH /0/                                                   0004252
C                                                                       0004253
C ......................................................................0004254
C       .  INITIALIZATION OF PROJECTION PARAMETERS (FILE INPUT)  .      0004255
C ......................................................................0004256
C                                                                       0004257
      ENTRY IF11Z0 (INFILE,data)                                        0004258
C                                                                       0004259
      IERROR = 0                                                        0004260
      READ (INFILE,END=060) ZONE,BUFF                                   0004261
      IF (SWITCH.NE.0 .AND. SWITCH.EQ.ZONE) RETURN                      0004262
  020 A = BUFF(1)                                                       0004263
      IF (A .LE. ZERO) A = AZZ                                          0004264
      LON0 = PAKRZ0 (BUFF(5))                                           0004265
      LAT0 = PAKRZ0 (BUFF(6))                                           0004266
      X0 = BUFF(7)                                                      0004267
      Y0 = BUFF(8)                                                      0004268
      SINPH0 = DSIN (LAT0)                                              0004269
      COSPH0 = DCOS (LAT0)                                              0004270
C                                                                       0004271
C LIST RESULTS OF PARAMETER INITIALIZATION.                             0004272
C                                                                       0004273
      CALL RADDZ0 (LON0,ANGS(1,1))                                      0004274
      CALL RADDZ0 (LAT0,ANGS(1,2))                                      0004275
c     IF (IPPARM .EQ. 0) PRINT 2000, A,ANGS,X0,Y0                       0004276
      IF (IPPARM .EQ. 0) PRINT 2000, A,angs(1,1),angs(2,1),
     .   angs(3,1),rangs1,angs(1,2),angs(2,2),angs(3,2),rangs2,
     .   X0,Y0
 2000 FORMAT (' INITIALIZATION PARAMETERS (LAMBERT AZIMUTHAL EQUAL-AREA'0004277
     .       ,' PROJECTION)'/                                           0004278
     .        ' RADIUS OF SPHERE             =',F12.2,' METERS'/        0004279
     .        ' LONGITUDE OF CENTER          = ',A1,2I3,F7.3/           0004280
     .        ' LATITUDE  OF CENTER          = ',A1,2I3,F7.3/           0004281
     .        ' FALSE EASTING                =',F12.2,' METERS'/        0004282
     .        ' FALSE NORTHING               =',F12.2,' METERS')        0004283
      DATA(1) = A                                                       0004284
      SWITCH = ZONE                                                     0004285
      RETURN                                                            0004286
  060 IF (IPEMSG .EQ. 0) PRINT 2010                                     0004287
 2010 FORMAT (' ERROR PJ11Z0'/                                          0004288
     .        ' MISSING PROJECTION PARAMETERS')                         0004289
      IERROR = 111                                                      0004290
      RETURN                                                            0004291
C                                                                       0004292
C ......................................................................0004293
C      .  INITIALIZATION OF PROJECTION PARAMETERS (ENTRY INPUT)  .      0004294
C ......................................................................0004295
C                                                                       0004296
      ENTRY IS11Z0 (ZZONE,DATA)                                         0004297
	zone = zzone
C                                                                       0004298
      IERROR = 0                                                        0004299
      IF (SWITCH.NE.0 .AND. SWITCH.EQ.ZONE) RETURN                      0004300
      DO 180 I = 1,8                                                    0004301
      BUFF(I) = DATA(I)                                                 0004302
  180 CONTINUE                                                          0004303
      GO TO 020                                                         0004304
C                                                                       0004305
C ......................................................................0004306
C                      .  FORWARD TRANSFORMATION  .                     0004307
C ......................................................................0004308
C                                                                       0004309
      ENTRY PF11Z0 (GEOG,PROJ)                                          0004310
C                                                                       0004311
      IERROR = 0                                                        0004312
      IF (SWITCH .NE. 0) GO TO 120                                      0004313
      IF (IPEMSG .EQ. 0) PRINT 2010                                     0004314
      IERROR = 112                                                      0004315
      RETURN                                                            0004316
  120 LON = ADJLZ0 (GEOG(1) - LON0)                                     0004317
      SINPHI = DSIN (GEOG(2))                                           0004318
      COSPHI = DCOS (GEOG(2))                                           0004319
      COSLON = DCOS (LON)                                               0004320
      G = SINPH0 * SINPHI + COSPH0 * COSPHI * COSLON                    0004321
      IF (G .NE. -ONE) GO TO 140                                        0004322
      CON = TWO * A                                                     0004323
      IF (IPEMSG .EQ. 0) PRINT 2020, CON                                0004324
 2020 FORMAT (' POINT PROJECTS INTO A CIRCLE OF RADIUS =',F12.2,        0004325
     .        ' METERS')                                                0004326
      IERROR = 113                                                      0004327
      RETURN                                                            0004328
  140 KSP = DSQRT (TWO / (ONE + G))                                     0004329
      PROJ(1) = X0 + A * KSP * COSPHI * DSIN (LON)                      0004330
      PROJ(2) = Y0 + A * KSP * (COSPH0 * SINPHI - SINPH0 * COSPHI *     0004331
     .          COSLON)                                                 0004332
      RETURN                                                            0004333
C                                                                       0004334
C ......................................................................0004335
C                      .  INVERSE TRANSFORMATION  .                     0004336
C ......................................................................0004337
C                                                                       0004338
      ENTRY PI11Z0 (PROJ,GEOG)                                          0004339
C                                                                       0004340
      IERROR = 0                                                        0004341
      IF (SWITCH .NE. 0) GO TO 220                                      0004342
      IF (IPEMSG .EQ. 0) PRINT 2010                                     0004343
      IERROR = 114                                                      0004344
      RETURN                                                            0004345
  220 X = PROJ(1) - X0                                                  0004346
      Y = PROJ(2) - Y0                                                  0004347
      RH = DSQRT (X * X + Y * Y)                                        0004348
      CON = RH / (TWO * A)                                              0004349
      IF (CON .LE. ONE) GO TO 230                                       0004350
      if (ipemsg .eq. 0) PRINT 2030
 2030 FORMAT (' ERROR PJ11Z0'/                                          0004352
     .        ' INPUT DATA ERROR')                                      0004353
      IERROR = 115                                                      0004354
      RETURN                                                            0004355
  230 Z = TWO *  DASIN (CON)                                            0004356
      SINZ = DSIN (Z)                                                   0004357
      COSZ = DCOS (Z)                                                   0004358
      GEOG(1) = LON0                                                    0004359
      IF (DABS(RH) .GT. EPSLN) GO TO 240                                0004360
      GEOG(2) = LAT0                                                    0004361
      RETURN                                                            0004362
  240 GEOG(2) =  DASIN (COSZ * SINPH0 + Y * SINZ * COSPH0 / RH)         0004363
      CON = DABS (LAT0) - HALFPI                                        0004364
      IF (DABS (CON) .GT. EPSLN) GO TO 260                              0004365
      IF (LAT0 .LT. ZERO) GO TO 250                                     0004366
      GEOG(1) = ADJLZ0 (LON0 + DATAN2 (X , -Y))                         0004367
      RETURN                                                            0004368
  250 GEOG(1) = ADJLZ0 (LON0 - DATAN2 (-X , Y))                         0004369
      RETURN                                                            0004370
  260 CON = COSZ - SINPH0 * DSIN (GEOG(2))                              0004371
      IF (CON .EQ. ZERO) RETURN                                         0004372
      GEOG(1) = ADJLZ0 (LON0 + DATAN2 ((X*SINZ*COSPH0) , (CON*RH)))     0004373
      RETURN                                                            0004374
C                                                                       0004375
      END                                                               0004376

C                   PJ12Z0                                                     
C **********************************************************************0004378
C ** U.S.G.S. GENERAL MAP PROJECTION PACKAGE ...... DR. A. A. ELASSAL **0004379
C **          MATHEMATICAL ANALYSIS BY JOHN SNYDER                    **0004380
C ** MODULE I                VERSION 1.0.0            NOVEMBER 1,1980 **0004381
C **********************************************************************0004382
C                      *  AZIMUTHAL EQUIDISTANT  *                      0004383
C **********************************************************************0004384
C                                                                       0004385
      SUBROUTINE PJ12Z0                                                 0004386
C                                                                       0004387
      IMPLICIT REAL*8 (A-Z)                                             0004388
	integer*4 zzone
      INTEGER*4 IERROR,IPEMSG,IPPARM                                    0004389
      INTEGER*4 SWITCH,I,ZONE,ANGS,INFILE                               0004390
      COMMON /SPHRZ0/ AZZ                                               0004391
C **** PARAMETERS **** A,LON0,LAT0,X0,Y0,SINPH0,COSPH0 *****************0004392
      COMMON /ERRMZ0/ IERROR                                            0004393
      COMMON /PRINZ0/ IPEMSG,IPPARM                                     0004394
c     COMMON /WORKZ0/ BUFF(15),ANGS(4,2)                                0004395
	COMMON /WORKZ0/ BUFF(15)
	COMMON /WK12Z0/ ANGS(4,2)
      real*4 rangs1,rangs2
      equivalence (rangs1,angs(4,1))
      equivalence (rangs2,angs(4,2))
      DIMENSION DATA(1),GEOG(1),PROJ(1)                                 0004396
      DATA HALFPI /1.5707963267948966D0/                                0004397
      DATA EPSLN /1.0D-10/                                              0004398
      DATA ZERO,ONE,TWO /0.0D0,1.0D0,2.0D0/                             0004399
      DATA SWITCH /0/                                                   0004400
C                                                                       0004401
C ......................................................................0004402
C       .  INITIALIZATION OF PROJECTION PARAMETERS (FILE INPUT)  .      0004403
C ......................................................................0004404
C                                                                       0004405
      ENTRY IF12Z0 (INFILE,data)                                        0004406
C                                                                       0004407
      IERROR = 0                                                        0004408
      READ (INFILE,END=060) ZONE,BUFF                                   0004409
      IF (SWITCH.NE.0 .AND. SWITCH.EQ.ZONE) RETURN                      0004410
  020 A = BUFF(1)                                                       0004411
      IF (A .LE. ZERO) A = AZZ                                          0004412
      LON0 = PAKRZ0 (BUFF(5))                                           0004413
      LAT0 = PAKRZ0 (BUFF(6))                                           0004414
      X0 = BUFF(7)                                                      0004415
      Y0 = BUFF(8)                                                      0004416
      SINPH0 = DSIN (LAT0)                                              0004417
      COSPH0 = DCOS (LAT0)                                              0004418
C                                                                       0004419
C LIST RESULTS OF PARAMETER INITIALIZATION.                             0004420
C                                                                       0004421
      CALL RADDZ0 (LON0,ANGS(1,1))                                      0004422
      CALL RADDZ0 (LAT0,ANGS(1,2))                                      0004423
c     IF (IPPARM .EQ. 0) PRINT 2000, A,ANGS,X0,Y0                       0004424
      IF (IPPARM .EQ. 0) PRINT 2000, A,angs(1,1),angs(2,1),
     .   angs(3,1),rangs1,angs(1,2),angs(2,2),angs(3,2),rangs2,
     .   X0,Y0
 2000 FORMAT (' INITIALIZATION PARAMETERS (AZIMUTHAL EQUIDISTANT',      0004425
     .        ' PROJECTION)'/                                           0004426
     .        ' RADIUS OF SPHERE             =',F12.2,' METERS'/        0004427
     .        ' LONGITUDE OF CENTER          = ',A1,2I3,F7.3/           0004428
     .        ' LATITUDE  OF CENTER          = ',A1,2I3,F7.3/           0004429
     .        ' FALSE EASTING                =',F12.2,' METERS'/        0004430
     .        ' FALSE NORTHING               =',F12.2,' METERS')        0004431
      DATA(1) = A                                                       0004432
      SWITCH = ZONE                                                     0004433
      RETURN                                                            0004434
  060 IF (IPEMSG .EQ. 0) PRINT 2010                                     0004435
 2010 FORMAT (' ERROR PJ12Z0'/                                          0004436
     .        ' MISSING PROJECTION PARAMETERS')                         0004437
      IERROR = 121                                                      0004438
      RETURN                                                            0004439
C                                                                       0004440
C ......................................................................0004441
C      .  INITIALIZATION OF PROJECTION PARAMETERS (ENTRY INPUT)  .      0004442
C ......................................................................0004443
C                                                                       0004444
      ENTRY IS12Z0 (ZZONE,DATA)                                         0004445
	zone = zzone
C                                                                       0004446
      IERROR = 0                                                        0004447
      IF (SWITCH.NE.0 .AND. SWITCH.EQ.ZONE) RETURN                      0004448
      DO 180 I = 1,8                                                    0004449
      BUFF(I) = DATA(I)                                                 0004450
  180 CONTINUE                                                          0004451
      GO TO 020                                                         0004452
C                                                                       0004453
C ......................................................................0004454
C                      .  FORWARD TRANSFORMATION  .                     0004455
C ......................................................................0004456
C                                                                       0004457
      ENTRY PF12Z0 (GEOG,PROJ)                                          0004458
C                                                                       0004459
      IERROR = 0                                                        0004460
      IF (SWITCH .NE. 0) GO TO 120                                      0004461
      IF (IPEMSG .EQ. 0) PRINT 2010                                     0004462
      IERROR = 122                                                      0004463
      RETURN                                                            0004464
  120 LON = ADJLZ0 (GEOG(1) - LON0)                                     0004465
      SINPHI = DSIN (GEOG(2))                                           0004466
      COSPHI = DCOS (GEOG(2))                                           0004467
      COSLON = DCOS (LON)                                               0004468
      G = SINPH0 * SINPHI + COSPH0 * COSPHI * COSLON                    0004469
      IF (DABS(DABS(G) - ONE) .GE. EPSLN) GO TO 140                     0004470
      KSP = ONE                                                         0004471
      IF (G .GE. ZERO) GO TO 160                                        0004472
      CON = TWO * HALFPI * A                                            0004473
      IF (IPEMSG .EQ. 0) PRINT 2020, CON                                0004474
 2020 FORMAT (' POINT PROJECTS INTO CIRCLE OF RADIUS =',F12.2,          0004475
     .        ' METERS')                                                0004476
      IERROR = 123                                                      0004477
      RETURN                                                            0004478
  140 Z =  DACOS (G)                                                    0004479
      KSP = Z / DSIN (Z)                                                0004480
  160 PROJ(1) = X0 + A * KSP * COSPHI * DSIN (LON)                      0004481
      PROJ(2) = Y0 + A * KSP * (COSPH0 * SINPHI - SINPH0 * COSPHI *     0004482
     .          COSLON)                                                 0004483
      RETURN                                                            0004484
C                                                                       0004485
C ......................................................................0004486
C                      .  INVERSE TRANSFORMATION  .                     0004487
C ......................................................................0004488
C                                                                       0004489
      ENTRY PI12Z0 (PROJ,GEOG)                                          0004490
C                                                                       0004491
      IERROR = 0                                                        0004492
      IF (SWITCH .NE. 0) GO TO 220                                      0004493
      IF (IPEMSG .EQ. 0) PRINT 2010                                     0004494
      IERROR = 124                                                      0004495
      RETURN                                                            0004496
  220 X = PROJ(1) - X0                                                  0004497
      Y = PROJ(2) - Y0                                                  0004498
      RH = DSQRT (X * X + Y * Y)                                        0004499
      IF (RH .LE. (TWO * HALFPI * A)) GO TO 230                         0004500
      IF (IPEMSG .EQ. 0) PRINT 2030                                     0004501
 2030 FORMAT (' ERROR PJ12Z0'/                                          0004502
     .        ' INPUT DATA ERROR')                                      0004503
      IERROR = 125                                                      0004504
      RETURN                                                            0004505
  230 Z = RH / A                                                        0004506
      SINZ = DSIN (Z)                                                   0004507
      COSZ = DCOS (Z)                                                   0004508
      GEOG(1) = LON0                                                    0004509
      IF (DABS(RH) .GT. EPSLN) GO TO 240                                0004510
      GEOG(2) = LAT0                                                    0004511
      RETURN                                                            0004512
  240 GEOG(2) =  DASIN (COSZ * SINPH0 + Y * SINZ * COSPH0 / RH)         0004513
      CON = DABS (LAT0) - HALFPI                                        0004514
      IF (DABS (CON) .GT. EPSLN) GO TO 260                              0004515
      IF (LAT0 .LT. ZERO) GO TO 250                                     0004516
      GEOG(1) = ADJLZ0 (LON0 + DATAN2 (X , -Y))                         0004517
      RETURN                                                            0004518
  250 GEOG(1) = ADJLZ0 (LON0 - DATAN2 (-X , Y))                         0004519
      RETURN                                                            0004520
  260 CON = COSZ - SINPH0 * DSIN (GEOG(2))                              0004521
      IF (CON .EQ. ZERO) RETURN                                         0004522
      GEOG(1) = ADJLZ0 (LON0 + DATAN2 ((X*SINZ*COSPH0) , (CON*RH)))     0004523
      RETURN                                                            0004524
C                                                                       0004525
      END                                                               0004526

C                   PJ13Z0                                                     
C **********************************************************************0004528
C ** U.S.G.S. GENERAL MAP PROJECTION PACKAGE ...... DR. A. A. ELASSAL **0004529
C **          MATHEMATICAL ANALYSIS BY JOHN SNYDER                    **0004530
C ** MODULE I                VERSION 1.0.0            NOVEMBER 1,1980 **0004531
C **********************************************************************0004532
C                            *  GNOMONIC  *                             0004533
C **********************************************************************0004534
C                                                                       0004535
      SUBROUTINE PJ13Z0                                                 0004536
C                                                                       0004537
      IMPLICIT REAL*8 (A-Z)                                             0004538
	integer*4 zzone
      INTEGER*4 IERROR,IPEMSG,IPPARM                                    0004539
      INTEGER*4 SWITCH,I,ZONE,ANGS,INFILE                               0004540
      COMMON /SPHRZ0/ AZZ                                               0004541
C **** PARAMETERS **** A,LON0,LAT0,X0,Y0,SINPH0,COSPH0 *****************0004542
      COMMON /ERRMZ0/ IERROR                                            0004543
      COMMON /PRINZ0/ IPEMSG,IPPARM                                     0004544
c     COMMON /WORKZ0/ BUFF(15),ANGS(4,2)                                0004545
	COMMON /WORKZ0/ BUFF(15)
	COMMON /WK13Z0/ ANGS(4,2)
      real*4 rangs1,rangs2
      equivalence (rangs1,angs(4,1))
      equivalence (rangs2,angs(4,2))
      DIMENSION DATA(1),GEOG(1),PROJ(1)                                 0004546
      DATA HALFPI /1.5707963267948966D0/                                0004547
      DATA EPSLN /1.0D-10/                                              0004548
      DATA ZERO,ONE /0.0D0,1.0D0/                                       0004549
      DATA SWITCH /0/                                                   0004550
C                                                                       0004551
C ......................................................................0004552
C       .  INITIALIZATION OF PROJECTION PARAMETERS (FILE INPUT)  .      0004553
C ......................................................................0004554
C                                                                       0004555
      ENTRY IF13Z0 (INFILE,data)                                        0004556
C                                                                       0004557
      IERROR = 0                                                        0004558
      READ (INFILE,END=060) ZONE,BUFF                                   0004559
      IF (SWITCH.NE.0 .AND. SWITCH.EQ.ZONE) RETURN                      0004560
  020 A = BUFF(1)                                                       0004561
      IF (A .LE. ZERO) A = AZZ                                          0004562
      LON0 = PAKRZ0 (BUFF(5))                                           0004563
      LAT0 = PAKRZ0 (BUFF(6))                                           0004564
      X0 = BUFF(7)                                                      0004565
      Y0 = BUFF(8)                                                      0004566
      SINPH0 = DSIN (LAT0)                                              0004567
      COSPH0 = DCOS (LAT0)                                              0004568
C                                                                       0004569
C LIST RESULTS OF PARAMETER INITIALIZATION.                             0004570
C                                                                       0004571
      CALL RADDZ0 (LON0,ANGS(1,1))                                      0004572
      CALL RADDZ0 (LAT0,ANGS(1,2))                                      0004573
c     IF (IPPARM .EQ. 0) PRINT 2000, A,ANGS,X0,Y0                       0004574
      IF (IPPARM .EQ. 0) PRINT 2000, A,angs(1,1),angs(2,1),
     .   angs(3,1),rangs1,angs(1,2),angs(2,2),angs(3,2),rangs2,
     .   X0,Y0
 2000 FORMAT (' INITIALIZATION PARAMETERS (GNOMONIC',                   0004575
     .        ' PROJECTION)'/                                           0004576
     .        ' RADIUS OF SPHERE             =',F12.2,' METERS'/        0004577
     .        ' LONGITUDE OF CENTER          = ',A1,2I3,F7.3/           0004578
     .        ' LATITUDE  OF CENTER          = ',A1,2I3,F7.3/           0004579
     .        ' FALSE EASTING                =',F12.2,' METERS'/        0004580
     .        ' FALSE NORTHING               =',F12.2,' METERS')        0004581
      DATA(1) = A                                                       0004582
      SWITCH = ZONE                                                     0004583
      RETURN                                                            0004584
  060 IF (IPEMSG .EQ. 0) PRINT 2010                                     0004585
 2010 FORMAT (' ERROR PJ13Z0'/                                          0004586
     .        ' MISSING PROJECTION PARAMETERS')                         0004587
      IERROR = 131                                                      0004588
      RETURN                                                            0004589
C                                                                       0004590
C ......................................................................0004591
C      .  INITIALIZATION OF PROJECTION PARAMETERS (ENTRY INPUT)  .      0004592
C ......................................................................0004593
C                                                                       0004594
      ENTRY IS13Z0 (ZZONE,DATA)                                         0004595
	zone = zzone
C                                                                       0004596
      IERROR = 0                                                        0004597
      IF (SWITCH.NE.0 .AND. SWITCH.EQ.ZONE) RETURN                      0004598
      DO 180 I = 1,8                                                    0004599
      BUFF(I) = DATA(I)                                                 0004600
  180 CONTINUE                                                          0004601
      GO TO 020                                                         0004602
C                                                                       0004603
C ......................................................................0004604
C                      .  FORWARD TRANSFORMATION  .                     0004605
C ......................................................................0004606
C                                                                       0004607
      ENTRY PF13Z0 (GEOG,PROJ)                                          0004608
C                                                                       0004609
      IERROR = 0                                                        0004610
      IF (SWITCH .NE. 0) GO TO 120                                      0004611
      IF (IPEMSG .EQ. 0) PRINT 2010                                     0004612
      IERROR = 132                                                      0004613
      RETURN                                                            0004614
  120 LON = ADJLZ0 (GEOG(1) - LON0)                                     0004615
      SINPHI = DSIN (GEOG(2))                                           0004616
      COSPHI = DCOS (GEOG(2))                                           0004617
      COSLON = DCOS (LON)                                               0004618
      G = SINPH0 * SINPHI + COSPH0 * COSPHI * COSLON                    0004619
      IF (G .GT. ZERO) GO TO 140                                        0004620
      IF (IPEMSG .EQ. 0) PRINT 2020                                     0004621
 2020 FORMAT (' POINT MAPS INTO INFINITY')                              0004622
      IERROR = 133                                                      0004623
      RETURN                                                            0004624
  140 KSP = ONE / G                                                     0004625
      PROJ(1) = X0 + A * KSP * COSPHI * DSIN (LON)                      0004626
      PROJ(2) = Y0 + A * KSP * (COSPH0 * SINPHI - SINPH0 * COSPHI *     0004627
     .          COSLON)                                                 0004628
      RETURN                                                            0004629
C                                                                       0004630
C ......................................................................0004631
C                      .  INVERSE TRANSFORMATION  .                     0004632
C ......................................................................0004633
C                                                                       0004634
      ENTRY PI13Z0 (PROJ,GEOG)                                          0004635
C                                                                       0004636
      IERROR = 0                                                        0004637
      IF (SWITCH .NE. 0) GO TO 220                                      0004638
      IF (IPEMSG .EQ. 0) PRINT 2010                                     0004639
      IERROR = 134                                                      0004640
      RETURN                                                            0004641
  220 X = PROJ(1) - X0                                                  0004642
      Y = PROJ(2) - Y0                                                  0004643
      RH = DSQRT (X * X + Y * Y)                                        0004644
      Z = DATAN (RH / A)                                                0004645
      SINZ = DSIN (Z)                                                   0004646
      COSZ = DCOS (Z)                                                   0004647
      GEOG(1) = LON0                                                    0004648
      IF (DABS(RH) .GT. EPSLN) GO TO 240                                0004649
      GEOG(2) = LAT0                                                    0004650
      RETURN                                                            0004651
  240 GEOG(2) =  DASIN (COSZ * SINPH0 + Y * SINZ * COSPH0 / RH)         0004652
      CON = DABS (LAT0) - HALFPI                                        0004653
      IF (DABS (CON) .GT. EPSLN) GO TO 260                              0004654
      IF (LAT0 .LT. ZERO) GO TO 250                                     0004655
      GEOG(1) = ADJLZ0 (LON0 + DATAN2 (X , -Y))                         0004656
      RETURN                                                            0004657
  250 GEOG(1) = ADJLZ0 (LON0 - DATAN2 (-X , Y))                         0004658
      RETURN                                                            0004659
  260 CON = COSZ - SINPH0 * DSIN (GEOG(2))                              0004660
      IF (CON .EQ. ZERO) RETURN                                         0004661
      GEOG(1) = ADJLZ0 (LON0 + DATAN2 ((X*SINZ*COSPH0) , (CON*RH)))     0004662
      RETURN                                                            0004663
C                                                                       0004664
      END                                                               0004665

C                   PJ14Z0                                                     
C **********************************************************************0004667
C ** U.S.G.S. GENERAL MAP PROJECTION PACKAGE ...... DR. A. A. ELASSAL **0004668
C **          MATHEMATICAL ANALYSIS BY JOHN SNYDER                    **0004669
C ** MODULE I                VERSION 1.0.0            NOVEMBER 1,1980 **0004670
C **********************************************************************0004671
C                          *  ORTHOGRAPHIC  *                           0004672
C **********************************************************************0004673
C                                                                       0004674
      SUBROUTINE PJ14Z0                                                 0004675
C                                                                       0004676
      IMPLICIT REAL*8 (A-Z)                                             0004677
	integer*4 zzone
      INTEGER*4 IERROR,IPEMSG,IPPARM                                    0004678
      INTEGER*4 SWITCH,I,ZONE,ANGS,INFILE                               0004679
      COMMON /SPHRZ0/ AZZ                                               0004680
C **** PARAMETERS **** A,LON0,LAT0,X0,Y0,SINPH0,COSPH0 *****************0004681
      COMMON /ERRMZ0/ IERROR                                            0004682
      COMMON /PRINZ0/ IPEMSG,IPPARM                                     0004683
c     COMMON /WORKZ0/ BUFF(15),ANGS(4,2)                                0004684
	COMMON /WORKZ0/ BUFF(15)
	COMMON /WK14Z0/ ANGS(4,2)
      real*4 rangs1,rangs2
      equivalence (rangs1,angs(4,1))
      equivalence (rangs2,angs(4,2))
      DIMENSION DATA(1),GEOG(1),PROJ(1)                                 0004685
      DATA HALFPI /1.5707963267948966D0/                                0004686
      DATA EPSLN /1.0D-10/                                              0004687
      DATA ZERO,ONE /0.0D0,1.0D0/                                       0004688
      DATA SWITCH /0/                                                   0004689
C                                                                       0004690
C ......................................................................0004691
C       .  INITIALIZATION OF PROJECTION PARAMETERS (FILE INPUT)  .      0004692
C ......................................................................0004693
C                                                                       0004694
      ENTRY IF14Z0 (INFILE,data)                                        0004695
C                                                                       0004696
      IERROR = 0                                                        0004697
      READ (INFILE,END=060) ZONE,BUFF                                   0004698
      IF (SWITCH.NE.0 .AND. SWITCH.EQ.ZONE) RETURN                      0004699
  020 A = BUFF(1)                                                       0004700
      IF (A .LE. ZERO) A = AZZ                                          0004701
      LON0 = PAKRZ0 (BUFF(5))                                           0004702
      LAT0 = PAKRZ0 (BUFF(6))                                           0004703
      X0 = BUFF(7)                                                      0004704
      Y0 = BUFF(8)                                                      0004705
      SINPH0 = DSIN (LAT0)                                              0004706
      COSPH0 = DCOS (LAT0)                                              0004707
C                                                                       0004708
C LIST RESULTS OF PARAMETER INITIALIZATION.                             0004709
C                                                                       0004710
      CALL RADDZ0 (LON0,ANGS(1,1))                                      0004711
      CALL RADDZ0 (LAT0,ANGS(1,2))                                      0004712
c     IF (IPPARM .EQ. 0) PRINT 2000, A,ANGS,X0,Y0                       0004713
      IF (IPPARM .EQ. 0) PRINT 2000, A,angs(1,1),angs(2,1),
     .   angs(3,1),rangs1,angs(1,2),angs(2,2),angs(3,2),rangs2,
     .   X0,Y0
 2000 FORMAT (' INITIALIZATION PARAMETERS (ORTHOGRAPHIC',               0004714
     .        ' PROJECTION)'/                                           0004715
     .        ' RADIUS OF SPHERE             =',F12.2,' METERS'/        0004716
     .        ' LONGITUDE OF CENTER          = ',A1,2I3,F7.3/           0004717
     .        ' LATITUDE  OF CENTER          = ',A1,2I3,F7.3/           0004718
     .        ' FALSE EASTING                =',F12.2,' METERS'/        0004719
     .        ' FALSE NORTHING               =',F12.2,' METERS')        0004720
      DATA(1) = A                                                       0004721
      SWITCH = ZONE                                                     0004722
      RETURN                                                            0004723
  060 IF (IPEMSG .EQ. 0) PRINT 2010                                     0004724
 2010 FORMAT (' ERROR PJ14Z0'/                                          0004725
     .        ' MISSING PROJECTION PARAMETERS')                         0004726
      IERROR = 141                                                      0004727
      RETURN                                                            0004728
C                                                                       0004729
C ......................................................................0004730
C      .  INITIALIZATION OF PROJECTION PARAMETERS (ENTRY INPUT)  .      0004731
C ......................................................................0004732
C                                                                       0004733
      ENTRY IS14Z0 (ZZONE,DATA)                                         0004734
	zone = zzone
C                                                                       0004735
      IERROR = 0                                                        0004736
      IF (SWITCH.NE.0 .AND. SWITCH.EQ.ZONE) RETURN                      0004737
      DO 180 I = 1,8                                                    0004738
      BUFF(I) = DATA(I)                                                 0004739
  180 CONTINUE                                                          0004740
      GO TO 020                                                         0004741
C                                                                       0004742
C ......................................................................0004743
C                      .  FORWARD TRANSFORMATION  .                     0004744
C ......................................................................0004745
C                                                                       0004746
      ENTRY PF14Z0 (GEOG,PROJ)                                          0004747
C                                                                       0004748
      IERROR = 0                                                        0004749
      IF (SWITCH .NE. 0) GO TO 120                                      0004750
      IF (IPEMSG .EQ. 0) PRINT 2010                                     0004751
      IERROR = 142                                                      0004752
      RETURN                                                            0004753
  120 LON = ADJLZ0 (GEOG(1) - LON0)                                     0004754
      SINPHI = DSIN (GEOG(2))                                           0004755
      COSPHI = DCOS (GEOG(2))                                           0004756
      COSLON = DCOS (LON)                                               0004757
      G = SINPH0 * SINPHI + COSPH0 * COSPHI * COSLON                    0004758
      KSP = ONE                                                         0004759
      IF (G.GT.ZERO .OR. DABS(G).LE.EPSLN) GO TO 140                    0004760
      IF (IPEMSG .EQ. 0) PRINT 2020                                     0004761
 2020 FORMAT (' POINT CANNOT BE PROJECTED')                             0004762
      IERROR = 143                                                      0004763
      RETURN                                                            0004764
  140 PROJ(1) = X0 + A * KSP * COSPHI * DSIN (LON)                      0004765
      PROJ(2) = Y0 + A * KSP * (COSPH0 * SINPHI - SINPH0 * COSPHI *     0004766
     .          COSLON)                                                 0004767
      RETURN                                                            0004768
C                                                                       0004769
C ......................................................................0004770
C                      .  INVERSE TRANSFORMATION  .                     0004771
C ......................................................................0004772
C                                                                       0004773
      ENTRY PI14Z0 (PROJ,GEOG)                                          0004774
C                                                                       0004775
      IERROR = 0                                                        0004776
      IF (SWITCH .NE. 0) GO TO 220                                      0004777
      IF (IPEMSG .EQ. 0) PRINT 2010                                     0004778
      IERROR = 144                                                      0004779
      RETURN                                                            0004780
  220 X = PROJ(1) - X0                                                  0004781
      Y = PROJ(2) - Y0                                                  0004782
      RH = DSQRT (X * X + Y * Y)                                        0004783
      IF (RH .LE. A) GO TO 230                                          0004784
      IF (IPEMSG .EQ. 0) PRINT 2030                                     0004785
 2030 FORMAT (' ERROR PJ14Z0'/                                          0004786
     .        ' INPUT DATA ERROR')                                      0004787
      IERROR = 145                                                      0004788
      RETURN                                                            0004789
  230 Z =  DASIN (RH / A)                                               0004790
      SINZ = DSIN (Z)                                                   0004791
      COSZ = DCOS (Z)                                                   0004792
      GEOG(1) = LON0                                                    0004793
      IF (DABS(RH) .GT. EPSLN) GO TO 240                                0004794
      GEOG(2) = LAT0                                                    0004795
      RETURN                                                            0004796
  240 GEOG(2) =  DASIN (COSZ * SINPH0 + Y * SINZ * COSPH0 / RH)         0004797
      CON = DABS (LAT0) - HALFPI                                        0004798
      IF (DABS (CON) .GT. EPSLN) GO TO 260                              0004799
      IF (LAT0 .LT. ZERO) GO TO 250                                     0004800
      GEOG(1) = ADJLZ0 (LON0 + DATAN2 (X , -Y))                         0004801
      RETURN                                                            0004802
  250 GEOG(1) = ADJLZ0 (LON0 - DATAN2 (-X , Y))                         0004803
      RETURN                                                            0004804
  260 CON = COSZ - SINPH0 * DSIN (GEOG(2))                              0004805
      IF (CON .EQ. ZERO) RETURN                                         0004806
      GEOG(1) = ADJLZ0 (LON0 + DATAN2 ((X*SINZ*COSPH0) , (CON*RH)))     0004807
      RETURN                                                            0004808
C                                                                       0004809
      END                                                               0004810

C                   PJ15Z0                                                     
C **********************************************************************0004812
C ** U.S.G.S. GENERAL MAP PROJECTION PACKAGE ...... DR. A. A. ELASSAL **0004813
C **          MATHEMATICAL ANALYSIS BY JOHN SNYDER                    **0004814
C ** MODULE I                VERSION 1.0.0            NOVEMBER 1,1980 **0004815
C **********************************************************************0004816
C              *  GENERAL VERTICAL NEAR-SIDE PERSPECTIVE  *             0004817
C **********************************************************************0004818
C                                                                       0004819
      SUBROUTINE PJ15Z0                                                 0004820
C                                                                       0004821
      IMPLICIT REAL*8 (A-Z)                                             0004822
	integer*4 zzone
      INTEGER*4 IERROR,IPEMSG,IPPARM                                    0004823
      INTEGER*4 SWITCH,I,ZONE,ANGS,INFILE                               0004824
      COMMON /SPHRZ0/ AZZ                                               0004825
C **** PARAMETERS **** A,P,LON0,LAT0,X0,Y0,SINPH0,COSPH0 ***************0004826
      COMMON /ERRMZ0/ IERROR                                            0004827
      COMMON /PRINZ0/ IPEMSG,IPPARM                                     0004828
c     COMMON /WORKZ0/ BUFF(15),ANGS(4,2)                                0004829
	COMMON /WORKZ0/ BUFF(15)
	COMMON /WK15Z0/ ANGS(4,2)
      real*4 rangs1,rangs2
      equivalence (rangs1,angs(4,1))
      equivalence (rangs2,angs(4,2))
      DIMENSION DATA(1),GEOG(1),PROJ(1)                                 0004830
      DATA HALFPI /1.5707963267948966D0/                                0004831
      DATA EPSLN /1.0D-10/                                              0004832
      DATA ZERO,ONE /0.0D0,1.0D0/                                       0004833
      DATA SWITCH /0/                                                   0004834
C                                                                       0004835
C ......................................................................0004836
C       .  INITIALIZATION OF PROJECTION PARAMETERS (FILE INPUT)  .      0004837
C ......................................................................0004838
C                                                                       0004839
      ENTRY IF15Z0 (INFILE,data)                                        0004840
C                                                                       0004841
      IERROR = 0                                                        0004842
      READ (INFILE,END=060) ZONE,BUFF                                   0004843
      IF (SWITCH.NE.0 .AND. SWITCH.EQ.ZONE) RETURN                      0004844
  020 A = BUFF(1)                                                       0004845
      IF (A .LE. ZERO) A = AZZ                                          0004846
      P = ONE + BUFF(3) / A                                             0004847
      LON0 = PAKRZ0 (BUFF(5))                                           0004848
      LAT0 = PAKRZ0 (BUFF(6))                                           0004849
      X0 = BUFF(7)                                                      0004850
      Y0 = BUFF(8)                                                      0004851
      SINPH0 = DSIN (LAT0)                                              0004852
      COSPH0 = DCOS (LAT0)                                              0004853
C                                                                       0004854
C LIST RESULTS OF PARAMETER INITIALIZATION.                             0004855
C                                                                       0004856
      CALL RADDZ0 (LON0,ANGS(1,1))                                      0004857
      CALL RADDZ0 (LAT0,ANGS(1,2))                                      0004858
c     IF (IPPARM .EQ. 0) PRINT 2000, A,BUFF(3),ANGS,X0,Y0               0004859
      IF (IPPARM .EQ. 0) PRINT 2000, A,buff(3),angs(1,1),angs(2,1),
     .   angs(3,1),rangs1,angs(1,2),angs(2,2),angs(3,2),rangs2,
     .   X0,Y0
 2000 FORMAT (' INITIALIZATION PARAMETERS (GENERAL VERTICAL NEAR-SIDE', 0004860
     .        ' PERSPECTIVE PROJECTION)'/                               0004861
     .        ' RADIUS OF SPHERE             =',F12.2,' METERS'/        0004862
     .        ' HEIGHT OF PERSPECTIVE POINT'/                           0004863
     .        ' ABOVE SPHERE                 =',F12.2,' METERS'/        0004864
     .        ' LONGITUDE OF CENTER          = ',A1,2I3,F7.3/           0004865
     .        ' LATITUDE  OF CENTER          = ',A1,2I3,F7.3/           0004866
     .        ' FALSE EASTING                =',F12.2,' METERS'/        0004867
     .        ' FALSE NORTHING               =',F12.2,' METERS')        0004868
      DATA(1) = A                                                       0004869
      SWITCH = ZONE                                                     0004870
      RETURN                                                            0004871
  060 IF (IPEMSG .EQ. 0) PRINT 2010                                     0004872
 2010 FORMAT (' ERROR PJ15Z0'/                                          0004873
     .        ' MISSING PROJECTION PARAMETERS')                         0004874
      IERROR = 151                                                      0004875
      RETURN                                                            0004876
C                                                                       0004877
C ......................................................................0004878
C      .  INITIALIZATION OF PROJECTION PARAMETERS (ENTRY INPUT)  .      0004879
C ......................................................................0004880
C                                                                       0004881
      ENTRY IS15Z0 (ZZONE,DATA)                                         0004882
	zone = zzone
C                                                                       0004883
      IERROR = 0                                                        0004884
      IF (SWITCH.NE.0 .AND. SWITCH.EQ.ZONE) RETURN                      0004885
      DO 180 I = 1,8                                                    0004886
      BUFF(I) = DATA(I)                                                 0004887
  180 CONTINUE                                                          0004888
      GO TO 020                                                         0004889
C                                                                       0004890
C ......................................................................0004891
C                      .  FORWARD TRANSFORMATION  .                     0004892
C ......................................................................0004893
C                                                                       0004894
      ENTRY PF15Z0 (GEOG,PROJ)                                          0004895
C                                                                       0004896
      IERROR = 0                                                        0004897
      IF (SWITCH .NE. 0) GO TO 120                                      0004898
      IF (IPEMSG .EQ. 0) PRINT 2010                                     0004899
      IERROR = 152                                                      0004900
      RETURN                                                            0004901
  120 LON = ADJLZ0 (GEOG(1) - LON0)                                     0004902
      SINPHI = DSIN (GEOG(2))                                           0004903
      COSPHI = DCOS (GEOG(2))                                           0004904
      COSLON = DCOS (LON)                                               0004905
      G = SINPH0 * SINPHI + COSPH0 * COSPHI * COSLON                    0004906
      IF (G .GE. (ONE / P)) GO TO 140                                   0004907
      IF (IPEMSG .EQ. 0) PRINT 2020                                     0004908
 2020 FORMAT (' POINT CANNOT BE PROJECTED')                             0004909
      IERROR = 153                                                      0004910
      RETURN                                                            0004911
  140 KSP = (P - ONE) / (P - G)                                         0004912
      PROJ(1) = X0 + A * KSP * COSPHI * DSIN (LON)                      0004913
      PROJ(2) = Y0 + A * KSP * (COSPH0 * SINPHI - SINPH0 * COSPHI *     0004914
     .          COSLON)                                                 0004915
      RETURN                                                            0004916
C                                                                       0004917
C ......................................................................0004918
C                      .  INVERSE TRANSFORMATION  .                     0004919
C ......................................................................0004920
C                                                                       0004921
      ENTRY PI15Z0 (PROJ,GEOG)                                          0004922
C                                                                       0004923
      IERROR = 0                                                        0004924
      IF (SWITCH .NE. 0) GO TO 220                                      0004925
      IF (IPEMSG .EQ. 0) PRINT 2010                                     0004926
      IERROR = 154                                                      0004927
      RETURN                                                            0004928
  220 X = PROJ(1) - X0                                                  0004929
      Y = PROJ(2) - Y0                                                  0004930
      RH = DSQRT (X * X + Y * Y)                                        0004931
      R = RH / A                                                        0004932
      CON = P - ONE                                                     0004933
      COM = P + ONE                                                     0004934
      IF (R .LE. DSQRT (CON / COM)) GO TO 230                           0004935
      IF (IPEMSG .EQ. 0) PRINT 2030                                     0004936
 2030 FORMAT (' ERROR PJ15Z0'/                                          0004937
     .        ' INPUT DATA ERROR')                                      0004938
      IERROR = 155                                                      0004939
      RETURN                                                            0004940
  230 SINZ = (P - DSQRT (ONE - R * R * COM / CON)) /                    0004941
     .       (CON / R + R / CON)                                        0004942
      Z =  DASIN (SINZ)                                                 0004943
      SINZ = DSIN (Z)                                                   0004944
      COSZ = DCOS (Z)                                                   0004945
      GEOG(1) = LON0                                                    0004946
      IF (DABS(RH) .GT. EPSLN) GO TO 240                                0004947
      GEOG(2) = LAT0                                                    0004948
      RETURN                                                            0004949
  240 GEOG(2) =  DASIN (COSZ * SINPH0 + Y * SINZ * COSPH0 / RH)         0004950
      CON = DABS (LAT0) - HALFPI                                        0004951
      IF (DABS (CON) .GT. EPSLN) GO TO 260                              0004952
      IF (LAT0 .LT. ZERO) GO TO 250                                     0004953
      GEOG(1) = ADJLZ0 (LON0 + DATAN2 (X , -Y))                         0004954
      RETURN                                                            0004955
  250 GEOG(1) = ADJLZ0 (LON0 - DATAN2 (-X , Y))                         0004956
      RETURN                                                            0004957
  260 CON = COSZ - SINPH0 * DSIN (GEOG(2))                              0004958
      IF (CON .EQ. ZERO) RETURN                                         0004959
      GEOG(1) = ADJLZ0 (LON0 + DATAN2 ((X*SINZ*COSPH0) , (CON*RH)))     0004960
      RETURN                                                            0004961
C                                                                       0004962
      END                                                               0004963

C                   PJ16Z0                                                     
C **********************************************************************0004965
C ** U.S.G.S. GENERAL MAP PROJECTION PACKAGE ...... DR. A. A. ELASSAL **0004966
C **          MATHEMATICAL ANALYSIS BY JOHN SNYDER                    **0004967
C ** MODULE I                VERSION 1.0.0            NOVEMBER 1,1980 **0004968
C **********************************************************************0004969
C                           *  SINUSOIDAL  *                            0004970
C **********************************************************************0004971
C                                                                       0004972
      SUBROUTINE PJ16Z0                                                 0004973
C                                                                       0004974
      IMPLICIT REAL*8 (A-Z)                                             0004975
	integer*4 zzone
      INTEGER*4 IERROR,IPEMSG,IPPARM                                    0004976
      INTEGER*4 SWITCH,I,ZONE,ANGS,INFILE                               0004977
      COMMON /SPHRZ0/ AZZ                                               0004978
C **** PARAMETERS **** A,LON0,X0,Y0 ************************************0004979
      COMMON /ERRMZ0/ IERROR                                            0004980
      COMMON /PRINZ0/ IPEMSG,IPPARM                                     0004981
c     COMMON /WORKZ0/ BUFF(15),ANGS(4)                                  0004982
	COMMON /WORKZ0/ BUFF(15)
	COMMON /WK16Z0/ ANGS(4)
      real*4 rangs
      equivalence (rangs,angs(4))
      DIMENSION DATA(1),GEOG(1),PROJ(1)                                 0004983
      DATA HALFPI /1.5707963267948966D0/                                0004984
      DATA EPSLN /1.0D-10/                                              0004985
      DATA ZERO /0.0D0/                                                 0004986
      DATA SWITCH /0/                                                   0004987
C                                                                       0004988
C ......................................................................0004989
C       .  INITIALIZATION OF PROJECTION PARAMETERS (FILE INPUT)  .      0004990
C ......................................................................0004991
C                                                                       0004992
      ENTRY IF16Z0 (INFILE,data)                                        0004993
C                                                                       0004994
      IERROR = 0                                                        0004995
      READ (INFILE,END=060) ZONE,BUFF                                   0004996
      IF (SWITCH.NE.0 .AND. SWITCH.EQ.ZONE) RETURN                      0004997
  020 A = BUFF(1)                                                       0004998
      IF (A .LE. ZERO) A = AZZ                                          0004999
      LON0 = PAKRZ0 (BUFF(5))                                           0005000
      X0 = BUFF(7)                                                      0005001
      Y0 = BUFF(8)                                                      0005002
C                                                                       0005003
C LIST RESULTS OF PARAMETER INITIALIZATION.                             0005004
C                                                                       0005005
      CALL RADDZ0 (LON0,ANGS)                                           0005006
c     IF (IPPARM .EQ. 0) PRINT 2000, A,ANGS,X0,Y0                       0005007
      IF (IPPARM .EQ. 0) PRINT 2000, A,angs(1),angs(2),
     .   angs(3),rangs,X0,Y0
 2000 FORMAT (' INITIALIZATION PARAMETERS (SINUSOIDAL',                 0005008
     .        ' PROJECTION)'/                                           0005009
     .        ' RADIUS OF SPHERE             =',F12.2,' METERS'/        0005010
     .        ' LONGITUDE OF C. MERIDIAN     = ',A1,2I3,F7.3/           0005011
     .        ' FALSE EASTING                =',F12.2,' METERS'/        0005012
     .        ' FALSE NORTHING               =',F12.2,' METERS')        0005013
      DATA(1) = A                                                       0005014
      SWITCH = ZONE                                                     0005015
      RETURN                                                            0005016
  060 IF (IPEMSG .EQ. 0) PRINT 2010                                     0005017
 2010 FORMAT (' ERROR PJ16Z0'/                                          0005018
     .        ' MISSING PROJECTION PARAMETERS')                         0005019
      IERROR = 161                                                      0005020
      RETURN                                                            0005021
C                                                                       0005022
C ......................................................................0005023
C      .  INITIALIZATION OF PROJECTION PARAMETERS (ENTRY INPUT)  .      0005024
C ......................................................................0005025
C                                                                       0005026
      ENTRY IS16Z0 (ZZONE,DATA)                                         0005027
	zone = zzone
C                                                                       0005028
      IERROR = 0                                                        0005029
      IF (SWITCH.NE.0 .AND. SWITCH.EQ.ZONE) RETURN                      0005030
      DO 080 I = 1,8                                                    0005031
      BUFF(I) = DATA(I)                                                 0005032
  080 CONTINUE                                                          0005033
      GO TO 020                                                         0005034
C                                                                       0005035
C ......................................................................0005036
C                      .  FORWARD TRANSFORMATION  .                     0005037
C ......................................................................0005038
C                                                                       0005039
      ENTRY PF16Z0 (GEOG,PROJ)                                          0005040
C                                                                       0005041
      IERROR = 0                                                        0005042
      IF (SWITCH .NE. 0) GO TO 120                                      0005043
      IF (IPEMSG .EQ. 0) PRINT 2010                                     0005044
      IERROR = 162                                                      0005045
      RETURN                                                            0005046
  120 LON = ADJLZ0 (GEOG(1) - LON0)                                     0005047
      PROJ(1) = X0 + A * LON * DCOS (GEOG(2))                           0005048
      PROJ(2) = Y0 + A * GEOG(2)                                        0005049
      RETURN                                                            0005050
C                                                                       0005051
C ......................................................................0005052
C                      .  INVERSE TRANSFORMATION  .                     0005053
C ......................................................................0005054
C                                                                       0005055
      ENTRY PI16Z0 (PROJ,GEOG)                                          0005056
C                                                                       0005057
      IERROR = 0                                                        0005058
      IF (SWITCH .NE. 0) GO TO 220                                      0005059
      IF (IPEMSG .EQ. 0) PRINT 2010                                     0005060
      IERROR = 163                                                      0005061
      RETURN                                                            0005062
  220 X = PROJ(1) - X0                                                  0005063
      Y = PROJ(2) - Y0                                                  0005064
      GEOG(2) = Y / A                                                   0005065
      IF (DABS(GEOG(2)) .LE. HALFPI) GO TO 230                          0005066
      IF (IPEMSG .EQ. 0) PRINT 2020                                     0005067
 2020 FORMAT (' ERROR PJ16Z0'/                                          0005068
     .        ' INPUT DATA ERROR')                                      0005069
      IERROR = 164                                                      0005070
      RETURN                                                            0005071
  230 CON = DABS (GEOG(2)) - HALFPI                                     0005072
      IF (DABS (CON) .GT. EPSLN) GO TO 240                              0005073
      GEOG(1) = LON0                                                    0005074
      RETURN                                                            0005075
  240 GEOG(1) = ADJLZ0 (LON0 + X / (A * DCOS (GEOG(2))))                0005076
      RETURN                                                            0005077
C                                                                       0005078
      END                                                               0005079

C                   PJ17Z0                                                     
C **********************************************************************0005081
C ** U.S.G.S. GENERAL MAP PROJECTION PACKAGE ...... DR. A. A. ELASSAL **0005082
C **          MATHEMATICAL ANALYSIS BY JOHN SNYDER                    **0005083
C ** MODULE I                VERSION 1.0.1            MAY 1 ,1981 ******0005084
C **********************************************************************0005085
C                  *  EQUIRECTANGULAR   *                               0005086
C **********************************************************************0005087
C                                                                       0005088
      SUBROUTINE PJ17Z0                                                 0005089
C                                                                       0005090
      IMPLICIT REAL*8 (A-Z)                                             0005091
	integer*4 zzone
      INTEGER*4 IERROR,IPEMSG,IPPARM                                    0005092
      INTEGER*4 SWITCH,I,ZONE,ANGS,INFILE                               0005093
      COMMON /SPHRZ0/ AZZ                                               0005094
C **** PARAMETERS **** A,LON0,X0,Y0,LAT1 *******************************0005095
      COMMON /ERRMZ0/ IERROR                                            0005096
      COMMON /PRINZ0/ IPEMSG,IPPARM                                     0005097
c     COMMON /WORKZ0/ BUFF(15),ANGS(4,2)                                0005098
	COMMON /WORKZ0/ BUFF(15)
	COMMON /WK17Z0/ ANGS(4,2)
      real*4 rangs1,rangs2
      equivalence (rangs1,angs(4,1))
      equivalence (rangs2,angs(4,2))
      DIMENSION DATA(1),GEOG(1),PROJ(1)                                 0005099
      DATA HALFPI /1.5707963267948966D0/                                0005100
      DATA ZERO /0.0D0/                                                 0005101
      DATA SWITCH /0/                                                   0005102
C                                                                       0005103
C ......................................................................0005104
C       .  INITIALIZATION OF PROJECTION PARAMETERS (FILE INPUT)  .      0005105
C ......................................................................0005106
C                                                                       0005107
      ENTRY IF17Z0 (INFILE,data)                                        0005108
C                                                                       0005109
      IERROR = 0                                                        0005110
      READ (INFILE,END=060) ZONE,BUFF                                   0005111
      IF (SWITCH.NE.0 .AND. SWITCH.EQ.ZONE) RETURN                      0005112
  020 A = BUFF(1)                                                       0005113
      IF (A .LE. ZERO) A = AZZ                                          0005114
      LAT1 = PAKRZ0 (BUFF(6))                                           0005114
      LON0 = PAKRZ0 (BUFF(5))                                           0005115
      X0 = BUFF(7)                                                      0005116
      Y0 = BUFF(8)                                                      0005117
C                                                                       0005118
C LIST RESULTS OF PARAMETER INITIALIZATION.                             0005119
C                                                                       0005120
      CALL RADDZ0 (LAT1,ANGS(1,1))                                      0005120
      CALL RADDZ0 (LON0,ANGS(1,2))                                      0005121
c     IF (IPPARM .EQ. 0) PRINT 2000, A,ANGS,X0,Y0                       0005122
      IF (IPPARM .EQ. 0) PRINT 2000, A,angs(1,1),angs(2,1),
     .   angs(3,1),rangs1,angs(1,2),angs(2,2),angs(3,2),rangs2,
     .   X0,Y0
 2000 FORMAT (' INITIALIZATION PARAMETERS (EQUIRECTANGULAR PROJECTION)'/0005123
     .        ' RADIUS OF SPHERE             =',F12.2,' METERS'/        0005125
     .        ' LATITUDE OF TRUE SCALE       = ',A1,2I2,F7.3/           0005125
     .        ' LONGITUDE OF C. MERIDIAN     = ',A1,2I3,F7.3/           0005126
     .        ' FALSE EASTING                =',F12.2,' METERS'/        0005127
     .        ' FALSE NORTHING               =',F12.2,' METERS')        0005128
      DATA(1) = A                                                       0005129
      SWITCH = ZONE                                                     0005130
      RETURN                                                            0005131
  060 IF (IPEMSG .EQ. 0) PRINT 2010                                     0005132
 2010 FORMAT (' ERROR PJ17Z0'/                                          0005133
     .        ' MISSING PROJECTION PARAMETERS')                         0005134
      IERROR = 171                                                      0005135
      RETURN                                                            0005136
C                                                                       0005137
C ......................................................................0005138
C      .  INITIALIZATION OF PROJECTION PARAMETERS (ENTRY INPUT)  .      0005139
C ......................................................................0005140
C                                                                       0005141
      ENTRY IS17Z0 (ZZONE,DATA)                                         0005142
	zone = zzone
C                                                                       0005143
      IERROR = 0                                                        0005144
      IF (SWITCH.NE.0 .AND. SWITCH.EQ.ZONE) RETURN                      0005145
      DO 080 I = 1,8                                                    0005146
      BUFF(I) = DATA(I)                                                 0005147
  080 CONTINUE                                                          0005148
      GO TO 020                                                         0005149
C                                                                       0005150
C ......................................................................0005151
C                      .  FORWARD TRANSFORMATION  .                     0005152
C ......................................................................0005153
C                                                                       0005154
      ENTRY PF17Z0 (GEOG,PROJ)                                          0005155
C                                                                       0005156
      IERROR = 0                                                        0005157
      IF (SWITCH .NE. 0) GO TO 120                                      0005158
      IF (IPEMSG .EQ. 0) PRINT 2010                                     0005159
      IERROR = 172                                                      0005160
      RETURN                                                            0005161
  120 LON = ADJLZ0 (GEOG(1) - LON0)                                     0005162
      PROJ(1) = X0 + A * LON * DCOS(LAT1)                               0005163
      PROJ(2) = Y0 + A * GEOG(2)                                        0005164
      RETURN                                                            0005165
C                                                                       0005166
C ......................................................................0005167
C                      .  INVERSE TRANSFORMATION  .                     0005168
C ......................................................................0005169
C                                                                       0005170
      ENTRY PI17Z0 (PROJ,GEOG)                                          0005171
C                                                                       0005172
      IERROR = 0                                                        0005173
      IF (SWITCH .NE. 0) GO TO 220                                      0005174
      IF (IPEMSG .EQ. 0) PRINT 2010                                     0005175
      IERROR = 173                                                      0005176
      RETURN                                                            0005177
  220 X = PROJ(1) - X0                                                  0005178
      Y = PROJ(2) - Y0                                                  0005179
      GEOG(2) = Y / A                                                   0005180
      IF (DABS(GEOG(2)) .LE. HALFPI) GO TO 240                          0005181
      IF (IPEMSG .EQ. 0) PRINT 2020                                     0005182
 2020 FORMAT (' ERROR PJ17Z0'/                                          0005183
     .        ' INPUT DATA ERROR')                                      0005184
      IERROR = 174                                                      0005185
      RETURN                                                            0005186
  240 GEOG(1) = ADJLZ0 (LON0 + X / (A * DCOS(LAT1) ))                   0005187
      RETURN                                                            0005188
C                                                                       0005189
      END                                                               0005190

C                   PJ18Z0                                                     
C **********************************************************************0005192
C ** U.S.G.S. GENERAL MAP PROJECTION PACKAGE ...... DR. A. A. ELASSAL **0005193
C **          MATHEMATICAL ANALYSIS BY JOHN SNYDER                    **0005194
C ** MODULE I                VERSION 1.0.0            NOVEMBER 1,1980 **0005195
C **********************************************************************0005196
C                       *  MILLER CYLINDRICAL  *                        0005197
C **********************************************************************0005198
C                                                                       0005199
      SUBROUTINE PJ18Z0                                                 0005200
C                                                                       0005201
      IMPLICIT REAL*8 (A-Z)                                             0005202
	integer*4 zzone
      INTEGER*4 IERROR,IPEMSG,IPPARM                                    0005203
      INTEGER*4 SWITCH,I,ZONE,ANGS,INFILE                               0005204
      COMMON /SPHRZ0/ AZZ                                               0005205
C **** PARAMETERS **** A,LON0,X0,Y0 ************************************0005206
      COMMON /ERRMZ0/ IERROR                                            0005207
      COMMON /PRINZ0/ IPEMSG,IPPARM                                     0005208
c     COMMON /WORKZ0/ BUFF(15),ANGS(4)                                  0005209
	COMMON /WORKZ0/ BUFF(15)
	COMMON /WK18Z0/ ANGS(4)
      real*4 rangs
      equivalence (rangs,angs(4))
      DIMENSION DATA(1),GEOG(1),PROJ(1)                                 0005210
      DATA FORTPI /0.78539816339744833D0/                               0005211
      DATA ZERO,ONEQ,TWOH /0.0D0,1.25D0,2.5D0/                          0005212
      DATA SWITCH /0/                                                   0005213
C                                                                       0005214
C ......................................................................0005215
C       .  INITIALIZATION OF PROJECTION PARAMETERS (FILE INPUT)  .      0005216
C ......................................................................0005217
C                                                                       0005218
      ENTRY IF18Z0 (INFILE,data)                                        0005219
C                                                                       0005220
      IERROR = 0                                                        0005221
      READ (INFILE,END=060) ZONE,BUFF                                   0005222
      IF (SWITCH.NE.0 .AND. SWITCH.EQ.ZONE) RETURN                      0005223
  020 A = BUFF(1)                                                       0005224
      IF (A .LE. ZERO) A = AZZ                                          0005225
      LON0 = PAKRZ0 (BUFF(5))                                           0005226
      X0 = BUFF(7)                                                      0005227
      Y0 = BUFF(8)                                                      0005228
C                                                                       0005229
C LIST RESULTS OF PARAMETER INITIALIZATION.                             0005230
C                                                                       0005231
      CALL RADDZ0 (LON0,ANGS)                                           0005232
c     IF (IPPARM .EQ. 0) PRINT 2000, A,ANGS,X0,Y0                       0005233
      IF (IPPARM .EQ. 0) PRINT 2000, A,angs(1),angs(2),
     .   angs(3),rangs,X0,Y0
 2000 FORMAT (' INITIALIZATION PARAMETERS (MILLER CYLINDRICAL',         0005234
     .        ' PROJECTION)'/                                           0005235
     .        ' RADIUS OF SPHERE             =',F12.2,' METERS'/        0005236
     .        ' LONGITUDE OF C. MERIDIAN     = ',A1,2I3,F7.3/           0005237
     .        ' FALSE EASTING                =',F12.2,' METERS'/        0005238
     .        ' FALSE NORTHING               =',F12.2,' METERS')        0005239
      DATA(1) = A                                                       0005240
      SWITCH = ZONE                                                     0005241
      RETURN                                                            0005242
  060 IF (IPEMSG .EQ. 0) PRINT 2010                                     0005243
 2010 FORMAT (' ERROR PJ18Z0'/                                          0005244
     .        ' MISSING PROJECTION PARAMETERS')                         0005245
      IERROR = 181                                                      0005246
      RETURN                                                            0005247
C                                                                       0005248
C ......................................................................0005249
C      .  INITIALIZATION OF PROJECTION PARAMETERS (ENTRY INPUT)  .      0005250
C ......................................................................0005251
C                                                                       0005252
      ENTRY IS18Z0 (ZZONE,DATA)                                         0005253
	zone = zzone
C                                                                       0005254
      IERROR = 0                                                        0005255
      IF (SWITCH.NE.0 .AND. SWITCH.EQ.ZONE) RETURN                      0005256
      DO 080 I = 1,8                                                    0005257
      BUFF(I) = DATA(I)                                                 0005258
  080 CONTINUE                                                          0005259
      GO TO 020                                                         0005260
C                                                                       0005261
C ......................................................................0005262
C                      .  FORWARD TRANSFORMATION  .                     0005263
C ......................................................................0005264
C                                                                       0005265
      ENTRY PF18Z0 (GEOG,PROJ)                                          0005266
C                                                                       0005267
      IERROR = 0                                                        0005268
      IF (SWITCH .NE. 0) GO TO 120                                      0005269
      IF (IPEMSG .EQ. 0) PRINT 2010                                     0005270
      IERROR = 182                                                      0005271
      RETURN                                                            0005272
  120 LON = ADJLZ0 (GEOG(1) - LON0)                                     0005273
      PROJ(1) = X0 + A * LON                                            0005274
      PROJ(2) = Y0 + A * DLOG (DTAN (FORTPI + GEOG(2) / TWOH)) * ONEQ   0005275
      RETURN                                                            0005276
C                                                                       0005277
C ......................................................................0005278
C                      .  INVERSE TRANSFORMATION  .                     0005279
C ......................................................................0005280
C                                                                       0005281
      ENTRY PI18Z0 (PROJ,GEOG)                                          0005282
C                                                                       0005283
      IERROR = 0                                                        0005284
      IF (SWITCH .NE. 0) GO TO 220                                      0005285
      IF (IPEMSG .EQ. 0) PRINT 2010                                     0005286
      IERROR = 183                                                      0005287
      RETURN                                                            0005288
  220 X = PROJ(1) - X0                                                  0005289
      Y = PROJ(2) - Y0                                                  0005290
      GEOG(1) = ADJLZ0 (LON0 + X / A)                                   0005291
      GEOG(2) = TWOH * DATAN (DEXP (Y / A / ONEQ)) - FORTPI * TWOH      0005292
      RETURN                                                            0005293
C                                                                       0005294
      END                                                               0005295

C                   PJ19Z0                                                     
C **********************************************************************0005297
C ** U.S.G.S. GENERAL MAP PROJECTION PACKAGE ...... DR. A. A. ELASSAL **0005298
C **          MATHEMATICAL ANALYSIS BY JOHN SNYDER                    **0005299
C ** MODULE I                VERSION 1.0.1            MAY 1, 1981 ******0005300
C **********************************************************************0005301
C                        *  VAN DER GRINTEN I  *                        0005302
C **********************************************************************0005303
C                                                                       0005304
      SUBROUTINE PJ19Z0                                                 0005305
C                                                                       0005306
      IMPLICIT REAL*8 (A-Z)                                             0005307
	integer*4 zzone
      INTEGER*4 IERROR,IPEMSG,IPPARM                                    0005308
      INTEGER*4 SWITCH,I,ZONE,ANGS,NIT,INFILE                           0005309
      COMMON /SPHRZ0/ AZZ                                               0005310
C **** PARAMETERS **** A,LON0,X0,Y0 ************************************0005311
      COMMON /ERRMZ0/ IERROR                                            0005312
      COMMON /PRINZ0/ IPEMSG,IPPARM                                     0005313
c     COMMON /WORKZ0/ BUFF(15),ANGS(4)                                  0005314
	COMMON /WORKZ0/ BUFF(15)
	COMMON /WK19Z0/ ANGS(4)
      real*4 rangs
      equivalence (rangs,angs(4))
      DIMENSION DATA(1),GEOG(1),PROJ(1)                                 0005315
      DATA PI /3.14159265358979323846D0/                                0005316
      DATA HALFPI /1.5707963267948966D0/                                0005317
      DATA EPSLN,TOL,NIT /1.0D-10,0.7D0,35/                             0005318
      DATA ZERO,HALF,ONE,TWO,FOUR /0.0D0,0.5D0,1.0D0,2.0D0,4.0D0/       0005319
      DATA SWITCH /0/                                                   0005320
C                                                                       0005321
C ......................................................................0005322
C       .  INITIALIZATION OF PROJECTION PARAMETERS (FILE INPUT)  .      0005323
C ......................................................................0005324
C                                                                       0005325
      ENTRY IF19Z0 (INFILE,data)                                        0005326
C                                                                       0005327
      IERROR = 0                                                        0005328
      READ (INFILE,END=060) ZONE,BUFF                                   0005329
      IF (SWITCH.NE.0 .AND. SWITCH.EQ.ZONE) RETURN                      0005330
  020 A = BUFF(1)                                                       0005331
      IF (A .LE. ZERO) A = AZZ                                          0005332
      LON0 = PAKRZ0 (BUFF(5))                                           0005333
      X0 = BUFF(7)                                                      0005334
      Y0 = BUFF(8)                                                      0005335
C                                                                       0005336
C LIST RESULTS OF PARAMETER INITIALIZATION.                             0005337
C                                                                       0005338
      CALL RADDZ0 (LON0,ANGS)                                           0005339
c     IF (IPPARM .EQ. 0) PRINT 2000, A,ANGS,X0,Y0                       0005340
      IF (IPPARM .EQ. 0) PRINT 2000, A,angs(1),angs(2),
     .   angs(3),rangs,X0,Y0
 2000 FORMAT (' INITIALIZATION PARAMETERS (VAN DER GRINTEN I',          0005341
     .        ' PROJECTION)'/                                           0005342
     .        ' RADIUS OF SPHERE             =',F12.2,' METERS'/        0005343
     .        ' LONGITUDE OF C. MERIDIAN     = ',A1,2I3,F7.3/           0005344
     .        ' FALSE EASTING                =',F12.2,' METERS'/        0005345
     .        ' FALSE NORTHING               =',F12.2,' METERS')        0005346
      DATA(1) = A                                                       0005347
      SWITCH = ZONE                                                     0005348
      RETURN                                                            0005349
  060 IF (IPEMSG .EQ. 0) PRINT 2010                                     0005350
 2010 FORMAT (' ERROR PJ19Z0'/                                          0005351
     .        ' MISSING PROJECTION PARAMETERS')                         0005352
      IERROR = 191                                                      0005353
      RETURN                                                            0005354
C                                                                       0005355
C ......................................................................0005356
C      .  INITIALIZATION OF PROJECTION PARAMETERS (ENTRY INPUT)  .      0005357
C ......................................................................0005358
C                                                                       0005359
      ENTRY IS19Z0 (ZZONE,DATA)                                         0005360
	zone = zzone
C                                                                       0005361
      IERROR = 0                                                        0005362
      IF (SWITCH.NE.0 .AND. SWITCH.EQ.ZONE) RETURN                      0005363
      DO 080 I = 1,8                                                    0005364
      BUFF(I) = DATA(I)                                                 0005365
  080 CONTINUE                                                          0005366
      GO TO 020                                                         0005367
C                                                                       0005368
C ......................................................................0005369
C                      .  FORWARD TRANSFORMATION  .                     0005370
C ......................................................................0005371
C                                                                       0005372
      ENTRY PF19Z0 (GEOG,PROJ)                                          0005373
C                                                                       0005374
      IERROR = 0                                                        0005375
      IF (SWITCH .NE. 0) GO TO 120                                      0005376
      IF (IPEMSG .EQ. 0) PRINT 2010                                     0005377
      IERROR = 192                                                      0005378
      RETURN                                                            0005379
  120 LON = ADJLZ0 (GEOG(1) - LON0)                                     0005380
      LAT = GEOG(2)                                                     0005381
      IF (DABS(LAT) .GT. EPSLN) GO TO 140                               0005382
      PROJ(1) = X0 + A * LON                                            0005383
      PROJ(2) = Y0                                                      0005384
      RETURN                                                            0005385
  140 THETA =  DASIN (DABS (LAT /HALFPI))                               0005386
      IF (DABS(LON) .GT. EPSLN) GO TO 160                               0005387
      PROJ(1) = X0                                                      0005388
      PROJ(2) = Y0 + PI * A * DSIGN( DTAN (HALF * THETA), LAT)          0005389
      RETURN                                                            0005390
  160 AL = HALF * DABS (PI / LON - LON / PI)                            0005391
      ASQ = AL * AL                                                     0005392
      SINTHT = DSIN (THETA)                                             0005393
      COSTHT = DCOS (THETA)                                             0005394
      G = COSTHT / (SINTHT + COSTHT - ONE)                              0005395
      GSQ = G * G                                                       0005396
      M = G * (TWO / SINTHT - ONE)                                      0005397
      MSQ = M * M                                                       0005398
      CON = PI * A * (AL * (G - MSQ) + DSQRT (ASQ * (G - MSQ)**2 -      0005399
     .      (MSQ + ASQ) * (GSQ - MSQ))) / (MSQ + ASQ)                   0005400
      CON = DSIGN (CON , LON)                                           0005401
      PROJ(1) = X0 + CON                                                0005402
      CON = DABS (CON / (PI * A))                                       0005403
      PROJ(2) = Y0 + DSIGN (PI * A * DSQRT (ONE - CON * CON -           0005404
     .          TWO * AL * CON) , LAT)                                  0005405
      RETURN                                                            0005406
C                                                                       0005407
C ......................................................................0005408
C                      .  INVERSE TRANSFORMATION  .                     0005409
C ......................................................................0005410
C                                                                       0005411
      ENTRY PI19Z0 (PROJ,GEOG)                                          0005412
C                                                                       0005413
      IERROR = 0                                                        0005414
      IF (SWITCH .NE. 0) GO TO 220                                      0005415
      IF (IPEMSG .EQ. 0) PRINT 2010                                     0005416
      IERROR = 193                                                      0005417
      RETURN                                                            0005418
  220 X = PROJ(1) - X0                                                  0005419
      Y = PROJ(2) - Y0                                                  0005420
      CON = DABS (Y / (PI * A))                                         0005421
      THETA = TWO * DATAN (CON)                                         0005422
      IF (DABS(X) .GT. EPSLN) GO TO 240                                 0005423
      GEOG(1) = LON0                                                    0005424
      GEOG(2) = HALFPI * DSIGN( DSIN (THETA), Y)                        0005425
      RETURN                                                            0005426
  240 IF (DABS(Y) .GT. EPSLN) GO TO 260                                 0005427
      GEOG(1) = ADJLZ0 (LON0 + X / A)                                   0005428
      GEOG(2) = ZERO                                                    0005429
      RETURN                                                            0005430
  260 IF (DSQRT(X*X+Y*Y) .LE. PI*A) GO TO 270                           0005431
      IF (IPEMSG .EQ. 0) PRINT 2020                                     0005432
 2020 FORMAT (' ERROR PI19Z0'/                                          0005433
     .        ' INPUT DATA ERROR')                                      0005434
      IERROR = 194                                                      0005435
      RETURN                                                            0005436
  270 CNN = CON * CON                                                   0005437
      COM = DABS (X / (PI * A))                                         0005438
      CMM = COM * COM                                                   0005439
      AL = (ONE - CMM - CNN) / (TWO * COM)                              0005440
      GEOG(1) = ADJLZ0 (LON0 + DSIGN (PI*(-AL + DSQRT (AL*AL+ONE)) , X))0005441
      PHI = THETA                                                       0005442
      IF (CON .GT. TOL) GO TO 320                                       0005443
C                                                                       0005444
C LOW LATITUDE CASE                                                     0005445
C                                                                       0005446
      DO 280 I = 1,NIT                                                  0005447
      THETA =  DASIN (PHI / HALFPI)                                     0005448
      SINTHT = DSIN (THETA)                                             0005449
      COSTHT = DCOS (THETA)                                             0005450
      G = COSTHT / (SINTHT + COSTHT - ONE)                              0005451
      D = CON / SINTHT - ONE / (ONE + COSTHT)                           0005452
      H = TWO - SINTHT                                                  0005453
      J = DTAN (HALF * THETA)                                           0005454
      DPHI = (CMM + CNN - TWO * D * G * H - J * J) * PI * COSTHT /      0005455
     .       (FOUR * (G * H * (CON * COSTHT / (ONE - COSTHT) + J) /     0005456
     .       (ONE + COSTHT) + D * G * ((ONE + TWO * COSTHT * COSTHT) /  0005457
     .       COSTHT + H * (COSTHT - SINTHT) / (SINTHT + COSTHT - ONE)) -0005458
     .       J * (J * J + ONE)))                                        0005459
      PHI = PHI - DPHI                                                  0005460
      IF (DABS(DPHI) .LT. EPSLN) GO TO 400                              0005461
  280 CONTINUE                                                          0005462
  300 IF (IPEMSG .EQ. 0) PRINT 2030, NIT                                0005463
 2030 FORMAT (' ERROR PI19Z0'/                                          0005464
     .        ' LATITUDE FAILED TO CONVERGE AFTER',I3,' ITERATIONS')    0005465
      IERROR = 195                                                      0005466
      RETURN                                                            0005467
C                                                                       0005468
C HIGH LATITUDE CASE.                                                   0005469
C                                                                       0005470
  320 LON = ADJLZ0 (GEOG(1) - LON0)                                     0005471
      DO 380 I = 1,NIT                                                  0005472
      IF (DABS(PHI) .GT. EPSLN) GO TO 330                               0005473
      Y1 = ZERO                                                         0005474
      GO TO 360                                                         0005475
  330 THETA =  DASIN (DABS (PHI /HALFPI))                               0005476
      IF (DABS(LON) .GT. EPSLN) GO TO 340                               0005477
      Y1 = PI * A * DTAN (HALF * THETA)                                 0005478
      GO TO 360                                                         0005479
  340 AL = HALF * DABS (PI / LON - LON / PI)                            0005480
      ASQ = AL * AL                                                     0005481
      SINTHT = DSIN (THETA)                                             0005482
      COSTHT = DCOS (THETA)                                             0005483
      G = COSTHT / (SINTHT + COSTHT - ONE)                              0005484
      GSQ = G * G                                                       0005485
      M = G * (TWO / SINTHT - ONE)                                      0005486
      MSQ = M * M                                                       0005487
      CON = DABS ((AL * (G - MSQ) + DSQRT (ASQ * (G - MSQ)**2 -         0005488
     .      (MSQ + ASQ) * (GSQ - MSQ))) / (MSQ + ASQ))                  0005489
      Y1 = DSIGN (PI * A * DSQRT (ONE - CON * CON -                     0005490
     .          TWO * AL * CON) , PHI)                                  0005491
  360 DPHI = ((DABS(Y) - Y1) / (PI * A - Y1)) * (HALFPI - PHI)          0005492
      PHI = PHI + DPHI                                                  0005493
      IF (DABS(DPHI) .LT. EPSLN) GO TO 400                              0005494
  380 CONTINUE                                                          0005495
      GO TO 300                                                         0005496
  400 GEOG(2) = DSIGN (PHI , Y)                                         0005497
      RETURN                                                            0005498
C                                                                       0005499
      END                                                               0005500

C                   PJ20Z0                                                     
C **********************************************************************0005502
C ** U.S.G.S. GENERAL MAP PROJECTION PACKAGE ...... DR. A. A. ELASSAL **0005503
C **          MATHEMATICAL ANALYSIS BY JOHN SNYDER                    **0005504
C ** MODULE I                VERSION 1.0.1            MAY 1, 1981 ******0005505
C **********************************************************************0005506
C                    *  OBLIQUE MERCATOR (HOTINE)  *                    0005507
C **********************************************************************0005508
C                                                                       0005509
      SUBROUTINE PJ20Z0                                                 0005510
C                                                                       0005511
      IMPLICIT REAL*8 (A-Z)                                             0005512
	integer*4 zzone
      INTEGER*4 IERROR,IPEMSG,IPPARM                                    0005513
      INTEGER*4 SWITCH,I,ZONE,ANGS1,ANGS2,MODE,INFILE                   0005514
      COMMON /ELLPZ0/ AZ,EZ,ESZ,E0Z,E1Z,E2Z,E3Z                         0005515
C **** PARAMETERS **** A,E,ES,KS0,ALPHA,LONC,LON1,LAT1,LON2,LAT2,LAT0 **0005516
C ********************** X0,Y0,GAMMA,LON0,AL,BL,EL *********************0005517
      COMMON /ERRMZ0/ IERROR                                            0005518
      COMMON /PRINZ0/ IPEMSG,IPPARM                                     0005519
c     COMMON /WORKZ0/ BUFF(15),ANGS1(4,5),ANGS2(4,3)                    0005520
	COMMON /WORKZ0/ BUFF(15)
	COMMON /WK20Z0/ ANGS1(4,5),ANGS2(4,3)
      real*4 rangs1,rangs2,rangs3,rangs11,rangs12,rangs13,
     . rangs14,rangs15
      equivalence (rangs1,angs2(4,1))
      equivalence (rangs2,angs2(4,2))
      equivalence (rangs3,angs2(4,3))
      equivalence (rangs11,angs1(4,1))
      equivalence (rangs12,angs1(4,2))
      equivalence (rangs13,angs1(4,3))
      equivalence (rangs14,angs1(4,4))
      equivalence (rangs15,angs1(4,5))
      DIMENSION DATA(1),GEOG(1),PROJ(1)                                 0005524
      DATA PI /3.14159265358979323846D0/                                0005521
      DATA HALFPI /1.5707963267948966D0/                                0005522
      DATA TOL,EPSLN /1.0D-7,1.0D-10/                                   0005523
      DATA ZERO,HALF,ONE /0.0D0,0.5D0,1.0D0/                            0005525
      DATA SWITCH /0/                                                   0005526
C                                                                       0005527
C ......................................................................0005528
C       .  INITIALIZATION OF PROJECTION PARAMETERS (FILE INPUT)  .      0005529
C ......................................................................0005530
C                                                                       0005531
      ENTRY IF20Z0 (INFILE,DATA)                                        0005532
C                                                                       0005533
      IERROR = 0                                                        0005534
      READ (INFILE,END=180) ZONE,BUFF                                   0005535
      IF (SWITCH.NE.0 .AND. SWITCH.EQ.ZONE) RETURN                      0005536
  020 MODE = 0                                                          0005537
      IF (BUFF(13) .NE. ZERO) MODE = 1                                  0005538
      IF (BUFF(1) .LE. ZERO) GO TO 100                                  0005539
      A = BUFF(1)                                                       0005540
      B = BUFF(2)                                                       0005541
      IF (B .GT. ZERO) GO TO 040                                        0005542
      E = ZERO                                                          0005543
      ES = ZERO                                                         0005544
      GO TO 120                                                         0005545
  040 IF (B .GT. ONE) GO TO 060                                         0005546
      E = DSQRT (B)                                                     0005547
      ES = B                                                            0005548
      GO TO 120                                                         0005549
  060 ES = ONE - (B / A) ** 2                                           0005550
      E = DSQRT (ES)                                                    0005551
      GO TO 120                                                         0005552
  100 A = AZ                                                            0005553
      E = EZ                                                            0005554
      ES = ESZ                                                          0005555
  120 KS0 = BUFF(3)                                                     0005556
      LAT0 = PAKRZ0 (BUFF(6))                                           0005557
      X0 = BUFF(7)                                                      0005558
      Y0 = BUFF(8)                                                      0005559
      SINPH0 = DSIN (LAT0)                                              0005560
      COSPH0 = DCOS (LAT0)                                              0005561
      CON = ONE - ES * SINPH0 * SINPH0                                  0005562
      COM = DSQRT (ONE - ES)                                            0005563
      BL = DSQRT (ONE + ES * COSPH0 ** 4 / (ONE - ES))                  0005564
      AL = A * BL * KS0 * COM / CON                                     0005565
      TS0 = TSFNZ0 (E,LAT0,SINPH0)                                      0005566
      CON = DSQRT (CON)                                                 0005567
      D = BL * COM / (COSPH0 * CON)                                     0005568
      F = D + DSIGN (DSQRT (DMAX1 ((D * D - ONE), 0.0D0)) , LAT0)       0005569
      EL = F * TS0 ** BL                                                0005570
      IF (IPPARM .EQ. 0) PRINT 2000, A,ES,KS0                           0005571
 2000 FORMAT (' INITIALIZATION PARAMETERS (OBLIQUE MERCATOR ''HOTINE''',0005572
     .        ' PROJECTION)'/                                           0005573
     .        ' SEMI-MAJOR AXIS OF ELLIPSOID =',F12.2,' METERS'/        0005574
     .        ' ECCENTRICITY SQUARED         =',F12.9/                  0005575
     .        ' SCALE AT CENTER              =',F12.9)                  0005576
      IF (MODE .EQ. 0) GO TO 140                                        0005577
      ALPHA = PAKRZ0 (BUFF(4))                                          0005578
      LONC = PAKRZ0 (BUFF(5))                                           0005579
      G = HALF * (F - ONE / F)                                          0005580
      GAMMA =  DASIN (DSIN (ALPHA) / D)                                 0005581
      LON0 = LONC -  DASIN (G * DTAN (GAMMA)) / BL                      0005582
C                                                                       0005583
C LIST INITIALIZATION PARAMETERS (CASE B).                              0005584
C                                                                       0005585
      CALL RADDZ0 (ALPHA,ANGS2(1,1))                                    0005586
      CALL RADDZ0 (LONC,ANGS2(1,2))                                     0005587
      CALL RADDZ0 (LAT0,ANGS2(1,3))                                     0005588
c     IF (IPPARM .EQ. 0) PRINT 2010, ANGS2                              0005589
      IF (IPPARM .EQ. 0) PRINT 2010, angs2(1,1),angs2(2,1),
     .   angs2(3,1),rangs1,angs2(1,2),angs2(2,2),angs2(3,2),rangs2,
     .   angs2(1,3),angs2(2,3),angs2(3,3),rangs3,X0,Y0
 2010 FORMAT (' AZIMUTH OF CENTRAL LINE      = ',A1,2I3,F7.3/           0005590
     .        ' LONGITUDE OF ORIGIN          = ',A1,2I3,F7.3/           0005591
     .        ' LATITUDE OF ORIGIN           = ',A1,2I3,F7.3)           0005592
      CON = DABS (LAT0)                                                 0005593
      IF (CON.GT.EPSLN .AND. DABS(CON - HALFPI).GT.EPSLN) GO TO 160     0005594
      IF (IPEMSG .EQ. 0) PRINT 2020                                     0005595
 2020 FORMAT (' ERROR PJ20Z0'/                                          0005596
     .        ' INPUT DATA ERROR')                                      0005597
      IERROR = 201                                                      0005598
      RETURN                                                            0005599
  140 LON1 = PAKRZ0 (BUFF(9))                                           0005600
      LAT1 = PAKRZ0 (BUFF(10))                                          0005601
      LON2 = PAKRZ0 (BUFF(11))                                          0005602
      LAT2 = PAKRZ0 (BUFF(12))                                          0005603
      SINPHI = DSIN (LAT1)                                              0005604
      TS1 = TSFNZ0 (E,LAT1,SINPHI)                                      0005605
      SINPHI = DSIN (LAT2)                                              0005606
      TS2 = TSFNZ0 (E,LAT2,SINPHI)                                      0005607
      H = TS1 ** BL                                                     0005608
      L = TS2 ** BL                                                     0005609
      F = EL / H                                                        0005610
      G = HALF * (F - ONE / F)                                          0005611
      J = (EL * EL - L * H) / (EL * EL + L * H)                         0005612
      P = (L - H) / (L + H)                                             0005613
      CALL RADDZ0 (LON2,ANGS1(1,3))                                     0005613
      DLON = LON1 - LON2                                                0005613
      IF (DLON .LT. -PI) LON2 = LON2 - 2.D0 * PI                        0005613
      IF (DLON .GT.  PI) LON2 = LON2 + 2.D0 * PI                        0005613
      DLON = LON1 - LON2                                                0005614
      LON0 = HALF * (LON1 + LON2) - DATAN (J * DTAN (HALF * BL *        0005615
     .       DLON) / P) / BL                                            0005616
      DLON = ADJLZ0 (LON1 - LON0)                                       0005617
      GAMMA = DATAN (DSIN (BL * DLON) / G)                              0005618
      ALPHA =  DASIN (D * DSIN (GAMMA))                                 0005619
      CALL RADDZ0 (LON1,ANGS1(1,1))                                     0005620
      CALL RADDZ0 (LAT1,ANGS1(1,2))                                     0005621
C     CALL RADDZ0 (LON2,ANGS1(1,3))                                     0005622
      CALL RADDZ0 (LAT2,ANGS1(1,4))                                     0005623
      CALL RADDZ0 (LAT0,ANGS1(1,5))                                     0005624
c     IF (IPPARM .EQ. 0) PRINT 2030, ANGS1                              0005625
      IF (IPPARM .EQ. 0) PRINT 2030,
     .   angs1(1,1),angs1(2,1),angs1(3,1),rangs11,
     .   angs1(1,2),angs1(2,2),angs1(3,2),rangs12,
     .   angs1(1,3),angs1(2,3),angs1(3,3),rangs13,
     .   angs1(1,4),angs1(2,4),angs1(3,4),rangs14,
     .   angs1(1,5),angs1(2,5),angs1(3,5),rangs15
 2030 FORMAT (' LONGITUDE OF 1ST POINT       = ',A1,2I3,F7.3/           0005626
     .        ' LATITUDE OF 1ST POINT        = ',A1,2I3,F7.3/           0005627
     .        ' LONGITUDE OF 2ND POINT       = ',A1,2I3,F7.3/           0005628
     .        ' LATITUDE OF 2ND POINT        = ',A1,2I3,F7.3/           0005629
     .        ' LATITUDE OF ORIGIN           = ',A1,2I3,F7.3)           0005630
      IF (DABS(LAT1 - LAT2) .LE. EPSLN) GO TO 150                       0005631
      CON = DABS (LAT1)                                                 0005632
      IF (CON.LE.EPSLN .OR. DABS(CON - HALFPI).LE.EPSLN) GO TO 150      0005633
      IF (DABS(DABS(LAT0) - HALFPI) .GT. EPSLN) GO TO 160               0005634
  150 IF (IPEMSG .EQ. 0) PRINT 2020                                     0005635
      IERROR = 202                                                      0005636
      RETURN                                                            0005637
  160 SINGAM = DSIN (GAMMA)                                             0005638
      COSGAM = DCOS (GAMMA)                                             0005639
      SINALF = DSIN (ALPHA)                                             0005640
      COSALF = DCOS (ALPHA)                                             0005641
      IF (IPEMSG .EQ. 0) PRINT 2040, X0,Y0                              0005642
 2040 FORMAT (' FALSE EASTING                =',F12.2,' METERS'/        0005643
     .        ' FALSE NORTHING               =',F12.2,' METERS')        0005644
      DATA(1) = A                                                       0005645
      DATA(2) = ES                                                      0005646
      SWITCH = ZONE                                                     0005647
      RETURN                                                            0005648
  180 IF (IPEMSG .EQ. 0) PRINT 2050                                     0005649
 2050 FORMAT (' ERROR PJ20Z0'/                                          0005650
     .        ' MISSING PROJECTION PARAMETERS')                         0005651
      IERROR = 203                                                      0005652
      RETURN                                                            0005653
C                                                                       0005654
C ......................................................................0005655
C      .  INITIALIZATION OF PROJECTION PARAMETERS (ENTRY INPUT)  .      0005656
C ......................................................................0005657
C                                                                       0005658
      ENTRY IS20Z0 (ZZONE,DATA)                                         0005659
	zone = zzone
C                                                                       0005660
      IERROR = 0                                                        0005661
      IF (SWITCH.NE.0 .AND. SWITCH.EQ.ZONE) RETURN                      0005662
      DO 200 I = 1,13                                                   0005663
      BUFF(I) = DATA(I)                                                 0005664
  200 CONTINUE                                                          0005665
      GO TO 020                                                         0005666
C                                                                       0005667
C ......................................................................0005668
C                      .  FORWARD TRANSFORMATION  .                     0005669
C ......................................................................0005670
C                                                                       0005671
      ENTRY PF20Z0 (GEOG,PROJ)                                          0005672
C                                                                       0005673
      IERROR = 0                                                        0005674
      IF (SWITCH .NE. 0) GO TO 220                                      0005675
      IF (IPEMSG .EQ. 0) PRINT 2050                                     0005676
      IERROR = 204                                                      0005677
      RETURN                                                            0005678
  220 SINPHI = DSIN (GEOG(2))                                           0005679
      DLON = ADJLZ0 (GEOG(1) - LON0)                                    0005680
      VL = DSIN (BL * DLON)                                             0005681
      IF (DABS(DABS(GEOG(2)) - HALFPI) .GT. EPSLN) GO TO 230            0005682
      UL = SINGAM * DSIGN (ONE , GEOG(2))                               0005683
      US = AL * GEOG(2) / BL                                            0005684
      GO TO 250                                                         0005685
  230 TS = TSFNZ0 (E,GEOG(2),SINPHI)                                    0005686
      Q = EL / TS ** BL                                                 0005687
      S = HALF * (Q - ONE / Q)                                          0005688
      T = HALF * (Q + ONE / Q)                                          0005689
      UL = (S * SINGAM - VL * COSGAM) / T                               0005690
      CON = DCOS (BL * DLON)                                            0005691
      IF (DABS(CON) .LT. TOL) GO TO 240                                 0005692
      US = AL * DATAN ((S * COSGAM + VL * SINGAM) / CON) / BL           0005693
      IF (CON .LT. ZERO) US = US + PI * AL / BL                         0005694
      GO TO 250                                                         0005695
  240 US = AL * BL * DLON                                               0005696
  250 IF (DABS(DABS(UL) - ONE) .GT. EPSLN) GO TO 255                    0005697
      IF (IPEMSG .EQ. 0) PRINT 2060                                     0005698
 2060 FORMAT (' ERROR PJ20Z0'/                                          0005699
     .        ' POINT PROJECTS INTO INFINITY')                          0005700
      IERROR = 205                                                      0005701
      RETURN                                                            0005702
  255 VS = HALF * AL * DLOG ((ONE - UL) / (ONE + UL)) / BL              0005703
  260 PROJ(1) = X0 + VS * COSALF + US * SINALF                          0005704
      PROJ(2) = Y0 + US * COSALF - VS * SINALF                          0005705
      RETURN                                                            0005706
C                                                                       0005707
C ......................................................................0005708
C                      .  INVERSE TRANSFORMATION  .                     0005709
C ......................................................................0005710
C                                                                       0005711
      ENTRY PI20Z0 (PROJ,GEOG)                                          0005712
C                                                                       0005713
      IERROR = 0                                                        0005714
      IF (SWITCH .NE. 0) GO TO 280                                      0005715
      IF (IPEMSG .EQ. 0) PRINT 2050                                     0005716
      IERROR = 206                                                      0005717
      RETURN                                                            0005718
  280 X = PROJ(1) - X0                                                  0005719
      Y = PROJ(2) - Y0                                                  0005720
      VS = X * COSALF - Y * SINALF                                      0005721
      US = Y * COSALF + X * SINALF                                      0005722
      Q = DEXP (- BL * VS / AL)                                         0005723
      S = HALF * (Q - ONE / Q)                                          0005724
      T = HALF * (Q + ONE / Q)                                          0005725
      VL = DSIN (BL * US / AL)                                          0005726
      UL = (VL * COSGAM + S * SINGAM) / T                               0005727
      IF (DABS (DABS (UL) - ONE) .GE. EPSLN) GO TO 300                  0005728
      GEOG(1) = LON0                                                    0005729
      GEOG(2) = DSIGN (HALFPI , UL)                                     0005730
      RETURN                                                            0005731
  300 CON = ONE / BL                                                    0005732
      TS = (EL / DSQRT ((ONE + UL) / (ONE - UL))) ** CON                0005733
      GEOG(2) = PHI2Z0 (E,TS)                                           0005734
      CON = DCOS (BL * US / AL)                                         0005735
      LON = LON0 - DATAN2 ((S * COSGAM - VL * SINGAM) , CON) / BL       0005736
      GEOG(1) = ADJLZ0 (LON)                                            0005737
      RETURN                                                            0005738
C                                                                       0005739
      END                                                               0005740

C **********************************************************************
C ** U.S.G.S. GENERAL MAP PROJECTION PACKAGE ...... EROS DATA CENTER **
C **          MATHEMATICAL ANALYSIS BY JOHN SNYDER                   **
C **                                                  JUNE,1982      **
C *********************************************************************
C                    *  SPACE OBLIQUE MERCATOR  * 
C *********************************************************************
C
      SUBROUTINE PJ21Z0
C
      IMPLICIT REAL*8 (A-Z)
      REAL RANGS1,RANGS2
      INTEGER*4 IERROR,IPEMSG,IPPARM,ANGS2
      INTEGER*4 SWITCH,I,ZONE,INFILE,N,L
	integer *4 zzone
      COMMON/PJ21/LON0,A,B,A2,A4,C1,C3,Q,T,U,W,XJ,P21,SA,CA,ES,S,START
      COMMON /ELLPZ0/ AZ,EZ,ESZ,E0Z,E1Z,E2Z,E3Z
C ********************** X0,Y0,GAMMA,LON0,AL,BL,EL *********************
      COMMON /ERRMZ0/ IERROR
      COMMON /PRINZ0/ IPEMSG,IPPARM 
c      COMMON /WORKZ0/ BUFF(15),ANGS2(4,2)
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
C
C ......................................................................
C       .  INITIALIZATION OF PROJECTION PARAMETERS (FILE INPUT)  .
C ......................................................................0
C
      ENTRY IF21Z0(INFILE,DATA)
C
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
 2000 FORMAT (' INITIALIZATION PARAMETERS (SPACE OBLIQUE MERCATOR', 
     .        ' PROJECTION)'/
     .        ' SEMI-MAJOR AXIS OF ELLIPSOID =',F12.2,' METERS'/
     .        ' ECCENTRICITY SQUARED         =',F12.9)
      ALF = PAKRZ0 (BUFF(4))
      LON0 = PAKRZ0 (BUFF(5))
      P21=BUFF(9)/1440.D0
      START=BUFF(11)
      CALL RADDZ0 (ALF,ANGS2(1,1))
      CALL RADDZ0 (LON0,ANGS2(1,2))
C     IF (IPPARM .EQ. 0) PRINT 2010, ANGS2,BUFF(9),BUFF(10),BUFF(11)
      IF (IPPARM .EQ. 0) PRINT 2010,
     .   ANGS2(1,1),ANGS2(2,1),ANGS2(3,1),RANGS1,
     .   ANGS2(1,2),ANGS2(2,2),ANGS2(3,2),RANGS2,
     .   BUFF(9),BUFF(10),BUFF(11)
 2010 FORMAT (' INCLINATION OF ORBIT           = ',A1,2I3,F7.3/ 
     .        ' LONGITUDE OF ASCENDING ORBIT   = ',A1,2I3,F7.3/
     .        ' PERIOD OF SATELLITE REVOLUTION = ',F15.10/
     .        ' LANDSAT RATIO                  = ',F15.10/
     .        ' LANDSAT END OF PATH FLAG       = ',F15.10)
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
 2050 FORMAT (' ERROR PJ21Z0'/
     .        ' MISSING PROJECTION PARAMETERS')
      IERROR = 211
      RETURN
C
C ......................................................................0
C      .  INITIALIZATION OF PROJECTION PARAMETERS (ENTRY INPUT)  .
C ......................................................................
C
      ENTRY IS21Z0(ZZONE,DATA)
      zone = zzone
C
      IERROR = 0
      IF (SWITCH.NE.0 .AND. SWITCH.EQ.ZONE) RETURN
      DO 200 I = 1,13
      BUFF(I) = DATA(I)
  200 CONTINUE
      GO TO 020
C
C ......................................................................
C                      .  FORWARD TRANSFORMATION  .
C ......................................................................0
C
      ENTRY PF21Z0 (GEOG,PROJ)
C
      IERROR = 0
      IF (SWITCH .NE. 0) GO TO 220
      IF (IPEMSG .EQ. 0) PRINT 2050
      IERROR = 212
      RETURN
  220 CONTINUE
      CONV=1.D-7
      DLAT=GEOG(2)
      DLON=GEOG(1)-LON0
C     TEST FOR LATITUDE AND LONGITUDE APPROACHING 90 DEGREES
C
      IF (DLAT.GT.1.570796D0) DLAT=1.570796D0
C
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
C
C     ADJUST FOR CONFUSION AT BEGINNING AND END OF LANDSAT ORBITS
C
  250 RLM=PI*BUFF(10)
      RLM2=RLM+2.0D0*PI
      N=N+1
      IF(N.GE.3) GO TO 300
      IF(TLAM.GT.RLM.AND.TLAM.LT.RLM2) GO TO 300
      IF(TLAM.LE.RLM)TLAMP=2.5D0*PI
      IF(TLAM.GE.RLM2) TLAMP=HALFPI
      GO TO 230
  260 IF (IPEMSG .EQ. 0) PRINT 50, SAV,TLAM,DLAT,DLON 
   50 FORMAT(1X,'  50 ITERATIONS WITHOUT CONV; SAV/TLAM ',5X,2F15.10/
     *2F25.10)
      IERROR = 214
      RETURN
  300 CONTINUE
C
C     TLAM COMPUTED - NOW COMPUTE TPHI
C
      DS=DSIN(TLAM)
      DD=DS*DS
      DP=DSIN(RADLT)
      TPHI=DASIN(((1.D0-ES)*CA*DP-SA*DCOS(RADLT)*DSIN(XLAMT))/DSQRT(1.D
     *0-ES*DP*DP))
C
C     COMPUTE X AND Y
C
      XTAN = (PI/4.0D0) + (TPHI/2.0D0)
      TANLG = DLOG(DTAN(XTAN))
      SD=DSIN(TLAM)
      SDSQ=SD*SD
      S=P21*SA*DCOS(TLAM)*DSQRT((1.D0+T*SDSQ)/((1.D0+W*SDSQ)*(1.D0+Q*SDS0
     *Q)))
      D=DSQRT(XJ*XJ+S*S)
      X=B*TLAM+A2*DSIN(2.D0*TLAM)+A4*DSIN(4.D0*TLAM)-TANLG*S/D
      PROJ(1)=A*X
      Y=C1*SD+C3*DSIN(3.D0*TLAM)+TANLG*XJ/D
      PROJ(2)=A*Y
      RETURN
C
C ......................................................................
C                      .  INVERSE TRANSFORMATION  .
C ......................................................................
C
      ENTRY PI21Z0 (PROJ,GEOG)
C
      IERROR = 0 
      IF (SWITCH .NE. 0) GO TO 280
      IF (IPEMSG .EQ. 0) PRINT 2050
      IERROR = 213
      RETURN
  280 CONTINUE
C
C     COMPUTES TRANSFORMED LAT/LON AND GEODETIC
C     LAT/LON GIVEN X-Y
C
C
C     BEGIN INVERSE COMPUTATION WITH APPROXIMATION FOR TLON. SOLVE
C     FOR TRANSFORMED LONG.
C
      X=PROJ(1)
      Y=PROJ(2)
      TLON= X/(A*B)
      CONV=1.D-9
      INUMB=0
  813 SAV=TLON
      SD=DSIN(TLON)
      SDSQ=SD*SD
      S=P21*SA*DCOS(TLON)*DSQRT((1.D0+T*SDSQ)/((1.D0+W*SDSQ)*(1.D0+Q*
     . SDSQ)))
      BLON=(X/A)+(Y/A)*S/XJ-A2*DSIN(2.D0*TLON)-A4*DSIN(4.D0*TLON)-(S/XJ)
     **(C1*DSIN(TLON)+C3*DSIN(3.D0*TLON)) 
      TLON=BLON/B
      DIF=TLON-SAV
      IF(DABS(DIF).LT.CONV) GO TO 814
      INUMB=INUMB+1
      IF(INUMB.LT.50) GO TO 813
      IERROR=215
      RETURN
C
C     COMPUTE TRANSFORMED LAT.
C
  814 CONTINUE
      ST=DSIN(TLON)
      DEFAC=DEXP(DSQRT(1.D0+S*S/XJ/XJ)*(Y/A-C1*ST-C3*DSIN(3.D0*TLON)))
      ACTAN=DATAN(DEFAC)
      TLAT=2.0D0*(ACTAN-(PI/4.0D0))
C
C     COMPUTE GEODETIC LONGITUDE
C
      DD=ST*ST
      IF(DABS(DCOS(TLON)).LT.1.D-7) TLON=TLON-1.D-7
      BIGK=DSIN(TLAT) 
      BIGK2=BIGK*BIGK
      XLAMT=DATAN(((1.D0-BIGK2/(1.D0-ES))*DTAN(TLON)*CA-BIGK*SA*DSQRT((1
     *.D0+Q*DD)*(1.D0-BIGK2)-BIGK2*U)/DCOS(TLON))/(1.D0-BIGK2*(1.D0+U)))
C
C     CORRECT INVERSE  QUADRANT
C
      IF(XLAMT.GE.0.D0) SL=1.D0
      IF(XLAMT.LT.0.D0) SL=-1.D0
      IF(DCOS(TLON).GE.0.D0) SCL=1.D0
      IF(DCOS(TLON).LT.0.D0) SCL=-1.0D0
      XLAMT=XLAMT-((PI/2.D0)*(1.D0-SCL)*SL)
      DLON=XLAMT-P21*TLON
C
C     COMPUTE GEODETIC LATITUDE
C
      IF(DABS(SA) .LT.1.D-7)DLAT=DASIN(BIGK/DSQRT((1.D0-ES)*(1.D0-ES) 
     *+ES*BIGK2))
      IF(DABS(SA) .LT.1.D-7)GO TO 1
      DLAT=DATAN((DTAN(TLON)*DCOS(XLAMT)-CA*DSIN(XLAMT))/((1.D0-ES)*SA))
    1 CONTINUE
      GEOG(1)=ADJLZ0(DLON+LON0)
      GEOG(2)=DLAT
      RETURN
      END

C                   QSFNZ0                                                     
C **********************************************************************0005742
C ** U.S.G.S. GENERAL MAP PROJECTION PACKAGE ...... DR. A. A. ELASSAL **0005743
C ** MODULE I                VERSION 1.0.0            NOVEMBER 1,1980 **0005744
C **********************************************************************0005745
      DOUBLE PRECISION FUNCTION QSFNZ0 (ECCENT,SINPHI,COSPHI)           0005746
C                                                                       0005747
C FUNCTION TO COMPUTE CONSTANT (SMALL Q).                               0005748
C                                                                       0005749
      IMPLICIT REAL*8 (A-Z)                                             0005750
      DATA HALF,ONE,TWO /0.5D0,1.0D0,2.0D0/                             0005751
      DATA EPSLN /1.0D-7/                                               0005752
C                                                                       0005753
      IF (ECCENT .LT. EPSLN) GO TO 020                                  0005754
      CON = ECCENT * SINPHI                                             0005755
      QSFNZ0 = (ONE - ECCENT * ECCENT) * (SINPHI / (ONE - CON * CON) -  0005756
     .         (HALF / ECCENT) * DLOG ((ONE - CON) / (ONE + CON)))      0005757
      RETURN                                                            0005758
C                                                                       0005759
  020 QSFNZ0 = TWO * SINPHI                                             0005760
      RETURN                                                            0005761
      END                                                               0005762

C                   RADDZ0                                                     
C **********************************************************************0005764
C ** U.S.G.S. GENERAL MAP PROJECTION PACKAGE ...... DR. A. A. ELASSAL **0005765
C ** MODULE I                VERSION 1.0.0            NOVEMBER 1,1980 **0005766
C **********************************************************************0005767
      SUBROUTINE RADDZ0 (RAD,IDMS)                                      0005768
C                                                                       0005769
C SUBROUTINE TO CONVERT ANGLE FROM RADIANS TO DMS.                      0005770
C SECONDS (IDMS) .                                                      0005771
C                                                                       0005772
      REAL*8 RAD,CON,RADSEC                                             0005773
      DIMENSION IDMS(1)                                                 0005774
      EQUIVALENCE (FLOT,INTG)                                           0005775
      DATA RADSEC,IBLANK,NEG /206264.80625D0,' ','-'/                   0005776
      DATA ZERO /0.0D0/                                                 0005777
C                                                                       0005778
C DETERMINE THE SIGN OF THE ANGLE.                                      0005779
C                                                                       0005780
      IDMS(1) = IBLANK                                                  0005781
      IF (RAD .LT. ZERO) IDMS(1) = NEG                                  0005782
C                                                                       0005783
C CONVERT THE ANGLE TO THOUSANDTH OF SECONDS.                           0005784
C                                                                       0005785
      CON = DABS(RAD) * RADSEC                                          0005786
      ISEC = CON                                                        0005787
C                                                                       0005788
C COMPUTE DEGREES PART OF THE ANGLE.                                    0005789
C                                                                       0005790
      INTG = ISEC / 3600                                                0005791
      IDMS(2) = INTG                                                    0005792
      ISEC = INTG * 3600                                                0005793
      CON = CON - DFLOAT(ISEC)                                          0005794
      ISEC = CON                                                        0005795
C                                                                       0005796
C COMPUTE MINUTES PART OF THE ANGLE.                                    0005797
C                                                                       0005798
      IDMS(3) = ISEC / 60                                               0005799
      ISEC = IDMS(3) * 60                                               0005800
      CON = CON - DFLOAT(ISEC)                                          0005801
C                                                                       0005802
C COMPUTE SECONDS PART OF THE ANGLE.                                    0005803
C                                                                       0005804
      FLOT = CON                                                        0005805
      IDMS(4) = INTG                                                    0005806
C                                                                       0005807
      RETURN                                                            0005808
      END                                                               0005809

C *********************************************************************
C ** U.S.G.S  GENERAL MAP PROJECTION PACKAGE ******* EROS DATA CENTER *
C ** MATHEMATICAL ANALYSIS BY JOHN SNYDER            JUNE,1982 ********
C ********************************************************************* 
C
C    SERIES TO CALCULATE A,B,C COEFFICIENTS TO CONVERT FROM
C    TRANSFORM LATITUDE,LONGITUDE TO SPACE OBLIQUE MERCATOR(SOM)
C    RECTANGULAR COORDINATES
C
C *********************************************************************
      SUBROUTINE SE21Z0(FB,FA2,FA4,FC1,FC3,DLAM)
      IMPLICIT REAL*8(A-Z)
      COMMON/PJ21/LON0,A,B,A2,A4,C1,C3,Q,T,U,W,XJ,P21,SA,CA,ES,S,START
C    CONVERT DLAM TO RADIANS
      DLAM=DLAM*0.0174532925D0
      SD=DSIN(DLAM) 
      SDSQ=SD*SD
      S=P21*SA*DCOS(DLAM)*DSQRT((1.D0+T*SDSQ)/((1.D0+W*SDSQ)*(1.D0+Q*SDS
     *Q)))
      H=DSQRT((1.D0+Q*SDSQ)/(1.D0+W*SDSQ))*(((1.D0+W*SDSQ)/((1.D0+Q*SDSQ
     *)*(1.D0+Q*SDSQ)))-P21*CA)
      SQ=DSQRT(XJ*XJ+S*S)
      FB=(H*XJ-S*S)/SQ
      FA2=FB*DCOS(2.D0*DLAM)
      FA4=FB*DCOS(4.D0*DLAM)
      FC=S*(H+XJ)/SQ
      FC1=FC*DCOS(DLAM)
      FC3=FC*DCOS(3.D0*DLAM)
      RETURN
      END

C****                                                              *****SPHD   
C ** U.S.G.S. GENERAL MAP PROJECTION PACKAGE ...... JOHN F. WAANANEN  **SPHD   
C ** MODULE I                VERSION 1.0.0             APRIL 2, 1981  **SPHD   
C ****                                                             *****SPHD   
C                                                                       SPHD   
      SUBROUTINE SPHDZ0(ISPH,PARM)                                      SPHD   
C     SUBROUTINE TO COMPUTE SPHEROID PARAMETERS                         SPHD   
C     SUBROUTINE SPHDZ0 IS REQUIRED FOR PROGRAMS NO. L176 AND NO. L177  SPHD   
C                                                                       SPHD   
C     ISPH IS THE SPHEROID CODE FROM THE FOLLOWING LIST:                SPHD  1
C     0 = CLARKE 1866           1 = CLARKE 1880                         SPHD  1
C     2 = BESSEL                3 = NEW INTERNATIONAL 1967              SPHD  1
C     4 = INTERNATIONAL 1909    5 = WGS 72                              SPHD  1
C     6 = EVEREST               7 = WGS 66                              SPHD  1
C     8 = GRS 1980              9 = AIRY                                SPHD  1
C    10 = MODIFIED EVEREST     11 = MODIFIED AIRY                       SPHD  1
C    12 = WALBECK              13 = SOUTHEAST ASIA                      SPHD  1
C    14 = AUSTRALIAN NATIONAL  15 = KRASSOVSKY                          SPHD  1
C    16 = HOUGH                17 = MERCURY 1960                        SPHD  1
C    18 = MODIFIED MERC 1968   19 = SPHERE OF RADIUS 6370997 M          SPHD  2
C                                                                       SPHD  2
C    PARM IS ARRAY OF PROJECTION PARAMETERS:                            SPHD  2
C       PARM(1) IS THE SEMI-MAJOR AXIS                                  SPHD  2
C       PARM(2) IS THE ECCENTRICITY SQUARED                             SPHD  2
C                                                                       SPHD  2
C     IF ISPH IS NEGATIVE, THE DEFAULT IS RESET FROM CLARKE 1866        SPHD  2
C     TO THE POSITIVE OF "ISPH" SPHEROID.                               SPHD  2
C                                                                       SPHD  2
C     IF ISPH = 0 , THE DEFAULT IS RESET TO CLARKE 1866                 SPHD  2
C                                                                       SPHD  3
C ****                                                             *****SPHD  3
C                                                                       SPHD  3
      IMPLICIT REAL*8 (A-H,O-Z)                                         SPHD  3
      DIMENSION PARM(2),AXIS(20),BXIS(20) !Alex, parm(15) to parm(2)   SPHD  3
C                                                                       SPHD  3
      COMMON/ELLPZ0/ AZ,EZ,ESZ,E0Z,E1Z,E2Z,E3Z                          SPHD  3
      COMMON/SPHRZ0/ AZZ                                                SPHD  3
      COMMON/ERRMZ0/ IERROR                                             SPHD  3
      COMMON/PRINZ0/ IPEMSG,IPPARM                                      SPHD  3
      COMMON/PROJZ0/ IPROJ                                              SPHD  4
C                                                                       SPHD  4
      DATA AXIS/6378206.4D0,6378249.145D0,6377397.155D0,6378157.5D0,    SPHD  4
     . 6378388.0D0,6378135.0D0,6377276.3452D0,6378145.0D0,6378137.0D0,  SPHD  4
     . 6377563.396D0,6377304.063D0,6377341.89D0,6376896.0D0,6378155.0D0,SPHD  4
     . 6378160.0D0,6378245.0D0,6378270.0D0,6378166.0D0,6378150.0D0,     SPHD  4
     . 6370997.0D0/                                                     SPHD  4
C                                                                       SPHD  4
      DATA BXIS/6356583.8D0,6356514.86955D0,6356078.96284D0,            SPHD  4
     . 6356772.2D0,6356911.94613D0,6356750.519915D0,6356075.4133D0,     SPHD  4
     . 6356759.769356D0,6356752.31414D0,6356256.91D0,6356103.039D0,     SPHD  5
     . 6356036.143D0,6355834.8467D0,6356773.3205D0,6356774.719D0,       SPHD  5
     . 6356863.0188D0,6356794.343479D0,6356784.283666D0,6356768.337303D0SPHD  5
     . ,6370997.0D0/                                                    SPHD  5
C                                                                       SPHD  5
      IF (ISPH.LT.0) GO TO 5                                            SPHD  5
      IF (PARM(1).NE.0.0D0.AND.IPROJ.NE.1) RETURN                       SPHD  5
    5 JSPH = IABS(ISPH) + 1                                             SPHD  5
      IF (JSPH.LE.20) GO TO 10                                          SPHD  5
      IERROR = 1211
      IF (IPEMSG.EQ.0) PRINT 1
    1 FORMAT(' Spheroid code of ',I5,' reset to 0')
      ISPH = 0                                                          SPHD  6
       JSPH = 1                                                         SPHD  6
   10 A = AXIS(JSPH)                                                    SPHD  6
      B = BXIS(JSPH)                                                    SPHD  6
      ES = (A*A-B*B)/(A*A)                                              SPHD  6
C     IF ISPH LE 0 THEN RESET DEFAULT                                   SPHD  6
      IF (ISPH.GT.0) GO TO 20                                           SPHD  6
      AZZ = 6370997.0D0                                                 SPHD  6
      EZ  = DSQRT(ES)                                                   SPHD  7
      E0Z = E0FNZ0(ES)                                                  SPHD  7
      E1Z = E1FNZ0(ES)                                                  SPHD  7
      E2Z = E2FNZ0(ES)                                                  SPHD  7
      E3Z = E3FNZ0(EZ)                                                  SPHD  7
      AZ  = A                                                           SPHD  7
      ESZ = ES                                                          SPHD  7
      IF (ES.EQ.0.0D0) AZZ=A                                            SPHD  7
      IF (IPPARM.NE.0) GO TO 20                                         SPHD  7
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
C                                                                       SPHD  8
   20 PARM(1) = A                                                       SPHD  9
      PARM(2) = ES                                                      SPHD  9
      RETURN                                                            SPHD  9
      END                                                               SPHD  9

C                   TSFNZ0                                                     
C **********************************************************************0005811
C ** U.S.G.S. GENERAL MAP PROJECTION PACKAGE ...... DR. A. A. ELASSAL **0005812
C ** MODULE I                VERSION 1.0.0            NOVEMBER 1,1980 **0005813
C **********************************************************************0005814
      DOUBLE PRECISION FUNCTION TSFNZ0 (ECCENT,PHI,SINPHI)              0005815
C                                                                       0005816
C FUNCTION TO COMPUTE CONSTANT (SMALL T).                               0005817
C                                                                       0005818
      IMPLICIT REAL*8 (A-Z)                                             0005819
      DATA HALF,ONE /0.5D0,1.0D0/                                       0005820
      DATA HALFPI /1.5707963267948966D0/                                0005821
C                                                                       0005822
      CON = ECCENT * SINPHI                                             0005823
      COM = HALF * ECCENT                                               0005824
      CON = ((ONE - CON) / (ONE + CON)) ** COM                          0005825
      TSFNZ0 = DTAN (HALF * (HALFPI - PHI)) / CON                       0005826
C                                                                       0005827
      RETURN                                                            0005828
      END                                                               0005829

C                   UNTFZ0                                                     
C **********************************************************************0005831
C ** U.S.G.S. GENERAL MAP PROJECTION PACKAGE ...... DR. A. A. ELASSAL **0005832
C ** MODULE I                VERSION 1.0.0            NOVEMBER 1,1980 **0005833
C **********************************************************************0005834
      SUBROUTINE UNTFZ0 (INUNIT,IOUNIT,FACTOR,IFLG)                     0005835
C                                                                       0005836
C SUBROUTINE TO DETERMINE CONVERGENCE FACTOR BETWEEN TWO LINEAL UNITS   0005837
C                                                                       0005838
C * INPUT ........                                                      0005839
C * INUNIT * UNIT CODE OF SOURCE.                                       0005840
C * IOUNIT * UNIT CODE OF TARGET.                                       0005841
C                                                                       0005842
C * OUTPUT .......                                                      0005843
C * FACTOR * CONVERGENCE FACTOR FROM SOURCE TO TARGET.                  0005844
C * IFLG   * RETURN FLAG = 0 , NORMAL RETURN.                           0005845
C                        = 1 , ABNORMAL RETURN.                         0005846
C                                                                       0005847
      IMPLICIT REAL*8 (A-H,O-Z)                                         0005848
      DIMENSION FACTRS(5,5)                                             0005849
      DATA ZERO,MAXUNT /0.0D0,5/                                        0005850
      DATA FACTRS /0.1000000000000000D01 , 0.0000000000000000D00 ,      0005851
     .             0.0000000000000000D00 , 0.2062648062470963D06 ,      0005852
     .             0.5729577951308231D02 ,                              0005853
     .             0.0000000000000000D00 , 0.1000000000000000D01 ,      0005854
     .             0.3048006096012192D00 , 0.0000000000000000D00 ,      0005855
     .             0.0000000000000000D00 ,                              0005856
     .             0.0000000000000000D00 , 0.3280833333333333D01 ,      0005857
     .             0.1000000000000000D01 , 0.0000000000000000D00 ,      0005858
     .             0.0000000000000000D00 ,                              0005859
     .             0.4848136811095360D-5 , 0.0000000000000000D00 ,      0005860
     .             0.0000000000000000D00 , 0.1000000000000000D01 ,      0005861
     .             0.2777777777777778D-3 ,                              0005862
     .             0.1745329251994330D-1 , 0.0000000000000000D00 ,      0005863
     .             0.0000000000000000D00 , 0.3600000000000000D04 ,      0005864
     .             0.1000000000000000D01 /                              0005865
C                                                                       0005866
      IF (INUNIT.LT.0 .OR. INUNIT.GE.MAXUNT) GO TO 020                  0005867
      IF (IOUNIT.GE.0 .AND. IOUNIT.LT.MAXUNT) GO TO 040                 0005868
  020 IFLG = 1                                                          0005869
      RETURN                                                            0005870
  040 FACTOR = FACTRS(IOUNIT+1 , INUNIT+1)                              0005871
      IF (FACTOR .EQ. ZERO) GO TO 020                                   0005872
      IFLG = 0                                                          0005873
      RETURN                                                            0005874
C                                                                       0005875
      END                                                               0005876
