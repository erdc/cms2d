    !*************************************************
    !MODULE CSHORE_VARS
    !*************************************************
    MODULE CSHORE_VARS
    USE GLOBAL,ONLY: IKIND
    real(ikind), allocatable :: QCSx(:),QCSy(:),Calpha(:),Salpha(:),DTLoc(:)
    real(ikind), allocatable :: BSLOPE(:),GSLOPERAW(:),GSLOPE(:),ASLOPERAW(:),ASLO,PE(:)
    real(ikind), allocatable :: CP(:),Re(:),RBETA(:),DBSTA(:),RQ(:),RQnew(:),RQf(:,:)
    real(ikind):: CTHETA, STHETA, SQR2, RBZERO, SQRG1,SQRG2, EFFF, EFFB
    real(ikind):: sgm1,sporo1,gsgm1,blp,shield, bld
    real(ikind):: GSLMAX,TANPHI,SLP,GAMMA_CS
    real(ikind):: SG, SPORO, SQR8, FB2
    logical CSHORE_ON
    END MODULE

    !*************************************************
    !subroutine CSHORE_point
    !************************************************
    subroutine CSHORE_point (J,Hrms,WT,H,WL,d50,WF,QRAW,Qonsh,Qoffsh)
    USE Global, only: i,grav,ikind,szparms,depth
    USE TRANS_VARS, only:twopi
    USE CSHORE_VARS
    implicit none
    real(ikind):: KP,WHRMS,CPt,SIGMA,SIGSTA,SIGT,US,IGT,VSIGT,USTD,UMEAN
    real(ikind):: RB,RS,US,PB,PS,VS,QBX,QSX,QRAW,ER,FCC,WFSGM1,GSD50S
    real(ikind):: Per,d50,WF,HRMS,WL,WT,H,GBX,GF,DF,STA,USTA,SIGMA2
    real(ikind):: Qonsh,Qoffsh
    integer:: j
    !INPUTS
    !Hrms Wave Height (m)
    !WT - wave period (s)
    !H - water depth (m)
    !WL - Wave length (m)
    !d50 - grain size (m)
    !WF - fall speed (m/s)
    !OUTPTUS
    !Qraw - Net trnsport (m3/m/s)
    !Qonsh - Onshore directed transport (m3/m/s)
    !Qoffsh- offshore directed trans. (m3/m/s)
    WKP = ((twopi/WT)**2)/Grav
    WFSGM1 = WF*SGM1
    GSD50S = GSGM1*D50*SHIELD
    WHRMS = Hrms
    CPt = WL/WT
    SIGMA = HRMS/SQR8
    SIGSTA = SIGMA/H
    SIGT = SIGSTA*CPt
    SIGMA2 = SIGMA**2.D0
    USIGT = -SIGSTA*GRAV*H/CPt/CPt
    USIGT =USIGT*(1.D0+(CPt/GRAV)*RQ(J)/SIGMA2)
    VSIGT = 0.0
    CALL GBXAGF(CTHETA,USIGT,STHETA,VSIGT,GBX,GF)
    DFSTA = FB2*GF*SIGT**3.D0/GRAV
    USTD = SIGT*CTHETA
    UMEAN= -USTD*SIGSTA*GRAV*H/CPt/CPt
    USTA = UMEAN/SIGT
    RB = DSQRT(GSD50S/FB2)/USTD
    RS = WF/USTD/FB2**0.3333D0
    US = USTA
    PB=0.5D0*(ERFCC((RB+US)/SQR2)+ERFCC((RBUS)/SQR2))
    PS=0.5D0*(ERFCC((RS+US)/SQR2)+ERFCC((RSUS)/SQR2))
    IF(PS.GT.PB) PS = PB
    VS = PS*(EFFF*DFSTA + EFFB*RBETA(J)*RQ(J))/WFSGM1
    VS = VS*DSQRT(1.D0+BSLOPE(J)*BSLOPE(J))
    QBX = BLD*PB*GSLOPE(J)*USTD**3
    QSX = ASLOPE(J)*UMEAN*VS
    Qonsh = QBX
    Qoffsh = QSX
    QRAW = (QBX + QSX) !/SPORO1
    end subroutine

    !*************************************************
    !subroutine CSHORE_spatial
    !*************************************************
    subroutine CSHORE_spatial()
    use global, only: ikind,location,ncells,pi,dx,&
      dy,ncellsD,depth, linktodummies,eta,grav,fuu,&
      gvv,advectx,advecty,x,y,szparms,rhow
    use TRANS_VARS, only: twopi
    use STRESS_VARS, only: Wdir,Wper,Whgt,Wlen,waves,wdiss
    use CSHORE_VARS
    implicit none
    real(ikind):: BSLOP1,BSLOP2,BSLOPEx,BSLOPEy,ABREAK,WKP,WHRMS,WT
    real(ikind):: QBREAK,QBOLD,H,HM,D,B,ERFCC,B1
    real(ikind):: flux1,flux2,sum1,sum2,term,dxt,dyt,dxww
    integer:: ncn,ncs,nce,ncw,j,k,ii,jj,i,kk

    BSLOP1 = -TANPHI*(GSLMAX-1.D0)/GSLMAX
    BSLOP2 = TANPHI*(GSLMAX+1.D0)/(GSLMAX+2.D0)
    !need this for slope and roller flux calcs
    do i=1,ncells
      Calpha(i) = cos(Wdir(I)*pi/180.)
      Salpha(i) = sin(Wdir(I)*pi/180.)
    enddo

    !store smoothed bathy in aslope (temporarily)
    call cshore_smooth(ncells,depth,aslope)
    ii=0
    do i=ncells+1,ncellsD
      ii=ii+1
      jj = linktodummies(ii)
      aslope(i) = aslope(jj)
    enddo
    
    !calcualte slope in wave direction - slope is positive when getting shallower
    DO I=1,NCELLS
      ncn = location(i,1)
      nce = location(i,2)
      ncs = location(i,3)
      ncw = location(i,4)
      !if(Calpha(i).ge.0) then
      BSLOPEx = 2*(aslope(ncw)-aslope(nce))/(DX(NCW)+2*DX(I)+DX(NCE))
      !else
      !BSLOPEx = 2*(aslope(ncw)-aslope(nce))/(DX(NCW)+DX(i))
      !endif
      !IF(Salpha(i).ge.0) then
      BSLOPEy = 2*(aslope(ncs)-aslope(ncn))/(DY(NCS)+2*DY(I)+DY(NCN))
      !else
      !BSLOPEy = 2*(aslope(ncs)-aslope(ncn))/(DY(NCS)+DY(I))
      !endif
      BSLOPE(I) = BSLOPEx*Calpha(i)+BSLOPEy*Salpha(i)
    enddo
    
    ii=0
    do i=ncells+1,ncellsD
      ii=ii+1
      jj = linktodummies(ii)
      BSLOPE(I) = BSLOPE(JJ)
    enddo
    
    DO 100 J=1,NCELLS
      IF(BSLOPE(J).LT.0.D0) THEN
        IF(BSLOPE(J).GT.BSLOP1) THEN
          GSLOPERAW(J) = TANPHI/(TANPHI + BSLOPE(J))
        ELSE
          GSLOPERAW(J) = GSLMAX
        ENDIF
      ELSE
        IF(BSLOPE(J).LT.BSLOP2) THEN
          GSLOPERAW(J) = (TANPHI - 2.D0*BSLOPE(J))/(TANPHI-BSLOPE(J))
        ELSE
          GSLOPERAW(J) = -GSLMAX
        ENDIF
      ENDIF
      ASLOPERAW(J) = SLP
      IF(BSLOPE(J).GT.0.D0) ASLOPERAW(J) = SLP + DSQRT(BSLOPE(J)/TANPHI)
100 CONTINUE

    ii=0
    do i=ncells+1,ncellsD
      ii=ii+1
      jj = linktodummies(ii)
      GSLOPERAW(I) = GSLOPERAW(JJ)
      ASLOPERAW(I) = ASLOPERAW(JJ)
    enddo
    CALL CSHORE_SMOOTH(NCELLS, GSLOPERAW, GSLOPE)
    CALL CSHORE_SMOOTH(NCELLS, ASLOPERAW, ASLOPE)
    DO J = 1,ncells
      CP(J) = Wlen(J)/Wper(J)
      RE(J)=CP(J)*CP(J)/GRAV
      RBETA(J)=RBZERO
      IF(BSLOPE(J).GT.0.D0) RBETA(J)=RBETA(J)+BSLOPE(J)*CTHETA
    ENDDO
    
    !Solve for breaker dissipation
    DO J=1,ncells
      H = eta(J)+depth(j)
      D = H
      WKP = twopi/Wlen(J)
      WT = Wper(J)
      WHRMS = Whgt(J)/sqrt(2.0)
      ABREAK = (TWOPI/WKP/H)*BSLOPE(J)*CTHETA/3.D0
      IF(ABREAK.LT.1.D0) ABREAK = 1.D0
      HM = 0.88D0/WKP*DTANH(GAMMA_CS*WKP*D/0.88D0)
      B = (WHRMS/HM)**2.D0
      B1= WHRMS/HM
      IF(B.LT.0.99999D0) THEN
        QBOLD = B/2.D0    
10      QBREAK = QBOLD - (1.D0-QBOLD + B*DLOG(QBOLD))/(B/QBOLD-1.D0)
        IF(QBREAK.LE.0.D0) QBREAK = QBOLD/2.D0
        IF(DABS(QBREAK-QBOLD).GT.1.D-6) THEN
          QBOLD = QBREAK
          GOTO 10
        ENDIF
      ELSE
        QBREAK = 1.D0
        HM=WHRMS
      ENDIF
      DBSTA(J) = 0.25D0*ABREAK*QBREAK*HM*HM/WT
    enddo
    
    !solve for Rq (roller flux)
    DO K=1,10000
      DO I = 1,NCELLS
        ncn = location(i,1)
        nce = location(i,2)
        ncs = location(i,3)
        ncw = location(i,4)
        Flux1 = (Re(i)*Calpha(i) + Re(ncw)*Calpha(ncw))/2.
        DTLoc(i) = 0.25*dx(i)/(Flux1+1.e-10)
        if(Flux1.gt.0) then
          Fuu(i) = Flux1*RQ(ncw)
        else
          Fuu(i) = Flux1*RQ(I)
        endif
        ncn = location(i,1)
        Flux2 = (Re(i)*Salpha(i) + Re(ncs)*Salpha(ncs))/2.
        if(Flux2.gt.0) then
          Gvv(i) = Gvv(i) + Flux2*RQ(ncs)
        else
          Gvv(i) = Gvv(i) + Flux2*RQ(i)
        endif
        DTLoc(i) = min(0.25*dx(i)/(abs(Flux1)+1.e-10),0.25*dy(i)/(abs(Flux2)+1.e-10))
      enddo

      do I=1,ncells
        DXT = (dx(i)+dx(location(i,4)))/2.0
        ADVECTX(i) = (Fuu(i)-Fuu(location(i,4)))/DXT
        DYT = (DY(i)+DY(location(i,3)))/2.0
        ADVECTY(i) = (Gvv(i)-Gvv(location(i,3)))/DYT
      enddo
      sum1 = 0.0
      sum2 = 0.0
      DO i=1,ncells
        term = abs(Rq(i)-Rqnew(i))
        sum1 = sum1 + term
        sum2 = max(sum2,term)
        RQ(i) = Rqnew(i)
      enddo
    enddo
    write(*,*)'wave roller flux solution did not converge'
    write(*,*)'sum1,sum2 =',sum1,sum2
!    enddo
    end

    !*************************************************
    subroutine ERFCC
    !*************************************************
    !FUNCTION ERFCC(X)
    use global, only: ikind
    real(ikind) X, Z, T, ERFCC
    Z=DABS(X)
    T=1.D0/(1.D0+0.5D0*Z)
    ERFCC=T*DEXP(-Z*Z-1.26551223D0+T*(1.00002368D0+T*(.37409196D0+&
      T*(.09678418D0+T*(-.18628806D0+T*(.27886807D0+&
      T*(-1.13520398D0+T*(1.48851587D0+&
      T*(-.82215223D0+T*.17087277D0)))))))))
    IF (X.LT.0.D0) ERFCC=2.D0-ERFCC
    RETURN
    END

    !*************************************************
    !subroutine CSHORE_SMOOTH
    !*************************************************
    SUBROUTINE CSHORE_SMOOTH(NUM,RAW,F)
    use global, only: location,ikind
    USE CSHORE_VARS, only: calpha,salpha
    implicit none
    real(ikind) raw(1),f(1),c2,s2
    integer ncn,nce,ncs,ncw,num,j
    DO 201 j = 1,NUM
      ncn = location(j,1)
      nce = location(j,2)
      ncs = location(j,3)
      ncw = location(j,4)
      c2 = calpha(j)**2
      s2 = salpha(j)**2
      F(j)=(RAW(J)+s2*RAW(ncn)+c2*RAW(NCE)+c2*raw(NCW)+s2*raw(ncs))/3.0
201 CONTINUE
    RETURN
    END

    !*************************************************
    !subroutine CSHORE_Defualts
    !*************************************************
    SUBROUTINE CSHORE_DEFAULTS()
    use global, only: PI,GRAV
    use cshore_vars
    implicit none
    !SOME OF THESE VALUES WILL BE UPDATED
    !IN THE CSHORE_INIT ROUTINE
    CTHETA = 1.0
    STHETA = 0.0
    SQR2 = sqrt(2.0)
    RBZERO = 0.1
    SQRG1 = sqrt(2.0/PI)
    SQRG2 = 2*SQRG1
    EFFF = 0.01
    EFFB = 0.005
    SG = 2.6
    SPORO = 0.4
    sporo1 = 1 - sporo
    SQR8 = sqrt(8.)
    blp = 0.002
    shield = 0.05
    GSLMAX = 10.0
    TANPHI = 0.63
    SLP = 0.2
    sgm1 = SG - 1.0
    GSGM1 = GRAV*SGM1
    GAMMA_CS = 0.80
    CSHORE_ON =.false.
    RETURN
    END

    !*************************************************
    subroutine cshore_cards(cardname,foundcard)
    !*************************************************
    use cshore_vars
    implicit none
    character*30 :: cardname,xcdum
    logical :: foundcard
    
    foundcard = .true.
    selectcase (cardname)
      case('CSHORE') !CWR
        !backspace(1)
        read(1,*) cardname, xcdum
        if(xcdum.eq.'ON ') CSHORE_ON = .true.
        write(*,*)'CSHORE WILL BE SIMULATED'
      case('CSHORE_RBZERO') !CWR
        !backspace(1)
        read(1,*) cardname, RBZERO
      case('CSHORE_EFFF') !CWR
        !backspace(1)
        read(1,*) cardname, EFFF
      case('CSHORE_EFFB') !CWR
        !backspace(1)
        read(1,*) cardname, EFFB
      case('CSHORE_BLP') !CWR
        !backspace(1)
        read(1,*) cardname, BLP
      case('CSHORE_SHIELD') !CWR
        !backspace(1)
        read(1,*) cardname, SHIELD
      case('CSHORE_GSLMAX') !CWR
        !backspace(1)
        read(1,*) cardname, GSLMAX
      case('CSHORE_TANPHI') !CWR
        !backspace(1)
        read(1,*) cardname, TANPHI
      case('CSHORE_SLP') !CWR
        !backspace(1)
        read(1,*) cardname, SLP
      case('CSHORE_GAMMA') !CWR
       !backspace(1)
       read(1,*) cardname, GAMMA_CS
      case default
        foundcard = .false.
    endselect
    
    return
    endsubroutine

    !*************************************************
    subroutine cshore_init()
    !*************************************************
    use global, only: ikind,ncellsd,grav,rhosed,rhow,poros
    use STRESS_VARS, only: waves
    use TRANS_VARS, only: twopi
    use cshore_vars
    implicit none
    if( .not. waves .and. CSHORE_ON) then
      write(*,*)"You are using CSHORE routines which"
      write(*,*)"require a wave simulation, to continue"
      write(*,*)"invoke a wave sim. - program stopped"
    endif
    if(CSHORE_ON)THEN
      allocate(QCSx(ncellsd),QCSy(ncellsd),Calpha(ncellsd),Salpha(ncellsd),DTLoc(ncellsd))
      allocate(BSLOPE(ncellsd),GSLOPERAW(ncellsd),GSLOPE(ncellsd),ASLOPERAW(ncellsd),ASLOPE(ncellsd))
      allocate(CP(ncellsd),Re(ncellsd),RBETA(ncellsd),DBSTA(ncellsd),RQ(ncellsd),RQnew(ncellsd),RQf(ncellsd,2))
      SG = rhosed/rhow
      SPORO = 1.0/(1.0-poros)
      sporo1 = 1 - sporo
      sgm1 = SG - 1.0
      GSGM1 = GRAV*SGM1
      BLD = BLP/GSGM1
      twopi = 8.0*atan(1.0)
      FB2 = 0.005
      qcsx = 0.0
      qcsy = 0.0
      Rq = 0.0
      rqnew = 0.0
      bslope = 0.0
      gsloperaw = 0.0
      gslope = 0.0
      asloperaw = 0.0
      aslope = 0.0
      dbsta = 0.0
      rbeta = 0.0
      re = 0.0
    ENDIF
    RETURN
    END SUBROUTINE

    !*************************************************
    !subroutine GBXAGF
    !*************************************************
    ! This subroutine computes GBX and GF for specified CTHETA, USIGT,
    ! STHETA and VSIGT for Gaussian variable R
    SUBROUTINE GBXAGF(CTHETA,USIGT,STHETA,VSIGT,GBX,GF)
    use global, only: ikind
    use CSHORE_vars, only: sqr2,sqrg1
    implicit none
    real(ikind) :: C1,C2,C3,GBX,GF,H,HM,B,WKP,WT,ERFCC
    real(ikind) :: CTHETA,USIGT,STHETA,VSIGT
    ! For normally incident waves, use analytical
    ! expresions involving complementary error function ERFCC below
    C1 = 1.D0-ERFCC(USIGT/SQR2)
    C2 = SQRG1*DEXP(-USIGT*USIGT/2.D0)
    C3 = 1.D0 + USIGT*USIGT
    GBX = C3*C1 + C2*USIGT
    GF = USIGT*(C3 + 2.D0)*C1 + (C3 + 1.D0)*C2
    RETURN
    END