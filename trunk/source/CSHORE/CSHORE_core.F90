!==================================================
module CSHORE_core
!==================================================
#include "CMS_cpp.h"
    implicit none

#ifdef CSHORE
contains
!*************************************************
    subroutine CSHORE_defaults()
! Sets the CSHORE default values
!*************************************************
    use const_def, only: pi
    use flow_def, only: grav
    use CSHORE_def
    implicit none
    
    !SOME OF THESE VALUES WILL BE UPDATED
    !IN THE CSHORE_INIT ROUTINE
    CTHETA = 1.0
    STHETA = 0.0
    SQR2 = sqrt(2.0)
    RBZERO = 0.1
    SQRG1 = sqrt(2.0/pi)
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
    GSGM1 = grav*SGM1
    GAMMA_CS = 0.80
    CSHORE_ON =.false.
    
    return
    endsubroutine CSHORE_defaults

!*************************************************
    subroutine cshore_cards(cardname,foundcard)
! Reads the cards related to CSHORE
!
! Author: Chris Reed, URS
! Modified by Alex Sanchez, USACE-CHL
!*************************************************
    use cshore_def
    implicit none
    integer :: ierr
    character(len=30) :: cardname,xcdum
    logical :: foundcard
    
    foundcard = .true.
    selectcase (cardname)
    case('CSHORE') !CWR
        call card_boolean(77,CSHORE_ON,ierr)
        
      case('CSHORE_RBZERO') !CWR
        backspace(77)
        read(77,*) cardname, RBZERO
        
      case('CSHORE_EFFF','SUSP_FRIC_EFF_COEF') !CWR
        backspace(77)
        read(77,*) cardname, EFFF
        
      case('CSHORE_EFFB','SUSP_BREAK_EFF_COEF') !CWR
        backspace(77)
        read(77,*) cardname, EFFB
        
      case('CSHORE_BLP') !CWR
        backspace(77)
        read(77,*) cardname, BLP
        
      case('CSHORE_SHIELD') !CWR
        backspace(77)
        read(77,*) cardname, SHIELD
        
      case('CSHORE_GSLMAX') !CWR
        backspace(77)
        read(77,*) cardname, GSLMAX
        
      case('CSHORE_TANPHI') !CWR
        backspace(77)
        read(77,*) cardname, TANPHI
        
      case('CSHORE_SLP') !CWR
        backspace(77)
        read(77,*) cardname, SLP
        
      case('CSHORE_GAMMA') !CWR
       backspace(77)
       read(77,*) cardname, GAMMA_CS
       
      case default
        foundcard = .false.
        
    endselect
    
    return
    endsubroutine CSHORE_cards

!**********************************************************
    subroutine cshore_init()
! Initializes the CSHORE variables
!**********************************************************
    use sed_def, only: rhosed,poros
    use size_def, only: ncellsD
    use prec_def
    use flow_def, only: grav,rhow
    use diag_def
    use cms_def, only: cmswave
    !use STRESS_VARS, only: waves
    use const_def, only: twopi
    use CSHORE_def
    implicit none
    
    if( .not. cmswave .and. CSHORE_ON) then
      call diag_print_error('You are using CSHORE routines which',&
       'requires a wave simulation, to continue')
    endif
    
    if(CSHORE_ON)then
      allocate(QCSx(ncellsD),QCSy(ncellsD),Calpha(ncellsD),Salpha(ncellsD),DTLoc(ncellsD))
      allocate(BSLOPE(ncellsD),GSLOPERAW(ncellsD),GSLOPE(ncellsD),ASLOPERAW(ncellsD))
      allocate(ASLOPE(ncellsD),ASLOPEX(ncellsD),ASLOPEY(ncellsD))
      allocate(CP(ncellsD),Re(ncellsD),RBETA(ncellsD),DBSTA(ncellsD),RQ(ncellsD),RQnew(ncellsD),RQf(ncellsD,2))
      SG = rhosed/rhow
      SPORO = 1.0/(1.0-poros)
      sporo1 = 1 - sporo
      sgm1 = SG - 1.0
      GSGM1 = grav*SGM1
      BLD = BLP/GSGM1
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
      aslopex = 0.0
      aslopey = 0.0
      dbsta = 0.0
      rbeta = 0.0
      re = 0.0
    endif
    
    !Note: This is also a good place to check variables and make sure they are consistent (Alex 10/24/14)
    
    return
    endsubroutine CSHORE_init

!*******************************************************************
    subroutine CSHORE_print()
! Prints the CSHORE settings to the screen and diagnostic file
! Author: Alex Sanchez, USACE-CHL
!*******************************************************************
    use CSHORE_def
    use diag_def
    implicit none
    integer :: i,iunit(2)
    
    iunit = (/6, dgunit/)
    open(dgunit,file=dgfile,access='append') 
    do i=1,2
      if(CSHORE_ON)then
        write(iunit(i),*)  'CSHORE                 ON'
      endif
    enddo
    
    return
    endsubroutine CSHORE_print
    
!!*************************************************
!    subroutine CSHORE_spatial()
!! Author: Brad Johnson, USACE-CHL
!! Modified by:
!!  Chris Reed, URS
!!  Alex Sanchez, USACE-CHL
!!*************************************************
!    use prec_def
!    use size_def, only: ncells,ncellsD
!    use geo_def, only: cell2cell,x,y,dx,dy,zb
!    use flow_def, only: grav,h,eta
!    use const_def, only: pi,twopi,deg2rad
!    use prec_def
!    use cms_def, only: cmswave
!    use wave_flowgrid_def, only: Wunitx,Wunity,Wper,Whgt,Wlen
!    use CSHORE_def
!    implicit none
!    real(ikind):: BSLOP1,BSLOP2,BSLOPEx,BSLOPEy,ABREAK,WKP,WHRMS,WT
!    real(ikind):: QBREAK,QBOLD,HM,B,ERFCC,B1
!    real(ikind):: flux1,flux2,sum1,sum2,term,dxt,dyt,dxww
!    integer:: ncn,ncs,nce,ncw,j,k,ii,jj,i,kk
!
!    BSLOP1 = -TANPHI*(GSLMAX-1.0_ikind)/GSLMAX
!    BSLOP2 = TANPHI*(GSLMAX+1.0_ikind)/(GSLMAX+2.0_ikind)
!
!! Calpha and Salpha are equivalent to Wunitx and Wunity (Alex 10/24/14)
!!    !need this for slope and roller flux calcs
!!    do i=1,ncells
!!      Calpha(i) = cos(Wang(I)*deg2rad)
!!      Salpha(i) = sin(Wang(I)*deg2rad)
!!    enddo
!
!    !store smoothed bathy in aslope (temporarily)
!    !call cshore_smooth(ncells,zb,aslope)
!    aslopex = dzbx; aslopey = dzby
!    call smooth_flowgrid_vec(aslopex,aslopey,1)
!    
!    !Calculate slope in wave direction - slope is positive when getting shallower
!    do i=1,ncells
!        BSLOPE(I) = aslopex(i)*Wunitx(i)+aslopey(i)*Wunity(i)!bed slope in wave direction, positive upslope
!    enddo    
!    do i=ncells+1,ncellsD
!      BSLOPE(I) = BSLOPE(cell2cell(1,i))
!    enddo
!    
!    do i=1,ncells
!      if(BSLOPE(i)<0.0_ikind) then
!        if(BSLOPE(i)>BSLOP1) then
!          GSLOPERAW(i) = TANPHI/(TANPHI + BSLOPE(i))
!        else
!          GSLOPERAW(i) = GSLMAX
!        endif
!      else
!        if(BSLOPE(i)<BSLOP2) then
!          GSLOPERAW(i) = (TANPHI - 2.0_ikind*BSLOPE(i))/(TANPHI-BSLOPE(i))
!        else
!          GSLOPERAW(i) = -GSLMAX
!        endif
!      endif
!      ASLOPERAW(i) = SLP
!      if(BSLOPE(i)>0.0_ikind) ASLOPERAW(i) = SLP + sqrt(BSLOPE(i)/TANPHI)
!    enddo
!
!    do i=ncells+1,ncellsD
!      jj = cell2cell(1,i)
!      GSLOPERAW(I) = GSLOPERAW(JJ)
!      ASLOPERAW(I) = ASLOPERAW(JJ)
!    enddo
!    GSLOPE = GSLOPERAW; ASLOPE = ASLOPERAW
!    call smooth_flowgrid_scal(GSLOPE,1)
!    call smooth_flowgrid_scal(ASLOPE,1)
!    !CALL CSHORE_SMOOTH(NCELLS, GSLOPERAW, GSLOPE)
!    !CALL CSHORE_SMOOTH(NCELLS, ASLOPERAW, ASLOPE)
!    
!    do J = 1,ncells
!      CP(J) = Wlen(J)/Wper(J)
!      RE(J)=CP(J)*CP(J)/grav
!      RBETA(J)=RBZERO
!      if(BSLOPE(J)>0.0_ikind) RBETA(J)=RBETA(J)+BSLOPE(J)*CTHETA
!    enddo
!    
!    !Solve for breaker dissipation
!    do J=1,ncells
!      D = h(i)
!      WKP = twopi/Wlen(J)
!      WT = Wper(J)
!      WHRMS = Whgt(J)/sqrt(2.0)
!      ABREAK = (TWOPI/WKP/D)*BSLOPE(J)*CTHETA/3.0_ikind
!      if(ABREAK<1.0_ikind) ABREAK = 1.0_ikind
!      HM = 0.88_ikind/WKP*tanh(GAMMA_CS*WKP*D/0.88_ikind)
!      B = (WHRMS/HM)**2.0_ikind
!      B1= WHRMS/HM
!      if(B<0.999990_ikind) then
!        QBOLD = B/2.0_ikind    
!10      QBREAK = QBOLD - (1.0_ikind-QBOLD + B*log(QBOLD))/(B/QBOLD-1.0_ikind)
!        if(QBREAK.LE.0.0_ikind) QBREAK = QBOLD/2.0_ikind
!        if(abs(QBREAK-QBOLD)>1.0e-6) then
!          QBOLD = QBREAK
!          GOTO 10
!        endif
!      else
!        QBREAK = 1.0_ikind
!        HM=WHRMS
!      endif
!      DBSTA(J) = 0.25_ikind*ABREAK*QBREAK*HM*HM/WT
!    enddo
!    
!    !solve for Rq (roller flux)
!    do K=1,10000  !Alex, this number seems excessive?
!      do i = 1,ncells
!        ncn = cell2cell(1,i)
!        nce = cell2cell(2,i)
!        ncs = cell2cell(3,i)
!        ncw = cell2cell(4,i)
!        Flux1 = (Re(i)*Calpha(i) + Re(ncw)*Calpha(ncw))/2.
!        DTLoc(i) = 0.25*dx(i)/(Flux1+1.e-10)
!        if(Flux1>0) then
!          Fuu(i) = Flux1*RQ(ncw)
!        else
!          Fuu(i) = Flux1*RQ(I)
!        endif
!        ncn = cell2cell(1,i)
!        Flux2 = (Re(i)*Salpha(i) + Re(ncs)*Salpha(ncs))/2.
!        if(Flux2>0) then
!          Gvv(i) = Gvv(i) + Flux2*RQ(ncs)
!        else
!          Gvv(i) = Gvv(i) + Flux2*RQ(i)
!        endif
!        DTLoc(i) = min(0.25*dx(i)/(abs(Flux1)+1.e-10),0.25*dy(i)/(abs(Flux2)+1.e-10))
!      enddo
!
!      do I=1,ncells
!        DXT = (dx(i)+dx(cell2cell(4,i)))/2.0
!        ADVECTX(i) = (Fuu(i)-Fuu(cell2cell(4,i)))/DXT
!        DYT = (DY(i)+DY(cell2cell(3,i)))/2.0
!        ADVECTY(i) = (Gvv(i)-Gvv(cell2cell(3,i)))/DYT
!      enddo
!      sum1 = 0.0
!      sum2 = 0.0
!      do i=1,ncells
!        term = abs(Rq(i)-Rqnew(i))
!        sum1 = sum1 + term
!        sum2 = max(sum2,term)
!        RQ(i) = Rqnew(i)
!      enddo
!    enddo
!
!    call diag_print_warning('wave roller flux solution did not converge',msg2)
!
!    return
!    endsubroutine CSHORE_spatial

!*****************************************************
    subroutine CSHORE_smooth(num,raw,F)
! Spatially smooths a variable based on a simple
! weighted average
!
! Input:
!   num - Number of cells to smooth
!   raw - Input variable
!
! Output:
!   F - smoothed variable
!
! Author: Brad Johnson, USACE-CHL
! Modified by:
!  Chris Reed, URS
!  Alex Sanchez, USACE-CHL
!
! Issue:
!  This routine needs to be modified for telescoping
!  and unstructured meshes (Alex Sanchez 10/23/14)
!*******************************************************
    use prec_def
    use geo_def, only: cell2cell
    USE CSHORE_def, only: calpha,salpha
    implicit none
    real(ikind) :: raw(1),f(1),c2,s2
    integer :: ncn,nce,ncs,ncw,num,j
    
    do j = 1,num
      ncn = cell2cell(1,j)
      nce = cell2cell(2,j)
      ncs = cell2cell(3,j)
      ncw = cell2cell(4,j)
      c2 = calpha(j)**2
      s2 = salpha(j)**2
      F(j)=(RAW(J)+s2*RAW(ncn)+c2*RAW(NCE)+c2*raw(NCW)+s2*raw(ncs))/3.0
    enddo
    
    return
    endsubroutine CSHORE_smooth

#endif
endmodule CSHORE_core