!************************************************************************    
    subroutine wave_default()
!************************************************************************    
    use wave_flowgrid_def
    implicit none
        
    !--- Constant wave parameters for testing and idealized cases --------
    constant_waves = .false.
    waveheight = 0.08
    waveperiod = 1.5
    wavedir = 0.0
    
    return
    endsubroutine wave_default
    
!***********************************************************************    
    subroutine wave_cards(cardname,foundcard)
!***********************************************************************    
    use wave_flowgrid_def
    use cms_def, only:  ignore_brk_restr
    use in_lib
    implicit none
    integer :: ierr
    character(len=*) :: cardname   
    logical :: foundcard
    
    foundcard = .true.
    selectcase(cardname)
    case('CONSTANT_WAVE_HEIGHT','WAVE_HEIGHT_CONSTANT','SIGNIFICANT_WAVE_HEIGHT_CONSTANT')
      call card_scalar(77,'m','m',waveheight,ierr)
      constant_waves = .true.
          
    case('CONSTANT_WAVE_PERIOD','WAVE_PERIOD_CONSTANT','PEAK_WAVE_PERIOD_CONSTANT')
      call card_scalar(77,'s','s',waveperiod,ierr)
      constant_waves = .true.
        
    case('CONSTANT_WAVE_DIRECTION','WAVE_DIRECTION_CONSTANT','MEAN_WAVE_DIRECTION_CONSTANT')
      call card_scalar(77,'deg','deg',wavedir,ierr) 
      constant_waves = .true.
    
    case('IGNORE_BREAKING_RESTRICTION')
      call card_boolean(77,ignore_brk_restr,ierr)

    case default 
      foundcard = .false.
      
    endselect
    
    return
    endsubroutine wave_cards
    
!***********************************************************************
    subroutine wave_init()
! Initializes the wave variables
! Author: Alex Sanchez, USACE-CHL
!***********************************************************************
#include "CMS_cpp.h"    
    use size_def, only: ncellsD
    use geo_def, only: grdfile,wgrdfile,azimuth_fl
    use cms_def
    use flow_def, only: iwet,rhow
    use comvarbl, only: timehrs,ctime,ramp,rampdur
    use const_def, only: deg2rad
    use diag_def
    use diag_lib
#ifdef XMDF_IO
    use in_xmdf_lib, only: readscalsteph5,readvecsteph5
#endif
    use wave_flowgrid_def
    implicit none
    integer :: i,ii,ierr
#ifdef DIAG_MODE
    logical :: isnankind
#endif
    
    if(constant_waves)then
      call wave_const_init
    endif

    if(noptset==3)then !Wave-Flow coupling
      nsteer=nsteer+1
      open(dgunit,file=dgfile,access='append') 
      call CMS_Wave_inline                      !call CMS_Wave   !modified MEB 10/15/2018
      close(dgunit)
      call steer_init !Needs to be called after CMS-Wave
      call interp_coef_flwav
      call interp_coef_wavfl
      !if(n2Dor3D==3) call allocate_wavestress3D2   !For 3D
      call rol_init      
      call freememory_fl_wav 
      call getwave        
      !if(n2Dor3D==3) call getwave3D    !For 3D
      call diag_print_message('*** Finished CMS-Wave Run ***',' ',&
            '*** Starting CMS-Flow Run ***' )      
      timehrs=ctime/3600.0
      tswave1=ctime
      tswave2=tswave1
      ramp=min(timehrs/rampdur,1.0)
      call tidevalue(tswave2,tide2)
      tide1=tide2
      do ii=1,ncellsD
        wavestrx1(ii)=wavestrx2(ii)
        wavestry1(ii)=wavestry2(ii)            
        Whgt1(ii)=Whgt2(ii)            
        Wper1(ii)=Wper2(ii)            
        Wunitx1(ii)=Wunitx2(ii)
        Wunity1(ii)=Wunity2(ii)
        !waveibr1(ii)=waveibr2(ii)
        wavediss1(ii)=wavediss2(ii) !Alex, Aug 3, 2009
        !Needed here to initialize time zero
        wavestrx(ii)=wavestrx2(ii)*ramp
        wavestry(ii)=wavestry2(ii)*ramp
        Whgt(ii)=Whgt2(ii)*sqrt(ramp) !because E~Hs^2
        Wper(ii)=Wper2(ii)          
        Wunitx(ii)=Wunitx2(ii)
        Wunity(ii)=Wunity2(ii)
        !waveibr(ii)=waveibr2(ii)
        wavediss(ii)=wavediss2(ii)*ramp
      enddo
      !do ii=1,ncells3DD
      !   wavestrx3D1(ii)=wavestrx3D2(ii)
      !   wavestry3D1(ii)=wavestry3D2(ii)
      !   wavestrx3D(ii)=wavestrx3D2(ii)*ramp
      !   wavestry3D(ii)=wavestry3D2(ii)*ramp
      ! enddo
    elseif(noptset==4)then !Constant wave conditions
      timehrs=ctime/3600.0
      ramp=min(timehrs/rampdur,1.0)
      
      nsteer=nsteer+1
      call wave_flgrid_init
#ifdef XMDF_IO
      call readscalsteph5(wgrdfile,wavpath, nsteer,tswave1,Whgt1,    ierr)           !Updated with 'wgrdfile' to use since there is no _grid.h5 file anymore.  MEB  06/10/2021
      call readscalsteph5(wgrdfile,perpath, nsteer,tswave1,Wper1,    ierr)              
      call readscalsteph5(wgrdfile,dirpath, nsteer,tswave1,Wang,     ierr)  
      call readscalsteph5(wgrdfile,disspath,nsteer,tswave1,wavediss1,ierr)  
      call readvecsteph5 (wgrdfile,radpath, nsteer,tswave1,wavestrx1,wavestry1,ierr)
#else
      call diag_print_error('Cannot read constant wave conditions without XMDF')
#endif
      do i=1,ncellsD
        Wunitx1(i) = cos((Wang(i)-azimuth_fl)*deg2rad)
        Wunity1(i) = sin((Wang(i)-azimuth_fl)*deg2rad)
        wavediss1(i) = -rhow*wavediss1(i)  !Flip sign and convert units [N/m/s]
        !if(wavediss1(i)>wavedisstol)then
        !  waveibr1(i) = 1.0
        !endif
      enddo
!      call rotate_vector(ncellsD,ncellsD,azimuth_fl,wavestrx1,wavestry1)      
!      call rotate_vector(ncellsD,ncellsD,azimuth_fl,Wunitx1,Wunity1) 
      tswave2 = max(tswave1,ctime)
      tswave1 = min(tswave1,ctime)  
      
      !Smoothing
      !do i=1,nbrksm
      !  where(iwet==0)
      !    waveibr1 = 1.0 !Set all dry cells as breaking
      !  endwhere  
      !  call smooth_flowgrid_scal(waveibr1,1)
      !enddo
      if(ndissm>0) call smooth_flowgrid_scal(wavediss1,ndissm)
      if(nradsm>0) call smooth_flowgrid_vec (wavestrx1,wavestry1,nradsm)
      if(npersm>0) call smooth_flowgrid_scal(Wper1,npersm) !Wave period smoothing
      do i=1,ncellsD
        if(wavediss1(i)<wavedisstol)then
          wavediss1(i) = 0.0
          !waveibr1(i) = 0.0
          !wavestrx1(i) = 0.0
          !wavestry1(i) = 0.0
        !else
        !  waveibr1(i) = 1.0
        endif
        !wavestrx1(i) = wavestrx1(i)*waveibr1(i) !Limit radiation stresses to breaking zone
        !wavestry1(i) = wavestry1(i)*waveibr1(i) !Limit radiation stresses to breaking zone
#ifdef DIAG_MODE
        if(isnankind(wavestrx1(i)) .or. isnankind(wavestry1(i)))then
          continue
        endif
#endif        
        !Needed here to initialize time zero
        Whgt(i)=Whgt1(i)*sqrt(ramp) !because E~Hs^2
        Wper(i)=Wper1(i)          
        Wunitx(i)=Wunitx1(i)
        Wunity(i)=Wunity1(i)
        wavestrx(i)=wavestrx1(i)*ramp
        wavestry(i)=wavestry1(i)*ramp
        wavediss(i)=wavediss1(i)*ramp
        !waveibr(i)=waveibr1(i)
      enddo
      
      nsteer=nsteer+1
#ifdef XMDF_IO
      !Updated with 'wgrdfile' to use since there is no _grid.h5 file anymore.  MEB  06/10/2021
      call readscalsteph5(wgrdfile,wavpath,nsteer,tswave2,Whgt2,ierr)      !ierr returns: -2 if File can't open, 3 if can't read timestep, 4 if can't read data record
#endif
      select case(ierr)
      case(:-1,1:)    !if(ierr<0 .or. ierr>0)then
        do i=1,ncellsD
          wavestrx2(i)=wavestrx1(i)
          wavestry2(i)=wavestry1(i)
          Whgt2(i)=Whgt1(i)
          Wper2(i)=Wper1(i)
          wavediss2(i)=wavediss1(i)
          !waveibr2(i)=waveibr1(i)
          Wunitx2(i)=Wunitx1(i)
          Wunity2(i)=Wunity1(i)
        enddo  
      case default
#ifdef XMDF_IO
        call readscalsteph5(wgrdfile,perpath, nsteer,tswave2,Wper2,              ierr) !Updated with 'wgrdfile' to use since there is no _grid.h5 file anymore.  MEB  06/10/2021
        call readscalsteph5(wgrdfile,dirpath, nsteer,tswave2,Wang,               ierr)
        call readscalsteph5(wgrdfile,disspath,nsteer,tswave2,wavediss2,          ierr)
        call readvecsteph5 (wgrdfile,radpath, nsteer,tswave2,wavestrx2,wavestry2,ierr)
#endif
        do i=1,ncellsD
          Wunitx2(i)=cos((Wang(i)-azimuth_fl)*deg2rad)
          Wunity2(i)=sin((Wang(i)-azimuth_fl)*deg2rad)
          wavediss2(i) = -rhow*wavediss2(i)  !Flip sign and convert units [N/m/s]
          !if(wavediss2(i)>wavedisstol)then
          !  waveibr2(i) = 1.0
          !endif
        enddo 
!        call rotate_vector(ncellsD,ncellsD,azimuth_fl,wavestrx2,wavestry2)   
!        call rotate_vector(ncellsD,ncellsD,azimuth_fl,Wunitx2,Wunity2)   
        
        !Smoothing
        !do i=1,nbrksm
        !  where(iwet==0)
        !    waveibr2=1.0 !Set all dry cells as breaking
        !  endwhere  
        !  call smooth_flowgrid_scal(waveibr2,1)
        !enddo
        if(ndissm>0) call smooth_flowgrid_scal(wavediss2,ndissm)
        if(nradsm>0) call smooth_flowgrid_vec (wavestrx2,wavestry2,nradsm)
        if(npersm>0) call smooth_flowgrid_scal(Wper2,npersm) !Wave period smoothing
        do i=1,ncellsD
          if(wavediss2(i)<wavedisstol)then
            wavediss2(i) = 0.0
            !waveibr2(i) = 0.0
            !wavestrx2(i) = 0.0
            !wavestry2(i) = 0.0
          !else
          !  waveibr2(i) = 1.0
          endif
          !wavestrx2(i) = wavestrx2(i)*waveibr2(i) !Limit radiation stresses to breaking zone
          !wavestry2(i) = wavestry2(i)*waveibr2(i) !Limit radiation stresses to breaking zone
#ifdef DIAG_MODE
          if(isnankind(wavestrx2(i)) .or. isnankind(wavestry2(i)))then
            continue
          endif
#endif
        enddo
      end select 

    endif
    
    return
    endsubroutine wave_init
    
!***************************************************************************
    subroutine wave_print()
!***************************************************************************
    use wave_flowgrid_def
    use diag_def
    implicit none
    integer :: i,iunit(2)
    
888 format(' ',A,T40,A)
222 format(' ',A,T40,F0.2,A)
    
    iunit = (/6, dgunit/)    
    open(dgunit,file=dgfile,access='append')     
    do i=1,2
      if(constant_waves)then
        write(iunit(i),*)        
        write(iunit(i),888) 'Constant Waves:'  
        write(iunit(i),222) '  Significant wave height:',waveheight,' m'
        write(iunit(i),222) '  Peak wave period:',waveperiod,' s'
        write(iunit(i),222) '  Mean wave direction:',wavedir,' deg'
      endif
    enddo
    close(dgunit)
    
    return
    endsubroutine wave_print
    
!!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!!  subroutine to obtain wave breaking energy dissipation term
!!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!    subroutine WVBRK2_alex(CAB,JBV,II,JM,JJ)
!!---------------------------------------------------------------
!!  energy dissipation term is Battjes and Janssen's(1978)
!!    bore model with a breaker parameter based on Miche's criterion of
!!    gama=0.73
!!  parameters used in this subroutine are calculated 
!!    in subroutine veloc
!!---------------------------------------------------------------
!    use global, only: NPF,MPD,IPMX,JPMX,IGPX,JGPX
!    use wave_lib, only: wavebreak_bj78,wave_Hmax
!    implicit none
!
!    integer :: ii,jm,jj,jbv
!    real :: cab,Hrms,wnk,Hmax,fp,Qb
!    integer :: kdate,idate
!    real :: DX,DY,DXX,DMESH,DTH,depmax,depmax0 !,aslope !Removed variable (Alex Sanchez 12/02/14)
!    common /DATD/DX,DY,DXX,DMESH,DTH,kdate,idate,depmax,depmax0 !,aslope(ipmx) !Removed variable (Alex Sanchez 12/02/14)
!    integer :: IBK
!    real :: DBAR,WL0,WCC,HSB,HSG,DCD
!    common /DATG/IBK,DBAR,WL0,WCC(JGPX),HSB(JGPX),HSG(JGPX),DCD(JGPX)
!    integer :: IWVBK,ibk3
!    real :: DEPM,DMNJ,SLF,wlmn,cmn,sigm
!    common /BREK/DEPM(JGPX),DMNJ(JGPX),SLF(JGPX),wlmn(JGPX),cmn(jgpx),  &
!                 sigm(jgpx),IWVBK,ibk3
!    integer :: imod,iprp,island,imd,iprpp,nonln,igrav,isolv,ixmdf,iproc,imud,iwnd
!    real :: PAI2,PAI,HPAI,RAD,akap,depmin0
!    common /VPAI/PAI2,PAI,HPAI,RAD,akap,imod,iprp,island,imd,iprpp,     &
!                 nonln,igrav,isolv,ixmdf,iproc,imud,iwnd,depmin0
!    COMMON /DATC/TP,PRD,IBND,VL,TIDE,KPEAK,IBACK,NBLOCK,IWIND
!    integer :: IBR,KPEAK,IBACK,NBLOCK,IWIND
!    real :: H13S,diss,TP,PRD,IBND,VL,TIDE
!    common /WAVS/H13S(IGPX,JGPX),IBR(IGPX,JGPX),DISS(IGPX,JGPX)             
!    integer :: ICK3,ICK4,ICHOICE,ibrij
!    real :: HS0,wd,ws,ph0,wdd,dissij
!    common /DATB/ICK3,ICK4,ICHOICE,HS0,wd,ws,ph0,wdd(mpd)
!
!    !Initialize
!    cab=0.0 
!    if(IBACK==0)then
!      diss(ii,jj)=0.0
!      ibr(ii,jj)=0    
!    endif 
!    
!! -- ibk=1 and jbv>=6 means no wave breaking
!    !if(IBK==1 .or. JBV>=6) return
!        
!    wnk=pai2/wlmn(jm)
!    if(wnk<1.0e-6) return        
!    
!    !Maximum waveheight
!    call wave_Hmax(depm(jm),slf(jm),wnk,Hmax) !slj is the bed slope
!!!    Hmax = 0.78*depm(jm)    
!    
!    !Check incipient breaking
!    Hrms=HSB(JJ)/1.414213562    
!!    if(Hsb(jj)<=0.72*Hmax .and. ibr(max(ii-1,1),jj)==0) return  !Alex, for random waves
!!    if(Hrms<=Hmax .and. ibr(max(ii-1,1),jj)==0) return  !Alex, only for regular waves
!!    i=max(ii-1,1)
!!    j=max(jj-1,1)
!!    j2=min(jj+1,JGPX)
!!    if(Hrms<=0.95*Hmax .and. (ibr(i,jj)+ibr(i,j)+ibr(i,j2))==0) return  !Alex, only for regular waves
!!    if(Hrms<=Hmax) return  !Alex, only for regular waves
!!    if((ibr(i,jj)+ibr(i,j)+ibr(i,j2))==0) return  !Alex, only for regular waves
!    
!    !Calculate terms
!    fp=sigm(jm)/pai2 !Frequency, 1/s
!    call wavebreak_bj78(Hmax,Hrms,fp,Qb,dissij,ibrij,cab)
!    
!    if(IBACK==0)then
!      diss(ii,jj) = dissij
!      ibr(ii,jj) = ibrij
!    endif
!    
!    return
!    endsubroutine WVBRK2_alex
!
!!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!!  subroutine to obtain wave breaking energy dissipation term
!!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!    subroutine WVBRK4_alex(CAB,JBV,II,JM,JJ)
!!---------------------------------------------------------------
!!  energy dissipation term is Battjes and Janssen's(1978)
!!    bore model with a breaker parameter based on Miche's criterion of
!!    gama=0.73
!!  parameters used in this subroutine are calculated 
!!    in subroutine veloc
!!---------------------------------------------------------------
!    use GLOBAL, only: NPF,MPD,IPMX,JPMX,IGPX,JGPX
!    use wave_lib, only: wavebreak_jb07,wave_Hmax
!    implicit none
!    integer :: ii,jm,jj,jbv
!    real :: cab,Hrms,wnk,Hmax,fp,Qb
!    integer :: kdate,idate
!    real :: DX,DY,DXX,DMESH,DTH,depmax,depmax0,aslop !,aslope !Removed variable (Alex Sanchez 12/02/14)
!    common /DATD/DX,DY,DXX,DMESH,DTH,kdate,idate,depmax,depmax0 !,aslop(ipmx) !Removed variable (Alex Sanchez 12/02/14)
!    integer :: IBK
!    real :: DBAR,WL0,WCC,HSB,HSG,DCD
!    common /DATG/IBK,DBAR,WL0,WCC(JGPX),HSB(JGPX),HSG(JGPX),DCD(JGPX)
!    integer :: IWVBK,ibk3
!    real :: DEPM,DMNJ,SLF,wlmn,cmn,sigm
!    common /BREK/DEPM(JGPX),DMNJ(JGPX),SLF(JGPX),wlmn(JGPX),cmn(jgpx),  &
!                 sigm(jgpx),IWVBK,ibk3
!    integer :: imod,iprp,island,imd,iprpp,nonln,igrav,isolv,ixmdf,iproc,imud,iwnd
!    real :: PAI2,PAI,HPAI,RAD,akap,depmin0
!    COMMON /VPAI/PAI2,PAI,HPAI,RAD,akap,imod,iprp,island,imd,iprpp,     &
!                 nonln,igrav,isolv,ixmdf,iproc,imud,iwnd,depmin0
!    COMMON /DATC/TP,PRD,IBND,VL,TIDE,KPEAK,IBACK,NBLOCK,IWIND
!    integer :: IBR,IBND,KPEAK,IBACK,NBLOCK,IWIND
!    real :: H13S,diss,TP,PRD,VL,TIDE
!    common /WAVS/H13S(IGPX,JGPX),IBR(IGPX,JGPX),DISS(IGPX,JGPX)             
!    integer :: ICK3,ICK4,ICHOICE,ibrij
!    real :: HS0,wd,ws,ph0,wdd,dissij
!    common /DATB/ICK3,ICK4,ICHOICE,HS0,wd,ws,ph0,wdd(mpd)
!    
!   !Initialize
!    cab=0.0 
!    if(IBACK==0)then
!      diss(ii,jj)=0.0
!      ibr(ii,jj)=0    
!    endif
!    
!! -- ibk=1 and jbv>=6 means no wave breaking
!    if(IBK==1 .or. JBV>=6) return
!        
!    wnk=pai2/wlmn(jm)
!    if(wnk<1.0e-6) return    
!            
!    !Maximum waveheight
!    call wave_Hmax(depm(jm),slf(jm),wnk,Hmax) !slj is the bed slope
!!!    Hmax = 0.78*depm(jm)    
!    
!    !Check incipient breaking    
!    Hrms=HSB(JJ)/1.414213562    
!!    if(Hsb(jj)<=0.7*Hmax .and. ibr(max(ii-1,1),jj)==0) return  !Alex, for random waves
!!    i=max(ii-1,1)
!!    j=max(jj-1,1)
!!    j2=min(jj+1,JGPX)
!!    ijbr=ibr(i,jj)+ibr(i,j)+ibr(i,j2)
!!    if(Hsb(jj)<=0.72*Hmax .and. ijbr==0) return  !Alex, for random waves
!!    if(Hsb(jj)<=Hmax .and. ijbr==0) return  !Alex, for regular waves
!!    if(Hsb(jj)<=0.7*Hmax .and. ijbr==0 .and. Hsb(jj)>=0.9*Hs0) return  !Alex, for random waves   
!!    if(Hrms<=Hmax .and. ibr(max(ii-1,1),jj)==0) return  !Alex, only for regular waves
!!    if(Hrms<=Hmax .and. ijbr==0) return  !Alex, only for regular waves
!!    if(Hrms<=0.95*Hmax .and. (ibr(i,jj)+ibr(i,j)+ibr(i,j2))==0) return  !Alex, only for regular waves
!!    if(Hrms<=Hmax) return  !Alex, only for regular waves
!!    if((ibr(i,jj)+ibr(i,j)+ibr(i,j2))==0) return  !Alex, only for regular waves
!    
!    !Calculate terms
!    fp=sigm(jm)/pai2 !Frequency, 1/s
!    call wavebreak_jb07(depm(jm),Hmax,Hrms,fp,Qb,dissij,ibrij,cab) 
!    
!    if(IBACK==0)then
!      diss(ii,jj) = dissij
!      ibr(ii,jj) = ibrij
!    endif
!!    if(jj==10)then
!!      write(*,*) II,JJ,depm(jm),slf(jm),Hmax,Hrms,fp
!!    endif  
!    
!    return
!    endsubroutine WVBRK4_alex
!    
!!****************************************************************    
!    function wave_diss(h,Hs,T,L) result(Dbr)
!!****************************************************************    
!    use const_def, only: twopi,pi,sqrttwo
!    use flow_def, only: rhow,grav
!    use prec_def
!    implicit none
!    real(ikind),intent(in) :: h,Hs,T,L
!    real(ikind) :: Hb,Hrms,kh,Hb2,Hrms2,Qb,val,Dbr
!    
!    kh=twopi*h/L
!    Hb=0.14*L*tanh(0.78*kh/0.88)
!    Hrms=Hs/sqrttwo
!    Hb2=Hb*Hb
!    Hrms2=Hrms*Hrms
!!    !Baldock
!!    Qb=exp(-Hb2/Hrms2)
!!    Dbr=0.25*Qb*rhow*grav/T*(Hb2+Hrms2)
!    !Battjes and Janssen
!    val=Hrms2/Hb2    
!    Qb=exp((0.9698*val-1.0)/(0.6574*val))    
!    if(Qb<0.5)then
!      Dbr=0.0
!      return
!    endif
!    Dbr=0.25*Qb*rhow*grav/T*Hb2  !N/m^2/s
!    
!    return
!    endfunction wave_diss
!
!!********************************************************************************
!    subroutine stokes_stress
!! Adds the Stokes stress term to the wave radiation stresses
!!
!! written by Alex Sanchez, USACE-CHL, May 2012
!!*********************************************************************************    
!    use flow_wavegrid_def, only: uwave,vwave,hwave
!    use wave_wavegrid_def, only: nwavei,nwavej,xwave,ywave,dxwav,dywav,wheight,wperiod,wdiss,wcos,wsin,wxrs1,wyrs1,wibr
!    use wave_flowgrid_def, only: wavestrx2,wavestry2
!    use rol_def, only: d2x,d2y
!    use flow_def, only: waveflux,rhow,hmin,grav
!    use const_def, only: twopi,pi,small
!    use prec_def
!    implicit none 
!
!    integer :: i,i1,i2,j,j1,j2
!    real(ikind) :: Sxx(nwavei,nwavej),Sxy(nwavei,nwavej),Syy(nwavei,nwavej)
!    real(ikind) :: E,wlength,c,hw !,hk,wavelength,val
!    real*4 :: q,d,uu,vv,om,cw,cg,sg,akk !,f,uvw !must be single 
!    
!    !Stokes velocities and stresses
!    do i=1,nwavei      
!      do j=1,nwavej        
!        hw=hwave(i,j)
!        if(hw<0.01 .or. wheight(i,j)<1.0e-4)then
!          Sxx(i,j)=0.0; Syy(i,j)=0.0; Sxy(i,j)=0.0
!          cycle
!        endif
!!        wlen=wavelength(wperiod(i,j),hw) !No wave-current interaction
!!        f = 1.0/wperiod(i,j)
!!        uvw = uwave(i,j)*wcos(i,j)+vwave(i,j)*wsin(i,j)
!!        call wkcgen(f,grav,pi,uvw,hw,iflag,akk) !STWAVE subroutine
!        q=atan2(wsin(i,j),wcos(i,j)); d=hw; uu=uwave(i,j); vv=vwave(i,j); om=twopi/wperiod(i,j)
!        call wccg(q,d,uu,vv,om,cw,cg,sg,akk) !CMS-Wave subroutine
!        !call wccg3(q,d,uu,vv,om,cw,cg,sg,akk) !Alex's subroutine
!        wlength=twopi/akk
!        c=wlength/wperiod(i,j)
!        E=0.0625*rhow*grav*wheight(i,j)*wheight(i,j)        
!        uu=E*wcos(i,j)/(rhow*hw*c)
!        vv=E*wsin(i,j)/(rhow*hw*c)
!        Sxx(i,j)=hw*uu*uu
!        Sxy(i,j)=hw*uu*vv 
!        Syy(i,j)=hw*vv*vv        
!      enddo
!    enddo
!    
!    !Add stokes stress gradients to wave forcing
!   ! do i=1,nwavei
!      !i1=min(i+1,nwavei); i2=max(i-1,1)
!      !do j=1,nwavej
!      !  j1=min(j+1,nwavej); j2=max(j-1,1)        
!   !     wxrs1(i,j)=wxrs1(i,j)-(Sxx(i1,j)-Sxx(i2,j))/d2x(i)
!   !     wxrs1(i,j)=wxrs1(i,j)-(Sxy(i,j1)-Sxy(i,j2))/d2y(j)
!   !     wyrs1(i,j)=wyrs1(i,j)-(Syy(i,j1)-Syy(i,j2))/d2y(j)
!   !     wyrs1(i,j)=wyrs1(i,j)-(Sxy(i1,j)-Sxy(i2,j))/d2x(i)
!   !   enddo
!   ! enddo
!   
!   ! do i=1,nwavei,nwavi-1
!      !i1=min(i+1,nwavei); i2=max(i-1,1)
!      !do j=1,nwavej,nwavej-1
!      !  j1=min(j+1,nwavej); j2=max(j-1,1)        
!   !     wxrs1(i,j)=wxrs1(i,j)-(Sxx(i1,j)-Sxx(i2,j))/delx(i)
!   !     wxrs1(i,j)=wxrs1(i,j)-(Sxy(i,j1)-Sxy(i,j2))/dely(j)
!   !     wyrs1(i,j)=wyrs1(i,j)-(Syy(i,j1)-Syy(i,j2))/dely(j)
!   !     wyrs1(i,j)=wyrs1(i,j)-(Sxy(i1,j)-Sxy(i2,j))/delx(i)
!   !   enddo
!   ! enddo
!   
!    do i=2,nwavei-1
!      i1=min(i+1,nwavei); i2=max(i-1,1)
!      i1=i+1; i2=i-1
!      do j=2,nwavej-1
!        j1=j+1; j2=j-1       
!        wxrs1(i,j)=wxrs1(i,j)-(Sxx(i1,j)-Sxx(i2,j))/d2x(i)
!        wxrs1(i,j)=wxrs1(i,j)-(Sxy(i,j1)-Sxy(i,j2))/d2y(j)
!        wyrs1(i,j)=wyrs1(i,j)-(Syy(i,j1)-Syy(i,j2))/d2y(j)
!        wyrs1(i,j)=wyrs1(i,j)-(Sxy(i1,j)-Sxy(i2,j))/d2x(i)
!      enddo
!    enddo   
!    
!    return 
!    endsubroutine stokes_stress
!    
!!********************************************************************************
!    subroutine wave_rad_stress
!! Over-rides the wave radiation stresses and is calculated based 
!! on the wave energy by uncommenting lines below.
!!
!! written by Alex Sanchez, USACE-CHL, May 2012
!!*********************************************************************************    
!    use flow_wavegrid_def, only: uwave,vwave,hwave
!    use wave_wavegrid_def, only: nwavei,nwavej,xwave,ywave,dxwav,dywav,&
!         wheight,wperiod,wdiss,wcos,wsin,wxrs1,wyrs1,wibr
!    use wave_flowgrid_def, only: wavestrx2,wavestry2
!    use rol_def, only: d2x,d2y
!    use flow_def, only: waveflux,rhow,hmin,grav
!    use const_def, only: twopi,pi,small
!    use prec_def
!    
!    implicit none 
!    integer :: i,i1,i2,j,j1,j2
!    real(ikind) :: Sxx(nwavei,nwavej),Sxy(nwavei,nwavej),Syy(nwavei,nwavej)
!    real(ikind) :: E,wlength,c,cn,hwkw,hw
!    real*4 :: q,d,uu,vv,om,cw,cg,sg,akk !must be single 
!    
!!$OMP PARALLEL DO PRIVATE(i,j,q,d,uu,vv,om,cw,cg,sg,akk,E,wlength,c,cn,hwkw,hw)
!    do i=1,nwavei      
!      do j=1,nwavej        
!        hw=hwave(i,j)
!        if(hw<0.01)then
!          Sxx(i,j)=0.0; Syy(i,j)=0.0; Sxy(i,j)=0.0
!          cycle
!        endif
!!        wlen=wavelength(wperiod(i,j),hw) !No wave-current interaction
!!        f = 1.0/wperiod(i,j)
!!        uvw = uwave(i,j)*wcos(i,j)+vwave(i,j)*wsin(i,j)
!!        call wkcgen(f,grav,pi,uvw,hw,iflag,akk) !STWAVE subroutine
!        q=atan2(wsin(i,j),wcos(i,j)); d=hw; uu=uwave(i,j); vv=vwave(i,j); om=twopi/wperiod(i,j)
!        call wccg(q,d,uu,vv,om,cw,cg,sg,akk) !CMS-Wave subroutine
!        !call wccg3(q,d,uu,vv,om,cw,cg,sg,akk) !Alex's subroutine
!        wlength=twopi/akk
!        c=wlength/wperiod(i,j)
!        hwkw=twopi/wlength*hw
!        hwkw=min(hwkw,10.0)
!        cn=0.5+hwkw/sinh(2.0*hwkw)
!        E=0.0625*rhow*grav*wheight(i,j)*wheight(i,j)
!        !Radiation stresses        
!        Sxx(i,j)=E*(cn*(1.0+wcos(i,j)*wcos(i,j))-0.5)
!        Syy(i,j)=E*(cn*(1.0+wsin(i,j)*wsin(i,j))-0.5)
!        Sxy(i,j)=E*cn*wcos(i,j)*wsin(i,j)        
!      enddo
!    enddo
!!$OMP END PARALLEL DO
!
!    !Over-ride wave radiation stress gradients
!!$OMP PARALLEL DO PRIVATE(i,i1,i2,j,j1,j2)
!    do i=1,nwavei
!      i1=min(i+1,nwavei); i2=max(i-1,1)
!      do j=1,nwavej
!        j1=min(j+1,nwavej); j2=max(j-1,1)        
!        wxrs1(i,j)=(-(Sxx(i1,j)-Sxx(i2,j))/d2x(i)-(Sxy(i,j1)-Sxy(i,j2))/d2y(j))/rhow
!        wyrs1(i,j)=(-(Syy(i,j1)-Syy(i,j2))/d2y(j)-(Sxy(i1,j)-Sxy(i2,j))/d2x(i))/rhow
!      enddo
!    enddo   
!!$OMP END PARALLEL DO
!    
!    return 
!    endsubroutine wave_rad_stress
!
!********************************************************************************
    subroutine wave_dissipation
! Adds the Stokes stress term to the wave radiation stresses
! and calculates the wave dissipation based on the divergence of the 
! wave energy flux.
!
! Note: The spectral wave radiation stresses from CMS-Wave can be overwritten 
!       and calculated here based on the wave energy by uncommenting lines below.
!
! written by Alex Sanchez, USACE-CHL, May 2012
!*********************************************************************************    
    use flow_wavegrid_def, only: uwave,vwave,hwave
    use wave_wavegrid_def, only: nwavei,nwavej,xwave,ywave,dxwav,dywav,&
          !wheight,wperiod,wdiss,wcos,wsin,wxrs1,wyrs1,wibr
          wheight,wperiod,wdiss,wcos,wsin,wxrs1,wyrs1
    use wave_flowgrid_def, only: wavestrx2,wavestry2
    use rol_def, only: d2x,d2y
    use flow_def, only: waveflux,rhow,hmin,grav
    use const_def, only: twopi,pi,small
    use prec_def
    implicit none 
    
    integer :: i,i1,i2,j,j1,j2
    real(ikind) :: ECg(nwavei,nwavej)
    real(ikind) :: E,wlength,c,cn,hk,hw,val
    real*4 :: q,d,uu,vv,om,cw,cg,sg,akk !must be single 
    
!$OMP PARALLEL DO PRIVATE(i,j,q,d,uu,vv,om,cw,cg,sg,akk,E,wlength,c,cn,hk,hw)
    do i=1,nwavei      
      do j=1,nwavej        
        hw=hwave(i,j)
        if(hw<0.01)then
          ECg(i,j)=0.0
          cycle
        endif
!        wlen=wavelength(wperiod(i,j),hw) !No wave-current interaction
!        f = 1.0/wperiod(i,j)
!        uvw = uwave(i,j)*wcos(i,j)+vwave(i,j)*wsin(i,j)
!        call wkcgen(f,grav,pi,uvw,hw,iflag,akk) !STWAVE subroutine
        q=atan2(wsin(i,j),wcos(i,j)); d=hw; uu=uwave(i,j); vv=vwave(i,j); om=twopi/wperiod(i,j)
        !call wccg(q,d,uu,vv,om,cw,cg,sg,akk) !CMS-Wave subroutine
        call wccg_inline(q,d,uu,vv,om,cw,cg,sg,akk) !CMS-Wave subroutine   !Lihwa  12/18/2020
        !call wccg3(q,d,uu,vv,om,cw,cg,sg,akk) !Alex's subroutine
        wlength=twopi/akk
        c=wlength/wperiod(i,j)
        hk=twopi/wlength*hw
        cn=0.5+hk/sinh(2.0*hk)
        E=0.0625*rhow*grav*wheight(i,j)*wheight(i,j)
        ECg(i,j)=E*cn*c    
      enddo
    enddo
!$OMP END PARALLEL DO

    !Modifies the wave dissipation
!$OMP PARALLEL DO PRIVATE(i,i1,i2,j,j1,j2)
    do i=1,nwavei
      i1=min(i+1,nwavei); i2=max(i-1,1)
      do j=1,nwavej
        j1=min(j+1,nwavej); j2=max(j-1,1)        
        val=-(ECg(i1,j)-ECg(i2,j))/d2x(i)-(ECg(i,j1)-ECg(i,j2))/d2y(j) !Central difference
        !wdiss(i,j)=max(wdiss(i,j),abs(val))*wibr(i,j)										 !
        wdiss(i,j)=max(wdiss(i,j),abs(val))	  !Lihwa 12/18/2020
        !wdiss(i,j)=abs(val)*wibr(i,j)
        !wdiss(i,j)=abs(val) !*wibr(i,j)
      enddo
    enddo   
!$OMP END PARALLEL DO

    return 
    endsubroutine wave_dissipation

!******************************************************
    subroutine wavebreaking(whgt,wper,wdir,h,u,v,Qb,cbr)
!******************************************************
    use size_def, only: ncells,ncellsD
    use geo_def, only: dzbx,dzby
    use comvarbl, only: ramp
    use const_def, only: sqrttwo,twopi
    use wave_lib, only: wavenumber,wavebreak_bj78,wave_Hmax
    use prec_def
    implicit none
    
    integer :: i,ib
    real(ikind),intent(in) :: h(ncellsD),u(ncellsD),v(ncellsD)
    real(ikind),intent(in) :: whgt(ncellsD),wper(ncellsD),wdir(ncellsD)
    real(ikind),intent(out):: Qb(ncellsD),cbr(ncellsD)
    real(ikind) :: wa
    real(4) :: h4,s,wk,Hmax,Hrms,val,fp,Db,Qb4,cbr4  !must be single

    s = 0.0_ikind
    val = sqrttwo*sqrt(max(ramp,1.0e-5))
    do i=1,ncells
      wa = twopi/wper(i)
      wk = wavenumber(wa,wdir(i),h(i),u(i),v(i))
      s = dzbx(i)*cos(wdir(i))+dzby(i)*sin(wdir(i))  !Bed slope in wave direction
      h4 = h(i)
      call wave_Hmax(h4,s,wk,Hmax)
      Hrms = whgt(i)/val
      fp = 1.0/wper(i)
      !call wavebreak_jb07(dep(i),Hmax,Hrms,fp,Qb(i),Db,ib,cbr(i)) !Jansen and Battjes (2007)
      call wavebreak_bj78(Hmax,Hrms,fp,Qb4,Db,ib,cbr4)         !Battjes and Jansen (1978)
      Qb(i) = Qb4
      cbr(i) = cbr4
    enddo

    return
    endsubroutine wavebreaking
    
!********************************************************
    subroutine wccg3(wd,h,u,v,wa,cw,cg,sig,wk)
! Description:
!  Solves the wave dispersion relation 
!   sig^2 = grav*wk*tanh(wk*h)
!  where 
!   sig = wa - k*cos(wd)*u - k*cos(wd)*v = wa - k*uk
!   sig = relative angular frequency [rad/s]
!   wk = wave number [rad/m]
!   wa = absolute angular frequency [rad/s]
!   wd = wave direction (Cartesian, going to) [rad]
!   u = current velocity in x direction [m/s]
!   v = current velocity in y direction [m/s]
!   uk = cos(q)*u + sin(q)*v [m/s]
!
! Uses the Newton-Raphson Method:
!   wk(n+1) = wk(n) - f(wk(n))/fp(wk(n))
! where
!   f = grav*wk*tanh(wk*h) - sig^2
!   fp = grav*(tanh(h*k)-h*k*(tanh(h*k)^2-1)) + 2*uk*sig
!
! Makes initial guess using 
!   Guo, J. 2002. Simple and explicit solution of wave
!      dispersion. Coastal Engineering, 46(2), 71-74.

!  when no solution and cg+u<0, sg=-1,....., are set
!********************************************************
    use diag_def
    use diag_lib
    implicit none    

    !Input/Output
    real,intent(in) :: wd,h,u,v,wa
    real,intent(out) :: cw,cg,sig,wk
    !Internal
    integer :: i
    real :: wki,uk,f,fp,tanhkh,xi,yi,wkh,err
    real, parameter :: g = 9.806
    real, parameter :: tol = 1.0e-4
    logical :: isnansingle
      
    !Initial wave number guess based on Guo (2002)
    xi = wa/sqrt(g/h) !=h*wa/sqrt(g*h) 
    yi = xi*xi/(1.0-exp(-xi**2.4908))**0.4015
    wki = yi/h !Initial guess without current-wave interaction 

    !Current velocity in wave direction    
    uk = cos(wd)*u + sin(wd)*v 
    
    !Simple correction for currents (tends to overcorrect, hense the relaxation)
    wki = 0.3*wki + 0.7*(wa-wki*uk)**2/g/tanh(wki*h)      

    !Newton-Raphson Iterations
    wk = wki   !save for iteration
    do i=1,10
      sig = wa-wk*uk
      tanhkh = tanh(wk*h)  
      f = g*wk*tanhkh - sig**2
      fp = g*(tanhkh - h*wk*(tanhkh**2-1.0)) + 2*uk*sig
      wk = wki - f/fp
      if(wk<1.0e-6 .or. wk>1.0e6)then !Did not converge, assume wave blocking                        
        !wk = 0.0; cw = 0.0; cg = 0.0; sig = -1.0 !Set to zero
        !wk = yi/h !Use initial guess
        wk = yi/h*2.5 !Use a small wave length (large wave number) to avoid divide by zeros        
        write(msg2,*) 'wa=',wa
        write(msg3,*) 'wd=',wd
        write(msg4,*) 'h=',h
        write(msg5,*) 'u=',u
        write(msg6,*) 'v=',v
        call diag_print_warning('Wave blocking detected in wavenumber',msg2,msg3,msg4,msg5,msg6)
        exit  
      elseif(isnansingle(wk))then
        wk = yi/h*2.5 !Use a small wave length (large wave number) to avoid divide by zeros
        write(msg2,*) 'wa=',wa
        write(msg3,*) 'wd=',wd
        write(msg4,*) 'h=',h
        write(msg5,*) 'u=',u
        write(msg6,*) 'v=',v
        call diag_print_error('Calculated NaN value for wk in wavenumber',msg2,msg3,msg4,msg5,msg6)
      endif    
      err = abs(wk-wki)
      if(err < tol) exit
      wki = wk
    enddo
    
    sig = wa-wk*uk
    cw=sig/wk !Wave celerity [m/s]
    wkh=min(wk*h,10.0) 
    cg=cw*(0.5+wk*h/sinh(2.0*wkh))  !Group velocity [m/s]

    return
    endsubroutine wccg3

!***************************************************************
    subroutine brkwavratio(dep,Hrms,L,gamma,Qb)
! Calculates the fraction of broken waves 
! according to Baldock et al. (1998)
!***************************************************************
    use const_def, only: twopi
    use prec_def
    implicit none
    real(ikind) :: Hrms,dep,Qb,Hb,wk,L,gamma
    
    wk = twopi/L !Wave number    
    Hb = 0.88/wk*tanh(gamma*wk*dep/0.88) !Miche (1951) criterion, maximum waveheight
    Qb = exp(-(Hb/Hrms)**2) !proportion of broken waves      
    
    return
    endsubroutine brkwavratio    
    
!***********************************************************************    
    subroutine wave_const_init
!***********************************************************************        
    use size_def
    use sed_def  
    use wave_flowgrid_def
    use geo_def, only: azimuth_fl
    use cms_def
    use const_def, only: deg2rad
    use prec_def
    implicit none
    integer :: i 
    
    !Constant wave parameters
    noptset = 5
    allocate(Whgt(ncellsD),Wper(ncellsD),wavediss(ncellsD))
    allocate(Wang(ncellsD),Wunitx(ncellsD),Wunity(ncellsD))    
    allocate(wavestrx(ncellsD),wavestry(ncellsD))
    allocate(Ssr(ncellsD))
    do i=1,ncellsD
      Whgt(i) = waveheight
      Wper(i) = waveperiod      
      Wang(i) = wavedir
      wavediss(i) = 0.0
      wavestrx(i) = 0.0
      wavestry(i) = 0.0
      Ssr(i) = 0.0
      Wunitx(i) = cos((Wang(i)-azimuth_fl)*deg2rad)
      Wunity(i) = sin((Wang(i)-azimuth_fl)*deg2rad)      
    enddo
    
    allocate(Worb(ncellsD),Worbrep(ncellsD),Wlen(ncellsD))
    call wave_const_update
    
    return
    endsubroutine wave_const_init
    
!***********************************************************************   
    subroutine wave_const_update()
! Prints sediment transport variables for debugging
! written by Alex Sanchez, USACE-CHL
!***********************************************************************
    use size_def
    use flow_def
    use sed_def  
    use wave_flowgrid_def
    use const_def, only: twopi,deg2rad
    use wave_lib
    use cms_def
    use prec_def
    implicit none
    integer :: i
    real(ikind) :: wa,wd,tol
     
    tol = 0.001_ikind
    do i=1,ncells
      Worb(i) = 0.0_ikind
      Worbrep(i) = 0.0_ikind
      Wlen(i) = 10.0_ikind
      if(waveheight>0.0)then
        wa = twopi/waveperiod
        wd = wavedir*deg2rad
        Wlen(i) = wavelength(wa,wd,h(i),u(i),v(i),tol)
        !Orbital velocities (includes ramp and wet/dry through Whgt)
         Worb(i) = waveorb_linear(h(i),Whgt(i),Wper(i),Wlen(i))
!!       Worbrep(i) = Worb(i)/sqrttwo !Old approach
         Worbrep(i) = waveorbrep_ss87(grav,h(i),Whgt(i),Wper(i)) !More accurate  
      endif    
    enddo
        
    return
    endsubroutine wave_const_update