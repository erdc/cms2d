!***************************************************************************   
    subroutine stat_default
! Output Defaults
! written by Alex Sanchez
!***************************************************************************        
    use stat_def
    use comvarbl, only: casename,flowpath
    implicit none
    
    calc_stats = .false.
    tstat(:) = -1.0
    nstatcount = 0   
    
    flowstats = .false.   
    tflowstat(1:3) = -1.0    
    nflowstatcount = 0 
        
    sedstats = .false. 
    tsedstat(:) = -1.0
    nsedstatcount = 0     
    
    salstats = .false.
    tsalstat(:) = -1.0
    nsalstatcount = 0
    
    wavestats = .false.
    twavstat(:) = -1.0
    nwavestatcount = 0  
    
    statfile = trim(flowpath) // trim(casename) // '_stat.h5'    
    
    return
    endsubroutine stat_default
    
!***************************************************************************   
    subroutine stat_cards(cardname,foundcard)
! Sets the default parameters for flow variables
! written by Alex Sanchez, USACE-CHL; Weiming Wu, NCCHE
!***************************************************************************        
    use stat_def
    use unitconv_lib, only: unitconv_vec
    implicit none 
    integer :: ierr
    character(len=37) :: cardname
    character(len=200) :: aline
    character(len=10) :: fromunits,tounits
    logical :: foundcard
    
    fromunits = 'hrs'
    tounits = 'sec'
    foundcard = .true.
    selectcase(cardname)  
      case('GLOBAL_STATISTICS')
        backspace(77)
        read(77,'(A)') aline
        read(aline,*,iostat=ierr) cardname,tstat(1),tstat(2),tstat(3),fromunits !last argument is optional internal
        calc_stats = .true.
        call unitconv_vec(fromunits,tounits,3,tstat)
        
      case('HYDRO_STATISTICS','FLOW_STATISTICS')
        backspace(77)
        read(77,'(A)') aline
        read(aline,*,iostat=ierr) cardname,tflowstat(1),tflowstat(2),tflowstat(3),fromunits !last argument is optional internal
        flowstats = .true.
        call unitconv_vec(fromunits,tounits,3,tflowstat)
        
      case('SEDIMENT_STATISTICS','SED_STATISTICS')
        backspace(77)
        read(77,'(A)') aline
        read(aline,*,iostat=ierr) cardname,tsedstat(1),tsedstat(2),tsedstat(3),fromunits !last argument is optional internal
        sedstats = .true. 
        call unitconv_vec(fromunits,tounits,3,tsedstat)
        
      case('SALINITY_STATISTICS','SALT_STATISTICS')
        backspace(77)
        read(77,'(A)') aline
        read(aline,*,iostat=ierr) cardname,tsalstat(1),tsalstat(2),tsalstat(3),fromunits !last argument is optional internal
        salstats = .true.
        call unitconv_vec(fromunits,tounits,3,tsalstat)
        
      case('WAVE_STATISTICS')
        backspace(77)
        read(77,'(A)') aline
        read(aline,*,iostat=ierr) cardname,twavstat(1),twavstat(2),twavstat(3),fromunits !last argument is optional internal
        wavestats = .true. 
        call unitconv_vec(fromunits,tounits,3,twavstat)
        
      case default
        foundcard = .false.
        
    endselect
    
    return
    endsubroutine stat_cards
    
!***********************************************************    
    subroutine stat_init
! Initialize global statistics variables
! There are several variables that contain the timing information for the calculations of statistics
! depending on the type of statistics.  These are: 
!   tstat     (global stats)
!   tflowstat (hydrodynamic stats)
!   tsedstat  (sediment stats)
!   tsalstat  (salinity stats)
!   twavstat  (wave stats)
! These variables are each composed as an array of three timing values as follows:
!   tstat(1) - Beginning time
!   tstat(2) - Ending time
!   tstat(3) - Update interval
!
! written by Alex Sanchez, USACE CHL
!***********************************************************    
    use size_def        
    use cms_def, only: noptset
    use geo_def, only: dx,dy
    use comvarbl, only: dtimebeg,dtime
    use flow_def, only: h,u,v,grav
    use out_def, only: simlabel
    use sed_def, only: sedtrans
    use sal_def, only: saltrans
    use der_lib, only: d2xy2d
    use stat_def
    use hot_def, only: coldstart
    use const_def, only: small
    implicit none          
    integer :: nn
    !real(ikind) :: d2varx2,d2vary2
    !real(ikind) :: dvarx,dvary
    logical :: ok
    
    nn = len_trim(simlabel)  
    statpath = simlabel(1:nn) // '/'
    if(coldstart)then
      inquire(file=statfile,exist=ok)
      if(ok)then
        open(45,file=statfile)
        close(45,status='delete')     
      endif
    endif
    
    !Global Statistics
    if(calc_stats)then
      flowstats = .true.
      if(noptset>=3) wavestats = .true.
      !Interval
      tstat(3) = max(tstat(3),dtimebeg)    
      !Sediment statistics
      if(sedtrans) sedstats = .true.            
      !Salinity statistics
      if(saltrans) salstats = .true.            
    endif
    
    !Flow Statistics
    if(flowstats) call stat_init_flow
    
    !Sediment Statistics
    if(.not.sedtrans) sedstats = .false.
    if(sedstats) call stat_init_sed
    
    !Salinity Statistics
    if(.not.saltrans) salstats = .false.
    if(salstats) call stat_init_sal

    !Waves
    if(noptset<3) wavestats = .false.
    if(wavestats) call stat_init_wave        
    
    return
    endsubroutine stat_init

!*********************************************************************************
    subroutine stat_init_flow
! Author: Alex Sanchez, USACE-CHL
!*********************************************************************************
    use size_def
    use geo_def
    use comvarbl, only: dtimebeg,dtime
    use flow_def, only: h,u,v,grav
    use out_def, only: simlabel
    use der_lib, only: d2xy2d
    use stat_def
    use hot_def, only: coldstart
    use const_def, only: small
    implicit none          
    integer :: i,ii,k,nck,idmax,idmin
    real(ikind) :: val
    !real(ikind) :: d2varx2,d2vary2
    !real(ikind) :: dvarx,dvary
    
    !Update Interval
    tflowstat(3) = max(tflowstat(3),dtimebeg)
      
    !Replace with global if not declared
    if(tflowstat(2)<0.0) tflowstat = tstat
    
    !Water levels
    allocate(etamax(ncellsD),etamin(ncellsD),etamean(ncellsD))
    etamax = 0.0; etamin = 0.0; etamean = 0.0
    allocate(tetamax(ncellsD),tetamin(ncellsD))
    tetamax = 0.0; tetamin = 0.0
    
    !Current velocities
    allocate(uvmax(ncellsD),umax(ncellsD),vmax(ncellsD))
    uvmax = 0.0; umax = 0.0; vmax = 0.0
    allocate(tuvmax(ncellsD))
    tuvmax = 0.0
    
    !Residual flow velocities
    allocate(uvr(ncellsD),ur(ncellsD),vr(ncellsD))
    uvr = 0.0; ur = 0.0;  vr = 0.0
    
    !Curvature
    allocate(uvcurv(ncellsD))
    uvcurv = 0.0
      
    !Hydro period
    allocate(hydper(ncellsD))
    hydper = 0.0    
      
    !Flow Gradients
    allocate(etamaxgrad(ncellsD),uvmaxgrad(ncellsD))
    etamaxgrad = 0.0; uvmaxgrad = 0.0      
      
    !Average Residuals
    allocate(rspavg(ncellsD),rsuavg(ncellsD),rsvavg(ncellsD))
    rspavg = 0.0; rsuavg = 0.0; rsvavg = 0.0
      
    !Maximum Residuals
    allocate(rspmax(ncellsD),rsumax(ncellsD),rsvmax(ncellsD))
    rspmax = 0.0; rsumax = 0.0; rsvmax = 0.0
      
    !CFL condition
    allocate(cflxmax(ncellsD),cflymax(ncellsD))
    cflxmax = 0.0; cflymax = 0.0
      
    !!Reynolds numbers
    !allocate(rexmax(ncellsD),reymax(ncellsD))
    !rexmax = 0.0; reymax = 0.0
      
    !Froude numbers
    allocate(frxmax(ncellsD),frymax(ncellsD))
    frxmax = 0.0; frymax = 0.0
      
    !Flow curvature
    allocate(umaxcurvx(ncellsD),umaxcurvy(ncellsD))
    umaxcurvx = 0.0; umaxcurvy = 0.0
    allocate(vmaxcurvx(ncellsD),vmaxcurvy(ncellsD))
    vmaxcurvx = 0.0; vmaxcurvy = 0.0
      
    if(ncellpoly==0)then
      !!Grid quality indicator based on long-wave celerity
      !allocate(qcelmax(ncellsD))
      !qcelmax = 0.0
      
      !!Grid quality indicators based on depth slopes
      !allocate(qhgradxmax(ncellsD),qhgradymax(ncellsD)) 
      !qhgradxmax = 0.0; qhgradymax = 0.0
      
      !!Grid quality indicators based on depth curvatures
      !allocate(qhcurvxmax(ncellsD),qhcurvymax(ncellsD))  
      !qhcurvxmax = 0.0; qhcurvymax = 0.0
      
      !Grid quality indicators based on courant numbers
      allocate(qcflxmax(ncellsD),qcflymax(ncellsD))  
      qcflxmax = 0.0; qcflymax = 0.0
      
      !Grid quality indicators based on Froude numbers
      allocate(qFrxmax(ncellsD),qFrymax(ncellsD))  
      qFrxmax = 0.0; qFrymax = 0.0
    endif
    
    !Grid Smoothness
    allocate(qareachg(ncellsD))    
    qareachg=0.0
    do i=1,ncells
      idmax=i; idmin=i
      do k=1,ncface(i)
        nck=cell2cell(k,i)  
        if(areap(nck)>areap(idmax))then
          idmax=nck
        elseif(areap(nck)<areap(idmin))then
          idmin=nck
        endif    
      enddo
      if(idmax/=idmin)then
        qareachg(i)=(areap(idmax)-areap(idmin))/sqrt(areap(idmin))&
          /sqrt((x(idmax)-x(idmin))**2+(y(idmax)-y(idmin))**2)
      endif
      !Too slow
      !do ii=1,ncells
      !  if(ii==i) cycle
      do k=1,ncface(i)        
        ii=cell2cell(k,i)
        !val=abs(sqrt(areap(i))-sqrt(areap(ii)))/sqrt((x(i)-x(ii))**2+(y(i)-y(ii))**2)
        !val=abs(areap(i)-areap(ii))/((x(i)-x(ii))**2+(y(i)-y(ii))**2)
        !val=abs(areap(i)**2-areap(ii)**2)/((x(i)-x(ii))**2+(y(i)-y(ii))**2)**2
        val=abs(areap(i)-areap(ii))/sqrt(min(areap(i),areap(ii)))&
           /sqrt((x(i)-x(ii))**2+(y(i)-y(ii))**2)
        qareachg(i)=max(qareachg(i),val)
      enddo
    enddo
    !allocate(temp(ncellsD))
    !do i=1,ncells    
    !  temp(i)=qareachg(i)
    !  nsum=0
    !  do k=1,ncface(i)
    !    if(qareachg(cell2cell(k,i))>qareachg(i))then
    !      nsum=nsum+1
    !    endif        
    !  enddo
    !  if(nsum>1)then
    !    do k=1,ncface(i)
    !      temp(i)=max(temp(i),qareachg(cell2cell(k,i)))
    !    enddo  
    !  endif
    !enddo        
    !qareachg=temp
    
    return
    endsubroutine stat_init_flow
  
!*********************************************************************
    subroutine stat_init_sed
!*********************************************************************
    use size_def
    use stat_def
    use comvarbl, only: dtimebeg
    implicit none
    
    !Update Interval
    tsedstat(3) = max(tsedstat(3),dtimebeg)
      
    !Replace with global if not declared
    if(tsedstat(2)<0.0) tsedstat = tstat

    !Maximum sediment transport
    allocate(qtxmax(ncellsD), qtymax(ncellsD), qtmax(ncellsD))
    qtxmax = 0.0; qtymax = 0.0; qtmax = 0.0
        
    !Net sediment transport    
    allocate(qtxn(ncellsD), qtyn(ncellsD), qtn(ncellsD))
    qtxn = 0.0; qtyn = 0.0; qtn = 0.0
    
    !Gross sediment transport    
    allocate(qtg(ncellsD))
    qtg = 0.0    
      
    !Residuals
    allocate(rsCtm(ncellsD))
    rsCtm = 0.0
    
    return
    endsubroutine stat_init_sed

!************************************************************
    subroutine stat_init_sal
!************************************************************
    use size_def
    use stat_def
    use comvarbl, only: dtimebeg
    implicit none
    
    !Update Interval
    tsalstat(3) = max(tsalstat(3),dtimebeg)
      
    !Replace with global if not declared
    if(tsalstat(2)<0.0) tsalstat = tstat
      
    !Maximum
    allocate(salmax(ncellsD))
    salmax = 0.0
    
    !Minimum
    allocate(salmin(ncellsD))
    salmin = 1.0e20
    
    !Average
    allocate(Salavg(ncellsD))
    Salavg = 0.0
      
    !Residuals
    allocate(rsSalm(ncellsD))
    rsSalm = 0.0
    
    return
    endsubroutine stat_init_sal
    
!************************************************************************
    subroutine stat_init_wave
!***********************************************************************
    use size_def
    use stat_def
    use comvarbl, only: dtimebeg
    implicit none
    
    !Interval
    if(twavstat(3)<0.0) twavstat(3) = dtimebeg 
    !Replace with global if not declared
    if(twavstat(2)<0.0) twavstat = tstat      
    !Statistics
    allocate(whgtmax(ncellsD),whgtmean(ncellsD),wpermean(ncellsD))
    allocate(whgtmeanx(ncellsD),whgtmeany(ncellsD))
    allocate(ursellmax(ncellsD),wdissmax(ncellsD))
    whgtmax = 0.0; whgtmean = 0.0; wpermean = 0.0
    whgtmeanx = 0.0; whgtmeany = 0.0
    ursellmax = 0.0
    wdissmax = 0.0
      
    return
    endsubroutine stat_init_wave
    
!***********************************************************    
    subroutine stat_print()
! Prints statistics settings
! written by Alex Sanchez, USACE-CHL
!***********************************************************    
    use stat_def
    use diag_def, only: dgfile,dgunit
    implicit none
    integer :: iunit(2),i

    iunit = (/6, dgunit/)    
    open(dgunit,file=dgfile,access='append') 
    do i=1,2
      if(flowstats) call stat_print_group(iunit(i),'Flow'    ,tflowstat)    
      if(sedstats)  call stat_print_group(iunit(i),'Sediment',tsedstat) 
      if(salstats)  call stat_print_group(iunit(i),'Salinity',tsalstat) 
      if(wavestats) call stat_print_group(iunit(i),'Waves'   ,twavstat)    
    enddo    
    close(dgunit)
    
    return
    contains
!---------------------------------------------------------------------
      subroutine stat_print_group(iuniti,groupname,tgroupstat)
!---------------------------------------------------------------------
      use prec_def
      use tool_def
      implicit none
      integer          :: iuniti,ierr
      real(ikind)      :: tgroupstat(3)
      character(len=*) :: groupname
    
354   format(' ',A,T40,A,A)    !Added for vstrlz function results
451   format(' ',A,T40,F0.2,A)
888   format(' ',A,T40,A)
      write(iuniti,*)    
      write(iuniti,888) trim(groupname)//' Statistics:','ON'
      write(iuniti,354) '  Starting at:',trim(vstrlz(tgroupstat(1)/3600.0,'(f0.3)')),' hrs'
      write(iuniti,451) '  Ending at:',tgroupstat(2)/3600.0,' hrs'
      write(iuniti,451) '  Duration:',(tgroupstat(2)-tgroupstat(1))/3600.0,' hrs'
      write(iuniti,451) '  Update Interval:',tgroupstat(3)/3600.0,' hrs'
    
      return
      endsubroutine stat_print_group 
    endsubroutine stat_print
    
!***************************************************************************
    subroutine stat_update()
!***************************************************************************    
    use stat_def
    use diag_lib
    implicit none
    logical :: write_extraLines = .false.
    
    if(flowstats.or.sedstats.or.salstats.or.wavestats) call diag_print_message('')    !Nicer output with empty line before stats are output. MEB  01/07/2020
    if(flowstats) call stat_update_flow  !Flow statistics
    if(sedstats)  call stat_update_sed   !Sediment statistics
    if(salstats)  call stat_update_sal   !Salinity statistics
    if(wavestats) call stat_update_wave  !Wave statistics
    
    return
    endsubroutine stat_update
    
!***************************************************************************
    subroutine stat_update_flow()
! Calculate the flow statistics
!
! Revision list:
!   2013-06-08 - Renamed average residuals and added maximum residuals
!   2014-08-22 - Added minimum and mean water levels
!
! written by Alex Sanchez, USACE=CHL
!***************************************************************************
#include "CMS_cpp.h"
    use comvarbl, only: ctime,timehrs,dtime,timesecs
    use der_def, only: gow
    use der_lib, only: d2xy2d
    use flow_def
    use flow_lib, only: streamwise_curvature
    use geo_def, only: x,y,dx,dy,areap,cell2cell,ncface,zb0
    use prec_def
    use sed_def, only: sedtrans
    use size_def
    use stat_def
    implicit none
    integer :: i
    real(ikind) :: rdninv,val,d2varx2,d2vary2,c,depavg
    !real(ikind),allocatable :: temp(:)
    !real(ikind),allocatable :: dzbdx(:),dzbdy(:)     !Initial bed slopes
    !real(ikind),allocatable :: d2zbdx2(:),d2zbdy2(:) !Initial bed curvatures
    
    if((ctime+1.0e-4)<tflowstat(1))then
      !if(abs(ctime-tflowstat(1))<1.0e-4) call stat_write_flow
      return
    endif
    !if(mod(ctime,tflowstat(3))>1.e-4) return
    
!    write(*,*) 'Calculating Flow Statistics'    
    nflowstatcount = nflowstatcount + 1

!!$OMP PARALLEL 
!!$OMP DO PRIVATE(i,val,c)
    !Hydrodynamics
    do i=1,ncells      
      if(iwet(i)==0) cycle
      
      !Maximum current and wse
      if(eta(i)>etamax(i))then  
        etamax(i) = eta(i)
        tetamax(i) = timehrs
      endif
      if(eta(i)<etamin(i))then
        etamin(i) = eta(i)
        tetamin = timehrs
      endif
      if(uv(i)>uvmax(i))then  
        uvmax(i) = uv(i)
        umax(i) = u(i)
        vmax(i) = v(i)
        tuvmax(i) = timehrs
      endif
      
      !Mean wse
      etamean(i) = etamean(i) + eta(i)
      
      !Residual currents
      ur(i) = ur(i) + u(i)
      vr(i) = vr(i) + v(i)    

      !Streamwise curvature
      uvcurv(i) = uvcurv(i) + streamwise_curvature(u(i),&
                      v(i),uv(i),dux(i),duy(i),dvx(i),dvy(i)) 
              
      !Hydroperiod
      hydper(i) = hydper(i) + iwet(i) 
      
      !Residuals
      rsuavg(i) = rsuavg(i) + rsu(i)
      rsvavg(i) = rsvavg(i) + rsv(i)
      rspavg(i) = rspavg(i) + rsp(i)
      rspmax(i) = max(rspmax(i),rsp(i))
      rsumax(i) = max(rsumax(i),rsu(i))
      rsvmax(i) = max(rsvmax(i),rsv(i))
      
      !Maximum spatial gradients
      val = sqrt(dpx(i)*dpx(i) + dpy(i)*dpy(i))*gravinv
      etamaxgrad(i) = max(etamaxgrad(i),val)  
      val = sqrt(dux(i)*dux(i) + duy(i)*duy(i))
      uvmaxgrad(i) = max(uvmaxgrad(i),val)       
      val = sqrt(dvx(i)*dvx(i) + dvy(i)*dvy(i))
      uvmaxgrad(i) = max(uvmaxgrad(i),val)     
      
      !Celerity
      c = sqrt(grav*h(i))
      
      !!Reynolds numbers
      !rexmax(i) = max(rexmax(i),h(i)*u(i)/vis(i))
      !reymax(i) = max(reymax(i),h(i)*v(i)/vis(i))
      
      !Froude numbers
      frxmax(i) = max(frxmax(i),u(i)/c)
      frymax(i) = max(frymax(i),v(i)/c)
      
      !Flow Curvature
      call d2xy2d(gow,i,dux,duy,d2varx2,d2vary2)
      umaxcurvx(i) = max(umaxcurvx(i),d2varx2)
      umaxcurvy(i) = max(umaxcurvy(i),d2vary2)
      call d2xy2d(gow,i,dvx,dvy,d2varx2,d2vary2) 
      vmaxcurvx(i) = max(vmaxcurvx(i),d2varx2)
      vmaxcurvy(i) = max(vmaxcurvy(i),d2vary2)
    enddo   
!!$OMP END DO

    !Courant numbers
    if(ncellpoly==0)then
      do i=1,ncells 
        c = sqrt(grav*h(i))
        val = (abs(u(i))+c)*dtime/dx(i)
        cflxmax(i) = max(cflxmax(i),val)
        val = (abs(v(i))+c)*dtime/dy(i)
        cflymax(i) = max(cflymax(i),val)      
      enddo
    endif
    
    !Grid Quality
    if(ncellpoly==0)then
!!$OMP DO PRIVATE(i,val)    
    do i=1,ncells  
      !!Grid quality indicator based on wave celerity
      !val=max(dx(i),dy(i))/sqrt(grav*h(i))
      !qcelmax(i)=max(qcelmax(i),val)
      !!Grid quality indicator based on flow depth gradients
      !val=abs(dhx(i))*dx(i) !/h(i)
      !qhgradxmax(i)=max(qhgradxmax(i),val)
      !val=abs(dhy(i))*dy(i) !/h(i)
      !qhgradymax(i)=max(qhgradymax(i),val)
      !!Grid quality indicator based on flow depth curvature
      !call d2xy2d(i,dhx,dhy,d2varx2,d2vary2)
      !val=abs(d2varx2)*dx(i)**2 !/h(i)
      !qhcurvxmax(i)=max(qhcurvxmax(i),val)
      !val=abs(d2vary2)*dy(i)**2 !/h(i)
      !qhcurvymax(i)=max(qhcurvymax(i),val)
      !Grid quality indicators based on courant numbers
      val=dx(i)/dtime/(abs(u(i))+sqrt(grav*h(i)))
      qcflxmax(i)=max(qcflxmax(i),val)
      val=dy(i)/dtime/(abs(v(i))+sqrt(grav*h(i)))
      qcflymax(i)=max(qcflymax(i),val)
      !Grid quality indicators based on Froude numbers
      val=dx(i)/dtime*abs(u(i))/sqrt(grav*h(i))
      qFrxmax(i)=max(qFrxmax(i),val)
      val=dy(i)/dtime*abs(v(i))/sqrt(grav*h(i))
      qFrymax(i)=max(qFrymax(i),val) 
    enddo    
!!$OMP END DO
    endif
!!$OMP END PARALLEL

    !If final count do averaging and unit conversions
    !if(abs(ctime-tflowstat(2))>1.0e-4) return
    if(mod(timesecs-tflowstat(1),tflowstat(2)-tflowstat(1))>1.0e-4) return
    
    !flowstats = .false. !Not needed anymore
    
    rdninv = 1.0/max(float(nflowstatcount),1.0)

    depavg = sum(h(1:ncells))/real(ncells)
    !Hydrodynamics
!$OMP PARALLEL DO PRIVATE(i)    
    do i=1,ncells
      !Mean wse
      etamean(i) = etamean(i)*rdninv ![m]
        
      !Residual currents
      ur(i) = ur(i)*rdninv  ![m/s]
      vr(i) = vr(i)*rdninv  ![m/s]
      uvr(i) = sqrt(ur(i)*ur(i) + vr(i)*vr(i)) ![m/s]    
      
      !Curvature
      uvcurv(i) = uvcurv(i)*rdninv
      
      !Hydroperiod  
      hydper(i) = hydper(i)*rdninv    
      
      !Residuals
      rsuavg(i) = rsuavg(i)*rdninv
      rsvavg(i) = rsvavg(i)*rdninv
      rspavg(i) = rspavg(i)*rdninv
      
      !!Scaling of grid quality variables
      !!Grid quality indicator based on flow depth gradients
      !qhgradxmax(i)=qhgradxmax(i)/depavg
      !qhgradymax(i)=qhgradymax(i)/depavg
      !!Grid quality indicator based on flow depth curvature
      !qhcurvxmax(i)=qhcurvxmax(i)/depavg
      !qhcurvymax(i)=qhcurvymax(i)/depavg
      
      !!Grid quality indicator based on flow depth gradients
      !qhgradxmax(i)=log10(1.0+qhgradxmax(i))
      !qhgradymax(i)=log10(1.0+qhgradymax(i))
      !!Grid quality indicator based on flow depth curvature
      !qhcurvxmax(i)=log10(1.0+qhcurvxmax(i))
      !qhcurvymax(i)=log10(1.0+qhcurvymax(i))
    enddo    
!$OMP END PARALLEL DO
    
    !Scale to make units more visible for display
    !rsuavg = log10(1.0+rsuavg)
    !rsvavg = log10(1.0+rsvavg)
    !rspavg = log10(1.0+rspavg)
    
    !Write Output
    call stat_write_flow
    
    !Re-initialize
    nflowstatcount = 0
!$OMP PARALLEL DO PRIVATE(i)
    do i=1,ncells
      !Maximum flow values  
      uvmax(i) = 0.0
      umax(i) = 0.0; vmax(i) = 0.0
      etamax(i) = 0.0
      etamin(i) = 1.0e20
    
      !Mean wse
      etamean(i) = 0.0
      
      !Residual flow velocities
      uvr(i) = 0.0
      ur(i) = 0.0; vr(i) = 0.0
    
      !Curvature
      uvcurv(i) = 0.0
      
      !Hydro period
      hydper(i) = 0.0    
      
      !Flow Gradients
      etamaxgrad(i) = 0.0
      uvmaxgrad(i) = 0.0      
      
      !Average Residuals
      rspavg(i) = 0.0
      rsuavg(i) = 0.0; rsvavg(i) = 0.0
      
      !Maximum Residuals
      rspmax(i) = 0.0
      rsumax(i) = 0.0; rsvmax(i) = 0.0
      
      !CFL condition
      cflxmax(i) = 0.0; cflymax(i) = 0.0
      
      !!Reynolds numbers
      !rexmax(i) = 0.0; reymax(i) = 0.0
      
      !Froude numbers
      frxmax(i) = 0.0; frymax(i) = 0.0

      !Flow curvature
      umaxcurvx(i) = 0.0; umaxcurvy(i) = 0.0
      vmaxcurvx(i) = 0.0; vmaxcurvy(i) = 0.0   
    enddo
!$OMP END PARALLEL DO

    if(ncellpoly==0)then
!$OMP PARALLEL DO PRIVATE(i)
      do i=1,ncells
        !!Grid quality indicator based on wave celerity  
        !qcelmax(i)=0.0
        !!Grid quality indicator based on flow depth gradients
        !qhgradxmax(i)=0.0
        !qhgradymax(i)=0.0
        !Grid quality indicator based on flow depth curvature
        !qhcurvxmax(i)=0.0
        !qhcurvymax(i)=0.0
        !Grid quality indicators based on courant numbers
        qcflxmax(i)=0.0
        qcflymax(i)=0.0
        !Grid quality indicators based on Froude numbers
        qFrxmax(i)=0.0
        qFrymax(i)=0.0
      enddo
!$OMP END PARALLEL DO
    endif
    
    return
    endsubroutine stat_update_flow
    
!***********************************************************    
    subroutine stat_update_sed()
! Calculate Sediment Statistics 
! written by Alex Sanchez, USACE=CHL
!***********************************************************        
#include "CMS_cpp.h"
    use size_def
    use stat_def
    use geo_def, only: zb
    use flow_def, only: iwet,grav
    use sed_def, only: sedtrans,rhosed,qtx,qty,rsCtkmax
    use comvarbl, only: ctime,timehrs,timesecs
    use prec_def
    implicit none
    integer :: i
    real(ikind) :: rdninv,fac
    !real(ikind) :: dzbdx(ncellsD),dzbdy(ncellsD)     !Bed slopes
    !real(ikind) :: d2zbdx2(ncellsD),d2zbdy2(ncellsD) !Bed curvatures
    
    if((ctime+1.0e-4)<=tsedstat(1))then
      !if(abs(ctime-tsedstat(1))<1.0e-4) call stat_write_sed  
      return 
    endif
    !if(mod(ctime,tsedstat(3))>1.e-4) return
    
    nsedstatcount = nsedstatcount + 1
        
    !Sediment transport
!$OMP PARALLEL DO PRIVATE(i)
    do i=1,ncells
      if(iwet(i)==0) cycle
      
      !Maximum transport rate, kg/m/s
      if(abs(qtx(i))>abs(qtxmax(i))) qtxmax(i) = qtx(i)
      if(abs(qty(i))>abs(qtymax(i))) qtymax(i) = qty(i)      
      
      !Net (Average) transport rates, kg/m/s
      qtxn(i) = qtxn(i) + qtx(i)
      qtyn(i) = qtyn(i) + qty(i)   
      
      !Gross transport rates
      qtg(i) = qtg(i) + sqrt(qtx(i)**2 + qty(i)**2)
      
      !Residual (mean of maximums)
      rsCtm(i) = rsCtm(i) + rsCtkmax(i)
    enddo   
!$OMP END PARALLEL DO

    !If final count do averaging and unit conversions
    !if(abs(ctime-tsedstat(2))>1.0e-5) return
    if(mod(timesecs-tsedstat(1),tsedstat(2)-tsedstat(1))>1.0e-4) return
    
    !sedstats = .false. !Not needed anymore
    
    rdninv = 1.0/max(float(nsedstatcount),1.0)
    
    !Averaging is done by multiplying by the time increment then dividing the time interval
    fac = tsedstat(3)/(tsedstat(2)-tsedstat(1))
       
    !Sediment Transport
!$OMP PARALLEL DO PRIVATE(i)    
    do i=1,ncells
      if(iwet(i)==0) cycle
      
      !Maximum transport rate [kg/m/s]
      qtmax(i) = sqrt(qtxmax(i)*qtxmax(i) + qtymax(i)*qtymax(i)) 
    
      !Average sediment transport [kg/m/s]
      qtxn(i) = fac*qtxn(i)
      qtyn(i) = fac*qtyn(i)
      qtn(i) = sqrt(qtxn(i)*qtxn(i) + qtyn(i)*qtyn(i))                  
      
      !Gross sediment transport [kg/m/s]
      qtg(i) =  fac*qtg(i)
      
      !Residual (mean of maximums)
      rsCtm(i) = rsCtm(i)*rdninv
    enddo
!$OMP END PARALLEL DO

    !Write sediment statistics
    call stat_write_sed

    !Re-initialize
    nsedstatcount = 0
!$OMP PARALLEL DO PRIVATE(i)
    do i=1,ncells
      qtmax(i) = 0.0  !Maximum transport rate [kg/m/s]
      qtxmax(i) = 0.0 !Maximum transport rate in x-direction [kg/m/s]
      qtymax(i) = 0.0 !Maximum transport rate in y-direction [kg/m/s]
      qtn(i) = 0.0    !Net transport rate [kg/m/s]
      qtxn(i) = 0.0   !Net transport rate in x-direction [kg/m/s]
      qtyn(i) = 0.0   !Net transport rate in y-direction [kg/m/s]
      qtg(i) = 0.0    !Gross transport rates
      rsCtm(i) = 0.0  !Residual (mean of maximums)
    enddo   
!$OMP END PARALLEL DO

    return
    endsubroutine stat_update_sed
    
!***********************************************************    
    subroutine stat_update_sal()
! Calculate Salinity Statistics 
! written by Alex Sanchez, USACE=CHL
!***********************************************************        
#include "CMS_cpp.h"
    use size_def, only: ncells,ncellsD
    use stat_def
    use flow_def, only: iwet,grav
    use sal_def
    use comvarbl, only: ctime,timehrs,timesecs
    use prec_def
    implicit none
    integer :: i
    real(ikind) :: rdninv
    
    if((ctime+1.0e-4)<tsalstat(1))then
      !if(abs(ctime-tsalstat(1))<1.0e-4) call stat_write_sal
      return
    endif
    !if(mod(ctime,tsalstat(3))>1.e-4) return
    
    nsalstatcount = nsalstatcount + 1
    
!$OMP PARALLEL DO PRIVATE(i)    
    do i=1,ncells
      salmax(i) = max(salmax(i),sal(i))
      salmin(i) = min(salmin(i),sal(i))
      salavg(i) = salavg(i) + sal(i)
      rssalm(i) = rssalm(i) + rssal(i) !Residual (mean of maximums)
    enddo
!$OMP END PARALLEL DO

    !If final count do averaging and unit conversions
    !if(abs(ctime-tsalstat(2))>1.0e-5) return
    if(mod(timesecs-tsalstat(1),tsalstat(2)-tsalstat(1))>1.0e-4) return
    
    !salstats = .false. !Not needed anymore
    
    rdninv = 1.0/max(float(nsalstatcount),1.0)
    
!$OMP PARALLEL DO PRIVATE(i)
    do i=1,ncells               
      salavg(i)=salavg(i)*rdninv
      rssalm(i)=rssalm(i)*rdninv
    enddo          
!$OMP END PARALLEL DO

    !Write output    
    call stat_write_sal  
    
    !Re-initializie
    nsalstatcount = 0
!$OMP PARALLEL DO PRIVATE(i)
    do i=1,ncells
      salmax(i) = 0.0
      salmin(i) = 1.0e20 
      salavg(i) = 0.0
      rssalm(i) = 0.0
    enddo
!$OMP END PARALLEL DO
    
    return
    endsubroutine stat_update_sal
    
!***********************************************************    
    subroutine stat_update_wave
! Calculate Wave Statistics 
! written by Alex Sanchez, USACE=CHL
!***********************************************************     
#include "CMS_cpp.h"
    use size_def
    use flow_def
    use stat_def
    use wave_flowgrid_def, only: whgt,wlen,wavediss,wper,wunitx,wunity
    use comvarbl, only: ctime,timehrs,timesecs
    use prec_def   
    implicit none
    integer :: i
    real(ikind) :: ursell,rdninv
    
    if((ctime+1.0e-4)<twavstat(1))then
      !if(abs(ctime-twavstat(1))<1.0e-4) call stat_write_wave  
      return
    endif
    !if(mod(ctime,twavstat(3))>1.e-4) return
    
!    write(*,*) 'Calculating Flow Statistics'    
    nwavestatcount = nwavestatcount + 1
        
!$OMP PARALLEL DO PRIVATE(i,Ursell)
    do i=1,ncells
      if(iwet(i)==0) cycle  
      Ursell = whgt(i)*wlen(i)**2/h(i)**3 !Ur = H*L^2/h^3  
      ursellmax(i) = max(ursellmax(i),Ursell) 
      whgtmax(i)   = max(whgtmax(i),whgt(i))
      wdissmax(i)  = max(wdissmax(i),wavediss(i))
      whgtmeanx(i) = whgtmeanx(i) + whgt(i)*wunitx(i)
      whgtmeany(i) = whgtmeany(i) + whgt(i)*wunity(i)
      wpermean(i)  = wpermean(i) + wper(i)
    enddo   
!$OMP END PARALLEL DO
    
    !If final count do averaging and unit conversions
    !if(abs(ctime-twavstat(2))>1.0e-4) return
    if(mod(timesecs-twavstat(1),twavstat(2)-twavstat(1))>1.0e-4) return
    
    !wavestats = .false. !Not needed anymore
    
    rdninv = 1.0/max(float(nwavestatcount),1.0)
    
!$OMP PARALLEL DO PRIVATE(i)
    do i=1,ncells
      if(iwet(i)==0) cycle  
      whgtmeanx(i) = whgtmeanx(i)*rdninv
      whgtmeany(i) = whgtmeany(i)*rdninv
      whgtmean(i) = sqrt(whgtmeanx(i)**2 + whgtmeany(i)**2)
      wpermean(i)  = wpermean(i)*rdninv
    enddo
!$OMP END PARALLEL DO

    !Write Output
    call stat_write_wave
    
    !Re-initialize
    nwavestatcount = 0
!$OMP PARALLEL DO PRIVATE(i)
    do i=1,ncells
      ursellmax(i) = 0.0
      whgtmax(i) = 0.0
      wdissmax(i) = 0.0
      whgtmeanx(i) = 0.0
      whgtmeany(i) = 0.0
      wpermean(i) = 0.0
      whgtmean(i) = 0.0
    enddo
!$OMP END PARALLEL DO

    return
    endsubroutine stat_update_wave    
    
!*****************************************************************************************
    subroutine stat_write_flow()
! Writes the flow statistics
! Author: Alex Sanchez, USACE-CHL
! Modified for ASCII output by Mitchell Brown, USACE-CHL - 02/27/2018
!*****************************************************************************************
#include "CMS_cpp.h"
    use size_def, only: ncellpoly
    use diag_def
    use diag_lib
    use stat_def
    use comvarbl, only: timehrs
    use prec_def
#ifdef XMDF_IO
    use out_lib, only: writescalh5,writevech5
#endif
    use out_lib, only: writescalTxt,writevecTxt
    use out_def, only: write_xmdf_output

    implicit none
    character(len=10) aext
    character(len=200) dirpath
    
787 format(A,T40,F7.2,A) 
    
    !call fileext(statfile,aext)
    if(write_xmdf_output)then          !If no XMDF output, then force to output the Stats in ASCII.
      aext='h5        '
      write(msg2,787) 'Writing XMDF Flow Statistics at: ',timehrs,' hrs'
    else
      aext='txt       '
      write(msg2,787) 'Writing ASCII Flow Statistics at: ',timehrs,' hrs'
    endif
    if(write_xmdf_output .or. timehrs .ne. 0.d0) then
    call diag_print_message(msg2)
    endif  
    
    select case(aext)
    case ('h5')
#ifdef XMDF_IO
    !Water levels
      call writescalh5(statfile,statpath,'WSE_Max',etamax,'m',timehrs,1)
      call writescalh5(statfile,statpath,'WSE_Min',etamin,'m',timehrs,1)
      call writescalh5(statfile,statpath,'WSE_Max_Time',tetamax,'hrs',timehrs,1)
      call writescalh5(statfile,statpath,'WSE_Min_Time',tetamin,'hrs',timehrs,1)
      call writescalh5(statfile,statpath,'WSE_Avg',etamean,'m',timehrs,1)
    call writescalh5(statfile,statpath,'Hydroperiod',hydper,'none',timehrs,1)
      call writescalh5(statfile,statpath,'WSE_Grad_Max_Mag',etamaxgrad,'none',timehrs,1)
    !Current Velocities
      call writevech5 (statfile,statpath,'Current_Vel_Max',umax,vmax,'m/s',timehrs,1)
      call writescalh5(statfile,statpath,'Current_Vel_Max_Mag',uvmax,'m/s',timehrs,1)
      call writescalh5(statfile,statpath,'Current_Vel_Max_Time',tuvmax,'m/s',timehrs,1)
      !call writescalh5(statfile,statpath,'Current_Vel_X_Max',umax,'m/s',timehrs,1)
      !call writescalh5(statfile,statpath,'Current_Vel_Y_Max',vmax,'m/s',timehrs,1)
      !call writescalh5(statfile,statpath,'Current_Vel_X_Max_Time',tumax,'m/s',timehrs,1)
      !call writescalh5(statfile,statpath,'Current_Vel_Y_Max_Time',tvmax,'m/s',timehrs,1)
    call writevech5(statfile,statpath,'Current_Residual',ur,vr,'m/s',timehrs,1)
    call writescalh5(statfile,statpath,'Current_Residual_Mag',uvr,'m/s',timehrs,1)     
    call writescalh5(statfile,statpath,'Flow_Curvature',uvcurv,'1/m',timehrs,1)
      call writescalh5(statfile,statpath,'Current_Grad_Max_Mag',uvmaxgrad,'1/s',timehrs,1)
    
    if(ncellpoly==0)then
      call writescalh5(statfile,statpath,'Grid_Qual_Smoothness',qareachg,'none',timehrs,1)
      !call writescalh5(statfile,statpath,'Grid_Qual_Celerity',qcelmax,'none',timehrs,1)
      call writescalh5(statfile,statpath,'Grid_Qual_Courant_X',qcflxmax,'none',timehrs,1)
      call writescalh5(statfile,statpath,'Grid_Qual_Courant_Y',qcflymax,'none',timehrs,1)
      call writescalh5(statfile,statpath,'Grid_Qual_Froude_X',qFrxmax,'none',timehrs,1)
      call writescalh5(statfile,statpath,'Grid_Qual_Froude_Y',qFrymax,'none',timehrs,1)
      !call writescalh5(statfile,statpath,'Grid_Qual_Flow_Depth_Grad_X',qhgradxmax,'none',timehrs,1)
      !call writescalh5(statfile,statpath,'Grid_Qual_Flow_Depth_Grad_Y',qhgradymax,'none',timehrs,1)
      !call writescalh5(statfile,statpath,'Grid_Qual_Flow_Depth_Curv_X',qhcurvxmax,'none',timehrs,1)
      !call writescalh5(statfile,statpath,'Grid_Qual_Flow_Depth_Curv_Y',qhcurvymax,'none',timehrs,1)
    endif
    
    call writescalh5(statfile,statpath,'Vx_Norm_Res_Avg',rsuavg,'none',timehrs,1)
    call writescalh5(statfile,statpath,'Vy_Norm_Res_Avg',rsvavg,'none',timehrs,1)
    call writescalh5(statfile,statpath,'Pres_Norm_Res_Avg',rspavg,'none',timehrs,1)
    call writescalh5(statfile,statpath,'Vx_Norm_Res_Max',rsumax,'none',timehrs,1)
    call writescalh5(statfile,statpath,'Vy_Norm_Res_Max',rsvmax,'none',timehrs,1)
    call writescalh5(statfile,statpath,'Pres_Norm_Res_Max',rspmax,'none',timehrs,1)
    !call writescalh5(statfile,statpath,'Vx_Norm_Res_Scaled',rsuavg,'none',timehrs,1)
    !call writescalh5(statfile,statpath,'Vy_Norm_Res_Scaled',rsvavg,'none',timehrs,1)
    !call writescalh5(statfile,statpath,'Pres_Norm_Res_Scaled',rspavg,'none',timehrs,1)
    
    if(ncellpoly==0)then
      call writescalh5(statfile,statpath,'Courant_Max_X',cflxmax,'none',timehrs,1)
      call writescalh5(statfile,statpath,'Courant_Max_Y',cflymax,'none',timehrs,1)
      !call writescalh5(statfile,statpath,'Reynolds_Max_X',rexmax,'none',timehrs,1)
      !call writescalh5(statfile,statpath,'Reynolds_Max_Y',reymax,'none',timehrs,1)
      call writescalh5(statfile,statpath,'Froude_Max_X',frxmax,'none',timehrs,1)    
      call writescalh5(statfile,statpath,'Froude_Max_Y',frymax,'none',timehrs,1)
    endif
    
    call writescalh5(statfile,statpath,'dudx_Max',umaxcurvx,'1/s',timehrs,1)
    call writescalh5(statfile,statpath,'dudy_Max',umaxcurvy,'1/s',timehrs,1)
    call writescalh5(statfile,statpath,'dvdx_Max',vmaxcurvx,'1/s',timehrs,1)
    call writescalh5(statfile,statpath,'dvdy_Max',vmaxcurvy,'1/s',timehrs,1)
#endif
    case('txt')
       dirpath='Statistics' 
      
      !Only write ASCII output at end, not the beginning. 
      if(timehrs.ne.0.d0) then
        !Water levels
        call writescalTxt(statfile,dirpath,'WSE_Max',etamax,1)
        call writescalTxt(statfile,dirpath,'WSE_Min',etamin,1)
        call writescalTxt(statfile,dirpath,'WSE_Max_Time',tetamax,1)
        call writescalTxt(statfile,dirpath,'WSE_Min_Time',tetamin,1)
        call writescalTxt(statfile,dirpath,'WSE_Avg',etamean,1)
        call writescalTxt(statfile,dirpath,'Hydroperiod',hydper,1)
        call writescalTxt(statfile,dirpath,'WSE_Grad_Max_Mag',etamaxgrad,1)
        !Current Velocities
        call writevecTxt (statfile,dirpath,'Current_Vel_Max',umax,vmax,1)
        call writescalTxt(statfile,dirpath,'Current_Vel_Max_Mag',uvmax,1)
        call writescalTxt(statfile,dirpath,'Current_Vel_Max_Time',tuvmax,1)
        call writevecTxt (statfile,dirpath,'Current_Residual',ur,vr,1)
        call writescalTxt(statfile,dirpath,'Current_Residual_Mag',uvr,1)     
        call writescalTxt(statfile,dirpath,'Flow_Curvature',uvcurv,1)
        call writescalTxt(statfile,dirpath,'Current_Grad_Max_Mag',uvmaxgrad,1)
      
        if(ncellpoly==0)then
          call writescalTxt(statfile,dirpath,'Grid_Qual_Smoothness',qareachg,1)
          call writescalTxt(statfile,dirpath,'Grid_Qual_Courant_X',qcflxmax,1)
          call writescalTxt(statfile,dirpath,'Grid_Qual_Courant_Y',qcflymax,1)
          call writescalTxt(statfile,dirpath,'Grid_Qual_Froude_X',qFrxmax,1)
          call writescalTxt(statfile,dirpath,'Grid_Qual_Froude_Y',qFrymax,1)
        endif
      
        call writescalTxt(statfile,dirpath,'Vx_Norm_Res_Avg',rsuavg,1)
        call writescalTxt(statfile,dirpath,'Vy_Norm_Res_Avg',rsvavg,1)
        call writescalTxt(statfile,dirpath,'Pres_Norm_Res_Avg',rspavg,1)
        call writescalTxt(statfile,dirpath,'Vx_Norm_Res_Max',rsumax,1)
        call writescalTxt(statfile,dirpath,'Vy_Norm_Res_Max',rsvmax,1)
        call writescalTxt(statfile,dirpath,'Pres_Norm_Res_Max',rspmax,1)
        
        if(ncellpoly==0)then
          call writescalTxt(statfile,dirpath,'Courant_Max_X',cflxmax,1)
          call writescalTxt(statfile,dirpath,'Courant_Max_Y',cflymax,1)
          call writescalTxt(statfile,dirpath,'Froude_Max_X',frxmax,1)    
          call writescalTxt(statfile,dirpath,'Froude_Max_Y',frymax,1)
        endif
    
        call writescalTxt(statfile,dirpath,'dudx_Max',umaxcurvx,1)
        call writescalTxt(statfile,dirpath,'dudy_Max',umaxcurvy,1)
        call writescalTxt(statfile,dirpath,'dvdx_Max',vmaxcurvx,1)
        call writescalTxt(statfile,dirpath,'dvdy_Max',vmaxcurvy,1)
      endif  
        
    end select
    
    return
    endsubroutine stat_write_flow
    
!********************************************************************************************
    subroutine stat_write_sed()
! Writes the sediment transport statistics
! Author: Alex Sanchez, USACE-CHL
! Modified for ASCII output by Mitchell Brown, USACE-CHL - 02/27/2018
!********************************************************************************************
#include "CMS_cpp.h"
    use size_def, only: ncellsD
    use stat_def
    use diag_def
    use diag_lib
    use comvarbl, only: timehrs
#ifdef XMDF_IO
    use out_lib, only: writescalh5,writevech5
#endif
    use out_lib, only: writescalTxt,writevecTxt
    use out_def, only: write_xmdf_output

    implicit none
    character(len=10) aext
    character(len=200) dirpath
    real(ikind) :: ones(ncellsD)

787 format(A,T40,F7.2,A) 
    
    !call fileext(statfile,aext)
    if(write_xmdf_output)then          !If no XMDF output, then force to output the Stats in ASCII.
      aext='h5        '
      write(msg2,787) 'Writing XMDF Sediment Statistics at: ',timehrs,' hrs'
    else
      aext='txt       '
      write(msg2,787) 'Writing ASCII Sediment Statistics at: ',timehrs,' hrs'
    endif
    if(write_xmdf_output .or. timehrs .ne. 0.d0) then
    call diag_print_message(msg2)
    endif  
    
    ones = 1.0
    
    select case(aext)
    case ('h5')
#ifdef XMDF_IO
    call writescalh5(statfile,statpath,'Ones',ones,'none',timehrs,1) !Used in SMS for flux calculations
    call writevech5(statfile,statpath,'Sed_Transp_Max',qtxmax,qtymax,'m/s',timehrs,1)
    call writescalh5(statfile,statpath,'Sed_Transp_Max_Mag',qtmax,'kg/m/s',timehrs,1)
    call writevech5(statfile,statpath,'Sed_Transp_Net',qtxn,qtyn,'kg/m/s',timehrs,1)
    call writescalh5(statfile,statpath,'Sed_Transp_Net_Mag',qtn,'kg/m/s',timehrs,1)
    call writescalh5(statfile,statpath,'Sed_Transp_Gross_Mag',qtg,'kg/m/s',timehrs,1)
    call writescalh5(statfile,statpath,'Conc_Norm_Res_Avg',rsCtm,'none',timehrs,1)
    !call writescalh5(statfile,statpath,'Conc_Norm_Res_Max',rsCtmax,'none',timehrs,1)
#endif
    case('txt')
      dirpath='Statistics' 
      
      !Only write ASCII output at end, not the beginning. 
      if(timehrs.ne.0.d0) then
        call writescalTxt(statfile,dirpath,'Ones',ones,1) !Used in SMS for flux calculations
        call writevecTxt (statfile,dirpath,'Sed_Transp_Max',qtxmax,qtymax,1)
        call writescalTxt(statfile,dirpath,'Sed_Transp_Max_Mag',qtmax,1)
        call writevecTxt (statfile,dirpath,'Sed_Transp_Net',qtxn,qtyn,1)
        call writescalTxt(statfile,dirpath,'Sed_Transp_Net_Mag',qtn,1)
        call writescalTxt(statfile,dirpath,'Sed_Transp_Gross_Mag',qtg,1)
        call writescalTxt(statfile,dirpath,'Conc_Norm_Res_Avg',rsCtm,1)
      endif    
    end select
    
    return
    endsubroutine stat_write_sed
    
!****************************************************************************************
    subroutine stat_write_sal()
! Writes the salinity transport statistics
! Author: Alex Sanchez, USACE-CHL
! Modified for ASCII output by Mitchell Brown, USACE-CHL - 02/27/2018
!***************************************************************************************
#include "CMS_cpp.h"
    use stat_def
    use diag_def
    use diag_lib
    use comvarbl, only: timehrs
#ifdef XMDF_IO
    use out_lib, only: writescalh5
#endif
    use out_lib, only: writescalTxt,writevecTxt
    use out_def, only: write_xmdf_output

    implicit none
    character(len=10) aext
    character(len=200) dirpath

787 format(A,T40,F7.2,A) 
    
    if(write_xmdf_output)then          !If no XMDF output, then force to output the Stats in ASCII.
      aext='h5        '
      write(msg2,787) 'Writing XMDF Salinity Statistics at: ',timehrs,' hrs'
    else
      aext='txt       '
      write(msg2,787) 'Writing ASCII Salinity Statistics at: ',timehrs,' hrs'
    endif
    if(write_xmdf_output .or. timehrs .ne. 0.d0) then
    call diag_print_message(msg2)
    endif  
    
    select case(aext)
    case ('h5')
#ifdef XMDF_IO
    call writescalh5(statfile,statpath,'Salinity_Max',salmax,'ppt',timehrs,1)
    call writescalh5(statfile,statpath,'Salinity_Min',salmin,'ppt',timehrs,1)
    call writescalh5(statfile,statpath,'Salinity_Avg',salavg,'ppt',timehrs,1)
    call writescalh5(statfile,statpath,'Salinity_Norm_Res',rssalm,'none',timehrs,1) 
#endif 
    case ('txt')
      dirpath='Statistics' 

      !Only write ASCII output at end, not the beginning. 
      if(timehrs.ne.0.d0) then
        call writescalTxt(statfile,dirpath,'Salinity_Max',salmax,1)
        call writescalTxt(statfile,dirpath,'Salinity_Min',salmin,1)
        call writescalTxt(statfile,dirpath,'Salinity_Avg',salavg,1)
        call writescalTxt(statfile,dirpath,'Salinity_Norm_Res',rssalm,1) 
      endif
    end select 
    return
    endsubroutine stat_write_sal
    
!*****************************************************************************************    
    subroutine stat_write_wave()
! Writes the wave statistics
! Author: Alex Sanchez, USACE-CHL
! Modified for ASCII output by Mitchell Brown, USACE-CHL - 02/27/2018
!*****************************************************************************************
#include "CMS_cpp.h"
    use stat_def
    use diag_def
    use diag_lib
    use comvarbl, only: timehrs
#ifdef XMDF_IO
    use out_lib, only: writescalh5,writevech5
#endif
    use out_lib, only: writescalTxt,writevecTxt
    use out_def, only: write_xmdf_output

    implicit none
    character(len=10) aext
    character(len=200) dirpath

787 format(A,T40,F7.2,A) 
    
    if(write_xmdf_output)then          !If no XMDF output, then force to output the Stats in ASCII.
      aext='h5        '
      write(msg2,787) 'Writing XMDF Wave Statistics at: ',timehrs,' hrs'
    else
      aext='txt       '
      write(msg2,787) 'Writing ASCII Wave Statistics at: ',timehrs,' hrs'
    endif
    if(write_xmdf_output .or. timehrs .ne. 0.d0) then
    call diag_print_message(msg2)
    endif  
    
    select case (aext)
    case('h5')
#ifdef XMDF_IO
    call writescalh5(statfile,statpath,'Wave_Height_Max',whgtmax,'m',timehrs,1)
    call writescalh5(statfile,statpath,'Wave_Breaking_Max',wdissmax,'none',timehrs,1)          
    call writescalh5(statfile,statpath,'Ursell_Max',ursellmax,'none',timehrs,1)    
    call writescalh5(statfile,statpath,'Wave_Height_Mean',whgtmean,'m',timehrs,1)
    call writevech5(statfile,statpath,'Wave_Height_Mean_Vec',whgtmeanx,whgtmeany,'m/s',timehrs,1)
    call writescalh5(statfile,statpath,'Wave_Period_Mean',wpermean,'sec',timehrs,1)
    !call writevech5(statfile,statpath,'Bed_Slope_Max',bedmaxslpx,bedmaxslpy,'none',timehrs,1)
#endif 
    case('txt')
      dirpath='Statistics' 
      
      !Only write ASCII output at end, not the beginning. 
      if(timehrs.ne.0.d0) then
        call writescalTxt(statfile,dirpath,'Wave_Height_Max',whgtmax,1)
        call writescalTxt(statfile,dirpath,'Wave_Breaking_Max',wdissmax,1)          
        call writescalTxt(statfile,dirpath,'Ursell_Max',ursellmax,1)    
        call writescalTxt(statfile,dirpath,'Wave_Height_Mean',whgtmean,1)
        call writevecTxt (statfile,dirpath,'Wave_Height_Mean_Vec',whgtmeanx,whgtmeany,1)
        call writescalTxt(statfile,dirpath,'Wave_Period_Mean',wpermean,1)
      endif        
    end select

    return
    endsubroutine stat_write_wave
    
!***********************************************************    
    subroutine stat_read()
! Reads the global statistics from the hotstart file
! written by Alex Sanchez, USACE-CHL
!***********************************************************    
#include "CMS_cpp.h"
    use size_def
    use stat_def
    use geo_def, only: zb
    use flow_def, only: eta
    use comvarbl, only: timehrs
    use sed_def, only: sedtrans
    use diag_lib
#ifdef XMDF_IO
    use in_xmdf_lib, only: readscalh5,readvech5
    use xmdf
#endif
    use in_lib, only: readscalTxt,readvecTxt
    use out_def, only: simlabel
    
    implicit none    
    integer :: nn,error,fid,gid,kunit
    character(len=200) :: apath
    character(len=10)  :: aext
     
    call fileext(trim(statfile),aext)
    select case (aext)
      case('h5')
#ifdef XMDF_IO  
    nn=len_trim(statpath)    
    call XF_OPEN_FILE(trim(statfile),READONLY,fid,error)   
    if(error<0)then  
      call diag_print_warning('Could not open simulation statistics file: ',statfile)
      return
    endif
    call XF_OPEN_GROUP(fid,trim(simlabel),gid,error)
    if(error<0)then  
      call diag_print_warning('Invalid dataset in simulation statistics file: ',statfile,&
         '  Path: ',apath)   
      return
    endif
    apath = statpath(1:nn) // 'Current_Max'   
    call readvech5(statfile,apath,umax,umax,error)
    
    apath = statpath(1:nn) // 'Current_Max_Mag'  
    call readscalh5(statfile,apath,uvmax,error)       
              
    apath = statpath(1:nn) // 'Water_Level_Max'  
    call readscalh5(statfile,apath,uvmax,error)      
    
    apath = statpath(1:nn) // 'Current_Residual'  
    call readvech5(statfile,apath,ur,vr,error)
    
    apath = statpath(1:nn) // 'Current_Residual_Mag'  
    call readscalh5(statfile,apath,uvr,error)
    
    apath = statpath(1:nn) // 'Hydroperiod'  
    call readscalh5(statfile,apath,hydper,error)
    
    apath = statpath(1:nn) // 'WSE_Max_Grad_Mag'  
    call readscalh5(statfile,apath,etamaxgrad,error)
    
    apath = statpath(1:nn) // 'Cur_Max_Grad_Mag'  
    call readscalh5(statfile,apath,uvmaxgrad,error)
    
    !apath = statpath(1:nn) // 'Bed_Max_Grad_mag'  
    !call readscalh5(statfile,apath,bedmaxgrad,error)
    
    if(.not.sedtrans) return   
    apath = statpath(1:nn) // 'Sed_Transp_Net'  
    call readvech5(statfile,apath,qtxn,qtyn,error)
    
    apath = statpath(1:nn) // 'Sed_Transp_Net_Mag'  
    call readscalh5(statfile,apath,qtn,error)
    
    apath = statpath(1:nn) // 'Sed_Transp_Gross_Mag'  
    call readscalh5(statfile,apath,qtg,error)       
    
    apath = statpath(1:nn) // 'Sed_Transp_Max'  
    call readvech5(statfile,apath,qtxmax,qtymax,error)
#endif
      case('txt')  
        kunit=600  
        open(kunit,FILE=statfile)
        
        call readvecTxt(statfile,umax,umax,error)      !'Current_Max'   
        call readscalTxt(statfile,uvmax,error)         !'Current_Max_Mag'  
        call readscalTxt(statfile,uvmax,error)         !'Water_Level_Max'  
        call readvecTxt(statfile,ur,vr,error)          !'Current_Residual'  
        call readscalTxt(statfile,uvr,error)           !'Current_Residual_Mag'  
        call readscalTxt(statfile,hydper,error)        !'Hydroperiod'  
        call readscalTxt(statfile,etamaxgrad,error)    !'WSE_Max_Grad_Mag'  
        call readscalTxt(statfile,uvmaxgrad,error)     !'Cur_Max_Grad_Mag'  
    
        if(.not.sedtrans) return   
        call readvecTxt(statfile,qtxn,qtyn,error)      !'Sed_Transp_Net'  
        call readscalTxt(statfile,qtn,error)           !'Sed_Transp_Net_Mag'  
        call readscalTxt(statfile,qtg,error)           !'Sed_Transp_Gross_Mag'  
        call readvecTxt(statfile,qtxmax,qtymax,error)  !'Sed_Transp_Max'  
        close(kunit)
      end select   
    
    return
    endsubroutine stat_read
