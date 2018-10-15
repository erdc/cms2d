!***********************************************************************
    subroutine CMS_Flow
! Main Subroutine of CMS2D Implicit Flow Solver
! Made by Weiming Wu, Oct. 2008 
! Modified by Alex Sanchez, USACE-CHL
!***********************************************************************
#include "CMS_cpp.h"
    use sed_def
    use sal_def
    use heat_def
    use hot_def
    use stat_def  
    use out_def
    use cms_def, only: cmswave
    use comvarbl
#ifdef DEV_MODE
    use wave_flowgrid_def, only: constant_waves
    use fric_def, only: mbedfric
    use flow_exp_mod
    use flow_semi_mod
    use geo_def, only: bathydata,zb
#endif    

#ifdef DREDGE
    use dredge_def
#endif    
    implicit none
    
    call prestart    
    
    call sim_start_print  !Start timer here

    !Output
    if(timehrs > timeout)then      
      call write_output !Global variable snapshots
      call stat_update  !Global variable statistics
    endif
    
    do while(ctime < stimet) !Note that stimet includes ramp period
       call timing   !ntime, mtime, dtime, ctime, timehrs, ramp         
       if(cmswave) call wave_eval !Run wave model, and interpolate in space and time
#ifdef DEV_MODE       
       if(constant_waves) call wave_const_update !Set constant wave parameters, for idealized cases and lab experiments
       if(mbedfric == 0) call fric_rough_eval
#endif 

!Temporary to match Chris' Files for choosing Explicit/Implicit
!!       if(nfsch==0)then !Implicit                     !Temporarily commented
         call flow_imp !u,v,p,eta,h,flux,vis,etc
!!       elseif(nfsch==1)then !Explicit                 !Temporarily commented
!!         call flow_exp !u,v,p,eta,h,flux,vis,etc      !Temporarily commented
!!       elseif(nfsch==2)then !Semi-implicit            !Temporarily commented
!!         call flow_semi !u,v,p,eta,h,flux,vis,etc     !Temporarily commented
!!       endif                                          !Temporarily commented
!#ifdef DEV_MODE
!       if(iFlow1D>0) call makeflow1D !Constant flow conditions in column direction, for testing only       
!#endif       
       if(sedtrans) call sed_imp !Sediment transport, implicit solution
       if(saltrans) call sal_imp !Salinity transport, implicit solution
       if(heattrans) call heat_imp  !Heat transfer, implicit solution

#ifdef DREDGE       
       if(dredging)then
         call dredge_eval
         if(sedtrans) call sed_concdepthchange !Correct concentrations for depth changes
       endif
#endif
#ifdef DEV_MODE       
       if(bathydata%ison)then
         call geo_bathy_update(zb,-1)
         call geo_flowdepth_update
       endif
#endif

       if(hot_out) call hot_write
       if(timehrs > timeout)then !Avoids writing before end of last simulation        
         call write_output     !Global variable snapshots
         call stat_update      !Gloval variable statistics
       endif
!#ifdef DEV_MODE
!       if(pred_corr) call flow_pred !Flow prediction based on simple explicit scheme  
!#endif
       call flow_update    !Reset variables for next iterations
    enddo
    if(save_point) call save_point_close
    call sim_end_print
    call cleanup
    
    return
    endsubroutine CMS_Flow
    
!*****************************************************************       
    subroutine timing
! Calculates the cumulative time, number of iterations and
! time step
! written by Alex Sanchez, USACE-CHL
!*****************************************************************
#include "CMS_cpp.h"
    use comvarbl
    use const_def, only: pi
    use time_lib
    use diag_def
    use diag_lib
    use cms_def, only: timestart,timenow
    use solv_def, only: iconv
    use prec_def
    implicit none
    integer :: ierr
    real(8) :: dtimetemp,rtime !,dtimetemp2,dtime2
    real(8) :: timedur,timerem,timelast,speed,err,timeint
    character(len=100) :: str
    
    ntime = ntime + 1  !Time step iteration counter
    
    if (etime > 0 .and. timesecs > 0.0) then      !This should be modified for implicit
      if(ntime==11 .or. mod(ntime,nprt)==0)then
        timelast = timenow
        timenow = time_jul()
        timeint = timenow - timelast  !Time interval between last speed check [sec]
        timedur = timenow - timestart !Total simulation clock time [sec]
        speed = dble(ctime-ctime1)/timeint !Note: computed using last speed check time interval
        ctime1 = ctime
        timerem = dble(stimet-ctime)/speed
        call time_sec2str(timesecs,str)
        write(msg,'(A,A)',iostat=ierr)     'Elapsed Simulation Time:  ',trim(str)
        call diag_print_message(' ',msg)
        call time_sec2str(timedur,str)
        write(msg,'(A,A)',iostat=ierr)     'Elapsed Clock Time:       ',trim(str)
        call diag_print_message(msg)
        write(msg,'(A,F10.3)',iostat=ierr) 'Computational Speed:      ',speed
        call diag_print_message(msg)
        call time_sec2str(timerem,str)
        write(msg,'(A,A)',iostat=ierr)     'Remaining Clock Time:     ',trim(str)
        call diag_print_message(msg)
      endif
    endif
    
!#ifdef DEV_MODE
!    if(ntime<=nspinup)then
!     deltime = dtimebeg*dble(ntime+1)/dble(nspinup+1) !Use double precision
!     dtime = real(deltime,kind=ikind) !Arbitrary precision for numerical compuations
!    else  
!#endif        
     if(iconv==3 .and. dtimebeg/deltime<=10)then
       write(6,*) 'Time Step Reduced'
       mtime = 0
       jtime = jtime+1
       deltime = dtimebeg/2**jtime !Avoids precision errors, Use double precision
       dtime = real(deltime,kind=ikind) !Arbitrary precision for numerical compuations
       !write(*,*) jtime,dtime,dtime2
     elseif(jtime>=1)then
       if((mtime>=2 .and. rmom(1)<rmomtargetp/100.0) .or. &
          (mtime>=4 .and. rmom(1)<rmomtargetp/10.0))then
         dtimetemp = dtimebeg/2**(jtime-1) !Avoids precision errors
         rtime = mod(timesecs+dtimetemp,dtimetemp)
         err = abs(rtime)/dtimetemp
#ifdef DIAG_MODE
         write(*,*) 'rtime: ',rtime,', err: ',err
#endif
         if(rtime<0.01 .or. err<0.01)then
           deltime = dtimetemp 
           dtime = real(deltime,kind=ikind)
           mtime = 0
           jtime = jtime-1 !>=0
           !write(*,*) jtime,dtimetemp,dtimetemp2
         endif
       endif        
     endif
     mtime = mtime + 1           
     timesecs = timesecs + deltime
     ctime = real(timesecs,kind=ikind)
     timehrs = ctime/3600.0_ikind
     ramp = ramp_func(timehrs,rampdur)
!#ifdef DEV_MODE
!    endif
!#endif

    if(ntime>1 .and. wtsch>1.0e-4 .and. abs(dtime-dtime1)<1.0e-4)then
    !if(ntime>1 .and. wtsch>1.0e-4)then
      ntsch = 2
    else
      ntsch = 1
    endif
    
    return
    endsubroutine timing    
    
!*************************************************
    subroutine flow_update
! update field values at different time steps
! by Weiming Wu, NCCHE, Oct. 2008
! modified by Alex Sanchez, USACE-CHL
!*************************************************
#include "CMS_cpp.h"
    use size_def
    use geo_def, only: zb
    use flow_def
    use comvarbl, only: dtime,dtime1,wtsch,ntime,nspinup,nfsch
    use sed_def
    use sal_def
    use heat_def
    use prec_def
    implicit none
    integer :: i

    dtime1 = dtime
    
!#ifdef DEV_MODE
!    if(ntime<=nspinup)then
!!$OMP PARALLEL
!!--- Hydro -------------------
!!!$OMP DO PRIVATE(i)
!!    do i=1,ncellsD
!!      !iwet(i)=iwet1(i)
!!      !flux(:,i)=flux1(:,i)
!!      !u(i)=u1(i) 
!!      !v(i)=v1(i)
!!      p(i)=p1(i)
!!      h(i)=h1(i)
!!      !pp(i)=0.0
!!      !dppx(i)=0.0
!!      !dppy(i)=0.0
!!    enddo
!!!$OMP END DO NOWAIT
!
!!--- Sediment -------------------
!    if(sedtrans)then
!!$OMP DO PRIVATE(i)
!      do i=1,ncellsD        
!        zb(i)=zb1(i)
!        Ctk(i,:)=Ctk1(i,:)
!        if(ibt>0)then
!          btk(i,:)=btk1(i,:)
!        endif  
!        if(.not.singlesize)then        
!          pbk(i,:,1)=pbk1(i,:) !Mixing layer
!          db(i,:)=db1(i,:)     !Bed layer thickness
!        endif      
!      enddo
!!$OMP END DO NOWAIT
!    endif
!    
!!--- Salinity -----------------------      
!    if(saltrans)then
!!$OMP DO PRIVATE(i)      
!      do i=1,ncellsD
!        sal(i)=sal1(i)      
!      enddo
!!$OMP END DO        
!    endif
!!$OMP END PARALLEL
!      return
!    endif
!#endif
    
!$OMP PARALLEL
!--- Hydro -------------------
    if(wtsch>1.0e-4)then
!$OMP DO PRIVATE(i)
      do i=1,ncellsD
        u2(i)=u1(i)
        v2(i)=v1(i)
        p2(i)=p1(i)
        h2(i)=h1(i)
      enddo
!$OMP END DO
    endif
    if(nfsch==0)then !Implicit
!$OMP DO PRIVATE(i)
      do i=1,ncellsD
        iwet1(i)=iwet(i)
        flux1(:,i)=flux(:,i)
        u1(i)=u(i) 
        v1(i)=v(i)
        p1(i)=p(i)
        h1(i)=h(i)
        pp(i)=0.0
        !!ppk(i,:)=0.0
        dppx(i)=0.0
        dppy(i)=0.0
      enddo
!$OMP END DO NOWAIT
    else !Explicit or semi-implicit
!$OMP DO PRIVATE(i)
      do i=1,ncellsD
        iwet1(i)=iwet(i)
        flux1(:,i)=flux(:,i)
        u1(i)=u(i) 
        v1(i)=v(i)
        p1(i)=p(i)
        h1(i)=h(i)
      enddo
!$OMP END DO NOWAIT
    endif
    
!--- Sediment -------------------
    if(sedtrans)then
!$OMP DO PRIVATE(i)
      do i=1,ncellsD        
        zb1(i)=zb(i)
        if(wtsch>1.0e-4) Ctk2(i,:)=Ctk1(i,:)
        Ctk1(i,:)=Ctk(i,:)
        if(ibt>0)then
          if(wtsch>1.0e-4) btk2(i,:)=btk1(i,:)
          btk1(i,:)=btk(i,:)
        endif  
        if(.not.singlesize)then        
          pbk1(i,:)=pbk(i,:,1) !Mixing layer
          db1(i,:)=db(i,:)     !Bed layer thickness
        endif      
      enddo
!$OMP END DO NOWAIT
    endif

!--- Salinity -----------------------      
    if(saltrans)then        
      if(wtsch>1.0e-4)then
!$OMP DO PRIVATE(i)          
        do i=1,ncellsD    
          sal2(i)=sal1(i)
        enddo
!$OMP END DO     
      endif  
!$OMP DO PRIVATE(i)      
      do i=1,ncellsD
        sal1(i)=sal(i)      
      enddo
!$OMP END DO        
    endif  

!--- Temperature -----------------------      
    if(heattrans)then        
      if(wtsch>1.0e-4)then
!$OMP DO PRIVATE(i)          
        do i=1,ncellsD    
          heat2(i)=heat1(i)
        enddo
!$OMP END DO     
      endif  

!$OMP DO PRIVATE(i)      
      do i=1,ncellsD
        heat1(i)=heat(i)      
      enddo
!$OMP END DO        
    endif  
!$OMP END PARALLEL

    return
    endsubroutine flow_update
   
 !*****************************************************
    subroutine makeflow1D
 !*****************************************************
    use geo_def, only: maxcol,maxrow,idmap
    use flow_def, only: u,v,p,uv
    use comvarbl, only: iflow1D
    use prec_def
    implicit none
    integer :: i,j,ii,id
    real(ikind) :: uavg,vavg,pavg,dmaxcol,dmaxrow    
    
    dmaxcol = dble(maxcol)
    dmaxrow = dble(maxrow)
    if(iflow1D==1)then
      do i=1,maxcol
        uavg = 0.0; vavg = 0.0; pavg = 0.0
        do j=1,maxrow        
          id = i + (j-1)*maxcol
          ii = idmap(id)
          uavg = uavg + u(ii)
          vavg = vavg + v(ii)
          pavg = pavg + p(ii)       
        enddo
        uavg = uavg/dmaxrow
        vavg = vavg/dmaxrow
        pavg = pavg/dmaxrow
        do j=1,maxrow   
          id = i + (j-1)*maxcol
          ii = idmap(id)
          u(ii) = uavg
          v(ii) = vavg
          uv(ii) = sqrt(uavg*uavg + vavg*vavg)
          p(ii) = pavg
        enddo
      enddo
    endif
     
    return
    endsubroutine

!***************************************************    
    subroutine cleanup
!***************************************************
    implicit none
    
    !call geo_cleanup
    call flow_cleanup
    !call sed_cleanup
    !call sal_cleanup
    
    return
    endsubroutine cleanup
    