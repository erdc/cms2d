!***********************************************************************
    subroutine CMS_Flow
! Main Subroutine of CMS2D Implicit Flow Solver
! Made by Weiming Wu, Oct. 2008 
! Modified by Alex Sanchez, USACE-CHL
!***********************************************************************
#include "CMS_cpp.h"
    use sed_def,  only: sedtrans
    use sal_def,  only: saltrans
    use heat_def, only: heattrans
    use hot_def,  only: timeout, hot_out
    use out_def,  only: simlabel, write_maxwse, write_xmdf_output, write_sup, outlist, save_point
    use cms_def,  only: cmswave
    use comvarbl, only: timehrs, ctime, stimet, flowpath, casename
    use dredge_def, only: dredging
    use wave_flowgrid_def, only: constant_waves
!! Added MEB  9/20/2021
    use flow_def, only: maxeta
    use out_lib,  only: write_scal_dat_file
#ifdef DEV_MODE
    use fric_def, only: mbedfric
    use flow_exp_mod
    use flow_semi_mod
    use geo_def, only: bathydata,zb
#endif    

#ifdef XMDF_IO
    use out_lib, only: writescalh5
#endif

    implicit none
    character(len=200) :: apath, aname
    
    call prestart    
    call sim_start_print  !Start timer here

    !Output
    if(timehrs > timeout)then      
      call write_output !Global variable snapshots
      call stat_update  !Global variable statistics
    endif
    
    do while(ctime < stimet) !Note that stimet includes ramp period
      call timing   !ntime, mtime, dtime, ctime, timehrs, ramp                        !added scalemorph_ramp - meb 03/11/2019
      if(cmswave) call wave_eval !Run wave model, and interpolate in space and time
      if(constant_waves) call wave_const_update !Set constant wave parameters, for idealized cases and lab experiments

#ifdef DEV_MODE       
      if(mbedfric == 0) call fric_rough_eval
#endif 

!Temporary to match Chris' Files for choosing Explicit/Implicit
!!      if(nfsch==0)then !Implicit                     !Temporarily commented
        call flow_imp !u,v,p,eta,h,flux,vis,etc
!!      elseif(nfsch==1)then !Explicit                 !Temporarily commented
!!        call flow_exp !u,v,p,eta,h,flux,vis,etc      !Temporarily commented
!!      elseif(nfsch==2)then !Semi-implicit            !Temporarily commented
!!        call flow_semi !u,v,p,eta,h,flux,vis,etc     !Temporarily commented
!!      endif                                          !Temporarily commented
!#ifdef DEV_MODE
!      if(iFlow1D>0) call makeflow1D !Constant flow conditions in column direction, for testing only       
!#endif       
      if(sedtrans) call sed_imp !Sediment transport, implicit solution
      if(saltrans) call sal_imp !Salinity transport, implicit solution
      if(heattrans) call heat_imp  !Heat transfer, implicit solution
      if(dredging)then
        call dredge_eval
        if(sedtrans) call sed_concdepthchange !Correct concentrations for depth changes
      endif
#ifdef DEV_MODE       
      if(bathydata%ison)then
        call geo_bathy_update(zb,-1)
        call geo_flowdepth_update
      endif
#endif

      if(hot_out) call hot_write
      if(timehrs > timeout)then !Avoids writing before end of last simulation        
        call write_output     !Global variable snapshots
        call stat_update      !Global variable statistics
      endif

!#ifdef DEV_MODE
!       if(pred_corr) call flow_pred !Flow prediction based on simple explicit scheme  
!#endif
      call flow_update    !Reset variables for next iterations
    enddo

    !added MEB 9/20/2021    Output the Maximum WSE to the correct file.
    if(write_maxwse) then !Maximum WSE 
          
#ifdef XMDF_IO
      apath = trim(simlabel)//'/'
      if(write_xmdf_output) call writescalh5(outlist(1)%afile,apath,'Maximum Water_Elevation',maxeta,'m',timehrs,0)
#endif   
      aname = trim(flowpath) // 'ASCII_Solutions' // '/' // trim(casename)
      if(write_sup) call write_scal_dat_file(aname,'Maximum_Water_Elevation','maxeta',maxeta) !SUPER ASCII File
                        !write_scal_dat_file(aname,'Water_Elevation','eta',eta) !SUPER ASCII File      
    endif
    
    if(save_point) call save_point_close
    call sim_end_print
    call cleanup
    
    return
    end subroutine CMS_Flow
    
!*****************************************************************       
    subroutine timing
! Calculates the cumulative time, number of iterations and
! time step
! written by Alex Sanchez, USACE-CHL
!*****************************************************************
#include "CMS_cpp.h"
    use comvarbl,  only: ntime, etime, timehrs, timesecs, nprt, ctime, ctime1, stimet, dtimebeg, deltime, mtime, jtime, dtime
    use comvarbl,  only: ramp, rampdur, wtsch, dtime1, ntsch, rmom, rmomtargetp
    use const_def, only: pi
    use time_lib,  only: time_jul, ramp_func, time_sec2str
    use diag_def,  only: msg, msg2
    use diag_lib,  only: diag_print_message
    use cms_def,   only: timestart,timenow
    use solv_def,  only: iconv
    use prec_def,  only: ikind
    use sed_def,   only: scalemorph
    use sed_def,   only: scalemorph_rampdur, scalemorph_ramp, scalemorph_orig  !added meb 03/11/2019
    use sed_def,   only: write_smorph_ramp,smorph_file
    
    implicit none
    integer     :: ierr
    real(8)     :: dtimetemp,rtime !,dtimetemp2,dtime2
    real(8)     :: timedur,timerem,timelast,speed,err,timeint
    logical     :: found
    character(len=100) :: str
    
    ntime = ntime + 1  !Time step iteration counter
    
    if (etime > 0 .and. timesecs > 0.0) then      !This should be modified for implicit
      if(ntime==10 .or. mod(ntime,nprt)==0)then
        timelast = timenow
        timenow = time_jul()
        timeint = timenow - timelast  !Time interval between last speed check [sec]
        timedur = timenow - timestart !Total simulation clock time [sec]
        speed = dble(ctime-ctime1)/timeint !Note: computed using last speed check time interval
        ctime1 = ctime
        timerem = dble(stimet-ctime)/speed
        msg=''; msg2=''
        call time_sec2str(timesecs,str,.false.)
        write(msg,'(A,A)',iostat=ierr)       'Elapsed Simulation Time:       ',trim(str)

        ! meb 02/06/2019
        ! Adding an output statement to write out the "Elapsed Morphologic Time" if Morphologic Acceleration Factor > 1
        if (scalemorph .gt. 1) then
          call time_sec2str(timesecs*scalemorph,str,.false.)
          write(msg2,'(A,A)',iostat=ierr)     '- Effective Morphologic Time:  ',trim(str)
        endif   
        !------------
        call diag_print_message(' ',msg,msg2)

        call time_sec2str(timedur,str,.false.)
        write(msg,'(A,A)',iostat=ierr)       'Elapsed Clock Time:            ',trim(str)
        call diag_print_message(msg)

        write(str,'(F0.3)',iostat=ierr)     speed
        write(msg,'(A,A)',iostat=ierr)       'Computational Speed:           ',trim(str)
        call diag_print_message(msg)

        call time_sec2str(timerem,str,.false.)
        write(msg,'(A,A)',iostat=ierr)       'Remaining Clock Time:          ',trim(str)
        call diag_print_message(msg)

        continue
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
      call diag_print_message('Time Step Reduced')
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

!added meb 03/11/2019     
1001 format (2x,f8.4,6x,f10.5)
    if (scalemorph_rampdur .gt. 0.0) then                                        !If there is a ramp, proceed
      scalemorph_ramp = ramp_func(timehrs,scalemorph_rampdur)
      if (scalemorph_orig .gt. 0.0) then                                         !If the original morph accel factor is not 0.0, proceed
        scalemorph = max(1.0 , scalemorph_orig*scalemorph_ramp)
        if(write_smorph_ramp .and. mod(timehrs,2.0)==0) then                         !If the user wants to see output, write to the file every half hour
          if (scalemorph .gt. 1.0 .and. scalemorph .lt. scalemorph_orig) then
            open (9,file=smorph_file,access='append')
            write(9,1001) scalemorph,timehrs   
            close (9)
          endif
        endif
      endif
    endif
     
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
    end subroutine timing    
    
!*************************************************
    subroutine flow_update
! update field values at different time steps
! by Weiming Wu, NCCHE, Oct. 2008
! modified by Alex Sanchez, USACE-CHL
!*************************************************
#include "CMS_cpp.h"
    use size_def, only: ncellsd
    use geo_def,  only: zb
    use flow_def, only: maxeta, u2, v2, p2, h2, u1, v1, p1, h1, iwet, iwet1,flux, flux1
    use flow_def, only: u, v, p, h, pp, dppx, dppy, eta
    use comvarbl, only: dtime,dtime1,wtsch,ntime,nspinup,nfsch
    use sed_def,  only: sedtrans, zb1, btk, btk1, btk2, ctk, ctk1, ctk2, ibt, singlesize, pbk, pbk1, db, db1
    use sal_def,  only: saltrans, sal, sal1, sal2
    use heat_def, only: heattrans, heat, heat1, heat2
    use prec_def, only: ikind
    
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
    
!Update maximum WSE values   !added MEB 9/20/2021
    do i=1,ncellsD
      maxeta(i) = max(eta(i),maxeta(i))
    enddo
    
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
    end subroutine flow_update
   
 !*****************************************************
    subroutine makeflow1D
 !*****************************************************
    use geo_def,  only: maxcol,maxrow,idmap
    use flow_def, only: u,v,p,uv
    use comvarbl, only: iflow1D
    use prec_def, only: ikind
    
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
    end subroutine

!***************************************************    
    subroutine cleanup
!***************************************************
    implicit none
    
    !call geo_cleanup
    call flow_cleanup
    !call sed_cleanup
    !call sal_cleanup
    
    return
    end subroutine cleanup
    