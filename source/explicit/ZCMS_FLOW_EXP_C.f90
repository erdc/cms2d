!***********************************************************************
      subroutine CMS_FLOW_EXP_C()
!***********************************************************************
!                                                                      |
!              CMS Model with Options of Numerical Schemes             |
!           Initial Development by Adele (Militello) Buttolph          |
!             Explicit Optimizations completed by Chris Reed           |
!                                                                      |
!                   Code Integration by Mitchell Brown                 |
!***********************************************************************
!   from CMS-Flow version 3.75                                         |
!***********************************************************************
#include "CMS_cpp.h"
    use EXP_Global_def,    only: maxunit, time, dt, ue, ve, linktodummies, sedtransexp, qxn, qyn, dtsed, rainfall, rain_time, outinterval
    USE EXP_transport_def, only: adss,tsed_elapse
    use NupdateMod, only: nupdateint,nupdatecnt,nupdate
    use sed_def,    only: sedtrans
    use flow_def,   only: h,eta,u,v,uv,mturbul
    use size_def,   only: ncells,ncellsd
    use geo_def,    only: zb,cell2cell
    use diag_def,   only: dgunit,dgfile
    use out_def,    only: save_point,obs_cell
    use sal_def,    only: saltrans
    use time_lib,   only: ramp_func
    use comvarbl,   only: tmax,ramp,timehrs,rampdur,ctime,stimet,ntime
    use met_def,    only: windconst,tauwindx,tauwindy,tauwx,tauwy,presconst,windvar,presvar,windsta,pressta
    use cms_def,    only: cmswave
    use dredge_def, only: dredging,dredge_time_lapse,dredge_interval
#ifdef PROFILE
    use watch_lib, only: watch_output, watch_destroy, watch_default, watch_stop, watch_start, watch_init, watch_print
#endif           
            
    implicit none
      
    integer :: i,id,mm,ll,ii,jj
    integer :: nmax,imod,imod6,itime
    !real    :: outinterval
    
#ifdef PROFILE
    call watch_default
    call watch_init
    call watch_print
#endif 
      
    !time step interval for updating focing etc. (speeds up simulation)
    !defined in Module NupdateMod
    NupdateInt = 20
    NupdateCnt = NupdateInt
    Nupdate = .true.
      
#ifdef PROFILE
    call watch_start('Initialize')
#endif       
    call prestart_EXP_shared  !set defaults and read card file'
    call sim_start_print  !Start timer here   !!Added MEB 4/14/2016

    maxunit = 10                                          
    call ReadCardFile_EXP                                    !read in control file 

    !call compute_reftime (IYR,IMO,IDAY,IHR,JDAYS,REFTIME)
    !call XF_CALENDAR_TO_JULIAN (0,IYR,IMO,IDAY,IHR,IMIN,ISEC,REFTIME,ERROR)

    call prestart_EXP      
      
    call initialize_flow_arrays          !read in grid file and allocate variables       
    !call PRINT_STATUS                    !print out status information about CMS parameters
    call Initialize_ActivityArray
    call INITIALIZE_BC                                         
    call INITIALIZE_WABC                 !initialize wind, rad, and wave stresses'
    call initialize_extrapolations_CWR()
    call INITIALIZE_TRANSPORT            !initialize salinity and sediment transport       
    call initialize_structures_CWR()
      
    !NEED TO call DEALLOCATE IMPLICIT VARIABLES
    call Zimplicit_var_deallocate()
      
!write mapping and activity arrays for review
#ifdef DEBUG
    open(unit=404,file='IDMAP.csv')
    write(*,*)'ncellsfull = ',ncellsfull
    do i=1,ncellsfull
      ID = IDMAP(i)
      if(ID.gt.0) then
        write(404,*)ID,i,(cell2cell(mm,ID),mm=1,4), (active(ID,LL),LL=1,3)
      else
        write(404,"(2i6)")ID,i
      endif
    enddo
    close(unit=404)
#endif

    NMAX = (TMAX-Time)*3600/DT + 1
    !allocate(szparms(25,ncellsd))  
      
#ifdef PROFILE
    call watch_stop('Initialize')
#endif      
      
    ramp=ramp_func(timehrs,rampdur) 
    call write_output

    !START OF TIMESTEPPING
    !DO N=1,NMAX
    do while (ctime<stimet) !Note that stimet includes ramp period
      call exp_timing
        
      if(cmswave) call wave_eval !Run wave model, and interpolate in space and time
#ifdef DEV_MODE       
      if(constant_waves) call wave_const_update !Set constant wave parameters, for idealized cases and lab experiments
      if(mbedfric == 0) call fric_rough_eval
#endif         
        
      Nupdate = .false.
      if(NupdateCnt .eq. NupdateInt) then
        Nupdate = .true.
        NupdateCnt = 0
      endif
      NupdateCnt = NupdateCnt + 1
           
      !set the ramp coefficient for this time step if maximum of 1.0 not already met.
      if (ramp.lt.1.0) ramp = ramp_func(timehrs,rampdur)
      call update_wse_bc()
      call update_Q_bc()
      call update_uv_from_qxqy()
           
      !Radiation stress and wave characteristics 
      !ctime=timesecs
      h=-zb+eta  !set h for wave interpoaltion and other functions   
        
      !if(cmswave) call waveinterpol         
        
      !update extrapolations
      call update_q_and_u_extrapolations()
         
#ifdef PROFILE
      call watch_start('CC_vel')
#endif            
         
      !get cell centered velocities for output       
!$OMP PARALLEL DO         
      do i=1,ncells
        u(i) = (uE(i)+uE(cell2cell(2,i)))/2.
        v(i) = (vE(i)+vE(cell2cell(1,i)))/2.
        uv(i) = sqrt(u(i)*u(i)+v(i)*v(i))
      enddo
!$OMP END PARALLEL DO

      ii=0
      do i=ncells+1,ncellsD
        ii=ii+1
        jj = linktodummies(ii)
        u(i) = u(jj)
        v(i) = v(jj)
        uv(i) = uv(jj)          
      enddo   
         
#ifdef PROFILE
      call watch_stop('CC_vel')
#endif            
         
      if(obs_cell)   call write_obs_cell !Observation cells
      if(save_point) call write_save_point !Save point cells, Mitch 5/8/2012
          
      !if(Nupdate) then 
      !  if(flowstats) call flow_stat    !Flow simulation statistics
      !  if(sedstats)  call sed_stat     !Sediment simulations statistics
      !  if(salstats)  call sal_stat     !Salinity simulation statistics       
      !endif
         
      !Wave Wetting and Drying
      if(cmswave .and. Nupdate) call wave_wetdry       
         
      !Wind and Atmospheric Pressure        
      if(Nupdate) then
        if(windconst) then
          call windcurve_eval    !Spatially constant wind field
          tauwindx = tauwx
          tauwindy = tauwy
        endif
        if(presconst) call prescurve_eval    !Spatially constant atmospheric pressure    
        if(windvar .or. presvar) call windpresfield_eval !Spatially variable wind and pressure fields   
        if(windsta .or. pressta) call metsta_eval   !Spatially variable wind and pressure stations
      endif
        
#ifdef PROFILE
      call watch_start('diff_fric')
#endif          
      if(Nupdate)then
        call fric_eval              
        if(mturbul .ge. 3) call derivativeEXP
        if(mturbul .eq. 4) call diswallEXP
        call flow_eddyvis         
      endif
        
#ifdef PROFILE
      call watch_stop('diff_fric')
      call watch_start('momentum')
#endif

      call update_momentum()
#ifdef PROFILE
      call watch_stop('momentum')
      call watch_start('wet_dry')
#endif 
      call update_wetdry()
#ifdef PROFILE
      call watch_stop('wet_dry')
      call watch_start('wse')
#endif
      call update_wse()
#ifdef PROFILE
      call watch_stop('wse')
#endif             
      call update_WABC()
      call update_sedtrans()  !this is explicit AD/TL scheme(s)
      if(sedtrans .and. .not. sedtransEXP) then  !this is the explicit NET scheme
        do i = 1,ncells
          ADSS(i)%qx = ADSS(i)%qx + qxn(i)*dt
          ADSS(i)%qy = ADSS(i)%qy + qyn(i)*dt
        enddo 
        tsed_elapse = tsed_elapse + dt
        if(tsed_elapse .ge. dtsed ) then   
          call sed_exp 
          do i = 1,ncells            
            ADSS(i)%qx = 0
            ADSS(i)%qy = 0     
          enddo
          tsed_elapse = 0.0
        endif
      endif
        
      if(saltrans) call update_salinity()
      IF(rainfall) call rain_evap()    
        
#ifdef PROFILE
      call watch_start('update_vars')
#endif 
      call update_variables()
#ifdef PROFILE
      call watch_stop('update_vars')
#endif 
      call update_culverts()
      call write_output       

      if(dredging) then
        dredge_time_lapse = dredge_time_lapse + dt
        if(dredge_time_lapse .ge. dredge_interval) then
          call dredge_Op
          dredge_time_lapse = 0
        endif
      endif     
       
      !time = time + dt/3600.d0
      !timesec = timesec + dt
      !time=timesec/3600.d0
      !timeHRS = time

      !rainfall
      if(rainfall) rain_time = rain_time+dt

      ITIME = NTIME+1
      
      if (outinterval .lt. 0.0) then           ! If not already set, set based on the DT
        if (dt .ge. 0.25) outinterval = 300.0  ! 5 minutes
        if (dt .lt. 0.25) outinterval =  30.0  ! 0.5 minutes
      endif

      IMOD  = NINT(OUTINTERVAL/DT)
      IMOD6 = NINT(OUTINTERVAL*6/DT)

      open(dgunit,file=dgfile,access='append') 
      
      IF ((IMOD.NE.0).AND.(IMOD6.NE.0)) THEN
        if((ITIME/IMOD)*IMOD.eq.ITIME)THEN                      !PRINT TIMESTEP COUNT EVERY DESIGNATED INTERVAL DEPENDENT ON TIMESTEP
          WRITE (*,7002)      ITIME,nmax-1
          WRITE (DGUNIT,7002) ITIME,nmax-1
        ENDIF
      ENDIF
        
      close(dgunit)
        
7001  FORMAT (' ',I0,' of ',I0,' time steps ',T40,'- ',F0.2,' Hrs')   !11/20/08
7002  FORMAT (' ',I0,' of ',I0,' time steps')  

    enddo
!END OF TIMESTEPPING
      
#ifdef PROFILE
    call watch_output
    call watch_destroy
#endif     

    call sim_end_print     !MEB- 04/14/2016

    return      
    end subroutine
