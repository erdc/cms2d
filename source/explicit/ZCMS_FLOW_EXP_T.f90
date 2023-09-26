!***********************************************************************
    subroutine CMS_FLOW_EXP_T()
!-----------------------------------------------------------------------
!                                                                      |
!              CMS Model with Options of Numerical Schemes             |
!           Initial Development by Adele (Militello) Buttolph          |
!             Explicit Optimizations completed by Chris Reed           |
!                                                                      |
!                   Code Integration by Mitchell Brown                 |
!-----------------------------------------------------------------------
!   from CMS-Flow version 3.75
!   Explicit Solver integrated with Implicit Code
!   and modified for using Telescoping Grid
!-----------------------------------------------------------------------
#include "CMS_cpp.h"
    use EXP_Global_def,    only: maxunit, time, dt, num_linktodummies, linktodummiestel, dtsed, rainfall, rain_time, outinterval
    USE EXP_transport_def, only: adss,tsed_elapse
    use EXP_TELESCOPING,   only: numregcells, numtbcells, regcells, tbcells, cellfaces, xface_qn, yface_qn, xface_length, yface_length
    use EXP_TELESCOPING,   only: xface_advdif_i, xface_advdif_c, yface_advdif_i, yface_advdif_c, numxfaces, numyfaces, xsedtransq, ysedtransq
    use NupdateMod, only: nupdateint,nupdatecnt,nupdate         
    use sed_def,    only: sedtrans
    use flow_def,   only: h,eta,u,v,uv,mturbul
    use size_def,   only: ncells, ncellsfull
    use geo_def,    only: zb
    use diag_def,   only: dgunit,dgfile
    use out_def,    only: save_point,obs_cell
    use sal_def,    only: saltrans
    use time_lib,   only: ramp_func
    use comvarbl,   only: tmax,ramp,timehrs,rampdur,ctime,stimet,ntime,dtime
    use met_def,    only: windconst,tauwindx,tauwindy,tauwx,tauwy,presconst,windvar,presvar,windsta,pressta
    use cms_def,    only: cmswave
    use dredge_def, only: dredging,dredge_time_lapse,dredge_interval
    use BalanceCheck_def    
#ifdef PROFILE
    use watch_lib
#endif
    implicit none
      
    integer :: i,id,mm,ll,ii
    integer :: nmax,imod,imod6,itime
    integer :: id1,id2,j
    real    :: totL,DepthT        
    !real    :: outinterval
    character*80 filename
    
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

    !CALL compute_reftime (IYR,IMO,IDAY,IHR,JDAYS,REFTIME)
    !CALL XF_CALENDAR_TO_JULIAN (0,IYR,IMO,IDAY,IHR,IMIN,ISEC,REFTIME,ERROR)

    Call prestart_EXP      
      
    Call initialize_flow_arrays_Tel          !read in grid file and allocate variables       
    !CALL PRINT_STATUS                        !print out status information about CMS parameters
    call Initialize_ActivityArray_Tel
    CALL INITIALIZE_BC_Tel                   !ALL Q AND SAL  - WILL NEED WORK FOR TELE                                    
    CALL INITIALIZE_WABC_Tel                 !initialize wind, rad, and wave stresses'
    call initialize_extrapolations_CWR_Tel()
    CALL INITIALIZE_TRANSPORT_tel            !WILL NEED WORK FOR TELE   !initialize salinity and sediment transport       
    call initialize_structures_CWR()
      
    !NEED TO CALL DEALLOCATE IMPLICIT VARIABLES
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

    if(balcheck) call Init_BalanceCheck_tel
    !write(*,*)'TEL REG = ',numREGXfaces,numREGYfaces    
    !write(*,*)'TEL     = ',numTBXfaces,numTBYfaces          
      
    !START OF TIMESTEPPING
    !DO N=1,NMAX
    do while (ctime<stimet) !Note that stimet includes ramp period
      call exp_timing

      !Added by Chris Reed - submitted 10/14/16        
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
      call update_Q_bc_tel()  
      call update_uv_from_qxqy_tel()
           
      !Radiation stress and wave characteristics 
      !ctime=timesecs
      h=-zb+eta  !set h for wave interpoaltion and other functions        
  
      !if(cmswave) call waveinterpol         
         
      !update extrapolations
      call update_q_and_u_extrapolations_tel()
        
#ifdef PROFILE
      call watch_start('CC_vel')
#endif         
         
!get cell centered velocities for output 
!$OMP PARALLEL DO private(i,DepthT)       
      do ii=1,numREGCells
        i=REGCells(ii)          
        depthT = eta(i) - zb(i)
        if(depthT .eq. 0) then
          depthT = max(depthT,0.0000001)
        endif
        U(i) = (xface_QN(cellfaces(3,i)) +xface_QN(cellfaces(7,i)))/2.0
        u(i) = u(i)/DepthT
        v(i) = (yface_QN(cellfaces(1,i)) +yface_QN(cellfaces(5,i)))/2.0
        v(i) = v(i)/DepthT
        uv(i) = sqrt(u(i)*u(i)+v(i)*v(i))
      enddo
!$OMP END PARALLEL DO
      
#ifdef PROFILE
      call watch_stop('CC_vel')
#endif       
  
!$omp parallel do private (i,totL,DepthT)     
      do ii=1,numTBCells
        i=TBCells(ii)          
        depthT = eta(i)-zb(i)
        U(i) =  xface_QN(cellfaces(3,i))*xface_length(cellfaces(3,i)) &
               +xface_QN(cellfaces(4,i))*xface_length(cellfaces(4,i)) &
               +xface_QN(cellfaces(7,i))*xface_length(cellfaces(7,i)) &
               +xface_QN(cellfaces(8,i))*xface_length(cellfaces(8,i)) 
        TotL =  xface_length(cellfaces(3,i)) &
               +xface_length(cellfaces(4,i)) &
               +xface_length(cellfaces(7,i)) &
               +xface_length(cellfaces(8,i))             
        u(i) = u(i)/(TotL*DepthT)
        v(i) =  yface_QN(cellfaces(1,i))*yface_length(cellfaces(1,i)) &
               +yface_QN(cellfaces(2,i))*yface_length(cellfaces(2,i)) &
               +yface_QN(cellfaces(5,i))*yface_length(cellfaces(5,i)) &
               +yface_QN(cellfaces(6,i))*yface_length(cellfaces(6,i)) 
        TotL =  yface_length(cellfaces(1,i)) &
               +yface_length(cellfaces(2,i)) &
               +yface_length(cellfaces(5,i)) &
               +yface_length(cellfaces(6,i))             
        v(i) = v(i)/(TotL*DepthT)
        uv(i) = sqrt(u(i)*u(i)+v(i)*v(i))
      enddo         
!$omp end parallel do

      do i=1,num_linktodummies
        id1 = linktodummiesTel(1,i) 
        id2 = linktodummiesTel(2,i)              
        u(id2) = u(id1)
        v(id2) = v(id1)
        uv(id2) = uv(id1)
      enddo
         
      if(obs_cell)   call write_obs_cell !Observation cells
      if(save_point) call write_save_point !Save point cells, Mitch 5/8/2012
                  
      !if(Nupdate) then 
      !  if(flowstats) call flow_stat   !Flow simulation statistics
      !  if(sedstats)  call sed_stat    !Sediment simulations statistics
      !  if(salstats)  call sal_stat    !Salinity simulation statistics       
      !endif
         
      !Added by Chris Reed - submitted 10/14/16        
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
#endif         

#ifdef PROFILE
      call watch_start('momentum')
#endif 
      !set advective/diffusive fluxes to zero     
      xface_advdif_I = 0.0    
      xface_advdif_C = 0.0
      yface_advdif_I = 0.0     
      yface_advdif_C = 0.0  
       
      !xface_q = 0.0
      !xface_vel = 0.0
      !yface_q = 10.0
      !yface_vel = 5.0 
      !eta = 0.0
      !do i=1,numYfaces
      !  if(yface_wall(i)) then
      !    yface_q(i) = 0.0
      !    yface_vel(i) = 0.0
      !  endif
      !enddo
      
      !Added by Chris Reed - submitted 10/14/16        
      call update_wetdry_tel_pre()
      call update_advection_tel()      
      call update_advection_tel_reg() 
       
      !yface_q = 0.0
      !yface_vel = 0.0      
       
      call update_momentum_tel()
      call update_momentum_tel_reg   

#ifdef PROFILE
      call watch_stop('momentum')
#endif 

#ifdef PROFILE
      call watch_start('wet_dry')
#endif 
      call update_wetdry_tel()
#ifdef PROFILE
      call watch_stop('wet_dry')
#endif         


#ifdef PROFILE
      call watch_start('wse')
#endif 

      call update_wse_tel()

#ifdef PROFILE
      call watch_stop('wse')
#endif         
      call update_WABC()
   
      !Modified by Chris Reed - submitted 10/14/16        
      !call update_sedtrans()  !this is explicit AD/TL scheme(s)
      if(sedtrans) then ! .and. .not. sedtransEXP) then  !this is the explicit NET scheme
        !do i = 1,ncells
        !  ADSS(i)%qx = ADSS(i)%qx + qxn(i)*dt
        !  ADSS(i)%qy = ADSS(i)%qy + qyn(i)*dt
        !enddo 
!$OMP PARALLEL DO 
        do i = 1,numxfaces
          xSedTransQ(i) = xSedTransQ(i) + xface_qn(i)*dt
        enddo 
!$OMP END PARALLEL DO
!$OMP PARALLEL DO 
        do i = 1,numyfaces
          ySedTransQ(i) = ySedTransQ(i) + yface_qn(i)*dt
        enddo
!$OMP END PARALLEL DO                   
        tsed_elapse = tsed_elapse + dt
        if(tsed_elapse .ge. dtsed ) then   
          call sed_exp_tel
          !do i = 1,ncells            
          !  ADSS(i)%qx = 0
          !  ADSS(i)%qy = 0     
          !enddo
!$OMP PARALLEL DO 
          do i = 1,numxfaces
            xSedTransQ(i) = 0
          enddo
!$OMP END PARALLEL DO
!$OMP PARALLEL DO 
          do i = 1,numyfaces
            ySedTransQ(i) = 0
          enddo
!$OMP END PARALLEL DO            
          tsed_elapse = 0.0
        endif
      endif
        
      if(saltrans) call update_salinity_tel()
      IF(rainfall) call rain_evap()  
        
#ifdef PROFILE
      call watch_start('update_vars')
#endif     
      call update_variables_tel()
         
#ifdef PROFILE
      call watch_stop('update_vars')
#endif 
      call update_culverts()
      call write_output             !This is a second "write_output" during the same iteration.  MEB- 041916

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
        
      !run mass balance checks
      if(balcheck) call balancecheck_tel()
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

#ifdef DEBUG
    !DEBUG-CHECK OUTUT
    open(unit=3050,file='xfaces2.csv')
    do i=1,numxfaces
      write(3050,"(3(i20,','),2(e15.7,','),e15.7)")i, &
           (xface_cells(j,i),j=1,2),xface_q(i),xface_length(i),xface_vel(i)
    enddo 
    close(3050)    
    open(unit=3050,file='yfaces2.csv')
    do i=1,numyfaces     
      write(3050,"(3(i20,','),2(e15.7,','),e15.7)")i, &
           (yface_cells(j,i),j=1,2),yface_q(i),yface_length(i),yface_vel(i)
    enddo 
    close(3050)      
      
    open(unit=9940,file='xfacestuff.csv')
    do i=1,numXfaces
      write(9940,"(i6,3(f15.9))")i,xface_q(i),xface_qn(i),xface_vel(i)
    enddo
#endif      

    return
    end subroutine