      subroutine CMS_FLOW_EXP_C()
!-----------------------------------------------------------------------
!                                                                      |
!              CMS Model with Options of Numerical Schemes             |
!           Initial Development by Adele (Militello) Buttolph          |
!             Explicit Optimizations completed by Chris Reed           |
!                                                                      |
!                   Code Integration by Mitchell Brown                 |
!-----------------------------------------------------------------------

!----------------------------------------------------------------------------------------------------
!   from CMS-Flow version 3.75
!   
!   BETA REVISION 
!     Brown   (08/13/2008) - Integrated Sanchez' NET scheme  
!     Reed    (08/20/2008) - Hardbottom/Avalanching changes 
!     Reed    (08/20/2008) - Symmetry corrections (momentum)
!     Brown   (08/28/2008) - Porosity correction and comments
!     Brown   (09/03/2008) - Minor fix for NET without waves
!     Brown   (09/05/2008) - Increased stack size to 20MB   
!   RELEASE REVISION 0
!     Brown   (09/05/2008) - Added check for GRDPATH and error message    
!     Brown   (09/16/2008) - Added AUTO and NORMAL (non-steering) HOTSTART ability
!     Brown   (09/22/2008) - Added code from Alex for different CAPACITY formulas and AVALANCHING for NET
!     Reed    (10/24/2008) - Symmetry fix (INITIALIZE_BC)
!     Sanchez (10/24/2008) - Minor changes in NET parameters - diffshear/diffsm
!     Sanchez (10/24/2008) - Change in Vector global output to Cell Faces (write_global_outputs-mb & open_create_dataset)
!     Brown   (11/11/2008) - Added SMS 10.1 card support for salinity parameters and BCs
!     Sanchez (12/05/2008) - Fix for SUSPLUND routine to WAVES AND CURRENTS (not WAVES OR CURRENTS)
!     Brown   (01/14/2009) - Activated OpenMP in Revision 0, instead of Revision 1 (with Variable D50)
!   RELEASE REVISION 1
!     Brown   (01/16/2009) - Added in Wall Friction (READ_GRID_XMDF)
!     Reed    (01/21/2009) - Added density difference due to salinity gradients into momentum equations
!     Reed    (01/26/2009) - Fix for non-steering hot start capability
!     Brown   (01/26/2009) - Added GUID information to global solution datasets.  Ensures solution is loaded to correct grid in SMS.
!     Brown   (01/17/2009) - Added Variable D50 routines - replaces SEDIMENT_GRAIN_SIZE
!     Brown   (02/05/2009) - Corrected some OpenMP calls and removed the statements for SEDTRANS.
!     Brown   (02/06/2009) - Fix for some ALLOCATE statements
!     Brown   (03/27/2009) - Fix for some H_Multi boundary condition statements
!   RELEASE REVISION 2
!     Reed    (05/27/2009) - Cohesive formulation updates
!     Brown   (05/27/2009) - Fix for write statement timing (tolerance)
!     Brown   (06/05/2009) - 32- and 64-bit version support.
!     Brown   (06/05/2009) - NaN improvement on global write:  -888 written to dataset where NaN exists.
!     Reed    (06/09/2009) - Fix for wave information passed from CMS-Wave - in SUBS-INIT.
!     Brown   (07/22/2009) - Additional output if Hotstart file cannot be written.
!     Brown   (07/22/2009) - Modification to handling of parallel processors/threads (print_status)
!     Sanchez (09/10/2009) - Updated NET subroutine that runs faster and has other improvements
!     Reed    (11/  /2009) - Added Wave bottom friction
!     Brown   (11/11/2009) - Fix for Extracted WL boundary condition.
!     Reed    (11/18/2009) - Fix for AD, Cohesive, and Salinity concentrations - volume needed to be recalculated after hot-start.
!----------------------------------------------------------------------------------------------------
	!PROGRAM CMS_Flow
#include "CMS_cpp.h"
	use EXP_Global_def
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
    
#ifdef DREDGE
    use dredge_def, only: dredging,dredge_time_lapse,dredge_interval
#endif      
#ifdef PROFILE
    use watch_lib
#endif           
            
    implicit none
      
    integer :: i,id,mm,ll,ii,jj
    integer :: nmax,imod,imod6,itime
    real    :: outinterval
    
!!     This section used for OpenMP - added 11/26/2008
!!$    integer omp_get_num_procs
!!$    integer omp_get_max_threads
!!$    integer omp_get_num_threads
!
!!$omp parallel private(numthreads)                                      !NLH 08/07
!!$    NCPU = omp_get_num_procs()
!!$    NTHR = omp_get_max_threads()
!!$    numthreads = omp_get_num_threads()                                !NLH 08/07
!!$omp end parallel                                                      !NLH 08/07

!================================================================================================	
! SET PROJECT PROPERTIES:  NDH 08/08
!     FORTRAN properties, preprocessor, Preprocess Source file: change from NO to YES
!     keep OpenMP Conditional Compilation = YES
!====================================================================================================

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
      call prestart  !set defaults and read card file'
            
      call sim_start_print  !Start timer here   !!Added MEB 4/14/2016

      maxunit = 10                                          
      call ReadCardFile_EXP                                    !read in control file 

      !CALL compute_reftime (IYR,IMO,IDAY,IHR,JDAYS,REFTIME)
      !CALL XF_CALENDAR_TO_JULIAN (0,IYR,IMO,IDAY,IHR,IMIN,ISEC,REFTIME,ERROR)

      Call prestart_EXP      
      
      Call initialize_flow_arrays          !read in grid file and allocate variables       
      CALL PRINT_STATUS                    !print out status information about CMS parameters
      call Initialize_ActivityArray
      CALL INITIALIZE_BC                                         
      CALL INITIALIZE_WABC                 !initialize wind, rad, and wave stresses'
      call initialize_extrapolations_CWR()
      CALL INITIALIZE_TRANSPORT            !initialize salinity and sediment transport       
      call initialize_structures_CWR()
      
      !NEED TO CALL DEALLOCATE IMPLICIT VARIABLES
      !call Zimplicit_var_deallocate()
      
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
          
        if(Nupdate) then 
          !if(flowstats) call flow_stat    !Flow simulation statistics
          !if(sedstats)  call sed_stat     !Sediment simulations statistics
          !if(salstats)  call sal_stat     !Salinity simulation statistics       
        endif
         
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

        call update_momentum()
#ifdef PROFILE
        call watch_stop('momentum')
#endif 

#ifdef PROFILE
        call watch_start('wet_dry')
#endif 
        call update_wetdry()
#ifdef PROFILE
        call watch_stop('wet_dry')
#endif         
        
#ifdef PROFILE
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

#ifdef DREDGE        
        if(dredging) then
          dredge_time_lapse = dredge_time_lapse + dt
          if(dredge_time_lapse .ge. dredge_interval) then
            call dredge_Op
            dredge_time_lapse = 0
          endif
        endif     
#endif
       
        !time = time + dt/3600.d0
        !timesec = timesec + dt
        !time=timesec/3600.d0
        !timeHRS = time

        !rainfall
        if(rainfall) rain_time = rain_time+dt

        ITIME = NTIME+1
        OUTINTERVAL = 300.0                  ! 5 MINUTES
        if (DT .LT. 0.25) OUTINTERVAL = 30.0 ! 0.5 MINUTES
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
        
 7001   FORMAT (' ',I0,' of ',I0,' time steps ',T40,'- ',F0.2,' Hrs')   !11/20/08
 7002   FORMAT (' ',I0,' of ',I0,' time steps')  

      ENDDO
!END OF TIMESTEPPING
      
#ifdef PROFILE
      call watch_output
      call watch_destroy
#endif     

      call sim_end_print     !MEB- 04/14/2016

      RETURN      
      END SUBROUTINE 
 
