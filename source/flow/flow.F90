!===============================================================================
! Flow (hydrodynamic) routines
! 
! Contains the following routines:
!   flow_default - Sets the default parameters for flow variables
!   flow_cards   - Reads the flow cards and overrides the default settings
!   flow_init    - Initializes the flow variables
!   flow_print   - Prints the flow parameters to the diagnostic file and screen
!   flow_alloc - Allocation of flow variables
!   flow_cleanup  - Deallcoates flow variables
!   flow_step_stat - Calculates flow time step statistics and prints them to 
!                    the screen and diagnostic file
!   flow_wetdry   - Wetting and drying algorithm
!   flow_pond     - Judges imponding and dead water
!   flow_convflux      - Initialising volume fluxes based on cell-center velocities
!
! written by Weiming Wu, NCCHE and Alex Sanchez, USACE-CHL
!===============================================================================

!***************************************************************************   
    subroutine flow_default
! Sets the default parameters for flow variables
! written by Alex Sanchez, USACE-CHL; Weiming Wu, NCCHE
!***************************************************************************
    use flow_def
    use comvarbl 
    use prec_def   
    implicit none

    !Gravity
    grav = 9.81_ikind  !Gravitational constant, m/s^2
    !grav = 9.80665_ikind !Standard Gravity [m/s^2] (established by the 3rd General Conference on Weights and Measures in 1901)
    gravinv = 1.0_ikind/grav
    sqrtgrav = sqrt(grav)
    
    !Water properties
    idensit = 2 !0-user specified, 1-Ekart (1958), 2-UNESCO (1981), 3-Fofonoff (1985)
    iviscos = 2 !0-user specified, 1-emp, 2-temp and salinity
    isalt = 0   !0-none(default), 1-user specified, 2-initial condition
    iheat = 0   !0-none(default), 1-user specified, 2-initial condition
    watertemp = 15.0  !Water temperature [ºC]
    watersalt = 35.0  !Water salinity [ppt]
    rhow = 1025.0   !Water density [kg/m^3] (approximate for 15ºC and 34 ppt)    
    viscos = 1.1812e-6 !Kinematic water viscosity [m^2/s] (approx for 15ºC and 35 ppt)    
    
    !Wetting and drying
    hmin = 0.05           !Minimum water depth [m]
    hdry = 0.05        !Wetting and drying criteria [m]
    hdry1 = 5.0*hdry   !Wetting and drying criteria for one-cell-wide channels [m]
    hdry2 = 2.0*hdry   !Wetting and drying criteria for two-cell-wide channels [m]
    ponding = .false.  !Allow water ponding. Otherwise isolated bodies of water are removed. Default off for stability
    narrowchannels = .false. !Allows narrow channels in wetting and drying. Default off for stability   
    
    !Coriolis
    icoriolisplane = 1 !1=f-plane,2=beta-plane
    fcoriolis = 0.0

    !Solution scheme
    nfsch = 0 !0-implicit, 1-explicit
    
    !Explicit flow solver settings    
    norder = 1           !Order of explicit solution
    nriem = 2            !Riemann solver

    !Time settings    
    stime = 0.0 
    etime = 21600.0 ![sec]
    nramp = 0
    rampdur = 0.0      
    ramp = 0.0
    nspinup = 0
    nprt = 20
    
    !Turblence model
    crossdiff = .false.
    mturbul = 5
    cvismax = 20.0
    cviscon = 1.0e-6 !constant (base value)
    cvisbot = 0.0667 !bottom current shear
    if(mturbul==5) cvisbot = 0.0667 !bottom current shear
    cvishor = 0.4    !horizontal current, OLD
    if(mturbul==5) cvishor = 0.2    !horizontal curren, NEW 
    if(mturbul==4) cvishor = 0.2    !horizontal current 
    cviswav = 0.5    !gamma in Larson and Kraus (1991) eddy viscosity for waves, OLD
    if(mturbul==5) cviswav = 0.5    !gamma in Larson and Kraus (1991) eddy viscosity for waves, NEW
    cviswavbrk = 2.5 !viscosity due to wave radiation stresses *IMPORTANT*
    if(mturbul==5) cviswavbrk = 0.1 !Viscosity due to wave breaking dissipation *IMPORTANT*

    !Wave mass flux
    waveflux = .false.

    !Numerics Methods
    nthr = 1  !Threads     
    
    !Skewness correction
    skewcor = .true. !0-Off, 1-On
    !facblend = 1.0  !Placeholder. Not implemented yet
    
    !Implicit scheme
    ntsch = 1   !1-First Order, 2-Second Order
    wtsch = -1.0 !0-First order, 1-Second order
    dtvar = .true. !Variable tim-stepping
    pred_corr = .false. !Predictor-corrector
    
    !Advection scheme
    ndsch = 4 !Exponential      
    !ndsch = 3  ! in Chris' code
    
    !Screen Output variable
    nprint = 10           
    
    !Convergence variables
    maxit = 0     !Maximum number of iterations for solver
    maxit0 = 0    !Outer-loop iterations
    nswp0 = 0     !Inner-loop iterations
    rmom = 0.0
    
    !Implicit Scheme Relaxation coefficients for pp, u , v , zs , te, ed, s , p ,den,vis
    relax = 0.6            !Wu
    relax(1) = 1.0       !pp     !This has to be 1.0   !Wu
    relax(2) = 0.8       !u
    relax(3) = relax(2)  !v
    relax(4) = 0.8       !c
    relax(5) = 0.8       !s
    relax(8) = 1.0       !p
    relaxsor = 1.8
    facpp = 0.2
    
    !Implicit Scheme Convergence thresholds
    rmommaxp     = 1.0e-3 !Divergence
    rmomtargetp  = 1.0e-4 !Reduce time step
    !rmominc4p    = 1.0e-5 !Increase time step after four iterations
    !rmominc2p    = 1.0e-6 !Increase time step after two iterations    
    rmomminp     = 1.0e-8 !Converged
    rmomratiop   = 0.5
    rmomabschgp  = 1.0e-8
    rmomrelchgp  = 0.05
    rmommaxuv    = 1.0e-2 !Divergence
    rmomtargetuv = 1.0e-3 !Reduce time step
    !rmominc4uv   = 1.0e-4 !Increase time step after four iterations
    !rmominc2uv   = 1.0e-5 !Increase time step after two iterations    
    rmomminuv    = 1.0e-7 !Converged
    rmomratiouv  = 1.0
    rmomabschguv = 1.0e-7
    rmomrelchguv = 0.01
    
    !Make flow 1D    
    iFlow1D = 0 
    
    !Volume conservation correction based on final inter-cell fluxes
    volcor = .false.
    flowvolbal = .false.
    
    !Limits
    presmax = 20.0*grav !m*(m/s^2)
    velmax = 10.0       !m/s
    
    !Advanced card default file name
    advfile = 'advanced.cmcards'
    
    return        
    endsubroutine flow_default
        
!*************************************************************   
    subroutine flow_cards(cardname,foundcard,doPrint)
! Reads wind data from Model Parameters file
! written by Alex Sanchez, USACE-CHL
!*************************************************************
    use diag_def
    use geo_def
    use comvarbl
    use diag_lib
    use flow_def
    use rol_def,  only: rolflux
    use flow_lib, only: water_viscosity_kinematic
    use solv_def, only: asolv
    use time_lib, only: julday,calendar2julian,julianday2calendarmonthday
    use heat_def, only: heatic
    use cms_def,  only: aValue
    use EXP_Global_def, only: outinterval
    implicit none
    
    character(len=*),intent(inout) :: cardname
    logical,intent(out)            :: foundcard
    logical,intent(in)             :: doPrint
    
    integer :: i,ierr,omp_get_max_threads
    character(len=37) :: cdum
    character(len=12) :: stringname
    character(len=200) :: apath,aname
    character(len=10) :: aext
    
    foundcard = .true.
    selectcase(cardname)
    case('PARAMS_FILE')
      backspace(77)
      read(77,*) cardname, mpfile
      call fileparts(mpfile,apath,aname,aext)  
      if(len_trim(apath)==0)then
        mpfile = trim(flowpath) // mpfile
      endif
      
    case('ADVANCED_FILE')
      backspace(77)
      read(77,*) cardname, advfile
      call fileparts(advfile,apath,aname,aext)
      if(len_trim(apath)==0)then
        advfile = trim(flowpath) // advfile
      endif

    case('2D_MODE')
     !DO NOTHING
     
    case('GRAVITY','GRAVITATIONAL_CONSTANT') 
      backspace(77)
      read(77,*) cardname, grav
      
    case('WATER_DENSITY')
      call card_scalar(77,'kg/m^3','kg/m^3',rhow,ierr)
      idensit = 0 !User specified
          
    case('WATER_TEMPERATURE')
      call card_scalar(77,'C','C',watertemp,ierr)
      heatic = watertemp
      
    case('WATER_VISCOSITY_KINEMATIC','WATER_KINMATIC_VISCOSITY')  
      backspace(77)
      read(77,*) cardname,viscos
      iviscos = 0 !User specified
      
    !=== Wetting and Drying ======================
    case('MIN_DEPTH','MINIMUM_DEPTH','MINIMUM_WATER_DEPTH')
      call card_scalar(77,'m','m',hmin,ierr)
      
    case('DRYING_DEPTH')
      call card_scalar(77,'m','m',hdry,ierr)    
    
    case('DRYING_DEPTH_ONE_CELL','ONE_CELL_DRYING_DEPTH')
      call card_scalar(77,'m','m',hdry1,ierr)
      
    case('DRYING_DEPTH_TWO_CELL','TWO_CELL_DRYING_DEPTH')
      call card_scalar(77,'m','m',hdry2,ierr)  
      
    case('WATER_PONDING')
      call card_boolean(77,ponding,ierr)
        
    case('ONE_CELL_WIDE_CHANNELS','ONE-CELL-WIDE_CHANNELS','NARROW_CHANNELS')
      call card_boolean(77,narrowchannels,ierr)
          
    !=== Coriolis =========================
    case('CORIOLIS_APPROXIMATION','CORIOLIS_PLANE')
      backspace(77)
      read(77,*) cardname, cdum
      call lowercase(cdum)
      if(cdum(1:1)=='b')then
        icoriolisplane = 2 !BETA
      else
        icoriolisplane = 1 !F
      endif
          
    !=== Timing ==========================
    case('HYDRO_TIMESTEP')
      call card_scalar(77,'sec','sec',dtime,ierr)
      deltime = dble(dtime)
      dtimebeg = dtime
      dtime1 = dtime
      if(dtime<1.0e-6 .and. doPrint)then
        write(msg2,*) dtime
        call diag_print_error('Invalid time step: ',msg2)
      endif  
      
    case('VARIABLE_TIMESTEP','ADAPTATIVE_TIMESTEP')  
      call card_boolean(77,dtvar,ierr)  
      
    case('STARTING_DATE_TIME')
      call card_datetime(77,iyr,imo,iday,ihr,imin,isec) !YYYY-MM-DD HH:MM:SS UTC
      call calendar2julian(iyr,imo,iday,ihr,imin,isec,tjulday0)
      tjulhr0 = tjulday0*24.0
      tjulhryr = julday(iyr,1,1)*24.0    !start of year, hours
      jday = (tjulhr0-tjulhryr)/24    !Julian day (1-366)      
      dtj = tjulhr0-tjulhryr
          
    case('STARTING_JDATE')
      backspace(77)
      read(77,*) cardname, stringname
      read(stringname(1:2),*) iyr
      if(iyr<100)then
        if(iyr>20)then
          iyr = iyr + 1900
        else
          iyr = iyr + 2000
        endif
      endif
      read(stringname(3:5),*) jday !julian day (1-366)
      call julianday2calendarmonthday(iyr,jday,imo,iday)
    
    case('STARTING_JDATE_HOUR')
      backspace(77)
      read(77,*) cardname, ihr      
      imin = 0; isec = 0
      call calendar2julian(iyr,imo,iday,ihr,imin,isec,tjulday0)
      tjulhr0 = tjulday0*24.0
      tjulhryr = julday(iyr,1,1)*24.0 !start of year, hours
      dtj = tjulhr0-tjulhryr
          
    case('DURATION_RUN','SIMULATION_DURATION')
      call card_scalar(77,'hrs','hrs',tmax,ierr)
      stimet = tmax*3600.0 !seconds
          
    case('DURATION_RAMP','RAMP_DURATION')
      call card_scalar(77,'days','hrs',rampdur,ierr)
          
    case('SPIN_UP_ITERATIONS','SPINUP_ITERATIONS')  
      backspace(77)
      read(77,*) cardname,nspinup !iterations                
      nspinup = max(nspinup,0)
      
    !=== Predictor-Corrector =================
    case('PREDICTOR_CORRECTOR_SCHEME','PREDICTOR_CORRECTOR')
      call card_boolean(77,pred_corr,ierr)
          
    !=== Advection ================  
    case('USE_ADVECTION_TERMS') !Outdated
      backspace(77)
      read(77,*) cardname, cdum
      if(cdum(1:3)=='OFF')then
        ndsch = 0
        if(doPrint) call diag_print_warning('Turning off advection',&
          '  Only recomended for testing or idealized cases')
      endif
          
    !=== Mixing ======================
    case('USE_MIXING_TERMS') !Outdated
      backspace(77)
      read(77,*) cardname, cdum
      if(cdum(1:3)=='OFF')then
        mturbul = 0   
        cviscon = 0.0
        if(doPrint) call diag_print_warning('Turning off diffusion',&
          '  Only recomended for testing or idealized cases')
      endif
        
    case('TURBULENCE_MODEL')
      backspace(77)
      read(77,*) cardname, cdum
      selectcase(cdum)
      case('CONSTANT')
        mturbul = 0
      case('FALCONER')
        mturbul = 1
      case('PARABOLIC')
        mturbul = 2
      case('SUBGRID-WU','SUBGRID_WU')
        mturbul = 3
      case('MIXING-LENGTH','MIXING_LENGTH','MIXING')
        mturbul = 4
      case('SUBGRID','SUBGRID-V2')
        mturbul = 5
      case default
        if(doPrint) call diag_print_warning('Invalid turbulence model',&
          '  Check CMS-Flow Card File','  Using Subgrid model')
        mturbul = 5
      endselect
         
    case('EDDY_VISCOSITY_CONSTANT')
      backspace(77)
      read(77,*) cardname, cviscon  
!      cviscon = max(cviscon,1.0e-15) !used to avoid divide by zero
            
    case('EDDY_VISCOSITY_BOTTOM')
      backspace(77)
      read(77,*) cardname, cvisbot            
         
    case('EDDY_VISCOSITY_HORIZONTAL','EDDY_VISCOSITY_HORIZ')
      backspace(77)
      read(77,*) cardname, cvishor    
        
    case('EDDY_VISCOSITY_WAVE')
      backspace(77)
      read(77,*) cardname, cviswav                                             
          
    case('EDDY_VISCOSITY_BREAKING','VISCOSITY_COEF_WAVE_BREAKING')
      backspace(77)
      read(77,*) cardname, cviswavbrk   
          
    case('EDDY_VISCOSITY_MAX_LIMIT','EDDY_VISCOSITY_MAX','EDDY_VISCOSITY_MAXIMUM')
      backspace(77)
      read(77,*) cardname, cvismax
      
    case('EDDY_VISCOSITY_CROSS_DIFFUSION','CROSS_DIFFUSION')  
      call card_boolean(77,crossdiff,ierr)
        
    !==== Wave Flux ======================================================
    case('WAVE_FLUX','WAVE_MASS_FLUX','STOKES_VELOCITY',&
         'STOKES_VELOCITIES','WAVE_VOLUME_FLUX','WAVE_TRANSPORT_VELOCITY')
      call card_boolean(77,waveflux,ierr)
        
    case('ROLLER_MASS_FLUX','ROLLER_VOLUME_FLUX','ROLLER_FLUX_VELOCITY')
      call card_boolean(77,rolflux,ierr)    
      
    case('EXPLICIT_PRINT_INTERVAL')
      call card_scalar(77,'sec','sec',outinterval,ierr)
      continue
          
    !==== Numerical Methods ====================================================
    case('SOLUTION_SCHEME')
      backspace(77)
      read(77,*) cardname, cdum 
      if(cdum(1:3)=='EXP')then
        nfsch = 1 !Explicit
        !call diag_print_warning('Invalid temporal solution scheme',&
        !  '  The inline model executable only contains the implicit temporal scheme.',&
        !  '  To run the explicit temporal scheme, use the explicit model executable.')            
      elseif(cdum(1:3)=='IMP')then
        nfsch = 0 !Implicit
      elseif(cdum(1:3)=='SEM')then  
        nfsch = 2 !Semi-implicits  
      endif
          
    case('RIEMANN_SOLVER')  
      backspace(77)
      read(77,*) cardname, cdum 
      do i=1,size(ariem)
        if(cdum==ariem(i))then
          nriem = i
          exit
        endif
      enddo          
        
    case('EXPLICIT_ORDER','EXPLICIT_SOLUTION_ORDER')  
      backspace(77)
      read(77,*) cardname, norder      
      norder = max(min(norder,2),1)            
          
    case('HYDRO_MAX_ITERATIONS','HYDRO_MAX_ITER','SOLVER_MAX_ITER')
      backspace(77)
      read(77,*) cardname, maxit
      maxit0=maxit          
        
    case('MATRIX_SOLVER','SOLVER_TYPE')
      backspace(77)
      read(77,*) cardname, cdum
      do i=0,size(asolv)-1
        if(cdum==asolv(i))then
          nsolv = i
          exit
        endif
      enddo         
    
    case('PRESSURE_CORRECTION_FACTOR','PRES_CORR_FAC','PRES_COR_FAC')
      backspace(77) 
      read(77,*) cardname, facpp  
      
    case('PRESSURE_RELAXATION','PRESSURE_RELAX','PRES_RELAX')
      backspace(77) 
      read(77,*) cardname, relax(8)
      
    case('VELOCITY_RELAXATION','VELOCITY_RELAX','VEL_RELAX')
      backspace(77) 
      read(77,*) cardname, relax(2)
      relax(3)=relax(2)    
          
    case('SUCCESSIVE_OVER_RELAXATION_CONSTANT','SOR_CONSTANT') 
      backspace(77) 
      read(77,*) cardname, relaxsor          
        
    case('PRESSURE_ITERATIONS','PRESSURE_ITER','PRES_ITER','PRESS_ITER') 
      backspace(77) 
      read(77,*) cardname, nswp(1)
      nswp0(1)=nswp(1)
          
    case('VELOCITY_ITERATIONS','VELOCITY_ITER','VEL_ITER') 
      backspace(77) 
      read(77,*) cardname, nswp(2)
      nswp(3)=nswp(2); nswp0(2:3)=nswp(2:3) 
          
    case('PRESSURE_MAX_RESIDUAL','PRES_MAX_RES')
      backspace(77) 
      read(77,*) cardname,rmommaxp
    
    case('PRESSURE_TARGET_RESIDUAL','PRES_TARGET_RES')
      backspace(77) 
      read(77,*) cardname,rmomtargetp
      
    case('PRESSURE_MIN_RESIDUAL','PRES_MIN_RES')
      backspace(77) 
      read(77,*) cardname,rmomminp
        
    case('VELOCITY_MAX_RESIDUAL','VEL_MAX_RES')
      backspace(77) 
      read(77,*) cardname,rmommaxuv
          
    case('VELOCITY_TARGET_RESIDUAL','VEL_TARGET_RES')
      backspace(77) 
      read(77,*) cardname,rmomtargetuv
          
    case('VELOCITY_MIN_RESIDUAL','VEL_MIN_RES')
      backspace(77) 
      read(77,*) cardname,rmomminuv
      
    case('VELOCITY_MAX','VELOCITY_MAXIMUM')
      call card_scalar(77,'m/s','m/s',velmax,ierr)  
      velmax = max(abs(velmax),0.0) 
    
    case('WSE_MAX','WATER_SURFACE_ELEVATION_MAXIMUM','WATER_LEVEL_MAX','WATER_LEVEL_MAXIMUM')
      call card_scalar(77,'m','m',presmax,ierr)  
      presmax = presmax*grav
      presmax = max(abs(presmax),0.0) 
      
    case('ADVECTION_SCHEME')
      backspace(77)
      read(77,*) cardname, cdum  
      do i=0,size(advsc)-1
        if(cdum==advsc(i))then
          ndsch = i
          exit
        endif
      enddo
      if(ndsch==0 .and. doPrint)then
        call diag_print_warning('Turning off advection',&
          '  Only recomended for testing or idealized cases')
      endif
          
    case('TEMPORAL_SCHEME','IMPLICIT_SCHEME')
      backspace(77)
      read(77,*) cardname, cdum  
      selectcase(cdum)
      case('THREE-LEVEL')
        wtsch = 1.0
      case default !('TWO-LEVEL')
        wtsch = 0.0
      endselect  
          
    case('IMPLICIT_WEIGHTING_FACTOR')
      backspace(77)
      read(77,*) cardname, wtsch
      wtsch = max(min(wtsch,1.0),0.0) !Apply limits        
          
    case('NUM_THREADS','OPENMP_THREADS') !Old
      backspace(77)
      read(77,*) cardname, nthr
      nthrmax = 1
!$      nthrmax = omp_get_max_threads()
      nthr = min(nthrmax,nthr)
        
    case('SKEWNESS_CORRECTION')
      call card_boolean(77,skewcor,ierr)
      
    case('VOLUME_CONSERVATION_CORRECTION','VOLUME_CORRECTION','WATER_VOLUME_CORRECTION')  
      call card_boolean(77,volcor,ierr)
      
    !=== Artificial hydrodynamics for testing =======
    case('MAKE_FLOW_1D')
      backspace(77)
      read(77,*) cardname,iFlow1D
        
    case('WATER_BALANCE','WATER_VOLUME_BALANCE')
      call card_boolean(77,flowvolbal,ierr)
    
    case default
      foundcard = .false.  
                          
    endselect
    
    return
    endsubroutine flow_cards

!***********************************************************************
    subroutine flow_init()
!   by Weiming Wu, NCCHE, Oct. 2008
!   modified by Alex Sanchez, USACE
!***********************************************************************
#include "CMS_cpp.h"
    use size_def
    use geo_def,  only: ncface,cell2cell,x,y,yc,avg_lat,areaavg,&
        xorigin,yorigin,azimuth_fl,lat,latfile,latpath
    use flow_def  
    use flow_lib
    use der_def,  only: nder
    use comvarbl
    use time_lib, only: ramp_func
    use const_def
    use secant_lib
    use prec_def
    use sed_def,  only: scalemorph_rampdur, scalemorph_ramp, scalemorph_orig, scalemorph
    implicit none
    
    integer :: iter    !omp_get_num_threads
    real(ikind) :: yavg,s(3),f(2),epssal,toldens,epsheat
        
    !Gravity-related variables
    gravinv = 1.0_ikind/grav
    sqrtgrav = sqrt(grav)
      
    !Water density
    selectcase(idensit) 
    case(0) !User-specified
      !Solve for salinity based on input temperature and density
      epssal = 0.001  !iterative change in salinity
      epsheat = 0.001 !iterative change in temperature
      toldens = 0.001 !error in density
      iter = 100      !maximum number of iterations
      s(1) = 25.0; s(2) = 35.0; !First guesses
      call secant(water_density_UNESCO,rhow,watertemp,s,f,epssal,toldens,iter)
      watersalt = s(3)
    case(1); rhow = water_density_ekart(watertemp,watersalt)
    case(2); rhow = water_density_UNESCO(watertemp,watersalt)
    case(3); rhow = water_density_fofonoff(watertemp,watersalt)    
    endselect
    
    !Kinematic viscosity, m^2/s
    selectcase(iviscos) 
    !case(0); User-specified
    case(1); viscos = water_viscosity_kinematic(watertemp)  !Temperature only
    case(2); viscos = water_viscosity_dynamic(watertemp,watersalt)/rhow !Temperature and salinity
    endselect
    
    !Timing variables  
!added meb 03/11/2019
    scalemorph_rampdur = max(scalemorph_rampdur, dtime/3600.0) !hours
    nramp = scalemorph_rampdur*3600.0/dtime  !Number of iterations
    scalemorph_ramp = 1.0
    if(nramp>1)then
      scalemorph_ramp = ramp_func(timehrs,scalemorph_rampdur)
    endif        
    if (scalemorph_orig .gt. 0.0) then 
      scalemorph = max(1.0 , scalemorph_orig*scalemorph_ramp)
    endif
!--------------------    
    
    rampdur = max(rampdur,dtime/3600.0) !hours
    nramp = rampdur*3600.0/dtime   !Number of iterations
    ramp = 1.0
    if(nramp>1)then
      ramp = ramp_func(timehrs,rampdur)
    endif        
    timehrs = stime/3600.0
    ctime = stime
    timesecs = dble(ctime)
    jtime = 0
    ntime = 0
    mtime = 0
    if(wtsch>1.0e-4  .and. nfsch==0)then
      ntsch=2
    else
      ntsch=1  
    endif                                

    !Print elapsed/projected time output
    if(nfsch==0)then
      nprt=20                   !Print elapsed output for Implicit scheme every XX timesteps
    elseif(nfsch==1)then
      nprt=4000                   !Print elapsed output for Explicit scheme every XX timesteps
    elseif(nfsch==2)then
      nprt=4000                 !Print elapsed output for Semi-implicit scheme every XX timesteps
    endif  
    
    call flow_alloc
    
    !--- Wetting and drying -------------------
    if(ncellpoly>0) narrowchannels = .true. !Always allowed for unstructured meshes
    if(.not.narrowchannels)then
      hdry1 = 5.0*hdry !One-cell-wide-channel
      hdry2 = 2.0*hdry !Two-cell-wide-channel
    endif
#ifdef DEV_MODE
    hmin = min(hmin,hdry)
#else
    hmin = hdry
#endif
    call flow_wetdry(0)
    
    !--- Time-stepping ---------------------------
    u1=u; v1=v; p1=p; h1=h
    !Implicit 3-level second-order scheme
    !This scheme is easier to impliment than the Crank-Nicolsen (CN) scheme.
    !It is also less prone to oscillations for large time steps. 
    !For small time steps however it is less accurate than the CN scheme.
    !Blending is allowed with the first order implicit Euler scheme to reduce 
    !oscillations for large time steps.
    if(ntsch==2)then
      u2 = u1; v2 = v1; p2 = p1; h2 = h1
      ctsch  = 1.0 + 0.5*wtsch
      ctsch1 = 1.0 + wtsch
      ctsch2 = 0.5*wtsch
    endif 

    !Latitudes  
    if(len_trim(latpath)>0 .and. trim(latpath)/='NONE')then
      allocate(lat(ncellsD))
      call read_latlon_dataset(lat,'Lats')  !Required because latitude dataset is saved as a group property
      avg_lat = sum(lat(1:ncells))/real(ncells,kind=ikind) !mean value
      deallocate(lat)
    endif

    !Coriolis    
    fcoriolis = 2.0*omega*sin(avg_lat*deg2rad) 
    allocate(fc(ncellsD)) 
    if(icoriolisplane==2)then !Beta-plane approximation      
      !Temporary variables
      betacoriolis = 2.0*omega*cos(avg_lat*deg2rad)/earthRadius
      !Average global y-coordinate
      if(ncellsimple>0)then
        yavg = sum(yc(1:ncells))/real(ncells,kind=ikind)
        fc(:) = fcoriolis + betacoriolis*(yc(:)-yavg)   !Beta-plane approximation
      else
        yavg = sum(y(1:ncells))/real(ncells,kind=ikind)  
        fc(:) = fcoriolis + betacoriolis*(y(:)-yavg)   !Beta-plane approximation
      endif      
    else
      fc(:) = fcoriolis    !f-plane approximation    
    endif

    !Turbulence model
    cvishor2areaavg = cvishor*cvishor*areaavg !Alex, March 31, 2011
    cvishor2areaavg2 = cvishor2areaavg**2
    if(mturbul>=3 .and. cvishor<1.0e-6) mturbul=2
    
    !Advection
    if(cviscon<1.0e-6 .and. mturbul==0 .and. ndsch>1 .and. ndsch<5)then
      ndsch = 1 !Upwind scheme
    endif
    
    !Settings for solver 
    if(ncelljoint>0 .and. (nsolv==6 .or. nsolv==7))then
      nsolv = 4
    endif
    selectcase(nsolv)
    case(0) !ADI
      if(nswp0(1)==0) nswp(1)=40        
      if(nswp0(2)==0) nswp(2)=10
      if(nswp0(3)==0) nswp(3)=10
      if(nswp0(4)==0) nswp(4)=10
      if(nswp0(5)==0) nswp(5)=10    
    case(1) !Gauss-Seidel
      if(nswp0(1)==0) nswp(1)=100
      if(nswp0(2)==0) nswp(2)=30
      if(nswp0(3)==0) nswp(3)=30
      if(nswp0(4)==0) nswp(4)=15
      if(nswp0(5)==0) nswp(5)=15
      if(maxit0==0)   maxit=30+ncells/5000
    case(2) !Gauss-Seidel-SOR
      if(nswp0(1)==0) nswp(1)=80        
      if(nswp0(2)==0) nswp(2)=20
      if(nswp0(3)==0) nswp(3)=20
      if(nswp0(4)==0) nswp(4)=15
      if(nswp0(5)==0) nswp(5)=15        
!      nswp(1)=100
!      nswp(2)=15
!      nswp(3)=15
!      nswp(4)=15
!      nswp(5)=10
      if(maxit0==0)   maxit=30+ncells/5000
    case(3) !BiCGSTAB
      if(nswp0(1)==0) nswp(1)=15
      if(nswp0(2)==0) nswp(2)=5
      if(nswp0(3)==0) nswp(3)=5
      if(nswp0(4)==0) nswp(4)=5
      if(nswp0(5)==0) nswp(5)=5
      if(maxit0==0)   maxit=20+ncells/10000
    case(4) !GMRES
      if(nswp0(1)==0) nswp(1)=20
      if(nswp0(2)==0) nswp(2)=5
      if(nswp0(3)==0) nswp(3)=5
      if(nswp0(4)==0) nswp(4)=3
!     if(nswp0(4)==0) nswp(4)=5
      if(nswp0(5)==0) nswp(5)=3
      if(maxit0==0)   maxit=20+ncells/10000 
    case(5) !Hybrid
      if(nswp0(1)==0) nswp(1)=20
      if(nswp0(2)==0) nswp(2)=5
      if(nswp0(3)==0) nswp(3)=5
      if(nswp0(4)==0) nswp(4)=3
!     if(nswp0(4)==0) nswp(4)=5
      if(nswp0(5)==0) nswp(5)=3
     if(maxit0==0)    maxit=20+ncells/10000  
    case(6) !SIP
      if(nswp0(1)==0) nswp(1)=30
      if(nswp0(2)==0) nswp(2)=7
      if(nswp0(3)==0) nswp(3)=7
      if(nswp0(4)==0) nswp(4)=5
      if(nswp0(5)==0) nswp(5)=5
      if(maxit0==0)   maxit=20+ncells/10000
    case(7) !ICCG
      if(nswp0(1)==0) nswp(1)=30
      if(nswp0(2)==0) nswp(2)=7
      if(nswp0(3)==0) nswp(3)=7
      if(nswp0(4)==0) nswp(4)=5
      if(nswp0(4)==0) nswp(5)=5
      if(maxit0==0)   maxit=20+ncells/10000
    case(8) !ICCGSTAB !Still under testing
      if(nswp0(1)==0) nswp(1)=30
      if(nswp0(2)==0) nswp(2)=7
      if(nswp0(3)==0) nswp(3)=7
      if(nswp0(4)==0) nswp(4)=5
      if(nswp0(5)==0) nswp(5)=5
      if(maxit0==0)   maxit=20+ncells/10000      
    endselect      
    nswp0=nswp      
    maxit0=maxit      
    rmom0=1.e10
    rmom=0.0
  
!$  if(nthr>=1) call omp_set_num_threads(nthr)
!!$  nthr = omp_get_num_threads() 

    !Skewness correction only for telescoping grids
    if(ncelljoint==0 .and. ncellpoly==0) skewcor=.false.

    return
    endsubroutine flow_init

!*****************************************************************
    subroutine flow_alloc()
! Allocates and initializes the flow variables
!
! written by Weiming Wu, NCCHE, and
!         Alex Sanchez, USACE-CHL
!*****************************************************************    
#include "CMS_cpp.h"
    use size_def
    use geo_def
    use flow_def
    use prec_def
    use comvarbl, only: ntsch,nfsch,pred_corr,norder
    implicit none
  
    !State
    allocate(iwet(ncellsD),iwet1(ncellsD))
    iwet=1; iwet1=1
    allocate(icorner(ncellsD))
    icorner=0
         
    !Shared implicit and explicit flow variables
    allocate(flux(nmaxfaces,ncellsD),flux1(nmaxfaces,ncellsD))
    allocate(p(ncellsD),p1(ncellsD),pk(nmaxfaces,ncellsD))
    allocate(dpx(ncellsD),dpy(ncellsD))
    allocate(eta(ncellsD),uv(ncellsD))
    allocate(h(ncellsD),h1(ncellsD)) 
    allocate(hk(nmaxfaces,ncellsD))
    allocate(u(ncellsD),u1(ncellsD),dux(ncellsD),duy(ncellsD))    
    allocate(v(ncellsD),v1(ncellsD),dvx(ncellsD),dvy(ncellsD))         
    flux=0.0; flux1=0.0
    p=0.0; p1=0.0; pk=0.0; dpx=0.0; dpy=0.0
    eta=0.0; uv=0.0
    h=0.0; h1=0.0; hk=0.0
    u=0.0; u1=0.0; dux=0.0; duy=0.0
    v=0.0; v1=0.0; dvx=0.0; dvy=0.0        
    
    allocate(su(ncellsD),sv(ncellsD),sp(ncellsD)) 
    su=0.0; sv=0.0; sp=0.0
    
    if(flowvolbal)then
      allocate(h0(ncellsD))
      h0 = 0.0
      volH2Ocumstg = 0.0
      volH2Ocumbnd = 0.0
      volH2Ocumrain = 0.0
      volH2Ocumevap = 0.0
      volH2Ocumflux = 0.0
    endif
    
    !-------------------------------------------------------------------
    !allocate(ppk(ncellsD,nmaxfaces))    
    !allocate(uk(ncellsD,nmaxfaces),vk(ncellsD,nmaxfaces))
    !ppk=0.0
    !uk=0.0; vk=0.0
    !-------------------------------------------------------------------
         
    !Implicit solver
    if(nfsch==0)then
      !Coefficient matrix
      allocate(acoef(nmaxfaces,ncellsD),ap(ncells))
      acoef=0.0; ap=0.0
      allocate(spu(ncellsD),spv(ncellsD)) 
      spu=0.0; spv=0.0
      allocate(ssu0(ncellsD),ssv0(ncellsD),sspp0(ncellsD)) 
      ssu0=0.0; ssv0=0.0; sspp0=0.0
      allocate(spuv0(ncellsD),sppp0(ncellsD))
      spuv0=0.0; sppp0=0.0
      
      !Pressure correction and momentum interpolation
      allocate(pp(ncellsD),dppx(ncellsD),dppy(ncellsD))
      pp=0.0; dppx=0.0; dppy=0.0
      allocate(sumu(ncellsD)) !dsumux,dsumuy are not used
      sumu=0.0
      allocate(apuareap(ncellsD),dapuareapx(ncellsD),dapuareapy(ncellsD)) 
      apuareap=0.0; dapuareapx=0.0; dapuareapy=0.0
      allocate(Hu(ncellsD),dHux(ncellsD),dHuy(ncellsD))
      Hu=0.0; dHux=0.0; dHuy=0.0
      allocate(Hv(ncellsD),dHvx(ncellsD),dHvy(ncellsD))
      Hv=0.0; dHvx=0.0; dHvy=0.0
      
      !Normalized residuals
      allocate(rsp(ncellsD),rsu(ncellsD),rsv(ncellsD))             
      rsp=0.0; rsu=0.0; rsv=0.0

      if(ntsch==2)then
        allocate(u2(ncellsD),v2(ncellsD),p2(ncellsD),h2(ncellsD))
        u2=0.0; v2=0.0; p2=0.0; h2=h
        !allocate(iwet2(ncellsD))
        !iwet2=1
      endif
    elseif(nfsch==1)then !Explicit solver
      !allocate(pk(nmaxfaces,ncellsD))
      !pk=0.0
      !allocate(uk(nmaxfaces,ncellsD),vk(nmaxfaces,ncellsD))      
      !uk=0.0;  vk=0.0      
      allocate(Huk(nmaxfaces,ncellsD),Hvk(nmaxfaces,ncellsD)) !Total Fluxes for u and v equations       
      Huk=0.0; Hvk=0.0
      allocate(detax(ncellsD),detay(ncellsD))
        detax=0.0; detay=0.0
      !if(norder==2)then
        allocate(dhx(ncellsD),dhy(ncellsD))
        dhx=0.0; dhy=0.0
        !allocate(limdhx(ncellsD),limdhy(ncellsD))
        !limdhx=0.0; limdhy=0.0
        !allocate(limdetax(ncellsD),limdetay(ncellsD))
        !limdetax=0.0; limdetay=0.0
        !allocate(limdux(ncellsD),limduy(ncellsD),limdvx(ncellsD),limdvy(ncellsD))
        !limdux=0.0; limduy=0.0; limdvx=0.0; limdvy=0.0
      !endif
    elseif(nfsch==2)then !Semi-implicit solver
      allocate(Huk(nmaxfaces,ncellsD),Hvk(nmaxfaces,ncellsD)) !Total Fluxes for u and v equations       
      Huk=0.0; Hvk=0.0
      allocate(detax(ncellsD),detay(ncellsD))
      detax=0.0; detay=0.0
      allocate(acoef(nmaxfaces,ncellsD),ap(ncells))
      acoef=0.0; ap=0.0
      allocate(hustar(ncellsD),dhustarx(ncellsD),dhustary(ncellsD))
      allocate(hvstar(ncellsD),dhvstarx(ncellsD),dhvstary(ncellsD))
      hustar=0.0; dhustarx=0.0; dhustary=0.0
      hvstar=0.0; dhvstarx=0.0; dhvstary=0.0
      allocate(spu(ncellsD),spv(ncellsD)) 
      spu=0.0; spv=0.0
      allocate(fluxstar(nmaxfaces,ncellsD))
      fluxstar=0.0
    endif 
         
    !Eddy viscosity
    allocate(vis(ncellsD),viskfl(nmaxfaces,ncellsD),visk(nmaxfaces,ncellsD))
    vis=0.0; visk=0.0; viskfl=0.0
    if(mturbul==4)then !Mixing-length model
      allocate(diswall(ncellsD))
      diswall=100000000.0  
    endif  
    
    !Stokes Velocities    
    allocate(us(ncellsD),vs(ncellsD))
    us=0.0; vs=0.0
    
#ifdef DEV_MODE
    !Predictor-Corrector Method (under development)    
    if(pred_corr)then
      allocate(hpred(ncellsD),upred(ncellsD),vpred(ncellsD))
      hpred=h; upred=0.0; vpred=0.0
    endif
#endif
    
    return
    endsubroutine flow_alloc
    
!***************************************************    
    subroutine flow_cleanup()
! Deallocates memory for flow module
! written by Alex Sanchez, USACE-CHL
!***************************************************
    use flow_def
    use comvarbl
    implicit none
    
    deallocate(iwet,iwet1)
    deallocate(icorner)
    deallocate(eta,uv,flux,flux1)    
    deallocate(h,h1,hk)
    if(allocated(dhx)) deallocate(dhx,dhy)
    if(allocated(ppk)) deallocate(ppk)
    deallocate(p,p1,dpx,dpy,pk)
    deallocate(u,u1,dux,duy)
    deallocate(v,v1,dvx,dvy)  
    
    !Implicit solver
    if(nfsch==0)then
      !Coefficient matrix
      deallocate(acoef,ap)
      deallocate(su,sv,sp,spu,spv)
      deallocate(ssu0,ssv0,sspp0) 
      deallocate(spuv0,sppp0)
      
      !Pressure correction and momentum interpolation
      deallocate(pp,dppx,dppy)
      deallocate(sumu) !dsumux,dsumuy are not used
      deallocate(apuareap,dapuareapx,dapuareapy) 
      deallocate(Hu,dHux,dHuy)
      deallocate(Hv,dHvx,dHvy)
      
      !Normalized residuals
      deallocate(rsp,rsu,rsv)

      !Second-order scheme
      if(ntsch==2)then
        deallocate(u2,v2,p2,h2)
        if(allocated(iwet2)) deallocate(iwet2)
      endif
    else !Explicit solver
      if(allocated(uk)) deallocate(uk)
      if(allocated(vk)) deallocate(vk)
      deallocate(Huk,Hvk) !Total Fluxes for u and v equations
    endif 
    
    !Eddy viscosity
    deallocate(vis,viskfl,visk)
    if(mturbul==4)then !Mixing-length model
      deallocate(diswall)
    endif 
    
    return
    endsubroutine flow_cleanup    
    
!**************************************************
    subroutine flow_print
! Prints the hydro setup parameters to the screen 
! and diagnositic file. 
! written by Alex Sanchez, USACE-CHL
!**************************************************   
#include "CMS_cpp.h"
    use size_def
    use geo_def
    use flow_def
    use der_def
    use diag_def
    use comvarbl
    use solv_def, only: asolv
    use rol_def,  only: rolflux
    use cms_def,  only: noptset
    use tool_def, only: vstrlz
    implicit none
    integer :: i,iunit(2)

111  format(' ',A)
222  format(' ',A,T40,A)    
342  format(' ',A,T40,F0.2,A)
353  format(' ',A,T40,F0.3,A)
354  format(' ',A,T40,A,A)    !Added for vstrlz function results
345  format(' ',A,T40,F0.5)
787  format(' ',A,T40,1x,A)
788  format(' ',A,T40,1x,A,A)
800  format(' ',A,T40,1pe10.3,A)    
801  format(' ',A,T40,I0)
5431 format(' ','    Simulation Start Time: 'T40,I4,'-',I2.2,'-',I2.2,' ',I2.2,':',I2.2,':',I2.2,' UTC')
    
!Note: There is no float format for a leading zero if value is <1.  
!  I am using the function 'vstrlz' to convert to a string for some of those cases.
     
    iunit = (/6, dgunit/)
    open(dgunit,file=dgfile,access='append') 
    
    do i=1,2    
      write(iunit(i),*)        
      write(iunit(i),111)       'Hydrodynamics'  
      write(iunit(i),111)       '  Water Properties:'
      write(iunit(i),342)       '    Temperature:',watertemp,' deg C'
      if(idensit>0 .or. iviscos==2)then
        write(iunit(i),342)     '    Salinity:',watersalt,' ppt'
      endif
      write(iunit(i),342)       '    Density:',rhow,' kg/m^3'
      write(iunit(i),354)       '    Kinematic Viscosity: ',trim(vstrlz(viscos,'(1pe10.3)')),' m^2/s'
      write(iunit(i),111)       '  Timing'
      write(iunit(i),5431)      iyr,imo,iday,ihr,imin,isec
      write(iunit(i),354)       '    Hydrodynamic time step: ',trim(vstrlz(dtime,'(f0.3)')),' sec'
      write(iunit(i),354)       '    Simulation Duration: ',trim(vstrlz(tmax,'(f0.3)')),' hours'
      write(iunit(i),354)       '    Ramp Duration: ',trim(vstrlz(rampdur,'(f0.3)')),' hours'    

      !write(iunit(i),*)    
      write(iunit(i),111)        '  Wetting and Drying'
#ifdef DEV_MODE
      write(iunit(i),354)       '    Minimum Depth: ',trim(vstrlz(hmin,'(1pe10.3)')),' m'
#endif
      write(iunit(i),354)       '    Drying Depth: ',trim(vstrlz(hdry,'(f0.3)')),' m'    
      if(ponding)then
        write(iunit(i),222)     '    Water Ponding: ','ON'
      else
        write(iunit(i),222)     '    Water Ponding: ','OFF'
      endif
      if(narrowchannels)then
        write(iunit(i),222)     '    Narrow Channels: ','ON'
        write(iunit(i),354)     '    One-cell Drying Depth: ',trim(vstrlz(hdry1,'(f0.3)')),' m'    
        write(iunit(i),354)     '    Two-cell Drying Depth: ',trim(vstrlz(hdry2,'(f0.3)')),' m'    
      else
        write(iunit(i),222)     '    Narrow Channels:','OFF'
      endif  
      
      if(noptset>=3)then
        if(waveflux)then
          write(iunit(i),222)   '  Wave Mass Flux:','ON'
          if(rolflux)then
            write(iunit(i),222) '  Roller Mass Flux:','ON'
          else
            write(iunit(i),222) '  Roller Mass Flux:','OFF'  
          endif
        else
          write(iunit(i),222)   '  Wave Mass Flux:','OFF'
        endif
      endif
    
      if(icoriolisplane==1)then
        write(iunit(i),222) '  Coriolis Approximation:','F-PLANE'  
        write(iunit(i),354) '    Average Latitude:',trim(vstrlz(avg_lat,'(f0.3)')),' deg'
        write(iunit(i),354) '    Constant Value:',trim(vstrlz(fcoriolis,'(1pe10.3)')),' 1/s'
      else
        write(iunit(i),222) '  Coriolis Approximation:','BETA-PLANE'
        write(iunit(i),354) '    Average Latitude:',trim(vstrlz(avg_lat,'(f0.3)')),' deg'
        write(iunit(i),354) '    Central Value:',trim(vstrlz(fcoriolis,'(1pe10.3)')),' 1/s'
      endif
    
      write(iunit(i),222)   '  Turbulence Model:',trim(aturb(mturbul))
      write(iunit(i),111)   '    Coefficients'
      write(iunit(i),354)   '    Constant:',trim(vstrlz(cviscon,'(1pe10.3)'))
      if(mturbul>=2)then
        write(iunit(i),354) '    Current Bottom Shear:',trim(vstrlz(cvisbot,'(f0.3)'))
      endif  
      if(mturbul>=3)then
        write(iunit(i),354) '    Current Horizontal Shear: ',trim(vstrlz(cvishor,'(f0.3)'))
      endif  
      if(noptset>=3)then
        write(iunit(i),354) '    Wave Bottom Shear:',trim(vstrlz(cviswav,'(f0.3)'))
        write(iunit(i),354) '    Wave Breaking:',trim(vstrlz(cviswavbrk,'(f0.3)'))
      endif    
    
      write(iunit(i),111)   '  Numerical Methods'
      if(nfsch==0)then !Implicit
        write(iunit(i),222) '    Solution Scheme:','IMPLICIT'
        if(dtime<1.0)then !Check time step
          write(iunit(i),*)     
          write(iunit(i),111) '      WARNING: Extremely small time step for implicit scheme'
        endif
#ifdef DEV_MODE
        if(pred_corr)then
          write(iunit(i),222) '    Predictor-Corrector Scheme:','ON'
        else
          write(iunit(i),222) '    Predictor-Corrector Scheme:','OFF'
        endif
#endif
        if(ntsch==1)then
          write(iunit(i),222) '    Temporal Scheme:','TWO-LEVEL'
        else
          write(iunit(i),222) '    Temporal Scheme:','THREE-LEVEL'
          write(iunit(i),354) '    Implicit Weighting Factor:',trim(vstrlz(wtsch,'(f0.3)'))
        endif
        write(iunit(i),222)   '    Advection Scheme:',trim(advsc(ndsch))
        write(iunit(i),222)   '    Matrix Solver:',trim(asolv(nsolv))
        if(nsolv==2 .or. nsolv==5)then
          write(iunit(i),354) '    SOR Constant:',trim(vstrlz(relaxsor,'(f0.3)'))
        endif   
        write(iunit(i),801)   '    Maximum Solver Iterations:',maxit
        write(iunit(i),801)   '    Pressure Iterations:',nswp(1)
        write(iunit(i),801)   '    Velocity Iterations:',nswp(2)   
        write(iunit(i),354)   '    Pressure Relaxation:',trim(vstrlz(relax(8),'(f0.3)'))
        write(iunit(i),354)   '    Velocity Relaxation:',trim(vstrlz(relax(2),'(f0.3)'))
        write(iunit(i),354)   '    Pressure Max Residual:',trim(vstrlz(rmommaxp,'(1pe10.3)'))
        write(iunit(i),354)   '    Pressure Target Residual:',trim(vstrlz(rmomtargetp,'(1pe10.3)'))
        write(iunit(i),354)   '    Pressure Min Residual:',trim(vstrlz(rmomminp,'(1pe10.3)'))
        write(iunit(i),354)   '    Velocity Max Residual:',trim(vstrlz(rmommaxuv,'(1pe10.3)'))
        write(iunit(i),354)   '    Velocity Target Residual:',trim(vstrlz(rmomtargetuv,'(1pe10.3)'))
        write(iunit(i),354)   '    Velocity Min Residual:',trim(vstrlz(rmomminuv,'(1pe10.3)'))
        write(iunit(i),354)   '    Water Level Max:',trim(vstrlz(presmax*gravinv,'(f0.3)')), ' m'
        write(iunit(i),354)   '    Velocity Max:',trim(vstrlz(velmax,'(f0.3)')),' m/s'
        if(skewcor)then
          write(iunit(i),222) '    Skewness Correction:','ON'
        else
          write(iunit(i),222) '    Skewness Correction:','OFF'
        endif
        write(iunit(i),222)   '    Slope Limiter:',trim(alim(nlim))
        write(iunit(i),222)   '    Spatial Derivative Scheme:',trim(ader(nder))
      elseif(nfsch==1)then    !Explicit
        write(iunit(i),222)   '  Solution Scheme:','EXPLICIT'
        write(iunit(i),801)   '    Order Accuracy:',norder
        write(iunit(i),222)   '    Riemann Solver:',trim(ariem(nriem))
        write(iunit(i),222)   '    Spatial Derivative Scheme:',trim(ader(nder))
        write(iunit(i),222)   '    Slope Limiter:',trim(alim(nlim))
      else !if(nfsch==1)then  !Semi-implicit  
        write(iunit(i),222)   '  Solution Scheme:','SEMI_IMPLICIT'
        write(iunit(i),801)   '    Order Accuracy:',norder  
        write(iunit(i),222)   '    Riemann Solver:',trim(ariem(nriem))
        write(iunit(i),222)   '    Spatial Derivative Scheme:',trim(ader(nder))
        write(iunit(i),222)   '    Slope Limiter:',trim(alim(nlim))
      endif   
      write(iunit(i),801)     '    Maximum Number of Threads:',nthrmax
      write(iunit(i),801)     '    Number of Threads used:',nthr     
    enddo
                
    close(dgunit)            
    
    return
    endsubroutine flow_print

!***********************************************************************
    subroutine flow_wetdry(iflagdry)
! Judge wetting and drying
! written by Weiming Wu, NCCHE, Oct. 2008
! modified by Alex Sanchez, USACE-CHL, 2013
!***********************************************************************
#include "CMS_cpp.h"
    use size_def
    use geo_def
    use flow_def
    use bnd_def
    use interp_def, only: fintp
    use comvarbl, only: dtime,pred_corr
    use const_def, only: small
    use diag_lib
    use prec_def
    implicit none
    integer :: i,j,k,kk,ksum,nck,iflagdry   !ibnd
    integer :: jcn,iwetface,iwetface1
    integer :: kew,keww,kns,knss,knns,keew,ndrysum,nsumi,nddryswitch
    real(ikind) :: usumi,vsumi,vissumi
    real(ikind) :: uni,unk,fp,fk
    
    !hdry=hmin  ! added by Wu temp
    
    if(iflagdry==0)then
!$OMP PARALLEL DO PRIVATE(i)         
      do i=1,ncellsD
        h(i)=p(i)*gravinv-zb(i)
        if(h(i)<hdry)then
          !h(i)=hdry
          h(i)=max(h(i),hmin)
          iwet(i)=0
        else
          iwet(i)=1
        endif
      enddo
!$OMP END PARALLEL DO         
    elseif(iflagdry==1)then
!$OMP PARALLEL DO PRIVATE(i)
      do i=1,ncellsD
        h(i)=p(i)*gravinv-zb(i)
        if(h(i)<hdry)then
          !h(i)=hdry
          h(i)=max(h(i),hmin)
          iwet(i)=0
        endif
      enddo
!$OMP END PARALLEL DO      
    endif
    
   ! !Cells are first assigned a wet/dry condition based on the 
   ! !cell average water depth. Then the neighboring cells are checked
   ! !for greater water levels 
   ! do i=1,ncells
      !do j=1,nxyface(i)
      !  k = kxyface(j,i)
      !  nck = cell2cell(k,i)
   !     if(nck>ncells) cycle
   !     if((h(i)>hdry   .and. h(nck)<hdry .and. &
   !         p(i)>p(nck) .and. (p(i)*gravinv-zb(nck))>hdry))then   !Right dry but lower than left
   !       iwet(nck) = 1
   !       h(nck) = max(h(nck),hmin)
   !     elseif((h(i)<hdry   .and. h(nck)>hdry .and. &
   !             p(i)<p(nck) .and. (p(nck)*gravinv-zb(i))>hdry))then   !Left dry but lower than right
   !       iwet(i) = 1
   !       h(i) = max(h(i),hmin)
   !     endif
   !   enddo !j-face
   ! enddo !i-cell

!--- Ghost cells -------------------------------------------    
    !Closed boundaries
    do j=1,W_str%ncells
      i=W_str%cells(j)
      k=W_str%faces(j)
      iwet(cell2cell(k,i))=0 !dry
    enddo     
    !!Open boundaries
    !do ibnd=1,nbndstr    !Closed by Wu  2/16/2014
    !  do j=1,bnd_str(ibnd)%ncells
    !    i=bnd_str(ibnd)%cells(j)
    !    k=bnd_str(ibnd)%faces(j)
    !    iwet(cell2cell(k,i))=iwet(i) !Copy to ghost cells
    !  enddo
    !enddo    

!--- Structures --------------------------------------    
    call struct_wetdry

    !narrowchannels=.false.  !Added by Wu
    
!--- Special cases -------------------------------------
    if(narrowchannels)then !Allow one-cell-wide channels (always allowed for unstructured meshes)
      nddryswitch=1
      do while(nddryswitch>0)
        nddryswitch=0
!$OMP PARALLEL DO PRIVATE(i,ksum) REDUCTION(+:nddryswitch)
        do i=1,ncells
          if(iwet(i)==0) cycle
          ksum=sum(iwet(cell2cell(1:ncface(i),i)))
          if(ksum==0)then
            iwet(i)=0
            nddryswitch=nddryswitch+1
          endif
        enddo !i-cells
!$OMP END PARALLEL DO      
        !continue
      enddo !while  
    else !Remove narrow channels (only for Cartesian grids)
      nddryswitch=1      
      do while(nddryswitch>0)
        nddryswitch=0
      if(ncellpoly==0)then  !Added by Wu
        do i=1,ncells
          if(iwet(i)==0) cycle
          !One wet cell surrounded by dry cells
          ksum=sum(iwet(cell2cell(1:ncface(i),i)))
          if(ksum<=1)then
            iwet(i)=0
            nddryswitch=nddryswitch+1
          endif
          !1 wet cell between dry cells treated as dry
          kns=0
          kew=0
          do k=1,ncface(i)
            nck=cell2cell(k,i)
            if(idirface(k,i)==1 .or. idirface(k,i)==3) kns=kns+iwet(nck)
            if(idirface(k,i)==2 .or. idirface(k,i)==4) kew=kew+iwet(nck)
          enddo
          !if(kns==0 .or. kew==0)then
          if((kns==0 .or. kew==0) .and. h(i)<hdry1)then    
            iwet(i)=0
            nddryswitch=nddryswitch+1
          endif
          !Two wet cells between dry cells are treated as dry cells
          knss=0
          do k=1,ncface(i)
            nck=cell2cell(k,i)
            if(idirface(k,i)==1) knss=knss+iwet(nck)
            if(idirface(k,i)==3)then
              do kk=1,ncface(nck)
                if(idirface(kk,nck)==3) knss=knss+iwet(cell2cell(kk,nck))
              enddo
            endif
          enddo
          !if(knss==0)then
          if(knss==0 .and. h(i)<hdry2)then    
            iwet(i)=0
            nddryswitch=nddryswitch+1
          endif
          knns=0
          do k=1,ncface(i)
            nck=cell2cell(k,i)
            if(idirface(k,i)==3) knns=knns+iwet(nck)
            if(idirface(k,i)==1)then
              do kk=1,ncface(nck)
                if(idirface(kk,nck)==1) knns=knns+iwet(cell2cell(kk,nck))
              enddo
            endif
          enddo
          !if(knns==0)then
          if(knss==0 .and. h(i)<hdry2)then     
            iwet(i)=0
            nddryswitch=nddryswitch+1
          endif
          keww=0
          do k=1,ncface(i)
            nck=cell2cell(k,i)
            if(idirface(k,i)==2) keww=keww+iwet(nck)
            if(idirface(k,i)==4)then
              do kk=1,ncface(nck)
                if(idirface(kk,nck)==4) keww=keww+iwet(cell2cell(kk,nck))
              enddo
            endif
          enddo
          !if(keww==0)then
          if(keww==0 .and. h(i)<hdry2)then    
            iwet(i)=0
            nddryswitch=nddryswitch+1
          endif             
          keew=0
          do k=1,ncface(i)
            nck=cell2cell(k,i)
            if(idirface(k,i)==4) keew=keew+iwet(nck)
            if(idirface(k,i)==2)then
              do kk=1,ncface(nck)
                if(idirface(kk,nck)==2) keew=keew+iwet(cell2cell(kk,nck))
              enddo
            endif
          enddo
          !if(keew==0)then
          if(keew==0 .and. h(i)<hdry2)then    
            iwet(i)=0
            nddryswitch=nddryswitch+1
          endif
        enddo !i-cells
      else  !if(ncellpoly>0) then   !Added by Wu
        do i=1,ncells
          if(iwet(i)==0) cycle
          !One wet cell surrounded by dry cells
          ksum=sum(iwet(cell2cell(1:ncface(i),i)))
          if(ksum<=1)then
            iwet(i)=0
            nddryswitch=nddryswitch+1
          endif
          !1 wet cell between dry cells treated as dry
          if(ncface(i)==4) then
            kns=0
            kew=0
            do k=1,ncface(i)
              nck=cell2cell(k,i)
              if(k==1 .or. k==3) kns=kns+iwet(nck)
              if(k==2 .or. k==4) kew=kew+iwet(nck)
            enddo
            if(kns==0 .or. kew==0)then
            !if((kns==0 .or. kew==0) .and. h(i)<hdry1)then    
              iwet(i)=0
              nddryswitch=nddryswitch+1
            endif
          endif
        enddo !i-cells
      endif               
      enddo !while
    endif !narrowchannels
   
    if(.not.ponding) call flow_pond

!--- Find corner wet nodes near dry nodes -------------------
    if(ncellpoly==0)then
!$OMP PARALLEL DO SHARED(icorner,iwet), PRIVATE(i,k,ndrysum)    
      do i=1,ncells
        icorner(i)=0
        if(iwet(i)==1) then
          ndrysum=0
          do k=1,ncface(i)
            if(iwet(cell2cell(k,i))==0) ndrysum=ndrysum+1
          enddo
          if(ndrysum>1) icorner(i)=1
        endif
      enddo
!$OMP END PARALLEL DO
    endif

!--- Initialize dry-to-wet nodes ---------------------------
    do j=1,5
!$OMP PARALLEL DO PRIVATE(i,usumi,vsumi,vissumi,nsumi)
      do i=1,ncells
        if(iwet(i)==1 .and. iwet1(i)==0)then
          usumi=sum(u(cell2cell(1:ncface(i),i)))
          vsumi=sum(v(cell2cell(1:ncface(i),i)))
          vissumi=sum(vis(cell2cell(1:ncface(i),i)))
          nsumi=sum(iwet(cell2cell(1:ncface(i),i)))            
          u(i)=usumi/(real(nsumi)+small)/1.25
          v(i)=vsumi/(real(nsumi)+small)/1.25
          vis(i)=vissumi/(real(nsumi)+small)
        endif
      enddo
!$OMP END PARALLEL DO
    enddo

!--- Wall distance (for mixing-length model) -----------------------------      
    if(mturbul==4) call distance2wall(diswall)    

!--- Apply new wet/dry criteria to fluxes and viscosity ----    
!Note: Because the fluxes and viscosities at the cell faces
! are calculated after one complete SIMPLEC cycles they should
! be initialized here for the cell faces which change wet/dry state
! from the previous iteration.
    do i=1,ncells
      do j=1,nxyface(i) !No repeat cell faces
        k=kxyface(j,i)      
        nck=cell2cell(k,i)
        jcn=llec2llec(k,i)
        iwetface=iwet(i)*iwet(nck)
        iwetface1=iwet1(i)*iwet1(nck)
        !if(iwetface/=iwetface1)then !Only change the faces that change state
           if(iwetface==0)then 
             flux(k,i)=0.0
             visk(k,i)=0.0
             flux(jcn,nck)=0.0
             visk(jcn,nck)=0.0
           elseif(iwetface==1 .and. iwetface==0)then !Initialize flux and viscosity
             uni=fnx(k,i)*u(i)+fny(k,i)*v(i)
             unk=fnx(k,i)*u(nck)+fny(k,i)*v(nck)
             fk=fintp(k,i); fp=1.0-fk
             flux(k,i)=ds(k,i)*hk(k,i)*(fk*unk+fp*uni) !Outward flux 
             visk(k,i)=fk*vis(nck)+fp*vis(i)
             flux(jcn,nck)=-flux(k,i)
             visk(jcn,nck)=visk(k,i)
           endif
         !endif
      enddo
    enddo
    
    !Check fluxes
#ifdef DIAG_MODE
    do i=1,ncells
      do k=1,ncface(i)
        if(iwet(i)*iwet(cell2cell(k,i))==0 .and. abs(flux(k,i))>1.0e-15)then
          call diag_print_error('Non-zero flux between wet and dry cells after flow_wetdry')
        endif
      enddo
    enddo
#endif

    return
    endsubroutine flow_wetdry
    
!***********************************************************************
    subroutine flow_wetdry0()
! Judge wetting and drying
! written by Weiming Wu, NCCHE, Oct. 2008
! modified by Alex Sanchez, USACE-CHL, 2013
!***********************************************************************
#include "CMS_cpp.h"
    use size_def
    use geo_def
    use flow_def
    use bnd_def
    use interp_def, only: fintp
    use comvarbl, only: dtime,pred_corr
    use const_def, only: small
    use diag_lib
    use prec_def
    implicit none
    integer :: i,j,k,kk,ksum,nck,ibnd
    integer :: jcn,iwetface,iwetface1
    integer :: kew,keww,kns,knss,knns,keew,ndrysum,nsumi,nddryswitch
    real(ikind) :: usumi,vsumi,vissumi
    real(ikind) :: uni,unk,fp,fk

#ifdef DEV_MODE
    if(pred_corr)then
!$OMP PARALLEL DO PRIVATE(i)
      do i=1,ncellsD
        h(i)=p(i)*gravinv-zb(i)
        if(hpred(i)<hdry)then !Use predicted water depth to judge wetting and drying
          h(i)=max(h(i),hmin)
          iwet(i)=0
        else
          iwet(i)=1
        endif
      enddo
!$OMP END PARALLEL DO
    else
#endif
!$OMP PARALLEL DO PRIVATE(i)
      do i=1,ncellsD
        !if(i==501)then
        !  continue
        !endif
        h(i)=p(i)*gravinv-zb(i)
        if(h(i)<hdry)then
          !h(i)=hdry
          h(i)=max(h(i),hmin)
          iwet(i)=0
        else
          iwet(i)=1
        endif
      enddo
!$OMP END PARALLEL DO
#ifdef DEV_MODE
    endif
#endif
    
!--- Ghost cells -------------------------------------------    
    !Closed boundaries
    do j=1,W_str%ncells
      i=W_str%cells(j)
      k=W_str%faces(j)
      iwet(cell2cell(k,i))=0 !dry
    enddo
    !Open boundaries
    do ibnd=1,nbndstr
      do j=1,bnd_str(ibnd)%ncells
        i=bnd_str(ibnd)%cells(j)
        k=bnd_str(ibnd)%faces(j)
        iwet(cell2cell(k,i))=iwet(i) !Copy to ghost cells
      enddo
    enddo

!--- Structures --------------------------------------    
    call struct_wetdry
      
!--- Special cases -------------------------------------
    if(narrowchannels)then !Allow one-cell-wide channels (always allowed for unstructured meshes)
      nddryswitch=1
      do while(nddryswitch>0)
        nddryswitch=0
!$OMP PARALLEL DO PRIVATE(i,ksum) REDUCTION(+:nddryswitch)
        do i=1,ncells
          if(iwet(i)==0) cycle
          ksum=sum(iwet(cell2cell(1:ncface(i),i)))
          if(ksum==0)then
            iwet(i)=0
            nddryswitch=nddryswitch+1
          endif
        enddo !i-cells
!$OMP END PARALLEL DO      
        !continue
      enddo !while  
    else !Remove narrow channels (only for Cartesian grids)
      nddryswitch=1      
      do while(nddryswitch>0)
        nddryswitch=0
        do i=1,ncells
          if(iwet(i)==0) cycle
          !One wet cell surrounded by dry cells
          ksum=sum(iwet(cell2cell(1:ncface(i),i)))
          if(ksum<=1)then
            iwet(i)=0
            nddryswitch=nddryswitch+1
          endif
          !1 wet cell between dry cells treated as dry
          kns=0
          kew=0
          do k=1,ncface(i)
            nck=cell2cell(k,i)
            if(idirface(k,i)==1 .or. idirface(k,i)==3) kns=kns+iwet(nck)
            if(idirface(k,i)==2 .or. idirface(k,i)==4) kew=kew+iwet(nck)
          enddo
          !if(kns==0 .or. kew==0)then
          if((kns==0 .or. kew==0) .and. h(i)<hdry1)then    
            iwet(i)=0
            nddryswitch=nddryswitch+1
          endif
          !Two wet cells between dry cells are treated as dry cells
          knss=0
          do k=1,ncface(i)
            nck=cell2cell(k,i)
            if(idirface(k,i)==1) knss=knss+iwet(nck)
            if(idirface(k,i)==3)then
              do kk=1,ncface(nck)
                if(idirface(kk,nck)==3) knss=knss+iwet(cell2cell(kk,nck))
              enddo
            endif
          enddo
          !if(knss==0)then
          if(knss==0 .and. h(i)<hdry2)then
            iwet(i)=0
            nddryswitch=nddryswitch+1
          endif
          knns=0
          do k=1,ncface(i)
            nck=cell2cell(k,i)
            if(idirface(k,i)==3) knns=knns+iwet(nck)
            if(idirface(k,i)==1)then
              do kk=1,ncface(nck)
                if(idirface(kk,nck)==1) knns=knns+iwet(cell2cell(kk,nck))
              enddo
            endif
          enddo
          !if(knns==0)then
          if(knss==0 .and. h(i)<hdry2)then
            iwet(i)=0
            nddryswitch=nddryswitch+1
          endif
          keww=0
          do k=1,ncface(i)
            nck=cell2cell(k,i)
            if(idirface(k,i)==2) keww=keww+iwet(nck)
            if(idirface(k,i)==4)then
              do kk=1,ncface(nck)
                if(idirface(kk,nck)==4) keww=keww+iwet(cell2cell(kk,nck))
              enddo
            endif
          enddo
          !if(keww==0)then
          if(keww==0 .and. h(i)<hdry2)then
            iwet(i)=0
            nddryswitch=nddryswitch+1
          endif             
          keew=0
          do k=1,ncface(i)
            nck=cell2cell(k,i)
            if(idirface(k,i)==4) keew=keew+iwet(nck)
            if(idirface(k,i)==2)then
              do kk=1,ncface(nck)
                if(idirface(kk,nck)==2) keew=keew+iwet(cell2cell(kk,nck))
              enddo
            endif
          enddo
          !if(keew==0)then
          if(keew==0 .and. h(i)<hdry2)then
            iwet(i)=0
            nddryswitch=nddryswitch+1
          endif
        enddo !i-cells
      enddo !while
    endif !narrowchannels

    if(.not.ponding) call flow_pond

!--- Find corner wet nodes near dry nodes -------------------
    if(ncellpoly==0)then
!$OMP PARALLEL DO SHARED(icorner,iwet), PRIVATE(i,k,ndrysum)    
      do i=1,ncells
        icorner(i)=0
        if(iwet(i)==1) then
          ndrysum=0
          do k=1,ncface(i)
            if(iwet(cell2cell(k,i))==0) ndrysum=ndrysum+1
          enddo
          if(ndrysum>1) icorner(i)=1
        endif
      enddo
!$OMP END PARALLEL DO
    endif

!--- Initialize dry-to-wet nodes ---------------------------
    do j=1,5
!$OMP PARALLEL DO PRIVATE(i,usumi,vsumi,vissumi,nsumi)
      do i=1,ncells
        if(iwet(i)==1 .and. iwet1(i)==0)then
          usumi=sum(u(cell2cell(1:ncface(i),i)))
          vsumi=sum(v(cell2cell(1:ncface(i),i)))
          vissumi=sum(vis(cell2cell(1:ncface(i),i)))
          nsumi=sum(iwet(cell2cell(1:ncface(i),i)))            
          u(i)=usumi/(real(nsumi)+small)/1.25
          v(i)=vsumi/(real(nsumi)+small)/1.25
          vis(i)=vissumi/(real(nsumi)+small)
        endif
      enddo
!$OMP END PARALLEL DO
    enddo

!--- Wall distance (for mixing-length model) -----------------------------      
    if(mturbul==4) call distance2wall(diswall)    

!--- Apply new wet/dry criteria to fluxes and viscosity ----    
!Note: Because the fluxes and viscosities at the cell faces
! are calculated after one complete SIMPLEC cycles they should
! be initialized here for the cell faces which change wet/dry state
! from the previous iteration.
    do i=1,ncells
      do j=1,nxyface(i) !No repeat cell faces
        k=kxyface(j,i)      
        nck=cell2cell(k,i)
        iwetface=iwet(i)*iwet(nck)
        iwetface1=iwet1(i)*iwet1(nck)
        if(iwetface/=iwetface1)then !Only change the faces that change state
           if(iwetface==0)then 
             flux(k,i)=0.0
             visk(k,i)=0.0 
           else !if(iwetface==1)then !Initialize flux and viscosity
             uni=fnx(k,i)*u(i)+fny(k,i)*v(i)
             unk=fnx(k,i)*u(nck)+fny(k,i)*v(nck)
             fk=fintp(k,i); fp=1.0-fk
             flux(k,i)=ds(k,i)*hk(k,i)*(fk*unk+fp*uni) !Outward flux 
             visk(k,i)=fk*vis(nck)+fp*vis(i)
           endif
           jcn=llec2llec(k,i)
           flux(jcn,nck)=-flux(k,i)
           visk(jcn,nck)=visk(k,i)
         endif
      enddo
    enddo
    
    !Check fluxes
#ifdef DIAG_MODE
    do i=1,ncells
      do k=1,ncface(i)
        if(iwet(i)*iwet(cell2cell(k,i))==0 .and. abs(flux(k,i))>1.0e-15)then
          call diag_print_error('Non-zero flux between wet and dry cells after flow_wetdry')
        endif
      enddo
    enddo
#endif

    return
    endsubroutine flow_wetdry0

!***********************************************************************
    subroutine flow_pond
! Judges imponding and dead water   
! written by Wu, Dec. 30, 2001
! modified by Alex Sanchez, USACE-CHL
!***********************************************************************
    use size_def, only: ncellsD,ncells
    use geo_def, only: ncface,cell2cell
    use flow_def, only: iwet
    use struct_def
    use bnd_def
    use prec_def
    implicit none
    integer :: i,j,k,ibnd,nck,npond
    integer :: ipond(ncellsD)

!$OMP PARALLEL DO PRIVATE(i)    
    do i=1,ncellsD
      ipond(i)=0
    enddo
!$OMP END PARALLEL DO

    !All Forcing Boundaries
    do ibnd=1,nbndstr
      do j=1,bnd_str(ibnd)%ncells
        i=bnd_str(ibnd)%cells(j)
        k=bnd_str(ibnd)%faces(j)
        ipond(i)=iwet(i)
        nck=cell2cell(k,i)
        ipond(nck)=iwet(nck)
      enddo
    enddo

    !Structures
    call struct_pond(iwet,ipond)

!!........ Check Internal Nodes
!!!      do ik=1,10   !Repeat 3 times --- Temporary
    npond=10
    do while(npond>0)
      npond=0
      do i=1,ncells
        if(ipond(i)==1)then
          do k=1,ncface(i)
            nck=cell2cell(k,i)
            npond=npond+max(0, iwet(nck)-ipond(nck))
            ipond(nck)=iwet(nck)
          enddo
        endif 
      enddo
      do i=ncells,1,-1
        if(ipond(i)==1)then
          do k=1,ncface(i)
            nck=cell2cell(k,i)
            npond=npond+max(0, iwet(nck)-ipond(nck))
            ipond(nck)=iwet(nck)
          enddo
        endif 
      enddo
    enddo
    
!$OMP PARALLEL DO PRIVATE(i)    
    do i=1,ncellsD
      if(ipond(i)==0 .and. iwet(i)==1)then
        iwet(i)=0
      endif
    enddo     
!$OMP END PARALLEL DO

    return
    endsubroutine flow_pond
    
!***********************************************************************
    subroutine flow_convflux()
! Initialising volume fluxes based on cell-center velocities
! by Weiming Wu, NCCHE, Oct. 2008
!***********************************************************************
    use size_def
    use geo_def, only: cell2cell,fnx,fny,llec2llec,nxyface,kxyface,ds
    use interp_def, only: fintp
    use flow_def, only: iwet,flux,u,v,hk
    use prec_def
    implicit none
    integer :: i,j,k,nck
    real(ikind) :: uni,unk

!$OMP PARALLEL DO PRIVATE(i,j,k,nck,uni,unk)  
    do i=1,ncells
       do j=1,nxyface(i)
         k=kxyface(j,i)
         nck=cell2cell(k,i)
         uni=fnx(k,i)*u(i)+fny(k,i)*v(i)
         unk=fnx(k,i)*u(nck)+fny(k,i)*v(nck)
         flux(k,i)=iwet(i)*iwet(nck)*ds(k,i)*hk(k,i)* &
                (fintp(k,i)*unk+(1.0-fintp(k,i))*uni) !Outward flux
         flux(llec2llec(k,i),nck)=-flux(k,i)         
       enddo
    enddo
!$OMP END PARALLEL DO

    call bndflux

    return
    endsubroutine flow_convflux
    
!*********************************************************
    subroutine flow_step_stat(iprintstat)
! Calculates flow time step statistics and prints them to 
! the screen and diagnostic file    
!
! written by Alex Sanchez, USACE-CHL
!*********************************************************
    use size_def
    use geo_def, only: mapid,dx,dy,areap
    use flow_def, only: u,v,p,h,iwet,gravinv,grav
    use comvarbl, only: uxtrm,vxtrm,pxtrm,crmax,&
        idu,idv,idp,idcr,nthr,nfsch,dtime
    use diag_lib
    use prec_def
    implicit none
    integer :: i,iprintstat,ierr
!    integer :: iugradmax,ivgradmax,ietagradmax
!    real(ikind) :: ugradmax,vgradmax,etagradmax
    integer :: idut,idvt,idpt,idcrt
    real(ikind) :: uxtrmt,vxtrmt,pxtrmt,crmaxt,c,cr
    character(len=200) :: msg
    
    !Initialize
    pxtrm=0.0; pxtrmt=0.0; idp=1; idpt=1
    uxtrm=0.0; uxtrmt=0.0; idu=1; idut=1
    vxtrm=0.0; vxtrmt=0.0; idv=1; idvt=1
    crmax=0.0; crmaxt=0.0; idcr=1; idcrt=1
!$omp parallel firstprivate(uxtrmt,vxtrmt,pxtrmt,crmaxt,idut,idvt,idpt,idcrt)
!$omp do private(i,c,cr)    
    do i=1,ncells
      if(iwet(i)==1)then
        if(abs(p(i))>abs(pxtrmt))then
          pxtrmt=p(i)
          idpt=i
        endif  
        if(abs(u(i))>abs(uxtrmt))then
          uxtrmt=u(i)
          idut=i
        endif
        if(abs(v(i))>abs(vxtrmt))then
          vxtrmt=v(i)
          idvt=i
        endif
        c=sqrt(grav*h(i)) !celerity
        if(ncellsimple>0)then
          cr=dtime*max((abs(u(i))+c)/dx(i),(abs(v(i))+c)/dy(i)) !Courant #
        else
          cr=dtime*max((abs(u(i))+c)/sqrt(areap(i)),(abs(v(i))+c)/sqrt(areap(i))) !Courant #  
        endif
        if(cr>crmaxt)then
          crmaxt=cr
          idcrt=i
        endif
      endif
    enddo   
!$omp end do
!$omp critical
    if(abs(pxtrmt)>abs(pxtrm))then
      pxtrm=pxtrmt
      idp=idpt
    endif
    if(abs(uxtrmt)>abs(uxtrm))then
      uxtrm=uxtrmt
      idu=idut
    endif
    if(abs(vxtrmt)>abs(vxtrm))then
      vxtrm=vxtrmt
      idv=idvt
    endif
    if(crmaxt>crmax)then
      crmax=crmaxt
      idcr=idcrt
    endif
!$omp end critical
!$omp end parallel
    if(allocated(mapid))then
      idp  = mapid(idp)
      idu  = mapid(idu)
      idv  = mapid(idv)
      idcr = mapid(idcr)
    endif
    
    if(iprintstat==0) return
    !if(abs(pxtrm)>1.0e5 .or. abs(uxtrm)>1.0e3 .or. abs(vxtrm)>1.0e3) return    
        
727 format('  WSE(',I6,')=',F9.4,', U(',I6,')=',F7.3,', V(',I6,')=',F7.3)
838 format('  WSE(',I6,')=',F9.4,', U(',I6,')=',F7.3,', V(',I6,')=',F7.3,', CR(',I6,')=',F6.3)
    
    if(nfsch==0)then
      write(msg,727,iostat=ierr) idp,pxtrm*gravinv,idu,uxtrm,idv,vxtrm
      call diag_print_message(msg)
    else
      write(msg,838,iostat=ierr) idp,pxtrm*gravinv,idu,uxtrm,idv,vxtrm,idcr,crmax
      call diag_print_message(msg)        
    endif
       
    return
    endsubroutine flow_step_stat
    