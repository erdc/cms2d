!==========================================================================
! CMS Sediment Transport routines
!
! Contains the following:
!   General
!     sed_default - Sets the sediment transport default parameters
!     sed_cards   - Reads the sediment transport cards from the control file
!     sed_init    - Initializes the sediment transport variables
!     sed_print   - Prints the sediment transport settings to the diagnositc
!                   file and the screen
!   Size class Input
!     sedclass_block - Reads the sediment size class block from the
!                            control file
!     sedclass_resize - Resizes a size class variable
!   Bed layer Input
!     bedlay_block  - Reads a bed layer block from the control file
!     bedlay_resize - Resizes a bed layer variable
!   Other
!     sed_balance - Computes a global sediment mass balance
!     sed_wave      - Calculates the source term in the bed change equation 
!                     due to waves asymmetry, roller and undertow
!     sed_total     - Updates several sediment transport variables once
!                     a converged solution is obtained
!     sed_step_stat - Computes the sediment transport extremal statistics
!                     and outputs them to the diagnostic file and screen
!
! written by Alex Sanchez, USACE-CHL
!            Weiming Wu, NCCHE
!===========================================================================    
    
!***************************************************************************   
    subroutine sed_default
! Sets the default parameters for sediment transport variables
!
! written by Alex Sanchez, USACE-CHL; Weiming Wu, NCCHE
!***************************************************************************    
    use sed_def
    !use wave_flowgrid_def, only: waveheight, waveperiod, wavedir
    implicit none
    
    integer :: i
    logical :: found
    
    !--- Boolean variables ------------------------------------------
    sedtrans = .false.     !Calculate total-load sediment transport          
    singlesize = .true.    !Is single sized sediment transport
    variableD50 = .false.  !Used for hiding and exposure in single sizes
    transD50 = .false.     !User has specified transport d50
    calcd50  = .false.     !CMS has calculated d50 from d50 dataset
    constd50 = .false.     !User has specified constant grain size
    sedbalance = .false.   !Global sediment balance
    sedcouple = .true.     !Coupling between transport, bed change and sorting    
    sedconstmix = .false.  !Use constant mixing coefficient
    calcmorph = .false.    !Calculate bed change
    calcbedcomp = .true.   !Calculate bed composition
    
    !--- Sediment transport model ------------------
    isedmodel = 1       
    
    !--- Transport Formula ------------------------------
    icapac = 1             !Transport capacity formula    
    Awatan = 0.5           !Watanabe coefficient

    !--- CSHORE Defaults ------------------------------ added 6/7/2019 bdj
    !CSeffb = 0.03          !CSHORE coefficient  !updated to 0.003 from Brad's Branch  05/15/2020
    CSeffb = 0.003         !CSHORE coefficient
    CSblp  = 0.002         !CSHORE coefficient
    CSslp  = 0.3           !CSHORE coefficient        
    
    !--- Horizontal Mixing -------------------------
    schmidt = 1.0          !Schmidt number ***************************************    
    cmixbedload = 0.0      !Coefficient used for bedload dispersion coefficient ~5.0 ***************
    
    !--- Fall velocity --------------------------------
    iws = 1                !Sediment fall velocity formula, 1-Soulsby, 2-Wu-Wang
    
    !--- Bed-slope effects ----------------------------------
    !Diffusion term
    do_bedslope = .true.   !Calculate bed-slope term
    dcoeff = 0.1           !Bed-slope diffusion coefficient   
    !Incipient motion and rate
    ibedslope = -1      !-1-Automatic, 0-None, 1-Dey (2001), 2-Bailard (1981), 2-Wu et al. (2000) 
    betaslope = 1.6     !Bed-load slope coefficient for Bailard (1981) method, default = 1.6 ~= 1/tan(a_repose)
    effslope  = 0.02    !Suspended-load efficiency coefficient for Bailard (1981) method, default = 0.02
    
    !--- Hiding and exposure --------------------------------
    ihidexpform = -1       !-1-Automatic,0-None,1-Egiazaroff,2-Parker,3-Wu
    mhe = -1.0             !Hiding & exposure coefficient, default function of transport formula, Lund-CIRP: 1.0, Wu: 0.6 
    varsigmamin = 0.01     !Minimum value for hiding and exposure correction
    varsigmamax = 100.0    !Maximum value for hiding and exposure correction
    
    !--- Adaptation coefficient options ---------------------------------
    !Total load adaptation
    iadapttot = 1          !Adaptation coefficient method for total load
    Ltot = 10.0            !Total-load adaptation length, m          
    Ttot = 2.0             !Total-load adaptation time, sec
    Lerotot = 10.0         !Total-load erosion adaptation length, m
    Ldeptot = 5.0          !Total load deposition adaptation length, m    
    
    !Suspended load adaptation
    iadaptsus = 1          !Adaptation coefficient method for suspended load     
    alphasus = 2.0         !Adaptation coefficient for suspended load 
    Lsus = 10.0            !Suspended-load adaptation length, m          
    Tsus = 2.0             !Suspended-load adaptation time, sec
    Lerosus = 10.0         !Suspended-load erosion adaptation length, m
    Ldepsus = 5.0          !Suspended load deposition adaptation length, m    
    
    !Bed load adaptation
    iadaptbed = 1          !Adaptation coefficient method for suspended load      
    Lbed = 10.0            !Bed-load adaptation length, m          
    Tbed = 2.0             !Bed-load adaptation time, sec
    fbed = 10.0            !Depth dependant factor
    Lerobed = 10.0         !Bed-load erosion adaptation length, m
    Ldepbed = 5.0          !Bed-load deposition adaptation length, m     
    
    !--- Sediment properties ---------------------------------------------------
    poros = 0.4            !Sediment porosity (constant for all sizes)
    rhosed = 2650.0        !Sediment density (constant for all sizes)    
    sedshape = 0.7         !Corey sediment shape factor (constant for all sizes)

    !--- Limiters ----------------------------------------------------------
    Cteqmax = 200.0        !Maximum mass concentration [kg/m^3]
    dzbmax = 0.5           !Maximum bed change per time step [m]
    !Cteqmax = (1.0-poros)*rhosed !Maximum mass concentration   !From Chris' code
    
    !--- Timing -----------------------------------------------
    tStartMorph = 0.0     !Morphology change starting time [hrs]
    tStartBedComp = 0.0   !Bed composition change starting time [hrs]
        
    !--- Solver  -------------------
    tolCtk = 1.0e-8        !Tolerance for total-load concentration 
    tolpbk = 1.0e-8        !Tolerance for bed composition 
    maxitersed = -1        !Maximum number of iteration of outer loop for sedimetns
    itermaxzb = 10         !Maximum number of iterations for bed change loop
    
    !--- Boundaries  ------------------------------------------------------------------------
    isedinflowbc = 1       !Sediment River Boundary condition, 1-Capacity, 2-Specified Influx
    facQtotin = 1.0        !Loading factor for inflow transport rate 
    nsedflux = 0           !Sediment flux boundaries
    nsedsource = 0         !Sediment sources
    nsedbc = 0             !nsedflux+nsedsource    
    sedfluxfile = ''
    sedsourcefile = ''
    
    !--- Scaling factors --------------------------------
    scalebed = 1.0           !Bed-load scaling factor
    scalesus = 1.0           !Suspended-load scaling factor
    scalemorph = 1.0         !Morphology scaling factor    
    scalemorph_orig = 1.0    !Original assigned value of Morphology scaling factor
    scalemorph_ramp = 1.0    !Time-dependent ramp value to apply
    scalemorph_rampdur = 0.0 !Ramp duration for morph scaling factor
    
    !--- Total-load (Beta) correction factor ----------------------------------------------------------
    ibt = 1                !0-constant btk, 1-exponential conc profile, 2-rouse conc profile
    betatot = 0.3          !Total-load correction factor, constant in time and for all sediment sizes   

    !--- Concentration profile --------------------------------------------
    iconcprof = 1  !1-Exponential, 2-Rouse, 3-Lund-CIRP, 4-Van Rijn
    
    !--- Hardbottom --------------------------------------------------
    hardbottom = .false.   !Calculate hard-bottom     
    nhard = 0
    
    !--- Avalanching ----------------------------------------------------
    do_aval = .false.      !Calculate avalanching
    a_repose = 32.0        !Angle of repose, degrees
    nmaxaval = 200         !Iterations for avalanching at every time step
    relax_aval = 0.1       !Relaxation coefficient for avalanching

    !--- Size classes ---------------------------------
    nsed = 0               !Number of sediment sizes classes
    
    !--- Bed Layering ----------------------------------------------------------------------------
    nlay = 0               !Max number of bed-material layers, at least one is needed for single-size
    nlayinp = 0            !Input number of bed-material layers
    mixlayconst = .false.  !Use a constant mixing layer thickness
    dmconst = 0.01         !Constant mixing layer thickness [m]
    dbmin = 0.05           !Minimum bed layer thickness [m]
    db1min = 1.0e-4        !Minimum active layer thickness [m]
    dbmax = 0.5            !Maximum bed layer thickness [m]

    !--- Wave-induced sediment transport -----------------------------
    !wavesedtrans = .true. !bdj
    wavesedtrans = .false.
    scaleoffshore = 1.0
    scaleonshore = 1.0
    scaleroller = 1.0
    scalewaveasym = 1.0
    scaleundertow = 1.0
    
    !--- Percentile diameter -----------
    nperinp = 0
    OutPerDiam = .false.
    ipd = 0
    do i=1,nperdiam
      ipd(iper(i)) = i
    enddo
    
    !--- Erosion of dry cells -------------------
    erosdry%calc = .false.
    erosdry%fac = 0.5  !fraction of erosion to move to neighboring dry cells
    erosdry%slopemin = 0.2
    erosdry%slopemax = 0.5
    
    write_smorph_ramp = .false.
    smorph_file = 'morfac_ramp.txt'
    inquire(file=smorph_file,exist=found)
    if (found) then 
      open (9,file=smorph_file)
      close (9,status='delete')
    endif

    return
    endsubroutine sed_default
    
!************************************************************
    subroutine sed_cards(cardname,foundcard)
! Reads the sediment transport cards from the control file
! Author: Alex Sanchez, USACE-CHL
!************************************************************
#include "CMS_cpp.h"
    use geo_def, only: grdfile
    use comvarbl, only: flowpath,rampdur,nswp,nswp0
    use math_lib, only: sortup
    use diag_lib
    use sed_def
    implicit none
    integer :: i,j,ks,ipr,ndiamlim,nn,ierr
    character(len=37) :: cardname,cdum    
    character(len=200) :: afile,apath
    character(len=120), allocatable :: temppath(:)
    real(ikind),allocatable :: temp(:)
    logical :: foundcard
    character(len=100) :: msg2
    
    foundcard = .true.
    selectcase(cardname)              
    !----- On Switch -----------------
    case('CALC_SEDIMENT_TRANSPORT')
      call card_boolean(77,sedtrans,ierr)
      
    !--- Model --------------------  
    case('SED_TRAN_FORMULATION','SED_TRANS_FORMULATION')
      backspace(77)
      read(77,*) cardname, cdum  
      do i=1,size(asedmodel)
        cardname=asedmodel(i)
        if(cdum(1:3)==cardname(1:3))then
          isedmodel = i
          !write(*,*)'bdj in sediment.F90 and isedmodel = ',isedmodel                                                                                                                                                                                                               
          exit
        endif
      enddo          
      if(isedmodel==4) icapac = 3
      if(isedmodel==5) icapac = 2
      isedmodel = min(isedmodel,3) !3,4, and 5 are the same model
      if(isedmodel==2)then
        call diag_print_error('A-D model is currently not supported in the inline model')     
      endif
!        if(cdum(1:3)/='NET')then
!          write(*,*) 'ERROR: Only the NET model is currently supported in the inline model'
!          open(dgunit,file=dgfile,access='append') 
!          write(dgunit,*) 'ERROR: Only the NET model is currently supported in the inline model'
!          close(dgunit)
!          write(*,*) 'Press any key to continue.'
!          read(*,*)
!          stop
!        endif
      
    !--- Time steps ----------------------
    case('SED_TRAN_CALC_INTERVAL')
      !DO NOTHING
      
    case('MORPH_UPDATE_INTERVAL')
      !DO NOTHING      
      
    !---- Hardbottom --------------------  
    case('HARDBOTTOM_DATASET')
      call card_dataset(77,grdfile,flowpath,hbfile,hbpath,1) 
      if(hbpath(1:4)/='NONE')then 
        hardbottom = .true.
      else
        hardbottom = .false.
      endif        
    
    !---- Scaling Factors ----------------
    case('BED_LOAD_SCALE_FACTOR')
      backspace(77)
      read(77,*) cardname, scalebed

    case('SUSP_LOAD_SCALE_FACTOR')
      backspace(77)
      read(77,*) cardname, scalesus

    case('A_COEFFICIENT_WATANABE')
      backspace(77)
      read(77,*) cardname, Awatan  

    case('CSHORE_EFFB','CSHORE_BREAK_EFFICIENCY')  !added 6/7/2019 bdj
      backspace(77)
      read(77,*) cardname, CSeffb

    case('CSHORE_BLP','CSHORE_BED_LOAD')           !added 6/7/2019 bdj
      backspace(77)
      read(77,*) cardname, CSblp

    case('CSHORE_SLP','CSHORE_SUSP_LOAD')          !added 6/7/2019 bdj
      backspace(77)
      read(77,*) cardname, CSslp
      
    case('MORPH_ACCEL_FACTOR')
      backspace(77)
      read(77,*) cardname, scalemorph 
      if(scalemorph<0.0)then
         call diag_print_error('Invalid value for morphologic acceleration factor',&
            '  Must be larger than 1.0 or equal to 0.0')
      elseif(scalemorph>1.0e-6 .and. scalemorph<(1.0-1.0e-6))then
        call diag_print_error('Invalid value for morphologic acceleration factor',&
            '  Must be larger than 1.0 or equal to 0.0')  
      elseif(scalemorph>30)then
        call diag_print_warning('Morphologic acceleration factor set above',&
          '   the recommended upper limit of 30')
      endif
      scalemorph = max(scalemorph,0.0)      
      scalemorph_orig = scalemorph
      
!meb added 3/11/2019
!-------------------------
      
    case('MORPH_ACCEL_RAMP_DURATION')
      !This sets the amount of time when the morphology acceleration faction starts increasing from one to the set value.
      call card_scalar(77,'days','hrs',scalemorph_rampdur,ierr)

    case('WRITE_ACCEL_RAMP_INFO')
      call card_boolean(77,write_smorph_ramp,ierr)  
!-------------------------
      
    !---- Transport Formula ---------------------
    case('NET_TRANSPORT_CAPACITY','TRANSPORT_FORMULA')
      backspace(77)
      read(77,*) cardname, cdum           
      do i=1,size(acapac)
        cardname=acapac(i)
        if(cdum(1:2)==cardname(1:2))then
          icapac = i
          !write(*,*)'bdj in sedimnet.F90, cdum, icapac = ',cdum,icapac                                                                       
          exit
        endif
      enddo
      
    case('CONCENTRATION_CAPACITY_MAX','CONCENTRATION_CAPACITY_MAX_TOTAL','CAPACITY_MAX_TOTAL')
      backspace(77)
      read(77,*) cardname, Cteqmax
      
    case('BED_CHANGE_MAX','BED_CHANGE_MAXIMUM','MAXIMUM_BED_CHANGE','MAX_BED_CHANGE')
      backspace(77)
      read(77,*) cardname, dzbmax  
      dzbmax = max(dzbmax,0.0)
    
    !---- Sediment Properties ------------------  
    case('SEDIMENT_DENSITY')
      call card_scalar(77,'kg/m^3','kg/m^3',rhosed,ierr)
        
    case('SEDIMENT_POROSITY')
      backspace(77)
      read(77,*) cardname, poros
      if(poros<0.3)then
        write(msg2,*) '  Porosity: ',poros
        call diag_print_warning('Extremely low sediment porosity specified',msg2,&
         '  Consider using a larger value',&
         '  Note: The sediment porosity is not a calibration parameter')  
      elseif(poros>0.46)then
        write(msg2,*) '  Porosity: ',poros
        call diag_print_warning('Extremely high sediment porosity specified',msg2,&
         '  Consider using a larger value',&
         '  Note: The sediment porosity is not a calibration parameter')      
      endif
      
    case('SEDIMENT_COREY_SHAPE_FACTOR')  
      backspace(77)
      read(77,*) cardname, sedshape      
      if(sedshape<0.3)then
        write(msg2,*) '  Corey Shape Factor: ',sedshape  
        call diag_print_warning('Relatively small Corey shape factor specified',msg2)
      elseif(sedshape>1.0)then
        call diag_print_warning('Corey shape factor cannot be greater than 1.0',&
          '  Setting to 1.0')
      endif
      sedshape = max(min(sedshape,1.0),0.1) !Cannot better greater than 1.0
      
    !--- Total Load correction factor -------------   
    case('TOTAL_LOAD_CORR_FACTOR_CONSTANT','TOTAL_LOAD_CORRECTION_FACTOR_CONSTANT')  
      backspace(77)
      read(77,*) cardname, betatot
      if(betatot>1.0)then
        call diag_print_warning('Total-load correction factor set larger than 1.0',&
           '  Limiting to 1.0')
        betatot = 1.0
      elseif(betatot<0.01)then
        call diag_print_warning('Total-load correction factor set to unrealistically ',&
          '  small value lower than 0.01','  Limiting to 0.01')
        betatot = 0.01
      endif
      ibt = 0
      
    case('CONCENTRATION_PROFILE')
      backspace(77)
      read(77,*) cardname, cdum   
      iconcprof = 1
      do i=1,size(aconcprof)
        cardname=aconcprof(i)
        if(cdum(1:3)==cardname(1:3))then            
          iconcprof = i
          exit
        endif
      enddo      
      if(iconcprof==1 .or. iconcprof==3)then
        ibt = 1 !Exponential
      else
        ibt = 2 !Rouse
      endif    
          
    !----- Avalanching ---------------- 
    case('USE_AVALANCHING')
      call card_boolean(77,do_aval,ierr)
        
    case('REPOSE_ANGLE','AVALANCHE_CRITICAL_BEDSLOPE')
      call card_scalar(77,'deg','deg',a_repose,ierr)
      if(a_repose+1.0e-5>=90.0)then
        do_aval = .false.
      elseif(a_repose<=5.0)then
        write(msg2,*)   'Repose angle: ', a_repose
        call diag_print_warning('Small input avalanching bed-slope angle:',msg2)
      else            
        do_aval = .true.
      endif        
        
    case('AVALANCHE_MAXIMUM_ITERATIONS','AVALANCHE_MAX_ITERATIONS','AVALANCHING_MAX_ITERATIONS')
      backspace(77)
      read(77,*) cardname, nmaxaval
    
    case('AVALANCHE_RELAX','AVALANCHE_RELAX_COEF','AVALANCHE_RELAX_COEFFICIENT','AVALANCHE_RELAXATION_COEFFICIENT')
      backspace(77)
      read(77,*) cardname, relax_aval     
     
    !---- Bedslope ---------------------------------    
    case('SLOPE_COEFFICIENT','BEDSLOPE_COEFFICIENT','BEDSLOPE_DIFFUSION_COEF','BEDSLOPE_DIFFUSION_COEFFICIENT')
      backspace(77)
      read(77,*) cardname, dcoeff
      if(dcoeff<0.0001) do_bedslope = .false.        
         
    !case('USE_BEDSLOPE')  !Outdated
    !  call card_boolean(77,do_bedslope,ierr)
    
    case('BEDSLOPE_TRANSPORT_FORMULATION','BEDSLOPE_CAPACITY_FORMULATION')
      backspace(77)  
      read(77,*) cardname, cdum
      if(cdum(1:3)=='DEY')then
        ibedslope = 1
      elseif(cdum(1:2)=='BA')then
        ibedslope = 2
      elseif(cdum(1:2)=='WU')then
        ibedslope = 3
      endif 
      
    case('BEDSLOPE_BEDLOAD_COEFFICIENT')
      backspace(77)
      read(77,*) cardname, betaslope
      betaslope = min(max(betaslope,0.0),4.0) !Limit value, note default = 1.6 
    
    case('BEDSLOPE_SUSPENDED_COEFFICIENT')
      backspace(77)
      read(77,*) cardname, effslope
      effslope = min(max(effslope,0.0),0.06) !Limit value, note default = 0.02
    
    !---- Fall Velocity ---------------------------        
    case('SEDIMENT_FALL_VELOCITY_FORMULA','SEDIMENT_FALL_VEL_FORM','FALL_VELOCITY_FORMULA')
      backspace(77)
      read(77,*) cardname, cdum !SOULSBY | WU-WANG
      if(cdum(1:2)=='SO')then				 !Previously, this was mishandling the assignment of type
        sedclass(:)%iws = 2
      elseif(cdum(1:2)=='WU')then
        sedclass(:)%iws = 3
      elseif(cdum(1:2)=='MA')then
        sedclass(:)%iws = 1
      endif   
      
    case('SEDIMENT_FALL_VELOCITY') !Only for single-size sediment transport
      if(nsed>1)then
        call diag_print_error('Use SEDIMENT_FALL_VELOCITIES for multiple sediment sizes')
      endif
      nlay = 1
      nsed = 1
      call sedclass_resize
      call card_scalar(77,'m/s','m/s',sedclass(1)%wsfall,ierr)
      sedclass(1)%iws = 1 !User-specified
        
    case('SEDIMENT_FALL_VELOCITIES')
      backspace(77)
      read(77,*) cardname, nsed
      call sedclass_resize
      backspace(77)
      read(77,*) cardname, nsed, (sedclass(ks)%wsfall,ks=1,nsed)    
      do ks=1,nsed
        sedclass(ks)%iws = 1 !User-Specified
      enddo
      
    !---- Mixing -----------------------------------  
    case('SCHMIDT_NUMBER')
      backspace(77)
      read(77,*) cardname, schmidt       
      sedconstmix = .false.      
        
    case('SEDIMENT_CONSTANT_MIX_COEFF','SEDIMENT_CONSTANT_MIX_COEFFICIENT')
      backspace(77)
      read(77,*) cardname, cmixsed    
      sedconstmix = .true.
      
    case('BEDLOAD_DISPERSION_PARAMETER')   
      backspace(77)
      read(77,*) cardname, cmixbedload
      
    !---- Boundaries --------------------------- 
    case('SEDIMENT_INFLOW_BC')
      backspace(77)
      read(77,*) cardname, cdum
      selectcase(cdum)
      case('CLEAR_WATER')
        isedinflowbc = 1
        facQtotin = 0.0
      case('INFLOW_RATE') 
        isedinflowbc = 2
      case('CAPACITY')
        isedinflowbc = 1  
      endselect  
      
!    case('SEDIMENT_CONC_CELLSTRING')             
      
    case('SEDIMENT_INFLOW_TRANSPORT_RATE')
      backspace(77)
      read(77,*) cardname, Qtotin
        
    case('SEDIMENT_INFLOW_LOADING_FACTOR')  
      backspace(77)
      read(77,*) cardname, facQtotin
      
    !----- Wave-Induced Sediment Transport ------------------    
    case('WAVE_INDUCED_SED_TRANS','WAVE-INDUCED_SED_TRANS','WAVESEDTRANS')
      call card_boolean(77,wavesedtrans,ierr)       
        
    case('ONSHORE_SED_TRANS_SCALE_FACTOR')
      backspace(77)
      read(77,*) cardname, scaleonshore
        
    case('OFFSHORE_SED_TRANS_SCALE_FACTOR')
      backspace(77)
      read(77,*) cardname, scaleoffshore   
      
    case('ROLLER_SED_TRANS_SCALE_FACTOR')
      backspace(77)
      read(77,*) cardname, scaleroller
        
    case('UNDERTOW_SED_TRANS_SCALE_FACTOR')
      backspace(77)
      read(77,*) cardname, scaleundertow
      
    case('WAVE_ASYM_SED_TRANS_SCALE_FACTOR')
      backspace(77)
      read(77,*) cardname, scalewaveasym          
      
    !------- Start of Morphology and Bed Composition Calculation -----------------      
    case('CALC_MORPH_DURING_RAMP','MORPH_DURING_RAMP_CALC')
      call card_boolean(77,calcmorph,ierr)
      if(.not.calcmorph) tStartMorph = RampDur
        
    case('MORPH_START_TIME','BED_CHANGE_START_TIME')
      call card_scalar(77,'days','hours',tStartMorph,ierr)
      tStartMorph = max(tStartMorph,0.0)
      
    case('CALC_BED_COMPOSITION','BED_COMPOSITION_CALC')
      call card_boolean(77,calcBedComp,ierr)
      
    case('BED_COMPOSITION_START_TIME','BED_COMP_START_TIME')
      call card_scalar(77,'days','hours',tStartBedComp,ierr)
      tStartBedComp = max(tStartBedComp,0.0)
      
    !---- Sediment Budget -------------------------
    case('GLOBAL_SEDIMENT_BUDGET','GLOBAL_SEDIMENT_BALANCE')
      call card_boolean(77,sedbalance,ierr)
      
    !---- Hiding and Exposurre --------------------------      
    case('HIDING_EXPOSURE_FORMULA','HIDING_EXPOSURE_FORMULATION')
      backspace(77)
      read(77,*) cardname,cdum
      if(cdum(1:2)=='EG')then !Egiazaroff
        ihidexpform = 1
      elseif(cdum(1:2)=='PA')then !Parker
        ihidexpform = 2
      elseif(cdum(1:2)=='WU')then !Wu
        ihidexpform = 3
      else
        call diag_print_warning('Invalid Sediment Hiding and Exposure Formula: ',cdum)  
        ihidexpform = 0 
      endif
        
    case('HIDING_EXPOSURE_COEFFICIENT','HIDING_EXPOSURE_COEFF','HIDING_EXPOSURE_COEF')
      backspace(77)  
      read(77,*) cardname, mhe
      mhe = max(min(mhe,10.0),0.0) !Limit values
      
    case('HIDING_EXPOSURE_MINIMUM','HIDING_EXPOSURE_CORRECTION_MINIMUM','HID_EXP_COR_MIN')  
      backspace(77)  
      read(77,*) cardname, varsigmamin
    
    case('HIDING_EXPOSURE_MAXIMUM','HIDING_EXPOSURE_CORRECTION_MAXIMUM','HID_EXP_COR_MAX')  
      call card_scalar(77,' ',' ',varsigmamax,ierr)
      
    !--- Sediment Size Class Block ---------------  
    case('SEDIMENT_SIZE_CLASS_BEGIN')
      call sedclass_block
      
    !--- Bed Layer Block ---------------  
    case('BED_LAYER_BEGIN')
      call bedlay_block
        
    !---- Since Grain Size -----------------------------   
    case('CONSTANT_GRAIN_SIZE')
      call card_scalar(77,'mm','mm',singleD50,ierr)
      nsed = 1; nlay = 1 
      call sedclass_resize
      call bedlay_resize
      sedclass(1)%diam = singleD50
      singlesize = .true.
      variableD50 = .false.
      constd50 = .true.
      calcd50 = .false.
            
    case('SEDIMENT_GRAIN_SIZE')  !For old input files
      call card_scalar(77,'mm','mm',singleD50,ierr)
      nsed = 1; nlay = 1
      call sedclass_resize
      call bedlay_resize
      sedclass(1)%diam = singleD50
      singlesize = .true.
      variableD50 = .false.            
           
    case('TRANSPORT_GRAIN_SIZE')
      call card_scalar(77,'mm','mm',singleD50,ierr)
      nsed = 1; nlay = 1
      call sedclass_resize
      call bedlay_resize
      sedclass(1)%diam = singleD50
      singlesize = .true.
      variableD50 = .true.              !Shouldn't this be False?  MEB 9/30/19
      transd50 = .true.
      calcd50 = .false.
      
    !---- Multiple Grain sizes ------------------      
    case('BED_COMPOSITION_INPUT') !used for surface layer only here
      backspace(77)
      read(77,*) cardname, cdum
      call bedlay_resize
      do i=1,size(apbkinp)
        cardname = apbkinp(i)
        if(cdum(1:3)==cardname(1:3))then
          bedlay(1)%ipbkinp = i
          exit
        endif
      enddo             
      !if(bedlay(1)%ipbkinp<=3 .or. bedlay(1)%ipbkinp==6)then              
      !  pbinitconst = .true. !composition uniform with depth
      !endif        
      singlesize = .false.
      bedlay(:)%ipbkinp = bedlay(1)%ipbkinp
      bedlay(:)%inppbk = .true.
      
    case('BED_FRACTIONAL_COMPOSITION_DATASET')
      backspace(77)
      read(77,*) cardname, pbkfile, pbkpath 
      call bedlay_resize
      bedlay(:)%ipbkinp = 5
      bedlay(:)%inppbk = .true.
      singlesize = .false.            
        
    case('SEDIMENT_STANDARD_DEVIATION')  
      call bedlay_resize   
      call card_scalar(77,'mm','mm',bedlay(1)%geostddev,ierr)
      bedlay(1)%geostddev = max(bedlay(1)%geostddev,1.1)
      bedlay(1)%geostddev = min(bedlay(1)%geostddev,10.0)
      bedlay(:)%geostddev = bedlay(1)%geostddev
      bedlay(:)%ipbkinp = 1 !D50_SIGMA
      bedlay(:)%inppbk = .true.
        
    case('SEDIMENT_SIZE_CLASS_NUMBER')  
      backspace(77)
      read(77,*) cardname, nsed 
      nsed = max(1,nsed)
        
    case('MULTIPLE_GRAIN_SIZES')
      backspace(77)
      read(77,*) cardname, nsed
      backspace(77)
      call sedclass_resize
      allocate(temp(nsed))
      read(77,*) cardname, nsed, (temp(ks),ks=1,nsed)
      call sortup(nsed,temp) !Make sure grain sizes are in ascending order
      sedclass(:)%diam = temp(:)
      sedclass(:)%idiam = 1 !User-specified diameter (Missed line? hli, 04/22/16)
      deallocate(temp)
      singlesize = .false.
        
    case('SEDIMENT_SIZE_CLASS_DIAMETERS')
      backspace(77)
      read(77,*) cardname, i
      backspace(77)
      if(i<nsed)then
        call diag_print_error('Number of sediment classes is defined incorrectly',&
          '  Check card: SEDIMENT_SIZE_CLASS_DIAMETER')
      endif
      nsed = max(nsed,i)
      call sedclass_resize
      allocate(temp(nsed))
      read(77,*) cardname, nsed, (temp(ks),ks=1,nsed)
      call sortup(nsed,temp) !Make sure grain sizes are in ascending order
      sedclass(:)%diam = temp(:)
      sedclass(:)%idiam = 1 !User-specified diameter
      deallocate(temp)
      singlesize = .false.    
        
    case('SEDIMENT_SIZE_CLASS_LIMITS')
      backspace(77)
      read(77,*) cardname, ndiamlim
      backspace(77)
      nsed = ndiamlim-1
      call sedclass_resize
      allocate(temp(ndiamlim))
      read(77,*) cardname, ndiamlim, (temp(i),i=1,ndiamlim)
      call sortup(ndiamlim,temp) !Make sure grain sizes are in ascending order
      do ks=1,nsed
        sedclass(ks)%diamlim(1:2) = temp(ks:ks+1)
      enddo
      deallocate(temp)
      singlesize = .false.
        
    case('SEDIMENT_SIZE_CLASS_FRACTIONS')
      backspace(77)
      read(77,*) cardname, i
      if(i<nsed)then
        call diag_print_error('Number of sediment classes is defined incorrectly',&
          '  Check card: SEDIMENT_SIZE_CLASS_DIAMETER')
      endif
      nsed = max(nsed,i)
      backspace(77)
      read(77,*) cardname, nsed, (bedlay(1)%pbconst(ks),ks=1,nsed)
      bedlay(:)%ipbkinp = 4
      bedlay(:)%inppbk = .true.
      
    case('BED_LAYER_THICKNESS_INPUT')  
      backspace(77)
      read(77,*) cardname, cdum
      do i=1,size(adbinp)
        if(cdum==adbinp(i))then
          bedlay(:)%idbinp = i
          exit
        endif
      enddo                
      singlesize = .false.
      
    case('BED_LAYER_THICKNESS_DATASET')
      call bedlay_resize  
      backspace(77)
      read(77,*) cardname, bedlay(1)%dbfile, bedlay(1)%dbpath 
      bedlay(:)%idbinp = 4
      bedlay(:)%dbfile = bedlay(1)%dbfile
      bedlay(:)%dbpath = bedlay(1)%dbpath 
      singlesize = .false.
      
    case('BED_LAYERS_MAX_NUMBER','BED_LAYER_MAX_NUMBER')
      backspace(77)
      read(77,*) cardname, nlay   
      call bedlay_resize
        
    case('BED_LAYERS_NUMBER','NUMBER_BED_LAYERS')
      backspace(77)
      read(77,*) cardname, nlayinp
      
    case('BED_LAYERS_CONSTANT_THICKNESS','BED_LAYER_CONSTANT_THICKNESS',&
         'BED_LAYERS_THICKNESS_CONSTANT','BED_LAYER_THICKNESS_CONSTANT')
      call bedlay_resize  
      backspace(77)
      read(77,*) cardname, bedlay(1)%dbconst
      bedlay(:)%dbconst = bedlay(1)%dbconst
      bedlay(:)%idbinp = 1
        
    case('BED_LAYERS_THICKNESS','BED_LAYER_THICKNESSES',&
         'BED_LAYER_THICKNESS','BED_LAYER_THICKNESS_LIST')
      backspace(77)
      read(77,*) cardname, nlayinp
      nlay=max(nlayinp,nlay)
      call bedlay_resize
      backspace(77)
      read(77,*) cardname, nlayinp, (bedlay(j)%dbconst,j=1,nlayinp)
      if(nlayinp<nlay) bedlay(nlayinp+1:nlay)%dbconst = bedlay(nlayinp)%dbconst !Fill-in unspecified values
      bedlay(nlayinp+1:nlay)%idbinp = 1      
      
    case('MIXING_LAYER_FORMULATION','ACTIVE_LAYER_FORMULATION')
      backspace(77)
      read(77,*) cardname, cdum
      if(cdum(1:3)=='CON')then
        mixlayconst = .true. !CONSTANT
      else
        mixlayconst = .false. !AUTOMATIC | VARIABLE
      endif
    
    case('MIXING_LAYER_CONSTANT_THICKNESS','MIXING_LAYER_THICKNESS',&
      'MIXING_LAYER_THICKNESS_CONSTANT',&
      'ACTIVE_LAYER_CONSTANT_THICKNESS','ACTIVE_LAYER_THICKNESS_CONSTANT')
      call card_scalar(77,'m','m',dmconst,ierr)
      mixlayconst = .true.          
        
    case('BED_LAYERS_MIN_THICKNESS','BED_LAYER_MIN_THICKNESS','INACTIVE_LAYER_MIN_THICKNESS') !minimum thickness is the same for mixing layer and all other layers  
      call card_scalar(77,'m','m',dbmin,ierr)
      
    case('ACTIVE_LAYER_MIN_THICKNESS') !minimum thickness is the same for mixing layer and all other layers  
      call card_scalar(77,'m','m',db1min,ierr)  
      
    case('BED_LAYERS_MAX_THICKNESS','BED_LAYER_MAX_THICKNESS') !minimum thickness is the same for mixing layer and all other layers
      call card_scalar(77,'m','m',dbmax,ierr)
        
    !---- Adaptation parameters ------------
    !Total load adaptation
    case('ADAPTATION_METHOD_TOTAL')
      backspace(77)
      read(77,*) cardname, cdum
      do i=1,size(atotm)
        if(cdum==atotm(i))then
          iadapttot = i
          exit
        endif
      enddo 
      
    case('ADAPTATION_LENGTH_TOTAL') !Alex, 09/04/2008 
      call card_scalar(77,'m','m',Ltot,ierr)
      iadapttot = 1
      if(abs(Ltot+1.0)<0.1)then
        iadapttot = 5 !Max of bed and suspended adaptation lengths
      elseif(Ltot<-1.9)then 
        iadapttot = 6 !Weighted average
      elseif(Ltot<0.1)then 
        iadapttot = 1
        Ltot = max(Ltot,0.1)
      endif
        
    case('ADAPTATION_TIME_TOTAL')
      call card_scalar(77,'s','s',Ttot,ierr)
      iadapttot = 2  
        
    case('ADAPTATION_EROS_LENGTH_TOTAL','ADAPTATION_LENGTH_TOTAL_EROSION')
      call card_scalar(77,'m','m',Lerotot,ierr)
        
    case('ADAPTATION_DEPO_LENGTH_TOTAL','ADAPTATION_DEPOS_LENGTH_TOTAL',&
         'ADAPTATION_DEPOSITION_LENGTH_TOTAL','ADAPTATION_LENGTH_TOTAL_DEPOSITION')
      call card_scalar(77,'m','m',Ldeptot,ierr)
        
    case('ADAPTATION_LENGTH_DATASET')  
      backspace(77)
      read(77,*) cardname, aLtotfile, aLtotpath  
      
    !Suspended load adaptation
    case('ADAPTATION_METHOD_SUSPENDED')
      backspace(77)
      read(77,*) cardname, cdum
      do i=1,size(asusm)
        if(cdum==asusm(i))then
          iadaptsus = i
          exit
        endif
      enddo
        
    case('ADAPTATION_COEFF_SUSPENDED','ADAPTATION_COEFF_METHOD_SUSP') !OLD format
      backspace(77)  
      read(77,*) cardname, i
      if(i==2)then
        iadaptsus = 5 !Armanini and diSilvio
      elseif(i==1)then
        iadaptsus = 6 !Lin
      endif
      
    case('ADAPTATION_LENGTH_SUSPENDED')
      call card_scalar(77,'m','m',Lsus,ierr)
      
    case('ADAPTATION_TIME_SUSPENDED')
      call card_scalar(77,'s','s',Tsus,ierr)
        
    case('ADAPTATION_COEFFICIENT_SUSPENDED')
      backspace(77)  
      read(77,*) cardname, alphasus
      alphasus = max(alphasus,0.0)  
      
    case('ADAPTATION_EROS_LENGTH_SUSPENDED')  
      call card_scalar(77,'m','m',Lerosus,ierr)
        
    case('ADAPTATION_DEPO_LENGTH_SUSPENDED','ADAPTATION_DEPOS_LENGTH_SUSPENDED')
      call card_scalar(77,'m','m',Ldepsus,ierr)
      
    !Bed load adaptation
    case('ADAPTATION_METHOD_BED')
      backspace(77)
      read(77,*) cardname, cdum
      do i=1,size(abedm)
        if(cdum==abedm(i))then
          iadaptbed = i
          exit
        endif
      enddo    
        
    case('ADAPTATION_LENGTH_BED')
      call card_scalar(77,'m','m',Lbed,ierr)          
      if(Lbed<0.0 .or. iadaptbed==5)then
        Lbed = abs(Lbed)
        fbed = max(Lbed,0.1) !Depth-dependant
        iadaptbed = 5
      else
        iadaptbed = 1  
      endif
      
    case('ADAPTATION_TIME_BED')
      call card_scalar(77,'s','s',Tbed,ierr)
        
    case('ADAPTATION_DEPTH_FACTOR_BED')
      backspace(77)
      read(77,*) cardname, fbed          
      fbed = max(abs(fbed),0.1) !Depth-dependant
      iadaptbed = 5
              
    case('ADAPTATION_EROS_LENGTH_BED')
      call card_scalar(77,'m','m',Lerobed,ierr)
         
    case('ADAPTATION_DEPO_LENGTH_BED')
      call card_scalar(77,'m','m',Ldepbed,ierr)
      
    !---- Numerical Methods ----------------        
    case('SEDIMENT_MAX_ITERATIONS','SEDIMENT_MAX_ITER') !Outer loop (for implicit temporal scheme)
      backspace(77)
      read(77,*) cardname, maxitersed
      maxitersed = max(maxitersed,5)
      
    case('SEDIMENT_CONCENTRATION_ITERATIONS','SEDIMENT_TRANSPORT_ITERATIONS') !Inner loop (matrix solver) iterations
      backspace(77)
      read(77,*) cardname, nswp(4)
      nswp(4) = max(nswp(4),5)
      nswp0(4)=nswp(4)
        
    case('SEDIMENT_COUPLING')
      call card_boolean(77,sedcouple,ierr)
        
    case('CONC_TOLERANCE') !Tolerance for total-load concentration (hli,05/10/16)
    tolpbk = 1.0e-8        
      backspace(77)
      read(77,*) cardname, tolCtk
      
    case('COMP_TOLERANCE') !Tolerance for bed composition (hli,05/10/16)
      backspace(77)
      read(77,*) cardname, tolpbk
        
!    case('SEDIMENT_INFLOW_BC') 

    !=== Erosion of dry cells ============================
#ifdef DEV_MODE
    case('EROSION_DRY_CELLS','ERODE_DRY_CELLS')
      call card_boolean(77,erosdry%calc,ierr)
      
    case('EROSION_FRACTION_DRY_NEIGHBORS','EROSION_TRANSFER_DRY_CELLS','EROSION_DRY_CELLS_FACTOR')
      backspace(77)
      read(77,*) cardname, erosdry%fac
      if(erosdry%fac>1.0e-4)then
        erosdry%calc = .true.
      else
        erosdry%calc = .false. 
      endif
      
    case('EROSION_DRY_CELLS_MIN_SLOPE','EROSION_DRY_CELLS_SLOPE_MIN')
      backspace(77)
      read(77,*) cardname, erosdry%slopemin
      
    case('EROSION_DRY_CELLS_MAX_SLOPE','EROSION_DRY_CELLS_SLOPE_MAX')
      backspace(77)
      read(77,*) cardname, erosdry%slopemax
#endif

    !=== Cohesive Sediments =================================
    case('COHESIVE_SEDIMENT')
      backspace(77)   
      read(77,*) cardname, ks   
      if(ks.eq.1) then 
        cohesivesed = .true.
        backspace(77)   
        read(77,*) cardname, ks, cohk1, cohk2, cohcp, cohn, cohr,   &
                   cohsalkmax, cohsalcp, cohsaln, cohtaubp, cohturbn1, cohturbn2, cohturbk1, &
                   cohdepmax0,cohdepmin0, coherodm0, coherodcr0,coherodn, &  
                   pcmax,pcmin  
        endif

    case('COHESIVE_SEDIMENT_SIMPLE')   !Less Input Parameters     
      backspace(77)   
      read(77,*) cardname, ks   
      if(ks.eq.1)then 
        cohesivesed = .true.
        backspace(77)   
        read(77,*) cardname, ks, cohk1, cohcp, cohsalkmax, cohtaubp, cohturbk1, &
                    cohdepmax0, coherodm0, coherodcr0  
        cohk2=0.008
        cohn=1.3
        cohr=4.65
        cohsalcp=30.0
        cohsaln=0.5
        cohturbn1=0.165 
        cohturbn2=0.165
        cohdepmin0=0.0
        !coherodn=1.0
        coherodn=2.5
        pcmax=0.6
        pcmin=0.1
      endif

    case('BED_CONSOLIDATION')
      backspace(77)   
      read(77,*) cardname, ks   
      if(ks.eq.1) then 
        consolidation = .true.
        allocate(tconsolid0(nlay))
        backspace(77)   
        read(77,*) cardname, ks, rhobedcoh0, arhobed, prhobed, rhobedcoh1yr, betarhobed, &  
                   rhobednoncoh,(tconsolid0(j), j=1,nlay), &
                   methcoherodcr,coherodcratrho0,coherodcrtau,rhobedcohercr0,coherodcrn
      endif
        
    !=== Boundary Conditions ==================  
    case('SEDIMENT_SOURCE_CELLSTRING')        
      nsedsource = nsedsource + 1
      allocate(temppath(nsedsource-1))
      do i=1,nsedsource-1
         temppath(i) = sedsourcepath(i)
      enddo
      deallocate(sedsourcepath)
      allocate(sedsourcepath(nsedsource))
      do i=1,nsedsource-1        
        sedsourcepath(i) = temppath(i)
      enddo
      deallocate(temppath)     
      backspace(77)   
      read(77,*) cardname, sedsourcefile, sedsourcepath(nsedsource)           
        
    case('SEDIMENT_FLUX_CELLSTRING')        
      nsedflux= nsedflux + 1  
      allocate(temppath(nsedflux-1))
      do i=1,nsedsource-1
         temppath(i) = sedfluxpath(i)
      enddo
      deallocate(sedfluxpath)
      allocate(sedfluxpath(nsedflux))
      do i=1,nsedflux-1        
        sedfluxpath(i) = temppath(i)
      enddo
      deallocate(temppath)   
      backspace(77)
      read(77,*) cardname, sedfluxfile, sedfluxpath(nsedflux)         
        
    case('CUSTOM_DATASET')
      call card_dataset(77,grdfile,flowpath,afile,apath,1)
      nn = len_trim(apath)
      read(apath(nn-1:nn),'(I2)') ipr
      do i=1,nperdiam
        if(ipr==iper(i)) exit
      enddo
      nlay = max(nlay,1) !At least one layer required
      call bedlay_resize
      bedlay(:)%perdiam(i)%file = afile
      bedlay(:)%perdiam(i)%path = apath
      bedlay(:)%perdiam(i)%inp = .true.
            
    case default  
      if(cardname(4:11)=='_DATASET')then                    
        call card_dataset(77,grdfile,flowpath,afile,apath,1)
        read(cardname(2:3),'(I2)') ipr
        do i=1,nperdiam
          if(ipr==iper(i))then
            if(iper(i)==50) variabled50 = .true.  
            exit
          endif
        enddo
        nlay = max(nlay,1) !At least one layer required
        call bedlay_resize
        bedlay(:)%perdiam(i)%file = afile
        bedlay(:)%perdiam(i)%path = apath
        bedlay(:)%perdiam(i)%inp = .true.
      else            
        foundcard = .false.
      endif
        
    endselect               
    
    return
    endsubroutine sed_cards
    
!*************************************************************
    subroutine sedclass_block()
! Reads a sediment size class block from the control file
!
! Change Log: 
!  03-24-14
!    New Feature: Added warning messages for invalid
!               formula specification.
!    Bug Fix: Invalid string check for WU_WANG formulas
!    Bug Fix: Invalid units for critical shields parameter
!    Bug Fix: Invalid units for critical shear stress
!
! Author: Alex Sanchez, USACE-CHL
!*************************************************************
    use geo_def, only: grdfile
    use comvarbl, only: rampdur
    use sed_def
    use diag_lib
    use prec_def
    implicit none
    integer :: ii,ks,ierr,iswap(1),iswap2
    character(len=37) :: cardname,cdum
    logical :: foundcard
    real(ikind),allocatable :: diamtemp(:)
    type(sed_size_class), allocatable :: sedclasstemp(:)
    
    nsed = nsed + 1    
    call sedclass_resize
    
d1: do ii=1,10
      foundcard = .true.
      read(77,*,iostat=ierr) cardname
      if(ierr/=0) exit d1
      if(cardname(1:1)=='!' .or. cardname(1:1)=='#') cycle
      selectcase(cardname)
        case('NAME')
          backspace(77)
          read(77,*) cardname,sedclass(nsed)%name
          
        case('DIAMETER','CHARACTERISTIC_DIAMETER')
          call card_scalar(77,'mm','mm',sedclass(nsed)%diam,ierr)
          sedclass(nsed)%idiam = 1
          
        case('DIAMETER_LIMITS','DIAMETER_BOUNDS')
          backspace(77)
          read(77,*) cardname,sedclass(nsed)%diamlim(1),sedclass(nsed)%diamlim(2) !mm
          sedclass(nsed)%idiam = 2
          
        case('COREY_SHAPE_FACTOR','SHAPE_FACTOR','SHAPE')  
          backspace(77)
          read(77,*) cardname,sedclass(nsed)%shape
          
        case('FALL_VELOCITY')
          backspace(77)  
          read(77,*,iostat=ierr) cardname, cdum
          if(cdum(1:2)=='SO')then
            sedclass(nsed)%iws = 2
          elseif(cdum(1:2)=='WU')then
            sedclass(nsed)%iws = 3
          else
            call card_scalar(77,'m/s','m/s',sedclass(nsed)%wsfall,ierr)  
            sedclass(nsed)%iws = 1 !User-specified  
          endif
          
        case('FALL_VELOCITY_FORMULA')
          backspace(77)  
          read(77,*,iostat=ierr) cardname, cdum
          select case (cdum(1:2))
          case ('SO')
            sedclass(nsed)%iws = 2
          case ('WU')  
            sedclass(nsed)%iws = 3
          case ('MA')
            sedclass(nsed)%iws = 1  
          case default    
            call diag_print_warning('Invalid Sediment Fall Velocity Formula',&
              '  Using SOULSBY')
            sedclass(nsed)%iws = 2
          end select
          
        case('SHIELDS_PARAMETER','SHIELDS','CRITICAL_SHIELDS','CRITICAL_SHIELDS_PARAMETER')
          backspace(77)  
          read(77,*,iostat=ierr) cardname, sedclass(nsed)%thetacr !Dimensionless  
          sedclass(nsed)%icr = 1
        
        case('CRITICAL_SHEAR_FORMULA','CRITICAL_STRESS_FORMULA')
          backspace(77)  
          read(77,*,iostat=ierr) cardname, cdum
          if(cdum(1:2)=='SO')then !Soulsby (1997)
            sedclass(nsed)%icr = 3
          elseif(cdum(1:2)=='WU')then !Wu and Wang (1999)
            sedclass(nsed)%icr = 4
          else
            call diag_print_warning('Invalid Sediment Critical Shear Stress Formula',&
              '  Using SOULSBY')
            sedclass(nsed)%icr = 3
          endif
          
        case('CRITICAL_SHEAR_STRESS','CRITICAL_STRESS','CRITICAL_SHEAR')
          backspace(77)
          read(77,*,iostat=ierr) cardname, cdum
          if(cdum(1:2)=='SO')then
            sedclass(nsed)%icr = 3
          elseif(cdum(1:2)=='WU')then
            sedclass(nsed)%icr = 4
          else
            call card_scalar(77,'Pa','Pa',sedclass(nsed)%taucr,ierr)    
            sedclass(nsed)%icr = 2
          endif
          
        !case('COVERAGE')
        !  call card_dataset(77,grdfile,flowpath,&
        !    sedclass(nsed)%covfile,sedclass(nsed)%covpath)
          
        case('SEDIMENT_SIZE_CLASS_END','END')
          exit d1
          
        case default
          foundcard = .false.
          call diag_print_warning('Invalid Card: ',cardname,&
            '  Found in block: ','    SEDIMENT_SIZE_CLASS_BEGIN')
        
      endselect
    enddo d1
    
    if(sedclass(nsed)%diam<0.0)then
      call diag_print_error('No characteristic grain size specified for sediment size class')
    endif

    if(nsed<2) return
    
    !Make sure characteristic diameters are sorted in ascending order
    allocate(sedclasstemp(1),diamtemp(nsed))
    do ks=1,nsed-1
      diamtemp=sedclass%diam
      iswap=minloc(diamtemp(ks:nsed)) !returns the index value for the minimum value
      iswap2=iswap(1)+ks-1  !global index
      if(iswap2/=ks)then
        sedclasstemp(1)=sedclass(ks)
        sedclass(ks)=sedclass(iswap2)
        sedclass(iswap2)=sedclasstemp(1)
      endif
    enddo   
    deallocate(sedclasstemp,diamtemp)
    
    return
    endsubroutine sedclass_block
    
!**************************************************    
    subroutine sedclass_resize()
! written by Alex Sanchez, USACE-CHL    
!**************************************************    
    use sed_def, only: sedclass,sed_size_class,nsed
    implicit none
    integer :: ks,nsedtemp
    type(sed_size_class),allocatable:: sedclasstemp(:)
    
    if(.not.allocated(sedclass))then
      allocate(sedclass(nsed))
      nsedtemp = 0
    else      
      nsedtemp = size(sedclass)
      if(nsed==nsedtemp) return
      allocate(sedclasstemp(nsedtemp))
      do ks=1,nsedtemp
        sedclasstemp(ks) = sedclass(ks)
      enddo
      deallocate(sedclass)
      allocate(sedclass(nsed))      
      do ks=1,nsedtemp
        sedclass(ks) = sedclasstemp(ks)
      enddo
      deallocate(sedclasstemp)
    endif
    
    !Set defaults
    do ks=nsedtemp+1,nsed
      sedclass(ks)%name = ' '
      sedclass(ks)%idiam   = 0    !Diameter method, 1-characteristic diam, 2-diameter limits
      sedclass(ks)%diam    = -1.0
      sedclass(ks)%diamlim(1:2) = -1.0
      sedclass(ks)%iws     = 0    !Fall velocity method, 0-none, 1-user, 2-Soulsby, 3-Wu and Wang
      sedclass(ks)%wsfall  = -1.0
      sedclass(ks)%icr     = 0    !Incipient motion, 0-none,1-User shields,2-user shear,3-Soulsby,4-Wu and Wang
      sedclass(ks)%thetacr = -1.0 !Shields parameter
      sedclass(ks)%taucr   = -1.0 !Critical shear stress
      sedclass(ks)%shape = 0.7    !Cory shape factor
      !sedclass(ks)%existcoverage = .false.
      !sedclass(ks)%covfile = ''
      !sedclass(ks)%covpath = ''
    enddo
    
    return
    endsubroutine sedclass_resize
    
!**************************************************
    subroutine bedlay_block()
! Reads a bed layer block from the input card file
! written by Alex Sanchez, USACE-CHL
!**************************************************
    use size_def
    use geo_def, only: grdfile
    use comvarbl, only: flowpath
    use sed_def
    use diag_lib
    use prec_def
    
    implicit none
    integer :: i,ii,ks,ierr,jlay,ipr,nsedchk
    real    :: value
    character(len=37) :: cardname,cdum
    character(len=200) :: file,path
    logical :: foundcard
    character(len=100) :: msg2,msg3,msg4

    nlayinp = nlayinp + 1    
    
    !apbkinp(0) = 'NONE'
    !apbkinp(1) = 'D50_SIGMA'            !Spatially variable, uniform in depth
    !apbkinp(2) = 'D16_D50_D84'          !Spatially variable, uniform in depth
    !apbkinp(3) = 'D35_D50_D90'          !Spatially variable, uniform in depth *** NOT WORKING *******
    !apbkinp(4) = 'SIZE_CLASS_FRACTIONS' !Specifies a constant distribution curve, Spatially uniform, and uniform in depth
    !apbkinp(6) = 'PERCENTILES'          !Specifies percentiles, uniform in depth
    
    !adbinp(0) = 'NONE'            !No information on bed layers, use default values
    !adbinp(1) = 'CONSTANT'        !Applies the same layer thickness to all layers, everywhere
    !adbinp(3) = 'LAYER_DATASET'   !Specifies the thickness of each layer and cell
    
    jlay = nlayinp !default value
    if(jlay>nlay)then
      nlay = jlay
      call bedlay_resize
    endif  
    
d1: do ii=1,30
      foundcard = .true.
      read(77,*,iostat=ierr) cardname
      if(ierr/=0) exit d1
      if(cardname(1:1)=='!' .or. cardname(1:1)=='#') cycle
      selectcase(cardname)
        case('BED_LAYER_END','END')
          exit d1
          
        case('BED_LAYER_NUMBER','LAYER_NUMBER','BED_LAYER_ID','LAYER','LAYER_ID')
          backspace(77)
          read(77,*,iostat=ierr) cardname, jlay
          if(jlay>nlay)then
            nlay = jlay  
            call bedlay_resize
          elseif (jlay==0) then
            jlay=1
          endif  
          !if(jlay>nlay)then
          !  write(*,*) 'ERROR: Increase maximum number of bed layers'
          !  write(*,*) ' before specifying the bed layer information'
          !  write(*,*) '   Press any key to continue'
          !  read(*,*)
          !  stop
          !endif          
          
        case('INITIAL_THICKNESS_VALUE','INITIAL_THICKNESS',&
           'THICKNESS_VALUE','THICKNESS','THICKNESS_CONSTANT')
          call card_scalar(77,'m','m',bedlay(jlay)%dbconst,ierr)  
          bedlay(jlay)%idbinp = 1
        
        case('INITIAL_THICKNESS_DATASET','THICKNESS_DATASET')
!          backspace(77)
!          read(77,*) cardname, bedlay(jlay)%dbfile, bedlay(jlay)%dbpath 

          call card_dataset(77,grdfile,flowpath,file,path,1)
          bedlay(jlay)%dbfile = file
          bedlay(jlay)%dbpath = path
          bedlay(jlay)%idbinp = 3
          
        case('BED_COMPOSITION_INPUT','COMPOSITION_INPUT','COMPOSITION')
          backspace(77)
          read(77,*) cardname, cdum
          call bedlay_resize
          do i=1,size(apbkinp)
            cardname = apbkinp(i)
            if(cdum(1:3)==cardname(1:3))then
              bedlay(jlay)%ipbkinp = i
              exit
            endif
          enddo
        
        case('SEDIMENT_STANDARD_DEVIATION','GEOMETRIC_STANDARD_DEVIATION',&
             'GEO_STD_DEV','STANDARD_DEVIATION','SIGMA')
          call card_scalar(77,'mm','mm',bedlay(jlay)%geostddev,ierr)
          bedlay(jlay)%geostddev = max(bedlay(jlay)%geostddev,1.1)
          bedlay(jlay)%geostddev = min(bedlay(jlay)%geostddev,10.0)           
          
        case('SEDIMENT_SIZE_CLASS_FRACTIONS','SIZE_CLASS_FRACTIONS','FRACTIONS')
          backspace(77)
          read(77,*) cardname, nsedchk
          if(nsedchk/=nsed)then
            call diag_print_error('Number of sediment classes is defined incorrectly',&
              'Check card: SEDIMENT_SIZE_CLASS_DIAMETER')
          endif
          backspace(77)
          allocate(bedlay(jlay)%pbconst(nsed))
          read(77,*,iostat=ierr) cardname, nsed, (bedlay(jlay)%pbconst(ks),ks=1,nsed)
          if(ierr/=0)then
            write(msg2,*) ' for bed layer ',jlay
            call diag_print_error('Could not read sediment size class fractions ',msg2)
          endif
          bedlay(jlay)%pbconst(:) = bedlay(jlay)%pbconst(:)/sum(bedlay(jlay)%pbconst(:)) !Normalize
          bedlay(jlay)%ipbkinp = 4
        
        case('FRACTIONS_DATASET','SEDIMENT_FRACTIONS_DATASET')
          call card_dataset(77,grdfile,flowpath,bedlay(jlay)%pbkfile,bedlay(jlay)%pbkpath,1)
          bedlay(jlay)%ipbkinp = 5
          
        case default
          if(cardname(4:11)=='_DATASET')then   
            call card_dataset(77,grdfile,flowpath,file,path,1)
            if(len_trim(file)>0 .and. len_trim(path)>0)then !Make sure card is not empty  
              read(cardname(2:3),'(I2)') ipr
              do i=1,nperdiam
                if(ipr==iper(i)) exit
              enddo
              bedlay(jlay)%perdiam(i)%file = file
              bedlay(jlay)%perdiam(i)%path = path
              if (path == 'Datasets/') then                 !This should handle SMS writing out empty Dxx_DATASET cards.
                bedlay(jlay)%perdiam(i)%inp = .false.
              else
                bedlay(jlay)%perdiam(i)%inp = .true.  
              endif  
            endif
          else            
            foundcard = .false.
            write(msg2,*) '  ',trim(cardname)
            write(msg3,*) 'Found in block:'
            write(msg4,*) '  BED_LAYER_BEGIN'
            call diag_print_warning('WARNING: Invalid Card: ',msg2,msg3,msg4)
          endif
      endselect
    enddo d1

    !Read in bed layer
    if(bedlay(jlay)%inppbk .and. jlay > 1)then
      write(msg2,*) 'Bed layer ',jlay,' already specified'
      call diag_print_error('Problem specify sediment bed layers',msg2)
    endif
    bedlay(jlay)%inppbk = .true.  
    
    !!Check bed composition input
    !if(bedlay(jlay)%ipbkinp/=0) return
    !
    !!Count number of input percentile datasets   
    !nperinp = 0    
    !do i=1,nperdiam
    !  if(bedlay(jlay)%perdiam(i)%inp)then
    !    nperinp = nperinp + 1
    !  endif
    !enddo
    !!Set input composition mode
    !selectcase(nperinp)
    !case(0)
    !   bedlay(jlay)%ipbkinp = 4 !SIZE_CLASS_FRACTIONS
    !case(1)
    !  bedlay(jlay)%ipbkinp = 1 !D50_SIGMA
    !case(2)
    !  write(*,*) 'ERROR: Invalid number of input percentile diameter datasets'
    !  write(*,*) '  for bed layer ',jlay
    !  write(*,*) '  Press any key to continue'
    !  read(*,*)
    !  stop
    !case(3)
    !  if(bedlay(jlay)%perdiam(ipd(16))%inp .and. &
    !     bedlay(jlay)%perdiam(ipd(84))%inp)then
    !    bedlay(jlay)%ipbkinp = 2 !D16_D50_D84  
    !  elseif(bedlay(jlay)%perdiam(ipd(35))%inp .and. &
    !         bedlay(jlay)%perdiam(ipd(90))%inp)then 
    !    bedlay(jlay)%ipbkinp = 3  !D35_D50_D90
    !  else
    !    bedlay(jlay)%ipbkinp = 6  !PERCENTILES 
    !  endif
    !case default
    !  bedlay(jlay)%ipbkinp = 6  !PERCENTILES   
    !endselect
    
    return
    endsubroutine bedlay_block
    
!**************************************************************    
    subroutine bedlay_resize()
! written by Alex Sanchez, USACE-CHL    
!**************************************************************
    use sed_def, only: bedlay,bed_layer,nlay,nsed,nperdiam
    implicit none
    integer :: i,j,nlaytemp
    type(bed_layer), allocatable :: bedlaytemp(:)
            
    if(nlay==0) nlay = 1
    if(.not.allocated(bedlay))then   
      allocate(bedlay(nlay))
      nlaytemp = 0
    else
      nlaytemp = size(bedlay)
      if(nlay==nlaytemp) return
      allocate(bedlaytemp(nlaytemp))
      do j=1,nlaytemp
        if(allocated(bedlay(j)%pbconst))then
          allocate(bedlaytemp(j)%pbconst(nsed))
        endif
        bedlaytemp(j) = bedlay(j)
      enddo  
      deallocate(bedlay)
      allocate(bedlay(nlay))
      do j=1,nlaytemp
        if(allocated(bedlaytemp(j)%pbconst))then
          allocate(bedlay(j)%pbconst(nsed))
        endif
        bedlay(j) = bedlaytemp(j)
      enddo  
      deallocate(bedlaytemp)
    endif
    
    !Set Defaults
    do j=nlaytemp+1,nlay      
      bedlay(j)%idbinp = 0         !Bed layer thickness specification mode       
      bedlay(j)%dbconst = -999.0   !Constant bed layer thickness
      bedlay(j)%dbfile = ''        !Bed layer thickness file
      bedlay(j)%dbpath = ''        !Bed layer thickness path 
      bedlay(j)%inppbk = .false.   !Input bed layer composition
      bedlay(j)%ipbkinp = 0        !Bed composition specified mode
      bedlay(j)%geostddev = -999.0 !Layer geometric standard deviation
      do i=1,nperdiam
        bedlay(j)%perdiam(i)%file = ''  
        bedlay(j)%perdiam(i)%path = ''          
      enddo
    enddo
    
    return
    endsubroutine bedlay_resize
    
!***************************************************************************   
    subroutine sed_init()
! Allocates and initializes the sediment transport variables
!
! written by Alex Sanchez, USACE-CHL; Weiming Wu, NCCHE
!***************************************************************************    
#include "CMS_cpp.h"
    use size_def, only: ncells,ncellsD
    use geo_def, only: zb,zb0,x,y
    use flow_def, only: rhow,viscos,grav,gravinv
    use comvarbl, only: ntsch
    use hot_def, only: coldstart
    use sed_def
    use beta_lib, only: bs_init
    use sed_lib
    use cms_def, only: noptset
    use wave_flowgrid_def    
    use const_def, only: deg2rad,small,eps
#ifdef XMDF_IO
    use in_xmdf_lib, only: readscalh5
#endif   
    use in_lib, only: readscalTxt
    use math_lib, only: avgval
    use diag_lib
    use prec_def
    
    implicit none
    integer :: i,j,jj,ks,ierr,k
    real(ikind), allocatable :: d16lay(:),d35lay(:),d50lay(:),d84lay(:),d90lay(:)
    character(len=200) :: file,path
    character(len=100) :: msg2
    character(len=10) :: aext
    
    !Internal variables
    solid = 1.0-poros
    specgrav = rhosed/rhow !Specific gravity
    s1grav = (specgrav-1.0)*grav  !(Specific gravity - 1)*gravity, internal variable  
    d2dstar = diam2dstar(grav,viscos,specgrav,1.0_ikind)
    
    !Surface percentile diameters
    allocate(d50(ncellsD),d90(ncellsD))  !Only ones required. Others only needed for output
    allocate(d50lay(ncellsD))   !(hli, 03/15/16)
    
    if(nsed>1) singlesize = .false.
    !==== Multiple-sizes sediment transport =======
    if(.not.singlesize)then        
      if(nsed<=1)then
        nsed = 5
      endif  
      if(nlay<5)then
        nlay = 5
        call bedlay_resize
      endif
      
      !--- Layer thickness thresholds ------------------------
      dbmax = max(dbmax,5.0*dbmin) !Maximum layer thickness for new layers. Existing layers are not affected
      if(mixlayconst)then
        dmconst = max(dmconst,db1min)  !for constant mixing layer thickness
      endif
      
      !--- Sediment size class diameters ------------
      allocate(diam(nsed),diamlim(nsed+1))
      !Calculate grain sizes if none specified
      if(sum(sedclass(:)%idiam)==0)then !No diameters specified 
        selectcase(bedlay(1)%ipbkinp)  !use surface layer
        case(1)
          call d50sigma2diamlim(nsed,d50lay,bedlay(1)%geostddev,diamlim) !Calculates grain size distribution from D50 and geometric standard deviation
          call diamlim2diam(nsed,diamlim,diam)         
        case(2)
          !bedlay(j)%geostddev=log(sqrt(maxval(d84lay(1:ncells))/minval(d16lay(1:ncells))))  !Estimate from extremal values  
          bedlay(1)%geostddev=log(sqrt(avgval(ncells,d84lay)/avgval(ncells,d16lay)))  !Estimate from extremal values  
          call d50sigma2diamlim(nsed,d50lay,bedlay(1)%geostddev,diamlim)
          call diamlim2diam(nsed,diamlim,diam)   
        case(3)
          !bedlay(j)%geostddev=log((maxval(d90lay(1:ncells))/minval(d35lay(1:ncells)))**0.61) !Estimate from extremal values 
          bedlay(1)%geostddev=log((avgval(ncells,d90lay)/avgval(ncells,d35lay))**0.61) !Estimate from extremal values 
          call d50sigma2diamlim(nsed,d50lay,bedlay(1)%geostddev,diamlim)
          call diamlim2diam(nsed,diamlim,diam)    
        case default
          call diag_print_error('Cannot determine sediment size class diameters')
        endselect
      endif      
      !Use input mode
      do ks=1,nsed
        selectcase(sedclass(ks)%idiam) 
        case(1)    
          diam(ks) = sedclass(ks)%diam !Characteristic diameter
        case(2)
          diam(ks) = sqrt(sedclass(ks)%diamlim(1)*sedclass(ks)%diamlim(2)) !Diameter limits
        case default
          write(msg2,*) '  Size class: ',ks
          call diag_print_error('Missing grain size for ',msg2)  
        endselect
      enddo
      
      !--- Sediment size class limits -----------------
      do ks=2,nsed
        selectcase(sedclass(ks)%idiam) 
        case(1)
          diamlim(ks)=sqrt(diam(ks)*diam(ks-1))
        case(2)
          diamlim(ks) = sedclass(ks)%diamlim(1)
          diamlim(ks+1) = sedclass(ks)%diamlim(2)
        endselect
      enddo
      selectcase(sedclass(1)%idiam) 
      case(1)
        diamlim(1) = diam(1)*diam(1)/diamlim(2)
      case(2)
        diamlim(1) = sedclass(1)%diamlim(1)
        diamlim(2) = sedclass(1)%diamlim(2)
      endselect
      selectcase(sedclass(nsed)%idiam) 
      case(1)
        diamlim(nsed+1) = diam(nsed)*diam(nsed)/diamlim(nsed)
      case(2)
        diamlim(nsed) = sedclass(nsed)%diamlim(1)
        diamlim(nsed+1) = sedclass(nsed)%diamlim(2)
      endselect
      
      !--- Bedchange for ks sediment size ------------
      allocate(dzbk(ncellsD,nsed))  !Only allocate for multiple grain sizes
      dzbk = 0.0
      
      !--- Bed material composition ------------------------------------
      allocate(pbk(ncellsD,nsed,nlay),pbk1(ncellsD,nsed))   !pbk1 contains only mixing layer to save memory *******  
      allocate(pbkstar(ncellsD,nsed))              
      do j=1,nlay        
        if(.not.bedlay(j)%inppbk) cycle        
        !Guess based on input datasets and parameters if no bed composition mode specified
        if(bedlay(j)%ipbkinp==0)then 
          !Count number of input percentile datasets   
          nperinp = 0    
          do i=1,nperdiam
            if(bedlay(j)%perdiam(i)%inp)then
              nperinp = nperinp + 1
            endif
          enddo
          !Set input composition mode
          selectcase(nperinp)
          case(0)
            if(allocated(bedlay(j)%pbconst))then
              bedlay(j)%ipbkinp = 4 !SIZE_CLASS_FRACTIONS
            else
              bedlay(j)%ipbkinp = 0 !None
            endif  
          case(1)            
            if(.not.bedlay(j)%perdiam(ipd(50))%inp .or. bedlay(j)%geostddev<1.0)then                
              write(msg2,*) '  Bed Layer: ',j  
              call diag_print_error('Could not determine bed layer composition mode for ',msg2)
            endif
            bedlay(j)%ipbkinp = 1 !D50_SIGMA
          case(2)
            write(msg2,*) '  Bed Layer: ',j    
            call diag_print_error('Invalid number of input percentile diameter datasets for ',msg2)
          case(3)
            if(bedlay(j)%perdiam(ipd(16))%inp .and. &
              bedlay(j)%perdiam(ipd(84))%inp)then
              bedlay(j)%ipbkinp = 2 !D16_D50_D84  
            elseif(bedlay(j)%perdiam(ipd(35))%inp .and. &
              bedlay(j)%perdiam(ipd(90))%inp)then 
              bedlay(j)%ipbkinp = 3 !D35_D50_D90
            else
              bedlay(j)%ipbkinp = 6 !PERCENTILES 
            endif
          case default
            bedlay(j)%ipbkinp = 6 !PERCENTILES   
          endselect  
        endif
        
        selectcase(bedlay(j)%ipbkinp)
        case(0) !None
          if(j==1)then
            call diag_print_error('Bed composition for first layer must be specified')
          endif
          
        case(1)  !D50 and sigma
          OutPerDiam(ipd(35)) = .true.; OutPerDiam(ipd(50)) = .true.; OutPerDiam(ipd(90)) = .true.
!          allocate(d50lay(ncellsD))
          file = bedlay(j)%perdiam(ipd(50))%file; path = bedlay(j)%perdiam(ipd(50))%path
          call fileext(trim(file),aext)      
          select case (aext)
          case('h5')
#ifdef XMDF_IO
            call readscalh5(file,path,d50lay,ierr) !Dataset should be in mm
#endif         
          case('txt')
            call readscalTxt(file,d50lay,ierr)
            
          case default
            write(msg2,*)"Unknown file"
            call diag_print_error(msg2)
          end select
          
          if(ierr<0) call dper_read_error_msg(file,path)          
          call bed_d50sigma(nsed,diam,diamlim,d50lay,bedlay(j)%geostddev,pbk(:,:,j)) !Note: assumes units of mm for d50, diam, sedsigma       
          deallocate(d50lay)
          
          !open(unit=4045,file='initialbedcomp.txt')
          !do i=1,ncells
          !    write(4045,*)x(i),y(i),(pbk(i,k,1),k=1,nsed)
          !enddo
          !close(4045)
          !stop
          
        case(2)  !D16, D50, D84
          OutPerDiam(ipd(16)) = .true.; OutPerDiam(ipd(50)) = .true.; OutPerDiam(ipd(84)) = .true.
          if(.not.allocated(d16lay)) allocate(d16lay(ncellsD))
          if(.not.allocated(d84lay)) allocate(d84lay(ncellsD))
          if(.not.allocated(d50lay)) allocate(d50lay(ncellsD))
          file = bedlay(j)%perdiam(ipd(16))%file; path = bedlay(j)%perdiam(ipd(16))%path

          call fileext(trim(file),aext)      
          select case (aext)
          case('h5')
#ifdef XMDF_IO
            call readscalh5(file,path,d16lay,ierr) !Dataset should be in mm    
#endif          
          case('txt')
            call readscalTxt(file,d16lay,ierr)
          end select
          
          if(ierr<0) call dper_read_error_msg(file,path)       
          file = bedlay(j)%perdiam(ipd(50))%file; path = bedlay(j)%perdiam(ipd(50))%path
          
          call fileext(trim(file),aext)      
          select case (aext)
          case('h5')
#ifdef XMDF_IO
            call readscalh5(file,path,d50lay,ierr)   
#endif          
          case('txt')
            call readscalTxt(file,d50lay,ierr)
          end select
          
          if(ierr<0) call dper_read_error_msg(file,path)     
          file = bedlay(j)%perdiam(ipd(84))%file; path = bedlay(j)%perdiam(ipd(84))%path

          call fileext(trim(file),aext)      
          select case (aext)
          case('h5')
#ifdef XMDF_IO
            call readscalh5(file,path,d84lay,ierr)    
#endif         
          case('txt')
            call readscalTxt(file,d84lay,ierr)
          end select
          
          if(ierr<0) call dper_read_error_msg(file,path)          
          call bed_d16d50d84(nsed,diam,diamlim,d16lay,d50lay,d84lay,pbk(:,:,j))
          deallocate(d16lay,d50lay,d84lay)
          
        case(3) !D35, D50, D90
          OutPerDiam(ipd(35)) = .true.; OutPerDiam(ipd(50)) = .true.; OutPerDiam(ipd(90)) = .true. 
          allocate(d35lay(ncellsD),d50lay(ncellsD),d90lay(ncellsD))
          file = bedlay(j)%perdiam(ipd(35))%file; path = bedlay(j)%perdiam(ipd(35))%path

          call fileext(trim(file),aext)      
          select case (aext)
          case('h5')
#ifdef XMDF_IO
          call readscalh5(file,path,d35lay,ierr)  
#endif          
          case('txt')
            call readscalTxt(file,d35lay,ierr)
          end select
          
          if(ierr<0) call dper_read_error_msg(file,path)       
          file = bedlay(j)%perdiam(ipd(50))%file; path = bedlay(j)%perdiam(ipd(50))%path
          
          call fileext(trim(file),aext)      
          select case (aext)
          case('h5')
#ifdef XMDF_IO
          call readscalh5(file,path,d50lay,ierr)   
#endif          
          case('txt')
            call readscalTxt(file,d50lay,ierr)
          end select
          if(ierr<0) call dper_read_error_msg(file,path)   
          file = bedlay(j)%perdiam(ipd(90))%file; path = bedlay(j)%perdiam(ipd(90))%path
          
          call fileext(trim(file),aext)      
          select case (aext)
          case('h5')
#ifdef XMDF_IO
          call readscalh5(file,path,d90lay,ierr)  
#endif         
          case('txt')
            call readscalTxt(file,d90lay,ierr)
          end select
          
          if(ierr<0) call dper_read_error_msg(file,path) 
          call bed_d35d50d90(nsed,diam,diamlim,d35lay,d50lay,d90lay,pbk(:,:,j)) !Calculates grain size distribution from D35,D50,D90
          deallocate(d35lay,d50lay,d90lay)
          
        case(4)  !Histogram
          OutPerDiam(ipd(35)) = .true.; OutPerDiam(ipd(50)) = .true.; OutPerDiam(ipd(90)) = .true.          
          bedlay(j)%pbconst = bedlay(j)%pbconst/sum(bedlay(j)%pbconst(1:nsed)) !Normalize bed fractions
          do ks=1,nsed
            pbk(:,ks,j) = bedlay(j)%pbconst(ks) !Initialize bed material
          enddo
          
        case(5)  !Dataset  (outdated)
          OutPerDiam(ipd(35)) = .true.; OutPerDiam(ipd(50)) = .true.; OutPerDiam(ipd(90)) = .true.  
          call read_pbk(j,bedlay(j)%pbkfile,bedlay(j)%pbkpath)
        
        case(6)  !Percentile Datasets
          call bed_perdiam(j)
          !Output any percentile diameter input for any bed layer
          do i=1,nperdiam
            if(bedlay(j)%perdiam(i)%inp) OutPerDiam(i) = .true.
          enddo        
          
        endselect
      enddo
      
      !--- Copy composition to unspecified layers from above downwards -------
      do j=2,nlay
        if(bedlay(j)%ipbkinp==0)then  
          do jj=j-1,1,-1
            if(bedlay(jj)%ipbkinp/=0)then !Bug fix, changed j to jj (Alex Sanchez 04-23-14)
              pbk(:,:,j) = pbk(:,:,jj)
              bedlay(j)%ipbkinp = bedlay(jj)%ipbkinp
              bedlay(j)%pbkfile = bedlay(jj)%pbkfile
              bedlay(j)%pbkpath = bedlay(jj)%pbkpath
              bedlay(j)%perdiam = bedlay(jj)%perdiam
              exit
            endif         
          enddo
        endif
      enddo
      
      !Convert diameters from mm to m
      diam = diam/1000.0       !convert from mm to m
      diamlim = diamlim/1000.0 !convert from mm to m 
      
      !Set previous time step composition
      pbk1 = pbk(:,:,1)        !Fraction of each sediment in first layer of previous time step  
      
      !Percentile Diamters
      !Used to calculate sediment percentiles
      allocate(logdiamlim(nsed+1))
      logdiamlim=log(diamlim)
      !Calculate necessary sediment percentiles for calculations
      !Others are calculated only as needed for output
      call sedpercentile(50,d50)
      call sedpercentile(90,d90)
      
      !Bed layer thickness
      allocate(db(ncellsD,nlay),db1(ncellsD,nlay))
      do j=1,nlay
        selectcase(bedlay(j)%idbinp)
        case(0) !None specified
          if(j==1)then !If first layer, use maximum thickness
            db(:,j) = dbmax
            bedlay(j)%dbconst = dbmax            
          else !Use thickness from above specified layer
            do jj=j-1,1,-1
              if(bedlay(jj)%idbinp/=0)then
                db(:,j) = db(:,jj)
                bedlay(j)%idbinp = bedlay(jj)%idbinp
                bedlay(j)%dbfile = bedlay(jj)%dbfile
                bedlay(j)%dbpath = bedlay(jj)%dbpath
                bedlay(j)%dbconst = bedlay(jj)%dbconst
                exit
              elseif(jj==1)then !If none found then use maximum thickness
                db(:,j) = dbmax  
                bedlay(j)%dbconst = dbmax
                bedlay(j)%idbinp = 1
              endif
            enddo   
          endif
        case(1) !Constant layer thickness
          db(:,j) = bedlay(j)%dbconst
        case(3) !Thickness datasets 
            
          call fileext(trim(bedlay(j)%dbfile),aext)      
          select case (aext)
          case('h5')
#ifdef XMDF_IO
          call readscalh5(bedlay(j)%dbfile,bedlay(j)%dbpath,db(:,j),ierr)
#endif         
          case('txt')
            call readscalTxt(bedlay(j)%dbfile,db(:,j),ierr)
          end select
          
          if(ierr/=0)then
            call diag_print_error('Problem reading bed layer thickness dataset')
          endif
        endselect
      enddo
      
      db1 = db !Previous iteration bed layer thickness
      !Insert Mixing layer (if possible)
      !Minimum mixing layer (1) thickness
      call mixing_layer !db(i,1)      
      !Recalculate layers
      do i=1,ncells
        !Determine if first layer is thick enough to split
        if(db1(i,1)+1.0e-6>=2.0*dbmin)then !Split first layer
          !First merge last two layers   
          db(i,nlay)=db1(i,nlay)+db1(i,nlay-1)     
          pbk(i,:,nlay)=(db1(i,nlay)*pbk(i,:,nlay)+db1(i,nlay-1)*pbk(i,:,nlay-1))/db(i,nlay)
          !Move index of layers 3 to nlay-1        
          db(i,3:nlay-1)=db1(i,2:nlay-2)
          pbk(i,:,3:nlay-1)=pbk(i,:,2:nlay-2)
          !Recalculate layers 1 and 2
          db(i,2)=db1(i,1)-db(i,1)          
          pbk(i,:,2)=pbk(i,:,1)          
          !Note: Mixing layer composition does not change   
        else !Layer too thin to split, revert back to initial value
          db(i,1)=db1(i,1)  
        endif
      enddo
      db1 = db !Make sure they are the same
      
    else !single size, Note: these are used even in single sized transport
      allocate(pbk(ncellsD,1,1))
      pbk = 1.0
      if(nlay==1)then
        if(len_trim(bedlay(1)%perdiam(ipd(50))%file)>0)then
          variableD50 = .true.  
        endif
      endif
      if(constd50)then
        variableD50 = .false.  
      endif
      if(variableD50)then
        if(.not.allocated(bedlay))then
          call diag_print_warning('Missing D50 dataset for variable D50',&
            '  Using transport grain size as D50')
          variableD50 = .false.            
        else
          if(len_trim(bedlay(1)%perdiam(ipd(50))%file)==0)then
            call diag_print_warning('Missing D50 dataset for variable D50',&
              '  Using transport grain size as D50')  
            variableD50 = .false.
          endif
        endif
      endif
      if(variableD50)then !read d50 dataset
        file = bedlay(1)%perdiam(ipd(50))%file
        path = bedlay(1)%perdiam(ipd(50))%path

        call fileext(trim(file),aext)      
        select case (aext)
        case('h5')
#ifdef XMDF_IO
        call readscalh5(file,path,d50,ierr)
#endif       
        case('txt')
          call readscalTxt(file,d50,ierr)
        end select
        
        if(ierr<0) call dper_read_error_msg(file,path)        
        d50 = d50/1000.0
      endif      
      if(.not.allocated(sedclass))then  !No grain size specified. Determine from d50
        calcd50 = .true.
        !variableD50 = .false.
        nsed = 1; nlay = 1
        call sedclass_resize
        call bedlay_resize
        !Calculate average d50
        sedclass(1)%diam=0.0
        do i=1,ncells
         sedclass(1)%diam=sedclass(1)%diam+d50(i)
         !if(abs(d50(i)-d50(1))>1.0e-6) variableD50 = .true. !Test for variable d50
        enddo
        sedclass(1)%diam=1000.0*sedclass(1)%diam/real(ncells,kind=ikind) !Average
        !if(.not.variableD50) sedclass(1)%diam=d50(1) !constant grain size
      endif
      allocate(diam(nsed))
      diam = sedclass(1)%diam/1000.0 !Convert from mm to m
      if(.not.variabled50) d50 = diam(1)
      d90 = 2.5*d50 !**********************
    endif

    !Sediment size class properties
    allocate(wsfall(nsed),dstar(nsed))
    allocate(thetacr(nsed),taucr(nsed),coreyshape(nsed))  
    wsfall = 0.0; dstar = 0.0; thetacr = 0.0; taucr = 0.0
    coreyshape(:) = sedshape !Initialize with default
    do ks=1,nsed    
      !Dimensionless diameter
      dstar(ks) = diam(ks)*d2dstar
      !Corey shape factor
      coreyshape(ks) = sedclass(ks)%shape
      !Fall Velocity
      selectcase(sedclass(ks)%iws)
      case(1) !User-specified    
        wsfall(ks) = sedclass(ks)%wsfall  !Fall velocity [m/s]      
      case(3) !Wu-Wang (2006)
        wsfall(ks) = fallvel_wu_wang(viscos,coreyshape(ks),diam(ks),dstar(ks))
      case default !(0,2) !None or Soulsby (1997)
        wsfall(ks) = fallvel_soulsby(viscos,diam(ks),dstar(ks),0.0_ikind) !Concentration not included  
      endselect
      if(wsfall(ks)<0.0)then
        call diag_print_error('Problem specifying sediment fall velocity')
      endif
      !Incipient Motion, Critical shields parameter and shear stress
      selectcase(sedclass(ks)%icr)
      case(1) !User-specified shields parameter
        thetacr(ks) = sedclass(ks)%thetacr 
        taucr(ks) = thetacr(ks)*(rhosed-rhow)*grav*diam(ks) !Critical shear stress
      case(2) !User-specified critical shear stress
        taucr(ks) = sedclass(ks)%taucr   !Critical shear stress [N/m^2]    
        thetacr(ks) = taucr(ks)/((rhosed-rhow)*grav*diam(ks)) !Critical shields parameter
      case(4) !Wu and Wang (1999)
        thetacr(ks) = shields_wuwang(dstar(ks))
        taucr(ks) = thetacr(ks)*(rhosed-rhow)*grav*diam(ks) !Critical shear stress  
      case default !(0,3) !None or Soulsby (1997)
        thetacr(ks) = shields_soulsby(dstar(ks))
        taucr(ks) = thetacr(ks)*(rhosed-rhow)*grav*diam(ks) !Critical shear stress
      endselect
      if(thetacr(ks).lt.0.0 .or. taucr(ks).lt.0.0)then
        call diag_print_error('Problem specifying incipient motion')
      endif
    enddo
    
    !Concentrations and related variables
    allocate(Ctk(ncellsD,nsed),Ctk1(ncellsD,nsed),Ct(ncellsD))
    Ctk = 0.0       !Fractional total-load depth-averaged concentration 
    Ctk1 = 0.0      !Previous time step Ctk
    Ct = 0.0        !Total-load depth-averaged concentration 
    if(ntsch==2)then
      allocate(Ctk2(ncellsD,nsed))      
      Ctk2=0.0
    endif    
    allocate(Ctkstar(ncellsD,nsed),CtstarP(ncellsD,nsed),Ctstar(ncellsD))    
    Ctkstar = 0.0   !Concentration capacity
    CtstarP = 0.0   !Concentration potentional capacity
    Ctstar = 0.0    !Total-load depth-averaged concentration capacity
    allocate(rsk(ncellsD,nsed),rs(ncellsD))
    rsk = 1.0       !Fraction of suspended sediment for each size class
    rs = 1.0        !Fraction of suspended sediment for all sediment    
    allocate(cak(ncellsD,nsed),epsvk(ncellsD,nsed))
    cak = 0.0       !Concentration at reference height for each size class
    epsvk = 0.0     !Vertical mixing coefficient    
    allocate(dCtkx(ncellsD,nsed),dCtky(ncellsD,nsed))
    dCtkx=0.0; dCtky=0.0 !Gradients
    
    !Normalized residuals
    allocate(rsCtk(ncellsD),rsCtkmax(ncellsD))
    rsCtk = 0.0
    rsCtkmax = 0.0
    
    !Total-load transport
    allocate(qtx(ncellsD),qty(ncellsD))  
    qtx = 0.0      !Total-load sediment transport rate in x-direction
    qty = 0.0      !Total-load sediment transport rate in y-direction    
    
    !Morphology    
    allocate(dzb(ncellsD),zb1(ncellsD))
    if(.not.allocated(zb0))then
      allocate(zb0(ncellsD))
    endif
    zb1 = zb       !Previous time step bathymetry
    zb0 = zb       !Initial bathymetry    
    dzb = 0.0      !Bedchange term    
    
    !Adaptation variables            
    allocate(alphat(ncellsD))
    alphat = 1.0    
    if(iadapttot>=5)then
      allocate(vLbed(ncellsD),vLsus(ncellsD))   
      vLbed = Lbed
      vLsus = Lsus         
    endif
    if(iadapttot>=5 .or. iadapttot==3)then     
      allocate(vLtot(ncellsD))
      vLtot = Ltot
    endif        
    
    if(iadapttot==3)then
      call fileext(trim(aLtotfile),aext)      
      select case (aext)
      case('h5')
#ifdef XMDF_IO
      call readscalh5(aLtotfile,aLtotpath,vLtot,ierr)
#endif     
      case('txt')
        call readscalTxt(aLtotfile,vLtot,ierr)
      end select
      
    endif  
    
    !Avalanching                         
    a_repose = tan(a_repose*deg2rad) !Convert repose angle to tan(phi)  
            
    !Bedslope term
    allocate(Sb(ncellsD,nsed))
    Sb = 0.0       !Bed-slope term
    if(ibedslope<0) ibedslope = 0 !if not specified use none
    !if(ibedslope<0)then !Not specified
    !  selectcase(icapac)
    !  case(1);     ibedslope = 2 !for Lund-CIRP use Dey (2001) which modifies the transport rate directly
    !  case(2,3,4); ibedslope = 1 !for Van Rijn, Watanabe, and Soulbsy-Van Rijn use Bailard (1981) which modifies the critical shear
    !  case(5);     ibedslope = 3 !for Wu et al. (2000) use Wu (2004) which add a shar the grain shear stress
    !  endselect  
    !endif
    
    !Bedload transport
    allocate(qbk(ncellsD,nsed))  
    qbk = 0.0      !current timestep    
    if(isedmodel==2)then
      allocate(qbk1(ncells,nsed))  
      qbk1 = 0.0     !previous timestep
    endif
    
    !Hiding and exposure correction factor
    if(ihidexpform<0)then !Not specified so use default
      if(nsed==1 .and. variableD50)then
        ihidexpform=2 !Parker
      else    
        ihidexpform=3 !Wu
      endif  
    endif
    if(ihidexpform>0)then
      allocate(varsigma(ncellsD,nsed))
      varsigma = 1.0  !Hiding and exposure correction factor
    endif
    if(mhe<0.0)then !Not specified so use transport formula specific default value
      selectcase(icapac)
      case(1); mhe = 1.0 !Lund-CIRP based on report by Wu 2011
      case(2); mhe = 0.3 !Van Rijn
      case(3); mhe = 0.5 !Watanabe
      case(4); mhe = 0.3 !Soulbsy-Van Rijn      
      case(5); mhe = 0.6 !Wu
      case(6); mhe = 1.0 !Temporary matching (1) until guidance from BDJ
      endselect
    endif
    
    !Hard-bottom cells
    allocate(hardzb(ncellsD)) !always use
    hardzb=-999.0 !Initialize
    if(hardbottom) call hardbottom_read
    call struct_hardzb

!Check total depth of layers
!    if(hardbottom .and. .not.singlesize)then
!      do ih=1,nhard
!        i=idhard(ih)
!        zbotlay=zb(i) !Bottom elevation of layer
!        do js=1,nlay
!          zbotlay=zbotlay-db(i,js)
!          if(zbotlay<=hardzb(i))then
!            db(i,js)=db(i,js)-hardzb(i)+zbotlay            
!            db(i,js)=max(small,db(i,js))
!            db1(i,js)=db(i,js)
!          endif
!        enddo !js
!      enddo !ih
!    endif !hardbottom

    !If user decides to write out morphology acceleration factor value (with ramp only), write out the header
1000 format ('Accel. Factor',4x,'Time(hrs)')    
1001 format (2x,f8.4,6x,f10.5)
     
    if(write_smorph_ramp .and. scalemorph_rampdur .ge. 0.0) then 
      open (9,file=smorph_file,status='new')
      write(9,1000) 
      write(9,1001) 1.0, 0.0
      close (9)
    endif
    
    !Maximum # of iterations
    if(maxitersed<0)then !If not specified, estimate based on size of grid, grain sizes, and morphologic scaling factor
      maxitersed = 18 + int(0.15*(real(ncells,kind=ikind)/1000.0)**0.5*nsed**1.2*scalemorph**0.8)  !changed scalemorph to computed 'val'
    endif
    maxitersed = max(maxitersed,5) !Minimum value
!!    if(singlesize)then         
!!      maxitersed = 1 !No need for iterations with single size
!!      nswp(4) = 20 
!!      tolCtk = 1.0e-7        !Tolerance for total-load concentration
!!    endif
    maxitersed0 = maxitersed  
    
    !Global Sediment Balance
    if(sedbalance)then
      allocate(sedvolcur(0:nsed),sedvolcum(0:nsed),sedvolnet(0:nsed)) !(cum-total,cum-sum,current;0-total,1:nsed-size classes)
      sedvolcum%storage    = 0.0
      sedvolcum%bedchange  = 0.0
      sedvolcum%erosion    = 0.0
      sedvolcum%deposition = 0.0
      sedvolcum%boundary   = 0.0
    endif

    !Wave induced sediment transport     
    if(noptset<3) wavesedtrans = .false.     
    if(wavesedtrans)then
      fsr = scaleSus*scaleOnshore*scaleRoller    !Suspended-load due to Roller
      fba = scaleBed*scaleOnshore*scaleWaveAsym  !Bed-load due to Wave assymetry
      fsm = scaleSus*scaleOffshore*scaleUndertow !Suspended-load due to Undertow
      fbm = scaleBed*scaleOffshore*scaleUndertow !Bed-load due to Undertow
      allocate(Qws(ncellsD,nsed),QwsP(ncellsD,nsed))
      Qws = 0.0; QwsP = 0.0
    endif        
    
    !For watanabe formula
    Awidg = Awatan/rhow*gravinv        
    
    !Boundary conditions
    call sedbnd_init

    !Total load correction factor
    allocate(btk(ncellsD,nsed),btk1(ncellsD,nsed))
    btk = betatot    !Total load correction factor
    btk1 = betatot   !Total load correction factor
    if(ntsch==2)then
      allocate(btk2(ncellsD,nsed))
      btk2 = betatot 
    endif
    if(ibt>0) call bs_init
    
    !Suspended load correction factor
    allocate(bsk(ncellsD,nsed))
    bsk=betatot

    !Bed load velocity
    allocate(ubk(ncellsD,nsed))
    ubk=0.0

    !Constant mixing coefficient
    if(sedconstmix)then
      schmidt = 1.0e25
    else
      cmixsed = 0.0
    endif  
    
    !Maximum sediment concentration
    Cteqmax = min(Cteqmax,(1.0-poros)*rhosed)
    
    return
    endsubroutine sed_init    
        
!**************************************************
    subroutine sed_print()
! Prints the sediment transport setup to the screen
! and the diagnositic file
! written by Alex Sanchez, USACE-CHL
!**************************************************   
#include "CMS_cpp.h"    
    use size_def    
    use sed_def
    use comvarbl,  only: dtime,nsolv,ndsch,ntsch,wtsch,advsc
    use diag_def
    use solv_def,  only: asolv
    use der_def
    use const_def, only: deg2rad
    use prec_def
    use tool_def,  only: vstrlz
    implicit none
    integer :: i,ii,j,jj,ks,iunit(2),ierr
    real(ikind) :: pbklow(nsed,nlay),pbkhigh(nsed,nlay),pbkmean(nsed,nlay)

    if(.not.singlesize)then
      pbklow = 1.0e6; pbkhigh = -1.0e6; pbkmean = 0.0
      do i=1,ncells
        do ks=1,nsed
          do j=1,nlay
            pbklow(ks,j) = min(pbklow(ks,j),pbk(i,ks,j))
            pbkhigh(ks,j) = max(pbkhigh(ks,j),pbk(i,ks,j))
            pbkmean(ks,j) = pbkmean(ks,j) + pbk(i,ks,j)
          enddo !j layer
        enddo !ks size
      enddo !i cell
      pbklow = pbklow*100.0
      pbkhigh = pbkhigh*100.0
      pbkmean = pbkmean*100.0/real(ncells,kind=ikind)
    endif

111 format(' ',A,T40,A)
233 format(' ',A,T40,I0)
354 format(' ',A,T40,A,A)    !Added for vstrlz function results
445 format(' ',A,T40,F0.2,A)
466 format(' ',A,T40,F0.3,A)
453 format(' ',A,T40,I2,A,I2)
955 format(' ',A,T40,F0.4,A)    
723 format(' ',3x,I3,1x,4(3x,F10.4),2x,F10.3)  
774 format(' ',5x,3(F8.3,F8.3,1x))
246 format(' ',3x,I5,3(2x,F8.3),4x,3(3x,F8.2))
945 format('    D',I02,' Dataset',A,A)
    
    iunit = (/6, dgunit/)
    open(dgunit,file=dgfile,access='append')
    
    do i=1,2
      write(iunit(i),*)     
      if(.not.sedtrans)then
        write(iunit(i),111)     'Sediment Transport:','OFF'
        if(i==1) cycle
        exit
      endif  
    
      write(iunit(i),111)       'Sediment Transport:','ON'
      write(iunit(i),111)       '  Transport Model:',trim(asedmodel(isedmodel))
      if(isedmodel==1)then
        write(iunit(i),354)     '  Max Total-load Conc. Capacity:',trim(vstrlz(Cteqmax,'(f0.2)')),' kg/m^3'
      endif
      write(iunit(i),111)       '  Transport Capacity Formula:',trim(acapac(icapac))
      write(iunit(i),445)       '  Timing'
      write(iunit(i),354)       '    Transport Time Step:',trim(vstrlz(dtime,'(f0.2)')),' sec'
      write(iunit(i),354)       '    Morphologic Time Step:',trim(vstrlz(dtime,'(f0.2)')),' sec'   
      write(iunit(i),354)       '    Morphology Starting Time:',trim(vstrlz(tStartMorph,'(f0.2)')),' hrs'
      if(nsed>1)then
        write(iunit(i),354)     '    Composition Starting Time:',trim(vstrlz(tStartBedComp,'(f0.2)')),' hrs'
      endif  
      
      write(iunit(i),445)       '  Sediment Properties'
      write(iunit(i),354)       '    Density:',trim(vstrlz(rhosed,'(f0.2)')),' kg/m^3'      
      write(iunit(i),354)       '    Porosity:',trim(vstrlz(poros,'(f0.2)'))
    
      if(singlesize)then !Single-size
        write(iunit(i),354)     '    Characteristic Diameter:',trim(vstrlz(diam(1)*1000.0,'(f0.3)')),'mm'
        write(iunit(i),354)     '    Fall Velocity:',trim(vstrlz(wsfall(1),'(f0.4)')),' m/s'
        write(iunit(i),354)     '    Corey Shape Factor:',trim(vstrlz(coreyshape(1),'(f0.4)'))
        if(icapac==1 .or. icapac==3)then
          write(iunit(i),354)   '    Shields Parameter:',trim(vstrlz(thetacr(1),'(f0.4)'))
          write(iunit(i),354)   '    Critical Shear Stress:',trim(vstrlz(taucr(1),'(f0.4)')),' Pa'
        endif        
        if(calcd50)then
          write(iunit(i),111)   '    Grain Size Determined from D50 Dataset'
        elseif(transd50)then
          write(iunit(i),111)   '    Transport Grain Size Specified'
        elseif(constd50)then
          write(iunit(i),111)   '    Constant Grain Size Specified'       
          if(variableD50)then
            write(iunit(i),111) '      WARNING: D50 Dataset will be Ignored'
          endif   
        endif
      else   !Multiple-sizes
        !Sediment size classess
        if(icapac==1 .or. icapac==3)then
          write(iunit(i),111)   '    Size   Characteristic    Fall        Shields    Critical   Corey Shape'
          write(iunit(i),111)   '    Class   Diameter,mm   Velocity,m/s  Parameter   Stress,Pa     Factor'
          do ks=1,nsed
            write(iunit(i),723)  ks,diam(ks)*1000.0,wsfall(ks),thetacr(ks),taucr(ks),coreyshape(ks)
          enddo
        else
          write(iunit(i),111)   '    Size   Characteristic    Fall      Corey Shape'
          write(iunit(i),111)   '    Class   Diameter,mm   Velocity,m/s   Factor'
          do ks=1,nsed
            write(iunit(i),723)  ks,diam(ks)*1000.0,wsfall(ks),coreyshape(ks)
          enddo  
        endif
      endif       

      !Total load adaptation length
      write(iunit(i),111)       '  Total-load Adaptation Coefficient'
      write(iunit(i),111)       '    Method:',trim(atotm(iadapttot))
      selectcase(iadapttot)
      case(1)
        write(iunit(i),354)     '    Adaptation length:',trim(vstrlz(Ltot,'(f0.2)')),' m'
      case(2)
        write(iunit(i),354)     '    Adaptation time:',trim(vstrlz(Ttot,'(f0.2)')),' sec'
      case(3)
        write(iunit(i),111)     '    Adaptation length file:',trim(aLtotfile)
        write(iunit(i),111)     '    Adaptation length path:',trim(aLtotpath)
      case(4)
        write(iunit(i),354)     '    Erosion length:',trim(vstrlz(Lerotot,'(f0.2)')),' m'
        write(iunit(i),354)     '    Deposition length:',trim(vstrlz(Ldeptot,'(f0.2)')),' m'      
      case(5,6)
        write(iunit(i),111)     '  Bed-load Adaptation Coefficient'
        write(iunit(i),111)     '    Method:',trim(abedm(iadaptbed))
        selectcase(iadaptbed)
        case(1)
          write(iunit(i),354)   '    Adaptation length:',trim(vstrlz(Lbed,'(f0.2)')),' m'
        case(2)
          write(iunit(i),354)   '    Adaptation time:',trim(vstrlz(Tbed,'(f0.2)')),' sec'  
        case(5)
          write(iunit(i),354)   '    Depth dependant factor:',trim(vstrlz(fbed,'(f0.2)'))
        endselect !iadaptbed 
        write(iunit(i),111)     '  Suspended-load Adaptation Coefficient '
        write(iunit(i),111)     '    Method:',trim(asusm(iadaptsus))
        selectcase(iadaptsus)
        case(1)
          write(iunit(i),354)   '    Adaptation length:',trim(vstrlz(Lsus,'(f0.2)')), ' m'
        case(2)
          write(iunit(i),354)   '    Adaptation time:',trim(vstrlz(Tsus,'(f0.2)')),' sec'         
        case(3)
          write(iunit(i),354)   '    Adaptation coefficient:',trim(vstrlz(alphasus,'(f0.2)'))
        endselect !iadaptsus
      endselect !iadapttot      
    
      !Hardbottom
      if(hardbottom)then
        write(iunit(i),111)     '  Hardbottom:','ON'
        write(iunit(i),111)     '    Hardbottom file:',trim(hbfile)
        if(len_trim(hbpath)>0)then
          write(iunit(i),111)   '    Hardbottom path:',trim(hbpath)
        endif
      else
        write(iunit(i),111)     '  Hardbottom:','OFF' 
      endif  
    
      !Avalanching
      if(do_aval)then
        write(iunit(i),111)     '  Avalanching:','ON'
        write(iunit(i),354)     '    Repose angle:',trim(vstrlz(atan(a_repose)/deg2rad,'(f0.2)')),' deg'
        write(iunit(i),354)     '    Relaxation coefficient:',trim(vstrlz(relax_aval,'(f0.2)'))
        write(iunit(i),233)     '    Maximum iterations:',nmaxaval  
      else
        write(iunit(i),111)     '  Avalanching:','OFF'  
      endif          
    
      write(iunit(i),111)       '  Scaling factors '    
      write(iunit(i),354)       '    Suspended Load:',trim(vstrlz(scalesus,'(f0.2)'))
      write(iunit(i),354)       '    Bed Load:',trim(vstrlz(scalebed,'(f0.2)'))
      write(iunit(i),354)       '    Morphologic:',trim(vstrlz(scalemorph_orig,'(f0.2)'))
      if (scalemorph_rampdur .gt. 0.0) then
        write(iunit(i),354)     '      Ramp duration (hrs):',trim(vstrlz(scalemorph_rampdur,'(f0.2)'))
      endif
     
      write(iunit(i),111)       '  Bed slope effects'   
      write(iunit(i),354)       '    Diffusion coefficient:',trim(vstrlz(dcoeff,'(f0.2)'))
      write(iunit(i),111)       '    Transport Correction:',trim(abedslope(ibedslope))
      if(ibedslope==2)then !Bailard (1981)
        write(iunit(i),354)     '    Bed-load coefficient:',trim(vstrlz(betaslope,'(f0.2)'))
        write(iunit(i),354)     '    Susp-load coefficient:',trim(vstrlz(effslope,'(f0.2)'))
      endif
    
      if(variableD50 .or. .not.singlesize)then
        write(iunit(i),111)     '  Hiding and Exposure'   
        write(iunit(i),111)     '    Formulation:',trim(ahidexpform(ihidexpform))
        if(ihidexpform==2 .or. ihidexpform==3)Then
          write(iunit(i),354)   '    Coefficient:',trim(vstrlz(mhe,'(f0.2)'))
        endif
      endif
      
      if(isedmodel<=2)then       
        if(sedconstmix)then
          write(iunit(i),354)   '  Constant Mixing Coefficient:',trim(vstrlz(cmixsed,'(f0.2)'))
        else
          write(iunit(i),354)   '  Schmidt Number:',trim(vstrlz(schmidt,'(f0.2)'))
        endif 
      endif    
     
      write(iunit(i),354)       '  Inflow Loading Factor:',trim(vstrlz(facQtotin,'(f0.2)'))  
      if(isedmodel<=2)then
        selectcase(ibt)
        case(0)
          write(iunit(i),354)   '  Constant Total Load Correction Factor: ',trim(vstrlz(betatot,'(f0.2)'))
        case(1)
          write(iunit(i),111)   '  Concentration profile:','EXPONENTIAL    (Used for Total Load Correction Factor)'
        case default
          write(iunit(i),111)   '  Concentration profile:','ROUSE    (Used for Total Load Correction Factor)'
        endselect 
      endif    
    
      !Bed Composition and layering
      if(singlesize)then
        if(variableD50)then
          write(iunit(i),111)   '  Variable D50:','ON'  
          write(iunit(i),111)   '    D50 Dataset File:',trim(bedlay(1)%perdiam(ipd(50))%file)
          write(iunit(i),111)   '    D50 Dataset Path:',trim(bedlay(1)%perdiam(ipd(50))%path)
        else
          write(iunit(i),111)   '  Variable D50:','OFF'
        endif  
      else
        write(iunit(i),111)     '  Mixing Layer'    
        if(mixlayconst)then          
          write(iunit(i),111)   '    Formulation:','CONSTANT'
          write(iunit(i),354)   '    Thickness:',trim(vstrlz(dmconst,'(f0.2)'))  !Constant mixing layer thickness [m]
        else
          write(iunit(i),111)   '    Formulation:','AUTOMATIC'  
        endif
        write(iunit(i),354)     '  Minimum Layer Thickness:',trim(vstrlz(dbmin,'(f0.2)'))   !Minimum layer thickness [m]
        write(iunit(i),354)     '  Maximum Layer Thickness:',trim(vstrlz(dbmax,'(f0.2)'))   !Maximum layer thickness for new layers. Existing layers are not affected [m]    
        write(iunit(i),233)     '  Number of Bed Layers:',nlay
        do j=1,nlay
          if(bedlay(j)%ipbkinp==0) cycle
          do jj=j+1,nlay
            if(bedlay(jj)%ipbkinp/=0) exit  
          enddo
          if(j==jj-1)then
            write(iunit(i),233) '   Bed Layer: ',j
          else
            write(iunit(i),453) '   Bed Layers: ',j,' - ',jj-1
          endif
          write(iunit(i),111)   '    Layer Thickness Method:',trim(adbinp(bedlay(j)%idbinp))
          selectcase(bedlay(j)%idbinp)
          case(0,1,2)
            write(iunit(i),354) '    Layer Thickness:',trim(vstrlz(bedlay(j)%dbconst,'(f0.2)')),' m'
          case(3)
            write(iunit(i),111) '    Layer Thickness Dataset File:',trim(bedlay(j)%dbfile)
            write(iunit(i),111) '    ----------------------- Path:',trim(bedlay(j)%dbpath)
          endselect
          write(iunit(i),111)   '    Composition Method:',trim(apbkinp(bedlay(j)%ipbkinp))
          selectcase(bedlay(j)%ipbkinp)
          case(1) !D50_SIGMA
            write(iunit(i),111) '    D50 Dataset File:',trim(bedlay(j)%perdiam(ipd(50))%file)
            write(iunit(i),111) '    ----------- Path:',trim(bedlay(j)%perdiam(ipd(50))%path)
            if (bedlay(j)%geostddev .ge. 0.0) then
              write(iunit(i),354) '    Geometric standard deviation:',trim(vstrlz(bedlay(j)%geostddev,'(f0.3)')),' mm'    
            endif
            write(iunit(i),111) '    Layer grain size distribution summary (all cells)'
            write(iunit(i),111) '            Lower     Upper  Characteristic   Minimum    Maximum     Mean'  
            write(iunit(i),111) '    Class  bound,mm  bound,mm  diameter,mm   fraction,% fraction,% fraction,%' 
            do ks=1,nsed
              write(iunit(i),246) ks, diamlim(ks)*1000.0,diamlim(ks+1)*1000.0,&
                 diam(ks)*1000.0,pbklow(ks,j),pbkhigh(ks,j),pbkmean(ks,j)
            enddo !ks
          case(2) !D16_D50_D84
            write(iunit(i),111) '    D16 Dataset File:',trim(bedlay(j)%perdiam(ipd(16))%file)
            write(iunit(i),111) '    ----------- Path:',trim(bedlay(j)%perdiam(ipd(16))%path)            
            write(iunit(i),111) '    D50 Dataset File:',trim(bedlay(j)%perdiam(ipd(50))%file)
            write(iunit(i),111) '    ----------- Path:',trim(bedlay(j)%perdiam(ipd(50))%path)
            write(iunit(i),111) '    D84 Dataset File:',trim(bedlay(j)%perdiam(ipd(84))%file)
            write(iunit(i),111) '    ----------- Path:',trim(bedlay(j)%perdiam(ipd(84))%path)
            write(iunit(i),111) '    Layer grain size distribution summary (all cells)'
            write(iunit(i),111) '            Lower     Upper  Characteristic   Minimum    Maximum     Mean'  
            write(iunit(i),111) '    Class  bound,mm  bound,mm  diameter,mm   fraction,% fraction,% fraction,%' 
            do ks=1,nsed
              write(iunit(i),246) ks, diamlim(ks)*1000.0,diamlim(ks+1)*1000.0,&
                 diam(ks)*1000.0,pbklow(ks,j),pbkhigh(ks,j),pbkmean(ks,j)
            enddo !ks
          case(3) !D35_D50_D90
            write(iunit(i),111) '    D35 Dataset File:',trim(bedlay(j)%perdiam(ipd(35))%file)
            write(iunit(i),111) '    ----------- Path:',trim(bedlay(j)%perdiam(ipd(35))%path)
            write(iunit(i),111) '    D50 Dataset File:',trim(bedlay(j)%perdiam(ipd(50))%file)
            write(iunit(i),111) '    ----------- Path:',trim(bedlay(j)%perdiam(ipd(50))%path)
            write(iunit(i),111) '    D90 Dataset File:',trim(bedlay(j)%perdiam(ipd(90))%file)
            write(iunit(i),111) '    ----------- Path:',trim(bedlay(j)%perdiam(ipd(90))%path)
            write(iunit(i),111) '    Layer grain size distribution summary (all cells)'
            write(iunit(i),111) '            Lower     Upper  Characteristic   Minimum    Maximum     Mean'  
            write(iunit(i),111) '    Class  bound,mm  bound,mm  diameter,mm   fraction,% fraction,% fraction,%' 
            do ks=1,nsed
              write(iunit(i),246) ks, diamlim(ks)*1000.0,diamlim(ks+1)*1000.0,&
                 diam(ks)*1000.0,pbklow(ks,j),pbkhigh(ks,j),pbkmean(ks,j)
            enddo !ks
          case(4) !SIZE_CLASS_FRACTIONS
            write(iunit(i),111) '            Lower     Upper  Characteristic  Size Class'
            write(iunit(i),111) '    Class  bound,mm  bound,mm  diameter,mm   Fraction,% '   
            do ks=1,nsed
              write(iunit(i),246) ks,diamlim(ks)*1000.0,diamlim(ks+1)*1000.0,diam(ks)*1000.0,bedlay(j)%pbconst(ks)*100.0
            enddo
          case(6) !PRECENTILES
            do ii=1,nperdiam
              if(bedlay(j)%perdiam(ii)%inp)then 
                write(iunit(i),945) iper(ii),' File:',trim(bedlay(j)%perdiam(ii)%file)
                write(iunit(i),945) iper(ii),' Path:',trim(bedlay(j)%perdiam(ii)%path)
              endif
            enddo   
            write(iunit(i),111) '    Layer grain size distribution summary (all cells)'
            write(iunit(i),111) '            Lower     Upper  Characteristic   Minimum    Maximum     Mean'  
            write(iunit(i),111) '    Class  bound,mm  bound,mm  diameter,mm   fraction,% fraction,% fraction,%' 
            do ks=1,nsed
              write(iunit(i),246) ks, diamlim(ks)*1000.0,diamlim(ks+1)*1000.0,&
                 diam(ks)*1000.0,pbklow(ks,j),pbkhigh(ks,j),pbkmean(ks,j)
            enddo !ks
          endselect          
        enddo !j layer        
      endif
      
      if(sedbalance)then
        write(iunit(i),111)     '  Global Sediment Balance Log:','ON'
      else
        write(iunit(i),111)     '  Global Sediment Balance Log:','OFF'  
      endif
      
#ifdef DEV_MODE
      if(erosdry%calc)then
        write(iunit(i),111)     '  Erosion of dry cells:','ON'
        write(iunit(i),354)     '    Erosion transfer:',trim(vstrlz(erosdry%fac,'(f0.2)'))
      else
        write(iunit(i),111)     '  Erosion of dry cells:','OFF'
      endif
#endif
    
      write(iunit(i),111)       '  Numerical Methods'
      write(iunit(i),111)       '    Matrix Solver:',trim(asolv(nsolv))
      if(ntsch==1)then
        write(iunit(i),111)     '    Temporal Scheme:','TWO-LEVEL'
      else
        write(iunit(i),111)     '    Temporal Scheme:','THREE-LEVEL'
        write(iunit(i),354)     '    Implicit Weighting Factor:',trim(vstrlz(wtsch,'(f0.2)')) 
      endif
      write(iunit(i),111)       '    Advection Scheme:',trim(advsc(ndsch))
      write(iunit(i),111)       '    Spatial Derivative Scheme:',trim(ader(nder))
      write(iunit(i),233)       '    Maximum transport iterations:',maxitersed    
      if(.not.singlesize)then      
        write(iunit(i),233)     '    Maximum bed iterations:',itermaxzb   
      endif    
    enddo
    close(dgunit)
    
    return
    end subroutine sed_print
    
!***********************************************************************
    subroutine hardbottom_read
! Reads the hardbottom cell id and depth from the model parameters file
! Creates hardbottom arrays idhard and hardbed, with length nhard
!
! written by Alex Sanchez, USACE-CHL
!***********************************************************************   
#include "CMS_cpp.h"
    use size_def, only: ncells,ncellsD,ncellsfull
    use geo_def
    use comvarbl
    use sed_def
#ifdef XMDF_IO
    use xmdf
    use in_xmdf_lib, only: readscalh5
#endif
    use in_lib, only: readscalTxt
    use diag_lib
    use diag_def, only: dgunit,dgfile
    use prec_def
    
    implicit none
    integer :: i,ih,ierr,idhardtemp(ncellsD)
    integer :: hbwarn(ncellsD), nhbwarn
    character(len=100) :: msg2,msg3,msg4,msg5
    character(len=10) :: aext

    call fileext(trim(hbfile),aext)      
    select case (aext)
    case('h5')
#ifdef XMDF_IO
      if (hbpath == 'Datasets/') then                 !This should handle SMS writing out empty HARDBOTTOM_DATASET card.
        call diag_print_warning('Hard Bottom dataset not specified - turning off hardbottom','')
        hardbottom=.false.
        return
      endif  
      call readscalh5(hbfile,hbpath,hardzb,ierr)
    
      if(ierr/=0)then
        write(msg2,*) '  File: ',trim(hbfile)
        write(msg3,*) '  Path: ',trim(hbpath)
        call diag_print_error('Could not open hardbottom dataset: ',msg2,msg3)
      endif
#endif
    case('txt')
      call readscalTxt(hbfile,hardzb,ierr)
      if(ierr/=0)then
        write(msg2,*) '  File: ',trim(hbfile)
        call diag_print_error('Could not open hardbottom dataset: ',msg2)
    endif
    end select
    
!Find number and id of hardbottom cells
    nhbwarn=0
    hbwarn=0
    do i=1,ncells
      if(abs(hardzb(i)+999.0)>1.0e-4)then
        nhard=nhard+1
        idhardtemp(nhard)=i
        hardzb(i) = -hardzb(i) !Note sign change from depths to elevations
        !Check elevations
        if(hardzb(i)>zb(i)+1.0e-3)then
          !Removing this message for now.  Adding IDs of cells to an array all to be printed at one time.  MEB 12/11/2018
          !write(msg2,*) '  Cell: ',mapid(i)
          !write(msg3,*) '  Hard bottom depth: ',-hardzb(i),' m'
          !write(msg4,*) '  Water Depth: ',-zb(i),' m'
          !write(msg5,*) '  Setting bed elevation as hard bottom'
          !call diag_print_warning('Specified hard bottom above bed elevation',msg2,msg3,msg4,msg5)
          nhbwarn = nhbwarn + 1
          hbwarn(nhbwarn) = i
        endif
        hardzb(i)=min(hardzb(i),zb(i)) !Note: Hard bottom must be at or below the bed elevation
      endif
    enddo
    
 99 format(15(i0,x))
100 format('Specified hard bottom above bed elevation for ',i0,' cells.')
101 format('Cell IDs written to file: "hb_warning.txt"')

    if(nhbwarn .gt. 0) then 
      write(msg2,100) nhbwarn
      write(msg3,101) 
      call diag_print_warning(msg2,msg3,'')
!      write(*,99)      (hbwarn(i),i=1,nhbwarn)
      
      open(200,file='hb_warning.txt',status='unknown') 
      write(200,99) (mapid(hbwarn(i)),i=1,nhbwarn)  !changed to the ID of the cell as in SMS.  01/29/2019
      close(200)
    endif        
    
    !Copy hardbottom info to smaller arrays
    if(nhard>0)then
      hardbottom=.true.
      allocate(hardbed(nhard),idhard(nhard))
      do ih=1,nhard
        idhard(ih)=idhardtemp(ih) !No mapping necessary
        if(hardzb(idhard(ih))<-900.0)then
          write(msg2,*) '  Cell: ',idhard(ih)
          write(msg3,*) '  Hard bottom: ',hardzb(idhard(ih))  
          call diag_print_error('Could not calculate hard bottom ID',msg2,msg3)
        endif
        hardbed(ih)=hardzb(idhardtemp(ih))  
      enddo   
    else
      hardbottom=.false.      
    endif
    
    return
    endsubroutine hardbottom_read
    
!******************************************************************
    subroutine sed_step_stat
! Calculates the sediment transport statistics for each time step
! Author: Alex Sanchez, USACE-CHL
!******************************************************************   
    use size_def
    use geo_def, only: mapid
    use flow_def, only: iwet
    use sed_def, only: Ct,dzb,nsed,pbk,pbk1
    use diag_def, only: msg
    use diag_lib
    use prec_def
    implicit none
    integer :: i,ks,idCt,idCt2,iddzb,iddzb2,iddpbm,iddpbm2,ierr
    real(ikind) :: Ctmax,Ctmax2,dzbxtr,dzbxtr2,dpbmxtr,dpbmxtr2,dpbm  
    
    !Initialize
    Ctmax =0.0;  Ctmax2 =0.0;  idCt=1;   idCt2=1
    dzbxtr=0.0;  dzbxtr2=0.0;  iddzb=1;  iddzb2=1
    dpbmxtr=0.0; dpbmxtr2=0.0; iddpbm=1; iddpbm2=1
!$OMP PARALLEL FIRSTPRIVATE(idCt2,Ctmax2,iddzb2,dzbxtr2,iddpbm2,dpbmxtr2,dpbm)
!$OMP DO PRIVATE(i)
    do i=1,ncells
      if(iwet(i)==1)then
        if(Ct(i)>abs(Ctmax2))then
          Ctmax2=Ct(i)
          idCt2=i
        endif
        if(abs(dzb(i))>abs(dzbxtr2))then
          dzbxtr2=dzb(i)
          iddzb2=i
        endif
      endif
    enddo
!$OMP END DO

!$OMP CRITICAL
    if(Ctmax2>Ctmax)then
      Ctmax=Ctmax2
      idCt=idCt2
    endif
    if(abs(dzbxtr2)>abs(dzbxtr))then
      dzbxtr=dzbxtr2
      iddzb=iddzb2
    endif
!$OMP END CRITICAL

    if(nsed>1)then
!$OMP DO PRIVATE(i,ks)
      do i=1,ncells
        if(iwet(i)==1)then
          do ks=1,nsed
            dpbm=pbk1(i,ks)-pbk(i,ks,1)  
            if(abs(dpbm)>abs(dpbmxtr2))then
              dpbmxtr2=dpbm
              iddpbm2=i
            endif
          enddo
        endif
      enddo
!$OMP END DO
!$OMP CRITICAL
      if(abs(dpbmxtr2)>abs(dpbmxtr))then
        dpbmxtr=dpbmxtr2
        iddpbm=iddpbm2
      endif
!$OMP END CRITICAL
    endif
!$OMP END PARALLEL

    if(ncellsimple>0)then
      idCt = mapid(idCt)
      iddzb = mapid(iddzb)
      if(nsed>1) iddpbm = mapid(iddpbm)
    endif
    
626 format('   Ct(',I6,')=',F9.5,', dzb(',I6,')=',F9.5)
727 format('   Ct(',I6,')=',F9.5,', dzb(',I6,')=',F9.5,', dpbk(',I6,')=',F9.5)
    !if(abs(Ctmax)>1.0e5 .or. abs(dzbxtr)>1.0e3) return
    if(nsed>1)then
      write(msg,727,iostat=ierr) idCt,Ctmax,iddzb,dzbxtr,iddpbm,dpbmxtr
    else
      !write(msg,626,iostat=ierr) idCt,Ctmax,iddzb,dzbxtr
    endif
    !call diag_print_message(msg)
    
    return
    endsubroutine sed_step_stat

!************************************************************************
    subroutine sed_total
! Updates cmb, Ct, Ctstar, rs, qtx and qty      
! written by Alex Sanchez, USACE-CHL;  Weiming Wu, NCCHE      
!************************************************************************      
    use size_def
    use geo_def, only: dzbx,dzby
    use flow_def, only: u,v,us,vs,vis,h,iwet !bdj added us, vs
    use wave_flowgrid_def, only: Wunitx,Wunity,Wang,Whgt,Wlen,Wper !bdj added Wang,Whgt,Wlen,Wper
    use fric_def, only: fricbedslope,cmb
    use const_def, only: small
    use sed_def
    use prec_def
    implicit none
    integer :: i
    real(ikind) :: fac
    real(ikind) :: CSsg,Hrms,sigT,CSPb,qb,qbx,qby,qsx,qsy !bdj    !removed CSSlp and CSBlp because they are initialized or user-specified 6/7/2019 bdj
!$OMP PARALLEL DO PRIVATE(i,fac)              
    do i=1,ncells
      !--- Total-load sediment concentrations 
      ! and fraction of suspended sediments ---  
      Ct(i)=sum(Ctk(i,:))
      Ctstar(i)=sum(Ctkstar(i,:))
      rs(i)=sum(Ctk(i,:)*rsk(i,:))/max(Ct(i),small)       
      !--- Total-load sediment transport vectors ---
      !Advective transports
      qtx(i)=u(i)*h(i)*Ct(i)
      qty(i)=v(i)*h(i)*Ct(i)         
      !Diffusive transport (approximate)
      fac=h(i)*vis(i)/schmidt  
      qtx(i)=qtx(i)-fac*sum(rsk(i,:)*dCtkx(i,:))
      qty(i)=qty(i)-fac*sum(rsk(i,:)*dCtky(i,:))
      !Bedslope transport (approximate)
      qtx(i)=qtx(i)-dcoeff*sum(qbk(i,:))*dzbx(i)
      qty(i)=qty(i)-dcoeff*sum(qbk(i,:))*dzby(i)  
    enddo !i
!$OMP END PARALLEL DO 
      
    !--- Wave-induced sediment transports ----------
     if(wavesedtrans)then
!$OMP PARALLEL DO PRIVATE(i)         
      do i=1,ncells
        qtx(i)=qtx(i)+sum(Qws(i,:))*Wunitx(i)
        qty(i)=qty(i)+sum(Qws(i,:))*Wunity(i)
      enddo
!$OMP END PARALLEL DO    
    endif      

! bdj some heavy-handed cshore mods
if (icapac.eq.6) then
   !CSslp = 0.5            !commented 6/7/2019 bdj
   !CSblp = 0.001          !commented 6/7/2019 bdj 
   CSsg = rhosed/1000.
!$OMP PARALLEL DO PRIVATE(qsx,qsy,Hrms,sigT,i,qb,qbx,qby,CSPb)   
   do i=1,ncells
      !--- Total-load sediment concentrations 
      ! and fraction of suspended sediments ---  
      Ct(i)=sum(Ctk(i,:))
      Ctstar(i)=sum(Ctkstar(i,:))
      rs(i)=sum(Ctk(i,:)*rsk(i,:))/max(Ct(i),small)       
      !--- Total-load sediment transport vectors ---
      !Suspended advective transport
      qsx=1.*(u(i)-(CSslp*us(i)))*h(i)*Ct(i)
      qsy=1.*(v(i)-(CSslp*vs(i)))*h(i)*Ct(i)         

      !Bedload transport
      Hrms = Whgt(i)/sqrt(2.)  
      sigT = (Hrms/sqrt(8.))*(Wlen(i)/Wper(i))/h(i)
      call prob_bedload(sigT,Wper(i),CSsg,diam(1),u(i),v(i),CSPb)
      !write(*,*),'bdj coming from prob_bedload,sigT,Wper(i),nsed,diam(1),Pb',sigT,Wper(i),nsed,diam(1),Pb
      qb = rhosed*(CSPb*CSblp*sigT**3.)/(9.81*(CSsg-1.))
      qbx = 1.*qb*wunitx(i)
      qby = 1.*qb*wunity(i)
      
      !Total = sus + bedload
      qtx(i) = qsx + qbx
      qty(i) = qsy + qby
      
      !write(1234,*)i,qbx,qsx,qtx(i) !bdj
      
      ! if(i.ge.1514.and.i.le.1524) then 
      !    !write(*,*),'bdj i wunitx(i) wunity(i) wang(i) cos sin',i,wunitx(i),wunity(i),wang(i),cos(wang(i)),sin(wang(i))
      !    write(*,*),'bdj i Hrms,sigT,qsx,qbx,qtx(i)',Hrms,sigT,qsx,qbx,qtx(i)
      ! endif
   enddo
!$OMP END PARALLEL DO   
endif

! bdj end some heavy-handed cshore mods
!$OMP PARALLEL DO PRIVATE(i)          
    do i=1,ncells
      qtx(i)=iwet(i)*qtx(i)  !kg/m/s
      qty(i)=iwet(i)*qty(i)  !kg/m/s
    enddo
!$OMP END PARALLEL DO
      
    return
    endsubroutine sed_total
    
    subroutine prob_bedload(sigT,Tp,CSsg,d50,u,v,CSPb)
! calculates the probability of bedload transport, Pb
! written by Brad Johnson, USACE-CHL;
! shear is on the basis of waves alone at this point	   
! updated from Brad's branch - 05/15/2020
!************************************************************************      
    use prec_def
    implicit none
    integer :: i,numsteps
    real(ikind) :: sigT,Tp,CSsg,d50,u,v,CSPb,dr
    real(ikind) :: fw,shields,tau_c,tau,dum
    real,DIMENSION(101) :: r,f,fdum,ta
 
    fw = 0.02
    shields = 0.05
    tau_c = 9810.*(CSsg-1.)*d50*shields
    r = (/ (I, I = -50,50,1) /)
    r = r/10.
    dr = r(2)-r(1)
    f = 1/SQRT(2.*3.14)*exp(-r**2./2)
    fdum = f
    ta = 1000*fw/2.*(r*sigT)**2.
    where(ta<tau_c)
       fdum = 0.
    end where
    CSPb = sum(fdum)*dr
    ! write(*,*) r
    ! write(*,*) ta
    ! write(*,*) fdum
    ! write(*,*) 'CSPb ',CSPb


    ! numsteps = 30 
    ! dum = 0.
    ! do i = 1,numsteps
    !   t = (float(i)-1.)/float(numsteps)*Tp
    !   u = sqrt(2.)*sigT*sin(t*2.*3.14/Tp)
    !   tau = 1000*fw/2.*u**2.
    !   if(tau.gt.tau_c) dum = dum+1.
    !   !write(*,*),'bdj sigT,t,u,tau,tau_c,dum',sigT,t,u,tau,tau_c,dum
    ! enddo
    ! CSPb = dum/float(numsteps)
    ! write(*,*)'bdj in prob_bedload',sigT,Tp,CSPb    

    return
    endsubroutine prob_bedload

    
!***********************************************************************   
    subroutine print_sedvar(i,ks)
! Prints sediment transport variables for debugging
!
! written by Alex Sanchez, USACE-CHL
!***********************************************************************
    use flow_def
    use geo_def, only: icol,irow,mapid,zb,igridtype
    use sed_def
    use wave_flowgrid_def
    use cms_def
    use comvarbl
    use diag_def
    implicit none
    integer :: i,ii,ks,iunit(2)
      
    iunit(1) = 6
    iunit(2) = dgunit
    open(dgunit,file=dgfile,access='append') 
    do ii=1,2
      write(iunit(ii),*)
      if(allocated(mapid))then
        write(iunit(ii),*) 'id =',mapid(i),'ks =',ks   
      else
        write(iunit(ii),*) 'id =',i,'ks =',ks     
      endif
      if(igridtype==0)then
        write(iunit(ii),*) ' Column =',icol(i),' Row =',irow(i)
      endif
      write(iunit(ii),*) 'Ctk1(i,ks) =',Ctk1(i,ks)
      write(iunit(ii),*) 'Ctk(i,ks) =',Ctk(i,ks)
      write(iunit(ii),*) 'Ctkstar(i,ks) =',Ctkstar(i,ks)
      write(iunit(ii),*) 'CtstarP(i,ks) =',CtstarP(i,ks)
      write(iunit(ii),*) 'pbk(i,ks,1) =',pbk(i,ks,1)
      write(iunit(ii),*) 'Sb(i,ks) =',Sb(i,ks)
      write(iunit(ii),*) 'alphat(i) =' ,alphat(i)
      write(iunit(ii),*) 'su(i) =',su(i)
      write(iunit(ii),*) 'sp(i) =',sp(i)
      write(iunit(ii),*) 'h(i) =',h(i)
      write(iunit(ii),*) 'eta(i) =',eta(i)
      write(iunit(ii),*) 'zb(i) =',zb(i)
      write(iunit(ii),*) 'dzb(i) =',dzb(i)
      write(iunit(ii),*) 'uv(i) =',uv(i)
      write(iunit(ii),*) 'u(i) =',u(i)
      write(iunit(ii),*) 'v(i) =',v(i)
      if(noptset>=3)then
        write(iunit(ii),*) 'Whgt(i) =',Whgt(i)
        write(iunit(ii),*) 'Worb(i) ',Worb(i)
        write(iunit(ii),*) 'Worbrep(i) =',Worbrep(i)
        write(iunit(ii),*) 'Wang(i) =',Wang(i)
        write(iunit(ii),*) 'Wper(i) =',Wper(i)    
        write(iunit(ii),*) 'Wlen(i) =',Wlen(i)
      endif
    enddo
    close(dgunit)
    
    return
    endsubroutine print_sedvar
    
!***********************************************************************   
    subroutine update_const_wave
! Prints sediment transport variables for debugging
! written by Alex Sanchez, USACE-CHL
!***********************************************************************
    use size_def
    use flow_def
    use sed_def  
    use wave_flowgrid_def
    use cms_def
    use const_def, only: pi,small
    use prec_def
    
    implicit none
    integer :: i
    real(ikind) :: wavelen, arg 
     
    do i=1,ncells
      Worb(i) = 0.0
      if(waveheight>0.0) then
        !Wavelength, Fenton (1990) Coatal Engineering 1.7% error, Alex
        wavelen = grav*waveperiod**2/2.0/pi+small !deepwater wavelength
        arg = tanh((2.0*pi*h(i)/wavelen)**0.75)**0.66667
        wavelen = wavelen*arg                    
        !Bottom wave orbital velocity
        arg = 2.0*pi*h(i)/wavelen !kh
        arg = min(arg,20.0)
        if(waveheight<=0.78*h(i))then
          Worb(i) = pi*waveheight/(waveperiod*sinh(arg))
        else
          Worb(i) = 0.78*pi*h(i)/(waveperiod*sinh(arg))
        endif
      endif    
    enddo
        
    return
    endsubroutine
    
!***********************************************************************
    subroutine sed_wave
! Calculates the source term in the bed change equation due to waves
! asymmetry, roller and undertow
! written by Alex Sanchez, USACE-ERDC-CHL
!***********************************************************************    
    use size_def
    use geo_def
    use flow_def
    use bnd_def
    use sed_def, only: nsed,nhard,idhard,hardzb,&
        Qws,QwsP,Sb,pbk,Ctk,Ctkstar,rhosed
    use wave_flowgrid_def, only: Wunitx,Wunity
    use interp_def
     use cms_def, only: cmswave 
     use const_def, only: small
     use prec_def
    implicit none
    integer :: i,j,k,ks,ih,ibnd,nck
    real(ikind) :: Sw,fk,fp
    real(ikind) :: awkx,awky,asum
        
    !--- Actual transport based on composition and loading ------
    do i=1,ncells
      do ks=1,nsed
!!        Qws(i,ks)=pbk(i,ks,1)*QwsP(i,ks)*max(min(Ctk(i,ks)/(Ctkstar(i,ks)+small),1.5),0.0) !m^2/sec, Bug fix, changed k to ks
        Qws(i,ks)=pbk(i,ks,1)*QwsP(i,ks)  !m^2/sec, Bug fix, changed k to ks
      enddo  
    enddo        
    
    !--- Calculate coefficients ------
    !Note that these coefficients are independant of the sediment size class ks
    !So they only need to be calculated once
    do i=1,ncells
      do k=1,ncface(i)     
        nck=cell2cell(k,i)        
        if(nck<=ncells)then 
          fk=fintp(k,i); fp=1.0-fk          
          awkx=fnx(k,i)*(Wunitx(nck)*fk+Wunitx(i)*fp) !cos(Wang(nck))
          awky=fny(k,i)*(Wunity(nck)*fk+Wunity(i)*fp) !sin(Wang(nck))
        else
          awkx=fnx(k,i)*Wunitx(i) !cos(Wang(nck))
          awky=fny(k,i)*Wunity(i) !sin(Wang(nck))
        endif
        acoef(k,i)=iwet(i)*iwet(nck)*ds(k,i)*(awkx+awky)
      enddo
    enddo
    
    !--- Boundary Conditions -----------
    !Apply open boundary conditions
    do ibnd=1,nbndstr
      do j=1,bnd_str(ibnd)%ncells
        i=bnd_str(ibnd)%cells(j)
        k=bnd_str(ibnd)%faces(j)
!        acoef(k,i)=0.0        
        Qws(cell2cell(k,i),:)=Qws(i,:)  !Zero-gradient BC
      enddo !im  
    enddo !ibnd       

    !Apply hardbottom
    do ih=1,nhard
      i=idhard(ih)
      if(zb(i)<=hardzb(i))then
        do k=1,ncface(i)
          nck=cell2cell(k,i)
          if(zb(nck)<=hardzb(nck) .or. zb(i)>zb(nck))then
            acoef(k,i)=0.0
            acoef(llec2llec(k,i),nck)=0.0   ! This is important.  Wu
          endif
        enddo !k    
      endif
    enddo !ih  

    !--- Solve for Sw --------------
    do ks=1,nsed 
    !Calculate bed-slope term    
      do i=1,ncells
        Sw=0.0; asum=0.0
        do k=1,ncface(i)
          Sw=Sw+min(acoef(k,i),0.0)*Qws(cell2cell(k,i),ks)
          asum=asum+max(acoef(k,i),0.0)
        enddo !k
        Sw=Sw+asum*Qws(i,ks)
        Sb(i,ks)=Sb(i,ks)+Sw/areap(i)
        if(Sb(i,ks)>50.0)then
          write(*,*) 'acoef(k,i) ',(acoef(k,i),k=1,ncface(i))
          write(*,*) 'pbk(i,ks,1)',pbk(i,ks,1)
          write(*,*) 'Qws(i,ks)',Qws(i,ks)
          write(*,*) 'Ctk(i,ks)',Ctk(i,ks)
          write(*,*) 'Ctkstar(i,ks)',Ctkstar(i,ks)
        endif
      enddo !i 
    enddo !ks 

    return
    endsubroutine

!*******************************************************************************
    subroutine sed_balance()
! Estimate total sediment budget
! written by Alex Sanchez, USACE-ERDC-CHL: Weiming Wu, NCCHE 
!*******************************************************************************         
    use size_def
    use geo_def
    use flow_def
    use comvarbl
    use diag_def
    use bnd_def
    use sed_def 
    use prec_def     
    implicit none
    integer :: i,j,k,ks,ibnd,nck,iunit(2)
    real(ikind) :: vol,fac,vel !,volabschg
    !real(ikind),dimension(nsed) :: volcumchgk,volchgk,volcumstgk,volstgk,volbalk,volbndk
    !real(ikind),dimension(nsed) :: volerosk,voldepok,volcumerosk,volcumdepok,volerrk
    
    !Current volumes
    sedvolcur(:)%storage    = 0.0
    sedvolcur(:)%bedchange  = 0.0
    sedvolcur(:)%erosion    = 0.0
    sedvolcur(:)%deposition = 0.0 
    sedvolcur(:)%boundary   = 0.0
    
    !Net Volumes (computed using total bed change and concentrations)
    sedvolnet(:)%storage    = 0.0
    sedvolnet(:)%bedchange  = 0.0
    sedvolnet(:)%erosion    = 0.0
    sedvolnet(:)%deposition = 0.0 
    sedvolnet(:)%boundary   = 0.0
    
  !=== Morphology, bed change and storage volumes ====      
    do i=1,ncells
      !Total Bedchange, deposition and erosion volumes
      vol = solid*(zb(i)-zb0(i))*areap(i)
      sedvolnet(0)%bedchange = sedvolnet(0)%bedchange + vol
      if(vol>0.0)then
        sedvolnet(0)%deposition = sedvolnet(0)%deposition + vol
      else
        sedvolnet(0)%erosion = sedvolnet(0)%erosion + vol
      endif
      !Fraction Bedchange
      do ks=1,nsed 
        if(nsed>1)then      
          vol = solid*dzbk(i,ks)*areap(i)
        else
          vol = solid*dzb(i)*areap(i)
        endif
        sedvolcur(0)%bedchange = sedvolcur(0)%bedchange + vol
        sedvolcur(ks)%bedchange = sedvolcur(ks)%bedchange + vol
        if(vol>0.0)then
          sedvolcur(0)%deposition = sedvolcur(0)%deposition + vol
          sedvolcur(ks)%deposition = sedvolcur(ks)%deposition + vol
        else
          sedvolcur(0)%erosion = sedvolcur(0)%erosion + vol
          sedvolcur(ks)%erosion = sedvolcur(ks)%erosion + vol
        endif
      enddo
      !Storage  
      fac=iwet(i)*areap(i)/rhosed
      do ks=1,nsed
        if(ntsch==1)then  
          vol = fac*(Ctk(i,ks)*h(i)/btk(i,ks) &
                   -Ctk1(i,ks)*h1(i)/btk1(i,ks))                         
        else
          vol = fac*(ctsch*Ctk(i,ks)*h(i)/btk(i,ks) &
                   -ctsch1*Ctk1(i,ks)*h1(i)/btk1(i,ks) &
                   +ctsch2*Ctk2(i,ks)*h2(i)/btk2(i,ks))    
        endif
        sedvolcur(0)%storage = sedvolcur(0)%storage + vol
        sedvolcur(ks)%storage = sedvolcur(ks)%storage + vol
        vol=fac*Ctk(i,ks)*h(i)/btk(i,ks) !Initial concentration assumed to be zero*****
        sedvolnet(0)%storage = sedvolnet(0)%storage + vol
        sedvolnet(ks)%storage = sedvolnet(ks)%storage + vol
      enddo
    enddo
   
    !All Boundaries
    do ibnd=1,nbndstr
      do j=1,bnd_str(ibnd)%ncells
        i=bnd_str(ibnd)%cells(j)
        k=bnd_str(ibnd)%faces(j)
        nck=cell2cell(k,i)
        vel=(fnx(k,i)*u(i)+fny(k,i)*v(i)) 
        fac=dtime*ds(k,i)*h(i)/rhosed
        do ks=1,nsed
          vol = fac*vel*Ctk(nck,ks)
          sedvolcur(0)%boundary = sedvolcur(0)%boundary + vol
          sedvolcur(ks)%boundary = sedvolcur(ks)%boundary + vol
        enddo
      enddo
    enddo        
    sedvolcum%boundary   = sedvolcum%boundary   + sedvolcur%boundary
    sedvolcum%storage    = sedvolcum%storage    + sedvolcur%storage
    sedvolcum%bedchange  = sedvolcum%bedchange  + sedvolcur%bedchange
    sedvolcum%erosion    = sedvolcum%erosion    + sedvolcur%erosion
    sedvolcum%deposition = sedvolcum%deposition + sedvolcur%deposition
    
    !Transport equation -> storage+transport=bedchange
    sedvolcur%balance = sedvolcur%storage + sedvolcur%boundary + sedvolcur%bedchange/max(scalemorph,1.0)  !Current time step
    sedvolcum%balance = sedvolcum%storage + sedvolcum%boundary + sedvolcum%bedchange/max(scalemorph,1.0)  !Cummulative
    !Note: sedvolnet%boundary cannot be calculated so sedvolcum%boundary is used
    sedvolnet%balance = sedvolnet%storage + sedvolcum%boundary + sedvolnet%bedchange/max(scalemorph,1.0)  !Net

    !!!Force mass balance
    !!volabschg = 0.0
    !!do i=1,ncells
    !!  vol = solid*dzb(i)*areap(i)
    !!  !vol = solid*(zb(i)-zb0(i))*areap(i)
    !!  volabschg = volabschg + abs(vol)
    !!enddo
    !!volabschg = max(volabschg,1.0e-20)
    !!do i=1,ncells
    !!  vol = solid*dzb(i)*areap(i)
    !!  !vol = solid*(zb(i)-zb0(i))*areap(i)
    !!  fac = abs(vol)/volabschg !Weight 0-1
    !!  zb(i) = zb(i) + edvolcur%balance*fac/(solid*areap(i))
    !!enddo
    !!!Recalculate Cumulative deposition and erosion volumes
    !!sedvol%bedchange = 0.0
    !!do i=1,ncells
    !!  vol = solid*(zb(i)-zb1(i))*areap(i)
    !!  sedvol%bedchange  =  sedvol%bedchange + vol
    !!enddo
    !!!Recalculate balance
    !!sedvol%balance = sedvol%storage + sedvol%boundary + sedvol%bedchange  !Current time step, total sediment
    
    !Percentage error
    do ks=0,nsed
      !Current
      sedvolcur(ks)%gross = sedvolcur(ks)%deposition - sedvolcur(ks)%erosion &
           + abs(sedvolcur(ks)%storage) + abs(sedvolcur(ks)%boundary) !normalizing volume
      sedvolcur(ks)%error = 100.0*sedvolcur(ks)%balance/max(sedvolcur(ks)%gross,1.0e-25)   
      !Cummulative
      sedvolcum(ks)%gross = sedvolcum(ks)%deposition - sedvolcum(ks)%erosion &
           + abs(sedvolcum(ks)%storage) + abs(sedvolcum(ks)%boundary) !normalizing volume
      sedvolcum(ks)%error = 100.0*sedvolcum(ks)%balance/max(sedvolcum(ks)%gross,1.0e-25)   
      !Net
      sedvolnet(ks)%gross = sedvolnet(ks)%deposition - sedvolnet(ks)%erosion &
           + abs(sedvolnet(ks)%storage) + abs(sedvolcum(ks)%boundary) !normalizing volume
      sedvolnet(ks)%error = 100.0*sedvolnet(ks)%balance/max(sedvolnet(ks)%gross,1.0e-25)
    enddo
                            
    !Print results               
751 format(' Sediment Balance Volumes, m^3 (Including Morphologic Acceleration Factor)')
752 format('                 Storage, Net Efflux, Erosion, Deposition, Balance, % Error ')
753 format(' Current       ',6(1x,1pe9.2),F5.2)
750 format(' Current_',I2.2,'    ',6(1x,1pe9.2),F5.2)
754 format(' Cumulative    ',6(1x,1pe9.2),F5.2)
755 format(' Cummulative_',I2.2,6(1x,1pe9.2),F5.2)
757 format(' Net           ',6(1x,1pe9.2),F5.2)
!844 format(' Current Sediment % Mass Error: ',1x,1pe9.2)    
!840 format(' Current Sediment',I2.2,' % Mass Error: ',1x,1pe9.2)   
!855 format(' Cumulative Sediment % Mass Error: ',1x,1pe9.2)
    
    iunit = (/6,dgunit/)
    open(dgunit,file=dgfile,access='append')  
    do i=1,2
      write(iunit(i),751)    
      write(iunit(i),752)
      write(iunit(i),753) sedvolcur(0)%storage,sedvolcur(0)%boundary,&
         sedvolcur(0)%erosion,sedvolcur(0)%deposition,sedvolcur(0)%balance,sedvolcur(0)%error
      if(nsed>1)then
        do ks=1,nsed
          write(iunit(i),750) ks,sedvolcur(ks)%storage,sedvolcur(ks)%boundary,&
           sedvolcur(ks)%erosion,sedvolcur(ks)%deposition,sedvolcur(ks)%balance,sedvolcur(ks)%error
        enddo
      endif
      write(iunit(i),754) sedvolcum(0)%storage,sedvolcum(0)%boundary,&
         sedvolcum(0)%erosion,sedvolcum(0)%deposition,sedvolcum(0)%balance,sedvolcum(0)%error
      if(nsed>1)then
        do ks=1,nsed
          write(iunit(i),755) ks,sedvolcum(ks)%storage,sedvolcum(ks)%boundary,&
           sedvolcum(ks)%erosion,sedvolcum(ks)%deposition,sedvolcum(ks)%balance,sedvolcum(ks)%error
        enddo
      endif
      write(iunit(i),757) sedvolnet(0)%storage,sedvolcum(0)%boundary,&
         sedvolnet(0)%erosion,sedvolnet(0)%deposition,sedvolnet(0)%balance,sedvolnet(0)%error
    enddo
    close(dgunit)
      
    return
    endsubroutine sed_balance
    
!!***************************************
!    subroutine fallvel_cohsed
!! Cohesive sediment fall velocity formula -  By Wu    
!!**************************************
!    use case_size
!    use comvarbl, only: viscos
!    use fl2d, only: bsxy
!    use sedmod   
!    use salmod, only: saltrans, sal
!    use precision
!    implicit none
!    integer :: i
!    !real(ikind) :: d(n),ds(n),ws(n)
!    real(ikind) :: correcsed, cohk,  correcsal, correcturb
! 
!    !wsfallmax=wsfall(1)*(0.000022/diam(1))**1.8*cohk1*cohcp**cohn*cohsalkmax*   &
!    !                           (1.0+cohturbk1*cohtaubp**cohturbn1)
!    !write(*,*) "Cohesive sediment wsfallmax=", wsfallmax
!
!    do i=1,ncells
!       wsfallcohsed(i)=wsfall(1)*(0.000022/diam(1))**1.8
!
!       !cohk1=wsfallmax/wsfallcohsed(i)/cohcp**cohn/cohsalkmax/   &
!       !               (1.0+cohturbk1*cohtaubp**cohturbn1)
!       if(Ctk(i,1).lt.cohcp) then ! Effect of Sed. Concentration on Settling
!          correcsed=cohk1*Ctk(i,1)**cohn
!       else
!          cohk=cohk1*cohcp**cohn/(1.0-cohk2*cohcp)**cohr
!          correcsed=cohk*(max(0.0,1.0-cohk2*Ctk(i,1)))**cohr
!       endif
!       correcsed=max(correcsed, wsfall(1)/wsfallcohsed(i))
!    
!       if(saltrans) then   ! Effect of Salinity on Settling
!          if(sal(i).lt.cohsalcp) then
!             correcsal=1.0+(cohsalkmax-1.0)*(sal(i)/cohsalcp)**cohsaln
!          else
!             correcsal=cohsalkmax
!          endif   
!       else            
!          correcsal=1.0
!       endif
!
!       if(cohbsxy(i).lt.cohtaubp) then   !effect of turbulence on settling
!          correcturb=1.0+cohturbk1*cohbsxy(i)**cohturbn1
!       else
!          correcturb=(1.0+cohturbk1*cohtaubp**cohturbn1)  &
!                         *(cohbsxy(i)/cohtaubp)**(-cohturbn2)
!       endif               
!            
!       wsfallcohsed(i)=wsfallcohsed(i)*correcsed*correcsal*correcturb
!       wsfallcohsed(i)=max(wsfallcohsed(i), wsfall(1))
!       !wsfallcohsed(i)=max(wsfallcohsed(i), 0.0001*wsfall(1))
!    enddo
!
!    return
!    endsubroutine 

!*******************************************************************    
    subroutine sed_concdepthchange
! Author: Alex Sanchez, USACE-CHL
!*******************************************************************
    use sed_def
    use flow_def
    use geo_def
    use size_def
    implicit none
    integer :: i
    real(ikind) ::val
    
!$OMP PARALLEL DO PRIVATE(i,val)
    do i=1,ncells
      val=h(i) !save old total water depth
      h(i)=max(hdry,p(i)*gravinv-zb(i))  !New total water depth
      !if(h(i)>2*hdry .and. val>2*hdry .and. abs(dzb(i))<0.25*hdry)then      
      if(h(i)>2*hdry .and. val>2*hdry)then     
        Ctk(i,:) = val*Ctk(i,:)/h(i)
      endif  
    enddo
!$OMP END PARALLEL DO

    return
    endsubroutine sed_concdepthchange
