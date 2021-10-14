!==========================================================================
! CMS Output routines
!
! Contains:
!   out_default - Sets the output module defaults
!   out_cards - Reads the output module cards from the control file
!   out_init - Initializes the output variables
!   out_print - Prints the output module settings to the screen and
!              the diagnostic file
!   write_output - writes the global variables to the solution files
!
! written by Alex Sanchez, USACE-CHL
!         Mitch Brown, USACE-CHL
!         Chris Reed, URS, R&R
!==========================================================================
    
!***************************************************************************   
    subroutine out_default
! Sets the output module defaults
! written by Alex Sanchez
!***************************************************************************        
#include "CMS_cpp.h"
    use out_def
    use comvarbl, only: casename,flowpath!,input_ver
    implicit none
    
    obs_cell = .false.           !Observation cells
    save_point = .false.         !Save point cells specified
    
    !Geometry
    write_areap     = .false.    !Turns on or off the current velocity magnitude = uv
    write_cell2cell = .false.    !Turns on or off the forward cell to cell connectivity
    write_llec2llec = .false.    !Turns on or off the backwards cell to cell connectivity
    write_cell2node = .false.    !Turns on or off the forward cell to node connectivity    
    write_node2cell = .false.    !Turns on or off the forward node to cell connectivity    
    write_geocells  = .false.    !Turns on or off the geometry file for all cells    
    
    !Hydro
    write_velmag = .true.        !Current velocity magnitude = uv
    write_veltotal = .false.     !Total flux velocity
    write_veltotalmag = .false.  !Total flux velocity magnitude
    write_pres = .false.         !WSE pressure = p
    write_totdep = .false.       !Total depth = h
    write_hpred = .false.        !Predicted total water depth
    write_wsepred = .false.      !Predicted water level
    write_velpred = .false.      !Predicted current velocity
    write_wsegrad = .false.      !WSE gradients
    write_velgrad = .false.      !Velocity gradients    
    write_prescor = .false.      !Pressure correction = pp
    write_prescorgrad = .false.  !Pressure correction gradients = (dppx,dppy)
    write_presgrad = .false.     !Pressure gradients = (dpx,dpy)
    write_presres = .false.      !Normalized residuals for p
    write_velres = .false.       !Normalized residuals for u,v
    write_maxwse = .false.       !Maximum water level output 
    
    !Sediment
    write_conc = .true.          !Concentration
    write_capac = .true.         !Concentration capacity
    write_fracsusp = .false.     !Fraction of suspended transport
    write_morph = .true.         !Morphology Change (bed change) 
    write_alphatot = .false.     !Adaptation coefficient for total-load
    write_lengthtot = .false.    !Adaptation length for total-load
    write_lengthsusp = .false.   !Adaptation length for suspended-load
    write_lengthbed = .false.    !Adaptation length for bed-load
    write_betatot = .false.      !Correction factor for suspended-load
    write_betasusp = .false.     !Correction factor for bed-load
    write_bedvel = .false.       !Bed-load velocity
       
    !Wave
    write_wave_detail = .false.  !Wave details
    write_wavstress = .false.    !Wave stress vectors
    write_wavbrkdiss = .false.   !Wave dissipation 
    write_wavbrkind = .false.    !Wave breaking index (0 or 1)
    write_fracbreak = .false.    !Fraction of breaking waves
    write_wavdisscoef = .false.  !Wave dissipation coefficient
    write_wavorbvel = .false.    !Wave bottom orbital velocity
    write_wavlen = .false.       !Wave length
    write_wavvel = .false.       !Wave mass transport velocity vector
    write_wavvelmag = .false.    !Wave mass transport velocity magnitude
    
    !Roller
    write_rollenergy = .false.   !Roller energy
    write_rollstress = .false.   !Roller stress   
    write_rolldiss = .false.     !Roller dissipation
     
    !Meteorological
    write_wndmag = .false.       !Wind velocity magnitude [m/s]
    write_wndstress = .false.    !Wind stress vector    [Pa]
    write_wndstressmag = .false. !Wind stress magnitude [Pa]
    write_atmpres  = .false.     !Atmospheric pressure  [Pa]
    write_atmpresgrad = .false.  !Atmospheric pressure gradient [Pa/m]
    
    !Roughness
    write_meanbedshear = .false.
    write_normrough = .false.
    write_normapprough = .false.
    
    !Compression
    ixmdfcomp = 0                !0-No XMDF compression, 1-XMDF Compression    
    
    !Input/Output
    write_ascii_input = .false.  !Turns on or off writing the ASCII input files
    write_tecplot     = .false.  !Tecplot files 
    write_netcdf      = .false.  !NetCDF files   !Just Preparing to add this feature  MEB 10/03/2019

#ifdef _WIN32
    write_xmdf_output = .true.   !By default, turns on writing XMDF output files on Windows
    write_sup = .false.          !  and turns off Super ASCII files
#else    
    write_xmdf_output = .false.  !By default, turns off writing XMDF output files on Linux
    write_sup = .true.           !  and turns on Super ASCII files
#endif
    
    !Simulation label
    simlabel = 'Simulation'
    
    !--- File names ------------------------------------------
    outpath = ''
    if(len_trim(flowpath)>0)then
      outpath = flowpath
    endif
    
    !--- Save points --------------
    maxout = 1000
    
    !Flow/Hydro
    savept(1)%group = 'HYDRO'
    savept(1)%time_inc = 0.0   !1 is Hydro group
    savept(1)%ncells = 0
    savept(1)%active = .false.   !make sure there are save points listed to make true.
    savept(1)%nvar = 4
    allocate(savept(1)%names(savept(1)%nvar))
    allocate(savept(1)%ounit(savept(1)%nvar))
    allocate(savept(1)%vals(savept(1)%nvar))
    savept(1)%names(1) = 'eta' ;   savept(1)%ounit(1) = 'm'       ; savept(1)%vals(1) = 1
    savept(1)%names(2) = 'uv'  ;   savept(1)%ounit(2) = 'm/sec'   ; savept(1)%vals(2) = 2
    savept(1)%names(3) = 'flux';   savept(1)%ounit(3) = 'm^2/sec' ; savept(1)%vals(3) = 2
    savept(1)%names(4) = 'totdep'; savept(1)%ounit(4) = 'm'       ; savept(1)%vals(4) = 1
    
    !Sediment transport
    savept(2)%group = 'SEDIMENT'
    savept(2)%time_inc = 0.0                   !2 is Sediment group
    savept(2)%ncells = 0
    savept(2)%active = .false.                     !make sure there are save points listed to make true.
    savept(2)%nvar = 4
    allocate(savept(2)%names(savept(2)%nvar))
    allocate(savept(2)%ounit(savept(2)%nvar))
    allocate(savept(2)%vals(savept(2)%nvar))
    savept(2)%names(1) = 'susp_conc'; savept(2)%ounit(1) = 'kg/m^3'  ; savept(2)%vals(1) = 1
    savept(2)%names(2) = 'trans'    ; savept(2)%ounit(2) = 'kg/m/sec'; savept(2)%vals(2) = 2
    savept(2)%names(3) = 'susp_epsv'; savept(2)%ounit(3) = 'm^2/sec' ; savept(2)%vals(3) = 1 !Only for advanced users
    savept(2)%names(4) = 'susp_beta'; savept(2)%ounit(4) = '-' ;       savept(2)%vals(4) = 1 !Only for advanced users
    
    !Salinity
    savept(3)%group = 'SALINITY'
    savept(3)%time_inc = 0.0                   !3 is Salinity group
    savept(3)%ncells = 0
    savept(3)%active = .false.                     !make sure there are save points listed to make true.
    savept(3)%nvar = 1
    allocate(savept(3)%names(1),savept(3)%ounit(1),savept(3)%vals(1))
    savept(3)%names(1) = 'salt' ; savept(3)%ounit(1) = 'ppt' ; savept(3)%vals(1) = 1
    
    !Waves
    savept(4)%group = 'WAVES'
    savept(4)%time_inc = 0.0                   !4 is Waves group
    savept(4)%ncells = 0
    savept(4)%active = .false.                     !make sure there are save points listed to make true.
    savept(4)%nvar = 4
    allocate(savept(4)%names(savept(4)%nvar))
    allocate(savept(4)%ounit(savept(4)%nvar))
    allocate(savept(4)%vals(savept(4)%nvar))
    savept(4)%names(1) = 'height'  ; savept(4)%ounit(1) = 'm'      ; savept(4)%vals(1) = 1
    savept(4)%names(2) = 'period'  ; savept(4)%ounit(2) = 'sec'    ; savept(4)%vals(2) = 1
    savept(4)%names(3) = 'hgt-vec' ; savept(4)%ounit(3) = 'm'      ; savept(4)%vals(3) = 2
    savept(4)%names(4) = 'diss'    ; savept(4)%ounit(4) = 'N/m/s'  ; savept(4)%vals(4) = 1
    
    !Morphology
    savept(5)%group = 'Morhology'
    savept(5)%time_inc = 0.0                   !5 is Morphology group
    savept(5)%ncells = 0
    savept(5)%active = .false.                     !make sure there are save points listed to make true.
    savept(5)%nvar = 2
    allocate(savept(5)%names(2),savept(5)%ounit(2),savept(5)%vals(2))
    savept(5)%names(1) = 'depth' ; savept(1)%ounit(1) = 'm' ; savept(5)%vals(1) = 1
    savept(5)%names(2) = 'morph' ; savept(2)%ounit(2) = 'm' ; savept(5)%vals(2) = 1
    
    return
    endsubroutine out_default

!***************************************************************************   
    subroutine out_cards(cardname,foundcard)
! Sets the default parameters for flow variables
! written by Alex Sanchez, USACE-CHL; Weiming Wu, NCCHE
!***************************************************************************
    use out_def
    use prec_def
    use comvarbl, only: etime
    use diag_lib, only: diag_print_warning
    implicit none
    integer :: i,j,cell,ierr
    character(len=37) :: cardname,cdum,ptname
    character(len=37) :: group(ngroups)
    character(len=120) :: saveargs,astring
    real(ikind) :: xsave,ysave
    logical :: foundcard
    
    interface
      function toUpper (astr)
        character(len=*),intent(in) :: astr
        character(len=len(astr)) :: toUpper
      end function
    end interface    
    
    do i=1,ngroups
      group(i) = ''
    enddo
    foundcard = .true.
    selectcase(cardname)    
      case('SIMULATION_LABEL')
        backspace(77)
        read(77,*) cardname, simlabel    
    
      case('OUTPUT_PATH')
        backspace(77)
        read(77,*) cardname, outpath
        
      !=== Time Series ====================  
      case('TIME_SERIES_1')
        backspace(77)
        read(77,*) cardname, outseries(1)%ntimes
        backspace(77)
        if(allocated(outseries(1)%times))then
          deallocate(outseries(1)%times) !override previous read
        endif          
        allocate(outseries(1)%times(outseries(1)%ntimes))
        read(77,*) cardname, outseries(1)%ntimes, (outseries(1)%times(j),j=1,outseries(1)%ntimes)          
        
      case('TIME_SERIES_2')
        backspace(77)
        read(77,*) cardname, outseries(2)%ntimes
        backspace(77)
        if(allocated(outseries(2)%times))then
          deallocate(outseries(2)%times) !override previous read
        endif  
        allocate(outseries(2)%times(outseries(2)%ntimes))
        read(77,*) cardname, outseries(2)%ntimes, (outseries(2)%times(j),j=1,outseries(2)%ntimes)

      case('TIME_SERIES_3')
        backspace(77)
        read(77,*) cardname, outseries(3)%ntimes
        backspace(77)
        if(allocated(outseries(3)%times))then
          deallocate(outseries(3)%times) !override previous read
        endif  
        allocate(outseries(3)%times(outseries(3)%ntimes))
        read(77,*) cardname, outseries(3)%ntimes, (outseries(3)%times(j),j=1,outseries(3)%ntimes)        

      case('TIME_SERIES_4')
        backspace(77)
        read(77,*) cardname, outseries(4)%ntimes
        backspace(77)
        if(allocated(outseries(4)%times))then
          deallocate(outseries(4)%times) !override previous read
        endif  
        allocate(outseries(4)%times(outseries(4)%ntimes))
        read(77,*) cardname, outseries(4)%ntimes, (outseries(4)%times(j),j=1,outseries(4)%ntimes)

      !==== Time Lists ==================================
      case('TIME_LIST_1')      !Alex, April 16, 2010    
        call read_time_list(1) 

      case('TIME_LIST_2')      !Alex, April 16, 2010      
        call read_time_list(2)
        
      case('TIME_LIST_3')      !Alex, April 16, 2010       
        call read_time_list(3)  
          
      case('TIME_LIST_4')      !Alex, April 16, 2010      
        call read_time_list(4)
      
      !==== Out Times List ===============================
      case('WSE_OUT_TIMES_LIST')
        backspace(77)
        read(77,*) cardname, outlist(1)%ilist
      
      case('VEL_OUT_TIMES_LIST')
        backspace(77)
        read(77,*) cardname, outlist(2)%ilist
          
      case('EDDY_OUT_TIMES_LIST','VISC_OUT_TIMES_LIST','EDDY_VISCOSITY_OUT_TIMES_LIST')
        backspace(77)
        read(77,*) cardname, outlist(3)%ilist
        
      case('SED_OUT_TIMES_LIST')
        backspace(77)
        read(77,*) cardname, outlist(4)%ilist
      
      case('MORPH_OUT_TIMES_LIST')
        backspace(77)
        read(77,*) cardname, outlist(5)%ilist
        
      case('TRANS_OUT_TIMES_LIST')
        backspace(77)
        read(77,*) cardname, outlist(9)%ilist
    
      case('BED_OUT_TIMES_LIST','BEDCOMP_OUT_TIMES_LIST')
        backspace(77)
        read(77,*) cardname, outlist(10)%ilist
        
      case('SEDMIX_OUT_TIMES_LIST','FRAC_OUT_TIMES_LIST')
        backspace(77)
        read(77,*) cardname, outlist(11)%ilist
          
      case('WAVE_OUT_TIMES_LIST','WAVES_OUT_TIMES_LIST')
        backspace(77)
        read(77,*) cardname, outlist(7)%ilist
      
      case('WAVE_DETAILS_OUT_TIMES_LIST')
        backspace(77)
        read(77,*) cardname, outlist(12)%ilist
        
      case('STRESS_OUT_TIMES_LIST','BED_SHEAR_STRESS_OUT_TIMES_LIST',&
           'ROUGHNESS_OUT_TIMES_LIST','ROUGH_OUT_TIMES_LIST','ROUGHNESS_STRESS_OUT_TIMES_LIST',&
           'FRICTION_OUT_TIMES_LIST','FRIC_OUT_TIMES_LIST')
        backspace(77)
        read(77,*) cardname, outlist(13)%ilist
                        
      !=== Observation Points =========================
      !Time series, u, v, and eta
      case('TIME_SERIES_INCREMENT') !Old
        backspace(77)        
        read(77,*) cardname, obs(1)%time_inc
      
      case('ELEV_OBS_CELLS_BEGIN','HYDRO_OBS_CELLS_BEGIN','TIME_SERIES_OBS_CELLS_BEGIN') !card should be changed just TIME_SERIES_OBS_CELLS_BEGIN
        i=1
        obs(i)%group = 'Time_Series'
        obs(i)%active = .true.
        obs(i)%nvar = 4 !u,v,eta,velpro
        allocate(obs(i)%names(obs(i)%nvar))
        obs(i)%names(1)='u'
        obs(i)%names(2)='v'
        obs(i)%names(3)='eta'
        obs(i)%names(4)='velpro'
        call read_obs_cells(i)
      
      case('ELEV_OBS_CELLS_END','END') !card should be changed just TIME_SERIES_OBS_CELLS_BEGIN
        !DO NOTHING
        
      !Flow, u*h, v*h
      case('FLOW_RATE_INCREMENT') !Old
        backspace(77)
        read(77,*) cardname, obs(2)%time_inc
        
      case('FLOW_OBS_CELLS_BEGIN')
        i=2        
        obs(i)%active = .true.
        obs(i)%group = 'Flow_Rate'
        obs(i)%nvar = 2 !u*h,v*h
        allocate(obs(i)%names(obs(i)%nvar))
        obs(i)%names(1)='Qx'
        obs(i)%names(2)='Qy'  
        call read_obs_cells(i)  
      
      case('FLOW_OBS_CELLS_END')
            
      !Transport, sediment and salinity
      case('Q_TRANS_RATE_INCREMENT') !Old
        backspace(77)
        read(77,*) cardname, obs(3)%time_inc
        
      case('Q_TRANS_OBS_CELLS_BEGIN')
        i=3
        obs(i)%group = 'Transport'
        obs(i)%active = .true.
        obs(i)%nvar = 5 !u*h,v*h
        allocate(obs(i)%names(obs(i)%nvar))
        obs(i)%names(1)='Csus'
        obs(i)%names(2)='qtx'
        obs(i)%names(3)='qty'
        obs(i)%names(4)='sal'
        obs(i)%names(5)='Csuspro'
        call read_obs_cells(i)
      
      case('Q_TRANS_OBS_CELLS_END')
        
      !Transport, sediment and salinity
      case('BED_COMPOSITION_INCREMENT')
        backspace(77)
        read(77,*) cardname, obs(4)%time_inc
        
      case('BED_OBS_CELLS_BEGIN')
        i=4
        obs(i)%group = 'Bed_Composition'
        obs(i)%active = .true.
        obs(i)%nvar = 1 !pbk(:,:,1)
        allocate(obs(i)%names(obs(i)%nvar))
        obs(i)%names(1)='bed'
        call read_obs_cells(i)
      
      case('BED_OBS_CELLS_END')
        !DO NOTHING  
        
      !=== Output Files and Options ==========================
      case('XMDF_COMPRESSION')
        backspace(77)
        read(77,*) cardname, cdum
        if(cdum(1:3)=='ON ')then
          ixmdfcomp = 1
        elseif(cdum(1:3)=='OFF')then
          ixmdfcomp = 0
        endif 
        
      case('ELAPSED_TIME_INTERVAL')   !added so that the Elapsed time output can be defined by user
        call card_scalar(77,'sec','sec',etime,ierr)
        
      case('GLOBAL_TECPLOT_FILES')
        call card_boolean(77,write_tecplot,ierr)
         
      case('GLOBAL_SUPER_FILES','GLOBAL_SUPER_ASCII')
        call card_boolean(77,write_sup,ierr)       
          
      case('WRITE_ASCII_INPUT_FILES','OUTPUT_ASCII_INPUT')
        call card_boolean(77,write_ascii_input,ierr)
          
      case('OUTPUT_WAVE_DETAILS','WAVE_OUTPUT_DETAILS')
        call card_boolean(77,write_wave_detail,ierr)
        
      case('USE_COMMON_SOLUTION_FILE')
        call card_boolean(77,use_common,ierr)
      
      case('GLOBAL_WATER_LEVEL_OUTPUT','WSE_OUT_FILE','WSE_SOL_FILE')
        backspace(77)
        read(77,*) cardname,outlist(1)%afile !,outlist(1)%apath
        outlist(1)%use_default=.false.
        !outlist(3)%afile = outlist(1)%afile                          !This was temporary and no longer necessary, MEB  02/16/2016
      
      case('GLOBAL_VELOCITY_OUTPUT','VEL_OUT_FILE','VEL_SOL_FILE')
        backspace(77)
        outlist(2)%use_default=.false.
        read(77,*) cardname,outlist(2)%afile !,outlist(2)%apath
      
      case('GLOBAL_EDDY_VISCOSITY_OUTPUT','EDDY_VISCOSITY_OUT_FILE',&
           'EDDY_OUT_FILE','VISC_OUT_FILE','VISC_SOL_FILE')
        backspace(77)
        read(77,*) cardname,outlist(3)%afile !,outlist(3)%apath
        outlist(3)%use_default=.false.
      
      case('GLBL_CONC_CAPAC_OUTPUT')   !Old format
!        backspace(77)
!        read(77,*) cardname,outlist(4)%afile !,outlist(4)%apath
        
      case('GLBL_CONCENTRATION_OUTPUT') !Old format
!        backspace(77)
!        read(77,*) cardname,outlist(4)%afile !,outlist(4)%apath

      case('SED_OUT_FILE','SED_SOL_FILE','SEDIMENT_OUTPUT_FILE') !New format
        backspace(77)
        read(77,*) cardname,outlist(4)%afile !,outlist(4)%apath      
        outlist(4)%use_default=.false.
      
      case('GLOBAL_MORPHOLOGY_OUTPUT','MORPH_OUT_FILE','MORPHOLOGY_OUTPUT_FILE')   !New
        backspace(77)
        read(77,*) cardname,outlist(5)%afile !,outlist(5)%apath      
        outlist(5)%use_default=.false.
      
      case('GLOBAL_SALINITY_OUTPUT','SAL_OUT_FILE','SALT_OUT_FILE','SALINITY_OUTPUT_FILE')
        backspace(77)
        read(77,*) cardname,outlist(6)%afile !,outlist(6)%apath      
        outlist(6)%use_default=.false.
      
      case('GLOBAL_TEMPERATURE_OUTPUT','HEAT_OUT_FILE','TEMPERATURE_OUTPUT_FILE')
        backspace(77)
        read(77,*) cardname,outlist(14)%afile !,outlist(14)%apath      
        outlist(14)%use_default=.false.
      
      case('GLOBAL_WAVES_OUTPUT','GLOBAL_WAVE_OUTPUT','WAVE_OUT_FILE','WAVES_OUT_FILE','WAVE_SOL_FILE')
        backspace(77)
        read(77,*) cardname,outlist(7)%afile !,outlist(7)%apath      
        outlist(7)%use_default=.false.
        
      case('GLOBAL_WIND_OUTPUT','WIND_OUT_FILE','WIND_SOL_FILE','WIND_OUTPUT_FILE')
        backspace(77)
        read(77,*) cardname,outlist(8)%afile !,outlist(8)%apath        
        outlist(8)%use_default=.false.
      
      case('GLOBAL_TRANS_RATE_OUTPUT','TRANS_OUT_FILE','TRANSPORT_OUTPUT_FILE')
        backspace(77)
        read(77,*) cardname,outlist(9)%afile !,outlist(9)%apath
        outlist(9)%use_default=.false.
      
      case('BED_OUT_FILE','BED_SOL_FILE','BEDCOMP_OUT_FILE','BEDCOMP_SOL_FILE')
        backspace(77)
        read(77,*) cardname,outlist(10)%afile !,outlist(10)%apath
        outlist(10)%use_default=.false.
      
      case('SEDMIX_OUT_FILE')
        backspace(77)
        read(77,*) cardname,outlist(11)%afile !,outlist(11)%apath
        outlist(11)%use_default=.false.
      
      case('WAVEADV_OUT_FILE','WAVEADV_OUTPUT_FILE')
        backspace(77)
        read(77,*) cardname,outlist(12)%afile !,outlist(12)%apath'
        outlist(12)%use_default=.false.
        
      case('BED_SHEAR_STRESS_OUT_FILE','SHEAR_OUT_FILE',&
           'ROUGHNESS_OUT_FILE','ROUGHNESS_OUTPUT_FILE')
        backspace(77)
        read(77,*) cardname,outlist(13)%afile !,outlist(13)%apath
        outlist(13)%use_default=.false.

      case('U_TIME_SERIES_OUTPUT')
      
      case('V_TIME_SERIES_OUTPUT')
      
      case('ETA_TIME_SERIES_OUTPUT')
      
      case('X_FLOW_RATE_OUTPUT')
      
      case('Y_FLOW_RATE_OUTPUT')
      
      case('QX_TRANS_RATE_OUTPUT')
      
      case('QY_TRANS_RATE_OUTPUT')

      !---- Booleans ----------------------------------    
      !Geometry    
      case('CELL_AREA_OUTPUT','OUTPUT_CELL_AREA','OUTPUT_AREA','AREA_OUTPUT')
        call card_boolean(77,write_areap,ierr)
        
      case('CELL_TO_CELL_OUTPUT','OUTPUT_CELL_TO_CELL')
        call card_boolean(77,write_cell2cell,ierr)
      
      case('CELL_FROM_CELL_OUTPUT','OUTPUT_CELL_FROM_CELL')
        call card_boolean(77,write_llec2llec,ierr)
        
      case('CELL_TO_NODE_OUTPUT','OUTPUT_CELL_TO_NODE')
        call card_boolean(77,write_cell2node,ierr)
        
      case('NODE_TO_CELL_OUTPUT','OUTPUT_NODE_TO_CELL')
        call card_boolean(77,write_node2cell,ierr)  
        
      case('GEOMETRY_OUTPUT','OUTPUT_GEOMETRY','GEO_CELLS_OUTPUT','OUTPUT_GEO_CELLS')
        call card_boolean(77,write_geocells,ierr)
        
      case('OUTPUT_GLOBAL_XMDF','GLOBAL_XMDF','GLOBAL_XMDF_OUTPUT')
        call card_boolean(77,write_xmdf_output,ierr)  
        write_sup = .false.
        write_xmdf_output = .true.
        
      case('WRITE_ASCII_OUTPUT_FILES')
        write_xmdf_output = .false.
        write_sup = .true.

      case('OUTPUT_FILE_TYPE')
        backspace(77)
        read(77,*) cardname,astring
        if    (toUpper(astring(1:4)) == "XMDF")  then   !XMDF, turn off other types
          write_xmdf_output = .true.
          write_sup         = .false.
        elseif(toUpper(astring(1:5)) == 'ASCII') then   !ASCII, turn off other types
          write_sup         = .true.
          write_xmdf_output = .false.
        elseif(toUpper(astring(1:4)) == 'BOTH')  then   !Other, allow for both XMDF and ASCII
          write_sup         = .true.
          write_xmdf_output = .true.
        else
          write_xmdf_output = .true.
          write_sup         = .false.
          call diag_print_warning('Unrecognized value for card, OUTPUT_FILE_TYPE','Defaulting to XMDF')
        endif
        
      !Hydrodynamics  
      case('VEL_MAG_OUTPUT','OUTPUT_VELOCITY_MAGNITUDE','VELOCITY_MAGNITUDE_OUTPUT')
        call card_boolean(77,write_velmag,ierr)
        
      case('TOTAL_VEL_OUTPUT','OUTPUT_TOTAL_VELOCITY',&
        'TOTAL_FLUX_VELOCITY_OUTPUT','TOTAL_FLUX_VEL_OUTPUT')
        call card_boolean(77,write_veltotal,ierr)  
        
      case('TOTAL_VEL_MAG_OUTPUT','OUTPUT_TOTAL_VELOCITY_MAGNITUDE',&
         'TOTAL_FLUX_VELOCITY_MAG_OUTPUT','TOTAL_FLUX_VEL_MAG_OUTPUT')
        call card_boolean(77,write_veltotalmag,ierr)   
      
      case('PRESSURE_OUTPUT','WATER_PRESSURE_OUTPUT','OUTPUT_PRESSURE','PRES_OUTPUT','OUTPUT_WATER_PRESSURE')
        call card_boolean(77,write_pres,ierr)
        
      case('PRESSURE_GRADIENT_OUTPUT','WATER_PRES_GRAD_OUTPUT',&
           'OUTPUT_PRESSURE_GRADIENT','PRES_GRAD_OUTPUT')
        call card_boolean(77,write_presgrad,ierr)  
      
      case('TOTAL_DEPTH_OUTPUT','WATER_DEPTH_OUTPUT',&
           'TOTAL_WATER_DEPTH_OUTPUT','OUTPUT_WATER_DEPTH')
        call card_boolean(77,write_totdep,ierr)
        
      case('PRED_WATER_DEPTH_OUTPUT','PREDICTED_TOTAL_WATER_DEPTH_OUTPUT',&
          'PREDICTED_WATER_DEPTH_OUTPUT')
        call card_boolean(77,write_hpred,ierr)
      
      case('PREDICTED_WATER_LEVEL_OUTPUT','PREDICTED_WATER_ELEVATION_OUTPUT',&
           'PRED_WATER_LEVEL_OUTPUT','PREDICTED_WATER_SURFACE_ELEVATION_OUTPUT')
        call card_boolean(77,write_wsepred,ierr)
        
      case('PRED_CUR_VEL_OUTPUT','PREDICTED_CURRENT_VELOCITY_OUTPUT')
        call card_boolean(77,write_velpred,ierr)
      
      case('VELOCITY_CURVATURE_OUTPUT','OUTPUT_VELOCITY_CURVATURE')
        call card_boolean(77,write_velcurv,ierr)
      
      case('VELOCITY_GRADIENT_OUTPUT','OUTPUT_VELOCITY_GRADIENT')
        call card_boolean(77,write_velgrad,ierr)
        
      case('STREAMWISE_CURVATURE_OUTPUT','OUTPUT_STREAMWISE_CURVATURE')
        call card_boolean(77,write_streamcurv,ierr)  
        
      case('WATER_LEVEL_GRADIENT_OUTPUT','OUTPUT_WATER_LEVEL_GRADIENT')
        call card_boolean(77,write_wsegrad,ierr)
        
      case('PRESSURE_CORRECTION_OUTPUT','PRES_COR_OUTPUT',&
          'OUTPUT_PRES_COR','OUTPUT_PRESSURE_CORRECTION')
        call card_boolean(77,write_prescor,ierr)
        
      case('PRESSURE_CORRECTION_GRADIENT_OUTPUT','PRES_COR_GRAD_OUTPUT',&
          'OUTPUT_PRES_COR_GRAD','OUTPUT_PRESSURE_CORRECTION_GRADIENT')
        call card_boolean(77,write_prescorgrad,ierr)  
        
      case('VOL_RES_OUTPUT','VOLUME_RESIDUAL_OUTPUT','OUTPUT_VOLUME_RESIDUAL',&
           'VOLUME_BALANCE_OUTPUT','OUTPUT_VOLUME_BALANCE')
        call card_boolean(77,write_volres,ierr)  
        
      case('PRESSURE_RESIDUALS_OUTPUT','OUTPUT_PRESSURE_RESIDUALS')
        call card_boolean(77,write_presres,ierr)
        
      case('VELOCITY_RESIDUALS_OUTPUT','OUTPUT_VELOCITY_RESIDUALS')
        call card_boolean(77,write_velres,ierr)  
        
      case('MAXIMUM_WATER_LEVEL_OUTPUT')            !Added MEB 9/20/2021  Provide a mechanism to turn this on
        call card_boolean(77,write_maxwse,ierr)
        
      !Sediment Transport  
      case('CONC_OUTPUT','CONCENTRATION_OUTPUT','OUTPUT_CONC','OUTPUT_CONCENTRATION') !Ct
        call card_boolean(77,write_conc,ierr)
        
      case('CONC_FRAC_OUTPUT','OUTPUT_CONC_FRAC','OUTPUT_CONCENCTRATION_FRACTION') !Ctk
        call card_boolean(77,write_concfrac,ierr)  
        
      case('CONC_SUSP_OUTPUT','OUTPUT_CONC_SUSP','OUTPUT_CONCENCTRATION_SUSPENDED') !Cs
        call card_boolean(77,write_concsusp,ierr)
      
      !case('CONC_SUSP_FRAC_OUTPUT','OUTPUT_CONC_SUSP_FRAC') !Csk
      !  call card_boolean(77,write_concsuspfrac,ierr)
        
      case('CAPAC_OUTPUT','CAPACITY_OUTPUT','OUTPUT_CAPAC','OUTPUT_CAPACITY') !Ctstar
        call card_boolean(77,write_capac,ierr)
        
      case('MORPH_OUTPUT','OUTPUT_MORPH','OUTPUT_MORPHOLOGY_CHANGE') !dzb
        call card_boolean(77,write_morph,ierr)
      
      case('FRAC_SUSP_OUTPUT','SUSP_FRAC_OUTPUT','OUTPUT_FRAC_SUSP','FRACTION_SUSPENDED_OUTPUT') !rs
        call card_boolean(77,write_fracsusp,ierr)
        
      case('BED_FRAC_OUTPUT','BED_FRACTION_OUTPUT','OUTPUT_BED_FRAC','OUTPUT_BED_FRACTION',&
           'SIZE_FRACTION_OUTPUT','OUTPUT_SIZE_FRACTION','SED_FRAC_OUTPUT','OUTPUT_SED_FRACTION') !pbk
        call card_boolean(77,write_sizefrac,ierr)
    
      case('THICKNESS_OUTPUT','OUTPUT_THICKNESS',&
           'LAYER_THICKNESS_OUTPUT','OUTPUT_LAYER_THICKNESS') !db
        call card_boolean(77,write_thickness,ierr)
        
      case('ADAPTATION_COEFFICIENT_TOTAL_OUTPUT','OUTPUT_ADAPTATION_COEFFICIENT_TOTAL',&
           'ADAPTATION_TOTAL_OUTPUT','OUTPUT_ADAPTATION_TOTAL') !alphatot
        call card_boolean(77,write_alphatot,ierr)
        
      case('ADAPTATION_LENGTH_TOTAL_OUTPUT','OUTPUT_ADAPTATION_LENGTH_TOTAL')
        call card_boolean(77,write_lengthtot,ierr) !Ltot
        
      case('ADAPTATION_LENGTH_SUSPENDED_OUTPUT','OUTPUT_ADAPTATION_SUSPENDED_TOTAL')
        call card_boolean(77,write_lengthsusp,ierr) !Lsus
      
      case('ADAPTATION_LENGTH_BED_OUTPUT','OUTPUT_ADAPTATION_BED_TOTAL')
        call card_boolean(77,write_lengthbed,ierr) !Lsus
        
      case('SEDIMENT_CORRECTION_TOTAL_OUTPUT','OUTPUT_SEDIMENT_CORRECTION_TOTAL')
        call card_boolean(77,write_betatot,ierr) !Betatot
        
      case('SEDIMENT_CORRECTION_SUSPENDED_OUTPUT','OUTPUT_SEDIMENT_CORRECTION_SUSPENDED')
        call card_boolean(77,write_betasusp,ierr) !Betasusp
        
      case('BEDLOAD_VELOCITY_OUTPUT','OUTPUT_BEDLOAD_VELOCITY',&
          'BED_LOAD_VELOCITY_OUTPUT','OUTPUT_BED_LOAD_VELOCITY')
        call card_boolean(77,write_bedvel,ierr) !ub
        
      !Waves  
      case('WAVE_DISS_OUTPUT','WAVE_DISSIP_OUTPUT',&
           'WAVE_BREAK_DISS_OUTPUT','WAVE_DISSIPATION_OUTPUT')
        call card_boolean(77,write_wavbrkdiss,ierr)
      
      case('WAVE_RAD_OUTPUT','WAVE_STRESS_OUTPUT','WAVE_RADIATION_STRESS_OUTPUT')
        call card_boolean(77,write_wavstress,ierr)
        
      case('WAVE_BREAK_INDEX_OUTPUT','WAVE_BREAKING_INDEX_OUTPUT')
        call card_boolean(77,write_wavbrkind,ierr)
      
      case('WAVE_FRAC_BREAK_OUTPUT')
        call card_boolean(77,write_fracbreak,ierr)
      
      case('WAVE_DISSIP_COEF_OUTPUT','WAVE_DISSIP_COEFF_OUTPUT')
        call card_boolean(77,write_wavdisscoef,ierr)
      
      case('WAVE_BOT_ORB_VEL_OUTPUT','WAVE_ORB_VEL_OUTPUT')
        call card_boolean(77,write_wavorbvel,ierr)
      
      case('WAVE_LENGTH_OUTPUT','OUTPUT_WAVE_LENGTH')
        call card_boolean(77,write_wavlen,ierr)
      
      case('WAVE_VELOCITY_OUTPUT','OUTPUT_WAVE_VELOCITY',&
        'WAVE_VEL_OUTPUT','WAVE_FLUX_VELOCITY_OUTPUT')
        call card_boolean(77,write_wavvel,ierr)
        
      case('WAVE_VELOCITY_MAGNITUDE_OUTPUT','OUTPUT_WAVE_VELOCITY_MAGNITUDE',&
        'WAVE_VEL_MAG_OUTPUT','WAVE_FLUX_VELOCITY_MAG_OUTPUT','WAVE_FLUX_VEL_MAG_OUTPUT')
        call card_boolean(77,write_wavvelmag,ierr)
        
      case('ROLLER_ENERGY_OUTPUT','OUTPUT_ROLLER_ENERGY')
        call card_boolean(77,write_rollenergy,ierr)
      
      case('ROLLER_STRESS_OUTPUT','OUTPUT_ROLLER_STRESS')
        call card_boolean(77,write_rollstress,ierr)
        
      case('ROLLER_DISSIPATION_OUTPUT','OUTPUT_ROLLER_DISSIPATION')
        call card_boolean(77,write_rolldiss,ierr)  
      
      case('OUTPUT_WIND_VELOCITY','WIND_VELOCITY_OUTPUT')
        call card_boolean(77,write_wndvel,ierr) 
        
      case('WIND_MAG_OUTPUT','OUTPUT_WIND_MAG','WIND_MAGNITUDE_OUTPUT','OUTPUT_WIND_MAGNITUDE')
        call card_boolean(77,write_wndmag,ierr) 
        
      case('WIND_STRESS_OUTPUT','OUTPUT_WIND_STRESS')
        call card_boolean(77,write_wndstress,ierr)   
      
      case('WIND_STRESS_MAG_OUTPUT','OUTPUT_WIND_STRESS_MAG')
        call card_boolean(77,write_wndstressmag,ierr)  
      
      case('ATM_PRESS_OUTPUT','OUTPUT_ATM_PRESS','ATM_PRES_OUTPUT','OUTPUT_ATM_PRES',&
          'OUTPUT_ATMOSPHERIC_PRESSURE')
        call card_boolean(77,write_atmpres,ierr)
        
      case('ATM_PRESS_GRAD_OUTPUT','OUTPUT_ATM_PRESS_GRAD','ATM_PRES_GRAD_OUTPUT','OUTPUT_ATM_PRES_GRAD',&
          'OUTPUT_ATMOSPHERIC_PRESSURE_GRADIENT')
        call card_boolean(77,write_atmpresgrad,ierr)  
        
      case('MEAN_BED_SHEAR_STRESS_OUTPUT','OUTPUT_MEAN_BED_SHEAR_STRESS')
        call card_boolean(77,write_meanbedshear,ierr)
        
      case('NORMALIZED_BED_ROUGHNESS_OUTPUT','NORMALIZED_ROUGHNESS_OUTPUT',&
           'NORM_BED_ROUGH_OUTPUT')
        call card_boolean(77,write_normrough,ierr)
        
      case('NORMALIZED_APPARENT_ROUGHNESS_OUTPUT','NORM_BED_APP_ROUGH_OUTPUT')
        call card_boolean(77,write_normapprough,ierr)  
        
      case('ROUGHNESS_HEIGHT_OUTPUT','OUTPUT_ROUGHNESS_HEIGHT')
        call card_boolean(77,write_rough,ierr) 
      
      case('RIPPLE_ROUGHNESS_HEIGHT_OUTPUT','OUTPUT_RIPPLE_ROUGHNESS_HEIGHT')
        call card_boolean(77,write_roughrip,ierr) 
      
      case('MEGARIPPLE_ROUGHNESS_HEIGHT_OUTPUT','OUTPUT_MEGARIPPLE_ROUGHNESS_HEIGHT')
        call card_boolean(77,write_roughmegarip,ierr) 
      
      case('DUNE_ROUGHNESS_HEIGHT_OUTPUT','OUTPUT_DUNE_ROUGHNESS_HEIGHT')
        call card_boolean(77,write_roughdune,ierr) 
      
      case('GRAIN_ROUGHNESS_HEIGHT_OUTPUT','OUTPUT_GRAIN_ROUGHNESS_HEIGHT')
        call card_boolean(77,write_roughgrain,ierr)
        
      case('TRANSPORT_ROUGHNESS_HEIGHT_OUTPUT','OUTPUT_TRANSPORT_ROUGHNESS_HEIGHT')
        call card_boolean(77,write_roughtrans,ierr)  
        
!  Save point cards - added 5/7/12 meb
      case('HYDRO_OUTPUT_INTERVAL','FLOW_OUTPUT_INTERVAL')
        call card_scalar(77,'hrs','sec',savept(1)%time_inc,ierr)

      case('SEDIMENT_OUTPUT_INTERVAL')
        call card_scalar(77,'hrs','sec',savept(2)%time_inc,ierr)
                
      case('SALINITY_OUTPUT_INTERVAL','SALT_OUTPUT_INTERVAL')  
        call card_scalar(77,'hrs','sec',savept(3)%time_inc,ierr)        
        
      case('WAVE_OUTPUT_INTERVAL','WAVES_OUTPUT_INTERVAL')
        call card_scalar(77,'hrs','sec',savept(4)%time_inc,ierr)        
      
      case('MORPHOLOGY_OUTPUT_INTERVAL','MORPH_OUTPUT_INTERVAL')
        call card_scalar(77,'hrs','sec',savept(5)%time_inc,ierr)     
        
      case('SAVE_POINT_LABEL')
        backspace(77)
        read(77,*) cardname, splabel

      case('SAVE_POINT_BEGIN')
        call savept_block
        
      case('SAVE_POINT')
        backspace(77)
        read(77,'(A120)') saveargs
        read(saveargs,*,iostat=ierr) cardname,ptname,cell,xsave,ysave,(group(i),i=1,ngroups)
        call savept_add(ptname,cell,xsave,ysave,ngroups,group)
    
      !case('SAVE_POINT_FILE') !This would be useful for having many save points
      !  backspace(77)
      !  read(77,'(A120)') saveargs
      !  read(saveargs,*,iostat=ierr) cardname,file,(group(i),i=1,ngroups)
        
      case default 
        foundcard = .false.   
        
    endselect
    
    return
    contains        
    
    !-----------------------------------------------------------------------------
        subroutine savept_add(ptname,cell,xsave,ysave,ngroups,group)
        !Adds (activates) a save point
        use prec_def
        use out_def, only: save_point
        implicit none
        !Input/Output
        integer :: cell,ngroups
        character(len=*) :: ptname,group(ngroups)
        real(ikind) :: xsave,ysave
        !Internal
        integer :: i,igroup
        
        do i=1,ngroups
          if(len_trim(group(i))>0) then
            igroup = savept_group(group(i))
            if(igroup/=0)then
              call savept_add_group(igroup, ptname, cell, xsave, ysave)
              group(i)=''
              if(.not.save_point) save_point = .true.
            endif
          endif
        enddo
        
        return
        endsubroutine savept_add
        
    !-------------------------------------------------------------
        subroutine savept_block
        use prec_def
        implicit none
        integer       :: cell,kk,ierr,i
        character(len=34) :: ptname
        character(len=100) :: groupstr
        real(ikind)   :: xsave,ysave
        logical :: foundcard
        
        do i=1,ngroups
          group(i) = ''
        enddo
        cell = 0
        
        do kk=1,10
          foundcard = .true.  
          read(77,*,iostat=ierr) cardname
          if(ierr/=0) exit
          if(cardname(1:1)=='!' .or. cardname(1:1)=='#') cycle      
          selectcase(cardname)
          case('NAME','LABEL')
            backspace(77)
            read(77,*) cardname,ptname 
            
          case('CELL')
            backspace(77)
            read(77,*) cardname,cell 
            
          case('COORDINATES')
            backspace(77)
            read(77,*) cardname,xsave,ysave 
            
          case('GROUPS')
            backspace(77)
            read(77,*) cardname,groupstr
            read(groupstr,*,iostat=ierr) (group(i),i=1,ngroups)
            
          case('SAVE_POINT_END','SAVEPOINT_END','END')
            exit
          
          case default
            foundcard = .false.
            write(*,*) 'WARNING: Card ',cardname,' not found in save point block'
            write(*,*) '  Check input CMS Card File. Skipping card.'
          
          endselect
        enddo
        
        call savept_add(ptname,cell,xsave,ysave,ngroups,group)
        
        return
        endsubroutine savept_block
        
    !-------------------------------------------------------------
        subroutine savept_add_group(igroup,name,cell,x,y)
        use out_def, only: save_pt_group, maxout
        use diag_lib, only: diag_print_error
        use diag_def, only: msg, msg2
        use prec_def
        implicit none
        integer       :: igroup,num_cells,cell
        character(len=*) :: name
        real(ikind)   :: x,y
        type(save_pt_group) :: tempsave
        
        num_cells = savept(igroup)%ncells
        if (num_cells==0) goto 100
        if (allocated(tempsave%id))   deallocate (tempsave%id)
        if (allocated(tempsave%cell)) deallocate (tempsave%cell)
        if (allocated(tempsave%x))    deallocate (tempsave%x)
        if (allocated(tempsave%y))    deallocate (tempsave%y)
        allocate (tempsave%id(num_cells),tempsave%x(num_cells),tempsave%y(num_cells),tempsave%cell(num_cells))
        do i=1,num_cells
          tempsave%id(i)   = savept(igroup)%id(i)
          tempsave%cell(i) = savept(igroup)%cell(i)
          tempsave%x(i)    = savept(igroup)%x(i)
          tempsave%y(i)    = savept(igroup)%y(i)
        enddo
        deallocate(savept(igroup)%id, savept(igroup)%x, savept(igroup)%y, savept(igroup)%cell)
100     allocate(savept(igroup)%id(num_cells+1), savept(igroup)%cell(num_cells+1))
        allocate(savept(igroup)%x(num_cells+1),  savept(igroup)%y(num_cells+1))
        do i=1,num_cells
          savept(igroup)%id(i)   = tempsave%id(i)
          savept(igroup)%cell(i) = tempsave%cell(i)
          savept(igroup)%x(i)    = tempsave%x(i)
          savept(igroup)%y(i)    = tempsave%y(i)
        enddo
        num_cells = num_cells+1
        savept(igroup)%active       = .true.
        savept(igroup)%id(num_cells)   = '"'//trim(name)//'"'
        savept(igroup)%cell(num_cells) = cell
        savept(igroup)%x(num_cells)    = x
        savept(igroup)%y(num_cells)    = y
        savept(igroup)%ncells       = num_cells
        
        if(num_cells > maxout) then    !MAXOUT defined in out_def
          write(msg,'(I0)') maxout
          msg2 = 'Number of Save Points exceeds the present limit of '//trim(msg)
          call diag_print_error(msg2)
        endif
        
        return
        endsubroutine savept_add_group
        
    !------------------------------------------------------------  
        function savept_group(group) result(idum)
        implicit none        
        character(len=*) :: group
        integer ::      idum
        
        selectcase(group)
        case('HYDRO','FLOW')
          idum = 1
        case('SEDIMENT','SED')
          idum = 2
        case('SALINITY','SAL','SALT')
          idum = 3
        case('WAVE','WAVES')
          idum = 4
        case('MORPH','MORPHOLOGY')
          idum = 5
        case default
          idum = 0
        end select
        
        return
        endfunction savept_group

    endsubroutine out_cards

!*************************************************************     
    subroutine out_init()
! Initializes the output variables
!
! written by Alex Sanchez, USACE-ERDC-CHL
!*************************************************************
#include "CMS_cpp.h"
    use flow_def, only: waveflux
    use comvarbl, only: iyr,imo,iday,ihr,imin,isec,timehrs,reftime,pred_corr,casename
    use geo_def, only: areap,bathydata
    use sed_def
    use sal_def
    use heat_def
    use cms_def
    use met_def
    use hot_def, only: coldstart,hotfile
    use bnd_def
    use rol_def, only: roller
    use out_def
    use out_lib
    use diag_lib
    use comvarbl, only: flowpath,nfsch,input_ver
    use diag_def, only: debug_mode
    use wave_flowgrid_def, only: constant_waves
#ifdef XMDF_IO
    use xmdf
    use IFPORT  !XMDF infers on Windows
#endif   

    implicit none      
    integer :: i,ierr,ncase,npath,nn
    character(len=200) :: apath,aname,outdirpath
    character(len=10) :: aext
    logical :: foundfile,found,created
    
#ifdef XMDF_IO
    call XF_CALENDAR_TO_JULIAN(0,iyr,imo,iday,ihr,imin,isec,reftime,ierr)       !Set reference time for XMDF from month, day, year, hour, minute and second.
    !call XF_JULIAN_TO_CALENDAR(0,iyr,imo,iday,ihr,imin,isec,reftime,error)
#endif

    npath = len_trim(outpath)
    ncase = len_trim(casename)  
    outprefix = outpath(1:npath) // casename(1:ncase)
    nn = npath+ncase
    hisfile           = outprefix(1:nn) // '_tecplot.his' !Tecplot History (Animation) File
    datfile           = outprefix(1:nn) // '_tecplot.dat' !Tecplot Snapshot file)
    goutfile          = outprefix(1:nn) // '_sol.h5'      !XMDF Global Solution File
    dgoutfile         = outprefix(1:nn) // '_diag.h5'     !XMDF Global Diagnostic File
    
    !If other filenames have been assigned with the use of a card, "USE_DEFAULT" for that type is set to FALSE and reassignment is skipped.   !MEB  02/16/2016
    do i=1,14
      if (outlist(i)%use_default .eqv. .false.) cycle     !Changed == to .eqv. for compatibility
      select case (i)
      case (1)
          outlist(i)%afile  = outprefix(1:nn) // '_wse.h5'       !Water surface elevation
      case (2)
          outlist(i)%afile  = outprefix(1:nn) // '_vel.h5'       !Current velocity
      case (3)
          outlist(i)%afile  = outprefix(1:nn) // '_visc.h5'      !Eddy viscosity    
      case (4)
          outlist(i)%afile  = outprefix(1:nn) // '_sed.h5'       !Sediment
      case (5)
          outlist(i)%afile  = outprefix(1:nn) // '_morph.h5'     !Morphology - Depth, Morphology Change
      case (6)
          outlist(i)%afile  = outprefix(1:nn) // '_salt.h5'      !Salinity - Concentrations
      case (7)
          outlist(i)%afile  = outprefix(1:nn) // '_wave.h5'      !Wave 
      case (8)
          outlist(i)%afile  = outprefix(1:nn) // '_met.h5'       !Meteorological 
      case (9)
          outlist(i)%afile  = outprefix(1:nn) // '_trans.h5'     !Transport 
      case (10)
          outlist(i)%afile  = outprefix(1:nn) // '_bedcomp.h5'   !Bed Composition
      case (11)
          outlist(i)%afile  = outprefix(1:nn) // '_sedfrac.h5'   !Sed Mix 
      case (12)
          outlist(i)%afile  = outprefix(1:nn) // '_wave.h5'      !Wave details
      case (13)
          outlist(i)%afile  = outprefix(1:nn) // '_fric.h5'      !Bed Friction/Roughness
      case (14)
          outlist(i)%afile  = outprefix(1:nn) // '_temp.h5'      !Temperature
      case default
          outlist(i)%afile  = outprefix(1:nn) // '_undefined.h5' !Undefined new file
      end select
    enddo
    
    !Override default names for SMS 11.0 and CMS version 3.75 - backward compatibility
    if(input_ver<3.8 .or. use_common)then
      do i=1,14
        outlist(i)%afile = goutfile
      enddo
    endif

    !--- Initialize output list, only for reference purposes ---
    outlist(1)%name  = 'WSE'         !Water surface elevation
    outlist(2)%name  = 'VEL'         !Current mangitude and direction
    outlist(3)%name  = 'EDDY'        !Eddy viscosity
    outlist(4)%name  = 'SED'         !Concentration concentration, capacity and total load transport vector
    outlist(5)%name  = 'MORPH'       !Water depth, and morphology change
    outlist(6)%name  = 'SALT'        !Salinity concentration
    outlist(7)%name  = 'WAVE'        !Wave height, period, and vector
    outlist(8)%name  = 'MET'         !Wind stress vector and magnitude
    outlist(9)%name  = 'TRANS'       !Sediment concentration, capacity, transport vector and salinity concentration
    outlist(10)%name = 'BEDCOMP'     !Bed layer thickness and fraction for nonuniform sediment
    outlist(11)%name = 'SEDFRAC'     !Sediment transport details for nonuniform sediment
    outlist(12)%name = 'WAVEADV'     !Advanced wave information
    outlist(13)%name = 'FRIC'        !Bed friction/roughness
    outlist(14)%name = 'TEMPERATURE' !Temperature
    
    do i=1,noutlist    
      outlist(i)%write_dat = .false.
      outlist(i)%inc = 1 
    enddo        

    !Hydrodynamics
    if(.not.pred_corr)then
      write_hpred = .false.
      write_wsepred = .false.
      write_velpred = .false.
    endif
    
    !--- Force output to be zero if no information is calculated ---    
    !Sediment
    if(.not.sedtrans)then
      outlist(4)%ilist = 0 !Sediment
      if(.not.bathydata%ison)then
        outlist(5)%ilist = 0 !Morphology
      endif
      outlist(10)%ilist = 0 !Bed thickness and fractions
      outlist(11)%ilist = 0 !Advanced sediment
    else
      if(singlesize)then
        outlist(10)%ilist = 0 !Bed thickness and fractions
        outlist(11)%ilist = 0 !Advanced sediment
      else !Multiple-grain sizes
        if(outlist(10)%ilist==0)then !If not given
          outlist(10)%ilist = outlist(5)%ilist !use morphology output list
          !outlist(10) = outlist(5) !use morphology output list
          outlist(10)%write_dat = .true.
          !if(input_ver<4.0) outlist(10)%afile = goutfile
        endif
      endif   
    endif
    
    !Salinity
    if(.not.saltrans)then
      outlist(6)%ilist = 0 !Salinity
    endif
    
    !Temperature
    if(.not.heattrans)then
      outlist(14)%ilist = 0 !Temperature
    else
      outlist(14)%ilist = 1                   !Otherwise, Temperature never got written out.  
      outlist(14)%write_dat = .true.
    endif
    
    !Transport
    if(.not.sedtrans .and. .not.saltrans .and. .not.heattrans)then  !Alex, added salinity
      outlist(9)%ilist = 0 !Trans      
    endif
    if(.not.sedtrans)then
      write_morph = .false.  !Morphology change 
      if(.not.bathydata%ison)then
        outlist(5)%ilist = 0 !Morphology group
      endif
    endif
    if(iadapttot/=3)then
      write_lengthtot = .false.
    endif
    if(iadapttot<5)then
      write_lengthsusp = .false.
      write_lengthbed = .false.
    endif
    
    !Waves
    if(noptset<=2 .and. .not.constant_waves) then    !Modified by Mitch to print out wave information if CONSTANT_WAVES_xxx are specified.  06/23/2017
      outlist(7)%ilist  = 0
      outlist(12)%ilist = 0        !Output details outdated
      write_wave_detail = .false.  !Wave details
      write_wavbrkdiss  = .false.  !Wave dissipation 
      write_wavstress   = .false.  !Wave stress vectors
      write_fracbreak   = .false.  !Fraction of breaking waves
      write_wavdisscoef = .false.  !Wave dissipation coefficient
      write_wavorbvel   = .false.  !Wave bottom orbital velocity
      write_wavlen      = .false.  !Wave length
      write_wavvel      = .false.  !Wave mass transport velocity
      write_wavvelmag   = .false.  !Wave mass transport velocity magnitude
      write_rollenergy  = .false.  !Roller energy
      write_rollstress  = .false.  !Roller stress   
    else
      if(write_wave_detail)then
        write_wavbrkdiss  = .true.
        write_wavbrkind  = .true.
        write_wavstress  = .true.
        write_fracbreak  = .true.
        write_wavorbvel  = .true.
        write_wavlen     = .true.
        write_rollenergy = .true.
        write_rollstress = .true.
      endif
      if(.not.roller)then
        write_rollenergy = .false.
        write_rollstress = .false.
        write_rolldiss   = .false.
      endif
      if(.not.waveflux)then
        write_wavvel = .false.
        write_wavvelmag = .false.
      endif
    endif        
    
    !Wind
    if(.not.windsta .and. .not.windconst .and. .not.windvar)then
      outlist(8)%ilist = 0                      !Alex, Bug Fix, Sept 1, 2009, changed 7 to 8
      write_wndmag = .false.
      write_wndstress = .false.
      write_wndstressmag = .false.
    endif
    
    !Atmospheric pressure
    if(.not.presconst .and. .not.presvar)then
      write_atmpres = .false.
      write_atmpresgrad = .false.
    endif
    
    !Implicit vs explicit schemes
    if(nfsch==0 .and. debug_mode)then
      write_presres = .true.
      write_velres = .true.
    elseif(nfsch>0)then
      write_presres = .false.
      write_velres = .false.
    endif
    
    !Friction/Roughness group
    if((write_meanbedshear .or. write_normrough .or. write_normapprough) &
      .and. outlist(13)%ilist==0)then
        outlist(13) = outlist(2)
    endif
    if(outlist(13)%ilist>0)then
      write_rough = .true.        !Turns on or off the normalized bed roughness height
      write_roughrip = .true.     !Turns on or off the roughness due to rippples
      write_roughmegarip = .true. !Turns on or off the roughness due to megarippples
      write_roughdune = .true.    !Turns on or off the roughness due to dunes
      write_roughgrain = .true.   !Turns on or off the roughness due to grains
      write_roughtrans = .true.   !Turns on or off the roughness due due to sediment transport
      write_mannings = .true.     !Turns on or off the Manning's roughness coefficient
    endif
    
    !--- Initialize each dataset time lists and variables ----
    npath = len_trim(flowpath) 
    do i=1,noutlist
      !sets up times for a generic dataset    
      if(outlist(i)%ilist==0)then
        outlist(i)%write_dat = .false.
        outlist(i)%ntimes = 0
      else
        outlist(i)%ntimes = outseries(outlist(i)%ilist)%ntimes
        if(outseries(outlist(i)%ilist)%ntimes==0)then
          outlist(i)%write_dat = .false.
        else
          if(allocated(outlist(i)%times)) deallocate(outlist(i)%times)
          allocate(outlist(i)%times(outlist(i)%ntimes))
          outlist(i)%times = outseries(outlist(i)%ilist)%times
          outlist(i)%inc = 1 
          outlist(i)%write_dat = .true.
        endif
      endif  

      if(outlist(i)%write_dat)then  !get increment and sim time synchronized
        do while(timehrs > (outlist(i)%times(outlist(i)%inc)+1.0e-4))
          outlist(i)%inc = outlist(i)%inc + 1
        enddo
        call fileparts(outlist(i)%afile,apath,aname,aext)
        if(len_trim(apath)==0)then
          outlist(i)%afile = flowpath(1:npath) // outlist(i)%afile
        endif        
      endif
    enddo    
    
    !Separate trans group into sediment and salinity groups
    if(outlist(9)%ilist>=1)then
      if(outlist(4)%ilist==0 .and. sedtrans)then
        outlist(4) = outlist(9) !Overide sediment with transport information
      endif
      if(outlist(6)%ilist==0 .and. saltrans)then
      !if(outlist(6)%ilist==0 .and. saltrans .and. input_ver<4.0)then    
        outlist(6) = outlist(9) !Overide salinity with transport information
      endif 
      if(outlist(14)%ilist==0 .and. heattrans)then
        outlist(14) = outlist(9) !Overide temperature with transport information
      endif
    endif  
    
    !Advanced wave information
    if(outlist(7)%ilist>=1)then
      if(outlist(12)%ilist==0 .and. write_wave_detail .and. cmswave)then
        outlist(12) = outlist(7) !Overide advanced wave with wave information
      endif  
    endif
        
    if(.not.waveflux)then
      write_veltotal = .false.
      write_veltotalmag = .false.
      write_wavvel = .false.
      write_wavvelmag = .false.
    endif        
    if(.not.write_veltotal) write_veltotalmag = .false.
    
    !Delete files if necessary 
    if(coldstart)then
      inquire(file=goutfile,exist=foundfile)
      if(foundfile)then
        open(100,file=goutfile)
        close(100,status='delete',iostat=ierr)    
        if (ierr < 0) then 
          call diag_print_error('Could not delete output file- '//trim(goutfile))
        endif          
      endif
      do i=1,14
        if(outlist(i)%write_dat)then
          inquire(file=outlist(i)%afile,exist=foundfile)
          if(foundfile)then
            open(100,file=outlist(i)%afile)
            close(100,status='delete',iostat=ierr)
            if (ierr < 0) then
              call diag_print_error('Could not delete output file- '//trim(outlist(i)%afile))
            endif
          endif
        endif  
      enddo
!     open(100,file=hotfile)
!     close(100,status='delete',err=943)      
    endif
    
    inquire(file=dgoutfile,exist=foundfile) 
    if(foundfile)then
      open(100,file=dgoutfile)
      close(100,status='delete',iostat=ierr)     
      if (ierr < 0) then
        call diag_print_error('Could not delete output file- '//trim(dgoutfile))
      endif
    endif
    
    !Initialize Super ASCII files
    if(write_sup)then
      !Save all these solution files to a subdirectory named "ASCII_Solutions"
      outdirpath='ASCII_Solutions'
#ifdef _WIN32      
      inquire(directory=trim(outdirpath), exist=found)
      if(.not.found) then
        created=MakeDirQQ(trim(outdirpath))
        if(.not.created)then
          call diag_print_error('Failed to create subdirectory- '//trim(outdirpath))
        endif
      endif
#else
      inquire(file=trim(outdirpath), exist=found)
      if(.not.found) then
        call system('mkdir '//(trim(outdirpath)))
      endif
#endif
      aname = trim(flowpath) // trim(outdirpath) // '/' //  trim(casename)
      call write_sup_file(aname) !Super File *.sup
      call write_xy_file(aname,casename)  !XY coordinate file *.xy  
    endif
    
    !Write static variables
    if(write_areap)then
      if(write_sup)then  
        aname = trim(outdirpath) // '/' // trim(flowpath) // trim(casename)  
        call write_scal_dat_file(aname,'area','area',areap)
      endif
#ifdef XMDF_IO
      apath = trim(simlabel) // '/'
      call writescalh5(goutfile,apath,'Cell_Area',areap,'m^2',timehrs,1)  !Cell Area
#endif     
    endif

    return
!-----------------------------------------------------------
!978 call diag_print_error('Could not access output file ')
!    return

    endsubroutine out_init
    
!*******************************************************************
    subroutine out_print()
! Displays the Output settings to the screen and the diagnostic file
! written by Alex Sanchez, USACE-CHL    
!*******************************************************************    
    use out_def
    use geo_def, only: mapid
    use sed_def, only: sedtrans
    use sal_def, only: saltrans
    use heat_def, only: heattrans
    use diag_def, only: dgunit,dgfile,debug_mode
    
    implicit none
    integer :: i,ii,j,k,iunit(2),ierr
    character(len=200) :: apath,aname,astring
    character(len=10) :: aext

787 format(' ',A,T40,A)    
485 format(' ',A,T40,F0.2,A)
262 format(' ',A,T40,I0,A)

    iunit = (/6,dgunit/)
    open(dgunit,file=dgfile,access='append')
    do i=1,2
      write(iunit(i),*)
      write(iunit(i),787)       'Global Output'
      write(iunit(i),787)       '  Simulation Label:',trim(simlabel)
      if(.not.write_xmdf_output) then
        write(iunit(i),787)     '  Output Path:','ASCII_Solutions/'
      else
        if(len_trim(outpath)>0)then
          write(iunit(i),787)   '  Output Path:',trim(outpath)
        endif
        call fileparts(outlist(1)%afile,apath,aname,aext)
        astring=trim(aname) // '.' // aext              
        write(iunit(i),787)     '  Water Level File:',trim(astring)
        call fileparts(outlist(2)%afile,apath,aname,aext)
        astring=trim(aname) // '.' // aext      
        write(iunit(i),787)     '  Current Velocity File:',trim(astring)
        if(outlist(3)%write_dat)then
          call fileparts(outlist(3)%afile,apath,aname,aext)
          astring=trim(aname) // '.' // aext     
          write(iunit(i),787)   '  Eddy Viscosity File:',trim(astring)
        endif
        if(sedtrans .or. saltrans)then
          call fileparts(outlist(9)%afile,apath,aname,aext)
          astring=trim(aname) // '.' // aext
          write(iunit(i),787)   '  Transport File:',trim(astring)
        endif   
        if(heattrans)then
          call fileparts(outlist(14)%afile,apath,aname,aext)
          astring=trim(aname) // '.' // aext
          write(iunit(i),787)   '  Temperature File:',trim(astring)
        endif   
        if(sedtrans)then
          call fileparts(outlist(5)%afile,apath,aname,aext)
          astring=trim(aname) // '.' // aext     
          write(iunit(i),787)   '  Morphology Change File:',trim(astring)
        endif
        if(outlist(7)%write_dat)then
          call fileparts(outlist(7)%afile,apath,aname,aext)
          astring=trim(aname) // '.' // aext 
          write(iunit(i),787)   '  Wave File:',trim(astring)
        endif
        if(outlist(8)%write_dat)then
          call fileparts(outlist(8)%afile,apath,aname,aext)
            astring=trim(aname) // '.' // aext 
          write(iunit(i),787)   '  Meteorological File:',trim(astring)
        endif
        if(outlist(10)%write_dat)then
          call fileparts(outlist(10)%afile,apath,aname,aext)
            astring=trim(aname) // '.' // aext
          write(iunit(i),787)   '  Bed Composition File:',trim(astring)
        endif  
        if(outlist(13)%write_dat)then
          call fileparts(outlist(13)%afile,apath,aname,aext)
          astring=trim(aname) // '.' // aext
          write(iunit(i),787)   '  Bed Friction/Roughness File:',trim(astring)
        endif
        if(write_tecplot)then
          write(iunit(i),787)   '  Tecplot Snapshot File:',trim(datfile)
        endif
        if(debug_mode)then
          write(iunit(i),787)   '  Diagnostic File:',trim(dgoutfile)
        endif      
      endif
      if(obs_cell)then
        !write(iunit(i),*)  
        write(iunit(i),787)     'Observation Cells:','ON' 
        if(obs(1)%active)then    
          write(iunit(i),485)   '  Time series increment:',obs(1)%time_inc,' sec'
        endif
        if(obs(2)%active)then    
          write(iunit(i),485)   '  Flow rate increment:',obs(2)%time_inc,' sec'
        endif
        if(obs(3)%active)then    
          write(iunit(i),485)   '  Transport increment:',obs(3)%time_inc,' sec'
        endif
        if(obs(4)%active)then
          write(iunit(i),485)   '  Bed increment:',obs(4)%time_inc,' sec'
        endif
      endif
      if(write_maxwse) then
        write(iunit(i),787)     '  Maximum Water Level Output:','ON'
      endif

234   format(' ',4x,A12,A40,'''',A,'''')
543   format(2x,I6,5x,A20,I6,6x,F12.3,8x,F12.3)
      
      if(save_point)then
        write(iunit(i),*)
        write(iunit(i),787)     'Save Points:','ON'   
        do k=1,ngroups
          if(.not.savept(k)%active) cycle
          write(iunit(i),787)   '  Group:',trim(savept(k)%group)
          write(iunit(i),485)   '   Output Interval:',savept(k)%time_inc,' sec'
          write(iunit(i),262)   '   Number of variables:',savept(k)%nvar
          write(iunit(i),787)   '    Variable,   File,                                   Units'
          do j=1,savept(k)%nvar
            call fileparts(savept(k)%files(j),apath,aname,aext)
              astring=trim(aname) // '.' // aext
            write(iunit(i),234,iostat=ierr) adjustl(savept(k)%names(j)),adjustl(astring),adjustl(trim(savept(k)%ounit(j)))
          enddo
          write(iunit(i),262)   '   Number of Points:',savept(k)%ncells
          write(iunit(i),787)   '  Point ID, Point Name,         Cell ID,   X-Coordinate (m),  Y-Coordinate (m)'
          do j=1,savept(k)%ncells
            ii = savept(k)%cell(j)
            if(allocated(mapid))then  
              ii = mapid(ii)
            endif
            write(iunit(i),543,iostat=ierr) j,adjustl(savept(k)%id(j)),ii,savept(k)%x(j),savept(k)%y(j)
          enddo
        enddo
      endif
    enddo    
    close(dgunit)
    
    return
    endsubroutine out_print

!**********************************************************************    
    subroutine write_output()
! Writes CMS-Flow Output
! Author: Alex Sanchez, USACE-ERDC-CHL    
!***********************************************************************
#include "CMS_cpp.h"
    use out_def, only: write_sup,obs_cell,save_point,&
       incXMDF_IO,write_tecplot,write_xmdf_output
    implicit none
    logical :: header
    
    header = .false.
    
    if(write_tecplot) call write_output_tecplot  !Tecplot files    
    
#ifdef XMDF_IO
    if(write_xmdf_output) call write_output_xmdf(header)
#endif   
    if(write_sup) call write_output_sup(header)    
    
    if(obs_cell)   call write_obs_cell   !Observation cells
    if(save_point) call write_save_point !Save point cells, Mitch 5/8/2012
    
    return
    endsubroutine write_output
    
!***********************************************************************    
    subroutine write_output_tecplot
! Writes the Tecplot output files
! Author: Weiming Wu, NCCHE
!***********************************************************************    
    use comvarbl, only: ntime
    use out_def, only: hisfile
    use out_lib, only: write_tecplot_dat,write_tecplot_his
    implicit none
    
    if(mod(ntime,2)==0)  call write_tecplot_dat  !Snapshot file
    if(ntime==0) open(30,file=hisfile)
    if(mod(ntime,12)==0) call write_tecplot_his  !History/animation file    
      
    return
    endsubroutine write_output_tecplot
    
!***********************************************************************    
    subroutine write_output_xmdf(header)
! Writes the XMDF output to each variable group    
!    
! Author: Alex Sanchez, USACE-CHL    
!***********************************************************************
#include "CMS_cpp.h"    
#ifdef XMDF_IO
    use der_def, only: nlim,rgx,rgy,gow
    use der_lib, only: curvxy
    use comvarbl, only: flowpath,casename,timehrs,nfsch,norder
    use diag_def, only: debug_mode
    use flow_def
    use flow_lib, only: streamwise_curvature
    use fric_def
    use fric_lib, only: fric_normapprough
    use geo_def, only: zb,zb0,areap
    use met_def
    use out_def
    use out_lib, only: writescalh5,writevech5
    use prec_def
    use rol_def, only: roller,rxrs,ryrs,roldiss
    use sal_def, only: sal
    use heat_def, only: heat
    use sed_def
    use size_def
    use wave_flowgrid_def, only: wunitx,wunity,whgt,wper,Wang,           &
        !wavediss,waveibr,Worb,wlen,wavestrx,wavestry,Ssr,wavstrx,wavstry
        wavediss,Worb,wlen,wavestrx,wavestry,Ssr,wavstrx,wavstry
    use dredge_def
#ifdef DEV_MODE
    use q3d_def, only: q3d,f3dxx,f3dxy,f3dyy,f3du,f3dv,wavcurint,        &
      udsx,udsy,uzplay,vzplay,nzplay,q3d_lay
#endif

implicit none    
    integer :: i,ii,j,ks,nn
    real(ikind) :: val(ncellsD),vecx(ncellsD),vecy(ncellsD),valdry
    character(len=5) :: apbk,alay
    character(len=3) :: aperdiam
    character(len=200) :: apath,aname
    logical :: header,check_time_list
    
62  format('_',I2.2)
71  format(1x,'(',I1,')')
72  format(1x,'(',I2,')')
    
    valdry = -999.0
    
    nn = len_trim(simlabel)
    apath = simlabel(1:nn)//'/'
    
    aname = trim(flowpath) // trim(casename)
    
    !WATER SURFACE ELEVATION GROUP
    if(check_time_list(1))then
      call print_output_header(header)
      call writescalh5(outlist(1)%afile,apath,'Water_Elevation',eta,'m',timehrs,0)
      if(write_presres)then !Pressure Residuals
        call writescalh5(outlist(1)%afile,apath,'Pres_Norm_Res',rsp,'none',timehrs,0)
        !val = log10(1.0+rsp)
        !call writescalh5(outlist(1)%afile,apath,'Pres_Norm_Res_Scaled',val,'none',timehrs,0)
      endif
      if(write_pres)then
        call writescalh5(outlist(1)%afile,apath,'Pressure',p,'m^2/s^2',timehrs,1)
#ifdef DEV_MODE
        call writescalh5(outlist(1)%afile,apath,'Pressure1',p1,'m^2/s^2',timehrs,1)  
#endif
      endif
      if(write_totdep)then
        call writescalh5(outlist(1)%afile,apath,'Total_Water_Depth',h,'m',timehrs,0)
      endif
      if(write_hpred)then
        call writescalh5(outlist(1)%afile,apath,'Pred_Water_Depth',hpred,'m',timehrs,0)
      endif
      if(write_wsepred)then
        val = hpred + zb !Predicted Water Surface Elevation
        call writescalh5(outlist(1)%afile,apath,'Pred_Water_Elevation',val,'m',timehrs,0)   
      endif
      if(write_wsegrad)then
        if(nfsch==0)then !Implicit
          vecx=dpx*gravinv
          vecy=dpy*gravinv
          call writescalh5(outlist(1)%afile,apath,'detax',vecx,'none',timehrs,0)
          call writescalh5(outlist(1)%afile,apath,'detay',vecy,'none',timehrs,0)
        else !Explicit or Semi-implicit
          call writescalh5(outlist(1)%afile,apath,'detax',detax,'none',timehrs,0)
          call writescalh5(outlist(1)%afile,apath,'detay',detay,'none',timehrs,0)
        endif
      endif
      if(write_presgrad)then
        call writescalh5(outlist(1)%afile,apath,'dpx',dpx,'none',timehrs,1)
        call writescalh5(outlist(1)%afile,apath,'dpy',dpy,'none',timehrs,1)
      endif
      if(write_prescor)then  
        call writescalh5(outlist(1)%afile,apath,'Pressure_Correction',pp,'m^2/s^2',timehrs,1)
      endif  
      if(write_prescorgrad)then        
        call writescalh5(outlist(1)%afile,apath,'dppx',dppx,'m/s^2',timehrs,1)
        call writescalh5(outlist(1)%afile,apath,'dppy',dppy,'m/s^2',timehrs,1)    
      endif
      if(write_volres)then  
        call flow_volres(val)
        val = 100.0*val/(areap*h)
        call writescalh5(outlist(1)%afile,apath,'Vol_Per_Err',val,'',timehrs,0)
      endif
#ifdef DEV_MODE
      if(nfsch==1 .and. norder==2)then !for explicit scheme
        call writescalh5(outlist(1)%afile,apath,'dhx',dhx,'',timehrs,0)
        call writescalh5(outlist(1)%afile,apath,'dhy',dhy,'',timehrs,0)
        if(nlim>0 .and. timehrs<1.0e-9 .and. ncellsimple>0)then
          call writescalh5(outlist(1)%afile,apath,'rgx',rgx,'',timehrs,0)
          call writescalh5(outlist(1)%afile,apath,'rgy',rgy,'',timehrs,0)
        endif
      endif
#endif
      !if(ncellpoly>0 .and. debug_mode)then
      !  call writescalh5(outlist(1)%afile,apath,'Hu',Hu,'m^2/s',timehrs,1)
      !  call writescalh5(outlist(1)%afile,apath,'Hv',Hv,'m^2/s',timehrs,1)    
      !  call writescalh5(outlist(1)%afile,apath,'apuareap',apuareap,'m/s^2',timehrs,1)
      !  call writescalh5(outlist(1)%afile,apath,'sumu',sumu,'m^2/s',timehrs,1)    
      !  call writescalh5(outlist(1)%afile,apath,'spu',spu,'m/s^2',timehrs,1)
      !  call writescalh5(outlist(1)%afile,apath,'spv',spv,'m/s^2',timehrs,1)
      !  call writescalh5(outlist(1)%afile,apath,'sp',sp,'m/s^2',timehrs,1)
      !  if(write_sup)then
      !    call write_scal_dat_file(aname,'Hu','Hu',Hu) !SUPER ASCII File
      !    call write_scal_dat_file(aname,'Hv','Hv',Hv) !SUPER ASCII File
      !    call write_scal_dat_file(aname,'apuareap','apuareap',apuareap) !SUPER ASCII File
      !    call write_scal_dat_file(aname,'sumu','sumu',sumu) !SUPER ASCII File
      !    call write_scal_dat_file(aname,'spu','spu',spu) !SUPER ASCII File
      !    call write_scal_dat_file(aname,'spv','spv',spv) !SUPER ASCII File
      !  endif  
      !endif
!      call writescalh5(outlist(1)%afile,apath,'limX',limx,'',timehrs,1) !Slope limiters
!      call writescalh5(outlist(1)%afile,apath,'limY',limy,'',timehrs,1) !Slope limiters

      if(.not. outlist(5)%write_dat) then
        if(dredging) then
          val=-zb
          call writescalh5(outlist(5)%afile,apath,'Depth',val,'m',timehrs,1)       
          val=zb-zb0
          call writescalh5(outlist(5)%afile,apath,'Morphology_Change',val,'m',timehrs,1)  
        endif
      endif
    endif
    
    !CURRENT VELOCITY GROUP
    if(check_time_list(2))then 
      call print_output_header(header)
      if(write_velmag)then
        call writescalh5(outlist(2)%afile,apath,'Current_Magnitude',uv,'m/s',timehrs,0)  !Eulerian current mangitude  
      endif      
      if(waveflux)then !****************************
        vecx=u-us; vecy=v-vs !Current (Eulerian) velocities    
        call writevech5(outlist(2)%afile,apath,'Current_Velocity',vecx,vecy,'m/s',timehrs,0)
        if(write_veltotal)then
          call writevech5(outlist(2)%afile,apath,'Total_Flux_Velocity',u,v,'m/s',timehrs,0)
          !if(write_sup) call write_vec_dat_file(aname,'Total_Velocity','totvel',u,v) !SUPER ASCII File
          if(write_veltotalmag)then
            do i=1,ncells
              val(i) = sqrt(u(i)*u(i) + v(i)*v(i))
            enddo
            call writescalh5(outlist(2)%afile,apath,'Total_Flux_Velocity_Magnitude',val,'m/s',timehrs,0)
          endif  
        endif          
      else
        call writevech5(outlist(2)%afile,apath,'Current_Velocity',u,v,'m/s',timehrs,0)    
      endif    
      if(write_velpred)then
        call writevech5(outlist(2)%afile,apath,'Predicted_Current_Velocity',upred,vpred,'m/s',timehrs,0)
      endif
      if(write_velgrad)then
        call writescalh5(outlist(2)%afile,apath,'dux',dux,'none',timehrs,0)     
        call writescalh5(outlist(2)%afile,apath,'duy',duy,'none',timehrs,0)
        call writescalh5(outlist(2)%afile,apath,'dvx',dvx,'none',timehrs,0)  
        call writescalh5(outlist(2)%afile,apath,'dvy',dvy,'none',timehrs,0)
!#ifdef DEV_MODE
!        if(nfsch==1 .and. nlim>0)then
!          call writescalh5(outlist(2)%afile,apath,'limdux',limdux,'none',timehrs,0)     
!          call writescalh5(outlist(2)%afile,apath,'limduy',limduy,'none',timehrs,0)
!          call writescalh5(outlist(2)%afile,apath,'limdvx',limdvx,'none',timehrs,0)  
!          call writescalh5(outlist(2)%afile,apath,'limdvy',limdvy,'none',timehrs,0)
!        endif
!#endif
      endif
      if(write_velcurv)then
        call curvxy(gow,dux,duy,vecx,vecy)
        call writescalh5(outlist(2)%afile,apath,'d2udx2',vecx,'1/m/s',timehrs,0)
        call writescalh5(outlist(2)%afile,apath,'d2udy2',vecy,'1/m/s',timehrs,0)
        call curvxy(gow,dvx,dvy,vecx,vecy)
        call writescalh5(outlist(2)%afile,apath,'d2vdx2',vecx,'1/m/s',timehrs,0)
        call writescalh5(outlist(2)%afile,apath,'d2vdy2',vecy,'1/m/s',timehrs,0)
      endif
      if(write_streamcurv)then
        do i=1,ncells
          vecx(i)=u(i)-us(i)
          vecy(i)=v(i)-vs(i)
          val(i) = streamwise_curvature(vecx(i),vecy(i),uv(i),duy(i),duy(i),dvx(i),dvy(i))
        enddo
        call writescalh5(outlist(2)%afile,apath,'Streamwise_Curvature',val,'1/m',timehrs,0)
      endif  
#ifdef DEV_MODE
      if(q3d)then
        call writescalh5(outlist(2)%afile,apath,'f3dxx',f3dxx,'m',timehrs,0)
        call writescalh5(outlist(2)%afile,apath,'f3dxy',f3dxy,'m',timehrs,0)
        call writescalh5(outlist(2)%afile,apath,'f3dyy',f3dyy,'m',timehrs,0)
        call writescalh5(outlist(2)%afile,apath,'f3du',f3du,'m',timehrs,0)
        call writescalh5(outlist(2)%afile,apath,'f3dv',f3dv,'m',timehrs,0)
        if(q3d_lay)then !Velocity layers
          call q3d_vel_lay !uzplay,vzplay
          do j=1,nzplay
            if(j<=9)then
              write(alay,'(I1)') j
            else
              write(alay,'(I2)') j
            endif 
            aname = 'Current_Velocity_Lay' // trim(alay)
            vecx(1:ncellsD) = uzplay(j,1:ncellsD)
            vecy(1:ncellsD) = vzplay(j,1:ncellsD)
            !vecx(1:ncellsD) = 0.0
            !vecy(1:ncellsD) = 0.0
            call writevech5(outlist(2)%afile,apath,aname,vecx,vecy,'m/s',timehrs,0)  
            do i=1,ncells
              val(i) = sqrt(vecx(i)*vecx(i) + vecy(i)*vecy(i))
            enddo
            aname = trim(aname) // '_mag'
            call writescalh5(outlist(2)%afile,apath,aname,val,'m/s',timehrs,0)
          enddo
        endif    
        if(wavcurint)then
          call writescalh5(outlist(2)%afile,apath,'udsx',udsx,'m',timehrs,0)
          call writescalh5(outlist(2)%afile,apath,'udsy',udsy,'m',timehrs,0)
        endif
      endif
#endif      
      if(write_velres)then !Velocity Residuals
        call writescalh5(outlist(2)%afile,apath,'Vx_Norm_Res',rsu,'none',timehrs,0)   
        call writescalh5(outlist(2)%afile,apath,'Vy_Norm_Res',rsv,'none',timehrs,0)  
        !vecx = log10(1.0+rsu)
        !vecy = log10(1.0+rsv)
        !call writescalh5(outlist(2)%afile,apath,'Vx_Norm_Res_Scaled',vecx,'none',timehrs,0)   
        !call writescalh5(outlist(2)%afile,apath,'Vy_Norm_Res_Scaled',vecy,'none',timehrs,0)
!!        call writevech5(outlist(2)%afile,apath,'Vx_Gradient',dux,duy,'1/s',timehrs,0)
!!        val=sqrt(dux**2+duy**2)
!!        call writescalh5(outlist(2)%afile,apath,'Vx_Gradient_Mag',val,'1/s',timehrs,0)
!!        call writevech5(outlist(2)%afile,apath,'Vy_Gradient',dvx,dvy,'1/s',timehrs,0)  
!!        val=sqrt(dvx**2+dvy**2)
!!        call writescalh5(outlist(2)%afile,apath,'Vy_Gradient_Mag',val,'1/s',timehrs,0) 
      endif
    endif
   
    !EDDY VISCOSITY GROUP  !Alex, BUG FIX, changed 4 to 3, Aug 24, 2009
    if(check_time_list(3))then 
      call print_output_header(header)
      call writescalh5(outlist(3)%afile,apath,'Eddy_Viscosity',vis,'m^2/s',timehrs,0)  
    endif    
    
    !SEDIMENT TRANSPORT GROUP
    if(check_time_list(4))then
      call print_output_header(header)      
      !Total-load sediment concentration, capacity and transports          
      if(write_conc)then       
        call writescalh5(outlist(4)%afile,apath,'Concentration',Ct,'kg/m^3',timehrs,0)  
      endif
      if(isedmodel/=3 .and. write_capac)then      
        call writescalh5(outlist(4)%afile,apath,'Capacity',Ctstar,'kg/m^3',timehrs,0)
      endif  
      if(write_concfrac)then  
        do ks=1,nsed
          write(apbk,62) ks
          aname='Concentration' // apbk
          call writescalh5(outlist(4)%afile,apath,aname,Ctk(:,ks),'kg/m^3',timehrs,0)            
          if(isedmodel/=3)then   
            aname='Capacity' // apbk
            call writescalh5(outlist(4)%afile,apath,aname,Ctkstar(:,ks),'kg/m^3',timehrs,0)
          endif
        enddo !ks
      endif
      call writevech5(outlist(4)%afile,apath,'Total_Sediment_Transport',qtx,qty,'kg/m/s',timehrs,0)
      if(write_fracsusp)then
        call writescalh5(outlist(4)%afile,apath,'Fraction_Suspended',rs,'none',timehrs,0)  
        !do ks=1,nsed
        !  name = 'Fraction_Suspended_' // aconc(ks)
        !  call writescalh5(outlist(11)%afile,apath,aname,rsk(:,ks),'none',timehrs,0)                      
        !  aname = 'Fraction_Mix_Layer_' // trim(aconc(ks))
        !  call writescalh5(outlist(11)%afile,apath,aname,pbk(:,ks,1),'none',timehrs,1) 
        !enddo !ks
      endif
      if(debug_mode)then
        call writescalh5(outlist(4)%afile,apath,'Conc_Res',rsCtkmax,'-',timehrs,0)  
      endif
      if(write_alphatot)then
        call writescalh5(outlist(4)%afile,apath,'alphat',alphat,'',timehrs,0) 
      endif
      
#ifdef DEV_MODE
      if(write_lengthsusp)then
        call writescalh5(outlist(4)%afile,apath,'vLsus',vLsus,'m',timehrs,0) 
      endif
      if(write_lengthtot)then
        call writescalh5(outlist(4)%afile,apath,'vLtot',vLtot,'m',timehrs,0) 
      endif
      if(write_betatot)then
        call writescalh5(outlist(4)%afile,apath,'Beta_t',btk(:,1),'',timehrs,0) 
        !do ks=1,nsed
          !aname='Beta_t'//aconc(ks)
          !call writescalh5(outlist(11)%afile,apath,aname,btk(:,ks),'',timehrs,0)
        !enddo !ks
      endif
      if(write_betasusp)then
        call writescalh5(outlist(4)%afile,apath,'Beta_s',bsk(:,1),'',timehrs,0) 
      endif
      if(write_bedvel)then
        call writescalh5(outlist(4)%afile,apath,'u_b',ubk(:,1),'m/s',timehrs,0)
      endif
#endif
      
      if(wavesedtrans)then
        do i=1,ncells           
          val(i)=sum(Qws(i,:))
          vecx(i)=Wunitx(i)*val(i)
          vecy(i)=Wunity(i)*val(i)
        enddo 
        call writescalh5(outlist(4)%afile,apath,'Wave_Sed_Trans_Onshore',val,'m^2/s',timehrs,0)
        call writevech5(outlist(4)%afile,apath,'Wave_Sed_Trans',vecx,vecy,'m^2/s',timehrs,0)
      endif  
    endif
    
    !MORPHOLOGY GROUP
    if(check_time_list(5))then 
      call print_output_header(header)
      val=-zb
      call writescalh5(outlist(5)%afile,apath,'Depth',val,'m',timehrs,1)
      if(write_morph)then   
        val=zb-zb0
        call writescalh5(outlist(5)%afile,apath,'Morphology_Change',val,'m',timehrs,1)
      endif   
    endif
    
    !SALINITY GROUP
    if(check_time_list(6))then       
      call writescalh5(outlist(6)%afile,apath,'Salinity',sal,'ppt',timehrs,0)
    endif
    
    !TEMPERATURE GROUP
    if(check_time_list(14))then
      call writescalh5(outlist(14)%afile,apath,'Temperature',heat,'c',timehrs,0)
    endif
    
    !WAVE GROUP
    if(check_time_list(7))then 
      call print_output_header(header)            
      call writescalh5(outlist(7)%afile,apath,'Wave_Height',Whgt,'m',timehrs,0)         
      call writescalh5(outlist(7)%afile,apath,'Wave_Period',Wper,'s',timehrs,0)      
      vecx=Whgt*Wunitx
      vecy=Whgt*Wunity
      call writevech5(outlist(7)%afile,apath,'Wave_Height_Vec',vecx,vecy,'m',timehrs,0)
      if(write_wavbrkdiss)then
        call writescalh5(outlist(7)%afile,apath,'Wave_Dissipation',wavediss,'none',timehrs,0)   
      endif
      !!****************** for testing ******************************************************
      !do i=1,ncells
      !  val(i) = h(i)*sqrt(dpx(i)*dpx(i)+dpy(i)*dpy(i))
      !enddo
      !call writescalh5(outlist(7)%afile,apath,'hdpxy',val,'s',timehrs,0)  
      !!****************** for testing ******************************************************
      !Fraction of broken waves and dissipation coefficient
      !if(write_fracbreak .or. write_wavdisscoef)then
        !call wavebreaking(Whgt,Wper,Wang,h,u,v,vecx,vecy) !Qb=vecx,cbr=vecy
      !endif  
      !if(write_fracbreak)then
      !  call writescalh5(outlist(7)%afile,apath,'Frac_Broken_Waves',vecx,'-',timehrs,0)
      !endif
      !if(write_wavbrkind)then
      !  call writescalh5(outlist(7)%afile,apath,'Wave_Breaking',waveibr,'none',timehrs,0)  
      !endif
      !if(write_wavdisscoef)then
      !  call writescalh5(outlist(7)%afile,apath,'Dissipation_Coeff',vecy,'-',timehrs,0)   
      !endif
      if(write_wavstress)then
        call writevech5(outlist(7)%afile,apath,'Wave_Rad_Str',wavestrx,wavestry,'none',timehrs,0)
        do i=1,ncells
          val(i) = sqrt(wavestrx(i)*wavestrx(i) + wavestry(i)*wavestry(i))
        enddo
        call writescalh5(outlist(7)%afile,apath,'Wave_Rad_Str_Mag',val,'none',timehrs,0)
        !!**************** For testing ******************************************************************
        !call writevech5(outlist(7)%afile,apath,'Wave_Rad_Str2',wavstrx,wavstry,'none',timehrs,0)
        !do i=1,ncells
        !  val(i) = sqrt(wavstrx(i)*wavstrx(i) + wavstry(i)*wavstry(i))
        !enddo
        !call writescalh5(outlist(7)%afile,apath,'Wave_Rad_Str_Mag2',val,'none',timehrs,0)   
        !!**************** For testing ******************************************************************
      endif      
      if(write_wavorbvel)then
        call writescalh5(outlist(7)%afile,apath,'Wave_Bot_Orb_Vel',Worb,'m/s',timehrs,0)  
      endif
      if(write_wavlen)then
        call writescalh5(outlist(7)%afile,apath,'Wave_Length',Wlen,'m',timehrs,0)
      endif        
      if(write_wavvel)then
        call writevech5(outlist(7)%afile,apath,'Wave_Flux_Velocity',us,vs,'m/s',timehrs,0)
        if(write_wavvelmag)then
          do i=1,ncells
            val(i) = sqrt(us(i)*us(i) + vs(i)*vs(i))
          enddo 
          call writescalh5(outlist(7)%afile,apath,'Wave_Flux_Velocity_Magnitude',val,'m/s',timehrs,0)
        endif    
      endif
      if(write_rollenergy)then
        val = Ssr/2.0
        call writescalh5(outlist(7)%afile,apath,'Roller_Energy',val,'',timehrs,0)
      endif
      if(write_rollstress)then
         call interp_vec_wavfl(rxrs,ryrs,vecx,vecy) !Extrapolate to zero
         call writevech5(outlist(7)%afile,apath,'Roller_Stress',vecx,vecy,'m/s',timehrs,0)   
      endif
      if(write_rolldiss)then
        call interp_scal_wavfl(roldiss,val) !Extrapolate to zero
        call writescalh5(outlist(7)%afile,apath,'Roller_Dissipation',val,'',timehrs,0)  
      endif
    endif      
    
    !METEO GROUP
    if(check_time_list(8))then 
      if(windvar .or. windsta)then
        call print_output_header(header)
        if(windformat/=7)then
          call writevech5(outlist(8)%afile,apath,'Wind_Velocity',uwind,vwind,'m/s',timehrs,1) 
          do i=1,ncells
            val(i) = sqrt(uwind(i)*uwind(i) + vwind(i)*vwind(i))
          enddo
          if(write_wndmag)then
            val(ncells+1:ncellsD)=0.0
            call writescalh5(outlist(8)%afile,apath,'Wind_Magnitude',val,'m/s',timehrs,1)
          endif
        endif
        if(write_wndstress)then
          call writevech5(outlist(8)%afile,apath,'Wind_Stress',tauwindx,tauwindy,'N/m^2',timehrs,1)
        endif
        if(write_wndstressmag)then
          do i=1,ncells
            val(i) = sqrt(tauwindx(i)*tauwindx(i) + tauwindy(i)*tauwindy(i))
          enddo
          call writescalh5(outlist(8)%afile,apath,'Wind_Stress_Magnitude',val,'N/m^2',timehrs,1)  
        endif
      elseif(windconst)then
        call print_output_header(header)
        vecx = wndx; vecy = wndy               
        call writevech5(outlist(8)%afile,apath,'Wind_Velocity',vecx,vecy,'m/s',timehrs,1)          
        val = sqrt(wndx*wndx + wndy*wndy)  !scalar to vector 
        call writescalh5(outlist(8)%afile,apath,'Wind_Magnitude',val,'m/s',timehrs,1)        
      endif
      if(write_atmpres)then
        call writescalh5(outlist(8)%afile,apath,'Atm_Pressure',pressatm,'Pa',timehrs,1)
      endif
      if(write_atmpresgrad)then
        call writescalh5(outlist(8)%afile,apath,'Atm_Pressure_GradX',pressatmdx,'Pa/m',timehrs,1)
        call writescalh5(outlist(8)%afile,apath,'Atm_Pressure_GradY',pressatmdy,'Pa/m',timehrs,1)    
      endif
    endif  
    
    !BED COMPOSITION GROUP
    if(check_time_list(10))then      
      !Sediment Percentiles   
      val=1000.0*d50 !Convert from m to mm
      call writescalh5(outlist(10)%afile,apath,'D50',val,'mm',timehrs,1)   
      if(outperdiam(ipd(90)))then
        val=1000.0*d90 !Convert from m to mm
        call writescalh5(outlist(10)%afile,apath,'D90',val,'mm',timehrs,1)     
      endif    
      do ii=1,nperdiam
        if(iper(ii)==50 .or. iper(ii)==90) cycle !Already output above  
        if(outperdiam(ii))then
          call sedpercentile(iper(ii),val)
          val=1000.0*val !Convert from m to mm
4324 format('D',I02)
          write(aperdiam,4324) iper(ii)
          call writescalh5(outlist(10)%afile,apath,aperdiam,val,'mm',timehrs,1)  
        endif
      enddo
      if(write_thickness .or. write_sizefrac)then
        do j=1,nlay
          if(j<=9)then
            write(alay,71) j
          else
            write(alay,72) j
          endif
          if(write_thickness)then
            aname = 'Thickness' // trim(alay)
            call writescalh5(outlist(10)%afile,apath,aname,db(:,j),'m',timehrs,1)
          endif
          if(write_sizefrac)then
            do ks=1,nsed
              write(apbk,62) ks  
              aname = 'Fraction' // trim(apbk) // trim(alay)
              call writescalh5(outlist(10)%afile,apath,aname,pbk(:,ks,j),'none',timehrs,1)              
            enddo !ks
          endif
        enddo !j
      endif
    endif
    
!    !FRACTIONAL SEDIMENT GROUP (Useful for hotstart)
!    if(check_time_list(11))then
!      if(write_fracsusp)then  
!        do ks=1,nsed
!          write(apbk,62) ks
!          aname='Concentration' // apbk
!          call writescalh5(outlist(11)%afile,apath,aname,Ctk(:,ks),'kg/m^3',timehrs,0)            
!          if(isedmodel/=3)then   
!            aname='Capacity' // apbk
!            call writescalh5(outlist(11)%afile,apath,aname,Ctkstar(:,ks),'kg/m^3',timehrs,0)
!          endif            
!!        aname = 'Fraction_Suspended_' // aconc(ks)
!!        call writescalh5(outlist(11)%afile,apath,aname,rsk(:,ks),'none',timehrs,0)                      
!!        aname = 'Fraction_Mix_Layer_' // trim(aconc(ks))
!!        call writescalh5(outlist(11)%afile,apath,aname,pbk(:,ks,1),'none',timehrs,1)                   
!!        aname='Beta_t'//aconc(ks)
!!        call writescalh5(outlist(11)%afile,apath,aname,btk(:,ks),'',timehrs,0)      
!        enddo !ks
!      endif
!      
!    endif
    
    !WAVE DETAILS GROUP (Outdated)
!    if(check_time_list(12))then      
!    endif  
    
    !ROUGHNESS GROUP
    if(check_time_list(13))then    
      if(write_meanbedshear)then
        call writescalh5(outlist(13)%afile,apath,'Bed_Shear_Stress_Mean',bsxy,'Pa',timehrs,0)
      endif
      if(write_normapprough)then
        do i=1,ncells
          val(i) = fric_normapprough(uelwc(i),uv(i),cfrict(i)) !Normalized apparent roughness length
        enddo  
        call writescalh5(outlist(13)%afile,apath,'Norm_App_Rough',val,'-',timehrs,0)  
      endif
      if(write_normrough)then
        do i=1,ncells
          val(i) = z0(i)/h(i) !Normalized roughness length
        enddo  
        call writescalh5(outlist(13)%afile,apath,'Norm_Rough',val,'-',timehrs,0)  
      endif
      if(write_rough)then !Total roughness height
        do i=1,ncells
          val(i) = z0(i)*30.0 !Roughness height
        enddo  
        call writescalh5(outlist(13)%afile,apath,'Roughness_Height_Total',val,'-',timehrs,0) 
      endif
      if(write_roughrip)then
        call writescalh5(outlist(13)%afile,apath,'Roughness_Height_Ripple',rksr,'-',timehrs,0) 
      endif
      if(write_roughmegarip)then
        call writescalh5(outlist(13)%afile,apath,'Roughness_Height_Megaripple',rksmr,'-',timehrs,0) 
      endif
      if(write_roughdune)then
        call writescalh5(outlist(13)%afile,apath,'Roughness_Height_Dune',rksd,'-',timehrs,0) 
      endif
      if(write_roughgrain)then
        call writescalh5(outlist(13)%afile,apath,'Roughness_Height_Grain',rksg,'-',timehrs,0) 
      endif
      if(write_roughtrans)then
        call writescalh5(outlist(13)%afile,apath,'Roughness_Height_Transport_Current',rkstc,'-',timehrs,0) 
      endif
      if(write_mannings)then
        call writescalh5(outlist(13)%afile,apath,'Mannings_Coefficient',coefman,'-',timehrs,0)   
      endif
    endif
#endif

    return
    end subroutine write_output_xmdf
    
!***********************************************************************    
    subroutine write_output_sup(header)
! writes the SMS Super ASCII Output Files
! written by Alex Sanchez, USACE-CHL
!***********************************************************************
#include "CMS_cpp.h"
    use comvarbl, only: flowpath,casename,timehrs,reftime
    use der_def, only: gow
    use der_lib, only: curvxy
    use diag_def, only: debug_mode
    use flow_def, only: iwet,h,eta,p,uv,u,v,uv,vis,us,vs,rsp,rsu,rsv,&
       hpred,upred,vpred,dux,duy,dvx,dvy,dpx,dpy,hmin,waveflux,dppx,dppy,pp,&
       Hu,Hv,apuareap,sp,spu,spv,sumu,gravinv
    use flow_lib, only: streamwise_curvature
    use fric_def, only: bsxy,cfrict,uelwc,z0
    use fric_lib, only: fric_normapprough
    use geo_def, only: zb,zb0
    use met_def, only: windvar,windconst,uwind,vwind,rain_evap,wndx,wndy,windsta,windformat,&
        nwnd_inc,wndvalsx,wndvalsy,wndtimes,presconst,presvar,pressatm,tauwindx,tauwindy
    use out_def
    use out_lib
    use prec_def
#ifdef DEV_MODE
    use q3d_def, only: q3d,f3dxx,f3dxy,f3dyy,f3du,f3dv,wavcurint,&
      udsx,udsy,uzplay,vzplay,nzplay,q3d_lay
#endif
    use rol_def, only: roller,rxrs,ryrs,roldiss
    use sal_def, only: sal
    use heat_def, only: heat
    use sed_def, only: nsed,nlay,isedmodel,singlesize,&
        d50,d90,Ct,Ctk,Ctstar,Ctkstar,aconc,Qtx,Qty,pbk,&
        db,wavesedtrans,Qws,rs,nperdiam,outperdiam,iper,ipd
    use size_def
    use wave_flowgrid_def, only: wunitx,wunity,whgt,wper,Wang,&
        !wavediss,waveibr,Worb,wlen,wavestrx,wavestry,Ssr,wavstrx,wavstry
        wavediss,Worb,wlen,wavestrx,wavestry,Ssr,wavstrx,wavstry
    implicit none    
    integer :: i,ii,j,ks,nn,nd,nc
    real(ikind) :: val(ncellsD),vecx(ncellsD),vecy(ncellsD)
    real(8) :: timed
    character(len=5) :: alay
    character(len=3) :: aperdiam
    character(len=200) :: apath,aname,outdirpath
    logical :: header,check_time_list
    
71  format(1x,'(',I1,')')
72  format(1x,'(',I2,')')
    
    nn = len_trim(simlabel)
    apath = simlabel(1:nn)//'/'
    outdirpath = 'ASCII_Solutions'
    
    !aname = trim(outdirpath) // '/' // trim(flowpath) // trim(casename)
    aname = trim(flowpath) // trim(outdirpath) // '/' //  trim(casename)
    
    if(debug_mode)then
      nd = ncellsD
    else
      nd = ncells
    endif
    nc = nd
    
    timed = timehrs !Double precision
    
    !WATER SURFACE ELEVATION GROUP
    if(check_time_list(1))then             
      call print_output_header(header)
      eta(ncells+1:ncellsD)=p(ncells+1:ncellsD)*gravinv
      call write_scal_dat_file(aname,'Water_Elevation','eta',eta) !SUPER ASCII File      
      if(write_presres)then !Pressure Residuals
        call write_scal_dat_file(aname,'Pres_Norm_Res','rsp',rsp) !SUPER ASCII File
      endif
      if(write_pres)then
        call write_scal_dat_file(aname,'Pressure','p',p) !SUPER ASCII File
      endif
      if(write_totdep)then
        call write_scal_dat_file(aname,'Total_Water_Depth','h',h) !SUPER ASCII File
      endif
      if(write_wsegrad)then
        vecx=dpx*gravinv
        vecy=dpy*gravinv
        call write_scal_dat_file(aname,'detax','detax',vecx) !SUPER ASCII File
        call write_scal_dat_file(aname,'detay','detay',vecy) !SUPER ASCII File
      endif
      if(write_presgrad)then
        call write_scal_dat_file(aname,'dpx','dpx',dpx) !SUPER ASCII File
        call write_scal_dat_file(aname,'dpy','dpy',dpy) !SUPER ASCII File  
      endif
      if(write_prescor)then  
        call write_scal_dat_file(aname,'Pressure_Correction','pp',pp) !SUPER ASCII File
      endif  
      if(write_prescorgrad)then 
        call write_scal_dat_file(aname,'dppx','dppx',dppx) !SUPER ASCII File
        call write_scal_dat_file(aname,'dppy','dppy',dppy) !SUPER ASCII File
      endif
      !if(ncellpoly>0 .and. debug_mode)then
      !    call write_scal_dat_file(aname,'Hu','Hu',Hu) !SUPER ASCII File
      !    call write_scal_dat_file(aname,'Hv','Hv',Hv) !SUPER ASCII File
      !    call write_scal_dat_file(aname,'apuareap','apuareap',apuareap) !SUPER ASCII File
      !    call write_scal_dat_file(aname,'sumu','sumu',sumu) !SUPER ASCII File
      !    call write_scal_dat_file(aname,'spu','spu',spu) !SUPER ASCII File
      !    call write_scal_dat_file(aname,'spv','spv',spv) !SUPER ASCII File
      !endif
    endif
    
    !CURRENT VELOCITY GROUP
    if(check_time_list(2))then 
      call print_output_header(header)
      if(waveflux)then !****************************
        vecx=u-us; vecy=v-vs !Current (Eulerian) velocities
        call write_vec_dat_file(aname,'Current_Velocity','vel',vecx,vecy) !SUPER ASCII File
        if(write_veltotal)then
          call write_vec_dat_file(aname,'Total_Velocity','totvel',u,v) !SUPER ASCII File
        endif          
      else  
        call write_vec_dat_file(aname,'Current_Velocity','vel',u,v) !SUPER ASCII File
      endif
      if(write_velgrad)then
        call write_scal_dat_file(aname,'dux','dux',dux) !SUPER ASCII File
        call write_scal_dat_file(aname,'duy','duy',duy) !SUPER ASCII File
        call write_scal_dat_file(aname,'dvx','dvx',dvx) !SUPER ASCII File
        call write_scal_dat_file(aname,'dvy','dvy',dvy) !SUPER ASCII File
      endif
      if(write_velcurv)then
        call curvxy(gow,dux,duy,vecx,vecy)
        call write_scal_dat_file(aname,'d2udx2','d2udx2',vecx) !SUPER ASCII File
        call write_scal_dat_file(aname,'d2udy2','d2udy2',vecy) !SUPER ASCII File
        call curvxy(gow,dvx,dvy,vecx,vecy)
        call write_scal_dat_file(aname,'d2vdx2','d2vdx2',vecx) !SUPER ASCII File
        call write_scal_dat_file(aname,'d2vdy2','d2vdy2',vecy) !SUPER ASCII File
      endif
      if(write_streamcurv)then
        do i=1,ncells
          vecx(i)=u(i)-us(i)
          vecy(i)=v(i)-vs(i)
          val(i) = streamwise_curvature(vecx(i),vecy(i),uv(i),duy(i),duy(i),dvx(i),dvy(i))
        enddo
        call write_scal_dat_file(aname,'curv','curv',val) !SUPER ASCII File
      endif
      if(write_velres)then !Velocity Residuals
        call write_scal_dat_file(aname,'Vx_Norm_Res','rsu',rsu) !SUPER ASCII File
        call write_scal_dat_file(aname,'Vy_Norm_Res','rsv',rsv) !SUPER ASCII File
        !vecx = log10(1.0+rsu)
        !vecy = log10(1.0+rsv)
        !call writescalh5(outlist(2)%afile,apath,'Vx_Norm_Res_Scaled',vecx,'none',timehrs,0)   
        !call writescalh5(outlist(2)%afile,apath,'Vy_Norm_Res_Scaled',vecy,'none',timehrs,0)
!!        call writevech5(outlist(2)%afile,apath,'Vx_Gradient',dux,duy,'1/s',timehrs,0)
!!        val=sqrt(dux**2+duy**2)
!!        call writescalh5(outlist(2)%afile,apath,'Vx_Gradient_Mag',val,'1/s',timehrs,0)
!!        call writevech5(outlist(2)%afile,apath,'Vy_Gradient',dvx,dvy,'1/s',timehrs,0)  
!!        val=sqrt(dvx**2+dvy**2)
!!        call writescalh5(outlist(2)%afile,apath,'Vy_Gradient_Mag',val,'1/s',timehrs,0) 
      endif
    endif
   
    !EDDY VISCOSITY GROUP  !Alex, BUG FIX, changed 4 to 3, Aug 24, 2009
    if(check_time_list(3))then 
      call print_output_header(header)
      call write_scal_dat_file(aname,'Eddy_Viscosity','vis',vis)
    endif    
    
    !SEDIMENT GROUP
    if(check_time_list(4))then
      call print_output_header(header)      
      !Total-load sediment concentration, capacity and transports          
      if(write_conc)then
        call write_scal_dat_file(aname,'Concentration','Ct',Ct)
      endif
      if(isedmodel/=3 .and. write_capac)then
        call write_scal_dat_file(aname,'Capacity','Ctstar',Ctstar)
      endif  
      call write_vec_dat_file(aname,'Total_Sediment_Transport','Qt',qtx,qty) !SUPER ASCII File
      if(write_fracsusp)then 
        call write_scal_dat_file(aname,'Fraction_Suspended','rs',rs)
      endif
!      call writescalh5(outlist(4)%afile,apath,'alphat',alphat,'',timehrs,0) 
!      call writescalh5(outlist(4)%afile,apath,'vLsus',vLsus,'',timehrs,0) 
!      call writescalh5(outlist(4)%afile,apath,'vLtot',vLtot,'',timehrs,0) 
      
!      val=0.0
!      do ks=1,nsed
!        val=val+btk(:,ks)
!      enddo
!      val=val/float(nsed)
!      call writescalh5(outlist(4)%afile,apath,'Beta_t',btk(:,1),'',timehrs,0) 
!      call writescalh5(outlist(4)%afile,apath,'Beta_s',bsk(:,1),'',timehrs,0) 
!      call writescalh5(outlist(4)%afile,apath,'u_b',ubk(:,1),'',timehrs,0)        
    endif
    
    !MORPHOLOGY GROUP
    if(check_time_list(5))then 
      call print_output_header(header)
      val=-zb
      call write_scal_dat_file(aname,'Depth','dep',val)
      if(write_morph)then
        val=zb-zb0  
        call write_scal_dat_file(aname,'Morphology_Change','dzb',val)           
      endif
    endif
    
    !SALINITY GROUP
    if(check_time_list(6))then
      call write_scal_dat_file(aname,'Salinity','sal',sal)
    endif
    
    !TEMPERATURE GROUP
    if(check_time_list(14))then
      call write_scal_dat_file(aname,'Temperature','temp',heat)
    endif
    
    !WAVE GROUP
    if(check_time_list(7))then 
      call print_output_header(header)
      call write_scal_dat_file(aname,'Wave_Height','Hs',Whgt)
      call write_scal_dat_file(aname,'Wave_Period','Tp',Wper)
      vecx=Whgt*Wunitx
      vecy=Whgt*Wunity
      call write_vec_dat_file(aname,'Wave_Vector','wav',vecx,vecy)           
      if(write_wavbrkdiss)then
        call write_scal_dat_file(aname,'Wave_Dissipation','wdiss',wavediss)
      endif
      !!****************** for testing ******************************************************
      !do i=1,ncells
      !  val(i) = h(i)*sqrt(dpx(i)*dpx(i)+dpy(i)*dpy(i))
      !enddo
      !call writescalh5(outlist(7)%afile,apath,'hdpxy',val,'s',timehrs,0)  
      !!****************** for testing ******************************************************
      !Fraction of broken waves and dissipation coefficient
      !if(write_fracbreak .or. write_wavdisscoef)then
        !call wavebreaking(Whgt,Wper,Wang,h,u,v,vecx,vecy) !Qb=vecx,cbr=vecy
      !endif  
      !if(write_fracbreak)then
      !  call writescalh5(outlist(7)%afile,apath,'Frac_Broken_Waves',vecx,'-',timehrs,0)
      !endif
      !if(write_wavbrkind)then
      !  call write_scal_dat_file(aname,'Wave_Breaking','wbrkind',waveibr)
      !endif
      !if(write_wavdisscoef)then
      !  call writescalh5(outlist(7)%afile,apath,'Dissipation_Coeff',vecy,'-',timehrs,0)   
      !endif
      if(write_wavstress)then
        call write_vec_dat_file(aname,'Wave_Stress','rad',wavestrx,wavestry) !SUPER ASCII File
      endif      
      if(write_wavorbvel)then
        call write_scal_dat_file(aname,'Wave_Bot_Orb_Vel','wOrb',Worb)
      endif
      if(write_wavlen)then
        call write_scal_dat_file(aname,'Wave_Length','wLen',Wlen)
      endif        
      if(write_wavvel)then
        call write_vec_dat_file(aname,'Wave_Flux_Velocity','wvel',us,vs) !SUPER ASCII File
      endif
      if(write_rollenergy)then
        val = Ssr/2.0
        call write_scal_dat_file(aname,'Roller_Energy','rolenergy',val)
      endif
      if(write_rollstress)then
         call interp_vec_wavfl(rxrs,ryrs,vecx,vecy) !Extrapolate to zero
         call write_vec_dat_file(aname,'Roller_Stress','rolstress',vecx,vecy) !SUPER ASCII File
      endif
      if(write_rolldiss)then
        call interp_scal_wavfl(roldiss,val) !Extrapolate to zero
        call write_scal_dat_file(aname,'Roller_Dissipation','roldiss',val)
      endif         
    endif      
    
    !METEO GROUP
    if(check_time_list(8))then 
      if(windvar .or. windsta)then
        call print_output_header(header)
        if(windformat/=7)then          
          call write_vec_dat_file(aname,'Wind_Velocity','wind',uwind,vwind) !SUPER ASCII File
        endif
        if(write_wndstress)then
          call write_vec_dat_file(aname,'Wind_Stress','wndstr',tauwindx,tauwindy) !SUPER ASCII File
        endif
        if(write_wndstressmag)then
          do i=1,ncells
            val(i) = sqrt(tauwindx(i)*tauwindx(i) + tauwindy(i)*tauwindy(i))
          enddo
          !call writescalh5(outlist(8)%afile,apath,'Wind_Stress_Magnitude',val,'N/m^2',timehrs,1)  
        endif
      elseif(windconst)then
        call print_output_header(header)
        vecx = wndx; vecy = wndy               
        !call writevech5(outlist(8)%afile,apath,'Wind_Velocity',vecx,vecy,'m/s',timehrs,1)          
        val = sqrt(wndx*wndx + wndy*wndy)  !scalar to vector 
        !call writescalh5(outlist(8)%afile,apath,'Wind_Magnitude',val,'m/s',timehrs,1)        
      endif
      if(write_atmpres)then
        val = pressatm/100.0 !Convert from Pascal to mbar
        !call writescalh5(outlist(8)%afile,apath,'Atm_Pressure',val,'mbar',timehrs,1)       
      endif
    endif  
    
    !BED GROUP                 !These were all previously commented out - meb 03/02/2018
    !if(check_time_list(10))then      
    !  !Sediment Percentiles   
    !  val=1000.0*d50 !Convert from m to mm
    !  !call writescalh5(outlist(10)%afile,apath,'D50',val,'mm',timehrs,1)   
    !  if(outperdiam(ipd(90)))then
    !    val=1000.0*d90 !Convert from m to mm
    !    !call writescalh5(outlist(10)%afile,apath,'D90',val,'mm',timehrs,1)     
    !  endif    
    !  do ii=1,nperdiam
    !    if(iper(ii)==50 .or. iper(ii)==90) cycle !Already output above  
    !    if(outperdiam(ii))then
    !      call sedpercentile(iper(ii),val)
    !      val=1000.0*val !Convert from m to mm
    !      4324 format('D',I02)
    !      write(aperdiam,4324) iper(ii)
    !      !call writescalh5(outlist(10)%afile,apath,aperdiam,val,'mm',timehrs,1)  
    !    endif
    !  enddo      
    !endif
    
    !MIXED SEDIMENT GROUP   !These were all previously commented out - meb 03/02/2018
    !if(check_time_list(11))then
      !do ks=1,nsed
        !aname='Concentration_'//aconc(ks)
        !call writescalh5(outlist(11)%afile,apath,aname,Ctk(:,ks),'kg/m^3',timehrs,0)            
        !if(isedmodel/=3)then   
          !aname='Capacity_' // aconc(ks)
          !call writescalh5(outlist(11)%afile,apath,aname,Ctkstar(:,ks),'kg/m^3',timehrs,0)
        !endif            
!        aname = 'Fraction_Suspended_' // aconc(ks)
!        call writescalh5(outlist(11)%afile,apath,aname,rsk(:,ks),'none',timehrs,0)                      
!        aname = 'Fraction_Mix_Layer_' // trim(aconc(ks))
!        call writescalh5(outlist(11)%afile,apath,aname,pbk(:,ks,1),'none',timehrs,1)                   
!        aname='Beta_t'//aconc(ks)
!        call writescalh5(outlist(11)%afile,apath,aname,btk(:,ks),'',timehrs,0)      
      !enddo !ks
      !do j=1,nlay
      !  if(j<=9)then
      !    write(alay,71) j
      !  else
      !    write(alay,72) j
      !  endif 
      !  aname = 'Thickness' // alay
        !call writescalh5(outlist(12)%afile,apath,aname,db(:,j),'m',timehrs,1)   
      !  do ks=1,nsed               
      !    aname = 'Fraction_' // trim(aconc(ks)) // alay
          !call writescalh5(outlist(12)%afile,apath,aname,pbk(:,ks,j),'none',timehrs,1)              
      !  enddo
      !enddo !j          
    !endif
    
    !WAVE DETAILS GROUP (Outdated)
!    if(check_time_list(12))then      
!    endif  
    
    !ROUGHNESS GROUP                 !These were all previously commented out - meb 03/02/2018
    !if(check_time_list(13))then    
    !  if(write_meanbedshear)then
        !call writescalh5(outlist(13)%afile,apath,'Bed_Shear_Stress_Mean',bsxy,'Pa',timehrs,0)
    !  endif
    !  if(write_normapprough)then
    !    do i=1,ncells
    !      val(i) = fric_normapprough(uelwc(i),uv(i),cfrict(i)) !Normalized apparent roughness
    !    enddo  
        !call writescalh5(outlist(13)%afile,apath,'Norm_App_Rough',val,'-',timehrs,0)  
    !  endif
    !  if(write_normrough)then
    !    do i=1,ncells
    !      val(i) = z0(i)/h(i) !Normalized roughness
    !    enddo  
        !call writescalh5(outlist(13)%afile,apath,'Norm_Rough',val,'-',timehrs,0)  
    !  endif
    !endif
    
    return
    endsubroutine write_output_sup
    
!***********************************************************************    
    logical function check_time_list(id)
! Checks if the current simulation time matches an output list time
!    
! Alex Sanchez, USACE-ERDC-CHL    
!***********************************************************************
    use comvarbl, only: timehrs,ctime,dtime,nfsch
    use out_def, only: outlist
    use prec_def
    implicit none
    integer :: id
    real :: listtime
    real(ikind) :: toltime 
  
    check_time_list = .false.    
    if(.not.outlist(id)%write_dat) return

    if(nfsch==0)then  !STACK:
      toltime = dtime/(3600.0*3.0)  !implicit
    else           !STACK:
      toltime = dtime/(3600.0*2.0)  !explicit
    endif          !STACK:   

    do while(timehrs>(outlist(id)%times(outlist(id)%inc))+toltime)
      outlist(id)%inc = outlist(id)%inc + 1        
    enddo
    listtime = outlist(id)%times(outlist(id)%inc)
    if(abs(timehrs-listtime)<toltime)then !write time     
      timehrs = listtime
      ctime = listtime*3600.0
      check_time_list = .true.
    endif    
   
    return
    endfunction check_time_list

!**********************************************************************    
    subroutine print_output_header(header)
!**********************************************************************   
    use diag_lib
    use comvarbl, only: timesecs
    use sed_def,  only: scalemorph
    use time_lib
    implicit none
    logical :: header    
    character(len=200) :: msg,msg2
    character(len=50) :: str
    
! meb 02/06/2019
! Added another write statement for effective morphologic time if Morphology Acceleration Factor > 1
    
701 format(' WRITING GLOBAL OUTPUT - ',A)
702 format(' WRITING GLOBAL OUTPUT -        ',A)
703 format(' - EFFECTIVE MORPHOLOGIC TIME - ',A)   
    
    if(.not.header)then
      if (scalemorph .eq. 1) then
        call time_sec2str(timesecs,str)  
        write(msg,701) trim(str)
        call diag_print_message(' ',msg)
      else
        call time_sec2str(timesecs,str)  
        write(msg,702) trim(str)
        call time_sec2str(timesecs*scalemorph,str)  
        write(msg2,703) trim(str)
        call diag_print_message(' ',msg,msg2)
      endif
      header = .true.
    endif
    
    return
    endsubroutine print_output_header

!**********************************************************************    
    subroutine write_debug
! Writes CMS-Flow Debug Information
! Alex Sanchez, USACE-ERDC-CHL    
!***********************************************************************
#include "CMS_cpp.h"
    use size_def, only: ncellpoly
    use out_def, only: write_sup
    implicit none
    
#ifdef XMDF_IO    
    call write_debug_xmdf !XMDF output diagnostic file
#endif
    if(ncellpoly>0 .and. write_sup)then
      call write_debug_sup  !Super ASCII output diagnostic file
    endif
    
    return
    endsubroutine write_debug
        
!**********************************************************************    
    subroutine write_debug_xmdf
! Writes CMS-Flow Debug Information
! Alex Sanchez, USACE-ERDC-CHL    
!***********************************************************************
#include "CMS_cpp.h"
#ifdef XMDF_IO
    use size_def
    use geo_def, only: zb
    use flow_def
    use comvarbl
    use met_def
    use bnd_def
    use wave_flowgrid_def
    use sed_def
    use cms_def
    use out_def, only: dgoutfile,simlabel,write_sup     
    use out_lib, only: writescalh5,writevech5 
    use prec_def
    implicit none   
    integer :: i,nn
    real(ikind) :: val(ncellsD),vecx(ncellsD),vecy(ncellsD)
    character :: apath*100
    
!!    timehrs = niter + maxit*ntime
    nn=len_trim(simlabel)
    apath=simlabel(1:nn)//'_diag/'
    
    val=0.0; vecx=0.0; vecy=0.0  !Initialize

    !Hydro  
    val = -zb
    call writescalh5(dgoutfile,apath,'Depth_Debug',val,'m',timehrs,1)

    val = iwet
    call writescalh5(dgoutfile,apath,'Wet_Debug',val,'m',timehrs,1)
    call writescalh5(dgoutfile,apath,'Total_Depth_Debug',h,'m',timehrs,1)

    eta=iwet*p*gravinv+(-999.0)*(1-iwet)
    call writescalh5(dgoutfile,apath,'Water_Elevation_Debug',eta,'m',timehrs,0)
    
    call writevech5 (dgoutfile,apath,'Current_Velocity_Debug',u,v,'m/s',timehrs,0)
    call writescalh5(dgoutfile,apath,'dux_Debug',dux,'m/s',timehrs,0)
    call writescalh5(dgoutfile,apath,'duy_Debug',duy,'m/s',timehrs,0)  
    call writescalh5(dgoutfile,apath,'dvx_Debug',dvx,'m/s',timehrs,0)  
    call writescalh5(dgoutfile,apath,'dvy_Debug',dvy,'m/s',timehrs,0)  
    call writescalh5(dgoutfile,apath,'Current_Magnitude_Debug',uv,'m/s',timehrs,0)   
    call writescalh5(dgoutfile,apath,'Eddy_Viscosity_Debug',vis,'m^2/s',timehrs,0)
    call writescalh5(dgoutfile,apath,'Pressure_Debug',p,'m^2/s^2',timehrs,1)
    call writescalh5(dgoutfile,apath,'dpx_Debug',dpx,'m/s^2',timehrs,1)
    call writescalh5(dgoutfile,apath,'dpy_Debug',dpy,'m/s^2',timehrs,1)
    call writescalh5(dgoutfile,apath,'Pressure_Correction_Debug',pp,'m^2/s^2',timehrs,1)
    call writescalh5(dgoutfile,apath,'dppx_Debug',dppx,'m/s^2',timehrs,1)
    call writescalh5(dgoutfile,apath,'dppy_Debug',dppy,'m/s^2',timehrs,1)
    
    !Residuals
    call writescalh5(dgoutfile,apath,'Vx_Norm_Res',rsu,'none',timehrs,1)   
    call writescalh5(dgoutfile,apath,'Vy_Norm_Res',rsv,'none',timehrs,1)   
    call writescalh5(dgoutfile,apath,'Pres_Norm_Res',rsp,'none',timehrs,1)
    
    !Sediment
    !Total-load sediment concentration, capacity and transports        
    if(sedtrans)then
      call writescalh5(dgoutfile,apath,'Concentration_Debug',Ct,'kg/m^3',timehrs,0)  
      if(isedmodel/=3)then      
        call writescalh5(dgoutfile,apath,'Capacity_Debug',Ctstar,'kg/m^3',timehrs,0)        
      endif  
      call writevech5(dgoutfile,apath,'Total_Sediment_Transport_Debug',qtx,qty,'m^2/s',timehrs,0)      
      call writescalh5(dgoutfile,apath,'alphat_Debug',alphat,'',timehrs,0) 
      if(iadapttot>=5)then
        call writescalh5(dgoutfile,apath,'vLsus_Debug',vLsus,'',timehrs,0) 
        call writescalh5(dgoutfile,apath,'vLtot_Debug',vLtot,'',timehrs,0) 
      endif
!      val=0.0
!      do ks=1,nsed
!        val=val+btk(:,ks)
!      enddo
!      val=val/float(nsed)
      call writescalh5(dgoutfile,apath,'Beta_t_Debug',btk(:,1),'',timehrs,0) 
      call writescalh5(dgoutfile,apath,'Beta_s_Debug',bsk(:,1),'',timehrs,0) 
      !call writescalh5(dgoutfile,apath,'u_b_Debug',ubk(:,1),'',timehrs,0)       
!      call writescalh5(dgoutfile,apath,'Fraction_Suspended',rsk(:,1),'none',timehrs,0)      
      if(.not.singlesize)then         
        !Sediment Percentiles   
        val=1000.0*d50 !Convert from m to mm
        call writescalh5(dgoutfile,apath,'D50',val,'mm',timehrs,1)    
        if(bedlay(1)%ipbkinp==2)then  
          call sedpercentile(16,val)
          val=1000.0*val !Convert from m to mm
          call writescalh5(dgoutfile,apath,'D16',val,'mm',timehrs,1)     
          call sedpercentile(84,val)  
          v=1000.0*val !Convert from m to mm
          call writescalh5(dgoutfile,apath,'D84',val,'mm',timehrs,1)               
        else
          call sedpercentile(35,val)
          val=1000.0*val !Convert from m to mm
          call writescalh5(dgoutfile,apath,'D35',val,'mm',timehrs,1)     
          call sedpercentile(90,val)  
          val=1000.0*val !Convert from m to mm
          call writescalh5(dgoutfile,apath,'D90',val,'mm',timehrs,1)      
        endif 
      endif 
    endif      
    
    !Waves  
    if(cmswave)then
      call writescalh5(dgoutfile,apath,'Wave_Height_Debug',Whgt,'m/s',timehrs,0)   
      call writescalh5(dgoutfile,apath,'Wave_Period_Debug',Wper,'s',timehrs,0) 
      do i=1,ncells
        vecx(i) = Whgt(i)*Wunitx(i)
        vecy(i) = Whgt(i)*Wunity(i)
      enddo
      call writevech5(dgoutfile,apath,'Wave_Height_Vec',vecx,vecy,'m',timehrs,0)
      call writevech5(dgoutfile,apath,'Wave_Rad_Stress_Debug',wavestrx,wavestry,'none',timehrs,0)   
      do i=1,ncellsD
        val(i)=sqrt(wavestrx(i)*wavestrx(i)+wavestry(i)*wavestry(i))
      enddo
      call writescalh5(dgoutfile,apath,'Wave_Rad_Str_Mag_Debug',val,'none',timehrs,0)         
    endif      

    !Wind
    if(windconst)then
      vecx=tauwx; vecy=tauwy; 
      call writevech5(dgoutfile,apath,'Wind_Stress_Debug',vecx,vecy,'none',timehrs,0)  
    elseif(windvar)then
      call writevech5(dgoutfile,apath,'Wind_Stress_Debug',tauwindx,tauwindy,'none',timehrs,0)      
    endif        
    
#endif

    return
    endsubroutine write_debug_xmdf

!**********************************************************************    
    subroutine write_debug_sup
! Writes CMS-Flow Debug Information
! Alex Sanchez, USACE-ERDC-CHL    
!***********************************************************************
    use size_def, only: ncellsd
    use geo_def, only: zb
    use flow_def, only: p,dpx,dpy,pp,h,vis,u,v,rsu,rsv,rsp
    use comvarbl, only: flowpath,casename
    use out_def, only: dgoutfile,simlabel,write_sup
    use out_lib, only: write_scal_dat_file,write_vec_dat_file
    use prec_def, only: ikind
    
    implicit none   
    integer :: nn
    real(ikind) :: val(ncellsD),vecx(ncellsD),vecy(ncellsD)
    character :: apath*100,aname*100
    character(len=100) :: outdirpath
    
!!    timehrs = niter + maxit*ntime
    nn=len_trim(simlabel)
    apath=simlabel(1:nn)//'/'
    outdirpath='ASCII_Solutions'
    
    val=0.0; vecx=0.0; vecy=0.0  !Initialize    

    aname = trim(outdirpath) // '/' // trim(flowpath) // trim(casename)      
    call write_scal_dat_file(aname,'Pressure_Debug','pres',p)
    call write_scal_dat_file(aname,'dpx_Debug','dpx',dpx)
    call write_scal_dat_file(aname,'dpy_Debug','dpy',dpy)
    call write_scal_dat_file(aname,'Pressure_Correction_Debug','pp',pp)
    call write_scal_dat_file(aname,'Depth_Debug','depth',-zb)
    call write_scal_dat_file(aname,'Total_Depth_Debug','depth',h)
    call write_scal_dat_file(aname,'Eddy_Viscosity_Debug','vis',vis)
    call write_vec_dat_file(aname,'Current_Velocity_Debug','vel',u,v)
    call write_scal_dat_file(aname,'Vx_Norm_Res','rsu',rsu)
    call write_scal_dat_file(aname,'Vy_Norm_Res','rsv',rsv)
    call write_scal_dat_file(aname,'Pres_Norm_Res','rsp',rsp)     
    
    return
    endsubroutine write_debug_sup        

!**********************************************    
    subroutine read_time_list(n)
! Reads the ouput time list
! written by Alex Sanchez    
!**********************************************    
    use out_def, only: outseries
    use diag_lib, only: diag_print_warning
    use prec_def
    use unitconv_lib, only: unitconv_var
    implicit none
    integer :: i,id,j,k,n,nlist,nsum,nmax,kth,ierr
    integer, allocatable :: nn(:),inc(:)
    real(8) :: tmin
    real(8), allocatable :: temp(:),tlist(:,:),tseries(:,:)   
    real(ikind) :: fac,con
    character(len=37) :: cardname
    character(len=10) :: fromunits,tounits
    character(len=1000) :: aline
    
    !Get number of list and time information    
    fromunits = 'hrs'

    backspace(77)
    read(77,*) cardname, nlist
    backspace(77)
    allocate(tlist(nlist,3))
    read(77,'(A)',iostat=ierr) aline
    if (cardname == 'TIME_LIST_1' .and. nlist .eq. 0) then   !Always turn on output to Time List 1 - MEB 03/14/2018
      call diag_print_warning('TIME_LIST_1 is turned off','Enabling for minimal hourly output for 720 hrs')
      nlist = 1
      if (allocated(tlist)) deallocate(tlist)
      allocate(tlist(nlist,3))
      tlist(1,1) = 0
      tlist(1,2) = 720
      tlist(1,3) = 1
    else
      read(aline,*,iostat=ierr) cardname, nsum, ((tlist(i,j),j=1,3),i=1,nlist), fromunits
    endif
    
    allocate(nn(nlist))
    allocate(inc(nlist))
    
    tounits = 'hrs'
    call unitconv_var(fromunits,tounits,fac,con)
    
    !Construct time series for each list
    nsum = 0; nmax = 0; nn = 1
    do i=1,nlist
      if(tlist(i,3)>0)then
        nn(i) = int((tlist(i,2)-tlist(i,1))/tlist(i,3)) + 1
      endif
      nsum = nsum + nn(i)       
      nmax = max(nmax,nn(i))      
    enddo    
    allocate(tseries(nlist,nmax))
    tseries = 0.0
    do i=1,nlist
      do j=1,nn(i)
        tseries(i,j) = tlist(i,1) + (j-1)*tlist(i,3)
      enddo
    enddo
       
    !Merge time series     
    if(allocated(temp)) deallocate(temp)
    allocate(temp(0:nsum))
    inc = 1; temp = -1.0; kth = 0        
    do k=1,nsum      
      !Find smallest value of all time series
      tmin = 1.0e25
      do i=1,nlist
        if(inc(i)<=nn(i))then
         if(tseries(i,inc(i))<=tmin)then
           id = i      
           tmin = tseries(i,inc(i))         
         endif   
        endif       
      enddo    
      !Add value to temp and move counter
      if(abs(tseries(id,inc(id))-temp(kth))>0.01*tlist(id,3))then !skip repeated values
        tmin = tseries(id,inc(id))
        kth = kth + 1
        temp(kth) = tseries(id,inc(id))                    
      endif
      inc(id) = inc(id) + 1    
    enddo        
    
    !Copy times to outseries
    outseries(n)%ntimes = kth
    if(allocated(outseries(n)%times)) deallocate(outseries(n)%times)
    allocate(outseries(n)%times(outseries(n)%ntimes))
    do k=1,outseries(n)%ntimes 
      outseries(n)%times(k) = fac*temp(k) + con
    enddo    
    
    deallocate(temp)
    
    return
    endsubroutine read_time_list
    
!**********************************************    
    subroutine read_obs_cells(i)
!**********************************************    
    use out_def, only: obs,obs_cell
    implicit none
    integer, parameter:: ntemp=1000
    integer :: i,n,itemp(ntemp), count
    character(len=37) :: adum(ntemp)
    
    count=0
    do n=1,ntemp
!      read(77,'(I,A)',err=881) itemp(n),adum(n)
!!      read(77,*,err=881) itemp(n),adum(n)
      read(77,*,err=881) itemp(n),adum(n)
      if(itemp(n)==0)then
        backspace(77)
        exit
      endif  
      count=count+1
    enddo       
881 obs(i)%ncells=count
    if(count==0) return
    obs_cell=.true. 
    allocate(obs(i)%cells(count))
    allocate(obs(i)%identifiers(count))
    do n=1,obs(i)%ncells
      obs(i)%cells(n) = itemp(n)
      obs(i)%identifiers(n) = adum(n)
    enddo
    
    return
    endsubroutine read_obs_cells
    
!********************************************************************************
    subroutine obs_cell_init()
! Initializes Observation cell Output
!
! written by Alex, USACE-ERDC-CHL
! last upated Nov 13, 2009 - added sediment and salinity
!********************************************************************************
#include "CMS_cpp.h"
    use geo_def,  only: idmap
    use comvarbl, only: casename,flowpath,dtimebeg
    use hot_def,  only: coldstart
    use sal_def,  only: saltrans
    use heat_def, only: heattrans
    use sed_def,  only: sedtrans
#ifdef DEV_MODE  
    use q3d_def,  only: q3d
#endif    
    use out_def,  only: obs
    implicit none
    
    integer :: i,j,k,nn,npath
    logical :: ok
    
    nn = len_trim(casename)
    npath = len_trim(flowpath) 
    do i=1,4 !groups
      if(obs(i)%active)then
        obs(i)%time_inc = max(obs(i)%time_inc,dtimebeg) !Minimum interval is the initial timestep
        allocate(obs(i)%units(obs(i)%nvar))
        allocate(obs(i)%files(obs(i)%nvar))        
        do j=1,obs(i)%nvar !variables          
          obs(i)%units(j)=50+i*j
          obs(i)%files(j)=casename(1:nn)//"_"//trim(obs(i)%names(j))//".txt"    
          obs(i)%files(j)=flowpath(1:npath)//obs(i)%files(j) !Bug fix, 1/18/11
        enddo
      else
        obs(i)%time_inc=-1.0
      endif
    enddo
      
753 format('% time(hrs)',255I10) !Observation File Header Format
      
    if(coldstart)then  !create data set (new XMDF file created already) 
      do i=1,4 !groups
        if(obs(i)%active)then  
          do j=1,obs(i)%nvar !variables
            if(i==3 .and. j<=3 .and. .not.sedtrans) cycle
            if(i==3 .and. j==4 .and. .not.saltrans) cycle
#ifdef DEV_MODE            
            if(i==1 .and. j==4 .and. .not.q3d) cycle
            if(i==3 .and. j==5 .and. .not.q3d) cycle
#endif
            open(obs(i)%units(j),file=obs(i)%files(j))
            write(obs(i)%units(j),753) (obs(i)%cells(k),k=1,obs(i)%ncells)
            close(obs(i)%units(j)) !closing the file is useful for viewing the results during the simulation
          enddo
        endif
      enddo
    else  !if hotstart - see if data set is created, an if not create it    
      do i=1,4 !groups
        if(obs(i)%active)then  
          do j=1,obs(i)%nvar !variables
            if(i==3 .and. j<=3 .and. .not.sedtrans) cycle
            if(i==3 .and. j==4 .and. .not.saltrans) cycle
#ifdef DEV_MODE           
            if(i==1 .and. j==4 .and. .not.q3d) cycle
#endif            
            inquire(file=obs(i)%files(j),exist=ok) 
            if(ok)then
              open(obs(i)%units(j),file=obs(i)%files(j),access='append')
            else  
              open(obs(i)%units(j),file=obs(i)%files(j))
            endif
            write(obs(i)%units(j),753) (obs(i)%cells(k),k=1,obs(i)%ncells)
            close(obs(i)%units(j)) !closing the file is useful for viewing the results during the simulation
          enddo
        endif
      enddo
    endif
    
    !Map cell id's to new grid
    do i=1,4 !groups
      if(obs(i)%active)then  
        do k=1,obs(i)%ncells
          obs(i)%cells(k)=idmap(obs(i)%cells(k))
        enddo
      endif
    enddo
    
    return
    endsubroutine obs_cell_init
    
!********************************************************************************
    subroutine write_obs_cell
! Writes data to all observation cell files    
! written by Alex Sanchez, USACE-ERDC-CHL
!********************************************************************************
    use size_def, only: ncellsD
    use out_def, only: obs
    use flow_def, only: u,v,eta,h,iwet    
    use sed_def, only: sedtrans,nsed,Ct,rs,qtx,qty,pbk
    use sal_def, only: saltrans,sal
    use heat_def, only: heattrans,heat
    use comvarbl, only: ctime,timesecs
    use const_def, only: eps
    use q3d_def, only: q3d
    !use hot_def, only: coldstart,hstarttime
    use prec_def
    implicit none
    integer :: ks
    real(ikind) :: Csus(ncellsD),qx(ncellsD),qy(ncellsD)

    !**** Needs testing ******************************
    !if(.not.coldstart)then
    !  if(hstarttime >= timehrs)then  !otherwise there could be duplicates written to this file.
    !    return
    !  endif
    !endif    
    
    !MEB - 05/23/2016
    !With very large values of CTIME (value in seconds), the MOD operation fails occasionally.
    !Replaced CTIME with TIMESECS, which is a double precision equivalent.
    
    !Time series
    
    if(obs(1)%active .and. mod(timesecs,obs(1)%time_inc)<0.5)then      
      call write_obs_var(1,1,u)
      call write_obs_var(1,2,v)
      call write_obs_var(1,3,eta)
      if(q3d) call write_obs_velpro !Write vertical velocity profiles
    endif
    
    !Flow
    if(obs(2)%active .and. mod(timesecs,obs(2)%time_inc)<0.5)then      
      qx = u*h
      qy = v*h
      call write_obs_var(2,1,qx)
      call write_obs_var(2,2,qy)
    endif
    
    !Transport
    if(obs(3)%active .and. mod(timesecs,obs(3)%time_inc)<0.5)then       
      if(sedtrans)then
        Csus = iwet*Ct*rs !Suspended sediment concentration
        call write_obs_var(3,1,Csus)
        call write_obs_var(3,2,qtx)
        call write_obs_var(3,3,qty)
        !if(q3d) call write_obs_sedpro !Write vertical velocity profiles
      endif
      if(saltrans) call write_obs_var(3,4,sal)
    endif
    
    !Bed Composition
    if(obs(4)%active .and. mod(timesecs,obs(4)%time_inc)<0.5)then      
      do ks=1,nsed
        call write_obs_var(4,1,pbk(:,ks,1))
      enddo      
    endif
    
    return
    endsubroutine write_obs_cell
        
!********************************************************************************   
    subroutine write_obs_var(i,j,var)
! Writes an observation cell variable to its corresponding file
! written by Alex Sanchez, USACE-ERDC-CHL
!********************************************************************************   
    use size_def, only: ncellsD
    use out_def, only: obs
    use comvarbl, only: timehrs
    use prec_def
    implicit none
    integer :: i,j,k
    real(ikind) :: var(ncellsD)
    
    open(obs(i)%units(j),file=obs(i)%files(j),access='append')
973 format(1x,F11.4,255(1X,1pe11.4))                                                  !Changed 10.4 to 11.4 to increase total time output possible from 11+ years to 114+ years.
    write(obs(i)%units(j),973) timehrs,(var(obs(i)%cells(k)),k=1,obs(i)%ncells)  
    close(obs(i)%units(j)) !closing the file is useful for viewing the results during the simulation
    
    return
    endsubroutine write_obs_var
    
!********************************************************************************   
    subroutine write_obs_velpro
! Writes an observayion cell profile to its corresponding file
! written by Alex Sanchez, USACE-ERDC-CHL
!********************************************************************************   
    use size_def, only: ncellsD
    use flow_def, only: h,u,v,us,vs
    use out_def, only: obs
    use comvarbl, only: timehrs
    use q3d_def, only: nzpsta,zpsta,uzpsta,vzpsta
    use prec_def
    implicit none
    !Internal variables
    integer :: i,j,k,iprofile
    real(ikind) :: uc,vc,vx,vy,wx,wy,zap,zwp
    real(ikind) :: Kc,visvm,visv0,visvslp
    
    open(obs(1)%units(4),file=obs(1)%files(4),access='append')
    !      timehrs,iprofile,h(i),u(i),v(i),us(i),vs(i),uc,vc,vx,vy,wx,wy,zap,zwp,Kc,visvm,visv0,visvslp,(uzpsta(j),vzpsta(j),j=1,nzpsta)
973 format(1x,F9.3,1x,I3,255(1x,1pe11.4)) !Note: nzpsta<=100, 17 + 2*nzpsta = 217 <= 255
    do k=1,obs(1)%ncells
      i = obs(1)%cells(k)
      call q3d_velpro(i,nzpsta,zpsta,iprofile,uc,vc,vx,vy,wx,wy,zap,zwp,&
        Kc,visvm,visv0,visvslp,uzpsta,vzpsta) !Compute velocity vertical profile
      write(obs(1)%units(4),973) timehrs,iprofile,h(i),u(i),v(i),us(i),vs(i),&
        uc,vc,vx,vy,wx,wy,zap,zwp,Kc,visvm,visv0,visvslp,(uzpsta(j),vzpsta(j),j=1,nzpsta)
    enddo  
    close(obs(1)%units(4)) !closing the file is useful for viewing the results during the simulation
    
    return
    endsubroutine write_obs_velpro
    
!!********************************************************************************   
!    subroutine write_obs_suspro
!! Writes an observayion cell profile to its corresponding file
!! written by Alex Sanchez, USACE-ERDC-CHL
!!********************************************************************************   
!    use size_def, only: ncellsD
!    use flow_def, only: h,u,v,us,vs
!    use out_def, only: obs
!    use comvarbl, only: timehrs
!    use q3d_def, only: nzpsta,zpsta,uzpsta,vzpsta
!    use prec_def
!    implicit none
!    !Internal variables
!    integer :: i,j,k,iprofile
!    real(ikind) :: uc,vc,vx,vy,wx,wy,zap,zwp
!    real(ikind) :: Kc,visvm,visv0,visvslp
!    
!    open(obs(1)%units(4),file=obs(1)%files(4),access='append')
!    !      timehrs,iprofile,h(i),u(i),v(i),us(i),vs(i),uc,vc,vx,vy,wx,wy,zap,zwp,Kc,visvm,visv0,visvslp,(uzpsta(j),vzpsta(j),j=1,nzpsta)
!973 format(1x,F9.3,1x,I3,255(1x,1pe11.4)) !Note: nzpsta<=100, 17 + 2*nzpsta = 217 <= 255
!    do k=1,obs(1)%ncells
!      i = obs(1)%cells(k)
!      !call q3d_suspro(i,nzpsta,zpsta,iprofile,uc,vc,vx,vy,wx,wy,zap,zwp,&
!      !  Kc,visvm,visv0,visvslp,uzpsta,vzpsta) !Compute velocity vertical profile
!      write(obs(1)%units(4),973) timehrs,iprofile,h(i),(rsk(i,ks)*Ctk(ks),ks=1,nsed),&
!        (visvk(ks),ks=1,nsed),((czpsta(j,ks),j=1,nzpsta),ks=1,nzpsta)
!    enddo  
!    close(obs(1)%units(4)) !closing the file is useful for viewing the results during the simulation
!    
!    return
!    endsubroutine write_obs_suspro    

!***********************************************************************
    subroutine save_point_init
! Initializes Save Point file output
! written by Mitchell Brown, USACE-ERDC-CHL, Lasted modified May 8, 2012
!***********************************************************************    
    use size_def, only: ncells,ncellpoly
    use out_def
    use diag_def
    use diag_lib
    use geo_def, only: idmap,hproj,vproj,xc,yc,x,y,bathydata,mapid
    use comvarbl, only: casename,flowpath,dtime
    use sal_def, only: saltrans
    use heat_def, only: heattrans
    use sed_def, only: sedtrans
    use hot_def, only: coldstart
    use cms_def, only: cmswave
    implicit none
    integer :: i,j,k,kk,nn,npath,npts,iunit,ksp
    logical :: ok = .false.
    character(len=200) :: spxfile
    real(ikind) :: distmin,dist
    
    nn  = len_trim(casename)
    npath = len_trim(flowpath)
    spxfile=flowpath(1:npath)//casename(1:nn)//".spx"
    do i=1,ngroups
      savept(i)%time_inc=max(savept(i)%time_inc,dtime) !********* IMPORTANT *****************
      allocate(savept(i)%funits(savept(i)%nvar))
      allocate(savept(i)%files(savept(i)%nvar))
      do j=1,savept(i)%nvar 
        !generate file names and unit numbers
        savept(i)%funits(j)=50+i*j
        savept(i)%files(j)=casename(1:nn)//"_"//trim(savept(i)%names(j))//".sp"
        savept(i)%files(j)=flowpath(1:npath)//savept(i)%files(j)  
      enddo
      if(.not.savept(i)%active) savept(i)%time_inc=-1.0
    enddo
    if (.not.sedtrans)then
      savept(2)%active = .false.  !if sediment transport off, do not process
      if(.not.bathydata%ison)then
        savept(5)%active = .false.  !if sediment transport off, do not process
      endif
    endif
    if (.not.saltrans) savept(3)%active = .false.  !if salinity transport off, do not process
    if (.not.cmswave)  savept(4)%active = .false.  !if not steering, do not process

    if(coldstart) then
      hproj=' ' ; vproj=' '       !Initially blank the projections.  
      call cleanup_savept_files   !Since this is coldstart, remove all previous save point files.
      open (100,file=spxfile,status='unknown')  !this file holds the names of all the SP files.
      write(100,'(A)') 'CMS_FLOW_SAVE_PTS_FILENAMES'
      write(100,'(A20,A)') "SAVE_POINT_LABEL    ", '"'//trim(splabel)//'"'
      do i=1,ngroups
        if (savept(i)%active) then
          do j=1,savept(i)%nvar
            !write header information to each file
            npts=savept(i)%ncells
            iunit = savept(i)%funits(j)
            write(100,'(A)') savept(i)%files(j)
            open (iunit,file=savept(i)%files(j),status='UNKNOWN')
            call write_savept_header(iunit,npts,i,j)
            close(iunit)
          enddo
        endif
      enddo
      close(100)
    else     !hotstarted
      do i=1,ngroups
        if (savept(i)%active) then
          do j=1,savept(i)%nvar
            inquire(file=savept(i)%files(j),exist=ok)
            npts=savept(i)%ncells
            iunit = savept(i)%funits(j)
            if (ok) then
              open(unit=iunit,file=savept(i)%files(j),access='APPEND')
            else
              open(unit=iunit,file=savept(i)%files(j),status='NEW')
              call write_savept_header(iunit,npts,i,j)
            endif
            close(iunit)
          enddo
        endif
      enddo
    endif
    
    !Calculate cell ID's
    do i=1,ngroups
      if(savept(i)%active)then
        do j=1,savept(i)%ncells
          !if(savept(i)%cell(j)<=0 .or. savept(i)%cell(j)>ncells)then !Cell ID not specified or invalid
            !Search for nearest cell  
            ksp = 1; distmin = 1.0e10
            if(ncellpoly==0)then
              do k=1,ncells  
                dist = sqrt((savept(i)%x(j)-xc(k))**2+(savept(i)%y(j)-yc(k))**2)
                if(dist<distmin)then
                  ksp = k
                  distmin = dist
                endif
              enddo
              kk = mapid(ksp) !Map to global grid
            else
               do k=1,ncells  
                dist = sqrt((savept(i)%x(j)-x(k))**2+(savept(i)%y(j)-y(k))**2)
                if(dist<distmin)then
                  ksp = k
                  distmin = dist
                endif
               enddo
               kk = ksp !Global index same as local index
            endif
            if(savept(i)%cell(j)>0 .and. savept(i)%cell(j)<=ncells .and. savept(i)%cell(j)/=kk)then
              !write(msg2,*) ' Save point name: ',savept(i)%names(j)
              write(msg3,*) ' Specified Cell ID: ',savept(i)%cell(j)
              write(msg4,*) ' Calculated Cell ID: ',kk
              !write(msg5,*) ' Specified Coordinates: ',savept(i)%x(j),savept(i)%y(j)' m' !Global coordiantes
              write(msg5,*) ' Calculated Coordinates: ',xc(ksp),yc(ksp),' m' !Global coordiantes
              call diag_print_warning('Specific cell ID for Save Point may be incorrect',msg3,msg4,msg5)
            endif
            savept(i)%cell(j) = ksp
          !endif
        enddo
      endif
    enddo
    
    !!map cell id's to new grid
    !do i=1,ngroups
    !  if(savept(i)%active)then
    !    do j=1,savept(i)%ncells
    !      savept(i)%cell(j) = idmap(savept(i)%cell(j))
    !    enddo
    !  endif
    !enddo
    
    return
    
    contains
!**********************************************************************    
      subroutine cleanup_savept_files
! Remove all previous save point files
!
! Author: Mitch Brown, USACE-CHL
!**********************************************************************      
      use out_def, only: savept
      implicit none
      integer :: i,j,ios
      logical :: ok
      
      do i=1,ngroups
        do j=1,savept(i)%nvar
          inquire(file=savept(i)%files(j),exist=ok, iostat=ios)
          if (ok .and. ios==0) then
             open (50,file=savept(i)%files(j))
             close(50,status='delete',err=999)
999          continue
          endif
        enddo
      enddo
      
      return
      endsubroutine cleanup_savept_files
      
!**********************************************************************     
      subroutine write_savept_header(iunit,npts,i,j)
!**********************************************************************      
      use out_def, only: savept
      use geo_def, only: hproj,vproj 
      use comvarbl, only: iyr,imo,iday,ihr,imin,version,revision
      implicit none
      integer i,j,k,iunit,npts
      integer idate(8)
      
780   format (A20,i4.4,"-",i2.2,"-",i2.2," ",i2.2,":",i2.2,":00 UTC")
781   format (2x,A)
782   format (2x,F12.4,x,F12.4)
783   format (A19,F5.2,'.',I2.2)
784   format (A20,i4.4,"-",i2.2,"-",i2.2," ",i2.2,":",i2.2)

      call date_and_time(values=idate)

      write(iunit,'(A20,A)')    "SAVE_POINT_OUTPUT   ", trim(savept(i)%names(j))
      write(iunit,780)          "REFERENCE_TIME      ", iyr,imo,iday,ihr,imin
      write(iunit,'(A20,A)')    "SAVE_POINT_LABEL    ", '"'//trim(splabel)//'"'
      write(iunit,784)          "CREATION_DATE       ", idate(1),idate(2),idate(3),idate(5),idate(6)
      write(iunit,783)          "CMS_VERSION         ", version,revision
      
      write(iunit,'(A20,A)')    "HORIZ_PROJECTION    ", trim(HProj)
      write(iunit,'(A20,A)')    "VERT_DATUM          ", trim(VProj)

      write(iunit,'(A,F0.2,A)') 'OUTPUT_INTERVAL     ', savept(i)%time_inc, ' sec'
      write(iunit,'(A)')        'TIME_UNITS          hrs'    !always HOURS, just in for clarity to user
      write(iunit,'(A20,A)')    "OUTPUT_UNITS        ", savept(i)%ounit(j)
      write(iunit,'(A20,I0)')   "NUMBER_POINTS       ", npts
      write(iunit,'(A)')        " "
      write(iunit,'(A)')        "NAME_BEGIN"  
      do k=1,npts    
        write(iunit,781) TRIM(savept(i)%id(k))
      enddo
      write(iunit,'(A)')     "NAME_END"      
      write(iunit,'(A)')     " "
      write(iunit,'(A)')     "XY_BEGIN"
      do k=1,npts
        write(iunit,782) savept(i)%x(k), savept(i)%y(k)
      enddo
      write(iunit,'(A)')     "XY_END"
      write(iunit,'(A)')     " "
      if (savept(i)%vals(j)==1) then
        write(iunit,'(A,T21,A,I0,A)')  "EXAMPLE_TS_OUTPUT:","'<time>  [<value(i)> ,i=1,",npts,"]'"
        write(iunit,'(A)') "SCALAR_TS_BEGIN"
      else
        write(iunit,'(A,T21,A,I0,A)')  "EXAMPLE_TS_OUTPUT:","'<time>  [<u_value(i)> <v_value(i)> ,i=1,",npts,"]'"
        write(iunit,'(A)') "VECTOR_TS_BEGIN    !Global space"
      endif
      !after this comes the actual data in the same order as POINT_NAMES above.

      return
      endsubroutine write_savept_header
    endsubroutine save_point_init
    
!***********************************************************************
    subroutine save_point_close
!***********************************************************************
    use out_def
    implicit none
    integer :: i,j,iunit
    
    do i=1,ngroups
      if(savept(i)%active)then
        do j=1,savept(i)%nvar
          iunit = savept(i)%funits(j)
          open(iunit,file=savept(i)%files(j),access='append')
          if(savept(i)%vals(j)==1)then
            write(iunit,'(A)') "SCALAR_TS_END"
          elseif(savept(i)%vals(j)==2)then
            write(iunit,'(A)') "VECTOR_TS_END"
          endif
          close(iunit)
        enddo
      endif
    enddo
    
    return
    endsubroutine save_point_close
    
!***********************************************************************
    subroutine write_save_point
! Writes data to all observation cell files
! 5/24/2012 - modify to output vector information in Global Space
! written by Mitchell Brown, USACE-ERDC-CHL, Lasted modified May 8, 2012
! modified by Alex Sanchez, USACE-ERDC-CHL, last modified June 17, 2013
!***********************************************************************    
#include "CMS_cpp.h"
    use size_def, only: ncellsD,ncells
    use geo_def, only: azimuth_fl,zb,zb0
    use hot_def, only: hstarttime,coldstart
    use out_def, only: savept,ngroups
    use flow_def, only: u,v,eta,h,iwet,us,vs
    use sed_def, only: sedtrans,nsed,Ct,rs,qtx,qty,pbk,epsvk,bsk
    use sal_def, only: saltrans,sal
    use heat_def, only: heattrans,heat
    use comvarbl, only: ctime,timehrs,timesecs
    use const_def, only: eps, deg2rad
    use wave_flowgrid_def, only: wunity, wunitx,whgt,wper,wavediss
    use prec_def
    use cms_def, only: cmswave
    implicit none
    integer :: i
    real(ikind) :: Csus(ncellsD),qx(ncellsD),qy(ncellsD),theta
    real(ikind) :: vecx(ncellsD),vecy(ncellsD),Ug(ncellsD),Vg(ncellsD),Xg(ncellsD),Yg(ncellsD)
    
    if(coldstart .or. ((hstarttime /= timehrs) .and. .not.coldstart))then  !otherwise there could be duplicates written to this file.
      theta = azimuth_fl*deg2rad  !For local -> global conversion
      do i=1,ngroups
        if(savept(i)%active .and. abs(mod(timesecs,savept(i)%time_inc))<0.5)then
          selectcase(i)
          case(1)    !Hydro
            call write_savept_scalar(1,1,eta)
            Ug = (u-us)*cos(theta)  - (v-vs)*sin(theta) !Convert to global  
            Vg = (u-us)*sin(theta)  + (v-vs)*cos(theta) !Convert to global  
            call write_savept_vector(1,2,Ug,Vg)
            qx = Ug*h  !Water current volume fluxes
            qy = Vg*h  !Water current volume fluxes
            call write_savept_vector(1,3,qx,qy)
            call write_savept_scalar(1,4,h) !Total water depth
          case(2)    !Sediment
            if (sedtrans) then
              Csus = iwet*Ct*rs  !Suspended sediment concentration
              Xg = qtx*cos(theta) - qty*sin(theta)
              Yg = qtx*sin(theta) + qty*cos(theta)     !Convert to global
              call write_savept_scalar(2,1,Csus)
              call write_savept_vector(2,2,Xg,Yg)
#ifdef DEV_MODE              
              do ii=1,ncells
                Xg(ii) = iwet(ii)*sum(pbk(ii,:,1)*epsvk(ii,:))  !Vertical mixing coefficient
                Yg(ii) = iwet(ii)*sum(pbk(ii,:,1)*bsk(ii,:))    !Suspended load correction factor
              enddo
#endif
              call write_savept_scalar(2,3,Xg)
              call write_savept_scalar(2,4,Yg)
            endif
          case (3)    !Salinity
            if (saltrans) then
              call write_savept_scalar(3,1,sal)
            endif
          case (4)    !Waves
            if (cmswave) then
              !Wave height vector
              vecx=whgt*wunitx
              vecy=whgt*wunity
              !Convert wave vector to global
              Xg = vecx*cos(theta) - vecy*sin(theta)
              Yg = vecx*sin(theta) + vecy*cos(theta) 
              call write_savept_scalar(4,1,Whgt)
              call write_savept_scalar(4,2,Wper)
              call write_savept_vector(4,3,Xg,Yg)
              call write_savept_scalar(4,4,wavediss)
            endif
          case (5) !Morphology
            vecx = -zb      !Depth
            call write_savept_scalar(5,1,vecx)
            if (sedtrans) then 
              vecy = zb - zb0 !Morphology change
              call write_savept_scalar(5,2,vecy)
            endif  
          case default
            continue
          end select
        endif
      enddo
    endif
      
    return
    endsubroutine write_save_point
    
!********************************************************************************   
    subroutine write_savept_scalar(i,j,var)
!> Writes a scalar savept cell variable to its corresponding file
!> written by Mitchell Brown, USACE-ERDC-CHL    
!********************************************************************************   
! Presently the maximum number of output values is 1000.  This value is checked upon cards being read.  MEB  10/22/2019
! This value is contained in the global output variable MAXOUT    
    
    use size_def, only: ncellsD
    use out_def, only: savept
    use comvarbl, only: timehrs
    use prec_def
    
    implicit none
    integer, intent(in)    :: i               !< Corresponds to the group number 
    integer, intent(in)    :: j               !< Corresponds to the member of the group
    integer                :: k
    real(ikind),intent(in) :: var(ncellsD)    !< scalar array of size 'ncellsD'
    
    open(savept(i)%funits(j),file=savept(i)%files(j),access='append')
973 format(1x,1pe11.4,2x,1000(1x,1pe11.4))
    write(savept(i)%funits(j),973) timehrs,(var(savept(i)%cell(k)),k=1,savept(i)%ncells)  
    close(savept(i)%funits(j)) !closing the file is useful for viewing the results during the simulation
    
    return
    endsubroutine write_savept_scalar

!********************************************************************************   
    subroutine write_savept_vector(i,j,var1,var2)
!> Writes a vector savept cell variable to its corresponding file
!> written by Mitchell Brown, USACE-ERDC-CHL    
!********************************************************************************   
! Presently the maximum number of output values is 1000.  This value is checked upon cards being read.  MEB  10/22/2019
! This value is contained in the global output variable MAXOUT    

    use size_def, only: ncellsD
    use out_def, only: savept
    use comvarbl, only: timehrs
    use prec_def
    
    implicit none
    integer, intent(in)    :: i               !< Corresponds to the group number 
    integer, intent(in)    :: j               !< Corresponds to the member of the group
    integer :: k
    real(ikind),intent(in) :: var1(ncellsD)   !< vector component of size 'ncellsD'
    real(ikind),intent(in) :: var2(ncellsD)   !< vector component of size 'ncellsD'
    
    open(savept(i)%funits(j),file=savept(i)%files(j),access='append')
973 format(1x,1pe11.4,2x,1000(1x,1pe11.4))
    write(savept(i)%funits(j),973) timehrs,  &
        (var1(savept(i)%cell(k)),var2(savept(i)%cell(k)),k=1,savept(i)%ncells)  
    close(savept(i)%funits(j)) !closing the file is useful for viewing the results during the simulation
    
    return
    endsubroutine write_savept_vector

!!********************************************************************************   
!    subroutine write_savept_velpro(i,j)
!!> Writes a scalar savept cell variable to its corresponding file
!!> written by Mitchell Brown, USACE-ERDC-CHL    
!!********************************************************************************   
!    use size_def, only: ncellsD
!    use out_def, only: savept
!    use comvarbl, only: timehrs
!    use prec_def
!    implicit none
!    integer, intent(in)    :: i               !< Corresponds to the group number 
!    integer, intent(in)    :: j               !< Corresponds to the member of the group
!    integer                :: k
!
!    open(savept(i)%funits(j),file=savept(i)%files(j),access='append')
!973 format(1x,F10.4,255(1X,1pe11.4))
!    do k=1,savept(i)%ncells  
!      i=savept(i)%cell(k)
!      call q3d_velpro(i,nzpsta,zpsta,iprofile,uc,vc,vx,vy,wx,wy,zap,zwp,&
!        Kc,visvm,visv0,visvslp,uzpsta,vzpsta) !Compute velocity vertical profile
!      write(savept(i)%funits(j),973) timehrs,(var(savept(i)%cell(k)),k=1,savept(i)%ncells)  
!    enddo
!    close(savept(i)%funits(j)) !closing the file is useful for viewing the results during the simulation
    !
    !return
    !endsubroutine write_savept_velpro

    
    
!Placing these NetCDF routines here for now so I don't have to recompile the 'out_lib' and all its dependencies.  MEB    
    
!!**************************************************************************
!    subroutine writescalnc(afile,apath,aname,var,aunits,timehr,iwritedry)
!! writes a scalar dataset to the NetCDF file with id number PID
!!
!! written by Mitchell Brown, USACE-ERDC-CHL   03/05/2020
!!**************************************************************************
!    use size_def, only: ncellsD,ncellsfull,ncellpoly
!    use netcdf 
!    use interp_lib, only: interp_scal_cell2node
!    use prec_def
!    implicit none
!    !Input/Output
!    character(len=*),intent(in) :: afile,apath,aname,aunits
!    real(ikind),intent(in) :: var(ncellsD),timehr
!    integer,intent(in) :: iwritedry
!    !Internal variables
!    integer :: fid,gid,ierr    
!    real(8) :: timed  !Must be double    
!    real(4) :: scalout(ncellsfull) !Must be single
!    character(len=200) :: afullpath
!    
!    if(ncellpoly>0)then
!      call interp_scal_cell2node(var,scalout,iwritedry)
!    else
!      call map_scal_active2full(var,scalout,iwritedry)  
!    endif
!    
!    afullpath=trim(apath)//trim(aname)    
!    ierr = nf90_open(trim(afile),ND90_WRITE,fid)            !Open existing NC file if it exists
!    if(ierr /= NF90_NOERR)then
!      ierr = nf90_create(trim(afile),NF90_CLOBBER,fid)       !Create new file if it didn't exist
!    endif                 
!    call OPEN_CREATE_NCDATASET(fid,trim(afullpath),gid,1,aunits,ierr)  !Open/create dataset          
!    timed = dble(timehr)
!    call XF_WRITE_SCALAR_TIMESTEP(gid,timed,ncellsfull,scalout,ierr) !Write data to XMDF file
!    call XF_CLOSE_GROUP(gid,ierr)  !Close dataset    
!    call XF_CLOSE_FILE(fid,ierr)   !Close XMDF file
!    
!    return
!    endsubroutine writescalnc
!
!!**************************************************************************
!    SUBROUTINE OPEN_CREATE_NCDATASET(OCID,STRING1,OCDID,DIM,OUNITS,OCERR)
!! THIS FUNCTION OPENS OR CREATES A SCALAR OR VECTOR DATASET AND 
!!   EXITS IF IT CAN'T CREATE ONE.
!!
!! - OCID:     parent group/file id
!! - STRING1:  name for dataset
!! - OCDID:    child group id
!! - DIM:      dataset dimensions (1- scalar, 2- vector)
!! - OUNITS:   Data units ('m','m/s','kg/m^3',etc.)
!! - REFTIME (in):  SMS Reference time for dataset
!!
!! written by MEB 03/05/20
!!  - fixed time units to hours 
!!  - added compression option
!!  - added output units for dataset
!!**************************************************************************
!    use comvarbl, only: reftime
!    use out_def, only: ixmdfcomp
!    use diag_lib
!    use netcdf
!    implicit none
!    integer OCERR,OCDID,DIM,OCDPID,OCID
!    character(LEN=*) STRING1,OUNITS    
!      
!    SELECT CASE (DIM)
!      CASE (1)
!        OCERR = nf90_inq_ncid (OCID, STRING1, OCDID) 
!        if (OCERR /= nf90_noerr) then
!          OCERR = nf90_def_grp (OCID, STRING1, OCDID)
!          if (OCERR /= nf90_noerr) then
!            call diag_print_error('Could not create scalar NC dataset: ',STRING1)
!          endif
!        endif
!
!!        CALL XF_OPEN_GROUP(OCID,STRING1,OCDID,OCERR)
!!        IF(OCERR<0)THEN
!!          !CALL XF_CREATE_SCALAR_DATASET(OCID,STRING1,'none',OUNITS,0,OCDID,OCERR)
!!          !Note: Time units always in hours (hard-wired from SMS)
!!          CALL XF_CREATE_SCALAR_DATASET(OCID,STRING1,OUNITS,TS_HOURS,ixmdfcomp,OCDID,OCERR)
!!          IF(OCERR<=0)THEN
!!            call diag_print_error('Could not create scalar dataset: ',STRING1)            
!!          ENDIF
!!          CALL XF_DATASET_REFTIME(OCDID,REFTIME,OCERR)
!!          CALL XF_CREATE_PROPERTY_GROUP(OCDID,OCDPID,OCERR)
!!          CALL XF_WRITE_PROPERTY_FLOAT(OCDPID,PROP_NULL_VALUE,1,-999.0,NONE,OCERR)
!!          CALL XF_SCALAR_DATA_LOCATION(OCDID,GRID_LOC_CENTER,OCERR)
!!        ENDIF
!      
!      CASE (2)
!        OCERR = nf90_inq_ncid (OCID, STRING1, OCDID) 
!        if (OCERR /= nf90_noerr) then
!          OCERR = nf90_def_grp (OCID, STRING1, OCDID)
!          if (OCERR /= nf90_noerr) then
!            call diag_print_error('Could not create vector NC dataset: ',STRING1)
!          endif
!        endif
!
!!        CALL XF_OPEN_GROUP(OCID,STRING1,OCDID,OCERR)
!!        IF(OCERR<0)THEN
!!          !CALL XF_CREATE_VECTOR_DATASET(OCID,STRING1,'none',OUNITS,0,OCDID,OCERR)
!!          !Note: Time units always in hours (hard-wired from SMS)
!!          CALL XF_CREATE_VECTOR_DATASET(OCID,STRING1,OUNITS,TS_HOURS,ixmdfcomp,OCDID,OCERR) 
!!          IF(OCERR<=0)THEN
!!            call diag_print_error('Could not create vector dataset: ',STRING1)   
!!          ENDIF
!!          CALL XF_DATASET_REFTIME(OCDID,REFTIME,OCERR)
!!          CALL XF_CREATE_PROPERTY_GROUP(OCDID,OCDPID,OCERR)
!!          CALL XF_WRITE_PROPERTY_FLOAT(OCDPID,PROP_NULL_VALUE,1,-999.0,NONE,OCERR)
!!          CALL XF_VECTORS_IN_LOCAL_COORDS(OCDID,OCERR)
!!!          CALL XF_VECTOR_2D_DATA_LOCS (OCDID,GRID_LOC_FACE_I,GRID_LOC_FACE_J,OCERR)
!!          CALL XF_VECTOR_2D_DATA_LOCS(OCDID,GRID_LOC_CENTER,GRID_LOC_CENTER,OCERR)
!!        ENDIF
!      
!      END SELECT
!      
!      RETURN
!    END SUBROUTINE OPEN_CREATE_NCDATASET

    