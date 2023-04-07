!================================================================================
! CMS Flow Boundary routines
!
! Contains:
!   bnd_default - Sets the default parameters for the flow boundary variables
!   bnd_cards - Reads the flow boundary cards from the card file
!   bnd_init - Initializes the flow boundary variables
!   bnd_print - Prints the flow boundary settings to the diagnostic file
!              and the screen
!   bnd_eval - Evaluates the boundary conditions
!   bound_uv - Applies the boundary conditions to the momentum equations
!   bound_pp - Applies the boundary conditions to the pressure correction equation
!   bnd_wse_adjust - Adjusts the water level at a boundary cell string to 
!                   include the wind and wave setup
!   bndflux - Updates the volume fluxes at all the boundaries
!   bndriver_vel - Updates the current velocity at flux boundaries
!   flux_alloc - Resizes the flux boundary condition variable 
!   flux_block - Reads a flux boundary block from the card file
!   wse_block - Reads an wse boundary condition block from the card file
!   vel_block - Reads a velocity block from the card file
!   singlewse_alloc - Resizes the single water level boundary condition variable 
!   harmonic_block - Reads a harmonic constituent block from the card file
!   tidal_block - Reads a tidal constituent block from the card file
!   parent_block - Reads a parent simulation block from the card file
!   tdb_block - Reads a tidal database block from the card file
!   multiwse_alloc - Resizes the multiple water level boundary condition variables   
!   multiwsevel_alloc - Resizes the multiple water level and velocity boundary
!                      condition variable
!   xshore_uv - Solves the cross-shore momentum equations for the current velocities
!   xshore_wse - Solves the cross-shore momentum equation for water levels
!   bnd_pp - Pressure correction boundary condition
!   bnd_vel_flux - Calculatest the velocity at a flux boundary condition
!   bnd_vel_openwse - Calculates the velocity at an open wse boundary condition
!   bnd_vel_openwsevel - Calculates the velocity at an open wse and velocity 
!                       boundary condition
!
! written by Alex Sanchez, USACE-CHL
!================================================================================

!***************************************************************************   
    subroutine bnd_default
! Sets the default parameters for the boundary variables
! written by Alex Sanchez, USACE-CHL; Weiming Wu, NCCHE
!***************************************************************************
    use bnd_def
    
    implicit none
    
    !Parameters
    veldamp = 0.0 !Damping factor for velocity at open wse boundaries
    
    !Boundary strings (contains cell/nodestring information)
    nbndstr = 0     
    nbcstr  = 0 !Temporary used for storing 2dm nodestrings
    
    !Boundaries (contain boundary information)
    nTHstr  = 0     !Tidal/harmonic BC, either 0 or 1
    nHstr   = 0     !Single WSE BC
    nQstr   = 0     !Flux BC
    nMHstr  = 0     !Multiple WSE BC
    nMHVstr = 0     !Multiple WSE and Velocity BC
    nCSstr  = 0     !Cross-shore boundaries
    nNHstr  = 0     !Nested WSE BC
    nNHVstr = 0     !Nested WSE and Velocity BC
    nNTHstr  = 0    !Nested Tidal Dabase WSE BC
    nNTHVstr  = 0   !Nested Tidal Dabase WSE and Velocity BC
    
    !Tidal constituent boundary
    Tread = .false.     !Tidal/harmonic cellstring
    tide_read = .false. !Tidal Constituents block
    harm_read = .false. !Harmonic Components block
    
    !Parent Simulations
    nParSim = 0
    
    return
    end subroutine bnd_default

!***************************************************************************   
    subroutine bnd_cards(cardname,foundcard)
! Reads the flow boundary cards from the card file
! written by Alex Sanchez, USACE-CHL; Weiming Wu, NCCHE
!***************************************************************************
    use bnd_def
    use flow_def, only: grav,viscos
    use const_def, only: pi
    use comvarbl, only: flowpath
    !use diag_def
    use prec_def, only: ikind
    
    implicit none    
    integer :: ierr  !Wu
    character(len=37) :: cardname
    logical :: foundcard    
    
    foundcard = .true.
    select case (cardname)
    case('BOUNDARY_VELOCITY_DAMPING','VELOCITY_DAMPING')
      backspace(77)  
      read(77,*) cardname,veldamp
      veldamp = min(max(veldamp,0.0),100.0)
      
    case('BOUNDARY_BEGIN')
      call bnd_block  
    
    case('QDRIVER_CELLSTRING')        !Flux BC
      !Note: when using this card, the id and data files are assumed to be same
      !To allow them to be different a block structure input format must be used
      call flux_alloc
      
      call card_bid(flowpath,Q_str(nQstr)%bidfile,Q_str(nQstr)%bidpath,Q_str(nQstr)%idnum)
      !Note: when this card is used, it is assumed that the id file is the same as the data file
      Q_str(nQstr)%fluxfile = Q_str(nQstr)%bidfile
      Q_str(nQstr)%fluxpath = Q_str(nQstr)%bidpath
      Q_str(nQstr)%istrtype = 1
      
    case('TIDAL_CELLSTRING')          !Tidal BC
      if(.not.tide_read) call tidal_alloc
      call card_bid(flowpath,TH_str(nTHstr)%bidfile,TH_str(nTHstr)%bidpath,TH_str(nTHstr)%idnum)
      Tread = .true.
      if(tide_read) Tread = .false.; tide_read = .false. !Both read in, prepare for next string
      !TH_str(nTHstr)%istidal = .true. !Default already true and this may overide user-input, Alex: 08-22-14
      TH_str(nTHstr)%istrtype = 1
     
    case('TIDAL_CONSTITUENTS_BEGIN')     
      tide_read = .true.   
      call read_tidal_const
      if(Tread)then !Both read in, prepare for next string
        Tread = .false.
        tide_read = .false.
     endif 
      
    case('TIDAL_CONSTITUENTS_END','END')
      
    case('HARMONIC_CELLSTRING')          !Harmonic BC
      !if(tide_read)then
      !  call diag_print_warning('Cannot specify both tidal and harmonic boundary condition')
      !endif
      if(.not.tide_read) call tidal_alloc
      call card_bid(flowpath,TH_str(nTHstr)%bidfile,TH_str(nTHstr)%bidpath,TH_str(nTHstr)%idnum)
      Tread = .true.          
      if(tide_read) Tread = .false.; tide_read = .false. !Both read in, prepare for next string
      TH_str(nTHstr)%istidal = .false.
      TH_str(nTHstr)%istrtype = 1
      
    case('HARMONIC_BEGIN','HARMONICS_BEGIN','HARMONIC_COMPONENTS_BEGIN')     
      tide_read = .true.   
      call read_harmonics
      if(Tread)then !Both read in, prepare for next string
        Tread = .false.
        tide_read = .false.
     endif   
      
    case('TIDAL_CORRECTIONS')
      call card_boolean(77,TH_str(nTHstr)%istidal,ierr)
      
    case('HDRIVER_CELLSTRING')        !Single Water Level BC        
      call singlewse_alloc
      call card_bid(flowpath,H_str(nHstr)%bidfile,H_str(nHstr)%bidpath,H_str(nHstr)%idnum)
      !Note: when this card is used, it is assumed that the id file is the same as the data file
      H_str(nHstr)%wsefile = H_str(nHstr)%bidfile
      H_str(nHstr)%wsepath = H_str(nHstr)%bidpath 
      H_str(nHstr)%istrtype = 1
        
    case('MULTI_HDRIVER_CELLSTRING')  !Multiple Water Level BC (Type 4)
      call multiwse_alloc
      call card_bid(flowpath,MH_str(nMHstr)%bidfile,MH_str(nMHstr)%bidpath,MH_str(nMHstr)%idnum)
      !Note: when this card is used, it is assumed that the id file is the same as the data file
      MH_str(nMHstr)%wsefile = MH_str(nMHstr)%bidfile
      MH_str(nMHstr)%wsepath = MH_str(nMHstr)%bidpath
      MH_str(nMHstr)%istrtype = 1
        
    case('MULTI_VDRIVER_CELLSTRING','MULTI_HVDRIVER_CELLSTRING') !Multiple Water Level and Velocity BC (Type 5)
      call multiwsevel_alloc
      call card_bid(flowpath,MHV_str(nMHVstr)%bidfile,MHV_str(nMHVstr)%bidpath,MHV_str(nMHVstr)%idnum)
      !Note: when this card is used, it is assumed that the id file is the same as the data file
      MHV_str(nMHVstr)%wsefile = MHV_str(nMHVstr)%bidfile
      MHV_str(nMHVstr)%wsepath = MHV_str(nMHVstr)%bidpath
      MHV_str(nMHVstr)%velfile = MHV_str(nMHVstr)%bidfile
      MHV_str(nMHVstr)%velpath = MHV_str(nMHVstr)%bidpath
      MHV_str(nMHVstr)%istrtype = 1
        
    case('CROSS-SHORE_CELLSTRING','CROSS_SHORE_CELLSTRING')    !Cross-shore BC (Type 6)
      call xshore_alloc
      call card_bid(flowpath,CS_str(nCSstr)%bidfile,CS_str(nCSstr)%bidpath,CS_str(nCSstr)%idnum)
      CS_str(nCSstr)%istrtype = 1
        
    case default
      foundcard = .false.
        
    end select

    return
    end subroutine bnd_cards

!****************************************************************    
    subroutine xshore_alloc
!****************************************************************
    use bnd_def, only: nCSstr,CS_str,CS_type
    
    implicit none
    integer :: i
    type(CS_type),allocatable :: CS_temp(:)
    
    nCSstr = nCSstr + 1
    if(nCSstr==1)then
      allocate(CS_str(nCSstr))
    else
      allocate(CS_temp(nCSstr-1))
      do i=1,nCSstr-1
        CS_temp(i) = CS_str(i)
      enddo
      deallocate(CS_str)
      allocate(CS_str(nCSstr))
      do i=1,nCSstr-1        
        CS_str(i) = CS_temp(i)  
      enddo
      deallocate(CS_temp)
    endif
      
    return
    end subroutine xshore_alloc
      
!**********************************************************    
    subroutine flux_alloc
! Resizes the flux boundary condition variable    
!**********************************************************
    use bnd_def, only: nQstr,Q_type,Q_str
    !use prec_def
    
    implicit none
    integer :: i
    type(Q_type),allocatable :: Q_temp(:)
    
    nQstr=nQstr+1
    if(nQstr==1)then
      allocate(Q_str(nQstr))
    else
      allocate(Q_temp(nQstr-1))
      do i=1,nQstr-1
        Q_temp(i)= Q_str(i)
      enddo
      deallocate(Q_str)
      allocate(Q_str(nQstr))      
      do i=1,nQstr-1
        Q_str(i) = Q_temp(i)
      enddo
      deallocate(Q_temp)
    endif
    
    !Initialize and set default values
    Q_str(nQstr)%idnum = 0    
    Q_str(nQstr)%ifluxmode = 2  !1-Constant flux, 2-Total-flux curve
    Q_str(nQstr)%bidfile = ''
    Q_str(nQstr)%bidpath = ''
    Q_str(nQstr)%fluxfile = ''
    Q_str(nQstr)%fluxpath = ''    
    Q_str(nQstr)%cmvel = 0.667  !Default value
    Q_str(nQstr)%angle = -999.0 !No default
    Q_str(nQstr)%qflux = 0.0
    Q_str(nQstr)%qfluxconst = 0.0
    Q_str(nQstr)%ifluxunits = 0  !Input flux units: 0-m^3/s/cell,1-m^3/s,2-ft^3/s
    Q_str(nQstr)%ncells = 0
    Q_str(nQstr)%ntimes = 0
    Q_str(nQstr)%inc = 1
    Q_str(nQstr)%nti = 1 !Order of temporal interpolation (Use first order for fluxes because they can be noisy)
    Q_str(nQstr)%nsw = 0 !Temporal smoothing width
    Q_str(nQstr)%nsi = 0 !Temporal smoothing iterations
    Q_str(nQstr)%nstages = 0 !Stage-Flow curve
    
    return
    end subroutine flux_alloc
    
!**************************************************************************************    
    subroutine flux_block(ibndtype,ifluxmode,fluxfile,fluxpath,qfluxconst,ifluxunits,&
                  angle_flux,cmvel,nti,nsi,nsw)
! Reads an a flux boundary condition block from the card file
!
! written by Alex Sanchez, USACE-CHL
!**************************************************************************************
    use geo_def, only: azimuth_fl
    use const_def, only: deg2rad
    use comvarbl, only: mpfile,flowpath
    use diag_lib, only: diag_print_error
    use prec_def, only: ikind
    
    implicit none
    !Input/Output
    character(len=*),intent(inout) :: fluxfile,fluxpath
    integer,intent(out) :: ibndtype,ifluxmode,ifluxunits,nti,nsi,nsw
    real(ikind),intent(out) :: qfluxconst,angle_flux,cmvel
    !Internal Variables
    integer :: kk,ierr
    character :: cardname*37,qunits*40,aline*200
    logical :: foundcard
    
    ibndtype = 0 !Q-1
    
    !Default units
    qunits = 'm^3/s/cell'
    
    do !  Removed the limit of only reading 10 cards in this block.  !kk=1,10  MEB 06/27/2018
      foundcard = .true.  
      cardname = ''
      read(77,*,iostat=ierr) cardname
      if(ierr/=0) exit
      if(cardname(1:1)=='!' .or. cardname(1:1)=='#' .or. cardname(1:1)=='*') cycle
      select case(cardname)
      case('RATING_CURVE','STAGE_FLOW_CURVE','STAGE_FLUX_CURVE')
        call card_dataset(77,mpfile,flowpath,fluxfile,fluxpath,1)
        ibndtype = 1 !Q-1
        ifluxmode = 3
        
      case('FLUX_CONSTANT','FLUX_VALUE','VALUE')  
        backspace(77)
        read(77,'(A)') aline
        read(aline,*,iostat=ierr) cardname,qfluxconst,qunits
        ibndtype = 1 !Q-1
        ifluxmode = 1
        
      case('FLUX_CURVE','FLUX_DATA','FLUX_TIMESERIES','CURVE')
        call card_dataset(77,mpfile,flowpath,fluxfile,fluxpath,1)
        ibndtype = 1 !Q-1
        ifluxmode = 2
        
      case('FLUX_UNITS','UNITS')     !For SMS 11.2 and previous - value is per cell, not total and so the FLUX_UNITS will be written to file from SMS 12.x and after.
        backspace(77)
        read(77,'(A)') aline
        read(aline,*,iostat=ierr) cardname,qunits
          
      case('INFLOW_DIRECTION','INFLOW_ANGLE','ANGLE','DIRECTION','POSITIVE_FLOW_DIRECTION')
        backspace(77)
        read(77,*) cardname,angle_flux              !degrees
        angle_flux = 90.0 - angle_flux - azimuth_fl !degrees, clockwise from North, local coordinate system
        if(angle_flux<0.0)then
          angle_flux = angle_flux + 360.0
        elseif(angle_flux>360.0)then
          angle_flux = angle_flux - 360.0
        endif 
        
      case('CONVEYANCE_COEFFICIENT','CONVEYANCE_COEF','CONVEYANCE')
        backspace(77)
        read(77,*) cardname,cmvel
        cmvel = max(min(cmvel,1.0),0.5) !Limit values
          
      case('TEMPORAL_INTERPOLATION_ORDER','ORDER','ORDER_INTERP')
        backspace(77)
        read(77,*) cardname,nti
        nti = min(max(nti,1),3) !Limit between 1 and 3
      
      case('TIME_SMOOTH_ITER')
        backspace(77)
        read(77,*) cardname,nsi
        nsi = min(max(nsi,0),3) !Limit between 0 and 10
        
      case('TIME_SMOOTH_WIDTH')
        backspace(77)
        read(77,*) cardname,nsw  
        nsw = min(max(nsw,0),9) !Limit between 0 and 9
        
      case('FLUX_BOUNDARY_END','FLUX_END','END')
        exit
          
      case default
        foundcard = .false.
          
      end select
    enddo
    
    if(ibndtype==0)then
      call diag_print_error('Input flux data must be specified for flux block')
    endif
    
    select case(qunits) !0-m^3/s/cell,1-m^3/s,2-ft^3/s
    case('m^3/s/cell','m^3/sec/cell','CUBIC_METERS_PER_SECOND_PER_CELL','!','#')     !Default
      ifluxunits = 0
    case('m^3/s','m^3/s/boundary','m^3/sec','m^3/sec/boundary','CUBIC_METERS_PER_SECOND_PER_BOUNDARY')
      ifluxunits = 1
    case('ft^3/s','ft^3/sec','ft^3/s/boundary','ft^3/sec/boundary','CUBIC_FEET_PER_SECOND_PER_BOUNDARY')
      ifluxunits = 2
    case('ft^3/s/cell','ft^3/sec/cell','CUBIC_FEET_PER_SECOND_PER_CELL')
      ifluxunits = 3  
    case default
      call diag_print_error('Unrecognized Flux Units',&
        ' Optional units are:',&
        '  CUBIC_METERS_PER_SECOND_PER_CELL',&
        '  CUBIC_METERS_PER_SECOND_PER_BOUNDARY',&
        '  CUBIC_FEET_PER_SECOND_PER_BOUNDARY')
    end select
    
    return
    end subroutine flux_block

!**************************************************************************    
    subroutine wse_block(ibndtype,istidal,wsefile,wsepath,wseconst,&
     ioffsetmode,offsetfile,offsetpath,wseoffset,&        !(hli,01/18/17)
     dwsex,dwsey,minterp,nti,nsi,nsw,nssi,nssw,wseout,wseadjust)
! Reads an wse boundary condition block from the card file
! written by Alex Sanchez, USACE-CHL
!**************************************************************************
    use geo_def, only: azimuth_fl
    use const_def, only: deg2rad
    use comvarbl, only: mpfile,flowpath
    use diag_lib, only: diag_print_error
    use prec_def, only: ikind
    
    implicit none
    !Input/Output
    character(len=*),intent(inout) :: wsefile,wsepath,offsetfile,offsetpath !(hli,01/18/17)
    integer,         intent(inout) :: ibndtype,minterp,nti,nsi,nsw,nssi,nssw
    integer,         intent(out)   :: ioffsetmode                           !(hli,01/18/17)
    real(ikind),     intent(inout) :: wseconst,wseoffset,dwsex,dwsey
    logical,         intent(inout) :: istidal,wseout,wseadjust
    !Internal Variables
    integer :: kk,ierr
    character(len=37) :: cardname
    character(len=40) :: cdum
    real(ikind) :: cosang,sinang,vecx,vecy
    logical :: foundcard
    
    ibndtype = 0  !Undefined
    
    do !  Removed the limit of only reading 20 cards in this block.  !kk=1,20  MEB 06/27/2018
      foundcard = .true. 
      cardname = ''
      read(77,*,iostat=ierr) cardname
      if(ierr/=0) exit
      if(cardname(1:1)=='!' .or. cardname(1:1)=='#' .or. cardname(1:1)=='*') cycle
      !call bnd_time_smooth_cards(cardname,nti,nsi,nsw,foundcard)
      !if(foundcard) cycle 
      !call bnd_space_smooth_cards(cardname,nssi,nssw,foundcard)
      !if(foundcard) cycle 
      select case(cardname)    
      case('WSE_BOUNDARY_END','WSE_END','WSE_FORCING_END','END')
        exit
          
      case('WSE_CONSTANT','WSE_VALUE','VALUE')
        call card_scalar(77,'m','m',wseconst,ierr)
        ibndtype = 3
        istidal = .false.
        nti = 0
        
      case('WSE_CURVE','CURVE')
        call card_dataset(77,mpfile,flowpath,wsefile,wsepath,1)
        ibndtype = 3
        istidal = .false.
        
      case('WSE_DATA','WSE_DATASET','DATASET','DATA','CURVES')      
        call card_dataset(77,mpfile,flowpath,wsefile,wsepath,1)
        ibndtype = 4
        istidal = .false.
        
      case('WSE_OUTPUT','OUT','OUTPUT') 
        call card_boolean(77,wseout,ierr)
        
      case('WSE_OUTPUT_FILE','OUTPUT_FILE') 
        backspace(77)
        read(77,*) cardname,wsefile
        wseout = .true.
        
      case('WSE_NEST','WSE_PARENT','SIMULATION','SOLUTION','WSE_SOLUTION')
        if(ibndtype/=8) ibndtype = 7
        istidal = .false.
        
      case('INTERPOLATION_METHOD')  
        backspace(77)
        read(77,*) cardname,cdum
        select case(cdum)
        case('PIECEWISE_POLYNOMIAL','PWP')    
          minterp = 1
        case('CUBIC_SPLINE','CS')
          minterp = 2
        case default
          minterp = 1
        end select
        
      case('TIME_INTERP_ORDER')
        backspace(77)
        read(77,*) cardname,nti
        nti = min(max(nti,1),3) !Limit between 1 and 3
        minterp = 1
        
      case('TIME_SMOOTH_ITER')
        backspace(77)
        read(77,*) cardname,nsi
        nsi = min(max(nsi,0),10) !Limit between 0 and 10
        
      case('TIME_SMOOTH_WIDTH')
        backspace(77)
        read(77,*) cardname,nsw  
        nsw = min(max(nsw,0),9) !Limit between 0 and 9
      
      case('SPACE_SMOOTH_ITER','SPACE_SMOOTH_ITERATIONS','SPACIAL_SMOOTHING_ITERATIONS',&
        'SMOOTHING_ITERATIONS','SMOOTH_ITERATIONS','SMOOTHING_ITER')
        backspace(77)
        read(77,*) cardname,nssi
        
      case('SPACE_SMOOTH_WIDTH','SPACAIAL_SMOOTHING_WIDTH','SMOOTHING_WIDTH','SMOOTH_WIDTH')
        backspace(77)
        read(77,*) cardname,nssw
        if(mod(nssw,2)/=0)then !Make sure it is odd
          nssw = nssw + 1
        endif
        
      case('WSE_OFFSET','OFFSET_WSE','OFFSET','WSE_OFFSET_CONSTANT','OFFSET_CONSTANT')
        call card_scalar(77,'m','m',wseoffset,ierr)  
        ioffsetmode=1
        
      case('WSE_OFFSET_CURVE','OFFSET_WSE_CURVE','OFFSET_CURVE')
        call card_dataset(77,mpfile,flowpath,offsetfile,offsetpath,1)   
        ioffsetmode=2
      
      case('WSE_ADJUSTMENT','ADJUSTMENT','WIND_WAVE_CORRECTION',&
         'WIND_WAVE_ADJUSTMENT','CORRECTION')  
        CALL card_boolean(77,wseadjust,ierr)
        
      case('WSE_GRADIENTS')
        if(ibndtype==4)then
          call diag_print_error('Cannot specify wse gradients for ',&
            '   nested boundary conditions')
        endif
        backspace(77)
        read(77,*) cardname,dwsex,dwsey
        !Rotate gradients to the local coordinate system
        cosang=cos(azimuth_fl*deg2rad)
        sinang=sin(azimuth_fl*deg2rad)
        vecx=dwsex
        vecy=dwsey
        dwsex= vecx*cosang+vecy*sinang
        dwsey=-vecx*sinang+vecy*cosang
      
      case('WSE_CONSTITUENTS')
        backspace(77)
        read(77,*) cardname,cdum
        select case(cdum)
        case('TIDAL_CONSTITUENTS','TIDAL')    
          ibndtype = 2  !2=TH 
          istidal = .true.
        case('HARMONIC_CONSTITUENTS','HARMONIC')
          ibndtype = 2  !2=TH 
          istidal = .false.
        case('TIDAL_DATABASE')
          if(ibndtype/=10) ibndtype = 9  !9=NTH 
        case default
          call diag_print_error('Invalid value for card WSE_CONSTITUENTS',&
            '  Valid options are:',&
            '    TIDAL_CONSTITUENTS',&
            '    HARMONIC_CONSTITUENTS',&
            '    TIDAL_DATABASE')
        end select
        
      case default
        foundcard = .false.
          
      end select
    enddo
    
    return
    end subroutine wse_block

!********************************************************************    
    subroutine vel_block(ibndtype,velfile,velpath,nti,nsi,nsw,&
      nssi,nssw,velout)
! Reads a velocity block from the card file
! written by Alex Sanchez, USACE-CHL
!********************************************************************
    use geo_def, only: azimuth_fl
    use const_def, only: deg2rad
    use comvarbl, only: mpfile,flowpath
    use diag_lib, only: diag_print_error
    !use prec_def
    
    implicit none
    !Input/Output
    integer,         intent(inout) :: ibndtype,nti,nsi,nsw,nssi,nssw
    character(len=*),intent(inout) :: velfile,velpath
    logical,         intent(inout) :: velout
    !Internal Variables
    integer :: kk,ierr
    character :: cardname*37,cdum*40
    logical :: foundcard
    
    ibndtype = 0  !Undefined
    
    do !  Removed the limit of only reading 10 cards in this block.  !kk=1,10  MEB 06/27/2018
      foundcard = .true.  
      read(77,*,iostat=ierr) cardname
      if(ierr/=0) exit
      if(cardname(1:1)=='!' .or. cardname(1:1)=='#' .or. cardname(1:1)=='*') cycle
      select case(cardname)        
      case('VEL_BOUNDARY_END','VEL_END','VEL_FORCING_END','END')
        exit
          
      case('VEL_DATA','VEL_DATASET','DATASET','DATA','CURVES')  
        call card_dataset(77,mpfile,flowpath,velfile,velpath,2)
        ibndtype = 5 !=MHV
        
      case('VEL_NEST','VEL_PARENT','SIMULATION','SOLUTION','VEL_SOLUTION')
        ibndtype = 8
        
      case('VEL_OUTPUT','OUT','OUTPUT')
        call card_boolean(77,velout,ierr)
        
      case('VEL_OUTPUT_FILE','OUTPUT_FILE') 
        backspace(77)
        read(77,*) cardname,velfile
        velout = .true.  
        
      case('TEMPORAL_INTERPOLATION_ORDER','ORDER','ORDER_INTERP')
        backspace(77)
        read(77,*) cardname,nti
        nti = min(max(nti,1),3) !Limit between 1 and 3
      
      case('TIME_SMOOTH_ITER')
        backspace(77)
        read(77,*) cardname,nsi
        nsi = min(max(nsi,0),10) !Limit between 0 and 10
        
      case('TIME_SMOOTH_WIDTH')
        backspace(77)
        read(77,*) cardname,nsw  
        nsw = min(max(nsw,0),9) !Limit between 0 and 9
      
      case('SPACE_SMOOTH_ITER','SPACE_SMOOTH_ITERATIONS','SPACIAL_SMOOTHING_ITERATIONS',&
        'SMOOTHING_ITERATIONS','SMOOTH_ITERATIONS','SMOOTHING_ITER')
        backspace(77)
        read(77,*) cardname,nssi
        
      case('SPACE_SMOOTH_WIDTH','SPACAIAL_SMOOTHING_WIDTH','SMOOTHING_WIDTH','SMOOTH_WIDTH')
        backspace(77)
        read(77,*) cardname,nssw
        if(mod(nssw,2)/=0)then !Make sure it is odd
          nssw = nssw + 1
        endif  
        
      case('VEL_CONSTITUENTS','CONSTITUENTS','VELOCITY_CONSTITUENTS')
        backspace(77)
        read(77,*) cardname,cdum
        select case(cdum)
        case('TIDAL_DATABASE')
          ibndtype = 9  !9=NTH 
          
        case default
          call diag_print_error('Invalid value for card VEL_CONSTITUENTS',&
            '  Only valid options is:','    TIDAL_DATABASE')
        end select
        
      case default
        foundcard = .false.
          
      end select
    enddo
    
    return
    end subroutine vel_block
     
!************************************************************************
    subroutine tidal_block(ibndtype,ntc,name,amp,phase,speed,f,vu,angle_wave,&
       ioffsetmode,offsetfile,offsetpath,wseoffset,nti)
! Reads a tidal constituent block from the card file
! Add "nti" for curve inperpolation
! written by Alex Sanchez, USACE-CHL
!************************************************************************
    use geo_def, only: azimuth_fl
    use const_def, only: deg2rad
    use comvarbl, only: mpfile,flowpath,iyr
    use bnd_def, only: ntf
    use tide_lib, only: tidal_data
    use unitconv_lib, only: unitconv_var
    use prec_def, only: ikind
    
    implicit none
    !Input/Output
    character(len=*),intent(inout) :: offsetfile,offsetpath !(hli,10/04/17)
    integer,         intent(out)   :: ioffsetmode           !(hli,10/04/17)
    real(ikind),     intent(inout) :: wseoffset
    
    integer,intent(inout) :: ibndtype
    integer,intent(out)  :: nti                   !(hli,10/06/17)
    integer              :: ntc         !Tidal constituents used
    real(ikind)          :: angle_wave  !Incident angle of tidal wave
    real(ikind), pointer :: amp(:)      !Amplitude [m] (constituent) 
    real(ikind), pointer :: speed(:)    !Speed [rad/hrs] (constituent)
    real(ikind), pointer :: phase(:)    !Phase [rad] (constituent)
    real(ikind), pointer :: f(:)        !Nodal factor [-] (constituent)
    real(ikind), pointer :: vu(:)       !Equilibrium argument [rad] (constituent)
    character(len=10), pointer :: name(:)     !Tidal Consitituent names (constituent)    
    !Internal Variables
    integer :: kk,ierr
    character(len=37) :: cardname
    character(len=6) :: ampunits,phaunits
    logical :: foundcard
    integer :: k,kth
    integer :: mcyctemp(ntf),ind(ntf)
    real(ikind) :: ampfac,phafac
    real(ikind) :: amptemp(ntf),phasetemp(ntf),speedtemp(ntf),ftemp(ntf),vutemp(ntf)
    character(len=10) :: nametemp(ntf)
    character(len=10) :: nameyr(ntf),astr  
    
!   Get constituent names, speeds, nodal corrections and equilibrium phases (relative to Greenwich)
    call tidal_data(iyr,nameyr,speedtemp,mcyctemp,ftemp,vutemp)

    ibndtype = 2 !=TH
    ampunits = 'm'
    phaunits = 'deg'
    ntc = 0 !Number of tidal constituents used
    do !  Removed the limit of only reading 10 cards in this block.  !kk=1,10  MEB 06/27/2018
      foundcard = .true.  
      read(77,*,iostat=ierr) cardname
      if(ierr/=0) exit
      if(cardname(1:1)=='!' .or. cardname(1:1)=='#' .or. cardname(1:1)=='*') cycle
      select case(cardname)
      case('TIDAL_CONSTITUENTS_END','TIDAL_END','END')
        exit
          
      case('AMPLITUDE_UNITS')
        backspace(77)
        read(77,*) cardname,ampunits
          
      case('PHASE_UNITS')
        backspace(77)
        read(77,*) cardname,phaunits
      
      case('INCIDENT_ANGLE','INCIDENT_DIRECTION')
        backspace(77)
        read(77,*) cardname,angle_wave  !degrees
        angle_wave = 90.0 - angle_wave - azimuth_fl !degrees, clockwise from North, local coordinate system
        if(angle_wave<0.0)then
          angle_wave = angle_wave + 360.0
        elseif(angle_wave>360.0)then
          angle_wave = angle_wave - 360.0
        endif
        angle_wave = angle_wave*deg2rad !Convert to rad
        
      case('CONSTITUENT')
        backspace(77)  
        read(77,*) cardname, astr  
        kth = 0
        do k=1,ntf
          if(astr==nameyr(k))then
            kth = k
            exit
          endif  
        enddo
        if(kth/=0)then
          ntc = ntc + 1
          ind(ntc) = kth
          backspace(77)   
          read(77,*) cardname,nametemp(kth),amptemp(kth),phasetemp(kth)
          nametemp(kth) = astr
        else
          foundcard = .false.
        endif
        ibndtype = 2
        
      case('WSE_OFFSET_CURVE','OFFSET_WSE_CURVE','OFFSET_CURVE')     !hli(10/04/17)
        call card_dataset(77,mpfile,flowpath,offsetfile,offsetpath,1)   
        ioffsetmode=2
        nti=2
 !     write(3000,*)'mpfile= ',mpfile,'flowpath = ',flowpath           !hli(10/04/17)
 !     write(3000,*)'offsetfile= ',offsetfile,'offsetpath = ',offsetpath  !hli(10/04/17)
          
      case default
        foundcard = .false.
          
      end select
    enddo
    
    call unitconv_var(ampunits,'m',ampfac)
    call unitconv_var(phaunits,'rad',phafac) 
    
    allocate(amp(ntc),speed(ntc))
    allocate(phase(ntc),name(ntc))
    allocate(f(ntc),vu(ntc))
    do k=1,ntc
      kth=ind(k)
      amp(k)   = amptemp(kth)*ampfac    !Convert units if necessary
      phase(k) = phasetemp(kth)*phafac  !Convert units if necessary
      speed(k) = speedtemp(kth)*deg2rad !Convert from degrees to radians
      f(k)     = ftemp(kth)
      vu(k)    = vutemp(kth)*deg2rad    !Convert from degrees to radians
      name(k)  = nametemp(kth)
    enddo
    
    return
    end subroutine tidal_block
    
!************************************************************************
    subroutine harmonic_block(ibndtype,ntc,amp,phase,speed,angle_wave)
! Reads an input harmonic constituent block from the card file
!
! written by Alex Sanchez, USACE-CHL
!************************************************************************
    use geo_def, only: azimuth_fl
    use const_def, only: deg2rad
    use comvarbl, only: mpfile,flowpath
    use bnd_def, only: ntf
    use unitconv_lib, only: unitconv_var
    use prec_def, only: ikind
    
    implicit none
    !Input/Output
    integer,intent(inout) :: ibndtype
    integer              :: ntc         !Tidal constituents used
    real(ikind)          :: angle_wave  !Incident angle of tidal wave
    real(ikind), pointer :: amp(:)      !Amplitude [m] (constituent) 
    real(ikind), pointer :: phase(:)    !Phase [rad] (constituent)
    real(ikind), pointer :: speed(:)    !Speed [rad/hrs] (constituent)
    !Internal Variables
    integer :: k,kk,ierr
    character(len=37) :: cardname
    character(len=6) :: ampunits,phaunits,spdunits
    real(ikind) :: ampfac,phafac,spdfac
    real(ikind) :: amptemp(ntf),phasetemp(ntf),speedtemp(ntf)
    logical :: foundcard
    
    ibndtype = 2 !=TH
    ampunits = 'm'
    phaunits = 'deg'
    spdunits = 'deg'
    ntc = 0 !Number of tidal constituents used
    do !  Removed the limit of only reading 20 cards in this block.  !kk=1,20  MEB 06/27/2018
      foundcard = .true.  
      read(77,*,iostat=ierr) cardname
      if(ierr/=0) exit
      if(cardname(1:1)=='!' .or. cardname(1:1)=='#' .or. cardname(1:1)=='*') cycle
      select case(cardname)
      case('HARMONIC_CONSTITUENTS_END','HARMONIC_END','END')
        exit
          
      case('AMPLITUDE_UNITS')
        backspace(77)
        read(77,*) cardname,ampunits
          
      case('PHASE_UNITS')
        backspace(77)
        read(77,*) cardname,phaunits
      
      case('SPEED_UNITS')
        backspace(77)
        read(77,*) cardname,spdunits
      
      case('INCIDENT_ANGLE','INCIDENT_DIRECTION')
        backspace(77)
        read(77,*) cardname,angle_wave              !degrees
        angle_wave = 90.0 - angle_wave - azimuth_fl !degrees, clockwise from North, local coordinate system
        if(angle_wave<0.0)then
          angle_wave = angle_wave + 360.0
        elseif(angle_wave>360.0)then
          angle_wave = angle_wave - 360.0
        endif
        angle_wave = angle_wave*deg2rad !Convert to rad
        
      case('HARMONIC','CONSTITUENT')
        ntc = ntc + 1
        backspace(77)   
        read(77,*) cardname,speedtemp(ntc),amptemp(ntc),phasetemp(ntc) !Preferred format, speed in deg/hr
        !read(77,*) cardname,amptemp(ntc),phasetemp(ntc),speedtemp(ntc)   !Current format in SMS 11.1, speed in cycles/hr
        !speedtemp(ntc) = speedtemp(ntc)*180.0 !convert from cycles/hr to deg/hr
        ibndtype = 2
        
      case default
        foundcard = .false.
          
      end select
    enddo
    
    call unitconv_var(ampunits,'m',ampfac)
    call unitconv_var(phaunits,'rad',phafac) 
    call unitconv_var(spdunits,'rad',spdfac) 
    
    allocate(amp(ntc),speed(ntc),phase(ntc))
    do k=1,ntc
      amp(k)   = amptemp(k)*ampfac    !Convert units if necessary
      phase(k) = phasetemp(k)*phafac  !Convert units if necessary
      speed(k) = speedtemp(k)*spdfac  !Convert units if necessary
    enddo
    
    return
    end subroutine harmonic_block

!*************************************************************************************
    subroutine tsta_block(ibndtype,station,ntc,name,amp,phase,speed,f,vu,angle_wave)
! Reads a tidal constituent block from the card file
!
! written by Alex Sanchez, USACE-CHL
!************************************************************************************
    use geo_def, only: azimuth_fl
    use const_def, only: deg2rad
    use comvarbl, only: mpfile,flowpath,iyr
    use bnd_def, only: ntf
    use tide_lib, only: tidal_data
    use unitconv_lib, only: unitconv_var
    use diag_lib, only: diag_print_error
    use prec_def, only: ikind
    
    implicit none
    !Input/Output
    integer,intent(inout) :: ibndtype
    character(len=*)     :: station
    integer              :: ntc         !Tidal constituents used
    real(ikind)          :: angle_wave  !Incident angle of tidal wave
    real(ikind), pointer :: amp(:)      !Amplitude [m] (constituent) 
    real(ikind), pointer :: speed(:)    !Speed [rad/hrs] (constituent)
    real(ikind), pointer :: phase(:)    !Phase [rad] (constituent)
    real(ikind), pointer :: f(:)        !Nodal factor [-] (constituent)
    real(ikind), pointer :: vu(:)       !Equilibrium argument [rad] (constituent)
    character(len=10), pointer :: name(:)     !Tidal Consitituent names (constituent)    
    !Internal Variables
    integer :: i,k,kk,ierr,id,kunit,ind,idsta,kerr
    integer :: mcyctemp(ntf)
    character(len=37) :: cardname
    character(len=100) :: namesta,tstafile,tstaname,tstapath,stationlower
    real(ikind) :: lonbnd,latbnd,lonsta,latsta,dist,distmin
    real(ikind) :: amptemp(ntf),phasetemp(ntf),speedtemp(ntf),ftemp(ntf),vutemp(ntf)
    real(ikind),parameter :: amptol = 0.001
    logical :: foundcard,instaname
    character(len=10) :: nametemp(ntf)

    !Get constituent names, speeds, nodal corrections and equilibrium phases (relative to Greenwich)
    call tidal_data(iyr,nametemp,speedtemp,mcyctemp,ftemp,vutemp)
    
    !Initialize
    ibndtype = 2 !=TH
    station = ''
    id = 0
    lonbnd = -999.0
    latbnd = -999.0
    tstaname = 'station_database.txt'
    tstapath = ''
    ntc = 0
    
    do !  Removed the limit of only reading 10 cards in this block.  !kk=1,10  MEB 06/27/2018
      foundcard = .true.  
      read(77,*,iostat=ierr) cardname
      if(ierr/=0) exit
      if(cardname(1:1)=='!' .or. cardname(1:1)=='#' .or. cardname(1:1)=='*') cycle
      select case(cardname)  
      case('TIDAL_STATION_END','STATION_END','END')
        exit
        
       case('DATABASE_NAME','STATION_DATABASE')
        backspace(77)
        read(77,*) cardname,tstaname
        
      case('DATABASE_PATH')
        backspace(77)
        read(77,*) cardname,tstapath
        
      case('DATABASE_FILE','DATABASE')
        backspace(77)
        read(77,*) cardname,tstafile 
        
      case('STATION_NAME','NAME')
        backspace(77)
        read(77,*) cardname,station
      
      case('STATION_ID','ID','ID_NUMBER')
        backspace(77)
        read(77,*) cardname,id
      
      case('LONGITUDE')
        backspace(77)
        read(77,*) cardname,lonbnd  
      
      case('LATITUDE')
        backspace(77)
        read(77,*) cardname,latbnd  
        
      case('CONSTITUENTS')
        backspace(77)
        read(77,*) cardname,ntc
        allocate(name(ntc))
        backspace(77)
        read(77,*) cardname,ntc,(name(k),k=1,ntc)
      
      case('INCIDENT_ANGLE','INCIDENT_DIRECTION')
        backspace(77)
        read(77,*) cardname,angle_wave  !degrees
        angle_wave = 90.0 - angle_wave - azimuth_fl !degrees, clockwise from North, local coordinate system
        if(angle_wave<0.0)then
          angle_wave = angle_wave + 360.0
        elseif(angle_wave>360.0)then
          angle_wave = angle_wave - 360.0
        endif
        angle_wave = angle_wave*deg2rad !Convert to rad  
        
      case default
        foundcard = .false.
            
      end select
    enddo
    
    tstafile = trim(tstapath) // trim(tstaname)
    stationlower = station
    call lowercase(stationlower)
    
    !Search file
    kunit = 932
    open(kunit,file=tstafile)
    !Skip header
    do k=1,3
      read(kunit,*) 
    enddo
    ierr = -1
    distmin = 1.0e6
    if(len_trim(station)>0)then
      instaname = .true.
    else
      instaname = .false.
    endif
    
    do
      read(kunit,*,iostat=kerr) ind,namesta,idsta,latsta,lonsta
      if(kerr/=0)then 
        exit
      endif
      if(id>0)then !ID
        if(id==idsta)then  
          backspace(kunit)
          read(kunit,*,iostat=ierr) ind,namesta,idsta,latsta,lonsta,(phasetemp(i),amptemp(i),i=1,ntf)
          station = namesta
          exit
        endif
      elseif(instaname)then !Name
        call lowercase(namesta)
        kk = index(namesta,stationlower)
        if(kk>0)then
          backspace(kunit)
          read(kunit,*,iostat=ierr) ind,namesta,idsta,latsta,lonsta,(phasetemp(i),amptemp(i),i=1,ntf)
          exit
        endif
      elseif(abs(lonbnd+999.0)>1.0e-3 .and. abs(latbnd+999.0)>1.0e-3)then !Specified Coordinates
        dist = sqrt((lonbnd-lonsta)**2+(latbnd-latsta)**2)
        if(dist<distmin)then
          distmin = dist
          station = namesta
          backspace(kunit)
          read(kunit,*,iostat=ierr) ind,namesta,idsta,latsta,lonsta,(phasetemp(i),amptemp(i),i=1,ntf)
        endif
      else !Use closest station to boundary   !Missing instructions - found 10/30/19 MEB
        continue  
      endif
    enddo
    
    if(ierr/=0)then
      call diag_print_error('Problem reading Station Constituents')
    endif
    
    if(ntc==0)then
      !Find and count significant amplitudes
      do k=1,ntf
        if(amptemp(k)>=amptol)then
          ntc = ntc + 1
        endif
      enddo
      allocate(name(ntc))
      i = 0
      do k=1,ntf
        if(amptemp(k)>=amptol)then
          i = i + 1
          name(i) = nametemp(k)
        endif
      enddo 
    endif  
    allocate(amp(ntc),phase(ntc),speed(ntc),f(ntc),vu(ntc))
      do k=1,ntc
        do i=1,ntf
          if(name(k)==nametemp(i))then
            amp(k) = amptemp(i)
            phase(k) = phasetemp(i)*deg2rad
            speed(k) = speedtemp(i)*deg2rad
            f(k) = ftemp(i)
            vu(k) = vutemp(i)*deg2rad
            exit
          endif     
        enddo
      enddo
    
    return
    end subroutine tsta_block
    
!*****************************************************************************    
    subroutine parent_block(ctlfilepar,grdfilepar,projpar,&
      wsefilepar,wsepathpar,velfilepar,velpathpar,tjuldaypar,timestarthr,nti)
! Reads a parent simulation block from the card file
!
! written by Alex Sanchez, USACE-CHL
!*****************************************************************************
    use geo_def, only: projection
    use const_def, only: deg2rad
    use comvarbl, only: mpfile,flowpath,tjulday0
    use time_lib, only: calendar2julian
    use prec_def, only: ikind
    
    implicit none
    !Input/Output
    integer,intent(inout) :: nti    
    real(ikind),intent(inout) :: tjuldaypar    !Parent simulation reference (starting) time in Julian days
    real(ikind),intent(inout) :: timestarthr   !Time relative to the CMS-Flow starting time (used for ADCIRC)
    character(len=200),intent(inout) :: wsefilepar,wsepathpar    !Water level data file and path
    character(len=200),intent(inout) :: velfilepar,velpathpar    !Velocity data file and path
    character(len=200),intent(inout) :: ctlfilepar               !Parent control file and path
    character(len=200),intent(inout) :: grdfilepar               !Parent grid file and path
    type(projection),intent(inout) :: projpar
    !Internal Variables
    integer :: kk,ierr
    integer :: iyrpar,imopar,idaypar,ihrpar,iminpar,isecpar
    character :: cardname*37,aext*10
    character(len=200) :: mpfilepar,flowpathpar,flownamepar
    logical :: foundcard
    
    do !  Removed the limit of only reading 10 cards in this block.  !kk=1,10  MEB 06/27/2018
      foundcard = .true.  
      read(77,*,iostat=ierr) cardname
      if(ierr/=0) exit
      if(cardname(1:1)=='!' .or. cardname(1:1)=='#' .or. cardname(1:1)=='*') cycle
      select case(cardname)        
      case('PARENT_END','PARENT_SIMULATION_END','END')
        exit
          
      case('CONTROL_FILE','PARENT_CONTROL_FILE','PARENT_CTL_FILE','CTL_FILE')
        backspace(77)
        read(77,*) cardname,ctlfilepar !,NH_str(nNHstr)%ctlpath
        call fileparts(ctlfilepar,flowpathpar,flownamepar,aext)
    
      case('GRID_FILE','PARENT_GRID_FILE')
        backspace(77)
        read(77,*) cardname,grdfilepar
        
      case('PARAMS_FILE','PARAMETERS_FILE','PARENT_PARAMETERS_FILE')
        backspace(77)
        read(77,*) cardname,mpfilepar
      
      case('SOL_FILE','SOLUTION_FILE')
        backspace(77)
        read(77,*) cardname,wsefilepar  
        velfilepar = wsefilepar  
        
      !case('WSE_FILE','PARENT_WSE_FILE','WSE_OUT_FILE','WSE_SOL_FILE')
      !  backspace(77)
      !  read(77,*) cardname,wsefilepar  
      !  
      !case('VEL_FILE','PARENT_VEL_FILE','VEL_OUT_FILE','VEL_SOL_FILE')
      !  backspace(77)
      !  read(77,*) cardname,velfilepar
        
      case('WSE_SOLUTION','WSE_DATASET','WSE_FILE','PARENT_WSE_FILE','WSE_OUT_FILE','WSE_SOL_FILE')
        backspace(77)                                 !added 6/10/21 MEB
        read(77,*) cardname, mpfilepar, flowpathpar   !added 6/10/21 MEB
        call card_dataset(77,mpfilepar,flowpathpar,wsefilepar,wsepathpar,1)
        
      case('VEL_SOLUTION','VEL_DATASET','VEL_FILE','PARENT_VEL_FILE','VEL_OUT_FILE','VEL_SOL_FILE')
        backspace(77)                                 !added 6/10/21 MEB
        read(77,*) cardname, mpfilepar, flowpathpar   !added 6/10/21 MEB
        call card_dataset(77,mpfilepar,flowpathpar,velfilepar,velpathpar,2)
        
      case('STARTING_DATE_TIME','START_DATE_TIME','PARENT_STARTING_DATE_TIME')
        call card_datetime(77,iyrpar,imopar,idaypar,ihrpar,iminpar,isecpar)  !YYYY-MM-DD HH:MM:SS UTC 
        call calendar2julian(iyrpar,imopar,idaypar,ihrpar,iminpar,isecpar,tjuldaypar)
        timestarthr = (tjuldaypar - tjulday0)*24.0
            
      case('RELATIVE_START_TIME','REFERENCE_TIME','PARENT_STARTING_TIME')
        call card_scalar(77,'hrs','hrs',timestarthr,ierr)
        tjuldaypar = tjulday0 + timestarthr/24.0
        
      case('TEMPORAL_INTERPOLATION_ORDER','ORDER','ORDER_INTERP')
        backspace(77)
        read(77,*) cardname,nti
        nti = min(max(nti,1),3) !Limit between 1 and 3
      
      case('HORIZONTAL_PROJECTION_BEGIN','HORIZ_PROJ_BEGIN')
        call proj_horiz_block(77,projpar)
        
      case default
        foundcard = .false.
          
      end select
    enddo
    
    return
    end subroutine parent_block    

!********************************************************************    
    subroutine tdb_block(ntcin,namein,tdbname,tdbpath,projtdb,nssi,nssw)
! Reads a tidal database block from the card file
! written by Alex Sanchez, USACE-CHL
!********************************************************************
    use geo_def, only: azimuth_fl,projection
    use const_def, only: deg2rad
    use comvarbl, only: mpfile,flowpath

    implicit none
    !Input/Output
    integer,intent(out)            :: ntcin    !Tidal constituents used     
    character(len=*),intent(inout), pointer :: namein(:)  !Input Tidal Consitituent names (constituent)
    character(len=*),intent(inout) :: tdbname  !Tidal Database Name, EC2001, ENPAC2003, LEPROVOST, 
    character(len=*),intent(inout) :: tdbpath  !Tidal Database file and path
    type(projection),intent(inout) :: projtdb  !Parent grid projection
    integer,intent(inout)          :: nssi     !Smoothing iterations (along string)
    integer,intent(inout)          :: nssw     !Smoothing window width (along string)
    
    !Internal Variables
    integer :: k,kk,ierr,nn
    character :: cardname*37
    logical :: foundcard
    
    do !  Removed the limit of only reading 10 cards in this block.  !kk=1,10  MEB 06/27/2018
      foundcard = .true.  
      read(77,*,iostat=ierr) cardname
      if(ierr/=0) exit
      if(cardname(1:1)=='!' .or. cardname(1:1)=='#' .or. cardname(1:1)=='*') cycle
      select case(cardname)
      case('TIDAL_DATABASE_END','DATABASE_END','END')
        exit
          
      case('CONSTITUENTS')
        backspace(77)
        read(77,*) cardname,ntcin
        allocate(namein(ntcin))
        backspace(77)
        read(77,*) cardname,ntcin,(namein(k),k=1,ntcin)
      
      case('TIDAL_DATABASE_NAME','DATABASE_NAME','NAME')
        backspace(77)
        read(77,*) cardname,tdbname         
        
      case('TIDAL_DATABASE_PATH','DATABASE_PATH','PATH')
        backspace(77)
        read(77,*) cardname,tdbpath
        nn = len_trim(tdbpath)
        if(tdbpath(nn:nn)/='\' .and. tdbpath(nn:nn)/='/')then
          tdbpath = trim(tdbpath) // '\'
        endif
        
       case('HORIZONTAL_PROJECTION_BEGIN','HORIZ_PROJ_BEGIN')
        call proj_horiz_block(77,projtdb)    
       
       case('SPACE_SMOOTH_ITER','SPACE_SMOOTH_ITERATIONS','SPATIAL_SMOOTHING_ITERATIONS',&
        'SMOOTHING_ITERATIONS','SMOOTH_ITERATIONS','SMOOTHING_ITER')
        backspace(77)
        read(77,*) cardname,nssi
        
      case('SPACE_SMOOTH_WIDTH','SPATIAL_SMOOTHING_WIDTH','SMOOTHING_WIDTH','SMOOTH_WIDTH')
        backspace(77)
        read(77,*) cardname,nssw
        if(mod(nssw,2)/=0)then !Make sure it is odd
          nssw = nssw + 1
        endif    
          
      case default
        foundcard = .false.
          
      end select
    enddo
    
    return
    end subroutine tdb_block
    
!**********************************************************    
    subroutine tidal_alloc
!**********************************************************
    use bnd_def, only: TH_str,nTHstr,ntf,TH_type,ioffsetmode !(hli 10/04/17)

    implicit none
    integer :: i,ntc
    type(TH_type), allocatable :: Tstrtemp(:)
    
    nTHstr = nTHstr + 1
    if(nTHstr==1)then
      if(.not.allocated(TH_str))then
        allocate(TH_str(nTHstr))
      endif
    else
      allocate(Tstrtemp(nTHstr-1))
      do i=1,nTHstr-1        
        ntc = TH_str(i)%ntc
        allocate(Tstrtemp(i)%amp(ntc),Tstrtemp(i)%speed(ntc))
        allocate(Tstrtemp(i)%phase(ntc),Tstrtemp(i)%name(ntc))
        allocate(Tstrtemp(i)%f(ntc),Tstrtemp(i)%vu(ntc))
        Tstrtemp(i) = TH_str(i)
      enddo
      deallocate(TH_str)
      allocate(TH_str(nTHstr))
      do i=1,nTHstr-1
        ntc = Tstrtemp(i)%ntc      
        allocate(TH_str(i)%amp(ntc),TH_str(i)%speed(ntc))
        allocate(TH_str(i)%phase(ntc),TH_str(i)%name(ntc))
        allocate(TH_str(i)%f(ntc),TH_str(i)%vu(ntc))
        TH_str(i) = Tstrtemp(i)
      enddo
      deallocate(Tstrtemp)
    endif
    
    !Initialize
    TH_str(nTHstr)%ioffsetmode = ioffsetmode  !1-Constant offset, 2-Offset curve (hli,01/18/17)
    TH_str(nTHstr)%offsetfile = ''
    TH_str(nTHstr)%offsetpath = ''    
    TH_str(nTHstr)%wsecurveoffset = 0.0
    TH_str(nTHstr)%ntc = 0
    TH_str(nTHstr)%nti = 2     !hli(10/06/17)
    TH_str(nTHstr)%inc = 1     !hli(10/06/17)
    TH_str(nTHstr)%angle = -999.0
    TH_str(nTHstr)%dwsex = 0.0
    TH_str(nTHstr)%dwsey = 0.0
    TH_str(nTHstr)%wseoffset = 0.0    
    TH_str(nTHstr)%istidal = .true.
    TH_str(nTHstr)%wseadjust = .true.
    TH_str(nTHstr)%station = ''
    
    return
    end subroutine tidal_alloc
    
!*************************************************************    
    subroutine singlewse_alloc
! Resizes the single water level boundary condition variable    
!*************************************************************
    use bnd_def, only: H_str,nHstr,H_type,ioffsetmode

    implicit none
    integer :: i
    type(H_type), allocatable :: Hstrtemp(:)    
    
    nHstr = nHstr + 1
    if(nHstr==1)then  
      allocate(H_str(nHstr))
    else
      allocate(Hstrtemp(nHstr-1))
      do i=1,nHstr-1
        Hstrtemp(i) = H_str(i)
      enddo
      deallocate(H_str)
      allocate(H_str(nHstr))
      do i=1,nHstr-1        
        H_str(i) = Hstrtemp(i)
      enddo
      deallocate(Hstrtemp)
    endif
    
    !Initialize and set default values
!    write(3000,*)'ioffsetmode (singlewse) = ',ioffsetmode 
    H_str(nHstr)%ioffsetmode = ioffsetmode  !1-Constant offset, 2-Offset curve (hli,01/18/17)
    H_str(nHstr)%offsetfile = ''
    H_str(nHstr)%offsetpath = ''    
    H_str(nHstr)%wsecurveoffset = 0.0
    H_str(nHstr)%wseoffset = 0.0
    H_str(nHstr)%dwsex = 0.0
    H_str(nHstr)%dwsey = 0.0
    H_str(nHstr)%bidfile = ''
    H_str(nHstr)%bidpath = ''
    H_str(nHstr)%wsefile = ''
    H_str(nHstr)%wsepath = ''
    H_str(nHstr)%ncells = 0
    H_str(nHstr)%ntimes = 0
    H_str(nHstr)%inc = 1
    H_str(nHstr)%nti = 2 !Order of temporal interpolation
    H_str(nHstr)%nsw = 0 !Temporal smoothing width
    H_str(nHstr)%nsi = 0 !Temporal smoothing iterations
    H_str(nHstr)%wseadjust = .true.
    !H_str(nHstr)%minterp    !Method for interpolation, 1-Piecewise polynomial, 2-cubic spline
    
    return
    end subroutine singlewse_alloc

!****************************************************************    
    subroutine multiwse_alloc
! Resizes the multiple water level boundary condition variables    
!****************************************************************
    use bnd_def, only: MH_str,nMHstr,MH_type

    implicit none
    integer :: i
    type(MH_type), allocatable :: MHstrtemp(:)    
    
    nMHstr = nMHstr + 1
    if(nMHstr==1)then
      allocate(MH_str(nMHstr))
    else
      allocate(MHstrtemp(nMHstr-1))
      do i=1,nMHstr-1
        MHstrtemp(i) = MH_str(i)
      enddo
      deallocate(MH_str)
      allocate(MH_str(nMHstr))
      do i=1,nMHstr-1        
        MH_str(i) = MHstrtemp(i) 
      enddo
      deallocate(MHstrtemp)
    endif    
    
    !Initialize values
    MH_str(nMHstr)%wseoffset = 0.0
    MH_str(nMHstr)%bidfile = ''
    MH_str(nMHstr)%bidpath = ''
    MH_str(nMHstr)%wsefile = ''
    MH_str(nMHstr)%wsepath = ''
    MH_str(nMHstr)%ncells = 0
    MH_str(nMHstr)%ntimes = 0
    MH_str(nMHstr)%inc = 1
    MH_str(nMHstr)%nti = 2  !Temporal interpolation order
    MH_str(nMHstr)%nsw = 0  !Temporal smoothing width
    MH_str(nMHstr)%nssi = 0  !Spatial smoothing iterations
    MH_str(nMHstr)%nssw = 0 !Spatial smoothing width
    
    return
    end subroutine multiwse_alloc

!**********************************************************    
    subroutine multiwsevel_alloc
! Resizes the multiple water level and velocity boundary
! condition variable
!**********************************************************
    use bnd_def, only: MHV_str,nMHVstr,MHV_type

    implicit none
    integer :: i
    type(MHV_type), allocatable :: MHVstrtemp(:)    
    
    nMHVstr = nMHVstr + 1
    if(nMHVstr==1)then
      allocate(MHV_str(nMHVstr))
    else
      allocate(MHVstrtemp(nMHVstr-1))
      do i=1,nMHVstr-1
        MHVstrtemp(i)= MHV_str(i)
      enddo
      deallocate(MHV_str)
      allocate(MHV_str(nMHVstr))
      do i=1,nMHVstr-1        
        MHV_str(i) = MHVstrtemp(i)        
      enddo
      deallocate(MHVstrtemp)
    endif    
    
    !Initialize values
    MHV_str(nMHVstr)%wseoffset = 0.0 !Can be useful to convert between datums
    MHV_str(nMHVstr)%bidfile = ''
    MHV_str(nMHVstr)%bidpath = ''
    MHV_str(nMHVstr)%wsefile = ''
    MHV_str(nMHVstr)%wsepath = ''
    MHV_str(nMHVstr)%velfile = ''
    MHV_str(nMHVstr)%velpath = ''
    MHV_str(nMHVstr)%ncells = 0
    !Water level options
    MHV_str(nMHVstr)%ntimeswse = 0    
    MHV_str(nMHVstr)%incwse = 1
    MHV_str(nMHVstr)%ntiwse = 2  !Temporal interpolation order
    MHV_str(nMHVstr)%nswwse = 0  !Temporal smoothing width
    MHV_str(nMHVstr)%nsiwse = 0  !Temporal smoothing iterations
    MHV_str(nMHVstr)%nssiwse = 0  !Spatial smoothing iterations
    MHV_str(nMHVstr)%nsswwse = 0 !Spatial smoothing width
    !Current velocity options
    MHV_str(nMHVstr)%ntimesvel = 0
    MHV_str(nMHVstr)%incvel = 1
    MHV_str(nMHVstr)%ntivel = 2  !Temporal interpolation order
    MHV_str(nMHVstr)%nswvel = 0  !Temporal smoothing width
    MHV_str(nMHVstr)%nsivel = 0  !Temporal smoothing iterations
    MHV_str(nMHVstr)%nssivel = 0  !Spatial smoothing iterations
    MHV_str(nMHVstr)%nsswvel = 0 !Spatial smoothing width
    
    return
    end subroutine multiwsevel_alloc

!****************************************************************
    subroutine read_tidal_const()
! Reads the tidal constituent information from the CMS control file
! written by Alex Sanchez, USACE-ERDC-CHL
!
! Incorporate tidal prediction program with the nodal factors.
! A total of 37 tidal contituents can be used for the prediction 
! (hli, 02/03/10).
!****************************************************************
    use bnd_def, only: TH_str,ntf,TH_type,nTHstr,Tread,nbndstr
    use comvarbl, only: iyr
    use tide_lib, only: tidal_data
    use const_def, only: deg2rad,twopi
    use prec_def, only: ikind
    
    implicit none
    integer :: i,k,nn,kth,ntc,ierr
    integer :: mcyctemp(ntf),ind(ntf)
    real(ikind) :: amptemp(ntf),phasetemp(ntf),speedtemp(ntf),ftemp(ntf),vutemp(ntf)
    character(len=10) :: nametemp(ntf)
    character(len=10) :: nameyr(ntf),astr
    character(len=32) :: cardname    
    
!   Get constituent names, speeds, nodal corrections and equilibrium phases (relative to Greenwich)
    call tidal_data(iyr,nameyr,speedtemp,mcyctemp,ftemp,vutemp)

    !backspace(77)   
    !read(77,*) cardname    !SKIP DOWN TO NEXT LINE   
             
    ntc = 0 !Number of tidal constituents used
d1: do i=1,ntf
      read(77,*,iostat=ierr) cardname
      if(ierr/=0)then
        return
      endif
      backspace(77)
      nn = len_trim(cardname)
      astr = cardname(19:nn)
d2:   do k=1,ntf
        if(astr==nameyr(k))then
          kth = k
          exit d2
        elseif(astr(1:4)=='_END')then
          exit d1
       endif  
      enddo d2
      ntc = ntc + 1
      ind(ntc) = kth      
      read(77,*) cardname, amptemp(kth), phasetemp(kth)
      nametemp(kth) = astr
    enddo d1        
    if(.not.Tread)then
      call tidal_alloc
    endif       
    TH_str(nTHstr)%ntc = ntc
    allocate(TH_str(nTHstr)%amp(ntc),TH_str(nTHstr)%speed(ntc))
    allocate(TH_str(nTHstr)%phase(ntc),TH_str(nTHstr)%name(ntc))
    allocate(TH_str(nTHstr)%f(ntc),TH_str(nTHstr)%vu(ntc))
    do k=1,TH_str(nTHstr)%ntc
      kth=ind(k)
      TH_str(nTHstr)%amp(k)   = amptemp(kth)      
      TH_str(nTHstr)%speed(k) = speedtemp(kth)*deg2rad !Convert from degrees to radians
      TH_str(nTHstr)%phase(k) = phasetemp(kth)*deg2rad !Convert from degrees to radians
      TH_str(nTHstr)%f(k)     = ftemp(kth)
      TH_str(nTHstr)%vu(k)    = vutemp(kth)*deg2rad    !Convert from degrees to radians
      TH_str(nTHstr)%name(k)  = nametemp(kth)
    enddo
    write(*,*) 'Tidal Constituents read'  
    
    return
    end subroutine read_tidal_const    
    
!****************************************************************
    subroutine read_harmonics
! Reads the tidal constituent information from the CMS control file
! written by Alex Sanchez, USACE-ERDC-CHL
!
! Incorporate tidal prediction program with the nodal factors.
! A total of 37 tidal contituents can be used for the prediction 
! (hli, 02/03/10).
!****************************************************************
    use bnd_def, only: TH_str,ntf,TH_type,nTHstr,Tread,nbndstr
    use const_def, only: deg2rad,twopi
    use prec_def, only: ikind
    
    implicit none
    integer :: i,k,ntc,ierr
    real(ikind) :: amptemp(ntf),phasetemp(ntf),speedtemp(ntf)
    character(len=32) :: cardname
    
    ntc = 0 !Number of tidal constituents used
d1: do i=1,ntf
      read(77,*,iostat=ierr) cardname
      backspace(77)
      if(ierr/=0) exit
      if(cardname(1:1)=='!' .or. cardname(1:1)=='#' .or. cardname(1:1)=='*') cycle
      select case(cardname)
      case('HARMONIC_END','HARMONICS_END','HARMONIC_COMPONENTS_END','END')  
        exit  
      case default      
        ntc = ntc + 1   
        if(cardname(1:4)=='HARM')then !SMS 11.1 format          
          read(77,*,iostat=ierr) cardname, amptemp(ntc),phasetemp(ntc),speedtemp(ntc)     
          if(ierr/=0)then
            ntc = ntc - 1
            exit
          endif  
          speedtemp(ntc) = speedtemp(ntc)*180.0 !Convert from cycles/hr to deg/hr
        else  !Preferred format 
          read(77,*,iostat=ierr) cardname, speedtemp(ntc),amptemp(ntc),phasetemp(ntc)   !Speed in deg/hr, phase in degrees
          if(ierr/=0)then
            ntc = ntc - 1
            exit
          endif  
        endif  
      end select
    enddo d1        
    if(.not.Tread)then
      call tidal_alloc
    endif      
    TH_str(nTHstr)%istidal = .false.
    TH_str(nTHstr)%ntc = ntc
    allocate(TH_str(nTHstr)%amp(ntc),TH_str(nTHstr)%speed(ntc))
    allocate(TH_str(nTHstr)%phase(ntc))
    do k=1,TH_str(nTHstr)%ntc
      TH_str(nTHstr)%amp(k)   = amptemp(k)      
      TH_str(nTHstr)%speed(k) = speedtemp(k)*deg2rad !Convert from degrees to radians
      TH_str(nTHstr)%phase(k) = phasetemp(k)*deg2rad !Convert from degrees to radians
    enddo
    write(*,*) 'Harmonic components read'  
    
    return
    end subroutine read_harmonics
    
!********************************************************************************
    subroutine bnd_init()
! Reads the boundary condition information from *_mp.h5 file and 
! calculates the spatial (grid) boundary condition variables needed for
! the model from the variables which are used to read from the xmdf files
!
! written by Alex Sanchez, USACE-CHL
!********************************************************************************
#include "CMS_cpp.h"
    use bnd_def  
    use bnd_lib, only: read_bndstr,read_fluxdata,read_snglwsedata,read_multiwsedata,read_multiveldata, read_offsetwsedata !(hli,01/20/17)
    use cms_def, only: cmswave
    use comvarbl, only: tjulday0
    use const_def, only: deg2rad,rad2deg
    use diag_def
    use diag_lib, only: diag_print_error, diag_print_warning
    use flow_def, only: grav
    use geo_def, only: ncface,cell2cell,idirface,x,y,zb,rx,ry,xc,yc,projfl,mapid
    use geo_lib, only: read_grid14,proj_horiz_conv
    use interp_lib, only: interp_coef_tel2pts,interp_coef_tri2pts
    use nest_lib
    use prec_def, only: ikind
    use size_def, only: ncellsd,ncellpoly,ncells
    !use struct_def
    !use tide_lib
    
    implicit none      
    integer :: i,ii,j,k,im,nbndcells,nck,mntp
    integer :: iriv,iwse,icsh,korient,id1,id2,ibnd,ipar,idpar,kk
    integer :: itide    !adding a descriptive loop control variable for the tidal and harmonic cell strings
    character(len=10) :: aext
    real(ikind) :: dep,cosang,sinang,val,distx,disty,tjuldaypar,xtrapdist
    integer :: ibndtemp(ncellsD)
    real(ikind), allocatable :: valtemp(:,:)
    type(MH_type),allocatable :: MH_temp(:)

!Remove repeated Multiple Water Level BC
    do i=1,nMHVstr
      j=1
      allocate(MH_temp(nMHstr))
      do while(j<=nMHstr)
        if(MHV_str(i)%bidpath==MH_str(j)%bidpath)then
          ii=0
          do k=1,nMHstr
            if(j/=k)then
              ii=ii+1
              MH_temp(ii) = MH_str(k)
            endif
          enddo
          deallocate(MH_str)
          nMHstr=nMHstr-1
          allocate(MH_str(nMHstr))
          do k=1,nMHstr
            MH_str(k) = MH_temp(k)
          enddo
        endif
        j=j+1
      enddo 
      deallocate(MH_temp)
    enddo
    
!Check tidal constituent boundary
    if(Tread .and. Tide_read)then
      nTHstr = 1
    elseif(Tread .and. .not.Tide_read)then
      write(msg2,*)  'Tread= ',Tread,'  Tide_read=',Tide_read        
      call diag_print_error('Missing tidal constituents for cell string',msg2)
    elseif(.not.Tread .and. Tide_read)then
      nTHstr = 0
      call diag_print_warning('Missing cell string for tidal constituents',&
        '  Tidal constituents will be ignored')
    endif
    
!    if(nTstrchk/=nTHstr)then
!      call diag_print_error('Unassigned tidal constituent information')
!    endif

!--- River/Flux BC (Type 1=Q) -------------------------------
    do iriv=1,nQstr
      !Read cell/node string  
      call read_bndstr(Q_str(iriv)%bidfile,Q_str(iriv)%bidpath,&
        Q_str(iriv)%istrtype,Q_str(iriv)%idnum,&
        Q_str(iriv)%ncells,Q_str(iriv)%cells,Q_str(iriv)%faces)
      
      !Determine orientation of river boundary
      if(Q_str(iriv)%angle<-900.0)then
        im = Q_str(iriv)%ncells/2  !Reference ID
        i = Q_str(iriv)%cells(im) !Reference cell
        k = Q_str(iriv)%faces(im)  !Reference face 
        if(Q_str(iriv)%istrtype==1)then
          korient = idirface(k,i) !Bug fix, Alex, June 8, 2010, changed i to iid
          !If inflow angle is not specified assume normal to boundary
          Q_str(iriv)%angle = float(4-korient)*90.0 !Angle is in local coordinates            
        else
          Q_str(iriv)%angle = atan2(ry(k,i),rx(k,i))*rad2deg-180.0 !Inflow Angle
        endif
      endif
      
      !Unit Conversion, ifluxunits = 0-m^3/s/cell, 1-m^3/s/boundary, 2-ft^3/s/boundary, default = 0      
      if(Q_str(iriv)%ifluxunits==0)then
        val = Q_str(iriv)%ncells !from m^3/s/cell to m^3/s/boundary
      elseif(Q_str(iriv)%ifluxunits==2)then
        val = 0.3048**3 !from ft^3/s to m^3/s/boundary 
      else
        val = 1.0  !m^3/s/boundary 
      endif
      
      if(Q_str(iriv)%ifluxmode==2)then !Time-series
        !Read flux data           
        call read_fluxdata(Q_str(iriv)%fluxfile,Q_str(iriv)%fluxpath,&
             Q_str(iriv)%ntimes,Q_str(iriv)%times,Q_str(iriv)%qcurv)
        
        !Check Interpolation order
        Q_str(iriv)%nti = max(min(Q_str(iriv)%nti,Q_str(iriv)%ntimes-1),1)
        
        !Temporal smoothing by moving average
        if(Q_str(iriv)%nsi>0 .and. Q_str(iriv)%nsw>0)then
          call moving_average(Q_str(iriv)%nsi,Q_str(iriv)%nsw,&
               Q_str(iriv)%ntimes,Q_str(iriv)%qcurv)
        endif  
        Q_str(iriv)%qcurv = val*Q_str(iriv)%qcurv  !Initial condition is set to zero [m^3/s/boundary]
      elseif(Q_str(iriv)%ifluxmode==1)then !Constant flux
        Q_str(iriv)%qfluxconst = val*Q_str(iriv)%qfluxconst  !Initial condition is set to zero [m^3/s/boundary] 
      elseif(Q_str(iriv)%ifluxmode==3)then !Stage-Flux curve
        call read_fluxdata(Q_str(iriv)%fluxfile,Q_str(iriv)%fluxpath,&
             Q_str(iriv)%nstages,Q_str(iriv)%stage,Q_str(iriv)%rflow)
        !Note: need to deal with units here ************
      endif
    enddo !i-str
    
!--- Tidal/Harmonic BC (Type 2=TH) ------------------------------------------------
    do itide=1,nTHstr                                                     !Switched loop control from 'iwse' to 'itide' for clarity - MEB 01/28/2020
      !Read cell/node string  
      call read_bndstr(TH_str(itide)%bidfile,TH_str(itide)%bidpath,&
        TH_str(itide)%istrtype,TH_str(itide)%idnum,&
        TH_str(itide)%ncells,TH_str(itide)%cells,TH_str(itide)%faces)
      
      !Water level adjustment
      if(TH_str(itide)%ncells<3) TH_str(itide)%wseadjust = .false.
      
      !Initialize variables
      allocate(TH_str(itide)%wsebnd0(TH_str(itide)%ncells)) !Initial water levels at boundary (not necessary the same as forcing)
      allocate(TH_str(itide)%wsebnd(TH_str(itide)%ncells)) !Water levels at boundary (not necessary the same as forcing)
      allocate(TH_str(itide)%wsevar(TH_str(itide)%ncells))         
      TH_str(itide)%wsebnd0(:) = 0.0 !May be overwritten in hot_read if initial condition is specified
      TH_str(itide)%wsebnd(:)  = 0.0
      TH_str(itide)%wsevar(:)  = 0.0
      if(TH_str(itide)%wseadjust)then
        allocate(TH_str(itide)%wseadj(TH_str(itide)%ncells))  !Adjusted water levels with wind and wave setup     
        TH_str(itide)%wseadj(:)  = 0.0      
      endif
      
      !Calculate phase difference due incident angle and regional wse gradients
      !allocate(TH_str(nTHstr)%psi(TH_str(nTHstr)%ncells,TH_str(nTHstr)%ntc))     !meb 1/31/2020  I think these should reference each cell string's PSI instead of always the max
      !TH_str(nTHstr)%psi = 0.0
      allocate(TH_str(itide)%psi(TH_str(itide)%ncells,TH_str(itide)%ntc))
      TH_str(itide)%psi = 0.0
      !Calculate characteristic depth
      if(TH_str(itide)%angle>-900.0)then
        dep = 0.0
        do j=1,TH_str(itide)%ncells
          i = TH_str(itide)%cells(j)
          dep = dep - zb(i) + TH_str(itide)%wseoffset  !h is undefined at this point
        enddo
        dep = dep/real(TH_str(itide)%ncells,kind=ikind)      
      endif      
      im = TH_str(itide)%ncells/2
      i = TH_str(itide)%cells(im) !Reference cell
      cosang = cos(TH_str(itide)%angle)
      sinang = sin(TH_str(itide)%angle)
      do j=1,TH_str(itide)%ncells
        ii=TH_str(itide)%cells(j)
        distx=x(ii)-x(i); disty=y(ii)-y(i)
        TH_str(itide)%wsevar(j) = distx*TH_str(itide)%dwsex + disty*TH_str(itide)%dwsey
        if(TH_str(itide)%angle<-900.0) cycle
        val = (distx*cosang+disty*sinang)/sqrt(grav*dep)/3600.0
        do k=1,TH_str(itide)%ntc
          TH_str(itide)%psi(j,k) = TH_str(itide)%psi(j,k) + TH_str(itide)%speed(k)*val
        enddo
      enddo
      if(TH_str(itide)%ioffsetmode==2)then !Time-series offset (hli,10/04/17)
        !Read offset data           
        call read_offsetwsedata(TH_str(itide)%offsetfile,TH_str(itide)%offsetpath,&
             TH_str(itide)%ntimesoffset,TH_str(itide)%offsettimes,TH_str(itide)%offsetcurve)
!        write(3000,*)'TH_str(itide)%offsetfile = ',TH_str(itide)%offsetfile
!        write(3000,*)'TH_str(itide)%offsettimes = ',TH_str(itide)%offsettimes
!        write(3000,*)'TH_str(itide)%offsetcurve = ',TH_str(itide)%offsetcurve
      endif
    enddo !i-str   

!--- Single Water Level BC (Type 3=H) -------------------------------------------
    do iwse=1,nHstr
      !Read cell/node string  
      call read_bndstr(H_str(iwse)%bidfile,H_str(iwse)%bidpath,&
        H_str(iwse)%istrtype,H_str(iwse)%idnum,&
        H_str(iwse)%ncells,H_str(iwse)%cells,H_str(iwse)%faces)
      
      !Water level adjustment
      if(H_str(iwse)%ncells<3) H_str(iwse)%wseadjust = .false.
      
      !Read water level data
      if(len_trim(H_str(iwse)%wsefile)>0)then
        call read_snglwsedata(H_str(iwse)%wsefile,H_str(iwse)%wsepath,&
             H_str(iwse)%ntimes,H_str(iwse)%times,H_str(iwse)%wsecurv)     
        
        !Check Interpolation order
        H_str(iwse)%nti = max(min(H_str(iwse)%nti,H_str(iwse)%ntimes-1),1)
        
        !Temporal smoothing by moving average
        if(H_str(iwse)%nsi>0 .and. H_str(iwse)%nsw>0)then
          call moving_average(H_str(iwse)%nsi,H_str(iwse)%nsw,&
               H_str(iwse)%ntimes,H_str(iwse)%wsecurv)
        endif
      endif
      
      !Initialize variables
      allocate(H_str(iwse)%wsebnd0(H_str(iwse)%ncells))
      allocate(H_str(iwse)%wsebnd(H_str(iwse)%ncells))
      allocate(H_str(iwse)%wsevar(H_str(iwse)%ncells))      
      H_str(iwse)%wsebnd0(:) = 0.0 !May be overwritten in hot_read if initial condition is specified
      H_str(iwse)%wsebnd(:)  = 0.0
      H_str(iwse)%wsevar(:)  = 0.0
      if(H_str(iwse)%wseadjust)then
        allocate(H_str(iwse)%wseadj(H_str(iwse)%ncells))      
        H_str(iwse)%wseadj(:)  = 0.0
      endif  
      
      !Add Temporally-Constant Surface gradients
      if(abs(H_str(iwse)%dwsex)>1.0e-20 .or. abs(H_str(iwse)%dwsey)>1.0e-20)then
        im = H_str(iwse)%ncells/2 !Reference id
        i = H_str(iwse)%cells(im) !Reference cell
        do j=1,H_str(iwse)%ncells
          ii=H_str(iwse)%cells(j)
          distx=x(ii)-x(i); disty=y(ii)-y(i)
          H_str(iwse)%wsevar(j) = distx*H_str(iwse)%dwsex + disty*H_str(iwse)%dwsey !Spatially variable and temporally constant water level
        enddo
      endif
      
      if(H_str(iwse)%ioffsetmode==2)then !Time-series offset (hli,01/18/17)
        !Read offset data           
        call read_offsetwsedata(H_str(iwse)%offsetfile,H_str(iwse)%offsetpath,&
             H_str(iwse)%ntimesoffset,H_str(iwse)%offsettimes,H_str(iwse)%offsetcurve)
      endif
    enddo !iwse

!--- Multiple Water Level BC (Type 4=MH) -----------------------------------------------
    do iwse=1,nMHstr
      !Read cell/node string  
      call read_bndstr(MH_str(iwse)%bidfile,MH_str(iwse)%bidpath,&
        MH_str(iwse)%istrtype,MH_str(iwse)%idnum,&
        MH_str(iwse)%ncells,MH_str(iwse)%cells,MH_str(iwse)%faces)
             
      !Read water level data             
      call read_multiwsedata(MH_str(iwse)%wsefile,MH_str(iwse)%wsepath,&
             MH_str(iwse)%ncells,MH_str(iwse)%cells,&
             MH_str(iwse)%ntimes,MH_str(iwse)%times,MH_str(iwse)%wsedata)

      !Spatial smoothing  
      if(MH_str(iwse)%nssi>0 .and. MH_str(iwse)%nssw>2)then
        allocate(valtemp(MH_str(iwse)%ncells,1))  
        do k=1,MH_str(iwse)%ntimes
          valtemp(:,1) = MH_str(iwse)%wsedata(k,:)
          call moving_average(MH_str(iwse)%nssi,MH_str(iwse)%nssw,&
               MH_str(iwse)%ncells,valtemp(:,1))
          MH_str(iwse)%wsedata(k,:) = valtemp(:,1)
        enddo
        deallocate(valtemp)
      endif
      
      !Temporal smoothing by moving average
      if(MH_str(iwse)%nsi>0 .and. MH_str(iwse)%nsw>0)then
        do j=1,MH_str(iwse)%ncells
          call moving_average(MH_str(iwse)%nsi,MH_str(iwse)%nsw,&
                 MH_str(iwse)%ntimes,MH_str(iwse)%wsedata(:,j))
        enddo
      endif
      
      !Initialize variables
      allocate(MH_str(iwse)%wsebnd(MH_str(iwse)%ncells))
      allocate(MH_str(iwse)%wsebnd0(MH_str(iwse)%ncells))
      MH_str(iwse)%wsebnd(:)  = 0.0
      MH_str(iwse)%wsebnd0(:) = 0.0 !May be overwritten in hot_read if initial condition is specified
    enddo !i-str
    
!--- Multiple Water Level and Velocity BC (Type 5=MHV) -----------------------------------
    do iwse=1,nMHVstr
      !Read cell/node string  
      call read_bndstr(MHV_str(iwse)%bidfile,MHV_str(iwse)%bidpath,&
        MHV_str(iwse)%istrtype,MHV_str(iwse)%idnum,&
        MHV_str(iwse)%ncells,MHV_str(iwse)%cells,MHV_str(iwse)%faces)
      
      !Read water level and velocity data
      call read_multiwsedata(MHV_str(iwse)%wsefile,MHV_str(iwse)%wsepath,&
             MHV_str(iwse)%ncells,MHV_str(iwse)%cells,&
             MHV_str(iwse)%ntimeswse,MHV_str(iwse)%timeswse,MHV_str(iwse)%wsedata)
      call read_multiveldata(MHV_str(iwse)%velfile,MHV_str(iwse)%velpath,&
             MHV_str(iwse)%ncells,MHV_str(iwse)%cells,&
             MHV_str(iwse)%ntimesvel,MHV_str(iwse)%timesvel,&
             MHV_str(iwse)%udata,MHV_str(iwse)%vdata)
      
      !Spatial smoothing by moving average
      if(MHV_str(iwse)%nssiwse>0 .and. MHV_str(iwse)%nsswwse>2)then
        allocate(valtemp(MHV_str(iwse)%ntimeswse,1))
        do k=1,MHV_str(iwse)%ntimeswse
          valtemp(:,1) = MHV_str(iwse)%wsedata(k,:)
          call moving_average(MHV_str(iwse)%nssiwse,MHV_str(iwse)%nsswwse,&
                 MHV_str(iwse)%ncells,valtemp(:,1))
          MHV_str(iwse)%wsedata(k,:) = valtemp(:,1)
        enddo
        deallocate(valtemp)
      endif  
      if(MHV_str(iwse)%nssivel>0 .and. MHV_str(iwse)%nsswvel>2)then  
        allocate(valtemp(MHV_str(iwse)%ntimesvel,2))
        do k=1,MHV_str(iwse)%ntimesvel
          valtemp(:,1) = MHV_str(iwse)%udata(k,:)
          valtemp(:,2) = MHV_str(iwse)%vdata(k,:)
          call moving_average(MHV_str(iwse)%nssivel,MHV_str(iwse)%nsswvel,&
                 MHV_str(iwse)%ncells,valtemp(:,1))  
          call moving_average(MHV_str(iwse)%nssivel,MHV_str(iwse)%nsswvel,&
                 MHV_str(iwse)%ncells,valtemp(:,2))  
          MHV_str(iwse)%udata(k,:)   = valtemp(:,1)
          MHV_str(iwse)%vdata(k,:)   = valtemp(:,2)
        enddo
        deallocate(valtemp)
      endif
      
      !Temporal smoothing by moving average
      if(MHV_str(iwse)%nsiwse>0 .and. MHV_str(iwse)%nswwse>0)then
        do j=1,MHV_str(iwse)%ncells
          call moving_average(MHV_str(iwse)%nsiwse,MHV_str(iwse)%nswwse,&
               MHV_str(iwse)%ntimeswse,MHV_str(iwse)%wsedata(:,j))
        enddo
      endif
      if(MHV_str(iwse)%nsivel>0 .and. MHV_str(iwse)%nswvel>0)then
        do j=1,MHV_str(iwse)%ncells
          call moving_average(MHV_str(iwse)%nsivel,MHV_str(iwse)%nswvel,&
               MHV_str(iwse)%ntimesvel,MHV_str(iwse)%udata(:,j))
          call moving_average(MHV_str(iwse)%nsivel,MHV_str(iwse)%nswvel,&
               MHV_str(iwse)%ntimesvel,MHV_str(iwse)%vdata(:,j))
        enddo
      endif
      
      !Initialize variables
      allocate(MHV_str(iwse)%wsebnd(MHV_str(iwse)%ncells))
      allocate(MHV_str(iwse)%wsebnd0(MHV_str(iwse)%ncells))
      allocate(MHV_str(iwse)%ubnd(MHV_str(iwse)%ncells))
      allocate(MHV_str(iwse)%ubnd0(MHV_str(iwse)%ncells))
      allocate(MHV_str(iwse)%vbnd(MHV_str(iwse)%ncells))
      allocate(MHV_str(iwse)%vbnd0(MHV_str(iwse)%ncells))
      MHV_str(iwse)%wsebnd(:)  = 0.0
      MHV_str(iwse)%wsebnd0(:) = 0.0
      MHV_str(iwse)%ubnd(:)    = 0.0
      MHV_str(iwse)%ubnd0(:)   = 0.0
      MHV_str(iwse)%vbnd(:)    = 0.0
      MHV_str(iwse)%vbnd0(:)   = 0.0                         
    enddo !i-str

!--- Cross-shore BC (Type 6=CS) -------------------------------------------------
    !if(nCSstr>0 .and. .not.cmswave)then
    !  call diag_print_error('CMS must be run in steering to use cross-shore boundary conditions)
    !endif
    if(ncellpoly>0 .and. nCSstr>0)then
      call diag_print_error('Cross-shore boundary condition not supported')
    endif
    ibndtemp = 0 !Temporary array for reordering CS_str cells
         
    do icsh=1,nCSstr    
      !Read cell/node string  
      call read_bndstr(CS_str(icsh)%bidfile,CS_str(icsh)%bidpath,&
        CS_str(icsh)%istrtype,CS_str(icsh)%idnum,&
        CS_str(icsh)%ncells,CS_str(icsh)%cells,CS_str(icsh)%faces)
                   
      !Use depths to determine whether to flip cell string 
      !so that the first cell is offshore
      nbndcells=CS_str(icsh)%ncells
      id1=CS_str(icsh)%cells(1)
      id2=CS_str(icsh)%cells(nbndcells)      
      if(zb(id1)>zb(id2))then !Flip cell string
        ibndtemp(1:nbndcells)=CS_str(icsh)%cells(1:nbndcells)
        CS_str(icsh)%cells(1:nbndcells)=ibndtemp(nbndcells:1:-1)
        ibndtemp(1:nbndcells)=CS_str(icsh)%faces(1:nbndcells)
        CS_str(icsh)%faces(1:nbndcells)=ibndtemp(nbndcells:1:-1)
      endif
      
      !Initialize variables
      allocate(CS_str(icsh)%wsecsh(CS_str(icsh)%ncells))
      allocate(CS_str(icsh)%ucsh(CS_str(icsh)%ncells))
      allocate(CS_str(icsh)%vcsh(CS_str(icsh)%ncells))
      CS_str(icsh)%wsecsh(:) = 0.0
      CS_str(icsh)%ucsh(:)   = 0.0
      CS_str(icsh)%vcsh(:)   = 0.0             
    enddo !icsh
    
!--- Parent Simulations ---------------------------------------------    
    do ipar=1,nParSim
      !Read parent grid and control file
      call fileext(ParSim(ipar)%ctlfilepar,aext)
      if(aext(1:7)=='cmcards')then !CMS Control File
        ParSim(ipar)%ipartype = 1 !CMS
        call read_parent_cmcards(ParSim(ipar)%ctlfilepar,tjuldaypar,&
          ParSim(ipar)%grdfilepar,ParSim(ipar)%projpar,ParSim(ipar)%typespathpar,&
          ParSim(ipar)%xoriginpar,ParSim(ipar)%yoriginpar,ParSim(ipar)%orientpar,ParSim(ipar)%telpargrd,&
          ParSim(ipar)%wsefilepar,ParSim(ipar)%wsepathpar,ParSim(ipar)%velfilepar,ParSim(ipar)%velpathpar)
        if(ParSim(ipar)%telpargrd)then !Telescoping grid
          ParSim(ipar)%nmaxfacespar = 6  
          call read_parent_grid_tel(ParSim(ipar)%grdfilepar,ParSim(ipar)%ncellspar,ParSim(ipar)%nptspar,&
            ParSim(ipar)%xpar,ParSim(ipar)%ypar,ParSim(ipar)%dxpar,ParSim(ipar)%dypar,&
            ParSim(ipar)%c2cpar,ParSim(ipar)%idfpar,ParSim(ipar)%ncfpar,ParSim(ipar)%activepar)
        else  !Cartesian grid
          ParSim(ipar)%nmaxfacespar = 4
#ifdef XMDF_IO
          call read_parent_grid_xmdf(ParSim(ipar)%grdfilepar,ParSim(ipar)%typespathpar,&
            ParSim(ipar)%ncellspar,ParSim(ipar)%nptspar,&
            ParSim(ipar)%xpar,ParSim(ipar)%ypar,ParSim(ipar)%dxpar,ParSim(ipar)%dypar,&
            ParSim(ipar)%c2cpar,ParSim(ipar)%idfpar,ParSim(ipar)%ncfpar,ParSim(ipar)%activepar)
#else
          call diag_print_error('Cannot read grid from *.h5 file without XMDF libraries')
#endif
        endif     
        !Reproject if necessary to child grid horizontal projection
        call proj_horiz_conv(ParSim(ipar)%projpar,projfl,ParSim(ipar)%nptspar,&
          ParSim(ipar)%xpar,ParSim(ipar)%ypar)
        ParSim(ipar)%t2hrs = 1.0 !hours to hours
        if(abs(ParSim(ipar)%tjuldaypar+999.0)<1.0e-4)then
          ParSim(ipar)%tjuldaypar = tjuldaypar   !Note: only overwrite if not specified
        endif  
        if(abs(ParSim(ipar)%timestarthr+999.0)<1.0e-4)then
          ParSim(ipar)%timestarthr = (ParSim(ipar)%tjuldaypar-tjulday0)*24.0 !Note: only overwrite if not specified
        endif
      else
        ParSim(ipar)%ipartype = 2 !ADCIRC  
        ParSim(ipar)%nmaxfacespar = 3
        call read_grid14(ParSim(ipar)%grdfilepar,ParSim(ipar)%nelemsfullpar,ParSim(ipar)%nptspar,&
          ParSim(ipar)%xpar,ParSim(ipar)%ypar,ParSim(ipar)%zpar,ParSim(ipar)%elem2node)    
        deallocate(ParSim(ipar)%zpar) !Not needed    
        !Check horizontal projection
        if(ParSim(ipar)%projpar%iHorizDatum==2 .or. ParSim(ipar)%projpar%iHorizCoordSystem==2)then !Local
          call diag_print_warning('No horizontal projection defined for ADCIRC grid',&
            '   Assuming Geographic, NAD83, Degrees')
          ParSim(ipar)%projpar%iHorizDatum = 1         !Horizontal Datum = NAD83
          ParSim(ipar)%projpar%iHorizCoordSystem = 0   !Horizontal Coordinate System = GEOGRAPHIC
          ParSim(ipar)%projpar%iHorizUnits = 4         !Horizontal Units = DEGREES
          ParSim(ipar)%projpar%iHorizZone = 0          !Horizontal Zone = NONE          !Added MEB 04/08/2021
          ParSim(ipar)%projpar%iVertDatum = 9          !Vertical Datum = LOCAL
          ParSim(ipar)%projpar%iVertUnits = 2          !Vertical units = METERS
          ParSim(ipar)%projpar%VertOffset = 0.0        !Vertical offset from datum
        endif    
        !Reproject parent grid to child grid projection
        call proj_horiz_conv(ParSim(ipar)%projpar,projfl,ParSim(ipar)%nptspar,&
          ParSim(ipar)%xpar,ParSim(ipar)%ypar)   
        ParSim(ipar)%t2hrs = 1.0/3600.0 !seconds to hours
      endif
      !Allocate variables
      allocate(ParSim(ipar)%timewsehrspar(ParSim(ipar)%ntiwsepar+1))
      ParSim(ipar)%timewsehrspar(:) = 0.0
      allocate(ParSim(ipar)%wsepar(ParSim(ipar)%nptspar,ParSim(ipar)%ntiwsepar+1))
      ParSim(ipar)%wsepar = 0.0
      if(ParSim(ipar)%velpar)then
        allocate(ParSim(ipar)%timevelhrspar(ParSim(ipar)%ntivelpar+1))  
        ParSim(ipar)%timevelhrspar(:) = 0.0
        allocate(ParSim(ipar)%upar(ParSim(ipar)%nptspar,ParSim(ipar)%ntivelpar+1))
        allocate(ParSim(ipar)%vpar(ParSim(ipar)%nptspar,ParSim(ipar)%ntivelpar+1))
        ParSim(ipar)%upar = 0.0
        ParSim(ipar)%vpar = 0.0
      endif
    enddo !ipar
    
!--- Nested Water Level BC (Type 7=NH) ---------------------------------------------
    do iwse=1,nNHstr
      !Read cell/node string  
      call read_bndstr(NH_str(iwse)%bidfile,NH_str(iwse)%bidpath,&
        NH_str(iwse)%istrtype,NH_str(iwse)%idnum,&
        NH_str(iwse)%ncells,NH_str(iwse)%cells,NH_str(iwse)%faces)  
        
      !Initialize variables
      allocate(NH_str(iwse)%xbnd(NH_str(iwse)%ncells))
      allocate(NH_str(iwse)%ybnd(NH_str(iwse)%ncells))
      NH_str(iwse)%xbnd(:) = xc(NH_str(iwse)%cells(:)) !Global x-coordinate
      NH_str(iwse)%ybnd(:) = yc(NH_str(iwse)%cells(:)) !Global y-coordinate
      allocate(NH_str(iwse)%timewsehrs(NH_str(iwse)%ntiwse+1))
      NH_str(iwse)%timewsehrs(:) = 0.0
      allocate(NH_str(iwse)%wsebnd (NH_str(iwse)%ncells))
      allocate(NH_str(iwse)%wsebnd0(NH_str(iwse)%ncells))
      allocate(NH_str(iwse)%wsedata(NH_str(iwse)%ncells,NH_str(iwse)%ntiwse+1))
      NH_str(iwse)%wsebnd(:)  = 0.0
      NH_str(iwse)%wsebnd0(:) = 0.0
      NH_str(iwse)%wsedata(:,:) = 0.0
     
      !Interpolation coefficients
      idpar = NH_str(iwse)%idpar
      if(ParSim(idpar)%ipartype==1)then !CMS
        mntp = 4  
      else !if(ParSim(idpar)%ipartype==2)then !ADCIRC
        mntp = 3  
      endif      
      NH_str(iwse)%mntp = mntp
      !Calculate interpolation coefficients    
      xtrapdist = 1.0e6 !Set large value to use nearest neighbor value      
      nbndcells = NH_str(iwse)%ncells
      allocate(NH_str(iwse)%intp(0:mntp,nbndcells),NH_str(iwse)%cntp(mntp,nbndcells))
      NH_str(iwse)%intp(:,:) = 0; NH_str(iwse)%cntp(:,:) = 0.0
      if(ParSim(idpar)%ipartype==1)then !CMS
        call interp_coef_tel2pts(ParSim(idpar)%ncellspar,&
          ParSim(idpar)%nptspar,ParSim(idpar)%nmaxfacespar,&
          ParSim(idpar)%xOriginpar,ParSim(idpar)%yOriginpar,ParSim(idpar)%orientpar, &
          ParSim(idpar)%xpar,ParSim(idpar)%ypar,ParSim(idpar)%dxpar,ParSim(idpar)%dypar,&
          ParSim(idpar)%c2cpar,ParSim(idpar)%idfpar,ParSim(idpar)%ncfpar, &            
          NH_str(iwse)%ncells,NH_str(iwse)%xbnd,NH_str(iwse)%ybnd,&
          xtrapdist,NH_str(iwse)%intp,NH_str(iwse)%cntp)     
      else
        call interp_coef_tri2pts(ParSim(idpar)%nelemsfullpar,ParSim(idpar)%nptspar,&
          ParSim(idpar)%xpar,ParSim(idpar)%ypar,ParSim(idpar)%elem2node, &            
          NH_str(iwse)%ncells,NH_str(iwse)%xbnd,NH_str(iwse)%ybnd,&
          xtrapdist,NH_str(iwse)%intp,NH_str(iwse)%cntp)  
      endif
    enddo !iwse-str
    
!--- Nested Water Level and Velocity BC (Type 8=NHV) -----------------------------------
    do iwse=1,nNHVstr
      !Read cell/node string  
      call read_bndstr(NHV_str(iwse)%bidfile,NHV_str(iwse)%bidpath,&
        NHV_str(iwse)%istrtype,NHV_str(iwse)%idnum,&
        NHV_str(iwse)%ncells,NHV_str(iwse)%cells,NHV_str(iwse)%faces)
           
      !Initialize variables
      allocate(NHV_str(iwse)%xbnd(NHV_str(iwse)%ncells))
      allocate(NHV_str(iwse)%ybnd(NHV_str(iwse)%ncells))
      NHV_str(iwse)%xbnd(:) = xc(NHV_str(iwse)%cells(:)) !Global x-coordinate
      NHV_str(iwse)%ybnd(:) = yc(NHV_str(iwse)%cells(:)) !Global y-coordinate
      allocate(NHV_str(iwse)%timewsehrs(NHV_str(iwse)%ntiwse+1))
      allocate(NHV_str(iwse)%timevelhrs(NHV_str(iwse)%ntiwse+1))
      NHV_str(iwse)%timewsehrs(:) = 0.0
      NHV_str(iwse)%timevelhrs(:) = 0.0
      allocate(NHV_str(iwse)%wsebnd(NHV_str(iwse)%ncells))
      allocate(NHV_str(iwse)%wsebnd0(NHV_str(iwse)%ncells))
      allocate(NHV_str(iwse)%wsedata(NHV_str(iwse)%ncells,NHV_str(iwse)%ntiwse+1))
      NHV_str(iwse)%wsebnd(:)  = 0.0
      NHV_str(iwse)%wsebnd0(:) = 0.0
      NHV_str(iwse)%wsedata(:,:) = 0.0
      allocate(NHV_str(iwse)%ubnd(NHV_str(iwse)%ncells))
      allocate(NHV_str(iwse)%ubnd0(NHV_str(iwse)%ncells))
      allocate(NHV_str(iwse)%udata(NHV_str(iwse)%ncells,NHV_str(iwse)%ntivel+1))
      NHV_str(iwse)%ubnd(:)  = 0.0
      NHV_str(iwse)%ubnd0(:) = 0.0
      NHV_str(iwse)%udata(:,:) = 0.0
      allocate(NHV_str(iwse)%vbnd(NHV_str(iwse)%ncells))
      allocate(NHV_str(iwse)%vbnd0(NHV_str(iwse)%ncells))
      allocate(NHV_str(iwse)%vdata(NHV_str(iwse)%ncells,NHV_str(iwse)%ntivel+1))
      NHV_str(iwse)%vbnd(:)  = 0.0
      NHV_str(iwse)%vbnd0(:) = 0.0
      NHV_str(iwse)%vdata(:,:) = 0.0
      
      !Interpolation coefficients
      idpar = NHV_str(iwse)%idpar
      if(ParSim(idpar)%ipartype==1)then !CMS
        mntp = 4  
      else !if(ParSim(idpar)%ipartype==2)then !ADCIRC
        mntp = 3  
      endif      
      NHV_str(iwse)%mntp = mntp
      !Calculate interpolation coefficients    
      xtrapdist = 1.0e6 !Set large value to use nearest neighbor value      
      nbndcells = NHV_str(iwse)%ncells
      allocate(NHV_str(iwse)%intp(0:mntp,nbndcells),NHV_str(iwse)%cntp(mntp,nbndcells))
      NHV_str(iwse)%intp(:,:) = 0; NHV_str(iwse)%cntp(:,:) = 0.0
      if(ParSim(idpar)%ipartype==1)then !CMS
        call interp_coef_tel2pts(ParSim(idpar)%ncellspar,&
          ParSim(idpar)%nptspar,ParSim(idpar)%nmaxfacespar,&
          ParSim(idpar)%xOriginpar,ParSim(idpar)%yOriginpar,ParSim(idpar)%orientpar, &
          ParSim(idpar)%xpar,ParSim(idpar)%ypar,ParSim(idpar)%dxpar,ParSim(idpar)%dypar,&
          ParSim(idpar)%c2cpar,ParSim(idpar)%idfpar,ParSim(idpar)%ncfpar, &            
          NHV_str(iwse)%ncells,NHV_str(iwse)%xbnd,NHV_str(iwse)%ybnd,&
          xtrapdist,NHV_str(iwse)%intp,NHV_str(iwse)%cntp)     
      else
        call interp_coef_tri2pts(ParSim(idpar)%nelemsfullpar,ParSim(idpar)%nptspar,&
          ParSim(idpar)%xpar,ParSim(idpar)%ypar,ParSim(idpar)%elem2node, &            
          NHV_str(iwse)%ncells,NHV_str(iwse)%xbnd,NHV_str(iwse)%ybnd,&
          xtrapdist,NHV_str(iwse)%intp,NHV_str(iwse)%cntp)  
      endif
    enddo !iwse-str
    
!--- Tidal Database WSE BC (Type 9=NTH) -----------------------------------
    do iwse=1,nNTHstr
      !Read cell/node string  
      call read_bndstr(NTH_str(iwse)%bidfile,NTH_str(iwse)%bidpath,&
        NTH_str(iwse)%istrtype,NTH_str(iwse)%idnum,&
        NTH_str(iwse)%ncells,NTH_str(iwse)%cells,NTH_str(iwse)%faces)
      
      !Initialize variables
      allocate(NTH_str(iwse)%xbnd(NTH_str(iwse)%ncells))
      allocate(NTH_str(iwse)%ybnd(NTH_str(iwse)%ncells))
      NTH_str(iwse)%xbnd(:) = xc(NTH_str(iwse)%cells(:)) !Global x-coordinate
      NTH_str(iwse)%ybnd(:) = yc(NTH_str(iwse)%cells(:)) !Global y-coordinate
      allocate(NTH_str(iwse)%wsebnd (NTH_str(iwse)%ncells))
      allocate(NTH_str(iwse)%wsebnd0(NTH_str(iwse)%ncells))
      allocate(NTH_str(iwse)%wseadj(NTH_str(iwse)%ncells))
      NTH_str(iwse)%wsebnd(:)  = 0.0
      NTH_str(iwse)%wsebnd0(:) = 0.0     
      NTH_str(iwse)%wseadj(:) = 0.0
      
      !Initialize constituent variables
      if(NTH_str(iwse)%ntc==0 .and. NTH_str(iwse)%ntcin==0) allocate(NTH_str(iwse)%namein(1))
      call tdb_init(NTH_str(iwse)%tdbname,NTH_str(iwse)%tdbpath,&
         NTH_str(iwse)%projtdb,NTH_str(iwse)%ncells,&
         NTH_str(iwse)%xbnd,NTH_str(iwse)%ybnd,&
         NTH_str(iwse)%nssi,NTH_str(iwse)%nssw,&
         NTH_str(iwse)%ntcin,NTH_str(iwse)%namein,&
         NTH_str(iwse)%ntc,NTH_str(iwse)%name,&
         NTH_str(iwse)%speed,NTH_str(iwse)%f,NTH_str(iwse)%vu,&
         NTH_str(iwse)%amp,NTH_str(iwse)%phase)
    enddo !iwse-str
    
!--- Tidal Database WSE and Velocity BC (Type 10=NTHV) --------------------------------
    do iwse=1,nNTHVstr
      !Read cell/node string  
      call read_bndstr(NTHV_str(iwse)%bidfile,NTHV_str(iwse)%bidpath,&
        NTHV_str(iwse)%istrtype,NTHV_str(iwse)%idnum,&
        NTHV_str(iwse)%ncells,NTHV_str(iwse)%cells,NTHV_str(iwse)%faces)
      
      !Initialize variables
      allocate(NTHV_str(iwse)%xbnd(NTHV_str(iwse)%ncells))
      allocate(NTHV_str(iwse)%ybnd(NTHV_str(iwse)%ncells))
      NTHV_str(iwse)%xbnd(:) = xc(NTHV_str(iwse)%cells(:)) !Global x-coordinate
      NTHV_str(iwse)%ybnd(:) = yc(NTHV_str(iwse)%cells(:)) !Global y-coordinate
      allocate(NTHV_str(iwse)%wsebnd (NTHV_str(iwse)%ncells))
      allocate(NTHV_str(iwse)%wsebnd0(NTHV_str(iwse)%ncells))
      allocate(NTHV_str(iwse)%wseadj(NTHV_str(iwse)%ncells))
      NTHV_str(iwse)%wsebnd(:)  = 0.0
      NTHV_str(iwse)%wsebnd0(:) = 0.0
      NTHV_str(iwse)%wseadj(:) = 0.0
      allocate(NTHV_str(iwse)%ubnd(NTHV_str(iwse)%ncells))
      allocate(NTHV_str(iwse)%ubnd0(NTHV_str(iwse)%ncells))
      NTHV_str(iwse)%ubnd(:)  = 0.0
      NTHV_str(iwse)%ubnd0(:) = 0.0
      allocate(NTHV_str(iwse)%vbnd(NTHV_str(iwse)%ncells))
      allocate(NTHV_str(iwse)%vbnd0(NTHV_str(iwse)%ncells))
      NTHV_str(iwse)%vbnd(:)  = 0.0
      NTHV_str(iwse)%vbnd0(:) = 0.0      
      
      !Initialize constituent variables
      if(NTHV_str(iwse)%ntc==0) allocate(NTHV_str(iwse)%namein(1))
      if(NTHV_str(iwse)%tdbname(1:9)=='LEPROVOST' .or. &
         NTHV_str(iwse)%tdbname(1:3)=='FES')then
        call diag_print_error('Invalid Boundary Specified',&
          '  Cannot use the LEPROVOST, FES952 or FES2004 tidal databases',&
          '  These databases do not contain velocities.',&
          '  Either change the boundary condition or select',&
          '  a the tidal database which contains current velocities',&
          '  such as EC2001 or ENPAC2003')
      endif
      call tdb_init(NTHV_str(iwse)%tdbname,NTHV_str(iwse)%tdbpath,&
         NTHV_str(iwse)%projtdb,NTHV_str(iwse)%ncells,&
         NTHV_str(iwse)%xbnd,NTHV_str(iwse)%ybnd,&
         NTHV_str(iwse)%nssi,NTHV_str(iwse)%nssw,&
         NTHV_str(iwse)%ntcin,NTHV_str(iwse)%namein,&
         NTHV_str(iwse)%ntc,NTHV_str(iwse)%name,&        
         NTHV_str(iwse)%speed,NTHV_str(iwse)%f,NTHV_str(iwse)%vu,&
         NTHV_str(iwse)%amp,NTHV_str(iwse)%phase,&
         NTHV_str(iwse)%ampu,NTHV_str(iwse)%phaseu,&
         NTHV_str(iwse)%ampv,NTHV_str(iwse)%phasev,NTHV_str(iwse)%angvel)        
    enddo !iwse-str

!!--- River BC - extrapolation ------------------       
!    allocate(iextrap(ncellsD))
!    iextrap = 0 !Initialize    
!    do iriv=1,nQstr
!      do j=1,Q_str(iriv)%ncells
!        i=Q_str(iriv)%cells(j)
!        k=Q_str(iriv)%faces(j)
!        iextrap(cell2cell(k,i))=1
!      enddo
!    enddo
    
!    call write_bnd_check

!--- Deallocate no longer used Parent Simulation Variables --------------------
    do ipar=1,nParSim
      if(ParSim(ipar)%ipartype==1)then
        deallocate(ParSim(ipar)%c2cpar,ParSim(ipar)%idfpar,ParSim(ipar)%ncfpar)
        deallocate(ParSim(ipar)%dxpar,ParSim(ipar)%dypar,ParSim(ipar)%activepar)
      else 
        deallocate(ParSim(ipar)%elem2node)
      endif
    enddo
    
!---- Store all boundary strings in one structure --------
    nbndstr = nQstr + nTHstr + nHstr + nMHstr + nMHVstr + nCSstr + nNHstr + nNHVstr + nNTHstr + nNTHVstr
    allocate(bnd_str(nbndstr))
    
    ibnd = 0 !Counter for all cell strings
    
    !River boundary
    do kk=1,nQstr
      call copy2bnd(ibnd,1,Q_str(kk)%ncells,Q_str(kk)%cells,Q_str(kk)%faces,Q_str(kk)%idnum)
    enddo 

    !Tidal/Harmonic boundary
    do kk=1,nTHstr  
      call copy2bnd(ibnd,2,TH_str(kk)%ncells,TH_str(kk)%cells,TH_str(kk)%faces,TH_str(kk)%idnum)
    enddo
    
    !Single Water Level BC  
    do kk=1,nHstr
      call copy2bnd(ibnd,3,H_str(kk)%ncells,H_str(kk)%cells,H_str(kk)%faces,H_str(kk)%idnum)
    enddo      

    !Multiple Water Surface Elevation BC  
    do kk=1,nMHstr
      call copy2bnd(ibnd,4,MH_str(kk)%ncells,MH_str(kk)%cells,MH_str(kk)%faces,MH_str(kk)%idnum)
    enddo

    !Multiple Water Surface Elevation and Velocity BC 
    do kk=1,nMHVstr
      call copy2bnd(ibnd,5,MHV_str(kk)%ncells,MHV_str(kk)%cells,MHV_str(kk)%faces,MHV_str(kk)%idnum)
    enddo
    
    !Cross-shore BC 
    do kk=1,nCSstr
      call copy2bnd(ibnd,6,CS_str(kk)%ncells,CS_str(kk)%cells,CS_str(kk)%faces, CS_str(kk)%idnum)
    enddo
    
    !Nested Water Surface Elevation BC  
    do kk=1,nNHstr
      call copy2bnd(ibnd,7,NH_str(kk)%ncells,NH_str(kk)%cells,NH_str(kk)%faces, NH_str(kk)%idnum)
    enddo
    
    !Nested Water Surface Elevation and Velocity BC 
    do kk=1,nNHVstr
      call copy2bnd(ibnd,8,NHV_str(kk)%ncells,NHV_str(kk)%cells,NHV_str(kk)%faces, NHV_str(kk)%idnum)
    enddo
    
    !Nested Tidal Database WSE Boundary
    do kk=1,nNTHstr
      call copy2bnd(ibnd,9,NTH_str(kk)%ncells,NTH_str(kk)%cells,NTH_str(kk)%faces, NTH_str(kk)%idnum)
    enddo
    
    !Nested Tidal Database WSE and Velocity Boundary
    do kk=1,nNTHVstr
      call copy2bnd(ibnd,10,NTHV_str(kk)%ncells,NTHV_str(kk)%cells,NTHV_str(kk)%faces, NTHV_str(kk)%idnum)
    enddo
    
    !Wall boundary
    ibndtemp = 1 !Initialize
    do ibnd=1,nbndstr  
      do j=1,bnd_str(ibnd)%ncells
        i = bnd_str(ibnd)%cells(j)
        k = bnd_str(ibnd)%faces(j)
        nck = cell2cell(k,i)
        if(nck==0)then
          if(allocated(mapid))then
            write(msg2,*) ' Problem at boundary cell: ',mapid(i)
          else
            write(msg2,*) ' Problem at boundary cell: ',i
          endif
          call diag_print_error(msg2)
        endif
        ibndtemp(nck) = 0
      enddo
    enddo
    !Count wall faces
    W_str%ncells = 0
    do i=1,ncells
      do k=1,ncface(i)
        nck = cell2cell(k,i)  
        if(nck>ncells .and. ibndtemp(nck)==1)then !Boundary with no assigned condition
          W_str%ncells = W_str%ncells+1 
        endif  
      enddo
    enddo
    !Initialize wall variables
    allocate(W_str%cells(W_str%ncells),W_str%faces(W_str%ncells))
    W_str%cells = 0; W_str%faces = 0    
    !Calculate wall cell ID's and faces
    j=0
    do i=1,ncells
      do k=1,ncface(i)
        nck = cell2cell(k,i)  
        if(nck>ncells .and. ibndtemp(nck)==1)then !Boundary with no assigned condition
          j=j+1  
          W_str%cells(j)=i
          W_str%faces(j)=k
        endif  
      enddo
    enddo
    
    return
    contains
!******************************************************************************    
    subroutine copy2bnd(ibnd,ibndtype,nbndcells,ibndcells,jbndfaces,idnum)
!******************************************************************************    
    use bnd_def, only: bnd_str
    
    implicit none
    integer,intent(inout) :: ibnd
    integer,intent(in) :: ibndtype,nbndcells,idnum
    integer,intent(in),dimension(nbndcells) :: ibndcells,jbndfaces
    
    ibnd = ibnd + 1
    bnd_str(ibnd)%ibndtype = ibndtype
    bnd_str(ibnd)%ncells = nbndcells
    allocate(bnd_str(ibnd)%cells(nbndcells))
    allocate(bnd_str(ibnd)%faces(nbndcells))
    bnd_str(ibnd)%cells(:) = ibndcells(:)
    bnd_str(ibnd)%faces(:) = jbndfaces(:)
    bnd_str(ibnd)%idnum = idnum
    
    return
    end subroutine copy2bnd
    
    end subroutine bnd_init
        
!*****************************************************************************************
    subroutine bnd_print()
! Prints the boundary condition setup to the screen and diagnostic file
! written by Alex Sanchez, USACE-CHL
!*****************************************************************************************    
    use bnd_def
    use geo_def,   only: azimuth_fl,mapid
    use struct_def
    use geo_def,   only: aHorizCoordSystem,aHorizDatum,aHorizUnits
    use diag_def,  only: dgunit, dgfile
    use time_lib,  only: julian2calendar
    use const_def, only: rad2deg,deg2rad
    use out_def,   only: write_ascii_input,outprefix
    use prec_def,  only: ikind
    use met_def,   only: wndspeed,wnddirection
    use tool_def,  only: vstrlz
    
    implicit none
    integer :: iyrpar,imopar,idaypar,ihrpar,iminpar,isecpar,ibnd
    integer :: i,j,k,iunit(2),iriv,iwse,icsh,ntc,nbndcells,idpar,kunit
    integer,allocatable :: ibndcells(:)
    integer :: nstrcells
    real(ikind) :: ampavg,phaseavg,ampmax,phasemax,ampmin,phasemin
    real(ikind) :: cosang,sinang,ang,valx,valy,val
    character(len=200) :: apath,aname,astring,bidoutfile
    character(len=10) :: aext

141 format(' ',A,T40,A)
241 format(' ',A,T40,I0)
261 format(' ',A,T40,I0)
262 format(' ',A,T40,I4.4)
341 format(' ',A,T40,A,A)     !Added for vstrlz function results
342 format(' ',A,T40,F0.2,A)
353 format(' ',A,T40,F0.3,A)
374 format(' ',A,T40,E10.3)
784 format(' ',A,T40,2(1x,E9.2))
450 format('       ID  Name   Speed     Amplitude  Phase   Nodal Factor  Eq. Arg.')
451 format('                 [deg/hr]      [m]     [deg]       [-]        [deg]')   
440 format('       ID    Speed     Amplitude  Phase')
441 format('            [deg/hr]      [m]     [deg]')     
452 format(' ',6x,I2,2x,A6,F9.5,2x,F7.3,2x,F8.2,3x,F7.3,4x,F8.2)
442 format(' ',6x,I2,2x,F9.5,2x,F7.3,2x,F8.2)
453 format(' ',6x,I2,2x,A6,F9.5,1x,3(F6.3),1x,3(F6.1),1x,F6.3,2x,F6.2)
541 format('       Summary of Water Level Constituents:')
531 format('       Summary of U-Velocity Constituents:')
532 format('       Summary of V-Velocity Constituents:')
542 format('                              Amplitude            Phase         Nodal     Eq.')
543 format('                   Speed         [m]               [deg]         Factor   Arg.')
533 format('                   Speed        [m/s]              [deg]         Factor   Arg.')
544 format('       ID  Name   [deg/hr]  Min   Max   Avg    Min   Max   Avg    [-]    [deg]')
431 format('          Start Date and Time:',T40,I4,'-',I2.2,'-',I2.2,' ',I2.2,':',I2.2,':',I2.2,' UTC')
    
    iunit = (/6, dgunit/)
    open(dgunit,file=dgfile,access='append') 
    do i=1,2
      write(iunit(i),*)  
      write(iunit(i),141)     '  Boundaries'
      if(veldamp>1.0e-6)then
        write(iunit(i),374)   '    Open WSE Boundary Velocity Damping: ',veldamp
      endif
      !--- River/Flux BC (Type 1=Q) -------------------------------      
      do iriv=1,nQstr
        write(iunit(i),241)   '    Flux Boundary:',iriv
        write(iunit(i),141)   '      Cellstring File:',trim(Q_str(iriv)%bidfile)        
        write(iunit(i),141)   '      Cellstring Path:',trim(Q_str(iriv)%bidpath)
        write(iunit(i),241)   '      Boundary Cells:',Q_str(iriv)%ncells
        if(Q_str(iriv)%ifluxmode==2)then
          write(iunit(i),141) '      Flux Time Series:'
          write(iunit(i),141) '        Data File:',trim(Q_str(iriv)%fluxfile)
          write(iunit(i),141) '        Data Path:',trim(Q_str(iriv)%fluxpath)
          write(iunit(i),261) '        Number of Points:',Q_str(iriv)%ntimes
        elseif(Q_str(iriv)%ifluxmode==1)then
          !Unit Conversion, ifluxunits = 0-m^3/s/cell, 1-m^3/s/boundary, 2-ft^3/s/boundary, default = 0      
          if(Q_str(iriv)%ifluxunits==0)then
            val = Q_str(iriv)%ncells !from m^3/s/cell to m^3/s/boundary
          elseif(Q_str(iriv)%ifluxunits==2)then
            val = 0.3048**3 !from ft^3/s to m^3/s/boundary 
          else
            val = 1.0  !m^3/s/boundary 
          endif  
          write(iunit(i),374) '      Flux Value:',Q_str(iriv)%qfluxconst/val
        elseif(Q_str(iriv)%ifluxmode==3)then
          write(iunit(i),141) '      Stage-Flux Curve:'
          write(iunit(i),141) '        Data File:',trim(Q_str(iriv)%fluxfile)
          write(iunit(i),141) '        Data Path:',trim(Q_str(iriv)%fluxpath)
          write(iunit(i),261) '        Number of Points:',Q_str(iriv)%ntimes
        endif
        if(Q_str(iriv)%ifluxunits==0)then
          write(iunit(i),141) '      Flux Input Units:','m^3/s/cell'
        elseif(Q_str(iriv)%ifluxunits==1)then
          write(iunit(i),141) '      Flux Input Units:','m^3/s/boundary'  
        else
          write(iunit(i),141) '      Flux Input Units:','ft^3/s/boundary'  
        endif
        ang = 90.0 - Q_str(iriv)%angle - azimuth_fl
        if(ang<0.0)then
          ang = ang + 360.0
        elseif(ang>360.0)then
          ang = ang - 360.0
        endif
        write(iunit(i),341)   '      Direction:',trim(vstrlz(ang,'(f0.3)')),' deg'
        write(iunit(i),341)   '      Conveyance Coefficient:',trim(vstrlz(Q_str(iriv)%cmvel,'(f0.3)'))
      enddo !i-str
    
      !--- Tidal/Harmonic BC (Type 2=TH) ------------------------------------------------
      do iwse=1,nTHstr
        if(TH_str(iwse)%istidal)then
          if(len_trim(TH_str(iwse)%station)>0)then
            write(iunit(i),241) '    Tidal Station Boundary:',iwse  
            write(iunit(i),141) '      Station Name:',trim(TH_str(iwse)%station)
          else
            write(iunit(i),241) '    Tidal Boundary:',iwse
          endif
        else
          write(iunit(i),241)   '    Harmonic Boundary:',iwse  
        endif
        write(iunit(i),141)     '      Cellstring File:',trim(TH_str(iwse)%bidfile)    
        write(iunit(i),141)     '      Cellstring Path:',trim(TH_str(iwse)%bidpath)
        write(iunit(i),261)     '      Boundary Cells:',TH_str(iwse)%ncells
        if (TH_str(iwse)%ioffsetmode.eq.1) then
          write(iunit(i),341)   '      Water/Sea Level Change Offset:',trim(vstrlz(TH_str(iwse)%wseoffset,'(f0.3)')),' m'
        elseif (TH_str(iwse)%ioffsetmode.eq.2)then
          write(iunit(i),141)   '      Water/Sea Level Change Curve File:',trim(TH_str(iwse)%offsetfile)    
        endif
        
        if(abs(TH_str(iwse)%dwsex)>1.0e-9 .or. abs(TH_str(iwse)%dwsey)>1.0e-9)then
          !Rotate gradients to the global coordinate system
          cosang=cos(-azimuth_fl*deg2rad); sinang=sin(-azimuth_fl*deg2rad)
          valx= TH_str(iwse)%dwsex*cosang+TH_str(iwse)%dwsey*sinang
          valy=-TH_str(iwse)%dwsex*sinang+TH_str(iwse)%dwsey*cosang
          if(abs(valx)<1.0e-10) valx=0.0
          if(abs(valy)<1.0e-10) valy=0.0  
          write(iunit(i),784)   '      Water Level Gradients:',TH_str(iwse)%dwsex,TH_str(iwse)%dwsey
        endif 
        if(TH_str(iwse)%angle>-900)then
          ang = 90.0 - TH_str(iwse)%angle*rad2deg - azimuth_fl
          if(ang<0.0)then
            ang = ang + 360.0
          elseif(ang>360.0)then
            ang = ang - 360.0
          endif
          write(iunit(i),341)   '      Incident Direction:',trim(vstrlz(ang,'(f0.3)')),' deg'
        endif               
        if(TH_str(iwse)%istidal)then
          write(iunit(i),261)   '      Number of Constituents:',TH_str(iwse)%ntc
          write(iunit(i),450)
          write(iunit(i),451)
          do k=1,TH_str(iwse)%ntc
            write(iunit(i),452) k,TH_str(iwse)%name(k),TH_str(iwse)%speed(k)*rad2deg,&
               TH_str(iwse)%amp(k),TH_str(iwse)%phase(k)*rad2deg,&
               TH_str(iwse)%f(k),TH_str(iwse)%vu(k)*rad2deg
          enddo
        else
          write(iunit(i),261)   '      Number of Harmonics:',TH_str(iwse)%ntc
          write(iunit(i),440)
          write(iunit(i),441)
          do k=1,TH_str(iwse)%ntc
            write(iunit(i),442) k,TH_str(iwse)%speed(k)*rad2deg,&
               TH_str(iwse)%amp(k),TH_str(iwse)%phase(k)*rad2deg
          enddo
        endif
        if(abs(TH_str(iwse)%dwsex)>1.0e-9 .or. abs(TH_str(iwse)%dwsey)>1.0e-9)then
          !Rotate gradients to the global coordinate system
          cosang=cos(-azimuth_fl*deg2rad); sinang=sin(-azimuth_fl*deg2rad)
          valx= TH_str(iwse)%dwsex*cosang+TH_str(iwse)%dwsey*sinang
          valy=-TH_str(iwse)%dwsex*sinang+TH_str(iwse)%dwsey*cosang
          if(abs(valx)<1.0e-10) valx=0.0
          if(abs(valy)<1.0e-10) valy=0.0
          write(iunit(i),784)   '      Water Level Gradients:',valx,valy
        endif
        if(TH_str(iwse)%wseadjust)then
          write(iunit(i),141)   '      WSE Forcing Correction:','ON'
        else
          write(iunit(i),141)   '      WSE Forcing Correction:','OFF'
        endif
      enddo !i-str   

      !--- Single Water Level BC (Type 3=H) -------------------------------------------
      do iwse=1,nHstr
        write(iunit(i),241)     '    Single Water Level Boundary:',iwse
        write(iunit(i),141)     '      Cellstring File:',trim(H_str(iwse)%bidfile)        
        write(iunit(i),141)     '      Cellstring Path:',trim(H_str(iwse)%bidpath)
        write(iunit(i),261)     '      Boundary Cells:',H_str(iwse)%ncells
        if(len_trim(H_str(iwse)%wsefile)>0)then
          write(iunit(i),141)   '      Water Level Data File:',trim(H_str(iwse)%wsefile)
          if(len_trim(H_str(iwse)%wsepath)>0)then
            write(iunit(i),141) '      Water Level Data Path:',trim(H_str(iwse)%wsepath)
          endif
          write(iunit(i),261)   '      Water Level Data Points:',H_str(iwse)%ntimes  
          write(iunit(i),261)   '      Temporal Interp Order:',H_str(iwse)%nti
        else !Constant
          write(iunit(i),353)   '      Water Level Value:',H_str(iwse)%wseconst,' m'  
        endif  
        if (H_str(iwse)%ioffsetmode.eq.1) then
          write(iunit(i),341)   '      Water/Sea Level Change Offset:',trim(vstrlz(H_str(iwse)%wseoffset,'(f0.3)')),' m'
        elseif (H_str(iwse)%ioffsetmode.eq.2)then
          write(iunit(i),141)   '      Water/Sea Level Change Curve File:',trim(H_str(iwse)%offsetfile)    
        endif
        if(abs(H_str(iwse)%dwsex)>1.0e-9 .or. abs(H_str(iwse)%dwsey)>1.0e-9)then
          !Rotate gradients to the global coordinate system
          cosang=cos(-azimuth_fl*deg2rad); sinang=sin(-azimuth_fl*deg2rad)
          valx= H_str(iwse)%dwsex*cosang+H_str(iwse)%dwsey*sinang
          valy=-H_str(iwse)%dwsex*sinang+H_str(iwse)%dwsey*cosang
          if(abs(valx)<1.0e-10) valx=0.0
          if(abs(valy)<1.0e-10) valy=0.0
          write(iunit(i),784)   '      Water Level Gradients:',valx,valy
        endif
        if(H_str(iwse)%wseadjust)then
          write(iunit(i),141)   '      WSE Forcing Correction:','ON'
        else
          write(iunit(i),141)   '      WSE Forcing Correction:','OFF'
        endif
      enddo !iwse

      !--- Multiple Water Level BC (Type 4=MH) -----------------------------------------------
      do iwse=1,nMHstr
        write(iunit(i),241)     '    Multi WSE Boundary: ',iwse
        write(iunit(i),141)     '      Cellstring File:',trim(MH_str(iwse)%bidfile)        
        write(iunit(i),141)     '      Cellstring Path:',trim(MH_str(iwse)%bidpath)
        write(iunit(i),261)     '      Boundary Cells:',MH_str(iwse)%ncells
        write(iunit(i),141)     '      WSE Data:'
        write(iunit(i),141)     '        File:',trim(MH_str(iwse)%wsefile)
        if(len_trim(MH_str(iwse)%wsepath)>0)then
          write(iunit(i),141)   '        Path:',trim(MH_str(iwse)%wsepath)
        endif
        write(iunit(i),261)     '        Data Points:',MH_str(iwse)%ntimes
        write(iunit(i),141)     '        Temporal Smoothing:'
        write(iunit(i),261)     '          Iterations:',MH_str(iwse)%nsi
        write(iunit(i),261)     '          Width:',MH_str(iwse)%nsw
        write(iunit(i),141)     '        Spatial Smoothing:'
        write(iunit(i),261)     '          Iterations:',MH_str(iwse)%nssi
        write(iunit(i),261)     '          Width:',MH_str(iwse)%nssw
      enddo !i-str
    
      !--- Multiple Water Level and Velocity BC (Type 5=MHV) -----------------------------------
      do iwse=1,nMHVstr
        write(iunit(i),241)     '    Multi WSE and Velocity Boundary: ',iwse
        write(iunit(i),141)     '      Cellstring File:',trim(MHV_str(iwse)%bidfile)        
        write(iunit(i),141)     '      Cellstring Path:',trim(MHV_str(iwse)%bidpath)
        write(iunit(i),261)     '      Boundary Cells:',MHV_str(iwse)%ncells
        write(iunit(i),141)     '      WSE Data:'
        write(iunit(i),141)     '        File:',trim(MHV_str(iwse)%wsefile)
        if(len_trim(MHV_str(iwse)%wsepath)>0)then
          write(iunit(i),141)   '        Path:',trim(MHV_str(iwse)%wsepath)
        endif
        write(iunit(i),261)     '        Time Steps:',MHV_str(iwse)%ntimeswse
        write(iunit(i),141)     '        Temporal Smoothing:'
        write(iunit(i),261)     '          Iterations:',MHV_str(iwse)%nsiwse
        write(iunit(i),261)     '          Width:',MHV_str(iwse)%nswwse
        write(iunit(i),141)     '        Spatial Smoothing:'
        write(iunit(i),261)     '          Iterations:',MHV_str(iwse)%nssiwse
        write(iunit(i),261)     '          Width:',MHV_str(iwse)%nsswwse
        write(iunit(i),141)     '      Velocity Data:'
        write(iunit(i),141)     '        File:',trim(MHV_str(iwse)%velfile)
        if(len_trim(MHV_str(iwse)%velpath)>0)then
          write(iunit(i),141)   '        Path:',trim(MHV_str(iwse)%velpath)
        endif
        write(iunit(i),261)     '        Time Steps:',MHV_str(iwse)%ntimesvel
        write(iunit(i),141)     '        Temporal Smoothing:'
        write(iunit(i),261)     '          Iterations:',MHV_str(iwse)%nsivel
        write(iunit(i),261)     '          Width:',MHV_str(iwse)%nswvel
        write(iunit(i),141)     '        Spatial Smoothing:'
        write(iunit(i),261)     '          Iterations:',MHV_str(iwse)%nssivel
        write(iunit(i),261)     '          Width:',MHV_str(iwse)%nsswvel
      enddo !i-str

      !--- Cross-shore BC (Type 6=CS) -------------------------------------------------
      do icsh=1,nCSstr    
        write(iunit(i),241)     '    Cross-shore Boundary:',icsh
        write(iunit(i),141)     '      Cellstring File:',trim(CS_str(icsh)%bidfile)        
        write(iunit(i),141)     '      Cellstring Path:',trim(CS_str(icsh)%bidpath)
        write(iunit(i),261)     '      Boundary Cells:',CS_str(iwse)%ncells
      enddo !icsh
    
      !--- Nested Water Level BC (Type 7=NH) ---------------------------------------------
      do iwse=1,nNHstr
        write(iunit(i),241)     '    Nested WSE Boundary: ',iwse
        write(iunit(i),141)     '      Cellstring File:',trim(NH_str(iwse)%bidfile)        
        write(iunit(i),141)     '      Cellstring Path:',trim(NH_str(iwse)%bidpath)
        write(iunit(i),261)     '      Boundary Cells:',NH_str(iwse)%ncells
        idpar=NH_str(iwse)%idpar
        write(iunit(i),241)     '      Parent Simulation:',idpar
        write(iunit(i),141)     '        Control File:',trim(ParSim(idpar)%ctlfilepar)
        write(iunit(i),141)     '        Grid File:',trim(ParSim(idpar)%grdfilepar)
        write(iunit(i),141)     '        WSE Input:'
        write(iunit(i),141)     '          File:',trim(ParSim(idpar)%wsefilepar)
        if(len_trim(ParSim(idpar)%wsepathpar)>0)then
          write(iunit(i),141)   '          Path:',trim(ParSim(idpar)%wsepathpar)   
        endif
        if(NH_str(iwse)%wseout)then
          write(iunit(i),141)   '        WSE Output:','ON'  
          call fileparts(NH_str(iwse)%wsefile,apath,aname,aext) 
          astring = trim(aname) // '.' // trim(aext)
          write(iunit(i),141)   '          File:',trim(astring)
          if(len_trim(apath)>0)then
            write(iunit(i),141) '          Path:',trim(apath)
          endif
        else
          write(iunit(i),141)   '        WSE Output:','OFF'  
        endif        
        call julian2calendar(ParSim(idpar)%tjuldaypar,iyrpar,imopar,idaypar,ihrpar,iminpar,isecpar)
        write(iunit(i),431)  iyrpar,imopar,idaypar,ihrpar,iminpar,isecpar 
        write(iunit(i),141)     '        Horizontal Projection'
        write(iunit(i),141)     '          Coordinate System:',trim(aHorizCoordSystem(ParSim(idpar)%projpar%iHorizCoordSystem))     
        if(ParSim(idpar)%projpar%iHorizCoordSystem/=22)then
          write(iunit(i),141)   '          Datum:',trim(aHorizDatum(ParSim(idpar)%projpar%iHorizDatum))
          if(ParSim(idpar)%projpar%iHorizZone/=0)then
            write(iunit(i),262) '          Zone:',ParSim(idpar)%projpar%iHorizZone
          endif
        endif      
        write(iunit(i),141)     '          Units:',trim(aHorizUnits(ParSim(idpar)%projpar%iHorizUnits))
      enddo !iwse-str
    
      !--- Nested Water Level and Velocity BC (Type 8=NHV) -----------------------------------
      do iwse=1,nNHVstr
        write(iunit(i),241)     '    Nested WSE and Velocity Boundary: ',iwse
        write(iunit(i),141)     '      Cellstring File:',trim(NHV_str(iwse)%bidfile)        
        write(iunit(i),141)     '      Cellstring Path:',trim(NHV_str(iwse)%bidpath)
        write(iunit(i),261)     '      Boundary Cells:',NHV_str(iwse)%ncells
        idpar = NHV_str(iwse)%idpar
        write(iunit(i),241)     '      Parent Simulation:       ',idpar
        write(iunit(i),141)     '        Control File:',trim(ParSim(idpar)%ctlfilepar)
        write(iunit(i),141)     '        Grid File:',trim(ParSim(idpar)%grdfilepar)
        write(iunit(i),141)     '        WSE Iput:'
        write(iunit(i),141)     '          File:',trim(ParSim(idpar)%wsefilepar)
        if(len_trim(ParSim(idpar)%wsepathpar)>0)then
          write(iunit(i),141)   '          Path:',trim(ParSim(idpar)%wsepathpar)
        endif
        write(iunit(i),141)     '        Vel Input: '
        write(iunit(i),141)     '          File:',trim(ParSim(idpar)%velfilepar)
        if(len_trim(ParSim(idpar)%velpathpar)>0)then
          write(iunit(i),141)   '          Path:',trim(ParSim(idpar)%velpathpar)   
        endif
        if(NHV_str(iwse)%wseout)then
          write(iunit(i),141)   '        WSE Output:','ON'  
          call fileparts(NHV_str(iwse)%wsefile,apath,aname,aext) 
          astring = trim(aname) // '.' // trim(aext)
          write(iunit(i),141)   '          File:',trim(astring)
          if(len_trim(apath)>0)then
            write(iunit(i),141) '          Path:',trim(apath)
          endif
        else
          write(iunit(i),141)   '        WSE Output:','OFF'
        endif
        if(NHV_str(iwse)%velout)then
          write(iunit(i),141)   '        Velocity Output:','ON'
          call fileparts(NHV_str(iwse)%velfile,apath,aname,aext)           
          astring = trim(aname) // '.' // trim(aext)
          write(iunit(i),141)   '          File:',trim(astring)
          if(len_trim(apath)>0)then
            write(iunit(i),141) '          Path:',trim(apath)
          endif
        else
          write(iunit(i),141)   '        Velocity Output:','OFF'  
        endif
        call julian2calendar(ParSim(idpar)%tjuldaypar,iyrpar,imopar,idaypar,ihrpar,iminpar,isecpar)
        write(iunit(i),431)  iyrpar,imopar,idaypar,ihrpar,iminpar,isecpar 
        write(iunit(i),141)     '        Horizontal Projection'
        write(iunit(i),141)     '          Coordinate System:',trim(aHorizCoordSystem(ParSim(idpar)%projpar%iHorizCoordSystem))     
        if(ParSim(idpar)%projpar%iHorizCoordSystem/=22)then
          write(iunit(i),141)   '          Datum:',trim(aHorizDatum(ParSim(idpar)%projpar%iHorizDatum))
          if(ParSim(idpar)%projpar%iHorizZone/=0)then
            write(iunit(i),262) '          Zone:',ParSim(idpar)%projpar%iHorizZone
          endif
        endif      
        write(iunit(i),141)     '          Units:',trim(aHorizUnits(ParSim(idpar)%projpar%iHorizUnits))
      enddo !iwse-str
    
      !--- Tidal Database WSE BC (Type 9=NTH) -----------------------------------
      do iwse=1,nNTHstr
        write(iunit(i),241)     '    Tidal DB WSE Boundary:',iwse
        write(iunit(i),141)     '      Cellstring File:',trim(NTH_str(iwse)%bidfile)        
        write(iunit(i),141)     '      Cellstring Path:',trim(NTH_str(iwse)%bidpath)        
        write(iunit(i),261)     '      Boundary Cells:',NTH_str(iwse)%ncells
        write(iunit(i),241)     '      Number of Constituents:',NTH_str(iwse)%ntc
        ntc = NTH_str(iwse)%ntc
        write(iunit(i),541)
        write(iunit(i),542)
        write(iunit(i),543)
        write(iunit(i),544)
        do k=1,ntc
          nbndcells = NTH_str(iwse)%ncells          
          ampmin = minval(NTH_str(iwse)%amp(1:nbndcells,k))
          ampmax = maxval(NTH_str(iwse)%amp(1:nbndcells,k))          
          phasemin = minval(NTH_str(iwse)%phase(1:nbndcells,k))*rad2deg          
          phasemax = maxval(NTH_str(iwse)%phase(1:nbndcells,k))*rad2deg
          call avgampdir(nbndcells,NTH_str(iwse)%amp(1:nbndcells,k),&
                 NTH_str(iwse)%phase(1:nbndcells,k),ampavg,phaseavg)
          phaseavg = phaseavg*rad2deg
          write(iunit(i),453) k,NTH_str(iwse)%name(k),NTH_str(iwse)%speed(k)*rad2deg,&
                ampmin,ampmax,ampavg,phasemin,phasemax,phaseavg,NTH_str(iwse)%f(k),NTH_str(iwse)%vu(k)*rad2deg
        enddo
        if(NTH_str(iwse)%wseadjust)then
          write(iunit(i),141)   '      WSE Forcing Correction:','ON'
        else
          write(iunit(i),141)   '      WSE Forcing Correction:','OFF'
        endif
        if(NTH_str(iwse)%wseout)then
          write(iunit(i),141)   '      WSE Output:','ON'  
        else
          write(iunit(i),141)   '      WSE Output:','OFF'  
        endif
        write(iunit(i),141)     '      Tidal Database Horizontal Projection'
        write(iunit(i),141)     '        Coordinate System:',trim(aHorizCoordSystem(NTH_str(iwse)%projtdb%iHorizCoordSystem))     
        if(NTH_str(iwse)%projtdb%iHorizCoordSystem/=22)then
          write(iunit(i),141)   '        Datum:',trim(aHorizDatum(NTH_str(iwse)%projtdb%iHorizDatum))
          if(NTH_str(iwse)%projtdb%iHorizZone>0)then
            write(iunit(i),262) '        Zone:',NTH_str(iwse)%projtdb%iHorizZone
          endif
        endif      
        write(iunit(i),141)     '        Units:',trim(aHorizUnits(NTH_str(iwse)%projtdb%iHorizUnits))
      enddo !iwse-str
    
      !--- Tidal Database WSE and Velocity BC (Type 10=NTHV) --------------------------------
      do iwse=1,nNTHVstr
        write(iunit(i),241)     '    Tidal DB WSE and Vel Boundary:',iwse
        write(iunit(i),141)     '      Cellstring File:',trim(NTHV_str(iwse)%bidfile)        
        write(iunit(i),141)     '      Cellstring Path:',trim(NTHV_str(iwse)%bidpath)
        write(iunit(i),261)     '      Boundary Cells:',NTHV_str(iwse)%ncells
        write(iunit(i),241)     '      Number of Constituents:',NTHV_str(iwse)%ntc
        ntc = NTHV_str(iwse)%ntc
        write(iunit(i),541)
        write(iunit(i),542)
        write(iunit(i),543)
        write(iunit(i),544)
        do k=1,ntc
          nbndcells = NTHV_str(iwse)%ncells
          ampmin = minval(NTHV_str(iwse)%amp(1:nbndcells,k))
          ampmax = maxval(NTHV_str(iwse)%amp(1:nbndcells,k))
          phasemin = minval(NTHV_str(iwse)%phase(1:nbndcells,k))*rad2deg
          phasemax = maxval(NTHV_str(iwse)%phase(1:nbndcells,k))*rad2deg
          call avgampdir(nbndcells,NTHV_str(iwse)%amp(1:nbndcells,k),&
                 NTHV_str(iwse)%phase(1:nbndcells,k),ampavg,phaseavg)
          phaseavg = phaseavg*rad2deg
          write(iunit(i),*) k,NTHV_str(iwse)%name(k),ampmin,ampmax,ampavg,&
             phasemin,phasemax,phaseavg,NTHV_str(iwse)%f(k),NTHV_str(iwse)%vu(k)*rad2deg
        enddo
        write(iunit(i),531)
        write(iunit(i),542)
        write(iunit(i),533)
        write(iunit(i),544)
        do k=1,ntc
          nbndcells = NTHV_str(iwse)%ncells
          ampmin = minval(NTHV_str(iwse)%ampu(1:nbndcells,k))
          ampmax = maxval(NTHV_str(iwse)%ampu(1:nbndcells,k))
          phasemin = minval(NTHV_str(iwse)%phaseu(1:nbndcells,k))*rad2deg
          phasemax = maxval(NTHV_str(iwse)%phaseu(1:nbndcells,k))*rad2deg
          call avgampdir(nbndcells,NTHV_str(iwse)%ampu(1:nbndcells,k),&
                 NTHV_str(iwse)%phaseu(1:nbndcells,k),ampavg,phaseavg)
          phaseavg = phaseavg*rad2deg
          write(iunit(i),453) k,NTHV_str(iwse)%name(k),NTHV_str(iwse)%speed(k)*rad2deg,&
                ampmin,ampmax,ampavg,phasemin,phasemax,phaseavg,&
                NTHV_str(iwse)%f(k),NTHV_str(iwse)%vu(k)*rad2deg
        enddo
        write(iunit(i),532)
        write(iunit(i),542)
        write(iunit(i),533)
        write(iunit(i),544)
        do k=1,ntc
          nbndcells = NTHV_str(iwse)%ncells
          ampmin = minval(NTHV_str(iwse)%ampv(1:nbndcells,k))
          ampmax = maxval(NTHV_str(iwse)%ampv(1:nbndcells,k))
          phasemin = minval(NTHV_str(iwse)%phasev(1:nbndcells,k))*rad2deg
          phasemax = maxval(NTHV_str(iwse)%phasev(1:nbndcells,k))*rad2deg
          call avgampdir(nbndcells,NTHV_str(iwse)%ampv(1:nbndcells,k),&
                 NTHV_str(iwse)%phasev(1:nbndcells,k),ampavg,phaseavg)
          phaseavg = phaseavg*rad2deg
          write(iunit(i),453) k,NTHV_str(iwse)%name(k),NTHV_str(iwse)%speed(k)*rad2deg,&
                ampmin,ampmax,ampavg,phasemin,phasemax,phaseavg,&
                NTHV_str(iwse)%f(k),NTHV_str(iwse)%vu(k)*rad2deg
        enddo
        if(NTHV_str(iwse)%wseadjust)then
          write(iunit(i),141)   '      WSE Forcing Correction:','ON'
        else
          write(iunit(i),141)   '      WSE Forcing Correction:','OFF'
        endif
        if(NTHV_str(iwse)%wseout)then
          write(iunit(i),141)   '      WSE Output:','ON'  
        else
          write(iunit(i),141)   '      WSE Output:','OFF'  
        endif
        if(NTHV_str(iwse)%velout)then
          write(iunit(i),141)   '      Velocity Output:','ON'
        else
          write(iunit(i),141)   '      Velocity Output:','OFF'  
        endif
        write(iunit(i),141)     '      Tidal Database Horizontal Projection'
        write(iunit(i),141)     '        Coordinate System:',trim(aHorizCoordSystem(NTHV_str(iwse)%projtdb%iHorizCoordSystem))     
        if(NTHV_str(iwse)%projtdb%iHorizCoordSystem/=22)then
          write(iunit(i),141)   '        Datum:',trim(aHorizDatum(NTHV_str(iwse)%projtdb%iHorizDatum))
          if(NTHV_str(iwse)%projtdb%iHorizZone/=0)then
            write(iunit(i),262) '        Zone:',NTHV_str(iwse)%projtdb%iHorizZone
          endif
        endif      
        write(iunit(i),141)     '        Units:',trim(aHorizUnits(NTHV_str(iwse)%projtdb%iHorizUnits))
      enddo !iwse-str
    enddo
    
    close(dgunit)
    
    !All boundary information has been read in, write the appropriate information out at this time.
    if(write_ascii_input)then
#ifdef _WIN32
      call write_datasets_to_ascii()   !added 02/21/18 meb        !Only possible on Windows
#endif      
      call write_boundary_to_ascii()   !added 02/15/18 meb
      call write_wind_to_ascii()       !added 02/28/18 meb
      if(allocated(wndspeed)) deallocate(wndspeed,wnddirection)   !Free up memory no longer needed
      !call write_wave_to_ascii()
    endif
    
    return
    end subroutine bnd_print

!****************************************************
    subroutine write_wind_to_ascii()
! Write out ASCII wind file  
! Completed by Mitch Brown, 02/28/2018
    use const_def, only: deg2rad
    use geo_def,   only: azimuth_fl
    use met_def
    use met_lib,   only: wind_heightcorr
    use prec_def,  only: ikind
#ifdef _WIN32
    use IFPORT
#endif
    use diag_lib,  only: diag_print_error

    implicit none
    integer            :: i,ival,velunit,dirunit
    character(len=100) :: windvelfile, winddirfile, astring, astring2, astring3, aname, apath, indirpath
    character(len=10)  :: aext
    logical            :: found,created
    
    if(windfile(1:1)==' ') then   !Windfile is empty, nothing to do here.
      return
    endif
    velunit = 700
    dirunit = 701
    call fileparts(windfile,apath,aname,aext)
    select case(aext)
    case('h5')
      ival=index(windfile,'_mp')-1
    case('xys')
      return        !Files were already read from this format.  Do not overwrite the files.
    end select
    
    !Save all these files to a subdirectory named "ASCII_Input"
    indirpath='ASCII_Input'
#ifdef _WIN32    
    inquire(directory=trim(indirpath), exist=found)
    if(.not.found) then
      created=MakeDirQQ(trim(indirpath))
      if(.not.created)then
        call diag_print_error('Failed to create subdirectory- '//trim(indirpath))
      endif
    endif
#else
    inquire(file=trim(indirpath), exist=found)
    if (.not.found) then 
      call system('mkdir '//(trim(indirpath)))
    endif
#endif
      
    windvelfile = trim(indirpath)// '/' //trim(aname(1:ival))// '_wind_vel.xys'
    winddirfile = trim(indirpath)// '/' //trim(aname(1:ival))// '_wind_dir.xys'
    open(velunit,file=windvelfile)
    open(dirunit,file=winddirfile)
    
    write(velunit,'(A3,x,i0,x,i0,x,A1,A1)') 'XYS',2,nwtimes,'"','"'
    write(dirunit,'(A3,x,i0,x,i0,x,A1,A1)') 'XYS',2,nwtimes,'"','"'

    do i=1,nwtimes
      write(astring, '(F10.1)') wndtimes(i)
      write(astring2,'(F10.1)') wndspeed(i)
      write(astring3,'(F10.1)') wnddirection(i)
      write(velunit,'(A)') trim(adjustl(astring))//' '//trim(adjustl(astring2))  !write speed
      write(dirunit,'(A)') trim(adjustl(astring))//' '//trim(adjustl(astring3))  !write direction
    enddo
    close(velunit)
    close(dirunit)
    
    return
    end subroutine write_wind_to_ascii    
    
#ifdef _WIN32
!only possible on Windows with XMDF
!****************************************************
    subroutine write_datasets_to_ascii()
! Write separate files containing grid-specific datasets that may be used in a CMS run.  
! An example is the bottom friction dataset that contains a Mannings or Roughness value 
!   for each cell in the grid.
! Completed by Mitch Brown, 02/21/2018
    
    use cms_def, only: dsetList,ndsets
    use size_def, only: ncells,ncellsD
    use prec_def, only: ikind
    use geo_def, only: lat,lon
    use out_lib, only: writescalTxt,writevecTxt
    use in_xmdf_lib, only: readscalh5,readvech5
    
    implicit none
    integer :: i, j, ibegin, iend, ierr
    character(len=200) :: apath,aname,astring,dsetname,indirpath
    character(len=10) :: aext
    real(ikind) :: var(ncellsD),vecx(ncellsD),vecy(ncellsD)
        
200 FORMAT (I0,x,I1,5x,'!Number of values, Number of dimensions')    
201 FORMAT (25(E10.3,2x))
202 FORMAT (25(E10.3,x,E10.3,2x))
    do i=1,ndsets
      call fileparts(trim(dsetList(i)%filename),apath,aname,aext)
      ibegin=index(dsetList(i)%path,'/',.true.)+1
      iend=len_trim(dsetList(i)%path)
      !Set datasetname for file
      dsetname=dsetList(i)%path
      dsetname=trim(dsetname(ibegin:iend))
      astring = trim(aname)//"_"//trim(dsetname)//'.txt'
      indirpath='ASCII_Input'
      
      !Determine dimension of dataset and write out appropriately
      select case (dsetList(i)%ndim)
      case(0)
        !Latitudes or Longitudes.  Write to scalar file.  Lat/Long section added 5/21/2018.
        if (dsetname=='Lats') then
          allocate(lat(ncellsD))
          call read_latlon_dataset(lat,'Lats')  !Required because latitude dataset is saved as a group property
          call writescalTxt(dsetList(i)%filename,trim(indirpath),trim(dsetname),lat,1)
          deallocate(lat)
        else
          allocate(lon(ncellsD))
          call read_latlon_dataset(lon,'Lons')  !Required because longitude dataset is saved as a group property
          call writescalTxt(dsetList(i)%filename,trim(indirpath),trim(dsetname),lon,1)
          deallocate(lon)
        endif  
      case(1)
        !Write to scalar file
        call readscalh5(dsetList(i)%filename,dsetList(i)%path,var,ierr)
        call writescalTxt(dsetList(i)%filename,trim(indirpath),trim(dsetname),var,1)
      case(2)
        !Write to vector file
        call readvech5(dsetList(i)%filename,dsetList(i)%path,vecx,vecy,ierr)
        call writevecTxt(dsetList(i)%filename,trim(indirpath),trim(dsetname),vecx,vecy,2)
      end select
    enddo
    
    return
    end subroutine write_datasets_to_ascii
#endif

!****************************************************
    subroutine write_boundary_to_ascii()
! Write files containing boundary forcing information so that it will be available to be read in later.
! Go in standard order of boundaries and write out a .bid file and .xys file containing the forcing information.
! Multiple WSE extracted boundaries will be more difficult and added later.
! Completed by Mitch Brown, 02/20/2018
    
    use bnd_def,  only: nbndstr,bnd_str,nqstr,q_str,nthstr,th_str,nhstr,h_str
    use out_def,  only: outprefix
    use out_lib,  only: write_xys
    use prec_def, only: ikind
    use geo_def,  only: mapid
    use diag_lib, only: diag_print_error
    use sal_def,  only: sal_str,nsalstr

#ifdef _WIN32
    use IFPORT
#endif
    
    implicit none
    integer              :: kunit,i,j, val, ibnd, nstrcells, correctID, id, ntimes
    integer, allocatable :: ibndcells(:)
    character(len=200)   :: bidoutfile, xysoutfile, astring,astring2, abnd, indirpath
    logical              :: found,created
    
    !Save all these files to a subdirectory named "ASCII_Input"
    indirpath='ASCII_Input'
#ifdef _WIN32    
    inquire(directory=trim(indirpath), exist=found)
    if(.not.found) then
      created=MakeDirQQ(trim(indirpath))
      if(.not.created)then
        call diag_print_error('Failed to create subdirectory- '//trim(indirpath))
      endif
    endif
#else
    inquire(file=trim(indirpath), exist=found)
    if (.not.found) then 
      call system('mkdir '//(trim(indirpath)))
    endif
#endif

    !Write Boundary ID file
    kunit=912
    bidoutfile = trim(indirpath)// '/' //trim(outprefix) // '.bid'   !ASCII Boundary ID File  
    open(kunit,file=bidoutfile)
    write(astring,'(I4)') nbndstr
    write(kunit,*) adjustl(trim(astring)),' !# of cellstrings'
    do ibnd=1,nbndstr
      !Determine cells without repeats
      allocate(ibndcells(0:bnd_str(ibnd)%ncells))
      ibndcells = 0
      nstrcells = 0
      do j=1,bnd_str(ibnd)%ncells
        i = mapid(bnd_str(ibnd)%cells(j))
        if (i .ne. ibndcells(nstrcells)) then
          nstrcells = nstrcells + 1
          ibndcells(nstrcells) = i
        endif
      enddo
      !Write cell ID's
      write(astring,'(I5)') nstrcells
      write(kunit,*) adjustl(trim(astring)), '!# of cells in this cellstring'
      do j=1,nstrcells
        write(astring,'(I7)') ibndcells(j)
        write(kunit,*) adjustl(trim(astring))
      enddo
      
      !Write forcing information
      select case(bnd_str(ibnd)%ibndtype)
      case(1) !Flux boundaries - ISTRTYPE = 1
        !Find correct boundary in list
        correctID = -1
        do i=1,nqstr
          if(q_str(i)%IDNUM == bnd_str(ibnd)%IDNUM) then
            correctID = i
            exit
          endif
        enddo
        if(correctID .lt. 0) call diag_print_error('Could not locate correct Boundary string')
        
        !Write curve information for correct boundary
        if (q_str(correctID)%ifluxmode == 2) then    !if iFluxMode==2(CURVE) instead of 1(CONSTANT)
          write(abnd,'(I0)') q_str(correctID)%idnum  !Use ID number from SMS (Boundary_#2 = 2) as the number to write out.          
          xysoutfile = trim(indirpath)// '/' //trim(outprefix) // '_q_'// trim(abnd) // '.xys'    !ASCII Boundary Forcing file
          open(kunit+1,file=xysoutfile)
          write(kunit+1,'(A3,x,i0,x,i0,x,A1,A1)') 'XYS',2,q_str(correctID)%ntimes,'"','"'
          do i=1,q_str(correctID)%ntimes
            if(q_str(correctID)%qcurv(i) == -0.0) then
              q_str(correctID)%qcurv(i) = 0.0
            endif
            write(astring, '(F10.2)') q_str(correctID)%times(i)
            write(astring2,'(F10.2)') q_str(correctID)%qcurv(i)
            write(kunit+1,'(A)') trim(adjustl(astring))//' '//trim(adjustl(astring2))
          enddo
          close(kunit+1)
        endif
        
        !Write out salinity information, if present
        if (nsalstr > 0) then
          id = q_str(correctID)%idnum
          ntimes = sal_str(correctID)%ntimes
          call write_xys (id, indirpath, outprefix, '_sal_', ntimes, sal_str(correctID)%timesal,sal_str(correctID)%val)
        endif     
                      
      case(2) !Tidal Harmonic boundaries - ISTRTYPE = 2 (Nothing to write here)
        continue
      case(3) !Single WSE boundaries - ISTRTYPE = 3
        !Find correct boundary in list
        do i=1,nhstr
          if(h_str(i)%IDNUM == bnd_str(ibnd)%IDNUM) then
            correctID = i
            exit
          endif
        enddo
        if(correctID .lt. 0) call diag_print_error('Error: Could not locate correct Boundary string')
        !Write curve information for correct boundary
        write(abnd,'(I0)') h_str(correctID)%idnum  !Use ID number from SMS (Boundary_#2 = 2) as the number to write out.          
        xysoutfile = trim(indirpath)// '/' //trim(outprefix) // '_h_'// trim(abnd) // '.xys'    !ASCII Boundary Forcing file
        open(kunit+1,file=xysoutfile)
        write(kunit+1,'(A3,x,i0,x,i0,x,A1,A1)') 'XYS',2,h_str(correctID)%ntimes,'"','"'
        do i=1,h_str(correctID)%ntimes
          if(h_str(correctID)%wsecurv(i) == -0.0) then
            h_str(correctID)%wsecurv(i) = 0.0
          endif
          write(astring, '(F10.2)') h_str(correctID)%times(i)
          write(astring2,'(F10.4)') h_str(correctID)%wsecurv(i)
          write(kunit+1,'(A)') trim(adjustl(astring))//' '//trim(adjustl(astring2))
        enddo
        close(kunit+1)
        
        !Write out salinity information, if present
        if (nsalstr > 0) then
          id = h_str(correctID)%idnum
          ntimes = sal_str(correctID)%ntimes
          call write_xys (id, indirpath, outprefix, '_sal_', ntimes, sal_str(correctID)%timesal,sal_str(correctID)%val)
        endif     

      end select  

      deallocate(ibndcells)
    enddo
    close(kunit)
    
    return
    end subroutine write_boundary_to_ascii

!***********************************************************************
    subroutine bound_uv()
! Applies the boundary conditions to the momentum equations
! written by Weiming Wu, NCCHE, Nov. 2008
! modified by Alex Sanchez, USACE-CHL
!***********************************************************************
#include "CMS_cpp.h"
    use size_def
    use geo_def
    use flow_def    
    use bnd_def
    use comvarbl, only: ntsch,ctsch1,ctsch2,relax
    use fric_def, only: wallfric,z0,wallfac
    use fric_lib, only: wall_coef
    use comp_lib, only: upwindcoef
    use const_def, only: small,pi
#ifdef DIAG_MODE
    use diag_lib, only: diag_print_error
#endif   
    use prec_def, only: ikind
    
    implicit none
    integer :: i,k,iriv,im,iwse,icsh,nck
    real(ikind) :: delta,upar,vpar,uvpar,z0wall,uvnorm,val,ddk

!--- River/Flux BC (Type 1=Q) ---------------------
    do iriv=1,nQstr
      call bnd_vel_flux(Q_str(iriv)%ncells,Q_str(iriv)%cells,Q_str(iriv)%faces)
    enddo

!--- Tidal/Harmonic BC (Type 2=TH) ----------------
    do iwse=1,nTHstr
      call bnd_vel_openwse(TH_str(iwse)%ncells,TH_str(iwse)%cells,TH_str(iwse)%faces)
    enddo
    
!--- Single Water Level BC (Type 3=H) ------------------------      
    do iwse=1,nHstr
      !call bnd_vel_openwse(H_str(iwse)%ncells,H_str(iwse)%cells,H_str(iwse)%faces)
      call bnd_vel_openwse_test(H_str(iwse)%ncells,H_str(iwse)%cells,H_str(iwse)%faces)  
    enddo
    
!--- Multiple Water Level BC (Type 4=MH) --------------------------      
    do iwse=1,nMHstr
      call bnd_vel_openwse(MH_str(iwse)%ncells,MH_str(iwse)%cells,MH_str(iwse)%faces)
    enddo    

!--- Multiple Water Level and Velocity BC (Type 5=MHV) ----------      
    do iwse=1,nMHVstr
      call bnd_vel_openwsevel(MHV_str(iwse)%ncells,MHV_str(iwse)%cells,MHV_str(iwse)%faces,MHV_str(iwse)%ubnd,MHV_str(iwse)%vbnd)
    enddo

!--- Nested Water Level BC (Type 7=NH) ----------------------------      
    do iwse=1,nNHstr
      call bnd_vel_openwse(NH_str(iwse)%ncells,NH_str(iwse)%cells,NH_str(iwse)%faces)
    enddo    

!--- Nested Water Level and Velocity BC (Type 8=NHV) -------------      
    do iwse=1,nNHVstr
      call bnd_vel_openwsevel(NHV_str(iwse)%ncells,NHV_str(iwse)%cells,NHV_str(iwse)%faces,NHV_str(iwse)%ubnd,NHV_str(iwse)%vbnd)
    enddo

!--- Tidal Database Water Level BC (Type 9=NTH) ----------------------------      
    do iwse=1,nNTHstr
      call bnd_vel_openwse(NTH_str(iwse)%ncells,NTH_str(iwse)%cells,NTH_str(iwse)%faces)
    enddo    

!--- Tidal Database Water Level and Velocity BC (Type 10=NTHV) -------------      
    do iwse=1,nNTHVstr
      call bnd_vel_openwsevel(NTHV_str(iwse)%ncells,NTHV_str(iwse)%cells,&
            NTHV_str(iwse)%faces,NTHV_str(iwse)%ubnd,NTHV_str(iwse)%vbnd)  
    enddo    

!********************************************************
!Note: Up to this point, sp is the same for u and v.
!Now sp is separated by the wall friction
!********************************************************

!--- Dry nodes --------------------------------------------------------
    if(wallfric)then
!$OMP PARALLEL DO PRIVATE(i,k,nck,delta,uvnorm,upar,vpar,uvpar,z0wall,val)            
      do i=1,ncells
        if(iwet(i)==1)then
          spu(i)=sp(i)
          spv(i)=sp(i) 
          do k=1,ncface(i)
            nck=cell2cell(k,i)
            if(iwet(nck)==0)then    !side is dry
              z0wall=z0(i)*wallfac
              if(abs(fny(k,i))>0.999)then   !north/south face
                upar=u(i) !+ry(k,i)*duy(i) !Wall parallel velocity                         
                delta=abs(ry(k,i)) !=0.5*dy(i)
                spu(i)=spu(i)-ds(k,i)*hk(k,i)*wall_coef(delta,z0wall)*abs(upar)
              elseif(abs(fnx(k,i))>0.999)then !east/west face
                vpar=v(i) !+rx(k,i)*dvx(i) !Wall parallel velocity
                delta=abs(rx(k,i)) !=0.5*dx(i)
                spv(i)=spv(i)-ds(k,i)*hk(k,i)*wall_coef(delta,z0wall)*abs(vpar)
              else  !Oblique face
                delta=sqrt(rx(k,i)*rx(k,i)+ry(k,i)*ry(k,i)) !Distance to wall
                uvnorm=u(i)*fnx(k,i)+v(i)*fny(k,i)  !Wall normal velocity magnitude
                upar=u(i)-uvnorm*fnx(k,i) !Wall parallel velocity x-component
                vpar=v(i)-uvnorm*fny(k,i) !Wall parallel velocity y-component
                uvpar=sqrt(upar*upar+vpar*vpar) !Magnitude of wall parallel velocity
                val=ds(k,i)*hk(k,i)*wall_coef(delta,z0wall)*abs(uvpar)
                spu(i)=spu(i)-val
                spv(i)=spv(i)-val
              endif
            endif
          enddo            
        else !if(iwet(i)==0)then
          u(i)=0.0; v(i)=0.0
          spu(i)=-1.0; spv(i)=-1.0; sp(i)=-1.0
          su(i)=0.0;   sv(i)=0.0
        endif
      enddo
!$OMP END PARALLEL DO
    else   !No wallfriction
!$OMP PARALLEL DO PRIVATE(i)       
      do i=1,ncells
        u(i)=iwet(i)*u(i)
        v(i)=iwet(i)*v(i)
        sp(i)=iwet(i)*sp(i)+real(iwet(i)-1,kind=ikind) !sp for wet, -1.0 for dry
        spu(i)=sp(i)
        spv(i)=spu(i)
        su(i)=iwet(i)*su(i)
        sv(i)=iwet(i)*sv(i)
      enddo    
!$OMP END PARALLEL DO
    endif
    
#ifdef DIAG_MODE
    do i=1,ncells
      do k=1,ncface(i)
        nck=cell2cell(k,i)  
        if(abs(acoef(k,i))>=1.0e-20 .and. iwet(i)*iwet(nck)==0)then
          call diag_print_error('Non-zero coefficient at dry boundary')
        endif
      enddo
    enddo
#endif

!--- Cross-shore BC (Type 6=CS) ---------------------------    
    call xshore_uv
    do icsh=1,nCSstr
      do im=1,CS_str(icsh)%ncells
        i=CS_str(icsh)%cells(im)
        if(iwet(i)==0) cycle
        k=CS_str(icsh)%faces(im)
        nck=cell2cell(k,i)
        ddk=visk(k,i)*hk(k,i)*dsxy(k,i)
        acoef(k,i)=upwindcoef(ddk,flux(k,i))
        if(idirface(k,i)==2 .or. idirface(k,i)==4)then
          u(nck)=CS_str(icsh)%ucsh(im)
          if(flux(k,i)<0.0)then !Inflow 
            su(i)=su(i)+acoef(k,i)*u(nck)
            spu(i)=spu(i)-acoef(k,i)
            !v(nck)=v(i)
            v(nck)=v(i)*h(i)/h(nck)*iwet(i)
          else
            u(nck)=u(i)*h(i)/h(nck)*iwet(i)
            v(nck)=v(i)*h(i)/h(nck)*iwet(i)
          endif
        else
          v(nck)=CS_str(icsh)%vcsh(im)
          if(flux(k,i)<0.0)then !Inflow
            sv(i)=sv(i)+acoef(k,i)*v(nck)   
            spv(i)=spv(i)-acoef(k,i)
            !u(nck)=u(i)
            u(nck)=u(i)*h(i)/h(nck)*iwet(i)
          else
            u(nck)=u(i)*h(i)/h(nck)*iwet(i)
            v(nck)=v(i)*h(i)/h(nck)*iwet(i)
          endif
        endif        
        acoef(k,i)=0.0
      enddo
    enddo    

!--- Structures -------------------    
    call struct_uv

    return
    end subroutine bound_uv

!***********************************************************************
    subroutine bound_pp()
! Applies the boundary conditions to the pressure correction equation
! written by Weiming Wu, NCCHE, NOV. 2008
!***********************************************************************
    use size_def, only: ncells
    use geo_def, only: ncface,cell2cell,dsxy,kkface,idirface
    use flow_def
    use bnd_def
    use prec_def, only: ikind
    
    implicit none
    integer :: i,j,k,kk,sumwet,nck,iriv,iwse,icsh,kkdf
    real(ikind) :: aploc,acoefik

!--- Dry nodes ---------------------------------------------------------------
!$OMP PARALLEL DO PRIVATE(i,k,nck,sumwet,aploc)
    do i=1,ncells
       !if(iwet(i)==1)then
       !  acoef(1:ncface(i),i)=iwet(cell2cell(1:ncface(i),i))*acoef(1:ncface(i),i)           
       !elseif(iwet(i)==0)then
       if(iwet(i)==0)then 
         sumwet=0.0
         do k=1,ncface(i)
           nck=cell2cell(k,i)
           if(nck<=ncells) sumwet=sumwet+iwet(nck)
         enddo             
         aploc=0.0
         su(i)=0.0
         if(sumwet>=1)then
           do k=1,ncface(i)
             nck=cell2cell(k,i)
             acoef(k,i)=iwet(nck)  
             aploc=aploc+acoef(k,i)
             su(i)=su(i)+acoef(k,i)*p(nck)
           enddo  
         else
           do k=1,ncface(i)
             nck=cell2cell(k,i)
             if(nck>ncells)then
               acoef(k,i)=0.0
             else
               acoef(k,i)=dsxy(k,i)
               aploc=aploc+acoef(k,i)
               su(i)=su(i)+acoef(k,i)*p(nck)
             endif             
           enddo
         endif
         !if(aploc<=0.000001) then
         !  write(*,*) 'WARNING: ap= ',aploc, 'at node ', i
         !  write(*,*) 'Press any key to continue'
         !  read(*,*)
         !  stop
         !endif
         su(i)=su(i)-aploc*p(i)
         sp(i)=0.0
       endif
    enddo
!$OMP END PARALLEL DO

!Note: dry faces are already included in coefsourcesink_pp
!!--- Wall BC ---------------------------------------------------------
!    do j=1,W_str%ncells
!      i=W_str%cells(j)
!      k=W_str%faces(j)
!      acoef(k,i)=0.0
!    enddo
    
!--- River BC (Type 1=Q) --------------------------------------------
    do iriv=1,nQstr
      do j=1,Q_str(iriv)%ncells
        i=Q_str(iriv)%cells(j)
        k=Q_str(iriv)%faces(j)
        acoef(k,i)=0.0
      enddo
    enddo

!--- Tidal/Harmonic BC (Type 2=TH) ---------------------------------------
!Water level is first assigend at the center cell of the cell string
!and then adjusted for wind and wave setup starting from the center cell
    do iwse=1,nTHstr      
      if(TH_str(iwse)%wseadjust)then
        call bnd_wse_adjust(TH_str(iwse)%ncells,TH_str(iwse)%cells,TH_str(iwse)%faces,TH_str(iwse)%wsebnd,TH_str(iwse)%wseadj)
        call bnd_pp(TH_str(iwse)%ncells,TH_str(iwse)%cells,TH_str(iwse)%faces,TH_str(iwse)%wseadj)
      else
        call bnd_pp(TH_str(iwse)%ncells,TH_str(iwse)%cells,TH_str(iwse)%faces,TH_str(iwse)%wsebnd)  
      endif
    enddo !iwse
    
!--- Single Water Surface Elevation BC (Type 3=H) -------------------------
!Water level is first assigend at the center cell of the cell string
!and then adjusted for wind and wave setup starting from the center cell      
    do iwse=1,nHstr        
      if(H_str(iwse)%wseadjust)then
        call bnd_wse_adjust(H_str(iwse)%ncells,H_str(iwse)%cells,H_str(iwse)%faces,H_str(iwse)%wsebnd,H_str(iwse)%wseadj)
        call bnd_pp(H_str(iwse)%ncells,H_str(iwse)%cells,H_str(iwse)%faces,H_str(iwse)%wseadj)
      else
        call bnd_pp(H_str(iwse)%ncells,H_str(iwse)%cells,H_str(iwse)%faces,H_str(iwse)%wsebnd)
      endif  
    enddo !iwse 

!--- Multiple Water Level BC (Type 4=MH) ---------------------------
    do iwse=1,nMHstr
      call bnd_pp(MH_str(iwse)%ncells,MH_str(iwse)%cells,MH_str(iwse)%faces,MH_str(iwse)%wsebnd)
    enddo !iwse
    
!--- Multiple Water Level and Velocity BC (Type 5=MHV) -----------------------
    do iwse=1,nMHVstr
      call bnd_pp(MHV_str(iwse)%ncells,MHV_str(iwse)%cells,MHV_str(iwse)%faces,MHV_str(iwse)%wsebnd)
    enddo !iwse        

!--- Cross-shore BC (Type 6=CS) ------------------------------------
    call xshore_wse
    do icsh=1,nCSstr
    !  call bnd_pp(CS_str(icsh)%ncells,CS_str(icsh)%cells,CS_str(icsh)%faces,CS_str(icsh)%wsecsh)
      do j=1,CS_str(icsh)%ncells
        i=CS_str(icsh)%cells(j)
        k=CS_str(icsh)%faces(j)
        nck=cell2cell(k,i)        
        if(flux(k,i)>0.0)then !Outflow
          pp(nck)=CS_str(icsh)%wsecsh(j)*grav-p(nck)   !Set given wse at the dummy node  
          if(ncface(i)==4)then
            acoefik=acoef(kkface(k),i)
          else
            acoefik=0.0; kkdf=kkface(idirface(k,i))
            do kk=1,ncface(i)    !Assume boundary face does not split
              if(idirface(kk,i)==kkdf) acoefik=acoefik+acoef(kk,i)
            enddo                   
          endif                        
          su(i)=su(i)+acoefik*pp(nck)  !using a_kk because a_k is not well defined  
          sp(i)=sp(i)-acoefik
        endif
        acoef(k,i)=0.0
      enddo !j  
    enddo
    
!--- Nested Water Level BC (Type 7=NH) ------------------------------
    do iwse=1,nNHstr
      call bnd_pp(NH_str(iwse)%ncells,NH_str(iwse)%cells,NH_str(iwse)%faces,NH_str(iwse)%wsebnd)
    enddo !iwse
    
!--- Nested Water Level and Velocity BC (Type 8=NHV) --------------------
    do iwse=1,nNHVstr
      call bnd_pp(NHV_str(iwse)%ncells,NHV_str(iwse)%cells,NHV_str(iwse)%faces,NHV_str(iwse)%wsebnd)  
    enddo !iwse

!--- Tidal Database Water Level BC (Type 9=NTH) ------------------------------
    do iwse=1,nNTHstr
      if(NTH_str(iwse)%wseadjust)then  
        call bnd_wse_adjust(NTH_str(iwse)%ncells,NTH_str(iwse)%cells,NTH_str(iwse)%faces,NTH_str(iwse)%wsebnd,NTH_str(iwse)%wseadj)
        call bnd_pp(NTH_str(iwse)%ncells,NTH_str(iwse)%cells,NTH_str(iwse)%faces,NTH_str(iwse)%wseadj)
      else
        call bnd_pp(NTH_str(iwse)%ncells,NTH_str(iwse)%cells,NTH_str(iwse)%faces,NTH_str(iwse)%wsebnd)  
      endif
    enddo !iwse
    
!--- Tidal Database Water Level and Velocity BC (Type 10=NTHV) --------------------
    do iwse=1,nNTHVstr
      if(NTHV_str(iwse)%wseadjust)then  
        call bnd_wse_adjust(NTHV_str(iwse)%ncells,NTHV_str(iwse)%cells,NTHV_str(iwse)%faces,NTHV_str(iwse)%wsebnd,NTHV_str(iwse)%wseadj)
        call bnd_pp(NTHV_str(iwse)%ncells,NTHV_str(iwse)%cells,NTHV_str(iwse)%faces,NTHV_str(iwse)%wseadj)
      else
        call bnd_pp(NTHV_str(iwse)%ncells,NTHV_str(iwse)%cells,NTHV_str(iwse)%faces,NTHV_str(iwse)%wsebnd)  
      endif
    enddo !iwse    
    
!--- Structures -------------------    
    call struct_pp

    return
    end subroutine bound_pp
    
!***********************************************************************
    subroutine bound_p()
! Applies the boundary conditions to the pressure correction equation
! written by Weiming Wu, NCCHE, NOV. 2008
!***********************************************************************
    use size_def, only: ncells
    use geo_def, only: ncface,cell2cell,dsxy,kkface,idirface
    use flow_def
    use bnd_def
    use prec_def, only: ikind
    
    implicit none
    integer :: i,j,k,kk,sumwet,nck,iriv,iwse,icsh,kkdf
    real(ikind) :: aploc,acoefik

!!--- Dry nodes ---------------------------------------------------------------
!!$OMP PARALLEL DO PRIVATE(i,k,nck,sumwet,aploc)
!    do i=1,ncells
!       !if(iwet(i)==1)then
!       !  acoef(1:ncface(i),i)=iwet(cell2cell(1:ncface(i),i))*acoef(1:ncface(i),i)           
!       !elseif(iwet(i)==0)then
!       if(iwet(i)==0)then 
!         sumwet=0.0
!         do k=1,ncface(i)
!           nck=cell2cell(k,i)
!           if(nck<=ncells) sumwet=sumwet+iwet(nck)
!         enddo             
!         aploc=0.0
!         su(i)=0.0
!         if(sumwet>=1)then
!           do k=1,ncface(i)
!             nck=cell2cell(k,i)
!             acoef(k,i)=iwet(nck)  
!             aploc=aploc+acoef(k,i)
!             su(i)=su(i)+acoef(k,i)*p(nck)
!           enddo  
!         else
!           do k=1,ncface(i)
!             nck=cell2cell(k,i)
!             if(nck>ncells)then
!               acoef(k,i)=0.0
!             else
!               acoef(k,i)=dsxy(k,i)
!               aploc=aploc+acoef(k,i)
!               su(i)=su(i)+acoef(k,i)*p(nck)
!             endif             
!           enddo
!         endif
!         !if(aploc<=0.000001) then
!         !  write(*,*) 'WARNING: ap= ',aploc, 'at node ', i
!         !  write(*,*) 'Press any key to continue'
!         !  read(*,*)
!         !  stop
!         !endif
!         su(i)=su(i)-aploc*p(i)
!         sp(i)=0.0
!       endif
!    enddo
!!$OMP END PARALLEL DO

!Note: dry faces are already included in coefsourcesink_pp
!!--- Wall BC ---------------------------------------------------------
!    do j=1,W_str%ncells
!      i=W_str%cells(j)
!      k=W_str%faces(j)
!      acoef(k,i)=0.0
!    enddo
    
!--- River BC (Type 1=Q) --------------------------------------------
    do iriv=1,nQstr
      do j=1,Q_str(iriv)%ncells
        i=Q_str(iriv)%cells(j)
        k=Q_str(iriv)%faces(j)
        acoef(k,i)=0.0
      enddo
    enddo

!--- Tidal/Harmonic BC (Type 2=TH) ---------------------------------------
!Water level is first assigend at the center cell of the cell string
!and then adjusted for wind and wave setup starting from the center cell
    do iwse=1,nTHstr      
      if(TH_str(iwse)%wseadjust)then
        call bnd_wse_adjust(TH_str(iwse)%ncells,TH_str(iwse)%cells,TH_str(iwse)%faces,TH_str(iwse)%wsebnd,TH_str(iwse)%wseadj)
        call bnd_p(TH_str(iwse)%ncells,TH_str(iwse)%cells,TH_str(iwse)%faces,TH_str(iwse)%wseadj)
      else
        call bnd_p(TH_str(iwse)%ncells,TH_str(iwse)%cells,TH_str(iwse)%faces,TH_str(iwse)%wsebnd)  
      endif
    enddo !iwse
    
!--- Single Water Surface Elevation BC (Type 3=H) -------------------------
!Water level is first assigend at the center cell of the cell string
!and then adjusted for wind and wave setup starting from the center cell      
    do iwse=1,nHstr        
      if(H_str(iwse)%wseadjust)then
        call bnd_wse_adjust(H_str(iwse)%ncells,H_str(iwse)%cells,H_str(iwse)%faces,H_str(iwse)%wsebnd,H_str(iwse)%wseadj)
        call bnd_p(H_str(iwse)%ncells,H_str(iwse)%cells,H_str(iwse)%faces,H_str(iwse)%wseadj)
      else
        call bnd_p(H_str(iwse)%ncells,H_str(iwse)%cells,H_str(iwse)%faces,H_str(iwse)%wsebnd)
      endif  
    enddo !iwse 

!--- Multiple Water Level BC (Type 4=MH) ---------------------------
    do iwse=1,nMHstr
      call bnd_p(MH_str(iwse)%ncells,MH_str(iwse)%cells,MH_str(iwse)%faces,MH_str(iwse)%wsebnd)
    enddo !iwse
    
!--- Multiple Water Level and Velocity BC (Type 5=MHV) -----------------------
    do iwse=1,nMHVstr
      call bnd_p(MHV_str(iwse)%ncells,MHV_str(iwse)%cells,MHV_str(iwse)%faces,MHV_str(iwse)%wsebnd)
    enddo !iwse        

!!--- Cross-shore BC (Type 6=CS) ------------------------------------
!    call xshore_wse
!    do icsh=1,nCSstr
!    !  call bnd_p(CS_str(icsh)%ncells,CS_str(icsh)%cells,CS_str(icsh)%faces,CS_str(icsh)%wsecsh)
!      do j=1,CS_str(icsh)%ncells
!        i=CS_str(icsh)%cells(j)
!        k=CS_str(icsh)%faces(j)
!        nck=cell2cell(k,i)        
!        if(flux(k,i)>0.0)then !Outflow
!          pp(nck)=CS_str(icsh)%wsecsh(j)*grav-p(nck)   !Set given wse at the dummy node  
!          if(ncface(i)==4)then
!            acoefik=acoef(kkface(k),i)
!          else
!            acoefik=0.0; kkdf=kkface(idirface(k,i))
!            do kk=1,ncface(i)    !Assume boundary face does not split
!              if(idirface(kk,i)==kkdf) acoefik=acoefik+acoef(kk,i)
!            enddo                   
!          endif                        
!          su(i)=su(i)+acoefik*pp(nck)  !using a_kk because a_k is not well defined  
!          sp(i)=sp(i)-acoefik
!        endif
!        acoef(k,i)=0.0
!      enddo !j  
!    enddo
    
!--- Nested Water Level BC (Type 7=NH) ------------------------------
    do iwse=1,nNHstr
      call bnd_p(NH_str(iwse)%ncells,NH_str(iwse)%cells,NH_str(iwse)%faces,NH_str(iwse)%wsebnd)
    enddo !iwse
    
!--- Nested Water Level and Velocity BC (Type 8=NHV) --------------------
    do iwse=1,nNHVstr
      call bnd_p(NHV_str(iwse)%ncells,NHV_str(iwse)%cells,NHV_str(iwse)%faces,NHV_str(iwse)%wsebnd)  
    enddo !iwse

!--- Tidal Database Water Level BC (Type 9=NTH) ------------------------------
    do iwse=1,nNTHstr
      if(NTH_str(iwse)%wseadjust)then  
        call bnd_wse_adjust(NTH_str(iwse)%ncells,NTH_str(iwse)%cells,NTH_str(iwse)%faces,NTH_str(iwse)%wsebnd,NTH_str(iwse)%wseadj)
        call bnd_p(NTH_str(iwse)%ncells,NTH_str(iwse)%cells,NTH_str(iwse)%faces,NTH_str(iwse)%wseadj)
      else
        call bnd_p(NTH_str(iwse)%ncells,NTH_str(iwse)%cells,NTH_str(iwse)%faces,NTH_str(iwse)%wsebnd)  
      endif
    enddo !iwse
    
!--- Tidal Database Water Level and Velocity BC (Type 10=NTHV) --------------------
    do iwse=1,nNTHVstr
      if(NTHV_str(iwse)%wseadjust)then  
        call bnd_wse_adjust(NTHV_str(iwse)%ncells,NTHV_str(iwse)%cells,NTHV_str(iwse)%faces,NTHV_str(iwse)%wsebnd,NTHV_str(iwse)%wseadj)
        call bnd_p(NTHV_str(iwse)%ncells,NTHV_str(iwse)%cells,NTHV_str(iwse)%faces,NTHV_str(iwse)%wseadj)
      else
        call bnd_p(NTHV_str(iwse)%ncells,NTHV_str(iwse)%cells,NTHV_str(iwse)%faces,NTHV_str(iwse)%wsebnd)  
      endif
    enddo !iwse    
    
!--- Structures -------------------    
    !call struct_pp

    return
    end subroutine bound_p    
    
!***********************************************************************
    subroutine bound_wse()
! Applies the boundary conditions to the pressure correction equation
! written by Weiming Wu, NCCHE, NOV. 2008
!***********************************************************************
    use size_def, only: ncells
    use geo_def, only: ncface,cell2cell,dsxy,kkface,idirface,zb
    use flow_def
    use bnd_def

    implicit none
    integer :: i,j,k,iwse,iriv,nck

!!--- Dry nodes ---------------------------------------------------------------
!!$OMP PARALLEL DO PRIVATE(i,k,nck,sumwet,aploc)
!    do i=1,ncells
!       !if(iwet(i)==1)then
!       !  acoef(1:ncface(i),i)=iwet(cell2cell(1:ncface(i),i))*acoef(1:ncface(i),i)           
!       !elseif(iwet(i)==0)then
!       if(iwet(i)==0)then 
!         sumwet=0.0
!         do k=1,ncface(i)
!           nck=cell2cell(k,i)
!           if(nck<=ncells) sumwet=sumwet+iwet(nck)
!         enddo             
!         aploc=0.0
!         su(i)=0.0
!         if(sumwet>=1)then
!           do k=1,ncface(i)
!             nck=cell2cell(k,i)
!             acoef(k,i)=iwet(nck)  
!             aploc=aploc+acoef(k,i)
!             su(i)=su(i)+acoef(k,i)*p(nck)
!           enddo  
!         else
!           do k=1,ncface(i)
!             nck=cell2cell(k,i)
!             if(nck>ncells)then
!               acoef(k,i)=0.0
!             else
!               acoef(k,i)=dsxy(k,i)
!               aploc=aploc+acoef(k,i)
!               su(i)=su(i)+acoef(k,i)*p(nck)
!             endif             
!           enddo
!         endif
!         !if(aploc<=0.000001) then
!         !  write(*,*) 'WARNING: ap= ',aploc, 'at node ', i
!         !  write(*,*) 'Press any key to continue'
!         !  read(*,*)
!         !  stop
!         !endif
!         su(i)=su(i)-aploc*p(i)
!         sp(i)=0.0
!       endif
!    enddo
!!$OMP END PARALLEL DO

!Note: dry faces are already included in coefsourcesink_pp
!!--- Wall BC ---------------------------------------------------------
!    do j=1,W_str%ncells
!      i=W_str%cells(j)
!      k=W_str%faces(j)
!      acoef(k,i)=0.0
!    enddo
    
!--- River BC (Type 1=Q) --------------------------------------------
    do iriv=1,nQstr
      do j=1,Q_str(iriv)%ncells
        i=Q_str(iriv)%cells(j)
        k=Q_str(iriv)%faces(j)
        nck=cell2cell(k,i)
        h(nck)=h(i)
        eta(nck)=h(nck)+zb(nck)
      enddo
    enddo

!--- Tidal/Harmonic BC (Type 2=TH) ---------------------------------------
!Water level is first assigend at the center cell of the cell string
!and then adjusted for wind and wave setup starting from the center cell
    do iwse=1,nTHstr      
      if(TH_str(iwse)%wseadjust)then
        call bnd_wse_adjust(TH_str(iwse)%ncells,TH_str(iwse)%cells,TH_str(iwse)%faces,TH_str(iwse)%wsebnd,TH_str(iwse)%wseadj)
        call bnd_wse(TH_str(iwse)%ncells,TH_str(iwse)%cells,TH_str(iwse)%faces,TH_str(iwse)%wseadj)
      else
        call bnd_wse(TH_str(iwse)%ncells,TH_str(iwse)%cells,TH_str(iwse)%faces,TH_str(iwse)%wsebnd)  
      endif
    enddo !iwse
    
!--- Single Water Surface Elevation BC (Type 3=H) -------------------------
!Water level is first assigend at the center cell of the cell string
!and then adjusted for wind and wave setup starting from the center cell      
    do iwse=1,nHstr        
      if(H_str(iwse)%wseadjust)then
        call bnd_wse_adjust(H_str(iwse)%ncells,H_str(iwse)%cells,H_str(iwse)%faces,H_str(iwse)%wsebnd,H_str(iwse)%wseadj)
        call bnd_wse(H_str(iwse)%ncells,H_str(iwse)%cells,H_str(iwse)%faces,H_str(iwse)%wseadj)
      else
        call bnd_wse(H_str(iwse)%ncells,H_str(iwse)%cells,H_str(iwse)%faces,H_str(iwse)%wsebnd)
      endif  
    enddo !iwse 

!--- Multiple Water Level BC (Type 4=MH) ---------------------------
    do iwse=1,nMHstr
      call bnd_wse(MH_str(iwse)%ncells,MH_str(iwse)%cells,MH_str(iwse)%faces,MH_str(iwse)%wsebnd)
    enddo !iwse
    
!--- Multiple Water Level and Velocity BC (Type 5=MHV) -----------------------
    do iwse=1,nMHVstr
      call bnd_wse(MHV_str(iwse)%ncells,MHV_str(iwse)%cells,MHV_str(iwse)%faces,MHV_str(iwse)%wsebnd)
    enddo !iwse        

!!--- Cross-shore BC (Type 6=CS) ------------------------------------
!    call xshore_wse
!    do icsh=1,nCSstr
!    !  call bnd_wse(CS_str(icsh)%ncells,CS_str(icsh)%cells,CS_str(icsh)%faces,CS_str(icsh)%wsecsh)
!      do j=1,CS_str(icsh)%ncells
!        i=CS_str(icsh)%cells(j)
!        k=CS_str(icsh)%faces(j)
!        nck=cell2cell(k,i)        
!        if(flux(k,i)>0.0)then !Outflow
!          pp(nck)=CS_str(icsh)%wsecsh(j)*grav-p(nck)   !Set given wse at the dummy node  
!          if(ncface(i)==4)then
!            acoefik=acoef(kkface(k),i)
!          else
!            acoefik=0.0; kkdf=kkface(idirface(k,i))
!            do kk=1,ncface(i)    !Assume boundary face does not split
!              if(idirface(kk,i)==kkdf) acoefik=acoefik+acoef(kk,i)
!            enddo                   
!          endif                        
!          su(i)=su(i)+acoefik*pp(nck)  !using a_kk because a_k is not well defined  
!          sp(i)=sp(i)-acoefik
!        endif
!        acoef(k,i)=0.0
!      enddo !j  
!    enddo
    
!--- Nested Water Level BC (Type 7=NH) ------------------------------
    do iwse=1,nNHstr
      call bnd_wse(NH_str(iwse)%ncells,NH_str(iwse)%cells,NH_str(iwse)%faces,NH_str(iwse)%wsebnd)
    enddo !iwse
    
!--- Nested Water Level and Velocity BC (Type 8=NHV) --------------------
    do iwse=1,nNHVstr
      call bnd_wse(NHV_str(iwse)%ncells,NHV_str(iwse)%cells,NHV_str(iwse)%faces,NHV_str(iwse)%wsebnd)  
    enddo !iwse

!--- Tidal Database Water Level BC (Type 9=NTH) ------------------------------
    do iwse=1,nNTHstr
      if(NTH_str(iwse)%wseadjust)then  
        call bnd_wse_adjust(NTH_str(iwse)%ncells,NTH_str(iwse)%cells,NTH_str(iwse)%faces,NTH_str(iwse)%wsebnd,NTH_str(iwse)%wseadj)
        call bnd_wse(NTH_str(iwse)%ncells,NTH_str(iwse)%cells,NTH_str(iwse)%faces,NTH_str(iwse)%wseadj)
      else
        call bnd_wse(NTH_str(iwse)%ncells,NTH_str(iwse)%cells,NTH_str(iwse)%faces,NTH_str(iwse)%wsebnd)  
      endif
    enddo !iwse
    
!--- Tidal Database Water Level and Velocity BC (Type 10=NTHV) --------------------
    do iwse=1,nNTHVstr
      if(NTHV_str(iwse)%wseadjust)then  
        call bnd_wse_adjust(NTHV_str(iwse)%ncells,NTHV_str(iwse)%cells,NTHV_str(iwse)%faces,NTHV_str(iwse)%wsebnd,NTHV_str(iwse)%wseadj)
        call bnd_wse(NTHV_str(iwse)%ncells,NTHV_str(iwse)%cells,NTHV_str(iwse)%faces,NTHV_str(iwse)%wseadj)
      else
        call bnd_wse(NTHV_str(iwse)%ncells,NTHV_str(iwse)%cells,NTHV_str(iwse)%faces,NTHV_str(iwse)%wsebnd)  
      endif
    enddo !iwse    
    
!--- Structures -------------------    
    !call struct_p

    return
    end subroutine bound_wse  
    
!*******************************************************************
    function wsesetupx(i,i1)
! Calculates local setup correction due to wind and wave forcing in x
!*******************************************************************
    use flow_def, only: h,gravinv
    use geo_def, only: x
    use met_def, only: windconst,windvar,tauwx,tauwindx
    use prec_def, only: ikind
    
    implicit none
    integer,intent(in) :: i,i1
    real(ikind):: wsesetupx
    
    interface
      function forcex(i)
        use prec_def
        integer,intent(in) :: i
        real(ikind) :: forcex
      end function
    endinterface
    
    wsesetupx=(x(i)-x(i1))*(forcex(i)+forcex(i1))/(h(i)+h(i1))*gravinv
    
    return
    end function wsesetupx
    
!*******************************************************************
    function wsesetupy(i,i1)
! Calculates local setup correction due to wind and wave forcing in y
!*******************************************************************   
    use flow_def, only: h,gravinv
    use geo_def, only: y
    use prec_def, only: ikind
    
    implicit none 
    integer,intent(in) :: i,i1
    real(ikind) :: wsesetupy
    
    interface
      function forcey(i)
        use prec_def
        integer,intent(in) :: i
        real(ikind) :: forcey
      end function
    endinterface
    
    wsesetupy=(y(i)-y(i1))*(forcey(i)+forcey(i1))/(h(i)+h(i1))*gravinv
       
    return
    end function wsesetupy
    
!*******************************************************************
    function forcex(i) result(fx)
! Calculates the sum of the wave and wind forces at cell i in the 
! x direction. Used at boundary cells
! written by Alex Sanchez, USACE-CHL
!*******************************************************************
    use cms_def, only: noptset
    use comvarbl, only: ramp
    use der_def, only: gow
    use der_lib, only: dx2d,dy2d
    use flow_def, only: u,us,h,rhow,waveflux
    use fric_def, only: bbl_stream,z0,cbcfuwcap
    use fric_lib, only: fric_streaming_stress
    use geo_def, only: areap
    use met_def, only: windconst,windvar,tauwx,tauwindx,&
        iwndlagr,cdWndareap,wndx,uwind,presvar,pressatmdx
    use prec_def, only: ikind
    use q3d_def, only: q3d,q3d_to_flow,f3dxx,f3dxy,f3dyy
    use wave_flowgrid_def, only: wavestrx,worbrep,wper,wlen,wunitx
    implicit none
    integer,intent(in) :: i
    real(ikind) :: fx,taustr,za,dvarx,dvary
    
    fx = 0.0
    if(noptset>=3)then
      !Wave and roller stresses  
      fx = fx + wavestrx(i)
      !Bottom boundary layer streaming (beta)
      if(bbl_stream)then
        !za = 1.6667e-05 !=0.2/1000/12
        za = z0(i)
        taustr = fric_streaming_stress(rhow,za,worbrep(i),wper(i),wlen(i)) 
        fx = fx+taustr*wunitx(i)/rhow 
      endif
      !Dispersion and wave-current interaction terms (beta)
      if(q3d .and. q3d_to_flow)then
        call dx2d(gow,i,f3dxx,dvarx)         
        call dy2d(gow,i,f3dxy,dvary)
        fx = fx - (dvarx+dvary)*ramp
      endif
      !Wave velocity forcing
      if(waveflux) fx = fx + cbcfuwcap(i)*us(i)/areap(i) !Note: cbcfuwcap may contain bed slope term
    endif

    !Wind
    if(iwndlagr==0)then !Eulerian
      if(windconst) fx = fx + tauwx
      if(windvar)   fx = fx + tauwindx(i)
    else !Lagrangian
      if(windconst) fx = fx + cdWndareap(i)*(wndx-u(i)+us(i))/areap(i)
      if(windvar)   fx = fx + cdWndareap(i)*(uwind(i)-u(i)+us(i))/areap(i)
    endif
    
    !Atmospheric pressure gradients
    if(presvar) fx = fx - pressatmdx(i)*h(i)/rhow
    
    
    return
    end function forcex
    
!*******************************************************************
    function forcey(i) result(fy)
! Calculates the sum of the wave and wind forces at cell i in the 
! y direction. Used at boundary cells
! written by Alex Sanchez, USACE-CHL
!*******************************************************************   
    use cms_def, only: noptset
    use comvarbl, only: ramp
    use der_lib, only: dx2d,dy2d
    use der_def, only: gow
    use flow_def, only: v,vs,h,rhow,waveflux
    use fric_def, only: bbl_stream,z0,cbcfuwcap
    use fric_lib, only: fric_streaming_stress
    use geo_def, only: areap
    use met_def, only: windconst,windvar,tauwy,tauwindy,&
        iwndlagr,cdWndareap,wndy,vwind,presvar,pressatmdy
    use prec_def, only: ikind
    use q3d_def, only: q3d,q3d_to_flow,f3dxx,f3dxy,f3dyy
    use wave_flowgrid_def, only: wavestry,worbrep,wper,wlen,wunity
    
    implicit none 
    integer,intent(in) :: i
    real(ikind) :: fy,za,taustr,dvarx,dvary
    
    fy = 0.0
    !Wave forcing
    if(noptset>=3)then
      !Wave and roller stresses
      fy = fy + wavestry(i)
      !Bottom boundary layer streaming (beta)
      if(bbl_stream)then
        !za = 1.6667e-05 !=0.2/1000/12
        za = z0(i)
        taustr = fric_streaming_stress(rhow,za,worbrep(i),wper(i),wlen(i)) 
        fy = fy + taustr*wunity(i)/rhow 
      endif
      !Dispersion and wave-current interaction terms (beta)
      if(q3d .and. q3d_to_flow)then
        call dx2d(gow,i,f3dxy,dvarx)
        call dy2d(gow,i,f3dyy,dvary)  
        fy = fy - (dvarx+dvary)*ramp
      endif 
      !Wave velocity forcing
      if(waveflux) fy = fy + cbcfuwcap(i)*vs(i)/areap(i) !Note: cbcfuwcap may contain bed slope term
    endif  
    
    !Wind forcing
    if(iwndlagr==0)then !Eulerian
      if(windconst) fy = fy + tauwy
      if(windvar)   fy = fy + tauwindy(i)
    else !Lagrangian
      if(windconst) fy = fy + cdWndareap(i)*(wndy-v(i)+vs(i))/areap(i)
      if(windvar)   fy = fy + cdWndareap(i)*(vwind(i)-v(i)+vs(i))/areap(i)
    endif
    
    !Atmospheric pressure gradients
    if(presvar) fy = fy - pressatmdy(i)*h(i)/rhow
    
    return
    end function forcey

!***********************************************************************
    subroutine bnd_eval()
!  Updates the flux and water level values at boundary cells
!  written by Alex Sanchez, USACE-ERDC-CHL        
!***********************************************************************
    use comvarbl, only: timehrs,ramp,dtj
    use geo_def, only: azimuth_fl
    use nest_lib, only: nest_wse_eval,nest_vel_eval
    use flow_def, only: h
    use bnd_def
    use plagr_lib
    use prec_def, only: ikind
    
    implicit none
    integer :: j,k,iriv,iwse,np,ipar,idpar
    real(ikind) :: wsebnd,stagebnd !,wsemax,wsemin
    integer,parameter :: nb = 4 !>=nti+1
    real(ikind) :: lb(nb),lb2(nb)
    
    lb = 0.0
    lb2 = 0.0

!--- Flux BC (Type 1=Q) -----------------------------------------------------
    do iriv=1,nQstr  !for each cell string
      if(Q_str(iriv)%ifluxmode==1)then !Constant
        Q_str(iriv)%qflux = Q_str(iriv)%qfluxconst
      elseif(Q_str(iriv)%ifluxmode==2)then !Time-series
        call plagr_fit(Q_str(iriv)%ntimes,Q_str(iriv)%times,timehrs,nb,lb,Q_str(iriv)%nti,np,Q_str(iriv)%inc)  
        !Q_STR(iriv)%qflux = plagr_eval(Q_STR(iriv)%ntimes,Q_STR(iriv)%qcurv,nb,lb,np,Q_STR(iriv)%inc)
        Q_str(iriv)%qflux = sum(lb(1:np+1)*Q_str(iriv)%qcurv(Q_str(iriv)%inc:Q_str(iriv)%inc+np))
      elseif(Q_str(iriv)%ifluxmode==3)then !Stage-Flow curve
        stagebnd = 0.0
        do j=1,Q_str(iriv)%ncells
          stagebnd = stagebnd + h(Q_str(iriv)%cells(j))
        enddo
        stagebnd = stagebnd/real(Q_str(iriv)%ncells,kind=ikind)
        k=0
        call plagr_fit(Q_str(iriv)%nstages,Q_str(iriv)%stage,stagebnd,nb,lb,1,np,k)
        Q_str(iriv)%qflux = plagr_eval(Q_str(iriv)%nstages,Q_str(iriv)%rflow,nb,lb,np,k)
      endif
      Q_str(iriv)%qflux = ramp*Q_str(iriv)%qflux  !Initial condition is set to zero
    enddo
    
!--- Tidal/Harmonic BC (Type 2=T) ---------------------------------------------------
    do iwse=1,nTHstr
      do j=1,TH_str(iwse)%ncells
        wsebnd = TH_str(iwse)%wseoffset
        if(TH_str(iwse)%istidal)then
          do k=1,TH_str(iwse)%ntc  
            wsebnd = wsebnd + TH_str(iwse)%f(k)*TH_str(iwse)%amp(k) &
                *cos(TH_str(iwse)%speed(k)*(timehrs+dtj) + TH_str(iwse)%vu(k) &
                - TH_str(iwse)%phase(k) + TH_str(iwse)%psi(j,k))
          enddo
        else
          do k=1,TH_str(iwse)%ntc  
            wsebnd = wsebnd + TH_str(iwse)%amp(k)*cos(TH_str(iwse)%speed(k)*timehrs &
                - TH_STR(iwse)%phase(k) + TH_str(iwse)%psi(j,k))
          enddo
        endif
 !       write(3000,*)'ioffsetmode (bnd_eval) = ',TH_str(iwse)%ioffsetmode 
 !       write(3000,*)'TH_str(iwse)%wseoffset = ',TH_str(iwse)%wseoffset 
        if(TH_str(iwse)%ioffsetmode==1)then !Constant (hli,10/05/17)
           wsebnd = wsebnd + TH_str(iwse)%wseoffset !Add Offset 
        elseif(TH_str(iwse)%ioffsetmode==2)then !Time-series
           call plagr_fit(TH_str(iwse)%ntimesoffset,TH_str(iwse)%offsettimes,timehrs,nb,lb,TH_str(iwse)%nti,np,TH_str(iwse)%inc)  !hli(10/06/17)
           TH_str(iwse)%wsecurveoffset = sum(lb(1:np+1)*TH_str(iwse)%offsetcurve(TH_str(iwse)%inc:TH_str(iwse)%inc+np))           !hli(10/06/17)
           wsebnd = wsebnd + TH_str(iwse)%wsecurveoffset
        endif
        TH_str(iwse)%wsebnd(j) = (1.0-ramp)*TH_str(iwse)%wsebnd0(j) &
                      + ramp*(wsebnd + TH_str(iwse)%wsevar(j))
!        write(3000,*)'j = ',j
!        write(3000,*)'TH_str(iwse)%wsebnd(j) = ',TH_str(iwse)%wsebnd(j)
      enddo !j-cell
    enddo !iwse-str

!--- Single Water Level BC (Type 3=H) ------------------------------------------
    do iwse=1,nHstr  !for each cell string
      if(H_STR(iwse)%ntimes>0)then !Time series
        call plagr_fit(H_str(iwse)%ntimes,H_str(iwse)%times,timehrs,nb,lb,H_str(iwse)%nti,np,H_str(iwse)%inc)  
        !wsebnd = plagr_eval(H_STR(iwse)%ntimes,H_STR(iwse)%wsecurv,nb,lb,np,H_STR(iwse)%inc)
        wsebnd = sum(lb(1:np+1)*H_str(iwse)%wsecurv(H_str(iwse)%inc:H_str(iwse)%inc+np))
        !!Limit value to be within neighboring values (higher order interpolations can produce extremal values)
        !wsemax = maxval(H_str(iwse)%wsecurv(H_str(iwse)%inc:H_str(iwse)%inc+np))
        !wsemin = minval(H_str(iwse)%wsecurv(H_str(iwse)%inc:H_str(iwse)%inc+np))
        !wsebnd = max(min(wsebnd,wsemax),wsemin)
      else !Constant
        wsebnd = H_str(iwse)%wseconst
      endif
!!745   format(2(F8.3,1x),I4,I3,5F12.4) !for testing
!!      write(23,745) timehrs,wsebnd,H_str(iwse)%inc,np,lb(1:np+1) !for testing
      !wsebnd = wsebnd + H_str(iwse)%wseoffset !Add Offset
      if(H_str(iwse)%ioffsetmode==1)then !Constant (hli,01/19/17)
      wsebnd = wsebnd + H_str(iwse)%wseoffset !Add Offset
      elseif(H_str(iwse)%ioffsetmode==2)then !Time-series
!        write(3000,*)'ntimesoffset = ',H_str(iwse)%ntimesoffset,'offsettimes =',H_str(iwse)%offsettimes 
!        call SINTER(H_str(iwse)%offsettimes,H_str(iwse)%offsetcurve,H_str(iwse)%times,H_str(iwse)%wsecurveoffset,H_str(iwse)%ntimesoffset,H_str(iwse)%ntimes)  
        call plagr_fit(H_str(iwse)%ntimesoffset,H_str(iwse)%offsettimes,timehrs,nb,lb,H_str(iwse)%nti,np,H_str(iwse)%inc)  
!        write(3000,*)'iwse =',iwse,'H_str(iwse)%nti = ',H_str(iwse)%nti,'H_str(iwse)%inc =',H_str(iwse)%inc
        H_str(iwse)%wsecurveoffset = sum(lb(1:np+1)*H_str(iwse)%offsetcurve(H_str(iwse)%inc:H_str(iwse)%inc+np))
        wsebnd = wsebnd + H_str(iwse)%wsecurveoffset
      endif
!      wsebnd = wsebnd + H_str(iwse)%wseoffset !Add Offset
      do j=1,H_str(iwse)%ncells
        H_str(iwse)%wsebnd(j) = (1.0-ramp)*H_str(iwse)%wsebnd0(j)+ramp*(wsebnd + H_str(iwse)%wsevar(j))    
      enddo
!      write(3000,*)'ioffsetmode (bnd_eval) = ',H_str(iwse)%ioffsetmode 
    enddo ! end of each cell string
    
!--- Multiple Water Level BC (Type 4=MH) ----------------------------------------
    do iwse=1,nMHstr  !for each cell string
      !Temporal interpolations  
      call plagr_fit(MH_str(iwse)%ntimes,MH_str(iwse)%times,timehrs,nb,lb,MH_str(iwse)%nti,np,MH_str(iwse)%inc)
      do j=1,MH_str(iwse)%ncells   !for each cell in string
        MH_str(iwse)%wsebnd(j) = sum(lb(1:np+1)*MH_str(iwse)%wsedata(MH_str(iwse)%inc:MH_str(iwse)%inc+np,j))
        MH_str(iwse)%wsebnd(j) = MH_str(iwse)%wsebnd(j) + MH_str(iwse)%wseoffset !Add offset here
        MH_str(iwse)%wsebnd(j) = ramp*MH_str(iwse)%wsebnd(j) + (1.0-ramp)*MH_str(iwse)%wsebnd0(j) !Ramp
      enddo !j
    enddo !iwse
    
!--- Multiple Water Level and Velocity BC (Type 5=MHV) --------------------------    
    do iwse=1,nMHVstr  !for each cell string
      !Temporal interpolations
      call plagr_fit(MHV_str(iwse)%ntimeswse,MHV_str(iwse)%timeswse,timehrs,nb,lb ,MHV_str(iwse)%ntivel,np,MHV_str(iwse)%incwse)
      call plagr_fit(MHV_str(iwse)%ntimesvel,MHV_str(iwse)%timesvel,timehrs,nb,lb2,MHV_str(iwse)%ntivel,np,MHV_str(iwse)%incvel)
      do j=1,MHV_str(iwse)%ncells   !for each cell in string
        MHV_str(iwse)%wsebnd(j) = sum(lb(1:np+1)*MHV_str(iwse)%wsedata(MHV_str(iwse)%incwse:MHV_str(iwse)%incwse+np,j))
        MHV_str(iwse)%ubnd(j) = sum(lb (1:np+1)*MHV_str(iwse)%udata(MHV_str(iwse)%incvel:MHV_str(iwse)%incvel+np,j))
        MHV_str(iwse)%vbnd(j) = sum(lb2(1:np+1)*MHV_str(iwse)%vdata(MHV_str(iwse)%incvel:MHV_str(iwse)%incvel+np,j))
        MHV_str(iwse)%wsebnd(j) = MHV_str(iwse)%wsebnd(j) + MHV_STR(iwse)%wseoffset !Add offset here
        MHV_str(iwse)%wsebnd(j) = ramp*MHV_str(iwse)%wsebnd(j) + (1.0-ramp)*MHV_str(iwse)%wsebnd0(j)
        MHV_str(iwse)%ubnd(j) = ramp*MHV_str(iwse)%ubnd(j) + (1.0-ramp)*MHV_str(iwse)%ubnd0(j)
        MHV_str(iwse)%vbnd(j) = ramp*MHV_str(iwse)%vbnd(j) + (1.0-ramp)*MHV_str(iwse)%vbnd0(j) 
      enddo
    enddo ! end of each cell string
    
!--- Parent Solutions --------------------------------------------------------------------
    do ipar=1,nParSim
      call par_wse_eval(ParSim(ipar)%wsefilepar,ParSim(ipar)%wsepathpar,ParSim(ipar)%nptspar,&
        ParSim(ipar)%t2hrs,ParSim(ipar)%timestarthr,ParSim(ipar)%ntiwsepar,&
        ParSim(ipar)%incwsepar,ParSim(ipar)%timewsehrspar,ParSim(ipar)%wsepar)  
      if(ParSim(ipar)%velpar)then
        call par_vel_eval(ParSim(ipar)%velfilepar,ParSim(ipar)%velpathpar,ParSim(ipar)%nptspar,&
          ParSim(ipar)%t2hrs,ParSim(ipar)%timestarthr,ParSim(ipar)%ntivelpar,&
          ParSim(ipar)%incvelpar,ParSim(ipar)%timevelhrspar,ParSim(ipar)%upar,ParSim(ipar)%vpar)  
      endif  
    enddo
    
!--- Nested Water Level BC (Type 7=NH) --------------------------    
    do iwse=1,nNHstr  !for each cell string
      idpar=NH_str(iwse)%idpar
      call nest_wse_eval(ParSim(idpar)%ntiwsepar,&
        ParSim(idpar)%timewsehrspar,ParSim(idpar)%nptspar,ParSim(idpar)%wsepar,&
        NH_str(iwse)%ncells,NH_str(iwse)%cells,NH_str(iwse)%mntp,NH_str(iwse)%intp,NH_str(iwse)%cntp,&
        NH_str(iwse)%ntiwse,NH_str(iwse)%timewsehrs,NH_str(iwse)%wsedata,NH_str(iwse)%wsebnd,&
        NH_str(iwse)%wseout,NH_str(iwse)%wsefile)
      NH_str(iwse)%wsebnd = NH_str(iwse)%wsebnd + NH_str(iwse)%wseoffset !Add offset here
      NH_str(iwse)%wsebnd = (1.0-ramp)*NH_str(iwse)%wsebnd0 + ramp*NH_str(iwse)%wsebnd   
    enddo
    
!--- Nested Water Level and Velocity BC (Type 8=NHV) --------------------------    
    do iwse=1,nNHVstr  !for each cell string
      idpar=NHV_str(iwse)%idpar
      call nest_wse_eval(ParSim(idpar)%ntiwsepar,&
        ParSim(idpar)%timewsehrspar,ParSim(idpar)%nptspar,ParSim(idpar)%wsepar,&
        NHV_str(iwse)%ncells,NHV_str(iwse)%cells,NHV_str(iwse)%mntp,NHV_str(iwse)%intp,NHV_str(iwse)%cntp,&
        NHV_str(iwse)%ntiwse,NHV_str(iwse)%timewsehrs,NHV_str(iwse)%wsedata,NHV_str(iwse)%wsebnd,&
        NHV_str(iwse)%wseout,NHV_str(iwse)%wsefile)  
      NHV_str(iwse)%wsebnd = NHV_str(iwse)%wsebnd + NHV_str(iwse)%wseoffset !Add offset here
      NHV_str(iwse)%wsebnd = (1.0-ramp)*NHV_str(iwse)%wsebnd0 + ramp*NHV_str(iwse)%wsebnd
      call nest_vel_eval(ParSim(idpar)%ntivelpar,&
        ParSim(idpar)%timevelhrspar,ParSim(idpar)%nptspar,ParSim(idpar)%upar,ParSim(idpar)%vpar,&
        NHV_str(iwse)%ncells,NHV_str(iwse)%cells,NHV_str(iwse)%mntp,NHV_str(iwse)%intp,NHV_str(iwse)%cntp,&
        NHV_str(iwse)%ntivel,NHV_str(iwse)%timevelhrs,NHV_str(iwse)%udata,NHV_str(iwse)%vdata,&
        NHV_str(iwse)%ubnd,NHV_str(iwse)%vbnd,NHV_str(iwse)%velout,NHV_str(iwse)%velfile)  
      NHV_str(iwse)%ubnd = (1.0-ramp)*NHV_STR(iwse)%ubnd0 + ramp*NHV_str(iwse)%ubnd
      NHV_str(iwse)%vbnd = (1.0-ramp)*NHV_STR(iwse)%vbnd0 + ramp*NHV_str(iwse)%vbnd              
    enddo
    
!--- Nested Tidal Database Water Level BC (Type 9-NTH) --------------------------------    
    do iwse=1,nNTHstr  !for each cell string
      NTH_str(iwse)%wsebnd = NTH_str(iwse)%wseoffset  
      do j=1,NTH_str(iwse)%ncells
        do k=1,NTH_str(iwse)%ntc
          NTH_str(iwse)%wsebnd(j) = NTH_str(iwse)%wsebnd(j) + NTH_str(iwse)%f(k)*NTH_str(iwse)%amp(j,k) &
             *cos(NTH_str(iwse)%speed(k)*(timehrs+dtj) + NTH_str(iwse)%vu(k) - NTH_str(iwse)%phase(j,k))
        enddo        
        NTH_str(iwse)%wsebnd(j) = (1.0-ramp)*NTH_str(iwse)%wsebnd0(j) + ramp*NTH_str(iwse)%wsebnd(j)
      enddo !j-cell
    enddo
    
!--- Nested Tidal Database WSE and Velocity BC (Type 10-NTHV) --------------------------------    
    do iwse=1,nNTHVstr  !for each cell string
      NTHV_str(iwse)%wsebnd = NTHV_str(iwse)%wseoffset
      NTHV_str(iwse)%ubnd = 0.0; NTHV_str(iwse)%vbnd = 0.0
      do j=1,NTHV_STR(iwse)%ncells
        do k=1,NTHV_str(iwse)%ntc  
          NTHV_str(iwse)%wsebnd(j) = NTHV_str(iwse)%wsebnd(j) + NTHV_str(iwse)%f(k)*NTHV_str(iwse)%amp(j,k) &
                   *cos(NTHV_STR(iwse)%speed(k)*(timehrs+dtj) + NTHV_str(iwse)%vu(k) - NTHV_str(iwse)%phase(j,k))
          NTHV_str(iwse)%ubnd(j) =   NTHV_str(iwse)%ubnd(j)   + NTHV_str(iwse)%f(k)*NTHV_str(iwse)%ampu(j,k) &
                   *cos(NTHV_STR(iwse)%speed(k)*(timehrs+dtj) + NTHV_str(iwse)%vu(k) - NTHV_str(iwse)%phaseu(j,k))
          NTHV_str(iwse)%vbnd(j) =   NTHV_str(iwse)%vbnd(j)   + NTHV_str(iwse)%f(k)*NTHV_str(iwse)%ampv(j,k) &
                   *cos(NTHV_STR(iwse)%speed(k)*(timehrs+dtj) + NTHV_str(iwse)%vu(k) - NTHV_str(iwse)%phasev(j,k))
        enddo
      enddo !j-cell
      call rotate_vector(NTHV_str(iwse)%ncells,NTHV_str(iwse)%ncells,&
            NTHV_str(iwse)%angvel,NTHV_str(iwse)%ubnd,NTHV_str(iwse)%vbnd)
      NTHV_str(iwse)%wsebnd = (1.0-ramp)*NTHV_str(iwse)%wsebnd0 + ramp*NTHV_str(iwse)%wsebnd
      NTHV_str(iwse)%ubnd   = (1.0-ramp)*NTHV_str(iwse)%ubnd0   + ramp*NTHV_str(iwse)%ubnd
      NTHV_str(iwse)%vbnd   = (1.0-ramp)*NTHV_str(iwse)%vbnd0   + ramp*NTHV_str(iwse)%vbnd
    enddo    

    return
    end subroutine bnd_eval
    
!***********************************************************************
    subroutine bndriver_vel()
! River in/outflow boundaries
!
! Variables:
!   qflux = total discharge across the boundary
!   coefman = manning's coefficient
!   h = total water depth
!   cmvel = conveyance coefficient
!
! written by Weiming Wu, NCCHE, Oct. 2008
! modified by Alex Sanchez, USACE-CHL
! last modified 02/14/2013
!***********************************************************************
#include "CMS_cpp.h"
    use size_def, only: ncells
    use geo_def, only: cell2cell,llec2llec,ds,fnx,fny,ncface,mapid
    use flow_def, only: iwet,hk,u,v,flux
    use bnd_def
    use comvarbl, only: ramp,niter
    use fric_def, only: coefman
    use const_def, only: deg2rad,small
    use diag_def
    use diag_lib, only: diag_print_error, diag_print_warning
    use prec_def, only: ikind
    
    implicit none
    integer :: i,j,k,iriv,nck,jcn
    real(ikind) :: sumw,val,cosang,sinang,qnorm,qfluxchk

    !Check fluxes
#ifdef DIAG_MODE
    do i=1,ncells
      do k=1,ncface(i)
        nck=cell2cell(k,i)  
        if(iwet(i)*iwet(nck)==0 .and. abs(flux(k,i))>1.0e-15)then
          write(msg2,*) '  Cell: ',mapid(i)
          write(msg3,*) '  Neighbor: ',mapid(nck)
          call diag_print_error('Non-zero flux between wet and dry cells before bndriver_vel',msg2,msg3)
        endif
      enddo
    enddo
#endif
    
    do iriv=1,nQstr
      sumw = 0.0 !Sum Weights
      cosang = cos(Q_str(iriv)%angle*deg2rad)
      sinang = sin(Q_str(iriv)%angle*deg2rad)
      !Compute weights and temporarily store in velocities at ghost cells
      do j=1,Q_str(iriv)%ncells
        i = Q_str(iriv)%cells(j)
        k = Q_str(iriv)%faces(j)
        nck = cell2cell(k,i)
        !!iwet(nck)=iwet(i)
        val = iwet(i)*iwet(nck)*hk(k,i)**Q_str(iriv)%cmvel/max(coefman(i),1.0e-5) !Note: h.^(1+r) is split to h^r*h, the h is below
        u(nck) = val*cosang !temporary storage of weights
        v(nck) = val*sinang !temporary storage of weights
        sumw = sumw + ds(k,i)*hk(k,i)*(u(nck)*fnx(k,i)+v(nck)*fny(k,i)) !Note: h.^(1+r) is split to h^r*h, the h is here
      enddo
      qnorm = Q_str(iriv)%qflux/max(abs(sumw),1.0e-6) !Note: sumq cannot be zero
      !Distribute total flux across boundary
      qfluxchk = 0.0
      do j=1,Q_str(iriv)%ncells
        i = Q_str(iriv)%cells(j)
        k = Q_str(iriv)%faces(j)
        nck = cell2cell(k,i)
        jcn = llec2llec(k,i)
        u(nck) = u(nck)*qnorm
        v(nck) = v(nck)*qnorm
        flux(k,i) = iwet(i)*iwet(nck)*ds(k,i)*hk(k,i)*(u(nck)*fnx(k,i)+v(nck)*fny(k,i))

        flux(jcn,nck) = -flux(k,i)
        qfluxchk = qfluxchk - flux(k,i)
      enddo
      if (Q_str(iriv)%qflux .ne. 0.0) then                      !Avoid divide by zero  MEB  04/08/2021
        val=abs((qfluxchk-Q_str(iriv)%qflux)/Q_str(iriv)%qflux) !Normalized error
        if(val>0.001)then                                      
          write(msg2,*) '  Flux Boundary:   ',Q_str(iriv)%idnum
          write(msg3,*) '  Specified Flux:  ',Q_str(iriv)%qflux,' m^3/s'
          write(msg4,*) '  Calculated Flux: ',qfluxchk,' m^3/s'  
          call diag_print_warning('Problem distributing flow at flux boundary ',msg2,msg3,msg4)
        endif
      endif
    enddo
    
    !Check fluxes
#ifdef DIAG_MODE
    do i=1,ncells
      do k=1,ncface(i)
        if(iwet(i)*iwet(cell2cell(k,i))==0 .and. abs(flux(k,i))>1.0e-15)then
          write(msg2,*) '  Cell: ',i
          write(msg3,*) '  Neighbor: ',cell2cell(k,i)
          call diag_print_error('Non-zero flux between wet and dry cells after bndriver_vel',msg2,msg3)
        endif
      enddo
    enddo
#endif

    return
    end subroutine bndriver_vel
    
!***********************************************************************
    subroutine bndflux
! Updatse volume fluxes at all the boundaries
! written by Weiming Wu, NCCHE, Oct. 2008
! modified by Alex Sanchez, USACE-CHL
!***********************************************************************
#include "CMS_cpp.h"
#ifdef DIAG_MODE
    use size_def, only: ncells
    use geo_def, only: cell2cell,ncface
    use flow_def, only: flux,iwet
#endif
    use bnd_def
    use diag_lib, only: diag_print_error

!!#ifdef DEV_MODE
!!    use comvarbl, only: dtime
!!    use flow_def, only: h,h1
!!    use geo_def, only: areap
!!#endif
    
    implicit none
    integer :: iriv,icsh
#ifdef DIAG_MODE
    integer :: i,k
#endif
!!#ifdef DEV_MODE
!!    integer :: nck
!!#endif

    !Check fluxes
#ifdef DIAG_MODE 
    do i=1,ncells
      do k=1,ncface(i)
        if(iwet(i)*iwet(cell2cell(k,i))==0 .and. abs(flux(k,i))>1.0e-15)then
          call diag_print_error('Non-zero flux between wet and dry cells')
        endif
      enddo
    enddo
#endif

!!#ifdef DEV_MODE
!!    call flow_negdepth
!!#endif

!--- Flux -------------------------------------
    do iriv=1,nQstr
      call fluxbnd(Q_str(iriv)%ncells,Q_str(iriv)%cells,Q_str(iriv)%faces)      
    enddo
    
!!--- Multiple WSE-Vel ------------------------------------    
!    do iwse=1,nMHVstr
!      call fluxbnd(MHV_str(iwse)%ncells,MHV_str(iwse)%cells,MHV_str(iwse)%faces)  
!    enddo
!    
!!--- Nested WSE-Vel ------------------------------------    
!    do iwse=1,nNHVstr
!      call fluxbnd(NHV_str(iwse)%ncells,NHV_str(iwse)%cells,NHV_str(iwse)%faces)  
!    enddo
!    
!!--- Tidal WSE-Vel ------------------------------------    
!    do iwse=1,nNTHVstr
!      call fluxbnd(NTHV_str(iwse)%ncells,NTHV_str(iwse)%cells,NTHV_str(iwse)%faces)  
!    enddo        
    
!--- Cross-shore ------------------------------------
    do icsh=1,nCSstr
      call fluxbnd(CS_str(icsh)%ncells,CS_str(icsh)%cells,CS_str(icsh)%faces)  
    enddo

!--- Structures ----
    call struct_flux

!!#ifdef DEV_MODE
!!!--- Open Boundaries ----------------------
!!! Computes fluxes at boundaries using the continuity equation
!!    do i=1,ncells
!!      if(iwet(i)==0) cycle
!!      do k=1,ncface(i)
!!        nck=cell2cell(k,i)  
!!        if(iwet(i)*iwet(nck)==0 .or. nck<=ncells) cycle
!!        flux(k,i) = h1(i) - h(i) - dtime*(sum(flux(1:ncface(i),i))-flux(k,i))/areap(i)   
!!      enddo
!!    enddo
!!#endif

!--- Dry nodes ------------    
!Note: dry faces already treated in coefsourcesink_pp but check here anyway
!#ifdef DIAG_MODE
!    do i=1,ncells
!      flux(1:ncface(i),i)=iwet(i)*iwet(cell2cell(1:ncface(i),i))*flux(1:ncface(i),i)
!    enddo
!#endif
    !Check fluxes
#ifdef DIAG_MODE
    do i=1,ncells
      do k=1,ncface(i)
        if(iwet(i)*iwet(cell2cell(k,i))==0 .and. abs(flux(k,i))>1.0e-15)then
          call diag_print_error('Non-zero flux between wet and dry cells')
        endif
      enddo
    enddo
#endif

    return
    end subroutine bndflux
    
!**********************************************************************
    subroutine bndcopy2ghost(val)
!**********************************************************************
    use size_def, only: ncellsD
    use geo_def, only: cell2cell
    use bnd_def
    use prec_def, only: ikind
    
    implicit none
    integer :: i,j,k,ibnd
    real(ikind) :: val(ncellsD)

    do ibnd=1,nbndstr
      do j=1,bnd_str(ibnd)%ncells
        i=bnd_str(ibnd)%cells(j)
        k=bnd_str(ibnd)%faces(j)
        val(cell2cell(k,i))=val(i)
      enddo
    enddo  

    return
    end subroutine bndcopy2ghost

!**********************************************************************
    subroutine xshore_uv
! Solves the cross-shore momentum equations for the current velocities
! written by Weiming Wu
! modified by Alex Sanchez, USACE-CHL
!**********************************************************************
    use geo_def, only: cell2cell,idirface,dx,dy
    use flow_def, only: iwet,h,u,v,us,vs,vis
    use bnd_def, only: nCSstr,CS_str
    use fric_def, only: cfrict,z0
    use wave_flowgrid_def, only: worb,worbrep,wper,wang    
    use comvarbl, only: niter
    use cms_def, only: noptset
    use const_def, only: small
    use prec_def, only: ikind
    
    implicit none
    integer :: i,k,kk,im,im1,icsh,nck
    integer :: iterxsh,idw,iup
    real(ikind) :: relaxcsh,relaxcsh1,cfuwc,fric_bed
    real(ikind) :: uhdelyyup,uhdelyydw,uhdelxxup,uhdelxxdw
    real(ikind) :: uxshim,vxshim,forcex,forcey
    
    !if(noptset<3)then !No waves
    !  do icsh=1,nCSstr
    !    im1=CS_str(icsh)%ncells/2        
    !    if(mod(idirface(im1,k),2)==0)then !East/West          
    !      CS_str(icsh)%ucsh(:)=0.0
    !    endif
    !  enddo
    !  return
    !endif
    relaxcsh=0.5
    relaxcsh1=1.0-relaxcsh
    iterxsh=max(5,30-niter**2)    
    do icsh=1,nCSstr
      im1=CS_str(icsh)%ncells/2 
      k  =CS_str(icsh)%faces(im1)
      !if(mod(idirface(im1,k),2)==0)then !East/West commented bdj 2021-02-25 
      if(mod(idirface(k,im1),2)==0)then !East/West 
        CS_str(icsh)%vcsh(:)=0.0
        do kk=1,iterxsh   
          im=1 !First node
          i  =CS_str(icsh)%cells(im)
          idw=CS_str(icsh)%cells(im+1)
          k  =CS_str(icsh)%faces(im)          
          nck=cell2cell(k,i)
          uhdelyydw=0.5*(vis(i)+vis(idw))*(h(i)+h(idw))/dy(i)/(dy(i)+dy(idw))*iwet(idw)
          if(noptset<3)then !No waves
            cfuwc=cfrict(i)
          else    
            cfuwc=fric_bed(h(i),cfrict(i),z0(i),CS_str(icsh)%ucsh(im),v(nck),&
                           us(i),vs(i),worb(i),worbrep(i),wper(i),Wang(i))   
          endif
          uxshim=(uhdelyydw*CS_str(icsh)%ucsh(im+1)+forcex(i))/(uhdelyydw+cfuwc+small)
          CS_str(icsh)%ucsh(im)=(relaxcsh1*CS_str(icsh)%ucsh(im)+relaxcsh*uxshim)*iwet(i)        
          do im=2,CS_str(icsh)%ncells-1
            i  =CS_str(icsh)%cells(im)
            iup=CS_str(icsh)%cells(im-1)
            idw=CS_str(icsh)%cells(im+1)
            k  =CS_str(icsh)%faces(im)
            nck=cell2cell(k,i)
            uhdelyyup=0.5*(vis(i)+vis(iup))*(h(i)+h(iup))/dy(i)/(dy(i)+dy(iup))*iwet(iup)
            uhdelyydw=0.5*(vis(i)+vis(idw))*(h(i)+h(idw))/dy(i)/(dy(i)+dy(idw))*iwet(idw)
            if(noptset<3)then !No waves
              cfuwc=cfrict(i)
            else
              cfuwc=fric_bed(h(i),cfrict(i),z0(i),CS_str(icsh)%ucsh(im),&
                             v(nck),us(i),vs(i),worb(i),worbrep(i),wper(i),Wang(i))
            endif
            uxshim=(uhdelyyup*CS_str(icsh)%ucsh(im-1)           &
                   +uhdelyydw*CS_str(icsh)%ucsh(im+1)+forcex(i)) &
                   /(uhdelyyup+uhdelyydw+cfuwc+small)
            CS_str(icsh)%ucsh(im)=(relaxcsh1*CS_str(icsh)%ucsh(im)+relaxcsh*uxshim)*iwet(i)   
          enddo
          im =CS_str(icsh)%ncells    !Last node
          i  =CS_str(icsh)%cells(im)
          iup=CS_str(icsh)%cells(im-1)
          k  =CS_str(icsh)%faces(im)
          nck=cell2cell(k,i)
          uhdelyyup=0.5*(vis(i)+vis(iup))*(h(i)+h(iup))/dy(i)/(dy(i)+dy(iup))*iwet(iup)
          if(noptset<3)then !No waves
            cfuwc=cfrict(i)
          else    
            cfuwc=fric_bed(h(i),cfrict(i),z0(i),CS_str(icsh)%ucsh(im),v(nck),&
                           us(i),vs(i),worb(i),worbrep(i),wper(i),Wang(i))  
          endif
          uxshim=(uhdelyyup*CS_str(icsh)%ucsh(im-1)+forcex(i))/(uhdelyyup+cfuwc+small)
          CS_str(icsh)%ucsh(im)=(relaxcsh1*CS_str(icsh)%ucsh(im)+relaxcsh*uxshim)*iwet(i)   
        enddo !kk-iteration  
      else  !Noth-South
        CS_str(icsh)%ucsh(:)=0.0
        do kk=1,iterxsh     
          im=1   !First node
          i  =CS_str(icsh)%cells(im)
          idw=CS_str(icsh)%cells(im+1)
          k  =CS_str(icsh)%faces(im)
          nck=cell2cell(k,i)
          uhdelxxdw=0.5*(vis(i)+vis(idw))*(h(i)+h(idw))/dx(i)/(dx(i)+dx(idw))*iwet(idw)   
          if(noptset<3)then !No waves
            cfuwc=cfrict(i)
          else
            cfuwc=fric_bed(h(i),cfrict(i),z0(i),u(nck),CS_str(icsh)%vcsh(im),&
                           us(i),vs(i),worb(i),worbrep(i),wper(i),Wang(i))
          endif
          vxshim=(uhdelxxdw*CS_str(icsh)%vcsh(im+1)+forcey(i))/(uhdelxxdw+cfuwc+small)
          CS_str(icsh)%vcsh(im)=(relaxcsh1*CS_str(icsh)%vcsh(im)+relaxcsh*vxshim)*iwet(i)
          do im=2,CS_str(icsh)%ncells-1
            i  =CS_str(icsh)%cells(im)
            iup=CS_str(icsh)%cells(im-1)
            idw=CS_str(icsh)%cells(im+1)
            k  =CS_str(icsh)%faces(im)
            nck=cell2cell(k,i)
            uhdelxxup=0.5*(vis(i)+vis(iup))*(h(i)+h(iup))/dx(i)/(dx(i)+dx(iup))*iwet(iup)
            uhdelxxdw=0.5*(vis(i)+vis(idw))*(h(i)+h(idw))/dx(i)/(dx(i)+dx(idw))*iwet(idw)
            if(noptset<3)then !No waves
              cfuwc=cfrict(i)
            else
              cfuwc=fric_bed(h(i),cfrict(i),z0(i),u(nck),CS_str(icsh)%vcsh(im),&
                             us(i),vs(i),worb(i),worbrep(i),wper(i),Wang(i))    
            endif
            vxshim=(uhdelxxup*CS_str(icsh)%vcsh(im-1)            &
                   +uhdelxxdw*CS_str(icsh)%vcsh(im+1)+forcey(i)) &
                   /(uhdelxxup+uhdelxxdw+cfuwc+small)
            CS_str(icsh)%vcsh(im)=(relaxcsh1*CS_str(icsh)%vcsh(im)+relaxcsh*vxshim)*iwet(i)
          enddo
          im =CS_str(icsh)%ncells    !Last node
          i  =CS_str(icsh)%cells(im)
          iup=CS_str(icsh)%cells(im-1)
          k  =CS_str(icsh)%faces(im)
          nck=cell2cell(k,i)
          uhdelxxup=0.5*(vis(i)+vis(iup))*(h(i)+h(iup))/dx(i)/(dx(i)+dx(iup))*iwet(iup)
          if(noptset<3)then !No waves
            cfuwc=cfrict(i)
          else    
            cfuwc=fric_bed(h(i),cfrict(i),z0(i),u(nck),CS_str(icsh)%vcsh(im),&
                           us(i),vs(i),worb(i),worbrep(i),wper(i),Wang(i))  
          endif
          vxshim=(uhdelxxup*CS_str(icsh)%vcsh(im-1)+forcey(i))/(uhdelxxup+cfuwc+small)
          CS_str(icsh)%vcsh(im)=(relaxcsh1*CS_str(icsh)%vcsh(im)+relaxcsh*vxshim)*iwet(i)
        enddo !kk-iteration  
      endif !East-West
    enddo !icsh
       
    return
    end subroutine xshore_uv
        
!****************************************************************          
    subroutine xshore_wse
! Solves the cross-shore momentum equation for water levels
! written by Weiming, Wu, NCCHE
! modified by Alex Sanchez, USACE-CHL
!****************************************************************  
    use geo_def, only: cell2cell, idirface
    use flow_def, only: p, gravinv
    use bnd_def
    use prec_def, only: ikind
    
    implicit none
    integer :: icsh,i,i1,k,im,im1,nck
    real(ikind) :: wsesetupx,wsesetupy
    do icsh=1,nCSstr
      im=1 !First node of each cross-shore bnd
      i =CS_str(icsh)%cells(im) 
      k =CS_str(icsh)%faces(im) 
      nck=cell2cell(k,i)
      !CS_str(icsh)%wsecsh(im) = p(nck)/grav !Old
      CS_str(icsh)%wsecsh(im) = p(i)*gravinv !New 09/10/13 
      im1=CS_str(icsh)%ncells/2
      k=CS_str(icsh)%faces(im1)
      !if(mod(idirface(im1,k),2)==1)then !North/South boundary commented bdj 2021-02-25 
      if(mod(idirface(k,im1),2)==1)then !North/South boundary
        do im=2,CS_str(icsh)%ncells   !Second node to end
          i1=CS_str(icsh)%cells(im-1)
          i =CS_str(icsh)%cells(im)          
          CS_str(icsh)%wsecsh(im)=CS_str(icsh)%wsecsh(im-1)+wsesetupx(i,i1)  
        enddo
      else                                !East/West boundary
        do im=2,CS_str(icsh)%ncells   !Second node to end
          i1=CS_str(icsh)%cells(im-1)
          i =CS_str(icsh)%cells(im)
          CS_str(icsh)%wsecsh(im)=CS_str(icsh)%wsecsh(im-1)+wsesetupy(i,i1)
        enddo
      endif
    enddo
    
    return
    end subroutine xshore_wse
    
!*********************************************************************
    subroutine bnd_wse_adjust(nbndcells,icells,kfaces,wsebnd,wseadj)
! Adjusts the water level at a boundary cell string to include 
! the wind and wave setup.
! written by Alex Sanchez, USACE-CHL    
!*********************************************************************
    use geo_def, only: rx,ry
    use prec_def, only: ikind
    
    implicit none
    !Input/Output    
    integer,intent(in)      :: nbndcells         !Boundary cells
    integer,intent(in)      :: icells(nbndcells) !Boundary cell id's
    integer,intent(in)      :: kfaces(nbndcells) !Boundary cell faces
    real(ikind),intent(in ) :: wsebnd(nbndcells) !Input water level at boundary
    real(ikind),intent(out) :: wseadj(nbndcells) !Output Water level at boundary
    !Internal
    integer :: i,i1,k,im,im0
    real(ikind) :: wsesetup(nbndcells) !Water level setup at boundary
    
    interface
      function wsesetupx(i,i1)
        use prec_def
        integer,intent(in) :: i,i1
        real(ikind) :: wsesetupx
      end function
    endinterface
    
    interface
      function wsesetupy(i,i1)
        use prec_def
        integer,intent(in) :: i,i1
        real(ikind) :: wsesetupy
      end function
    endinterface
    
    im0=nbndcells/2 !Center cell, starting or reference point
    wsesetup(im0) = 0.0 !Only need to initialize reference point
    !Upper side from reference cell
    do im=im0+1,nbndcells
      i1=icells(im-1)
      i =icells(im)
      k =kfaces(im)
      if(abs(ry(k,i))>=0.999)then !North/South boundary
        wsesetup(im)=wsesetup(im-1)+wsesetupx(i,i1)
      elseif(abs(rx(k,i))>=0.999)then !East/West boundary
        wsesetup(im)=wsesetup(im-1)+wsesetupy(i,i1)  
      else                            !Oblique boundary
        wsesetup(im)=wsesetup(im-1)+wsesetupx(i,i1)+wsesetupy(i,i1)  
      endif
    enddo
    !Lower side from reference cell
    do im=im0-1,1,-1
      i1=icells(im+1)
      i =icells(im)
      k =kfaces(im)
      if(abs(ry(k,i))>=0.999)then !North/South boundary
        wsesetup(im)=wsesetup(im+1)+wsesetupx(i,i1) 
      elseif(abs(rx(k,i))>=0.999)then !East/West boundary
        wsesetup(im)=wsesetup(im+1)+wsesetupy(i,i1)
      else                            !Oblique boundary
        wsesetup(im)=wsesetup(im+1)+wsesetupx(i,i1)+wsesetupy(i,i1)
      endif
    enddo
    wseadj = wsebnd + wsesetup
    
    return
    end subroutine bnd_wse_adjust
    
!*********************************************************
    subroutine bndpath2id(apath,id)
!*********************************************************    
    implicit none
    integer :: i,j,nn,ierr
    integer,intent(out) :: id    
    character(len=*),intent(in) :: apath
    
    nn = len_trim(apath)
    id = 1
    ierr = -1  !added 5/21/2018 meb
    if(nn==0) return    
    
    !Test to see if apath is only a number
    select case(nn)
    case (1)
      read(apath,'(I1)',iostat=ierr) id
    case (2)
      read(apath,'(I2)',iostat=ierr) id
    case (3)
      read(apath,'(I3)',iostat=ierr) id
    end select 

    if(ierr==0) then 
      return
    endif
        
    do i=nn,1,-1
      select case(apath(i:i))
      case('#')        
        read(apath(i+1:nn),*) id
        return
          
      case('(')
        do j=i+1,nn
          if(apath(j+1:j+1)==')')then
            read(apath(i+1:j),*) id
            return
          endif
        enddo
          
      case('/')
        id = 1
        return

      end select  
    enddo
    
    return
    end subroutine bndpath2id

!**************************************************************************
    subroutine fluxbnd(nbndcells,icells,kfaces)
! Flux boundary condition
! written by Alex Sanchez, USACE-CHL
!*************************************************************************
    use geo_def, only: ds,cell2cell,llec2llec,fnx,fny
    use flow_def, only: flux,iwet,u,v,h,hk
    
    implicit none
    !Input/Output
    integer,intent(in) :: nbndcells
    integer,intent(in) :: icells(nbndcells),kfaces(nbndcells)
    !Internal Variables
    integer :: i,j,k,nck,jcn
    
    do j=1,nbndcells
      i=icells(j)
      k=kfaces(j)
      nck=cell2cell(k,i) !Forward connectivity 
      flux(k,i)=iwet(i)*ds(k,i)*hk(k,i)*(fnx(k,i)*u(nck)+fny(k,i)*v(nck))
      jcn=llec2llec(k,i) !Backwards connectivity
      flux(jcn,nck)=-flux(k,i)
    enddo
    
    return
    end subroutine fluxbnd    
    
!**************************************************************************
    subroutine bnd_pp(nbndcells,icells,kfaces,wsebnd)
! Pressure correction boundary condition
!*************************************************************************
    use size_def, only: ncellsimple,ncelljoint,ncellpoly,ncells
    use geo_def, only: cell2cell,idirface,ncface,kkface,fnx,fny
    use flow_def, only: flux,iwet,p,pp,su,sp,acoef,grav
    use prec_def, only: ikind
    
    implicit none
    !Input/Output
    integer,intent(in) :: nbndcells
    integer,intent(in) :: icells(nbndcells),kfaces(nbndcells)
    real(ikind),intent(in) :: wsebnd(nbndcells)
    !Internal Variables
    integer :: i,j,k,nck !,kk,kkdf    
    !real(ikind) :: acoefik
    !real(ikind) :: fac,sumfac
    
    do j=1,nbndcells
      i=icells(j)
      k=kfaces(j)
      nck=cell2cell(k,i)
      pp(nck)=wsebnd(j)*grav-p(nck)   !Set given wse at the dummy node
      su(i)=su(i)+acoef(k,i)*pp(nck)  !using a_kk because a_k is not well defined  
      sp(i)=sp(i)-acoef(k,i)
      acoef(k,i)=0.0
    enddo !j
    
    !if(ncellsimple==ncells)then !Nonuniform Cartesian grid
    !  do j=1,nbndcells
    !    i=icells(j)
    !    k=kfaces(j)
    !    nck=cell2cell(k,i)
    !    pp(nck)=wsebnd(j)*grav-p(nck)   !Set given wse at the dummy node
    !    acoefik=acoef(kkface(k),i)
    !    su(i)=su(i)+acoefik*pp(nck)  !using a_kk because a_k is not well defined  
    !    sp(i)=sp(i)-acoefik
    !    acoef(k,i)=0.0
    !  enddo !j
    !elseif(ncelljoint>0)then !Telescoping grid
    !  do j=1,nbndcells
    !    i=icells(j)
    !    k=kfaces(j)
    !    nck=cell2cell(k,i)
    !    pp(nck)=wsebnd(j)*grav-p(nck)   !Set given wse at the dummy node
    !    if(ncface(i)==4)then
    !      acoefik=acoef(kkface(k),i)
    !    else
    !      acoefik=0.0; kkdf=kkface(idirface(k,i))
    !      do kk=1,ncface(i)    !Assume boundary face does not split
    !        if(idirface(kk,i)==kkdf) acoefik=acoefik+acoef(kk,i)
    !      enddo                   
    !    endif                        
    !    su(i)=su(i)+acoefik*pp(nck)  !using a_kk because a_k is not well defined  
    !    sp(i)=sp(i)-acoefik
    !    acoef(k,i)=0.0
    !  enddo !j  
    !elseif(ncellpoly>0)then !Unstructured Polyhedral grid
    !  do j=1,nbndcells
    !    i=icells(j)
    !    k=kfaces(j)
    !    nck=cell2cell(k,i)
    !    pp(nck)=wsebnd(j)*grav-p(nck)   !Set given wse at the dummy node
    !    !!-- set coefficient as average of opposite face(s) -----
    !    !acoefik=0.0; sumfac=0.0
    !    !do kk=1,ncface(i)    !Assume boundary face does not split
    !    !  if(kk==k) cycle
    !    !  fac=max(-fnx(k,i)*fnx(i,kk),0.0)+max(-fny(k,i)*fny(i,kk),0.0) !Determine if face is opposite facing
    !    !  acoefik=acoefik+acoef(kk,i)*fac
    !    !  sumfac=sumfac+fac
    !    !enddo 
    !    !acoefik=acoefik/max(sumfac,1.0e-15)
    !    !--- Use acoef(k,i) since acoef(kk,i) is not well defined ---
    !    acoefik=acoef(k,i)
    !    su(i)=su(i)+acoefik*pp(nck)  !using a_kk because a_k is not well defined  
    !    sp(i)=sp(i)-acoefik
    !    acoef(k,i)=0.0
    !  enddo !j
    !endif 
    
    return
    end subroutine bnd_pp
    
!**************************************************************************
    subroutine bnd_p(nbndcells,icells,kfaces,wsebnd)
! Pressure boundary condition
!*************************************************************************
    use size_def, only: ncellsimple,ncelljoint,ncellpoly,ncells
    use geo_def, only: cell2cell,idirface,ncface,kkface,fnx,fny,zb
    use flow_def, only: flux,iwet,p,h,su,sp,acoef,grav
    use prec_def, only: ikind
    
    implicit none
    !Input/Output
    integer,intent(in) :: nbndcells
    integer,intent(in) :: icells(nbndcells),kfaces(nbndcells)
    real(ikind),intent(in) :: wsebnd(nbndcells)
    !Internal Variables
    integer :: i,j,k,nck !,kk,kkdf    
    !real(ikind) :: acoefik
    !real(ikind) :: fac,sumfac
    
    do j=1,nbndcells
      i=icells(j)
      k=kfaces(j)
      nck=cell2cell(k,i)
      p(nck)=wsebnd(j)*grav
      h(nck)=wsebnd(j)-zb(nck) !eta(nck)=h(nck)+zb(nck)
      su(i)=su(i)+acoef(k,i)*p(nck) !using a_kk because a_k is not well defined  
      sp(i)=sp(i)-acoef(k,i)
      acoef(k,i)=0.0
    enddo !j
    
    return
    end subroutine bnd_p    
    
!**************************************************************************
    subroutine bnd_wse(nbndcells,icells,kfaces,wsebnd)
! Pressure boundary condition
!*************************************************************************
    use geo_def, only: cell2cell,zb,ds,fnx,fny
    use flow_def, only: p,h,eta,u,v,hk,flux,grav
    use prec_def, only: ikind
    
    implicit none
    !Input/Output
    integer,intent(in) :: nbndcells
    integer,intent(in) :: icells(nbndcells),kfaces(nbndcells)
    real(ikind),intent(in) :: wsebnd(nbndcells)
    !Internal Variables
    integer :: i,j,k,nck !,kk,kkdf
    
    do j=1,nbndcells
      i=icells(j)
      k=kfaces(j)
      nck=cell2cell(k,i)
      eta(nck)=wsebnd(j)  !Set given wse at the dummy node
      p(nck)=eta(nck)*grav 
      h(nck)=eta(nck)-zb(nck)
      !u(nck)=u(i)*h(i)/h(nck)
      !v(nck)=v(i)*h(i)/h(nck)
      u(nck)=fnx(k,i)*flux(k,i)/(ds(k,i)*hk(k,i))
      v(nck)=fny(k,i)*flux(k,i)/(ds(k,i)*hk(k,i))
    enddo !j
    
    return
    end subroutine bnd_wse
    
!********************************************************    
    subroutine bndvelxtrap(nbndcells,icells,kfaces)
! Extrapolates the boundary velocity using
! a simple flux formulation
!********************************************************    
    use geo_def, only: cell2cell
    use flow_def, only: h,u,v,iwet
    use bnd_def, only: bnd_str
    
    implicit none
    !Input
    integer,intent(in) :: nbndcells
    integer,intent(in),dimension(nbndcells) :: icells,kfaces
    !Internal Variables
    integer :: i,j,k,nck
    
    do j=1,nbndcells
      i=icells(j)
      k=kfaces(j)
      nck=cell2cell(k,i)
      u(nck)=u(i)*h(i)/h(nck)*iwet(i)
      v(nck)=v(i)*h(i)/h(nck)*iwet(i)
    enddo
    
    return
    end subroutine bndvelxtrap

!********************************************************    
    subroutine bndpresdepxtrap(nbndcells,icells,kfaces)
! Extrapolates the boundary water depth and 
! pressure
!********************************************************    
    use geo_def, only: cell2cell,zb
    use flow_def, only: h,p,iwet,hmin,gravinv
    use bnd_def, only: bnd_str
    
    implicit none
    !Input
    integer,intent(in) :: nbndcells
    integer,intent(in),dimension(nbndcells) :: icells,kfaces
    !Internal Variables
    integer :: i,j,k,nck
    
    do j=1,nbndcells
      i=icells(j)
      k=kfaces(j)
      nck=cell2cell(k,i)
      p(nck)=p(i)
      !p(nck)=p(i)+2.0*(rx(k,i)*dpx(i)+ry(k,i)*dpy(i)) !Linear extrapolation
      !p(cell2cell(k,i))=2.0*p(i)-p(cell2cell(llec2llec(k,i),i)) !Linear extrapolation
      h(nck)=max(hmin,p(nck)*gravinv-zb(nck))
    enddo
    
    return
    end subroutine bndpresdepxtrap
    
