!************************************************************
    subroutine bnd_block
! Reads a boundary block:
! BOUNDARY_BEGIN
!   <CARDS>
! BOUNDARY_END
!
! written by Alex Sanchez, USACE-CHL
!************************************************************
    use bnd_def, only: nparsim,parsim,q_str,nqstr,th_str,nthstr,h_str,nhstr,mh_str,nmhstr,mhv_str,nmhvstr,nh_str,nnhstr,nhv_str,nnhvstr
    use bnd_def, only: nth_str,nnthstr,nthv_str,nnthvstr,ioffsetmode !1-Constant offset, 2-Offset curve (hli,01/18/17)
    use sal_def, only: nsalstr,sal_str
    use geo_def, only: azimuth_fl,projection
    use geo_lib, only: proj_default
    use const_def, only: deg2rad
    use comvarbl, only: mpfile,flowpath,tjulday0,iyr
    use time_lib, only: calendar2julian
    use tide_lib, only: tidal_data
    !use unitconv_lib
    use diag_lib, only: diag_print_error, diag_print_warning
    use prec_def, only: ikind
    use cms_def,  only: aValue
    implicit none    

    real(ikind), parameter :: undef = -999.0
    integer, parameter     :: iundef = -999
    integer            :: i,k,kk,ierr
    logical            :: foundcard,foundfile
    character(len=20)  :: qunits
    character(len=37)  :: cardname
    
    !Boundary String
    character(len=200) :: bidfile,bidpath    !Boundary ID file and path
    character(len=100) :: bndname
    integer :: idnum       !Boundary id number
    integer :: istrtype    !Input string type, 1-cellstring, 2-nodestring
    integer :: ibndtype    !Boundary type,1=Q,2=T,3=H,4=MH,5=MHV,6=CS,7=NH,8=NHV,9=NTH,10-NTHV
    
    !Flux Boundary
    integer     :: ifluxmode   !1-Constant flux, 2-Total Flux curve
!    integer     :: ioffsetmode   !1-Constant offset, 2-Offset curve (hli,01/18/17)
    logical     :: fluxblockread
    real(ikind) :: cmvel       !Coefficient [-]
    real(ikind) :: angle_flux  !Angle of flow clockwise from north [rad]
    real(ikind) :: qfluxconst
    integer     :: ifluxunits  !Input flux units: 0-m^3/s/cell,1-m^3/s,2-ft^3/s

    character(len=200) :: fluxfile,fluxpath  !Flux data file and path
    character(len=200) :: offsetfile,offsetpath  !offset data file and path (hli,01/18/17)
    integer    :: ntiq 
    
    !Time-series Smoothing
    integer :: nsw         !Temporal smoothing width
    integer :: nsi         !Temporal smoothing iterations
    integer :: nswwse      !Temporal smoothing width
    integer :: nsiwse      !Temporal smoothing iterations
    integer :: nswvel      !Temporal smoothing width
    integer :: nsivel      !Temporal smoothing iterations
    
    !Tidal Boundary Condition
    integer              :: ntc         !Tidal constituents used
    logical              :: istidal     !true for tidal, false for harmonic
    real(ikind)          :: angle_wave  !Incident angle of tidal wave
    real(ikind), pointer :: amp(:)      !Amplitude [m] (constituent) 
    real(ikind), pointer :: speed(:)    !Speed [rad/hrs] (constituent)
    real(ikind), pointer :: phase(:)    !Phase [rad] (constituent)
    real(ikind), pointer :: f(:)        !Nodal factor [-] (constituent)
    real(ikind), pointer :: vu(:)       !Equilibrium argument [rad] (constituent)
    character(len=10), pointer :: name(:) !Tidal Consitituent names (constituent)
    character(len=100) :: station
    
    !WSE Block
    logical :: wseblockread
    integer :: minterp  !Method for interpolation, 1-Piecewise polynomial, 2-cubic spline
    integer :: ntiwse   !Interpolation order
    character(len=200) :: wsefile,wsepath !Water level data file and path
    real(ikind) :: wseoffset !wse offset
    real(ikind) :: dwsex     !Regional steady water level gradinet
    real(ikind) :: dwsey     !Regional steady water level gradinet
    logical :: wseout        !Water level output file
    logical :: wseadjust      !Turns on or off the wse adjustment/correction due to wind and waves 
    
    !Vel Block
    logical :: velblockread
    integer :: ntivel
    real(ikind) :: wseconst
    character(len=200) :: velfile,velpath  !Velocity data file and path
    logical :: velout           !Velocity output file
    
    !Nested Boundary Conditions
    logical :: nestblockread
    integer :: idpar
    real(ikind) :: tjuldaypar    !Parent simulation reference (starting) time in Julian days
    real(ikind) :: timestarthr   !Time relative to the CMS-Flow starting time (used for ADCIRC)
    character(len=200) :: wsefilepar,wsepathpar !Water level data file and path
    character(len=200) :: velfilepar,velpathpar !Velocity data file and path
    character(len=200) :: ctlfilepar            !Parent control file and path
    character(len=200) :: grdfilepar            !Parent grid file and path
    type(projection) :: projpar
    
    !Tidal Database Boundary Conditions
    logical :: tdbblockread
    integer :: ntcin  !Tidal constituents used     
    character(len=10), pointer :: namein(:) !Input Tidal Consitituent names (constituent)
    character(len=10)  :: tdbname !Tidal Database Name, EC2001, ENPAC2003, LEPROVOST, 
    character(len=200) :: tdbpath !Tidal Database file and path
    type(projection)   :: projtdb !Parent grid projection
    
    !Spatial Smoothing
    integer :: nssi  !Smoothing iterations (along string)
    integer :: nssw  !Smoothing window width (along string)
    integer :: nssiwse  !Smoothing iterations (along string)
    integer :: nsswwse  !Smoothing window width (along string)
    integer :: nssivel  !Smoothing iterations (along string)
    integer :: nsswvel  !Smoothing window width (along string)
    
    !Salinity
    integer :: isaltype
    character(len=200) :: salfile,salpath
    real(ikind) :: salbnd
    
    !Sediment
    integer :: isedtype
    character(len=200) :: sedfile,sedpath
    real(ikind) :: sedbnd
    
    interface
      subroutine tidal_block(ibndtype,ntc,name,amp,phase,speed,f,vu,angle_wave,ioffsetmode,offsetfile,offsetpath,wseoffset,nti)   !hli(10/04/17)
        use prec_def
        implicit none
        integer,intent(inout) :: ibndtype
        integer              :: ntc          !Tidal constituents used
        real(ikind)          :: angle_wave   !Incident angle of tidal wave
        real(ikind), pointer :: amp(:)       !Amplitude [m] (constituent) 
        real(ikind), pointer :: speed(:)     !Speed [rad/hrs] (constituent)
        real(ikind), pointer :: phase(:)     !Phase [rad] (constituent)
        real(ikind), pointer :: f(:)         !Nodal factor [-] (constituent)
        real(ikind), pointer :: vu(:)        !Equilibrium argument [rad] (constituent)
        character(len=10), pointer :: name(:) !Tidal Consitituent names (constituent)    
        character(len=*),intent(inout) :: offsetfile,offsetpath !(hli,10/04/17)
        integer,         intent(out)   :: ioffsetmode           !(hli,10/04/17)
        integer,         intent(out)   :: nti                   !(hli,10/04/17)
        real(ikind),     intent(inout) :: wseoffset             !(hli,10/04/17)
      end subroutine
    endinterface
    
    interface
      subroutine tsta_block(ibndtype,station,ntc,name,amp,phase,speed,f,vu,angle_wave)
        use prec_def
        implicit none
        integer,intent(inout) :: ibndtype
        character(len=*)     :: station
        integer              :: ntc          !Tidal constituents used
        real(ikind)          :: angle_wave   !Incident angle of tidal wave
        real(ikind), pointer :: amp(:)       !Amplitude [m] (constituent) 
        real(ikind), pointer :: speed(:)     !Speed [rad/hrs] (constituent)
        real(ikind), pointer :: phase(:)     !Phase [rad] (constituent)
        real(ikind), pointer :: f(:)         !Nodal factor [-] (constituent)
        real(ikind), pointer :: vu(:)        !Equilibrium argument [rad] (constituent)
        character(len=10), pointer :: name(:) !Tidal Consitituent names (constituent)    
      end subroutine
    endinterface
    
    interface
      subroutine harmonic_block(ibndtype,ntc,amp,phase,speed,angle_wave)
        use prec_def
        implicit none
        integer,intent(inout) :: ibndtype
        integer              :: ntc         !Tidal constituents used
        real(ikind)          :: angle_wave  !Incident angle of tidal wave
        real(ikind), pointer :: amp(:)      !Amplitude [m] (constituent) 
        real(ikind), pointer :: phase(:)    !Phase [rad] (constituent)
        real(ikind), pointer :: speed(:)    !Speed [rad/hrs] (constituent)
      end subroutine
    endinterface
    
    interface
      subroutine tdb_block(ntcin,namein,tdbname,tdbpath,projtdb,nssi,nssw)
        use geo_def, only: projection
        use prec_def
        implicit none
        integer,intent(out) :: ntcin !Tidal constituents used     
        character(len=*),intent(inout),pointer :: namein(:) !Input Tidal Consitituent names (constituent)
        character(len=*),intent(inout) :: tdbname !Tidal Database Name, EC2001, ENPAC2003, LEPROVOST, 
        character(len=*),intent(inout) :: tdbpath !Tidal Database file and path
        type(projection),intent(inout) :: projtdb !Parent grid projection
        integer,intent(inout) :: nssi      !Smoothing iterations (along string)
        integer,intent(inout) :: nssw     !Smoothing window width (along string)
      end subroutine 
    endinterface
    
    !--- Initialize -----
    bidfile = ''
    bidpath = ''
    bndname = ''
    idnum = iundef
    istrtype = iundef
    ibndtype = iundef
    
    !Flux
    ifluxmode = 2   !1-Constant flux, 2-Total-flux curve
    fluxblockread = .false.
    cmvel= 0.667    !Coefficient [-]
    angle_flux = undef  !Angle of flow clockwise from north [rad]
    qfluxconst = 0.0   
    ifluxunits = 0  !Input flux units: 0-m^3/s/cell,1-m^3/s,2-ft^3/s
    fluxfile = ''
    fluxpath = ''
    qunits = 'm^3/s/cell'
    ntiq = 1 !Linear
    
    !Temporal Smoothing (same for all boundaries)
    nsi = 0
    nsw = 0
    nswwse = 0  !Smoothing width
    nsiwse = 0  !Smoothing iterations for wse
    nswvel = 0  !Smoothing width for wse
    nsivel = 0  !Smoothing iterations for velocities
    
    !Spatial smoothing 
    nssi = 0     !Smoothing iterations (along string)
    nssw = 0     !Smoothing window width (along string)
    nssiwse = 0  !Smoothing iterations (along string)
    nsswwse = 0  !Smoothing window width (along string)
    nssivel = 0  !Smoothing iterations (along string)
    nsswvel = 0  !Smoothing window width (along string)
    
    !Tidal Boundary Condition
    ntc = 0            !Tidal constituents used
    istidal = .true.   !true for tidal, false for harmonic
    angle_wave = undef !Incident angle of tidal wave
    station = ''
    
    !WSE Boundary Conditions
    wseblockread = .false.
    wseconst  = 0.0
    wseoffset = 0.0  !wse offset
    dwsex     = 0.0  !Regional steady water level gradinet
    dwsey     = 0.0  !Regional steady water level gradinet
    minterp = 1  !Piecewise polynomial interpolation
    ntiwse = 2   !second order
    wsefile = '' !Water level data file
    wsepath = '' !Water level data path
    wseout = .false.    
    wseadjust = .true.
    
    !Vel Boundary Conditions
    velblockread = .false.
    ntivel = 1
    velfile = '' !Velocity data file
    velpath = '' !Velocity data path
    velout = .false.
    
    !Nested Boundary Conditions
    nestblockread = .false.
    idpar = 0
    timestarthr = undef
    tjuldaypar = undef
    wsefilepar = ''   !Water level data file
    wsepathpar = ''   !Water level data path
    velfilepar = ''   !Velocity data file
    velpathpar = ''   !Velocity data path
    ctlfilepar = ''   !Parent control file and path
    grdfilepar = ''   !Parent grid file and path
    call proj_default(projpar) !Parent grid projection
    
    !Tidal Database Boundary Conditions
    tdbblockread = .false.
    ntcin  = 0  !Tidal constituents used
    tdbname = ''  !Tidal Database Name, EC2001, ENPAC2003, LEPROVOST, 
    tdbpath = ''  !Tidal Database file and path
    nssi = iundef  !Smoothing iterations (along string)
    nssw = iundef !Smoothing window width (along string)
    call proj_default(projtdb) !Parent grid projection
    
    !Salinity
    isaltype = iundef
    salfile = ''
    salpath = ''
    salbnd = undef
    
    !Sediment
    isedtype = iundef
    sedfile = ''
    sedpath = ''
    sedbnd = undef
    
    do !kk=1,60
      foundcard = .true.
      read(77,*,iostat=ierr) cardname
      if(ierr/=0) exit
      if(cardname(1:1)=='!' .or. cardname(1:1)=='#' .or. cardname(1:1)=="*") cycle  
      select case(cardname)
      case('BOUNDARY_END','END')
        exit
          
      case('NAME')
        backspace(77)
        read(77,*) cardname,bndname
        
      case('CELLSTRING')
        call card_bid(flowpath,bidfile,bidpath,idnum)
        istrtype = 1
          
      case('NODESTRING')  
        call card_bid(flowpath,bidfile,bidpath,idnum)
        istrtype = 2 
      
      !case('COORDINATES') !Not implimented yet
      !  backspace(77)
      !  read(77,*) cardname,bxyfile
          
      case('FLUX_BEGIN','FLUX_BOUNDARY_BEGIN')
        call flux_block(ibndtype,ifluxmode,fluxfile,fluxpath,qfluxconst,ifluxunits,angle_flux,cmvel,ntiq,nsi,nsw)
        ibndtype = 1
        fluxblockread = .true.
        
      case('WSE_BEGIN','WSE_BOUNDARY_BEGIN','WSE_FORCING_BEGIN')
        call wse_block(ibndtype,istidal,wsefile,wsepath,wseconst,ioffsetmode,offsetfile,offsetpath,wseoffset, &  !(hli)
           dwsex,dwsey,minterp,ntiwse,nsiwse,nswwse,nssiwse,nsswwse,wseout,wseadjust)
        wseblockread = .true.
        
      case('VEL_BEGIN')
        call vel_block(ibndtype,velfile,velpath,ntivel,nsivel,nswvel,nssivel,nsswvel,velout) !Note: same interpolation order for wse and vel
        velblockread = .true.
        
      case('TIDAL_CONSTITUENTS_BEGIN','TIDAL_BEGIN')     
        call tidal_block(ibndtype,ntc,name,amp,phase,speed,f,vu,angle_wave,ioffsetmode,offsetfile,offsetpath,wseoffset,ntiwse)       !(hli 10/04/17)
        istidal = .true.
        
      case('HARMONIC_CONSTITUENTS_BEGIN','HARMONIC_BEGIN')
        call harmonic_block(ibndtype,ntc,amp,phase,speed,angle_wave)
        istidal = .false.
      
      case('TIDAL_STATION_BEGIN','STATION_BEGIN')
        call tsta_block(ibndtype,station,ntc,name,amp,phase,speed,f,vu,angle_wave)
        istidal = .true.
        
      case('PARENT_BEGIN','PARENT_SIMULATION_BEGIN')
        if(ibndtype==0) ibndtype=7  
        if(ibndtype==7 .or. ibndtype==8)then
          call parent_block(ctlfilepar,grdfilepar,projpar,wsefilepar,wsepathpar,velfilepar,velpathpar,tjuldaypar,timestarthr,ntiwse)
        else
          call diag_print_error('Found conflicting boundary specifications','  Parent block specification conflicts with the boundary type')
        endif
        nestblockread = .true.
        
      case('TIDAL_DATABASE_BEGIN','DATABASE_BEGIN')  
        call tdb_block(ntcin,namein,tdbname,tdbpath,projtdb,nssi,nssw)
        tdbblockread = .true.
       
      case('SALINITY_BLOCK','SALT_BLOCK','SAL_BLOCK')
        !call sal_block
        
      !case('SEDIMENT_TOTAL_FLUX_CURVE')
      !  call card_dataset(77,mpfile,flowpath,sedfile,sedpath)
      !  isedtype = 1 !Total sediemnt transport rate [kg/s]
      !
      !case('SEDIMENT_FRAC_FLUX_CURVE')
      !  call card_dataset(77,mpfile,flowpath,sedfile,sedpath)
      !  isedtype = 2 !Total sediemnt transport rate [kg/s]
      !  
      !case('SEDIMENT_TOTAL_FLUX_CONSTANT')
      !  backspace(77)
      !  read(77,*) cardname,sedbnd
      !  isedtype = 1
      !    
      !case('SEDIMENT_FRAC_FLUX_CONSTANT')
      !  backspace(77)
      !  read(77,*) cardname,sedbnd
      !  isedtype = 2
      
      case default
        foundcard = .false.
        call diag_print_warning('Unrecognized Card: '//trim(cardname))
      end select
    enddo
    
    if(len_trim(bidfile)==0)then
      call diag_print_error('Missing Boundary ID file')
    endif
    
    if(fluxblockread)then
      if(wseblockread .or. velblockread)then
        call diag_print_error('Invalid boundary specification',&
          '  Cannot specify both flux and wse/vel conditions',&
          '  at the same boundary')
      endif
    endif
    
    if(nestblockread)then
      if(.not.wseblockread)then
        call diag_print_warning('WSE block not specified for nested BC')
      endif
      if(velblockread)then
        ibndtype = 8 !NHV Nested wse and vel
      else
        ibndtype = 7 !NH Nested wse
      endif
      if(nParSim==0)then
        call parsim_init
        idpar = nParSim
      else
        call check_parent_files(grdfilepar,ctlfilepar) !wsefilepar,wsepathpar
        do i=1,nParSim
          if(grdfilepar==ParSim(i)%grdfilepar)then
            idpar = i
            exit
          endif
        enddo
        if(idpar==0)then
          call parsim_init
          idpar = nParSim
        endif
      endif
      ParSim(idpar)%ctlfilepar  = ctlfilepar
      ParSim(idpar)%grdfilepar  = grdfilepar
      ParSim(idpar)%projpar     = projpar
      ParSim(idpar)%tjuldaypar  = tjuldaypar
      ParSim(idpar)%timestarthr = timestarthr
      ParSim(idpar)%wsefilepar  = wsefilepar
      ParSim(idpar)%wsepathpar  = wsepathpar
      ParSim(idpar)%ntiwsepar   = max(ParSim(idpar)%ntiwsepar,ntiwse)
      if(ibndtype==8)then
        ParSim(idpar)%velpar     = .true.  
        ParSim(idpar)%ntivelpar  = max(ParSim(idpar)%ntivelpar,ntivel)
        ParSim(idpar)%velfilepar = velfilepar
        ParSim(idpar)%velpathpar = velpathpar
      endif
    endif
    
    if(tdbblockread)then
      if(.not.wseblockread)then
        call diag_print_warning('WARNING: WSE block not specified for tidal database BC')
      endif
      if(velblockread)then
        ibndtype = 10 !NTHV Nested wse and vel
      else
        ibndtype = 9 !NTH Nested wse
      endif  
    endif
    
    !Flow boundary type
    select case(ibndtype)
    case(1) !Q - Flux
      call flux_alloc 
      Q_str(nQstr)%bidfile = bidfile
      Q_str(nQstr)%bidpath = bidpath
      Q_str(nQstr)%idnum = idnum
      Q_str(nQstr)%istrtype = istrtype
      Q_str(nQstr)%ifluxmode = ifluxmode
      if(Q_str(nQstr)%ifluxmode==1)then !Constant
        Q_str(nQstr)%qfluxconst = qfluxconst  
      else !Time-series or stage-flow curve
        if(len_trim(fluxfile)==0)then
          fluxfile = bidfile
          fluxpath = bidpath
        endif 
        Q_str(nQstr)%fluxfile = fluxfile
        Q_str(nQstr)%fluxpath = fluxpath
      endif
      Q_str(nQstr)%angle = angle_flux
      Q_str(nQstr)%cmvel = cmvel
      Q_str(nQstr)%ifluxunits = ifluxunits
      Q_str(nQstr)%nti = ntiq
      Q_str(nQstr)%nsi = nsi
      Q_str(nQstr)%nsw = nsw
      
    case(2) !TH - Tidal
      call tidal_alloc  
      TH_str(nTHstr)%bidfile = bidfile
      TH_str(nTHstr)%bidpath = bidpath
      TH_str(nTHstr)%idnum = idnum
      TH_str(nTHstr)%istrtype = istrtype
      TH_str(nTHstr)%angle = angle_wave
      TH_str(nTHstr)%wseoffset = wseoffset      
      TH_str(nTHstr)%ntc = ntc
      TH_str(nTHstr)%wseadjust = wseadjust
      TH_str(nTHstr)%station = station
      TH_str(nTHstr)%offsetfile = offsetfile !hli
      TH_str(nTHstr)%offsetpath = offsetpath !hli
      TH_str(nTHstr)%nti = ntiwse !hli
      if(ntc==0) return
      allocate(TH_str(nTHstr)%amp(ntc),TH_str(nTHstr)%speed(ntc))
      allocate(TH_str(nTHstr)%phase(ntc),TH_str(nTHstr)%name(ntc))
      allocate(TH_str(nTHstr)%f(ntc),TH_str(nTHstr)%vu(ntc))
      TH_str(nTHstr)%istidal = istidal 
      do k=1,ntc
        TH_str(nTHstr)%amp(k)   = amp(k)
        TH_str(nTHstr)%phase(k) = phase(k)
        TH_str(nTHstr)%speed(k) = speed(k)
      enddo
      deallocate(amp,phase,speed)
      if(istidal)then
        do k=1,ntc
          TH_str(nTHstr)%f(k)     = f(k)
          TH_str(nTHstr)%vu(k)    = vu(k)
          TH_str(nTHstr)%name(k)  = name(k)
        enddo
        deallocate(f,vu,name)
      endif
      
    case(3) !H - Single WSE
      call singlewse_alloc
      H_str(nHstr)%bidfile = bidfile
      H_str(nHstr)%bidpath = bidpath
      H_str(nHstr)%idnum = idnum
      H_str(nHstr)%istrtype = istrtype
      H_str(nHstr)%wseconst = wseconst
      H_str(nHstr)%wsefile = wsefile
      H_str(nHstr)%wsepath = wsepath
      H_str(nHstr)%dwsex = dwsex
      H_str(nHstr)%dwsey = dwsey
      H_str(nHstr)%wseoffset = wseoffset
      H_str(nHstr)%minterp = minterp
      H_str(nHstr)%nti = ntiwse
      H_str(nHstr)%nsi = nsiwse
      H_str(nHstr)%nsw = nswwse
      H_str(nHstr)%wseadjust = wseadjust
      H_str(nHstr)%offsetfile = offsetfile !hli
      H_str(nHstr)%offsetpath = offsetpath !hli
      
      
    case(4) !MH - Multiple WSE
      call multiwse_alloc
      MH_str(nMHstr)%bidfile = bidfile
      MH_str(nMHstr)%bidpath = bidpath
      MH_str(nMHstr)%idnum = idnum
      MH_str(nMHstr)%istrtype = istrtype
      if(len_trim(wsefile)==0)then
        wsefile = bidfile  
        wsepath = bidpath
      endif
      MH_str(nMHstr)%wsefile = wsefile
      MH_str(nMHstr)%wsepath = wsepath
      MH_str(nMHstr)%nti = ntiwse
      MH_str(nMHstr)%nsw = nswwse     !Temporal smoothing width
      MH_str(nMHstr)%nsi = nsiwse     !Temporal smoothing iterations
      MH_str(nMHstr)%nssi = nssiwse  !Spatial smoothing iterations
      MH_str(nMHstr)%nssw = nsswwse !Spatial smoothing window width
      MH_str(nMHstr)%wseoffset = wseoffset
    
    case(5) !MHV - Multiple WSE and Vel
      call multiwsevel_alloc
      MHV_str(nMHVstr)%bidfile = bidfile
      MHV_str(nMHVstr)%bidpath = bidpath
      MHV_str(nMHVstr)%idnum = idnum
      MHV_str(nMHVstr)%istrtype = istrtype
      if(len_trim(wsefile)==0)then
        wsefile = bidfile  
        wsepath = bidpath
      endif
      MHV_str(nMHVstr)%wsefile = wsefile
      MHV_str(nMHVstr)%wsepath = wsepath
      MHV_str(nMHVstr)%ntiwse = ntiwse
      MHV_str(nMHVstr)%nswwse = nswwse          !Temporal smoothing width
      MHV_str(nMHVstr)%nsiwse = nsiwse          !Temporal smoothing iterations
      MHV_str(nMHVstr)%ntiwse = ntiwse
      MHV_str(nMHVstr)%nssiwse = nssiwse  !Spatial smoothing iterations
      MHV_str(nMHVstr)%nsswwse = nsswwse !Spatial smoothing window width
      MHV_str(nMHVstr)%wseoffset = wseoffset
      if(len_trim(velfile)==0)then
        velfile = bidfile  
        velpath = bidpath
      endif
      MHV_str(nMHVstr)%velfile = velfile
      MHV_str(nMHVstr)%velpath = velpath
      MHV_str(nMHVstr)%nswvel = nswvel          !Temporal smoothing width
      MHV_str(nMHVstr)%nsivel = nsivel          !Temporal smoothing iterations
      MHV_str(nMHVstr)%nssivel = nssivel  !Spatial smoothing iterations
      MHV_str(nMHVstr)%nsswvel = nsswvel !Spatial smoothing window width

    case(7) !NH - Nested WSE
      call nestwse_init
      NH_str(nNHstr)%bidfile = bidfile
      NH_str(nNHstr)%bidpath = bidpath
      NH_str(nNHstr)%idnum = idnum
      NH_str(nNHstr)%istrtype = istrtype
      NH_str(nNHstr)%wseoffset = wseoffset
      NH_str(nNHstr)%ntiwse = ntiwse
      NH_str(nNHstr)%wseout = wseout
      if(wseout)then
        if(len_trim(wsefile)==0)then  
          write(wsefile,'(A,I1,A)')  'Nest_bnd_wse',idnum,'.tsd'
        endif  
        inquire(file=wsefile,exist=foundfile)  
        if(.not.foundfile)then
          wsefile = trim(flowpath) // wsefile
        endif        
        NH_str(nNHstr)%wsefile = wsefile
      endif
      NH_str(nNHstr)%idpar = idpar
      
    case(8) !Nested WSE and Vel
      call nestwsevel_init
      NHV_str(nNHVstr)%bidfile = bidfile
      NHV_str(nNHVstr)%bidpath = bidpath
      NHV_str(nNHVstr)%idnum = idnum
      NHV_str(nNHVstr)%istrtype = istrtype
      NHV_str(nNHVstr)%wseoffset = wseoffset
      NHV_str(nNHVstr)%ntiwse = ntiwse
      NHV_str(nNHVstr)%ntivel = ntivel
      NHV_str(nNHVstr)%wseout = wseout
      NHV_str(nNHVstr)%velout = velout
      if(wseout)then
        if(len_trim(wsefile)==0)then  
          write(wsefile,'(A,I1,A)')  'Nest_bnd_wse',idnum,'.tsd'
        endif  
        inquire(file=wsefile,exist=foundfile)  
        if(.not.foundfile)then
          wsefile = trim(flowpath) // wsefile
        endif
        NHV_str(nNHVstr)%wsefile = wsefile
      endif
      if(velout)then
        if(len_trim(velfile)==0)then  
          write(velfile,'(A,I1,A)')  'Nest_bnd_vel',idnum,'.tsd'
        endif  
        inquire(file=velfile,exist=foundfile)  
        if(.not.foundfile)then
          velfile = trim(flowpath) // velfile
        endif
        NHV_str(nNHVstr)%velfile = velfile
      endif
      NHV_str(nNHVstr)%idpar = idpar
      
    case(9) !Tidal database WSE
      call tidalwse_init
      NTH_str(nNTHstr)%bidfile = bidfile
      NTH_str(nNTHstr)%bidpath = bidpath
      NTH_str(nNTHstr)%idnum = idnum
      NTH_str(nNTHstr)%istrtype = istrtype
      NTH_str(nNTHstr)%tdbname = tdbname
      NTH_str(nNTHstr)%tdbpath = tdbpath
      NTH_str(nNTHstr)%wseoffset = wseoffset
      NTH_str(nNTHstr)%ntcin = ntcin
      if(ntcin>0)then
        allocate(NTH_str(nNTHstr)%namein(ntcin))
        NTH_str(nNTHstr)%namein(:) = namein(:)
      endif
      NTH_str(nNTHstr)%projtdb = projtdb
      NTH_str(nNTHstr)%wseout = wseout
      NTH_str(nNTHstr)%wseadjust = wseadjust
      NTH_str(nNTHstr)%nssi = nssi  !Spatial smoothing iterations
      NTH_str(nNTHstr)%nssw = nssw  !Spatial smoothing window width
      
    case(10) !Tidal database WSE and Vel
      call tidalwsevel_init
      NTHV_str(nNTHVstr)%bidfile = bidfile
      NTHV_str(nNTHVstr)%bidpath = bidpath
      NTHV_str(nNTHVstr)%idnum = idnum
      NTHV_str(nNTHVstr)%istrtype = istrtype
      NTHV_str(nNTHVstr)%tdbname = tdbname
      NTHV_str(nNTHVstr)%tdbpath = tdbpath
      NTHV_str(nNTHVstr)%wseoffset = wseoffset
      NTHV_str(nNTHVstr)%ntcin = ntcin
      if(ntcin>0)then
        allocate(NTHV_str(nNTHVstr)%namein(ntcin))
        NTHV_str(nNTHVstr)%namein(:) = namein(:)
      endif
      NTHV_str(nNTHVstr)%nssi = nssi
      NTHV_str(nNTHVstr)%nssw = nssw
      NTHV_str(nNTHVstr)%projtdb = projtdb
      NTHV_str(nNTHVstr)%wseout = wseout
      NTHV_str(nNTHVstr)%velout = velout
      NTHV_str(nNTHVstr)%wseadjust = wseadjust
      !NTHV_str(nNTHVstr)%nssi = nssi  !Spatial smoothing iterations   !!Repeated
      !NTHV_str(nNTHVstr)%nssw = nssw !Spatial smoothing window width  !!Repeated
      
    case default !Cross-shore boundary condition
          
    end select
    
    !Salinity boundary type
    select case(isaltype)
    case(1) !Constant along boundary
      call salstr_resize
      sal_str(nsalstr)%bidfile = bidfile
      sal_str(nsalstr)%bidpath = bidpath
      !sal_str(nsalstr)%idnum = idnum
      !sal_str(nsalstr)%istrtype = istrtype
      if(len_trim(salfile)==0)then
        salfile = bidfile  
        salpath = bidpath
      endif
      sal_str(nsalstr)%salfile = salfile
      sal_str(nsalstr)%salpath = salpath
      
    !case(2) !Nested (not implimented yet)
    end select
    
    !!Sediment boundary type
    !select case(isedtype)
    !case(1) !Constant along boundary
    !  call sedstr_resize
    !  sed_str(nsedstr)%bidfile = bidfile
    !  sed_str(nsedstr)%bidpath = bidpath
    !  sed_str(nsedstr)%idnum = idnum
    !  !sal_str(nsalstr)%istrtype = istrtype
    !  if(len_trim(salfile)==0)then
    !    sedfile = bidfile  
    !    sedpath = bidpath
    !  endif
    !  sed_str(nsedstr)%sedfile = sedfile
    !  sed_str(nsedstr)%sedpath = sedpath
    !!case(2) !Nested (not implimented yet)
    !end select
    
    return
    end subroutine bnd_block
        
!*******************************************************************************    
    subroutine check_parent_files(grdfilepar,ctlfilepar)
         !wsefilepar,wsepathpar,velfilepar,velpathpar)
!*******************************************************************************  
    use comvarbl, only: flowpath
    use diag_lib, only: diag_print_error
    
    implicit none
    character(len=*),intent(inout) :: grdfilepar,ctlfilepar
    !character(len=*),intent(inout) :: wsefilepar,wsepathpar
    !character(len=*),intent(inout),optional :: velfilepar,velpathpar
    !Internal variables
    logical :: foundfile
    character :: apath*200,aname*100,aext*10
    
    !Grid File
    foundfile = .false.
    if(len_trim(grdfilepar)==0 .and. len_trim(ctlfilepar)==0)then
      call diag_print_error('Must specify at the parent grid or control file')
    endif
    
    if(len_trim(grdfilepar)==0)then !None specified
      call fileparts(ctlfilepar,apath,aname,aext) 
      if(aext(1:2)=='15' .or. aext(1:3)=='ctl')then        
        grdfilepar = trim(apath) // trim(aname) // '.14'  
      elseif(aext(1:7)=='cmcards')then
        grdfilepar = trim(apath) // trim(aname) // '_grid.h5'
      else
        call diag_print_error('Invalid Parent Control File Name: ',ctlfilepar)
      endif
      inquire(file=grdfilepar,exist=foundfile)  
      if(.not.foundfile .and. (aext(1:2)=='15' .or. aext(1:3)=='ctl'))then
        grdfilepar = apath // aname // '.grd'
      endif  
    endif
    inquire(file=grdfilepar,exist=foundfile)
    if(.not.foundfile)then
      call diag_print_error('Could not find Parent Grid File:', grdfilepar)    
    endif
    
    !Control file
    foundfile = .false.  
    if(len_trim(ctlfilepar)==0)then !None specified
      call fileparts(grdfilepar,apath,aname,aext) 
      if(aext(1:2)=='14' .or. aext(1:3)=='grd')then        
        ctlfilepar = trim(apath) // trim(aname) // '.15'  
      elseif(aext(1:2)=='h5')then
        ctlfilepar = trim(apath) // trim(aname) // '.cmcards'
      else
        call diag_print_error('Invalid Parent Grid File Name', grdfilepar)
      endif
      inquire(file=ctlfilepar,exist=foundfile)  
      if(.not.foundfile .and. (aext(1:2)=='14' .or. aext(1:3)=='grd'))then
        ctlfilepar = trim(apath) // trim(aname) // '.ctl'
      endif  
    endif
    inquire(file=ctlfilepar,exist=foundfile) 
    if(.not.foundfile)then
      call diag_print_error('Could not find Parent Control File:', ctlfilepar)
    endif
      
    return
    end subroutine check_parent_files