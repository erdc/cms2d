!===================================================================
module bnd_def
! Flow Boundary Variable Definitions Module
!
! Description:
!   Contains the variable definitions for the hydrodynamic
!   boundary conditions.
! 
! written by Alex Sanchez, USACE-CHL
!===================================================================
    use geo_def,  only: projection
    use prec_def, only: ikind
    implicit none
    save 
    
    !Parameterss
    real(ikind) :: veldamp !Damping factor for velocity at open wse boundaries
    
!--- All Forcing Boundary Strings --------------------------------------------
    integer :: nbndstr                   !Boundary strings
    type bnd_type !Node string driver
      integer          :: idnum          !Boundary id number
      integer          :: ibndtype       !Boundary type,1=Q,2=T,3=H,4=MH,5=MHV,6=CS,7=NH,8=NHV,9=NTH,10-NTHV
      integer          :: ncells         !Cells in string    
      integer          :: nnodes         !Nodes in string
      integer, pointer :: cells(:)       !Cell id's  
      integer, pointer :: nodes(:)       !Node id's
      integer, pointer :: faces(:)       !Boundary face
    endtype bnd_type
    type(bnd_type), allocatable :: bnd_str(:)   

!--- Mesh Boundary Strings --------------------------------------------
!Used to store the boundary strings in the *.2dm file
!Because the node strings are saved at the end of the file in a weird format
!it is easier to read them once, store them, and then reference them later
!Note: Not all of the nodes strings in the 2dm file are necessarily forcing
!strings and they may not be necessarilly all of the strings either.
    integer :: nbcstr                   !Boundary strings
    type bc_type !Node string driver
      integer          :: ncells         !Cells in string    
      integer          :: nnodes         !Nodes in string
      integer, pointer :: cells(:)       !Cell id's  
      integer, pointer :: nodes(:)       !Node id's
      integer, pointer :: faces(:)       !Boundary face
    endtype bc_type
    type(bc_type), allocatable :: bc_str(:)          
    
!--- Wall (Closed) BC ------------------------------------------------------    
! This wall boundary condition is only applied at "outer" wall boundaries
! which are dry faces neighboring internal and outer (ghost or dummy) cells.
! The boundary condition is necessary for the pressure correction equation
    type W_type
      integer              :: ncells      !Cells in string  
      integer, allocatable :: cells(:)    !Cell id's
      integer, allocatable :: faces(:)    !Boundary face
    endtype W_type
    type(W_type) :: W_str !Only one string necessary. Cells do not have to be contiguous
    
!--- Flux (River) BC (Type 1=Q) ------------------------------------------------
    integer :: nQstr                      !Flux cell strings, >=0
    type Q_type
      integer              :: idnum       !Boundary id number
      integer              :: istrtype    !Input string type, 1-cellstring, 2-nodestring
      integer              :: ifluxmode   !1-Constant flux, 2-Total-flux curve, 3-Rating curve, 4-Rating equation
      character(len=200)   :: bidfile     !Boundary ID file
      character(len=200)   :: bidpath     !Boundary ID path
      integer              :: ncells      !Cells in string
      integer,     pointer :: cells(:)    !Cell id's
      integer,     pointer :: faces(:)    !Boundary face
      integer              :: ntimes      !Times in boundary data
      integer              :: inc         !Current time step of boundary data      
      real(ikind)          :: cmvel       !Coefficient [-]
      real(ikind)          :: angle       !Angle of flow clockwise from north [rad]
      logical              :: specified   !True when a user specifies the angle in the parameter file.
      real(ikind)          :: qflux       !Interpolated river flux [m^3/s]
      real(ikind)          :: qfluxconst  !Constant flux [m^3/s]
      real(ikind), pointer :: times(:)    !Times [hrs]
      integer              :: nstages     !Number of stages
      real(ikind), pointer :: stage(:)    !Stage (water level) [m] (used for flow-rating curve)
      real(ikind), pointer :: qcurv(:)    !Data Values (time)
      real(ikind), pointer :: rflow(:)    !Flow-rate [m^3/s] (used for flow-rating curve)
      integer              :: ifluxunits  !Input flux units: 0-m^3/s/cell,1-m^3/s/boundary,2-ft^3/s/boundary
      integer              :: nti         !Order of temporal interpolation (default=2)
      integer              :: nsw         !Smoothing width
      integer              :: nsi         !Smoothing iterations      
      character(len=200) :: fluxfile,fluxpath  !Flux data file and path
    endtype Q_type 
    type(Q_type), allocatable :: Q_str(:)    

!--- Tidal/Harmonic BC (Type 2=TH) --------------------------------------------------------
    integer, parameter :: ntf = 37    !Total number of possible tidal constituents    
    logical Tread,tide_read,harm_read
    integer :: nTHstr                 !Tidal cell strings, either 0 or 1    
    type TH_type       
      integer                       :: idnum       !Boundary id number 
      integer                       :: istrtype    !Input id type, 1-cellstring, 2-nodestring
      integer                       :: ncells      !Cells in string
      character(len=200)            :: bidfile     !Boundary ID file
      character(len=200)            :: bidpath     !Boundary ID path
      integer,              pointer :: cells(:)    !Cell id's (cell)
      integer,              pointer :: faces(:)    !Boundary face (cell)
      logical                       :: istidal     !true for tidal, false for harmonic
      character(len=100)            :: station     !Name of station
      integer                       :: ntc         !Tidal constituents used
      integer                       :: nti         !Order of temporal interpolation (default=2) !hli(10/6/17)
      integer                       :: inc         !Current time step of boundary data  !hli(10/6/17)
      real(ikind),      allocatable :: amp(:)      !Amplitude [m] (constituent) 
      real(ikind),      allocatable :: speed(:)    !Speed [rad/hrs] (constituent)
      real(ikind),      allocatable :: phase(:)    !Phase [rad] (constituent)
      real(ikind),      allocatable :: f(:)        !Nodal factor [-] (constituent)
      real(ikind),      allocatable :: vu(:)       !Equilibrium argument [rad] (constituent)
      character(len=10),allocatable :: name(:)     !Tidal Consitituent names (constituent)
      real(ikind)                   :: angle       !Incident angle of tidal wave
      logical                       :: specified   !True when a user specifies the angle in the parameter file.
      real(ikind),      allocatable :: psi(:,:)    !Phase difference due to incident wave angle (cell,constituent) 
      real(ikind)                   :: dwsex       !Regional steady water level gradinet
      real(ikind)                   :: dwsey       !Regional steady water level gradinet
      real(ikind)                   :: wseoffset   !wse offset [m]
      real(ikind),      allocatable :: wsebnd0(:)  !Initial wse at each cell (may be different from forcing value) [m] (cell)      
      real(ikind),      allocatable :: wsebnd(:)   !Wse at each cell (may be different from forcing value) [m]
      real(ikind),      allocatable :: wsevar(:)   !Temporally constant spatially variation along boundary [m] (cell)
      real(ikind),      allocatable :: wseadj(:)   !Current wse adjusted for wave and wind setup [m] (cell)
      logical                       :: wseadjust   !Turns on or off the wse adjustment/correction due to wind and waves 
      character(len=200)            :: offsetfile  !Offset file           (hli,10/04/17)
      character(len=200)            :: offsetpath  !Offset path
      integer                       :: ioffsetmode !1-Constant offset, 2-Offset curve 
      integer                       :: ntimesoffset!Times in offset curve 
      real(ikind)                   :: wsecurveoffset !Interpolated offset Values [m] 
      real(ikind),          pointer :: offsettimes(:) !Offset times [hrs] 
      real(ikind),          pointer :: offsetcurve(:) !Input offset Values [m]           
    endtype TH_type   
    type(TH_type), allocatable :: TH_str(:)
   
!--- Single Water Level BC (Type 3=H) -------------------------------------------------
    integer :: nHstr
    type H_type 
      integer                  :: idnum      !Boundary id number
      integer                  :: istrtype   !Input id type, 1-cellstring, 2-nodestring
      character(len=200)       :: bidfile    !Boundary ID file
      character(len=200)       :: bidpath    !Boundary ID path
      integer                  :: ncells     !Cells in string
      integer,         pointer :: cells(:)   !Cell id's
      integer,         pointer :: faces(:)   !Boundary face
      integer                  :: ntimes     !Times in boundary data
      integer                  :: inc        !Current time step of boundary data
      real(ikind)              :: wseconst   !Constant wse [m]
      real(ikind)              :: wseoffset  !wse offset
      real(ikind)              :: dwsex      !Regional steady water level gradinet
      real(ikind)              :: dwsey      !Regional steady water level gradinet
      real(ikind), allocatable :: wsebnd0(:) !Initial wse at each cell (may be different from forcing value) [m]
      real(ikind), allocatable :: wsebnd(:)  !Wse at each cell (may be different from forcing value) [m]
      real(ikind), allocatable :: wsevar(:)  !Temporally constant spatially variation along boundary [m] (cell)
      real(ikind), allocatable :: wseadj(:)  !Current wse adjusted for wave and wind setup [m]      
      real(ikind),     pointer :: times(:)   !Times [hrs]
      real(ikind),     pointer :: wsecurv(:) !Data Values (time) [m]      
      integer                  :: minterp    !Method for interpolation, 1-Piecewise polynomial, 2-cubic spline
      real(ikind),     pointer :: d2wse(:)   !Second-derivative used for spline interpolation
      integer                  :: nti        !Order of temporal interpolation (default=2)
      integer                  :: nsw        !Temporal smoothing width
      integer                  :: nsi        !Temporal smoothing iterations
      character(len=200)       :: wsefile    !Water level data file
      character(len=200)       :: wsepath    !Water level data file
      logical                  :: wseadjust  !Turns on or off the wse adjustment/correction due to wind and waves
      character(len=200)       :: offsetfile !Offset file (hli,01/19/17)
      character(len=200)       :: offsetpath !Offset path     
      integer                  :: ioffsetmode !1-Constant offset, 2-Offset curve 
      integer                  :: ntimesoffset!Times in offset curve 
      real(ikind)              :: wsecurveoffset !Interpolated offset Values [m] 
      real(ikind),     pointer :: offsettimes(:) !Offset times [hrs] 
      real(ikind),     pointer :: offsetcurve(:) !Input offset Values [m]           
    endtype H_type    
    type(H_type), allocatable :: H_str(:)    
           
!--- Multiple Water Level BC (Type 4=MH) -------------------------------------
    integer :: nMHstr
    type MH_type
      integer                  :: idnum        !Boundary id number
      integer                  :: istrtype     !Input id type, 1-cellstring, 2-nodestring
      character(len=200)       :: bidfile      !Boundary ID file
      character(len=200)       :: bidpath      !Boundary ID path
      integer                  :: ncells       !Cells in string
      integer,         pointer :: cells(:)     !Cell id's
      integer,         pointer :: faces(:)     !Boundary face
      integer                  :: ntimes       !Times in boundary data
      integer                  :: inc          !Current time step of boundary data      
      real(ikind)              :: wseoffset    !Water level offset which can be used to convert between vertical datums
      real(ikind), allocatable :: wsebnd(:)    !Interpolated eta at each cell at current time [m]
      real(ikind), allocatable :: wsebnd0(:)   !Initial wse at each cell (may be different from forcing) [m]
      real(ikind),     pointer :: times(:)     !Times [hrs]
      real(ikind),     pointer :: wsedata(:,:) !Data Values (time,cell)      
      integer                  :: nti          !Order of temporal interpolation (default=2)
      integer                  :: nsi          !Temporal smoothing iterations
      integer                  :: nsw          !Temporal smoothing width
      integer                  :: nssi         !Spatial smoothing iterations
      integer                  :: nssw         !Spatial smoothing window width
      character(len=200)       :: wsefile,wsepath !Water level data file and path
    endtype MH_type
    type(MH_type), allocatable :: MH_str(:)
    
!--- Multiple Water Level and Velocity BC (Type 5=MHV) -------------------------
    integer :: nMHVstr   
    type MHV_type !Multiple eta, u and v values, one for each cell in string
      integer                  :: idnum        !Boundary id number
      integer                  :: istrtype     !Input id type, 1-cellstring, 2-nodestring
      character(len=200)       :: bidfile      !Boundary ID file
      character(len=200)       :: bidpath      !Boundary ID path
      integer                  :: ncells       !Cells in string
      integer,         pointer :: cells(:)     !Cell id's
      integer,         pointer :: faces(:)     !Boundary face
      integer                  :: ntimeswse       !Times in boundary data
      integer                  :: incwse          !Current time step of boundary data      
      real(ikind)              :: wseoffset    !Water level offset which can be used to convert between vertical datums
      real(ikind), allocatable :: wsebnd(:)    !Interpolated eta at each cell at current time [m]
      real(ikind), allocatable :: wsebnd0(:)   !Initial wse at each cell (may be different from forcing) [m]
      integer                  :: ntimesvel       !Times in boundary data
      integer                  :: incvel          !Current time step of boundary data 
      real(ikind), allocatable :: ubnd(:)      !Interpolated u at each cell at current time [m/s]
      real(ikind), allocatable :: ubnd0(:)     !Initial u at each cell (may be different from forcing) [m/s]
      real(ikind), allocatable :: vbnd(:)      !Interpolated u at each cell at current time [m/s]
      real(ikind), allocatable :: vbnd0(:)     !Initial v at each cell (may be different from forcing) [m/s]
      real(ikind),     pointer :: timeswse(:)     !Times [hrs]
      real(ikind),     pointer :: wsedata(:,:) !Data Values (time,cell)
      real(ikind),     pointer :: timesvel(:)     !Times [hrs]
      real(ikind),     pointer :: udata(:,:)   !Data Values (time,cell)
      real(ikind),     pointer :: vdata(:,:)   !Data Values (time,cell)
      integer                  :: ntiwse          !Order of temporal interpolation (default=2) (same for vel and wse)
      integer                  :: ntivel          !Order of temporal interpolation (default=2) (same for vel and wse)
      integer                  :: nswwse          !Temporal smoothing width
      integer                  :: nsiwse          !Temporal smoothing iterations
      integer                  :: nswvel          !Temporal smoothing width
      integer                  :: nsivel          !Temporal smoothing iterations
      integer                  :: nssiwse         !Spatial smoothing iterations
      integer                  :: nsswwse         !Spatial smoothing window width
      integer                  :: nssivel         !Spatial smoothing iterations
      integer                  :: nsswvel         !Spatial smoothing window width
      character(len=200)       :: wsefile,wsepath !Input Water level data file and path
      character(len=200)       :: velfile,velpath !Input Velocity data file and path
    endtype MHV_type
    type(MHV_type), allocatable :: MHV_str(:) 

!--- Cross-shore BC (Type 6=CS) ----------------------------------------
    integer :: nCSstr
    type CS_type 
      integer                  :: idnum     !Boundary id number
      integer                  :: istrtype  !Input id type, 1-cellstring, 2-nodestring
      character(len=200)       :: bidfile   !Boundary ID file
      character(len=200)       :: bidpath   !Boundary ID path
      integer                  :: ncells    !Cells in string
      integer,         pointer :: cells(:)  !Cell id's
      integer,         pointer :: faces(:)  !Boundary face     
      real(ikind), allocatable :: wsecsh(:) !wse along cell string
      real(ikind), allocatable :: ucsh(:)   !u-velocity along cell string
      real(ikind), allocatable :: vcsh(:)   !v-velocity along cell string
    endtype CS_type
    type(CS_type), allocatable :: CS_str(:)    
    
!--- Parent Simulation ---------------------------    
    integer nParSim
    type ParSim_type
      integer                  :: ipartype      !0-none, 1-CMS, 2-ADCIRC   
      logical                  :: velpar        !Turns on and off the velocity  
      integer                  :: ntiwsepar     !Highest temporal interpolation order (1-3)
      integer                  :: ntivelpar     !Highest temporal interpolation order (1-3)
      integer                  :: incwsepar     !Current time step of boundary data
      integer                  :: incvelpar     !Current time step of boundary data
      real(ikind)              :: t2hrs         !Conversion factor from input time series to hours
      real(ikind)              :: tjuldaypar    !Parent simulation reference (starting) time in Julian days
      real(ikind)              :: timestarthr   !Time relative to the CMS-Flow starting time (used for ADCIRC)
      real(ikind), allocatable :: timewsehrspar(:) !wse times
      real(ikind), allocatable :: timevelhrspar(:) !velocity times
      real(ikind), allocatable :: wsepar(:,:)   !wse data
      real(ikind), allocatable :: upar(:,:)     !U-velocity data (pt,time)
      real(ikind), allocatable :: vpar(:,:)     !V-velocity data (pt,time)
      integer                  :: nptspar       !Size of parent grid output points    
      integer,         pointer :: activepar(:)  !Size of parent grid output points
      character(len=200)       :: ctlfilepar    !Parent control file and path
      character(len=200)       :: grdfilepar    !Parent grid file and path
      character(len=200)       :: wsefilepar,wsepathpar  !Water level data file and path
      character(len=200)       :: velfilepar,velpathpar  !Velocity data file and path
      character(len=200)       :: typespathpar  !Types path for XMDF grid  
      logical                  :: telpargrd     !Telescoping grid
      type(projection)         :: projpar       !Parent grid projection
      integer                  :: ncellspar     !Active grid size of parent grid   
      integer                  :: nmaxfacespar  !Maximum number of faces
      real(ikind)              :: xoriginpar    !Global x-origin of parent grid (if Cartesian)
      real(ikind)              :: yoriginpar    !Global y-origin of parent grid (if Cartesian)   
      real(ikind)              :: orientpar     !Orientation (grid angle) of parent grid
      integer,         pointer :: c2cpar(:,:)   !Cell-to-cells connectivity
      integer,         pointer :: idfpar(:,:)   !Direction of cell face (1-N,2-E,3-South,4-West)
      integer,         pointer :: ncfpar(:)     !Number of cell faces on each cell
      integer                  :: nelemsfullpar !Active grid size of parent grid    
      integer,         pointer :: elem2node(:,:) !Element to node connectivity
      real(ikind),     pointer :: dxpar(:),dypar(:)  !Cell dimensiones of full parent grid 
      real(ikind),     pointer :: xpar(:),ypar(:),zpar(:)   !Global coordinates and elevation of full parent grid       
    endtype ParSim_type
    type(ParSim_type), allocatable :: ParSim(:)
    
!--- Nested Water Level BC (Type 7=NH) -----------------------------------------
    integer :: nNHstr
    type NH_type !Nested eta values, one for each cell in string
      integer                  :: idnum         !Boundary id number
      integer                  :: istrtype      !Input id type, 1-cellstring, 2-nodestring
      character(len=200)       :: bidfile       !Boundary ID file
      character(len=200)       :: bidpath       !Boundary ID path
      integer                  :: ncells        !Cells in string      
      integer,         pointer :: cells(:)      !Cell id's
      integer,         pointer :: faces(:)      !Boundary face  
      real(ikind), allocatable :: xbnd(:)       !Global x-coordinate of boundary cell centroid
      real(ikind), allocatable :: ybnd(:)       !Global x-coordinate of boundary cell centroid
      integer                  :: ntiwse        !Temporal interpolation order (1-3)
      real(ikind), allocatable :: timewsehrs(:) !Times for data [hrs] (time)
      real(ikind)              :: wseoffset     !WSE offset [m]
      real(ikind), allocatable :: wsebnd0(:)    !Interpolated eta at each cell at initial time [m] (cell)
      real(ikind), allocatable :: wsebnd(:)     !Interpolated eta at each cell at current time [m] (cell)
      real(ikind), allocatable :: wsedata(:,:)  !WSE data [m] (time,cell)      
      integer                  :: mntp          !Maximum number of interpolation points (used for allocation)
      integer,         pointer :: intp(:,:)     !Interpolation id's (local grid id, interp mesh id's)
      real(ikind),     pointer :: cntp(:,:)     !Interpolation coefficients
      integer                  :: nssiwse       !Spatial smoothing iterations
      integer                  :: nsswwse       !Spatial smoothing window width
      logical                  :: wseout        !Output wse time-series 
      character(len=200)       :: wsefile       !Output wse file
      integer                  :: idpar         !Parent ID number
    endtype NH_type
    type(NH_type), allocatable :: NH_str(:)
    
!--- Nested Water Level and Velocity BC (Type 8=NHV) ---------------------------------
    integer :: nNHVstr
    type NHV_type !Nested eta, u and v values, one for each cell in string
      integer                  :: idnum         !Boundary id number
      integer                  :: istrtype      !Input id type, 1-cellstring, 2-nodestring
      character(len=200)       :: bidfile       !Boundary ID file
      character(len=200)       :: bidpath       !Boundary ID path
      integer                  :: ncells        !Cells in string      
      integer,         pointer :: cells(:)      !Cell id's
      integer,         pointer :: faces(:)      !Boundary face
      real(ikind), allocatable :: xbnd(:)       !Global x-coordinate of boundary cell centroid
      real(ikind), allocatable :: ybnd(:)       !Global x-coordinate of boundary cell centroid
      integer                  :: ntiwse        !Temporal interpolation order (1-3)
      integer                  :: ntivel        !Temporal interpolation order (1-3)
      real(ikind)              :: wseoffset     !WSE offset [m]
      real(ikind), allocatable :: wsebnd0(:)    !Interpolated wse at each cell at initial time [m]
      real(ikind), allocatable :: wsebnd(:)     !Interpolated wse at each cell at current time [m]
      real(ikind), allocatable :: wsedata(:,:)  !WSE data [m] (time,cell)   
      real(ikind), allocatable :: ubnd0(:)      !Interpolated u at each cell at initial time [m/s]
      real(ikind), allocatable :: ubnd(:)       !Interpolated u at each cell at current time [m/s]
      real(ikind), allocatable :: udata(:,:)    !U data [m] (time,cell)
      real(ikind), allocatable :: vbnd0(:)      !v at each cell at initial time [m/s]
      real(ikind), allocatable :: vbnd(:)       !Interpolated u at each cell at current time [m/s]      
      real(ikind), allocatable :: vdata(:,:)    !V data [m] (time,cell)      
      real(ikind), allocatable :: timewsehrs(:) !Times for wse data [hrs] (time)
      real(ikind), allocatable :: timevelhrs(:) !Times for velocity data [hrs] (time)
      integer                  :: mntp          !Maximum number of interpolation points (used for allocation)
      integer,         pointer :: intp(:,:)     !Interpolation id's (local grid id, interp mesh id's)
      real(ikind),     pointer :: cntp(:,:)     !Interpolation coefficients      
      integer                  :: nssiwse       !Smoothing iterations
      integer                  :: nsswwse       !Smoothing window width
      integer                  :: nssivel       !Smoothing iterations
      integer                  :: nsswvel       !Smoothing window width
      real(ikind)              :: angvel        !Rotation angle for current velocities from parent to child grid [rad]
      logical                  :: wseout        !Output wse time-series      
      logical                  :: velout        !Output Velocity time-series 
      character(len=200)       :: wsefile       !Output wse file
      character(len=200)       :: velfile       !Output velocity file
      integer                  :: idpar         !Parent ID number
    endtype NHV_type  
    type(NHV_type), allocatable :: NHV_str(:)          
    
!--- Nested Tidal Database WSE BC (Type 9=NTH) -------------------------------------------
    integer :: nNTHstr
    type NTH_type !Nested eta values, one for each cell in string
      integer                   :: idnum         !Boundary id number
      integer                   :: istrtype      !Input id type, 1-cellstring, 2-nodestring
      character(len=200)        :: bidfile       !Boundary ID file
      character(len=200)        :: bidpath       !Boundary ID path
      integer                   :: ncells        !Cells in string      
      integer,         pointer  :: cells(:)      !Cell id's
      integer,         pointer  :: faces(:)      !Boundary face
      real(ikind), allocatable  :: xbnd(:)       !Global x-coordinate of boundary cell centroid
      real(ikind), allocatable  :: ybnd(:)       !Global x-coordinate of boundary cell centroid
      real(ikind)               :: wseoffset     !wse offset
      real(ikind), allocatable  :: wsebnd(:)     !Interpolated eta at each cell at current time [m/s]
      real(ikind), allocatable  :: wseadj(:)     !Adjusted eta for wave and wind setup at each cell at current time [m/s]
      real(ikind), allocatable  :: wsebnd0(:)    !Interpolated eta at each cell at initial time [m/s]
      integer                   :: ntcin         !Tidal constituents used     
      character(len=10),pointer :: namein(:)     !Input Tidal Consitituent names (constituent)
      integer                   :: ntc           !Tidal constituents used      
      real(ikind),     pointer  :: amp(:,:)      !Amplitude [m] (cell,constituent)       
      real(ikind),     pointer  :: phase(:,:)    !Phase [rad] (cell,constituent)
      real(ikind),     pointer  :: f(:)          !Nodal factor [-] (constituent)
      real(ikind),     pointer  :: vu(:)         !Equilibrium argument [rad] (constituent)
      real(ikind),     pointer  :: speed(:)      !Speed [rad/hrs] (constituent)
      character(len=10),pointer :: name(:)       !Tidal Consitituent names (constituent) 
      character(len=10)         :: tdbname       !Tidal Database Name, EC2001, ENPAC2003, LEPROVOST
      character(len=200)        :: tdbpath       !Tidal Database file and path      
      type(projection)          :: projtdb       !Parent grid projection
      logical                   :: wseout        !Output wse time-series 
      character(len=200)        :: wsefile       !Output wse file
      integer                   :: nssi          !Spatial smoothing iterations
      integer                   :: nssw          !Spatial smoothing window width
      logical                   :: wseadjust     !Turns on or off the wse adjustment/correction due to wind and waves
    endtype NTH_type
    type(NTH_type), allocatable :: NTH_str(:)    

!--- Nested Tidal Database WSE and Velocity BC (Type 10=NTHV) ----------------------------------
    integer :: nNTHVstr
    type NTHV_type !Nested eta values, one for each cell in string
      integer                   :: idnum         !Boundary id number
      integer                   :: istrtype      !Input id type, 1-cellstring, 2-nodestring
      character(len=200)        :: bidfile       !Boundary ID file
      character(len=200)        :: bidpath       !Boundary ID path
      integer                   :: ncells        !Cells in string      
      integer,         pointer  :: cells(:)      !Cell id's
      integer,         pointer  :: faces(:)      !Boundary face
      real(ikind), allocatable  :: xbnd(:)       !Global x-coordinate of boundary cell centroid
      real(ikind), allocatable  :: ybnd(:)       !Global x-coordinate of boundary cell centroid
      real(ikind)               :: wseoffset     !wse offset
      real(ikind), allocatable  :: wsebnd0(:)    !Interpolated eta at each cell at initial time [m/s]
      real(ikind), allocatable  :: wseadj(:)     !Adjusted eta for wave and wind setup at each cell at current time [m/s]
      real(ikind), allocatable  :: wsebnd(:)     !Interpolated eta at each cell at current time [m/s]            
      real(ikind), allocatable  :: ubnd0(:)      !v at each cell at initial time [m/s]
      real(ikind), allocatable  :: ubnd(:)       !Interpolated eta at each cell at current time [m/s]
      real(ikind), allocatable  :: vbnd0(:)      !v at each cell at initial time [m/s]
      real(ikind), allocatable  :: vbnd(:)       !Interpolated eta at each cell at initial time [m/s]
      integer                   :: ntcin         !Tidal constituents used     
      character(len=10),pointer :: namein(:)     !Input Tidal Consitituent names (constituent)
      integer                   :: ntc           !Tidal constituents used      
      real(ikind),     pointer  :: amp(:,:)      !Amplitude [m | m/s] (cell,constituent)
      real(ikind),     pointer  :: ampu(:,:)     !Amplitude [m | m/s] (cell,constituent)
      real(ikind),     pointer  :: ampv(:,:)     !Amplitude [m | m/s] (cell,constituent)      
      real(ikind),     pointer  :: phase(:,:)    !Phase [rad] (cell,constituent)
      real(ikind),     pointer  :: phaseu(:,:)   !Phase [rad] (cell,constituent)
      real(ikind),     pointer  :: phasev(:,:)   !Phase [rad] (cell,constituent)
      real(ikind),     pointer  :: f(:)          !Nodal factor [-] (constituent)
      real(ikind),     pointer  :: vu(:)         !Equilibrium argument [rad] (constituent)
      real(ikind),     pointer  :: speed(:)      !Speed [rad/hrs] (cell,constituent)
      character(len=10),pointer :: name(:)       !Tidal Consitituent names (constituent)
      character(len=10)         :: tdbname       !Tidal Database Name, EC2001, ENPAC2003, LEPROVOST, 
      character(len=200)        :: tdbpath       !Tidal Database file and path
      type(projection)          :: projtdb       !Parent grid projection
      logical                   :: wseout        !Output wse time-series       
      logical                   :: velout        !Output velocity time-series 
      character(len=200)        :: wsefile       !Output wse file
      character(len=200)        :: velfile       !Output velocity file
      integer                   :: nssi          !Smoothing iterations for constituents
      integer                   :: nssw          !Smoothing window width for constituents
      real(ikind)               :: angvel        !Rotation angle for current velocities from parent to child grid [rad]
      logical                   :: wseadjust     !Turns on or off the wse adjustment/correction due to wind and waves
    endtype NTHV_type
    type(NTHV_type), allocatable :: NTHV_str(:)
    
    integer, allocatable :: bndcorner(:)
    integer              :: nbndcorner,ioffsetmode !hli
    
end module bnd_def