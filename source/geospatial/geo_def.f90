!===================================================================
module size_def
! Grid size module
!===================================================================
    implicit none
    save
    integer:: ncellsfull   !Full cells or elements (Total grid cells active and inactive, but not ghost cells)
    integer:: ncells       !Active cells or elements
    integer:: ncellsD      !Active cells plus dummy or ghost cells
    integer:: nnodes       !# of nodes
    integer:: ncellsimple  !Regular Cartesian cells
    integer:: ncelljoint   !Telescoping or joint cells
    integer:: ncellpoly    !Polygonal cells
    integer:: nmaxfaces    !Maximum # of faces per cell
    !integer:: ncmaxcells  !Maximum # of neighbor cells per cell (including extended stencil
    integer:: ndmaxfaces   !Maximum # of faces per dimension per cell
    integer:: nmaxcells    !Maximum # of cells per node
    integer:: nmaxextcells !Maximum # of cells in extended stencil
endmodule size_def

!===================================================================
module geo_def
! CMS Geopatial variable definitions
! Author: Alex Sanchez, USACE-CHL
!===================================================================
    use prec_def     
    implicit none    
    save  
    
    !Files
    character(len=200) :: grdfile = " ",grdpath = " "      !Grid file and path
    character(len=200) :: telfile              !Telescoping grid file
    character(len=200) :: latfile,latpath      !Latitude file name and path
    character(len=200) :: lonfile,lonpath      !Longitude file name and path
    character(len=200) :: typespath,proppath,rootpath
    character(len=200) :: meshfile             !From Chris' code 10/26/2015
    !Time-varying depth datasets (at every grid cell)
    type bathydatatype
      logical                :: ison
      real(ikind)            :: tjulday !Depth reference (starting) time in Julian days
      integer                :: inc2        !Index of times corresponding to depth2 (i.e. depth2->times(inc2))
      integer                :: ntimes      !# of dataset times
      !real(ikind),pointer    :: times(:)
      character(len=200)     :: file
      character(len=200)     :: path
      real(ikind)            :: time0
      real(ikind)            :: time1
      real(ikind)            :: time2
      real(ikind)            :: timen
      !real(ikind), allocatable :: depth0(:)           !Water depths at cell centers (positive is downwards)
      real(ikind), allocatable :: depth1(:)           !Water depths at cell centers (positive is downwards)
      real(ikind), allocatable :: depth2(:)           !Water depths at cell centers (positive is downwards)
    endtype
    type(bathydatatype):: bathydata
    
    !Connectivity (general)
    integer,allocatable:: ncface(:)                !# of cell faces
    integer,allocatable:: cell2cell(:,:)           !Forward connectivity table to neighest neighbors
    integer,allocatable:: llec2llec(:,:)           !Backward connectivity table 
    integer,allocatable:: nxface(:),kxface(:,:)    !No repeat cell faces in x-direction only
    integer,allocatable:: nyface(:),kyface(:,:)    !No repeat cell faces in y-direction only
    integer,allocatable:: nxyface(:),kxyface(:,:)  !No repeat cell faces in x and y directions
    integer,allocatable:: idirface(:,:)            !Direction of cell face 1-North,2-East,3-South,4-West
    integer,allocatable:: idcellsimple(:)          !Non-joint or regular Cartesian cells
    integer,allocatable:: idcelljoint(:)           !Joint (telescoping) Cartesian cells
    integer,allocatable:: nncell(:)                !Number of cells neighboring each node
    integer,allocatable:: cell2node(:,:)           !Cell to node connectivity
    integer,allocatable:: node2cell(:,:)           !Node to cell connectivity
    integer,allocatable:: cell2upwdcell(:,:)       !Second upwind connectivity exists
    integer,allocatable:: ncnode(:)
    integer,allocatable:: cell2extcell(:,:)        !Forward connectivity for extended stencil (node
         
    !--- Geometry ---------------
    !Coordinates
    real(ikind),allocatable:: x(:),y(:)    !Cell centroid local coordinates [m]
    real(ikind),allocatable:: xc(:),yc(:)  !Cell centroid x and y coordinates of flow grid on the global coordinate system [m]
    real(ikind),allocatable:: xn(:),yn(:)  !Node x and y coordinates of flow grid on the global coordinate system [m]
    real(ikind) :: avg_lat     !Average grid latitude (degrees) (used for Coriolis calculation)
    
    !Cell-dimensions
    real(ikind),allocatable:: dx(:),dy(:)           !Grid resolution for Cartesian grids [m]
    real(ikind),allocatable:: areap(:)               !Cell area [m^2]
    real(ikind) :: areaavg !Average cell area
    real(ikind) :: areamin !Minimum cell area
    real(ikind) :: areamax !Maximum cell area
    
    !Unit Vectors
    real(ikind),allocatable:: fnx(:,:),fny(:,:) !Cell face normal unit vectors
    
    !Distance Vectors
    real(ikind),allocatable:: ds(:,:)           !Cell face/edge length [m]
    real(ikind),allocatable:: dc(:,:)           !Distance from centroid to centroid =sqrt((x(i)-x(nck))**2+(y(i)-y(nck))**2) [m]
    real(ikind),allocatable:: dn(:,:)           !Normal distance from centroid to centroid
    real(ikind),allocatable:: dnx(:,:),dny(:,:) !Normal distance vectors from centroid to centroid
    real(ikind),allocatable:: dsxy(:,:)         !used for diffusion terms
    real(ikind),allocatable:: dsx(:,:),dsy(:,:) !ds(k,i)*[fnx(k,i),fny(k,i)] [m]
    real(ikind),allocatable:: rnx(:,:),rny(:,:) !Normal distance vector from O to f [m]
    real(ikind),allocatable:: rpx(:,:),rpy(:,:) !Parallel distance vector from P to O [m]
    real(ikind),allocatable:: rx(:,:), ry(:,:)  !Short distance vector from P to f  [m]
    real(ikind),allocatable:: dpara(:,:),dnorm(:,:)         !From Chris' code 10/26/2015
    
    !Bathymetry
    real(ikind), allocatable :: zb0(:)          !Initial bed elevations at cell centers (negative is downwards)
    real(ikind), allocatable :: zb(:)           !Bed elevations at cell centers (negative is downwards)
    real(ikind), allocatable :: zbk(:,:)        !Bed elevation at cell faces (negative is downwards)
    real(ikind), allocatable :: zbn(:)          !Bed elevation at nodes (negative is downwards)
    real(ikind), allocatable :: dzbx(:),dzby(:) !Bedslopes
    
    !Grid properties
    integer :: nflowgrd       !Flow grid modification number
    integer :: maxrow,maxcol  !Maximum number of rows and columns for regular and nonuniform Cartesian grids
    integer :: igridtype      !0-Non-telescoping, 1-Telescoping, 2-Polygonal, 3-Curvilinear    
    real(ikind) :: azimuth_fl  !Grid orientation (degrees)
    real(ikind) :: xOrigin,yOrigin !Cartesian grid origin (m)
    integer, allocatable :: icol(:),irow(:) 
    real*4, allocatable :: dyy(:),Z(:),dxx(:),lat(:),lon(:) !temporary variables, must be single
    integer, allocatable :: idmap(:),mapid(:)
    integer, allocatable :: cell_type(:)
    integer :: kkface(4)
    real(ikind):: signface(4)
    
    !Projection
    character(len=5)  :: aHorizDatum(0:2)        !Horizontal Datum
    character(len=40) :: aHorizCoordSystem(0:22) !Horizontal Coordinate System
    character(len=9)  :: aHorizUnits(0:6)        !Horizontal Units
    character(len=6)  :: aVertDatum(0:9)         !Vertical Datum
    character(len=6)  :: aVertUnits(1:4)         !Vertical Units
    type projection
      integer :: iHorizDatum        !0-NAD27,1-NAD83,2-Local      
      integer :: iHorizCoordSystem  !Code number for input horizontal coordinate system
      integer :: iHorizZone         !Code number for input horizontal coordinate zone
      integer :: iHorizUnits        !Code number for input horizontal coordinate units
      integer :: iVertDatum         !Code number for input vertical coordinate datum
      integer :: iVertUnits         !Code number for input vertical coordinate units
      real(ikind) :: VertOffset     !Vertical offset from datum
    endtype projection
    type(projection) :: projfl    
    integer :: nzones
    type zonetype
      character(len=200) :: name  !Zone name
      character(len=4)   :: fip   !Zone code
      character(len=10)  :: utm   !UTM zone
    endtype zonetype    
    type(zonetype), allocatable :: zones(:)
    character(len=200) :: HProj,VProj     !added to keep track on projection
    
endmodule geo_def 

