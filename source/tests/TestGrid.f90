MODULE TestGrid

USE ErrorDefinitions
USE XmdfDefs
USE Xmdf

CHARACTER(LEN=*),PARAMETER :: GRID_CART2D_GROUP_NAME = 'Grid Cart2D Group'
CHARACTER(LEN=*),PARAMETER :: GRID_CURV2D_GROUP_NAME = 'Grid Curv2D Group'
CHARACTER(LEN=*),PARAMETER :: GRID_CART3D_GROUP_NAME = 'Grid Cart3D Group'

CONTAINS

!****************************
! ---------------------------------------------------------------------------
! FUNCTION  TG_READ_GRID
! PURPOSE   Read a grid and write data to a text file
! NOTES     
! ---------------------------------------------------------------------------
SUBROUTINE  TG_READ_GRID(a_Id, a_Outfile, error)
INTEGER(XID), INTENT(IN) :: a_Id
INTEGER, INTENT(IN)        :: a_Outfile
INTEGER, INTENT(OUT)       :: error
INTEGER nGridType, nExtrudeType, nDims, nCellsI, nCellsJ
INTEGER nCellsK, nLayers, nOrientation, nValsI, nValsJ
INTEGER nValsK, nExtrudeVals, nCompOrigin, nUDir
CHARACTER(265) strGridType, strExtrudeType
LOGICAL*2 bDefined
REAL(DOUBLE) dOrigin(3), dBearing, dDip, dRoll
REAL(DOUBLE), ALLOCATABLE :: dExtrudeVals(:), dCoordI(:), dCoordJ(:)
REAL(DOUBLE), ALLOCATABLE :: dCoordK(:)
INTEGER i

nGridType = 0
nExtrudeType = 0
nDims = 0
nCellsI = 0
nCellsJ = 0
nCellsK = 0
nLayers = 0
nOrientation = 0
nValsI = 0
nValsJ = 0
nValsK = 0
nExtrudeVals = 0
nCompOrigin = 1
nUDir = 1
bDefined = .FALSE.
dBearing = 0.0
dDip = 0.0
dRoll =0.0
i = 0
error = 1

  ! Grid type
call XF_GET_GRID_TYPE(a_Id, nGridType, error)
if (error < 0) then
  return
endif
SELECT CASE (nGridType) 
  CASE (GRID_TYPE_CARTESIAN)
    strGridType = 'Cartesian'
  CASE (GRID_TYPE_CURVILINEAR)
    strGridType = 'Curvilinear'
  CASE (GRID_TYPE_CARTESIAN_EXTRUDED)
    strGridType = 'Cartesian extruded'
  CASE (GRID_TYPE_CURVILINEAR_EXTRUDED)
    strGridType = 'Curvilinear extruded'
  CASE DEFAULT
    WRITE(*,*) 'Invalid grid type'
    error = -1
END SELECT
WRITE(a_Outfile,*) 'The grid type is: ', strGridType(1:LEN_TRIM(strGridType))

  ! Number of dimensions
call XF_GET_NUMBER_OF_DIMENSIONS(a_Id, nDims, error)
if (error .LT. 0) then
  return
endif
if (nDims == 2) then
  WRITE(a_Outfile,*) 'The grid is two-dimensional'
elseif (nDims == 3) then
  WRITE(a_Outfile,*) 'The grid is three-dimensional'
else
  WRITE(*,*) 'The grid dimensions are invalid'
  error = -1
endif

  ! Extrusion type if applicable
if (nGridType .EQ. GRID_TYPE_CARTESIAN_EXTRUDED .OR. &
    nGridType .EQ. GRID_TYPE_CURVILINEAR_EXTRUDED) then
  call XF_GET_EXTRUSION_TYPE(a_Id, nExtrudeType, error)
  if (error < 0) then
    return
  endif
  SELECT CASE (nExtrudeType)
    case (EXTRUDE_SIGMA)
      strExtrudeType = 'Sigma stretch'
    case (EXTRUDE_CARTESIAN)
      strExtrudeType = 'Cartesian'
    case (EXTRUDE_CURV_AT_CORNERS)
      strExtrudeType = 'Curvilinear at Corners'
    case (EXTRUDE_CURV_AT_CELLS)
      strExtrudeType = 'Curvilinear at Cells'
  END SELECT
  WRITE(a_Outfile,*) 'The grid is extruding using: ', &
                     strExtrudeType(1:LEN_TRIM(strExtrudeType))
endif

  ! Origin
call XF_ORIGIN_DEFINED(a_Id, bDefined, error)
if (error < 0) then
  return
endif
if (bDefined) then
  call XF_GET_ORIGIN(a_Id, dOrigin(1), dOrigin(2), dOrigin(3), error)
  if (error < 0) then
    return
  endif
  WRITE(a_Outfile,*) 'The grid origin is ', dOrigin(1), ' ',&
                      dOrigin(2), ' ', dOrigin(3)
endif
  
  ! Orientation
call XF_GET_ORIENTATION(a_Id, nOrientation, error)
if (error < 0) then
  return
endif
if (nOrientation == ORIENTATION_RIGHT_HAND) then
  WRITE(a_Outfile,*) 'The grid has a right hand orientation'
elseif (nOrientation == ORIENTATION_LEFT_HAND) then
  WRITE(a_Outfile,*) 'The grid has a left hand orientation'
else 
  WRITE(*,*) 'Invalid grid orientation';
  error = -1
  return
endif

  ! Bearing
call XF_BEARING_DEFINED(a_Id, bDefined, error)
if (error < 0) then
  return
endif
if (bDefined) then
  call XF_GET_BEARING(a_Id, dBearing, error)
  if (error < 0) then
    return
  endif
  WRITE(a_Outfile,*) 'The grid bearing is ', dBearing
endif

  ! Dip
call XF_DIP_DEFINED(a_Id, bDefined, error)
if (error < 0) then
  return
endif
if (bDefined) then
  call XF_GET_DIP(a_Id, dDip, error)
  if (error < 0) then
    return
  endif
  WRITE(a_Outfile,*) 'The grid Dip is ', dDip
endif

if (nDims == 3) then
  ! Roll
  call XF_ROLL_DEFINED(a_Id, bDefined, error)
  if (error < 0) then
    return
  endif
  if (bDefined) then
    call XF_GET_ROLL(a_Id, dRoll, error)
    if (error < 0) then
      return
    endif
    WRITE(a_Outfile,*) 'The grid Roll is ', dRoll
  endif
endif

  ! Computational origin
call XF_COMPUTATIONAL_ORIGIN_DEFINED(a_Id, bDefined, error)
if (error < 0) then
  return
endif
if (bDefined) then
  call XF_GET_COMPUTATIONAL_ORIGIN(a_Id, nCompOrigin, error)
  if (error < 0) then
    return
  endif
  WRITE(a_Outfile,*) 'The grid Computational Origin is ', nCompOrigin
else 
  WRITE(a_Outfile,*) 'The grid Computational Origin is not defined';
endif


  ! U Direction
call XF_GET_U_DIRECTION_DEFINED(a_Id, bDefined, error)
if (error < 0) then
  return 
endif
if (bDefined) then
  call XF_GET_U_DIRECTION(a_Id, nUDir, error)
  if (error < 0) then
    return
  endif
    WRITE(a_Outfile,*) 'The grid U Direction is ', nUDir
else 
  WRITE(a_Outfile,*) 'The grid U Direction is not defined'
endif

  ! number of cells in each direction
call XF_GET_NUMBER_CELLS_IN_I(a_Id, nCellsI, error)
if (error >= 0) then
  call XF_GET_NUMBER_CELLS_IN_J(a_Id, nCellsJ, error)
  if ((error >= 0) .AND. (nDims == 3)) then 
    call XF_GET_NUMBER_CELLS_IN_K(a_Id, nCellsK, error)
  endif
endif
if (error < 0) then
  return
endif
WRITE(a_Outfile,*) 'Number of cells in I ', nCellsI
WRITE(a_Outfile,*) 'Number of cells in J ', nCellsJ
if (nDims == 3) then
  WRITE(a_Outfile,*) 'Number of cells in K ', nCellsK
endif

  ! Grid coordinates 
if (nGridType == GRID_TYPE_CARTESIAN .OR. &
    nGridType == GRID_TYPE_CARTESIAN_EXTRUDED) then
  nValsI = nCellsI
  nValsJ = nCellsJ
  if (nDims == 3) then
    nValsK = nCellsK
  endif
elseif (nGridType == GRID_TYPE_CURVILINEAR .OR. &
        nGridType == GRID_TYPE_CURVILINEAR_EXTRUDED) then
  if (nDims == 3) then
      ! three dimensions
    nValsK = (nCellsI + 1) * (nCellsJ + 1) * (nCellsK + 1)
    nValsJ = nValsK
    nValsI = nValsJ
  else
      ! two dimensions
    nValsJ = (nCellsI + 1) * (nCellsJ + 1)
    nValsI = nValsJ
  endif
else
  WRITE(*,*) 'Invalid grid type'
  error = -1
  return
endif

ALLOCATE(dCoordI(nValsI))
ALLOCATE(dCoordJ(nValsJ))
if (nDims == 3) then
  ALLOCATE(dCoordK(nValsK))
endif

call XF_GET_GRID_COORDS_I(a_Id, nValsI, dCoordI, error)
if (error >= 0) then
  call XF_GET_GRID_COORDS_J(a_Id, nValsJ, dCoordJ, error)
  if ((error >= 0) .AND. (nDims == 3)) then
    call XF_GET_GRID_COORDS_K(a_Id, nValsK, dCoordK, error)
  endif
endif
if (error < 0) then
  WRITE(*,*) 'Error reading coordinates'
  error = -1
  return
endif

WRITE(a_Outfile,*) 'The Coordinates in direction I:'
do i = 1, nValsI 
  if (mod(i,5) == 0) then
    WRITE(a_Outfile,*) ''
  endif
  WRITE(a_Outfile,*) dCoordI(i)
enddo
WRITE(a_Outfile,*) ''

WRITE(a_Outfile,*) 'The Coordinates in direction J:'
do i = 1, nValsJ
  if (mod(i,5) == 0) then
    WRITE(a_Outfile,*) ''
  endif
  WRITE(a_Outfile,*) dCoordJ(i)
enddo
WRITE(a_Outfile,*) ''

if (nDims == 3) then
  WRITE(a_Outfile,*) 'The Coordinates in direction K:'
  do i = 1, nValsK
    if (mod(i,5) == 0) then
      WRITE(a_Outfile,*) ''
    endif
    WRITE(a_Outfile,*) dCoordK(i)
  enddo
endif
WRITE(a_Outfile,*) ''

if (ALLOCATED(dCoordI)) DEALLOCATE(dCoordI)
if (ALLOCATED(dCoordJ)) DEALLOCATE(dCoordJ)
if (ALLOCATED(dCoordK)) DEALLOCATE(dCoordK)

!  // Extrude data
if (nGridType .EQ. GRID_TYPE_CARTESIAN_EXTRUDED .OR. &
    nGridType .EQ. GRID_TYPE_CURVILINEAR_EXTRUDED) then
  call XF_GET_EXTRUDE_NUM_LAYERS(a_Id, nLayers, error)
  if (error < 0) then
    return
  endif

  SELECT CASE(nExtrudeType)
    case (EXTRUDE_SIGMA)
      nExtrudeVals = nLayers
    case (EXTRUDE_CURV_AT_CORNERS)
      nExtrudeVals = (nCellsI + 1) * (nCellsJ + 1) * nLayers
    case (EXTRUDE_CURV_AT_CELLS)
      nExtrudeVals = nCellsI * nCellsJ * nLayers
  END SELECT

  ALLOCATE(dExtrudeVals(nExtrudeVals))

  call XF_GET_EXTRUDE_VALUES(a_Id, nExtrudeVals, dExtrudeVals, error)
  if (error < 0) then
    return
  endif

  WRITE(*,*) 'The extrude values are:'
  do i = 1, nExtrudeVals
    if (mod(i,5) == 0) then
      WRITE(a_Outfile,*) ''
    endif
    WRITE(a_Outfile,*) dExtrudeVals(i)
  enddo
  if (ALLOCATED(dExtrudeVals)) DEALLOCATE(dExtrudeVals)
endif

return

END SUBROUTINE TG_READ_GRID

!****************************
!----------------------------------------------------------------------------
! FUNCTION   TG_WRITE_TEST_GRID_CART_2D
! PURPOSE    Write a file that contains data for a 2D Cartesian Grid
! NOTES      A picture of the grid is in the file (TestGridCart2D.gif)
!            returns TRUE on success and FALSE on failure
!----------------------------------------------------------------------------
SUBROUTINE TG_WRITE_TEST_GRID_CART_2D(Filename, error)
use const_def, only: READWRITE, READONLY, OVERWRITE
  
CHARACTER(LEN=*), INTENT(IN) :: Filename
INTEGER, INTENT(OUT) :: error
INTEGER        nDimensions
INTEGER        nCellsI, nCellsJ
INTEGER        nGridType
INTEGER        nCompOrigin, nUDir
REAL(DOUBLE)   dOriginX, dOriginY, dOriginZ
INTEGER        nOrientation
REAL(DOUBLE)   dBearing
REAL(DOUBLE)   PlanesI(5), PlanesJ(5)
INTEGER        i, j, iSpcZone
INTEGER(XID) xFileId, xGridId, xCoordId
INTEGER        status
INTEGER        tmpOut1, tmpOut2

  nDimensions = 2
  nCellsI = 5
  nCellsJ = 5
  nGridType = GRID_TYPE_CARTESIAN
  nCompOrigin = 4
  nUDir = -2
  dOriginX = 10.0
  dOriginY = 10.0
  dOriginZ = 0.0
  nOrientation = ORIENTATION_RIGHT_HAND
  dBearing = 45.0
  xFileId = NONE
  xGridId = NONE

    ! Fill in the grid plane data with a constant size of 30
  do i = 1, nCellsI
    PlanesI(i) = i*30.0
  enddo
  do j = 1, nCellsJ
    PlanesJ(j) = j*30.0
  enddo

    ! create the file
  call XF_CREATE_FILE(Filename, OVERWRITE, xFileId, status)
  if (status < 0) then
    error = -1  
    return
  endif

    ! create the group to store the grid
  call XF_CREATE_GROUP_FOR_GRID(xFileId, GRID_CART2D_GROUP_NAME, xGridId, status)
  if (status < 0) then
    call XF_CLOSE_FILE(xFileId, error)
    error = -1
    return
  endif

    ! Write the grid information to the file
  call XF_SET_GRID_TYPE(xGridId, nGridType, tmpOut1)
  call XF_SET_NUMBER_OF_DIMENSIONS(xGridId, nDimensions, tmpOut2)

  if ((tmpOut1 < 0) .OR. (tmpOut2 < 0)) then
    call XF_CLOSE_GROUP(xGridId, error)
    call XF_CLOSE_FILE(xFileId, error)
    error = -1
    return
  endif

    ! set origin and orientation
  call XF_SET_ORIGIN(xGridId, dOriginX, dOriginY, dOriginZ, tmpOut1) 
  call XF_SET_ORIENTATION(xGridId, nOrientation, tmpOut2)

  if ((tmpOut1 < 0) .OR. (tmpOut2 .LT. 0)) then
    call XF_CLOSE_GROUP(xGridId, error)
    call XF_CLOSE_FILE(xFileId, error)
    error = -1
    return
  endif

    ! Set bearing
  call XF_SET_BEARING(xGridId, dBearing, tmpOut1)
  if (tmpOut1 < 0) then
    call XF_CLOSE_GROUP(xGridId, error)
    call XF_CLOSE_FILE(xFileId, error)
    error = -1
    return
  endif

   ! Set computational origin
  call XF_SET_COMPUTATIONAL_ORIGIN(xGridId, nCompOrigin, tmpOut1)
  if (tmpOut1 < 0) then
    call XF_CLOSE_GROUP(xGridId, error)
    call XF_CLOSE_FILE(xFileId, error)
    error = -1
    return
  endif

    ! Set u direction
  call XF_SET_U_DIRECTION(xGridId, nUDir, tmpOut1)
  if (tmpOut1 < 0) then
    call XF_CLOSE_GROUP(xGridId, error)
    call XF_CLOSE_FILE(xFileId, error)
    error = -1
    return
  endif

    ! Write the grid geometry to the file
    ! Set the number of cells in each direction
  call XF_SET_NUMBER_CELLS_IN_I(xGridId, nCellsI, tmpOut1)
  call XF_SET_NUMBER_CELLS_IN_J(xGridId, nCellsJ, tmpOut2)
  if ((tmpOut1 < 0) .OR. (tmpOut2 < 0)) then
    call XF_CLOSE_GROUP(xGridId, error)
    call XF_CLOSE_FILE(xFileId, error)
    error = -1
    return
  endif

    ! Set the grid plane locations
  call XF_SET_GRID_COORDS_I(xGridId, nCellsI, PlanesI, tmpOut1)
  call XF_SET_GRID_COORDS_J(xGridId, nCellsJ, PlanesJ, tmpOut2)  
  if ((tmpOut1 < 0) .OR. (tmpOut2 < 0)) then
    call XF_CLOSE_GROUP(xGridId, error)
    call XF_CLOSE_FILE(xFileId, error)
    error = -1
    return
  endif

  ! Write Coordinate file - for GridCart2D, we will set the coordinate system
  !   to be State Plane NAD27.
  call XF_CREATE_COORDINATE_GROUP(xFileId, xCoordId, status)
  if (status < 0) then
    call XF_CLOSE_GROUP(xGridId, error)
    call XF_CLOSE_FILE(xFileId, error)
  error = status
  return
  endif

  iSpcZone = 3601  ! Oregon North

  call XF_SET_HORIZ_DATUM(xCoordId, HORIZ_DATUM_STATE_PLANE_NAD27, error)
  call XF_SET_HORIZ_UNITS(xCoordId, COORD_UNITS_US_FEET, error)

  call XF_SET_VERT_DATUM(xCoordId, VERT_DATUM_LOCAL, error)
  call XF_SET_VERT_UNITS(xCoordId, COORD_UNITS_US_FEET, error)

    ! write additional information
  call XF_SET_SPC_ZONE(xCoordId, iSpcZone, error)

  call XF_CLOSE_GROUP(xCoordId, error)
  xCoordId = 0

  ! release memory
call XF_CLOSE_GROUP(xGridId, error)
call XF_CLOSE_FILE(xFileId, error)
return

END SUBROUTINE TG_WRITE_TEST_GRID_CART_2D

!****************************
!------------------------------------------------------------------------------
! FUNCTION   tgWriteTestGridCurv2D
! PURPOSE    Write a file that contains data for a 2D Curvilinear Grid
! NOTES      A picture of the grid is TestGridCurv2D.gif
!            returns TRUE on success and FALSE on failure
!------------------------------------------------------------------------------
SUBROUTINE TG_WRITE_TEST_GRID_CURV_2D(Filename, Compression, error)
use const_def, only: READWRITE, READONLY, OVERWRITE
  
CHARACTER(LEN=*), INTENT(IN) :: Filename
INTEGER, INTENT(IN) :: Compression
INTEGER, INTENT(OUT) :: error
INTEGER        nDimensions
INTEGER        nCompOrigin, nUDir
INTEGER        nCellsI, nCellsJ
INTEGER        nCells
INTEGER        nCorners
INTEGER        nGridType
REAL(DOUBLE)   xVals(16), yVals(16) !There are 16 corners
INTEGER        i
INTEGER(XID) xFileId, xGridId, xPropId, xDatasetsId, xScalarId
INTEGER(XID) xCoordId
REAL(DOUBLE)   dNullValue(1)
INTEGER        nOrientation
REAL           fDsetCellVals(6) !For cell-centered dataset
REAL           fDsetCornerVals(12) !For corner-centered dataset
REAL(DOUBLE)   tempdouble, dCppLat, dCppLon
INTEGER*1      bDsetCellActive(6)
INTEGER*1      bDsetCornerActive(12)
INTEGER        status
INTEGER        tmpOut1, tmpOut2

  nDimensions = 2
  nCompOrigin = 1
  nUDir = 1;
  nCellsI = 2
  nCellsJ = 3
  nCells = nCellsI*nCellsJ
  nCorners = (nCellsI + 1)*(nCellsJ + 1)
  nGridType = GRID_TYPE_CURVILINEAR
  xFileId = NONE
  xGridId = NONE
  xPropId = NONE
  xDatasetsId = NONE
  xScalarId = NONE
  dNullValue(1) = -999.0
  nOrientation = ORIENTATION_RIGHT_HAND

    ! There is no cell in the top right corner so we have a NullValue for
    ! the top right corner

    ! xValues row by row
  xVals(1) = 0.0
  xVals(2) = 7.5
  xVals(3) = 15.0
  xVals(4) = 2.5
  xVals(5) = 10.0
  xVals(6) = 17.5
  xVals(7) = 3.5
  xVals(8) = 11.0
  xVals(9) = 18.5
  xVals(10) = 0.0
  xVals(11) = 7.5
  xVals(12) = dNullValue(1)

    ! yValues row by row
  yVals(1) = 0.0
  yVals(2) = 0.0
  yVals(3) = 0.0
  yVals(4) = 10.0
  yVals(5) = 10.0
  yVals(6) = 10.0
  yVals(7) = 20.0
  yVals(8) = 20.0
  yVals(9) = 20.0
  yVals(10) = 30.0
  yVals(11) = 30.0
  yVals(12) = dNullValue(1)

    ! cell centered velocity dataset values
  fDsetCellVals(1) = 2.1
  fDsetCellVals(2) = 2.0
  fDsetCellVals(3) = 1.9
  fDsetCellVals(4) = 2.3
  fDsetCellVals(5) = 2.5
  fDsetCellVals(6) = dNullValue(1)

    ! all are active except the last value
  do i = 1, nCells
    bDsetCellActive(i) = 1
  enddo
  bDsetCellActive(nCells) = 0

    ! corner centered elevation dataset values
  fDsetCornerVals(1) = 1.0
  fDsetCornerVals(2) = 0.8 
  fDsetCornerVals(3) = 1.2
  fDsetCornerVals(4) = 1.4
  fDsetCornerVals(5) = 1.8
  fDsetCornerVals(6) = 2.2
  fDsetCornerVals(7) = 1.8
  fDsetCornerVals(8) = 1.4
  fDsetCornerVals(9) = 2.0
  fDsetCornerVals(10) = 1.0
  fDsetCornerVals(11) = 1.8
  fDsetCornerVals(12) = 2.2

    ! all are active except the last value
  do i = 1, nCorners
    bDsetCornerActive(i) = 1
  enddo
  bDsetCornerActive(nCorners) = 0

    ! create the file
  call XF_CREATE_FILE(Filename, OVERWRITE, xFileId, status)
  if (status < 0) then
    error = -1
    return
  endif

    ! create the group to store the grid
  call XF_CREATE_GROUP_FOR_GRID(xFileId, GRID_CURV2D_GROUP_NAME, xGridId, status)
  if (status < 0) then
    call XF_CLOSE_FILE(xFileId, error)
    error = -1
    return
  endif

    ! Write the grid information to the file
  call XF_SET_GRID_TYPE(xGridId, nGridType, tmpOut1)
  call XF_SET_NUMBER_OF_DIMENSIONS(xGridId, nDimensions, tmpOut2)
  if ((tmpOut1 < 0) .OR. (tmpOut2 < 0)) then
    call XF_CLOSE_GROUP(xGridId, error)
    call XF_CLOSE_FILE(xFileId, error)
    error = -1
    return
  endif

    ! set orientation
  call XF_SET_ORIENTATION(xGridId, nOrientation, tmpOut1)  
  if (tmpOut1 < 0 ) then
    call XF_CLOSE_GROUP(xGridId, error)
    call XF_CLOSE_FILE(xFileId, error)
    error = -1
    return
  endif

    ! Set computational origin
  call XF_SET_COMPUTATIONAL_ORIGIN(xGridId, nCompOrigin, tmpOut1)  
  if (tmpOut1 < 0 ) then
    call XF_CLOSE_GROUP(xGridId, error)
    call XF_CLOSE_FILE(xFileId, error)
    error = -1
    return
  endif

    ! Set u direction
  call XF_SET_U_DIRECTION(xGridId, nUDir, tmpOut1)  
  if (tmpOut1 < 0 ) then
    call XF_CLOSE_GROUP(xGridId, error)
    call XF_CLOSE_FILE(xFileId, error)
    error = -1
    return
  endif

    ! Write the grid geometry to the file
    ! Set the number of cells in each direction
  call XF_SET_NUMBER_CELLS_IN_I(xGridId, nCellsI, tmpOut1)
  call XF_SET_NUMBER_CELLS_IN_J(xGridId, nCellsJ, tmpOut2)
  if ((tmpOut1 < 0) .OR. (tmpOut2 < 0)) then
    call XF_CLOSE_GROUP(xGridId, error)
    call XF_CLOSE_FILE(xFileId, error)
    error = -1
    return
  endif
    ! Set a NullValue.  This is used to identify locations in the grid that are
    ! not being used.  In our case no geometry is defined for the top right
    ! corner.
  call XF_CREATE_GRID_PROPERTY_GROUP(xGridId, xPropId, tmpOut1)
  if (xPropId < 0) then
    call XF_CLOSE_GROUP(xGridId, error)
    call XF_CLOSE_FILE(xFileId, error)
    error = -1
    return
  endif

  call XF_WRITE_PROPERTY_DOUBLE(xPropId, PROP_NULL_VALUE, 1, dNullValue, NONE, &
                                tmpOut1)  
  if (tmpOut1 < 0) then
    call XF_CLOSE_GROUP(xPropId, error)
    call XF_CLOSE_GROUP(xGridId, error)
    call XF_CLOSE_FILE(xFileId, error)
    error = -1
    return
  endif
  call XF_CLOSE_GROUP(xPropId, error)

    ! Set the grid plane locations
  call XF_SET_GRID_COORDS_I(xGridId, nCorners, xVals, tmpOut1)
  call XF_SET_GRID_COORDS_J(xGridId, nCorners, yVals, tmpOut2)
  if ((tmpOut1 < 0) .OR. (tmpOut2 < 0)) then
    call XF_CLOSE_GROUP(xGridId, error)
    call XF_CLOSE_FILE(xFileId, error)
    error = -1
    return
  endif

    ! Create the datasets group
  call XF_CREATE_GENERIC_GROUP(xGridId, 'Datasets', xDatasetsId, tmpOut1)  
  if (tmpOut1 < 0) then
    call XF_CLOSE_GROUP(xGridId, error)
    call XF_CLOSE_FILE(xFileId, error)
    error = -1
    return
  endif

    ! Create the cell-centered dataset
  call XF_CREATE_SCALAR_DATASET(xDatasetsId, 'Velocity Mag', 'ft/s', TS_MINUTES, &
                                Compression, xScalarId, tmpOut1)  
  if (tmpOut1 < 0) then
    call XF_CLOSE_GROUP(xDatasetsId, error)
    call XF_CLOSE_GROUP(xGridId, error)
    call XF_CLOSE_FILE(xFileId, error)
    error = -1
    return
  endif

    ! specify that the dataset is cell-centered
  call XF_SCALAR_DATA_LOCATION(xScalarId, GRID_LOC_CENTER, tmpOut1)
  if (tmpOut1 < 0) then
    call XF_CLOSE_GROUP(xScalarId, error)
    call XF_CLOSE_GROUP(xDatasetsId, error)
    call XF_CLOSE_GROUP(xGridId, error)
    call XF_CLOSE_FILE(xFileId, error)
    error = -1
    return
  endif

    ! Write the data
    ! set a temporary variable to pass into the function
  tempdouble = 0.0
  call XF_WRITE_SCALAR_TIMESTEP(xScalarId, tempdouble, nCells, fDsetCellVals, tmpOut1)
  call XF_WRITE_ACTIVITY_TIMESTEP(xScalarId, nCells, bDsetCellActive, tmpOut2)
  if ((tmpOut1 < 0) .OR. (tmpOut1 < 0)) then
    call XF_CLOSE_GROUP(xScalarId, error)
    call XF_CLOSE_GROUP(xDatasetsId, error)
    call XF_CLOSE_GROUP(xGridId, error)
    call XF_CLOSE_FILE(xFileId, error)
    error = -1
    return
  endif

    ! close the cell-centered dataset
  call XF_CLOSE_GROUP(xScalarId, error)

    ! Create the corner-centered dataset
  call XF_CREATE_SCALAR_DATASET(xDatasetsId, 'elevation', 'ft', TS_MINUTES, &
                                Compression, xScalarId, tmpOut1)  
  if (tmpOut1 < 0) then
    call XF_CLOSE_GROUP(xDatasetsId, error)
    call XF_CLOSE_GROUP(xGridId, error)
    call XF_CLOSE_FILE(xFileId, error)
    error = -1
    return
  endif

    ! specify that the dataset is corner-centered
  call XF_SCALAR_DATA_LOCATION(xScalarId, GRID_LOC_CORNER, tmpOut1)  
  if (tmpOut1 < 0) then
    call XF_CLOSE_GROUP(xScalarId, error)
    call XF_CLOSE_GROUP(xDatasetsId, error)
    call XF_CLOSE_GROUP(xGridId, error)
    call XF_CLOSE_FILE(xFileId, error)
    error = -1
    return
  endif

    ! Write the data
    ! set a temporary variable to pass into the function
  tempdouble = 0.0
  call XF_WRITE_SCALAR_TIMESTEP(xScalarId, tempdouble, nCorners, fDsetCornerVals, tmpOut1)
  call XF_WRITE_ACTIVITY_TIMESTEP(xScalarId, nCorners, bDsetCornerActive, tmpOut2)
  if ((tmpOut1 < 0) .OR. (tmpOut2 < 0)) then
    call XF_CLOSE_GROUP(xScalarId, error)
    call XF_CLOSE_GROUP(xDatasetsId, error)
    call XF_CLOSE_GROUP(xGridId, error)
    call XF_CLOSE_FILE(xFileId, error)
    error = -1
    return
  endif

  ! Write Coordinate file - for GridCurv2D, we will set the coordinate system
  !   to be CPP, with CPP Latitude and CPP Longitude settings written 
  !   to the file.
  call XF_CREATE_COORDINATE_GROUP(xFileId, xCoordId, status)
  if (status < 0) then
    call XF_CLOSE_GROUP(xScalarId, error)
    call XF_CLOSE_GROUP(xDatasetsId, error)
    call XF_CLOSE_GROUP(xGridId, error)
    call XF_CLOSE_FILE(xFileId, error)
    error = status
    return
  endif

  dCppLat = 56.0;  ! Made-up value
  dCppLon = 23.0;  ! Made-up value

  call XF_SET_HORIZ_DATUM(xCoordId, HORIZ_DATUM_CPP, error)
  call XF_SET_HORIZ_UNITS(xCoordId, COORD_UNITS_METERS, error)

  call XF_SET_VERT_DATUM(xCoordId, VERT_DATUM_LOCAL, error)
  call XF_SET_VERT_UNITS(xCoordId, COORD_UNITS_METERS, error)

    ! write additional information
  call XF_SET_CPP_LAT(xCoordId, dCppLat, error)
  call XF_SET_CPP_LON(xCoordId, dCppLon, error)

  call XF_CLOSE_GROUP(xCoordId, error)
  xCoordId = 0

  ! release memory
call XF_CLOSE_GROUP(xScalarId, error)
call XF_CLOSE_GROUP(xDatasetsId, error)
call XF_CLOSE_GROUP(xGridId, error)
call XF_CLOSE_FILE(xFileId, error)
return

END SUBROUTINE TG_WRITE_TEST_GRID_CURV_2D

!****************************
!------------------------------------------------------------------------------
! FUNCTION   tgWriteTestGridCart3D
! PURPOSE    Write a file that contains data for a 2D Cartesian Grid
! NOTES      A picture of the grid is in the file (TestGridCart2D.gif)
!            returns TRUE on success and FALSE on failure
!------------------------------------------------------------------------------
SUBROUTINE TG_WRITE_TEST_GRID_CART_3D(Filename, Compression, error)
use const_def, only: READWRITE, READONLY, OVERWRITE
  
CHARACTER(LEN=*), INTENT(IN) :: Filename
INTEGER, INTENT(IN) :: Compression
INTEGER, INTENT(OUT) :: error
INTEGER        nDimensions
INTEGER        nCompOrigin, nUDir
INTEGER        nCellsI, nCellsJ, nCellsK
INTEGER        nGridType
REAL(DOUBLE)   dOriginX, dOriginY, dOriginZ
INTEGER        nOrientation
REAL(DOUBLE)   dBearing, dDip, dRoll
REAL(DOUBLE)   PlanesI(5), PlanesJ(5), PlanesK(3)
INTEGER        i, j, status, iSpcZone
INTEGER(XID) xFileId, xGridId, xPropId, xCoordId
INTEGER        nCells
INTEGER        Active(75)
INTEGER        tmpOut1, tmpOut2, tmpOut3

nDimensions = 3
nCompOrigin = 8
nUDir = -2
nCellsI = 5
nCellsJ = 5
nCellsK = 3
nGridType = GRID_TYPE_CARTESIAN
dOriginX = 10.0
dOriginY = 10.0
dOriginZ = 0.0
nOrientation = ORIENTATION_RIGHT_HAND
dBearing = 45.0
dDip = 0.0
dRoll = 0.0
xFileId = NONE
xGridId = NONE
xPropId = NONE
nCells = nCellsI * nCellsJ * nCellsK

  ! Fill in the grid plane data with a constant size of 30
do i = 1, nCellsI
  PlanesI(i) = i*30.0
enddo
do j = 1, nCellsJ
  PlanesJ(j) = j*30.0
enddo
do j = 1,  nCellsK
  PlanesK(j) = j*30.0
enddo

  ! fill in the activity array
  ! default array to active
do i = 1, nCells
  Active(i) = 1
enddo

  ! two cells are inactive (identified by array index)
  ! i = 0, j = 0, k = 0  and i = 4, j = 4, k = 0
Active(1) = 0
Active(4*nCellsJ*nCellsK+4*nCellsK+1) = 0

  ! create the file
call XF_CREATE_FILE(Filename, OVERWRITE, xFileId, tmpOut1)  
if (tmpOut1 < 0) then
  error = -1
  return
endif

  ! create the group to store the grid
call XF_CREATE_GROUP_FOR_GRID(xFileId, GRID_CART3D_GROUP_NAME, xGridId, tmpOut1)
if (tmpOut1 < 0) then
  call XF_CLOSE_FILE(xFileId, error)
  error = -1
  return
endif

  ! Write the grid information to the file
call XF_SET_GRID_TYPE(xGridId, nGridType, tmpOut1)
call XF_SET_NUMBER_OF_DIMENSIONS(xGridId, nDimensions, tmpOut2)
if ((tmpOut1 < 0) .OR. (tmpOut2 < 0)) then
  call XF_CLOSE_GROUP(xGridId, error)
  call XF_CLOSE_FILE(xFileId, error)
  error = -1
  return
endif

  ! set origin and orientation
call XF_SET_ORIGIN(xGridId, dOriginX, dOriginY, dOriginZ, tmpOut1)
call XF_SET_ORIENTATION(xGridId, nOrientation, tmpOut2)
if ((tmpOut1 < 0) .OR. (tmpOut2 < 0)) then
  call XF_CLOSE_GROUP(xGridId, error)
  call XF_CLOSE_FILE(xFileId, error)
  error = -1
  return
endif
 
  ! Set bearing, dip and roll
call XF_SET_BEARING(xGridId, dBearing, tmpOut1)
call XF_SET_DIP(xGridId, dDip, tmpOut2)
call XF_SET_ROLL(xGridId, dRoll, tmpOut3)
if ((tmpOut1 < 0) .OR. (tmpOut2 < 0) .OR. (tmpOut3 < 0)) then
  call XF_CLOSE_GROUP(xGridId, error)
  call XF_CLOSE_FILE(xFileId, error)
  error = -1
  return
endif

 ! Set computational origin
call XF_SET_COMPUTATIONAL_ORIGIN(xGridId, nCompOrigin, tmpOut1)
if (tmpOut1 < 0) then
  call XF_CLOSE_GROUP(xGridId, error)
  call XF_CLOSE_FILE(xFileId, error)
  error = -1
  return
endif

  ! Set u direction
call XF_SET_U_DIRECTION(xGridId, nUDir, tmpOut1)
if (tmpOut1 < 0) then
  call XF_CLOSE_GROUP(xGridId, error)
  call XF_CLOSE_FILE(xFileId, error)
  error = -1
  return
endif

  ! Write the grid geometry to the file
  ! Set the number of cells in each direction
call XF_SET_NUMBER_CELLS_IN_I(xGridId, nCellsI, tmpOut1)
call XF_SET_NUMBER_CELLS_IN_J(xGridId, nCellsJ, tmpOut2)
call XF_SET_NUMBER_CELLS_IN_K(xGridId, nCellsK, tmpOut3)  
if ((tmpOut1 < 0) .OR. (tmpOut2 < 0) .OR. (tmpOut3 < 0)) then
  call XF_CLOSE_GROUP(xGridId, error)
  call XF_CLOSE_FILE(xFileId, error)
  error = -1
  return
endif

  ! Set the grid plane locations
call XF_SET_GRID_COORDS_I(xGridId, nCellsI, PlanesI, tmpOut1)
call XF_SET_GRID_COORDS_J(xGridId, nCellsJ, PlanesJ, tmpOut2)
call XF_SET_GRID_COORDS_K(xGridId, nCellsK, PlanesK, tmpOut3)  

if ((tmpOut1 < 0) .OR. (tmpOut2 < 0) .OR. (tmpOut3 < 0)) then
  call XF_CLOSE_GROUP(xGridId, error)
  call XF_CLOSE_FILE(xFileId, error)
  error = -1
  return
endif

  ! Write the activity array
call XF_CREATE_GRID_CELL_PROP_GRP(xGridId, xPropId, tmpOut1)  
if (xPropId < 0) then
  call XF_CLOSE_GROUP(xGridId, error)
  call XF_CLOSE_FILE(xFileId, error)
  error = -1
  return
endif

call XF_WRITE_PROPERTY_INT(xPropId, PROP_ACTIVITY, nCells, Active, &
                           Compression, tmpOut1)

if (tmpOut1 < 0) then
  call XF_CLOSE_GROUP(xPropId, error)
  call XF_CLOSE_GROUP(xGridId, error)
  call XF_CLOSE_FILE(xFileId, error)
  error = -1
  return
endif
  
call XF_CLOSE_GROUP(xPropId, error)

  ! Write Coordinate file - for GridCart3D, we will set the coordinate system
  !   to be State Plane NAD27.
  call XF_CREATE_COORDINATE_GROUP(xFileId, xCoordId, status)
  if (status < 0) then
    call XF_CLOSE_GROUP(xGridId, error)
    call XF_CLOSE_FILE(xFileId, error)
    error = status
  endif

  iSpcZone = 3601; ! Oregon North

  call XF_SET_HORIZ_DATUM(xCoordId, HORIZ_DATUM_STATE_PLANE_NAD27, error)
  call XF_SET_HORIZ_UNITS(xCoordId, COORD_UNITS_US_FEET, error)

  call XF_SET_VERT_DATUM(xCoordId, VERT_DATUM_LOCAL, error)
  call XF_SET_VERT_UNITS(xCoordId, COORD_UNITS_US_FEET, error)

    ! write additional information
  call XF_SET_SPC_ZONE(xCoordId, iSpcZone, error)

  call XF_CLOSE_GROUP(xCoordId, error)
  xCoordId = 0

  ! release memory
call XF_CLOSE_GROUP(xGridId, error)
call XF_CLOSE_FILE(xFileId, error)
return

END SUBROUTINE TG_WRITE_TEST_GRID_CART_3D

END MODULE TestGrid