MODULE TestTimestep

USE Xmdf

CHARACTER(LEN=*), PARAMETER :: DATASETS_LOCATION = 'Datasets'
CHARACTER(LEN=*), PARAMETER :: SCALAR_A_LOCATION = 'Scalars/ScalarA'
CHARACTER(LEN=*), PARAMETER :: SCALAR_B_LOCATION = 'Scalars/ScalarB'
CHARACTER(LEN=*), PARAMETER :: VECTOR2D_A_LOCATION = 'Vectors/Vector2D_A'
CHARACTER(LEN=*), PARAMETER :: VECTOR2D_B_LOCATION = 'Vectors/Vector2D_B'

CONTAINS

! --------------------------------------------------------------------------
! FUNCTION ttiTestNumTimes
! PURPOSE  Change the NumTimes to truncate timesteps
! NOTES    
! --------------------------------------------------------------------------
RECURSIVE SUBROUTINE TTI_Test_Num_Times (a_DatasetId, a_Itimestep, error)
INTEGER(XID), INTENT(IN) :: a_DatasetId
INTEGER, INTENT(IN)        :: a_Itimestep
INTEGER, INTENT(OUT)       :: error
INTEGER                    NumTimes
INTEGER                    Itimestep

error = 1

! Fortran loops are 1 to 3 but C is 0 to 2, etc
Itimestep = a_Itimestep - 1

  ! truncate just written timestep and test error conditions
if (1 == Itimestep .OR. 3 == Itimestep .OR. 5 == Itimestep) then
      ! Test setting NumTimes after end of dataset
      call XF_SET_DATASET_NUM_TIMES( a_DatasetId, Itimestep + 2, error )
      if (error >= 0) then
         WRITE(*,*) 'ERROR1: XF_SET_DATASET_NUM_TIMES must return ERROR.'
      endif

      if (1 == Itimestep) then
         Itimestep = 1;
      endif
      if (3 == Itimestep) then
         Itimestep = 2;
      endif
      if (5 == Itimestep) then
         Itimestep = 3;
      endif

      ! Write actual NumTimes
      call XF_SET_DATASET_NUM_TIMES( a_DatasetId, Itimestep, error )
      if (error < 0) then
         WRITE(*,*) 'ERROR2: xfSetDatasetNumTimes must NOT return error.'
      endif

      ! Test setting NumTimes after end step.
      call XF_SET_DATASET_NUM_TIMES( a_DatasetId, Itimestep + 1, error )
      if (error >= 0) then
         WRITE(*,*) 'ERROR3: xfSetDatasetNumTimes must return ERROR.'
      endif

      ! Test reading NumTimes
      call XF_GET_DATASET_NUM_TIMES( a_DatasetId, NumTimes, error )
      if (error < 0) then
         WRITE(*,*) 'ERROR4: xfSetDatasetNumTimes must NOT return error.'
      endif
      if (NumTimes .NE. Itimestep) then
         WRITE(*,*) 'ERROR5: xfGetDatasetNumTimes must return CORRECT NumTimes.'
      endif
endif

return

END SUBROUTINE
!ttiTestNumTimes
! --------------------------------------------------------------------------
! FUNCTION ttReadDatasets
! PURPOSE  Read a dataset group from an XMDF file and output information to
!          to a text file
! NOTES    
! --------------------------------------------------------------------------
RECURSIVE SUBROUTINE TT_READ_DATASETS (a_xGroupId, a_FileUnit, error)
INTEGER(XID), INTENT(IN) :: a_xGroupId
INTEGER, INTENT(IN)        :: a_FileUnit
INTEGER, INTENT(OUT)       :: error
INTEGER                   nPaths, nMaxPathLength, j
CHARACTER, ALLOCATABLE, DIMENSION(:) :: Paths
CHARACTER(LEN=500)       IndividualPath
INTEGER                   nStatus, i
INTEGER(XID)            xScalarId, xVectorId, xMultiId
INTEGER                   nMultiDatasets

xScalarId = NONE
xVectorId = NONE

nMultiDatasets = 0
nPaths = 0
nMaxPathLength = 0

  ! Look for scalar datasets
call XF_GET_SCALAR_DATASETS_INFO(a_xGroupId, nPaths, nMaxPathLength, nStatus)
if (nStatus >= 0 .AND. nPaths > 0) then
  allocate(Paths(nPaths*nMaxPathLength))
  call XF_GET_SCALAR_DATASET_PATHS(a_xGroupId, nPaths, nMaxPathLength, Paths, &
                                                                         error)
endif
if (nStatus < 0) then
  error = -1
  return
endif

  ! Output number and paths to scalar datasets
WRITE(a_FileUnit,*) 'Number of Scalars ', nPaths
do i=2, nPaths
  IndividualPath = ''
  do j=1, nMaxPathLength-1
    IndividualPath(j:j) = Paths((i-1)*nMaxPathLength+j)
  enddo
  WRITE(a_FileUnit,*) 'Reading scalar: ', IndividualPath(1:nMaxPathLength-1)
  call XF_OPEN_GROUP(a_xGroupId, IndividualPath(1:nMaxPathLength-1), &
                                                           xScalarId, nStatus)
  if (nStatus < 0) then
    error = -1
  return
  endif

  call TTI_READ_SCALAR(xScalarId, a_FileUnit, nStatus)
  call XF_CLOSE_GROUP(xScalarId, error)
  if (nStatus < 0) then
    WRITE(*,*) 'Error reading scalar dataset.'
    error = -1
  return
  endif
enddo

if (allocated(Paths)) deallocate(Paths)
  ! Look for vector datasets
call XF_GET_VECTOR_DATASETS_INFO(a_xGroupId, nPaths, nMaxPathLength, nStatus)
if (nStatus >= 0 .AND. nPaths > 0) then
  allocate(Paths(nPaths*nMaxPathLength))
  call XF_GET_VECTOR_DATASET_PATHS(a_xGroupId, nPaths, nMaxPathLength, Paths, error)
endif
if (nStatus < 0) then
  error = -1
  return
endif

  ! Output number and paths to scalar datasets
WRITE(a_FileUnit,*) 'Number of Vectors ', nPaths
do i=2, nPaths
  do j=1, nMaxPathLength-1
    IndividualPath(j:j) = Paths((i-1)*nMaxPathLength+j)
  enddo
  WRITE(a_FileUnit,*) 'Reading Vector: ', &
                      IndividualPath(1:nMaxPathLength-1)
  call XF_OPEN_GROUP(a_xGroupId, IndividualPath(1:nMaxPathLength-1), &
                                                          xVectorId, nStatus)
  if (nStatus < 0) then
    error = -1
  return
  endif
  call TTI_READ_VECTOR(xVectorId, a_FileUnit, nStatus)
  call XF_CLOSE_GROUP(xVectorId, error)
  if (nStatus < 0) then
    WRITE(*,*) 'Error reading vector dataset.'
    error = -1
  return
  endif
enddo

if (allocated(Paths)) deallocate(Paths)

! find multidataset folders
call XF_GET_GRP_PTHS_SZ_MLT_DSETS(a_xGroupId, nMultiDatasets, &
                                                      nMaxPathLength, nStatus)
if (nStatus >= 0 .AND. nMultiDatasets > 0) then
  allocate(Paths(nMultiDatasets*nMaxPathLength))
  call XF_GET_ALL_GRP_PATHS_MLT_DSETS(a_xGroupId, nMultiDatasets, &
                                                 nMaxPathLength, Paths, error)
  if (nStatus < 0) then
    error = -1
    return
  endif

  ! Output number and paths to multidatasets
  WRITE(a_FileUnit,*) 'Number of Multidatasets ', nMultiDatasets
  do i=2, nMultiDatasets
    IndividualPath = ''
    do j=1, nMaxPathLength-1
      IndividualPath(j:j) = Paths((i-1)*nMaxPathLength+j)
    enddo
    WRITE(a_FileUnit,*) 'Reading multidataset: ', &
                                             IndividualPath(1:nMaxPathLength-1)
    call XF_OPEN_GROUP(a_xGroupId, IndividualPath(1:nMaxPathLength-1), &
                                                           xMultiId, nStatus)
    if (nStatus < 0) then
      error = -1
    return
    endif

    call TT_READ_DATASETS(xMultiId, a_FileUnit, nStatus)
    call XF_CLOSE_GROUP(xMultiId, error)
    if (nStatus < 0) then
      WRITE(*,*) 'Error reading multidatasets.'
      error = -1
    return
    endif
  enddo
endif
if (allocated(Paths)) deallocate(Paths)

error = 1
return

END SUBROUTINE
!ttReadDatasets
! --------------------------------------------------------------------------
! FUNCTION ttReadActivityScalarAIndex
! PURPOSE  Read all timestep values for a particular index
! NOTES    
! --------------------------------------------------------------------------
SUBROUTINE TT_READ_ACTIVITY_SCALAR_A_INDEX(a_Filename, a_Index, error)
use const_def, only: READWRITE, READONLY, OVERWRITE
  
CHARACTER(LEN=*), INTENT(IN) :: a_Filename
INTEGER, INTENT(IN)   :: a_Index
INTEGER, INTENT(OUT)  :: error
INTEGER                  status
INTEGER(XID)           xFileId, xDsetsId, xScalarAId
INTEGER                  nTimesteps, i
INTEGER, ALLOCATABLE  :: bActive(:)

xFileId = NONE
xDsetsId = NONE
xScalarAId = NONE

  ! open the file
call XF_OPEN_FILE(a_Filename, READONLY, xFileId, status)
if (status < 0) then
  error = -1
  return
endif

  ! open the dataset group
call XF_OPEN_GROUP(xFileId, DATASETS_LOCATION, xDsetsId, status)
if (status >= 0) then
  call XF_OPEN_GROUP(xDsetsId, SCALAR_A_LOCATION, xScalarAId, status)
endif
if (status < 0) then
  error = status
  return
endif

  ! Find out the number of timesteps in the file
CALL XF_GET_DATASET_NUM_TIMES(xScalarAId, nTimesteps, status)
if (status < 0) then
  error = status
  return
endif

if (nTimesteps < 1) then
  error = -1
  return
endif

  ! Read the values for the index
allocate(bActive(nTimesteps))
call XF_READ_ACTIVE_VALS_AT_INDEX(xScalarAId, a_Index, 1, nTimesteps, &
                                       bActive, status)
  ! output the data
WRITE(*,*) ''
WRITE(*,*) 'Reading activity for scalar A slice at index: ', a_Index
do i=1, nTimesteps
  WRITE(*,*) bActive(i), ' '
enddo

deallocate(bActive)

error = status
return

END SUBROUTINE
! ttReadActivityScalarAIndex

! --------------------------------------------------------------------------
! FUNCTION ttReadScalarAIndex
! PURPOSE  Read all timestep values for a particular index
! NOTES    
! --------------------------------------------------------------------------
SUBROUTINE TT_READ_SCALAR_A_INDEX (a_Filename, a_Index, error)
use const_def, only: READWRITE, READONLY, OVERWRITE
  
CHARACTER(LEN=*), INTENT(IN) :: a_Filename
INTEGER, INTENT(IN)   :: a_Index
INTEGER, INTENT(OUT)  :: error
INTEGER              status
INTEGER(XID)       xFileId, xDsetsId, xScalarAId
INTEGER              nTimesteps, i
REAL, ALLOCATABLE :: fValues(:)

xFileId = NONE
xDsetsId = NONE
xScalarAId = NONE

  ! open the file
call XF_OPEN_FILE(a_Filename, READONLY, xFileId, status)
if (status < 0) then
  error = -1
  return
endif

  ! open the dataset group
call XF_OPEN_GROUP(xFileId, DATASETS_LOCATION, xDsetsId, status)
if (status >= 0) then
  call XF_OPEN_GROUP(xDsetsId, SCALAR_A_LOCATION, xScalarAId, status)
endif
if (status < 0) then
  error = status
  return
endif

  ! Find out the number of timesteps in the file
call XF_GET_DATASET_NUM_TIMES(xScalarAId, nTimesteps, status)
if (status < 0) then
  error = status
  return
endif

if (nTimesteps < 1) then
  error = -1
  return
endif

  ! Read the values for the index
allocate (fValues(nTimesteps))
call XF_READ_SCALAR_VALUES_AT_INDEX(xScalarAId, a_Index, 1, nTimesteps, &
                                     fValues, status)

  ! output the data
WRITE(*,*) ''
WRITE(*,*) 'Reading scalar A slice at index: ', a_Index
do i=1, nTimesteps
  WRITE(*,*) fValues(i), ' '
enddo

deallocate(fValues)

error = status
return

END SUBROUTINE
! ttReadScalarAtIndex

! --------------------------------------------------------------------------
! FUNCTION ttWriteScalarA
! PURPOSE  Write scalar Dataset to an HDF5 File
! NOTES    This tests dynamic data sets, and activity
!          This dataset is dynamic concentrations (mg/L) with output times
!          in minutes.
!          Dataset is for a mesh and so nActive is the number of elements
!          which is not the same as the nValues which would be number of nodes
!          reads/writes a reference time in julian days
! --------------------------------------------------------------------------
SUBROUTINE TT_WRITE_SCALAR_A (a_Filename, a_Compression, error)
  use const_def, only: READWRITE, READONLY, OVERWRITE
  
  CHARACTER(LEN=*), INTENT(IN) :: a_Filename
  INTEGER, INTENT(IN)          :: a_Compression
  INTEGER, INTENT(OUT)         :: error
  INTEGER(XID)      xFileId, xDsetsId, xScalarAId, xCoordId
  INTEGER      nValues, nTimes, nActive
  REAL(DOUBLE) dTime, dJulianReftime
  INTEGER      iTimestep, iActive, iHpgnZone
  REAL         fValues(10) ! nValues
  INTEGER*1    bActivity(10) ! activity
  INTEGER      i, status

  ! initialize the data
  nValues = 10
  nTimes = 3
  nActive = 8
  dTime = 0.0

  ! 5th item in data set is always inactive, others active
  do iActive = 1, nActive
    bActivity(iActive) = 1
  enddo 
  bActivity(6) = 0


  ! create the file
  call XF_CREATE_FILE(a_Filename, OVERWRITE, xFileId, error)
  if (error .LT. 0) then
      ! close the file
    call XF_CLOSE_FILE(xFileId, error)
    return
  endif

  ! create the group where we will put all the datasets 
  call XF_CREATE_GENERIC_GROUP(xFileId, DATASETS_LOCATION, xDsetsId, status)
  if (status < 0) then
    call XF_CLOSE_FILE(xFileId, error)
    error = -1
  return
  endif

  ! Create the scalar A dataset group
  call XF_CREATE_SCALAR_DATASET(xDsetsId, SCALAR_A_LOCATION, 'mg/L', &
              TS_HOURS, a_Compression, xScalarAId, status)
  if (status .LT. 0) then
      ! close the dataset
    call XF_CLOSE_GROUP(xScalarAId, error)
    call XF_CLOSE_GROUP(xDsetsId, error)
    call XF_CLOSE_FILE(xFileId, error)
    error = status
    return 
  endif

  ! Add in a reftime.  This is a julian day for:
  ! noon July 1, 2003
  dJulianReftime = 2452822.0;
  call XF_WRITE_REFTIME(xScalarAId, dJulianReftime, status)
  if (status < 0) then
    call XF_CLOSE_GROUP(xScalarAId, error)
    call XF_CLOSE_GROUP(xDsetsId, error)
    call XF_CLOSE_FILE(xFileId, error)
  endif

  ! Loop through timesteps adding them to the file
  do iTimestep = 1, nTimes
    ! We will have an 0.5 hour timestep
    dTime = iTimestep * 0.5

    fValues(1) = dTime
    do i = 2, nValues
      fValues(i) = fValues(i-1)*2.5
    end do

    ! write the dataset array values
    call XF_WRITE_SCALAR_TIMESTEP(xScalarAId, dTime, nValues, fValues, error)
    if (error .GE. 0) then
      ! write activity array
      call XF_WRITE_ACTIVITY_TIMESTEP(xScalarAId, nActive, bActivity, error)
    end if 

    call TTI_Test_Num_Times(xScalarAid, iTimestep, error)
  enddo

  ! Write Coordinate file - for ScalarA, we will set the coordinate system
  !   to be Geographic HPGN, with HPGN settings written to the file.
  call XF_CREATE_COORDINATE_GROUP(xFileId, xCoordId, status)
  if (status < 0) then
    call XF_CLOSE_GROUP(xScalarAId, error)
    call XF_CLOSE_GROUP(xDsetsId, error)
    call XF_CLOSE_FILE(xFileId, error)
    error = -1
  return
  endif

    ! set HPGN Zone for test
  iHpgnZone = 29   ! Utah
    ! Write Coordinate Information to file
  call XF_SET_HORIZ_DATUM(xCoordId, HORIZ_DATUM_GEOGRAPHIC_HPGN, error)
  call XF_SET_HORIZ_UNITS(xCoordId, COORD_UNITS_METERS, error)
  call XF_SET_VERT_DATUM(xCoordId, VERT_DATUM_LOCAL, error)
  call XF_SET_VERT_UNITS(xCoordId, COORD_UNITS_METERS, error)

    ! write additional information
  call XF_SET_HPGN_AREA(xCoordId, iHpgnZone, error)

  call XF_CLOSE_GROUP(xCoordId, error)
  xCoordId = 0;

  ! close the dataset
  call XF_CLOSE_GROUP(xScalarAId, error)
  call XF_CLOSE_GROUP(xDsetsId, error)
  call XF_CLOSE_FILE(xFileId, error)

  return
END SUBROUTINE
! ttWriteScalarA

! --------------------------------------------------------------------------
! FUNCTION TT_WRITE_SCALAR_B
! PURPOSE  Write scalar Dataset to an HDF5 File
! NOTES    This tests dynamic data sets, and activity
!          This dataset is dynamic concentrations (mg/L) with output times
!          in minutes.
!          Dataset is for a mesh and so nActive is the number of elements
!          which is not the same as the nValues which would be number of nodes
!          reads/writes a reference time in julian days
! --------------------------------------------------------------------------
SUBROUTINE TT_WRITE_SCALAR_B (a_Filename, a_Compression, a_Overwrite, error)
  use const_def, only: READWRITE, READONLY, OVERWRITE
  
  CHARACTER(LEN=*), INTENT(IN) :: a_Filename
  INTEGER, INTENT(IN)          :: a_Compression
  LOGICAL, INTENT(IN)          :: a_Overwrite
  INTEGER, INTENT(OUT)         :: error
  INTEGER(XID)      xFileId, xDsetsId, xScalarBId, xCoordId
  INTEGER      nValues, nTimes, nActive
  REAL(DOUBLE) dTime, dJulianReftime
  INTEGER      iTimestep, iActive
  REAL         fValues(10) ! nValues
  INTEGER*1    bActivity(10) ! activity
  INTEGER      i, status

  ! initialize the data
  nValues = 10
  nTimes = 3
  nActive = 8
  dTime = 0.0
  i = 0

  ! 5th item in data set is always inactive, others active
  do iActive = 1, nActive
    bActivity(iActive) = 1
  enddo 
  bActivity(6) = 0

  if (a_Overwrite) then
      ! open the already-existing file
    call XF_OPEN_FILE(a_Filename, READWRITE, xFileId, status)
    if (status < 0) then
      error = -1
      return
    endif
      ! open the group where we have all the datasets
    call XF_OPEN_GROUP(xFileId, DATASETS_LOCATION, xDsetsId, status)
    if (status < 0) then
      call XF_CLOSE_FILE(xFileId, error)
      error = -1
      return
    endif
  else
      ! create the file
    call XF_CREATE_FILE(a_Filename, OVERWRITE, xFileId, error)
    if (error .LT. 0) then
        ! close the file
      call XF_CLOSE_FILE(xFileId, error)
      return
    endif

      ! create the group where we will put all the datasets 
    call XF_CREATE_GENERIC_GROUP(xFileId, DATASETS_LOCATION, xDsetsId, status)
    if (status < 0) then
      call XF_CLOSE_FILE(xFileId, error)
      error = -1
    return
    endif
  endif

  ! Create/Overwrite the scalar B dataset group
  call XF_CREATE_SCALAR_DATASET(xDsetsId, SCALAR_B_LOCATION, 'mg/L', &
              TS_HOURS, a_Compression, xScalarBId, status)
  if (status < 0) then
      ! close the dataset
    call XF_CLOSE_GROUP(xScalarBId, error)
    call XF_CLOSE_GROUP(xDsetsId, error)
    call XF_CLOSE_FILE(xFileId, error)
    error = status
    return 
  endif

  ! Add in a reftime.  This is a julian day for:
  ! noon July 1, 2003
  dJulianReftime = 2452822.0;
  call XF_WRITE_REFTIME(xScalarBId, dJulianReftime, status)
  if (status < 0) then
    call XF_CLOSE_GROUP(xScalarBId, error)
    call XF_CLOSE_GROUP(xDsetsId, error)
    call XF_CLOSE_FILE(xFileId, error)
  endif

  if (.NOT. a_Overwrite) then
      ! Loop through timesteps adding them to the file
    do iTimestep = 1, nTimes
        ! We will have an 0.5 hour timestep
      dTime = iTimestep * 0.5

      fValues(1) = dTime
      do i = 2, nValues
        fValues(i) = fValues(i-1)*2.5
      end do

        ! write the dataset array values
      call XF_WRITE_SCALAR_TIMESTEP(xScalarBId, dTime, nValues, fValues, error)
      if (error .GE. 0) then
          ! write activity array
        call XF_WRITE_ACTIVITY_TIMESTEP(xScalarBId, nActive, bActivity, error)
      end if
      if (error < 0) then
          call XF_CLOSE_GROUP(xScalarBId, error)
          call XF_CLOSE_GROUP(xDsetsId, error)
          call XF_CLOSE_FILE(xFileId, error)
      endif
    enddo
  else
      ! Loop through timesteps adding them to the file
    do iTimestep = 1, nTimes
        ! We will have an 1.5 hour timestep
      dTime = iTimestep * 1.5

      fValues(1) = dTime
      do i = 2, nValues
        fValues(i) = fValues(i-1)*1.5
      end do

        ! write the dataset array values
      call XF_WRITE_SCALAR_TIMESTEP(xScalarBId, dTime, nValues, fValues, error)
      if (error .GE. 0) then
          ! write activity array
        call XF_WRITE_ACTIVITY_TIMESTEP(xScalarBId, nActive, bActivity, error)
      end if
      if (error < 0) then
          call XF_CLOSE_GROUP(xScalarBId, error)
          call XF_CLOSE_GROUP(xDsetsId, error)
          call XF_CLOSE_FILE(xFileId, error)
      endif
    enddo
  endif

  if (.NOT. a_Overwrite) then
    ! Write Coordinate file - for ScalarB, we will set the coordinate system
    !   to be UTM, with UTM Zone settings written to the file.
    call XF_CREATE_COORDINATE_GROUP(xFileId, xCoordId, status)
    if (status < 0) then
      call XF_CLOSE_GROUP(xScalarBId, error)
      call XF_CLOSE_GROUP(xDsetsId, error)
      call XF_CLOSE_FILE(xFileId, error)
    error = -1
    return
    endif

     ! Write Coord Info to file
    call XF_SET_HORIZ_DATUM(xCoordId, HORIZ_DATUM_UTM, error)
    call XF_SET_HORIZ_UNITS(xCoordId, COORD_UNITS_METERS, error)

    call XF_SET_VERT_DATUM(xCoordId, VERT_DATUM_LOCAL, error)
    call XF_SET_VERT_UNITS(xCoordId, COORD_UNITS_METERS, error)

      ! write additional information - we'll use the max value for this test
    call XF_SET_UTM_ZONE(xCoordId, UTM_ZONE_MAX, error)

    call XF_CLOSE_GROUP(xCoordId, error)
    xCoordId = 0
  endif

  ! close the dataset
  call XF_CLOSE_GROUP(xScalarBId, error)
  call XF_CLOSE_GROUP(xDsetsId, error)
  call XF_CLOSE_FILE(xFileId, error)

  error = 1
  return
END SUBROUTINE
! ttWriteScalarB
!------------------------------------------------------------------------------
!  FUNCTION TT_WRITE_COORDS_TO_MULTI
!  PURPOSE  Write coordinate system to a multidataset file
!  NOTES
!------------------------------------------------------------------------------
SUBROUTINE TT_WRITE_COORDS_TO_MULTI (a_xFileId, error)
INTEGER(XID), INTENT(IN) :: a_xFileId
INTEGER, INTENT(OUT)       :: error
INTEGER(XID)    xCoordId
INTEGER           status

  ! Write Coordinate file - for Multidatasets, we will set the coordinate system
  !   to be UTM, with UTM Zone settings written to the file.
  call XF_CREATE_COORDINATE_GROUP(a_xFileId, xCoordId, status)
  if (status < 0) then
    error = status
  return
  endif

    ! Write Coord Info to file
  call XF_SET_HORIZ_DATUM(xCoordId, HORIZ_DATUM_UTM, error)
  call XF_SET_HORIZ_UNITS(xCoordId, COORD_UNITS_METERS, error)

  call XF_SET_VERT_DATUM(xCoordId, VERT_DATUM_LOCAL, error)
  call XF_SET_VERT_UNITS(xCoordId, COORD_UNITS_METERS, error)

    ! write additional information - we'll use the max value for this test
  call XF_SET_UTM_ZONE(xCoordId, UTM_ZONE_MAX, error)

  call XF_CLOSE_GROUP(xCoordId, error)
  xCoordId = 0

  return
END SUBROUTINE

! --------------------------------------------------------------------------
! FUNCTION ttWriteScalarAToMulti
! PURPOSE  Write scalar Dataset to an HDF5 File
! NOTES    This tests dynamic data sets, and activity
!          This dataset is dynamic concentrations (mg/L) with output times
!          in minutes.
!          Dataset is for a mesh and so nActive is the number of elements
!          which is not the same as the nValues which would be number of nodes
!          reads/writes a reference time in julian days
! --------------------------------------------------------------------------
SUBROUTINE TT_WRITE_SCALAR_A_TO_MULTI (a_GroupID, status)
 ! CHARACTER(LEN=*), INTENT(IN) :: a_Filename
 ! INTEGER, INTENT(IN)          :: a_Compression
 ! INTEGER, INTENT(OUT)         :: error
  INTEGER(XID)      xFileId, xDsetsId, xScalarAId
  INTEGER(XID)      a_GroupID
  INTEGER      nValues, nTimes, nActive
  REAL(DOUBLE) dTime, dJulianReftime
  INTEGER      iTimestep, iActive
  REAL         fValues(10) ! nValues
  INTEGER*1    bActivity(10) ! activity
  INTEGER      i, status

  ! initialize the data
  nValues = 10
  nTimes  = 3
  nActive = 8
  dTime   = 0.0

  ! 5th item in data set is always inactive, others active
  do iActive = 1, nActive
    bActivity(iActive) = 1
  enddo 
  bActivity(6) = 0

  ! Create the scalar A dataset group
  call XF_CREATE_SCALAR_DATASET(a_GroupID, SCALAR_A_LOCATION, 'mg/L', &
              TS_HOURS, NONE, xScalarAId, status)
  if (status .LT. 0) then
      ! close the dataset
    call XF_CLOSE_GROUP(xScalarAId, status)
    call XF_CLOSE_GROUP(xDsetsId, status)
    call XF_CLOSE_FILE(xFileId, status)
    return 
  endif

  ! Add in a reftime.  This is a julian day for:
  ! noon July 1, 2003
  dJulianReftime = 2452822.0;
  call XF_WRITE_REFTIME(xScalarAId, dJulianReftime, status)
  if (status < 0) then
    call XF_CLOSE_GROUP(xScalarAId, status)
    call XF_CLOSE_GROUP(xDsetsId, status)
    call XF_CLOSE_FILE(xFileId, status)
  endif

  ! Loop through timesteps adding them to the file
  do iTimestep = 1, nTimes
    ! We will have an 0.5 hour timestep
    dTime = iTimestep * 0.5

    fValues(1) = dTime
    do i = 2, nValues
      fValues(i) = fValues(i-1)*2.5
    end do

    ! write the dataset array values
    call XF_WRITE_SCALAR_TIMESTEP(xScalarAId, dTime, nValues, fValues, status)
    if (status .GE. 0) then
      ! write activity array
      call XF_WRITE_ACTIVITY_TIMESTEP(xScalarAId, nActive, bActivity, status)
    end if 

    call TTI_Test_Num_Times(xScalarAId, iTimestep, status)
  enddo

  ! close the dataset
  call XF_CLOSE_GROUP(xScalarAId, status)
  !call XF_CLOSE_GROUP(a_GroupID, status)
  !call XF_CLOSE_FILE(a_FileID, status)

  return
END SUBROUTINE
! ttWriteScalarAToMulti
! --------------------------------------------------------------------------
! FUNCTION ttReadVector2DAIndex
! PURPOSE  Read all timestep values for a particular index
! NOTES    
! --------------------------------------------------------------------------
SUBROUTINE TT_READ_VECTOR2D_A_INDEX (a_Filename, a_Index, error)
use const_def, only: READWRITE, READONLY, OVERWRITE
  
CHARACTER(LEN=*), INTENT(IN) :: a_Filename
INTEGER, INTENT(IN)   :: a_Index
INTEGER, INTENT(OUT)  :: error
INTEGER           status
INTEGER(XID)    xFileId, xDsetsId, xVector2DA
INTEGER           nTimesteps, i
REAL, ALLOCATABLE     :: fValues(:)

xFileId = NONE
xDsetsId = NONE
xVector2DA = NONE

  ! open the file
call XF_OPEN_FILE(a_Filename, READONLY, xFileId, status)
if (status < 0) then
  error = -1
  return
endif

  ! open the dataset group
call XF_OPEN_GROUP(xFileId, DATASETS_LOCATION, xDsetsId, status)
if (status >= 0) then
  call XF_OPEN_GROUP(xDsetsId, VECTOR2D_A_LOCATION, xVector2DA, status)
endif
if (status < 0) then
  error = status
  return
endif

  ! Find out the number of timesteps in the file
call XF_GET_DATASET_NUM_TIMES(xVector2DA, nTimesteps, status)
if (status < 0) then
  error = status
  return
endif

if (nTimesteps < 1) then
  error = -1
  return
endif

  ! Read the values for the index
allocate(fValues(nTimesteps*2))
call XF_READ_VECTOR_VALUES_AT_INDEX(xVector2DA, a_Index, 1, nTimesteps, 2, &
                                     fValues, status)

  ! output the data
WRITE(*,*) ''
WRITE(*,*) 'Reading vector 2D A slice at index: ', a_Index
do i=1, nTimesteps
  WRITE(*,*) fValues(i*2-1), ' ', fValues(i*2)
enddo
WRITE(*,*) ''

deallocate(fValues)

error = status
return

END SUBROUTINE
!ttReadVector2DAIndex

! --------------------------------------------------------------------------
! FUNCTION ttWriteVector2D_A
! PURPOSE  Write scalar Dataset to an HDF5 File
! NOTES    This tests dynamic data sets, and activity
!          This dataset is dynamic concentrations (mg/L) with output times
!          in minutes.
!          Dataset is for a mesh and so nActive is the number of elements
!          which is not the same as the nValues which would be number of nodes
!          reads/writes a reference time in julian days
! --------------------------------------------------------------------------
SUBROUTINE TT_WRITE_VECTOR2D_A (a_Filename, a_Compression, error)
  use const_def, only: READWRITE, READONLY, OVERWRITE
  
  CHARACTER(LEN=*), INTENT(IN) :: a_Filename
  INTEGER, INTENT(IN)          :: a_Compression
  INTEGER, INTENT(OUT)         :: error
  INTEGER(XID)      xFileId, xDsetsId, xVector2D_A, xCoordId
  INTEGER      nValues, nTimes, nComponents, nActive
  REAL(DOUBLE) dTime
  INTEGER      iTimestep, iActive
  REAL, DIMENSION(2, 100) :: fValues ! nComponents, nValues
  INTEGER*1    bActivity(100) ! activity
  INTEGER      i, j, status
  INTEGER      iHpgnZone

  ! initialize the data
  nComponents = 2
  nValues = 100
  nTimes = 6
  nActive = 75
  dTime = 0.0

  ! 5th item in data set is always inactive, others active
  bActivity(1) = 0
  do iActive = 2, nActive
    if (mod(iActive-1, 3) == 0) then
      bActivity(iActive) = 0
    else
    bActivity(iActive) = 1
  endif
  enddo

  ! create the file
  call XF_CREATE_FILE(a_Filename, OVERWRITE, xFileId, error)
  if (error .LT. 0) then
    ! close the dataset
    call XF_CLOSE_FILE(xFileId, error)
    return
  endif

  ! create the group where we will put all the datasets 
  call XF_CREATE_GENERIC_GROUP(xFileId, DATASETS_LOCATION, xDsetsId, status)
  if (status < 0) then
    call XF_CLOSE_FILE(xFileId, error)
    error = -1
  return
  endif

  ! Create the vector dataset group
  call XF_CREATE_VECTOR_DATASET(xDsetsId, VECTOR2D_A_LOCATION, 'ft/s', &
              TS_SECONDS, a_Compression, xVector2D_A, status)
  if (status .LT. 0) then
      ! close the dataset
    call XF_CLOSE_GROUP(xVector2D_A, error)
    call XF_CLOSE_GROUP(xDsetsId, error)
    call XF_CLOSE_FILE(xFileId, error)
    error = status
    return 
  endif

  ! Loop through timesteps adding them to the file
  do iTimestep = 1, nTimes
    ! We will have an 0.5 hour timestep
    dTime = iTimestep * 0.5

    do i = 1, nValues
      do j = 1, nComponents
        fValues(j,i) = ((i-1)*nComponents + j)*dTime
      end do
    end do

    ! write the dataset array values
    call XF_WRITE_VECTOR_TIMESTEP(xVector2D_A, dTime, nValues, nComponents, &
                                  fValues, error)
    if (error .GE. 0) then
      ! write activity array
      call XF_WRITE_ACTIVITY_TIMESTEP(xVector2D_A, nActive, bActivity, error)
    end if 

    call TTI_Test_Num_Times(xVector2d_A, iTimestep, error)
  enddo

  ! Write Coordinate file - for Vector2D_A, we will set the coordinate system
  !   to be Geographic HPGN, with HPGN settings written to the file.
  call XF_CREATE_COORDINATE_GROUP(xFileId, xCoordId, status)
  if (status < 0) then
    call XF_CLOSE_GROUP(xVector2D_A, error)
    call XF_CLOSE_GROUP(xDsetsId, error)
    call XF_CLOSE_FILE(xFileId, error)
  error = -1
  return
  endif

    ! set HPGN info for test
  iHpgnZone = 29   ! Utah

  call XF_SET_HORIZ_DATUM(xCoordId, HORIZ_DATUM_GEOGRAPHIC_HPGN, error)
  call XF_SET_HORIZ_UNITS(xCoordId, COORD_UNITS_METERS, error)
  call XF_SET_VERT_DATUM(xCoordId, VERT_DATUM_LOCAL, error)
  call XF_SET_VERT_UNITS(xCoordId, COORD_UNITS_METERS, error)

    ! write additional information
  call XF_SET_HPGN_AREA(xCoordId, iHpgnZone, error)

  call XF_CLOSE_GROUP(xCoordId, error)
  xCoordId = 0

  ! close the dataset
  call XF_CLOSE_GROUP(xVector2D_A, error)
  call XF_CLOSE_GROUP(xDsetsId, error)
  call XF_CLOSE_FILE(xFileId, error)

  return
END SUBROUTINE
! ttWriteVector2D_A

! --------------------------------------------------------------------------
! FUNCTION TT_WRITE_VECTOR2D_B
! PURPOSE  Write scalar Dataset to an HDF5 File
! NOTES    This tests dynamic data sets, and activity
!          This dataset is dynamic concentrations (mg/L) with output times
!          in minutes.
!          Dataset is for a mesh and so nActive is the number of elements
!          which is not the same as the nValues which would be number of nodes
!          reads/writes a reference time in julian days
! --------------------------------------------------------------------------
SUBROUTINE TT_WRITE_VECTOR2D_B (a_Filename, a_Compression, a_Overwrite, error)
  use const_def, only: READWRITE, READONLY, OVERWRITE
  
  CHARACTER(LEN=*), INTENT(IN) :: a_Filename
  INTEGER, INTENT(IN)          :: a_Compression
  LOGICAL, INTENT(IN)          :: a_Overwrite
  INTEGER, INTENT(OUT)         :: error
  INTEGER(XID)      xFileId, xDsetsId, xVector2D_B, xCoordId
  INTEGER      nValues, nTimes, nComponents, nActive
  REAL(DOUBLE) dTime
  INTEGER      iTimestep, iActive
  REAL, DIMENSION(2, 100) :: fValues
  INTEGER*1    bActivity(100)
  INTEGER      i, j, status

    ! initialize the data
  nComponents = 2
  nValues = 100
  nTimes = 6
  nActive = 75
  dTime = 0.0

    ! 5th item in data set is always inactive, others active
  bActivity(1) = 0
  do iActive = 2, nActive
    if (mod(iActive-1, 3) == 0) then
      bActivity(iActive) = 0
    else
      bActivity(iActive) = 1
    endif
  enddo

  if (a_Overwrite) then
      ! open the already-existing file
    call XF_OPEN_FILE(a_Filename, READWRITE, xFileId, status)
    if (status < 0) then
      error = -1
      return
    endif
      ! open the group where we have all the datasets
    call XF_OPEN_GROUP(xFileId, DATASETS_LOCATION, xDsetsId, status)
    if (status < 0) then
      call XF_CLOSE_FILE(xFileId, error)
      error = -1
      return
    endif
  else
      ! create the file
    call XF_CREATE_FILE(a_Filename, OVERWRITE, xFileId, error)
    if (error .LT. 0) then
        ! close the dataset
      call XF_CLOSE_FILE(xFileId, error)
      return
    endif

      ! create the group where we will put all the datasets 
    call XF_CREATE_GENERIC_GROUP(xFileId, DATASETS_LOCATION, xDsetsId, status)
    if (status < 0) then
      call XF_CLOSE_FILE(xFileId, error)
      error = -1
      return
    endif
  endif

    ! Create/Overwrite the vector dataset group
  call XF_CREATE_VECTOR_DATASET(xDsetsId, VECTOR2D_B_LOCATION, 'ft/s', &
                                TS_SECONDS, a_Compression, xVector2D_B, status)
  if (status .LT. 0) then
      ! close the dataset
    call XF_CLOSE_GROUP(xVector2D_B, error)
    call XF_CLOSE_GROUP(xDsetsId, error)
    call XF_CLOSE_FILE(xFileId, error)
    error = status
    return 
  endif

  if (.NOT. a_Overwrite) then
      ! Loop through timesteps adding them to the file
    do iTimestep = 1, nTimes
        ! We will have an 0.5 hour timestep
      dTime = iTimestep * 0.5
      do i = 1, nValues
        do j = 1, nComponents
          fValues(j,i) = ((i-1)*nComponents + j)*dTime
        end do
      end do
        ! write the dataset array values
      call XF_WRITE_VECTOR_TIMESTEP(xVector2D_B, dTime, nValues, nComponents, &
                                    fValues, error)
      if (error .GE. 0) then
          ! write activity array
        call XF_WRITE_ACTIVITY_TIMESTEP(xVector2D_B, nActive, bActivity, error)
      end if
      if (error < 0) then
        call XF_CLOSE_GROUP(xVector2D_B, error)
        call XF_CLOSE_GROUP(xDsetsId, error)
        call XF_CLOSE_FILE(xFileId, error)
      endif
    enddo
  else
      ! Loop through timesteps adding them to the file
    do iTimestep = 1, nTimes
        ! We will have an 1.5 hour timestep
      dTime = iTimestep * 1.5
      do i = 1, nValues
        do j = 1, nComponents
          fValues(j,i) = ((i-1)*nComponents + j)*dTime
        end do
      end do
        ! write the dataset array values
      call XF_WRITE_VECTOR_TIMESTEP(xVector2D_B, dTime, nValues, nComponents, &
                                    fValues, error)
      if (error .GE. 0) then
          ! write activity array
        call XF_WRITE_ACTIVITY_TIMESTEP(xVector2D_B, nActive, bActivity, error)
      end if
      if (error < 0) then
        call XF_CLOSE_GROUP(xVector2D_B, error)
        call XF_CLOSE_GROUP(xDsetsId, error)
        call XF_CLOSE_FILE(xFileId, error)
      endif
    enddo
  endif

  if (.NOT. a_Overwrite) then
    ! Write Coordinate file - for ScalarB, we will set the coordinate system
    !   to be UTM, with UTM Zone settings written to the file.
    call XF_CREATE_COORDINATE_GROUP(xFileId, xCoordId, status)
    if (status < 0) then
      call XF_CLOSE_GROUP(xVector2D_B, error)
      call XF_CLOSE_GROUP(xDsetsId, error)
      call XF_CLOSE_FILE(xFileId, error)
    error = -1
    return
    endif

      ! write the coordinate data to the file
    call XF_SET_HORIZ_DATUM(xCoordId, HORIZ_DATUM_UTM, error)
    call XF_SET_HORIZ_UNITS(xCoordId, COORD_UNITS_METERS, error)
    call XF_SET_VERT_DATUM(xCoordId, VERT_DATUM_LOCAL, error)
    call XF_SET_VERT_UNITS(xCoordId, COORD_UNITS_METERS, error)

      ! write additional information - we'll use the max UTM zone for the test
    call XF_SET_UTM_ZONE(xCoordId, UTM_ZONE_MAX, error)

    call XF_CLOSE_GROUP(xCoordId, error)
    xCoordId = 0
  endif

  ! close the dataset
  call XF_CLOSE_GROUP(xVector2D_B, error)
  call XF_CLOSE_GROUP(xDsetsId, error)
  call XF_CLOSE_FILE(xFileId, error)

  return
END SUBROUTINE
! ttWriteVector2D_B

! --------------------------------------------------------------------------
! FUNCTION ttWriteVector2D_AToMulti
! PURPOSE  Write scalar Dataset to an HDF5 File
! NOTES    This tests dynamic data sets, and activity
!          This dataset is dynamic concentrations (mg/L) with output times
!          in minutes.
!          Dataset is for a mesh and so nActive is the number of elements
!          which is not the same as the nValues which would be number of nodes
!          reads/writes a reference time in julian days
! --------------------------------------------------------------------------
SUBROUTINE TT_WRITE_VECTOR2D_A_TO_MULTI (a_FileID, a_GroupID, status)
  INTEGER(XID)      xVector2D_A
  INTEGER(XID)      a_FileID, a_GroupID
  INTEGER      nValues, nTimes, nComponents, nActive
  REAL(DOUBLE) dTime
  INTEGER      iTimestep, iActive
  REAL, DIMENSION(2, 100) :: fValues ! nComponents, nValues
  INTEGER*1    bActivity(100) ! activity
  INTEGER      i, j, status

  ! initialize the data
  nComponents = 2
  nValues = 100
  nTimes = 6
  nActive = 75
  dTime = 0.0

  ! 5th item in data set is always inactive, others active
  bActivity(1) = 0
  do iActive = 2, nActive
    if (mod(iActive-1, 3) == 0) then
      bActivity(iActive) = 0
    else
    bActivity(iActive) = 1
  endif
  enddo

  ! Create the vector dataset group
  call XF_CREATE_VECTOR_DATASET(a_GroupID, VECTOR2D_A_LOCATION, 'ft/s', &
              TS_SECONDS, NONE, xVector2D_A, status)
  if (status .LT. 0) then
      ! close the dataset
    call XF_CLOSE_GROUP(xVector2D_A, status)
    call XF_CLOSE_GROUP(a_GroupID, status)
    call XF_CLOSE_FILE(a_FileID, status)
    return 
  endif

  ! Loop through timesteps adding them to the file
  do iTimestep = 1, nTimes
    ! We will have an 0.5 hour timestep
    dTime = iTimestep * 0.5

    do i = 1, nValues
      do j = 1, nComponents
        fValues(j,i) = ((i-1)*nComponents + j)*dTime
      end do
    end do

    ! write the dataset array values
    call XF_WRITE_VECTOR_TIMESTEP(xVector2D_A, dTime, nValues, nComponents, &
                                  fValues, status)
    if (status .GE. 0) then
      ! write activity array
      call XF_WRITE_ACTIVITY_TIMESTEP(xVector2D_A, nActive, bActivity, status)
    end if 

    call TTI_Test_Num_Times(xVector2D_A, iTimestep, status)
  enddo

  ! close the dataset
  call XF_CLOSE_GROUP(xVector2D_A, status)
  return
END SUBROUTINE
! ttWriteVector2D_AToMulti
! --------------------------------------------------------------------------
! FUNCTION ttiReadScalar
! PURPOSE  Read a scalar from an XMDF file and output information to
!          to a text file
! NOTES    
! --------------------------------------------------------------------------
SUBROUTINE TTI_READ_SCALAR (a_xScalarId, FileUnit, error)
  INTEGER(XID), INTENT(IN) ::  a_xScalarId
  INTEGER, INTENT(IN) ::         FileUnit
  INTEGER, INTENT(OUT) :: error
  INTEGER             nTimes, nValues, nActive
  LOGICAL*2             bUseReftime
  INTEGER             iTime
  CHARACTER(LEN=100)   TimeUnits
  REAL(DOUBLE), ALLOCATABLE :: Times(:)
  REAL, ALLOCATABLE         :: Values(:), Minimums(:), Maximums(:)
  INTEGER, ALLOCATABLE      :: Active(:)
  REAL(DOUBLE)                 Reftime
nTimes = NONE
nValues = NONE
nActive = None

  ! read the time units
  call XF_GET_DATASET_TIME_UNITS(a_xScalarId, TimeUnits, error)
  if (error < 0) return

  WRITE(FileUnit,*) 'Time units: ', TimeUnits(1:LEN_TRIM(TimeUnits))

  ! see if we are using a reftime
  call XF_USE_REFTIME (a_xScalarId, bUseReftime, error)
  if (error < 0) then
    return
  endif
  if (bUseReftime) then
    call XF_READ_REFTIME (a_xScalarId, Reftime, error)
    if (error < 0) then
      return
  endif
    WRITE(FileUnit,*) 'Reftime: ', Reftime
  endif

  ! read in the number of values and number of active values
  call XF_GET_DATASET_NUMVALS(a_xScalarId, nValues, error)
  if (error .GE. 0) then
    call XF_GET_DATASET_NUMACTIVE(a_xScalarId, nActive, error)
  endif
  if (error .LT. 0) return 

  if (nValues <= 0) then
    WRITE(FileUnit, *) 'No data to read in.'
    error = -1
    return 
  endif

  ! read in the number of times
  call XF_GET_DATASET_NUM_TIMES(a_xScalarId, nTimes, error)
  if (error < 0) then
    return 
  endif

  ! Read in the individual time values
  allocate(Times(nTimes))

  call XF_GET_DATASET_TIMES(a_xScalarId, nTimes, Times, error)
  if (error < 0) return 

  ! Read in the minimum and maximum values
  allocate(Minimums(nTimes))
  allocate(Maximums(nTimes))

  call XF_GET_DATASET_MINS(a_xScalarId, nTimes, Minimums, error)
  if (error >= 0) then
    call XF_GET_DATASET_MAXS(a_xScalarId, nTimes, Maximums, error)
  endif
  if (error < 0) then
    deallocate(Times)
    deallocate(Minimums)
    deallocate(Maximums)
    return
  endif

  allocate(Values(nValues))
  if (nActive .GT. 0) then
    allocate(Active(nActive))
  endif

  WRITE(FileUnit,*) 'Number Timesteps: ', nTimes
  WRITE(FileUnit,*) 'Number Values: ', nValues
  WRITE(FileUnit,*) 'Number Active: ', nActive
  WRITE(FileUnit,*) ''

  ! loop through the timesteps, read the values and active values and write
  ! them to the text file
  do iTime = 1, nTimes
    call XF_READ_SCALAR_VALUES_TIMESTEP(a_xScalarId, iTime, nValues, Values, error)
    if (error >= 0 .AND. nActive > 0) then
      call XF_READ_ACTIVITY_TIMESTEP(a_xScalarId, iTime, nActive, Active, error)
    endif

    ! Write the time, min, max, values and active values to the text output
    ! file.
    WRITE(FileUnit,*) 'Timestep at  ', Times(iTime)
    WRITE(FileUnit,*) 'Min: ', Minimums(iTime)
    WRITE(FileUnit,*) 'Max: ', Maximums(iTime)

    WRITE(FileUnit,*) 'Values:'
    WRITE(FileUnit,*) Values(1:nValues)
    WRITE(FileUnit,*) ''

    WRITE(FileUnit,*) 'Activity:'
    WRITE(FileUnit,*) Active(1:nActive)
    WRITE(FileUnit,*) ''
  end do

  if (allocated(Times)) then
    deallocate(Times)
  endif
  
  if (allocated(Minimums)) then
    deallocate(Minimums)
  endif

  if (allocated(Maximums)) then
    deallocate(Maximums)
  endif

  if (allocated(Values)) then
    deallocate(Values)
  endif

  if (allocated(Active)) then
    deallocate(Active)
  endif

  return
END SUBROUTINE
! ttiReadScalar

! --------------------------------------------------------------------------
! FUNCTION TTI_READ_VECTOR
! PURPOSE  Read a vector from an XMDF file and output information to
!          to a text file
! NOTES    
! --------------------------------------------------------------------------
SUBROUTINE TTI_READ_VECTOR (a_xVectorId, FileUnit, error)
  INTEGER(XID), INTENT(IN) ::  a_xVectorId
  INTEGER, INTENT(IN) ::         FileUnit
  INTEGER, INTENT(OUT) :: error
  INTEGER             nTimes, nValues, nActive, nComponents
  INTEGER             iTime, i
  LOGICAL*2            bUseReftime
  CHARACTER(LEN=100)   TimeUnits
  REAL(DOUBLE), ALLOCATABLE :: Times(:)
  REAL, ALLOCATABLE, DIMENSION (:, :) :: Values
  REAL, ALLOCATABLE         :: Minimums(:), Maximums(:)
  INTEGER, ALLOCATABLE      :: Active(:)
  REAL(DOUBLE)                 Reftime

nTimes = NONE
nValues = NONE
nActive = NONE
nComponents = NONE

  ! read the time units
  call XF_GET_DATASET_TIME_UNITS(a_xVectorId, TimeUnits, error)
  if (error < 0) return

  WRITE(FileUnit,*) 'Time units: ', TimeUnits(1:LEN_TRIM(TimeUnits))

  ! see if we are using a reftime
  call XF_USE_REFTIME (a_xVectorId, bUseReftime, error)
  if (error < 0) then
    return
  endif
  if (bUseReftime) then
    call XF_READ_REFTIME (a_xVectorId, Reftime, error)
    if (error < 0) then
      return
  endif
    WRITE(FileUnit,*) 'Reftime: ', Reftime
  endif

  ! read in the number of values and number of active values
  call XF_GET_DATASET_NUMVALS(a_xVectorId, nValues, error)
  if (error .GE. 0) then
    call XF_GET_DATASET_NUMCOMPONENTS(a_xVectorId, nComponents, error)
    if (error .GE. 0) then
      call XF_GET_DATASET_NUMACTIVE(a_xVectorId, nActive, error)
    endif
  endif
  if (error .LT. 0) return 

  if (nValues <= 0) then
    WRITE(FileUnit, *) 'No data to read in.'
    error = -1
    return 
  endif

  ! read in the number of times
  call XF_GET_DATASET_NUM_TIMES(a_xVectorId, nTimes, error)
  if (error < 0) then
    return 
  endif

  ! Read in the individual time values
  allocate(Times(nTimes))

  call XF_GET_DATASET_TIMES(a_xVectorId, nTimes, Times, error)
  if (error < 0) return 

  ! Read in the minimum and maximum values
  allocate(Minimums(nTimes))
  allocate(Maximums(nTimes))

  call XF_GET_DATASET_MINS(a_xVectorId, nTimes, Minimums, error)
  if (error >= 0) then
    call XF_GET_DATASET_MAXS(a_xVectorId, nTimes, Maximums, error)
  endif
  if (error < 0) then
    deallocate(Times)
    deallocate(Minimums)
    deallocate(Maximums)
    return
  endif

  allocate(Values(nComponents, nValues))
  if (nActive .GT. 0) then
    allocate(Active(nActive))
  endif

  WRITE(FileUnit,*) 'Number Timesteps: ', nTimes
  WRITE(FileUnit,*) 'Number Values: ', nValues
  WRITE(FileUnit,*) 'Number Components: ', nComponents
  WRITE(FileUnit,*) 'Number Active: ', nActive

  ! loop through the timesteps, read the values and active values and write
  ! them to the text file
  do iTime = 1, nTimes
    call XF_READ_VECTOR_VALUES_TIMESTEP(a_xVectorId, iTime, nValues, &
                                        nComponents, Values, error)
    if (error >= 0 .AND. nActive > 0) then
      call XF_READ_ACTIVITY_TIMESTEP(a_xVectorId, iTime, nActive, Active, error)
    endif

    ! Write the time, min, max, values and active values to the text output
    ! file.
  WRITE(FileUnit,*) ''
    WRITE(FileUnit,*) 'Timestep at  ', Times(iTime)
    WRITE(FileUnit,*) 'Min: ', Minimums(iTime)
    WRITE(FileUnit,*) 'Max: ', Maximums(iTime)

    WRITE(FileUnit,*) 'Values:'
    do i=1, nValues
      WRITE(FileUnit,*) Values(1:nComponents,i:i)
    enddo
    WRITE(FileUnit,*) ''

    WRITE(FileUnit,*) 'Activity:'
    WRITE(FileUnit,*) Active(1:nActive)
    WRITE(FileUnit,*) ''    
  WRITE(FileUnit,*) ''    

  end do

  if (allocated(Times)) then
    deallocate(Times)
  endif
  
  if (allocated(Minimums)) then
    deallocate(Minimums)
  endif

  if (allocated(Maximums)) then
    deallocate(Maximums)
  endif

  if (allocated(Values)) then
    deallocate(Values)
  endif

  if (allocated(Active)) then
    deallocate(Active)
  endif

  return
END SUBROUTINE
! ttiReadVector

END MODULE TestTimestep
