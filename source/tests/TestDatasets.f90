MODULE TestDatasets

   USE Xmdf

   CHARACTER(LEN=*), PARAMETER :: DATASETS_LOCATION = 'Datasets'
   CHARACTER(LEN=*), PARAMETER :: SCALAR_A_LOCATION = 'Scalars/ScalarA'
   CHARACTER(LEN=*), PARAMETER :: SCALAR_A_LOCATION_FULL = 'Datasets/Scalars/ScalarA'
   CHARACTER(LEN=*), PARAMETER :: SCALAR_B_LOCATION = 'Scalars/ScalarB'
   CHARACTER(LEN=*), PARAMETER :: VECTOR2D_A_LOCATION = 'Vectors/Vector2D_A'
   CHARACTER(LEN=*), PARAMETER :: VECTOR2D_B_LOCATION = 'Vectors/Vector2D_B'

CONTAINS
! ---------------------------------------------------------------------------
! FUNCTION  tdEditScalarAValues
! PURPOSE
! NOTES
! ---------------------------------------------------------------------------
   SUBROUTINE TD_EDIT_SCALAR_A_VALUES(a_Filename, a_Compression, error)
      use const_def, only: READWRITE, READONLY, OVERWRITE
      CHARACTER(LEN=*), INTENT(IN) :: a_Filename
      INTEGER, INTENT(IN) :: a_Compression
      INTEGER, INTENT(OUT) :: error
      INTEGER(XID) xFileId, xScalarId
      INTEGER, PARAMETER :: editNumValues = 3
      INTEGER editTimestep
      INTEGER, DIMENSION(editNumValues) :: indices
      REAL, DIMENSION(editNumValues) :: new_values

      CALL TD_WRITE_SCALAR_A(a_Filename, a_Compression, error)
      if (error < 0) then
         return
      end if

      ! open the file and edit the values
      CALL XF_OPEN_FILE(a_Filename, READWRITE, xFileId, error)
      if (error < 0) then
         return
      end if

      CALL XF_OPEN_GROUP(xFileId, SCALAR_A_LOCATION_FULL, xScalarId, error); 
      if (error < 0) then
         CALL XF_CLOSE_FILE(xFileId, error)
         return
      end if

      ! Edit values in timestep 1, make index 1 = 4, index 5 = 40,
      ! and index 10 = 400
      editTimestep = 1
      indices(1) = 1; 
      indices(2) = 5; 
      indices(3) = 10; 
      new_values(1) = 4.0; 
      new_values(2) = 40.0; 
      new_values(3) = 400.0; 
      CALL XF_CHANGE_SCALAR_VALUES_TIMESTEP_FLOAT(xScalarId, editTimestep, editNumValues, &
                                                  indices, new_values, error)
      if (error < 0) then
         CALL XF_CLOSE_GROUP(xScalarId, error)
         CALL XF_CLOSE_FILE(xFileId, error)
         return
      end if

!  RDJ - Technically we shouldn't have to close and reopen a file that we are
!        editing but it seems to fix a crash that the tests are having.
!  CALL XF_CLOSE_GROUP(xScalarId, error)
!  CALL XF_CLOSE_FILE(xFileId, error)
!    ! open the file and edit the values
!  CALL XF_OPEN_FILE(a_Filename, .FALSE., xFileId, error)
!  if (error < 0) then
!    return
!  endif
!
!  CALL XF_OPEN_GROUP(xFileId, SCALAR_A_LOCATION_FULL, xScalarId, error);
!  if (error < 0) then
!    CALL XF_CLOSE_FILE(xFileId, error)
!    return
!  endif

      ! Edit values in timestep 2, make index 2 = 6, index 3 = 60, and
      ! index 9 = 600
      editTimestep = 2
      indices(1) = 2
      indices(2) = 3
      indices(3) = 9
      new_values(1) = -6.0
      new_values(2) = 60.0
      new_values(3) = 6000.0

      CALL XF_CHANGE_SCALAR_VALUES_TIMESTEP_FLOAT(xScalarId, editTimestep, editNumValues, &
                                                  indices, new_values, error)
      if (error < 0) then
         CALL XF_CLOSE_GROUP(xScalarId, error)
         CALL XF_CLOSE_FILE(xFileId, error)
         return
      end if

      CALL XF_CLOSE_GROUP(xScalarId, error)
      CALL XF_CLOSE_FILE(xFileId, error)

      return
   END SUBROUTINE
! tdEditScalarAValues

! ---------------------------------------------------------------------------
! FUNCTION  tdReadScalarAIndices
! PURPOSE
! NOTES
! ---------------------------------------------------------------------------
   SUBROUTINE TD_READ_SCALAR_A_INDICES(a_Filename, a_nIndices, a_indices, error)
      use const_def, only: READWRITE, READONLY, OVERWRITE

      CHARACTER(LEN=*), INTENT(IN) :: a_Filename
      INTEGER, INTENT(IN) :: a_nIndices
      INTEGER, DIMENSION(*), INTENT(IN) :: a_indices
      INTEGER, INTENT(OUT) :: error
      INTEGER(XID) xFileId, xDsetsId, xScalarAId
      INTEGER nTimesteps
      REAL, ALLOCATABLE, DIMENSION(:) ::  fValues
      INTEGER nValues
      INTEGER id, i, j

      ! open the file
      CALL XF_OPEN_FILE(a_Filename, READONLY, xFileId, error)
      if (error < 0) then
         return
      end if

      ! open the dataset group
      CALL XF_OPEN_GROUP(xFileId, DATASETS_LOCATION, xDsetsId, error)
      if (error >= 0) then
         CALL XF_OPEN_GROUP(xDsetsId, SCALAR_A_LOCATION, xScalarAId, error)
      end if
      if (error < 0) then
         return
      end if

      ! Find out the number of timesteps in the file
      CALL XF_GET_DATASET_NUM_TIMES(xScalarAId, nTimesteps, error)
      if (error < 0) then
         return
      end if
      if (nTimesteps < 1) then
         error = -1
         return
      end if

      ! Read the values for the index
      nValues = nTimesteps*a_nIndices
      allocate (fValues(nValues))
      CALL XF_READ_SCALAR_VALUES_AT_INDICES_FLOAT(xScalarAId, a_nIndices, a_indices, 1, &
                                                  nTimesteps, fValues, error)
      if (error < 0) then
         return
      end if

      ! output the data
      WRITE (*, *) ''
      WRITE (*, *) 'Reading scalar A indices'
      id = 1; 
      do i = 1, nTimesteps
         WRITE (*, *) 'Timestep: ', i
         do j = 1, a_nIndices
            WRITE (*, *) 'index:', a_indices(j), ' value: ', fValues(id)
            id = id + 1
         end do
      end do
      WRITE (*, *) ''

      deallocate (fValues)

      return

   END SUBROUTINE
! TD_READ_SCALARA_INDICES

! --------------------------------------------------------------------------
! FUNCTION tdReadDatasets
! PURPOSE  Read a dataset group from an XMDF file and output information to
!          to a text file
! NOTES
! --------------------------------------------------------------------------
   RECURSIVE SUBROUTINE TD_READ_DATASETS(a_xGroupId, a_FileUnit, error)
      INTEGER(XID), INTENT(IN) :: a_xGroupId
      INTEGER, INTENT(IN)        :: a_FileUnit
      INTEGER, INTENT(OUT)       :: error
      INTEGER nPaths, nMaxPathLength, j
      CHARACTER, ALLOCATABLE, DIMENSION(:) :: Paths
      CHARACTER(LEN=500) IndividualPath
      INTEGER nStatus, i
      INTEGER(XID) xScalarId, xVectorId, xMultiId
      INTEGER nMultiDatasets

      xScalarId = NONE
      xVectorId = NONE

      nMultiDatasets = 0
      nPaths = 0
      nMaxPathLength = 0

      ! Look for scalar datasets
      call XF_GET_SCALAR_DATASETS_INFO(a_xGroupId, nPaths, nMaxPathLength, nStatus)
      if (nStatus >= 0 .AND. nPaths > 0) then
         allocate (Paths(nPaths*nMaxPathLength))
         call XF_GET_SCALAR_DATASET_PATHS(a_xGroupId, nPaths, nMaxPathLength, Paths, &
                                          error)
      end if
      if (nStatus < 0) then
         error = -1
         return
      end if

      ! Output number and paths to scalar datasets
      WRITE (a_FileUnit, *) 'Number of Scalars ', nPaths
      do i = 2, nPaths
         IndividualPath = ''
         do j = 1, nMaxPathLength - 1
            IndividualPath(j:j) = Paths((i - 1)*nMaxPathLength + j)
         end do
         WRITE (a_FileUnit, *) 'Reading scalar: ', IndividualPath(1:nMaxPathLength - 1)
         call XF_OPEN_GROUP(a_xGroupId, IndividualPath(1:nMaxPathLength - 1), &
                            xScalarId, nStatus)
         if (nStatus < 0) then
            error = -1
            return
         end if

         call TDI_READ_SCALAR(xScalarId, a_FileUnit, nStatus)
         call XF_CLOSE_GROUP(xScalarId, error)
         if (nStatus < 0) then
            WRITE (*, *) 'Error reading scalar dataset.'
            error = -1
            return
         end if
      end do

      if (allocated(Paths)) deallocate (Paths)
      ! Look for vector datasets
      call XF_GET_VECTOR_DATASETS_INFO(a_xGroupId, nPaths, nMaxPathLength, nStatus)
      if (nStatus >= 0 .AND. nPaths > 0) then
         allocate (Paths(nPaths*nMaxPathLength))
         call XF_GET_VECTOR_DATASET_PATHS(a_xGroupId, nPaths, nMaxPathLength, Paths, error)
      end if
      if (nStatus < 0) then
         error = -1
         return
      end if

      ! Output number and paths to scalar datasets
      WRITE (a_FileUnit, *) 'Number of Vectors ', nPaths
      do i = 2, nPaths
         do j = 1, nMaxPathLength - 1
            IndividualPath(j:j) = Paths((i - 1)*nMaxPathLength + j)
         end do
         WRITE (a_FileUnit, *) 'Reading Vector: ', &
            IndividualPath(1:nMaxPathLength - 1)
         call XF_OPEN_GROUP(a_xGroupId, IndividualPath(1:nMaxPathLength - 1), &
                            xVectorId, nStatus)
         if (nStatus < 0) then
            error = -1
            return
         end if
         call TDI_READ_VECTOR(xVectorId, a_FileUnit, nStatus)
         call XF_CLOSE_GROUP(xVectorId, error)
         if (nStatus < 0) then
            WRITE (*, *) 'Error reading vector dataset.'
            error = -1
            return
         end if
      end do

      if (allocated(Paths)) deallocate (Paths)

! find multidataset folders
      call XF_GET_GRP_PTHS_SZ_MLT_DSETS(a_xGroupId, nMultiDatasets, &
                                        nMaxPathLength, nStatus)
      if (nStatus >= 0 .AND. nMultiDatasets > 0) then
         allocate (Paths(nMultiDatasets*nMaxPathLength))
         call XF_GET_ALL_GRP_PATHS_MLT_DSETS(a_xGroupId, nMultiDatasets, &
                                             nMaxPathLength, Paths, error)
         if (nStatus < 0) then
            error = -1
            return
         end if

         ! Output number and paths to multidatasets
         WRITE (a_FileUnit, *) 'Number of Multidatasets ', nMultiDatasets
         do i = 2, nMultiDatasets
            IndividualPath = ''
            do j = 1, nMaxPathLength - 1
               IndividualPath(j:j) = Paths((i - 1)*nMaxPathLength + j)
            end do
            WRITE (a_FileUnit, *) 'Reading multidataset: ', &
               IndividualPath(1:nMaxPathLength - 1)
            call XF_OPEN_GROUP(a_xGroupId, IndividualPath(1:nMaxPathLength - 1), &
                               xMultiId, nStatus)
            if (nStatus < 0) then
               error = -1
               return
            end if

            call TD_READ_DATASETS(xMultiId, a_FileUnit, nStatus)
            call XF_CLOSE_GROUP(xMultiId, error)
            if (nStatus < 0) then
               WRITE (*, *) 'Error reading multidatasets.'
               error = -1
               return
            end if
         end do
      end if
      if (allocated(Paths)) deallocate (Paths)

      error = 1
      return

   END SUBROUTINE
!tdReadDatasets
! --------------------------------------------------------------------------
! FUNCTION tdReadActivityScalarAIndex
! PURPOSE  Read all timestep values for a particular index
! NOTES
! --------------------------------------------------------------------------
   SUBROUTINE TD_READ_ACTIVITY_SCALAR_A_INDEX(a_Filename, a_Index, error)
      use const_def, only: READWRITE, READONLY, OVERWRITE

      CHARACTER(LEN=*), INTENT(IN) :: a_Filename
      INTEGER, INTENT(IN)   :: a_Index
      INTEGER, INTENT(OUT)  :: error
      INTEGER status
      INTEGER(XID) xFileId, xDsetsId, xScalarAId
      INTEGER nTimesteps, i
      INTEGER, ALLOCATABLE  :: bActive(:)

      xFileId = NONE
      xDsetsId = NONE
      xScalarAId = NONE

      ! open the file
      call XF_OPEN_FILE(a_Filename, READONLY, xFileId, status)
      if (status < 0) then
         error = -1
         return
      end if

      ! open the dataset group
      call XF_OPEN_GROUP(xFileId, DATASETS_LOCATION, xDsetsId, status)
      if (status >= 0) then
         call XF_OPEN_GROUP(xDsetsId, SCALAR_A_LOCATION, xScalarAId, status)
      end if
      if (status < 0) then
         error = status
         return
      end if

      ! Find out the number of timesteps in the file
      CALL XF_GET_DATASET_NUM_TIMES(xScalarAId, nTimesteps, status)
      if (status < 0) then
         error = status
         return
      end if

      if (nTimesteps < 1) then
         error = -1
         return
      end if

      ! Read the values for the index
      allocate (bActive(nTimesteps))
      call XF_READ_ACTIVE_VALS_AT_INDEX(xScalarAId, a_Index, 1, nTimesteps, &
                                        bActive, status)
      ! output the data
      WRITE (*, *) ''
      WRITE (*, *) 'Reading activity for scalar A slice at index: ', a_Index
      do i = 1, nTimesteps
         WRITE (*, *) bActive(i), ' '
      end do

      deallocate (bActive)

      error = status
      return

   END SUBROUTINE
! tdReadActivityScalarAIndex

! --------------------------------------------------------------------------
! FUNCTION tdReadScalarAIndex
! PURPOSE  Read all timestep values for a particular index
! NOTES
! --------------------------------------------------------------------------
   SUBROUTINE TD_READ_SCALAR_A_INDEX(a_Filename, a_Index, error)
      use const_def, only: READWRITE, READONLY, OVERWRITE

      CHARACTER(LEN=*), INTENT(IN) :: a_Filename
      INTEGER, INTENT(IN)   :: a_Index
      INTEGER, INTENT(OUT)  :: error
      INTEGER status
      INTEGER(XID) xFileId, xDsetsId, xScalarAId
      INTEGER nTimesteps, i
      REAL, ALLOCATABLE :: fValues(:)

      xFileId = NONE
      xDsetsId = NONE
      xScalarAId = NONE

      ! open the file
      call XF_OPEN_FILE(a_Filename, READONLY, xFileId, status)
      if (status < 0) then
         error = -1
         return
      end if

      ! open the dataset group
      call XF_OPEN_GROUP(xFileId, DATASETS_LOCATION, xDsetsId, status)
      if (status >= 0) then
         call XF_OPEN_GROUP(xDsetsId, SCALAR_A_LOCATION, xScalarAId, status)
      end if
      if (status < 0) then
         error = status
         return
      end if

      ! Find out the number of timesteps in the file
      call XF_GET_DATASET_NUM_TIMES(xScalarAId, nTimesteps, status)
      if (status < 0) then
         error = status
         return
      end if

      if (nTimesteps < 1) then
         error = -1
         return
      end if

      ! Read the values for the index
      allocate (fValues(nTimesteps))
      call XF_READ_SCALAR_VALUES_AT_INDEX(xScalarAId, a_Index, 1, nTimesteps, &
                                          fValues, status)

      ! output the data
      WRITE (*, *) ''
      WRITE (*, *) 'Reading scalar A slice at index: ', a_Index
      do i = 1, nTimesteps
         WRITE (*, *) fValues(i), ' '
      end do

      deallocate (fValues)

      error = status
      return

   END SUBROUTINE
! tdReadScalarAtIndex

! --------------------------------------------------------------------------
! FUNCTION tdWriteScalarA
! PURPOSE  Write scalar Dataset to an HDF5 File
! NOTES    This tests dynamic data sets, and activity
!          This dataset is dynamic concentrations (mg/L) with output times
!          in minutes.
!          Dataset is for a mesh and so nActive is the number of elements
!          which is not the same as the nValues which would be number of nodes
!          reads/writes a reference time in julian days
! --------------------------------------------------------------------------
   SUBROUTINE TD_WRITE_SCALAR_A(a_Filename, a_Compression, error)
      use const_def, only: READWRITE, READONLY, OVERWRITE

      CHARACTER(LEN=*), INTENT(IN) :: a_Filename
      INTEGER, INTENT(IN)          :: a_Compression
      INTEGER, INTENT(OUT)         :: error
      INTEGER(XID) xFileId, xDsetsId, xScalarAId, xCoordId
      INTEGER nValues, nTimes, nActive
      REAL(DOUBLE) dTime, dJulianReftime
      INTEGER iTimestep, iActive, iHpgnZone
      REAL fValues(10) ! nValues
      INTEGER*1 bActivity(10) ! activity
      INTEGER i, status

      ! initialize the data
      nValues = 10
      nTimes = 3
      nActive = 8
      dTime = 0.0

      ! 5th item in data set is always inactive, others active
      do iActive = 1, nActive
         bActivity(iActive) = 1
      end do
      bActivity(6) = 0

      ! create the file
      call XF_CREATE_FILE(a_Filename, OVERWRITE, xFileId, error)
      if (error .LT. 0) then
         ! close the file
         call XF_CLOSE_FILE(xFileId, error)
         return
      end if

      ! create the group where we will put all the datasets
      call XF_CREATE_GENERIC_GROUP(xFileId, DATASETS_LOCATION, xDsetsId, status)
      if (status < 0) then
         call XF_CLOSE_FILE(xFileId, error)
         error = -1
         return
      end if

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
      end if

      ! Add in a reftime.  This is a julian day for:
      ! noon July 1, 2003
      dJulianReftime = 2452822.0; 
      call XF_WRITE_REFTIME(xScalarAId, dJulianReftime, status)
      if (status < 0) then
         call XF_CLOSE_GROUP(xScalarAId, error)
         call XF_CLOSE_GROUP(xDsetsId, error)
         call XF_CLOSE_FILE(xFileId, error)
      end if

      ! Loop through timesteps adding them to the file
      do iTimestep = 1, nTimes
         ! We will have an 0.5 hour timestep
         dTime = iTimestep*0.5

         fValues(1) = dTime
         do i = 2, nValues
            fValues(i) = fValues(i - 1)*2.5
         end do

         ! write the dataset array values
         call XF_WRITE_SCALAR_TIMESTEP(xScalarAId, dTime, nValues, fValues, error)
         if (error .GE. 0) then
            ! write activity array
            call XF_WRITE_ACTIVITY_TIMESTEP(xScalarAId, nActive, bActivity, error)
         end if
      end do

      ! Write Coordinate file - for ScalarA, we will set the coordinate system
      !   to be Geographic HPGN, with HPGN settings written to the file.
      call XF_CREATE_COORDINATE_GROUP(xFileId, xCoordId, status)
      if (status < 0) then
         call XF_CLOSE_GROUP(xScalarAId, error)
         call XF_CLOSE_GROUP(xDsetsId, error)
         call XF_CLOSE_FILE(xFileId, error)
         error = -1
         return
      end if

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
! tdWriteScalarA

! --------------------------------------------------------------------------
! FUNCTION tdWriteScalarAPIECES
! PURPOSE  Write scalar Dataset to an HDF5 File
! NOTES    This tests dynamic data sets, and activity
!          This dataset is dynamic concentrations (mg/L) with output times
!          in minutes.
!          Dataset is for a mesh and so nActive is the number of elements
!          which is not the same as the nValues which would be number of nodes
!          reads/writes a reference time in julian days
! --------------------------------------------------------------------------
   SUBROUTINE TD_WRITE_SCALAR_A_PIECES(a_Filename, a_Compression, error)
      use const_def, only: READWRITE, READONLY, OVERWRITE

      CHARACTER(LEN=*), INTENT(IN) :: a_Filename
      INTEGER, INTENT(IN)          :: a_Compression
      INTEGER, INTENT(OUT)         :: error
      INTEGER(XID) xFileId, xDsetsId, xScalarAId, xCoordId
      INTEGER nValues, nTimes, nActive
      REAL(DOUBLE) dTime, dJulianReftime
      INTEGER iTimestep, iActive, iHpgnZone
      REAL fValues(10) ! nValues
      INTEGER*1 bActivity(10) ! activity
      INTEGER i, status
      REAL minvalue, maxvalue
      INTEGER(XID) timestepId
      INTEGER(XID) activeTs

      ! initialize the data
      nValues = 10
      nTimes = 3
      nActive = 8
      dTime = 0.0

      ! 5th item in data set is always inactive, others active
      do iActive = 1, nActive
         bActivity(iActive) = 1
      end do
      bActivity(6) = 0

      ! create the file
      call XF_CREATE_FILE(a_Filename, OVERWRITE, xFileId, error)
      if (error .LT. 0) then
         ! close the file
         call XF_CLOSE_FILE(xFileId, error)
         return
      end if

      ! create the group where we will put all the datasets
      call XF_CREATE_GENERIC_GROUP(xFileId, DATASETS_LOCATION, xDsetsId, status)
      if (status < 0) then
         call XF_CLOSE_FILE(xFileId, error)
         error = -1
         return
      end if

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
      end if

      ! Add in a reftime.  This is a julian day for:
      ! noon July 1, 2003
      dJulianReftime = 2452822.0; 
      call XF_WRITE_REFTIME(xScalarAId, dJulianReftime, status)
      if (status < 0) then
         call XF_CLOSE_GROUP(xScalarAId, error)
         call XF_CLOSE_GROUP(xDsetsId, error)
         call XF_CLOSE_FILE(xFileId, error)
      end if

      ! Loop through timesteps adding them to the file
      do iTimestep = 1, nTimes
         ! We will have an 0.5 hour timestep
         dTime = iTimestep*0.5

         fValues(1) = dTime
         minvalue = fValues(1)
         maxvalue = fValues(1)
         do i = 2, nValues
            fValues(i) = fValues(i - 1)*2.5
            minvalue = min(minvalue, fValues(i))
            maxvalue = max(maxvalue, fValues(i))
         end do

         ! write the dataset array values
         call XF_INITIALIZE_SCALAR_TIMESTEP(xScalarAId, dTime, nValues, minvalue, &
                                            maxvalue, timestepId, error)

         ! write data in pairs
         do i = 1, nValues, +2
            call XF_WRITE_SCALAR_TIMESTEP_PORTION(xScalarAId, timestepId, 2, i, &
                                                  fValues(i), error)
         end do

         if (error .GE. 0) then
            ! write activity array
            call XF_INITIALIZE_ACTIVITY_TIMESTEP(xScalarAId, nActive, activeTs, error)

            do i = 1, nActive, +2
               call XF_WRITE_ACTIVITY_TIMESTEP_PORTION(xScalarAId, activeTs, 2, &
                                                       i, bActivity(i), error)
            end do

         end if
      end do

      ! Write Coordinate file - for ScalarA, we will set the coordinate system
      !   to be Geographic HPGN, with HPGN settings written to the file.
      call XF_CREATE_COORDINATE_GROUP(xFileId, xCoordId, status)
      if (status < 0) then
         call XF_CLOSE_GROUP(xScalarAId, error)
         call XF_CLOSE_GROUP(xDsetsId, error)
         call XF_CLOSE_FILE(xFileId, error)
         error = -1
         return
      end if

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
! tdWriteScalarAPieces

! --------------------------------------------------------------------------
! FUNCTION tdWriteScalarAPIECESAltMinMax
! PURPOSE  Write scalar Dataset to an HDF5 File
! NOTES    This tests dynamic data sets, and activity
!          This dataset is dynamic concentrations (mg/L) with output times
!          in minutes.
!          Dataset is for a mesh and so nActive is the number of elements
!          which is not the same as the nValues which would be number of nodes
!          reads/writes a reference time in julian days
! --------------------------------------------------------------------------
   SUBROUTINE TD_WRITE_SCALAR_A_PIECES_ALT_MIN_MAX(a_Filename, a_Compression, error)
      use const_def, only: READWRITE, READONLY, OVERWRITE

      CHARACTER(LEN=*), INTENT(IN) :: a_Filename
      INTEGER, INTENT(IN)          :: a_Compression
      INTEGER, INTENT(OUT)         :: error
      INTEGER(XID) xFileId, xDsetsId, xScalarAId, xCoordId
      INTEGER nValues, nTimes, nActive
      REAL(DOUBLE) dTime, dJulianReftime
      INTEGER iTimestep, iActive, iHpgnZone
      REAL fValues(10) ! nValues
      INTEGER*1 bActivity(10) ! activity
      INTEGER i, status
      REAL minvalue, maxvalue
      INTEGER(XID) timestepId
      INTEGER(XID) activeTs

      ! initialize the data
      nValues = 10
      nTimes = 3
      nActive = 8
      dTime = 0.0

      ! 5th item in data set is always inactive, others active
      do iActive = 1, nActive
         bActivity(iActive) = 1
      end do
      bActivity(6) = 0

      ! create the file
      call XF_CREATE_FILE(a_Filename, OVERWRITE, xFileId, error)
      if (error .LT. 0) then
         ! close the file
         call XF_CLOSE_FILE(xFileId, error)
         return
      end if

      ! create the group where we will put all the datasets
      call XF_CREATE_GENERIC_GROUP(xFileId, DATASETS_LOCATION, xDsetsId, status)
      if (status < 0) then
         call XF_CLOSE_FILE(xFileId, error)
         error = -1
         return
      end if

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
      end if

      ! Add in a reftime.  This is a julian day for:
      ! noon July 1, 2003
      dJulianReftime = 2452822.0; 
      call XF_WRITE_REFTIME(xScalarAId, dJulianReftime, status)
      if (status < 0) then
         call XF_CLOSE_GROUP(xScalarAId, error)
         call XF_CLOSE_GROUP(xDsetsId, error)
         call XF_CLOSE_FILE(xFileId, error)
      end if

      ! Loop through timesteps adding them to the file
      do iTimestep = 1, nTimes
         ! We will have an 0.5 hour timestep
         dTime = iTimestep*0.5

         fValues(1) = dTime
         minvalue = fValues(1)
         maxvalue = fValues(1)
         do i = 2, nValues
            fValues(i) = fValues(i - 1)*2.5
            minvalue = min(minvalue, fValues(i))
            maxvalue = max(maxvalue, fValues(i))
         end do

         ! write the dataset array values
         call XF_INITIALIZE_SCALAR_TIMESTEP(xScalarAId, dTime, nValues, minvalue, &
                                            maxvalue, timestepId, error)

         ! write data in pairs
         do i = 1, nValues, +2
            call XF_WRITE_SCALAR_TIMESTEP_PORTION(xScalarAId, timestepId, 2, i, &
                                                  fValues(i), error)
         end do

         minvalue = 0.1111*timestepId
         maxvalue = 1111*timestepId
         call XF_SET_DATASET_TIMESTEP_MIN_MAX(xScalarAId, timestepId, minvalue, &
                                              maxvalue, error)

         if (error .GE. 0) then
            ! write activity array
            call XF_INITIALIZE_ACTIVITY_TIMESTEP(xScalarAId, nActive, activeTs, error)

            do i = 1, nActive, +2
               call XF_WRITE_ACTIVITY_TIMESTEP_PORTION(xScalarAId, activeTs, 2, &
                                                       i, bActivity(i), error)
            end do

         end if
      end do

      ! Write Coordinate file - for ScalarA, we will set the coordinate system
      !   to be Geographic HPGN, with HPGN settings written to the file.
      call XF_CREATE_COORDINATE_GROUP(xFileId, xCoordId, status)
      if (status < 0) then
         call XF_CLOSE_GROUP(xScalarAId, error)
         call XF_CLOSE_GROUP(xDsetsId, error)
         call XF_CLOSE_FILE(xFileId, error)
         error = -1
         return
      end if

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
! TD_WRITE_SCALAR_A_PIECES_ALT_MIN_MAX

! --------------------------------------------------------------------------
! FUNCTION TD_WRITE_SCALAR_B
! PURPOSE  Write scalar Dataset to an HDF5 File
! NOTES    This tests dynamic data sets, and activity
!          This dataset is dynamic concentrations (mg/L) with output times
!          in minutes.
!          Dataset is for a mesh and so nActive is the number of elements
!          which is not the same as the nValues which would be number of nodes
!          reads/writes a reference time in julian days
! --------------------------------------------------------------------------
   SUBROUTINE TD_WRITE_SCALAR_B(a_Filename, a_Compression, a_Overwrite, error)
      use const_def, only: READWRITE, READONLY, OVERWRITE

      CHARACTER(LEN=*), INTENT(IN) :: a_Filename
      INTEGER, INTENT(IN)          :: a_Compression
      LOGICAL, INTENT(IN)          :: a_Overwrite
      INTEGER, INTENT(OUT)         :: error
      INTEGER(XID) xFileId, xDsetsId, xScalarBId, xCoordId
      INTEGER nValues, nTimes, nActive
      REAL(DOUBLE) dTime, dJulianReftime
      INTEGER iTimestep, iActive
      REAL fValues(10) ! nValues
      INTEGER*1 bActivity(10) ! activity
      INTEGER i, status

      ! initialize the data
      nValues = 10
      nTimes = 3
      nActive = 8
      dTime = 0.0
      i = 0

      ! 5th item in data set is always inactive, others active
      do iActive = 1, nActive
         bActivity(iActive) = 1
      end do
      bActivity(6) = 0

      if (a_Overwrite) then
         ! open the already-existing file
         call XF_OPEN_FILE(a_Filename, READWRITE, xFileId, status)
         if (status < 0) then
            error = -1
            return
         end if
         ! open the group where we have all the datasets
         call XF_OPEN_GROUP(xFileId, DATASETS_LOCATION, xDsetsId, status)
         if (status < 0) then
            call XF_CLOSE_FILE(xFileId, error)
            error = -1
            return
         end if
      else
         ! create the file
         call XF_CREATE_FILE(a_Filename, OVERWRITE, xFileId, error)
         if (error .LT. 0) then
            ! close the file
            call XF_CLOSE_FILE(xFileId, error)
            return
         end if

         ! create the group where we will put all the datasets
         call XF_CREATE_GENERIC_GROUP(xFileId, DATASETS_LOCATION, xDsetsId, status)
         if (status < 0) then
            call XF_CLOSE_FILE(xFileId, error)
            error = -1
            return
         end if
      end if

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
      end if

      ! Add in a reftime.  This is a julian day for:
      ! noon July 1, 2003
      dJulianReftime = 2452822.0; 
      call XF_WRITE_REFTIME(xScalarBId, dJulianReftime, status)
      if (status < 0) then
         call XF_CLOSE_GROUP(xScalarBId, error)
         call XF_CLOSE_GROUP(xDsetsId, error)
         call XF_CLOSE_FILE(xFileId, error)
      end if

      if (.NOT. a_Overwrite) then
         ! Loop through timesteps adding them to the file
         do iTimestep = 1, nTimes
            ! We will have an 0.5 hour timestep
            dTime = iTimestep*0.5

            fValues(1) = dTime
            do i = 2, nValues
               fValues(i) = fValues(i - 1)*2.5
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
            end if
         end do
      else
         ! Loop through timesteps adding them to the file
         do iTimestep = 1, nTimes
            ! We will have an 1.5 hour timestep
            dTime = iTimestep*1.5

            fValues(1) = dTime
            do i = 2, nValues
               fValues(i) = fValues(i - 1)*1.5
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
            end if
         end do
      end if

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
         end if

         ! Write Coord Info to file
         call XF_SET_HORIZ_DATUM(xCoordId, HORIZ_DATUM_UTM, error)
         call XF_SET_HORIZ_UNITS(xCoordId, COORD_UNITS_METERS, error)

         call XF_SET_VERT_DATUM(xCoordId, VERT_DATUM_LOCAL, error)
         call XF_SET_VERT_UNITS(xCoordId, COORD_UNITS_METERS, error)

         ! write additional information - we'll use the max value for this test
         call XF_SET_UTM_ZONE(xCoordId, UTM_ZONE_MAX, error)

         call XF_CLOSE_GROUP(xCoordId, error)
         xCoordId = 0
      end if

      ! close the dataset
      call XF_CLOSE_GROUP(xScalarBId, error)
      call XF_CLOSE_GROUP(xDsetsId, error)
      call XF_CLOSE_FILE(xFileId, error)

      error = 1
      return
   END SUBROUTINE
! tdWriteScalarB
!------------------------------------------------------------------------------
!  FUNCTION TD_WRITE_COORDS_TO_MULTI
!  PURPOSE  Write coordinate system to a multidataset file
!  NOTES
!------------------------------------------------------------------------------
   SUBROUTINE TD_WRITE_COORDS_TO_MULTI(a_xFileId, error)
      INTEGER(XID), INTENT(IN) :: a_xFileId
      INTEGER, INTENT(OUT)       :: error
      INTEGER(XID) xCoordId
      INTEGER status

      ! Write Coordinate file - for Multidatasets, we will set the coordinate system
      !   to be UTM, with UTM Zone settings written to the file.
      call XF_CREATE_COORDINATE_GROUP(a_xFileId, xCoordId, status)
      if (status < 0) then
         error = status
         return
      end if

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
! FUNCTION tdWriteScalarAToMulti
! PURPOSE  Write scalar Dataset to an HDF5 File
! NOTES    This tests dynamic data sets, and activity
!          This dataset is dynamic concentrations (mg/L) with output times
!          in minutes.
!          Dataset is for a mesh and so nActive is the number of elements
!          which is not the same as the nValues which would be number of nodes
!          reads/writes a reference time in julian days
! --------------------------------------------------------------------------
   SUBROUTINE TD_WRITE_SCALAR_A_TO_MULTI(a_GroupID, status)
      ! CHARACTER(LEN=*), INTENT(IN) :: a_Filename
      ! INTEGER, INTENT(IN)          :: a_Compression
      ! INTEGER, INTENT(OUT)         :: error
      INTEGER(XID) xFileId, xDsetsId, xScalarAId
      INTEGER(XID) a_GroupID
      INTEGER nValues, nTimes, nActive
      REAL(DOUBLE) dTime, dJulianReftime
      INTEGER iTimestep, iActive
      REAL fValues(10) ! nValues
      INTEGER*1 bActivity(10) ! activity
      INTEGER i, status

      ! initialize the data
      nValues = 10
      nTimes = 3
      nActive = 8
      dTime = 0.0

      ! 5th item in data set is always inactive, others active
      do iActive = 1, nActive
         bActivity(iActive) = 1
      end do
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
      end if

      ! Add in a reftime.  This is a julian day for:
      ! noon July 1, 2003
      dJulianReftime = 2452822.0; 
      call XF_WRITE_REFTIME(xScalarAId, dJulianReftime, status)
      if (status < 0) then
         call XF_CLOSE_GROUP(xScalarAId, status)
         call XF_CLOSE_GROUP(xDsetsId, status)
         call XF_CLOSE_FILE(xFileId, status)
      end if

      ! Loop through timesteps adding them to the file
      do iTimestep = 1, nTimes
         ! We will have an 0.5 hour timestep
         dTime = iTimestep*0.5

         fValues(1) = dTime
         do i = 2, nValues
            fValues(i) = fValues(i - 1)*2.5
         end do

         ! write the dataset array values
         call XF_WRITE_SCALAR_TIMESTEP(xScalarAId, dTime, nValues, fValues, status)
         if (status .GE. 0) then
            ! write activity array
            call XF_WRITE_ACTIVITY_TIMESTEP(xScalarAId, nActive, bActivity, status)
         end if
      end do

      ! close the dataset
      call XF_CLOSE_GROUP(xScalarAId, status)
      !call XF_CLOSE_GROUP(a_GroupID, status)
      !call XF_CLOSE_FILE(a_FileID, status)

      return
   END SUBROUTINE
! tdWriteScalarAToMulti
! --------------------------------------------------------------------------
! FUNCTION tdReadVector2DAIndex
! PURPOSE  Read all timestep values for a particular index
! NOTES
! --------------------------------------------------------------------------
   SUBROUTINE TD_READ_VECTOR2D_A_INDEX(a_Filename, a_Index, error)
      use const_def, only: READWRITE, READONLY, OVERWRITE

      CHARACTER(LEN=*), INTENT(IN) :: a_Filename
      INTEGER, INTENT(IN)   :: a_Index
      INTEGER, INTENT(OUT)  :: error
      INTEGER status
      INTEGER(XID) xFileId, xDsetsId, xVector2DA
      INTEGER nTimesteps, i
      REAL, ALLOCATABLE     :: fValues(:)

      xFileId = NONE
      xDsetsId = NONE
      xVector2DA = NONE

      ! open the file
      call XF_OPEN_FILE(a_Filename, READONLY, xFileId, status)
      if (status < 0) then
         error = -1
         return
      end if

      ! open the dataset group
      call XF_OPEN_GROUP(xFileId, DATASETS_LOCATION, xDsetsId, status)
      if (status >= 0) then
         call XF_OPEN_GROUP(xDsetsId, VECTOR2D_A_LOCATION, xVector2DA, status)
      end if
      if (status < 0) then
         error = status
         return
      end if

      ! Find out the number of timesteps in the file
      call XF_GET_DATASET_NUM_TIMES(xVector2DA, nTimesteps, status)
      if (status < 0) then
         error = status
         return
      end if

      if (nTimesteps < 1) then
         error = -1
         return
      end if

      ! Read the values for the index
      allocate (fValues(nTimesteps*2))
      call XF_READ_VECTOR_VALUES_AT_INDEX(xVector2DA, a_Index, 1, nTimesteps, 2, &
                                          fValues, status)

      ! output the data
      WRITE (*, *) ''
      WRITE (*, *) 'Reading vector 2D A slice at index: ', a_Index
      do i = 1, nTimesteps
         WRITE (*, *) fValues(i*2 - 1), ' ', fValues(i*2)
      end do
      WRITE (*, *) ''

      deallocate (fValues)

      error = status
      return

   END SUBROUTINE
!tdReadVector2DAIndex

! --------------------------------------------------------------------------
! FUNCTION tdWriteVector2D_A
! PURPOSE  Write scalar Dataset to an HDF5 File
! NOTES    This tests dynamic data sets, and activity
!          This dataset is dynamic concentrations (mg/L) with output times
!          in minutes.
!          Dataset is for a mesh and so nActive is the number of elements
!          which is not the same as the nValues which would be number of nodes
!          reads/writes a reference time in julian days
! --------------------------------------------------------------------------
   SUBROUTINE TD_WRITE_VECTOR2D_A(a_Filename, a_Compression, error)
      use const_def, only: READWRITE, READONLY, OVERWRITE

      CHARACTER(LEN=*), INTENT(IN) :: a_Filename
      INTEGER, INTENT(IN)          :: a_Compression
      INTEGER, INTENT(OUT)         :: error
      INTEGER(XID) xFileId, xDsetsId, xVector2D_A, xCoordId
      INTEGER nValues, nTimes, nComponents, nActive
      REAL(DOUBLE) dTime
      INTEGER iTimestep, iActive
      REAL, DIMENSION(2, 100) :: fValues ! nComponents, nValues
      INTEGER*1 bActivity(100) ! activity
      INTEGER i, j, status
      INTEGER iHpgnZone

      ! initialize the data
      nComponents = 2
      nValues = 100
      nTimes = 6
      nActive = 75
      dTime = 0.0

      ! 5th item in data set is always inactive, others active
      bActivity(1) = 0
      do iActive = 2, nActive
         if (mod(iActive - 1, 3) == 0) then
            bActivity(iActive) = 0
         else
            bActivity(iActive) = 1
         end if
      end do

      ! create the file
      call XF_CREATE_FILE(a_Filename, OVERWRITE, xFileId, error)
      if (error .LT. 0) then
         ! close the dataset
         call XF_CLOSE_FILE(xFileId, error)
         return
      end if

      ! create the group where we will put all the datasets
      call XF_CREATE_GENERIC_GROUP(xFileId, DATASETS_LOCATION, xDsetsId, status)
      if (status < 0) then
         call XF_CLOSE_FILE(xFileId, error)
         error = -1
         return
      end if

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
      end if

      ! Loop through timesteps adding them to the file
      do iTimestep = 1, nTimes
         ! We will have an 0.5 hour timestep
         dTime = iTimestep*0.5

         do i = 1, nValues
            do j = 1, nComponents
               fValues(j, i) = ((i - 1)*nComponents + j)*dTime
            end do
         end do

         ! write the dataset array values
         call XF_WRITE_VECTOR_TIMESTEP(xVector2D_A, dTime, nValues, nComponents, &
                                       fValues, error)
         if (error .GE. 0) then
            ! write activity array
            call XF_WRITE_ACTIVITY_TIMESTEP(xVector2D_A, nActive, bActivity, error)
         end if
      end do

      ! Write Coordinate file - for Vector2D_A, we will set the coordinate system
      !   to be Geographic HPGN, with HPGN settings written to the file.
      call XF_CREATE_COORDINATE_GROUP(xFileId, xCoordId, status)
      if (status < 0) then
         call XF_CLOSE_GROUP(xVector2D_A, error)
         call XF_CLOSE_GROUP(xDsetsId, error)
         call XF_CLOSE_FILE(xFileId, error)
         error = -1
         return
      end if

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
! tdWriteVector2D_A

! --------------------------------------------------------------------------
! FUNCTION tdWriteVector2D_A_PIECES
! PURPOSE  Write scalar Dataset to an HDF5 File
! NOTES    This tests dynamic data sets, and activity
!          This dataset is dynamic concentrations (mg/L) with output times
!          in minutes.
!          Dataset is for a mesh and so nActive is the number of elements
!          which is not the same as the nValues which would be number of nodes
!          reads/writes a reference time in julian days
! --------------------------------------------------------------------------
   SUBROUTINE TD_WRITE_VECTOR2D_A_PIECES(a_Filename, a_Compression, error)
      use const_def, only: READWRITE, READONLY, OVERWRITE

      CHARACTER(LEN=*), INTENT(IN) :: a_Filename
      INTEGER, INTENT(IN)          :: a_Compression
      INTEGER, INTENT(OUT)         :: error
      INTEGER(XID) xFileId, xDsetsId, xVector2D_A, xCoordId
      INTEGER nValues, nTimes, nComponents, nActive
      REAL(DOUBLE) dTime
      INTEGER iTimestep, iActive
      REAL*4, DIMENSION(2, 100) :: fValues ! nComponents, nValues
      INTEGER*1 bActivity(100) ! activity
      INTEGER i, j, status
      INTEGER iHpgnZone
      INTEGER nValuesToWrite, nComponentsToWrite, startComponent
      REAL*4 minvalue, maxvalue
      INTEGER(XID) timeId
      REAL*8 mag

      ! initialize the data
      nComponents = 2
      nValues = 100
      nTimes = 6
      nActive = 75
      dTime = 0.0

      ! 5th item in data set is always inactive, others active
      bActivity(1) = 0
      do iActive = 2, nActive
         if (mod(iActive - 1, 3) == 0) then
            bActivity(iActive) = 0
         else
            bActivity(iActive) = 1
         end if
      end do

      ! create the file
      call XF_CREATE_FILE(a_Filename, OVERWRITE, xFileId, error)
      if (error .LT. 0) then
         ! close the dataset
         call XF_CLOSE_FILE(xFileId, error)
         return
      end if

      ! create the group where we will put all the datasets
      call XF_CREATE_GENERIC_GROUP(xFileId, DATASETS_LOCATION, xDsetsId, status)
      if (status < 0) then
         call XF_CLOSE_FILE(xFileId, error)
         error = -1
         return
      end if

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
      end if

      ! Loop through timesteps adding them to the file
      do iTimestep = 1, nTimes
         ! We will have an 0.5 hour timestep
         dTime = iTimestep*0.5

         do i = 1, nValues
            do j = 1, nComponents
               fValues(j, i) = ((i - 1)*nComponents + j)*dTime
            end do
         end do

         minvalue = 99999.0
         maxvalue = 0.0
         do i = 1, nValues
            mag = 0.0
            do j = 1, nComponents
               mag = mag + fValues(j, i)**2
            end do
            mag = mag**0.5

            minvalue = min(minvalue, mag)
            maxvalue = max(maxvalue, mag)
         end do

         ! write the dataset array values
         call XF_INITIALIZE_VECTOR_TIMESTEP(xVector2D_A, dTime, nValues, nComponents, &
                                            minvalue, maxvalue, timeId, error)

         nValuesToWrite = 2
         nComponentsToWrite = 2
         startComponent = 1

         do i = 1, nValues, +2
            call XF_WRITE_VECTOR_TIMESTEP_PORTION(xVector2D_A, timeId, nValuesToWrite, &
                                                  nComponentsToWrite, i, startComponent, fValues(1, i), error)
         end do

         if (error .GE. 0) then
            ! write activity array
            call XF_WRITE_ACTIVITY_TIMESTEP(xVector2D_A, nActive, bActivity, error)
         end if
      end do

      ! Write Coordinate file - for Vector2D_A, we will set the coordinate system
      !   to be Geographic HPGN, with HPGN settings written to the file.
      call XF_CREATE_COORDINATE_GROUP(xFileId, xCoordId, status)
      if (status < 0) then
         call XF_CLOSE_GROUP(xVector2D_A, error)
         call XF_CLOSE_GROUP(xDsetsId, error)
         call XF_CLOSE_FILE(xFileId, error)
         error = -1
         return
      end if

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
! tdWriteVector2D_A_PIECES

! --------------------------------------------------------------------------
! FUNCTION TD_WRITE_VECTOR2D_B
! PURPOSE  Write scalar Dataset to an HDF5 File
! NOTES    This tests dynamic data sets, and activity
!          This dataset is dynamic concentrations (mg/L) with output times
!          in minutes.
!          Dataset is for a mesh and so nActive is the number of elements
!          which is not the same as the nValues which would be number of nodes
!          reads/writes a reference time in julian days
! --------------------------------------------------------------------------
   SUBROUTINE TD_WRITE_VECTOR2D_B(a_Filename, a_Compression, a_Overwrite, error)
      use const_def, only: READONLY, READWRITE, OVERWRITE

      CHARACTER(LEN=*), INTENT(IN) :: a_Filename
      INTEGER, INTENT(IN)          :: a_Compression
      LOGICAL, INTENT(IN)          :: a_Overwrite
      INTEGER, INTENT(OUT)         :: error
      INTEGER(XID) xFileId, xDsetsId, xVector2D_B, xCoordId
      INTEGER nValues, nTimes, nComponents, nActive
      REAL(DOUBLE) dTime
      INTEGER iTimestep, iActive
      REAL, DIMENSION(2, 100) :: fValues
      INTEGER*1 bActivity(100)
      INTEGER i, j, status

      ! initialize the data
      nComponents = 2
      nValues = 100
      nTimes = 6
      nActive = 75
      dTime = 0.0

      ! 5th item in data set is always inactive, others active
      bActivity(1) = 0
      do iActive = 2, nActive
         if (mod(iActive - 1, 3) == 0) then
            bActivity(iActive) = 0
         else
            bActivity(iActive) = 1
         end if
      end do

      if (a_Overwrite) then
         ! open the already-existing file
         call XF_OPEN_FILE(a_Filename, READWRITE, xFileId, status)
         if (status < 0) then
            error = -1
            return
         end if
         ! open the group where we have all the datasets
         call XF_OPEN_GROUP(xFileId, DATASETS_LOCATION, xDsetsId, status)
         if (status < 0) then
            call XF_CLOSE_FILE(xFileId, error)
            error = -1
            return
         end if
      else
         ! create the file
         call XF_CREATE_FILE(a_Filename, OVERWRITE, xFileId, error)
         if (error .LT. 0) then
            ! close the dataset
            call XF_CLOSE_FILE(xFileId, error)
            return
         end if

         ! create the group where we will put all the datasets
         call XF_CREATE_GENERIC_GROUP(xFileId, DATASETS_LOCATION, xDsetsId, status)
         if (status < 0) then
            call XF_CLOSE_FILE(xFileId, error)
            error = -1
            return
         end if
      end if

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
      end if

      if (.NOT. a_Overwrite) then
         ! Loop through timesteps adding them to the file
         do iTimestep = 1, nTimes
            ! We will have an 0.5 hour timestep
            dTime = iTimestep*0.5
            do i = 1, nValues
               do j = 1, nComponents
                  fValues(j, i) = ((i - 1)*nComponents + j)*dTime
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
            end if
         end do
      else
         ! Loop through timesteps adding them to the file
         do iTimestep = 1, nTimes
            ! We will have an 1.5 hour timestep
            dTime = iTimestep*1.5
            do i = 1, nValues
               do j = 1, nComponents
                  fValues(j, i) = ((i - 1)*nComponents + j)*dTime
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
            end if
         end do
      end if

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
         end if

         ! write the coordinate data to the file
         call XF_SET_HORIZ_DATUM(xCoordId, HORIZ_DATUM_UTM, error)
         call XF_SET_HORIZ_UNITS(xCoordId, COORD_UNITS_METERS, error)
         call XF_SET_VERT_DATUM(xCoordId, VERT_DATUM_LOCAL, error)
         call XF_SET_VERT_UNITS(xCoordId, COORD_UNITS_METERS, error)

         ! write additional information - we'll use the max UTM zone for the test
         call XF_SET_UTM_ZONE(xCoordId, UTM_ZONE_MAX, error)

         call XF_CLOSE_GROUP(xCoordId, error)
         xCoordId = 0
      end if

      ! close the dataset
      call XF_CLOSE_GROUP(xVector2D_B, error)
      call XF_CLOSE_GROUP(xDsetsId, error)
      call XF_CLOSE_FILE(xFileId, error)

      return
   END SUBROUTINE
! tdWriteVector2D_B

! --------------------------------------------------------------------------
! FUNCTION tdWriteVector2D_AToMulti
! PURPOSE  Write scalar Dataset to an HDF5 File
! NOTES    This tests dynamic data sets, and activity
!          This dataset is dynamic concentrations (mg/L) with output times
!          in minutes.
!          Dataset is for a mesh and so nActive is the number of elements
!          which is not the same as the nValues which would be number of nodes
!          reads/writes a reference time in julian days
! --------------------------------------------------------------------------
   SUBROUTINE TD_WRITE_VECTOR2D_A_TO_MULTI(a_FileID, a_GroupID, status)
      INTEGER(XID) xVector2D_A
      INTEGER(XID) a_FileID, a_GroupID
      INTEGER nValues, nTimes, nComponents, nActive
      REAL(DOUBLE) dTime
      INTEGER iTimestep, iActive
      REAL, DIMENSION(2, 100) :: fValues ! nComponents, nValues
      INTEGER*1 bActivity(100) ! activity
      INTEGER i, j, status

      ! initialize the data
      nComponents = 2
      nValues = 100
      nTimes = 6
      nActive = 75
      dTime = 0.0

      ! 5th item in data set is always inactive, others active
      bActivity(1) = 0
      do iActive = 2, nActive
         if (mod(iActive - 1, 3) == 0) then
            bActivity(iActive) = 0
         else
            bActivity(iActive) = 1
         end if
      end do

      ! Create the vector dataset group
      call XF_CREATE_VECTOR_DATASET(a_GroupID, VECTOR2D_A_LOCATION, 'ft/s', &
                                    TS_SECONDS, NONE, xVector2D_A, status)
      if (status .LT. 0) then
         ! close the dataset
         call XF_CLOSE_GROUP(xVector2D_A, status)
         call XF_CLOSE_GROUP(a_GroupID, status)
         call XF_CLOSE_FILE(a_FileID, status)
         return
      end if

      ! Loop through timesteps adding them to the file
      do iTimestep = 1, nTimes
         ! We will have an 0.5 hour timestep
         dTime = iTimestep*0.5

         do i = 1, nValues
            do j = 1, nComponents
               fValues(j, i) = ((i - 1)*nComponents + j)*dTime
            end do
         end do

         ! write the dataset array values
         call XF_WRITE_VECTOR_TIMESTEP(xVector2D_A, dTime, nValues, nComponents, &
                                       fValues, status)
         if (status .GE. 0) then
            ! write activity array
            call XF_WRITE_ACTIVITY_TIMESTEP(xVector2D_A, nActive, bActivity, status)
         end if
      end do

      ! close the dataset
      call XF_CLOSE_GROUP(xVector2D_A, status)
      return
   END SUBROUTINE
! tdWriteVector2D_AToMulti
! --------------------------------------------------------------------------
! FUNCTION tdiReadScalar
! PURPOSE  Read a scalar from an XMDF file and output information to
!          to a text file
! NOTES
! --------------------------------------------------------------------------
   SUBROUTINE TDI_READ_SCALAR(a_xScalarId, FileUnit, error)
      INTEGER(XID), INTENT(IN) ::  a_xScalarId
      INTEGER, INTENT(IN) ::         FileUnit
      INTEGER, INTENT(OUT) :: error
      INTEGER nTimes, nValues, nActive
      LOGICAL*2 bUseReftime
      INTEGER iTime
      CHARACTER(LEN=100) TimeUnits
      REAL(DOUBLE), ALLOCATABLE :: Times(:)
      REAL, ALLOCATABLE         :: Values(:), Minimums(:), Maximums(:)
      INTEGER, ALLOCATABLE      :: Active(:)
      REAL(DOUBLE) Reftime
      nTimes = NONE
      nValues = NONE
      nActive = None

      ! read the time units
      call XF_GET_DATASET_TIME_UNITS(a_xScalarId, TimeUnits, error)
      if (error < 0) return

      WRITE (FileUnit, *) 'Time units: ', TimeUnits(1:LEN_TRIM(TimeUnits))

      ! see if we are using a reftime
      call XF_USE_REFTIME(a_xScalarId, bUseReftime, error)
      if (error < 0) then
         return
      end if
      if (bUseReftime) then
         call XF_READ_REFTIME(a_xScalarId, Reftime, error)
         if (error < 0) then
            return
         end if
         WRITE (FileUnit, *) 'Reftime: ', Reftime
      end if

      ! read in the number of values and number of active values
      call XF_GET_DATASET_NUMVALS(a_xScalarId, nValues, error)
      if (error .GE. 0) then
         call XF_GET_DATASET_NUMACTIVE(a_xScalarId, nActive, error)
      end if
      if (error .LT. 0) return

      if (nValues <= 0) then
         WRITE (FileUnit, *) 'No data to read in.'
         error = -1
         return
      end if

      ! read in the number of times
      call XF_GET_DATASET_NUM_TIMES(a_xScalarId, nTimes, error)
      if (error < 0) then
         return
      end if

      ! Read in the individual time values
      allocate (Times(nTimes))

      call XF_GET_DATASET_TIMES(a_xScalarId, nTimes, Times, error)
      if (error < 0) return

      ! Read in the minimum and maximum values
      allocate (Minimums(nTimes))
      allocate (Maximums(nTimes))

      call XF_GET_DATASET_MINS(a_xScalarId, nTimes, Minimums, error)
      if (error >= 0) then
         call XF_GET_DATASET_MAXS(a_xScalarId, nTimes, Maximums, error)
      end if
      if (error < 0) then
         deallocate (Times)
         deallocate (Minimums)
         deallocate (Maximums)
         return
      end if

      allocate (Values(nValues))
      if (nActive .GT. 0) then
         allocate (Active(nActive))
      end if

      WRITE (FileUnit, *) 'Number Timesteps: ', nTimes
      WRITE (FileUnit, *) 'Number Values: ', nValues
      WRITE (FileUnit, *) 'Number Active: ', nActive
      WRITE (FileUnit, *) ''

      ! loop through the timesteps, read the values and active values and write
      ! them to the text file
      do iTime = 1, nTimes
         call XF_READ_SCALAR_VALUES_TIMESTEP(a_xScalarId, iTime, nValues, Values, error)
         if (error >= 0 .AND. nActive > 0) then
            call XF_READ_ACTIVITY_TIMESTEP(a_xScalarId, iTime, nActive, Active, error)
         end if

         ! Write the time, min, max, values and active values to the text output
         ! file.
         WRITE (FileUnit, *) 'Timestep at  ', Times(iTime)
         WRITE (FileUnit, *) 'Min: ', Minimums(iTime)
         WRITE (FileUnit, *) 'Max: ', Maximums(iTime)

         WRITE (FileUnit, *) 'Values:'
         WRITE (FileUnit, *) Values(1:nValues)
         WRITE (FileUnit, *) ''

         WRITE (FileUnit, *) 'Activity:'
         WRITE (FileUnit, *) Active(1:nActive)
         WRITE (FileUnit, *) ''
      end do

      if (allocated(Times)) then
         deallocate (Times)
      end if

      if (allocated(Minimums)) then
         deallocate (Minimums)
      end if

      if (allocated(Maximums)) then
         deallocate (Maximums)
      end if

      if (allocated(Values)) then
         deallocate (Values)
      end if

      if (allocated(Active)) then
         deallocate (Active)
      end if

      return
   END SUBROUTINE
! tdiReadScalar

! --------------------------------------------------------------------------
! FUNCTION TDI_READ_VECTOR
! PURPOSE  Read a vector from an XMDF file and output information to
!          to a text file
! NOTES
! --------------------------------------------------------------------------
   SUBROUTINE TDI_READ_VECTOR(a_xVectorId, FileUnit, error)
      INTEGER(XID), INTENT(IN) ::  a_xVectorId
      INTEGER, INTENT(IN) ::  FileUnit
      INTEGER, INTENT(OUT) :: error
      INTEGER nTimes, nValues, nActive, nComponents
      INTEGER iTime, i
      LOGICAL*2 bUseReftime
      CHARACTER(LEN=100) TimeUnits
      REAL(DOUBLE), ALLOCATABLE :: Times(:)
      REAL, ALLOCATABLE, DIMENSION(:, :) :: Values
      REAL, ALLOCATABLE         :: Minimums(:), Maximums(:)
      INTEGER, ALLOCATABLE      :: Active(:)
      REAL(DOUBLE) Reftime

      nTimes = NONE
      nValues = NONE
      nActive = NONE
      nComponents = NONE

      ! read the time units
      call XF_GET_DATASET_TIME_UNITS(a_xVectorId, TimeUnits, error)
      if (error < 0) return

      WRITE (FileUnit, *) 'Time units: ', TimeUnits(1:LEN_TRIM(TimeUnits))

      ! see if we are using a reftime
      call XF_USE_REFTIME(a_xVectorId, bUseReftime, error)
      if (error < 0) then
         return
      end if
      if (bUseReftime) then
         call XF_READ_REFTIME(a_xVectorId, Reftime, error)
         if (error < 0) then
            return
         end if
         WRITE (FileUnit, *) 'Reftime: ', Reftime
      end if

      ! read in the number of values and number of active values
      call XF_GET_DATASET_NUMVALS(a_xVectorId, nValues, error)
      if (error .GE. 0) then
         call XF_GET_DATASET_NUMCOMPONENTS(a_xVectorId, nComponents, error)
         if (error .GE. 0) then
            call XF_GET_DATASET_NUMACTIVE(a_xVectorId, nActive, error)
         end if
      end if
      if (error .LT. 0) return

      if (nValues <= 0) then
         WRITE (FileUnit, *) 'No data to read in.'
         error = -1
         return
      end if

      ! read in the number of times
      call XF_GET_DATASET_NUM_TIMES(a_xVectorId, nTimes, error)
      if (error < 0) then
         return
      end if

      ! Read in the individual time values
      allocate (Times(nTimes))

      call XF_GET_DATASET_TIMES(a_xVectorId, nTimes, Times, error)
      if (error < 0) return

      ! Read in the minimum and maximum values
      allocate (Minimums(nTimes))
      allocate (Maximums(nTimes))

      call XF_GET_DATASET_MINS(a_xVectorId, nTimes, Minimums, error)
      if (error >= 0) then
         call XF_GET_DATASET_MAXS(a_xVectorId, nTimes, Maximums, error)
      end if
      if (error < 0) then
         deallocate (Times)
         deallocate (Minimums)
         deallocate (Maximums)
         return
      end if

      allocate (Values(nComponents, nValues))
      if (nActive .GT. 0) then
         allocate (Active(nActive))
      end if

      WRITE (FileUnit, *) 'Number Timesteps: ', nTimes
      WRITE (FileUnit, *) 'Number Values: ', nValues
      WRITE (FileUnit, *) 'Number Components: ', nComponents
      WRITE (FileUnit, *) 'Number Active: ', nActive

      ! loop through the timesteps, read the values and active values and write
      ! them to the text file
      do iTime = 1, nTimes
         call XF_READ_VECTOR_VALUES_TIMESTEP(a_xVectorId, iTime, nValues, &
                                             nComponents, Values, error)
         if (error >= 0 .AND. nActive > 0) then
            call XF_READ_ACTIVITY_TIMESTEP(a_xVectorId, iTime, nActive, Active, error)
         end if

         ! Write the time, min, max, values and active values to the text output
         ! file.
         WRITE (FileUnit, *) ''
         WRITE (FileUnit, *) 'Timestep at  ', Times(iTime)
         WRITE (FileUnit, *) 'Min: ', Minimums(iTime)
         WRITE (FileUnit, *) 'Max: ', Maximums(iTime)

         WRITE (FileUnit, *) 'Values:'
         do i = 1, nValues
            WRITE (FileUnit, *) Values(1:nComponents, i:i)
         end do
         WRITE (FileUnit, *) ''

         WRITE (FileUnit, *) 'Activity:'
         WRITE (FileUnit, *) Active(1:nActive)
         WRITE (FileUnit, *) ''
         WRITE (FileUnit, *) ''

      end do

      if (allocated(Times)) then
         deallocate (Times)
      end if

      if (allocated(Minimums)) then
         deallocate (Minimums)
      end if

      if (allocated(Maximums)) then
         deallocate (Maximums)
      end if

      if (allocated(Values)) then
         deallocate (Values)
      end if

      if (allocated(Active)) then
         deallocate (Active)
      end if

      return
   END SUBROUTINE
! tdiReadVector

END MODULE TestDatasets
