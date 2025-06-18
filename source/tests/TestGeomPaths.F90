MODULE TestGeomPaths

   USE TestDatasets
   USE Xmdf
   USE ErrorDefinitions
   USE XmdfDefs

CONTAINS

   SUBROUTINE TM_READ_TEST_PATHS(Filename, OutFilename, error)
      use const_def, only: READWRITE, READONLY, OVERWRITE

      CHARACTER(LEN=*), INTENT(IN) :: Filename, OutFilename
      INTEGER, INTENT(OUT)         :: error
      INTEGER OutUnit
      INTEGER(XID)               :: xFileId, xGroupId
      INTEGER nGroups, nMaxPathLength, nDims, nPaths, nTimes
      INTEGER i, j, nStatus
      CHARACTER, ALLOCATABLE        :: Paths(:)
      REAL(DOUBLE), ALLOCATABLE    :: times(:), locs(:)
      REAL(DOUBLE) Nullval
      CHARACTER(LEN=BIG_STRING_SIZE) :: IndividualPath
      INTEGER StartLoc
      INTEGER, DIMENSION(2)        :: PathIndices
      INTEGER iTime, iPath

      ! open the XMDF file
      CALL XF_OPEN_FILE(Filename, READONLY, xFileId, nStatus)
      if (nStatus < 0) then
         WRITE (*, *) 'Unable to open XMDF file TM_READ_TEST_PATHS.'
         error = nStatus; 
         return
      end if

      ! open the Output file
      OutUnit = 79
      OPEN (UNIT=OutUnit, FILE=OutFilename, STATUS='REPLACE', ACTION='WRITE', &
            IOSTAT=error)
      if (OutUnit == 0) then
         call XF_CLOSE_FILE(xFileId, error)
         error = -1
         return
      end if

      ! find the geomotric path groups
      ! Get the number and paths of datasets in the file.
      CALL XF_GRP_PTHS_SZ_FOR_GEOM_PTHS(xFileId, nGroups, &
                                        nMaxPathLength, nStatus)
      if (nStatus >= 0 .AND. nGroups > 0) then
         allocate (Paths(nMaxPathLength*nGroups))
         CALL XF_GRP_PTHS_FOR_GEOM_PTHS(xFileId, nGroups, &
                                        nMaxPathLength, Paths, nStatus)
      end if
      if (nStatus < 0) then
         CALL XF_CLOSE_FILE(xFileId, nStatus)
         error = -1
         return
      end if
      ! Report the number and paths to individual geom paths groups in the file.
      WRITE (OutUnit, *) 'Number of geometric paths in file: ', nGroups

      do i = 1, nGroups
         StartLoc = (i - 1)*nMaxPathLength + 1
         IndividualPath = ''
         do j = 1, nMaxPathLength - 1
            IndividualPath(j:j) = Paths(StartLoc + j - 1)
         end do
         WRITE (OutUnit, *) 'Reading particles in group: ', &
            IndividualPath(1:LEN_TRIM(IndividualPath))

         CALL XF_OPEN_GROUP(xFileId, IndividualPath, xGroupId, nStatus)
         if (nStatus >= 0) then
            ! read the dimensionality of the paths
            CALL XF_GET_PATH_DIMENSIONALITY(xGroupId, nDims, nStatus)
            if (nStatus >= 0) then
               WRITE (OutUnit, *) 'Group dimensionality:', nDims
               ! read the number of paths
               CALL XF_GET_NUMBER_OF_PATHS(xGroupId, nPaths, nStatus)
               if (nStatus >= 0) then
                  WRITE (OutUnit, *) 'Number of paths in group:', nPaths
                  ! read the number of timesteps
                  CALL XF_GET_NUMBER_OF_TIMES(xGroupId, nTimes, nStatus)
                  if (nStatus >= 0) then
                     WRITE (OutUnit, *) 'Number of timesteps in group:', nTimes
                     ! allocate times array
                     allocate (times(nTimes))
                     CALL XF_GET_PATH_TIMES_ARRAY(xGroupId, nTimes, times, nStatus)
                     if (nStatus >= 0) then
                        ! get space for the data location values
                        allocate (locs(nPaths*nDims))
                        ! read null value
                        CALL XF_GET_PATH_NULL_VAL(xGroupId, NullVal, nStatus)
                        if (nStatus >= 0) then
                           WRITE (OutUnit, *) 'Null value:', NullVal
                           do iTime = 1, nTimes
                              WRITE (OutUnit, *) 'Timestep: ', times(iTime)
                              ! read the data for this timestep
                              CALL XF_READ_PATH_LOCATIONS_AT_TIME(xGroupId, iTime, &
                                                                  1, nPaths, locs, nStatus)
                              if (nStatus >= 0) then
                                 ! write  the data for this timestep
                                 WRITE (OutUnit, *) '  X        Y'
                                 if (nDims == 3) then
                                    WRITE (OutUnit, *) '       Z'
                                 end if
                                 WRITE (OutUnit, *) ''
                                 do j = 1, nPaths
                                    if (locs(j*nDims) == NullVal) then
                                       WRITE (OutUnit, *) 'Particle not active yet'
                                    else
                                       WRITE (OutUnit, *) locs((j - 1)*nDims + 1), ' ', locs((j - 1)*nDims + 2)
                                       if (nDims == 3) then
                                          WRITE (OutUnit, *) ' ', locs((j - 1)*nDims + 3)
                                       end if
                                       WRITE (OutUnit, *) ''
                                    end if
                                 end do ! Loop through the paths
                              end if ! if timestep read
                           end do ! Loop through timesteps
                        end if ! if null value read
                        if (allocated(locs)) deallocate (locs)
                     end if ! if we get get the times array
                     ! get space for the data location values - 1 particle
                     ! all times
                     allocate (locs(nTimes*nDims))
                     do iPath = 1, nPaths
                        ! read the values for this particle for all timesteps
                        CALL XF_READ_PATH_LOCS_FOR_PART(xGroupId, iPath, &
                                                        1, nTimes, locs, nStatus)
                        if (nStatus >= 0) then
                           ! write  the data for this path
                           WRITE (OutUnit, *) 'Time       X        Y'
                           if (nDims == 3) then
                              WRITE (OutUnit, *) '        Z'
                           end if
                           WRITE (OutUnit, *) ''
                           do j = 1, nTimes
                              if (locs((j - 1)*nDims + 1) .NE. NullVal) then
                                 WRITE (OutUnit, *) times(j), ' ', locs((j - 1)*nDims + 1), ' ', &
                                    locs((j - 1)*nDims + 2)
                                 if (nDims == 3) then
                                    WRITE (OutUnit, *) ' ', locs((j - 1)*nDims + 3)
                                 end if
                                 WRITE (OutUnit, *) ''
                              end if
                           end do ! Loop through the times
                        end if ! if read path locations for particle
                     end do ! Loop through the paths
                     if (allocated(locs)) deallocate (locs)
                  end if ! if allocation of locs suceeded
                  ! get space for the data location values - 2 particle
                  ! all times
                  allocate (locs(nTimes*nDims*2))
                  PathIndices(1) = 1
                  PathIndices(2) = nPaths
                  CALL XF_READ_PATH_LOCS_FOR_PARTS(xGroupId, 2, &
                                                   PathIndices, 1, nTimes, locs, nStatus)
                  if (nStatus >= 0) then
                     ! write  the data for these 2 paths
                     if (nDims == 3) then
                        WRITE (OutUnit, *) 'Timestep       X1        Y1        Z1        Xn        Yn        Zn'
                     else
                        WRITE (OutUnit, *) 'Timestep       X1        Y1        Xn        Yn'
                     end if
                     do j = 1, nTimes
                        if (nDims == 3) then
                           WRITE (OutUnit, *) times(j), ' ', locs((j - 1)*2*nDims + 1), ' ', &
                              locs((j - 1)*2*nDims + 2), ' ', locs((j - 1)*2*nDims + 3), ' ', &
                              locs((j - 1)*2*nDims + 4), locs((j - 1)*2*nDims + 5), ' ', locs((j - 1)*2*nDims + 6)
                        else
                           WRITE (OutUnit, *) times(j), ' ', locs((j - 1)*2*nDims + 1), ' ', &
                              locs((j - 1)*2*nDims + 2), ' ', locs((j - 1)*2*nDims + 3)
                        end if
                     end do
                  else
                     WRITE (*, *) 'Error reading path locations for multiple particles'
                  end if
               end if
               if (allocated(locs)) deallocate (locs)
               if (allocated(times)) deallocate (times)
               CALL XF_CLOSE_GROUP(xGroupId, nStatus)
               if (nStatus < 0) then
                  WRITE (*, *) 'Error reading geometric paths..'
               end if
            end if
         end if
      end do ! loop through groups
      ! free the paths
      if (allocated(Paths)) deallocate (Paths)
      ! close the files
      close (OutUnit)
      CALL XF_CLOSE_FILE(xFileId, nStatus)
      error = nStatus; 
      return

   END SUBROUTINE
! tmReadTestPaths

   SUBROUTINE TM_WRITE_TEST_PATHS(Filename, Compression, error)
      use const_def, only: READWRITE, READONLY, OVERWRITE

      CHARACTER(LEN=*), INTENT(IN) :: Filename
      INTEGER, INTENT(IN)         :: Compression
      INTEGER, INTENT(OUT)         :: error
      INTEGER nPaths
      REAL(DOUBLE), DIMENSION(6)   :: pLocs
      REAL, DIMENSION(2)           :: pSpeeds
      REAL(DOUBLE) NullVal

      REAL NullValSinglePrec
      INTEGER(XID) xFileId, xPathGroupId, xSpeedId, xPropGroupId
      INTEGER status
      REAL(DOUBLE) Timeval

      nPaths = 0
      NullVal = -99999.9d0
      NullValSinglePrec = -99999.9

      ! create the file
      call XF_CREATE_FILE(Filename, OVERWRITE, xFileId, status)
      if (status < 0) then
         error = -1
         return
      end if

      ! create the group to store the particle paths
      call XF_CREATE_GEOMETRIC_PATH_GROUP(xFileId, "particles", &
                                          'abcdefglasdfjoaieur', Compression, xPathGroupId, NullVal, status)
      if (status < 0) then
         error = -1
         return
      end if

      ! create the data set to store the speed
      call XF_CREATE_SCALAR_DSET_EXTNDBL(xPathGroupId, 'Vmag', 'm/s', &
                                         TS_SECONDS, NullValSinglePrec, Compression, xSpeedId, status)
      if (status < 0) then
         error = -1
         return
      end if

      call XF_CREATE_PROPERTY_GROUP(xSpeedId, xPropGroupId, status)
      if (status < 0) then
         call XF_CLOSE_GROUP(xSpeedId, status)
         error = -1
         return
      end if
      call XF_WRITE_PROPERTY_FLOAT(xPropGroupId, PROP_NULL_VALUE, 1, &
                                   NullValSinglePrec, -1, status)
      call XF_CLOSE_GROUP(xPropGroupId, status)

      ! Setup the arrays for the path group at timestep 0
      ! particle location at the first timestep
      nPaths = 1
      pLocs(1) = 1.0
      pLocs(2) = 2.0
      pLocs(3) = 3.0

      ! store the particles for the first timestep
      Timeval = 0.0
      call XF_WRITE_PARTICLE_TIMESTEP(xPathGroupId, 3, Timeval, nPaths, pLocs, status)
      if (status < 0) then
         error = -1
         return
      end if
      ! set up and store the speed at timestep 0
      pSpeeds(1) = 1.1
      Timeval = 0.0
      call XF_WRITE_SCALAR_TIMESTEP(xSpeedId, Timeval, 1, pSpeeds, status)
      if (status < 0) then
         error = -1
         return
      end if

      ! Setup the arrays for the path group at timestep 1
      ! particle location at the first timestep
      pLocs(1) = 4.0
      pLocs(2) = 5.0
      pLocs(3) = 6.0

      ! store the particles for the second timestep
      Timeval = 1.0
      call XF_WRITE_PARTICLE_TIMESTEP(xPathGroupId, 3, Timeval, nPaths, pLocs, status)
      if (status < 0) then
         error = -1
         return
      end if

      ! set up and store the speed at timestep 2
      pSpeeds(1) = 1.2
      call XF_WRITE_SCALAR_TIMESTEP(xSpeedId, Timeval, 1, pSpeeds, status)
      if (status < 0) then
         error = -1
         return
      end if

      ! Setup the arrays for the path group at timestep 3-add a particle
      ! particle location at the first timestep
      nPaths = 2
      pLocs(1) = 7.0
      pLocs(2) = 8.0
      pLocs(3) = 9.0
      pLocs(4) = -1.0
      pLocs(5) = -2.0
      pLocs(6) = -3.0

      ! store the particles for the timestep 3
      Timeval = 2.0
      call XF_WRITE_PARTICLE_TIMESTEP(xPathGroupId, 3, Timeval, nPaths, pLocs, status)
      if (status < 0) then
         error = -1
         return
      end if

      ! extend the data set for speed
      call XF_EXTEND_SCALAR_DATASET(xSpeedId, 2, status)
      if (status < 0) then
         error = -1
         return
      end if

      ! set up and store the speed at timestep 4
      pSpeeds(1) = 1.3
      pSpeeds(2) = 2.1
      call XF_WRITE_SCALAR_TIMESTEP(xSpeedId, Timeval, 2, pSpeeds, status)
      if (status < 0) then
         error = -1
         return
      end if

      ! Setup the arrays for the path group at timestep 3-inactive particle(static)
      ! particle location at the first timestep
      pLocs(1) = 7.0
      pLocs(2) = 8.0
      pLocs(3) = 9.0
      pLocs(4) = -4.0
      pLocs(5) = -5.0
      pLocs(6) = -6.0

      ! store the particles for timestep 4
      Timeval = 3.0
      call XF_WRITE_PARTICLE_TIMESTEP(xPathGroupId, 3, Timeval, nPaths, pLocs, status)
      if (status < 0) then
         error = -1
         return
      end if

      ! set up and store the speed at timestep 4
      pSpeeds(1) = NullVal
      pSpeeds(2) = 2.2
      call XF_WRITE_SCALAR_TIMESTEP(xSpeedId, Timeval, 2, pSpeeds, status)
      if (status < 0) then
         error = -1
         return
      end if

      ! close the resources
      call XF_CLOSE_GROUP(xSpeedId, status)
      call XF_CLOSE_GROUP(xPathGroupId, status)
      call XF_CLOSE_FILE(xFileId, status)

      error = 1
      return
   END SUBROUTINE

END MODULE TestGeomPaths
