MODULE TestMesh

   USE TestDatasets
   USE Xmdf
   USE ErrorDefinitions
   USE XmdfDefs

   CHARACTER(LEN=*), PARAMETER :: MESH_A_GROUP_NAME = 'MeshA Group'
   CHARACTER(LEN=*), PARAMETER :: MESH_B_GROUP_NAME = 'MeshB Group'

CONTAINS

!****************************
! ---------------------------------------------------------------------------
! FUNCTION  TM_READ_MESH
! PURPOSE   Read a mesh file and write a text output file
! NOTES
! ---------------------------------------------------------------------------
   SUBROUTINE TM_READ_MESH(xGroupId, a_OutFile, error)
      INTEGER(XID), INTENT(IN) :: xGroupId
      INTEGER, INTENT(IN)        :: a_OutFile
      INTEGER, INTENT(OUT)       :: error
      INTEGER nElems, nNodes, nNodesPerElem, nElemType, nNodeId
      LOGICAL*2 ElementsOneType
      INTEGER status, i, j
      INTEGER StrType, UIntType, IntType, DblType, FloatType
      INTEGER(XID) xPropGroupId
      INTEGER(XID), ALLOCATABLE      :: NodesInElem(:)
      INTEGER, ALLOCATABLE :: ElemTypes(:)
      REAL(DOUBLE), ALLOCATABLE :: XNodeLocs(:), YNodeLocs(:), ZNodeLocs(:)

      ! Get the number of elements, nodes, and Maximum number of nodes per element
      call XF_GET_NUMBER_OF_ELEMENTS(xGroupId, nElems, status)
      if (status >= 0) then
         call XF_GET_NUMBER_OF_NODES(xGroupId, nNodes, status)
         if (status >= 0) then
            call XF_GET_MAX_NODES_IN_ELEM(xGroupId, nNodesPerElem, status)
         end if
      end if

      if (status < 0) then
         error = -1
         return
      end if

      ! Do Element information first
      WRITE (a_OutFile, *) 'Number of Elements: ', nElems

      ! Element types
      call XF_ARE_ALL_ELEMS_SAME_TYPE(xGroupId, ElementsOneType, status)
      if (status < 0) then
         error = -1
         return
      end if

      if (ElementsOneType) then
         call XF_READ_ELEM_TYPES_SINGLE_VALUE(xGroupId, nElemType, status)
         WRITE (a_OutFile, *) 'All elements are type ', nElemType
      else
         allocate (ElemTypes(nElems))
         call XF_READ_ELEM_TYPES(xGroupId, nElems, ElemTypes, status)
         if (status < 0) then
            error = -1
            return
         end if
         WRITE (a_OutFile, *) 'Element Types:'
         do i = 1, nElems
            WRITE (a_OutFile, *) 'Elem ', i, ', ', 'Type ', ElemTypes(i)
         end do
         deallocate (ElemTypes)
      end if

      ! Nodes in each element
      allocate (NodesInElem(nElems*nNodesPerElem))
      call XF_READ_ELEM_NODE_IDS(xGroupId, nElems, nNodesPerElem, NodesInElem, error)
      if (error < 0) then
         WRITE (*, *) 'Error reading mesh'
         return
      end if

      do i = 1, nElems
         WRITE (a_OutFile, *) 'Elem: ', i, ' - '
         do j = 1, nNodesPerElem
            nNodeId = NodesInElem((i - 1)*nNodesPerElem + j)
            if (nNodeId > 0) then  ! -1 is for unused array locations
               WRITE (a_OutFile, *) nNodeId, ' '
            end if
         end do
         WRITE (a_OutFile, *) ''
      end do

      if (allocated(NodesInElem)) deallocate (NodesInElem)

      ! NodeLocations
      allocate (XNodeLocs(nNodes))
      allocate (YNodeLocs(nNodes))
      allocate (ZNodeLocs(nNodes))

      call XF_READ_X_NODE_LOCATIONS(xGroupId, nNodes, XNodeLocs, status)
      if (status >= 0) then
         call XF_READ_Y_NODE_LOCATIONS(xGroupId, nNodes, YNodeLocs, status)
         if (status >= 0) then
            call XF_READ_Z_NODE_LOCATIONS(xGroupId, nNodes, ZNodeLocs, status)
         end if
      end if

      WRITE (a_OutFile, *) 'Node Locations:'
      do i = 1, nNodes
         WRITE (a_OutFile, *) 'Node: ', i, ' Location: ', XNodeLocs(i), ' ', &
            YNodeLocs(i), ' ', ZNodeLocs(i)
      end do

      deallocate (XNodeLocs)
      deallocate (YNodeLocs)
      deallocate (ZNodeLocs)

      ! Open Properties Group
      call XF_OPEN_GROUP(xGroupId, 'PROPERTIES', xPropGroupId, status)
      if (status < 0) then
         WRITE (a_OutFile, *) ''
         WRITE (a_OutFile, *) 'Properties Group not found'
         WRITE (a_OutFile, *) ''
         error = -1
         return
      end if

      call XF_GET_PROPERTY_TYPE(xPropGroupId, 'String', StrType, status)
      call XF_GET_PROPERTY_TYPE(xPropGroupId, 'UInt', UIntType, status)
      call XF_GET_PROPERTY_TYPE(xPropGroupId, 'Int', IntType, status)
      call XF_GET_PROPERTY_TYPE(xPropGroupId, 'Double', DblType, status)
      call XF_GET_PROPERTY_TYPE(xPropGroupId, 'Float', FloatType, status)

      ! Property Types:
      WRITE (a_OutFile, *) ''
      if (StrType == XF_TYPE_STRING) then
         WRITE (a_OutFile, *) 'String Property Type Read Correctly'
      else
         WRITE (a_OutFile, *) 'Error in Getting String Property Type'
      end if
      if (UIntType == XF_TYPE_UINT) then
         WRITE (a_OutFile, *) 'Unsigned Integer Property Type Read Correctly'
      else
         WRITE (a_OutFile, *) 'Error in Getting Unsigned Integer Property Type'
      end if
      if (IntType == XF_TYPE_INT) then
         WRITE (a_OutFile, *) 'Integer Property Type Read Correctly'
      else
         WRITE (a_OutFile, *) 'Error in Getting Integer Property Type'
      end if
      if (DblType == XF_TYPE_DOUBLE) then
         WRITE (a_OutFile, *) 'Double Property Type Read Correctly'
      else
         WRITE (a_OutFile, *) 'Error in Getting Double Property Type'
      end if
      if (FloatType == XF_TYPE_FLOAT) then
         WRITE (a_OutFile, *) 'Float Property Type Read Correctly'
      else
         WRITE (a_OutFile, *) 'Error in Getting Float Property Type'
      end if
      WRITE (a_OutFile, *) ''

      error = 0
      return

   END SUBROUTINE

!****************************
!------------------------------------------------------------------------------
! FUNCTION   TM_WRITE_TEST_MESH_A
! PURPOSE    Write a file that contains data for an all tri mesh
! NOTES      A picture of the mesh is in the file (TestMeshA.gif)
!            error equals TRUE on success and FALSE on failure
!------------------------------------------------------------------------------
   SUBROUTINE TM_WRITE_TEST_MESH_A(Filename, Compression, error)
      use const_def, only: READWRITE, READONLY, OVERWRITE

      CHARACTER(LEN=*), INTENT(IN) :: Filename
      INTEGER, INTENT(IN)          :: Compression
      INTEGER, INTENT(OUT)         :: error
      INTEGER nElements, nNodes
      INTEGER(XID) xFileId, xMeshId, xPropGroupId, xCoordId
      REAL(DOUBLE), DIMENSION(5) :: dNodeLocsX
      REAL(DOUBLE), DIMENSION(5) :: dNodeLocsY
      REAL(DOUBLE), DIMENSION(5) :: dNodeLocsZ
      INTEGER(XID), DIMENSION(3, 3)    :: iElementNodes
      INTEGER status, ElemType, propint(1), iEllipse
      CHARACTER(LEN=BIG_STRING_SIZE) propstring
      INTEGER propuint(1)
      REAL(DOUBLE) propdouble(1), dMajorR, dMinorR
      REAL propfloat

      ! set Element type to EL_TYPE_TRI_LINEAR
      ElemType = EL_TYPE_TRI_LINEAR

      nElements = 3
      nNodes = 5
      xFileId = NONE
      xMeshId = NONE

      ! Setup the arrays for the mesh data
      ! nodes
      dNodeLocsX(1) = 0.0
      dNodeLocsX(2) = 5.0
      dNodeLocsX(3) = 0.0
      dNodeLocsX(4) = 5.0
      dNodeLocsX(5) = 7.5

      dNodeLocsY(1) = 5.0
      dNodeLocsY(2) = 5.0
      dNodeLocsY(3) = 0.0
      dNodeLocsY(4) = 0.0
      dNodeLocsY(5) = 2.5

      dNodeLocsZ(1) = 0.0
      dNodeLocsZ(2) = 0.0
      dNodeLocsZ(3) = 0.0
      dNodeLocsZ(4) = 0.0
      dNodeLocsZ(5) = 0.0

      ! nodes for each element must be counter-clockwize
      iElementNodes(1, 1) = 1
      iElementNodes(2, 1) = 3
      iElementNodes(3, 1) = 2

      iElementNodes(1, 2) = 2
      iElementNodes(2, 2) = 3
      iElementNodes(3, 2) = 4

      iElementNodes(1, 3) = 5
      iElementNodes(2, 3) = 2
      iElementNodes(3, 3) = 4

      ! create the file
      call XF_CREATE_FILE(Filename, OVERWRITE, xFileId, status)
      if (status < 0) then
         ! close the resources
         call XF_CLOSE_FILE(xFileId, error)
         error = -1
         return
      end if

      ! create the group to store the mesh
      call XF_CREATE_GROUP_FOR_MESH(xFileId, MESH_A_GROUP_NAME, xMeshId, status)
      if (status < 0) then
         ! close the resources
         call XF_CLOSE_GROUP(xMeshId, error)
         call XF_CLOSE_FILE(xFileId, error)
         error = -1
         return
      end if

      ! Element types - all are linear triangles
      call XF_SET_ALL_ELEMS_SAME_TYPE(xMeshId, ElemType, status)
      if (status < 0) then
         ! close the resources
         call XF_CLOSE_GROUP(xMeshId, error)
         call XF_CLOSE_FILE(xFileId, error)
         error = -1
         return
      end if

      ! node information
      call XF_SET_NUMBER_OF_NODES(xMeshId, nNodes, status)
      if (status < 0) then
         ! close the resources
         call XF_CLOSE_GROUP(xMeshId, error)
         call XF_CLOSE_FILE(xFileId, error)
         error = -1
         return
      end if

      call XF_WRITE_X_NODE_LOCATIONS(xMeshId, nNodes, dNodeLocsX, Compression, status)
      if (status < 0) then
         ! close the resources
         call XF_CLOSE_GROUP(xMeshId, error)
         call XF_CLOSE_FILE(xFileId, error)
         error = -1
         return
      end if

      call XF_WRITE_Y_NODE_LOCATIONS(xMeshId, nNodes, dNodeLocsY, status)
      if (status < 0) then
         ! close the resources
         call XF_CLOSE_GROUP(xMeshId, error)
         call XF_CLOSE_FILE(xFileId, error)
         error = -1
         return
      end if

      call XF_WRITE_Z_NODE_LOCATIONS(xMeshId, nNodes, dNodeLocsZ, status)
      if (status < 0) then
         ! close the resources
         call XF_CLOSE_GROUP(xMeshId, error)
         call XF_CLOSE_FILE(xFileId, error)
         error = -1
         return
      end if

      ! element information
      call XF_SET_NUMBER_OF_ELEMENTS(xMeshId, nElements, status)
      if (status < 0) then
         ! close the resources
         call XF_CLOSE_GROUP(xMeshId, error)
         call XF_CLOSE_FILE(xFileId, error)
         error = -1
         return
      end if

      ! Write the node array ids for the nodes in each element
      call XF_WRITE_ELEM_NODE_IDS(xMeshId, nElements, 3, iElementNodes, &
                                  Compression, status)
      if (status < 0) then
         ! close the resources
         call XF_CLOSE_GROUP(xMeshId, error)
         call XF_CLOSE_FILE(xFileId, error)
         error = -1
         return
      end if

      ! Write the Property File
      call XF_CREATE_MESH_PROPERTY_GROUP(xMeshId, xPropGroupId, status)
      if (status < 0) then
         call XF_CLOSE_GROUP(xMeshId, error)
         call XF_CLOSE_FILE(xFileId, error)
         error = -1
         return
      end if

      propstring = 'Property String'
      propuint = 5
      propint = -5
      propdouble = 5.6789012345d0
      propfloat = 5.6789

      call XF_WRITE_PROPERTY_STRING(xPropGroupId, 'String', propstring, status)
      call XF_WRITE_PROPERTY_UINT(xPropGroupId, 'UInt', 1, propuint, NONE, status)
      call XF_WRITE_PROPERTY_INT(xPropGroupId, 'Int', 1, propint, NONE, status)
      call XF_WRITE_PROPERTY_DOUBLE(xPropGroupId, 'Double', 1, &
                                    propdouble, NONE, status)
      call XF_WRITE_PROPERTY_FLOAT(xPropGroupId, 'Float', 1, &
                                   propfloat, NONE, status)

      ! Write Coordinate file - for MeshA, we will set the coordinate system to be
      !   Geogrpahic, with Latitude, Longitude, and user-defined ellipsoid settings
      !   written to the file.
      call XF_CREATE_COORDINATE_GROUP(xFileId, xCoordId, status)
      if (status < 0) then
         call XF_CLOSE_GROUP(xPropGroupId, error)
         call XF_CLOSE_GROUP(xMeshId, error)
         call XF_CLOSE_FILE(xFileId, error)
         error = status
         return
      end if

      ! set coordinate values
      iEllipse = 32    ! User defined
      dMajorR = 45.0   ! Made up
      dMinorR = 32.0   ! Made up

      call XF_SET_HORIZ_DATUM(xCoordId, HORIZ_DATUM_GEOGRAPHIC, error)
      call XF_SET_HORIZ_UNITS(xCoordId, COORD_UNITS_US_FEET, error)

      call XF_SET_VERT_DATUM(xCoordId, VERT_DATUM_LOCAL, error)
      call XF_SET_VERT_UNITS(xCoordId, COORD_UNITS_US_FEET, error)

      ! write additional information
      call XF_SET_ELLIPSE(xCoordId, iEllipse, error)
      call XF_SET_LAT(xCoordId, LATITUDE_NORTH, error)
      call XF_SET_LON(xCoordId, LONGITUDE_EAST, error)
      call XF_SET_MAJOR_R(xCoordId, dMajorR, error)
      call XF_SET_MINOR_R(xCoordId, dMinorR, error)

      call XF_CLOSE_GROUP(xCoordId, error)
      xCoordId = 0

      ! close the resources
      call XF_CLOSE_GROUP(xPropGroupId, error)
      call XF_CLOSE_GROUP(xMeshId, error)
      call XF_CLOSE_FILE(xFileId, error)

      return

   END SUBROUTINE

!****************************
!------------------------------------------------------------------------------
! FUNCTION   TM_WRITE_TEST_MESH_B
! PURPOSE    Write a file that contains data for an mixed quad/tri linear mesh
! NOTES      A picture of the mesh is in the file (TestMeshB.gif)
!            error equals TRUE on success and FALSE on failure
!------------------------------------------------------------------------------
   SUBROUTINE TM_WRITE_TEST_MESH_B(Filename, Compression, error)
      use const_def, only: READWRITE, READONLY, OVERWRITE

      CHARACTER(LEN=*), INTENT(IN) :: Filename
      INTEGER, INTENT(IN)          :: Compression
      INTEGER, INTENT(OUT)         :: error
      INTEGER nElements, nNodes, nMaxNodePerElem
      INTEGER(XID) xFileId, xMeshId, xPropGroupId, xCoordId
      REAL(DOUBLE), DIMENSION(5)   :: dNodeLocsX
      REAL(DOUBLE), DIMENSION(5)   :: dNodeLocsY
      REAL(DOUBLE), DIMENSION(5)   :: dNodeLocsZ
      INTEGER(XID), DIMENSION(4, 2)      :: iElementNodes
      INTEGER, DIMENSION(2)        :: iElementTypes
      INTEGER status, propint(1), iEllipse
      CHARACTER(LEN=BIG_STRING_SIZE) propstring
      INTEGER propuint(1)
      REAL(DOUBLE) propdouble(1)
      REAL propfloat

      nElements = 2
      nNodes = 5
      nMaxNodePerElem = 4
      xFileId = NONE
      xMeshId = NONE

      ! Setup the arrays for the mesh data
      ! nodes
      dNodeLocsX(1) = 0.0
      dNodeLocsX(2) = 5.0
      dNodeLocsX(3) = 0.0
      dNodeLocsX(4) = 5.0
      dNodeLocsX(5) = 7.5

      dNodeLocsY(1) = 5.0
      dNodeLocsY(2) = 5.0
      dNodeLocsY(3) = 0.0
      dNodeLocsY(4) = 0.0
      dNodeLocsY(5) = 2.5

      dNodeLocsZ(1) = 0.0
      dNodeLocsZ(2) = 0.0
      dNodeLocsZ(3) = 0.0
      dNodeLocsZ(4) = 0.0
      dNodeLocsZ(5) = 0.0

      ! nodes for each element must be counter-clockwize
      iElementNodes(1, 1) = 1
      iElementNodes(2, 1) = 3
      iElementNodes(3, 1) = 4
      iElementNodes(4, 1) = 2

      iElementNodes(1, 2) = 2
      iElementNodes(2, 2) = 4
      iElementNodes(3, 2) = 5
      iElementNodes(4, 2) = NONE; 
      iElementTypes(1) = EL_TYPE_QUADRILATERAL_LINEAR; 
      iElementTypes(2) = EL_TYPE_TRI_LINEAR; 
      ! create the file
      call XF_CREATE_FILE(Filename, OVERWRITE, xFileId, status)
      if (status < 0) then
         ! close the resources
         call XF_CLOSE_FILE(xFileId, error)
         error = -1
         return
      end if

      ! create the group to store the mesh
      call XF_CREATE_GROUP_FOR_MESH(xFileId, MESH_B_GROUP_NAME, xMeshId, status)
      if (status < 0) then
         ! close the resources
         call XF_CLOSE_GROUP(xMeshId, error)
         call XF_CLOSE_FILE(xFileId, error)
         error = -1
         return
      end if

      ! node information
      call XF_SET_NUMBER_OF_NODES(xMeshId, nNodes, status)
      if (status < 0) then
         ! close the resources
         call XF_CLOSE_GROUP(xMeshId, error)
         call XF_CLOSE_FILE(xFileId, error)
         error = -1
         return
      end if

      call XF_WRITE_X_NODE_LOCATIONS(xMeshId, nNodes, dNodeLocsX, &
                                     Compression, status)
      if (status < 0) then
         ! close the resources
         call XF_CLOSE_GROUP(xMeshId, error)
         call XF_CLOSE_FILE(xFileId, error)
         error = -1
         return
      end if

      call XF_WRITE_Y_NODE_LOCATIONS(xMeshId, nNodes, dNodeLocsY, status)
      if (status < 0) then
         ! close the resources
         call XF_CLOSE_GROUP(xMeshId, error)
         call XF_CLOSE_FILE(xFileId, error)
         error = -1
         return
      end if

      call XF_WRITE_Z_NODE_LOCATIONS(xMeshId, nNodes, dNodeLocsZ, status)
      if (status < 0) then
         ! close the resources
         call XF_CLOSE_GROUP(xMeshId, error)
         call XF_CLOSE_FILE(xFileId, error)
         error = -1
         return
      end if

      ! element information
      call XF_SET_NUMBER_OF_ELEMENTS(xMeshId, nElements, status)
      if (status < 0) then
         ! close the resources
         call XF_CLOSE_GROUP(xMeshId, error)
         call XF_CLOSE_FILE(xFileId, error)
         error = -1
         return
      end if

      ! Element types
      call XF_WRITE_ELEM_TYPES(xMeshId, nElements, iElementTypes, &
                               Compression, status)
      if (status < 0) then
         ! close the resources
         call XF_CLOSE_GROUP(xMeshId, error)
         call XF_CLOSE_FILE(xFileId, error)
         error = -1
         return
      end if

      ! Write the node array ids for the nodes in each element
      call XF_WRITE_ELEM_NODE_IDS(xMeshId, nElements, nMaxNodePerElem, &
                                  iElementNodes, Compression, status)
      if (status < 0) then
         ! close the resources
         call XF_CLOSE_GROUP(xMeshId, error)
         call XF_CLOSE_FILE(xFileId, error)
         error = -1
         return
      end if

      ! Write the Property File
      call XF_CREATE_MESH_PROPERTY_GROUP(xMeshId, xPropGroupId, status)
      if (status < 0) then
         call XF_CLOSE_GROUP(xMeshId, error)
         call XF_CLOSE_GROUP(xFileId, error)
         error = -1
         return
      end if

      propstring = 'String Property'
      propuint = 2
      propint = -2
      propdouble = 2.3456789012d0
      propfloat = 2.3456

      call XF_WRITE_PROPERTY_STRING(xPropGroupId, 'String', propstring, status)
      call XF_WRITE_PROPERTY_UINT(xPropGroupId, 'UInt', 1, propuint, NONE, status)
      call XF_WRITE_PROPERTY_INT(xPropGroupId, 'Int', 1, propint, NONE, status)
      call XF_WRITE_PROPERTY_DOUBLE(xPropGroupId, 'Double', 1, &
                                    propdouble, NONE, status)
      call XF_WRITE_PROPERTY_FLOAT(xPropGroupId, 'Float', 1, &
                                   propfloat, NONE, status)

      ! Write Coordinate file - for MeshB, we will set the coordinate system to be
      !   Geogrpahic, with Latitude, Longitude, and standard ellipsoid settings
      !   written to the file.
      call XF_CREATE_COORDINATE_GROUP(xFileId, xCoordId, status)
      if (status < 0) then
         call XF_CLOSE_GROUP(xPropGroupId, error)
         call XF_CLOSE_GROUP(xMeshId, error)
         call XF_CLOSE_FILE(xFileId, error)
         error = status
         return
      end if

      ! set coordinate values
      iEllipse = 21   ! International 1924

      call XF_SET_HORIZ_DATUM(xCoordId, HORIZ_DATUM_GEOGRAPHIC, error)
      call XF_SET_HORIZ_UNITS(xCoordId, COORD_UNITS_METERS, error)

      call XF_SET_VERT_DATUM(xCoordId, VERT_DATUM_NGVD_88, error)
      call XF_SET_VERT_UNITS(xCoordId, COORD_UNITS_METERS, error)

      ! write additional information
      call XF_SET_ELLIPSE(xCoordId, iEllipse, error)
      call XF_SET_LAT(xCoordId, LATITUDE_SOUTH, error)
      call XF_SET_LON(xCoordId, LONGITUDE_WEST, error)

      call XF_CLOSE_GROUP(xCoordId, error)
      xCoordId = 0

      ! close the resources
      call XF_CLOSE_GROUP(xPropGroupId, error)
      call XF_CLOSE_GROUP(xMeshId, error)
      call XF_CLOSE_FILE(xFileId, error)

      return

   END SUBROUTINE

END MODULE TestMesh
