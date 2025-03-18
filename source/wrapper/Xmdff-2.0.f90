! XMDF License & Copyright Notice Copyright Notice and Statement for the
! eXtensible Model Data Format (XMDF) Software Library and API.  XMDF
! Software Library and Utilities Copyright 2004, 2005 by the
! Enivironmental Modeling Research Laboratory and Brigham Young
! University.  All rights reserved.  Redistribution and use in source
! and binary forms, with or without modification, are permitted for any
! purpose (including commercial purposes) provided that the following
! conditions are met: 1.  Redistributions of source code must retain the
! above copyright notice, this list of conditions, and the following
! disclaimer.  2.  Redistributions in binary form must reproduce the
! above copyright notice, this list of conditions, and the following
! disclaimer in the documentation and/or materials provided with the
! distribution.  3.  In addition, redistributions of modified forms of
! the source or binary code must carry prominent notices stating that
! the original code was changed and the date of the change.  4.  All
! publications or advertising materials mentioning features or use of
! this software are asked, but not required, to acknowledge that it was
! developed by the EMRL at Brigham Young University in partnership with
! the US Army Engineer Research and Development Center and to credit the
! contributors.  5.  Neither the name of the University nor the names of
! the Contributors maybe used to endorse or promote products derived
! from this software without specific prior written permission from the
! University or the Contributors, as appropriate for the name(s) to be
! used.  6.  THIS SOFTWARE IS PROVIDED BY THE UNIVERSITY AND THE
! CONTRIBUTORS "AS IS"WITH NO WARRANTY OF ANY KIND, EITHER EXPRESSED OR
! IMPLIED.  In no event shall the University or the Contributors be
! liable for any damages suffered by the users arising out of the use of
! this software, even if advised of the possibility of such damage.
! --------------------------------------------------------------------------
! Portions of XMDF utilize the HDF5 software library.  The following
! statement applies to HDF5: Copyright Notice and Statement for NCSA
! Hierarchical Data Format (HDF) Software Library and Utilities NCSA
! HDF5 (Hierarchical Data Format 5) Software Library and Utilities
! Copyright 1998, 1999, 2000, 2001, 2002, 2003, 2004, 2005 by the Board
! of Trustees of the University of Illinois.  All rights reserved.
! Contributors: National Center for Supercomputing Applications (NCSA)
! at the University of Illinois at Urbana-Champaign (UIUC), Lawrence
! Livermore National Laboratory (LLNL), Sandia National Laboratories
! (SNL), Los Alamos National Laboratory (LANL), Jean-loup Gailly and
! Mark Adler (gzip library).  Redistribution and use in source and
! binary forms, with or without modification, are permitted for any
! purpose (including commercial purposes) provided that the following
! conditions are met: 1.  Redistributions of source code must retain the
! above copyright notice, this list of conditions, and the following
! disclaimer.  2.  Redistributions in binary form must reproduce the
! above copyright notice, this list of conditions, and the following
! disclaimer in the documentation and/or materials provided with the
! distribution.  3.  In addition, redistributions of modified forms of
! the source or binary code must carry prominent notices stating that
! the original code was changed and the date of the change.  4.  All
! publications or advertising materials mentioning features or use of
! this software are asked, but not required, to acknowledge that it was
! developed by the National Center for Supercomputing Applications at
! the University of Illinois at Urbana-Champaign and to credit the
! contributors.  5.  Neither the name of the University nor the names of
! the Contributors may be used to endorse or promote products derived
! from this software without specific prior written permission from the
! University or the Contributors, as appropriate for the name(s) to be
! used.  6.  THIS SOFTWARE IS PROVIDED BY THE UNIVERSITY AND THE
! CONTRIBUTORS "AS IS" WITH NO WARRANTY OF ANY KIND, EITHER EXPRESSED OR
! IMPLIED.  In no event shall the University or the Contributors be
! liable for any damages suffered by the users arising out of the use of
! this software, even if advised of the possibility of such damage.
! --------------------------------------------------------------------------
! Portions of HDF5 were developed with support from the University of
! California, Lawrence Livermore National Laboratory (UC LLNL).  The
! following statement applies to those portions of the product and must
! be retained in any redistribution of source code, binaries,
! documentation, and/or accompanying materials: This work was partially
! produced at the University of California, Lawrence Livermore National
! Laboratory (UC LLNL) under contract no. W-7405-ENG-48 (Contract 48)
! between the U.S. Department of Energy (DOE) and The Regents of the
! University of California (University) for the operation of UC LLNL.
! DISCLAIMER: This work was prepared as an account of work sponsored by
! an agency of the United States Government.  Neither the United States
! Government nor the University of California nor any of their
! employees, makes any warranty, express or implied, or assumes any
! liability or responsibility for the accuracy, completeness, or
! usefulness of any information, apparatus, product, or process
! disclosed, or represents that its use would not infringe privately-
! owned rights.  Reference herein to any specific commercial products,
! process, or service by trade name, trademark, manufacturer, or
! otherwise, does not necessarily constitute or imply its endorsement,
! recommendation, or favoring by the United States Government or the
! University of California.  The views and opinions of authors expressed
! herein do not necessarily state or reflect those of the United States
! Government or the University of California, and shall not be used for
! advertising or product endorsement purposes.
! --------------------------------------------------------------------------
! This work was partially produced at the University of California,
! Lawrence Livermore National Laboratory (UC LLNL) under contract
! no.W-7405-ENG-48 (Contract 48) between the U.S. Department of Energy
! (DOE) and The Regents of the University of California (University) for
! the operation of UC LLNL.  DISCLAIMER:This work was prepared as an
! account of work sponsored by an agency of the United States
! Government.  Neither the United States Government nor the University
! of California nor any of their employees, makes any warranty, express
! or implied, or assumes any liability or responsibility for the
! accuracy, completeness, or usefulness of any information, apparatus,
! product, or process disclosed, or represents that its use would not
! infringe privately-owned rights.  Reference herein to any specific
! commercial products, process, or service by trade name, trademark,
! manufacturer, or otherwise, does not necessarily constitute or imply
! its endorsement, recommendation, or favoring by the United States
! Government or the University of California.  The views and opinions of
! authors expressed herein do not necessarily state or reflect those of
! the United States Government or the University of California, and
! shall not be used for advertising or product endorsement purposes.
! --------------------------------------------------------------------------


MODULE ERRORDEFINITIONS
  !/ This file contains the error codes for errors thrown only from XMDF
  !/ If the error came from HDF5, you have to use their error routines

  !/ File errors -40xx

  INTEGER,PARAMETER :: ERROR_FILE_NOT_HDF5   =-4001 
  INTEGER,PARAMETER :: ERROR_FILE_NOT_XMDF=-4002 

  !/ Attribute errors -41xx

  INTEGER,PARAMETER :: ERROR_ATTRIBUTE_NOT_SUPPORTED=-4101

  !/ Datatype errors -42xx

  INTEGER,PARAMETER :: ERROR_INCORRECT_DATATYPE=-4201

  !/ Dataset errors -43xx

  INTEGER,PARAMETER :: ERROR_DATASET_SIZE_INCORRECT=-4301
  INTEGER,PARAMETER :: ERROR_DATASET_NO_DATA       =-4302
  INTEGER,PARAMETER :: ERROR_DATASET_DOES_NOT_EXIST=-4303
  INTEGER,PARAMETER :: ERROR_DATASET_INVALID       =-4304

  !/ Group errors -44xx

  INTEGER,PARAMETER :: ERROR_GROUP_TYPE_INCONSISTENT=-4401

  !/ Mesh errors -45xx
  INTEGER,PARAMETER :: ERROR_ELEMENT_NUM_INCORRECT= -4501   !inconsistent element number
  INTEGER,PARAMETER :: ERROR_NODE_NUM_INCORRECT   = -4502
  INTEGER,PARAMETER :: ERROR_NOT_MESH_GROUP       = -4503
  INTEGER,PARAMETER :: ERROR_MESH_INCOMPLETE      = -4504
  INTEGER,PARAMETER :: ERROR_MESH_INVALID         = -4505

  ! Grid errors
  INTEGER,PARAMETER :: ERROR_GRID_TYPE_INVALID           = -4601
  INTEGER,PARAMETER :: ERROR_GRID_NUM_DIMS = -4602
  INTEGER,PARAMETER :: ERROR_GRID_EXTRUDE_TYPE_INVALID   = -4603
  INTEGER,PARAMETER :: ERROR_GRID_NUMVALS_INCORRECT      = -4604

  ! Errors that I don't want to name 

  INTEGER,PARAMETER :: ERROR_OTHER                = -9901

    ! Guid names and lengths
  CHARACTER(LEN=*),PARAMETER :: XF_GUID         = 'Guid'
  INTEGER,PARAMETER :: XF_GUID_STRINGLENGTH     = 37

  ! Overwrite options in the function xfSetupToWriteDatasets
  INTEGER,PARAMETER :: XF_OVERWRITE_CLEAR_FILE          = 1
  ! XF_OVERWRITE_CLEAR_DATASET_GROUP
  INTEGER,PARAMETER :: XF_OVERWRITE_CLEAR_DATASET_GRP   = 2
  INTEGER,PARAMETER :: XF_OVERWRITE_NONE                = 3



END MODULE ERRORDEFINITIONS

!****************************************

MODULE XMDFDEFS

  !/ type defines

  !/ XMDF version number
  !/ Single precision float
  REAL*4, PARAMETER :: XMDF_VERSION = 1.20
  CHARACTER(LEN=*),PARAMETER :: FILE_VERSION= 'File Version'
  LOGICAL*2 :: READONLY=.true., READWRITE=.false.

  INTEGER, PARAMETER :: NULL = 0
  INTEGER, PARAMETER :: NONE   = -1

  ! When we write out data we want to write it out the same type as the PC for
  ! faster reading during post-processing
  INTEGER      XF_FLOAT_WRITE_TYPE
  INTEGER      XF_DOUBLE_WRITE_TYPE
  LOGICAL             FLOAT_TYPES_INITIALIZED
  !INTEGER             ,PARAMETER   :: XTYPE_WRITE_FLOAT = H5T_IEEE_F32LE
  !INTEGER             ,PARAMETER   :: XTYPE_WRITE_DOUBLE = H5T_IEEE_F64LE

  CHARACTER(LEN=*),PARAMETER :: FILE_TYPE_ID = 'File Type'
  CHARACTER(LEN=*),PARAMETER :: FILE_TYPE_XMDF='Xmdf'

  CHARACTER(LEN=*),PARAMETER :: FILE_CONTENTS='Contents'
  CHARACTER(LEN=*),PARAMETER :: FILE_CONTENTS_DATASETS='Datasets'

  CHARACTER(LEN=*),PARAMETER :: GROUP_TYPE='Grouptype'
  CHARACTER(LEN=*),PARAMETER :: GROUP_TYPE_GENERIC='Generic'
  CHARACTER(LEN=*),PARAMETER :: GROUP_TYPE_MESH='MESH'
  CHARACTER(LEN=*),PARAMETER :: GROUP_TYPE_MESH_ELEMS='ElemGroup'
  CHARACTER(LEN=*),PARAMETER :: GROUP_TYPE_MESH_NODES='NodeGroup'
  CHARACTER(LEN=*),PARAMETER :: GROUP_TYPE_GRID='GRID'
  CHARACTER(LEN=*),PARAMETER :: GROUP_TYPE_XSECS='XSECS'
  CHARACTER(LEN=*),PARAMETER :: GROUP_TYPE_DATASETS='DATASETS'
  CHARACTER(LEN=*),PARAMETER :: GROUP_TYPE_PROPERTIES='PROPERTIES'
  CHARACTER(LEN=*),PARAMETER :: GROUP_TYPE_MULTI_DATASETS='MULTI DATASETS'
  CHARACTER(LEN=*),PARAMETER :: GROUP_TYPE_DATASET_SCALAR='DATASET SCALAR'
  CHARACTER(LEN=*),PARAMETER :: GROUP_TYPE_DATASET_VECTOR='DATASET VECTOR'
  CHARACTER(LEN=*),PARAMETER :: GROUP_TYPE_COORDS='Coordinates'
  CHARACTER(LEN=*),PARAMETER :: GROUP_TYPE_GEOMETRIC_PATH='GEOMETRIC PATHS'

  CHARACTER(LEN=*),PARAMETER :: TIME_UNITS='TimeUnits'
  CHARACTER(LEN=*),PARAMETER :: TS_SECONDS ='Seconds'
  CHARACTER(LEN=*),PARAMETER :: TS_MINUTES ='Minutes'
  CHARACTER(LEN=*),PARAMETER :: TS_HOURS ='Hours'
  CHARACTER(LEN=*),PARAMETER :: TS_DAYS ='Days'
  CHARACTER(LEN=*),PARAMETER :: TS_NOT_APPLICABLE = 'None'
  INTEGER,PARAMETER :: TIME_UNITS_MAXLENGTH=500

  INTEGER,PARAMETER :: BIG_STRING_SIZE=1000
  INTEGER,PARAMETER :: DEFAULT_GROUP_SIZE=10
  INTEGER,PARAMETER :: MAX_ID_NAME=512

  ! Property Types
  INTEGER,PARAMETER :: XF_TYPE_INT    = 1
  INTEGER,PARAMETER :: XF_TYPE_FLOAT  = 2
  INTEGER,PARAMETER :: XF_TYPE_DOUBLE = 3
  INTEGER,PARAMETER :: XF_TYPE_STRING = 4
  INTEGER,PARAMETER :: XF_TYPE_UINT   = 5
  INTEGER,PARAMETER :: XF_TYPE_OTHER  = 11

  ! Dataset Arrays
  CHARACTER(LEN=*),PARAMETER :: DATASET_DSET_VALUES ='Values'
  CHARACTER(LEN=*),PARAMETER :: DATASET_DSET_TIMES  ='Times'
  CHARACTER(LEN=*),PARAMETER :: DATASET_DSET_MINS   ='Mins'
  CHARACTER(LEN=*),PARAMETER :: DATASET_DSET_MAXS   ='Maxs'
  CHARACTER(LEN=*),PARAMETER :: DATASET_REFTIME      ='Reftime'

  CHARACTER(LEN=*),PARAMETER :: DATASET_DSET_ACTIVE ='Active'
  CHARACTER(LEN=*),PARAMETER :: DATASET_DSET_REFTIME='Reftime'

  INTEGER,PARAMETER :: DATASET_TIME_CHUNK_SIZE=10
  INTEGER,PARAMETER :: DOUBLE=8

  ! name for multiple datasets folder under a mesh, grid, or scattered dataset
  CHARACTER(LEN=*),PARAMETER :: MULTI_DATASET_LOCATION='Datasets'

  CHARACTER(LEN=*),PARAMETER :: MESH_NODE_GROUP='Nodes'
  CHARACTER(LEN=*),PARAMETER :: MESH_ELEM_GROUP='Elements'
  CHARACTER(LEN=*),PARAMETER :: MESH_ATT_NUMELEMS='NumElems'
  CHARACTER(LEN=*),PARAMETER :: MESH_ATT_NUMNODES='NumNodes'
  CHARACTER(LEN=*),PARAMETER :: MESH_DSET_NODE_LOCS='NodeLocs'
  CHARACTER(LEN=*),PARAMETER :: MESH_DSET_ELEM_TYPES='Types'
  CHARACTER(LEN=*),PARAMETER :: MESH_DSET_ELEM_NODE_IDS='Nodeids'

  CHARACTER(LEN=*),PARAMETER :: MESH_PROP_GROUP='PROPERTIES'
  CHARACTER(LEN=*),PARAMETER :: MESH_NODE_PROP_GROUP='Nodes\\PROPERTIES'
  CHARACTER(LEN=*),PARAMETER :: MESH_ELEM_PROP_GROUP='Nodes\\PROPERTIES'

    ! Attributes & datasets used for grids
  CHARACTER(LEN=*),PARAMETER :: GRID_ATT_TYPE='Gridtype'
  CHARACTER(LEN=*),PARAMETER :: GRID_ATT_NUM_DIMENSIONS='Numdimensions'
  CHARACTER(LEN=*),PARAMETER :: GRID_ATT_EXTRUDE_TYPE='Extrudetype'
  CHARACTER(LEN=*),PARAMETER :: GRID_ATT_ORIGIN='Origin'
  CHARACTER(LEN=*),PARAMETER :: GRID_ATT_ORIENTATION='Orientation'
  CHARACTER(LEN=*),PARAMETER :: GRID_ATT_BEARING='Bearing'
  CHARACTER(LEN=*),PARAMETER :: GRID_ATT_DIP='Dip'
  CHARACTER(LEN=*),PARAMETER :: GRID_ATT_ROLL='Roll'
  CHARACTER(LEN=*),PARAMETER :: GRID_ATT_COMP_ORIG='CompOrigin'
  CHARACTER(LEN=*),PARAMETER :: GRID_ATT_U_DIR='Direction'
  CHARACTER(LEN=*),PARAMETER :: GRID_ATT_NUM_I='NumI'
  CHARACTER(LEN=*),PARAMETER :: GRID_ATT_NUM_J='NumJ'
  CHARACTER(LEN=*),PARAMETER :: GRID_ATT_NUM_K='NumK'
  CHARACTER(LEN=*),PARAMETER :: GRID_ATT_EXTRUDE_LAYERS='Extrudelayers'
  CHARACTER(LEN=*),PARAMETER :: GRID_DSET_COORDS_I='CoordsI'
  CHARACTER(LEN=*),PARAMETER :: GRID_DSET_COORDS_J='CoordsJ'
  CHARACTER(LEN=*),PARAMETER :: GRID_DSET_COORDS_K='CoordsK'
  CHARACTER(LEN=*),PARAMETER :: GRID_DSET_EXTRUDE_VALS='Extrudevals'

  CHARACTER(LEN=*),PARAMETER :: GRID_PROP_GROUP='PROPERTIES'
  CHARACTER(LEN=*),PARAMETER :: GRID_CELL_PROP_GROUP='GridCellProps'
  CHARACTER(LEN=*),PARAMETER :: GRID_NODE_PROP_GROUP='GridNodeProps'

  ! Reserved property names
  CHARACTER(LEN=*),PARAMETER :: PROP_ACTIVITY='activity'
  CHARACTER(LEN=*),PARAMETER :: PROP_NULL_VALUE='nullvalue'

  ! Dataset value locations for a grid
  INTEGER, PARAMETER :: GRID_LOC_CENTER = 0
  INTEGER, PARAMETER :: GRID_LOC_CORNER = 1
  INTEGER, PARAMETER :: GRID_LOC_FACES  = 2
  INTEGER, PARAMETER :: GRID_LOC_FACE_I = 3
  INTEGER, PARAMETER :: GRID_LOC_FACE_J = 4
  INTEGER, PARAMETER :: GRID_LOC_FACE_K = 5
  INTEGER, PARAMETER :: GRID_LOC_COLUMN = 6

  ! extrusion types for grids
  INTEGER, PARAMETER :: EXTRUDE_SIGMA           = 0
  INTEGER, PARAMETER :: EXTRUDE_CARTESIAN       = 1
  INTEGER, PARAMETER :: EXTRUDE_CURV_AT_CORNERS = 2
  INTEGER, PARAMETER :: EXTRUDE_CURV_AT_CELLS   = 3
  INTEGER, PARAMETER :: EXTRUDE_MIN             = 0
  INTEGER, PARAMETER :: EXTRUDE_MAX             = 3

  ! Grid orientation
  INTEGER, PARAMETER :: ORIENTATION_RIGHT_HAND = 0
  INTEGER, PARAMETER :: ORIENTATION_LEFT_HAND  = 1

  ! Grid types
  INTEGER, PARAMETER :: GRID_TYPE_CARTESIAN            = 0
  INTEGER, PARAMETER :: GRID_TYPE_CURVILINEAR          = 1
  INTEGER, PARAMETER :: GRID_TYPE_CARTESIAN_EXTRUDED   = 2
  INTEGER, PARAMETER :: GRID_TYPE_CURVILINEAR_EXTRUDED = 3
  INTEGER, PARAMETER :: GRID_TYPE_MIN                  = 0
  INTEGER, PARAMETER :: GRID_TYPE_MAX                  = 3

! Attributes used for geometric paths
  CHARACTER(LEN=*),PARAMETER :: GEOMPATH_ATT_COMPRESSION = "GeompathCompression"
  CHARACTER(LEN=*),PARAMETER :: GEOMPATH_ATT_NULLVALUE = "GeompathNullvalue"
  CHARACTER(LEN=*),PARAMETER :: GEOMPATH_ATT_MINS = "GeompathMins"
  CHARACTER(LEN=*),PARAMETER :: GEOMPATH_ATT_MAXS = "GeompathMaxs"
  CHARACTER(LEN=*),PARAMETER :: GEOMPATH_DSET_LOCS =  "Locations"
  CHARACTER(LEN=*),PARAMETER :: GEOMPATH_DSET_TIMES =  "Times"
  INTEGER, PARAMETER :: GEOMPATH_TIME_CHUNK_SIZE =  10

! type defines for coordinate conversions
  ! Horizontal datum
  INTEGER, PARAMETER :: HORIZ_DATUM_LOCAL              = 0
  INTEGER, PARAMETER :: HORIZ_DATUM_GEOGRAPHIC         = 1
  INTEGER, PARAMETER :: HORIZ_DATUM_GEOGRAPHIC_NAD27   = 2
  INTEGER, PARAMETER :: HORIZ_DATUM_GEOGRAPHIC_NAD83   = 3
  INTEGER, PARAMETER :: HORIZ_DATUM_GEOGRAPHIC_HPGN    = 4
  INTEGER, PARAMETER :: HORIZ_DATUM_UTM                = 5
  INTEGER, PARAMETER :: HORIZ_DATUM_UTM_NAD27          = 6
  INTEGER, PARAMETER :: HORIZ_DATUM_UTM_NAD83          = 7
  INTEGER, PARAMETER :: HORIZ_DATUM_UTM_HPGN           = 8
  INTEGER, PARAMETER :: HORIZ_DATUM_STATE_PLANE_NAD27  = 9
  INTEGER, PARAMETER :: HORIZ_DATUM_STATE_PLANE_NAD83  =10
  INTEGER, PARAMETER :: HORIZ_DATUM_STATE_PLANE_HPGN   =11
  INTEGER, PARAMETER :: HORIZ_DATUM_CPP                =12

  ! Horizontal or vertical units
  INTEGER, PARAMETER :: COORD_UNITS_US_FEET            = 0
  INTEGER, PARAMETER :: COORD_UNITS_INTERNATIONAL_FEET = 1
  INTEGER, PARAMETER :: COORD_UNITS_METERS             = 2

  ! vertical datum
  INTEGER, PARAMETER :: VERT_DATUM_LOCAL     = 0
  INTEGER, PARAMETER :: VERT_DATUM_NGVD_29   = 1
  INTEGER, PARAMETER :: VERT_DATUM_NGVD_88   = 2

  ! latitude
  INTEGER, PARAMETER :: LATITUDE_NORTH       = 0
  INTEGER, PARAMETER :: LATITUDE_SOUTH       = 1

  ! longitude
  INTEGER, PARAMETER :: LONGITUDE_EAST       = 0
  INTEGER, PARAMETER :: LONGITUDE_WEST       = 1

  ! valid UTM zones
  INTEGER, PARAMETER :: UTM_ZONE_MIN         = 1
  INTEGER, PARAMETER :: UTM_ZONE_MAX         =60

  ! Attributes used for coordinate groups
  CHARACTER(LEN=*),PARAMETER :: COORD_ATT_HORIZ_DATUM   = 'HorizontalDatum'
  CHARACTER(LEN=*),PARAMETER :: COORD_ATT_HORIZ_UNITS   = 'HorizontalUnits'
  CHARACTER(LEN=*),PARAMETER :: COORD_ATT_VERT_DATUM    = 'VerticalDatum'
  CHARACTER(LEN=*),PARAMETER :: COORD_ATT_VERT_UNITS    = 'VerticalUnits'
  CHARACTER(LEN=*),PARAMETER :: COORD_ATT_LATITUDE      = 'Latitude'
  CHARACTER(LEN=*),PARAMETER :: COORD_ATT_LONGITUDE     = 'Longitude'
  CHARACTER(LEN=*),PARAMETER :: COORD_ATT_UTM_ZONE      = 'UtmZone'
  CHARACTER(LEN=*),PARAMETER :: COORD_ATT_SPC_ZONE      = 'SpcZone'
  CHARACTER(LEN=*),PARAMETER :: COORD_ATT_HPGN_AREA     = 'HpgnArea'
  CHARACTER(LEN=*),PARAMETER :: COORD_ATT_CPP_LATITUDE  = 'Latitude'
  CHARACTER(LEN=*),PARAMETER :: COORD_ATT_CPP_LONGITUDE = 'Longitude'
  CHARACTER(LEN=*),PARAMETER :: COORD_ATT_ELLIPSE       = 'Ellipse'
  CHARACTER(LEN=*),PARAMETER :: COORD_ATT_MAJOR_R       = 'MajorR'
  CHARACTER(LEN=*),PARAMETER :: COORD_ATT_MINOR_R       = 'MinorR'
  CHARACTER(LEN=*),PARAMETER :: COORD_ATT_VERSION       = 'Version'

! Calendar Eras
  INTEGER, PARAMETER :: ERA_IS_BCE = 0
  INTEGER, PARAMETER :: ERA_IS_CE  = 1

! Datasets
  CHARACTER(LEN=*),PARAMETER :: DATASET_UNITS  = 'DatasetUnits'
  CHARACTER(LEN=*),PARAMETER :: DATASET_COMPRESSION = 'DatasetCompression'
  CHARACTER(LEN=*),PARAMETER :: DATASET_NULLVALUE = 'DatasetNullvalue'

  ! Element-Related Parameters
  INTEGER, PARAMETER :: EL_TYPE_TRI_LINEAR = 200
  INTEGER, PARAMETER :: EL_TYPE_QUADRILATERAL_LINEAR = 210

  ! Used to determine chunking for properties
  INTEGER, PARAMETER :: XF_PROP_CHUNK_NUMBER = 100
  INTEGER, PARAMETER :: XF_PROP_CHUNK_MIN = 10
  INTEGER, PARAMETER :: XF_PROP_CHUNK_MAX = 1000

  ! The number below is there because I couldn't get unlimited to work
  INTEGER, PARAMETER :: REALLY_BIG_NUMBER = 214783647

  ! Attributes for datasets
  CHARACTER(LEN=*), PARAMETER :: DATASET_ATT_LOCALCOORDS = 'LocalCoords'
  CHARACTER(LEN=*), PARAMETER :: DATASET_ATT_DATALOCATION = 'DatasetLocation'
  CHARACTER(LEN=*), PARAMETER :: DATASET_ATT_EXTENDFILL   = &
                                                       'DatasetExtendFillValue'

  ! Attribute of integer - used in writing
  INTEGER     XF_INTEGER_WRITE_TYPE
  INTEGER     XF_UINT_WRITE_TYPE

  ! Used to pass 1 as a size_t
  INTEGER, PARAMETER :: N1_SIZE_T = 1

  INTEGER, PARAMETER :: MAX_NUMBER_DSET_DIMS = 7

  ! binary flags for HDF5 functions
  INTEGER, PARAMETER :: H5F_OBJ_FILE_F     = 1
  INTEGER, PARAMETER :: H5F_OBJ_DATASET_F  = 2
  INTEGER, PARAMETER :: H5F_OBJ_GROUP_F    = 4
  INTEGER, PARAMETER :: H5F_OBJ_DATATYPE_F = 8
  INTEGER, PARAMETER :: H5F_OBJ_ATTR_F     = 10
  INTEGER, PARAMETER :: H5F_OBJ_ALL_F      = 31


END MODULE XMDFDEFS

!************************************************************************************

MODULE XMDF

USE ERRORDEFINITIONS
USE XMDFDEFS

CONTAINS

SUBROUTINE XF_INITIALIZE (error)
INTEGER, INTENT(OUT) :: error


!            INTEGER, EXTERNAL :: _xfpinitialize_f
!  MS FORTRAN needs explicit interface for C functions called here.
!
  INTERFACE
    INTEGER FUNCTION xfpinitialize_f()
	  !MS$ATTRIBUTES C,reference::xfpinitialize_f

    END FUNCTION xfpinitialize_f
  END INTERFACE

  error = xfpinitialize_f()

  return
END SUBROUTINE
  
!*************************

SUBROUTINE XF_CLOSE (error)
INTEGER, INTENT(OUT) :: error


!            INTEGER, EXTERNAL :: xfpclose_f
!  MS FORTRAN needs explicit interface for C functions called here.
!
  INTERFACE
    INTEGER FUNCTION xfpclose_f()
	  !MS$ATTRIBUTES C,reference::xfpclose_f

    END FUNCTION xfpclose_f
  END INTERFACE

  error = xfpclose_f()

  return
END SUBROUTINE

!*************************

SUBROUTINE XF_CREATE_FILE (a_Filename, a_Overwrite, Id, error)

  CHARACTER(LEN=*), INTENT(IN) :: a_Filename
	LOGICAL*2, INTENT(IN)        :: a_Overwrite
	INTEGER, INTENT(OUT)         :: Id
	INTEGER, INTENT(OUT)         :: error
	INTEGER                         namelen ! Length of the name


!            INTEGER, EXTERNAL :: _xfcreatefile_f
!  MS FORTRAN needs explicit interface for C functions called here.
!
  INTERFACE
    INTEGER FUNCTION xfcreatefile_f(a_Filename, namelen, Id, a_Overwrite)
	  !MS$ATTRIBUTES C,reference::xfcreatefile_f
	  !DEC$ATTRIBUTES reference :: a_Filename
    CHARACTER(LEN=*), INTENT(IN) :: a_Filename
	  INTEGER, INTENT(OUT)         :: Id
	  LOGICAL*2, INTENT(IN)        :: a_Overwrite
	  INTEGER, INTENT(IN)          :: namelen

    END FUNCTION xfcreatefile_f
  END INTERFACE

  namelen = LEN_TRIM(a_Filename)
  error = xfcreatefile_f(a_Filename, namelen, Id, a_Overwrite)

return
END SUBROUTINE

!*************************

SUBROUTINE XF_OPEN_FILE (a_Filename, a_ReadOnly, Id, error)

  CHARACTER(LEN=*), INTENT(IN) :: a_Filename
	LOGICAL*2, INTENT(IN)        :: a_ReadOnly
	INTEGER, INTENT(OUT)         :: Id
	INTEGER, INTENT(OUT)         :: error
	INTEGER                         namelen ! Length of the name

!            INTEGER, EXTERNAL :: _xfopenfile_f
!  MS FORTRAN needs explicit interface for C functions called here.
!
  INTERFACE
    INTEGER FUNCTION xfopenfile_f(a_Filename, namelen, Id, a_ReadOnly)
	  !MS$ATTRIBUTES C,reference::xfopenfile_f
	  !DEC$ATTRIBUTES reference :: a_Filename
    CHARACTER(LEN=*), INTENT(IN) :: a_Filename
	  INTEGER, INTENT(OUT)         :: Id
	  LOGICAL*2, INTENT(IN)          :: a_ReadOnly
	  INTEGER, INTENT(IN)          :: namelen

    END FUNCTION xfopenfile_f
  END INTERFACE

  namelen = LEN_TRIM(a_Filename)
  error = xfopenfile_f(a_Filename, namelen, Id, a_ReadOnly)
  

return

END SUBROUTINE

!****************************

SUBROUTINE XF_CLOSE_FILE (a_Id, error)

	INTEGER, INTENT(IN)  :: a_Id
	INTEGER, INTENT(OUT)         :: error

!            INTEGER, EXTERNAL :: _xfclosefile_f
!  MS FORTRAN needs explicit interface for C functions called here.
!
  INTERFACE
    INTEGER FUNCTION xfclosefile_f(a_Id)
	  !MS$ATTRIBUTES C,reference::xfclosefile_f
	  !DEC$ATTRIBUTES reference :: a_Id
	  INTEGER, INTENT(IN)       :: a_Id

    END FUNCTION xfclosefile_f
  END INTERFACE

  error = xfclosefile_f(a_Id)

return

END SUBROUTINE

!------------------------------------------------------------------------------
! SUBROUTINE  XF_GET_LIBRARY_VERSION
! PURPOSE     Returns the current version of XMDF
! NOTES     
!------------------------------------------------------------------------------
SUBROUTINE XF_GET_LIBRARY_VERSION (a_Version, error)

REAL*4, INTENT(OUT) :: a_Version
INTEGER, INTENT(OUT) :: error

!            INTEGER, EXTERNAL :: _xfgetlibraryversion_f
!  MS FORTRAN needs explicit interface for C functions called here.
!
  INTERFACE
    INTEGER FUNCTION xfgetlibraryversion_f(a_Version)
	  !MS$ATTRIBUTES C,reference::xfgetlibraryversion_f
	  REAL*4, INTENT(OUT)         :: a_Version

    END FUNCTION xfgetlibraryversion_f
  END INTERFACE

  error = xfgetlibraryversion_f(a_Version)

return

END SUBROUTINE

!------------------------------------------------------------------------------
! SUBROUTINE  XF_GET_LIBRARY_VERSION_FILE
! PURPOSE     Obtain the version of XMDF from the library file
! NOTES     
!------------------------------------------------------------------------------
SUBROUTINE XF_GET_LIBRARY_VERSION_FILE (a_File, a_Version, error)

INTEGER, INTENT(IN) :: a_File
REAL*4, INTENT(OUT)          :: a_Version
INTEGER, INTENT(OUT)       :: error


!            INTEGER, EXTERNAL :: _xfGetLibraryVersionFile_f
!  MS FORTRAN needs explicit interface for C functions called here.
!
  INTERFACE
    INTEGER FUNCTION xfgetlibraryversionfile_f(a_File, a_Version)
	  !MS$ATTRIBUTES C,reference::xfgetlibraryversionfile_f
	  INTEGER, INTENT(IN)        :: a_File
	  REAL*4, INTENT(OUT)          :: a_Version

    END FUNCTION xfgetlibraryversionfile_f
  END INTERFACE

  error = xfgetlibraryversionfile_f(a_File, a_Version)

return

END SUBROUTINE
!------------------------------------------------------------------------------
! SUBROUTINE XF_CREATE_COORDINATE_GROUP
! PURPOSE    Creates a group to store coordinate information in.
! NOTES
!------------------------------------------------------------------------------
SUBROUTINE XF_CREATE_COORDINATE_GROUP (a_ParentId, a_ChildId, error)
INTEGER, INTENT(IN)  :: a_ParentId
INTEGER, INTENT(OUT) :: a_ChildId
INTEGER, INTENT(OUT)        :: error

!            INTEGER, EXTERNAL :: _xfCreateCoordinateGroup_f
!  MS FORTRAN needs explicit interface for C functions called here.
!
  INTERFACE
    INTEGER FUNCTION xfcreatecoordinategroup_f(a_ParentId, a_ChildId)
	  !MS$ATTRIBUTES C,reference::xfcreatecoordinategroup_f
    INTEGER, INTENT(IN)  :: a_ParentId
    INTEGER, INTENT(OUT) :: a_ChildId

    END FUNCTION xfcreatecoordinategroup_f
  END INTERFACE

  error = xfcreatecoordinategroup_f(a_ParentId, a_ChildId)

  return
END SUBROUTINE
!------------------------------------------------------------------------------
! SUBROUTINE XF_OPEN_COORDINATE_GROUP
! PURPOSE    Opens a coordinate group
! NOTES
!------------------------------------------------------------------------------
SUBROUTINE XF_OPEN_COORDINATE_GROUP (a_ParentId, a_ChildId, error)
INTEGER, INTENT(IN)  :: a_ParentId
INTEGER, INTENT(OUT) :: a_ChildId
INTEGER, INTENT(OUT)        :: error

!            INTEGER, EXTERNAL :: _xfOpenCoordinateGroup_f
!  MS FORTRAN needs explicit interface for C functions called here.
!
  INTERFACE
    INTEGER FUNCTION xfopencoordinategroup_f(a_ParentId, a_ChildId)
	  !MS$ATTRIBUTES C,reference::xfopencoordinategroup_f
    INTEGER, INTENT(IN)  :: a_ParentId
    INTEGER, INTENT(OUT) :: a_ChildId

    END FUNCTION xfopencoordinategroup_f
  END INTERFACE

  error = xfopencoordinategroup_f(a_ParentId, a_ChildId)
  
  return
END SUBROUTINE

!******************************

SUBROUTINE XF_CREATE_PROPERTY_GROUP (a_ParentId, Id, error)

INTEGER, INTENT(IN)  :: a_ParentId
INTEGER, INTENT(OUT) :: Id
INTEGER, INTENT(OUT) :: error

!            INTEGER, EXTERNAL :: _xfCreatePropertyGroup_f
!  MS FORTRAN needs explicit interface for C functions called here.
!
  INTERFACE
    INTEGER FUNCTION xfcreatepropertygroup_f(a_ParentId, Id)
	  !MS$ATTRIBUTES C,reference::xfcreatepropertygroup_f
    INTEGER, INTENT(IN)  :: a_ParentId
    INTEGER, INTENT(OUT) :: Id

    END FUNCTION xfcreatepropertygroup_f
  END INTERFACE

  error = xfcreatepropertygroup_f(a_ParentId, Id)

return

END SUBROUTINE XF_CREATE_PROPERTY_GROUP

!******************************

SUBROUTINE XF_WRITE_PROPERTY_STRING (a_Id, a_Name, a_Attributes, error)

INTEGER, INTENT(IN)   :: a_Id
CHARACTER(LEN=*), INTENT(IN) :: a_Name, a_Attributes
INTEGER          strlen, Number, tmplen 
INTEGER          error

!            INTEGER, EXTERNAL :: _xfWritePropertyString_f
!  MS FORTRAN needs explicit interface for C functions called here.
!
  INTERFACE
    INTEGER FUNCTION xfwritepropertystring_f(a_Id, a_Name, tmplen, &
                                             Number, a_Attributes, strlen)
	!MS$ATTRIBUTES C,reference::xfwritepropertystring_f
	!DEC$ATTRIBUTES reference :: a_Name
	!DEC$ATTRIBUTES reference :: a_Attributes
    INTEGER, INTENT(IN)          :: a_Id
    CHARACTER(LEN=*), INTENT(IN) :: a_Name, a_Attributes
	INTEGER, INTENT(IN)          :: strlen, Number, tmplen

    END FUNCTION xfwritepropertystring_f
  END INTERFACE

  strlen = LEN_TRIM(a_Name)
  Number = 1
  tmplen = LEN_TRIM(a_Attributes)

  error = xfwritepropertystring_f(a_Id, a_Name, strlen, &
                                  Number, a_Attributes, tmplen)

return

END SUBROUTINE

!*****************************

SUBROUTINE XF_WRITE_PROPERTY_INT (a_Id, a_Name, a_Number, a_Properties, &
                                  a_Compression, error)
INTEGER,                INTENT(IN)  :: a_Id
CHARACTER(LEN=*),              INTENT(IN)  :: a_Name
INTEGER,                       INTENT(IN)  :: a_Number
INTEGER, DIMENSION(*),         INTENT(IN)  :: a_Properties
INTEGER,                       INTENT(IN)  :: a_Compression
INTEGER,                       INTENT(OUT) :: error
INTEGER                                       namelen
  

!            INTEGER, EXTERNAL :: _xfWritePropertyInt_f
!  MS FORTRAN needs explicit interface for C functions called here.
!
  INTERFACE
    INTEGER FUNCTION xfwritepropertyint_f(a_Id, a_Name, namelen, a_Number, & 
						                              a_Properties, a_Compression)
	  !MS$ATTRIBUTES C,reference::xfwritepropertyint_f
    !DEC$ATTRIBUTES reference :: a_Name
    INTEGER,               INTENT(IN)  :: a_Id
    CHARACTER(LEN=*),      INTENT(IN)  :: a_Name
    INTEGER,               INTENT(IN)  :: a_Number
    INTEGER, DIMENSION(*), INTENT(IN)  :: a_Properties 
    INTEGER,               INTENT(IN)  :: a_Compression
    INTEGER,               INTENT(IN)  :: namelen

    END FUNCTION xfwritepropertyint_f
  END INTERFACE

  namelen = LEN_TRIM(a_Name)
  error = xfwritepropertyint_f(a_Id, a_Name, namelen, a_Number, &
						                   a_Properties, a_Compression)

return

END SUBROUTINE

!------------------------------------------------------------------------------
! SUBROUTINE  XF_WRITE_PROPERTY_UINT
! PURPOSE     writes a unsigned integer attribute (really a dataset) to the 
!             folder
! NOTES
!------------------------------------------------------------------------------
SUBROUTINE XF_WRITE_PROPERTY_UINT (a_Id, a_Name, a_Number, a_Properties, &
                                  a_Compression, error)
INTEGER,        INTENT(IN)          :: a_Id
CHARACTER(LEN=*),      INTENT(IN)   :: a_Name
INTEGER,               INTENT(IN)   :: a_Number
INTEGER, DIMENSION(*), INTENT(IN)   :: a_Properties
INTEGER,               INTENT(IN)   :: a_Compression
INTEGER,               INTENT(OUT)  :: error
INTEGER                                namelen


!            INTEGER, EXTERNAL :: _xfWritePropertyUnsignedInt_f
!  MS FORTRAN needs explicit interface for C functions called here.
!
  INTERFACE
    INTEGER FUNCTION xfwritepropertyunsignedint_f(a_Id, a_Name, namelen, a_Number, &
                                                  a_Properties, a_Compression)
	  !MS$ATTRIBUTES C,reference::xfwritepropertyunsignedint_f
    !DEC$ATTRIBUTES reference :: a_Name
    INTEGER,               INTENT(IN)  :: a_Id
    CHARACTER(LEN=*),      INTENT(IN)  :: a_Name
    INTEGER,               INTENT(IN)  :: a_Number
    INTEGER, DIMENSION(*), INTENT(IN)  :: a_Properties
    INTEGER,               INTENT(IN)  :: a_Compression
    INTEGER,               INTENT(IN)  :: namelen

    END FUNCTION xfwritepropertyunsignedint_f
  END INTERFACE

  namelen = LEN_TRIM(a_Name)
  error = xfwritepropertyunsignedint_f(a_Id, a_Name, namelen, a_Number, a_Properties, & 
                                       a_Compression)

return

END SUBROUTINE

!*****************************

SUBROUTINE XF_WRITE_PROPERTY_DOUBLE (a_Id, a_Name, a_Number, a_Properties, &
                                     a_Compression, error)
INTEGER, INTENT(IN)                 :: a_Id
CHARACTER(LEN=*), INTENT(IN)        :: a_Name
INTEGER, INTENT(IN)                 :: a_Number
REAL*8, DIMENSION(*), INTENT(IN)    :: a_Properties 
INTEGER, INTENT(IN)                 :: a_Compression
INTEGER                             :: namelen
INTEGER, INTENT(OUT)                :: error


!            INTEGER, EXTERNAL :: _xfWritePropertyDouble_f
!  MS FORTRAN needs explicit interface for C functions called here.
!
  INTERFACE
    INTEGER FUNCTION xfwritepropertydouble_f(a_Id, a_Name, namelen, a_Number, a_Properties, & 
                                             a_Compression)
	  !MS$ATTRIBUTES C,reference::xfwritepropertydouble_f
    !DEC$ATTRIBUTES reference :: a_Name
    INTEGER,              INTENT(IN)  :: a_Id
    CHARACTER(LEN=*),     INTENT(IN)  :: a_Name
    INTEGER,              INTENT(IN)  :: a_Number
    REAL*8, DIMENSION(*), INTENT(IN)   :: a_Properties 
    INTEGER,              INTENT(IN)  :: a_Compression
    INTEGER,              INTENT(IN)  :: namelen

    END FUNCTION xfwritepropertydouble_f
  END INTERFACE

  namelen = LEN_TRIM(a_Name)
  error = xfwritepropertydouble_f(a_Id, a_Name, namelen, a_Number, a_Properties, &
                                  a_Compression)

return

END SUBROUTINE

!*****************************

!-----------------------------------------------------------------------------
!  SUBROUTINE  xfwritepropertyfloat
!  PURPOSE   
!  NOTES     
!-----------------------------------------------------------------------------
SUBROUTINE XF_WRITE_PROPERTY_FLOAT (a_Id, a_Name, a_Number, a_Properties, &
                                    a_Compression, error)
INTEGER, INTENT(IN)          :: a_Id
CHARACTER(LEN=*), INTENT(IN) :: a_Name
INTEGER, INTENT(IN)          :: a_Number
REAL*4, INTENT(IN)             :: a_Properties
INTEGER, INTENT(IN)          :: a_Compression
INTEGER, INTENT(OUT)         :: error
INTEGER                         namelen


!            INTEGER, EXTERNAL :: _xfWritePropertyFloat_f
!  MS FORTRAN needs explicit interface for C functions called here.
!
  INTERFACE
    INTEGER FUNCTION xfwritepropertyfloat_f(a_Id, a_Name, namelen, a_Number, a_Properties, & 
                                            a_Compression)
	  !MS$ATTRIBUTES C,reference::xfwritepropertyfloat_f
    !DEC$ATTRIBUTES reference :: a_Name
    INTEGER,          INTENT(IN)  :: a_Id
    CHARACTER(LEN=*), INTENT(IN)  :: a_Name
    INTEGER,          INTENT(IN)  :: a_Number
    REAL*4,             INTENT(IN)  :: a_Properties
    INTEGER,          INTENT(IN)  :: a_Compression
    INTEGER,          INTENT(IN)  :: namelen

    END FUNCTION xfwritepropertyfloat_f
  END INTERFACE

  namelen = LEN_TRIM(a_Name)
  error = xfwritepropertyfloat_f(a_Id, a_Name, namelen, a_Number, a_Properties, &
                                 a_Compression)

return

END SUBROUTINE

!*****************************

SUBROUTINE XF_DOES_PROP_W_NAME_EXIST (a_Id, a_Name, a_Exists, error)

LOGICAL*2        a_Exists
INTEGER        a_Id
INTEGER        error
CHARACTER      a_Name
INTEGER        namelen


!            INTEGER, EXTERNAL :: _xfDoesPropertyWithNameExist_f
!  MS FORTRAN needs explicit interface for C functions called here.
!
  INTERFACE
    INTEGER FUNCTION xfdoespropertywithnameexist_f(a_Id, a_Name, namelen, a_Exists)
	  !MS$ATTRIBUTES C,reference::xfdoespropertywithnameexist_f
    !DEC$ATTRIBUTES reference :: a_Name
    INTEGER           a_Id
    CHARACTER         a_Name
    LOGICAL*2           a_Exists
    INTEGER           namelen

    END FUNCTION xfdoespropertywithnameexist_f
  END INTERFACE

  namelen = LEN_TRIM(a_Name)
  error = xfdoespropertywithnameexist_f(a_Id, a_Name, namelen, a_Exists)


return

END SUBROUTINE

!*****************************

SUBROUTINE XF_GET_PROPERTY_NUMBER (a_Id, a_Name, a_Number, error)

INTEGER,          INTENT(IN)  :: a_Id
INTEGER,          INTENT(OUT) :: a_Number
CHARACTER(LEN=*), INTENT(IN)  :: a_Name
INTEGER                          namelen
INTEGER,          INTENT(OUT) :: error


!            INTEGER, EXTERNAL :: _xfGetPropertyNumber_f
!  MS FORTRAN needs explicit interface for C functions called here.
!
  INTERFACE
    INTEGER FUNCTION xfgetpropertynumber_f(a_Id, a_Name, namelen, a_Number)
	  !MS$ATTRIBUTES C,reference::xfgetpropertynumber_f
    !DEC$ATTRIBUTES reference :: a_Name
    INTEGER,          INTENT(IN)  :: a_Id
    INTEGER,          INTENT(OUT) :: a_Number
    CHARACTER(LEN=*), INTENT(IN)  :: a_Name
    INTEGER,          INTENT(IN)  :: namelen

    END FUNCTION xfgetpropertynumber_f
  END INTERFACE

  namelen = LEN_TRIM(a_Name)
  error = xfgetpropertynumber_f (a_Id, a_Name, namelen, a_Number)


return

END SUBROUTINE

!*****************************

SUBROUTINE XF_READ_PROPERTY_INT (a_Id, a_Name, a_Number, a_Properties, error)

INTEGER,          INTENT(IN)  :: a_Id
INTEGER,          INTENT(IN)  :: a_Number
INTEGER,          INTENT(OUT) :: a_Properties
CHARACTER(LEN=*), INTENT(IN)  :: a_Name
INTEGER                          namelen
INTEGER,          INTENT(OUT) :: error


!            INTEGER, EXTERNAL :: _xfReadPropertyInt_f
!  MS FORTRAN needs explicit interface for C functions called here.
!
  INTERFACE
    INTEGER FUNCTION xfreadpropertyint_f(a_Id, a_Name, namelen, a_Number, a_Properties)
	  !MS$ATTRIBUTES C,reference::xfreadpropertyint_f
    !DEC$ATTRIBUTES reference :: a_Name
    INTEGER,          INTENT(IN)  :: a_Id, a_Number
    INTEGER,          INTENT(OUT) :: a_Properties
    CHARACTER(LEN=*), INTENT(IN)  :: a_Name
    INTEGER,          INTENT(IN)  :: namelen

    END FUNCTION xfreadpropertyint_f
  END INTERFACE

  namelen = LEN_TRIM(a_Name)
  error = xfreadpropertyint_f (a_Id, a_Name, namelen, a_Number, a_Properties)


return

END SUBROUTINE

!*********************************

SUBROUTINE XF_GET_PROPERTY_STRING_LENGTH (a_Id, a_Name, a_Number, a_MaxLength, error)

INTEGER          a_Id
CHARACTER(LEN=*) a_Name
INTEGER          a_Number
INTEGER          a_MaxLength
INTEGER          namelen
INTEGER          error


!            INTEGER, EXTERNAL :: _xfGetPropertyStringLength_f
!  MS FORTRAN needs explicit interface for C functions called here.
!
  INTERFACE
    INTEGER FUNCTION xfgetpropertystringlength_f(a_Id, a_Name, namelen, a_Number, &
                                                 a_MaxLength)
	  !MS$ATTRIBUTES C,reference::xfgetpropertystringlength_f
    !DEC$ATTRIBUTES reference :: a_Name
    INTEGER           a_Id
    CHARACTER         a_Name
    INTEGER           a_Number
    INTEGER           a_MaxLength
    INTEGER           namelen

    END FUNCTION xfgetpropertystringlength_f
  END INTERFACE

  namelen = LEN_TRIM(a_Name)
  error = xfgetpropertystringlength_f(a_Id, a_Name, namelen, a_Number, a_MaxLength)


return

END SUBROUTINE

!------------------------------------------------------------------------------
!  SUBROUTINE  XF_GET_PROPERTY_TYPE
!  PURPOSE     Gets the property type from a dataset
!  NOTES
!------------------------------------------------------------------------------
SUBROUTINE XF_GET_PROPERTY_TYPE (a_GroupId, a_Name, a_Type, error)
INTEGER,          INTENT(IN)  :: a_GroupId
CHARACTER(LEN=*), INTENT(IN)  :: a_Name
INTEGER,          INTENT(OUT) :: a_Type
INTEGER,          INTENT(OUT) :: error
INTEGER                          namelen
    
    
!            INTEGER, EXTERNAL :: _xfGetPropertyType_f
!  MS FORTRAN needs explicit interface for C functions called here.
!
  INTERFACE
    INTEGER FUNCTION xfgetpropertytype_f(a_GroupId, a_Name, namelen, a_Type)
	  !MS$ATTRIBUTES C,reference::xfgetpropertytype_f
    !DEC$ATTRIBUTES reference :: a_Name
    INTEGER,          INTENT(IN)  :: a_GroupId
    CHARACTER(LEN=*), INTENT(IN)  :: a_Name
    INTEGER,          INTENT(OUT) :: a_Type
    INTEGER,          INTENT(IN)  :: namelen

    END FUNCTION xfgetpropertytype_f
  END INTERFACE

  namelen = LEN_TRIM(a_Name)
  error = xfgetpropertytype_f(a_GroupId, a_Name, namelen, a_Type)
    
    
  return

END SUBROUTINE

!******************

SUBROUTINE XF_READ_PROPERTY_STRING (a_Id, a_Name, a_Number, a_MaxLength, &
                                    a_Attributes, error)

INTEGER          a_Id
INTEGER          a_MaxLength, a_Number
CHARACTER(LEN=*) a_Name, a_Attributes
INTEGER          namelen, attrlen
INTEGER          error


!            INTEGER, EXTERNAL :: _xfReadPropertyString_f
!  MS FORTRAN needs explicit interface for C functions called here.
!
  INTERFACE
    INTEGER FUNCTION xfreadpropertystring_f(a_Id, a_Name, namelen, a_Number, &
                                         a_MaxLength, a_Attributes, attrlen)
	  !MS$ATTRIBUTES C,reference::xfreadpropertystring_f
    !DEC$ATTRIBUTES reference :: a_Name
    INTEGER           a_Id
    CHARACTER(LEN=*)  a_Name, a_Attributes
    INTEGER           a_Number
    INTEGER           a_MaxLength
    INTEGER           namelen, attrlen

    END FUNCTION xfreadpropertystring_f
  END INTERFACE

  namelen = LEN_TRIM(a_Name)
  attrlen = LEN_TRIM(a_Attributes)
  error = xfreadpropertystring_f(a_Id, a_Name, namelen, a_Number, a_MaxLength, &
                              a_Attributes, attrlen)

  return

END SUBROUTINE

!******************

!------------------------------------------------------------------------------
! SUBROUTINE  XF_READ_PROPERTY_FLOAT
! PURPOSE     Reads the float dataset from the properties directory
! NOTES       the variable properties must already be allocated to a_Number
!------------------------------------------------------------------------------
SUBROUTINE XF_READ_PROPERTY_FLOAT (a_Id, a_Name, a_Number, a_Properties, error)
INTEGER, INTENT(IN)          :: a_Id
CHARACTER(LEN=*), INTENT(IN) :: a_Name
INTEGER, INTENT(IN)          :: a_Number
REAL*4, INTENT(OUT)            :: a_Properties
INTEGER, INTENT(OUT)         :: error
INTEGER                         namelen


!            INTEGER, EXTERNAL :: _xfReadPropertyFloat_f
!  MS FORTRAN needs explicit interface for C functions called here.
!
  INTERFACE
    INTEGER FUNCTION xfreadpropertyfloat_f(a_Id, a_Name, namelen, a_Number, &
                                           a_Properties)
	  !MS$ATTRIBUTES C,reference::xfreadpropertyfloat_f
    !DEC$ATTRIBUTES reference :: a_Name
    INTEGER,          INTENT(IN)  :: a_Id
    CHARACTER(LEN=*), INTENT(IN)  :: a_Name
    INTEGER,          INTENT(IN)  :: a_Number
    REAL*4,             INTENT(OUT) :: a_Properties
    INTEGER,          INTENT(IN)  :: namelen

    END FUNCTION xfreadpropertyfloat_f
  END INTERFACE

  namelen = LEN_TRIM(a_Name)
  error = xfreadpropertyfloat_f(a_Id, a_Name, namelen, a_Number, a_Properties)


return

END SUBROUTINE

!------------------------------------------------------------------------------
! SUBROUTINE  XF_READ_PROPERTY_DOUBLE
! PURPOSE     Reads the double dataset from the properties directory
! NOTES       the variable properties must already be allocated to a_Number
!------------------------------------------------------------------------------
SUBROUTINE XF_READ_PROPERTY_DOUBLE (a_Id, a_Name, a_Number, a_Properties, error)
INTEGER, INTENT(IN)          :: a_Id
CHARACTER(LEN=*), INTENT(IN) :: a_Name
INTEGER, INTENT(IN)          :: a_Number
REAL*8, DIMENSION(*), INTENT(OUT) :: a_Properties
INTEGER, INTENT(OUT)         :: error
INTEGER                         namelen


!            INTEGER, EXTERNAL :: _xfReadPropertyDouble_f
!  MS FORTRAN needs explicit interface for C functions called here.
!
  INTERFACE
    INTEGER FUNCTION xfreadpropertydouble_f(a_Id, a_Name, namelen, a_Number, &
                                           a_Properties)
	  !MS$ATTRIBUTES C,reference::xfreadpropertydouble_f
    !DEC$ATTRIBUTES reference :: a_Name
    INTEGER,          INTENT(IN)  :: a_Id
    CHARACTER(LEN=*), INTENT(IN)  :: a_Name
    INTEGER,          INTENT(IN)  :: a_Number
    REAL*8, DIMENSION(*), INTENT(OUT) :: a_Properties
    INTEGER,          INTENT(IN)  :: namelen

    END FUNCTION xfreadpropertydouble_f
  END INTERFACE

  namelen = LEN_TRIM(a_Name)
  error = xfreadpropertydouble_f(a_Id, a_Name, namelen, a_Number, a_Properties)


return

END SUBROUTINE

!**************************

SUBROUTINE XF_CREATE_GENERIC_GROUP (a_FileId, a_Path, a_GroupId, error)
INTEGER,          INTENT(IN)  :: a_FileId
CHARACTER(LEN=*), INTENT(IN)  :: a_Path
INTEGER,          INTENT(OUT) :: a_GroupId
INTEGER,          INTENT(OUT) :: error
INTEGER                          pathlen

!            INTEGER, EXTERNAL :: _xfCreateGenericGroup_f
!  MS FORTRAN needs explicit interface for C functions called here.
!
  INTERFACE
    INTEGER FUNCTION xfcreategenericgroup_f(a_FileId, a_Path, pathlen, a_GroupId)
	  !MS$ATTRIBUTES C,reference::xfcreategenericgroup_f
    !DEC$ATTRIBUTES reference :: a_Path
    INTEGER,          INTENT(IN)  :: a_FileId
    CHARACTER(LEN=*), INTENT(IN)  :: a_Path
    INTEGER,          INTENT(OUT) :: a_GroupId
    INTEGER,          INTENT(IN)  :: pathlen

    END FUNCTION xfcreategenericgroup_f
  END INTERFACE

  pathlen = LEN_TRIM(a_Path)
  error = xfcreategenericgroup_f(a_FileId, a_Path, pathlen, a_GroupId)

return

END SUBROUTINE

!**************************

SUBROUTINE XF_CREATE_GROUP_FOR_MESH (a_FileId, a_Path, a_GroupId, error)
INTEGER,          INTENT(IN)   :: a_FileId
CHARACTER(LEN=*), INTENT(IN)   :: a_Path
INTEGER,          INTENT(OUT)  :: a_GroupId
INTEGER,          INTENT(OUT)  :: error
INTEGER                           pathlen


!            INTEGER, EXTERNAL :: _xfCreateGroupForMesh_f
!  MS FORTRAN needs explicit interface for C functions called here.
!
  INTERFACE
    INTEGER FUNCTION xfcreategroupformesh_f(a_FileId, a_Path, pathlen, a_GroupId)
	  !MS$ATTRIBUTES C,reference::xfcreategroupformesh_f
    !DEC$ATTRIBUTES reference :: a_Path
    INTEGER,          INTENT(IN)  :: a_FileId
    CHARACTER(LEN=*), INTENT(IN)  :: a_Path
    INTEGER,          INTENT(OUT) :: a_GroupId
    INTEGER,          INTENT(IN)  :: pathlen

    END FUNCTION xfcreategroupformesh_f
  END INTERFACE

  pathlen = LEN_TRIM(a_Path)
  error = xfcreategroupformesh_f(a_FileId, a_Path, pathlen, a_GroupId)


return

END SUBROUTINE

!**************************

SUBROUTINE XF_CREATE_GROUP_FOR_GRID (a_FileId, a_Path, a_GroupId, error)
INTEGER, INTENT(IN)   :: a_FileId
CHARACTER(LEN=*), INTENT(IN) :: a_Path
INTEGER, INTENT(OUT)  :: a_GroupId
INTEGER, INTENT(OUT)         :: error
INTEGER                         pathlen

!            INTEGER, EXTERNAL :: _xfCreateGroupForGrid_f
!  MS FORTRAN needs explicit interface for C functions called here.
!
  INTERFACE
    INTEGER FUNCTION xfcreategroupforgrid_f(a_FileId, a_Path, pathlen, a_GroupId)
	  !MS$ATTRIBUTES C,reference::xfcreategroupforgrid_f
    !DEC$ATTRIBUTES reference :: a_Path
    INTEGER,          INTENT(IN)  :: a_FileId
    CHARACTER(LEN=*), INTENT(IN)  :: a_Path
    INTEGER,          INTENT(OUT) :: a_GroupId
    INTEGER,          INTENT(IN)  :: pathlen

    END FUNCTION xfcreategroupforgrid_f
  END INTERFACE

  pathlen = LEN_TRIM(a_Path)
  error = xfcreategroupforgrid_f(a_FileId, a_Path, pathlen, a_GroupId)


return

END SUBROUTINE

!**************************

SUBROUTINE XF_CREATE_GROUP_FOR_XSEC (a_FileId, a_Path, a_GroupId, error)
INTEGER, INTENT(IN)   :: a_FileId
CHARACTER(LEN=*), INTENT(IN) :: a_Path
INTEGER, INTENT(OUT)  :: a_GroupId
INTEGER, INTENT(OUT)         :: error
INTEGER                         pathlen

!            INTEGER, EXTERNAL :: _xfCreateGroupForXsec_f
!  MS FORTRAN needs explicit interface for C functions called here.
!
  INTERFACE
    INTEGER FUNCTION xfcreategroupforxsec_f(a_FileId, a_Path, pathlen, a_GroupId)
	  !MS$ATTRIBUTES C,reference::xfcreategroupforxsec_f
    !DEC$ATTRIBUTES reference :: a_Path
    INTEGER,          INTENT(IN)  :: a_FileId
    CHARACTER(LEN=*), INTENT(IN)  :: a_Path
    INTEGER,          INTENT(OUT) :: a_GroupId
    INTEGER,          INTENT(IN)  :: pathlen

    END FUNCTION xfcreategroupforxsec_f
  END INTERFACE

  pathlen = LEN_TRIM(a_Path)
  error = xfcreategroupforxsec_f(a_FileId, a_Path, pathlen, a_GroupId)


return

END SUBROUTINE

!**************************

SUBROUTINE XF_OPEN_GROUP (a_ParentId, a_Path, a_GroupId, error)

INTEGER, INTENT(IN)          :: a_ParentId
CHARACTER(LEN=*), INTENT(IN) :: a_Path
INTEGER, INTENT(OUT)         :: a_GroupId
INTEGER, INTENT(OUT)         :: error
INTEGER    pathlen

!            INTEGER, EXTERNAL :: _xfOpenGroup_f
!  MS FORTRAN needs explicit interface for C functions called here.
!
  INTERFACE
    INTEGER FUNCTION xfopengroup_f(a_ParentId, a_Path, pathlen, a_GroupId)
	  !MS$ATTRIBUTES C,reference::xfopengroup_f
	  !DEC$ATTRIBUTES reference :: a_Path
    INTEGER,          INTENT(IN)  :: a_ParentId
	CHARACTER(LEN=*), INTENT(IN)  :: a_Path
	INTEGER,          INTENT(OUT) :: a_GroupId
	INTEGER,          INTENT(IN)  :: pathlen

    END FUNCTION xfopengroup_f
  END INTERFACE

  pathlen = LEN_TRIM(a_Path)
  error = xfopengroup_f(a_ParentId, a_Path, pathlen, a_GroupId)

return

END SUBROUTINE

!******************************

SUBROUTINE XF_CLOSE_GROUP (a_GroupId, error)

INTEGER        a_GroupId
INTEGER        error

!            INTEGER, EXTERNAL :: _xfCloseGroup_f
!  MS FORTRAN needs explicit interface for C functions called here.
!
  INTERFACE
    INTEGER FUNCTION xfclosegroup_f(a_GroupId)
	  !MS$ATTRIBUTES C,reference::xfclosegroup_f
    INTEGER,          INTENT(OUT) :: a_GroupId

    END FUNCTION xfclosegroup_f
  END INTERFACE

  error = xfclosegroup_f(a_GroupId)

return

END SUBROUTINE

!**************************

SUBROUTINE XF_GRP_PTHS_SZ_FOR_MESHES (a_FileId, a_Number, a_Maxsize, error)
INTEGER, INTENT(IN) :: a_FileId
INTEGER                       a_Number
INTEGER, INTENT(OUT)       :: a_Maxsize
INTEGER, INTENT(OUT)       :: error

!            INTEGER, EXTERNAL :: _xfGetGroupPathsSizeForMeshes_f
!  MS FORTRAN needs explicit interface for C functions called here.
!
  INTERFACE
    INTEGER FUNCTION xfgetgrouppathssizeformeshes_f(a_FileId, a_Number, a_Maxsize)
	  !MS$ATTRIBUTES C,reference::xfgetgrouppathssizeformeshes_f
    INTEGER,          INTENT(IN)  :: a_FileId
    INTEGER,          INTENT(IN)  :: a_Number
    INTEGER,          INTENT(OUT) :: a_Maxsize

    END FUNCTION xfgetgrouppathssizeformeshes_f
  END INTERFACE

  error = xfgetgrouppathssizeformeshes_f(a_FileId, a_Number, a_Maxsize)

return

END SUBROUTINE

!**************************

SUBROUTINE XF_GET_GROUP_PATHS_FOR_MESHES (a_FileId, a_Num, a_Maxsize, a_Paths, error)
INTEGER, INTENT(IN)      :: a_FileId
INTEGER, INTENT(IN)             :: a_Num
INTEGER, INTENT(IN)             :: a_Maxsize
INTEGER, INTENT(OUT)            :: error
CHARACTER,  DIMENSION(*)        :: a_Paths
INTEGER                            pathlen

!            INTEGER, EXTERNAL :: _xfGetGroupPathsForMeshes_f
!  MS FORTRAN needs explicit interface for C functions called here.
!
  INTERFACE
    INTEGER FUNCTION xfgetgrouppathsformeshes_f(a_FileId, a_Num, a_Maxsize, a_Paths, &
                                            pathlen)
	  !MS$ATTRIBUTES C,reference::xfgetgrouppathsformeshes_f
    !DEC$ATTRIBUTES reference :: a_Paths
    INTEGER,          INTENT(IN)      :: a_FileId
    INTEGER,          INTENT(IN)      :: a_Num
    INTEGER,          INTENT(IN)      :: a_Maxsize
    CHARACTER,        DIMENSION(*)    :: a_Paths
    INTEGER,          INTENT(IN)      :: pathlen

    END FUNCTION xfgetgrouppathsformeshes_f
  END INTERFACE

  pathlen = LEN(a_Paths)
  error = xfgetgrouppathsformeshes_f(a_FileId, a_Num, a_Maxsize, a_Paths, pathlen)

return

END SUBROUTINE

!**************************

SUBROUTINE XF_GRP_PTHS_SZ_FOR_GRIDS (a_FileId, a_Number, a_Maxsize, error)
INTEGER, INTENT(IN) :: a_FileId
INTEGER                       a_Number
INTEGER, INTENT(OUT)       :: a_Maxsize
INTEGER, INTENT(OUT)       :: error

!            INTEGER, EXTERNAL :: _xfGetGroupPathsSizeForGrids_f
!  MS FORTRAN needs explicit interface for C functions called here.
!
  INTERFACE
    INTEGER FUNCTION xfgetgrouppathssizeforgrids_f(a_FileId, a_Number, a_Maxsize)
	!MS$ATTRIBUTES C,reference::xfgetgrouppathssizeforgrids_f
    INTEGER,          INTENT(IN)      :: a_FileId
    INTEGER,          INTENT(IN)      :: a_Number
    INTEGER,          INTENT(OUT)     :: a_Maxsize

    END FUNCTION xfgetgrouppathssizeforgrids_f
  END INTERFACE

  error = xfgetgrouppathssizeforgrids_f(a_FileId, a_Number, a_Maxsize)

return

END SUBROUTINE

!**************************

SUBROUTINE XF_GET_GROUP_PATHS_FOR_GRIDS(a_FileId, a_Num, a_Maxsize, a_Paths, error)
INTEGER, INTENT(IN)      :: a_FileId
INTEGER, INTENT(IN)             :: a_Num
INTEGER, INTENT(IN)             :: a_Maxsize
CHARACTER,  DIMENSION(*)        :: a_Paths
INTEGER, INTENT(OUT)            :: error
INTEGER                         :: pathlen

!            INTEGER, EXTERNAL :: _xfGetGroupPathsForGrids_f
!  MS FORTRAN needs explicit interface for C functions called here.
!
  INTERFACE
    INTEGER FUNCTION xfgetgrouppathsforgrids_f(a_FileId, a_Num, a_Maxsize, a_Paths, &
                                               pathlen)
	  !MS$ATTRIBUTES C,reference::xfgetgrouppathsforgrids_f
    !DEC$ATTRIBUTES reference :: a_Paths
    INTEGER,          INTENT(IN)      :: a_FileId
    INTEGER,          INTENT(IN)      :: a_Num
    CHARACTER,        DIMENSION(*)    :: a_Paths
    INTEGER,          INTENT(IN)      :: a_Maxsize
    INTEGER,          INTENT(IN)      :: pathlen

    END FUNCTION xfgetgrouppathsforgrids_f
  END INTERFACE

  pathlen = LEN(a_Paths)
  error = xfgetgrouppathsforgrids_f (a_FileId, a_Num, a_Maxsize, a_Paths, pathlen)


return

END SUBROUTINE

!**************************

SUBROUTINE XF_GRP_PTHS_SZ_FOR_XSECS (a_FileId, a_Num, a_Maxsize, error)
INTEGER, INTENT(IN)      :: a_FileId
INTEGER, INTENT(IN)             :: a_Num
INTEGER, INTENT(IN)             :: a_Maxsize
INTEGER, INTENT(OUT)            :: error

!            INTEGER, EXTERNAL :: _xfGetGroupPathsSizeForXsecs_f
!  MS FORTRAN needs explicit interface for C functions called here.
!
  INTERFACE
    INTEGER FUNCTION xfgetgrouppathssizeforxsecs_f(a_FileId, a_Num, a_Maxsize)
	  !MS$ATTRIBUTES C,reference::xfgetgrouppathssizeforxsecs_f
    INTEGER,          INTENT(IN)      :: a_FileId
    INTEGER,          INTENT(IN)      :: a_Num
    INTEGER,          INTENT(IN)      :: a_Maxsize

    END FUNCTION xfgetgrouppathssizeforxsecs_f
  END INTERFACE

  error = xfgetgrouppathssizeforxsecs_f(a_FileId, a_Num, a_Maxsize)


return

END SUBROUTINE

!**************************

SUBROUTINE XF_GET_GROUP_PATHS_FOR_XSECS (a_FileId, a_Num, a_Maxsize, a_Paths, error)
INTEGER, INTENT(IN)      :: a_FileId
INTEGER, INTENT(IN)             :: a_Num
INTEGER, INTENT(IN)             :: a_Maxsize
CHARACTER,  DIMENSION(*)        :: a_Paths
INTEGER, INTENT(OUT)            :: error
INTEGER                         :: pathlen

!            INTEGER, EXTERNAL :: _xfGetGroupPathsForXsecs_f
!  MS FORTRAN needs explicit interface for C functions called here.
!
  INTERFACE
    INTEGER FUNCTION xfgetgrouppathsforxsecs_f(a_FileId, a_Num, a_Maxsize, a_Paths, &
                                               pathlen)
	  !MS$ATTRIBUTES C,reference::xfgetgrouppathsforxsecs_f
    !DEC$ATTRIBUTES reference :: a_Paths
    INTEGER,          INTENT(IN)      :: a_FileId
    INTEGER,          INTENT(IN)      :: a_Num
    CHARACTER,        DIMENSION(*)    :: a_Paths
    INTEGER,          INTENT(IN)      :: a_Maxsize
    INTEGER,          INTENT(IN)      :: pathlen

    END FUNCTION xfgetgrouppathsforxsecs_f
  END INTERFACE

  pathlen = LEN(a_Paths)
  error = xfgetgrouppathsforxsecs_f (a_FileId, a_Num, a_Maxsize, a_Paths, pathlen)


return

END SUBROUTINE

!**************************

SUBROUTINE XF_GRP_PTHS_SZ_FOR_GEOM_PTHS (a_FileId, a_Number, &
                                                   a_Maxsize, error)
INTEGER, INTENT(IN) :: a_FileId
INTEGER                       a_Number
INTEGER, INTENT(OUT)       :: a_Maxsize
INTEGER, INTENT(OUT)       :: error

!            INTEGER, EXTERNAL :: _xfGetPathSizeForGeomPaths_f
!  MS FORTRAN needs explicit interface for C functions called here.
!
  INTERFACE
    INTEGER FUNCTION xfgetpathsizeforgeompaths_f(a_FileId, a_Number, a_Maxsize)
	  !MS$ATTRIBUTES C,reference::xfgetpathsizeforgeompaths_f
    INTEGER,          INTENT(IN)      :: a_FileId
    INTEGER,          INTENT(IN)      :: a_Number
    INTEGER,          INTENT(OUT)     :: a_Maxsize

    END FUNCTION xfgetpathsizeforgeompaths_f
  END INTERFACE

  error = xfgetpathsizeforgeompaths_f (a_FileId, a_Number, a_Maxsize)

return

END SUBROUTINE

!**************************

SUBROUTINE XF_GRP_PTHS_FOR_GEOM_PTHS(a_FileId, a_Num, a_Maxsize, &
                                             a_Paths, error)
INTEGER, INTENT(IN)      :: a_FileId
INTEGER, INTENT(IN)             :: a_Num
INTEGER, INTENT(IN)             :: a_Maxsize
CHARACTER,  DIMENSION(*)        :: a_Paths
INTEGER, INTENT(OUT)            :: error
INTEGER                            pathlen

!            INTEGER, EXTERNAL :: _xfGetGroupPathsForGeomPaths_f
!  MS FORTRAN needs explicit interface for C functions called here.
!
  INTERFACE
    INTEGER FUNCTION xfgetgrouppathsforgeompaths_f(a_FileId, a_Num, a_Maxsize, a_Paths, &
                                                   pathlen)
	  !MS$ATTRIBUTES C,reference::xfgetgrouppathsforgeompaths_f
    !DEC$ATTRIBUTES reference :: a_Paths
    INTEGER,          INTENT(IN)      :: a_FileId
    INTEGER,          INTENT(IN)      :: a_Num
    CHARACTER,        DIMENSION(*)    :: a_Paths
    INTEGER,          INTENT(IN)      :: a_Maxsize
    INTEGER,          INTENT(IN)      :: pathlen

    END FUNCTION xfgetgrouppathsforgeompaths_f
  END INTERFACE

  pathlen = LEN(a_Paths)
  error = xfgetgrouppathsforgeompaths_f (a_FileId, a_Num, a_Maxsize, a_Paths, pathlen)

return

END SUBROUTINE

!**************************

SUBROUTINE XF_SET_NUMBER_OF_ELEMENTS (a_Id, a_nElems, error)
INTEGER, INTENT(IN) :: a_Id
INTEGER, INTENT(IN)        :: a_nElems
INTEGER, INTENT(OUT)       :: error

!            INTEGER, EXTERNAL :: _xfSetNumberOfElements_f
!  MS FORTRAN needs explicit interface for C functions called here.
!
  INTERFACE
    INTEGER FUNCTION xfsetnumberofelements_f (a_Id, a_nElems)
	  !MS$ATTRIBUTES C,reference::xfsetnumberofelements_f
    INTEGER,          INTENT(IN)      :: a_Id
    INTEGER,          INTENT(IN)      :: a_nElems

    END FUNCTION xfsetnumberofelements_f
  END INTERFACE

  error = xfsetnumberofelements_f (a_Id, a_nElems)

return

END SUBROUTINE

!***********

SUBROUTINE XF_OPEN_PROPERTY_GROUP (a_ParentId, a_GroupId, error)
  INTEGER, INTENT(IN)  :: a_ParentId
  INTEGER, INTENT(OUT) :: a_GroupId
  INTEGER, INTENT(OUT)        :: error

!            INTEGER, EXTERNAL :: _xfOpenPropertyGroup_f
!  MS FORTRAN needs explicit interface for C functions called here.
!
  INTERFACE
    INTEGER FUNCTION xfopenpropertygroup_f (a_ParentId, a_GroupId)
	  !MS$ATTRIBUTES C,reference::xfopenpropertygroup_f
    INTEGER,          INTENT(IN)      :: a_ParentId
    INTEGER,          INTENT(OUT)     :: a_GroupId

    END FUNCTION xfopenpropertygroup_f
  END INTERFACE

  error = xfopenpropertygroup_f (a_ParentId, a_GroupId)

  return 
END SUBROUTINE ! xfopenpropertygroup*/

!------------------------------------------------------------------------------
! SUBROUTINE  XF_GET_GROUP_ABSOLUTE_PATH_SIZE
! PURPOSE     returns the size of the absolute path of the group with a given
!             ID.
! NOTES
!------------------------------------------------------------------------------
SUBROUTINE XF_GET_GROUP_ABSOLUTE_PATH_SIZE (a_GroupId, a_PathLength, error)
INTEGER, INTENT(IN)    :: a_GroupId
INTEGER , INTENT(OUT)  :: a_PathLength
INTEGER, INTENT(OUT)   :: error

!            INTEGER, EXTERNAL :: _xfGetGroupAbsolutePathSize_f
!  MS FORTRAN needs explicit interface for C functions called here.
!
  INTERFACE
    INTEGER FUNCTION xfgetgroupabsolutepathsize_f (a_GroupId, a_PathLength)
	  !MS$ATTRIBUTES C,reference::xfgetgroupabsolutepathsize_f
    INTEGER,          INTENT(IN)      :: a_GroupId
    INTEGER,          INTENT(OUT)     :: a_PathLength

    END FUNCTION xfgetgroupabsolutepathsize_f
  END INTERFACE

  error = xfgetgroupabsolutepathsize_f (a_GroupId, a_PathLength)

  return

END SUBROUTINE

!------------------------------------------------------------------------------
! SUBROUTINE  XF_GET_GROUP_ABSOLUTE_PATH
! PURPOSE     returns the absolute path of the group with a given ID.
! NOTES
!------------------------------------------------------------------------------
SUBROUTINE XF_GET_GROUP_ABSOLUTE_PATH (a_GroupId, a_PathLength, a_Path, error)
INTEGER, INTENT(IN)    :: a_GroupId
INTEGER, INTENT(IN)    :: a_PathLength
CHARACTER(LEN=*)              :: a_Path
INTEGER, INTENT(OUT)          :: error
INTEGER                       :: pathlen

!            INTEGER, EXTERNAL :: _xfGetGroupAbsolutePath_f
!  MS FORTRAN needs explicit interface for C functions called here.
!
  INTERFACE
    INTEGER FUNCTION xfgetgroupabsolutepath_f(a_GroupId, a_PathLength, a_Path, pathlen)
	  !MS$ATTRIBUTES C,reference::xfgetgroupabsolutepath_f
    !DEC$ATTRIBUTES reference :: a_Path
    INTEGER,          INTENT(IN)      :: a_GroupId
    INTEGER,          INTENT(IN)      :: a_PathLength
    CHARACTER(LEN=*)                  :: a_Path
    INTEGER,          INTENT(IN)      :: pathlen

    END FUNCTION xfgetgroupabsolutepath_f
  END INTERFACE

  pathlen = LEN_TRIM(a_Path)
  error = xfgetgroupabsolutepath_f(a_GroupId, a_PathLength, a_Path, pathlen)

return

END SUBROUTINE

!**************************

SUBROUTINE XF_SET_ALL_ELEMS_SAME_TYPE (a_Id, a_Type, error)
INTEGER, INTENT(IN) :: a_Id
INTEGER, INTENT(IN)        :: a_Type
INTEGER, INTENT(OUT)       :: error


!            INTEGER, EXTERNAL :: _xfSetAllElemsSameType_f
!  MS FORTRAN needs explicit interface for C functions called here.
!
  INTERFACE
    INTEGER FUNCTION xfsetallelemssametype_f(a_Id, a_Type)
	  !MS$ATTRIBUTES C,reference::xfsetallelemssametype_f
    INTEGER,          INTENT(IN)      :: a_Id
    INTEGER,          INTENT(IN)      :: a_Type

    END FUNCTION xfsetallelemssametype_f
  END INTERFACE

  error = xfsetallelemssametype_f (a_Id, a_Type)

return

END SUBROUTINE

!**************************

SUBROUTINE XF_WRITE_ELEM_TYPES (a_Id, a_nElems, a_Type, a_Compression, error)
INTEGER,         INTENT(IN)  :: a_Id
INTEGER,                INTENT(IN)  :: a_nElems, a_Compression
INTEGER, DIMENSION(*),  INTENT(IN)  :: a_Type
INTEGER,                INTENT(OUT) :: error

!            INTEGER, EXTERNAL :: _xfWriteElemTypes_f
!  MS FORTRAN needs explicit interface for C functions called here.
!
  INTERFACE
    INTEGER FUNCTION xfwriteelemtypes_f(a_Id, a_nElems, a_Type, a_Compression)
	  !MS$ATTRIBUTES C,reference::xfwriteelemtypes_f
    INTEGER,                 INTENT(IN)      :: a_Id
    INTEGER,                 INTENT(IN)      :: a_nElems, a_Compression
    INTEGER,  DIMENSION(*),  INTENT(IN)      :: a_Type

    END FUNCTION xfwriteelemtypes_f
  END INTERFACE

  error = xfwriteelemtypes_f (a_Id, a_nElems, a_Type, a_Compression)


return

END SUBROUTINE

!******************************

SUBROUTINE XF_WRITE_ELEM_NODE_IDS (a_Id, a_nElems, a_nMaxNodes, a_Ids, a_Compression, error)
INTEGER,         INTENT(IN)        :: a_Id
INTEGER,                INTENT(IN)        :: a_nElems, a_nMaxNodes, a_Compression
INTEGER, DIMENSION(*),  INTENT(IN)        :: a_Ids
INTEGER,                INTENT(OUT)       :: error

!            INTEGER, EXTERNAL :: _xfWriteElemNodeIds_f
!  MS FORTRAN needs explicit interface for C functions called here.
!
  INTERFACE
    INTEGER FUNCTION xfwriteelemnodeids_f(a_Id, a_nElems, a_nMaxNodes, a_Ids, a_Compression)
	  !MS$ATTRIBUTES C,reference::xfwriteelemnodeids_f
    INTEGER,              INTENT(IN)      :: a_Id
    INTEGER,              INTENT(IN)      :: a_nElems, a_nMaxNodes, a_Compression
    INTEGER, DIMENSION(*),INTENT(IN)      :: a_Ids

    END FUNCTION xfwriteelemnodeids_f
  END INTERFACE

  error = xfwriteelemnodeids_f (a_Id, a_nElems, a_nMaxNodes, a_Ids, a_Compression)

return

END SUBROUTINE

!******************************

SUBROUTINE XF_SET_NUMBER_OF_NODES (a_Id, a_nNodes, error)
INTEGER, INTENT(IN) :: a_Id
INTEGER, INTENT(IN)        :: a_nNodes
INTEGER, INTENT(OUT)       :: error


!            INTEGER, EXTERNAL :: _xfSetNumberOfNodes_f
!  MS FORTRAN needs explicit interface for C functions called here.
!
  INTERFACE
    INTEGER FUNCTION xfsetnumberofnodes_f(a_Id, a_nNodes)
	  !MS$ATTRIBUTES C,reference::xfsetnumberofnodes_f
    INTEGER,          INTENT(IN)      :: a_Id
    INTEGER,          INTENT(IN)      :: a_nNodes

    END FUNCTION xfsetnumberofnodes_f
  END INTERFACE

  error = xfsetnumberofnodes_f (a_Id, a_nNodes)

return

END SUBROUTINE

!******************************

SUBROUTINE XF_WRITE_X_NODE_LOCATIONS (a_Id, a_nNodes, &
                                      a_Locs, a_Compression, error)
INTEGER, INTENT(IN)          :: a_Id
INTEGER, INTENT(IN)                 :: a_nNodes, a_Compression
REAL*8, DIMENSION(*), INTENT(OUT)   :: a_Locs
INTEGER, INTENT(OUT)                :: error


!            INTEGER, EXTERNAL :: _xfWriteXNodeLocations_f
!  MS FORTRAN needs explicit interface for C functions called here.
!
  INTERFACE
    INTEGER FUNCTION xfwritexnodelocations_f(a_Id, a_nNodes, a_Locs, a_Compression)
	  !MS$ATTRIBUTES C,reference::xfwritexnodelocations_f
    INTEGER,              INTENT(IN)    :: a_Id
    INTEGER,              INTENT(IN)    :: a_nNodes, a_Compression
    REAL*8, DIMENSION(*), INTENT(OUT)   :: a_Locs


    END FUNCTION xfwritexnodelocations_f
  END INTERFACE

  error = xfwritexnodelocations_f (a_Id, a_nNodes, a_Locs, a_Compression)

return

END SUBROUTINE

!******************************

SUBROUTINE XF_WRITE_Y_NODE_LOCATIONS (a_Id, a_nNodes, a_Locs, error)
INTEGER, INTENT(IN) :: a_Id
INTEGER, INTENT(IN)        :: a_nNodes
REAL*8, DIMENSION(*), INTENT(OUT)  :: a_Locs
INTEGER, INTENT(OUT)       :: error


!            INTEGER, EXTERNAL :: _xfWriteYNodeLocations_f
!  MS FORTRAN needs explicit interface for C functions called here.
!
  INTERFACE
    INTEGER FUNCTION xfwriteynodelocations_f(a_Id, a_nNodes, a_Locs)
	  !MS$ATTRIBUTES C,reference::xfwriteynodelocations_f
    INTEGER,              INTENT(IN)    :: a_Id
    INTEGER,              INTENT(IN)    :: a_nNodes
    REAL*8, DIMENSION(*), INTENT(OUT)   :: a_Locs


    END FUNCTION xfwriteynodelocations_f
  END INTERFACE

  error = xfwriteynodelocations_f (a_Id, a_nNodes, a_Locs)

return

END SUBROUTINE

!******************************

SUBROUTINE XF_WRITE_Z_NODE_LOCATIONS (a_Id, a_nNodes, a_Locs, error)
INTEGER, INTENT(IN) :: a_Id
INTEGER, INTENT(IN)        :: a_nNodes
REAL*8, DIMENSION(*), INTENT(OUT)  :: a_Locs
INTEGER, INTENT(OUT)       :: error


!            INTEGER, EXTERNAL :: _xfWriteZNodeLocations_f
!  MS FORTRAN needs explicit interface for C functions called here.
!
  INTERFACE
    INTEGER FUNCTION xfwriteznodelocations_f(a_Id, a_nNodes, a_Locs)
	  !MS$ATTRIBUTES C,reference::xfwriteznodelocations_f
    INTEGER,              INTENT(IN)    :: a_Id
    INTEGER,              INTENT(IN)    :: a_nNodes
    REAL*8, DIMENSION(*), INTENT(OUT)   :: a_Locs


    END FUNCTION xfwriteznodelocations_f
  END INTERFACE

  error = xfwriteznodelocations_f (a_Id, a_nNodes, a_Locs)

return

END SUBROUTINE

!******************************

SUBROUTINE XF_GET_NUMBER_OF_ELEMENTS (a_Id, a_nElems, error)
INTEGER, INTENT(IN) :: a_Id
INTEGER, INTENT(OUT)       :: a_nElems, error


!            INTEGER, EXTERNAL :: _xfGetNumberOfElements_f
!  MS FORTRAN needs explicit interface for C functions called here.
!
  INTERFACE
    INTEGER FUNCTION xfgetnumberofelements_f(a_Id, a_nElems)
	  !MS$ATTRIBUTES C,reference::xfgetnumberofelements_f
    INTEGER,              INTENT(IN)    :: a_Id
    INTEGER,              INTENT(OUT)   :: a_nElems


    END FUNCTION xfgetnumberofelements_f
  END INTERFACE

  error = xfgetnumberofelements_f (a_Id, a_nElems)

return

END SUBROUTINE

!******************************

SUBROUTINE XF_ARE_ALL_ELEMS_SAME_TYPE (a_Id, a_Same, error)
INTEGER, INTENT(IN)        :: a_Id
LOGICAL*2, INTENT(OUT)       :: a_Same
INTEGER, INTENT(OUT)       :: error


!            INTEGER, EXTERNAL :: _xfAreAllElemsSameType_f
!  MS FORTRAN needs explicit interface for C functions called here.
!
  INTERFACE
    INTEGER FUNCTION xfareallelemssametype_f(a_Id, a_Same)
	  !MS$ATTRIBUTES C,reference::xfareallelemssametype_f
    INTEGER,       INTENT(IN)    :: a_Id
    LOGICAL*2,       INTENT(OUT)   :: a_Same


    END FUNCTION xfareallelemssametype_f
  END INTERFACE

  error = xfareallelemssametype_f (a_Id, a_Same)

return

END SUBROUTINE

!******************************

SUBROUTINE XF_READ_ELEM_TYPES_SINGLE_VALUE (a_Id, a_Type, error)
INTEGER, INTENT(IN) :: a_Id
INTEGER, INTENT(OUT)       :: a_Type
INTEGER, INTENT(OUT)       :: error


!            INTEGER, EXTERNAL :: _xfReadElemTypesSingleValue_f
!  MS FORTRAN needs explicit interface for C functions called here.
!
  INTERFACE
    INTEGER FUNCTION xfreadelemtypessinglevalue_f(a_Id, a_Type)
	  !MS$ATTRIBUTES C,reference::xfreadelemtypessinglevalue_f
    INTEGER,       INTENT(IN)    :: a_Id
    INTEGER,       INTENT(OUT)   :: a_Type


    END FUNCTION xfreadelemtypessinglevalue_f
  END INTERFACE

  error = xfreadelemtypessinglevalue_f (a_Id, a_Type)


return

END SUBROUTINE

!******************************

SUBROUTINE XF_READ_ELEM_TYPES (a_Id, a_nElems, a_Type, error)
INTEGER,         INTENT(IN)    :: a_Id
INTEGER,                INTENT(IN)    :: a_nElems
INTEGER, DIMENSION(*),  INTENT(OUT)   :: a_Type
INTEGER,                INTENT(OUT)   :: error


!            INTEGER, EXTERNAL :: _xfReadElemTypes_f
!  MS FORTRAN needs explicit interface for C functions called here.
!
  INTERFACE
    INTEGER FUNCTION xfreadelemtypes_f(a_Id, a_nElems, a_Type)
	  !MS$ATTRIBUTES C,reference::xfreadelemtypes_f
    INTEGER,               INTENT(IN)    :: a_Id
    INTEGER,               INTENT(IN)    :: a_nElems
    INTEGER, DIMENSION(*), INTENT(OUT)   :: a_Type


    END FUNCTION xfreadelemtypes_f
  END INTERFACE

  error = xfreadelemtypes_f (a_Id, a_nElems, a_Type)


return

END SUBROUTINE

!******************************

SUBROUTINE XF_GET_MAX_NODES_IN_ELEM (a_Id, a_nMaxNodes, error)
INTEGER, INTENT(IN) :: a_Id
INTEGER, INTENT(OUT)       :: a_nMaxNodes, error


!            INTEGER, EXTERNAL :: _xfGetMaxNodesInElem_f
!  MS FORTRAN needs explicit interface for C functions called here.
!
  INTERFACE
    INTEGER FUNCTION xfgetmaxnodesinelem_f(a_Id, a_nMaxNodes)
	  !MS$ATTRIBUTES C,reference::xfgetmaxnodesinelem_f
    INTEGER,       INTENT(IN)    :: a_Id
    INTEGER,       INTENT(OUT)   :: a_nMaxNodes


    END FUNCTION xfgetmaxnodesinelem_f
  END INTERFACE

  error = xfgetmaxnodesinelem_f (a_Id, a_nMaxNodes)


return

END SUBROUTINE

!******************************

SUBROUTINE XF_READ_ELEM_NODE_IDS (a_Id, a_nElems, a_nMaxNodes, a_Ids, error)
INTEGER, INTENT(IN) :: a_Id
INTEGER, INTENT(IN)        :: a_nElems, a_nMaxNodes
INTEGER, INTENT(OUT)       :: a_Ids(a_nMaxNodes,*), error


!            INTEGER, EXTERNAL :: _xfReadElemNodeIds_f
!  MS FORTRAN needs explicit interface for C functions called here.
!
  INTERFACE
    INTEGER FUNCTION xfreadelemnodeids_f(a_Id, a_nElems, a_nMaxNodes, a_Ids)
	  !MS$ATTRIBUTES C,reference::xfreadelemnodeids_f
    INTEGER,       INTENT(IN)    :: a_Id
    INTEGER,       INTENT(IN)    :: a_nElems, a_nMaxNodes
    INTEGER,       INTENT(OUT)   :: a_Ids(a_nMaxNodes, *)


    END FUNCTION xfreadelemnodeids_f
  END INTERFACE

  error = xfreadelemnodeids_f (a_Id, a_nElems, a_nMaxNodes, a_Ids)

return

END SUBROUTINE

!******************************

SUBROUTINE XF_GET_NUMBER_OF_NODES (a_Id, a_nNodes, error)
INTEGER, INTENT(IN) :: a_Id
INTEGER, INTENT(OUT)       :: a_nNodes, error


!            INTEGER, EXTERNAL :: _xfGetNumberOfNodes_f
!  MS FORTRAN needs explicit interface for C functions called here.
!
  INTERFACE
    INTEGER FUNCTION xfgetnumberofnodes_f(a_Id, a_nNodes)
	  !MS$ATTRIBUTES C,reference::xfgetnumberofnodes_f
    INTEGER,       INTENT(IN)    :: a_Id
    INTEGER,       INTENT(OUT)   :: a_nNodes


    END FUNCTION xfgetnumberofnodes_f
  END INTERFACE

  error = xfgetnumberofnodes_f (a_Id, a_nNodes)

return

END SUBROUTINE

!******************************

SUBROUTINE XF_READ_X_NODE_LOCATIONS (a_Id, a_nNodes, a_Locs, error)
INTEGER, INTENT(IN) :: a_Id
INTEGER, INTENT(IN)        :: a_nNodes
REAL*8, DIMENSION(*), INTENT(OUT)  :: a_Locs
INTEGER, INTENT(OUT)       :: error


!            INTEGER, EXTERNAL :: _xfReadXNodeLocations_f
!  MS FORTRAN needs explicit interface for C functions called here.
!
  INTERFACE
    INTEGER FUNCTION xfreadxnodelocations_f(a_Id, a_nNodes, a_Locs)
	  !MS$ATTRIBUTES C,reference::xfreadxnodelocations_f
    INTEGER,              INTENT(IN)    :: a_Id
    INTEGER,              INTENT(IN)    :: a_nNodes
    REAL*8, DIMENSION(*), INTENT(OUT)   :: a_Locs


    END FUNCTION xfreadxnodelocations_f
  END INTERFACE

  error = xfreadxnodelocations_f (a_Id, a_nNodes, a_Locs)


return

END SUBROUTINE

!******************************

SUBROUTINE XF_READ_Y_NODE_LOCATIONS (a_Id, a_nNodes, a_Locs, error)
INTEGER, INTENT(IN) :: a_Id
INTEGER, INTENT(IN)        :: a_nNodes
REAL*8, DIMENSION(*), INTENT(OUT)  :: a_Locs
INTEGER, INTENT(OUT)       :: error


!            INTEGER, EXTERNAL :: _xfReadYNodeLocations_f
!  MS FORTRAN needs explicit interface for C functions called here.
!
  INTERFACE
    INTEGER FUNCTION xfreadynodelocations_f(a_Id, a_nNodes, a_Locs)
	  !MS$ATTRIBUTES C,reference::xfreadynodelocations_f
    INTEGER,              INTENT(IN)    :: a_Id
    INTEGER,              INTENT(IN)    :: a_nNodes
    REAL*8, DIMENSION(*), INTENT(OUT)   :: a_Locs


    END FUNCTION xfreadynodelocations_f
  END INTERFACE

  error = xfreadynodelocations_f (a_Id, a_nNodes, a_Locs)


return

END SUBROUTINE

!******************************

SUBROUTINE XF_READ_Z_NODE_LOCATIONS (a_Id, a_nNodes, a_Locs, error)
INTEGER, INTENT(IN) :: a_Id
INTEGER, INTENT(IN)        :: a_nNodes
REAL*8, DIMENSION(*), INTENT(OUT)  :: a_Locs
INTEGER, INTENT(OUT)       :: error


!            INTEGER, EXTERNAL :: _xfReadZNodeLocations_f
!  MS FORTRAN needs explicit interface for C functions called here.
!
  INTERFACE
    INTEGER FUNCTION xfreadznodelocations_f(a_Id, a_nNodes, a_Locs)
	  !MS$ATTRIBUTES C,reference::xfreadznodelocations_f
    INTEGER,              INTENT(IN)    :: a_Id
    INTEGER,              INTENT(IN)    :: a_nNodes
    REAL*8, DIMENSION(*), INTENT(OUT)   :: a_Locs


    END FUNCTION xfreadznodelocations_f
  END INTERFACE

  error = xfreadznodelocations_f (a_Id, a_nNodes, a_Locs)


return

END SUBROUTINE

!******************************

SUBROUTINE XF_CREATE_MESH_PROPERTY_GROUP (a_Id, a_PropId, error)
INTEGER, INTENT(IN)  :: a_Id
INTEGER, INTENT(OUT) :: a_PropId
INTEGER, INTENT(OUT)        :: error


!            INTEGER, EXTERNAL :: _xfCreateMeshPropertyGroup_f
!  MS FORTRAN needs explicit interface for C functions called here.
!
  INTERFACE
    INTEGER FUNCTION xfcreatemeshpropertygroup_f(a_Id, a_PropId)
	  !MS$ATTRIBUTES C,reference::xfcreatemeshpropertygroup_f
    INTEGER,              INTENT(IN)    :: a_Id
    INTEGER,              INTENT(OUT)   :: a_PropId


    END FUNCTION xfcreatemeshpropertygroup_f
  END INTERFACE

  error = xfcreatemeshpropertygroup_f (a_Id, a_PropId)


return

END SUBROUTINE

!******************************

SUBROUTINE XF_GET_MESH_PROPERTY_GROUP (a_Id, a_PropId, error)
INTEGER, INTENT(IN)  :: a_Id
INTEGER, INTENT(OUT) :: a_PropId
INTEGER, INTENT(OUT)        :: error


!            INTEGER, EXTERNAL :: _xfGetMeshPropertyGroup_f
!  MS FORTRAN needs explicit interface for C functions called here.
!
  INTERFACE
    INTEGER FUNCTION xfgetmeshpropertygroup_f(a_Id, a_PropId)
	  !MS$ATTRIBUTES C,reference::xfgetmeshpropertygroup_f
    INTEGER,              INTENT(IN)    :: a_Id
    INTEGER,              INTENT(OUT)   :: a_PropId


    END FUNCTION xfgetmeshpropertygroup_f
  END INTERFACE

  error = xfgetmeshpropertygroup_f (a_Id, a_PropId)


return

END SUBROUTINE

!******************************

SUBROUTINE XF_CREATE_MESH_NODE_PROP_GRP (a_Id, a_PropId, error)
INTEGER, INTENT(IN)  :: a_Id
INTEGER, INTENT(OUT) :: a_PropId
INTEGER, INTENT(OUT)        :: error


!            INTEGER, EXTERNAL :: _xfCreateMeshNodePropertyGroup_f
!  MS FORTRAN needs explicit interface for C functions called here.
!
  INTERFACE
    INTEGER FUNCTION xfcreatemeshnodepropertygroup_f(a_Id, a_PropId)
	  !MS$ATTRIBUTES C,reference::xfcreatemeshnodepropertygroup_f
    INTEGER,              INTENT(IN)    :: a_Id
    INTEGER,              INTENT(OUT)   :: a_PropId


    END FUNCTION xfcreatemeshnodepropertygroup_f
  END INTERFACE

  error = xfcreatemeshnodepropertygroup_f (a_Id, a_PropId)


return

END SUBROUTINE

!******************************

SUBROUTINE XF_GET_MESH_NODE_PROPERTY_GROUP (a_Id, a_PropId, error)
INTEGER, INTENT(IN)  :: a_Id
INTEGER, INTENT(OUT) :: a_PropId
INTEGER, INTENT(OUT)        :: error


!            INTEGER, EXTERNAL :: _xfGetMeshNodePropertyGroup_f
!  MS FORTRAN needs explicit interface for C functions called here.
!
  INTERFACE
    INTEGER FUNCTION xfgetmeshnodepropertygroup_f(a_Id, a_PropId)
	  !MS$ATTRIBUTES C,reference::xfgetmeshnodepropertygroup_f
    INTEGER,              INTENT(IN)    :: a_Id
    INTEGER,              INTENT(OUT)   :: a_PropId


    END FUNCTION xfgetmeshnodepropertygroup_f
  END INTERFACE

  error = xfgetmeshnodepropertygroup_f (a_Id, a_PropId)


return

END SUBROUTINE

!******************************

SUBROUTINE XF_CREATE_MESH_ELEM_PROP_GRP (a_Id, a_PropId, error)
INTEGER, INTENT(IN)  :: a_Id
INTEGER, INTENT(OUT) :: a_PropId
INTEGER, INTENT(OUT)        :: error


!            INTEGER, EXTERNAL :: _xfCreateMeshElemPropGrp_f
!  MS FORTRAN needs explicit interface for C functions called here.
!
  INTERFACE
    INTEGER FUNCTION xfcreatemeshelempropgrp_f(a_Id, a_PropId)
	  !MS$ATTRIBUTES C,reference::xfcreatemeshelempropgrp_f
    INTEGER,              INTENT(IN)    :: a_Id
    INTEGER,              INTENT(OUT)   :: a_PropId


    END FUNCTION xfcreatemeshelempropgrp_f
  END INTERFACE

  error = xfcreatemeshelempropgrp_f (a_Id, a_PropId)


return

END SUBROUTINE

!******************************

SUBROUTINE XF_GET_MESH_ELEM_PROP_GRP (a_Id, a_PropId, error)
INTEGER, INTENT(IN)  :: a_Id
INTEGER, INTENT(OUT) :: a_PropId
INTEGER, INTENT(OUT)        :: error


!            INTEGER, EXTERNAL :: _xfGetMeshElementPropertyGroup_f
!  MS FORTRAN needs explicit interface for C functions called here.
!
  INTERFACE
    INTEGER FUNCTION xfgetmeshelementpropertygroup_f(a_Id, a_PropId)
	  !MS$ATTRIBUTES C,reference::xfgetmeshelementpropertygroup_f
    INTEGER,              INTENT(IN)    :: a_Id
    INTEGER,              INTENT(OUT)   :: a_PropId


    END FUNCTION xfgetmeshelementpropertygroup_f
  END INTERFACE

  error = xfgetmeshelementpropertygroup_f (a_Id, a_PropId)


return

END SUBROUTINE

!******************************

SUBROUTINE XF_SET_GRID_TYPE (a_Id, a_GridType, error)
INTEGER, INTENT(IN) :: a_Id
INTEGER, INTENT(IN)        :: a_GridType
INTEGER, INTENT(OUT)       :: error


!            INTEGER, EXTERNAL :: _xfSetGridType_f
!  MS FORTRAN needs explicit interface for C functions called here.
!
  INTERFACE
    INTEGER FUNCTION xfsetgridtype_f(a_Id, a_GridType)
	  !MS$ATTRIBUTES C,reference::xfsetgridtype_f
    INTEGER,              INTENT(IN)   :: a_Id
    INTEGER,              INTENT(IN)   :: a_GridType


    END FUNCTION xfsetgridtype_f
  END INTERFACE

  error = xfsetgridtype_f (a_Id, a_GridType)


return

END SUBROUTINE

!******************************

SUBROUTINE XF_SET_NUMBER_OF_DIMENSIONS (a_Id, a_NumDimensions, error)
INTEGER, INTENT(IN) :: a_Id
INTEGER, INTENT(IN)        :: a_NumDimensions
INTEGER, INTENT(OUT)       :: error

!            INTEGER, EXTERNAL :: _xfSetNumberOfDimensions_f
!  MS FORTRAN needs explicit interface for C functions called here.
!
  INTERFACE
    INTEGER FUNCTION xfsetnumberofdimensions_f(a_Id, a_NumDimensions)
	  !MS$ATTRIBUTES C,reference::xfsetnumberofdimensions_f
    INTEGER,              INTENT(IN)   :: a_Id
    INTEGER,              INTENT(IN)   :: a_NumDimensions


    END FUNCTION xfsetnumberofdimensions_f
  END INTERFACE

  error = xfsetnumberofdimensions_f (a_Id, a_NumDimensions)


return

END SUBROUTINE

!******************************

SUBROUTINE XF_SET_EXTRUSION_TYPE (a_Id, a_ExtrudeType, error)
INTEGER, INTENT(IN) :: a_Id
INTEGER, INTENT(IN)        :: a_ExtrudeType
INTEGER, INTENT(OUT)       :: error

!            INTEGER, EXTERNAL :: _xfSetExtrusionType_f
!  MS FORTRAN needs explicit interface for C functions called here.
!
  INTERFACE
    INTEGER FUNCTION xfsetextrusiontype_f(a_Id, a_ExtrudeType)
	  !MS$ATTRIBUTES C,reference::xfsetextrusiontype_f
    INTEGER,              INTENT(IN)   :: a_Id
    INTEGER,              INTENT(IN)   :: a_ExtrudeType


    END FUNCTION xfsetextrusiontype_f
  END INTERFACE

  error = xfsetextrusiontype_f (a_Id, a_ExtrudeType)

return

END SUBROUTINE

!******************************

SUBROUTINE XF_SET_ORIGIN (a_Id, a_x, a_y, a_z, error)
INTEGER, INTENT(IN) :: a_Id
REAL*8,  INTENT(IN)   :: a_x
REAL*8,  INTENT(IN)   :: a_y
REAL*8,  INTENT(IN)   :: a_z
INTEGER, INTENT(OUT)       :: error

!            INTEGER, EXTERNAL :: _xfSetOrigin_f
!  MS FORTRAN needs explicit interface for C functions called here.
!
  INTERFACE
    INTEGER FUNCTION xfsetorigin_f(a_Id, a_x, a_y, a_z)
	  !MS$ATTRIBUTES C,reference::xfsetorigin_f
    INTEGER, INTENT(IN)    :: a_Id
    REAL*8,  INTENT(IN)    :: a_x
    REAL*8,  INTENT(IN)    :: a_y
    REAL*8,  INTENT(IN)    :: a_z



    END FUNCTION xfsetorigin_f
  END INTERFACE

  error = xfsetorigin_f (a_Id, a_x, a_y, a_z)

return

END SUBROUTINE

!******************************

SUBROUTINE XF_SET_ORIENTATION (a_Id, a_Orientation, error)
INTEGER, INTENT(IN) :: a_Id
INTEGER, INTENT(IN)        :: a_Orientation
INTEGER, INTENT(OUT)       :: error

!            INTEGER, EXTERNAL :: _xfSetOrientation_f
!  MS FORTRAN needs explicit interface for C functions called here.
!
  INTERFACE
    INTEGER FUNCTION xfsetorientation_f(a_Id, a_Orientation)
	  !MS$ATTRIBUTES C,reference::xfsetorientation_f
    INTEGER,              INTENT(IN)    :: a_Id
    INTEGER,              INTENT(IN)    :: a_Orientation



    END FUNCTION xfsetorientation_f
  END INTERFACE

  error = xfsetorientation_f (a_Id, a_Orientation)

return

END SUBROUTINE

!******************************

SUBROUTINE XF_SET_BEARING (a_Id, a_Bearing, error)
INTEGER, INTENT(IN)   :: a_Id
REAL*8,  INTENT(IN)   :: a_Bearing
INTEGER, INTENT(OUT)  :: error


!            INTEGER, EXTERNAL :: _xfSetBearing_f
!  MS FORTRAN needs explicit interface for C functions called here.
!
  INTERFACE
    INTEGER FUNCTION xfsetbearing_f(a_Id, a_Bearing)
	  !MS$ATTRIBUTES C,reference::xfsetbearing_f
    INTEGER,  INTENT(IN)    :: a_Id
    REAL*8,   INTENT(IN)    :: a_Bearing



    END FUNCTION xfsetbearing_f
  END INTERFACE

  error = xfsetbearing_f (a_Id, a_Bearing)

return

END SUBROUTINE

!******************************

SUBROUTINE XF_SET_DIP (a_Id, a_Dip, error)
INTEGER, INTENT(IN)   :: a_Id
REAL*8,  INTENT(IN)   :: a_Dip
INTEGER, INTENT(OUT)  :: error


!            INTEGER, EXTERNAL :: _xfSetDip_f
!  MS FORTRAN needs explicit interface for C functions called here.
!
  INTERFACE
    INTEGER FUNCTION xfsetdip_f(a_Id, a_Dip)
	  !MS$ATTRIBUTES C,reference::xfsetdip_f
    INTEGER, INTENT(IN)    :: a_Id
    REAL*8,  INTENT(IN)    :: a_Dip



    END FUNCTION xfsetdip_f
  END INTERFACE

  error = xfsetdip_f (a_Id, a_Dip)

return

END SUBROUTINE

!******************************

SUBROUTINE XF_SET_ROLL (a_Id, a_Roll, error)
INTEGER, INTENT(IN)    :: a_Id
REAL*8,  INTENT(IN)    :: a_Roll
INTEGER, INTENT(OUT)   :: error


!            INTEGER, EXTERNAL :: _xfSetRoll_f
!  MS FORTRAN needs explicit interface for C functions called here.
!
  INTERFACE
    INTEGER FUNCTION xfsetroll_f(a_Id, a_Roll)
	  !MS$ATTRIBUTES C,reference::xfsetroll_f
    INTEGER,  INTENT(IN)    :: a_Id
    REAL*8,   INTENT(IN)    :: a_Roll



    END FUNCTION xfsetroll_f
  END INTERFACE

  error = xfsetroll_f (a_Id, a_Roll)

return

END SUBROUTINE

!*******************************************

SUBROUTINE XF_SET_COMPUTATIONAL_ORIGIN (a_Id, a_Origin, error)
INTEGER, INTENT(IN)  :: a_Id
INTEGER, INTENT(IN)         :: a_Origin
INTEGER, INTENT(OUT)        :: error


!            INTEGER, EXTERNAL :: _xfSetComputationalOrigin_f
!  MS FORTRAN needs explicit interface for C functions called here.
!
  INTERFACE
    INTEGER FUNCTION xfsetcomputationalorigin_f(a_Id, a_Origin)
	  !MS$ATTRIBUTES C,reference::xfsetcomputationalorigin_f
    INTEGER,              INTENT(IN)    :: a_Id
    INTEGER,              INTENT(IN)    :: a_Origin



    END FUNCTION xfsetcomputationalorigin_f
  END INTERFACE

  error = xfsetcomputationalorigin_f (a_Id, a_Origin)

return

END SUBROUTINE

!********************************************

SUBROUTINE XF_SET_U_DIRECTION (a_Id, a_Direction, error)
INTEGER, INTENT(IN)  :: a_Id
INTEGER, INTENT(IN)         :: a_Direction
INTEGER                        error

!            INTEGER, EXTERNAL :: _xfSetUDirection_f
!  MS FORTRAN needs explicit interface for C functions called here.
!
  INTERFACE
    INTEGER FUNCTION xfsetudirection_f(a_Id, a_Direction)
	  !MS$ATTRIBUTES C,reference::xfsetudirection_f
    INTEGER,              INTENT(IN)    :: a_Id
    INTEGER,              INTENT(IN)    :: a_Direction



    END FUNCTION xfsetudirection_f
  END INTERFACE

  error = xfsetudirection_f (a_Id, a_Direction)

return

END SUBROUTINE

!******************************

SUBROUTINE XF_SET_NUMBER_CELLS_IN_I (a_Id, a_NumI, error)
INTEGER, INTENT(IN)  :: a_Id
INTEGER, INTENT(IN)         :: a_NumI
INTEGER, INTENT(OUT)        :: error

!            INTEGER, EXTERNAL :: _xfSetNumberCellsInI_f
!  MS FORTRAN needs explicit interface for C functions called here.
!
  INTERFACE
    INTEGER FUNCTION xfsetnumbercellsini_f(a_Id, a_NumI)
	  !MS$ATTRIBUTES C,reference::xfsetnumbercellsini_f
    INTEGER,              INTENT(IN)    :: a_Id
    INTEGER,              INTENT(IN)    :: a_NumI



    END FUNCTION xfsetnumbercellsini_f
  END INTERFACE

  error = xfsetnumbercellsini_f (a_Id, a_NumI)

return

END SUBROUTINE

!******************************

SUBROUTINE XF_SET_NUMBER_CELLS_IN_J (a_Id, a_NumJ, error)
INTEGER             , INTENT(IN)  :: a_Id
INTEGER, INTENT(IN)         :: a_NumJ
INTEGER, INTENT(OUT)        :: error


!            INTEGER, EXTERNAL :: _xfSetNumberCellsInJ_f
!  MS FORTRAN needs explicit interface for C functions called here.
!
  INTERFACE
    INTEGER FUNCTION xfsetnumbercellsinj_f(a_Id, a_NumJ)
	  !MS$ATTRIBUTES C,reference::xfsetnumbercellsinj_f
    INTEGER,              INTENT(IN)    :: a_Id
    INTEGER,              INTENT(IN)    :: a_NumJ



    END FUNCTION xfsetnumbercellsinj_f
  END INTERFACE

  error = xfsetnumbercellsinj_f (a_Id, a_NumJ)

return

END SUBROUTINE

!******************************

SUBROUTINE XF_SET_NUMBER_CELLS_IN_K (a_Id, a_NumK, error)
INTEGER             , INTENT(IN)  :: a_Id
INTEGER, INTENT(IN)         :: a_NumK
INTEGER, INTENT(OUT)        :: error

!            INTEGER, EXTERNAL :: _xfSetNumberCellsInK_f
!  MS FORTRAN needs explicit interface for C functions called here.
!
  INTERFACE
    INTEGER FUNCTION xfsetnumbercellsink_f(a_Id, a_NumK)
	  !MS$ATTRIBUTES C,reference::xfsetnumbercellsink_f
    INTEGER,              INTENT(IN)    :: a_Id
    INTEGER,              INTENT(IN)    :: a_NumK



    END FUNCTION xfsetnumbercellsink_f
  END INTERFACE

  error = xfsetnumbercellsink_f (a_Id, a_NumK)

return

END SUBROUTINE

!******************************

SUBROUTINE XF_SET_GRID_COORDS_I (a_Id, a_NumVals, a_iValues, error)
INTEGER             , INTENT(IN) :: a_Id
INTEGER, INTENT(IN)        :: a_NumVals
REAL*8, DIMENSION(*), INTENT(IN)   :: a_iValues
INTEGER, INTENT(OUT)       :: error


!            INTEGER, EXTERNAL :: _xfSetGridCoordsI_f
!  MS FORTRAN needs explicit interface for C functions called here.
!
  INTERFACE
    INTEGER FUNCTION xfsetgridcoordsi_f(a_Id, a_NumVals, a_iValues)
	  !MS$ATTRIBUTES C,reference::xfsetgridcoordsi_f
    INTEGER,              INTENT(IN)    :: a_Id
    INTEGER,              INTENT(IN)    :: a_NumVals
    REAL*8, DIMENSION(*), INTENT(IN)    :: a_iValues


    END FUNCTION xfsetgridcoordsi_f
  END INTERFACE

  error = xfsetgridcoordsi_f (a_Id, a_NumVals, a_iValues)

return

END SUBROUTINE

!******************************

SUBROUTINE XF_SET_GRID_COORDS_J (a_Id, a_NumVals, a_jValues, error)
INTEGER             , INTENT(IN)  :: a_Id
INTEGER, INTENT(IN)         :: a_NumVals
REAL*8, DIMENSION(*), INTENT(IN)    :: a_jValues
INTEGER, INTENT(OUT)        :: error

!            INTEGER, EXTERNAL :: _xfSetGridCoordsJ_f
!  MS FORTRAN needs explicit interface for C functions called here.
!
  INTERFACE
    INTEGER FUNCTION xfsetgridcoordsj_f(a_Id, a_NumVals, a_jValues)
	  !MS$ATTRIBUTES C,reference::xfsetgridcoordsj_f
    INTEGER,              INTENT(IN)    :: a_Id
    INTEGER,              INTENT(IN)    :: a_NumVals
    REAL*8, DIMENSION(*), INTENT(IN)    :: a_jValues


    END FUNCTION xfsetgridcoordsj_f
  END INTERFACE

  error = xfsetgridcoordsj_f (a_Id, a_NumVals, a_jValues)

return

END SUBROUTINE

!******************************

SUBROUTINE XF_SET_GRID_COORDS_K (a_Id, a_NumVals, a_kValues, error)
INTEGER             , INTENT(IN)  :: a_Id
INTEGER, INTENT(IN)         :: a_NumVals
REAL*8, DIMENSION(*), INTENT(IN)    :: a_kValues
INTEGER, INTENT(OUT)        :: error

!            INTEGER, EXTERNAL :: _xfSetGridCoordsK_f
!  MS FORTRAN needs explicit interface for C functions called here.
!
  INTERFACE
    INTEGER FUNCTION xfsetgridcoordsk_f(a_Id, a_NumVals, a_kValues)
	  !MS$ATTRIBUTES C,reference::xfsetgridcoordsk_f
    INTEGER,              INTENT(IN)    :: a_Id
    INTEGER,              INTENT(IN)    :: a_NumVals
    REAL*8, DIMENSION(*), INTENT(IN)    :: a_kValues


    END FUNCTION xfsetgridcoordsk_f
  END INTERFACE

  error = xfsetgridcoordsk_f (a_Id, a_NumVals, a_kValues)

return

END SUBROUTINE

!******************************

SUBROUTINE XF_WRITE_EXTRUDE_LAYER_DATA (a_Id, a_NumLayers, a_NumVals, a_Values, error)
INTEGER             , INTENT(IN) :: a_Id
INTEGER, INTENT(IN)        :: a_NumLayers
INTEGER, INTENT(IN)        :: a_NumVals
REAL*8, DIMENSION(*), INTENT(IN)   :: a_Values
INTEGER, INTENT(OUT)       :: error

!            INTEGER, EXTERNAL :: _xfWriteExtrudeLayerData_f
!  MS FORTRAN needs explicit interface for C functions called here.
!
  INTERFACE
    INTEGER FUNCTION xfwriteextrudelayerdata_f(a_Id, a_NumLayers, a_NumVals, a_Values)
	  !MS$ATTRIBUTES C,reference::xfwriteextrudelayerdata_f
    INTEGER,              INTENT(IN)    :: a_Id
    INTEGER,              INTENT(IN)    :: a_NumLayers
    INTEGER,              INTENT(IN)    :: a_NumVals
    REAL*8, DIMENSION(*), INTENT(IN)    :: a_Values


    END FUNCTION xfwriteextrudelayerdata_f
  END INTERFACE

  error = xfwriteextrudelayerdata_f (a_Id, a_NumLayers, a_NumVals, a_Values)

return

END SUBROUTINE

!******************************

SUBROUTINE XF_GET_GRID_TYPE (a_Id, a_GridType, error)
INTEGER             , INTENT(IN) :: a_Id
INTEGER, INTENT(OUT)       :: a_GridType
INTEGER, INTENT(OUT)       :: error

!            INTEGER, EXTERNAL :: _xfGetGridType_f
!  MS FORTRAN needs explicit interface for C functions called here.
!
  INTERFACE
    INTEGER FUNCTION xfgetgridtype_f(a_Id, a_GridType)
	  !MS$ATTRIBUTES C,reference::xfgetgridtype_f
    INTEGER,              INTENT(IN)    :: a_Id
    INTEGER,              INTENT(OUT)   :: a_GridType


    END FUNCTION xfgetgridtype_f
  END INTERFACE

  error = xfgetgridtype_f (a_Id, a_GridType)

return

END SUBROUTINE

!******************************

SUBROUTINE XF_GET_EXTRUSION_TYPE (a_Id, a_ExtrudeType, error)
INTEGER             , INTENT(IN) :: a_Id
INTEGER, INTENT(OUT)       :: a_ExtrudeType
INTEGER, INTENT(OUT)       :: error

!            INTEGER, EXTERNAL :: _xfGetExtrusionType_f
!  MS FORTRAN needs explicit interface for C functions called here.
!
  INTERFACE
    INTEGER FUNCTION xfgetextrusiontype_f(a_Id, a_ExtrudeType)
	  !MS$ATTRIBUTES C,reference::xfgetextrusiontype_f
    INTEGER,              INTENT(IN)    :: a_Id
    INTEGER,              INTENT(OUT)   :: a_ExtrudeType


    END FUNCTION xfgetextrusiontype_f
  END INTERFACE

  error = xfgetextrusiontype_f (a_Id, a_ExtrudeType)

return

END SUBROUTINE

!******************************

SUBROUTINE XF_GET_NUMBER_OF_DIMENSIONS(a_Id, a_NumDimensions, error)
INTEGER             , INTENT(IN) :: a_Id
INTEGER, INTENT(OUT)       :: a_NumDimensions
INTEGER, INTENT(OUT)       :: error

!            INTEGER, EXTERNAL :: _xfGetNumberOfDimensions_f
!  MS FORTRAN needs explicit interface for C functions called here.
!
  INTERFACE
    INTEGER FUNCTION xfgetnumberofdimensions_f(a_Id, a_NumDimensions)
	  !MS$ATTRIBUTES C,reference::xfgetnumberofdimensions_f
    INTEGER,              INTENT(IN)    :: a_Id
    INTEGER,              INTENT(OUT)   :: a_NumDimensions


    END FUNCTION xfgetnumberofdimensions_f
  END INTERFACE

  error = xfgetnumberofdimensions_f (a_Id, a_NumDimensions)

return

END SUBROUTINE

!******************************

SUBROUTINE XF_ORIGIN_DEFINED (a_Id, a_bDefined, error)
INTEGER             , INTENT(IN) :: a_Id
LOGICAL*2, INTENT(OUT)       :: a_bDefined
INTEGER, INTENT(OUT)       :: error


!            INTEGER, EXTERNAL :: _xfOriginDefined_f
!  MS FORTRAN needs explicit interface for C functions called here.
!
  INTERFACE
    INTEGER FUNCTION xforigindefined_f(a_Id, a_bDefined)
	  !MS$ATTRIBUTES C,reference::xforigindefined_f
    INTEGER,              INTENT(IN)    :: a_Id
    LOGICAL*2,              INTENT(OUT)   :: a_bDefined


    END FUNCTION xforigindefined_f
  END INTERFACE

  error = xforigindefined_f (a_Id, a_bDefined)


return

END SUBROUTINE

!******************************
SUBROUTINE XF_GET_ORIGIN (a_Id, a_x, a_y, a_z, error)
INTEGER, INTENT(IN)   :: a_Id
REAL*8,  INTENT(OUT)  :: a_x
REAL*8,  INTENT(OUT)  :: a_y
REAL*8,  INTENT(OUT)  :: a_z
INTEGER, INTENT(OUT)  :: error

!            INTEGER, EXTERNAL :: _xfGetOrigin_f
!  MS FORTRAN needs explicit interface for C functions called here.
!
  INTERFACE
    INTEGER FUNCTION xfgetorigin_f(a_Id, a_x, a_y, a_z)
	  !MS$ATTRIBUTES C,reference::xfgetorigin_f
    INTEGER, INTENT(IN)    :: a_Id
    REAL*8,  INTENT(OUT)   :: a_x, a_y, a_z


    END FUNCTION xfgetorigin_f
  END INTERFACE

  error = xfgetorigin_f (a_Id, a_x, a_y, a_z)


return

END SUBROUTINE

!******************************

SUBROUTINE XF_GET_ORIENTATION (a_Id, a_Orientation, error)
INTEGER             , INTENT(IN) :: a_Id
INTEGER, INTENT(OUT)       :: a_Orientation
INTEGER, INTENT(OUT)       :: error


!            INTEGER, EXTERNAL :: _xfGetOrientation_f
!  MS FORTRAN needs explicit interface for C functions called here.
!
  INTERFACE
    INTEGER FUNCTION xfgetorientation_f(a_Id, a_Orientation)
	  !MS$ATTRIBUTES C,reference::xfgetorientation_f
    INTEGER,              INTENT(IN)    :: a_Id
    INTEGER,              INTENT(OUT)   :: a_Orientation


    END FUNCTION xfgetorientation_f
  END INTERFACE

  error = xfgetorientation_f (a_Id, a_Orientation)

return

END SUBROUTINE

!******************************

SUBROUTINE XF_BEARING_DEFINED (a_Id, a_bDefined, error)
INTEGER             , INTENT(IN)  :: a_Id
LOGICAL*2, INTENT(OUT)        :: a_bDefined
INTEGER, INTENT(OUT)        :: error


!            INTEGER, EXTERNAL :: _xfBearingDefined_f
!  MS FORTRAN needs explicit interface for C functions called here.
!
  INTERFACE
    INTEGER FUNCTION xfbearingdefined_f(a_Id, a_bDefined)
	  !MS$ATTRIBUTES C,reference::xfbearingdefined_f
    INTEGER,              INTENT(IN)    :: a_Id
    LOGICAL*2,              INTENT(OUT)   :: a_bDefined


    END FUNCTION xfbearingdefined_f
  END INTERFACE

  error = xfbearingdefined_f (a_Id, a_bDefined)

return

END SUBROUTINE

!******************************

SUBROUTINE XF_GET_BEARING (a_Id, a_bearing, error)
INTEGER, INTENT(IN)    :: a_Id
REAL*8,  INTENT(OUT)   :: a_bearing
INTEGER, INTENT(OUT)   :: error


!            INTEGER, EXTERNAL :: _xfGetBearing_f
!  MS FORTRAN needs explicit interface for C functions called here.
!
  INTERFACE
    INTEGER FUNCTION xfgetbearing_f(a_Id, a_bearing)
	  !MS$ATTRIBUTES C,reference::xfgetbearing_f
    INTEGER,              INTENT(IN)    :: a_Id
    REAL*8,               INTENT(OUT)   :: a_bearing


    END FUNCTION xfgetbearing_f
  END INTERFACE

  error = xfgetbearing_f (a_Id, a_bearing)

return

END SUBROUTINE

!******************************

SUBROUTINE XF_DIP_DEFINED (a_Id, a_bDefined, error)
INTEGER             , INTENT(IN) :: a_Id
LOGICAL*2, INTENT(OUT)       :: a_bDefined
INTEGER, INTENT(OUT)       :: error

!            INTEGER, EXTERNAL :: _xfDipDefined_f
!  MS FORTRAN needs explicit interface for C functions called here.
!
  INTERFACE
    INTEGER FUNCTION xfdipdefined_f(a_Id, a_bDefined)
	  !MS$ATTRIBUTES C,reference::xfdipdefined_f
    INTEGER,              INTENT(IN)    :: a_Id
    LOGICAL*2,              INTENT(OUT)   :: a_bDefined


    END FUNCTION xfdipdefined_f
  END INTERFACE

  error = xfdipdefined_f (a_Id, a_bDefined)

return

END SUBROUTINE

!******************************

SUBROUTINE XF_GET_DIP (a_Id, a_dip, error)
INTEGER, INTENT(IN)   :: a_Id
REAL*8,  INTENT(OUT)  :: a_dip
INTEGER, INTENT(OUT)  :: error


!            INTEGER, EXTERNAL :: _xfGetDip_f
!  MS FORTRAN needs explicit interface for C functions called here.
!
  INTERFACE
    INTEGER FUNCTION xfgetdip_f(a_Id, a_dip)
	  !MS$ATTRIBUTES C,reference::xfgetdip_f
    INTEGER,  INTENT(IN)    :: a_Id
    REAL*8,   INTENT(OUT)   :: a_dip


    END FUNCTION xfgetdip_f
  END INTERFACE

  error = xfgetdip_f (a_Id, a_dip)

return

END SUBROUTINE

!******************************

SUBROUTINE XF_ROLL_DEFINED (a_Id, a_bDefined, error)
INTEGER             , INTENT(IN) :: a_Id
LOGICAL*2, INTENT(OUT)       :: a_bDefined
INTEGER, INTENT(OUT)       :: error


!            INTEGER, EXTERNAL :: _xfRollDefined_f
!  MS FORTRAN needs explicit interface for C functions called here.
!
  INTERFACE
    INTEGER FUNCTION xfrolldefined_f(a_Id, a_bDefined)
	  !MS$ATTRIBUTES C,reference::xfrolldefined_f
    INTEGER,              INTENT(IN)    :: a_Id
    LOGICAL*2,              INTENT(OUT)   :: a_bDefined


    END FUNCTION xfrolldefined_f
  END INTERFACE

  error = xfrolldefined_f (a_Id, a_bDefined)

return

END SUBROUTINE

!******************************

SUBROUTINE XF_GET_ROLL (a_Id, a_Roll, error)
INTEGER, INTENT(IN)   :: a_Id
REAL*8,  INTENT(OUT)  :: a_Roll
INTEGER, INTENT(OUT)  :: error

!            INTEGER, EXTERNAL :: _xfGetRoll_f
!  MS FORTRAN needs explicit interface for C functions called here.
!
  INTERFACE
    INTEGER FUNCTION xfgetroll_f(a_Id, a_Roll)
	  !MS$ATTRIBUTES C,reference::xfgetroll_f
    INTEGER, INTENT(IN)    :: a_Id
    REAL*8,  INTENT(OUT)   :: a_Roll


    END FUNCTION xfgetroll_f
  END INTERFACE

  error = xfgetroll_f (a_Id, a_Roll)

return

END SUBROUTINE

!********************************

SUBROUTINE XF_COMPUTATIONAL_ORIGIN_DEFINED (a_Id, a_bDefined, error)
INTEGER             , INTENT(IN) :: a_Id
LOGICAL*2, INTENT(OUT)       :: a_bDefined
INTEGER                       error

!            INTEGER, EXTERNAL :: _xfComputationalOriginDefined_f
!  MS FORTRAN needs explicit interface for C functions called here.
!
  INTERFACE
    INTEGER FUNCTION xfcomputationalorigindefined_f(a_Id, a_bDefined)
	  !MS$ATTRIBUTES C,reference::xfcomputationalorigindefined_f
    INTEGER,              INTENT(IN)    :: a_Id
    LOGICAL*2,              INTENT(OUT)   :: a_bDefined


    END FUNCTION xfcomputationalorigindefined_f
  END INTERFACE

  error = xfcomputationalorigindefined_f (a_Id, a_bDefined)
return

END SUBROUTINE 

!********************************

SUBROUTINE XF_GET_COMPUTATIONAL_ORIGIN (a_Id, a_Origin, error)
INTEGER             , INTENT(IN) :: a_Id
INTEGER, INTENT(OUT)       :: a_Origin
INTEGER                       error

!            INTEGER, EXTERNAL :: _xfGetComputationalOrigin_f
!  MS FORTRAN needs explicit interface for C functions called here.
!
  INTERFACE
    INTEGER FUNCTION xfgetcomputationalorigin_f(a_Id, a_Origin)
	  !MS$ATTRIBUTES C,reference::xfgetcomputationalorigin_f
    INTEGER,              INTENT(IN)    :: a_Id
    INTEGER,              INTENT(OUT)   :: a_Origin


    END FUNCTION xfgetcomputationalorigin_f
  END INTERFACE

  error = xfgetcomputationalorigin_f (a_Id, a_Origin)

return

END SUBROUTINE

!**********************************

SUBROUTINE XF_GET_U_DIRECTION_DEFINED (a_Id, a_bDefined, error)
INTEGER             , INTENT(IN) :: a_Id
LOGICAL*2, INTENT(OUT)       :: a_bDefined
INTEGER                       error

!            INTEGER, EXTERNAL :: _xfGetUDirectionDefined_f
!  MS FORTRAN needs explicit interface for C functions called here.
!
  INTERFACE
    INTEGER FUNCTION xfgetudirectiondefined_f(a_Id, a_bDefined)
	  !MS$ATTRIBUTES C,reference::xfgetudirectiondefined_f
    INTEGER,              INTENT(IN)    :: a_Id
    LOGICAL*2,              INTENT(OUT)   :: a_bDefined


    END FUNCTION xfgetudirectiondefined_f
  END INTERFACE

  error = xfgetudirectiondefined_f (a_Id, a_bDefined)
return

END SUBROUTINE

!**********************************

SUBROUTINE XF_GET_U_DIRECTION (a_Id, a_Direction, error)
INTEGER             , INTENT(IN) :: a_Id
INTEGER, INTENT(OUT)       :: a_Direction
INTEGER                       error

!            INTEGER, EXTERNAL :: _xfGetUDirection_f
!  MS FORTRAN needs explicit interface for C functions called here.
!
  INTERFACE
    INTEGER FUNCTION xfgetudirection_f(a_Id, a_Direction)
	  !MS$ATTRIBUTES C,reference::xfgetudirection_f
    INTEGER,              INTENT(IN)    :: a_Id
    INTEGER,              INTENT(OUT)   :: a_Direction


    END FUNCTION xfgetudirection_f
  END INTERFACE

  error = xfgetudirection_f (a_Id, a_Direction)

return

END SUBROUTINE

!******************************

SUBROUTINE XF_GET_NUMBER_CELLS_IN_I (a_Id, a_NumI, error)
INTEGER             , INTENT(IN) :: a_Id
INTEGER, INTENT(OUT)       :: a_NumI
INTEGER, INTENT(OUT)       :: error


!            INTEGER, EXTERNAL :: _xfGetNumberCellsInI_f
!  MS FORTRAN needs explicit interface for C functions called here.
!
  INTERFACE
    INTEGER FUNCTION xfgetnumbercellsini_f(a_Id, a_NumI)
	  !MS$ATTRIBUTES C,reference::xfgetnumbercellsini_f
    INTEGER,              INTENT(IN)    :: a_Id
    INTEGER,              INTENT(OUT)   :: a_NumI


    END FUNCTION xfgetnumbercellsini_f
  END INTERFACE

  error = xfgetnumbercellsini_f (a_Id, a_NumI)

return

END SUBROUTINE

!******************************

SUBROUTINE XF_GET_NUMBER_CELLS_IN_J (a_Id, a_NumJ, error)
INTEGER             , INTENT(IN) :: a_Id
INTEGER, INTENT(OUT)       :: a_NumJ
INTEGER, INTENT(OUT)       :: error


!            INTEGER, EXTERNAL :: _xfGetNumberCellsInJ_f
!  MS FORTRAN needs explicit interface for C functions called here.
!
  INTERFACE
    INTEGER FUNCTION xfgetnumbercellsinj_f(a_Id, a_NumJ)
	  !MS$ATTRIBUTES C,reference::xfgetnumbercellsinj_f
    INTEGER,              INTENT(IN)    :: a_Id
    INTEGER,              INTENT(OUT)   :: a_NumJ


    END FUNCTION xfgetnumbercellsinj_f
  END INTERFACE

  error = xfgetnumbercellsinj_f (a_Id, a_NumJ)

return

END SUBROUTINE

!******************************

SUBROUTINE XF_GET_NUMBER_CELLS_IN_K (a_Id, a_NumK, error)
INTEGER             , INTENT(IN) :: a_Id
INTEGER, INTENT(OUT)       :: a_NumK
INTEGER, INTENT(OUT)       :: error


!            INTEGER, EXTERNAL :: _xfGetNumberCellsInK_f
!  MS FORTRAN needs explicit interface for C functions called here.
!
  INTERFACE
    INTEGER FUNCTION xfgetnumbercellsink_f(a_Id, a_NumK)
	  !MS$ATTRIBUTES C,reference::xfgetnumbercellsink_f
    INTEGER,              INTENT(IN)    :: a_Id
    INTEGER,              INTENT(OUT)   :: a_NumK


    END FUNCTION xfgetnumbercellsink_f
  END INTERFACE

  error = xfgetnumbercellsink_f (a_Id, a_NumK)

return

END SUBROUTINE

!******************************

SUBROUTINE XF_GET_GRID_COORDS_I (a_Id, a_NumVals, a_iValues, error)
INTEGER             , INTENT(IN) :: a_Id
INTEGER, INTENT(IN)        :: a_NumVals
REAL*8, DIMENSION(*), INTENT(OUT)  :: a_iValues
INTEGER, INTENT(OUT)       :: error


!            INTEGER, EXTERNAL :: _xfGetGridCoordsI_f
!  MS FORTRAN needs explicit interface for C functions called here.
!
  INTERFACE
    INTEGER FUNCTION xfgetgridcoordsi_f(a_Id, a_NumVals, a_iValues)
	  !MS$ATTRIBUTES C,reference::xfgetgridcoordsi_f
    INTEGER,              INTENT(IN)    :: a_Id
    INTEGER,              INTENT(IN)    :: a_NumVals
    REAL*8, DIMENSION(*), INTENT(OUT)   :: a_iValues


    END FUNCTION xfgetgridcoordsi_f 
  END INTERFACE

  error = xfgetgridcoordsi_f (a_Id, a_NumVals, a_iValues)

return

END SUBROUTINE

!******************************

SUBROUTINE XF_GET_GRID_COORDS_J (a_Id, a_NumVals, a_jValues, error)
INTEGER             , INTENT(IN) :: a_Id
INTEGER, INTENT(IN)        :: a_NumVals
REAL*8, DIMENSION(*), INTENT(OUT)  :: a_jValues
INTEGER, INTENT(OUT)       :: error


!            INTEGER, EXTERNAL :: _xfGetGridCoordsJ_f
!  MS FORTRAN needs explicit interface for C functions called here.
!
  INTERFACE
    INTEGER FUNCTION xfgetgridcoordsj_f(a_Id, a_NumVals, a_jValues)
	  !MS$ATTRIBUTES C,reference::xfgetgridcoordsj_f
    INTEGER,              INTENT(IN)    :: a_Id
    INTEGER,              INTENT(IN)    :: a_NumVals
    REAL*8, DIMENSION(*), INTENT(OUT)   :: a_jValues


    END FUNCTION xfgetgridcoordsj_f
  END INTERFACE

  error = xfgetgridcoordsj_f (a_Id, a_NumVals, a_jValues)

return

END SUBROUTINE

!******************************

SUBROUTINE XF_GET_GRID_COORDS_K (a_Id, a_NumVals, a_kValues, error)
INTEGER             , INTENT(IN) :: a_Id
INTEGER, INTENT(IN)        :: a_NumVals
REAL*8, DIMENSION(*), INTENT(OUT)  :: a_kValues
INTEGER, INTENT(OUT)       :: error


!            INTEGER, EXTERNAL :: _xfGetGridCoordsK_f
!  MS FORTRAN needs explicit interface for C functions called here.
!
  INTERFACE
    INTEGER FUNCTION xfgetgridcoordsk_f(a_Id, a_NumVals, a_kValues)
	  !MS$ATTRIBUTES C,reference::xfgetgridcoordsk_f
    INTEGER,              INTENT(IN)    :: a_Id
    INTEGER,              INTENT(IN)    :: a_NumVals
    REAL*8, DIMENSION(*), INTENT(OUT)   :: a_kValues


    END FUNCTION xfgetgridcoordsk_f
  END INTERFACE

  error = xfgetgridcoordsk_f(a_Id, a_NumVals, a_kValues)

return

END SUBROUTINE

!******************************

SUBROUTINE XF_GET_EXTRUDE_NUM_LAYERS (a_Id, a_NumLayers, error)
INTEGER             , INTENT(IN) :: a_Id
INTEGER, INTENT(OUT)       :: a_NumLayers
INTEGER, INTENT(OUT)       :: error

!            INTEGER, EXTERNAL :: _xfGetExtrudeNumLayers_f
!  MS FORTRAN needs explicit interface for C functions called here.
!
  INTERFACE
    INTEGER FUNCTION xfgetextrudenumlayers_f(a_Id, a_NumLayers)
	  !MS$ATTRIBUTES C,reference::xfgetextrudenumlayers_f
    INTEGER,              INTENT(IN)    :: a_Id
    INTEGER,              INTENT(OUT)   :: a_NumLayers


    END FUNCTION xfgetextrudenumlayers_f
  END INTERFACE

  error = xfgetextrudenumlayers_f (a_Id, a_NumLayers)

return

END SUBROUTINE

!******************************

SUBROUTINE XF_GET_EXTRUDE_VALUES (a_Id, a_NumVals, a_Values, error)
INTEGER             , INTENT(IN) :: a_Id
INTEGER, INTENT(IN)        :: a_NumVals
REAL*8, DIMENSION(*), INTENT(OUT)  :: a_Values
INTEGER, INTENT(OUT)       :: error


!            INTEGER, EXTERNAL :: _xfGetExtrudeValues_f
!  MS FORTRAN needs explicit interface for C functions called here.
!
  INTERFACE
    INTEGER FUNCTION xfgetextrudevalues_f(a_Id, a_NumVals, a_Values)
	  !MS$ATTRIBUTES C,reference::xfgetextrudevalues_f
    INTEGER,              INTENT(IN)    :: a_Id
    INTEGER,              INTENT(IN)    :: a_NumVals
    REAL*8, DIMENSION(*), INTENT(OUT)   :: a_Values


    END FUNCTION xfgetextrudevalues_f
  END INTERFACE

  error = xfgetextrudevalues_f (a_Id, a_NumVals, a_Values)

return

END SUBROUTINE

!******************************

SUBROUTINE XF_CREATE_GRID_PROPERTY_GROUP (a_Id, a_PropId, error)
INTEGER             , INTENT(IN)  :: a_Id
INTEGER             , INTENT(OUT) :: a_PropId
INTEGER, INTENT(OUT)        :: error


!            INTEGER, EXTERNAL :: _xfCreateGridPropertyGroup_f
!  MS FORTRAN needs explicit interface for C functions called here.
!
  INTERFACE
    INTEGER FUNCTION xfcreategridpropertygroup_f(a_Id, a_PropId)
	  !MS$ATTRIBUTES C,reference::xfcreategridpropertygroup_f
    INTEGER,              INTENT(IN)    :: a_Id
    INTEGER,              INTENT(OUT)   :: a_PropId


    END FUNCTION xfcreategridpropertygroup_f
  END INTERFACE

  error = xfcreategridpropertygroup_f (a_Id, a_PropId)

return

END SUBROUTINE

!******************************

SUBROUTINE XF_GET_GRID_PROPERTY_GROUP (a_Id, a_PropId, error)
INTEGER             , INTENT(IN)  :: a_Id
INTEGER             , INTENT(OUT) :: a_PropId
INTEGER, INTENT(OUT)        :: error


!            INTEGER, EXTERNAL :: _xfGetGridPropertyGroup_f
!  MS FORTRAN needs explicit interface for C functions called here.
!
  INTERFACE
    INTEGER FUNCTION xfgetgridpropertygroup_f(a_Id, a_PropId)
	  !MS$ATTRIBUTES C,reference::xfgetgridpropertygroup_f
    INTEGER,              INTENT(IN)    :: a_Id
    INTEGER,              INTENT(OUT)   :: a_PropId


    END FUNCTION xfgetgridpropertygroup_f
  END INTERFACE

  error = xfgetgridpropertygroup_f (a_Id, a_PropId)

return

END SUBROUTINE

!******************************

SUBROUTINE XF_CREATE_GRID_CELL_PROP_GRP (a_Id, a_PropId, error)
INTEGER             , INTENT(IN)  :: a_Id
INTEGER             , INTENT(OUT) :: a_PropId
INTEGER, INTENT(OUT)        :: error


!            INTEGER, EXTERNAL :: _xfCreateGridCellPropertyGroup_f
!  MS FORTRAN needs explicit interface for C functions called here.
!
  INTERFACE
    INTEGER FUNCTION xfcreategridcellpropertygroup_f(a_Id, a_PropId)
	  !MS$ATTRIBUTES C,reference::xfcreategridcellpropertygroup_f
    INTEGER,              INTENT(IN)    :: a_Id
    INTEGER,              INTENT(OUT)   :: a_PropId


    END FUNCTION xfcreategridcellpropertygroup_f
  END INTERFACE

  error = xfcreategridcellpropertygroup_f (a_Id, a_PropId)

return

END SUBROUTINE

!******************************

SUBROUTINE XF_GET_GRID_CELL_PROPERTY_GROUP (a_Id, a_PropId, error)
INTEGER             , INTENT(IN)  :: a_Id
INTEGER             , INTENT(OUT) :: a_PropId
INTEGER, INTENT(OUT)        :: error


!            INTEGER, EXTERNAL :: _xfGetGridCellPropertyGroup_f
!  MS FORTRAN needs explicit interface for C functions called here.
!
  INTERFACE
    INTEGER FUNCTION xfgetgridcellpropertygroup_f(a_Id, a_PropId)
	  !MS$ATTRIBUTES C,reference::xfgetgridcellpropertygroup_f
    INTEGER,              INTENT(IN)    :: a_Id
    INTEGER,              INTENT(OUT)   :: a_PropId

    END FUNCTION xfgetgridcellpropertygroup_f
  END INTERFACE

  error = xfgetgridcellpropertygroup_f (a_Id, a_PropId)

return

END SUBROUTINE

!******************************

SUBROUTINE XF_CREATE_GRID_NODE_PROP_GRP (a_Id, a_PropId, error)
INTEGER             , INTENT(IN)  :: a_Id
INTEGER             , INTENT(OUT) :: a_PropId
INTEGER, INTENT(OUT)        :: error


!            INTEGER, EXTERNAL :: _xfCreateGridNodePropertyGroup_f
!  MS FORTRAN needs explicit interface for C functions called here.
!
  INTERFACE
    INTEGER FUNCTION xfcreategridnodepropertygroup_f(a_Id, a_PropId)
	  !MS$ATTRIBUTES C,reference::xfcreategridnodepropertygroup_f
    INTEGER,              INTENT(IN)    :: a_Id
    INTEGER,              INTENT(OUT)   :: a_PropId


    END FUNCTION xfcreategridnodepropertygroup_f
  END INTERFACE

  error = xfcreategridnodepropertygroup_f (a_Id, a_PropId)

return

END SUBROUTINE

!******************************

SUBROUTINE XF_GET_GRID_NODE_PROPERTY_GROUP (a_Id, a_PropId, error)
INTEGER             , INTENT(IN)  :: a_Id
INTEGER             , INTENT(OUT) :: a_PropId
INTEGER, INTENT(OUT)        :: error


!            INTEGER, EXTERNAL :: _xfGetGridNodePropertyGroup_f
!  MS FORTRAN needs explicit interface for C functions called here.
!
  INTERFACE
    INTEGER FUNCTION xfgetgridnodepropertygroup_f(a_Id, a_PropId)
	  !MS$ATTRIBUTES C,reference::xfgetgridnodepropertygroup_f
    INTEGER,              INTENT(IN)    :: a_Id
    INTEGER,              INTENT(OUT)   :: a_PropId


    END FUNCTION xfgetgridnodepropertygroup_f
  END INTERFACE

  error = xfgetgridnodepropertygroup_f (a_Id, a_PropId)

return

END SUBROUTINE

! --------------------------------------------------------------------------- */
! FUNCTION  XF_CREATE_GEOMETRIC_PATH_GROUP 
! PURPOSE    
! NOTES     
! --------------------------------------------------------------------------- */
SUBROUTINE XF_CREATE_GEOMETRIC_PATH_GROUP(a_ParentId, a_Path, a_Guid,        &
                               a_Compression, a_PathGroup, a_NullVal, error)
  INTEGER,              INTENT(IN)   :: a_ParentId
  CHARACTER(LEN=*),     INTENT(IN)   :: a_Path, a_Guid
  INTEGER,              INTENT(IN)   :: a_Compression
  INTEGER,              INTENT(OUT)  :: a_PathGroup
  REAL*8,               INTENT(IN)   :: a_NullVal
  INTEGER,              INTENT(OUT)  :: error
  INTEGER                          pathlen, guidlen

 !            INTEGER, EXTERNAL :: _xfCreateGeometricPathGroup_f
!  MS FORTRAN needs explicit interface for C functions called here.
!
  INTERFACE
    INTEGER FUNCTION xfcreategeometricpathgroup_f(a_ParentId, a_Path, pathlen, &
                                                  a_Guid, guidlen, a_Compression, &
                                                  a_PathGroup, a_NullVal)
	  !MS$ATTRIBUTES C,reference::xfcreategeometricpathgroup_f
    !DEC$ATTRIBUTES reference :: a_Path
    !DEC$ATTRIBUTES reference :: a_Guid
    INTEGER,              INTENT(IN)    :: a_ParentId, a_Compression
    INTEGER,              INTENT(OUT)   :: a_PathGroup
    CHARACTER(LEN=*),     INTENT(IN)    :: a_Path, a_Guid
    REAL*8,               INTENT(IN)    :: a_NullVal
    INTEGER,              INTENT(IN)    :: pathlen, guidlen


    END FUNCTION xfcreategeometricpathgroup_f
  END INTERFACE

  pathlen = LEN_TRIM(a_Path)
  guidlen = LEN_TRIM(a_Guid)
  error = xfcreategeometricpathgroup_f (a_ParentId, a_Path, pathlen, &
                                        a_Guid, guidlen, a_Compression, &
                                        a_PathGroup, a_NullVal)


  return
END SUBROUTINE ! XF_CREATE_GEOMETRIC_PATH_GROUP

! --------------------------------------------------------------------------- */
! FUNCTION  XF_WRITE_PARTICLE_TIMESTEP */
! PURPOSE    */
! NOTES      */
! --------------------------------------------------------------------------- */
SUBROUTINE XF_WRITE_PARTICLE_TIMESTEP(a_Id, a_nDim, a_Time, a_nPaths, a_Locs, &
                                      error)
  INTEGER,              INTENT(IN)   :: a_Id
  INTEGER,              INTENT(IN)   :: a_nDim
  REAL*8,               INTENT(IN)   :: a_Time
  INTEGER,              INTENT(IN)   :: a_nPaths
  REAL*8, DIMENSION(*), INTENT(IN)   :: a_Locs
  INTEGER,              INTENT(OUT)  :: error


!            INTEGER, EXTERNAL :: _xfWriteParticleTimestep_f
!  MS FORTRAN needs explicit interface for C functions called here.
!
  INTERFACE
    INTEGER FUNCTION xfwriteparticletimestep_f(a_Id, a_nDim, a_Time, a_nPaths, a_Locs)
	  !MS$ATTRIBUTES C,reference::xfwriteparticletimestep_f
    INTEGER,              INTENT(IN)    :: a_Id
    INTEGER,              INTENT(IN)    :: a_nDim, a_nPaths
	REAL*8,               INTENT(IN)    :: a_Time
    REAL*8, DIMENSION(*), INTENT(IN)    :: a_Locs


    END FUNCTION xfwriteparticletimestep_f
  END INTERFACE

  error = xfwriteparticletimestep_f (a_Id, a_nDim, a_Time, a_nPaths, a_Locs)

  return 
END SUBROUTINE ! XF_WRITE_PARTICLE_TIMESTEP

! --------------------------------------------------------------------------- */
! FUNCTION  XF_GET_PATH_NULL_VAL
! PURPOSE    
! NOTES      
! --------------------------------------------------------------------------- */
SUBROUTINE XF_GET_PATH_NULL_VAL(GroupId, NullVal, error)
  INTEGER, INTENT(IN)   ::  GroupId
  REAL*8,  INTENT(OUT)  ::  NullVal
  INTEGER, INTENT(OUT)  ::  error

  !            INTEGER, EXTERNAL :: _xfGetPathNullVal_f
!  MS FORTRAN needs explicit interface for C functions called here.
!
  INTERFACE
    INTEGER FUNCTION xfgetpathnullval_f(GroupId, NullVal)
	  !MS$ATTRIBUTES C,reference::xfgetpathnullval_f
    INTEGER,              INTENT(IN)    :: GroupId
    REAL*8,               INTENT(OUT)   :: NullVal


    END FUNCTION xfgetpathnullval_f
  END INTERFACE

  error = xfgetpathnullval_f (GroupId, NullVal)


END SUBROUTINE ! XF_GET_PATH_NULL_VAL 

! --------------------------------------------------------------------------- */
! FUNCTION  XF_GET_NUMBER_OF_PATHS
! PURPOSE    
! NOTES      
! --------------------------------------------------------------------------- */
SUBROUTINE XF_GET_NUMBER_OF_PATHS(GroupId, NumPaths, error)
  INTEGER             , INTENT(IN)   :: GroupId
  INTEGER, INTENT(OUT)         :: NumPaths
  INTEGER, INTENT(OUT)         :: error

  
  !            INTEGER, EXTERNAL :: _xfGetNumberOfPaths_f
!  MS FORTRAN needs explicit interface for C functions called here.
!
  INTERFACE
    INTEGER FUNCTION xfgetnumberofpaths_f(GroupId, NumPaths)
	  !MS$ATTRIBUTES C,reference::xfgetnumberofpaths_f
    INTEGER,              INTENT(IN)    :: GroupId
    INTEGER,              INTENT(OUT)   :: NumPaths


    END FUNCTION xfgetnumberofpaths_f
  END INTERFACE

  error = xfgetnumberofpaths_f (GroupId, NumPaths)

  return
END SUBROUTINE ! XF_GET_NUMBER_OF_PATHS */

! --------------------------------------------------------------------------- */
! FUNCTION  XF_GET_NUMBER_OF_TIMES */
! PURPOSE    */
! NOTES      */
! --------------------------------------------------------------------------- */
SUBROUTINE XF_GET_NUMBER_OF_TIMES(GroupId, NumTimes, error)
  INTEGER             , INTENT(IN)   ::  GroupId
  INTEGER, INTENT(OUT)         ::  NumTimes
  INTEGER, INTENT(OUT)         ::  error
  

  !            INTEGER, EXTERNAL :: _xfGetNumberOfTimes_f
!  MS FORTRAN needs explicit interface for C functions called here.
!
  INTERFACE
    INTEGER FUNCTION xfgetnumberoftimes_f(GroupId, NumTimes)
	  !MS$ATTRIBUTES C,reference::xfgetnumberoftimes_f
    INTEGER,              INTENT(IN)    :: GroupId
    INTEGER,              INTENT(OUT)   :: NumTimes


    END FUNCTION xfgetnumberoftimes_f
  END INTERFACE

  error = xfgetnumberoftimes_f (GroupId, NumTimes)


  return
END SUBROUTINE ! XF_GET_NUMBER_OF_TIMES 

! --------------------------------------------------------------------------- */
! FUNCTION  XF_GET_PATH_DIMENSIONALITY
! PURPOSE   
! NOTES     
! --------------------------------------------------------------------------- */
SUBROUTINE XF_GET_PATH_DIMENSIONALITY(GroupId, NumDims, error)
  INTEGER             , INTENT(IN)  ::  GroupId
  INTEGER, INTENT(OUT)        ::  NumDims, error
  

  !            INTEGER, EXTERNAL :: _xfGetPathDimensionality_f
!  MS FORTRAN needs explicit interface for C functions called here.
!
  INTERFACE
    INTEGER FUNCTION xfgetpathdimensionality_f(GroupId, NumDims)
	  !MS$ATTRIBUTES C,reference::xfgetpathdimensionality_f
    INTEGER,              INTENT(IN)    :: GroupId
    INTEGER,              INTENT(OUT)   :: NumDims


    END FUNCTION xfgetpathdimensionality_f
  END INTERFACE

  error = xfgetpathdimensionality_f (GroupId, NumDims)


  return
END SUBROUTINE ! xfgetpathdimensionality */

! --------------------------------------------------------------------------- */
! FUNCTION  XF_GET_PATH_TIMES_ARRAY
! PURPOSE   
! NOTES     
! --------------------------------------------------------------------------- */
SUBROUTINE XF_GET_PATH_TIMES_ARRAY(GroupId, NumTimes, Times, error)
  INTEGER,              INTENT(IN)    :: GroupId
  INTEGER,              INTENT(IN)    :: NumTimes
  REAL*8, DIMENSION(*), INTENT(INOUT) :: Times
  INTEGER                                error

  
!            INTEGER, EXTERNAL :: _xfGetPathTimesArray_f
!  MS FORTRAN needs explicit interface for C functions called here.
!
  INTERFACE
    INTEGER FUNCTION xfgetpathtimesarray_f(GroupId, NumTimes, Times)
	  !MS$ATTRIBUTES C,reference::xfgetpathtimesarray_f
    INTEGER,              INTENT(IN)    :: GroupId
    INTEGER,              INTENT(IN)    :: NumTimes
    REAL*8, DIMENSION(*), INTENT(INOUT) :: Times


    END FUNCTION xfgetpathtimesarray_f
  END INTERFACE

  error = xfgetpathtimesarray_f (GroupId, NumTimes, Times)

  
END SUBROUTINE ! xfgetpathtimesarray */

! --------------------------------------------------------------------------- */
! FUNCTION  XF_READ_PATH_LOCATIONS_AT_TIME
! PURPOSE   
! NOTES     
! --------------------------------------------------------------------------- */
SUBROUTINE XF_READ_PATH_LOCATIONS_AT_TIME(GroupId, TimeIndex, FirstPathIndex, &
                                          NumIndicies, Locs, error)
  INTEGER             , INTENT(IN)  :: GroupId
  INTEGER, INTENT(IN)         :: TimeIndex, FirstPathIndex, NumIndicies
  REAL*8, DIMENSION(*), INTENT(INOUT) :: Locs
  INTEGER, INTENT(OUT)        :: error


!            INTEGER, EXTERNAL :: _xfReadPathLocationsAtTime_f
!  MS FORTRAN needs explicit interface for C functions called here.
!
  INTERFACE
    INTEGER FUNCTION xfreadpathlocationsattime_f(GroupId, TimeIndex, FirstPathIndex, &
                                                 NumIndicies, Locs)
	  !MS$ATTRIBUTES C,reference::xfreadpathlocationsattime_f
    INTEGER,              INTENT(IN)    :: GroupId
    INTEGER,              INTENT(IN)    :: TimeIndex, FirstPathIndex, NumIndicies
    REAL*8, DIMENSION(*), INTENT(INOUT) :: Locs


    END FUNCTION xfreadpathlocationsattime_f
  END INTERFACE

  error = xfreadpathlocationsattime_f (GroupId, TimeIndex, FirstPathIndex, &
                                       NumIndicies, Locs)

END SUBROUTINE ! XF_READ_PATH_LOCATIONS_AT_TIME */

! --------------------------------------------------------------------------- */
! FUNCTION  XF_READ_PATH_LOCS_FOR_PART
! PURPOSE   
! NOTES     
! --------------------------------------------------------------------------- */
SUBROUTINE XF_READ_PATH_LOCS_FOR_PART(GroupId, PathIndex,        &
                                   FirstTimeIndex, NumTimes, Locs, error)
  INTEGER             , INTENT(IN)  :: GroupId
  INTEGER, INTENT(IN)         :: PathIndex, FirstTimeIndex, NumTimes
  REAL*8, DIMENSION(*), INTENT(INOUT) :: Locs
  INTEGER, INTENT(OUT)        :: error
  
  
!            INTEGER, EXTERNAL :: _xfReadPathLocsForParticle_f
!  MS FORTRAN needs explicit interface for C functions called here.
!
  INTERFACE
    INTEGER FUNCTION xfreadpathlocsforparticle_f(GroupId, PathIndex, &
                                                      FirstTimeIndex, &
                                                      NumTimes, Locs)
	  !MS$ATTRIBUTES C,reference::xfreadpathlocsforparticle_f
    INTEGER,              INTENT(IN)    :: GroupId
    INTEGER,              INTENT(IN)    :: PathIndex, FirstTimeIndex, NumTimes
    REAL*8, DIMENSION(*), INTENT(INOUT) :: Locs


    END FUNCTION xfreadpathlocsforparticle_f
  END INTERFACE

  error = xfreadpathlocsforparticle_f (GroupId, PathIndex, &
                                            FirstTimeIndex, &
                                            NumTimes, Locs)


END SUBROUTINE ! XF_READ_PATH_LOCS_FOR_PART

! --------------------------------------------------------------------------- */
! FUNCTION  XF_READ_PATH_LOCS_FOR_PARTS
! PURPOSE   
! NOTES     
! --------------------------------------------------------------------------- */
SUBROUTINE XF_READ_PATH_LOCS_FOR_PARTS(GroupId, NumPaths,  &
                           PathIndices, FirstTimeIndex, NumTimes, Locs, error)
  INTEGER             ,        INTENT(IN)          :: GroupId
  INTEGER,               INTENT(IN)          :: NumPaths
  INTEGER, DIMENSION(*), INTENT(IN)          :: PathIndices
  INTEGER,               INTENT(IN)          :: FirstTimeIndex
  INTEGER,               INTENT(IN)          :: NumTimes
  REAL*8, DIMENSION(*),  INTENT(INOUT)       :: Locs
  INTEGER,               INTENT(OUT)         :: error
  

!            INTEGER, EXTERNAL :: _xfReadPathLocsForParticles_f
!  MS FORTRAN needs explicit interface for C functions called here.
!
  INTERFACE
    INTEGER FUNCTION xfreadpathlocsforparticles_f(GroupId, NumPaths, PathIndices, &
                                                       FirstTimeIndex, NumTimes, Locs)
	  !MS$ATTRIBUTES C,reference::xfreadpathlocsforparticles_f
    INTEGER,               INTENT(IN)    :: GroupId
    INTEGER,               INTENT(IN)    :: NumPaths, FirstTimeIndex, NumTimes
    INTEGER, DIMENSION(*), INTENT(IN)    :: PathIndices
    REAL*8, DIMENSION(*),  INTENT(INOUT) :: Locs


    END FUNCTION xfreadpathlocsforparticles_f
  END INTERFACE

  error = xfreadpathlocsforparticles_f (GroupId, NumPaths, PathIndices, &
                                             FirstTimeIndex, NumTimes, Locs)

  
  return
END SUBROUTINE ! XF_READ_PATH_LOCS_FOR_PARTS */

! ---------------------------------------------------------------------------
! FUNCTION  XF_SETUP_TO_WRITE_DATASETS
! PURPOSE   
! NOTES     
! ---------------------------------------------------------------------------
SUBROUTINE XF_SETUP_TO_WRITE_DATASETS (a_Filename, a_MultiDatasetsGroupPath, &
                        a_PathInMultiDatasetsGroup, a_SpatialDataObjectGuid, &
                        a_OverwriteOptions, a_FileId, a_GroupId, error)
CHARACTER(LEN=*), INTENT(IN)    :: a_Filename
CHARACTER(LEN=*), INTENT(IN)    :: a_MultiDatasetsGroupPath 
CHARACTER(LEN=*)                   a_PathInMultiDatasetsGroup !please don't specify an intent.
CHARACTER(LEN=*), INTENT(IN)    :: a_SpatialDataObjectGuid
INTEGER,          INTENT(IN)    :: a_OverwriteOptions
INTEGER             ,   INTENT(OUT)   :: a_FileId
INTEGER             ,   INTENT(OUT)   :: a_GroupId
INTEGER,          INTENT(OUT)   :: error
INTEGER                         :: filelen, multilen, pathlen, spatlen


!            INTEGER, EXTERNAL :: _xfSetupToWriteDatasets_f
!  MS FORTRAN needs explicit interface for C functions called here.
!
  INTERFACE
    INTEGER FUNCTION xfsetuptowritedatasets_f(a_Filename, filelen, a_MultiDatasetsGroupPath, &
                                              multilen, a_PathInMultiDatasetsGroup, pathlen, &
                                              a_SpatialDataObjectGuid, spatlen, &
                                              a_OverwriteOptions, a_FileId, a_GroupId)

	  !MS$ATTRIBUTES C,reference::xfsetuptowritedatasets_f
    !DEC$ATTRIBUTES reference :: a_Filename
    !DEC$ATTRIBUTES reference :: a_MultiDatasetsGroupPath
    !DEC$ATTRIBUTES reference :: a_PathInMultiDatasetsGroup
    !DEC$ATTRIBUTES reference :: a_SpatialDataObjectGuid
    INTEGER,              INTENT(OUT)   :: a_GroupId, a_FileId
    INTEGER,              INTENT(IN)    :: a_OverwriteOptions

    CHARACTER(LEN=*),     INTENT(IN)    :: a_Filename, a_MultiDatasetsGroupPath
    CHARACTER(LEN=*),     INTENT(IN)    :: a_SpatialDataObjectGuid
    CHARACTER(LEN=*)                    :: a_PathInMultiDatasetsGroup 
    INTEGER,              INTENT(IN)    :: filelen, multilen, pathlen, spatlen


    END FUNCTION xfsetuptowritedatasets_f
  END INTERFACE


  filelen  = LEN_TRIM(a_Filename)
  multilen = LEN_TRIM(a_MultiDatasetsGroupPath)
  pathlen  = LEN_TRIM(a_PathInMultiDatasetsGroup)
  spatlen  = LEN_TRIM(a_SpatialDataObjectGuid)
  
  error = xfsetuptowritedatasets_f (a_Filename, filelen, a_MultiDatasetsGroupPath, multilen, &
                                    a_PathInMultiDatasetsGroup, pathlen, &
                                    a_SpatialDataObjectGuid, spatlen, &
                                    a_OverwriteOptions, a_FileId, a_GroupId)


return

END SUBROUTINE ! xfsetuptowritedatasets
! ---------------------------------------------------------------------------
! FUNCTION  xfcreatemultidatasetsgroup
! PURPOSE   
! NOTES     
! ---------------------------------------------------------------------------
SUBROUTINE XF_CREATE_MULTI_DATASETS_GROUP (a_Id, a_Path, a_Guid, a_MultiId, error)
INTEGER             ,   INTENT(IN)  :: a_Id
CHARACTER(LEN=*), INTENT(IN)  :: a_Path, a_Guid
INTEGER             ,   INTENT(OUT) :: a_MultiId
INTEGER, INTENT(OUT)          :: error
INTEGER                          pathlen, guidlen


!            INTEGER, EXTERNAL :: _xfCreateMultiDatasetsGroup_f
!  MS FORTRAN needs explicit interface for C functions called here.
!
  INTERFACE
    INTEGER FUNCTION xfcreatemultidatasetsgroup_f(a_Id, a_Path, pathlen, a_Guid, &
                                                  guidlen, a_MultiId)
	  !MS$ATTRIBUTES C,reference::xfcreatemultidatasetsgroup_f
    !DEC$ATTRIBUTES reference :: a_Path
    !DEC$ATTRIBUTES reference :: a_Guid
    INTEGER,              INTENT(IN)    :: a_Id
    INTEGER,              INTENT(OUT)   :: a_MultiId
    CHARACTER(LEN=*),     INTENT(IN)    :: a_Path, a_Guid
    INTEGER                             :: pathlen, guidlen
    


    END FUNCTION xfcreatemultidatasetsgroup_f
  END INTERFACE

  pathlen  = LEN_TRIM(a_Path)
  guidlen  = LEN_TRIM(a_Guid)
  error = xfcreatemultidatasetsgroup_f (a_Id, a_Path, pathlen, a_Guid, &
                                        guidlen, a_MultiId)


return

END SUBROUTINE ! xfcreatemultidatasetsgroup
!******************************
! ---------------------------------------------------------------------------
! FUNCTION  xfgetgrouppathssizeformultidatasets
! PURPOSE   
! NOTES     
! ---------------------------------------------------------------------------
SUBROUTINE XF_GET_GRP_PTHS_SZ_MLT_DSETS (a_Id, Num, Maxsize, error)
INTEGER             , INTENT(IN)    :: a_Id
INTEGER                          Num
INTEGER                          Maxsize
INTEGER,          INTENT(OUT) :: error

!            INTEGER, EXTERNAL :: _xfGetGrpPathsSizeForMltDsets_f
!  MS FORTRAN needs explicit interface for C functions called here.
!
  INTERFACE
    INTEGER FUNCTION xfgetgrppathssizeformltdsets_f(a_Id, Num, Maxsize)

	  !MS$ATTRIBUTES C,reference::xfgetgrppathssizeformltdsets_f
    INTEGER,              INTENT(IN)    :: a_Id
    INTEGER,              INTENT(IN)    :: Num, Maxsize
    


    END FUNCTION xfgetgrppathssizeformltdsets_f
  END INTERFACE


  error = xfgetgrppathssizeformltdsets_f (a_Id, Num, Maxsize)


return 
END SUBROUTINE ! xfgetgrouppathssizeformultidatasets
!******************************
! ---------------------------------------------------------------------------
! FUNCTION  xfgetallgrouppathsformultidatasets
! PURPOSE   
! NOTES     
! ---------------------------------------------------------------------------
SUBROUTINE XF_GET_ALL_GRP_PATHS_MLT_DSETS (a_Id, a_Num, a_Maxsize, &
                                           a_Paths, error)
INTEGER             , INTENT(IN)    :: a_Id
INTEGER                          a_Num
INTEGER                          a_Maxsize
CHARACTER,     DIMENSION(*)   :: a_Paths 
INTEGER,       INTENT(OUT)    :: error
INTEGER                          pathlen


!            INTEGER, EXTERNAL :: _xfGetAllGrpPathsForMltDsets_f
!  MS FORTRAN needs explicit interface for C functions called here.
!
  INTERFACE
    INTEGER FUNCTION xfgetallgrppathsformltdsets_f(a_Id, a_Num, a_Maxsize, &
                                                          a_Paths, pathlen)

	  !MS$ATTRIBUTES C,reference::xfgetallgrppathsformltdsets_f
    !DEC$ATTRIBUTES reference :: a_Paths
    INTEGER,              INTENT(IN)    ::  a_Id
    INTEGER,              INTENT(IN)    ::  a_Num, a_Maxsize
    CHARACTER,            DIMENSION(*)  ::  a_Paths
    INTEGER,              INTENT(IN)    ::  pathlen
    


    END FUNCTION xfgetallgrppathsformltdsets_f
  END INTERFACE


  pathlen = LEN(a_Paths)
  error = xfgetallgrppathsformltdsets_f (a_Id, a_Num, a_Maxsize, &
                                                a_Paths, pathlen)



return
END SUBROUTINE ! xfgetallgrouppathsformultidatasets
! ---------------------------------------------------------------------------
! FUNCTION  xfgetdatasetssdoguid
! PURPOSE   
! NOTES     
! ---------------------------------------------------------------------------
SUBROUTINE XF_GET_DATASETS_SDO_GUID (a_MultiDatasetsGroup, a_GUID, error)
INTEGER             , INTENT(IN) :: a_MultiDatasetsGroup
CHARACTER(LEN=*)              a_GUID
INTEGER, INTENT(OUT)       :: error
INTEGER                       a_guidlen

!            INTEGER, EXTERNAL :: _xfGetDatasetsSdoGuid_f
!  MS FORTRAN needs explicit interface for C functions called here.
!
  INTERFACE
    INTEGER FUNCTION xfgetdatasetssdoguid_f(a_MultiDatasetsGroup, a_GUID, a_guidlen)

	  !MS$ATTRIBUTES C,reference::xfgetdatasetssdoguid_f
    !DEC$ATTRIBUTES reference :: a_GUID
    INTEGER,              INTENT(IN)    :: a_MultiDatasetsGroup
    CHARACTER(LEN=*)                    :: a_GUID
    INTEGER,              INTENT(IN)    :: a_guidlen
    


    END FUNCTION xfgetdatasetssdoguid_f
  END INTERFACE


  a_guidlen = LEN_TRIM(a_GUID)
  error = xfgetdatasetssdoguid_f (a_MultiDatasetsGroup, a_GUID, a_guidlen)


return 
END SUBROUTINE ! xfgetdatasetssdoguid
! ---------------------------------------------------------------------------
! FUNCTION  XF_OPEN_MULTI_DATASETS_GROUP
! PURPOSE   Open or create a multi-datasets group inside a mesh or grid
! NOTES     
! ---------------------------------------------------------------------------
SUBROUTINE XF_OPEN_MULTI_DATASETS_GROUP (a_Id, DatasetsGroupId, error)
INTEGER             , INTENT(IN)  :: a_Id
INTEGER             , INTENT(OUT) :: DatasetsGroupId
INTEGER, INTENT(OUT)        :: error


!            INTEGER, EXTERNAL :: _xfOpenMultiDatasetsGroup_f
!  MS FORTRAN needs explicit interface for C functions called here.
!
  INTERFACE
    INTEGER FUNCTION xfopenmultidatasetsgroup_f(a_Id, DatasetsGroupId)

	  !MS$ATTRIBUTES C,reference::xfopenmultidatasetsgroup_f
    INTEGER,              INTENT(IN)    :: a_Id
    INTEGER,              INTENT(OUT)   :: DatasetsGroupId
    

    END FUNCTION xfopenmultidatasetsgroup_f
  END INTERFACE


  error = xfopenmultidatasetsgroup_f (a_Id, DatasetsGroupId)


return 
END SUBROUTINE ! xfopenmultidatasetsgroup

!******************************

SUBROUTINE XF_CREATE_SCALAR_DATASET (a_DatasetsGroupId, a_Path, a_Units,&
                                 a_TimeUnits, a_Compression, a_DatasetId, error)
INTEGER, INTENT(IN)          :: a_DatasetsGroupId
CHARACTER(LEN=*), INTENT(IN) :: a_Path, a_Units, a_TimeUnits
INTEGER, INTENT(IN)          :: a_Compression
INTEGER             , INTENT(OUT)  :: a_DatasetId
INTEGER, INTENT(OUT)         :: error
INTEGER                         namelen1, namelen2, namelen3

!            INTEGER, EXTERNAL :: _xfCreateFile_f
!  MS FORTRAN needs explicit interface for C functions called here.

  INTERFACE
    INTEGER FUNCTION xfcreatescalardataset_f(a_DatasetsGroupId, a_Path, namelen1, a_Units,&
	                                         namelen2, a_TimeUnits, namelen3, a_Compression,&
											 a_DatasetId)
	  !MS$ATTRIBUTES C,reference::xfcreatescalardataset_f
	  !DEC$ATTRIBUTES reference :: a_Path
	  !DEC$ATTRIBUTES reference :: a_Units
	  !DEC$ATTRIBUTES reference :: a_TimeUnits
	  INTEGER, INTENT(IN)          :: a_DatasetsGroupId
     CHARACTER(LEN=*), INTENT(IN) :: a_Path, a_Units, a_TimeUnits
	  INTEGER, INTENT(IN)          :: a_Compression
	  INTEGER, INTENT(OUT)         :: a_DatasetId
	  INTEGER, INTENT(IN)          :: namelen1, namelen2, namelen3 

    END FUNCTION xfcreatescalardataset_f
  END INTERFACE

  namelen1 = LEN_TRIM(a_Path)
  namelen2 = LEN_TRIM(a_Units)
  namelen3 = LEN_TRIM(a_TimeUnits)

  error = xfcreatescalardataset_f(a_DatasetsGroupId, a_Path, namelen1, a_Units,&
	                              namelen2, a_TimeUnits, namelen3, a_Compression,&
								  a_DatasetId)
  return

END SUBROUTINE


!********************************

SUBROUTINE XF_CREATE_VECTOR_DATASET (a_DatasetsGroupId, a_Path, a_Units,&
                                     a_TimeUnits, a_Compression, a_DatasetId, error)

INTEGER, INTENT(IN)          :: a_DatasetsGroupId
CHARACTER(LEN=*), INTENT(IN) :: a_Path, a_Units, a_TimeUnits
INTEGER, INTENT(IN)          :: a_Compression
INTEGER             , INTENT(OUT)  :: a_DatasetId
INTEGER, INTENT(OUT)         :: error
INTEGER                      namelen1, namelen2, namelen3

!            INTEGER, EXTERNAL :: _xfCreateFile_f
!  MS FORTRAN needs explicit interface for C functions called here.

  INTERFACE
    INTEGER FUNCTION xfcreatevectordataset_f(a_DatasetsGroupId, a_Path, namelen1, a_Units,&
	                                         namelen2, a_TimeUnits, namelen3, a_Compression,&
											 a_DatasetId)
	  !MS$ATTRIBUTES C,reference::xfcreatevectordataset_f
	  !DEC$ATTRIBUTES reference :: a_Path
	  !DEC$ATTRIBUTES reference :: a_Units
	  !DEC$ATTRIBUTES reference :: a_TimeUnits
	  INTEGER, INTENT(IN)          :: a_DatasetsGroupId
     CHARACTER(LEN=*), INTENT(IN) :: a_Path, a_Units, a_TimeUnits
	  INTEGER, INTENT(IN)          :: a_Compression
	  INTEGER, INTENT(OUT)         :: a_DatasetId
	  INTEGER, INTENT(IN)          :: namelen1, namelen2, namelen3 

    END FUNCTION xfcreatevectordataset_f
  END INTERFACE

  namelen1 = LEN_TRIM(a_Path)
  namelen2 = LEN_TRIM(a_Units)
  namelen3 = LEN_TRIM(a_TimeUnits)

  error = xfcreatevectordataset_f(a_DatasetsGroupId, a_Path, namelen1, a_Units,&
	                              namelen2, a_TimeUnits, namelen3, a_Compression,&
								                a_DatasetId)
return


END SUBROUTINE

!* FUNCTION  XF_CREATE_SCALAR_DSET_EXTNDBL
!* PURPOSE   Create a scalar dataset that can be extended
!* NOTES     The intermediate groups in the path may or may not be created
!******************2***********************************************************/
SUBROUTINE XF_CREATE_SCALAR_DSET_EXTNDBL(a_DatasetsGroupId, a_Path, &
                            a_Units, a_TimeUnits, a_FillVal, a_Compression, &
                            a_DatasetId, error)
 
  INTEGER             ,   INTENT(IN)    :: a_DatasetsGroupId
  CHARACTER(LEN=*), INTENT(IN)    :: a_Path, a_Units, a_TimeUnits
  REAL*4,             INTENT(IN)    :: a_FillVal
  INTEGER,          INTENT(IN)    :: a_Compression
  INTEGER             ,   INTENT(OUT)   :: a_DatasetId 
  INTEGER,          INTENT(OUT)   :: error
  INTEGER                            pathlen, unitlen, timelen

  
!            INTEGER, EXTERNAL :: _xfCreateScalarDsetExtndbl_f
!  MS FORTRAN needs explicit interface for C functions called here.

  INTERFACE
    INTEGER FUNCTION xfcreatescalardsetextndbl_f(a_DatasetsGroupId, a_Path, pathlen, &
                                                      a_Units, unitlen, a_TimeUnits, &
                                                      timelen, a_FillVal, a_Compression, & 
                                                      a_DatasetId)
	  !MS$ATTRIBUTES C,reference::xfcreatescalardsetextndbl_f
	  !DEC$ATTRIBUTES reference :: a_Path
	  !DEC$ATTRIBUTES reference :: a_Units
	  !DEC$ATTRIBUTES reference :: a_TimeUnits
	  INTEGER, INTENT(IN)           :: a_DatasetsGroupId
    CHARACTER(LEN=*), INTENT(IN)  :: a_Path, a_Units, a_TimeUnits
	  INTEGER, INTENT(IN)           :: a_Compression
	  INTEGER, INTENT(OUT)          :: a_DatasetId
	  INTEGER, INTENT(IN)           :: pathlen, unitlen, timelen
    REAL*4,    INTENT(IN)           :: a_FillVal

    END FUNCTION xfcreatescalardsetextndbl_f
  END INTERFACE

  pathlen = LEN_TRIM(a_Path)
  unitlen = LEN_TRIM(a_Units)
  timelen = LEN_TRIM(a_TimeUnits)

  error = xfcreatescalardsetextndbl_f (a_DatasetsGroupId, a_Path, pathlen, &
                                             a_Units, unitlen, a_TimeUnits, &
                                             timelen, a_FillVal, a_Compression, & 
                                             a_DatasetId)
  
  return
END SUBROUTINE

! FUNCTION  xfextendscalardataset
! PURPOSE   Create a scalar dataset that can be extended
! NOTES     The intermediate groups in the path may or may not be created
SUBROUTINE XF_EXTEND_SCALAR_DATASET (a_Id, a_NewSize, error)
  INTEGER             , INTENT(IN)  :: a_Id
  INTEGER             , INTENT(IN)  :: a_NewSize
  INTEGER             , INTENT(OUT) :: error


!            INTEGER, EXTERNAL :: _xfExtendScalarDataset_f
!  MS FORTRAN needs explicit interface for C functions called here.

  INTERFACE
    INTEGER FUNCTION xfextendscalardataset_f(a_Id, a_NewSize)
	  !MS$ATTRIBUTES C,reference::xfextendscalardataset_f
	  INTEGER, INTENT(IN)           :: a_Id, a_NewSize

    END FUNCTION xfextendscalardataset_f
  END INTERFACE

  
  error = xfextendscalardataset_f (a_Id, a_NewSize)
  
  return
END SUBROUTINE

!********************************

SUBROUTINE XF_WRITE_SCALAR_TIMESTEP (a_Id, a_Time, a_NumValues, a_Values, error)


INTEGER,              INTENT(IN)   :: a_Id
REAL*8,               INTENT(IN)   :: a_Time
INTEGER,              INTENT(IN)   :: a_NumValues
REAL*4, DIMENSION(*),   INTENT(IN)   :: a_Values
INTEGER,              INTENT(OUT)  :: error

!            INTEGER, EXTERNAL :: _xfCreateFile_f
!  MS FORTRAN needs explicit interface for C functions called here.

  INTERFACE
    INTEGER FUNCTION xfwritescalartimestep_f(a_Id, a_Time, a_NumValues, a_Values)
	  !MS$ATTRIBUTES C,reference::xfwritescalartimestep_f
	  INTEGER,               INTENT(IN)  :: a_Id
      REAL*8,                INTENT(IN)  :: a_Time
      INTEGER,               INTENT(IN)  :: a_NumValues
	  REAL*4,    DIMENSION(*), INTENT(IN)  :: a_Values

    END FUNCTION xfwritescalartimestep_f
  END INTERFACE


  error = xfwritescalartimestep_f(a_Id, a_Time, a_NumValues, a_Values)

return

END SUBROUTINE

!********************************

SUBROUTINE XF_SET_DATASET_TIMESTEP_MIN_MAX (a_Id, a_TimeId, a_Minvalue, &
                                            a_Maxvalue, error)


INTEGER,              INTENT(IN)   :: a_Id
INTEGER,              INTENT(IN)   :: a_TimeId
REAL*4,               INTENT(IN)   :: a_Minvalue, a_Maxvalue
INTEGER,              INTENT(OUT)  :: error

!            INTEGER, EXTERNAL :: _xfCreateFile_f
!  MS FORTRAN needs explicit interface for C functions called here.

  INTERFACE
    INTEGER FUNCTION xfsetdatasettimestepminmax_f(a_Id, a_TimeId, &
                                                  a_Minvalue, a_Maxvalue)
	  !MS$ATTRIBUTES C,reference::xfsetdatasettimestepminmax_f
	  INTEGER,               INTENT(IN)  :: a_Id
      INTEGER,               INTENT(IN)  :: a_TimeId
      REAL*4,               INTENT(IN)   :: a_Minvalue, a_Maxvalue

    END FUNCTION xfsetdatasettimestepminmax_f
  END INTERFACE


  error = xfsetdatasettimestepminmax_f(a_Id, a_TimeId, a_Minvalue, a_Maxvalue)

return

END SUBROUTINE

!********************************

SUBROUTINE XF_WRITE_SCALAR_TS_MIN_MAX(a_Id, a_Time, a_NumValues, &
                                      a_Values, a_Minimum, a_Maximum, error)
INTEGER ,                 INTENT(IN)    ::  a_Id
REAL*8,                   INTENT(IN)    ::  a_Time
INTEGER,                  INTENT(IN)    ::  a_NumValues
REAL*4,    DIMENSION(*),    INTENT(IN)    ::  a_Values
REAL*4,                     INTENT(IN)    ::  a_Minimum, a_Maximum
INTEGER,                  INTENT(OUT)   ::  error


!            INTEGER, EXTERNAL :: _xfWriteScalarTimestepMinMax_f
!  MS FORTRAN needs explicit interface for C functions called here.

  INTERFACE
    INTEGER FUNCTION xfwritescalartimestepminmax_f(a_Id, a_Time, a_NumValues, &
                                                   a_Values, a_Minimum, a_Maximum)
	  !MS$ATTRIBUTES C,reference::xfwritescalartimestepminmax_f
	  INTEGER,               INTENT(IN)   :: a_Id, a_NumValues
      REAL*8,                INTENT(IN)   :: a_Time
	  REAL*4,    DIMENSION(*), INTENT(IN)   :: a_Values
      REAL*4,                  INTENT(IN)   :: a_Minimum, a_Maximum


    END FUNCTION xfwritescalartimestepminmax_f
  END INTERFACE


  error = xfwritescalartimestepminmax_f(a_Id, a_Time, a_NumValues, &
                                        a_Values, a_Minimum, a_Maximum)


return

END SUBROUTINE

!********************************

SUBROUTINE XF_INITIALIZE_SCALAR_TIMESTEP (a_Id, a_Time, a_NumValues, a_Minimum,&
                                          a_Maximum, a_timestepId, error)

INTEGER,              INTENT(IN)   :: a_Id
REAL*8,               INTENT(IN)   :: a_Time
INTEGER,              INTENT(IN)   :: a_NumValues
REAL*4,               INTENT(IN)   :: a_Minimum, a_Maximum
INTEGER,              INTENT(OUT)  :: a_timestepId
INTEGER,              INTENT(OUT)  :: error

!            INTEGER, EXTERNAL :: _xfInitializeScalarTimestep_f
!  MS FORTRAN needs explicit interface for C functions called here.

  INTERFACE
    INTEGER FUNCTION xfinitializescalartimestep_f(a_Id, a_Time, a_NumValues, &
                       a_Minimum, a_Maximum, a_timestepId)
	  !MS$ATTRIBUTES C,reference::xfinitializescalartimestep_f
	  INTEGER,               INTENT(IN)  :: a_Id
      REAL*8,                INTENT(IN)  :: a_Time
      INTEGER,               INTENT(IN)  :: a_NumValues
      REAL*4,               INTENT(IN)   :: a_Minimum, a_Maximum
      INTEGER,              INTENT(OUT)  :: a_timestepId

    END FUNCTION xfinitializescalartimestep_f
  END INTERFACE


  error = xfinitializescalartimestep_f(a_Id, a_Time, a_NumValues, a_Minimum, &
                                  a_Maximum, a_timestepId)

return

END SUBROUTINE

!********************************

SUBROUTINE XF_WRITE_SCALAR_TIMESTEP_PORTION (a_Id, a_timeStepId, &
                    a_NumValuesToWrite, a_startIndex, a_Values, error)

INTEGER,              INTENT(IN)   :: a_Id
INTEGER,              INTENT(IN)   :: a_timeStepId
INTEGER,              INTENT(IN)   :: a_NumValuesToWrite
INTEGER,              INTENT(IN)   :: a_startIndex
REAL*4, DIMENSION(*),   INTENT(IN)   :: a_Values
INTEGER,              INTENT(OUT)  :: error

!            INTEGER, EXTERNAL :: _xfWriteScalarTimestepPortion_f
!  MS FORTRAN needs explicit interface for C functions called here.

  INTERFACE
    INTEGER FUNCTION xfwritescalartimestepportion_f(a_Id, a_timeStepId, &
                    a_NumValuesToWrite, a_startIndex, a_Values)
	  !MS$ATTRIBUTES C,reference::xfwritescalartimestepportion_f
	  INTEGER,               INTENT(IN)  :: a_Id
      INTEGER,              INTENT(IN)   :: a_timeStepId
      INTEGER,              INTENT(IN)   :: a_NumValuesToWrite
      INTEGER,              INTENT(IN)   :: a_startIndex
      REAL*4, DIMENSION(*),   INTENT(IN)   :: a_Values

    END FUNCTION xfwritescalartimestepportion_f
  END INTERFACE


  error = xfwritescalartimestepportion_f(a_Id, a_timeStepId, &
                    a_NumValuesToWrite, a_startIndex, a_Values)

return

END SUBROUTINE

!********************************

SUBROUTINE XF_WRITE_VECTOR_TIMESTEP (a_Id, a_Time, a_NumValues, a_NumComponents,&
                                     a_Values, error)
INTEGER,        INTENT(IN)                                          :: a_Id
REAL*8,         INTENT(IN)                                          :: a_Time
INTEGER,        INTENT(IN)                                          :: a_NumValues
INTEGER,        INTENT(IN)                                          :: a_NumComponents
REAL*4,           INTENT(IN), DIMENSION(a_NumComponents, a_NumValues) :: a_Values
INTEGER,        INTENT(OUT)                                         :: error


!            INTEGER, EXTERNAL :: _xfWriteVectorTimestep_f
!  MS FORTRAN needs explicit interface for C functions called here.

  INTERFACE
    INTEGER FUNCTION xfwritevectortimestep_f(a_Id, a_Time, a_NumValues, a_NumComponents, &
                                             a_Values)
	  !MS$ATTRIBUTES C,reference::xfwritevectortimestep_f
	  INTEGER, INTENT(IN)   :: a_Id, a_NumValues, a_NumComponents
      REAL*8,  INTENT(IN)   :: a_Time
      REAL*4,    INTENT(IN), DIMENSION(a_NumComponents, a_NumValues) :: a_Values


    END FUNCTION xfwritevectortimestep_f
  END INTERFACE


  error = xfwritevectortimestep_f (a_Id, a_Time, a_NumValues, a_NumComponents, a_Values)



return
END SUBROUTINE

!********************************

SUBROUTINE XF_WRITE_VECTOR_TS_MIN_MAX(a_Id, a_Time, a_NumValues, a_NumComponents, &
                                            a_Values, a_Minimum, a_Maximum, error)
INTEGER,                                        INTENT(IN)  :: a_Id
REAL*8,                                         INTENT(IN)  :: a_Time
INTEGER,                                        INTENT(IN)  :: a_NumValues, a_NumComponents
REAL*4, DIMENSION(a_NumComponents, a_NumValues),  INTENT(IN)  :: a_Values
REAL*4,                                           INTENT(IN)  :: a_Minimum, a_Maximum
INTEGER,                                        INTENT(OUT) :: error


!            INTEGER, EXTERNAL :: _xfWriteVectorTimestepMinMax_f
!  MS FORTRAN needs explicit interface for C functions called here.

  INTERFACE
    INTEGER FUNCTION xfwritevectortimestepminmax_f(a_Id, a_Time, a_NumValues, a_NumComponents, &
                                                   a_Values, a_Minimum, a_Maximum)
	  !MS$ATTRIBUTES C,reference::xfwritevectortimestepminmax_f
	  INTEGER,              INTENT(IN)  :: a_Id, a_NumValues, a_NumComponents
      REAL*8,               INTENT(IN)  :: a_Time
      REAL*4, DIMENSION(a_NumComponents, a_NumValues),  INTENT(IN)  :: a_Values
      REAL*4,                 INTENT(IN)  :: a_Minimum, a_Maximum


    END FUNCTION xfwritevectortimestepminmax_f
  END INTERFACE


  error = xfwritevectortimestepminmax_f (a_Id, a_Time, a_NumValues, a_NumComponents, &
                                         a_Values, a_Minimum, a_Maximum)


return

END SUBROUTINE
       
       !********************************

SUBROUTINE XF_INITIALIZE_VECTOR_TIMESTEP(a_Id, a_Time, a_NumValues, a_NumComponents, &
                                         a_Minimum, a_Maximum, timeId, error)
INTEGER,                                        INTENT(IN)  :: a_Id
REAL*8,                                         INTENT(IN)  :: a_Time
INTEGER,                                        INTENT(IN)  :: a_NumValues, a_NumComponents
REAL*4,                                           INTENT(IN)  :: a_Minimum, a_Maximum
INTEGER,                                        INTENT(OUT) :: timeId
INTEGER,                                        INTENT(OUT) :: error


!            INTEGER, EXTERNAL :: _xfInitializeVectorTimestep_f
!  MS FORTRAN needs explicit interface for C functions called here.

  INTERFACE
    INTEGER FUNCTION xfinitializevectortimestep_f(a_Id, a_Time, a_NumValues, a_NumComponents, &
                                                  a_Minimum, a_Maximum, timeId)
	  !MS$ATTRIBUTES C,reference::xfinitializevectortimestep_f
	  INTEGER,              INTENT(IN)  :: a_Id, a_NumValues, a_NumComponents
      REAL*8,               INTENT(IN)  :: a_Time
      REAL*4,               INTENT(IN)  :: a_Minimum, a_Maximum
      INTEGER,              INTENT(OUT) :: timeId


    END FUNCTION xfinitializevectortimestep_f
  END INTERFACE


  error = xfinitializevectortimestep_f (a_Id, a_Time, a_NumValues, a_NumComponents, &
                                        a_Minimum, a_Maximum, timeId)

return

END SUBROUTINE

!********************************

SUBROUTINE XF_WRITE_VECTOR_TIMESTEP_PORTION(a_Id, a_TimeId, a_NumValuesToWrite, &
              a_NumComponentsToWrite, a_startIndex, a_startComponent, a_Values, error)
INTEGER,                                        INTENT(IN)  :: a_Id
INTEGER,                                         INTENT(IN)  :: a_TimeId
INTEGER,                                        INTENT(IN)  :: a_NumValuesToWrite, a_NumComponentsToWrite
INTEGER,                                        INTENT(IN)  :: a_startIndex, a_startComponent
REAL*4, DIMENSION(*),  INTENT(IN)                           :: a_Values
INTEGER,                                        INTENT(OUT) :: error


!            INTEGER, EXTERNAL :: _xfWriteVectorTimestepPortion_f
!  MS FORTRAN needs explicit interface for C functions called here.

  INTERFACE
    INTEGER FUNCTION xfwritevectortimestepportion_f(a_Id, a_TimeId, a_NumValuesToWrite, &
              a_NumComponentsToWrite, a_startIndex, a_startComponent, a_Values)
	  !MS$ATTRIBUTES C,reference::xfwritevectortimestepportion_f
        INTEGER,                                        INTENT(IN)  :: a_Id
        INTEGER,                                         INTENT(IN)  :: a_TimeId
        INTEGER,                                        INTENT(IN)  :: a_NumValuesToWrite, a_NumComponentsToWrite
        INTEGER,                                        INTENT(IN) :: a_startIndex, a_startComponent
        REAL*4, DIMENSION(*),  INTENT(IN)                           :: a_Values


    END FUNCTION xfwritevectortimestepportion_f
  END INTERFACE


  error = xfwritevectortimestepportion_f (a_Id, a_TimeId, a_NumValuesToWrite, &
              a_NumComponentsToWrite, a_startIndex, a_startComponent, a_Values)


return

END SUBROUTINE

       
!***********************

SUBROUTINE XF_WRITE_ACTIVITY_TIMESTEP (a_Id, a_NumActive, a_Active, error)
INTEGER             ,        INTENT(IN)  :: a_Id
INTEGER,               INTENT(IN)  :: a_NumActive
!INTEGER, DIMENSION(*), INTENT(IN)  :: a_Active
INTEGER*1, dimension(a_numactive), INTENT(IN)  :: a_Active
INTEGER,               INTENT(OUT) :: error

!            INTEGER, EXTERNAL :: _xfWriteActivityTimestep_f
!  MS FORTRAN needs explicit interface for C functions called here.

  INTERFACE
    INTEGER FUNCTION xfwriteactivitytimestep_f(a_Id, a_NumActive, a_Active)

	  !MS$ATTRIBUTES C,reference::xfwriteactivitytimestep_f
	  INTEGER,               INTENT(IN)   :: a_Id
    INTEGER,               INTENT(IN)   :: a_NumActive
    !INTEGER, DIMENSION(*), INTENT(IN)   :: a_Active
    INTEGER*1, dimension(a_numactive), INTENT(IN)   :: a_Active


    END FUNCTION xfwriteactivitytimestep_f
  END INTERFACE


  error = xfwriteactivitytimestep_f (a_Id, a_NumActive, a_Active)


return

END SUBROUTINE

!********************************

SUBROUTINE XF_INITIALIZE_ACTIVITY_TIMESTEP (a_Id, a_NumActive, a_timestepId, error)

INTEGER,              INTENT(IN)   :: a_Id
INTEGER,              INTENT(IN)   :: a_NumActive
INTEGER,              INTENT(OUT)  :: a_timestepId
INTEGER,              INTENT(OUT)  :: error

!            INTEGER, EXTERNAL :: _xfInitializeActivityTimestep_f
!  MS FORTRAN needs explicit interface for C functions called here.

  INTERFACE
    INTEGER FUNCTION xfinitializeactivitytimestep_f(a_Id, a_NumActive, &
                       a_timestepId)
	  !MS$ATTRIBUTES C,reference::xfinitializeactivitytimestep_f
	  INTEGER,               INTENT(IN)  :: a_Id
      INTEGER,               INTENT(IN)  :: a_NumActive
      INTEGER,              INTENT(OUT)  :: a_timestepId

    END FUNCTION xfinitializeactivitytimestep_f
  END INTERFACE


  error = xfinitializeactivitytimestep_f(a_Id, a_NumActive, a_timestepId)

return

END SUBROUTINE

!********************************

SUBROUTINE XF_WRITE_ACTIVITY_TIMESTEP_PORTION (a_Id, a_timeStepId, &
                    a_NumValuesToWrite, a_startIndex, a_ActiveValues, error)
INTEGER,              INTENT(IN)   :: a_Id
INTEGER,              INTENT(IN)   :: a_timeStepId
INTEGER,              INTENT(IN)   :: a_NumValuesToWrite
INTEGER,              INTENT(IN)   :: a_startIndex
INTEGER*1, dimension(*), INTENT(IN)  :: a_ActiveValues
INTEGER,              INTENT(OUT)  :: error

!            INTEGER, EXTERNAL :: _xfwriteactivitytimestepportion_f
!  MS FORTRAN needs explicit interface for C functions called here.

  INTERFACE
    INTEGER FUNCTION xfwriteactivitytimestepportion_f(a_Id, a_timeStepId, &
                    a_NumValuesToWrite, a_startIndex, a_ActiveValues)
	  !MS$ATTRIBUTES C,reference::xfwriteactivitytimestepportion_f
	  INTEGER,               INTENT(IN)  :: a_Id
      INTEGER,              INTENT(IN)   :: a_timeStepId
      INTEGER,              INTENT(IN)   :: a_NumValuesToWrite
      INTEGER,              INTENT(IN)   :: a_startIndex
      INTEGER*1, dimension(*), INTENT(IN)  :: a_ActiveValues

    END FUNCTION xfwriteactivitytimestepportion_f
  END INTERFACE


  error = xfwriteactivitytimestepportion_f(a_Id, a_timeStepId, &
                    a_NumValuesToWrite, a_startIndex, a_ActiveValues)

return

END SUBROUTINE

!***********************

SUBROUTINE XF_SCALAR_DATASET_GROUP_ID (a_Id, error)
INTEGER,        INTENT(IN)  :: a_Id
INTEGER,        INTENT(OUT) :: error

!            INTEGER, EXTERNAL :: _xfGetScalarDatasetGroupId_f
!  MS FORTRAN needs explicit interface for C functions called here.

  INTERFACE
    INTEGER FUNCTION xfgetscalardatasetgroupid_f(a_Id)

	  !MS$ATTRIBUTES C,reference::xfgetscalardatasetgroupid_f
	  INTEGER,            INTENT(IN)   :: a_Id


    END FUNCTION xfgetscalardatasetgroupid_f
  END INTERFACE


  error = xfgetscalardatasetgroupid_f (a_Id)


return

END SUBROUTINE

!***********************

SUBROUTINE XF_VECTOR_DATASET_GROUP_ID (a_Id, error)
INTEGER             , INTENT(IN)  :: a_Id
INTEGER,        INTENT(OUT) :: error

!            INTEGER, EXTERNAL :: _xfGetVectorDatasetGroupId_f
!  MS FORTRAN needs explicit interface for C functions called here.

  INTERFACE
    INTEGER FUNCTION xfgetvectordatasetgroupid_f(a_Id)

	  !MS$ATTRIBUTES C,reference::xfgetvectordatasetgroupid_f
	  INTEGER,            INTENT(IN)   :: a_Id


    END FUNCTION xfgetvectordatasetgroupid_f
  END INTERFACE


  error = xfgetvectordatasetgroupid_f (a_Id)


return

END SUBROUTINE

!***********************

SUBROUTINE XF_WRITE_REFTIME (a_Id, a_Reftime, error)
INTEGER,  INTENT(IN)   :: a_Id
REAL*8,   INTENT(IN)   :: a_Reftime
INTEGER,  INTENT(OUT)  :: error


!            INTEGER, EXTERNAL :: _xfWriteReftime_f
!  MS FORTRAN needs explicit interface for C functions called here.

  INTERFACE
    INTEGER FUNCTION xfwritereftime_f(a_Id, a_Reftime)
	  !MS$ATTRIBUTES C,reference::xfwritereftime_f
	  INTEGER,            INTENT(IN)   :: a_Id
      REAL*8, INTENT(IN)   :: a_Reftime


    END FUNCTION xfwritereftime_f
  END INTERFACE


  error = xfwritereftime_f (a_Id, a_Reftime)


return

END SUBROUTINE

!***********************

SUBROUTINE XF_USE_REFTIME (a_Id, a_bUseReftime, error)
INTEGER             , INTENT(IN)    :: a_Id
LOGICAL*2,        INTENT(OUT)   :: a_bUseReftime
INTEGER,        INTENT(OUT)   :: error


!            INTEGER, EXTERNAL :: _xfUseReftime_f
!  MS FORTRAN needs explicit interface for C functions called here.

  INTERFACE
    INTEGER FUNCTION xfusereftime_f(a_Id, a_bUseReftime)
	  !MS$ATTRIBUTES C,reference::xfusereftime_f
	  INTEGER,      INTENT(IN)   :: a_Id
    LOGICAL*2,      INTENT(OUT)  :: a_bUseReftime


    END FUNCTION xfusereftime_f
  END INTERFACE

  error = xfusereftime_f (a_Id, a_bUseReftime)


return

END SUBROUTINE

!***********************

SUBROUTINE XF_READ_REFTIME (a_Id, a_dReftime, error)
INTEGER,        INTENT(IN)   :: a_Id
REAL*8,         INTENT(OUT)  :: a_dReftime
INTEGER,        INTENT(OUT)  :: error


!            INTEGER, EXTERNAL :: _xfReadReftime_f
!  MS FORTRAN needs explicit interface for C functions called here.

  INTERFACE
    INTEGER FUNCTION xfreadreftime_f(a_Id, a_dReftime)
	  !MS$ATTRIBUTES C,reference::xfreadreftime_f
	  INTEGER,              INTENT(IN)   :: a_Id
      REAL*8,               INTENT(OUT)  :: a_dReftime


    END FUNCTION xfreadreftime_f
  END INTERFACE


  error = xfreadreftime_f (a_Id, a_dReftime)


return

END SUBROUTINE

!***********************

SUBROUTINE XF_GET_SCALAR_DATASETS_INFO (a_Id, a_Number, a_MaxPathLength, error)
INTEGER             , INTENT(IN)  :: a_Id
INTEGER                     :: a_Number
INTEGER,        INTENT(OUT) :: a_MaxPathLength
INTEGER,        INTENT(OUT) :: error


!            INTEGER, EXTERNAL :: _xfGetScalarDatasetsInfo_f
!  MS FORTRAN needs explicit interface for C functions called here.

  INTERFACE
    INTEGER FUNCTION xfgetscalardatasetsinfo_f(a_Id, a_Number, a_MaxPathLength)
	  !MS$ATTRIBUTES C,reference::xfgetscalardatasetsinfo_f
	  INTEGER,            INTENT(IN)     :: a_Id
    INTEGER                            :: a_Number
    INTEGER,            INTENT(OUT)    :: a_MaxPathLength


    END FUNCTION xfgetscalardatasetsinfo_f
  END INTERFACE


  error = xfgetscalardatasetsinfo_f (a_Id, a_Number, a_MaxPathLength)


return

END SUBROUTINE

!***********************

SUBROUTINE XF_GET_SCALAR_DATASET_PATHS (a_Id, a_Number, a_MaxPathLength, a_Paths, error)
INTEGER,          INTENT(IN)   :: a_Id
INTEGER,          INTENT(IN)   :: a_Number
INTEGER,          INTENT(IN)   :: a_MaxPathLength
CHARACTER,        DIMENSION(*) :: a_Paths
INTEGER,          INTENT(OUT)  :: error
INTEGER                          pathlen


!            INTEGER, EXTERNAL :: _xfGetScalarDatasetPaths_f
!  MS FORTRAN needs explicit interface for C functions called here.

  INTERFACE
    INTEGER FUNCTION xfgetscalardatasetpaths_f(a_Id, a_Number, a_MaxPathLength, &
                                               a_Paths, pathlen)
	!MS$ATTRIBUTES C,reference::xfgetscalardatasetpaths_f
    !DEC$ATTRIBUTES reference :: a_Paths
	INTEGER,            INTENT(IN)       :: a_Id
    INTEGER,            INTENT(IN)       :: a_Number, a_MaxPathLength
    CHARACTER,          DIMENSION(*)     ::  a_Paths
    INTEGER,            INTENT(IN)       :: pathlen
    


    END FUNCTION xfgetscalardatasetpaths_f
  END INTERFACE


  pathlen = LEN(a_Paths)
  error = xfgetscalardatasetpaths_f (a_Id, a_Number, a_MaxPathLength, &
                                     a_Paths, pathlen)


return

END SUBROUTINE

!***********************

SUBROUTINE XF_GET_VECTOR_DATASETS_INFO (a_Id, a_Number, a_MaxPathLength, error)
INTEGER             , INTENT(IN)    :: a_Id
INTEGER                       :: a_Number
INTEGER,        INTENT(OUT)   :: a_MaxPathLength
INTEGER,        INTENT(OUT)   :: error

!            INTEGER, EXTERNAL :: _xfGetVectorDatasetsInfo_f
!  MS FORTRAN needs explicit interface for C functions called here.

  INTERFACE
    INTEGER FUNCTION xfgetvectordatasetsinfo_f(a_Id, a_Number, a_MaxPathLength)

	  !MS$ATTRIBUTES C,reference::xfgetvectordatasetsinfo_f
	  INTEGER,            INTENT(IN)     :: a_Id
    INTEGER                            :: a_Number
    INTEGER,            INTENT(OUT)    :: a_MaxPathLength


    END FUNCTION xfgetvectordatasetsinfo_f
  END INTERFACE


  error = xfgetvectordatasetsinfo_f (a_Id, a_Number, a_MaxPathLength)

return

END SUBROUTINE

!***********************

SUBROUTINE XF_GET_VECTOR_DATASET_PATHS (a_Id, a_Number, a_MaxPathLength, a_Paths, error)
INTEGER,          INTENT(IN)    :: a_Id
INTEGER,          INTENT(IN)    :: a_Number
INTEGER,          INTENT(IN)    :: a_MaxPathLength
CHARACTER,        DIMENSION(*)  :: a_Paths
INTEGER,          INTENT(OUT)   :: error
INTEGER                            pathlen

!            INTEGER, EXTERNAL :: _xfGetVectorDatasetPaths_f
!  MS FORTRAN needs explicit interface for C functions called here.

  INTERFACE
    INTEGER FUNCTION xfgetvectordatasetpaths_f(a_Id, a_Number, a_MaxPathLength, &
                                               a_Paths, pathlen)

	!MS$ATTRIBUTES C,reference::xfgetvectordatasetpaths_f
    !DEC$ATTRIBUTES reference :: a_Paths
	INTEGER,            INTENT(IN)     :: a_Id
    INTEGER                            :: a_Number
    INTEGER,            INTENT(IN)     :: a_MaxPathLength
    CHARACTER,          DIMENSION(*)   :: a_Paths
    INTEGER,            INTENT(IN)     :: pathlen



    END FUNCTION xfgetvectordatasetpaths_f
  END INTERFACE

  pathlen = LEN(a_Paths)
  error = xfgetvectordatasetpaths_f (a_Id, a_Number, a_MaxPathLength, &
                                     a_Paths, pathlen)

return

END SUBROUTINE

!***********************

SUBROUTINE XF_GET_DATASET_REFTIME (a_Id, a_dReftime, error)
INTEGER,              INTENT(IN)    :: a_Id
REAL*8,               INTENT(OUT)   :: a_dReftime
INTEGER,              INTENT(OUT)   :: error

!            INTEGER, EXTERNAL :: _xfReadDatasetReftime_f
!  MS FORTRAN needs explicit interface for C functions called here.

  INTERFACE
    INTEGER FUNCTION xfreaddatasetreftime_f(a_Id, a_dReftime)

	  !MS$ATTRIBUTES C,reference::xfreaddatasetreftime_f
	  INTEGER,              INTENT(IN)   :: a_Id
      REAL*8,               INTENT(OUT)  :: a_dReftime



    END FUNCTION xfreaddatasetreftime_f
  END INTERFACE

  error = xfreaddatasetreftime_f (a_Id, a_dReftime)
return

END SUBROUTINE

!**************************

SUBROUTINE XF_SET_DATASET_NUM_TIMES (a_Id, a_NumTimes, error)
INTEGER, INTENT(IN) :: a_Id
INTEGER, INTENT(IN)        :: a_NumTimes
INTEGER, INTENT(OUT)       :: error

!            INTEGER, EXTERNAL :: _xfSetDatasetNumTimes_f
!  MS FORTRAN needs explicit interface for C functions called here.
!
  INTERFACE
    INTEGER FUNCTION xfsetdatasetnumtimes_f (a_Id, a_NumTimes)
	  !MS$ATTRIBUTES C,reference::xfsetdatasetnumtimes_f
    INTEGER,          INTENT(IN)      :: a_Id
    INTEGER,          INTENT(IN)      :: a_NumTimes

    END FUNCTION xfsetdatasetnumtimes_f
  END INTERFACE

  error = xfsetdatasetnumtimes_f (a_Id, a_NumTimes)

return

END SUBROUTINE

! ******

SUBROUTINE XF_CHANGE_SCALAR_VALUES_TIMESTEP_FLOAT(FileId, TimestepIndex, NumValsToEdit, &
                                                  Indices, NewValues, error)
INTEGER, INTENT(OUT) ::   FileId
INTEGER, INTENT(IN) :: TimestepIndex
INTEGER, INTENT(IN)       ::   NumValsToEdit
INTEGER, DIMENSION(*), INTENT(IN) ::   Indices
REAL*4, DIMENSION(*), INTENT(IN) :: NewValues
INTEGER, INTENT(OUT)      ::   error

!            INTEGER, EXTERNAL :: _xfChangeScalarValuesTimestepFloat_f
!  MS FORTRAN needs explicit interface for C functions called here.

  INTERFACE
    INTEGER FUNCTION xfchangescalarvaluestimestepfloat_f(FileId, TimestepIndex, NumValsToEdit, &
                                                  Indices, NewValues)
	  !MS$ATTRIBUTES C,reference::xfchangescalarvaluestimestepfloat_f
    INTEGER, INTENT(OUT) ::   FileId
    INTEGER, INTENT(IN) :: TimestepIndex
    INTEGER, INTENT(IN)       ::   NumValsToEdit
    INTEGER, DIMENSION(*), INTENT(IN) ::   Indices
    REAL*4, DIMENSION(*), INTENT(IN) :: NewValues

    END FUNCTION xfchangescalarvaluestimestepfloat_f
  END INTERFACE

  error = xfchangescalarvaluestimestepfloat_f (FileId, TimestepIndex, NumValsToEdit, &
                                               Indices, NewValues)

return

END SUBROUTINE


!***********************

SUBROUTINE XF_GET_DATASET_NUM_TIMES (a_Id, a_NumTimes, error)
INTEGER,  INTENT(IN)   :: a_Id
INTEGER,  INTENT(OUT)  :: a_NumTimes, error


!            INTEGER, EXTERNAL :: _xfGetDatasetNumTimes_f
!  MS FORTRAN needs explicit interface for C functions called here.

  INTERFACE
    INTEGER FUNCTION xfgetdatasetnumtimes_f(a_Id, a_Numtimes)

	  !MS$ATTRIBUTES C,reference::xfgetdatasetnumtimes_f
	  INTEGER,            INTENT(IN)     :: a_Id
      INTEGER,            INTENT(OUT)    :: a_NumTimes


    END FUNCTION xfgetdatasetnumtimes_f
  END INTERFACE

  error = xfgetdatasetnumtimes_f (a_Id, a_Numtimes)

return

END SUBROUTINE

!***************************

SUBROUTINE XF_GET_DATASET_NUMVALS (a_Id, a_Numvals, error)
INTEGER,        INTENT(IN)    :: a_Id
INTEGER,        INTENT(OUT)   :: a_Numvals
INTEGER,        INTENT(OUT)   :: error


!            INTEGER, EXTERNAL :: _xfGetDatasetNumVals_f
!  MS FORTRAN needs explicit interface for C functions called here.

  INTERFACE
    INTEGER FUNCTION xfgetdatasetnumvals_f(a_Id, a_Numvals)

	  !MS$ATTRIBUTES C,reference::xfgetdatasetnumvals_f
	  INTEGER,            INTENT(IN)     :: a_Id
      INTEGER,            INTENT(OUT)    :: a_Numvals


    END FUNCTION xfgetdatasetnumvals_f
  END INTERFACE

  error = xfgetdatasetnumvals_f (a_Id, a_Numvals)


return

END SUBROUTINE

!***************************

SUBROUTINE XF_GET_DATASET_NUMACTIVE (a_Id, a_NumActivevals, error)
INTEGER             , INTENT(IN)    :: a_Id
INTEGER,        INTENT(OUT)   :: a_NumActiveVals
INTEGER,        INTENT(OUT)   :: error


!            INTEGER, EXTERNAL :: _xfGetDatasetNumActive_f
!  MS FORTRAN needs explicit interface for C functions called here.

  INTERFACE
    INTEGER FUNCTION xfgetdatasetnumactive_f(a_Id, a_NumActiveVals)

	  !MS$ATTRIBUTES C,reference::xfgetdatasetnumactive_f
	  INTEGER,            INTENT(IN)     :: a_Id
    INTEGER,            INTENT(OUT)    :: a_NumActiveVals


    END FUNCTION xfgetdatasetnumactive_f
  END INTERFACE

  error = xfgetdatasetnumactive_f (a_Id, a_NumActiveVals)


return

END SUBROUTINE

!******************************

SUBROUTINE XF_GET_DATASET_NUMCOMPONENTS (a_Id, a_NumComponents, error)
INTEGER             , INTENT(IN)    :: a_Id
INTEGER,        INTENT(OUT)   :: a_NumComponents, error


!            INTEGER, EXTERNAL :: _xfGetDatasetVecNumComponents_f
!  MS FORTRAN needs explicit interface for C functions called here.

  INTERFACE
    INTEGER FUNCTION xfgetdatasetvecnumcomponents_f(a_Id, a_NumComponents)

	  !MS$ATTRIBUTES C,reference::xfgetdatasetvecnumcomponents_f
	  INTEGER,            INTENT(IN)     :: a_Id
    INTEGER,            INTENT(OUT)    :: a_NumComponents


    END FUNCTION xfgetdatasetvecnumcomponents_f
  END INTERFACE

  error = xfgetdatasetvecnumcomponents_f (a_Id, a_NumComponents)

return

END SUBROUTINE

!************************

SUBROUTINE XF_GET_DATASET_UNITS (a_Id, Units, error)
INTEGER             ,   INTENT(IN)  :: a_Id
CHARACTER(LEN=*), INTENT(OUT) :: Units
INTEGER,          INTENT(OUT) :: error
INTEGER                          unitlen


!            INTEGER, EXTERNAL :: _xfGetDatasetUnits_f
!  MS FORTRAN needs explicit interface for C functions called here.

  INTERFACE
    INTEGER FUNCTION xfgetdatasetunits_f(a_Id, Units, unitlen)

	  !MS$ATTRIBUTES C,reference::xfgetdatasetunits_f
    !DEC$ATTRIBUTES reference :: Units
	  INTEGER,            INTENT(IN)     :: a_Id
    CHARACTER(LEN=*),   INTENT(OUT)    :: Units
    INTEGER,            INTENT(IN)     :: unitlen



    END FUNCTION xfgetdatasetunits_f
  END INTERFACE

  error = xfgetdatasetunits_f (a_Id, Units, unitlen)

return

END SUBROUTINE

!************************

SUBROUTINE XF_GET_DATASET_TIME_UNITS (a_Id, Units, error)
INTEGER             ,   INTENT(IN)  :: a_Id
CHARACTER(LEN=*), INTENT(OUT) :: Units
INTEGER,          INTENT(OUT) :: error
INTEGER                          unitlen


!            INTEGER, EXTERNAL :: _xfGetDatasetTimeUnits_f
!  MS FORTRAN needs explicit interface for C functions called here.

  INTERFACE
    INTEGER FUNCTION xfgetdatasettimeunits_f(a_Id, Units, unitlen)

	  !MS$ATTRIBUTES C,reference::xfgetdatasettimeunits_f
    !DEC$ATTRIBUTES reference :: Units
	  INTEGER,            INTENT(IN)     :: a_Id
    CHARACTER(LEN=*),   INTENT(OUT)    :: Units
    INTEGER,            INTENT(IN)     :: unitlen



    END FUNCTION xfgetdatasettimeunits_f
  END INTERFACE

  error = xfgetdatasettimeunits_f (a_Id, Units, unitlen)

return

END SUBROUTINE

!*************************

SUBROUTINE XF_GET_DATASET_TIMES (a_Id, a_nTimes, a_Times, error)
INTEGER             , INTENT(IN) :: a_Id
INTEGER, INTENT(IN)        :: a_nTimes
REAL*8, DIMENSION(*), INTENT(OUT)  :: a_Times
INTEGER, INTENT(OUT)       :: error


!            INTEGER, EXTERNAL :: _xfGetDatasetTimes_f
!  MS FORTRAN needs explicit interface for C functions called here.

  INTERFACE
    INTEGER FUNCTION xfgetdatasettimes_f(a_Id, a_nTimes, a_Times)
    
	  !MS$ATTRIBUTES C,reference::xfgetdatasettimes_f
	  INTEGER,              INTENT(IN)   :: a_Id
    INTEGER,              INTENT(IN)   :: a_nTimes
    REAL*8, DIMENSION(*), INTENT(OUT)  :: a_Times


    END FUNCTION xfgetdatasettimes_f
  END INTERFACE


  error = xfgetdatasettimes_f (a_Id, a_nTimes, a_Times)



return

END SUBROUTINE

!****************************

SUBROUTINE XF_GET_DATASET_MINS (a_Id, a_nTimes, a_Mins, error)
INTEGER             ,     INTENT(IN)    :: a_Id
INTEGER,            INTENT(IN)    ::  a_nTimes
REAL*4, DIMENSION(*), INTENT(OUT)   :: a_Mins
INTEGER,            INTENT(OUT)   ::  error


!            INTEGER, EXTERNAL :: _xfGetDatasetMins_f
!  MS FORTRAN needs explicit interface for C functions called here.

  INTERFACE
    INTEGER FUNCTION xfgetdatasetmins_f(a_Id, a_nTimes, a_Mins)
    
	  !MS$ATTRIBUTES C,reference::xfgetdatasetmins_f
	  INTEGER,            INTENT(IN)   :: a_Id
    INTEGER,            INTENT(IN)   :: a_nTimes
    REAL*4, DIMENSION(*), INTENT(OUT)  :: a_Mins


    END FUNCTION xfgetdatasetmins_f
  END INTERFACE


  error = xfgetdatasetmins_f (a_Id, a_nTimes, a_Mins)


return

END SUBROUTINE

!*****************************

SUBROUTINE XF_GET_DATASET_MAXS (a_Id, a_nTimes, a_Maxs, error)
INTEGER             ,     INTENT(IN)    :: a_Id
INTEGER,            INTENT(IN)    :: a_nTimes
REAL*4, DIMENSION(*), INTENT(OUT)   :: a_Maxs
INTEGER,            INTENT(OUT)   :: error

!            INTEGER, EXTERNAL :: _xfGetDatasetMaxs_f
!  MS FORTRAN needs explicit interface for C functions called here.

  INTERFACE
    INTEGER FUNCTION xfgetdatasetmaxs_f(a_Id, a_nTimes, a_Maxs)
    
	  !MS$ATTRIBUTES C,reference::xfgetdatasetmaxs_f
	  INTEGER,            INTENT(IN)   :: a_Id
    INTEGER,            INTENT(IN)   :: a_nTimes
    REAL*4, DIMENSION(*), INTENT(OUT)  :: a_Maxs


    END FUNCTION xfgetdatasetmaxs_f
  END INTERFACE


  error = xfgetdatasetmaxs_f (a_Id, a_nTimes, a_Maxs)


return

END SUBROUTINE

!******************************

SUBROUTINE XF_READ_ACTIVITY_TIMESTEP(a_Id, a_TimestepIndex, a_NumVals,  &
                                     a_Activity, error)
INTEGER,               INTENT(IN)  :: a_Id
INTEGER,               INTENT(IN)  :: a_TimestepIndex, a_NumVals
INTEGER, DIMENSION(*), INTENT(OUT) :: a_Activity
INTEGER,               INTENT(OUT) :: error


!            INTEGER, EXTERNAL :: _xfReadActivityTimestep_f
!  MS FORTRAN needs explicit interface for C functions called here.

  INTERFACE
    INTEGER FUNCTION xfreadactivitytimestep_f(a_Id, a_TimestepIndex, a_NumVals, &
                                              a_Activity)
    
	  !MS$ATTRIBUTES C,reference::xfreadactivitytimestep_f
	  INTEGER,                 INTENT(IN)   :: a_Id
      INTEGER,               INTENT(IN)   :: a_TimestepIndex, a_NumVals
      INTEGER, DIMENSION(*), INTENT(OUT) :: a_Activity


    END FUNCTION xfreadactivitytimestep_f
  END INTERFACE


  error = xfreadactivitytimestep_f (a_Id, a_TimestepIndex, a_NumVals, a_Activity)

return
END SUBROUTINE

! I don't remember why we have this function
! Remove it until we need it
!SUBROUTINE XF_GET_DATASET_2D_VALUES_FLOAT (a_Id, a_Values, error)

!INTEGER        status, error
!INTEGER              DatasetId, a_Id
!INTEGER        Scalar, Vector
!REAL*4           a_Values(:)
!
!Scalar=XFI_IS_GROUP_OF_TYPE (a_Id, GROUP_TYPE_DATASET_SCALAR, error)
!Vector=XFI_IS_GROUP_OF_TYPE (a_Id, GROUP_TYPE_DATASET_VECTOR, error)
!if ((Scalar .NE. 1).AND.(Vector .NE. 1)) then
!  error = ERROR_DATASET_INVALID
!  return
!endif
!
!call H5Dopen_f (a_Id, DATASET_DSET_VALUES, DatasetId, error)
!if (DatasetId .LT. 0) then
!  error = DatasetId
!  return
!endif
!
!call H5dread_real_2 (DatasetId, H5T_NATIVE_REAL, a_Values, 2, status)
!
!error = status
!return
!
!END SUBROUTINE

!******************************

SUBROUTINE XF_READ_ACTIVE_VALS_AT_INDEX (a_Id, a_Index, a_FirstTime, a_NumTimes, a_Values, error)
INTEGER             ,        INTENT(IN)        :: a_Id
INTEGER,               INTENT(IN)        :: a_Index
INTEGER,               INTENT(IN)        :: a_FirstTime
INTEGER,               INTENT(IN)        :: a_NumTimes
INTEGER, DIMENSION(*), INTENT(OUT)       :: a_Values
INTEGER,               INTENT(OUT)       :: error


!            INTEGER, EXTERNAL :: _xfReadActivityValuesAtIndex_f
!  MS FORTRAN needs explicit interface for C functions called here.

  INTERFACE
    INTEGER FUNCTION xfreadactivityvaluesatindex_f(a_Id, a_Index, a_FirstTime, &
                                                   a_NumTimes, a_Values)
    
	  !MS$ATTRIBUTES C,reference::xfreadactivityvaluesatindex_f
	  INTEGER,               INTENT(IN)   :: a_Id
    INTEGER,               INTENT(IN)   :: a_Index, a_FirstTime, a_NumTimes
    INTEGER, DIMENSION(*), INTENT(OUT)  :: a_Values


    END FUNCTION xfreadactivityvaluesatindex_f
  END INTERFACE


  error = xfreadactivityvaluesatindex_f (a_Id, a_Index, a_FirstTime, &
                                         a_NumTimes, a_Values)


return

END SUBROUTINE

!******************************

SUBROUTINE XF_READ_SCALAR_VALUES_TIMESTEP(a_Id, a_TimestepIndex, a_NumVals, a_Values, &
                                          error)
INTEGER             ,        INTENT(IN)         :: a_Id
INTEGER,               INTENT(IN)         :: a_TimestepIndex, a_NumVals
REAL*4,    DIMENSION(*), INTENT(OUT)        :: a_Values
INTEGER,               INTENT(OUT)        :: error


!            INTEGER, EXTERNAL :: _xfReadScalarValuesTimestep_f
!  MS FORTRAN needs explicit interface for C functions called here.

  INTERFACE
    INTEGER FUNCTION xfreadscalarvaluestimestep_f(a_Id, a_TimestepIndex, a_NumVals, &
                                                  a_Values)
    
	  !MS$ATTRIBUTES C,reference::xfreadscalarvaluestimestep_f
	  INTEGER,            INTENT(IN)   :: a_Id
    INTEGER,            INTENT(IN)   :: a_TimestepIndex, a_NumVals
    REAL*4, DIMENSION(*), INTENT(OUT)  :: a_Values


    END FUNCTION xfreadscalarvaluestimestep_f
  END INTERFACE


  error = xfreadscalarvaluestimestep_f (a_Id, a_TimestepIndex, a_NumVals, &
                                        a_Values)

return
END SUBROUTINE

!******************************

SUBROUTINE XF_READ_SCALAR_VALUES_AT_INDEX (a_Id, a_Index, a_FirstTime, &
                                           a_NumTimes, a_Values, error)
INTEGER             , INTENT(IN) :: a_Id
INTEGER,            INTENT(IN)        :: a_Index
INTEGER,            INTENT(IN)        :: a_FirstTime
INTEGER,            INTENT(IN)        :: a_NumTimes
REAL*4, DIMENSION(*), INTENT(OUT)       :: a_Values
INTEGER,            INTENT(OUT)       :: error


!            INTEGER, EXTERNAL :: _xfReadScalarValuesAtIndex_f
!  MS FORTRAN needs explicit interface for C functions called here.

  INTERFACE
    INTEGER FUNCTION xfreadscalarvaluesatindex_f(a_Id, a_Index, a_FirstTime, &
                                                 a_NumTimes, a_Values)
    
	  !MS$ATTRIBUTES C,reference::xfreadscalarvaluesatindex_f
	  INTEGER,            INTENT(IN)   :: a_Id
    INTEGER,            INTENT(IN)   :: a_Index, a_FirstTime, a_NumTimes
    REAL*4, DIMENSION(*), INTENT(OUT)  :: a_Values


    END FUNCTION xfreadscalarvaluesatindex_f
  END INTERFACE


  error = xfreadscalarvaluesatindex_f (a_Id, a_Index, a_FirstTime, &
                                       a_NumTimes, a_Values)

return

END SUBROUTINE

SUBROUTINE XF_READ_SCALAR_VALUES_AT_INDICES_FLOAT (a_Id, a_nIndices, a_Indices, a_FirstTime, &
                                           a_NumTimes, a_Values, error)
INTEGER             , INTENT(IN) :: a_Id
INTEGER,           INTENT(IN)        :: a_nIndices
INTEGER, DIMENSION(*), INTENT(IN)    :: a_Indices
INTEGER,            INTENT(IN)        :: a_FirstTime
INTEGER,            INTENT(IN)        :: a_NumTimes
REAL*4, DIMENSION(*), INTENT(OUT)       :: a_Values
INTEGER,            INTENT(OUT)       :: error


!            INTEGER, EXTERNAL :: _xfReadScalarValuesAtIndex_f
!  MS FORTRAN needs explicit interface for C functions called here.

  INTERFACE
    INTEGER FUNCTION xfreadscalarvaluesatindicesfloat_f(a_Id, a_nIndices, a_Indices, a_FirstTime, &
                                                 a_NumTimes, a_Values)
    
	  !MS$ATTRIBUTES C,reference::xfreadscalarvaluesatindicesfloat_f
    INTEGER             , INTENT(IN) :: a_Id
    INTEGER,           INTENT(IN)        :: a_nIndices
    INTEGER, DIMENSION(*), INTENT(IN)    :: a_Indices
    INTEGER,            INTENT(IN)        :: a_FirstTime
    INTEGER,            INTENT(IN)        :: a_NumTimes
    REAL*4, DIMENSION(*), INTENT(OUT)       :: a_Values

    END FUNCTION xfreadscalarvaluesatindicesfloat_f
  END INTERFACE


  error = xfreadscalarvaluesatindicesfloat_f (a_Id, a_nIndices, a_Indices, a_FirstTime, &
                                       a_NumTimes, a_Values)

return

END SUBROUTINE

!******************************

SUBROUTINE XF_READ_VECTOR_VALUES_TIMESTEP(a_Id, a_TimestepIndex, a_NumVals, &
                                          a_NumComponents, a_Values, error)
INTEGER             , INTENT(IN)  :: a_Id
INTEGER, INTENT(IN)         :: a_TimestepIndex, a_NumVals, a_NumComponents
REAL*4, INTENT(OUT), DIMENSION(a_NumComponents, a_NumVals) :: a_Values
INTEGER, INTENT(OUT)        :: error


!            INTEGER, EXTERNAL :: _xfReadVectorValuesTimestep_f
!  MS FORTRAN needs explicit interface for C functions called here.

  INTERFACE
    INTEGER FUNCTION xfreadvectorvaluestimestep_f(a_Id, a_TimestepIndex, a_NumVals, &
                                                  a_NumComponents, a_Values)
    
	  !MS$ATTRIBUTES C,reference::xfreadvectorvaluestimestep_f
	  INTEGER,            INTENT(IN)   :: a_Id
    INTEGER,            INTENT(IN)   :: a_TimestepIndex, a_NumVals, a_NumComponents
    REAL*4, INTENT(OUT), DIMENSION(a_NumComponents, a_NumVals) :: a_Values


    END FUNCTION xfreadvectorvaluestimestep_f
  END INTERFACE


  error = xfreadvectorvaluestimestep_f (a_Id, a_TimestepIndex, a_NumVals, &
                                        a_NumComponents, a_Values)



return
END SUBROUTINE

!******************************

SUBROUTINE XF_READ_VECTOR_VALUES_AT_INDEX (a_Id, a_Index, a_FirstTime, &
                                           a_NumTimes, a_NumComponents, a_Values, error)
INTEGER             ,     INTENT(IN)        :: a_Id
INTEGER,            INTENT(IN)        :: a_Index
INTEGER,            INTENT(IN)        :: a_FirstTime
INTEGER,            INTENT(IN)        :: a_NumTimes
INTEGER,            INTENT(IN)        :: a_NumComponents
REAL*4, DIMENSION(*), INTENT(OUT)       :: a_Values
INTEGER,            INTENT(OUT)       :: error


!            INTEGER, EXTERNAL :: _xfReadVectorValuesAtIndex_f
!  MS FORTRAN needs explicit interface for C functions called here.

  INTERFACE
    INTEGER FUNCTION xfreadvectorvaluesatindex_f(a_Id, a_Index, a_FirstTime, a_NumTimes, &
                                                 a_NumComponents, a_Values)
    
	  !MS$ATTRIBUTES C,reference::xfreadvectorvaluesatindex_f
	  INTEGER,            INTENT(IN)   :: a_Id
    INTEGER,            INTENT(IN)   :: a_Index, a_FirstTime, a_NumTimes, a_NumComponents
    REAL*4, DIMENSION(*), INTENT(OUT)  :: a_Values


    END FUNCTION xfreadvectorvaluesatindex_f
  END INTERFACE


  error = xfreadvectorvaluesatindex_f (a_Id, a_Index, a_FirstTime, a_NumTimes, &
                                       a_NumComponents, a_Values)

return

END SUBROUTINE

!******************************

SUBROUTINE XF_SCALAR_DATA_LOCATION (a_Id, a_DataLoc, error)
INTEGER             , INTENT(IN) :: a_Id
INTEGER, INTENT(IN)        :: a_DataLoc
INTEGER, INTENT(OUT)       :: error


!            INTEGER, EXTERNAL :: _xfScalarDataLocation_f
!  MS FORTRAN needs explicit interface for C functions called here.

  INTERFACE
    INTEGER FUNCTION xfscalardatalocation_f(a_Id, a_DataLoc)
    
	  !MS$ATTRIBUTES C,reference::xfscalardatalocation_f
	  INTEGER,            INTENT(IN)   :: a_Id
    INTEGER,            INTENT(IN)   :: a_DataLoc


    END FUNCTION xfscalardatalocation_f
  END INTERFACE


  error = xfscalardatalocation_f (a_Id, a_DataLoc)


return

END SUBROUTINE

!******************************

SUBROUTINE XF_VECTOR_2D_DATA_LOCS (a_Id, a_DataLocI, a_DataLocJ, error)
INTEGER             , INTENT(IN) :: a_Id
INTEGER, INTENT(IN)        :: a_DataLocI, a_DataLocJ
INTEGER, INTENT(OUT)       :: error


!            INTEGER, EXTERNAL :: _xfVector2DDataLocations_f
!  MS FORTRAN needs explicit interface for C functions called here.

  INTERFACE
    INTEGER FUNCTION xfvector2ddatalocations_f(a_Id, a_DataLocI, a_DataLocJ)
    
	  !MS$ATTRIBUTES C,reference::xfvector2ddatalocations_f
	  INTEGER,            INTENT(IN)   :: a_Id
    INTEGER,            INTENT(IN)   :: a_DataLocI, a_DataLocJ


    END FUNCTION xfvector2ddatalocations_f
  END INTERFACE


  error = xfvector2ddatalocations_f (a_Id, a_DataLocI, a_DataLocJ)


return

END SUBROUTINE

!******************************

SUBROUTINE XF_VECTOR_3D_DATA_LOCS (a_Id, a_DataLocI, a_DataLocJ, a_DataLock, error)
INTEGER             , INTENT(IN) :: a_Id
INTEGER, INTENT(IN)        :: a_DataLocI, a_DataLocJ, a_DataLock
INTEGER, INTENT(OUT)       :: error


!            INTEGER, EXTERNAL :: _xfVector3DDataLocations_f
!  MS FORTRAN needs explicit interface for C functions called here.

  INTERFACE
    INTEGER FUNCTION xfvector3ddatalocations_f(a_Id, a_DataLocI, a_DataLocJ, a_DataLock)
    
	  !MS$ATTRIBUTES C,reference::xfvector3ddatalocations_f
	  INTEGER,            INTENT(IN)   :: a_Id
    INTEGER,            INTENT(IN)   :: a_DataLocI, a_DataLocJ, a_DataLock


    END FUNCTION xfvector3ddatalocations_f
  END INTERFACE


  error = xfvector3ddatalocations_f (a_Id, a_DataLocI, a_DataLocJ, a_DataLocK)


return

END SUBROUTINE

!******************************

SUBROUTINE XF_GET_SCALAR_DATA_LOC (a_Id, a_DataLoc, error)
INTEGER             , INTENT(IN) :: a_Id
INTEGER, INTENT(IN)        :: a_DataLoc
INTEGER, INTENT(OUT)       :: error


!            INTEGER, EXTERNAL :: _xfGetScalarDataLocation_f
!  MS FORTRAN needs explicit interface for C functions called here.

  INTERFACE
    INTEGER FUNCTION xfgetscalardatalocation_f(a_Id, a_DataLoc)
    
	  !MS$ATTRIBUTES C,reference::xfgetscalardatalocation_f
	  INTEGER,            INTENT(IN)   :: a_Id
    INTEGER,            INTENT(IN)   :: a_DataLoc


    END FUNCTION xfgetscalardatalocation_f
  END INTERFACE


  error = xfgetscalardatalocation_f (a_Id, a_DataLoc)


return

END SUBROUTINE

!******************************

SUBROUTINE XF_GET_VECTOR_2D_DATA_LOCS (a_Id, a_DataLocI, a_DataLocJ, error)
INTEGER             , INTENT(IN) :: a_Id
INTEGER, INTENT(IN)        :: a_DataLocI, a_DataLocJ
INTEGER, INTENT(OUT)       :: error


!            INTEGER, EXTERNAL :: _xfGetVector2DDataLocations_f
!  MS FORTRAN needs explicit interface for C functions called here.

  INTERFACE
    INTEGER FUNCTION xfgetvector2ddatalocations_f(a_Id, a_DataLocI, a_DataLocJ)
    
	  !MS$ATTRIBUTES C,reference::xfgetvector2ddatalocations_f
	  INTEGER,            INTENT(IN)   :: a_Id
    INTEGER,            INTENT(IN)   :: a_DataLocI, a_DataLocJ


    END FUNCTION xfgetvector2ddatalocations_f
  END INTERFACE


  error = xfgetvector2ddatalocations_f (a_Id, a_DataLocI, a_DataLocJ)


return

END SUBROUTINE

!******************************

SUBROUTINE XF_GET_VECTOR_3D_DATA_LOCS (a_Id, a_DataLocI, a_DataLocJ, a_DataLocK, error)
INTEGER             , INTENT(IN) :: a_Id
INTEGER, INTENT(IN)        :: a_DataLocI, a_DataLocJ, a_DataLocK
INTEGER, INTENT(OUT)       :: error


!            INTEGER, EXTERNAL :: _xfGetVector3DDataLocations_f
!  MS FORTRAN needs explicit interface for C functions called here.

  INTERFACE
    INTEGER FUNCTION xfgetvector3ddatalocations_f(a_Id, a_DataLocI, a_DataLocJ, a_DataLocK)
    
	  !MS$ATTRIBUTES C,reference::xfgetvector3ddatalocations_f
	  INTEGER,            INTENT(IN)   :: a_Id
    INTEGER,            INTENT(IN)   :: a_DataLocI, a_DataLocJ, a_DataLocK


    END FUNCTION xfgetvector3ddatalocations_f
  END INTERFACE


  error = xfgetvector3ddatalocations_f (a_Id, a_DataLocI, a_DataLocJ, a_DataLocK)


return

END SUBROUTINE

!******************************

SUBROUTINE XF_VECTORS_IN_LOCAL_COORDS (a_Id, error)

  INTEGER             , INTENT(IN) :: a_Id
  INTEGER, INTENT(OUT)       :: error

  !            INTEGER, EXTERNAL :: _xfVectorsInLocalCoords_f
!  MS FORTRAN needs explicit interface for C functions called here.

  INTERFACE
    INTEGER FUNCTION xfvectorsinlocalcoords_f(a_Id)
    
	  !MS$ATTRIBUTES C,reference::xfvectorsinlocalcoords_f
	  INTEGER,            INTENT(IN)   :: a_Id


    END FUNCTION xfvectorsinlocalcoords_f
  END INTERFACE


  error = xfvectorsinlocalcoords_f (a_Id)

return

END SUBROUTINE

!******************************

SUBROUTINE XF_ARE_VECTORS_IN_LOCAL_COORDS (a_Id, a_LocalCoords, error)

  INTEGER             , INTENT(IN) :: a_Id
  INTEGER, INTENT(OUT)       :: a_LocalCoords, error


!            INTEGER, EXTERNAL :: _xfAreVectorsInLocalCoords_f
!  MS FORTRAN needs explicit interface for C functions called here.

  INTERFACE
    INTEGER FUNCTION xfarevectorsinlocalcoords_f(a_Id, a_LocalCoords)
    
	  !MS$ATTRIBUTES C,reference::xfarevectorsinlocalcoords_f
	  INTEGER,            INTENT(IN)   :: a_Id
    INTEGER,            INTENT(OUT)  :: a_LocalCoords


    END FUNCTION xfarevectorsinlocalcoords_f
  END INTERFACE


  error = xfarevectorsinlocalcoords_f (a_Id, a_LocalCoords)
 
return
END SUBROUTINE XF_ARE_VECTORS_IN_LOCAL_COORDS

!------------------------------------------------------------------------------
! SUBROUTINE  XF_GET_HORIZ_DATUM
! PURPOSE     Read the horizontal datum (type defines at top of xmdf.f90)
! NOTES
!------------------------------------------------------------------------------
SUBROUTINE XF_GET_HORIZ_DATUM (a_CoordId, a_val, error)
INTEGER             , INTENT(IN) :: a_CoordId
INTEGER, INTENT(OUT)       :: a_val
INTEGER, INTENT(OUT)       :: error

!            INTEGER, EXTERNAL :: _xfGetHorizDatum_f
!  MS FORTRAN needs explicit interface for C functions called here.

  INTERFACE
    INTEGER FUNCTION xfgethorizdatum_f(a_CoordId, a_val)
    
	  !MS$ATTRIBUTES C,reference::xfgethorizdatum_f
	  INTEGER,            INTENT(IN)   :: a_CoordId
    INTEGER,            INTENT(OUT)  :: a_val


    END FUNCTION xfgethorizdatum_f
  END INTERFACE


  error = xfgethorizdatum_f (a_CoordId, a_val)

  return
END SUBROUTINE
!------------------------------------------------------------------------------
! SUBROUTINE  XF_GET_HORIZ_UNITS
! PURPOSE     Read the horizontal units (type defines at top of xmdf.f90)
! NOTES
!------------------------------------------------------------------------------
SUBROUTINE XF_GET_HORIZ_UNITS (a_CoordId, a_val, error)
INTEGER             , INTENT(IN) :: a_CoordId
INTEGER, INTENT(OUT)       :: a_val
INTEGER, INTENT(OUT)       :: error


!            INTEGER, EXTERNAL :: _xfGetHorizUnits_f
!  MS FORTRAN needs explicit interface for C functions called here.

  INTERFACE
    INTEGER FUNCTION xfgethorizunits_f(a_CoordId, a_val)
    
	  !MS$ATTRIBUTES C,reference::xfgethorizunits_f
	  INTEGER,            INTENT(IN)   :: a_CoordId
    INTEGER,            INTENT(OUT)  :: a_val


    END FUNCTION xfgethorizunits_f
  END INTERFACE


  error = xfgethorizunits_f (a_CoordId, a_val)


  return
END SUBROUTINE
!------------------------------------------------------------------------------
! SUBROUTINE  XF_GET_VERT_DATUM
! PURPOSE     Read the vertical datum (type defines at top of xmdf.f90)
! NOTES
!------------------------------------------------------------------------------
SUBROUTINE XF_GET_VERT_DATUM (a_CoordId, a_val, error)
INTEGER             , INTENT(IN) :: a_CoordId
INTEGER, INTENT(OUT)       :: a_val
INTEGER, INTENT(OUT)       :: error


!            INTEGER, EXTERNAL :: _xfGetVertDatum_f
!  MS FORTRAN needs explicit interface for C functions called here.

  INTERFACE
    INTEGER FUNCTION xfgetvertdatum_f(a_CoordId, a_val)
    
	  !MS$ATTRIBUTES C,reference::xfgetvertdatum_f
	  INTEGER,            INTENT(IN)   :: a_CoordId
    INTEGER,            INTENT(OUT)  :: a_val


    END FUNCTION xfgetvertdatum_f
  END INTERFACE


  error = xfgetvertdatum_f (a_CoordId, a_val)


  return
END SUBROUTINE
!------------------------------------------------------------------------------
! SUBROUTINE  XF_GET_VERT_UNITS
! PURPOSE     Read the vertical units (type defines at top of xmdf.f90)
! NOTES
!------------------------------------------------------------------------------
SUBROUTINE XF_GET_VERT_UNITS (a_CoordId, a_val, error)
INTEGER             , INTENT(IN) :: a_CoordId
INTEGER, INTENT(OUT)       :: a_val
INTEGER, INTENT(OUT)       :: error

!            INTEGER, EXTERNAL :: _xfGetVertUnits_f
!  MS FORTRAN needs explicit interface for C functions called here.

  INTERFACE
    INTEGER FUNCTION xfgetvertunits_f(a_CoordId, a_val)
    
	  !MS$ATTRIBUTES C,reference::xfgetvertunits_f
	  INTEGER,            INTENT(IN)   :: a_CoordId
    INTEGER,            INTENT(OUT)  :: a_val


    END FUNCTION xfgetvertunits_f
  END INTERFACE


  error = xfgetvertunits_f (a_CoordId, a_val)

  return
END SUBROUTINE
!------------------------------------------------------------------------------
! SUBROUTINE  XF_GET_LAT
! PURPOSE     Read whether the lattitude is North or South latitude
!              (type defines at top of xmdf.f90)
! NOTES
!------------------------------------------------------------------------------
SUBROUTINE XF_GET_LAT (a_CoordId, a_val, error)
INTEGER             , INTENT(IN) :: a_CoordId
INTEGER, INTENT(OUT)       :: a_val
INTEGER, INTENT(OUT)       :: error

!            INTEGER, EXTERNAL :: _xfGetLat_f
!  MS FORTRAN needs explicit interface for C functions called here.

  INTERFACE
    INTEGER FUNCTION xfgetlat_f(a_CoordId, a_val)
    
	  !MS$ATTRIBUTES C,reference::xfgetlat_f
	  INTEGER,            INTENT(IN)   :: a_CoordId
    INTEGER,            INTENT(OUT)  :: a_val


    END FUNCTION xfgetlat_f
  END INTERFACE


  error = xfgetlat_f (a_CoordId, a_val)

  return
END SUBROUTINE
!------------------------------------------------------------------------------
! SUBROUTINE  XF_GET_LON
! PURPOSE     Read whether the longitude is East or West
!              (type defines at top of xmdf.f90)
! NOTES
!------------------------------------------------------------------------------
SUBROUTINE XF_GET_LON (a_CoordId, a_val, error)
INTEGER             , INTENT(IN) :: a_CoordId
INTEGER, INTENT(OUT)       :: a_val
INTEGER, INTENT(OUT)       :: error

!            INTEGER, EXTERNAL :: _xfGetLon_f
!  MS FORTRAN needs explicit interface for C functions called here.

  INTERFACE
    INTEGER FUNCTION xfgetlon_f(a_CoordId, a_val)
    
	  !MS$ATTRIBUTES C,reference::xfgetlon_f
	  INTEGER,            INTENT(IN)   :: a_CoordId
    INTEGER,            INTENT(OUT)  :: a_val


    END FUNCTION xfgetlon_f
  END INTERFACE


  error = xfgetlon_f (a_CoordId, a_val)
  
  return
END SUBROUTINE
!------------------------------------------------------------------------------
! SUBROUTINE  XF_GET_UTM_ZONE
! PURPOSE     Read the UTM zone (should be a number between 1 and 60)
! NOTES
!------------------------------------------------------------------------------
SUBROUTINE XF_GET_UTM_ZONE (a_CoordId, a_val, error)
INTEGER             , INTENT(IN) :: a_CoordId
INTEGER, INTENT(OUT)       :: a_val
INTEGER, INTENT(OUT)       :: error

!            INTEGER, EXTERNAL :: _xfGetUTMZone_f
!  MS FORTRAN needs explicit interface for C functions called here.

  INTERFACE
    INTEGER FUNCTION xfgetutmzone_f(a_CoordId, a_val)
    
	  !MS$ATTRIBUTES C,reference::xfgetutmzone_f
	  INTEGER,            INTENT(IN)   :: a_CoordId
    INTEGER,            INTENT(OUT)  :: a_val


    END FUNCTION xfgetutmzone_f
  END INTERFACE


  error = xfgetutmzone_f (a_CoordId, a_val)


  return
END SUBROUTINE
!------------------------------------------------------------------------------
! SUBROUTINE  XF_GET_SPC_ZONE
! PURPOSE     Read the SPC zone (Lookup numbers in XMDF documentation)
! NOTES
!------------------------------------------------------------------------------
SUBROUTINE XF_GET_SPC_ZONE (a_CoordId, a_val, error)
INTEGER             , INTENT(IN) :: a_CoordId
INTEGER, INTENT(OUT)       :: a_val
INTEGER, INTENT(OUT)       :: error

!            INTEGER, EXTERNAL :: _xfGetSPCZone_f
!  MS FORTRAN needs explicit interface for C functions called here.

  INTERFACE
    INTEGER FUNCTION xfgetspczone_f(a_CoordId, a_val)
    
	  !MS$ATTRIBUTES C,reference::xfgetspczone_f
	  INTEGER,            INTENT(IN)   :: a_CoordId
    INTEGER,            INTENT(OUT)  :: a_val


    END FUNCTION xfgetspczone_f
  END INTERFACE


  error = xfgetspczone_f (a_CoordId, a_val)


  return
END SUBROUTINE
!------------------------------------------------------------------------------
! SUBROUTINE  XF_GET_HPGN_AREA
! PURPOSE     Read the HPGN area (Lookup numbers in XMDF documentation)
! NOTES
!------------------------------------------------------------------------------
SUBROUTINE XF_GET_HPGN_AREA (a_CoordId, a_val, error)
INTEGER, INTENT(IN)        :: a_CoordId
INTEGER, INTENT(OUT)       :: a_val
INTEGER, INTENT(OUT)       :: error


!            INTEGER, EXTERNAL :: _xfGetHPGNArea_f
!  MS FORTRAN needs explicit interface for C functions called here.

  INTERFACE
    INTEGER FUNCTION xfgethpgnarea_f(a_CoordId, a_val)
    
	  !MS$ATTRIBUTES C,reference::xfgethpgnarea_f
	  INTEGER,            INTENT(IN)   :: a_CoordId
    INTEGER,            INTENT(OUT)  :: a_val


    END FUNCTION xfgethpgnarea_f
  END INTERFACE


  error = xfgethpgnarea_f (a_CoordId, a_val)


  return
END SUBROUTINE
!------------------------------------------------------------------------------
! SUBROUTINE  XF_GET_CPP_LAT
! PURPOSE     Read the Carte Prallelo Grammatique Projetion Factor for
!             converting the latitude
! NOTES
!------------------------------------------------------------------------------
SUBROUTINE XF_GET_CPP_LAT (a_CoordId, a_val, error)
INTEGER, INTENT(IN)   :: a_CoordId
REAL*8,  INTENT(OUT)  :: a_val
INTEGER, INTENT(OUT)  :: error


!            INTEGER, EXTERNAL :: _xfGetCPPLat_f
!  MS FORTRAN needs explicit interface for C functions called here.

  INTERFACE
    INTEGER FUNCTION xfgetcpplat_f(a_CoordId, a_val)
    
	  !MS$ATTRIBUTES C,reference::xfgetcpplat_f
	  INTEGER, INTENT(IN)   :: a_CoordId
      REAL*8,  INTENT(OUT)  :: a_val


    END FUNCTION xfgetcpplat_f
  END INTERFACE


  error = xfgetcpplat_f (a_CoordId, a_val)

  return
END SUBROUTINE
!------------------------------------------------------------------------------
! SUBROUTINE  XF_GET_CPP_LON
! PURPOSE     Read the Carte Prallelo Grammatique Projetion Factor for
!             converting the longitude
! NOTES
!------------------------------------------------------------------------------
SUBROUTINE XF_GET_CPP_LON (a_CoordId, a_val, error)
INTEGER, INTENT(IN)   :: a_CoordId
REAL*8,  INTENT(OUT)  :: a_val
INTEGER, INTENT(OUT)  :: error


!            INTEGER, EXTERNAL :: _xfGetCPPLon_f
!  MS FORTRAN needs explicit interface for C functions called here.

  INTERFACE
    INTEGER FUNCTION xfgetcpplon_f(a_CoordId, a_val)
    
	  !MS$ATTRIBUTES C,reference::xfgetcpplon_f
	  INTEGER, INTENT(IN)   :: a_CoordId
      REAL*8,  INTENT(OUT)  :: a_val


    END FUNCTION xfgetcpplon_f
  END INTERFACE


  error = xfgetcpplon_f (a_CoordId, a_val)

  return
END SUBROUTINE
!------------------------------------------------------------------------------
! SUBROUTINE  XF_GET_ELLIPSE
! PURPOSE     Read the Ellipse number based upon XMDF documentation
! NOTES
!------------------------------------------------------------------------------
SUBROUTINE XF_GET_ELLIPSE (a_CoordId, a_val, error)
INTEGER, INTENT(IN)        :: a_CoordId
INTEGER, INTENT(OUT)       :: a_val
INTEGER, INTENT(OUT)       :: error

!            INTEGER, EXTERNAL :: _xfGetEllipse_f
!  MS FORTRAN needs explicit interface for C functions called here.

  INTERFACE
    INTEGER FUNCTION xfgetellipse_f(a_CoordId, a_val)
    
	  !MS$ATTRIBUTES C,reference::xfgetellipse_f
	  INTEGER,            INTENT(IN)   :: a_CoordId
      INTEGER,            INTENT(OUT)  :: a_val


    END FUNCTION xfgetellipse_f
  END INTERFACE


  error = xfgetellipse_f (a_CoordId, a_val)


  return
END SUBROUTINE
!------------------------------------------------------------------------------
! SUBROUTINE  XF_GET_MAJOR_R
! PURPOSE     Read a user-defined ellipse major radius
! NOTES
!------------------------------------------------------------------------------
SUBROUTINE XF_GET_MAJOR_R (a_CoordId, a_val, error)
INTEGER, INTENT(IN)   :: a_CoordId
REAL*8,  INTENT(OUT)  :: a_val
INTEGER, INTENT(OUT)  :: error


!            INTEGER, EXTERNAL :: _xfGetMajorR_f
!  MS FORTRAN needs explicit interface for C functions called here.

  INTERFACE
    INTEGER FUNCTION xfgetmajorr_f(a_CoordId, a_val)
    
	  !MS$ATTRIBUTES C,reference::xfgetmajorr_f
	  INTEGER, INTENT(IN)   :: a_CoordId
      REAL*8,  INTENT(OUT)  :: a_val


    END FUNCTION xfgetmajorr_f
  END INTERFACE


  error = xfgetmajorr_f (a_CoordId, a_val)

  return
END SUBROUTINE
!------------------------------------------------------------------------------
! SUBROUTINE  XF_GET_MINOR_R
! PURPOSE     Read a user-defined ellipse minor radius
! NOTES
!------------------------------------------------------------------------------
SUBROUTINE XF_GET_MINOR_R (a_CoordId, a_val, error)
INTEGER, INTENT(IN)   :: a_CoordId
REAL*8,  INTENT(OUT)  :: a_val
INTEGER, INTENT(OUT)  :: error


!            INTEGER, EXTERNAL :: _xfGetMinorR_f
!  MS FORTRAN needs explicit interface for C functions called here.

  INTERFACE
    INTEGER FUNCTION xfgetminorr_f(a_CoordId, a_val)
    
	  !MS$ATTRIBUTES C,reference::xfgetminorr_f
	  INTEGER, INTENT(IN)   :: a_CoordId
      REAL*8,  INTENT(OUT)  :: a_val


    END FUNCTION xfgetminorr_f
  END INTERFACE


  error = xfgetminorr_f (a_CoordId, a_val)

  return
END SUBROUTINE
!------------------------------------------------------------------------------
! SUBROUTINE  XF_SET_HORIZ_DATUM
! PURPOSE
! NOTES
!------------------------------------------------------------------------------
SUBROUTINE XF_SET_HORIZ_DATUM (a_CoordId, a_val, error)
INTEGER, INTENT(IN)        :: a_CoordId
INTEGER, INTENT(IN)        :: a_val
INTEGER, INTENT(OUT)       :: error

!            INTEGER, EXTERNAL :: _xfSetHorizDatum_f
!  MS FORTRAN needs explicit interface for C functions called here.

  INTERFACE
    INTEGER FUNCTION xfsethorizdatum_f(a_CoordId, a_val)
    
	  !MS$ATTRIBUTES C,reference::xfsethorizdatum_f
	  INTEGER,            INTENT(IN)   :: a_CoordId
      INTEGER,            INTENT(IN)   :: a_val


    END FUNCTION xfsethorizdatum_f
  END INTERFACE


  error = xfsethorizdatum_f (a_CoordId, a_val)


  return
END SUBROUTINE
!------------------------------------------------------------------------------
! SUBROUTINE  XF_SET_HORIZ_UNITS
! PURPOSE
! NOTES
!------------------------------------------------------------------------------
SUBROUTINE XF_SET_HORIZ_UNITS (a_CoordId, a_val, error)
INTEGER, INTENT(IN)        :: a_CoordId
INTEGER, INTENT(IN)        :: a_val
INTEGER, INTENT(OUT)       :: error

!            INTEGER, EXTERNAL :: _xfSetHorizUnits_f
!  MS FORTRAN needs explicit interface for C functions called here.

  INTERFACE
    INTEGER FUNCTION xfsethorizunits_f(a_CoordId, a_val)
    
	  !MS$ATTRIBUTES C,reference::xfsethorizunits_f
	  INTEGER,            INTENT(IN)   :: a_CoordId
      INTEGER,            INTENT(IN)   :: a_val


    END FUNCTION xfsethorizunits_f
  END INTERFACE


  error = xfsethorizunits_f (a_CoordId, a_val)

  return
END SUBROUTINE
!------------------------------------------------------------------------------
! SUBROUTINE  XF_SET_VERT_DATUM
! PURPOSE
! NOTES
!------------------------------------------------------------------------------
SUBROUTINE XF_SET_VERT_DATUM (a_CoordId, a_val, error)
INTEGER, INTENT(IN)        :: a_CoordId
INTEGER, INTENT(IN)        :: a_val
INTEGER, INTENT(OUT)       :: error

!            INTEGER, EXTERNAL :: _xfSetVertDatum_f
!  MS FORTRAN needs explicit interface for C functions called here.

  INTERFACE
    INTEGER FUNCTION xfsetvertdatum_f(a_CoordId, a_val)
    
	  !MS$ATTRIBUTES C,reference::xfsetvertdatum_f
	  INTEGER,            INTENT(IN)   :: a_CoordId
      INTEGER,            INTENT(IN)   :: a_val


    END FUNCTION xfsetvertdatum_f
  END INTERFACE


  error = xfsetvertdatum_f (a_CoordId, a_val)

  return
END SUBROUTINE
!------------------------------------------------------------------------------
! SUBROUTINE  XF_SET_VERT_UNITS
! PURPOSE
! NOTES
!------------------------------------------------------------------------------
SUBROUTINE XF_SET_VERT_UNITS (a_CoordId, a_val, error)
INTEGER             , INTENT(IN) :: a_CoordId
INTEGER, INTENT(IN)        :: a_val
INTEGER, INTENT(OUT)       :: error

!            INTEGER, EXTERNAL :: _xfSetVertUnits_f
!  MS FORTRAN needs explicit interface for C functions called here.

  INTERFACE
    INTEGER FUNCTION xfsetvertunits_f(a_CoordId, a_val)
    
	  !MS$ATTRIBUTES C,reference::xfsetvertunits_f
	  INTEGER,            INTENT(IN)   :: a_CoordId
    INTEGER,            INTENT(IN)   :: a_val


    END FUNCTION xfsetvertunits_f
  END INTERFACE


  error = xfsetvertunits_f (a_CoordId, a_val)


  return
END SUBROUTINE
!------------------------------------------------------------------------------
! SUBROUTINE  XF_SET_LAT
! PURPOSE
! NOTES
!------------------------------------------------------------------------------
SUBROUTINE XF_SET_LAT (a_CoordId, a_val, error)
INTEGER             , INTENT(IN) :: a_CoordId
INTEGER, INTENT(IN)        :: a_val
INTEGER, INTENT(OUT)       :: error

!            INTEGER, EXTERNAL :: _xfSetLat_f
!  MS FORTRAN needs explicit interface for C functions called here.

  INTERFACE
    INTEGER FUNCTION xfsetlat_f(a_CoordId, a_val)
    
	  !MS$ATTRIBUTES C,reference::xfsetlat_f
	  INTEGER,            INTENT(IN)   :: a_CoordId
    INTEGER,            INTENT(IN)   :: a_val


    END FUNCTION xfsetlat_f
  END INTERFACE


  error = xfsetlat_f (a_CoordId, a_val)

  return
END SUBROUTINE
!------------------------------------------------------------------------------
! SUBROUTINE  XF_SET_LON
! PURPOSE
! NOTES
!------------------------------------------------------------------------------
SUBROUTINE XF_SET_LON (a_CoordId, a_val, error)
INTEGER             , INTENT(IN) :: a_CoordId
INTEGER, INTENT(IN)        :: a_val
INTEGER, INTENT(OUT)       :: error

!            INTEGER, EXTERNAL :: _xfSetLon_f
!  MS FORTRAN needs explicit interface for C functions called here.

  INTERFACE
    INTEGER FUNCTION xfsetlon_f(a_CoordId, a_val)
    
	  !MS$ATTRIBUTES C,reference::xfsetlon_f
	  INTEGER,            INTENT(IN)   :: a_CoordId
    INTEGER,            INTENT(IN)   :: a_val


    END FUNCTION xfsetlon_f
  END INTERFACE


  error = xfsetlon_f (a_CoordId, a_val)

  return
END SUBROUTINE
!------------------------------------------------------------------------------
! SUBROUTINE  XF_SET_UTM_ZONE
! PURPOSE
! NOTES
!------------------------------------------------------------------------------
SUBROUTINE XF_SET_UTM_ZONE (a_CoordId, a_val, error)
INTEGER             , INTENT(IN) :: a_CoordId
INTEGER, INTENT(IN)        :: a_val
INTEGER, INTENT(OUT)       :: error

!            INTEGER, EXTERNAL :: _xfSetUTMZone_f
!  MS FORTRAN needs explicit interface for C functions called here.

  INTERFACE
    INTEGER FUNCTION xfsetutmzone_f(a_CoordId, a_val)
    
	  !MS$ATTRIBUTES C,reference::xfsetutmzone_f
	  INTEGER,            INTENT(IN)   :: a_CoordId
    INTEGER,            INTENT(IN)   :: a_val


    END FUNCTION xfsetutmzone_f
  END INTERFACE


  error = xfsetutmzone_f (a_CoordId, a_val)

  return
END SUBROUTINE
!------------------------------------------------------------------------------
! SUBROUTINE  XF_SET_SPC_ZONE
! PURPOSE
! NOTES
!------------------------------------------------------------------------------
SUBROUTINE XF_SET_SPC_ZONE (a_CoordId, a_val, error)
INTEGER             , INTENT(IN) :: a_CoordId
INTEGER, INTENT(IN)        :: a_val
INTEGER, INTENT(OUT)       :: error

!            INTEGER, EXTERNAL :: _xfSetSPCZone_f
!  MS FORTRAN needs explicit interface for C functions called here.

  INTERFACE
    INTEGER FUNCTION xfsetspczone_f(a_CoordId, a_val)
    
	  !MS$ATTRIBUTES C,reference::xfsetspczone_f
	  INTEGER,            INTENT(IN)   :: a_CoordId
    INTEGER,            INTENT(IN)   :: a_val


    END FUNCTION xfsetspczone_f
  END INTERFACE


  error = xfsetspczone_f (a_CoordId, a_val)

  return
END SUBROUTINE
!------------------------------------------------------------------------------
! SUBROUTINE  XF_SET_HPGN_AREA
! PURPOSE
! NOTES
!------------------------------------------------------------------------------
SUBROUTINE XF_SET_HPGN_AREA(a_CoordId, a_val, error)
INTEGER             , INTENT(IN) :: a_CoordId
INTEGER, INTENT(IN)        :: a_val
INTEGER, INTENT(OUT)       :: error

!            INTEGER, EXTERNAL :: _xfSetHPGNArea_f
!  MS FORTRAN needs explicit interface for C functions called here.

  INTERFACE
    INTEGER FUNCTION xfsethpgnarea_f(a_CoordId, a_val)
    
	  !MS$ATTRIBUTES C,reference::xfsethpgnarea_f
	  INTEGER,            INTENT(IN)   :: a_CoordId
    INTEGER,            INTENT(IN)   :: a_val


    END FUNCTION xfsethpgnarea_f
  END INTERFACE


  error = xfsethpgnarea_f (a_CoordId, a_val)

  return
END SUBROUTINE
!------------------------------------------------------------------------------
! SUBROUTINE  XF_SET_CPP_LAT
! PURPOSE
! NOTES
!------------------------------------------------------------------------------
SUBROUTINE XF_SET_CPP_LAT (a_CoordId, a_val, error)
INTEGER, INTENT(IN)   :: a_CoordId
REAL*8,  INTENT(IN)   :: a_val
INTEGER, INTENT(OUT)  :: error


!            INTEGER, EXTERNAL :: _xfSetCPPLat_f
!  MS FORTRAN needs explicit interface for C functions called here.

  INTERFACE
    INTEGER FUNCTION xfsetcpplat_f(a_CoordId, a_val)
    
	  !MS$ATTRIBUTES C,reference::xfsetcpplat_f
	  INTEGER,  INTENT(IN)   :: a_CoordId
      REAL*8,   INTENT(IN)   :: a_val


    END FUNCTION xfsetcpplat_f
  END INTERFACE


  error = xfsetcpplat_f (a_CoordId, a_val)

  return
END SUBROUTINE
!------------------------------------------------------------------------------
! SUBROUTINE  XF_SET_CPP_LON
! PURPOSE
! NOTES
!------------------------------------------------------------------------------
SUBROUTINE XF_SET_CPP_LON (a_CoordId, a_val, error)
INTEGER, INTENT(IN)   :: a_CoordId
REAL*8,  INTENT(IN)   :: a_val
INTEGER, INTENT(OUT)  :: error


!            INTEGER, EXTERNAL :: _xfSetCPPLon_f
!  MS FORTRAN needs explicit interface for C functions called here.

  INTERFACE
    INTEGER FUNCTION xfsetcpplon_f(a_CoordId, a_val)
    
	  !MS$ATTRIBUTES C,reference::xfsetcpplon_f
	  INTEGER, INTENT(IN)   :: a_CoordId
      REAL*8,  INTENT(IN)   :: a_val


    END FUNCTION xfsetcpplon_f
  END INTERFACE


  error = xfsetcpplon_f (a_CoordId, a_val)

  return
END SUBROUTINE
!------------------------------------------------------------------------------
! SUBROUTINE  XF_SET_ELLIPSE
! PURPOSE
! NOTES
!------------------------------------------------------------------------------
SUBROUTINE XF_SET_ELLIPSE (a_CoordId, a_val, error)
INTEGER             , INTENT(IN) :: a_CoordId
INTEGER, INTENT(IN)        :: a_val
INTEGER, INTENT(OUT)       :: error

!            INTEGER, EXTERNAL :: _xfSetEllipse_f
!  MS FORTRAN needs explicit interface for C functions called here.

  INTERFACE
    INTEGER FUNCTION xfsetellipse_f(a_CoordId, a_val)
    
	  !MS$ATTRIBUTES C,reference::xfsetellipse_f
	  INTEGER,            INTENT(IN)   :: a_CoordId
    INTEGER,            INTENT(IN)   :: a_val


    END FUNCTION xfsetellipse_f
  END INTERFACE


  error = xfsetellipse_f (a_CoordId, a_val)


  return
END SUBROUTINE
!------------------------------------------------------------------------------
! SUBROUTINE  XF_SET_MAJOR_R
! PURPOSE
! NOTES
!------------------------------------------------------------------------------
SUBROUTINE XF_SET_MAJOR_R (a_CoordId, a_val, error)
INTEGER, INTENT(IN) :: a_CoordId
REAL*8,  INTENT(IN)   :: a_val
INTEGER, INTENT(OUT)       :: error


!            INTEGER, EXTERNAL :: _xfSetMajorR_f
!  MS FORTRAN needs explicit interface for C functions called here.

  INTERFACE
    INTEGER FUNCTION xfsetmajorr_f(a_CoordId, a_val)
    
	  !MS$ATTRIBUTES C,reference::xfsetmajorr_f
	  INTEGER, INTENT(IN)   :: a_CoordId
      REAL*8,  INTENT(IN)   :: a_val


    END FUNCTION xfsetmajorr_f
  END INTERFACE


  error = xfsetmajorr_f (a_CoordId, a_val)



  return
END SUBROUTINE
!------------------------------------------------------------------------------
! SUBROUTINE  XF_SET_MINOR_R
! PURPOSE
! NOTES
!------------------------------------------------------------------------------
SUBROUTINE XF_SET_MINOR_R (a_CoordId, a_val, error)
INTEGER, INTENT(IN)   :: a_CoordId
REAL*8,  INTENT(IN)   :: a_val
INTEGER, INTENT(OUT)  :: error


!            INTEGER, EXTERNAL :: _xfSetMinorR_f
!  MS FORTRAN needs explicit interface for C functions called here.

  INTERFACE
    INTEGER FUNCTION xfsetminorr_f(a_CoordId, a_val)
    
	  !MS$ATTRIBUTES C,reference::xfsetminorr_f
	  INTEGER, INTENT(IN)   :: a_CoordId
      REAL*8,  INTENT(IN)   :: a_val


    END FUNCTION xfsetminorr_f
  END INTERFACE


  error = xfsetminorr_f (a_CoordId, a_val)

  return
END SUBROUTINE

!-----------------------------------------------------------------------------
! SUBROUTINE XF_CALENDAR_TO_JULIAN
! PURPOSE    Convert Calendar to Julian Date
! NOTES      (Julian day number algorithm adopted from Press et al.)
!            era is #defined (use #defines): ERA_IS_BCE (BC), ERA_IS_CE (AD)
!           -Taken from JavaScript Code found at website:
!                 http://aa.usno.navy.mil/data/docs/JulianDate.html
!            Contact info. provided by website:
!              Marc A. Murison
!              Astronomical Applications Dept.
!              U.S. Naval Observatory
!              3450 Massachusetts Ave, NW
!              Washington, DC  20392-5420
!-----------------------------------------------------------------------------
SUBROUTINE XF_CALENDAR_TO_JULIAN (era, yr, mo, day, hr, &
                                  min, sec, julian, error)
  INTEGER, INTENT(IN)       :: era
  INTEGER, INTENT(IN)       :: yr
  INTEGER, INTENT(IN)       :: mo
  INTEGER, INTENT(IN)       :: day
  INTEGER, INTENT(IN)       :: hr
  INTEGER, INTENT(IN)       :: min
  INTEGER, INTENT(IN)       :: sec
  REAL*8,  INTENT(OUT)      :: julian
  INTEGER, INTENT(OUT)      :: error
  

!            INTEGER, EXTERNAL :: _xfCalendarToJulian_f
!  MS FORTRAN needs explicit interface for C functions called here.

  INTERFACE
    INTEGER FUNCTION xfcalendartojulian_f(era, yr, mo, day, hr, & 
                                          min, sec, julian)
    
	  !MS$ATTRIBUTES C,reference::xfcalendartojulian_f
    REAL*8,  INTENT(OUT)  :: julian
    INTEGER, INTENT(IN)   :: era, yr, mo, day, hr, min, sec


    END FUNCTION xfcalendartojulian_f
  END INTERFACE


  error = xfcalendartojulian_f (era, yr, mo, day, hr, min, &
                                sec, julian) 


  return

END SUBROUTINE

!-----------------------------------------------------------------------------
! SUBROUTINE XF_JULIAN_TO_CALENDAR
! PURPOSE    Convert Julian Date To Calendar
! NOTES      (algorithm adopted from Press et al.)
!            era is #defined (use #defines): ERA_IS_BCE (BC), ERA_IS_CE (AD)
!           -Taken from JavaScript Code found at website:
!                 http://aa.usno.navy.mil/data/docs/JulianDate.html
!            Contact info. provided by website:
!              Marc A. Murison
!              Astronomical Applications Dept.
!              U.S. Naval Observatory
!              3450 Massachusetts Ave, NW
!              Washington, DC  20392-5420
!-----------------------------------------------------------------------------
SUBROUTINE XF_JULIAN_TO_CALENDAR (era, yr, mo, day, hr, &
                                  min, sec, julian, error)
  INTEGER, INTENT(OUT) :: era
  INTEGER, INTENT(OUT) :: yr
  INTEGER, INTENT(OUT) :: mo
  INTEGER, INTENT(OUT) :: day
  INTEGER, INTENT(OUT) :: hr
  INTEGER, INTENT(OUT) :: min
  INTEGER, INTENT(OUT) :: sec
  REAL*8,  INTENT(IN)  :: julian
  INTEGER, INTENT(OUT) :: error
  

!            INTEGER, EXTERNAL :: _xfJulianToCalendar_f
!  MS FORTRAN needs explicit interface for C functions called here.

  INTERFACE
    INTEGER FUNCTION xfjuliantocalendar_f(era, yr, mo, day, hr, & 
                                          min, sec, julian)
    
	  !MS$ATTRIBUTES C,reference::xfjuliantocalendar_f
    REAL*8,   INTENT(IN)    :: julian
    INTEGER,  INTENT(OUT)   :: era, yr, mo, day, hr, min, sec


    END FUNCTION xfjuliantocalendar_f
  END INTERFACE


  error = xfjuliantocalendar_f (era, yr, mo, day, hr, min, &
                                sec, julian) 

  return

END SUBROUTINE

!************************************************************************************

SUBROUTINE XF_DATASET_REFTIME (a_Id, a_Reftime, error)

  INTEGER, INTENT(IN)  :: a_Id
  REAL*8,  INTENT(IN)  :: a_Reftime
  INTEGER, INTENT(OUT) :: error
  

!            INTEGER, EXTERNAL :: _xfDatasetReftime_f
!  MS FORTRAN needs explicit interface for C functions called here.

  INTERFACE
    INTEGER FUNCTION xfdatasetreftime_f(a_Id, a_Reftime)
    
	  !MS$ATTRIBUTES C,reference::xfdatasetreftime_f
    INTEGER, INTENT(IN)    :: a_Id
    REAL*8,  INTENT(IN)    :: a_Reftime
    


    END FUNCTION xfdatasetreftime_f
  END INTERFACE


  error = xfdatasetreftime_f (a_Id, a_Reftime)

  return

END SUBROUTINE

!************************************************************************************

SUBROUTINE XF_GET_NUM_OPEN_IDENTIFIERS (a_Id, a_Num, error)

INTEGER,     INTENT(IN)  :: a_Id
INTEGER,     INTENT(OUT) :: a_Num
INTEGER,     INTENT(OUT) :: error
  

!            INTEGER, EXTERNAL :: _xfGetNumOpenIdentifiers_f
!  MS FORTRAN needs explicit interface for C functions called here.

  INTERFACE
    INTEGER FUNCTION xfgetnumopenidentifiers_f (a_Id, a_Num)
    
	  !MS$ATTRIBUTES C,reference::xfgetnumopenidentifiers_f
    INTEGER,     INTENT(IN)  :: a_Id, a_Num


    END FUNCTION xfgetnumopenidentifiers_f
  END INTERFACE
  
  a_Num = 0
  error = xfgetnumopenidentifiers_f (a_Id, a_Num)
  return

END SUBROUTINE


!******************************************************************************
!* FUNCTION  XFI_CLOSE_OPEN_IDENTIFIERS
!* PURPOSE   Closes all open identifiers for a file 
!* NOTES     
!******************2***********************************************************
SUBROUTINE XFI_CLOSE_OPEN_IDENTIFIERS (a_Id, error)
INTEGER             ,       INTENT(IN) :: a_Id
INTEGER,              INTENT(OUT):: error


!            INTEGER, EXTERNAL :: _xfpCloseOpenIdentifiers_f
!  MS FORTRAN needs explicit interface for C functions called here.

  INTERFACE
    INTEGER FUNCTION xfpcloseopenidentifiers_f (a_Id)
    
	  !MS$ATTRIBUTES C,reference::xfpcloseopenidentifiers_f
    INTEGER,     INTENT(IN)  :: a_Id


    END FUNCTION xfpcloseopenidentifiers_f 
  END INTERFACE


  error = xfpcloseopenidentifiers_f (a_Id)

return
END SUBROUTINE ! XFICloseOpenIdentifiers




!******************************************************************************
!* FUNCTION  XF_SET_NUM_OF_XSECTS
!* PURPOSE   Set the number of Cross Sections
!* NOTES  
!************************************************************************************

SUBROUTINE XF_SET_NUM_OF_XSECTS (a_Id, a_nXSects, a_compression, error)
INTEGER,          INTENT(IN)  :: a_Id
INTEGER,          INTENT(IN)  :: a_nXSects, a_compression
INTEGER,          INTENT(OUT) :: error
  

!            INTEGER, EXTERNAL :: _xfsetnumberofxsects_f
!  MS FORTRAN needs explicit interface for C functions called here.

  INTERFACE
    INTEGER FUNCTION xfsetnumberofxsects_f (a_Id, a_nXSects, a_compression)
    
	  !MS$ATTRIBUTES C,reference::xfsetnumberofxsects_f
    INTEGER,          INTENT(IN)  :: a_Id
    INTEGER,          INTENT(IN)  :: a_nXSects, a_compression


    END FUNCTION xfsetnumberofxsects_f
  END INTERFACE


  error = xfsetnumberofxsects_f (a_Id, a_nXSects, a_compression)

  return

END SUBROUTINE ! XF_SET_NUM_OF_XSECTS

!******************************************************************************
!* FUNCTION  XF_GET_NUM_OF_XSECTS
!* PURPOSE   Set the number of Cross Sections
!* NOTES  
!************************************************************************************

SUBROUTINE XF_GET_NUM_OF_XSECTS (a_Id, a_nXSects, error)
INTEGER,          INTENT(IN)  :: a_Id
INTEGER,          INTENT(IN)  :: a_nXSects
INTEGER,          INTENT(OUT) :: error
  

!            INTEGER, EXTERNAL :: _xfgetnumberofxsects_f
!  MS FORTRAN needs explicit interface for C functions called here.

  INTERFACE
    INTEGER FUNCTION xfgetnumberofxsects_f (a_Id, a_nXSects)
    
	  !MS$ATTRIBUTES C,reference::xfgetnumberofxsects_f
    INTEGER,          INTENT(IN)  :: a_Id
    INTEGER,          INTENT(IN)  :: a_nXSects


    END FUNCTION xfgetnumberofxsects_f
  END INTERFACE


  error = xfgetnumberofxsects_f (a_Id, a_nXSects)

  return

END SUBROUTINE ! XF_GET_NUM_OF_XSECTS


!******************************************************************************
!* FUNCTION  XF_SET_CS_ID
!* PURPOSE   Set the Cross Sections ID
!* NOTES  
!************************************************************************************

SUBROUTINE XF_SET_CS_ID (a_Id, a_NumVals, a_PropId, a_compression, error)
INTEGER,          INTENT(IN)  :: a_Id
INTEGER,          INTENT(IN)  :: a_NumVals, a_PropId, a_compression
INTEGER,          INTENT(OUT) :: error
  

!            INTEGER, EXTERNAL :: _xfsetcsid_f
!  MS FORTRAN needs explicit interface for C functions called here.

  INTERFACE
    INTEGER FUNCTION xfsetcsid_f (a_Id, a_NumVals, a_PropId, a_compression)
    
	  !MS$ATTRIBUTES C,reference::xfsetcsid_f
    INTEGER,          INTENT(IN)  :: a_Id
    INTEGER,          INTENT(IN)  :: a_NumVals, a_PropId, a_compression


    END FUNCTION xfsetcsid_f
  END INTERFACE


  error = xfsetcsid_f (a_Id, a_NumVals, a_PropId, a_compression)

  return

END SUBROUTINE ! XF_SET_CS_ID

!******************************************************************************
!* FUNCTION  XF_GET_CS_ID
!* PURPOSE   Get the Cross Sections ID
!* NOTES  
!************************************************************************************

SUBROUTINE XF_GET_CS_ID (a_Id, a_NumVals, a_PropId, error)
INTEGER,          INTENT(IN)  :: a_Id
INTEGER,          INTENT(IN)  :: a_NumVals, a_PropId
INTEGER,          INTENT(OUT) :: error
  

!            INTEGER, EXTERNAL :: _xfgetcsid_f
!  MS FORTRAN needs explicit interface for C functions called here.

  INTERFACE
    INTEGER FUNCTION xfgetcsid_f (a_Id, a_NumVals, a_PropId)
    
	  !MS$ATTRIBUTES C,reference::xfgetcsid_f
    INTEGER,          INTENT(IN)  :: a_Id
    INTEGER,          INTENT(IN)  :: a_NumVals, a_PropId


    END FUNCTION xfgetcsid_f
  END INTERFACE


  error = xfgetcsid_f (a_Id, a_NumVals, a_PropId)

  return

END SUBROUTINE ! XF_GET_CS_ID


!******************************************************************************
!* FUNCTION  XF_SET_CS_NAME
!* PURPOSE   Set the Cross Sections Name
!* NOTES  
!************************************************************************************

SUBROUTINE XF_SET_CS_NAME (a_Id, a_NumVals, a_PropId, a_StrLen, a_compression, error)
INTEGER,          INTENT(IN)  :: a_Id
INTEGER,          INTENT(IN)  :: a_NumVals, a_StrLen, a_compression
CHARACTER(LEN=*), INTENT(IN)  :: a_PropId
INTEGER,          INTENT(OUT) :: error
  

!            INTEGER, EXTERNAL :: _xfsetcsname_f
!  MS FORTRAN needs explicit interface for C functions called here.

  INTERFACE
    INTEGER FUNCTION xfsetcsname_f (a_Id, a_NumVals, a_PropId, a_StrLen, a_compression)
    
	  !MS$ATTRIBUTES C,reference::xfsetcsname_f
    !DEC$ATTRIBUTES reference :: a_PropId
    INTEGER,          INTENT(IN)  :: a_Id
    INTEGER,          INTENT(IN)  :: a_NumVals, a_StrLen, a_compression
    CHARACTER(LEN=*), INTENT(IN)  :: a_PropId


    END FUNCTION xfsetcsname_f
  END INTERFACE


  error = xfsetcsname_f (a_Id, a_NumVals, a_PropId, a_StrLen, a_compression)

  return

END SUBROUTINE ! XF_SET_CS_NAME

!******************************************************************************
!* FUNCTION  XF_GET_CS_NAME
!* PURPOSE   Get the Cross Sections Name
!* NOTES  
!************************************************************************************

SUBROUTINE XF_GET_CS_NAME (a_Id, a_NumVals, a_PropId, a_StrLen, error)
INTEGER,          INTENT(IN)  :: a_Id
INTEGER,          INTENT(IN)  :: a_NumVals, a_StrLen
CHARACTER(LEN=*), INTENT(IN)  :: a_PropId
INTEGER,          INTENT(OUT) :: error
  

!            INTEGER, EXTERNAL :: _xfgetcsid_f
!  MS FORTRAN needs explicit interface for C functions called here.

  INTERFACE
    INTEGER FUNCTION xfgetcsname_f (a_Id, a_NumVals, a_PropId, a_StrLen)
    
	  !MS$ATTRIBUTES C,reference::xfgetcsname_f
    !DEC$ATTRIBUTES reference :: a_PropId
    INTEGER,          INTENT(IN)  :: a_Id
    INTEGER,          INTENT(IN)  :: a_NumVals, a_StrLen
    CHARACTER(LEN=*), INTENT(IN)  :: a_PropId


    END FUNCTION xfgetcsname_f
  END INTERFACE


  error = xfgetcsname_f (a_Id, a_NumVals, a_PropId, a_StrLen)

  return

END SUBROUTINE ! XF_GET_CS_NAME


!******************************************************************************
!* FUNCTION  XF_GET_CS_NAME_LEN
!* PURPOSE   Get the Cross Sections Name Length
!* NOTES  
!************************************************************************************

SUBROUTINE XF_GET_CS_NAME_LEN (a_Id, a_NumVals, a_StrLen, error)
INTEGER,          INTENT(IN)  :: a_Id
INTEGER,          INTENT(IN)  :: a_NumVals, a_StrLen
INTEGER,          INTENT(OUT) :: error
  

!            INTEGER, EXTERNAL :: _xfgetcsnamelen_f
!  MS FORTRAN needs explicit interface for C functions called here.

  INTERFACE
    INTEGER FUNCTION xfgetcsnamelen_f (a_Id, a_NumVals, a_StrLen)
    
	  !MS$ATTRIBUTES C,reference::xfgetcsnamelen_f
    INTEGER,          INTENT(IN)  :: a_Id
    INTEGER,          INTENT(IN)  :: a_NumVals, a_StrLen


    END FUNCTION xfgetcsnamelen_f
  END INTERFACE


  error = xfgetcsnamelen_f (a_Id, a_NumVals, a_StrLen)

  return

END SUBROUTINE ! XF_GET_CS_NAME_LEN


!******************************************************************************
!* FUNCTION  XF_SET_REACH_NAME
!* PURPOSE   Set the Reach Name
!* NOTES  
!************************************************************************************

SUBROUTINE XF_SET_REACH_NAME (a_Id, a_NumVals, a_PropId, a_StrLen, error)
INTEGER,          INTENT(IN)  :: a_Id
INTEGER,          INTENT(IN)  :: a_NumVals, a_StrLen
CHARACTER(LEN=*), INTENT(IN)  :: a_PropId
INTEGER,          INTENT(OUT) :: error
  

!            INTEGER, EXTERNAL :: _xfsetreachname_f
!  MS FORTRAN needs explicit interface for C functions called here.

  INTERFACE
    INTEGER FUNCTION xfsetreachname_f (a_Id, a_NumVals, a_PropId, a_StrLen)
    
	  !MS$ATTRIBUTES C,reference::xfsetreachname_f
    !DEC$ATTRIBUTES reference :: a_PropId
    INTEGER,          INTENT(IN)  :: a_Id
    INTEGER,          INTENT(IN)  :: a_NumVals, a_StrLen
    CHARACTER(LEN=*), INTENT(IN)  :: a_PropId


    END FUNCTION xfsetreachname_f
  END INTERFACE


  error = xfsetreachname_f (a_Id, a_NumVals, a_PropId, a_StrLen)

  return

END SUBROUTINE ! XF_SET_REACH_NAME

!******************************************************************************
!* FUNCTION  XF_GET_REACH_NAME
!* PURPOSE   Get the Reach Name
!* NOTES  
!************************************************************************************

SUBROUTINE XF_GET_REACH_NAME (a_Id, a_NumVals, a_PropId, a_StrLen, error)
INTEGER,          INTENT(IN)  :: a_Id
INTEGER,          INTENT(IN)  :: a_NumVals, a_StrLen
CHARACTER(LEN=*), INTENT(IN)  :: a_PropId
INTEGER,          INTENT(OUT) :: error
  

!            INTEGER, EXTERNAL :: _xfsetreachname_f
!  MS FORTRAN needs explicit interface for C functions called here.

  INTERFACE
    INTEGER FUNCTION xfsetreachname_f (a_Id, a_NumVals, a_PropId, a_StrLen)
    
	  !MS$ATTRIBUTES C,reference::xfsetreachname_f
    !DEC$ATTRIBUTES reference :: a_PropId
    INTEGER,          INTENT(IN)  :: a_Id
    INTEGER,          INTENT(IN)  :: a_NumVals, a_StrLen
    CHARACTER(LEN=*), INTENT(IN)  :: a_PropId


    END FUNCTION xfsetreachname_f
  END INTERFACE


  error = xfsetreachname_f (a_Id, a_NumVals, a_PropId, a_StrLen)

  return

END SUBROUTINE ! XF_GET_REACH_NAME


!******************************************************************************
!* FUNCTION  XF_GET_REACH_NAME_LEN
!* PURPOSE   Get the Reach Name Length
!* NOTES  
!************************************************************************************

SUBROUTINE XF_GET_REACH_NAME_LEN (a_Id, a_NumVals, a_StrLen, error)
INTEGER,          INTENT(IN)  :: a_Id
INTEGER,          INTENT(IN)  :: a_NumVals, a_StrLen
INTEGER,          INTENT(OUT) :: error
  

!            INTEGER, EXTERNAL :: _xfGetReachNameLen_f
!  MS FORTRAN needs explicit interface for C functions called here.

  INTERFACE
    INTEGER FUNCTION xfgetreachnamelen_f (a_Id, a_NumVals, a_StrLen)
    
	  !MS$ATTRIBUTES C,reference::xfgetreachnamelen_f
    INTEGER,          INTENT(IN)  :: a_Id
    INTEGER,          INTENT(IN)  :: a_NumVals, a_StrLen


    END FUNCTION xfgetreachnamelen_f
  END INTERFACE


  error = xfgetreachnamelen_f (a_Id, a_NumVals, a_StrLen)

  return

END SUBROUTINE ! XF_GET_REACH_NAME_LEN


!******************************************************************************
!* FUNCTION  XF_SET_TOPO_ID
!* PURPOSE   Set the Topo Id
!* NOTES  
!************************************************************************************

SUBROUTINE XF_SET_TOPO_ID (a_Id, a_NumVals, a_PropId, a_StrLen, a_compression, error)
INTEGER,          INTENT(IN)  :: a_Id
INTEGER,          INTENT(IN)  :: a_NumVals, a_StrLen, a_compression
CHARACTER(LEN=*), INTENT(IN)  :: a_PropId
INTEGER,          INTENT(OUT) :: error
  

!            INTEGER, EXTERNAL :: _xfsettopoid_f
!  MS FORTRAN needs explicit interface for C functions called here.

  INTERFACE
    INTEGER FUNCTION xfsettopoid_f (a_Id, a_NumVals, a_PropId, a_StrLen, a_compression)
    
	  !MS$ATTRIBUTES C,reference::xfsettopoid_f
    !DEC$ATTRIBUTES reference :: a_propid
    INTEGER,          INTENT(IN)  :: a_Id
    INTEGER,          INTENT(IN)  :: a_NumVals, a_StrLen, a_compression
    CHARACTER(LEN=*), INTENT(IN)  :: a_PropId


    END FUNCTION xfsettopoid_f
  END INTERFACE


  error = xfsettopoid_f (a_Id, a_NumVals, a_PropId, a_StrLen, a_compression)

  return

END SUBROUTINE ! XF_SET_TOPO_ID

!******************************************************************************
!* FUNCTION  XF_GET_TOPO_ID
!* PURPOSE   Get the Topo Id
!* NOTES  
!************************************************************************************

SUBROUTINE XF_GET_TOPO_ID (a_Id, a_NumVals, a_PropId, a_StrLen, error)
INTEGER,          INTENT(IN)  :: a_Id
INTEGER,          INTENT(IN)  :: a_NumVals, a_StrLen
CHARACTER(LEN=*), INTENT(IN)  :: a_PropId
INTEGER,          INTENT(OUT) :: error
  

!            INTEGER, EXTERNAL :: _xfgettopoid_f
!  MS FORTRAN needs explicit interface for C functions called here.

  INTERFACE
    INTEGER FUNCTION xfgettopoid_f (a_Id, a_NumVals, a_StrLen, a_PropId)
    
	  !MS$ATTRIBUTES C,reference::xfgettopoid_f
    !DEC$ATTRIBUTES reference :: a_PropId
    INTEGER,          INTENT(IN)  :: a_Id
    INTEGER,          INTENT(IN)  :: a_NumVals, a_StrLen
    CHARACTER(LEN=*), INTENT(IN)  :: a_PropId


    END FUNCTION xfgettopoid_f
  END INTERFACE


  error = xfgettopoid_f (a_Id, a_NumVals, a_StrLen, a_PropId)

  return

END SUBROUTINE ! XF_GET_TOPO_ID


!******************************************************************************
!* FUNCTION  XF_GET_TOPO_ID_LEN
!* PURPOSE   Get the Topo Id String Length
!* NOTES  
!************************************************************************************

SUBROUTINE XF_GET_TOPO_ID_LEN (a_Id, a_NumVals, a_StrLen, error)
INTEGER,          INTENT(IN)  :: a_Id
INTEGER,          INTENT(IN)  :: a_NumVals, a_StrLen
INTEGER,          INTENT(OUT) :: error
  

!            INTEGER, EXTERNAL :: _xfgettopoidlen_f
!  MS FORTRAN needs explicit interface for C functions called here.

  INTERFACE
    INTEGER FUNCTION xfgettopoidlen_f (a_Id, a_NumVals, a_StrLen)
    
	  !MS$ATTRIBUTES C,reference::xfgettopoidlen_f
    INTEGER,          INTENT(IN)  :: a_Id
    INTEGER,          INTENT(IN)  :: a_NumVals, a_StrLen


    END FUNCTION xfgettopoidlen_f
  END INTERFACE


  error = xfgettopoidlen_f (a_Id, a_NumVals, a_StrLen)

  return

END SUBROUTINE ! XF_GET_TOPO_ID_LEN


!******************************************************************************
!* FUNCTION  XF_SET_STATION
!* PURPOSE   Set the Station
!* NOTES  
!************************************************************************************

SUBROUTINE XF_SET_STATION (a_Id, a_NumVals, a_PropId, a_compression, error)
INTEGER,          INTENT(IN)  :: a_Id
INTEGER,          INTENT(IN)  :: a_NumVals, a_compression
REAL*4,             INTENT(IN)  :: a_PropId
INTEGER,          INTENT(OUT) :: error
  

!            INTEGER, EXTERNAL :: _xfsetstation_f
!  MS FORTRAN needs explicit interface for C functions called here.

  INTERFACE
    INTEGER FUNCTION xfsetstation_f (a_Id, a_NumVals, a_PropId, a_compression)
    
	  !MS$ATTRIBUTES C,reference::xfsetstation_f
    INTEGER,          INTENT(IN)  :: a_Id
    INTEGER,          INTENT(IN)  :: a_NumVals, a_compression
    REAL*4,             INTENT(IN)  :: a_PropId


    END FUNCTION xfsetstation_f
  END INTERFACE


  error = xfsetstation_f (a_Id, a_NumVals, a_PropId, a_compression)

  return

END SUBROUTINE ! XF_SET_STATION

!******************************************************************************
!* FUNCTION  XF_GET_STATION
!* PURPOSE   Get the Station
!* NOTES  
!************************************************************************************

SUBROUTINE XF_GET_STATION (a_Id, a_NumVals, a_PropId, error)
INTEGER,          INTENT(IN)  :: a_Id
INTEGER,          INTENT(IN)  :: a_NumVals
REAL*4,             INTENT(IN)  :: a_PropId
INTEGER,          INTENT(OUT) :: error
  

!            INTEGER, EXTERNAL :: _xfgetstation_f
!  MS FORTRAN needs explicit interface for C functions called here.

  INTERFACE
    INTEGER FUNCTION xfgetstation_f (a_Id, a_NumVals, a_PropId)
    
	  !MS$ATTRIBUTES C,reference::xfgetstation_f
    INTEGER,          INTENT(IN)  :: a_Id
    INTEGER,          INTENT(IN)  :: a_NumVals
    REAL*4,             INTENT(IN)  :: a_PropId


    END FUNCTION xfgetstation_f
  END INTERFACE


  error = xfgetstation_f (a_Id, a_NumVals, a_PropId)

  return

END SUBROUTINE ! XF_GET_STATION

!******************************************************************************
!* FUNCTION  XF_SET_TYPE
!* PURPOSE   Set the Type
!* NOTES  
!************************************************************************************

SUBROUTINE XF_SET_TYPE (a_Id, a_NumVals, a_PropId, a_compression, error)
INTEGER,          INTENT(IN)  :: a_Id
INTEGER,          INTENT(IN)  :: a_NumVals, a_PropId, a_compression
INTEGER,          INTENT(OUT) :: error
  

!            INTEGER, EXTERNAL :: _xfsettype_f
!  MS FORTRAN needs explicit interface for C functions called here.

  INTERFACE
    INTEGER FUNCTION xfsettype_f (a_Id, a_NumVals, a_PropId, a_compression)
    
	  !MS$ATTRIBUTES C,reference::xfsettype_f
    INTEGER,          INTENT(IN)  :: a_Id
    INTEGER,          INTENT(IN)  :: a_NumVals, a_PropId, a_compression


    END FUNCTION xfsettype_f
  END INTERFACE


  error = xfsettype_f (a_Id, a_NumVals, a_PropId, a_compression)

  return

END SUBROUTINE ! XF_SET_TYPE

!******************************************************************************
!* FUNCTION  XF_GET_TYPE
!* PURPOSE   Get the Type
!* NOTES  
!************************************************************************************

SUBROUTINE XF_GET_TYPE (a_Id, a_NumVals, a_PropId, error)
INTEGER,          INTENT(IN)  :: a_Id
INTEGER,          INTENT(IN)  :: a_NumVals, a_PropId
INTEGER,          INTENT(OUT) :: error
  

!            INTEGER, EXTERNAL :: _xfgettype_f
!  MS FORTRAN needs explicit interface for C functions called here.

  INTERFACE
    INTEGER FUNCTION xfgettype_f (a_Id, a_NumVals, a_PropId)
    
	  !MS$ATTRIBUTES C,reference::xfgettype_f
    INTEGER,          INTENT(IN)  :: a_Id
    INTEGER,          INTENT(IN)  :: a_NumVals, a_PropId


    END FUNCTION xfgettype_f
  END INTERFACE


  error = xfgettype_f (a_Id, a_NumVals, a_PropId)

  return

END SUBROUTINE ! XF_GET_TYPE

!******************************************************************************
!* FUNCTION  XF_SET_PTYPE
!* PURPOSE   Set the Type
!* NOTES  
!************************************************************************************

SUBROUTINE XF_SET_PTYPE (a_Id, a_NumVals, a_PropId, a_compression, error)
INTEGER,          INTENT(IN)  :: a_Id
INTEGER,          INTENT(IN)  :: a_NumVals, a_PropId, a_compression
INTEGER,          INTENT(OUT) :: error
  

!            INTEGER, EXTERNAL :: _xfsetptype_f
!  MS FORTRAN needs explicit interface for C functions called here.

  INTERFACE
    INTEGER FUNCTION xfsetptype_f (a_Id, a_NumVals, a_PropId, a_compression)
    
	  !MS$ATTRIBUTES C,reference::xfsetptype_f
    INTEGER,          INTENT(IN)  :: a_Id
    INTEGER,          INTENT(IN)  :: a_NumVals, a_PropId, a_compression


    END FUNCTION xfsetptype_f
  END INTERFACE


  error = xfsetptype_f (a_Id, a_NumVals, a_PropId, a_compression)

  return

END SUBROUTINE ! XF_SET_PTYPE

!******************************************************************************
!* FUNCTION  XF_GET_PTYPE
!* PURPOSE   Get the Type
!* NOTES  
!************************************************************************************

SUBROUTINE XF_GET_PTYPE (a_Id, a_NumVals, a_PropId, error)
INTEGER,          INTENT(IN)  :: a_Id
INTEGER,          INTENT(IN)  :: a_NumVals, a_PropId
INTEGER,          INTENT(OUT) :: error
  

!            INTEGER, EXTERNAL :: _xfgettype_f
!  MS FORTRAN needs explicit interface for C functions called here.

  INTERFACE
    INTEGER FUNCTION xfgettype_f (a_Id, a_NumVals, a_PropId)
    
	  !MS$ATTRIBUTES C,reference::xfgettype_f
    INTEGER,          INTENT(IN)  :: a_Id
    INTEGER,          INTENT(IN)  :: a_NumVals, a_PropId


    END FUNCTION xfgettype_f
  END INTERFACE


  error = xfgettype_f (a_Id, a_NumVals, a_PropId)

  return

END SUBROUTINE ! XF_GET_PTYPE


!******************************************************************************
!* FUNCTION  XF_SET_P_CS_DB_LINK
!* PURPOSE   Set the Cross Section Database Link
!* NOTES  
!************************************************************************************

SUBROUTINE XF_SET_P_CS_DB_LINK (a_Id, a_NumVals, a_PropId, a_compression, error)
INTEGER,          INTENT(IN)  :: a_Id
INTEGER,          INTENT(IN)  :: a_NumVals, a_PropId, a_compression
INTEGER,          INTENT(OUT) :: error
  

!            INTEGER, EXTERNAL :: _xfsetpcsdblink_f
!  MS FORTRAN needs explicit interface for C functions called here.

  INTERFACE
    INTEGER FUNCTION xfsetpcsdblink_f (a_Id, a_NumVals, a_PropId, a_compression)
    
	  !MS$ATTRIBUTES C,reference::xfsetpcsdblink_f
    INTEGER,          INTENT(IN)  :: a_Id
    INTEGER,          INTENT(IN)  :: a_NumVals, a_PropId, a_compression


    END FUNCTION xfsetpcsdblink_f
  END INTERFACE


  error = xfsetpcsdblink_f (a_Id, a_NumVals, a_PropId, a_compression)

  return

END SUBROUTINE ! XF_SET_P_CS_DB_LINK

!******************************************************************************
!* FUNCTION  XF_GET_P_CS_DB_LINK
!* PURPOSE   Get the Cross Section Database Link
!* NOTES  
!************************************************************************************

SUBROUTINE XF_GET_P_CS_DB_LINK (a_Id, a_NumVals, a_PropId, error)
INTEGER,          INTENT(IN)  :: a_Id
INTEGER,          INTENT(IN)  :: a_NumVals, a_PropId
INTEGER,          INTENT(OUT) :: error
  

!            INTEGER, EXTERNAL :: _xfgetpcsdblink_f
!  MS FORTRAN needs explicit interface for C functions called here.

  INTERFACE
    INTEGER FUNCTION xfgetpcsdblink_f (a_Id, a_NumVals, a_PropId)
    
	  !MS$ATTRIBUTES C,reference::xfgetpcsdblink_f
    INTEGER,          INTENT(IN)  :: a_Id
    INTEGER,          INTENT(IN)  :: a_NumVals, a_PropId


    END FUNCTION xfgetpcsdblink_f
  END INTERFACE


  error = xfgetpcsdblink_f (a_Id, a_NumVals, a_PropId)

  return

END SUBROUTINE ! XF_GET_P_CS_DB_LINK


!******************************************************************************
!* FUNCTION  XF_SET_NOTE
!* PURPOSE   Set the Cross Sections Note
!* NOTES  
!************************************************************************************

SUBROUTINE XF_SET_NOTE (a_Id, a_NumVals, a_PropId, a_StrLen, a_compression, error)
INTEGER,          INTENT(IN)  :: a_Id
INTEGER,          INTENT(IN)  :: a_NumVals, a_StrLen, a_compression
CHARACTER(LEN=*), INTENT(IN)  :: a_PropId
INTEGER,          INTENT(OUT) :: error
  

!            INTEGER, EXTERNAL :: _xfsetnote_f
!  MS FORTRAN needs explicit interface for C functions called here.

  INTERFACE
    INTEGER FUNCTION xfsetnote_f (a_Id, a_NumVals, a_PropId, a_StrLen, a_compression)
    
	  !MS$ATTRIBUTES C,reference::xfsetnote_f
    !DEC$ATTRIBUTES reference :: a_PropId
    INTEGER,          INTENT(IN)  :: a_Id
    INTEGER,          INTENT(IN)  :: a_NumVals, a_StrLen, a_compression
    CHARACTER(LEN=*), INTENT(IN)  :: a_PropId


    END FUNCTION xfsetnote_f
  END INTERFACE


  error = xfsetnote_f (a_Id, a_NumVals, a_PropId, a_StrLen, a_compression)

  return

END SUBROUTINE ! XF_SET_NOTE

!******************************************************************************
!* FUNCTION  XF_GET_NOTE
!* PURPOSE   Get the Cross Sections Note
!* NOTES  
!************************************************************************************

SUBROUTINE XF_GET_NOTE (a_Id, a_NumVals, a_PropId, a_StrLen, error)
INTEGER,          INTENT(IN)  :: a_Id
INTEGER,          INTENT(IN)  :: a_NumVals, a_StrLen
CHARACTER(LEN=*), INTENT(IN)  :: a_PropId
INTEGER,          INTENT(OUT) :: error
  

!            INTEGER, EXTERNAL :: _xfgetnote_f
!  MS FORTRAN needs explicit interface for C functions called here.

  INTERFACE
    INTEGER FUNCTION xfgetnote_f (a_Id, a_NumVals, a_PropId, a_StrLen)
    
	  !MS$ATTRIBUTES C,reference::xfgetnote_f
    !DEC$ATTRIBUTES reference :: a_PropId
    INTEGER,          INTENT(IN)  :: a_Id
    INTEGER,          INTENT(IN)  :: a_NumVals, a_StrLen
    CHARACTER(LEN=*), INTENT(IN)  :: a_PropId


    END FUNCTION xfgetnote_f
  END INTERFACE


  error = xfgetnote_f (a_Id, a_NumVals, a_PropId, a_StrLen)

  return

END SUBROUTINE ! XF_GET_NOTE


!******************************************************************************
!* FUNCTION  XF_GET_NOTE_LEN
!* PURPOSE   Get the Cross Sections Note Length
!* NOTES  
!************************************************************************************

SUBROUTINE XF_GET_NOTE_LEN (a_Id, a_NumVals, a_StrLen, error)
INTEGER,          INTENT(IN)  :: a_Id
INTEGER,          INTENT(IN)  :: a_NumVals, a_StrLen
INTEGER,          INTENT(OUT) :: error
  

!            INTEGER, EXTERNAL :: _xfgetnotelen_f
!  MS FORTRAN needs explicit interface for C functions called here.

  INTERFACE
    INTEGER FUNCTION xfgetnotelen_f (a_Id, a_NumVals, a_StrLen)
    
	  !MS$ATTRIBUTES C,reference::xfgetnotelen_f
    INTEGER,          INTENT(IN)  :: a_Id
    INTEGER,          INTENT(IN)  :: a_NumVals, a_StrLen


    END FUNCTION xfgetnotelen_f
  END INTERFACE


  error = xfgetnotelen_f (a_Id, a_NumVals, a_StrLen)

  return

END SUBROUTINE ! XF_GET_NOTE_LEN

!******************************************************************************
!* FUNCTION  XF_SET_XSECT_GEOM_X
!* PURPOSE   Set the Cross Section Geometry X
!* NOTES  
!************************************************************************************

SUBROUTINE XF_SET_XSECT_GEOM_X (a_Id, a_index, a_NumVals, a_iValues, a_compression, error)
INTEGER,          INTENT(IN)  :: a_Id
INTEGER,          INTENT(IN)  :: a_NumVals, a_index, a_compression
REAL*4,             INTENT(IN)  :: a_iValues
INTEGER,          INTENT(OUT) :: error
  

!            INTEGER, EXTERNAL :: _xfsetxsectgeomx_f
!  MS FORTRAN needs explicit interface for C functions called here.

  INTERFACE
    INTEGER FUNCTION xfsetxsectgeomx_f (a_Id, a_index, a_NumVals, a_iValues, a_compression)
    
	  !MS$ATTRIBUTES C,reference::xfsetxsectgeomx_f
    INTEGER,          INTENT(IN)  :: a_Id
    INTEGER,          INTENT(IN)  :: a_NumVals, a_index, a_compression
    REAL*4,             INTENT(IN)  :: a_iValues


    END FUNCTION xfsetxsectgeomx_f
  END INTERFACE


  error = xfsetxsectgeomx_f (a_Id, a_index, a_NumVals, a_iValues, a_compression)

  return

END SUBROUTINE ! XF_SET_XSECT_GEOM_X

!******************************************************************************
!* FUNCTION  XF_GET_XSECT_GEOM_X
!* PURPOSE   Get the Cross Section Geometry X
!* NOTES  
!************************************************************************************

SUBROUTINE XF_GET_XSECT_GEOM_X (a_Id, a_index, a_NumVals, a_iValues, error)
INTEGER,          INTENT(IN)  :: a_Id
INTEGER,          INTENT(IN)  :: a_NumVals, a_index
REAL*4,             INTENT(IN)  :: a_iValues
INTEGER,          INTENT(OUT) :: error
  

!            INTEGER, EXTERNAL :: _xfgetxsectgeomx_f
!  MS FORTRAN needs explicit interface for C functions called here.

  INTERFACE
    INTEGER FUNCTION xfgetxsectgeomx_f (a_Id, a_index, a_NumVals, a_iValues)
    
	  !MS$ATTRIBUTES C,reference::xfgetxsectgeomx_f
    INTEGER,          INTENT(IN)  :: a_Id
    INTEGER,          INTENT(IN)  :: a_NumVals, a_index
    REAL*4,             INTENT(IN)  :: a_iValues


    END FUNCTION xfgetxsectgeomx_f
  END INTERFACE


  error = xfgetxsectgeomx_f (a_Id, a_index, a_NumVals, a_iValues)

  return

END SUBROUTINE ! XF_GET_XSECT_GEOM_X

!******************************************************************************
!* FUNCTION  XF_SET_XSECT_GEOM_Y
!* PURPOSE   Set the Cross Section Geometry Y
!* NOTES  
!************************************************************************************

SUBROUTINE XF_SET_XSECT_GEOM_Y (a_Id, a_index, a_NumVals, a_iValues, a_compression, error)
INTEGER,          INTENT(IN)  :: a_Id
INTEGER,          INTENT(IN)  :: a_NumVals, a_index, a_compression
REAL*4,             INTENT(IN)  :: a_iValues
INTEGER,          INTENT(OUT) :: error
  

!            INTEGER, EXTERNAL :: _xfsetxsectgeomy_f
!  MS FORTRAN needs explicit interface for C functions called here.

  INTERFACE
    INTEGER FUNCTION xfsetxsectgeomy_f (a_Id, a_index, a_NumVals, a_iValues, a_compression)
    
	  !MS$ATTRIBUTES C,reference::xfsetxsectgeomy_f
    INTEGER,          INTENT(IN)  :: a_Id
    INTEGER,          INTENT(IN)  :: a_NumVals, a_index, a_compression
    REAL*4,             INTENT(IN)  :: a_iValues


    END FUNCTION xfsetxsectgeomy_f
  END INTERFACE


  error = xfsetxsectgeomy_f (a_Id, a_index, a_NumVals, a_iValues, a_compression)

  return

END SUBROUTINE ! XF_SET_XSECT_GEOM_Y

!******************************************************************************
!* FUNCTION  XF_GET_XSECT_GEOM_Y
!* PURPOSE   Get the Cross Section Geometry Y
!* NOTES  
!************************************************************************************

SUBROUTINE XF_GET_XSECT_GEOM_Y (a_Id, a_index, a_NumVals, a_iValues, error)
INTEGER,          INTENT(IN)  :: a_Id
INTEGER,          INTENT(IN)  :: a_NumVals, a_index
REAL*4,             INTENT(IN)  :: a_iValues
INTEGER,          INTENT(OUT) :: error
  

!            INTEGER, EXTERNAL :: _xfgetxsectgeomy_f
!  MS FORTRAN needs explicit interface for C functions called here.

  INTERFACE
    INTEGER FUNCTION xfgetxsectgeomy_f (a_Id, a_index, a_NumVals, a_iValues)
    
	  !MS$ATTRIBUTES C,reference::xfgetxsectgeomy_f
    INTEGER,          INTENT(IN)  :: a_Id
    INTEGER,          INTENT(IN)  :: a_NumVals, a_index
    REAL*4,             INTENT(IN)  :: a_iValues


    END FUNCTION xfgetxsectgeomy_f
  END INTERFACE


  error = xfgetxsectgeomy_f (a_Id, a_index, a_NumVals, a_iValues)

  return

END SUBROUTINE ! XF_GET_XSECT_GEOM_Y


!******************************************************************************
!* FUNCTION  XF_SET_XSECT_GEOM_D
!* PURPOSE   Set the Cross Section Geometry D
!* NOTES  
!************************************************************************************

SUBROUTINE XF_SET_XSECT_GEOM_D (a_Id, a_index, a_NumVals, a_iValues, a_compression, error)
INTEGER,          INTENT(IN)  :: a_Id
INTEGER,          INTENT(IN)  :: a_NumVals, a_index, a_compression
REAL*4,             INTENT(IN)  :: a_iValues
INTEGER,          INTENT(OUT) :: error
  

!            INTEGER, EXTERNAL :: _xfsetxsectgeomd_f
!  MS FORTRAN needs explicit interface for C functions called here.

  INTERFACE
    INTEGER FUNCTION xfsetxsectgeomd_f (a_Id, a_index, a_NumVals, a_iValues, a_compression)
    
	  !MS$ATTRIBUTES C,reference::xfsetxsectgeomd_f
    INTEGER,          INTENT(IN)  :: a_Id
    INTEGER,          INTENT(IN)  :: a_NumVals, a_index, a_compression
    REAL*4,             INTENT(IN)  :: a_iValues


    END FUNCTION xfsetxsectgeomd_f
  END INTERFACE


  error = xfsetxsectgeomd_f (a_Id, a_index, a_NumVals, a_iValues, a_compression)

  return

END SUBROUTINE ! XF_SET_XSECT_GEOM_D

!******************************************************************************
!* FUNCTION  XF_GET_XSECT_GEOM_D
!* PURPOSE   Get the Cross Section Geometry D
!* NOTES  
!************************************************************************************

SUBROUTINE XF_GET_XSECT_GEOM_D (a_Id, a_index, a_NumVals, a_iValues, error)
INTEGER,          INTENT(IN)  :: a_Id
INTEGER,          INTENT(IN)  :: a_NumVals, a_index
REAL*4,             INTENT(IN)  :: a_iValues
INTEGER,          INTENT(OUT) :: error
  

!            INTEGER, EXTERNAL :: _xfgetxsectgeomd_f
!  MS FORTRAN needs explicit interface for C functions called here.

  INTERFACE
    INTEGER FUNCTION xfgetxsectgeomd_f (a_Id, a_index, a_NumVals, a_iValues)
    
	  !MS$ATTRIBUTES C,reference::xfgetxsectgeomd_f
    INTEGER,          INTENT(IN)  :: a_Id
    INTEGER,          INTENT(IN)  :: a_NumVals, a_index
    REAL*4,             INTENT(IN)  :: a_iValues


    END FUNCTION xfgetxsectgeomd_f
  END INTERFACE


  error = xfgetxsectgeomd_f (a_Id, a_index, a_NumVals, a_iValues)

  return

END SUBROUTINE ! XF_GET_XSECT_GEOM_D


!******************************************************************************
!* FUNCTION  XF_SET_XSECT_GEOM_Z
!* PURPOSE   Set the Cross Section Geometry Z
!* NOTES  
!************************************************************************************

SUBROUTINE XF_SET_XSECT_GEOM_Z (a_Id, a_index, a_NumVals, a_iValues, a_compression, error)
INTEGER,          INTENT(IN)  :: a_Id
INTEGER,          INTENT(IN)  :: a_NumVals, a_index, a_compression
REAL*4,             INTENT(IN)  :: a_iValues
INTEGER,          INTENT(OUT) :: error
  

!            INTEGER, EXTERNAL :: _xfsetxsectgeomz_f
!  MS FORTRAN needs explicit interface for C functions called here.

  INTERFACE
    INTEGER FUNCTION xfsetxsectgeomz_f (a_Id, a_index, a_NumVals, a_iValues, a_compression)
    
	  !MS$ATTRIBUTES C,reference::xfsetxsectgeomz_f
    INTEGER,          INTENT(IN)  :: a_Id
    INTEGER,          INTENT(IN)  :: a_NumVals, a_index, a_compression
    REAL*4,             INTENT(IN)  :: a_iValues


    END FUNCTION xfsetxsectgeomz_f
  END INTERFACE


  error = xfsetxsectgeomz_f (a_Id, a_index, a_NumVals, a_iValues, a_compression)

  return

END SUBROUTINE ! XF_SET_XSECT_GEOM_Z

!******************************************************************************
!* FUNCTION  XF_GET_XSECT_GEOM_Z
!* PURPOSE   Get the Cross Section Geometry Z
!* NOTES  
!************************************************************************************

SUBROUTINE XF_GET_XSECT_GEOM_Z (a_Id, a_index, a_NumVals, a_iValues, error)
INTEGER,          INTENT(IN)  :: a_Id
INTEGER,          INTENT(IN)  :: a_NumVals, a_index
REAL*4,             INTENT(IN)  :: a_iValues
INTEGER,          INTENT(OUT) :: error
  

!            INTEGER, EXTERNAL :: _xfgetxsectgeomz_f
!  MS FORTRAN needs explicit interface for C functions called here.

  INTERFACE
    INTEGER FUNCTION xfgetxsectgeomz_f (a_Id, a_index, a_NumVals, a_iValues)
    
	  !MS$ATTRIBUTES C,reference::xfgetxsectgeomz_f
    INTEGER,          INTENT(IN)  :: a_Id
    INTEGER,          INTENT(IN)  :: a_NumVals, a_index
    REAL*4,             INTENT(IN)  :: a_iValues


    END FUNCTION xfgetxsectgeomz_f
  END INTERFACE


  error = xfgetxsectgeomz_f (a_Id, a_index, a_NumVals, a_iValues)

  return

END SUBROUTINE ! XF_GET_XSECT_GEOM_Z


!******************************************************************************
!* FUNCTION  XF_SET_XSECT_LINE_FROM
!* PURPOSE   Set the Cross Section Line Property From
!* NOTES  
!************************************************************************************

SUBROUTINE XF_SET_XSECT_LINE_FROM (a_Id, a_index, a_NumVals, a_PropId, a_compression, error)
INTEGER,          INTENT(IN)  :: a_Id
INTEGER,          INTENT(IN)  :: a_NumVals, a_index, a_compression
REAL*4,             INTENT(IN)  :: a_PropId
INTEGER,          INTENT(OUT) :: error
  

!            INTEGER, EXTERNAL :: _xfsetxsectlinepropfrom_f
!  MS FORTRAN needs explicit interface for C functions called here.

  INTERFACE
    INTEGER FUNCTION xfsetxsectlinepropfrom_f (a_Id, a_index, a_NumVals, a_PropId, a_compression)
    
	  !MS$ATTRIBUTES C,reference::xfsetxsectlinepropfrom_f
    INTEGER,          INTENT(IN)  :: a_Id
    INTEGER,          INTENT(IN)  :: a_NumVals, a_index, a_compression
    REAL*4,             INTENT(IN)  :: a_PropId


    END FUNCTION xfsetxsectlinepropfrom_f
  END INTERFACE


  error = xfsetxsectlinepropfrom_f (a_Id, a_index, a_NumVals, a_PropId, a_compression)

  return

END SUBROUTINE ! XF_SET_XSECT_LINE_FROM

!******************************************************************************
!* FUNCTION  XF_GET_XSECT_LINE_FROM
!* PURPOSE   Get the Cross Section Line Property From
!* NOTES  
!************************************************************************************

SUBROUTINE XF_GET_XSECT_LINE_FROM (a_Id, a_index, a_NumVals, a_PropId, error)
INTEGER,          INTENT(IN)  :: a_Id
INTEGER,          INTENT(IN)  :: a_NumVals, a_index
REAL*4,             INTENT(IN)  :: a_PropId
INTEGER,          INTENT(OUT) :: error
  

!            INTEGER, EXTERNAL :: _xfgetxsectlinepropfrom_f
!  MS FORTRAN needs explicit interface for C functions called here.

  INTERFACE
    INTEGER FUNCTION xfgetxsectlinepropfrom_f (a_Id, a_index, a_NumVals, a_PropId)
    
	  !MS$ATTRIBUTES C,reference::xfgetxsectlinepropfrom_f
    INTEGER,          INTENT(IN)  :: a_Id
    INTEGER,          INTENT(IN)  :: a_NumVals, a_index
    REAL*4,             INTENT(IN)  :: a_PropId


    END FUNCTION xfgetxsectlinepropfrom_f
  END INTERFACE


  error = xfgetxsectlinepropfrom_f (a_Id, a_index, a_NumVals, a_PropId)

  return

END SUBROUTINE ! XF_GET_XSECT_LINE_FROM


!******************************************************************************
!* FUNCTION  XF_SET_XSECT_LINE_TO
!* PURPOSE   Set the Cross Section Line Property To
!* NOTES  
!************************************************************************************

SUBROUTINE XF_SET_XSECT_LINE_TO (a_Id, a_index, a_NumVals, a_PropId, a_compression, error)
INTEGER,          INTENT(IN)  :: a_Id
INTEGER,          INTENT(IN)  :: a_NumVals, a_index, a_compression
REAL*4,             INTENT(IN)  :: a_PropId
INTEGER,          INTENT(OUT) :: error
  

!            INTEGER, EXTERNAL :: _xfsetxsectlinepropto_f
!  MS FORTRAN needs explicit interface for C functions called here.

  INTERFACE
    INTEGER FUNCTION xfsetxsectlinepropto_f (a_Id, a_index, a_NumVals, a_PropId, a_compression)
    
	  !MS$ATTRIBUTES C,reference::xfsetxsectlinepropto_f
    INTEGER,          INTENT(IN)  :: a_Id
    INTEGER,          INTENT(IN)  :: a_NumVals, a_index, a_compression
    REAL*4,             INTENT(IN)  :: a_PropId


    END FUNCTION xfsetxsectlinepropto_f
  END INTERFACE


  error = xfsetxsectlinepropto_f (a_Id, a_index, a_NumVals, a_PropId, a_compression)

  return

END SUBROUTINE ! XF_SET_XSECT_LINE_TO

!******************************************************************************
!* FUNCTION  XF_GET_XSECT_LINE_TO
!* PURPOSE   Get the Cross Section Line Property To
!* NOTES  
!************************************************************************************

SUBROUTINE XF_GET_XSECT_LINE_TO (a_Id, a_index, a_NumVals, a_PropId, error)
INTEGER,          INTENT(IN)  :: a_Id
INTEGER,          INTENT(IN)  :: a_NumVals, a_index
REAL*4,             INTENT(IN)  :: a_PropId
INTEGER,          INTENT(OUT) :: error
  

!            INTEGER, EXTERNAL :: _xfgetxsectlinepropto_f
!  MS FORTRAN needs explicit interface for C functions called here.

  INTERFACE
    INTEGER FUNCTION xfgetxsectlinepropto_f (a_Id, a_index, a_NumVals, a_PropId)
    
	  !MS$ATTRIBUTES C,reference::xfgetxsectlinepropto_f
    INTEGER,          INTENT(IN)  :: a_Id
    INTEGER,          INTENT(IN)  :: a_NumVals, a_index
    REAL*4,             INTENT(IN)  :: a_PropId


    END FUNCTION xfgetxsectlinepropto_f
  END INTERFACE


  error = xfgetxsectlinepropto_f (a_Id, a_index, a_NumVals, a_PropId)

  return

END SUBROUTINE ! XF_GET_XSECT_LINE_TO

!******************************************************************************
!* FUNCTION  XF_SET_XSECT_LINE_TYPE
!* PURPOSE   Set the Cross Section Line Property Type
!* NOTES  
!************************************************************************************

SUBROUTINE XF_SET_XSECT_LINE_TYPE (a_Id, a_index, a_NumVals, a_PropId, a_compression, error)
INTEGER,          INTENT(IN)  :: a_Id
INTEGER,          INTENT(IN)  :: a_NumVals, a_index, a_compression
REAL*4,             INTENT(IN)  :: a_PropId
INTEGER,          INTENT(OUT) :: error
  

!            INTEGER, EXTERNAL :: _xfsetxsectlineproptype_f
!  MS FORTRAN needs explicit interface for C functions called here.

  INTERFACE
    INTEGER FUNCTION xfsetxsectlineproptype_f (a_Id, a_index, a_NumVals, a_PropId, a_compression)
    
	  !MS$ATTRIBUTES C,reference::xfsetxsectlineproptype_f
    INTEGER,          INTENT(IN)  :: a_Id
    INTEGER,          INTENT(IN)  :: a_NumVals, a_index, a_compression
    REAL*4,             INTENT(IN)  :: a_PropId


    END FUNCTION xfsetxsectlineproptype_f
  END INTERFACE


  error = xfsetxsectlineproptype_f (a_Id, a_index, a_NumVals, a_PropId, a_compression)

  return

END SUBROUTINE ! XF_SET_XSECT_LINE_TYPE

!******************************************************************************
!* FUNCTION  XF_GET_XSECT_LINE_TYPE
!* PURPOSE   Get the Cross Section Line Property Type
!* NOTES  
!************************************************************************************

SUBROUTINE XF_GET_XSECT_LINE_TYPE (a_Id, a_index, a_NumVals, a_PropId, error)
INTEGER,          INTENT(IN)  :: a_Id
INTEGER,          INTENT(IN)  :: a_NumVals, a_index
REAL*4,             INTENT(IN)  :: a_PropId
INTEGER,          INTENT(OUT) :: error
  

!            INTEGER, EXTERNAL :: _xfgetxsectlineproptype_f
!  MS FORTRAN needs explicit interface for C functions called here.

  INTERFACE
    INTEGER FUNCTION xfgetxsectlineproptype_f (a_Id, a_index, a_NumVals, a_PropId)
    
	  !MS$ATTRIBUTES C,reference::xfgetxsectlineproptype_f
    INTEGER,          INTENT(IN)  :: a_Id
    INTEGER,          INTENT(IN)  :: a_NumVals, a_index
    REAL*4,             INTENT(IN)  :: a_PropId


    END FUNCTION xfgetxsectlineproptype_f
  END INTERFACE


  error = xfgetxsectlineproptype_f (a_Id, a_index, a_NumVals, a_PropId)

  return

END SUBROUTINE ! XF_GET_XSECT_LINE_TYPE


!******************************************************************************
!* FUNCTION  XF_SET_XSECT_LINE_IVALUE
!* PURPOSE   Set the Cross Section Line Property I Value
!* NOTES  
!************************************************************************************

SUBROUTINE XF_SET_XSECT_LINE_IVALUE (a_Id, a_index, a_NumVals, a_PropId, a_compression, error)
INTEGER,          INTENT(IN)  :: a_Id
INTEGER,          INTENT(IN)  :: a_NumVals, a_index, a_compression
REAL*4,             INTENT(IN)  :: a_PropId
INTEGER,          INTENT(OUT) :: error
  

!            INTEGER, EXTERNAL :: _xfsetxsectlinepropivalue_f
!  MS FORTRAN needs explicit interface for C functions called here.

  INTERFACE
    INTEGER FUNCTION xfsetxsectlinepropivalue_f (a_Id, a_index, a_NumVals, a_PropId, a_compression)
    
	  !MS$ATTRIBUTES C,reference::xfsetxsectlinepropivalue_f
    INTEGER,          INTENT(IN)  :: a_Id
    INTEGER,          INTENT(IN)  :: a_NumVals, a_index, a_compression
    REAL*4,             INTENT(IN)  :: a_PropId


    END FUNCTION xfsetxsectlinepropivalue_f
  END INTERFACE


  error = xfsetxsectlinepropivalue_f (a_Id, a_index, a_NumVals, a_PropId, a_compression)

  return

END SUBROUTINE ! XF_SET_XSECT_LINE_IVALUE

!******************************************************************************
!* FUNCTION  XF_GET_XSECT_LINE_IVALUE
!* PURPOSE   Get the Cross Section Line Property I Value
!* NOTES  
!************************************************************************************

SUBROUTINE XF_GET_XSECT_LINE_IVALUE (a_Id, a_index, a_NumVals, a_PropId, error)
INTEGER,          INTENT(IN)  :: a_Id
INTEGER,          INTENT(IN)  :: a_NumVals, a_index
REAL*4,             INTENT(IN)  :: a_PropId
INTEGER,          INTENT(OUT) :: error
  

!            INTEGER, EXTERNAL :: _xfgetxsectlinepropivalue_f
!  MS FORTRAN needs explicit interface for C functions called here.

  INTERFACE
    INTEGER FUNCTION xfgetxsectlinepropivalue_f (a_Id, a_index, a_NumVals, a_PropId)
    
	  !MS$ATTRIBUTES C,reference::xfgetxsectlinepropivalue_f
    INTEGER,          INTENT(IN)  :: a_Id
    INTEGER,          INTENT(IN)  :: a_NumVals, a_index
    REAL*4,             INTENT(IN)  :: a_PropId


    END FUNCTION xfgetxsectlinepropivalue_f
  END INTERFACE


  error = xfgetxsectlinepropivalue_f (a_Id, a_index, a_NumVals, a_PropId)

  return

END SUBROUTINE ! XF_GET_XSECT_LINE_IVALUE


!******************************************************************************
!* FUNCTION  XF_SET_XSECT_LINE_FVALUE
!* PURPOSE   Set the Cross Section Line Property F Value
!* NOTES  
!************************************************************************************

SUBROUTINE XF_SET_XSECT_LINE_FVALUE (a_Id, a_index, a_NumVals, a_PropId, a_compression, error)
INTEGER,          INTENT(IN)  :: a_Id
INTEGER,          INTENT(IN)  :: a_NumVals, a_index, a_compression
REAL*4,             INTENT(IN)  :: a_PropId
INTEGER,          INTENT(OUT) :: error
  

!            INTEGER, EXTERNAL :: _xfsetxsectlinepropfvalue_f
!  MS FORTRAN needs explicit interface for C functions called here.

  INTERFACE
    INTEGER FUNCTION xfsetxsectlinepropfvalue_f (a_Id, a_index, a_NumVals, a_PropId, a_compression)
    
	  !MS$ATTRIBUTES C,reference::xfsetxsectlinepropfvalue_f
    INTEGER,          INTENT(IN)  :: a_Id
    INTEGER,          INTENT(IN)  :: a_NumVals, a_index, a_compression
    REAL*4,             INTENT(IN)  :: a_PropId


    END FUNCTION xfsetxsectlinepropfvalue_f
  END INTERFACE


  error = xfsetxsectlinepropfvalue_f (a_Id, a_index, a_NumVals, a_PropId, a_compression)

  return

END SUBROUTINE ! XF_SET_XSECT_LINE_FVALUE

!******************************************************************************
!* FUNCTION  XF_GET_XSECT_LINE_FVALUE
!* PURPOSE   Get the Cross Section Line Property F Value
!* NOTES  
!************************************************************************************

SUBROUTINE XF_GET_XSECT_LINE_FVALUE (a_Id, a_index, a_NumVals, a_PropId, error)
INTEGER,          INTENT(IN)  :: a_Id
INTEGER,          INTENT(IN)  :: a_NumVals, a_index
REAL*4,             INTENT(IN)  :: a_PropId
INTEGER,          INTENT(OUT) :: error
  

!            INTEGER, EXTERNAL :: _xfgetxsectlinepropfvalue_f
!  MS FORTRAN needs explicit interface for C functions called here.

  INTERFACE
    INTEGER FUNCTION xfgetxsectlinepropfvalue_f (a_Id, a_index, a_NumVals, a_PropId)
    
	  !MS$ATTRIBUTES C,reference::xfgetxsectlinepropfvalue_f
    INTEGER,          INTENT(IN)  :: a_Id
    INTEGER,          INTENT(IN)  :: a_NumVals, a_index
    REAL*4,             INTENT(IN)  :: a_PropId


    END FUNCTION xfgetxsectlinepropfvalue_f
  END INTERFACE


  error = xfgetxsectlinepropfvalue_f (a_Id, a_index, a_NumVals, a_PropId)

  return

END SUBROUTINE ! XF_GET_XSECT_LINE_FVALUE


!******************************************************************************
!* FUNCTION  XF_SET_XSECT_LINE_ID
!* PURPOSE   Set the Cross Section Line Property ID
!* NOTES  
!************************************************************************************

SUBROUTINE XF_SET_XSECT_LINE_ID (a_Id, a_NumVals, a_PropId, a_compression, error)
INTEGER,          INTENT(IN)  :: a_Id
INTEGER,          INTENT(IN)  :: a_NumVals, a_compression, a_PropId
INTEGER,          INTENT(OUT) :: error
  

!            INTEGER, EXTERNAL :: _xfsetxsectlinepropid_f
!  MS FORTRAN needs explicit interface for C functions called here.

  INTERFACE
    INTEGER FUNCTION xfsetxsectlinepropid_f (a_Id, a_NumVals, a_PropId, a_compression)
    
	  !MS$ATTRIBUTES C,reference::xfsetxsectlinepropid_f
    INTEGER,          INTENT(IN)  :: a_Id
    INTEGER,          INTENT(IN)  :: a_NumVals, a_compression, a_PropId


    END FUNCTION xfsetxsectlinepropid_f
  END INTERFACE


  error = xfsetxsectlinepropid_f (a_Id, a_NumVals, a_PropId, a_compression)

  return

END SUBROUTINE ! XF_SET_XSECT_LINE_ID

!******************************************************************************
!* FUNCTION  XF_GET_XSECT_LINE_ID
!* PURPOSE   Get the Cross Section Line Property ID
!* NOTES  
!************************************************************************************

SUBROUTINE XF_GET_XSECT_LINE_ID (a_Id, a_NumVals, a_PropId, error)
INTEGER,          INTENT(IN)  :: a_Id
INTEGER,          INTENT(IN)  :: a_NumVals, a_PropId
INTEGER,          INTENT(OUT) :: error
  

!            INTEGER, EXTERNAL :: _xfgetxsectlinepropid_f
!  MS FORTRAN needs explicit interface for C functions called here.

  INTERFACE
    INTEGER FUNCTION xfgetxsectlinepropid_f (a_Id, a_NumVals, a_PropId)
    
	  !MS$ATTRIBUTES C,reference::xfgetxsectlinepropid_f
    INTEGER,          INTENT(IN)  :: a_Id
    INTEGER,          INTENT(IN)  :: a_NumVals, a_PropId


    END FUNCTION xfgetxsectlinepropid_f
  END INTERFACE


  error = xfgetxsectlinepropid_f (a_Id, a_NumVals, a_PropId)

  return

END SUBROUTINE ! XF_GET_XSECT_LINE_ID


!******************************************************************************
!* FUNCTION  XF_SET_XSECT_LINE_NAME
!* PURPOSE   Set the Cross Sections Line Properties Name
!* NOTES  
!************************************************************************************

SUBROUTINE XF_SET_XSECT_LINE_NAME (a_Id, a_NumVals, a_PropId, a_StrLen, a_compression, error)
INTEGER,          INTENT(IN)  :: a_Id
INTEGER,          INTENT(IN)  :: a_NumVals, a_StrLen, a_compression
CHARACTER(LEN=*), INTENT(IN)  :: a_PropId
INTEGER,          INTENT(OUT) :: error
  

!            INTEGER, EXTERNAL :: _xfsetxsectlinepropname_f
!  MS FORTRAN needs explicit interface for C functions called here.

  INTERFACE
    INTEGER FUNCTION xfsetxsectlinepropname_f (a_Id, a_NumVals, a_PropId, a_StrLen, a_compression)
    
	  !MS$ATTRIBUTES C,reference::xfsetxsectlinepropname_f
    !DEC$ATTRIBUTES reference :: a_PropId
    INTEGER,          INTENT(IN)  :: a_Id
    INTEGER,          INTENT(IN)  :: a_NumVals, a_StrLen, a_compression
    CHARACTER(LEN=*), INTENT(IN)  :: a_PropId


    END FUNCTION xfsetxsectlinepropname_f
  END INTERFACE


  error = xfsetxsectlinepropname_f (a_Id, a_NumVals, a_PropId, a_StrLen, a_compression)

  return

END SUBROUTINE ! XF_SET_XSECT_LINE_NAME

!******************************************************************************
!* FUNCTION  XF_GET_XSECT_LINE_NAME
!* PURPOSE   Get the Cross Sections Line Properties Name
!* NOTES  
!************************************************************************************

SUBROUTINE XF_GET_XSECT_LINE_NAME (a_Id, a_NumVals, a_PropId, a_StrLen, error)
INTEGER,          INTENT(IN)  :: a_Id
INTEGER,          INTENT(IN)  :: a_NumVals, a_StrLen
CHARACTER(LEN=*), INTENT(IN)  :: a_PropId
INTEGER,          INTENT(OUT) :: error
  

!            INTEGER, EXTERNAL :: _xfgetxsectlinepropname_f
!  MS FORTRAN needs explicit interface for C functions called here.

  INTERFACE
    INTEGER FUNCTION xfgetxsectlinepropname_f (a_Id, a_NumVals, a_PropId, a_StrLen)
    
	  !MS$ATTRIBUTES C,reference::xfgetxsectlinepropname_f
    !DEC$ATTRIBUTES reference :: a_PropId
    INTEGER,          INTENT(IN)  :: a_Id
    INTEGER,          INTENT(IN)  :: a_NumVals, a_StrLen
    CHARACTER(LEN=*), INTENT(IN)  :: a_PropId


    END FUNCTION xfgetxsectlinepropname_f
  END INTERFACE


  error = xfgetxsectlinepropname_f (a_Id, a_NumVals, a_PropId, a_StrLen)

  return

END SUBROUTINE ! XF_GET_XSECT_LINE_NAME


!******************************************************************************
!* FUNCTION  XF_GET_XSECT_LINE_NAME_LEN
!* PURPOSE   Get the Cross Sections Line Properties Name Length
!* NOTES  
!************************************************************************************

SUBROUTINE XF_GET_XSECT_LINE_NAME_LEN (a_Id, a_NumVals, a_StrLen, error)
INTEGER,          INTENT(IN)  :: a_Id
INTEGER,          INTENT(IN)  :: a_NumVals, a_StrLen
INTEGER,          INTENT(OUT) :: error
  

!            INTEGER, EXTERNAL :: _xfgetxsectlinepropnamelen_f
!  MS FORTRAN needs explicit interface for C functions called here.

  INTERFACE
    INTEGER FUNCTION xfgetxsectlinepropnamelen_f (a_Id, a_NumVals, a_StrLen)
    
	  !MS$ATTRIBUTES C,reference::xfgetxsectlinepropnamelen_f
    INTEGER,          INTENT(IN)  :: a_Id
    INTEGER,          INTENT(IN)  :: a_NumVals, a_StrLen


    END FUNCTION xfgetxsectlinepropnamelen_f
  END INTERFACE


  error = xfgetxsectlinepropnamelen_f (a_Id, a_NumVals, a_StrLen)

  return

END SUBROUTINE ! XF_GET_XSECT_LINE_NAME_LEN


!******************************************************************************
!* FUNCTION  XF_SET_XSECT_LINE_DESC
!* PURPOSE   Set the Cross Sections Line Properties Description
!* NOTES  
!************************************************************************************

SUBROUTINE XF_SET_XSECT_LINE_DESC (a_Id, a_NumVals, a_PropId, a_StrLen, a_compression, error)
INTEGER,          INTENT(IN)  :: a_Id
INTEGER,          INTENT(IN)  :: a_NumVals, a_StrLen, a_compression
CHARACTER(LEN=*), INTENT(IN)  :: a_PropId
INTEGER,          INTENT(OUT) :: error
  

!            INTEGER, EXTERNAL :: _xfsetxsectlinepropdesc_f
!  MS FORTRAN needs explicit interface for C functions called here.

  INTERFACE
    INTEGER FUNCTION xfsetxsectlinepropdesc_f (a_Id, a_NumVals, a_PropId, a_StrLen, a_compression)
    
	  !MS$ATTRIBUTES C,reference::xfsetxsectlinepropdesc_f
    !DEC$ATTRIBUTES reference :: a_PropId
    INTEGER,          INTENT(IN)  :: a_Id
    INTEGER,          INTENT(IN)  :: a_NumVals, a_StrLen, a_compression
    CHARACTER(LEN=*), INTENT(IN)  :: a_PropId


    END FUNCTION xfsetxsectlinepropdesc_f
  END INTERFACE


  error = xfsetxsectlinepropdesc_f (a_Id, a_NumVals, a_PropId, a_StrLen, a_compression)

  return

END SUBROUTINE ! XF_SET_XSECT_LINE_DESC

!******************************************************************************
!* FUNCTION  XF_GET_XSECT_LINE_DESC
!* PURPOSE   Get the Cross Sections Line Properties Description
!* NOTES  
!************************************************************************************

SUBROUTINE XF_GET_XSECT_LINE_DESC (a_Id, a_NumVals, a_PropId, a_StrLen, error)
INTEGER,          INTENT(IN)  :: a_Id
INTEGER,          INTENT(IN)  :: a_NumVals, a_StrLen
CHARACTER(LEN=*), INTENT(IN)  :: a_PropId
INTEGER,          INTENT(OUT) :: error
  

!            INTEGER, EXTERNAL :: _xfgetxsectlinepropdesc_f
!  MS FORTRAN needs explicit interface for C functions called here.

  INTERFACE
    INTEGER FUNCTION xfgetxsectlinepropdesc_f (a_Id, a_NumVals, a_PropId, a_StrLen)
    
	  !MS$ATTRIBUTES C,reference::xfgetxsectlinepropdesc_f
    !DEC$ATTRIBUTES reference :: a_PropId
    INTEGER,          INTENT(IN)  :: a_Id
    INTEGER,          INTENT(IN)  :: a_NumVals, a_StrLen
    CHARACTER(LEN=*), INTENT(IN)  :: a_PropId


    END FUNCTION xfgetxsectlinepropdesc_f
  END INTERFACE


  error = xfgetxsectlinepropdesc_f (a_Id, a_NumVals, a_PropId, a_StrLen)

  return

END SUBROUTINE ! XF_GET_XSECT_LINE_DESC


!******************************************************************************
!* FUNCTION  XF_GET_XSECT_LINE_DESC_LEN
!* PURPOSE   Get the Cross Sections Description Length
!* NOTES  
!************************************************************************************

SUBROUTINE XF_GET_XSECT_LINE_DESC_LEN (a_Id, a_NumVals, a_StrLen, error)
INTEGER,          INTENT(IN)  :: a_Id
INTEGER,          INTENT(IN)  :: a_NumVals, a_StrLen
INTEGER,          INTENT(OUT) :: error
  

!            INTEGER, EXTERNAL :: _xfgetxsectlinepropdesclen_f
!  MS FORTRAN needs explicit interface for C functions called here.

  INTERFACE
    INTEGER FUNCTION xfgetxsectlinepropdesclen_f (a_Id, a_NumVals, a_StrLen)
    
	  !MS$ATTRIBUTES C,reference::xfgetxsectlinepropdesclen_f
    INTEGER,          INTENT(IN)  :: a_Id
    INTEGER,          INTENT(IN)  :: a_NumVals, a_StrLen


    END FUNCTION xfgetxsectlinepropdesclen_f
  END INTERFACE


  error = xfgetxsectlinepropdesclen_f (a_Id, a_NumVals, a_StrLen)

  return

END SUBROUTINE ! XF_GET_XSECT_LINE_DESC_LEN


!******************************************************************************
!* FUNCTION  XF_SET_XSECT_LINE_CATEGORY
!* PURPOSE   Set the Cross Section Line Property Category
!* NOTES  
!************************************************************************************

SUBROUTINE XF_SET_XSECT_LINE_CATEGORY (a_Id, a_NumVals, a_PropId, a_compression, error)
INTEGER,          INTENT(IN)  :: a_Id
INTEGER,          INTENT(IN)  :: a_NumVals, a_compression, a_PropId
INTEGER,          INTENT(OUT) :: error
  

!            INTEGER, EXTERNAL :: _xfsetxsectlinepropcategory_f
!  MS FORTRAN needs explicit interface for C functions called here.

  INTERFACE
    INTEGER FUNCTION xfsetxsectlinepropcategory_f (a_Id, a_NumVals, a_PropId, a_compression)
    
	  !MS$ATTRIBUTES C,reference::xfsetxsectlinepropcategory_f
    INTEGER,          INTENT(IN)  :: a_Id
    INTEGER,          INTENT(IN)  :: a_NumVals, a_compression, a_PropId


    END FUNCTION xfsetxsectlinepropcategory_f
  END INTERFACE


  error = xfsetxsectlinepropcategory_f (a_Id, a_NumVals, a_PropId, a_compression)

  return

END SUBROUTINE ! XF_SET_XSECT_LINE_CATEGORY

!******************************************************************************
!* FUNCTION  XF_GET_XSECT_LINE_CATEGORY
!* PURPOSE   Get the Cross Section Line Property Category
!* NOTES  
!************************************************************************************

SUBROUTINE XF_GET_XSECT_LINE_CATEGORY (a_Id, a_NumVals, a_PropId, error)
INTEGER,          INTENT(IN)  :: a_Id
INTEGER,          INTENT(IN)  :: a_NumVals, a_PropId
INTEGER,          INTENT(OUT) :: error
  

!            INTEGER, EXTERNAL :: _xfgetxsectlinepropcategory_f
!  MS FORTRAN needs explicit interface for C functions called here.

  INTERFACE
    INTEGER FUNCTION xfgetxsectlinepropcategory_f (a_Id, a_NumVals, a_PropId)
    
	  !MS$ATTRIBUTES C,reference::xfgetxsectlinepropcategory_f
    INTEGER,          INTENT(IN)  :: a_Id
    INTEGER,          INTENT(IN)  :: a_NumVals, a_PropId


    END FUNCTION xfgetxsectlinepropcategory_f
  END INTERFACE


  error = xfgetxsectlinepropcategory_f (a_Id, a_NumVals, a_PropId)

  return

END SUBROUTINE ! XF_GET_XSECT_LINE_CATEGORY


!******************************************************************************
!* FUNCTION  XF_SET_XSECT_LINE_DEFAULT
!* PURPOSE   Set the Cross Section Line Property Default
!* NOTES  
!************************************************************************************

SUBROUTINE XF_SET_XSECT_LINE_DEFAULT (a_Id, a_NumVals, a_PropId, a_compression, error)
INTEGER,          INTENT(IN)  :: a_Id
INTEGER,          INTENT(IN)  :: a_NumVals, a_compression
REAL*4,             INTENT(IN)  :: a_PropId
INTEGER,          INTENT(OUT) :: error
  

!            INTEGER, EXTERNAL :: _xfsetxsectlinepropfdefault_f
!  MS FORTRAN needs explicit interface for C functions called here.

  INTERFACE
    INTEGER FUNCTION xfsetxsectlinepropfdefault_f (a_Id, a_NumVals, a_PropId, a_compression)
    
	  !MS$ATTRIBUTES C,reference::xfsetxsectlinepropfdefault_f
    INTEGER,          INTENT(IN)  :: a_Id
    INTEGER,          INTENT(IN)  :: a_NumVals, a_compression
    REAL*4,             INTENT(IN)  :: a_PropId


    END FUNCTION xfsetxsectlinepropfdefault_f
  END INTERFACE


  error = xfsetxsectlinepropfdefault_f (a_Id, a_NumVals, a_PropId, a_compression)

  return

END SUBROUTINE ! XF_SET_XSECT_LINE_DEFAULT

!******************************************************************************
!* FUNCTION  XF_GET_XSECT_LINE_DEFAULT
!* PURPOSE   Get the Cross Section Line Property Default
!* NOTES  
!************************************************************************************

SUBROUTINE XF_GET_XSECT_LINE_DEFAULT (a_Id, a_NumVals, a_PropId, error)
INTEGER,          INTENT(IN)  :: a_Id
INTEGER,          INTENT(IN)  :: a_NumVals
REAL*4,             INTENT(IN)  :: a_PropId
INTEGER,          INTENT(OUT) :: error
  

!            INTEGER, EXTERNAL :: _xfgetxsectlinepropfdefault_f
!  MS FORTRAN needs explicit interface for C functions called here.

  INTERFACE
    INTEGER FUNCTION xfgetxsectlinepropfdefault_f (a_Id, a_NumVals, a_PropId)
    
	  !MS$ATTRIBUTES C,reference::xfgetxsectlinepropfdefault_f
    INTEGER,          INTENT(IN)  :: a_Id
    INTEGER,          INTENT(IN)  :: a_NumVals
    REAL*4,             INTENT(IN)  :: a_PropId


    END FUNCTION xfgetxsectlinepropfdefault_f
  END INTERFACE


  error = xfgetxsectlinepropfdefault_f (a_Id, a_NumVals, a_PropId)

  return

END SUBROUTINE ! XF_GET_XSECT_LINE_DEFAULT


!******************************************************************************
!* FUNCTION  XF_SET_XSECT_LINE_EXCLUSIVE
!* PURPOSE   Set the Cross Section Line Property Exclusive
!* NOTES  
!************************************************************************************

SUBROUTINE XF_SET_XSECT_LINE_EXCLUSIVE (a_Id, a_NumVals, a_PropId, a_compression, error)
INTEGER,          INTENT(IN)  :: a_Id
INTEGER,          INTENT(IN)  :: a_NumVals, a_compression, a_PropId
INTEGER,          INTENT(OUT) :: error
  

!            INTEGER, EXTERNAL :: _xfsetxsectlinepropexclusive_f
!  MS FORTRAN needs explicit interface for C functions called here.

  INTERFACE
    INTEGER FUNCTION xfsetxsectlinepropexclusive_f (a_Id, a_NumVals, a_PropId, a_compression)
    
	  !MS$ATTRIBUTES C,reference::xfsetxsectlinepropexclusive_f
    INTEGER,          INTENT(IN)  :: a_Id
    INTEGER,          INTENT(IN)  :: a_NumVals, a_compression, a_PropId


    END FUNCTION xfsetxsectlinepropexclusive_f
  END INTERFACE


  error = xfsetxsectlinepropexclusive_f (a_Id, a_NumVals, a_PropId, a_compression)

  return

END SUBROUTINE ! XF_SET_XSECT_LINE_EXCLUSIVE

!******************************************************************************
!* FUNCTION  XF_GET_XSECT_LINE_EXCLUSIVE
!* PURPOSE   Get the Cross Section Line Property Exclusive
!* NOTES  
!************************************************************************************

SUBROUTINE XF_GET_XSECT_LINE_EXCLUSIVE (a_Id, a_NumVals, a_PropId, error)
INTEGER,          INTENT(IN)  :: a_Id
INTEGER,          INTENT(IN)  :: a_NumVals, a_PropId
INTEGER,          INTENT(OUT) :: error
  

!            INTEGER, EXTERNAL :: _xfgetxsectlinepropexclusive_f
!  MS FORTRAN needs explicit interface for C functions called here.

  INTERFACE
    INTEGER FUNCTION xfgetxsectlinepropexclusive_f (a_Id, a_NumVals, a_PropId)
    
	  !MS$ATTRIBUTES C,reference::xfgetxsectlinepropexclusive_f
    INTEGER,          INTENT(IN)  :: a_Id
    INTEGER,          INTENT(IN)  :: a_NumVals, a_PropId


    END FUNCTION xfgetxsectlinepropexclusive_f
  END INTERFACE


  error = xfgetxsectlinepropexclusive_f (a_Id, a_NumVals, a_PropId)

  return

END SUBROUTINE ! XF_GET_XSECT_LINE_EXCLUSIVE


!******************************************************************************
!* FUNCTION  XF_SET_XSECT_LINE_ENUM_GROUP
!* PURPOSE   Set the Cross Section Line Property Enumeration Group Number
!* NOTES  
!************************************************************************************

SUBROUTINE XF_SET_XSECT_LINE_ENUM_GROUP (a_Id, a_PropId, a_compression, error)
INTEGER,          INTENT(IN)  :: a_Id
INTEGER,          INTENT(IN)  :: a_compression, a_PropId
INTEGER,          INTENT(OUT) :: error
  

!            INTEGER, EXTERNAL :: _xfSetNumberOfLnPropEnumGroup_f
!  MS FORTRAN needs explicit interface for C functions called here.

  INTERFACE
    INTEGER FUNCTION xfsetnumberoflnpropenumgroup_f (a_Id, a_PropId, a_compression)
    
	  !MS$ATTRIBUTES C,reference::xfsetnumberoflnpropenumgroup_f
    INTEGER,          INTENT(IN)  :: a_Id
    INTEGER,          INTENT(IN)  :: a_compression, a_PropId


    END FUNCTION xfsetnumberoflnpropenumgroup_f
  END INTERFACE


  error = xfsetnumberoflnpropenumgroup_f (a_Id, a_PropId, a_compression)

  return

END SUBROUTINE ! XF_SET_XSECT_LINE_ENUM_GROUP

!******************************************************************************
!* FUNCTION  XF_GET_XSECT_LINE_ENUM_GROUP
!* PURPOSE   Get the Cross Section Line Property Enumeration Group Number
!* NOTES  
!************************************************************************************

SUBROUTINE XF_GET_XSECT_LINE_ENUM_GROUP (a_Id, a_PropId, error)
INTEGER,          INTENT(IN)  :: a_Id
INTEGER,          INTENT(IN)  :: a_PropId
INTEGER,          INTENT(OUT) :: error
  

!            INTEGER, EXTERNAL :: _xfgetnumberoflnpropenumgroup_f
!  MS FORTRAN needs explicit interface for C functions called here.

  INTERFACE
    INTEGER FUNCTION xfgetnumberoflnpropenumgroup_f (a_Id, a_PropId)
    
	  !MS$ATTRIBUTES C,reference::xfgetnumberoflnpropenumgroup_f
    INTEGER,          INTENT(IN)  :: a_Id
    INTEGER,          INTENT(IN)  :: a_PropId


    END FUNCTION xfgetnumberoflnpropenumgroup_f
  END INTERFACE


  error = xfgetnumberoflnpropenumgroup_f (a_Id, a_PropId)

  return

END SUBROUTINE ! XF_GET_XSECT_LINE_ENUM_GROUP


!******************************************************************************
!* FUNCTION  XF_SET_XSECT_LINE_ENUM
!* PURPOSE   Set the Cross Section Line Property Enumeration Number
!* NOTES  
!************************************************************************************

SUBROUTINE XF_SET_XSECT_LINE_ENUM (a_Id, a_PropId, a_compression, error)
INTEGER,          INTENT(IN)  :: a_Id
INTEGER,          INTENT(IN)  :: a_compression, a_PropId
INTEGER,          INTENT(OUT) :: error
  

!            INTEGER, EXTERNAL :: _xfsetnumberoflinepropenum
!  MS FORTRAN needs explicit interface for C functions called here.

  INTERFACE
    INTEGER FUNCTION xfsetnumberoflinepropenum_f (a_Id, a_PropId, a_compression)
    
	  !MS$ATTRIBUTES C,reference::xfsetnumberoflinepropenum_f
    INTEGER,          INTENT(IN)  :: a_Id
    INTEGER,          INTENT(IN)  :: a_compression, a_PropId


    END FUNCTION xfsetnumberoflinepropenum_f
  END INTERFACE


  error = xfsetnumberoflinepropenum_f (a_Id, a_PropId, a_compression)

  return

END SUBROUTINE ! XF_SET_XSECT_LINE_ENUM

!******************************************************************************
!* FUNCTION  XF_GET_XSECT_LINE_ENUM
!* PURPOSE   Get the Cross Section Line Property Enumeration Number
!* NOTES  
!************************************************************************************

SUBROUTINE XF_GET_XSECT_LINE_ENUM (a_Id, a_PropId, error)
INTEGER,          INTENT(IN)  :: a_Id
INTEGER,          INTENT(IN)  :: a_PropId
INTEGER,          INTENT(OUT) :: error
  

!            INTEGER, EXTERNAL :: _xfgetnumberoflinepropenum_f
!  MS FORTRAN needs explicit interface for C functions called here.

  INTERFACE
    INTEGER FUNCTION xfgetnumberoflinepropenum_f (a_Id, a_PropId)
    
	  !MS$ATTRIBUTES C,reference::xfgetnumberoflinepropenum_f
    INTEGER,          INTENT(IN)  :: a_Id
    INTEGER,          INTENT(IN)  :: a_PropId


    END FUNCTION xfgetnumberoflinepropenum_f
  END INTERFACE


  error = xfgetnumberoflinepropenum_f (a_Id, a_PropId)

  return

END SUBROUTINE ! XF_GET_XSECT_LINE_ENUM


!******************************************************************************
!* FUNCTION  XF_SET_XSECT_LINE_ENUM_ID
!* PURPOSE   Set the Cross Section Line Property Exclusive
!* NOTES  
!************************************************************************************

SUBROUTINE XF_SET_XSECT_LINE_ENUM_ID (a_Id, a_index, a_NumVals, a_PropId, a_compression, error)
INTEGER,          INTENT(IN)  :: a_Id
INTEGER,          INTENT(IN)  :: a_index, a_NumVals, a_compression, a_PropId
INTEGER,          INTENT(OUT) :: error
  

!            INTEGER, EXTERNAL :: _xfsetxsectlinepropenumid_f
!  MS FORTRAN needs explicit interface for C functions called here.

  INTERFACE
    INTEGER FUNCTION xfsetxsectlinepropenumid_f (a_Id, a_index, a_NumVals, a_PropId, a_compression)
    
	  !MS$ATTRIBUTES C,reference::xfsetxsectlinepropenumid_f
    INTEGER,          INTENT(IN)  :: a_Id
    INTEGER,          INTENT(IN)  :: a_index, a_NumVals, a_compression, a_PropId


    END FUNCTION xfsetxsectlinepropenumid_f
  END INTERFACE


  error = xfsetxsectlinepropenumid_f (a_Id, a_index, a_NumVals, a_PropId, a_compression)

  return

END SUBROUTINE ! XF_SET_XSECT_LINE_ENUM_ID

!******************************************************************************
!* FUNCTION  XF_GET_XSECT_LINE_ENUM_ID
!* PURPOSE   Get the Cross Section Line Property Enumeration ID
!* NOTES  
!************************************************************************************

SUBROUTINE XF_GET_XSECT_LINE_ENUM_ID (a_Id, a_index, a_NumVals, a_PropId, error)
INTEGER,          INTENT(IN)  :: a_Id
INTEGER,          INTENT(IN)  :: a_index, a_NumVals, a_PropId
INTEGER,          INTENT(OUT) :: error
  

!            INTEGER, EXTERNAL :: _xfgetxsectlinepropenumid_f
!  MS FORTRAN needs explicit interface for C functions called here.

  INTERFACE
    INTEGER FUNCTION xfgetxsectlinepropenumid_f (a_Id, a_index, a_NumVals, a_PropId)
    
	  !MS$ATTRIBUTES C,reference::xfgetxsectlinepropenumid_f
    INTEGER,          INTENT(IN)  :: a_Id
    INTEGER,          INTENT(IN)  :: a_index, a_NumVals, a_PropId


    END FUNCTION xfgetxsectlinepropenumid_f
  END INTERFACE


  error = xfgetxsectlinepropenumid_f (a_Id, a_index, a_NumVals, a_PropId)

  return

END SUBROUTINE ! XF_GET_XSECT_LINE_ENUM_ID


!******************************************************************************
!* FUNCTION  XF_SET_XSECT_LINE_ENUM_MAT_ID
!* PURPOSE   Set the Cross Section Line Property Enumeration Material ID
!* NOTES  
!************************************************************************************

SUBROUTINE XF_SET_XSECT_LINE_ENUM_MAT_ID (a_Id, a_index, a_NumVals, a_PropId, a_compression, error)
INTEGER,          INTENT(IN)  :: a_Id
INTEGER,          INTENT(IN)  :: a_index, a_NumVals, a_compression, a_PropId
INTEGER,          INTENT(OUT) :: error
  

!            INTEGER, EXTERNAL :: _xfsetxsectlinepropenummatid_f
!  MS FORTRAN needs explicit interface for C functions called here.

  INTERFACE
    INTEGER FUNCTION xfsetxsectlinepropenummatid_f (a_Id, a_index, a_NumVals, a_PropId, a_compression)
    
	  !MS$ATTRIBUTES C,reference::xfsetxsectlinepropenummatid_f
    INTEGER,          INTENT(IN)  :: a_Id
    INTEGER,          INTENT(IN)  :: a_index, a_NumVals, a_compression, a_PropId


    END FUNCTION xfsetxsectlinepropenummatid_f
  END INTERFACE


  error = xfsetxsectlinepropenummatid_f (a_Id, a_index, a_NumVals, a_PropId, a_compression)

  return

END SUBROUTINE ! XF_SET_XSECT_LINE_ENUM_MAT_ID

!******************************************************************************
!* FUNCTION  XF_GET_XSECT_LINE_ENUM_MAT_ID
!* PURPOSE   Get the Cross Section Line Property Enumeration Material ID
!* NOTES  
!************************************************************************************

SUBROUTINE XF_GET_XSECT_LINE_ENUM_MAT_ID (a_Id, a_index, a_NumVals, a_PropId, error)
INTEGER,          INTENT(IN)  :: a_Id
INTEGER,          INTENT(IN)  :: a_index, a_NumVals, a_PropId
INTEGER,          INTENT(OUT) :: error
  

!            INTEGER, EXTERNAL :: _xfgetxsectlinepropenummatid_f
!  MS FORTRAN needs explicit interface for C functions called here.

  INTERFACE
    INTEGER FUNCTION xfgetxsectlinepropenummatid_f (a_Id, a_index, a_NumVals, a_PropId)
    
	  !MS$ATTRIBUTES C,reference::xfgetxsectlinepropenummatid_f
    INTEGER,          INTENT(IN)  :: a_Id
    INTEGER,          INTENT(IN)  :: a_index, a_NumVals, a_PropId


    END FUNCTION xfgetxsectlinepropenummatid_f
  END INTERFACE


  error = xfgetxsectlinepropenummatid_f (a_Id, a_index, a_NumVals, a_PropId)

  return

END SUBROUTINE ! XF_GET_XSECT_LINE_ENUM_MAT_ID


!******************************************************************************
!* FUNCTION  XF_SET_XSECT_LINE_ENUM_NAME
!* PURPOSE   Set the Cross Sections Name
!* NOTES  
!************************************************************************************

SUBROUTINE XF_SET_XSECT_LINE_ENUM_NAME (a_Id, a_NumVals, a_PropId, a_StrLen, a_compression, error)
INTEGER,          INTENT(IN)  :: a_Id
INTEGER,          INTENT(IN)  :: a_NumVals, a_StrLen, a_compression
CHARACTER(LEN=*), INTENT(IN)  :: a_PropId
INTEGER,          INTENT(OUT) :: error
  

!            INTEGER, EXTERNAL :: _xfsetxsectlinepropenumname_f
!  MS FORTRAN needs explicit interface for C functions called here.

  INTERFACE
    INTEGER FUNCTION xfsetxsectlinepropenumname_f (a_Id, a_NumVals, a_PropId, a_StrLen, a_compression)
    
	  !MS$ATTRIBUTES C,reference::xfsetxsectlinepropenumname_f
    !DEC$ATTRIBUTES reference :: a_PropId
    INTEGER,          INTENT(IN)  :: a_Id
    INTEGER,          INTENT(IN)  :: a_NumVals, a_StrLen, a_compression
    CHARACTER(LEN=*), INTENT(IN)  :: a_PropId


    END FUNCTION xfsetxsectlinepropenumname_f
  END INTERFACE


  error = xfsetxsectlinepropenumname_f (a_Id, a_NumVals, a_PropId, a_StrLen, a_compression)

  return

END SUBROUTINE ! XF_SET_XSECT_LINE_ENUM_NAME

!******************************************************************************
!* FUNCTION  XF_GET_XSECT_LINE_ENUM_NAME
!* PURPOSE   Get the Cross Section Line Property Enumeration Name
!* NOTES  
!************************************************************************************

SUBROUTINE XF_GET_XSECT_LINE_ENUM_NAME (a_Id, a_NumVals, a_PropId, a_StrLen, error)
INTEGER,          INTENT(IN)  :: a_Id
INTEGER,          INTENT(IN)  :: a_NumVals, a_StrLen
CHARACTER(LEN=*), INTENT(IN)  :: a_PropId
INTEGER,          INTENT(OUT) :: error
  

!            INTEGER, EXTERNAL :: _xfgetxsectlinepropenumname_f
!  MS FORTRAN needs explicit interface for C functions called here.

  INTERFACE
    INTEGER FUNCTION xfgetxsectlinepropenumname_f (a_Id, a_NumVals, a_PropId, a_StrLen)
    
	  !MS$ATTRIBUTES C,reference::xfgetxsectlinepropenumname_f
    !DEC$ATTRIBUTES reference :: a_PropId
    INTEGER,          INTENT(IN)  :: a_Id
    INTEGER,          INTENT(IN)  :: a_NumVals, a_StrLen
    CHARACTER(LEN=*), INTENT(IN)  :: a_PropId


    END FUNCTION xfgetxsectlinepropenumname_f
  END INTERFACE


  error = xfgetxsectlinepropenumname_f (a_Id, a_NumVals, a_PropId, a_StrLen)

  return

END SUBROUTINE ! XF_GET_XSECT_LINE_ENUM_NAME


!******************************************************************************
!* FUNCTION  XF_GET_XSECT_LINE_ENUM_NAME_LEN
!* PURPOSE   Get the Cross Sections Name Length
!* NOTES  
!************************************************************************************

SUBROUTINE XF_GET_XSECT_LINE_ENUM_NAME_LEN (a_Id, a_NumVals, a_StrLen, error)
INTEGER,          INTENT(IN)  :: a_Id
INTEGER,          INTENT(IN)  :: a_NumVals, a_StrLen
INTEGER,          INTENT(OUT) :: error
  

!            INTEGER, EXTERNAL :: _xfgetxsectlnpropenumnamelen_f
!  MS FORTRAN needs explicit interface for C functions called here.

  INTERFACE
    INTEGER FUNCTION xfgetxsectlnpropenumnamelen_f (a_Id, a_NumVals, a_StrLen)
    
	  !MS$ATTRIBUTES C,reference::xfgetxsectlnpropenumnamelen_f
    INTEGER,          INTENT(IN)  :: a_Id
    INTEGER,          INTENT(IN)  :: a_NumVals, a_StrLen


    END FUNCTION xfgetxsectlnpropenumnamelen_f
  END INTERFACE


  error = xfgetxsectlnpropenumnamelen_f (a_Id, a_NumVals, a_StrLen)

  return

END SUBROUTINE ! XF_GET_XSECT_LINE_ENUM_NAME_LEN


!******************************************************************************
!* FUNCTION  XF_SET_XSECT_POINT_MEASURE
!* PURPOSE   Set the Cross Section Point Property Measure
!* NOTES  
!************************************************************************************

SUBROUTINE XF_SET_XSECT_POINT_MEASURE (a_Id, a_index, a_NumVals, a_PropId, a_compression, error)
INTEGER,          INTENT(IN)  :: a_Id
INTEGER,          INTENT(IN)  ::  a_index,a_NumVals, a_compression
REAL*4,             INTENT(IN)  :: a_PropId
INTEGER,          INTENT(OUT) :: error
  

!            INTEGER, EXTERNAL :: _xfsetxsectpointpropmeasure_f
!  MS FORTRAN needs explicit interface for C functions called here.

  INTERFACE
    INTEGER FUNCTION xfsetxsectpointpropmeasure_f (a_Id, a_index, a_NumVals, a_PropId, a_compression)
    
	  !MS$ATTRIBUTES C,reference::xfsetxsectpointpropmeasure_f
    INTEGER,          INTENT(IN)  :: a_Id
    INTEGER,          INTENT(IN)  :: a_index, a_NumVals, a_compression
    REAL*4,             INTENT(IN)  :: a_PropId


    END FUNCTION xfsetxsectpointpropmeasure_f
  END INTERFACE


  error = xfsetxsectpointpropmeasure_f (a_Id, a_index, a_NumVals, a_PropId, a_compression)

  return

END SUBROUTINE ! XF_SET_XSECT_POINT_MEASURE

!******************************************************************************
!* FUNCTION  XF_GET_XSECT_POINT_MEASURE
!* PURPOSE   Get the Cross Section Point Property Measure
!* NOTES  
!************************************************************************************

SUBROUTINE XF_GET_XSECT_POINT_MEASURE (a_Id, a_index, a_NumVals, a_PropId, error)
INTEGER,          INTENT(IN)  :: a_Id
INTEGER,          INTENT(IN)  :: a_index, a_NumVals
REAL*8,           INTENT(IN)  :: a_PropId
INTEGER,          INTENT(OUT) :: error
  

!            INTEGER, EXTERNAL :: _xfgetxsectpointpropmeasure_f
!  MS FORTRAN needs explicit interface for C functions called here.

  INTERFACE
    INTEGER FUNCTION xfgetxsectpointpropmeasure_f (a_Id, a_index, a_NumVals, a_PropId)
    
	  !MS$ATTRIBUTES C,reference::xfgetxsectpointpropmeasure_f
    INTEGER,          INTENT(IN)  :: a_Id
    INTEGER,          INTENT(IN)  :: a_index, a_NumVals
    REAL*8,           INTENT(IN)  :: a_PropId


    END FUNCTION xfgetxsectpointpropmeasure_f
  END INTERFACE


  error = xfgetxsectpointpropmeasure_f (a_Id, a_index, a_NumVals, a_PropId)

  return

END SUBROUTINE ! XF_GET_XSECT_POINT_MEASURE


!******************************************************************************
!* FUNCTION  XF_SET_XSECT_POINT_TYPE
!* PURPOSE   Set the Cross Section Point Property Measure
!* NOTES  
!************************************************************************************

SUBROUTINE XF_SET_XSECT_POINT_TYPE (a_Id, a_index, a_NumVals, a_PropId, a_compression, error)
INTEGER,          INTENT(IN)  :: a_Id
INTEGER,          INTENT(IN)  :: a_index,a_NumVals, a_compression, a_PropId
INTEGER,          INTENT(OUT) :: error
  

!            INTEGER, EXTERNAL :: _xfsetxsectpointproptype_f
!  MS FORTRAN needs explicit interface for C functions called here.

  INTERFACE
    INTEGER FUNCTION xfsetxsectpointproptype_f (a_Id, a_index, a_NumVals, a_PropId, a_compression)
    
	  !MS$ATTRIBUTES C,reference::xfsetxsectpointproptype_f
    INTEGER,          INTENT(IN)  :: a_Id
    INTEGER,          INTENT(IN)  :: a_index, a_NumVals, a_compression, a_PropId


    END FUNCTION xfsetxsectpointproptype_f
  END INTERFACE


  error = xfsetxsectpointproptype_f (a_Id, a_index, a_NumVals, a_PropId, a_compression)

  return

END SUBROUTINE ! XF_SET_XSECT_POINT_TYPE

!******************************************************************************
!* FUNCTION  XF_GET_XSECT_POINT_TYPE
!* PURPOSE   Get the Cross Section Point Property Measure
!* NOTES  
!************************************************************************************

SUBROUTINE XF_GET_XSECT_POINT_TYPE (a_Id, a_index, a_NumVals, a_PropId, error)
INTEGER,          INTENT(IN)  :: a_Id
INTEGER,          INTENT(IN)  :: a_index, a_NumVals, a_PropId
INTEGER,          INTENT(OUT) :: error
  

!            INTEGER, EXTERNAL :: xfgetxsectpointproptype_f
!  MS FORTRAN needs explicit interface for C functions called here.

  INTERFACE
    INTEGER FUNCTION xfgetxsectpointproptype_f (a_Id, a_index, a_NumVals, a_PropId)
    
	  !MS$ATTRIBUTES C,reference::xfgetxsectpointproptype_f
    INTEGER,          INTENT(IN)  :: a_Id
    INTEGER,          INTENT(IN)  :: a_index, a_NumVals, a_PropId


    END FUNCTION xfgetxsectpointproptype_f
  END INTERFACE


  error = xfgetxsectpointproptype_f (a_Id, a_index, a_NumVals, a_PropId)

  return

END SUBROUTINE ! XF_GET_XSECT_POINT_TYPE


!******************************************************************************
!* FUNCTION  XF_SET_XSECT_POINT_ID
!* PURPOSE   Set the Cross Section Point Property Measure
!* NOTES  
!************************************************************************************

SUBROUTINE XF_SET_XSECT_POINT_ID (a_Id, a_NumVals, a_PropId, a_compression, error)
INTEGER,          INTENT(IN)  :: a_Id
INTEGER,          INTENT(IN)  :: a_NumVals, a_compression, a_PropId
INTEGER,          INTENT(OUT) :: error
  

!            INTEGER, EXTERNAL :: _xfsetxsectpointpropid_f
!  MS FORTRAN needs explicit interface for C functions called here.

  INTERFACE
    INTEGER FUNCTION xfsetxsectpointpropid_f (a_Id, a_NumVals, a_PropId, a_compression)
    
	  !MS$ATTRIBUTES C,reference::xfsetxsectpointpropid_f
    INTEGER,          INTENT(IN)  :: a_Id
    INTEGER,          INTENT(IN)  :: a_NumVals, a_compression, a_PropId


    END FUNCTION xfsetxsectpointpropid_f
  END INTERFACE


  error = xfsetxsectpointpropid_f (a_Id, a_NumVals, a_PropId, a_compression)

  return

END SUBROUTINE ! XF_SET_XSECT_POINT_ID

!******************************************************************************
!* FUNCTION  XF_GET_XSECT_POINT_ID
!* PURPOSE   Get the Cross Section Point Property Measure
!* NOTES  
!************************************************************************************

SUBROUTINE XF_GET_XSECT_POINT_ID (a_Id, a_NumVals, a_PropId, error)
INTEGER,          INTENT(IN)  :: a_Id
INTEGER,          INTENT(IN)  :: a_NumVals, a_PropId
INTEGER,          INTENT(OUT) :: error
  

!            INTEGER, EXTERNAL :: _xfgetxsectpointpropid_f
!  MS FORTRAN needs explicit interface for C functions called here.

  INTERFACE
    INTEGER FUNCTION xfgetxsectpointpropid_f (a_Id, a_NumVals, a_PropId)
    
	  !MS$ATTRIBUTES C,reference::xfgetxsectpointpropid_f
    INTEGER,          INTENT(IN)  :: a_Id
    INTEGER,          INTENT(IN)  :: a_NumVals, a_PropId


    END FUNCTION xfgetxsectpointpropid_f
  END INTERFACE


  error = xfgetxsectpointpropid_f (a_Id, a_NumVals, a_PropId)

  return

END SUBROUTINE ! XF_GET_XSECT_POINT_ID

!******************************************************************************
!* FUNCTION  XF_SET_XSECT_POINT_NAME
!* PURPOSE   Set the Cross Sections Point Property Name
!* NOTES  
!************************************************************************************

SUBROUTINE XF_SET_XSECT_POINT_NAME (a_Id, a_NumVals, a_PropId, a_StrLen, a_compression, error)
INTEGER,          INTENT(IN)  :: a_Id
INTEGER,          INTENT(IN)  :: a_NumVals, a_StrLen, a_compression
CHARACTER(LEN=*), INTENT(IN)  :: a_PropId
INTEGER,          INTENT(OUT) :: error
  

!            INTEGER, EXTERNAL :: _xfsetxsectpointpropname_f
!  MS FORTRAN needs explicit interface for C functions called here.

  INTERFACE
    INTEGER FUNCTION xfsetxsectpointpropname_f (a_Id, a_NumVals, a_PropId, a_StrLen, a_compression)
    
	  !MS$ATTRIBUTES C,reference::xfsetxsectpointpropname_f
    !DEC$ATTRIBUTES reference :: a_PropId
    INTEGER,          INTENT(IN)  :: a_Id
    INTEGER,          INTENT(IN)  :: a_NumVals, a_StrLen, a_compression
    CHARACTER(LEN=*), INTENT(IN)  :: a_PropId


    END FUNCTION xfsetxsectpointpropname_f
  END INTERFACE


  error = xfsetxsectpointpropname_f (a_Id, a_NumVals, a_PropId, a_StrLen, a_compression)

  return

END SUBROUTINE ! XF_SET_XSECT_POINT_NAME

!******************************************************************************
!* FUNCTION  XF_GET_XSECT_POINT_NAME
!* PURPOSE   Get the Cross Sections Name
!* NOTES  
!************************************************************************************

SUBROUTINE XF_GET_XSECT_POINT_NAME (a_Id, a_NumVals, a_PropId, a_StrLen, error)
INTEGER,          INTENT(IN)  :: a_Id
INTEGER,          INTENT(IN)  :: a_NumVals, a_StrLen
CHARACTER(LEN=*), INTENT(IN)  :: a_PropId
INTEGER,          INTENT(OUT) :: error
  

!            INTEGER, EXTERNAL :: _xfgetxsectpointpropname_f
!  MS FORTRAN needs explicit interface for C functions called here.

  INTERFACE
    INTEGER FUNCTION xfgetxsectpointpropname_f (a_Id, a_NumVals, a_PropId, a_StrLen)
    
	  !MS$ATTRIBUTES C,reference::xfgetxsectpointpropname_f
    !DEC$ATTRIBUTES reference :: a_PropId
    INTEGER,          INTENT(IN)  :: a_Id
    INTEGER,          INTENT(IN)  :: a_NumVals, a_StrLen
    CHARACTER(LEN=*), INTENT(IN)  :: a_PropId


    END FUNCTION xfgetxsectpointpropname_f
  END INTERFACE


  error = xfgetxsectpointpropname_f (a_Id, a_NumVals, a_PropId, a_StrLen)

  return

END SUBROUTINE ! XF_GET_XSECT_POINT_NAME


!******************************************************************************
!* FUNCTION  XF_GET_XSECT_POINT_NAME_LEN
!* PURPOSE   Get the Cross Sections Name Length
!* NOTES  
!************************************************************************************

SUBROUTINE XF_GET_XSECT_POINT_NAME_LEN (a_Id, a_NumVals, a_StrLen, error)
INTEGER,          INTENT(IN)  :: a_Id
INTEGER,          INTENT(IN)  :: a_NumVals, a_StrLen
INTEGER,          INTENT(OUT) :: error
  

!            INTEGER, EXTERNAL :: _xfgetxsectpointpropnamelen_f
!  MS FORTRAN needs explicit interface for C functions called here.

  INTERFACE
    INTEGER FUNCTION xfgetxsectpointpropnamelen_f (a_Id, a_NumVals, a_StrLen)
    
	  !MS$ATTRIBUTES C,reference::xfgetxsectpointpropnamelen_f
    INTEGER,          INTENT(IN)  :: a_Id
    INTEGER,          INTENT(IN)  :: a_NumVals, a_StrLen


    END FUNCTION xfgetxsectpointpropnamelen_f
  END INTERFACE


  error = xfgetxsectpointpropnamelen_f (a_Id, a_NumVals, a_StrLen)

  return

END SUBROUTINE ! XF_GET_XSECT_POINT_NAME_LEN

!******************************************************************************
!* FUNCTION  XF_SET_XSECT_POINT_DESC
!* PURPOSE   Set the Cross Sections Point Property Description
!* NOTES  
!************************************************************************************

SUBROUTINE XF_SET_XSECT_POINT_DESC (a_Id, a_NumVals, a_PropId, a_StrLen, a_compression, error)
INTEGER,          INTENT(IN)  :: a_Id
INTEGER,          INTENT(IN)  :: a_NumVals, a_StrLen, a_compression
CHARACTER(LEN=*), INTENT(IN)  :: a_PropId
INTEGER,          INTENT(OUT) :: error
  

!            INTEGER, EXTERNAL :: _xfsetxsectpointpropdesc_f
!  MS FORTRAN needs explicit interface for C functions called here.

  INTERFACE
    INTEGER FUNCTION xfsetxsectpointpropdesc_f (a_Id, a_NumVals, a_PropId, a_StrLen, a_compression)
    
	  !MS$ATTRIBUTES C,reference::xfsetxsectpointpropdesc_f
    !DEC$ATTRIBUTES reference :: a_PropId
    INTEGER,          INTENT(IN)  :: a_Id
    INTEGER,          INTENT(IN)  :: a_NumVals, a_StrLen, a_compression
    CHARACTER(LEN=*), INTENT(IN)  :: a_PropId


    END FUNCTION xfsetxsectpointpropdesc_f
  END INTERFACE


  error = xfsetxsectpointpropdesc_f (a_Id, a_NumVals, a_PropId, a_StrLen, a_compression)

  return

END SUBROUTINE ! XF_SET_XSECT_POINT_DESC

!******************************************************************************
!* FUNCTION  XF_GET_XSECT_POINT_DESC
!* PURPOSE   Get the Cross Sections Description
!* NOTES  
!************************************************************************************

SUBROUTINE XF_GET_XSECT_POINT_DESC (a_Id, a_NumVals, a_PropId, a_StrLen, error)
INTEGER,          INTENT(IN)  :: a_Id
INTEGER,          INTENT(IN)  :: a_NumVals, a_StrLen
CHARACTER(LEN=*), INTENT(IN)  :: a_PropId
INTEGER,          INTENT(OUT) :: error
  

!            INTEGER, EXTERNAL :: _xfgetxsectpointpropdesc_f
!  MS FORTRAN needs explicit interface for C functions called here.

  INTERFACE
    INTEGER FUNCTION xfgetxsectpointpropdesc_f (a_Id, a_NumVals, a_PropId, a_StrLen)
    
	  !MS$ATTRIBUTES C,reference::xfgetxsectpointpropdesc_f
    !DEC$ATTRIBUTES reference :: a_PropId
    INTEGER,          INTENT(IN)  :: a_Id
    INTEGER,          INTENT(IN)  :: a_NumVals, a_StrLen
    CHARACTER(LEN=*), INTENT(IN)  :: a_PropId


    END FUNCTION xfgetxsectpointpropdesc_f
  END INTERFACE


  error = xfgetxsectpointpropdesc_f (a_Id, a_NumVals, a_PropId, a_StrLen)

  return

END SUBROUTINE ! XF_GET_XSECT_POINT_DESC


!******************************************************************************
!* FUNCTION  XF_GET_XSECT_POINT_DESC_LEN
!* PURPOSE   Get the Cross Sections Description Length
!* NOTES  
!************************************************************************************

SUBROUTINE XF_GET_XSECT_POINT_DESC_LEN (a_Id, a_NumVals, a_StrLen, error)
INTEGER,          INTENT(IN)  :: a_Id
INTEGER,          INTENT(IN)  :: a_NumVals, a_StrLen
INTEGER,          INTENT(OUT) :: error
  

!            INTEGER, EXTERNAL :: _xfgetxsectpointpropdesclen_f
!  MS FORTRAN needs explicit interface for C functions called here.

  INTERFACE
    INTEGER FUNCTION xfgetxsectpointpropdesclen_f (a_Id, a_NumVals, a_StrLen)
    
	  !MS$ATTRIBUTES C,reference::xfgetxsectpointpropdesclen_f
    INTEGER,          INTENT(IN)  :: a_Id
    INTEGER,          INTENT(IN)  :: a_NumVals, a_StrLen


    END FUNCTION xfgetxsectpointpropdesclen_f
  END INTERFACE


  error = xfgetxsectpointpropdesclen_f (a_Id, a_NumVals, a_StrLen)

  return

END SUBROUTINE ! XF_GET_XSECT_POINT_NAME_LEN


!******************************************************************************
!* FUNCTION  XF_SET_XSECT_POINT_EXCLUSIVE
!* PURPOSE   Set the Cross Section Point Property Exclusive
!* NOTES  
!************************************************************************************

SUBROUTINE XF_SET_XSECT_POINT_EXCLUSIVE (a_Id, a_NumVals, a_PropId, a_compression, error)
INTEGER,          INTENT(IN)  :: a_Id
INTEGER,          INTENT(IN)  :: a_NumVals, a_compression, a_PropId
INTEGER,          INTENT(OUT) :: error
  

!            INTEGER, EXTERNAL :: _xfSetXSectPointPropExclusive_f
!  MS FORTRAN needs explicit interface for C functions called here.

  INTERFACE
    INTEGER FUNCTION xfsetxsectpointpropexclusive_f (a_Id, a_NumVals, a_PropId, a_compression)
    
	  !MS$ATTRIBUTES C,reference::xfsetxsectpointpropexclusive_f
    INTEGER,          INTENT(IN)  :: a_Id
    INTEGER,          INTENT(IN)  :: a_NumVals, a_compression, a_PropId


    END FUNCTION xfsetxsectpointpropexclusive_f
  END INTERFACE


  error = xfsetxsectpointpropexclusive_f (a_Id, a_NumVals, a_PropId, a_compression)

  return

END SUBROUTINE ! XF_SET_XSECT_POINT_EXCLUSIVE

!******************************************************************************
!* FUNCTION  XF_GET_XSECT_POINT_EXCLUSIVE
!* PURPOSE   Get the Cross Section Point Property Exclusive
!* NOTES  
!************************************************************************************

SUBROUTINE XF_GET_XSECT_POINT_EXCLUSIVE (a_Id, a_NumVals, a_PropId, error)
INTEGER,          INTENT(IN)  :: a_Id
INTEGER,          INTENT(IN)  :: a_NumVals, a_PropId
INTEGER,          INTENT(OUT) :: error
  

!            INTEGER, EXTERNAL :: _xfgetxsectpointpropexclusive_f
!  MS FORTRAN needs explicit interface for C functions called here.

  INTERFACE
    INTEGER FUNCTION xfgetxsectpointpropexclusive_f (a_Id, a_NumVals, a_PropId)
    
	  !MS$ATTRIBUTES C,reference::xfgetxsectpointpropexclusive_f
    INTEGER,          INTENT(IN)  :: a_Id
    INTEGER,          INTENT(IN)  :: a_NumVals, a_PropId


    END FUNCTION xfgetxsectpointpropexclusive_f
  END INTERFACE


  error = xfgetxsectpointpropexclusive_f (a_Id, a_NumVals, a_PropId)

  return

END SUBROUTINE ! XF_GET_XSECT_POINT_EXCLUSIVE


!******************************************************************************
!* FUNCTION  XF_SET_XSECT_TOPO_NAME
!* PURPOSE   Set the Cross Sections Topo Name
!* NOTES  
!************************************************************************************

SUBROUTINE XF_SET_XSECT_TOPO_NAME (a_Id, a_NumVals, a_PropId, a_StrLen, a_compression, error)
INTEGER,          INTENT(IN)  :: a_Id
INTEGER,          INTENT(IN)  :: a_NumVals, a_StrLen, a_compression
CHARACTER(LEN=*), INTENT(IN)  :: a_PropId
INTEGER,          INTENT(OUT) :: error
  

!            INTEGER, EXTERNAL :: _xfsetxsecttoponame_f
!  MS FORTRAN needs explicit interface for C functions called here.

  INTERFACE
    INTEGER FUNCTION xfsetxsecttoponame_f (a_Id, a_NumVals, a_PropId, a_StrLen, a_compression)
    
	  !MS$ATTRIBUTES C,reference::xfsetxsecttoponame_f
    !DEC$ATTRIBUTES reference :: a_PropId
    INTEGER,          INTENT(IN)  :: a_Id
    INTEGER,          INTENT(IN)  :: a_NumVals, a_StrLen, a_compression
    CHARACTER(LEN=*), INTENT(IN)  :: a_PropId


    END FUNCTION xfsetxsecttoponame_f
  END INTERFACE


  error = xfsetxsecttoponame_f (a_Id, a_NumVals, a_PropId, a_StrLen, a_compression)

  return

END SUBROUTINE ! XF_SET_XSECT_TOPO_NAME

!******************************************************************************
!* FUNCTION  XF_GET_XSECT_TOPO_NAME
!* PURPOSE   Get the Cross Sections Topo Name
!* NOTES  
!************************************************************************************

SUBROUTINE XF_GET_XSECT_TOPO_NAME (a_Id, a_NumVals, a_PropId, a_StrLen, error)
INTEGER,          INTENT(IN)  :: a_Id
INTEGER,          INTENT(IN)  :: a_NumVals, a_StrLen
CHARACTER(LEN=*), INTENT(IN)  :: a_PropId
INTEGER,          INTENT(OUT) :: error
  

!            INTEGER, EXTERNAL :: _xfgetxsecttoponame_f
!  MS FORTRAN needs explicit interface for C functions called here.

  INTERFACE
    INTEGER FUNCTION xfgetxsecttoponame_f (a_Id, a_NumVals, a_PropId, a_StrLen)
    
	  !MS$ATTRIBUTES C,reference::xfgetxsecttoponame_f
    !DEC$ATTRIBUTES reference :: a_PropId
    INTEGER,          INTENT(IN)  :: a_Id
    INTEGER,          INTENT(IN)  :: a_NumVals, a_StrLen
    CHARACTER(LEN=*), INTENT(IN)  :: a_PropId


    END FUNCTION xfgetxsecttoponame_f
  END INTERFACE


  error = xfgetxsecttoponame_f (a_Id, a_NumVals, a_PropId, a_StrLen)

  return

END SUBROUTINE ! XF_GET_XSECT_TOPO_NAME


!******************************************************************************
!* FUNCTION  XF_GET_XSECT_TOPO_NAME_LEN
!* PURPOSE   Get the Cross Sections Topo Name Length
!* NOTES  
!************************************************************************************

SUBROUTINE XF_GET_XSECT_TOPO_NAME_LEN (a_Id, a_NumVals, a_StrLen, error)
INTEGER,          INTENT(IN)  :: a_Id
INTEGER,          INTENT(IN)  :: a_NumVals, a_StrLen
INTEGER,          INTENT(OUT) :: error
  

!            INTEGER, EXTERNAL :: _xfgetxsecttoponamelen_f
!  MS FORTRAN needs explicit interface for C functions called here.

  INTERFACE
    INTEGER FUNCTION xfgetxsecttoponamelen_f (a_Id, a_NumVals, a_StrLen)
    
	  !MS$ATTRIBUTES C,reference::xfgetxsecttoponamelen_f
    INTEGER,          INTENT(IN)  :: a_Id
    INTEGER,          INTENT(IN)  :: a_NumVals, a_StrLen


    END FUNCTION xfgetxsecttoponamelen_f
  END INTERFACE


  error = xfgetxsecttoponamelen_f (a_Id, a_NumVals, a_StrLen)

  return

END SUBROUTINE ! XF_GET_XSECT_TOPO_NAME_LEN


!******************************************************************************
!* FUNCTION  XF_SET_XSECT_TOPO_DESC
!* PURPOSE   Set the Cross Sections Topo Name
!* NOTES  
!************************************************************************************

SUBROUTINE XF_SET_XSECT_TOPO_DESC (a_Id, a_NumVals, a_PropId, a_StrLen, a_compression, error)
INTEGER,          INTENT(IN)  :: a_Id
INTEGER,          INTENT(IN)  :: a_NumVals, a_StrLen, a_compression
CHARACTER(LEN=*), INTENT(IN)  :: a_PropId
INTEGER,          INTENT(OUT) :: error
  

!            INTEGER, EXTERNAL :: _xfsetxsecttopodesc_f
!  MS FORTRAN needs explicit interface for C functions called here.

  INTERFACE
    INTEGER FUNCTION xfsetxsecttopodesc_f (a_Id, a_NumVals, a_PropId, a_StrLen, a_compression)
    
	  !MS$ATTRIBUTES C,reference::xfsetxsecttopodesc_f
    !DEC$ATTRIBUTES reference :: a_PropId
    INTEGER,          INTENT(IN)  :: a_Id
    INTEGER,          INTENT(IN)  :: a_NumVals, a_StrLen, a_compression
    CHARACTER(LEN=*), INTENT(IN)  :: a_PropId


    END FUNCTION xfsetxsecttopodesc_f
  END INTERFACE


  error = xfsetxsecttopodesc_f (a_Id, a_NumVals, a_PropId, a_StrLen, a_compression)

  return

END SUBROUTINE ! XF_SET_XSECT_TOPO_DESC

!******************************************************************************
!* FUNCTION  XF_GET_XSECT_TOPO_DESC
!* PURPOSE   Get the Cross Sections Topo Name
!* NOTES  
!************************************************************************************

SUBROUTINE XF_GET_XSECT_TOPO_DESC (a_Id, a_NumVals, a_PropId, a_StrLen, error)
INTEGER,          INTENT(IN)  :: a_Id
INTEGER,          INTENT(IN)  :: a_NumVals, a_StrLen
CHARACTER(LEN=*), INTENT(IN)  :: a_PropId
INTEGER,          INTENT(OUT) :: error
  

!            INTEGER, EXTERNAL :: _xfgetxsecttopodesc_f
!  MS FORTRAN needs explicit interface for C functions called here.

  INTERFACE
    INTEGER FUNCTION xfgetxsecttopodesc_f (a_Id, a_NumVals, a_PropId, a_StrLen)
    
	  !MS$ATTRIBUTES C,reference::xfgetxsecttopodesc_f
    !DEC$ATTRIBUTES reference :: a_PropId
    INTEGER,          INTENT(IN)  :: a_Id
    INTEGER,          INTENT(IN)  :: a_NumVals, a_StrLen
    CHARACTER(LEN=*), INTENT(IN)  :: a_PropId


    END FUNCTION xfgetxsecttopodesc_f
  END INTERFACE


  error = xfgetxsecttopodesc_f (a_Id, a_NumVals, a_PropId, a_StrLen)

  return

END SUBROUTINE ! XF_GET_XSECT_TOPO_DESC


!******************************************************************************
!* FUNCTION  XF_GET_XSECT_TOPO_DESC_LEN
!* PURPOSE   Get the Cross Sections Topo Name Length
!* NOTES  
!************************************************************************************

SUBROUTINE XF_GET_XSECT_TOPO_DESC_LEN (a_Id, a_NumVals, a_StrLen, error)
INTEGER,          INTENT(IN)  :: a_Id
INTEGER,          INTENT(IN)  :: a_NumVals, a_StrLen
INTEGER,          INTENT(OUT) :: error
  

!            INTEGER, EXTERNAL :: _xfgetxsecttopodesclen_f
!  MS FORTRAN needs explicit interface for C functions called here.

  INTERFACE
    INTEGER FUNCTION xfgetxsecttopodesclen_f (a_Id, a_NumVals, a_StrLen)
    
	  !MS$ATTRIBUTES C,reference::xfgetxsecttopodesclen_f
    INTEGER,          INTENT(IN)  :: a_Id
    INTEGER,          INTENT(IN)  :: a_NumVals, a_StrLen


    END FUNCTION xfgetxsecttopodesclen_f
  END INTERFACE


  error = xfgetxsecttopodesclen_f (a_Id, a_NumVals, a_StrLen)

  return

END SUBROUTINE ! XF_GET_XSECT_TOPO_DESC_LEN

!******************************************************************************
!* FUNCTION  XF_SET_WKT
!* PURPOSE   Set the WKT string 
!* NOTES     
!************************************************************************************

SUBROUTINE XF_SET_WKT(a_Id, a_WKT, error)
INTEGER,          INTENT(IN)  :: a_Id
CHARACTER(LEN=*), INTENT(IN)  :: a_WKT
INTEGER,          INTENT(OUT) :: error
INTEGER                       :: StrLen

!            INTEGER, EXTERNAL :: _xfsetwkt_f
!  MS FORTRAN needs explicit interface for C functions called here.

  INTERFACE
    INTEGER FUNCTION xfsetwkt_f (a_Id, a_WKT, a_StrLen)
    
	  !MS$ATTRIBUTES C,reference::xfsetwkt_f
    INTEGER,          INTENT(IN)  :: a_Id
    CHARACTER(LEN=*), INTENT(IN)  :: a_WKT
    INTEGER,          INTENT(IN)  :: a_StrLen

    END FUNCTION xfsetwkt_f
  END INTERFACE
  
  StrLen = LEN_TRIM(a_WKT)
  error = xfsetwkt_f(a_Id, a_WKT, StrLen)
  
  return
  
END SUBROUTINE

!******************************************************************************
!* FUNCTION  XF_SET_ATTRIBUTE_STRING
!* PURPOSE   Write the string as an attribute
!* NOTES     
!************************************************************************************

SUBROUTINE XF_SET_ATTRIBUTE_STRING(a_Id, a_Name, a_String, error)
INTEGER,          INTENT(IN)  :: a_Id
CHARACTER(LEN=*), INTENT(IN)  :: a_Name, a_String
INTEGER,          INTENT(OUT) :: error
INTEGER                       :: name_len, str_len

!            INTEGER, EXTERNAL :: _xfsetattributestring_f
!  MS FORTRAN needs explicit interface for C functions called here.

  INTERFACE
    INTEGER FUNCTION xfsetattributestring_f (a_Id, a_Name, name_len, a_String, str_len)
    
	  !MS$ATTRIBUTES C,reference::xfsetattributestring_f
    INTEGER,          INTENT(IN)  :: a_Id
    CHARACTER(LEN=*), INTENT(IN)  :: a_Name, a_String
    INTEGER,          INTENT(IN)  :: name_len, str_len

    END FUNCTION xfsetattributestring_f
  END INTERFACE
  
  name_len = LEN_TRIM(a_Name)
  str_len = LEN_TRIM(a_String)
  error = xfsetattributestring_f(a_Id, a_Name, name_len, a_String, str_len)
  
  return
  
END SUBROUTINE

!******************************************************************************
!* FUNCTION  XF_GET_WKT_STRING_SIZE
!* PURPOSE   Get the size of the WKT String 
!* NOTES     
!************************************************************************************

SUBROUTINE XF_GET_WKT_STRING_SIZE (a_Id, a_StrLen, error)
INTEGER,          INTENT(IN)  :: a_Id
INTEGER,          INTENT(IN)  :: a_StrLen
INTEGER,          INTENT(OUT) :: error
  

!            INTEGER, EXTERNAL :: _xfgetwktstringsize_f
!  MS FORTRAN needs explicit interface for C functions called here.

  INTERFACE
    INTEGER FUNCTION xfgetwktstringsize_f (a_Id, a_StrLen)
    
	  !MS$ATTRIBUTES C,reference::xfgetwktstringsize_f
    INTEGER,          INTENT(IN)  :: a_Id
    INTEGER,          INTENT(IN)  :: a_StrLen


    END FUNCTION xfgetwktstringsize_f
  END INTERFACE


  error = xfgetwktstringsize_f (a_Id, a_StrLen)

  return

END SUBROUTINE ! XF_GET_WKT_STRING_SIZE

!******************************************************************************
!* FUNCTION  XF_GET_WKT
!* PURPOSE   Get the WKT String 
!* NOTES     
!************************************************************************************

SUBROUTINE XF_GET_WKT (a_Id, a_String, error)
INTEGER,          INTENT(IN)  :: a_Id
CHARACTER(LEN=*), INTENT(IN)  :: a_String
INTEGER,          INTENT(OUT) :: error
  

!            INTEGER, EXTERNAL :: _xfgetwkt_f
!  MS FORTRAN needs explicit interface for C functions called here.

  INTERFACE
    INTEGER FUNCTION xfgetwkt_f (a_Id, a_String)
    
	  !MS$ATTRIBUTES C,reference::xfgetwkt_f
    INTEGER,          INTENT(IN)  :: a_Id
    CHARACTER(LEN=*), INTENT(IN)  :: a_String


    END FUNCTION xfgetwkt_f
  END INTERFACE


  error = xfgetwkt_f (a_Id, a_String)

  return

END SUBROUTINE ! XF_GET_WKT

!******************************************************************************
!* FUNCTION  XF_GET_COORD_VERSION
!* PURPOSE   Get the Coordinate Version
!* NOTES     
!************************************************************************************

SUBROUTINE XF_GET_COORD_VERSION (a_Id, a_Version, error)
INTEGER,  INTENT(IN)  :: a_Id
INTEGER,  INTENT(IN)  :: a_Version
INTEGER,  INTENT(OUT) :: error
  

!            INTEGER, EXTERNAL :: _xfgetcoordversion_f
!  MS FORTRAN needs explicit interface for C functions called here.

  INTERFACE
    INTEGER FUNCTION xfgetcoordversion_f (a_Id, a_Version)
    
	  !MS$ATTRIBUTES C,reference::xfgetcoordversion_f
    INTEGER, INTENT(IN)  :: a_Id
    INTEGER, INTENT(IN)  :: a_Version


    END FUNCTION xfgetcoordversion_f
  END INTERFACE


  error = xfgetcoordversion_f (a_Id, a_Version)

  return

END SUBROUTINE ! XF_GET_COORD_VERSION


END MODULE XMDF
