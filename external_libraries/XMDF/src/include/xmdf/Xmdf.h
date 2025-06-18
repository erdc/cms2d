#ifndef XMDF_DOT_H
#define XMDF_DOT_H

/*! \file xmdf.h */
/*
XMDF License & Copyright Notice
Copyright Notice and Statement for the eXtensible Model Data Format (XMDF) Software Library and API.
XMDF Software Library and Utilities Copyright 2004, 2005 by the Enivironmental Modeling Research Laboratory and Brigham Young University.  All rights reserved.
Redistribution and use in source and binary forms, with or without modification, are permitted for any purpose (including commercial purposes) provided that the following conditions are met:
1.  Redistributions of source code must retain the above copyright notice, this list of conditions, and the following disclaimer.
2.  Redistributions in binary form must reproduce the above copyright notice, this list of conditions, and the following disclaimer in the documentation and/or materials provided with the distribution.
3.  In addition, redistributions of modified forms of the source or binary code must carry prominent notices stating that the original code was changed and the date of the change.
4.  All publications or advertising materials mentioning features or use of this software are asked, but not required, to acknowledge that it was developed by the EMRL at Brigham Young University in partnership with the US Army Engineer Research and Development Center and to credit the contributors.
5.  Neither the name of the University nor the names of the Contributors maybe used to endorse or promote products derived from this software without specific prior written permission from the University or the Contributors, as appropriate for the name(s) to be used.
6.  THIS SOFTWARE IS PROVIDED BY THE UNIVERSITY AND THE CONTRIBUTORS "AS IS"WITH NO WARRANTY OF ANY KIND, EITHER EXPRESSED OR IMPLIED.  In no event shall the University or the Contributors be liable for any damages suffered by the users arising out of the use of this software, even if advised of the possibility of such damage.
--------------------------------------------------------------------------
Portions of XMDF utilize the HDF5 software library.  The following statement applies to HDF5:
Copyright Notice and Statement for NCSA Hierarchical Data Format (HDF) Software Library and Utilities
NCSA HDF5 (Hierarchical Data Format 5) Software Library and Utilities Copyright 1998, 1999, 2000, 2001, 2002, 2003, 2004, 2005 by the Board of Trustees of the University of Illinois.  All rights reserved.
Contributors: National Center for Supercomputing Applications (NCSA) at the University of Illinois at Urbana-Champaign (UIUC), Lawrence Livermore National Laboratory (LLNL), Sandia National Laboratories (SNL), Los Alamos National Laboratory (LANL), Jean-loup Gailly and Mark Adler (gzip library).
Redistribution and use in source and binary forms, with or without modification, are permitted for any purpose (including commercial purposes) provided that the following conditions are met:
1.  Redistributions of source code must retain the above copyright notice, this list of conditions, and the following disclaimer.
2.  Redistributions in binary form must reproduce the above copyright notice, this list of conditions, and the following disclaimer in the documentation and/or materials provided with the distribution.
3.  In addition, redistributions of modified forms of the source or binary code must carry prominent notices stating that the original code was changed and the date of the change.
4.  All publications or advertising materials mentioning features or use of this software are asked, but not required, to acknowledge that it was developed by the National Center for Supercomputing Applications at the University of Illinois at Urbana-Champaign and to credit the contributors.
5.  Neither the name of the University nor the names of the Contributors may be used to endorse or promote products derived from this software without specific prior written permission from the University or the Contributors, as appropriate for the name(s) to be used.
6.  THIS SOFTWARE IS PROVIDED BY THE UNIVERSITY AND THE CONTRIBUTORS "AS IS" WITH NO WARRANTY OF ANY KIND, EITHER EXPRESSED OR IMPLIED.  In no event shall the University or the Contributors be liable for any damages suffered by the users arising out of the use of this software, even if advised of the possibility of such damage.
--------------------------------------------------------------------------
Portions of HDF5 were developed with support from the University of California, Lawrence Livermore National Laboratory (UC LLNL).
The following statement applies to those portions of the product and must be retained in any redistribution of source code, binaries, documentation, and/or accompanying materials:
This work was partially produced at the University of California, Lawrence Livermore National Laboratory (UC LLNL) under contract no. W-7405-ENG-48 (Contract 48) between the U.S. Department of Energy (DOE) and The Regents of the University of California (University) for the operation of UC LLNL.
DISCLAIMER:
This work was prepared as an account of work sponsored by an agency of the United States Government.  Neither the United States Government nor the University of California nor any of their employees, makes any warranty, express or implied, or assumes any liability or responsibility for the accuracy, completeness, or usefulness of any information, apparatus, product, or process  disclosed, or represents that its use would not infringe privately- owned rights.  Reference herein to any specific commercial products,  process, or service by trade name, trademark, manufacturer, or  otherwise, does not necessarily constitute or imply its endorsement,  recommendation, or favoring by the United States Government or the  University of California.  The views and opinions of authors expressed herein do not necessarily state or reflect those of the United States Government or the University of California, and shall not be used for advertising or product endorsement purposes.
--------------------------------------------------------------------------
This work was partially produced at the University of California, Lawrence Livermore National Laboratory (UC LLNL) under contract no.W-7405-ENG-48 (Contract 48) between the U.S. Department of Energy (DOE) and The Regents of the University of California (University) for the operation of UC LLNL.
DISCLAIMER:This work was prepared as an account of work sponsored by an agency of the United States Government.  Neither the United States Government nor the University of California nor any of their employees, makes any warranty, express or implied, or assumes any liability or responsibility for the accuracy, completeness, or usefulness of any information, apparatus, product, or process disclosed, or represents that its use would not infringe privately-owned rights.  Reference herein to any specific commercial products, process, or service by trade name, trademark, manufacturer, or otherwise, does not necessarily constitute or imply its endorsement, recommendation, or favoring by the United States Government or the University of California.  The views and opinions of authors expressed herein do not necessarily state or reflect those of the United States Government or the University of California, and shall not be used for advertising or product endorsement purposes.
--------------------------------------------------------------------------
*/

/* The following ifdef block is the standard way of creating macros which make
   exporting from a DLL simpler. All files within this DLL are compiled with
   the XMDF_EXPORTS symbol defined on the command line. this symbol should
   not be defined on any project that uses this DLL. This way any other project
   whose source files include this file see  functions as being
   imported from a DLL, wheras this DLL sees symbols defined with this macro
   as being exported. */
#ifdef NO_XMDF_API
#define XMDF_API  
#else
#ifdef _WIN32
#ifdef XMDF_EXPORTS
#define XMDF_API __declspec(dllexport) 
#else
#define XMDF_API __declspec(dllimport)
#endif
#else
#define XMDF_API  // Shared library symbols are exported by default on non-windows platforms
#endif
#endif

/* XMDF version number */
/* Single precision float */ 
#define XMDF_VERSION  99.99


#include <hdf5.h>
#include <float.h>
#include <limits.h>
#include <math.h>

/* type defines */
#define xmbool unsigned char
#define xid hid_t

#define XFALSE 0
#define XTRUE  !XFALSE

#define MAX_ID_NAME 512

typedef struct {
  hid_t       id;
  H5I_type_t  type;
  char        name[MAX_ID_NAME];
} IdentifierInfo;


/* #defines */
#define NONE -1

#ifndef MAX_FLOAT
#define MAX_FLOAT FLT_MAX
#endif

#ifndef MIN_FLOAT
#define MIN_FLOAT (-FLT_MAX)
#endif

  /* File Version */
#define FILE_VERSION "File Version"

  /* Property Types */
#define XF_TYPE_INT     1
#define XF_TYPE_FLOAT   2
#define XF_TYPE_DOUBLE  3
#define XF_TYPE_STRING  4
#define XF_TYPE_UINT    5
#define XF_TYPE_OTHER   11

  /* Time Variables */
#define TIME_UNITS "TimeUnits"
#define TS_SECONDS "Seconds"
#define TS_MINUTES "Minutes"
#define TS_HOURS   "Hours"
#define TS_DAYS    "Days"
#define TS_NOT_APPLICABLE "None"
#define TIME_UNITS_MAXLENGTH 25

#define NO_UNITS "None"
#define UNITS_MAXLENGTH 100

/* Element types */
#define EL_TYPE_JUNCTION3              3
#define EL_TYPE_JUNCTION4              4
#define EL_TYPE_JUNCTION5              5
#define EL_TYPE_JUNCTION6              6
#define EL_TYPE_JUNCTION7              7
#define EL_TYPE_JUNCTION8              8
  /* 1D elements  */
#define EL_TYPE_1D_LINEAR            100 
#define EL_TYPE_1D_QUAD              101 
#define EL_TYPE_TRANS_1D_2D          110
  /* 2D elements  */
#define EL_TYPE_TRI_LINEAR           200
#define EL_TYPE_TRI_QUAD             201
#define EL_TYPE_QUADRILATERAL_LINEAR 210
#define EL_TYPE_QUADRILATERAL_QUAD8  211
#define EL_TYPE_QUADRILATERAL_QUAD9  212
  /* 3D elements */
#define EL_TYPE_TET_LINEAR           300
#define EL_TYPE_PRISM_LINEAR         310
#define EL_TYPE_HEX_LINEAR           320
#define EL_TYPE_PYRAMID_LINEAR       330

  /* Grid types */
#define GRID_TYPE_CARTESIAN            0
#define GRID_TYPE_CURVILINEAR          1
#define GRID_TYPE_CARTESIAN_EXTRUDED   2
#define GRID_TYPE_CURVILINEAR_EXTRUDED 3
#define GRID_TYPE_MIN                  0
#define GRID_TYPE_MAX                  3

  /* extrusion types for grids  */
#define EXTRUDE_SIGMA                  0
#define EXTRUDE_CARTESIAN              1
#define EXTRUDE_CURV_AT_CORNERS        2
#define EXTRUDE_CURV_AT_CELLS          3
#define EXTRUDE_MIN                    0
#define EXTRUDE_MAX                    3

 /* Grid orientation  */
#define ORIENTATION_RIGHT_HAND         0
#define ORIENTATION_LEFT_HAND          1

 /* Reserved property names */
#define PROP_ACTIVITY                 "activity"
#define PROP_NULL_VALUE               "nullvalue"

/* Dataset value locations for a grid */
#define GRID_LOC_CENTER                0
#define GRID_LOC_CORNER                1
#define GRID_LOC_FACES                 2
#define GRID_LOC_FACE_I                3
#define GRID_LOC_FACE_J                4
#define GRID_LOC_FACE_K                5
#define GRID_LOC_COLUMN                6
#define GRID_LOC_CENTER_CORNERS_FACES  7 /* Centers, corners, and faces */

/* #defines for coordinate conversions  */
  /* Horizontal datum  */
#define HORIZ_DATUM_LOCAL               0
#define HORIZ_DATUM_GEOGRAPHIC          1
#define HORIZ_DATUM_GEOGRAPHIC_NAD27    2
#define HORIZ_DATUM_GEOGRAPHIC_NAD83    3
#define HORIZ_DATUM_GEOGRAPHIC_HPGN     4
#define HORIZ_DATUM_UTM                 5
#define HORIZ_DATUM_UTM_NAD27           6
#define HORIZ_DATUM_UTM_NAD83           7
#define HORIZ_DATUM_UTM_HPGN            8
#define HORIZ_DATUM_STATE_PLANE_NAD27   9
#define HORIZ_DATUM_STATE_PLANE_NAD83  10
#define HORIZ_DATUM_STATE_PLANE_HPGN   11
#define HORIZ_DATUM_CPP                12

  /* Horizontal or vertical units  */
#define COORD_UNITS_US_FEET             0
#define COORD_UNITS_INTERNATIONAL_FEET  1
#define COORD_UNITS_METERS              2

  /* vertical datum */
#define VERT_DATUM_LOCAL                0
#define VERT_DATUM_NGVD_29              1
#define VERT_DATUM_NGVD_88              2

  /* latitude  */
#define LATITUDE_NORTH                  0
#define LATITUDE_SOUTH                  1

  /* longitude */
#define LONGITUDE_EAST                  0
#define LONGITUDE_WEST                  1

  /* valid UTM zones */
#define UTM_ZONE_MIN                    1
#define UTM_ZONE_MAX                   60

  /* Guid names and lengths */
#define XF_GUID                    "Guid"
#define XF_GUID_STRINGLENGTH           37

  /* Material Mapping File for Classifying */
#define XMDF_MAT_MAPING_TARGET_IDS     "Target Ids"
#define XMDF_MAT_MAPING_BACKGROUND_IDS "Background Ids"
#define XMDF_MAT_MAPING_IDS            "Mapping Ids"

/* Calendar Eras */
#define ERA_IS_BCE   0
#define ERA_IS_CE    1

/* Macros */
#define Xmax(a,b)    (((a) >= (b)) ? (a) : (b))
#define Xmin(a,b)    (((a) >= (b)) ? (b) : (a))

#define MAX_DOUBLE 1.7976931348623157e+308 
#define MIN_DOUBLE -1.7976931348623157e+308

/* Property for node ids */
#define XF_IDS                      "Ids"

/* Error Retrieval */
#define XF_MAX_ERROR_MSG_SIZE         256
#define XF_MAX_ERROR_STACK            256

/* Overwrite options in the function xfSetupToWriteDatasets */
#define XF_OVERWRITE_CLEAR_FILE           1
#define XF_OVERWRITE_CLEAR_DATASET_GROUP  2
#define XF_OVERWRITE_NONE                 3

/* File Functions */
#ifdef __cplusplus
extern "C" {
#endif

/* Library versions */
XMDF_API xid xfGetLibraryVersion(float *a_Version);
XMDF_API xid xfGetLibraryVersionFile(xid a_File, float *a_Version);
/* Coordinates */
XMDF_API xid xfCreateCoordinateGroup (xid a_ParentId, xid *a_ChildId);
XMDF_API xid xfOpenCoordinateGroup (xid a_ParentId, xid *a_ChildId);

XMDF_API xid xfCreateFile(const char *a_File, xid *Id, xmbool a_Overwrite);
XMDF_API xid xfCreateInMemoryFile(const char *a_File, xid *Id);

XMDF_API xid xfOpenFile(const char *a_File, xid *Id, xmbool a_ReadOnly);
XMDF_API xid xfCloseFile(xid a_FileId);

/* retreiving error messages */
/* error messages are cleared whenever a file is opened or closed */
XMDF_API xid xfGetNumErrorMessages(int *a_Num);
/* a_Errors must be allocated to a_Num * XF_MAX_ERROR_MSG_SIZE */
XMDF_API xid xfGetErrorMessages(int a_Num, char **a_Errors);

/* functions to determine and close identifiers not yet closed */
XMDF_API xid xfGetNumOpenIdentifiers(xid a_FileId, int *a_Num);
XMDF_API xid xfGetOpenIdentifiersInfo(xid a_FileId, int a_Num,
                                      IdentifierInfo *a_Info);

/* Functions used to work with Properties */
XMDF_API xid xfCreatePropertyGroup(xid a_ParentId, xid *Id);

XMDF_API xid xfWritePropertyString(xid a_Id, const char *a_Name, 
                   int a_Number, int a_StringLength, const char *a_Properties);
XMDF_API xid xfWritePropertyUnsignedInt(xid a_Id, const char *a_Name, 
                   int a_Number, const unsigned int *a_Properties, 
                   int a_Compression);
XMDF_API xid xfWritePropertyInt(xid a_Id, const char *a_Name, 
                   int a_Number, const int *a_Properties,
                   int a_Compression);
XMDF_API xid xfWritePropertyDouble(xid a_Id, const char *a_Name, 
                   int a_Number, const double *a_Properties,
                   int a_Compression);
XMDF_API xid xfWritePropertyFloat(xid a_Id, const char *a_Name, 
                   int a_Number, const float *a_Properties,
                   int a_Compression);

XMDF_API xid xfDoesPropertyWithNameExist(xid a_Id, const char *a_Name,
                                             xmbool *Exists);
XMDF_API xid xfGetPropertyNumber(xid a_Id, const char *a_Name,
                                    int *a_Number);
XMDF_API xid xfGetPropertyStringLength(xid a_Id, const char *a_Name, 
                                   int *a_Number, int *a_MaxLength);
XMDF_API xid xfGetPropertyType(xid a_GroupId, const char *a_Name, int *a_Type);
XMDF_API xid xfAllocateReadPropertyString (xid a_Id, const char *a_Name,
                                           int *a_Number, int *a_MaxLength,
                                           char **a_Properties);
XMDF_API xid xfReadPropertyString(xid a_Id, const char *a_Name,
                    int a_Number, int a_MaxLength, char *a_Properties);
XMDF_API xid xfReadPropertyUnsignedInt(xid a_Id, const char *a_Name,
                    int a_Number, unsigned int *a_Properties);
XMDF_API xid xfReadPropertyInt(xid a_Id, const char *a_Name,
                    int a_Number, int *a_Properties);
XMDF_API xid xfReadPropertyDouble(xid a_Id, const char *a_Name,
                    int a_Number, double *a_Properties);
XMDF_API xid xfReadPropertyFloat(xid a_Id, const char *a_Name,
                    int a_Number, float *a_Properties);

  /* group functions */
XMDF_API xid xfCreateGenericGroup(xid a_FileId, const char *a_Path,
                                     xid *a_GroupId);
XMDF_API xid xfCreateGroupForMesh(xid a_FileId, const char *a_Path,
                                         xid *a_GroupId);
XMDF_API xid xfCreateGroupForGrid(xid a_FileId, const char *a_Path,
                                         xid *a_GroupId);
XMDF_API xid xfCreateStationForGrid(xid a_FileId, const char *a_Path,
                                         xid *a_GroupId);
XMDF_API xid xfCreateGroupForXsec(xid a_FileId, const char *a_Path,
                                         xid *a_GroupId);
XMDF_API xid xfCreateGroupForMatSim(xid a_FileId, const char *a_Path,
                                         xid *a_GroupId);

XMDF_API xid xfOpenGroup (xid a_ParentId, const char *a_Path, 
                             xid *a_GroupId);
XMDF_API xid xfCloseGroup(xid GroupId);

XMDF_API xid xfGetGroupPathsSizeForMeshes(xid FileId, int *Num,int *Maxsize);
XMDF_API xid xfGetGroupPathsForMeshes(xid a_FileId, int a_Num, int a_Maxsize,
                                         char *a_Paths);

XMDF_API xid xfGetGroupPathsSizeForGrids(xid FileId, int *Num,int *Maxsize);
XMDF_API xid xfGetGroupPathsForGrids(xid a_FileId, int a_Num, int a_Maxsize,
                                     char *a_Paths);

XMDF_API xid xfGetGroupPathsSizeForXsecs(xid a_FileId, int *a_Num,int *Maxsize);
XMDF_API xid xfGetGroupPathsForXsecs(xid a_FileId, int a_Num, int a_Maxsize,
                                     char *a_Paths);

XMDF_API xid xfGetGroupPathsSizeForGeomPaths(xid FileId, int *Num,int *Maxsize);
XMDF_API xid xfGetGroupPathsForGeomPaths(xid a_FileId, int a_Num, int a_Maxsize,
                                         char *a_Paths);

XMDF_API xid xfOpenPropertyGroup(xid a_ParentId, xid *a_GroupId);

XMDF_API xid xfGetGroupAbsolutePathSize(xid a_GroupId, int *a_PathLength);
XMDF_API xid xfGetGroupAbsolutePath(xid a_GroupId, int a_PathLength,
                                          char *a_Path);

  /* float and double write types */
XMDF_API hid_t xfGetDoubleType (void);
XMDF_API hid_t xfGetFloatType (void);
XMDF_API void xfSetFloatType (int a_BigEndian);

/* writing/reading reftimes */
XMDF_API xid xfWriteReftime(xid a_Id, double a_Reftime);
XMDF_API xid xfUseReftime(xid a_Id, xmbool *a_bUseReftime);
XMDF_API xid xfReadReftime(xid a_Id, double *a_Reftime);

  /* functions to convert between calendar dates and julian dates */
XMDF_API xid xfCalendarToJulian(xmbool a_bEra, int a_yr, int a_mo, int a_day, 
                        int a_hr, int a_min, int a_sec, double *a_julian);
XMDF_API xid xfJulianToCalendar (xmbool *a_bEra, int *a_yr, int *a_mo, int *a_day,
                        int *a_hr, int *a_min, int *a_sec, double a_julian);

  /* mesh functions */
    /* writing */
XMDF_API xid xfSetNumberOfElements(xid a_Id, int a_nElems);
XMDF_API xid xfSetAllElemsSameType(xid a_Id, int a_Type);
XMDF_API xid xfWriteElemTypes(xid a_Id, int a_nElems, const int *a_Type,
                                 int Compression);

XMDF_API xid xfWriteElemNodeIds(xid a_Id, int a_nElems, int a_nMaxNodes,
                                   int *a_Ids, int Compression);
XMDF_API xid xfSetNumberOfNodes(xid a_Id, int a_nNodes);
  /* x locations must be written first */
XMDF_API xid xfWriteXNodeLocations(xid a_Id, int a_nNodes, double *a_Locs,
                                      int Compression);
XMDF_API xid xfWriteYNodeLocations(xid a_Id, int a_nNodes, double *a_Locs);
XMDF_API xid xfWriteZNodeLocations(xid a_Id, int a_nNodes, double *a_Locs);

  /* reading */
XMDF_API xid xfGetNumberOfElements(xid a_Id, int *a_nElems);
XMDF_API xid xfAreAllElemsSameType(xid a_Id, xmbool *a_Same);
XMDF_API xid xfReadElemTypesSingleValue(xid a_Id, int *a_Type);
XMDF_API xid xfReadElemTypes(xid a_Id, int a_nElems, int *a_Type);

XMDF_API xid xfGetMaxNodesInElem(xid a_Id, int *a_nMaxNodes);
XMDF_API xid xfReadElemNodeIds(xid a_Id, int a_nElems, int a_nMaxNodes,
                                  int *a_Ids);

XMDF_API xid xfGetNumberOfNodes(xid a_Id, int *a_nNodes);
XMDF_API xid xfReadXNodeLocations(xid a_Id, int a_nNodes, double *a_Locs);
XMDF_API xid xfReadYNodeLocations(xid a_Id, int a_nNodes, double *a_Locs);
XMDF_API xid xfReadZNodeLocations(xid a_Id, int a_nNodes, double *a_Locs);

  /* Mesh property groups */
XMDF_API xid xfCreateMeshPropertyGroup (xid a_Id, xid *a_PropId);
XMDF_API xid xfGetMeshPropertyGroup(xid a_Id, xid *a_PropId);
XMDF_API xid xfCreateMeshNodePropertyGroup (xid a_Id, xid *a_PropId);
XMDF_API xid xfGetMeshNodePropertyGroup(xid a_Id, xid *a_PropId);
XMDF_API xid xfCreateMeshElementPropertyGroup (xid a_Id, xid *a_PropId);
XMDF_API xid xfGetMeshElementPropertyGroup(xid a_Id, xid *a_PropId);

  /* Grid functions */
XMDF_API xid xfSetGridType(xid a_Id, int a_GridType);

XMDF_API xid xfSetNumberOfDimensions(xid a_Id, int a_NumDimensions);
XMDF_API xid xfSetExtrusionType(xid a_Id, int a_ExtrudeType);
XMDF_API xid xfSetOrigin(xid a_Id, double a_x, double a_y, double a_z);
XMDF_API xid xfSetOrientation(xid a_Id, int a_Orientation);
XMDF_API xid xfSetBearing(xid a_Id, double a_Bearing);
XMDF_API xid xfSetDip(xid a_Id, double a_Dip);
XMDF_API xid xfSetRoll(xid a_Id, double a_Roll);
XMDF_API xid xfSetComputationalOrigin(xid a_Id, int a_origin);
XMDF_API xid xfSetUDirection(xid a_Id, int a_direction);

XMDF_API xid xfSetNumberCellsInI(xid a_Id, int a_NumI);
XMDF_API xid xfSetNumberCellsInJ(xid a_Id, int a_NumJ);
XMDF_API xid xfSetNumberCellsInK(xid a_Id, int a_NumK);

XMDF_API xid xfSetGridCoordsI(xid a_Id, int a_NumVals, double *a_iValues);
XMDF_API xid xfSetGridCoordsJ(xid a_Id, int a_NumVals, double *a_jValues);
XMDF_API xid xfSetGridCoordsK(xid a_Id, int a_NumVals, double *a_kValues);

XMDF_API xid xfWriteExtrudeLayerData(xid a_Id, int a_NumLayers, int a_NumVals,
                                     double *a_Values);

XMDF_API xid xfGetGridType(xid a_Id, int *a_GridType);
XMDF_API xid xfGetExtrusionType(xid a_Id, int *a_ExtrudeType);
XMDF_API xid xfGetNumberOfDimensions(xid a_Id, int *a_NumDimensions);
XMDF_API xid xfOriginDefined(xid a_Id, xmbool *a_bDefined);
XMDF_API xid xfGetOrigin(xid a_Id, double *a_x, double *a_y, double *a_z);
XMDF_API xid xfGetOrientation(xid a_Id, int *a_Orientation);
XMDF_API xid xfBearingDefined(xid a_Id, xmbool *a_bDefined);
XMDF_API xid xfGetBearing(xid a_Id, double *a_bearing);
XMDF_API xid xfDipDefined(xid a_Id, xmbool *a_bDefined);
XMDF_API xid xfGetDip(xid a_Id, double *a_dip);
XMDF_API xid xfRollDefined(xid a_Id, xmbool *a_bDefined);
XMDF_API xid xfGetRoll(xid a_Id, double *a_Roll);
XMDF_API xid xfComputationalOriginDefined(xid GroupId, xmbool *bDefined);
XMDF_API xid xfGetComputationalOrigin(xid GroupId, int *Origin);
XMDF_API xid xfGetUDirectionDefined(xid GroupId, xmbool *bDefined);
XMDF_API xid xfGetUDirection (xid GroupId, int *Direction);
XMDF_API xid xfGetNumberCellsInI(xid a_Id, int *a_NumI);
XMDF_API xid xfGetNumberCellsInJ(xid a_Id, int *a_NumJ);
XMDF_API xid xfGetNumberCellsInK(xid a_Id, int *a_NumK);

XMDF_API xid xfGetGridCoordsI(xid a_Id, int a_NumVals, double *a_iValues);
XMDF_API xid xfGetGridCoordsJ(xid a_Id, int a_NumVals, double *a_jValues);
XMDF_API xid xfGetGridCoordsK(xid a_Id, int a_NumVals, double *a_kValues);

XMDF_API xid xfGetExtrudeNumLayers(xid a_Id, int *a_NumLayers);
XMDF_API xid xfGetExtrudeValues(xid a_Id, int a_NumVals, double *a_Values);

    /* Grid Property group */
XMDF_API xid xfCreateGridPropertyGroup(xid a_Id, xid *a_PropId);
XMDF_API xid xfGetGridPropertyGroup(xid a_Id, xid *a_PropId);
XMDF_API xid xfCreateGridCellPropertyGroup(xid a_Id, xid *a_PropId);
XMDF_API xid xfGetGridCellPropertyGroup(xid a_Id, xid *a_PropId);
XMDF_API xid xfCreateGridNodePropertyGroup(xid a_Id, xid *a_PropId);
XMDF_API xid xfGetGridNodePropertyGroup(xid a_Id, xid *a_PropId);






  /* Xsec Functions */

XMDF_API xid xfSetNumberOfXSects (xid a_Id, int *a_nXSects, 
                                  int a_compression);
XMDF_API xid xfGetNumberOfXSects (xid a_Id, int *a_nXSects);
XMDF_API xid xfSetCSID (xid a_Id, int a_NumVals, int *a_PropId, 
                        int a_compression);
XMDF_API xid xfGetCSID (xid a_Id,  int *a_NumVals, int *a_PropId);
XMDF_API xid xfSetCSName (xid a_Id, int a_NumVals, int a_StrLen, char *a_PropId);
XMDF_API xid xfGetCSName (xid a_Id, int *a_NumVals, int *a_StrLen, char *a_PropId);
XMDF_API xid xfGetCSNameLen (xid a_Id, int *a_NumVals, int *a_StrLen);
XMDF_API xid xfSetReachName (xid a_Id, int a_NumVals, int a_StrLen, char *a_PropId);
XMDF_API xid xfGetReachName (xid a_Id, int *a_NumVals, int *a_StrLen, char *a_PropId);
XMDF_API xid xfGetReachNameLen (xid a_Id, int *a_NumVals, int *a_StrLen);
XMDF_API xid xfSetTopoID (xid a_Id, int a_NumVals, int a_StrLen, char *a_PropId);
XMDF_API xid xfGetTopoID (xid a_Id, int *a_NumVals, int *a_StrLen, char *a_PropId);
XMDF_API xid xfGetTopoIDLen (xid a_Id, int *a_NumVals, int *a_StrLen);
XMDF_API xid xfSetStation (xid a_Id, int a_NumVals, double *a_PropId, 
                           int a_compression);
XMDF_API xid xfGetStation (xid a_Id, int *a_NumVals, double *a_PropId);
XMDF_API xid xfSetType (xid a_Id, int a_NumVals, int *a_PropId, 
                        int a_compression);
XMDF_API xid xfGetType (xid a_Id, int *a_NumVals, int *a_PropId);
XMDF_API xid xfSetpType (xid a_Id, int a_NumVals, int *a_PropId, 
                         int a_compression);
XMDF_API xid xfGetpType (xid a_Id, int *a_NumVals, int *a_PropId);
XMDF_API xid xfSetpCSDBLink (xid a_Id, int a_NumVals, int *a_PropId, 
                             int a_compression);
XMDF_API xid xfGetpCSDBLink (xid a_Id, int *a_NumVals, int *a_PropId);
XMDF_API xid xfSetNote (xid a_Id, int a_NumVals, int a_StrLen, char *a_PropId);
XMDF_API xid xfGetNote (xid a_Id, int *a_NumVals, int *a_StrLen, char *a_PropId);
XMDF_API xid xfGetNoteLen(xid a_Id, int *a_NumVals, int *a_StrLen);

  /*Geometery */
XMDF_API xid xfSetXSectGeomX (xid a_Id, int a_index, int a_NumVals, 
                             double *a_iValues, int a_compression);
XMDF_API xid xfGetXSectGeomX (xid a_Id, int a_index, int *a_NumVals, 
                             double *a_iValues);
XMDF_API xid xfSetXSectGeomY (xid a_Id, int a_index, int a_NumVals, 
                             double *a_iValues, int a_compression);
XMDF_API xid xfGetXSectGeomY (xid a_Id, int a_index, int *a_NumVals, 
                             double *a_iValues);
XMDF_API xid xfSetXSectGeomD (xid a_Id, int a_index, int a_NumVals, 
                             double *a_iValues, int a_compression);
XMDF_API xid xfGetXSectGeomD (xid a_Id, int a_index, int *a_NumVals, 
                             double *a_iValues);
XMDF_API xid xfSetXSectGeomZ (xid a_Id, int a_index, int a_NumVals, 
                             double *a_iValues, int a_compression);
XMDF_API xid xfGetXSectGeomZ (xid a_Id, int a_index, int *a_NumVals, 
                             double *a_iValues);

  /*XSecs--Line properties */
XMDF_API xid xfSetXSectLinePropFrom (xid a_Id, int a_index, int a_NumVals, 
                                     double *a_PropId, int a_compression);
XMDF_API xid xfGetXSectLinePropFrom (xid a_Id, int a_index, int *a_NumVals, 
                                     double *a_PropId);
XMDF_API xid xfSetXSectLinePropTo (xid a_Id, int a_index, int a_NumVals, 
                                   double *a_PropId, int a_compression);
XMDF_API xid xfGetXSectLinePropTo (xid a_Id, int a_index, int *a_NumVals, 
                                   double *a_PropId);
XMDF_API xid xfSetXSectLinePropType (xid a_Id, int a_index, int a_NumVals, 
                                     int *a_PropId, int a_compression);
XMDF_API xid xfGetXSectLinePropType (xid a_Id, int a_index, int *a_NumVals, 
                                     int *a_PropId);
XMDF_API xid xfSetXSectLinePropIValue (xid a_Id, int a_index, int a_NumVals, 
                                       int *a_PropId, int a_compression);
XMDF_API xid xfGetXSectLinePropIValue (xid a_Id, int a_index, int *a_NumVals, 
                                       int *a_PropId);
XMDF_API xid xfSetXSectLinePropFValue (xid a_Id, int a_index, int a_NumVals, 
                                       double *a_PropId, int a_compression);
XMDF_API xid xfGetXSectLinePropFValue (xid a_Id, int a_index, int *a_NumVals, 
                                       double *a_PropId);

  /*XSec--Line properties--2 */
XMDF_API xid xfSetXSectLinePropID (xid a_Id, int a_NumVals, int *a_PropId, 
                                   int a_compression);
XMDF_API xid xfGetXSectLinePropID (xid a_Id, int *a_NumVals, int *a_PropId);
XMDF_API xid xfSetXSectLinePropName (xid a_Id, int a_NumVals, int a_StrLen, 
                                     char *a_PropId);
XMDF_API xid xfGetXSectLinePropName (xid a_Id, int *a_NumVals, int *a_StrLen, 
                                     char *a_PropId);
XMDF_API xid xfGetXSectLinePropNameLen (xid a_Id, int *a_NumVals, int *a_StrLen);
XMDF_API xid xfSetXSectLinePropDesc (xid a_Id, int a_NumVals, int a_StrLen, 
                                     char *a_PropId);
XMDF_API xid xfGetXSectLinePropDesc (xid a_Id, int *a_NumVals, int *a_StrLen, 
                                     char *a_PropId);
XMDF_API xid xfGetXSectLinePropDescLen (xid a_Id, int *a_NumVals, int *a_StrLen);
XMDF_API xid xfSetXSectLinePropCategory (xid a_Id, int a_NumVals, int *a_PropId, 
                                         int a_compression);
XMDF_API xid xfGetXSectLinePropCategory (xid a_Id, int *a_NumVals, int *a_PropId);
XMDF_API xid xfSetXSectLinePropFDefault (xid a_Id, int a_NumVals, double *a_PropId, 
                                         int a_compression);
XMDF_API xid xfGetXSectLinePropFDefault (xid a_Id, int *a_NumVals, double *a_PropId);
XMDF_API xid xfSetXSectLinePropExclusive (xid a_Id, int a_NumVals, int *a_PropId, 
                                          int a_compression);
XMDF_API xid xfGetXSectLinePropExclusive (xid a_Id, int *a_NumVals, int *a_PropId);

XMDF_API xid xfSetNumberOfLinePropEnumGroup (xid a_Id, int *a_nPropNum, 
                                             int a_compression);
XMDF_API xid xfGetNumberOfLinePropEnumGroup (xid a_Id, int *a_nPropNum);
XMDF_API xid xfSetNumberOfLinePropEnum (xid a_Id, int a_index, int *a_nPropNum, 
                                             int a_compression);
XMDF_API xid xfGetNumberOfLinePropEnum (xid a_Id, int a_index, int *a_nPropNum);

XMDF_API xid xfSetXSectLinePropEnumID (xid a_Id, int a_index, int a_NumVals, 
                                       int *a_PropId, int a_compression);
XMDF_API xid xfGetXSectLinePropEnumID (xid a_Id, int a_index, int *a_NumVals, 
                                       int *a_PropId);
XMDF_API xid xfSetXSectLinePropEnumMatID (xid a_Id, int a_index, int a_NumVals, 
                                          int *a_PropId, int a_compression);
XMDF_API xid xfGetXSectLinePropEnumMatID (xid a_Id, int a_index, int *a_NumVals, 
                                          int *a_PropId);
XMDF_API xid xfSetXSectLinePropEnumName (xid a_Id, int a_index, int a_NumVals, 
                                         int a_StrLen, char *a_PropId);
XMDF_API xid xfGetXSectLinePropEnumName (xid a_Id, int a_index, int *a_NumVals, 
                                         int *a_StrLen, char *a_PropId);
XMDF_API xid xfGetXSectLinePropEnumNameLen (xid a_Id, int a_index, int *a_NumVals, 
                                            int *a_StrLen);
  /*XSec--Point properties--1 */
XMDF_API xid xfSetXSectPointPropMeasure (xid a_Id, int a_index, int a_NumVals, 
                                         double *a_PropId, int a_compression);
XMDF_API xid xfGetXSectPointPropMeasure (xid a_Id, int a_index, int *a_NumVals, 
                                         double *a_PropId);
XMDF_API xid xfSetXSectPointPropType (xid a_Id, int a_index, int a_NumVals, int *a_PropId, 
                                      int a_compression);
XMDF_API xid xfGetXSectPointPropType (xid a_Id, int a_index, int *a_NumVals, int *a_PropId);

  /*XSec--Point Properties--2 */
XMDF_API xid xfSetXSectPointPropID (xid a_Id, int a_NumVals, 
                                          int *a_PropId, int a_compression);
XMDF_API xid xfGetXSectPointPropID (xid a_Id, int *a_NumVals, 
                                          int *a_PropId);
XMDF_API xid xfSetXSectPointPropName (xid a_Id, int a_NumVals, int a_StrLen, 
                                          char *a_PropId);
XMDF_API xid xfGetXSectPointPropName (xid a_Id, int *a_NumVals, int *a_StrLen, 
                                          char *a_PropId);
XMDF_API xid xfGetXSectPointPropNameLen (xid a_Id, int *a_NumVals, int *a_StrLen);
XMDF_API xid xfSetXSectPointPropDesc (xid a_Id, int a_NumVals, int a_StrLen, 
                                          char *a_PropId);
XMDF_API xid xfGetXSectPointPropDesc (xid a_Id, int *a_NumVals, int *a_StrLen, 
                                          char *a_PropId);
XMDF_API xid xfGetXSectPointPropDescLen (xid a_Id, int *a_NumVals, int *a_StrLen);
XMDF_API xid xfSetXSectPointPropExclusive (xid a_Id, int a_NumVals, 
                                           int *a_PropId, int a_compression);
XMDF_API xid xfGetXSectPointPropExclusive (xid a_Id, int *a_NumVals, 
                                           int *a_PropId);
  /*XSec--Profiles */
XMDF_API xid xfSetXSectTopoName (xid a_Id, int a_NumVals, int a_StrLen, char *a_PropId);
XMDF_API xid xfGetXSectTopoName (xid a_Id, int *a_NumVals, int *a_StrLen, char *a_PropId);
XMDF_API xid xfGetXSectTopoNameLen (xid a_Id, int *a_NumVals, int *a_StrLen);
XMDF_API xid xfSetXSectTopoDesc (xid a_Id, int a_NumVals, int a_StrLen, char *a_PropId);
XMDF_API xid xfGetXSectTopoDesc (xid a_Id, int *a_NumVals, int *a_StrLen, char *a_PropId);
XMDF_API xid xfGetXSectTopoDescLen (xid a_Id, int *a_NumVals, int *a_StrLen);

  /* Metadata to XML */
XMDF_API xid xfGetMetadataXML(xid a_Id, char **a_xml);
XMDF_API xid xfGetMetadataProfileXML (xid a_Id, char **a_xml);
XMDF_API xid xfGetMetadataSourceXML (xid a_Id, char **a_xml);
XMDF_API xid xfGetMetadataSpatialXML (xid a_Id, char **a_xml);

  /* Geometric path group */
XMDF_API xid xfCreateGeometricPathGroup(xid a_ParentId, const char *a_Path,
               const char *a_Guid, int a_Compression,
               xid *a_PathGroup, double a_NullVal);
XMDF_API xid xfWriteParticleTimestep(xid a_PathGroup, int a_nDim, double a_Time,
               int a_nPaths, double *a_Locs);

XMDF_API xid xfGetPathNullVal(xid GroupId, double *NullVal);
XMDF_API xid xfGetNumberOfPaths(xid GroupId, int *NumPaths);
XMDF_API xid xfGetNumberOfTimes(xid GroupId, int *NumTimes);
XMDF_API xid xfGetPathDimensionality(xid GroupId, int *NumDims);
XMDF_API xid xfGetPathTimesArray(xid GroupId, int NumTimes, double *Times);
XMDF_API xid xfReadPathLocationsAtTime(xid GroupId, int TimeIndex,
                                      int FirstPathIndex,
                                      int NumIndicies, double *Locs);
XMDF_API xid xfReadPathLocationsForParticle(xid GroupId, int PathIndex, 
                                   int FirstTimeIndex, int NumTimes,
                                   double *Locs);
XMDF_API xid xfReadPathLocationsForParticles(xid GroupId, int NumPaths,
                                    const int *PathIndices,
                                    int FirstTimeIndex, int NumTimes,
                                    double *Locs);

  /* Function to setup evertying necessary to start writing datasets */
XMDF_API xid xfSetupToWriteDatasets(const char *a_Filename,
                           const char *MultiDatasetsGroupPath, 
                           const char *PathInMultiDatasetsGroup,
                           const char *SpatialDataObjectGuid,
                           int OverwriteOptions, xid *FileId, xid *GroupId);
  /* Overload for in-memory file */
XMDF_API xid xfSetupToWriteDatasets2(const char *a_Filename,
                           const char *MultiDatasetsGroupPath, 
                           const char *PathInMultiDatasetsGroup,
                           const char *SpatialDataObjectGuid,
                           int OverwriteOptions, xid *FileId, xid *GroupId,
                           xmbool a_inMemory);

  /* Multi-datasets group */
XMDF_API xid xfCreateMultiDatasetsGroup(xid a_Id, const char *a_Path,
                                        const char *a_Guid, xid *a_MultiId);
XMDF_API xid xfGetGroupPathsSizeForMultiDatasets(xid a_Id, int *Num,
                                                 int *Maxsize);
XMDF_API xid xfGetAllGroupPathsForMultiDatasets(xid a_Id, int a_Num,
                                                int a_Maxsize, char *a_Paths);
XMDF_API xid xfGetDatasetsSdoGuid(xid a_MultiDatasetsGroup, char *a_GUID);
XMDF_API xid xfOpenMultiDatasetsGroup(xid a_Id, xid *DatasetsGroupId);

  /* writing dataset functions */
  /* scalar */
XMDF_API xid xfCreateScalarDataset(xid DatasetGroupId, const char *Path,
                    const char *Units, const char *TimeUnits, int Compression,
                    xid *DatasetId);
XMDF_API xid xfCreateScalarDatasetExtendable(xid DatasetGroupId, const char *Path,
                    const char *Units, const char *TimeUnits, float aFillVal,
                    int Compression, xid *DatasetId);
XMDF_API xid xfExtendScalarDataset(xid DatasetId, int aNewSize);

XMDF_API xid xfWriteScalarTimestep(xid Id, double Time, 
                    int NumValues, const float *Values);
XMDF_API xid xfWriteScalarTimestepMinMax(xid Id, double Time, int NumValues,
                    const float *Values, float Min, float Max);
XMDF_API xid xfSetDatasetNumTimes(xid a_Id, int a_NumTimes);

  /* The following functions can be used together to write scalar timesteps
     in parts. This is useful if you have the data in multiple arrays */
XMDF_API xid xfInitializeScalarTimestep(xid xScalarAId, double dTime,
                    int nValues, float minvalue, float maxvalue, 
                    hsize_t *timestepId);
XMDF_API xid xfWriteScalarTimestepPortion(xid xScalarAId, int timestepId,
                   int numValuesToWrite, int startIndex, const float *a_values);

  /* This can be used to overwrite specified min and max values for the timestep
     which is necessary if you don't know them up front and writing
     datasets in portions. Applies to both scalar and vector timesteps*/
XMDF_API xid xfSetDatasetTimestepMinMax(xid xDatasetId, int timestepId,
                                        float minvalue, float maxvalue);

  /* vector */
XMDF_API xid xfCreateVectorDataset(xid DatasetGroupId, const char *Path,
                    const char *Units, const char *TimeUnits, int Compression,
                    xid *DatasetId);

XMDF_API xid xfWriteVectorTimestep(xid Id, double Time,
                    int NumValues, int NumComponents, const float *Values);
XMDF_API xid xfWriteVectorTimestepMinMax(xid Id, double Time,
                    int NumValues, int NumComponents, const float *Values,
                    float Min, float Max);

  /* The following functions can be used together to write vector timesteps
     in parts. This is useful if you have the data in multiple arrays */
XMDF_API xid xfInitializeVectorTimestep(xid a_Id, double dTime,
                   int nValues, int nComponents, float minvalue, float maxvalue, 
                   hsize_t *timestepId);
XMDF_API xid xfWriteVectorTimestepPortion(xid a_Id, int timestepId,
                    int numValuesToWrite, int nComponentsToWrite,
                    int startIndex, int startComponent, const float *a_values);

  /* activity */
XMDF_API xid xfWriteActivityTimestep (xid a_Id, int a_NumActive,
                                         const unsigned char *a_Active);
XMDF_API xid xfInitializeActivityTimestep (xid a_Id, int a_NumActive,
                                           hsize_t *a_timestepId);
XMDF_API xid xfWriteActivityTimestepPortion(xid a_Id, hsize_t a_timestepId,
               int a_NumValuesToWrite, int a_startIndex,
               const unsigned char *a_activityValues);
/* writing/reading reftimes */
XMDF_API xid xfDatasetReftime(xid a_Id, double a_Reftime);
XMDF_API xid xfUseDatasetReftime(xid a_Id, xmbool *a_bUseReftime);
XMDF_API xid xfReadDatasetReftime(xid a_Id, double *a_Reftime);

/* finding paths to datasets */
XMDF_API xid xfGetScalarDatasetsInfo(xid a_Id, int *a_Number,
                                        int *a_MaxPathLength);
XMDF_API xid xfGetScalarDatasetPaths(xid a_Id, int a_Number,
                                        int a_MaxLength, char *Paths);
XMDF_API xid xfGetStationInfo(xid a_Id, int *a_Number,
                                        int *a_MaxPathLength);
XMDF_API xid xfGetStationPaths(xid a_Id, int a_Number,
                                        int a_MaxLength, char *Paths);
XMDF_API xid xfGetVectorDatasetsInfo(xid a_Id, int *a_Number,
                                        int *a_MaxPathLength);
XMDF_API xid xfGetVectorDatasetPaths(xid a_Id, int a_Number,
                                        int a_MaxLength, char *Paths);

  /* used to change scalar values in an existing dataset */
  /* The file has to be closed and reopened before reading these values */
XMDF_API xid xfChangeScalarValuesTimestepFloat(xid a_Id, int a_TimestepIndex, 
               int a_NumValsToEdit, int *a_Indices, float *a_NewValues);

  /* reading dataset functions */
XMDF_API xid xfGetDatasetNumTimes(xid a_Id, int *a_Numtimes);
XMDF_API xid xfGetDatasetNumVals (xid a_Id, int *a_Numvals);
XMDF_API xid xfGetDatasetNumActive(xid a_Id, int *a_NumActivevals);
/* for vector datasets */
XMDF_API xid xfGetDatasetVecNumComponents(xid a_Id, int *a_NumComponents);

XMDF_API xid xfGetDatasetTimeUnits (xid a_Id, char *Units); /* units 100 long */
XMDF_API xid xfGetDatasetUnits     (xid a_Id, char *Units); /* units 100 long */

XMDF_API xid xfGetVectorDatasetGroupId (xid a_Id);
XMDF_API xid xfGetScalarDatasetGroupId (xid a_Id);

  /* Functions to get entire arrays.  Arrays must already be allocated */
XMDF_API xid xfGetDatasetTimes(xid a_DatasetId, int a_NumTimes,
                               double *a_Times);
XMDF_API xid xfGetDatasetMins(xid a_DatasetId, int a_NumTimes, float *a_Mins);
XMDF_API xid xfGetDatasetMaxs(xid a_DatasetId, int a_NumTimes, float *a_Maxs);
XMDF_API xid xfGetDatasetMinsFloat(xid a_Id,   int a_NumTimes, float  *a_Mins);
XMDF_API xid xfGetDatasetMaxsFloat(xid a_Id,   int a_NumTimes, float  *a_Maxs);
XMDF_API xid xfGetDatasetMinsDouble(xid a_Id,  int a_NumTimes, double *a_Mins);
XMDF_API xid xfGetDatasetMaxsDouble(xid a_Id,  int a_NumTimes, double *a_Maxs);
/* XMDF_API xid xfGetDatasetActivity(xid a_DatasetId, unsigned char *a_Active); */

  /* Functions to get portions of datasets.  Arrays must already be allocated */
XMDF_API xid xfReadScalarValuesTimestep(xid a_DatasetId, 
                    int a_TimestepIndex, int a_NumVals, float *a_Values);
XMDF_API xid xfReadScalarValuesTimestepFloat(xid a_DatasetId, 
                    int a_TimestepIndex, int a_NumVals, float *a_Values);
XMDF_API xid xfReadScalarValuesTimestepFloatPortion (xid a_Id, 
                    int a_TimestepIndex, int a_Start, int a_NumVals,
                    float *a_Values);
XMDF_API xid xfReadScalarValuesTimestepDouble(xid a_Id, 
                  int a_TimestepIndex, int a_NumVals, double *a_Values);
XMDF_API xid xfReadScalarValuesTimestepDoublePortion(xid a_Id, 
                  int a_TimestepIndex, int a_Start, int a_NumVals,
                  double *a_Values);
XMDF_API xid xfReadScalarValuesTimestepInt(xid a_Id, 
                  int a_TimestepIndex, int a_NumVals, int *a_Values);
XMDF_API xid xfReadScalarValuesAtIndex(xid a_DatasetId,
                    int a_Index, int a_FirstTimestep, int a_NumTimes,
                    float *a_Values);
XMDF_API xid xfReadScalarValuesAtIndexFloat(xid a_DatasetId,
                    int a_Index, int a_FirstTimestep, int a_NumTimes,
                    float *a_Values);
XMDF_API xid xfReadScalarValuesAtIndexDouble(xid a_Id, 
                    int a_Index, int a_FirstTime, int a_NumTimes,
                    double *a_Values);
XMDF_API xid xfReadScalarValuesAtIndexInt(xid a_Id,
                    int a_Index, int a_FirstTime, int a_NumTimes,
                    int *a_Values);
XMDF_API xid xfReadScalarValuesAtIndices(xid a_Id, 
                    int a_nIndices, const int *a_Indices, int a_FirstTime,
                    int a_NumTimes, float *a_Values);
XMDF_API xid xfReadScalarValuesAtIndicesFloat(xid a_Id, 
                    int a_nIndices, const int *a_Indices, int a_FirstTime,
                    int a_NumTimes, float *a_Values);
XMDF_API xid xfReadVectorValuesTimestep(xid a_DatasetId,
                    int a_TimestepIndex, int a_NumVals, int a_NumComponents,
                    float *a_Values);
XMDF_API xid xfReadVectorValuesTimestepFloat(xid a_Id, 
                    int a_TimestepIndex, int a_NumVals, int a_NumComponents, 
                    float *a_Values);
XMDF_API xid xfReadVectorValuesTimestepFloatPortion(xid a_Id, 
                    int a_TimestepIndex, int a_Start, int a_NumVals,
                    int a_NumComponents, float *a_Values);
XMDF_API xid xfReadVectorValuesAtIndex(xid a_Id, int a_Index, 
                    int a_FirstTime, int a_NumTimes, int a_NumComponents,
                    float *a_Values);
XMDF_API xid xfReadVectorValuesTimestepDoublePortion (xid a_Id, 
                    int a_TimestepIndex, int a_Start, int a_NumVals,
                    int a_NumComponents, double *a_Values);
XMDF_API xid xfReadVectorValuesAtIndexFloat(xid a_Id, int a_Index,
                    int a_FirstTime, int a_NumTimes, int a_NumComponents,
                    float *a_Values);
XMDF_API xid xfReadVectorValuesTimestepDouble(xid a_Id, 
                    int a_TimestepIndex, int a_NumVals, int a_NumComponents, 
                    double *a_Values);
XMDF_API xid xfReadVectorValuesAtIndices(xid a_Id, 
                    int a_nIndices, const int *a_Indices, int a_FirstTime,
                    int a_NumTimes, int a_NumComponents, float *a_Values);
XMDF_API xid xfReadVectorValuesAtIndicesFloat(xid a_Id, 
                    int a_nIndices, const int *a_Indices, int a_FirstTime,
                    int a_NumTimes, int a_NumComponents, float *a_Values);

XMDF_API xid xfReadActivityTimestep(xid a_Id, int Index, int NumActive,
                    xmbool *Active);
XMDF_API xid xfReadActivityValuesAtIndex(xid a_Id, int a_Index,
                    int a_FirstTime, int a_NumTimes, xmbool *a_Values);
 
XMDF_API hssize_t xfGetParallelNumValuesToRead();
XMDF_API void xfSetParallelNumValuesToRead(hssize_t a);

  /* reserved dataset property functions */
XMDF_API xid xfScalarDataLocation(xid a_Id, int a_DataLoc);
XMDF_API xid xfVector2DDataLocations(xid a_Id, int a_DataLocI, int a_DataLocJ);
XMDF_API xid xfVector3DDataLocations(xid a_Id, int a_DataLocI, int a_DataLocJ,
                                     int a_DataLocK);
XMDF_API xid xfGetScalarDataLocation(xid a_Id, int *a_DataLoc);
XMDF_API xid xfGetVector2DDataLocations(xid a_Id, int *a_DataLocI,
                               int *a_DataLocJ);
XMDF_API xid xfGetVector3DDataLocations(xid a_Id, int *a_DataLocI,
                               int *a_DataLocJ, int *a_DataLocK);
XMDF_API xid xfVectorsInLocalCoords (xid a_Id);
XMDF_API xid xfAreVectorsInLocalCoords (xid a_Id, int *a_LocalCoords);

/* Coordinate conversions */

  // The method of storing projections changed after version 1.5 of XMDF.
  // This function will return 1 if reading a file in 1.5,
  // and 2 if reading a file after 1.5
XMDF_API xid xfGetCoordVersion(xid a_CoordId, int *a_version);

  // These functions are used for all versions of XMDF
XMDF_API xid xfGetHorizUnits(xid a_CoordId, int *a_val);
XMDF_API xid xfGetVertDatum(xid a_CoordId, int *a_val);
XMDF_API xid xfGetVertUnits(xid a_CoordId, int *a_val);

XMDF_API xid xfSetHorizUnits(xid a_CoordId, int a_val);
XMDF_API xid xfSetVertDatum(xid a_CoordId, int a_val);
XMDF_API xid xfSetVertUnits(xid a_CoordId, int a_val);

  // These apply to the new system for storing projection information
  // XMDF versions after 1.5

  // Read whether we are using local coordinates or a defined projection
XMDF_API xid xfGetUsingLocalCoords(xid a_CoordId, xmbool *a_localCoords);
XMDF_API xid xfSetUsingLocalCoords(xid a_CoordId, xmbool *a_localCoords);

  // If we are using a defined projection it is stored in the 
  // well-known text format here (don't call if using local)
XMDF_API xid xfSetWKT (xid a_CoordId, const char * a_wkt);
XMDF_API xid xfGetWKTStringSize(xid a_CoordId, int *a_size);
  // a_wkt should already be allocated to the string size
XMDF_API xid xfGetWKT (xid a_CoordId, char *a_wkt);

  // These apply to the old system for storing projection information
  // XMDF version 1.5 and before.
XMDF_API xid xfGetHorizDatum(xid a_CoordId, int *a_val);
XMDF_API xid xfGetHorizUnits(xid a_CoordId, int *a_val);
XMDF_API xid xfGetVertDatum(xid a_CoordId, int *a_val);
XMDF_API xid xfGetVertUnits(xid a_CoordId, int *a_val);
XMDF_API xid xfGetLat(xid a_CoordId, int *a_val);
XMDF_API xid xfGetLon(xid a_CoordId, int *a_val);
XMDF_API xid xfGetUTMZone(xid a_CoordId, int *a_val);
XMDF_API xid xfGetSPCZone(xid a_CoordId, int *a_val);
XMDF_API xid xfGetHPGNArea(xid a_CoordId, int *a_val);
XMDF_API xid xfGetCPPLat(xid a_CoordId, double *a_val);
XMDF_API xid xfGetCPPLon(xid a_CoordId, double *a_val);
XMDF_API xid xfGetEllipse(xid a_CoordId, int *a_val);
XMDF_API xid xfGetMajorR(xid a_CoordId, double *a_val);
XMDF_API xid xfGetMinorR(xid a_CoordId, double *a_val);

XMDF_API xid xfSetHorizDatum(xid a_CoordId, int a_val);
XMDF_API xid xfSetLat(xid a_CoordId, int a_val);
XMDF_API xid xfSetLon(xid a_CoordId, int a_val);
XMDF_API xid xfSetUTMZone(xid a_CoordId, int a_val);
XMDF_API xid xfSetSPCZone(xid a_CoordId, int a_val);
XMDF_API xid xfSetHPGNArea(xid a_CoordId, int a_val);
XMDF_API xid xfSetCPPLat(xid a_CoordId, double a_val);
XMDF_API xid xfSetCPPLon(xid a_CoordId, double a_val);
XMDF_API xid xfSetEllipse(xid a_CoordId, int a_val);
XMDF_API xid xfSetMajorR(xid a_CoordId, double a_val);
XMDF_API xid xfSetMinorR(xid a_CoordId, double a_val);


#ifdef __cplusplus
}

#define EMRL_XMDF

/******************************************************************************
 * CLASS     ECreateGroup
 * PURPOSE   calls xfiCreateGroup, for EMRL only!
 * NOTES     
 ******************************************************************************/
class ECreateGroup
{
  public:
    ECreateGroup() {;}
    ~ECreateGroup() {;}

    int CallXfiCreateGroup(xid a_Id, const char *Path, 
                           const char *a_GroupType);
};
#endif


#endif
