
/* Xmdf.c : Defines the entry point for the DLL application. */ 

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

#include "xmdf/Xmdf.h"
#include "xmdf/xmdf_private.h"
#include "xmdf/xmdf_timestep.h"
#include "hdf5.h"
#include "xmdf/ErrorDefinitions.h"
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#ifdef EMRL
#include <windows.h>
#else
#ifdef PCAPI
#include <windows.h>
BOOL WINAPI DllMain( HANDLE hModule, 
                       DWORD  ul_reason_for_call, 
                       LPVOID lpReserved
           )
{
    switch (ul_reason_for_call)
  {
    case DLL_PROCESS_ATTACH:
    case DLL_THREAD_ATTACH:
    case DLL_THREAD_DETACH:
    case DLL_PROCESS_DETACH:
      break;
    }
    return XTRUE;
}
#endif
#endif

#define DBG Dbg(__FILE__,__LINE__);
#extern void Dbg();

static xid xfiCreateFile (const char * a_Filename, xid *Id, xmbool a_Overwrite,
                          xmbool a_inMemory);
static xid xfiInitializeScalarTimestep(xid xScalarAId, double dTime,
              int nValues, int a_DataType, double minvalue, double maxvalue,
              hsize_t *timestepId);
static xid xfiWriteScalarTimestepPortion (xid a_id, hsize_t a_timeId,
              int a_NumValues, int a_startIndex, int a_DataType,
              const void *a_values);
static xid xfiSetDatasetTimestepMinMax(xid xDatasetId, int timestepId,
                  int a_DataType, double minvalue, double maxvalue);
static xid xfiInitializeVectorTimestep (xid a_Id, double a_Time, 
              int a_NumValues, int a_NumComponents, int a_DataType,
              double a_Min, double a_Max, hsize_t *a_timestepId);
static xid xfiWriteVectorTimestepPortion (xid a_id, hsize_t a_timeId,
              int a_NumValuesToWrite, int a_NumComponentsToWrite, 
              int a_startIndex, int a_startComponent, int a_DataType,
              const void *a_values);

static xid xfiWriteScalarTimestep(xid a_Id, double a_Time, int a_NumValues,
                                    const void *a_Values, int a_DataType);
static xid xfiWriteScalarTimestepMinMax(xid a_Id, double a_Time, int a_NumValues,
                                          const void *a_Values, int a_DataType,
                                          double a_Min, double a_Max);
static xid xfiWriteVectorTimestep(xid a_Id, double a_Time, 
                                    int a_NumValues, int a_NumComponents, 
                                    const void *a_Values, int a_DataType);
static xid xfiWriteVectorTimestepMinMax(xid a_Id, double a_Time, 
                                          int a_NumValues, int a_NumComponents,
                                          const void *a_Values, int a_DataType,
                                          double a_Min, double a_Max);
static xid xfiReadScalarValuesTimestepPortion(xid a_Id, int a_TimestepIndex,
                                                int a_Start, int a_NumVals,
                                                double *a_dValues, float *a_fValues);
static xid xfiReadVectorValuesTimestepPortion(xid a_Id, int a_TimestepIndex, 
                                                int a_Start, int a_NumVals,
                                                int a_NumComponents, 
                                                double *a_dValues, float *a_fValues);
static xid xfiReadWriteScalarValuesFloatAtIndices (xid a_Id, 
                    ReadWrite_enum a_readWrite,
                    int a_nIndices, const int *a_Indices, int a_FirstTime,
                    int a_NumTimes, float *a_Values);
static void xfiGetParIndexFromIndex(xid a_Id, int a_TimestepIndex,
                             int a_index, int *a_parindex);
static xid xfiReadWriteVectorValuesFloatAtIndices(xid a_Id, 
                         ReadWrite_enum a_readWrite,
                         int a_nIndices, const int *a_Indices, int a_FirstTime,
                         int a_NumTimes, int a_nComponents, float *a_Values);
xid xfpCreateDset (xid a_Id, const char *a_Name,
                const XDatasetParams *a_Params, hid_t a_hdfType, xid *a_DsetId);
xid xftExtendDset2DFirstDim (xid a_Id, const char *a_Name, 
                 hsize_t a_NumNew, hsize_t a_NumDim2, hid_t a_incomingHdfType);
xid xfpWriteDset2DPortion (xid a_Id, const char *a_Name, 
                hid_t hdfTypeInMemory, hsize_t a_Start1, int a_Num1, int a_Start2,
                int a_Num2, const void *a_Array);
xid xftExtendDset3DFirstDim (xid a_Id, const char *a_Name, 
                                   hsize_t a_NumNew, hsize_t a_NumDim2,
                                   hsize_t a_NumDim3, hid_t a_hdfType);
xid xfpWriteDset3DPortion (xid a_Id, const char *a_Name, 
                hid_t hdfTypeInMemory, hsize_t a_Start1, int a_Num1, int a_Start2,
                int a_Num2, int a_Start3, int a_Num3, const void *a_Array);

void Dbg(char const*const a_file, int a_line)
{
  FILE* fp = fopen("c:\\dbg.txt", "a");
  fprintf(fp, "%s:%d\n", a_file, a_line);
  fclose(fp);
}

/*!
 * \ingroup files
 * xfCreateFile()
 * \brief Creates a file to read using the Xmdf library.
 * \param a_Filename Takes the path and name of the XMDF file, e.g. "input.h5"
 * \param Id Returns a handle to the file.
 * \param a_Overwrite Takes !0 to overwrite the file or 0 otherwise.
 * \return A negative value upon failure.
 *
 * FORTRAN: SUBROUTINE XF_CREATE_FILE(Name, Overwrite, FileId, Error)
 * - CHARACTER(LEN=*), INTENT(IN):: Name
 * - LOGICAL, INTENT(IN)       ::   Overwrite
 * - INTEGER(XID), INTENT(OUT) ::   FileId
 * - INTEGER(8), INTENT(OUT)      ::   Error
 **********************************************************************/
 XMDF_API xid xfCreateFile (const char * a_Filename, xid *Id, xmbool a_Overwrite)
{
  return xfiCreateFile(a_Filename, Id, a_Overwrite, XFALSE);
} /* xfCreateFile*/
/*!
 * \ingroup files
 * xfCreateFile()
 * \brief Creates an "in memory file" to read using the Xmdf library.
 * \brief This is useful for creating/working with data using Xmdf that is fast
 * \brief but temporary
 * \param a_Filename Takes the path and name of the XMDF file, e.g. "input.h5"
 * \param Id Returns a handle to the file.
 * \return A negative value upon failure.
 *
 * FORTRAN: SUBROUTINE XF_CREATE_FILE(Name, Overwrite, FileId, Error)
 * - CHARACTER(LEN=*), INTENT(IN):: Name
 * - LOGICAL, INTENT(IN)       ::   Overwrite
 * - INTEGER(XID), INTENT(OUT) ::   FileId
 * - INTEGER(8), INTENT(OUT)      ::   Error
 **********************************************************************/
XMDF_API xid xfCreateInMemoryFile (const char *a_File, xid *Id)
{
  return xfiCreateFile(a_File, Id, XTRUE, XTRUE);
} /* xfCreateInMemoryFile */
/*-----------------------------------------------------------------------------*/
/* FUNCTION  xfOpenFile*/
/*! PURPOSE:   Opens a file to read using the Xmdf library.
 * NOTES     Also checks to make sure the file is an HDF5 file.
 *           If you attempt to open a file that is already open, you must open
 *           it as read only. */
/*-----------------------------------------------------------------------------*/
XMDF_API xid xfOpenFile(const char * a_Filename, xid *Id, xmbool a_ReadOnly)
{
  htri_t    IsHdf5;
  char *  Filetype = NULL;
  int       status;
  unsigned  majnum, minnum, relnum;
  hid_t     accessProperties = 0;

  *Id = -1;
  /* Initialize the double and float types*/
  if (!xfpGetInitialized()) {
    xfpInitialize();
  }

  /* clear the error stack*/
  xfpClearErrors();
  /* set the XMDF error handler*/
  H5Eset_auto1(xfpHDF5ErrorHandler, NULL);

  status = H5get_libversion(&majnum, &minnum, &relnum);

  /* make sure file is an HDF5 file*/
  IsHdf5 = H5Fis_hdf5(a_Filename);
  if (IsHdf5 <= 0) {
    /* if we are not an HDF5 file return a negative number*/
    return ERROR_FILE_NOT_HDF5;
  }

  accessProperties = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fclose_degree(accessProperties, H5F_CLOSE_STRONG);

  if (a_ReadOnly) {
    *Id = H5Fopen(a_Filename, H5F_ACC_RDONLY, accessProperties);
  }
  else {
    *Id = H5Fopen(a_Filename, H5F_ACC_RDWR, accessProperties);
  }
  H5Pclose(accessProperties);


  /* check to make sure file is an Xmdf*/
  status = xfpReadDatasetString(*Id, FILE_TYPE_ID, &Filetype);
  if (status < 0 || strcmp(Filetype, FILE_TYPE_XMDF) != 0) {
    H5Fclose(*Id);
    free(Filetype);
    return ERROR_FILE_NOT_XMDF;
  }

  if (Filetype)
    free(Filetype);
  return *Id;
} /* xfOpenFile*/
/*-----------------------------------------------------------------------------*/
/* FUNCTION  xfCloseFile*/
/*! PURPOSE:   Closes a file
 * NOTES     */
/*-----------------------------------------------------------------------------*/
XMDF_API hid_t xfCloseFile(hid_t a_Id)
{ 
  return H5Fclose(a_Id);
} /* xfCloseFile*/
/*-----------------------------------------------------------------------------*/
/* FUNCTION  xfGetNumErrorMessages*/
/*! PURPOSE:   Get the number of messages on the error stack
 * NOTES     */
/*-----------------------------------------------------------------------------*/
XMDF_API xid xfGetNumErrorMessages (int *a_Num)
{
  return xfpGetNumErrorMessages(a_Num);
} /* xfGetNumErrorMessages*/
/*-----------------------------------------------------------------------------*/
/* FUNCTION  xfGetErrorMessages*/
/*! PURPOSE:
 * NOTES     */
/*-----------------------------------------------------------------------------*/
XMDF_API xid xfGetErrorMessages (int a_Num, char **a_Errors)
{
  return xfpGetErrorMessages(a_Num, a_Errors);
} /* xfGetErrorMessages*/
/******************************************************************************/
/* FUNCTION  xfGetNumOpenIdentifiers */
/*! PURPOSE:   Get the number of open identifiers in a file.
 *  - NOTES: */
/******************************************************************************/
XMDF_API xid xfGetNumOpenIdentifiers(xid a_Id, int *a_Num)
{ 
  *a_Num = H5Fget_obj_count(a_Id, H5F_OBJ_ALL);
  return *a_Num;
} /* xfGetNumOpenIdentifiers */
/******************************************************************************
 * FUNCTION  xfGetOpenIdentifiersInfo
 *//*! PURPOSE:   Get the information for open identifiers in a file.
 * - NOTES:     a_Info should already be allocated to the correct number. */
/*******************2***********************************************************/
XMDF_API xid xfGetOpenIdentifiersInfo(xid a_Id, int a_Num,
                                      IdentifierInfo *a_Info)
{
  int    NumOpen;
  hid_t  *Ids;
  int    i, status;

  xfGetNumOpenIdentifiers(a_Id, &NumOpen);
  if (NumOpen != a_Num) {
    /* ERROR number is not correct. */
    return -2;
  }

  /* get the object ids. */
  Ids = (hid_t*)malloc(a_Num*sizeof(hid_t));
  status = H5Fget_obj_ids(a_Id, H5F_OBJ_ALL, a_Num, Ids);
  if (status < 0) {
    free(Ids);
    return status;
  }

  for (i = 0; i < a_Num; i++) {
    /* copy id into structure */
    a_Info[i].id = Ids[i];
    /* fill in the object names and type */
    status = H5Iget_name(a_Info[i].id, a_Info[i].name, MAX_ID_NAME);
    if (status < 0) {
      free(Ids);
      return status;
    }
    a_Info[i].type = H5Iget_type(a_Info[i].id);
    if (a_Info[i].type == H5I_BADID) {
      free(Ids);
      return -3;
    }
  }

  /* free the ids */
  free(Ids);
  return 1;
} /* xfGetOpenIdentifiersInfo */

/* ----------------------------------------------------------------------------- */
/* FUNCTION  xfGetLibraryVersion */
/*! PURPOSE:   Returns the current version of XMDF
 * NOTES      */
/*------------------------------------------------------------------------------- */
XMDF_API xid xfGetLibraryVersion(float *a_Version)
{
  *a_Version = (float)XMDF_VERSION;
  return 1;
} /* xfGetLibraryVersion */
/*-----------------------------------------------------------------------------*/
/* FUNCTION  xfGetLibraryVersionFile */
/*! PURPOSE:   Obtain the version of XMDF from the library file
 * NOTES     */
/*------------------------------------------------------------------------------- */
XMDF_API xid xfGetLibraryVersionFile(xid a_File, float *a_Version)
{
  int  status;

  status = xfReadPropertyFloat(a_File, FILE_VERSION, 1, a_Version);

  return status;
} /* xfGetLibraryVersionFile */
/*-----------------------------------------------------------------------------*/
/* FUNCTION  xfCreateCoordinateGroup*/
/*! PURPOSE:   Creates a group to store coordinate information in
 * NOTES*/
/*-----------------------------------------------------------------------------*/
XMDF_API xid xfCreateCoordinateGroup (xid a_ParentId, xid *a_ChildId)
{
  int version = 2;
  int error;

  /* Create group with name "Coordinates"*/
  *a_ChildId = xfpCreateGroup(a_ParentId, GROUP_TYPE_COORDS, GROUP_TYPE_COORDS);

  error = xfpWriteAttributeInt(*a_ChildId, COORD_ATT_VERSION, 1, &version);

  if (*a_ChildId <= 0)
  {
    error = 0;
  }
  else
    if (error == 0)
    {
      error = 1;
    }
  
  return error;
} /* xfCreateCoordinateGroup*/
/*-----------------------------------------------------------------------------*/
/* FUNCTION  xfOpenCoordinateGroup*/
/*! PURPOSE:   Opens a coordinate grou
 * NOTES*/
/*-----------------------------------------------------------------------------*/
XMDF_API xid xfOpenCoordinateGroup (xid a_ParentId, xid *a_ChildId)
{
  /* Open "Coordinates" Group*/
  xfOpenGroup(a_ParentId, GROUP_TYPE_COORDS, a_ChildId);

  return *a_ChildId;
} /* xfOpenCoordinateGroup*/
/******************************************************************************
 * FUNCTION  xfCreatePropertyGroup
 *//*! PURPOSE:   Creates a group to store attributes in
 * NOTES     *//*
 ******************2***********************************************************/
XMDF_API xid xfCreatePropertyGroup(xid a_ParentId, xid *Id)
{ 
  *Id = xfpCreateGroup(a_ParentId,GROUP_TYPE_PROPERTIES, GROUP_TYPE_PROPERTIES);

  return *Id;
} /* xfCreatePropertyGroup */
/******************************************************************************
 * FUNCTION  xfWritePropertyString */
 /*! PURPOSE:   writes a string attribute (really a dataset) to the folder
  * - NOTES:     */
/******************************************************************************/
XMDF_API xid xfWritePropertyString (xid a_Id, const char *a_Name,
                  int a_Number, int a_StringLength, const char *a_Attributes)
{ 
  return xfpWriteDatasetStrings(a_Id, a_Name, a_Number, a_StringLength,
                                 a_Attributes);
} /* xfWritePropertyString */
/*! --------------------------------------------------------------------------*/
/* FUNCTION  xfWritePropertyUnsignedInt */
/*! PURPOSE:   writes a unsigned INTEGER(8) attribute (really a dataset) to the
 *            folder
 *  - NOTES:     */
/*------------------------------------------------------------------------------ */
XMDF_API xid xfWritePropertyUnsignedInt (xid a_Id, const char *a_Name, 
                   int a_Number, const unsigned int *a_Properties, 
                   int a_Compression)
{
  XDatasetParams   Params;
  /* Chunk if we are using compression */
  xmbool            bChunked = (a_Compression > 0);
  int              nChunksize;
  int              error;

  if (a_Compression > 0 && a_Number < XF_DATASET_MIN_COMPRESS_SIZE) {
    bChunked = XFALSE;
    a_Compression = NONE;
  }

  /* only do rank 1 arrays for now */
  xfpDsetParamsInit(&Params, 1, bChunked, a_Compression);

  /* if chunking chunk into XF_PROP_CHUNK_NUMBER pieces */
  /* But require chunksize to be between XF_PROP_CHUNK_MIN and XF_PROP_CHUNK_MAX */
  nChunksize = a_Number / XF_PROP_CHUNK_NUMBER;
  nChunksize = Xmax(nChunksize, XF_PROP_CHUNK_MIN);
  nChunksize = Xmin(nChunksize, XF_PROP_CHUNK_MAX);
  nChunksize = Xmin(nChunksize, a_Number);
  xfpDsetParamsSetSizes(&Params, 0, a_Number, a_Number, nChunksize);

  /* write the array */
  error = xfpWriteDsetUInt(a_Id, a_Name, H5T_NATIVE_UINT, &Params,
                            a_Properties);
  xfpDsetParamsDestroy(&Params);

  return error;
} /* xfWritePropertyUnsignedInt */
/* --------------------------------------------------------------------------*/
/* FUNCTION  xfWritePropertyInt */
/*! PURPOSE:
 * - NOTES:     */
/*------------------------------------------------------------------------------ */
XMDF_API xid xfWritePropertyInt (xid a_Id, const char *a_Name, 
                   int a_Number, const int *a_Properties, 
                   int a_Compression)
{
  XDatasetParams   Params;
  /* Chunk if we are using compression */
  xmbool            bChunked = (a_Compression > 0);
  int              nChunksize;
  int              error;

  if (a_Compression > 0 && a_Number < XF_DATASET_MIN_COMPRESS_SIZE) {
    bChunked = XFALSE;
    a_Compression = NONE;
  }

  /* only do rank 1 arrays for now */
  xfpDsetParamsInit(&Params, 1, bChunked, a_Compression);

  /* if chunking chunk into XF_PROP_CHUNK_NUMBER pieces */
  /* But require chunksize to be between XF_PROP_CHUNK_MIN and XF_PROP_CHUNK_MAX */
  nChunksize = a_Number / XF_PROP_CHUNK_NUMBER;
  nChunksize = Xmax(nChunksize, XF_PROP_CHUNK_MIN);
  nChunksize = Xmin(nChunksize, XF_PROP_CHUNK_MAX);
  nChunksize = Xmin(nChunksize, a_Number);
  xfpDsetParamsSetSizes(&Params, 0, a_Number, a_Number, nChunksize);

  /* write the array */
  error = xfpWriteDsetInt(a_Id, a_Name, H5T_NATIVE_INT, &Params,
                            a_Properties);
  xfpDsetParamsDestroy(&Params);

  return error;
} /* xfWritePropertyInt */
/* --------------------------------------------------------------------------*/
/* FUNCTION  xfWritePropertyDouble */
/*! PURPOSE:
 * - NOTES:     */
/*------------------------------------------------------------------------------ */
XMDF_API xid xfWritePropertyDouble (xid a_Id, const char *a_Name, 
                   int a_Number, const double *a_Properties, 
                   int a_Compression)
{
  XDatasetParams   Params;
  /* Chunk if we are using compression */
  xmbool            bChunked = (a_Compression > 0);
  int              nChunksize;
  int              error;

  if (a_Compression > 0 && a_Number < XF_DATASET_MIN_COMPRESS_SIZE) {
    bChunked = XFALSE;
    a_Compression = NONE;
  }

  /* only do rank 1 arrays for now */
  xfpDsetParamsInit(&Params, 1, bChunked, a_Compression);

  /* if chunking chunk into XF_PROP_CHUNK_NUMBER pieces */
  /* But require chunksize to be between XF_PROP_CHUNK_MIN and XF_PROP_CHUNK_MAX */
  nChunksize = a_Number / XF_PROP_CHUNK_NUMBER;
  nChunksize = Xmax(nChunksize, XF_PROP_CHUNK_MIN);
  nChunksize = Xmin(nChunksize, XF_PROP_CHUNK_MAX);
  nChunksize = Xmin(nChunksize, a_Number);
  xfpDsetParamsSetSizes(&Params, 0, a_Number, a_Number, nChunksize);

  /* write the array */
  error = xfpWriteDsetDouble(a_Id, a_Name, &Params, a_Properties);
  xfpDsetParamsDestroy(&Params);

  return error;
} /* xfWritePropertyDouble */
/* -------------------------------------------------------------------------- */
/* FUNCTION  xfWritePropertyFloat */
/*! PURPOSE:
 * - NOTES:     */
/*------------------------------------------------------------------------------ */
XMDF_API xid xfWritePropertyFloat (xid a_Id, const char *a_Name, 
                   int a_Number, const float *a_Properties, 
                   int a_Compression)
{
  XDatasetParams   Params;
  /* Chunk if we are using compression */
  xmbool            bChunked = (a_Compression > 0);
  int              nChunksize;
  int              error;

  if (a_Compression > 0 && a_Number < XF_DATASET_MIN_COMPRESS_SIZE) {
    bChunked = XFALSE;
    a_Compression = NONE;
  }

  /* only do rank 1 arrays for now */
  xfpDsetParamsInit(&Params, 1, bChunked, a_Compression);

  /* if chunking chunk into XF_PROP_CHUNK_NUMBER pieces */
  /* But require chunksize to be between XF_PROP_CHUNK_MIN and XF_PROP_CHUNK_MAX */
  nChunksize = a_Number / XF_PROP_CHUNK_NUMBER;
  nChunksize = Xmax(nChunksize, XF_PROP_CHUNK_MIN);
  nChunksize = Xmin(nChunksize, XF_PROP_CHUNK_MAX);
  nChunksize = Xmin(nChunksize, a_Number);
  xfpDsetParamsSetSizes(&Params, 0, a_Number, a_Number, nChunksize);

  /* write the array */
  error = xfpWriteDsetFloat(a_Id, a_Name, &Params, a_Properties);
  xfpDsetParamsDestroy(&Params);

  return error;
} /* xfWritePropertyFloat */
/******************************************************************************/
/*   FUNCTION  xfDoesPropertyWithNameExist */
/*! PURPOSE:   Looks to see if an Property with a given name exists
 *   - NOTES:     */
/******************************************************************************/
XMDF_API xid xfDoesPropertyWithNameExist(xid a_Id, const char *a_Name,
                                             xmbool *a_Exists)
{
  xid  Id;
  herr_t (*old_func)(void*);
  void    *old_client_data;

  H5Eget_auto1(&old_func, &old_client_data);

  /* Turn off error handling */
  H5Eset_auto1(NULL, NULL);

  /* try to open the dataset */
  Id = H5Dopen1(a_Id, a_Name);

  /* Restore previous error handler */
  H5Eset_auto1(old_func, old_client_data);

  if (Id > 0) {
    H5Dclose(Id);
    *a_Exists = XTRUE;
  }
  else {
    *a_Exists = XFALSE;
  }
  return XTRUE;
} /* xfDoesPropertyWithNameExist */
/******************************************************************************/
/*   FUNCTION  xfGetPropertyNumber */
/*! PURPOSE:   Gets the number of items in a property array
 *   - NOTES:     */
/******************************************************************************/
XMDF_API xid xfGetPropertyNumber(xid a_Id, const char *a_Name, int *a_Number)
{
  hid_t   DsetId, DataspaceId;
  herr_t  status;
  htri_t  IsSimple;
  int     Rank;
  hsize_t *Dims = NULL, *MaxDims = NULL;

  DsetId = H5Dopen1(a_Id, a_Name);
  if (DsetId < 0) {
   return DsetId;
  }

  /* Get the dataspace */
  DataspaceId = H5Dget_space(DsetId);

  /* Dataspace must be simply defined */
  IsSimple = H5Sis_simple(DataspaceId);
  if (!IsSimple) {
    H5Sclose(DataspaceId);
    H5Dclose(DsetId);
    return ERROR_ATTRIBUTE_NOT_SUPPORTED;
  }

  /* Get the rank and dimensions for simple dataspace */
  status = xfpGetSimpleDataspaceInfo(DataspaceId, &Rank, &Dims, &MaxDims);
  if (status < 0) {
    H5Sclose(DataspaceId);
    H5Dclose(DsetId);
    return ERROR_ATTRIBUTE_NOT_SUPPORTED;
  }

  /* See how many items we have */
  *a_Number = (int)Dims[0];

  free(Dims);
  free(MaxDims);
  H5Sclose(DataspaceId);
  H5Dclose(DsetId);

  return XTRUE;
} /* xfGetPropertyNumber */
/******************************************************************************/
/*   FUNCTION  xfGetPropertyStringLength */
/*! PURPOSE:   Gets the string length from a dataset
 *   - NOTES:     */
/******************************************************************************/
XMDF_API xid xfGetPropertyStringLength(xid a_Id, const char *a_Name,
                                     int *a_Number, int *a_MaxLength)
{
  hid_t   DsetId, DataspaceId, DataTypeId;
  herr_t  status;
  htri_t  IsSimple;
  int     Rank;
  hsize_t *Dims = NULL, *MaxDims = NULL;
  H5T_class_t DatasetClass;

  DsetId = H5Dopen1(a_Id, a_Name);
  if (DsetId < 0) {
   return DsetId;
  }

  /* open the datatype and make sure it is a string type */
  DataTypeId = H5Dget_type(DsetId);
  if (DataTypeId < 0) {
    return DataTypeId;
  }

  /* Make sure it is a string type */
  DatasetClass = H5Tget_class(DataTypeId);
  if (DatasetClass != H5T_STRING) {
    return ERROR_INCORRECT_DATATYPE;
  }

  /* Get the size of the string */
  *a_MaxLength = H5Tget_size(DataTypeId);

  H5Tclose(DataTypeId);

  /* Get the dataspace */
  DataspaceId = H5Dget_space(DsetId);

  /* Dataspace must be simply defined */
  IsSimple = H5Sis_simple(DataspaceId);
  if (!IsSimple) {
    H5Sclose(DataspaceId);
    H5Dclose(DsetId);
    return ERROR_ATTRIBUTE_NOT_SUPPORTED;
  }

  /* Get the rank and dimensions for simple dataspace */
  status = xfpGetSimpleDataspaceInfo(DataspaceId, &Rank, &Dims, &MaxDims);
  if (status < 0 || Rank != 1) {
    H5Sclose(DataspaceId);
    H5Dclose(DsetId);
    return ERROR_ATTRIBUTE_NOT_SUPPORTED;
  }

  /* See how many strings we have */
  *a_Number = (int)Dims[0];
 
  free(Dims);
  free(MaxDims);
  H5Sclose(DataspaceId);
  H5Dclose(DsetId);

  return XTRUE;
} /* xfGetPropertyStringLength */
/*-----------------------------------------------------------------------------*/
/* FUNCTION  xfGetPropertyType*/
/*! PURPOSE:   Gets the property type from a datase
 * - NOTES:*/
/*-----------------------------------------------------------------------------*/
XMDF_API xid xfGetPropertyType(xid a_GroupId, const char *a_Name, int *a_Type)
{
  hid_t        DsetId, DtypeId;
  H5T_class_t  Class;
  int          Sign;
  size_t       Size;

    /* open the dataset*/
  DsetId = H5Dopen1(a_GroupId, a_Name);
  if (DsetId < 0) {
   return DsetId;
  }

    /* get the datatype Id*/
  DtypeId = H5Dget_type(DsetId);
  if (DtypeId < 0) {
    H5Dclose(DsetId);
    return DtypeId;
  }

    /* get the class type*/
  Class = H5Tget_class(DtypeId);

    /* set a_type to appropriate type*/
  switch (Class) {
    case H5T_INTEGER:
      Sign = H5Tget_sign(DtypeId);
      if (Sign == H5T_SGN_NONE) {
        *a_Type = XF_TYPE_UINT;
      }
      else if (Sign == H5T_SGN_2) {
        *a_Type = XF_TYPE_INT;
      }
      else {
        *a_Type = XF_TYPE_OTHER;
      }
      break;
    case H5T_STRING:
      *a_Type = XF_TYPE_STRING;
      break;
    case H5T_FLOAT:
        /* find out if the type is a float or a double*/
      Size = H5Tget_size(DtypeId);
        if (Size == 4) {
          *a_Type = XF_TYPE_FLOAT;
        }
        else if (Size == 8) {
          *a_Type = XF_TYPE_DOUBLE;
        }
        else {
          *a_Type = XF_TYPE_OTHER;
        }
      break;
    default:
        /* if the type isn't an int, string, or float, set as other*/
      *a_Type = XF_TYPE_OTHER;
      break;
  }

  H5Tclose(DtypeId);
  H5Dclose(DsetId);

  return XTRUE;
} /* xfGetPropertyType*/

/******************************************************************************/
/*   FUNCTION  xfAllocateReadPropertyString */
/*! PURPOSE:   Reads the string dataset with the property ID
 *   - NOTES:     the variable does not need to be allocated before calling
 *             this function.  After this function is called, the variable needs
 *             to be free'd. */
/******************************************************************************/
XMDF_API xid xfAllocateReadPropertyString (xid a_Id, const char *a_Name,
                                           int *a_Number, int *a_MaxLength,
                                           char **a_Properties)
{
  int status;

  if (status == xfGetPropertyStringLength(a_Id, a_Name, a_Number,
      a_MaxLength)) {
    if (*a_Number > 0 && *a_MaxLength > 0) {
      *a_Properties = (char *)malloc(*a_Number*(*a_MaxLength+1)*
        sizeof(char));
      if (*a_Properties) {
        return(xfReadPropertyString(a_Id, a_Name, *a_Number, *a_MaxLength,
          *a_Properties));
      }
      else {
        status = ERROR_OTHER;
      }
    }
  }
  return status;
} /* xfAllocateReadPropertyString */

/******************************************************************************/
/*   FUNCTION  xfReadPropertyString */
/*! PURPOSE:   Reads the string dataset from the attributes directory
 *   - NOTES:     the variable attributes must already be allocated
 *             to a_Number*a_MaxLength*/
/******************************************************************************/
XMDF_API xid xfReadPropertyString(xid a_Id, const char *a_Name, 
                  int a_Number, int a_MaxLength, char *a_Properties)
{
  hid_t   DataTypeId, StringTypeId;
  herr_t  status = 1;
  int  StrSize;
  H5T_class_t DatasetClass;
  /* try putting all in same function */
  hid_t   DsetId, DataspaceId;
  htri_t  IsSimple;
  int     Rank;
  hsize_t *Dims = NULL, *MaxDims = NULL;

  /* make our own string type */
  /* Create the string type */
  StringTypeId = H5Tcopy(H5T_C_S1);
  if (StringTypeId < 0) {
    return StringTypeId;
  }
  H5Tset_strpad(StringTypeId, H5T_STR_NULLTERM);


#ifdef DATASETSTRING_WORKING
  H5Tclose(DataTypeId);
  return xfpReadDset1D(a_Id, a_Name, a_Number, H5T_STRING, StringTypeId,
                          (void*)a_Properties);
#else

  DsetId = H5Dopen1(a_Id, a_Name);
  if (DsetId < 0) {
   return DsetId;
  }
  /* open the datatype and make sure it is the correct type */
  DataTypeId = H5Dget_type(DsetId);
  if (DataTypeId < 0) {
    return DataTypeId;
  }
  /* Make sure it is a string type */
  DatasetClass = H5Tget_class(DataTypeId);
  if (DatasetClass != H5T_STRING) {
    H5Tclose(DataTypeId);
    return ERROR_INCORRECT_DATATYPE;
  }

  /* Make sure the maximum length is correct */
  StrSize = (int)H5Tget_size(DataTypeId);
  if (StrSize > a_MaxLength) {
    H5Tclose(DataTypeId);
    return ERROR_DATASET_INVALID;
  }
  H5Tclose(DataTypeId);

  /* Get the dataspace */
  DataspaceId = H5Dget_space(DsetId);

  /* Dataspace must be simply defined */
  IsSimple = H5Sis_simple(DataspaceId);
  if (!IsSimple) {
    H5Sclose(DataspaceId);
    H5Dclose(DsetId);
    H5Tclose(DataTypeId);
    return ERROR_DATASET_INVALID;
  }

  /* Get the rank and dimensions for simple dataspace */
  status = xfpGetSimpleDataspaceInfo(DataspaceId, &Rank, &Dims, &MaxDims);
  if (status < 0 || Rank != 1) {
    H5Sclose(DataspaceId);
    H5Aclose(DsetId);
    H5Tclose(DataTypeId);
    return ERROR_DATASET_INVALID;
  }

  /* Make sure the number of strings matches */
  if (a_Number != Dims[0]) {
    H5Sclose(DataspaceId);
    H5Aclose(DsetId);
    H5Tclose(DataTypeId);
    return ERROR_DATASET_INVALID;
  }
  /* set the size on the string type */
  H5Tset_size(StringTypeId, StrSize); 

  status = H5Dread(DsetId, StringTypeId, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                   a_Properties);

  H5Tclose(DataTypeId);
  H5Sclose(DataspaceId);
  H5Dclose(DsetId);
 
  free(Dims);
  free(MaxDims);
#endif

  H5Tclose(StringTypeId);

  return status;
} /* xfReadPropertyString */
/******************************************************************************/
/*   FUNCTION  xfReadPropertyUnsignedInt */
/*! PURPOSE:   Reads the unsigned int dataset from the properties directory
 *   - NOTES:     the variable properties must already be allocated to a_Number.*/
/******************************************************************************/
XMDF_API xid xfReadPropertyUnsignedInt(xid a_Id, const char *a_Name, 
                                       int a_Number, unsigned int *a_Properties)
{
  int status;

  status = xfpReadDset1D(a_Id, a_Name, a_Number, H5T_INTEGER,
              H5T_NATIVE_UINT, (void*)a_Properties);

  return status;
} /* xfReadPropertyUnsignedInt */
/******************************************************************************/
/*   FUNCTION  xfReadPropertyInt */
/*! PURPOSE:   Reads the int dataset from the properties directory
 *   - NOTES:     the variable properties must already be allocated to a_Number.*/
/******************************************************************************/
XMDF_API xid xfReadPropertyInt(xid a_Id, const char *a_Name, 
                                       int a_Number, int *a_Properties)
{
  int status;

  status = xfpReadDset1D(a_Id, a_Name, a_Number, H5T_INTEGER, H5T_NATIVE_INT,
              (void*)a_Properties);
  return status;
} /* xfReadPropertyInt */
/******************************************************************************/
/*   FUNCTION  xfReadPropertyDouble */
/*! PURPOSE:   Reads the double dataset from the properties directory
 *   - NOTES:     the variable properties must already be allocated to a_Number.*/
/******************************************************************************/
XMDF_API xid xfReadPropertyDouble(xid a_Id, const char *a_Name, 
                                       int a_Number, double *a_Properties)
{
  int status;

  status = xfpReadDset1D(a_Id, a_Name, a_Number, H5T_FLOAT, 
              H5T_NATIVE_DOUBLE, (void*)a_Properties);
  return status;
} /* xfReadPropertyDouble */
/******************************************************************************/
/*   FUNCTION  xfReadPropertyFloat */
/*! PURPOSE:   Reads the float dataset from the properties directory
 *   - NOTES:     the variable properties must already be allocated to a_Number. */
/******************************************************************************/
XMDF_API xid xfReadPropertyFloat(xid a_Id, const char *a_Name, 
                                       int a_Number, float *a_Properties)
{
  int status;

  status = xfpReadDset1D(a_Id, a_Name, a_Number, H5T_FLOAT, 
                         H5T_NATIVE_FLOAT, (void*)a_Properties);
  return status;
} /* xfReadPropertyFloat */
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfCreateGenericGroup */
/*! PURPOSE:   Create a group of generic type (anything may be stored below)
 * - NOTES:     These groups are just created for organization.  */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfCreateGenericGroup (xid a_FileId, const char *a_Path,
                                     xid *a_GroupId)
{
  *a_GroupId = xfpCreateGroup(a_FileId, a_Path, GROUP_TYPE_GENERIC); 

  return *a_GroupId;
} /* xfCreateGenericGroup */
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfCreateGroupForMatSim */
/*! PURPOSE:   Create a group of Material Simulation
 * - NOTES:     These groups are just created for organization.  */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfCreateGroupForMatSim (xid a_FileId, const char *a_Path,
                                     xid *a_GroupId)
{
  *a_GroupId = xfpCreateGroup(a_FileId, a_Path, GROUP_TYPE_MAT_SIM); 

  return *a_GroupId;
} /* xfCreateGenericGroup */
/******************************************************************************/
/*   FUNCTION  xfCreateGroupForMesh */
/*! PURPOSE:   Create an Xformat group to store a mesh
 *   - NOTES:    */
/******************************************************************************/
XMDF_API xid xfCreateGroupForMesh(xid a_FileId, const char *a_Path,
                                     xid *a_GroupId)
{
  xid  NodeGroup, ElemGroup;

  *a_GroupId = xfpCreateGroup(a_FileId, a_Path, GROUP_TYPE_MESH); 

  NodeGroup = xfpCreateGroup(*a_GroupId, MESH_NODE_GROUP,GROUP_TYPE_MESH_NODES);
  ElemGroup = xfpCreateGroup(*a_GroupId, MESH_ELEM_GROUP,GROUP_TYPE_MESH_ELEMS);

  /* see if there was an error */
  if (NodeGroup < 0 || ElemGroup < 0) {
    if (*a_GroupId > 0) {
      xfCloseGroup(*a_GroupId);
    }
    if (NodeGroup > 0) {
      xfCloseGroup(NodeGroup);
    }
    if (ElemGroup > 0) {
      xfCloseGroup(ElemGroup);
    }
    return -1;
  }
  xfCloseGroup(NodeGroup);
  xfCloseGroup(ElemGroup);

  return *a_GroupId;
} /* xfCreateGroupForMesh */
/******************************************************************************/
/*   FUNCTION  xfCreateGroupForGrid */
/*! PURPOSE:   Create an Xformat group to store a grid
 *   - NOTES:    */
/******************************************************************************/
XMDF_API xid xfCreateGroupForGrid(xid a_FileId, const char *a_Path, xid *a_GroupId)
{
  *a_GroupId = xfpCreateGroup(a_FileId, a_Path, GROUP_TYPE_GRID); 

  return *a_GroupId;
} /* xfCreateGroupForGrid */
/******************************************************************************/
/*   FUNCTION  xfCreateStationForGrid */
/*! PURPOSE:   Create an Xformat group to store a station for a grid
 *   - NOTES:    */
/******************************************************************************/
XMDF_API xid xfCreateStationForGrid(xid a_FileId, const char *a_Path, xid *a_GroupId)
{
  *a_GroupId = xfpCreateGroup(a_FileId, a_Path, GROUP_TYPE_STATION); 

  return *a_GroupId;
} /* xfCreateStationForGrid */

/******************************************************************************/
/*   FUNCTION  xfCreateGroupForXsec */
/*! PURPOSE:   Create an Xformat group to store a xsec
 *   - NOTES:    */
/******************************************************************************/
XMDF_API xid xfCreateGroupForXsec(xid a_FileId, const char *a_Path,
                                     xid *a_GroupId)
{
  xid  CrossSectionGeometryGroup;
  xid  LinePropertiesGroup;
  xid  PointPropertiesGroup;
  xid  ProfilesGroup;
  xid  XSECSGroup;

  *a_GroupId = xfpCreateGroup(a_FileId, a_Path, GROUP_TYPE_MESH); 

  CrossSectionGeometryGroup = xfpCreateGroup(*a_GroupId, XSEC_GEOMETRY_GROUP,
                                             GROUP_TYPE_XSEC_GEOM);
  LinePropertiesGroup = xfpCreateGroup(*a_GroupId, XSEC_LINE_GROUP,
                                       GROUP_TYPE_XSEC_LINE);
  PointPropertiesGroup = xfpCreateGroup(*a_GroupId, XSEC_POINT_GROUP,
                                        GROUP_TYPE_XSEC_POINT);
  ProfilesGroup = xfpCreateGroup(*a_GroupId, XSEC_PROFILES_GROUP,
                                 GROUP_TYPE_XSEC_PROFILES);
  XSECSGroup = xfpCreateGroup(*a_GroupId, XSEC_XSECS_GROUP,
                              GROUP_TYPE_XSEC_XSECS);

  /* see if there was an error */
  if (CrossSectionGeometryGroup < 0 || LinePropertiesGroup < 0 ||
      PointPropertiesGroup < 0 || ProfilesGroup < 0 ||
      XSECSGroup < 0) {
    if (*a_GroupId > 0) {
      xfCloseGroup(*a_GroupId);
    }
    if (CrossSectionGeometryGroup > 0) {
      xfCloseGroup(CrossSectionGeometryGroup);
    }
    if (LinePropertiesGroup > 0) {
      xfCloseGroup(LinePropertiesGroup);
    }
    if (PointPropertiesGroup > 0) {
      xfCloseGroup(PointPropertiesGroup);
    }
    if (ProfilesGroup > 0) {
      xfCloseGroup(ProfilesGroup);
    }
    if (XSECSGroup > 0) {
      xfCloseGroup(XSECSGroup);
    }
    return -1;
  }
  xfCloseGroup(CrossSectionGeometryGroup);
  xfCloseGroup(LinePropertiesGroup);
  xfCloseGroup(PointPropertiesGroup);
  xfCloseGroup(ProfilesGroup);
  xfCloseGroup(XSECSGroup);

  return *a_GroupId;
} /* xfCreateGroupForXsec */

/******************************************************************************/
/*   FUNCTION  xfOpenGroup */
/*! PURPOSE:   Open an Xformat group
 *   - NOTES:    */
/******************************************************************************/
XMDF_API xid xfOpenGroup (xid a_ParentId, const char *a_Path, xid *a_GroupId)
{
  *a_GroupId = H5Gopen1(a_ParentId, a_Path);
  return *a_GroupId;
} /* xfOpenGroup */
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfCloseGroup */
/*! PURPOSE:   Close an Xformat group
 * - NOTES:     */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfCloseGroup(xid a_GroupId)
{
  return H5Gclose(a_GroupId);
} /* xfCloseGroup */
/******************************************************************************/
/*   FUNCTION  xfGetGroupPathsSizeForMeshes */
/*! PURPOSE:   Find the number and maximum size for paths to datasets in an
 *             Xmdf file
 *   - NOTES:    */
/******************************************************************************/
XMDF_API xid xfGetGroupPathsSizeForMeshes(xid a_FileId, int *a_Number,
                                               int *a_Maxsize)
{
  return xfpNumGroupsOfType(a_FileId, GROUP_TYPE_MESH, a_Number, a_Maxsize, 
                     XFALSE);
} /* xfGetGroupPathsSizeForMeshes */
/******************************************************************************/
/*   FUNCTION  xfGetGroupPathsSizeForMeshes */
/*! PURPOSE:   Find the number and maximum size for paths to datasets in an
 *             Xmdf file
 *   - NOTES:    */
/******************************************************************************/
XMDF_API xid xfGetGroupPathsForMeshes (xid a_FileId, int a_Num,
                                              int a_Maxsize, char *a_Paths)
{
  return xfpPathGroupsOfType(a_FileId, GROUP_TYPE_MESH, 
                             a_Num, a_Maxsize, a_Paths, XFALSE);
} /* xfGetGroupPathsSizeForMeshes */
/******************************************************************************/
/*   FUNCTION  xfGetGroupPathsSizeForGrids */
/*! PURPOSE:   Find the number and maximum size for paths to Grids in an
 *             Xmdf file
 *   - NOTES:    */
/******************************************************************************/
XMDF_API xid xfGetGroupPathsSizeForGrids(xid a_FileId, int *a_Number,
                                               int *a_Maxsize)
{
  return xfpNumGroupsOfType(a_FileId, GROUP_TYPE_GRID, a_Number, a_Maxsize, 
                     XFALSE);
} /* xfGetGroupPathsSizeForGrids */
/******************************************************************************/
/*   FUNCTION  xfGetGroupPathsSizeForGrids */
/*! PURPOSE:   Find the number and maximum size for paths to datasets in an
 *             Xmdf file
 *   - NOTES:    */
/******************************************************************************/
XMDF_API xid xfGetGroupPathsForGrids (xid a_FileId, int a_Num,
                                              int a_Maxsize, char *a_Paths)
{
  return xfpPathGroupsOfType(a_FileId, GROUP_TYPE_GRID, a_Num, a_Maxsize,
                             a_Paths, XFALSE);
} /* xfGetGroupPathsSizeForGrids */
/******************************************************************************/
/*   FUNCTION  xfGetGroupPathsSizeForXsecs */
/*! PURPOSE:   Find the number and maximum size for paths to Xsecs in an
 *             Xmdf file
 *   - NOTES:    */
/******************************************************************************/
XMDF_API xid xfGetGroupPathsSizeForXsecs (xid a_FileId, int *a_Num,
                                          int *Maxsize)
{
  return xfpNumGroupsOfType(a_FileId, GROUP_TYPE_XSECS, a_Num, Maxsize, 
                            XFALSE);
} /* xfGetGroupPathsSizeForXsecs */
/******************************************************************************/
/*   FUNCTION  xfGetGroupPathsForXsecs */
/*! PURPOSE:   Find the number and maximum size for paths to Xsecs in an
 *             Xmdf file
 *   - NOTES:    */
/******************************************************************************/
XMDF_API xid xfGetGroupPathsForXsecs (xid a_FileId, int a_Num, int a_Maxsize,
                                      char *a_Paths)
{
  return xfpPathGroupsOfType(a_FileId, GROUP_TYPE_XSECS, a_Num, a_Maxsize,
                             a_Paths, XFALSE);
} /* xfGetGroupPathsForXsecs */
/******************************************************************************/
/*   FUNCTION  xfGetGroupPathsSizeForGeomPaths */
/*! PURPOSE:   Find the number and maximum size for paths to datasets in an
 *             Xmdf file
 *   - NOTES:    */
/******************************************************************************/
XMDF_API xid xfGetGroupPathsSizeForGeomPaths(xid a_FileId, int *a_Number,
                                               int *a_Maxsize)
{
  return xfpNumGroupsOfType(a_FileId, GROUP_TYPE_GEOMETRIC_PATH, a_Number, a_Maxsize, 
                     XFALSE);
} /* xfGetGroupPathsSizeForGeomPaths */
/******************************************************************************/
/*   FUNCTION  xfGetGroupPathsForGeomPaths */
/*! PURPOSE:   Find the number and maximum size for paths to datasets in an
 *             Xmdf file
 *   - NOTES:    */
/******************************************************************************/
XMDF_API xid xfGetGroupPathsForGeomPaths (xid a_FileId, int a_Num,
                                              int a_Maxsize, char *a_Paths)
{
  return xfpPathGroupsOfType(a_FileId, GROUP_TYPE_GEOMETRIC_PATH, 
                             a_Num, a_Maxsize, a_Paths, XFALSE);
} /* xfGetGroupPathsForGeomPaths */
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfOpenPropertyGroup */
/*! PURPOSE:   Open a property Group
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfOpenPropertyGroup (xid a_ParentId, xid *a_GroupId)
{
  *a_GroupId = H5Gopen1(a_ParentId, GROUP_TYPE_PROPERTIES);
  return *a_GroupId;
} /* xfOpenPropertyGroup*/
/* ---------------------------------------------------------------------------*/
/* FUNCTION  xfGetGroupAbsolutePathSize*/
/*! PURPOSE:   returns the size of the absolute path of the group with a given
 *           ID.
 * - NOTES:*/
/* ---------------------------------------------------------------------------*/
XMDF_API xid xfGetGroupAbsolutePathSize (xid a_GroupId, int *a_PathLength)
{
  size_t      size = 0;

  *a_PathLength = H5Iget_name(a_GroupId, NULL, size);

  return *a_PathLength;
} /* xfGetGroupAbsolutePathSize*/
/* ---------------------------------------------------------------------------*/
/* FUNCTION  xfGetGroupAbsolutePath*/
/*! PURPOSE:   returns the absolute path of the group with a given ID
 * - NOTES:*/
/* ---------------------------------------------------------------------------*/
XMDF_API xid xfGetGroupAbsolutePath (xid a_GroupId, int a_PathLength,
                                       char *a_Path)
{
  int         status = -1, i = 0;
  char        *temppath = NULL, *mypath = NULL;
  int         tmppathlen = 0;

  if (a_PathLength > 1) {
    tmppathlen = a_PathLength;
    tmppathlen++;
    temppath = (char *)malloc(tmppathlen*sizeof(char));
    mypath = (char *)malloc((tmppathlen-1)*sizeof(char));
    status = H5Iget_name(a_GroupId, temppath, tmppathlen);
    for (i=0; i<(tmppathlen-1); i++) {
      a_Path[i]=temppath[i+1];
    }
    status = 0;
  }
  if (temppath) {
    free(temppath);
  }
  if (mypath) {
    free(mypath);
  }
  return status;

} /* xfGetGroupAbsolutePath*/
/*-----------------------------------------------------------------------------*/
/* FUNCTION xfGetDoubleType*/
/*! PURPOSE:  Gets the current double type (either big endian
 *          or little endian)

 * - NOTES:*/
/*-----------------------------------------------------------------------------*/
 XMDF_API hid_t xfGetDoubleType ()
{
  return xfpGetDefaultDoubleType();
} /* xfGetDoubleType*/
/*-----------------------------------------------------------------------------*/
/* FUNCTION xfGetFloatType*/
/*! PURPOSE:  Gets the current float type (either big endian
 *          or little endian)

 * - NOTES:*/
/*-----------------------------------------------------------------------------*/
 XMDF_API hid_t xfGetFloatType ()
{
  return xfpGetDefaultFloatType();
} /* xfGetFloatType*/
/*-----------------------------------------------------------------------------*/
/* FUNCTION  xfSetFloatType*/
/*! PURPOSE:   Set the Float AND Double Types to be written to an XMDF file
 * - NOTES:     Pass in "1" for big endian and "0" for little endian.  */
/*-----------------------------------------------------------------------------*/
 XMDF_API void xfSetFloatType (int a_BigEndian)
{
  xfpDefaultWriteFloatType(a_BigEndian);
} /* xfSetFloatType*/
/* ---------------------------------------------------------------------------*/
/* FUNCTION  xfSetNumberOfElements*/
/*! PURPOSE:   Set the number of elements for the mes
 * - NOTES:     This number is written as an attribute to a mesh. */
/* ---------------------------------------------------------------------------*/
XMDF_API xid xfSetNumberOfElements (xid a_Id, int a_nElems)
{
  xid ElemGroup;
  int status;

  /* open the element group */
  status = xfOpenGroup(a_Id, MESH_ELEM_GROUP, &ElemGroup);
  if (status < 0) {
    return status;
  }
  status = xfpWriteAttributeInt(ElemGroup, MESH_ATT_NUMELEMS, 1, &a_nElems);
  xfCloseGroup(ElemGroup);

  return status;
} /* xfSetNumberOfElements */
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfSetAllElemsSameType */
/*! PURPOSE:   Write the element type dataset as a single type which is applied
 *           to all elements

 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfSetAllElemsSameType (xid a_Id, int a_Type)
{
  xid             ElemGroup;
  int             status;
  XDatasetParams  Params;

  /* open the element group */
  status = xfOpenGroup(a_Id, MESH_ELEM_GROUP, &ElemGroup);
  if (status < 0) {
    return status;
  }

  /* Setup the dataset parameters */
  xfpDsetParamsInit(&Params, 1, XFALSE, NONE); 
  xfpDsetParamsSetSizes(&Params, 0, 1, 1, 1);

  status = xfpWriteDsetInt(ElemGroup, MESH_DSET_ELEM_TYPES, H5T_NATIVE_INT,
                            &Params, &a_Type);
  xfpDsetParamsDestroy(&Params);
  
  xfCloseGroup(ElemGroup);

  return status;
} /* xfSetAllElemsSameType */
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfWriteElemTypes */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfWriteElemTypes (xid a_Id, int a_nElems, const int *a_Type,
                                 int a_Compression)
{
  xid             ElemGroup;
  int             status, nElemsInFile;
  XDatasetParams  Params;
  xmbool           Chunk = (xmbool)(a_Compression >= 0); /* chunk if using compression */
  int             Chunksize = a_nElems;

  if (a_Compression > 0 && a_nElems < XF_DATASET_MIN_COMPRESS_SIZE) {
    Chunk = XFALSE;
    a_Compression = NONE;
  }

  /* If we are chunking set the chunksize  */
  if (a_Compression >= 0) {
    Chunksize = a_nElems / XF_PROP_CHUNK_NUMBER;
    Chunksize = Xmin(Chunksize, XF_PROP_CHUNK_MIN);
    Chunksize = Xmax(Chunksize, XF_PROP_CHUNK_MAX);
    Chunksize = Xmin(Chunksize, a_nElems);
  }

  /* open the element group */
  status = xfOpenGroup(a_Id, MESH_ELEM_GROUP, &ElemGroup);
  if (status < 0) {
    return status;
  }
  
  /* make sure the number of elements is correct */
  status = xfpReadAttributeInt(ElemGroup, MESH_ATT_NUMELEMS, 1, &nElemsInFile);
  if (status < 0 || nElemsInFile != a_nElems) {
    xfCloseGroup(ElemGroup);
    return ERROR_ELEMENT_NUM_INCORRECT;
  }

  /* Setup the dataset parameters */
  xfpDsetParamsInit(&Params, 1, Chunk, a_Compression); 
  xfpDsetParamsSetSizes(&Params, 0, a_nElems, a_nElems, Chunksize);

  status = xfpWriteDsetInt(ElemGroup, MESH_DSET_ELEM_TYPES, H5T_NATIVE_INT,
                            &Params, a_Type);
  xfpDsetParamsDestroy(&Params);
  
  xfCloseGroup(ElemGroup);

  return status;
} /* xfWriteElemTypes */
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfWriteElemNodeIds */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfWriteElemNodeIds (xid a_Id, int a_nElems, int a_nMaxNodes,
                                   int *a_Ids, int a_Compression)
{
  xid             ElemGroup;
  int             status, nElemsInFile;
  XDatasetParams  Params;
  xmbool           Chunk = (a_Compression >= 0); /* chunk if using compression */
  int             Chunksize = a_nElems;

  if (a_Compression > 0 && a_nElems < XF_DATASET_MIN_COMPRESS_SIZE) {
    Chunk = XFALSE;
    a_Compression = NONE;
  }

  /* If we are chunking set the chunksize  */
  if (a_Compression >= 0) {
    Chunksize = a_nElems / XF_PROP_CHUNK_NUMBER;
    Chunksize = Xmin(Chunksize, XF_PROP_CHUNK_MIN); 
    Chunksize = Xmax(Chunksize, XF_PROP_CHUNK_MAX);
    Chunksize = Xmin(Chunksize, a_nElems);
  }

  /* open the element group */
  status = xfOpenGroup(a_Id, MESH_ELEM_GROUP, &ElemGroup);
  if (status < 0) {
    return status;
  }
  
  /* make sure the number of elements is correct */
  status = xfpReadAttributeInt(ElemGroup, MESH_ATT_NUMELEMS, 1, &nElemsInFile);
  if (status < 0 || nElemsInFile != a_nElems) {
    xfCloseGroup(ElemGroup);
    return ERROR_ELEMENT_NUM_INCORRECT;
  }

  /* Setup the dataset parameters */
  xfpDsetParamsInit(&Params, 2, Chunk, a_Compression); 
  xfpDsetParamsSetSizes(&Params, 0, a_nElems, a_nElems, Chunksize);
  xfpDsetParamsSetSizes(&Params, 1, a_nMaxNodes, a_nMaxNodes, a_nMaxNodes);

/* MKI BUGFIX Aug 9 2004 - Issue #26 - check for existing elements*/
  /* write element node ids, if elements exist*/
  if (nElemsInFile) {
    status = xfpWriteDsetInt(ElemGroup, MESH_DSET_ELEM_NODE_IDS, H5T_NATIVE_INT,
                              &Params, a_Ids);
  }
  xfpDsetParamsDestroy(&Params);
  
  xfCloseGroup(ElemGroup);

  return status;
} /* xfWriteElemNodeIds */
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfSetNumberOfNodes */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfSetNumberOfNodes (xid a_Id, int a_nNodes)
{
  xid NodeGroup;
  int status;

  /* open the node group */
  status = xfOpenGroup(a_Id, MESH_NODE_GROUP, &NodeGroup);
  if (status < 0) {
    return status;
  }
  status = xfpWriteAttributeInt(NodeGroup, MESH_ATT_NUMNODES, 1, &a_nNodes);
  xfCloseGroup(NodeGroup);

  return status;
} /* xfSetNumberOfNodes */
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfWriteXNodeLocations */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfWriteXNodeLocations (xid a_Id, int a_nNodes, double *a_Locs,
                                      int a_Compression)
{
  xid             NodeGroup, LocId, DsetSpaceId, MemSpaceId, DsetCreateParamsId;
  int             status, nNodesInFile;
  XDatasetParams  Params;
  xmbool           Chunk = (a_Compression >= 0); /* chunk if using compression */
  hsize_t         Dims[2];
/*   hssize_t        NumSelected; */
  int             Chunksize = a_nNodes;

  if (a_Compression > 0 && a_nNodes < XF_DATASET_MIN_COMPRESS_SIZE) {
    Chunk = XFALSE;
    a_Compression = NONE;
  }

  /* If we are chunking set the chunksize  */
  if (a_Compression >= 0) {
    Chunksize = a_nNodes / XF_PROP_CHUNK_NUMBER;
    Chunksize = Xmin(Chunksize, XF_PROP_CHUNK_MIN); 
    Chunksize = Xmax(Chunksize, XF_PROP_CHUNK_MAX);
    Chunksize = Xmin(Chunksize, a_nNodes);
  }

  /* open the element group */
  status = xfOpenGroup(a_Id, MESH_NODE_GROUP, &NodeGroup);
  if (status < 0) {
    return status;
  }
  
  /* make sure the number of elements is correct */
  status = xfpReadAttributeInt(NodeGroup, MESH_ATT_NUMNODES, 1, &nNodesInFile);
  if (status < 0 || nNodesInFile != a_nNodes) {
    xfCloseGroup(NodeGroup);
    return ERROR_NODE_NUM_INCORRECT;
  }

  /* Setup the dataset parameters */
  xfpDsetParamsInit(&Params, 2, Chunk, a_Compression); 
  xfpDsetParamsSetSizes(&Params, 0, a_nNodes, a_nNodes, Chunksize);
  xfpDsetParamsSetSizes(&Params, 1, 3, 3, 3);

  /* Create the dataset dataspace */
  DsetSpaceId = xfpDsetParamsCreateSpace(&Params);  

  /* Create the dataset creation properties */
  DsetCreateParamsId = xfpDsetParamsCreateProps(&Params);

  if (DsetSpaceId < 0 || DsetCreateParamsId < 0) {
    if (DsetSpaceId >= 0) {
      H5Sclose(DsetSpaceId);
    }
    if (DsetCreateParamsId >= 0) {
      H5Pclose(DsetCreateParamsId);
    }
  }
  xfpDsetParamsDestroy(&Params);
  

  /* Create the dataset */
  LocId = H5Dcreate1(NodeGroup, MESH_DSET_NODE_LOCS, xfpGetDefaultDoubleType(),
                    DsetSpaceId, DsetCreateParamsId);
/*  if (LocId < 0) { */
/*    H5Epush(__FILE__, __FUNCTION__, __LINE__,  */
/*  } */

  /* Close the dataspace and dataset params */
  DsetSpaceId = H5Sclose(DsetSpaceId);
  DsetCreateParamsId = H5Pclose(DsetCreateParamsId);
  if (DsetSpaceId < 0 || DsetCreateParamsId < 0) {
    xfCloseGroup(NodeGroup);
    return -1;
  }

  /* Select the dataspace for the X values of the dataset */
  DsetSpaceId = H5Dget_space(LocId);


  {
    const hsize_t        start[2] = {0,0};
    hsize_t         count[2], block[2];

/*     start[0] = 0; */
/*     start[1] = 0; */
    count[0] = 1;
    count[1] = 1;
    block[0] = a_nNodes;
    block[1] = 1;
    status = H5Sselect_hyperslab(DsetSpaceId, H5S_SELECT_SET, start, NULL,
                                 count, block);
  }
  if (status < 0) {
    H5Dclose(LocId);
    xfCloseGroup(NodeGroup);
    return status;
  }
  
  status = H5Sselect_valid(DsetSpaceId);
  if (status < 0) {
    H5Dclose(LocId);
    xfCloseGroup(NodeGroup);
    return status;
  }

/*   NumSelected = H5Sget_select_npoints(DsetSpaceId); */

  /* Create the memory dataspace */
  Dims[0] = a_nNodes;
  MemSpaceId = H5Screate_simple(1, Dims, NULL);

  /* Write the X locations */
  status = H5Dwrite(LocId, H5T_NATIVE_DOUBLE, MemSpaceId, DsetSpaceId,
                    H5P_DEFAULT, a_Locs);

  H5Sclose(MemSpaceId);
  H5Sclose(DsetSpaceId);
  H5Dclose(LocId);
  xfCloseGroup(NodeGroup);

  return status;
} /* xfWriteXNodeLocations */
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfWriteYNodeLocations */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfWriteYNodeLocations (xid a_Id, int a_nNodes, double *a_Locs)
{
  xid             NodeGroup, LocId, DsetSpaceId, MemSpaceId;
  int             status, nNodesInFile;
  hsize_t         Dims[2];
/*   hssize_t NumSelected; */

  /* open the element group */
  status = xfOpenGroup(a_Id, MESH_NODE_GROUP, &NodeGroup);
  if (status < 0) {
    return status;
  }
  
  /* make sure the number of elements is correct */
  status = xfpReadAttributeInt(NodeGroup, MESH_ATT_NUMNODES, 1, &nNodesInFile);
  if (status < 0 || nNodesInFile != a_nNodes) {
    xfCloseGroup(NodeGroup);
    return ERROR_NODE_NUM_INCORRECT;
  }

  /* Open the Dataset */
  LocId = H5Dopen1(NodeGroup, MESH_DSET_NODE_LOCS); 

  /* Select the dataspace for the Y values of the dataset */
  DsetSpaceId = H5Dget_space(LocId);

  {
    const hsize_t        start[2] = {0,1};
    hsize_t         count[2], block[2];

    count[0] = 1;
    count[1] = 1;
    block[0] = a_nNodes;
    block[1] = 1;
    status = H5Sselect_hyperslab(DsetSpaceId, H5S_SELECT_SET, start, NULL,
                                 count, block);
  }
  if (status < 0) {
    H5Dclose(LocId);
    xfCloseGroup(NodeGroup);
    return status;
  }
  
  status = H5Sselect_valid(DsetSpaceId);
  if (status < 0) {
    H5Dclose(LocId);
    xfCloseGroup(NodeGroup);
    return status;
  }

/*   NumSelected = H5Sget_select_npoints(DsetSpaceId); */

  /* Create the memory dataspace */
  Dims[0] = a_nNodes;
  MemSpaceId = H5Screate_simple(1, Dims, NULL);

  /* Write the X locations */
  status = H5Dwrite(LocId, H5T_NATIVE_DOUBLE, MemSpaceId, DsetSpaceId,
                    H5P_DEFAULT, a_Locs);

  H5Sclose(MemSpaceId);
  H5Sclose(DsetSpaceId);
  H5Dclose(LocId);
  xfCloseGroup(NodeGroup);

  return status;
} /* xfWriteYNodeLocations */
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfWriteZNodeLocations */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfWriteZNodeLocations (xid a_Id, int a_nNodes, double *a_Locs)
{
  xid             NodeGroup, LocId, DsetSpaceId, MemSpaceId;
  int             status, nNodesInFile;
  hsize_t         Dims[2];
  hsize_t         start[2];
  hsize_t         count[2], block[2];
/*   hssize_t NumSelected; */

  /* open the element group */
  status = xfOpenGroup(a_Id, MESH_NODE_GROUP, &NodeGroup);
  if (status < 0) {
    return status;
  }
  
  /* make sure the number of elements is correct */
  status = xfpReadAttributeInt(NodeGroup, MESH_ATT_NUMNODES, 1, &nNodesInFile);
  if (status < 0 || nNodesInFile != a_nNodes) {
    xfCloseGroup(NodeGroup);
    return ERROR_NODE_NUM_INCORRECT;
  }

  /* Open the Dataset */
  LocId = H5Dopen1(NodeGroup, MESH_DSET_NODE_LOCS); 

  /* Select the dataspace for the Y values of the dataset */
  DsetSpaceId = H5Dget_space(LocId);

  start[0] = 0;
  start[1] = 2;
  count[0] = 1;
  count[1] = 1;
  block[0] = a_nNodes;
  block[1] = 1;
  status = H5Sselect_hyperslab(DsetSpaceId, H5S_SELECT_SET, start, NULL,
                               count, block);
  if (status < 0) {
    H5Dclose(LocId);
    xfCloseGroup(NodeGroup);
    return status;
  }
  
  status = H5Sselect_valid(DsetSpaceId);
  if (status < 0) {
    H5Dclose(LocId);
    xfCloseGroup(NodeGroup);
    return status;
  }

/*   NumSelected = H5Sget_select_npoints(DsetSpaceId); */

  /* Create the memory dataspace */
  Dims[0] = a_nNodes;
  MemSpaceId = H5Screate_simple(1, Dims, NULL);

  /* Write the X locations */
  status = H5Dwrite(LocId, H5T_NATIVE_DOUBLE, MemSpaceId, DsetSpaceId,
                    H5P_DEFAULT, a_Locs);

  H5Sclose(MemSpaceId);
  H5Sclose(DsetSpaceId);
  H5Dclose(LocId);
  xfCloseGroup(NodeGroup);

  return status;
} /* xfWriteZNodeLocations */
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfGetNumberOfElements */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfGetNumberOfElements (xid a_Id, int *a_nElems)
{
  xid ElemGroup;
  int status;

  if (!xfpIsGroupOfType(a_Id, GROUP_TYPE_MESH)) {
    return ERROR_NOT_MESH_GROUP;
  }
  /* open the element group */
  status = xfOpenGroup(a_Id, MESH_ELEM_GROUP, &ElemGroup);
  if (status < 0) {
    return status;
  }
  status = xfpReadAttributeInt(ElemGroup, MESH_ATT_NUMELEMS, 1, a_nElems);
  xfCloseGroup(ElemGroup);

  return status;
} /* xfGetNumberOfElements */
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfAreAllElemsSameType */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfAreAllElemsSameType (xid a_Id, xmbool *a_Same)
{
  xid      ElemGroup, ElemTypeId, SpaceId;
  herr_t   Simple;
  int      Rank, status;
  hsize_t  *Dims = NULL, *Maxdims = NULL;

  if (!xfpIsGroupOfType(a_Id, GROUP_TYPE_MESH)) {
    return ERROR_NOT_MESH_GROUP;
  }

  /* open the element group */
  status = xfOpenGroup(a_Id, MESH_ELEM_GROUP, &ElemGroup);
  if (status < 0) {
    return status;
  }
  /* open the element type dataset */
  ElemTypeId = H5Dopen1(ElemGroup, MESH_DSET_ELEM_TYPES);
  if (ElemTypeId < 0) {
    return ERROR_MESH_INVALID;
  }

  /* Get the dataspace and make sure it is a simple dataset of rank 1 */
  SpaceId = H5Dget_space(ElemTypeId);
  if (SpaceId < 0) {
    H5Dclose(ElemTypeId);
    return ERROR_MESH_INVALID;
  }

  Simple = H5Sis_simple(SpaceId);
  if (!Simple) {
    H5Sclose(SpaceId);
    H5Dclose(ElemTypeId);
    return ERROR_MESH_INVALID;
  }

  xfpGetSimpleDataspaceInfo(SpaceId, &Rank, &Dims, &Maxdims);
  if (Rank != 1) {
    H5Sclose(SpaceId);
    H5Dclose(ElemTypeId);
    free(Dims);
    free(Maxdims);
    return ERROR_MESH_INVALID;
  }

  if (Dims[0] == 1) {
    *a_Same = XTRUE;
  }
  else {
    *a_Same = XFALSE;
  }

  free(Dims);
  free(Maxdims);

  H5Sclose(SpaceId);
  H5Dclose(ElemTypeId);
  return status;
} /* xfAreAllElemsSameType */
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfReadElemTypesSingleValue */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfReadElemTypesSingleValue (xid a_Id, int *a_Type)
{
  xid      ElemGroup, ElemTypeId, SpaceId;
  herr_t   Simple;
  int      Rank, status;
  hsize_t  *Dims = NULL, *Maxdims = NULL;

  if (!xfpIsGroupOfType(a_Id, GROUP_TYPE_MESH)) {
    return ERROR_NOT_MESH_GROUP;
  }
  /* open the element group */
  status = xfOpenGroup(a_Id, MESH_ELEM_GROUP, &ElemGroup);
  if (status < 0) {
    return status;
  }

  /* open the element type dataset */
  ElemTypeId = H5Dopen1(ElemGroup, MESH_DSET_ELEM_TYPES);
  if (ElemTypeId < 0) {
    return ERROR_MESH_INVALID;
  }

  /* Get the dataspace and make sure it is a simple dataset of rank 1 */
  SpaceId = H5Dget_space(ElemTypeId);
  if (SpaceId < 0) {
    H5Dclose(ElemTypeId);
    return ERROR_MESH_INVALID;
  }

  Simple = H5Sis_simple(SpaceId);
  if (!Simple) {
    H5Sclose(SpaceId);
    H5Dclose(ElemTypeId);
    return ERROR_MESH_INVALID;
  }

  xfpGetSimpleDataspaceInfo(SpaceId, &Rank, &Dims, &Maxdims);
  if (Rank != 1) {
    H5Sclose(SpaceId);
    H5Dclose(ElemTypeId);
    free(Dims);
    free(Maxdims);
    return ERROR_MESH_INVALID;
  }

  if (Dims[0] != 1) {
    H5Sclose(SpaceId);
    H5Dclose(ElemTypeId);
    free(Dims);
    free(Maxdims);
    return ERROR_MESH_INVALID;
  }

  status = H5Dread(ElemTypeId, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                   a_Type); 

  free(Dims);
  free(Maxdims);

  H5Sclose(SpaceId);
  H5Dclose(ElemTypeId);
  return status;
} /* xfReadElemTypesSingleValue */
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfReadElemTypes */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfReadElemTypes (xid a_Id, int a_nElems, int *a_Type)
{
  xid      ElemGroup, ElemTypeId, SpaceId;
  herr_t   Simple;
  int      Rank, i, status;
  int      SingleType;
  hsize_t  *Dims = NULL, *Maxdims = NULL;

  if (!xfpIsGroupOfType(a_Id, GROUP_TYPE_MESH)) {
    return ERROR_NOT_MESH_GROUP;
  }

  /* open the element group */
  status = xfOpenGroup(a_Id, MESH_ELEM_GROUP, &ElemGroup);
  if (status < 0) {
    return status;
  }
  /* open the element type dataset */
  ElemTypeId = H5Dopen1(ElemGroup, MESH_DSET_ELEM_TYPES);
  if (ElemTypeId < 0) {
    return ERROR_MESH_INVALID;
  }

  /* Get the dataspace and make sure it is a simple dataset of rank 1 */
  SpaceId = H5Dget_space(ElemTypeId);
  if (SpaceId < 0) {
    H5Dclose(ElemTypeId);
    return ERROR_MESH_INVALID;
  }

  Simple = H5Sis_simple(SpaceId);
  if (!Simple) {
    H5Sclose(SpaceId);
    H5Dclose(ElemTypeId);
    return ERROR_MESH_INVALID;
  }

  xfpGetSimpleDataspaceInfo(SpaceId, &Rank, &Dims, &Maxdims);
  if (Rank != 1) {
    H5Sclose(SpaceId);
    H5Dclose(ElemTypeId);
    free(Dims);
    free(Maxdims);
    return ERROR_MESH_INVALID;
  }

  if (Dims[0] == 1) {
    /* The file was saved only specifying a single type, set whole array to be */
    /* this type */
    status = H5Dread(ElemTypeId, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                     &SingleType);
    for (i = 0; i < a_nElems; i++) {
      a_Type[i] = SingleType;
    }
  }
  else {
    status = H5Dread(ElemTypeId, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                     a_Type);
  }

  free(Dims);
  free(Maxdims);

  H5Sclose(SpaceId);
  H5Dclose(ElemTypeId);
  return status;
}
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfGetMaxNodesInElem */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfGetMaxNodesInElem (xid a_Id, int *a_nMaxNodes)
{
  xid      ElemGroup, ElemNodeId, SpaceId;
  herr_t   Simple;
  int      Rank, status;
  hsize_t  *Dims = NULL, *Maxdims = NULL;

  if (!xfpIsGroupOfType(a_Id, GROUP_TYPE_MESH)) {
    return ERROR_NOT_MESH_GROUP;
  }
  /* open the element group */
  status = xfOpenGroup(a_Id, MESH_ELEM_GROUP, &ElemGroup);
  if (status < 0) {
    return status;
  }

  /* open the element node id dataset */
  ElemNodeId = H5Dopen1(ElemGroup, MESH_DSET_ELEM_NODE_IDS);
  if (ElemNodeId < 0) {
    return ERROR_MESH_INVALID;
  }

  /* Get the dataspace and make sure it is a simple dataset of rank 1 */
  SpaceId = H5Dget_space(ElemNodeId);
  if (SpaceId < 0) {
    H5Dclose(ElemNodeId);
    return ERROR_MESH_INVALID;
  }

  Simple = H5Sis_simple(SpaceId);
  if (!Simple) {
    H5Sclose(SpaceId);
    H5Dclose(ElemNodeId);
    return ERROR_MESH_INVALID;
  }

  xfpGetSimpleDataspaceInfo(SpaceId, &Rank, &Dims, &Maxdims);
  if (Rank != 2) {
    H5Sclose(SpaceId);
    H5Dclose(ElemNodeId);
    free(Dims);
    free(Maxdims);
    return ERROR_MESH_INVALID;
  }

  /* the second dimension is the maximum number of nodes in an element */
  *a_nMaxNodes = (int)Dims[1];

  free(Dims);
  free(Maxdims);

  H5Sclose(SpaceId);
  H5Dclose(ElemNodeId);
  return status;
} /* xfGetMaxNodesInElem */
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfReadElemNodeIds */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfReadElemNodeIds (xid a_Id, int a_nElems, int a_nMaxNodes,
                                  int *a_Ids)
{
  xid      ElemGroup, ElemNodeId, SpaceId;
  herr_t   Simple;
  int      Rank, status;
  hsize_t  *Dims = NULL, *Maxdims = NULL;

  if (!xfpIsGroupOfType(a_Id, GROUP_TYPE_MESH)) {
    return ERROR_NOT_MESH_GROUP;
  }

  /* open the element group */
  status = xfOpenGroup(a_Id, MESH_ELEM_GROUP, &ElemGroup);
  if (status < 0) {
    return status;
  }

/* MKI BUGFIX Aug 9 2004 - Issue #26 - check for existing elements*/
  if (a_nElems > 0) {
    /* open the element node id dataset*/
    ElemNodeId = H5Dopen1(ElemGroup, MESH_DSET_ELEM_NODE_IDS);
    if (ElemNodeId < 0) {
      return ERROR_MESH_INVALID;
    }

    /* Get the dataspace and make sure it is a simple dataset of rank 1 */
    SpaceId = H5Dget_space(ElemNodeId);
    if (SpaceId < 0) {
      H5Dclose(ElemNodeId);
      return ERROR_MESH_INVALID;
    }

    Simple = H5Sis_simple(SpaceId);
    if (!Simple) {
      H5Sclose(SpaceId);
      H5Dclose(ElemNodeId);
      return ERROR_MESH_INVALID;
    }

    xfpGetSimpleDataspaceInfo(SpaceId, &Rank, &Dims, &Maxdims);
    if (Rank != 2) {
      H5Sclose(SpaceId);
      H5Dclose(ElemNodeId);
      free(Dims);
      free(Maxdims);
      return ERROR_MESH_INVALID;
    }
  
    /* Make sure the dimensions are consistent*/
    if (Dims[0] != a_nElems || Dims[1] != a_nMaxNodes) {
      H5Sclose(SpaceId);
      H5Dclose(ElemNodeId);
      free(Dims);
      free(Maxdims);
      return ERROR_MESH_INVALID;
    }

    /* Read in the element node ids*/
    status = H5Dread(ElemNodeId, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                     a_Ids);

    free(Dims);
    free(Maxdims);

    H5Sclose(SpaceId);
    H5Dclose(ElemNodeId);
  }
  return status;
} /* xfReadElemNodeIds */
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfGetNumberOfNodes */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfGetNumberOfNodes (xid a_Id, int *a_nNodes)
{
  xid NodeGroup;
  int status;

  if (!xfpIsGroupOfType(a_Id, GROUP_TYPE_MESH)) {
    return ERROR_NOT_MESH_GROUP;
  }

  /* open the node group */
  status = xfOpenGroup(a_Id, MESH_NODE_GROUP, &NodeGroup);
  if (status < 0) {
    return status;
  }
  status = xfpReadAttributeInt(NodeGroup, MESH_ATT_NUMNODES, 1, a_nNodes);
  xfCloseGroup(NodeGroup);

  return status;
} /* xfGetNumberOfNodes */
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfReadXNodeLocations */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfReadXNodeLocations (xid a_Id, int a_nNodes, double *a_Locs)
{
  xid      NodeGroup, NodeLocId, SpaceId, MemSpaceId;
  herr_t   Simple;
  int      Rank, status;
  hsize_t  *Dims = NULL, *Maxdims = NULL;
  hsize_t  start[2];
  hsize_t  count[2], block[2];

  if (!xfpIsGroupOfType(a_Id, GROUP_TYPE_MESH)) {
    return ERROR_NOT_MESH_GROUP;
  }

  /* open the node group */
  status = xfOpenGroup(a_Id, MESH_NODE_GROUP, &NodeGroup);
  if (status < 0) {
    return status;
  }

  /* open the element node id dataset */
  NodeLocId = H5Dopen1(NodeGroup, MESH_DSET_NODE_LOCS);
  if (NodeLocId < 0) {
    return ERROR_MESH_INVALID;
  }

  /* Get the dataspace and make sure it is a simple dataset of rank 1 */
  SpaceId = H5Dget_space(NodeLocId);
  if (SpaceId < 0) {
    H5Dclose(NodeLocId);
    return ERROR_MESH_INVALID;
  }

  Simple = H5Sis_simple(SpaceId);
  if (!Simple) {
    H5Sclose(SpaceId);
    H5Dclose(NodeLocId);
    return ERROR_MESH_INVALID;
  }

  xfpGetSimpleDataspaceInfo(SpaceId, &Rank, &Dims, &Maxdims);
  if (Rank != 2) {
    H5Sclose(SpaceId);
    H5Dclose(NodeLocId);
    free(Dims);
    free(Maxdims);
    return ERROR_MESH_INVALID;
  }
  
  /* Make sure the dimensions are consistent */
  if (Dims[0] != a_nNodes || Dims[1] != 3) {
    H5Sclose(SpaceId);
    H5Dclose(NodeLocId);
    free(Dims);
    free(Maxdims);
    return ERROR_MESH_INVALID;
  }

  /* Select the X locations hyperslab */
  start[0] = 0;
  start[1] = 0;
  count[0] = 1;
  count[1] = 1;
  block[0] = a_nNodes;
  block[1] = 1;
  status = H5Sselect_hyperslab(SpaceId, H5S_SELECT_SET, start, NULL,
                               count, block);
  if (status < 0) {
    H5Dclose(NodeLocId);
    xfCloseGroup(NodeGroup);
    return status;
  }
  
  status = H5Sselect_valid(SpaceId);
  if (status < 0) {
    H5Dclose(NodeLocId);
    xfCloseGroup(NodeGroup);
    return status;
  }

  /* Create the memory dataspace */
  Dims[0] = a_nNodes;
  MemSpaceId = H5Screate_simple(1, Dims, NULL);
 
  status = H5Dread(NodeLocId, H5T_NATIVE_DOUBLE, MemSpaceId, SpaceId,
                   H5P_DEFAULT, a_Locs);

  free(Dims);
  free(Maxdims);

  H5Sclose(SpaceId);
  H5Dclose(NodeLocId);
  return status;
} /* xfReadXNodeLocations */
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfReadYNodeLocations */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfReadYNodeLocations (xid a_Id, int a_nNodes, double *a_Locs)
{
  xid      NodeGroup, NodeLocId, SpaceId, MemSpaceId;
  herr_t   Simple;
  int      Rank, status;
  hsize_t  *Dims = NULL, *Maxdims = NULL;
  hsize_t  start[2];
  hsize_t  count[2], block[2];

  if (!xfpIsGroupOfType(a_Id, GROUP_TYPE_MESH)) {
    return ERROR_NOT_MESH_GROUP;
  }

  /* open the node group */
  status = xfOpenGroup(a_Id, MESH_NODE_GROUP, &NodeGroup);
  if (status < 0) {
    return status;
  }

  /* open the element node id dataset */
  NodeLocId = H5Dopen1(NodeGroup, MESH_DSET_NODE_LOCS);
  if (NodeLocId < 0) {
    return ERROR_MESH_INVALID;
  }

  /* Get the dataspace and make sure it is a simple dataset of rank 1 */
  SpaceId = H5Dget_space(NodeLocId);
  if (SpaceId < 0) {
    H5Dclose(NodeLocId);
    return ERROR_MESH_INVALID;
  }

  Simple = H5Sis_simple(SpaceId);
  if (!Simple) {
    H5Sclose(SpaceId);
    H5Dclose(NodeLocId);
    return ERROR_MESH_INVALID;
  }

  xfpGetSimpleDataspaceInfo(SpaceId, &Rank, &Dims, &Maxdims);
  if (Rank != 2) {
    H5Sclose(SpaceId);
    H5Dclose(NodeLocId);
    free(Dims);
    free(Maxdims);
    return ERROR_MESH_INVALID;
  }
  
  /* Make sure the dimensions are consistent */
  if (Dims[0] != a_nNodes || Dims[1] != 3) {
    H5Sclose(SpaceId);
    H5Dclose(NodeLocId);
    free(Dims);
    free(Maxdims);
    return ERROR_MESH_INVALID;
  }

  /* Select the Y locations hyperslab */
  start[0] = 0;
  start[1] = 1;
  count[0] = 1;
  count[1] = 1;
  block[0] = a_nNodes;
  block[1] = 1;
  status = H5Sselect_hyperslab(SpaceId, H5S_SELECT_SET, start, NULL,
                               count, block);
  if (status < 0) {
    H5Dclose(NodeLocId);
    xfCloseGroup(NodeGroup);
    return status;
  }
  
  status = H5Sselect_valid(SpaceId);
  if (status < 0) {
    H5Dclose(NodeLocId);
    xfCloseGroup(NodeGroup);
    return status;
  }

  /* Create the memory dataspace */
  Dims[0] = a_nNodes;
  MemSpaceId = H5Screate_simple(1, Dims, NULL);
 
  status = H5Dread(NodeLocId, H5T_NATIVE_DOUBLE, MemSpaceId, SpaceId,
                   H5P_DEFAULT, a_Locs);

  free(Dims);
  free(Maxdims);

  H5Sclose(SpaceId);
  H5Dclose(NodeLocId);
  return status;
} /* xfReadYNodeLocations */
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfReadZNodeLocations */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfReadZNodeLocations (xid a_Id, int a_nNodes, double *a_Locs)
{
  xid      NodeGroup, NodeLocId, SpaceId, MemSpaceId;
  herr_t   Simple;
  int      Rank, status;
  hsize_t  *Dims = NULL, *Maxdims = NULL;
  hsize_t  start[2];
  hsize_t  count[2], block[2];

  if (!xfpIsGroupOfType(a_Id, GROUP_TYPE_MESH)) {
    return ERROR_NOT_MESH_GROUP;
  }

  /* open the node group */
  status = xfOpenGroup(a_Id, MESH_NODE_GROUP, &NodeGroup);
  if (status < 0) {
    return status;
  }

  /* open the element node id dataset */
  NodeLocId = H5Dopen1(NodeGroup, MESH_DSET_NODE_LOCS);
  if (NodeLocId < 0) {
    return ERROR_MESH_INVALID;
  }

  /* Get the dataspace and make sure it is a simple dataset of rank 1 */
  SpaceId = H5Dget_space(NodeLocId);
  if (SpaceId < 0) {
    H5Dclose(NodeLocId);
    return ERROR_MESH_INVALID;
  }

  Simple = H5Sis_simple(SpaceId);
  if (!Simple) {
    H5Sclose(SpaceId);
    H5Dclose(NodeLocId);
    return ERROR_MESH_INVALID;
  }

  xfpGetSimpleDataspaceInfo(SpaceId, &Rank, &Dims, &Maxdims);
  if (Rank != 2) {
    H5Sclose(SpaceId);
    H5Dclose(NodeLocId);
    free(Dims);
    free(Maxdims);
    return ERROR_MESH_INVALID;
  }
  
  /* Make sure the dimensions are consistent */
  if (Dims[0] != a_nNodes || Dims[1] != 3) {
    H5Sclose(SpaceId);
    H5Dclose(NodeLocId);
    free(Dims);
    free(Maxdims);
    return ERROR_MESH_INVALID;
  }

  /* Select the Z locations hyperslab */
  start[0] = 0;
  start[1] = 2;
  count[0] = 1;
  count[1] = 1;
  block[0] = a_nNodes;
  block[1] = 1;
  status = H5Sselect_hyperslab(SpaceId, H5S_SELECT_SET, start, NULL,
                               count, block);
  if (status < 0) {
    H5Dclose(NodeLocId);
    xfCloseGroup(NodeGroup);
    return status;
  }
  
  status = H5Sselect_valid(SpaceId);
  if (status < 0) {
    H5Dclose(NodeLocId);
    xfCloseGroup(NodeGroup);
    return status;
  }

  /* Create the memory dataspace */
  Dims[0] = a_nNodes;
  MemSpaceId = H5Screate_simple(1, Dims, NULL);
 
  status = H5Dread(NodeLocId, H5T_NATIVE_DOUBLE, MemSpaceId, SpaceId,
                   H5P_DEFAULT, a_Locs);

  free(Dims);
  free(Maxdims);

  H5Sclose(SpaceId);
  H5Dclose(NodeLocId);

  return status;
} /* xfReadZNodeLocations */
/* --------------------------------------------------------------------------*/
/* FUNCTION  xfCreateMeshPropertyGroup                                       */
/*! PURPOSE:
 * - NOTES:                                                                     */
/* --------------------------------------------------------------------------*/
XMDF_API xid xfCreateMeshPropertyGroup (xid a_Id, xid *a_PropId)
{
  int     error = 1;
  herr_t (*old_func)(void*);
  void    *old_client_data;

  H5Eget_auto1(&old_func, &old_client_data);

  /* Turn off error handling */
  H5Eset_auto1(NULL, NULL);

  /* Open group */
  error = xfOpenGroup(a_Id, MESH_PROP_GROUP, a_PropId);

  /* Restore previous error handler */
  H5Eset_auto1(old_func, old_client_data);

  /* If group does not exist create it */
  if (error < 0) {
    *a_PropId = xfpCreateGroup(a_Id, MESH_PROP_GROUP, GROUP_TYPE_PROPERTIES);
    error = *a_PropId;
  }

  return error;
} /* xfCreateMeshPropertyGroup */
/* --------------------------------------------------------------------------*/
/* FUNCTION  xfGetMeshPropertyGroup                                          */
/*! PURPOSE:
 * - NOTES:                                                                     */
/* --------------------------------------------------------------------------*/
XMDF_API xid xfGetMeshPropertyGroup (xid a_Id, xid *a_PropId)
{
  int     error = 1;
  herr_t (*old_func)(void*);
  void    *old_client_data;

  H5Eget_auto1(&old_func, &old_client_data);

  /* Turn off error handling */
  H5Eset_auto1(NULL, NULL);

  /* Open group */
  error = xfOpenGroup(a_Id, MESH_PROP_GROUP, a_PropId);

  /* Restore previous error handler */
  H5Eset_auto1(old_func, old_client_data);

#ifndef XMDF_NO_CREATE_GROUP_ON_OPEN
  /* If group does not exist create it */
  if (error < 0) {
    *a_PropId = xfpCreateGroup(a_Id, MESH_PROP_GROUP, GROUP_TYPE_PROPERTIES);
    error = *a_PropId;
  }
#endif

  return error;
} /* xfGetMeshPropertyGroup */
/* -------------------------------------------------------------------------- */
/* FUNCTION  xfCreateMeshNodePropertyGroup */
/*! PURPOSE:
 * - NOTES:      */
/* -------------------------------------------------------------------------- */
XMDF_API xid xfCreateMeshNodePropertyGroup (xid a_Id, xid *a_PropId)
{
  int     error = 1;
  herr_t (*old_func)(void*);
  void    *old_client_data;

  H5Eget_auto1(&old_func, &old_client_data);

  /* Turn off error handling */
  H5Eset_auto1(NULL, NULL);

  /* Open group */
  error = xfOpenGroup(a_Id, MESH_NODE_PROP_GROUP, a_PropId);

  /* Restore previous error handler */
  H5Eset_auto1(old_func, old_client_data);

  /* If group does not exist create it */
  if (error < 0) {
    *a_PropId = xfpCreateGroup(a_Id, MESH_NODE_PROP_GROUP, GROUP_TYPE_PROPERTIES);
    error = *a_PropId;
  }

  return error;
} /* xfCreateMeshNodePropertyGroup */
/* -------------------------------------------------------------------------- */
/* FUNCTION  xfGetMeshNodePropertyGroup */
/*! PURPOSE:
 * - NOTES:      */
/* -------------------------------------------------------------------------- */
XMDF_API xid xfGetMeshNodePropertyGroup (xid a_Id, xid *a_PropId)
{
  int     error = 1;
  herr_t (*old_func)(void*);
  void    *old_client_data;

  H5Eget_auto1(&old_func, &old_client_data);

  /* Turn off error handling */
  H5Eset_auto1(NULL, NULL);

  /* Open group */
  error = xfOpenGroup(a_Id, MESH_NODE_PROP_GROUP, a_PropId);

  /* Restore previous error handler */
  H5Eset_auto1(old_func, old_client_data);

#ifndef XMDF_NO_CREATE_GROUP_ON_OPEN
  /* If group does not exist create it */
  if (error < 0) {
    *a_PropId = xfpCreateGroup(a_Id, MESH_NODE_PROP_GROUP, GROUP_TYPE_PROPERTIES);
    error = *a_PropId;
  }
#endif

  return error;
} /* xfGetMeshNodePropertyGroup */
/* -------------------------------------------------------------------------- */
/* FUNCTION  xfCreateMeshElementPropertyGroup */
/*! PURPOSE:
 * - NOTES:      */
/* -------------------------------------------------------------------------- */
XMDF_API xid xfCreateMeshElementPropertyGroup (xid a_Id, xid *a_PropId)
{
  int     error = 1;
  herr_t (*old_func)(void*);
  void    *old_client_data;

  H5Eget_auto1(&old_func, &old_client_data);

  /* Turn off error handling */
  H5Eset_auto1(NULL, NULL);

  /* Open group */
  error = xfOpenGroup(a_Id, MESH_ELEM_PROP_GROUP, a_PropId);
 
  /* Restore previous error handler */
  H5Eset_auto1(old_func, old_client_data);

  /* If group does not exist create it */
  if (error < 0) {
    *a_PropId = xfpCreateGroup(a_Id, MESH_ELEM_PROP_GROUP, GROUP_TYPE_PROPERTIES);
    error = *a_PropId;
  }

  return error;
} /* xfCreateGridMeshElementPropertyGroup */
/* -------------------------------------------------------------------------- */
/* FUNCTION  xfGetMeshElementPropertyGroup */
/*! PURPOSE:
 * - NOTES:      */
/* -------------------------------------------------------------------------- */
XMDF_API xid xfGetMeshElementPropertyGroup (xid a_Id, xid *a_PropId)
{
  int     error = 1;
  herr_t (*old_func)(void*);
  void    *old_client_data;

  H5Eget_auto1(&old_func, &old_client_data);

  /* Turn off error handling */
  H5Eset_auto1(NULL, NULL);

  /* Open group */
  error = xfOpenGroup(a_Id, MESH_ELEM_PROP_GROUP, a_PropId);
 
  /* Restore previous error handler */
  H5Eset_auto1(old_func, old_client_data);

#ifndef XMDF_NO_CREATE_GROUP_ON_OPEN
  /* If group does not exist create it */
  if (error < 0) {
    *a_PropId = xfpCreateGroup(a_Id, MESH_ELEM_PROP_GROUP, GROUP_TYPE_PROPERTIES);
    error = *a_PropId;
  }
#endif

  return error;
} /* xfGetGridMeshElementPropertyGroup */
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfSetGridType */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfSetGridType (xid a_Id, int a_GridType)
{
  int error = 1;

  if (a_GridType < GRID_TYPE_MIN || a_GridType > GRID_TYPE_MAX) {
    return ERROR_GRID_TYPE_INVALID;
  }

  error = xfpWriteAttributeInt(a_Id, GRID_ATT_TYPE, 1, &a_GridType);

  return error;
} /* xfSetGridType */
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfSetNumberOfDimensions */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfSetNumberOfDimensions (xid a_Id, int a_NumDimensions)
{
  int error = 1;

  if (a_NumDimensions < 2 || a_NumDimensions > 3) {
    return ERROR_GRID_NUM_DIMENSIONS_INVALID;
  }

  error = xfpWriteAttributeInt(a_Id, GRID_ATT_NUM_DIMENSIONS, 1, &a_NumDimensions);

  return error;
} /* xfSetNumberOfDimensions */
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfSetExtrusionType */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfSetExtrusionType (xid a_Id, int a_ExtrudeType)
{
  int error = 1;

  if (a_ExtrudeType < EXTRUDE_MIN || a_ExtrudeType > EXTRUDE_MAX) {
    return ERROR_GRID_TYPE_INVALID;
  }

  error = xfpWriteAttributeInt(a_Id, GRID_ATT_EXTRUDE_TYPE, 1, &a_ExtrudeType);

  return error;
} /* xfSetExtrusionType */
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfSetOrigin */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfSetOrigin (xid a_Id, double a_x, double a_y, double a_z)
{
  double    loc[3];
  int       error;

  loc[0] = a_x;
  loc[1] = a_y;
  loc[2] = a_z;

  error = xfpWriteAttributeDouble(a_Id, GRID_ATT_ORIGIN, 3, loc);

  return error;
} /* xfSetOrigin */
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfSetOrientation */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfSetOrientation (xid a_Id, int a_Orientation)
{
  int error = 1;

  if (a_Orientation != ORIENTATION_RIGHT_HAND &&
      a_Orientation != ORIENTATION_LEFT_HAND) {
    return ERROR_GRID_TYPE_INVALID;
  }

  error = xfpWriteAttributeInt(a_Id, GRID_ATT_ORIENTATION, 1, &a_Orientation);

  return error;
} /* xfSetOrientation */
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfSetBearing */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfSetBearing (xid a_Id, double a_Bearing)
{
  int error = 1;

  error = xfpWriteAttributeDouble(a_Id, GRID_ATT_BEARING, 1, &a_Bearing);

  return error;
} /* xfSetBearing */
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfSetDip */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfSetDip (xid a_Id, double a_Dip)
{
  int error = 1;

  error = xfpWriteAttributeDouble(a_Id, GRID_ATT_DIP, 1, &a_Dip);

  return error;
} /* xfSetDip */
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfSetRoll */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfSetRoll (xid a_Id, double a_Roll)
{
  int error = 1;

  error = xfpWriteAttributeDouble(a_Id, GRID_ATT_ROLL, 1, &a_Roll);

  return error;
} /* xfSetRoll */
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfSetComputationalOrigin */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfSetComputationalOrigin (xid a_Id, int a_Origin)
{
  int error = 1;

  error = xfpWriteAttributeInt(a_Id, GRID_ATT_COMP_ORIG, 1, &a_Origin);

  return error;
} /* xfSetComputationalOrigin */
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfSetUDirection */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfSetUDirection (xid a_Id, int a_Direction)
{
  int error = 1;

  error = xfpWriteAttributeInt(a_Id, GRID_ATT_U_DIR, 1, &a_Direction);

  return error;
} /* xfSetUDirection */
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfSetNumberCellsInI */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfSetNumberCellsInI (xid a_Id, int a_NumI)
{
  int error = 1;

  error = xfpWriteAttributeInt(a_Id, GRID_ATT_NUM_I, 1, &a_NumI);

  return error;
} /* xfSetNumberCellsInI */
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfSetNumberCellsInJ */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfSetNumberCellsInJ (xid a_Id, int a_NumJ)
{
  int error = 1;

  error = xfpWriteAttributeInt(a_Id, GRID_ATT_NUM_J, 1, &a_NumJ);

  return error;
} /* xfSetNumberCellsInJ */
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfSetNumberCellsInK */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfSetNumberCellsInK (xid a_Id, int a_NumK)
{
  int error = 1;

  error = xfpWriteAttributeInt(a_Id, GRID_ATT_NUM_K, 1, &a_NumK);

  return error;
} /* xfSetNumberCellsInK */
/*-----------------------------------------------------------------------------*/
/* FUNCTION  xfSetGridCoordsI*/
/*! PURPOSE:
 * - NOTES:*/
/*-----------------------------------------------------------------------------*/
XMDF_API xid xfSetGridCoordsI (xid a_Id, int a_NumVals, double *a_iValues)
{
  int         error = 1;
  XDatasetParams Params;

  xfpDsetParamsInit(&Params, 1, 0, NONE);
  xfpDsetParamsSetSizes(&Params, 0, a_NumVals, a_NumVals, a_NumVals);

  error = xfpWriteDsetDouble(a_Id, GRID_DSET_COORDS_I, &Params, a_iValues);
  xfpDsetParamsDestroy(&Params);

  return error;
} /* xfSetGridCoordsI */
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfSetGridCoordsJ */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfSetGridCoordsJ (xid a_Id, int a_NumVals, double *a_jValues)
{
  int         error = 1;
  XDatasetParams Params;

  xfpDsetParamsInit(&Params, 1, 0, NONE);
  xfpDsetParamsSetSizes(&Params, 0, a_NumVals, a_NumVals, a_NumVals);

  error = xfpWriteDsetDouble(a_Id, GRID_DSET_COORDS_J, &Params, a_jValues);
  xfpDsetParamsDestroy(&Params);

  return error;
} /* xfSetGridCoordsJ */
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfSetGridCoordsK */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfSetGridCoordsK (xid a_Id, int a_NumVals, double *a_kValues)
{
  int         error = 1;
  XDatasetParams Params;
  int            nDims;

  /* The number of grid dimensions should already have been set to 3 */
  error = xfGetNumberOfDimensions(a_Id, &nDims);
  if (error < 0) {
    return error;
  }
  if (nDims != 3) {
    return ERROR_GRID_NUM_DIMENSIONS_INVALID;
  }

  xfpDsetParamsInit(&Params, 1, 0, NONE);
  xfpDsetParamsSetSizes(&Params, 0, a_NumVals, a_NumVals, a_NumVals);

  error = xfpWriteDsetDouble(a_Id, GRID_DSET_COORDS_K, &Params, a_kValues);
  xfpDsetParamsDestroy(&Params);

  return error;
} /* xfSetGridCoordsK */
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfWriteExtrudeLayerData */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfWriteExtrudeLayerData (xid a_Id, int a_NumLayers, int a_NumVals,
                                     double *a_Values)
{
  int    error;
  int    NumI = 0, NumJ = 0, nValsCorrect = 0;
  XDatasetParams Params;
  int    ExtrudeType = 0;

  error = xfGetExtrusionType(a_Id, &ExtrudeType);
  if (error < 0) {
    return error;
  }

  /* First make sure that the number of layers and number of values are correct */
  switch (ExtrudeType) {
    case EXTRUDE_SIGMA:
      /* The number of layers and number of values should be the same for SIGMA extruded */
      /* grids */
      if (a_NumLayers != a_NumVals) {
        return ERROR_GRID_NUMVALS_INCORRECT;
      }
      break;
    case EXTRUDE_CURV_AT_CORNERS:
      error = xfGetNumberCellsInI(a_Id, &NumI);
      if (error <= 0) {
        error = xfGetNumberCellsInJ(a_Id, &NumJ);
      }
      if (error < 0) {
        return error;
      }
      /* Compute the correct number of values */
      nValsCorrect = (NumI + 1) * (NumJ + 1) * (a_NumLayers * 1); 
      if (nValsCorrect != a_NumVals) {
        return ERROR_GRID_NUMVALS_INCORRECT;
      }
      break;
    case EXTRUDE_CURV_AT_CELLS:
      error = xfGetNumberCellsInI(a_Id, &NumI);
      if (error <= 0) {
        error = xfGetNumberCellsInJ(a_Id, &NumJ);
      }
      if (error < 0) {
        return error;
      }
      /* Compute the correct number of values */
      nValsCorrect = NumI * NumJ * a_NumLayers; 
      if (nValsCorrect != a_NumVals) {
        return ERROR_GRID_NUMVALS_INCORRECT;
      }
      break;
    default:
      return ERROR_GRID_EXTRUDE_TYPE_INVALID;
      break;
  }

  /* Write the number of extrude layer attribute */
  xfpWriteAttributeInt(a_Id, GRID_ATT_EXTRUDE_LAYERS, 1, &a_NumLayers);

  /* Write the values */
  xfpDsetParamsInit(&Params, 1, 0, NONE);
  xfpDsetParamsSetSizes(&Params, 0, a_NumVals, a_NumVals, a_NumVals);

  error = xfpWriteDsetDouble(a_Id, GRID_DSET_EXTRUDE_VALS, &Params, a_Values);
  xfpDsetParamsDestroy(&Params);

  return error;
} /* xfWriteExtrudeLayerData */
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfGetGridType */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfGetGridType (xid a_Id, int *a_GridType)
{
  int error = 1;

  error = xfpReadAttributeInt(a_Id, GRID_ATT_TYPE, 1, a_GridType);

  return error;
} /* xfGetGridType */
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfGetExtrusionType */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfGetExtrusionType (xid a_Id, int *a_ExtrudeType)
{
  int error = 1;

  error = xfpReadAttributeInt(a_Id, GRID_ATT_EXTRUDE_TYPE, 1, a_ExtrudeType);

  return error;
} /* xfGetExtrusionType */
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfGetNumberOfDimensions */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfGetNumberOfDimensions (xid a_Id, int *a_NumDimensions)
{
  int error = 1;

  error = xfpReadAttributeInt(a_Id, GRID_ATT_NUM_DIMENSIONS, 1, a_NumDimensions);

  return error;
} /* xfGetNumberOfDimensions */
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfOriginDefined */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfOriginDefined (xid a_Id, xmbool *a_bDefined)
{
  int  error = 1;

  error = xfpDoesAttributeExist(a_Id, GRID_ATT_ORIGIN, a_bDefined);

  return error;
} /* xfOriginDefined */
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfGetOrigin */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfGetOrigin (xid a_Id, double *a_x, double *a_y, double *a_z)
{
  int  error = 1;
  double loc[3];

  error = xfpReadAttributeDouble(a_Id, GRID_ATT_ORIGIN, 3, loc);
  if (error < 0) {
    return error;
  }

  *a_x = loc[0];
  *a_y = loc[1];
  *a_z = loc[2];

  return error;
} /* xfGetOrigin */
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfGetOrientation */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfGetOrientation (xid a_Id, int *a_Orientation)
{
  int error = 1;

  error = xfpReadAttributeInt(a_Id, GRID_ATT_ORIENTATION, 1, a_Orientation);

  return error;
} /* xfGetOrientation */
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfBearingDefined */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfBearingDefined (xid a_Id, xmbool *a_bDefined)
{
  int error = 1;

  error = xfpDoesAttributeExist(a_Id, GRID_ATT_BEARING, a_bDefined);

  return error;
} /* xfBearingDefined */
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfGetBearing */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfGetBearing (xid a_Id, double *a_bearing)
{
  int error = 1;

  error = xfpReadAttributeDouble(a_Id, GRID_ATT_BEARING, 1, a_bearing);

  return error;
} /* xfGetBearing */
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfDipDefined */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfDipDefined (xid a_Id, xmbool *a_bDefined)
{
  int error = 1;

  error = xfpDoesAttributeExist(a_Id, GRID_ATT_DIP, a_bDefined);

  return error;
} /* xfDipDefined */
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfGetDip */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfGetDip (xid a_Id, double *a_dip)
{
  int error = 1;

  error = xfpReadAttributeDouble(a_Id, GRID_ATT_DIP, 1, a_dip);

  return error;
} /* xfGetDip */
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfRollDefined */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfRollDefined (xid a_Id, xmbool *a_bDefined)
{
  int error = 1;

  error = xfpDoesAttributeExist(a_Id, GRID_ATT_ROLL, a_bDefined);

  return error;
} /* xfRollDefined */
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfGetRoll */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfGetRoll (xid a_Id, double *a_Roll)
{
  int error = 1;

  error = xfpReadAttributeDouble(a_Id, GRID_ATT_ROLL, 1, a_Roll);

  return error;
} /* xfGetRoll */
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfComputationalOriginDefined */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfComputationalOriginDefined (xid GroupId, xmbool *bDefined)
{
  int error = 1;

  error = xfpDoesAttributeExist(GroupId, GRID_ATT_COMP_ORIG, bDefined);

  return error;
} /* xfComputationalOriginDefined */
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfGetComputationalOrigin */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfGetComputationalOrigin (xid GroupId, int *Origin)
{
  int error = 1;

  error = xfpReadAttributeInt(GroupId, GRID_ATT_COMP_ORIG, 1, Origin);

  return error;
} /* xfGetComputationalOrigin */
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfGetUDirectionDefined */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfGetUDirectionDefined (xid GroupId, xmbool *bDefined)
{
  int error = 1;

  error = xfpDoesAttributeExist(GroupId, GRID_ATT_U_DIR, bDefined);

  return error;
} /* xfGetUDirectionDefined */
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfGetUDirection */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfGetUDirection (xid GroupId, int *Direction)
{
  int error = 1;

  error = xfpReadAttributeInt(GroupId, GRID_ATT_U_DIR, 1, Direction);

  return error;
} /* xfGetUDirection */
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfGetNumberCellsInI */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfGetNumberCellsInI (xid a_Id, int *a_NumI)
{
  int error = 1;

  error = xfpReadAttributeInt(a_Id, GRID_ATT_NUM_I, 1, a_NumI);

  return error;
} /* xfGetNumberCellsInI */
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfGetNumberCellsInJ */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfGetNumberCellsInJ (xid a_Id, int *a_NumJ)
{
  int error = 1;

  error = xfpReadAttributeInt(a_Id, GRID_ATT_NUM_J, 1, a_NumJ);

  return error;
} /* xfGetNumberCellsInJ */
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfGetNumberCellsInK */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfGetNumberCellsInK (xid a_Id, int *a_NumK)
{
  int error = 1;

  error = xfpReadAttributeInt(a_Id, GRID_ATT_NUM_K, 1, a_NumK);

  return error;
} /* xfGetNumberCellsInK */
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfGetGridCoordsI */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfGetGridCoordsI (xid a_Id, int a_NumVals, double *a_iValues)
{
  int error = 1;

  error = xfpReadDset1DDouble(a_Id, GRID_DSET_COORDS_I, a_NumVals, a_iValues);

  return error;
} /* xfGetGridCoordsI */
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfGetGridCoordsJ */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfGetGridCoordsJ (xid a_Id, int a_NumVals, double *a_jValues)
{
  int error = 1;

  error = xfpReadDset1DDouble(a_Id, GRID_DSET_COORDS_J, a_NumVals, a_jValues);

  return error;
} /* xfGetGridCoordsJ */
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfGetGridCoordsK */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfGetGridCoordsK (xid a_Id, int a_NumVals, double *a_kValues)
{
  int error = 1;

  error = xfpReadDset1DDouble(a_Id, GRID_DSET_COORDS_K, a_NumVals, a_kValues);

  return error;
} /* xfGetGridCoordsK */
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfGetExtrudeNumLayers */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfGetExtrudeNumLayers (xid a_Id, int *a_NumLayers)
{
  int  error = 1;

  error = xfpReadAttributeInt(a_Id, GRID_ATT_EXTRUDE_LAYERS, 1, a_NumLayers);

  return error;
} /* xfGetExtrudeNumLayers */
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfGetExtrudeValues */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfGetExtrudeValues (xid a_Id, int a_NumVals, double *a_Values)
{
  int error = 1;

  error = xfpReadDset1DDouble(a_Id, GRID_DSET_EXTRUDE_VALS, a_NumVals, a_Values);

  return error;
} /* xfGetExtrudeValues */
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfCreateGridPropertyGroup                                         */
/*! PURPOSE:
 * - NOTES:                                                                       */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfCreateGridPropertyGroup (xid a_Id, xid *a_PropId)
{
  int     error = 1;
  herr_t (*old_func)(void*);
  void    *old_client_data;

  H5Eget_auto1(&old_func, &old_client_data);

  /* Turn off error handling */
  H5Eset_auto1(NULL, NULL);

  /* Open group */
  error = xfOpenGroup(a_Id, GRID_PROP_GROUP, a_PropId);

  /* Restore previous error handler */
  H5Eset_auto1(old_func, old_client_data);

  /* If group does not exist create it */
  if (error < 0) {
    *a_PropId = xfpCreateGroup(a_Id, GRID_PROP_GROUP, GROUP_TYPE_PROPERTIES);
    error = *a_PropId;
  }

  return error;
} /* xfCreateGridPropertyGroup */
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfGetGridPropertyGroup                                            */
/*! PURPOSE:
 * - NOTES:                                                                       */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfGetGridPropertyGroup (xid a_Id, xid *a_PropId)
{
  int     error = 1;
  herr_t (*old_func)(void*);
  void    *old_client_data;

  H5Eget_auto1(&old_func, &old_client_data);

  /* Turn off error handling */
  H5Eset_auto1(NULL, NULL);

  /* Open group */
  error = xfOpenGroup(a_Id, GRID_PROP_GROUP, a_PropId);

  /* Restore previous error handler */
  H5Eset_auto1(old_func, old_client_data);

#ifndef XMDF_NO_CREATE_GROUP_ON_OPEN
  /* If group does not exist create it */
  if (error < 0) {
    *a_PropId = xfpCreateGroup(a_Id, GRID_PROP_GROUP, GROUP_TYPE_PROPERTIES);
    error = *a_PropId;
  }
#endif

  return error;
} /* xfGetGridPropertyGroup */
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfCreateGridCellPropertyGroup                                     */
/*! PURPOSE:
 * - NOTES:                                                                       */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfCreateGridCellPropertyGroup (xid a_Id, xid *a_PropId)
{
  int     error = 1;
  herr_t (*old_func)(void*);
  void    *old_client_data;

  H5Eget_auto1(&old_func, &old_client_data);

  /* Turn off error handling */
  H5Eset_auto1(NULL, NULL);

  /* Open group */
  error = xfOpenGroup(a_Id, GRID_CELL_PROP_GROUP, a_PropId);

  /* Restore previous error handler */
  H5Eset_auto1(old_func, old_client_data);

  /* If group does not exist create it */
  if (error < 0) {
    *a_PropId = xfpCreateGroup(a_Id, GRID_CELL_PROP_GROUP, GROUP_TYPE_PROPERTIES);
    error = *a_PropId;
  }

  return error;
} /* xfCreateGridCellPropertyGroup */
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfGetGridCellPropertyGroup                                        */
/*! PURPOSE:
 * - NOTES:                                                                       */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfGetGridCellPropertyGroup (xid a_Id, xid *a_PropId)
{
  int     error = 1;
  herr_t (*old_func)(void*);
  void    *old_client_data;

  H5Eget_auto1(&old_func, &old_client_data);

  /* Turn off error handling */
  H5Eset_auto1(NULL, NULL);

  /* Open group */
  error = xfOpenGroup(a_Id, GRID_CELL_PROP_GROUP, a_PropId);

  /* Restore previous error handler */
  H5Eset_auto1(old_func, old_client_data);

#ifndef XMDF_NO_CREATE_GROUP_ON_OPEN
  /* If group does not exist create it */
  if (error < 0) {
    *a_PropId = xfpCreateGroup(a_Id, GRID_CELL_PROP_GROUP, GROUP_TYPE_PROPERTIES);
    error = *a_PropId;
  }
#endif

  return error;
} /* xfGetGridCellPropertyGroup */
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfCreateGridNodePropertyGroup                                     */
/*! PURPOSE:
 * - NOTES:                                                                       */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfCreateGridNodePropertyGroup (xid a_Id, xid *a_PropId)
{
  int     error = 1;
  herr_t (*old_func)(void*);
  void    *old_client_data;

  H5Eget_auto1(&old_func, &old_client_data);

  /* Turn off error handling */
  H5Eset_auto1(NULL, NULL);

  /* Open group */
  error = xfOpenGroup(a_Id, GRID_NODE_PROP_GROUP, a_PropId);
 
  /* Restore previous error handler */
  H5Eset_auto1(old_func, old_client_data);

  /* If group does not exist create it */
  if (error < 0) {
    *a_PropId = xfpCreateGroup(a_Id, GRID_NODE_PROP_GROUP, GROUP_TYPE_PROPERTIES);
    error = *a_PropId;
  }

  return error;
} /* xfCreateGridNodePropertyGroup */
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfGetGridNodePropertyGroup                                        */
/*! PURPOSE:
 * - NOTES:                                                                       */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfGetGridNodePropertyGroup (xid a_Id, xid *a_PropId)
{
  int     error = 1;
  herr_t (*old_func)(void*);
  void    *old_client_data;

  H5Eget_auto1(&old_func, &old_client_data);

  /* Turn off error handling */
  H5Eset_auto1(NULL, NULL);

  /* Open group */
  error = xfOpenGroup(a_Id, GRID_NODE_PROP_GROUP, a_PropId);
 
  /* Restore previous error handler */
  H5Eset_auto1(old_func, old_client_data);

#ifndef XMDF_NO_CREATE_GROUP_ON_OPEN
  /* If group does not exist create it */
  if (error < 0) {
    *a_PropId = xfpCreateGroup(a_Id, GRID_NODE_PROP_GROUP, GROUP_TYPE_PROPERTIES);
    error = *a_PropId;
  }
#endif

  return error;
} /* xfGetGridNodePropertyGroup */

/* --------------------------------------------------------------------------- */
/* FUNCTION  xfSetNumberOfXSects */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfSetNumberOfXSects (xid a_Id, int *a_nXSects, int a_compression)
{
  xid XSectGroup, XSecGeometryGroup, LnPropGroup, PtPropGroup;
  xid XSecGeometrySubGroup, LnPropSubGroup, PtPropSubGroup;
  int status, i;
  char grpname[30];

  /* open the XSect group */
  status = xfOpenGroup(a_Id, XSEC_XSECS_GROUP, &XSectGroup);
  if (status < 0) {
    return status;
  }
  status = xfWritePropertyInt(XSectGroup, XSEC_XSECS_GROUP, 1, a_nXSects, 
                              a_compression);
  xfCloseGroup(XSectGroup);

  /* open the Geometry group */
  status = xfOpenGroup(a_Id, XSEC_GEOMETRY_GROUP, &XSecGeometryGroup);
  for (i=0; i<(*a_nXSects); i++) {
    sprintf(grpname, "%s%d", XSEC_GEOMETRY_GROUP, i);
    XSecGeometrySubGroup = xfpCreateGroup(XSecGeometryGroup, grpname, 
                                          GROUP_TYPE_XSEC_GEOM);
    xfCloseGroup(XSecGeometrySubGroup);
  }
  xfCloseGroup(XSecGeometryGroup);

  /* open the Line group */
  status = xfOpenGroup(a_Id, XSEC_LINE_GROUP, &LnPropGroup);
  for (i=0; i<(*a_nXSects); i++) {
    sprintf(grpname, "%s%d", XSEC_LINE_GROUP, i);
    LnPropSubGroup = xfpCreateGroup(LnPropGroup, grpname, 
            GROUP_TYPE_XSEC_LINE);
    xfCloseGroup(LnPropSubGroup);
    if (grpname == NULL) {
      free(grpname);
    }
  }
  xfCloseGroup(LnPropGroup);

  /* open the Point group */
  status = xfOpenGroup(a_Id, XSEC_POINT_GROUP, &PtPropGroup);
  for (i=0; i<(*a_nXSects); i++) {
    sprintf(grpname, "%s%d", XSEC_POINT_GROUP, i);
    PtPropSubGroup = xfpCreateGroup(PtPropGroup, grpname, 
          GROUP_TYPE_XSEC_POINT);
    xfCloseGroup(PtPropSubGroup);
  }
  xfCloseGroup(PtPropGroup);

  return status;
} /* xfSetNumberOfXSects */


/* --------------------------------------------------------------------------- */
/* FUNCTION  xfGetNumberOfXSects */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfGetNumberOfXSects (xid a_Id, int *a_nXSects)
{
  xid XSectGroup;
  int status;
  xmbool exists;

  /* open the XSect group */
  status = xfOpenGroup(a_Id, XSEC_XSECS_GROUP, &XSectGroup);
  if (status < 0) {
    return status;
  }
  xfDoesPropertyWithNameExist(XSectGroup, XSEC_XSECS_GROUP, &exists);
  if (exists) {
    xfReadPropertyInt(XSectGroup, XSEC_XSECS_GROUP, 1, a_nXSects);
  }

  xfCloseGroup(XSectGroup);

  return status;
} /* xfGetNumberOfXSects */

/* --------------------------------------------------------------------------- */
/* FUNCTION  xfSetCSID */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfSetCSID (xid a_Id, int a_NumVals, int *a_PropId,
                        int a_compression)
{
  xid XSectGroup;
  int status;

  /* open the XSect group */
  status = xfOpenGroup(a_Id, XSEC_XSECS_GROUP, &XSectGroup);
  if (status < 0) {
    return status;
  }
  status = xfWritePropertyInt(XSectGroup, XMDF_CSID, a_NumVals, a_PropId, 
                              a_compression);
  xfCloseGroup(XSectGroup);

  return status;
} /* xfSetCSID */


/* --------------------------------------------------------------------------- */
/* FUNCTION  xfGetCSID */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfGetCSID (xid a_Id,  int *a_NumVals, int *a_PropId)
{
  xid XSectGroup;
  int status;
  xmbool exists;

  /* open the XSect group */
  status = xfOpenGroup(a_Id, XSEC_XSECS_GROUP, &XSectGroup);
  if (status < 0) {
    return status;
  }
  xfDoesPropertyWithNameExist(XSectGroup, XMDF_CSID, &exists);
  if (exists) {
    xfGetPropertyNumber(XSectGroup, XMDF_CSID, a_NumVals);
    xfReadPropertyInt(XSectGroup, XMDF_CSID, (*a_NumVals), a_PropId);
  }

  xfCloseGroup(XSectGroup);

  return status;
} /* xfGetCSID */


/* --------------------------------------------------------------------------- */
/* FUNCTION  xfSetCSName */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfSetCSName (xid a_Id, int a_NumVals, int a_StrLen, 
                          char *a_PropId)
{
  xid XSectGroup;
  int status;

  /* open the XSect group */
  status = xfOpenGroup(a_Id, XSEC_XSECS_GROUP, &XSectGroup);
  if (status < 0) {
    return status;
  }
  status = xfWritePropertyString(XSectGroup, XMDF_CSNAME, a_NumVals, 
                                 a_StrLen, a_PropId);

  xfCloseGroup(XSectGroup);

  return status;
} /* xfSetCSName */


/* --------------------------------------------------------------------------- */
/* FUNCTION  xfGetCSName */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfGetCSName (xid a_Id, int *a_NumVals, int *a_StrLen, 
                          char *a_PropId)
{
  xid XSectGroup;
  int status;
  xmbool exists;

  /* open the XSect group */
  status = xfOpenGroup(a_Id, XSEC_XSECS_GROUP, &XSectGroup);
  if (status < 0) {
    return status;
  }
  xfDoesPropertyWithNameExist(XSectGroup, XMDF_CSNAME, &exists);
  if (exists) {
    xfGetPropertyNumber(XSectGroup, XMDF_CSNAME, a_NumVals);
    xfGetPropertyStringLength(XSectGroup, XMDF_CSNAME, a_NumVals, a_StrLen);
    xfReadPropertyString(XSectGroup, XMDF_CSNAME, (*a_NumVals), (*a_StrLen),
                         a_PropId);
  }

  xfCloseGroup(XSectGroup);

  return status;
} /* xfGetCSName */


/* --------------------------------------------------------------------------- */
/* FUNCTION  xfGetCSNameLen */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfGetCSNameLen (xid a_Id, int *a_NumVals, int *a_StrLen)
{
  xid XSectGroup;
  int status;
  xmbool exists;

  /* open the XSect group */
  status = xfOpenGroup(a_Id, XSEC_XSECS_GROUP, &XSectGroup);
  if (status < 0) {
    return status;
  }
  xfDoesPropertyWithNameExist(XSectGroup, XMDF_CSNAME, &exists);
  if (exists) {
    xfGetPropertyNumber(XSectGroup, XMDF_CSNAME, a_NumVals);
    xfGetPropertyStringLength(XSectGroup, XMDF_CSNAME, a_NumVals, a_StrLen);
  }

  xfCloseGroup(XSectGroup);

  return status;
} /* xfGetCSNameLen */

/* --------------------------------------------------------------------------- */
/* FUNCTION  xfSetReachName */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfSetReachName (xid a_Id, int a_NumVals, int a_StrLen, 
                             char *a_PropId)
{
  xid XSectGroup;
  int status;

  /* open the XSect group */
  status = xfOpenGroup(a_Id, XSEC_XSECS_GROUP, &XSectGroup);
  if (status < 0) {
    return status;
  }
  status = xfWritePropertyString(XSectGroup, XMDF_REACHNAME, a_NumVals, 
                                 a_StrLen, a_PropId);

  xfCloseGroup(XSectGroup);

  return status;
} /* xfSetReachName */


/* --------------------------------------------------------------------------- */
/* FUNCTION  xfGetReachName */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfGetReachName (xid a_Id, int *a_NumVals, int *a_StrLen, 
                             char *a_PropId)
{
  xid XSectGroup;
  int status;
  xmbool exists;

  /* open the XSect group */
  status = xfOpenGroup(a_Id, XSEC_XSECS_GROUP, &XSectGroup);
  if (status < 0) {
    return status;
  }
  xfDoesPropertyWithNameExist(XSectGroup, XMDF_REACHNAME, &exists);
  if (exists) {
    xfGetPropertyNumber(XSectGroup, XMDF_REACHNAME, a_NumVals);
    xfGetPropertyStringLength(XSectGroup, XMDF_REACHNAME, a_NumVals, a_StrLen);
    xfReadPropertyString(XSectGroup, XMDF_REACHNAME, (*a_NumVals), (*a_StrLen),
                         a_PropId);
  }

  xfCloseGroup(XSectGroup);

  return status;
} /* xfGetReachName */

/* --------------------------------------------------------------------------- */
/* FUNCTION  xfGetReachNameLen */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfGetReachNameLen (xid a_Id, int *a_NumVals, int *a_StrLen)
{
  xid XSectGroup;
  int status;
  xmbool exists;

  /* open the XSect group */
  status = xfOpenGroup(a_Id, XSEC_XSECS_GROUP, &XSectGroup);
  if (status < 0) {
    return status;
  }
  xfDoesPropertyWithNameExist(XSectGroup, XMDF_REACHNAME, &exists);
  if (exists) {
    xfGetPropertyNumber(XSectGroup, XMDF_REACHNAME, a_NumVals);
    xfGetPropertyStringLength(XSectGroup, XMDF_REACHNAME, a_NumVals, a_StrLen);
  }

  xfCloseGroup(XSectGroup);

  return status;
} /* xfGetReachNameLen */

/* --------------------------------------------------------------------------- */
/* FUNCTION  xfSetTopoID */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfSetTopoID (xid a_Id, int a_NumVals, int a_StrLen, 
                          char *a_PropId)
{
  xid XSectGroup;
  int status;

  /* open the XSect group */
  status = xfOpenGroup(a_Id, XSEC_XSECS_GROUP, &XSectGroup);
  if (status < 0) {
    return status;
  }
  status = xfWritePropertyString(XSectGroup, XMDF_TOPO_ID, a_NumVals, a_StrLen,
                                 a_PropId);

  xfCloseGroup(XSectGroup);

  return status;
} /* xfSetTopoID */


/* --------------------------------------------------------------------------- */
/* FUNCTION  xfGetTopoID */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfGetTopoID (xid a_Id, int *a_NumVals, int *a_StrLen, 
                          char *a_PropId)
{
  xid XSectGroup;
  int status;
  xmbool exists;

  /* open the XSect group */
  status = xfOpenGroup(a_Id, XSEC_XSECS_GROUP, &XSectGroup);
  if (status < 0) {
    return status;
  }
  xfDoesPropertyWithNameExist(XSectGroup, XMDF_TOPO_ID, &exists);
  if (exists) {
    xfGetPropertyNumber(XSectGroup, XMDF_TOPO_ID, a_NumVals);
    xfGetPropertyStringLength(XSectGroup, XMDF_TOPO_ID, a_NumVals, a_StrLen);
    xfReadPropertyString(XSectGroup, XMDF_TOPO_ID, (*a_NumVals), (*a_StrLen), 
                         a_PropId);
  }

  xfCloseGroup(XSectGroup);

  return status;
} /* xfGetTopoID */

/* --------------------------------------------------------------------------- */
/* FUNCTION  xfGetTopoIDLen */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfGetTopoIDLen (xid a_Id, int *a_NumVals, int *a_StrLen)
{
  xid XSectGroup;
  int status;
  xmbool exists;

  /* open the XSect group */
  status = xfOpenGroup(a_Id, XSEC_XSECS_GROUP, &XSectGroup);
  if (status < 0) {
    return status;
  }
  xfDoesPropertyWithNameExist(XSectGroup, XMDF_TOPO_ID, &exists);
  if (exists) {
    xfGetPropertyNumber(XSectGroup, XMDF_TOPO_ID, a_NumVals);
    xfGetPropertyStringLength(XSectGroup, XMDF_TOPO_ID, a_NumVals, a_StrLen);
  }

  xfCloseGroup(XSectGroup);

  return status;
} /* xfGetTopoIDLen */

/* --------------------------------------------------------------------------- */
/* FUNCTION  xfSetStation */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfSetStation (xid a_Id, int a_NumVals, double *a_PropId, 
                           int a_compression)
{
  xid XSectGroup;
  int status;

  /* open the XSect group */
  status = xfOpenGroup(a_Id, XSEC_XSECS_GROUP, &XSectGroup);
  if (status < 0) {
    return status;
  }
  status = xfWritePropertyDouble(XSectGroup, XMDF_STATION, a_NumVals, a_PropId,
                                 a_compression);
  xfCloseGroup(XSectGroup);

  return status;
} /* xfSetStation */


/* --------------------------------------------------------------------------- */
/* FUNCTION  xfGetStation */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfGetStation (xid a_Id, int *a_NumVals, double *a_PropId)
{
  xid XSectGroup;
  int status;
  xmbool exists;

  /* open the XSect group */
  status = xfOpenGroup(a_Id, XSEC_XSECS_GROUP, &XSectGroup);
  if (status < 0) {
    return status;
  }
  xfDoesPropertyWithNameExist(XSectGroup, XMDF_STATION, &exists);
  if (exists) {
    xfGetPropertyNumber(XSectGroup, XMDF_STATION, a_NumVals);
    xfReadPropertyDouble(XSectGroup, XMDF_STATION, (*a_NumVals), a_PropId);
  }

  xfCloseGroup(XSectGroup);

  return status;
} /* xfGetStation */

/* --------------------------------------------------------------------------- */
/* FUNCTION  xfSetType */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfSetType (xid a_Id, int a_NumVals, int *a_PropId, 
                        int a_compression)
{
  xid XSectGroup;
  int status;

  /* open the XSect group */
  status = xfOpenGroup(a_Id, XSEC_XSECS_GROUP, &XSectGroup);
  if (status < 0) {
    return status;
  }
  status = xfWritePropertyInt(XSectGroup, XMDF_TYPE, a_NumVals, a_PropId, 
                              a_compression);
  xfCloseGroup(XSectGroup);

  return status;
} /* xfSetType */


/* --------------------------------------------------------------------------- */
/* FUNCTION  xfGetType */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfGetType (xid a_Id, int *a_NumVals, int *a_PropId)
{
  xid XSectGroup;
  int status;
  xmbool exists;

  /* open the XSect group */
  status = xfOpenGroup(a_Id, XSEC_XSECS_GROUP, &XSectGroup);
  if (status < 0) {
    return status;
  }
  xfDoesPropertyWithNameExist(XSectGroup, XMDF_TYPE, &exists);
  if (exists) {
    xfGetPropertyNumber(XSectGroup, XMDF_TYPE, a_NumVals);
    xfReadPropertyInt(XSectGroup, XMDF_TYPE, (*a_NumVals), a_PropId);
  }

  xfCloseGroup(XSectGroup);

  return status;
} /* xfGetType */

/* --------------------------------------------------------------------------- */
/* FUNCTION  xfSetpType */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfSetpType (xid a_Id, int a_NumVals, int *a_PropId, 
                         int a_compression)
{
  xid XSectGroup;
  int status;

  /* open the XSect group */
  status = xfOpenGroup(a_Id, XSEC_XSECS_GROUP, &XSectGroup);
  if (status < 0) {
    return status;
  }
  status = xfWritePropertyInt(XSectGroup, XMDF_PTYPE, a_NumVals, a_PropId, 
                              a_compression);
  xfCloseGroup(XSectGroup);

  return status;
} /* xfSetpType */


/* --------------------------------------------------------------------------- */
/* FUNCTION  xfGetpType */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfGetpType (xid a_Id, int *a_NumVals, int *a_PropId)
{
  xid XSectGroup;
  int status;
  xmbool exists;

  /* open the XSect group */
  status = xfOpenGroup(a_Id, XSEC_XSECS_GROUP, &XSectGroup);
  if (status < 0) {
    return status;
  }
  xfDoesPropertyWithNameExist(XSectGroup, XMDF_PTYPE, &exists);
  if (exists) {
    xfGetPropertyNumber(XSectGroup, XMDF_PTYPE, a_NumVals);
    xfReadPropertyInt(XSectGroup, XMDF_PTYPE, (*a_NumVals), a_PropId);
  }

  xfCloseGroup(XSectGroup);

  return status;
} /* xfGetpType */

/* --------------------------------------------------------------------------- */
/* FUNCTION  xfSetpCSDBLink */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfSetpCSDBLink (xid a_Id, int a_NumVals, int *a_PropId, 
                             int a_compression)
{
  xid XSectGroup;
  int status;

  /* open the XSect group */
  status = xfOpenGroup(a_Id, XSEC_XSECS_GROUP, &XSectGroup);
  if (status < 0) {
    return status;
  }
  status = xfWritePropertyInt(XSectGroup, XMDF_CSDBLINK, a_NumVals, a_PropId, 
                              a_compression);
  xfCloseGroup(XSectGroup);

  return status;
} /* xfSetpCSDBLink */


/* --------------------------------------------------------------------------- */
/* FUNCTION  xfGetpCSDBLink */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfGetpCSDBLink (xid a_Id, int *a_NumVals, int *a_PropId)
{
  xid XSectGroup;
  int status;
  xmbool exists;

  /* open the XSect group */
  status = xfOpenGroup(a_Id, XSEC_XSECS_GROUP, &XSectGroup);
  if (status < 0) {
    return status;
  }
  xfDoesPropertyWithNameExist(XSectGroup, XSEC_XSECS_GROUP, &exists);
  if (exists) {
    xfGetPropertyNumber(XSectGroup, XMDF_CSNAME, a_NumVals);
    xfReadPropertyInt(XSectGroup, XMDF_CSNAME, (*a_NumVals), a_PropId);
  }

  xfCloseGroup(XSectGroup);

  return status;
} /* xfGetpCSDBLink */

/* --------------------------------------------------------------------------- */
/* FUNCTION  xfSetNote */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfSetNote (xid a_Id, int a_NumVals, int a_StrLen, char *a_PropId)
{
  xid XSectGroup;
  int status;

  /* open the XSect group */
  status = xfOpenGroup(a_Id, XSEC_XSECS_GROUP, &XSectGroup);
  if (status < 0) {
    return status;
  }
  
  status = xfWritePropertyString(XSectGroup, XMDF_NOTE, a_NumVals, a_StrLen, 
                                 a_PropId);

  xfCloseGroup(XSectGroup);

  return status;
} /* xfSetNote */


/* --------------------------------------------------------------------------- */
/* FUNCTION  xfGetNote */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfGetNote (xid a_Id, int *a_NumVals, int *a_StrLen, 
                        char *a_PropId)
{
  xid XSectGroup;
  int status;
  xmbool exists;

  /* open the XSect group */
  status = xfOpenGroup(a_Id, XSEC_XSECS_GROUP, &XSectGroup);
  if (status < 0) {
    return status;
  }
  xfDoesPropertyWithNameExist(XSectGroup, XMDF_NOTE, &exists);
  if (exists) {
    xfGetPropertyNumber(XSectGroup, XMDF_NOTE, a_NumVals);
    xfGetPropertyStringLength(XSectGroup, XMDF_NOTE, a_NumVals, a_StrLen);
    xfReadPropertyString(XSectGroup, XMDF_NOTE, (*a_NumVals), (*a_StrLen), 
                         a_PropId);
  }

  xfCloseGroup(XSectGroup);

  return status;
} /* xfGetNote */

/* --------------------------------------------------------------------------- */
/* FUNCTION  xfGetNoteLen */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfGetNoteLen(xid a_Id, int *a_NumVals, int *a_StrLen) 
{
  xid XSectGroup;
  int status;
  xmbool exists;

  /* open the XSect group */
  status = xfOpenGroup(a_Id, XSEC_XSECS_GROUP, &XSectGroup);
  if (status < 0) {
    return status;
  }
  xfDoesPropertyWithNameExist(XSectGroup, XMDF_NOTE, &exists);
  if (exists) {
    xfGetPropertyNumber(XSectGroup, XMDF_NOTE, a_NumVals);
    xfGetPropertyStringLength(XSectGroup, XMDF_NOTE, a_NumVals, a_StrLen);
  }

  xfCloseGroup(XSectGroup);

  return status;
}


/* --------------------------------------------------------------------------- */
/* FUNCTION  xfSetXSectGeomX */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfSetXSectGeomX (xid a_Id, int a_index, int a_NumVals, 
                             double *a_iValues, int a_compression)
{
  xid XSecGeometryGroup, XSecGeometrySubGroup;
  int status;
  char grpname[30];

  /* open the XSect group */
  status = xfOpenGroup(a_Id, XSEC_GEOMETRY_GROUP, &XSecGeometryGroup);

  /* open the Geometry group */
  sprintf(grpname, "%s%d", XSEC_GEOMETRY_GROUP, a_index);
  status = xfOpenGroup(XSecGeometryGroup, grpname, &XSecGeometrySubGroup);
  if (status < 0) {
    return status;
  }
  status = xfWritePropertyDouble(XSecGeometrySubGroup, XMDF_GEOMETRY_X, 
                                   a_NumVals, a_iValues, a_compression);
  xfCloseGroup(XSecGeometrySubGroup);
  xfCloseGroup(XSecGeometryGroup);

  return status;
} /* xfSetXSectGeomX */


/* --------------------------------------------------------------------------- */
/* FUNCTION  xfGetXSectGeomX */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfGetXSectGeomX (xid a_Id, int a_index, int *a_NumVals, 
                             double *a_iValues) 
{
  xid XSecGeometryGroup, XSecGeometrySubGroup;
  int status;
  xmbool exists = XFALSE;
  char grpname[30];

  /* open the XSect group */
  status = xfOpenGroup(a_Id, XSEC_GEOMETRY_GROUP, &XSecGeometryGroup);

  /* open the Geometry group */
  sprintf(grpname, "%s%d", XSEC_GEOMETRY_GROUP, a_index);
  status = xfOpenGroup(XSecGeometryGroup, grpname, &XSecGeometrySubGroup);
  if (status < 0) {
    return status;
  }
  xfDoesPropertyWithNameExist(XSecGeometrySubGroup, XMDF_GEOMETRY_X, &exists);
  if (exists) {
    xfGetPropertyNumber(XSecGeometrySubGroup, XMDF_GEOMETRY_X, a_NumVals);
    xfReadPropertyDouble(XSecGeometrySubGroup, XMDF_GEOMETRY_X, (*a_NumVals), 
                         a_iValues);
  }

  xfCloseGroup(XSecGeometrySubGroup);
  xfCloseGroup(XSecGeometryGroup);

  return status;
}

/* --------------------------------------------------------------------------- */
/* FUNCTION  xfSetXSectGeomY */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfSetXSectGeomY (xid a_Id, int a_index, int a_NumVals, 
                             double *a_iValues, int a_compression)
{
  xid XSecGeometryGroup, XSecGeometrySubGroup;
  int status;
  char grpname[30];

  /* open the XSect group */
  status = xfOpenGroup(a_Id, XSEC_GEOMETRY_GROUP, &XSecGeometryGroup);

  /* open the Geometry group */
  sprintf(grpname, "%s%d", XSEC_GEOMETRY_GROUP, a_index);
  status = xfOpenGroup(XSecGeometryGroup, grpname, &XSecGeometrySubGroup);
  if (status < 0) {
    return status;
  }
  status = xfWritePropertyDouble(XSecGeometrySubGroup, XMDF_GEOMETRY_Y, 
                                   a_NumVals, a_iValues, a_compression);
  xfCloseGroup(XSecGeometrySubGroup);
  xfCloseGroup(XSecGeometryGroup);

  return status;
}


/* --------------------------------------------------------------------------- */
/* FUNCTION  xfGetXSectGeomY */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfGetXSectGeomY (xid a_Id, int a_index, int *a_NumVals, 
                             double *a_iValues)
{
  xid XSecGeometryGroup, XSecGeometrySubGroup;
  int status;
  xmbool  exists = XFALSE;
  char grpname[30];

  /* open the XSect group */
  status = xfOpenGroup(a_Id, XSEC_GEOMETRY_GROUP, &XSecGeometryGroup);

  /* open the Geometry group */
  sprintf(grpname, "%s%d", XSEC_GEOMETRY_GROUP, a_index);
  status = xfOpenGroup(XSecGeometryGroup, grpname, &XSecGeometrySubGroup);
  if (status < 0) {
    return status;
  }
  xfDoesPropertyWithNameExist(XSecGeometrySubGroup, XMDF_GEOMETRY_Y, &exists);
  if (exists) {
    xfGetPropertyNumber(XSecGeometrySubGroup, XMDF_GEOMETRY_Y, a_NumVals);
    xfReadPropertyDouble(XSecGeometrySubGroup, XMDF_GEOMETRY_Y, (*a_NumVals), 
                         a_iValues);
  }

  xfCloseGroup(XSecGeometrySubGroup);
  xfCloseGroup(XSecGeometryGroup);

  return status;
}


/* --------------------------------------------------------------------------- */
/* FUNCTION  xfSetXSectGeomD */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfSetXSectGeomD (xid a_Id, int a_index, int a_NumVals, 
                             double *a_iValues, int a_compression)
{
  xid XSecGeometryGroup, XSecGeometrySubGroup;
  int status;
  char grpname[30];

  /* open the XSect group */
  status = xfOpenGroup(a_Id, XSEC_GEOMETRY_GROUP, &XSecGeometryGroup);

  /* open the Geometry group */
  sprintf(grpname, "%s%d", XSEC_GEOMETRY_GROUP, a_index);
  status = xfOpenGroup(XSecGeometryGroup, grpname, &XSecGeometrySubGroup);
  if (status < 0) {
    return status;
  }
  status = xfWritePropertyDouble(XSecGeometrySubGroup, XMDF_GEOMETRY_D, 
                                   a_NumVals, a_iValues, a_compression);
  xfCloseGroup(XSecGeometrySubGroup);
  xfCloseGroup(XSecGeometryGroup);

  return status;
}


/* --------------------------------------------------------------------------- */
/* FUNCTION  xfGetXSectGeomD */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfGetXSectGeomD (xid a_Id, int a_index, int *a_NumVals, 
                             double *a_iValues)
{
  xid XSecGeometryGroup, XSecGeometrySubGroup;
  int status;
  xmbool  exists = XFALSE;
  char grpname[30];

  /* open the XSect group */
  status = xfOpenGroup(a_Id, XSEC_GEOMETRY_GROUP, &XSecGeometryGroup);

  /* open the Geometry group */
  sprintf(grpname, "%s%d", XSEC_GEOMETRY_GROUP, a_index);
  status = xfOpenGroup(XSecGeometryGroup, grpname, &XSecGeometrySubGroup);
  if (status < 0) {
    return status;
  }
  xfDoesPropertyWithNameExist(XSecGeometrySubGroup, XMDF_GEOMETRY_D, &exists);
  if (exists) {
    xfGetPropertyNumber(XSecGeometrySubGroup, XMDF_GEOMETRY_D, a_NumVals);
    xfReadPropertyDouble(XSecGeometrySubGroup, XMDF_GEOMETRY_D, (*a_NumVals), 
                         a_iValues);
  }

  xfCloseGroup(XSecGeometrySubGroup);
  xfCloseGroup(XSecGeometryGroup);

  return status;
}


/* --------------------------------------------------------------------------- */
/* FUNCTION  xfSetXSectGeomZ */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfSetXSectGeomZ (xid a_Id, int a_index, int a_NumVals, 
                             double *a_iValues, int a_compression)
{
  xid XSecGeometryGroup, XSecGeometrySubGroup;
  int status;
  char grpname[30];

  /* open the XSect group */
  status = xfOpenGroup(a_Id, XSEC_GEOMETRY_GROUP, &XSecGeometryGroup);

  /* open the Geometry group */
  sprintf(grpname, "%s%d", XSEC_GEOMETRY_GROUP, a_index);
  status = xfOpenGroup(XSecGeometryGroup, grpname, &XSecGeometrySubGroup);
  if (status < 0) {
    return status;
  }
  status = xfWritePropertyDouble(XSecGeometrySubGroup, XMDF_GEOMETRY_Z, 
                                   a_NumVals, a_iValues, a_compression);
  xfCloseGroup(XSecGeometrySubGroup);
  xfCloseGroup(XSecGeometryGroup);

  return status;
}


/* --------------------------------------------------------------------------- */
/* FUNCTION  xfGetXSectGeomZ */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfGetXSectGeomZ (xid a_Id, int a_index, int *a_NumVals, 
                             double *a_iValues)
{
  xid XSecGeometryGroup, XSecGeometrySubGroup;
  int status;
  xmbool  exists = XFALSE;
  char grpname[30];

  /* open the XSect group */
  status = xfOpenGroup(a_Id, XSEC_GEOMETRY_GROUP, &XSecGeometryGroup);

  /* open the Geometry group */
  sprintf(grpname, "%s%d", XSEC_GEOMETRY_GROUP, a_index);
  status = xfOpenGroup(XSecGeometryGroup, grpname, &XSecGeometrySubGroup);
  if (status < 0) {
    return status;
  }
  xfDoesPropertyWithNameExist(XSecGeometrySubGroup, XMDF_GEOMETRY_Z, &exists);
  if (exists) {
    xfGetPropertyNumber(XSecGeometrySubGroup, XMDF_GEOMETRY_Z, a_NumVals);
    xfReadPropertyDouble(XSecGeometrySubGroup, XMDF_GEOMETRY_Z, (*a_NumVals), 
                         a_iValues);
  }

  xfCloseGroup(XSecGeometrySubGroup);
  xfCloseGroup(XSecGeometryGroup);

  return status;
}


/*Line Properties */

/* --------------------------------------------------------------------------- */
/* FUNCTION  xfSetXSectLinePropFrom */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfSetXSectLinePropFrom (xid a_Id, int a_index, int a_NumVals, 
                                     double *a_PropId, int a_compression)
{
  xid XSecLinePropertyGroup, XSecLinePropertySubGroup;
  int status;
  char grpname[30];

  /* open the XSect group */
  status = xfOpenGroup(a_Id, XSEC_LINE_GROUP, &XSecLinePropertyGroup);

  /* open the Property group */
  sprintf(grpname, "%s%d", XSEC_LINE_GROUP, a_index);
  status = xfOpenGroup(XSecLinePropertyGroup, grpname, 
        &XSecLinePropertySubGroup);
  if (status < 0) {
    return status;
  }
  status = xfWritePropertyDouble(XSecLinePropertySubGroup, XMDF_LP_FROM, 
                                   a_NumVals, a_PropId, a_compression);
  xfCloseGroup(XSecLinePropertySubGroup);
  xfCloseGroup(XSecLinePropertyGroup);

  return status;
}


/* --------------------------------------------------------------------------- */
/* FUNCTION  xfGetXSectLinePropFrom */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfGetXSectLinePropFrom (xid a_Id, int a_index, int *a_NumVals, 
                                     double *a_PropId)
{
  xid XSecLinePropertyGroup, XSecLinePropertySubGroup;
  int status;
  xmbool  exists = XFALSE;
  char grpname[30];

  /* open the XSect group */
  status = xfOpenGroup(a_Id, XSEC_LINE_GROUP, &XSecLinePropertyGroup);

  /* open the Property group */
  sprintf(grpname, "%s%d", XSEC_LINE_GROUP, a_index);
  status = xfOpenGroup(XSecLinePropertyGroup, grpname, 
        &XSecLinePropertySubGroup);
  if (status < 0) {
    return status;
  }
  xfDoesPropertyWithNameExist(XSecLinePropertySubGroup, XMDF_LP_FROM, &exists);
  if (exists) {
    xfGetPropertyNumber(XSecLinePropertySubGroup, XMDF_LP_FROM, a_NumVals);
    xfReadPropertyDouble(XSecLinePropertySubGroup, XMDF_LP_FROM, (*a_NumVals), 
          a_PropId);
  }

  xfCloseGroup(XSecLinePropertySubGroup);
  xfCloseGroup(XSecLinePropertyGroup);

  return status;
}

/* --------------------------------------------------------------------------- */
/* FUNCTION  xfSetXSectLinePropTo */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfSetXSectLinePropTo (xid a_Id, int a_index, int a_NumVals, 
                                   double *a_PropId, int a_compression) 
{
  xid XSecLinePropertyGroup, XSecLinePropertySubGroup;
  int status;
  char grpname[30];

  /* open the XSect group */
  status = xfOpenGroup(a_Id, XSEC_LINE_GROUP, &XSecLinePropertyGroup);

  /* open the Line Property group */
  sprintf(grpname, "%s%d", XSEC_LINE_GROUP, a_index);
  status = xfOpenGroup(XSecLinePropertyGroup, grpname, 
        &XSecLinePropertySubGroup);
  if (status < 0) {
    return status;
  }
  status = xfWritePropertyDouble(XSecLinePropertySubGroup, XMDF_LP_TO, 
                                   a_NumVals, a_PropId, a_compression);
  xfCloseGroup(XSecLinePropertySubGroup);
  xfCloseGroup(XSecLinePropertyGroup);

  return status;
}

/* --------------------------------------------------------------------------- */
/* FUNCTION  xfGetXSectLinePropTo */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfGetXSectLinePropTo (xid a_Id, int a_index, int *a_NumVals, 
                                   double *a_PropId)
{
  xid XSecLinePropertyGroup, XSecLinePropertySubGroup;
  int status;
  xmbool  exists = XFALSE;
  char grpname[30];

  /* open the XSect group */
  status = xfOpenGroup(a_Id, XSEC_LINE_GROUP, &XSecLinePropertyGroup);

  /* open the Property group */
  sprintf(grpname, "%s%d", XSEC_LINE_GROUP, a_index);
  status = xfOpenGroup(XSecLinePropertyGroup, grpname, 
        &XSecLinePropertySubGroup);
  if (status < 0) {
    return status;
  }
  xfDoesPropertyWithNameExist(XSecLinePropertySubGroup, XMDF_LP_TO, &exists);
  if (exists) {
    xfGetPropertyNumber(XSecLinePropertySubGroup, XMDF_LP_TO, a_NumVals);
    xfReadPropertyDouble(XSecLinePropertySubGroup, XMDF_LP_TO, (*a_NumVals), 
                         a_PropId);
  }

  xfCloseGroup(XSecLinePropertySubGroup);
  xfCloseGroup(XSecLinePropertyGroup);

  return status;
}

/* --------------------------------------------------------------------------- */
/* FUNCTION  xfSetXSectLinePropType */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfSetXSectLinePropType (xid a_Id, int a_index, int a_NumVals, 
                                     int *a_PropId, int a_compression)
{
  xid XSecLinePropertyGroup, XSecLinePropertySubGroup;
  int status;
  char grpname[30];

  /* open the XSect group */
  status = xfOpenGroup(a_Id, XSEC_LINE_GROUP, &XSecLinePropertyGroup);

  /* open the Line Property group */
  sprintf(grpname, "%s%d", XSEC_LINE_GROUP, a_index);
  status = xfOpenGroup(XSecLinePropertyGroup, grpname, 
        &XSecLinePropertySubGroup);
  if (status < 0) {
    return status;
  }
  status = xfWritePropertyInt(XSecLinePropertySubGroup, XMDF_LP_TYPE, 
                                   a_NumVals, a_PropId, a_compression);
  xfCloseGroup(XSecLinePropertySubGroup);
  xfCloseGroup(XSecLinePropertyGroup);

  return status;
}

/* --------------------------------------------------------------------------- */
/* FUNCTION  xfGetXSectLinePropType */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfGetXSectLinePropType (xid a_Id, int a_index, int *a_NumVals, 
                                     int *a_PropId) 
{
  xid XSecLinePropertyGroup, XSecLinePropertySubGroup;
  int status;
  xmbool  exists = XFALSE;
  char grpname[30];

  /* open the XSect group */
  status = xfOpenGroup(a_Id, XSEC_LINE_GROUP, &XSecLinePropertyGroup);

  /* open the Property group */
  sprintf(grpname, "%s%d", XSEC_LINE_GROUP, a_index);
  status = xfOpenGroup(XSecLinePropertyGroup, grpname, 
        &XSecLinePropertySubGroup);
  if (status < 0) {
    return status;
  }
  xfDoesPropertyWithNameExist(XSecLinePropertySubGroup, XMDF_LP_TYPE, &exists);
  if (exists) {
    xfGetPropertyNumber(XSecLinePropertySubGroup, XMDF_LP_TYPE, a_NumVals);
    xfReadPropertyInt(XSecLinePropertySubGroup, XMDF_LP_TYPE, (*a_NumVals), 
                      a_PropId);
  }

  xfCloseGroup(XSecLinePropertySubGroup);
  xfCloseGroup(XSecLinePropertyGroup);

  return status;
}

/* --------------------------------------------------------------------------- */
/* FUNCTION  xfSetXSectLinePropIValue */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfSetXSectLinePropIValue (xid a_Id, int a_index, int a_NumVals, 
                                       int *a_PropId, int a_compression)
{
  xid XSecLinePropertyGroup, XSecLinePropertySubGroup;
  int status;
  char grpname[30];

  /* open the XSect group */
  status = xfOpenGroup(a_Id, XSEC_LINE_GROUP, &XSecLinePropertyGroup);

  /* open the Line Property group */
  sprintf(grpname, "%s%d", XSEC_LINE_GROUP, a_index);
  status = xfOpenGroup(XSecLinePropertyGroup, grpname, 
        &XSecLinePropertySubGroup);
  if (status < 0) {
    return status;
  }
  status = xfWritePropertyInt(XSecLinePropertySubGroup, XMDF_LP_IVALUE, 
                                   a_NumVals, a_PropId, a_compression);
  xfCloseGroup(XSecLinePropertySubGroup);
  xfCloseGroup(XSecLinePropertyGroup);

  return status;
}

/* --------------------------------------------------------------------------- */
/* FUNCTION  xfGetXSectLinePropIValue */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfGetXSectLinePropIValue (xid a_Id, int a_index, int *a_NumVals, 
                                       int *a_PropId)
{
  xid XSecLinePropertyGroup, XSecLinePropertySubGroup;
  int status;
  xmbool  exists = XFALSE;
  char grpname[30];

  /* open the XSect group */
  status = xfOpenGroup(a_Id, XSEC_LINE_GROUP, &XSecLinePropertyGroup);

  /* open the Property group */
  sprintf(grpname, "%s%d", XSEC_LINE_GROUP, a_index);
  status = xfOpenGroup(XSecLinePropertyGroup, grpname, 
        &XSecLinePropertySubGroup);
  if (status < 0) {
    return status;
  }
  xfDoesPropertyWithNameExist(XSecLinePropertySubGroup, XMDF_LP_IVALUE, 
        &exists);
  if (exists) {
    xfGetPropertyNumber(XSecLinePropertySubGroup, XMDF_LP_IVALUE, a_NumVals);
    xfReadPropertyInt(XSecLinePropertySubGroup, XMDF_LP_IVALUE, (*a_NumVals), 
                      a_PropId);
  }

  xfCloseGroup(XSecLinePropertySubGroup);
  xfCloseGroup(XSecLinePropertyGroup);

  return status;
}

/* --------------------------------------------------------------------------- */
/* FUNCTION  xfSetXSectLinePropFValue */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfSetXSectLinePropFValue (xid a_Id, int a_index, int a_NumVals, 
                                       double *a_PropId, int a_compression)
{
  xid XSecLinePropertyGroup, XSecLinePropertySubGroup;
  int status;
  char grpname[30];

  /* open the XSect group */
  status = xfOpenGroup(a_Id, XSEC_LINE_GROUP, &XSecLinePropertyGroup);

  /* open the Line Property group */
  sprintf(grpname, "%s%d", XSEC_LINE_GROUP, a_index);
  status = xfOpenGroup(XSecLinePropertyGroup, grpname, 
        &XSecLinePropertySubGroup);
  if (status < 0) {
    return status;
  }
  status = xfWritePropertyDouble(XSecLinePropertySubGroup, XMDF_LP_FVALUE, 
                                   a_NumVals, a_PropId, a_compression);
  xfCloseGroup(XSecLinePropertySubGroup);
  xfCloseGroup(XSecLinePropertyGroup);

  return status;
}

/* --------------------------------------------------------------------------- */
/* FUNCTION  xfGetXSectLinePropFValue */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfGetXSectLinePropFValue (xid a_Id, int a_index, int *a_NumVals, 
                                       double *a_PropId)
{
  xid XSecLinePropertyGroup, XSecLinePropertySubGroup;
  int status;
  xmbool  exists = XFALSE;
  char grpname[30];

  /* open the XSect group */
  status = xfOpenGroup(a_Id, XSEC_LINE_GROUP, &XSecLinePropertyGroup);

  /* open the Property group */
  sprintf(grpname, "%s%d", XSEC_LINE_GROUP, a_index);
  status = xfOpenGroup(XSecLinePropertyGroup, grpname, 
        &XSecLinePropertySubGroup);
  if (status < 0) {
    return status;
  }
  xfDoesPropertyWithNameExist(XSecLinePropertySubGroup, XMDF_LP_FVALUE, 
        &exists);
  if (exists) {
    xfGetPropertyNumber(XSecLinePropertySubGroup, XMDF_LP_FVALUE, a_NumVals);
    xfReadPropertyDouble(XSecLinePropertySubGroup, XMDF_LP_FVALUE, 
          (*a_NumVals), a_PropId);
  }

  xfCloseGroup(XSecLinePropertySubGroup);
  xfCloseGroup(XSecLinePropertyGroup);

  return status;
}

/*wmsejj */
/*line properties--2 */
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfSetXSectLinePropID */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfSetXSectLinePropID (xid a_Id, int a_NumVals, 
                                   int *a_PropId, int a_compression)
{
  xid XSecLinePropertyGroup;
  int status;

  /* open the XSect group */
  status = xfOpenGroup(a_Id, XSEC_LINE_GROUP, &XSecLinePropertyGroup);

  if (status < 0) {
    return status;
  }
  status = xfWritePropertyInt(XSecLinePropertyGroup, XMDF_LNID, 
                                   a_NumVals, a_PropId, a_compression);
  xfCloseGroup(XSecLinePropertyGroup);

  return status;
}
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfGetXSectLinePropID */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfGetXSectLinePropID (xid a_Id, int *a_NumVals, 
                                   int *a_PropId)
{
  xid XSecLinePropertyGroup;
  int status;
  xmbool exists = XFALSE;

  /* open the XSect group */
  status = xfOpenGroup(a_Id, XSEC_LINE_GROUP, &XSecLinePropertyGroup);

  if (status < 0) {
    return status;
  }
  xfDoesPropertyWithNameExist(XSecLinePropertyGroup, XMDF_LNID, &exists);
  if (exists) {
    xfGetPropertyNumber(XSecLinePropertyGroup, XMDF_LNID, a_NumVals);
    xfReadPropertyInt(XSecLinePropertyGroup, XMDF_LNID, (*a_NumVals), 
          a_PropId);
  }

  xfCloseGroup(XSecLinePropertyGroup);

  return status;
}
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfSetXSectLinePropName */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfSetXSectLinePropName (xid a_Id, int a_NumVals, int a_StrLen, 
                                     char *a_PropId)
{
  xid XSecLinePropertyGroup;
  int status;

  /* open the XSect group */
  status = xfOpenGroup(a_Id, XSEC_LINE_GROUP, &XSecLinePropertyGroup);

  if (status < 0) {
    return status;
  }

  status = xfWritePropertyString(XSecLinePropertyGroup, XMDF_LN_NAME, 
          a_NumVals, a_StrLen, a_PropId);

  xfCloseGroup(XSecLinePropertyGroup);

  return status;
}
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfGetXSectLinePropName */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfGetXSectLinePropName (xid a_Id, int *a_NumVals, int *a_StrLen, 
                                     char *a_PropId)
{
  xid XSecLinePropertyGroup;
  int status;
  xmbool exists = XFALSE;

  /* open the XSect group */
  status = xfOpenGroup(a_Id, XSEC_LINE_GROUP, &XSecLinePropertyGroup);

  if (status < 0) {
    return status;
  }
  xfDoesPropertyWithNameExist(XSecLinePropertyGroup, XMDF_LN_NAME, &exists);
  if (exists) {
    xfGetPropertyNumber(XSecLinePropertyGroup, XMDF_LN_NAME, a_NumVals);
    xfGetPropertyStringLength(XSecLinePropertyGroup, XMDF_LN_NAME, a_NumVals, 
                              a_StrLen);
    xfReadPropertyString(XSecLinePropertyGroup, XMDF_LN_NAME, (*a_NumVals), 
                        (*a_StrLen), a_PropId);
  }
  xfCloseGroup(XSecLinePropertyGroup);

  return status;
}

/* --------------------------------------------------------------------------- */
/* FUNCTION  xfGetXSectLinePropNameLen */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfGetXSectLinePropNameLen (xid a_Id, int *a_NumVals, 
                                        int *a_StrLen)
{
  xid XSecLinePropertyGroup;
  int status;
  xmbool exists = XFALSE;

  /* open the XSect group */
  status = xfOpenGroup(a_Id, XSEC_LINE_GROUP, &XSecLinePropertyGroup);

  if (status < 0) {
    return status;
  }
  xfDoesPropertyWithNameExist(XSecLinePropertyGroup, XMDF_LN_NAME, &exists);
  if (exists) {
    xfGetPropertyNumber(XSecLinePropertyGroup, XMDF_LN_NAME, a_NumVals);
    xfGetPropertyStringLength(XSecLinePropertyGroup, XMDF_LN_NAME, a_NumVals,
                              a_StrLen);
  }
  xfCloseGroup(XSecLinePropertyGroup);

  return status;
}

/* --------------------------------------------------------------------------- */
/* FUNCTION  xfSetXSectLinePropDesc */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfSetXSectLinePropDesc (xid a_Id, int a_NumVals, int a_StrLen, 
                                     char *a_PropId)
{
  xid XSecLinePropertyGroup;
  int status;

  /* open the XSect group */
  status = xfOpenGroup(a_Id, XSEC_LINE_GROUP, &XSecLinePropertyGroup);

  if (status < 0) {
    return status;
  }

  status = xfWritePropertyString(XSecLinePropertyGroup, XMDF_LN_DESC, 
                                 a_NumVals, a_StrLen, a_PropId);

  xfCloseGroup(XSecLinePropertyGroup);

  return status;
}
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfGetXSectLinePropDesc */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfGetXSectLinePropDesc (xid a_Id, int *a_NumVals, int *a_StrLen, 
                                     char *a_PropId)
{
  xid XSecLinePropertyGroup;
  int status;
  xmbool exists = XFALSE;

  /* open the XSect group */
  status = xfOpenGroup(a_Id, XSEC_LINE_GROUP, &XSecLinePropertyGroup);

  if (status < 0) {
    return status;
  }
  xfDoesPropertyWithNameExist(XSecLinePropertyGroup, XMDF_LN_DESC, &exists);
  if (exists) {
    xfGetPropertyNumber(XSecLinePropertyGroup, XMDF_LN_DESC, a_NumVals);
    xfGetPropertyStringLength(XSecLinePropertyGroup, XMDF_LN_DESC, a_NumVals, 
                              a_StrLen);
    xfReadPropertyString(XSecLinePropertyGroup, XMDF_LN_DESC, (*a_NumVals), 
                        (*a_StrLen), a_PropId);
  }

  xfCloseGroup(XSecLinePropertyGroup);

  return status;
}

/* --------------------------------------------------------------------------- */
/* FUNCTION  xfGetXSectLinePropDescLen */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfGetXSectLinePropDescLen (xid a_Id, int *a_NumVals, 
                                        int *a_StrLen)
{
  xid XSecLinePropertyGroup;
  int status;
  xmbool exists = XFALSE;

  /* open the XSect group */
  status = xfOpenGroup(a_Id, XSEC_LINE_GROUP, &XSecLinePropertyGroup);

  if (status < 0) {
    return status;
  }
  xfDoesPropertyWithNameExist(XSecLinePropertyGroup, XMDF_LN_DESC, &exists);
  if (exists) {
    xfGetPropertyNumber(XSecLinePropertyGroup, XMDF_LN_DESC, a_NumVals);
    xfGetPropertyStringLength(XSecLinePropertyGroup, XMDF_LN_DESC, a_NumVals,
                              a_StrLen);
  }

  xfCloseGroup(XSecLinePropertyGroup);

  return status;
}

/* --------------------------------------------------------------------------- */
/* FUNCTION  xfSetXSectLinePropCategory */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfSetXSectLinePropCategory (xid a_Id, int a_NumVals, 
                                         int *a_PropId, int a_compression)
{
  xid XSecLinePropertyGroup;
  int status;

  /* open the XSect group */
  status = xfOpenGroup(a_Id, XSEC_LINE_GROUP, &XSecLinePropertyGroup);

  if (status < 0) {
    return status;
  }
  status = xfWritePropertyInt(XSecLinePropertyGroup, XMDF_LN_CATEGORY, 
                                   a_NumVals, a_PropId, a_compression);
  xfCloseGroup(XSecLinePropertyGroup);

  return status;
}
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfGetXSectLinePropCategory */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfGetXSectLinePropCategory (xid a_Id, int *a_NumVals, 
                                         int *a_PropId)
{
  xid XSecLinePropertyGroup;
  int status;
  xmbool exists = XFALSE;

  /* open the XSect group */
  status = xfOpenGroup(a_Id, XSEC_LINE_GROUP, &XSecLinePropertyGroup);

  if (status < 0) {
    return status;
  }
  xfDoesPropertyWithNameExist(XSecLinePropertyGroup, XMDF_LN_CATEGORY, 
        &exists);
  if (exists) {
    xfGetPropertyNumber(XSecLinePropertyGroup, XMDF_LN_CATEGORY, a_NumVals);
    xfReadPropertyInt(XSecLinePropertyGroup, XMDF_LN_CATEGORY, (*a_NumVals), 
                      a_PropId);
  }

  xfCloseGroup(XSecLinePropertyGroup);

  return status;
}
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfSetXSectLinePropFDefault */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfSetXSectLinePropFDefault (xid a_Id, int a_NumVals, 
                                         double *a_PropId, int a_compression)
{
  xid XSecLinePropertyGroup;
  int status;

  /* open the XSect group */
  status = xfOpenGroup(a_Id, XSEC_LINE_GROUP, &XSecLinePropertyGroup);

  if (status < 0) {
    return status;
  }
  status = xfWritePropertyDouble(XSecLinePropertyGroup, XMDF_LN_FDEFAULT, 
                                   a_NumVals, a_PropId, a_compression);
  xfCloseGroup(XSecLinePropertyGroup);

  return status;
}
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfGetXSectLinePropFDefault */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfGetXSectLinePropFDefault (xid a_Id, int *a_NumVals, 
                                         double *a_PropId)
{
  xid XSecLinePropertyGroup;
  int status;
  xmbool exists = XFALSE;

  /* open the XSect group */
  status = xfOpenGroup(a_Id, XSEC_LINE_GROUP, &XSecLinePropertyGroup);

  if (status < 0) {
    return status;
  }
  xfDoesPropertyWithNameExist(XSecLinePropertyGroup, XMDF_LN_FDEFAULT, 
        &exists);
  if (exists) {
    xfGetPropertyNumber(XSecLinePropertyGroup, XMDF_LN_FDEFAULT, a_NumVals);
    xfReadPropertyDouble(XSecLinePropertyGroup, XMDF_LN_FDEFAULT, (*a_NumVals),
                         a_PropId);
  }

  xfCloseGroup(XSecLinePropertyGroup);

  return status;
}
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfSetXSectLinePropExclusive */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfSetXSectLinePropExclusive (xid a_Id, int a_NumVals, 
                                          int *a_PropId, int a_compression)
{
  xid XSecLinePropertyGroup;
  int status;

  /* open the XSect group */
  status = xfOpenGroup(a_Id, XSEC_LINE_GROUP, &XSecLinePropertyGroup);

  if (status < 0) {
    return status;
  }
  status = xfWritePropertyInt(XSecLinePropertyGroup, XMDF_LN_EXCLUSIVE, 
                                   a_NumVals, a_PropId, a_compression);
  xfCloseGroup(XSecLinePropertyGroup);

  return status;
}
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfGetXSectLinePropExclusive */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfGetXSectLinePropExclusive (xid a_Id, int *a_NumVals, 
                                          int *a_PropId)
{
  xid XSecLinePropertyGroup;
  int status;
  xmbool exists = XFALSE;

  /* open the XSect group */
  status = xfOpenGroup(a_Id, XSEC_LINE_GROUP, &XSecLinePropertyGroup);

  if (status < 0) {
    return status;
  }
  xfDoesPropertyWithNameExist(XSecLinePropertyGroup, XMDF_LN_EXCLUSIVE, 
        &exists);
  if (exists) {
    xfGetPropertyNumber(XSecLinePropertyGroup, XMDF_LN_EXCLUSIVE, a_NumVals);
    xfReadPropertyInt(XSecLinePropertyGroup, XMDF_LN_EXCLUSIVE, (*a_NumVals), 
                      a_PropId);
  }

  xfCloseGroup(XSecLinePropertyGroup);

  return status;
}


/* --------------------------------------------------------------------------- */
/* FUNCTION  xfSetNumberOfLinePropEnumGroup */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfSetNumberOfLinePropEnumGroup (xid a_Id, int *a_nPropNum, 
                                        int a_compression)
{
  xid XSecLinePropertyGroup, XSecLinePropertySubGroup;
  int status, i;
  char grpname[30];

  /* open the XSect group */
  status = xfOpenGroup(a_Id, XSEC_LINE_GROUP, &XSecLinePropertyGroup);
  if (status < 0) {
    return status;
  }

  status = xfWritePropertyInt(XSecLinePropertyGroup, XMDF_XSECT_PROPNUM, 
                                1, a_nPropNum, a_compression);

  /* open the Line Property group */
  for (i=0; i<(*a_nPropNum); i++) {
    sprintf(grpname, "%s%d", XMDF_XSECT_PROPNUM2, i);
    XSecLinePropertySubGroup = xfpCreateGroup(XSecLinePropertyGroup, 
                                   grpname, GROUP_TYPE_XSEC_LINE);
    xfCloseGroup(XSecLinePropertySubGroup);
  }
  xfCloseGroup(XSecLinePropertyGroup);
  return status;
}
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfGetNumberOfLinePropEnumGroup */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfGetNumberOfLinePropEnumGroup (xid a_Id, int *a_nXSects)
{
  xid XSecLinePropertyGroup;
  int status;
  xmbool exists;

  /* open the XSect group */
  status = xfOpenGroup(a_Id, XSEC_LINE_GROUP, &XSecLinePropertyGroup);
  if (status < 0) {
    return status;
  }
  xfDoesPropertyWithNameExist(XSecLinePropertyGroup, XMDF_XSECT_PROPNUM, 
        &exists);
  if (exists) {
    xfReadPropertyInt(XSecLinePropertyGroup, XMDF_XSECT_PROPNUM, 1, a_nXSects);
  }

  xfCloseGroup(XSecLinePropertyGroup);

  return status;
}


/* --------------------------------------------------------------------------- */
/* FUNCTION  xfSetNumberOfLinePropEnum */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfSetNumberOfLinePropEnum (xid a_Id, int a_index, 
                                        int *a_nPropNum, int a_compression)
{
  xid XSecLinePropertyGroup, XSecLinePropertySubGroup;
  int status;
  char grpname[30];

  /* open the XSect group */
  status = xfOpenGroup(a_Id, XSEC_LINE_GROUP, &XSecLinePropertyGroup);
  if (status < 0) {
    return status;
  }

  /* open the Line Property group */
  sprintf(grpname, "%s%d", XMDF_XSECT_PROPNUM2, a_index);
  status = xfOpenGroup(XSecLinePropertyGroup, grpname, 
                      &XSecLinePropertySubGroup);
  if (status < 0) {
    return status;
  }
  status = xfWritePropertyInt(XSecLinePropertySubGroup, XMDF_XSECT_PROPENUMNUM,
                                1, a_nPropNum, a_compression);

  xfCloseGroup(XSecLinePropertySubGroup);
  xfCloseGroup(XSecLinePropertyGroup);
  return status;
}
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfGetNumberOfLinePropEnum */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfGetNumberOfLinePropEnum (xid a_Id, int a_index, int *a_nXSects)
{
  xid XSecLinePropertyGroup, XSecLinePropertySubGroup;
  int status;
  xmbool exists;
  char grpname[30];

  /* open the XSect group */
  status = xfOpenGroup(a_Id, XSEC_LINE_GROUP, &XSecLinePropertyGroup);
  if (status < 0) {
    return status;
  }
  /* open the Line Property group */
  sprintf(grpname, "%s%d", XMDF_XSECT_PROPNUM2, a_index);
  status = xfOpenGroup(XSecLinePropertyGroup, grpname, 
                      &XSecLinePropertySubGroup);
  if (status < 0) {
    return status;
  }
  xfDoesPropertyWithNameExist(XSecLinePropertySubGroup, XMDF_XSECT_PROPENUMNUM,
        &exists);
  if (exists) {
    xfReadPropertyInt(XSecLinePropertySubGroup, XMDF_XSECT_PROPENUMNUM, 1, 
          a_nXSects);
  }

  xfCloseGroup(XSecLinePropertyGroup);

  return status;
}

/* --------------------------------------------------------------------------- */
/* FUNCTION  xfSetXSectLinePropEnumID */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfSetXSectLinePropEnumID (xid a_Id, int a_index, int a_NumVals, 
                                       int *a_PropId, int a_compression)
{
  xid XSecLinePropertyGroup, XSecLinePropertySubGroup;
  int status;
  char grpname[30];

  /* open the XSect group */
  status = xfOpenGroup(a_Id, XSEC_LINE_GROUP, &XSecLinePropertyGroup);
  if (status < 0) {
    return status;
  }
  /* open the Line Property group */
  sprintf(grpname, "%s%d", XMDF_XSECT_PROPNUM2, a_index);
  status = xfOpenGroup(XSecLinePropertyGroup, grpname, 
                      &XSecLinePropertySubGroup);
  if (status < 0) {
    return status;
  }
  status = xfWritePropertyInt(XSecLinePropertySubGroup, XMDF_LN_PROPID, 
                                   a_NumVals, a_PropId, a_compression);
  xfCloseGroup(XSecLinePropertySubGroup);
  xfCloseGroup(XSecLinePropertyGroup);
  return status;
}
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfGetXSectLinePropEnumID */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfGetXSectLinePropEnumID (xid a_Id, int a_index, int *a_NumVals, 
                                       int *a_PropId)
{
  xid XSecLinePropertyGroup, XSecLinePropertySubGroup;
  int status;
  xmbool exists;
  char grpname[30];

  /* open the XSect group */
  status = xfOpenGroup(a_Id, XSEC_LINE_GROUP, &XSecLinePropertyGroup);
  if (status < 0) {
    return status;
  }
  /* open the Line Property group */
  sprintf(grpname, "%s%d", XMDF_XSECT_PROPNUM2, a_index);
  status = xfOpenGroup(XSecLinePropertyGroup, grpname, 
        &XSecLinePropertySubGroup);
  
  xfDoesPropertyWithNameExist(XSecLinePropertySubGroup, XMDF_LN_PROPID, 
        &exists);
  if (exists) {
    xfGetPropertyNumber(XSecLinePropertySubGroup, XMDF_LN_PROPID, a_NumVals);
    xfReadPropertyInt(XSecLinePropertySubGroup, XMDF_LN_PROPID, (*a_NumVals), 
                      a_PropId);
  }
  xfCloseGroup(XSecLinePropertySubGroup);
  xfCloseGroup(XSecLinePropertyGroup);

  return status;
}
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfSetXSectLinePropEnumMatID */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfSetXSectLinePropEnumMatID (xid a_Id, int a_index, int a_NumVals,
                                          int *a_PropId, int a_compression)
{
  xid XSecLinePropertyGroup, XSecLinePropertySubGroup;
  int status;
  char grpname[30];

  /* open the XSect group */
  status = xfOpenGroup(a_Id, XSEC_LINE_GROUP, &XSecLinePropertyGroup);
  if (status < 0) {
    return status;
  }
  /* open the Line Property group */
  sprintf(grpname, "%s%d", XMDF_XSECT_PROPNUM2, a_index);
  status = xfOpenGroup(XSecLinePropertyGroup, grpname, 
                      &XSecLinePropertySubGroup);
  if (status < 0) {
    return status;
  }
  status = xfWritePropertyInt(XSecLinePropertySubGroup, XMDF_LN_PROPENUMID, 
                                   a_NumVals, a_PropId, a_compression);
  xfCloseGroup(XSecLinePropertySubGroup);
  xfCloseGroup(XSecLinePropertyGroup);
  return status;
}
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfGetXSectLinePropEnumMatID */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfGetXSectLinePropEnumMatID (xid a_Id, int a_index, 
                                          int *a_NumVals, int *a_PropId)
{
  xid XSecLinePropertyGroup, XSecLinePropertySubGroup;
  int status;
  xmbool exists;
  char grpname[30];

  /* open the XSect group */
  status = xfOpenGroup(a_Id, XSEC_LINE_GROUP, &XSecLinePropertyGroup);
  if (status < 0) {
    return status;
  }
  /* open the Line Property group */
  sprintf(grpname, "%s%d", XMDF_XSECT_PROPNUM2, a_index);
  status = xfOpenGroup(XSecLinePropertyGroup, grpname, 
                      &XSecLinePropertySubGroup);
  
  xfDoesPropertyWithNameExist(XSecLinePropertySubGroup, XMDF_LN_PROPENUMID,
                              &exists);
  if (exists) {
    xfGetPropertyNumber(XSecLinePropertySubGroup, XMDF_LN_PROPENUMID, 
                        a_NumVals);
    xfReadPropertyInt(XSecLinePropertySubGroup, XMDF_LN_PROPENUMID, 
                    (*a_NumVals), a_PropId);
  }
  xfCloseGroup(XSecLinePropertySubGroup);
  xfCloseGroup(XSecLinePropertyGroup);

  return status;
}

/* --------------------------------------------------------------------------- */
/* FUNCTION  xfSetXSectLinePropEnumName */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfSetXSectLinePropEnumName (xid a_Id, int a_index, int a_NumVals, 
                                         int a_StrLen, char *a_PropId)
{
  xid XSecLinePropertyGroup, XSecLinePropertySubGroup;
  int status;
  char grpname[30];

  /* open the XSect group */
  status = xfOpenGroup(a_Id, XSEC_LINE_GROUP, &XSecLinePropertyGroup);
  if (status < 0) {
    return status;
  }
  /* open the Line Property group */
  sprintf(grpname, "%s%d", XMDF_XSECT_PROPNUM2, a_index);
  status = xfOpenGroup(XSecLinePropertyGroup, grpname, 
                      &XSecLinePropertySubGroup);
  if (status < 0) {
    return status;
  }

  status = xfWritePropertyString(XSecLinePropertySubGroup, XMDF_LN_PROPNAME,  
                                 a_NumVals, a_StrLen, a_PropId);

  xfCloseGroup(XSecLinePropertySubGroup);
  xfCloseGroup(XSecLinePropertyGroup);
  return status;
}
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfGetXSectLinePropEnumName */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfGetXSectLinePropEnumName (xid a_Id, int a_index, int *a_NumVals,
                                         int *a_StrLen, char *a_PropId)
{
  xid XSecLinePropertyGroup, XSecLinePropertySubGroup;
  int status;
  xmbool exists;
  char grpname[30];

  /* open the XSect group */
  status = xfOpenGroup(a_Id, XSEC_LINE_GROUP, &XSecLinePropertyGroup);
  if (status < 0) {
    return status;
  }
  /* open the Line Property group */
  sprintf(grpname, "%s%d", XMDF_XSECT_PROPNUM2, a_index);
  status = xfOpenGroup(XSecLinePropertyGroup, grpname, 
                      &XSecLinePropertySubGroup);
  if (status < 0) {
    return status;
  }
  
  xfDoesPropertyWithNameExist(XSecLinePropertySubGroup, XMDF_LN_PROPNAME, 
                             &exists);
  if (exists) {
    xfGetPropertyNumber(XSecLinePropertySubGroup, XMDF_LN_PROPNAME, a_NumVals);
    xfGetPropertyStringLength(XSecLinePropertySubGroup, XMDF_LN_PROPNAME, 
                              a_NumVals, a_StrLen);
    xfReadPropertyString(XSecLinePropertySubGroup, XMDF_LN_PROPNAME, 
                        (*a_NumVals), (*a_StrLen), a_PropId);
  }
  xfCloseGroup(XSecLinePropertySubGroup);
  xfCloseGroup(XSecLinePropertyGroup);

  return status;
}

/* --------------------------------------------------------------------------- */
/* FUNCTION  xfGetXSectLinePropEnumNameLen */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfGetXSectLinePropEnumNameLen (xid a_Id, int a_index, 
                                            int *a_NumVals, int *a_StrLen)
{
  xid XSecLinePropertyGroup, XSecLinePropertySubGroup;
  int status;
  xmbool exists;
  char grpname[30];

  /* open the XSect group */
  status = xfOpenGroup(a_Id, XSEC_LINE_GROUP, &XSecLinePropertyGroup);
  if (status < 0) {
    return status;
  }
  /* open the Line Property group */
  sprintf(grpname, "%s%d", XMDF_XSECT_PROPNUM2, a_index);
  status = xfOpenGroup(XSecLinePropertyGroup, grpname, 
                      &XSecLinePropertySubGroup);
  if (status < 0) {
    return status;
  }
  
  xfDoesPropertyWithNameExist(XSecLinePropertySubGroup, XMDF_LN_PROPNAME, 
                              &exists);
  if (exists) {
    xfGetPropertyNumber(XSecLinePropertySubGroup, XMDF_LN_PROPNAME, a_NumVals);
    xfGetPropertyStringLength(XSecLinePropertySubGroup, XMDF_LN_PROPNAME, 
                              a_NumVals, a_StrLen);
  }
  xfCloseGroup(XSecLinePropertySubGroup);
  xfCloseGroup(XSecLinePropertyGroup);

  return status;
}

/*wmsejj */
/*Point Properties */
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfSetXSectPointPropMeasure */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfSetXSectPointPropMeasure (xid a_Id, int a_index, int a_NumVals, 
                                         double *a_PropId, int a_compression)
{
  xid XSecPointPropertyGroup, XSecPointPropertySubGroup;
  int status;
  char grpname[30];

  /* open the XSect group */
  status = xfOpenGroup(a_Id, XSEC_POINT_GROUP, &XSecPointPropertyGroup);
  if (status < 0) {
    return status;
  }
  /* open the Point Property group */
  sprintf(grpname, "%s%d", XSEC_POINT_GROUP, a_index);
  status = xfOpenGroup(XSecPointPropertyGroup, grpname, 
        &XSecPointPropertySubGroup);
  if (status < 0) {
    return status;
  }

  status = xfWritePropertyDouble(XSecPointPropertySubGroup, XMDF_PP_MEASURE, 
                                   a_NumVals, a_PropId, a_compression);
  xfCloseGroup(XSecPointPropertySubGroup);
  xfCloseGroup(XSecPointPropertyGroup);

  return status;
}

/* --------------------------------------------------------------------------- */
/* FUNCTION  xfGetXSectPointPropMeasure */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfGetXSectPointPropMeasure (xid a_Id, int a_index, int *a_NumVals,
                                         double *a_PropId)
{
  xid XSecPointPropertyGroup, XSecPointPropertySubGroup;
  int status;
  xmbool exists = XFALSE;
  char grpname[30];

  /* open the XSect group */
  status = xfOpenGroup(a_Id, XSEC_POINT_GROUP, &XSecPointPropertyGroup);
  if (status < 0) {
    return status;
  }
  /* open the Point Property group */
  sprintf(grpname, "%s%d", XSEC_POINT_GROUP, a_index);
  status = xfOpenGroup(XSecPointPropertyGroup, grpname, 
        &XSecPointPropertySubGroup);
  if (status < 0) {
    return status;
  }

  xfDoesPropertyWithNameExist(XSecPointPropertySubGroup, XMDF_PP_MEASURE, 
        &exists);
  if (exists) {
    xfGetPropertyNumber(XSecPointPropertySubGroup, XMDF_PP_MEASURE, a_NumVals);
    xfReadPropertyDouble(XSecPointPropertySubGroup, XMDF_PP_MEASURE, 
            (*a_NumVals), a_PropId);
  }

  xfCloseGroup(XSecPointPropertySubGroup);
  xfCloseGroup(XSecPointPropertyGroup);

  return status;
}

/* --------------------------------------------------------------------------- */
/* FUNCTION  xfSetXSectPointPropType */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfSetXSectPointPropType (xid a_Id, int a_index, int a_NumVals, 
                                      int *a_PropId, 
                                      int a_compression)
{
  xid XSecPointPropertyGroup, XSecPointPropertySubGroup;
  int status;
  char grpname[30];

  /* open the XSect group */
  status = xfOpenGroup(a_Id, XSEC_POINT_GROUP, &XSecPointPropertyGroup);

  if (status < 0) {
    return status;
  }
  /* open the Point Property group */
  sprintf(grpname, "%s%d", XSEC_POINT_GROUP, a_index);
  status = xfOpenGroup(XSecPointPropertyGroup, grpname, 
                      &XSecPointPropertySubGroup);
  if (status < 0) {
    return status;
  }

  status = xfWritePropertyInt(XSecPointPropertySubGroup, XMDF_PP_TYPE, 
                                   a_NumVals, a_PropId, a_compression);
  xfCloseGroup(XSecPointPropertySubGroup);
  xfCloseGroup(XSecPointPropertyGroup);

  return status;
}
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfGetXSectPointPropType */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfGetXSectPointPropType (xid a_Id, int a_index, int *a_NumVals, 
                                      int *a_PropId)
{
  xid XSecPointPropertyGroup, XSecPointPropertySubGroup;
  int status;
  xmbool exists = XFALSE;
  char grpname[30];

  /* open the XSect group */
  status = xfOpenGroup(a_Id, XSEC_POINT_GROUP, &XSecPointPropertyGroup);

  if (status < 0) {
    return status;
  }
  /* open the Point Property group */
  sprintf(grpname, "%s%d", XSEC_POINT_GROUP, a_index);
  status = xfOpenGroup(XSecPointPropertyGroup, grpname, 
                      &XSecPointPropertySubGroup);
  if (status < 0) {
    return status;
  }

  xfDoesPropertyWithNameExist(XSecPointPropertySubGroup, XMDF_PP_TYPE, 
                              &exists);
  if (exists) {
    xfGetPropertyNumber(XSecPointPropertySubGroup, XMDF_PP_TYPE, a_NumVals);
    xfReadPropertyInt(XSecPointPropertySubGroup, XMDF_PP_TYPE, (*a_NumVals), 
                      a_PropId);
  }

  xfCloseGroup(XSecPointPropertySubGroup);
  xfCloseGroup(XSecPointPropertyGroup);

  return status;
}
/*wmsejj */

/*XSec--Point Properties--2 */
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfSetXSectPointPropID */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfSetXSectPointPropID (xid a_Id, int a_NumVals, 
                                          int *a_PropId, int a_compression)
{
  xid XSecPointPropGroup;
  int status;

  /* open the XSect group */
  status = xfOpenGroup(a_Id, XSEC_POINT_GROUP, &XSecPointPropGroup);

  status = xfWritePropertyInt(XSecPointPropGroup, XMDF_PTID, 
                              a_NumVals, a_PropId, a_compression);
  xfCloseGroup(XSecPointPropGroup);

  return status;
}
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfGetXSectPointPropID */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfGetXSectPointPropID (xid a_Id, int *a_NumVals, 
                                          int *a_PropId)
{
  xid XSecPointPropGroup;
  int status;
  xmbool exists = XFALSE;

  /* open the XSect group */
  status = xfOpenGroup(a_Id, XSEC_POINT_GROUP, &XSecPointPropGroup);

  xfDoesPropertyWithNameExist(XSecPointPropGroup, XMDF_PTID, &exists);
  if (exists) {
    xfGetPropertyNumber(XSecPointPropGroup, XMDF_PTID, a_NumVals);
    xfReadPropertyInt(XSecPointPropGroup, XMDF_PTID, (*a_NumVals), a_PropId);
  }

  xfCloseGroup(XSecPointPropGroup);

  return status;
}
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfSetXSectPointPropName */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfSetXSectPointPropName (xid a_Id, int a_NumVals, int a_StrLen, 
                                          char *a_PropId)
{
  xid XSecPointPropGroup;
  int status;

  /* open the XSect group */
  status = xfOpenGroup(a_Id, XSEC_POINT_GROUP, &XSecPointPropGroup);

  status = xfWritePropertyString(XSecPointPropGroup, XMDF_PP_NAME, a_NumVals, 
                                 a_StrLen, a_PropId);
  
  xfCloseGroup(XSecPointPropGroup);

  return status;
}

/* --------------------------------------------------------------------------- */
/* FUNCTION  xfGetXSectPointPropName */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfGetXSectPointPropName (xid a_Id, int *a_NumVals, int *a_StrLen, 
                                          char *a_PropId)
{
  xid XSecPointPropGroup;
  int status;
  xmbool exists = XFALSE;

  /* open the XSect group */
  status = xfOpenGroup(a_Id, XSEC_POINT_GROUP, &XSecPointPropGroup);

  xfDoesPropertyWithNameExist(XSecPointPropGroup, XMDF_PP_NAME, &exists);
  if (exists) {
    xfGetPropertyNumber(XSecPointPropGroup, XMDF_PP_NAME, a_NumVals);
    xfGetPropertyStringLength(XSecPointPropGroup, XMDF_PP_NAME, a_NumVals, 
          a_StrLen);
    xfReadPropertyString(XSecPointPropGroup, XMDF_PP_NAME, (*a_NumVals), 
          (*a_StrLen), a_PropId);
  }

  xfCloseGroup(XSecPointPropGroup);

  return status;
}

/* --------------------------------------------------------------------------- */
/* FUNCTION  xfGetXSectPointPropNameLen */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfGetXSectPointPropNameLen (xid a_Id, int *a_NumVals, 
                                         int *a_StrLen)
{
  xid XSecPointPropGroup;
  int status;
  xmbool exists = XFALSE;

  /* open the XSect group */
  status = xfOpenGroup(a_Id, XSEC_POINT_GROUP, &XSecPointPropGroup);

  xfDoesPropertyWithNameExist(XSecPointPropGroup, XMDF_PP_NAME, &exists);
  if (exists) {
    xfGetPropertyNumber(XSecPointPropGroup, XMDF_PP_NAME, a_NumVals);
    xfGetPropertyStringLength(XSecPointPropGroup, XMDF_PP_NAME, a_NumVals, 
          a_StrLen);
  }

  xfCloseGroup(XSecPointPropGroup);

  return status;
}

/* --------------------------------------------------------------------------- */
/* FUNCTION  xfSetXSectPointPropDesc */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfSetXSectPointPropDesc (xid a_Id, int a_NumVals, int a_StrLen, 
                                          char *a_PropId)
{
  xid XSecPointPropGroup;
  int status;

  /* open the XSect group */
  status = xfOpenGroup(a_Id, XSEC_POINT_GROUP, &XSecPointPropGroup);

  status = xfWritePropertyString(XSecPointPropGroup, XMDF_PP_DESC, a_NumVals, 
                                 a_StrLen, a_PropId);

  xfCloseGroup(XSecPointPropGroup);

  return status;
}

/* --------------------------------------------------------------------------- */
/* FUNCTION  xfGetXSectPointPropDesc */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfGetXSectPointPropDesc (xid a_Id, int *a_NumVals, int *a_StrLen, 
                                          char *a_PropId)
{
  xid XSecPointPropGroup;
  int status;
  xmbool exists = XFALSE;

  /* open the XSect group */
  status = xfOpenGroup(a_Id, XSEC_POINT_GROUP, &XSecPointPropGroup);

  xfDoesPropertyWithNameExist(XSecPointPropGroup, XMDF_PP_DESC, &exists);
  if (exists) {
    xfGetPropertyNumber(XSecPointPropGroup, XMDF_PP_DESC, a_NumVals);
    xfGetPropertyStringLength(XSecPointPropGroup, XMDF_NOTE, a_NumVals, 
          a_StrLen);
    xfReadPropertyString(XSecPointPropGroup, XMDF_PP_DESC, (*a_NumVals), 
          (*a_StrLen), a_PropId);
  }

  xfCloseGroup(XSecPointPropGroup);

  return status;
}

/* --------------------------------------------------------------------------- */
/* FUNCTION  xfGetXSectPointPropDescLen */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfGetXSectPointPropDescLen (xid a_Id, int *a_NumVals, 
                                         int *a_StrLen)
{
  xid XSecPointPropGroup;
  int status;
  xmbool exists = XFALSE;

  /* open the XSect group */
  status = xfOpenGroup(a_Id, XSEC_POINT_GROUP, &XSecPointPropGroup);

  xfDoesPropertyWithNameExist(XSecPointPropGroup, XMDF_PP_DESC, &exists);
  if (exists) {
    xfGetPropertyNumber(XSecPointPropGroup, XMDF_PP_DESC, a_NumVals);
    xfGetPropertyStringLength(XSecPointPropGroup, XMDF_NOTE, a_NumVals, 
          a_StrLen);
  }

  xfCloseGroup(XSecPointPropGroup);

  return status;
}

/* --------------------------------------------------------------------------- */
/* FUNCTION  xfSetXSectPointPropExclusive */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfSetXSectPointPropExclusive (xid a_Id, int a_NumVals, 
                                           int *a_PropId, int a_compression)
{
  xid XSecPointPropGroup;
  int status;

  /* open the XSect group */
  status = xfOpenGroup(a_Id, XSEC_POINT_GROUP, &XSecPointPropGroup);

  status = xfWritePropertyInt(XSecPointPropGroup, XMDF_PP_EXCLUSIVE, 
                                   a_NumVals, a_PropId, a_compression);
  xfCloseGroup(XSecPointPropGroup);

  return status;
}

/* --------------------------------------------------------------------------- */
/* FUNCTION  xfGetXSectPointPropExclusive */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfGetXSectPointPropExclusive (xid a_Id, int *a_NumVals, 
                                           int *a_PropId)
{
  xid XSecPointPropGroup;
  int status;

  /* open the XSect group */
  status = xfOpenGroup(a_Id, XSEC_POINT_GROUP, &XSecPointPropGroup);

  xfCloseGroup(XSecPointPropGroup);

  return status;
}

/* --------------------------------------------------------------------------- */
/* FUNCTION  xfSetXSectTopoName */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfSetXSectTopoName (xid a_Id, int a_NumVals, int a_StrLen, 
                                 char *a_PropId)
{
  xid XSecProfileGroup;
  int status;

  /* open the XSect group */
  status = xfOpenGroup(a_Id, XSEC_PROFILES_GROUP, &XSecProfileGroup);

  if (status < 0) {
    return status;
  }

  status = xfWritePropertyString(XSecProfileGroup, XMDF_TOPO_NAME, a_NumVals, 
                                 a_StrLen, a_PropId);
  xfCloseGroup(XSecProfileGroup);

  return status;
}

/* --------------------------------------------------------------------------- */
/* FUNCTION  xfGetXSectTopoName */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfGetXSectTopoName (xid a_Id, int *a_NumVals, int *a_StrLen, 
                                 char *a_PropId)
{
  xid XSecProfileGroup;
  int status;
  xmbool exists = XFALSE;

  /* open the XSect group */
  status = xfOpenGroup(a_Id, XSEC_PROFILES_GROUP, &XSecProfileGroup);

  if (status < 0) {
    return status;
  }
  xfDoesPropertyWithNameExist(XSecProfileGroup, XMDF_TOPO_NAME, &exists);
  if (exists) {
    xfGetPropertyNumber(XSecProfileGroup, XMDF_TOPO_NAME, a_NumVals);
    xfGetPropertyStringLength(XSecProfileGroup, XMDF_TOPO_NAME, a_NumVals, a_StrLen);
    xfReadPropertyString(XSecProfileGroup, XMDF_TOPO_NAME, (*a_NumVals), (*a_StrLen),
                         a_PropId);
  }

  xfCloseGroup(XSecProfileGroup);

  return status;
}

/* --------------------------------------------------------------------------- */
/* FUNCTION  xfGetXSectTopoNameLen */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfGetXSectTopoNameLen (xid a_Id, int *a_NumVals, int *a_StrLen)
{
  xid XSecProfileGroup;
  int status;
  xmbool exists = XFALSE;

  /* open the XSect group */
  status = xfOpenGroup(a_Id, XSEC_PROFILES_GROUP, &XSecProfileGroup);

  if (status < 0) {
    return status;
  }
  xfDoesPropertyWithNameExist(XSecProfileGroup, XMDF_TOPO_NAME, &exists);
  if (exists) {
    xfGetPropertyNumber(XSecProfileGroup, XMDF_TOPO_NAME, a_NumVals);
    xfGetPropertyStringLength(XSecProfileGroup, XMDF_TOPO_NAME, a_NumVals, 
          a_StrLen);
  }

  xfCloseGroup(XSecProfileGroup);

  return status;
}

/* --------------------------------------------------------------------------- */
/* FUNCTION  xfSetXSectTopoDesc */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfSetXSectTopoDesc (xid a_Id, int a_NumVals, int a_StrLen, 
                                 char *a_PropId)
{
  xid XSecProfileGroup;
  int status;

  /* open the XSect group */
  status = xfOpenGroup(a_Id, XSEC_PROFILES_GROUP, &XSecProfileGroup);

  if (status < 0) {
    return status;
  }

  status = xfWritePropertyString(XSecProfileGroup, XMDF_TOPO_DESC, a_NumVals, 
                                 a_StrLen, a_PropId);

  xfCloseGroup(XSecProfileGroup);

  return status;
}
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfGetXSectTopoDesc */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfGetXSectTopoDesc (xid a_Id, int *a_NumVals, int *a_StrLen, 
                                 char *a_PropId)
{
  xid XSecProfileGroup;
  int status;
  xmbool exists = XFALSE;

  /* open the XSect group */
  status = xfOpenGroup(a_Id, XSEC_PROFILES_GROUP, &XSecProfileGroup);

  if (status < 0) {
    return status;
  }
  xfDoesPropertyWithNameExist(XSecProfileGroup, XMDF_TOPO_DESC, &exists);
  if (exists) {
    xfGetPropertyNumber(XSecProfileGroup, XMDF_TOPO_DESC, a_NumVals);
    xfGetPropertyStringLength(XSecProfileGroup, XMDF_TOPO_DESC, a_NumVals,
          a_StrLen);
    xfReadPropertyString(XSecProfileGroup, XMDF_TOPO_DESC, (*a_NumVals), 
          (*a_StrLen), a_PropId);
  }

  xfCloseGroup(XSecProfileGroup);

  return status;
}

/* --------------------------------------------------------------------------- */
/* FUNCTION  xfGetXSectTopoDescLen */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfGetXSectTopoDescLen (xid a_Id, int *a_NumVals, int *a_StrLen)
{
  xid XSecProfileGroup;
  int status;
  xmbool exists = XFALSE;

  /* open the XSect group */
  status = xfOpenGroup(a_Id, XSEC_PROFILES_GROUP, &XSecProfileGroup);

  if (status < 0) {
    return status;
  }
  xfDoesPropertyWithNameExist(XSecProfileGroup, XMDF_TOPO_DESC, &exists);
  if (exists) {
    xfGetPropertyNumber(XSecProfileGroup, XMDF_TOPO_DESC, a_NumVals);
    xfGetPropertyStringLength(XSecProfileGroup, XMDF_TOPO_DESC, a_NumVals, 
          a_StrLen);
  }

  xfCloseGroup(XSecProfileGroup);

  return status;
}



/*wmsejj */

/* --------------------------------------------------------------------------- */
/* FUNCTION  xfGetMetadataXML */
/*! PURPOSE:   Gets a single XML string containing all the metadata for the
 *           specified object's ID in the XMDF file

 * - NOTES:     Be sure to free the XML string returned from this function.  */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfGetMetadataXML (xid a_Id, char **a_xml)
{
  xid MetadataGroup;
  int status, NumTopics, NumComments, StrLenComments, StrLenTopics,
    year, month, day, hour, minute, second, i, NumProperties, PropertyLen;
  size_t len;
  xmbool exists;
  char *abstract = NULL, *comments = NULL, *topics = NULL, *ProfileXML = NULL,
    *purpose = NULL, *SourceXML = NULL, *SpatialXML = NULL, *title = NULL,
    datestring[11];
  double startdate = 0.0;
  unsigned char era;

  /* open the Metadata group */
  status = xfOpenGroup(a_Id, METADATA_MAIN, &MetadataGroup);
  if (status < 0) {
    return status;
  }
  /* Read Abstract */
  xfAllocateReadPropertyString(MetadataGroup, METADATA_MAIN_ABSTRACT,
    &NumProperties, &PropertyLen, &abstract);
  /* Read Generic Comments */
  xfDoesPropertyWithNameExist(MetadataGroup, METADATA_MAIN_COMMENTS, &exists);
  if (exists) {
    xfGetPropertyStringLength(MetadataGroup, METADATA_MAIN_COMMENTS,
      &NumComments, &StrLenComments);
    if (NumComments > 0 && StrLenComments > 0) {
      comments = (char *)malloc(NumComments*(StrLenComments+1)*sizeof(char));
      xfReadPropertyString(MetadataGroup, METADATA_MAIN_COMMENTS, NumComments,
        StrLenComments, comments);
    }
  }
  /* Read Generic Topics */
  xfDoesPropertyWithNameExist(MetadataGroup, METADATA_MAIN_TOPICS, &exists);
  if (exists) {
    xfGetPropertyStringLength(MetadataGroup, METADATA_MAIN_TOPICS, &NumTopics,
      &StrLenTopics);
    if (NumTopics > 0 && StrLenTopics > 0) {
      topics = (char *)malloc(NumTopics*(StrLenTopics+1)*sizeof(char));
      xfReadPropertyString(MetadataGroup, METADATA_MAIN_TOPICS, NumTopics,
        StrLenTopics, topics);
    }
  }

  /* Read Metadata Profile Information */
  xfGetMetadataProfileXML(MetadataGroup, &ProfileXML);
  /* Read Purpose */
  xfAllocateReadPropertyString(MetadataGroup, METADATA_MAIN_PURPOSE,
    &NumProperties, &PropertyLen, &purpose);
  /* Read Metadata Source Information */
  xfGetMetadataSourceXML(MetadataGroup, &SourceXML);
  /* Read Metadata Spatial Information */
  xfGetMetadataSpatialXML(MetadataGroup, &SpatialXML);
  /* Read Start Date */
  xfReadPropertyDouble(MetadataGroup, METADATA_MAIN_STARTDATE, 1, &startdate);
  /* Read Title */
  xfAllocateReadPropertyString(MetadataGroup, METADATA_MAIN_TITLE,
    &NumProperties, &PropertyLen, &title);

  xfCloseGroup(MetadataGroup);

  /* Create the XML string */
  len = 0;

  len += strlen(METATAG_MAIN);
  *a_xml = (char *)malloc((len+1)*sizeof(char));
  strcpy(*a_xml, METATAG_MAIN);

  len += strlen(METATAG_MAIN_TITLE);
  *a_xml = (char *)realloc(*a_xml, (len+1)*sizeof(char));
  strcat(*a_xml, METATAG_MAIN_TITLE);

  if (title) {
    len += strlen(title);
    *a_xml = (char *)realloc(*a_xml, (len+1)*sizeof(char));
    strcat(*a_xml, title);
  }

  len += strlen(METATAG_MAIN_TITLE_END);
  *a_xml = (char *)realloc(*a_xml, (len+1)*sizeof(char));
  strcat(*a_xml, METATAG_MAIN_TITLE_END);

  len += strlen(METATAG_MAIN_ABSTRACT);
  *a_xml = (char *)realloc(*a_xml, (len+1)*sizeof(char));
  strcat(*a_xml, METATAG_MAIN_ABSTRACT);

  if (abstract) {
    len += strlen(abstract);
    *a_xml = (char *)realloc(*a_xml, (len+1)*sizeof(char));
    strcat(*a_xml, abstract);
  }

  len += strlen(METATAG_MAIN_ABSTRACT_END);
  *a_xml = (char *)realloc(*a_xml, (len+1)*sizeof(char));
  strcat(*a_xml, METATAG_MAIN_ABSTRACT_END);

  len += strlen(METATAG_MAIN_PURPOSE);
  *a_xml = (char *)realloc(*a_xml, (len+1)*sizeof(char));
  strcat(*a_xml, METATAG_MAIN_PURPOSE);

  if (purpose) {
    len += strlen(purpose);
    *a_xml = (char *)realloc(*a_xml, (len+1)*sizeof(char));
    strcat(*a_xml, purpose);
  }

  len += strlen(METATAG_MAIN_PURPOSE_END);
  *a_xml = (char *)realloc(*a_xml, (len+1)*sizeof(char));
  strcat(*a_xml, METATAG_MAIN_PURPOSE_END);

  len += strlen(METATAG_MAIN_STARTDATE);
  *a_xml = (char *)realloc(*a_xml, (len+1)*sizeof(char));
  strcat(*a_xml, METATAG_MAIN_STARTDATE);

  xfJulianToCalendar(&era, &year, &month, &day, &hour, &minute, &second,
    startdate);
  sprintf(datestring, "%2d/%2d/%4d", month, day, year);
  len += strlen(datestring);
  *a_xml = (char *)realloc(*a_xml, (len+1)*sizeof(char));
  strcat(*a_xml, datestring);

  len += strlen(METATAG_MAIN_STARTDATE_END);
  *a_xml = (char *)realloc(*a_xml, (len+1)*sizeof(char));
  strcat(*a_xml, METATAG_MAIN_STARTDATE_END);

  if (topics) {
    for (i=0; i<NumTopics; i++) {
      len += strlen(METATAG_MAIN_TOPIC);
      *a_xml = (char *)realloc(*a_xml, (len+1)*sizeof(char));
      strcat(*a_xml, METATAG_MAIN_TOPIC);

      len += strlen(&topics[i*StrLenTopics]);
      *a_xml = (char *)realloc(*a_xml, (len+1)*sizeof(char));
      strcat(*a_xml, &topics[i*StrLenTopics]);

      len += strlen(METATAG_MAIN_TOPIC_END);
      *a_xml = (char *)realloc(*a_xml, (len+1)*sizeof(char));
      strcat(*a_xml, METATAG_MAIN_TOPIC_END);
    }
  }

  if (comments) {
    for (i=0; i<NumComments; i++) {
      len += strlen(METATAG_MAIN_COMMENT);
      *a_xml = (char *)realloc(*a_xml, (len+1)*sizeof(char));
      strcat(*a_xml, METATAG_MAIN_COMMENT);

      len += strlen(&comments[i*StrLenComments]);
      *a_xml = (char *)realloc(*a_xml, (len+1)*sizeof(char));
      strcat(*a_xml, &comments[i*StrLenComments]);

      len += strlen(METATAG_MAIN_COMMENT_END);
      *a_xml = (char *)realloc(*a_xml, (len+1)*sizeof(char));
      strcat(*a_xml, METATAG_MAIN_COMMENT_END);
    }
  }

  if (ProfileXML) {
    len += strlen(ProfileXML);
    *a_xml = (char *)realloc(*a_xml, (len+1)*sizeof(char));
    strcat(*a_xml, ProfileXML);
  }

  if (SourceXML) {
    len += strlen(SourceXML);
    *a_xml = (char *)realloc(*a_xml, (len+1)*sizeof(char));
    strcat(*a_xml, SourceXML);
  }

  if (SpatialXML) {
    len += strlen(SpatialXML);
    *a_xml = (char *)realloc(*a_xml, (len+1)*sizeof(char));
    strcat(*a_xml, SpatialXML);
  }

  len += strlen(METATAG_MAIN_END);
  *a_xml = (char *)realloc(*a_xml, (len+1)*sizeof(char));
  strcat(*a_xml, METATAG_MAIN_END);

  /* Clean up all the allocated memory */
  if (abstract) {
    free(abstract);
  }
  if (comments) {
    free(comments);
  }
  if (topics) {
    free(topics);
  }
  if (purpose) {
    free(purpose);
  }
  if (title) {
    free(title);
  }
  if (ProfileXML) {
    free(ProfileXML);
  }
  if (SourceXML) {
    free(SourceXML);
  }
  if (SpatialXML) {
    free(SpatialXML);
  }

  return status;
} /* xfGetMetadataXML */

/* --------------------------------------------------------------------------- */
/* FUNCTION  xfGetMetadataProfileXML */
/*! PURPOSE:   Gets an XML string containing the profile associated with the
 *           metadata in the XMDF file

 * - NOTES:     Be sure to free the XML string returned from this function.  */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfGetMetadataProfileXML (xid a_Id, char **a_xml)
{
  xid ProfileGroup;
  int status, NumProperties, PropertyLen;
  size_t len;
  char *address = NULL, *addresstype = NULL, *contactperson = NULL,
    *contactposition = NULL, *email = NULL, *fax = NULL, *organization = NULL,
    *phone = NULL;

  /* open the Profile group */
  status = xfOpenGroup(a_Id, METADATA_MAIN_PROFILE, &ProfileGroup);
  if (status < 0) {
    return status;
  }
  /* Read Address */
  xfAllocateReadPropertyString(ProfileGroup, METADATA_PROFILE_ADDRESS,
    &NumProperties, &PropertyLen, &address);
  /* Read Address Type */
  xfAllocateReadPropertyString(ProfileGroup, METADATA_PROFILE_ADDRESSTYPE,
    &NumProperties, &PropertyLen, &addresstype);
  /* Read Contact Person */
  xfAllocateReadPropertyString(ProfileGroup, METADATA_PROFILE_CONTACTPERSON,
    &NumProperties, &PropertyLen, &contactperson);
  /* Read Contact Position */
  xfAllocateReadPropertyString(ProfileGroup, METADATA_PROFILE_CONTACTPOSITION,
    &NumProperties, &PropertyLen, &contactposition);
  /* Read Email */
  xfAllocateReadPropertyString(ProfileGroup, METADATA_PROFILE_EMAIL,
    &NumProperties, &PropertyLen, &email);
  /* Read Fax */
  xfAllocateReadPropertyString(ProfileGroup, METADATA_PROFILE_FAX,
    &NumProperties, &PropertyLen, &fax);
  /* Read Organization */
  xfAllocateReadPropertyString(ProfileGroup, METADATA_PROFILE_ORGANIZATION,
    &NumProperties, &PropertyLen, &organization);
  /* Read Phone */
  xfAllocateReadPropertyString(ProfileGroup, METADATA_PROFILE_PHONE,
    &NumProperties, &PropertyLen, &phone);

  xfCloseGroup(ProfileGroup);

  /* Create the XML string */
  len = 0;

  len += strlen(METATAG_MAIN_PROFILE);
  *a_xml = (char *)malloc((len+1)*sizeof(char));
  strcpy(*a_xml, METATAG_MAIN_PROFILE);

  len += strlen(METATAG_PROFILE_ADDRESS);
  *a_xml = (char *)realloc(*a_xml, (len+1)*sizeof(char));
  strcat(*a_xml, METATAG_PROFILE_ADDRESS);

  if (address) {
    len += strlen(address);
    *a_xml = (char *)realloc(*a_xml, (len+1)*sizeof(char));
    strcat(*a_xml, address);
  }

  len += strlen(METATAG_PROFILE_ADDRESS_END);
  *a_xml = (char *)realloc(*a_xml, (len+1)*sizeof(char));
  strcat(*a_xml, METATAG_PROFILE_ADDRESS_END);

  len += strlen(METATAG_PROFILE_ADDRESSTYPE);
  *a_xml = (char *)realloc(*a_xml, (len+1)*sizeof(char));
  strcat(*a_xml, METATAG_PROFILE_ADDRESSTYPE);

  if (addresstype) {
    len += strlen(addresstype);
    *a_xml = (char *)realloc(*a_xml, (len+1)*sizeof(char));
    strcat(*a_xml, addresstype);
  }

  len += strlen(METATAG_PROFILE_ADDRESSTYPE_END);
  *a_xml = (char *)realloc(*a_xml, (len+1)*sizeof(char));
  strcat(*a_xml, METATAG_PROFILE_ADDRESSTYPE_END);

  len += strlen(METATAG_PROFILE_CONTACTPERSON);
  *a_xml = (char *)realloc(*a_xml, (len+1)*sizeof(char));
  strcat(*a_xml, METATAG_PROFILE_CONTACTPERSON);

  if (contactperson) {
    len += strlen(contactperson);
    *a_xml = (char *)realloc(*a_xml, (len+1)*sizeof(char));
    strcat(*a_xml, contactperson);
  }

  len += strlen(METATAG_PROFILE_CONTACTPERSON_END);
  *a_xml = (char *)realloc(*a_xml, (len+1)*sizeof(char));
  strcat(*a_xml, METATAG_PROFILE_CONTACTPERSON_END);

  len += strlen(METATAG_PROFILE_CONTACTPOSITION);
  *a_xml = (char *)realloc(*a_xml, (len+1)*sizeof(char));
  strcat(*a_xml, METATAG_PROFILE_CONTACTPOSITION);

  if (contactposition) {
    len += strlen(contactposition);
    *a_xml = (char *)realloc(*a_xml, (len+1)*sizeof(char));
    strcat(*a_xml, contactposition);
  }

  len += strlen(METATAG_PROFILE_CONTACTPOSITION_END);
  *a_xml = (char *)realloc(*a_xml, (len+1)*sizeof(char));
  strcat(*a_xml, METATAG_PROFILE_CONTACTPOSITION_END);

  len += strlen(METATAG_PROFILE_EMAIL);
  *a_xml = (char *)realloc(*a_xml, (len+1)*sizeof(char));
  strcat(*a_xml, METATAG_PROFILE_EMAIL);

  if (email) {
    len += strlen(email);
    *a_xml = (char *)realloc(*a_xml, (len+1)*sizeof(char));
    strcat(*a_xml, email);
  }

  len += strlen(METATAG_PROFILE_EMAIL_END);
  *a_xml = (char *)realloc(*a_xml, (len+1)*sizeof(char));
  strcat(*a_xml, METATAG_PROFILE_EMAIL_END);

  len += strlen(METATAG_PROFILE_FAX);
  *a_xml = (char *)realloc(*a_xml, (len+1)*sizeof(char));
  strcat(*a_xml, METATAG_PROFILE_FAX);

  if (fax) {
    len += strlen(fax);
    *a_xml = (char *)realloc(*a_xml, (len+1)*sizeof(char));
    strcat(*a_xml, fax);
  }

  len += strlen(METATAG_PROFILE_FAX_END);
  *a_xml = (char *)realloc(*a_xml, (len+1)*sizeof(char));
  strcat(*a_xml, METATAG_PROFILE_FAX_END);

  len += strlen(METATAG_PROFILE_ORGANIZATION);
  *a_xml = (char *)realloc(*a_xml, (len+1)*sizeof(char));
  strcat(*a_xml, METATAG_PROFILE_ORGANIZATION);

  if (organization) {
    len += strlen(organization);
    *a_xml = (char *)realloc(*a_xml, (len+1)*sizeof(char));
    strcat(*a_xml, organization);
  }

  len += strlen(METATAG_PROFILE_ORGANIZATION_END);
  *a_xml = (char *)realloc(*a_xml, (len+1)*sizeof(char));
  strcat(*a_xml, METATAG_PROFILE_ORGANIZATION_END);

  len += strlen(METATAG_PROFILE_PHONE);
  *a_xml = (char *)realloc(*a_xml, (len+1)*sizeof(char));
  strcat(*a_xml, METATAG_PROFILE_PHONE);

  if (phone) {
    len += strlen(phone);
    *a_xml = (char *)realloc(*a_xml, (len+1)*sizeof(char));
    strcat(*a_xml, phone);
  }

  len += strlen(METATAG_PROFILE_PHONE_END);
  *a_xml = (char *)realloc(*a_xml, (len+1)*sizeof(char));
  strcat(*a_xml, METATAG_PROFILE_PHONE_END);

  len += strlen(METATAG_MAIN_PROFILE_END);
  *a_xml = (char *)realloc(*a_xml, (len+1)*sizeof(char));
  strcat(*a_xml, METATAG_MAIN_PROFILE_END);

  /* Clean up all the allocated memory */
  if (address) {
    free(address);
  }
  if (addresstype) {
    free(addresstype);
  }
  if (contactperson) {
    free(contactperson);
  }
  if (contactposition) {
    free(contactposition);
  }
  if (email) {
    free(email);
  }
  if (fax) {
    free(fax);
  }
  if (organization) {
    free(organization);
  }
  if (phone) {
    free(phone);
  }

  return status;
} /* xfGetMetadataProfileXML */

/* --------------------------------------------------------------------------- */
/* FUNCTION  xfGetMetadataSourceXML */
/*! PURPOSE:   Gets an XML string containing the source data associated with the
 *           metadata in the XMDF file

 * - NOTES:     Be sure to free the XML string returned from this function.  */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfGetMetadataSourceXML (xid a_Id, char **a_xml)
{
  xid SourceGroup;
  int status, NumXMSGenerated = 0, StrLenXMSGenerated = 0, i,
    NumProperties, PropertyLen;
  size_t len;
  char *userdefined = NULL, *XMSGenerated = NULL;
  xmbool exists;

  /* open the Profile group */
  status = xfOpenGroup(a_Id, METADATA_MAIN_SOURCE, &SourceGroup);
  if (status < 0) {
    return status;
  }
  /* Read Generic Topics */
  xfDoesPropertyWithNameExist(SourceGroup, METADATA_SOURCE_XMSGENERATED,
    &exists);
  if (exists) {
    xfGetPropertyStringLength(SourceGroup, METADATA_SOURCE_XMSGENERATED,
      &NumXMSGenerated, &StrLenXMSGenerated);
    if (NumXMSGenerated > 0 && StrLenXMSGenerated > 0) {
      XMSGenerated = (char *)malloc(NumXMSGenerated*(StrLenXMSGenerated+1)*
        sizeof(char));
      xfReadPropertyString(SourceGroup, METADATA_SOURCE_XMSGENERATED,
        NumXMSGenerated, StrLenXMSGenerated, XMSGenerated);
    }
  }
  /* Read User Defined */
  xfAllocateReadPropertyString(SourceGroup, METADATA_SOURCE_USERDEFINED,
    &NumProperties, &PropertyLen, &userdefined);

  xfCloseGroup(SourceGroup);

  /* Create the XML string */
  len = 0;

  len += strlen(METATAG_MAIN_SOURCE);
  *a_xml = (char *)malloc((len+1)*sizeof(char));
  strcpy(*a_xml, METATAG_MAIN_SOURCE);

  if (XMSGenerated) {
    for (i=0; i<NumXMSGenerated; i++) {
      len += strlen(METATAG_SOURCE_XMSGENERATED);
      *a_xml = (char *)realloc(*a_xml, (len+1)*sizeof(char));
      strcat(*a_xml, METATAG_SOURCE_XMSGENERATED);

      len += strlen(&XMSGenerated[i*StrLenXMSGenerated]);
      *a_xml = (char *)realloc(*a_xml, (len+1)*sizeof(char));
      strcat(*a_xml, &XMSGenerated[i*StrLenXMSGenerated]);

      len += strlen(METATAG_SOURCE_XMSGENERATED_END);
      *a_xml = (char *)realloc(*a_xml, (len+1)*sizeof(char));
      strcat(*a_xml, METATAG_SOURCE_XMSGENERATED_END);
    }
  }

  len += strlen(METATAG_SOURCE_USERDEFINED);
  *a_xml = (char *)realloc(*a_xml, (len+1)*sizeof(char));
  strcat(*a_xml, METATAG_SOURCE_USERDEFINED);

  if (userdefined) {
    len += strlen(userdefined);
    *a_xml = (char *)realloc(*a_xml, (len+1)*sizeof(char));
    strcat(*a_xml, userdefined);
  }

  len += strlen(METATAG_SOURCE_USERDEFINED_END);
  *a_xml = (char *)realloc(*a_xml, (len+1)*sizeof(char));
  strcat(*a_xml, METATAG_SOURCE_USERDEFINED_END);

  len += strlen(METATAG_MAIN_SOURCE_END);
  *a_xml = (char *)realloc(*a_xml, (len+1)*sizeof(char));
  strcat(*a_xml, METATAG_MAIN_SOURCE_END);

  /* Clean up all the allocated memory */
  if (XMSGenerated) {
    free(XMSGenerated);
  }
  if (userdefined) {
    free(userdefined);
  }

  return status;
} /* xfGetMetadataSourceXML */

/* ------------------------------------------------------------------------- */
/* FUNCTION  xfGetMetadataSpatialXML */
/*! PURPOSE:   Gets an XML string containing the spatial data associated with the
 *           metadata in the XMDF file

 * - NOTES:     Be sure to free the XML string returned from this function.  */
/* ------------------------------------------------------------------------- */
XMDF_API xid xfGetMetadataSpatialXML (xid a_Id, char **a_xml)
{
  xid SpatialGroup;
  int status;
  size_t len;
  char coord[25];
  double bottombound, leftbound, rightbound, topbound;

  /* open the Profile group */
  status = xfOpenGroup(a_Id, METADATA_MAIN_SPATIAL, &SpatialGroup);
  if (status < 0) {
    return status;
  }
  /* Read Bottom Bound */
  xfReadPropertyDouble(SpatialGroup, METADATA_SPATIAL_BOTTOMBOUND, 1, &bottombound);
  /* Read Bottom Bound */
  xfReadPropertyDouble(SpatialGroup, METADATA_SPATIAL_LEFTBOUND, 1, &leftbound);
  /* Read Bottom Bound */
  xfReadPropertyDouble(SpatialGroup, METADATA_SPATIAL_RIGHTBOUND, 1, &rightbound);
  /* Read Bottom Bound */
  xfReadPropertyDouble(SpatialGroup, METADATA_SPATIAL_TOPBOUND, 1, &topbound);

  xfCloseGroup(SpatialGroup);

  /* Create the XML string */
  len = 0;

  len += strlen(METATAG_MAIN_SPATIAL);
  *a_xml = (char *)malloc((len+1)*sizeof(char));
  strcpy(*a_xml, METATAG_MAIN_SPATIAL);

  len += strlen(METATAG_SPATIAL_BOTTOMBOUND);
  *a_xml = (char *)realloc(*a_xml, (len+1)*sizeof(char));
  strcat(*a_xml, METATAG_SPATIAL_BOTTOMBOUND);

  sprintf(coord, "%.6lf", bottombound);
  len += strlen(coord);
  *a_xml = (char *)realloc(*a_xml, (len+1)*sizeof(char));
  strcat(*a_xml, coord);

  len += strlen(METATAG_SPATIAL_BOTTOMBOUND_END);
  *a_xml = (char *)realloc(*a_xml, (len+1)*sizeof(char));
  strcat(*a_xml, METATAG_SPATIAL_BOTTOMBOUND_END);

  len += strlen(METATAG_SPATIAL_LEFTBOUND);
  *a_xml = (char *)realloc(*a_xml, (len+1)*sizeof(char));
  strcat(*a_xml, METATAG_SPATIAL_LEFTBOUND);

  sprintf(coord, "%.6lf", leftbound);
  len += strlen(coord);
  *a_xml = (char *)realloc(*a_xml, (len+1)*sizeof(char));
  strcat(*a_xml, coord);

  len += strlen(METATAG_SPATIAL_LEFTBOUND_END);
  *a_xml = (char *)realloc(*a_xml, (len+1)*sizeof(char));
  strcat(*a_xml, METATAG_SPATIAL_LEFTBOUND_END);

  len += strlen(METATAG_SPATIAL_RIGHTBOUND);
  *a_xml = (char *)realloc(*a_xml, (len+1)*sizeof(char));
  strcat(*a_xml, METATAG_SPATIAL_RIGHTBOUND);

  sprintf(coord, "%.6lf", rightbound);
  len += strlen(coord);
  *a_xml = (char *)realloc(*a_xml, (len+1)*sizeof(char));
  strcat(*a_xml, coord);

  len += strlen(METATAG_SPATIAL_RIGHTBOUND_END);
  *a_xml = (char *)realloc(*a_xml, (len+1)*sizeof(char));
  strcat(*a_xml, METATAG_SPATIAL_RIGHTBOUND_END);

  len += strlen(METATAG_SPATIAL_TOPBOUND);
  *a_xml = (char *)realloc(*a_xml, (len+1)*sizeof(char));
  strcat(*a_xml, METATAG_SPATIAL_TOPBOUND);

  sprintf(coord, "%.6lf", topbound);
  len += strlen(coord);
  *a_xml = (char *)realloc(*a_xml, (len+1)*sizeof(char));
  strcat(*a_xml, coord);

  len += strlen(METATAG_SPATIAL_TOPBOUND_END);
  *a_xml = (char *)realloc(*a_xml, (len+1)*sizeof(char));
  strcat(*a_xml, METATAG_SPATIAL_TOPBOUND_END);

  len += strlen(METATAG_MAIN_SPATIAL_END);
  *a_xml = (char *)realloc(*a_xml, (len+1)*sizeof(char));
  strcat(*a_xml, METATAG_MAIN_SPATIAL_END);

  return status;
} /* xfGetMetadataSpatialXML */

/* --------------------------------------------------------------------------- */
/* FUNCTION  xfCreateGeometricPathGroup */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfCreateGeometricPathGroup(xid a_ParentId, const char *a_Path,
               const char *a_Guid, int a_Compression,
               xid *a_PathGroup, double a_NullVal)
{
  int     error = 1;
    /* create the path group*/
  *a_PathGroup = xfpCreateGroup(a_ParentId, a_Path, GROUP_TYPE_GEOMETRIC_PATH);
  if (*a_PathGroup < 0) {
    xfpAddXMDFError("Failure creating Geometric Path Group");
    return *a_PathGroup;
  }
    /* write the guid*/
  error = xfpWriteDatasetString(*a_PathGroup, XF_GUID, a_Guid);
  if (error < 0) {
    xfpAddXMDFError("Failure writing GUID for Geometric Path Group");
    return error;
  }
    /* save the compression attribute*/
  error = xfpWriteAttributeInt(*a_PathGroup, GEOMPATH_ATT_COMPRESSION, 1, &a_Compression);
  if (error < 0) {
    xfpAddXMDFError("Failure writing compression type for Geometric Path Group");
    return error;
  }
    /* save the null value attribute*/
  error = xfpWriteAttributeDouble(*a_PathGroup, GEOMPATH_ATT_NULLVALUE, 1, &a_NullVal);
  if (error < 0) {
    xfpAddXMDFError("Failure writing NULLVALUE for Geometric Path Group");
    return error;
  }

  return error;
} /* xfCreateGeometricPathGroup */

/* --------------------------------------------------------------------------- */
/* FUNCTION  xfWriteParticleTimestep */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfWriteParticleTimestep(xid a_Id, int a_nDim, double a_Time,
               int a_nPaths, double *a_Locs)
{
  xid     LocsId, TimesId, minsId, maxsId;
  int     status = 1, i, j;
  hid_t   DataspaceId;
  herr_t (*old_func)(void*);
  void    *old_client_data;
  XDatasetParams Params;
  int      Compression;
  htri_t  IsSimple;
  int     Rank;
  hsize_t *Dims = NULL, newDims[3], Dim, *MaxDims = NULL;
  double  NullValue, mins[3], maxs[3];

  if (a_nPaths < 1) {
    xfpAddXMDFError("Error-Trying to save a timestep with 0 Geometric Paths");
    return -1;
  }

  H5Eget_auto1(&old_func, &old_client_data);

  /* Turn off error handling */
  H5Eset_auto1(NULL, NULL);

  /* Delete the time, values, mins, and maxs, activity datasets if they exist */
  LocsId = H5Dopen1(a_Id, GEOMPATH_DSET_LOCS);
  TimesId = H5Dopen1(a_Id, GEOMPATH_DSET_TIMES);

  /* Restore previous error handler */
  H5Eset_auto1(old_func, old_client_data);

  /* Get the compression level for the dataset */
  status = xfpReadAttributeInt(a_Id, GEOMPATH_ATT_COMPRESSION, 1, &Compression);
  if (status < 0) {
    if (LocsId > 0) {
      xfCloseGroup(LocsId);
    }
    if (TimesId > 0) {
      H5Dclose(TimesId);
    }
    return status;
  }

  /* Dataset time values */
  if (TimesId < 0) {
    /* Create the time dataset */
    xfpDsetParamsInit(&Params, 1, XTRUE, Compression);
    xfpDsetParamsSetSizes(&Params, 0, 1, H5S_UNLIMITED, 
                          GEOMPATH_TIME_CHUNK_SIZE);
    status = xfpWriteDsetDouble(a_Id, GEOMPATH_DSET_TIMES, &Params, &a_Time);
    if (status < 0) {
      xfpDsetParamsDestroy(&Params);
      xfpAddXMDFError("Error-Unable to write Geometric Path times");
      return status;
    }
  }
  else {
    H5Dclose(TimesId);
    /* Time dataset exists append to it */
    status = xfpAppendDset1DDouble(a_Id, GEOMPATH_DSET_TIMES, 1, &a_Time);
  }

    /* Dataset of point locations-First Timestep */
  if (LocsId < 0) {
    /* Create the values dataset */
    xfpDsetParamsReset(&Params, 3, XTRUE, Compression);
    /* Time dimension */
    xfpDsetParamsSetSizes(&Params, 0, 1, H5S_UNLIMITED, 1);
    /* Particle id dimension */
    xfpDsetParamsSetSizes(&Params, 1, a_nPaths, H5S_UNLIMITED, 20);
    /* Particle location dimension */
    xfpDsetParamsSetSizes(&Params, 2, a_nDim, a_nDim, a_nDim);
    /* get the null value to set as a fill value */
    status = xfpReadAttributeDouble(a_Id, GEOMPATH_ATT_NULLVALUE, 1, &NullValue);
    if (status < 0) {
      xfpDsetParamsDestroy(&Params);
      xfpAddXMDFError("Error-Unable to read Geometric Path NullValue");
      return status;
    }
    xfpDsetParamsUseFillValue(&Params, H5T_NATIVE_DOUBLE, 
                              (void *)(&NullValue));
    status = xfpWriteDsetDouble(a_Id, GEOMPATH_DSET_LOCS, &Params, a_Locs);
    xfpDsetParamsDestroy(&Params);
    if (status < 0) {
      return status;
    }
  }
    /* Dataset of point locations-Subsequent Timestep */
  else {
    /* Get the dataspace */
    DataspaceId = H5Dget_space(LocsId);
    /* Dataspace must be simply defined */
    IsSimple = H5Sis_simple(DataspaceId);
    if (!IsSimple) {
      H5Sclose(DataspaceId);
      H5Dclose(LocsId);
      xfpAddXMDFError("Error-Geometric Paths data is not simple");
      return -1;
    }
      /* Get the rank and dimensions for simple dataspace */
    status = xfpGetSimpleDataspaceInfo(DataspaceId, &Rank, &Dims, &MaxDims);
    if (status < 0 || Rank != a_nDim) {
      H5Sclose(DataspaceId);
      H5Dclose(LocsId);
      xfpAddXMDFError("Error-Can't get Geometric Paths data info");
      return -1;
    }
      /* see if we need to extend in the particle direction */
    if (Dims[1] != a_nPaths) {
        /* don't allow it to get smaller */
      if (Dims[1] > a_nPaths) {
        H5Sclose(DataspaceId);
        H5Dclose(LocsId);
        free(Dims);
        free(MaxDims);
        xfpAddXMDFError("Error-Can't reduce the number of Geometric Paths");
        return -1;
      }
        /* make the space bigger to accomodate the new particles */
      newDims[0] = Dims[0];
      newDims[1] = a_nPaths;
      newDims[2] = Dims[2];
      status = H5Dextend(LocsId, newDims);
      if (status < 0) {
        H5Sclose(DataspaceId);
        H5Dclose(LocsId);
        free(Dims);
        free(MaxDims);
        xfpAddXMDFError("Error-Can't extend the space of the Geometric Paths");
        return -1;
      }
  
    }
    free(Dims);
    free(MaxDims);
    H5Sclose(DataspaceId);
    H5Dclose(LocsId);
    /* values array exists append to it */
    xfpAppendDset3DDoubleFirstDim(a_Id, GEOMPATH_DSET_LOCS, 1, a_nPaths, a_nDim,
                                  a_Locs);
    if (status < 0) {
      return status;
    }
  }
    /* do the mins and maxs */
  H5Eget_auto1(&old_func, &old_client_data);

  /* Turn off error handling */
  H5Eset_auto1(NULL, NULL);

  minsId = H5Aopen_name(a_Id, GEOMPATH_ATT_MINS);
  maxsId = H5Aopen_name(a_Id, GEOMPATH_ATT_MAXS);

   /* Restore previous error handler */
  H5Eset_auto1(old_func, old_client_data);
    /* add the mins and maxes */
  if (minsId < 0 && maxsId < 0) {
    Dim = a_nDim;
    DataspaceId = H5Screate_simple(1, &Dim, &Dim);
    minsId = H5Acreate1(a_Id, GEOMPATH_ATT_MINS, xfpGetDefaultDoubleType(), 
                       DataspaceId, H5P_DEFAULT);
    maxsId = H5Acreate1(a_Id, GEOMPATH_ATT_MAXS, xfpGetDefaultDoubleType(), 
                       DataspaceId, H5P_DEFAULT);
    if (minsId < 0 || maxsId < 0) {
      H5Sclose(DataspaceId);
      xfpAddXMDFError("Error-Can't setup mins and maxes of Geometric Paths");
      return -1;
    }
    mins[0] = mins[1] = mins[2] = MAX_DOUBLE;
    maxs[0] = maxs[1] = maxs[2] = MIN_DOUBLE;
  }
    /* inconsistent mins and maxes - this is a problem */
  else if (minsId < 0 || maxsId < 0) {
    if (minsId > 0) {
      H5Aclose(minsId);
    }
    if (maxsId > 0) {
      H5Aclose(maxsId);
    }
    xfpAddXMDFError("Error-Inconsistend mins and maxes of Geometric Paths");
    return -1;
  }
    /* get the current mins and maxes */
  else {
    H5Aread(minsId, H5T_NATIVE_DOUBLE, mins);
    H5Aread(maxsId, H5T_NATIVE_DOUBLE, maxs);
  }
    /* enlarge the mins and maxes */
  for (i=0; i<a_nPaths*3; ++i) {
    j = i%3;
    if (a_Locs[i] < mins[j]) {
      mins[j] = a_Locs[i];
    } 
    if (a_Locs[i] > maxs[j]) {
      maxs[j] = a_Locs[i];
    } 
  }
  H5Awrite(minsId, H5T_NATIVE_DOUBLE, mins);
  H5Awrite(maxsId, H5T_NATIVE_DOUBLE, maxs);
  H5Aclose(minsId);
  H5Aclose(maxsId);

  return status;
} /* xfWriteParticleTimestep */

/* --------------------------------------------------------------------------- */
/* FUNCTION  xfGetPathNullVal */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfGetPathNullVal(xid GroupId, double *NullVal)
{
  return xfpReadAttributeDouble(GroupId, GEOMPATH_ATT_NULLVALUE, 1, NullVal);
} /* xfGetPathNullVal */

/* --------------------------------------------------------------------------- */
/* FUNCTION  xfGetNumberOfPaths */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfGetNumberOfPaths(xid GroupId, int *NumPaths)
{
  int     Rank;
  hsize_t *Dims = NULL, *Maxdims = NULL;
  hid_t   DsetId, DataspaceId;
  herr_t  status;
  htri_t  IsSimple;

  DsetId = H5Dopen1(GroupId, GEOMPATH_DSET_LOCS);
  if (DsetId < 0) {
   return DsetId;
  }

  /* Get the dataspace */
  DataspaceId = H5Dget_space(DsetId);

  /* Dataspace must be simply defined */
  IsSimple = H5Sis_simple(DataspaceId);
  if (!IsSimple) {
    H5Sclose(DataspaceId);
    H5Dclose(DsetId);
    return ERROR_ATTRIBUTE_NOT_SUPPORTED;
  }

  /* Get the rank and dimensions for simple dataspace */
  status = xfpGetSimpleDataspaceInfo(DataspaceId, &Rank, &Dims, &Maxdims);
  if (status > 0 && Rank == 3) {
    *NumPaths = (int)Dims[1];
    free(Dims);
    free(Maxdims);
    H5Sclose(DataspaceId);
    H5Dclose(DsetId);
    return 1;
  }
  if (Dims) {
    free(Dims);
  }
  if (Maxdims) {
    free(Maxdims);
  }
  H5Sclose(DataspaceId);
  H5Dclose(DsetId);
  return -1;
} /* xfGetNumberOfPaths */

/* --------------------------------------------------------------------------- */
/* FUNCTION  xfGetNumberOfTimes */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfGetNumberOfTimes(xid GroupId, int *NumTimes)
{
  int     Rank;
  hsize_t *Dims = NULL, *Maxdims = NULL;
  hid_t   DsetId, DataspaceId;
  herr_t  status;
  htri_t  IsSimple;

  DsetId = H5Dopen1(GroupId, GEOMPATH_DSET_TIMES);
  if (DsetId < 0) {
   return DsetId;
  }

  /* Get the dataspace */
  DataspaceId = H5Dget_space(DsetId);

  /* Dataspace must be simply defined */
  IsSimple = H5Sis_simple(DataspaceId);
  if (!IsSimple) {
    H5Sclose(DataspaceId);
    H5Dclose(DsetId);
    return ERROR_ATTRIBUTE_NOT_SUPPORTED;
  }

  /* Get the rank and dimensions for simple dataspace */
  status = xfpGetSimpleDataspaceInfo(DataspaceId, &Rank, &Dims, &Maxdims);
  if (status > 0 && Rank == 1) {
    *NumTimes = (int)Dims[0];
    free(Dims);
    free(Maxdims);
    H5Sclose(DataspaceId);
    H5Dclose(DsetId);
    return 1;
  }
  if (Dims) {
    free(Dims);
  }
  if (Maxdims) {
    free(Maxdims);
  }
  H5Sclose(DataspaceId);
  H5Dclose(DsetId);
  return -1;
} /* xfGetNumberOfTimes */

/* --------------------------------------------------------------------------- */
/* FUNCTION  xfGetPathDimensionality */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfGetPathDimensionality(xid GroupId, int *NumDims)
{
  int     Rank;
  hsize_t *Dims = NULL, *Maxdims = NULL;
  hid_t   DsetId, DataspaceId;
  herr_t  status;
  htri_t  IsSimple;

  DsetId = H5Dopen1(GroupId, GEOMPATH_DSET_LOCS);
  if (DsetId < 0) {
   return DsetId;
  }

  /* Get the dataspace */
  DataspaceId = H5Dget_space(DsetId);

  /* Dataspace must be simply defined */
  IsSimple = H5Sis_simple(DataspaceId);
  if (!IsSimple) {
    H5Sclose(DataspaceId);
    H5Dclose(DsetId);
    return ERROR_ATTRIBUTE_NOT_SUPPORTED;
  }

  /* Get the rank and dimensions for simple dataspace */
  status = xfpGetSimpleDataspaceInfo(DataspaceId, &Rank, &Dims, &Maxdims);
  if (status > 0 && Rank == 3) {
    *NumDims = (int)Dims[2];
    free(Dims);
    free(Maxdims);
    H5Sclose(DataspaceId);
    H5Dclose(DsetId);
    return 1;
  }
  if (Dims) {
    free(Dims);
  }
  if (Maxdims) {
    free(Maxdims);
  }
  H5Sclose(DataspaceId);
  H5Dclose(DsetId);
  return -1;
} /* xfGetPathDimensionality */

/* --------------------------------------------------------------------------- */
/* FUNCTION  xfGetPathTimesArray */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfGetPathTimesArray(xid GroupId, int NumTimes, double *Times)
{
  return xfpReadDset1DDouble(GroupId, GEOMPATH_DSET_TIMES, NumTimes, Times);
} /* xfGetPathTimesArray */

/* --------------------------------------------------------------------------- */
/* FUNCTION  xfReadPathLocationsAtTime */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfReadPathLocationsAtTime(xid GroupId, int TimeIndex,
                                      int FirstPathIndex,
                                      int NumIndicies, double *Locs)
{
  int nDims = 0;

  if (xfGetPathDimensionality(GroupId, &nDims) < 0) {
    xfpAddXMDFError("Error-Unable to get particle dimensionality");
    return -1;
  }
  return xfpReadDset3DDoublePortion(GroupId, GEOMPATH_DSET_LOCS, TimeIndex-1, 1,
            FirstPathIndex-1, NumIndicies, 0, nDims, Locs);
} /* xfReadPathLocationsAtTime */

/* --------------------------------------------------------------------------- */
/* FUNCTION  xfReadPathLocationsForParticle */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfReadPathLocationsForParticle(xid GroupId, int PathIndex, 
                                   int FirstTimeIndex, int NumTimes,
                                   double *Locs)
{
  int nDims = 0;

  if (xfGetPathDimensionality(GroupId, &nDims) < 0) {
    xfpAddXMDFError("Error-Unable to get particle dimensionality");
    return -1;
  }
  return xfpReadDset3DDoublePortion(GroupId, GEOMPATH_DSET_LOCS, FirstTimeIndex-1,
            NumTimes, PathIndex-1, 1, 0, nDims, Locs);
} /* xfReadPathLocationsForParticle */

/* --------------------------------------------------------------------------- */
/* FUNCTION  xfReadPathLocationsForParticles */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfReadPathLocationsForParticles(xid GroupId, int NumPaths,
                                    const int *PathIndices,
                                    int FirstTimeIndex, int NumTimes,
                                    double *Locs)
{
  hssize_t     *nStart1, *nStart2, *nStart3;
  hsize_t      *nNum1, *nNum2, *nNum3; 
  int           i, nStatus, nDims = 0;

  if (xfGetPathDimensionality(GroupId, &nDims) < 0) {
    xfpAddXMDFError("Error-Unable to get particle dimensionality");
    return -1;
  }
  /* allocate the bounds for the particle limits */
  nStart1 = (hssize_t *)malloc(NumPaths*sizeof(hssize_t));
  nStart2 = (hssize_t *)malloc(NumPaths*sizeof(hssize_t));
  nStart3 = (hssize_t *)malloc(NumPaths*sizeof(hssize_t));
  nNum1 = (hsize_t *)malloc(NumPaths*sizeof(hsize_t));
  nNum2 = (hsize_t *)malloc(NumPaths*sizeof(hsize_t));
  nNum3 = (hsize_t *)malloc(NumPaths*sizeof(hsize_t));
    /* fill in the bounds */
  if (nStart1 && nStart2 && nStart3 && nNum1 && nNum2 && nNum3) {
    for (i=0; i<NumPaths; ++i) {
      nStart1[i] = FirstTimeIndex-1;
      nNum1[i]   = NumTimes;
      nStart2[i] = PathIndices[i]-1;
      nNum2[i]   = 1;
      nStart3[i] = 0;
      nNum3[i]   = nDims;
    }
    nStatus = xfpReadDset3DDoublePortions(GroupId, GEOMPATH_DSET_LOCS, NumPaths,
        nStart1, nNum1, nStart2, nNum2, nStart3, nNum3, Locs);
    free(nStart1);
    free(nStart2);
    free(nStart3);
    free(nNum1);
    free(nNum2);
    free(nNum3);
    return nStatus;
  }
    /* clean up after malloc failure */
  if (nStart1) {
    free(nStart1);
  }
  if (nStart2) {
    free(nStart2);
  }
  if (nStart3) {
    free(nStart3);
  }
  if (nNum1) {
    free(nNum1);
  }
  if (nNum2) {
    free(nNum2);
  }
  if (nNum3) {
    free(nNum3);
  }
  return -1;
} /* xfReadPathLocationsForParticles */
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfSetupToWriteDatasets */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfSetupToWriteDatasets (const char *a_Filename,
                         const char *a_MultiDatasetsGroupPath, 
                         const char *a_PathInMultiDatasetsGroup,
                         const char *a_SpatialDataObjectGuid,
                         int a_OverwriteOptions, xid *a_FileId, xid *a_GroupId)
{
  return xfSetupToWriteDatasets2(a_Filename, a_MultiDatasetsGroupPath,
           a_PathInMultiDatasetsGroup, a_SpatialDataObjectGuid,
           a_OverwriteOptions, a_FileId, a_GroupId, XFALSE);
} /* xfSetupToWriteDatasets */
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfSetupToWriteDatasets */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfSetupToWriteDatasets2 (const char *a_Filename,
                         const char *a_MultiDatasetsGroupPath, 
                         const char *a_PathInMultiDatasetsGroup,
                         const char *a_SpatialDataObjectGuid,
                         int a_OverwriteOptions, xid *a_FileId, xid *a_GroupId,
                         xmbool a_inMemory)
{
  xid    xFileId = 0, xMultiDatasetsId = 0, xStartDatasets = 0;
  xmbool  bFileNew = XFALSE, bMultiDatasetsNew = XFALSE, bFail = XFALSE;
  char   TestGuid[XF_GUID_STRINGLENGTH];
  int    status;

  /* if overwrite file is not on try opening the file */
  bFileNew = XFALSE;
  if (a_OverwriteOptions != XF_OVERWRITE_CLEAR_FILE) {
    if (xfOpenFile(a_Filename, &xFileId, XFALSE) < 0) {
      /* File could not be opened create */
      bFileNew = XTRUE;
      status = xfiCreateFile(a_Filename, &xFileId, XTRUE, a_inMemory);
      if (status < 0) {
        return status;
      }
    }
  }
  else {
    /* if overwrite file option is on create the file with overwrite on */
    bFileNew = XTRUE;
    status = xfiCreateFile(a_Filename, &xFileId, XTRUE, a_inMemory);
    if (status < 0) {
      return status;
    }
  }

  /* if the file is not new see if the path to the multidatasets group exists */
  if (bFileNew == XFALSE) {
    status = xfOpenGroup(xFileId, a_MultiDatasetsGroupPath, &xMultiDatasetsId);
    if (status >= 0 && xMultiDatasetsId != 0) {
      /* if the path to multidatasets exists make sure it is correct type */
      /* and that the GUIDs match, if not give an error */
      if (xfpIsGroupOfType(xMultiDatasetsId, GROUP_TYPE_MULTI_DATASETS)
                                                                 == XFALSE) {
        bFail = XTRUE;
      }
      else {
        if (xfGetDatasetsSdoGuid(xMultiDatasetsId, TestGuid) < 0 ||
            strcmp(TestGuid, a_SpatialDataObjectGuid) != 0) {
          bFail = XTRUE;
        }
      }
      if (bFail) {
        xfCloseFile(xFileId);
        return -1;
      }
    }
    else {
      /* if the path to multidatasets does not exist create the group */
      bMultiDatasetsNew = XTRUE;
      if (xfCreateMultiDatasetsGroup(xFileId, a_MultiDatasetsGroupPath,
                                     a_SpatialDataObjectGuid,
                                     &xMultiDatasetsId) < 0) {
        xfCloseFile(xFileId);
        return -1;
      }
    }
  }
  else {
    /* if the file is new the multi datasets group needs to be created */
    bMultiDatasetsNew = XTRUE;
    if (xfCreateMultiDatasetsGroup(xFileId, a_MultiDatasetsGroupPath,
                                   a_SpatialDataObjectGuid,
                                   &xMultiDatasetsId) < 0) {
      xfCloseFile(xFileId);
      return -1;
    }
  }

  /* if the multidatasets group and the start path are the same */
  if (strcmp(a_PathInMultiDatasetsGroup, "") == 0) {
    /* if we are supposed to clear the group do it */
    if (a_OverwriteOptions == XF_OVERWRITE_CLEAR_DATASET_GROUP) {
      xfpClearGroup(xMultiDatasetsId);
      /* when we clear the group the Guid gets removed.  Replace the Guid. */
      xfpWriteDatasetString(xMultiDatasetsId, XF_GUID, a_SpatialDataObjectGuid);
    }
    *a_FileId = xFileId;
    *a_GroupId = xMultiDatasetsId;
    return XTRUE;
  }

  /* if the file is not new and the path to the multidatasets is not new, */
  /* see if the path to start writing datasets exists */
  if (bFileNew == XFALSE && bMultiDatasetsNew == XFALSE &&
      xfOpenGroup(xMultiDatasetsId, a_PathInMultiDatasetsGroup, &xStartDatasets)
                                                                         >= 0) {
      /* if the path to start writing datasets exists and the overwrite */
      /* everything in path is selected clear the folder */
    if (a_OverwriteOptions == XF_OVERWRITE_CLEAR_DATASET_GROUP) {
      xfpClearGroup(xStartDatasets);
    }
  }
  else { 
    /* if the path to start writing datasets does not exist then create it */
    if (xfCreateGenericGroup(xMultiDatasetsId, a_PathInMultiDatasetsGroup,
                                                         &xStartDatasets) < 0) {
      xfCloseGroup(xMultiDatasetsId);
      xfCloseFile(xFileId);
      return -1;
    }
  }

  /* close the multidatasets group we no longer need it */
  xfCloseGroup(xMultiDatasetsId);

  /* return the id to the file, and the group to start writing in */
  *a_FileId = xFileId;
  *a_GroupId = xStartDatasets;
  return XTRUE;
} /* xfSetupToWriteDatasets */
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfCreateMultiDatasetsGroup */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfCreateMultiDatasetsGroup (xid a_Id, const char *a_Path,
                                        const char *a_Guid, xid *a_MultiId)
{
  *a_MultiId = xfpCreateGroup(a_Id, a_Path, GROUP_TYPE_MULTI_DATASETS);
  if (*a_MultiId < 0) {
    return *a_MultiId;
  }
  return xfpWriteDatasetString(*a_MultiId, XF_GUID, a_Guid);
} /* xfCreateMultiDatasetsGroup */
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfGetGroupPathsSizeForMultiDatasets */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfGetGroupPathsSizeForMultiDatasets (xid a_Id, int *Num,
                                                 int *Maxsize)
{
  return xfpNumGroupsOfType(a_Id, GROUP_TYPE_MULTI_DATASETS, Num, Maxsize,
                            XFALSE);
} /* xfGetGroupPathsSizeForMultiDatasets */
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfGetAllGroupPathsForMultiDatasets */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfGetAllGroupPathsForMultiDatasets (xid a_Id, int a_Num,
                                                int a_Maxsize, char *a_Paths)
{
  return xfpPathGroupsOfType(a_Id, GROUP_TYPE_MULTI_DATASETS, a_Num, a_Maxsize,
                             a_Paths, XFALSE);
} /* xfGetAllGroupPathsForMultiDatasets */
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfGetDatasetsSdoGuid */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfGetDatasetsSdoGuid (xid a_MultiDatasetsGroup, char *a_GUID)
{
  return xfReadPropertyString(a_MultiDatasetsGroup, XF_GUID, 1,
                              XF_GUID_STRINGLENGTH, a_GUID);
} /* xfGetDatasetsSdoGuid */
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfOpenMultiDatasetsGroup */
/*! PURPOSE:   Open or create a multi-datasets group inside a mesh or grid
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfOpenMultiDatasetsGroup (xid a_Id, xid *DatasetsGroupId)
{
  char  strGuid[XF_GUID_STRINGLENGTH];

  /* make sure that the group is a mesh or grid */
  if (!xfpIsGroupOfType(a_Id, GROUP_TYPE_MESH) &&
      !xfpIsGroupOfType(a_Id, GROUP_TYPE_GRID)) {
    return -1;
  }

  /* see if the group exists */
  if (xfOpenGroup(a_Id, MULTI_DATASET_LOCATION, DatasetsGroupId) >= 0 &&
      *DatasetsGroupId > 0) {
    return *DatasetsGroupId;
  }

  /* the group does not exist, see if SpatialDataObject has a GUID */
  if (xfReadPropertyString(a_Id, "PROPERTIES/Guid", 1, XF_GUID_STRINGLENGTH, strGuid)
                                                                         >= 0) {
    /* create the group with the GUID */
    xfCreateMultiDatasetsGroup(a_Id, MULTI_DATASET_LOCATION, strGuid,
                               DatasetsGroupId);
  }
  else {
    *DatasetsGroupId = xfpCreateGroup(a_Id, MULTI_DATASET_LOCATION,
                                      GROUP_TYPE_MULTI_DATASETS);
  }
  return *DatasetsGroupId;
} /* xfOpenMultiDatasetsGroup */
/*-----------------------------------------------------------------------------*/
/* FUNCTION  xfCreateScalarDataset*/
/*! PURPOSE:   Create a scalar dataset
 * - NOTES:     The intermediate groups in the path may or may not be created. */
/*-----------------------------------------------------------------------------*/
XMDF_API xid xfCreateScalarDataset(xid a_DatasetsGroupId, const char *a_Path,
                    const char *a_Units, const char *a_TimeUnits,
                    int a_Compression, xid *a_DatasetId)
{
  xid         DatasetGroup = NONE, ParentId = NONE;
  int         status = -1, i=0;
  int         PathLen;
  const char *strGroups = NULL;
  size_t      EndOfGroup;
  char        *ChildPath = NULL, *ParentPath = NULL;
  herr_t     (*old_func)(void*);
  void        *old_client_data;

  H5Eget_auto1(&old_func, &old_client_data);

    /* Turn off error handling*/
  H5Eset_auto1(NULL, NULL);

    /* Check to see if group exists*/
  status = xfOpenGroup(a_DatasetsGroupId, a_Path, &DatasetGroup);
  xfCloseGroup(DatasetGroup);
  
    /* Restore previous error handler*/
  H5Eset_auto1(old_func, old_client_data);

    /* if group already exists, delete it and rewrite it.*/
    /* (in effect, we're overwriting it)*/
  if (status >= 0) {
    strGroups = a_Path;
    PathLen = strlen(strGroups);
    ParentPath = (char *)malloc(strlen(a_Path)+1*sizeof(char));
    ChildPath = (char *)malloc(strlen(a_Path)+1*sizeof(char));

      /* get the first substring until the end or until we hit a slash "/"*/
    EndOfGroup = strcspn(strGroups, "/");
    strncpy(ParentPath, strGroups, EndOfGroup);
    ParentPath[EndOfGroup] = '\0';

      /* get the substring after the slash "/"*/
    for (i=(EndOfGroup+1); i < PathLen; i++) {
      ChildPath[i-EndOfGroup-1] = strGroups[i];
    }
    ChildPath[PathLen-EndOfGroup-1] = '\0';

      /* delete the group so it can be rewritten*/
    status = xfOpenGroup(a_DatasetsGroupId, ParentPath, &ParentId);
    H5Gunlink(ParentId, ChildPath);
    status = xfCloseGroup(ParentId);
    if (ChildPath) {
      free(ChildPath);
      ChildPath = NULL;
    }
    if (ParentPath) {
      free(ParentPath);
      ParentPath = NULL;
    }
  } /* if (status >= 0)*/

    /* create the HDF5 group*/
  DatasetGroup = xfpCreateGroup(a_DatasetsGroupId, a_Path, 
                                                    GROUP_TYPE_DATASET_SCALAR);
  if (DatasetGroup < 0) {
    return DatasetGroup;
  }
  *a_DatasetId = DatasetGroup;

    /* set the time units*/
  status = xfpWriteAttributeString(DatasetGroup, TIME_UNITS, a_TimeUnits);
  if (status < 0) {
    xfCloseGroup(DatasetGroup);
    return status;
  }

    /* set the units*/
  status = xfpWriteAttributeString(DatasetGroup, DATASET_ATT_UNITS, a_Units);
  if (status < 0) {
    xfCloseGroup(DatasetGroup);
    return status;
  }

    /* store the compression level*/
  status = xfpWriteAttributeInt(DatasetGroup, DATASET_ATT_COMPRESSION, 1,
                                &a_Compression); 
  if (status < 0) {
    xfCloseGroup(DatasetGroup);
    return status;
  }

  return DatasetGroup;
} /* xfCreateScalarDataset */

/******************************************************************************/
/*   FUNCTION  xfCreateScalarDatasetExtendable */
/*! PURPOSE:   Create a scalar dataset that can be extended
 *   - NOTES:     The intermediate groups in the path may or may not be created. */
/******************************************************************************/
XMDF_API xid xfCreateScalarDatasetExtendable(xid a_DatasetsGroupId, const char *a_Path,
                    const char *a_Units, const char *a_TimeUnits, float a_FillVal,
                    int a_Compression, xid *a_DatasetId)
{
    /* create the data set */
  if (xfCreateScalarDataset(a_DatasetsGroupId, a_Path, a_Units, a_TimeUnits,
    a_Compression, a_DatasetId) < 0) {
    return -1;
  }
    /* make it extendable */
  if (xfpWriteAttributeFloat(*a_DatasetId, DATASET_ATT_EXTENDFILL, 1, 
         &a_FillVal) < 0) {
    xfpAddXMDFError("Error: couldn't add the fill value to this stupid data set");
    return -1;
  }

  return 1;
} /* xfCreateScalarDatasetExtendable */

/******************************************************************************/
/*   FUNCTION  xfExtendScalarDataset */
/*! PURPOSE:   Create a scalar dataset that can be extended
 *   - NOTES:     The intermediate groups in the path may or may not be created. */
/******************************************************************************/
XMDF_API xid xfExtendScalarDataset (xid a_Id, int aNewSize)
{
  xid      ValsId;
  hsize_t *Dims = NULL, newDims[2], *MaxDims = NULL;
  hid_t    DataspaceId;
  int      Rank, status;
  htri_t   IsSimple;

  ValsId = H5Dopen1(a_Id, DATASET_DSET_VALUES);
  if (ValsId < 0) {
    xfpAddXMDFError("Error: couldn't extend this stupid data set");
    return -1;
  }

    /* Get the dataspace */
  DataspaceId = H5Dget_space(ValsId);
    /* Dataspace must be simply defined */
  IsSimple = H5Sis_simple(DataspaceId);
  if (!IsSimple) {
    H5Sclose(DataspaceId);
    H5Dclose(ValsId);
    xfpAddXMDFError("Error-Dataset Values are not simple");
    return -1;
  }
      /* Get the rank and dimensions for simple dataspace */
  status = xfpGetSimpleDataspaceInfo(DataspaceId, &Rank, &Dims, &MaxDims);
  if (status < 0 || Rank != 2) {
    H5Sclose(DataspaceId);
    H5Dclose(ValsId);
    xfpAddXMDFError("Error-Can't get dataset info or rank is wrong");
    return -1;
  }
      /* make the dataset bigger to accomodate the new size */
  newDims[0] = Dims[0];
  newDims[1] = aNewSize;
  status = H5Dextend(ValsId, newDims);
      /* close and free the stuff now */
  H5Sclose(DataspaceId);
  H5Dclose(ValsId);
  free(Dims);
  free(MaxDims);
  if (status < 0) {
    xfpAddXMDFError("Error-Can't extend the space of the Data Set");
  }
  return status;
} /* xfExtendScalarDataset */

/******************************************************************************/
/*   FUNCTION  xfCreateVectorDataset */
/*! PURPOSE:   Create a vector dataset
 *   - NOTES:     The intermediate groups in the path may or may not be created. */
/******************************************************************************/
XMDF_API xid xfCreateVectorDataset(xid a_DatasetsGroupId, const char *a_Path,
                    const char *a_Units, const char *a_TimeUnits,
                    int a_Compression, xid *a_DatasetId)
{
  xid         DatasetGroup = NONE, ParentId = NONE;
  int         status = -1, i=0;
  int         PathLen;
  const char *strGroups = NULL;
  size_t      EndOfGroup;
  char        *ChildPath, *ParentPath; /*, *PPath;*/
  herr_t     (*old_func)(void*);
  void        *old_client_data;

  ChildPath = "";
  ParentPath = "";
/*   PPath = ""; */

  H5Eget_auto1(&old_func, &old_client_data);

    /* Turn off error handling*/
  H5Eset_auto1(NULL, NULL);

    /* Check to see if group exists*/
  status = xfOpenGroup(a_DatasetsGroupId, a_Path, &DatasetGroup);
  xfCloseGroup(DatasetGroup);
  
    /* Restore previous error handler*/
  H5Eset_auto1(old_func, old_client_data);

    /* if group already exists, delete it and rewrite it.*/
    /* (in effect, we're overwriting it)*/
  if (status >= 0) {
    strGroups = a_Path;
    PathLen = strlen(strGroups);
    ParentPath = (char *)malloc(strlen(a_Path)+1*sizeof(char));
    ChildPath = (char *)malloc(strlen(a_Path)+1*sizeof(char));

      /* get the first substring until the end or until we hit a slash "/"*/
    EndOfGroup = strcspn(strGroups, "/");
    strncpy(ParentPath, strGroups, EndOfGroup);
    ParentPath[EndOfGroup] = '\0';

      /* get the substring after the slash "/"*/
    for (i=(EndOfGroup+1); i < PathLen; i++) {
      ChildPath[i-EndOfGroup-1] = strGroups[i];
    }
    ChildPath[PathLen-EndOfGroup-1] = '\0';

      /* delete the group so it can be rewritten*/
    status = xfOpenGroup(a_DatasetsGroupId, ParentPath, &ParentId);
    H5Gunlink(ParentId, ChildPath);
    status = xfCloseGroup(ParentId);
    if (ChildPath) {
      free(ChildPath);
    }
    if (ParentPath) {
      free(ParentPath);
    }
  } /* if (status >= 0)*/

  /* create the HDF5 group   */
  DatasetGroup = xfpCreateGroup(a_DatasetsGroupId, a_Path, 
                                                    GROUP_TYPE_DATASET_VECTOR);
  if (DatasetGroup < 0) {
    return DatasetGroup;
  }
  *a_DatasetId = DatasetGroup;

  /* set the time units  */
  status = xfpWriteAttributeString(DatasetGroup, TIME_UNITS, a_TimeUnits);
  if (status < 0) {
    xfCloseGroup(DatasetGroup);
    return status;
  }

  /* set the units */
  status = xfpWriteAttributeString(DatasetGroup, DATASET_ATT_UNITS, a_Units);
  if (status < 0) {
    xfCloseGroup(DatasetGroup);
    return status;
  }

  /* store the compression level */
  status = xfpWriteAttributeInt(DatasetGroup, DATASET_ATT_COMPRESSION, 
                                1, &a_Compression); 
  if (status < 0) {
    xfCloseGroup(DatasetGroup);
    return status;
  }

  return DatasetGroup;
} /* xfCreateVectorDataset */
/******************************************************************************/
/*   FUNCTION  xfWriteScalarTimestep */
/*! PURPOSE:
 *   - NOTES:     Append a timestep to a float scalar dataset. */
/******************************************************************************/
XMDF_API xid xfWriteScalarTimestep(xid a_Id, double a_Time, int 
                    a_NumValues, const float *a_Values)
{
  
  return xfiWriteScalarTimestep(a_Id, a_Time, a_NumValues, 
                         (void*)a_Values, DATASET_TYPE_FLOAT);

} /* xfWriteScalarTimestep */
/******************************************************************************/
/*   FUNCTION  xfWriteScalarTimestepMinMax */
/*! PURPOSE:
 *   - NOTES:     Append a timestep to a float scalar dataset. 
 *             The min and max should be set to the min or max of the timestep*/
/******************************************************************************/
XMDF_API xid xfWriteScalarTimestepMinMax (xid a_Id, double a_Time, int 
                    a_NumValues, const float *a_Values, float a_Min, 
                    float a_Max)
{
  return xfiWriteScalarTimestepMinMax(a_Id, a_Time, a_NumValues,
                               (void*)a_Values, DATASET_TYPE_FLOAT,
                               (double)a_Min, (double)a_Max);

} /* xfWriteScalarTimestepMinMax */
/******************************************************************************/
/*   FUNCTION  xfWriteScalarTimestepFloat */
/*! PURPOSE:
 *   - NOTES:     Append a timestep to a float scalar dataset. */
/******************************************************************************/
XMDF_API xid xfWriteScalarTimestepFloat(xid a_Id, double a_Time, int 
                    a_NumValues, const float *a_Values)
{
  
  return xfiWriteScalarTimestep(a_Id, a_Time, a_NumValues, 
                         (void*)a_Values, DATASET_TYPE_FLOAT);

} /* xfWriteScalarTimestepFloat */
/******************************************************************************/
/*   FUNCTION  xfWriteScalarTimestepMinMaxFloat */
/*! PURPOSE:
 *   - NOTES:     Append a timestep to a float scalar dataset. 
 *             The min and max should be set to the min or max of the timestep */
/******************************************************************************/
XMDF_API xid xfWriteScalarTimestepMinMaxFloat (xid a_Id, double a_Time, int 
                    a_NumValues, const float *a_Values, float a_Min, 
                    float a_Max)
{
  return xfiWriteScalarTimestepMinMax(a_Id, a_Time, a_NumValues,
                               (void*)a_Values, DATASET_TYPE_FLOAT,
                               (double)a_Min, (double)a_Max);

} /* xfWriteScalarTimestepMinMaxFloat */
/******************************************************************************/
/*   FUNCTION  xfWriteScalarTimestepDouble */
/*! PURPOSE:
 *   - NOTES:     Append a timestep to a double scalar dataset. */
/******************************************************************************/
XMDF_API xid xfWriteScalarTimestepDouble(xid a_Id, double a_Time, int 
                    a_NumValues, const double *a_Values)
{
  
  return xfiWriteScalarTimestep(a_Id, a_Time, a_NumValues, 
                         (void*)a_Values, DATASET_TYPE_DOUBLE);

} /* xfWriteScalarTimestepFloat */
/******************************************************************************/
/*   FUNCTION  xfWriteScalarTimestepMinMaxDouble */
/*! PURPOSE:
 *   - NOTES:     Append a timestep to a double scalar dataset. 
 *             The min and max should be set to the min or max of the timestep */
/******************************************************************************/
XMDF_API xid xfWriteScalarTimestepMinMaxDouble (xid a_Id, double a_Time, int 
                    a_NumValues, const double *a_Values, double a_Min, 
                    double a_Max)
{
  return xfiWriteScalarTimestepMinMax(a_Id, a_Time, a_NumValues,
                               (void*)a_Values, DATASET_TYPE_DOUBLE,
                               a_Min, a_Max);

} /* xfWriteScalarTimestepMinMaxDouble */
/******************************************************************************/
/*   FUNCTION  xfiCreateFile */
/*! PURPOSE: */
/******************************************************************************/
xid xfiCreateFile (const char * a_Filename, xid *Id, xmbool a_Overwrite,
                   xmbool a_inMemory)
{
  htri_t    IsHdf5;
  int       status;
  float     version = (float)XMDF_VERSION;
  hid_t     accessProperties;
  size_t    memory_increment = 10000;
  hbool_t   backing_store = XFALSE;

  *Id = -1;

  /* Initialize the double and float types*/
  if (!xfpGetInitialized()) {
    xfpInitialize();
  }

  /* clear the error stack in case of previous errors*/
  xfpClearErrors();
  
  H5Eset_auto1(xfpHDF5ErrorHandler, NULL);

  /* For some reason adding this makes everything work */
  IsHdf5 = H5Fis_hdf5(a_Filename);
  IsHdf5 = IsHdf5;

  accessProperties = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fclose_degree(accessProperties, H5F_CLOSE_STRONG);
  if (a_inMemory) {
    status = H5Pset_fapl_core(accessProperties, memory_increment,
                              backing_store);
    if (status < 0) {
      return status;
    }
  }

  /* Create the file*/
  if (a_Overwrite) {
     *Id = H5Fcreate(a_Filename, H5F_ACC_TRUNC, H5P_DEFAULT, accessProperties);
  }
  else {
    *Id = H5Fcreate(a_Filename, H5F_ACC_EXCL, H5P_DEFAULT, accessProperties);
  }
  H5Pclose(accessProperties);

  if (*Id >= 0) {
    /* write a string attribute specifying the file is an Xmdf file*/
    status = xfpWriteDatasetString(*Id, FILE_TYPE_ID, FILE_TYPE_XMDF);
    if (status < 0) {
      H5Fclose(*Id);
      *Id = status;
    }

    status = xfWritePropertyFloat(*Id, FILE_VERSION, 1, &version, NONE);
    if (status < 0) {
      H5Fclose(*Id);
      *Id = status;
    }
  }

  return *Id;
} // xfiCreateFile
/******************************************************************************/
/*   FUNCTION  xfiDetermineTimestepMinsAndMaxs */
/*! PURPOSE:
 *   - NOTES:     Determine the minimum and maximum values for a dataset
                  Pass NULL for a_NullValue if not using a NULL value */
/******************************************************************************/
static xid xfiDetermineTimestepMinsAndMaxs (xid a_Id, int a_NumValues, 
            const void *a_Values, int a_DataType, double *a_min, double *a_max)
{
  float            fNullValue = 0.0;
  double           dNullValue = 0.0;
  int              iValue = 0;
  xmbool           bUseNullValue = XFALSE;
  herr_t           (*old_func)(void*);
  void            *old_client_data = NULL;
  xid              propGroupId;
  int              status = -1;

  if (a_NumValues < 1) {
    return ERROR_DATASET_NO_DATA;
  }

  H5Eget_auto1(&old_func, &old_client_data);

  /* Turn off error handling */
  H5Eset_auto1(NULL, NULL);

    /* get the null value if there is one */
  if (xfOpenPropertyGroup(a_Id, &propGroupId) >= 0) {
    if (a_DataType == DATASET_TYPE_FLOAT) {
      status = xfReadPropertyFloat(propGroupId, PROP_NULL_VALUE, 1, &fNullValue);
    }
    else if (a_DataType == DATASET_TYPE_DOUBLE) {
      status = xfReadPropertyDouble(propGroupId, PROP_NULL_VALUE, 1, &dNullValue);
    }
    if (status >= 0) {
      bUseNullValue = XTRUE;
    }
    xfCloseGroup(propGroupId);
  }

  /* Restore previous error handler */
  H5Eset_auto1(old_func, old_client_data);

  *a_min = MAX_DOUBLE;
  *a_max = MIN_DOUBLE;
  for (iValue = 0; iValue < a_NumValues; iValue++) {      
    /* get min and max values */
    if (a_DataType == DATASET_TYPE_FLOAT) {
      if (bUseNullValue == XFALSE || ((float*)a_Values)[iValue] != fNullValue) {
        *a_min = Xmin(*a_min, ((float*)a_Values)[iValue]);
        *a_max = Xmax(*a_max, ((float*)a_Values)[iValue]);
      }
    }
    else if (a_DataType == DATASET_TYPE_DOUBLE) {
      if (bUseNullValue == XFALSE || ((double*)a_Values)[iValue] != dNullValue) {
        *a_min = Xmin(*a_min, ((double*)a_Values)[iValue]);
        *a_max = Xmax(*a_max, ((double*)a_Values)[iValue]);
      }
    }
  }
  /* No values.  Use Null value for extremes */
  if (*a_max < *a_min) {
    if (a_DataType == DATASET_TYPE_FLOAT) {
      *a_min = *a_max = fNullValue;
    }
    else if (a_DataType == DATASET_TYPE_DOUBLE) {
      *a_min = *a_max = dNullValue;
    }
  }
  return 1;
} /* xfiDetermineTimestepMinsAndMaxs */
/******************************************************************************/
/*   FUNCTION  xfWriteScalarTimestep */
/*! PURPOSE:
 *   - NOTES:     Append a timestep to a scalar dataset. */
/******************************************************************************/
static xid xfiWriteScalarTimestep (xid a_Id, double a_Time, int a_NumValues,
                                     const void *a_Values, int a_DataType)
{
  double           Min, Max;

  if (a_NumValues < 1) {
    return ERROR_DATASET_NO_DATA;
  }

  xfiDetermineTimestepMinsAndMaxs (a_Id, a_NumValues, a_Values, a_DataType,
                                   &Min, &Max);

  return xfiWriteScalarTimestepMinMax(a_Id, a_Time, a_NumValues, a_Values,
                                      a_DataType, Min, Max);
} /* xfiWriteScalarTimestep */
/******************************************************************************/
/*   FUNCTION  xfiInitializeScalarTimestep */
/*! PURPOSE:
 *   - NOTES:     Append a timestep to a scalar dataset but not write values
 *             The min and max should be set to the min or max of the timestep */
/******************************************************************************/
static xid xfiInitializeScalarTimestep (xid a_Id, double a_Time,int a_NumValues,
              int a_DataType, double a_Min, double a_Max, hsize_t *a_timestepId)
{
  xmbool           bExtendable=XFALSE;
  float            FillVal;
  xid              ValuesId, TimesId, MinsId, MaxsId;
  int              status = 1;
  herr_t           (*old_func)(void*);
  void            *old_client_data;
  XDatasetParams   Params;
  int              Compression;
  hid_t            hdfType = H5T_NATIVE_FLOAT;
  xid              dsetId;
  int              numTimes = 0;
  hsize_t         *Dims = NULL, *MaxDims = NULL;
  xid              DatasetId, DataspaceId;
  int              Rank;

  if (a_NumValues < 1) {
    return ERROR_DATASET_NO_DATA;
  }

  H5Eget_auto1(&old_func, &old_client_data);

  /* Turn off error handling */
  H5Eset_auto1(NULL, NULL);

  /* See if the time, values, mins, maxs, and activity datasets if they exist */
  ValuesId = H5Dopen1(a_Id, DATASET_DSET_VALUES);
  TimesId = H5Dopen1(a_Id, DATASET_DSET_TIMES);
  MinsId = H5Dopen1(a_Id, DATASET_DSET_MINS); 
  MaxsId = H5Dopen1(a_Id, DATASET_DSET_MAXS);

  /* Restore previous error handler */
  H5Eset_auto1(old_func, old_client_data);

    // Set the HDF data type
  if (a_DataType == DATASET_TYPE_DOUBLE) {
    hdfType = xfpGetDefaultDoubleType();
  }

  /* Get the compression level for the dataset */
  status = xfpReadAttributeInt(a_Id, DATASET_ATT_COMPRESSION, 1, &Compression);
  if (status < 0) {
    if (ValuesId > 0) {
      xfCloseGroup(ValuesId);
    }
    if (TimesId > 0) {
      H5Dclose(TimesId);
    }
    if (MinsId > 0) {
      H5Dclose(MinsId);
    }
    if (MaxsId > 0) {
      H5Dclose(MaxsId);
    }
    return status;
  }

  status = xftIncNumTimes( a_Id );
  if (status < 0) {
    return status;
  }

  *a_timestepId = NONE;
  status = xftGetNumTimes( a_Id, &numTimes );
  if (status >= 0) {
    *a_timestepId = numTimes;
  }

  /* Dataset time values */
  if (TimesId < 0) {
    /* Create the time dataset */
    xfpDsetParamsInit(&Params, 1, XTRUE, Compression);
    xfpDsetParamsSetSizes(&Params, 0, 1, H5S_UNLIMITED, 
                          DATASET_TIME_CHUNK_SIZE);
    status = xfpWriteDsetDouble(a_Id, DATASET_DSET_TIMES, &Params, &a_Time);
    if (status < 0) {
      xfpDsetParamsDestroy(&Params);
      return status;
    }
    if (*a_timestepId == NONE) {
      *a_timestepId = 1;
    }
  }
  else {
    H5Dclose(TimesId);
    /* Time dataset exists append to it */
    status = xftAppendDset1DDouble(a_Id, DATASET_DSET_TIMES, 1, &a_Time);

    DatasetId = H5Dopen1(a_Id, DATASET_DSET_TIMES);
    if (DatasetId < 0) {
      return DatasetId;
    }
    DataspaceId = H5Dget_space(DatasetId);
    if (DataspaceId < 0) {
      H5Dclose(DatasetId);
      return DataspaceId;
    }
    /* Dataspace must be simple and have a rank of 1 */
    if (H5Sis_simple(DataspaceId) <= 0) {
      H5Dclose(DatasetId);
      H5Sclose(DataspaceId);
      return ERROR_DATASET_INVALID;
    }
    status = xftGetSimpleDataspaceInfoFromName(a_Id, DATASET_DSET_TIMES, !0, 
                                               &Rank, &Dims, &MaxDims);
    H5Dclose(DatasetId);
    H5Sclose(DataspaceId);
    if (status < 0) {
      return status;
    }
    if (Rank != 1) {
      return -1;
    }
    if (*a_timestepId == NONE) {
      *a_timestepId = Dims[0];
    }
    free(Dims);
    free(MaxDims);
  }

  /* Dataset values */
  if (ValuesId < 0) {
    if (xfpReadAttributeFloat(a_Id, DATASET_ATT_EXTENDFILL, 1, &FillVal) >= 0) {
      bExtendable = XTRUE;
    }
    
    /* Create the values dataset */
    xfpDsetParamsReset(&Params, 2, XTRUE, Compression);
    /* Time dimension */
    xfpDsetParamsSetSizes(&Params, 0, 1, H5S_UNLIMITED, 1);
    /* value dimension */
    if (bExtendable) {
      xfpDsetParamsSetSizes(&Params, 1, a_NumValues, H5S_UNLIMITED,
         Xmin(50, a_NumValues));
      if (a_DataType == DATASET_TYPE_FLOAT) {
        xfpDsetParamsUseFillValue(&Params, H5T_NATIVE_FLOAT, &FillVal);
      }
      else if (a_DataType == DATASET_TYPE_DOUBLE) { 
        xfpDsetParamsUseFillValue(&Params, H5T_NATIVE_DOUBLE, &FillVal);
      }
    }
    else {
      xfpDsetParamsSetSizes(&Params, 1, a_NumValues, a_NumValues, a_NumValues);
    }
    status = xfpCreateDset(a_Id, DATASET_DSET_VALUES, &Params, hdfType,
                           &dsetId);
    xfCloseGroup(dsetId);

    /* write out he data type because it maybe useful to know when reading */
    xfpWriteAttributeInt(a_Id, DATASET_TYPE, 1, &a_DataType);

    if (status < 0) {
      xfpDsetParamsDestroy(&Params);
      return status;
    }
  }
  else {
    H5Dclose(ValuesId);

    status = xftExtendDset2DFirstDim (a_Id, DATASET_DSET_VALUES, 
                 1, a_NumValues, hdfType);
    if (status < 0) {
      return status;
    }
  }

  if (MinsId < 0 && MaxsId < 0) {
    /* create the min and max datasets */
    xfpDsetParamsReset(&Params, 1, XTRUE, Compression);
    xfpDsetParamsSetSizes(&Params, 0, 1, H5S_UNLIMITED, 
                          DATASET_TIME_CHUNK_SIZE);
    
    if (a_DataType == DATASET_TYPE_FLOAT) {
      /* write min float */
      float fMin = (float)(a_Min);
      status = xfpWriteDsetFloat(a_Id, DATASET_DSET_MINS, &Params, &fMin);
      if (status >= 0) {
        /* write max float */
        float fMax = (float)(a_Max);
        status = xfpWriteDsetFloat(a_Id, DATASET_DSET_MAXS, &Params, &fMax);
      }
    }
    else if (a_DataType == DATASET_TYPE_DOUBLE) {
      /* write min double */
      status = xfpWriteDsetDouble(a_Id, DATASET_DSET_MINS, &Params, &a_Min);
      if (status >= 0) {
        /* write max double */
        status = xfpWriteDsetDouble(a_Id, DATASET_DSET_MAXS, &Params, &a_Max);
      }
    }
   
    xfpDsetParamsDestroy(&Params);
  }
  else if (MinsId > 0 && MaxsId > 0) {
    H5Dclose(MinsId);
    H5Dclose(MaxsId);
    /* The min and max datasets exist, append */
    if (a_DataType == DATASET_TYPE_FLOAT) {
      /* write min float */
       float fMin = (float)(a_Min);
      status = xftAppendDset1DFloat(a_Id, DATASET_DSET_MINS, 1, &fMin);
      if (status >= 0) {
        /* write max float */
        float fMax = (float)(a_Max);
        status = xftAppendDset1DFloat(a_Id, DATASET_DSET_MAXS, 1, &fMax);
      }
    }
    else if (a_DataType == DATASET_TYPE_DOUBLE) {
      /* write min double */
      status = xftAppendDset1DDouble(a_Id, DATASET_DSET_MINS, 1, &a_Min);
      if (status >= 0) {
        /* write max double */
        status = xftAppendDset1DDouble(a_Id, DATASET_DSET_MAXS, 1, &a_Max);
      }
    }
  }
  else {
    /* something is wrong.  Either both datasets should exist or neither */
    return ERROR_OTHER;
  }

  return status;
} /* xfiInitializeScalarTimestep */
/******************************************************************************/
/*   FUNCTION  xfiWriteScalarTimestepPortion */
/*! PURPOSE:
 *   - NOTES:     Append a timestep to a scalar dataset but not write values
 *             The min and max should be set to the min or max of the timestep */
/******************************************************************************/
static xid xfiWriteScalarTimestepPortion (xid a_id, hsize_t a_timeId,
              int a_NumValues, int a_startIndex, int a_DataType,
              const void *a_values)
{
  hid_t  hdfType;
  int    status;

  if (a_DataType == DATASET_TYPE_FLOAT) {
    hdfType = xfpGetDefaultFloatType();
  }
  else if (a_DataType == DATASET_TYPE_DOUBLE) {
    hdfType = xfpGetDefaultDoubleType();
  } 

  status = xfpWriteDset2DPortion(a_id, DATASET_DSET_VALUES, hdfType,
              a_timeId - 1, 1, a_startIndex - 1, a_NumValues, a_values);

  return status;
} /* xfiWriteScalarTimestepPortion */
/******************************************************************************/
/*   FUNCTION  xfiSetDatasetTimestepMinMax */
/*! PURPOSE:   overwrite min/max values for a specific timestep
 *   - NOTES:     A negative error is returned if the write fails
 *                timestepId is a 1 based index into the dataset and is returned
 *                from the Initialize function */
/******************************************************************************/
static xid xfiSetDatasetTimestepMinMax(xid xDatasetId, int timestepId,
                         int a_DataType, double minvalue, double maxvalue)
{ 
  float minFloat = (float)minvalue;
  float maxFloat = (float)maxvalue;
  int status = 0;

  if (a_DataType == DATASET_TYPE_FLOAT) {
    status = xfpWriteDset1DFloatPortion(xDatasetId, DATASET_DSET_MINS,
                     timestepId - 1, 1,  &minFloat);
    if (status < 0) {
      return status;
    }
    status = xfpWriteDset1DFloatPortion(xDatasetId, DATASET_DSET_MAXS,
                     timestepId - 1, 1, &maxFloat);
    if (status < 0) {
      return status;
    }
  }
  else if (a_DataType == DATASET_TYPE_DOUBLE) {
    status = xfpWriteDset1DDoublePortion(xDatasetId, DATASET_DSET_MINS,
                     timestepId - 1, 1, &minvalue);
    if (status < 0) {
      return status;
    }
    status = xfpWriteDset1DDoublePortion(xDatasetId, DATASET_DSET_MAXS,
                     timestepId - 1, 1, &maxvalue);
    if (status < 0) {
      return status;
    }
  }
  return status;
} /* xfiSetDatasetTimestepMinMax */
/******************************************************************************/
/*   FUNCTION  xfiWriteScalarTimestepMinMax */
/*! PURPOSE:
 *   - NOTES:     Append a timestep to a scalar dataset. 
 *             The min and max should be set to the min or max of the timestep */
/******************************************************************************/
static xid xfiWriteScalarTimestepMinMax (xid a_Id, double a_Time, int a_NumValues,
                                           const void *a_Values, int a_DataType,
                                           double a_Min, double a_Max)
{
  int              status = 1;
  hsize_t          timeId;
 
  status = xfiInitializeScalarTimestep(a_Id, a_Time, a_NumValues, a_DataType,
              a_Min, a_Max, &timeId);
  if (status < 0) {
    return status;
  }

  /* Dataset values */
    // We want to write the dataset values, we know the id of the timestep
  status = xfiWriteScalarTimestepPortion(a_Id, timeId, a_NumValues, 1,
              a_DataType, a_Values);
  return status;
} /* xfiWriteScalarTimestepMinMax */
/******************************************************************************/
/*   FUNCTION  xfWriteVectorTimestep */
/*! PURPOSE:
 *   - NOTES:     Append a timestep to a vector dataset. */
/******************************************************************************/
XMDF_API xid xfWriteVectorTimestep(xid a_Id, double a_Time, int 
                    a_NumValues, int a_NumComponents, const float *a_Values)
{
  return xfiWriteVectorTimestep(a_Id, a_Time, a_NumValues, a_NumComponents,
                         (void*)a_Values, DATASET_TYPE_FLOAT);
} /* xfWriteVectorTimestep */
/******************************************************************************/
/*   FUNCTION  xfWriteVectorTimestepMinMax */
/*! PURPOSE:
 *   - NOTES:     Append a timestep to a vector dataset. 
 *             The min and max should be set to the min or max (magnitude)
 *             of the timestep */
/******************************************************************************/
XMDF_API xid xfWriteVectorTimestepMinMax(xid a_Id, double a_Time, int 
                    a_NumValues, int a_NumComponents, const float *a_Values,
                    float a_Min, float a_Max)
{
  return xfiWriteVectorTimestepMinMax(a_Id, a_Time, a_NumValues, a_NumComponents,
                               (void*)a_Values, DATASET_TYPE_FLOAT,
                               (double)a_Min, (double)a_Max);
} /* xfWriteVectorTimestepMinMax */
/******************************************************************************/
/*   FUNCTION  xfWriteVectorTimestepFloat */
/*! PURPOSE:
 *   - NOTES:     Append a timestep to a float Vector dataset. */
/******************************************************************************/
XMDF_API xid xfWriteVectorTimestepFloat(xid a_Id, double a_Time,
                                        int a_NumValues, int a_NumComponents,
                                        const float *a_Values)
{
  
  return xfiWriteVectorTimestep(a_Id, a_Time, a_NumValues, a_NumComponents,
                         (void*)a_Values, DATASET_TYPE_FLOAT);

} /* xfWriteVectorTimestepFloat */
/******************************************************************************/
/*   FUNCTION  xfWriteVectorTimestepMinMaxFloat */
/*! PURPOSE:
 *   - NOTES:     Append a timestep to a float Vector dataset. 
 *             The min and max should be set to the min or max of the timestep */
/******************************************************************************/
XMDF_API xid xfWriteVectorTimestepMinMaxFloat (xid a_Id, double a_Time,
                                               int a_NumValues, int a_NumComponents,
                                               const float *a_Values,
                                               double a_Min, double a_Max)
{
  
  return xfiWriteVectorTimestepMinMax(a_Id, a_Time, a_NumValues, a_NumComponents,
                               (void*)a_Values, DATASET_TYPE_FLOAT,
                               a_Min, a_Max);

} /* xfWriteVectorTimestepMinMaxFloat */
/******************************************************************************/
/*   FUNCTION  xfWriteVectorTimestepDouble */
/*! PURPOSE:
 *   - NOTES:     Append a timestep to a double Vector dataset. */
/******************************************************************************/
XMDF_API xid xfWriteVectorTimestepDouble(xid a_Id, double a_Time, int a_NumValues,
                                         int a_NumComponents, const double *a_Values)
{
  
  return xfiWriteVectorTimestep(a_Id, a_Time, a_NumValues, a_NumComponents,
                         (void*)a_Values, DATASET_TYPE_DOUBLE);

} /* xfWriteVectorTimestepDouble */
/******************************************************************************/
/*   FUNCTION  xfWriteVectorTimestepMinMaxDouble */
/*! PURPOSE:
 *   - NOTES:     Append a timestep to a double Vector dataset. 
 *             The min and max should be set to the min or max of the timestep */
/******************************************************************************/
XMDF_API xid xfWriteVectorTimestepMinMaxDouble (xid a_Id, double a_Time, 
                                                int  a_NumValues, int a_NumComponents,
                                                const double *a_Values,
                                                double a_Min, double a_Max)
{
  return xfiWriteVectorTimestepMinMax(a_Id, a_Time, a_NumValues, a_NumComponents,
                               (void*)a_Values, DATASET_TYPE_DOUBLE,
                               a_Min, a_Max);

} /* xfWriteVectorTimestepMinMaxDouble */
/******************************************************************************/
/*   FUNCTION  xfInitializeVectorTimestep */
/*! PURPOSE:
 *   - NOTES:  Initialize and setup the datasets for a new vector timestep
 *             The min and max should be set to the min or max of the timestep */
/******************************************************************************/
XMDF_API xid xfInitializeVectorTimestep (xid a_Id, double a_Time, 
               int  a_NumValues, int a_NumComponents, float a_Min, float a_Max,
               hsize_t *a_timestepId)
{
  return xfiInitializeVectorTimestep(a_Id, a_Time, a_NumValues, a_NumComponents,
                                     DATASET_TYPE_FLOAT, (double)a_Min,
                                     (double)a_Max,
                                     a_timestepId);
} /* xfInitializeVectorTimestep */
/******************************************************************************/
/*   FUNCTION  xfWriteVectorTimestepPortion */
/*! PURPOSE:
 *   - NOTES:     Append a timestep to a double Vector dataset. 
 *             The min and max should be set to the min or max of the timestep */
/******************************************************************************/
XMDF_API xid xfWriteVectorTimestepPortion (xid a_Id, int timestepId,
                    int numValuesToWrite, int nComponentsToWrite,
                    int startIndex, int startComponent, const float *a_values)
{
  return xfiWriteVectorTimestepPortion(a_Id, timestepId, numValuesToWrite, 
            nComponentsToWrite, startIndex, startComponent, DATASET_TYPE_FLOAT,
            (void*)a_values);
} /* xfWriteVectorTimestepPortion */

/******************************************************************************/
/*   FUNCTION  xfWriteVectorTimestep */
/*! PURPOSE:
 *   - NOTES:     Append a timestep to a vector dataset. */
/******************************************************************************/
static xid xfiWriteVectorTimestep(xid a_Id, double a_Time, 
                                    int a_NumValues, int a_NumComponents, 
                                    const void *a_Values, int a_DataType)
{
  xid              propGroupId;
  double           Min, Max;
  xmbool           bUseNullValue = XFALSE;
  float            fNullValue = 0.0;
  double           dNullValue = 0.0;
  herr_t         (*old_func)(void*);
  void            *old_client_data;
  int              iValue, iComp, status;
  double           mag_sqr;


  status = -1;

  if (a_NumValues < 1) {
    return ERROR_DATASET_NO_DATA;
  }

  H5Eget_auto1(&old_func, &old_client_data);

  /* Turn off error handling */
  H5Eset_auto1(NULL, NULL);

    /* get the null value if there is one */
  if (xfOpenPropertyGroup(a_Id, &propGroupId) >= 0) {
    if (a_DataType == DATASET_TYPE_FLOAT) {
      status = xfReadPropertyFloat(propGroupId, PROP_NULL_VALUE, 1, &fNullValue);
    }
    else if (a_DataType == DATASET_TYPE_DOUBLE) {
      status = xfReadPropertyDouble(propGroupId, PROP_NULL_VALUE, 1, &dNullValue);
    }
    if (status >= 0) {
      bUseNullValue = XTRUE;
    }
    xfCloseGroup(propGroupId);
  }

  /* Restore previous error handler */
  H5Eset_auto1(old_func, old_client_data);

  Min = MAX_DOUBLE;
  Max = MIN_DOUBLE;
  for (iValue = 0; iValue < a_NumValues; iValue++) {      
    mag_sqr = 0.0;
    for (iComp = 0; iComp < a_NumComponents; iComp++) {
      if (a_DataType == DATASET_TYPE_FLOAT) {
        if (bUseNullValue == XFALSE || 
            ((float*)a_Values)[iValue*a_NumComponents + iComp] != fNullValue) {
          mag_sqr += pow(((float*)a_Values)[iValue*a_NumComponents + iComp], 2);
        }
      }
      else if (a_DataType == DATASET_TYPE_DOUBLE) {
        if (bUseNullValue == XFALSE || 
            ((double*)a_Values)[iValue*a_NumComponents + iComp] != dNullValue) {
          mag_sqr += pow(((double*)a_Values)[iValue*a_NumComponents + iComp], 2);
        }
      }
    }
    Min = Xmin(Min, mag_sqr);
    Max = Xmax(Max, mag_sqr);
  }
  Min = (double)pow(Min, (double)1.0/(double)a_NumComponents);
  Max = (double)pow(Max, (double)1.0/(double)a_NumComponents);

  return xfiWriteVectorTimestepMinMax(a_Id, a_Time, a_NumValues, a_NumComponents,
                                     a_Values, a_DataType, Min, Max);

} /* xfiWriteVectorTimestep */
/******************************************************************************/
/*   FUNCTION  xfiWriteVectorTimestepMinMax */
/*! PURPOSE:
 *   - NOTES:     Append a timestep to a vector dataset. 
 *             The min and max should be set to the min or max (magnitude)
 *             of the timestep */
/******************************************************************************/
static xid xfiWriteVectorTimestepMinMax (xid a_Id, double a_Time, 
                                          int a_NumValues, int a_NumComponents,
                                          const void *a_Values, int a_DataType,
                                          double a_Min, double a_Max)
{
  int              status = 1;
  hsize_t timeId;

  status = xfiInitializeVectorTimestep(a_Id, a_Time, a_NumValues,
              a_NumComponents, a_DataType, a_Min, a_Max,
              &timeId);
  if (status < 0) {
    return status;
  }

  status = xfiWriteVectorTimestepPortion(a_Id, timeId, a_NumValues, 
              a_NumComponents, 1, 1, a_DataType, a_Values);

  return status;
} /* xfiWriteVectorTimestepMinMax */
/******************************************************************************/
/*   FUNCTION  xfiInitializeVectorTimestep */
/*! PURPOSE:
 *   - NOTES:     Append a timestep to a vector dataset. 
 *             The min and max should be set to the min or max (magnitude)
 *             of the timestep */
/******************************************************************************/
static xid xfiInitializeVectorTimestep (xid a_Id, double a_Time, 
              int a_NumValues, int a_NumComponents, int a_DataType,
              double a_Min, double a_Max, hsize_t *a_timestepId)
{
  xid              ValuesId, TimesId, MinsId, MaxsId;
  int              status = 1;
  herr_t         (*old_func)(void*);
  void            *old_client_data;
  XDatasetParams   Params;
  int              Compression;
  hid_t            hdfType = H5T_NATIVE_FLOAT;
  xid              dsetId;
  int              numTimes = 0;
  hsize_t         *Dims = NULL, *MaxDims = NULL;
  xid              DatasetId, DataspaceId;
  int              Rank;

  if (a_NumValues < 1) {
    return ERROR_DATASET_NO_DATA;
  }

  H5Eget_auto1(&old_func, &old_client_data);

  /* Turn off error handling */
  H5Eset_auto1(NULL, NULL);

  /* Delete the time, values, mins, and maxs, activity datasets if they exist */
  ValuesId = H5Dopen1(a_Id, DATASET_DSET_VALUES);
  TimesId = H5Dopen1(a_Id, DATASET_DSET_TIMES);
  MinsId = H5Dopen1(a_Id, DATASET_DSET_MINS); 
  MaxsId = H5Dopen1(a_Id, DATASET_DSET_MAXS);

  /* Restore previous error handler */
  H5Eset_auto1(old_func, old_client_data);

    // Set the HDF data type
  if (a_DataType == DATASET_TYPE_DOUBLE) {
    hdfType = xfpGetDefaultDoubleType();
  }

  /* Get the compression level for the dataset */
  status = xfpReadAttributeInt(a_Id, DATASET_ATT_COMPRESSION, 1, &Compression);
  if (status < 0) {
    if (ValuesId > 0) {
      H5Dclose(ValuesId);
    }
    if (TimesId > 0) {
      H5Dclose(TimesId);
    }
    if (MinsId > 0) {
      H5Dclose(MinsId);
    }
    if (MaxsId > 0) {
      H5Dclose(MaxsId);
    }
    return status;
  }

  status = xftIncNumTimes( a_Id );
  if (status < 0) {
    return status;
  }

  *a_timestepId = NONE;
  status = xftGetNumTimes(a_Id, &numTimes);
  if (status >= 0) {
    *a_timestepId = numTimes;
  }

  /* Dataset time values */
  if (TimesId < 0) {
    /* Create the time dataset */
    xfpDsetParamsInit(&Params, 1, XTRUE, Compression);
    xfpDsetParamsSetSizes(&Params, 0, 1, H5S_UNLIMITED, 
                          DATASET_TIME_CHUNK_SIZE);
    status = xfpWriteDsetDouble(a_Id, DATASET_DSET_TIMES, &Params, &a_Time);
    if (status < 0) {
      xfpDsetParamsDestroy(&Params);
      return status;
    }
    if (*a_timestepId == NONE) {
      *a_timestepId = 1;
    }
  }
  else {
    H5Dclose(TimesId);
    /* Time dataset exists append to it */
    status = xftAppendDset1DDouble(a_Id, DATASET_DSET_TIMES, 1, &a_Time);

    DatasetId = H5Dopen1(a_Id, DATASET_DSET_TIMES);
    if (DatasetId < 0) {
      return DatasetId;
    }
    DataspaceId = H5Dget_space(DatasetId);
    if (DataspaceId < 0) {
      H5Dclose(DatasetId);
      return DataspaceId;
    }
    /* Dataspace must be simple and have a rank of 1 */
    if (H5Sis_simple(DataspaceId) <= 0) {
      H5Dclose(DatasetId);
      H5Sclose(DataspaceId);
      return ERROR_DATASET_INVALID;
    }
    status = xftGetSimpleDataspaceInfoFromName(a_Id, DATASET_DSET_TIMES, !0, 
                                               &Rank, &Dims, &MaxDims);
    H5Dclose(DatasetId);
    H5Sclose(DataspaceId);
    if (status < 0) {
      return status;
    }
    if (Rank != 1) {
      return -1;
    }
    if (*a_timestepId == NONE) {
      *a_timestepId = Dims[0];
    }
    free(Dims);
    free(MaxDims);
  }

  /* Dataset values */
  if (ValuesId < 0) {
    /* Create the values dataset */
    xfpDsetParamsReset(&Params, 3, XTRUE, Compression);
    /* Time dimension */
    xfpDsetParamsSetSizes(&Params, 0, 1, H5S_UNLIMITED, 1);
    /* value dimension */
    xfpDsetParamsSetSizes(&Params, 1, a_NumValues, a_NumValues, a_NumValues);
    /* Components dimension */
    xfpDsetParamsSetSizes(&Params, 2, a_NumComponents, a_NumComponents,
                          a_NumComponents);

    status = xfpCreateDset(a_Id, DATASET_DSET_VALUES, &Params, hdfType, &dsetId);
    xfCloseGroup(dsetId);

     /* write out he data type because it maybe useful to know when reading */
    xfpWriteAttributeInt(a_Id, DATASET_TYPE, 1, &a_DataType);

    if (status < 0) {
      xfpDsetParamsDestroy(&Params);
      return status;
    }
  }
  else {
    H5Dclose(ValuesId);

    /* values array exists append to it */
    status = xftExtendDset3DFirstDim(a_Id, DATASET_DSET_VALUES, 1, a_NumValues,
                                     a_NumComponents, hdfType);    
    if (status < 0) {
      return status;
    }
  }

  if (MinsId < 0 && MaxsId < 0) {
    /* create the min and max datasets */
    xfpDsetParamsReset(&Params, 1, XTRUE, Compression);
    xfpDsetParamsSetSizes(&Params, 0, 1, H5S_UNLIMITED, 
                          DATASET_TIME_CHUNK_SIZE);
    if (a_DataType == DATASET_TYPE_FLOAT) {
      /* write min float */
      float fMin = (float)(a_Min);
      status = xfpWriteDsetFloat(a_Id, DATASET_DSET_MINS, &Params, &fMin);
      if (status >= 0) {
        /* write max float */
        float fMax = (float)(a_Max);
        status = xfpWriteDsetFloat(a_Id, DATASET_DSET_MAXS, &Params, &fMax);
      }
    }
    else if (a_DataType == DATASET_TYPE_DOUBLE) {
      /* write min double */
      status = xfpWriteDsetDouble(a_Id, DATASET_DSET_MINS, &Params, &a_Min);
      if (status >= 0) {
        /* write max double */
        status = xfpWriteDsetDouble(a_Id, DATASET_DSET_MAXS, &Params, &a_Max);
      }
    }
    xfpDsetParamsDestroy(&Params);
  }
  else if (MinsId > 0 && MaxsId > 0) {
    H5Dclose(MinsId);
    H5Dclose(MaxsId);
    /* The min and max datasets exist, append */
    if (a_DataType == DATASET_TYPE_FLOAT) {
      /* write min float */
       float fMin = (float)(a_Min);
      status = xftAppendDset1DFloat(a_Id, DATASET_DSET_MINS, 1, &fMin);
      if (status >= 0) {
        /* write max float */
        float fMax = (float)(a_Max);
        status = xftAppendDset1DFloat(a_Id, DATASET_DSET_MAXS, 1, &fMax);
      }
    }
    else if (a_DataType == DATASET_TYPE_DOUBLE) {
      /* write min double */
      status = xftAppendDset1DDouble(a_Id, DATASET_DSET_MINS, 1, &a_Min);
      if (status >= 0) {
        /* write max double */
        status = xftAppendDset1DDouble(a_Id, DATASET_DSET_MAXS, 1, &a_Max);
      }
    }
  }
  else {
    /* something is wrong.  Either both datasets should exist or neither */
    return ERROR_OTHER;
  }

  return status;

} /* xfiInitializeVectorTimestep */
/******************************************************************************/
/*   FUNCTION  xfiWriteVectorTimestepPortion */
/*! PURPOSE:
 *   - NOTES:     Append a timestep to a scalar dataset but not write values
 *             The min and max should be set to the min or max of the timestep */
/******************************************************************************/
static xid xfiWriteVectorTimestepPortion (xid a_id, hsize_t a_timeId,
              int a_NumValuesToWrite, int a_NumComponentsToWrite, 
              int a_startIndex, int a_startComponent, int a_DataType,
              const void *a_values)
{
  hid_t  hdfType;
  int    status;

  if (a_DataType == DATASET_TYPE_FLOAT) {
    hdfType = xfpGetDefaultFloatType();
  }
  else if (a_DataType == DATASET_TYPE_DOUBLE) {
    hdfType = xfpGetDefaultDoubleType();
  } 

  status = xfpWriteDset3DPortion(a_id, DATASET_DSET_VALUES, hdfType,
              a_timeId - 1, 1, a_startIndex - 1, a_NumValuesToWrite,
              a_startComponent - 1, a_NumComponentsToWrite, a_values);

  return status;
} /* xfiWriteVectorTimestepPortion */

/******************************************************************************/
/*   FUNCTION  xfWriteActivityTimestep */
/*! PURPOSE:
 *   - NOTES:     Write dataset activity information. */
/******************************************************************************/
XMDF_API xid xfWriteActivityTimestep (xid a_Id, int a_NumActive, 
                                      const unsigned char *a_Active)
{
  int            status;
  hsize_t            timestepId;

  status = xfInitializeActivityTimestep(a_Id, a_NumActive, &timestepId);
  if (status < 0) {
    return status;
  }

  status = xfWriteActivityTimestepPortion(a_Id, timestepId, a_NumActive, 1,
                                          a_Active);

  return status;
} /* xfWriteActivityTimestep */
/******************************************************************************/
/*   FUNCTION  xfInitializeActivityTimestep */
/*! PURPOSE:
 *   - NOTES:     Write dataset activity information. */
/******************************************************************************/
XMDF_API xid xfInitializeActivityTimestep (xid a_Id, int a_NumActive,
                                           hsize_t *a_timestepId)
{
  xid            ActiveId;
  int            status;
  XDatasetParams Params;
  int      Compression;
  herr_t (*old_func)(void*);
  void    *old_client_data;
  hid_t          dsetId;
  hid_t DatasetId, DataspaceId;
  int              Rank;
  hsize_t         *Dims = NULL, *MaxDims = NULL;

  if (a_NumActive < 1) {
    return ERROR_DATASET_NO_DATA;
  }

  H5Eget_auto1(&old_func, &old_client_data);

  /* Turn off error handling */
  H5Eset_auto1(NULL, NULL);

  /* Open the activity information if exists */
  ActiveId = H5Dopen1(a_Id, DATASET_DSET_ACTIVE);

  /* Restore previous error handler */
  H5Eset_auto1(old_func, old_client_data);

  /* Get the compression level for the dataset */
  status = xfpReadAttributeInt(a_Id, DATASET_ATT_COMPRESSION, 1, &Compression);

  if (ActiveId < 0) {
    /* create the activity dataset */
    xfpDsetParamsInit(&Params, 2, XTRUE, Compression);
    /* Time dimension */
    xfpDsetParamsSetSizes(&Params, 0, 1, H5S_UNLIMITED, 1);
    /* value dimension */
    xfpDsetParamsSetSizes(&Params, 1, a_NumActive, a_NumActive, a_NumActive);

    status = xfpCreateDset(a_Id, DATASET_DSET_ACTIVE, &Params, 
                             H5T_NATIVE_UCHAR, &dsetId);
    xfCloseGroup(dsetId);                             

    xfpDsetParamsDestroy(&Params);

    *a_timestepId = 1;
  }
  else {
    H5Dclose(ActiveId);


    /* Activity datasets exists append */
    status = xftExtendDset2DFirstDim(a_Id, DATASET_DSET_ACTIVE, 1, a_NumActive,
                                     H5T_NATIVE_UCHAR);
    if (status >= 0) {
      DatasetId = H5Dopen1(a_Id, DATASET_DSET_ACTIVE);
      if (DatasetId < 0) {
        return DatasetId;
      }
      DataspaceId = H5Dget_space(DatasetId);
      if (DataspaceId < 0) {
        H5Dclose(DatasetId);
        return DataspaceId;
      }
      /* Dataspace must be simple and have a rank of 1 */
      if (H5Sis_simple(DataspaceId) <= 0) {
        H5Dclose(DatasetId);
        H5Sclose(DataspaceId);
        return ERROR_DATASET_INVALID;
      }
      status = xftGetSimpleDataspaceInfoFromName(a_Id, DATASET_DSET_ACTIVE, 0, 
                                                 &Rank, &Dims, &MaxDims);
      H5Dclose(DatasetId);
      H5Sclose(DataspaceId);
      if (status < 0) {
        return status;
      }
      if (Rank != 2) {
        return -1;
      }
      *a_timestepId = Dims[0];
      free(Dims);
      free(MaxDims);
    }
  }
  return status;
} /* xfInitializeActivityTimestep */
/******************************************************************************/
/*   FUNCTION  xfWriteActivityTimestepPortion */
/*! PURPOSE:
 *   - NOTES:     Write dataset activity information. */
/******************************************************************************/
XMDF_API xid xfWriteActivityTimestepPortion (xid a_Id, hsize_t a_timestepId,
               int a_NumValuesToWrite, int a_startIndex, 
               const unsigned char *a_activityValues)
{
  int            status;

  status = xfpWriteDset2DPortion(a_Id, DATASET_DSET_ACTIVE, H5T_NATIVE_UCHAR,
              a_timestepId - 1, 1, a_startIndex - 1, a_NumValuesToWrite,
              a_activityValues);
  
  return status;
} /* xfWriteActivityTimestepPortion */

/******************************************************************************/
/*   FUNCTION  xfDatasetReftime */
/*! PURPOSE:   Set a reference time for a dataset
 *   - NOTES:     */
/******************************************************************************/
XMDF_API xid xfDatasetReftime (xid a_Id, double a_Reftime)
{
  int    status;

  status = xfpWriteAttributeDouble(a_Id, DATASET_REFTIME, 1, &a_Reftime);

  return status;
} /* xfDatasetReftime */
/******************************************************************************/
/*   FUNCTION  xfUseDatasetReftime */
/*! PURPOSE:   See if a reference time exists for a dataset
 *   - NOTES:     */
/******************************************************************************/
XMDF_API xid xfUseDatasetReftime (xid a_Id, xmbool *a_bUseReftime)
{
  int    status;

  status = xfpDoesAttributeExist(a_Id, DATASET_REFTIME, a_bUseReftime);

  return status;
} /* xfDatasetReftime */
/******************************************************************************/
/*   FUNCTION  xfReadDatasetReftime */
/*! PURPOSE:   Read a reference time for a dataset
 *   - NOTES:     */
/******************************************************************************/
XMDF_API xid xfReadDatasetReftime (xid a_Id, double *a_dReftime)
{
  int    status;

  status = xfpReadAttributeDouble(a_Id, DATASET_REFTIME, 1, a_dReftime);

  return status;
} /* xfReadDatasetReftime */
/******************************************************************************/
/*   FUNCTION  xfGetScalarDatasetsInfo */
/*! PURPOSE:
 *   - NOTES:     Get the number and max path length of scalar datasets in the path. */
/******************************************************************************/
XMDF_API xid xfGetScalarDatasetsInfo (xid a_Id, int *a_Number,
                                         int *a_MaxPathLength)
{
  int    status;

  status = xfpNumGroupsOfType(a_Id, GROUP_TYPE_DATASET_SCALAR, a_Number, 
                              a_MaxPathLength, XFALSE);

  return status;
} /* xfGetScalarDatasetsInfo */
/******************************************************************************/
/*   FUNCTION  xfGetScalarDatasetPaths */
/*! PURPOSE:   Get the paths to scalar datasets under a starting group
 *   - NOTES:     The Path array must already be allocated to a size Number by Maxlength. */
/******************************************************************************/
XMDF_API xid xfGetScalarDatasetPaths (xid a_Id, int a_Number, 
                                         int a_MaxPathLength, char *a_Paths)
{
  int    status;

  status = xfpPathGroupsOfType(a_Id, GROUP_TYPE_DATASET_SCALAR, a_Number, 
                              a_MaxPathLength, a_Paths, XFALSE);

  return status;
} /* xfGetScalarDatasetPaths */
/******************************************************************************/
/*   FUNCTION  xfGetStationInfo */
/*! PURPOSE:
 *   - NOTES:     Get the number and max path length of station datasets in the path. */
/******************************************************************************/
XMDF_API xid xfGetStationInfo (xid a_Id, int *a_Number, 
                                         int *a_MaxPathLength)
{
  int    status;

  status = xfpNumGroupsOfType(a_Id, GROUP_TYPE_STATION, a_Number, 
                              a_MaxPathLength, XFALSE);

  return status;
} /* xfGetStationInfo */
/******************************************************************************/
/*   FUNCTION  xfGetStationPaths */
/*! PURPOSE:   Get the paths to station datasets (3D Grid) under a starting group
 *   - NOTES:     The Path array must already be allocated to a size Number by Maxlength. */
/******************************************************************************/
XMDF_API xid xfGetStationPaths (xid a_Id, int a_Number, 
                                         int a_MaxPathLength, char *a_Paths)
{
  int    status;

  status = xfpPathGroupsOfType(a_Id, GROUP_TYPE_STATION, a_Number, 
                              a_MaxPathLength, a_Paths, XFALSE);

  return status;
} /* xfGetStationPaths */
/******************************************************************************/
/*   FUNCTION  xfGetVectorDatasetsInfo */
/*! PURPOSE:
 *   - NOTES:     Get the number and max path length of scalar datasets in the path. */
/******************************************************************************/
XMDF_API xid xfGetVectorDatasetsInfo (xid a_Id, int *a_Number, 
                                         int *a_MaxPathLength)
{
  int    status;

  status = xfpNumGroupsOfType(a_Id, GROUP_TYPE_DATASET_VECTOR, a_Number, 
                              a_MaxPathLength, XFALSE);

  return status;
} /* xfGetVectorDatasetsInfo */
/******************************************************************************/
/*   FUNCTION  xfGetVectorDatasetPaths */
/*! PURPOSE:   Get the paths to scalar datasets under a starting group
 *   - NOTES:     The Path array must already be allocated to a size Number by Maxlength. */
/******************************************************************************/
XMDF_API xid xfGetVectorDatasetPaths (xid a_Id, int a_Number, 
                                         int a_MaxPathLength, char *a_Paths)
{
  int    status;

  status = xfpPathGroupsOfType(a_Id, GROUP_TYPE_DATASET_VECTOR, a_Number, 
                              a_MaxPathLength, a_Paths, XFALSE);

  return status;
} /* xfGetVectorDatasetPaths */
/******************************************************************************/
/*   FUNCTION  xfGetScalarDatasetGroupId */
/*! PURPOSE:
 *   - NOTES:     Open the Scalar datasets for a Mesh or Grid. */
/******************************************************************************/
XMDF_API xid xfGetScalarDatasetGroupId (xid a_Id)
{
/*   int    status; */
  xid    Id;

  /*status =*/ xfOpenGroup(a_Id, GROUP_TYPE_DATASET_SCALAR, &Id);

  return Id;
} /* xfGetScalarDatasetGroupId */
/******************************************************************************/
/*   FUNCTION  xfGetVectorDatasetGroupId */
/*! PURPOSE:
 *   - NOTES:     Open the Scalar datasets for a Mesh or Grid. */
/******************************************************************************/
XMDF_API xid xfGetVectorDatasetGroupId (xid a_Id)
{
/*   int    status; */
  xid    Id;

  /*status =*/ xfOpenGroup(a_Id, GROUP_TYPE_DATASET_VECTOR, &Id);

  return Id;
} /* xfGetVectorDatasetGroupId */
/*!
 * \ingroup files
 * xfChangeScalarValuesTimestepFloat()
 * \brief Changes float dataset values at one or more indices in a timestep
 * \param a_Id The id of the scalar dataset
 * \param a_TimestepIndex One based timestep index to modify
 * \param a_NumValsToEdit The number of values to change
 * \param a_Indices A one-dimensional array of indices to modify must be of size a_NumValsToEdit
 * \param a_NewValues A one-dimensional array of values must be of size a_NumValsToEdit
 * \return A negative value upon failure.
 *
 * FORTRAN: SUBROUTINE XF_CHANGE_SCALAR_VALUES_TIMESTEP_FLOAT(Name, Overwrite, FileId, Error)
 * - INTEGER(XID), INTENT(OUT) ::   FileId
 * - INTEGER(8), INTENT(IN) :: TimestepIndex
 * - INTEGER(8), INTENT(IN)       ::   NumValsToEdit
 * - INTEGER(8), DIMENSION(*), INTENT(IN) ::   Indices
 * - INTEGER(8), DIMENSION(*), INTENT(INT) :: NewValues
 * - INTEGER(8), INTENT(OUT)      ::   Error
 **********************************************************************/
XMDF_API xid xfChangeScalarValuesTimestepFloat (xid a_Id, int a_TimestepIndex,
               int a_NumValsToEdit, int *a_Indices, float *a_NewValues)
{
  int status = 0;
  int numValues = 0;
  int numTimesteps = 0;
  float *values = NULL;
  double myMin = 0.0, myMax = 0.0;
  //int flushStatus = 0;
  hid_t  MinsId, MaxsId;
  float myMinFloat = 0.0, myMaxFloat = 0.0;

  status = xfiReadWriteScalarValuesFloatAtIndices (a_Id, RW_WRITE, 
            a_NumValsToEdit, a_Indices, a_TimestepIndex, 1, a_NewValues);
  //flushStatus = H5Fflush(H5Iget_file_id(a_Id), H5F_SCOPE_GLOBAL);
  // this caused problems - don't do
  // H5Fflush(H5Iget_file_id(a_Id), H5F_SCOPE_GLOBAL);
  if (status >= 0) {
    // compute new minimum and maximum values
    status = xfGetDatasetNumTimes(a_Id, &numTimesteps);
    if (status < 0) {
      return status;
    }
    status = xfGetDatasetNumVals(a_Id, &numValues);
    if (status < 0) {
      return status;
    }
    values = (float*)malloc(numValues*sizeof(float));
    status = xfReadScalarValuesTimestep(a_Id, a_TimestepIndex, numValues, 
                                        values);
    if (status < 0) {
      free (values);
      return status;
    }
    status = xfiDetermineTimestepMinsAndMaxs(a_Id, numValues, (void*)values, 
                DATASET_TYPE_FLOAT, &myMin, &myMax);
    free (values);
    if (status < 0) {
      return status;
    }
      // correct the minimum and maximum values
    MinsId = H5Dopen1(a_Id, DATASET_DSET_MINS); 
    MaxsId = H5Dopen1(a_Id, DATASET_DSET_MAXS);
    if (MinsId > 0 && MaxsId > 0) {
      myMinFloat = (float)myMin;
      myMaxFloat = (float)myMax;
      xfpWriteDset1DFloatPortion(a_Id, DATASET_DSET_MINS, a_TimestepIndex - 1,
                                 1, &myMinFloat);
      xfpWriteDset1DFloatPortion(a_Id, DATASET_DSET_MAXS, a_TimestepIndex - 1,
                                 1, &myMaxFloat);
    }
    if (MinsId > 0) {
      H5Dclose(MinsId);
    }
    if (MaxsId > 0) {
      H5Dclose(MaxsId);
    }

  }
  return status;
} /* xfChangeScalarValuesTimestep */
/******************************************************************************/
/*   FUNCTION  xfGetDatasetNumTimes */
/*! PURPOSE:
 *   - NOTES:     */
/******************************************************************************/
XMDF_API xid xfGetDatasetNumTimes (xid a_Id, int *a_Numtimes)
{
  int     status = XTRUE;
  xid     TimesId, SpaceId;
  htri_t  Simple;
  int     Rank;
  hsize_t *Dims = NULL, *Maxdims = NULL;

  /* make sure group is a dataset group */
  if (!xfpIsGroupOfType(a_Id, GROUP_TYPE_DATASET_SCALAR) &&
      !xfpIsGroupOfType(a_Id, GROUP_TYPE_DATASET_VECTOR)) {
    return ERROR_DATASET_INVALID;
  }

  /* open the timestep dataset and see how big it is */
  TimesId = H5Dopen1(a_Id, DATASET_DSET_TIMES);
  if (TimesId < 0) {
    return ERROR_DATASET_INVALID;
  }

  /* Get the dataspace and make sure it is a simple dataset of rank 1 */
  SpaceId = H5Dget_space(TimesId);
  if (SpaceId < 0) {
    H5Dclose(TimesId);
    return ERROR_DATASET_INVALID;
  }

  Simple = H5Sis_simple(SpaceId);
  if (!Simple) {
    H5Sclose(SpaceId);
    H5Dclose(TimesId);
    return ERROR_DATASET_INVALID;
  }

  status = xfpGetSimpleDataspaceInfo(SpaceId, &Rank, &Dims, &Maxdims);
  if (status < 0 || Rank != 1) {
    if (Dims) {
      free(Dims);
    }
    if (Maxdims) {
      free(Maxdims);
    }
    H5Sclose(SpaceId);
    H5Dclose(TimesId);
    return ERROR_DATASET_INVALID;
  }

  /* attempt to get NumTimes */
  status = xftGetNumTimes(a_Id, a_Numtimes);
  if (status < 0)
  {
    H5Sclose(SpaceId);
    H5Dclose(TimesId);
    free(Dims);
    free(Maxdims);
    return status;
  }

  /* if NumTimes is not used then use end of dataset */
  if (*a_Numtimes < 0)
    {
      *a_Numtimes = (int)Dims[0];
    }

  H5Sclose(SpaceId);
  H5Dclose(TimesId);
  free(Dims);
  free(Maxdims);
  return status;
} /* xfGetDatasetNumTimes */
/******************************************************************************/
/*   FUNCTION  xfSetDatasetNumTimes */
/*! PURPOSE:   Reduce the number of timesteps in the dataset
 *   - NOTES:     A negative error is returned if NumTimes is greater than
 *    the number of time steps in the dataset.
 *    If NumTimes does not exist in the file or is negative, then
 *    all stored time steps are used by default. */
/******************************************************************************/
XMDF_API xid xfSetDatasetNumTimes(xid a_Id, int a_NumTimes)
{
  int        status;
  xid        DatasetId, DataspaceId;
  int        numTimes;

  status = xfGetDatasetNumTimes( a_Id, &numTimes );
  if (status < 0 || numTimes < 0) /* no NumTimes so compare Dataset size */
  {
    DatasetId = H5Dopen1(a_Id, DATASET_DSET_VALUES);
    if (DatasetId < 0) {
     return DatasetId;
    }

    DataspaceId = H5Dget_space(DatasetId);
    if (DataspaceId < 0) {
      H5Dclose(DatasetId);
      return DataspaceId;
    }

      /* Dataspace must be simple */
    if (H5Sis_simple(DataspaceId) <= 0) {
      H5Dclose(DatasetId);
      H5Sclose(DataspaceId);
      return ERROR_DATASET_INVALID;
    }

    {
      hsize_t   *OldDims = NULL, *OldMaxDims = NULL;
      int        Rank;
      status = xfpGetSimpleDataspaceInfo(DataspaceId, &Rank, &OldDims, &OldMaxDims);
      H5Dclose(DatasetId);
      H5Sclose(DataspaceId);
      if (status < 0) {
        return status;
      }

      if (OldDims == NULL ||  /* old dims must be non-null */
          OldDims[0] < a_NumTimes) /* num times must be within range of OldDims */
      {
        return -1;
      }
    }
  }
  else if (numTimes < a_NumTimes) /* compare NumTimes value */
    return -1;

  status = xfpWriteAttributeInt(a_Id, DATASET_NUMTIMES, 1, &a_NumTimes);

  return status;
} /* xfSetDatasetNumTimes */
/******************************************************************************/
/*   FUNCTION  xfInitializeScalarTimestep */
/*! PURPOSE:   Setup to write a scalar dataset in parts
 *   - NOTES:     A negative error is returned if the timestep cannot be created
 *                timestepId is returned and can be used to write data in the
 *                timestep
 *                minvalue and maxvalue need to be set to the minimum and 
 *                maximum values of the timestep */
/******************************************************************************/
XMDF_API xid xfInitializeScalarTimestep(xid xScalarAId, double dTime,
                    int nValues, float minvalue, float maxvalue, 
                    hsize_t *timestepId)
{
  return xfiInitializeScalarTimestep(
    xScalarAId, 
    dTime, 
    nValues, 
    DATASET_TYPE_FLOAT, 
    minvalue, 
    maxvalue, 
    timestepId);
} /* xfInitializeScalarTimestep */
/******************************************************************************/
/*   FUNCTION  xfWriteScalarTimestepPortion */
/*! PURPOSE:   Write values into a precreated xmdf scalar dataset
 *   - NOTES:     A negative error is returned if the write fails
 *                timestepId is a 1 based index into the dataset and is returned
 *                from the Initialize function */
/******************************************************************************/
XMDF_API xid xfWriteScalarTimestepPortion(xid xScalarId, int timestepId,
                    int numValuesToWrite, int startIndex, const float *a_values)
{
  return xfiWriteScalarTimestepPortion(xScalarId, timestepId, numValuesToWrite,
            startIndex, DATASET_TYPE_FLOAT, a_values);
} /* xfWriteScalarTimestepPortion */
/******************************************************************************/
/*   FUNCTION  xfSetDatasetTimestepMinMax */
/*! PURPOSE:   overwrite min/max values for a specific timestep
 *   - NOTES:     A negative error is returned if the write fails
 *                timestepId is a 1 based index into the dataset and is returned
 *                from the Initialize function */
/******************************************************************************/
XMDF_API xid xfSetDatasetTimestepMinMax(xid xDatasetId, int timestepId,
                                        float minvalue, float maxvalue)
{ 
  return xfiSetDatasetTimestepMinMax(xDatasetId, timestepId, DATASET_TYPE_FLOAT,
                                     minvalue, maxvalue);
} /* xfSetDatasetTimestepMinMax */
/******************************************************************************/
/*   FUNCTION  xfGetDatasetNumVals */
/*! PURPOSE:
 *   - NOTES:     */
/******************************************************************************/
XMDF_API xid xfGetDatasetNumVals (xid a_Id, int *a_Numvals)
{
  int       status = XTRUE;
  xid       ValuesId, SpaceId;
  herr_t    Simple;
  int       Rank;
  hsize_t  *Dims = NULL, *Maxdims = NULL;
  int       i, k;
  xid       ISpaceId, IndicesId, IMemspaceId;
  hsize_t   Count[2];
  hsize_t   Offset[2], Memsize;
  int       *Indices;
  hssize_t   numvalues;
  hssize_t   maxvalues;
  hsize_t    numloops;
  hsize_t    start;

  /* see if it is a scalar dataset */
  if (xfpIsGroupOfType(a_Id, GROUP_TYPE_DATASET_SCALAR)) {
    /* open the values dataset */
    ValuesId = H5Dopen1(a_Id, DATASET_DSET_VALUES);
    if (ValuesId < 0) {
      return ERROR_DATASET_INVALID;
    }
  
    /* Get the dataspace and make sure it is a simple dataset of rank 2 */
    SpaceId = H5Dget_space(ValuesId);
    if (SpaceId < 0) {
      H5Dclose(ValuesId);
      return ERROR_DATASET_INVALID;
    }

    Simple = H5Sis_simple(SpaceId);
    if (!Simple) {
      H5Sclose(SpaceId);
      H5Dclose(ValuesId);
      return ERROR_DATASET_INVALID;
    }

    xfpGetSimpleDataspaceInfo(SpaceId, &Rank, &Dims, &Maxdims);
    if (Rank != 2) {
      H5Sclose(SpaceId);
      H5Dclose(ValuesId);
      return ERROR_DATASET_INVALID;
    }
    /* See if the Dataset was written out in parallel */
    /* Because if it was then it was written differently and to get the data we */
    /* do something different. We know it was written out in parallel if the indices */
    /* array exists. */
    IndicesId = H5Dopen1(a_Id, DATASET_DSET_INDICES);
    if (IndicesId > -1) {
       /* read the indices array to read the values in correctly */
      /* It is defaulted at 2,000,000 values */
      maxvalues = xfpGetParallelNumValuesToRead();
      *a_Numvals=0;
      numloops = (Dims[1] / maxvalues + 1);
      for (i=0; i<numloops; i=i+1) {
        
        /* get how many values to read */
        if ((int)(Dims[1]) > (i+1) * maxvalues) {
          numvalues = maxvalues;
        }
        else {
          numvalues = (Dims[1] - i * maxvalues);
        }

        /* get where to start reading the values from the dataset */
        start = 0 + i * maxvalues;

        /* create the array to hold the data */
        Indices = (int*) malloc((size_t)numvalues*sizeof(int));
        
        /*Read the array*/
        ISpaceId = H5Dget_space(IndicesId);
        if (ISpaceId < 0) {
          H5Dclose(SpaceId);
          H5Dclose(ValuesId);
          return ERROR_DATASET_INVALID;
        }


        /* select the hyperslab */
        Count[0] = 1;
        Count[1] = numvalues;

        Offset[0] = 0;
        Offset[1] = start;

        status = H5Sselect_hyperslab(ISpaceId, H5S_SELECT_SET, Offset, NULL, Count, 
                                    NULL);
        if (status < 0) {
          H5Sclose(ISpaceId);
          H5Sclose(SpaceId);
          H5Dclose(SpaceId);
          H5Dclose(ValuesId);
          return status;
        }

        /* Make sure the portion is valid */
        if (H5Sselect_valid(ISpaceId) <= 0) {
          return ERROR_DATASET_SIZE_INCORRECT;
        }

        /* Create a memory dataspace for the data */
        Memsize = H5Sget_select_npoints(ISpaceId);
        IMemspaceId = H5Screate_simple(1, &Memsize, NULL);
        if (IMemspaceId < 0) {
          H5Sclose(ISpaceId);
          H5Dclose(IndicesId);
          H5Sclose(SpaceId);
          H5Dclose(ValuesId);
          return IMemspaceId;
        }

        /* Read the data */
        status = H5Dread(IndicesId, H5T_NATIVE_INT, IMemspaceId, ISpaceId,
                        H5P_DEFAULT, Indices);

        /*count the values that do not have a -1 as a value */
        for (k=0; k<numvalues; k++) {
          if (Indices[k] != -1) {
            *a_Numvals=*a_Numvals+1;
          }
        }
        H5Sclose(ISpaceId);
        H5Sclose(IMemspaceId);
        if (Indices) {
          free(Indices);
        }
      }
      H5Dclose(IndicesId);
    }
    else {
      *a_Numvals = (int)Dims[1];
    }
    H5Sclose(SpaceId);
    H5Dclose(ValuesId);
  }
  /* vector dataset */
  else if (xfpIsGroupOfType(a_Id, GROUP_TYPE_DATASET_VECTOR)) {
    /* open the values dataset */
    ValuesId = H5Dopen1(a_Id, DATASET_DSET_VALUES);
    if (ValuesId < 0) {
      return ERROR_DATASET_INVALID;
    }
  
    /* Get the dataspace and make sure it is a simple dataset of rank 3 */
    SpaceId = H5Dget_space(ValuesId);
    if (SpaceId < 0) {
      H5Dclose(ValuesId);
      return ERROR_DATASET_INVALID;
    }

    Simple = H5Sis_simple(SpaceId);
    if (!Simple) {
      H5Sclose(SpaceId);
      H5Dclose(ValuesId);
      return ERROR_DATASET_INVALID;
    }

    xfpGetSimpleDataspaceInfo(SpaceId, &Rank, &Dims, &Maxdims);
    if (Rank != 3) {
      H5Sclose(SpaceId);
      H5Dclose(ValuesId);
      return ERROR_DATASET_INVALID;
    }

    *a_Numvals = (int)Dims[1];
    H5Sclose(SpaceId);
    H5Dclose(ValuesId);
  }
  else {
    /* group is not a valid dataset type */
    return ERROR_DATASET_INVALID;
  }

  if (Dims) {
    free(Dims);
    Dims = NULL;
  }
  if (Maxdims) {
    free(Maxdims);
    Maxdims = NULL;
  }

  return status;
} /* xfGetDatasetNumVals */
/******************************************************************************/
/*   FUNCTION  xfGetDatasetNumActive */
/*! PURPOSE:
 *   - NOTES:     */
/******************************************************************************/
XMDF_API xid xfGetDatasetNumActive (xid a_Id, int *a_NumActivevals)
{
  int    status = XTRUE;
  xid    ActiveId, SpaceId;
  hsize_t *Dims = NULL, *Maxdims = NULL;
  herr_t  Simple;
  int     Rank;
  herr_t (*old_func)(void*);
  void  *old_client_data;

  /* make sure group is a dataset group */
  if (!xfpIsGroupOfType(a_Id, GROUP_TYPE_DATASET_SCALAR) &&
      !xfpIsGroupOfType(a_Id, GROUP_TYPE_DATASET_VECTOR)) {
    return ERROR_DATASET_INVALID;
  }

  /* open the activity dataset */

  H5Eget_auto1(&old_func, &old_client_data);

  /* Turn off error handling */
  H5Eset_auto1(NULL, NULL);

  ActiveId = H5Dopen1(a_Id, DATASET_DSET_ACTIVE);

  /* Restore previous error handler */
  H5Eset_auto1(old_func, old_client_data);

  if (ActiveId < 0) {
    /* If there is no activity dataset everything is active, set numvals to 0 */
    *a_NumActivevals = 0;
    return XTRUE;
  }

  /* Get the dataspace and make sure it is a simple dataset of rank 2 */
  SpaceId = H5Dget_space(ActiveId);
  if (SpaceId < 0) {
    H5Dclose(ActiveId);
    return ERROR_DATASET_INVALID;
  }

  Simple = H5Sis_simple(SpaceId);
  if (!Simple) {
    H5Sclose(SpaceId);
    H5Dclose(ActiveId);
    return ERROR_DATASET_INVALID;
  }

  xfpGetSimpleDataspaceInfo(SpaceId, &Rank, &Dims, &Maxdims);
  if (Rank != 2) {
    H5Sclose(SpaceId);
    H5Dclose(ActiveId);
    return ERROR_DATASET_INVALID;
  }

  *a_NumActivevals = (int)Dims[1];
  H5Sclose(SpaceId);
  H5Dclose(ActiveId);

  free(Dims);
  free(Maxdims);

  return status;
} /* xfGetDatasetNumActive */
/******************************************************************************/
/*   FUNCTION  xfGetDatasetVecNumComponents */
/*! PURPOSE:
 *   - NOTES:     */
/******************************************************************************/
XMDF_API xid xfGetDatasetVecNumComponents (xid a_Id, int *a_NumComponents)
{
  int    status = XTRUE;
  xid    ValuesId, SpaceId;
  herr_t Simple;
  int    Rank;
  hsize_t *Dims = NULL, *Maxdims = NULL;

  /* Make sure it is a vector dataset */
  if (!xfpIsGroupOfType(a_Id, GROUP_TYPE_DATASET_VECTOR)) {
    /* group is not a valid dataset type */
    return ERROR_DATASET_INVALID;
  }

  /* open the values dataset */
  ValuesId = H5Dopen1(a_Id, DATASET_DSET_VALUES);
  if (ValuesId < 0) {
    return ERROR_DATASET_INVALID;
  }

  /* Get the dataspace and make sure it is a simple dataset of rank 3 */
  SpaceId = H5Dget_space(ValuesId);
  if (SpaceId < 0) {
    H5Dclose(ValuesId);
    return ERROR_DATASET_INVALID;
  }

  Simple = H5Sis_simple(SpaceId);
  if (!Simple) {
    H5Sclose(SpaceId);
    H5Dclose(ValuesId);
    return ERROR_DATASET_INVALID;
  }

  xfpGetSimpleDataspaceInfo(SpaceId, &Rank, &Dims, &Maxdims);
  if (Rank != 3) {
    H5Sclose(SpaceId);
    H5Dclose(ValuesId);
    return ERROR_DATASET_INVALID;
  }

  *a_NumComponents = (int)Dims[2];
  H5Sclose(SpaceId);
  H5Dclose(ValuesId);

  if (Dims) {
    free(Dims);
    Dims = NULL;
  }
  if (Maxdims) {
    free(Maxdims);
    Maxdims = NULL;
  }

  return status;
} /* xfGetDatasetVecNumComponents */

/******************************************************************************/
/*   FUNCTION  xfGetDatasetTimeUnits */
/*! PURPOSE:
 *   - NOTES:     The units variable should arleady be allocated. 
 *             size >= TIME_UNITS_MAXLENGTH */
/******************************************************************************/
XMDF_API xid xfGetDatasetTimeUnits (xid a_Id, char *Units)
{
  int    status;
  char  *ReadUnits = NULL;

  /* make sure group is a dataset group */
  if (!xfpIsGroupOfType(a_Id, GROUP_TYPE_DATASET_SCALAR) &&
      !xfpIsGroupOfType(a_Id, GROUP_TYPE_DATASET_VECTOR)) {
    return ERROR_DATASET_INVALID;
  }

  /* open the Time attribute */
  status = xfpReadAttributeString(a_Id, TIME_UNITS, 1, &ReadUnits);
  if (status < 0 || strlen(ReadUnits) > TIME_UNITS_MAXLENGTH) {
    if (ReadUnits) {
      free(ReadUnits);
    }
    return ERROR_DATASET_INVALID;
  }

  strcpy(Units, ReadUnits);

  free(ReadUnits);
  return status;
} /* xfGetDatasetTimeUnits */

/******************************************************************************/
/*   FUNCTION  xfGetDatasetUnits */
/*! PURPOSE:
 *   - NOTES:     The units variable should arleady be allocated. 
 *             size >= UNITS_MAXLENGTH */
/******************************************************************************/
XMDF_API xid xfGetDatasetUnits (xid a_Id, char *Units)
{
  int    status;
  char  *ReadUnits = NULL;

  /* make sure group is a dataset group */
  if (!xfpIsGroupOfType(a_Id, GROUP_TYPE_DATASET_SCALAR) &&
      !xfpIsGroupOfType(a_Id, GROUP_TYPE_DATASET_VECTOR)) {
    return ERROR_DATASET_INVALID;
  }

  /* open the Time attribute */
  status = xfpReadAttributeString(a_Id, DATASET_ATT_UNITS, 1, &ReadUnits);
  if (status < 0 || strlen(ReadUnits) > UNITS_MAXLENGTH) {
    if (ReadUnits) {
      free(ReadUnits);
    }
    return ERROR_DATASET_INVALID;
  }

  strcpy(Units, ReadUnits);

  free(ReadUnits);
  return status;
} /* xfGetDatasetUnits */
/******************************************************************************/
/*   FUNCTION  xfGetDatasetTimes */
/*! PURPOSE:
 *   - NOTES:     The times array must already be allocated to the correct number. */
/******************************************************************************/
XMDF_API xid xfGetDatasetTimes (xid a_Id, int a_NumTimes, double *a_Times)
{
  int    status;

  /* make sure group is a dataset group */
  if (!xfpIsGroupOfType(a_Id, GROUP_TYPE_DATASET_SCALAR) &&
      !xfpIsGroupOfType(a_Id, GROUP_TYPE_DATASET_VECTOR)) {
    return ERROR_DATASET_INVALID;
  }

  status = xftReadDset1DDouble(a_Id, DATASET_DSET_TIMES, a_NumTimes, a_Times);

  return status;
} /* xfGetDatasetTimes */
/******************************************************************************/
/*   FUNCTION  xfGetDatasetMins */
/*! PURPOSE:
 *   - NOTES:     */
/******************************************************************************/
XMDF_API xid xfGetDatasetMins(xid a_Id, int a_NumTimes, float *a_Mins)
{
  return xfGetDatasetMinsFloat(a_Id, a_NumTimes, a_Mins);

} /* xfGetDatasetMins */
/******************************************************************************/
/*   FUNCTION  xfGetDatasetMaxs */
/*! PURPOSE:
 *   - NOTES:     */
/******************************************************************************/
XMDF_API xid xfGetDatasetMaxs (xid a_Id, int a_NumTimes, float *a_Maxs)
{
  return xfGetDatasetMaxsFloat(a_Id, a_NumTimes, a_Maxs);

} /* xfGetDatasetMaxs */
/******************************************************************************/
/*   FUNCTION  xfGetDatasetMinsFloat */
/*! PURPOSE:
 *   - NOTES:     */
/******************************************************************************/
XMDF_API xid xfGetDatasetMinsFloat (xid a_Id, int a_NumTimes, float *a_Mins)
{
  int    status;

  /* make sure group is a dataset group */
  if (!xfpIsGroupOfType(a_Id, GROUP_TYPE_DATASET_SCALAR) &&
      !xfpIsGroupOfType(a_Id, GROUP_TYPE_DATASET_VECTOR)) {
    return ERROR_DATASET_INVALID;
  }

  /* read the data */
  status = xftReadDset1DFloat(a_Id, DATASET_DSET_MINS, a_NumTimes, a_Mins);

  return status;
} /* xfGetDatasetMinsFloat */
/******************************************************************************/
/*   FUNCTION  xfGetDatasetMaxsFloat */
/*! PURPOSE:
 *   - NOTES:     */
/******************************************************************************/
XMDF_API xid xfGetDatasetMaxsFloat (xid a_Id, int a_NumTimes, float *a_Maxs)
{
  int    status;

  /* make sure group is a dataset group */
  if (!xfpIsGroupOfType(a_Id, GROUP_TYPE_DATASET_SCALAR) &&
      !xfpIsGroupOfType(a_Id, GROUP_TYPE_DATASET_VECTOR)) {
    return ERROR_DATASET_INVALID;
  }

  /* read the data */
  status = xftReadDset1DFloat(a_Id, DATASET_DSET_MAXS, a_NumTimes, a_Maxs);

  return status;
} /* xfGetDatasetMaxsFloat */
/******************************************************************************/
/*   FUNCTION  xfGetDatasetMinsDouble */
/*! PURPOSE:
 *   - NOTES:     */
/******************************************************************************/
XMDF_API xid xfGetDatasetMinsDouble (xid a_Id, int a_NumTimes, double *a_Mins)
{
  int    status;

  /* make sure group is a dataset group */
  if (!xfpIsGroupOfType(a_Id, GROUP_TYPE_DATASET_SCALAR) &&
      !xfpIsGroupOfType(a_Id, GROUP_TYPE_DATASET_VECTOR)) {
    return ERROR_DATASET_INVALID;
  }

  /* read the data */
  status = xftReadDset1DDouble(a_Id, DATASET_DSET_MINS, a_NumTimes, a_Mins);

  return status;
} /* xfGetDatasetMinsDouble */
/******************************************************************************/
/*   FUNCTION  xfGetDatasetMaxsDouble */
/*! PURPOSE:
 *   - NOTES:     */
/******************************************************************************/
XMDF_API xid xfGetDatasetMaxsDouble (xid a_Id, int a_NumTimes, double *a_Maxs)
{
  int    status;

  /* make sure group is a dataset group */
  if (!xfpIsGroupOfType(a_Id, GROUP_TYPE_DATASET_SCALAR) &&
      !xfpIsGroupOfType(a_Id, GROUP_TYPE_DATASET_VECTOR)) {
    return ERROR_DATASET_INVALID;
  }

  /* read the data */
  status = xftReadDset1DDouble(a_Id, DATASET_DSET_MAXS, a_NumTimes, a_Maxs);

  return status;
} /* xfGetDatasetMaxsDouble */
/******************************************************************************/
/*   FUNCTION  xfGetDatasetActivityTimestep */
/*! PURPOSE:   Read the activity values for a specific timestep
 *   - NOTES:     a_Values must already be allocated of size a_NumVals.      */
/******************************************************************************/
XMDF_API xid xfReadActivityTimestep (xid a_Id, 
                  int a_TimestepIndex, int a_NumVals,
                  xmbool *a_Values)
{
  int    status;

  /* make sure group is a dataset group */
  if (!xfpIsGroupOfType(a_Id, GROUP_TYPE_DATASET_SCALAR) &&
      !xfpIsGroupOfType(a_Id, GROUP_TYPE_DATASET_VECTOR)) {
    return ERROR_DATASET_INVALID;
  }

  status = xftReadDset2DUCharPortion(a_Id, DATASET_DSET_ACTIVE, 
                                      a_TimestepIndex-1, 1, 0,
                                      a_NumVals, a_Values);

  return status;
} /* xfReadActivityTimestep */
/******************************************************************************/
/*   FUNCTION  xfReadActivityValuesAtIndex */
/*! PURPOSE:   Read the activity values for specific indices
 *   - NOTES:     a_Values must already be allocated of size a_NumTimes.      */
/******************************************************************************/
XMDF_API xid xfReadActivityValuesAtIndex (xid a_Id, 
                  int a_Index, int a_FirstTime, int a_NumTimes,
                  xmbool *a_Values)
{
  int    status = 1;

  /* make sure group is a dataset group */
  if (!xfpIsGroupOfType(a_Id, GROUP_TYPE_DATASET_SCALAR) &&
      !xfpIsGroupOfType(a_Id, GROUP_TYPE_DATASET_VECTOR)) {
    return ERROR_DATASET_INVALID;
  }

  status = xftReadDset2DUCharPortion(a_Id, DATASET_DSET_ACTIVE, 
              a_FirstTime-1, a_NumTimes, a_Index-1, 1, a_Values);

  return status;
} /* xfReadActivityValuesAtIndex */
/******************************************************************************/
/*   FUNCTION  xfReadScalarValuesTimestep */
/*! PURPOSE:   Read the scalar values for a specific timestep
 *   - NOTES:     a_Values must already be allocated of size a_NumVals.      */
/******************************************************************************/
XMDF_API xid xfReadScalarValuesTimestep (xid a_Id, 
                  int a_TimestepIndex, int a_NumVals, float *a_Values)
{
  return xfReadScalarValuesTimestepFloat(a_Id, a_TimestepIndex, a_NumVals, a_Values);

} /* xfReadScalarValuesTimestep */
/*!***** (OVERLOAD) ************************************************************/
XMDF_API xid xfReadScalarValuesTimestepFloat (xid a_Id, 
                  int a_TimestepIndex, int a_NumVals, float *a_Values)
{
  int    status;

  /* make sure group is a dataset group */
  if (!xfpIsGroupOfType(a_Id, GROUP_TYPE_DATASET_SCALAR)) {
    return ERROR_DATASET_INVALID;
  }

  status = xfReadScalarValuesTimestepFloatPortion(a_Id, a_TimestepIndex, 1,
                                                  a_NumVals, a_Values);

  return status;
} /* xfReadScalarValuesTimestepFloat */
/******************************************************************************/
/*   FUNCTION  xfReadScalarValuesTimestepFloatPortion */
/*! PURPOSE:   Read the scalar values for a specific timestep
 *   - NOTES:    */
/******************************************************************************/
XMDF_API xid xfReadScalarValuesTimestepFloatPortion (xid a_Id, 
                  int a_TimestepIndex, int a_Start, int a_NumVals,
                  float *a_Values)
{
  return xfiReadScalarValuesTimestepPortion(a_Id, a_TimestepIndex, a_Start, 
                                             a_NumVals, NULL, a_Values);

}
/*!***** (OVERLOAD) ************************************************************/
XMDF_API xid xfReadScalarValuesTimestepDouble (xid a_Id, 
                  int a_TimestepIndex, int a_NumVals, double *a_Values)
{
  return xfReadScalarValuesTimestepDoublePortion(a_Id, a_TimestepIndex, 1, 
                                                 a_NumVals, a_Values);
} /* xfReadScalarValuesTimestepDouble */
/******************************************************************************/
/*   FUNCTION  xfReadScalarValuesTimestepFloatPortion */
/*! PURPOSE:   Read the scalar values for a specific timestep
 *   - NOTES:    */
/******************************************************************************/
XMDF_API xid xfReadScalarValuesTimestepDoublePortion (xid a_Id, 
                  int a_TimestepIndex, int a_Start, int a_NumVals,
                  double *a_Values)
{
  return xfiReadScalarValuesTimestepPortion(a_Id, a_TimestepIndex, 1, 
                                             a_NumVals, a_Values, NULL);

}
/******************************************************************************/
/*   FUNCTION  xfiReadScalarValuesTimestepPortion */
/*! PURPOSE:   Read the scalar values for a specific timestep
 *   - NOTES:    */
/******************************************************************************/
static xid xfiReadScalarValuesTimestepPortion (xid a_Id, 
                  int a_TimestepIndex, int a_Start, int a_NumVals,
                  double *a_dValues, float *a_fValues)
{
  int        status;
  xid        IndicesId;
  double    *dArray;
  float     *fArray;
  int       *Indices;
  hssize_t   numvalues;
  hssize_t   maxvalues;
  hsize_t    numloops;
  hsize_t    start;
  int        Rank; 
  hsize_t   *Dims, *MaxDims;
  int        i, j;

  status = -1;

  /* make sure group is a dataset group */
  if (!xfpIsGroupOfType(a_Id, GROUP_TYPE_DATASET_SCALAR)) {
    return ERROR_DATASET_INVALID;
  }

  /* See if the Dataset was written out in parallel */
  /* Because if it was then it was written differently and to get the data we */
  /* do something different. We know it was written out in parallel if the indices */
  /* array exists. */
  IndicesId = H5Dopen1(a_Id, DATASET_DSET_INDICES);
  if (IndicesId > -1) {
    H5Dclose(IndicesId);
    fArray = NULL;
    dArray = NULL;
    /* For simpilcity and to not overshoot the ram */
    /* I will read only so many values of the array */
    /* and then loop throught them the values */
    /* It is defaulted at 2,000,000 values */
    maxvalues = xfpGetParallelNumValuesToRead();

    xftGetSimpleDataspaceInfoFromName(a_Id, DATASET_DSET_VALUES, 0,
          &Rank, &Dims, &MaxDims);

    numloops = (Dims[1] / maxvalues + 1);
    for (i=0; i<numloops; i=i+1) {
      
      /* get how many values to read */
      if (Dims[1] > (int)((i+1) * maxvalues)) {
        numvalues = maxvalues;
      }
      else {
        numvalues = (Dims[1] - i * maxvalues);
        if (!numvalues) break;
      }
      
      /* get where to start reading the values from the dataset */
      start = 0 + i * maxvalues;

      /* read the indices array to read the values in correctly */
      /* create the array to hold the data */
      Indices = (int*) malloc((size_t)numvalues*sizeof(int));/* malloc outside loop !!! */

      xftReadDset2DIntPortion(a_Id, DATASET_DSET_INDICES, a_TimestepIndex-1, 1,
                              start, numvalues, Indices);    

      if (a_dValues) {
        /* create the array to hold the data */
        dArray = (double*) malloc((size_t)numvalues*sizeof(double));
        status = xftReadDset2DDoublePortion(a_Id, DATASET_DSET_VALUES,
                                            a_TimestepIndex-1, 1, start,
                                            numvalues, dArray);
      }
      else {
        /* create the array to hold the data */
        fArray = (float*) malloc((size_t)numvalues*sizeof(float));
        status = xftReadDset2DFloatPortion(a_Id, DATASET_DSET_VALUES,
                                           a_TimestepIndex-1, 1, start,
                                           numvalues, fArray);
      }
      if (status < 0)
        {
          if (Indices) free(Indices);
          if (dArray) free(dArray);
          if (fArray) free(fArray);
          return status;
        }

      for (j = 0; j < numvalues; j = j+ 1) {
        if (Indices[j] != -1 && Indices[j] > a_Start-1 &&
            Indices[j]-1 < a_Start-1 + a_NumVals) {
          if (a_dValues) {
            a_dValues[Indices[j]-1] = dArray[j];
          }
          else {
            a_fValues[Indices[j]-1] = fArray[j];
          }
        }
      }

      /* free the Arrays */
      free(Indices);
      if (dArray) {
        free(dArray);
      }
      if (fArray) {
        free(fArray);
      }
    }
  }
  else {
    if (a_dValues) {
      status = xftReadDset2DDoublePortion(a_Id, DATASET_DSET_VALUES, 
                                          a_TimestepIndex-1, 1, a_Start-1,
                                          a_NumVals, a_dValues);
    }
    else {
      status = xftReadDset2DFloatPortion(a_Id, DATASET_DSET_VALUES, 
                                          a_TimestepIndex-1, 1, a_Start-1,
                                          a_NumVals, a_fValues);
    }
  }

  return status;
} /* xfReadScalarValuesTimestepDouble */

/*!***** (OVERLOAD) ************************************************************/
XMDF_API xid xfReadScalarValuesTimestepInt (xid a_Id, 
                  int a_TimestepIndex, int a_NumVals, int *a_Values)
{
  int    status;

  /* make sure group is a dataset group */
  if (!xfpIsGroupOfType(a_Id, GROUP_TYPE_DATASET_SCALAR)) {
    return ERROR_DATASET_INVALID;
  }

  status = xftReadDset2DIntPortion(a_Id, DATASET_DSET_VALUES, 
                                      a_TimestepIndex-1, 1, 0,
                                      a_NumVals, a_Values);

  return status;
} /* xfReadScalarValuesTimestepDouble */
/******************************************************************************/
/*   FUNCTION  xfReadScalarValuesAtIndex */
/*! PURPOSE:   Read the scalar values for a specific index
 *   - NOTES:     a_Values must already be allocated of size a_NumTimes.      */
/******************************************************************************/
XMDF_API xid xfReadScalarValuesAtIndex (xid a_Id, 
                    int a_Index, int a_FirstTime, int a_NumTimes,
                    float *a_Values)
{
  int    status = 1;

  status = xfReadScalarValuesAtIndexFloat(a_Id, a_Index, a_FirstTime,
                                          a_NumTimes, a_Values);

  return status;
} /* xfReadScalarValuesAtIndex */
/*!***** (OVERLOAD) ************************************************************/
XMDF_API xid xfReadScalarValuesAtIndexFloat (xid a_Id, 
                    int a_Index, int a_FirstTime, int a_NumTimes,
                    float *a_Values)
{
  int        status;
  int        newIndex;
  xid        IndicesId;

  status = 1;

  /* make sure group is a dataset group */
  if (!xfpIsGroupOfType(a_Id, GROUP_TYPE_DATASET_SCALAR)) {
    return ERROR_DATASET_INVALID;
  }

  /* See if the Dataset was written out in parallel */
  /* Because if it was then it was written differently and to get the data we */
  /* do something different. We know it was written out in parallel if the indices */
  /* array exists. */
  IndicesId = H5Dopen1(a_Id, DATASET_DSET_INDICES);
  if (IndicesId > -1) {
    H5Dclose(IndicesId);
    
    xfiGetParIndexFromIndex(a_Id, a_FirstTime-1, a_Index-1, &newIndex);

    status = xftReadDset2DFloatPortion(a_Id, DATASET_DSET_VALUES, a_FirstTime-1,
                                      a_NumTimes, newIndex, 1, a_Values);  
  }
  else {
    status = xftReadDset2DFloatPortion(a_Id, DATASET_DSET_VALUES, a_FirstTime-1,
                                      a_NumTimes, a_Index - 1, 1, a_Values);
  }

  return status;
} /* xfReadScalarValuesAtIndexFloat */
/*!***** (OVERLOAD) ************************************************************/
XMDF_API xid xfReadScalarValuesAtIndexDouble (xid a_Id, 
                    int a_Index, int a_FirstTime, int a_NumTimes,
                    double *a_Values)
{
  int        status;
  int        newIndex;
  xid        IndicesId;

  /* make sure group is a dataset group */
  if (!xfpIsGroupOfType(a_Id, GROUP_TYPE_DATASET_SCALAR)) {
    return ERROR_DATASET_INVALID;
  }

  /* See if the Dataset was written out in parallel */
  /* Because if it was then it was written differently and to get the data we */
  /* do something different. We know it was written out in parallel if the indices */
  /* array exists. */
  IndicesId = H5Dopen1(a_Id, DATASET_DSET_INDICES);
  if (IndicesId > -1) {
    H5Dclose(IndicesId);
    
    xfiGetParIndexFromIndex(a_Id, a_FirstTime-1, a_Index-1, &newIndex);

    status = xftReadDset2DDoublePortion(a_Id, DATASET_DSET_VALUES, a_FirstTime-1,
                                      a_NumTimes, newIndex, 1, a_Values);     
  }
  else {
    status = xftReadDset2DDoublePortion(a_Id, DATASET_DSET_VALUES, a_FirstTime-1,
                                      a_NumTimes, a_Index - 1, 1, a_Values);
  }

  return status;
} /* xfReadScalarValuesAtIndexDouble */
/*!***** (OVERLOAD) ************************************************************/
XMDF_API xid xfReadScalarValuesAtIndexInt (xid a_Id,
                    int a_Index, int a_FirstTime, int a_NumTimes,
                    int *a_Values)
{
  int        status;
  int        newIndex;
  xid        IndicesId;

  /* make sure group is a dataset group */
  if (!xfpIsGroupOfType(a_Id, GROUP_TYPE_DATASET_SCALAR)) {
    return ERROR_DATASET_INVALID;
  }

  /* See if the Dataset was written out in parallel */
  /* Because if it was then it was written differently and to get the data we */
  /* do something different. We know it was written out in parallel if the indices */
  /* array exists. */
  IndicesId = H5Dopen1(a_Id, DATASET_DSET_INDICES);
  if (IndicesId > -1) {
    H5Dclose(IndicesId);
    
    xfiGetParIndexFromIndex(a_Id, a_FirstTime-1, a_Index-1, &newIndex);

    status = xftReadDset2DIntPortion(a_Id, DATASET_DSET_VALUES, a_FirstTime-1,
                                      a_NumTimes, newIndex, 1, a_Values);     
  }
  else {
    status = xftReadDset2DIntPortion(a_Id, DATASET_DSET_VALUES, a_FirstTime-1,
                                      a_NumTimes, a_Index - 1, 1, a_Values);
  }

  return status;
} /* xfReadScalarValuesAtIndexInt */
/******************************************************************************/
/*   FUNCTION  xfiReadWriteScalarValuesFloatAtIndices */
/*! PURPOSE:   Read the scalar values for a specific index
 *   - NOTES:     a_Values must already be allocated of size a_NumTimes.      */
/******************************************************************************/
static xid xfiReadWriteScalarValuesFloatAtIndices (xid a_Id, 
                    ReadWrite_enum a_readWrite,
                    int a_nIndices, const int *a_Indices, int a_FirstTime,
                    int a_NumTimes, float *a_Values)
{
  xid status = -1;
  size_t numIndices2D;
  hsize_t *indices2D=NULL;
  size_t pos2D;
  int iTime;
  int index;
  int timeIndex;

  /* make sure group is a dataset group */
  if (!xfpIsGroupOfType(a_Id, GROUP_TYPE_DATASET_SCALAR)) {
    return ERROR_DATASET_INVALID;
  }

  numIndices2D = a_nIndices*a_NumTimes; 
    // 2 for the # of indices for dataset
  indices2D = (hsize_t*)malloc(numIndices2D*2*sizeof(hsize_t));
  pos2D = 0;
  for (iTime = 0; iTime < a_NumTimes; ++iTime) {
    timeIndex = iTime + a_FirstTime - 1;
    for (index = 0; index < a_nIndices; ++index) {
      indices2D[pos2D] = timeIndex;
      ++pos2D;
      indices2D[pos2D] = a_Indices[index] - 1;
      ++pos2D;
    }
  }
  status = xftReadWriteDsetFloatIndices (a_Id, a_readWrite,
              DATASET_DSET_VALUES, numIndices2D, indices2D, 2, a_Values);

  free(indices2D);

  return status;
} /* xfiReadWriteScalarValuesFloatAtIndices */
/******************************************************************************/
/*   FUNCTION  xfiGetParIndexFromIndex */
/*! PURPOSE:   Gets the index for a value from the index array that was written
 *             out if it came from the parallel world

 *   - NOTES:         */
/******************************************************************************/
static void xfiGetParIndexFromIndex (xid a_Id, int a_TimestepIndex,
                             int a_index, int *a_parindex)
{
  xmbool     indexfound;
  int       *Indices;
  hssize_t   numvalues;
  hssize_t   maxvalues;
  hsize_t    numloops;
  hsize_t    start;
  int        Rank; 
  hsize_t   *Dims, *MaxDims;
  int        i, j;
    

  indexfound = XFALSE;

  /* For simpilcity and to not overshoot the ram */
  /* I will read only so many values of the array */
  /* and then loop throught them the values */
  /* It is defaulted at 2,000,000 values */
  maxvalues = xfpGetParallelNumValuesToRead();

  xfpGetSimpleDataspaceInfoFromName(a_Id, DATASET_DSET_VALUES,
                                     &Rank, &Dims, &MaxDims);
  numloops = (Dims[1] / maxvalues + 1);
  for (i=0; i<numloops && !indexfound; i=i+1) {
  
    /* get how many values to read */
    if ((int)(Dims[1]) > (int)((i+1) * maxvalues)) {
      numvalues = maxvalues;
    }
    else {
      numvalues = (Dims[1] - i * maxvalues);
    }
  
    /* get where to start reading the values from the dataset */
    start = 0 + i * maxvalues;

    /* read the indices array to read the values in correctly */
    /* create the array to hold the data */
    Indices = (int*) malloc((size_t)numvalues*sizeof(int));

    xfpReadDset2DIntPortion(a_Id, DATASET_DSET_INDICES, a_TimestepIndex, 1,
                            start, numvalues, Indices);
  
    /* get the right index from the array */

    for (j = 0; j < numvalues && !indexfound; j = j+ 1) {
      if (Indices[j] != -1 && Indices[j]-1 == a_index) {
        *a_parindex = (int)(start + j);
        indexfound = XTRUE;
      }
    }

    /* free the Arrays */
    free(Indices);
  }
}
/******************************************************************************/
/*   FUNCTION  xfReadScalarValuesAtIndices */
/*! PURPOSE:   Read the scalar values for a specific index
 *   - NOTES:     a_Values must already be allocated of size a_NumTimes.      */
/******************************************************************************/
XMDF_API xid xfReadScalarValuesAtIndices (xid a_Id, 
                    int a_nIndices, const int *a_Indices, int a_FirstTime,
                    int a_NumTimes, float *a_Values)
{
  return xfiReadWriteScalarValuesFloatAtIndices (a_Id, RW_READ,
                    a_nIndices, a_Indices, a_FirstTime, a_NumTimes, a_Values);
} /* xfReadScalarValuesAtIndices */
/*!***** (OVERLOAD) ************************************************************/
XMDF_API xid xfReadScalarValuesAtIndicesFloat (xid a_Id, 
                    int a_nIndices, const int *a_Indices, int a_FirstTime,
                    int a_NumTimes, float *a_Values)
{
  return xfReadScalarValuesAtIndices(a_Id, a_nIndices, a_Indices, a_FirstTime,
                                     a_NumTimes, a_Values);
}
/******************************************************************************/
/*   FUNCTION  xfReadVectorValuesTimestep */
/*! PURPOSE:   Read the vector values for a specific timestep
 *   - NOTES:     a_Values must already be allocated to size
 *             a_NumVals * a_NumComponents */
/******************************************************************************/
XMDF_API xid xfReadVectorValuesTimestep (xid a_Id, 
                  int a_TimestepIndex, int a_NumVals, int a_NumComponents, 
                  float *a_Values)
{
  return xfReadVectorValuesTimestepFloat(a_Id, a_TimestepIndex, a_NumVals,
                                         a_NumComponents, a_Values);
} /* xfGetDatasetVectorValuesFloatTimestep */
/*!***** (OVERLOAD) ************************************************************/
XMDF_API xid xfReadVectorValuesTimestepFloat (xid a_Id, 
                  int a_TimestepIndex, int a_NumVals, int a_NumComponents, 
                  float *a_Values)
{
  int    status = 1;

  /* make sure group is a vector dataset group */
  if (!xfpIsGroupOfType(a_Id, GROUP_TYPE_DATASET_VECTOR)) {
    return ERROR_DATASET_INVALID;
  }

  status = xfReadVectorValuesTimestepFloatPortion(a_Id, a_TimestepIndex, 1, 
                                                  a_NumVals, a_NumComponents,
                                                  a_Values);

  return status;
} /* xfGetDatasetVectorValuesFloatTimestepFloat */
/******************************************************************************/
/*   FUNCTION  xfReadVectorValuesTimestepFloatPortion */
/*! PURPOSE:   Read the vector values for a specific timestep
 *   - NOTES:    */
/******************************************************************************/
XMDF_API xid xfReadVectorValuesTimestepFloatPortion (xid a_Id, 
                  int a_TimestepIndex, int a_Start, int a_NumVals,
                  int a_NumComponents, float *a_Values)
{
  return xfiReadVectorValuesTimestepPortion(a_Id, a_TimestepIndex, a_Start, 
                                            a_NumVals, a_NumComponents,
                                            NULL, a_Values);
}

/*!***** (OVERLOAD) ************************************************************/
XMDF_API xid xfReadVectorValuesTimestepDouble (xid a_Id, 
                  int a_TimestepIndex, int a_NumVals, int a_NumComponents, 
                  double *a_Values)
{

  return xfReadVectorValuesTimestepDoublePortion(a_Id, a_TimestepIndex, 1, 
                                                 a_NumVals, a_NumComponents,
                                                 a_Values);

} /* xfReadScalarValuesTimestepDouble */
/******************************************************************************/
/*   FUNCTION  xfReadVectorValuesTimestepDoublePortion */
/*! PURPOSE:   Read the vector values for a specific timestep
 *   - NOTES:    */
/******************************************************************************/
XMDF_API xid xfReadVectorValuesTimestepDoublePortion (xid a_Id, 
                  int a_TimestepIndex, int a_Start, int a_NumVals,
                  int a_NumComponents, double *a_Values)
{
  return xfiReadVectorValuesTimestepPortion(a_Id, a_TimestepIndex, a_Start, 
                                            a_NumVals, a_NumComponents,
                                            a_Values, NULL);
}
/******************************************************************************/
/*   FUNCTION  xfReadVectorValuesTimestepDoublePortion */
/*! PURPOSE:   Read the vector values for a specific timestep
 *   - NOTES:    */
/******************************************************************************/
static xid xfiReadVectorValuesTimestepPortion (xid a_Id, 
                  int a_TimestepIndex, int a_Start, int a_NumVals,
                  int a_NumComponents, double *a_dValues, float *a_fValues)
{
  int        status = 1;
  xid        IndicesId;
  double    *dArray;
  float     *fArray;
  int       *Indices;
  hssize_t   numvalues;
  hssize_t   maxvalues;
  hsize_t    numloops;
  hsize_t    start;
  int        Rank; 
  hsize_t   *Dims, *MaxDims;
  int        i, k, j;

  /* make sure group is a vector dataset group */
  if (!xfpIsGroupOfType(a_Id, GROUP_TYPE_DATASET_VECTOR)) {
    return ERROR_DATASET_INVALID;
  }

  /* See if the Dataset was written out in parallel */
  /* Because if it was then it was written differently and to get the data we */
  /* do something different. We know it was written out in parallel if the indices */
  /* array exists. */
  IndicesId = H5Dopen1(a_Id, DATASET_DSET_INDICES);
  if (IndicesId > -1) {
    H5Dclose(IndicesId);
    dArray = NULL;
    fArray = NULL;
    /* For simpilcity and to not overshoot the ram */
    /* I will read only so many values of the array */
    /* and then loop throught them the values */
    /* It is defaulted at 2,000,000 values */
    maxvalues = xfpGetParallelNumValuesToRead();

    xftGetSimpleDataspaceInfoFromName(a_Id, DATASET_DSET_VALUES, 0,
                                       &Rank, &Dims, &MaxDims);
    numloops = (Dims[1] / maxvalues + 1);
    for (i=0; i<numloops; i=i+1) {
      
      /* get how many values to read */
      if ((int)Dims[1] > (int)((i+1) * maxvalues)) {
        numvalues = maxvalues;
      }
      else {
        numvalues = (Dims[1] - i * maxvalues);
      }
      
      /* get where to start reading the values from the dataset */
      start = 0 + i * maxvalues;

      /* read the indices array to read the values in correctly */
      /* create the array to hold the data */
      Indices = (int*) malloc((size_t)numvalues*sizeof(int));

      xftReadDset2DIntPortion(a_Id, DATASET_DSET_INDICES, a_TimestepIndex-1, 1,
                              start, numvalues, Indices);    

      /* create the array to hold the data */
      if (a_dValues) {
        dArray = (double*) malloc((size_t)(numvalues*a_NumComponents)*sizeof(double));
        xftReadDset3DDoublePortion(a_Id, DATASET_DSET_VALUES, 
                                   a_TimestepIndex-1, 1, start, numvalues,
                                   0, a_NumComponents, dArray);
      }
      else {
        fArray = (float*) malloc((size_t)(numvalues*a_NumComponents)*sizeof(float));
        xftReadDset3DFloatPortion(a_Id, DATASET_DSET_VALUES, 
                                   a_TimestepIndex-1, 1, start, numvalues,
                                   0, a_NumComponents, fArray);
      }

      /* reset the data */
      for (j=0; j < numvalues; j=j+1) {
        if (Indices[j] != -1 && Indices[j] > a_Start-1 &&
            Indices[j]-1 < a_Start-1 + a_NumVals) {
          /* do each component */
          for (k=0; k<a_NumComponents; k=k+1) {
            if (a_dValues) {
              a_dValues[(Indices[j]-1)*a_NumComponents + k] = dArray[(j*a_NumComponents)+k];
            }
            else {
              a_fValues[(Indices[j]-1)*a_NumComponents + k] = fArray[(j*a_NumComponents)+k];
            }
          }
        }
      }

      /* free the Arrays */
      free(Indices);
      if (dArray) {
        free(dArray);
      }
      if (fArray) {
        free(fArray);
      }
    }
  }
  else {
    if (a_dValues) {
      status = xftReadDset3DDoublePortion(a_Id, DATASET_DSET_VALUES, 
                                          a_TimestepIndex-1, 1, a_Start-1, a_NumVals,
                                          0, a_NumComponents, a_dValues);
    }
    else {
      status = xftReadDset3DFloatPortion(a_Id, DATASET_DSET_VALUES, 
                                         a_TimestepIndex-1, 1, a_Start-1, a_NumVals,
                                         0, a_NumComponents, a_fValues);
    }
  }

  return status;
} /* xfGetDatasetVectorValuesFloatTimestepDouble */

/******************************************************************************/
/*   FUNCTION  xfReadVectorValuesAtIndex */
/*! PURPOSE:   Read the vector values for specific indices
 *   - NOTES:     a_Values must already be allocated of size
 *             a_NumTimes * a_NumComponents      */
/******************************************************************************/
XMDF_API xid xfReadVectorValuesAtIndex (xid a_Id, int a_Index,
                    int a_FirstTime, int a_NumTimes, int a_NumComponents,
                    float *a_Values)
{

  return xfReadVectorValuesAtIndexFloat(a_Id, a_Index, a_FirstTime, a_NumTimes,
                                        a_NumComponents, a_Values);
} /* xfReadVectorValuesAtIndex */
/*!***** (OVERLOAD) ************************************************************/
XMDF_API xid xfReadVectorValuesAtIndexFloat (xid a_Id, int a_Index,
                    int a_FirstTime, int a_NumTimes, int a_NumComponents,
                    float *a_Values)
{
  int    status = 1;
  int    newIndex;
  xid    IndicesId;

  /* make sure group is a dataset group */
  if (!xfpIsGroupOfType(a_Id, GROUP_TYPE_DATASET_VECTOR)) {
    return ERROR_DATASET_INVALID;
  }

  /* See if the Dataset was written out in parallel */
  /* Because if it was then it was written differently and to get the data we */
  /* do something different. We know it was written out in parallel if the indices */
  /* array exists. */
  IndicesId = H5Dopen1(a_Id, DATASET_DSET_INDICES);
  if (IndicesId > -1) {
    H5Dclose(IndicesId);
    
    xfiGetParIndexFromIndex(a_Id, a_FirstTime-1, a_Index-1, &newIndex);

    status = xftReadDset3DFloatPortion(a_Id, DATASET_DSET_VALUES, a_FirstTime-1,
               a_NumTimes, newIndex, 1, 0, a_NumComponents, a_Values); 
  }
  else {

    status = xftReadDset3DFloatPortion(a_Id, DATASET_DSET_VALUES, a_FirstTime-1,
               a_NumTimes, a_Index-1, 1, 0, a_NumComponents, a_Values);
  }


  return status;
} /* xfReadVectorValuesAtIndexFloat */
/*!***** (OVERLOAD) ************************************************************/
XMDF_API xid xfReadVectorValuesAtIndexDouble (xid a_Id, int a_Index,
                    int a_FirstTime, int a_NumTimes, int a_NumComponents,
                    double *a_Values)
{
  int    status = 1;
  int    newIndex;
  xid    IndicesId;

  /* make sure group is a dataset group */
  if (!xfpIsGroupOfType(a_Id, GROUP_TYPE_DATASET_VECTOR)) {
    return ERROR_DATASET_INVALID;
  }

  /* See if the Dataset was written out in parallel */
  /* Because if it was then it was written differently and to get the data we */
  /* do something different. We know it was written out in parallel if the indices */
  /* array exists. */
  IndicesId = H5Dopen1(a_Id, DATASET_DSET_INDICES);
  if (IndicesId > -1) {
    H5Dclose(IndicesId);
    
    xfiGetParIndexFromIndex(a_Id, a_FirstTime-1, a_Index-1, &newIndex);

    status = xftReadDset3DDoublePortion(a_Id, DATASET_DSET_VALUES, a_FirstTime-1,
               a_NumTimes, newIndex, 1, 0, a_NumComponents, a_Values);   
  }
  else {

    status = xftReadDset3DDoublePortion(a_Id, DATASET_DSET_VALUES, a_FirstTime-1,
                 a_NumTimes, a_Index-1, 1, 0, a_NumComponents, a_Values);
  }

  return status;
} /* xfReadVectorValuesAtIndexDouble */
/******************************************************************************/
/*   FUNCTION  xfiReadWriteVectorValuesFloatAtIndices */
/*! PURPOSE:   Read the scalar values for a specific index
 *   - NOTES:     a_Values must already be allocated of size a_NumTimes.      */
/******************************************************************************/
static xid xfiReadWriteVectorValuesFloatAtIndices (xid a_Id, 
                         ReadWrite_enum a_readWrite,
                         int a_nIndices, const int *a_Indices, int a_FirstTime,
                         int a_NumTimes, int a_nComponents, float *a_Values)
{
  xid status = -1;
  size_t numIndices3D;
  hsize_t *indices3D=NULL;
  size_t pos3D;
  int iTime;
  int index;
  int timeIndex;
  int component;

  /* make sure group is a dataset group */
  if (!xfpIsGroupOfType(a_Id, GROUP_TYPE_DATASET_VECTOR)) {
    return ERROR_DATASET_INVALID;
  }

  numIndices3D = a_nIndices*a_NumTimes*a_nComponents; 
    // 2 for the # of indices for dataset
  indices3D = (hsize_t*)malloc(numIndices3D*a_nComponents*2*sizeof(hsize_t));
  pos3D = 0;
  for (iTime = 0; iTime < a_NumTimes; ++iTime) {
    timeIndex = iTime + a_FirstTime - 1;
    for (index = 0; index < a_nIndices; ++index) {
      for (component = 0; component < a_nComponents; ++component) {
        indices3D[pos3D] = timeIndex;
        ++pos3D;
        indices3D[pos3D] = a_Indices[index] - 1;
        ++pos3D;
        indices3D[pos3D] = component;
        ++pos3D;
      }
    }
  }
  status = xftReadWriteDsetFloatIndices (a_Id, a_readWrite,
        DATASET_DSET_VALUES, numIndices3D, indices3D, 3, a_Values);

  free(indices3D);

  return status;
} /* xfiReadWriteVectorValuesFloatAtIndices */
/******************************************************************************/
/*   FUNCTION  xfReadVectorValuesAtIndices */
/*! PURPOSE:   Read the vector values for a specific index
 *   - NOTES:     a_Values must already be allocated of size a_NumTimes.      */
/******************************************************************************/
XMDF_API xid xfReadVectorValuesAtIndices (xid a_Id, 
                    int a_nIndices, const int *a_Indices, int a_FirstTime,
                    int a_NumTimes, int a_NumComponents, float *a_Values)
{
  return xfiReadWriteVectorValuesFloatAtIndices (a_Id, RW_READ,
    a_nIndices, a_Indices, a_FirstTime, a_NumTimes, a_NumComponents, a_Values);
} /* xfReadVectorValuesAtIndices */
/*!***** (OVERLOAD) ************************************************************/
XMDF_API xid xfReadVectorValuesAtIndicesFloat (xid a_Id, 
                    int a_nIndices, const int *a_Indices, int a_FirstTime,
                    int a_NumTimes, int a_NumComponents, float *a_Values)
{
  return xfReadVectorValuesAtIndices(a_Id, a_nIndices, a_Indices, a_FirstTime,
                                     a_NumTimes, a_NumComponents, a_Values);
}
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfScalarDataLocation */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfScalarDataLocation (xid a_Id, int a_DataLoc)
{
  int     error = 1;

  error = xfpWriteAttributeInt(a_Id, DATASET_ATT_DATALOCATION, 1, &a_DataLoc);

  return error;
} /* xfScalarDataLocation */
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfVector2DDataLocations */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfVector2DDataLocations (xid a_Id, int a_DataLocI, int a_DataLocJ)
{
  int     error = 1;

  error = xfpWriteAttributeInt(a_Id, DATASET_ATT_DATALOCATIONI, 1, &a_DataLocI);
  if (error >= 0) {
    error = xfpWriteAttributeInt(a_Id, DATASET_ATT_DATALOCATIONJ, 1, &a_DataLocJ);
  }

  return error;
} /* xfVector2DDataLocations */
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfVector3DDataLocations */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfVector3DDataLocations (xid a_Id, int a_DataLocI, int a_DataLocJ,
                                      int a_DataLocK)
{
  int     error = 1;

  error = xfpWriteAttributeInt(a_Id, DATASET_ATT_DATALOCATIONI, 1, &a_DataLocI);
  if (error >= 0) {
    error = xfpWriteAttributeInt(a_Id, DATASET_ATT_DATALOCATIONJ, 1, &a_DataLocJ);
    if (error >= 0) {
      error = xfpWriteAttributeInt(a_Id, DATASET_ATT_DATALOCATIONK, 1, &a_DataLocK);
    }
  }

  return error;
} /* xfVector3DDataLocations */
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfGetScalarDataLocation */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfGetScalarDataLocation (xid a_Id, int *a_DataLoc)
{
  int error = 1;

  error = xfpReadAttributeInt(a_Id, DATASET_ATT_DATALOCATION, 1, a_DataLoc);

  return error;
} /* xfGetScalarDataLocation */
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfGetVector2DDataLocations */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfGetVector2DDataLocations (xid a_Id, int *a_DataLocI,
                               int *a_DataLocJ)
{
  int error = 1;

  error = xfpReadAttributeInt(a_Id, DATASET_ATT_DATALOCATIONI, 1, a_DataLocI);
  if (error >= 0) {
    error = xfpReadAttributeInt(a_Id, DATASET_ATT_DATALOCATIONJ, 1, a_DataLocJ);
  }

  return error;
} /* xfGetVector2DDataLocations */
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfGetVector3DDataLocations */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfGetVector3DDataLocations (xid a_Id, int *a_DataLocI,
                               int *a_DataLocJ, int *a_DataLocK)
{
  int error = 1;

  error = xfpReadAttributeInt(a_Id, DATASET_ATT_DATALOCATIONI, 1, a_DataLocI);
  if (error >= 0) {
    error = xfpReadAttributeInt(a_Id, DATASET_ATT_DATALOCATIONJ, 1, a_DataLocJ);
    if (error >= 0) {
      error = xfpReadAttributeInt(a_Id, DATASET_ATT_DATALOCATIONK, 1, a_DataLocK);
    }
  }

  return error;
} /* xfGetVector3DDataLocations */
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfVectorsInLocalCoords */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfVectorsInLocalCoords (xid a_Id)
{
  int     error = 1, on = 1;

  error = xfpWriteAttributeInt(a_Id, DATASET_ATT_LOCALCOORDS, 1, &on);
  
  return error;
} /* xfVectorsInLocalCoords */
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfAreVectorsInLocalCoords */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfAreVectorsInLocalCoords (xid a_Id, int *a_LocalCoords)
{
  int     error = 1, on = 1;

  *a_LocalCoords = 0;
  error = xfpReadAttributeInt(a_Id, DATASET_ATT_LOCALCOORDS, 1, &on);
  if (error >= 0 && on == 1) {
    *a_LocalCoords = 1;
  }
  
  return error;
} /* xfAreVectorsInLocalCoords */

/* --------------------------------------------------------------------------- */
/* FUNCTION  xfGetCoordVersion */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfGetCoordVersion (xid a_CoordId, int *a_version)
{
  int    error = 1;

  error = xfpReadAttributeInt(a_CoordId, COORD_ATT_VERSION, 1, a_version);
  if (error < 0) {
    *a_version = 1;
    error = 0;
  }

  return error;
} /* xfGetCoordVersion */
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfGetUsingLocalCoords */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfGetUsingLocalCoords (xid a_CoordId, xmbool *a_localCoords)
{
  int    error = 1;
  
  int ival = 0;
  error = xfpReadAttributeInt(a_CoordId, COORD_ATT_LOCAL, 1, &ival);

  *a_localCoords = (ival == 1);

  return error;
} /* xfGetUsingLocalCoords */
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfSetUsingLocalCoords */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfSetUsingLocalCoords (xid a_CoordId, xmbool *a_localCoords)
{
  int    error = 1;

  int ival = 0;
  if (*a_localCoords) {
    ival = 1;
  }
  error = xfpWriteAttributeInt(a_CoordId, COORD_ATT_LOCAL, 1, &ival);

  return error;
} /* xfSetUsingLocalCoords */

/* --------------------------------------------------------------------------- */
/* FUNCTION  xfSetWKT */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfSetWKT (xid a_CoordId, const char * a_wkt)
{
  int    error = 1;

  error = xfpWriteAttributeString(a_CoordId, COORD_ATT_WKT, a_wkt);

  return error;
} /* xfSetWKT */

/* --------------------------------------------------------------------------- */
/* FUNCTION  xfSetAttributeString */
/*! PURPOSE:
 * - NOTES: Added to write of CF compliance attributes to objects.   */
 /* --------------------------------------------------------------------------- */
XMDF_API xid xfSetAttributeString(xid a_Id, const char *a_Name, const char *a_String)
{
    int    error = 1;
    
    error = xfpWriteAttributeString(a_Id, a_Name, a_String);  /* set the name attribute to the value in a_str */

    return error;
} /* xfSetAttributeString */

/* --------------------------------------------------------------------------- */
/* FUNCTION  xfGetWKTStringSize */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfGetWKTStringSize (xid a_CoordId, int *a_size)
{
  int    error = 1;

  char *wkt = NULL;
  error = xfpReadAttributeString(a_CoordId, COORD_ATT_WKT, 1, &wkt);
  if (wkt != NULL)
  {
    *a_size = strlen(wkt);
    free(wkt);
  }

  return error;
} /* xfGetWKTStringSize */

/* --------------------------------------------------------------------------- */
/* FUNCTION  xfGetWKT */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfGetWKT (xid a_CoordId, char *a_wkt)
{
  int    error = 1;

  char *wkt = NULL;
  error = xfpReadAttributeString(a_CoordId, COORD_ATT_WKT, 1, &wkt);
  strcpy(a_wkt, wkt);
  free(wkt);

  return error;
} /* xfGetWKT */

/* --------------------------------------------------------------------------- */
/* FUNCTION  xfGetHorizDatum */
/*! PURPOSE:   Read the horizontal datum (\#defines in xmdf.h)
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfGetHorizDatum (xid a_CoordId, int *a_val)
{
  int    error = 1;

  error = xfpReadAttributeInt(a_CoordId, COORD_ATT_HORIZ_DATUM, 1, a_val);

  return error;
} /* xfGetHorizDatum */
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfGetHorizUnits */
/*! PURPOSE:   Read the horizontal units (\#defines in xmdf.h)
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfGetHorizUnits (xid a_CoordId, int *a_val)
{
  int    error = 1;

  error = xfpReadAttributeInt(a_CoordId, COORD_ATT_HORIZ_UNITS, 1, a_val);

  return error;
} /* xfGetHorizUnits */
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfGetVertDatum */
/*! PURPOSE:   Read the vertical datum (\#defines in xmdf.h)
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfGetVertDatum (xid a_CoordId, int *a_val)
{
  int    error = 1;

  error = xfpReadAttributeInt(a_CoordId, COORD_ATT_VERT_DATUM, 1, a_val);

  return error;
} /* xfGetVertDatum */
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfGetVertUnits */
/*! PURPOSE:   Read the vertical units (\#defines in xmdf.h)
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfGetVertUnits (xid a_CoordId, int *a_val)
{
  int    error = 1;

  error = xfpReadAttributeInt(a_CoordId, COORD_ATT_VERT_UNITS, 1, a_val);

  return error;
} /* xfGetVertUnits */
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfGetLat */
/*! PURPOSE:   Read whether the lattitude is North or South latitude
 *           (\#defines in xmdf.h)

 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfGetLat (xid a_CoordId, int *a_val)
{
  int    error = 1;

  error = xfpReadAttributeInt(a_CoordId, COORD_ATT_LATITUDE, 1, a_val);

  return error;
} /* xfGetLat */
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfGetLon */
/*! PURPOSE:   Read whether the longitude is East or West
 *           (\#defines in xmdf.h)
 */
/* PURPOSE
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfGetLon (xid a_CoordId, int *a_val)
{
  int    error = 1;

  error = xfpReadAttributeInt(a_CoordId, COORD_ATT_LONGITUDE, 1, a_val);

  return error;
} /* xfGetLon */
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfGetUTMZone */
/*! PURPOSE:   Read the UTM zone (should be a number between 1 and 60)
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfGetUTMZone (xid a_CoordId, int *a_val)
{
  int    error = 1;

  error = xfpReadAttributeInt(a_CoordId, COORD_ATT_UTM_ZONE, 1, a_val);

  return error;
} /* xfGetUTMZone */
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfGetSPCZone */
/*! PURPOSE:   Read the SPC zone (Lookup numbers in XMDF documentation)
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfGetSPCZone (xid a_CoordId, int *a_val)
{
  int    error = 1;

  error = xfpReadAttributeInt(a_CoordId, COORD_ATT_SPC_ZONE, 1, a_val);

  return error;
} /* xfGetSPCZone */
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfGetHPGNArea */
/*! PURPOSE:   Read the HPGN area (Lookup numbers in XMDF documentation)
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfGetHPGNArea (xid a_CoordId, int *a_val)
{
  int    error = 1;

  error = xfpReadAttributeInt(a_CoordId, COORD_ATT_HPGN_AREA, 1, a_val);

  return error;
} /* xfGetHPGNArea */
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfGetCPPLat */
/*! PURPOSE:   Read the Carte Prallelo Grammatique Projetion Factor for
 *           converting the latitude

 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfGetCPPLat (xid a_CoordId, double *a_val)
{
  int    error = 1;

  error = xfpReadAttributeDouble(a_CoordId, COORD_ATT_CPP_LATITUDE, 1, a_val);

  return error;
} /* xfGetCPPLat */
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfGetCPPLon */
/*! PURPOSE:   Read the Carte Prallelo Grammatique Projetion Factor for
 *  converting the longitude

 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfGetCPPLon (xid a_CoordId, double *a_val)
{
  int    error = 1;

  error = xfpReadAttributeDouble(a_CoordId, COORD_ATT_CPP_LONGITUDE, 1, a_val);

  return error;
} /* xfGetCPPLon */
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfGetEllipse */
/*! PURPOSE:   Read the Ellipse number based upon XMDF documentation
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfGetEllipse (xid a_CoordId, int *a_val)
{
  int    error = 1;

  error = xfpReadAttributeInt(a_CoordId, COORD_ATT_ELLIPSE, 1, a_val);

  return error;
} /* xfGetEllipse */
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfGetMajorR */
/*! PURPOSE:   Read a user-defined ellipse major radius
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfGetMajorR (xid a_CoordId, double *a_val)
{
  int    error = 1;

  error = xfpReadAttributeDouble(a_CoordId, COORD_ATT_MAJOR_R, 1, a_val);

  return error;
} /* xfGetMajorR */
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfGetMinorR */
/*! PURPOSE:   Read a user-defined ellipse minor radius
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfGetMinorR (xid a_CoordId, double *a_val)
{
  int    error = 1;

  error = xfpReadAttributeDouble(a_CoordId, COORD_ATT_MINOR_R, 1, a_val);

  return error;
} /* xfGetMinorR */

/* --------------------------------------------------------------------------- */
/* FUNCTION  xfSetHorizDatum */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfSetHorizDatum (xid a_CoordId, int a_val)
{
  int    error = 1;

  error = xfpWriteAttributeInt(a_CoordId, COORD_ATT_HORIZ_DATUM, 1, &a_val);

  return error;
} /* xfSetHorizDatum */
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfSetHorizUnits */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfSetHorizUnits (xid a_CoordId, int a_val)
{
  int    error = 1;

  error = xfpWriteAttributeInt(a_CoordId, COORD_ATT_HORIZ_UNITS, 1, &a_val);

  return error;
} /* xfSetHorizUnits */
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfSetVertDatum */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfSetVertDatum (xid a_CoordId, int a_val)
{
  int    error = 1;

  error = xfpWriteAttributeInt(a_CoordId, COORD_ATT_VERT_DATUM, 1, &a_val);

  return error;
} /* xfSetVertDatum */
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfSetVertUnits */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfSetVertUnits (xid a_CoordId, int a_val)
{
  int    error = 1;

  error = xfpWriteAttributeInt(a_CoordId, COORD_ATT_VERT_UNITS, 1, &a_val);

  return error;
} /* xfSetVertUnits */
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfSetLat */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfSetLat (xid a_CoordId, int a_val)
{
  int    error = 1;

  error = xfpWriteAttributeInt(a_CoordId, COORD_ATT_LATITUDE, 1, &a_val);

  return error;
} /* xfSetLat */
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfSetLon */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfSetLon (xid a_CoordId, int a_val)
{
  int    error = 1;

  error = xfpWriteAttributeInt(a_CoordId, COORD_ATT_LONGITUDE, 1, &a_val);

  return error;
} /* xfSetLon */
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfSetUTMZone */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfSetUTMZone (xid a_CoordId, int a_val)
{
  int    error = 1;

  error = xfpWriteAttributeInt(a_CoordId, COORD_ATT_UTM_ZONE, 1, &a_val);

  return error;
} /* xfSetUTMZone */
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfSetSPCZone */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfSetSPCZone (xid a_CoordId, int a_val)
{
  int    error = 1;

  error = xfpWriteAttributeInt(a_CoordId, COORD_ATT_SPC_ZONE, 1, &a_val);

  return error;
} /* xfSetSPCZone */
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfSetHPGNArea */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfSetHPGNArea (xid a_CoordId, int a_val)
{
  int    error = 1;

  error = xfpWriteAttributeInt(a_CoordId, COORD_ATT_HPGN_AREA, 1, &a_val);

  return error;
} /* xfSetHPGNArea */
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfSetCPPLat */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfSetCPPLat (xid a_CoordId, double a_val)
{
  int    error = 1;

  error = xfpWriteAttributeDouble(a_CoordId, COORD_ATT_CPP_LATITUDE, 1, &a_val);

  return error;
} /* xfSetCPPLat */
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfSetCPPLon */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfSetCPPLon (xid a_CoordId, double a_val)
{
  int    error = 1;

  error = xfpWriteAttributeDouble(a_CoordId, COORD_ATT_CPP_LONGITUDE, 1, &a_val);

  return error;
} /* xfSetCPPLon */
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfSetEllipse */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfSetEllipse (xid a_CoordId, int a_val)
{
  int    error = 1;

  error = xfpWriteAttributeInt(a_CoordId, COORD_ATT_ELLIPSE, 1, &a_val);

  return error;
} /* xfSetEllipse */
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfSetMajorR */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfSetMajorR (xid a_CoordId, double a_val)
{
  int    error = 1;

  error = xfpWriteAttributeDouble(a_CoordId, COORD_ATT_MAJOR_R, 1, &a_val);

  return error;
} /* xfSetMajorR */
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfSetMinorR */
/*! PURPOSE:
 * - NOTES:      */
/* --------------------------------------------------------------------------- */
XMDF_API xid xfSetMinorR (xid a_CoordId, double a_val)
{
  int    error = 1;

  error = xfpWriteAttributeDouble(a_CoordId, COORD_ATT_MINOR_R, 1, &a_val);

  return error;
} /* xfSetMinorR */
/******************************************************************************/
/*   FUNCTION  xfWriteReftime */
/*! PURPOSE:   Set a reference time for a dataset
 *   - NOTES:     */
/******************************************************************************/
XMDF_API xid xfWriteReftime (xid a_Id, double a_Reftime)
{
  int    status;

  status = xfpWriteAttributeDouble(a_Id, DATASET_REFTIME, 1, &a_Reftime);

  return status;
} /* xfWriteReftime */
/******************************************************************************/
/*   FUNCTION  xfUseReftime */
/*! PURPOSE:   See if a reference time exists for a dataset
 *   - NOTES:     */
/******************************************************************************/
XMDF_API xid xfUseReftime (xid a_Id, xmbool *a_bUseReftime)
{
  int    status;

  status = xfpDoesAttributeExist(a_Id, DATASET_REFTIME, a_bUseReftime);

  return status;
} /* xfUseReftime */
/******************************************************************************/
/*   FUNCTION  xfReadReftime */
/*! PURPOSE:   Read a reference time for a dataset
 *   - NOTES:     */
/******************************************************************************/
XMDF_API xid xfReadReftime (xid a_Id, double *a_dReftime)
{
  int    status;

  status = xfpReadAttributeDouble(a_Id, DATASET_REFTIME, 1, a_dReftime);

  return status;
} /* xfReadReftime */
/*----------------------------------------------------------------------------- */
/* FUNCTION  xfCalendarToJulian */
/*! PURPOSE:   Convert Calendar to Julian Date
 * - NOTES:     Returns 1 if successful and -1 if unsuccessful. 
 *           (Julian day number algorithm adopted from Press et al.)
 *           era = 0 represents BCE (BC), and era = 1 represents CE (AD).
 *          \verbatim
 *          -Taken from JavaScript Code found at website:
 *                http://aa.usno.navy.mil/data/docs/JulianDate.html
 *           Contact info. provided by website:
 *             Marc A. Murison
 *             Astronomical Applications Dept.
 *             U.S. Naval Observatory
 *             3450 Massachusetts Ave, NW
 *             Washington, DC  20392-5420    \endverbatim */
/*---------------------------------------------------------------------------- */
XMDF_API xid xfCalendarToJulian(xmbool a_bEra, int a_yr, int a_mo, int a_day, 
                       int a_hr, int a_min, int a_sec, double *a_julian) 
{
  double  jy, ja, jm, /*jd, jd0,*/ tmpyr; 
  double  intgr, gregcal, dayfrac, frac;

  if (a_yr==0) {
      /* There is no year 0 in the Julian system! */
    return -1;
  }
  if ((a_yr==1582) && (a_mo==10) && (a_day>4) && 
                      (a_day<15) && (a_bEra == ERA_IS_BCE)) {
      /* The dates 5 through 14 October, 1582 */
        /* do not exist in the Gregorian system! */
    return -1;
  }

  tmpyr = a_yr;

  if (a_bEra == ERA_IS_BCE) {
    tmpyr = 1 - tmpyr;
  }
  if (a_mo > 2) {
    jy = tmpyr;
    jm = a_mo + 1;
  }
  else {
    jy = tmpyr - 1;
    jm = a_mo + 13;
  }

  intgr = floor(floor(365.25*jy) + floor(30.6001*jm) + a_day + 1720995);

    /* check for switch to Gregorian calendar */
  gregcal = 15 + 31*(10 + 12*1582);
  if (a_day + 31*(a_mo + 12*tmpyr) >= gregcal) {
    ja = floor(0.01*jy);
    intgr += 2 - ja + floor(0.25*ja);
  }

    /* correct for half-day offset */
  dayfrac = a_hr/24.0 - 0.5;
  if (dayfrac < 0.0) {
    dayfrac += 1.0;
    --intgr;
  }

    /* now set the fraction of a day */
  frac = dayfrac + (a_min + a_sec/60.0)/60.0/24.0;


    /* round to nearest second */  
   // we don't want to round


  *a_julian = intgr + frac;
  return 1;
} /* xfCalendarToJulian */

/*---------------------------------------------------------------------------- */
/* FUNCTION  xfJulianToCalendar */
/*! PURPOSE:   Convert Julian Date To Calendar
 * - NOTES:     Returns 1 if successful and -1 if unsuccessful. 
 *           (algorithm adopted from Press et al.)
 *           era = 0 represents BCE (BC), and era = 1 represents CE (AD).
 *          \verbatim
 *          -Taken from JavaScript Code found at website:
 *                http://aa.usno.navy.mil/data/docs/JulianDate.html
 *           Contact info. provided by website:
 *             Marc A. Murison
 *             Astronomical Applications Dept.
 *             U.S. Naval Observatory
 *             3450 Massachusetts Ave, NW
 *             Washington, DC  20392-5420  \endverbatim */
/*---------------------------------------------------------------------------- */
XMDF_API xid xfJulianToCalendar (xmbool *a_bEra, int *a_yr, int *a_mo, int *a_day,
                        int *a_hr, int *a_min, int *a_sec, double a_julian)
{
  double  j1, j2, j3, j4, j5;
  double  intgr, frac, gregjd, tmp;
  double  dayfrac, f;

    /* get the date from the Julian day number */
  intgr  = floor(a_julian);
  frac   = a_julian - intgr;
  gregjd = 2299161;
    /* Gregorian Calendar Correction */
  if (intgr >= gregjd) {
    tmp = floor(((intgr-1867216) - 0.25) / 36524.25);
    j1 = intgr + 1 + tmp - floor(0.25*tmp);
  }
  else {
    j1 = intgr;
  }

    /*correction for half day offset */
  dayfrac = frac + 0.5;
  if (dayfrac >= 1.0) {
    dayfrac -= 1.0;
    ++j1;
  }

  j2 = j1 + 1524;
  j3 = floor(6680.0 + ((j2 - 2439870) - 122.1)/365.25);
  j4 = floor(j3*365.25);
  j5 = floor((j2-j4)/30.6001);

  *a_day = (int)floor(j2 - j4 - floor(j5*30.6001));
  *a_mo  = (int)floor(j5 - 1);
  if(*a_mo > 12 ) {
    *a_mo -= 12;
  }
  *a_yr = (int)floor(j3 - 4715);
  if (*a_mo > 2) {
    --*a_yr;
  }
  if (*a_yr <= 0) {
    --*a_yr;
  }
  
    /* get time of day from day fraction */
  *a_hr  = (int)floor(dayfrac * 24.0);
  *a_min = (int)floor((dayfrac*24.0 - *a_hr)*60.0);
  f   = ((dayfrac*24.0 - *a_hr)*60.0 - *a_min)*60.0;
  *a_sec = (int)floor(f);
  f   -= *a_sec;
  if (f>0.5) {
    ++*a_sec;
  }
  if (*a_sec==60) {
    *a_sec = 0;
    ++*a_min;
  }
  if (*a_min == 60)  {
    *a_min = 0;
    ++*a_hr;
  }
  if (*a_hr == 24) {
    *a_hr = 0;
    ++*a_day;
  }

  if (*a_yr < 0) {
    *a_yr = -*a_yr;
    *a_bEra = XFALSE;
  }
  else {
    *a_bEra = XTRUE;
  }

  return 1;
} /* xfJulianToCalendar */
/******************************************************************************/
/*   FUNCTION  xfGetParallelNumValuesToRead */
/*! PURPOSE:   This function is to return a number of values to read in for a
 *             datset at a time so that we dont use all of the RAM
 *   - NOTES:      */
/******************************************************************************/
 XMDF_API hssize_t xfGetParallelNumValuesToRead ()
{
  return xfpGetParallelNumValuesToRead();

} /* xfGetParallelNumValuesToRead */
/******************************************************************************/
/*   FUNCTION  xfSetParallelNumValuesToRead */
/*! PURPOSE:   This function allows you to set a number of values to read in for a
 *             datset at a time so that we dont use all of the RAM
 *   - NOTES:      */
/******************************************************************************/
 XMDF_API void xfSetParallelNumValuesToRead (hssize_t a)
{
  xfpSetParallelNumValuesToRead(a);

} /* xfSetParallelNumValuesToRead */
