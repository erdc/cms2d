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

/*****************************************************************************
 * file       Xmdf_f.c.cpp
 *
 * purpose    Used as a wrapper to call C functions from the FORTRAN code in XMDF
 *
 * notes  
 *
 * coded by
 *
 ******************************************************************************/


#include "xmdf/Xmdf.h"
#include "xmdf/xmdf_private.h"
#include <stdlib.h>
#include <string.h>

#define DBG Dbg(__FILE__,__LINE__);
extern void Dbg();

/*-----------------------------------------------------------------------------*/
/* FUNCTION xfpInitialize_f_                                                   */
/* PURPOSE                                                                     */
/* NOTES                                                                       */
/*                                                                             */
/*-----------------------------------------------------------------------------*/
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER
#define func_name(XFPINITIALIZE_F,xfpinitialize_f,xfpinitialize_f_) XFPINITIALIZE_F
#elif defined FLOWER
#define func_name(XFPINITIALIZE_F,xfpinitialize_f,xfpinitialize_f_) xfpinitialize_f
#else
#define func_name(XFPINITIALIZE_F,xfpinitialize_f,xfpinitialize_f_) xfpinitialize_f_
#endif
XMDF_API  int func_name(XFPINITIALIZE_F,xfpinitialize_f,xfpinitialize_f_) ()
{
  /*In FORTRAN the library must be initialized*/
  H5open();

  /*call C version*/
  xfpInitialize();

  return 0;
 }/* xfpInitialize_f_*/
 /*-----------------------------------------------------------------------------*/
/* FUNCTION xfpclose_f_                                                        */
/* PURPOSE                                                                      */
/* NOTES                                                                        */
/*                                                                              */
/*----------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER
#define func_name(XFPCLOSE_F,xfpclose_f,xfpclose_f_) XFPCLOSE_F
#elif defined FLOWER
#define func_name(XFPCLOSE_F,xfpclose_f,xfpclose_f_) xfpclose_f
#else
#define func_name(XFPCLOSE_F,xfpclose_f,xfpclose_f_) xfpclose_f_
#endif
XMDF_API  int func_name(XFPCLOSE_F,xfpclose_f,xfpclose_f_) ()
{
  int nOpen = 0;
  xfGetNumOpenIdentifiers(H5F_OBJ_ALL, &nOpen);

  if(nOpen > 0)
  {
    IdentifierInfo* ii = malloc(sizeof(IdentifierInfo)*nOpen);

    xfGetOpenIdentifiersInfo(H5F_OBJ_ALL, nOpen, ii);

    free(ii);
  }

  /*In FORTRAN the library must be closed*/
  H5close();

  return 0;

 }/* xfpClose_f_*/
/*-----------------------------------------------------------------------------*/
/* FUNCTION  xfCreateFile_f_*/
/* PURPOSE   Creates a file to read using the Xmdf library.*/
/* NOTES     Calls C version*/
/*-----------------------------------------------------------------------------*/
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFCREATEFILE_F,xfcreatefile_f,xfcreatefile_f_) XFCREATEFILE_F
#elif defined FLOWER
#define func_name(XFCREATEFILE_F,xfcreatefile_f,xfcreatefile_f_) xfcreatefile_f
#else
#define func_name(XFCREATEFILE_F,xfcreatefile_f,xfcreatefile_f_) xfcreatefile_f_
#endif
XMDF_API  int func_name(XFCREATEFILE_F,xfcreatefile_f,xfcreatefile_f_) (const char * a_Filename,
                                                              int* a_namelen,
                                                              xid *a_Id,
                                                              short *a_Overwrite)
{

  xmbool Overwrite; 
  int	error;
  char* filename_copy = NULL;

  Overwrite = (*a_Overwrite) ? XTRUE : XFALSE; 

  filename_copy = (char*)malloc((*a_namelen+1)*sizeof(char));
  strncpy(filename_copy, a_Filename, *a_namelen);
  filename_copy[*a_namelen] = '\0';
  
  error=xfCreateFile(filename_copy, a_Id, Overwrite);

  free(filename_copy);
  return error;
} /* xfCreateFile_f_*/
/*-----------------------------------------------------------------------------*/
/* FUNCTION  xfOpenFile_f_*/
/* PURPOSE   Opens a file to read using the Xmdf library.*/
/* NOTES     Calls C version*/
/*-----------------------------------------------------------------------------*/
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFOPENFILE_F,xfopenfile_f,xfopenfile_f_) XFOPENFILE_F
#elif defined FLOWER
#define func_name(XFOPENFILE_F,xfopenfile_f,xfopenfile_f_) xfopenfile_f
#else
#define func_name(XFOPENFILE_F,xfopenfile_f,xfopenfile_f_) xfopenfile_f_
#endif
XMDF_API  int func_name(XFOPENFILE_F,xfopenfile_f,xfopenfile_f_) (const char * a_Filename,
                                                        int* a_namelen, xid *a_Id,
                                                        short *a_ReadOnly)
{

  xmbool ReadOnly;
  int error;
  char * filename_copy = NULL;

  ReadOnly = (*a_ReadOnly) ? XTRUE : XFALSE;

  filename_copy = (char*)malloc((*a_namelen+1)*sizeof(char));
  strncpy(filename_copy, a_Filename, *a_namelen);
  filename_copy[*a_namelen] = '\0';

  error=xfOpenFile(filename_copy, a_Id, ReadOnly);

  free(filename_copy);
  return error;
} /* xfOpenFile_f_*/
/*-----------------------------------------------------------------------------*/
/* FUNCTION  xfclosefile_f_*/
/* PURPOSE   Opens a file to read using the Xmdf library.*/
/* NOTES     Calls C version*/
/*-----------------------------------------------------------------------------*/
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFCLOSEFILE_F,xfclosefile_f,xfclosefile_f_) XFCLOSEFILE_F
#elif defined FLOWER
#define func_name(XFCLOSEFILE_F,xfclosefile_f,xfclosefile_f_) xfclosefile_f
#else
#define func_name(XFCLOSEFILE_F,xfclosefile_f,xfclosefile_f_) xfclosefile_f_
#endif
XMDF_API  int func_name(XFCLOSEFILE_F,xfclosefile_f,xfclosefile_f_) (hid_t *a_Id)
{

  int error;
  error=xfCloseFile(*a_Id);

  return error;
} /* xfclosefile_f_*/
/*-----------------------------------------------------------------------------
FUNCTION  xfGetLibraryVersionFile_f_
PURPOSE   Obtain the version of XMDF from the library file
NOTES     Calls C version
------------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFGETLIBRARYVERSIONFILE_F,xfgetlibraryversionfile_f,xfgetlibraryversionfile_f_) XFGETLIBRARYVERSIONFILE_F
#elif defined FLOWER
#define func_name(XFGETLIBRARYVERSIONFILE_F,xfgetlibraryversionfile_f,xfgetlibraryversionfile_f_) xfgetlibraryversionfile_f
#else
#define func_name(XFGETLIBRARYVERSIONFILE_F,xfgetlibraryversionfile_f,xfgetlibraryversionfile_f_) xfgetlibraryversionfile_f_
#endif
XMDF_API  int func_name(XFGETLIBRARYVERSIONFILE_F,xfgetlibraryversionfile_f,xfgetlibraryversionfile_f_) 
           (xid* a_File, float *a_Version)
{
  int  status;

  status = xfGetLibraryVersionFile(*a_File, a_Version);

  return status;
} /* xfGetLibraryVersionFile */
/*-----------------------------------------------------------------------------
FUNCTION  xfgetlibraryversion_f_
PURPOSE   Obtain the version of XMDF from the library file
NOTES     Calls C version
------------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFGETLIBRARYVERSION_F,xfgetlibraryversion_f,xfgetlibraryversion_f_) XFGETLIBRARYVERSION_F
#elif defined FLOWER
#define func_name(XFGETLIBRARYVERSION_F,xfgetlibraryversion_f,xfgetlibraryversion_f_) xfgetlibraryversion_f
#else
#define func_name(XFGETLIBRARYVERSION_F,xfgetlibraryversion_f,xfgetlibraryversion_f_) xfgetlibraryversion_f_
#endif
XMDF_API  int func_name(XFGETLIBRARYVERSION_F,xfgetlibraryversion_f,xfgetlibraryversion_f_)
              (float *a_Version)
{
  int  status;

  status = xfGetLibraryVersion(a_Version);

  return status;
} /* xfgetlibraryversion_f_*/
/*-----------------------------------------------------------------------------*/
/* FUNCTION  xfCreateCoordinateGroup_f_*/
/* PURPOSE   Create a Vector dataset using the Xmdf library.*/
/* NOTES     Calls C version*/
/*-----------------------------------------------------------------------------*/
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFCREATECOORDINATEGROUP_F,xfcreatecoordinategroup_f,xfcreatecoordinategroup_f_) XFCREATECOORDINATEGROUP_F
#elif defined FLOWER
#define func_name(XFCREATECOORDINATEGROUP_F,xfcreatecoordinategroup_f,xfcreatecoordinategroup_f_) xfcreatecoordinategroup_f
#else
#define func_name(XFCREATECOORDINATEGROUP_F,xfcreatecoordinategroup_f,xfcreatecoordinategroup_f_) xfcreatecoordinategroup_f_
#endif
XMDF_API  int func_name(XFCREATECOORDINATEGROUP_F,xfcreatecoordinategroup_f,xfcreatecoordinategroup_f_)
              (xid *a_ParentId, xid *a_ChildId)
{

  int    error;

  error=xfCreateCoordinateGroup(*a_ParentId, a_ChildId);

  return error;
} /* xfCreateCoordinateGroup_f_*/
/*-----------------------------------------------------------------------------*/
/* FUNCTION  xfOpenCoordinateGroup_f_*/
/* PURPOSE   Create a Vector dataset using the Xmdf library.*/
/* NOTES     Calls C version*/
/*-----------------------------------------------------------------------------*/
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFOPENCOORDINATEGROUP_F,xfopencoordinategroup_f,xfopencoordinategroup_f_) XFOPENCOORDINATEGROUP_F
#elif defined FLOWER
#define func_name(XFOPENCOORDINATEGROUP_F,xfopencoordinategroup_f,xfopencoordinategroup_f_) xfopencoordinategroup_f
#else
#define func_name(XFOPENCOORDINATEGROUP_F,xfopencoordinategroup_f,xfopencoordinategroup_f_) xfopencoordinategroup_f_
#endif
XMDF_API  int func_name(XFOPENCOORDINATEGROUP_F,xfopencoordinategroup_f,xfopencoordinategroup_f_)
              (xid *a_ParentId, xid *a_ChildId)
{

  int    error;

  error=xfOpenCoordinateGroup(*a_ParentId, a_ChildId);

  return error;
} /* xfOpenCoordinateGroup_f_*/

/*-----------------------------------------------------------------------------*/
/* FUNCTION  xfCreatePropertyGroup_f_*/
/* PURPOSE   Create a Vector dataset using the Xmdf library.*/
/* NOTES     Calls C version*/
/*-----------------------------------------------------------------------------*/
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFCREATEPROPERTYGROUP_F,xfcreatepropertygroup_f,xfcreatepropertygroup_f_) XFCREATEPROPERTYGROUP_F
#elif defined FLOWER
#define func_name(XFCREATEPROPERTYGROUP_F,xfcreatepropertygroup_f,xfcreatepropertygroup_f_) xfcreatepropertygroup_f
#else
#define func_name(XFCREATEPROPERTYGROUP_F,xfcreatepropertygroup_f,xfcreatepropertygroup_f_) xfcreatepropertygroup_f_
#endif
XMDF_API  int func_name(XFCREATEPROPERTYGROUP_F,xfcreatepropertygroup_f,xfcreatepropertygroup_f_)
              (xid *a_ParentId, xid *Id)
{

  int    error;

  error=xfCreatePropertyGroup(*a_ParentId, Id);

  return error;

} /* xfCreatePropertyGroup_f_*/
 /* --------------------------------------------------------------------------- */
/* FUNCTION  xfWritePropertyInt */
/* PURPOSE   */
/* NOTES     */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFWRITEPROPERTYINT_F,xfwritepropertyint_f,xfwritepropertyint_f_) XFWRITEPROPERTYINT_F
#elif defined FLOWER
#define func_name(XFWRITEPROPERTYINT_F,xfwritepropertyint_f,xfwritepropertyint_f_) xfwritepropertyint_f
#else
#define func_name(XFWRITEPROPERTYINT_F,xfwritepropertyint_f,xfwritepropertyint_f_) xfwritepropertyint_f_
#endif
XMDF_API  int func_name(XFWRITEPROPERTYINT_F,xfwritepropertyint_f,xfwritepropertyint_f_)
              (xid *a_Id, const char *a_Name, int *a_namelen, 
                          int *a_Number, const int *a_Properties, 
                          int *a_Compression)
 {

  int	error;
  char* name_copy = NULL;

  name_copy = (char*)malloc((*a_namelen+1)*sizeof(char));
  strncpy(name_copy, a_Name, *a_namelen);
  name_copy[*a_namelen] = '\0';

  error = xfWritePropertyInt(*a_Id, name_copy, *a_Number, a_Properties, *a_Compression);

  free(name_copy);
  return error;

 } /*xfWritePropertyInt_f_*/
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfWritePropertyUnsignedInt_f_ */
/* PURPOSE   */
/* NOTES     */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFWRITEPROPERTYUNSIGNEDINT_F,xfwritepropertyunsignedint_f,xfwritepropertyunsignedint_f_) XFWRITEPROPERTYUNSIGNEDINT_F
#elif defined FLOWER
#define func_name(XFWRITEPROPERTYUNSIGNEDINT_F,xfwritepropertyunsignedint_f,xfwritepropertyunsignedint_f_) xfwritepropertyunsignedint_f
#else
#define func_name(XFWRITEPROPERTYUNSIGNEDINT_F,xfwritepropertyunsignedint_f,xfwritepropertyunsignedint_f_) xfwritepropertyunsignedint_f_
#endif
XMDF_API  int func_name(XFWRITEPROPERTYUNSIGNEDINT_F,xfwritepropertyunsignedint_f,xfwritepropertyunsignedint_f_)
              (xid *a_Id, const char *a_Name, int *a_length, 
                                  int *a_Number, const unsigned int *a_Properties, 
                                  int *a_Compression)
 {

  int	error;
  char* name_copy = NULL;

  name_copy = (char*)malloc((*a_length+1)*sizeof(char));
  strncpy(name_copy, a_Name, *a_length);
  name_copy[*a_length] = '\0';

  error = xfWritePropertyUnsignedInt(*a_Id, name_copy, *a_Number, a_Properties, 
                                     *a_Compression);

  free(name_copy);
  return error;
 } /*xfWritePropertyUnsignedInt_f_*/
 /* --------------------------------------------------------------------------- */
/* FUNCTION  xfWritePropertyDouble_f_ */
/* PURPOSE   */
/* NOTES     */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFWRITEPROPERTYDOUBLE_F,xfwritepropertydouble_f,xfwritepropertydouble_f_) XFWRITEPROPERTYDOUBLE_F
#elif defined FLOWER
#define func_name(XFWRITEPROPERTYDOUBLE_F,xfwritepropertydouble_f,xfwritepropertydouble_f_) xfwritepropertydouble_f
#else
#define func_name(XFWRITEPROPERTYDOUBLE_F,xfwritepropertydouble_f,xfwritepropertydouble_f_) xfwritepropertydouble_f_
#endif
XMDF_API  int func_name(XFWRITEPROPERTYDOUBLE_F,xfwritepropertydouble_f,xfwritepropertydouble_f_)
              (xid *a_Id, const char *a_Name, int *a_length, 
                             int *a_Number, const double *a_Properties, 
                             int *a_Compression)
 {

  int	error;
  char* name_copy = NULL;

  name_copy = (char*)malloc((*a_length+1)*sizeof(char));
  strncpy(name_copy, a_Name, *a_length);
  name_copy[*a_length] = '\0';

  error = xfWritePropertyDouble(*a_Id, name_copy, *a_Number, a_Properties,
                                *a_Compression);

  free(name_copy);
  return error;

 } /*xfWritePropertyDouble_f_*/
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfWritePropertyFloat_f_ */
/* PURPOSE   */
/* NOTES     */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFWRITEPROPERTYFLOAT_F,xfwritepropertyfloat_f,xfwritepropertyfloat_f_) XFWRITEPROPERTYFLOAT_F
#elif defined FLOWER
#define func_name(XFWRITEPROPERTYFLOAT_F,xfwritepropertyfloat_f,xfwritepropertyfloat_f_) xfwritepropertyfloat_f
#else
#define func_name(XFWRITEPROPERTYFLOAT_F,xfwritepropertyfloat_f,xfwritepropertyfloat_f_) xfwritepropertyfloat_f_
#endif
XMDF_API  int func_name(XFWRITEPROPERTYFLOAT_F,xfwritepropertyfloat_f,xfwritepropertyfloat_f_)
              (xid *a_Id, const char *a_Name, int *a_length,  
                            int *a_Number, float *a_Properties, 
                            int *a_Compression)
 {

  int	error;
  char* name_copy = NULL;

  name_copy = (char*)malloc((*a_length+1)*sizeof(char));
  strncpy(name_copy, a_Name, *a_length);
  name_copy[*a_length] = '\0';

  error = xfWritePropertyFloat(*a_Id, name_copy, *a_Number, a_Properties,
                               *a_Compression);

  free(name_copy);
  return error;

 } /*xfWritePropertyFloat_f_*/
 /*---------------------------------------------------------------------------- */
 /* FUNCTION  xfDoesPropertyWithNameExist */
 /* PURPOSE   Looks to see if an Property with a given name exists */
 /* NOTES     */
 /*---------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFDOESPROPERTYWITHNAMEEXIST_F,xfdoespropertywithnameexist_f,xfdoespropertywithnameexist_f_) XFDOESPROPERTYWITHNAMEEXIST_F
#elif defined FLOWER
#define func_name(XFDOESPROPERTYWITHNAMEEXIST_F,xfdoespropertywithnameexist_f,xfdoespropertywithnameexist_f_) xfdoespropertywithnameexist_f
#else
#define func_name(XFDOESPROPERTYWITHNAMEEXIST_F,xfdoespropertywithnameexist_f,xfdoespropertywithnameexist_f_) xfdoespropertywithnameexist_f_
#endif
XMDF_API  int func_name(XFDOESPROPERTYWITHNAMEEXIST_F,xfdoespropertywithnameexist_f,xfdoespropertywithnameexist_f_)
              (xid *a_Id, const char *a_Name, int *a_length,
                                   short *a_Exists)
 { 

  xmbool  Exists;
  int     error;
  char   *name_copy = NULL;

  Exists = (*a_Exists) ? XTRUE : XFALSE;

  name_copy = (char*)malloc((*a_length+1)*sizeof(char));
  strncpy(name_copy, a_Name, *a_length);
  name_copy[*a_length] = '\0';

  error = xfDoesPropertyWithNameExist(*a_Id, name_copy, &Exists);

  free(name_copy);
  return error;


 } /*xfDoesPropertyWithNameExist_f_*/
 /******************************************************************************
 * FUNCTION  xfGetPropertyNumber_f_
 * PURPOSE   Gets the number of items in a property array
 * NOTES     
 ******************2***********************************************************/
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFGETPROPERTYNUMBER_F,xfgetpropertynumber_f,xfgetpropertynumber_f_) XFGETPROPERTYNUMBER_F
#elif defined FLOWER
#define func_name(XFGETPROPERTYNUMBER_F,xfgetpropertynumber_f,xfgetpropertynumber_f_) xfgetpropertynumber_f
#else
#define func_name(XFGETPROPERTYNUMBER_F,xfgetpropertynumber_f,xfgetpropertynumber_f_) xfgetpropertynumber_f_
#endif
XMDF_API  int func_name(XFGETPROPERTYNUMBER_F,xfgetpropertynumber_f,xfgetpropertynumber_f_)
              (xid *a_Id, const char *a_Name, int *a_namelen,
                            int *a_Number)
{
  
  int	error;
  char* name_copy = NULL;

  name_copy = (char*)malloc((*a_namelen+1)*sizeof(char));
  strncpy(name_copy, a_Name, *a_namelen);
  name_copy[*a_namelen] = '\0';

  error = xfGetPropertyNumber(*a_Id, name_copy, a_Number);

  free(name_copy);
  return error;

 } /*xfGetPropertyNumber_f_*/
/******************************************************************************
 * FUNCTION  xfreadpropertyint_f_
 * PURPOSE   Reads the int dataset from the properties directory
 * NOTES     the variable properties must already be allocated to a_Number
 ******************2***********************************************************/
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFREADPROPERTYINT_F,xfreadpropertyint_f,xfreadpropertyint_f_) XFREADPROPERTYINT_F
#elif defined FLOWER
#define func_name(XFREADPROPERTYINT_F,xfreadpropertyint_f,xfreadpropertyint_f_) xfreadpropertyint_f
#else
#define func_name(XFREADPROPERTYINT_F,xfreadpropertyint_f,xfreadpropertyint_f_) xfreadpropertyint_f_
#endif
XMDF_API  int func_name(XFREADPROPERTYINT_F,xfreadpropertyint_f,xfreadpropertyint_f_)
              (xid *a_Id, const char *a_Name, int *a_namelen, 
                          int *a_Number, int *a_Properties)
{

  int	error;
  char* name_copy = NULL;

  name_copy = (char*)malloc((*a_namelen+1)*sizeof(char));
  strncpy(name_copy, a_Name, *a_namelen);
  name_copy[*a_namelen] = '\0';

  error = xfReadPropertyInt(*a_Id, name_copy, *a_Number, a_Properties);

  free(name_copy);
  return error;

 } /*xfReadPropertyInt_f_*/
/*-----------------------------------------------------------------------------*/
/* FUNCTION  xfGetPropertyStringLength_f_ */
/* PURPOSE   */
/* NOTES */
/*-----------------------------------------------------------------------------*/
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFGETPROPERTYSTRINGLENGTH_F,xfgetpropertystringlength_f,xfgetpropertystringlength_f_) XFGETPROPERTYSTRINGLENGTH_F
#elif defined FLOWER
#define func_name(XFGETPROPERTYSTRINGLENGTH_F,xfgetpropertystringlength_f,xfgetpropertystringlength_f_) xfgetpropertystringlength_f
#else
#define func_name(XFGETPROPERTYSTRINGLENGTH_F,xfgetpropertystringlength_f,xfgetpropertystringlength_f_) xfgetpropertystringlength_f_
#endif
XMDF_API  int func_name(XFGETPROPERTYSTRINGLENGTH_F,xfgetpropertystringlength_f,xfgetpropertystringlength_f_)
              (xid *a_Id, const char *a_Name, int *a_length,
                                 int *a_Number, int *a_MaxLength)
 {
  
  int	error;
  char* name_copy = NULL;

  name_copy = (char*)malloc((*a_length+1)*sizeof(char));
  strncpy(name_copy, a_Name, *a_length);
  name_copy[*a_length] = '\0';

  error = xfGetPropertyStringLength(*a_Id, name_copy, a_Number, a_MaxLength);

  free(name_copy);
  return error;

 } /*xfGetPropertyStringLength_f_*/
/*-----------------------------------------------------------------------------*/
/* FUNCTION  xfGetPropertyType_f_ */
/* PURPOSE   Gets the property type from a dataset */
/* NOTES */
/*-----------------------------------------------------------------------------*/
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFGETPROPERTYTYPE_F,xfgetpropertytype_f,xfgetpropertytype_f_) XFGETPROPERTYTYPE_F
#elif defined FLOWER
#define func_name(XFGETPROPERTYTYPE_F,xfgetpropertytype_f,xfgetpropertytype_f_) xfgetpropertytype_f
#else
#define func_name(XFGETPROPERTYTYPE_F,xfgetpropertytype_f,xfgetpropertytype_f_) xfgetpropertytype_f_
#endif
XMDF_API  int func_name(XFGETPROPERTYTYPE_F,xfgetpropertytype_f,xfgetpropertytype_f_)
              (xid *a_GroupId, const char *a_Name, int *a_length,
                         int *a_Type)
 {

  int	error;
  char* name_copy = NULL;

  name_copy = (char*)malloc((*a_length+1)*sizeof(char));
  strncpy(name_copy, a_Name, *a_length);
  name_copy[*a_length] = '\0';

  error = xfGetPropertyType(*a_GroupId, name_copy, a_Type);

  free(name_copy);
  return error;

 } /*xfGetPropertyType_f_*/
 /*-----------------------------------------------------------------------------*/
 /* FUNCTION  xfReadPropertyString_f_ */
 /* PURPOSE   Reads the string dataset from the attributes directory */
 /* NOTES     the variable attributes must already be allocated */
 /*           to a_Number*a_MaxLength */
 /*-----------------------------------------------------------------------------*/
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFREADPROPERTYSTRING_F,xfreadpropertystring_f,xfreadpropertystring_f_) XFREADPROPERTYSTRING_F
#elif defined FLOWER
#define func_name(XFREADPROPERTYSTRING_F,xfreadpropertystring_f,xfreadpropertystring_f_) xfreadpropertystring_f
#else
#define func_name(XFREADPROPERTYSTRING_F,xfreadpropertystring_f,xfreadpropertystring_f_) xfreadpropertystring_f_
#endif
XMDF_API  int func_name(XFREADPROPERTYSTRING_F,xfreadpropertystring_f,xfreadpropertystring_f_)
              (xid *a_Id, const char *a_Name, int *a_length, 
                            int *a_Number, int *a_MaxLength, char *a_Properties,
                            int *a_attrlength)
{

  int	error;
  char* name_copy = NULL;

  name_copy = (char*)malloc((*a_length+1)*sizeof(char));
  strncpy(name_copy, a_Name, *a_length);
  name_copy[*a_length] = '\0';

  error = xfReadPropertyString(*a_Id, name_copy, *a_Number, *a_MaxLength, a_Properties);

  free(name_copy);
  return error;

 } /*xfReadPropertyString_f_*/
 /*-----------------------------------------------------------------------------*/
 /* FUNCTION  xfAllocateReadPropertyString_f */
 /* PURPOSE   Reads the string dataset from the attributes directory */
 /* NOTES     This allocates the variable automatically */
 /*-----------------------------------------------------------------------------*/
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFALLOCATEREADPROPERTYSTRING_,xfAllocateReadPropertyString_,xfAllocateReadPropertyString_f) XFALLOCATEREADPROPERTYSTRING_
#elif defined FLOWER
#define func_name(XFALLOCATEREADPROPERTYSTRING_,xfAllocateReadPropertyString_,xfAllocateReadPropertyString_f) xfAllocateReadPropertyString_
#else
#define func_name(XFALLOCATEREADPROPERTYSTRING_,xfAllocateReadPropertyString_,xfAllocateReadPropertyString_f) xfAllocateReadPropertyString_f
#endif
XMDF_API  int func_name(XFALLOCATEREADPROPERTYSTRING_,xfAllocateReadPropertyString_,xfAllocateReadPropertyString_f)
              (xid *a_Id, const char *a_Name, int *a_length,
                                    int *a_Number, int *a_MaxLength,
                                    char **a_Properties, int *a_attrlength)
{

  int	error;
  char* name_copy = NULL;

  name_copy = (char*)malloc((*a_length+1)*sizeof(char));
  strncpy(name_copy, a_Name, *a_length);
  name_copy[*a_length] = '\0';

  error = xfAllocateReadPropertyString(*a_Id, name_copy, a_Number, a_MaxLength,
									    a_Properties);
  free(name_copy);
  return error;

 } /*xfAllocateReadPropertyString_f*/
 /*-----------------------------------------------------------------------------*/
 /* FUNCTION  xfReadPropertyFloat_f */
 /* FUNCTION  xfReadPropertyFloat_f_ */
 /* PURPOSE   Reads the float dataset from the properties directory */
 /* NOTES     the variable properties must already be allocated to a_Number */
 /*-----------------------------------------------------------------------------*/
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFREADPROPERTYFLOAT_F,xfreadpropertyfloat_f,xfreadpropertyfloat_f_) XFREADPROPERTYFLOAT_F
#elif defined FLOWER
#define func_name(XFREADPROPERTYFLOAT_F,xfreadpropertyfloat_f,xfreadpropertyfloat_f_) xfreadpropertyfloat_f
#else
#define func_name(XFREADPROPERTYFLOAT_F,xfreadpropertyfloat_f,xfreadpropertyfloat_f_) xfreadpropertyfloat_f_
#endif
XMDF_API  int func_name(XFREADPROPERTYFLOAT_F,xfreadpropertyfloat_f,xfreadpropertyfloat_f_)
              (xid *a_Id, const char *a_Name, int *a_length, 
                           int *a_Number, float *a_Properties)
{

  int	error;
  char* name_copy = NULL;

  name_copy = (char*)malloc((*a_length+1)*sizeof(char));
  strncpy(name_copy, a_Name, *a_length);
  name_copy[*a_length] = '\0';

  error = xfReadPropertyFloat(*a_Id, name_copy, *a_Number, a_Properties);

  free(name_copy);
  return error;

 } /*xfReadPropertyFloat_f_*/
 /******************************************************************************
 * FUNCTION  xfReadPropertyDouble_f_
 * PURPOSE   Reads the double dataset from the properties directory
 * NOTES     the variable properties must already be allocated to a_Number
 ******************2***********************************************************/
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFREADPROPERTYDOUBLE_F,xfreadpropertydouble_f,xfreadpropertydouble_f_) XFREADPROPERTYDOUBLE_F
#elif defined FLOWER
#define func_name(XFREADPROPERTYDOUBLE_F,xfreadpropertydouble_f,xfreadpropertydouble_f_) xfreadpropertydouble_f
#else
#define func_name(XFREADPROPERTYDOUBLE_F,xfreadpropertydouble_f,xfreadpropertydouble_f_) xfreadpropertydouble_f_
#endif
XMDF_API  int func_name(XFREADPROPERTYDOUBLE_F,xfreadpropertydouble_f,xfreadpropertydouble_f_)
              (xid *a_Id, const char *a_Name, int *a_namelen,
                             int *a_Number, double *a_Properties)
{
  int	error;
  char* name_copy = NULL;

  name_copy = (char*)malloc((*a_namelen+1)*sizeof(char));
  strncpy(name_copy, a_Name, *a_namelen);
  name_copy[*a_namelen] = '\0';

  error = xfReadPropertyDouble(*a_Id, name_copy, *a_Number, a_Properties);

  free(name_copy);
  return error;

 } /*xfReadPropertyDouble_f_*/
 /* --------------------------------------------------------------------------- */
 /* FUNCTION  xfCreateGenericGroup_f_ */
 /* PURPOSE   Create a group of generic type (anything may be stored below) */
 /* NOTES     These groups are just created for organization */
 /* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFCREATEGENERICGROUP_F,xfcreategenericgroup_f,xfcreategenericgroup_f_) XFCREATEGENERICGROUP_F
#elif defined FLOWER
#define func_name(XFCREATEGENERICGROUP_F,xfcreategenericgroup_f,xfcreategenericgroup_f_) xfcreategenericgroup_f
#else
#define func_name(XFCREATEGENERICGROUP_F,xfcreategenericgroup_f,xfcreategenericgroup_f_) xfcreategenericgroup_f_
#endif
XMDF_API  int func_name(XFCREATEGENERICGROUP_F,xfcreategenericgroup_f,xfcreategenericgroup_f_)
              (xid *a_FileId, const char *a_Path, int *a_length,
                             xid *a_GroupId)
 {
  int	error;
  char *path_copy = NULL;

  path_copy = (char*)malloc((*a_length+1)*sizeof(char));
  strncpy(path_copy, a_Path, *a_length);
  path_copy[*a_length] = '\0';

  error = xfCreateGenericGroup(*a_FileId, path_copy, a_GroupId);

  free(path_copy);
  return error;

 } /*xfCreateGenericGroup_f_*/
 /* --------------------------------------------------------------------------- */
 /* FUNCTION  xfCreateGroupForMesh_f_ */
 /* PURPOSE   Create an Xformat group to store a mesh */
 /* NOTES    */
 /* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFCREATEGROUPFORMESH_F,xfcreategroupformesh_f,xfcreategroupformesh_f_) XFCREATEGROUPFORMESH_F
#elif defined FLOWER
#define func_name(XFCREATEGROUPFORMESH_F,xfcreategroupformesh_f,xfcreategroupformesh_f_) xfcreategroupformesh_f
#else
#define func_name(XFCREATEGROUPFORMESH_F,xfcreategroupformesh_f,xfcreategroupformesh_f_) xfcreategroupformesh_f_
#endif
XMDF_API  int func_name(XFCREATEGROUPFORMESH_F,xfcreategroupformesh_f,xfcreategroupformesh_f_)
              (xid *a_FileId, const char *a_Path, int *a_length,
                            xid *a_GroupId)
 {

  int error;
  char *path_copy = NULL;

  path_copy = (char*)malloc((*a_length+1)*sizeof(char));
  strncpy(path_copy, a_Path, *a_length);
  path_copy[*a_length] = '\0';

  error = xfCreateGroupForMesh(*a_FileId, path_copy, a_GroupId);

  free(path_copy);
  return error;

 } /*xfCreateGroupForMesh_f_*/
 /* --------------------------------------------------------------------------- */
 /* FUNCTION  xfCreateGroupForGrid_f_ */
 /* PURPOSE   Create an Xformat group to store a mesh */
 /* NOTES    */
 /* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFCREATEGROUPFORGRID_F,xfcreategroupforgrid_f,xfcreategroupforgrid_f_) XFCREATEGROUPFORGRID_F
#elif defined FLOWER
#define func_name(XFCREATEGROUPFORGRID_F,xfcreategroupforgrid_f,xfcreategroupforgrid_f_) xfcreategroupforgrid_f
#else
#define func_name(XFCREATEGROUPFORGRID_F,xfcreategroupforgrid_f,xfcreategroupforgrid_f_) xfcreategroupforgrid_f_
#endif
XMDF_API  int func_name(XFCREATEGROUPFORGRID_F,xfcreategroupforgrid_f,xfcreategroupforgrid_f_)
              (xid *a_FileId, const char *a_Path, int *a_length, 
                            xid *a_GroupId)
 {

  int error;
  char *path_copy = NULL;

  path_copy = (char*)malloc((*a_length+1)*sizeof(char));
  strncpy(path_copy, a_Path, *a_length);
  path_copy[*a_length] = '\0';

  error = xfCreateGroupForGrid(*a_FileId, path_copy, a_GroupId);

  free(path_copy);
  return error;

 } /*xfCreateGroupForGrid_f_*/
 /******************************************************************************
 * FUNCTION  xfCreateGroupForXsec_f_
 * PURPOSE   Create an Xformat group to store a xsec
 * NOTES    
 ******************2***********************************************************/
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFCREATEGROUPFORXSEC_F,xfcreategroupforxsec_f,xfcreategroupforxsec_f_) XFCREATEGROUPFORXSEC_F
#elif defined FLOWER
#define func_name(XFCREATEGROUPFORXSEC_F,xfcreategroupforxsec_f,xfcreategroupforxsec_f_) xfcreategroupforxsec_f
#else
#define func_name(XFCREATEGROUPFORXSEC_F,xfcreategroupforxsec_f,xfcreategroupforxsec_f_) xfcreategroupforxsec_f_
#endif
XMDF_API  int func_name(XFCREATEGROUPFORXSEC_F,xfcreategroupforxsec_f,xfcreategroupforxsec_f_)
              (xid *a_FileId, const char *a_Path, int *a_pathlen,
                             xid *a_GroupId)
{
  int error;
  char *path_copy = NULL;

  path_copy = (char*)malloc((*a_pathlen+1)*sizeof(char));
  strncpy(path_copy, a_Path, *a_pathlen);
  path_copy[*a_pathlen] = '\0';

  error = xfCreateGroupForXsec(*a_FileId, path_copy, a_GroupId);

  free(path_copy);
  return error;

 } /*xfCreateGroupForXsec_f_*/
 /* --------------------------------------------------------------------------- */
 /* FUNCTION  xfOpenGroup_f_ */
 /* PURPOSE   Open an Xformat group */
 /* NOTES    */
 /* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFOPENGROUP_F,xfopengroup_f,xfopengroup_f_) XFOPENGROUP_F
#elif defined FLOWER
#define func_name(XFOPENGROUP_F,xfopengroup_f,xfopengroup_f_) xfopengroup_f
#else
#define func_name(XFOPENGROUP_F,xfopengroup_f,xfopengroup_f_) xfopengroup_f_
#endif
XMDF_API  int func_name(XFOPENGROUP_F,xfopengroup_f,xfopengroup_f_)
              (xid *a_ParentId, const char *a_Path, int *a_length, 
                    xid *a_GroupId)
 {
    
  int   error;
  char *path_copy = NULL;

  path_copy = (char*)malloc((*a_length+1)*sizeof(char));
  strncpy(path_copy, a_Path, *a_length);
  path_copy[*a_length] = '\0';

  error = xfOpenGroup(*a_ParentId, path_copy, a_GroupId);
 
  free(path_copy);
  return error;

 } /*xfOpenGroup_f_*/
 /* --------------------------------------------------------------------------- */
 /* FUNCTION  xfCloseGroup_f_ */
 /* PURPOSE   Close an Xformat group */
 /* NOTES     */
 /* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFCLOSEGROUP_F,xfclosegroup_f,xfclosegroup_f_) XFCLOSEGROUP_F
#elif defined FLOWER
#define func_name(XFCLOSEGROUP_F,xfclosegroup_f,xfclosegroup_f_) xfclosegroup_f
#else
#define func_name(XFCLOSEGROUP_F,xfclosegroup_f,xfclosegroup_f_) xfclosegroup_f_
#endif
XMDF_API  int func_name(XFCLOSEGROUP_F,xfclosegroup_f,xfclosegroup_f_)
              (xid *a_GroupId)
 {

  int error;

  error = xfCloseGroup(*a_GroupId);

  return error;

 } /*xfCloseGroup_f_*/
 /* --------------------------------------------------------------------------- */
 /* FUNCTION  xfGetGroupPathsSizeForMeshes_f_ */
 /* PURPOSE   Find the number and maximum size for paths to datasets in an */
 /*           Xmdf file */
 /* NOTES    */
 /* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFGETGROUPPATHSSIZEFORMESHES_F,xfgetgrouppathssizeformeshes_f,xfgetgrouppathssizeformeshes_f_) XFGETGROUPPATHSSIZEFORMESHES_F
#elif defined FLOWER
#define func_name(XFGETGROUPPATHSSIZEFORMESHES_F,xfgetgrouppathssizeformeshes_f,xfgetgrouppathssizeformeshes_f_) xfgetgrouppathssizeformeshes_f
#else
#define func_name(XFGETGROUPPATHSSIZEFORMESHES_F,xfgetgrouppathssizeformeshes_f,xfgetgrouppathssizeformeshes_f_) xfgetgrouppathssizeformeshes_f_
#endif
XMDF_API  int func_name(XFGETGROUPPATHSSIZEFORMESHES_F,xfgetgrouppathssizeformeshes_f,xfgetgrouppathssizeformeshes_f_)
              (xid *a_FileId, int *a_Number, int *a_Maxsize)
 {

  int error;

  error = xfGetGroupPathsSizeForMeshes(*a_FileId, a_Number, a_Maxsize);

  return error;

 } /*xfGetGroupPathsSizeForMeshes_f_*/
 /* --------------------------------------------------------------------------- */
 /* FUNCTION  xfGetGroupPathsSizeForMeshes_f_ */
 /* PURPOSE   Find the number and maximum size for paths to datasets in an */
 /*           Xmdf file */
 /* NOTES    */
 /* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFGETGROUPPATHSFORMESHES_F,xfgetgrouppathsformeshes_f,xfgetgrouppathsformeshes_f_) XFGETGROUPPATHSFORMESHES_F
#elif defined FLOWER
#define func_name(XFGETGROUPPATHSFORMESHES_F,xfgetgrouppathsformeshes_f,xfgetgrouppathsformeshes_f_) xfgetgrouppathsformeshes_f
#else
#define func_name(XFGETGROUPPATHSFORMESHES_F,xfgetgrouppathsformeshes_f,xfgetgrouppathsformeshes_f_) xfgetgrouppathsformeshes_f_
#endif
XMDF_API  int func_name(XFGETGROUPPATHSFORMESHES_F,xfgetgrouppathsformeshes_f,xfgetgrouppathsformeshes_f_)
              (xid *a_FileId, int *a_Num, int *a_Maxsize, 
                                 char *a_Paths, int *a_length)
{

  int   error;
  error = xfGetGroupPathsForMeshes(*a_FileId, *a_Num, *a_Maxsize, a_Paths);
  return error;

 } /*xfGetGroupPathsForMeshes_f_*/
 /* --------------------------------------------------------------------------- */
 /* FUNCTION  xfGetGroupPathsSizeForGrids_f_ */
 /* PURPOSE   Find the number and maximum size for paths to Grids in an */
 /*           Xmdf file */
 /* NOTES   */ 
 /* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFGETGROUPPATHSSIZEFORGRIDS_F,xfgetgrouppathssizeforgrids_f,xfgetgrouppathssizeforgrids_f_) XFGETGROUPPATHSSIZEFORGRIDS_F
#elif defined FLOWER
#define func_name(XFGETGROUPPATHSSIZEFORGRIDS_F,xfgetgrouppathssizeforgrids_f,xfgetgrouppathssizeforgrids_f_) xfgetgrouppathssizeforgrids_f
#else
#define func_name(XFGETGROUPPATHSSIZEFORGRIDS_F,xfgetgrouppathssizeforgrids_f,xfgetgrouppathssizeforgrids_f_) xfgetgrouppathssizeforgrids_f_
#endif
XMDF_API  int func_name(XFGETGROUPPATHSSIZEFORGRIDS_F,xfgetgrouppathssizeforgrids_f,xfgetgrouppathssizeforgrids_f_)
              (xid *a_FileId, int *a_Number, int *a_Maxsize)
 {

  int	error;

  error = xfGetGroupPathsSizeForGrids(*a_FileId, a_Number, a_Maxsize);

  return error;

 } /*xfGetGroupPathsSizeForGrids_f_*/
 /* --------------------------------------------------------------------------- 
 * FUNCTION  xfGetGroupPathsSizeForGrids_f_
 * PURPOSE   Find the number and maximum size for paths to datasets in an 
 *           Xmdf file
 * NOTES    
 * --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFGETGROUPPATHSFORGRIDS_F,xfgetgrouppathsforgrids_f,xfgetgrouppathsforgrids_f_) XFGETGROUPPATHSFORGRIDS_F
#elif defined FLOWER
#define func_name(XFGETGROUPPATHSFORGRIDS_F,xfgetgrouppathsforgrids_f,xfgetgrouppathsforgrids_f_) xfgetgrouppathsforgrids_f
#else
#define func_name(XFGETGROUPPATHSFORGRIDS_F,xfgetgrouppathsforgrids_f,xfgetgrouppathsforgrids_f_) xfgetgrouppathsforgrids_f_
#endif
XMDF_API  int func_name(XFGETGROUPPATHSFORGRIDS_F,xfgetgrouppathsforgrids_f,xfgetgrouppathsforgrids_f_)
              (xid *a_FileId, int *a_Num, int *a_Maxsize, 
                                char *a_Paths, int *a_length)
{

  int	error;
  error = xfGetGroupPathsForGrids(*a_FileId, *a_Num, *a_Maxsize, a_Paths);
  return error;

 } /*xfGetGroupPathsForGrids_f_*/
 /******************************************************************************
 * FUNCTION  xfGetGroupPathsSizeForXsecs_f_
 * PURPOSE   Find the number and maximum size for paths to Xsecs in an 
 *           Xmdf file
 * NOTES    
 ******************2***********************************************************/
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFGETGROUPPATHSSIZEFORXSECS_F,xfgetgrouppathssizeforxsecs_f,xfgetgrouppathssizeforxsecs_f_) XFGETGROUPPATHSSIZEFORXSECS_F
#elif defined FLOWER
#define func_name(XFGETGROUPPATHSSIZEFORXSECS_F,xfgetgrouppathssizeforxsecs_f,xfgetgrouppathssizeforxsecs_f_) xfgetgrouppathssizeforxsecs_f
#else
#define func_name(XFGETGROUPPATHSSIZEFORXSECS_F,xfgetgrouppathssizeforxsecs_f,xfgetgrouppathssizeforxsecs_f_) xfgetgrouppathssizeforxsecs_f_
#endif
XMDF_API  int func_name(XFGETGROUPPATHSSIZEFORXSECS_F,xfgetgrouppathssizeforxsecs_f,xfgetgrouppathssizeforxsecs_f_)
              (xid *a_FileId, int *a_Num, int *a_Maxsize)
{
  int error;

  error = xfGetGroupPathsSizeForXsecs(*a_FileId, a_Num, a_Maxsize);

  return error;

 } /*xfGetGroupPathsSizeForXsecs_f_*/
 /******************************************************************************
 * FUNCTION  xfGetGroupPathsForXsecs_f_
 * PURPOSE   Find the number and maximum size for paths to Xsecs in an 
 *           Xmdf file
 * NOTES    
 ******************2***********************************************************/
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFGETGROUPPATHSFORXSECS_F,xfgetgrouppathsforxsecs_f,xfgetgrouppathsforxsecs_f_) XFGETGROUPPATHSFORXSECS_F
#elif defined FLOWER
#define func_name(XFGETGROUPPATHSFORXSECS_F,xfgetgrouppathsforxsecs_f,xfgetgrouppathsforxsecs_f_) xfgetgrouppathsforxsecs_f
#else
#define func_name(XFGETGROUPPATHSFORXSECS_F,xfgetgrouppathsforxsecs_f,xfgetgrouppathsforxsecs_f_) xfgetgrouppathsforxsecs_f_
#endif
XMDF_API  int func_name(XFGETGROUPPATHSFORXSECS_F,xfgetgrouppathsforxsecs_f,xfgetgrouppathsforxsecs_f_)
              (xid *a_FileId, int *a_Num, int *a_Maxsize,
                                char *a_Paths, int *a_pathlen)
{
  int	error;
  error = xfGetGroupPathsForXsecs(*a_FileId, *a_Num, *a_Maxsize, a_Paths);
  return error;

 } /*xfGetGroupPathsForXsecs_f_*/
 /* --------------------------------------------------------------------------- 
 * FUNCTION  xfGetPathSizeForGeomPaths_f_
 * PURPOSE   Find the number and maximum size for paths to datasets in an 
 *           Xmdf file
 * NOTES    
 * --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFGETPATHSIZEFORGEOMPATHS_F,xfgetpathsizeforgeompaths_f,xfgetpathsizeforgeompaths_f_) XFGETPATHSIZEFORGEOMPATHS_F
#elif defined FLOWER
#define func_name(XFGETPATHSIZEFORGEOMPATHS_F,xfgetpathsizeforgeompaths_f,xfgetpathsizeforgeompaths_f_) xfgetpathsizeforgeompaths_f
#else
#define func_name(XFGETPATHSIZEFORGEOMPATHS_F,xfgetpathsizeforgeompaths_f,xfgetpathsizeforgeompaths_f_) xfgetpathsizeforgeompaths_f_
#endif
XMDF_API  int func_name(XFGETPATHSIZEFORGEOMPATHS_F,xfgetpathsizeforgeompaths_f,xfgetpathsizeforgeompaths_f_)
              (xid *a_FileId, int *a_Number, int *a_Maxsize)
{

  int error;

  error = xfGetGroupPathsSizeForGeomPaths(*a_FileId, a_Number, a_Maxsize);

  return error;

 } /*xfGetPathSizeForGeomPaths_f_*/
 /* --------------------------------------------------------------------------- 
 * FUNCTION  xfGetGroupPathsForGeomPaths_f_
 * PURPOSE   Find the number and maximum size for paths to datasets in an 
 *           Xmdf file
 * NOTES    
 * --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFGETGROUPPATHSFORGEOMPATHS_F,xfgetgrouppathsforgeompaths_f,xfgetgrouppathsforgeompaths_f_) XFGETGROUPPATHSFORGEOMPATHS_F
#elif defined FLOWER
#define func_name(XFGETGROUPPATHSFORGEOMPATHS_F,xfgetgrouppathsforgeompaths_f,xfgetgrouppathsforgeompaths_f_) xfgetgrouppathsforgeompaths_f
#else
#define func_name(XFGETGROUPPATHSFORGEOMPATHS_F,xfgetgrouppathsforgeompaths_f,xfgetgrouppathsforgeompaths_f_) xfgetgrouppathsforgeompaths_f_
#endif
XMDF_API  int func_name(XFGETGROUPPATHSFORGEOMPATHS_F,xfgetgrouppathsforgeompaths_f,xfgetgrouppathsforgeompaths_f_)
              (xid *a_FileId, int *a_Num, int *a_Maxsize, 
                                    char *a_Paths, int *a_length)
 {

  int   error;
  error = xfGetGroupPathsForGeomPaths(*a_FileId, *a_Num, *a_Maxsize, a_Paths);
  return error;

 } /*xfGetGroupPathsForGeomPaths_f_*/
 /* ---------------------------------------------------------------------------*/
/* FUNCTION  xfSetNumberOfElements_f_ */
/* PURPOSE   Set the number of elements for the mesh*/
/* NOTES     This number is written as an attribute to a mesh*/
/* ---------------------------------------------------------------------------*/
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFSETNUMBEROFELEMENTS_F,xfsetnumberofelements_f,xfsetnumberofelements_f_) XFSETNUMBEROFELEMENTS_F
#elif defined FLOWER
#define func_name(XFSETNUMBEROFELEMENTS_F,xfsetnumberofelements_f,xfsetnumberofelements_f_) xfsetnumberofelements_f
#else
#define func_name(XFSETNUMBEROFELEMENTS_F,xfsetnumberofelements_f,xfsetnumberofelements_f_) xfsetnumberofelements_f_
#endif
XMDF_API  int func_name(XFSETNUMBEROFELEMENTS_F,xfsetnumberofelements_f,xfsetnumberofelements_f_)
              (xid *a_Id, int *a_nElems)
 {

  int error;

  error = xfSetNumberOfElements(*a_Id, *a_nElems);

  return error;

 } /*xfSetNumberOfElements_f_*/
 /* --------------------------------------------------------------------------- */
 /* FUNCTION  xfOpenPropertyGroup_f_*/
 /* PURPOSE   Open a property Group */
 /* NOTES      */
 /* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFOPENPROPERTYGROUP_F,xfopenpropertygroup_f,xfopenpropertygroup_f_) XFOPENPROPERTYGROUP_F
#elif defined FLOWER
#define func_name(XFOPENPROPERTYGROUP_F,xfopenpropertygroup_f,xfopenpropertygroup_f_) xfopenpropertygroup_f
#else
#define func_name(XFOPENPROPERTYGROUP_F,xfopenpropertygroup_f,xfopenpropertygroup_f_) xfopenpropertygroup_f_
#endif
XMDF_API  int func_name(XFOPENPROPERTYGROUP_F,xfopenpropertygroup_f,xfopenpropertygroup_f_)
              (xid *a_ParentId, xid *a_GroupId)
 {

  int error;

  error = xfOpenPropertyGroup(*a_ParentId, a_GroupId);

  return error;

 } /*xfOpenPropertyGroup_f_*/
 /* ---------------------------------------------------------------------------*/
/* FUNCTION  xfGetGroupAbsolutePathSize_f_ */
/* PURPOSE   returns the size of the absolute path of the group with a given*/
/*           ID.*/
/* NOTES*/
/* ---------------------------------------------------------------------------*/
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFGETGROUPABSOLUTEPATHSIZE_F,xfgetgroupabsolutepathsize_f,xfgetgroupabsolutepathsize_f_) XFGETGROUPABSOLUTEPATHSIZE_F
#elif defined FLOWER
#define func_name(XFGETGROUPABSOLUTEPATHSIZE_F,xfgetgroupabsolutepathsize_f,xfgetgroupabsolutepathsize_f_) xfgetgroupabsolutepathsize_f
#else
#define func_name(XFGETGROUPABSOLUTEPATHSIZE_F,xfgetgroupabsolutepathsize_f,xfgetgroupabsolutepathsize_f_) xfgetgroupabsolutepathsize_f_
#endif
XMDF_API  int func_name(XFGETGROUPABSOLUTEPATHSIZE_F,xfgetgroupabsolutepathsize_f,xfgetgroupabsolutepathsize_f_)
              (xid *a_GroupId, int *a_PathLength)
{

  int error;

  error = xfGetGroupAbsolutePathSize(*a_GroupId, a_PathLength);

  return error;

 } /*xfGetGroupAbsolutePathSize_f_*/
 /* ---------------------------------------------------------------------------*/
/* FUNCTION  xfGetGroupAbsolutePath_f_ */
/* PURPOSE   returns the absolute path of the group with a given ID.*/
/* NOTES*/
/* ---------------------------------------------------------------------------*/
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFGETGROUPABSOLUTEPATH_F,xfgetgroupabsolutepath_f,xfgetgroupabsolutepath_f_) XFGETGROUPABSOLUTEPATH_F
#elif defined FLOWER
#define func_name(XFGETGROUPABSOLUTEPATH_F,xfgetgroupabsolutepath_f,xfgetgroupabsolutepath_f_) xfgetgroupabsolutepath_f
#else
#define func_name(XFGETGROUPABSOLUTEPATH_F,xfgetgroupabsolutepath_f,xfgetgroupabsolutepath_f_) xfgetgroupabsolutepath_f_
#endif
XMDF_API  int func_name(XFGETGROUPABSOLUTEPATH_F,xfgetgroupabsolutepath_f,xfgetgroupabsolutepath_f_)
              (xid *a_GroupId, int *a_PathLength,
                               char *a_Path, int *a_length)
{
  int	error;
  error = xfGetGroupAbsolutePath(*a_GroupId, *a_PathLength, a_Path);
  return error;

 } /*xfGetGroupAbsolutePath_f_*/
 /* --------------------------------------------------------------------------- */
/* FUNCTION  xfSetAllElemsSameType_f_ */
/* PURPOSE   Write the element type dataset as a single type which is applied */
/*           to all elements */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFSETALLELEMSSAMETYPE_F,xfsetallelemssametype_f,xfsetallelemssametype_f_) XFSETALLELEMSSAMETYPE_F
#elif defined FLOWER
#define func_name(XFSETALLELEMSSAMETYPE_F,xfsetallelemssametype_f,xfsetallelemssametype_f_) xfsetallelemssametype_f
#else
#define func_name(XFSETALLELEMSSAMETYPE_F,xfsetallelemssametype_f,xfsetallelemssametype_f_) xfsetallelemssametype_f_
#endif
XMDF_API  int func_name(XFSETALLELEMSSAMETYPE_F,xfsetallelemssametype_f,xfsetallelemssametype_f_)
              (xid *a_Id, int *a_Type)
{

  int error;

  error = xfSetAllElemsSameType(*a_Id, *a_Type);

  return error;

 } /*xfSetAllElemsSameType_f_*/ 
 /* --------------------------------------------------------------------------- */
/* FUNCTION  xfWriteElemTypes_f_ */
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFWRITEELEMTYPES_F,xfwriteelemtypes_f,xfwriteelemtypes_f_) XFWRITEELEMTYPES_F
#elif defined FLOWER
#define func_name(XFWRITEELEMTYPES_F,xfwriteelemtypes_f,xfwriteelemtypes_f_) xfwriteelemtypes_f
#else
#define func_name(XFWRITEELEMTYPES_F,xfwriteelemtypes_f,xfwriteelemtypes_f_) xfwriteelemtypes_f_
#endif
XMDF_API  int func_name(XFWRITEELEMTYPES_F,xfwriteelemtypes_f,xfwriteelemtypes_f_)
              (xid *a_Id, int *a_nElems, const int *a_Type,
                         int *a_Compression)
{

  int error;

  error = xfWriteElemTypes(*a_Id, *a_nElems, a_Type, *a_Compression);

  return error;

 } /*xfWriteElemTypes_f_*/
 /* --------------------------------------------------------------------------- */
/* FUNCTION  xfWriteElemNodeIds_f_ */
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFWRITEELEMNODEIDS_F,xfwriteelemnodeids_f,xfwriteelemnodeids_f_) XFWRITEELEMNODEIDS_F
#elif defined FLOWER
#define func_name(XFWRITEELEMNODEIDS_F,xfwriteelemnodeids_f,xfwriteelemnodeids_f_) xfwriteelemnodeids_f
#else
#define func_name(XFWRITEELEMNODEIDS_F,xfwriteelemnodeids_f,xfwriteelemnodeids_f_) xfwriteelemnodeids_f_
#endif
XMDF_API  int func_name(XFWRITEELEMNODEIDS_F,xfwriteelemnodeids_f,xfwriteelemnodeids_f_)
              (xid *a_Id, int *a_nElems, int *a_nMaxNodes,
                           int *a_Ids, int *a_Compression)
{

  int error;

  error = xfWriteElemNodeIds(*a_Id, *a_nElems, *a_nMaxNodes, a_Ids, *a_Compression);

  return error;

 } /*xfWriteElemNodeIds_f_*/
 /* --------------------------------------------------------------------------- */
/* FUNCTION  xfsetnumberofnodes_f_ */
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFSETNUMBEROFNODES_F,xfsetnumberofnodes_f,xfsetnumberofnodes_f_) XFSETNUMBEROFNODES_F
#elif defined FLOWER
#define func_name(XFSETNUMBEROFNODES_F,xfsetnumberofnodes_f,xfsetnumberofnodes_f_) xfsetnumberofnodes_f
#else
#define func_name(XFSETNUMBEROFNODES_F,xfsetnumberofnodes_f,xfsetnumberofnodes_f_) xfsetnumberofnodes_f_
#endif
XMDF_API  int func_name(XFSETNUMBEROFNODES_F,xfsetnumberofnodes_f,xfsetnumberofnodes_f_)
              (xid *a_Id, int *a_nNodes)
{

  int error;

  error = xfSetNumberOfNodes(*a_Id, *a_nNodes);

  return error;

 } /*xfSetNumberOfNodes_f_*/
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfWriteXNodeLocations_f_ */
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFWRITEXNODELOCATIONS_F,xfwritexnodelocations_f,xfwritexnodelocations_f_) XFWRITEXNODELOCATIONS_F
#elif defined FLOWER
#define func_name(XFWRITEXNODELOCATIONS_F,xfwritexnodelocations_f,xfwritexnodelocations_f_) xfwritexnodelocations_f
#else
#define func_name(XFWRITEXNODELOCATIONS_F,xfwritexnodelocations_f,xfwritexnodelocations_f_) xfwritexnodelocations_f_
#endif
XMDF_API  int func_name(XFWRITEXNODELOCATIONS_F,xfwritexnodelocations_f,xfwritexnodelocations_f_)
              (xid *a_Id, int *a_nNodes, double *a_Locs,
                              int *a_Compression)
{

  int error;

  error = xfWriteXNodeLocations(*a_Id, *a_nNodes, a_Locs, *a_Compression);

  return error;

 } /*xfWriteXNodeLocations_f_*/
 /* --------------------------------------------------------------------------- */
/* FUNCTION  xfWriteYNodeLocations_f_ */
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFWRITEYNODELOCATIONS_F,xfwriteynodelocations_f,xfwriteynodelocations_f_) XFWRITEYNODELOCATIONS_F
#elif defined FLOWER
#define func_name(XFWRITEYNODELOCATIONS_F,xfwriteynodelocations_f,xfwriteynodelocations_f_) xfwriteynodelocations_f
#else
#define func_name(XFWRITEYNODELOCATIONS_F,xfwriteynodelocations_f,xfwriteynodelocations_f_) xfwriteynodelocations_f_
#endif
XMDF_API  int func_name(XFWRITEYNODELOCATIONS_F,xfwriteynodelocations_f,xfwriteynodelocations_f_)
              (xid *a_Id, int *a_nNodes, double *a_Locs)
{

  int error;

  error = xfWriteYNodeLocations(*a_Id, *a_nNodes, a_Locs);

  return error;

 } /*xfWriteYNodeLocations_f_*/
 /* --------------------------------------------------------------------------- */
/* FUNCTION  xfWriteZNodeLocations_f_ */
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFWRITEZNODELOCATIONS_F,xfwriteznodelocations_f,xfwriteznodelocations_f_) XFWRITEZNODELOCATIONS_F
#elif defined FLOWER
#define func_name(XFWRITEZNODELOCATIONS_F,xfwriteznodelocations_f,xfwriteznodelocations_f_) xfwriteznodelocations_f
#else
#define func_name(XFWRITEZNODELOCATIONS_F,xfwriteznodelocations_f,xfwriteznodelocations_f_) xfwriteznodelocations_f_
#endif
XMDF_API  int func_name(XFWRITEZNODELOCATIONS_F,xfwriteznodelocations_f,xfwriteznodelocations_f_)
              (xid *a_Id, int *a_nNodes, double *a_Locs)
{

  int error;

  error = xfWriteZNodeLocations(*a_Id, *a_nNodes, a_Locs);

  return error;

 } /*xfWriteZNodeLocations_f_*/
 /* --------------------------------------------------------------------------- */
/* FUNCTION  xfGetNumberOfElements_f_ */
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFGETNUMBEROFELEMENTS_F,xfgetnumberofelements_f,xfgetnumberofelements_f_) XFGETNUMBEROFELEMENTS_F
#elif defined FLOWER
#define func_name(XFGETNUMBEROFELEMENTS_F,xfgetnumberofelements_f,xfgetnumberofelements_f_) xfgetnumberofelements_f
#else
#define func_name(XFGETNUMBEROFELEMENTS_F,xfgetnumberofelements_f,xfgetnumberofelements_f_) xfgetnumberofelements_f_
#endif
XMDF_API  int func_name(XFGETNUMBEROFELEMENTS_F,xfgetnumberofelements_f,xfgetnumberofelements_f_)
              (xid *a_Id, int *a_nElems)
{

  int error;

  error = xfGetNumberOfElements(*a_Id, a_nElems);

  return error;

 } /*xfGetNumberOfElements_f_*/
 /* --------------------------------------------------------------------------- */
/* FUNCTION  xfAreAllElemsSameType_f_ */
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFAREALLELEMSSAMETYPE_F,xfareallelemssametype_f,xfareallelemssametype_f_) XFAREALLELEMSSAMETYPE_F
#elif defined FLOWER
#define func_name(XFAREALLELEMSSAMETYPE_F,xfareallelemssametype_f,xfareallelemssametype_f_) xfareallelemssametype_f
#else
#define func_name(XFAREALLELEMSSAMETYPE_F,xfareallelemssametype_f,xfareallelemssametype_f_) xfareallelemssametype_f_
#endif
XMDF_API  int func_name(XFAREALLELEMSSAMETYPE_F,xfareallelemssametype_f,xfareallelemssametype_f_)
              (xid *a_Id, short *a_Same)
{

  xmbool  Same;
  int     error;

  Same = (*a_Same) ? XTRUE : XFALSE;
  error = xfAreAllElemsSameType(*a_Id, &Same);
  *a_Same = Same == XTRUE;

  return error;

 } /*xfAreAllElemsSameType_f_*/
 /* --------------------------------------------------------------------------- */
/* FUNCTION  xfReadElemTypesSingleValue_f_ */
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFREADELEMTYPESSINGLEVALUE_F,xfreadelemtypessinglevalue_f,xfreadelemtypessinglevalue_f_) XFREADELEMTYPESSINGLEVALUE_F
#elif defined FLOWER
#define func_name(XFREADELEMTYPESSINGLEVALUE_F,xfreadelemtypessinglevalue_f,xfreadelemtypessinglevalue_f_) xfreadelemtypessinglevalue_f
#else
#define func_name(XFREADELEMTYPESSINGLEVALUE_F,xfreadelemtypessinglevalue_f,xfreadelemtypessinglevalue_f_) xfreadelemtypessinglevalue_f_
#endif
XMDF_API  int func_name(XFREADELEMTYPESSINGLEVALUE_F,xfreadelemtypessinglevalue_f,xfreadelemtypessinglevalue_f_)
              (xid *a_Id, int *a_Type)
{

  int error;

  error = xfReadElemTypesSingleValue(*a_Id, a_Type);

  return error;

 } /*xfReadElemTypesSingleValue_f_*/
 /* --------------------------------------------------------------------------- */
/* FUNCTION  xfReadElemTypes_f_ */
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFREADELEMTYPES_F,xfreadelemtypes_f,xfreadelemtypes_f_) XFREADELEMTYPES_F
#elif defined FLOWER
#define func_name(XFREADELEMTYPES_F,xfreadelemtypes_f,xfreadelemtypes_f_) xfreadelemtypes_f
#else
#define func_name(XFREADELEMTYPES_F,xfreadelemtypes_f,xfreadelemtypes_f_) xfreadelemtypes_f_
#endif
XMDF_API  int func_name(XFREADELEMTYPES_F,xfreadelemtypes_f,xfreadelemtypes_f_)
              (xid *a_Id, int *a_nElems, int *a_Type)
{

  int error;

  error = xfReadElemTypes(*a_Id, *a_nElems, a_Type);

  return error;

 } /*xfReadElemTypes_f_*/
 /* --------------------------------------------------------------------------- */
/* FUNCTION  xfGetMaxNodesInElem_f_ */
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFGETMAXNODESINELEM_F,xfgetmaxnodesinelem_f,xfgetmaxnodesinelem_f_) XFGETMAXNODESINELEM_F
#elif defined FLOWER
#define func_name(XFGETMAXNODESINELEM_F,xfgetmaxnodesinelem_f,xfgetmaxnodesinelem_f_) xfgetmaxnodesinelem_f
#else
#define func_name(XFGETMAXNODESINELEM_F,xfgetmaxnodesinelem_f,xfgetmaxnodesinelem_f_) xfgetmaxnodesinelem_f_
#endif
XMDF_API  int func_name(XFGETMAXNODESINELEM_F,xfgetmaxnodesinelem_f,xfgetmaxnodesinelem_f_)
              (xid *a_Id, int *a_nMaxNodes)
{

  int error;

  error = xfGetMaxNodesInElem(*a_Id, a_nMaxNodes);

  return error;

 } /*xfGetMaxNodesInElem_f_*/
 /* --------------------------------------------------------------------------- */
/* FUNCTION  xfReadElemNodeIds_f_ */
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFREADELEMNODEIDS_F,xfreadelemnodeids_f,xfreadelemnodeids_f_) XFREADELEMNODEIDS_F
#elif defined FLOWER
#define func_name(XFREADELEMNODEIDS_F,xfreadelemnodeids_f,xfreadelemnodeids_f_) xfreadelemnodeids_f
#else
#define func_name(XFREADELEMNODEIDS_F,xfreadelemnodeids_f,xfreadelemnodeids_f_) xfreadelemnodeids_f_
#endif
XMDF_API  int func_name(XFREADELEMNODEIDS_F,xfreadelemnodeids_f,xfreadelemnodeids_f_)
              (xid *a_Id, int *a_nElems, int *a_nMaxNodes,
                          int *a_Ids)
{

  int error;

  error = xfReadElemNodeIds(*a_Id, *a_nElems, *a_nMaxNodes, a_Ids);

  return error;

 } /*xfReadElemNodeIds_f_*/
 /* --------------------------------------------------------------------------- */
/* FUNCTION  xfGetNumberOfNodes_f_ */
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFGETNUMBEROFNODES_F,xfgetnumberofnodes_f,xfgetnumberofnodes_f_) XFGETNUMBEROFNODES_F
#elif defined FLOWER
#define func_name(XFGETNUMBEROFNODES_F,xfgetnumberofnodes_f,xfgetnumberofnodes_f_) xfgetnumberofnodes_f
#else
#define func_name(XFGETNUMBEROFNODES_F,xfgetnumberofnodes_f,xfgetnumberofnodes_f_) xfgetnumberofnodes_f_
#endif
XMDF_API  int func_name(XFGETNUMBEROFNODES_F,xfgetnumberofnodes_f,xfgetnumberofnodes_f_)
              (xid *a_Id, int *a_nNodes)
{

  int error;

  error = xfGetNumberOfNodes(*a_Id, a_nNodes);

  return error;

 } /*xfGetNumberOfNodes_f_*/
 /* --------------------------------------------------------------------------- */
/* FUNCTION  xfReadXNodeLocations_f_ */
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFREADXNODELOCATIONS_F,xfreadxnodelocations_f,xfreadxnodelocations_f_) XFREADXNODELOCATIONS_F
#elif defined FLOWER
#define func_name(XFREADXNODELOCATIONS_F,xfreadxnodelocations_f,xfreadxnodelocations_f_) xfreadxnodelocations_f
#else
#define func_name(XFREADXNODELOCATIONS_F,xfreadxnodelocations_f,xfreadxnodelocations_f_) xfreadxnodelocations_f_
#endif
XMDF_API  int func_name(XFREADXNODELOCATIONS_F,xfreadxnodelocations_f,xfreadxnodelocations_f_)
              (xid *a_Id, int *a_nNodes, double *a_Locs)
{

  int error;

  error = xfReadXNodeLocations(*a_Id, *a_nNodes, a_Locs);

  return error;

 } /*xfReadXNodeLocations_f_*/
 /* --------------------------------------------------------------------------- */
/* FUNCTION  xfReadYNodeLocations_f_ */
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFREADYNODELOCATIONS_F,xfreadynodelocations_f,xfreadynodelocations_f_) XFREADYNODELOCATIONS_F
#elif defined FLOWER
#define func_name(XFREADYNODELOCATIONS_F,xfreadynodelocations_f,xfreadynodelocations_f_) xfreadynodelocations_f
#else
#define func_name(XFREADYNODELOCATIONS_F,xfreadynodelocations_f,xfreadynodelocations_f_) xfreadynodelocations_f_
#endif
XMDF_API  int func_name(XFREADYNODELOCATIONS_F,xfreadynodelocations_f,xfreadynodelocations_f_)
              (xid *a_Id, int *a_nNodes, double *a_Locs)
{

  int error;

  error = xfReadYNodeLocations(*a_Id, *a_nNodes, a_Locs);

  return error;

 } /*xfReadYNodeLocations_f_*/
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfReadZNodeLocations_f_ */
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFREADZNODELOCATIONS_F,xfreadznodelocations_f,xfreadznodelocations_f_) XFREADZNODELOCATIONS_F
#elif defined FLOWER
#define func_name(XFREADZNODELOCATIONS_F,xfreadznodelocations_f,xfreadznodelocations_f_) xfreadznodelocations_f
#else
#define func_name(XFREADZNODELOCATIONS_F,xfreadznodelocations_f,xfreadznodelocations_f_) xfreadznodelocations_f_
#endif
XMDF_API  int func_name(XFREADZNODELOCATIONS_F,xfreadznodelocations_f,xfreadznodelocations_f_)
              (xid *a_Id, int *a_nNodes, double *a_Locs)
{ 

  int error;

  error = xfReadZNodeLocations(*a_Id, *a_nNodes, a_Locs);

  return error;

 } /*xfReadZNodeLocations_f_*/
 /* --------------------------------------------------------------------------*/
/* FUNCTION  xfCreateMeshPropertyGroup_f_                                    */
/* PURPOSE                                                                   */
/* NOTES                                                                     */
/* --------------------------------------------------------------------------*/
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFCREATEMESHPROPERTYGROUP_F,xfcreatemeshpropertygroup_f,xfcreatemeshpropertygroup_f_) XFCREATEMESHPROPERTYGROUP_F
#elif defined FLOWER
#define func_name(XFCREATEMESHPROPERTYGROUP_F,xfcreatemeshpropertygroup_f,xfcreatemeshpropertygroup_f_) xfcreatemeshpropertygroup_f
#else
#define func_name(XFCREATEMESHPROPERTYGROUP_F,xfcreatemeshpropertygroup_f,xfcreatemeshpropertygroup_f_) xfcreatemeshpropertygroup_f_
#endif
XMDF_API  int func_name(XFCREATEMESHPROPERTYGROUP_F,xfcreatemeshpropertygroup_f,xfcreatemeshpropertygroup_f_)
              (xid *a_Id, xid *a_PropId)
{

  int error;

  error = xfCreateMeshPropertyGroup(*a_Id, a_PropId);

  return error;

 } /*xfCreateMeshPropertyGroup_f_*/
 /* --------------------------------------------------------------------------*/
/* FUNCTION  xfGetMeshPropertyGroup_f_                                         */
/* PURPOSE                                                                   */
/* NOTES                                                                     */
/* --------------------------------------------------------------------------*/
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFGETMESHPROPERTYGROUP_F,xfgetmeshpropertygroup_f,xfgetmeshpropertygroup_f_) XFGETMESHPROPERTYGROUP_F
#elif defined FLOWER
#define func_name(XFGETMESHPROPERTYGROUP_F,xfgetmeshpropertygroup_f,xfgetmeshpropertygroup_f_) xfgetmeshpropertygroup_f
#else
#define func_name(XFGETMESHPROPERTYGROUP_F,xfgetmeshpropertygroup_f,xfgetmeshpropertygroup_f_) xfgetmeshpropertygroup_f_
#endif
XMDF_API  int func_name(XFGETMESHPROPERTYGROUP_F,xfgetmeshpropertygroup_f,xfgetmeshpropertygroup_f_)
              (xid *a_Id, xid *a_PropId)
{

  int error;

  error = xfGetMeshPropertyGroup(*a_Id, a_PropId);

  return error;

 } /*xfGetMeshPropertyGroup_f_*/
 /* -------------------------------------------------------------------------- */
/* FUNCTION  xfCreateMeshNodePropertyGroup_f_ */
/* PURPOSE    */
/* NOTES      */
/* -------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFCREATEMESHNODEPROPERTYGROUP_F,xfcreatemeshnodepropertygroup_f,xfcreatemeshnodepropertygroup_f_) XFCREATEMESHNODEPROPERTYGROUP_F
#elif defined FLOWER
#define func_name(XFCREATEMESHNODEPROPERTYGROUP_F,xfcreatemeshnodepropertygroup_f,xfcreatemeshnodepropertygroup_f_) xfcreatemeshnodepropertygroup_f
#else
#define func_name(XFCREATEMESHNODEPROPERTYGROUP_F,xfcreatemeshnodepropertygroup_f,xfcreatemeshnodepropertygroup_f_) xfcreatemeshnodepropertygroup_f_
#endif
XMDF_API  int func_name(XFCREATEMESHNODEPROPERTYGROUP_F,xfcreatemeshnodepropertygroup_f,xfcreatemeshnodepropertygroup_f_)
              (xid *a_Id, xid *a_PropId)
{

  int error;

  error = xfCreateMeshNodePropertyGroup(*a_Id, a_PropId);

  return error;

 } /*xfCreateMeshNodePropertyGroup_f_*/
 /* -------------------------------------------------------------------------- */
/* FUNCTION  xfGetMeshNodePropertyGroup_f_ */
/* PURPOSE    */
/* NOTES      */
/* -------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFGETMESHNODEPROPERTYGROUP_F,xfgetmeshnodepropertygroup_f,xfgetmeshnodepropertygroup_f_) XFGETMESHNODEPROPERTYGROUP_F
#elif defined FLOWER
#define func_name(XFGETMESHNODEPROPERTYGROUP_F,xfgetmeshnodepropertygroup_f,xfgetmeshnodepropertygroup_f_) xfgetmeshnodepropertygroup_f
#else
#define func_name(XFGETMESHNODEPROPERTYGROUP_F,xfgetmeshnodepropertygroup_f,xfgetmeshnodepropertygroup_f_) xfgetmeshnodepropertygroup_f_
#endif
XMDF_API  int func_name(XFGETMESHNODEPROPERTYGROUP_F,xfgetmeshnodepropertygroup_f,xfgetmeshnodepropertygroup_f_)
              (xid *a_Id, xid *a_PropId)
{

  int error;

  error = xfGetMeshNodePropertyGroup(*a_Id, a_PropId);

  return error;

 } /*xfGetMeshNodePropertyGroup_f_*/
 /* -------------------------------------------------------------------------- */
/* FUNCTION  xfCreateMeshElemPropGrp_f_ */
/* PURPOSE    */
/* NOTES      */
/* -------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFCREATEMESHELEMPROPGRP_F,xfcreatemeshelempropgrp_f,xfcreatemeshelempropgrp_f_) XFCREATEMESHELEMPROPGRP_F
#elif defined FLOWER
#define func_name(XFCREATEMESHELEMPROPGRP_F,xfcreatemeshelempropgrp_f,xfcreatemeshelempropgrp_f_) xfcreatemeshelempropgrp_f
#else
#define func_name(XFCREATEMESHELEMPROPGRP_F,xfcreatemeshelempropgrp_f,xfcreatemeshelempropgrp_f_) xfcreatemeshelempropgrp_f_
#endif
XMDF_API  int func_name(XFCREATEMESHELEMPROPGRP_F,xfcreatemeshelempropgrp_f,xfcreatemeshelempropgrp_f_)
              (xid *a_Id, xid *a_PropId)
{

  int error;

  error = xfCreateMeshElementPropertyGroup(*a_Id, a_PropId);

  return error;

 } /*xfCreateMeshElemPropGrp_f_*/
/* -------------------------------------------------------------------------- */
/* FUNCTION  xfGetMeshElementPropertyGroup_f_ */
/* PURPOSE    */
/* NOTES      */
/* -------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFGETMESHELEMENTPROPERTYGROUP_F,xfgetmeshelementpropertygroup_f,xfgetmeshelementpropertygroup_f_) XFGETMESHELEMENTPROPERTYGROUP_F
#elif defined FLOWER
#define func_name(XFGETMESHELEMENTPROPERTYGROUP_F,xfgetmeshelementpropertygroup_f,xfgetmeshelementpropertygroup_f_) xfgetmeshelementpropertygroup_f
#else
#define func_name(XFGETMESHELEMENTPROPERTYGROUP_F,xfgetmeshelementpropertygroup_f,xfgetmeshelementpropertygroup_f_) xfgetmeshelementpropertygroup_f_
#endif
XMDF_API  int func_name(XFGETMESHELEMENTPROPERTYGROUP_F,xfgetmeshelementpropertygroup_f,xfgetmeshelementpropertygroup_f_)
              (xid *a_Id, xid *a_PropId)
{

  int error;

  error = xfGetMeshElementPropertyGroup(*a_Id, a_PropId);

  return error;

 } /*xfGetMeshElementPropertyGroup_f_*/
 /* --------------------------------------------------------------------------- */
/* FUNCTION  xfSetGridType_f_ */
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFSETGRIDTYPE_F,xfsetgridtype_f,xfsetgridtype_f_) XFSETGRIDTYPE_F
#elif defined FLOWER
#define func_name(XFSETGRIDTYPE_F,xfsetgridtype_f,xfsetgridtype_f_) xfsetgridtype_f
#else
#define func_name(XFSETGRIDTYPE_F,xfsetgridtype_f,xfsetgridtype_f_) xfsetgridtype_f_
#endif
XMDF_API xid func_name(XFSETGRIDTYPE_F,xfsetgridtype_f,xfsetgridtype_f_)
              (xid *a_Id, int *a_GridType)
{

  int error;

  error = xfSetGridType(*a_Id, *a_GridType);

  return error;

 } /*xfSetGridType_f_*/
 /* --------------------------------------------------------------------------- */
/* FUNCTION  xfSetNumberOfDimensions_f_ */
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFSETNUMBEROFDIMENSIONS_F,xfsetnumberofdimensions_f,xfsetnumberofdimensions_f_) XFSETNUMBEROFDIMENSIONS_F
#elif defined FLOWER
#define func_name(XFSETNUMBEROFDIMENSIONS_F,xfsetnumberofdimensions_f,xfsetnumberofdimensions_f_) xfsetnumberofdimensions_f
#else
#define func_name(XFSETNUMBEROFDIMENSIONS_F,xfsetnumberofdimensions_f,xfsetnumberofdimensions_f_) xfsetnumberofdimensions_f_
#endif
XMDF_API xid func_name(XFSETNUMBEROFDIMENSIONS_F,xfsetnumberofdimensions_f,xfsetnumberofdimensions_f_)
              (xid *a_Id, int *a_NumDimensions)
{

  int error;

  error = xfSetNumberOfDimensions(*a_Id, *a_NumDimensions);

  return error;

 } /*xfSetNumberOfDimensions_f_*/
 /* --------------------------------------------------------------------------- */
/* FUNCTION  xfSetExtrusionType_f_ */
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFSETEXTRUSIONTYPE_F,xfsetextrusiontype_f,xfsetextrusiontype_f_) XFSETEXTRUSIONTYPE_F
#elif defined FLOWER
#define func_name(XFSETEXTRUSIONTYPE_F,xfsetextrusiontype_f,xfsetextrusiontype_f_) xfsetextrusiontype_f
#else
#define func_name(XFSETEXTRUSIONTYPE_F,xfsetextrusiontype_f,xfsetextrusiontype_f_) xfsetextrusiontype_f_
#endif
XMDF_API xid func_name(XFSETEXTRUSIONTYPE_F,xfsetextrusiontype_f,xfsetextrusiontype_f_)
              (xid *a_Id, int *a_ExtrudeType)
{
  
  int error;

  error = xfSetExtrusionType(*a_Id, *a_ExtrudeType);

  return error;

 } /*xfSetExtrusionType_f_*/
 /* --------------------------------------------------------------------------- */
/* FUNCTION  xfSetOrigin_f_ */
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFSETORIGIN_F,xfsetorigin_f,xfsetorigin_f_) XFSETORIGIN_F
#elif defined FLOWER
#define func_name(XFSETORIGIN_F,xfsetorigin_f,xfsetorigin_f_) xfsetorigin_f
#else
#define func_name(XFSETORIGIN_F,xfsetorigin_f,xfsetorigin_f_) xfsetorigin_f_
#endif
XMDF_API xid func_name(XFSETORIGIN_F,xfsetorigin_f,xfsetorigin_f_)
              (xid *a_Id, double *a_x, double *a_y, double *a_z)
{

  int error;

  error = xfSetOrigin(*a_Id, *a_x, *a_y, *a_z);

  return error;

 } /*xfSetOrigin_f_*/
 /* --------------------------------------------------------------------------- */
/* FUNCTION  xfSetOrientation_f_ */
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFSETORIENTATION_F,xfsetorientation_f,xfsetorientation_f_) XFSETORIENTATION_F
#elif defined FLOWER
#define func_name(XFSETORIENTATION_F,xfsetorientation_f,xfsetorientation_f_) xfsetorientation_f
#else
#define func_name(XFSETORIENTATION_F,xfsetorientation_f,xfsetorientation_f_) xfsetorientation_f_
#endif
XMDF_API xid func_name(XFSETORIENTATION_F,xfsetorientation_f,xfsetorientation_f_)
              (xid *a_Id, int *a_Orientation)
{

  int error;

  error = xfSetOrientation(*a_Id, *a_Orientation);

  return error;

 } /*xfSetOrientation_f_*/
 /* --------------------------------------------------------------------------- */
/* FUNCTION  xfSetBearing_f_ */
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFSETBEARING_F,xfsetbearing_f,xfsetbearing_f_) XFSETBEARING_F
#elif defined FLOWER
#define func_name(XFSETBEARING_F,xfsetbearing_f,xfsetbearing_f_) xfsetbearing_f
#else
#define func_name(XFSETBEARING_F,xfsetbearing_f,xfsetbearing_f_) xfsetbearing_f_
#endif
XMDF_API xid func_name(XFSETBEARING_F,xfsetbearing_f,xfsetbearing_f_)
              (xid *a_Id, double *a_Bearing)
{

  int error;

  error = xfSetBearing(*a_Id, *a_Bearing);

  return error;

 } /*xfSetBearing_f_*/
 /* --------------------------------------------------------------------------- */
/* FUNCTION  xfSetDip_f_ */
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFSETDIP_F,xfsetdip_f,xfsetdip_f_) XFSETDIP_F
#elif defined FLOWER
#define func_name(XFSETDIP_F,xfsetdip_f,xfsetdip_f_) xfsetdip_f
#else
#define func_name(XFSETDIP_F,xfsetdip_f,xfsetdip_f_) xfsetdip_f_
#endif
XMDF_API xid func_name(XFSETDIP_F,xfsetdip_f,xfsetdip_f_)
              (xid *a_Id, double *a_Dip)
{
  
  int error;

  error = xfSetDip(*a_Id, *a_Dip);

  return error;

 } /*xfSetDip_f_*/
 /* --------------------------------------------------------------------------- */
/* FUNCTION  xfSetRoll_f_ */
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFSETROLL_F,xfsetroll_f,xfsetroll_f_) XFSETROLL_F
#elif defined FLOWER
#define func_name(XFSETROLL_F,xfsetroll_f,xfsetroll_f_) xfsetroll_f
#else
#define func_name(XFSETROLL_F,xfsetroll_f,xfsetroll_f_) xfsetroll_f_
#endif
XMDF_API xid func_name(XFSETROLL_F,xfsetroll_f,xfsetroll_f_)
              (xid *a_Id, double *a_Roll)
{

  int error;

  error = xfSetRoll(*a_Id, *a_Roll);

  return error;

 } /*xfSetRoll_f_*/
 /* --------------------------------------------------------------------------- */
/* FUNCTION  xfSetComputationalOrigin_f_ */
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFSETCOMPUTATIONALORIGIN_F,xfsetcomputationalorigin_f,xfsetcomputationalorigin_f_) XFSETCOMPUTATIONALORIGIN_F
#elif defined FLOWER
#define func_name(XFSETCOMPUTATIONALORIGIN_F,xfsetcomputationalorigin_f,xfsetcomputationalorigin_f_) xfsetcomputationalorigin_f
#else
#define func_name(XFSETCOMPUTATIONALORIGIN_F,xfsetcomputationalorigin_f,xfsetcomputationalorigin_f_) xfsetcomputationalorigin_f_
#endif
XMDF_API xid func_name(XFSETCOMPUTATIONALORIGIN_F,xfsetcomputationalorigin_f,xfsetcomputationalorigin_f_)
              (xid *a_Id, int *a_Origin)
{

  int error;

  error = xfSetComputationalOrigin(*a_Id, *a_Origin);

  return error;

 } /*xfSetComputationalOrigin_f_*/
 /* --------------------------------------------------------------------------- */
/* FUNCTION  xfSetUDirection_f_ */
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFSETUDIRECTION_F,xfsetudirection_f,xfsetudirection_f_) XFSETUDIRECTION_F
#elif defined FLOWER
#define func_name(XFSETUDIRECTION_F,xfsetudirection_f,xfsetudirection_f_) xfsetudirection_f
#else
#define func_name(XFSETUDIRECTION_F,xfsetudirection_f,xfsetudirection_f_) xfsetudirection_f_
#endif
XMDF_API xid func_name(XFSETUDIRECTION_F,xfsetudirection_f,xfsetudirection_f_)
              (xid *a_Id, int *a_Direction)
{

  int error;

  error = xfSetUDirection(*a_Id, *a_Direction);

  return error;

 } /*xfSetUDirection_f_*/
 /* --------------------------------------------------------------------------- */
/* FUNCTION  xfSetNumberCellsInI_f_ */
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFSETNUMBERCELLSINI_F,xfsetnumbercellsini_f,xfsetnumbercellsini_f_) XFSETNUMBERCELLSINI_F
#elif defined FLOWER
#define func_name(XFSETNUMBERCELLSINI_F,xfsetnumbercellsini_f,xfsetnumbercellsini_f_) xfsetnumbercellsini_f
#else
#define func_name(XFSETNUMBERCELLSINI_F,xfsetnumbercellsini_f,xfsetnumbercellsini_f_) xfsetnumbercellsini_f_
#endif
XMDF_API xid func_name(XFSETNUMBERCELLSINI_F,xfsetnumbercellsini_f,xfsetnumbercellsini_f_)
              (xid *a_Id, int *a_NumI)
{

  int error;

  error = xfSetNumberCellsInI(*a_Id, *a_NumI);

  return error;

 } /*xfSetNumberCellsInI_f_*/
 /* --------------------------------------------------------------------------- */
/* FUNCTION  xfSetNumberCellsInJ_f_ */
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFSETNUMBERCELLSINJ_F,xfsetnumbercellsinj_f,xfsetnumbercellsinj_f_) XFSETNUMBERCELLSINJ_F
#elif defined FLOWER
#define func_name(XFSETNUMBERCELLSINJ_F,xfsetnumbercellsinj_f,xfsetnumbercellsinj_f_) xfsetnumbercellsinj_f
#else
#define func_name(XFSETNUMBERCELLSINJ_F,xfsetnumbercellsinj_f,xfsetnumbercellsinj_f_) xfsetnumbercellsinj_f_
#endif
XMDF_API xid func_name(XFSETNUMBERCELLSINJ_F,xfsetnumbercellsinj_f,xfsetnumbercellsinj_f_)
              (xid *a_Id, int *a_NumJ)
{

  int error;

  error = xfSetNumberCellsInJ(*a_Id, *a_NumJ);

  return error;

 } /*xfSetNumberCellsInJ_f_*/
 /* --------------------------------------------------------------------------- */
/* FUNCTION  xfSetNumberCellsInK_f_ */
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFSETNUMBERCELLSINK_F,xfsetnumbercellsink_f,xfsetnumbercellsink_f_) XFSETNUMBERCELLSINK_F
#elif defined FLOWER
#define func_name(XFSETNUMBERCELLSINK_F,xfsetnumbercellsink_f,xfsetnumbercellsink_f_) xfsetnumbercellsink_f
#else
#define func_name(XFSETNUMBERCELLSINK_F,xfsetnumbercellsink_f,xfsetnumbercellsink_f_) xfsetnumbercellsink_f_
#endif
XMDF_API xid func_name(XFSETNUMBERCELLSINK_F,xfsetnumbercellsink_f,xfsetnumbercellsink_f_)
              (xid *a_Id, int *a_NumK)
{

  int error;

  error = xfSetNumberCellsInK(*a_Id, *a_NumK);

  return error;

 } /*xfSetNumberCellsInK_f_*/
 /*-----------------------------------------------------------------------------*/
/* FUNCTION  xfSetGridCoordsI_f_ */
/* PURPOSE*/
/* NOTES*/
/*-----------------------------------------------------------------------------*/
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFSETGRIDCOORDSI_F,xfsetgridcoordsi_f,xfsetgridcoordsi_f_) XFSETGRIDCOORDSI_F
#elif defined FLOWER
#define func_name(XFSETGRIDCOORDSI_F,xfsetgridcoordsi_f,xfsetgridcoordsi_f_) xfsetgridcoordsi_f
#else
#define func_name(XFSETGRIDCOORDSI_F,xfsetgridcoordsi_f,xfsetgridcoordsi_f_) xfsetgridcoordsi_f_
#endif
XMDF_API xid func_name(XFSETGRIDCOORDSI_F,xfsetgridcoordsi_f,xfsetgridcoordsi_f_)
              (xid *a_Id, int *a_NumVals, double *a_iValues)
{

  int error;

  error = xfSetGridCoordsI(*a_Id, *a_NumVals, a_iValues);

  return error;

 } /*xfSetGridCoordsI_f_*/
 /* --------------------------------------------------------------------------- */
/* FUNCTION  xfSetGridCoordsJ_f_ */
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFSETGRIDCOORDSJ_F,xfsetgridcoordsj_f,xfsetgridcoordsj_f_) XFSETGRIDCOORDSJ_F
#elif defined FLOWER
#define func_name(XFSETGRIDCOORDSJ_F,xfsetgridcoordsj_f,xfsetgridcoordsj_f_) xfsetgridcoordsj_f
#else
#define func_name(XFSETGRIDCOORDSJ_F,xfsetgridcoordsj_f,xfsetgridcoordsj_f_) xfsetgridcoordsj_f_
#endif
XMDF_API xid func_name(XFSETGRIDCOORDSJ_F,xfsetgridcoordsj_f,xfsetgridcoordsj_f_)
              (xid *a_Id, int *a_NumVals, double *a_jValues)
{

  int error;

  error = xfSetGridCoordsJ(*a_Id, *a_NumVals, a_jValues);

  return error;

 } /*xfSetGridCoordsJ_f_*/
 /* --------------------------------------------------------------------------- */
/* FUNCTION  xfSetGridCoordsK_f_ */
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFSETGRIDCOORDSK_F,xfsetgridcoordsk_f,xfsetgridcoordsk_f_) XFSETGRIDCOORDSK_F
#elif defined FLOWER
#define func_name(XFSETGRIDCOORDSK_F,xfsetgridcoordsk_f,xfsetgridcoordsk_f_) xfsetgridcoordsk_f
#else
#define func_name(XFSETGRIDCOORDSK_F,xfsetgridcoordsk_f,xfsetgridcoordsk_f_) xfsetgridcoordsk_f_
#endif
XMDF_API xid func_name(XFSETGRIDCOORDSK_F,xfsetgridcoordsk_f,xfsetgridcoordsk_f_)
              (xid *a_Id, int *a_NumVals, double *a_kValues)
{

  int error;

  error = xfSetGridCoordsK(*a_Id, *a_NumVals, a_kValues);

  return error;

 } /*xfSetGridCoordsK_f_*/
 /* --------------------------------------------------------------------------- */
/* FUNCTION  xfWriteExtrudeLayerData_f_ */
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFWRITEEXTRUDELAYERDATA_F,xfwriteextrudelayerdata_f,xfwriteextrudelayerdata_f_) XFWRITEEXTRUDELAYERDATA_F
#elif defined FLOWER
#define func_name(XFWRITEEXTRUDELAYERDATA_F,xfwriteextrudelayerdata_f,xfwriteextrudelayerdata_f_) xfwriteextrudelayerdata_f
#else
#define func_name(XFWRITEEXTRUDELAYERDATA_F,xfwriteextrudelayerdata_f,xfwriteextrudelayerdata_f_) xfwriteextrudelayerdata_f_
#endif
XMDF_API xid func_name(XFWRITEEXTRUDELAYERDATA_F,xfwriteextrudelayerdata_f,xfwriteextrudelayerdata_f_)
              (xid *a_Id, int *a_NumLayers, int *a_NumVals,
                                double *a_Values)
{

  int error;

  error = xfWriteExtrudeLayerData(*a_Id, *a_NumLayers, *a_NumVals, a_Values);

  return error;

 } /*xfWriteExtrudeLayerData_f_*/
 /* --------------------------------------------------------------------------- */
/* FUNCTION  xfGetGridType_f_ */
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFGETGRIDTYPE_F,xfgetgridtype_f,xfgetgridtype_f_) XFGETGRIDTYPE_F
#elif defined FLOWER
#define func_name(XFGETGRIDTYPE_F,xfgetgridtype_f,xfgetgridtype_f_) xfgetgridtype_f
#else
#define func_name(XFGETGRIDTYPE_F,xfgetgridtype_f,xfgetgridtype_f_) xfgetgridtype_f_
#endif
XMDF_API xid func_name(XFGETGRIDTYPE_F,xfgetgridtype_f,xfgetgridtype_f_)
              (xid *a_Id, int *a_GridType)
{

  int error;

  error = xfGetGridType(*a_Id, a_GridType);

  return error;

 } /*xfGetGridType_f_*/
 /* --------------------------------------------------------------------------- */
/* FUNCTION  xfGetExtrusionType_f_ */
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFGETEXTRUSIONTYPE_F,xfgetextrusiontype_f,xfgetextrusiontype_f_) XFGETEXTRUSIONTYPE_F
#elif defined FLOWER
#define func_name(XFGETEXTRUSIONTYPE_F,xfgetextrusiontype_f,xfgetextrusiontype_f_) xfgetextrusiontype_f
#else
#define func_name(XFGETEXTRUSIONTYPE_F,xfgetextrusiontype_f,xfgetextrusiontype_f_) xfgetextrusiontype_f_
#endif
XMDF_API xid func_name(XFGETEXTRUSIONTYPE_F,xfgetextrusiontype_f,xfgetextrusiontype_f_)
              (xid *a_Id, int *a_ExtrudeType)
{

  int error;

  error = xfGetExtrusionType(*a_Id, a_ExtrudeType);

  return error;

 } /*xfGetExtrusionType_f_*/
 /* --------------------------------------------------------------------------- */
/* FUNCTION  xfGetNumberOfDimensions_f_ */
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFGETNUMBEROFDIMENSIONS_F,xfgetnumberofdimensions_f,xfgetnumberofdimensions_f_) XFGETNUMBEROFDIMENSIONS_F
#elif defined FLOWER
#define func_name(XFGETNUMBEROFDIMENSIONS_F,xfgetnumberofdimensions_f,xfgetnumberofdimensions_f_) xfgetnumberofdimensions_f
#else
#define func_name(XFGETNUMBEROFDIMENSIONS_F,xfgetnumberofdimensions_f,xfgetnumberofdimensions_f_) xfgetnumberofdimensions_f_
#endif
XMDF_API xid func_name(XFGETNUMBEROFDIMENSIONS_F,xfgetnumberofdimensions_f,xfgetnumberofdimensions_f_)
              (xid *a_Id, int *a_NumDimensions)
{

  int error;

  error = xfGetNumberOfDimensions(*a_Id, a_NumDimensions);

  return error;

 } /*xfGetNumberOfDimensions_f_*/
 /* --------------------------------------------------------------------------- */
/* FUNCTION  xfOriginDefined_f_ */
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFORIGINDEFINED_F,xforigindefined_f,xforigindefined_f_) XFORIGINDEFINED_F
#elif defined FLOWER
#define func_name(XFORIGINDEFINED_F,xforigindefined_f,xforigindefined_f_) xforigindefined_f
#else
#define func_name(XFORIGINDEFINED_F,xforigindefined_f,xforigindefined_f_) xforigindefined_f_
#endif
XMDF_API xid func_name(XFORIGINDEFINED_F,xforigindefined_f,xforigindefined_f_)
              (xid *a_Id, short *a_bDefined)
{

  xmbool  bDefined;
  int     error;

  bDefined = (*a_bDefined) ? XTRUE : XFALSE;
  error = xfOriginDefined(*a_Id, &bDefined);

  return error;

 } /*xfOriginDefined_f_*/
 /* --------------------------------------------------------------------------- */
/* FUNCTION  xfGetOrigin_f_ */
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFGETORIGIN_F,xfgetorigin_f,xfgetorigin_f_) XFGETORIGIN_F
#elif defined FLOWER
#define func_name(XFGETORIGIN_F,xfgetorigin_f,xfgetorigin_f_) xfgetorigin_f
#else
#define func_name(XFGETORIGIN_F,xfgetorigin_f,xfgetorigin_f_) xfgetorigin_f_
#endif
XMDF_API xid func_name(XFGETORIGIN_F,xfgetorigin_f,xfgetorigin_f_)
              (xid *a_Id, double *a_x, double *a_y, double *a_z)
{

  int error;

  error = xfGetOrigin(*a_Id, a_x, a_y, a_z);

  return error;

 } /*xfGetOrigin_f_*/
 /* --------------------------------------------------------------------------- */
/* FUNCTION  xfGetOrientation_f_ */
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFGETORIENTATION_F,xfgetorientation_f,xfgetorientation_f_) XFGETORIENTATION_F
#elif defined FLOWER
#define func_name(XFGETORIENTATION_F,xfgetorientation_f,xfgetorientation_f_) xfgetorientation_f
#else
#define func_name(XFGETORIENTATION_F,xfgetorientation_f,xfgetorientation_f_) xfgetorientation_f_
#endif
XMDF_API xid func_name(XFGETORIENTATION_F,xfgetorientation_f,xfgetorientation_f_)
              (xid *a_Id, int *a_Orientation)
{

  int error;

  error = xfGetOrientation(*a_Id, a_Orientation);

  return error;

 } /*xfGetOrientation_f_*/
 /* --------------------------------------------------------------------------- */
/* FUNCTION  xfBearingDefined_f_ */
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFBEARINGDEFINED_F,xfbearingdefined_f,xfbearingdefined_f_) XFBEARINGDEFINED_F
#elif defined FLOWER
#define func_name(XFBEARINGDEFINED_F,xfbearingdefined_f,xfbearingdefined_f_) xfbearingdefined_f
#else
#define func_name(XFBEARINGDEFINED_F,xfbearingdefined_f,xfbearingdefined_f_) xfbearingdefined_f_
#endif
XMDF_API xid func_name(XFBEARINGDEFINED_F,xfbearingdefined_f,xfbearingdefined_f_)
              (xid *a_Id, short *a_bDefined)
{

  xmbool  bDefined;
  int     error;

  bDefined = (*a_bDefined) ? XTRUE : XFALSE;
  error = xfBearingDefined(*a_Id, &bDefined);

  return error;

 } /*xfBearingDefined_f_*/
 /* --------------------------------------------------------------------------- */
/* FUNCTION  xfGetBearing_f_ */
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFGETBEARING_F,xfgetbearing_f,xfgetbearing_f_) XFGETBEARING_F
#elif defined FLOWER
#define func_name(XFGETBEARING_F,xfgetbearing_f,xfgetbearing_f_) xfgetbearing_f
#else
#define func_name(XFGETBEARING_F,xfgetbearing_f,xfgetbearing_f_) xfgetbearing_f_
#endif
XMDF_API xid func_name(XFGETBEARING_F,xfgetbearing_f,xfgetbearing_f_)
              (xid *a_Id, double *a_bearing)
{

  int error;

  error = xfGetBearing(*a_Id, a_bearing);

  return error;

 } /*xfGetBearing_f_*/
 /* --------------------------------------------------------------------------- */
/* FUNCTION  xfDipDefined_f_ */
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFDIPDEFINED_F,xfdipdefined_f,xfdipdefined_f_) XFDIPDEFINED_F
#elif defined FLOWER
#define func_name(XFDIPDEFINED_F,xfdipdefined_f,xfdipdefined_f_) xfdipdefined_f
#else
#define func_name(XFDIPDEFINED_F,xfdipdefined_f,xfdipdefined_f_) xfdipdefined_f_
#endif
XMDF_API xid func_name(XFDIPDEFINED_F,xfdipdefined_f,xfdipdefined_f_)
              (xid *a_Id, short *a_bDefined)
{

  xmbool  bDefined;
  int     error;

  bDefined = (*a_bDefined) ? XTRUE : XFALSE;
  error = xfDipDefined(*a_Id, &bDefined);

  return error;

 } /*xfDipDefined_f_*/
 /* --------------------------------------------------------------------------- */
/* FUNCTION  xfGetDip_f_ */
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFGETDIP_F,xfgetdip_f,xfgetdip_f_) XFGETDIP_F
#elif defined FLOWER
#define func_name(XFGETDIP_F,xfgetdip_f,xfgetdip_f_) xfgetdip_f
#else
#define func_name(XFGETDIP_F,xfgetdip_f,xfgetdip_f_) xfgetdip_f_
#endif
XMDF_API xid func_name(XFGETDIP_F,xfgetdip_f,xfgetdip_f_)
              (xid *a_Id, double *a_dip)
{
  int error;

  error = xfGetDip(*a_Id, a_dip);

  return error;

} /* xfGetDip_f_*/
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfRollDefined_f_*/
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFROLLDEFINED_F,xfrolldefined_f,xfrolldefined_f_) XFROLLDEFINED_F
#elif defined FLOWER
#define func_name(XFROLLDEFINED_F,xfrolldefined_f,xfrolldefined_f_) xfrolldefined_f
#else
#define func_name(XFROLLDEFINED_F,xfrolldefined_f,xfrolldefined_f_) xfrolldefined_f_
#endif
XMDF_API xid func_name(XFROLLDEFINED_F,xfrolldefined_f,xfrolldefined_f_)
              (xid *a_Id, short *a_bDefined)
{
  xmbool  bDefined;  
  int     error;

  bDefined = (*a_bDefined) ? XTRUE : XFALSE;
  error = xfRollDefined(*a_Id, &bDefined);

  return error;
} /* xfRollDefined_f_*/
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfGetRoll_f_*/
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFGETROLL_F,xfgetroll_f,xfgetroll_f_) XFGETROLL_F
#elif defined FLOWER
#define func_name(XFGETROLL_F,xfgetroll_f,xfgetroll_f_) xfgetroll_f
#else
#define func_name(XFGETROLL_F,xfgetroll_f,xfgetroll_f_) xfgetroll_f_
#endif
XMDF_API xid func_name(XFGETROLL_F,xfgetroll_f,xfgetroll_f_)
              (xid *a_Id, double *a_Roll)
{
  int error;

  error = xfGetRoll(*a_Id, a_Roll);

  return error;

} /* xfGetRoll_f_*/
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfComputationalOriginDefined_f_*/
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFCOMPUTATIONALORIGINDEFINED_F,xfcomputationalorigindefined_f,xfcomputationalorigindefined_f_) XFCOMPUTATIONALORIGINDEFINED_F
#elif defined FLOWER
#define func_name(XFCOMPUTATIONALORIGINDEFINED_F,xfcomputationalorigindefined_f,xfcomputationalorigindefined_f_) xfcomputationalorigindefined_f
#else
#define func_name(XFCOMPUTATIONALORIGINDEFINED_F,xfcomputationalorigindefined_f,xfcomputationalorigindefined_f_) xfcomputationalorigindefined_f_
#endif
XMDF_API xid func_name(XFCOMPUTATIONALORIGINDEFINED_F,xfcomputationalorigindefined_f,xfcomputationalorigindefined_f_)
              (xid *GroupId, short *bDefined)
{
  xmbool  Defined;
  int     error;

  Defined = (*bDefined) ? XTRUE : XFALSE;
  error = xfComputationalOriginDefined(*GroupId, &Defined);

  return error;

} /* xfComputationalOriginDefined_f_*/
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfGetComputationalOrigin_f_*/
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFGETCOMPUTATIONALORIGIN_F,xfgetcomputationalorigin_f,xfgetcomputationalorigin_f_) XFGETCOMPUTATIONALORIGIN_F
#elif defined FLOWER
#define func_name(XFGETCOMPUTATIONALORIGIN_F,xfgetcomputationalorigin_f,xfgetcomputationalorigin_f_) xfgetcomputationalorigin_f
#else
#define func_name(XFGETCOMPUTATIONALORIGIN_F,xfgetcomputationalorigin_f,xfgetcomputationalorigin_f_) xfgetcomputationalorigin_f_
#endif
XMDF_API xid func_name(XFGETCOMPUTATIONALORIGIN_F,xfgetcomputationalorigin_f,xfgetcomputationalorigin_f_)
              (xid *GroupId, int *Origin)
{
  int error;

  error = xfGetComputationalOrigin(*GroupId, Origin);

  return error;

} /* xfGetComputationalOrigin_f_*/
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfGetUDirectionDefined_f_*/
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFGETUDIRECTIONDEFINED_F,xfgetudirectiondefined_f,xfgetudirectiondefined_f_) XFGETUDIRECTIONDEFINED_F
#elif defined FLOWER
#define func_name(XFGETUDIRECTIONDEFINED_F,xfgetudirectiondefined_f,xfgetudirectiondefined_f_) xfgetudirectiondefined_f
#else
#define func_name(XFGETUDIRECTIONDEFINED_F,xfgetudirectiondefined_f,xfgetudirectiondefined_f_) xfgetudirectiondefined_f_
#endif
XMDF_API xid func_name(XFGETUDIRECTIONDEFINED_F,xfgetudirectiondefined_f,xfgetudirectiondefined_f_)
              (xid *GroupId, short *bDefined)
{
  
  xmbool    Defined;
  int       error;

  Defined = (*bDefined) ? XTRUE : XFALSE;
  error = xfGetUDirectionDefined(*GroupId, &Defined);

  return error;

} /* xfGetUDirectionDefined_f_*/
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfGetUDirection_f_*/
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFGETUDIRECTION_F,xfgetudirection_f,xfgetudirection_f_) XFGETUDIRECTION_F
#elif defined FLOWER
#define func_name(XFGETUDIRECTION_F,xfgetudirection_f,xfgetudirection_f_) xfgetudirection_f
#else
#define func_name(XFGETUDIRECTION_F,xfgetudirection_f,xfgetudirection_f_) xfgetudirection_f_
#endif
XMDF_API xid func_name(XFGETUDIRECTION_F,xfgetudirection_f,xfgetudirection_f_)
              (xid *GroupId, int *Direction)
{
  int error;

  error = xfGetUDirection(*GroupId, Direction);

  return error;

} /* xfGetUDirection_f_*/
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfGetNumberCellsInI_f_*/
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFGETNUMBERCELLSINI_F,xfgetnumbercellsini_f,xfgetnumbercellsini_f_) XFGETNUMBERCELLSINI_F
#elif defined FLOWER
#define func_name(XFGETNUMBERCELLSINI_F,xfgetnumbercellsini_f,xfgetnumbercellsini_f_) xfgetnumbercellsini_f
#else
#define func_name(XFGETNUMBERCELLSINI_F,xfgetnumbercellsini_f,xfgetnumbercellsini_f_) xfgetnumbercellsini_f_
#endif
XMDF_API xid func_name(XFGETNUMBERCELLSINI_F,xfgetnumbercellsini_f,xfgetnumbercellsini_f_)
              (xid *a_Id, int *a_NumI)
{
  int error;

  error = xfGetNumberCellsInI(*a_Id, a_NumI);

  return error;

} /* xfGetNumberCellsInI_f_*/
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfGetNumberCellsInJ_f_*/
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFGETNUMBERCELLSINJ_F,xfgetnumbercellsinj_f,xfgetnumbercellsinj_f_) XFGETNUMBERCELLSINJ_F
#elif defined FLOWER
#define func_name(XFGETNUMBERCELLSINJ_F,xfgetnumbercellsinj_f,xfgetnumbercellsinj_f_) xfgetnumbercellsinj_f
#else
#define func_name(XFGETNUMBERCELLSINJ_F,xfgetnumbercellsinj_f,xfgetnumbercellsinj_f_) xfgetnumbercellsinj_f_
#endif
XMDF_API xid func_name(XFGETNUMBERCELLSINJ_F,xfgetnumbercellsinj_f,xfgetnumbercellsinj_f_)
              (xid *a_Id, int *a_NumJ)
{
  int error;

  error = xfGetNumberCellsInJ(*a_Id, a_NumJ);

  return error;

} /* xfGetNumberCellsInJ_f_*/
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfGetNumberCellsInK_f_*/
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFGETNUMBERCELLSINK_F,xfgetnumbercellsink_f,xfgetnumbercellsink_f_) XFGETNUMBERCELLSINK_F
#elif defined FLOWER
#define func_name(XFGETNUMBERCELLSINK_F,xfgetnumbercellsink_f,xfgetnumbercellsink_f_) xfgetnumbercellsink_f
#else
#define func_name(XFGETNUMBERCELLSINK_F,xfgetnumbercellsink_f,xfgetnumbercellsink_f_) xfgetnumbercellsink_f_
#endif
XMDF_API xid func_name(XFGETNUMBERCELLSINK_F,xfgetnumbercellsink_f,xfgetnumbercellsink_f_)
              (xid *a_Id, int *a_NumK)
{
  int error;

  error = xfGetNumberCellsInK(*a_Id, a_NumK);

  return error;

} /* xfGetNumberCellsInK_f_*/
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfGetGridCoordsI_f_*/
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFGETGRIDCOORDSI_F,xfgetgridcoordsi_f,xfgetgridcoordsi_f_) XFGETGRIDCOORDSI_F
#elif defined FLOWER
#define func_name(XFGETGRIDCOORDSI_F,xfgetgridcoordsi_f,xfgetgridcoordsi_f_) xfgetgridcoordsi_f
#else
#define func_name(XFGETGRIDCOORDSI_F,xfgetgridcoordsi_f,xfgetgridcoordsi_f_) xfgetgridcoordsi_f_
#endif
XMDF_API xid func_name(XFGETGRIDCOORDSI_F,xfgetgridcoordsi_f,xfgetgridcoordsi_f_)
              (xid *a_Id, int *a_NumVals, double *a_iValues)
{
  int error;

  error = xfGetGridCoordsI(*a_Id, *a_NumVals, a_iValues);

  return error;

} /* xfGetGridCoordsI_f_*/
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfGetGridCoordsJ_f_*/
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFGETGRIDCOORDSJ_F,xfgetgridcoordsj_f,xfgetgridcoordsj_f_) XFGETGRIDCOORDSJ_F
#elif defined FLOWER
#define func_name(XFGETGRIDCOORDSJ_F,xfgetgridcoordsj_f,xfgetgridcoordsj_f_) xfgetgridcoordsj_f
#else
#define func_name(XFGETGRIDCOORDSJ_F,xfgetgridcoordsj_f,xfgetgridcoordsj_f_) xfgetgridcoordsj_f_
#endif
XMDF_API xid func_name(XFGETGRIDCOORDSJ_F,xfgetgridcoordsj_f,xfgetgridcoordsj_f_)
              (xid *a_Id, int *a_NumVals, double *a_jValues)
{
  int error;

  error = xfGetGridCoordsJ(*a_Id, *a_NumVals, a_jValues);

  return error;

} /* xfGetGridCoordsJ_f_*/
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfGetGridCoordsK_f_*/
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFGETGRIDCOORDSK_F,xfgetgridcoordsk_f,xfgetgridcoordsk_f_) XFGETGRIDCOORDSK_F
#elif defined FLOWER
#define func_name(XFGETGRIDCOORDSK_F,xfgetgridcoordsk_f,xfgetgridcoordsk_f_) xfgetgridcoordsk_f
#else
#define func_name(XFGETGRIDCOORDSK_F,xfgetgridcoordsk_f,xfgetgridcoordsk_f_) xfgetgridcoordsk_f_
#endif
XMDF_API xid func_name(XFGETGRIDCOORDSK_F,xfgetgridcoordsk_f,xfgetgridcoordsk_f_)
              (xid *a_Id, int *a_NumVals, double *a_kValues)
{
  int error;

  error = xfGetGridCoordsK(*a_Id, *a_NumVals, a_kValues);

  return error;

} /* xfGetGridCoordsK_f_*/
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfGetExtrudeNumLayers_f_*/
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFGETEXTRUDENUMLAYERS_F,xfgetextrudenumlayers_f,xfgetextrudenumlayers_f_) XFGETEXTRUDENUMLAYERS_F
#elif defined FLOWER
#define func_name(XFGETEXTRUDENUMLAYERS_F,xfgetextrudenumlayers_f,xfgetextrudenumlayers_f_) xfgetextrudenumlayers_f
#else
#define func_name(XFGETEXTRUDENUMLAYERS_F,xfgetextrudenumlayers_f,xfgetextrudenumlayers_f_) xfgetextrudenumlayers_f_
#endif
XMDF_API xid func_name(XFGETEXTRUDENUMLAYERS_F,xfgetextrudenumlayers_f,xfgetextrudenumlayers_f_)
              (xid *a_Id, int *a_NumLayers)
{
  int error;

  error = xfGetExtrudeNumLayers(*a_Id, a_NumLayers);

  return error;

} /* xfGetExtrudeNumLayers_f_*/
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfgetextrudevalues_f_*/
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFGETEXTRUDEVALUES_F,xfgetextrudevalues_f,xfgetextrudevalues_f_) XFGETEXTRUDEVALUES_F
#elif defined FLOWER
#define func_name(XFGETEXTRUDEVALUES_F,xfgetextrudevalues_f,xfgetextrudevalues_f_) xfgetextrudevalues_f
#else
#define func_name(XFGETEXTRUDEVALUES_F,xfgetextrudevalues_f,xfgetextrudevalues_f_) xfgetextrudevalues_f_
#endif
XMDF_API xid func_name(XFGETEXTRUDEVALUES_F,xfgetextrudevalues_f,xfgetextrudevalues_f_)
              (xid *a_Id, int *a_NumVals, double *a_Values)
{
  int error;

  error = xfGetExtrudeValues(*a_Id, *a_NumVals, a_Values);

  return error;

} /* xfGetExtrudeValues_f_*/
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfCreateGridPropertyGroup_f_                                        */
/* PURPOSE                                                                     */
/* NOTES                                                                       */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFCREATEGRIDPROPERTYGROUP_F,xfcreategridpropertygroup_f,xfcreategridpropertygroup_f_) XFCREATEGRIDPROPERTYGROUP_F
#elif defined FLOWER
#define func_name(XFCREATEGRIDPROPERTYGROUP_F,xfcreategridpropertygroup_f,xfcreategridpropertygroup_f_) xfcreategridpropertygroup_f
#else
#define func_name(XFCREATEGRIDPROPERTYGROUP_F,xfcreategridpropertygroup_f,xfcreategridpropertygroup_f_) xfcreategridpropertygroup_f_
#endif
XMDF_API xid func_name(XFCREATEGRIDPROPERTYGROUP_F,xfcreategridpropertygroup_f,xfcreategridpropertygroup_f_)
              (xid *a_Id, xid *a_PropId)
{
  int error;

  error = xfCreateGridPropertyGroup(*a_Id, a_PropId);

  return error;

} /* xfCreateGridPropertyGroup_f_*/
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfGetGridPropertyGroup_f_                                           */
/* PURPOSE                                                                     */
/* NOTES                                                                       */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFGETGRIDPROPERTYGROUP_F,xfgetgridpropertygroup_f,xfgetgridpropertygroup_f_) XFGETGRIDPROPERTYGROUP_F
#elif defined FLOWER
#define func_name(XFGETGRIDPROPERTYGROUP_F,xfgetgridpropertygroup_f,xfgetgridpropertygroup_f_) xfgetgridpropertygroup_f
#else
#define func_name(XFGETGRIDPROPERTYGROUP_F,xfgetgridpropertygroup_f,xfgetgridpropertygroup_f_) xfgetgridpropertygroup_f_
#endif
XMDF_API xid func_name(XFGETGRIDPROPERTYGROUP_F,xfgetgridpropertygroup_f,xfgetgridpropertygroup_f_)
              (xid *a_Id, xid *a_PropId)
{
  int error;

  error = xfGetGridPropertyGroup(*a_Id, a_PropId);

  return error;

} /* xfGetGridPropertyGroup_f_*/
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfCreateGridCellPropertyGroup_f_                                   */
/* PURPOSE                                                                     */
/* NOTES                                                                       */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFCREATEGRIDCELLPROPERTYGROUP_F,xfcreategridcellpropertygroup_f,xfcreategridcellpropertygroup_f_) XFCREATEGRIDCELLPROPERTYGROUP_F
#elif defined FLOWER
#define func_name(XFCREATEGRIDCELLPROPERTYGROUP_F,xfcreategridcellpropertygroup_f,xfcreategridcellpropertygroup_f_) xfcreategridcellpropertygroup_f
#else
#define func_name(XFCREATEGRIDCELLPROPERTYGROUP_F,xfcreategridcellpropertygroup_f,xfcreategridcellpropertygroup_f_) xfcreategridcellpropertygroup_f_
#endif
XMDF_API xid func_name(XFCREATEGRIDCELLPROPERTYGROUP_F,xfcreategridcellpropertygroup_f,xfcreategridcellpropertygroup_f_)
              (xid *a_Id, xid *a_PropId)
{
  int error;

  error = xfCreateGridCellPropertyGroup(*a_Id, a_PropId);

  return error;

} /* xfCreateGridCellPropertyGroup_f_*/
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfGetGridCellPropertyGroup_f_                                       */
/* PURPOSE                                                                     */
/* NOTES                                                                       */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFGETGRIDCELLPROPERTYGROUP_F,xfgetgridcellpropertygroup_f,xfgetgridcellpropertygroup_f_) XFGETGRIDCELLPROPERTYGROUP_F
#elif defined FLOWER
#define func_name(XFGETGRIDCELLPROPERTYGROUP_F,xfgetgridcellpropertygroup_f,xfgetgridcellpropertygroup_f_) xfgetgridcellpropertygroup_f
#else
#define func_name(XFGETGRIDCELLPROPERTYGROUP_F,xfgetgridcellpropertygroup_f,xfgetgridcellpropertygroup_f_) xfgetgridcellpropertygroup_f_
#endif
XMDF_API xid func_name(XFGETGRIDCELLPROPERTYGROUP_F,xfgetgridcellpropertygroup_f,xfgetgridcellpropertygroup_f_)
              (xid *a_Id, xid *a_PropId)
{
  int error;

  error = xfGetGridCellPropertyGroup(*a_Id, a_PropId);

  return error;

} /* xfGetGridCellPropertyGroup_f_*/
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfCreateGridNodePropertyGroup_f_                                    */
/* PURPOSE                                                                     */
/* NOTES                                                                       */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFCREATEGRIDNODEPROPERTYGROUP_F,xfcreategridnodepropertygroup_f,xfcreategridnodepropertygroup_f_) XFCREATEGRIDNODEPROPERTYGROUP_F
#elif defined FLOWER
#define func_name(XFCREATEGRIDNODEPROPERTYGROUP_F,xfcreategridnodepropertygroup_f,xfcreategridnodepropertygroup_f_) xfcreategridnodepropertygroup_f
#else
#define func_name(XFCREATEGRIDNODEPROPERTYGROUP_F,xfcreategridnodepropertygroup_f,xfcreategridnodepropertygroup_f_) xfcreategridnodepropertygroup_f_
#endif
XMDF_API xid func_name(XFCREATEGRIDNODEPROPERTYGROUP_F,xfcreategridnodepropertygroup_f,xfcreategridnodepropertygroup_f_)
              (xid *a_Id, xid *a_PropId)
{
  int error;

  error = xfCreateGridNodePropertyGroup(*a_Id, a_PropId);

  return error;

} /* xfCreateGridNodePropertyGroup_f_*/
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfGetGridNodePropertyGroup_f_                                       */
/* PURPOSE                                                                     */
/* NOTES                                                                       */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFGETGRIDNODEPROPERTYGROUP_F,xfgetgridnodepropertygroup_f,xfgetgridnodepropertygroup_f_) XFGETGRIDNODEPROPERTYGROUP_F
#elif defined FLOWER
#define func_name(XFGETGRIDNODEPROPERTYGROUP_F,xfgetgridnodepropertygroup_f,xfgetgridnodepropertygroup_f_) xfgetgridnodepropertygroup_f
#else
#define func_name(XFGETGRIDNODEPROPERTYGROUP_F,xfgetgridnodepropertygroup_f,xfgetgridnodepropertygroup_f_) xfgetgridnodepropertygroup_f_
#endif
XMDF_API xid func_name(XFGETGRIDNODEPROPERTYGROUP_F,xfgetgridnodepropertygroup_f,xfgetgridnodepropertygroup_f_)
              (xid *a_Id, xid *a_PropId)
{
  int error;

  error = xfGetGridNodePropertyGroup(*a_Id, a_PropId);

  return error;

} /* xfGetGridNodePropertyGroup_f_*/
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfCreateGeometricPathGroup_f_*/
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFCREATEGEOMETRICPATHGROUP_F,xfcreategeometricpathgroup_f,xfcreategeometricpathgroup_f_) XFCREATEGEOMETRICPATHGROUP_F
#elif defined FLOWER
#define func_name(XFCREATEGEOMETRICPATHGROUP_F,xfcreategeometricpathgroup_f,xfcreategeometricpathgroup_f_) xfcreategeometricpathgroup_f
#else
#define func_name(XFCREATEGEOMETRICPATHGROUP_F,xfcreategeometricpathgroup_f,xfcreategeometricpathgroup_f_) xfcreategeometricpathgroup_f_
#endif
XMDF_API xid func_name(XFCREATEGEOMETRICPATHGROUP_F,xfcreategeometricpathgroup_f,xfcreategeometricpathgroup_f_)
              (xid *a_ParentId, const char *a_Path, int *a_pathlen,
                                   const char *a_Guid, int *a_guidlen, int *a_Compression,
                                   xid *a_PathGroup, double *a_NullVal)
{
  int	error;
  char *path_copy = NULL, *guid_copy = NULL;

  path_copy = (char*)malloc((*a_pathlen+1)*sizeof(char));
  guid_copy = (char*)malloc((*a_guidlen+1)*sizeof(char));
  strncpy(path_copy, a_Path, *a_pathlen);
  strncpy(guid_copy, a_Guid, *a_guidlen);
  path_copy[*a_pathlen] = '\0';
  guid_copy[*a_guidlen] = '\0';

  error = xfCreateGeometricPathGroup(*a_ParentId, path_copy, guid_copy, *a_Compression,
                                     a_PathGroup, *a_NullVal);

  free(path_copy);
  free(guid_copy);
  return error;

} /* xfCreateGeometricPathGroup_f_*/

/* --------------------------------------------------------------------------- */
/* FUNCTION  xfWriteParticleTimestep_f_*/
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFWRITEPARTICLETIMESTEP_F,xfwriteparticletimestep_f,xfwriteparticletimestep_f_) XFWRITEPARTICLETIMESTEP_F
#elif defined FLOWER
#define func_name(XFWRITEPARTICLETIMESTEP_F,xfwriteparticletimestep_f,xfwriteparticletimestep_f_) xfwriteparticletimestep_f
#else
#define func_name(XFWRITEPARTICLETIMESTEP_F,xfwriteparticletimestep_f,xfwriteparticletimestep_f_) xfwriteparticletimestep_f_
#endif
XMDF_API xid func_name(XFWRITEPARTICLETIMESTEP_F,xfwriteparticletimestep_f,xfwriteparticletimestep_f_)
              (xid *a_Id, int *a_nDim, double *a_Time,
                                int *a_nPaths, double *a_Locs)
{
  int error;

  error = xfWriteParticleTimestep(*a_Id, *a_nDim, *a_Time, *a_nPaths, a_Locs);

  return error;

} /* xfWriteParticleTimestep_f_*/

/* --------------------------------------------------------------------------- */
/* FUNCTION  xfGetPathNullVal_f_*/
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFGETPATHNULLVAL_F,xfgetpathnullval_f,xfgetpathnullval_f_) XFGETPATHNULLVAL_F
#elif defined FLOWER
#define func_name(XFGETPATHNULLVAL_F,xfgetpathnullval_f,xfgetpathnullval_f_) xfgetpathnullval_f
#else
#define func_name(XFGETPATHNULLVAL_F,xfgetpathnullval_f,xfgetpathnullval_f_) xfgetpathnullval_f_
#endif
XMDF_API xid func_name(XFGETPATHNULLVAL_F,xfgetpathnullval_f,xfgetpathnullval_f_)
              (xid *GroupId, double *NullVal)
{
  int error;

  error = xfGetPathNullVal(*GroupId, NullVal);

  return error;
} /* xfGetPathNullVal_f_*/

/* --------------------------------------------------------------------------- */
/* FUNCTION  xfGetNumberOfPaths_f_*/
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFGETNUMBEROFPATHS_F,xfgetnumberofpaths_f,xfgetnumberofpaths_f_) XFGETNUMBEROFPATHS_F
#elif defined FLOWER
#define func_name(XFGETNUMBEROFPATHS_F,xfgetnumberofpaths_f,xfgetnumberofpaths_f_) xfgetnumberofpaths_f
#else
#define func_name(XFGETNUMBEROFPATHS_F,xfgetnumberofpaths_f,xfgetnumberofpaths_f_) xfgetnumberofpaths_f_
#endif
XMDF_API xid func_name(XFGETNUMBEROFPATHS_F,xfgetnumberofpaths_f,xfgetnumberofpaths_f_)
              (xid *GroupId, int *NumPaths)
{
  int error;

  error = xfGetNumberOfPaths(*GroupId, NumPaths);

  return error;

} /* xfGetNumberOfPaths_f_*/

/* --------------------------------------------------------------------------- */
/* FUNCTION  xfGetNumberOfTimes_f_*/
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFGETNUMBEROFTIMES_F,xfgetnumberoftimes_f,xfgetnumberoftimes_f_) XFGETNUMBEROFTIMES_F
#elif defined FLOWER
#define func_name(XFGETNUMBEROFTIMES_F,xfgetnumberoftimes_f,xfgetnumberoftimes_f_) xfgetnumberoftimes_f
#else
#define func_name(XFGETNUMBEROFTIMES_F,xfgetnumberoftimes_f,xfgetnumberoftimes_f_) xfgetnumberoftimes_f_
#endif
XMDF_API xid func_name(XFGETNUMBEROFTIMES_F,xfgetnumberoftimes_f,xfgetnumberoftimes_f_)
              (xid *GroupId, int *NumTimes)
{
  int error;

  error = xfGetNumberOfTimes(*GroupId, NumTimes);

  return error;

} /* xfGetNumberOfTimes_f_*/

/* --------------------------------------------------------------------------- */
/* FUNCTION  xfGetPathDimensionality_f_*/
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFGETPATHDIMENSIONALITY_F,xfgetpathdimensionality_f,xfgetpathdimensionality_f_) XFGETPATHDIMENSIONALITY_F
#elif defined FLOWER
#define func_name(XFGETPATHDIMENSIONALITY_F,xfgetpathdimensionality_f,xfgetpathdimensionality_f_) xfgetpathdimensionality_f
#else
#define func_name(XFGETPATHDIMENSIONALITY_F,xfgetpathdimensionality_f,xfgetpathdimensionality_f_) xfgetpathdimensionality_f_
#endif
XMDF_API xid func_name(XFGETPATHDIMENSIONALITY_F,xfgetpathdimensionality_f,xfgetpathdimensionality_f_)
              (xid *GroupId, int *NumDims)
{
  int error;

  error = xfGetPathDimensionality(*GroupId, NumDims);

  return error;

} /* xfGetPathDimensionality_f_*/

/* --------------------------------------------------------------------------- */
/* FUNCTION  xfGetPathTimesArray_f_*/
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFGETPATHTIMESARRAY_F,xfgetpathtimesarray_f,xfgetpathtimesarray_f_) XFGETPATHTIMESARRAY_F
#elif defined FLOWER
#define func_name(XFGETPATHTIMESARRAY_F,xfgetpathtimesarray_f,xfgetpathtimesarray_f_) xfgetpathtimesarray_f
#else
#define func_name(XFGETPATHTIMESARRAY_F,xfgetpathtimesarray_f,xfgetpathtimesarray_f_) xfgetpathtimesarray_f_
#endif
XMDF_API xid func_name(XFGETPATHTIMESARRAY_F,xfgetpathtimesarray_f,xfgetpathtimesarray_f_)
              (xid *GroupId, int *NumTimes, double *Times)
{
  int error;

  error = xfGetPathTimesArray(*GroupId, *NumTimes, Times);

  return error;

} /* xfGetPathTimesArray_f_*/

/* --------------------------------------------------------------------------- */
/* FUNCTION  xfReadPathLocationsAtTime_f_*/
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFREADPATHLOCATIONSATTIME_F,xfreadpathlocationsattime_f,xfreadpathlocationsattime_f_) XFREADPATHLOCATIONSATTIME_F
#elif defined FLOWER
#define func_name(XFREADPATHLOCATIONSATTIME_F,xfreadpathlocationsattime_f,xfreadpathlocationsattime_f_) xfreadpathlocationsattime_f
#else
#define func_name(XFREADPATHLOCATIONSATTIME_F,xfreadpathlocationsattime_f,xfreadpathlocationsattime_f_) xfreadpathlocationsattime_f_
#endif
XMDF_API xid func_name(XFREADPATHLOCATIONSATTIME_F,xfreadpathlocationsattime_f,xfreadpathlocationsattime_f_)
              (xid *GroupId, int *TimeIndex,
                                  int *FirstPathIndex,
                                  int *NumIndicies, double *Locs)
{
  int error;

  error = xfReadPathLocationsAtTime(*GroupId, *TimeIndex, *FirstPathIndex, 
                                    *NumIndicies, Locs);

  return error;

} /* xfReadPathLocationsAtTime_f_*/

/* --------------------------------------------------------------------------- */
/* FUNCTION  xfReadPathLocsForParticle_f_*/
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFREADPATHLOCSFORPARTICLE_F,xfreadpathlocsforparticle_f,xfreadpathlocsforparticle_f_) XFREADPATHLOCSFORPARTICLE_F
#elif defined FLOWER
#define func_name(XFREADPATHLOCSFORPARTICLE_F,xfreadpathlocsforparticle_f,xfreadpathlocsforparticle_f_) xfreadpathlocsforparticle_f
#else
#define func_name(XFREADPATHLOCSFORPARTICLE_F,xfreadpathlocsforparticle_f,xfreadpathlocsforparticle_f_) xfreadpathlocsforparticle_f_
#endif
XMDF_API xid func_name(XFREADPATHLOCSFORPARTICLE_F,xfreadpathlocsforparticle_f,xfreadpathlocsforparticle_f_)
              (xid *GroupId, int *PathIndex, 
                                       int *FirstTimeIndex, int *NumTimes,
                                       double *Locs)
{
  int error;

  error = xfReadPathLocationsForParticle(*GroupId, *PathIndex, *FirstTimeIndex,
                                         *NumTimes, Locs);

  return error;

} /* xfReadPathLocsForParticle_f_*/

/* --------------------------------------------------------------------------- */
/* FUNCTION  xfReadPathLocsForParticles_f_*/
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFREADPATHLOCSFORPARTICLES_F,xfreadpathlocsforparticles_f,xfreadpathlocsforparticles_f_) XFREADPATHLOCSFORPARTICLES_F
#elif defined FLOWER
#define func_name(XFREADPATHLOCSFORPARTICLES_F,xfreadpathlocsforparticles_f,xfreadpathlocsforparticles_f_) xfreadpathlocsforparticles_f
#else
#define func_name(XFREADPATHLOCSFORPARTICLES_F,xfreadpathlocsforparticles_f,xfreadpathlocsforparticles_f_) xfreadpathlocsforparticles_f_
#endif
XMDF_API xid func_name(XFREADPATHLOCSFORPARTICLES_F,xfreadpathlocsforparticles_f,xfreadpathlocsforparticles_f_)
              (xid *GroupId, int *NumPaths,
                                        const int *PathIndices,
                                        int *FirstTimeIndex, int *NumTimes,
                                        double *Locs)
{
  int error;

  error = xfReadPathLocationsForParticles(*GroupId, *NumPaths, PathIndices, 
                                          *FirstTimeIndex, *NumTimes, Locs);

  return error;

} /* xfReadPathLocsForParticles_f_*/

/* --------------------------------------------------------------------------- */
/* FUNCTION  xfSetupToWriteDatasets_f_*/
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFSETUPTOWRITEDATASETS_F,xfsetuptowritedatasets_f,xfsetuptowritedatasets_f_) XFSETUPTOWRITEDATASETS_F
#elif defined FLOWER
#define func_name(XFSETUPTOWRITEDATASETS_F,xfsetuptowritedatasets_f,xfsetuptowritedatasets_f_) xfsetuptowritedatasets_f
#else
#define func_name(XFSETUPTOWRITEDATASETS_F,xfsetuptowritedatasets_f,xfsetuptowritedatasets_f_) xfsetuptowritedatasets_f_
#endif
XMDF_API xid func_name(XFSETUPTOWRITEDATASETS_F,xfsetuptowritedatasets_f,xfsetuptowritedatasets_f_)
              (const char *a_Filename, int *a_f_ilelen,
                               const char *a_MultiDatasetsGroupPath, int * a_multilen,
                               const char *a_PathInMultiDatasetsGroup, int * a_pathlen, 
                               const char *a_SpatialDataObjectGuid, int *a_spatlen,
                               int *a_OverwriteOptions, hid_t *a_FileId, hid_t *a_GroupId)
{
  int	error;
  char *file_copy = NULL, *multi_copy = NULL, *path_copy = NULL, *spat_copy = NULL;

  file_copy = (char*)malloc((*a_f_ilelen+1)*sizeof(char));
  multi_copy = (char*)malloc((*a_multilen+1)*sizeof(char));
  path_copy = (char*)malloc((*a_pathlen+1)*sizeof(char));
  spat_copy = (char*)malloc((*a_spatlen+1)*sizeof(char));
  strncpy(file_copy, a_Filename, *a_f_ilelen);
  strncpy(multi_copy, a_MultiDatasetsGroupPath, *a_multilen);
  strncpy(path_copy, a_PathInMultiDatasetsGroup, *a_pathlen);
  strncpy(spat_copy, a_SpatialDataObjectGuid, *a_spatlen);
  file_copy[*a_f_ilelen] = '\0';
  multi_copy[*a_multilen] = '\0';
  path_copy[*a_pathlen] = '\0';
  spat_copy[*a_spatlen] = '\0';

  error = xfSetupToWriteDatasets(file_copy, multi_copy, path_copy, spat_copy,
                                 *a_OverwriteOptions, a_FileId, a_GroupId);

  free(file_copy);
  free(multi_copy);
  free(path_copy);
  free(spat_copy);
  return error;

} /* xfSetupToWriteDatasets_f_*/
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfCreateMultiDatasetsGroup_f_*/
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFCREATEMULTIDATASETSGROUP_F,xfcreatemultidatasetsgroup_f,xfcreatemultidatasetsgroup_f_) XFCREATEMULTIDATASETSGROUP_F
#elif defined FLOWER
#define func_name(XFCREATEMULTIDATASETSGROUP_F,xfcreatemultidatasetsgroup_f,xfcreatemultidatasetsgroup_f_) xfcreatemultidatasetsgroup_f
#else
#define func_name(XFCREATEMULTIDATASETSGROUP_F,xfcreatemultidatasetsgroup_f,xfcreatemultidatasetsgroup_f_) xfcreatemultidatasetsgroup_f_
#endif
XMDF_API xid func_name(XFCREATEMULTIDATASETSGROUP_F,xfcreatemultidatasetsgroup_f,xfcreatemultidatasetsgroup_f_)
              (xid *a_Id, const char *a_Path, int *a_pathlen,
               const char *a_Guid, int *a_guidlen, xid *a_MultiId)
{
  int error;
  char *path_copy = NULL, *guid_copy = NULL;

  path_copy = (char*)malloc((*a_pathlen+1)*sizeof(char));
  guid_copy = (char*)malloc((*a_guidlen+1)*sizeof(char));
  strncpy(path_copy, a_Path, *a_pathlen);
  strncpy(guid_copy, a_Guid, *a_guidlen);
  path_copy[*a_pathlen] = '\0';
  guid_copy[*a_guidlen] = '\0';

  error = xfCreateMultiDatasetsGroup(*a_Id, path_copy, guid_copy, a_MultiId);

  free(path_copy);
  free(guid_copy);
  return error;

} /* xfCreateMultiDatasetsGroup_f_*/
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfGetGrpPathsSizeForMltDsets_f_*/
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFGETGRPPATHSSIZEFORMLTDSETS_F,xfgetgrppathssizeformltdsets_f,xfgetgrppathssizeformltdsets_f_) XFGETGRPPATHSSIZEFORMLTDSETS_F
#elif defined FLOWER
#define func_name(XFGETGRPPATHSSIZEFORMLTDSETS_F,xfgetgrppathssizeformltdsets_f,xfgetgrppathssizeformltdsets_f_) xfgetgrppathssizeformltdsets_f
#else
#define func_name(XFGETGRPPATHSSIZEFORMLTDSETS_F,xfgetgrppathssizeformltdsets_f,xfgetgrppathssizeformltdsets_f_) xfgetgrppathssizeformltdsets_f_
#endif
XMDF_API xid func_name(XFGETGRPPATHSSIZEFORMLTDSETS_F,xfgetgrppathssizeformltdsets_f,xfgetgrppathssizeformltdsets_f_)
              (xid *a_Id, int *Num, int *Maxsize)
{
  int error;

  error = xfGetGroupPathsSizeForMultiDatasets(*a_Id, Num, Maxsize);

  return error;

} /* xfGetGrpPathsSizeForMltDsets_f_*/
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfGetAllGrpPathsForMltDsets_f_*/
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFGETALLGRPPATHSFORMLTDSETS_F,xfgetallgrppathsformltdsets_f,xfgetallgrppathsformltdsets_f_) XFGETALLGRPPATHSFORMLTDSETS_F
#elif defined FLOWER
#define func_name(XFGETALLGRPPATHSFORMLTDSETS_F,xfgetallgrppathsformltdsets_f,xfgetallgrppathsformltdsets_f_) xfgetallgrppathsformltdsets_f
#else
#define func_name(XFGETALLGRPPATHSFORMLTDSETS_F,xfgetallgrppathsformltdsets_f,xfgetallgrppathsformltdsets_f_) xfgetallgrppathsformltdsets_f_
#endif
XMDF_API xid func_name(XFGETALLGRPPATHSFORMLTDSETS_F,xfgetallgrppathsformltdsets_f,xfgetallgrppathsformltdsets_f_)
              (xid *a_Id, int *a_Num, int *a_Maxsize, char *a_Paths, 
               int *a_pathlen)
{
  int	error;
  error = xfGetAllGroupPathsForMultiDatasets(*a_Id, *a_Num, *a_Maxsize, a_Paths);
  return error;

} /* xfGetAllGrpPathsForMltDsets_f_*/
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfGetDatasetsSdoGuid_f_*/
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFGETDATASETSSDOGUID_F,xfgetdatasetssdoguid_f,xfgetdatasetssdoguid_f_) XFGETDATASETSSDOGUID_F
#elif defined FLOWER
#define func_name(XFGETDATASETSSDOGUID_F,xfgetdatasetssdoguid_f,xfgetdatasetssdoguid_f_) xfgetdatasetssdoguid_f
#else
#define func_name(XFGETDATASETSSDOGUID_F,xfgetdatasetssdoguid_f,xfgetdatasetssdoguid_f_) xfgetdatasetssdoguid_f_
#endif
XMDF_API xid func_name(XFGETDATASETSSDOGUID_F,xfgetdatasetssdoguid_f,xfgetdatasetssdoguid_f_)
              (xid *a_MultiDatasetsGroup, char *a_GUID, int *a_guidlen)
{
  int	error;
  error = xfGetDatasetsSdoGuid(*a_MultiDatasetsGroup, a_GUID);
  return error;

} /* xfGetDatasetsSdoGuid_f_*/
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfOpenMultiDatasetsGroup_f_*/
/* PURPOSE   Open or create a multi-datasets group inside a mesh or grid */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFOPENMULTIDATASETSGROUP_F,xfopenmultidatasetsgroup_f,xfopenmultidatasetsgroup_f_) XFOPENMULTIDATASETSGROUP_F
#elif defined FLOWER
#define func_name(XFOPENMULTIDATASETSGROUP_F,xfopenmultidatasetsgroup_f,xfopenmultidatasetsgroup_f_) xfopenmultidatasetsgroup_f
#else
#define func_name(XFOPENMULTIDATASETSGROUP_F,xfopenmultidatasetsgroup_f,xfopenmultidatasetsgroup_f_) xfopenmultidatasetsgroup_f_
#endif
XMDF_API xid func_name(XFOPENMULTIDATASETSGROUP_F,xfopenmultidatasetsgroup_f,xfopenmultidatasetsgroup_f_)
              (xid *a_Id, xid *DatasetsGroupId)
{
  int error;

  error = xfOpenMultiDatasetsGroup(*a_Id, DatasetsGroupId);

  return error;

} /* xfOpenMultiDatasetsGroup_f_*/
 /*-----------------------------------------------------------------------------*/
/* FUNCTION  xfCreateScalarDataset_f_*/
/* PURPOSE   Create a Scalar dataset using the Xmdf library.*/
/* NOTES     Calls C version*/
/*-----------------------------------------------------------------------------*/
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFCREATESCALARDATASET_F,xfcreatescalardataset_f,xfcreatescalardataset_f_) XFCREATESCALARDATASET_F
#elif defined FLOWER
#define func_name(XFCREATESCALARDATASET_F,xfcreatescalardataset_f,xfcreatescalardataset_f_) xfcreatescalardataset_f
#else
#define func_name(XFCREATESCALARDATASET_F,xfcreatescalardataset_f,xfcreatescalardataset_f_) xfcreatescalardataset_f_
#endif
XMDF_API xid func_name(XFCREATESCALARDATASET_F,xfcreatescalardataset_f,xfcreatescalardataset_f_)
              (xid *a_DatasetsGroupId, const char *a_Path,
                              int *a_namelen1, const char *a_Units, int *a_namelen2,
                              const char *a_TimeUnits, int *a_namelen3,
                              int *a_Compression, xid *a_DatasetId)
{


  int	compression;
  xid	DatasetsGroupId;
  int	error;
  char *path_copy = NULL, *unit_copy = NULL, *time_copy = NULL;

  path_copy = (char*)malloc((*a_namelen1+1)*sizeof(char));
  unit_copy = (char*)malloc((*a_namelen2+1)*sizeof(char));
  time_copy = (char*)malloc((*a_namelen3+1)*sizeof(char));
  strncpy(path_copy, a_Path, *a_namelen1);
  strncpy(unit_copy, a_Units, *a_namelen2);
  strncpy(time_copy, a_TimeUnits, *a_namelen3);
  path_copy[*a_namelen1] = '\0';
  unit_copy[*a_namelen2] = '\0';
  time_copy[*a_namelen3] = '\0';

  compression = *a_Compression;
  DatasetsGroupId = *a_DatasetsGroupId;
  
  error=xfCreateScalarDataset(DatasetsGroupId, path_copy, unit_copy, time_copy,
                              compression, a_DatasetId);

  free(path_copy);
  free(unit_copy);
  free(time_copy);
  return error;
} /* xfCreateScalarDataset_f_*/
/*-----------------------------------------------------------------------------*/
/* FUNCTION  xfVectorScalarDataset_f_*/
/* PURPOSE   Create a Vector dataset using the Xmdf library.*/
/* NOTES     Calls C version*/
/*-----------------------------------------------------------------------------*/
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFCREATEVECTORDATASET_F,xfcreatevectordataset_f,xfcreatevectordataset_f_) XFCREATEVECTORDATASET_F
#elif defined FLOWER
#define func_name(XFCREATEVECTORDATASET_F,xfcreatevectordataset_f,xfcreatevectordataset_f_) xfcreatevectordataset_f
#else
#define func_name(XFCREATEVECTORDATASET_F,xfcreatevectordataset_f,xfcreatevectordataset_f_) xfcreatevectordataset_f_
#endif
XMDF_API xid func_name(XFCREATEVECTORDATASET_F,xfcreatevectordataset_f,xfcreatevectordataset_f_)
              (xid *a_DatasetsGroupId, const char *a_Path,
                              int *a_namelen1, const char *a_Units, int *a_namelen2,
                              const char *a_TimeUnits, int *a_namelen3,
                              int *a_Compression, xid *a_DatasetId)
{


  int	compression;
  xid	DatasetsGroupId;
  int	error;
  char *path_copy = NULL, *unit_copy = NULL, *time_copy = NULL;

  path_copy = (char*)malloc((*a_namelen1+1)*sizeof(char));
  unit_copy = (char*)malloc((*a_namelen2+1)*sizeof(char));
  time_copy = (char*)malloc((*a_namelen3+1)*sizeof(char));
  strncpy(path_copy, a_Path, *a_namelen1);
  strncpy(unit_copy, a_Units, *a_namelen2);
  strncpy(time_copy, a_TimeUnits, *a_namelen3);
  path_copy[*a_namelen1] = '\0';
  unit_copy[*a_namelen2] = '\0';
  time_copy[*a_namelen3] = '\0';

  compression = *a_Compression;
  DatasetsGroupId = *a_DatasetsGroupId;
  
  error=xfCreateVectorDataset(DatasetsGroupId, path_copy, unit_copy, time_copy,
                              compression, a_DatasetId);

  free(path_copy);
  free(unit_copy);
  free(time_copy);
  return error;
} /* xfCreateVectorDataset_f_*/
/******************************************************************************
 * FUNCTION  xfCreateScalarDsetExtndbl_f_
 * PURPOSE   Create a scalar dataset that can be extended
 * NOTES     The intermediate groups in the path may or may not be created
 ******************2***********************************************************/
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFCREATESCALARDSETEXTNDBL_F,xfcreatescalardsetextndbl_f,xfcreatescalardsetextndbl_f_) XFCREATESCALARDSETEXTNDBL_F
#elif defined FLOWER
#define func_name(XFCREATESCALARDSETEXTNDBL_F,xfcreatescalardsetextndbl_f,xfcreatescalardsetextndbl_f_) xfcreatescalardsetextndbl_f
#else
#define func_name(XFCREATESCALARDSETEXTNDBL_F,xfcreatescalardsetextndbl_f,xfcreatescalardsetextndbl_f_) xfcreatescalardsetextndbl_f_
#endif
XMDF_API xid func_name(XFCREATESCALARDSETEXTNDBL_F,xfcreatescalardsetextndbl_f,xfcreatescalardsetextndbl_f_)
              (xid *a_DatasetsGroupId, const char *a_Path,
                                        int *a_pathlen, const char *a_Units, 
                                        int *a_unitlen, const char *a_TimeUnits, 
                                        int *a_timelen, float *a_FillVal, int *a_Compression, 
                                        xid *a_DatasetId)
{
  int	error;
  char *path_copy = NULL, *unit_copy = NULL, *time_copy = NULL;

  path_copy = (char*)malloc((*a_pathlen+1)*sizeof(char));
  unit_copy = (char*)malloc((*a_unitlen+1)*sizeof(char));
  time_copy = (char*)malloc((*a_timelen+1)*sizeof(char));
  strncpy(path_copy, a_Path, *a_pathlen);
  strncpy(unit_copy, a_Units, *a_unitlen);
  strncpy(time_copy, a_TimeUnits, *a_timelen);
  path_copy[*a_pathlen] = '\0';
  unit_copy[*a_unitlen] = '\0';
  time_copy[*a_timelen] = '\0';

  error = xfCreateScalarDatasetExtendable(*a_DatasetsGroupId, path_copy, unit_copy, time_copy, 
                                          *a_FillVal, *a_Compression, a_DatasetId);

  free(path_copy);
  free(unit_copy);
  free(time_copy);
  return error;

} /* xfCreateScalarDsetExtndbl_f_*/

/******************************************************************************
 * FUNCTION  xfExtendScalarDataset_f_
 * PURPOSE   Create a scalar dataset that can be extended
 * NOTES     The intermediate groups in the path may or may not be created
 ******************2***********************************************************/
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFEXTENDSCALARDATASET_F,xfextendscalardataset_f,xfextendscalardataset_f_) XFEXTENDSCALARDATASET_F
#elif defined FLOWER
#define func_name(XFEXTENDSCALARDATASET_F,xfextendscalardataset_f,xfextendscalardataset_f_) xfextendscalardataset_f
#else
#define func_name(XFEXTENDSCALARDATASET_F,xfextendscalardataset_f,xfextendscalardataset_f_) xfextendscalardataset_f_
#endif
XMDF_API xid func_name(XFEXTENDSCALARDATASET_F,xfextendscalardataset_f,xfextendscalardataset_f_)
              (xid *a_Id, int *aNewSize)
{
  int error;

  error = xfExtendScalarDataset(*a_Id, *aNewSize);

  return error;

} /* xfExtendScalarDataset_f_*/
 /*-----------------------------------------------------------------------------*/
/* FUNCTION  xfWriteScalarTimestep_f_*/
/* PURPOSE   Create a Vector dataset using the Xmdf library.*/
/* NOTES     Calls C version*/
/*-----------------------------------------------------------------------------*/
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFWRITESCALARTIMESTEP_F,xfwritescalartimestep_f,xfwritescalartimestep_f_) XFWRITESCALARTIMESTEP_F
#elif defined FLOWER
#define func_name(XFWRITESCALARTIMESTEP_F,xfwritescalartimestep_f,xfwritescalartimestep_f_) xfwritescalartimestep_f
#else
#define func_name(XFWRITESCALARTIMESTEP_F,xfwritescalartimestep_f,xfwritescalartimestep_f_) xfwritescalartimestep_f_
#endif
XMDF_API xid func_name(XFWRITESCALARTIMESTEP_F,xfwritescalartimestep_f,xfwritescalartimestep_f_)
              (xid *a_Id, double *a_Time, int *a_NumValues,
                              const float *a_Values)
{


  xid    Id;
  double Time;
  int    NumValues;
  int    error;

  Id = *a_Id;
  Time = *a_Time;
  NumValues = *a_NumValues;
  
  error=xfWriteScalarTimestep(Id, Time, NumValues, a_Values);
  return error;
} /* xfWriteScalarTimestep_f_*/

/******************************************************************************
 * FUNCTION  xfWriteScalarTimestepMinMax_f_
 * PURPOSE   
 * NOTES     Append a timestep to a scalar dataset
 *           The min and max should be set to the min or max of the timestep
 ******************2***********************************************************/
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFWRITESCALARTIMESTEPMINMAX_F,xfwritescalartimestepminmax_f,xfwritescalartimestepminmax_f_) XFWRITESCALARTIMESTEPMINMAX_F
#elif defined FLOWER
#define func_name(XFWRITESCALARTIMESTEPMINMAX_F,xfwritescalartimestepminmax_f,xfwritescalartimestepminmax_f_) xfwritescalartimestepminmax_f
#else
#define func_name(XFWRITESCALARTIMESTEPMINMAX_F,xfwritescalartimestepminmax_f,xfwritescalartimestepminmax_f_) xfwritescalartimestepminmax_f_
#endif
XMDF_API xid func_name(XFWRITESCALARTIMESTEPMINMAX_F,xfwritescalartimestepminmax_f,xfwritescalartimestepminmax_f_)
              (xid *a_Id, double *a_Time, int *a_NumValues, 
                                    const float *a_Values, float *a_Min, 
                                    float *a_Max)
{
  int error;

  error = xfWriteScalarTimestepMinMax(*a_Id, *a_Time, *a_NumValues, a_Values, *a_Min, 
                                      *a_Max);

  return error;

} /* xfWriteScalarTimestepMinMax_f_*/
 /*-----------------------------------------------------------------------------*/
/* FUNCTION  xfInitializeScalarTimestep_f */
/*-----------------------------------------------------------------------------*/
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFINITIALIZESCALARTIMESTEP_F,xfinitializescalartimestep_f,xfinitializescalartimestep_f_) XFINITIALIZESCALARTIMESTEP_F
#elif defined FLOWER
#define func_name(XFINITIALIZESCALARTIMESTEP_F,xfinitializescalartimestep_f,xfinitializescalartimestep_f_) xfinitializescalartimestep_f
#else
#define func_name(XFINITIALIZESCALARTIMESTEP_F,xfinitializescalartimestep_f,xfinitializescalartimestep_f_) xfinitializescalartimestep_f_
#endif
XMDF_API xid func_name(XFINITIALIZESCALARTIMESTEP_F,xfinitializescalartimestep_f,xfinitializescalartimestep_f_)
              (xid *a_Id, double *a_Time, int *a_NumValues,
               const float *minvalue, const float *maxvalue,
               int *a_timestepId)
{
  xid    Id;
  double Time;
  int    NumValues;
  int    error;
  hsize_t timestep;

  Id = *a_Id;
  Time = *a_Time;
  NumValues = *a_NumValues;
  timestep = *a_timestepId;
  
  error=xfInitializeScalarTimestep(Id, Time, NumValues, *minvalue, *maxvalue,
                                   &timestep);
  *a_timestepId = (int)timestep;

  return error;
} /* xfInitializeScalarTimestep_f*/
 /*-----------------------------------------------------------------------------*/
/* FUNCTION  xfWriteScalarTimestepPortion_f*/
/* PURPOSE   Create a Vector dataset using the Xmdf library.*/
/* NOTES     Calls C version*/
/*-----------------------------------------------------------------------------*/
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFWRITESCALARTIMESTEPPORTION_F,xfwritescalartimestepportion_f,xfwritescalartimestepportion_f_) XFWRITESCALARTIMESTEPPORTION_F
#elif defined FLOWER
#define func_name(XFWRITESCALARTIMESTEPPORTION_F,xfwritescalartimestepportion_f,xfwritescalartimestepportion_f_) xfwritescalartimestepportion_f
#else
#define func_name(XFWRITESCALARTIMESTEPPORTION_F,xfwritescalartimestepportion_f,xfwritescalartimestepportion_f_) xfwritescalartimestepportion_f_
#endif
XMDF_API xid func_name(XFWRITESCALARTIMESTEPPORTION_F,xfwritescalartimestepportion_f,xfwritescalartimestepportion_f_)
              (xid *a_Id, int *a_timestepId, int *a_NumValuesToWrite,
               int *a_startIndex, const float *a_Values)
{


  xid    Id;
  int    timestepId;
  int    NumValuesToWrite;
  int    StartIndex;
  int    error;

  Id = *a_Id;
  timestepId = *a_timestepId;
  NumValuesToWrite = *a_NumValuesToWrite;
  StartIndex = *a_startIndex;
  
  error=xfWriteScalarTimestepPortion(Id, timestepId, NumValuesToWrite,
          StartIndex, a_Values);
  return error;
} /* xfWriteScalarTimestepPortion_f*/
/*-----------------------------------------------------------------------------*/
/* FUNCTION  xfSetDatasetTimestepMinMax_f*/
/* PURPOSE   Create a Vector dataset using the Xmdf library.*/
/* NOTES     Calls C version*/
/*-----------------------------------------------------------------------------*/
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFSETDATASETTIMESTEPMINMAX_F,xfsetdatasettimestepminmax_f,xfsetdatasettimestepminmax_f_) XFSETDATASETTIMESTEPMINMAX_F
#elif defined FLOWER
#define func_name(XFSETDATASETTIMESTEPMINMAX_F,xfsetdatasettimestepminmax_f,xfsetdatasettimestepminmax_f_) xfsetdatasettimestepminmax_f
#else
#define func_name(XFSETDATASETTIMESTEPMINMAX_F,xfsetdatasettimestepminmax_f,xfsetdatasettimestepminmax_f_) xfsetdatasettimestepminmax_f_
#endif
XMDF_API xid func_name(XFSETDATASETTIMESTEPMINMAX_F,xfsetdatasettimestepminmax_f,xfsetdatasettimestepminmax_f_)
              (xid *a_Id, int *a_timestepId, const float *a_minvalue,
               const float *a_maxvalue)
{


  xid    Id;
  int    timestepId;
  int    error;
  float  minvalue, maxvalue;

  Id = *a_Id;
  timestepId = *a_timestepId;
  minvalue = *a_minvalue;
  maxvalue = *a_maxvalue;
  
  error = xfSetDatasetTimestepMinMax(Id, timestepId, minvalue, maxvalue);

  return error;
} /* xfSetDatasetTimestepMinMax_f*/

/******************************************************************************
 * FUNCTION  xfWriteVectorTimestep_f_
 * PURPOSE   
 * NOTES     Append a timestep to a vector dataset
 ******************2***********************************************************/
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFWRITEVECTORTIMESTEP_F,xfwritevectortimestep_f,xfwritevectortimestep_f_) XFWRITEVECTORTIMESTEP_F
#elif defined FLOWER
#define func_name(XFWRITEVECTORTIMESTEP_F,xfwritevectortimestep_f,xfwritevectortimestep_f_) xfwritevectortimestep_f
#else
#define func_name(XFWRITEVECTORTIMESTEP_F,xfwritevectortimestep_f,xfwritevectortimestep_f_) xfwritevectortimestep_f_
#endif
XMDF_API xid func_name(XFWRITEVECTORTIMESTEP_F,xfwritevectortimestep_f,xfwritevectortimestep_f_)
              (xid *a_Id, double *a_Time, int *a_NumValues, 
                              int *a_NumComponents, const float *a_Values)
{
  int error;

  error = xfWriteVectorTimestep(*a_Id, *a_Time, *a_NumValues, *a_NumComponents, a_Values);

  return error;

} /* xfWriteVectorTimestep_f_*/
/******************************************************************************
 * FUNCTION  xfWriteVectorTimestepMinMax_f_
 * PURPOSE   
 * NOTES     Append a timestep to a vector dataset
 *           The min and max should be set to the min or max (magnitude)
 *           of the timestep
 ******************2***********************************************************/
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFWRITEVECTORTIMESTEPMINMAX_F,xfwritevectortimestepminmax_f,xfwritevectortimestepminmax_f_) XFWRITEVECTORTIMESTEPMINMAX_F
#elif defined FLOWER
#define func_name(XFWRITEVECTORTIMESTEPMINMAX_F,xfwritevectortimestepminmax_f,xfwritevectortimestepminmax_f_) xfwritevectortimestepminmax_f
#else
#define func_name(XFWRITEVECTORTIMESTEPMINMAX_F,xfwritevectortimestepminmax_f,xfwritevectortimestepminmax_f_) xfwritevectortimestepminmax_f_
#endif
XMDF_API xid func_name(XFWRITEVECTORTIMESTEPMINMAX_F,xfwritevectortimestepminmax_f,xfwritevectortimestepminmax_f_)
              (xid *a_Id, double *a_Time, int *a_NumValues, 
                                    int *a_NumComponents, const float *a_Values,
                                    float *a_Min, float *a_Max)
{
  int error;

  error = xfWriteVectorTimestepMinMax(*a_Id, *a_Time, *a_NumValues, *a_NumComponents, 
                                      a_Values, *a_Min, *a_Max);

  return error;

} /* xfWriteVectorTimestepMinMax_f_*/
/******************************************************************************
 * FUNCTION  xfInitializeVectorTimestep_f
 * PURPOSE   
 * NOTES     Append a timestep to a vector dataset
 *           The min and max should be set to the min or max (magnitude)
 *           of the timestep
 ******************2***********************************************************/
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFINITIALIZEVECTORTIMESTEP_F,xfinitializevectortimestep_f,xfinitializevectortimestep_f_) XFINITIALIZEVECTORTIMESTEP_F
#elif defined FLOWER
#define func_name(XFINITIALIZEVECTORTIMESTEP_F,xfinitializevectortimestep_f,xfinitializevectortimestep_f_) xfinitializevectortimestep_f
#else
#define func_name(XFINITIALIZEVECTORTIMESTEP_F,xfinitializevectortimestep_f,xfinitializevectortimestep_f_) xfinitializevectortimestep_f_
#endif
XMDF_API xid func_name(XFINITIALIZEVECTORTIMESTEP_F,xfinitializevectortimestep_f,xfinitializevectortimestep_f_)
              (xid *a_Id, double *a_Time, int *a_NumValues, 
               int *a_NumComponents, float *a_Min, float *a_Max, int *a_timeId)
{
  int error;

  hsize_t timeId = *a_timeId;

  error = xfInitializeVectorTimestep(*a_Id, *a_Time, *a_NumValues, *a_NumComponents, 
                                     *a_Min, *a_Max, &timeId);

  *a_timeId = (int)timeId;

  return error;
} /* xfInitializeVectorTimestep_f */
/******************************************************************************
 * FUNCTION  xfWriteVectorTimestepPortion_f
 * PURPOSE   
 ******************2***********************************************************/
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFWRITEVECTORTIMESTEPPORTION_F,xfwritevectortimestepportion_f,xfwritevectortimestepportion_f_) XFWRITEVECTORTIMESTEPPORTION_F
#elif defined FLOWER
#define func_name(XFWRITEVECTORTIMESTEPPORTION_F,xfwritevectortimestepportion_f,xfwritevectortimestepportion_f_) xfwritevectortimestepportion_f
#else
#define func_name(XFWRITEVECTORTIMESTEPPORTION_F,xfwritevectortimestepportion_f,xfwritevectortimestepportion_f_) xfwritevectortimestepportion_f_
#endif
XMDF_API xid func_name(XFWRITEVECTORTIMESTEPPORTION_F,xfwritevectortimestepportion_f,xfwritevectortimestepportion_f_)
              (xid *a_Id, int *a_TimeId, int *a_NumValuesToWrite, 
               int *a_NumComponentsToWrite, int *startIndex,
               int *startComponent, const float *a_Values)
{
  int error;

  error = xfWriteVectorTimestepPortion(*a_Id, *a_TimeId, *a_NumValuesToWrite,
           *a_NumComponentsToWrite, *startIndex, *startComponent, a_Values);

  return error;
} /* xfWriteVectorTimestepPortion_f*/

/******************************************************************************
 * FUNCTION  xfWriteActivityTimestep_f_
 * PURPOSE   
 * NOTES     Write dataset activity information
 ******************2***********************************************************/
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFWRITEACTIVITYTIMESTEP_F,xfwriteactivitytimestep_f,xfwriteactivitytimestep_f_) XFWRITEACTIVITYTIMESTEP_F
#elif defined FLOWER
#define func_name(XFWRITEACTIVITYTIMESTEP_F,xfwriteactivitytimestep_f,xfwriteactivitytimestep_f_) xfwriteactivitytimestep_f
#else
#define func_name(XFWRITEACTIVITYTIMESTEP_F,xfwriteactivitytimestep_f,xfwriteactivitytimestep_f_) xfwriteactivitytimestep_f_
#endif
XMDF_API xid func_name(XFWRITEACTIVITYTIMESTEP_F,xfwriteactivitytimestep_f,xfwriteactivitytimestep_f_)
              (xid *a_Id, int *a_NumActive, const unsigned char *a_Active)
{
  int error;

  error = xfWriteActivityTimestep(*a_Id, *a_NumActive, a_Active);

  return error;

} /* xfWriteActivityTimestep_f_*/
 /*-----------------------------------------------------------------------------*/
/* FUNCTION  xfInitializeActivityTimestep_f */
/*-----------------------------------------------------------------------------*/
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFINITIALIZEACTIVITYTIMESTEP_F,xfinitializeactivitytimestep_f,xfinitializeactivitytimestep_f_) XFINITIALIZEACTIVITYTIMESTEP_F
#elif defined FLOWER
#define func_name(XFINITIALIZEACTIVITYTIMESTEP_F,xfinitializeactivitytimestep_f,xfinitializeactivitytimestep_f_) xfinitializeactivitytimestep_f
#else
#define func_name(XFINITIALIZEACTIVITYTIMESTEP_F,xfinitializeactivitytimestep_f,xfinitializeactivitytimestep_f_) xfinitializeactivitytimestep_f_
#endif
XMDF_API xid func_name(XFINITIALIZESCALARTIMESTEP_F,xfinitializeactivitytimestep_f,xfinitializeactivitytimestep_f_)
              (xid *a_Id, int *a_NumActive, xid *a_timestepId)
{
  xid    Id;
  int    NumActive;
  int    error;
  hsize_t timestepId = *a_timestepId;

  Id = *a_Id;
  NumActive = *a_NumActive;
  
  error = xfInitializeActivityTimestep(Id, NumActive, &timestepId);

  *a_timestepId = (xid)timestepId;

  return error;
} /* xfInitializeActivityTimestep_f*/
 /*-----------------------------------------------------------------------------*/
/* FUNCTION  xfWriteActivityTimestepPortion_f*/
/* NOTES     Calls C version*/
/*-----------------------------------------------------------------------------*/
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFWRITEACTIVITYTIMESTEPPORTION_F,xfwriteactivitytimestepportion_f,xfwriteactivitytimestepportion_f_) XFWRITEACTIVITYTIMESTEPPORTION_F
#elif defined FLOWER
#define func_name(XFWRITEACTIVITYTIMESTEPPORTION_F,xfwriteactivitytimestepportion_f,xfwriteactivitytimestepportion_f_) xfwriteactivitytimestepportion_f
#else
#define func_name(XFWRITEACTIVITYTIMESTEPPORTION_F,xfwriteactivitytimestepportion_f,xfwriteactivitytimestepportion_f_) xfwriteactivitytimestepportion_f_
#endif
XMDF_API xid func_name(XFWRITEACTIVITYTIMESTEPPORTION_F,xfwriteactivitytimestepportion_f,xfwriteactivitytimestepportion_f_)
              (xid *a_Id, int *a_timestepId, int *a_NumValuesToWrite,
               int *a_startIndex, const unsigned char *a_ActivityValues)
{


  xid    Id;
  int    timestepId;
  int    NumValuesToWrite;
  int    StartIndex;
  int    error;

  Id = *a_Id;
  timestepId = *a_timestepId;
  NumValuesToWrite = *a_NumValuesToWrite;
  StartIndex = *a_startIndex;
  
  error=xfWriteActivityTimestepPortion(Id, timestepId, NumValuesToWrite,
          StartIndex, a_ActivityValues);

  return error;
} /* xfWriteActivityTimestepPortion_f*/
 /******************************************************************************
 * FUNCTION  xfGetScalarDatasetGroupId_f_
 * PURPOSE   
 * NOTES     Open the Scalar datasets for a Mesh or Grid
 ******************2***********************************************************/
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFGETSCALARDATASETGROUPID_F,xfgetscalardatasetgroupid_f,xfgetscalardatasetgroupid_f_) XFGETSCALARDATASETGROUPID_F
#elif defined FLOWER
#define func_name(XFGETSCALARDATASETGROUPID_F,xfgetscalardatasetgroupid_f,xfgetscalardatasetgroupid_f_) xfgetscalardatasetgroupid_f
#else
#define func_name(XFGETSCALARDATASETGROUPID_F,xfgetscalardatasetgroupid_f,xfgetscalardatasetgroupid_f_) xfgetscalardatasetgroupid_f_
#endif
XMDF_API xid func_name(XFGETSCALARDATASETGROUPID_F,xfgetscalardatasetgroupid_f,xfgetscalardatasetgroupid_f_)
              (xid *a_Id)
{
  int error;

  error = xfGetScalarDatasetGroupId(*a_Id);

  return error;

} /* xfGetScalarDatasetGroupId_f_*/
 /******************************************************************************
 * FUNCTION  xfGetVectorDatasetGroupId_f_
 * PURPOSE   
 * NOTES     Open the Scalar datasets for a Mesh or Grid
 ******************2***********************************************************/
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFGETVECTORDATASETGROUPID_F,xfgetvectordatasetgroupid_f,xfgetvectordatasetgroupid_f_) XFGETVECTORDATASETGROUPID_F
#elif defined FLOWER
#define func_name(XFGETVECTORDATASETGROUPID_F,xfgetvectordatasetgroupid_f,xfgetvectordatasetgroupid_f_) xfgetvectordatasetgroupid_f
#else
#define func_name(XFGETVECTORDATASETGROUPID_F,xfgetvectordatasetgroupid_f,xfgetvectordatasetgroupid_f_) xfgetvectordatasetgroupid_f_
#endif
XMDF_API xid func_name(XFGETVECTORDATASETGROUPID_F,xfgetvectordatasetgroupid_f,xfgetvectordatasetgroupid_f_)
              (xid *a_Id)
{
  int error;

  error = xfGetVectorDatasetGroupId(*a_Id);

  return error;

} /* xfGetVectorDatasetGroupId_f_*/
 /******************************************************************************
 * FUNCTION  xfWriteReftime_f_
 * PURPOSE   Set a reference time for a dataset
 * NOTES     
 ******************2***********************************************************/
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFWRITEREFTIME_F,xfwritereftime_f,xfwritereftime_f_) XFWRITEREFTIME_F
#elif defined FLOWER
#define func_name(XFWRITEREFTIME_F,xfwritereftime_f,xfwritereftime_f_) xfwritereftime_f
#else
#define func_name(XFWRITEREFTIME_F,xfwritereftime_f,xfwritereftime_f_) xfwritereftime_f_
#endif
XMDF_API xid func_name(XFWRITEREFTIME_F,xfwritereftime_f,xfwritereftime_f_)
              (xid *a_Id, double *a_Reftime)
{
  int error;

  error = xfWriteReftime(*a_Id, *a_Reftime);

  return error;

} /* xfWriteReftime_f_*/
/******************************************************************************
 * FUNCTION  xfUseReftime_f_
 * PURPOSE   See if a reference time exists for a dataset
 * NOTES     
 ******************2***********************************************************/
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFUSEREFTIME_F,xfusereftime_f,xfusereftime_f_) XFUSEREFTIME_F
#elif defined FLOWER
#define func_name(XFUSEREFTIME_F,xfusereftime_f,xfusereftime_f_) xfusereftime_f
#else
#define func_name(XFUSEREFTIME_F,xfusereftime_f,xfusereftime_f_) xfusereftime_f_
#endif
XMDF_API xid func_name(XFUSEREFTIME_F,xfusereftime_f,xfusereftime_f_)
              (xid *a_Id, short *a_bUseReftime)
{
  
  xmbool  bUseReftime;
  int     error;

  bUseReftime = (*a_bUseReftime) ? XTRUE : XFALSE;
  error = xfUseReftime(*a_Id, &bUseReftime);

  return error;

} /* xfUseReftime_f_*/
/******************************************************************************
 * FUNCTION  xfReadReftime_f_
 * PURPOSE   Read a reference time for a dataset
 * NOTES     
 ******************2***********************************************************/
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFREADREFTIME_F,xfreadreftime_f,xfreadreftime_f_) XFREADREFTIME_F
#elif defined FLOWER
#define func_name(XFREADREFTIME_F,xfreadreftime_f,xfreadreftime_f_) xfreadreftime_f
#else
#define func_name(XFREADREFTIME_F,xfreadreftime_f,xfreadreftime_f_) xfreadreftime_f_
#endif
XMDF_API xid func_name(XFREADREFTIME_F,xfreadreftime_f,xfreadreftime_f_)
              (xid *a_Id, double *a_dReftime)
{
  int error;

  error = xfReadReftime(*a_Id, a_dReftime);

  return error;

} /* xfReadReftime_f_*/
 /******************************************************************************
 * FUNCTION  xfGetScalarDatasetsInfo_f_
 * PURPOSE   
 * NOTES     Get the number and max path length of scalar datasets in the path
 ******************2***********************************************************/
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFGETSCALARDATASETSINFO_F,xfgetscalardatasetsinfo_f,xfgetscalardatasetsinfo_f_) XFGETSCALARDATASETSINFO_F
#elif defined FLOWER
#define func_name(XFGETSCALARDATASETSINFO_F,xfgetscalardatasetsinfo_f,xfgetscalardatasetsinfo_f_) xfgetscalardatasetsinfo_f
#else
#define func_name(XFGETSCALARDATASETSINFO_F,xfgetscalardatasetsinfo_f,xfgetscalardatasetsinfo_f_) xfgetscalardatasetsinfo_f_
#endif
XMDF_API xid func_name(XFGETSCALARDATASETSINFO_F,xfgetscalardatasetsinfo_f,xfgetscalardatasetsinfo_f_)
              (xid *a_Id, int *a_Number, int *a_MaxPathLength)
{
  int error;

  error = xfGetScalarDatasetsInfo(*a_Id, a_Number, a_MaxPathLength);

  return error;

} /* xfGetScalarDatasetsInfo_f_*/
/******************************************************************************
 * FUNCTION  xfGetScalarDatasetPaths_f_
 * PURPOSE   Get the paths to scalar datasets under a starting group
 * NOTES     The Path array must already be allocated to a size Number by Maxlength
 ******************2***********************************************************/
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFGETSCALARDATASETPATHS_F,xfgetscalardatasetpaths_f,xfgetscalardatasetpaths_f_) XFGETSCALARDATASETPATHS_F
#elif defined FLOWER
#define func_name(XFGETSCALARDATASETPATHS_F,xfgetscalardatasetpaths_f,xfgetscalardatasetpaths_f_) xfgetscalardatasetpaths_f
#else
#define func_name(XFGETSCALARDATASETPATHS_F,xfgetscalardatasetpaths_f,xfgetscalardatasetpaths_f_) xfgetscalardatasetpaths_f_
#endif
XMDF_API xid func_name(XFGETSCALARDATASETPATHS_F,xfgetscalardatasetpaths_f,xfgetscalardatasetpaths_f_)
              (xid *a_Id, int *a_Number, int *a_MaxPathLength, 
                                char *a_Paths, int *a_pathlen)
{
  int error;
  error = xfGetScalarDatasetPaths(*a_Id, *a_Number, *a_MaxPathLength, a_Paths);
  return error;

} /* xfGetScalarDatasetPaths_f_*/
/******************************************************************************
 * FUNCTION  xfGetVectorDatasetsInfo_f_
 * PURPOSE   
 * NOTES     Get the number and max path length of scalar datasets in the path
 ******************2***********************************************************/
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFGETVECTORDATASETSINFO_F,xfgetvectordatasetsinfo_f,xfgetvectordatasetsinfo_f_) XFGETVECTORDATASETSINFO_F
#elif defined FLOWER
#define func_name(XFGETVECTORDATASETSINFO_F,xfgetvectordatasetsinfo_f,xfgetvectordatasetsinfo_f_) xfgetvectordatasetsinfo_f
#else
#define func_name(XFGETVECTORDATASETSINFO_F,xfgetvectordatasetsinfo_f,xfgetvectordatasetsinfo_f_) xfgetvectordatasetsinfo_f_
#endif
XMDF_API xid func_name(XFGETVECTORDATASETSINFO_F,xfgetvectordatasetsinfo_f,xfgetvectordatasetsinfo_f_)
              (xid *a_Id, int *a_Number, int *a_MaxPathLength)
{
  int error;

  error = xfGetVectorDatasetsInfo(*a_Id, a_Number, a_MaxPathLength);

  return error;

} /* xfGetVectorDatasetsInfo_f_*/
/******************************************************************************
 * FUNCTION  xfGetVectorDatasetPaths_f_
 * PURPOSE   Get the paths to scalar datasets under a starting group
 * NOTES     The Path array must already be allocated to a size Number by Maxlength
 ******************2***********************************************************/
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFGETVECTORDATASETPATHS_F,xfgetvectordatasetpaths_f,xfgetvectordatasetpaths_f_) XFGETVECTORDATASETPATHS_F
#elif defined FLOWER
#define func_name(XFGETVECTORDATASETPATHS_F,xfgetvectordatasetpaths_f,xfgetvectordatasetpaths_f_) xfgetvectordatasetpaths_f
#else
#define func_name(XFGETVECTORDATASETPATHS_F,xfgetvectordatasetpaths_f,xfgetvectordatasetpaths_f_) xfgetvectordatasetpaths_f_
#endif
XMDF_API xid func_name(XFGETVECTORDATASETPATHS_F,xfgetvectordatasetpaths_f,xfgetvectordatasetpaths_f_)
              (xid *a_Id, int *a_Number, int *a_MaxPathLength, 
                                char *a_Paths, int *a_pathlen)
{
  int	error;
  error = xfGetVectorDatasetPaths(*a_Id, *a_Number, *a_MaxPathLength, a_Paths);
  return error;

} /* xfGetVectorDatasetPaths_f_*/
 /******************************************************************************
 * FUNCTION  xfReadDatasetReftime_f_
 * PURPOSE   Read a reference time for a dataset
 * NOTES     
 ******************2***********************************************************/
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFREADDATASETREFTIME_F,xfreaddatasetreftime_f,xfreaddatasetreftime_f_) XFREADDATASETREFTIME_F
#elif defined FLOWER
#define func_name(XFREADDATASETREFTIME_F,xfreaddatasetreftime_f,xfreaddatasetreftime_f_) xfreaddatasetreftime_f
#else
#define func_name(XFREADDATASETREFTIME_F,xfreaddatasetreftime_f,xfreaddatasetreftime_f_) xfreaddatasetreftime_f_
#endif
XMDF_API xid func_name(XFREADDATASETREFTIME_F,xfreaddatasetreftime_f,xfreaddatasetreftime_f_)
              (xid *a_Id, double *a_dReftime)
{
  int error;

  error = xfReadDatasetReftime(*a_Id, a_dReftime);

  return error;

} /* xfReadDatasetReftime_f_*/
/******************************************************************************
 * FUNCTION  xfSetDatasetNumTimes_f_
 * PURPOSE   
 * NOTES     
 ******************2***********************************************************/
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFSETDATASETNUMTIMES_F,xfsetdatasetnumtimes_f,xfsetdatasetnumtimes_f_) XFSETDATASETNUMTIMES_F
#elif defined FLOWER
#define func_name(XFSETDATASETNUMTIMES_F,xfsetdatasetnumtimes_f,xfsetdatasetnumtimes_f_) xfsetdatasetnumtimes_f
#else
#define func_name(XFSETDATASETNUMTIMES_F,xfsetdatasetnumtimes_f,xfsetdatasetnumtimes_f_) xfsetdatasetnumtimes_f_
#endif
XMDF_API xid func_name(XFSETDATASETNUMTIMES_F,xfsetdatasetnumtimes_f,xfsetdatasetnumtimes_f_)
              (xid *a_Id, int *a_Numtimes)
{
  int error;

  error = xfSetDatasetNumTimes(*a_Id, *a_Numtimes);

  return error;

} /* xfSetDatasetNumTimes_f_*/
/******************************************************************************
 * FUNCTION  xfChangeScalarValuesTimestepFloat_f
 * PURPOSE   
 * NOTES     
 ******************2***********************************************************/
#ifdef func_name
#undef func_name
#endif

#if defined FUPPER
#define func_name XFCHANGESCALARVALUESTIMESTEPFLOAT_F
#elif defined FLOWER
#define func_name xfchangescalarvaluestimestepfloat_f
#else
#define func_name xfchangescalarvaluestimestepfloat_f_
#endif

XMDF_API xid func_name (xid *a_Id, int *a_TimestepIndex,
                        int *a_NumValsToEdit, int *a_Indices, float *a_NewValues)
{
  int error;

  error = xfChangeScalarValuesTimestepFloat(*a_Id, *a_TimestepIndex, *a_NumValsToEdit, a_Indices, a_NewValues);

  return error;

} /* xfChangeScalarValuesTimestepFloat_f*/
/******************************************************************************
 * FUNCTION  xfGetDatasetNumTimes_f_
 * PURPOSE   
 * NOTES     
 ******************2***********************************************************/
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFGETDATASETNUMTIMES_F,xfgetdatasetnumtimes_f,xfgetdatasetnumtimes_f_) XFGETDATASETNUMTIMES_F
#elif defined FLOWER
#define func_name(XFGETDATASETNUMTIMES_F,xfgetdatasetnumtimes_f,xfgetdatasetnumtimes_f_) xfgetdatasetnumtimes_f
#else
#define func_name(XFGETDATASETNUMTIMES_F,xfgetdatasetnumtimes_f,xfgetdatasetnumtimes_f_) xfgetdatasetnumtimes_f_
#endif
XMDF_API xid func_name(XFGETDATASETNUMTIMES_F,xfgetdatasetnumtimes_f,xfgetdatasetnumtimes_f_)
              (xid *a_Id, int *a_Numtimes)
{
  int error;

  error = xfGetDatasetNumTimes(*a_Id, a_Numtimes);

  return error;

} /* xfGetDatasetNumTimes_f_*/
/******************************************************************************
 * FUNCTION  xfGetDatasetNumVals_f_
 * PURPOSE   
 * NOTES     
 ******************2***********************************************************/
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFGETDATASETNUMVALS_F,xfgetdatasetnumvals_f,xfgetdatasetnumvals_f_) XFGETDATASETNUMVALS_F
#elif defined FLOWER
#define func_name(XFGETDATASETNUMVALS_F,xfgetdatasetnumvals_f,xfgetdatasetnumvals_f_) xfgetdatasetnumvals_f
#else
#define func_name(XFGETDATASETNUMVALS_F,xfgetdatasetnumvals_f,xfgetdatasetnumvals_f_) xfgetdatasetnumvals_f_
#endif
XMDF_API xid func_name(XFGETDATASETNUMVALS_F,xfgetdatasetnumvals_f,xfgetdatasetnumvals_f_)
              (xid *a_Id, int *a_Numvals)
{
  int error;

  error = xfGetDatasetNumVals(*a_Id, a_Numvals);

  return error;

} /* xfGetDatasetNumVals_f_*/
/******************************************************************************
 * FUNCTION  xfGetDatasetNumActive_f_
 * PURPOSE   
 * NOTES     
 ******************2***********************************************************/
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFGETDATASETNUMACTIVE_F,xfgetdatasetnumactive_f,xfgetdatasetnumactive_f_) XFGETDATASETNUMACTIVE_F
#elif defined FLOWER
#define func_name(XFGETDATASETNUMACTIVE_F,xfgetdatasetnumactive_f,xfgetdatasetnumactive_f_) xfgetdatasetnumactive_f
#else
#define func_name(XFGETDATASETNUMACTIVE_F,xfgetdatasetnumactive_f,xfgetdatasetnumactive_f_) xfgetdatasetnumactive_f_
#endif
XMDF_API xid func_name(XFGETDATASETNUMACTIVE_F,xfgetdatasetnumactive_f,xfgetdatasetnumactive_f_)
              (xid *a_Id, int *a_NumActivevals)
{
  int error;

  error = xfGetDatasetNumActive(*a_Id, a_NumActivevals);

  return error;

} /* xfGetDatasetNumActive_f_*/
/******************************************************************************
 * FUNCTION  xfGetDatasetVecNumComponents_f_
 * PURPOSE   
 * NOTES     
 ******************2***********************************************************/
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFGETDATASETVECNUMCOMPONENTS_F,xfgetdatasetvecnumcomponents_f,xfgetdatasetvecnumcomponents_f_) XFGETDATASETVECNUMCOMPONENTS_F
#elif defined FLOWER
#define func_name(XFGETDATASETVECNUMCOMPONENTS_F,xfgetdatasetvecnumcomponents_f,xfgetdatasetvecnumcomponents_f_) xfgetdatasetvecnumcomponents_f
#else
#define func_name(XFGETDATASETVECNUMCOMPONENTS_F,xfgetdatasetvecnumcomponents_f,xfgetdatasetvecnumcomponents_f_) xfgetdatasetvecnumcomponents_f_
#endif
XMDF_API xid func_name(XFGETDATASETVECNUMCOMPONENTS_F,xfgetdatasetvecnumcomponents_f,xfgetdatasetvecnumcomponents_f_)
              (xid *a_Id, int *a_NumComponents)
{
  int error;

  error = xfGetDatasetVecNumComponents(*a_Id, a_NumComponents);

  return error;

} /* xfGetDatasetVecNumComponents_f_*/
/******************************************************************************
 * FUNCTION  xfGetDatasetTimeUnits_f_
 * PURPOSE   
 * NOTES     The units variable should arleady be allocated 
 *           size >= TIME_UNITS_MAXLENGTH
 ******************2***********************************************************/
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFGETDATASETTIMEUNITS_F,xfgetdatasettimeunits_f,xfgetdatasettimeunits_f_) XFGETDATASETTIMEUNITS_F
#elif defined FLOWER
#define func_name(XFGETDATASETTIMEUNITS_F,xfgetdatasettimeunits_f,xfgetdatasettimeunits_f_) xfgetdatasettimeunits_f
#else
#define func_name(XFGETDATASETTIMEUNITS_F,xfgetdatasettimeunits_f,xfgetdatasettimeunits_f_) xfgetdatasettimeunits_f_
#endif
XMDF_API xid func_name(XFGETDATASETTIMEUNITS_F,xfgetdatasettimeunits_f,xfgetdatasettimeunits_f_)
              (xid *a_Id, char *Units, int *a_unitlen)
{
  int	error;
  error = xfGetDatasetTimeUnits(*a_Id, Units);
  return error;

} /* xfGetDatasetTimeUnits_f_*/
/******************************************************************************
 * FUNCTION  xfgetdatasetunits_f_
 * PURPOSE   Get the units under a starting group
 * NOTES     The a_Units array must already be allocated to a size Number by a_MaxUnitsLength
 ******************2***********************************************************/
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFGETDATASETUNITS_F,xfgetdatasetunits_f,xfgetdatasetunits_f_) XFGETDATASETUNITS_F
#elif defined FLOWER
#define func_name(XFGETDATASETUNITS_F,xfgetdatasetunits_f,xfgetdatasetunits_f_) xfgetdatasetunits_f
#else
#define func_name(XFGETDATASETUNITS_F,xfgetdatasetunits_f,xfgetdatasetunits_f_) xfgetdatasetunits_f_
#endif
XMDF_API xid func_name(XFGETDATASETUNITS_F,xfgetdatasetunits_f,xfgetdatasetunits_f_)
              (xid *a_Id, char *a_Units, int *a_unitslen)
{
  int error;
  error = xfGetDatasetUnits(*a_Id, a_Units);
  return error;

} /* xfgetdatasetunits_f_ */
  /******************************************************************************
 * FUNCTION  xfGetDatasetTimes_f_
 * PURPOSE   
 * NOTES     The times array must already be allocated to the correct number
 ******************2***********************************************************/
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFGETDATASETTIMES_F,xfgetdatasettimes_f,xfgetdatasettimes_f_) XFGETDATASETTIMES_F
#elif defined FLOWER
#define func_name(XFGETDATASETTIMES_F,xfgetdatasettimes_f,xfgetdatasettimes_f_) xfgetdatasettimes_f
#else
#define func_name(XFGETDATASETTIMES_F,xfgetdatasettimes_f,xfgetdatasettimes_f_) xfgetdatasettimes_f_
#endif
XMDF_API xid func_name(XFGETDATASETTIMES_F,xfgetdatasettimes_f,xfgetdatasettimes_f_)
              (xid *a_Id, int *a_NumTimes, double *a_Times)
{
  int error;

  error = xfGetDatasetTimes(*a_Id, *a_NumTimes, a_Times);

  return error;

} /* xfGetDatasetTimes_f_*/
/******************************************************************************
 * FUNCTION  xfGetDatasetMinsFloat_f_
 * PURPOSE   
 * NOTES     
 ******************2***********************************************************/
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFGETDATASETMINS_F,xfgetdatasetmins_f,xfgetdatasetmins_f_) XFGETDATASETMINS_F
#elif defined FLOWER
#define func_name(XFGETDATASETMINS_F,xfgetdatasetmins_f,xfgetdatasetmins_f_) xfgetdatasetmins_f
#else
#define func_name(XFGETDATASETMINS_F,xfgetdatasetmins_f,xfgetdatasetmins_f_) xfgetdatasetmins_f_
#endif
XMDF_API xid func_name(XFGETDATASETMINS_F,xfgetdatasetmins_f,xfgetdatasetmins_f_)
              (xid *a_Id, int *a_NumTimes, float *a_Mins)
{
  int error;

  error = xfGetDatasetMins(*a_Id, *a_NumTimes, a_Mins);

  return error;

} /* xfGetDatasetMins_f_*/
/******************************************************************************
 * FUNCTION  xfGetDatasetMaxs_f_
 * PURPOSE   
 * NOTES     
 ******************2***********************************************************/
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFGETDATASETMAXS_F,xfgetdatasetmaxs_f,xfgetdatasetmaxs_f_) XFGETDATASETMAXS_F
#elif defined FLOWER
#define func_name(XFGETDATASETMAXS_F,xfgetdatasetmaxs_f,xfgetdatasetmaxs_f_) xfgetdatasetmaxs_f
#else
#define func_name(XFGETDATASETMAXS_F,xfgetdatasetmaxs_f,xfgetdatasetmaxs_f_) xfgetdatasetmaxs_f_
#endif
XMDF_API xid func_name(XFGETDATASETMAXS_F,xfgetdatasetmaxs_f,xfgetdatasetmaxs_f_)
              (xid *a_Id, int *a_NumTimes, float *a_Maxs)
{
  int error;

  error = xfGetDatasetMaxs(*a_Id, *a_NumTimes, a_Maxs);

  return error;

} /* xfGetDatasetMaxs_f_*/
/******************************************************************************
 * FUNCTION  xfGetDatasetActivityTimestep_f_
 * PURPOSE   Read the activity values for a specific timestep
 * NOTES     a_Values must already be allocated of size a_NumVals     
 ******************2***********************************************************/
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFREADACTIVITYTIMESTEP_F,xfreadactivitytimestep_f,xfreadactivitytimestep_f_) XFREADACTIVITYTIMESTEP_F
#elif defined FLOWER
#define func_name(XFREADACTIVITYTIMESTEP_F,xfreadactivitytimestep_f,xfreadactivitytimestep_f_) xfreadactivitytimestep_f
#else
#define func_name(XFREADACTIVITYTIMESTEP_F,xfreadactivitytimestep_f,xfreadactivitytimestep_f_) xfreadactivitytimestep_f_
#endif
XMDF_API xid func_name(XFREADACTIVITYTIMESTEP_F,xfreadactivitytimestep_f,xfreadactivitytimestep_f_)
              (xid *a_Id, int *a_TimestepIndex, int *a_NumVals,
                                int *a_Values)
{
  
  int     error, i;

  xmbool *tmpvals;
  tmpvals = (xmbool *)malloc(*a_NumVals*sizeof(xmbool));

  /*  error = xfReadActivityTimestep(*a_Id, (*a_TimestepIndex)-1, *a_NumVals, tmpvals);*/
  error = xfReadActivityTimestep(*a_Id, *a_TimestepIndex, *a_NumVals, tmpvals);

  for (i=0; i<*a_NumVals; i=i+1) {
    a_Values[i] = (tmpvals[i]) ? XTRUE : XFALSE;
  }

  return error;

} /* xfReadActivityTimestep_f_*/
/******************************************************************************
 * FUNCTION  xfReadActivityValuesAtIndex_f_
 * PURPOSE   Read the activity values for specific indices
 * NOTES     a_Values must already be allocated of size a_NumTimes     
 ******************2***********************************************************/
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFREADACTIVITYVALUESATINDEX_F,xfreadactivityvaluesatindex_f,xfreadactivityvaluesatindex_f_) XFREADACTIVITYVALUESATINDEX_F
#elif defined FLOWER
#define func_name(XFREADACTIVITYVALUESATINDEX_F,xfreadactivityvaluesatindex_f,xfreadactivityvaluesatindex_f_) xfreadactivityvaluesatindex_f
#else
#define func_name(XFREADACTIVITYVALUESATINDEX_F,xfreadactivityvaluesatindex_f,xfreadactivityvaluesatindex_f_) xfreadactivityvaluesatindex_f_
#endif
XMDF_API xid func_name(XFREADACTIVITYVALUESATINDEX_F,xfreadactivityvaluesatindex_f,xfreadactivityvaluesatindex_f_)
              (xid *a_Id, int *a_Index, int *a_FirstTime, 
                                    int *a_NumTimes, int *a_Values)
{
  int error, i;

  xmbool *tmpvals;
  tmpvals = (xmbool *)malloc(*a_NumTimes*sizeof(xmbool));
  
  /*  error = xfReadActivityValuesAtIndex(*a_Id, (*a_Index)-1, (*a_FirstTime)-1, *a_NumTimes, */
  /*                                      tmpvals);*/
  error = xfReadActivityValuesAtIndex(*a_Id, *a_Index, *a_FirstTime, *a_NumTimes, 
                                      tmpvals);

  for (i=0; i<*a_NumTimes; i=i+1) {
    a_Values[i] = (tmpvals[i]) ? XTRUE : XFALSE;
  }


  return error;

} /* xfReadActivityValuesAtIndex_f_*/
/******************************************************************************
 * FUNCTION  xfReadScalarValuesTimestep_f_
 * PURPOSE   Read the scalar values for a specific timestep
 * NOTES     a_Values must already be allocated of size a_NumVals     
 ******************2***********************************************************/
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFREADSCALARVALUESTIMESTEP_F,xfreadscalarvaluestimestep_f,xfreadscalarvaluestimestep_f_) XFREADSCALARVALUESTIMESTEP_F
#elif defined FLOWER
#define func_name(XFREADSCALARVALUESTIMESTEP_F,xfreadscalarvaluestimestep_f,xfreadscalarvaluestimestep_f_) xfreadscalarvaluestimestep_f
#else
#define func_name(XFREADSCALARVALUESTIMESTEP_F,xfreadscalarvaluestimestep_f,xfreadscalarvaluestimestep_f_) xfreadscalarvaluestimestep_f_
#endif
XMDF_API xid func_name(XFREADSCALARVALUESTIMESTEP_F,xfreadscalarvaluestimestep_f,xfreadscalarvaluestimestep_f_)
              (xid *a_Id, int *a_TimestepIndex, int *a_NumVals, 
                                   float *a_Values)
{
  int error;

  error = xfReadScalarValuesTimestep(*a_Id, *a_TimestepIndex, *a_NumVals, a_Values);

  return error;

} /* xfReadScalarValuesTimestep_f_*/
 /******************************************************************************
 * FUNCTION  xfReadScalarValuesAtIndex_f_
 * PURPOSE   Read the scalar values for a specific index
 * NOTES     a_Values must already be allocated of size a_NumTimes     
 ******************2***********************************************************/
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFREADSCALARVALUESATINDEX_F,xfreadscalarvaluesatindex_f,xfreadscalarvaluesatindex_f_) XFREADSCALARVALUESATINDEX_F
#elif defined FLOWER
#define func_name(XFREADSCALARVALUESATINDEX_F,xfreadscalarvaluesatindex_f,xfreadscalarvaluesatindex_f_) xfreadscalarvaluesatindex_f
#else
#define func_name(XFREADSCALARVALUESATINDEX_F,xfreadscalarvaluesatindex_f,xfreadscalarvaluesatindex_f_) xfreadscalarvaluesatindex_f_
#endif
XMDF_API xid func_name(XFREADSCALARVALUESATINDEX_F,xfreadscalarvaluesatindex_f,xfreadscalarvaluesatindex_f_)
              (xid *a_Id, int *a_Index, int *a_FirstTime, 
                                  int *a_NumTimes, float *a_Values)
{
   int error;

  error = xfReadScalarValuesAtIndex(*a_Id, *a_Index, *a_FirstTime, *a_NumTimes, 
                                    a_Values);

  return error;
  
} /* xfReadScalarValuesAtIndex_f_*/
 /******************************************************************************
 * FUNCTION  xfReadScalarValuesAtIndex_f_
 * PURPOSE   Read the scalar values for a specific index
 * NOTES     a_Values must already be allocated of size a_NumTimes     
 ******************2***********************************************************/
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFREADSCALARVALUESATINDICESFLOAT_F,xfreadscalarvaluesatindicesfloat_f,xfreadscalarvaluesatindicesfloat_f_) XFREADSCALARVALUESATINDICESFLOAT_F
#elif defined FLOWER
#define func_name(XFREADSCALARVALUESATINDICESFLOAT_F,xfreadscalarvaluesatindicesfloat_f,xfreadscalarvaluesatindicesfloat_f_) xfreadscalarvaluesatindicesfloat_f
#else
#define func_name(XFREADSCALARVALUESATINDICESFLOAT_F,xfreadscalarvaluesatindicesfloat_f,xfreadscalarvaluesatindicesfloat_f_) xfreadscalarvaluesatindicesfloat_f_
#endif
XMDF_API xid func_name(XFREADSCALARVALUESATINDICESFLOAT_F,xfreadscalarvaluesatindicesfloat_f,xfreadscalarvaluesatindicesfloat_f_)
              (xid *a_Id, 
                    int *a_nIndices, const int *a_Indices, int *a_FirstTime,
                    int *a_NumTimes, float *a_Values)
{
   int error;

  error = xfReadScalarValuesAtIndicesFloat(*a_Id, *a_nIndices, a_Indices, *a_FirstTime,
                                           *a_NumTimes, a_Values);

  return error;
  
} /* xfReadScalarValuesAtIndex_f_*/
  /******************************************************************************
 * FUNCTION  xfReadVectorValuesTimestep_f_
 * PURPOSE   Read the vector values for a specific timestep
 * NOTES     a_Values must already be allocated to size 
 *           a_NumVals * a_NumComponents
 ******************2***********************************************************/
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFREADVECTORVALUESTIMESTEP_F,xfreadvectorvaluestimestep_f,xfreadvectorvaluestimestep_f_) XFREADVECTORVALUESTIMESTEP_F
#elif defined FLOWER
#define func_name(XFREADVECTORVALUESTIMESTEP_F,xfreadvectorvaluestimestep_f,xfreadvectorvaluestimestep_f_) xfreadvectorvaluestimestep_f
#else
#define func_name(XFREADVECTORVALUESTIMESTEP_F,xfreadvectorvaluestimestep_f,xfreadvectorvaluestimestep_f_) xfreadvectorvaluestimestep_f_
#endif
XMDF_API xid func_name(XFREADVECTORVALUESTIMESTEP_F,xfreadvectorvaluestimestep_f,xfreadvectorvaluestimestep_f_)
              (xid *a_Id, int *a_TimestepIndex, int *a_NumVals, 
                                   int *a_NumComponents, float *a_Values)
{
  int error;

  error = xfReadVectorValuesTimestep(*a_Id, *a_TimestepIndex, *a_NumVals, 
                                     *a_NumComponents, a_Values);

  return error;

} /* xfReadVectorValuesTimestep_f_*/
/******************************************************************************
 * FUNCTION  xfReadVectorValuesAtIndex_f_
 * PURPOSE   Read the vector values for specific indices
 * NOTES     a_Values must already be allocated of size 
 *           a_NumTimes * a_NumComponents     
 ******************2***********************************************************/
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFREADVECTORVALUESATINDEX_F,xfreadvectorvaluesatindex_f,xfreadvectorvaluesatindex_f_) XFREADVECTORVALUESATINDEX_F
#elif defined FLOWER
#define func_name(XFREADVECTORVALUESATINDEX_F,xfreadvectorvaluesatindex_f,xfreadvectorvaluesatindex_f_) xfreadvectorvaluesatindex_f
#else
#define func_name(XFREADVECTORVALUESATINDEX_F,xfreadvectorvaluesatindex_f,xfreadvectorvaluesatindex_f_) xfreadvectorvaluesatindex_f_
#endif
XMDF_API xid func_name(XFREADVECTORVALUESATINDEX_F,xfreadvectorvaluesatindex_f,xfreadvectorvaluesatindex_f_)
              (xid *a_Id, int *a_Index, int *a_FirstTime, 
                                  int *a_NumTimes, int *a_NumComponents,
                                  float *a_Values)
{
  int error;

  /*  error = xfReadVectorValuesAtIndex(*a_Id, (*a_Index)-1, (*a_FirstTime)-1, *a_NumTimes, */
  /*                                    *a_NumComponents, a_Values);*/
  error = xfReadVectorValuesAtIndex(*a_Id, *a_Index, *a_FirstTime, *a_NumTimes, 
                                    *a_NumComponents, a_Values);

  return error;

} /* xfReadVectorValuesAtIndex_f_*/
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfScalarDataLocation_f_*/
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFSCALARDATALOCATION_F,xfscalardatalocation_f,xfscalardatalocation_f_) XFSCALARDATALOCATION_F
#elif defined FLOWER
#define func_name(XFSCALARDATALOCATION_F,xfscalardatalocation_f,xfscalardatalocation_f_) xfscalardatalocation_f
#else
#define func_name(XFSCALARDATALOCATION_F,xfscalardatalocation_f,xfscalardatalocation_f_) xfscalardatalocation_f_
#endif
XMDF_API xid func_name(XFSCALARDATALOCATION_F,xfscalardatalocation_f,xfscalardatalocation_f_)
              (xid *a_Id, int *a_DataLoc)
{
  int error;

  error = xfScalarDataLocation(*a_Id, *a_DataLoc);

  return error;

 }  /* xfScalarDataLocation_f_*/
  /* --------------------------------------------------------------------------- */
/* FUNCTION  xfVector2DDataLocations_f_*/
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFVECTOR2DDATALOCATIONS_F,xfvector2ddatalocations_f,xfvector2ddatalocations_f_) XFVECTOR2DDATALOCATIONS_F
#elif defined FLOWER
#define func_name(XFVECTOR2DDATALOCATIONS_F,xfvector2ddatalocations_f,xfvector2ddatalocations_f_) xfvector2ddatalocations_f
#else
#define func_name(XFVECTOR2DDATALOCATIONS_F,xfvector2ddatalocations_f,xfvector2ddatalocations_f_) xfvector2ddatalocations_f_
#endif
XMDF_API xid func_name(XFVECTOR2DDATALOCATIONS_F,xfvector2ddatalocations_f,xfvector2ddatalocations_f_)
              (xid *a_Id, int *a_DataLocI, int *a_DataLocJ)
{
  int error;

  error = xfVector2DDataLocations(*a_Id, *a_DataLocI, *a_DataLocJ);

  return error;

} /* xfVector2DDataLocations_f_*/
 /* --------------------------------------------------------------------------- */
/* FUNCTION  xfVector3DDataLocations_f_*/
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFVECTOR3DDATALOCATIONS_F,xfvector3ddatalocations_f,xfvector3ddatalocations_f_) XFVECTOR3DDATALOCATIONS_F
#elif defined FLOWER
#define func_name(XFVECTOR3DDATALOCATIONS_F,xfvector3ddatalocations_f,xfvector3ddatalocations_f_) xfvector3ddatalocations_f
#else
#define func_name(XFVECTOR3DDATALOCATIONS_F,xfvector3ddatalocations_f,xfvector3ddatalocations_f_) xfvector3ddatalocations_f_
#endif
XMDF_API xid func_name(XFVECTOR3DDATALOCATIONS_F,xfvector3ddatalocations_f,xfvector3ddatalocations_f_)
              (xid *a_Id, int *a_DataLocI, int *a_DataLocJ,
                                int *a_DataLocK)
{
  int error;

  error = xfVector3DDataLocations(*a_Id, *a_DataLocI, *a_DataLocJ, *a_DataLocJ);

  return error;

} /* xfVector3DDataLocations_f_*/
 /* --------------------------------------------------------------------------- */
/* FUNCTION  xfGetScalarDataLocation_f_*/
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFGETSCALARDATALOCATION_F,xfgetscalardatalocation_f,xfgetscalardatalocation_f_) XFGETSCALARDATALOCATION_F
#elif defined FLOWER
#define func_name(XFGETSCALARDATALOCATION_F,xfgetscalardatalocation_f,xfgetscalardatalocation_f_) xfgetscalardatalocation_f
#else
#define func_name(XFGETSCALARDATALOCATION_F,xfgetscalardatalocation_f,xfgetscalardatalocation_f_) xfgetscalardatalocation_f_
#endif
XMDF_API xid func_name(XFGETSCALARDATALOCATION_F,xfgetscalardatalocation_f,xfgetscalardatalocation_f_)
              (xid *a_Id, int *a_DataLoc)
{
  int error;

  error = xfGetScalarDataLocation(*a_Id, a_DataLoc);

  return error;

} /* xfGetScalarDataLocation_f_*/
 /* --------------------------------------------------------------------------- */
/* FUNCTION  xfGetVector2DDataLocations_f_*/
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFGETVECTOR2DDATALOCATIONS_F,xfgetvector2ddatalocations_f,xfgetvector2ddatalocations_f_) XFGETVECTOR2DDATALOCATIONS_F
#elif defined FLOWER
#define func_name(XFGETVECTOR2DDATALOCATIONS_F,xfgetvector2ddatalocations_f,xfgetvector2ddatalocations_f_) xfgetvector2ddatalocations_f
#else
#define func_name(XFGETVECTOR2DDATALOCATIONS_F,xfgetvector2ddatalocations_f,xfgetvector2ddatalocations_f_) xfgetvector2ddatalocations_f_
#endif
XMDF_API xid func_name(XFGETVECTOR2DDATALOCATIONS_F,xfgetvector2ddatalocations_f,xfgetvector2ddatalocations_f_)
              (xid *a_Id, int *a_DataLocI, int *a_DataLocJ)
{
  int error;

  error = xfGetVector2DDataLocations(*a_Id, a_DataLocI, a_DataLocJ);

  return error;

} /* xfGetVector2DDataLocations_f_*/
 /* --------------------------------------------------------------------------- */
/* FUNCTION  xfGetVector3DDataLocations_f_*/
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFGETVECTOR3DDATALOCATIONS_F,xfgetvector3ddatalocations_f,xfgetvector3ddatalocations_f_) XFGETVECTOR3DDATALOCATIONS_F
#elif defined FLOWER
#define func_name(XFGETVECTOR3DDATALOCATIONS_F,xfgetvector3ddatalocations_f,xfgetvector3ddatalocations_f_) xfgetvector3ddatalocations_f
#else
#define func_name(XFGETVECTOR3DDATALOCATIONS_F,xfgetvector3ddatalocations_f,xfgetvector3ddatalocations_f_) xfgetvector3ddatalocations_f_
#endif
XMDF_API xid func_name(XFGETVECTOR3DDATALOCATIONS_F,xfgetvector3ddatalocations_f,xfgetvector3ddatalocations_f_)
              (xid *a_Id, int *a_DataLocI,
                                   int *a_DataLocJ, int *a_DataLocK)
{
  int error;

  error = xfGetVector3DDataLocations(*a_Id, a_DataLocI, a_DataLocJ, a_DataLocK);

  return error;

} /* xfGetVector3DDataLocations_f_*/
 /* --------------------------------------------------------------------------- */
/* FUNCTION  xfVectorsInLocalCoords_f_*/
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFVECTORSINLOCALCOORDS_F,xfvectorsinlocalcoords_f,xfvectorsinlocalcoords_f_) XFVECTORSINLOCALCOORDS_F
#elif defined FLOWER
#define func_name(XFVECTORSINLOCALCOORDS_F,xfvectorsinlocalcoords_f,xfvectorsinlocalcoords_f_) xfvectorsinlocalcoords_f
#else
#define func_name(XFVECTORSINLOCALCOORDS_F,xfvectorsinlocalcoords_f,xfvectorsinlocalcoords_f_) xfvectorsinlocalcoords_f_
#endif
XMDF_API xid func_name(XFVECTORSINLOCALCOORDS_F,xfvectorsinlocalcoords_f,xfvectorsinlocalcoords_f_)
              (xid *a_Id)
{
  int error;

  error = xfVectorsInLocalCoords(*a_Id);

  return error;

} /* xfVectorsInLocalCoords_f_*/
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfAreVectorsInLocalCoords_f_*/
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFAREVECTORSINLOCALCOORDS_F,xfarevectorsinlocalcoords_f,xfarevectorsinlocalcoords_f_) XFAREVECTORSINLOCALCOORDS_F
#elif defined FLOWER
#define func_name(XFAREVECTORSINLOCALCOORDS_F,xfarevectorsinlocalcoords_f,xfarevectorsinlocalcoords_f_) xfarevectorsinlocalcoords_f
#else
#define func_name(XFAREVECTORSINLOCALCOORDS_F,xfarevectorsinlocalcoords_f,xfarevectorsinlocalcoords_f_) xfarevectorsinlocalcoords_f_
#endif
XMDF_API xid func_name(XFAREVECTORSINLOCALCOORDS_F,xfarevectorsinlocalcoords_f,xfarevectorsinlocalcoords_f_)
              (xid *a_Id, int *a_LocalCoords)
{
  int error;

  error = xfAreVectorsInLocalCoords(*a_Id, a_LocalCoords);

  return error;

} /* xfAreVectorsInLocalCoords_f_*/

/* --------------------------------------------------------------------------- */
/* FUNCTION  xfGetHorizDatum_f_*/
/* PURPOSE   Read the horizontal datum (#defines in xmdf.h) */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFGETHORIZDATUM_F,xfgethorizdatum_f,xfgethorizdatum_f_) XFGETHORIZDATUM_F
#elif defined FLOWER
#define func_name(XFGETHORIZDATUM_F,xfgethorizdatum_f,xfgethorizdatum_f_) xfgethorizdatum_f
#else
#define func_name(XFGETHORIZDATUM_F,xfgethorizdatum_f,xfgethorizdatum_f_) xfgethorizdatum_f_
#endif
XMDF_API xid func_name(XFGETHORIZDATUM_F,xfgethorizdatum_f,xfgethorizdatum_f_)
              (xid *a_CoordId, int *a_val)
{
  int error;

  error = xfGetHorizDatum(*a_CoordId, a_val);

  return error;

} /* xfGetHorizDatum_f_*/
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfGetHorizUnits_f_*/
/* PURPOSE   Read the horizontal units (#defines in xmdf.h) */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFGETHORIZUNITS_F,xfgethorizunits_f,xfgethorizunits_f_) XFGETHORIZUNITS_F
#elif defined FLOWER
#define func_name(XFGETHORIZUNITS_F,xfgethorizunits_f,xfgethorizunits_f_) xfgethorizunits_f
#else
#define func_name(XFGETHORIZUNITS_F,xfgethorizunits_f,xfgethorizunits_f_) xfgethorizunits_f_
#endif
XMDF_API xid func_name(XFGETHORIZUNITS_F,xfgethorizunits_f,xfgethorizunits_f_)
              (xid *a_CoordId, int *a_val)
{
  int error;

  error = xfGetHorizUnits(*a_CoordId, a_val);

  return error;

} /* xfGetHorizUnits_f_*/
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfGetVertDatum_f_*/
/* PURPOSE   Read the vertical datum (#defines in xmdf.h) */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFGETVERTDATUM_F,xfgetvertdatum_f,xfgetvertdatum_f_) XFGETVERTDATUM_F
#elif defined FLOWER
#define func_name(XFGETVERTDATUM_F,xfgetvertdatum_f,xfgetvertdatum_f_) xfgetvertdatum_f
#else
#define func_name(XFGETVERTDATUM_F,xfgetvertdatum_f,xfgetvertdatum_f_) xfgetvertdatum_f_
#endif
XMDF_API xid func_name(XFGETVERTDATUM_F,xfgetvertdatum_f,xfgetvertdatum_f_)
              (xid *a_CoordId, int *a_val)
{
  int error;

  error = xfGetVertDatum(*a_CoordId, a_val);

  return error;

} /* xfGetVertDatum_f_*/
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfGetVertUnits_f_*/
/* PURPOSE   Read the vertical units (#defines in xmdf.h) */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFGETVERTUNITS_F,xfgetvertunits_f,xfgetvertunits_f_) XFGETVERTUNITS_F
#elif defined FLOWER
#define func_name(XFGETVERTUNITS_F,xfgetvertunits_f,xfgetvertunits_f_) xfgetvertunits_f
#else
#define func_name(XFGETVERTUNITS_F,xfgetvertunits_f,xfgetvertunits_f_) xfgetvertunits_f_
#endif
XMDF_API xid func_name(XFGETVERTUNITS_F,xfgetvertunits_f,xfgetvertunits_f_)
              (xid *a_CoordId, int *a_val)
{
  int error;

  error = xfGetVertUnits(*a_CoordId, a_val);

  return error;

} /* xfGetVertUnits_f_*/
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfGetLat_f_*/
/* PURPOSE   Read whether the lattitude is North or South latitude */
/*           (#defines in xmdf.h) */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFGETLAT_F,xfgetlat_f,xfgetlat_f_) XFGETLAT_F
#elif defined FLOWER
#define func_name(XFGETLAT_F,xfgetlat_f,xfgetlat_f_) xfgetlat_f
#else
#define func_name(XFGETLAT_F,xfgetlat_f,xfgetlat_f_) xfgetlat_f_
#endif
XMDF_API xid func_name(XFGETLAT_F,xfgetlat_f,xfgetlat_f_)
              (xid *a_CoordId, int *a_val)
{
  int error;

  error = xfGetLat(*a_CoordId, a_val);

  return error;

} /* xfGetLat_f_*/
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfGetLon_f_*/
/* PURPOSE   Read whether the longitude is East or West  */
/*           (#defines in xmdf.h) */
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFGETLON_F,xfgetlon_f,xfgetlon_f_) XFGETLON_F
#elif defined FLOWER
#define func_name(XFGETLON_F,xfgetlon_f,xfgetlon_f_) xfgetlon_f
#else
#define func_name(XFGETLON_F,xfgetlon_f,xfgetlon_f_) xfgetlon_f_
#endif
XMDF_API xid func_name(XFGETLON_F,xfgetlon_f,xfgetlon_f_)
              (xid *a_CoordId, int *a_val)
{
  int error;

  error = xfGetLon(*a_CoordId, a_val);

  return error;

} /* xfGetLon_f_*/
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfGetUTMZone_f_*/
/* PURPOSE   Read the UTM zone (should be a number between 1 and 60) */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFGETUTMZONE_F,xfgetutmzone_f,xfgetutmzone_f_) XFGETUTMZONE_F
#elif defined FLOWER
#define func_name(XFGETUTMZONE_F,xfgetutmzone_f,xfgetutmzone_f_) xfgetutmzone_f
#else
#define func_name(XFGETUTMZONE_F,xfgetutmzone_f,xfgetutmzone_f_) xfgetutmzone_f_
#endif
XMDF_API xid func_name(XFGETUTMZONE_F,xfgetutmzone_f,xfgetutmzone_f_)
              (xid *a_CoordId, int *a_val)
{
  int error;

  error = xfGetUTMZone(*a_CoordId, a_val);

  return error;

} /* xfGetUTMZone_f_*/
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfGetSPCZone_f_*/
/* PURPOSE   Read the SPC zone (Lookup numbers in XMDF documentation) */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFGETSPCZONE_F,xfgetspczone_f,xfgetspczone_f_) XFGETSPCZONE_F
#elif defined FLOWER
#define func_name(XFGETSPCZONE_F,xfgetspczone_f,xfgetspczone_f_) xfgetspczone_f
#else
#define func_name(XFGETSPCZONE_F,xfgetspczone_f,xfgetspczone_f_) xfgetspczone_f_
#endif
XMDF_API xid func_name(XFGETSPCZONE_F,xfgetspczone_f,xfgetspczone_f_)
              (xid *a_CoordId, int *a_val)
{
  int error;

  error = xfGetSPCZone(*a_CoordId, a_val);

  return error;

} /* xfGetSPCZone_f_*/
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfGetHPGNArea_f_*/
/* PURPOSE   Read the HPGN area (Lookup numbers in XMDF documentation) */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFGETHPGNAREA_F,xfgethpgnarea_f,xfgethpgnarea_f_) XFGETHPGNAREA_F
#elif defined FLOWER
#define func_name(XFGETHPGNAREA_F,xfgethpgnarea_f,xfgethpgnarea_f_) xfgethpgnarea_f
#else
#define func_name(XFGETHPGNAREA_F,xfgethpgnarea_f,xfgethpgnarea_f_) xfgethpgnarea_f_
#endif
XMDF_API xid func_name(XFGETHPGNAREA_F,xfgethpgnarea_f,xfgethpgnarea_f_)
              (xid *a_CoordId, int *a_val)
{
  int error;

  error = xfGetHPGNArea(*a_CoordId, a_val);

  return error;

} /* xfGetHPGNArea_f_*/
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfGetCPPLat_f_*/
/* PURPOSE   Read the Carte Prallelo Grammatique Projetion Factor for  */
/*           converting the latitude */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFGETCPPLAT_F,xfgetcpplat_f,xfgetcpplat_f_) XFGETCPPLAT_F
#elif defined FLOWER
#define func_name(XFGETCPPLAT_F,xfgetcpplat_f,xfgetcpplat_f_) xfgetcpplat_f
#else
#define func_name(XFGETCPPLAT_F,xfgetcpplat_f,xfgetcpplat_f_) xfgetcpplat_f_
#endif
XMDF_API xid func_name(XFGETCPPLAT_F,xfgetcpplat_f,xfgetcpplat_f_)
              (xid *a_CoordId, double *a_val)
{
  int error;

  error = xfGetCPPLat(*a_CoordId, a_val);

  return error;

} /* xfGetCPPLat_f_*/
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfGetCPPLon_f_*/
/* PURPOSE   Read the Carte Prallelo Grammatique Projetion Factor for  */
/*           converting the longitude */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFGETCPPLON_F,xfgetcpplon_f,xfgetcpplon_f_) XFGETCPPLON_F
#elif defined FLOWER
#define func_name(XFGETCPPLON_F,xfgetcpplon_f,xfgetcpplon_f_) xfgetcpplon_f
#else
#define func_name(XFGETCPPLON_F,xfgetcpplon_f,xfgetcpplon_f_) xfgetcpplon_f_
#endif
XMDF_API xid func_name(XFGETCPPLON_F,xfgetcpplon_f,xfgetcpplon_f_)
              (xid *a_CoordId, double *a_val)
{
  int error;

  error = xfGetCPPLon(*a_CoordId, a_val);

  return error;

} /* xfGetCPPLon_f_*/
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfGetEllipse_f_*/
/* PURPOSE   Read the Ellipse number based upon XMDF documentation */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFGETELLIPSE_F,xfgetellipse_f,xfgetellipse_f_) XFGETELLIPSE_F
#elif defined FLOWER
#define func_name(XFGETELLIPSE_F,xfgetellipse_f,xfgetellipse_f_) xfgetellipse_f
#else
#define func_name(XFGETELLIPSE_F,xfgetellipse_f,xfgetellipse_f_) xfgetellipse_f_
#endif
XMDF_API xid func_name(XFGETELLIPSE_F,xfgetellipse_f,xfgetellipse_f_)
              (xid *a_CoordId, int *a_val)
{
  int error;

  error = xfGetEllipse(*a_CoordId, a_val);

  return error;

} /* xfGetEllipse_f_*/
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfGetMajorR_f_*/
/* PURPOSE   Read a user-defined ellipse major radius */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFGETMAJORR_F,xfgetmajorr_f,xfgetmajorr_f_) XFGETMAJORR_F
#elif defined FLOWER
#define func_name(XFGETMAJORR_F,xfgetmajorr_f,xfgetmajorr_f_) xfgetmajorr_f
#else
#define func_name(XFGETMAJORR_F,xfgetmajorr_f,xfgetmajorr_f_) xfgetmajorr_f_
#endif
XMDF_API xid func_name(XFGETMAJORR_F,xfgetmajorr_f,xfgetmajorr_f_)
              (xid *a_CoordId, double *a_val)
{
  int error;

  error = xfGetMajorR(*a_CoordId, a_val);

  return error;

} /* xfGetMajorR_f_*/
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfGetMinorR_f_*/
/* PURPOSE   Read a user-defined ellipse minor radius */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFGETMINORR_F,xfgetminorr_f,xfgetminorr_f_) XFGETMINORR_F
#elif defined FLOWER
#define func_name(XFGETMINORR_F,xfgetminorr_f,xfgetminorr_f_) xfgetminorr_f
#else
#define func_name(XFGETMINORR_F,xfgetminorr_f,xfgetminorr_f_) xfgetminorr_f_
#endif
XMDF_API xid func_name(XFGETMINORR_F,xfgetminorr_f,xfgetminorr_f_)
              (xid *a_CoordId, double *a_val)
{
  int error;

  error = xfGetMinorR(*a_CoordId, a_val);

  return error;

} /* xfGetMinorR_f_*/
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfSetHorizDatum_f_*/
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFSETHORIZDATUM_F,xfsethorizdatum_f,xfsethorizdatum_f_) XFSETHORIZDATUM_F
#elif defined FLOWER
#define func_name(XFSETHORIZDATUM_F,xfsethorizdatum_f,xfsethorizdatum_f_) xfsethorizdatum_f
#else
#define func_name(XFSETHORIZDATUM_F,xfsethorizdatum_f,xfsethorizdatum_f_) xfsethorizdatum_f_
#endif
XMDF_API xid func_name(XFSETHORIZDATUM_F,xfsethorizdatum_f,xfsethorizdatum_f_)
              (xid *a_CoordId, int *a_val)
{
  int error;

  error = xfSetHorizDatum(*a_CoordId, *a_val);

  return error;

} /* xfSetHorizDatum_f_*/
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfSetHorizUnits_f_*/
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFSETHORIZUNITS_F,xfsethorizunits_f,xfsethorizunits_f_) XFSETHORIZUNITS_F
#elif defined FLOWER
#define func_name(XFSETHORIZUNITS_F,xfsethorizunits_f,xfsethorizunits_f_) xfsethorizunits_f
#else
#define func_name(XFSETHORIZUNITS_F,xfsethorizunits_f,xfsethorizunits_f_) xfsethorizunits_f_
#endif
XMDF_API xid func_name(XFSETHORIZUNITS_F,xfsethorizunits_f,xfsethorizunits_f_)
              (xid *a_CoordId, int *a_val)
{
  int error;

  error = xfSetHorizUnits(*a_CoordId, *a_val);

  return error;

} /* xfSetHorizUnits_f_*/
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfSetVertDatum_f_*/
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFSETVERTDATUM_F,xfsetvertdatum_f,xfsetvertdatum_f_) XFSETVERTDATUM_F
#elif defined FLOWER
#define func_name(XFSETVERTDATUM_F,xfsetvertdatum_f,xfsetvertdatum_f_) xfsetvertdatum_f
#else
#define func_name(XFSETVERTDATUM_F,xfsetvertdatum_f,xfsetvertdatum_f_) xfsetvertdatum_f_
#endif
XMDF_API xid func_name(XFSETVERTDATUM_F,xfsetvertdatum_f,xfsetvertdatum_f_)
              (xid *a_CoordId, int *a_val)
{
  int error;

  error = xfSetVertDatum(*a_CoordId, *a_val);

  return error;

} /* xfSetVertDatum_f_*/
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfSetVertUnits */
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFSETVERTUNITS_F,xfsetvertunits_f,xfsetvertunits_f_) XFSETVERTUNITS_F
#elif defined FLOWER
#define func_name(XFSETVERTUNITS_F,xfsetvertunits_f,xfsetvertunits_f_) xfsetvertunits_f
#else
#define func_name(XFSETVERTUNITS_F,xfsetvertunits_f,xfsetvertunits_f_) xfsetvertunits_f_
#endif
XMDF_API xid func_name(XFSETVERTUNITS_F,xfsetvertunits_f,xfsetvertunits_f_)
              (xid *a_CoordId, int *a_val)
{
  int error;

  error = xfSetVertUnits(*a_CoordId, *a_val);

  return error;

} /* xfSetVertUnits_f_*/
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfSetLat_f_*/
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFSETLAT_F,xfsetlat_f,xfsetlat_f_) XFSETLAT_F
#elif defined FLOWER
#define func_name(XFSETLAT_F,xfsetlat_f,xfsetlat_f_) xfsetlat_f
#else
#define func_name(XFSETLAT_F,xfsetlat_f,xfsetlat_f_) xfsetlat_f_
#endif
XMDF_API xid func_name(XFSETLAT_F,xfsetlat_f,xfsetlat_f_)
              (xid *a_CoordId, int *a_val)
{
  int error;

  error = xfSetLat(*a_CoordId, *a_val);

  return error;

} /* xfSetLat_f_*/
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfSetLon_f_*/
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFSETLON_F,xfsetlon_f,xfsetlon_f_) XFSETLON_F
#elif defined FLOWER
#define func_name(XFSETLON_F,xfsetlon_f,xfsetlon_f_) xfsetlon_f
#else
#define func_name(XFSETLON_F,xfsetlon_f,xfsetlon_f_) xfsetlon_f_
#endif
XMDF_API xid func_name(XFSETLON_F,xfsetlon_f,xfsetlon_f_)
              (xid *a_CoordId, int *a_val)
{
  int error;

  error = xfSetLon(*a_CoordId, *a_val);

  return error;

} /* xfSetLon_f_*/
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfSetUTMZone_f_*/
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFSETUTMZONE_F,xfsetutmzone_f,xfsetutmzone_f_) XFSETUTMZONE_F
#elif defined FLOWER
#define func_name(XFSETUTMZONE_F,xfsetutmzone_f,xfsetutmzone_f_) xfsetutmzone_f
#else
#define func_name(XFSETUTMZONE_F,xfsetutmzone_f,xfsetutmzone_f_) xfsetutmzone_f_
#endif
XMDF_API xid func_name(XFSETUTMZONE_F,xfsetutmzone_f,xfsetutmzone_f_)
              (xid *a_CoordId, int *a_val)
{
  int error;

  error = xfSetUTMZone(*a_CoordId, *a_val);

  return error;

} /* xfSetUTMZone_f_*/
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfSetSPCZone_f_*/
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFSETSPCZONE_F,xfsetspczone_f,xfsetspczone_f_) XFSETSPCZONE_F
#elif defined FLOWER
#define func_name(XFSETSPCZONE_F,xfsetspczone_f,xfsetspczone_f_) xfsetspczone_f
#else
#define func_name(XFSETSPCZONE_F,xfsetspczone_f,xfsetspczone_f_) xfsetspczone_f_
#endif
XMDF_API xid func_name(XFSETSPCZONE_F,xfsetspczone_f,xfsetspczone_f_)
              (xid *a_CoordId, int *a_val)
{
  int error;

  error = xfSetSPCZone(*a_CoordId, *a_val);

  return error;

} /* xfSetSPCZone_f_*/
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfSetHPGNArea_f_*/
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFSETHPGNAREA_F,xfsethpgnarea_f,xfsethpgnarea_f_) XFSETHPGNAREA_F
#elif defined FLOWER
#define func_name(XFSETHPGNAREA_F,xfsethpgnarea_f,xfsethpgnarea_f_) xfsethpgnarea_f
#else
#define func_name(XFSETHPGNAREA_F,xfsethpgnarea_f,xfsethpgnarea_f_) xfsethpgnarea_f_
#endif
XMDF_API xid func_name(XFSETHPGNAREA_F,xfsethpgnarea_f,xfsethpgnarea_f_)
              (xid *a_CoordId, int *a_val)
{
  int error;

  error = xfSetHPGNArea(*a_CoordId, *a_val);

  return error;

} /* xfSetHPGNArea_f_*/
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfSetCPPLat_f_*/
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFSETCPPLAT_F,xfsetcpplat_f,xfsetcpplat_f_) XFSETCPPLAT_F
#elif defined FLOWER
#define func_name(XFSETCPPLAT_F,xfsetcpplat_f,xfsetcpplat_f_) xfsetcpplat_f
#else
#define func_name(XFSETCPPLAT_F,xfsetcpplat_f,xfsetcpplat_f_) xfsetcpplat_f_
#endif
XMDF_API xid func_name(XFSETCPPLAT_F,xfsetcpplat_f,xfsetcpplat_f_)
              (xid *a_CoordId, double *a_val)
{
  int error;

  error = xfSetCPPLat(*a_CoordId, *a_val);

  return error;

} /* xfSetCPPLat_f_*/
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfSetCPPLon_f_*/
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFSETCPPLON_F,xfsetcpplon_f,xfsetcpplon_f_) XFSETCPPLON_F
#elif defined FLOWER
#define func_name(XFSETCPPLON_F,xfsetcpplon_f,xfsetcpplon_f_) xfsetcpplon_f
#else
#define func_name(XFSETCPPLON_F,xfsetcpplon_f,xfsetcpplon_f_) xfsetcpplon_f_
#endif
XMDF_API xid func_name(XFSETCPPLON_F,xfsetcpplon_f,xfsetcpplon_f_)
              (xid *a_CoordId, double *a_val)
{
  int error;

  error = xfSetCPPLon(*a_CoordId, *a_val);

  return error;

} /* xfSetCPPLon_f_*/
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfSetEllipse_f_*/
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFSETELLIPSE_F,xfsetellipse_f,xfsetellipse_f_) XFSETELLIPSE_F
#elif defined FLOWER
#define func_name(XFSETELLIPSE_F,xfsetellipse_f,xfsetellipse_f_) xfsetellipse_f
#else
#define func_name(XFSETELLIPSE_F,xfsetellipse_f,xfsetellipse_f_) xfsetellipse_f_
#endif
XMDF_API xid func_name(XFSETELLIPSE_F,xfsetellipse_f,xfsetellipse_f_)
              (xid *a_CoordId, int *a_val)
{
  int error;

  error = xfSetEllipse(*a_CoordId, *a_val);

  return error;

} /* xfSetEllipse_f_*/
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfSetMajorR_f_*/
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFSETMAJORR_F,xfsetmajorr_f,xfsetmajorr_f_) XFSETMAJORR_F
#elif defined FLOWER
#define func_name(XFSETMAJORR_F,xfsetmajorr_f,xfsetmajorr_f_) xfsetmajorr_f
#else
#define func_name(XFSETMAJORR_F,xfsetmajorr_f,xfsetmajorr_f_) xfsetmajorr_f_
#endif
XMDF_API xid func_name(XFSETMAJORR_F,xfsetmajorr_f,xfsetmajorr_f_)
              (xid *a_CoordId, double *a_val)
{
  int error;

  error = xfSetMajorR(*a_CoordId, *a_val);

  return error;

} /* xfSetMajorR_f_*/
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfSetMinorR_f_*/
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFSETMINORR_F,xfsetminorr_f,xfsetminorr_f_) XFSETMINORR_F
#elif defined FLOWER
#define func_name(XFSETMINORR_F,xfsetminorr_f,xfsetminorr_f_) xfsetminorr_f
#else
#define func_name(XFSETMINORR_F,xfsetminorr_f,xfsetminorr_f_) xfsetminorr_f_
#endif
XMDF_API xid func_name(XFSETMINORR_F,xfsetminorr_f,xfsetminorr_f_)
              (xid *a_CoordId, double *a_val)
{
  int error;

  error = xfSetMinorR(*a_CoordId, *a_val);

  return error;

} /* xfSetMinorR_f_*/
 /*----------------------------------------------------------------------------- */
/* FUNCTION  xfCalendarToJulian_f_*/
/* PURPOSE   Convert Calendar to Julian Date */
/* NOTES     Returns 1 if successful and -1 if unsuccessful */
/*           (Julian day number algorithm adopted from Press et al.) */
/*           era = 0 represents BCE (BC), and era = 1 represents CE (AD). */
/*          -Taken from JavaScript Code found at website: */
/*                http://aa.usno.navy.mil/data/docs/JulianDate.html */
/*           Contact info. provided by website: */
/*             Marc A. Murison */
/*             Astronomical Applications Dept. */
/*             U.S. Naval Observatory */
/*             3450 Massachusetts Ave, NW */
/*             Washington, DC  20392-5420 */
/*---------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFCALENDARTOJULIAN_F,xfcalendartojulian_f,xfcalendartojulian_f_) XFCALENDARTOJULIAN_F
#elif defined FLOWER
#define func_name(XFCALENDARTOJULIAN_F,xfcalendartojulian_f,xfcalendartojulian_f_) xfcalendartojulian_f
#else
#define func_name(XFCALENDARTOJULIAN_F,xfcalendartojulian_f,xfcalendartojulian_f_) xfcalendartojulian_f_
#endif
XMDF_API xid func_name(XFCALENDARTOJULIAN_F,xfcalendartojulian_f,xfcalendartojulian_f_)
              (int *a_bEra, int *a_yr, int *a_mo, int *a_day, 
                          int *a_hr, int *a_min, int *a_sec, double *a_julian) 
{
  int error;
  xmbool bEra = (*a_bEra) ? XTRUE : XFALSE;
  error = xfCalendarToJulian(bEra, *a_yr, *a_mo, *a_day, *a_hr, *a_min, *a_sec, 
                             a_julian);
  return error;

} /* xfCalendarToJulian_f_*/

/*---------------------------------------------------------------------------- */
/* FUNCTION  xfJulianToCalendar_f_*/
/* PURPOSE   Convert Julian Date To Calendar */
/* NOTES     Returns 1 if successful and -1 if unsuccessful */
/*           (algorithm adopted from Press et al.) */
/*           era = 0 represents BCE (BC), and era = 1 represents CE (AD). */
/*          -Taken from JavaScript Code found at website: */
/*                http://aa.usno.navy.mil/data/docs/JulianDate.html */
/*           Contact info. provided by website: */
/*             Marc A. Murison */
/*             Astronomical Applications Dept. */
/*             U.S. Naval Observatory */
/*             3450 Massachusetts Ave, NW */
/*             Washington, DC  20392-5420 */
/*---------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFJULIANTOCALENDAR_F,xfjuliantocalendar_f,xfjuliantocalendar_f_) XFJULIANTOCALENDAR_F
#elif defined FLOWER
#define func_name(XFJULIANTOCALENDAR_F,xfjuliantocalendar_f,xfjuliantocalendar_f_) xfjuliantocalendar_f
#else
#define func_name(XFJULIANTOCALENDAR_F,xfjuliantocalendar_f,xfjuliantocalendar_f_) xfjuliantocalendar_f_
#endif
XMDF_API xid func_name(XFJULIANTOCALENDAR_F,xfjuliantocalendar_f,xfjuliantocalendar_f_)
              (int *a_bEra, int *a_yr, int *a_mo, int *a_day,
                          int *a_hr, int *a_min, int *a_sec, double *a_julian)
{
  int error;
  xmbool bEra;

  error = xfJulianToCalendar(&bEra, a_yr, a_mo, a_day, a_hr, a_min, a_sec, 
                             *a_julian);
  
  *a_bEra = (bEra) ? XTRUE : XFALSE;

  return error;

} /* xfJulianToCalendar_f_*/
/******************************************************************************
 * FUNCTION  xfDatasetReftime_f_
 * PURPOSE   Set a reference time for a dataset
 * NOTES     
 ******************2***********************************************************/
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFDATASETREFTIME_F,xfdatasetreftime_f,xfdatasetreftime_f_) XFDATASETREFTIME_F
#elif defined FLOWER
#define func_name(XFDATASETREFTIME_F,xfdatasetreftime_f,xfdatasetreftime_f_) xfdatasetreftime_f
#else
#define func_name(XFDATASETREFTIME_F,xfdatasetreftime_f,xfdatasetreftime_f_) xfdatasetreftime_f_
#endif
XMDF_API xid func_name(XFDATASETREFTIME_F,xfdatasetreftime_f,xfdatasetreftime_f_)
              (xid *a_Id, double *a_Reftime)
{

  int error;

  error = xfDatasetReftime (*a_Id, *a_Reftime);

  return error;

} /* xfDatasetReftime_f_*/

 /******************************************************************************
 * FUNCTION  xfGetNumOpenIdentifiers_f_
 * PURPOSE   Get the number of open identifiers in a file.
 * NOTES     
 ******************2***********************************************************/
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFGETNUMOPENIDENTIFIERS_F,xfgetnumopenidentifiers_f,xfgetnumopenidentifiers_f_) XFGETNUMOPENIDENTIFIERS_F
#elif defined FLOWER
#define func_name(XFGETNUMOPENIDENTIFIERS_F,xfgetnumopenidentifiers_f,xfgetnumopenidentifiers_f_) xfgetnumopenidentifiers_f
#else
#define func_name(XFGETNUMOPENIDENTIFIERS_F,xfgetnumopenidentifiers_f,xfgetnumopenidentifiers_f_) xfgetnumopenidentifiers_f_
#endif
XMDF_API xid func_name(XFGETNUMOPENIDENTIFIERS_F,xfgetnumopenidentifiers_f,xfgetnumopenidentifiers_f_)
              (xid *a_Id, int *a_Num)
{ 

  int error;
//int ii;
//DBG;
//for (ii = -5; ii <= 5; ++ii) Dbg("before",(&error)[ii]);
  error = xfGetNumOpenIdentifiers(*a_Id, a_Num);
//for (ii = -5; ii <= 5; ++ii) Dbg("after",(&error)[ii]);
//Dbg("error",error);
//DBG;
//__asm int 3
  return error;

} /* xfGetNumOpenIdentifiers_f_*/
 /******************************************************************************
 * FUNCTION  xfpCloseOpenIdentifiers_f_
 * PURPOSE   Closes all open identifiers for a file 
 * NOTES     
 ******************2***********************************************************/
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFPCLOSEOPENIDENTIFIERS_F,xfpcloseopenidentifiers_f,xfpcloseopenidentifiers_f_) XFPCLOSEOPENIDENTIFIERS_F
#elif defined FLOWER
#define func_name(XFPCLOSEOPENIDENTIFIERS_F,xfpcloseopenidentifiers_f,xfpcloseopenidentifiers_f_) xfpcloseopenidentifiers_f
#else
#define func_name(XFPCLOSEOPENIDENTIFIERS_F,xfpcloseopenidentifiers_f,xfpcloseopenidentifiers_f_) xfpcloseopenidentifiers_f_
#endif
XMDF_API xid func_name(XFPCLOSEOPENIDENTIFIERS_F,xfpcloseopenidentifiers_f,xfpcloseopenidentifiers_f_)
              (xid *a_Id)
{

  int error;

  error = xfpCloseOpenIdentifiers(*a_Id);

  return error;

} /* xfpCloseOpenIdentifiers_f_*/
 /******************************************************************************
 * FUNCTION  xfWritePropertyString
 * PURPOSE   writes a string attribute (really a dataset) to the folder
 * NOTES     
 ******************2***********************************************************/
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFWRITEPROPERTYSTRING_F,xfwritepropertystring_f,xfwritepropertystring_f_) XFWRITEPROPERTYSTRING_F
#elif defined FLOWER
#define func_name(XFWRITEPROPERTYSTRING_F,xfwritepropertystring_f,xfwritepropertystring_f_) xfwritepropertystring_f
#else
#define func_name(XFWRITEPROPERTYSTRING_F,xfwritepropertystring_f,xfwritepropertystring_f_) xfwritepropertystring_f_
#endif
XMDF_API xid func_name(XFWRITEPROPERTYSTRING_F,xfwritepropertystring_f,xfwritepropertystring_f_)
              (xid *a_Id, const char *a_Name, int* a_StringLength,
                              int* a_Number, const char *a_Attributes,
                              int *a_tmplen)
{ 
  int error;
  char *name_copy = NULL, *att_copy = NULL;

  name_copy = (char*)malloc((*a_StringLength+1)*sizeof(char));
  att_copy = (char*)malloc((*a_tmplen+1)*sizeof(char));
  strncpy(name_copy, a_Name, *a_StringLength);
  strncpy(att_copy, a_Attributes, *a_tmplen);
  name_copy[*a_StringLength] = '\0';
  att_copy[*a_tmplen] = '\0';

  error = xfWritePropertyString(*a_Id, name_copy, *a_Number, *a_tmplen, att_copy);
  
  free(name_copy);
  free(att_copy);
  return error;

} /* xfWritePropertyString */

/* --------------------------------------------------------------------------- */
/* FUNCTION  xfGetMetadataXML_f */
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFGETMETADATAXML_F,xfGetMetadataXML_f,xfGetMetadataXML_f_ ) XFGETMETADATAXML_F
#elif defined FLOWER
#define func_name(XFGETMETADATAXML_F,xfGetMetadataXML_f,xfGetMetadataXML_f_ ) xfGetMetadataXML_f
#else
#define func_name(XFGETMETADATAXML_F,xfGetMetadataXML_f,xfGetMetadataXML_f_ ) xfGetMetadataXML_f_ 
#endif
XMDF_API xid func_name(XFGETMETADATAXML_F,xfGetMetadataXML_f,xfGetMetadataXML_f )
              (xid *a_Id, char **a_xml, int *pathlength)
{
  int	error;

  error = xfGetMetadataXML(*a_Id, a_xml);

  return error;

} /* xfGetMetadataXML_f */

/* --------------------------------------------------------------------------- */
/* FUNCTION  xfGetMetadataProfileXML_f */
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFGETMETADATAPROFILEXML_F,xfGetMetadataProfileXML_f,xfGetMetadataProfileXML_f_ ) XFGETMETADATAPROFILEXML_F
#elif defined FLOWER
#define func_name(XFGETMETADATAPROFILEXML_F,xfGetMetadataProfileXML_f,xfGetMetadataProfileXML_f_ ) xfGetMetadataProfileXML_f
#else
#define func_name(XFGETMETADATAPROFILEXML_F,xfGetMetadataProfileXML_f,xfGetMetadataProfileXML_f_ ) xfGetMetadataProfileXML_f_
#endif
XMDF_API xid func_name(XFGETMETADATAPROFILEXML_F,xfGetMetadataProfileXML_f,xfGetMetadataProfileXML_f_ )
              (xid *a_Id, char **a_xml, int *pathlength)
{
  int	 error;

  error = xfGetMetadataProfileXML(*a_Id, a_xml);
  
  return error;

} /* xfGetMetadataProfileXML_f */

/* --------------------------------------------------------------------------- */
/* FUNCTION  xfGetMetadataSourceXML_f */
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFGETMETADATASOURCEXML_F,xfGetMetadataSourceXML_f,xfGetMetadataSourceXML_f_ ) XFGETMETADATASOURCEXML_F
#elif defined FLOWER
#define func_name(XFGETMETADATASOURCEXML_F,xfGetMetadataSourceXML_f,xfGetMetadataSourceXML_f_ ) xfGetMetadataSourceXML_f
#else
#define func_name(XFGETMETADATASOURCEXML_F,xfGetMetadataSourceXML_f,xfGetMetadataSourceXML_f_ ) xfGetMetadataSourceXML_f_
#endif
XMDF_API xid func_name(XFGETMETADATASOURCEXML_F,xfGetMetadataSourceXML_f,xfGetMetadataSourceXML_f_)
              (xid *a_Id, char **a_xml, int *pathlength)
{
  int	 error;

  error = xfGetMetadataSourceXML(*a_Id, a_xml);
  
  return error;

} /* xfGetMetadataSourceXML_f */

/* --------------------------------------------------------------------------- */
/* FUNCTION  xfGetMetadataSpatialXML_f */
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFGETMETADATASPATIALXML_F,xfGetMetadataSpatialXML_f,xfGetMetadataSpatialXML_f_) XFGETMETADATASPATIALXML_F
#elif defined FLOWER
#define func_name(XFGETMETADATASPATIALXML_F,xfGetMetadataSpatialXML_f,xfGetMetadataSpatialXML_f_) xfGetMetadataSpatialXML_f
#else
#define func_name(XFGETMETADATASPATIALXML_F,xfGetMetadataSpatialXML_f,xfGetMetadataSpatialXML_f_) xfGetMetadataSpatialXML_f_
#endif
XMDF_API xid func_name(XFGETMETADATASPATIALXML_F,xfGetMetadataSpatialXML_f,xfGetMetadataSpatialXML_f_)
              (xid *a_Id, char **a_xml, int *pathlength)
{
  int	 error;

  error = xfGetMetadataSpatialXML(*a_Id, a_xml);
  
  return error;

} /* xfGetMetadataSpatialXML_f */
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfsetnumberofxsects_f_*/
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFSETNUMBEROFXSECTS_F,xfsetnumberofxsects_f,xfsetnumberofxsects_f_) XFSETNUMBEROFXSECTS_F
#elif defined FLOWER
#define func_name(XFSETNUMBEROFXSECTS_F,xfsetnumberofxsects_f,xfsetnumberofxsects_f_) xfsetnumberofxsects_f
#else
#define func_name(XFSETNUMBEROFXSECTS_F,xfsetnumberofxsects_f,xfsetnumberofxsects_f_) xfsetnumberofxsects_f_
#endif
XMDF_API xid func_name(XFSETNUMBEROFXSECTS_F,xfsetnumberofxsects_f,xfsetnumberofxsects_f_)
              (xid *a_Id, int *a_nXSects, int *a_compression)
{
  int error;

  error = xfSetNumberOfXSects(*a_Id, a_nXSects, *a_compression);
  
  return error;

} /* xfsetnumberofxsects_f_*/

/* --------------------------------------------------------------------------- */
/* FUNCTION  xfgetnumberofxsects_f_*/
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFGETNUMBEROFXSECTS_F,xfgetnumberofxsects_f,xfgetnumberofxsects_f_) XFGETNUMBEROFXSECTS_F
#elif defined FLOWER
#define func_name(XFGETNUMBEROFXSECTS_F,xfgetnumberofxsects_f,xfgetnumberofxsects_f_) xfgetnumberofxsects_f
#else
#define func_name(XFGETNUMBEROFXSECTS_F,xfgetnumberofxsects_f,xfgetnumberofxsects_f_) xfgetnumberofxsects_f_
#endif
XMDF_API xid func_name(XFGETNUMBEROFXSECTS_F,xfgetnumberofxsects_f,xfgetnumberofxsects_f_)
              (xid *a_Id, int *a_nXSects)
{
  int error;

  error = xfGetNumberOfXSects(*a_Id, a_nXSects);
  
  return error;

} /* xfgetnumberofxsects_f_*/

/* --------------------------------------------------------------------------- */
/* FUNCTION  xfsetcsid_f_*/
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFSETCSID_F,xfsetcsid_f,xfsetcsid_f_) XFSETCSID_F
#elif defined FLOWER
#define func_name(XFSETCSID_F,xfsetcsid_f,xfsetcsid_f_) xfsetcsid_f
#else
#define func_name(XFSETCSID_F,xfsetcsid_f,xfsetcsid_f_) xfsetcsid_f_
#endif
XMDF_API xid func_name(XFSETCSID_F,xfsetcsid_f,xfsetcsid_f_)
              (xid *a_Id, int *a_NumVals, int *a_PropId,
                        int *a_compression)
{
  int error;

  error = xfSetCSID(*a_Id, *a_NumVals, a_PropId,
                        *a_compression);
  
  return error;

} /* xfsetcsid_f_*/

/* --------------------------------------------------------------------------- */
/* FUNCTION  xfgetcsid_f_*/
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFGETCSID_F,xfgetcsid_f,xfgetcsid_f_) XFGETCSID_F
#elif defined FLOWER
#define func_name(XFGETCSID_F,xfgetcsid_f,xfgetcsid_f_) xfgetcsid_f
#else
#define func_name(XFGETCSID_F,xfgetcsid_f,xfgetcsid_f_) xfgetcsid_f_
#endif
XMDF_API xid func_name(XFGETCSID_F,xfgetcsid_f,xfgetcsid_f_)
              (xid *a_Id, int *a_NumVals, int *a_PropId)
{
  int error;

  error = xfGetCSID(*a_Id, a_NumVals, a_PropId);
  
  return error;

} /* xfgetcsid_f_*/


/* --------------------------------------------------------------------------- */
/* FUNCTION  xfsetcsname_f_*/
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFSETCSNAME_F,xfsetcsname_f,xfsetcsname_f_) XFSETCSNAME_F
#elif defined FLOWER
#define func_name(XFSETCSNAME_F,xfsetcsname_f,xfsetcsname_f_) xfsetcsname_f
#else
#define func_name(XFSETCSNAME_F,xfsetcsname_f,xfsetcsname_f_) xfsetcsname_f_
#endif
XMDF_API xid func_name(XFSETCSNAME_F,xfsetcsname_f,xfsetcsname_f_)
              (xid *a_Id, int *a_NumVals, const char *a_PropId, int *a_StrLen)
{
  int	error;
  char *prop_copy = NULL;

  prop_copy = (char*)malloc((*a_StrLen+1)*sizeof(char));
  strncpy(prop_copy, a_PropId, *a_StrLen);
  prop_copy[*a_StrLen] = '\0';

  error = xfSetCSName(*a_Id, *a_NumVals, *a_StrLen, prop_copy);
  
  free(prop_copy);
  return error;

} /* xfsetcsname_f_*/

/* --------------------------------------------------------------------------- */
/* FUNCTION  xfgetcsname_f_*/
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFGETCSNAME_F,xfgetcsname_f,xfgetcsname_f_) XFGETCSNAME_F
#elif defined FLOWER
#define func_name(XFGETCSNAME_F,xfgetcsname_f,xfgetcsname_f_) xfgetcsname_f
#else
#define func_name(XFGETCSNAME_F,xfgetcsname_f,xfgetcsname_f_) xfgetcsname_f_
#endif
XMDF_API xid func_name(XFGETCSNAME_F,xfgetcsname_f,xfgetcsname_f_)
              (xid *a_Id, int *a_NumVals, char *a_PropId, int *a_StrLen)
{
  int	error;
  error = xfSetCSName(*a_Id, *a_NumVals, *a_StrLen, a_PropId);
  return error;

} /* xfgetcsname_f_*/


/* --------------------------------------------------------------------------- */
/* FUNCTION  xfgetcsnamelen_f_*/
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFGETCSNAMELEN_F,xfgetcsnamelen_f,xfgetcsnamelen_f_) XFGETCSNAMELEN_F
#elif defined FLOWER
#define func_name(XFGETCSNAMELEN_F,xfgetcsnamelen_f,xfgetcsnamelen_f_) xfgetcsnamelen_f
#else
#define func_name(XFGETCSNAMELEN_F,xfgetcsnamelen_f,xfgetcsnamelen_f_) xfgetcsnamelen_f_
#endif
XMDF_API xid func_name(XFGETCSNAMELEN_F,xfgetcsnamelen_f,xfgetcsnamelen_f_)
              (xid *a_Id, int *a_NumVals, int *a_StrLen)
{
  int error;

  error = xfGetCSNameLen(*a_Id, a_NumVals, a_StrLen);
  
  return error;

} /* xfgetcsnamelen_f_*/


/* --------------------------------------------------------------------------- */
/* FUNCTION  xfsetreachname_f_*/
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFSETREACHNAME_F,xfsetreachname_f,xfsetreachname_f_) XFSETREACHNAME_F
#elif defined FLOWER
#define func_name(XFSETREACHNAME_F,xfsetreachname_f,xfsetreachname_f_) xfsetreachname_f
#else
#define func_name(XFSETREACHNAME_F,xfsetreachname_f,xfsetreachname_f_) xfsetreachname_f_
#endif
XMDF_API xid func_name(XFSETREACHNAME_F,xfsetreachname_f,xfsetreachname_f_)
              (xid *a_Id, int *a_NumVals, const char *a_PropId, int *a_StrLen)
{
  int	error;
  char *prop_copy = NULL;

  prop_copy = (char*)malloc((*a_StrLen+1)*sizeof(char));
  strncpy(prop_copy, a_PropId, *a_StrLen);
  prop_copy[*a_StrLen] = '\0';

  error = xfSetReachName(*a_Id, *a_NumVals, *a_StrLen, prop_copy);
  
  free(prop_copy);
  return error;

} /* xfsetreachname_f_*/

 /* --------------------------------------------------------------------------- */
/* FUNCTION  xfgetreachname_f_*/
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFGETREACHNAME_F,xfgetreachname_f,xfgetreachname_f_) XFGETREACHNAME_F
#elif defined FLOWER
#define func_name(XFGETREACHNAME_F,xfgetreachname_f,xfgetreachname_f_) xfgetreachname_f
#else
#define func_name(XFGETREACHNAME_F,xfgetreachname_f,xfgetreachname_f_) xfgetreachname_f_
#endif
XMDF_API xid func_name(XFGETREACHNAME_F,xfgetreachname_f,xfgetreachname_f_)
              (xid *a_Id, int *a_NumVals, char *a_PropId, int *a_StrLen)
{
  int	error;
  error = xfGetReachName(*a_Id, a_NumVals, a_StrLen, a_PropId);
  return error;

} /* xfgetreachname_f_*/

/* --------------------------------------------------------------------------- */
/* FUNCTION  xfgetreachnamelen_f_*/
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFGETREACHNAMELEN_F,xfgetreachnamelen_f,xfgetreachnamelen_f_) XFGETREACHNAMELEN_F
#elif defined FLOWER
#define func_name(XFGETREACHNAMELEN_F,xfgetreachnamelen_f,xfgetreachnamelen_f_) xfgetreachnamelen_f
#else
#define func_name(XFGETREACHNAMELEN_F,xfgetreachnamelen_f,xfgetreachnamelen_f_) xfgetreachnamelen_f_
#endif
XMDF_API xid func_name(XFGETREACHNAMELEN_F,xfgetreachnamelen_f,xfgetreachnamelen_f_)
              (xid *a_Id, int *a_NumVals, int *a_StrLen)
{
  int error;

  error = xfGetReachNameLen(*a_Id, a_NumVals, a_StrLen);
  
  return error;

} /* xfgetreachnamelen_f_*/

/* --------------------------------------------------------------------------- */
/* FUNCTION  xfsettopoid_f_*/
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFSETTOPOID_F,xfsettopoid_f,xfsettopoid_f_) XFSETTOPOID_F
#elif defined FLOWER
#define func_name(XFSETTOPOID_F,xfsettopoid_f,xfsettopoid_f_) xfsettopoid_f
#else
#define func_name(XFSETTOPOID_F,xfsettopoid_f,xfsettopoid_f_) xfsettopoid_f_
#endif
XMDF_API xid func_name(XFSETTOPOID_F,xfsettopoid_f,xfsettopoid_f_)
              (xid *a_Id, int *a_NumVals, const char *a_PropId, int *a_StrLen)
{
  int	error;
  char *prop_copy = NULL;

  prop_copy = (char*)malloc((*a_StrLen+1)*sizeof(char));
  strncpy(prop_copy, a_PropId, *a_StrLen);
  prop_copy[*a_StrLen] = '\0';

  error = xfSetTopoID(*a_Id, *a_NumVals, *a_StrLen, prop_copy);

  free(prop_copy);
  return error;

} /* xfsettopoid_f_*/

/* --------------------------------------------------------------------------- */
/* FUNCTION  xfgettopoid_f_*/
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFGETTOPOID_F,xfgettopoid_f,xfgettopoid_f_) XFGETTOPOID_F
#elif defined FLOWER
#define func_name(XFGETTOPOID_F,xfgettopoid_f,xfgettopoid_f_) xfgettopoid_f
#else
#define func_name(XFGETTOPOID_F,xfgettopoid_f,xfgettopoid_f_) xfgettopoid_f_
#endif
XMDF_API xid func_name(XFGETTOPOID_F,xfgettopoid_f,xfgettopoid_f_)
              (xid *a_Id, int *a_NumVals, char *a_PropId, int *a_StrLen)
{
  int	error;
  error = xfGetTopoID(*a_Id, a_NumVals, a_StrLen, a_PropId);
  return error;

} /* xfgettopoid_f_*/

/* --------------------------------------------------------------------------- */
/* FUNCTION  xfgettopoidlen_f_*/
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFGETTOPOIDLEN_F,xfgettopoidlen_f,xfgettopoidlen_f_) XFGETTOPOIDLEN_F
#elif defined FLOWER
#define func_name(XFGETTOPOIDLEN_F,xfgettopoidlen_f,xfgettopoidlen_f_) xfgettopoidlen_f
#else
#define func_name(XFGETTOPOIDLEN_F,xfgettopoidlen_f,xfgettopoidlen_f_) xfgettopoidlen_f_
#endif
XMDF_API xid func_name(XFGETTOPOIDLEN_F,xfgettopoidlen_f,xfgettopoidlen_f_)
              (xid *a_Id, int *a_NumVals, int *a_StrLen)
{
  int error;

  error = xfGetTopoIDLen(*a_Id, a_NumVals, a_StrLen);
  
  return error;

} /* xfgettopoidlen_f_*

 */
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfsetstation_f_*/
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFSETSTATION_F,xfsetstation_f,xfsetstation_f_) XFSETSTATION_F
#elif defined FLOWER
#define func_name(XFSETSTATION_F,xfsetstation_f,xfsetstation_f_) xfsetstation_f
#else
#define func_name(XFSETSTATION_F,xfsetstation_f,xfsetstation_f_) xfsetstation_f_
#endif
XMDF_API xid func_name(XFSETSTATION_F,xfsetstation_f,xfsetstation_f_)
              (xid *a_Id, int *a_NumVals, double *a_PropId, int *a_compression)
{
  int error;

  error = xfSetStation(*a_Id, *a_NumVals, a_PropId, *a_compression);
  
  return error;

} /* xfsetstation_f_* */


/* --------------------------------------------------------------------------- */
/* FUNCTION  xfgetstation_f_*/
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFGETSTATION_F,xfgetstation_f,xfgetstation_f_) XFGETSTATION_F
#elif defined FLOWER
#define func_name(XFGETSTATION_F,xfgetstation_f,xfgetstation_f_) xfgetstation_f
#else
#define func_name(XFGETSTATION_F,xfgetstation_f,xfgetstation_f_) xfgetstation_f_
#endif
XMDF_API xid func_name(XFGETSTATION_F,xfgetstation_f,xfgetstation_f_)
              (xid *a_Id, int *a_NumVals, double *a_PropId)
{
  int error;

  error = xfGetStation(*a_Id, a_NumVals, a_PropId);
  
  return error;

} /* xfgetstation_f_* */


/* --------------------------------------------------------------------------- */
/* FUNCTION  xfsettype_f_*/
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFSETTYPE_F,xfsettype_f,xfsettype_f_) XFSETTYPE_F
#elif defined FLOWER
#define func_name(XFSETTYPE_F,xfsettype_f,xfsettype_f_) xfsettype_f
#else
#define func_name(XFSETTYPE_F,xfsettype_f,xfsettype_f_) xfsettype_f_
#endif
XMDF_API xid func_name(XFSETTYPE_F,xfsettype_f,xfsettype_f_)
              (xid *a_Id, int *a_NumVals, int *a_PropId, int *a_compression)
{
  int error;

  error = xfSetType(*a_Id, *a_NumVals, a_PropId, *a_compression);
  
  return error;

} /* xfsettype_f_*/


/* --------------------------------------------------------------------------- */
/* FUNCTION  xfgettype_f_*/
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFGETTYPE_F,xfgettype_f,xfgettype_f_) XFGETTYPE_F
#elif defined FLOWER
#define func_name(XFGETTYPE_F,xfgettype_f,xfgettype_f_) xfgettype_f
#else
#define func_name(XFGETTYPE_F,xfgettype_f,xfgettype_f_) xfgettype_f_
#endif
XMDF_API xid func_name(XFGETTYPE_F,xfgettype_f,xfgettype_f_)
              (xid *a_Id, int *a_NumVals, int *a_PropId)
{
  int error;

  error = xfGetType(*a_Id, a_NumVals, a_PropId);
  
  return error;

} /* xfgettype_f_*/


/* --------------------------------------------------------------------------- */
/* FUNCTION  xfsetptype_f_*/
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFSETPTYPE_F,xfsetptype_f,xfsetptype_f_) XFSETPTYPE_F
#elif defined FLOWER
#define func_name(XFSETPTYPE_F,xfsetptype_f,xfsetptype_f_) xfsetptype_f
#else
#define func_name(XFSETPTYPE_F,xfsetptype_f,xfsetptype_f_) xfsetptype_f_
#endif
XMDF_API xid func_name(XFSETPTYPE_F,xfsetptype_f,xfsetptype_f_)
              (xid *a_Id, int *a_NumVals, int *a_PropId, 
                         int *a_compression)
{
  int error;

  error = xfSetpType(*a_Id, *a_NumVals, a_PropId, *a_compression);
  
  return error;

} /* xfsetptype_f_*/

/* --------------------------------------------------------------------------- */
/* FUNCTION  xfgetptype_f_*/
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFGETPTYPE_F,xfgetptype_f,xfgetptype_f_) XFGETPTYPE_F
#elif defined FLOWER
#define func_name(XFGETPTYPE_F,xfgetptype_f,xfgetptype_f_) xfgetptype_f
#else
#define func_name(XFGETPTYPE_F,xfgetptype_f,xfgetptype_f_) xfgetptype_f_
#endif
XMDF_API xid func_name(XFGETPTYPE_F,xfgetptype_f,xfgetptype_f_)
              (xid *a_Id, int *a_NumVals, int *a_PropId)
{
  int error;

  error = xfGetpType(*a_Id, a_NumVals, a_PropId);
  
  return error;

} /* xfgetptype_f_*/


/* --------------------------------------------------------------------------- */
/* FUNCTION  xfsetpcsdblink_f_*/
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFSETPCSDBLINK_F,xfsetpcsdblink_f,xfsetpcsdblink_f_) XFSETPCSDBLINK_F
#elif defined FLOWER
#define func_name(XFSETPCSDBLINK_F,xfsetpcsdblink_f,xfsetpcsdblink_f_) xfsetpcsdblink_f
#else
#define func_name(XFSETPCSDBLINK_F,xfsetpcsdblink_f,xfsetpcsdblink_f_) xfsetpcsdblink_f_
#endif
XMDF_API xid func_name(XFSETPCSDBLINK_F,xfsetpcsdblink_f,xfsetpcsdblink_f_)
              (xid *a_Id, int *a_NumVals, int *a_PropId, int *a_compression)
{
  int error;

  error = xfSetpCSDBLink(*a_Id, *a_NumVals, a_PropId, *a_compression);
  
  return error;

} /* xfsetpcsdblink_f_*/

/* --------------------------------------------------------------------------- */
/* FUNCTION  xfgetpcsdblink_f_*/
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFGETPCSDBLINK_F,xfgetpcsdblink_f,xfgetpcsdblink_f_) XFGETPCSDBLINK_F
#elif defined FLOWER
#define func_name(XFGETPCSDBLINK_F,xfgetpcsdblink_f,xfgetpcsdblink_f_) xfgetpcsdblink_f
#else
#define func_name(XFGETPCSDBLINK_F,xfgetpcsdblink_f,xfgetpcsdblink_f_) xfgetpcsdblink_f_
#endif
XMDF_API xid func_name(XFGETPCSDBLINK_F,xfgetpcsdblink_f,xfgetpcsdblink_f_)
              (xid *a_Id, int *a_NumVals, int *a_PropId)
{
  int error;

  error = xfGetpCSDBLink(*a_Id, a_NumVals, a_PropId);
  
  return error;

} /* xfgetpcsdblink_f_*/


/* --------------------------------------------------------------------------- */
/* FUNCTION  xfsetnote_f_*/
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFSETNOTE_F,xfsetnote_f,xfsetnote_f_) XFSETNOTE_F
#elif defined FLOWER
#define func_name(XFSETNOTE_F,xfsetnote_f,xfsetnote_f_) xfsetnote_f
#else
#define func_name(XFSETNOTE_F,xfsetnote_f,xfsetnote_f_) xfsetnote_f_
#endif
XMDF_API xid func_name(XFSETNOTE_F,xfsetnote_f,xfsetnote_f_)
              (xid *a_Id, int *a_NumVals, const char *a_PropId, int *a_StrLen)
{
  int	error;
  char *prop_copy = NULL;

  prop_copy = (char*)malloc((*a_StrLen+1)*sizeof(char));
  strncpy(prop_copy, a_PropId, *a_StrLen);
  prop_copy[*a_StrLen] = '\0';

  error = xfSetNote(*a_Id, *a_NumVals, *a_StrLen, prop_copy);
  
  free(prop_copy);
  return error;

} /* xfsetnote_f_*/


/* --------------------------------------------------------------------------- */
/* FUNCTION  xfgetnote_f_*/
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFGETNOTE_F,xfgetnote_f,xfgetnote_f_) XFGETNOTE_F
#elif defined FLOWER
#define func_name(XFGETNOTE_F,xfgetnote_f,xfgetnote_f_) xfgetnote_f
#else
#define func_name(XFGETNOTE_F,xfgetnote_f,xfgetnote_f_) xfgetnote_f_
#endif
XMDF_API xid func_name(XFGETNOTE_F,xfgetnote_f,xfgetnote_f_)
              (xid *a_Id, int *a_NumVals, char *a_PropId, int *a_StrLen)
{
  int	error;
  error = xfGetNote(*a_Id, a_NumVals, a_StrLen, a_PropId);
  return error;

} /* xfgetnote_f_*/


/* --------------------------------------------------------------------------- */
/* FUNCTION  xfgetnotelen_f_*/
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFGETNOTELEN_F,xfgetnotelen_f,xfgetnotelen_f_) XFGETNOTELEN_F
#elif defined FLOWER
#define func_name(XFGETNOTELEN_F,xfgetnotelen_f,xfgetnotelen_f_) xfgetnotelen_f
#else
#define func_name(XFGETNOTELEN_F,xfgetnotelen_f,xfgetnotelen_f_) xfgetnotelen_f_
#endif
XMDF_API xid func_name(XFGETNOTELEN_F,xfgetnotelen_f,xfgetnotelen_f_)
              (xid *a_Id, int *a_NumVals, int *a_StrLen) 
{
  int error;

  error = xfGetNoteLen(*a_Id, a_NumVals, a_StrLen);
  
  return error;

} /* xfgetnotelen_f_*/


/* --------------------------------------------------------------------------- */
/* FUNCTION  xfsetxsectgeomx_f_*/
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFSETXSECTGEOMX_F,xfsetxsectgeomx_f,xfsetxsectgeomx_f_) XFSETXSECTGEOMX_F
#elif defined FLOWER
#define func_name(XFSETXSECTGEOMX_F,xfsetxsectgeomx_f,xfsetxsectgeomx_f_) xfsetxsectgeomx_f
#else
#define func_name(XFSETXSECTGEOMX_F,xfsetxsectgeomx_f,xfsetxsectgeomx_f_) xfsetxsectgeomx_f_
#endif
XMDF_API xid func_name(XFSETXSECTGEOMX_F,xfsetxsectgeomx_f,xfsetxsectgeomx_f_)
              (xid *a_Id, int *a_index, int *a_NumVals, 
                             double *a_iValues, int *a_compression)
{
  int error;

  error = xfSetXSectGeomX(*a_Id, *a_index, *a_NumVals, a_iValues, *a_compression);
  
  return error;

} /* xfsetxsectgeomx_f_*/

/* --------------------------------------------------------------------------- */
/* FUNCTION  xfgetxsectgeomx_f_*/
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFGETXSECTGEOMX_F,xfgetxsectgeomx_f,xfgetxsectgeomx_f_) XFGETXSECTGEOMX_F
#elif defined FLOWER
#define func_name(XFGETXSECTGEOMX_F,xfgetxsectgeomx_f,xfgetxsectgeomx_f_) xfgetxsectgeomx_f
#else
#define func_name(XFGETXSECTGEOMX_F,xfgetxsectgeomx_f,xfgetxsectgeomx_f_) xfgetxsectgeomx_f_
#endif
XMDF_API xid func_name(XFGETXSECTGEOMX_F,xfgetxsectgeomx_f,xfgetxsectgeomx_f_)
              (xid *a_Id, int *a_index, int *a_NumVals, 
                             double *a_iValues)
{
  int error;

  error = xfGetXSectGeomX(*a_Id, *a_index, a_NumVals, a_iValues);
  
  return error;

} /* xfgetxsectgeomx_f_*/


/* --------------------------------------------------------------------------- */
/* FUNCTION  xfsetxsectgeomy_f_*/
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFSETXSECTGEOMY_F,xfsetxsectgeomy_f,xfsetxsectgeomy_f_) XFSETXSECTGEOMY_F
#elif defined FLOWER
#define func_name(XFSETXSECTGEOMY_F,xfsetxsectgeomy_f,xfsetxsectgeomy_f_) xfsetxsectgeomy_f
#else
#define func_name(XFSETXSECTGEOMY_F,xfsetxsectgeomy_f,xfsetxsectgeomy_f_) xfsetxsectgeomy_f_
#endif
XMDF_API xid func_name(XFSETXSECTGEOMY_F,xfsetxsectgeomy_f,xfsetxsectgeomy_f_)
              (xid *a_Id, int *a_index, int *a_NumVals, 
                             double *a_iValues, int *a_compression)
{
  int error;

  error = xfSetXSectGeomY(*a_Id, *a_index, *a_NumVals, a_iValues, *a_compression);
  
  return error;

} /* xfsetxsectgeomy_f_*/


/* --------------------------------------------------------------------------- */
/* FUNCTION  xfgetxsectgeomy_f_*/
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFGETXSECTGEOMY_F,xfgetxsectgeomy_f,xfgetxsectgeomy_f_) XFGETXSECTGEOMY_F
#elif defined FLOWER
#define func_name(XFGETXSECTGEOMY_F,xfgetxsectgeomy_f,xfgetxsectgeomy_f_) xfgetxsectgeomy_f
#else
#define func_name(XFGETXSECTGEOMY_F,xfgetxsectgeomy_f,xfgetxsectgeomy_f_) xfgetxsectgeomy_f_
#endif
XMDF_API xid func_name(XFGETXSECTGEOMY_F,xfgetxsectgeomy_f,xfgetxsectgeomy_f_)
              (xid *a_Id, int *a_index, int *a_NumVals, 
                             double *a_iValues)
{
  int error;

  error = xfGetXSectGeomY(*a_Id, *a_index, a_NumVals, a_iValues);
  
  return error;

} /* xfgetxsectgeomy_f_*/


/* --------------------------------------------------------------------------- */
/* FUNCTION  xfsetxsectgeomd_f_*/
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFSETXSECTGEOMD_F,xfsetxsectgeomd_f,xfsetxsectgeomd_f_) XFSETXSECTGEOMD_F
#elif defined FLOWER
#define func_name(XFSETXSECTGEOMD_F,xfsetxsectgeomd_f,xfsetxsectgeomd_f_) xfsetxsectgeomd_f
#else
#define func_name(XFSETXSECTGEOMD_F,xfsetxsectgeomd_f,xfsetxsectgeomd_f_) xfsetxsectgeomd_f_
#endif
XMDF_API xid func_name(XFSETXSECTGEOMD_F,xfsetxsectgeomd_f,xfsetxsectgeomd_f_)
              (xid *a_Id, int *a_index, int *a_NumVals, 
                             double *a_iValues, int *a_compression)
{
  int error;

  error = xfSetXSectGeomD(*a_Id, *a_index, *a_NumVals, a_iValues, *a_compression);
  
  return error;

} /* xfsetxsectgeomd_f_*/


/* --------------------------------------------------------------------------- */
/* FUNCTION  xfgetxsectgeomd_f_*/
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFGETXSECTGEOMD_F,xfgetxsectgeomd_f,xfgetxsectgeomd_f_) XFGETXSECTGEOMD_F
#elif defined FLOWER
#define func_name(XFGETXSECTGEOMD_F,xfgetxsectgeomd_f,xfgetxsectgeomd_f_) xfgetxsectgeomd_f
#else
#define func_name(XFGETXSECTGEOMD_F,xfgetxsectgeomd_f,xfgetxsectgeomd_f_) xfgetxsectgeomd_f_
#endif
XMDF_API xid func_name(XFGETXSECTGEOMD_F,xfgetxsectgeomd_f,xfgetxsectgeomd_f_)
              (xid *a_Id, int *a_index, int *a_NumVals, 
                             double *a_iValues)
{
  int error;

  error = xfGetXSectGeomD(*a_Id, *a_index, a_NumVals, a_iValues);
  
  return error;

} /* xfgetxsectgeomd_f_*/


/* --------------------------------------------------------------------------- */
/* FUNCTION  xfsetxsectgeomz_f_*/
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFSETXSECTGEOMZ_F,xfsetxsectgeomz_f,xfsetxsectgeomz_f_) XFSETXSECTGEOMZ_F
#elif defined FLOWER
#define func_name(XFSETXSECTGEOMZ_F,xfsetxsectgeomz_f,xfsetxsectgeomz_f_) xfsetxsectgeomz_f
#else
#define func_name(XFSETXSECTGEOMZ_F,xfsetxsectgeomz_f,xfsetxsectgeomz_f_) xfsetxsectgeomz_f_
#endif
XMDF_API xid func_name(XFSETXSECTGEOMZ_F,xfsetxsectgeomz_f,xfsetxsectgeomz_f_)
              (xid *a_Id, int *a_index, int *a_NumVals, 
                             double *a_iValues, int *a_compression)
{
  int error;

  error = xfSetXSectGeomZ(*a_Id, *a_index, *a_NumVals, a_iValues, *a_compression);
  
  return error;

} /* xfsetxsectgeomz_f_*/


/* --------------------------------------------------------------------------- */
/* FUNCTION  xfgetxsectgeomz_f_*/
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFGETXSECTGEOMZ_F,xfgetxsectgeomz_f,xfgetxsectgeomz_f_) XFGETXSECTGEOMZ_F
#elif defined FLOWER
#define func_name(XFGETXSECTGEOMZ_F,xfgetxsectgeomz_f,xfgetxsectgeomz_f_) xfgetxsectgeomz_f
#else
#define func_name(XFGETXSECTGEOMZ_F,xfgetxsectgeomz_f,xfgetxsectgeomz_f_) xfgetxsectgeomz_f_
#endif
XMDF_API xid func_name(XFGETXSECTGEOMZ_F,xfgetxsectgeomz_f,xfgetxsectgeomz_f_)
              (xid *a_Id, int *a_index, int *a_NumVals, 
                             double *a_iValues)
{
  int error;

  error = xfGetXSectGeomZ(*a_Id, *a_index, a_NumVals, a_iValues);
  
  return error;

} /* xfgetxsectgeomz_f_*/


/*Line Properties*/

/* --------------------------------------------------------------------------- */
/* FUNCTION  xfsetxsectlinepropfrom_f_*/
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFSETXSECTLINEPROPFROM_F,xfsetxsectlinepropfrom_f,xfsetxsectlinepropfrom_f_) XFSETXSECTLINEPROPFROM_F
#elif defined FLOWER
#define func_name(XFSETXSECTLINEPROPFROM_F,xfsetxsectlinepropfrom_f,xfsetxsectlinepropfrom_f_) xfsetxsectlinepropfrom_f
#else
#define func_name(XFSETXSECTLINEPROPFROM_F,xfsetxsectlinepropfrom_f,xfsetxsectlinepropfrom_f_) xfsetxsectlinepropfrom_f_
#endif
XMDF_API xid func_name(XFSETXSECTLINEPROPFROM_F,xfsetxsectlinepropfrom_f,xfsetxsectlinepropfrom_f_)
              (xid *a_Id, int *a_index, int *a_NumVals, 
                                     double *a_PropId, int *a_compression)
{
  int error;

  error = xfSetXSectLinePropFrom(*a_Id, *a_index, *a_NumVals, a_PropId, *a_compression);
  
  return error;

} /* xfsetxsectlinepropfrom_f_*/


/* --------------------------------------------------------------------------- */
/* FUNCTION  xfgetxsectlinepropfrom_f_*/
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFGETXSECTLINEPROPFROM_F,xfgetxsectlinepropfrom_f,xfgetxsectlinepropfrom_f_) XFGETXSECTLINEPROPFROM_F
#elif defined FLOWER
#define func_name(XFGETXSECTLINEPROPFROM_F,xfgetxsectlinepropfrom_f,xfgetxsectlinepropfrom_f_) xfgetxsectlinepropfrom_f
#else
#define func_name(XFGETXSECTLINEPROPFROM_F,xfgetxsectlinepropfrom_f,xfgetxsectlinepropfrom_f_) xfgetxsectlinepropfrom_f_
#endif
XMDF_API xid func_name(XFGETXSECTLINEPROPFROM_F,xfgetxsectlinepropfrom_f,xfgetxsectlinepropfrom_f_)
              (xid *a_Id, int *a_index, int *a_NumVals, 
                                     double *a_PropId)
{
  int error;

  error = xfGetXSectLinePropFrom(*a_Id, *a_index, a_NumVals, a_PropId);
  
  return error;

} /* xfgetxsectlinepropfrom_f_*/


/* --------------------------------------------------------------------------- */
/* FUNCTION  xfsetxsectlinepropto_f_*/
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFSETXSECTLINEPROPTO_F,xfsetxsectlinepropto_f,xfsetxsectlinepropto_f_) XFSETXSECTLINEPROPTO_F
#elif defined FLOWER
#define func_name(XFSETXSECTLINEPROPTO_F,xfsetxsectlinepropto_f,xfsetxsectlinepropto_f_) xfsetxsectlinepropto_f
#else
#define func_name(XFSETXSECTLINEPROPTO_F,xfsetxsectlinepropto_f,xfsetxsectlinepropto_f_) xfsetxsectlinepropto_f_
#endif
XMDF_API xid func_name(XFSETXSECTLINEPROPTO_F,xfsetxsectlinepropto_f,xfsetxsectlinepropto_f_)
              (xid *a_Id, int *a_index, int *a_NumVals, 
                                   double *a_PropId, int *a_compression) 
{
  int error;

  error = xfSetXSectLinePropTo(*a_Id, *a_index, *a_NumVals, a_PropId, *a_compression);
  
  return error;

} /* xfsetxsectlinepropto_f_*/


/* --------------------------------------------------------------------------- */
/* FUNCTION  xfgetxsectlinepropto_f_*/
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFGETXSECTLINEPROPTO_F,xfgetxsectlinepropto_f,xfgetxsectlinepropto_f_) XFGETXSECTLINEPROPTO_F
#elif defined FLOWER
#define func_name(XFGETXSECTLINEPROPTO_F,xfgetxsectlinepropto_f,xfgetxsectlinepropto_f_) xfgetxsectlinepropto_f
#else
#define func_name(XFGETXSECTLINEPROPTO_F,xfgetxsectlinepropto_f,xfgetxsectlinepropto_f_) xfgetxsectlinepropto_f_
#endif
XMDF_API xid func_name(XFGETXSECTLINEPROPTO_F,xfgetxsectlinepropto_f,xfgetxsectlinepropto_f_)
              (xid *a_Id, int *a_index, int *a_NumVals, 
                                   double *a_PropId, int *a_compression) 
{
  int error;

  error = xfGetXSectLinePropTo(*a_Id, *a_index, a_NumVals, a_PropId);
  
  return error;

} /* xfgetxsectlinepropto_f_*/


/* --------------------------------------------------------------------------- */
/* FUNCTION  xfsetxsectlineproptype_f_*/
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFSETXSECTLINEPROPTYPE_F,xfsetxsectlineproptype_f,xfsetxsectlineproptype_f_) XFSETXSECTLINEPROPTYPE_F
#elif defined FLOWER
#define func_name(XFSETXSECTLINEPROPTYPE_F,xfsetxsectlineproptype_f,xfsetxsectlineproptype_f_) xfsetxsectlineproptype_f
#else
#define func_name(XFSETXSECTLINEPROPTYPE_F,xfsetxsectlineproptype_f,xfsetxsectlineproptype_f_) xfsetxsectlineproptype_f_
#endif
XMDF_API xid func_name(XFSETXSECTLINEPROPTYPE_F,xfsetxsectlineproptype_f,xfsetxsectlineproptype_f_)
              (xid *a_Id, int *a_index, int *a_NumVals, 
                                     int *a_PropId, int *a_compression)
{
  int error;

  error = xfSetXSectLinePropType(*a_Id, *a_index, *a_NumVals, a_PropId, *a_compression);
  
  return error;

} /* xfsetxsectlineproptype_f_*/


/* --------------------------------------------------------------------------- */
/* FUNCTION  xfgetxsectlineproptype_f_*/
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFGETXSECTLINEPROPTYPE_F,xfgetxsectlineproptype_f,xfgetxsectlineproptype_f_) XFGETXSECTLINEPROPTYPE_F
#elif defined FLOWER
#define func_name(XFGETXSECTLINEPROPTYPE_F,xfgetxsectlineproptype_f,xfgetxsectlineproptype_f_) xfgetxsectlineproptype_f
#else
#define func_name(XFGETXSECTLINEPROPTYPE_F,xfgetxsectlineproptype_f,xfgetxsectlineproptype_f_) xfgetxsectlineproptype_f_
#endif
XMDF_API xid func_name(XFGETXSECTLINEPROPTYPE_F,xfgetxsectlineproptype_f,xfgetxsectlineproptype_f_)
              (xid *a_Id, int *a_index, int *a_NumVals, 
                                     int *a_PropId)
{
  int error;

  error = xfGetXSectLinePropType(*a_Id, *a_index, a_NumVals, a_PropId);
  
  return error;

} /* xfgetxsectlineproptype_f_*/


/* --------------------------------------------------------------------------- */
/* FUNCTION  xfsetxsectlinepropivalue_f_*/
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFSETXSECTLINEPROPIVALUE_F,xfsetxsectlinepropivalue_f,xfsetxsectlinepropivalue_f_) XFSETXSECTLINEPROPIVALUE_F
#elif defined FLOWER
#define func_name(XFSETXSECTLINEPROPIVALUE_F,xfsetxsectlinepropivalue_f,xfsetxsectlinepropivalue_f_) xfsetxsectlinepropivalue_f
#else
#define func_name(XFSETXSECTLINEPROPIVALUE_F,xfsetxsectlinepropivalue_f,xfsetxsectlinepropivalue_f_) xfsetxsectlinepropivalue_f_
#endif
XMDF_API xid func_name(XFSETXSECTLINEPROPIVALUE_F,xfsetxsectlinepropivalue_f,xfsetxsectlinepropivalue_f_)
              (xid *a_Id, int *a_index, int *a_NumVals, 
                                       int *a_PropId, int *a_compression)
{
  int error;

  error = xfSetXSectLinePropIValue(*a_Id, *a_index, *a_NumVals, a_PropId, *a_compression);
  
  return error;

} /* xfsetxsectlinepropivalue_f_*/


/* --------------------------------------------------------------------------- */
/* FUNCTION  xfgetxsectlinepropivalue_f_*/
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFGETXSECTLINEPROPIVALUE_F,xfgetxsectlinepropivalue_f,xfgetxsectlinepropivalue_f_) XFGETXSECTLINEPROPIVALUE_F
#elif defined FLOWER
#define func_name(XFGETXSECTLINEPROPIVALUE_F,xfgetxsectlinepropivalue_f,xfgetxsectlinepropivalue_f_) xfgetxsectlinepropivalue_f
#else
#define func_name(XFGETXSECTLINEPROPIVALUE_F,xfgetxsectlinepropivalue_f,xfgetxsectlinepropivalue_f_) xfgetxsectlinepropivalue_f_
#endif
XMDF_API xid func_name(XFGETXSECTLINEPROPIVALUE_F,xfgetxsectlinepropivalue_f,xfgetxsectlinepropivalue_f_)
              (xid *a_Id, int *a_index, int *a_NumVals, 
                                     int *a_PropId)
{
  int error;

  error = xfGetXSectLinePropIValue(*a_Id, *a_index, a_NumVals, a_PropId);
  
  return error;

} /* xfgetxsectlinepropivalue_f_*/


/* --------------------------------------------------------------------------- */
/* FUNCTION  xfsetxsectlinepropfvalue_f_*/
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFSETXSECTLINEPROPFVALUE_F,xfsetxsectlinepropfvalue_f,xfsetxsectlinepropfvalue_f_) XFSETXSECTLINEPROPFVALUE_F
#elif defined FLOWER
#define func_name(XFSETXSECTLINEPROPFVALUE_F,xfsetxsectlinepropfvalue_f,xfsetxsectlinepropfvalue_f_) xfsetxsectlinepropfvalue_f
#else
#define func_name(XFSETXSECTLINEPROPFVALUE_F,xfsetxsectlinepropfvalue_f,xfsetxsectlinepropfvalue_f_) xfsetxsectlinepropfvalue_f_
#endif
XMDF_API xid func_name(XFSETXSECTLINEPROPFVALUE_F,xfsetxsectlinepropfvalue_f,xfsetxsectlinepropfvalue_f_)
              (xid *a_Id, int *a_index, int *a_NumVals, 
                                       double *a_PropId, int *a_compression)
{
  int error;

  error = xfSetXSectLinePropFValue(*a_Id, *a_index, *a_NumVals, a_PropId, *a_compression);
  
  return error;

} /* xfsetxsectlinepropfvalue_f_*/


/* --------------------------------------------------------------------------- */
/* FUNCTION  xfgetxsectlinepropfvalue_f_*/
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFGETXSECTLINEPROPFVALUE_F,xfgetxsectlinepropfvalue_f,xfgetxsectlinepropfvalue_f_) XFGETXSECTLINEPROPFVALUE_F
#elif defined FLOWER
#define func_name(XFGETXSECTLINEPROPFVALUE_F,xfgetxsectlinepropfvalue_f,xfgetxsectlinepropfvalue_f_) xfgetxsectlinepropfvalue_f
#else
#define func_name(XFGETXSECTLINEPROPFVALUE_F,xfgetxsectlinepropfvalue_f,xfgetxsectlinepropfvalue_f_) xfgetxsectlinepropfvalue_f_
#endif
XMDF_API xid func_name(XFGETXSECTLINEPROPFVALUE_F,xfgetxsectlinepropfvalue_f,xfgetxsectlinepropfvalue_f_)
              (xid *a_Id, int *a_index, int *a_NumVals, 
                                       double *a_PropId)
{
  int error;

  error = xfGetXSectLinePropFValue(*a_Id, *a_index, a_NumVals, a_PropId);
  
  return error;

} /* xfgetxsectlinepropfvalue_f_*/


/* --------------------------------------------------------------------------- */
/* FUNCTION  xfsetxsectlinepropid_f_*/
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFSETXSECTLINEPROPID_F,xfsetxsectlinepropid_f,xfsetxsectlinepropid_f_) XFSETXSECTLINEPROPID_F
#elif defined FLOWER
#define func_name(XFSETXSECTLINEPROPID_F,xfsetxsectlinepropid_f,xfsetxsectlinepropid_f_) xfsetxsectlinepropid_f
#else
#define func_name(XFSETXSECTLINEPROPID_F,xfsetxsectlinepropid_f,xfsetxsectlinepropid_f_) xfsetxsectlinepropid_f_
#endif
XMDF_API xid func_name(XFSETXSECTLINEPROPID_F,xfsetxsectlinepropid_f,xfsetxsectlinepropid_f_)
              (xid *a_Id, int *a_NumVals, 
                                   int *a_PropId, int *a_compression)
{
  int error;

  error = xfSetXSectLinePropID(*a_Id, *a_NumVals, a_PropId, *a_compression);
  
  return error;

} /* xfsetxsectlinepropid_f_*/


/* --------------------------------------------------------------------------- */
/* FUNCTION  xfgetxsectlinepropid_f_*/
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFGETXSECTLINEPROPID_F,xfgetxsectlinepropid_f,xfgetxsectlinepropid_f_) XFGETXSECTLINEPROPID_F
#elif defined FLOWER
#define func_name(XFGETXSECTLINEPROPID_F,xfgetxsectlinepropid_f,xfgetxsectlinepropid_f_) xfgetxsectlinepropid_f
#else
#define func_name(XFGETXSECTLINEPROPID_F,xfgetxsectlinepropid_f,xfgetxsectlinepropid_f_) xfgetxsectlinepropid_f_
#endif
XMDF_API xid func_name(XFGETXSECTLINEPROPID_F,xfgetxsectlinepropid_f,xfgetxsectlinepropid_f_)
              (xid *a_Id, int *a_NumVals, 
                                   int *a_PropId)
{
  int error;

  error = xfGetXSectLinePropID(*a_Id, a_NumVals, a_PropId);
  
  return error;

} /* xfgetxsectlinepropid_f_*/


/* --------------------------------------------------------------------------- */
/* FUNCTION  xfsetxsectlinepropname_f_*/
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFSETXSECTLINEPROPNAME_F,xfsetxsectlinepropname_f,xfsetxsectlinepropname_f_) XFSETXSECTLINEPROPNAME_F
#elif defined FLOWER
#define func_name(XFSETXSECTLINEPROPNAME_F,xfsetxsectlinepropname_f,xfsetxsectlinepropname_f_) xfsetxsectlinepropname_f
#else
#define func_name(XFSETXSECTLINEPROPNAME_F,xfsetxsectlinepropname_f,xfsetxsectlinepropname_f_) xfsetxsectlinepropname_f_
#endif
XMDF_API xid func_name(XFSETXSECTLINEPROPNAME_F,xfsetxsectlinepropname_f,xfsetxsectlinepropname_f_)
              (xid *a_Id, int *a_NumVals, const char *a_PropId, int *a_StrLen )
{
  int   error;
  char *prop_copy = NULL;

  prop_copy = (char*)malloc((*a_StrLen+1)*sizeof(char));
  strncpy(prop_copy, a_PropId, *a_StrLen);
  prop_copy[*a_StrLen] = '\0';

  error = xfSetXSectLinePropName(*a_Id, *a_NumVals, *a_StrLen, prop_copy);
  
  free(prop_copy);
  return error;

} /* xfsetxsectlinepropname_f_*/


/* --------------------------------------------------------------------------- */
/* FUNCTION  xfgetxsectlinepropname_f_*/
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFGETXSECTLINEPROPNAME_F,xfgetxsectlinepropname_f,xfgetxsectlinepropname_f_) XFGETXSECTLINEPROPNAME_F
#elif defined FLOWER
#define func_name(XFGETXSECTLINEPROPNAME_F,xfgetxsectlinepropname_f,xfgetxsectlinepropname_f_) xfgetxsectlinepropname_f
#else
#define func_name(XFGETXSECTLINEPROPNAME_F,xfgetxsectlinepropname_f,xfgetxsectlinepropname_f_) xfgetxsectlinepropname_f_
#endif
XMDF_API xid func_name(XFGETXSECTLINEPROPNAME_F,xfgetxsectlinepropname_f,xfgetxsectlinepropname_f_)
              (xid *a_Id, int *a_NumVals, char *a_PropId, int *a_StrLen)
{
  int	error;
  error = xfGetXSectLinePropName(*a_Id, a_NumVals, a_StrLen, a_PropId);
  return error;

} /* xfgetxsectlinepropname_f_*/


/* --------------------------------------------------------------------------- */
/* FUNCTION  xfgetxsectlinepropnamelen_f_*/
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFGETXSECTLINEPROPNAMELEN_F,xfgetxsectlinepropnamelen_f,xfgetxsectlinepropnamelen_f_) XFGETXSECTLINEPROPNAMELEN_F
#elif defined FLOWER
#define func_name(XFGETXSECTLINEPROPNAMELEN_F,xfgetxsectlinepropnamelen_f,xfgetxsectlinepropnamelen_f_) xfgetxsectlinepropnamelen_f
#else
#define func_name(XFGETXSECTLINEPROPNAMELEN_F,xfgetxsectlinepropnamelen_f,xfgetxsectlinepropnamelen_f_) xfgetxsectlinepropnamelen_f_
#endif
XMDF_API xid func_name(XFGETXSECTLINEPROPNAMELEN_F,xfgetxsectlinepropnamelen_f,xfgetxsectlinepropnamelen_f_)
              (xid *a_Id, int *a_NumVals, int *a_StrLen)
{
  int error;

  error = xfGetXSectLinePropNameLen(*a_Id, a_NumVals, a_StrLen);
  
  return error;

} /* xfgetxsectlinepropnamelen_f_*/


/* --------------------------------------------------------------------------- */
/* FUNCTION  xfsetxsectlinepropdesc_f_*/
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFSETXSECTLINEPROPDESC_F,xfsetxsectlinepropdesc_f,xfsetxsectlinepropdesc_f_) XFSETXSECTLINEPROPDESC_F
#elif defined FLOWER
#define func_name(XFSETXSECTLINEPROPDESC_F,xfsetxsectlinepropdesc_f,xfsetxsectlinepropdesc_f_) xfsetxsectlinepropdesc_f
#else
#define func_name(XFSETXSECTLINEPROPDESC_F,xfsetxsectlinepropdesc_f,xfsetxsectlinepropdesc_f_) xfsetxsectlinepropdesc_f_
#endif
XMDF_API xid func_name(XFSETXSECTLINEPROPDESC_F,xfsetxsectlinepropdesc_f,xfsetxsectlinepropdesc_f_)
              (xid *a_Id, int *a_NumVals, const char *a_PropId, int *a_StrLen)
{
  int	error;
  char *prop_copy = NULL;

  prop_copy = (char*)malloc((*a_StrLen+1)*sizeof(char));
  strncpy(prop_copy, a_PropId, *a_StrLen);
  prop_copy[*a_StrLen] = '\0';

  error = xfSetXSectLinePropDesc(*a_Id, *a_NumVals, *a_StrLen, prop_copy);
  
  free(prop_copy);
  return error;

} /* xfsetxsectlinepropdesc_f_*/


/* --------------------------------------------------------------------------- */
/* FUNCTION  xfgetxsectlinepropdesc_f_*/
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFGETXSECTLINEPROPDESC_F,xfgetxsectlinepropdesc_f,xfgetxsectlinepropdesc_f_) XFGETXSECTLINEPROPDESC_F
#elif defined FLOWER
#define func_name(XFGETXSECTLINEPROPDESC_F,xfgetxsectlinepropdesc_f,xfgetxsectlinepropdesc_f_) xfgetxsectlinepropdesc_f
#else
#define func_name(XFGETXSECTLINEPROPDESC_F,xfgetxsectlinepropdesc_f,xfgetxsectlinepropdesc_f_) xfgetxsectlinepropdesc_f_
#endif
XMDF_API xid func_name(XFGETXSECTLINEPROPDESC_F,xfgetxsectlinepropdesc_f,xfgetxsectlinepropdesc_f_)
              (xid *a_Id, int *a_NumVals, char *a_PropId, int *a_StrLen)
{
  int	error;
  error = xfGetXSectLinePropDesc(*a_Id, a_NumVals, a_StrLen, a_PropId);
  return error;

} /* xfgetxsectlinepropdesc_f_*/

 
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfgetxsectlinepropdesclen_f_*/
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFGETXSECTLINEPROPDESCLEN_F,xfgetxsectlinepropdesclen_f,xfgetxsectlinepropdesclen_f_) XFGETXSECTLINEPROPDESCLEN_F
#elif defined FLOWER
#define func_name(XFGETXSECTLINEPROPDESCLEN_F,xfgetxsectlinepropdesclen_f,xfgetxsectlinepropdesclen_f_) xfgetxsectlinepropdesclen_f
#else
#define func_name(XFGETXSECTLINEPROPDESCLEN_F,xfgetxsectlinepropdesclen_f,xfgetxsectlinepropdesclen_f_) xfgetxsectlinepropdesclen_f_
#endif
XMDF_API xid func_name(XFGETXSECTLINEPROPDESCLEN_F,xfgetxsectlinepropdesclen_f,xfgetxsectlinepropdesclen_f_)
              (xid *a_Id, int *a_NumVals, int *a_StrLen)
{
  int error;

  error = xfGetXSectLinePropDescLen(*a_Id, a_NumVals, a_StrLen);
  
  return error;

} /* xfgetxsectlinepropdesclen_f_*/


/* --------------------------------------------------------------------------- */
/* FUNCTION  xfsetxsectlinepropcategory_f_*/
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFSETXSECTLINEPROPCATEGORY_F,xfsetxsectlinepropcategory_f,xfsetxsectlinepropcategory_f_) XFSETXSECTLINEPROPCATEGORY_F
#elif defined FLOWER
#define func_name(XFSETXSECTLINEPROPCATEGORY_F,xfsetxsectlinepropcategory_f,xfsetxsectlinepropcategory_f_) xfsetxsectlinepropcategory_f
#else
#define func_name(XFSETXSECTLINEPROPCATEGORY_F,xfsetxsectlinepropcategory_f,xfsetxsectlinepropcategory_f_) xfsetxsectlinepropcategory_f_
#endif
XMDF_API xid func_name(XFSETXSECTLINEPROPCATEGORY_F,xfsetxsectlinepropcategory_f,xfsetxsectlinepropcategory_f_)
              (xid *a_Id, int *a_NumVals, 
                                         int *a_PropId, int *a_compression)
{
  int error;

  error = xfSetXSectLinePropCategory(*a_Id, *a_NumVals, a_PropId, *a_compression);
  
  return error;

} /* xfsetxsectlinepropcategory_f_*/


/* --------------------------------------------------------------------------- */
/* FUNCTION  xfgetxsectlinepropcategory_f_*/
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFGETXSECTLINEPROPCATEGORY_F,xfgetxsectlinepropcategory_f,xfgetxsectlinepropcategory_f_) XFGETXSECTLINEPROPCATEGORY_F
#elif defined FLOWER
#define func_name(XFGETXSECTLINEPROPCATEGORY_F,xfgetxsectlinepropcategory_f,xfgetxsectlinepropcategory_f_) xfgetxsectlinepropcategory_f
#else
#define func_name(XFGETXSECTLINEPROPCATEGORY_F,xfgetxsectlinepropcategory_f,xfgetxsectlinepropcategory_f_) xfgetxsectlinepropcategory_f_
#endif
XMDF_API xid func_name(XFGETXSECTLINEPROPCATEGORY_F,xfgetxsectlinepropcategory_f,xfgetxsectlinepropcategory_f_)
              (xid *a_Id, int *a_NumVals, 
                                         int *a_PropId)
{
  int error;

  error = xfGetXSectLinePropCategory(*a_Id, a_NumVals, a_PropId);
  
  return error;

} /* xfgetxsectlinepropcategory_f_*/


/* --------------------------------------------------------------------------- */
/* FUNCTION  xfsetxsectlinepropfdefault_f_*/
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFSETXSECTLINEPROPFDEFAULT_F,xfsetxsectlinepropfdefault_f,xfsetxsectlinepropfdefault_f_) XFSETXSECTLINEPROPFDEFAULT_F
#elif defined FLOWER
#define func_name(XFSETXSECTLINEPROPFDEFAULT_F,xfsetxsectlinepropfdefault_f,xfsetxsectlinepropfdefault_f_) xfsetxsectlinepropfdefault_f
#else
#define func_name(XFSETXSECTLINEPROPFDEFAULT_F,xfsetxsectlinepropfdefault_f,xfsetxsectlinepropfdefault_f_) xfsetxsectlinepropfdefault_f_
#endif
XMDF_API xid func_name(XFSETXSECTLINEPROPFDEFAULT_F,xfsetxsectlinepropfdefault_f,xfsetxsectlinepropfdefault_f_)
              (xid *a_Id, int *a_NumVals, 
                                         double *a_PropId, int *a_compression)
{
  int error;

  error = xfSetXSectLinePropFDefault(*a_Id, *a_NumVals, a_PropId, *a_compression);
  
  return error;

} /* xfsetxsectlinepropfdefault_f_*/


/* --------------------------------------------------------------------------- */
/* FUNCTION  xfgetxsectlinepropfdefault_f_*/
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFGETXSECTLINEPROPFDEFAULT_F,xfgetxsectlinepropfdefault_f,xfgetxsectlinepropfdefault_f_) XFGETXSECTLINEPROPFDEFAULT_F
#elif defined FLOWER
#define func_name(XFGETXSECTLINEPROPFDEFAULT_F,xfgetxsectlinepropfdefault_f,xfgetxsectlinepropfdefault_f_) xfgetxsectlinepropfdefault_f
#else
#define func_name(XFGETXSECTLINEPROPFDEFAULT_F,xfgetxsectlinepropfdefault_f,xfgetxsectlinepropfdefault_f_) xfgetxsectlinepropfdefault_f_
#endif
XMDF_API xid func_name(XFGETXSECTLINEPROPFDEFAULT_F,xfgetxsectlinepropfdefault_f,xfgetxsectlinepropfdefault_f_)
              (xid *a_Id, int *a_NumVals, 
                                         double *a_PropId)
{
  int error;

  error = xfGetXSectLinePropFDefault(*a_Id, a_NumVals, a_PropId);
  
  return error;

} /* xfgetxsectlinepropfdefault_f_*/


/* --------------------------------------------------------------------------- */
/* FUNCTION  xfsetxsectlinepropexclusive_f_*/
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFSETXSECTLINEPROPEXCLUSIVE_F,xfsetxsectlinepropexclusive_f,xfsetxsectlinepropexclusive_f_) XFSETXSECTLINEPROPEXCLUSIVE_F
#elif defined FLOWER
#define func_name(XFSETXSECTLINEPROPEXCLUSIVE_F,xfsetxsectlinepropexclusive_f,xfsetxsectlinepropexclusive_f_) xfsetxsectlinepropexclusive_f
#else
#define func_name(XFSETXSECTLINEPROPEXCLUSIVE_F,xfsetxsectlinepropexclusive_f,xfsetxsectlinepropexclusive_f_) xfsetxsectlinepropexclusive_f_
#endif
XMDF_API xid func_name(XFSETXSECTLINEPROPEXCLUSIVE_F,xfsetxsectlinepropexclusive_f,xfsetxsectlinepropexclusive_f_)
              (xid *a_Id, int *a_NumVals, 
                                          int *a_PropId, int *a_compression)
{
  int error;

  error = xfSetXSectLinePropExclusive(*a_Id, *a_NumVals, a_PropId, *a_compression);
  
  return error;

} /* xfsetxsectlinepropexclusive_f_*/


/* --------------------------------------------------------------------------- */
/* FUNCTION  xfgetxsectlinepropexclusive_f_*/
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFGETXSECTLINEPROPEXCLUSIVE_F,xfgetxsectlinepropexclusive_f,xfgetxsectlinepropexclusive_f_) XFGETXSECTLINEPROPEXCLUSIVE_F
#elif defined FLOWER
#define func_name(XFGETXSECTLINEPROPEXCLUSIVE_F,xfgetxsectlinepropexclusive_f,xfgetxsectlinepropexclusive_f_) xfgetxsectlinepropexclusive_f
#else
#define func_name(XFGETXSECTLINEPROPEXCLUSIVE_F,xfgetxsectlinepropexclusive_f,xfgetxsectlinepropexclusive_f_) xfgetxsectlinepropexclusive_f_
#endif
XMDF_API xid func_name(XFGETXSECTLINEPROPEXCLUSIVE_F,xfgetxsectlinepropexclusive_f,xfgetxsectlinepropexclusive_f_)
              (xid *a_Id, int *a_NumVals, 
                                          int *a_PropId)
{
  int error;

  error = xfGetXSectLinePropExclusive(*a_Id, a_NumVals, a_PropId);
  
  return error;

} /* xfgetxsectlinepropexclusive_f_*/


/* --------------------------------------------------------------------------- */
/* FUNCTION  xfsetnumberoflnpropenumgroup_f_*/
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFSETNUMBEROFLNPROPENUMGROUP_F,xfsetnumberoflnpropenumgroup_f,xfsetnumberoflnpropenumgroup_f_) XFSETNUMBEROFLNPROPENUMGROUP_F
#elif defined FLOWER
#define func_name(XFSETNUMBEROFLNPROPENUMGROUP_F,xfsetnumberoflnpropenumgroup_f,xfsetnumberoflnpropenumgroup_f_) xfsetnumberoflnpropenumgroup_f
#else
#define func_name(XFSETNUMBEROFLNPROPENUMGROUP_F,xfsetnumberoflnpropenumgroup_f,xfsetnumberoflnpropenumgroup_f_) xfsetnumberoflnpropenumgroup_f_
#endif
XMDF_API xid func_name(XFSETNUMBEROFLNPROPENUMGROUP_F,xfsetnumberoflnpropenumgroup_f,xfsetnumberoflnpropenumgroup_f_)
              (xid *a_Id, int *a_nPropNum, 
                                        int *a_compression)
{
  int error;

  error = xfSetNumberOfLinePropEnumGroup(*a_Id, a_nPropNum, *a_compression);
  
  return error;

} /* xfsetnumberoflnpropenumgroup_f_*/


/* --------------------------------------------------------------------------- */
/* FUNCTION  xfgetnumberoflnpropenumgroup_f_*/
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFGETNUMBEROFLNPROPENUMGROUP_F,xfgetnumberoflnpropenumgroup_f,xfgetnumberoflnpropenumgroup_f_) XFGETNUMBEROFLNPROPENUMGROUP_F
#elif defined FLOWER
#define func_name(XFGETNUMBEROFLNPROPENUMGROUP_F,xfgetnumberoflnpropenumgroup_f,xfgetnumberoflnpropenumgroup_f_) xfgetnumberoflnpropenumgroup_f
#else
#define func_name(XFGETNUMBEROFLNPROPENUMGROUP_F,xfgetnumberoflnpropenumgroup_f,xfgetnumberoflnpropenumgroup_f_) xfgetnumberoflnpropenumgroup_f_
#endif
XMDF_API xid func_name(XFGETNUMBEROFLNPROPENUMGROUP_F,xfgetnumberoflnpropenumgroup_f,xfgetnumberoflnpropenumgroup_f_)
              (xid *a_Id, int *a_nXSects)
{
  int error;

  error = xfGetNumberOfLinePropEnumGroup(*a_Id, a_nXSects);
  
  return error;

} /* xfgetnumberoflnpropenumgroup_f_*/


/* --------------------------------------------------------------------------- */
/* FUNCTION  xfsetnumberoflinepropenum_f_*/
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFSETNUMBEROFLINEPROPENUM_F,xfsetnumberoflinepropenum_f,xfsetnumberoflinepropenum_f_) XFSETNUMBEROFLINEPROPENUM_F
#elif defined FLOWER
#define func_name(XFSETNUMBEROFLINEPROPENUM_F,xfsetnumberoflinepropenum_f,xfsetnumberoflinepropenum_f_) xfsetnumberoflinepropenum_f
#else
#define func_name(XFSETNUMBEROFLINEPROPENUM_F,xfsetnumberoflinepropenum_f,xfsetnumberoflinepropenum_f_) xfsetnumberoflinepropenum_f_
#endif
XMDF_API xid func_name(XFSETNUMBEROFLINEPROPENUM_F,xfsetnumberoflinepropenum_f,xfsetnumberoflinepropenum_f_)
              (xid *a_Id, int *a_index, int *a_nXSects, 
                                        int *a_compression)
{
  int error;

  error = xfSetNumberOfLinePropEnum(*a_Id, *a_index, a_nXSects, *a_compression);
  
  return error;

} /* xfsetnumberoflinepropenum_f_*/


/* --------------------------------------------------------------------------- */
/* FUNCTION  xfgetnumberoflinepropenum_f_*/
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFGETNUMBEROFLINEPROPENUM_F,xfgetnumberoflinepropenum_f,xfgetnumberoflinepropenum_f_) XFGETNUMBEROFLINEPROPENUM_F
#elif defined FLOWER
#define func_name(XFGETNUMBEROFLINEPROPENUM_F,xfgetnumberoflinepropenum_f,xfgetnumberoflinepropenum_f_) xfgetnumberoflinepropenum_f
#else
#define func_name(XFGETNUMBEROFLINEPROPENUM_F,xfgetnumberoflinepropenum_f,xfgetnumberoflinepropenum_f_) xfgetnumberoflinepropenum_f_
#endif
XMDF_API xid func_name(XFGETNUMBEROFLINEPROPENUM_F,xfgetnumberoflinepropenum_f,xfgetnumberoflinepropenum_f_)
              (xid *a_Id, int *a_index, int *a_nXSects)
{
  int error;

  error = xfGetNumberOfLinePropEnum(*a_Id, *a_index, a_nXSects);
  
  return error;

} /* xfgetnumberoflinepropenum_f_*/



/* --------------------------------------------------------------------------- */
/* FUNCTION  xfsetxsectlinepropenumid_f_*/
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFSETXSECTLINEPROPENUMID_F,xfsetxsectlinepropenumid_f,xfsetxsectlinepropenumid_f_) XFSETXSECTLINEPROPENUMID_F
#elif defined FLOWER
#define func_name(XFSETXSECTLINEPROPENUMID_F,xfsetxsectlinepropenumid_f,xfsetxsectlinepropenumid_f_) xfsetxsectlinepropenumid_f
#else
#define func_name(XFSETXSECTLINEPROPENUMID_F,xfsetxsectlinepropenumid_f,xfsetxsectlinepropenumid_f_) xfsetxsectlinepropenumid_f_
#endif
XMDF_API xid func_name(XFSETXSECTLINEPROPENUMID_F,xfsetxsectlinepropenumid_f,xfsetxsectlinepropenumid_f_)
              (xid *a_Id, int *a_index, int *a_NumVals, 
                                       int *a_PropId, int *a_compression)
{
  int error;

  error = xfSetXSectLinePropEnumID(*a_Id, *a_index, *a_NumVals, a_PropId, *a_compression);
  
  return error;

} /* xfsetxsectlinepropenumid_f_*/


/* --------------------------------------------------------------------------- */
/* FUNCTION  xfgetxsectlinepropenumid_f_*/
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFGETXSECTLINEPROPENUMID_F,xfgetxsectlinepropenumid_f,xfgetxsectlinepropenumid_f_) XFGETXSECTLINEPROPENUMID_F
#elif defined FLOWER
#define func_name(XFGETXSECTLINEPROPENUMID_F,xfgetxsectlinepropenumid_f,xfgetxsectlinepropenumid_f_) xfgetxsectlinepropenumid_f
#else
#define func_name(XFGETXSECTLINEPROPENUMID_F,xfgetxsectlinepropenumid_f,xfgetxsectlinepropenumid_f_) xfgetxsectlinepropenumid_f_
#endif
XMDF_API xid func_name(XFGETXSECTLINEPROPENUMID_F,xfgetxsectlinepropenumid_f,xfgetxsectlinepropenumid_f_)
              (xid *a_Id, int *a_index, int *a_NumVals, 
                                       int *a_PropId)
{
  int error;

  error = xfGetXSectLinePropEnumID(*a_Id, *a_index, a_NumVals, a_PropId);
  
  return error;

} /* xfgetxsectlinepropenumid_f_*/


/* --------------------------------------------------------------------------- */
/* FUNCTION  xfsetxsectlinepropenummatid_f_*/
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFSETXSECTLINEPROPENUMMATID_F,xfsetxsectlinepropenummatid_f,xfsetxsectlinepropenummatid_f_) XFSETXSECTLINEPROPENUMMATID_F
#elif defined FLOWER
#define func_name(XFSETXSECTLINEPROPENUMMATID_F,xfsetxsectlinepropenummatid_f,xfsetxsectlinepropenummatid_f_) xfsetxsectlinepropenummatid_f
#else
#define func_name(XFSETXSECTLINEPROPENUMMATID_F,xfsetxsectlinepropenummatid_f,xfsetxsectlinepropenummatid_f_) xfsetxsectlinepropenummatid_f_
#endif
XMDF_API xid func_name(XFSETXSECTLINEPROPENUMMATID_F,xfsetxsectlinepropenummatid_f,xfsetxsectlinepropenummatid_f_)
              (xid *a_Id, int *a_index, int *a_NumVals, 
                                          int *a_PropId, int *a_compression)
{
  int error;

  error = xfSetXSectLinePropEnumMatID(*a_Id, *a_index, *a_NumVals, a_PropId, *a_compression);
  
  return error;

} /* xfsetxsectlinepropenummatid_f_*/


/* --------------------------------------------------------------------------- */
/* FUNCTION  xfgetxsectlinepropenummatid_f_*/
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFGETXSECTLINEPROPENUMMATID_F,xfgetxsectlinepropenummatid_f,xfgetxsectlinepropenummatid_f_) XFGETXSECTLINEPROPENUMMATID_F
#elif defined FLOWER
#define func_name(XFGETXSECTLINEPROPENUMMATID_F,xfgetxsectlinepropenummatid_f,xfgetxsectlinepropenummatid_f_) xfgetxsectlinepropenummatid_f
#else
#define func_name(XFGETXSECTLINEPROPENUMMATID_F,xfgetxsectlinepropenummatid_f,xfgetxsectlinepropenummatid_f_) xfgetxsectlinepropenummatid_f_
#endif
XMDF_API xid func_name(XFGETXSECTLINEPROPENUMMATID_F,xfgetxsectlinepropenummatid_f,xfgetxsectlinepropenummatid_f_)
              (xid *a_Id, int *a_index, int *a_NumVals, 
                                          int *a_PropId)
{
  int error;

  error = xfGetXSectLinePropEnumMatID(*a_Id, *a_index, a_NumVals, a_PropId);
  
  return error;

} /* xfgetxsectlinepropenummatid_f_*/


/* --------------------------------------------------------------------------- */
/* FUNCTION  xfsetxsectlinepropenumname_f_*/
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFSETXSECTLINEPROPENUMNAME_F,xfsetxsectlinepropenumname_f,xfsetxsectlinepropenumname_f_) XFSETXSECTLINEPROPENUMNAME_F
#elif defined FLOWER
#define func_name(XFSETXSECTLINEPROPENUMNAME_F,xfsetxsectlinepropenumname_f,xfsetxsectlinepropenumname_f_) xfsetxsectlinepropenumname_f
#else
#define func_name(XFSETXSECTLINEPROPENUMNAME_F,xfsetxsectlinepropenumname_f,xfsetxsectlinepropenumname_f_) xfsetxsectlinepropenumname_f_
#endif
XMDF_API xid func_name(XFSETXSECTLINEPROPENUMNAME_F,xfsetxsectlinepropenumname_f,xfsetxsectlinepropenumname_f_)
              (xid *a_Id, int *a_index, int *a_NumVals, 
                                     const char *a_PropId, int *a_StrLen)
{
  int	error;
  char *prop_copy = NULL;

  prop_copy = (char*)malloc((*a_StrLen+1)*sizeof(char));
  strncpy(prop_copy, a_PropId, *a_StrLen);
  prop_copy[*a_StrLen] = '\0';

  error = xfSetXSectLinePropEnumName(*a_Id, *a_index, *a_NumVals, *a_StrLen, prop_copy);
  
  free(prop_copy);
  return error;

} /* xfsetxsectlinepropenumname_f_*/


/* --------------------------------------------------------------------------- */
/* FUNCTION  xfgetxsectlinepropenumname_f_*/
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFGETXSECTLINEPROPENUMNAME_F,xfgetxsectlinepropenumname_f,xfgetxsectlinepropenumname_f_) XFGETXSECTLINEPROPENUMNAME_F
#elif defined FLOWER
#define func_name(XFGETXSECTLINEPROPENUMNAME_F,xfgetxsectlinepropenumname_f,xfgetxsectlinepropenumname_f_) xfgetxsectlinepropenumname_f
#else
#define func_name(XFGETXSECTLINEPROPENUMNAME_F,xfgetxsectlinepropenumname_f,xfgetxsectlinepropenumname_f_) xfgetxsectlinepropenumname_f_
#endif
XMDF_API xid func_name(XFGETXSECTLINEPROPENUMNAME_F,xfgetxsectlinepropenumname_f,xfgetxsectlinepropenumname_f_)
              (xid *a_Id, int *a_index, int *a_NumVals, 
                                    char *a_PropId, int *a_StrLen)
{
  int	error;
  error = xfGetXSectLinePropEnumName(*a_Id, *a_index, a_NumVals, a_StrLen, a_PropId);
  return error;

} /* xfgetxsectlinepropenumname_f_*/


/* --------------------------------------------------------------------------- */
/* FUNCTION  xfgetxsectlnpropenumnamelen_f_*/
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFGETXSECTLNPROPENUMNAMELEN_F,xfgetxsectlnpropenumnamelen_f,xfgetxsectlnpropenumnamelen_f_) XFGETXSECTLNPROPENUMNAMELEN_F
#elif defined FLOWER
#define func_name(XFGETXSECTLNPROPENUMNAMELEN_F,xfgetxsectlnpropenumnamelen_f,xfgetxsectlnpropenumnamelen_f_) xfgetxsectlnpropenumnamelen_f
#else
#define func_name(XFGETXSECTLNPROPENUMNAMELEN_F,xfgetxsectlnpropenumnamelen_f,xfgetxsectlnpropenumnamelen_f_) xfgetxsectlnpropenumnamelen_f_
#endif
XMDF_API xid func_name(XFGETXSECTLNPROPENUMNAMELEN_F,xfgetxsectlnpropenumnamelen_f,xfgetxsectlnpropenumnamelen_f_)
              (xid *a_Id, int *a_index, int *a_NumVals,
                                            int *a_StrLen)
{
  int error;

  error = xfGetXSectLinePropEnumNameLen(*a_Id, *a_index, a_NumVals, a_StrLen);
  
  return error;

} /* xfgetxsectlnpropenumnamelen_f_*/


/* --------------------------------------------------------------------------- */
/* FUNCTION  xfsetxsectpointpropmeasure_f_*/
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFSETXSECTPOINTPROPMEASURE_F,xfsetxsectpointpropmeasure_f,xfsetxsectpointpropmeasure_f_) XFSETXSECTPOINTPROPMEASURE_F
#elif defined FLOWER
#define func_name(XFSETXSECTPOINTPROPMEASURE_F,xfsetxsectpointpropmeasure_f,xfsetxsectpointpropmeasure_f_) xfsetxsectpointpropmeasure_f
#else
#define func_name(XFSETXSECTPOINTPROPMEASURE_F,xfsetxsectpointpropmeasure_f,xfsetxsectpointpropmeasure_f_) xfsetxsectpointpropmeasure_f_
#endif
XMDF_API xid func_name(XFSETXSECTPOINTPROPMEASURE_F,xfsetxsectpointpropmeasure_f,xfsetxsectpointpropmeasure_f_)
              (xid *a_Id, int *a_index, int *a_NumVals, 
                                         double *a_PropId, int *a_compression)
{
  int error;

  error = xfSetXSectPointPropMeasure(*a_Id, *a_index, *a_NumVals, a_PropId, *a_compression);
  
  return error;

} /* xfsetxsectpointpropmeasure_f_*/


/* --------------------------------------------------------------------------- */
/* FUNCTION  xfgetxsectpointpropmeasure_f_*/
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFGETXSECTPOINTPROPMEASURE_F,xfgetxsectpointpropmeasure_f,xfgetxsectpointpropmeasure_f_) XFGETXSECTPOINTPROPMEASURE_F
#elif defined FLOWER
#define func_name(XFGETXSECTPOINTPROPMEASURE_F,xfgetxsectpointpropmeasure_f,xfgetxsectpointpropmeasure_f_) xfgetxsectpointpropmeasure_f
#else
#define func_name(XFGETXSECTPOINTPROPMEASURE_F,xfgetxsectpointpropmeasure_f,xfgetxsectpointpropmeasure_f_) xfgetxsectpointpropmeasure_f_
#endif
XMDF_API xid func_name(XFGETXSECTPOINTPROPMEASURE_F,xfgetxsectpointpropmeasure_f,xfgetxsectpointpropmeasure_f_)
              (xid *a_Id, int *a_index, int *a_NumVals, 
                                         double *a_PropId)
{
  int error;

  error = xfGetXSectPointPropMeasure(*a_Id, *a_index, a_NumVals, a_PropId);
  
  return error;

} /* xfgetxsectpointpropmeasure_f_*/


/* --------------------------------------------------------------------------- */
/* FUNCTION  xfsetxsectpointproptype_f_*/
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFSETXSECTPOINTPROPTYPE_F,xfsetxsectpointproptype_f,xfsetxsectpointproptype_f_) XFSETXSECTPOINTPROPTYPE_F
#elif defined FLOWER
#define func_name(XFSETXSECTPOINTPROPTYPE_F,xfsetxsectpointproptype_f,xfsetxsectpointproptype_f_) xfsetxsectpointproptype_f
#else
#define func_name(XFSETXSECTPOINTPROPTYPE_F,xfsetxsectpointproptype_f,xfsetxsectpointproptype_f_) xfsetxsectpointproptype_f_
#endif
XMDF_API xid func_name(XFSETXSECTPOINTPROPTYPE_F,xfsetxsectpointproptype_f,xfsetxsectpointproptype_f_)
              (xid *a_Id, int *a_index, int *a_NumVals, 
                                      int *a_PropId, 
                                      int *a_compression)
{
  int error;

  error = xfSetXSectPointPropType(*a_Id, *a_index, *a_NumVals, a_PropId, *a_compression);
  
  return error;

} /* xfsetxsectpointproptype_f_*/


/* --------------------------------------------------------------------------- */
/* FUNCTION  xfgetxsectpointproptype_f_*/
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFGETXSECTPOINTPROPTYPE_F,xfgetxsectpointproptype_f,xfgetxsectpointproptype_f_) XFGETXSECTPOINTPROPTYPE_F
#elif defined FLOWER
#define func_name(XFGETXSECTPOINTPROPTYPE_F,xfgetxsectpointproptype_f,xfgetxsectpointproptype_f_) xfgetxsectpointproptype_f
#else
#define func_name(XFGETXSECTPOINTPROPTYPE_F,xfgetxsectpointproptype_f,xfgetxsectpointproptype_f_) xfgetxsectpointproptype_f_
#endif
XMDF_API xid func_name(XFGETXSECTPOINTPROPTYPE_F,xfgetxsectpointproptype_f,xfgetxsectpointproptype_f_)
              (xid *a_Id, int *a_index, int *a_NumVals, 
                                      int *a_PropId)
{
  int error;

  error = xfGetXSectPointPropType(*a_Id, *a_index, a_NumVals, a_PropId);
  
  return error;

} /* xfgetxsectpointpropmeasure_f_*/


/* --------------------------------------------------------------------------- */
/* FUNCTION  xfsetxsectpointpropid_f_*/
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFSETXSECTPOINTPROPID_F,xfsetxsectpointpropid_f,xfsetxsectpointpropid_f_) XFSETXSECTPOINTPROPID_F
#elif defined FLOWER
#define func_name(XFSETXSECTPOINTPROPID_F,xfsetxsectpointpropid_f,xfsetxsectpointpropid_f_) xfsetxsectpointpropid_f
#else
#define func_name(XFSETXSECTPOINTPROPID_F,xfsetxsectpointpropid_f,xfsetxsectpointpropid_f_) xfsetxsectpointpropid_f_
#endif
XMDF_API xid func_name(XFSETXSECTPOINTPROPID_F,xfsetxsectpointpropid_f,xfsetxsectpointpropid_f_)
              (xid *a_Id, int *a_NumVals, 
                                          int *a_PropId, int *a_compression)
{
  int error;

  error = xfSetXSectPointPropID(*a_Id, *a_NumVals, a_PropId, *a_compression);
  
  return error;

} /* xfsetxsectpointpropid_f_*/
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfgetxsectpointpropid_f_*/
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFGETXSECTPOINTPROPID_F,xfgetxsectpointpropid_f,xfgetxsectpointpropid_f_) XFGETXSECTPOINTPROPID_F
#elif defined FLOWER
#define func_name(XFGETXSECTPOINTPROPID_F,xfgetxsectpointpropid_f,xfgetxsectpointpropid_f_) xfgetxsectpointpropid_f
#else
#define func_name(XFGETXSECTPOINTPROPID_F,xfgetxsectpointpropid_f,xfgetxsectpointpropid_f_) xfgetxsectpointpropid_f_
#endif
XMDF_API xid func_name(XFGETXSECTPOINTPROPID_F,xfgetxsectpointpropid_f,xfgetxsectpointpropid_f_)
              (xid *a_Id, int *a_NumVals, int *a_PropId)
{
  int error;

  error = xfGetXSectPointPropID(*a_Id, a_NumVals, a_PropId);
  
  return error;

} /* xfgetxsectpointpropid_f_*/
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfsetxsectpointpropname_f_*/
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFSETXSECTPOINTPROPNAME_F,xfsetxsectpointpropname_f,xfsetxsectpointpropname_f_) XFSETXSECTPOINTPROPNAME_F
#elif defined FLOWER
#define func_name(XFSETXSECTPOINTPROPNAME_F,xfsetxsectpointpropname_f,xfsetxsectpointpropname_f_) xfsetxsectpointpropname_f
#else
#define func_name(XFSETXSECTPOINTPROPNAME_F,xfsetxsectpointpropname_f,xfsetxsectpointpropname_f_) xfsetxsectpointpropname_f_
#endif
XMDF_API xid func_name(XFSETXSECTPOINTPROPNAME_F,xfsetxsectpointpropname_f,xfsetxsectpointpropname_f_)
              (xid *a_Id, int *a_NumVals, const char *a_PropId, int *a_StrLen)
{
  int	error;
  char *prop_copy = NULL;

  prop_copy = (char*)malloc((*a_StrLen+1)*sizeof(char));
  strncpy(prop_copy, a_PropId, *a_StrLen);
  prop_copy[*a_StrLen] = '\0';

  error = xfSetXSectPointPropName(*a_Id, *a_NumVals, *a_StrLen, prop_copy);
  
  free(prop_copy);
  return error;

} /* xfsetxsectpointpropname_f_*/
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfgetxsectpointpropname_f_*/
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFGETXSECTPOINTPROPNAME_F,xfgetxsectpointpropname_f,xfgetxsectpointpropname_f_) XFGETXSECTPOINTPROPNAME_F
#elif defined FLOWER
#define func_name(XFGETXSECTPOINTPROPNAME_F,xfgetxsectpointpropname_f,xfgetxsectpointpropname_f_) xfgetxsectpointpropname_f
#else
#define func_name(XFGETXSECTPOINTPROPNAME_F,xfgetxsectpointpropname_f,xfgetxsectpointpropname_f_) xfgetxsectpointpropname_f_
#endif
XMDF_API xid func_name(XFGETXSECTPOINTPROPNAME_F,xfgetxsectpointpropname_f,xfgetxsectpointpropname_f_)
              (xid *a_Id, int *a_NumVals, char *a_PropId, int *a_StrLen)
{
  int	error;
  error = xfGetXSectPointPropName(*a_Id, a_NumVals, a_StrLen, a_PropId);
  return error;

} /* xfgetxsectpointpropname_f_*/
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfgetxsectpointpropnamelen_f_*/
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFGETXSECTPOINTPROPNAMELEN_F,xfgetxsectpointpropnamelen_f,xfgetxsectpointpropnamelen_f_) XFGETXSECTPOINTPROPNAMELEN_F
#elif defined FLOWER
#define func_name(XFGETXSECTPOINTPROPNAMELEN_F,xfgetxsectpointpropnamelen_f,xfgetxsectpointpropnamelen_f_) xfgetxsectpointpropnamelen_f
#else
#define func_name(XFGETXSECTPOINTPROPNAMELEN_F,xfgetxsectpointpropnamelen_f,xfgetxsectpointpropnamelen_f_) xfgetxsectpointpropnamelen_f_
#endif
XMDF_API xid func_name(XFGETXSECTPOINTPROPNAMELEN_F,xfgetxsectpointpropnamelen_f,xfgetxsectpointpropnamelen_f_)
              (xid *a_Id, int *a_NumVals, int *a_StrLen)
{
  int error;

  error = xfGetXSectPointPropNameLen(*a_Id, a_NumVals, a_StrLen);
  
  return error;

} /* xfgetxsectpointpropnamelen_f_*/
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfsetxsectpointpropdesc_f_*/
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFSETXSECTPOINTPROPDESC_F,xfsetxsectpointpropdesc_f,xfsetxsectpointpropdesc_f_) XFSETXSECTPOINTPROPDESC_F
#elif defined FLOWER
#define func_name(XFSETXSECTPOINTPROPDESC_F,xfsetxsectpointpropdesc_f,xfsetxsectpointpropdesc_f_) xfsetxsectpointpropdesc_f
#else
#define func_name(XFSETXSECTPOINTPROPDESC_F,xfsetxsectpointpropdesc_f,xfsetxsectpointpropdesc_f_) xfsetxsectpointpropdesc_f_
#endif
XMDF_API xid func_name(XFSETXSECTPOINTPROPDESC_F,xfsetxsectpointpropdesc_f,xfsetxsectpointpropdesc_f_)
              (xid *a_Id, int *a_NumVals, const char *a_PropId, int *a_StrLen)
{
  int	error;
  char *prop_copy = NULL;

  prop_copy = (char*)malloc((*a_StrLen+1)*sizeof(char));
  strncpy(prop_copy, a_PropId, *a_StrLen);
  prop_copy[*a_StrLen] = '\0';

  error = xfSetXSectPointPropDesc(*a_Id, *a_NumVals, *a_StrLen, prop_copy);
  
  free(prop_copy);
  return error;

} /* xfsetxsectpointpropdesc_f_*/
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfgetxsectpointpropdesc_f_*/
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFGETXSECTPOINTPROPDESC_F,xfgetxsectpointpropdesc_f,xfgetxsectpointpropdesc_f_) XFGETXSECTPOINTPROPDESC_F
#elif defined FLOWER
#define func_name(XFGETXSECTPOINTPROPDESC_F,xfgetxsectpointpropdesc_f,xfgetxsectpointpropdesc_f_) xfgetxsectpointpropdesc_f
#else
#define func_name(XFGETXSECTPOINTPROPDESC_F,xfgetxsectpointpropdesc_f,xfgetxsectpointpropdesc_f_) xfgetxsectpointpropdesc_f_
#endif
XMDF_API xid func_name(XFGETXSECTPOINTPROPDESC_F,xfgetxsectpointpropdesc_f,xfgetxsectpointpropdesc_f_)
              (xid *a_Id, int *a_NumVals, char *a_PropId, int *a_StrLen)
{
  int	error;
  error = xfGetXSectPointPropDesc(*a_Id, a_NumVals, a_StrLen, a_PropId);
  return error;

} /* xfgetxsectpointpropdesc_f_*/
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfgetxsectpointpropdesclen_f_*/
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFGETXSECTPOINTPROPDESCLEN_F,xfgetxsectpointpropdesclen_f,xfgetxsectpointpropdesclen_f_) XFGETXSECTPOINTPROPDESCLEN_F
#elif defined FLOWER
#define func_name(XFGETXSECTPOINTPROPDESCLEN_F,xfgetxsectpointpropdesclen_f,xfgetxsectpointpropdesclen_f_) xfgetxsectpointpropdesclen_f
#else
#define func_name(XFGETXSECTPOINTPROPDESCLEN_F,xfgetxsectpointpropdesclen_f,xfgetxsectpointpropdesclen_f_) xfgetxsectpointpropdesclen_f_
#endif
XMDF_API xid func_name(XFGETXSECTPOINTPROPDESCLEN_F,xfgetxsectpointpropdesclen_f,xfgetxsectpointpropdesclen_f_)
              (xid *a_Id, int *a_NumVals, int *a_StrLen)
{
  int error;

  error = xfGetXSectPointPropDescLen(*a_Id, a_NumVals, a_StrLen);
  
  return error;

} /* xfgetxsectpointpropdesclen_f_*/
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfsetxsectpointpropexclusive_f_*/
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFSETXSECTPOINTPROPEXCLUSIVE_F,xfsetxsectpointpropexclusive_f,xfsetxsectpointpropexclusive_f_) XFSETXSECTPOINTPROPEXCLUSIVE_F
#elif defined FLOWER
#define func_name(XFSETXSECTPOINTPROPEXCLUSIVE_F,xfsetxsectpointpropexclusive_f,xfsetxsectpointpropexclusive_f_) xfsetxsectpointpropexclusive_f
#else
#define func_name(XFSETXSECTPOINTPROPEXCLUSIVE_F,xfsetxsectpointpropexclusive_f,xfsetxsectpointpropexclusive_f_) xfsetxsectpointpropexclusive_f_
#endif
XMDF_API xid func_name(XFSETXSECTPOINTPROPEXCLUSIVE_F,xfsetxsectpointpropexclusive_f,xfsetxsectpointpropexclusive_f_)
              (xid *a_Id, int *a_NumVals, 
                                           int *a_PropId, int *a_compression)
{
  int error;

  error = xfSetXSectPointPropExclusive(*a_Id, *a_NumVals, a_PropId, *a_compression);
  
  return error;

} /* xfsetxsectpointpropexclusive_f_*/
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfgetxsectpointpropexclusive_f_*/
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFGETXSECTPOINTPROPEXCLUSIVE_F,xfgetxsectpointpropexclusive_f,xfgetxsectpointpropexclusive_f_) XFGETXSECTPOINTPROPEXCLUSIVE_F
#elif defined FLOWER
#define func_name(XFGETXSECTPOINTPROPEXCLUSIVE_F,xfgetxsectpointpropexclusive_f,xfgetxsectpointpropexclusive_f_) xfgetxsectpointpropexclusive_f
#else
#define func_name(XFGETXSECTPOINTPROPEXCLUSIVE_F,xfgetxsectpointpropexclusive_f,xfgetxsectpointpropexclusive_f_) xfgetxsectpointpropexclusive_f_
#endif
XMDF_API xid func_name(XFGETXSECTPOINTPROPEXCLUSIVE_F,xfgetxsectpointpropexclusive_f,xfgetxsectpointpropexclusive_f_)
              (xid *a_Id, int *a_NumVals, int *a_PropId)
{
  int error;

  error = xfGetXSectPointPropExclusive(*a_Id, a_NumVals, a_PropId);
  
  return error;

} /* xfgetxsectpointpropexclusive_f_*/
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfsetxsecttoponame_f_*/
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFSETXSECTTOPONAME_F,xfsetxsecttoponame_f,xfsetxsecttoponame_f_) XFSETXSECTTOPONAME_F
#elif defined FLOWER
#define func_name(XFSETXSECTTOPONAME_F,xfsetxsecttoponame_f,xfsetxsecttoponame_f_) xfsetxsecttoponame_f
#else
#define func_name(XFSETXSECTTOPONAME_F,xfsetxsecttoponame_f,xfsetxsecttoponame_f_) xfsetxsecttoponame_f_
#endif
XMDF_API xid func_name(XFSETXSECTTOPONAME_F,xfsetxsecttoponame_f,xfsetxsecttoponame_f_)
              (xid *a_Id, int *a_NumVals, 
	       const char *a_PropId, int *a_StrLen)
{
  int	error;
  char *prop_copy = NULL;

  prop_copy = (char*)malloc((*a_StrLen+1)*sizeof(char));
  strncpy(prop_copy, a_PropId, *a_StrLen);
  prop_copy[*a_StrLen] = '\0';

  error = xfSetXSectTopoName(*a_Id, *a_NumVals, *a_StrLen, prop_copy);
  
  free(prop_copy);
  return error;

} /* xfsetxsecttoponame_f_*/
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfgetxsecttoponame_f_*/
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFGETXSECTTOPONAME_F,xfgetxsecttoponame_f,xfgetxsecttoponame_f_) XFGETXSECTTOPONAME_F
#elif defined FLOWER
#define func_name(XFGETXSECTTOPONAME_F,xfgetxsecttoponame_f,xfgetxsecttoponame_f_) xfgetxsecttoponame_f
#else
#define func_name(XFGETXSECTTOPONAME_F,xfgetxsecttoponame_f,xfgetxsecttoponame_f_) xfgetxsecttoponame_f_
#endif
XMDF_API xid func_name(XFGETXSECTTOPONAME_F,xfgetxsecttoponame_f,xfgetxsecttoponame_f_)
              (xid *a_Id, int *a_NumVals, char *a_PropId, int *a_StrLen)
{
  int	error;
  error = xfGetXSectTopoName(*a_Id, a_NumVals, a_StrLen, a_PropId);
  return error;

} /* xfgetxsecttoponame_f_*/
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfgetxsecttoponamelen_f_*/
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFGETXSECTTOPONAMELEN_F,xfgetxsecttoponamelen_f,xfgetxsecttoponamelen_f_) XFGETXSECTTOPONAMELEN_F
#elif defined FLOWER
#define func_name(XFGETXSECTTOPONAMELEN_F,xfgetxsecttoponamelen_f,xfgetxsecttoponamelen_f_) xfgetxsecttoponamelen_f
#else
#define func_name(XFGETXSECTTOPONAMELEN_F,xfgetxsecttoponamelen_f,xfgetxsecttoponamelen_f_) xfgetxsecttoponamelen_f_
#endif
XMDF_API xid func_name(XFGETXSECTTOPONAMELEN_F,xfgetxsecttoponamelen_f,xfgetxsecttoponamelen_f_)
              (xid *a_Id, int *a_NumVals, int *a_StrLen)
{
  int error;

  error = xfGetXSectTopoNameLen(*a_Id, a_NumVals, a_StrLen);
  
  return error;

} /* xfgetxsecttoponamelen_f_*/
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfsetxsecttopodesc_f_*/
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFSETXSECTTOPODESC_F,xfsetxsecttopodesc_f,xfsetxsecttopodesc_f_) XFSETXSECTTOPODESC_F
#elif defined FLOWER
#define func_name(XFSETXSECTTOPODESC_F,xfsetxsecttopodesc_f,xfsetxsecttopodesc_f_) xfsetxsecttopodesc_f
#else
#define func_name(XFSETXSECTTOPODESC_F,xfsetxsecttopodesc_f,xfsetxsecttopodesc_f_) xfsetxsecttopodesc_f_
#endif
XMDF_API xid func_name(XFSETXSECTTOPODESC_F,xfsetxsecttopodesc_f,xfsetxsecttopodesc_f_)
              (xid *a_Id, int *a_NumVals, const char *a_PropId, int *a_StrLen)
{
  int	error;
  char *prop_copy = NULL;

  prop_copy = (char*)malloc((*a_StrLen+1)*sizeof(char));
  strncpy(prop_copy, a_PropId, *a_StrLen);
  prop_copy[*a_StrLen] = '\0';

  error = xfSetXSectTopoDesc(*a_Id, *a_NumVals, *a_StrLen, prop_copy);
  
  free(prop_copy);
  return error;

} /* xfsetxsecttopodesc_f_*/
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfgetxsecttopodesc_f_*/
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFGETXSECTTOPODESC_F,xfgetxsecttopodesc_f,xfgetxsecttopodesc_f_) XFGETXSECTTOPODESC_F
#elif defined FLOWER
#define func_name(XFGETXSECTTOPODESC_F,xfgetxsecttopodesc_f,xfgetxsecttopodesc_f_) xfgetxsecttopodesc_f
#else
#define func_name(XFGETXSECTTOPODESC_F,xfgetxsecttopodesc_f,xfgetxsecttopodesc_f_) xfgetxsecttopodesc_f_
#endif
XMDF_API xid func_name(XFGETXSECTTOPODESC_F,xfgetxsecttopodesc_f,xfgetxsecttopodesc_f_)
              (xid *a_Id, int *a_NumVals, char *a_PropId, int *a_StrLen)
{
  int	error;
  error = xfGetXSectTopoDesc(*a_Id, a_NumVals, a_StrLen, a_PropId);
  return error;

} /* xfgetxsecttopodesc_f_*/
/* --------------------------------------------------------------------------- */
/* FUNCTION  xfgetxsecttopodesclen_f_*/
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFGETXSECTTOPODESCLEN_F,xfgetxsecttopodesclen_f,xfgetxsecttopodesclen_f_) XFGETXSECTTOPODESCLEN_F
#elif defined FLOWER
#define func_name(XFGETXSECTTOPODESCLEN_F,xfgetxsecttopodesclen_f,xfgetxsecttopodesclen_f_) xfgetxsecttopodesclen_f
#else
#define func_name(XFGETXSECTTOPODESCLEN_F,xfgetxsecttopodesclen_f,xfgetxsecttopodesclen_f_) xfgetxsecttopodesclen_f_
#endif
XMDF_API xid func_name(XFGETXSECTTOPODESCLEN_F,xfgetxsecttopodesclen_f,xfgetxsecttopodesclen_f_)
              (xid *a_Id, int *a_NumVals, int *a_StrLen)
{
  int error;

  error = xfGetXSectTopoDescLen(*a_Id, a_NumVals, a_StrLen);
  
  return error;

} /* xfgetxsecttopodesclen_f_*/

/* --------------------------------------------------------------------------- */
/* FUNCTION  xfsetwkt_f_*/
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFSETWKT_F,xfsetwkt_f,xfsetwkt_f_) XFSETWKT_F 
#elif defined FLOWER
#define func_name(XFSETWKT_F,xfsetwkt_f,xfsetwkt_f_) xfsetwkt_f
#else
#define func_name(XFSETWKT_F,xfsetwkt_f,xfsetwkt_f_) xfsetwkt_f_
#endif
XMDF_API xid func_name(XFSETWKT_F, xfsetwkt_f, xfsetwkt_f_)
              (xid* a_Id, char* a_wkt, int *a_StrLen)
{
  int error;
  char* str_copy = NULL;

  str_copy = (char*)malloc((*a_StrLen + 1) * sizeof(char));
  strncpy(str_copy, a_wkt, *a_StrLen);
  str_copy[*a_StrLen] = '\0';
  error = xfSetWKT(*a_Id, a_wkt);
  free(str_copy);

  return error;

} /* xfsetwkt_f_*/

/* --------------------------------------------------------------------------- */
/* FUNCTION  xfSetAttributeString_f_*/
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFSETATTRIBUTESTRING_F,xfsetattributestring_f,xfsetattributestring_f_) XFSETATTRIBUTESTRING_F 
#elif defined FLOWER
#define func_name(XFSETATTRIBUTESTRING_F,xfsetattributestring_f,xfsetattributestring_f_) xfsetattributestring_f
#else
#define func_name(XFSETATTRIBUTESTRING_F,xfsetattributestring_f,xfsetattributestring_f_) xfsetattributestring_f_
#endif
XMDF_API xid func_name(XFSETATTRIBUTESTRING_F, xfsetattributestring_f, xfsetattributestring_f_)
              (xid* a_ID, char* a_Name, int *name_len, char* a_String, int *string_len)
{

  int error;
  char* name_copy = NULL;
  char* str_copy = NULL;

  name_copy = (char*)malloc((*name_len + 1) * sizeof(char));
  str_copy = (char*)malloc((*string_len + 1) * sizeof(char));
  strncpy(name_copy, a_Name, *name_len);
  strncpy(str_copy, a_String, *string_len);
  name_copy[*name_len] = '\0';
  str_copy[*string_len] = '\0';

  error = xfSetAttributeString(*a_ID, name_copy, str_copy);
  free(name_copy);
  return error;
} /* xfSetAttributeString_f_ */


/* --------------------------------------------------------------------------- */
/* FUNCTION  xfgetwktstringsize_f_*/
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFGETWKTSTRINGSIZE_F,xfgetwktstringsize_f,xfgetwktstringsize_f_) XFGETWKTSTRINGSIZE_F 
#elif defined FLOWER
#define func_name(XFGETWKTSTRINGSIZE_F,xfgetwktstringsize_f,xfgetwktstringsize_f_) xfgetwktstringsize_f
#else
#define func_name(XFGETWKTSTRINGSIZE_F,xfgetwktstringsize_f,xfgetwktstringsize_f_) xfgetwktstringsize_f_
#endif
XMDF_API xid func_name(XFGETWKTSTRINGSIZE_F,xfgetwktstringsize_f,xfgetwktstringsize_f_)
              (xid *a_Id, int *a_StrLen)
{
  int error;

  error = xfGetWKTStringSize(*a_Id, a_StrLen);
  
  return error;

} /* xfgetwktstringsize_f_*/

/* --------------------------------------------------------------------------- */
/* FUNCTION  xfgetwkt_f_*/
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFGETWKT_F,xfgetwkt_f,xfgetwkt_f_) XFGETWKT_F
#elif defined FLOWER
#define func_name(XFGETWKT_F,xfgetwkt_f,xfgetwkt_f_) xfgetwkt_f
#else
#define func_name(XFGETWKT_F,xfgetwkt_f,xfgetwkt_f_) xfgetwkt_f_
#endif
XMDF_API xid func_name(XFGETWKT_F,xfgetwkt_f,xfgetwkt_f_)
              (xid *a_Id, char *a_String)
{
  int error;

  error = xfGetWKT (*a_Id, a_String);
  
  return error;

} /* xfgetwkt_f_*/

/* --------------------------------------------------------------------------- */
/* FUNCTION  xfGetCoordVersion_f_*/
/* PURPOSE    */
/* NOTES      */
/* --------------------------------------------------------------------------- */
#ifdef func_name
#undef func_name
#endif
#if defined FUPPER

#define func_name(XFGETCOORDVERSION_F,xfgetcoordversion_f,xfgetcoordversion_f_) XFGETCOORDVERSION_F
#elif defined FLOWER
#define func_name(XFGETCOORDVERSION_F,xfgetcoordversion_f,xfgetcoordversion_f_) xfgetcoordversion_f
#else
#define func_name(XFGETCOORDVERSION_F,xfgetcoordversion_f,xfgetcoordversion_f_) xfgetcoordversion_f_
#endif
XMDF_API xid func_name(XFGETCOORDVERSION_F,xfgetcoordversion_f,xfgetcoordversion_f_)
              (xid *a_Id, int *a_Version)
{
  int error;

  error = xfGetCoordVersion (*a_Id, a_Version);
  
  return error;

} /* xfGetCoordVersion_f_*/
