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
 * file       Xmdf_f.h.cpp
 *
 * purpose    Used as a wrapper to call C functions from the FORTRAN code in XMDF
 *
 * notes  
 *
 * coded by
 *
 ******************************************************************************/






xid xfCreateFile_f(const char * a_Filename, int* namelen, xid *a_Id, short *a_Overwrite);
xid xfOpenFile_f(const char * a_Filename, int* namelen, xid *a_Id, short *a_ReadOnly);
