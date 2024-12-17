## NOTICE:
CMS Version 5.4.  The CMS code is supplied by the U.S. Army Corps of Engineers, Engineer Research and Development Center (USACE-ERDC). Some portions of this software were written by employees of ERDC-CHL, a laboratory of the U.S. Army Corps of Engineers, and are not under copyright. Other portions are copyright their creators and licensed under licenses listed in the included documentation. USACE-ERDC makes no guarantees about the results, or appropriateness of outputs, obtained from CMS.  

## LIST OF CONDITIONS:
Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
- Redistributions of source code must retain the above notice, this list of conditions, and the following disclaimer.
- Redistributions in binary form must reproduce the above notice, this list of conditions, and the following disclaimer in the documentation and/or other materials provided with the distribution.
- The names of the U.S. Government, the U.S. Army Corps of Engineers, Engineer Research and Development Center, the Coastal and Hydraulics Laboratory may not be used to endorse or promote products derived from this software without specific prior written permission.  Nor may the names of its contributors be used to endorse or promote products derived from this software without specific prior written permission.

## DISCLAIMER:
THIS SOFTWARE IS PROVIDED BY THE U.S. ARMY ENGINEER RESEARCH AND DEVELOPMENT CENTER, COASTAL AND HYDRAULICS LABORATORY (ERDC-CHL) "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL ERDC-CHL BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

## THIRD-PARTY USE NOTICE:
This version of CMS uses the following libraries and source code:

### XMDF Binary Library [Aquaveo]
This library allows for input/output using platform-independent binary files to save disk space and for better communication with the Surface-water Modeling System GUI.

Copyright 2024, Aquaveo, LLC
	
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

### General Cartographic Coordinates Transformation Package (GCTP) [USGS]
This project contains library code needed to map information between horizontal projections.

Documentation - https://www.cmascenter.org/ioapi/documentation/all_versions/html/GCTP.pdf

The General Cartographic Transformation Package (GCTP) is a system of computer subroutines, written in FORTRAN, designed to permit the transformation of coordinate pairs from one map projection to another. GCTP should not be used to transform coordinates between spheroids, because a datum shift should be applied to the geographic coordinates in most cases. It is a subroutine package that must be linked to and called by other FORTRAN programs. The GCTP is the standard computer software used by the National Mapping Division for map projection computations. The mathematical algorithms used in GCTP meet the accuracy specifications of USGS Circular 878-B, "Representation of Geographic Point Locations for Information Interchange," which has been proposed as Federal Information Processing Standards Publication 70-1. This software is approved for use with all products of the National Mapping Program.

Created by USGS. No license or copyright.

### State Plane Coordinate System (SPCS) [NOAA NGS]
This module contains all the subroutines needed to convert NAD 83 state plane coordinates to NAD 83 geographic positions and conversely. Includes defining constants for NAD 83 coordinate zones. State plane coordinates are entered or computed to 1 mm accuracy, while the latitudes and longitudes entered or computed correspond to approximately 0.3 mm accuracy.

Original Source: https://www.ngs.noaa.gov/PC_PROD/SPCS83/

This program and supporting information is furnished by the Government of the United States of America, and is accepted and used by the recipient with the understanding that the United States Government makes no warranties, express or implied, concerning the accuracy, completeness, reliability, or suitability of this program, of its constituent parts, or of any supporting data. The Government of the United States of America shall be under no liability whatsoever resulting from any use of this program. This program should not be relied upon as the sole basis for solving a problem whose incorrect solution could result in injury to person or property. This program is property of the Government of the United States of America. Therefore, the recipient further agrees not to assert proprietary rights therein and not to represent this program to anyone as being other than a Government program.

Created by NOAA. No license or copyright.
   
### UTM to GEO [MIT]
These routines convert UTM to Lat/Longitude and vice-versa, using the WGS-84 (GPS standard) or Clarke 1866 Datums.

Copyright (C) 1998, Massachusetts Institute of Technology, All Rights Reserved.

### SZIP
Precompiled binary library included in the XMDF library mentioned above.

Copyright (C) 2001 Science & Technology Corporation @ UNM. All rights released.  
Copyright (C) 2003 Lowell H. Miles and Jack A. Venbrux.  Licensed to ICs, LLC, for distribution by the University of Illinois' National Center for Supercomputing Applications as a part of the HDF data storage and retrieval file format and software library products package.  All rights reserved.  Do not modify or use for other purposes.

Revocable, royalty-free, nonexclusive sublicense to use SZIP compression software routines and underlying patents for non-commercial, scientific use only is hereby granted by ICs, LLC, to users of and in conjunction with HDF data storage and retrieval file format and software library products.