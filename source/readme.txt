              U.S. Army Corps of Engineers
            Coastal Inlets Research Program
                Coastal Modeling System
  Coupled Hydrodynamic, Wave, and Sediment Transport Model
       For the latest executable for CMS please visit
        https://cirp.usace.army.mil/products/cms.php

    By using this software the user has agreed to the
    terms and conditions of CMS license agreement.
    A copy of the license can be obtained from
    https://cirpwiki.info/wiki/CMS_License.
     
    This software is distributed on an "AS IS" basis
    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND,
    either express or implied.


 The command line accepts multiple arguments, the following formats are valid:
    Wave model only: 
        > CMS2D_V5.4.exe Wave.sim
    Flow model only:
        > CMS2D_V5.4.exe Flow.cmcards
    Flow and Wave models in steering (w/default time interval of 3 hours)
        > CMS2D_V5.4.exe Wave.sim Flow.cmcards
   If there are both a Wave and Flow simulation specified, user may specify time interval 
   as 4th argument.

    Flow and Wave models in steering (with specified time interval of 1 hour)
        > CMS2D_V5.4.exe Wave.sim Flow.cmcards 1.0

 There are some tools you can use as well.  To access the tools, use the following format:
   > CMS2D_V5.4.exe  tools           !Case doesn't matter


 Compiling in Microsoft Visual Studio (2019 or greater)
    Release Configuration Properties: 
       Project Settings:
          Linker  | Input - Additional Dependencies 
             => "..\external_libraries\XMDF\xmdf2.2x64.lib"
          Linker  | System - Stack Reserve Size => 92428800
          Fortran | Optimization - Optimization => Maximum Speed plus Higher Level Optimizations
          Fortran | Optimization - Global Optimizations => Yes (/Og) 
          Fortran | Language - Process OpenMP Directives => Generate Parallel Code (/Qopenmp)
          Fortran | Floating Point | Floating Point Exception Handling =>
                     Underflow gives 0.0; Abort on other IEEE exceptions (/fpe:0) (giving problems with 
                     ICCG solver)
                  or Produce NaN, signed infinities, and denormal results (more stable for some reason)
          Fortran | Run-time - Generate Traceback Information  => Yes
          Fortran | Run-time - Check Array and Sting Bounds => Yes
          Fortran | Run-time - Check Uninitialized Variables => Yes
          Fortran | Preprocessor - Preprocess source file => Yes
          Fortran | Preprocessor - Additional include directories => ..\source\main\

    Debug Configuration Properties:
       Project Settings:
          Linker  | Input - Additional Dependencies 
             => "..\external_libraries\XMDF\xmdf2.2dx64.lib"
          Linker  | System - Stack Reserve Size => 92428800
          Linker  | Debugging - Generate Debug Infor =>Yes (/DEBUG)
          Linker  | Input - Ignore specific library => LIBCMT.lib (bug in Visual Studio 2010+)
          Fortran | Debugging - Debug Information Format => Full(/debug:full)
          Fortran | Debugging - Information for PARAMETER Constants => All(/debug-parameters:all)
          Fortran | Diagnostics - Compile time diagnostics => Show All (/warn:all)
          Fortran | Preprocessor - Preprocess source file => Yes
          Fortran | Preprocessor - Additional include directories => ..\source\main\
       Specific file Properties:
	     for the following files:
           - CMS-Wave_v3.4W_Jan2025.f90
           Fortran | Diagnostics - Compile time diagnostics => Disable All (/warn:none) 
     
 CMake Compiling (mainly for Linux)
   To build CMS On linux:
     - You must have CMAKE software 3.5 or above
     - You must have the HDF5 (and Zlib) libraries installed on your system.
        
   Instructions (updated 06/18/25):
     1. Clone repository or unzip CMS Source into a folder.
     2. Go to the root directory of the repo. 
     3. Type `cmake .` at the prompt    !(do not enter the ` characters)
     4. Type `make clean`
     5. Type `make`
        - This will make the XMDF library and the CMS executable using the XMDF library.
        - There may be a few warnings during XMDF compilation, but it should still create a library.
     6. If no errors, an exectuable named 'cms' will be in the 'source' directory. Move and rename the executable as needed.
