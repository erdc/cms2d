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


 Compiling in Microsoft Visual Studio (2010 or greater)
    Release Configuration Properties: 
       Project Settings:
          Linker | Input - Additional Dependencies (replace <###> with working version)
		     => "..\external_libraries\XMDF\xmdf<###>.lib"    (for 32 bit) 
             => "..\external_libraries\XMDF\xmdf<###>x64.lib" (for 64 bit)
          Linker | System - Stack Reserve Size => 46214400
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
          Linker | Input - Additional Dependencies (replace <###> with working version)
		     => "..\libraries\XMDF\xmdf<###>.lib"  (for 32 bit) 
             => "..\libraries\XMDF\xmdf<###>x64.lib" (for 64 bit)
          Linker | System - Stack Reserve Size => 46214400
          Linker | Debugging - Generate Debug Infor =>Yes (/DEBUG)
          Linker | Input - Ignore specific library => LIBCMT.lib (bug in Visual Studio 2010+)
          Fortran | Debugging - Debug Information Format => Full(/debug:full)
          Fortran | Debugging - Information for PARAMETER Constants => All(/debug-parameters:all)
          Fortran | Diagnostics - Compile time diagnostics => Show All (/warn:all)
          Fortran | Preprocessor - Preprocess source file => Yes
          Fortran | Preprocessor - Additional include directories => ..\source\main\
       Specific file Properties:
	     for the following files:
           - CMS-Wave_v3.3W_30Apr2020.f90
           Fortran | Diagnostics - Compile time diagnostics => Disable All (/warn:none) 
     
 CMake Compiling (mainly for Linux)
   To build CMS On linux you must have CMAKE software above 3.4
   
   Instructions:
     1. Go to the directory containing the "source" subdirectory (generally "source"). 
     2. Type `cmake .` at the prompt    !(do not enter the ` character)
     3. Type `make clean`
     4. Type `make`
     5. If no errors, an exectuable named 'cms' will be in the working directory. Rename the executable as needed.
