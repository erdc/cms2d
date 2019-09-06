              U.S. Army Corps of Engineers
            Coastal Inlets Research Program
                Coastal Modeling System
 Coupled Hydrodynamic, Wave, and Sediment Transport Model
       For the latest version of CMS please visit
          http://cirp.usace.army.mil/products/

    By using this software the user has agreed to the
    terms and conditions of CMS license agreement.
    A copy of the license can be obtained from the
    website shown above.
     
    This software is distributed on an "AS IS" basis
    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND,
    either express or implied.


 The command line accepts two arguments
 CMS-Wave requires that its input file (*.sim)
 be the first file.
 Therefore, the following formats are valid
   1. Wave model only: 
        > cms2d_v4p1r35-x32p.exe Wave.sim
   2. Flow model only:
        > cms2d_v4p1r35-x32p.exe Flow.cmcards
   3. Flow and Wave models in steering
        > cms2d_v4p1r35-x32p.exe Wave.sim Flow.cmcards
   4. Flow and Wave models in steering
        > cms2d_v4p1r35-x32p.exe Wave.sim Flow.cmcards 3.0 1


 Compiling in Microsoft Visual Studio (2010)
    Release Configuration Properties: 
       Project Settings:
          Linker | Input - Additional Dependencies => "..\external_libraries\XMDF\xmdf2.0.lib"  (for 32 bit) 
                                                 => "..\external_libraries\XMDF\xmdf2.0x64.lib" (for 64 bit)
          Linker | System - Stack Reserve Size => 46214400
          Fortran | Optimization - Optimization => Maximum Speed plus Higher Level Optimizations
          Fortran | Optimization - Global Optimizations => Yes (/Og) 
          Fortran | Language - Process OpenMP Directives => Generate Parallel Code (/Qopenmp)
          Fortran | Floating Point | Floating Point Exception Handling =>
                     Underflow gives 0.0; Abort on other IEEE exceptions (/fpe:0) (giving problems with ICCG solver)
                  or Produce NaN, signed infinities, and denormal results (more stable for some reason)
          Fortran | Run-time - Generate Traceback Information  => Yes
          Fortran | Run-time - Check Array and Sting Bounds => Yes
          Fortran | Run-time - Check Uninitialized Variables => Yes
          Fortran | Preprocessor - Preprocess source file => Yes
          Fortran | Preprocessor - Additional include directories => ..\source\main\

    Debug Configuration Properties:
       Project Settings:
          Linker | Input - Additional Dependencies => "..\libraries\XMDF\xmdf2.0.lib"  (for 32 bit) 
                                                   => "..\libraries\XMDF\xmdf2.0x64.lib" (for 64 bit)
          Linker | System - Stack Reserve Size => 46214400
          Linker | Debugging - Generate Debug Infor =>Yes (/DEBUG)
          Linker | Input - Ignore specific library => LIBCMT.lib (bug in Visual Studio 2010)
          Fortran | Debugging - Debug Information Format => Full(/debug:full)
          Fortran | Debugging - Information for PARAMETER Constants => All(/debug-parameters:all)
          Fortran | Diagnostics - Compile time diagnostics => Show All (/warn:all)
          Fortran | Preprocessor - Preprocess source file => Yes
          Fortran | Preprocessor - Additional include directories => ..\source\main\
       Specific file Properties:
          CMS-Wave_v3-2W_10Oct2011.f90
             Fortran | Diagnostics - Compile time diagnostics => Disable All (/warn:none)
          gctp.f
             Fortran | Diagnostics - Compile time diagnostics => Disable All (/warn:none)
          spcs83_combined.F
             Fortran | Diagnostics - Compile time diagnostics => Disable All (/warn:none) 
          utm2geo.F
             Fortran | Language - Fixed Form Line Length => 132 Columns (/extend_source:132)
     
     
 CMake Compiling (mainly for Linux)
   To build CMS On linux you must have CMAKE software above 3.4
   
   Instructions:
     1. Go to the directory containing the "source" subdirectory (generally "source"). 
     2. Type cmake . at the prompt
     3. Type make clean
     4. Type make
     5. If no errors, an exectuable named 'cms' will be in the working directory. Rename the executable as needed.
