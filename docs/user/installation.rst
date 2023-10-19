.. _installation:

Installation and Usage
======================

CMS is a Fortran package that can be compiled on Windows or Linux from source. For the source code, see the `GitHub repository <https://github.com/erdc/cms2d>`_.

Windows
-------
- Microsoft Visual Studio 2019 (Community version is fine) 
- Intel oneAPI Base and HPC kits
- `Get latest CMS 5.3 64-bit executable <https://cirpwiki.info/wiki/CMS_Releases>`_

1) Install MS Visual Studio 2019
	* Make sure **"Desktop development with C++"** is selected at a minimum
2) Install Intel oneAPI Base kit (C++ language only)
3) Install Intel oneAPI HPC kit (adds Fortran)

Note: It has been mentioned there may be an Intel oneAPI Fortran only build which may be possibly used but we have not yet tested it.

Linux
-----
- Gnu compiler (tested with gcc 4.8.5)
- CMake 2.8 or higher
- `Download latest CMS source <https://github.com/erdc/cms2d>`_

1) Clone a copy of the CMS source from the repository.
2) Change directories to the "source" directory
3) Run the command **`cmake  .`**
4) Run the command **`make cms`**

This will leave an executable named 'cms' in the local working directory.

To clean all of the intermediate files and old executable, type 'make clean' from the 
"source" directory.

Note: Linux builds have not been fully tested on HPC platform and some extra configuration may be necessary.

Running CMS
===========
CMS can be run on Windows from within the SMS interface or from the command prompt. On Linux, CMS is run from the command prompt. 

The following options are available from the command prompt:

* Wave Model only:
	> <executable>  <wave>.sim
* Flow Model only: 
	> <executable>  <flow>.cmcards
	
If steering is desired with both Flow and Wave, there are two options:

1) Include necessary Wave steering information within the <flow>.cmcard file to include the cards:
	* CMS-WAVE_SIM_FILE	 <wave>.sim  (include a path if the files aren't in the same folder)
	* STEERING_INTERVAL  <value>     (<value> is in hours)
	* WAVE_WATER_LEVEL   <value>     ('LAST TIME STEP', 'TIDAL', or 'TIDAL_PLUS_VARIATION')
	Then run the Flow model only as indicated above.
2) Specify both the flow and wave paramter files and give additional arguments at the end.
	> <executable>  <flow>.cmcards  <wave>.sim  <interval>  
