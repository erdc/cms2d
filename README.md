# Coastal Modeling System (CMS)
CMS is a coastal modeling system that couples a wave, circulation, and morphology model together to get better predictions in the near-shore.

CMS is developed by the USACE - Coastal and Hydraulics Laboratory

Documentation for the CMS is available from the [CIRP Wiki](https://cirpwiki.info/wiki/CMS).

# Building/Editing
## Obtain the source code using the OpenSource GitHub repository.
- Go the location you wish to place the source.
- If using a command prompt with Git, run `git clone https://github.com/erdc/cms2d.git`
- If using a browser, retrieve a zip file with the source the following URL and unzip.
  - https://github.com/erdc/cms2d/archive/refs/heads/master.zip

## Windows platform
Building CMS on Windows requires use of the Microsoft Visual Studio GUI. This repository was created using the Intel Fortran compiler. The versions of the GUI and compiler primarily being used are as follows:
- MS Visual Studio 2019 or 2022
- Intel oneAPI 2021.4 (vs2019) or 2024.1 (vs2022)

Once the sofware has been installed, the user should load the appropriate **solution file** for the UI.

## Linux platform
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
    - There may be a few warnings during XMDF compilation, but it should create a working library.
	
- If no errors, an executable named 'cms' will be in the 'source' directory. Move and rename the executable as needed.
- This will leave you with an executable named 'cms' in your the "source" directory. 
Note: Linux builds have not been fully tested on HPC platform and some extra configuration may be necessary.

# File/Folder structure
Below is a list of folders and descriptions of the contents.
- **docs** - this folder contains the documentation files needed for the [CMS Documentation page](https://cms2d.readthedocs.io/) that is in progress.
- **external_libraries** - subfolders containing third-party source and binary libraries to link with.
- **Intel_vs2019** - this folder contains the solution file needed to load into Visual Studio 2019.
- **Intel_vs2022** - this folder contains the solution file needed to load into Visual Studio 2022.
- **source** - contains subfolders for all the main CMS code functionality separated into process type.
- **testing** - contains CI/CD tests

Additional files in the main folder and their descriptions:
- [README.md](README.md) - this file.
- [UNLICENSE.md](UNLICENSE.md) - contains information regarding the use of this repository and the third-party licenses.
- [CMS Terms and Conditions.txt](<CMS Terms and Conditions.txt>) - additional terms and conditions.
- [.gitlab-ci.yml](.gitlab-ci.yml) - contains needed continuous integration build information.
- [pyproject.toml](pyproject.toml) - contains needed information for the ReadTheDocs documentation pages.

Additional files in the "source" folder and their descriptions:
- [CMakeLists.txt](source/CMakeLists.txt) - contains information for compiling on Linux using CMake.
- [logsheet.txt](source/logsheet.txt) - contains information about changes between versions of CMS.
- [readme.txt](source/readme.txt) - contains information about configuration settings for Visual Studio on Windows and proper usages of the CMS executable from the command line.
