## Coastal Modeling System (CMS)

CMS is a coastal modeling system that couples a wave, circulation, and morphology model together to 
get better predictions in the near-shore

CMS is developed by the USACE - Coastal and Hydraulics Laboratory


To build CMS On linux you must have CMAKE softward above 3.4
1) In file `"source/main/CMS_cpp.h"`, change `#def WIN_OS` to `#undef WIN_OS`

2) `cd trunk/source`

3) `cmake .`

4) `make`
This will leave you with a cms executable in your local directory 

To clean all of the intermediate files and old executable, type 'make clean' from the 
"trunk/source" directory.
