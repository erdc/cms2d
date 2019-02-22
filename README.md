## Coastal Modeling System (CMS)

CMS is a coastal modeling system that couples a wave, circulation, and morphology model together to get better predictions in the nea-shore

CMS is developed by the USACE - Coastal and Hydraulics Laboratory


to build CMS On linux you must have CMAKE softward above 3.4
in file `"source/main/CMS_cpp.h"`
change `#def WIN_OS` to `#undef WIN_OS`

`cd trunk/source`

`cmake`

`make .`
