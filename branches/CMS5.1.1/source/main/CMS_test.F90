!=====================================================
module CMS_test
!=====================================================
#include "CMS_cpp.h"
    implicit none
    
contains   
#ifdef UNIT_TEST
!**************************************************    
    subroutine CMS_test_run
! Performs CMS tests of individual routines and
! outputs the results to the diagnostic file
! written by Alex Sanchez, USACE-CHL
!**************************************************
    use comvarbl
    use watch_test
    use geo_test
    use math_test
    use time_test
    !use wave_test
    implicit none
    
    dgfile = 'CMS_DIAG.txt' !Diagnostic file is always in flow path
    dgunit = 9
    
    call watch_test_all
    call math_test_all
    call geo_test_all
    call time_test_all
    !call wave_test_all

    write(*,*) ' Press any key to continue'
    read(*,*)
    
    return
    endsubroutine CMS_test_run
#endif

endmodule CMS_test    
    
    
    