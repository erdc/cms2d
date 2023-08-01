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
    use diag_def,   only: dgfile, dgunit
    use geo_test,   only: geo_test_all
    use math_test,  only: math_test_all 
    use time_test,  only: time_test_all
    use watch_test, only: watch_test_all

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
    end subroutine CMS_test_run
#endif

end module CMS_test    
    
    
    