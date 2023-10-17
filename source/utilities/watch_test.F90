!==========================================================
module watch_test
!==========================================================    
#include "CMS_cpp.h"
    implicit none
    
contains    

#ifdef UNIT_TEST
!**********************************************************
    subroutine watch_test_all()
! Destroys the watches    
! written by Alex Sanchez, USACE-CHL        
!**********************************************************
    use watch_lib, only: watch_stop, watch_start, watch_default, watch_init, watch_print, watch_output, watch_destroy
    
    implicit none    
    integer, parameter :: n = 100000
    integer :: i,j,iwatch
    real, dimension(n) :: a,b,c
    
    !Initialize
    call watch_default
    call watch_init
    call watch_print
    
    !--- Test 1 --------------------------
    call watch_start('Watch_test 1')    
    !Do some work
    do j=1,100
      a = 1.0 + j
      b = 3.0 - j + a
      do i=1,n
        c(i) = sin(a(i)) + cos(b(i)) + (i*3.0)**3
      enddo
    enddo    
    call watch_stop('Watch_test 1')

    !--- Test 2 --------------------------
    call watch_start('Watch_test 2',iwatch)    
    !Do some work
    do j=1,100
      a = 1.0 + j
      b = 3.0 - j + a
      do i=1,n
        c(i) = sin(a(i)) + cos(b(i)) + (i*3.0)**3
      enddo
    enddo    
    call watch_stop(watch_index=iwatch)
    
    !Finalize
    call watch_output
    call watch_destroy
    
    return
    end subroutine watch_test_all  
#endif

end module watch_test    