!==========================================================
module math_test
! Math library unit tests
! written by Alex Sanchez, USACE-CHL
!==========================================================    
#include "CMS_cpp.h"    
    implicit none
    
contains
#ifdef UNIT_TEST
!**********************************************************
    subroutine math_test_all()  
! written by Alex Sanchez, USACE-CHL        
!**********************************************************
    use math_lib, only: lognpdf, erfinv, logninv, sortup
    use diag_lib, only: diag_print_message
    use prec_def, only: ikind
    
    implicit none
    integer, parameter :: nsort = 4
    integer :: indsort(nsort)
    real(ikind) :: xsort(nsort),xsort2(nsort)
    real(ikind) :: val
    character(len=200) :: msg
    
    call diag_print_message('')
    call diag_print_message('--- Math Unit Tests ----')
    
    !--- Sorting -------------------------
    xsort = (/3.0, 2.0, 5.0, 4.0/)
    xsort2 = xsort
    call sortup(nsort,xsort,indsort)
    xsort2 = xsort2(indsort)
    write(msg,*) 'xsort = ',xsort
    call diag_print_message(msg)
    write(msg,*) 'xsort2 = ',xsort2
    call diag_print_message(msg)

    !--- Functions ----------------------------
    val = lognpdf(1.0,0.5,1.5)
    write(msg,'(A,F6.4,A)') 'lognpdf(1.0) :: Calculated = ',val,', Matlab = 0.2516'
    call diag_print_message(msg)
    
    val = erfinv(0.5)
    write(msg,'(A,F6.4,A)') 'erfinv(0.5) :: Calculated = ',val,', Matlab = 0.4769'
    call diag_print_message(msg)
    
    val = logninv(0.2,0.5,1.0)
    write(msg,'(A,F6.4,A)') 'logninv(0.5) :: Calculated = ',val,', Matlab = 0.7106'
    call diag_print_message(msg)
    
    val = erf(1.0)
    write(msg,'(A,F6.4,A)') 'erf(1.0) :: Calculated = ',val,', Matlab = 0.8427'
    call diag_print_message(msg)
    
    return
    end subroutine math_test_all
#endif

end module math_test