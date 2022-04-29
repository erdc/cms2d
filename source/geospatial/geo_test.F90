!==========================================================
module geo_test
! Geospatial library unit tests
! written by Alex Sanchez, USACE-CHL
!==========================================================    
#include "CMS_cpp.h"    
    implicit none
    
contains
#ifdef UNIT_TEST
!**********************************************************
    subroutine geo_test_all()  
! written by Alex Sanchez, USACE-CHL        
!**********************************************************
    use geo_lib
    use diag_lib
    use prec_def
    implicit none
    integer :: ierr
    integer, parameter :: np = 4
    real(ikind) :: xp(np),yp(np)
    real(ikind) :: x1,x2,x3,x4,y1,y2,y3,y4,px,py
    character(len=200) :: msg,msg2
    
    !--- Polynomial vertex sorting ------
    xp = (/319810.4, 320738.6, 321213.7, 321965.6/)
    yp = (/286610.4, 284836.1, 287794.6, 286171.8/)
    call poly_sort(np,np,xp,yp)
    
    !--- Intercept between two lines -----------------
    x1 = 0.0; y1 = 2.0; x2 = 2.0; y2 = 0.0;
    x3 = 0.0; y3 = 0.0; x4 = 2.0; y4 = 2.0;    
    call line_line_intercept(x1,y1,x2,y2,x3,y3,x4,y4,px,py,ierr)
    write(msg,*) 'Calculated: px = ',px,', py = ',py,' ierr = ',ierr
    write(msg2,*) 'Correct:    px = ',1.0,', py = ',1.0,' ierr = ',0
    
    x1 = 1.5; y1 = 0.6; x2 = 2.0; y2 = 0.0;
    x3 = 0.0; y3 = 0.0; x4 = 2.0; y4 = 2.0;
    call line_line_intercept(x1,y1,x2,y2,x3,y3,x4,y4,px,py,ierr)
    write(msg,*) 'Calculated: px = ',px,', py = ',py,' ierr = ',ierr
    write(msg2,*) 'Correct:    px = ',1.0909,', py = ',1.0909,' ierr = ',1
    
    return
    end subroutine geo_test_all
#endif

end module geo_test
