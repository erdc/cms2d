!==========================================================
module time_test
! Time library unit tests
! written by Alex Sanchez, USACE-CHL
!==========================================================    
#include "CMS_cpp.h"    
    implicit none
    
contains
#ifdef UNIT_TEST
!**********************************************************
    subroutine time_test_all
! written by Alex Sanchez, USACE-CHL        
!**********************************************************
    use time_lib
    use diag_lib
    use prec_def
    implicit none
    integer :: iyr,imo,iday,ihr,imin,isec,imilsec
    real(ikind) :: tjul,tjul2
    real(8) :: timesec
    character(len=200) :: msg
    
645 format(A,1x,I4,1x,5(I2,1x),I3)    
    
    call diag_print_message(' ','--- Timing Unit Tests ---')
    tjul=367.5
    write(msg,*) 'tjul = ',tjul
    call diag_print_message(msg)
    
    call julian2calendar(tjul,iyr,imo,iday,ihr,imin,isec) ! The origin of the Julian date is such that tjul=367.0 on 1901-01-01 00:00:00
    write(msg,645) 'iyr,imo,iday,ihr,imin,isec = ',iyr,imo,iday,ihr,imin,isec
    call diag_print_message(' ','call julian2calendar(tjul,iyr,imo,iday,ihr,imin,isec)','msg = ',msg)
    
    call calendar2julian(iyr,imo,iday,ihr,imin,isec,tjul2)
    write(msg,*) 'tjul2 = ',tjul2
    call diag_print_message(' ','call calendar2julian(iyr,imo,iday,ihr,imin,isec,tjul2)','msg = ',msg)
    
    iyr = 2014; imo = 6; iday = 24; ihr = 22; imin = 30; isec = 30; imilsec = 212
    write(msg,645) 'iyr,imo,iday,ihr,imin,isec,imilsec = ',iyr,imo,iday,ihr,imin,isec,imilsec
    call diag_print_message(' ',msg)
        
    call time_cal2str(msg,iyr,imo,iday)
    call diag_print_message(' ','call time_cal2str(msg,iyr,imo,iday)','msg = ',msg)
    
    call time_cal2str(msg,iyr,imo,iday,ihr,imin,isec)
    call diag_print_message(' ','call time_cal2str(msg,iyr,imo,iday,ihr,imin,isec)','msg = ',msg)
    
    call time_cal2str(msg,iyr,imo,iday,ihr,imin,isec,imilsec)
    call diag_print_message(' ','call time_cal2str(msg,iyr,imo,iday,ihr,imin,isec,imilsec)','msg = ',msg)
    
    !           days            hrs         min      sec        milsec
    timesec = 4*86400.0d0 + 2*3600.0d0 + 1*60.0d0 + 3*1.0d0 + 355*0.001d0
    msg = 'timesec = 4*86400.0d0 + 2*3600.0d0 + 1*60.0d0 + 3*1.0d0 + 355*0.001d0'
    call diag_print_message('',msg)
    call time_sec2str(timesec,msg)
    call diag_print_message('','call time_sec2str(timesec,msg) ',msg)
    
    return
    end subroutine time_test_all
#endif

end module time_test

!call waveorbrms_jonswap(3.0,10.2480,10.0,ubr)