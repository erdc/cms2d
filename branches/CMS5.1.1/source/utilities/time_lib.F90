!=================================================================================
module time_lib
! Timing library
!
! Contains:
!   time_cpu - Returns the the CPU time in seconds 
!   time_clock - Returns the the wall clock time in seconds
!   time_jul - Returns the Julian time in seconds
!   julday - Converts a Gregorian calendar date to a Julian date
!   gregdate - Converts a Julian date to a Gregorian calendar date
!   julian2calendar - Converts a Julian time to a Gregorian calendar time
!   calendar2julian - Converts a Gregorian calendar time to a Julian time
!   julianday2calendarmonthday - Converts a julian date to a calendar month and day
!   index2calendar - Converts a time index to a calendar time
!   ramp_func - Computes the ramp function
!   time_cal2str - Converts a calendar date and time to a nicely formatted string
!   time_sec2str - Converts an elapsed time in seconds to nicely formated string using
!                  days, hours, minutes, and seconds
!   time_s2dhs - Converts an elapsed time in seconds to days, hours, minutes, 
!                seconds, and miliseconds
!
! Authors: Alex Sanchez, USACE-CHL
!          Mitch Brown, USACE-CHL
!=================================================================================
    implicit none
    
contains

!**********************************************************
    function time_cpu() result(t)
! Returns the the CPU time in seconds    
! written by Alex Sanchez, USACE-CHL    
!**********************************************************
    implicit none
    real(8) :: t
    
    call cpu_time(t)
   
    return
    endfunction time_cpu
    
!*******************************************************************
    function time_clock() result(t)
! Returns the the wall clock time in seconds
! Note: The wall clock is a 24-hr clock
! written by Alex Sanchez, USACE-CHL    
!*******************************************************************
    implicit none
    real(8) :: t
    integer :: clock_count,clock_rate,clock_max
    
    call system_clock(clock_count,clock_rate,clock_max)    
    t = dble(clock_count)/dble(clock_rate)
   
    return
    endfunction time_clock   
    
!************************************************************************    
    function time_jul(itain) result(t)
! Returns the Julian time in seconds
! The origin of the Julian date is such that 
! time_jul=367.0d0*86400.0d0 on 1901-01-01 00:00:00  
! written by Alex Sanchez, USACE-CHL
!************************************************************************    
    implicit none
    integer,intent(in),optional :: itain(8)
    integer :: ita(8)
    real(8) :: t
    
    if(present(itain))then
      ita = itain 
    else  
      call date_and_time(values=ita)    
    endif  
    t = dble(julday(ita(1),ita(2),ita(3)))*86400.0d0 + dble(ita(5))*3600.0d0 &
           + dble(ita(6))*60.0d0 + dble(ita(7)) + 0.001d0*dble(ita(8))  !julian seconds  
    
    return
    endfunction time_jul

!************************************************************************      
    function julday(iyr,imo,iday)
! DATE ROUTINE julday(iyr,imo,iday) CONVERTS CALENDAR DATE TO JULIAN DATE. 
! SEE CACM 1968 11(10):657, LETTER TO THE
! EDITOR BY HENRY F. FLIEGEL AND THOMAS C. VAN FLANDERN.
! THE ORIGIN IS SUCH THAT julday=367 ON 1/1/1901
!************************************************************************      
    implicit none
    integer,intent(in) :: iyr,imo,iday
    integer :: julday
    
    julday=iday-32075+1461*(iyr+4800+(imo-14)/12)/4 &
           +367*(imo-2-(imo-14)/12*12)/12 &
           -3*((iyr+4900+(imo-14)/12)/100)/4         
	julday = julday - 2415019  !ADDED TO MAKE CONSISTENT WITH EXCELL JULIAN DATE
    
    return
    endfunction julday
      
!************************************************************************       
    subroutine gregdate(jday,iyr,imo,iday)
! COMPUTES THE GREGORIAN CALENDAR DATE (iyr,imo,iday)
! GIVEN THE JULIAN DATE (julday)
! THE ORIGIN OF IS SUCH THAT julday=367 ON 1/1/1901
!************************************************************************ 
    implicit none
    integer,intent(in) :: jday
    integer,intent(out) :: iyr,imo,iday
    integer :: n,h
    
    h=jday+68569+2415019 !ADDED TO MAKE CONSISTENT WITH EXCELL JULIAN DATE      
    n=4*h/146097
    h=h-(146097*n+3)/4
    iyr=4000*(h+1)/1461001
    h=h-1461*iyr/4+31
    imo=80*h/2447
    iday=h-2447*imo/80
    h=imo/11
    imo=imo+2-12*h
    iyr=100*(n-49)+iyr+h      
    
    return
    endsubroutine gregdate

!************************************************************************       
    subroutine julian2calendar(tjul,iyr,imo,iday,ihr,imin,isec)
! The origin of the Julian date is such that 
! tjul=367.0 on 1901-01-01 00:00:00  
!************************************************************************
    use prec_def
    implicit none
    !Input/Output
    real(ikind),intent(in) :: tjul
    integer,intent(out) :: iyr,imo,iday,ihr,imin,isec
    !Internal
    integer :: jday
    real(ikind) :: hr,min,sec
    
    jday = floor(tjul) !Julian days
    call gregdate(jday,iyr,imo,iday)    
    hr = (tjul-floor(tjul))*24.0
    min = (hr-floor(hr))*60.0
    sec = (min-floor(min))*60.0
    ihr = int(hr)
    imin = int(min)
    isec = int(sec)
    
    return
    endsubroutine julian2calendar
    
!************************************************************************       
    subroutine calendar2julian(iyr,imo,iday,ihr,imin,isec,tjul)
! The origin of the Julian date is such that 
! tjul=367.0 on 1901-01-01 00:00:00  
!************************************************************************
    use prec_def
    implicit none
    !Input/Output
    integer,intent(in) :: iyr,imo,iday,ihr,imin,isec
    real(ikind),intent(out) :: tjul
    
    tjul = real(julday(iyr,imo,iday),kind=ikind) &
         + real(ihr,kind=ikind)/24.0 &
         + real(imin,kind=ikind)/1440.0 & 
         + real(isec,kind=ikind)/86400.0 !julian days
    
    return
    endsubroutine calendar2julian

!*****************************************************************************
    subroutine julianday2calendarmonthday(iyr,jday,imo,iday)
! Converts a julian date to a calendar month and day
!
!Input
! iyr - calendar year
! jday - julian day (1-366)
!
!Output
! imo - calendar month
! iday - calendar day
!
! written by Mitch Brown, USACE-CHL
! modified by Alex Sanchez, USACE-CHL
!*****************************************************************************
    use diag_lib
    implicit none      
    integer,intent(in)  :: jday,iyr
    integer,intent(out) :: imo,iday
    integer :: jdays(12)
    character(len=100) :: msg
    
!    real(Double) :: reftime0     !(hli, 02/04/10), Alex: Not used 
      
    ! CALCULATE REFERENCE TIME FOR BINARY DATASETS - 04/25/07 MEB
    if(mod(iyr,4)==0)then
      jdays=(/31,60,91,121,152,182,213,244,274,305,335,366/)
    else
      jdays=(/31,59,90,120,151,181,212,243,273,304,334,365/)
    endif
    if(jday==0)then
      write(msg,*) 'Year: ',iyr,', Julein day: ',jday
      call diag_print_error('Invalid Julian Day: ',msg)
    endif    
    if(jday<=jdays(1))then 
      imo = 1;  iday = jday
    elseif(jday<=jdays(2))then
      imo = 2 ; iday = jday-jdays(imo-1)
    elseif(jday<=jdays(3))then
      imo = 3 ; iday = jday-jdays(imo-1)
    elseif(jday<=jdays(4))then
      imo = 4 ; iday = jday-jdays(imo-1)
    elseif(jday<=jdays(5))then
      imo = 5 ; iday = jday-jdays(imo-1)
    elseif(jday<=jdays(6))then
      imo = 6 ; iday = jday-jdays(imo-1)
    elseif(jday<=jdays(7))then
      imo = 7 ; iday = jday-jdays(imo-1)
    elseif(jday<=jdays(8))then
      imo = 8 ; iday = jday-jdays(imo-1)
    elseif(jday<=jdays(9))then
      imo = 9 ; iday = jday-jdays(imo-1)
    elseif(jday<=jdays(10))then
      imo = 10 ; iday = jday-jdays(imo-1)
    elseif(jday<=jdays(11))then
      imo = 11 ; iday = jday-jdays(imo-1)
    else
      imo = 12 ; iday = jday-jdays(imo-1)
    endif    
    
    return
    endsubroutine julianday2calendarmonthday
    
!**************************************************************************
    function ramp_func(timehrs,rampdur) result(ramp)
! Computes the ramp function used to slowly increase the model forcing
! from zero.
!
! Input:
!   timehrs - Time in hours from start of simulation
!   rampdur - Duration of ramp period in hours
!
! Output:
!  ramp_func - ramp scalar  
!
! written by Alex Sanchez, USACE-CHL
!**************************************************************************
    use const_def, only: pi
    use prec_def
    implicit none    
    real(ikind),intent(in):: timehrs,rampdur
    real(ikind):: ramp    
    
    ramp = 0.5-0.5*cos(pi*min(timehrs,rampdur)/rampdur)
    
    return
    endfunction ramp_func

!********************************************************************
    subroutine index2calendar(index,iyr,imo,iday,ihr,imin,isec,ierr)
! Converts a CMS-Wave time index to a calendar date and time
! 
! Input: 
!   index - CMS-Wave index
! 
! Output:
!  iyr,imo,iday,ihr,imin,isec - Calendar date and time
! 
! Author: Alex Sanchez, USACE-CHL
!********************************************************************
    integer,intent(in) :: index
    integer,intent(out) :: iyr,imo,iday,ihr,imin,isec,ierr
    character(len=10) :: adatestr
    
    write(adatestr,'(I10)') index
    adatestr = adjustl(adatestr)
    imin = 0
    isec = 0
    ierr = 0
    if(index<9999)then !index
      ierr = -1
      iyr = index
      imo = 0; iday = 0; imin = 0; isec = 0
    elseif(index<99999)then !mddhh
      iyr = 2000   
      read(adatestr,'(I1,2I2)') imo,iday,ihr
    elseif(index<999999)then !mmddhh
      iyr = 2000
      read(adatestr,'(3I2)') imo,iday,ihr
    elseif(index<9999999)then !ymmddhh
      read(adatestr,'(I1,3I2)') iyr,imo,iday,ihr
      iyr = 2000 + iyr
    elseif(index<99999999)then !yymmddhh
      read(adatestr,'(4I2)') iyr,imo,iday,ihr
      if(iyr>2030) iyr = iyr - 100
      iyr = 2000 + iyr
    else !(ind>99999999)
      read(adatestr,'(I4,3I2)') iyr,imo,iday,ihr
    endif
    
    return
    endsubroutine index2calendar
    
!**************************************************************************
    subroutine time_cal2str(str,iyr,imo,iday,ihr,imin,isec,imilsec)
! Converts a calendar date and time to a nicely formatted string
!
! Author: Alex Sanchez, USACE-CHL    
!**************************************************************************
    use diag_lib
    implicit none
    !Input/Output
    integer, intent(in) :: iyr,imo,iday
    integer, intent(in), optional :: ihr,imin,isec,imilsec
    character(len=*), intent(inout) :: str
    !Internal
    integer :: ierr

640 format(2x,I4,'-',I2.2,'-',I2.2)    
740 format(2x,I4,'-',I2.2,'-',I2.2,1x,I2.2,':',I2.2,':',I2.2)    
840 format(2x,I4,'-',I2.2,'-',I2.2,1x,I2.2,':',I2.2,':',I2.2,'.',I3.3)
    
    if(present(imilsec))then
      write(str,840,iostat=ierr) iyr,imo,iday,ihr,imin,isec,imilsec
    elseif(present(isec))then
      write(str,740,iostat=ierr) iyr,imo,iday,ihr,imin,isec
    else
      write(str,640,iostat=ierr) iyr,imo,iday
    endif
    
    if(ierr/=0)then 
      call diag_print_warning('Could not convert calendar to string in time_cal2str()')
      str = ''
    endif
    
    return
    endsubroutine time_cal2str
    
!**************************************************************************
    subroutine time_sec2str(timesec,str)
! Converts an elapsed time in seconds to nicely formated string using
! days, hours, minutes, and seconds
!
! Author: Alex Sanchez, USACE-CHL    
!**************************************************************************
    use diag_lib
    implicit none
    !Input/Output
    real(8), intent(in) :: timesec
    character(len=*), intent(inout) :: str
    !Internal
    integer :: idays,ihrs,imin,isec,imilsec,ierr
    real(8) :: sec

850 format(I0,' days, ',I0,' hrs, ',I0,' min, ',F7.3,' s')
750 format(I0,' hrs, ',I0,' min, ',F7.3,' s')
650 format(I0,' min, ',F7.3,' s')
550 format(F7.3,' s')
    
    call time_s2dhs(timesec,idays,ihrs,imin,isec,imilsec)    
    sec = dble(isec) + dble(imilsec)*0.001d0
    
    if(idays>=1)then
      ihrs=ihrs+(idays*24)
    endif
    
    if(ihrs>=1)then
      write(str,750,iostat=ierr) ihrs,imin,sec
    elseif(imin>=1)then
      write(str,650,iostat=ierr) imin,sec
    else
      write(str,550,iostat=ierr) sec
    endif
    
    if(ierr/=0)then
      call diag_print_warning('Could not convert time to string in time_sec2str()')
      str = ''
    endif
    
    return
    endsubroutine time_sec2str    
    
!**************************************************************************
    subroutine time_s2dhs(timesec,idays,ihrs,imin,isec,imilsec)
! Converts an elapsed time in seconds to days, hours, minutes, 
! seconds, and miliseconds
!    
! written by Alex Sanchez, USACE-CHL    
!**************************************************************************
    implicit none
    !Input/Output
    real(8), intent(in) :: timesec
    integer, intent(out) :: idays,ihrs,imin,isec,imilsec    
    !Internal
    real(8) :: sec
    
    idays = int(timesec/86400.0d0)
    ihrs = int(timesec/3600.0d0-idays*24.0d0)
    imin = int(timesec/60.0d0-ihrs*60.0d0-idays*1440.0d0)
    sec = timesec - idays*86400.0d0 - ihrs*3600.0d0 - imin*60.0d0   
    isec = int(sec)
    imilsec = nint((sec-dble(isec))*1000.d0)    
    
    return
    endsubroutine time_s2dhs
    
endmodule time_lib

