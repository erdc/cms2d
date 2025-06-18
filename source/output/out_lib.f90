!==============================================================================
module out_lib
! Output Library
!
! Contains the following:
!   SMS Related
!     Time-Series Data Files
!       open_tsd    - Opens an SMS Time-Series Data file
!       append_tsd  - Appends an SMS Time-Series Data file
!       write_tsd   - Writes an SMS Time-Series Data file
!     Super ASCII
!       write_sup_file - Writes an SMS Super ASCII file *.sup
!       write_vec_dat_file - Writes a vector array to an SMS ASCII data file
!       write_scal_dat_file - Writes a scalar array to an SMS ASCII data file
!   Tecplot
!     write_tecplot_dat - Writes a tecplot data file *.dat
!     write_tecplot_his - Writes a tecplot history file *.his 
!
! written by Alex Sanchez, USACE-CHL
!            Weiming Wu, NCCHE
!            Mitchell Brown, USACE-CHL
!==============================================================================    
#include "CMS_cpp.h"
    implicit none
    
contains    

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! SMS End
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++    
    
!**************************************************************************
    subroutine write_tsd(afile,aname,imode,atype,ndat,nt,tstart,t,dat)
! Writes an SMS tsd file  
! imode = 0 - Open and write header
! imode = 1 - Open, write header and data time(s)
! imode = 2 - Open and append a data time(s)
! 
! Times are converted from seconds to hours for output.  
! Two examples of tsd files are:
!----------------------------------------------------------------------------
!TIME_SERIES
!"Current Buoy 1"  "Velocity - Mag. & Dir."  3  5  "01/01/2008 00:00:00"
!      0.0  2.0   25.0
!   3600.0  2.5   30.0
!   7200.0  3.0   35.0
!  10800.0  3.5   45.0
!  14400.0  4.6   50.0
!-------------------------------------------------------------------
!TIME_SERIES
!"Blind Pass"  "Unassigned"   2  4  "01/11/2001 07:40:00"
!      0.0   -0.539
!   1200.0   -0.506
!   2400.0   -0.451
!   3600.0   -0.377
!-------------------------------------------------------------------
! written by Alex Sanchez, USACE-CHL
!**************************************************************************
    use time_lib, only: julian2calendar
    use diag_lib
    use prec_def
    implicit none
    !Input/Output
    integer,intent(in) :: nt,ndat,imode
    character(len=*),intent(in) :: afile
    character(len=*),intent(in) :: aname,atype
    real(ikind), intent(in) :: tstart
    real(ikind), intent(in), optional :: t(:),dat(:,:)
    !Internal variables
    integer :: i,j,iyr,imo,iday,ihr,imin,isec
    character(len=19) :: astarttime
    logical :: found
        
    inquire(file=afile,exist=found)
    if(.not.found)then
      call diag_print_error('Could not open file: ',trim(afile),'  File does not exist')
    endif     

432 format(I2,'/',I2,'/',I4,1x,I2,':',I2,':',I2)  !01/01/2008 00:00:00
    call julian2calendar(tstart,iyr,imo,iday,ihr,imin,isec)
    
    !Header
    if(imode<=1)then
      open(454,file=afile)  
      write(astarttime,432) imo,iday,iyr,ihr,imin,isec
      write(454,*) aname,atype,ndat,nt,astarttime
    else
      open(454,file=afile,access='append')  
    endif
    !Append data
    if(imode>=1)then
      do i=1,max(nt,1)
        write(454,*) t(i),(dat(i,j),j=1,ndat)    
      enddo
    endif
    close(454)
    
    return
    end subroutine write_tsd
    
!**************************************************************************
    subroutine open_tsd(afile,aname,atype,ndat,nt,tstart)
! Opens an SMS Time-Series Data (*.tsd) file
! written by Alex Sanchez, USACE-CHL
!**************************************************************************
    use time_lib, only: julian2calendar
    use diag_lib
    use prec_def
    implicit none
    !Input/Output
    integer,intent(in) :: nt,ndat
    character(len=*),intent(in) :: afile
    character(len=*),intent(in) :: aname,atype
    real(ikind), intent(in) :: tstart
    !Internal variables
    integer :: iyr,imo,iday,ihr,imin,isec,ncol
    
    call julian2calendar(tstart,iyr,imo,iday,ihr,imin,isec)    
    ncol = ndat + 1
    open(454,file=afile)
    write(454,'(A)') 'TIME_SERIES' !first line of header
223 format('"',A,'"',2x,'"',A,'"',2x,I6,2x,I6,2x,&
       '"',I2.2,'/',I2.2,'/',I4,1x,I2.2,':',I2.2,':',I2.2,'"')    
    write(454,223) aname,atype,ncol,nt,imo,iday,iyr,ihr,imin,isec !second line of header
    close(454)
    
    return
    end subroutine open_tsd
    
!**************************************************************************
    subroutine append_tsd(afile,ndat,t,dat)
! Appends data to an SMS Time-Series Data (*.tsd) File.
! The file is opened and closed so it can be opened during the simulation.
! written by Alex Sanchez, USACE-CHL
! Last modified 01/09/14
!**************************************************************************
    use diag_lib
    use prec_def
    implicit none
    !Input/Output
    integer,intent(in) :: ndat
    real(ikind), intent(in) :: t,dat(ndat)
    character(len=*),intent(in) :: afile
    !Internal variables
    integer :: j
    logical :: foundfile
        
    inquire(file=afile,exist=foundfile)
    if(.not.foundfile)then
      call diag_print_error('Could not find output file: ',trim(afile))
    endif     

675 format(F13.3,1x,10000(F11.4))    
    open(454,file=afile,access='append')
    write(454,675) t,(dat(j),j=1,ndat)
    close(454)
    
    return
    end subroutine append_tsd

!***********************************************************************
    subroutine write_sup_file(aname)
! Writes a *.sup file 
! written by Alex Sanchez, USACE-ERDC-CHL, 
! Last modified April 8, 2011
!***********************************************************************
    use diag_lib, only: diag_print_error
#ifdef _WIN32
    use IFPORT
#endif    
    
    implicit none
    integer :: nunit
    character(len=*) :: aname
    character(len=200) :: filesup,filexy
    
101 format('SUPER')
102 format('SCAT2D    "',A,'"')       

    filesup = trim(aname) // '_sol.sup'    
    filexy = trim(aname) // '.xy'
    
    nunit = 46
    open(nunit,file=filesup)
    write(nunit,101)
    write(nunit,102) trim(filexy)    
    close(nunit)
    
    return
    end subroutine write_sup_file

!***********************************************************************
    subroutine write_xy_file(aname,acase)
! Writes a *.xy file 
! written by Alex Sanchez, USACE-ERDC-CHL, Lasted modified April 8, 2011
!***********************************************************************   
    use size_def
    use geo_def, only: xorigin,yorigin,azimuth_fl,x,y
    use diag_def, only: debug_mode
    use const_def, only: deg2rad
    use prec_def
    implicit none
    integer :: i,nunit,id,nc
    real(ikind) :: xtemp, ytemp, cosang, sinang
    character(len=*) :: aname,acase
    character(len=200) :: filexy
    
201 format('SCAT2D')
202 format('BEGSET')
203 format('NAME  "',A,'"')
204 format('ID',1x,I5)
205 format('DELEV 0.0')
206 format('IXY',1x,I6)
! 207 format(I8,1x,F12.5,1x,F12.5)   !for whatever reason, this was causing a run-time error - meb
! 207 format(I8,1x,F12.4,1x,F12.4)   !commented by bdj as roundoff error was printed to xy file
207 format(I8,1x,F12.3,1x,F12.3)    
208 format('ENDSET')

    filexy = trim(aname) // '.xy'

    nunit = 32; id = 17255
    if(debug_mode)then 
      nc = ncellsD !Include ghost cells
    else
      nc = ncells
    endif 
    
    open(nunit,file=filexy)
    write(nunit,201)
    write(nunit,202)
    write(nunit,203) trim(acase)
    write(nunit,204) id
    write(nunit,205)
    write(nunit,206) nc    
!    write(nunit,207) (i,x(i),y(i),i=ncells,1,-1)
    cosAng = cos(azimuth_fl*deg2rad)
    sinAng = sin(azimuth_fl*deg2rad)
    do i=nc,1,-1
     xtemp = xOrigin + x(i)*cosAng - y(i)*sinAng     !This writes the .xy files out with the wrong angle convention  !meb 01/13/2017   !switched back 05/22/2018
     ytemp = yOrigin + x(i)*sinAng + y(i)*cosAng
!      xtemp = xOrigin + x(i)*cosAng + y(i)*sinAng     !This write the file out correctly.                             !meb 01/13/2017
!      ytemp = yOrigin - x(i)*sinAng + y(i)*cosAng
      write(nunit,207) i,xtemp,ytemp
    enddo
    write(nunit,208)
    close(nunit)
    
    return
    end subroutine write_xy_file

!***********************************************************************
    subroutine write_scal_dat_file(aname,avarname,asuffix,var)
! writes scalar data to a *.dat file 
! written by Alex Sanchez, USACE-ERDC-CHL, Nov. 17, 2008
!***********************************************************************  
    use size_def
    use comvarbl, only: timehrs,reftime
    use diag_def, only: debug_mode
    use prec_def
    
    implicit none  
    !Input/Output
    character(len=*),intent(in) :: aname,avarname,asuffix
    real(ikind),intent(in) :: var(ncellsD)
    !Internal Variables
    integer :: i,nunit,id,nc,k,ierr
    character(len=200) :: filedat,filesup,aline,cardname
    logical :: isnankind, found, foundfile
    logical :: isHotStart=.false.

103 format('DATA    "',A,'"')
301 format('DATASET')
302 format('OBJTYPE "scat2d"')
304 format('BEGSCL')
305 format('OBJID ',I5)
306 format('ND',1x,I6)
307 format('NC',1x,I6)
308 format('NAME',1x,'"',A,'"')
309 format('RT_JULIAN ',F15.3) ! 18250.0
310 format('TIMEUNITS hours')
311 format('TS',1x,I1,1x,F12.4) 
312 format(F12.0) 
313 format(1pe12.5) 
314 format('ENDDS')    

    filedat = trim(aname) // '_' // trim(asuffix) // '.dat'

    !Test to see if this is hotstart dat file.
    i=index(filedat,'HotStart')
    if (i .gt. 0) then
      isHotStart = .true.
      !if HotStart file, delete file previously written to'
      open(100,file=filedat)
      close(100,status='DELETE')
    endif  
    
    nunit = 124
    id = 17255
    if(debug_mode)then 
      nc = ncellsD !Include ghost cells
    else
      nc = ncells
    endif  
    
    !Testing new way
    inquire(file=filedat,exist=found)
    if(found) then
      !File exists. but this is a cold start and at the beginning. Delete the file and start over so you don't get duplicate records.  05/21/2018
      if (timehrs .eq. 0.d0) then
        open(nunit,file=filedat)
        close(nunit,status='DELETE')
        found = .false.  !reset this so it goes through the next section properly.
      endif
    endif
    
    if(.not.found) then   !file doesn't exist, so write header information
    !if(timehrs<1.e-6)then       !Old way
      open(nunit,file=filedat,status='unknown',action='write')
      write(nunit,301)
      write(nunit,302)  
      write(nunit,304)
      write(nunit,305) id
      write(nunit,306) nc
      write(nunit,307) nc      
      write(nunit,308) trim(avarname)
      write(nunit,309) reftime
      write(nunit,310)   
      
      if(isHotStart) then
        filesup = trim(aname) // '.sup'  
      else
        filesup = trim(aname) // '_sol.sup'  
      endif

      !Only add the filename the first time.
      open(46,file=filesup,status='old')
      foundfile=.false.
      do
        read(46,'(A100)',iostat=ierr) aline
        if(ierr/=0) exit
        i = index(aline,trim(filedat))
        if (i.gt.0) then
          foundfile=.true.
          exit
        endif
      enddo
      close(46)
      i = 0
      if(.not.foundfile) then
        open(46,file=filesup,status='old',position='append')
        !Ensure name written has the 'ASCII_Solutions' prefix like the others  MEB 9/30/2021
        i = index(filedat,'ASCII_Solutions')
        if (i == 0) then
          write(46,103) aline(10:25)//trim(filedat)
        else
          write(46,103) trim(filedat)
        endif
        close(46)
      endif  
    else
      open(nunit,file=filedat,status='old',position='append')
      backspace(nunit)
    endif       
    write(nunit,311) 0,timehrs    
    do i=1,nc
      if(isnankind(var(i)))then
        write(nunit,312) -999.0
      else
        write(nunit,313) var(i)
      endif
    enddo   
    write(nunit,314)
    close(nunit)    
    
    return
    end subroutine write_scal_dat_file
    
!***********************************************************************
    subroutine sms_write_dat_scal(aprefix,asuffix,aname, RT_Julian,time,atimeunits,nd,nc,ns,var,stat)
! writes scalar data to a *.dat file 
!
! Description:
!  Writes an SMS ASCII Data File (*.dat)
!
! Input:
!  aprefix - File name prefix including path
!  asuffix - File name suffix without extension
!  aname - Name of the data set 
!  RT_Julian - Reference time as a Julian number
!  time - Time with respect to the reference time
!  atimeunits -  time units (e.g. hours, minutes, seconds)
!  nd - Number of data values that will be listed per time step. 
!       This number should correspond to the total number of vertices, 
!       nodes, cells centers (cell-centered grid), cell corners 
!       (mesh-centered grid), maximum node id (meshes) or scatter points. 
!  nc - Maximum element id (meshes) or the number of cells (grids).
!  var - Output scalar field of size (1:nd)
!  stat - Status variable. Optional. 1 - Active, 0 - Inactive
!
! Output:
!   SMS ASCII data file
!
! Author: Alex Sanchez, USACE-ERDC-CHL
!***********************************************************************
    use diag_def, only: debug_mode
    use prec_def
    implicit none  
    !Input/Output
    character(len=*),intent(in) :: aprefix,asuffix,aname
    real(8),intent(in) :: RT_Julian,time
    character(len=*),intent(in) :: atimeunits
    integer,intent(in) :: nd,nc,ns
    real(ikind),intent(in) :: var(ns)
    integer,intent(in),optional :: stat(ns)
    !Internal Variables
    integer :: i,nunit,id
    character(len=200) :: filedat,filesup
    logical :: isnankind

103 format('DATA    "',A,'"')
301 format('DATASET')
302 format('OBJTYPE "scat2d"')
304 format('BEGSCL')
305 format('OBJID ',I5)
306 format('ND',1x,I6)
307 format('NC',1x,I6)
308 format('NAME',1x,'"',A,'"')
309 format('RT_JULIAN ',F15.3) ! 18250.0
310 format('TIMEUNITS',1x,A)
311 format('TS',1x,I1,1x,F12.4) 
312 format(F12.0)
313 format(1pe12.5) 
314 format('ENDDS') 
315 format(I1) 

    filedat = trim(aname) // '_' // trim(asuffix) // '.dat'
    
    nunit = 124
    id = 17255 
    if(time<1.0e-6)then !First time step
      open(nunit,file=filedat,status='unknown',action='write')
      write(nunit,301)
      write(nunit,302)  
      write(nunit,304)
      write(nunit,305) id
      write(nunit,306) nc
      write(nunit,307) nc      
      write(nunit,308) trim(aname)
      write(nunit,309) RT_Julian
      write(nunit,310)
      filesup = trim(aprefix) // '_sol.sup'  
      open(46,file=filesup,status='old',position='append')
      write(46,103) trim(filedat)
      close(46)
    else
      open(nunit,file=filedat,status='old',position='append')
      backspace(nunit)
    endif         
    if(present(stat))then
      write(nunit,311) 1,time
      do i=1,nd
        write(nunit,315) stat(i)
      enddo 
    else
      write(nunit,311) 0,time   
    endif
    do i=1,nd
      if(isnankind(var(i)))then
        write(nunit,312) -999.0
      else
        write(nunit,313) var(i)
      endif
    enddo   
    write(nunit,314)
    close(nunit)    
    
    return
    end subroutine sms_write_dat_scal

!***********************************************************************
    subroutine write_vec_dat_file(aname,avarname,asuffix,varx,vary)
! writes vector data to a *.dat file 
! written by Alex Sanchez, USACE-ERDC-CHL, Nov. 17, 2008
!***********************************************************************  
    use size_def
    use geo_def, only: azimuth_fl    
    use comvarbl, only: timehrs,reftime
    use diag_def, only: debug_mode
    use const_def, only: deg2rad
    use prec_def
    implicit none 
    !Input/Output
    character(len=*),intent(in) :: aname,avarname,asuffix
    real(ikind),intent(in) :: varx(ncellsD),vary(ncellsD)
    !Internal Variables
    integer :: i,nunit,id,nc,ierr,k
    real(ikind) :: cosAng,sinAng,xvec,yvec
    character(len=200) :: filedat,filesup,aline,cardname
    logical :: isnankind, found, foundfile
    logical :: isHotStart = .false.

103 format('DATA    "',A,'"')
301 format('DATASET')
302 format('OBJTYPE "scat2d"')
304 format('BEGVEC')
303 format('VECTYPE 0')     !missing from file format - added 031618 MEB
305 format('OBJID',1x,I5)
306 format('ND',1x,I6)
307 format('NC',1x,I6)
308 format('NAME',1x,'"',A,'"')
309 format('RT_JULIAN ',F15.3) ! 18250.0
310 format('TIMEUNITS hours')
311 format('TS',1x,I1,1x,F12.4) 
312 format(2F12.0)
313 format(1pe12.5,1x,1pe12.5)
314 format('ENDDS')    
  
    filedat = trim(aname) // '_' // trim(asuffix) // '.dat'

    !Test to see if this is hotstart dat file.
    i=index(filedat,'HotStart')
    if (i .gt. 0) then
      isHotStart = .true.

      !if HotStart file, delete file previously written to'
      open(100,file=filedat)
      close(100,status='DELETE')
    endif  

    nunit = 124
    id = 17255
    if(debug_mode)then 
      nc = ncellsD !Include ghost cells
    else
      nc = ncells
    endif 

    !Testing new way
    inquire(file=filedat,exist=found)
    if(found) then
      !File exists. but this is a cold start and at the beginning. Delete the file and start over so you don't get duplicate records.  05/21/2018
      if (timehrs .eq. 0.d0) then
        open(nunit,file=filedat)
        close(nunit,status='DELETE')
        found = .false.  !reset this so it goes through the next section properly.
      endif
    endif
   
    if(.not.found) then   !file doesn't exist, so write header information
    !if(timehrs<1.e-6)then
      open(nunit,file=filedat,status='unknown',action='write')
      write(nunit,301)
      write(nunit,302)  
      write(nunit,304)
      write(nunit,303)
      write(nunit,305) id
      write(nunit,306) nc
      write(nunit,307) nc
      write(nunit,308) trim(avarname)
      write(nunit,309) reftime
      write(nunit,310)
      if(isHotStart) then
        filesup = trim(aname) // '.sup'  
      else
        filesup = trim(aname) // '_sol.sup'  
      endif
      !Search for filename and only add the filename the first time
      open(46,file=filesup,status='old')
      foundfile=.false.
      do
        read(46,'(A100)',iostat=ierr) aline
        if(ierr/=0) exit
        i = index(aline,trim(filedat))
        if (i.gt.0) then
          foundfile=.true.
          exit
        endif
      enddo
      close(46)
      if(.not.foundfile) then
        open(46,file=filesup,status='old',position='append')
        write(46,103) trim(filedat)
        close(46)
      endif  
    else
      open(nunit,file=filedat,status='old',position='append')
      backspace(nunit)
    endif     
    write(nunit,311) 0,timehrs    
    sinAng = sin(azimuth_fl*deg2rad)
    cosAng = cos(azimuth_fl*deg2rad)
    do i=1,nc
      if(isnankind(varx(i)))then
        write(nunit,312) -999.0,-999.0
      else
        xvec = varx(i)*cosAng - vary(i)*sinAng
        yvec = varx(i)*sinAng + vary(i)*cosAng
        write(nunit,313) xvec,yvec
      endif
    enddo
    write(nunit,314)
    close(nunit) 
    
    return
    end subroutine write_vec_dat_file

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! SMS End
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++    

!**************************************************************************
    subroutine writescalTxt(afile,apath,dname,var,iwritedry)
! writes a scalar dataset to a TXT dataset file with id number PID
! - APATH differs in this from the XMDF format to allow for a subdirectory
!   to hold the many statistics ASCII files in one location.
! written by Mitchell Brown, USACE-ERDC-CHL  
! afile = filename
! apath = subdirectory to place files (will create if it doesn't exist)
! dname = dataset name for this variable
!**************************************************************************
    use size_def, only: ncellsD,ncellsfull,ncellpoly
    use interp_lib, only: interp_scal_cell2node
    use prec_def
    use diag_lib, only: diag_print_error
#ifdef _WIN32    
    use IFPORT
#endif
    implicit none
    !Input/Output
    character(len=*),intent(in) :: afile,apath,dname
    real(ikind),intent(in) :: var(ncellsD)
    integer,intent(in) :: iwritedry
    !Internal variables
    integer :: kunit,j
    real(8) :: timed  !Must be double    
    real(4) :: scalout(ncellsfull) !Must be single
    character(len=200) :: afullpath, thefile, ppath, pname
    character(len=10)  :: pext
    logical :: found, created
    
    call fileparts(afile,ppath,pname,pext)  !Break 'afile' into <directory path>, <name>, and <extension>.
    
    kunit = 601
    if(ncellpoly>0)then
      call interp_scal_cell2node(var,scalout,iwritedry)
    else
      call map_scal_active2full(var,scalout,iwritedry)  
    endif
    
!    afullpath=trim(apath)//trim(aname)    
    if (apath==' ') then
      thefile=trim(pname)//"_"//trim(dname)//'.txt'  !Don't create a directory and put into the current directory path
    else
#ifdef _WIN32        
      inquire(directory=trim(apath), exist=found)
      if (.not.found) then 
        created=MakeDirQQ(trim(apath))
        if (.not.created) then
          call diag_print_error('Failed to create subdirectory- '//trim(apath))
        endif
      endif
#else
      inquire(file=trim(apath), exist=found)
      if (.not.found) then 
        call system('mkdir '//(trim(apath)))
      endif
#endif
      thefile=trim(apath)//'/'//trim(pname)//"_"//trim(dname)//'.txt'
    endif  

200 FORMAT (I0,x,I1,5x,'!Number of cells, Number of dimensions')    
201 FORMAT (25(E10.3,2x))
    
    open(kunit,FILE=thefile)
    write(kunit,200) ncellsfull, 1  !Always 1 for Scalar
    write(kunit,201) (scalout(j),j=1,ncellsfull)
    close(kunit)
    
    return
    end subroutine writescalTxt
    
!**************************************************************************
    subroutine writevecTxt(afile,apath,dname,varx,vary,iwritedry)
! writes a vector dataset to the TXT dataset file with id ncellsfull PID
!
! written by Mitchell Brown, USACE-ERDC-CHL  
! afile = filename
! apath = subdirectory to place files (will create if it doesn't exist)
! dname = dataset name for this variable    
!**************************************************************************
    use size_def, only: ncellsD,ncellsfull,ncellpoly
    use diag_lib, only: diag_print_error
    use interp_lib, only: interp_vec_cell2node
    use prec_def
#ifdef _WIN32
    use IFPORT
#endif
    implicit none
    !Input/Output    
    character(len=*),intent(in) :: afile,apath,dname
    real(ikind),     intent(in) :: varx(ncellsD),vary(ncellsD)
    integer,         intent(in) :: iwritedry    
    !Internal Variables
    integer :: kunit,j
    real*8 :: timed !Must be double for XMDF subs    
    real*4 :: vecout(ncellsfull*2) !Must be single 
    character(len=200) :: afullpath, thefile, ppath, pname
    character(len=10)  :: pext
    logical :: found, created
    
    call fileparts(afile,ppath,pname,pext)  !Break 'afile' into <directory path>, <name>, and <extension>.
    
    kunit = 602
    if(ncellpoly>0)then
      call interp_vec_cell2node(varx,vary,vecout,iwritedry)
    else
      call map_vec_active2full(varx,vary,vecout,iwritedry)  
    endif
    
!    afullpath=trim(apath)//trim(aname)    
    if (apath==' ') then
      thefile=trim(pname)//"_"//trim(dname)//'.txt'  !Don't create a directory and put into the current directory path
    else
#ifdef _WIN32        
      inquire(directory=trim(apath), exist=found)
      if (.not.found) then 
        created=MakeDirQQ(trim(apath))
        if (.not.created) then
          call diag_print_error('Failed to create subdirectory- '//trim(apath))
        endif
      endif
#else
      inquire(file=trim(apath), exist=found)
      if (.not.found) then 
        call system('mkdir '//(trim(apath)))
      endif
#endif
      thefile=trim(apath)//'/'//trim(pname)//"_"//trim(dname)//'.txt'
    endif  

200 FORMAT (I0,x,I1,5x,'!Number of cells, Number of dimensions')    
202 FORMAT (25(E10.3,x,E10.3,2x))
    
    open(kunit,FILE=thefile)
    write(kunit,200) ncellsfull, 2  !Always 2 for Vector
    write(kunit,202) (vecout(j),j=1,ncellsfull*2)
    close(kunit)    
          
    return
    end subroutine writevecTxt
    
!**************************************************************************
    subroutine write_xys (id,apath,outprefix,aprefix,ntimes,times,val)
    
    implicit none
    integer            :: id,ntimes,i
    character(len=*)   :: apath,aprefix,outprefix
    character(len=200) :: xysoutfile,abnd,astring,astring2
    real(4)            :: times(ntimes),val(ntimes)
    
    write(abnd,'(I0)') id
    xysoutfile = trim(apath)// '/' //trim(outprefix) // trim(aprefix) // trim(abnd) // '.xys'
    open(913,file=xysoutfile)
    write(913,'(A3,x,i0,x,i0,x,A1,A1)') 'XYS',2,ntimes,'"','"'
    do i=1,ntimes
      if (val(i) .eq. -0.0) val(i) = 0.0
      write(astring, '(F10.2)') times(i)
      write(astring2,'(F10.2)') val(i)
      write(913,'(A)') trim(adjustl(astring))//' '//trim(adjustl(astring2))
    enddo
    close(913)
    
    return
    end subroutine write_xys
    
    
#ifdef XMDF_IO
!**************************************************************************
    subroutine writescalh5(afile,apath,aname,var,aunits,timehr,iwritedry,writecf_var)
! writes a scalar dataset to the xmdf file with id number PID
!
! written by Alex Sanchez, USACE-ERDC-CHL  
! added optional argument, Mitchell Brown, USACE-ERDC-CHL  10/02/2023
!**************************************************************************
    use xmdf
    use prec_def
    use size_def, only: ncellsD,ncellsfull,ncellpoly
    use const_def, only: READWRITE
    use interp_lib, only: interp_scal_cell2node

    implicit none
    
    !Input/Output
    character(len=*),intent(in) :: afile,apath,aname,aunits
    real(ikind),intent(in) :: var(ncellsD),timehr
    integer,intent(in) :: iwritedry
    logical,intent(in), optional :: writecf_var   !This logical indicates whether to pull from the 'cf_vals' list and writes extra attributes to the dataset.
    
    
    !Internal variables
    integer(XID) :: fid,gid
    integer ::ierr    
    real(8) :: timed  !Must be double    
    real(4) :: scalout(ncellsfull) !Must be single
    character(len=200) :: afullpath
    
    if(ncellpoly>0)then
      call interp_scal_cell2node(var,scalout,iwritedry)
    else
      call map_scal_active2full(var,scalout,iwritedry)  
    endif
    
    afullpath=trim(apath)//trim(aname)    
    call XF_OPEN_FILE(trim(afile),readwrite,fid,ierr) !Open XMDF file
    if(ierr<0)then
      call XF_CREATE_FILE(trim(afile),readwrite,fid,ierr)
    endif                         
    call OPEN_CREATE_DATASET(fid,trim(afullpath),gid,1,aunits,ierr)  !Open/create dataset
    if(present(writecf_var) .and. ierr == -666) call write_h5_general_cf_info(gid, trim(aname)) !Write general CF attributes on new file creation  MEB 02/21/2025

    timed = dble(timehr)
    call XF_WRITE_SCALAR_TIMESTEP(gid,timed,ncellsfull,scalout,ierr) !Write data to XMDF file
    call XF_CLOSE_GROUP(gid,ierr)  !Close dataset    
    call XF_CLOSE_FILE(fid,ierr)   !Close XMDF file
    
    return
    end subroutine writescalh5
    
!**************************************************************************
    subroutine write_h5_general_cf_info (a_Id, aname)
! when a new XMDF file is create, writes the appropriate CF compliance attributes
!
! Mitchell Brown, USACE-ERDC-CHL  02/21/2025
!**************************************************************************
    use xmdf,     only: XF_SET_ATTRIBUTE_STRING, XF_CLOSE_GROUP
    use XMDFDEFS, only: XID
    use out_def,  only: cf_vars,ncf_vars
    use diag_lib, only: diag_print_warning
    implicit none
    
    integer(XID), intent(in) :: a_Id
    character(len=*), intent(in) :: aname
    
    integer(XID) :: pid
    integer             :: error, i, list_item
    character(len=50)   :: output_name
    character(len=100)  :: long_name
    character(len=100)  :: standard_name
    character(len=20)   :: units
    character(len=10)   :: positive
    
    list_item = -1
    do i=1,ncf_vars
      if (trim(cf_vars(i)%output_name) == aname) then
        list_item = i
        exit
      endif
    enddo
    
    if (list_item < 0) then
      call diag_print_warning('Could not find CF match to: '//trim(aname))
      return
    endif
    
    output_name = cf_vars(list_item)%output_name
    long_name = cf_vars(list_item)%long_name
    standard_name = cf_vars(list_item)%standard_name
    units = cf_vars(list_item)%units
    positive = cf_vars(list_item)%positive
    
    call OPEN_CREATE_DATASET(a_Id,'PROPERTIES',pid,1,'',error)                   !Open dataset (this dataset should already have been created).
    call XF_SET_ATTRIBUTE_STRING (pid, 'short_name', output_name, error)
    call XF_SET_ATTRIBUTE_STRING (pid, 'long_name', long_name, error)
    call XF_SET_ATTRIBUTE_STRING (pid, 'standard_name', standard_name, error)
    call XF_SET_ATTRIBUTE_STRING (pid, 'units', units, error)
    call XF_SET_ATTRIBUTE_STRING (pid, 'positive', positive, error)
    call XF_CLOSE_GROUP(pid,error)  !Close dataset    
    
    return
    end subroutine write_h5_general_cf_info
    
    
!**************************************************************************
    subroutine writevech5(afile,apath,aname,varx,vary,aunits,timehr,iwritedry,writecf_var)
! writes a vector dataset to the xmdf file with id ncellsfull PID
!
! written by Alex Sanchez, USACE-ERDC-CHL  
!**************************************************************************
    use size_def, only: ncellsD,ncellsfull,ncellpoly
    use xmdf
    use interp_lib, only: interp_vec_cell2node
    use prec_def
    use const_def, only: READWRITE
    implicit none
    
    !Input/Output    
    character(len=*),intent(in)  :: afile,apath,aname,aunits
    real(ikind),     intent(in)  :: varx(ncellsD),vary(ncellsD),timehr
    integer,         intent(in)  :: iwritedry
    logical,intent(in), optional :: writecf_var   !This logical indicates whether to pull from the 'cf_vals' list and writes extra attributes to the dataset.

    !Internal Variables
    integer(XID) :: pid,did
    integer :: ierr
    real*8 :: timed !Must be double for XMDF subs    
    real*4 :: vecout(ncellsfull*2) !Must be single 
    character :: afullpath*100    
    
    if(ncellpoly>0)then
      call interp_vec_cell2node(varx,vary,vecout,iwritedry)
    else
      call map_vec_active2full(varx,vary,vecout,iwritedry)  
    endif
    
    afullpath=trim(apath)//trim(aname)    
    call XF_OPEN_FILE(trim(afile),readwrite,pid,ierr) !Open XMDF file
    if(ierr<0)then
      call XF_CREATE_FILE(trim(afile),readwrite,pid,ierr)
    endif 
    call OPEN_CREATE_DATASET(PID,trim(afullpath),did,2,aunits,ierr) !Open/create dataset
    if(present(writecf_var) .and. ierr == -666) call write_h5_general_cf_info(did, trim(aname)) !Write general CF attributes on new file creation  MEB 02/21/2025

    timed = dble(timehr)
    call XF_WRITE_VECTOR_TIMESTEP(did,timed,ncellsfull,2,vecout,ierr) !Write data to XMDF file    
    call XF_CLOSE_GROUP(did,ierr) !Close dataset   
    call XF_CLOSE_FILE(pid,ierr) !Close XMDF file      
          
    return
    end subroutine writevech5
    
!**************************************************************************
    SUBROUTINE OPEN_CREATE_DATASET(OCID,STRING1,OCDID,DIM,OUNITS,OCERR)
! THIS FUNCTION OPENS OR CREATES A SCALAR OR VECTOR DATASET AND 
!   EXITS IF IT CAN'T CREATE ONE.
!
! - OCID:     parent group/file id
! - STRING1:  name for dataset
! - OCDID:    child group id
! - DIM:      dataset dimensions (1- scalar, 2- vector)
! - OUNITS:   Data units ('m','m/s','kg/m^3',etc.)
! - REFTIME (in):  SMS Reference time for dataset
!
! written by MEB 04/01/07
! Updates by Alex Sanchez 7/27/2010
!        - fixed time units to hours 
!        - added compression option
!        - added output units for dataset
!**************************************************************************
    use comvarbl, only: reftime
    use out_def, only: ixmdfcomp
    use diag_lib
    use XMDF
    implicit none
    
    integer(XID), intent(in) :: OCID !Parent group/file id
    integer(XID), intent(out) :: OCDID  !Child group id
    integer, intent(in)    :: DIM
    integer, intent(inout) :: OCERR
    character(LEN=*), intent(in) :: STRING1,OUNITS    

    integer(XID)          :: OCDPID  !Property group id

      
    SELECT CASE (DIM)
      CASE (1)
        CALL XF_OPEN_GROUP(OCID,STRING1,OCDID,OCERR)
        IF(OCERR<0)THEN
          !Note: Time units always in hours (hard-wired from SMS)
          CALL XF_CREATE_SCALAR_DATASET(OCID,STRING1,OUNITS,TS_HOURS,ixmdfcomp,OCDID,OCERR)
          IF(OCERR<=0)THEN
            call diag_print_error('Could not create scalar dataset: ',STRING1)            
          ENDIF
          CALL XF_DATASET_REFTIME(OCDID,REFTIME,OCERR)
          CALL XF_CREATE_PROPERTY_GROUP(OCDID,OCDPID,OCERR)
          CALL XF_WRITE_PROPERTY_FLOAT(OCDPID,PROP_NULL_VALUE,1,-999.0,NONE,OCERR)
          CALL XF_SCALAR_DATA_LOCATION(OCDID,GRID_LOC_CENTER,OCERR)
          OCERR = -666  !Something specific to test if the dataset was newly created  !MEB 02/24/2025
        ENDIF
      
      CASE (2)
        CALL XF_OPEN_GROUP(OCID,STRING1,OCDID,OCERR)
        IF(OCERR<0)THEN
          !Note: Time units always in hours (hard-wired from SMS)
          CALL XF_CREATE_VECTOR_DATASET(OCID,STRING1,OUNITS,TS_HOURS,ixmdfcomp,OCDID,OCERR) 
          IF(OCERR<=0)THEN
            call diag_print_error('Could not create vector dataset: ',STRING1)   
          ENDIF
          CALL XF_DATASET_REFTIME(OCDID,REFTIME,OCERR)
          CALL XF_CREATE_PROPERTY_GROUP(OCDID,OCDPID,OCERR)
          CALL XF_WRITE_PROPERTY_FLOAT(OCDPID,PROP_NULL_VALUE,1,-999.0,NONE,OCERR)
          CALL XF_VECTORS_IN_LOCAL_COORDS(OCDID,OCERR)
          CALL XF_VECTOR_2D_DATA_LOCS(OCDID,GRID_LOC_CENTER,GRID_LOC_CENTER,OCERR)
          OCERR = -666  !Something specific to test if the dataset was newly created  !MEB 02/24/2025
        ENDIF
      
      END SELECT
      
     RETURN
    END SUBROUTINE OPEN_CREATE_DATASET
#endif


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! XMDF End
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Tecplot Begin
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            
!***********************************************************************
    subroutine write_tecplot_dat
!   write results for tecplot                                                                
!   by Weiming Wu, NCCHE, Oct. 2008
!***********************************************************************
    use size_def
    use geo_def
    use flow_def
    use fric_def, only: bsxy
    use out_def, only: datfile
    use comvarbl
    use sed_def
    use sal_def
    use heat_def
    use wave_flowgrid_def
    use met_def, only: windconst,windvar,tauwindx,tauwindy,presvar,pressatm
    use cms_def
    implicit none
    integer :: i, ictec
      
    ictec = 31

    open (unit=ictec,file=datfile)
    rewind ictec

    write (ictec,*) 'title= "CMS2D Results"'
    write (ictec,*) 'variables = "x","y","u","v",'
    write (ictec,*) '"uv","h","zs","zb","taub","vis","iwet"'
    if(sedtrans) write (ictec,*) ',"sed"'
    if(saltrans) write (ictec,*) ',"sal"'
    if(heattrans) write (ictec,*) ',"temperature"'
    if(noptset==3) write (ictec,*) ',"wavestrx","wavestry","waveh","wavet","waveob","waveangle","wavediss"'
    if(windconst .or. windvar) write(ictec,*) ',"tauwindx","tauwindy"'
    if(presvar) write(ictec,*) ',"pressatm"'
    write (ictec,*) 'zone f=block,'
    write (ictec,*) 'i =',ncellsD

    write (ictec,*) (x(i),i=1,ncellsD)
    write (ictec,*) (y(i),i=1,ncellsD)
    write (ictec,*) (u(i),i=1,ncellsD)
    write (ictec,*) (v(i),i=1,ncellsD)
    write (ictec,*) (uv(i),i=1,ncellsD)
    write (ictec,*) (h(i),i=1,ncellsD)
    write (ictec,*) (p(i)*gravinv,i=1,ncellsD)
    write (ictec,*) (zb(i),i=1,ncellsD)
    write (ictec,*) (bsxy(i),i=1,ncellsD)
    write (ictec,*) (vis(i),i=1,ncellsD)
    write (ictec,*) (iwet(i),i=1,ncellsD)
    if(sedtrans) write (ictec,*) (ct(i),i=1,ncellsD)
    if(saltrans) write (ictec,*) (sal(i),i=1,ncellsD)
    if(heattrans) write (ictec,*) (heat(i),i=1,ncellsD)
      if(noptset==3) then
      write (ictec,*) (wavestrx(i),i=1,ncellsD)
      write (ictec,*) (wavestry(i),i=1,ncellsD)
      write (ictec,*) (Whgt(i),i=1,ncellsD)
      write (ictec,*) (Wper(i),i=1,ncellsD)
      write (ictec,*) (Worb(i),i=1,ncellsD)
      write (ictec,*) (Wang(i),i=1,ncellsD)
      write (ictec,*) (wavediss(i),i=1,ncellsD)
    endif
    if(windconst .or. windvar) then
      write (ictec,*) (tauwindx(i),i=1,ncellsD)
      write (ictec,*) (tauwindy(i),i=1,ncellsD)
    endif
    if(presvar) write (ictec,*) (pressatm(i),i=1,ncellsD)
    close(ictec)

    return
    end subroutine write_tecplot_dat
    
!***********************************************************************
    subroutine write_tecplot_his
!   write history results for tecplot                                                                
!   by Weiming Wu, NCCHE, Oct. 2008
!***********************************************************************
    use size_def
    use geo_def
    use flow_def
    use fric_def, only: bsxy
    use struct_def
    use comvarbl
    use sed_def
    use sal_def
    use heat_def
    use wave_flowgrid_def
    use met_def, only: windconst,windvar,tauwindx,tauwindy,presvar,pressatm
    use cms_def
    use prec_def
    implicit none
    integer :: i,ictec
            
    ictec = 30

    write (ictec,*) 'title= "CMS2D Results"'
    write (ictec,*) 'variables = "x","y","u","v",'
    write (ictec,*) '"uv","h","zs","zb","taub","vis","iwet"'
    if(sedtrans) write (ictec,*) ',"sed"'
    if(saltrans) write (ictec,*) ',"sal"'
    if(heattrans) write (ictec,*) ',"temperature"'
    if(noptset==3) write (ictec,*) ',"wavestrx","wavestry","waveh","wavet","waveob","wavangle","wavediss"'
    if(windconst .or. windvar) write(ictec,*) ',"tauwindx","tauwindy"'
    if(presvar) write(ictec,*) ',"pressatm"'
    write (ictec,*) 'zone t="',ctime/3600.0, '",    f=block,'
    write (ictec,*) 'i =',ncellsD

    write (ictec,*) (x(i),i=1,ncellsD)
    write (ictec,*) (y(i),i=1,ncellsD)
    write (ictec,*) (u(i),i=1,ncellsD)
    write (ictec,*) (v(i),i=1,ncellsD)
    write (ictec,*) (uv(i),i=1,ncellsD)
    write (ictec,*) (h(i),i=1,ncellsD)
    write (ictec,*) (p(i)*gravinv,i=1,ncellsD)
    write (ictec,*) (zb(i),i=1,ncellsD)
    write (ictec,*) (bsxy(i),i=1,ncellsD)
    write (ictec,*) (vis(i),i=1,ncellsD)
    write (ictec,*) (iwet(i),i=1,ncellsD)
    if(sedtrans) write (ictec,*) (ct(i),i=1,ncellsD)
    if(saltrans) write (ictec,*) (sal(i),i=1,ncellsD)
    if(heattrans) write (ictec,*) (heat(i),i=1,ncellsD)
    if(noptset==3) then
      write (ictec,*) (wavestrx(i),i=1,ncellsD)
      write (ictec,*) (wavestry(i),i=1,ncellsD)
      write (ictec,*) (Whgt(i),i=1,ncellsD)
      write (ictec,*) (Wper(i),i=1,ncellsD)
      write (ictec,*) (Worb(i),i=1,ncellsD)
      write (ictec,*) (Wang(i),i=1,ncellsD)
      write (ictec,*) (wavediss(i),i=1,ncellsD)
    endif
    if(windconst .or. windvar) then
      write (ictec,*) (tauwindx(i),i=1,ncellsD)
      write (ictec,*) (tauwindy(i),i=1,ncellsD)
    endif
    if(presvar) write (ictec,*) (pressatm(i),i=1,ncellsD)

    return
    end subroutine write_tecplot_his

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Tecplot End
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++    

!*******************************************************************
    subroutine add_item_to_cf_list(name, long, standard, units, positive)
! Adds items to the cf compliance list of variables
! written by Mitchell Brown, USACE-CHL  10/02/2023
!*******************************************************************
    use out_def, only: cf_vars,cf_var_type,ncf_vars
    implicit none
    
    integer icf
    
    character(len=*),intent(in) :: name, long, standard, units, positive
    type(cf_var_type), allocatable :: temp(:)
    
    ncf_vars = ncf_vars + 1
    if(ncf_vars == 1) then              !Add first element to list
      allocate(cf_vars(1))              
    else                                !Extend list, move old elements to larger list using a temp.
      allocate(temp(ncf_vars - 1))      
      do icf=1,ncf_vars - 1
        temp(icf) = cf_vars(icf)
      enddo
      deallocate(cf_vars)
      allocate(cf_vars(ncf_vars))
      do icf=1,ncf_vars-1
        cf_vars(icf) = temp(icf)
      enddo
      deallocate(temp)
    endif
    
    cf_vars(ncf_vars)%units = units
    cf_vars(ncf_vars)%long_name = long
    cf_vars(ncf_vars)%output_name = name
    cf_vars(ncf_vars)%positive = positive
    cf_vars(ncf_vars)%standard_name = standard
    
    return
    end subroutine add_item_to_cf_list
    
!*******************************************************************
    subroutine init_cf_var_list()
! Initializes the cf compliance list of variables
! written by Mitchell Brown, USACE-CHL  10/02/2023
!*******************************************************************
    use out_def, only: cf_vars, ncf_vars
    implicit none
    
    ncf_vars = 0
    call add_item_to_cf_list('Water_Elevation','sea surface elevation','surface_elevation','m','up')
    call add_item_to_cf_list('Total_Water_Depth','bathymetry plus surface elevation','total_water_depth','m','down')
    call add_item_to_cf_list('Current_Velocity','vertical averaged velocity vector','u/v_velocity','m/s','x/y_direction')
    call add_item_to_cf_list('Current_Magnitude','magnitude of vertical averaged velocity','velocity_magnitude','m/s','up')
    call add_item_to_cf_list('Depth','bathymetry plus change due to sediment transport','depth_through_time','m','down')
    call add_item_to_cf_list('Morphology_Change','change due to sediment transport','depth_change_through_time','m','up')
    call add_item_to_cf_list('Eddy_Viscosity','diffusivity due to eddy advection','ocean_tracer_diffusivity_due_to_parameterized_mesoscale_eddy_advection','m^2/s','up')
    call add_item_to_cf_list('Concentration','concentration of suspended sediment','mass_concentration_of_suspended_sediment_in_sea_water','kg/m^3','up')
    call add_item_to_cf_list('Capacity','maximum capacity of sea water to hold sediment','maximum_capacity_of_sea_water_to_hold_suspended_sediment','kg/m^3','up')
    call add_item_to_cf_list('Total_Sediment_Transport','total sediment transport across unit distance in ocean', &
                             'total_sediment_transport_across_unit_distance_in_ocean','kg/m/s','x/y_direction')
    call add_item_to_cf_list('Fraction_Suspended','fraction of suspended sediment of specific grain size', &
                             'fraction_suspended_sediment_for_specific_grain_size_in_sea_water','nondimensional','up')
    call add_item_to_cf_list('Salinity','concentration of salinity','mass_concentration_of_salinity_in_sea_water','ppt','up')
    call add_item_to_cf_list('Temperature','sea water conservative temperature','sea_water_redistributed_conservative_temperature','degree C','up')
    call add_item_to_cf_list('Wave_Height','significant wave height','sea_surface_wave_significant_height','m','up')
    call add_item_to_cf_list('Wave_Period','parabolic fit to the period of the peak of the energy', &
                             'sea_surface_wave_period_at_variance_spectral_density_maximum','s','up')
    call add_item_to_cf_list('Wave_Height_Vec','significant wave height vectors','sea_surface_wave_significant_height_vectors','m','up')
    call add_item_to_cf_list('Wave_Dissipation','reduction in stress due to dissipation of waves', &
                             'sea_surface_downward_stress_due_to_dissipation_of_sea_surface_waves','m^2/s','down')
    call add_item_to_cf_list('Wave_Rad_Str','wave xy radiation stress vectors','sea_surface_wave_xy_radiation_stress','m^2/s^2','up')
    call add_item_to_cf_list('Wave_Rad_Str_Mag','wave xy radiation stress magnitude','sea_surface_wave_xy_radiation_stress','m^2/s^2','up')
    call add_item_to_cf_list('Pressure','sea water pressure','sea_water_pressure_due_to_sea_water','m^2/s^2','up')
    call add_item_to_cf_list('Wind_Velocity','wind_velocity_vector','wind_u/v_velocity','m/s','u/v_direction')
    call add_item_to_cf_list('Wind_Magnitude','magnitude of wind u/v velocity','velocity_magnitude','m/s','up')
    call add_item_to_cf_list('Wind_Stress','wind xy stress vectors','wind_xy_stress','N/m^2','x/y_direction')
    call add_item_to_cf_list('Wind_Stress_Magnitude','magnitude of wind xy stress','wind_xy_stress_magnitude','N/m^2','up')
    call add_item_to_cf_list('Atm_Pressure','atmospheric pressure at sea level','air_pressure_at_mean_sea_level','Pa','up')
    
    return
    end subroutine init_cf_var_list

end module out_lib    
