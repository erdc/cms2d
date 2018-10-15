!======================================================================
module in_lib
! Input Library
!
! Contains the following routines:
!   SMS file formats:
!     read_xys       - Reads an SMS XY Series (*.xys) file
!     read_tsd       - Reads an SMS Time-Series Data (*.tsd) file
!     read_dat       - Reads an SMS Dataset (*.dat) file
!   XMDF file formats:
!     read_dataseth5 - Reads a dataset from an XMDF file
!     readscalh5     - Reads a scalar from an XMDF file 
!     readscalsteph5 - Reads a time step of a scalar from an XMDF file 
!     readscallasth5 - Reads the last time step of a scalar from an XMDF file 
!     readvech5      - Reads a vector from an XMDF file 
!     readvecsteph5  - Reads a time step of a vector from an XMDF file 
!     readveclasth5  - Reads the last time step of a vector from an XMDF file
!   ADCIRC solution files:
!     readstep63 - Reads an ADCIRC fort.63 solution step
!     readstep64 - Reads an ADCIRC fort.64 solution step
!
! written by Alex Sanchez, USACE-CHL
!======================================================================
#include "CMS_cpp.h"
    implicit none

contains
    
!*****************************************************************
    subroutine read_xys(afile,ns,xs,ys)
! Reads an SMS XYS File
!
!  Note: fortran has problems reading ASCII files in Unicode or 
!  other encodings. It is recommended to use ANSI encoding for all
!  ASCII files.
!
! written by Alex Sanchez, USACE-CHL
!*****************************************************************
    use diag_lib
    use prec_def
    implicit none    
    !Input/Output
    integer,intent(out) :: ns
    character(len=*),intent(in) :: afile
    real(ikind),intent(inout),pointer :: xs(:),ys(:)
    !Internal
    integer :: i,ierr
    character :: adum*3,msg*200
    logical :: foundfile
    
    inquire(file=afile,exist=foundfile)
    if(.not.foundfile)then
      write(msg,*,iostat=ierr) 'Could not find file: ',trim(afile)
      call diag_print_error(msg)
    endif     
    open(454,file=afile)
    read(454,*,iostat=ierr) adum,i,ns
    if(ierr/=0)then
      call diag_print_error('Problem reading header of XYS file: ',afile)  
    endif
    !if(allocated(xs)) deallocate(xs)
    !if(allocated(ys)) deallocate(ys)
    allocate(xs(ns),ys(ns))
    read(454,*,iostat=ierr) (xs(i),ys(i),i=1,ns)
    if(ierr/=0)then
      call diag_print_error('Problem reading data in XYS file: ',afile)  
    endif
    close(454)
    !write(*,*) afile,' read succesfully'
    
    return
    endsubroutine read_xys
    
!**************************************************************************
    subroutine read_tsd(afile,aname,atype,ndat,nt,tjuldaybeg,t,dat)
! Reads an SMS tsd file  
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
!
!  Note: fortran has problems reading ASCII files in Unicode or 
!  other encodings. It is recommended to use ANSI encoding for all
!  ASCII files.
!
! written by Alex Sanchez, USACE-CHL
!**************************************************************************
    use time_lib, only: calendar2julian
    use diag_lib
    use prec_def
    implicit none
    !Input/Output
    integer,intent(out) :: nt,ndat
    character(len=*),intent(in) :: afile
    character(len=*),intent(inout) :: aname,atype
    real(ikind), intent(out) :: tjuldaybeg
    real(ikind), intent(inout), pointer :: t(:),dat(:,:)
    !Internal variables
    integer :: i,j,iyear,imonth,iday,ihour,imin,isec,ncol,ierr
    character(len=19) :: astarttime
    character(len=200) :: msg
    logical :: found
        
    inquire(file=afile,exist=found)
    if(.not.found)then
      write(msg,*,iostat=ierr) 'Could not find file: ',trim(afile)
      call diag_print_error(msg)
    endif     
    open(454,file=afile)
    read(454,*,iostat=ierr) aname
    if(ierr/=0)then
      call diag_print_error('Problem reading file: ',trim(afile))  
    endif
    if(aname(1:11)/='TIME_SERIES')then
      backspace(454) !Check for missing first line
    endif
    read(454,*,iostat=ierr) aname,atype,ncol,nt,astarttime
    if(ierr/=0)then
      call diag_print_error('Problem reading header in TSD file: ',trim(afile))  
    endif
    if(nt<=0)then !Determine length of file        
      do i=1,1000000
        read(454,*,iostat=ierr)
        if(ierr/=0) exit
        nt = nt + 1
      enddo
      close(454)
      open(454,file=afile)
      read(454,*,iostat=ierr) msg !skip
      read(454,*,iostat=ierr) msg
    endif    
    if(ncol<2)then
      call diag_print_error('Invalid Number of Colums in TSD File: ',trim(afile),&
       '  Must be equal or greater than 2')
    endif
    ndat = ncol - 1 !Subtract time column to get number of data columns    
432 format(I2,'/',I2,'/',I4,1x,I2,':',I2,':',I2)  !01/01/2008 00:00:00
433 format(I2,'/',I2,'/',I4,2x,I2,':',I2,':',I2)  !01/01/2008  00:00:00  !Note extra space between date and time
    read(astarttime,432,iostat=ierr) imonth,iday,iyear,ihour,imin,isec
    if (ierr .ne. 0) read(astarttime,433,iostat=ierr) imonth,iday,iyear,ihour,imin,isec
    call calendar2julian(iyear,imonth,iday,ihour,imin,isec,tjuldaybeg)
    allocate(t(nt),dat(nt,ndat))
    do i=1,nt
      read(454,*,iostat=ierr) t(i),(dat(i,j),j=1,ndat)
      if(ierr/=0)then
        write(msg,*,iostat=ierr) '  at line',i+2
        call diag_print_error('Problem reading data in TSD file: ',trim(afile),msg)  
      endif
    enddo    
    close(454)
    !write(*,*) trim(afile),' read succesfully'
    
    return
    endsubroutine read_tsd
    
!Note: Other useful routines usages would be something like this
!  call read_dat(afile,nscal=ns,scaldat=sdat,name='Water_Elevation',step='Last')  
!  call read_dat(afile,nscal=ns,scaldat=sdat,name='Water_Elevation',time=6.0)      
!  call read_dat(afile,nvec=nv,vecdat=vdat,name='Current_Velocity',step='Last')  
!  call read_dat(afile,nvec=nv,vecdat=vdat,name='Current_Velocity',time=6.0)   

!  call read_scal(afile,nscal,scaldat,name='Water_Elevation',step='Last')  
!  call read_scal(afile,nscal,scaldat,name='Water_Elevation',time=6.0)      
!  call read_vec(afile,nvec,vecdat,name='Current_Velocity',step='Last')  
!  call read_vec(afile,nvec,vecdat,name='Current_Velocity',time=6.0)    
    
!*************************************************************    
    subroutine read_dat(afile,nscal,scaldat,nvec,vecdat)
! Reads an SMS Dataset (*.dat) file
! Reads all of the scalar and vector datasets and stores them
! in the variables scaldat and vecdat of length nscal and nvec,
! respectively. The datasets can then be parsed by name. 
!
! written by Alex Sanchez, USACE-CHL    
!*************************************************************    
    use prec_def
    use in_def
    implicit none
    !Input/Output
    character(len=*),  intent(in) :: afile
    integer, intent(out) :: nscal,nvec
    type(scaldattype), intent(inout), pointer :: scaldat(:)
    type(vecdattype),  intent(inout), pointer :: vecdat(:)    
    !Internal variables
    integer :: iunit,ierr
    character(len=37) :: cardname,objtype
    logical :: foundcard,founddataset
    
    foundcard = .true.
    founddataset = .false.
    nscal = 0
    nvec = 0
    iunit = 778
    open(iunit,file=afile)
    ierr = 0
    do while(ierr==0)
	  read(iunit,*,iostat=ierr) cardname
      if(ierr==-1) exit !End of File
      selectcase(cardname)           
      case('DATASET')
	    founddataset = .true.
      
      case('BEGSCL')  
        call read_scal(iunit,nscal,scaldat)
      
      case('BEGVEC')  
        call read_vec(iunit,nvec,vecdat)
        
      case('OBJTYPE')
	    backspace(iunit)
        read(iunit,*,iostat=ierr) cardname,objtype
        
      endselect    
    enddo
    close(iunit)
    
    return
    endsubroutine read_dat 

!*************************************************************    
    subroutine read_datnew(afile,ndat,dat)
! Reads an SMS Dataset (*.dat) file
! Reads all of the scalar and vector datasets and stores them
! in the variables scaldat and vecdat of length nscal and nvec,
! respectively. The datasets can then be parsed by name. 
!
! written by Alex Sanchez, USACE-CHL    
!*************************************************************    
    use prec_def
    use in_def
    implicit none
    !Input/Output
    character(len=*),  intent(in) :: afile
    integer, intent(out) :: ndat
    type(dat_type), intent(inout), pointer :: dat(:)
    !Internal variables
    integer :: iunit,ierr
    character(len=37) :: cardname,objtype
    logical :: foundcard,founddataset
    
    foundcard = .true.
    founddataset = .false.
    ndat = 0
    iunit = 778
    open(iunit,file=afile)
    ierr = 0
    do while(ierr==0)
	  read(iunit,*,iostat=ierr) cardname
      if(ierr==-1) exit !End of File
      selectcase(cardname)           
      case('DATASET')
	    founddataset = .true.
      
      case('BEGSCL')  
        call read_scalvec(iunit,1,ndat,dat)
      
      case('BEGVEC')   
        call read_scalvec(iunit,2,ndat,dat)
        
      case('OBJTYPE')
	    backspace(iunit)
        read(iunit,*,iostat=ierr) cardname,objtype
        
      endselect    
    enddo
    close(iunit)
    
    return
    endsubroutine read_datnew 

!**************************************************************************
    subroutine read_scalvec(iunit,ndim,ndat,dat)
! Reads a scalar dataset from an SMS ASCII data file (*.dat)    
! written by Alex Sanchez, USACE-CHL    
!**************************************************************************
    use prec_def
    use diag_lib
    use in_def
    implicit none
    !Input/Output
    integer,intent(in) :: iunit,ndim
    integer,intent(inout) :: ndat
    type(dat_type), intent(inout), pointer :: dat(:)    
    !Internal variables
    integer :: i,j,ierr
    type(dat_type), allocatable :: dat_temp(:)
    character(len=37) :: cardname
    
    !Allocate temporary variable  
    if(ndat>0)then
      allocate(dat_temp(ndat))     
      do i=1,ndat
        allocate(dat_temp(i)%time(dat(i)%nt))
        allocate(dat_temp(i)%val(dat(i)%nd,dat(i)%ndim,dat(i)%nt))
        if(dat(i)%istatus==1)then
          allocate(dat_temp(1)%stat(dat(i)%nd,dat(i)%nt))
        endif
      enddo
    
      !Copy information to temporary variable
      dat_temp = dat
    
      !Deallocate variable
      deallocate(dat)
    endif
    
    !Resize
    ndat = ndat + 1
    
    !Re-allocate variable
    allocate(dat(ndat)) 
    do i=1,ndat
      allocate(dat(i)%time(dat(i)%nt))
      allocate(dat(i)%val(dat(i)%nd,dat(i)%ndim,dat(i)%nt))
      if(dat(i)%istatus==1)then
        allocate(dat(i)%stat(dat(i)%nd,dat(i)%nt))
      endif
    enddo
    
    if(ndat>1)then
      !Copy information back from temporary variable
      dat(1:ndat-1) = dat_temp(1:ndat-1)
    
      !Dellocate temporary variable
      deallocate(dat_temp)
    endif
    
    !Set defaults
    dat(ndat)%name = ''
    dat(ndat)%id = 0 !ID #
    dat(ndat)%nd = 0 !# of data values
    dat(ndat)%nc = 0 !# of cells or elements
    dat(ndat)%nt = 0 !# of time steps
    dat(ndat)%usereftime = .false. !Reference time specified by user
    dat(ndat)%reftime = -1.0       !Reference time
    dat(ndat)%time_units = 'seconds'
    dat(ndat)%itype = 0   !0-applied to nodes, 1-applied to elements/cells    
    dat(ndat)%ndim  = ndim   !1-Scalar, 2-2D Vector
    
    !Read information from file
    ierr = 0
    do while(ierr==0)
	  read(iunit,*,iostat=ierr) cardname      
      selectcase(cardname)        
      case('OBJTYPE')
	    backspace(iunit)
        read(iunit,*,iostat=ierr) cardname,dat(ndat)%id
        
      case('ND')
	    backspace(iunit)
        read(iunit,*,iostat=ierr) cardname,dat(ndat)%nd
        
      case('NC')
	    backspace(iunit)
        read(iunit,*,iostat=ierr) cardname,dat(ndat)%nc
      
      case('NAME')
	    backspace(iunit)
        read(iunit,*,iostat=ierr) cardname,dat(ndat)%name  
        
      case('RT_JULIAN')
        backspace(iunit)
        read(iunit,*,iostat=ierr) cardname,dat(ndat)%reftime
        dat(ndat)%usereftime = .true.
        
      case('TIMEUNITS')
	    backspace(iunit)
        read(iunit,*,iostat=ierr) cardname,dat(ndat)%time_units    
      
      case('TS')  
        !Increase time by one
        dat(ndat)%nt = dat(ndat)%nt + 1 
        
        !Read new time step
        allocate(dat_temp(1))                
        allocate(dat_temp(1)%time(dat(ndat)%nt))
        backspace(iunit)
        read(iunit,*,iostat=ierr) cardname,dat(ndat)%istatus,dat_temp(1)%time(dat(ndat)%nt)
        if(dat(ndat)%istatus==1)then !Read activity or status
          allocate(dat_temp(1)%stat(dat(ndat)%nd,dat(ndat)%nt)) 
          read(iunit,*,iostat=ierr) (dat_temp(1)%stat(i,dat(ndat)%nt),i=1,dat(ndat)%nd)
        endif
        allocate(dat_temp(1)%val(dat(ndat)%nd,dat(ndat)%ndim,dat(ndat)%nt))
        read(iunit,*,iostat=ierr) ((dat_temp(1)%val(i,j,dat(ndat)%nt), &
                                   j=1,dat(ndat)%ndim),i=1,dat(ndat)%nd)
        
        !Copy previous time steps to temporary variable
        if(dat(ndat)%nt>1)then
          dat_temp(1)%time(1:dat(ndat)%nt-1) = dat(ndat)%time(1:dat(ndat)%nt-1)
          dat_temp(1)%val(:,:,1:dat(ndat)%nt-1) = dat(ndat)%val(:,:,1:dat(ndat)%nt-1)
          if(dat(ndat)%istatus==1)then
            dat_temp(1)%stat(:,1:dat(ndat)%nt-1) = dat(ndat)%stat(:,1:dat(ndat)%nt-1)
          endif
        endif
      
        !Deallocate variables
        deallocate(dat(ndat)%time)
        deallocate(dat(ndat)%val)
        if(dat(ndat)%istatus==1 .and. allocated(dat(ndat)%stat))then
          deallocate(dat(ndat)%stat)
        endif  
        
        dat(ndat)%ndim = dat_temp(1)%ndim
        
        !Re-allocate variables
        allocate(dat(ndat)%time(dat(ndat)%nt))
        allocate(dat(ndat)%val(dat(ndat)%nd,dat(ndat)%ndim,dat(ndat)%nt))
        if(dat(ndat)%istatus==1)then
          allocate(dat(ndat)%stat(dat(ndat)%nd,dat(ndat)%nt))
        endif
        
        !Copy back information from temporary variable
        dat(ndat)%time(:) = dat_temp(1)%time(:)
        dat(ndat)%val(:,:,:) = dat_temp(1)%val(:,:,:)
        if(dat(ndat)%istatus==1)then
          dat(ndat)%stat(:,:) = dat_temp(1)%stat(:,:)
        endif
                        
        !Deallocate temporary variable
        deallocate(dat_temp)        
        
      case('ENDDS')
        exit
        
      endselect    
    enddo
    
    if(ierr/=0)then
      call diag_print_error('Problem reading dataset: ',dat(ndat)%name)    
    endif        

    return
    endsubroutine read_scalvec     
    
!**************************************************************************
    subroutine read_scal(iunit,nscal,scaldat)
! Reads a scalar dataset from an SMS ASCII data file (*.dat)    
! written by Alex Sanchez, USACE-CHL    
!**************************************************************************
    use prec_def
    use diag_lib
    use in_def
    implicit none
    !Input/Output
    integer,intent(in) :: iunit
    integer,intent(inout) :: nscal
    type(scaldattype), intent(inout), pointer :: scaldat(:)    
    !Internal variables
    integer :: i,ierr
    type(scaldattype), allocatable :: scaldattemp(:)
    character(len=37) :: cardname
    
    !Allocate temporary variable  
    if(nscal>0)then
      allocate(scaldattemp(nscal))     
      do i=1,nscal
        allocate(scaldattemp(i)%time(scaldat(i)%nt))
        allocate(scaldattemp(i)%val(scaldat(i)%nd,scaldat(i)%nt))
        if(scaldat(i)%istatus==1)then
          allocate(scaldattemp(1)%stat(scaldat(i)%nd,scaldat(i)%nt))
        endif
      enddo
    
      !Copy information to temporary variable
      scaldattemp = scaldat
    
      !Deallocate variable
      deallocate(scaldat)
    endif
    
    !Resize
    nscal = nscal + 1
    
    !Re-allocate variable
    allocate(scaldat(nscal)) 
    do i=1,nscal
      allocate(scaldat(i)%time(scaldat(i)%nt),scaldat(i)%val(scaldat(i)%nd,scaldat(i)%nt))
      if(scaldat(i)%istatus==1)then
        allocate(scaldat(i)%stat(scaldat(i)%nd,scaldat(i)%nt))
      endif
    enddo
    
    if(nscal>1)then
      !Copy information back from temporary variable
      scaldat(1:nscal-1) = scaldattemp(1:nscal-1)
    
      !Dellocate temporary variable
      deallocate(scaldattemp)
    endif
    
    !Set defaults
    scaldat(nscal)%name = ''
    scaldat(nscal)%id = 0
    scaldat(nscal)%nd = 0
    scaldat(nscal)%nc = 0
    scaldat(nscal)%nt = 0
    scaldat(nscal)%reftime = -1
    scaldat(nscal)%time_units = 'seconds'
    
    !Read information from file
    ierr = 0
    do while(ierr==0)
	  read(iunit,*,iostat=ierr) cardname      
      selectcase(cardname)        
      case('OBJTYPE')
	    backspace(iunit)
        read(iunit,*,iostat=ierr) cardname,scaldat(nscal)%id
        
      case('ND')
	    backspace(iunit)
        read(iunit,*,iostat=ierr) cardname,scaldat(nscal)%nd
        
      case('NC')
	    backspace(iunit)
        read(iunit,*,iostat=ierr) cardname,scaldat(nscal)%nc
      
      case('NAME')
	    backspace(iunit)
        read(iunit,*,iostat=ierr) cardname,scaldat(nscal)%name  
        
      case('RT_JULIAN')
        backspace(iunit)
        read(iunit,*,iostat=ierr) cardname,scaldat(nscal)%reftime
        scaldat(nscal)%usereftime = .true.
        
      case('TIMEUNITS')
	    backspace(iunit)
        read(iunit,*,iostat=ierr) cardname,scaldat(nscal)%time_units    
      
      case('TS')  
        !Increase time by one
        scaldat(nscal)%nt = scaldat(nscal)%nt + 1 
        
        !Read new time step
        allocate(scaldattemp(1))                
        allocate(scaldattemp(1)%time(scaldat(nscal)%nt))
        backspace(iunit)
        read(iunit,*,iostat=ierr) cardname,scaldat(nscal)%istatus,scaldattemp(1)%time(scaldat(nscal)%nt)
        if(scaldat(nscal)%istatus==1)then !Read activity or status
          allocate(scaldattemp(1)%stat(scaldat(nscal)%nd,scaldat(nscal)%nt))  
          read(iunit,*,iostat=ierr) (scaldattemp(1)%stat(i,scaldat(nscal)%nt),i=1,scaldat(nscal)%nd)
        endif
        allocate(scaldattemp(1)%val(scaldat(nscal)%nd,scaldat(nscal)%nt))
        read(iunit,*,iostat=ierr) (scaldattemp(1)%val(i,scaldat(nscal)%nt),i=1,scaldat(nscal)%nd)
        
        !Copy previous time steps to temporary variable
        if(scaldat(nscal)%nt>1)then
          scaldattemp(1)%time(1:scaldat(nscal)%nt-1) = scaldat(nscal)%time(1:scaldat(nscal)%nt-1)
          scaldattemp(1)%val(:,1:scaldat(nscal)%nt-1) = scaldat(nscal)%val(:,1:scaldat(nscal)%nt-1)
          if(scaldat(nscal)%istatus==1)then
            scaldattemp(1)%stat(:,1:scaldat(nscal)%nt-1) = scaldat(nscal)%stat(:,1:scaldat(nscal)%nt-1)
          endif
        endif
      
        !Deallocate variables
        deallocate(scaldat(nscal)%time)
        deallocate(scaldat(nscal)%val)
        if(scaldat(nscal)%istatus==1 .and. allocated(scaldat(nscal)%stat))then
          deallocate(scaldat(nscal)%stat)
        endif        
        
        !Re-allocate variables
        allocate(scaldat(nscal)%time(scaldat(nscal)%nt))
        allocate(scaldat(nscal)%val(scaldat(nscal)%nd,scaldat(nscal)%nt))
        if(scaldat(nscal)%istatus==1)then
          allocate(scaldat(nscal)%stat(scaldat(nscal)%nd,scaldat(nscal)%nt))
        endif
        
        !Copy back information from temporary variable
        scaldat(nscal)%time(:) = scaldattemp(1)%time(:)
        scaldat(nscal)%val(:,:) = scaldattemp(1)%val(:,:)
        if(scaldat(nscal)%istatus==1)then
          scaldat(nscal)%stat(:,:) = scaldattemp(1)%stat(:,:)
        endif
                        
        !Deallocate temporary variable
        deallocate(scaldattemp)        
        
      case('ENDDS')
        exit
        
      endselect    
    enddo
    
    if(ierr/=0)then
      call diag_print_error('Problem reading scalar dataset: ',scaldat(nscal)%name)    
    endif        

    return
    endsubroutine read_scal 
    
!*************************************************************
    subroutine read_vec(iunit,nvec,vecdat) !name,time,step
! written by Alex Sanchez, USACE-CHL    
!*************************************************************
    use prec_def
    use diag_lib
    use in_def
    implicit none
    !Input/Output
    integer,intent(in) :: iunit
    integer,intent(inout) :: nvec
    type(vecdattype), intent(inout), pointer :: vecdat(:)    
    !Internal variables
    integer :: i,j,ierr
    type(vecdattype), allocatable :: vecdattemp(:)
    character(len=37) :: cardname
    
    if(nvec>0)then
      !Allocate temporary variable   
      allocate(vecdattemp(nvec))     
      do i=1,nvec
        allocate(vecdattemp(i)%time(vecdat(i)%nt))
        allocate(vecdattemp(i)%val(vecdat(i)%nd,2,vecdat(i)%nt))
        if(vecdat(i)%istatus==1)then
          allocate(vecdattemp(1)%stat(vecdat(i)%nd,vecdat(i)%nt))
        endif
      enddo
    
      !Copy information to temporary variable
      vecdattemp = vecdat
    
      !Deallocate variable
      deallocate(vecdat)
    endif
    
    !Resize
    nvec = nvec + 1
    
    !Re-allocate variable
    allocate(vecdat(nvec)) 
    do i=1,nvec
      allocate(vecdat(i)%time(vecdat(i)%nt))
      allocate(vecdat(i)%val(vecdat(i)%nd,2,vecdat(i)%nt))
      if(vecdat(i)%istatus==1)then
        allocate(vecdat(i)%stat(vecdat(i)%nd,vecdat(i)%nt))
      endif
    enddo
    
    if(nvec>1)then
      !Copy information back from temporary variable
      vecdat(1:nvec-1) = vecdattemp(1:nvec-1)
    
      !Dellocate temporary variable
      deallocate(vecdattemp)
    endif
    
    
    !Set defaults
    vecdat(nvec)%name = ''
    vecdat(nvec)%id = 0
    vecdat(nvec)%nd = 0
    vecdat(nvec)%nc = 0
    vecdat(nvec)%nt = 0
    vecdat(nvec)%reftime = -1
    vecdat(nvec)%time_units = 'seconds'
    
    !Read information from file
    ierr = 0
    do while(ierr==0)
	  read(iunit,*,iostat=ierr) cardname      
      selectcase(cardname)        
      case('OBJTYPE')
	    backspace(iunit)
        read(iunit,*,iostat=ierr) cardname,vecdat(nvec)%id
        
      case('ND')
	    backspace(iunit)
        read(iunit,*,iostat=ierr) cardname,vecdat(nvec)%nd
        
      case('NC')
	    backspace(iunit)
        read(iunit,*,iostat=ierr) cardname,vecdat(nvec)%nc
      
      case('NAME')
	    backspace(iunit)
        read(iunit,*,iostat=ierr) cardname,vecdat(nvec)%name  
        
      case('RT_JULIAN')
        backspace(iunit)
        read(iunit,*,iostat=ierr) cardname,vecdat(nvec)%reftime
        vecdat(nvec)%usereftime = .true.
        
      case('TIMEUNITS')
	    backspace(iunit)
        read(iunit,*,iostat=ierr) cardname,vecdat(nvec)%time_units    
      
      case('TS')
	    !Increase time by one
        vecdat(nvec)%nt = vecdat(nvec)%nt + 1 
        
        !Read new time step
        allocate(vecdattemp(1))                
        allocate(vecdattemp(1)%time(vecdat(nvec)%nt))
        backspace(iunit)
        read(iunit,*,iostat=ierr) cardname,vecdat(nvec)%istatus,vecdattemp(1)%time(vecdat(nvec)%nt)
        if(vecdat(nvec)%istatus==1)then !Read activity or status
          allocate(vecdattemp(1)%stat(vecdat(nvec)%nd,vecdat(nvec)%nt))  
          read(iunit,*,iostat=ierr) (vecdattemp(1)%stat(i,vecdat(nvec)%nt),i=1,vecdat(nvec)%nd)
        endif
        allocate(vecdattemp(1)%val(vecdat(nvec)%nd,2,vecdat(nvec)%nt))
        read(iunit,*,iostat=ierr) ((vecdattemp(1)%val(i,j,vecdat(nvec)%nt),j=1,2),i=1,vecdat(nvec)%nd)
        
        !Copy previous time steps to temporary variable
        if(vecdat(nvec)%nt>1)then
          vecdattemp(1)%time(1:vecdat(nvec)%nt-1) = vecdat(nvec)%time(1:vecdat(nvec)%nt-1)
          vecdattemp(1)%val(:,:,1:vecdat(nvec)%nt-1) = vecdat(nvec)%val(:,:,1:vecdat(nvec)%nt-1)
          if(vecdat(nvec)%istatus==1)then
            vecdattemp(1)%stat(:,1:vecdat(nvec)%nt-1) = vecdat(nvec)%stat(:,1:vecdat(nvec)%nt-1)
          endif
        endif
      
        !Deallocate variables
        deallocate(vecdat(nvec)%time)
        deallocate(vecdat(nvec)%val)
        if(vecdat(nvec)%istatus==1 .and. allocated(vecdat(nvec)%stat))then
          deallocate(vecdat(nvec)%stat)
        endif        
        
        !Re-allocate variables
        allocate(vecdat(nvec)%time(vecdat(nvec)%nt))
        allocate(vecdat(nvec)%val(vecdat(nvec)%nd,2,vecdat(nvec)%nt))
        if(vecdat(nvec)%istatus==1)then
          allocate(vecdat(nvec)%stat(vecdat(nvec)%nd,vecdat(nvec)%nt))
        endif
        
        !Copy back information from temporary variable
        vecdat(nvec)%time(:) = vecdattemp(1)%time(:)
        vecdat(nvec)%val(:,:,:) = vecdattemp(1)%val(:,:,:)
        if(vecdat(nvec)%istatus==1)then
          vecdat(nvec)%stat(:,:) = vecdattemp(1)%stat(:,:)
        endif
                        
        !Deallocate temporary variable
        deallocate(vecdattemp)  
        
      case('ENDDS')
        exit
        
      endselect    
    enddo
    
    if(ierr/=0)then
      call diag_print_error('Problem reading vector dataset: ',vecdat(nvec)%name)
    endif  

    return
    endsubroutine read_vec    
    
#ifdef XMDF_IO
!*************************************************************      
    subroutine read_dataseth5(afile,apath,aname,nn,vec)
! Reads a dataset from an XMDF file
! written by Alex Sanchez, USACE-CHL
!*************************************************************  
    use xmdf
    use prec_def
    use diag_lib, only: diag_print_error
    implicit none
    
    !Input/Output
    character(len=*),intent(in) :: afile,apath,aname
    integer,intent(out) :: nn
    real(ikind),intent(inout),pointer :: vec(:)
    
    !Internal variables
    integer :: pid,gid,ierr    
    real(4),allocatable :: ftemp(:) !Should be single    
       
    call XF_OPEN_FILE(afile,READONLY,pid,ierr)
    if (ierr < 0) call diag_print_error ('Unable to open model parameters file.')
    
    call XF_OPEN_GROUP(pid,apath,gid,ierr)
    if (ierr < 0) call diag_print_error ('Unable to open dataset group.','Check model parameters file for correct input.')

    call XF_GET_PROPERTY_NUMBER(gid,trim(aname),nn,ierr)
    if (ierr < 0) call diag_print_error ('Unable to obtain number of times.')

    allocate(vec(nn),ftemp(nn))
    call XF_READ_PROPERTY_FLOAT(gid,trim(aname),nn,ftemp(1),ierr)
    if (ierr < 0) call diag_print_error ('Unable to get list of times.')

    vec = ftemp !Useful for converting precision    
    call XF_CLOSE_GROUP(gid,ierr)
    call XF_CLOSE_FILE(pid,ierr)    
    deallocate(ftemp)  
    
    return
    endsubroutine read_dataseth5
    
!************************************************************************
    subroutine readscalh5(afile,apath,var,ierr)
! Reads a scalar dataset located in apath from afile XMDF file 
! written by Alex Sanchez, USACE-CHL
!************************************************************************
    use size_def, only: ncellsD,ncellsfull,ncellpoly
    use geo_def, only: idmap
    use interp_lib, only: interp_scal_node2cell
    use diag_lib
    use xmdf
    use prec_def
    implicit none
    !Input/Output
    character(len=*), intent(in) :: afile,apath
    real(ikind), intent(out) :: var(ncellsD)
    integer, intent(out) :: ierr
    !Internal Variables
    integer :: fid,gid
    real(4) :: vtemp(ncellsfull)
    character(len=200) :: msg2,msg3

    call XF_OPEN_FILE(trim(afile),READONLY,fid,ierr)        
    if(ierr<0) call diag_print_error('Could not open file: ',trim(afile))
    
    call XF_OPEN_GROUP(fid,trim(apath),gid,ierr)
    if(ierr<0)then
      call XF_CLOSE_FILE(fid,ierr)
      write(msg2,*,iostat=ierr) '  File: ',trim(afile)
      write(msg3,*,iostat=ierr) '  Path: ',trim(apath)
      call diag_print_error('Could not open dataset from',msg2,msg3)
      return
    endif
       
    call XF_READ_SCALAR_VALUES_TIMESTEP(gid,1,ncellsfull,vtemp,ierr)
    if(ierr<0)then
      !call XF_CLOSE_GROUP(gid,ierr)
      !call XF_CLOSE_FILE(fid,ierr)  
      !write(msg2,*,iostat=ierr) '  File: ',trim(afile)
      !write(msg3,*,iostat=ierr) '  Path: ',trim(apath)
      !call diag_print_error('Could not read scalar dataset from',msg2,msg3)
      return
    endif
    
    call XF_CLOSE_GROUP(gid,ierr)
    call XF_CLOSE_FILE(fid,ierr)
                
    if(ncellpoly>0)then
      call interp_scal_node2cell(vtemp,var) !Interpolate node to cell centers
    else
      call map_scal_full2active(vtemp,var) !Convert from full to active grid 
    endif
      
    return
    endsubroutine readscalh5
      
!************************************************************************
    subroutine readscalsteph5(afile,apath,itsind,thrs,var,ierr)
! Reads the a time step of a scalar located in apath from afile XMDF file 
! written by Alex Sanchez, USACE-CHL
!************************************************************************
    use size_def, only: ncellsD,ncellsfull,ncellpoly
    use geo_def, only: idmap
    use interp_lib, only: interp_scal_node2cell
    use diag_def
    use diag_lib
    use xmdf
    use prec_def
    implicit none
    !Input/Output
    character(len=*), intent(in) :: afile,apath
    integer, intent(in) :: itsind !Time step index
    real(ikind), intent(out):: var(ncellsD),thrs
    integer, intent(out) :: ierr    
    !Internal Variables
    integer :: fid,gid,ntimes
    real(8),allocatable :: timesd(:) !Output times
    real(4) :: vtemp(ncellsfull) !Must be single
    
    call XF_OPEN_FILE(trim(afile),READONLY,fid,ierr)        
    if(ierr<0)then
      call diag_print_error('Could not open file: ',afile)
    endif
          
    call XF_OPEN_GROUP(fid,trim(apath),gid,ierr)
    if(ierr<0)then
      call XF_CLOSE_FILE(fid,ierr)
      !write(msg2,*,iostat=ierr) '  File: ',trim(afile)
      !write(msg3,*,iostat=ierr) '  Path: ',trim(apath)
      !call diag_print_warning('Could not open dataset from',msg2,msg3)
      ierr = -2 !Could not open group
      return
    endif     
       
    call XF_READ_SCALAR_VALUES_TIMESTEP(gid,itsind,ncellsfull,vtemp,ierr)
    if(ierr<0)then
      call XF_CLOSE_GROUP(gid,ierr)
      call XF_CLOSE_FILE(fid,ierr)  
      !write(msg2,*,iostat=ierr) '  File: ',trim(afile)
      !write(msg3,*,iostat=ierr) '  Path: ',trim(apath)
      !write(msg4,*,iostat=ierr) '  Time Step: ',itsind
      !call diag_print_warning('Could not read scalar time step',msg2,msg3,msg4)
      ierr = 3 !Could not read timestep value
      return
    endif 

    call XF_GET_DATASET_NUM_TIMES(gid,ntimes,ierr)
    allocate(timesd(ntimes))
    call XF_GET_DATASET_TIMES(gid,ntimes,timesd,ierr)
    if(ierr<0)then
      thrs = -999.0
      write(msg2,*,iostat=ierr) '  File: ',trim(afile)
      write(msg3,*,iostat=ierr) '  Path: ',trim(apath)
      call diag_print_warning('Could not read time stamp from',msg2,msg3)
      ierr = 4 !Could not read time stamp
    else
      thrs = timesd(itsind)
    endif
    
    call XF_CLOSE_GROUP(gid,ierr)
    call XF_CLOSE_FILE(fid,ierr)
           
    if(ncellpoly>0)then
      call interp_scal_node2cell(vtemp,var) !Interpolate node to cell centers
    else
      call map_scal_full2active(vtemp,var) !Convert from full to active grid 
    endif
      
    return
    endsubroutine readscalsteph5

!************************************************************************
    subroutine readscallasth5(afile,apath,ntimes,var,reftimed,thrs,ierr)
! Reads a scalar located in apath from afile XMDF file 
! written by Alex Sanchez, USACE-CHL
!************************************************************************
!#include "CMS_cpp.h"
    use size_def, only: ncellsD,ncellsfull,ncellpoly
    use geo_def, only: idmap
    use interp_lib, only: interp_scal_node2cell
    use diag_def
    use diag_lib
    use xmdf
    use prec_def
    implicit none
    !Input/Output
    character(len=*), intent(in) :: afile,apath
    real(ikind), intent(out) :: var(ncellsD),thrs    
    real(8),intent(out) :: reftimed          !Reference time 
    integer, intent(out) :: ierr,ntimes
    !Internal Variables
    integer :: fid,gid    
    real(8), allocatable :: timed(:) !Output time
    real(4) :: vtemp(ncellsfull) !Must be single
    
    call XF_OPEN_FILE(trim(afile),READONLY,fid,ierr)   
    if(ierr<0) call diag_print_error('Could not open file: ',afile)
      
    call XF_OPEN_GROUP(fid,trim(apath),gid,ierr)
    if(ierr<0)then
      call XF_CLOSE_FILE(fid,ierr)
!      write(msg2,*,iostat=ierr) '  File: ',trim(afile)
!      write(msg3,*,iostat=ierr) '  Path: ',trim(apath)
!      call diag_print_warning('Could not open dataset from ',msg2,msg3)
      ierr = -2
      return
    endif

    call XF_GET_DATASET_NUM_TIMES(gid,ntimes,ierr)
    if(ierr<0)then
      call XF_CLOSE_GROUP(gid,ierr)
      call XF_CLOSE_FILE(fid,ierr)
      write(msg2,*,iostat=ierr) '  File: ',trim(afile)
      write(msg3,*,iostat=ierr) '  Path: ',trim(apath)
      call diag_print_error('Could not read number of times from',msg2,msg3)
    endif
    
    call XF_READ_SCALAR_VALUES_TIMESTEP(gid,ntimes,ncellsfull,vtemp,ierr)
    if(ierr<0)then
      call XF_CLOSE_GROUP(gid,ierr)
      call XF_CLOSE_FILE(fid,ierr)  
      write(msg2,*,iostat=ierr) '  File: ',trim(afile)
      write(msg3,*,iostat=ierr) '  Path: ',trim(apath)
      write(msg4,*,iostat=ierr) '  Time step: ',ntimes
      call diag_print_error('Could not read last time step scalar values from',msg2,msg3,msg4)
    endif
    
    allocate(timed(ntimes))
    call XF_GET_DATASET_TIMES(gid,ntimes,timed,ierr)
    if(ierr<0)then
      thrs = -999.0
      write(msg2,*,iostat=ierr) '  File: ',trim(afile)
      write(msg3,*,iostat=ierr) '  Path: ',trim(apath)
      write(msg4,*,iostat=ierr) '  Time step: ',ntimes
      call diag_print_warning('Could not read times from',msg2,msg3,msg4)
      ierr = 4 !Could not read time stamp
    else
      thrs = timed(ntimes)
    endif
    deallocate(timed)
        
    call XF_GET_DATASET_REFTIME(gid,reftimed,ierr)
    if(ierr<0)then
      reftimed = -999.0d0
      write(msg2,*,iostat=ierr) '  File: ',trim(afile)
      write(msg3,*,iostat=ierr) '  Path: ',trim(apath)
      call diag_print_warning('Could not read reference time from',msg2,msg3)
      ierr = 5 !Could not read time stamp
    endif  

    call XF_CLOSE_GROUP(gid,ierr)
    call XF_CLOSE_FILE(fid,ierr)
            
    if(ncellpoly>0)then
      call interp_scal_node2cell(vtemp,var) !Interpolate node to cell centers
    else
      call map_scal_full2active(vtemp,var) !Convert from full to active grid 
    endif
      
    return
    endsubroutine readscallasth5
    
!************************************************************************
    subroutine readscaltimeh5(afile,apath,thrs,var,reftimed,ierr)
! Reads a scalar located in apath from afile XMDF file 
! written by Alex Sanchez, USACE-CHL
!************************************************************************
!#include "CMS_cpp.h"
    use size_def, only: ncellsD,ncellsfull,ncellpoly
    use geo_def, only: idmap
    use interp_lib, only: interp_scal_node2cell
    use diag_def
    use diag_lib
    use xmdf
    use prec_def
    implicit none
    !Input
    character(len=*), intent(in) :: afile,apath
    real(ikind), intent(in) :: thrs
    !Output
    real(ikind), intent(out) :: var(ncellsD)    
    real(8),intent(out) :: reftimed          !Reference time 
    integer, intent(out) :: ierr
    !Internal Variables
    integer :: i,fid,gid,nstep,ntimes
    real(8), allocatable :: timed(:) !Output time
    real(8) :: thrsd,terrd,terrdmin
    real(4) :: vtemp(ncellsfull) !Must be single
    
    call XF_OPEN_FILE(trim(afile),READONLY,fid,ierr)   
    if(ierr<0) call diag_print_error('Could not open file: ',afile)
      
    call XF_OPEN_GROUP(fid,trim(apath),gid,ierr)
    if(ierr<0)then
      call XF_CLOSE_FILE(fid,ierr)
!#ifdef DIAG_MODE
!      write(msg2,*,iostat=ierr) '  File: ',trim(afile)
!      write(msg3,*,iostat=ierr) '  Path: ',trim(apath)
!      call diag_print_warning('Could not open dataset from ',msg2,msg3)
!#endif
      ierr = -2
      return
    endif

    call XF_GET_DATASET_NUM_TIMES(gid,ntimes,ierr)
    if(ierr<0)then
      call XF_CLOSE_GROUP(gid,ierr)
      call XF_CLOSE_FILE(fid,ierr)
      write(msg2,*,iostat=ierr) '  File: ',trim(afile)
      write(msg3,*,iostat=ierr) '  Path: ',trim(apath)
      call diag_print_error('Could not read number of times from',msg2,msg3)
    endif
    
    allocate(timed(ntimes))
    call XF_GET_DATASET_TIMES(gid,ntimes,timed,ierr)
    if(ierr<0)then
      write(msg2,*,iostat=ierr) '  File: ',trim(afile)
      write(msg3,*,iostat=ierr) '  Path: ',trim(apath)
      write(msg4,*,iostat=ierr) '  Time step: ',ntimes
      call diag_print_warning('Could not read times from',msg2,msg3,msg4)
      ierr = 4 !Could not read time stamp
    else
      nstep = 1
      terrdmin = 1.0e9
      thrsd = dble(thrs)
      do i=1,ntimes
        terrd = abs(thrsd-timed(i))
        if(terrd<terrdmin)then
          nstep = i
          terrdmin = terrd
        endif
      enddo
      if(terrdmin>0.001)then
        write(msg2,*,iostat=ierr) '  File: ',trim(afile)
        write(msg3,*,iostat=ierr) '  Path: ',trim(apath)
        write(msg4,*,iostat=ierr) '  Time: ',thrs,' hrs'
        call diag_print_warning('Could not find time: ',msg2,msg3,msg4)
      endif
    endif
    deallocate(timed)
    
    call XF_READ_SCALAR_VALUES_TIMESTEP(gid,nstep,ncellsfull,vtemp,ierr)
    if(ierr<0)then
      call XF_CLOSE_GROUP(gid,ierr)
      call XF_CLOSE_FILE(fid,ierr)  
      write(msg2,*,iostat=ierr) '  File: ',trim(afile)
      write(msg3,*,iostat=ierr) '  Path: ',trim(apath)
      write(msg4,*,iostat=ierr) '  Time step: ',ntimes
      call diag_print_error('Could not read last time step scalar values from',msg2,msg3,msg4)
    endif
    
    call XF_GET_DATASET_REFTIME(gid,reftimed,ierr)
    if(ierr<0)then
      reftimed = -999.0d0
      write(msg2,*,iostat=ierr) '  File: ',trim(afile)
      write(msg3,*,iostat=ierr) '  Path: ',trim(apath)
      call diag_print_warning('Could not read reference time from',msg2,msg3)
      ierr = 5 !Could not read time stamp
    endif  

    call XF_CLOSE_GROUP(gid,ierr)
    call XF_CLOSE_FILE(fid,ierr)
            
    if(ncellpoly>0)then
      call interp_scal_node2cell(vtemp,var) !Interpolate node to cell centers
    else
      call map_scal_full2active(vtemp,var) !Convert from full to active grid 
    endif
      
    return
    endsubroutine readscaltimeh5  

!************************************************************************
    subroutine readvech5(afile,apath,vecx,vecy,ierr)
! Reads a scalar located in apath from afile XMDF file 
! written by Alex Sanchez, USACE-CHL
!************************************************************************
    use size_def
    use geo_def, only: idmap
    use interp_lib, only: interp_vec_node2cell
    use diag_lib
    use xmdf
    use prec_def
    implicit none
    !Input/Output
    character(len=200), intent(in) :: afile,apath
    real(ikind), intent(out) :: vecx(ncellsD),vecy(ncellsD)
    integer, intent(out) :: ierr
    !Internal
    integer :: fid,gid
    real(4) :: vtemp(ncellsfull*2) !Must be single
    character(len=200) :: msg2,msg3

    call XF_OPEN_FILE(trim(afile),READONLY,fid,ierr)        
    if(ierr<0) call diag_print_error('Could not open file: ',afile)
      
    call XF_OPEN_GROUP(fid,trim(apath),gid,ierr)
    if(ierr<0)then
      call XF_CLOSE_FILE(fid,ierr)  
      write(msg2,*,iostat=ierr) '  File: ',trim(afile)
      write(msg3,*,iostat=ierr) '  Path: ',trim(apath)
      call diag_print_error('Could not open dataset from',msg2,msg3)
      return
    endif
    
    call XF_READ_VECTOR_VALUES_TIMESTEP(gid,1,ncellsfull,2,vtemp,ierr)     
    if(ierr<0)then
      call XF_CLOSE_GROUP(gid,ierr)
      call XF_CLOSE_FILE(fid,ierr)
      write(msg2,*,iostat=ierr) '  File: ',trim(afile)
      write(msg3,*,iostat=ierr) '  Path: ',trim(apath)
      call diag_print_error('Could not read vector dataset from',msg2,msg3) 
      return
    endif

    call XF_CLOSE_GROUP(gid,ierr)
    call XF_CLOSE_FILE(fid,ierr)
    
    if(ncellpoly>0)then
      call interp_vec_node2cell(vtemp,vecx,vecy)  !Map from nodes to cell-centers
    else
      call map_vec_full2active(vtemp,vecx,vecy) !Convert from full to active grid 
    endif    
      
    return
    endsubroutine readvech5   
     
!************************************************************************
    subroutine readvecsteph5(afile,apath,itsind,thrs,vecx,vecy,ierr)
! Reads a scalar located in apath from afile XMDF file 
! written by Alex Sanchez, USACE-CHL
!************************************************************************
    use size_def,only: ncellsD,ncellsfull,ncellpoly
    use geo_def, only: idmap
    use interp_lib, only: interp_vec_node2cell
    use diag_lib
    use xmdf
    use prec_def
    implicit none
    !Input/Output
    character(len=*), intent(in) :: afile,apath
    integer, intent(in) :: itsind !Time step index
    real(ikind), intent(out) :: vecx(ncellsD),vecy(ncellsD),thrs
    integer, intent(out) :: ierr
    !Internal Variables
    integer :: fid,gid,ntimes
    real(8),allocatable :: timesd(:) !Output times
    real(4) :: vtemp(ncellsfull*2) !Must be single
    character(len=200) :: msg2,msg3,msg4

    call XF_OPEN_FILE(trim(afile),READONLY,fid,ierr)        
    if(ierr<0)then
      call diag_print_error('Could not open file: ',afile)
    endif
              
    call XF_OPEN_GROUP(fid,trim(apath),gid,ierr)
    if(ierr<0)then
      call XF_CLOSE_FILE(fid,ierr)
      !write(msg2,*,iostat=ierr) '  File: ',trim(afile)
      !write(msg3,*,iostat=ierr) '  Path: ',trim(apath)
      !call diag_print_warning('Could not open dataset from',msg2,msg3)
      ierr = -2 !Could not open group
      return
    endif
               
    call XF_READ_VECTOR_VALUES_TIMESTEP(gid,itsind,ncellsfull,2,vtemp,ierr)
    if(ierr<0)then
      call XF_CLOSE_GROUP(gid,ierr)
      call XF_CLOSE_FILE(fid,ierr)
      write(msg2,*,iostat=ierr) '  File: ',trim(afile)
      write(msg3,*,iostat=ierr) '  Path: ',trim(apath)
      write(msg4,*,iostat=ierr) '  Time step: ',itsind
      call diag_print_warning('Could not vector time step values',msg2,msg3)
      ierr = 3 !Could not read timestep
      return
    endif
    
    call XF_GET_DATASET_NUM_TIMES(gid,ntimes,ierr)
    allocate(timesd(ntimes))
    call XF_GET_DATASET_TIMES(gid,ntimes,timesd,ierr)
    if(ierr<0)then
      thrs = -999.0
      write(msg2,*,iostat=ierr) '  File: ',trim(afile)
      write(msg3,*,iostat=ierr) '  Path: ',trim(apath)
      write(msg4,*,iostat=ierr) '  Time step: ',itsind
      call diag_print_warning('Could not read time stamp from',msg2,msg3,msg4)
      ierr = 4
    else
      thrs = timesd(itsind)
    endif
    
    call XF_CLOSE_GROUP(gid,ierr)
    call XF_CLOSE_FILE(fid,ierr)
            
    if(ncellpoly>0)then
      call interp_vec_node2cell(vtemp,vecx,vecy)  !Map from nodes to cell-centers
    else
      call map_vec_full2active(vtemp,vecx,vecy) !Convert from full to active grid 
    endif
      
    return
    endsubroutine readvecsteph5

!********************************************************************************
      subroutine readveclasth5(afile,apath,ntimes,vecx,vecy,reftimed,thrs,ierr)
! Reads a scalar located in apath from afile XMDF file 
! written by Alex Sanchez, USACE-CHL
!********************************************************************************
    use size_def, only: ncellsD,ncellsfull,ncellpoly
    use geo_def, only: idmap
    use interp_lib, only: interp_vec_node2cell
    use diag_def
    use diag_lib
    use prec_def
    use xmdf
    implicit none
    !Input/Output
    character(len=*), intent(in) :: afile,apath
    real(ikind), intent(out) :: vecx(ncellsD),vecy(ncellsD),thrs
    real(8), intent(out) :: reftimed            !Reference time
    integer, intent(out) :: ierr,ntimes
    !Internal Variables
    integer :: fid,gid    
    real(8), allocatable :: timed(:)  !Output times
    real(4) :: vtemp(ncellsfull*2) !Must be single

    call XF_OPEN_FILE(trim(afile),READONLY,fid,ierr)        
    if(ierr<0)then
      ierr = -1 !Could not open file
      call diag_print_error('Could not open file: ',afile)
    endif
      
    call XF_OPEN_GROUP(fid,trim(apath),gid,ierr)
    if(ierr<0)then
      call XF_CLOSE_FILE(fid,ierr)  
      !write(msg2,*,iostat=ierr) '  File: ',trim(afile)
      !write(msg3,*,iostat=ierr) '  Path: ',trim(apath)
      !call diag_print_warning('Could not open dataset from',msg2,msg3)
      ierr = -2
      return
    endif
      
    call XF_GET_DATASET_NUM_TIMES(gid,ntimes,ierr)
    if(ierr<0)then
      call XF_CLOSE_GROUP(gid,ierr)
      call XF_CLOSE_FILE(fid,ierr)  
      write(msg2,*,iostat=ierr) '  File: ',trim(afile)
      write(msg3,*,iostat=ierr) '  Path: ',trim(apath)
      call diag_print_error('Could not read number of times from ',msg2,msg3)
    endif
    
    call XF_READ_VECTOR_VALUES_TIMESTEP(gid,ntimes,ncellsfull,2,vtemp,ierr)
    if(ierr<0)then
      call XF_CLOSE_GROUP(gid,ierr)
      call XF_CLOSE_FILE(fid,ierr)  
      return
    endif
      
!    string = trim(apath) // 'TIME'
!    call XF_OPEN_GROUP(fid,trim(string),gid,ierr)
!    if(ierr<0)then
!      call diag_print_error('Invalid dataset path: ',apath)
!    endif    
!    call XF_READ_SCALAR_VALUES_TIMESTEP(gid,1,1,thrs,ierr)
        
    allocate(timed(ntimes))
    call XF_GET_DATASET_TIMES(gid,ntimes,timed,ierr)   
    if(ierr<0)then
      thrs = -999.0
      write(msg2,*,iostat=ierr) '  File: ',trim(afile)
      write(msg3,*,iostat=ierr) '  Path: ',trim(apath)
      write(msg4,*,iostat=ierr) '  Time step: ',ntimes
      call diag_print_warning('Could not read times from ',msg2,msg3,msg4)
      ierr = 4
    else
      thrs = timed(ntimes)
    endif
    deallocate(timed)
    
    call XF_GET_DATASET_REFTIME(gid,reftimed,ierr)
    if(ierr<0)then
      reftimed = -999.0d0
      write(msg2,*,iostat=ierr) '  File: ',trim(afile)
      write(msg3,*,iostat=ierr) '  Path: ',trim(apath)
      call diag_print_warning('Could not read reference time from',msg2,msg3)
      ierr = 5
    endif
      
    call XF_CLOSE_GROUP(gid,ierr)
    call XF_CLOSE_FILE(fid,ierr)
            
    if(ncellpoly>0)then
      call interp_vec_node2cell(vtemp,vecx,vecy)  !Map from nodes to cell-centers
    else
      call map_vec_full2active(vtemp,vecx,vecy) !Convert from full to active grid 
    endif
      
    return
    endsubroutine readveclasth5

!********************************************************************************
      subroutine readvectimeh5(afile,apath,thrs,vecx,vecy,reftimed,ierr)
! Reads a scalar located in apath from afile XMDF file 
! written by Alex Sanchez, USACE-CHL
!********************************************************************************
    use size_def, only: ncellsD,ncellsfull,ncellpoly
    use geo_def, only: idmap
    use interp_lib, only: interp_vec_node2cell
    use diag_def
    use diag_lib
    use prec_def
    use xmdf
    implicit none
    !Input
    character(len=*), intent(in) :: afile,apath
    real(ikind), intent(in) :: thrs
    !Output
    real(ikind), intent(out) :: vecx(ncellsD),vecy(ncellsD)
    real(8), intent(out) :: reftimed            !Reference time
    integer, intent(out) :: ierr
    !Internal Variables
    integer :: i,fid,gid,nstep,ntimes
    real(8), allocatable :: timed(:)  !Output times
    real(8) :: thrsd,terrd,terrdmin
    real(4) :: vtemp(ncellsfull*2) !Must be single

    call XF_OPEN_FILE(trim(afile),READONLY,fid,ierr)        
    if(ierr<0)then
      ierr = -1 !Could not open file
      call diag_print_error('Could not open file: ',afile)
    endif
      
    call XF_OPEN_GROUP(fid,trim(apath),gid,ierr)
    if(ierr<0)then
      call XF_CLOSE_FILE(fid,ierr)  
      !write(msg2,*,iostat=ierr) '  File: ',trim(afile)
      !write(msg3,*,iostat=ierr) '  Path: ',trim(apath)
      !call diag_print_warning('Could not open dataset from',msg2,msg3)
      ierr = -2
      return
    endif
      
    call XF_GET_DATASET_NUM_TIMES(gid,ntimes,ierr)
    if(ierr<0)then
      call XF_CLOSE_GROUP(gid,ierr)
      call XF_CLOSE_FILE(fid,ierr)  
      write(msg2,*,iostat=ierr) '  File: ',trim(afile)
      write(msg3,*,iostat=ierr) '  Path: ',trim(apath)
      call diag_print_error('Could not read number of times from ',msg2,msg3)
    endif
      
!    string = trim(apath) // 'TIME'
!    call XF_OPEN_GROUP(fid,trim(string),gid,ierr)
!    if(ierr<0)then
!      call diag_print_error('Invalid dataset path: ',apath)
!    endif    
!    call XF_READ_SCALAR_VALUES_TIMESTEP(gid,1,1,thrs,ierr)
        
    allocate(timed(ntimes))
    call XF_GET_DATASET_TIMES(gid,ntimes,timed,ierr)   
    if(ierr<0)then
      write(msg2,*,iostat=ierr) '  File: ',trim(afile)
      write(msg3,*,iostat=ierr) '  Path: ',trim(apath)
      write(msg4,*,iostat=ierr) '  Time step: ',ntimes
      call diag_print_warning('Could not read times from ',msg2,msg3,msg4)
      ierr = 4
    else
      nstep = 1
      terrdmin = 1.0e9
      thrsd = dble(thrs)
      do i=1,ntimes
        terrd = abs(thrsd-timed(i))
        if(terrd<terrdmin)then
          nstep = i
          terrdmin = terrd
        endif
      enddo
      if(terrdmin>0.001)then
        write(msg2,*,iostat=ierr) '  File: ',trim(afile)
        write(msg3,*,iostat=ierr) '  Path: ',trim(apath)
        write(msg4,*,iostat=ierr) '  Time: ',thrs,' hrs'
        call diag_print_warning('Could not find time: ',msg2,msg3,msg4)
      endif
    endif
    deallocate(timed)
    
    call XF_READ_VECTOR_VALUES_TIMESTEP(gid,nstep,ncellsfull,2,vtemp,ierr)
    if(ierr<0)then
      call XF_CLOSE_GROUP(gid,ierr)
      call XF_CLOSE_FILE(fid,ierr)  
      return
    endif
    
    call XF_GET_DATASET_REFTIME(gid,reftimed,ierr)
    if(ierr<0)then
      reftimed = -999.0d0
      write(msg2,*,iostat=ierr) '  File: ',trim(afile)
      write(msg3,*,iostat=ierr) '  Path: ',trim(apath)
      call diag_print_warning('Could not read reference time from',msg2,msg3)
      ierr = 5
    endif
      
    call XF_CLOSE_GROUP(gid,ierr)
    call XF_CLOSE_FILE(fid,ierr)
            
    if(ncellpoly>0)then
      call interp_vec_node2cell(vtemp,vecx,vecy)  !Map from nodes to cell-centers
    else
      call map_vec_full2active(vtemp,vecx,vecy) !Convert from full to active grid 
    endif
      
    return
    endsubroutine readvectimeh5

#endif     
    
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>      
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      
!********************************************************************    
    subroutine readstep63(nnodes,timesec,wsen,ierr)
! Reads the ADCIRC Water Surface Elevation ASCII file (fort.63)
! written by Alex Sanchez, USACE-CHL
!********************************************************************
    use diag_lib
    use prec_def
    implicit none
    !Input/Output
    integer,        intent(in)  :: nnodes  !# of nodes  
    integer,        intent(out) :: ierr    !error message
    real(ikind),    intent(out) :: timesec !time in seconds
    real(ikind),    intent(out) :: wsen(nnodes) !wse at each node
    !Internal variables
    integer :: k,kdum 
    character(len=100) :: msg,msg2
    
    read(63,*,iostat=ierr) timesec
    if(ierr/=0)then
      write(msg,*,iostat=ierr)  '  timesec=',timesec
      write(msg2,*,iostat=ierr) '  ierr=',ierr  
      call diag_print_error('Could not read time stamp information in',&
          '  ADCIRC fort.63 file',msg,msg2)
    endif
    
    read(63,*,iostat=ierr) (kdum,wsen(k),k=1,nnodes)
    if(ierr/=0)then
      call diag_print_error('Could not read water level data in',&
         '  ADCIRC fort.63 file')
    endif
    
    return
    endsubroutine readstep63
    
!********************************************************************    
    subroutine readstep64(nnodes,timesec,un,vn,ierr)
! Reads the ADCIRC Water Surface Elevation ASCII file (fort.63)
! written by Alex Sanchez, USACE-CHL
!********************************************************************
    use diag_lib
    use prec_def
    implicit none
    !Input/Output
    integer,        intent(in)  :: nnodes  !Number of nodes
    integer,        intent(out) :: ierr    !error message
    real(ikind),    intent(out) :: timesec !time in seconds 
    real(ikind),    intent(out) :: un(nnodes),vn(nnodes)  !Current velocities   
    !Internal variables
    integer :: k,kdum
    character(len=100) :: msg,msg2,msg3
    
    read(64,*,iostat=ierr) timesec
    if(ierr/=0)then
      write(msg,*,iostat=ierr)  '  timesec=',timesec
      write(msg2,*,iostat=ierr) '  ierr=',ierr  
      call diag_print_error('Could not read time stamp information in',&
          '  ADCIRC fort.64 file',msg,msg2)
    endif
    
    do k=1,nnodes
      !Read velocities
      read(64,*,iostat=ierr) kdum,un(k),vn(k)
      if(ierr/=0)then
        write(msg,*,iostat=ierr)  '  k=',k    
        write(msg2,*,iostat=ierr) '  ierr=',ierr  
        call diag_print_error('Could not read velocity data in',&
          '  ADCIRC fort.64 file',msg,msg2)  
      endif
      !Check values
      if(abs(un(k))>5.0 .or. abs(vn(k))>5.0)then
        write(msg,*,iostat=ierr)  '  k=',k
        write(msg2,*,iostat=ierr) '  un(k)=',un(k)
        write(msg3,*,iostat=ierr) '  vn(k)=',vn(k)
        call diag_print_warning('Large velocities found in ',&
          '  ADCIRC fort.64 file',msg,msg2,msg3)
        !continue
      endif
    enddo
    
    return
    endsubroutine readstep64      
      
endmodule in_lib