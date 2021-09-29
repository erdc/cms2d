!******************************************************
    subroutine CMS_Tools_Dialog    
! Perfom one of a selection of internal tools
! written by Mitchell Brown, USACE-CHL
! 04/14/2021
!******************************************************    
    use diag_lib, only: diag_print_message
    implicit none
    integer :: ichoice 
    
50  continue
#ifdef _WIN32
    call system('cls')
#else
    call system('clear')
#endif
    write(*,*) ''
    write(*,*) '*****************************************************************'
    write(*,*) 'Make a choice from the following CMS Tools'
    write(*,*) '1 - Print RT_JULIAN/REFTIME corresponding to actual date'
    write(*,*) '2 - Print actual date corresponding to RT_JULIAN/REFTIME'
    write(*,*) '3 - Read WSE Dataset and output a Max WSE dataset for all times'
    write(*,*) '9 - Exit'
    write(*,*) '*****************************************************************'
    write(*,*) ''
100 read (*,*) ichoice
    
    select case (ichoice)
    case (1) 
      call CMS_date2reftime
    case (2) 
      call CMS_reftime2date
    case (3)
      call CMS_MaxWSE
    case (9)
      STOP
    case default
      call diag_print_message ('Selection not available, make another selection')
      goto 100
    end select
    
    goto 50
    
    contains
      subroutine CMS_MaxWSE
        use XMDF
        use out_lib, only: OPEN_CREATE_DATASET
        use EXP_Global_def, only: RROUND
        use comvarbl, only: reftime
        use diag_lib, only: diag_print_error
        
        implicit none
        integer :: fid, gid, ierr, a_Number, a_MaxPathLength, iloc
        integer :: nfid, ndid
        integer :: i, j, a_NumTimes, nCells, complete
        real(4) :: Num
        real(4), allocatable :: Values(:), max_array(:)
        real(8), allocatable :: Times(:)
        character(100) :: filename, newfile, a_Path
        character(100), allocatable :: Paths(:)
        logical :: exists

        do
          write(*,'(A)') 'Enter filename of WSE solution file (*_wse.h5, *_sol.h5, etc)'
          read(*,*) filename
          inquire(FILE=trim(filename),EXIST=exists)
          if (.not.exists) then
            write(*,'(A)') 'File not found, try again'
          else
            exit
          endif
        enddo

        call XF_OPEN_FILE(trim(filename),readwrite,fid,ierr) !Open XMDF file
        if(ierr<0)then
          write(*,'(A)') 'Error opening file'
        endif          
        call XF_GET_SCALAR_DATASETS_INFO (fid, a_Number, a_MaxPathLength, ierr)
        allocate(Paths(a_Number))
        do i=1,a_Number
          call XF_GET_SCALAR_DATASET_PATHS (fid, i, a_MaxPathLength, Paths, ierr)
          a_Path = Paths(i)
          iloc = index(a_Path,'Water_Elevation')
          if (iloc <= 0) then
            cycle
          else
            a_Path(iloc+15:100) = ' '
            exit
          endif
        enddo
        deallocate(Paths)
        
        if (i > a_Number) then
          write (*,'(A)') 'No Water Elevation Dataset found'
          stop
        endif 
        
        call OPEN_CREATE_DATASET(fid,trim(a_Path),gid,1,'m',ierr)
        
        !Parse through the datasets in the file and get max value
        call XF_READ_REFTIME(gid, Reftime, ierr)
        call XF_GET_DATASET_NUM_TIMES(gid, a_NumTimes, ierr)
        allocate(Times(a_NumTimes))
        times = 0.0 ;
        
        call XF_GET_DATASET_TIMES(gid, a_NumTimes, Times, ierr)
        call XF_GET_DATASET_NUMVALS(gid, nCells, ierr)
        allocate(max_array(nCells),Values(nCells))
        max_array = -1000.0 ; values = 0.0
        complete = 0
        write(*,*) ' '
        write(*,'(A)') ' Parsing existing Water Elevations'
        do i=1,a_NumTimes
          call XF_READ_SCALAR_VALUES_TIMESTEP(gid, i, nCells, Values, ierr)
          do j=1,nCells
            max_array(j) = max(Values(j),max_array(j))
          enddo
          num = real(i)/real(a_NumTimes)*100
          if (num .ge. complete) then           
            write(*,'(2x,I3,x,A)') complete, 'Percent Complete'
            complete = complete + 10
          endif
        enddo
        call XF_CLOSE_GROUP(gid,ierr)  !Close old dataset    
        call XF_CLOSE_FILE(fid,ierr)   !Close old file
        
        !Create new file for the Max Water Elevations based on original filename
        iloc = index(filename,'_wse')
        newfile = filename(1:iloc-1)//'_maxwse.h5'
        iloc = index(a_Path,'Water_Elevation')
        a_Path = a_Path(1:iloc-1)//'Max_Water_Elevation'

        inquire(FILE=trim(newfile), EXIST=exists)
        if (exists) then
          open(100,file=newfile,iostat=ierr)
          close(100,status='delete',iostat=ierr)
        endif            
        call XF_CREATE_FILE(trim(newfile),readwrite,nfid,ierr)
        call OPEN_CREATE_DATASET(nfid,trim(a_Path),ndid,1,'m',ierr)
        if (ierr.ge.0) then
          write(*,*) ' '
          write(*,'(A,A)')" Writing Maximum elevations to file: '"//trim(newfile)//"'"
          write(*,*) ' '
        else
          call diag_print_error('Error creating dataset')
        endif
        call XF_WRITE_SCALAR_TIMESTEP(ndid,0.d0,nCells,max_array,ierr)
        
        call XF_CLOSE_GROUP(ndid,ierr) !Close new dataset
        call XF_CLOSE_FILE(nfid,ierr)  !Close new file 
        write(*,*) 'Press the Enter key to continue'
        read(*,*) 
        
        return
      end subroutine CMS_MaxWSE
    
      subroutine CMS_date2reftime
        use XMDF
        use comvarbl, only: iyr, imo, iday, ihr, imin, isec, reftime
        implicit none
        integer :: ierr
        
        
        write(*,'(A)') 'Enter date (yyyy-mm-dd hh:mm:ss)'
        call card_datetime(5,iyr,imo,iday,ihr,imin,isec) 
        call XF_CALENDAR_TO_JULIAN(0,iyr,imo,iday,ihr,imin,isec,reftime,ierr)
        write(*,*)''
        write(*,'(A,F0.6)') 'RT_JULIAN date is: ',reftime
        write(*,'(A)') 'Press any key to continue.'
        read(*,*)
        write(*,*)''
        
        return
      end subroutine CMS_date2reftime
    
      subroutine CMS_reftime2date
        use XMDF
        use comvarbl, only: iyr, imo, iday, ihr, imin, isec, reftime
        implicit none
        integer :: ierr, era = 0

100     format('Calendar date is: ',i4.4,'-',i2.2,'-',i2.2, 1x ,i2.2,':',i2.2,':',i2.2)
        
        write(*,*) 'Enter reference time:'
        read (*,*) reftime 
        call XF_JULIAN_TO_CALENDAR(era,iyr,imo,iday,ihr,imin,isec,reftime,ierr)
        write(*,*)''
        write(*,100) iyr, imo, iday, ihr, imin, isec
        write(*,'(A)') 'Press any key to continue.'
        read(*,*)
        write(*,*)''
        
        return
      end subroutine CMS_reftime2date

    end subroutine CMS_Tools_Dialog
!******************************************************    
    
!******************************************************
    logical function findCard(aFile,aCard,aValue)    
! Find a given card in the cardfile and return the rest of the selected line in the variable, aValue.
! written by Mitch Brown, USACE-CHL 06/25/2021
!
! The correct interface is given below
    !interface  
    ! function findCard(aFile,aCard,aValue)
    !    character(len=*),intent(in)    :: aFile
    !    character(len=*),intent(in)    :: aCard
    !    character(len=100),intent(out) :: aValue
    !    logical :: findCard
    !  end function
    !end interface
!******************************************************    
    implicit none
    character(len=*), intent(in)    :: aFile
    character(len=*), intent(in)    :: aCard
    character(len=100), intent(out) :: aValue
    
    character(len=80)  :: testCard
    character(len=100) :: aLine
    logical            :: foundCard = .false.
    integer            :: iloc
    
    open(77,file=aFile,status='unknown')
    do 
      read(77,*,end=100) testCard
      if (trim(testCard) == trim(aCard)) then
        backspace(77)
        read(77,'(A100)') aLine
        iloc = index(aLine,' ')
        aLine = adjustL(aLine(iloc:))
        read(aLine, '(A100)') aValue
        findCard = .true.
        return
        exit
      else
        findCard = .false.
      endif
    enddo

100 findCard = .false.

    return
    end function findcard
!******************************************************    
   
    
!******************************************************
    subroutine fileparts(astr,apath,aname,aext)    
! Determines the parts of a file name
! written by Alex Sanchez, USACE-CHL
!******************************************************    
    implicit none
    !Input/Output
    character(len=*),intent(in) :: astr
    character(len=*),intent(inout) :: apath,aname,aext
    !Internal
    integer :: i,k,nn    
    
    nn = len_trim(astr)    
    !Determine path
    apath = ''   
    do k=nn,1,-1
      if(astr(k:k)=='\' .or. astr(k:k)=='/')then
        apath = astr(1:k)     
        exit
      endif
    enddo
    
    !Determine name and extension
    aname = astr(max(k+1,1):nn)
    aext = ''    
    do i=nn,k+1,-1
      if(astr(i:i)=='.')then
        aext = astr(i+1:nn)    
        aname = trim(astr(max(k+1,1):i-1))
        return
      endif
    enddo          

    return
    endsubroutine fileparts
    
!!******************************************************
!    subroutine fileparts(astr,apath,aname,aext)    
!! Determines the parts of a file name
!! written by Alex Sanchez, USACE-CHL
!!******************************************************    
!    implicit none
!    !Input/Output
!    character(len=*),intent(in) :: astr
!    character(len=*),intent(inout) :: apath
!    character(len=*),intent(inout),optional :: aname,aext
!    !Internal
!    integer :: i,k,nn    
!    
!    nn = len_trim(astr)    
!    !Determmin path
!    apath = ''   
!    do k=nn,1,-1
!      if(astr(k:k)=='\' .or. astr(k:k)=='/')then
!        apath = astr(1:k)     
!        exit
!      endif
!    enddo    
!    
!    if(.not.present(aname)) return
!    
!    !Determine name and extension
!    aname = astr(max(k+1,1):nn)
!    if(present(aext)) aext = ''    
!    do i=nn,k+1,-1
!      if(astr(i:i)=='.')then
!        if(present(aext)) aext = astr(i+1:nn)    
!        aname = trim(astr(max(k+1,1):i-1))
!        return
!      endif
!    enddo          
!
!    return
!    endsubroutine fileparts
    
!******************************************************
    subroutine filepath(astr,apath)  
! Determines the path of a file
! written by Alex Sanchez, USACE-CHL
!******************************************************    
    implicit none
    !Input/Output
    character(len=*),intent(in) :: astr
    character(len=*),intent(inout) :: apath
    !Internal
    integer :: k,nn    
    
    nn = len_trim(astr)
    apath = ''   
    do k=nn,1,-1
      if(astr(k:k)=='\' .or. astr(k:k)=='/')then
        apath = astr(1:k)     
        return
      endif
    enddo       

    return
    endsubroutine filepath
    
!******************************************************
    subroutine fileext(astr,aext)  
! Determines the extension of a file
! written by Alex Sanchez, USACE-CHL    
!******************************************************    
    implicit none
    !Input/Output
    character(len=*),intent(in) :: astr
    character(len=*),intent(inout) :: aext
    !Internal
    integer :: i,k,nn
        
    !Determine extension
    nn = len_trim(astr)
    k = index(astr,'\')
    if(k==0) k = index(astr,'/')
    aext = ' '
    do i=nn,max(1,k+1),-1
      if(astr(i:i)=='.')then
        aext = astr(i+1:nn)        
        exit
      endif
    enddo    
    
    return
    endsubroutine fileext
    
!******************************************************
    subroutine remove_underscores(astr)    
!******************************************************    
    implicit none
    character(len=*) :: astr
    integer :: i,nn
    
    nn = len_trim(astr)
    do i=1,nn
      if(astr(i:i)=='_')then
        astr(i:i) = ' '
      endif
    enddo    
    
    return
    endsubroutine remove_underscores

!********************************************************
    subroutine uppercase(str)
!********************************************************
    implicit none
    character(len=*), intent(inout) :: str
    integer :: i,del

    del = iachar('a') - iachar('A')
    do i=1,len_trim(str)
      if(lge(str(i:i),'a') .and. lle(str(i:i),'z'))then
        str(i:i) = achar(iachar(str(i:i)) - del)
      endif
    enddo

    return
    endsubroutine uppercase
    
!********************************************************
    function toUpper(str) result(aString)
!********************************************************
    implicit none
    character(len=*), intent(in)  :: str
    character(len=len(str))       :: aString
    integer :: i,del

    del = iachar('a') - iachar('A')
    do i=1,len_trim(str)
      if(lge(str(i:i),'a') .and. lle(str(i:i),'z'))then
        aString(i:i) = achar(iachar(str(i:i)) - del)
      else
        aString(i:i) = str(i:i)
     endif
    enddo

    return
    end function toUpper


!********************************************************
    subroutine lowercase(str)
!********************************************************
    implicit none
    character(len=*), intent(inout) :: str
    integer :: i,del

    del = iachar('a') - iachar('A')
    do i=1,len_trim(str)
      if(lge(str(i:i),'A') .and. lle(str(i:i),'Z'))then
        str(i:i) = achar(iachar(str(i:i)) + del)
      endif
    enddo

    return
    endsubroutine lowercase

!********************************************************
    function toLower(str) result(aString)
!********************************************************
    implicit none
    character(len=*), intent(in)  :: str
    character(len=len(str))       :: aString
    integer :: i,del

    del = iachar('a') - iachar('A')
    do i=1,len_trim(str)
      if(lge(str(i:i),'A') .and. lle(str(i:i),'Z'))then
        aString(i:i) = achar(iachar(str(i:i)) + del)
      else
        aString(i:i) = str(i:i)
      endif
    enddo

    return
    end function toLower

    
!*************************************************************
    function trimspace(str)
! Return index of first space
!*************************************************************    
    integer :: nn,trimspace
    character(len=*),intent(in):: str
    
    nn = len(str)      
    do while(str(nn:nn)==' ' .and. nn>0)
      nn = nn - 1
    enddo
    trimspace = nn      
    
    return
    endfunction trimspace

!*****************************************************************
    subroutine removequotes(string)
!*****************************************************************
    integer :: i,j,js,je,nn
    character(len=*),intent(inout) :: string
    
    nn=len_trim(string)
    js=1; je=nn
    do i=1,nn
      if(string(i:i)=='"')then
        do j=i+1,nn
          if(string(j:j)=='"')then
            string = string(i+1:j-1)
            return
          endif
        enddo
      endif
    enddo

    return
    endsubroutine removequotes
    
!*****************************************************************
    subroutine countquotes(string,count)
!*****************************************************************
    integer :: i,j,js,je,nn
    character(len=*),intent(in) :: string
    integer, intent(out) :: count
    
    count=0
    nn=len_trim(string)
    js=1; je=nn
    do i=1,nn
      if(string(i:i)=='"') then
        count=count+1  
      endif
    enddo

    return
    endsubroutine countquotes

!!***********************************************************************    
!    subroutine copy_file_to_temp(file,tempfile)
!!***********************************************************************    
!    !use ifport, only: system
!    use DFLIB, only: systemqq
!    use diag_def
!    use diag_lib
!    implicit none
!    !Input
!    character(len=*),intent(in) :: file
!    !Output
!    character(len=*),intent(out) :: tempfile
!    !Internal
!    integer :: ierr
!    character(len=200) :: thepath
!    character(len=200) :: thename
!    character(len=10) :: theext
!    logical :: foundfile
!    
!    !Name of temporary file
!    call fileparts(file,thepath,thename,theext)
!    tempfile = trim(thepath)//trim(thename)//'_temporary.h5'
!    !call filepath(file,thepath)
!    !tempfile = trim(thepath)//'temp.h5'
!    
!    ! Check if temp file exists, if so, delete it
!    inquire(file=tempfile,exist=foundfile)
!    if(foundfile)then
!      open(100,file=tempfile,iostat=ierr)
!      close(100,status='delete',iostat=ierr)
!    endif
!
!    ! Make a copy of ICFILE to temp.h5 and read from that.
!    !res = system('copy '//trim(file)//' '//trim(tempfile)//' > trash.txt')
!    ierr = systemqq('copy '//trim(file)//' '//trim(tempfile)//' > trash.txt')
!    inquire(file='trash.txt',exist=foundfile)
!    if(foundfile)then  !remove scratch file
!      open(100,file='trash.txt',iostat=ierr)
!      close(100,status='delete',iostat=ierr)
!    endif
!    
!    return
!    endsubroutine copy_file_to_temp
    
!************************************************************
    subroutine delete_file(file)
!************************************************************    
    implicit none
    integer :: ierr
    character(len=*),intent(in) :: file
    logical :: foundfile
    
    inquire(file=file,exist=foundfile)
    if(foundfile)then
      open(777,file=file,iostat=ierr)
      close(777,status='delete',iostat=ierr)
      !if(ierr/=0)then
        !res = systemqq('del '//trim(tfile))  
      !endif
    endif    
    
    return
    endsubroutine delete_file
    
!-----------------------------------------------------------------------
!   converts a real*4 value with n significant digits to a character 
!   string
!-----------------------------------------------------------------------
      subroutine split_real_to_integers (x, n, first, second) 

      integer*4      :: i, idot, n, lstring, first, second, loc
      real*4         :: x
      logical        :: significant_number
      character (len=20) :: string, part1, part2

      if (n < 2) n = 2
      if (n > 7) n = 7

      if (n == 2) then
          write (string, '(f0.2)') x
      else if (n == 3) then
          write (string, '(f0.3)') x
      else if (n == 4) then
          write (string, '(f0.4)') x
      else if (n == 5) then
          write (string, '(f0.5)') x
      else if (n == 6) then
          write (string, '(f0.6)') x
      else 
          write (string, '(f0.7)') x
      end if

      lstring = len_trim(string)
!
!     remove trailing zeros 
!
      idot = index (string(1:lstring),'.')
 
      do i = lstring, idot+2, -1
         if (string(i:i) == '0') then
             string(i:i) = ' '
             significant_number = .false.
         else
             significant_number = .true.
         end if
         if (significant_number) exit 
      end do

      loc=index(string,'.')
      part1 = string(:loc-1)
      part2 = string(loc+1:)

      read(part1,*) first
      read(part2,*) second
      
      return      
      end subroutine split_real_to_integers    
    
!************************************************************
    character(len=20) function Int2Str(k)
!   "Convert an integer to string."
!************************************************************
    implicit none
    integer, intent(in) :: k
    
    write (Int2Str, *) k
    Int2Str = adjustl(Int2Str)
    
    end function Int2Str
    
    module tool_def
      implicit none
  
      contains
!************************************************************
      function vstrlz(flt,gfmt) result(gbuf)
! This function will take a float variable and a format declaration, then convert it to a string. 
! The string will add a leading 0, +0, or -0 if necessary.  Also, it will remove a leading space from the PE format
      character(len=20)                      :: gbuf
      real, intent(in)                       :: flt
      character(len=*), optional, intent(in) :: gfmt
      character(len=40)                      :: gtmp
      integer                                :: istat

      gtmp = ''
      if(present(gfmt) ) then   ! specified format
        write(gtmp, gfmt, iostat=istat ) flt  
      else                      ! generic format
        write(gtmp, '(g0)', iostat=istat) flt  
      endif
      
      if( istat /= 0 ) then
        gbuf='****'
        return
      endif
      
      if    (gtmp(1:1) == '.' ) then
        gbuf = '0'//trim(gtmp)
      elseif(gtmp(1:2) == '-.' ) then
        gbuf = '-0.'//trim(gtmp(3:))
      elseif(gtmp(1:2) == '+.' ) then ! S format adds a +
        gbuf = '+0.'//trim(gtmp(3:))
      else
        gbuf = trim(adjustl(gtmp))
      endif
      
      return
      end function vstrlz
!************************************************************
    end module tool_def      