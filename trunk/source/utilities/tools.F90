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
    !Determmin path
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