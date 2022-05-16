!**********************************************************************
    subroutine card_scalar(inunit,defunits,tounits,scalar,ierr)
! Reads a scalar from the card file
! written by Alex Sanchez, USACE-CHL
!**********************************************************************
    use prec_def
    use diag_lib
    use unitconv_lib
    implicit none
    !Input/Output
    integer,         intent(in) :: inunit   !Input card file unit
    character(len=*),intent(in) :: defunits !Default unit string
    character(len=*),intent(in) :: tounits  !Output unit string    
    real(ikind),     intent(out) :: scalar  !Scalar value
    integer,         intent(out) :: ierr    !Error code
    !Internal
    character(len=10) :: fromunits,temp
    character         :: extra_attr
    character(len=37) :: cardname
    character(len=200) :: aline    
    real(ikind) :: scalarin
    
    scalarin = -9999.0    
    backspace(inunit)
    read(inunit,'(A)') aline
    fromunits  = ' '    
    extra_attr = ' '
!    read(aline,*,iostat=ierr) cardname,scalarin,fromunits,extra_attr
    read(aline,*,end=50) cardname,scalarin,fromunits,extra_attr
50  if(fromunits==' ' .and. abs(scalarin+9999.0)<1.0e-5)then
      call diag_print_error('Invalid Card Specification: ',&
        aline,'  Missing Card Value')
      ierr = -1
    else
      scalar = scalarin
      if (extra_attr .ne. ' ' .and. extra_attr .ne. '!') fromunits = trim(fromunits)//' '//extra_attr
      ierr = 0
    endif
    if(ierr==-1 .or. fromunits(1:1)==' ' .or. &
      fromunits(1:1)=='!' .or. fromunits(1:1)=='#')then
      fromunits = defunits !None specified in file, use default value
    endif

    !General linear unit conversion
    call unitconv_scal(fromunits,tounits,scalar)
    
    return
    end subroutine card_scalar

!*****************************************************************
    subroutine card_boolean(inunit,bool,ierr)
! Reads a boolean (true or false) from the card file
! written by Alex Sanchez, USACE-CHL
!*****************************************************************    
    use diag_lib
    implicit none
    !Input/Output    
    integer, intent(in) :: inunit
    logical, intent(inout) :: bool
    integer, intent(out) :: ierr
    !Internal
    character(len=37) :: cardname,cdum    
    
    ierr = 0
    backspace(inunit)
    read(inunit,*,iostat=ierr) cardname, cdum     
    if(ierr/=0)then
      ierr = -1
      call diag_print_warning('Could not read CMS-Flow Control file card: ',cardname)
    endif
    call uppercase(cdum)
    if(cdum(1:2)=='ON')then
      bool = .true.
    elseif(cdum(1:3)=='OFF')then
      bool = .false.
    else
      ierr = -2
      call diag_print_warning('Invalid value for card: ',cardname,&
        '  Value must be either ON or OFF')
    endif
        
    return
    end subroutine card_boolean
    
!!*****************************************************************************
!    subroutine card_file(inunit,defaultfile,defaultpath,filename)
!!*****************************************************************************
!    use diag_lib
!    implicit none
!    !Input/Output
!    integer,intent(in) :: inunit
!    character(len=*),intent(in) :: defaultfile,defaultpath
!    character(len=*),intent(inout) :: filename
!    !Internal Variables
!    integer :: ierr
!    character(len=37)  :: cardname
!    character(len=200) :: apath,aname
!    character(len=10)  :: aext
!    logical :: foundfile
!    
!    backspace(inunit)
!    read(inunit,*,iostat=ierr) cardname, filename
!    
!    if(ierr/=0)then
!      filename = trim(defaultpath) // trim(defaultfile)
!    endif
!    
!    call fileparts(filename,apath,aname,aext)  
!        
!    if(apath(1:2)=='.')then !Relative path
!        
!    if(apath(1:1)=='.')then !Relative path
!        
!    elseif(len_trim(apath)==0)then
!      filename = trim(defaultpath) // filename  
!    endif
!    
!    return
!    end subroutine card_file
    
!*****************************************************************************
    subroutine card_dataset(inunit,defaultfile,defaultpath,datafile,datapath,ndim)
! Reads a dataset card from the CMS conrol/card file
!
! written by Alex Sanchez, USACE-CHL
!*****************************************************************************
    use diag_lib
    implicit none
    !Input/Output
    integer,intent(in) :: inunit
    character(len=*),intent(in) :: defaultfile,defaultpath
    character(len=*),intent(inout) :: datafile,datapath
    integer, intent(in) :: ndim
    !Internal Variables
    character(len=10)  :: aext
    character(len=37)  :: cardname
    character(len=200) :: apath,aname
    character(len=400) :: aline,afile
    integer :: ierr, ival,nquotes,aNumber
    logical :: foundfile
! added 5/21/2018
    integer :: ival1,ival2,ilen
    logical :: isMP, isLatLong
    
    backspace(inunit)
    read(inunit,'(A)') aline
    read(aline,*,iostat=ierr) cardname, datafile    
    if(ierr/=0)then
      call diag_print_error('Could not read dataset card: ',cardname)
    endif
    call fileparts(datafile,apath,aname,aext)

    ilen=len_trim(aPath)-1
    if (ilen.gt.0) then
#ifdef _WIN32
      aFile = trim(aPath(1:ilen))//'\'//trim(aName)//'.'//trim(aExt)
#else
      aFile = trim(aPath(1:ilen))//'/'//trim(aName)//'.'//trim(aExt)
#endif
      datafile = trim(aFile)
    endif

    if(len_trim(aext)==0)then !Is NOT a file. It is a path
      datapath = datafile
      datafile = defaultfile !Assume the file is the grid file
      return  
    endif            
    call lowercase(aext)
    select case(aext)
    case('h5')
      call countquotes(aline,nquotes)        !Account for lacking final quote on Bed Layer cards
      if(mod(nquotes,2).eq.1) then
        aline=trim(aline)//' "'
      endif
      read(aline,*,iostat=ierr) cardname, datafile, datapath 
      
      !Ignore hotstarting file  MEB  1/15/2021
      if(ierr/=0 .and. cardname(1:20)/='INITIAL_STARTUP_FILE' .and. cardname(1:22) /= 'INITIAL_CONDITION_FILE')then
        call diag_print_warning('Path not specified for dataset card: ',cardname)
      endif
      
      !If boundary condition dataset, do not add to list
      isMP = .false. ; isLatLong = .false.   !some modifications to this area 5/21/2018
      ival=-1
      ival = INDEX(datafile,"mp")
      if (ival .gt. 0) then
        isMP = .true.
        ival1 = -1; ival2= -1
        ival1 = INDEX(datapath,"Lats")
        ival2 = INDEX(datapath,"Lons")
        if (ival1 .gt. 0 .or. ival2 .gt. 0) isLatLong=.true.
      endif
      
      if (ival .le. 0) then
        if(datapath(1:10).eq.'Datasets/ ') then                           
          continue                                               !If this happens, most likely it is a blank Dxx_DATASET from the Bed Layer Block that SMS writes.  Ignore and continue.
        else
          call addXMDFDataset2List (datafile,datapath,ndim)      !Add file and path to a running list of all XMDF files for output to ASCII later, if requested.
        endif  
!added 5/21/2018
      elseif (isLatLong) then
          call addXMDFDataset2List (datafile,datapath,ndim)      !Add file and path to a running list of all XMDF files for output to ASCII later, if requested.
!
      endif
    case('bid')  
      read(aline,*,iostat=ierr) cardname, datafile, datapath 
      if(len_trim(datapath)==0) then  !Path may need to be read as an integer into a string
        read(aline,*,iostat=ierr) cardname, datafile, aNumber
        write(datapath,'(I1)') aNumber      !Temporarily one digit limitation.  If needed, can do more.
      endif
      
    case('txt')
      !Set some text to the path to make known there is a latitude file
      datapath='ASCII_Input'
    end select
    
    inquire(file=datafile,exist=foundfile)       
    if(.not.foundfile)then !If not found, try adding path
      datafile  = trim(defaultpath)//datafile
      inquire(file=datafile,exist=foundfile)
      if(.not.foundfile)then
        call diag_print_error('Could not find file: ',datafile)  
      endif  
    endif
    
    return
    end subroutine card_dataset
    
!*****************************************************************************
    subroutine addUnknownCard2List (iFile)
! Adds unknown cards to a list for viewing later.
!
! written by Mitchell Brown, USACE-CHL - 09/04/2019
!*****************************************************************************
    use cms_def, only: cardList,nCards
    implicit none
    integer,intent(in) :: iFile
    
    integer :: i,ilen
    logical :: exists = .false.
    character(len=200),allocatable :: oldCards(:), oldVals(:)
    character(len=200) :: aCard,aValue
    
    backspace(iFile)
    read(iFile,'(A200)') aValue
    ilen=index(aValue,' ')
    aCard=aValue(1:ilen)
    do i=1,ilen
      aValue(i:i)=' '
    enddo
    aValue = adjustl(aValue)
    
    !increase total number of cards in list by 1 and move old information to new array
    if(nCards==0)then
      allocate(cardList(1))
      nCards=1
    else
      !check to see if card is already in list
      do i=1,nCards
        if(cardList(i)%cardname == trim(aCard)) exists = .true.
        if(exists) return
      enddo
      
      !Move old information to temp array
      allocate(oldCards(ncards),oldVals(ncards))
      do i=1,nCards
        oldCards(i) = trim(cardList(i)%cardname)
        oldVals(i)  = trim(cardList(i)%value)
      enddo
      !Increase allocation by 1
      deallocate(cardList) ; allocate(cardList(nCards+1))
      !Move information back from temp array
      do i=1,nCards
        cardList(i)%cardname = trim(oldCards(i))
        cardList(i)%value    = trim(oldVals(i))
      enddo
      nCards=nCards+1
    endif
    cardList(nCards)%cardname = trim(aCard)   !Add new card
    cardList(nCards)%value    = trim(aValue)  !Add new value
    if (allocated(oldCards)) deallocate(oldCards,oldVals)
    
    return
    end subroutine addUnknownCard2List

!*****************************************************************************
    subroutine addXMDFDataset2List (datafile,datapath,ndim)
! Reads a dataset card from the CMS conrol/card file
!
! written by Mitchell Brown, USACE-CHL - 02/22/2018
!*****************************************************************************
    use cms_def, only: dsetList, ndsets
    implicit none
    character(len=*),intent(in) :: datafile,datapath
    integer, intent(in) :: ndim

    integer :: i
    character(len=200),allocatable :: oldfiles(:),oldpaths(:)
    integer, allocatable :: olddims(:)
!  added 5/21/2018
    logical :: exists
    exists = .false.

    !increase total number of datafiles in list by 1 and move old information to new array
    if (ndsets==0) then
      allocate(dsetList(1))
      ndsets = 1
    else
      !check to see if dataset is already in list - 05/21/2018
      do i=1,ndsets
        if (dsetList(i)%path == trim(datapath)) exists = .true.
        if (exists) return
      enddo
        
      !Move old information to temp array     
      allocate(oldfiles(ndsets),oldpaths(ndsets),olddims(ndsets))
      do i=1,ndsets                     
        oldfiles(i) = trim(dsetList(i)%filename)
        oldpaths(i) = trim(dsetList(i)%path)
        olddims(i)  = dsetList(i)%ndim
      enddo
      !Increase allocation by 1
      deallocate (dsetList)
      allocate (dsetList(ndsets+1))
      !Move information back from temp array     
      do i=1,ndsets                     
        dsetList(i)%filename = trim(oldfiles(i))
        dsetList(i)%path     = trim(oldpaths(i))
        dsetList(i)%ndim     = olddims(i)
      enddo
      ndsets = ndsets + 1
    endif
    dsetList(ndsets)%filename = trim(datafile)  !Add new file
    dsetList(ndsets)%path     = trim(datapath)  !Add new path
    dsetList(ndsets)%ndim     = ndim            !Add new dimension
    if (allocated(oldfiles)) deallocate(oldfiles,oldpaths)
    
    return
    end subroutine addXMDFDataset2List

!***************************************************************************    
    subroutine card_datetime(inunit,iyr,imo,iday,ihr,imin,isec)
! Reads a date and time card using the format (ISO Standard)
!   STARTING_DATE_TIME  YYYY-MM-DD HH:MM:SS UTC
!
! written by Alex Sanchez, USACE-CHL
!***************************************************************************   
    use diag_lib
    use prec_def
    implicit none
    !Input/Output
    integer,intent(in) :: inunit
    integer,intent(out) :: iyr,imo,iday,ihr,imin,isec
    !Internal variables
    integer :: iloc,ierr
    character(len=37) :: cdum
    character(len=200) :: aline
    
    ierr = 0
    backspace(inunit)
    read(inunit,'(A)') aline
    
784 format(I4,'-',I2,'-',I2) !YYYY-MM-DD
    iloc = index(aline,'-',back=.false.)
    cdum = aline(iloc-4:)
    read(cdum,784,iostat=ierr) iyr,imo,iday
    
    if(ierr/=0)then
      call diag_print_error('Could not read date: ',aline)
    endif
    
673 format(I2,':',I2,':',I2) !hh:mm:ss    
    iloc = index(aline,':',back=.false.)
    cdum = aline(iloc-2:)
    read(cdum,673,iostat=ierr) ihr,imin,isec
    
    if(ierr/=0)then
      call diag_print_error('Could not read date: ',cdum)
    endif
        
    return
    end subroutine card_datetime
    
!*****************************************************************    
    subroutine card_bid(defpath,bidfile,bidpath,idnum)
! Reads a boundary ID card for example:
! CELLSTRING "flow_mp.h5" " "PROPERTIES/Model Params/Boundary_#1" 
! written by Alex Sanchez, USACE-CHL
!*****************************************************************
    use diag_lib
    implicit none
    !Input
    character(len=*),intent(in)    :: defpath  !Default input path
    character(len=*),intent(inout) :: bidfile  !Boundary ID File (*.h5, *.2dm, *,bid) (including directory path)
    character(len=*),intent(inout) :: bidpath  !Boundary ID Path (for XMDF files)
    integer,         intent(out)   :: idnum    !Boundary ID number within bid file
    !Internal Variables
    character(len=10)  :: aext
    character(len=37)  :: cardname
    character(len=200) :: aname,apath
    character(len=300) :: aFile
    integer            :: ilen
    logical            :: foundfile
    
    !Read file and path
    backspace(77)
    read(77,*) cardname,bidfile,bidpath
    call fileparts(bidfile,apath,aname,aext)
    
    ilen=len_trim(aPath)-1
    if(ilen .gt. 0) then
#ifdef _WIN32
      aFile = trim(aPath(1:ilen))//'\'//trim(aName)//'.'//trim(aExt)
#else
      aFile = trim(aPath(1:ilen))//'/'//trim(aName)//'.'//trim(aExt)
#endif
      bidfile = trim(aFile)
    endif

    !Check file exists and guess path if necessary
    inquire(file=bidfile,exist=foundfile)  
    if(.not.foundfile .and. len_trim(apath)==0)then
      bidfile = trim(defpath) // trim(bidfile)
      inquire(file=bidfile,exist=foundfile) 
    endif
    if(.not.foundfile)then
      call diag_print_error('Could not find Boundary ID file: ',bidfile)
    endif  
    
    call lowercase(aext)
    select case(aext)
    case('h5')
      call bndpath2id(bidpath,idnum)
    case('2dm','bid')
      call bndpath2id(bidpath,idnum)
      !read(bidpath,*) idnum
    case default
      call diag_print_error('Unsupported Boundary ID File type: ',bidpath,&
        '  File must be *.h5, *.2dm, or *.bid')
    end select
    
    return
    end subroutine card_bid    