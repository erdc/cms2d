!=====================================================================
module watch_lib
! CMS Watch Library
! 
! Description:
!   Used to time sections of the code and profile the code efficiency
!
! Routines:
!   watch_default - Sets the watch defaults
!   watch_init - Initializes the watches
!   watch_print - Prints the initial parameters 
!   watch_start - Starts a watch
!   watch_stop - Stops a watch
!   watch_name2index - Returns the index of a watch
!   watch_output - Prints output 
!   watch_destroy - Deallocates and resets variables
!
! Authors:
!   Alex Sanchez, UASCE-CHL
!=====================================================================    
    implicit none
    
    private    
    public :: watch_default,watch_init,watch_start,watch_stop,&
      watch_print,watch_output,watch_destroy
    
contains   

!**********************************************************
    subroutine watch_default
! Initializes the watches    
! written by Alex Sanchez, USACE-CHL    
!**********************************************************
    use watch_def
    implicit none
    
    nwatches = 0    
    watchdef%active  =  .false.
    watchdef%running =  .false.   
    watchdef%name   = ''
    watchdef%cpu%lap =  0.0d0
    watchdef%cpu%split =  0.0d0
    watchdef%clock%lap =  0.0d0
    watchdef%clock%split = 0.0d0
    
    return
    end subroutine watch_default

!**********************************************************
    subroutine watch_init
! Initializes the watches    
! written by Alex Sanchez, USACE-CHL    
!**********************************************************
    use watch_def
    implicit none
    integer :: istat,i
    
    allocate(watch(nmaxwatches),stat=istat)
    
    do i=1,nmaxwatches
      watch(i) = watchdef
    enddo
    
    return
    end subroutine watch_init
    
!**********************************************************
    subroutine watch_print() 
! written by Alex Sanchez, USACE-CHL        
!**********************************************************
    use watch_def
    use diag_def, only: dgunit,dgfile

    implicit none
    integer :: iunit(2),i
    
222 format(' ',A,T40,A)    
    
    iunit = (/6, dgunit/)
    open(dgunit,file=dgfile,access='append') 
    
    do i=1,2        
      write(iunit(i),*) ''   
      write(iunit(i),222) ' Watches:','ON'
    enddo
    
    return
    end subroutine watch_print    
   
!**********************************************************
    subroutine watch_start(watch_name,watch_index)
!
! Usage:
!   call watch_start(watch_name)    
!   call watch_start(watch_name,watch_index)
!
! written by Alex Sanchez
!**********************************************************
    use watch_def
    use diag_lib
    use time_lib
    implicit none
    character(len=*), intent(in) :: watch_name
    integer, intent(out), optional :: watch_index
    !Internal
    integer :: iw
    character(len=200) :: msg
      
    iw = watch_name2index(watch_name)
    
    if(present(watch_index))then
      watch_index = iw  
    endif
    
    if(watch(iw)%running)then
      write(msg,*) 'Trying to starts running watch: '
      call diag_print_warning(msg,watch_name)
      return
    endif
    
    watch(iw)%cpu%lap = time_cpu()     !CPU time of start of lap
    !watch(iw)%clock%lap = time_clock() !Wall clock time of start of lap
    watch(iw)%clock%lap = time_jul() !Wall clock time of start of lap
    watch(iw)%active  =  .true.
    watch(iw)%running =  .true.   
    
    return
    end subroutine watch_start

!********************************************************************************
    subroutine watch_stop(watch_name,watch_index)
! Stops a timing watch
!
! Usage:
!   call watch_stop('MyWatchName')
!   call watch_stop('MyWatchName',iwatch)
!   call watch_stop(watch_name='MyWatchName')
!   call watch_stop(watch_index=iwatch)
!
! written by Alex Sanchez, USACE-CHL    
!********************************************************************************
    use watch_def
    use diag_lib
    use time_lib
    implicit none
    !Input/Output
    character(len=*), intent(in), optional :: watch_name
    integer, intent(in), optional :: watch_index
    !Internal
    integer :: iw
    character(len=200) :: msg

    if(present(watch_index))then
      iw = watch_index
    elseif(present(watch_name))then
      iw = watch_name2index(watch_name)
    else
      call diag_print_error('Call to watch_stop requires a watch name or index')
    endif
    if(iw>nmaxwatches)then
      write(msg,*) 'Invalid watch index: ',iw,' >',nmaxwatches
      call diag_print_error(msg)
    endif
    if(.not.watch(iw)%active)then
      write(msg,*)  'Trying to stop inactive watch: ',iw
      call diag_print_error(msg)
    endif
    
    !CPU Time    
    watch(iw)%cpu%split = watch(iw)%cpu%split + time_cpu() - watch(iw)%cpu%lap 
    
    !Wall clock Time
    !watch(iw)%clock%split = watch(iw)%clock%split + time_clock() - watch(iw)%clock%lap 
    watch(iw)%clock%split = watch(iw)%clock%split + time_jul() - watch(iw)%clock%lap 
    
    !Stop running
    watch(iw)%running = .false.

    return
    end subroutine watch_stop
    
!***********************************************************
    function watch_name2index(watch_namein) result(iw)
! written by Alex Sanchez, USACE-CHL    
!***********************************************************    
    use watch_def
    use diag_lib
    implicit none
    !Input/Output
    character(len=*), intent(in) :: watch_namein    
    integer :: iw
    !Internal    
    character(len=watchcharlen) :: watch_name
    character(len=200) :: msg
    logical :: watchfound    
    
    !Takes care of character length differences
    watch_name = watch_namein
    
    !Find matching watch if it exist
    watchfound = .false.
    do iw=1,nwatches
      if(watch(iw)%name==watch_name)then        
        watchfound = .true. !already known watch
        return
      endif  
    enddo
    
    nwatches = nwatches + 1
    if(nwatches>nmaxwatches)then
      write(msg,*)  'Maximum size of watches',nmaxwatches,'exceeded'
      call diag_print_error(msg)
    endif
    iw = nwatches
    watch(iw)%name = watch_name
    
    return        
    end function watch_name2index    
    
!**********************************************************
    subroutine watch_output() 
! Outputs the watch information to the screen 
! and diagnostic file.
! written by Alex Sanchez, USACE-CHL        
!**********************************************************
    use watch_def
    use diag_lib
    implicit none
    integer :: i
    character(len=200) :: msg
    
444 format(1x,A,F10.4,6x,F10.4)
    
    call diag_print_message(' ') 
    call diag_print_message(' Watch Name                    Clock Time (sec)   CPU Time (sec)') 
    do i=1,nwatches
      write(msg,444) watch(i)%name,watch(i)%clock%split,watch(i)%cpu%split
      call diag_print_message(msg)
    enddo
    
    return
    end subroutine watch_output    
    
!**********************************************************
    subroutine watch_destroy()
! Destroys the watches    
! written by Alex Sanchez, USACE-CHL        
!**********************************************************
    use watch_def
    implicit none    
    
    if(allocated(watch)) deallocate(watch)

    nwatches = 0
    
    return
    end subroutine watch_destroy  
   
end module watch_lib
