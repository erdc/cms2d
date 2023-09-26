!***************************************************************
    subroutine diag_default
!***************************************************************    
    use diag_def, only: msg, msg2, msg3, msg4, msg5, msg6, msg7, msg8, debug_mode
    implicit none
    
    !Debug model    
    debug_mode = .false.
    
    msg  = ''
    msg2 = ''
    msg3 = ''
    msg4 = ''
    msg5 = ''
    msg6 = ''
    msg7 = ''
    msg8 = ''
    
    return
    end subroutine diag_default
    
!*************************************************************   
    subroutine diag_cards(cardname,foundcard)
! Reads diagnositic variable
! written by Alex Sanchez, USACE-CHL
!*************************************************************
    use diag_def, only: debug_mode, dgfile, debug_level
    
    implicit none
    integer :: ierr
    character(len=37) :: cardname
    logical :: foundcard
    
    foundcard = .true.
    select case(cardname)
    case('DEBUG_MODE','DIAGNOSTIC_MODE')
      call card_boolean(77,debug_mode,ierr)
     
    case('DIAGNOSTIC_FILE')
      backspace(77)
      read(77,*,iostat=ierr) cardname,dgfile
      
    case('DIAGNOSTIC_LEVEL')
      backspace(77)
      read(77,*,iostat=ierr) cardname,debug_level
      
    case default
      foundcard = .false.        
      
    end select
    
    return
    end subroutine diag_cards
    
 !********************************************************************************
    subroutine diag_print
! Prints the general CMS settings to the screen and diagnostic file
! written by Alex Sanchez, USACE-CHL
!********************************************************************************  
    use diag_def, only: dgunit, dgfile, debug_mode
    
    implicit none
    integer :: i,iunit(2)
    
    open(dgunit,file=dgfile,access='append') 
    iunit = (/6,dgunit/)
    do i=1,2        
      if(debug_mode)then
        write(iunit(i),*)
        write(iunit(i),'(A,T40,A)') 'Debug Mode:','ON'
      endif
    enddo
    close(dgunit)
    
    return
    end subroutine diag_print