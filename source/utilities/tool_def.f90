!***********************************************************************************************************************    
    module tool_def
      implicit none
      
      interface vstrlz                        !Overload the function so that it works for both real(4) and real(8) variables  MEB  03/29/2022
        module procedure vstrlz4, vstrlz8
      end interface

      interface rround                        !Overload the function so that it works for both real(4) and real(8) variables  MEB  04/05/2022
        module procedure rround4, rround8
      end interface
  
  contains
  
!************************************************************************ 
    pure function center(str,width) result(out) 
! Returns str centered within a field of 'width' characters. 
! Leading/trailing blanks in str are ignored; if str is longer than 
! width it is truncated.  Extra odd space goes to the right. 
!************************************************************************ 
    character(len=*), intent(in)  :: str 
    integer,          intent(in)  :: width 
    character(len=:), allocatable :: out 
    integer :: n, lpad 
 
    n = len_trim(adjustl(str)) 
    if (n >= width) then 
      out = adjustl(str) 
      out = out(1:width) 
    else 
      lpad = (width - n)/2 
      out = repeat(' ',lpad)//trim(adjustl(str)) 
    endif 
 
    end function center 
    
!************************************************************
      function vstrlz4(flt,gfmt) result(gbuf)
! This function will take a single precision float (4-bit) variable and a format declaration, then convert it to a string. 
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
      end function vstrlz4

!************************************************************
      function vstrlz8(flt,gfmt) result(gbuf)
! This function will take a double precision float (8-bit) variable and a format declaration, then convert it to a string. 
! The string will add a leading 0, +0, or -0 if necessary.  Also, it will remove a leading space from the PE format
      character(len=20)                      :: gbuf
      real(8), intent(in)                    :: flt
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
      end function vstrlz8
      
!****************************************************
! double precision rounding function - 12/06/07 meb
!   two inputs - value, and #digits to round off to
!   returns - rounded value to specified precision
!
!   This and single precision version - combined into 
!   overloaded function 'rround' - MEB  04/05/22
!****************************************************
    double precision function rround8 (X,P)    
    use prec_def
    implicit none
    real(8) X    !double precision value passed in
    integer K,P  !digits of precision
    real(8) PR,R
    
    PR=10.d0**P      
    K = NINT(PR*(X-AINT(X))) 
    R = AINT(X) + K/PR
    RROUND8 = R
    
    end function rround8
  
!****************************************************
! single precision rounding function - 12/06/07 meb
!   two inputs - value, and #digits to round off to
!   returns - rounded value to specified precision
!
!   This and double precision version - combined into 
!   overloaded function 'rround' - MEB  04/05/22
!****************************************************
    real function rround4 (X,P)
    implicit none
    real(4) X    !real value passed in
    integer K,P  !digits of precision
    real(4) PR,R
    
    PR=10.0**P      
    K = NINT(PR*(X-AINT(X))) 
    R = AINT(X) + K/PR
    RROUND4 = R
    
    end function rround4
    
        
!************************************************************
  pure function addQuotes(str) result(delimited)
! "Returns the delivered string surrounded by single quotes."
! Added by Mitchell Brown, 04/25/2022
!************************************************************
    character(*), intent(in) :: str
    character(len(str)+2) delimited

    delimited = "'"//str//"'"
  end function addQuotes
  
  
!************************************************************
  pure function adjustc(string,length)
! DESCRIPTION center text using implicit or explicit length
!************************************************************
    character(len=*),intent(in)  :: string         ! input string to trim and center
    integer,intent(in),optional  :: length         ! line length to center text in
    character(len=:),allocatable :: adjustc        ! output string
    integer                      :: inlen
    integer                      :: ileft          ! left edge of string if it is centered
    
    if(present(length))then                        ! optional length
      inlen=length                                 ! length will be requested length
      if(inlen.le.0)then                           ! bad input length
        inlen=len(string)                          ! could not use input value, fall back to length of input string
      endif
    else                                           ! output length was not explicitly specified, use input string length
      inlen=len(string)
    endif
    allocate(character(len=inlen):: adjustc)       ! create output at requested length
    adjustc(1:inlen)=' '                           ! initialize output string to all blanks

    ileft =(inlen-len_trim(adjustl(string)))/2     ! find starting point to start input string to center it
    if(ileft.gt.0)then                             ! if string will fit centered in output
      adjustc(ileft+1:inlen)=adjustl(string)       ! center the input text in the output string
    else                                           ! input string will not fit centered in output string
      adjustc(1:inlen)=adjustl(string)             ! copy as much of input to output as can
    endif
  end function adjustc  

!************************************************************
  pure logical function is_digit(c)
! Input:   single character
! Returns: .true. or .false.
!************************************************************
    character(len=1), intent(in) :: c
    is_digit = (c >= '0' .and. c <= '9')
  end function is_digit
  
!**************************************************************    
  subroutine check_percentile_file(filename, pathname)
! Checks each percentile file to make sure that there are no negative values.
! Added by M. Brown - 5/11/2026
!**************************************************************    
    use XMDF
    use diag_def,  only: msg, msg2
    use diag_lib,  only: diag_print_error
    use const_def, only: READWRITE
    implicit none
  
    character(len=*),intent(in) :: filename,pathname
  
    integer :: ERROR, NTIMES
    real(4) :: A_MIN(1)
    integer(XID) :: PID, DGID
  
    CALL XF_OPEN_FILE (TRIM(filename),READWRITE,PID,ERROR)
    IF (ERROR.LT.0) THEN
      write(msg,*) 'CANNOT OPEN FILE: ', TRIM(filename)
      call DIAG_PRINT_ERROR(msg)   
    ENDIF
  
    CALL XF_OPEN_GROUP (PID,pathname,DGID,ERROR)
    IF (ERROR.LT.0) THEN 
      call XF_CLOSE_FILE(PID, ERROR)
      write(msg,*) 'CANNOT OPEN DATASET: ', TRIM(pathname)
      call DIAG_PRINT_ERROR(msg)
    ENDIF
  
    CALL XF_GET_DATASET_MINS(DGID, 1, A_MIN, ERROR)
    if (A_MIN(1) < 0.0) then
      call XF_CLOSE_FILE(PID, ERROR)
      write(MSG,*) 'Minimum dataset value in dataset: ',TRIM(pathname),' is negative.'
      write(msg2,*) 'Sediment percentile datasets must contain only positive values.'
      call DIAG_PRINT_ERROR(msg, msg2)
    endif
  
  end subroutine
    
FUNCTION xmdf_error(code) RESULT(name) 
  use ERRORDEFINITIONS 
  implicit none 
   
  integer, intent(in) :: code 
  character(len=50)   :: name 
   
  select case(code) 
  ! File errors -40xx 
    case(ERROR_FILE_NOT_HDF5)             ; name = 'ERROR_FILE_NOT_HDF5 (-4001)' 
    case(ERROR_FILE_NOT_XMDF)             ; name = 'ERROR_FILE_NOT_XMDF (-4002)' 
  ! Attribute errors -41xx 
    case(ERROR_ATTRIBUTE_NOT_SUPPORTED)   ; name = 'ERROR_ATTRIBUTE_NOT_SUPPORTED (-4101)' 
  ! Datatype errors -42xx 
    case(ERROR_INCORRECT_DATATYPE)        ; name = 'ERROR_INCORRECT_DATATYPE (-4201)' 
  ! Dataset errors -43xx 
    case(ERROR_DATASET_SIZE_INCORRECT)    ; name = 'ERROR_DATASET_SIZE_INCORRECT (-4301)' 
    case(ERROR_DATASET_NO_DATA)           ; name = 'ERROR_DATASET_NO_DATA (-4302)' 
    case(ERROR_DATASET_DOES_NOT_EXIST)    ; name = 'ERROR_DATASET_DOES_NOT_EXIST (-4303)' 
    case(ERROR_DATASET_INVALID)           ; name = 'ERROR_DATASET_INVALID (-4304)' 
  ! Group errors -44xx 
    case(ERROR_GROUP_TYPE_INCONSISTENT)   ; name = 'ERROR_GROUP_TYPE_INCONSISTENT (-4401)' 
  ! Mesh errors -45xx 
    case(ERROR_ELEMENT_NUM_INCORRECT)     ; name = 'ERROR_ELEMENT_NUM_INCORRECT (-4501)'   !inconsistent element number 
    case(ERROR_NODE_NUM_INCORRECT)        ; name = 'ERROR_NODE_NUM_INCORRECT (-4502)' 
    case(ERROR_NOT_MESH_GROUP)            ; name = 'ERROR_NOT_MESH_GROUP (-4503)' 
    case(ERROR_MESH_INCOMPLETE)           ; name = 'ERROR_MESH_INCOMPLETE (-4504)' 
    case(ERROR_MESH_INVALID)              ; name = 'ERROR_MESH_INVALID (-4505)' 
  ! Grid errors 
    case(ERROR_GRID_TYPE_INVALID)         ; name = 'ERROR_GRID_TYPE_INVALID (-4601)' 
    case(ERROR_GRID_NUM_DIMS)             ; name = 'ERROR_GRID_NUM_DIMS (-4602)' 
    case(ERROR_GRID_EXTRUDE_TYPE_INVALID) ; name = 'ERROR_GRID_EXTRUDE_TYPE_INVALID (-4603)' 
    case(ERROR_GRID_NUMVALS_INCORRECT)    ; name = 'ERROR_GRID_NUMVALS_INCORRECT (-4604)' 
  ! Others 
    case(ERROR_OTHER)                     ; name = 'ERROR_OTHER (-9901)' 
  end select 
   
  end FUNCTION xmdf_error 

!************************************************************
    end module tool_def      