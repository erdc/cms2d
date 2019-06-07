!=================================================================
! CMS Dredging Routines
!
! Contains the following:
!   dredge_default - Sets default variables for dredge operations
!   dredge_cards   - Reads the dredging cards from the 
!                    control file
!   dredge_init    - Initializes the dredging variables
!   dredge_print   - Prints the salinity transport settings to the
!                    screen and diagnositic file
!   dredge_stat    - computes dredging statistics (not operational)
! written by Chris Reed
!=================================================================
   
!**************************************************************
    subroutine dredge_default
! Sets default variables for dredging simulation
! written by Chris Reed
!**************************************************************    
    use dredge_def
    implicit none
    
    dredging = .false.      
    ndredge_operations = 0
    dredge_interval = 5.0  !minutes
    write_dredge_diag = .false.
    
    return
    endsubroutine dredge_default

!**************************************************************
    subroutine dredge_cards(cardname,foundcard)
!   written by Chris Reed
!   minor modifications by Mitchell Brown 6/5/2019
!   - Added #ifdef statement, because dredging only works with XMDF files at present
!   - Added OUTPUT_DREDGE_DIAGNOSTICS card to only write extra information if desired.
!**************************************************************    
#include "CMS_cpp.h"
    use dredge_def
    use diag_lib, only: diag_print_error
    implicit none

    character cardname*37
	logical :: foundcard
    integer :: ierr
        
    foundcard = .true.
    selectcase(cardname)    
    case('DREDGE_UPDATE_INTERVAL')         
      backspace(77)
      read(77,*)cardname,dredge_interval

    case('DREDGE_OPERATION_BEGIN')
#ifdef WIN_OS        
      ndredge_operations = ndredge_operations + 1 
      dredging = .true.
      call dredge_operation_increment()         
      call dredge_op_block()
#else
      call diag_print_error('Dredge Operation not presently possible with ASCII-only input')
#endif
      
    case('OUTPUT_DREDGE_DIAGNOSTICS')
      call card_boolean(77,write_dredge_diag,ierr)        
    
    case default
      foundcard = .false.
!      write(*,*)trim(cardname)," not recognized"
       
    endselect
    
    return
    endsubroutine dredge_cards        
        
    
!********************************************************************    
    subroutine dredge_op_block()
! Reads an wse boundary condition block from the card file
! written by Chris Reed
!********************************************************************
    use dredge_def
    implicit none

    integer :: kk,ierr
    character cardname*37,operation_name*32
	logical :: foundcard

    do kk=1,20
      foundcard = .true.  
      read(77,*,iostat=ierr) cardname
      if(ierr/=0) exit
      if(cardname(1:1)=='!' .or. cardname(1:1)=='#') cycle      
      selectcase(cardname)
      case('NAME')    
        backspace(77)
        read(77,*)cardname,operation_name
        dredge_operations(ndredge_operations)%name = operation_Name  
      
      case('DREDGE_BEGIN')   
        call dredge_source_block()        
          
      case('PLACEMENT_BEGIN')
        call placement_block()        
        
      case('OUTPUT_DREDGE_DIAGNOSTICS')
        call card_boolean(77,write_dredge_diag,ierr)        
        
      case('DREDGE_OPERATION_END')
          return
          
      case default
        foundcard = .false.
        write(*,*)trim(cardname)," not recognized"          
      endselect
    enddo
    
    return
    endsubroutine dredge_op_block    
    

!********************************************************************    
    subroutine dredge_source_block()
! Reads an wse boundary condition block from the card file
! written by Chris Reed
!********************************************************************
    use dredge_def
    use diag_def
    use diag_lib
    use unitconv_lib, only: unitconv_var
    implicit none

    integer :: kk,ierr,m
    real(ikind) :: scalar,fac,con
    character*200 File,Path
	logical :: foundcard
    !Internal
    character :: cardname*37,aline*200,fromunits*10,tounits*10,defunits*10    

    do kk=1,20
      foundcard = .true.  
      read(77,*,iostat=ierr) cardname
      if(ierr/=0) exit
      if(cardname(1:1)=='!' .or. cardname(1:1)=='#') cycle      
      selectcase(cardname)
      case('DEPTH_DATASET')
        backspace(77)
        read(77,*)cardname,File,Path             
        dredge_operations(ndredge_operations)%DredgeSourceAreaFile = trim(File)
        dredge_operations(ndredge_operations)%DredgeSourceAreaPath = trim(Path)
!        write(msg2,*)'  Dredge file: ',trim(dredge_operations(ndredge_operations)%DredgeSourceAreaFile)
!        write(msg3,*)'  Dredge path: ',trim(dredge_operations(ndredge_operations)%DredgeSourceAreaPath)
!        call diag_print_message(msg2,msg3)
      
      case('START_METHOD')
        backspace(77)          
        read(77,*)cardname,METHOD 
        selectcase(trim(METHOD))
        case('CELL', 'SPECIFIED CELL')               !meb added 'SPECIFIED CELL' 6/5/2019 to be consistent between interface and CMS input
          dredge_operations(ndredge_operations)%Dredge_Approach = 2
        case('SHALLOW')
          dredge_operations(ndredge_operations)%Dredge_Approach = 1          
        case default
          write(msg,*)'Incorrect value entered for Start method ',TRIM(METHOD)
          call diag_print_error(msg)
        endselect
      
      case('START_CELL')
        backspace(77)          
        read(77,*)cardname, dredge_operations(ndredge_operations)%dredge_start_cell  
        !write(*,*)'for dredge oepartion number ',ndredge_operations
        !write(*,*)'start cell = ',dredge_operations(ndredge_operations)%dredge_start_cell
          
      case('DREDGE_RATE')
        backspace(77)          
        read(77,*)cardname, dredge_operations(ndredge_operations)%rate  
      
      case('TRIGGER_METHOD')
      !only used to check that proper trigger info is supplied
        backspace(77)          
        read(77,*)cardname,METHOD 
      
      case('TRIGGER_DEPTH')
        if(METHOD(1:5) == 'DEPTH' .or. METHOD(1:7) == 'PERCENT' )then
          call card_scalar(77,'m','m',scalar,ierr)
          dredge_operations(ndredge_operations)%Trigger_Depth  = scalar
          dredge_operations(ndredge_operations)%trigger_approach = 1
        else
          write(msg2,*)'Incorrect Card entered for Trigger method: ',TRIM(METHOD)
          call diag_print_error(msg2)
        endif 
        
      case('TRIGGER_VOLUME')
        if(METHOD(1:6) == 'VOLUME')then
          call card_scalar(77,'m^3','m^3',scalar,ierr)     
          dredge_operations(ndredge_operations)%trigger_approach = 2 
          dredge_operations(ndredge_operations)%Trigger_Vol = scalar     
        else
          write(msg2,*)'Incorrect Card entered for Trigger method: ',TRIM(METHOD)
          call diag_print_error(msg2)
        endif 
         
      case('TRIGGER_PERCENT')
        if(METHOD(1:7) == 'PERCENT') then        
          backspace(77)
          read(77,*)cardname,dredge_operations(ndredge_operations)%Trigger_percentage 
          dredge_operations(ndredge_operations)%trigger_approach = 3 
        else
          write(msg2,*)'Incorrect Card entered for Trigger method: ',TRIM(METHOD)
          call diag_print_error(msg2)
         endif 
         
      case('TRIGGER_TIME_PERIODS') 
        if(METHOD(1:12) == 'TIME_PERIODS') then
          backspace(77)
          read(77,*)cardname,num_trigger_TS
          if(allocated(trigger_start)) deallocate(trigger_start)
          if(allocated(trigger_finish)) deallocate(trigger_finish)     
          allocate(trigger_start(num_trigger_TS),trigger_finish(num_trigger_TS))
          dredge_operations(ndredge_operations)%trigger_approach = 4     
          !backspace(77)
          !read(77,*)cardname,num_trigger_TS,((trigger_start(m),trigger_finish(m)),m=1,num_trigger_TS) 
          backspace(77)
          read(77,'(A)') aline
          read(aline,*,iostat=ierr) cardname,num_trigger_TS
          dredge_operations(ndredge_operations)%num_trigger_intervals = num_trigger_TS      
          allocate(dredge_operations(ndredge_operations)%Trigger_start(num_trigger_TS))
          allocate(dredge_operations(ndredge_operations)%Trigger_finish(num_trigger_TS))        
          backspace(77)
          read(77,'(A)') aline 
          defunits = 'hr'
          tounits = 'hr'
          fromunits = ' '      

          read(aline,*,iostat=ierr) cardname, num_trigger_TS, (trigger_start(m),trigger_finish(m),m=1,num_trigger_TS), fromunits
          if(ierr==-1 .or. fromunits(1:1)==' ' .or. &
             fromunits(1:1)=='!' .or. fromunits(1:1)=='#')then
            fromunits = defunits !None specified in file, use default value
          endif 
          call unitconv_var(fromunits,tounits,fac,con) 
          do m=1,dredge_operations(ndredge_operations)%num_trigger_intervals
            dredge_operations(ndredge_operations)%Trigger_start(m) = fac*trigger_start(m)+con 
            dredge_operations(ndredge_operations)%Trigger_finish(m) = fac*trigger_finish(m)+con
          enddo 
        else
          write(msg2,*)'Incorrect Card entered for Trigger method: ',TRIM(METHOD)
          call diag_print_error(msg2)
        endif 
          
      case('DISTRIBUTION')
        backspace(77)          
        read(77,*)cardname,METHOD           
        selectcase(trim(METHOD))
        case('SEQUENTIAL')
          dredge_operations(ndredge_operations)%PAmethod = 1    
        case('PERCENT')
          dredge_operations(ndredge_operations)%PAmethod = 2       
        case default
          write(msg2,*)'DISTRIBUTION CARD INPUT NOT RECOGNIZED'
          call diag_print_error(msg2)
        endselect  
      
      case('DREDGE_END')
        return
       
      case default
        foundcard = .false.
        write(msg2,*) trim(cardname)," not recognized"
        call diag_print_warning(msg2)
        
      endselect
    enddo
    
    return
    endsubroutine dredge_source_block            
    
!********************************************************************    
    subroutine placement_block()
! Reads an wse boundary condition block from the card file
! written by Chris Reed
!********************************************************************
    use dredge_def
    implicit none

    integer :: kk,ierr
    real(ikind) :: scalar
    character*32 cardname
    character*200 file,path
	logical :: foundcard
    
    if (write_dredge_diag) then 
      open(unit=2056,file='dredge_module_diagnostics.txt',status='OLD')   !Moved default write statements to a diagnostic type file.  MEB 04/24/2017
      write(2056,*)'IN placement block num areas = ',num_place_areas
    endif
    
    num_place_areas = num_place_areas + 1
    dredge_operations(ndredge_operations)%NumPlacementAreas = num_place_areas
    call placement_area_increment()
    
    if (write_dredge_diag) then
      write(2056,*)'and updated to num areas = ',num_place_areas   
      write(2056,*)'num ar set to: ', dredge_operations(ndredge_operations)%NumPlacementAreas
    endif  
     
    do kk=1,20
      foundcard = .true.  
      read(77,*,iostat=ierr) cardname
      if(ierr/=0) exit
      if(cardname(1:1)=='!' .or. cardname(1:1)=='#') cycle      
      selectcase(cardname) 

      case('AREA','AREA_DATASET')
      backspace(77)
      read(77,*)cardname,File,Path             
      dredge_operations(ndredge_operations)%DredgePlaceAreaFile(num_place_areas) = trim(File)
      dredge_operations(ndredge_operations)%DredgePlaceAreaPath(num_place_areas) = trim(Path)
      
      case('START_METHOD','PLACEMENT_METHOD')   
        backspace(77)
        read(77,*)cardname,METHOD
        selectcase(trim(METHOD))
        case('UNIFORM')
          dredge_operations(ndredge_operations)%Placement_Approach(num_place_areas)=1
        
        case('CELL','SPECIFIED CELL')                !meb added 'SPECIFIED CELL' 6/5/2019 to be consistent between interface and CMS input
          dredge_operations(ndredge_operations)%Placement_Approach(num_place_areas)=2      
          
        case default
          write(*,*)'PLACEMENT_METHOD CARD INPUT NOT RECOGNIZED'
          STOP
        endselect  
       
      case('DISTRIBUTION_PERCENTAGE')         
        backspace(77)
        !write(*,*)"DP2",ndredge_operations,num_place_areas,dredge_operations(ndredge_operations)%DredgePlacementAllocation(num_place_areas)
        read(77,*)cardname,dredge_operations(ndredge_operations)%DredgePlacementAllocation(num_place_areas)
        !write(*,*)"DP3",ndredge_operations,num_place_areas,dredge_operations(ndredge_operations)%DredgePlacementAllocation(num_place_areas)
       
      case('START_CELL')         
        backspace(77)          
        read(77,*)cardname, dredge_operations(ndredge_operations)%Placement_Start_Cell(num_place_areas)  
      
      case('THICKNESS_LIMIT')
        call card_scalar(77,'m','m',scalar,ierr) 
        dredge_operations(ndredge_operations)%Placement_Thickness(num_place_areas) = scalar
      
      case('DEPTH_LIMIT')
        call card_scalar(77,'m','m',scalar,ierr) 
        dredge_operations(ndredge_operations)%Placement_Depth(num_place_areas) = scalar
      
      case('PLACEMENT_END')
        return
          
      case default
        foundcard = .false.
        write(*,*)trim(cardname)," not recognized"
        
      endselect
    enddo
    
    if (write_dredge_diag) close(2056)
    
    return
    endsubroutine placement_block   
    
!**************************************************************
    subroutine dredge_operation_increment()
!   subroutine checks to see if the operationName has already been defined, and
!   if not it will add a new operation name
!   it returns the OperationNum associated with the operationName, the OperationNum is tha array index * in the 
!   type Dredge_Operations(*)
!   written by Chris Reed
!**************************************************************    
    use dredge_def
    implicit none

    integer      :: i
    logical      :: found
    character*50 :: dredge_diag_file
    type(dredge_operations_type), allocatable :: dredge_operations_TMP(:)
    
    !First time opening this file
    dredge_diag_file='dredge_module_diagnostics.txt'
    
    INQUIRE(FILE=dredge_diag_file,EXIST=FOUND)
    IF (FOUND) THEN
      OPEN(2056,FILE=dredge_diag_file)
      CLOSE(2056,STATUS='DELETE')
    ENDIF

    if (write_dredge_diag) then
      open(unit=2056,file=dredge_diag_file,status='NEW')
      write(2056,*)
      write(2056,*) 'Number of Dredge Operations= ',ndredge_operations
    endif  
          
    if(ndredge_operations > 1) then  !then adding new operation to array already allocated            
      allocate(dredge_operations_TMP(ndredge_operations-1))
      do i=1,ndredge_operations-1
        allocate(dredge_operations_TMP(i)%DredgePlaceAreaFile(dredge_operations(i)%NumPlacementAreas))
        allocate(dredge_operations_TMP(i)%DredgePlaceAreaPath(dredge_operations(i)%NumPlacementAreas))          
        allocate(dredge_operations_TMP(i)%Placement_Approach(dredge_operations(i)%NumPlacementAreas)) 
        allocate(dredge_operations_TMP(i)%Placement_Start_Cell(dredge_operations(i)%NumPlacementAreas))          
        allocate(dredge_operations_TMP(i)%Placement_Limit(dredge_operations(i)%NumPlacementAreas,1)) 
        allocate(dredge_operations_TMP(i)%Placement_Thickness(dredge_operations(i)%NumPlacementAreas))          
        allocate(dredge_operations_TMP(i)%Placement_Depth(dredge_operations(i)%NumPlacementAreas))           
        allocate(dredge_operations_TMP(i)%DredgePlacementAllocation(dredge_operations(i)%NumPlacementAreas))          
        dredge_operations_TMP(i) = dredge_operations(i)
      enddo
          
      deallocate(dredge_operations)
      allocate(dredge_operations(ndredge_operations))
      do i=1,ndredge_operations-1
        dredge_operations(i) = dredge_operations_TMP(i)
      enddo          
    else            
      allocate(dredge_operations(ndredge_operations))         
    endif
            
    dredge_operations(ndredge_operations)%name = "no_name"
    dredge_operations(ndredge_operations)%active = .false.          
    dredge_operations(ndredge_operations)%NumPlacementAreas = 0
    allocate(dredge_operations(ndredge_operations)%DredgePlaceAreaFile(1))
    allocate(dredge_operations(ndredge_operations)%DredgePlaceAreaPath(1))                  
    dredge_operations(ndredge_operations)%Trigger_Depth = 0.0 
    dredge_operations(ndredge_operations)%Trigger_Vol = 0.0 
    dredge_operations(ndredge_operations)%Trigger_Percentage = 0.0            
    dredge_operations(ndredge_operations)%Trigger_Approach = 1 
!    dredge_operations(ndredge_operations)%Dredge_Depth(1) = 0.0   !m               !This isn't allocated now and is reallocated later.  Commenting out.  MEB 052516
    dredge_operations(ndredge_operations)%Rate = 0  !m3/day
    dredge_operations(ndredge_operations)%Dredge_Approach = 1
    dredge_operations(ndredge_operations)%dredge_start_cell = 1
    dredge_operations(ndredge_operations)%cell_inc = 1   
    dredge_operations(ndredge_operations)%PAmethod = 1           
    allocate(dredge_operations(ndredge_operations)%Placement_Approach(1))          
    dredge_operations(ndredge_operations)%Placement_Approach(1) = 1  
    allocate(dredge_operations(ndredge_operations)%Placement_Start_cell(1))          
    dredge_operations(ndredge_operations)%Placement_Start_cell(1) = 1                  
    allocate(dredge_operations(ndredge_operations)%Placement_Limit(1,1))          
    dredge_operations(ndredge_operations)%Placement_Limit(1,1) = -999.0 !m 
    allocate(dredge_operations(ndredge_operations)%Placement_Thickness(1))          
    dredge_operations(ndredge_operations)%Placement_Thickness(1) = -999 !m            
    allocate(dredge_operations(ndredge_operations)%Placement_Depth(1))          
    dredge_operations(ndredge_operations)%Placement_Depth(1) = -999 !m                   
    allocate(dredge_operations(ndredge_operations)%DredgePlacementAllocation(1))
    dredge_operations(ndredge_operations)%DredgePlacementAllocation(1)  = -999 !percent
    !write(*,*)"DP0",ndredge_operations,num_place_areas,dredge_operations(ndredge_operations)%DredgePlacementAllocation(1)           
    num_place_areas = 0
           
    if (write_dredge_diag) close(2056)
    
    return
    endsubroutine dredge_operation_increment    

!**************************************************************
    subroutine placement_area_increment()
!   subroutine checks to see if the operationName has already been defined, and
!   if not it will add a new operation name
!   it returns the OperationNum associated with the operationName, the OperationNum is tha array index * in the 
!   type Dredge_Operations(*)
!   written by Chris Reed
!**************************************************************    
    use dredge_def
    implicit none

    integer i
    character*200, allocatable :: DredgePlaceAreaFile(:)
    character*200, allocatable :: DredgePlaceAreaPath(:)          
    integer, allocatable :: Placement_Approach(:) 
    integer, allocatable :: Placement_Start_Cell(:)          
    !real, allocatable :: Placement_Limit(:) 
    real, allocatable :: Placement_Thickness(:)          
    real, allocatable :: Placement_Depth(:)           
    real, allocatable :: DredgePlacementAllocation(:) 
          
    if(num_place_areas > 1) then  !then adding new area to array already allocated
      allocate(DredgePlaceAreaFile(num_place_areas-1))
      allocate(DredgePlaceAreaPath(num_place_areas-1))          
      allocate(Placement_Approach(num_place_areas-1)) 
      allocate(Placement_Start_Cell(num_place_areas-1))          
      !allocate(Placement_Limit(num_place_areas-1)) 
      allocate(Placement_Thickness(num_place_areas-1))          
      allocate(Placement_Depth(num_place_areas-1))           
      allocate(DredgePlacementAllocation(num_place_areas-1))
          
      do i=1,num_place_areas-1
        DredgePlaceAreaFile(i)=dredge_operations(ndredge_operations)%DredgePlaceAreaFile(i)
        DredgePlaceAreaPath(i)=dredge_operations(ndredge_operations)%DredgePlaceAreaPath(i)          
        Placement_Approach(i)=dredge_operations(ndredge_operations)%Placement_Approach(i) 
        Placement_Start_Cell(i)=dredge_operations(ndredge_operations)%Placement_Start_Cell(i)          
        !Placement_Limit(i)=dredge_operations(ndredge_operations)%Placement_Limit(i) 
        Placement_Thickness(i)=dredge_operations(ndredge_operations)%Placement_Thickness(i)          
        Placement_Depth(i)=dredge_operations(ndredge_operations)%Placement_Depth(i)           
        DredgePlacementAllocation(i)=dredge_operations(ndredge_operations)%DredgePlacementAllocation(i)   
      enddo
         
      deallocate(dredge_operations(ndredge_operations)%DredgePlaceAreaFile)
      deallocate(dredge_operations(ndredge_operations)%DredgePlaceAreaPath)          
      deallocate(dredge_operations(ndredge_operations)%Placement_Approach) 
      deallocate(dredge_operations(ndredge_operations)%Placement_Start_Cell)          
      !deallocate(dredge_operations(ndredge_operations)%Placement_Limit) 
      deallocate(dredge_operations(ndredge_operations)%Placement_Thickness)          
      deallocate(dredge_operations(ndredge_operations)%Placement_Depth)           
      deallocate(dredge_operations(ndredge_operations)%DredgePlacementAllocation)  
       
      allocate(dredge_operations(ndredge_operations)%DredgePlaceAreaFile(num_place_areas))
      allocate(dredge_operations(ndredge_operations)%DredgePlaceAreaPath(num_place_areas))          
      allocate(dredge_operations(ndredge_operations)%Placement_Approach(num_place_areas)) 
      allocate(dredge_operations(ndredge_operations)%Placement_Start_Cell(num_place_areas))          
      !allocate(dredge_operations(ndredge_operations)%Placement_Limit(num_place_areas)) 
      allocate(dredge_operations(ndredge_operations)%Placement_Thickness(num_place_areas))          
      allocate(dredge_operations(ndredge_operations)%Placement_Depth(num_place_areas))           
      allocate(dredge_operations(ndredge_operations)%DredgePlacementAllocation(num_place_areas))          
        
      do i=1,num_place_areas-1
        dredge_operations(ndredge_operations)%DredgePlaceAreaFile(i)= DredgePlaceAreaFile(i)
        dredge_operations(ndredge_operations)%DredgePlaceAreaPath(i)= DredgePlaceAreaPath(i)      
        dredge_operations(ndredge_operations)%Placement_Approach(i)= Placement_Approach(i)
        dredge_operations(ndredge_operations)%Placement_Start_Cell(i)= Placement_Start_Cell(i)    
        !dredge_operations(ndredge_operations)%Placement_Limit(i)= Placement_Limit(i)
        dredge_operations(ndredge_operations)%Placement_Thickness(i)= Placement_Thickness(i)      
        dredge_operations(ndredge_operations)%Placement_Depth(i)= Placement_Depth(i)     
        dredge_operations(ndredge_operations)%DredgePlacementAllocation(i)=  DredgePlacementAllocation(i)
      enddo
         
      !allocate(dredge_operations(ndredge_operations)%DredgePlaceAreaFile(num_place_areas))
      !allocate(dredge_operations(ndredge_operations)%DredgePlaceAreaPath(num_place_areas))                            
      !allocate(dredge_operations(ndredge_operations)%Placement_Approach(num_place_areas))          
      dredge_operations(ndredge_operations)%Placement_Approach(num_place_areas) = 1  
      !allocate(dredge_operations(ndredge_operations)%Placement_Start_cell(num_place_areas))          
      dredge_operations(ndredge_operations)%Placement_Start_cell(num_place_areas) = 1                  
      !allocate(dredge_operations(ndredge_operations)%Placement_Limit(num_place_areas))          
      !dredge_operations(ndredge_operations)%Placement_Limit(1) = 99999.0 !m 
      !allocate(dredge_operations(ndredge_operations)%Placement_Thickness(num_place_areas))          
      dredge_operations(ndredge_operations)%Placement_Thickness(num_place_areas) = -999 !m            
      !allocate(dredge_operations(ndredge_operations)%Placement_Depth(num_place_areas))          
      dredge_operations(ndredge_operations)%Placement_Depth(num_place_areas) = -999 !m                   
      !allocate(dredge_operations(ndredge_operations)%DredgePlacementAllocation(num_place_areas))
      dredge_operations(ndredge_operations)%DredgePlacementAllocation(num_place_areas)  = -999  !percent
      !write(*,*)"DP1",ndredge_operations,num_place_areas,dredge_operations(ndredge_operations)%DredgePlacementAllocation(num_place_areas)           
    endif
           
    return
    endsubroutine placement_area_increment        
    
    
    
!**************************************************************
    subroutine dredge_print
!   subroutine writes summary of dredge operations setup
!   written by Chris Reed
!**************************************************************    
    use dredge_def
    use prec_def
    use diag_def, only: dgunit,dgfile
    implicit none

    integer i,j,k,m, iunit(2)
    real(ikind) :: sum

111 format(' ',A)
222 format(' ',A,1x,A)    
!342 format(' ',A,F5.2,A)
!353 format(' ',A,F6.3,A)
!163 format(' ',A,F8.3,A)
133 format(' ',A,F10.2)
134 format(' ',A,F12.4)
135 format(' ',A,I0,2(1x,F10.2))
345 format(' ',A,F8.5)
445 format(' ',A,1x,I0)    
787 format(' ',A,1x,A)
!788 format(' ',A,1x,A,A)
!799 format(' ',A,1pe9.2)
!845 format(' ',A,F7.3,A)    


    iunit = (/6,dgunit/)
    open(dgunit,file=dgfile,access='append')   !Moved dredge setup output the CMS-type diagnostic output file.  MEB 03/18/2019
    
    do i=1,2
      write(iunit(i),*) 
	  write(iunit(i),111)            'Dredge Module Setup'  
      write(iunit(i),445)            '  Number of dredge ops       ',ndredge_operations
      do j=1,ndredge_operations
        write(iunit(i),787)          '  Operation name:              ',trim(dredge_operations(j)%name)
        write(iunit(i),787)          '    Source area file:          ',trim(dredge_operations(j)%DredgeSourceAreaFile)
        write(iunit(i),787)          '    Source area path:          ',trim(dredge_operations(j)%DredgeSourceAreaPath) 
        write(iunit(i),445)          '    Number source area cells:  ',dredge_operations(j)%NumDredgeAreaCells
      
        sum=0.0
        do k=1,dredge_operations(j)%NumDredgeAreaCells
          sum=sum + dredge_operations(j)%dredge_depth(k)
        enddo
        sum=sum/dredge_operations(j)%NumDredgeAreaCells

        write(iunit(i),345)          '    Average dredge depth:      ',sum
        write(iunit(i),134)          '    Dredge rate (m^3/sec):     ',dredge_operations(j)%Rate
         
        if(dredge_operations(j)%dredge_approach == 2) then
          write(iunit(i),222)        '    Dredge Start method:       ','SPECIFIED CELL'          !Previously this was 'CELL'.  I changed to make terminology consistent between interface and CMS input.
          write(iunit(i),445)        '    Starting cell number:      ',dredge_operations(j)%dredge_start_cell
        elseif(dredge_operations(j)%dredge_approach == 1) then
          write(iunit(i),222)        '    Dredge Start method:       ','SHALLOW'                 !Previously this was 'UNIFORM', but that is not a Dredge Option, it is a Placement Option.
        endif
            
        if(dredge_operations(j)%Trigger_Approach == 1) then
          write(iunit(i),222)        '    Trigger method:            ','DEPTH'
          write(iunit(i),345)        '    Trigger depth (m):         ',dredge_operations(j)%Trigger_Depth
        elseif(dredge_operations(j)%Trigger_Approach == 2) then
          write(iunit(i),222)        '    Trigger method:            ','VOLUME'
          write(iunit(i),345)        '    Trigger volume (m^3):      ',dredge_operations(j)%Trigger_Vol
        elseif(dredge_operations(j)%Trigger_Approach == 2) then
          write(iunit(i),222)        '    Trigger method:            ','PERCENT'
          write(iunit(i),345)        '    Trigger percent (%):       ',dredge_operations(j)%Trigger_percentage  
        elseif(dredge_operations(j)%Trigger_Approach == 4) then
          write(iunit(i),222)        '    Trigger method:            ','TIME PERIODS'
          write(iunit(i),445)        '    Num. trigger time periods: ',dredge_operations(j)%num_trigger_intervals
          do m=1,dredge_operations(j)%num_trigger_intervals
            write(iunit(i),135)      '      Times (hrs):             ',m,dredge_operations(j)%Trigger_Start(m),dredge_operations(j)%Trigger_Finish(m)
          enddo
        endif
         
        if(dredge_operations(j)%PAMETHOD == 1) then    
          write(iunit(i),222)        '    Distribution method        ','SEQUENTIAL'
        elseif(dredge_operations(j)%PAMETHOD == 2) then           
          write(iunit(i),222)        '    Distribution method        ','PERCENT'
        endif
            
        write(iunit(i),445)          '    Number of placement areas: ',dredge_operations(j)%NumPlacementAreas 
        do k=1,dredge_operations(j)%NumPlacementAreas
          write(iunit(i),445)        '      Placement Area (PA)',k
          write(iunit(i),787)        '        Placement area file:   ',trim(dredge_operations(j)%DredgePlaceAreaFile(k))
          write(iunit(i),787)        '        Placement area path:   ',trim(dredge_operations(j)%DredgePlaceAreaPath(k))
          write(iunit(i),445)        '        Number PA cells:       ',dredge_operations(j)%NumPlacementAreaCells(k)
          
          if(dredge_operations(j)%Placement_Approach(k) == 1) then
            write(iunit(i),222)      '        Placement method:      ','UNIFORM'
          elseif(dredge_operations(j)%Placement_Approach(k) == 2) then
            write(iunit(i),222)      '        Placement method:      ','SPECIFIED CELL'
            write(iunit(i),445)      '        Starting cell:         ',dredge_operations(j)%Placement_start_cell(k) 
          endif        
          if(dredge_operations(j)%Placement_Depth(k) > -999.0) &
            write(iunit(i),133)      '        Depth Limit (m):       ',dredge_operations(j)%Placement_Depth(k)
          if(dredge_operations(j)%Placement_thickness(k) > -999.0) & 
            write(iunit(i),133)      '        Thickness Limit (m):   ',dredge_operations(j)%Placement_thickness(k)
          if(dredge_operations(j)%PAMETHOD == 2) & 
            write(iunit(i),133)      '        Distribution % (%)     ',dredge_operations(j)%DredgePlacementAllocation(k)
        enddo         
      enddo
    enddo  

    close(dgunit)
    
    return
    endsubroutine dredge_print  

!**************************************************************
    subroutine dredge_init
! Written By Chris Reed 
!**************************************************************       
#include "CMS_cpp.h"
    use size_def
    use geo_def,     only: idmap,x,y,dx,dy
    use geo_def,     only: zb,zb0
    !use comvarbl,    only:  ntsch
    use dredge_def
    use sed_def,     only: singlesize,nsed    
#ifdef XMDF_IO
    use xmdf
    use in_xmdf_lib, only: readscalh5
#endif
    use in_lib,      only: readscalTxt
    use diag_lib,    only: diag_print_error
    use diag_def,    only: msg
    implicit none

    integer :: error,dfile_id,dcell_id,sumcells,k,i,j,summax,kk,ii,jj,ncnt
    integer :: n,m,mm
    integer, allocatable :: cells(:)
    real(ikind) value_1,value_2
    real(ikind), allocatable :: TMPARRAY(:),cell_dist(:)
    character*200 :: datapath,datafile
    character*7   :: STipt 
    character*50  :: filedetails,dredge_diag_file
    character*300 :: line  
    character*1   :: ext1
    character*19, allocatable :: linePA(:)
    character(len=10) :: aext
    logical       :: found

    allocate (TMPARRAY(ncellsfull))
    
    !variables that holds the time series output values
    allocate(dredgeTS_Vars(ndredge_operations,10))
    allocate(DredgeUnit(ndredge_operations))   

!    open(unit=2056,file='dredge_module_diagnostics.txt',status='OLD',access='APPEND') !Moved default write statements to a diagnostic type file.  MEB 04/24/2017

    !
    !!First time opening this file
    !dredge_diag_file='dredge_module_diagnostics.txt'
    !INQUIRE(FILE=dredge_diag_file,EXIST=FOUND)
    !IF (FOUND) THEN
    !  OPEN(2056,FILE=dredge_diag_file)
    !  CLOSE(2056,STATUS='DELETE')
    !ENDIF
    !open(unit=2056,file=dredge_diag_file,status='NEW')

    do k=1,ndredge_operations
      if (write_dredge_diag) then  
        write(2056,*) ' '        
        write(2056,*)'k = ',k,'  in dredge init'
        write(2056,*)'name= ',adjustL(trim(dredge_operations(k)%DredgeSourceAreaFile))
        write(2056,*)'path= ',adjustL(trim(dredge_operations(k)%DredgeSourceAreaPath))
      endif
      
      !holdover from original code - it is allocated in card reader section and needs to be deallocated before being properly allocated in this section
      if (allocated(dredge_operations(k)%Placement_limit)) deallocate(dredge_operations(k)%Placement_limit)  
        
      ! get a list of source area cells and convert to active grid
      datafile = dredge_operations(k)%DredgeSourceAreaFile
      inquire(file=datafile,exist=found)
      if(.not.found)then
        msg='File not found: '//trim(datafile)
        call diag_print_error(msg)
      endif
      
      datapath = dredge_operations(k)%DredgeSourceAreaPath
      error = -999
      
      call fileext(trim(datafile),aext)      
      select case (aext)
      case('h5')
#ifdef XMDF_IO
      call XF_OPEN_FILE(datafile,READONLY,DFILE_ID,error)  
      if (write_dredge_diag) write(2056,*)'Internal dredge data file number = ',error
      
      call XF_OPEN_GROUP(DFILE_ID,trim(dataPATH),DCELL_ID,error)
      if (write_dredge_diag) write(2056,*)'Internal dredge data path number = ',error      
      
      call XF_READ_SCALAR_VALUES_TIMESTEP(DCELL_ID,1,ncellsfull,TMPARRAY,error)
      if (write_dredge_diag) write(2056,*)'Internal dredge data read error  = ',error      
#endif
      case('txt')
        call readscalTxt(datafile,TMPARRAY,error)
      end select
      
      sumcells = 0
      do i = 1,ncellsfull
        if(TMPARRAY(i) > 1.0e-20) then
          sumcells = sumcells + 1
        endif
      enddo
      dredge_operations(k)%NumDredgeAreaCells = sumcells
      allocate (dredge_operations(k)%DredgeAreaCells(sumcells),dredge_operations(k)%Dredge_depth(sumcells))
      dredge_operations(k)%dredge_depth = 0.0
      kk=0
      do i = 1,ncellsfull
        if(TMPARRAY(i) > 1.0e-20) then
          kk=kk+1
          if (write_dredge_diag) write(2056,*)k,i,kk,idmap(i)
          dredge_operations(k)%DredgeAreaCells(kk) = idmap(i) 
          dredge_operations(k)%dredge_depth(kk) = TMPARRAY(i)
        endif
      enddo  
#ifdef XMDF_IO
      call XF_CLOSE_FILE(DFILE_ID,error)
#endif      

      ! get list of cells in each placement area and convert to active grid
      SumMax = 0
      allocate(dredge_operations(k)%NumPlacementAreaCells(dredge_operations(k)%NumPlacementAreas))     
      DO j=1,dredge_operations(k)%NumPlacementAreas
        if (write_dredge_diag) then
          write(2056,*) ' '        
          write(2056,*)'j = ',j,'  in dredge init (placement areas)'
          write(2056,*)'name= ',adjustL(trim(dredge_operations(k)%DredgePlaceAreaFile(j)))
          write(2056,*)'path= ',adjustL(trim(dredge_operations(k)%DredgePlaceAreaPath(j)))
        endif  
         
        datafile = dredge_operations(k)%DredgePlaceAreaFile(j)        
        inquire(file=datafile,exist=found)
        if(.not.found)then
          msg='File not found: '//trim(datafile)
          call diag_print_error(msg)
        endif
        datapath = dredge_operations(k)%DredgePlaceAreaPath(j)

#ifdef XMDF_IO
        call XF_OPEN_FILE(datafile,READONLY,DFILE_ID,error) 
        if (write_dredge_diag) write(2056,*)'Internal placement data file number = ',error        
        
        call XF_OPEN_GROUP(DFILE_ID,trim(dataPATH),DCELL_ID,error)
        if (write_dredge_diag) write(2056,*)'Internal placement data path number = ',error         
        
        call XF_READ_SCALAR_VALUES_TIMESTEP(DCELL_ID,1,ncellsfull,TMPARRAY,error)
        if (write_dredge_diag) write(2056,*)'Internal placement data read error  = ',error          
#endif
        
        sumcells = 0
        do i = 1,ncellsfull
          if(TMPARRAY(i) > 1.0e-20) sumcells = sumcells + 1
        enddo
        dredge_operations(k)%NumPlacementAreaCells(j) = sumcells
        if (write_dredge_diag) write(2056,*)'PA # of cells: ',k,j,dredge_operations(k)%NumPlacementAreaCells(j)
        SumMax = max(SumMax,sumcells)
#ifdef XMDF_IO
        call XF_CLOSE_FILE(DFILE_ID,error)  
#endif
        if (write_dredge_diag) write(2056,*)'sumMax = ',sumMax
      ENDDO
      allocate(dredge_operations(k)%PlacementAreaCells(dredge_operations(k)%NumPlacementAreas,SumMax)) 
      allocate(dredge_operations(k)%Placement_Limit(dredge_operations(k)%NumPlacementAreas,SumMax))        
      allocate(dredge_operations(k)%Placement_Area(dredge_operations(k)%NumPlacementAreas))
      DO j=1,dredge_operations(k)%NumPlacementAreas 
        datafile = dredge_operations(k)%DredgePlaceAreaFile(j)        
        datapath = dredge_operations(k)%DredgePlaceAreaPath(j)
#ifdef XMDF_IO
        call XF_OPEN_FILE(datafile,READONLY,DFILE_ID,error)         
        call XF_OPEN_GROUP(DFILE_ID,trim(dataPATH),DCELL_ID,error)
        call XF_READ_SCALAR_VALUES_TIMESTEP(DCELL_ID,1,ncellsfull,TMPARRAY,error)   
#endif        
        kk=0
        do i = 1,ncellsfull
          if(TMPARRAY(i) > 1.0e-20) then
            kk=kk+1
            if (write_dredge_diag) write(2056,*)j,i,kk,idmap(i)
            dredge_operations(k)%PlacementAreaCells(j,kk) = idmap(i)  
          endif
        enddo 
        !calculate footprint area for placement area
        dredge_operations(k)%Placement_Area(j)=0.
        if (write_dredge_diag) write(2056,*)'k,j= ',k,j
        if (write_dredge_diag) write(2056,*)'PA # of cells: ',dredge_operations(k)%NumPlacementAreaCells(j)
        do i=1,dredge_operations(k)%NumPlacementAreaCells(j)
          if (write_dredge_diag) write(2056,*)'j,i = ',j,i
          ii=dredge_operations(k)%PlacementAreaCells(j,i)
          dredge_operations(k)%Placement_Area(j)=dredge_operations(k)%Placement_Area(j)+dx(ii)*dy(ii)
        enddo
#ifdef XMDF_IO
        call XF_CLOSE_FILE(DFILE_ID,error)      
#endif
      ENDDO 
    
      !calculate placement limit for each placment area
      DO j=1,dredge_operations(k)%NumPlacementAreas
        value_1 = dredge_operations(k)%placement_depth(j)
        if(value_1 .gt. -999) then
          do i=1,dredge_operations(k)%NumPlacementAreaCells(j)
            ii=dredge_operations(k)%PlacementAreaCells(j,i)    
            value_2 = -zb(ii) - dredge_operations(k)%placement_thickness(j)
            if(value_2 .lt. 999) then
              dredge_operations(k)%placement_limit(j,i)=min(value_1,value_2)
            else
              dredge_operations(k)%placement_limit(j,i)=value_1
            endif
            if (write_dredge_diag) write(2056,*)k,j,i,dredge_operations(k)%placement_limit(j,i),dredge_operations(k)%placement_thickness(j),zb(ii)
          enddo
        else
          do i=1,dredge_operations(k)%NumPlacementAreaCells(j)
            ii=dredge_operations(k)%PlacementAreaCells(j,i)    
            dredge_operations(k)%placement_limit(j,i)=-zb(ii) - dredge_operations(k)%placement_thickness(j)
            if (write_dredge_diag) write(2056,*)k,j,i,dredge_operations(k)%placement_limit(j,i),dredge_operations(k)%placement_thickness(j),zb(ii)
          enddo            
        endif
      enddo
      !pause   !removed by MEB 2/23/2017
        
      !for tracking volumes
      dredge_operations(k)%total_dredged_vol = 0.
      dredge_operations(k)%total_placed_vol = 0.
      allocate(dredge_operations(k)%total_placed_by_area(dredge_operations(k)%NumPlacementAreas))
      allocate(dredge_operations(k)%placed_vol_by_area(dredge_operations(k)%NumPlacementAreas))
      dredge_operations(k)%total_placed_by_area = 0.
      dredge_operations(k)%placed_vol_by_area = 0.
    
      !processing for dredge approach 2  (sort cells from closest to farthest from start cell)
      if(dredge_operations(k)%dredge_approach == 2) then  !sort cells by proximity to start cell
        if (write_dredge_diag) write(2056,*)"srt cell ",k, dredge_operations(k)%dredge_start_cell
        dredge_operations(k)%dredge_start_cell = idmap(dredge_operations(k)%dredge_start_cell)
        ii=dredge_operations(k)%dredge_start_cell 
        ncnt=dredge_operations(k)%NumDredgeAreaCells
        allocate(cells(ncnt),cell_dist(ncnt))
        do i=1,ncnt
          jj=dredge_operations(k)%DredgeAreaCells(i)
          cell_dist(i) = sqrt((x(ii)-x(jj))**2+(y(ii)-y(jj))**2)
          cells(i) = jj
        enddo 
        call sort2(cell_dist,cells,ncnt)  
        dredge_operations(k)%DredgeAreaCells = cells
      endif
 
      !processing for placement approach 2  (sort cells from closest to farthest from start cell)
      allocate(dredge_operations(k)%placement_cell_inc(dredge_operations(k)%NumPlacementAreas))
      DO j=1,dredge_operations(k)%NumPlacementAreas
        if(dredge_operations(k)%placement_approach(j) == 2) then  !sort cells by proximity to start cell
          dredge_operations(k)%placement_start_cell(j) = idmap(dredge_operations(k)%placement_start_cell(j))
          ii=dredge_operations(k)%placement_start_cell(j) 
          ncnt=dredge_operations(k)%NumPlacementAreaCells(j)
          if(allocated(cells)) deallocate(cells)
          if(allocated(cell_dist)) deallocate(cell_dist)
          allocate(cells(ncnt),cell_dist(ncnt))
          do i=1,ncnt
            jj=dredge_operations(k)%PlacementAreaCells(j,i)
            cell_dist(i) = sqrt((x(ii)-x(jj))**2+(y(ii)-y(jj))**2)
            cells(i) = jj
          enddo 
          call sort2(cell_dist,cells,ncnt)  
          do i=1,ncnt   
            dredge_operations(k)%PlacementAreaCells(j,i) = cells(i)
            !write(*,*)k,j,cells(i),dredge_operations(k)%PlacementAreaCells(j,i)
          enddo
          dredge_operations(k)%placement_cell_inc(j) = 1
        endif 
      ENDDO
        
      !processing for trigger approach 3  (convert trigger_percentage form percent to fraction)
      !if(dredge_operations(k)%dredge_approach == 3) dredge_operations(k)%trigger_percentage = dredge_operations(k)%trigger_percentage/100.0
      !convert dredge rate for each operation from m3/day to m3/sec
      dredge_operations(k)%rate =  dredge_operations(k)%rate / (3600*24)
      
      !open time series output files for dredge tracking
      write(STipt,"(i7)") k
      N = index(STipt," ",back=.true.)       

      if (write_dredge_diag) then
        filedetails = 'DredgeDetailsOpNum_'//trim(STipt(N+1:7))//'.csv'
        write(*,*)'Writing data into file: ',trim(filedetails)
        DredgeUnit(k) = 850+k
        open(unit=DredgeUnit(k),file=filedetails)
      
        !line = 'time, DVol_Tem_Tot, DVol_Planned, DVol_Dredged, Vol_Placed, Vol_Offgrd, Vol_Offgrd_Tot, Vol_dredged_Tot, Vol_placed_Tot, Dredge_vol_all'
        line = 'Time, Dredged Vol., Placed Vol., Excess Vol., Tot. Dredge Vol., Total Placed Vol., Total Excess Vol., Dredge Vol. Remaining, Vol. in Source Area'
        line=adjustL(line)
        if(dredge_operations(k)%NumPlacementAreas > 1) then            
          !write(*,*) 'in op multi line ',k
          mm=dredge_operations(k)%NumPlacementAreas
          allocate(linePA(mm))
          do m=1,mm
            write(ext1,"(i1)")m
            linePA(m) = ",VolPA_"//ext1//",TotVolPA_"//ext1   
          enddo
          write(DredgeUnit(k),999) trim(line), (linePA(m),m=1,mm)  !MEB 04/24/2017 Switched from using the FORMAT backslash character to concatenate for issues with GNU compiler.
        else            
          write(DredgeUnit(k),"(A)")trim(line)
        endif
      endif  
    enddo  !end of dredge operations
999 FORMAT (A,1000(',',a19))
    
    
    !convert dredge time interval from minutes to seconds
    dredge_interval = dredge_interval * 60
    
    !in case sediment transport has not been activated then:
    if(.not. allocated(zb0) ) then
      allocate(zb0(ncellsD))
      zb0 = zb
    endif
    
    
    !find max array dimension needed for bed_change and allocate
    summax = 0
    do k=1,ndredge_operations
      summax = max(summax,dredge_operations(k)%NumDredgeAreaCells)
      DO j=1,dredge_operations(k)%NumPlacementAreas
        summax = max(summax,dredge_operations(k)%NumPlacementAreaCells(j))
      enddo
    enddo
    allocate(bed_change(summax))
    if(.not. singlesize) allocate(dredge_mat_gradation(nsed))
    
    if (write_dredge_diag) close(2056)
    
    return
    endsubroutine dredge_init
