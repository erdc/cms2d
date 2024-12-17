!*****************************************************************
    subroutine tg_alloc
! Resizes the tidegate structure variables
!*****************************************************************
    use struct_def
    implicit none

    integer :: itg
    type(TG_type), allocatable :: TG_temp(:)
    
    itidegate = itidegate + 1
    if (itidegate == 1) then
      allocate(TG_struct(1))
    else
      allocate(TG_temp(itidegate-1))
      do itg=1,itidegate-1
        TG_temp(itg) = TG_struct(itg)
      enddo
      deallocate(TG_struct)
      allocate(TG_struct(itidegate))
      do itg=1,itidegate-1
        TG_struct(itg) = TG_temp(itg)
      enddo
      deallocate(TG_temp)
    endif      
    
    !Initialize and set default values
    TG_struct(itidegate)%coefTGlateral = 0.0
    TG_struct(itidegate)%orientTGsea = 1
    TG_struct(itidegate)%coefTG_b2s = 0.0
    TG_struct(itidegate)%coefTG_s2b = 0.0
    TG_struct(itidegate)%open_height = 0.0
    TG_struct(itidegate)%bottom_elev = 0.0
    TG_struct(itidegate)%method = 1
    TG_struct(itidegate)%mtidegateopt = 1
    TG_struct(itidegate)%reg_start_time = 0
    TG_struct(itidegate)%reg_open_freq = 0
    TG_struct(itidegate)%reg_open_dur = 0.0
    TG_struct(itidegate)%ntidegateopt = 0
    
    return
    end subroutine tg_alloc

!************************************************************************************
    subroutine new_tg_block()
! Reads the multiple block structure for tidegates (meb, 10/23/2024)
!   Based on new_culvert_block
!************************************************************************************
    use const_def, only: deg2rad
    use struct_def
    use diag_lib
    implicit none
    
    integer :: i,itg,ierr,maxtidegate
    integer :: ival,numids
    real    :: rval
    character(len=34) :: cardname
    character(len=5)  :: cdum
    logical :: foundcard

    foundcard = .true.
    
    call tg_alloc  !increment and initialize multiple weir structure
    
    do
      read(77,*,iostat=ierr) cardname
      select case(cardname)
        case('TIDEGATE_STRUCT_END','END')
          exit
          
        case('CELL_IDS')
          backspace(77)
          read(77,*) cardname,numids    !No. of IDs
          TG_struct(itidegate)%ncells = numids
          allocate(TG_struct(itidegate)%cells(numids))
          backspace(77)
          read(77,*) cardname,numids,(TG_struct(itidegate)%cells(itg),itg=1,numids)
          
        case('DISTRIBUTION_COEFFICIENT')
          backspace(77)
          read(77,*) cardname, TG_struct(itidegate)%coefTGlateral

        case('ORIENTATION_SEA')
          backspace(77)
          read(77,*) cardname, cdum
          selectcase(cdum(1:4))
          case('NORT')
            TG_struct(itidegate)%orientTGsea = 1
          case('EAST')
            TG_struct(itidegate)%orientTGsea = 2
          case('SOUT')
            TG_struct(itidegate)%orientTGsea = 3
          case('WEST')
            TG_struct(itidegate)%orientTGsea = 4
          case default
            call diag_print_error('Orientation value must be one of the four main cardinal directions')
          end select
            
        case('FLOW_COEFFICIENT_FROM_BAY','FLOW_COEFF_FROM_BAY')
          backspace(77)
          read(77,*) cardname, rval
          TG_struct(itidegate)%coefTG_b2s = rval
          
        case('FLOW_COEFFICIENT_FROM_SEA','FLOW_COEFF_FROM_SEA')
          backspace(77)
          read(77,*) cardname, rval
          TG_struct(itidegate)%coefTG_s2b = rval
          
        case('OPENING_HEIGHT')
          backspace(77)
          read(77,*) cardname, rval
          TG_struct(itidegate)%open_height = rval

        case('BOTTOM_ELEVATION')
          backspace(77)
          read(77,*) cardname, rval
          TG_struct(itidegate)%bottom_elev = rval

        case('METH','METHOD')
          backspace(77)
          read(77,*) cardname, ival 
          if (ival .lt. 1 .or. ival .gt. 2) then
            call diag_print_error('Method must be either 1 or 2')
          endif 
          TG_struct(itidegate)%method = ival
          
        case('SCHED_OPERATION_TYPE')
          backspace(77)
          read(77,*) cardname, cdum
          selectcase(cdum(1:3))
            case ('REG')
              TG_struct(itidegate)%mtidegateopt = 1
              TG_struct(itidegate)%ntidegateopt = 3
            case ('DES')
              TG_struct(itidegate)%mtidegateopt = 2
              !will get 'ntidegateopt' below
            case ('EBB')
              TG_struct(itidegate)%mtidegateopt = 3
              !no options added for this type
            case ('UCG')
              TG_struct(itidegate)%mtidegateopt = 4
              !no options added for this type
            case default
              call diag_print_error('Unknown operation schedule type')
          endselect
            
        case('REG_START_TIME')
          backspace(77)
          read(77,*) cardname, TG_struct(itidegate)%reg_start_time
          
        case('REG_OPEN_FREQUENCY')
          backspace(77)
          read(77,*) cardname, TG_struct(itidegate)%reg_open_freq
          
        case('REG_OPEN_DURATION')
          backspace(77)
          read(77,*) cardname, TG_struct(itidegate)%reg_open_dur
          
        case('NUM_DESIGNATED_TIMES')
          backspace(77)
          read(77,*) cardname, ival
          TG_struct(itidegate)%ntidegateopt = ival*2
          allocate(TG_struct(itidegate)%des_start_time(ival))
          allocate(TG_struct(itidegate)%des_open_dur(ival))
          
        case('DES_START_TIME')
          backspace(77)
          ival = int(TG_struct(itidegate)%ntidegateopt/2)
          read(77,*) cardname,(TG_struct(itidegate)%des_start_time(itg),itg=1,ival)
          
        case('DES_OPEN_DURATION')
          backspace(77)
          ival = int(TG_struct(itidegate)%ntidegateopt/2)
          read(77,*) cardname,(TG_struct(itidegate)%des_open_dur(itg),itg=1,ival)
          
        case default
          write(*,*) 'WARNING: Card ',cardname,' not found'
          foundcard = .false.
         
        end select
    enddo
    
    return
    end subroutine new_tg_block
  
!****************************************************    
    subroutine tg_init
!****************************************************
#include "CMS_cpp.h"
    use struct_def
    use const_def, only: deg2rad    
    implicit none
    
    integer itg, i, j, k, ival, maxtidegate, maxtidegateopt, cell_inc, opts_inc
    
    numtidegate = itidegate
    allocate(ntidegate(0:numtidegate),orienttidegate(numtidegate),methtidegate(numtidegate),  &
             coeftidegate(numtidegate,2),elevtidegate(numtidegate),orienttgbay(numtidegate),  &
             orienttgsea(numtidegate),openhgttidegate(numtidegate),Qtottidegate(numtidegate))
    allocate(mtidegateopt(numtidegate),ntidegateopt(0:numtidegate),aschedtype(numtidegate),   &
             aorienttgsea(numtidegate))

    ntidegate(0) = 0
    ntidegateopt(0) = 0
    maxtidegate = 0      !Maximum number of cells associated with tide gates
    maxtidegateopt = 0   !Maximum number of option times associated with tide gates
    do i=1,numtidegate
      maxtidegate     = maxtidegate + TG_struct(i)%ncells
      maxtidegateopt  = maxtidegateopt + TG_struct(i)%ntidegateopt
      ntidegate(i)    = maxtidegate
      ntidegateopt(i) = maxtidegateopt 
    enddo
    allocate(idtidegate(maxtidegate),opentidegate(maxtidegate),qtidegate(maxtidegate), &
             coeftglateral(maxtidegate),dqtgdzdown(maxtidegate),dqtgdzup(maxtidegate))
    allocate(tidegateopt(maxtidegateopt))
    dqtgdzdown = 0.0;  dqtgdzup = 0.0  
     
    cell_inc = 0
    opts_inc = 0
    do itg = 1,numtidegate
      do j = 1,TG_struct(itg)%ncells
        cell_inc=cell_inc+1
        idtidegate(cell_inc)    = TG_struct(itg)%cells(j)
        coeftglateral(cell_inc) = TG_struct(itg)%coefTGLateral
      enddo
      orienttidegate(itg)  = TG_struct(itg)%orientTGsea
      coeftidegate(itg,1)  = TG_struct(itg)%coefTG_b2s
      coeftidegate(itg,2)  = TG_struct(itg)%coefTG_s2b
      openhgttidegate(itg) = TG_struct(itg)%open_height
      elevtidegate(itg)    = TG_struct(itg)%bottom_elev
      methtidegate(itg)    = TG_struct(itg)%method
      orienttgsea(itg)     = TG_struct(itg)%orientTGsea   !Define direction of sea side at local coordinate

      if(orienttidegate(itg).eq.1) then
        orienttgbay(itg) = 3
        aorienttgsea(itg) = 'North'
      elseif(orienttidegate(itg).eq.2) then
        orienttgbay(itg) = 4
        aorienttgsea(itg) = 'East'
      elseif(orienttidegate(itg).eq.3) then
        orienttgbay(itg) = 1
        aorienttgsea(itg) = 'South'
      elseif(orienttidegate(itg).eq.4) then
        orienttgbay(itg) = 2
        aorienttgsea(itg) = 'West'
      endif
      ival = TG_struct(itg)%mtidegateopt
      mtidegateopt(itg) = ival
      if (ival .eq. 1) then
        tidegateopt(opts_inc+1) = TG_struct(itg)%reg_start_time
        tidegateopt(opts_inc+2) = TG_struct(itg)%reg_open_freq
        tidegateopt(opts_inc+3) = TG_struct(itg)%reg_open_dur
        aschedtype(itg) = 'Regular'
        opts_inc = opts_inc + 3
      elseif (ival .eq. 2) then
        do k=1,int(TG_struct(itg)%ntidegateopt/2)
          tidegateopt(opts_inc+1) = TG_struct(itg)%des_start_time(k)
          tidegateopt(opts_inc+2) = TG_struct(itg)%des_open_dur(k)
          opts_inc = opts_inc + 2
        enddo
        aschedtype(itg) = 'Designated'
      elseif (ival .eq. 3) then
        aschedtype(itg) = 'Open on Ebb'
      elseif (ival .eq. 4) then
        aschedtype(itg) = 'Uncontrolled'
      endif
    enddo
      
    return
    end subroutine tg_init
         
!************************************************************************************
    subroutine tide_gate_block
! Reads tidal gate structures (hli, 02/28/13) 
!************************************************************************************
    use const_def, only: deg2rad
    use struct_def
    use diag_lib, only: diag_print_error                                                                    !MEB 06/21  better descriptive output without application error
    implicit none
    integer :: i,ii,ierr,maxtidegate
    character(len=34) :: cardname
    logical :: foundcard

    foundcard = .true.
    do
      read(77,*,iostat=ierr) cardname
      select case(cardname)
        case('TIDE_GATE_END','END')
          exit
          
        case('NUMBER_TIDE_GATE')
          backspace(77)
          read(77,*) cardname,numtidegate
          allocate(ntidegate(0:numtidegate),orienttidegate(numtidegate),methtidegate(numtidegate),  &
                   coeftidegate(numtidegate,2),elevtidegate(numtidegate),orienttgbay(numtidegate),  &
                   orienttgsea(numtidegate),openhgttidegate(numtidegate),Qtottidegate(numtidegate) )

        case('NUM_CELL_TIDE_GATE')
          backspace(77)
          read(77,*,iostat=ierr) cardname,(ntidegate(ii),ii=1,numtidegate)              
          if(ierr/=0) call diag_print_error('Must specify the number of cells for each tide gate')           !MEB 06/21  better descriptive output without application error
          
          ntidegate(0) = 0
          maxtidegate = 0
          do ii=1,numtidegate
            maxtidegate = maxtidegate + ntidegate(ii)
          enddo
          do ii=1,numtidegate
            ntidegate(ii) = ntidegate(ii-1) + ntidegate(ii)
          enddo

          allocate(idtidegate(maxtidegate),opentidegate(maxtidegate),qtidegate(maxtidegate), &
                   coeftglateral(maxtidegate),dqtgdzdown(maxtidegate),dqtgdzup(maxtidegate))

        case('CELLS')
          backspace(77)
          read(77,*,iostat=ierr) cardname,(idtidegate(ii),ii=1,maxtidegate)     
          if(ierr/=0) call diag_print_error('Must specify the Cell IDs for each tide gate')
          
        case('DISTRIBUTION_COEFFICIENT')
          backspace(77)
          read(77,*,iostat=ierr) cardname,(coeftglateral(ii),ii=1,maxtidegate)     
          if(ierr/=0) call diag_print_error('Must specify the distribution coefficient for each tide gate')          

        case('ORIENTATION')
          backspace(77)
          read(77,*,iostat=ierr) cardname,(orienttidegate(ii),ii=1,numtidegate)     
          if(ierr/=0) call diag_print_error('Must specify the orientation for each tide gate')
          
        case('FLOW_COEFFICIENT')
          backspace(77)
          read(77,*,iostat=ierr) cardname,(coeftidegate(ii,1),coeftidegate(ii,2),ii=1,numtidegate) 
          if(ierr/=0) call diag_print_error('Must specify the flow coefficient for each tide gate')
          
        case('OPEN_HEIGHT')
          backspace(77)
          read(77,*,iostat=ierr) cardname,(openhgttidegate(ii),ii=1,numtidegate) 
          if(ierr/=0) call diag_print_error('Must specify the gate opening height for each tide gate')
          
        case('BOTTOM_ELEVATION')
          backspace(77)
          read(77,*,iostat=ierr) cardname,(elevtidegate(ii),ii=1,numtidegate)  
          if(ierr/=0) call diag_print_error('Must specify the bottom elevation for each tide gate')

        case('METH')
          backspace(77)
          read(77,*,iostat=ierr) cardname,(methtidegate(ii),ii=1,numtidegate)     
          if(ierr/=0) call diag_print_error('Must specify the flux calculation method for each tide gate')

        !Initialize       
        dqtgdzdown = 0.0;  dqtgdzup = 0.0  
     
        do ii=1,numtidegate
          orienttgsea(ii) = orienttidegate(ii)   !Define direction of sea side at local coordinate
          if(orienttidegate(ii).eq.1) then
            orienttgbay(ii) = 3
          elseif(orienttidegate(ii).eq.2) then
            orienttgbay(ii) = 4
          elseif(orienttidegate(ii).eq.3) then
            orienttgbay(ii) = 1
          elseif(orienttidegate(ii).eq.4) then
            orienttgbay(ii) = 2
          endif
        enddo

      case('SCHEDULE_BEGIN')
        call tg_schedule_block

      case default
        write(*,*) 'WARNING: Card ',cardname,' not found'
        foundcard = .false.
         
    end select
    enddo
    
    return
    end subroutine tide_gate_block

!************************************************************************************
    subroutine tg_schedule_block
! Reads tidal gate schedule (hli, 02/28/13) 
!************************************************************************************
    use const_def, only: deg2rad
    use struct_def
    use diag_lib, only: diag_print_error
    implicit none
    integer :: i,ii,ierr,k,maxtidegateopt
    character(len=34) :: cardname
    character(len=100), allocatable :: operation_type(:)
    logical :: foundcard

    foundcard = .true.
    
    allocate(mtidegateopt(numtidegate),ntidegateopt(0:numtidegate),operation_type(numtidegate))

    do
      read(77,*,iostat=ierr) cardname
      select case(cardname)
        case('SCHEDULE_END','END')
          exit
        
        case('OPERATION_TYPE')
          backspace(77)
          read(77,*,iostat=ierr) cardname,(operation_type(ii),ii=1,numtidegate)
          if(ierr/=0) call diag_print_error('Must specify the operation schedule type for each tide gate')
          do ii=1,numtidegate
            if(operation_type(ii)(1:3).eq.'REG') then
              mtidegateopt(ii)=1
            elseif(operation_type(ii)(1:3).eq.'DES') then
              mtidegateopt(ii)=2
            elseif(operation_type(ii)(1:3).eq.'EBB') then
              mtidegateopt(ii)=3
            else
              mtidegateopt(ii)=4
            endif
          enddo

        case('NUM_CONTROL_ELEMENT')
          backspace(77)
          read(77,*,iostat=ierr) cardname,(ntidegateopt(ii),ii=1,numtidegate)  
          if(ierr/=0) call diag_print_error('Must specify the number of controlling elements for each tide gate')
          ntidegateopt(0)=0
          do ii=1,numtidegate
            ntidegateopt(ii)=ntidegateopt(ii-1)+ntidegateopt(ii)
          enddo
          maxtidegateopt=ntidegateopt(numtidegate)
          allocate(tidegateopt(maxtidegateopt))
          k=0
          
        case('REG_START_TIME')
          k=k+1
          backspace(77)
          read(77,*,iostat=ierr) cardname,tidegateopt(k)    
          if(ierr/=0) call diag_print_error('Must specify the regular start time (hr) for each tide gate')
          
        case('REG_OPEN_FREQUENCY')
          k=k+1
          backspace(77)
          read(77,*,iostat=ierr) cardname,tidegateopt(k) 
          if(ierr/=0) call diag_print_error('Must specify the regular opening frequency (hr) for each tide gate')
          
        case('REG_OPEN_DURATION')
          k=k+1
          backspace(77)
          read(77,*,iostat=ierr) cardname,tidegateopt(k)   
          if(ierr/=0) call diag_print_error('Must specify the regular opening duration (hr) for each tide gate')
          
        case('DES_START_TIME')
          k=k+1
          backspace(77)
          read(77,*,iostat=ierr) cardname,(tidegateopt(ii),ii=k,maxtidegateopt,2)   
          if(ierr/=0) call diag_print_error('Must specify the designated start time (hr) for each tide gate')
          
        case('DES_OPEN_DURATION')
          k=k+1
          backspace(77)
          read(77,*,iostat=ierr) cardname,(tidegateopt(ii),ii=k,maxtidegateopt,2)  
          if(ierr/=0) call diag_print_error('Must specify the designated opening duration (hr) for each tide gate')
          
      case default
        write(*,*) 'WARNING: Card ',cardname,' not found'
        foundcard = .false.
         
    end select
    enddo
    
    return
    end subroutine tg_schedule_block
    
