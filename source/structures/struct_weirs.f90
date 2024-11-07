!*****************************************************************
    subroutine weir_alloc
! Resizes the weir structure variables
!*****************************************************************
    use struct_def
    implicit none

    integer :: iwr
    type(WR_type), allocatable :: WR_temp(:)
    
    iweir = iweir + 1
    if (iweir == 1) then
      allocate(WR_struct(1))
    else
      allocate(WR_temp(iweir-1))
      do iwr=1,iweir-1
        WR_temp(iwr) = WR_struct(iwr)
      enddo
      deallocate(WR_struct)
      allocate(WR_struct(iweir))
      do iwr=1,iweir-1
        WR_struct(iwr) = WR_temp(iwr)
      enddo
      deallocate(WR_temp)
    endif      
    
    !Initialize and set default values
    WR_struct(iweir)%coefWeirLateral = 0.0
    WR_struct(iweir)%ncells       = 0
    WR_struct(iweir)%orientWeir   = 1
    WR_struct(iweir)%weirType     = 1
    WR_struct(iweir)%coefWeir_b2s = 0.0
    WR_struct(iweir)%coefWeir_s2b = 0.0
    WR_struct(iweir)%elevWeir     = 0.0
    WR_struct(iweir)%methWeir     = 1
    
    return
    end subroutine weir_alloc
    
!************************************************************************************
    subroutine new_weir_block()
! Reads the multiple block structure for weirs (meb, 01/28/2019)
!   Based on weir_block (hli, 02/08/13) 
!************************************************************************************
    use const_def, only: deg2rad
    use struct_def
    use diag_lib
    implicit none
    
    integer :: i,iwr,ierr,maxweir
    integer :: ival,numids
    real    :: rval
    character(len=34) :: cardname
    character(len=5)  :: cdum
    logical :: foundcard

    foundcard = .true.
    
    call weir_alloc  !increment and initialize multiple weir structure
    
    do
      read(77,*,iostat=ierr) cardname
      select case(cardname)
        case('WEIR_STRUCT_END','END')
          exit
          
        case('NUM_CELL_WEIRS','NUM_CELL_WEIR')
          backspace(77)
          read(77,*) cardname,ival
          WR_struct(iweir)%ncells = ival
          allocate(WR_struct(iweir)%cells(ival))

        case('CELLS')
          backspace(77)
          ival = WR_struct(iweir)%ncells 
          read(77,*,iostat=ierr) cardname, (WR_struct(iweir)%cells(iwr),iwr=1,ival)
          if(ierr/=0) call diag_print_error('Must specify the number of cells for each weir')            !MEB 06/21  better descriptive output without application error

        case('CELL_IDS')
          backspace(77)
          read(77,*) cardname,numids    !No. of IDs
          WR_struct(iweir)%ncells = numids
          allocate(WR_struct(iweir)%cells(numids))
          backspace(77)
          read(77,*) cardname,numids,(WR_struct(iweir)%cells(iwr),iwr=1,numids)
          
        case('DISTRIBUTION_COEFFICIENT')
          backspace(77)
          read(77,*) cardname, WR_struct(iweir)%coefweirlateral

        case('ORIENTATION')
          backspace(77)
          read(77,*) cardname, ival
          if (ival .lt. 1 .or. ival .gt. 4) then
            call diag_print_error('Orientation value must be between 1 and 4, inclusively')
          endif 
          WR_struct(iweir)%orientWeir = ival

        case('ORIENTATION_SEA')
          backspace(77)
          read(77,*) cardname, cdum
          selectcase(cdum(1:4))
          case('NORT')
            WR_struct(iweir)%orientWeir = 1
          case('EAST')
            WR_struct(iweir)%orientWeir = 2
          case('SOUT')
            WR_struct(iweir)%orientWeir = 3
          case('WEST')
            WR_struct(iweir)%orientWeir = 4
          case default
            call diag_print_error('Orientation value must be one of the four main cardinal directions')
          end select
            
        case('TYPE')
          backspace(77)
          read(77,*) cardname, ival 
          if (ival .lt. 1 .or. ival .gt. 2) then
            call diag_print_error('Weir type must be either 1 or 2')
          endif 
          WR_struct(iweir)%weirType = ival
          
        case('CREST_TYPE')
          backspace(77)
          read(77,*) cardname, cdum
          ival = 1
          if (cdum(1:5) == 'BROAD') ival = 2
          WR_struct(iweir)%weirType = ival          

        case('FLOW_COEFFICIENT_FROM_BAY','FLOW_COEFF_FROM_BAY')
          backspace(77)
          read(77,*) cardname, rval
          WR_struct(iweir)%coefWeir_b2s = rval
          
        case('FLOW_COEFFICIENT_FROM_SEA','FLOW_COEFF_FROM_SEA')
          backspace(77)
          read(77,*) cardname, rval
          WR_struct(iweir)%coefWeir_s2b = rval
          
        case('CREST_ELEVATION')
          backspace(77)
          read(77,*) cardname, rval
          WR_struct(iweir)%elevWeir = rval

        case('METH','METHOD')
          backspace(77)
          read(77,*) cardname, ival 
          if (ival .lt. 1 .or. ival .gt. 2) then
            call diag_print_error('Method must be either 1 or 2')
          endif 
          WR_struct(iweir)%methWeir = ival
          
        case default
          write(*,*) 'WARNING: Card ',cardname,' not found'
          foundcard = .false.
         
        end select
    enddo
    
    return
    end subroutine new_weir_block
    
    
!************************************************************************************
    subroutine weir_block()
! Reads the block structure for weirs (hli, 02/08/13) 
!************************************************************************************
    use const_def, only: deg2rad
    use struct_def
    use diag_lib, only: diag_print_error                                                                 !MEB 06/21  better descriptive output without application error
    implicit none
    integer :: i,iwr,ierr,maxweir
    character(len=34) :: cardname
    logical :: foundcard

    foundcard = .true.
    do i=1,30
      read(77,*,iostat=ierr) cardname
      select case(cardname)
        case('NUMBER_WEIRS','NUMBER_WEIR')
          backspace(77)
          read(77,*) cardname,numweir
          allocate(nweir(0:numweir),orientweir(numweir),iweirtype(numweir),methweir(numweir),  &
                   coefweir(numweir,2),elevweir(numweir),orientweirbay(numweir),  &
                   orientweirsea(numweir),Qtotweir(numweir) )

        case('NUM_CELL_WEIRS','NUM_CELL_WEIR')
          backspace(77)
          read(77,*,iostat=ierr) cardname,(nweir(iwr),iwr=1,numweir)              
          if(ierr/=0) call diag_print_error('Must specify the number of cells for each weir')            !MEB 06/21  better descriptive output without application error
          nweir(0) = 0
          maxweir = 0
          do iwr=1,numweir
            maxweir = maxweir + nweir(iwr)
          enddo
          do iwr=1,numweir
            nweir(iwr) = nweir(iwr-1) + nweir(iwr)
          enddo

          allocate(idweir(maxweir),qweir(maxweir),coefweirlateral(maxweir),dqweirdzdown(maxweir),dqweirdzup(maxweir))

        case('CELLS')
          backspace(77)
          read(77,*,iostat=ierr) cardname,(idweir(iwr),iwr=1,maxweir)     
          if(ierr/=0) call diag_print_error('Must specify the Cell IDs in each weir')

        case('DISTRIBUTION_COEFFICIENT')
          backspace(77)
          read(77,*,iostat=ierr) cardname,(coefweirlateral(iwr),iwr=1,maxweir)     
          if(ierr/=0) call diag_print_error('Must specify the distribution coefficient values for each weir')

        case('ORIENTATION')
          backspace(77)
          read(77,*,iostat=ierr) cardname,(orientweir(iwr),iwr=1,numweir)     
          if(ierr/=0) call diag_print_error('Must specify the orientation for each weir')

        case('TYPE')
          backspace(77)
          read(77,*,iostat=ierr) cardname,(iweirtype(iwr),iwr=1,numweir)     
          if(ierr/=0) call diag_print_error('Must specify the type for each weir')

        case('FLOW_COEFFICIENT')
          backspace(77)
          read(77,*,iostat=ierr) cardname,(coefweir(iwr,1),coefweir(iwr,2),iwr=1,numweir)     
          if(ierr/=0) call diag_print_error('Must specify the flow coefficient for each weir')

        case('CREST_ELEVATION')
          backspace(77)
          read(77,*,iostat=ierr) cardname,(elevweir(iwr),iwr=1,numweir)     
          if(ierr/=0) call diag_print_error('Must specify the crest elevation for each weir')

        case('METH')
          backspace(77)
          read(77,*,iostat=ierr) cardname,(methweir(iwr),iwr=1,numweir)     
          if(ierr/=0) call diag_print_error('Must specify the flux calculation method for each weir')

        !Initialize       
        dqweirdzdown = 0.0;  dqweirdzup = 0.0  
     
        do iwr=1,numweir
          orientweirsea(iwr) = orientweir(iwr)   !Define direction of sea side at local coordinate
          if(orientweir(iwr).eq.1) then
            orientweirbay(iwr) = 3
          elseif(orientweir(iwr).eq.2) then
            orientweirbay(iwr) = 4
          elseif(orientweir(iwr).eq.3) then
            orientweirbay(iwr) = 1
          elseif(orientweir(iwr).eq.4) then
            orientweirbay(iwr) = 2
          endif
        enddo

      case('WEIR_END','END')
        exit
          
      case default
        write(*,*) 'WARNING: Card ',cardname,' not found'
        foundcard = .false.
         
      end select
    enddo
    
    return
  end subroutine weir_block

!****************************************************    
    subroutine weir_init
!****************************************************
#include "CMS_cpp.h"
    use struct_def
    implicit none
    
    integer iwr, j, kk, maxweir
    
    numweir = iweir
    allocate(nweir(0:numweir),orientweir(numweir),iweirtype(numweir),methweir(numweir),  &
             coefweir(numweir,2),elevweir(numweir),orientweirbay(numweir),  &
             orientweirsea(numweir),Qtotweir(numweir) )
    allocate(aweirtype(numweir),aorientweirsea(numweir))
    
    nweir(0) = 0
    maxweir = 0
    do iwr=1,numweir
      maxweir = maxweir + WR_struct(iwr)%ncells
      nweir(iwr) = maxweir
    enddo
    !do iwr=1,numweir
    !  maxweir = maxweir + nweir(iwr)           !Combined this into loop above
    !enddo
    !do iwr=1,numweir
    !  nweir(iwr) = nweir(iwr-1) + nweir(iwr)   !Combined this into loop above
    !enddo

    allocate(idweir(maxweir),qweir(maxweir),coefweirlateral(maxweir),dqweirdzdown(maxweir),dqweirdzup(maxweir))

    kk=0
    do iwr=1,numweir
      do j = 1,WR_struct(iwr)%ncells
        kk=kk+1
        idweir(kk) = WR_struct(iwr)%cells(j)
        coefweirlateral(kk) = WR_struct(iwr)%coefWeirLateral
      enddo
      orientweir(iwr) = WR_struct(iwr)%orientWeir
      iweirtype(iwr)  = WR_struct(iwr)%weirType
      coefweir(iwr,1) = WR_struct(iwr)%coefWeir_b2s
      coefweir(iwr,2) = WR_struct(iwr)%coefWeir_s2b
      elevweir(iwr)   = WR_struct(iwr)%elevWeir
      methweir(iwr)   = WR_struct(iwr)%methWeir
      selectcase(iweirtype(iwr))
      case(1)
        aweirtype(iwr) = 'SHARP'
      case(2)
        aweirtype(iwr) = 'BROAD'
      end select
    enddo

    !Initialize       
    dqweirdzdown = 0.0;  dqweirdzup = 0.0  
     
    do iwr=1,numweir
      orientweirsea(iwr) = orientweir(iwr)   !Define direction of sea side at local coordinate
      if(orientweir(iwr).eq.1) then
        orientweirbay(iwr) = 3
        aorientweirsea(iwr) = 'NORTH'
      elseif(orientweir(iwr).eq.2) then
        orientweirbay(iwr) = 4
        aorientweirsea(iwr) = 'EAST '
      elseif(orientweir(iwr).eq.3) then
        orientweirbay(iwr) = 1
        aorientweirsea(iwr) = 'SOUTH'
      elseif(orientweir(iwr).eq.4) then
        orientweirbay(iwr) = 2
        aorientweirsea(iwr) = 'WEST '
      endif
    enddo
        
    return
    end subroutine weir_init
  
  