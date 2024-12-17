!*****************************************************************
    subroutine rubble_mound_alloc
! Resizes the rubble mound structure variables
!*****************************************************************
    use struct_def
    implicit none

    integer :: irm
    type(RM_type), allocatable :: RM_temp(:)
    
    irubmound = irubmound + 1
    if (irubmound == 1) then
      allocate(RM_struct(1))
    else
      allocate(RM_temp(irubmound-1))
      do irm=1,irubmound-1
        RM_temp(irm) = RM_struct(irm)
      enddo
      deallocate(RM_struct)
      allocate(RM_struct(irubmound))
      do irm=1,irubmound-1
        RM_struct(irm) = RM_temp(irm)
      enddo
      deallocate(RM_temp)
    endif      
    
    !Initialize and set default values
    RM_struct(irubmound)%rockdia_const       = 0.0
    RM_struct(irubmound)%structporo_const    = 0.0
    RM_struct(irubmound)%structbaseD_const   = 0.0
    RM_struct(irubmound)%rubmoundmeth        = 0
    RM_struct(irubmound)%ncells              = 0
    RM_struct(irubmound)%name                = ''
    RM_struct(irubmound)%rockdia_dset(2)     = ''
    RM_struct(irubmound)%structporo_dset(2)  = ''
    RM_struct(irubmound)%structbaseD_dset(2) = ''
    RM_struct(irubmound)%structmeth_dset(2)  = ''
    RM_struct(irubmound)%inc = 0
    
    return
    end subroutine rubble_mound_alloc

!*****************************************************************
    subroutine read_rubble_mound(cardname,foundcard)
!*****************************************************************
    use bnd_def
    use const_def, only: pi
    use flow_def, only: grav
    use struct_def
    implicit none
    integer :: k   !hli(11/28/12)
    integer :: numids, irm  !meb 01/28/2019
    character(len=34) :: cardname
    character(len=200) :: aline,adum
    logical :: foundcard
        
    foundcard = .true.
    rmblock=1

    call rubble_mound_alloc         !added meb 1/28/2019

    do
      read(77,*) cardname
      if(cardname(1:1)=='!' .or. cardname(1:1)=='#' .or. cardname(1:1)=='*') cycle
      
      select case (cardname) 
        case('RUBBLE_MOUND_END','END')
          exit
          
        case('RUBBLE_MOUND_DATASET')  !Modify input format for rubble mound structure (hli,11/26/12)
          backspace(77)
          read(77,*) cardname, arubmoundfile, arubmoundpath

        case('ROCK_DIAMETER_DATASET')
          backspace(77)
          read(77,*) cardname, arockdiamfile, arockdiampath
          RM_struct(irubmound)%rockdia_dset(1) = arockdiamfile
          RM_struct(irubmound)%rockdia_dset(2) = arockdiampath

        case('STRUCTURE_POROSITY_DATASET')
          backspace(77)
          read(77,*) cardname, astructporofile, astructporopath    
          RM_struct(irubmound)%structporo_dset(1) = astructporofile
          RM_struct(irubmound)%structporo_dset(2) = astructporopath

        case('STRUCTURE_BASE_DEPTH_DATASET')
          backspace(77)
          read(77,*) cardname, astructbaseDfile, astructbaseDpath    
          RM_struct(irubmound)%structbaseD_dset(1) = astructbaseDfile
          RM_struct(irubmound)%structbaseD_dset(2) = astructbaseDpath

        case('FORCHHEIMER_COEFF_METHOD_DATASET')
          backspace(77)
          read(77,*) cardname, astructmethfile, astructmethpath    
          RM_struct(irubmound)%structmeth_dset(1) = astructmethfile
          RM_struct(irubmound)%structmeth_dset(2) = astructmethpath

        case('CELL_IDS')
          backspace(77)
          read(77,*) cardname,numids    !No. of IDs
          RM_struct(irubmound)%ncells = numids
          allocate(RM_struct(irubmound)%cells(numids))
          backspace(77)
          read(77,*) cardname,numids,(RM_struct(irubmound)%cells(irm),irm=1,numids)
          
        case('NAME')
          backspace(77)
          read(77,'(A100)') aline
          read(aline,'(A,A100)') cardname,adum
          RM_struct(irubmound)%name = adjustl(adum)                    !meb added index 3/17/2020
          
        case('ROCK_DIAMETER_CONSTANT')
          backspace(77)
          read(77,*) cardname,RM_struct(irubmound)%rockdia_const       !meb added index 1/28/2019
    
        case('STRUCTURE_POROSITY_CONSTANT')
          backspace(77)
          read(77,*) cardname,RM_struct(irubmound)%structporo_const    !meb added index 1/28/2019
           
        case('STRUCTURE_BASE_DEPTH_CONSTANT')
          backspace(77)
          read(77,*) cardname,RM_struct(irubmound)%structbaseD_const   !meb added index 1/28/2019
  
        case('FORCHHEIMER_COEFF_METHOD')  
          backspace(77)
          read(77,*) cardname,RM_struct(irubmound)%rubmoundmeth        !meb added index 1/28/2019   
          
        case default
          foundcard = .false.
        end select              
      enddo
    end subroutine read_rubble_mound
    
