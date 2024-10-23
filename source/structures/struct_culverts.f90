!************************************************************************************
    subroutine culvert_block()
! Reads the block structure for a single culvert
! written by Alex Sanchez, USACE-CHL    
!************************************************************************************
    use const_def, only: deg2rad
    use diag_lib, only: diag_print_error                                                                !MEB 06/21  better descriptive output without application error
    use struct_def
    implicit none
    integer :: i,icv,ierr
    character(len=37) :: cardname,cdum
    character(len=37),allocatable :: strings(:)
    
    do i=1,20
      read(77,*,iostat=ierr) cardname      
      if(ierr/=0) exit
      if(cardname(1:1)=='!' .or. cardname(1:1)=='#') cycle
      select case(cardname)
        case('NUMBER_CULVERTS')
          backspace(77)
          read(77,*) cardname,numculvert
          allocate(strings(numculvert))
          allocate(idculvert(numculvert,2),aculverttype(numculvert),iculverttype(numculvert),iculvertflap(numculvert),  &
                   culvertrad(numculvert),culvertwidth(numculvert),culvertheight(numculvert),  &
                   culvertelevbay(numculvert),culvertelevsea(numculvert),culvertlength(numculvert),  &
                   cvheadlossbayentr(numculvert),cvheadlossbayexit(numculvert),cvheadlossseaentr(numculvert),  &
                   cvheadlossseaexit(numculvert),culvertfrict(numculvert),culvertmann(numculvert),  &
                   qculvert(numculvert),dqcvdzdown(numculvert),dqcvdzup(numculvert),  &
                   uvculvert(numculvert),angleculvertbay(numculvert),angleculvertsea(numculvert) ) 
          !Initialize
          qculvert=0.0   !by Wu
          dqcvdzdown=0.0
          dqcvdzup=0.0
          iculverttype=2
          aculverttype='BOX'
          iculvertflap=0
          culvertfrict=0.04
          culvertmann=0.02
          cvheadlossbayentr=0.5
          cvheadlossseaentr=0.5
          cvheadlossbayexit=0.5
          cvheadlossseaexit=0.5
          
        case('CELLS')
          backspace(77)
          read(77,*,iostat=ierr) cardname,(idculvert(icv,1),idculvert(icv,2),icv=1,numculvert)  
          if(ierr/=0) call diag_print_error('Must specify bay and sea side cell IDs for each culvert')   !MEB 06/21  better descriptive output without application error

        case('TYPE')
          backspace(77)
          read(77,*,iostat=ierr) cardname,(strings(icv),icv=1,numculvert)
          if(ierr/=0) call diag_print_error('Must specify a BOX or CIR type for each culvert')
          do icv=1,numculvert
            cdum = strings(icv)  
            if(cdum(1:3)=='CIR')then
              iculverttype(icv) = 1
            else 
              iculverttype(icv) = 2
            endif 
            aculverttype(icv) = trim(cdum)
          enddo   
          
        case('FLAP_GATE')
          backspace(77)
          read(77,*,iostat=ierr) cardname,(strings(icv),icv=1,numculvert)
          if(ierr/=0) call diag_print_error('Must specify Flap Gate ON or OFF value for each culvert')
          do icv=1,numculvert
            cdum = strings(icv)  
            if(cdum(1:2)=='ON')then
              iculvertflap(icv) = 1
            else
              iculvertflap(icv) = 0
            endif
          enddo          
        
        case('RADIUS')
          backspace(77)
          read(77,*,iostat=ierr) cardname,(culvertrad(icv),icv=1,numculvert)
          if(ierr/=0) call diag_print_error('Must specify a radius value for each culvert')
          do icv=1,numculvert
            culvertwidth(icv)=2.0*culvertrad(icv)
          enddo

        case('WIDTH')
          backspace(77)
          read(77,*,iostat=ierr) cardname,(culvertwidth(icv),icv=1,numculvert)
          if(ierr/=0) call diag_print_error('Must specify a width value for each culvert')
          do icv=1,numculvert
            culvertrad(icv)=0.5*culvertwidth(icv)
          enddo
            
        case('HEIGHT')
          backspace(77)
          read(77,*,iostat=ierr) cardname,(culvertheight(icv),icv=1,numculvert)
          if(ierr/=0) call diag_print_error('Must specify a height value for each culvert')
          
        case('LENGTH')
          backspace(77)
          read(77,*,iostat=ierr) cardname,(culvertlength(icv),icv=1,numculvert)
          if(ierr/=0) call diag_print_error('Must specify a length value for each culvert')
          
        case('DARCY_FRICTION_FACTOR')
          backspace(77)
          read(77,*,iostat=ierr) cardname,(culvertfrict(icv),icv=1,numculvert)
          if(ierr/=0) call diag_print_error('Must specify a Darcy friction value for each culvert')
          do icv=1,numculvert
            culvertfrict(icv)=max(culvertfrict(icv),0.001)
          enddo 
          
        case('MANNINGS_COEFFICIENT')
          backspace(77)
          read(77,*,iostat=ierr) cardname,(culvertmann(icv),icv=1,numculvert)
          if(ierr/=0) call diag_print_error('Must specify a Mannings coefficient value for each culvert')
          do icv=1,numculvert
            culvertmann(icv)=max(culvertmann(icv),0.001)
          enddo
            
        case('INVERT_ELEVATIONS')
          backspace(77)
          read(77,*,iostat=ierr) cardname,(culvertelevbay(icv),culvertelevsea(icv),icv=1,numculvert)
          if(ierr/=0) call diag_print_error('Must specify a Bay and Sea side invert elevation value for each culvert')
          
        case('ENTRY_HEAD_LOSSES')
          backspace(77)
          read(77,*,iostat=ierr) cardname,(cvheadlossbayentr(icv),cvheadlossseaentr(icv),icv=1,numculvert)
          if(ierr/=0) call diag_print_error('Must specify a Bay and Sea side entry head loss value for each culvert')
          do icv=1,numculvert
            cvheadlossbayentr(icv)=max(cvheadlossbayentr(icv),0.0)
            cvheadlossseaentr(icv)=max(cvheadlossseaentr(icv),0.0)
          enddo
            
        case('EXIT_HEAD_LOSSES')
          backspace(77)
          read(77,*,iostat=ierr) cardname,(cvheadlossbayexit(icv),cvheadlossseaexit(icv),icv=1,numculvert)
          if(ierr/=0) call diag_print_error('Must specify a Bay and Sea side exit head loss value for each culvert')
          do icv=1,numculvert
            cvheadlossbayexit(icv)=max(cvheadlossbayexit(icv),0.0)
            cvheadlossseaexit(icv)=max(cvheadlossseaexit(icv),0.0)
          enddo
          
        case('OUTFLOW_ANGLES')
          backspace(77)
          read(77,*,iostat=ierr) cardname,(angleculvertbay(icv),angleculvertsea(icv),icv=1,numculvert)
          if(ierr/=0) call diag_print_error('Must specify Bay and Sea outflow angle values for each culvert')
          do icv=1,numculvert
            angleculvertbay(icv)=angleculvertbay(icv)*deg2rad
            angleculvertsea(icv)=angleculvertsea(icv)*deg2rad
          enddo
          
        case('CULVERT_END','END')
          exit
          
        case default
          write(*,*) 'WARNING: Card ',trim(cardname),' not found'
          
      end select
    enddo
    
    return
    end subroutine culvert_block
  
!*****************************************************************
    subroutine culvert_alloc
! Resizes the culvert structure variables
!*****************************************************************
    use struct_def
    implicit none

    integer :: icv
    type(CV_type), allocatable :: CV_temp(:)
    
    iculvert = iculvert + 1
    if (iculvert == 1) then
      allocate(CV_struct(1))
    else
      allocate(CV_temp(iculvert-1))
      do icv=1,iculvert-1
        CV_temp(icv) = CV_struct(icv)
      enddo
      deallocate(CV_struct)
      allocate(CV_struct(iculvert))
      do icv=1,iculvert-1
        CV_struct(icv) = CV_temp(icv)
      enddo
      deallocate(CV_temp)
    endif      
    
    !Initialize and set default values
    CV_struct(iculvert)%cell_bay          = 0
    CV_struct(iculvert)%cell_sea          = 0
    CV_struct(iculvert)%culvertType       = 2
    CV_struct(iculvert)%radius            = 0.0
    CV_struct(iculvert)%height            = 0.0
    CV_struct(iculvert)%width             = 0.0
    CV_struct(iculvert)%length            = 0.0
    CV_struct(iculvert)%darcyCoef         = 0.04
    CV_struct(iculvert)%manningsCoef      = 0.025
    CV_struct(iculvert)%invertElevBay     = 0.0
    CV_struct(iculvert)%invertElevSea     = 0.0
    CV_struct(iculvert)%entryCoefWeir_b2s = 0.0
    CV_struct(iculvert)%entryCoefWeir_s2b = 0.0
    CV_struct(iculvert)%exitCoefWeir_b2s  = 0.0
    CV_struct(iculvert)%exitCoefWeir_s2b  = 0.0
    CV_struct(iculvert)%outflowAngle_bay  = 0.0
    CV_struct(iculvert)%outflowAngle_sea  = 0.0
    CV_struct(iculvert)%culvertFlap       = .false.
    
    return
    end subroutine culvert_alloc
    
!************************************************************************************
    subroutine new_culvert_block()
! Reads the multiple block structure for culverts (meb, 10/22/2024)
!   Based on new_weir_block
!************************************************************************************
    use const_def, only: deg2rad
    use struct_def
    use diag_lib
    implicit none
    
    integer :: i,icv,ierr,maxculvert
    integer :: ival,numids
    real    :: rval
    character(len=34) :: cardname
    character(len=10) :: cdum
    logical :: foundcard

    foundcard = .true.
    
    call culvert_alloc  !increment and initialize multiple culvert structures
    
    do
      read(77,*,iostat=ierr) cardname
      select case(cardname)
        case('CULVERT_STRUCT_END','END')
          exit
          
        case('CELLS')
          backspace(77)
          read(77,*,iostat=ierr) cardname, CV_struct(iculvert)%cell_bay, CV_struct(iculvert)%cell_sea
          if(ierr/=0) call diag_print_error('There must be two cells IDs specified.')  

        case('TYPE')
          backspace(77)
          read(77,*) cardname, cdum
          if (cdum(1:3) .eq. 'BOX') then
            CV_struct(iculvert)%culvertType = 2
          elseif (cdum(1:3) .eq. 'CIR') then
            CV_struct(iculvert)%culvertType = 1
          else
            call diag_print_error('Culvert type must be either BOX or CIRCLE.')
          endif 

        case('FLAP_GATE')
          call card_boolean(77,CV_struct(iculvert)%culvertFlap,ierr)

        case('RADIUS')
          if (CV_struct(iculvert)%culvertType == 2) call diag_print_error('Culvert type is BOX. RADIUS card encountered.')
          backspace(77)
          read(77,*) cardname, CV_struct(iculvert)%radius

        case('HEIGHT')
          if (CV_struct(iculvert)%culvertType == 1) call diag_print_error('Culvert type is CIRCLE. HEIGHT card encountered.')
          backspace(77)
          read(77,*) cardname, CV_struct(iculvert)%height

        case('WIDTH')
          if (CV_struct(iculvert)%culvertType == 1) call diag_print_error('Culvert type is CIRCLE. WIDTH card encountered.')
          backspace(77)
          read(77,*) cardname, CV_struct(iculvert)%width

        case('LENGTH')
          backspace(77)
          read(77,*) cardname, CV_struct(iculvert)%length

        case('DARCY_FRICTION_FACTOR')
          backspace(77)
          read(77,*) cardname, CV_struct(iculvert)%darcyCoef

        case('MANNINGS_COEFFICIENT')
          backspace(77)
          read(77,*) cardname, CV_struct(iculvert)%manningsCoef

        case('INVERT_ELEVATIONS')
          backspace(77)
          read(77,*,iostat=ierr) cardname, CV_struct(iculvert)%invertElevBay, CV_struct(iculvert)%invertElevSea
          if(ierr/=0) call diag_print_error('There must be two invert elevation values specified.')  
          
        case('ENTRY_HEAD_LOSS_BAY')
          backspace(77)
          read(77,*) cardname, CV_struct(iculvert)%entryCoefWeir_b2s

        case('ENTRY_HEAD_LOSS_SEA')
          backspace(77)
          read(77,*) cardname, CV_struct(iculvert)%entryCoefWeir_s2b

        case('EXIT_HEAD_LOSS_BAY')
          backspace(77)
          read(77,*) cardname, CV_struct(iculvert)%exitCoefWeir_b2s

        case('EXIT_HEAD_LOSS_SEA')
          backspace(77)
          read(77,*) cardname, CV_struct(iculvert)%exitCoefWeir_s2b

        case('OUTFLOW_ANGLES')
          backspace(77)
          read(77,*,iostat=ierr) cardname, CV_struct(iculvert)%outflowAngle_Bay, CV_struct(iculvert)%outflowAngle_Sea
          if(ierr/=0) call diag_print_error('There must be two outflow angles values specified.')  
          
        case default
          write(*,*) 'WARNING: Card ',cardname,' not found'
          foundcard = .false.
         
        end select
    enddo
    
    return
    end subroutine new_culvert_block

!****************************************************    
    subroutine culvert_init
!****************************************************
#include "CMS_cpp.h"
    use struct_def
    use const_def, only: deg2rad    
    implicit none
    
    integer icv, j, kk, maxculvert
    
    numculvert = iculvert
    allocate(idculvert(numculvert,2),aculverttype(numculvert),iculverttype(numculvert),             &
             iculvertflap(numculvert),culvertrad(numculvert),culvertwidth(numculvert),              &
             culvertheight(numculvert),culvertelevbay(numculvert),culvertelevsea(numculvert),       &
             culvertlength(numculvert),cvheadlossbayentr(numculvert),cvheadlossbayexit(numculvert), &
             cvheadlossseaentr(numculvert),cvheadlossseaexit(numculvert),culvertfrict(numculvert),  &
             culvertmann(numculvert),qculvert(numculvert),dqcvdzdown(numculvert),                   &
             dqcvdzup(numculvert),uvculvert(numculvert),angleculvertbay(numculvert),                &
             angleculvertsea(numculvert) )
    
    qculvert=0.0   !by Wu
    dqcvdzdown=0.0
    dqcvdzup=0.0
    iculverttype=2
    aculverttype='BOX'
    iculvertflap=0
    culvertfrict=0.04
    culvertmann=0.02
    cvheadlossbayentr=0.5
    cvheadlossseaentr=0.5
    cvheadlossbayexit=0.5
    cvheadlossseaexit=0.5
    
    do icv=1,numculvert
      idculvert(icv,1) = CV_struct(icv)%cell_bay
      idculvert(icv,2) = CV_struct(icv)%cell_sea
      if (CV_struct(icv)%culvertType == 1) then
        iculverttype(icv)  = 1
        aculverttype(icv)  = 'CIR'
        culvertrad(icv)    = CV_struct(icv)%radius
        culvertwidth(icv)  = 2.0 * CV_struct(icv)%radius
        culvertheight(icv) = 2.0 * CV_struct(icv)%radius
      else
        iculverttype(icv)  = 2
        aculverttype(icv)  = 'BOX'
        culvertrad(icv)    = 0.5 * CV_struct(icv)%width
        culvertwidth(icv)  = CV_struct(icv)%width
        culvertheight(icv) = CV_struct(icv)%height
      endif
      if (CV_struct(icv)%culvertFlap) then
        iculvertflap(icv) = 1
      else 
        iculvertflap(icv) = 0
      endif
      culvertlength(icv) = CV_struct(icv)%length
      culvertfrict(icv) = max(CV_struct(icv)%darcyCoef,0.001)
      culvertmann(icv)  = max(CV_struct(icv)%manningsCoef, 0.001)
      culvertelevbay(icv) = CV_struct(icv)%invertElevBay
      culvertelevsea(icv) = CV_struct(icv)%invertElevSea
      cvheadlossbayentr(icv) = max(CV_struct(icv)%entryCoefWeir_b2s,0.0)
      cvheadlossseaentr(icv) = max(CV_struct(icv)%entryCoefWeir_s2b,0.0)
      cvheadlossbayexit(icv) = max(CV_struct(icv)%exitCoefWeir_b2s,0.0)
      cvheadlossseaexit(icv) = max(CV_struct(icv)%exitCoefWeir_s2b,0.0)
      angleculvertbay(icv) = CV_struct(icv)%outflowAngle_bay * deg2rad
      angleculvertsea(icv) = CV_struct(icv)%outflowAngle_sea * deg2rad
    enddo
        
    return
    end subroutine culvert_init
       
