!***********************************************************************    
    subroutine veg_default
!***********************************************************************       
    use veg_def
    implicit none
    
    !Vegetation Parameters
    hv = 1.0      !Height of vegetation
    Cdv = 1.0     !Drag coefficient
    Dv = 0.0064   !Stem diameter in meters
    alphac = 0.75
    
    end subroutine veg_default

!***********************************************************************    
    subroutine veg_cards(cardname,foundcard)
!***********************************************************************       
    use veg_def
    implicit none
    character(len=32) :: cardname,cdum
    logical :: foundcard
    
    foundcard = .true.
    select case (cardname)          
      case('INCLUDE_VEGETATION') !Alex
        backspace(77)
        read(77,*) cardname, cdum
        if(cdum.eq.'ON ') veg = .true.              
          
      case('VEGETATION_STEM_DENSITY')  !Alex        
        backspace(77)
        read(77,*) cardname, vegdensfile, vegdenspath      
            
      case('VEGETATION_DIAMETER') !Alex
        backspace(77)
        read(77,*) cardname, Dv    
  
      case('VEGETATION_HEIGHT') !Alex
        backspace(77)
        read(77,*) cardname, Hv 
        
      case default
        foundcard = .false.  

    end select
    
    end subroutine veg_cards
    
!***********************************************************************
    subroutine veg_init
!***********************************************************************    
#include "CMS_cpp.h"    
    use veg_def
    use size_def
    use diag_lib, only: diag_print_error
#ifdef XMDF_IO  
    use in_xmdf_lib, only: readscalh5
#endif   
    implicit none
    integer :: ierr
    
    allocate(nv(ncellsD))    
    
#ifdef XMDF_IO   
    call readscalh5(vegdensfile,vegdenspath,nv,ierr)
#else
    call diag_print_error('Cannot read vegetation density without XMDF')
#endif
    
    end subroutine veg_init
    
!***********************************************************************
    subroutine veg_print()
!***********************************************************************    
    use veg_def
    use diag_def, only: dgunit,dgfile
    implicit none
    integer :: i,iunit(2)
    
    iunit = (/6, dgunit/)
    
888 format(' ',A,T40,A)
    
    open(dgunit,file=dgfile,access='append')     
    do i=1,2
      write(iunit(i),*)     
      if(veg)then
        write(iunit(i),888) ' Vegetation Flow Drag:','ON'
      else    
        write(iunit(i),888) ' Vegetation Flow Drag:','OFF'
      endif  
    enddo       
    close(dgunit)
    
    end subroutine veg_print    
    
!***********************************************************************
    subroutine veg_drag
!***********************************************************************
    use size_def
    use geo_def, only: areap
    use flow_def, only: uv,h,sp
    use veg_def
    use prec_def
    implicit none
    integer :: i
    real(ikind) :: cd,val
    
    val = alphac*0.5*Cdv*Dv
    do i=1,ncells
      cd=val*nv(i)*(min(hv,h(i))**2)/h(i)
      sp(i)=sp(i)-cd*uv(i)*areap(i)
    enddo
    
    end subroutine veg_drag
