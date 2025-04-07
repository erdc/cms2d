!======================================================================
module in_xmdf_lib
! Input Library
!
! Contains the following routines:
!   SMS file formats:
!     read_xys       - Reads an SMS XY Series (*.xys) file
!     read_tsd       - Reads an SMS Time-Series Data (*.tsd) file
!     read_dat       - Reads an SMS Dataset (*.dat) file
!   ADCIRC solution files:
!     readstep63 - Reads an ADCIRC fort.63 solution step
!     readstep64 - Reads an ADCIRC fort.64 solution step
!
! written by Alex Sanchez, USACE-CHL
!======================================================================
#include "CMS_cpp.h"
    implicit none

contains
    
#ifdef XMDF_IO
!*************************************************************      
    subroutine read_dataseth5(afile,apath,aname,nn,vec)
! Reads a dataset from an XMDF file
! written by Alex Sanchez, USACE-CHL
!*************************************************************  
    use xmdf
    use prec_def
    use diag_lib, only: diag_print_error
    use diag_def, only: msg
    implicit none
    
    !Input/Output
    character(len=*),intent(in) :: afile,apath,aname
    integer,intent(out) :: nn
    real(ikind),intent(inout),pointer :: vec(:)
    
    !Internal variables
    integer :: pid,gid,ierr,iloc
    real(4),allocatable :: ftemp(:) !Should be single    
    character(100) :: thepath
       
    call XF_OPEN_FILE(afile,READONLY,pid,ierr)
    msg = "Unable to open file: '" // trim(afile) // "'"
    if (ierr < 0) call diag_print_error (msg)
    
    thepath = trim(apath)
    iloc=index(thepath,'\')
    if (iloc.gt.0) thepath(iloc:iloc)='/'
    call XF_OPEN_GROUP(pid,trim(thepath),gid,ierr)
    msg = "Unable to open dataset: '" // trim(thepath) // "from file: '" // trim(afile) // "'"
    if (ierr < 0) call diag_print_error (msg)

    call XF_GET_PROPERTY_NUMBER(gid,trim(aname),nn,ierr)
    msg = "Unable to obtain number of times from file: '" // trim(afile) // "'"
    if (ierr < 0) call diag_print_error (msg)

    allocate(vec(nn),ftemp(nn))
    call XF_READ_PROPERTY_FLOAT(gid,trim(aname),nn,ftemp(1),ierr)
    msg = "Unable to get list of times from file: '" // trim(afile) // "'"
    if (ierr < 0) call diag_print_error (msg)

    vec = ftemp !Useful for converting precision    
    call XF_CLOSE_GROUP(gid,ierr)
    call XF_CLOSE_FILE(pid,ierr)    
    deallocate(ftemp)  
    
    return
    end subroutine read_dataseth5
    

!************************************************************************
    subroutine readscalh5(afile,apath,var,ierr)
! Reads a scalar dataset located in apath from afile XMDF file 
! written by Alex Sanchez, USACE-CHL
!************************************************************************
    use size_def, only: ncellsD,ncellsfull,ncellpoly
    use geo_def, only: idmap
    use interp_lib, only: interp_scal_node2cell
    use diag_lib
    use xmdf
    use prec_def
    implicit none
    
    !Input/Output
    character(len=*), intent(in) :: afile,apath
    real(ikind), intent(out) :: var(ncellsD)
    integer, intent(out) :: ierr
    !Internal Variables
    integer :: fid,gid,iloc
    real(4) :: vtemp(ncellsfull)
    character(len=200) :: msg2,msg3,thepath

    call XF_OPEN_FILE(trim(afile),READONLY,fid,ierr)        
    if(ierr<0) call diag_print_error('Could not open file: ',trim(afile))
    
    thepath = trim(apath)
    iloc=index(thepath,'\')
    if (iloc.gt.0) thepath(iloc:iloc)='/'
    call XF_OPEN_GROUP(fid,trim(thepath),gid,ierr)
    if(ierr<0)then
      call XF_CLOSE_FILE(fid,ierr)
      write(msg2,*,iostat=ierr) '  File: ',trim(afile)
      write(msg3,*,iostat=ierr) '  Path: ',trim(thepath)
      call diag_print_error('Could not open dataset from',msg2,msg3)
      return
    endif
       
    call XF_READ_SCALAR_VALUES_TIMESTEP(gid,1,ncellsfull,vtemp,ierr)
    if(ierr<0)then
      !call XF_CLOSE_GROUP(gid,ierr)
      !call XF_CLOSE_FILE(fid,ierr)  
      !write(msg2,*,iostat=ierr) '  File: ',trim(afile)
      !write(msg3,*,iostat=ierr) '  Path: ',trim(thepath)
      !call diag_print_error('Could not read scalar dataset from',msg2,msg3)
      return
    endif
    
    call XF_CLOSE_GROUP(gid,ierr)
    call XF_CLOSE_FILE(fid,ierr)
                
    if(ncellpoly>0)then
      call interp_scal_node2cell(vtemp,var) !Interpolate node to cell centers
    else
      call map_scal_full2active(vtemp,var) !Convert from full to active grid 
    endif
      
    return
    end subroutine readscalh5
      
!************************************************************************
    subroutine readscalsteph5(afile,apath,itsind,thrs,var,ierr)
! Reads the a time step of a scalar located in apath from afile XMDF file 
! written by Alex Sanchez, USACE-CHL
!************************************************************************
    use size_def, only: ncellsD,ncellsfull,ncellpoly
    use geo_def, only: idmap
    use interp_lib, only: interp_scal_node2cell
    use diag_def
    use diag_lib
    use xmdf
    use prec_def
    implicit none
    !Input/Output
    character(len=*), intent(in) :: afile,apath
    integer, intent(in) :: itsind !Time step index
    real(ikind), intent(out):: var(ncellsD),thrs
    integer, intent(out) :: ierr    
    !Internal Variables
    integer :: fid,gid,ntimes,iloc
    real(8),allocatable :: timesd(:) !Output times
    real(4) :: vtemp(ncellsfull) !Must be single
    character(100) :: thepath
    
    call XF_OPEN_FILE(trim(afile),READONLY,fid,ierr)        
    if(ierr<0)then
      call diag_print_error('Could not open file: ',afile)
    endif
          
    thepath = trim(apath)
    iloc=index(thepath,'\')
    if (iloc.gt.0) thepath(iloc:iloc)='/'
    call XF_OPEN_GROUP(fid,trim(thepath),gid,ierr)
    if(ierr<0)then
      call XF_CLOSE_FILE(fid,ierr)
      !write(msg2,*,iostat=ierr) '  File: ',trim(afile)
      !write(msg3,*,iostat=ierr) '  Path: ',trim(thepath)
      !call diag_print_warning('Could not open dataset from',msg2,msg3)
      ierr = -2 !Could not open group
      return
    endif     
       
    call XF_READ_SCALAR_VALUES_TIMESTEP(gid,itsind,ncellsfull,vtemp,ierr)
    if(ierr<0)then
      call XF_CLOSE_GROUP(gid,ierr)
      call XF_CLOSE_FILE(fid,ierr)  
      !write(msg2,*,iostat=ierr) '  File: ',trim(afile)
      !write(msg3,*,iostat=ierr) '  Path: ',trim(thepath)
      !write(msg4,*,iostat=ierr) '  Time Step: ',itsind
      !call diag_print_warning('Could not read scalar time step',msg2,msg3,msg4)
      ierr = 3 !Could not read timestep value
      return
    endif 

    call XF_GET_DATASET_NUM_TIMES(gid,ntimes,ierr)
    allocate(timesd(ntimes))
    call XF_GET_DATASET_TIMES(gid,ntimes,timesd,ierr)
    if(ierr<0)then
      thrs = -999.0
      write(msg2,*,iostat=ierr) '  File: ',trim(afile)
      write(msg3,*,iostat=ierr) '  Path: ',trim(thepath)
      call diag_print_warning('Could not read time stamp from',msg2,msg3)
      ierr = 4 !Could not read time stamp
    else
      thrs = timesd(itsind)
    endif
    
    call XF_CLOSE_GROUP(gid,ierr)
    call XF_CLOSE_FILE(fid,ierr)
           
    if(ncellpoly>0)then
      call interp_scal_node2cell(vtemp,var) !Interpolate node to cell centers
    else
      call map_scal_full2active(vtemp,var) !Convert from full to active grid 
    endif
      
    return
    end subroutine readscalsteph5

!************************************************************************
    subroutine readscallasth5(afile,apath,ntimes,var,reftimed,thrs,ierr)
! Reads a scalar located in apath from afile XMDF file 
! written by Alex Sanchez, USACE-CHL
!************************************************************************
!#include "CMS_cpp.h"
    use size_def, only: ncellsD,ncellsfull,ncellpoly
    use geo_def, only: idmap
    use interp_lib, only: interp_scal_node2cell
    use diag_def
    use diag_lib
    use xmdf
    use prec_def
    implicit none
    !Input/Output
    character(len=*), intent(in) :: afile,apath
    real(ikind), intent(out) :: var(ncellsD),thrs    
    real(8),intent(out) :: reftimed          !Reference time 
    integer, intent(out) :: ierr,ntimes
    !Internal Variables
    integer :: fid,gid,iloc
    real(8), allocatable :: timed(:) !Output time
    real(4) :: vtemp(ncellsfull) !Must be single
    character(100) :: thepath
    
    call XF_OPEN_FILE(trim(afile),READONLY,fid,ierr)   
    if(ierr<0) call diag_print_error('Could not open file: ',afile)
      
    thepath = trim(apath)
    iloc=index(thepath,'\')
    if (iloc.gt.0) thepath(iloc:iloc)='/'
    call XF_OPEN_GROUP(fid,trim(thepath),gid,ierr)
    if(ierr<0)then
      call XF_CLOSE_FILE(fid,ierr)
!      write(msg2,*,iostat=ierr) '  File: ',trim(afile)
!      write(msg3,*,iostat=ierr) '  Path: ',trim(thepath)
!      call diag_print_warning('Could not open dataset from ',msg2,msg3)
      ierr = -2
      return
    endif

    call XF_GET_DATASET_NUM_TIMES(gid,ntimes,ierr)
    if(ierr<0)then
      call XF_CLOSE_GROUP(gid,ierr)
      call XF_CLOSE_FILE(fid,ierr)
      write(msg2,*,iostat=ierr) '  File: ',trim(afile)
      write(msg3,*,iostat=ierr) '  Path: ',trim(thepath)
      call diag_print_error('Could not read number of times from',msg2,msg3)
    endif
    
    call XF_READ_SCALAR_VALUES_TIMESTEP(gid,ntimes,ncellsfull,vtemp,ierr)
    if(ierr<0)then
      call XF_CLOSE_GROUP(gid,ierr)
      call XF_CLOSE_FILE(fid,ierr)  
      write(msg2,*,iostat=ierr) '  File: ',trim(afile)
      write(msg3,*,iostat=ierr) '  Path: ',trim(thepath)
      write(msg4,*,iostat=ierr) '  Time step: ',ntimes
      call diag_print_error('Could not read last time step scalar values from',msg2,msg3,msg4)
    endif
    
    allocate(timed(ntimes))
    call XF_GET_DATASET_TIMES(gid,ntimes,timed,ierr)
    if(ierr<0)then
      thrs = -999.0
      write(msg2,*,iostat=ierr) '  File: ',trim(afile)
      write(msg3,*,iostat=ierr) '  Path: ',trim(thepath)
      write(msg4,*,iostat=ierr) '  Time step: ',ntimes
      call diag_print_warning('Could not read times from',msg2,msg3,msg4)
      ierr = 4 !Could not read time stamp
    else
      thrs = timed(ntimes)
    endif
    deallocate(timed)
        
    call XF_GET_DATASET_REFTIME(gid,reftimed,ierr)
    if(ierr<0)then
      reftimed = -999.0d0
      write(msg2,*,iostat=ierr) '  File: ',trim(afile)
      write(msg3,*,iostat=ierr) '  Path: ',trim(thepath)
      call diag_print_warning('Could not read reference time from',msg2,msg3)
      ierr = 5 !Could not read time stamp
    endif  

    call XF_CLOSE_GROUP(gid,ierr)
    call XF_CLOSE_FILE(fid,ierr)
            
    if(ncellpoly>0)then
      call interp_scal_node2cell(vtemp,var) !Interpolate node to cell centers
    else
      call map_scal_full2active(vtemp,var) !Convert from full to active grid 
    endif
      
    return
    end subroutine readscallasth5
    
!************************************************************************
    subroutine readscaltimeh5(afile,apath,thrs,var,reftimed,ierr)
! Reads a scalar located in apath from afile XMDF file 
! written by Alex Sanchez, USACE-CHL
!************************************************************************
!#include "CMS_cpp.h"
    use size_def, only: ncellsD,ncellsfull,ncellpoly
    use geo_def, only: idmap
    use interp_lib, only: interp_scal_node2cell
    use diag_def
    use diag_lib
    use xmdf
    use prec_def
    implicit none
    !Input
    character(len=*), intent(in) :: afile,apath
    real(ikind), intent(in) :: thrs
    !Output
    real(ikind), intent(out) :: var(ncellsD)    
    real(8),intent(out) :: reftimed          !Reference time 
    integer, intent(out) :: ierr
    !Internal Variables
    integer :: i,fid,gid,nstep,ntimes,iloc
    real(8), allocatable :: timed(:) !Output time
    real(8) :: thrsd,terrd,terrdmin
    real(4) :: vtemp(ncellsfull) !Must be single
    character(100) :: thepath
    
    call XF_OPEN_FILE(trim(afile),READONLY,fid,ierr)   
    if(ierr<0) call diag_print_error('Could not open file: ',afile)
      
    thepath = trim(apath)
    iloc=index(thepath,'\')
    if (iloc.gt.0) thepath(iloc:iloc)='/'
    call XF_OPEN_GROUP(fid,trim(thepath),gid,ierr)
    if(ierr<0)then
      call XF_CLOSE_FILE(fid,ierr)
!#ifdef DIAG_MODE
!      write(msg2,*,iostat=ierr) '  File: ',trim(afile)
!      write(msg3,*,iostat=ierr) '  Path: ',trim(thepath)
!      call diag_print_warning('Could not open dataset from ',msg2,msg3)
!#endif
      ierr = -2
      return
    endif

    call XF_GET_DATASET_NUM_TIMES(gid,ntimes,ierr)
    if(ierr<0)then
      call XF_CLOSE_GROUP(gid,ierr)
      call XF_CLOSE_FILE(fid,ierr)
      write(msg2,*,iostat=ierr) '  File: ',trim(afile)
      write(msg3,*,iostat=ierr) '  Path: ',trim(thepath)
      call diag_print_error('Could not read number of times from',msg2,msg3)
    endif
    
    allocate(timed(ntimes))
    call XF_GET_DATASET_TIMES(gid,ntimes,timed,ierr)
    if(ierr<0)then
      write(msg2,*,iostat=ierr) '  File: ',trim(afile)
      write(msg3,*,iostat=ierr) '  Path: ',trim(thepath)
      write(msg4,*,iostat=ierr) '  Time step: ',ntimes
      call diag_print_warning('Could not read times from',msg2,msg3,msg4)
      ierr = 4 !Could not read time stamp
    else
      nstep = 1
      terrdmin = 1.0e9
      thrsd = dble(thrs)
      do i=1,ntimes
        terrd = abs(thrsd-timed(i))
        if(terrd<terrdmin)then
          nstep = i
          terrdmin = terrd
        endif
      enddo
      if(terrdmin>0.001)then
        write(msg2,*,iostat=ierr) '  File: ',trim(afile)
        write(msg3,*,iostat=ierr) '  Path: ',trim(thepath)
        write(msg4,*,iostat=ierr) '  Time: ',thrs,' hrs'
        call diag_print_warning('Could not find time: ',msg2,msg3,msg4)
      endif
    endif
    deallocate(timed)
    
    call XF_READ_SCALAR_VALUES_TIMESTEP(gid,nstep,ncellsfull,vtemp,ierr)
    if(ierr<0)then
      call XF_CLOSE_GROUP(gid,ierr)
      call XF_CLOSE_FILE(fid,ierr)  
      write(msg2,*,iostat=ierr) '  File: ',trim(afile)
      write(msg3,*,iostat=ierr) '  Path: ',trim(thepath)
      write(msg4,*,iostat=ierr) '  Time step: ',ntimes
      call diag_print_error('Could not read last time step scalar values from',msg2,msg3,msg4)
    endif
    
    call XF_GET_DATASET_REFTIME(gid,reftimed,ierr)
    if(ierr<0)then
      reftimed = -999.0d0
      write(msg2,*,iostat=ierr) '  File: ',trim(afile)
      write(msg3,*,iostat=ierr) '  Path: ',trim(thepath)
      call diag_print_warning('Could not read reference time from',msg2,msg3)
      ierr = 5 !Could not read time stamp
    endif  

    call XF_CLOSE_GROUP(gid,ierr)
    call XF_CLOSE_FILE(fid,ierr)
            
    if(ncellpoly>0)then
      call interp_scal_node2cell(vtemp,var) !Interpolate node to cell centers
    else
      call map_scal_full2active(vtemp,var) !Convert from full to active grid 
    endif
      
    return
    end subroutine readscaltimeh5  

!************************************************************************
    subroutine readvech5(afile,apath,vecx,vecy,ierr)
! Reads a vector located in apath from afile XMDF file 
! written by Alex Sanchez, USACE-CHL
!************************************************************************
    use size_def
    use geo_def, only: idmap
    use interp_lib, only: interp_vec_node2cell
    use diag_lib
    use xmdf
    use prec_def
    implicit none
    !Input/Output
    character(len=200), intent(in) :: afile,apath
    real(ikind), intent(out) :: vecx(ncellsD),vecy(ncellsD)
    integer, intent(out) :: ierr
    !Internal
    integer :: fid,gid,iloc
    real(4) :: vtemp(ncellsfull*2) !Must be single
    character(len=200) :: msg2,msg3,thepath

    call XF_OPEN_FILE(trim(afile),READONLY,fid,ierr)        
    if(ierr<0) call diag_print_error('Could not open file: ',afile)

    thepath = trim(apath)
    iloc=index(thepath,'\')
    if (iloc.gt.0) thepath(iloc:iloc)='/'    
    call XF_OPEN_GROUP(fid,trim(thepath),gid,ierr)
    if(ierr<0)then
      call XF_CLOSE_FILE(fid,ierr)  
      write(msg2,*,iostat=ierr) '  File: ',trim(afile)
      write(msg3,*,iostat=ierr) '  Path: ',trim(thepath)
      call diag_print_error('Could not open dataset from',msg2,msg3)
      return
    endif
    
    call XF_READ_VECTOR_VALUES_TIMESTEP(gid,1,ncellsfull,2,vtemp,ierr)     
    if(ierr<0)then
      call XF_CLOSE_GROUP(gid,ierr)
      call XF_CLOSE_FILE(fid,ierr)
      write(msg2,*,iostat=ierr) '  File: ',trim(afile)
      write(msg3,*,iostat=ierr) '  Path: ',trim(thepath)
      call diag_print_error('Could not read vector dataset from',msg2,msg3) 
      return
    endif

    call XF_CLOSE_GROUP(gid,ierr)
    call XF_CLOSE_FILE(fid,ierr)
    
    if(ncellpoly>0)then
      call interp_vec_node2cell(vtemp,vecx,vecy)  !Map from nodes to cell-centers
    else
      call map_vec_full2active(vtemp,vecx,vecy) !Convert from full to active grid 
    endif    
      
    return
    end subroutine readvech5   
     
!************************************************************************
    subroutine readvecsteph5(afile,apath,itsind,thrs,vecx,vecy,ierr)
! Reads a scalar located in apath from afile XMDF file 
! written by Alex Sanchez, USACE-CHL
!************************************************************************
    use size_def,only: ncellsD,ncellsfull,ncellpoly
    use geo_def, only: idmap
    use interp_lib, only: interp_vec_node2cell
    use diag_lib
    use xmdf
    use prec_def
    implicit none
    !Input/Output
    character(len=*), intent(in) :: afile,apath
    integer, intent(in) :: itsind !Time step index
    real(ikind), intent(out) :: vecx(ncellsD),vecy(ncellsD),thrs
    integer, intent(out) :: ierr
    !Internal Variables
    integer :: fid,gid,ntimes,iloc
    real(8),allocatable :: timesd(:) !Output times
    real(4) :: vtemp(ncellsfull*2) !Must be single
    character(len=200) :: msg2,msg3,msg4,thepath

    call XF_OPEN_FILE(trim(afile),READONLY,fid,ierr)        
    if(ierr<0)then
      call diag_print_error('Could not open file: ',afile)
    endif
              
    thepath = trim(apath)
    iloc=index(thepath,'\')
    if (iloc.gt.0) thepath(iloc:iloc)='/'
    call XF_OPEN_GROUP(fid,trim(thepath),gid,ierr)
    if(ierr<0)then
      call XF_CLOSE_FILE(fid,ierr)
      !write(msg2,*,iostat=ierr) '  File: ',trim(afile)
      !write(msg3,*,iostat=ierr) '  Path: ',trim(thepath)
      !call diag_print_warning('Could not open dataset from',msg2,msg3)
      ierr = -2 !Could not open group
      return
    endif
               
    call XF_READ_VECTOR_VALUES_TIMESTEP(gid,itsind,ncellsfull,2,vtemp,ierr)
    if(ierr<0)then
      call XF_CLOSE_GROUP(gid,ierr)
      call XF_CLOSE_FILE(fid,ierr)
      !write(msg2,*,iostat=ierr) '  File: ',trim(afile)
      !write(msg3,*,iostat=ierr) '  Path: ',trim(thepath)
      !write(msg4,*,iostat=ierr) '  Time step: ',itsind
      !call diag_print_warning('Could not read vector time step values',msg2,msg3)
      ierr = 3 !Could not read timestep
      return
    endif
    
    call XF_GET_DATASET_NUM_TIMES(gid,ntimes,ierr)
    allocate(timesd(ntimes))
    call XF_GET_DATASET_TIMES(gid,ntimes,timesd,ierr)
    if(ierr<0)then
      thrs = -999.0
      write(msg2,*,iostat=ierr) '  File: ',trim(afile)
      write(msg3,*,iostat=ierr) '  Path: ',trim(thepath)
      write(msg4,*,iostat=ierr) '  Time step: ',itsind
      call diag_print_warning('Could not read time stamp from',msg2,msg3,msg4)
      ierr = 4
    else
      thrs = timesd(itsind)
    endif
    
    call XF_CLOSE_GROUP(gid,ierr)
    call XF_CLOSE_FILE(fid,ierr)
            
    if(ncellpoly>0)then
      call interp_vec_node2cell(vtemp,vecx,vecy)  !Map from nodes to cell-centers
    else
      call map_vec_full2active(vtemp,vecx,vecy) !Convert from full to active grid 
    endif
      
    return
    end subroutine readvecsteph5

!********************************************************************************
      subroutine readveclasth5(afile,apath,ntimes,vecx,vecy,reftimed,thrs,ierr)
! Reads a scalar located in apath from afile XMDF file 
! written by Alex Sanchez, USACE-CHL
!********************************************************************************
    use size_def, only: ncellsD,ncellsfull,ncellpoly
    use geo_def, only: idmap
    use interp_lib, only: interp_vec_node2cell
    use diag_def
    use diag_lib
    use prec_def
    use xmdf
    implicit none
    !Input/Output
    character(len=*), intent(in) :: afile,apath
    real(ikind), intent(out) :: vecx(ncellsD),vecy(ncellsD),thrs
    real(8), intent(out) :: reftimed            !Reference time
    integer, intent(out) :: ierr,ntimes
    !Internal Variables
    integer :: fid,gid,iloc
    real(8), allocatable :: timed(:)  !Output times
    real(4) :: vtemp(ncellsfull*2) !Must be single
    character(100) :: thepath

    call XF_OPEN_FILE(trim(afile),READONLY,fid,ierr)        
    if(ierr<0)then
      ierr = -1 !Could not open file
      call diag_print_error('Could not open file: ',afile)
    endif
      
    thepath = trim(apath)
    iloc=index(thepath,'\')
    if (iloc.gt.0) thepath(iloc:iloc)='/'
    call XF_OPEN_GROUP(fid,trim(thepath),gid,ierr)
    if(ierr<0)then
      call XF_CLOSE_FILE(fid,ierr)  
      !write(msg2,*,iostat=ierr) '  File: ',trim(afile)
      !write(msg3,*,iostat=ierr) '  Path: ',trim(thepath)
      !call diag_print_warning('Could not open dataset from',msg2,msg3)
      ierr = -2
      return
    endif
      
    call XF_GET_DATASET_NUM_TIMES(gid,ntimes,ierr)
    if(ierr<0)then
      call XF_CLOSE_GROUP(gid,ierr)
      call XF_CLOSE_FILE(fid,ierr)  
      write(msg2,*,iostat=ierr) '  File: ',trim(afile)
      write(msg3,*,iostat=ierr) '  Path: ',trim(thepath)
      call diag_print_error('Could not read number of times from ',msg2,msg3)
    endif
    
    call XF_READ_VECTOR_VALUES_TIMESTEP(gid,ntimes,ncellsfull,2,vtemp,ierr)
    if(ierr<0)then
      call XF_CLOSE_GROUP(gid,ierr)
      call XF_CLOSE_FILE(fid,ierr)  
      return
    endif
      
!    string = trim(apath) // 'TIME'
!    call XF_OPEN_GROUP(fid,trim(string),gid,ierr)
!    if(ierr<0)then
!      call diag_print_error('Invalid dataset path: ',thepath)
!    endif    
!    call XF_READ_SCALAR_VALUES_TIMESTEP(gid,1,1,thrs,ierr)
        
    allocate(timed(ntimes))
    call XF_GET_DATASET_TIMES(gid,ntimes,timed,ierr)   
    if(ierr<0)then
      thrs = -999.0
      write(msg2,*,iostat=ierr) '  File: ',trim(afile)
      write(msg3,*,iostat=ierr) '  Path: ',trim(thepath)
      write(msg4,*,iostat=ierr) '  Time step: ',ntimes
      call diag_print_warning('Could not read times from ',msg2,msg3,msg4)
      ierr = 4
    else
      thrs = timed(ntimes)
    endif
    deallocate(timed)
    
    call XF_GET_DATASET_REFTIME(gid,reftimed,ierr)
    if(ierr<0)then
      reftimed = -999.0d0
      write(msg2,*,iostat=ierr) '  File: ',trim(afile)
      write(msg3,*,iostat=ierr) '  Path: ',trim(thepath)
      call diag_print_warning('Could not read reference time from',msg2,msg3)
      ierr = 5
    endif
      
    call XF_CLOSE_GROUP(gid,ierr)
    call XF_CLOSE_FILE(fid,ierr)
            
    if(ncellpoly>0)then
      call interp_vec_node2cell(vtemp,vecx,vecy)  !Map from nodes to cell-centers
    else
      call map_vec_full2active(vtemp,vecx,vecy) !Convert from full to active grid 
    endif
      
    return
    end subroutine readveclasth5

!********************************************************************************
      subroutine readvectimeh5(afile,apath,thrs,vecx,vecy,reftimed,ierr)
! Reads a scalar located in apath from afile XMDF file 
! written by Alex Sanchez, USACE-CHL
!********************************************************************************
    use size_def, only: ncellsD,ncellsfull,ncellpoly
    use geo_def, only: idmap
    use interp_lib, only: interp_vec_node2cell
    use diag_def
    use diag_lib
    use prec_def
    use xmdf
    implicit none
    !Input
    character(len=*), intent(in) :: afile,apath
    real(ikind), intent(in) :: thrs
    !Output
    real(ikind), intent(out) :: vecx(ncellsD),vecy(ncellsD)
    real(8), intent(out) :: reftimed            !Reference time
    integer, intent(out) :: ierr
    !Internal Variables
    integer :: i,fid,gid,nstep,ntimes,iloc
    real(8), allocatable :: timed(:)  !Output times
    real(8) :: thrsd,terrd,terrdmin
    real(4) :: vtemp(ncellsfull*2) !Must be single
    character(100) :: thepath

    call XF_OPEN_FILE(trim(afile),READONLY,fid,ierr)        
    if(ierr<0)then
      ierr = -1 !Could not open file
      call diag_print_error('Could not open file: ',afile)
    endif
      
    thepath = trim(apath)
    iloc=index(thepath,'\')
    if (iloc.gt.0) thepath(iloc:iloc)='/'
    call XF_OPEN_GROUP(fid,trim(thepath),gid,ierr)
    if(ierr<0)then
      call XF_CLOSE_FILE(fid,ierr)  
      !write(msg2,*,iostat=ierr) '  File: ',trim(afile)
      !write(msg3,*,iostat=ierr) '  Path: ',trim(thepath)
      !call diag_print_warning('Could not open dataset from',msg2,msg3)
      ierr = -2
      return
    endif
      
    call XF_GET_DATASET_NUM_TIMES(gid,ntimes,ierr)
    if(ierr<0)then
      call XF_CLOSE_GROUP(gid,ierr)
      call XF_CLOSE_FILE(fid,ierr)  
      write(msg2,*,iostat=ierr) '  File: ',trim(afile)
      write(msg3,*,iostat=ierr) '  Path: ',trim(thepath)
      call diag_print_error('Could not read number of times from ',msg2,msg3)
    endif
      
!    string = trim(apath) // 'TIME'
!    call XF_OPEN_GROUP(fid,trim(string),gid,ierr)
!    if(ierr<0)then
!      call diag_print_error('Invalid dataset path: ',thepath)
!    endif    
!    call XF_READ_SCALAR_VALUES_TIMESTEP(gid,1,1,thrs,ierr)
        
    allocate(timed(ntimes))
    call XF_GET_DATASET_TIMES(gid,ntimes,timed,ierr)   
    if(ierr<0)then
      write(msg2,*,iostat=ierr) '  File: ',trim(afile)
      write(msg3,*,iostat=ierr) '  Path: ',trim(thepath)
      write(msg4,*,iostat=ierr) '  Time step: ',ntimes
      call diag_print_warning('Could not read times from ',msg2,msg3,msg4)
      ierr = 4
    else
      nstep = 1
      terrdmin = 1.0e9
      thrsd = dble(thrs)
      do i=1,ntimes
        terrd = abs(thrsd-timed(i))
        if(terrd<terrdmin)then
          nstep = i
          terrdmin = terrd
        endif
      enddo
      if(terrdmin>0.001)then
        write(msg2,*,iostat=ierr) '  File: ',trim(afile)
        write(msg3,*,iostat=ierr) '  Path: ',trim(thepath)
        write(msg4,*,iostat=ierr) '  Time: ',thrs,' hrs'
        call diag_print_warning('Could not find time: ',msg2,msg3,msg4)
      endif
    endif
    deallocate(timed)
    
    call XF_READ_VECTOR_VALUES_TIMESTEP(gid,nstep,ncellsfull,2,vtemp,ierr)
    if(ierr<0)then
      call XF_CLOSE_GROUP(gid,ierr)
      call XF_CLOSE_FILE(fid,ierr)  
      return
    endif
    
    call XF_GET_DATASET_REFTIME(gid,reftimed,ierr)
    if(ierr<0)then
      reftimed = -999.0d0
      write(msg2,*,iostat=ierr) '  File: ',trim(afile)
      write(msg3,*,iostat=ierr) '  Path: ',trim(thepath)
      call diag_print_warning('Could not read reference time from',msg2,msg3)
      ierr = 5
    endif
      
    call XF_CLOSE_GROUP(gid,ierr)
    call XF_CLOSE_FILE(fid,ierr)
            
    if(ncellpoly>0)then
      call interp_vec_node2cell(vtemp,vecx,vecy)  !Map from nodes to cell-centers
    else
      call map_vec_full2active(vtemp,vecx,vecy) !Convert from full to active grid 
    endif
      
    return
    end subroutine readvectimeh5
#endif     

end module