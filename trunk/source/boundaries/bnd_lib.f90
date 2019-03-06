!===============================================================================
module bnd_lib
! Boundary Condition Library
!
! Contains the following routines
!  Boundary Strings:
!    read_bndstr - Reads a cellstring from an ASCII or XMDF file,
!                  alllocates and initializes and cell id array
!  Input Formats:
!    read_fluxdata - Reads the boundary flux data from the appropriate file
!    read_snglwsedata - Reads a single water level curve from the from the appropriate file   
!    read_multiwsedata - Reads multiple water level curves from the appropriate file
!    read_multiwseh5   -  Reads the multiple water curve data from an XMDF file
!    read_multiveldata - Reads multiple current velocity curves from the appropriate file
!    read_multivelh5 - Reads multiple current velocity curves from an XMDF file
!
! written by Alex Sanchez, USACE-CHL
!===============================================================================
#include "CMS_cpp.h"
    implicit none
    
contains

!******************************************************************************************
    subroutine read_bndstr(bidfile,bidpath,iBndStrType,idbnd,nstrcells,icells,kfaces)
! Reads a cellstring from an ASCII or XMDF file,
! alllocates and initializes and cell id array
! written by Alex Sanchez, USACE-CHL
!******************************************************************************************
#include "CMS_cpp.h"
    use diag_lib, only: diag_print_error
    
    implicit none
    !Input/Output
    character(len=*),intent(in) :: bidfile             !Boundary ID file
    character(len=*),intent(in) :: bidpath             !Boundary ID path 
    integer,         intent(in) :: iBndStrType         !Boundary string type, 1-cellstring, 2-nodestring
    integer,         intent(inout) :: idbnd            !Boundary ID number (for *.2dm and *.bid files)
    integer,         intent(out):: nstrcells           !Number of cells in boundary string
    integer,         intent(inout),pointer:: icells(:) !Boundary cells
    integer,         intent(inout),pointer:: kfaces(:) !Boundary faces  
    !Internal Variables
    integer :: nstrnodes
    character*10 :: aext
    integer, pointer :: jnodes(:)
    
    call fileext(bidfile,aext)
    selectcase(aext)
    case('h5')        
      !Only cellstrings stored in h5 file
      if(iBndStrType==2)then
        call diag_print_error('Only cellstrings allowed in *.h5 files')
      endif
#ifdef XMDF_IO
      call read_cellstr_h5(bidfile,bidpath,nstrcells,icells) 
#endif
    case('bid')
      if(idbnd==0) read(bidpath,*) idbnd  
      if(iBndStrType==1)then
        call bnd_read_str_bid(bidfile,idbnd,nstrcells,icells)
      else
        call bnd_read_str_bid(bidfile,idbnd,nstrnodes,jnodes)       
      endif
    case('2dm')
      !Only nodestrings stored in m2d file  
      if(iBndStrType==1)then
        call diag_print_error('Only nodestrings allowed in *.2dm files')
      endif
      if(idbnd==0) read(bidpath,*) idbnd
      call bnd_read_nodestr_2dm(idbnd,nstrnodes,jnodes) 
    endselect
    
    if(iBndStrType==1)then
      !Map id's from full to active grid
      call map_cell_full2active(nstrcells,icells)  
      !Check cell string corners and repeat boundary cells if necessary  
      call cellstr_bnd(nstrcells,icells,kfaces) 
    elseif(iBndStrType==2)then
      !Convert nodestring to cellstring and determine cell faces
      call nodestr_bnd(nstrnodes,jnodes,nstrcells,icells,kfaces)  
    else
      call diag_print_error('Invalid boundary string type')
    endif
    
    return
    endsubroutine read_bndstr    
    
#ifdef XMDF_IO
!*****************************************************************************
    subroutine read_cellstr_h5(bidfile,bidpath,nstrcells,icells)
!*****************************************************************************
    use size_def, only: ncells    
    use diag_lib, only: diag_print_error
    use xmdf
    
    implicit none
    !Input/Output
    character(len=*),intent(in)  :: bidfile,bidpath
    integer,         intent(out) :: nstrcells
    integer,         intent(inout),pointer :: icells(:)    
    !Internal variables
    integer, allocatable :: idtemp(:)
    integer :: pid,gid,ierr    
    character(len=200) :: str1,str2
    
    allocate(idtemp(10000))
    call XF_OPEN_FILE(bidfile,READONLY,pid,ierr)
    if(ierr<0)then
      call diag_print_error('Could not open XMDF Boundary ID File: ',&
        bidfile,'  Check the input and rerun')
    endif
    call XF_OPEN_GROUP(pid,bidpath,gid,ierr)
    if(ierr<0)then
      str1 = '  File ' // trim(bidfile)
      str2 = '  Path ' // trim(bidpath)  
      call diag_print_error('Incorrect Boundary ID Path',&
        str1,str2,'  Check the input and rerun')
    endif
    call XF_GET_PROPERTY_NUMBER(gid,'Cells',nstrcells,ierr)
    call XF_READ_PROPERTY_INT(gid,'Cells',nstrcells,idtemp(1),ierr)    
    call XF_CLOSE_GROUP(gid,ierr)
    call XF_CLOSE_FILE(pid,ierr)
    allocate(icells(nstrcells))
    icells = idtemp(1:nstrcells)
    deallocate(idtemp)
        
    return
    endsubroutine read_cellstr_h5
#endif

!*****************************************************************************
    subroutine bnd_read_str_bid(bidfile,idbnd,nstrelem,ielem)
! Reads a boundary string ID information from a Boundary ID file (*.bid)
! The string in the BID file may be a cellstring or a nodestring
! written by Alex Sanchez, USACE-CHL    
!*****************************************************************************
    use diag_lib, only: diag_print_error
    
    implicit none
    !Input/Output
    character(len=*),intent(in)  :: bidfile !Boundary ID file
    integer, intent(in)  :: idbnd    !Boundary ID within bidFile
    integer, intent(out) :: nstrelem !# of boundary string elements (cells or nodes)
    integer, intent(inout),pointer :: ielem(:) !Boundary string element (cell or node) ID's
    !Internal Variables
    integer :: i,k,nbndstr,idum
    logical :: foundfile
    
    inquire(file=bidfile,exist=foundfile)
    if(.not.foundfile)then
      call diag_print_error('Could not find the Boundary ID file: ',bidfile)
    endif
    open(38,file=bidfile,status='old')
    read(38,*) nbndstr       !# of boundary strings
    if(idbnd>nbndstr)then                                   !Removed this check because each type of boundary has it's own BID file now.  MEB 02/21/2018
      call diag_print_error('Boundary ID number is larger than the number of ',&
        '  strings in the Boundary ID File: ',bidfile)
    endif
    do k=1,nbndstr
      read(38,*) nstrelem  
      if(k/=idbnd)then                                      !No need to skip elements at present time.  MEB 02/21/2018 
        read(38,*) (idum,i=1,nstrelem) !Skip elements
      else
        allocate(ielem(nstrelem)) 
        read(38,*) (ielem(i),i=1,nstrelem)
        close(38)
        return
      endif
    enddo
    call diag_print_error('Could not read Boundary ID File: ',bidfile)
    
    return
    endsubroutine bnd_read_str_bid
    
!*****************************************************************************
    subroutine bnd_read_nodestr_2dm(idbnd,nstrnodes,jnodes)   
! Extracts the boundary cellstring from the bc_str variable.
!
! Because the node strings are saved at the end of the file in a weird format
! it is easier to read them once, store them, and then reference them later
!
! Note: Not all of the nodes strings in the 2dm file are necessarily forcing
! strings and they may not be necessarilly all of the strings either.
!
! written by Alex Sanchez, USACE-CHL
!*****************************************************************************
    use bnd_def, only: bc_str,nbcstr
    use diag_lib, only: diag_print_error
    use prec_def, only: ikind
    
    implicit none
    !Input/Output
    integer,intent(in)  :: idbnd
    integer,intent(out) :: nstrnodes
    integer,intent(inout),pointer :: jnodes(:)

    !Initialize
    if(idbnd>nbcstr)then
      call diag_print_error('Boundary ID does not exist')
    endif
    nstrnodes = bc_str(idbnd)%nnodes
    !if(allocated(jnodes)) deallocate(jnodes)
    allocate(jnodes(nstrnodes))
    jnodes=bc_str(idbnd)%nodes
    
    return
    endsubroutine bnd_read_nodestr_2dm
    
!*************************************************************
    subroutine cellstr_bnd(nstrcells,icells,kfaces)
! Determines the cellstring boundary faces and if necessary
! repeats corner cells to account for multiple boundary faces.
! Only used for Cartesian grids. Boundaries for unstructured
! grids are specified using node strings which avoids confusion
! in determining the boundary faces.
!
! Change Log:
!  Bug fix (Alex Sanchez 5/17/14):  Corrected problem with 
!    telescoping refinement at cell second and second to 
!    last cells. 
!
! Author: Alex Sanchez, USACE-CHL
!*************************************************************
    use size_def, only: ncells,ncellsD
    use geo_def, only: ncface,cell2cell,idirface   
    use diag_lib, only: diag_print_error
    use bnd_def, only: nbndcorner, bndcorner
    
    implicit none
    !Input/Output
    integer,intent(inout)         :: nstrcells
    integer,intent(inout),pointer :: icells(:)
    integer,intent(inout),pointer :: kfaces(:)
    !Internal Variables
    integer :: i,j,k,id1,id2,id,kk
    integer :: nbndtemp,ibndtemp(ncellsD),kfacetemp(ncellsD)

    nbndtemp = 0    !Counter for cell-string, cells may be repeated if at corners   
            
    !Use second cell in string to determine outward face of first cell in string
    id1 = icells(1)
    id2 = icells(2)
    do k=1,ncface(id2) !Determine the outward face (boundary)
      if(cell2cell(k,id2)>ncells)then
        nbndtemp = nbndtemp + 1
	    ibndtemp(nbndtemp) = id1
        do kk=1,ncface(id1)
          if(idirface(kk,id1)==idirface(k,id2))then
            kfacetemp(nbndtemp) = kk
            exit
          endif
        enddo
	  endif   
    enddo
	    
    do j=2,nstrcells-1   !Wu
      id = icells(j)
      do k=1,ncface(id) !Determine the outward face (boundary)
        if(cell2cell(k,id)>ncells)then 
	      nbndtemp = nbndtemp + 1
	      ibndtemp(nbndtemp) = id       
	      kfacetemp(nbndtemp) = k 	          
	    endif
      enddo !k        
    enddo !j   
	    
    !Use second-to-last cell to determine outward face of last cell in string
    id1 = icells(nstrcells)
    id2 = icells(nstrcells-1)
    do k=1,ncface(id2) !Determine the outward face (boundary)
      if(cell2cell(k,id2)>ncells)then
        nbndtemp = nbndtemp + 1
	    ibndtemp(nbndtemp) = id1
        do kk=1,ncface(id1)
          if(idirface(kk,id1)==idirface(k,id2))then
            kfacetemp(nbndtemp) = kk
            exit
          endif
        enddo
      endif
    enddo !k-face
    
    if(nbndtemp==0) then
        !there is an error in the cellstring specification somewhere
        call diag_print_error('Cellstring has no cells at boundary edge')
    endif
    
    
    !Note: The variable 'bndcorner' is never used except here.  Do we even need all the extra code here for it?
    if(nstrcells/=nbndtemp)then 
      nbndcorner = nbndtemp-nstrcells
      if (allocated(bndcorner)) deallocate (bndcorner)    !meb 03/06/2019  If another WSE cellstring has a corner, you have to deallocate first.
      allocate(bndcorner(nbndcorner))
      
      deallocate(icells)
      nstrcells = nbndtemp
      allocate(icells(nstrcells))      
      icells(1:nstrcells) =  ibndtemp(1:nstrcells)      
      j=1
      do i=2,nstrcells
        if(icells(i-1)==icells(i)) then
          bndcorner(j)=icells(i)
          j=j+1
        endif
      enddo
    endif
    allocate(kfaces(nstrcells))
    kfaces(1:nstrcells) = kfacetemp(1:nstrcells)
    
    return
    endsubroutine cellstr_bnd 
    
!*************************************************************
    subroutine cellstr_bnd_old(nstrcells,icells,kfaces)
!*************************************************************
    use size_def, only: ncells,ncellsD
    use geo_def, only: ncface,cell2cell    
    
    implicit none
    !Input/Output
    integer,intent(inout)         :: nstrcells
    integer,intent(inout),pointer :: icells(:)
    integer,intent(inout),pointer :: kfaces(:)
    !Internal Variables
    integer :: j,k,id1,id2,id
    integer :: nbndtemp,ibndtemp(ncellsD),kfacetemp(ncellsD)

    nbndtemp = 0    !Counter for cell-string, cells may be repeated if at corners   
            
    !Use second cell in string to determine outward face of first cell in string
    id1 = icells(1)
    id2 = icells(2)
    do k=1,ncface(id2) !Determine the outward face (boundary)
      if(cell2cell(k,id2)>ncells)then
        nbndtemp = nbndtemp + 1
	    ibndtemp(nbndtemp) = id1
	    kfacetemp(nbndtemp) = k
	  endif   
    enddo
	    
    do j=2,nstrcells-1   !   !Wu
      id = icells(j)
      do k=1,ncface(id) !Determine the outward face (boundary)
        if(cell2cell(k,id)>ncells)then 
	      nbndtemp = nbndtemp + 1
	      ibndtemp(nbndtemp) = id       
	      kfacetemp(nbndtemp) = k 	          
	    endif
      enddo !k        
    enddo !j   
	    
    !Use second-to-last cell to determine outward face of last cell in string
    id1 = icells(nstrcells)
    id2 = icells(nstrcells-1)
    do k=1,ncface(id2) !Determine the outward face (boundary)
      if(cell2cell(k,id2)>ncells)then
        nbndtemp = nbndtemp + 1
	    ibndtemp(nbndtemp) = id1
	    kfacetemp(nbndtemp) = k
      endif
    enddo !k-face
    
    if(nstrcells/=nbndtemp)then 
      deallocate(icells)
      nstrcells = nbndtemp
      allocate(icells(nstrcells))      
      icells(1:nstrcells) =  ibndtemp(1:nstrcells)      
    endif
    allocate(kfaces(nstrcells))
    kfaces(1:nstrcells) = kfacetemp(1:nstrcells)
    
    return
    endsubroutine cellstr_bnd_old
    
!*****************************************************************************
    subroutine nodestr_bnd(nstrnodes,jnodes,nstrcells,icells,kfaces)   
! Determine boundary cells and faces from node strings
! written by Alex Sanchez, USACE-CHL
!*****************************************************************************
    use size_def, only: ncells
    use geo_def, only: nncell,node2cell,ncface,cell2node
    use bnd_def, only: bc_str,nbcstr
    use prec_def, only: ikind
    
    implicit none
    !Input/Output
    integer,intent(in) :: nstrnodes           !Number of nodes on boundary string
    integer,intent(in) :: jnodes(nstrnodes)   !Node ID's on boundary string
    integer,intent(out):: nstrcells           !Number of cells on boundary string
    integer,intent(inout),pointer:: icells(:) !Boundary cells
    integer,intent(inout),pointer:: kfaces(:) !Boundary faces
    !Internal variables
    integer :: ii,j,k,jj,jj2,nd1,nd2

    nstrcells = nstrnodes - 1
    !if(allocated(icells)) deallocate(icells)
    !if(allocated(kfaces)) deallocate(kfaces)
    allocate(icells(nstrcells),kfaces(nstrcells))
    
    icells=0 
    do j=1,nstrnodes-1
      nd1=jnodes(j)
      nd2=jnodes(j+1)
      !Search neighboring cells for one connected to nd1 and nd2
 dok: do k=1,nncell(nd1)
        ii=node2cell(k,nd1)
        if(ii>ncells) cycle
        do jj=1,ncface(ii)
          if(jj==ncface(ii))then
            jj2=1
          else
            jj2=jj+1  
          endif
          if((cell2node(jj,ii)==nd1 .and. cell2node(jj2,ii)==nd2) .or. &
             (cell2node(jj,ii)==nd2 .and. cell2node(jj2,ii)==nd1))then
            icells(j)=ii
            kfaces(j)=jj
            exit dok
          endif
        enddo !jj neighboring node
      enddo dok !k neighboring cell
    enddo !j node pair
    
    return
    endsubroutine nodestr_bnd

!*************************************************************************
    subroutine read_fluxdata(datfile,datpath,ntimes,times,fluxes)
! Reads the boundary flux data from the appropriate file
! written by Alex Sanchez, USACE-CHL
!*****************************************************************************    
#include "CMS_cpp.h"
    use prec_def, only: ikind
    use diag_lib
    use comvarbl, only: tjulday0
    use in_lib, only: read_xys,read_tsd

#ifdef XMDF_IO    
    use in_xmdf_lib, only: read_dataseth5
#endif
    implicit none

    !Input/Output
    integer,intent(out):: ntimes
    real(ikind),intent(inout),pointer:: times(:)
    real(ikind),intent(inout),pointer:: fluxes(:)
    character(len=*),intent(in) :: datfile,datpath
    !Internal variables
    character(len=10) :: aext
    character(len=200) :: aname,atype
    integer :: ndat
    real(ikind) :: tjuldaybegflux
    real(ikind), pointer :: dat(:,:)
    
    call fileext(datfile,aext)      
    if(aext(1:2)=='h5')then
#ifdef XMDF_IO
      call read_dataseth5(datfile,datpath,'Times',ntimes,times)
      call read_dataseth5(datfile,datpath,'Flow',ntimes,fluxes)
#else
      call diag_print_error('Cannot read flux time series from *.h5 file without XMDF libraries')
#endif
    elseif(aext(1:3)=='xys')then
      call read_xys(datfile,ntimes,times,fluxes)
    elseif(aext(1:3)=='tsd')then 
      call read_tsd(datfile,aname,atype,ndat,ntimes,tjuldaybegflux,times,dat)
      allocate(fluxes(ntimes))
      fluxes(:) = dat(:,1)
      deallocate(dat)
      times = times/3600.0 + 24.0*(tjuldaybegflux - tjulday0) !Note conversion to hours
    endif 

    return
    endsubroutine read_fluxdata

!*****************************************************************************
    subroutine read_offsetwsedata(datfile,datpath,ntimes,times,wseoffset)
! Reads a wse offset curve
! written by Honghai Li, USACE-CHL (01/19/17)
!*****************************************************************************    
#include "CMS_cpp.h"
    use in_lib, only: read_xys,read_tsd
    use comvarbl, only: tjulday0
    use diag_lib, only: diag_print_error
    use prec_def, only: ikind
#ifdef XMDF_IO    
    use in_xmdf_lib, only: read_dataseth5
#endif
    
    implicit none
    !Input/Output
    integer,intent(out):: ntimes
    real(ikind),intent(inout),pointer:: times(:)
    real(ikind),intent(inout),pointer:: wseoffset(:)
    character(len=*),intent(in) :: datfile,datpath
    !Internal Variables
    integer :: ndat
    character(len=10) :: aext
    character(len=100) :: aname,atype
    real(ikind) :: tjuldaybegwse
    real(ikind), pointer :: dat(:,:)

    call fileext(datfile,aext)
    if(aext(1:2)=='h5')then
#ifdef XMDF_IO
      call read_dataseth5(datfile,datpath,'Offset_Times',ntimes,times)
      call read_dataseth5(datfile,datpath,'Offset',ntimes,wseoffset)
#else
      call diag_print_error('Cannot read wse offset time series from *.h5 file without XMDF libraries')
#endif
    elseif(aext(1:3)=='xys')then
      call read_xys(datfile,ntimes,times,wseoffset)
    elseif(aext(1:3)=='tsd')then
      call read_tsd(datfile,aname,atype,ndat,ntimes,tjuldaybegwse,times,dat)
      allocate(wseoffset(ntimes))
      wseoffset = dat(:,1)
      deallocate(dat)
!      write(2000,*)'times = ',times,'tjuldaybegwse =',tjuldaybegwse,'tjulday0 =',tjulday0 
      times = times/3600.0 + 24.0*(tjuldaybegwse - tjulday0) !Note conversions to hours
!      write(2000,*)'times = ',times
    endif
    
    return
    endsubroutine read_offsetwsedata

!*****************************************************************************
    subroutine read_snglwsedata(datfile,datpath,ntimes,times,wse)
! Reads a single water level boundary condition
! written by Alex Sanchez, USACE-CHL
!*****************************************************************************    
#include "CMS_cpp.h"
    use in_lib,   only: read_xys,read_tsd
    use comvarbl, only: tjulday0
    use diag_lib, only: diag_print_error
    use prec_def, only: ikind
#ifdef XMDF_IO
    use in_xmdf_lib,   only: read_dataseth5
#endif
    
    implicit none
    !Input/Output
    integer,intent(out):: ntimes
    real(ikind),intent(inout),pointer:: times(:)
    real(ikind),intent(inout),pointer:: wse(:)
    character(len=*),intent(in) :: datfile,datpath
    !Internal Variables
    integer :: ndat
    character(len=10) :: aext
    character(len=100) :: aname,atype
    real(ikind) :: tjuldaybegwse
    real(ikind), pointer :: dat(:,:)

    call fileext(datfile,aext)
    if(aext(1:2)=='h5')then
#ifdef XMDF_IO
      call read_dataseth5(datfile,datpath,'Times',ntimes,times)
      call read_dataseth5(datfile,datpath,'WaterLevel',ntimes,wse)
#else
      call diag_print_error('Cannot read single wse time series from *.h5 file without XMDF libraries')
#endif
    elseif(aext(1:3)=='xys')then
      call read_xys(datfile,ntimes,times,wse)
    elseif(aext(1:3)=='tsd')then
      call read_tsd(datfile,aname,atype,ndat,ntimes,tjuldaybegwse,times,dat)
      allocate(wse(ntimes))
      wse = dat(:,1)
      deallocate(dat)
      times = times/3600.0 + 24.0*(tjuldaybegwse - tjulday0) !Note conversions to hours
    endif
    
    return
    endsubroutine read_snglwsedata

!***************************************************************************************
    subroutine read_multiwsedata(datfile,datpath,nstrcells,icells,ntimes,times,wsedata)
! Reads a multiple wse boundary condition data from an input file
! written by Alex Sanchez, USACE-CHL     
!***************************************************************************************
#include "CMS_cpp.h"    
    use diag_lib, only: diag_print_error
    use comvarbl, only: tjulday0
    use in_lib, only: read_tsd
    use prec_def, only: ikind
    
    implicit none
    !Input/Output
    integer,intent(in):: nstrcells
    integer,intent(inout),pointer:: icells(:)
    integer,intent(out):: ntimes
    real(ikind),intent(inout),pointer:: times(:)
    real(ikind),intent(inout),pointer:: wsedata(:,:)
    !Internal Variables
    integer :: ndat
    real(ikind) :: tjuldaybegwse
    character(len=*),intent(in) :: datfile,datpath
    character(len=100) :: aname,atype
    character(len=10) :: aext

    call fileext(datfile,aext)
    if(aext(1:2)=='h5')then
#ifdef XMDF_IO
      call read_multiwseh5(datfile,datpath,nstrcells,icells,ntimes,times,wsedata)
#else
      call diag_print_error('Cannot read muliple wse ',&
        '  time series from an *.h5 file without XMDF libraries')
#endif
    elseif(aext(1:3)=='tsd')then
      call read_tsd(datfile,aname,atype,ndat,ntimes,tjuldaybegwse,times,wsedata)
      if(nstrcells/=ndat)then
        call diag_print_error('Problem reading multiple wse data',&
          '  Number of time-series curves must be equal to the number of boundary cells')
      endif
      times = times/3600.0 + 24.0*(tjuldaybegwse - tjulday0) !Note conversions to hours
    endif
    
    return
    endsubroutine read_multiwsedata

!********************************************************************************************
    subroutine read_multiveldata(datfile,datpath,nstrcells,icells,ntimes,times,udata,vdata)
! Reads a multiple current velocity boundary condition data from an input file
! written by Alex Sanchez, USACE-CHL     
!********************************************************************************************
#include "CMS_cpp.h"
    use diag_lib, only: diag_print_error
    use comvarbl, only: tjulday0
    use in_lib, only: read_tsd
    use prec_def, only: ikind
    
    implicit none
    !Input/Output
    integer,intent(in):: nstrcells
    integer,intent(inout),pointer:: icells(:)
    integer,intent(out):: ntimes
    real(ikind),intent(inout),pointer:: times(:)
    real(ikind),intent(inout),pointer:: udata(:,:),vdata(:,:) !(time,cell)
    character(len=*),intent(in) :: datfile,datpath
    !Internal Variables
    integer :: i,j,j2,j2m1
    integer :: ndat,ncurves
    real(ikind) :: tjuldaybegvel
    character(len=10) :: aext
    character(len=100) :: aname,atype
    real(ikind), pointer :: veldata(:,:)
    
    call fileext(datfile,aext)
    if(aext(1:2)=='h5')then
#ifdef XMDF_IO
      call read_multivelh5(datfile,datpath,nstrcells,icells,ntimes,times,udata,vdata)
#else
      call diag_print_error('Cannot read muliple velocity ',&
        '  time series from an *.h5 file without XMDF libraries')
#endif
    elseif(aext(1:3)=='tsd')then  
      call read_tsd(datfile,aname,atype,ndat,ntimes,tjuldaybegvel,times,veldata)
      ncurves = ndat/2
      if(nstrcells/=ncurves)then
        call diag_print_error('Problem reading multiple velocity data',&
          '  Number of time-series curves must be equal to the number of boundary cells')
      endif
      allocate(udata(ntimes,nstrcells),vdata(ntimes,nstrcells))
      do i=1,ntimes
        do j=1,ncurves
          j2 = j*2
          j2m1  = j2 - 1       
          udata(i,j) = veldata(i,j2m1)
          vdata(i,j) = veldata(i,j2)
        enddo
      enddo
      deallocate(veldata)
      times = times/3600.0 + 24.0*(tjuldaybegvel - tjulday0) !Note conversions to hours
    endif    
      
    return
    endsubroutine read_multiveldata
    
#ifdef XMDF_IO
!**********************************************************************************
    subroutine read_multiwseh5(datfile,datpath,nstrcells,icells,ntimes,times,wsedata)
!**********************************************************************************
    use geo_def, only: mapid
    use xmdf
    use prec_def, only: ikind
    use comvarbl, only: input_ver
    
    implicit none
    !Input/Output
    integer,intent(in):: nstrcells
    integer,intent(inout),pointer:: icells(:)
    integer,intent(out):: ntimes
    real(ikind),intent(inout),pointer:: times(:)    
    real(ikind),intent(inout),pointer:: wsedata(:,:) !(time,cell)
    character(len=*),intent(in) :: datfile,datpath
    !Internal Variables
    character :: astring*50,card*80, time_dset*10
    integer :: j,error,PID,GID,WID
    real(4),allocatable :: ftemp(:) !Should be single precision
    
    time_dset = 'Times'
    if (input_ver >= 5.0) time_dset = 'WSE_Times'
    call XF_OPEN_FILE(datfile,READONLY,PID,error)
    call XF_OPEN_GROUP(PID,datpath,GID,error)
    call XF_GET_PROPERTY_NUMBER(GID,time_dset,ntimes,error)             !call XF_GET_PROPERTY_NUMBER(GID,'Times',ntimes,error)           !MEB change from 'Times' 6/28/2016
    allocate(ftemp(ntimes),times(ntimes),wsedata(ntimes,nstrcells))
    call XF_READ_PROPERTY_FLOAT(GID,time_dset,ntimes,ftemp(1),error)    !call XF_READ_PROPERTY_FLOAT(GID,'Times',ntimes,ftemp(1),error)  !MEB change from 'Times' 6/28/2016
    times = ftemp
    call XF_OPEN_GROUP(GID,"WSE",WID,error)
    do j=1,nstrcells
      write(unit=astring,fmt='(I0)') mapid(icells(j))
      card = 'WaterLevel_'//trim(astring)
      call XF_READ_PROPERTY_FLOAT(WID,trim(card),ntimes,ftemp(1),error)
      if(error<0)then
        call error_invalid_dataset(datfile,datpath,card)            
      endif
      wsedata(:,j) = ftemp(:)
    enddo
    call XF_CLOSE_GROUP(GID,error)
    call XF_CLOSE_FILE(PID,error)
    deallocate(ftemp)

    return
    endsubroutine read_multiwseh5
 
!************************************************************************************************
    subroutine read_multivelh5(datfile,datpath,nstrcells,icells,ntimes,times,udata,vdata)
!************************************************************************************************
    use geo_def, only: mapid
    use diag_lib, only: diag_print_error
    use xmdf
    use prec_def, only: ikind
    use comvarbl, only: input_ver
    
    implicit none
    !Input/Output
    integer,intent(in):: nstrcells
    integer,intent(inout),pointer:: icells(:)
    integer,intent(out):: ntimes
    real(ikind),intent(inout),pointer:: times(:)    
    real(ikind),intent(inout),pointer:: udata(:,:),vdata(:,:)
    character(len=*),intent(in) :: datfile,datpath
    !Internal Variables
    character :: astring*50,card*80, time_dset*10
    integer :: j,error,PID,GID,LID,RID,TID,BID
    real(4),allocatable :: ftemp(:),ftemp2(:) !Should be single precision
      
    time_dset='Times'
    if(input_ver >= 5.0) time_dset='Vel_Times'
    
    call XF_OPEN_FILE(datfile,READONLY,PID,error)
    call XF_OPEN_GROUP(PID,datpath,GID,error)
    call XF_GET_PROPERTY_NUMBER(GID,time_dset,ntimes,error)            !call XF_GET_PROPERTY_NUMBER(GID,'Times',ntimes,error)           !MEB change from 'Times' 6/28/2016
    allocate(times(ntimes),udata(ntimes,nstrcells),vdata(ntimes,nstrcells))
    allocate(ftemp(ntimes),ftemp2(ntimes))
    call XF_READ_PROPERTY_FLOAT(GID,time_dset,ntimes,ftemp(1),error)   !call XF_READ_PROPERTY_FLOAT(GID,'Times',ntimes,ftemp(1),error)  !MEB change from 'Times' 6/28/2016
    times = ftemp
    call XF_OPEN_GROUP(GID,"Left",LID,error)
    call XF_OPEN_GROUP(GID,"Right",RID,error)
    call XF_OPEN_GROUP(GID,"Top",TID,error)
    call XF_OPEN_GROUP(GID,"Bottom",BID,error)
    do j=1,nstrcells
      write(unit=astring,fmt='(I0)') mapid(icells(j))
      card = 'Velocity_'//trim(astring)
      call XF_READ_PROPERTY_FLOAT(LID,trim(card),ntimes,ftemp(1),error) !Left face
      if(error<0)then
        call error_invalid_dataset(datfile,datpath,card)
      endif
      if(error<0) call diag_print_error('Problem reading Multiple WSE-Vel BC')
      udata(:,j) = ftemp(:)
      call XF_READ_PROPERTY_FLOAT(RID,trim(card),ntimes,ftemp(1),error) !Right face      
      if(error<0) call diag_print_error('Problem reading Multiple WSE-Vel BC')         
      call XF_READ_PROPERTY_FLOAT(TID,trim(card),ntimes,ftemp2(1),error) !Top face
      if(error<0) call diag_print_error('Problem reading Multiple WSE-Vel BC')
      vdata(:,j) = ftemp2(:)
      call XF_READ_PROPERTY_FLOAT(BID,trim(card),ntimes,ftemp2(1),error) !Bot face  
      if(error<0) call diag_print_error('Problem reading Multiple WSE-Vel BC')
                           
      !Average Left and Right faces
      udata(:,j)=0.5*(ftemp(:)+udata(:,j))            
      !Average Top and Bottom faces
      vdata(:,j)=0.5*(ftemp2(:)+vdata(:,j))     
    enddo !j     
    call XF_CLOSE_GROUP(LID,error)
    call XF_CLOSE_GROUP(RID,error)
    call XF_CLOSE_GROUP(TID,error)
    call XF_CLOSE_GROUP(BID,error)
    call XF_CLOSE_GROUP(GID,error)
    call XF_CLOSE_FILE(PID,error)
    deallocate(ftemp,ftemp2)
     
    return
    endsubroutine read_multivelh5
  
!************************************************************************************************
    subroutine read_multiwsevelh5(datfile,datpath,nstrcells,icells,ntimes,times,wsevel)
! Reads multiple water level and velocity data from an XMDF file
! The values are specified at each boundary cell
! written by Alex Sanchez, USACE-CHL
!************************************************************************************************
    use geo_def, only: mapid
    use diag_lib, only: diag_print_error
    use xmdf
    use prec_def, only: ikind
    use comvarbl, only: input_ver
    
    implicit none
    !Input/Output
    integer,intent(in):: nstrcells
    integer,intent(inout),pointer:: icells(:)
    integer,intent(out):: ntimes
    real(ikind),intent(inout),pointer:: times(:)    
    real(ikind),intent(inout),pointer:: wsevel(:,:,:)
    character(len=*),intent(in) :: datfile,datpath
    !Internal Variables
    character :: astring*50,card*80, time_dset*10
    integer :: j,error,PID,GID,LID,RID,TID,BID,WID
    real(4),allocatable :: ftemp(:),ftemp2(:) !Should be single precision
      
    time_dset='Times'
    if (input_ver >= 5.0) time_dset='WSE_Times'
    
    call XF_OPEN_FILE(datfile,READONLY,PID,error)
    call XF_OPEN_GROUP(PID,datpath,GID,error)
    call XF_GET_PROPERTY_NUMBER(GID,time_dset,ntimes,error)            !call XF_GET_PROPERTY_NUMBER(GID,'Times',ntimes,error)            !MEB change from 'Times' 6/28/2016
    allocate(times(ntimes),wsevel(ntimes,nstrcells,3))
    allocate(ftemp(ntimes),ftemp2(ntimes))
    call XF_READ_PROPERTY_FLOAT(GID,time_dset,ntimes,ftemp(1),error)   !call XF_READ_PROPERTY_FLOAT(GID,'Times',ntimes,ftemp(1),error)   !MEB change from 'Times' 6/28/2016
    times = ftemp
    call XF_OPEN_GROUP(GID,"Left",LID,error)
    call XF_OPEN_GROUP(GID,"Right",RID,error)
    call XF_OPEN_GROUP(GID,"Top",TID,error)
    call XF_OPEN_GROUP(GID,"Bottom",BID,error) 
    call XF_OPEN_GROUP(GID,"WSE",WID,error)
    do j=1,nstrcells
      write(unit=astring,fmt='(I0)') mapid(icells(j))
      card = 'Velocity_'//trim(astring)
      call XF_READ_PROPERTY_FLOAT(LID,trim(card),ntimes,ftemp(1),error) !Left face
      if(error<0)then
        call error_invalid_dataset(datfile,datpath,card)
      endif
      if(error<0) call diag_print_error('Problem reading Multiple WSE-Vel BC')
      wsevel(:,j,1) = ftemp(:)
      call XF_READ_PROPERTY_FLOAT(RID,trim(card),ntimes,ftemp(1),error) !Right face      
      if(error<0) call diag_print_error('Problem reading Multiple WSE-Vel BC')         
      call XF_READ_PROPERTY_FLOAT(TID,trim(card),ntimes,ftemp2(1),error) !Top face
      if(error<0) call diag_print_error('Problem reading Multiple WSE-Vel BC')
      wsevel(:,j,2) = ftemp2(:)
      call XF_READ_PROPERTY_FLOAT(BID,trim(card),ntimes,ftemp2(1),error) !Bot face  
      if(error<0) call diag_print_error('Problem reading Multiple WSE-Vel BC')
                           
      !Average Left and Right faces
      wsevel(:,j,1)=0.5*(ftemp(:)+wsevel(:,j,1))            
      !Average Top and Bottom faces
      wsevel(:,j,2)=0.5*(ftemp2(:)+wsevel(:,j,2))         
                   
      card = 'WaterLevel_'//trim(astring) 
      call XF_READ_PROPERTY_FLOAT(WID,trim(card),ntimes,ftemp(1),error) !WSE
      wsevel(:,j,3) = ftemp(:) !WSE     
    enddo !j     
    call XF_CLOSE_GROUP(LID,error)
    call XF_CLOSE_GROUP(RID,error)
    call XF_CLOSE_GROUP(TID,error)
    call XF_CLOSE_GROUP(BID,error)
    call XF_CLOSE_GROUP(WID,error)
    call XF_CLOSE_GROUP(GID,error)
    call XF_CLOSE_FILE(PID,error)
    deallocate(ftemp,ftemp2)
     
    return
    endsubroutine read_multiwsevelh5
#endif

endmodule bnd_lib