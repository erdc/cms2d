!====================================================================
module nest_lib
! Nesting Library
!
! written by Alex Sanchez, USACE-CHL
!====================================================================
#include "CMS_cpp.h"
    implicit none

contains

!***********************************************************************************
    subroutine read_parent_cmcards(ctlfilepar,tjuldaypar,&
        grdfilepar,projpar,typespathpar,xoriginpar,yoriginpar,orientpar,telpargrd,&
        wsefilepar,wsepathpar,velfilepar,velpathpar)
! Reads the parent CMS-Flow Control (Card) File.
! A description of the input/output variables is provided below.
! written by Alex Sanchez, USACE-CHL
!***********************************************************************************
    use geo_def, only: projection,projfl
    use time_lib, only: calendar2julian,julianday2calendarmonthday
    use diag_lib
    use prec_def
    implicit none
    !Input/Output
    character(len=*),intent(in)    :: ctlfilepar              !Parent Control File
    real(ikind),     intent(out)   :: tjuldaypar              !Parent simulation reference (start) time Julian days
    character(len=*),intent(inout) :: grdfilepar,typespathpar !Parent Grid File and cell types path
    type(projection),intent(out)   :: projpar                 !Parent grid projection
    real(ikind),     intent(out)   :: xoriginpar,yoriginpar,orientpar !Parent grid origin and orientation
    logical,         intent(out)   :: telpargrd                       !Is parent grid telescoping
    character(len=*),intent(inout) :: wsefilepar,wsepathpar   !Parent water level solution file and path
    character(len=*),intent(inout),optional :: velfilepar,velpathpar  !Parent current velocity solution file and path    
    !Internal Variabl
    integer :: k,ierr
    integer :: iyrpar,imopar,idaypar,ihrpar,iminpar,isecpar,jdaypar
    character(len=200) :: flowpathpar,flownamepar,apath,aname
    character :: cardname*37,aext*10,simlabelpar*100,cdum*34
    logical :: foundfile
    
    xoriginpar=0.0
    yoriginpar=0.0
    orientpar=0.0
    
    write(*,*) 'Reading parent CMS grid control file:'
    write(*,*) trim(ctlfilepar)
    
    call fileparts(ctlfilepar,flowpathpar,flownamepar,aext)
    
    simlabelpar = 'Simulation'    
    
    inquire(file=ctlfilepar,exist=foundfile)
    if(.not.foundfile)then
      call diag_print_error('Could not find parent control file: ',ctlfilepar)
    endif
    open(88,file=ctlfilepar,status='unknown')
    telpargrd = .false.
    do
      read(88,*,iostat=ierr) cardname
      if(ierr/=0) exit
      selectcase(cardname)
        case('SIMULATION_LABEL')
          backspace(88)
          read(88,*) cardname, simlabelpar
        
        case('GRID_CELL_TYPES')
          backspace(88)
          read(88,*) cardname, typespathpar
        
        case('GRID_FILE')
          backspace(88)
          read(88,*) cardname, grdfilepar
            call fileparts(grdfilepar,apath,aname,aext)
          if(len_trim(apath)==0)then
            grdfilepar = trim(flowpathpar) // grdfilepar
          endif
        
        case('TELESCOPING','TELESCOPING_FILE')
          backspace(88)
          read(88,*) cardname, grdfilepar
          telpargrd = .true.
          call fileparts(grdfilepar,apath,aname,aext)
          if(len_trim(apath)==0)then
            grdfilepar = trim(flowpathpar) // grdfilepar
          endif
        
        case('GRID_ANGLE','GRID_ORIENTATION')
          call card_scalar(88,'deg','deg',orientpar,ierr)
          orientpar = 360.0 - orientpar  !Alex
          
        case('GRID_ORIGIN_X')
          call card_scalar(88,'m','m',xoriginpar,ierr)
         
        case('GRID_ORIGIN_Y')
          call card_scalar(88,'m','m',yoriginpar,ierr)
        
        case('STARTING_DATE_TIME','START_DATE_TIME')
          call card_datetime(88,iyrpar,imopar,idaypar,ihrpar,iminpar,isecpar) !YYYY-MM-DD HH:MM:SS UTC
          call calendar2julian(iyrpar,imopar,idaypar,ihrpar,iminpar,isecpar,tjuldaypar)
          
        case('STARTING_JDATE')
          backspace(88)
          read(88,*) cardname, cdum
          read(cdum(1:2),*) iyrpar
          if(iyrpar<100)then
            if(iyrpar>50)then       !if(iyrpar>20)then    MEB Modified - 02/22/2021   This changed all dates in 2021 to 1921.  SMS is now using the card STARTING_DATE_TIME, so this shouldn't be an issue much longer.
              iyrpar = iyrpar + 1900
            else
              iyrpar = iyrpar + 2000
            endif
          endif
          read(cdum(3:5),*) jdaypar !julian day (1-366)
          call julianday2calendarmonthday(iyrpar,jdaypar,imopar,idaypar)
    
        case('STARTING_JDATE_HOUR')
          backspace(88)
          read(88,*) cardname, ihrpar
          iminpar = 0; isecpar = 0
          call calendar2julian(iyrpar,imopar,idaypar,ihrpar,iminpar,isecpar,tjuldaypar)
        
        case('HORIZONTAL_PROJECTION_BEGIN','HORIZ_PROJ_BEGIN')
          call proj_horiz_block(88,projpar)
        
        case('VERTICAL_PROJECTION_BEGIN','VERT_PROJ_BEGIN')
          call proj_vert_block(88,projpar)
          
        case('GLOBAL_WATER_LEVEL_OUTPUT')
          backspace(88)
          read(88,*) cardname, wsefilepar !,wsepathpar  !Path written by SMS is incorrect
          call fileparts(wsefilepar,apath,aname,aext)
          if(len_trim(apath)==0)then
            wsefilepar = trim(flowpathpar) // wsefilepar
          endif
        
        case('WSE_OUT_FILE')
          backspace(88)
          read(88,*) cardname, wsefilepar  
          call fileparts(wsefilepar,apath,aname,aext)
          if(len_trim(apath)==0)then
            wsefilepar = trim(flowpathpar) // wsefilepar
          endif        
          
        case('GLOBAL_VELOCITY_OUTPUT')
          if(present(velfilepar))then
            backspace(88)          
            read(88,*) cardname, velfilepar !,velpathpar  !Path written by SMS is incorrect
            call fileparts(velfilepar,apath,aname,aext)
            if(len_trim(apath)==0)then
              velfilepar = trim(flowpathpar) // velfilepar
            endif
          endif
          
        case('VEL_OUT_FILE')
          if(present(velfilepar))then  
            backspace(88)
            read(88,*) cardname, velfilepar     
            call fileparts(velfilepar,apath,aname,aext)
            if(len_trim(apath)==0)then
              velfilepar = trim(flowpathpar) // velfilepar
            endif
          endif
          
      endselect
    enddo
741 close(88)
    
    if(projfl%iHorizDatum==2 .or. projfl%iHorizCoordSystem==22)then
      if(projpar%iHorizDatum/=2 .or. projpar%iHorizCoordSystem/=22)then  !LOCAL
        call diag_print_error('If the child CMS-Flow simulation has a local projection',&
           '  then the Parent CMS-Flow Simulatation must also have a local projection',&
           '  either in the parent or child CMS-Flow card files.')
      endif
    else
      if(projpar%iHorizDatum==2 .or. projpar%iHorizCoordSystem==22)then  !LOCAL
        call diag_print_error('If the child CMS-Flow simulation does not have a local projection',&
           '  then the parent CMS-Flow simulation must also have a non-local projection.')
      endif
    endif    
    
    !Check if paths are unspecified and set defaults
    if(len_trim(wsefilepar)==0)then
      wsefilepar = trim(flowpathpar) // trim(flownamepar) // '_wse.h5'
      inquire(file=wsefilepar,exist=foundfile) 
      if(.not.foundfile)then
        wsefilepar = trim(flowpathpar) // trim(flownamepar) // '_sol.h5'
        inquire(file=wsefilepar,exist=foundfile) 
        if(.not.foundfile)then
          call diag_print_error('Could not find parent wse solution file: ',wsefilepar)
        endif
      endif
    endif
    if(len_trim(wsepathpar)==0)then
      wsepathpar = trim(simlabelpar) // '/Water_Elevation'
    endif
    
    if(len_trim(velfilepar)==0)then
      velfilepar = trim(flowpathpar) // trim(flownamepar) // '_vel.h5'
      inquire(file=velfilepar,exist=foundfile) 
      if(.not.foundfile)then
        velfilepar = trim(flowpathpar) // trim(flownamepar) // '_sol.h5'
        inquire(file=velfilepar,exist=foundfile) 
        if(.not.foundfile)then
          call diag_print_error('Could not find parent velocity solution file: ',velfilepar)
        endif
      endif
    endif
    if(len_trim(velpathpar)==0)then  
      velpathpar = trim(simlabelpar) // '/Current_Velocity'    
    endif
    
    return
    end subroutine read_parent_cmcards

!**************************************************************************    
    subroutine read_parent_grid_tel(grdfilepar,ncellspar,ncellsfullpar,&
                  xpar,ypar,dxpar,dypar,c2cpar,idfpar,ncfpar,activepar)
! Reads a parent telescoping grid. A description of the input/output variables
! is provided below.
! written by Alex Sanchez, USACE-CHL
!**************************************************************************
    use const_def, only: deg2rad
    use diag_lib
    use prec_def
    implicit none    
    !Input/Output
    integer,      intent(out)         :: ncellspar         !Grid size of active cells on parent grid
    integer,      intent(out)         :: ncellsfullpar     !Grid size of full grid
    integer,      intent(out),pointer :: activepar(:)      !1-active,0-inactive
    integer,      intent(out),pointer :: c2cpar(:,:)       !Cell-to-cell connectivity
    integer,      intent(out),pointer :: idfpar(:,:)       !Direction of cell face (1-N,2-E,3-South,4-West)
    integer,      intent(out),pointer :: ncfpar(:)         !Number of cell faces at each cell
    real(ikind),  intent(out),pointer :: xpar(:),ypar(:)   !Local coordinates of full parent grid     
    real(ikind),  intent(out),pointer :: dxpar(:),dypar(:) !Cell-dimensiones of full parent grid     
    character(len=*),intent(in)          :: grdfilepar        !Grid file of parent grid
    !Internal variables
    integer :: i,id,k,kk,loctemp(8)
    real(ikind) :: val,ztemp
    logical :: ok
    
    inquire(file=grdfilepar,exist=ok)
    if(.not.ok)then
      call diag_print_error('Could not find parent grid: ',grdfilepar,&
         '   Make sure path is assigned to parent grid.')
    endif
    
    open(44,file=grdfilepar,status='unknown')
    read(44,*) !Skip header
    read(44,*) val,val,val,ncellsfullpar !Do not read grid angle or coordinates since they may be incorrect 
    allocate(xpar(ncellsfullpar),ypar(ncellsfullpar),dxpar(ncellsfullpar),dypar(ncellsfullpar))
    allocate(c2cpar(ncellsfullpar,6),idfpar(ncellsfullpar,6),ncfpar(ncellsfullpar),activepar(ncellsfullpar))
    do i=1,ncellsfullpar
      read(44,*) id,xpar(i),ypar(i),dxpar(i),dypar(i),(loctemp(k),k=1,8),ztemp
      ncfpar(i) = 0
      do k=1,8
        if(loctemp(k)/=0)then
          ncfpar(i) = ncfpar(i) + 1
          kk=ncfpar(i)
          c2cpar(i,kk) = loctemp(k)
          selectcase(k)
            case(1,2); idfpar(i,kk)=1
            case(3,4); idfpar(i,kk)=2
            case(5,6); idfpar(i,kk)=3
            case(7,8); idfpar(i,kk)=4
          endselect
        endif        
      enddo
      if(abs(ztemp+999.0)>1.0e-4)then
        activepar(i)=1
      else
        activepar(i)=0  
      endif
    enddo    
    close(44)
    
    ncellspar = sum(activepar)
    
    return
    end subroutine read_parent_grid_tel
    
#ifdef XMDF_IO
!****************************************************************************************    
    subroutine read_parent_grid_xmdf(grdfilepar,typespathpar,ncellspar,ncellsfullpar,&
                  xpar,ypar,dxpar,dypar,c2cpar,idfpar,ncfpar,activepar)
! Reads the Cartesian from the XMDF file.
! The XMDF file is only used for regular and nonuniform Cartesian grids. 
! The telescoping Cartesian grids are stored in the *.tel file.
! written by Alex Sanchez, USACE-CHL
!****************************************************************************************
    use const_def, only: deg2rad
    use xmdf
    use prec_def
    implicit none
    !Input/Output
    integer,      intent(out)         :: ncellspar                       !Grid size of active cells on parent grid
    integer,      intent(out)         :: ncellsfullpar                   !Grid size of full grid
    real(ikind),  intent(out),pointer :: xpar(:),ypar(:)                 !Global coordinates of full parent grid 
    real(ikind),  intent(out),pointer :: dxpar(:),dypar(:)               !Global coordinates of full parent grid 
    integer,      intent(out),pointer :: activepar(:)                    !1-active,0-inactive
    integer,      intent(out),pointer :: c2cpar(:,:)                     !Cell-to-cell connectivity
    integer,      intent(out),pointer :: idfpar(:,:)                     !Direction of cell face (1-N,2-E,3-South,4-West)
    integer,      intent(out),pointer :: ncfpar(:)                       !Number of cell faces at each cell
    character(len=*),intent(in)          :: grdfilepar,typespathpar         !Grid file, and path to activity dataset
    !Internal variables
    integer :: i,j,id,iloc,k,FID,RID,TID,ierr
    integer :: maxcolpar,maxrowpar   !Maximum number of columns and rows on parent grid
    real(4),allocatable :: dxxpar(:),dyypar(:)     !Must be single for XMDF libraries
    character(len=200) :: proppathpar,rootpathpar
    
    !Read XMDF grid file
    call XF_OPEN_FILE(grdfilepar,READONLY,FID,ierr)    
    iloc=index(typespathpar,'/',BACK=.TRUE.)
    proppathpar = typespathpar(1:iloc)                  
    call XF_OPEN_GROUP(FID,trim(proppathpar),TID,ierr)
    call XF_GET_PROPERTY_NUMBER(TID,'CellTypes',ncellsfullpar,ierr)
    allocate(activepar(ncellsfullpar))
    call XF_READ_PROPERTY_INT(TID,'CellTypes',ncellsfullpar,activepar(1),ierr)
    iloc=index(grdfilepar,'/')
    rootpathpar = grdfilepar(1:iloc-1)
    call XF_OPEN_GROUP(FID,trim(rootpathpar),RID,ierr)
    call XF_GET_PROPERTY_NUMBER(RID,'CoordsI',maxcolpar,ierr)
    call XF_GET_PROPERTY_NUMBER(RID,'CoordsJ',maxrowpar,ierr)
    allocate(dxxpar(maxcolpar),dyypar(maxrowpar))
    call XF_READ_PROPERTY_FLOAT(RID,'CoordsI',maxcolpar,dxxpar(1),ierr)
    call XF_READ_PROPERTY_FLOAT(RID,'CoordsJ',maxrowpar,dyypar(1),ierr)    
    call XF_CLOSE_GROUP(RID,ierr)
    call XF_CLOSE_FILE(FID,ierr)

    allocate(xpar(ncellsfullpar),ypar(ncellsfullpar),dxpar(ncellsfullpar),dypar(ncellsfullpar))
    allocate(c2cpar(ncellsfullpar,4),idfpar(ncellsfullpar,4),ncfpar(ncellsfullpar),activepar(ncellsfullpar))
    
    !Parent flow grid global coordinates
    do j=1,maxrowpar        
      do i=1,maxcolpar    
        id = i + (j-1)*maxcolpar  
        c2cpar(id,1) = i     + j*maxcolpar     !North        
        c2cpar(id,2) = i + 1 + (j-1)*maxcolpar !East
        c2cpar(id,3) = i +     (j-2)*maxcolpar !South  
        c2cpar(id,4) = i - 1 + (j-1)*maxcolpar !West
        do k=1,4
          if(c2cpar(id,k)>ncellsfullpar .or. c2cpar(id,k)<1) c2cpar(id,k) = 0  
        enddo
        if(i==1)then
          dxpar(id) = dxxpar(i)
        else
          dxpar(id) = (dxxpar(i) - dxxpar(i-1)) 
        endif
        if(j==1)then
          dypar(id) = dyypar(j)
        else
          dypar(id) = (dyypar(j) - dyypar(j-1))
        endif    
        xpar(id) = dxxpar(i) - 0.5*dxpar(id)  !Local x-coordinate 
        ypar(id) = dyypar(j) - 0.5*dypar(id)  !Local y-coordinate
        activepar(id)=min(activepar(id),1)
        idfpar(id,1:4)=(/1:4/)
        ncfpar(id)=4
      enddo
    enddo
    ncellspar = sum(activepar)
    
    deallocate(dxxpar,dyypar)
    
    return
    end subroutine read_parent_grid_xmdf

!********************************************************************************
    subroutine readparscalsteph5(afile,apath,nptspar,inc,thrs,var,ierr)
! Reads a vector from a parent CMS XMDF solution file.
! Because the whole dataset has to be read in, there is no point in mapping the
! data from the full grid (SMS grid with all cells) to the active grid 
! (CMS grid with only active ocean cells).
! written by Alex Sanchez, USACE-CHL
!********************************************************************************
    use xmdf
    use prec_def
    implicit none
    !Input/Output
    integer, intent(in) :: nptspar
    integer, intent(in) :: inc !Time step index
    integer, intent(out) :: ierr
    integer :: fid,gid,ntimes
    real(ikind), intent(inout) :: var(nptspar),thrs
    character(len=*),intent(in) :: afile,apath
    !Internal
    real*8,allocatable :: timesd(:)
    real*4 :: vtemp(nptspar) !Must be single    

    call XF_OPEN_FILE(trim(afile),READONLY,fid,ierr)        
4   if(ierr<0) return
          
    call XF_OPEN_GROUP(fid,trim(apath),gid,ierr)
    if(ierr<0) return    
    
    call XF_GET_PROPERTY_NUMBER(gid,'Times',ntimes,ierr)
    if(ierr<0) return
    if(inc>ntimes)then
      ierr = -999
      return
    endif
    
    allocate(timesd(ntimes))
    call XF_GET_DATASET_TIMES(gid,ntimes,timesd,ierr)
    if(ierr<0) return
    thrs = real(timesd(inc)) !Hours
    deallocate(timesd)
    
    call XF_READ_SCALAR_VALUES_TIMESTEP(gid,inc,nptspar,vtemp,ierr)
    if(ierr<0) return
    var = vtemp
    
    call XF_CLOSE_GROUP(gid,ierr)
    call XF_CLOSE_FILE(fid,ierr)
      
    return
    end subroutine readparscalsteph5

!********************************************************************************
    subroutine readparvecsteph5(afile,apath,inc,nptspar,thrs,vecx,vecy,ierr)
! Reads a vector from a parent CMS XMDF solution file
! Because the whole dataset has to be read in, there is no point in mapping the
! data from the full grid (SMS grid with all cells) to the active grid 
! (CMS grid with only active ocean cells).
!
! written by Alex Sanchez, USACE-CHL
!********************************************************************************
    use xmdf
    use prec_def
    implicit none
    !Input/Output
    integer, intent(in) :: nptspar
    integer, intent(in) :: inc !Time step index
    integer, intent(out) :: ierr
    real(ikind), intent(inout) :: vecx(nptspar),vecy(nptspar),thrs
    character(len=*),intent(in) :: afile,apath
    !Internal
    integer :: fid,gid,ntimes,i
    real*8,allocatable :: timesd(:)
    real*4 :: vtemp(nptspar*2) !Must be single

    call XF_OPEN_FILE(trim(afile),READONLY,fid,ierr)        
    if(ierr<0) return
          
    call XF_OPEN_GROUP(fid,trim(apath),gid,ierr)
    if(ierr<0) return    
    
    call XF_GET_PROPERTY_NUMBER(gid,'Times',ntimes,ierr)
    if(ierr<0) return
    if(inc>ntimes)then
      ierr = -999
      return
    endif
    
    allocate(timesd(ntimes))
    call XF_GET_DATASET_TIMES(gid,ntimes,timesd,ierr)
    if(ierr<0) return
    thrs = real(timesd(inc)) !Hours
    deallocate(timesd)
    
    call XF_READ_VECTOR_VALUES_TIMESTEP(gid,inc,nptspar,2,vtemp,ierr)   
    if(ierr<0) return
    do i=1,nptspar
      vecx(i) = vtemp(2*i-1)
      vecy(i) = vtemp(2*i)
    enddo
    
    call XF_CLOSE_GROUP(gid,ierr)
    call XF_CLOSE_FILE(fid,ierr)
      
    return
    end subroutine readparvecsteph5
#endif
    
!******************************************************************************************
    subroutine nestParentCMS_init(ctlfilepar,tjuldaypar,grdfilepar,ncellsfullpar,projpar,&
        nbnd,xbnd,ybnd,mntp,intp,cntp,wsefilepar,wsepathpar,velfilepar,velpathpar,angvel)
! Initializes a nested boundary condition within a parent CMS simulation.
! A description of the input/output variables is provided below.
! written by Alex Sanchez, USACE-CHL
!******************************************************************************************
#include "CMS_cpp.h"
    use geo_def, only: projection,projfl,azimuth_fl
    use geo_lib, only: proj_horiz_conv
    use interp_lib, only: interp_coef_tel2pts
    use const_def, only: deg2rad
    use diag_lib
    use prec_def
    implicit none
    !Input/Output
    character(len=*),intent(in)    :: ctlfilepar    !Parent CMS-Flow Control (Card) File
    real(ikind),     intent(out)   :: tjuldaypar    !Parent simulation reference time in Julian Days
    character(len=*),intent(inout) :: grdfilepar    !Parent CMS-Flow Grid File
    integer,         intent(out)   :: ncellsfullpar !Number of output points (cells) on parent grid
    type(projection),intent(out)   :: projpar       !Parent simulation projection
    integer,         intent(in)    :: nbnd          !Number of boundary cells/nodes
    real(ikind),     intent(in)    :: xbnd(nbnd)    !x-coordinate (global) of boundary point (cell/nodes)
    real(ikind),     intent(in)    :: ybnd(nbnd)    !y-coordinate (global) of boundary point (cell/nodes) 
    integer,         intent(out)   :: mntp          !Maximum number of interpolation cells
    integer,         intent(inout),pointer :: intp(:,:) !Interpolation id's (local grid id, interp mesh id's)
    real(ikind),     intent(inout),pointer :: cntp(:,:) !Interpolation coefficients    
    character(len=*),intent(inout) :: wsefilepar,wsepathpar !Parent CMS-Flow wse solution file
    character(len=*),intent(inout),optional :: velfilepar,velpathpar !Parent CMS-Flow velocity solution file
    real(ikind),     intent(out),  optional :: angvel       !Rotation angle to convert from parent to the child local coordinate systems [rad]
    !Internal variables
    integer              :: ncellspar          !Active grid size of parent grid   
    integer              :: nmaxfacespar       !Maximum number of faces
    real(ikind)          :: xoriginpar         !Global x-origin of parent grid (if Cartesian)
    real(ikind)          :: yoriginpar         !Global y-origin of parent grid (if Cartesian)   
    real(ikind)          :: orientpar          !Orientation (grid angle) of parent grid
    integer,     pointer :: activepar(:)       !Global x-coordinate of full parent grid 
    integer,     pointer :: c2cpar(:,:)        !Cell-to-cells connectivity
    integer,     pointer :: idfpar(:,:)        !Direction of cell face (1-N,2-E,3-South,4-West)
    integer,     pointer :: ncfpar(:)          !Number of cell faces on each cell
    real(ikind), pointer :: xpar(:),ypar(:)    !Global coordinates of full parent grid 
    real(ikind), pointer :: dxpar(:),dypar(:)  !Cell dimensiones of full parent grid 
    real(ikind)          :: xtrapdist
    character(len=200)   :: typespathpar
    logical              :: telpargrd
    
    !Read Card file
    if(present(velfilepar) .and. present(velpathpar))then
      call read_parent_cmcards(ctlfilepar,tjuldaypar,grdfilepar,projpar,&
        typespathpar,xoriginpar,yoriginpar,orientpar,telpargrd,&
        wsefilepar,wsepathpar,velfilepar,velpathpar)
    else
      call read_parent_cmcards(ctlfilepar,tjuldaypar,grdfilepar,projpar,&
        typespathpar,xoriginpar,yoriginpar,orientpar,telpargrd,&
        wsefilepar,wsepathpar)
    endif
    
    if(present(angvel)) angvel = (azimuth_fl - orientpar)*deg2rad
    
    !Read Grid
    if(telpargrd)then !Telescoping grid
      nmaxfacespar = 6  
      call read_parent_grid_tel(grdfilepar,ncellspar,ncellsfullpar,&
             xpar,ypar,dxpar,dypar,c2cpar,idfpar,ncfpar,activepar)
    else  !Cartesian grid
      nmaxfacespar = 4
#ifdef XMDF_IO
      call read_parent_grid_xmdf(grdfilepar,typespathpar,ncellspar,ncellsfullpar,&
             xpar,ypar,dxpar,dypar,c2cpar,idfpar,ncfpar,activepar)
#else
      call diag_print_error('Cannot read grid from *.h5 file without XMDF libraries')
#endif
    endif
    
    !Reproject if necessary
    call proj_horiz_conv(projpar,projfl,ncellsfullpar,xpar,ypar)
    
    !Calculate interpolation coefficients    
    xtrapdist = 1.0e6 !Set large value to use nearest neighbor value
    mntp = 4  
    allocate(intp(0:mntp,nbnd),cntp(mntp,nbnd))
    intp(:,:) = 0; cntp(:,:) = 0.0
    call interp_coef_tel2pts(ncellspar,ncellsfullpar,nmaxfacespar,&
            xOriginpar,yOriginpar,orientpar, &
            xpar,ypar,dxpar,dypar,c2cpar,idfpar,ncfpar, &            
            nbnd,xbnd,ybnd,xtrapdist,intp,cntp)
    
    !if(debug_mode)then
    !  open(57,file='IntpCoefNestCMS.txt')  
    !  do i=1,nbnd
    !    write(57,*) i,(intp(k,i),k=0,mntp),(cntp(k,i),k=1,mntp)
    !  enddo
    !  close(57)
    !endif
    
    return
    end subroutine nestParentCMS_init
    
!***********************************************************************************
    subroutine nestParentADCIRC_init(grdfilepar,nnodesfullpar,projpar,&
       nbnd,xbnd,ybnd,mntp,intp,cntp,wsefilepar,wsepathpar,velfilepar,velpathpar,angvel)
! Initializes the nesting of a CMS grid within a larger ADCIRC grid
! Calculates the interpolation coefficients
!***********************************************************************************
    use geo_def, only: projection,projfl,azimuth_fl
    use geo_lib, only: read_grid14,proj_horiz_conv
    use const_def, only: deg2rad
    use interp_lib, only: interp_coef_tri2pts
    use diag_lib
    use prec_def
    implicit none
    !Input/Output
    character(len=*),intent(in)  :: grdfilepar     !Grid file  
    type(projection),intent(out) :: projpar        !Parent grid projection
    integer,         intent(in)  :: nbnd           !Number of boundary points
    real(ikind),     intent(in)  :: xbnd(nbnd)     !x-coordinate (global) of boundary points (cells/nodes)
    real(ikind),     intent(in)  :: ybnd(nbnd)     !y-coordinate (global) of boundary points (cells/nodes)
    integer,         intent(out) :: nnodesfullpar  !Number of output points (cells) on parent grid
    integer,         intent(out)  :: mntp          !Maximum number of interpolation cells
    integer,         intent(inout),pointer :: intp(:,:)   !Interpolation id's (local grid id, interp mesh id's)
    real(ikind),     intent(inout),pointer :: cntp(:,:)   !Interpolation coefficients
    character(len=*),intent(inout):: wsefilepar    !Parent water level solution file
    character(len=*),intent(inout):: wsepathpar    !Parent water level solution path
    character(len=*),intent(inout),optional :: velfilepar !Parent current velocity solution file
    character(len=*),intent(inout),optional :: velpathpar !Parent current velocity solution path
    real(ikind),     intent(out),  optional :: angvel         !Roation angle needed for rotate velocities onto CMS grid
    !Internal variables
    integer              :: nelemsfullpar            !Active grid size of parent grid    
    integer,     pointer :: elem2node(:,:)           !Element to node connectivity
    real(ikind), pointer :: xpar(:),ypar(:),zpar(:)  !Global coordinates of full parent grid
    real(ikind) :: xtrapdist
    character :: aext*10,aline*100
    logical :: foundfile
       
    !Read grid
    call read_grid14(grdfilepar,nelemsfullpar,nnodesfullpar,xpar,ypar,zpar,elem2node)    
    deallocate(zpar) !Not needed
    
    !Projection Convertion
    if(projpar%iHorizDatum==2 .or. projpar%iHorizCoordSystem==2)then !Local
      call diag_print_warning('No horizontal projection defined for ADCIRC grid',&
        '   Assuming Geographic, NAD83, Degrees')
      projpar%iHorizDatum = 1         !Horizontal Datum = NAD83
      projpar%iHorizCoordSystem = 0   !Horizontal Coordinate System = GEOGRAPHIC
      projpar%iHorizUnits = 4         !Horizontal Units = DEGREES
      projpar%iVertDatum = 9          !Vertical Datum = LOCAL
      projpar%iVertUnits = 2          !Vertical units = METERS
      projpar%VertOffset = 0.0       !Vertical offset from datum
    endif
    
    !Reproject parent grid to child grid projection
    call proj_horiz_conv(projpar,projfl,nnodesfullpar,xpar,ypar)    
    
    !Calculate interpolation coefficients
    xtrapdist = 1.0e6 !Set large value to use nearest neighbor value
    mntp = 3
    allocate(intp(0:mntp,nbnd),cntp(mntp,nbnd))
    intp = 0; cntp = 0.0  
    call interp_coef_tri2pts(nelemsfullpar,nnodesfullpar,xpar,ypar,elem2node, &            
             nbnd,xbnd,ybnd,xtrapdist,intp,cntp)
    
    inquire(file=wsefilepar,exist=foundfile)  
    if(.not.foundfile)then
      call diag_print_error('Could not find parent water level solution file: ',wsefilepar)
    endif  
    call fileext(wsefilepar,aext)
    if(aext(1:2)=='63')then
      wsepathpar = '' !Not used  
      open(63,file=wsefilepar)
      read(63,*) aline !Skip first line
      read(63,*) aline !Skip second line     
    elseif(aext(1:2)=='h5')then
      if(len_trim(wsepathpar)==0)then
        wsepathpar = 'Datasets/Water Surface Elevation (63)'  
        call diag_print_warning('Path not specified for depth-averaged current velocity in XMDF file',&
          '  Setting path to default: ',wsepathpar)      
      endif
    endif
    
    if(.not.present(velfilepar)) return
    
    inquire(file=velfilepar,exist=foundfile)  
    if(.not.foundfile)then
      call diag_print_error('Could not find parent current velocity solution file: ',velfilepar)
    endif   
    call fileext(velfilepar,aext)
    if(aext(1:2)=='64')then
      velpathpar = '' !Not used  
      open(64,file=velfilepar,status='old')
      read(64,*) !Skip first line
      read(64,*) !Skip second line      
    elseif(aext(1:2)=='h5')then
      if(len_trim(velpathpar)==0)then
        velpathpar = 'Datasets/Depth-averaged Velocity (64)'  
        call diag_print_warning('Path not specified for depth-averaged current velocity in XMDF file',&
          '  Setting path to default: ',velpathpar)
      endif
    endif
    
    if(present(angvel)) angvel = azimuth_fl*deg2rad
    
    return
    end subroutine nestParentADCIRC_init
       
!**********************************************************************************
    subroutine tdb_init(tdbname,tdbpath,projtdb,nbnd,xbnd,ybnd,niter,mwin,&
        ntcin,namein,ntc,name,speed,f,vu,etamp,etpha,utamp,utpha,vtamp,vtpha,angvel)
! Initializes the nested tidal database boundary condition
! Supports the ADCIRC EC2001 and ENPAC2003 databases as well as the 
! LeProvost, FES952 and FES2004 databases.
! The option is provided to smooth (running mean) the consituents alon the boundary
! since sharp gradients can occur for coarse databases due to linear interpolation.
! A description of the input and output variables is provided below.
! written by Alex Sanchez, USACE-CHL
!**********************************************************************************
    use geo_def, only: azimuth_fl
    use geo_def, only: projection,projfl
    use geo_lib, only: proj_horiz_conv
    use comvarbl, only: tjulhr0,tjulhryr,iyr
    use const_def, only: deg2rad,twopi
    use tide_lib, only: tide_adcirc,tide_fes,tidal_data
    use bnd_def, only: ntf
    use diag_lib
    use prec_def
    implicit none
    !Input/Output
    integer,         intent(in)    :: nbnd           !Number of boundary points
    real(ikind),     intent(in)    :: xbnd(nbnd)     !x-coordinte (global) of points (usually boundary cell centroids or nodes)
    real(ikind),     intent(in)    :: ybnd(nbnd)     !y-coordinte (global) of points (usually boundary cell centroids or nodes)
    integer,         intent(in)    :: niter          !Spatial Smoothing iterations
    integer,         intent(inout) :: mwin           !Spatial Smoothing window width
    character(len=*),intent(in)    :: tdbname        !Tidal database file
    character(len=*),intent(in)    :: tdbpath        !Tidal database path
    type(projection),intent(inout) :: projtdb        !Tidal database projection
    integer,         intent(in)    :: ntcin          !Number of input constituents
    character(len=6),     intent(in)    :: namein(ntcin)  !Names of input constituents
    integer,         intent(inout) :: ntc            !Number of output constituents
    character(len=6),     intent(out), pointer :: name(:) !Names of output constituents
    real(ikind),     intent(out), pointer :: speed(:),f(:),vu(:) !Constituent speed, nodal factor and equilibrium argument
    real(ikind),     intent(out), pointer :: etamp(:,:),etpha(:,:) !Water level amplitude and phase
    real(ikind),     intent(inout), pointer, optional :: utamp(:,:),utpha(:,:) !U-velocity amplitude (m) and phase (deg)
    real(ikind),     intent(inout), pointer, optional :: vtamp(:,:),vtpha(:,:) !V-Velocity amplitude (m) and phase (deg)
    real(ikind),     intent(out),   optional :: angvel
    !Internal variables
    integer :: i,ii,j,jj,k,ind(ntf),mcyctemp(ntf)
    real(ikind) :: distmin,dist,xout(nbnd),yout(nbnd)
    real(ikind) :: speedtemp(ntf),ftemp(ntf),vutemp(ntf)
    character(len=6) :: nametemp(ntf)
    character(len=100) :: msg2,msg3,msg4,msg5
    
    !Get corrections for all constituents correspondint that year
    call tidal_data(iyr,nametemp,speedtemp,mcyctemp,ftemp,vutemp)
    
    !If Projection is not specified (local), then assume Geographic, NAD83, degrees
    if(projtdb%iHorizDatum==2)then !2=Local
      call diag_print_warning('Tidal Database Horizontal Projection not specified',&
        ' Assuming Projection: Geographic, NAD83, degrees')
      projtdb%iHorizDatum = 1         !Horizontal Datum = NAD83
      projtdb%iHorizCoordSystem = 0   !Horizontal Coordinate System = GEOGRAPHIC
      projtdb%iHorizUnits = 4         !Horizontal Units = DEGREES
      projtdb%iVertDatum = 9          !Vertical Datum = LOCAL
      projtdb%iVertUnits = 2          !Vertical units = METERS
      projtdb%VertOffset = 0.0        !Vertical offset from datum
    endif
    
    !Reproject input coordinates to tidal database horizontal projection
    xout(:) = xbnd(:); yout(:) = ybnd(:)
    call proj_horiz_conv(projfl,projtdb,nbnd,xout,yout)
    
    !Read parent grid and control file
    if(tdbname(1:6)=='EC2001' .or. tdbname(1:9)=='ENPAC2003')then !ADCIRC Tidal Database   
      if(present(vtpha))then
        call tide_adcirc(tdbname,tdbpath,nbnd,xout,yout,&
           ntcin,namein,ntc,name,speed,etamp,etpha,utamp,utpha,vtamp,vtpha)      
      else
        call tide_adcirc(tdbname,tdbpath,nbnd,xout,yout,&
           ntcin,namein,ntc,name,speed,etamp,etpha)     
      endif
    elseif(tdbname(1:9)=='LEPROVOST' .or. tdbname(1:3)=='FES')then !LeProvost Tidal Database
      call tide_fes(tdbname,tdbpath,nbnd,xout,yout,ntcin,namein,ntc,name,etamp,etpha)  
      allocate(speed(ntc))
      speed = -1.0
      !Check for undefined points
      do ii=1,nbnd
        if(etamp(ii,1)<-900.0)then !if undefined use neighest defined point
          distmin = 1.0e20
          do jj=1,nbnd
            if(jj==ii .or. etamp(jj,1)<-900.0) cycle
            dist = sqrt((xout(ii)-xout(jj))**2+(yout(ii)-yout(jj))**2)
            if(dist<distmin) j = jj
          enddo
          etamp(ii,:)=etamp(j,:)
          etpha(ii,:)=etpha(j,:)
        endif        
      enddo
      !allocate(utamp(nbnd,ntc),utpha(nbnd,ntc))
      !allocate(vtamp(nbnd,ntc),vtpha(nbnd,ntc))
      !utamp = 0.0; utpha = 0.0
      !vtamp = 0.0; vtpha = 0.0
      
    elseif(tdbname(1:4)=='TPXO')then !TPXO Tidal Database (all versions)
        !Not implimented yet
    else
      call diag_print_error('Invalid Parent Grid Control File Extention')
    endif
    
    !Calculate index for names
    ind = 0
    do k=1,ntc
      do j=1,ntf
        if(name(k)==nametemp(j))then
          ind(k) = j
          exit
        elseif(name(k)=='STEADY')then
          ind(k) = 0
          exit
        elseif(j==ntf)then
          ind(k) = -1  
          write(msg2,*) '  Name: ',name(k)
          write(msg3,*) '  Speed: ',speed(k)
          write(msg4,*) '  Setting nodal factor to 1.0'
          write(msg5,*) '  Setting equilibrium argument to 0.0'
          call diag_print_warning('Tidal constituent not found',msg2,msg3,msg4,msg5)
        endif
      enddo !j
    enddo !k
    
    !Check corrections and speeds
    allocate(f(ntc),vu(ntc))
    do k=1,ntc
      i = ind(k)  
      if(i<=0)then      
        f(k) = 1.0
        vu(k) = 0.0
      else
        f(k) = ftemp(i)
        vu(k) = vutemp(i)       
        if(speed(k)<0.0)then
          speed(k) = speedtemp(i)
        elseif(abs(speed(k)-speedtemp(i))>1.0e-4)then
          write(msg2,*) '  Name: ',name(k)
          write(msg3,*) '  Speed: ',speed(k)
          write(msg4,*) '  Setting nodal factor to 1.0'
          write(msg5,*) '  Setting equilibrium argument to 0.0'
          call diag_print_warning('Invalid tidal constituent ',&
            '  speed found in tidal database',msg2,msg3,msg4,msg5)
          read(*,*)
          stop
        endif
      endif
    enddo
    
    !Smoothing along string. 
    if(niter>0 .and. mwin>=3)then
      do k=1,ntc
        call smoothampdir(niter,mwin,nbnd,etamp(:,k),etpha(:,k))
      enddo        
    endif
    
    !Corrections to amplitudes and phases
    do k=1,ntc
      vu(k)      = vu(k)*deg2rad  
      speed(k)   = speed(k)*deg2rad      
      etpha(:,k) = etpha(:,k)*deg2rad
    enddo
    
    if(.not.present(angvel)) return
    
    angvel = azimuth_fl*deg2rad !Radians
    
    !Smoothing along cellstring
    if(niter>0 .and. mwin>=3)then
      do k=1,ntc
        call smoothampdir(niter,mwin,nbnd,utamp(:,k),utpha(:,k))
        call smoothampdir(niter,mwin,nbnd,vtamp(:,k),vtpha(:,k))   
      enddo     
    endif
    
    do k=1,ntc
      utpha(:,k) = utpha(:,k)*deg2rad
      vtpha(:,k) = vtpha(:,k)*deg2rad 
    enddo
    
    return
    end subroutine tdb_init
        
!***************************************************************************
    subroutine interp_nest(nptspar,nbnd,mntp,intp,cntp,valdry,varpar,varnest)
! Interpolates a variable from a parent grid to a boundary cellstring of
! the child grid.
! written by Alex Sanchez, USACE-CHL
!***************************************************************************
    use prec_def
    implicit none
    !Input/Output
    integer,    intent(in) :: nptspar
    integer,    intent(in) :: nbnd
    integer,    intent(in) :: mntp
    integer,    intent(in) :: intp(0:mntp,nbnd)
    real(ikind),intent(in) :: cntp(mntp,nbnd)
    real(ikind),intent(in) :: valdry
    real(ikind),intent(in) :: varpar(nptspar)
    real(ikind),intent(inout) :: varnest(nbnd)
    !Internal Variables
    integer :: j,k
    real(ikind) :: csum
            
    do j=1,nbnd
      varnest(j) = 0.0
      csum=0.0
      do k=1,intp(0,j)
        !Check for wetting and drying, CMS uses NaN = -999.0, and ADCIRC uses NaN = -9999.0  
        if(varpar(intp(k,j))>-90.0)then !Valid NaN values are all <=90.0
          csum = csum + cntp(k,j)
          varnest(j) = varnest(j) + cntp(k,j)*varpar(intp(k,j))
        endif             
      enddo
      if(csum<1.0e-4)then !All points NaN, make cell dry
        varnest(j) = valdry  !Dry
      else
        varnest(j) = varnest(j)/csum !Wet    
      endif   
      !if(abs(varnest(j))>1.0e15)then  
      !  write(*,*) j,varnest(j)  
      !endif
    enddo
    
    return
    end subroutine interp_nest
    
!************************************************************************************
    subroutine nest_wse_eval(ntipar,timewsehrspar,nptspar,wsepar,&
      nbnd,ibndcells,mntp,intp,cntp,nti,timewsehrs,wsedata,wsebnd,wseout,wsefile)
! Reads the parent solution file and performs spatial and temporal interpolations. 
! The wse is first spatially interpolated and then interpolated in time.
! A description of the input/output variables is provided below
! written by Alex Sanchez, USACE-CHL
!************************************************************************************
    use comvarbl, only: timehrs,tjulday0,flowpath
    use flow_def, only: iwet
    use plagr_lib
    use out_lib, only: open_tsd,append_tsd
    use prec_def
    implicit none
    !Input/Output
    integer,         intent(in)    :: ntipar                   !Maximum order of interpolation for parent grid
    real(ikind),     intent(in)    :: timewsehrspar(ntipar+1)  !Current velocity solution times used for temporal interpolation
    integer,         intent(in)    :: nptspar                  !# of solution points of the parent grid
    real(ikind),     intent(in)    :: wsepar(nptspar,ntipar+1) !WSE data for parent simulation [m]
    integer,         intent(in)    :: nbnd                !# of boundary cells/nodes
    integer,         intent(in)    :: ibndcells(nbnd)     !Boundary cell ID's
    integer,         intent(in)    :: mntp                !Maximum # of interpolation weights per boundary element (cell/node)
    integer,         intent(in)    :: intp(0:mntp,nbnd)   !Interpolation indeces from the parent to child grids
    real(ikind),     intent(in)    :: cntp(mntp,nbnd)     !Interpolation coefficients from the parent to child grids    
    integer,         intent(in)    :: nti                 !Order of interpolation of child grid
    !integer,         intent(inout) :: inc                 !Time increment of the current velocity solution
    real(ikind),     intent(inout) :: timewsehrs(nti+1)   !Current velocity solution times used for temporal interpolation
    real(ikind),     intent(inout) :: wsedata(nbnd,nti+1) !Water levels at times 1 and 2 used for interpolation
    real(ikind),     intent(out)   :: wsebnd(nbnd)        !Interpolated water level at timehrs
    logical,         intent(in)    :: wseout              !Write wse output
    character(len=*),intent(in)    :: wsefile             !Output wse file
    !Internal variables
    integer :: j,jj,ierr,ntimes,kth,nti1,ipar
    integer :: indgap(nbnd)
    integer :: i,k,np
    integer,parameter :: nb = 3 !Must be greater than largest nti
    real(ikind) :: lb(nb),wsum,timewsesec,valdry
    
    !if(wseout)then
    !  write(wsefile,'(A,I1,A)')  'Nest_bnd_wse',idnum,'.tsd'
    !  wsefile = trim(flowpath) // wsefile
    !endif
    nti1 = nti + 1
        
    valdry = -999.0_ikind
    
    !Read dataset if necessary
    if(timewsehrs(2)<1.0e-5)then !First time
      !Spatial interpolation from parent grid to nested grid
      do k=1,nti1
        timewsehrs(k) = timewsehrspar(k)  
        call interp_nest(nptspar,nbnd,mntp,intp,cntp,valdry,wsepar(:,k),wsedata(:,k)) !Dry cells are assigned -999.0
      enddo
      if(wseout)then
        call open_tsd(wsefile,'Boundary WSE','Unassigned',nbnd,0,tjulday0)
        do k=1,nti1
          timewsesec = timewsehrs(k)*3600.0_ikind !Convert from hrs to seconds
          call append_tsd(wsefile,nbnd,timewsesec,wsedata(:,k))
        enddo
      endif
    endif
    
    !Read new time step(s) if necesary
    kth = max(nti,2) !Better than using the last value. Note nti is at least 1
    ierr = 0  
    do while(timehrs>timewsehrs(kth) .and. ierr>=0)
      do k=1,nti
        timewsehrs(k) = timewsehrs(k+1)
        do i=1,nbnd
          wsedata(i,k) = wsedata(i,k+1)
        enddo
      enddo
      ipar=1
      do i=1,ntipar+1
        if(timewsehrspar(i)>timewsehrs(nti1))then
          ipar=i
          exit
        endif
      enddo
      timewsehrs(nti1) = timewsehrspar(ipar)
      !Spatial interpolation from parent grid to nested grid
      call interp_nest(nptspar,nbnd,mntp,intp,cntp,valdry,&
        wsepar(:,ipar),wsedata(:,nti1)) !Dry cells are assigned -999.0
      if(wseout)then
        timewsesec = timewsehrs(nti1)*3600.0_ikind !Convert from hrs to seconds
        call append_tsd(wsefile,nbnd,timewsesec,wsedata(:,nti1))
      endif
    enddo      
    
    !Temporal interpolation
    if(timehrs<=timewsehrs(1)+1.0e-5)then
      do i=1,nbnd  
        wsebnd(i) = wsedata(i,1)
      enddo
    elseif(timehrs>=timewsehrs(nti1))then
      do i=1,nbnd  
        wsebnd(i) = wsedata(i,nti1)
      enddo
    else
      np = nti; ntimes = nti1; k = 1
      call plagr_fit(ntimes,timewsehrs,timehrs,nb,lb,nti,np,k)  
      do i=1,nbnd
        wsebnd(i) = 0.0  
        wsum = 0.0
        do j=1,np+1
          jj=k+j-1  
          if(wsedata(i,jj)>-99.0)then !Not dry
            wsebnd(i) = wsebnd(i) + lb(j)*wsedata(i,jj)
            wsum = wsum + lb(j)
          endif
        enddo
        if(abs(wsum)<1.0e-15)then
          wsebnd(i) = -999.0 !Dry
        else  
          wsebnd(i) = wsebnd(i)/wsum !Wet
        endif  
      enddo
    endif
    
    !Check wetting and drying of boundary cells    
    do j=1,nbnd
      if(iwet(ibndcells(j))==1 .and. wsebnd(j)==-999.0)then !Inconsistent boundary
        indgap(j)=1 !Dry
      else
        indgap(j)=0 !Wet
      endif
    enddo        
    call interp_gap(nbnd,indgap,wsebnd)
    
    !if(timehrs<=timehr1+1.0e-5)then
    !  wsebnd = wsebnd1
    !elseif(timehrs>=timehr2-1.0e-5)then
    !  wsebnd = wsebnd2
    !else !if(timehrs>timehr1 .and. timehrs<timehr2)then
    !  fac=(timehrs-timehr1)/(timehr2-timehr1)
    !  fac=max(min(fac,1.0),0.0) !Avoids extrapolation
    !  wsebnd = (1.0-fac)*wsebnd1 + fac*wsebnd2   
    !endif

    return
    end subroutine nest_wse_eval
    
!*************************************************************************************
    subroutine nest_vel_eval(ntipar,timevelhrspar,nptspar,upar,vpar,&
      nbnd,ibndcells,mntp,intp,cntp,nti,timevelhrs,udata,vdata,ubnd,vbnd,velout,velfile)
! Reads the parent solution file and performs spatial and temporal interpolations.
! The wse is first spatially interpolated and then interpolated in time.
! A description of the input/output variables is provided below
! written by Alex Sanchez, USACE-CHL
!*************************************************************************************
    use comvarbl, only: timehrs,tjulday0,flowpath
    use flow_def, only: iwet
    use plagr_lib
    use out_lib, only: open_tsd,append_tsd
    use diag_lib
    use prec_def
    implicit none
    !Input/Output
    integer,         intent(in)    :: ntipar                  !Order of temporal interpolation [-]
    real(ikind),     intent(in)    :: timevelhrspar(ntipar+1) !Current velocitiy solution times used for temporal interpolation
    integer,         intent(in)    :: nptspar                 !# of solution points of the parent grid
    real(ikind),     intent(in)    :: upar(nptspar,ntipar+1)  !U-Current velocity data for parent simulation [m/s]
    real(ikind),     intent(inout) :: vpar(nptspar,ntipar+1)  !V-Current velocity data for parent simulation [m/s]
    integer,         intent(in)    :: nbnd              !# of boundary cells
    integer,         intent(in)    :: ibndcells(nbnd)   !Boundary cell ID's
    integer,         intent(in)    :: mntp              !Maximum # of interpolation weights per boundary element (cell/node)
    integer,         intent(in)    :: intp(0:mntp,nbnd) !Interpolation indeces from the parent to child grids
    real(ikind),     intent(in)    :: cntp(mntp,nbnd)   !Interpolation coefficients from the parent to child grids
    integer,         intent(in)    :: nti               !Order of temporal interpolation [-]
    !integer,         intent(inout) :: inc               !Time increment of the parent current velocity solution
    real(ikind),     intent(inout) :: timevelhrs(nti+1) !Current velocitiy solution times used for temporal interpolation
    real(ikind),     intent(inout) :: udata(nbnd,nti+1) !U-Current velocity data [m/s]
    real(ikind),     intent(inout) :: vdata(nbnd,nti+1) !V-Current velocity data [m/s]
    real(ikind),     intent(out)   :: ubnd(nbnd)        !Interpolated current velocites at timehrs [m/s]
    real(ikind),     intent(out)   :: vbnd(nbnd)        !Interpolated current velocites at timehrs [m/s]    
    logical,         intent(in)    :: velout            !Output velocity
    character(len=*),intent(in)    :: velfile           !Output velocity file
    !Internal variables
    integer :: i,j,jj,k,ierr,kth,np,ntimes,nbnd2,nti1,ipar
    !integer :: indgap(nbnd)
    real(ikind) :: veldata(nbnd*2)
    integer,parameter :: nb = 3 !Must be greater than largest nti
    real(ikind) ::lb(nb),timevelsec,valdry
    character(len=100) :: msg2,msg3,msg4,msg5
    
    if(velout)then
      !write(velfile,'(A,I1,A)')  'Nest_bnd_vel',idnum,'.tsd'
      !velfile = trim(flowpath) // velfile
      nbnd2 = nbnd*2
    endif
    
    nti1 = nti + 1
    valdry = 0.0_ikind
    
    !Read dataset if necessary
    if(timevelhrs(2)<1.0e-5)then      
      !Spatial interpolation from parent grid to nested grid
      do k=1,nti1
        timevelhrs(k) = timevelhrspar(k)
        call interp_nest(nptspar,nbnd,mntp,intp,cntp,valdry,upar(:,k),udata(:,k)) !Dry cells are set to 0.0
        call interp_nest(nptspar,nbnd,mntp,intp,cntp,valdry,vpar(:,k),vdata(:,k)) !Dry cells are set to 0.0
      enddo    
      if(velout)then
        call open_tsd(velfile,'Boundary Vel','Unassigned',nbnd2,0,tjulday0)  
        do k=1,nti1
          do j=1,nbnd
            jj=j*2  
            veldata(jj-1) = udata(j,k)
            veldata(jj)   = vdata(j,k)   
          enddo
          timevelsec = timevelhrs(k)*3600.0 !Convert from hrs to seconds
          call append_tsd(velfile,nbnd2,timevelsec,veldata)
        enddo  
      endif
    endif
    
    !Read new time(s) step
    kth = max(nti,2)
    ierr = 0  
    do while(timehrs>timevelhrs(kth) .and. ierr>=0)
      do k=1,nti
        timevelhrs(k) = timevelhrs(k+1)
        do i=1,nbnd
          udata(i,k) = udata(i,k+1)
          vdata(i,k) = vdata(i,k+1)
        enddo
      enddo
      ipar=1
      do i=1,ntipar+1
        if(timevelhrspar(i)>timevelhrs(nti1))then
          ipar=i
          exit  
        endif
      enddo
      timevelhrs(nti1) = timevelhrspar(ipar)
      !Spatial interpolation from parent grid to nested grid
      call interp_nest(nptspar,nbnd,mntp,intp,cntp,valdry,upar(:,ipar),udata(:,nti1)) !Dry cells are set to 0.0
      call interp_nest(nptspar,nbnd,mntp,intp,cntp,valdry,vpar(:,ipar),vdata(:,nti1)) !Dry cells are set to 0.0
      if(velout)then
        do j=1,nbnd
          jj=j*2  
          veldata(jj-1) = udata(j,k)
          veldata(jj)   = vdata(j,k)   
        enddo
        timevelsec = timevelhrs(nti1)*3600.0_ikind !Convert from hrs to seconds
        call append_tsd(velfile,nbnd2,timevelsec,veldata)
      endif
    enddo
    
    !Temporal interpolation
    if(timehrs<=timevelhrs(1))then
      do i=1,nbnd    
        ubnd(i) = udata(i,1)
        vbnd(i) = vdata(i,1)
        if(abs(ubnd(i))>3.0_ikind .or. abs(vbnd(i))>3.0_ikind)then
          write(msg2,*) '   Cell: ',i
          write(msg3,*) '   Velocities: ',ubnd(i),vbnd(i)
          write(msg4,*) '   Interpolation Index: ',intp(:,i)
          write(msg5,*) '   Interpolation Coeff: ',cntp(:,i)
          call diag_print_warning('Large forcing velocities',msg2,msg3,msg4,msg5)
        endif
      enddo
    elseif(timehrs>=timevelhrs(nti1))then
      do i=1,nbnd  
        ubnd(i) = udata(i,nti1)
        vbnd(i) = vdata(i,nti1)
        if(abs(ubnd(i))>3.0_ikind .or. abs(vbnd(i))>3.0_ikind)then
          write(msg2,*) '   Cell: ',i
          write(msg3,*) '   Velocities: ',ubnd(i),vbnd(i)
          write(msg4,*) '   Interpolation Index: ',intp(:,i)
          write(msg5,*) '   Interpolation Coeff: ',cntp(:,i)
          call diag_print_warning('Large forcing velocities found',msg2,msg3,msg4,msg5) 
        endif
      enddo
    else
      np = nti; ntimes=nti1; k = 1
      call plagr_fit(ntimes,timevelhrs,timehrs,nb,lb,nti,np,k)  
      do i=1,nbnd
        !Note: dry cells have zero velocity and are included in the temporal interpolation
        ubnd(i)=valdry !==0.0
        vbnd(i)=valdry !==0.0
        do j=1,np+1
          jj=k+j-1  
          ubnd(i) = ubnd(i) + lb(j)*udata(i,jj)
          vbnd(i) = vbnd(i) + lb(j)*vdata(i,jj)
        enddo
        if(abs(ubnd(i))>3.0_ikind .or. abs(vbnd(i))>3.0_ikind)then
          write(msg2,*) '   Cell: ',i
          write(msg3,*) '   Velocities: ',ubnd(i),vbnd(i)
          write(msg4,*) '   Interpolation Index: ',intp(:,i)
          write(msg5,*) '   Interpolation Coeff: ',cntp(:,i)
          call diag_print_warning('Large forcing velocities found',msg2,msg3,msg4,msg5)   
        endif
      enddo
    endif
    
    !if(timehrs<=timehr1+1.0e-5)then
    !  ubnd = ubnd1
    !  vbnd = vbnd1
    !elseif(timehrs>=timehr2-1.0e-5)then
    !  ubnd = ubnd2
    !  vbnd = vbnd2
    !else !if(timehrs>timehr1 .and. timehrs<timehr2)then
    !  fac=(timehrs-timehr1)/(timehr2-timehr1)
    !  fac=max(min(fac,1.0),0.0) !Avoids extrapolation
    !  ubnd = (1.0-fac)*ubnd1 + fac*ubnd2      
    !  vbnd = (1.0-fac)*vbnd1 + fac*vbnd2       
    !endif    
    
    !!Check wetting and drying of boundary cells    
    !do j=1,nbnd
    !  if(iwet(ibndcells(j))==1 .and. ubnd(j)==0.0 .and. vbnd(j)==0.0)then !Inconsistent boundary
    !    indgap(j)=1
    !  else
    !    indgap(j)=0  
    !  endif
    !enddo        
    !call interp_gap(nbnd,indgap,ubnd)
    !call interp_gap(nbnd,indgap,vbnd)
    
    !Check for large velocities
    do i=1,nbnd
      if(abs(ubnd(i))>3.0_ikind .or. abs(vbnd(i))>3.0_ikind)then
        write(msg2,*) '   Cell: ',i
        write(msg3,*) '   Velocities: ',ubnd(i),vbnd(i)
        write(msg4,*) '   Interpolation Index: ',intp(:,i)
        write(msg5,*) '   Interpolation Coeff: ',cntp(:,i)
        call diag_print_warning('Large forcing velocities found',msg2,msg3,msg4,msg5)       
      endif
    enddo
    
    return
    end subroutine nest_vel_eval

!***************************************************    
    subroutine interp_gap(n,igap,a)
!***************************************************
    use diag_lib
    use prec_def
    implicit none      
    !Input/Output
    integer :: n
    integer :: igap(n) !1-gap, 0-not gap    
    real(ikind) :: a(n)
    !Internal variables
    integer :: i,j,j1,j2
    real(ikind) :: fac
    character(len=100) :: msg2
    
    do j=1,n !Index loop
      if(igap(j)==1)then !is gap          
        j1=0 !Left index
        do i=max(j-1,1),1,-1
          if(igap(i)==0)then !Find first non-gap value
            j1=i
            exit
          endif   
        enddo
        j2=0 !Right index
        do i=min(j+1,n),n
          if(igap(i)==0)then !Find first non-gap value
            j2=i 
            exit
          endif  
        enddo        
        !Fill-in gap
        if(j1>0 .and. j2>0)then !Interpolate linearly
          fac=real(j-j1,kind=ikind)/real(j2-j1,kind=ikind) 
          fac=max(min(fac,1.0),0.0) !Avoids extrapolation
          a(j)=(1.0-fac)*a(j1)+fac*a(j2)
        elseif(j1>0)then !Left value
          a(j)=a(j1)
        elseif(j2>0)then !Right value
          a(j)=a(j2)
        else
          write(msg2,*) '  Boundary String Cell:',j
          call diag_print_warning('Could not fill-in gap',msg2)
        endif
      endif !is gap
    enddo !index loop
    
    return
    end subroutine interp_gap

end module nest_lib
