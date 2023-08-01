!=================================================================================
! CMS Geopatial routines
!
! Contains the following:
!   General
!     geo_default   - Sets the default geospatial parameters
!     geo_cards     - Reads the geospatial cards from the control file 
!     geo_init_cart - Initializes geospatial variables for Cartesian grids 
!     geo_init_poly - Initializes geospatial variables for Polyhedral meshes
!     geo_print     - Prints the geospatial settings to the diagnostic file 
!                     and the screen
!   Grids I/O
!     read_grid_xmdf - Reads the SMS XMDF grid file *_grid.h5 for Nonuniform Cartesian grids
!     read_tel       - Reads an SMS ASCII telescoping grid file *.tel
!     read_2dm       - Reads an SMS ASCII 2D mesh file *.2dm
!   Mapping (indexing)
!     map_cell_full2active - Maps cell ID's from the full grid to 
!                            the active grid
!     map_scal_full2active - Maps a scalar array from the full grid 
!                            to the active grid
!     map_vec_full2active - Maps a vector array from the full grid 
!                            to the active grid
!   Horizontal and Vertical Projections
!     proj_horiz_block - Reads a horizontal projection block
!     proj_vert_block  - Reads a vertical projection block
!     get_coord_info - Reads the horizontal projection information from the
!                      CMS-Flow XMDF grid file *_grid.h5
!
! written by Alex Sanchez, USACE-CHL
!=================================================================================

!*************************************************************
    subroutine geo_default()
! Sets the default values to geospatial variables
! written by Alex Sanchez, USACE-CHL
!*************************************************************
    use size_def
    use geo_def
    use geo_lib, only: proj_default
    use prec_def
    implicit none
    
    !Grid size
    ncellsfull = 0   !Full cells
    ncells = 0       !# of active cells
    ncellsD = 0      !# of cells including dummy (ghost) cells
    ncellsimple = 0  !# of simple (regular) cells
    ncelljoint = 0   !# of joint (telescoping) cells
    ncellpoly = 0    !# of  polygonal cells
    nnodes = 0       !# of nodes
    nmaxfaces = 4    !Maximum # of faces in x and y directions    
    ndmaxfaces = 2   !Maximum # of faces in each direction
    nmaxcells = 4    !Maximum # of cells per node
    
    !Internal CONSTANT variables    
    kkface(1)=3     !index for opposit face 
    kkface(2)=4
    kkface(3)=1
    kkface(4)=2
    signface(1)=1.0    !sign at face for flux divergence term in pp equation
    signface(2)=1.0
    signface(3)=-1.0
    signface(4)=-1.0
    
    !Grid parameters    
    igridtype = 0    !0-Non-telescoping, 1-Telescoping, 2-Polygonal, 3-Curvilinear
    avg_lat = 0.0
    azimuth_fl = 0.0
    xorigin = 0.0
    yorigin = 0.0
    nflowgrd = 0
    latpath = ''
    lonpath = ''
    
    !Projection defaults
    call proj_default(projfl)
    HProj = ''
    VProj = ''
    
    !Horizontal Coordinate datum    
    aHorizDatum(0) = 'NAD27'
    aHorizDatum(1) = 'NAD83'
    aHorizDatum(2) = 'LOCAL'
    
    !Horizontal Coordinate system    
    aHorizCoordSystem(0)  = 'GEOGRAPHIC'
    aHorizCoordSystem(1)  = 'UTM'
    aHorizCoordSystem(2)  = 'STATE_PLANE'
    aHorizCoordSystem(3)  = 'ALBERS_CONICAL_EQUAL-AREA'
    aHorizCoordSystem(4)  = 'LAMBERT_CONFORMAL_CONIC'
    aHorizCoordSystem(5)  = 'MERCATOR'
    aHorizCoordSystem(6)  = 'POLAR_STEREOGRAPHIC'
    aHorizCoordSystem(7)  = 'POLYCONIC'
    aHorizCoordSystem(8)  = 'EQUIDISTANT_CONIC'
    aHorizCoordSystem(9)  = 'TRANSVERSE_MERCATOR'
    aHorizCoordSystem(10) = 'STEREOGRAPHIC'
    aHorizCoordSystem(11) = 'LAMBERT_AZIMUTHAL_EQUAL-AREA'
    aHorizCoordSystem(12) = 'AZIMUTHAL_EQUIDISTANT'
    aHorizCoordSystem(13) = 'GNOMONIC'
    aHorizCoordSystem(14) = 'ORTHOGRAPHIC'
    aHorizCoordSystem(15) = 'GENERAL_VERTICAL_NEAR-SIDE_PERSPECTIVE'
    aHorizCoordSystem(16) = 'SINUSOIDAL'
    aHorizCoordSystem(17) = 'EQUIRECTANGULAR' !(PLATE CARREE)
    aHorizCoordSystem(18) = 'MILLER_CYLINDRICAL'
    aHorizCoordSystem(19) = 'VAN_DER_GRINTEN_I'
    aHorizCoordSystem(20) = 'OBLIQUE_MERCATOR' !(HOTINE)
    aHorizCoordSystem(21) = 'SPACE_OBLIQUE_MERCATOR'
    aHorizCoordSystem(22) = 'LOCAL'

    !Horizontal Coordinate Units    
    aHorizUnits(0) = 'RADIANS'
    aHorizUnits(1) = 'FEET'
    aHorizUnits(2) = 'METERS'
    aHorizUnits(3) = 'SECONDS'
    aHorizUnits(4) = 'DEGREES'
    aHorizUnits(5) = 'US_FEET'
    aHorizUnits(6) = 'INTL_FEET'

    !Vertical Coordinate Datum    
    aVertDatum(0) = 'SWL'
    aVertDatum(1) = 'MSL'
    aVertDatum(2) = 'MTL'
    aVertDatum(3) = 'MLLW'
    aVertDatum(4) = 'MLW'
    aVertDatum(5) = 'MHHW'
    aVertDatum(6) = 'MHW'
    aVertDatum(7) = 'NAVD88'
    aVertDatum(8) = 'NGVD29'
    aVertDatum(9) = 'LOCAL'
    
    !Vertical Coordinate Units    
    projfl%iVertUnits = 2
    aVertUnits(1) = 'FEET'
    aVertUnits(2) = 'METERS'
    aVertUnits(3) = 'US_FEET'
    aVertUnits(4) = 'INTL_FEET'
    projfl%VertOffset = 0.0  !Vertical offset from datum
    
    return
    end subroutine geo_default

!*************************************************************
    subroutine geo_cards(cardname,foundcard,doPrint)
! Reads the input vards from the CMS-Flow Control File and
! overrides default values
! written by Alex Sanchez, USACE-CHL
!*************************************************************    
    use geo_def
    use comvarbl, only: mpfile, flowpath, casename, SMS_ver
    use prec_def
    implicit none
    
    !Input/Output
    character(len=*),intent(inout) :: cardname
    logical,intent(out) :: foundcard
    logical,intent(in)  :: doPrint
    
    !Internal Variables
    integer :: ierr
    character(len=200) :: apath,aname
    character(len=10) :: aext
    character(len=37) :: cdum
    
    foundcard = .true.
    select case(cardname)        
      case('GRID_MODIFICATION_NUMBER')
        backspace(77)
        read(77,*) cardname, nflowgrd
      
      case('GRID_FILE')
        if (telfile(1:1) .eq. ' ') then
          backspace(77)
          read(77,*) cardname, grdfile    
          call fileparts(grdfile,apath,aname,aext)  
          if(len_trim(apath)==0) grdfile = trim(flowpath) // grdfile
          call lowercase(aext)
          if(aext(1:3)=='2dm') igridtype = 2  !0-Non-telescoping, 1-Telescoping, 2-Polygonal, 3-Curvilinear
          if(aext(1:3)=='tel') igridtype = 1
        endif
        
      case('TELESCOPING','TELESCOPING_FILE')
        backspace(77)
        read(77,*) cardname, telfile
        call fileparts(telfile,apath,aname,aext)  
        if(len_trim(apath)==0) telfile = trim(flowpath) // telfile
          
      case('BATHYMETRY_DATASET')
        backspace(77)   
        read(77,*) cardname, cdum, grdpath
        
      case('BATHYMETRY_UPDATE_DATASET','BATHYMETRY_DATASETS')
        backspace(77)   
        read(77,*) cardname, bathydata%file, bathydata%path
        bathydata%ison = .true.
          
      case('GRID_CELL_TYPES')
        backspace(77)
        read(77,*) cardname, typespath
          
      case('GRID_ANGLE')
        call card_scalar(77,'deg','deg',azimuth_fl,ierr)
        if (SMS_ver .le. 13.0) then
          azimuth_fl = 360.0 - azimuth_fl  !Alex
        endif
              
      case('GRID_ORIGIN_X')
        call card_scalar(77,'m','m',xorigin,ierr)
          
      case('GRID_ORIGIN_Y')
        call card_scalar(77,'m','m',yorigin,ierr)
          
      case('CELL_LATITUDES')
        call card_dataset(77,mpfile,flowpath,latfile,latpath,0)  !0 for Lats and Lons - 05/21/2018
        
      case('CELL_LONGITUDES')
        call card_dataset(77,mpfile,flowpath,lonfile,lonpath,0)  !0 for Lats and Lons - 05/21/2018 
          
      case('AVERAGE_LATITUDE','LATITUDE_AVERAGE')
        call card_scalar(77,'deg','deg',avg_lat,ierr)     
        
      case('AVERAGE_LONGITUDE')
        backspace(77)
!!        read(77,*) cardname, avg_lon
        read(77,*) cardname  
      
      case('HORIZONTAL_PROJECTION_BEGIN','HORIZ_PROJ_BEGIN')
        if(doPrint) call proj_horiz_block(77,projfl)            !Only process if in main section not the Solution Scheme part
        
      case('VERTICAL_PROJECTION_BEGIN','VERT_PROJ_BEGIN')
        if(doPrint) call proj_vert_block(77,projfl)             !Only process if in main section not the Solution Scheme part
          
      case default
          foundcard = .false.  
                          
    end select
    
    return
    end subroutine geo_cards

!**************************************************
    subroutine geo_print()
! Prints the geospatial parameters
! written by Alex Sanchez, USACE-CHL
!**************************************************   
    use size_def
    use diag_lib
    use diag_def, only: dgunit,dgfile
    use out_def, only: write_ascii_input
    use geo_def    
    use tool_def, only: vstrlz
    implicit none
    integer :: i,iunit(2)
    real    :: val
    character(len=200) :: apath,aname,astring
    character(len=10) :: aext

133 format(' ',A,T40,F0.3,A)
144 format(' ',A,T40,A,A)
222 format(' ',A,T40,I0)
787 format(' ',A,T40,A)
888 format(' ',A)
    
    iunit = (/6, dgunit/)
    
    open(dgunit,file=dgfile,access='append') 
    
    val = 360.0 - azimuth_fl
    do i=1,2    
      write(iunit(i),*)
      write(iunit(i),888)    'Grid '         
      
      if (telfile .eq. " ") then
        call fileparts(grdfile,apath,aname,aext)
      else
        call fileparts(telfile,apath,aname,aext)
      endif
      astring=trim(aname) // '.' // aext
      write(iunit(i),787)    '  Grid File:',trim(astring)
      
      if(igridtype<=1)then
        if(grdpath(1:1)/=' ') write(iunit(i),787)  '  Depth Dataset:',trim(grdpath)  
        write(iunit(i),144)  '  Orientation:',trim(vstrlz(val,'(f0.3)')),' deg'    !Writes out the same value that was written in the .cmcards file - MEB 01/13/2017
        write(iunit(i),144)  '  x-Origin:',trim(vstrlz(xorigin,'(f0.3)')),' m'
        write(iunit(i),144)  '  y-Origin:',trim(vstrlz(yorigin,'(f0.3)')),' m'
      endif
      write(iunit(i),144)    '  Average Latitude:',trim(vstrlz(avg_lat,'(f0.3)')),' deg'
      write(iunit(i),222)    '  Total cells:',ncellsD
      write(iunit(i),222)    '  Active cells:',ncells
      if(igridtype==1)then
        write(iunit(i),787)  '  Telescoping:','ON'
        write(iunit(i),222)  '  Regular cells:',ncellsimple
        write(iunit(i),222)  '  Joint cells:',ncelljoint
      else
        write(iunit(i),787)  '  Telescoping:','OFF'
      endif 
!    write(*,*) 'CMS-Flow Latitudes Dataset: ',LATPATH
!    write(dgunit,*) 'CMS-Flow Latitudes Dataset: ',LATPATH
      write(iunit(i),888)    '  Horizontal Projection'
      write(iunit(i),787)    '    Coordinate System:',trim(aHorizCoordSystem(projfl%iHorizCoordSystem))     
      if(projfl%iHorizCoordSystem/=22)then
        write(iunit(i),787)  '    Datum:',trim(aHorizDatum(projfl%iHorizDatum))
        write(iunit(i),222)  '    Zone:',projfl%iHorizZone
      endif      
      write(iunit(i),787)    '    Units:',trim(aHorizUnits(projfl%iHorizUnits))
      
    enddo

    close(dgunit)
    
    if(write_ascii_input)then
      if(igridtype<1)then !Nonuniform Cartesian grid
        call diag_print_message('Writing ASCII grid file')
        call grid_cart_write_ascii  
      !elseif
      endif
    endif
    
    return
    end subroutine geo_print
     
!*************************************************************
    subroutine read_grid_xmdf()
! Reads grid from XMDF file
!
! writen by Alex Sanchez, USACE-ERDC-CHL
! last modified 120811
!*************************************************************          
#include "CMS_cpp.h"
    use size_def
    use geo_def
    use diag_lib
    use comvarbl
#ifdef XMDF_IO  
    use xmdf
#endif
    use fric_def, only: constbotfric,cbotfric
    use bnd_def      
    implicit none
    integer :: ILOC,Icheck,ierr
    integer :: FILE_ID,TYPES_ID,DSET_ID,ROOT_ID
    
    nmaxfaces = 4  !Maximum # of faces in all directions
    ndmaxfaces = 2 !Maximum # of faces in each direction

!---- Read XMDF file ---------------------------------------------
#ifdef XMDF_IO
    call XF_OPEN_FILE(grdfile,READONLY,FILE_ID,ierr)
      
    !Grid parameters
    iloc = index(typespath,'/',back=.true.)
    proppath = typespath(1:iloc)                  
    call XF_OPEN_GROUP(FILE_ID,trim(proppath),TYPES_ID,ierr)
    call XF_GET_PROPERTY_NUMBER(TYPES_ID,'CellTypes',ncellsfull,ierr)
      
    iloc = index(grdpath,'/')
    rootpath = grdpath(1:iloc)
    call XF_OPEN_GROUP(FILE_ID,trim(rootpath),ROOT_ID,ierr)
    call XF_GET_PROPERTY_NUMBER(ROOT_ID,'CoordsI',maxcol,ierr)
    call XF_GET_PROPERTY_NUMBER(ROOT_ID,'CoordsJ',maxrow,ierr)

    !check to see if grid strucgture is valid
    Icheck = maxrow*maxcol - ncellsfull
    if(Icheck/=0) then
      call diag_print_error('Grid structure is not understood by this code')
    endif
    
    allocate(z(ncellsfull),dxx(maxcol),dyy(maxrow))
    allocate(cell_type(ncellsfull)) 

    !GET DX/DY COORDINATES
    call XF_READ_PROPERTY_FLOAT(ROOT_ID,'CoordsI',maxcol,dxx(1),ierr)
    call XF_READ_PROPERTY_FLOAT(ROOT_ID,'CoordsJ',maxrow,dyy(1),ierr)

    !GET DEPTHS
    call XF_OPEN_GROUP(FILE_ID,trim(GRDPATH),DSET_ID,ierr)
    call XF_READ_SCALAR_VALUES_TIMESTEP(DSET_ID,1,ncellsfull,Z,ierr)      
      
    !GET CELLTYPES
    call XF_READ_PROPERTY_INT(TYPES_ID,'CellTypes',ncellsfull,cell_type(1),ierr)
    call XF_CLOSE_FILE(FILE_ID,ierr)

#endif

    call geo_cart

    return
    end subroutine read_grid_xmdf
    
!**************************************************************************
    subroutine read_tel()
! Reads the telescoping grid file
!
! written by Alex Sanchez, USACE-ERDC-CHL  
!**************************************************************************
    use size_def 
    use geo_def
    use flow_def
    use struct_def
    use comvarbl
    use const_def, only: deg2rad,small
    use diag_lib
    use diag_def, only: msg,msg2,msg3
    use prec_def
    implicit none
    
    integer :: i,ii,j,k,id,nck,nck2,ndum
    integer, allocatable :: loctemp(:,:)    
    integer, parameter :: inan = -99
    real(ikind), parameter :: fnan = -999.0          
    real(ikind) :: xg,yg,val,cosAng,sinAng   
    real(ikind), allocatable :: xtemp(:),ytemp(:),dxtemp(:),dytemp(:),ztemp(:)
    
    integer :: nIsolated
    type isoType 
      integer     :: id, col, row
      real(ikind) :: xg,yg
    endtype isoType
    type(isoType), allocatable :: isolated(:)

    nmaxfaces = 6  !Maximum # of faces in all directions
    ndmaxfaces = 4 !Maximum # of faces in each direction

!=== Load temporary grid (full grid including inactive land cells) ===
    open(55,file=telfile)
    read(55,*) !Skip header
    read(55,*) val,val,val,ncellsfull !Do not read grid angle or coordinates since they be incorrect 
    allocate(idmap(0:ncellsfull))       
    allocate(loctemp(ncellsfull,8))
    allocate(xtemp(ncellsfull),ytemp(ncellsfull))
    allocate(dxtemp(ncellsfull),dytemp(ncellsfull))
    allocate(ztemp(0:ncellsfull)) !needed for below    
    ztemp(0) = 0
    idmap(0) = 0
    do ii=1,ncellsfull
      continue
      read(55,*) id,xtemp(ii),ytemp(ii),dxtemp(ii),dytemp(ii),&
          loctemp(ii,1),loctemp(ii,2),loctemp(ii,3),loctemp(ii,4),&
          loctemp(ii,5),loctemp(ii,6),loctemp(ii,7),loctemp(ii,8),&
          ztemp(ii)   
    enddo
    close(55)    

!== Save geometry infomation ======
    cosAng = cos(azimuth_fl*deg2rad)
    sinAng = sin(azimuth_fl*deg2rad)

!=== Determine case size and mapping index ====
    ncells = 0
    ncellsD = 0  
    do ii=1,ncellsfull                
      if(abs(ztemp(ii)-fnan)<small)then
        idmap(ii) = 0
        cycle !Skip land cells
      endif      
      ncells = ncells + 1 !Active cells
      idmap(ii) = ncells    
      
!      if(abs(xtemp(ii)-4555.0+xOrigin)<small .and. abs(ytemp(ii)-2355.0+yOrigin)<small)then
!        continue
!      endif
!      if(abs(xtemp(ii)-2695.0+xOrigin)<small .and. abs(ytemp(ii)-1485.0+yOrigin)<small)then
!        continue
!      endif
!      if(ii==28055)then
!        continue
!      endif
            
      do k=1,7,2
        nck = loctemp(ii,k)
        nck2 = loctemp(ii,k+1)
        if(nck==0 .and. nck2==0)then !At boundary of grid
          ncellsD = ncellsD + 1
          loctemp(ii,k) = inan 
        elseif((nck==0 .or. abs(ztemp(nck)-fnan)<small) .and. &
           (nck2==0 .or. abs(ztemp(nck2)-fnan)<small))then !Both are inactive cells, replace with on
          ncellsD = ncellsD + 1
          loctemp(ii,k) = inan
          loctemp(ii,k+1) = 0
        elseif(nck>0 .and. abs(ztemp(nck)-fnan)<small)then !Only one is inactive
          ncellsD = ncellsD + 1
          loctemp(ii,k) = inan   !Remove id. Will add a new one later
        elseif(nck2>0 .and. abs(ztemp(nck2)-fnan)<small)then !Only one is inactive
          ncellsD = ncellsD + 1
          loctemp(ii,k+1) = inan   !Remove id. Will add a new one later
        endif
      enddo  
!      !Add extra dummy cells for corners      
!      ncf = 0          
!      do j=1,8
!        if(ia(j)>0) ncf = ncf + 1
!      enddo  
!      if(ncf==2)then
!        ncellsD = ncellsD + 1    
!      endif
    enddo
 
    ncellsD = ncellsD + ncells
    
!!    write(*,'(A,I10)') ' Active cells = ',ncells
!!    write(*,'(A,I10)') ' Ghost cells =  ',ncellsD-ncells
!!    write(*,'(A,I10)') ' Total cells =  ',ncellsD 

!=== Map from full grid to active grid ====
    allocate(mapid(ncellsD))
    allocate(ncface(ncellsD))
    allocate(cell2cell(nmaxfaces,ncellsD))
    allocate(idirface(nmaxfaces,ncellsD))
    allocate(zb(ncellsD))
    allocate(zbk(nmaxfaces,ncellsD))
    allocate(x(ncellsD),y(ncellsD))
    allocate(dx(ncellsD),dy(ncellsD))
    ncface = 0
    cell2cell = 0
    idirface = 0
    zb = 0.0
    zbk=0.0
    x = 0.0; y = 0.0
    dx = 0.0; dy = 0.0
    ndum = ncells
    do ii=1,ncellsfull         
      if(idmap(ii)==0) then
        cycle
      endif
      i = idmap(ii) 
      mapid(i) = ii
      x(i) = xtemp(ii)
      y(i) = ytemp(ii)
      dx(i) = dxtemp(ii)
      dy(i) = dytemp(ii)
      zb(i) = -ztemp(ii)  !Note sign change ****
      
      !Calculate ncface, loconnect and idirface      
      k = 0
      do j=1,8        
        nck = loctemp(ii,j)                        
        if(nck==inan)then !Dummy cells
          k = k + 1
          ncface(i) = k
          ndum = ndum + 1
          cell2cell(k,i) = ndum
          ncface(ndum) = 4
          if(j==1)then  !North     !Wu
            idirface(k,i) = 1
            cell2cell(3,ndum) = i
            if(loctemp(ii,2)>0)then
              x(ndum) = x(i) - 0.25 * dx(i)
              y(ndum) = y(i) + 0.75 * dy(i) 
              dx(ndum) = 0.5 * dx(i)
              dy(ndum) = 0.5 * dy(i)
            else
              x(ndum) = x(i)
              y(ndum) = y(i) + dy(i) 
              dx(ndum) = dx(i)
              dy(ndum) = dy(i)
            endif
          endif
          if(j==2)then  !North
            idirface(k,i) = 1
            cell2cell(3,ndum) = i
            if(loctemp(ii,1)>0)then
              x(ndum) = x(i) + 0.25 * dx(i)
              y(ndum) = y(i) + 0.75 * dy(i) 
              dx(ndum) = 0.5 * dx(i)
              dy(ndum) = 0.5 * dy(i)
            else
              x(ndum) = x(i)
              y(ndum) = y(i) + dy(i) 
              dx(ndum) = dx(i)
              dy(ndum) = dy(i)
            endif
          endif
          if(j==3)then !East
            idirface(k,i) = 2
            cell2cell(4,ndum) = i
            if(loctemp(ii,4)>0)then
              x(ndum) = x(i) + 0.75 * dx(i)
              y(ndum) = y(i) + 0.25 * dy(i) 
              dx(ndum) = 0.5 * dx(i)
              dy(ndum) = 0.5 * dy(i)
            else
              x(ndum) = x(i) + dx(i)
              y(ndum) = y(i) 
              dx(ndum) = dx(i)
              dy(ndum) = dy(i)
            endif
          endif
          if(j==4)then !East
            idirface(k,i) = 2
            cell2cell(4,ndum) = i
            if(loctemp(ii,3)>0)then
              x(ndum) = x(i) + 0.75 * dx(i)
              y(ndum) = y(i) - 0.25 * dy(i) 
              dx(ndum) = 0.5 * dx(i)
              dy(ndum) = 0.5 * dy(i)
            else
              x(ndum) = x(i) + dx(i)
              y(ndum) = y(i) 
              dx(ndum) = dx(i)
              dy(ndum) = dy(i)
            endif
          endif
          if(j==5)then !South
            idirface(k,i) = 3
            cell2cell(1,ndum) = i
            if(loctemp(ii,6)>0)then
              x(ndum) = x(i) + 0.25 * dx(i)
              y(ndum) = y(i) - 0.75 * dy(i) 
              dx(ndum) = 0.5 * dx(i)
              dy(ndum) = 0.5 * dy(i)
            else
              x(ndum) = x(i)
              y(ndum) = y(i) - dy(i)
              dx(ndum) = dx(i)
              dy(ndum) = dy(i)
            endif
          endif
          if(j==6)then !South
            idirface(k,i) = 3
            cell2cell(1,ndum) = i
            if(loctemp(ii,5)>0)then
              x(ndum) = x(i) - 0.25 * dx(i)
              y(ndum) = y(i) - 0.75 * dy(i) 
              dx(ndum) = 0.5 * dx(i)
              dy(ndum) = 0.5 * dy(i)
            else
              x(ndum) = x(i)
              y(ndum) = y(i) - dy(i)
              dx(ndum) = dx(i)
              dy(ndum) = dy(i)
            endif
          endif
          if(j==7)then  !West
            idirface(k,i) = 4
            cell2cell(2,ndum) = i
            if(loctemp(ii,8)>0)then
              x(ndum) = x(i) - 0.75 * dx(i)
              y(ndum) = y(i) - 0.25 * dy(i) 
              dx(ndum) = 0.5 * dx(i)
              dy(ndum) = 0.5 * dy(i)
            else
              x(ndum) = x(i) - dx(i)
              y(ndum) = y(i)     
              dx(ndum) = dx(i)
              dy(ndum) = dy(i)
            endif
          endif
          if(j==8)then  !West  
            idirface(k,i) = 4
            cell2cell(2,ndum) = i
            if(loctemp(ii,7)>0)then
              x(ndum) = x(i) - 0.75 * dx(i)
              y(ndum) = y(i) + 0.25 * dy(i) 
              dx(ndum) = 0.5 * dx(i)
              dy(ndum) = 0.5 * dy(i)
            else
              x(ndum) = x(i) - dx(i)
              y(ndum) = y(i)     
              dx(ndum) = dx(i)
              dy(ndum) = dy(i)
            endif
          endif
          zb(ndum) = zb(i)       !Wu
        elseif(nck>0 .and. idmap(nck)>0)then !Internal cells
          k = k + 1
          ncface(i) = k
          cell2cell(k,i) = idmap(nck)
          if(j<=2)then
            idirface(k,i) = 1 !North
          elseif(j<=4)then
            idirface(k,i) = 2 !East
          elseif(j<=6)then
            idirface(k,i) = 3 !South
          else
            idirface(k,i) = 4 !West          
          endif
        endif    
      enddo !j
    enddo !ii
    
!Check
    if(ndum/=ncellsD)then
      call diag_print_error('Problem reading telescoping file')
    endif        

!**** TEMPORARY *******************
!Set constant manning, surface elevation
!    manncont = 0.025
!    do i=1,ncellsD
!      coefman(i) = manncont   !Wu
!      eta(i)=0.0   !Wu
!      p(i) = 0.0 !Initial water surface elevation
!    enddo
!**** TEMPORARY *******************

    allocate( isolated(ncells) )
    nIsolated = 0  
!--- Check for single cell lakes --------------------
    do i=1,ncells
      if (cell2cell(1,i)>ncells .and. & 
          cell2cell(2,i)>ncells .and. &
          cell2cell(3,i)>ncells .and. &
          cell2cell(4,i)>ncells) then
        !Convert from local to global coordinates
        xg=xOrigin+x(i)*cosAng-y(i)*sinAng
        yg=yOrigin+x(i)*sinAng+y(i)*cosAng        

8615 format('  i=' ,I6,' j=',I6)          
8617 format('  id=' ,I6)       
        
    !Modified the code so that the user can get a list of all cell IDs to fix before CMS stops, MEB 01/31/2020
        if(igridtype==0)then
          nIsolated = nIsolated + 1
          isolated(nIsolated)%col = icol(i)
          isolated(nIsolated)%row = irow(i)
          isolated(nIsolated)%xg  = xg
          isolated(nIsolated)%yg  = yg
        else
          nIsolated = nIsolated + 1
          isolated(nIsolated)%id = mapid(i)
          isolated(nIsolated)%xg = xg
          isolated(nIsolated)%yg = yg
        endif 
      endif      
    enddo
    
100 FORMAT ('ERROR: ',i0,1x,'isolated ocean cells identified')
101 FORMAT ('  Column: ',i0,3x,'Row: ',i0,3x,'X: ',f0.5,3x,'Y: ',f0.5)
102 FORMAT ('  ID:',i0,3x,'X:',f0.5,3x,'Y:',f0.5)    
103 FORMAT ('Correct isolated cells before continuing')
    
    if (nIsolated .gt. 0) then
      write(msg,100) nIsolated
      call diag_print_message(' ',msg)
      do i=1,nIsolated
        if(igridtype == 0) then
          write(msg2,101) isolated(i)%col, isolated(i)%row, isolated(i)%xg, isolated(i)%yg
        else
          write(msg2,102) isolated(i)%id, isolated(i)%xg, isolated(i)%yg
        endif
        call diag_print_message(msg2)
      enddo
      write(msg,103)
      call diag_print_message(msg)
      write(*,*) 'Press any key to continue.'
      read(*,*)
      stop
    endif
        
!    if(debug_mode)then
!    !Cell Mapping
!      icasen = index(casename,' ')-1 
!      filename = casename(1:icasen)//'.idmap'
!      open(58,file=filename) 
!      write(58,*) 'ii, idmap(ii)'
!818   format(I8,1x,I8)      
!      do ii=1,ncellsfull
!        write(58,818) ii,idmap(ii)
!      enddo
!      close(58)
!      
!      filename = casename(1:icasen)//'.dbug'
!      open(78,file=filename) 
!      write(78,*) ncells, ncellsD
!      write(78,*) 'i,mapid(i),xg,yg,-zb,ncface(i),(cell2cell(k,i),k=1,ncface(i))'      
!714   format(I7,1x,I7,1x,2(F12.3,1x),F9.3,1x,I6,6(I7,1x))  
!      do i=1,ncells       
!        xg=xOrigin+x(i)*cosAng-y(i)*sinAng
!        yg=yOrigin+x(i)*sinAng+y(i)*cosAng
!        write(78,714) i,mapid(i),xg,yg,-zb(i),ncface(i),(cell2cell(k,i),k=1,ncface(i))    
!        write(78,714) i,mapid(i),dx(i),dy(i),-zb(i),ncface(i),(idirface(k,i),k=1,ncface(i))        
!      enddo
!715   format(I7,1x,I7,1x,2(F12.3,1x),F9.3)        
!      do i=ncells+1,ncellsD
!        xg=xOrigin+x(i)*cosAng-y(i)*sinAng
!        yg=yOrigin+x(i)*sinAng+y(i)*cosAng
!        write(78,715) i,0,xg,yg,-zb(i)
!      enddo
!      close(78)      
!    endif        
    
    deallocate(xtemp,ytemp,dxtemp,dytemp,ztemp,loctemp)     
    deallocate( isolated )
    
    return
    end subroutine read_tel    
    
!*************************************************************
    subroutine read_2dm()
! written by Alex, USACE-ERDC-CHL
!*************************************************************        
    use size_def
    use geo_def
    use diag_lib
    use comvarbl, only: casename,flowpath
    use bnd_def, only: bc_type,bc_str,nbcstr
    use prec_def
    implicit none
    integer :: i,j,id,istart,iend,nseg,nmod,nns,ncs,ierr,nbn,nbc
    integer :: nd,idnstrtemp(0:1000)
    character(len=10) :: crd
    type(bc_type), allocatable :: bc_temp(:)

    !Initialize variables
    nmaxfaces = 3   !Maximum number of faces per cell (Initialized as minimum value here)
    nmaxcells = 8   !Maximum number of cells per node ******** Hard Coded **************
    nbc = 0         !Total number of boundary cells
    nbn = 0         !Total number of boundary nodes
    
    !Get size of grid and read node strings
    open(85,file=grdfile,iostat=ierr)
    if(ierr/=0)then
      call diag_print_error('Could not open 2dm file: ',grdfile)
    endif
    do
      read(85,*,iostat=ierr) crd
      if(ierr/=0) exit
      select case(crd)
        case('E4Q') !Quads are treated as general polygons
          ncellpoly = ncellpoly + 1 
          nmaxfaces = max(nmaxfaces,4)
        case('E3T') !Triangles are treated as general polygons
          ncellpoly = ncellpoly + 1
          nmaxfaces = max(nmaxfaces,3)
        case('E4R') !Regular Cartesian cells
          ncellsimple = ncellsimple + 1
          nmaxfaces = max(nmaxfaces,4)
        case('E5J') !Joint Cartesian cells
          ncelljoint = ncelljoint + 1     
          nmaxfaces = max(nmaxfaces,6)
        case('ND') !Nodes
          nnodes = nnodes + 1
        case('NS') !Node strings        
          nbcstr = nbcstr + 1  !Node string count   
          if(nbcstr==1)then
            allocate(bc_str(nbcstr))
          else
            allocate(BC_TEMP(nbcstr-1))
            BC_TEMP(1:nbcstr-1) = bc_str(1:nbcstr-1)
            deallocate(bc_str)
            allocate(bc_str(nbcstr))
            bc_str(1:nbcstr-1)=BC_TEMP(1:nbcstr-1)
            deallocate(BC_TEMP)
          endif
          idnstrtemp = 0
          nns = 1  !length of node string
          iend=0
          nseg=1          
          do while(idnstrtemp(iend)>=0)
            nns = nns + 1
            nmod = mod(nns-1,10)
            if(nmod/=0)then
              backspace(85) !Not new line
            else
              nseg = nseg + 1   !New line
            endif  
            istart = 10*(nseg-1) + 1
            iend = istart + nmod
            read(85,*) crd,idnstrtemp(istart:iend)
          enddo
          ncs = nns - 1
          bc_str(nbcstr)%nnodes = nns
          bc_str(nbcstr)%ncells = ncs
          allocate(bc_str(nbcstr)%nodes(nns))
          allocate(bc_str(nbcstr)%cells(bc_str(nbcstr)%ncells))
          allocate(bc_str(nbcstr)%faces(bc_str(nbcstr)%ncells))
          idnstrtemp(nns) = abs(idnstrtemp(nns))
          bc_str(nbcstr)%nodes = idnstrtemp(1:nns)
          !bc_str(nbcstr)%idnum = nbcstr
          nbn = nbn + nns    !Total number of boundary nodes
          nbc = nbc + ncs    !Total number of boundary cells
      end select  
    enddo
    close(85)
    
    !Total number of cells or elements
    ncells = ncellsimple + ncelljoint + ncellpoly
    !Total number of cells or elements including ghost cells
    ncellsD = ncells + nbc !*3 !Temporary
    
    !Load element and nodes information
    allocate(cell2node(nmaxfaces,ncellsD),ncface(ncellsD))
    cell2node = 0; ncface = 0
    allocate(xn(nnodes),yn(nnodes),zbn(nnodes))
    xn = 0.0; yn = 0.0; zbn = 0.0
    nd = 0
    open(85,file=grdfile,iostat=ierr)
    do
      read(85,*,iostat=ierr) crd
      if(ierr/=0) exit
      select case(crd)
        case('E4Q')
          backspace(85)
          read(85,*,iostat=ierr) crd,id,(cell2node(j,id),j=1,4)
          if(ierr/=0) exit
          ncface(id) = 4
        case('E3T')
          backspace(85)
          read(85,*,iostat=ierr) crd,id,(cell2node(j,id),j=1,3)
          if(ierr/=0) exit
          ncface(id) = 3
        case('ND')
          backspace(85)
          read(85,*,iostat=ierr) crd,id,xn(id),yn(id),zbn(id)
          if(ierr/=0) exit
          zbn(id)=-zbn(id) !Depth converted to elevation
          nd = nd + 1
          if(id/=nd)then
            write(*,*) id, nd 
            read(*,*) 
          endif
      end select  
    enddo
    close(85)
    
    !Allocate variables
    call assigned_bnd           !Make boundary ghost cells and connect to nodes, cell2node,ncface
    call connect_node_to_cell   !nncell,nod2cell
    call forward_cell_to_cell   !cell2cell
    call unassigned_bnd         !Adds a wall at unasigned boundaries and connects to nodes,cell2node,ncface,cell2cell
    call connect_node_to_cell   !nncell, node2cell
    call backwards_cell_to_cell !llec2llec

    !Total nodes including inactive ones (variable used for output)
    ncellsfull = nnodes
    
    !!allocate(idmap(ncellsD),mapid(ncellsD))
    !!do i=1,ncellsD
    !!  idmap(i)=i
    !!  mapid(i)=i
    !!enddo
    
    return
    end subroutine read_2dm
    
!****************************************************    
    subroutine assigned_bnd
!Node to cell connection
!****************************************************  
    use size_def, only: ncells,ncellsD,nmaxcells,nnodes
    use geo_def, only: ncface,node2cell,nncell,cell2node
    use bnd_def, only: BC_str,nbcstr
    implicit none
    integer :: nck,j,k
    
    nck=ncells
    ncface(ncells+1:ncellsD)=0
    do j=1,nbcstr
      do k=1,BC_str(j)%nnodes-1
        nck=nck+1
        ncface(nck)=ncface(nck)+1
        cell2node(ncface(nck),nck) = BC_str(j)%nodes(k)        
        cell2node(ncface(nck)+1,nck) = BC_str(j)%nodes(k+1)
        if(nck>ncells) ncface(nck)=1
      enddo
    enddo
    
    return
    end subroutine assigned_bnd
    
!****************************************************    
    subroutine connect_node_to_cell
!Node to cell connection
!****************************************************    
    use size_def, only: ncells,ncellsD,nmaxcells,nnodes
    use geo_def, only: ncface,node2cell,nncell,cell2node
    use comvarbl, only: casename,flowpath
    use out_def, only: write_node2cell
    use diag_lib
    implicit none
    integer :: numnode,i,j,nd
    character(len=200) :: filename,msg2,msg3,msg4
    
    if(allocated(node2cell)) deallocate(node2cell)
    if(allocated(nncell)) deallocate(nncell)
    allocate(node2cell(nmaxcells,nnodes),nncell(nnodes))
    node2cell = 0; nncell = 0 
    
    do i=1,ncellsD
      numnode=ncface(i)
      if(i>ncells) numnode=2
      do j=1,numnode
        nd=cell2node(j,i)          
        nncell(nd)=nncell(nd)+1
        if(nncell(nd)>nmaxcells)then
          write(msg2,*) '  Node:                           ',nd
          write(msg2,*) '  Neighboring cells:              ',nncell(nd)
          write(msg3,*) '  Maximum # of neighboring cells: ',nmaxcells
          call diag_print_error('Problem calculating node to cell connectivity',msg2,msg3,msg4)
        endif
        node2cell(nncell(nd),nd)=i
      enddo
    enddo
    
    if(write_node2cell)then
      filename = trim(flowpath) // trim(casename) // '_node2cell.txt'
      open(67,file=filename)
      write(67,'(A)')'j,(node2cell(i,j),i=1,nncell(j))'
      do j=1,nnodes
        write(67,'(I6,8(I8))') j,(node2cell(i,j),i=1,nncell(j))
      enddo
      close(67)
    endif
    
    return
    end subroutine connect_node_to_cell

!***************************************************
    subroutine forward_cell_to_cell
!***************************************************    
    use size_def, only: ncells,ncellsD,nmaxfaces
    use geo_def, only: ncface,cell2cell,cell2node,node2cell,nncell
    implicit none
    integer :: i,ii,j,jj,k,numnode,numnode2,nd,nd2,nd3,ncf
    
    if(allocated(cell2cell))then
      deallocate(cell2cell)
    endif
    allocate(cell2cell(nmaxfaces,ncellsD))   
    cell2cell=0
        
!--- Cell to cell connection ---------------
    do i=1,ncellsD
      numnode=ncface(i)
      ncf=0
      if(i>ncells) numnode=2 !Ghost cells only have 2 nodes
loopj: do j=1,numnode  !number of faces
        nd=cell2node(j,i)
        if(j<numnode)then
          nd2=cell2node(j+1,i)        
        else
          nd2=cell2node(1,i) 
        endif
        !The current face is identified by nodes nd and nd2
        !the forward connectivity is found by identifying the
        !same 2 nodes from the cells neighboring node nd
        do k=1,nncell(nd)
          ii=node2cell(k,nd)
          if(ii==i) cycle
          numnode2=ncface(ii)
          if(ii>ncells) numnode2=2
          do jj=1,numnode2
            nd3=cell2node(jj,ii)
            if(nd2==nd3)then
              ncf=ncf+1
              cell2cell(j,i)=ii
              !!exit loopj
            endif
          enddo !jj
        enddo !k
      enddo loopj !j
    enddo !i
    
    return
    end subroutine forward_cell_to_cell

!***************************************************************************
    subroutine unassigned_bnd
!***************************************************************************    
    use size_def, only: ncells,ncellsD,nmaxfaces
    use geo_def, only: ncface,cell2cell,cell2node
    use comvarbl, only: flowpath,casename
    use out_def, only: write_cell2node,write_cell2cell
    implicit none
    integer:: ndum,i,j,k,kk,ncellsDold,nck
    integer,allocatable:: cell2cell2(:,:),ncface2(:),cell2node2(:,:)
    character(len=200) :: filename
    
!--- Check for empty missing cell to cell conections and add ghost cells -----
    ndum=0
    do i=1,ncells
      do k=1,ncface(i)
        if(cell2cell(k,i)==0)then
          ndum=ndum+1
        endif
      enddo
    enddo    
    if(ndum>0)then
      ncellsDold=ncellsD
      ncellsD=ncellsD+ndum
      allocate(cell2cell2(nmaxfaces,ncellsD),ncface2(ncellsD))       
      allocate(cell2node2(nmaxfaces,ncellsD))
      cell2cell2(:,1:ncellsDold)=cell2cell(:,1:ncellsDold)
      ncface2(1:ncellsDold)=ncface(1:ncellsDold)
      cell2node2(:,1:ncellsDold)=cell2node(:,1:ncellsDold)
      deallocate(cell2cell,ncface,cell2node)      
      allocate(cell2cell(nmaxfaces,ncellsD),ncface(ncellsD))       
      allocate(cell2node(nmaxfaces,ncellsD))
      cell2cell(:,1:ncellsDold)=cell2cell2(:,1:ncellsDold)
      ncface(1:ncellsDold)=ncface2(1:ncellsDold)
      cell2node(:,1:ncellsDold)=cell2node2(:,1:ncellsDold)
      deallocate(cell2cell2,ncface2,cell2node2)
      
      nck=ncellsDold
      do i=1,ncells
        do k=1,ncface(i)
           if(cell2cell(k,i)==0)then
            nck=nck+1
            ncface(nck)=1
            cell2cell(k,i)=nck
            cell2cell(1,nck)=i
            cell2cell(2:,nck)=0
            cell2node(1,nck)=cell2node(k,i)
            kk=k+1
            if(kk>ncface(i)) kk=1
            cell2node(2,nck)=cell2node(kk,i)
          endif
        enddo
      enddo
    endif
    
    if(write_cell2node)then
      filename = trim(flowpath) // trim(casename) // '_cell2node.txt'
      open(67,file=filename)
      write(67,'(A)')'i,(cell2node(j,i),j=1,ncface(i))'
      do i=1,ncells
        write(67,*) i,(cell2node(j,i),j=1,ncface(i))    
      enddo
      do i=ncells+1,ncellsD
        write(67,*) i,(cell2node(j,i),j=1,2)
      enddo
      close(67)
    endif
    
    if(write_cell2cell)then
      filename = trim(flowpath) // trim(casename) // '_cell2cell.txt'
      open(67,file=filename)
      write(67,'(A)') 'i,(cell2cell(j,i),j=1,ncface(i))'
      do i=1,ncellsD
        write(67,*) i,(cell2cell(j,i),j=1,ncface(i))    
      enddo
      close(67)    
    endif
    
    return
    end subroutine unassigned_bnd

!***********************************************************
    subroutine backwards_cell_to_cell
!***********************************************************    
    use size_def, only: ncells,ncellsD,nmaxfaces
    use geo_def, only: ncface,llec2llec,cell2cell
    use diag_def
    use diag_lib
    implicit none
    integer :: i,k,kk,nck,nckk
    
!--- Backward cell to cell connectivity ------------------------------------------
    allocate(llec2llec(nmaxfaces,ncellsD))  !Backward connectivity  
    llec2llec=0
    do i=1,ncellsD      !Derive llec2llec from cell2cell 
      do k=1,ncface(i)
        nck=cell2cell(k,i)        
        llec2llec(k,i)=0
!        if(nck<=ncells)then
          do kk=1,ncface(nck)
            nckk=cell2cell(kk,nck)
            if(nckk==i) llec2llec(k,i)=kk
          enddo
!        endif
        if(llec2llec(k,i)==0)then
          write(msg2,*) '  Cell: ',i
          write(msg3,*) '  Neighbor: ',cell2cell(k,i)
          call diag_print_error('Problem calculating backwards cell to cell connectivity',msg2,msg3)
        endif
      enddo
    enddo
    
    return
    end subroutine backwards_cell_to_cell
    
!*************************************************************
    subroutine geo_init()
! Initializes the Geospatial variables
! Author: Alex Sanchez, USACE-CHL
!*************************************************************
    use geo_def
    implicit none
    
    if(igridtype<=1)then
      call geo_init_cart     !Geospatial variable initialization
    else
      call geo_init_poly     !Geospatial variable initialization  
    endif 
    !if(n2Dor3D==3) then  
    !   call geometry3D      !For 3D
    !   call allocate_fl3D
    !   call geo_var3D
    !   if(noptset==3) call allocate_wavestress3D  !For 3D   
    !endif
    
    if(bathydata%ison) call geo_init_bathydata
    
    return
    end subroutine geo_init

!*********************************************************************
    subroutine geo_init_bathydata()
! Initializes the depth update variable
! Author: Alex Sanchez, USACE-CHL
!*********************************************************************    
#include "CMS_cpp.h"
    use size_def, only: ncellsD
    use geo_def
    use diag_lib, only: diag_print_error
    use prec_def
#ifdef XMDF_IO
    use in_xmdf_lib, only: read_dataseth5,readscalsteph5
#endif
    implicit none
    
    integer :: ierr
    
!#ifdef XMDF_IO
!    !call readscalh5(bathydata%file,bathydata%path,hardzb,ierr)
!    call read_dataseth5(bathydata%file,bathydata%path,'Times',bathydata%ntimes,bathydata%times)
!#else
!    call diag_print_error('Cannot read depth datasets from *.h5 file without XMDF libraries')
!#endif
!    bathydata%timen = bathydata%times(bathydata%ntimes)
!    deallocate(bathydata%times)

    !allocate(bathydata%depth0(ncellsD))
    if(.not.allocated(zb0))then
      allocate(zb0(ncellsD))
    endif
    zb0 = zb !Initialize
    allocate(bathydata%depth1(ncellsD))
    allocate(bathydata%depth2(ncellsD))
    bathydata%depth1 = -zb !Initialize
    bathydata%depth2 = -zb !Initialize

    bathydata%inc2 = 1
#ifdef XMDF_IO
    call readscalsteph5(bathydata%file,bathydata%path,bathydata%inc2,bathydata%time1,bathydata%depth1,ierr)
#else
    call diag_print_error('Cannot read depth update datasets from *.h5 file without XMDF libraries')
#endif
    bathydata%timen = huge(1.0)
    bathydata%time1 = bathydata%time1/3600.0
    bathydata%time0 = bathydata%time1
    bathydata%tjulday = bathydata%time1/24.0 !Convert from hours to days **** Hard coded for now *****     
    bathydata%time1 = 0.0 !Assumes first time is at 0.0 *** Hard coded for now *****
    bathydata%time2 = 0.0
    bathydata%depth2 = bathydata%depth1

    call geo_bathy_update(zb,-1)
    
    return
    end subroutine geo_init_bathydata
    
!************************************************************
    subroutine geo_bathy_update(zbed,msignbed)
! Updates the grid bathymetry (bed elevations) by interpolating
! from a temporally variable dataset
! Author: Alex Sanchez, USACE-CHL
!************************************************************
#include "CMS_cpp.h"
    use size_def, only: ncells,ncellsD
    use geo_def
    use prec_def
    use diag_lib, only: diag_print_error
    use comvarbl, only: timehrs,ramp
#ifdef XMDF_IO    
    use in_xmdf_lib, only: read_dataseth5,readscalsteph5
#endif
    use bnd_def, only: nbndstr,bnd_str
    use prec_def
    implicit none
    integer :: i,inc1,ierr,msignbed
    integer :: ibnd,j,k,nck
    real(ikind):: fac,fac0,fac1,fac2
    real(ikind):: zbed(ncellsD)
    !integer :: np
    !integer, parameter :: nb = 2 !>=np+1
    !real(ikind) :: lb(nb)
    
    !nti = 1 !Input polynomial order
    !call plagr_fit(bathydata%ntimes,bathydata%times,timehrs,nb,lb,nti,np,bathydata%inc2)  
    !yi = sum(lb(1:np+1)*y(k:k+np))
    
    !call interp_linear_fit(bathydata%ntimes,bathydata%times,timehrs,inc,fac1,fac2)
    
    do while(bathydata%time1>timehrs)
      bathydata%inc2 = bathydata%inc2 - 1
      inc1 = bathydata%inc2 - 1 
      bathydata%time2 = bathydata%time1
      bathydata%depth2 = bathydata%depth1
#ifdef XMDF_IO
      call readscalsteph5(bathydata%file,bathydata%path,inc1,bathydata%time1,bathydata%depth1,ierr)
      if(ierr/=0)then
        exit
      endif
#else
      call diag_print_error('Cannot read depth update datasets from *.h5 file without XMDF libraries')
#endif
      bathydata%time1 = bathydata%time1/3600.0
      bathydata%time1 = bathydata%time1 - bathydata%time0
      continue
    enddo
    
    do while(bathydata%time2<=timehrs .and. bathydata%time2<bathydata%timen)
      bathydata%inc2 = bathydata%inc2 + 1
      bathydata%time1 = bathydata%time2
      bathydata%depth1 = bathydata%depth2
#ifdef XMDF_IO
      call readscalsteph5(bathydata%file,bathydata%path,bathydata%inc2,bathydata%time2,bathydata%depth2,ierr)
      if(ierr==3)then !Reached in of time series
        bathydata%ntimes = bathydata%inc2 - 1
        bathydata%timen = bathydata%time2
        exit
      endif
#else
      call diag_print_error('Cannot read depth update datasets from *.h5 file without XMDF libraries')
#endif
      bathydata%time2 = bathydata%time2/3600.0
      bathydata%time2 = bathydata%time2 - bathydata%time0
      !call read_dataseth5(bathydata%file,bathydata%path,'Values',bathydata%ntimes,bathydata%depth2)
      continue
    enddo
    
    !Temporal interpolation
    fac = (timehrs-bathydata%time1)/(bathydata%time2-bathydata%time1)
    fac = max(min(fac,1.0),0.0) !Avoids extrapolation
    fac0 = -msignbed*(1.0-ramp)  !Note sign
    fac1 = msignbed*ramp*(1.0-fac)  !Apply ramp here
    fac2 = msignbed*ramp*fac        !Apply ramp here
!$OMP PARALLEL DO PRIVATE(i)
    do i=1,ncells
      zbed(i) = fac0*zb0(i) + fac1*bathydata%depth1(i) + fac2*bathydata%depth2(i)
    enddo
!$OMP END PARALLEL DO
    
    !Copy to ghost cells
    do ibnd=1,nbndstr
      do j=1,bnd_str(ibnd)%ncells
        i=bnd_str(ibnd)%cells(j)
        k=bnd_str(ibnd)%faces(j)
        nck=cell2cell(k,i)
        zbed(nck)=zbed(i)
      enddo
    enddo

    return
    end subroutine geo_bathy_update
    
!***********************************************************
    subroutine geo_flowdepth_update
!***********************************************************
    use size_def
    use bnd_def
    use flow_def
    use geo_def
    use der_def
    use der_lib
    use interp_lib
    use prec_def
    implicit none
    integer :: i,j,k,ibnd,nck
    real(ikind) :: val
    
    !=== Bed-slopes ====--------------===========
    call der_grad_eval(goa,0,zb,dzbx,dzby) !Bed-slope

    !=== Update bed elevation at cell faces =====
    call interp_scal_cell2face(zb,0,zbk,dzbx,dzby)
    
!$OMP PARALLEL DO PRIVATE(i,val)
    do i=1,ncells
      val = h(i)
      h(i) = max(hdry,p(i)*gravinv-zb(i))  !New total water depth
      u(i) = val*u(i)/h(i)
      v(i) = val*v(i)/h(i)
    enddo
!$OMP END PARALLEL DO
    
    !Copy to ghost cells
    do ibnd=1,nbndstr
      do j=1,bnd_str(ibnd)%ncells
        i=bnd_str(ibnd)%cells(j)
        k=bnd_str(ibnd)%faces(j)
        nck=cell2cell(k,i)
        h(nck)=max(hmin,p(nck)*gravinv-zb(nck))
      enddo
    enddo

    return
    end subroutine geo_flowdepth_update
    
!*************************************************************     
    subroutine geo_init_cart()
! Calculates geometrical quantities
!
! Input:
!  ncells, ncellsD, x, y, dx, dy, 
!  cell2cell, ncface, idirface,
!
! Output:
!  dxymin, dsxy, fintp, areap, areaavg, llec2llec, ds, 
!  ncellsimple, ncelljoint, idcellsimple, idcelljoint
!
! written by Weiming Wu, NCCHE
! modified by Alex Sanchez, USACE-CHL
!*************************************************************         
#include "CMS_cpp.h"
    use size_def
    use geo_def    
    use flow_def
    use struct_def
    use comvarbl
    use const_def, only: deg2rad
    use diag_def
    use diag_lib
    use prec_def
    implicit none
    integer :: i,k,kk,nck,nckk,ncelljointtemp,ncellsimpletemp  !j, jcn, jjcn
    !integer :: inorth,isouth,ieast,iwest,icell,icnt
#ifdef DIAG_MODE
    integer :: ii
#endif
    integer,allocatable :: indictface(:,:),isjoint(:) !Temporary arrays
    real(ikind) :: costheta,sintheta
    
    if(projfl%iHorizUnits/=2)then
      call diag_print_error('Input horizontal units must be meters')
    endif
    
!--- Cell areas ------------------------------------------------------
    allocate(areap(ncellsD))
    do i=1,ncellsD
      areap(i)=dx(i)*dy(i)
    enddo
    areaavg=sum(areap(1:ncells))/real(ncells)
    areamax=maxval(areap(1:ncells))
    areamin=minval(areap(1:ncells))

!--- Backwards connectivity ------------------------------------------
    allocate(llec2llec(nmaxfaces,ncellsD))
    llec2llec=0
    do i=1,ncells      !Derive llec2llec from cell2cell 
      do k=1,ncface(i)
        nck=cell2cell(k,i)
        if(nck <= ncells)then
          do kk=1,ncface(nck)
            nckk=cell2cell(kk,nck)
            if(nckk == i)then
              llec2llec(k,i) = kk
              llec2llec(k,i) = kk     !I have no idea why, but without some re-evaluation of this variable, my machine locks.  A simple assignment is faster than an IF that does nothing.  Mitch
            endif
          enddo
        else
          if(idirface(k,i) == 1) llec2llec(k,i)=3
          if(idirface(k,i) == 3) llec2llec(k,i)=1
          if(idirface(k,i) == 2) llec2llec(k,i)=4
          if(idirface(k,i) == 4) llec2llec(k,i)=2
        endif
        if(llec2llec(k,i)==0)then
          call diag_print_error('Problem in calculating backwards cell to cell connectivity')
        endif
      enddo
    enddo

!--- Face lengths ------------------------------------------
    allocate(ds(nmaxfaces,ncellsD))
    ds=0.0
    do i=1,ncells
      do k=1,ncface(i)
        nck=cell2cell(k,i)
        if(idirface(k,i)==1 .or. idirface(k,i)==3) ds(k,i)=min(dx(i),dx(nck))
        if(idirface(k,i)==2 .or. idirface(k,i)==4) ds(k,i)=min(dy(i),dy(nck))
        if(ds(k,i)<1.0e-20)then
          write(msg2,*) '  Cell: ',mapid(i)
          write(msg3,*) '  Face: ',k
          write(msg4,*) '  Neighbor: ',mapid(nck)
          call diag_print_error('Problem in calculating cell face length at: ',msg2,msg3,msg4)
        endif
      enddo
    enddo

!--- distance from cell center to cell center ------- used in avalanching and diffusion terms
    allocate(dc(nmaxfaces,ncellsD))
    dc=0.0
    do i=1,ncells !Alex, July 20, 2009
      do k=1,ncface(i)
         nck=cell2cell(k,i)
         dc(k,i)=sqrt((x(i)-x(nck))**2+(y(i)-y(nck))**2)                  
      enddo
    enddo

!--- Group Cells as simple and joint -------------------------------  
    allocate(isjoint(ncellsD))
    isjoint=0 !Initialize
    do i=1,ncells
      do k=1,ncface(i)
        nck=cell2cell(k,i)
        if(nck>ncells) cycle !All ghost/dummy cells are regular cells
        if(idirface(k,i)==2 .or. idirface(k,i)==4)then !X-direction
          if(dy(i)/dy(nck)>1.2 .or. dy(i)/dy(nck)<0.8)then !Joint cell
            isjoint(i)=1; isjoint(nck)=1
          endif
        else  !Y-direction
          if(dx(i)/dx(nck)>1.2 .or. dx(i)/dx(nck)<0.8)then !Joint cell
            isjoint(i)=1; isjoint(nck)=1
          endif
        endif
      enddo
    enddo
    ncelljoint=sum(isjoint(1:ncells))
    ncellsimple=ncells-ncelljoint
    ncelljointtemp=ncelljoint   !Save
    ncellsimpletemp=ncellsimple !Save
    allocate(idcellsimple(ncellsimple),idcelljoint(ncelljoint))
    ncellsimple=0; ncelljoint=0  
    do i=1,ncells
      if(isjoint(i)==1)then
        ncelljoint=ncelljoint+1
        idcelljoint(ncelljoint)=i
      else
        ncellsimple=ncellsimple+1
        idcellsimple(ncellsimple)=i  
      endif
    enddo
    if(ncellsimple+ncelljoint/=ncells .or. &
       ncellsimple/=ncellsimpletemp .or. &
       ncelljoint/=ncelljointtemp)then
      write(msg2,*) '  Simple cells: ',ncellsimple 
      write(msg3,*) '  Joint cells: ',ncelljoint
      write(msg4,*) '  Total cells: ',ncells
      call diag_print_error('Problem simple and joint cells',msg2,msg3,msg4)
    endif
    
!--- Calculate no-repeat cell connectivity ---- 
    allocate(indictface(nmaxfaces,ncellsD))
    indictface=1 !indicator for visiting cell face (0 = avoid repeating)
    allocate(nxyface(ncellsD),kxyface(nmaxfaces,ncellsD))
    allocate(nxface(ncellsD),kxface(ndmaxfaces,ncellsD))
    allocate(nyface(ncellsD),kyface(ndmaxfaces,ncellsD))
    nxyface=0; kxyface=0
    nxface=0; kxface=0
    nyface=0; kyface=0
    do i=1,ncells
      nxface(i)=0; nyface(i)=0; nxyface(i)=0
      do k=1,ncface(i)
        if(indictface(k,i)==1)then
          nck=cell2cell(k,i)
          indictface(k,i)=0
          indictface(llec2llec(k,i),nck)=0
          nxyface(i)=nxyface(i)+1
          kxyface(nxyface(i),i)=k
          if(idirface(k,i)==2 .or. idirface(k,i)==4)then !x-direction
            nxface(i)=nxface(i)+1
            kxface(nxface(i),i)=k
          else                                               !y-direction
            nyface(i)=nyface(i)+1
            kyface(nyface(i),i)=k
          endif
        endif  
      enddo
    enddo 
    
!--- Vectors ------------------------------
    allocate(fnx(nmaxfaces,ncellsD),fny(nmaxfaces,ncellsD))
    allocate(dn(nmaxfaces,ncellsD))
    allocate(rx(nmaxfaces,ncellsD),ry(nmaxfaces,ncellsD))
    do i=1,ncells
      do k=1,ncface(i)
        nck=cell2cell(k,i)
        select case(idirface(k,i))
          case(1) !North
            fnx(k,i)= 0.0
            fny(k,i)= 1.0 
            dn(k,i)=abs(y(nck)-y(i)) !Normal
            rx(k,i)=x(nck)-x(i)
            ry(k,i)=dy(i)/2.0
          case(2) !East
            fnx(k,i)= 1.0
            fny(k,i)= 0.0         
            dn(k,i)=abs(x(nck)-x(i)) !Normal
            rx(k,i)=dx(i)/2.0 
            ry(k,i)=y(nck)-y(i)
          case(3) !South
            fnx(k,i)= 0.0
            fny(k,i)=-1.0 
            dn(k,i)=abs(y(nck)-y(i)) !Normal
            rx(k,i)=x(nck)-x(i)
            ry(k,i)=-dy(i)/2.0
          case(4) !West
            fnx(k,i)=-1.0
            fny(k,i)= 0.0 
            dn(k,i)=abs(x(nck)-x(i)) !Normal
            rx(k,i)=-dx(i)/2.0 
            ry(k,i)=y(nck)-y(i)
          case default
            call diag_print_error('Problem calculating unit vectors for Cartesian grid')
        end select
      enddo  
    enddo
    
!--- General Geometry variables -------------------------    
    allocate(rpx(nmaxfaces,ncellsD),rpy(nmaxfaces,ncellsD))
    allocate(rnx(nmaxfaces,ncellsD),rny(nmaxfaces,ncellsD))
    allocate(dnx(nmaxfaces,ncellsD),dny(nmaxfaces,ncellsD))    
    allocate(dsx(nmaxfaces,ncellsD),dsy(nmaxfaces,ncellsD))
    allocate(dsxy(nmaxfaces,ncellsD))
    do i=1,ncells
      do k=1,ncface(i)
        !Normal centroid to face vectors
        rnx(k,i)=fnx(k,i)*abs(rx(k,i))
        rny(k,i)=fny(k,i)*abs(ry(k,i))
        !Normal centroid to centroid vectors
        dnx(k,i)=fnx(k,i)*dn(k,i)
        dny(k,i)=fny(k,i)*dn(k,i)
        !Perpendicular centroid to face vectors
        rpx(k,i)=rx(k,i)-rnx(k,i)
        rpy(k,i)=ry(k,i)-rny(k,i)        
        !Face/edge vectors
        dsx(k,i)=fnx(k,i)*ds(k,i)
        dsy(k,i)=fny(k,i)*ds(k,i)
        !Diffusion term variable
        dsxy(k,i)=ds(k,i)/dn(k,i)
      enddo
    enddo
    
!--- Global Cartesian Coordinates ----------------------
    allocate(xc(ncellsD),yc(ncellsD))
    costheta = cos((azimuth_fl)*deg2rad)  
    sintheta = sin((azimuth_fl)*deg2rad)
    do i=1,ncellsD 
      xc(i)=xOrigin+x(i)*costheta-y(i)*sintheta
      yc(i)=yOrigin+x(i)*sintheta+y(i)*costheta
    enddo

!--- Projection ----------------------------------------------
    !if(len_trim(hproj)==0) call get_coord_info     
    
    !Construct strings
757 format('"',A,', ',A,', ZONE ',I4,', ',A,'"')
151 format('"',A,', ',A,'"')    
    write(HPROJ,757)trim(aHorizDatum(projfl%iHorizDatum)),&
        trim(aHorizCoordSystem(projfl%iHorizCoordSystem)),&
        projfl%iHorizZone,&
        trim(aHorizUnits(projfl%iHorizUnits))
    write(VPROJ,151) trim(aVertDatum(projfl%iVertDatum)),&
        trim(aVertUnits(projfl%iVertUnits))

!--- Bed Slope ----------------------
    allocate(dzbx(ncellsD),dzby(ncellsD))
    dzbx=0.0; dzby=0.0
    
#ifdef DEV_MODE
!--- Second Upwind Connectivity -----------------------------------
!The variable is setup so that we have:
!                              Downstream(d)  First-Upwind(c)     Second-Upwind(u)
! Outflow (flux(k,i)>0.0)    cell2cell(k,i)        i            cell2upwdcell(i,k)
! Inflow  (flux(k,i)<0.0)          i          cell2cell(k,i)   cell2upwdcell(nck,jcn)
!The variable is only active (>0) at faces where d,c,u are inline with each other
    allocate(cell2upwdcell(ncellsD,nmaxfaces))
    cell2upwdcell = 0
    do i=1,ncells
      do k=1,ncface(i)
        nck=cell2cell(k,i)  !Foward connectivity
        jcn=llec2llec(k,i)  !Backwards connectivity  
        if(abs(rpx(k,i))+abs(rpy(k,i))>1.0e-6 .or. &
           abs(rpx(jcn,nck))+abs(rpy(jcn,nck))>1.0e-6) cycle !skip refinement faces
        !Find opposing face
        do kk=1,ncface(i)
          if(kk==k) cycle !Skip if same face
          if(kkface(idirface(kk,i))/=idirface(k,i)) cycle !skip if not opposing face
          !if(abs(rpx(nck,kk))+abs(rpy(nck,kk))>1.0e-6) cycle !skip if refinement face
          nckk=cell2cell(kk,i)  !Foward connectivity
          jjcn=llec2llec(kk,i)  !Backwards connectivity  
          if(abs(rpx(kk,i))+abs(rpy(kk,i))>1.0e-6 .or. &
           abs(rpx(jjcn,nckk))+abs(rpy(jjcn,nckk))>1.0e-6) cycle !skip refinement faces
          cell2upwdcell(i,k) = nckk !Found apposing nonrefined neighbor cell
          exit          
        enddo
      enddo
    enddo
#endif    
    
#ifdef DIAG_MODE
!---- Check for joint cells ----
    isjoint = 0
    indictface = 0
    do ii=1,ncelljoint
      i=idcelljoint(ii)
      isjoint(i) = 1
      do k=1,ncface(i)
        nck=cell2cell(k,i) !Foward connectivity
        isjoint(nck) = 1
        indictface(k,i) = 1
      enddo
    enddo
    do ii=1,ncelljoint
      i=idcelljoint(ii)  
      do j=1,nxyface(i) !no repeat for all sides
        k=kxyface(j,i)
        nck=cell2cell(k,i) !Foward connectivity
        if(nck>ncells) cycle
        if(isjoint(nck)==0 .or. indictface(k,i)==0)then
          write(msg2,*) '  Cell = ',mapid(i)
          write(msg3,*) '  Neighbor = ',mapid(nck)
          write(msg4,*) '  isjoint(i)   = ',isjoint(i)
          write(msg5,*) '  isjoint(nck) = ',isjoint(nck)
          write(msg6,*) '  indictface(1:ncface(i),i) = ',indictface(1:ncface(i),i)          
          call diag_print_error('Problem calculating cell to cell connectivity',&
            '  Inconsistent neighboring cells ',msg2,msg3,msg4,msg5,msg6)
        endif 
      enddo
      do j=1,nxface(i) !no repeat east and west sides
        k=kxface(j,i)
        nck=cell2cell(k,i) !Foward connectivity
        if(nck>ncells) cycle
        if(isjoint(nck)==0 .or. indictface(k,i)==0)then
          write(msg2,*) '  Cell = ',mapid(i)
          write(msg3,*) '  Neighbor = ',mapid(nck)
          write(msg4,*) '  isjoint(i)   = ',isjoint(i)
          write(msg5,*) '  isjoint(nck) = ',isjoint(nck)
          write(msg6,*) '  indictface(1:ncface(i),i) = ',indictface(1:ncface(i),i)            
          call diag_print_error('Problem calculating cell to cell connectivity',&
            '  Inconsistent neighboring cells ',msg2,msg3,msg4,msg5,msg6)
        endif 
      enddo
      do j=1,nyface(i) !no repeat north and south sides
        k=kyface(j,i)
        nck=cell2cell(k,i) !Foward connectivity
        if(nck>ncells) cycle
        if(isjoint(nck)==0 .or. indictface(k,i)==0)then
          write(msg2,*) '  Cell = ',mapid(i)
          write(msg3,*) '  Neighbor = ',mapid(nck)
          write(msg4,*) '  isjoint(i)   = ',isjoint(i)
          write(msg5,*) '  isjoint(nck) = ',isjoint(nck)
          write(msg6,*) '  indictface(1:ncface(i),i) = ',indictface(1:ncface(i),i)           
          call diag_print_error('Problem calculating cell to cell connectivity',&
            '  Inconsistent neighboring cells ',msg2,msg3,msg4,msg5,msg6)
        endif 
      enddo
    enddo
    
!--- Check No-repeat cell connectivity -----
    !XY No-Repeat
    indictface=1 !indicator for visiting cell face (0 = avoid repeating)
    do i=1,ncells
      do j=1,nxyface(i) !No Repeat of cell faces
        k=kxyface(j,i)
        nck=cell2cell(k,i) !Forward connectivity 
        jcn=llec2llec(k,i) !Backward connectivity
        indictface(k,i)=0
        indictface(jcn,nck)=0
      enddo
    enddo
    do i=1,ncells
      do k=1,ncface(i)
        if(indictface(k,i)==1)then
          call diag_print_error('Problem calculating no-repeat cell connetivity')
        endif
      enddo
    enddo
    
    !X No-Repeat
    indictface=1 !indicator for visiting cell face (0 = avoid repeating)
    do i=1,ncells
      do j=1,nxface(i) !No Repeat of cell faces
        k=kxface(j,i)
        nck=cell2cell(k,i) !Forward connectivity 
        jcn=llec2llec(k,i) !Backward connectivity
        indictface(k,i)=0
        indictface(jcn,nck)=0
      enddo
    enddo
    do i=1,ncells
      do k=1,ncface(i)
        if(idirface(k,i)==2 .or. idirface(k,i)==4)then !X-direction        
          if(indictface(k,i)==1)then        
            call diag_print_error('Problem calculating X-direction no-repeat cell connetivity')
          endif
        endif
      enddo
    enddo
    
    !Y No-Repeat
    indictface=1 !indicator for visiting cell face (0 = avoid repeating)
    do i=1,ncells
      do j=1,nyface(i) !No Repeat of cell faces
        k=kyface(j,i)
        nck=cell2cell(k,i) !Forward connectivity 
        jcn=llec2llec(k,i) !Backward connectivity
        indictface(k,i)=0
        indictface(jcn,nck)=0
      enddo
    enddo
    do i=1,ncells
      do k=1,ncface(i)
        if(idirface(k,i)==1 .or. idirface(k,i)==3)then !Y-direction
          if(indictface(k,i)==1)then          
            call diag_print_error('Problem calculating Y-direction no-repeat cell connetivity')
          endif
        endif
      enddo
    enddo
    
!!--- Check face count ---    
!    indictface = 1.0
!    do i=1,ncells
!      do j=1,nxyface(i) !No repeat faces
!        k=kxyface(j,i)
!        nck=cell2cell(k,i) !Forward connectivity 
!        jcn=llec2llec(k,i) !Backward connectivity
!        indictface(k,i)=0
!        indictface(jcn,nck)=0
!      enddo
!    enddo    
!    do j=1,W_str%ncells
!      i=W_str%cells(j)
!      k=W_str%faces(j)
!      if(indictface(k,i)>0.5)then
!        call diag_print_error('Problem with face count')
!      endif
!    enddo    
    
#endif

    deallocate(indictface,isjoint)
    
    return
    end subroutine geo_init_cart
    
!***************************************************************************************
    subroutine geo_init_poly()
! Calculates geometrical quantities
! Author: Alex Sanchez
! Last Modified: 01/22/14
!***************************************************************************************
#include "CMS_cpp.h"
    use size_def
    use geo_def
    use geo_lib, only: polyvar,inpoly
    use comvarbl, only: flowpath,casename
    use out_def, only: write_geocells
    use diag_def
    use diag_lib
    use prec_def
    implicit none
    integer :: i,ii,j,k,kk,nck,jcn,nd,nck2
    real(ikind) :: dof,dif
    !real(ikind) :: xr,yr
    integer,dimension(nmaxfaces,ncellsD) :: indictface
    real(ikind),dimension(nmaxfaces) :: xni,yni,fnxi,fnyi,dsi,dsxi,dsyi
    real(ikind),dimension(nmaxfaces) :: rxi,ryi,rnxi,rnyi,rpxi,rpyi
    character(len=200) :: filename
    
!--- Geometry variables ----------------
    allocate(zb(ncellsD),zbk(nmaxfaces,ncellsD))
    allocate(dzbx(ncellsD),dzby(ncellsD))    
    allocate(x(ncellsD),y(ncellsD),areap(ncellsD))
    allocate(fnx(nmaxfaces,ncellsD),fny(nmaxfaces,ncellsD))
    allocate(ds(nmaxfaces,ncellsD),dsx(nmaxfaces,ncellsD),dsy(nmaxfaces,ncellsD))     
    allocate(rx(nmaxfaces,ncellsD),ry(nmaxfaces,ncellsD))    
    allocate(rnx(nmaxfaces,ncellsD),rny(nmaxfaces,ncellsD))
    allocate(rpx(nmaxfaces,ncellsD),rpy(nmaxfaces,ncellsD))    
    zb=0.0; zbk=0.0
    dzbx=0.0; dzby=0.0
    x=0.0; y=0.0; areap=0.0    
    fnx = 0.0; fny = 0.0    
    ds=0.0; dsx = 0.0; dsy = 0.0
    rx=0.0;   ry=0.0
    rnx=0.0; rny=0.0  
    rpx=0.0; rpy=0.0  
    do i=1,ncells
      zb(i)=sum(zbn(cell2node(1:ncface(i),i)))/real(ncface(i),kind=ikind) !Set as average of nodes ******* TEMPORARY **********
      !Copy information into temporary variables
      do k=1,ncface(i)
        xni(k)=xn(cell2node(k,i)); yni(k)=yn(cell2node(k,i))        
      enddo
      !Calculate polygonal variables
      call polyvar(nmaxfaces,ncface(i),xni,yni,&
           areap(i),x(i),y(i),fnxi,fnyi,dsi,dsxi,dsyi,&
           rxi,ryi,rnxi,rnyi,rpxi,rpyi)
      !Copy variables back from temporary variables
      do k=1,ncface(i)
        xni(k)=xn(cell2node(k,i)); yni(k)=yn(cell2node(k,i))        
        fnx(k,i)=fnxi(k); fny(k,i)=fnyi(k)
        ds(k,i)=dsi(k); dsx(k,i)=dsxi(k); dsy(k,i)=dsyi(k)        
        rx(k,i)=rxi(k); ry(k,i)=ryi(k)        
        rnx(k,i)=rnxi(k); rny(k,i)=rnyi(k)
        rpx(k,i)=rpxi(k); rpy(k,i)=rpyi(k)
        !!!Check if recontruction point is within polygon
        !!xr = x(i)+rpx(k,i); yr = y(i)+rpy(k,i)
        !!if(.not.inpoly(xr,yr,nmaxfaces,ncface(i),xni,yni))then
        !!    write(*,*)
        !!    write(msg,*)  '  Cell: ',i
        !!    write(msg2,*) '  Face: ',k
        !!    write(msg3,*) '  fnx(k,i),fny(k,i): ',fnx(k,i),fny(k,i)
        !!    write(msg4,*) '  x(i),y(i): ',x(i),y(i)
        !!    write(msg5,*) '  xn(:): ', xni(1:ncface(i))
        !!    write(msg6,*) '  yn(:): ', yni(1:ncface(i))
        !!    write(msg7,*) '  xr,yr: ',xr,yr
        !!    call diag_print_warning('Reconstruction point outside of cell',&
        !!       msg,msg2,msg3,msg4,msg5,msg6,'Moving point to cell boundary')
        !!    continue   
        !!endif
      enddo
    enddo    

!-- Ghost cell geometry -----
    do i=ncells+1,ncellsD
      nck=cell2cell(1,i)
      jcn=llec2llec(1,i)
      areap(i) = areap(nck)
      zb(i) = zb(nck)
      ds(1,i) = ds(jcn,nck)
      x(i) = x(nck) + 2.0*rnx(jcn,nck)
      y(i) = y(nck) + 2.0*rny(jcn,nck)
      !rx(1,i) = -rx(jcn,nck)
      !ry(1,i) = -ry(jcn,nck)
      rx(1,i) = x(nck) + rx(jcn,nck) - x(i)
      ry(1,i) = y(nck) + ry(jcn,nck) - y(i)
      fnx(1,i) = -fnx(jcn,nck)
      fny(1,i) = -fny(jcn,nck)
      dsx(1,i) = -dsx(jcn,nck)
      dsy(1,i) = -dsy(jcn,nck)
      rnx(1,i) = -rnx(jcn,nck)
      rny(1,i) = -rny(jcn,nck)
      !rpx(1,i) = 0.0
      !rpy(1,i) = 0.0
      rpx(1,i) = -rpx(jcn,nck)
      rpy(1,i) = -rpy(jcn,nck)
    enddo    

!--- Average cell area --------------------
    areaavg=sum(areap(1:ncells))/real(ncells)
    areamax=maxval(areap(1:ncells))
    areamin=minval(areap(1:ncells))

!--- Interpolation, diffusion and cell-center distances ------------
    allocate(dc(nmaxfaces,ncellsD),dsxy(nmaxfaces,ncellsD))
    allocate(dn(nmaxfaces,ncellsD))
    allocate(dnx(nmaxfaces,ncellsD),dny(nmaxfaces,ncellsD))
    dc=0.0; dsxy=0.0
    dn=0.0
    dnx=0.0; dny=0.0
    do i=1,ncellsD
      do k=1,ncface(i)
        nck=cell2cell(k,i)
        jcn=llec2llec(k,i)
        !Normal distance between i and nck
        dof=sqrt(rnx(k,i)**2+rny(k,i)**2)        
        dif=sqrt(rnx(jcn,nck)**2+rny(jcn,nck)**2)
        dn(k,i)=dof+dif             
        !Normal centroid to centroid vectors
        dnx(k,i)=fnx(k,i)*dn(k,i)
        dny(k,i)=fny(k,i)*dn(k,i)
        !Used in diffusion terms
        dsxy(k,i)=ds(k,i)/dn(k,i)
        !Cell-to-cell distance
        dc(k,i)=sqrt((x(i)-x(nck))**2+(y(i)-y(nck))**2)
      enddo
    enddo   
    
!--- No Grouping of Cells -------------------------------  
    ncellpoly=ncells 
    
!--- No-Repeat (face/cell) connectivity ----       
    allocate(nxyface(ncellsD),kxyface(nmaxfaces,ncellsD))           
    nxyface=0; kxyface=0
    indictface=1 !indicator for visiting cell face (0 = avoid repeating)
    do i=1,ncells
      nxyface(i)=0
      do k=1,ncface(i)
        nck=cell2cell(k,i)
        if(indictface(k,i)==1)then          
          jcn=llec2llec(k,i)  
          indictface(jcn,nck)=0
          nxyface(i)=nxyface(i)+1
          kxyface(nxyface(i),i)=k
        endif  
      enddo
    enddo 
    
!---- Extended computational stencil -------------------------
    nmaxextcells = 12
    allocate(cell2extcell(nmaxextcells,ncellsD))
    allocate(ncnode(ncellsD))
    ncnode = 0
    cell2extcell = 0
di: do i=1,ncellsD
  dj: do j=1,ncface(i)
        nd = cell2node(j,i)
   dok: do k=1,nncell(nd)
          nck = node2cell(k,nd)
          !if(nck==i .or. nck>ncells .or. nck==0) cycle dok
          if(nck==i .or. nck==0) cycle dok
          !Search list
          ii = 0
    dokk: do kk=1,ncface(i)
            nck2 = cell2cell(kk,i)
            !if(nck2==i .or. nck2>ncells) cycle dokk
            if(nck2==i) cycle dokk
            if(nck==nck2)then
              ii = 1
              cycle dok !Found match
            endif
          enddo dokk
          !nck not in list cell2cell
          if(ii==0)then
            if(ncnode(i)==nmaxextcells)then
              call diag_print_warning('Reached limit for extended stencil')
              cycle di !No more space
            endif  
            ncnode(i) = ncnode(i) + 1
            cell2extcell(ncnode(i),i) = nck
          endif
        enddo dok
      enddo dj
    enddo di
    
#ifdef DIAG_MODE
    !Check No-repeat cell connectivity
    indictface=1 !indicator for visiting cell face (0 = avoid repeating)
    do i=1,ncellpoly
      do j=1,nxyface(i) !No Repeat of cell faces
        k=kxyface(j,i)
        nck=cell2cell(k,i) !Forward connectivity 
        jcn=llec2llec(k,i) !Backward connectivity
        indictface(k,i)=0
        indictface(jcn,nck)=0
      enddo
    enddo
    do i=1,ncells
      do k=1,ncface(i)
        if(indictface(k,i)==1)then     
          call diag_print_error('Problem calculating no-repeat cell connetivity.',&
            '  Check the input mesh.','  Try reordering the nodes in SMS.')
        endif
      enddo
    enddo
#endif
    
    if(write_geocells)then
      !Geometry file
410   format(i6,3(1pe12.4),i4)
432   format(8x,i2,i8,2x,2(F10.5),20(1pe12.4))
      filename = trim(flowpath) // trim(casename) // '_geo.txt'
      open(67,file=filename)
      write(67,*) ncellsD,nmaxfaces
      do i=1,ncellsD
!        x(i),y(i),fnx(i,:),fny(i,:),ds(i,:),dsx(i,:),dsy(i,:),&
!             rx(i,:),ry(i,:),rnx(i,:),rny(i,:),rpx(i,:),rpy(i,:)
        write(67,'(A)') '     i,       x,          y,        areap,   ncface'
        write(67,410) i,x(i),y(i),areap(i),ncface(i)
        write(67,'(A,A)') '         k, cell2cell,   fnx,      fny,        ds,        dsx,       dsy,      ',&
          '   rx,         ry,        rnx,        rny,        rpx,        rpy,         dsxy,       dn'
        do k=1,ncface(i)
          write(67,432) k,cell2cell(k,i),fnx(k,i),fny(k,i),ds(k,i),dsx(k,i),dsy(k,i),&
            rx(k,i),ry(k,i),rnx(k,i),rny(k,i),rpx(k,i),rpy(k,i),dsxy(k,i),dn(k,i)     
        enddo
      enddo
      close(67)
    endif
    
    return
    end subroutine geo_init_poly

!*************************************************************
    subroutine map_cell_full2active(nstrcells,cells)
! Maps the cell ID's from the full grid to the active grid
! written by Alex Sanchez, USACE-CHL
!*************************************************************
    use geo_def, only: idmap
    use diag_lib
    implicit none
    !Input/Output
    integer,intent(in) :: nstrcells
    integer,intent(inout) :: cells(nstrcells)
    !Internal
    integer :: i,j,k
    character(len=100) :: msg2
    
    do j=1,nstrcells
      i = cells(j)
      k = idmap(i)
      if(k/=0)then          
        cells(j) = k
      else
        write(msg2,*) '  Cell: ',i  
        call diag_print_error('Problem mapping cell from full grid to active',msg2,&
          '  Check all cell references including cell strings and structures.')
      endif
    enddo
    
    return
    end subroutine map_cell_full2active
    
!********************************************************    
    subroutine map_scal_active2full(var,scalout,iwritedry)
! Map scalar to full grid
! written by Alex Sanchez, USACE-CHL
!********************************************************    
    use size_def, only: ncellsD,ncellsfull
    use geo_def, only: idmap
    use flow_def, only: iwet
    use prec_def
    implicit none
    !Input/Output
    integer,intent(in) :: iwritedry
    real(ikind),intent(in) :: var(ncellsD)
    real(4),intent(out) :: scalout(ncellsfull) !Must be single for XMDF libraries
    !Internal Variables
    integer :: i,k
    
    scalout = -999.0 !Initialize
    if(iwritedry==1)then
        do i=1,ncellsfull
          k=idmap(i)  
          if(k/=0)then
            scalout(i) = var(k)
          endif  
        enddo 
    else
        do i=1,ncellsfull
          k=idmap(i)                      
          if(k/=0)then
            if(iwet(k)==1)then
              scalout(i) = var(k)  
            endif  
          endif
        enddo  
    endif
    
    return
    end subroutine map_scal_active2full

!********************************************************    
    subroutine map_vec_active2full(varx,vary,vecout,iwritedry)
!Map scalar to full grid        
!********************************************************    
    use size_def, only: ncellsD,ncellsfull
    use geo_def, only: idmap
    use flow_def, only: iwet
    use prec_def
    implicit none
    !Input/Output
    integer,intent(in) :: iwritedry
    real(ikind),intent(in) :: varx(ncellsD),vary(ncellsD)
    real(4),intent(out) :: vecout(ncellsfull*2)  !Must be single for XMDF Libraries
    !Internal Variables
    integer :: i,i2,k
    
    vecout = -999.0 !Intialize
    if(iwritedry==1)then !Map vector to full grid   
      do i=1,ncellsfull
        k = idmap(i)        
        if(k/=0)then
           i2=i*2 
           vecout(i2-1) = varx(k)
           vecout(i2)   = vary(k) 
        endif
      enddo
    else
      do i=1,ncellsfull
        k = idmap(i)
        if(k/=0)then
           if(iwet(k)==1)then
              i2=i*2 
              vecout(i2-1) = varx(k)
              vecout(i2)   = vary(k)
           endif
        endif
      enddo
    endif 
    
    return
    end subroutine map_vec_active2full    
    
!********************************************************    
    subroutine map_scal_full2active(vtemp,var)
!Map scalar from full 2 active grid    
!********************************************************    
    use size_def, only: ncellsD,ncellsfull
    use geo_def, only: idmap
    use prec_def
    implicit none
    !Input/Output
    real(4),    intent(in) :: vtemp(ncellsfull) !Must be single
    real(ikind),intent(out):: var(ncellsD)
    !Internal Variables
    integer :: i,k
    
    do i=1,ncellsfull
      k=idmap(i)  
      if(k/=0)then
        var(k) = vtemp(i)
      endif
    enddo

    return
    end subroutine map_scal_full2active
    
!********************************************************    
    subroutine map_vec_full2active(vtemp,vecx,vecy)
!Map scalar from full 2 active grid    
!********************************************************    
    use size_def, only: ncellsD,ncellsfull
    use geo_def, only: idmap
    use prec_def
    implicit none
    !Input/Output
    real*4,     intent(in) :: vtemp(ncellsfull*2) !Must be single
    real(ikind),intent(out):: vecx(ncellsD),vecy(ncellsD)
    !Internal Variables
    integer :: i,i2,k
    
    do i=1,ncellsfull
      k=idmap(i)
      if(k/=0)then
        i2=i*2  
        vecx(k)=vtemp(i2-1)
        vecy(k)=vtemp(i2)
      endif
    enddo

    return
    end subroutine map_vec_full2active
        
!**************************************************    
    subroutine proj_horiz_block(kunit,proj)
! Reads a horizontal projection block from the 
! CMS Card file. 
! written by Alex Sanchez, USACE-CHL
!**************************************************
    use geo_def, only: projection,aHorizDatum,&
       aHorizCoordSystem,aHorizUnits
    use diag_lib
    implicit none
    
    !Input/Output
    integer,intent(in) :: kunit
    type(projection),intent(inout) :: proj
    !Internal
    integer :: i,k,ierr,ilen1,ilen2
    integer :: i27=0, i83=0, inad=0
    character(len=37) :: cardname
    character(len=200) :: aline,dtype,azone
    logical :: foundcard,matched
    
d1: do k=1,10
      foundcard = .true.
      read(kunit,*,iostat=ierr) cardname
      if(ierr/=0) exit
      if(cardname(1:1)=='!' .or. cardname(1:1)=='#') cycle
      select case(cardname)  
      case('HORIZONTAL_PROJECTION_END','HORIZ_PROJ_END','END')
        exit d1        
      
      case('COORDINATE_DATUM','HORIZONTAL_DATUM','DATUM')
        backspace(kunit)
        read(kunit,'(A200)') aline
        read(aline,*) cardname, dtype

        matched=.false.
        do i=0,2    !0 = NAD27, 1 = NAD83, 2 = LOCAL
          if(dtype(1:5)==aHorizDatum(i))then
            proj%iHorizDatum = i
            matched=.true.
            exit
          endif
        enddo
        if(.not.matched)then
          inad = index(dtype,'NAD')  
          i27  = index(dtype,'27')
          i83  = index(dtype,'83')
          if (inad .gt. 0) then
            if (i27 > 0) then
              proj%iHorizDatum = 0
              matched=.true.
            elseif (i83 > 0) then
              proj%iHorizDatum = 1
              matched=.true.
            else
              write(*,*)  
              write(*,*) 'ERROR: Invalid Input Horizontal Coordinate Datum: '
              write(*,*) trim(aline)
              write(*,*) '- Will run as LOCAL Datum'
              write(*,*)
              proj%iHorizDatum = 2
            endif
          endif
        endif
       
      case('COORDINATE_SYSTEM','COORD_SYSTEM','SYSTEM','HORIZONTAL_COORDINATE_SYSTEM')
        backspace(kunit)
        read(kunit,*) cardname, aline
        do i=0,22
          if(aline(1:40)==aHorizCoordSystem(i))then
            proj%iHorizCoordSystem = i
            exit
          endif
          !if(i==21)then
          !  write(*,*)  
          !  write(*,*) 'ERROR: Invalid Input Horizontal Coordinate System: '
          !  write(*,*) trim(aline)
          !  read(*,*)
          !  stop
          !endif
        enddo
        
      case('COORDINATE_ZONE','ZONE')
        backspace(kunit)
        read(kunit,'(A200)') aline
        aline=adjustl(aline)
        ilen1=len_trim(aline)
        ilen2=len_trim(cardname)
        if (aline(ilen1:ilen1)==',') then
          aline(ilen1:ilen1)=' '
          ilen1=len_trim(aline)
        endif
        if (ilen1==ilen2) then
          !no zone specified, skip this card
          continue
        else
          read(aline,*) cardname, azone
          read(azone,*,iostat=ierr) proj%iHorizZone
          if(ierr/=0)then
            proj%iHorizZone = 0
          endif
        endif  
        
      case('COORDINATE_UNITS','HORIZONTAL_UNITS','UNITS')
        backspace(kunit)
        read(kunit,*) cardname, aline
        do i=0,4
          if(aline(1:15)==aHorizUnits(i))then
            proj%iHorizUnits = i
            exit
          endif
          !if(i==4)then
          !  write(*,*)  
          !  write(*,*) 'ERROR: Invalid Input Horizontal Coordinate Units: '
          !  write(*,*) trim(aline)
          !  read(*,*)
          !  stop
          !endif
        enddo
        
      case('HORIZONTAL_PROJECTION_BEGIN','HORIZ_PROJ_BEGIN')
        call diag_print_error('Found a horizontal projection block within another ',&
          '  horizontal projection block')
          
      case default
        foundcard = .false.
        
      end select
    enddo d1

    if(proj%iHorizCoordSystem==2 .and. proj%iHorizZone==0)then
      call diag_print_error('Invalid horizontal projection specified',&
        '  Zone must be specified for State Plane Coordinate System')   
    endif
    
    return
    end subroutine proj_horiz_block

!**************************************************    
    subroutine proj_vert_block(kunit,proj)
! Reads a vertical projection block from the 
! CMS Card file. 
! written by Alex Sanchez, USACE-CHL    
!**************************************************
    use geo_def, only: projection,aVertDatum,aVertUnits
    use diag_lib
    implicit none
    !Input/Output
    integer,intent(in) :: kunit
    type(projection),intent(inout) :: proj
    !Internal
    integer :: i,k,ierr
    character(len=37) :: cardname
    character(len=200) :: aline
    logical :: foundcard
    
    
d1: do k=1,10
      foundcard = .true.
      read(kunit,*,iostat=ierr) cardname
      if(ierr/=0) exit
      if(cardname(1:1)=='!' .or. cardname(1:1)=='#') cycle
      select case(cardname)  
      case('DATUM','VERTICAL_DATUM')
        backspace(kunit)
        read(kunit,*) cardname, aline
        ierr = -1
        if (aline(1:7)=='NAVD   ') aline(1:7)='NAVD88 '
        if (aline(1:7)=='NGVD   ') aline(1:7)='NGVD27 '
        do i=0,size(aVertDatum)-1
          if(aline(1:6)==aVertDatum(i))then
            proj%iVertDatum = i
            ierr = 0
            exit
          endif
        enddo
        if(ierr==-1)then
          call diag_print_error('Invalid Input Vertical Coordinate Datum: ',aline)
        endif
        
      case('UNITS','VERTICAL_UNITS')
        backspace(kunit)
        read(kunit,*) cardname, aline
        ierr = -1
        do i=1,size(aVertUnits)
          if(aline(1:6)==aVertUnits(i))then
            proj%iVertUnits = i
            ierr = 0
            exit
          endif          
        enddo
        if(ierr==-1)then
          call diag_print_error('Invalid Input Vertical Coordinate Units: ',aline)
        endif
      
      case('OFFSET','VERTICAL_OFFSET')
        backspace(kunit)
        read(kunit,*) cardname, proj%VertOffset    
      
      case('VERTICAL_PROJECTION_END','VERT_PROJ_END','END')
          exit d1          
      
      case('VERTICAL_PROJECTION_BEGIN','VERT_PROJ_BEGIN')
        call diag_print_error('Found a vertical projection block within another ',&
           '  vertical projection block')
      
      case default
        foundcard = .false.
        
      end select
    enddo d1
    
    return
    end subroutine proj_vert_block

!**********************************************************************
    subroutine get_coord_info
! Written by Mitch Brown, USACE-CHL    
! modified by Alex Sanchez, USACE-CHL
!**********************************************************************
#include "CMS_cpp.h"
#ifdef XMDF_IO
    use geo_def, only: projfl,aHorizDatum,grdfile,proppath,typespath,HProj,VProj
    use geo_lib, only: zone2fip
    use diag_lib
    use xmdf
    
    implicit none
    integer             :: COORD_ID, CID, FILE_ID, ILOC
    integer             :: coordVersion, nn, istart, iend, ierr
    integer             :: HunitsNo, VdatumNo, VunitsNo
    character(len=200)  :: COORDPATH, ProjCS, Hunits, Vdatum, Vunits
    character(len=100)  :: Hdatum, HZone,aHfip
    character(len=1000) :: WKT
      
    call XF_OPEN_FILE(GRDFILE,READONLY,FILE_ID,ierr)
    ILOC=INDEX(TYPESPATH,'/',BACK=.TRUE.)
    PROPPATH = TYPESPATH(1:ILOC-1)
    ILOC=INDEX(PROPPATH,'/',BACK=.TRUE.)
    COORDPATH = PROPPATH(1:ILOC-1)
    call XF_OPEN_GROUP(FILE_ID,TRIM(COORDPATH),COORD_ID,ierr)
    call XF_OPEN_COORDINATE_GROUP(COORD_ID,CID,ierr)
    if(ierr<0)then     !NO COORDINATE GROUP EXISTS, ASSUME LOCAL, m, LOCAL, m
      HPROJ  = '"Local, m"'
      VDATUM = "Local"
      VUNITS = "m"
    else        
      call XF_GET_HORIZ_UNITS(CID,HunitsNo,ierr)  !0 = US Feet, 1 = Intl. Feet, 2 = meters
      select case(HUNITSNO)
        case(0)
          HUNITS="US ft"
          projfl%iHorizUnits=1
        case(1)
          HUNITS="Intl ft"  
          projfl%iHorizUnits=1
        case default !(2)
          HUNITS="m"  
          projfl%iHorizUnits=2
      end select
    
      call XF_GET_VERT_UNITS(CID,VunitsNo,ierr)   !2 = meters
      select case(VunitsNo)
        case(0)
          VUNITS="US ft"
          projfl%iVertUnits=1
        case(1)
          VUNITS="Intl ft"  
          projfl%iVertUnits=1
        case(2)
          VUNITS="m"  
          projfl%iVertUnits=2
      end select
      
      call XF_GET_VERT_DATUM(CID,VdatumNo,ierr)   !0 = Local, 1 = NGVD29, 2 = NAVD88
      select case(VdatumNo)
        case(0)
          VDATUM="Local"
          projfl%iVertDatum=9
        case(1)
          VDATUM="NGVD29"  
          projfl%iVertDatum=8
        case(2)
          VDATUM="NAVD88"  
          projfl%iVertDatum=7
      end select
      
      !call XF_GET_HORIZ_DATUM(CID,coordVersion,ierr)
      coordVersion = 1
      call xf_Get_Coord_Version (CID, coordVersion, ierr)
      call XF_GET_WKT_STRING_SIZE(CID,nn,ierr)
      if(ierr==0)then
        call XF_GET_WKT(CID,WKT(1:nn),ierr)
        ISTART = INDEX(WKT,'PROJCS')
        IEND = INDEX(WKT,'GEOGCS')
        PROJCS = WKT(ISTART+8:IEND-3)
        IEND = INDEX(WKT,'_')
        PROJCS = WKT(ISTART+8:IEND-1)
        ISTART = IEND        
        HZone = WKT(ISTART+1:IEND-3)
        if(PROJCS(1:3)=='UTM')then
          projfl%iHorizCoordSystem = 1
          IEND = INDEX(HZone(6:),'_')
          HZone = HZone(6:IEND+4)
          read(HZone,'(I4)') projfl%iHorizZone
        elseif(PROJCS(1:5)=='NAD83' .or. PROJCS(1:5)=='NAD27')then
          PROJCS = 'State Plane'  
          projfl%iHorizCoordSystem = 2 !State Plane          
          call zone2fip(HZone,aHfip)
          read(aHfip(1:4),'(I4)') projfl%iHorizZone
        endif
        ISTART = INDEX(WKT,'DATUM')
        IEND = INDEX(WKT,'SPHEROID')
        HDatum = WKT(ISTART+7:IEND-3)
        select case(HDatum)
          case('NAD83','D_NORTH_AMERICAN_1983')
            projfl%iHorizDatum=1
            HDatum='NAD83'          
          case('NAD27','D_NORTH_AMERICAN_1927')
            projfl%iHorizDatum=0
            HDatum='NAD27'  
        end select
      else
        call diag_print_error('Error reading WKT_STRING_SIZE')
      endif
    endif    
    call XF_CLOSE_FILE(FILE_ID,ierr)    
    
    !Construct strings
    HPROJ = '"'//trim(HDatum)//', '//trim(PROJCS)//', '//trim(Hzone)//', '//trim(HUNITS)//'"'
    VPROJ = '"'//trim(VDATUM)//', '//trim(VUNITS)//'"'

#endif
    return
    end subroutine get_coord_info
    
!***************************************************************************************
    subroutine grid_cart_write_ascii()
! Writes the CMS-Flow regular and nonuniform Cartesian grid ASCII files
!
! Author: Alex Sanchez, USACE-CHL
!***************************************************************************************    
    use prec_def
    use size_def
    use geo_def
    use out_def, only: outprefix
    implicit none
    
    integer :: i,j,id,kunit,ii
    real(ikind) :: val
    logical :: isrect
    character(len=200)::grdoutfile
    
    !Determine of regular or rectilinear
    isrect = .false.
    do i=1,ncells
      if(abs(dx(1)-dx(i))>1.0e-3 .or. abs(dy(1)-dy(i))>1.0e-3)then
        isrect = .true.
      endif
    enddo
    
222 format(A)
555 format(1pe13.4)
333 format(I5)
444 format(A,A)
342 format(A,1pe13.4,A)
    
    kunit = 784
    grdoutfile = trim(outprefix) // '_grid.cart'   !ASCII Grid
    open(kunit,file=grdoutfile)
    if(isrect)then
      write(kunit,444)  'GRID_TYPE                ','RECTILINEAR'
    else
      write(kunit,444)  'GRID_TYPE                ','REGULAR'
    endif
    write(kunit,'(A,F4.1)') 'GRID_VERSION             ',1.0
    write(kunit,'(A,I5)')   'GRID_MODIFICATION_NUMBER   ',0 
    !Projections
    !call proj_horiz_write(kunit)
    !call proj_vert_write(kunit)
    val = 360.0 - azimuth_fl

    write(kunit,342) 'GRID_ANGLE           ',val," 'deg'"
    write(kunit,342) 'GRID_ORIGIN_X        ',xOrigin," 'm'"
    write(kunit,342) 'GRID_ORIGIN_Y        ',yOrigin," 'm'"
    write(kunit,'(A,I6)')  'GRID_NUMBER_CELLS_X  ',maxcol
    write(kunit,'(A,I6)')  'GRID_NUMBER_CELLS_Y  ',maxrow
    
    if(isrect)then
      !write(kunit,222)    'GRID_COORDSI_BEGIN'
      !do i=1,maxcol
      !  write(kunit,555) dxx(i)
      !enddo
      !write(kunit,222)    'GRID_COORDSI_END'
      !write(kunit,222)    'GRID_COORDSJ_BEGIN'
      !do j=1,maxrow
      !  write(kunit,555) dyy(j)
      !enddo
      !write(kunit,222)    'GRID_COORDSJ_END'
      write(kunit,222)    'GRID_RESOLUTION_X_BEGIN'
      do i=1,maxcol
        if(i==1)then
          val = dxx(1)  
        else
          val = dxx(i)-dxx(i-1)
        endif
        write(kunit,555) val
      enddo
      write(kunit,222)    'GRID_RESOLUTION_X_END'
      write(kunit,222)    'GRID_RESOLUTION_Y_BEGIN'
      do j=1,maxrow 
       if(j==1)then
          val = dyy(1)  
        else
          val = dyy(j)-dyy(j-1)
        endif
        write(kunit,555) val
      enddo
      write(kunit,222)    'GRID_RESOLUTION_Y_END'
    else
      write(kunit,342) 'GRID_RESOLUTION_X    ',dxx(1),' m'
      write(kunit,342) 'GRID_RESOLUTION_Y    ',dyy(1),' m'
    endif
    

    write(kunit,222)      'GRID_DEPTHS_BEGIN'
    do j=1,maxrow
      do i=1,maxcol
        id = i + (j-1)*maxcol
        ii = idmap(id)
        if(ii>0)then
          write(kunit,555) -zb(ii)
        else
          write(kunit,555) -999.0
        endif
      enddo
    enddo        
    write(kunit,222)     'GRID_DEPTHS_END'
    write(kunit,222)      'GRID_ACTIVITY_BEGIN'
    do j=1,maxrow
      do i=1,maxcol
        id = i + (j-1)*maxcol
        ii = idmap(id)
        if(ii>0)then
          write(kunit,333) 1
        else
          write(kunit,333) 0
        endif
      enddo
    enddo
    write(kunit,222)      'GRID_ACTIVITY_END'
    close(kunit)
    
    return
    end subroutine grid_cart_write_ascii
    
!***********************************************************************************    
    subroutine read_grid_cart()
! Reads an ASCII Cartesian (non-telescoping) grid file
! Author: Alex Sanchez, USACE-CHL
!***********************************************************************************    
#include "CMS_cpp.h"
    use size_def
    use geo_def
    use diag_def
    use diag_lib
    use comvarbl
    use fric_def, only: constbotfric,cbotfric
    use bnd_def      
    implicit none
    integer :: ierr,icart,kunit
    integer :: i,ii,j,k
    real(ikind) :: gridversion
    character(len=37) :: cardname,cdum
    logical :: foundcard,foundfile
    
    nmaxfaces = 4  !Maximum # of faces in all directions
    ndmaxfaces = 2 !Maximum # of faces in each direction

    icart = 0
    
!---- Read file ---------------------------------------------
    kunit = 784
    inquire(file=grdfile,exist=foundfile)
    if(.not.foundfile)then
      call diag_print_error('Grid file: ',grdfile,'  not found.')
    endif
    open(kunit,file=grdfile,status='unknown')
    do
      read(kunit,*,iostat=ierr) cardname
      if(ierr/=0) exit
      foundcard = .true.
      select case(cardname)
      case('GRID_TYPE')
        backspace(kunit)
        read(kunit,*) cardname,cdum
        select case(cdum)
        case('RECTILINEAR')
          icart = 2
        case('REGULAR')
          icart = 1  
        end select
        
      case('GRID_MODIFICATION_NUMBER')
        backspace(kunit)
        read(kunit,*,iostat=ierr) cardname,gridversion
        
      case('GRID_ANGLE')
        call card_scalar(kunit,'deg','deg',azimuth_fl,ierr)
        if (SMS_ver .le. 13.0) then
          azimuth_fl = 360.0 - azimuth_fl  !Alex
        endif
              
      case('GRID_ORIGIN_X')
        call card_scalar(kunit,'m','m',xorigin,ierr)
          
      case('GRID_ORIGIN_Y')
        call card_scalar(kunit,'m','m',yorigin,ierr)
          
      case('GRID_NUMBER_CELLS_X')
        backspace(kunit)
        read(kunit,*,iostat=ierr) cardname,maxcol
      
      case('GRID_NUMBER_CELLS_Y')
        backspace(kunit)
        read(kunit,*,iostat=ierr) cardname,maxrow
        
      case('GRID_RESOLUTION_X_BEGIN')
        allocate(dxx(maxcol))
        do i=1,maxcol
          read(kunit,*,iostat=ierr) dxx(i)
          if(ierr/=0)then
            call diag_print_error('Problem reading grid resolution in x-direction')
            exit
          endif
          if(i>1) dxx(i) = dxx(i) + dxx(i-1)
        enddo
        read(kunit,*,iostat=ierr) !Skip end card
      
      case('GRID_RESOLUTION_Y_BEGIN')
        allocate(dyy(maxrow))
        do j=1,maxrow
          read(kunit,*,iostat=ierr) dyy(j)
          if(ierr/=0)then
            call diag_print_error('Problem reading grid resolution in y-direction')
            exit
          endif
          if(j>1) dyy(j) = dyy(j) + dyy(j-1)
        enddo
        read(kunit,*,iostat=ierr) !Skip end card
        
      case('GRID_COORDSI_BEGIN','GRID_COORDS_I_BEGIN') !Reads coordinates in i-direction as written by SMS in the *_grid.h5 file
        allocate(dxx(maxcol))
        read(kunit,*,iostat=ierr) (dxx(i),i=1,maxcol)
        if(ierr/=0)then
          call diag_print_error('Problem reading grid resolution in x-direction')
          exit
        endif
        read(kunit,*,iostat=ierr) !Skip end card
      
      case('GRID_COORDSJ_BEGIN','GRID_COORDS_J_BEGIN') !Reads coordinates in j-direction as written by SMS in the *_grid.h5 file
        allocate(dyy(maxrow))
        read(kunit,*,iostat=ierr) (dyy(j),j=1,maxrow)
        if(ierr/=0)then
          call diag_print_error('Problem reading grid resolution in y-direction')
          exit
        endif        
        read(kunit,*,iostat=ierr) !Skip end card  
        
      case('GRID_DEPTHS_BEGIN')
        ncellsfull = maxcol*maxrow
        allocate(z(ncellsfull))
        read(kunit,*,iostat=ierr) (z(ii),ii=1,ncellsfull)
        if(ierr/=0)then
          call diag_print_error('Problem reading water depths')
          exit
        endif
        read(kunit,*,iostat=ierr) !Skip end card
      
      case('GRID_ACTIVITY_BEGIN','GRID_ACTIVE_BEGIN')
        ncellsfull = maxcol*maxrow
        allocate(cell_type(ncellsfull))
        read(kunit,*,iostat=ierr) (cell_type(ii),ii=1,ncellsfull)
        if(ierr/=0)then
          call diag_print_error('Problem reading cell activity')
          exit
        endif
        read(kunit,*,iostat=ierr) !Skip end card
        
      case default
        foundcard = .false.
        
      end select
    enddo
    
    ncellsfull = maxcol*maxrow
    if(.not.allocated(cell_type))then
        allocate(cell_type(ncellsfull))
        do ii=1,ncellsfull
          if(abs(z(ii)+999.0)<1.0)then
            cell_type = 0
          else
            cell_type = 1
          endif
        enddo
    endif
    
    call geo_cart
    
    return
    end subroutine read_grid_cart
    
!*******************************************************************    
    subroutine geo_cart
!*******************************************************************
    use size_def
    use geo_def
    use prec_def
    use diag_def
    use diag_lib
    implicit none
    integer :: i,j,k,id,icell,icnt,inorth,ieast,isouth,iwest
    
    nmaxfaces = 4  !Maximum # of faces in all directions
    ndmaxfaces = 2 !Maximum # of faces in each direction
    
    !Allocation
    allocate(idmap(ncellsfull),mapid(ncellsfull))      
 
!---- Sum active cells --------------------------------------------
    !get ncells - count active cells which are "1" in cell_type
    !which hold attributes "XMDF: CellTypes"
    ncells = 0
    do j=1,maxrow
      do i=1,maxcol
        id = i + (j-1)*maxcol
        if(cell_type(id)>=1) ncells = ncells+1
      enddo
    enddo

!--- Calculate number of dummy cells -----------------------------------
    !see how many dummy row of cells on boundary cells
    !this is done by searching throught the "full" grid
    !and for each active cell, looking to all sides of it 
    ncellsD = ncells
    do j=1,maxrow
      do i=1,maxcol
        id = i + (j-1)*maxcol
        if(cell_type(id)>=1)then
          !if on edge of full grid, then need dummy cell
          !two dummy cells are assigned for active corner cells 
          if(j==1)then
            ncellsD = ncellsD+1
            inorth = i + j*maxcol
            if(cell_type(inorth)==0)then
                 ncellsD = ncellsD + 1
            endif
          endif
          if(j==maxrow)then
            ncellsD = ncellsD + 1
            isouth = i + (j-2)*maxcol
            if(cell_type(isouth)==0)then
            ncellsD = ncellsD + 1
            endif
          endif
          if(i==1)then
            ncellsD = ncellsD+1
            ieast = i+1 + (j-1)*maxcol
            if(cell_type(ieast)==0)then
              ncellsD = ncellsD + 1
            endif
          endif
          if(i==maxcol) then
            ncellsD = ncellsD + 1
            iwest = i-1 + (j-1)*maxcol
            if(cell_type(iwest)==0)then
              ncellsD = ncellsD + 1
            endif
          endif
          !for interior active cells - look to both sides
          !and see if neighbor cell is active - if not add dummy
          if(j>1.and.j<maxrow)then
            isouth = i + (j-2)*maxcol
            if(cell_type(isouth)==0)then
              ncellsD = ncellsD + 1
            endif
            inorth = i + j*maxcol
            if(cell_type(inorth)==0)then
              ncellsD = ncellsD + 1
            endif
          endif
          if(i>1.and.i<maxcol) then
            iwest = i-1 + (j-1)*maxcol
            if(cell_type(iwest)==0)then
              ncellsD = ncellsD + 1
            endif
            ieast = i+1 + (j-1)*maxcol
            if(cell_type(ieast)==0)then
              ncellsD = ncellsD + 1
            endif
          endif
        endif
      enddo
    enddo

!---- Map from full grid to active grid --------------------------    
    allocate(irow(ncellsD),icol(ncellsD)) 
    irow=0; icol=0
    allocate(zb(ncellsD))
    zb=0.0
    allocate(zbk(nmaxfaces,ncellsD))
    zbk=0.0
    allocate(x(ncellsD),y(ncellsD),dx(ncellsD),dy(ncellsD))
    x=0.0; y=0.0; dx=0.0; dy=0.0        
    !assign values to arrays - mapping from "full" grid to "active" grid
    ncells = 0
    do j=1,maxrow
      do i=1,maxcol
        id = i + (j-1)*maxcol
        idmap(id) = 0        
        if(cell_type(id)>=1)then
          ncells = ncells+1
          idmap(id) = ncells
          mapid(ncells) = id
          zb(ncells) = -Z(id) !Note sign change
          icol(ncells) = i
          irow(ncells) = j
          if(i==1)then
            dx(ncells) = dxx(i)
          else
            dx(ncells) = (dxx(i) - dxx(i - 1)) 
          endif
          if(j==1)then
            dy(ncells) = dyy(j)
          else
            dy(ncells) = (dyy(j) - dyy(j - 1))
          endif
          x(ncells) = dxx(i) - 0.5*dx(ncells)
          y(ncells) = dyy(j) - 0.5*dy(ncells)       
        endif
      enddo
    enddo

!--- Calculate cell connectivity -----------------------
    allocate(cell2cell(4,ncellsD))
    cell2cell = 0
    !determine the location array for the active grid
    do i=1,maxcol
      do j=1,maxrow
        id = i + (j-1)*maxcol
        icell = idmap(id) 
        if(icell>0)then
          isouth = 0
          if(j>1)  isouth = i + (j-2)*maxcol
          inorth = 0
          if(j<maxrow) inorth = i + j*maxcol
          iwest = 0
          if(i>1) iwest = i-1 + (j-1)*maxcol
          ieast = 0
          if(i<maxcol) ieast = i+1 + (j-1)*maxcol
          if(isouth>0)then
            if(cell_type(isouth)>=1) cell2cell(3,icell) = idmap(isouth)
          else
            cell2cell(3,icell) = 0
          endif
          if(inorth>0)then
            if(cell_type(inorth)>=1) cell2cell(1,icell) = idmap(inorth)
          else
            cell2cell(1,icell) = 0
          endif
          if(iwest>0)then
            if(cell_type(iwest)>=1) cell2cell(4,icell) = idmap(iwest)
          else
            cell2cell(4,icell) = 0
          endif
          if(ieast>0)then
            if(cell_type(ieast)>=1) cell2cell(2,icell) = idmap(ieast)
          else
            cell2cell(2,icell) = 0
          endif
        endif
      enddo
    enddo

    !deallocate(cell_type,Z,DXX,DYY)
    deallocate(cell_type,Z)
    
!---- Copy values to dummy cells --------------
    !put dummy row of cells on boundary cells
    ncellsD = ncells
    icnt = 0
    do i=1,ncells    
      if(cell2cell(1,i)<=0)then
        ncellsD = ncellsD + 1
        cell2cell(1,i) = ncellsD
        cell2cell(3,ncellsD) = i
        dx(ncellsD) = dx(i)
        dy(ncellsD) = dy(i)
        x(ncellsD) = x(i)
        y(ncellsD) = y(i) + dy(ncellsD)
        zb(ncellsD) = zb(i)
        icnt = icnt + 1
      endif
      if(cell2cell(2,i)<=0)then
        ncellsD = ncellsD + 1
        cell2cell(2,i) = ncellsD
        cell2cell(4,ncellsD) = i
        dx(ncellsD) = dx(i)
        dy(ncellsD) = dy(i)
        x(ncellsD) = x(i) + dx(ncellsD)
        y(ncellsD) = y(i)
        zb(ncellsD) = zb(i)
        icnt = icnt + 1
      endif
      if(cell2cell(3,i)<=0)then
        ncellsD = ncellsD + 1
        cell2cell(3,i) = ncellsD
        cell2cell(1,ncellsD) = i
        dx(ncellsD) = dx(i)
        dy(ncellsD) = dy(i)
        x(ncellsD) = x(i)
        y(ncellsD) = y(i) - dy(ncellsD)
        zb(ncellsD) = zb(i)
        icnt = icnt + 1
      endif
      if(cell2cell(4,i)<=0)then
        ncellsD = ncellsD + 1
        cell2cell(4,i) = ncellsD
        cell2cell(2,ncellsD) = i
        dx(ncellsD) = dx(i)
        dy(ncellsD) = dy(i)
        x(ncellsD) = x(i) - dx(ncellsD)
        y(ncellsD) = y(i)
        zb(ncellsD) = zb(i)
        icnt = icnt + 1
      endif
    enddo  
           
!--- fix number of faces to 4 --------------------     
    allocate(ncface(ncellsD))
    ncface = 4

!--- Connectivity directions -------------------------------  
    allocate(idirface(4,ncellsD))
    idirface = 0
    do i=1,ncells
      do k=1,4    
        idirface(k,i) = k  !Note: directions are always in the same order *************************        
      enddo
    enddo    
    
!--- Check for single cell lakes --------------------
    do i=1,ncells
      if(cell2cell(1,i)>ncells .and. & 
         cell2cell(2,i)>ncells .and. &
         cell2cell(3,i)>ncells .and. &
         cell2cell(4,i)>ncells)then
8615 format('i=' ,I6,' j=',I6)          
        write(msg2,8615) icol(i),irow(i)  
        call diag_print_error('Invalid ocean cell at:',msg2)
      endif      
    enddo
    
    return
    end subroutine geo_cart