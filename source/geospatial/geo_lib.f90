!=============================================================================
module geo_lib
! Geospatial library
!
! Contains the following:
!   Polygonal variables
!     polyvar - Calculates general polygonal variables
!   Element
!     intriangle - Determines if a point is within a triangle
!     triangle_area - Calculates the area of a triangle
!     inquad - Determines if a point is within a quad
!     quad_area - Calculates the area of a quad
!     quad_split_area - Splits the area of a quad into four triangles 
!                    based on a given point
!     lin_line_intercept - Determines the intercept between two lines
!
!   Grid I/O
!     read_grid14 - Reads an ADCIRC unstructured triangular mesh
!     write_grid14 - Writes an ADCIRC unstructured triangular mesh
!   Grid extraction
!     trisubrect - Extracts a submesh from an unstructured triangular mesh
!                  overlapping a rectangular domain specified by 4 points
!   Horizontal and Vertical Projections
!     proj_default - Sets the default values to a projection variables
!     proj_undef - Sets a projection variable to undefined values
!     proj_horiz_conv - Horizontal projection conversion 
!     reproject27 - Horizontal projection conversion in NAD27
!     reproject83 - Horizontal projection conversion in NAD83
!     zone2fip - Converts a zone name to a FIP
!     zones_create - Creates a database of zones with their
!                   name, FIP, and UTM numbers
!
! written by Alex Sanchez, USACE-CHL
!=============================================================================
#include "CMS_cpp.h"
    implicit none

contains

!****************************************************************
    subroutine polyvar(nmf,nf,xn,yn,area,xc,yc,cnx,cny,&
                       s,sx,sy,rx,ry,rnx,rny,rpx,rpy)
! Calculates polygon information including:
! area - polygon area
! xc,yc - centroid coordinates
! fnx,fny - outward unit vectors
! rnx,rny - normal distance vector from O to f
! rpx,rpy - parallel distance vector from P to O
!
! where P is the cell centroid, f is the cell face center, 
! and O is a ghost point connecting P to f through rn* and rp*
!
! written by Alex Sanchez
!****************************************************************
    use prec_def
    implicit none    
    integer,    intent(in) :: nmf,nf
    real(ikind),intent(in),dimension(nmf) :: xn,yn
    real(ikind),intent(out):: xc,yc
    real(ikind),intent(out),dimension(nmf) :: cnx,cny,s,sx,sy
    real(ikind),intent(out),dimension(nmf) :: rx,rnx,rny,ry,rpx,rpy
    !Internal variables
    integer :: j,k,jj,kk,ierr
    real(ikind):: val,x0,y0,rm,rn,area,xr,yr,px,py
    real(ikind),dimension(nmf) :: xf,yf,xn0,yn0
    
    !Initialize variables
    area = 0.0_ikind
    xc = 0.0_ikind
    yc = 0.0_ikind
    
    !Using a different reference frame avoids precision problems
    x0 = minval(xn(1:nf))
    y0 = minval(yn(1:nf))
    xn0 = xn - x0
    yn0 = yn - y0
    do k=1,nf
      if(k<nf)then
        j=k+1
      else
        j=1
      endif
      val = xn0(k)*yn0(j)-xn0(j)*yn0(k)
      area = area + val
      xf(k) = 0.5_ikind*(xn0(k)+xn0(j))
      yf(k) = 0.5_ikind*(yn0(k)+yn0(j))
      xc = xc + xf(k)*val
      yc = yc + yf(k)*val
    enddo
    area = area/2.0_ikind
    xc = xc/area/3.0_ikind + x0
    yc = yc/area/3.0_ikind + y0
    xf(:) = xf(:) + x0
    yf(:) = yf(:) + y0
    area = abs(area) 

    do k=1,nf
      if(k<nf)then
        j=k+1
      else
        j=1
      endif

      !Normal unit vectors
      sx(k) = xn(j)-xn(k)
      sy(k) = yn(j)-yn(k)
      s(k) = sqrt(sx(k)*sx(k)+sy(k)*sy(k))
      cnx(k) = -sy(k)/s(k)
      cny(k) = sx(k)/s(k)
      
      !Distance vector from P to f
      rx(k) = xf(k)-xc
      ry(k) = yf(k)-yc
      rm = sqrt(rx(k)*rx(k)+ry(k)*ry(k))
      
      !Determine to flip vector
      val=rx(k)*cnx(k)+ry(k)*cny(k)
      if(val<1.0e-7)then
        cnx(k) = -cnx(k)
        cny(k) = -cny(k)
      endif
      
      !Side length vectors
      sx(k) = s(k)*cnx(k)
      sy(k) = s(k)*cny(k)
            
      !Normal and Perpendicular components
      rn = cnx(k)*rx(k)+cny(k)*ry(k)
      rnx(k) = rn*cnx(k)
      rny(k) = rn*cny(k)
      rpx(k) = rx(k)-rnx(k)
      rpy(k) = ry(k)-rny(k)
      
      !Check if recontruction point is within polygon
      !If reconstruction point is not within polygon
      !move to boundary along normal to face.
      !There probably is a fancier way of doing this but this works.
      xr = xc+rpx(k); yr = yc+rpy(k)
      if(.not.inpoly(xr,yr,nmf,nf,xn,yn))then
        do kk=1,nf
          if(kk<nf)then
            jj=kk+1
          else
            jj=1
          endif
          if(kk==k .and. jj==j) cycle !Skip if same face
          call line_line_intercept(xf(k),yf(k),xr,yr,&
                 xn(jj),yn(jj),xn(kk),yn(kk),px,py,ierr)
          if(ierr==0)then !The face normal line and face intercept
            !Reset the distance vector
            rpx(k) = px - xc !no longer perpendicular
            rpy(k) = py - yc !no longer perpendicular
            rnx(k) = xf(k) - px
            rny(k) = yf(k) - py
            exit
          endif
        enddo  
      endif
    enddo
    
    return
    end subroutine polyvar
    
!********************************************************
    function intriangle(xi,yi,xt,yt,at) result(in)
! Determines whether the point (xi,yi) is within the
! triangle (xt(3),yt(3))
! The triangle points can be in any order.
! written by Alex Sanchez, USACE-CHL
!********************************************************    
    use prec_def
    implicit none
    real(ikind),intent(in) :: xi,yi,xt(3),yt(3)
    real(ikind),intent(in),optional :: at
    real(ikind) :: x2(3),y2(3),abc,pab,pbc,pac
    real(ikind) :: suma,err
    logical :: in
    
    !Area of ABC triangle
    if(present(at))then
      abc = at
    else  
      abc = triangle_area(xt,yt)
    endif
    
    if(abc<=1.0e-15)then
      in = .false.
      return
    endif
    
    !Initialize
    x2 = xt; y2 = yt
    
    !Area of PBC triangle
    x2(1) = xi; y2(1) = yi
    pbc = triangle_area(x2,y2)
    
    !Area of PAC triangle        
    x2(2) = xt(1); y2(2) = yt(1)
    pac = triangle_area(x2,y2)
    
    !Area of PAB triangle        
    x2(3) = xt(2);  y2(3) = yt(2)
    pab = triangle_area(x2,y2)
    
    !Total area
    suma = pbc+pac+pab
    
    !Normalized error
    err = abs(suma-abc)/abc  
    if(err<=0.0001)then
      in = .true.
    else
      in = .false.
    endif
    
    return
    end function intriangle
    
!*************************************************
    function triangle_area(xt,yt) result(at)
! Calculates the area of a triangle (x(3),y(3))
! The triangle points can be in any order.
! written by Alex Sanchez, USACE-CHL
!*************************************************      
    use prec_def
    implicit none
    real(ikind),intent(in) :: xt(3),yt(3)
    real(ikind) :: xn(3),yn(3),at
    
    !Subtract min to reduce precision error
    xn = xt-minval(xt)
    yn = yt-minval(yt)    
               
    at = 0.5*abs(xn(1)*(yn(2)-yn(3)) &
           +xn(2)*(yn(3)-yn(1)) &
           +xn(3)*(yn(1)-yn(2)))
       
    return
    end function triangle_area
    
!****************************************************
    function inquad(xi,yi,xq,yq,aq) result(isin)
! Determines if the point (xi,yi) is within
! the quadrilateral (xq(4),yq(4)).
! The order of the points does NOT matter.
! Splits a quad into 4 triangles using the 
! point (xi,yi). If the sum of triangle areas
! is equal to the quad area than the point is
! within the quad.
! The quad is given by 
!  A--B
!  :  :
!  D--C
!
! written by Alex Sanchez, USACE-CHL
!****************************************************    
    use prec_def
    implicit none
    !Input
    real(ikind),intent(in) :: xi,yi         !Point 
    real(ikind),intent(in) :: xq(4),yq(4) !Quadrilateral
    real(ikind),intent(in),optional :: aq   !Area of the quad  
    !Output
    logical :: isin   !true if point is within quadrilateral
    !Internal
    real(ikind) :: aquad
    real(ikind) :: xt(3),yt(3) !Temporary triangle coordinates
    real(ikind) :: at(4)       !Triangle areas
    
    !Area of triangle PAB
    xt = (/xi,xq(1),xq(2)/)
    yt = (/yi,yq(1),yq(2)/)
    at(1) = triangle_area(xt,yt)
    
    !Area of triangle PCB (replace A with C)
    xt(2) = xq(3); yt(2) = yq(3)
    at(2) = triangle_area(xt,yt)
    
    !Area of triangle PCD (replace B with D)
    xt(3) = xq(4); yt(3) = yq(4)
    at(3) = triangle_area(xt,yt)
    
    !Area of triangle PAD (replace C with A)
    xt(2) = xq(1);  yt(2) = yq(1)
    at(4) = triangle_area(xt,yt)
    
    !Area of quad. Optional argument avoids having to recompute in some cases
    if(present(aq))then
      aquad = aq
    else
      aquad = quad_area(xq,yq)
    endif
    
    if(abs(sum(at)-aquad)<0.001*aquad)then
      isin = .true.
    else
      isin = .false.
    endif
    
    return
    end function inquad
    
!**********************************************
    function quad_area(xq,yq) result(aq)
! Calculates the area of a quadrilateral    
! The order of the points does NOT matter.
! written by Alex Sanchez, USACE-CHL
!**********************************************      
    use prec_def
    implicit none
    real(ikind),intent(in) :: xq(4),yq(4)
    real(ikind) :: xn(4),yn(4),aq
    
    !Subtract min to reduce precision error
    xn = xq - minval(xq)
    yn = yq - minval(yq)
    
    aq = xn(1)*(yn(2)-yn(4)) &
       + xn(2)*(yn(3)-yn(1)) &
       + xn(3)*(yn(4)-yn(2)) &
       + xn(4)*(yn(1)-yn(3))
    aq = 0.5*abs(aq)
    
    return
    end function quad_area
    
!***************************************************  
    function poly_area(ns,np,xp,yp) result(ap)
! Calculates the area of a polygon (xp(np),yp(np))
! The triangle points can be in any order.
! written by Alex Sanchez, USACE-CHL
!***************************************************
    use prec_def
    implicit none
    !Input
    integer,    intent(in) :: ns,np
    real(ikind),intent(in) :: xp(ns),yp(ns) !Note: May have a size larger then np
    !Output
    real(ikind) :: ap
    !Internal variables
    integer :: j,k
    real(ikind) :: xp0(np),yp0(np)
    
    !Subtract min value to reduce precision errors
    xp0 = xp - minval(xp(1:np)) 
    yp0 = yp - minval(yp(1:np))
    
    ap = 0.0 !Initialize
    do k=1,np
    if(k<np)then
        j=k+1
      else
        j=1
      endif
      ap = ap + xp0(k)*yp0(j)-xp0(j)*yp0(k)
    enddo
    ap = 0.5*abs(ap)

    return    
    end function poly_area

!*********************************************************
    function inpolyold(xi,yi,ns,np,xp,yp,ap) result(isin)
! Determines if the point (xi,yi) is within
! the quadrilateral (xp(np),yp(np)).
! The order of the points does NOT matter.
! Splits a quad into 4 triangles using the
! point (xi,yi). If the sum of triangle areas
! is equal to the quad area than the point is
! within the polygon.
!
! written by Alex Sanchez, USACE-CHL
!*********************************************************
    use prec_def
    implicit none
    !Input
    real(ikind),intent(in) :: xi,yi
    integer,    intent(in) :: ns,np
    real(ikind),intent(in) :: xp(ns),yp(ns)
    real(ikind),intent(in),optional :: ap
    !Output
    logical :: isin
    !Internal variables
    integer :: j,k
    real(ikind) :: xt(3),yt(3),sumat,apoly,err
    
    !Area of polygon. Optional argument avoids having to recompute in some cases
    if(present(ap))then
      apoly = ap
    else
      apoly = poly_area(ns,np,xp,yp)
    endif
    
    !Compute area of triangles
    sumat = 0.0
    do k=1,np
      if(k<np)then
        j=k+1
      else
        j=1
      endif
      xt = (/xi,xp(k),xp(j)/)
      yt = (/yi,yp(k),yp(j)/)
      sumat = sumat + triangle_area(xt,yt)
      err = (sumat-apoly)/apoly !Normalized error
      if(err>0.005)then
        isin = .false.
        return
      endif
    enddo

    !Determine is point is inside polygon
    err = (sumat-apoly)/apoly !Normalized error
    if(abs(err)<0.001)then
      isin = .true.
    else
      isin = .false.
    endif

    return
    end function inpolyold

!********************************************************
    function inpoly(xi,yi,ns,np,xp,yp) result(isin) 
! Determines if the point (xi,yi) is within
! the polygon (xp(np),yp(np)).
! The size of the arrays xp and yp is ns.
! The order of the points does NOT matter.
! Splits the polygon into triangles using the
! point (xi,yi). If the sum of triangle areas
! is equal to the polygon area than the point is
! within the polygon.
!
! Author: Alex Sanchez, USACE-CHL
! Last Modified: 01/22/14
!********************************************************
    use prec_def
    implicit none 
    !Input
    real(ikind),intent(in) :: xi,yi
    integer,    intent(in) :: ns,np
    real(ikind),intent(in) :: xp(ns),yp(ns)
    !Output
    logical :: isin
    !Internal variables
    integer:: j,k
    real(ikind) :: xv(ns),yv(ns),xa
    
    !Make sorted vertices
    xv = xp; yv = yp
    call poly_sort(ns,np,xv,yv)
    
    isin = .false. 
    do k=1,np 
      if(k<np)then
        j=k+1
      else
        j=1
      endif 
      xa=(xv(j)-xv(k))*(yi-yv(k))/(yv(j)-yv(k))+xv(k)
      if((((yv(k)<=yi) .and. (yi<yv(j))) .or. &
          ((yv(j)<=yi) .and. (yi<yv(k)))) .and. (xi<xa))then
        isin = .not.isin 
      endif
    enddo
    
    return 
    end function inpoly
    
!********************************************************************************    
    subroutine poly_sort(ns,np,xp,yp)
! Sorts the vertices of a polygon
! Author: Alex Sanchez, USACE-CHL
! Last Modified: 01/22/14
!********************************************************************************    
    use prec_def
    use math_lib, only: sortup
    implicit none 
    !Input/Output
    integer,    intent(in) :: ns,np
    real(ikind),intent(inout) :: xp(ns),yp(ns)
    !Internal
    integer :: k,ind(np)
    real(ikind) :: xpm,ypm,b(np)
    
    xpm = sum(xp(1:np))/real(np,kind=ikind)
    ypm = sum(yp(1:np))/real(np,kind=ikind)
    
    do k=1,np
      b(k) = atan2(xp(k)-xpm,yp(k)-ypm)
    enddo
    
    call sortup(np,b,ind)
    xp(1:np) = xp(ind)
    yp(1:np) = yp(ind)
    
    return
    end subroutine poly_sort
    
!******************************************************************************    
    subroutine line_line_intercept(x1,y1,x2,y2,x3,y3,x4,y4,px,py,ierr)
! Determines the intercept between two lines
! Author: Alex Sanchez, USACE-CHL
!******************************************************************************
    use prec_def
    implicit none
    !Input
    real(ikind),intent(in) :: x1,y1,x2,y2,x3,y3,x4,y4
    !Output
    real(ikind),intent(out) :: px,py
    integer, intent(out) :: ierr
    !Internal variables
    real(ikind) :: a,b,d

    d = (x1-x2)*(y3-y4)-(y1-y2)*(x3-x4)
    if(abs(d)<1.0e-10)then !lines are parallel
      px = 0.0; py = 0.0;
      ierr = -1 !Lines do not intercept
      return
    endif
    
    a = x1*y2-y1*x2
    b = x3*y4-y3*x4
    px = (a*(x3-x4)-(x1-x2)*b)/d
    py = (a*(y3-y4)-(y1-y2)*b)/d
    
    if(px<max(x1,x2) .and. px>min(y1,y2) .and. &
       py<max(y1,y2) .and. py>min(y1,y2) .and. &
       px<max(x3,x4) .and. px>min(y3,y4) .and. &
       py<max(y3,y4) .and. py>min(y3,y4))then
      ierr = 0  !Point on both lines
    else
      ierr = 1  !Point not on both lines
    endif

    return
    end subroutine line_line_intercept

!********************************************************************************
    subroutine read_grid14(grd14file,numelems,numnodes,xn,yn,zn,elem2node)
! Reads the parent ADCIRC grid file    
!********************************************************************************
    use diag_def
    use diag_lib
    use prec_def
    implicit none
    !Input/Output
    integer,      intent(out)            :: numelems       !Number of elements
    integer,      intent(out)            :: numnodes       !Number of nodes
    integer,      intent(inout), pointer :: elem2node(:,:) !Element to node connectivity
    real(ikind),  intent(inout), pointer :: xn(:)          !Nodal point global coordinates (output points)
    real(ikind),  intent(inout), pointer :: yn(:)          !Nodal point global coordinates (output points)
    real(ikind),  intent(inout), pointer :: zn(:)          !Nodal point global coordinates (output points)
    character(len=*),intent(in)          :: grd14file      !ADCIRC grid file
    !Internal variables
    integer :: i,k,id,numedges
    logical :: found

    inquire(file=grd14file,exist=found)
    if(.not.found)then
      call diag_print_error('Could not find ADCIRC grid file: ',grd14file)
    endif
    open(unit=14,file=grd14file)
    read(14,*) !Skip first line
    read(14,*) numelems,numnodes
    
    allocate(xn(numnodes),yn(numnodes),zn(numnodes))
    do i=1,numnodes
      read(14,*) id,xn(i),yn(i),zn(i)
    enddo
    allocate(elem2node(3,numelems))
    do i=1,numelems
      read(14,*) id,numedges,(elem2node(k,i),k=1,3)
      if(id/=i)then
        write(msg2,*) '  Element: ',i
        call diag_print_error('Problem reading ADCIRC grid connectivity at ',msg2)
      endif
    enddo
    close(14)
    
    return
    end subroutine read_grid14
    
!********************************************************************
    subroutine write_grid14(grdfile,grdname,ne,nn,xn,yn,zn,e2n)
! Writes an ADCIRC Grid File    
! Author: Alex Sanchez, USACE-CHL
!********************************************************************
    use prec_def
    implicit none
    !Input/Output
    integer,         intent(in) :: ne,nn     !# of elements and nodes
    integer,         intent(in) :: e2n(3,ne) !element to node connectivity
    real(ikind),     intent(in) :: xn(nn),yn(nn),zn(nn) !node coordinates
    character(len=*),intent(in) :: grdfile !Grid file
    character(len=*),intent(in) :: grdname !Grid name
    !Internal variables
    integer :: j,k

222 format(I6,I6)
333 format(I6,3(1x,F20.10))
444 format(I6,1x,I3,3(1x,I6))
    
    open(144,file=grdfile)
    write(144,'(A)') grdname
    write(144,222) ne,nn
    do k=1,nn
      write(144,333) k,xn(k),yn(k),zn(k)
    enddo
    do j=1,ne
      write(144,444) j,3,(e2n(k,j),k=1,3)
    enddo    
    close(144)
    
    return
    end subroutine write_grid14

!********************************************************************    
    subroutine trisubrect(ne,nn,xn,yn,zn,e2n,xmin,xmax,ymin,ymax,nes,nns,xns,yns,zns,e2ns,kns)
! Extracts a submesh from an unstructured triangular overlapping a rectangular domain
! Author: Alex Sanchez, USACE-CHL
!********************************************************************
    use diag_lib
    use prec_def
    implicit none        
    !---- Input/Output -----------------------------------------------------
    !Parent triangular mesh
    integer,       intent(in) :: ne            !Number of elements
    integer,       intent(in) :: nn            !Number of nodes
    integer,       intent(in) :: e2n(3,ne)     !Element to node connectivity
    real(ikind),   intent(in) :: xn(nn),yn(nn),zn(nn) !Node global geometry
    !Rectangular subdomain
    real(ikind),   intent(in) :: xmin,xmax,ymin,ymax !Rectangular domain covering points
    !Child triangular mesh
    integer,     intent(out) :: nes !Number of elements on sub grid
    integer,     intent(out) :: nns !Number of nodes on sub grid
    integer,     intent(out), pointer :: e2ns(:,:) !Element to node mappin on sub grid    
    integer,     intent(out), pointer :: kns(:)  !Node mapping from sub to full grid
    real(ikind), intent(out), pointer :: xns(:),yns(:),zns(:) !Subgrid nodal geometry    
    !---- Internal variables ------------------------------------------------------
    integer :: i,j,k,n
    integer, allocatable :: ies(:)  !Element in (1) or out (0) of rectangular domain
    integer, allocatable :: ins(:)  !Node in (1) or out (0) of rectangular domain
    integer, allocatable :: jes(:)  !Element mapping from sub to full grid    
    integer, allocatable :: jns(:)  !Node mapping from full to sub grid    
    real(ikind) :: xnk,ynk
    
    !--- Subgrid Element sublist -----------------
    !Find elements that overlap rectangular domain
    !and tag elements and nodes
    allocate(ies(ne),ins(nn))  
    ies = 0
    ins = 0    
    do j=1,ne
      do k=1,3
        n = e2n(k,j)  !Element to node on full grid
        xnk = xn(n); ynk = yn(n)
        if(xnk>=xmin .and. xnk<=xmax .and. &
           ynk>=ymin .and. ynk<=ymax)then
          ies(j) = 1    !Add element to list
          ins(e2n(1:3,j)) = 1      !Set nodes as within the rectangular domain
          cycle
        endif   
      enddo
    enddo
    
    !----- Make Element and Node Sublist -------------------
    !Element sublist
    nes = sum(ies) !Number of elements on subgrid
    allocate(jes(nes))
    i=0
    do j=1,ne
      if(ies(j)==1)then
        i = i + 1
        jes(i) = j
      endif
    enddo
    if(i/=nes)then
      call diag_print_error('Problem counting elements on subgrid')
    endif
    
    !Node sublist
    nns = sum(ins)  !Number of nodes on subgrid
    allocate(kns(nns)) !Over-sized, node mapping from sub to full grid
    allocate(jns(nn))  !Node mapping from full to sub grid
    k = 0
    do n=1,nn
      if(ins(n)==1)then
        k = k + 1
        kns(k) = n !Node sub list
        jns(n) = k
      endif
    enddo
    if(k/=nns)then
      call diag_print_error('Problem counting nodes on subgrid')
    endif    
        
    !--- Element to node mapping on subgrid -----------------
    allocate(e2ns(3,nes))
    do i=1,nes
      j = jes(i)
      do k=1,3
        e2ns(k,i) = jns(e2n(k,j))
      enddo
    enddo
    
    !---- Extract subnode locations ------------------------
    allocate(xns(nns),yns(nns),zns(nns))
    do j=1,nns
      k=kns(j)
      xns(j)=xn(k)
      yns(j)=yn(k)
      zns(j)=zn(k)
    enddo    
    
    deallocate(ies,ins,jes,jns)
    
    return
    end subroutine trisubrect
    
!***********************************************************************
    subroutine proj_default(proj)
! Sets the default values for the horizontal and vertical projections
! The default values are local and meters for both
! written by Alex Sanchez, USACE-CHL
!***********************************************************************
    use geo_def, only: projection
    implicit none
    type(projection),intent(inout) :: proj
    
    proj%iHorizDatum = 2         !Horizontal Datum = LOCAL
    proj%iHorizCoordSystem = 22  !Horizontal Coordinate System = LOCAL
    proj%iHorizUnits = 2         !Horizontal Units = METERS
    proj%iVertDatum = 9          !Vertical Datum = LOCAL
    proj%iVertUnits = 2          !Vertical units = METERS
    proj%VertOffset = 0.0        !Vertical Offset from Datum
    
    return
    end subroutine proj_default

!***********************************************************************
    subroutine proj_undef(iundef,proj)
! Sets the default values for the horizontal and vertical projections
! The default values are local and meters for both
! written by Alex Sanchez, USACE-CHL
!***********************************************************************
    use geo_def, only: projection
    implicit none
    integer,intent(in) :: iundef
    type(projection),intent(inout) :: proj
    
    proj%iHorizDatum = iundef         !Horizontal Datum = LOCAL
    proj%iHorizCoordSystem = iundef   !Horizontal Coordinate System = LOCAL
    proj%iHorizUnits = iundef         !Horizontal Units = METERS
    proj%iVertDatum = iundef          !Vertical Datum = LOCAL
    proj%iVertUnits = iundef          !Vertical units = METERS
    proj%VertOffset = real(iundef)    !Vertical Offset from Datum
    
    return
    end subroutine proj_undef

!************************************************************************
    subroutine proj_horiz_conv(projfrom,projto,ncalc,xpts,ypts)
! Horizontal Projection Conversion
!************************************************************************    
#include "CMS_cpp.h"
    use geo_def, only: projection
    use diag_lib
    use prec_def
    implicit none
    integer,     intent(in) :: ncalc
    real(ikind), intent(inout) :: xpts(ncalc),ypts(ncalc)
    type(projection) :: projfrom,projto
      
    if(projfrom%iHorizDatum/=projto%iHorizDatum)then
      call diag_print_error('Cannot convert between different Horizontal Projections')
    endif
    
!   if(projfrom%iHorizDatum==projto%iHorizDatum .and. &                !We know these are the same, otherwise the program would have stopped already
    if(projfrom%iHorizCoordSystem==projto%iHorizCoordSystem .and. &
      projfrom%iHorizUnits==projto%iHorizUnits)then
      !write(*,*) ' Projections are the same'  
      return 
    endif
      
    if(projfrom%iHorizDatum==2 .and. projto%iHorizDatum/=2)then !from local
      call diag_print_warning('Cannot convert from a local horizontal projection',&
        '   No horizontal projection conversion applied')
      return
    elseif(projfrom%iHorizDatum/=2 .and. projto%iHorizDatum==2)then !to local
      call diag_print_warning('Cannot convert to a local horizontal projection',&
        '   No horizontal projection conversion applied')
      return  
    endif
      
#ifdef PROJ_CONV
    if(projfrom%iHorizDatum==0)then !NAD27
      call reproject27(projfrom,projto,ncalc,xpts,ypts)  
    elseif(projfrom%iHorizDatum==1)then !NAD83
      call reproject83(projfrom,projto,ncalc,xpts,ypts)  
    endif    
#else
    call diag_print_error('Cannot Convert horizontal projection without projection libraries')
#endif

    return
    end subroutine proj_horiz_conv

#ifdef PROJ_CONV
!************************************************************************
    subroutine reproject27(projfrom,projto,ncalc,xpts,ypts)
! Inteface to gctp.f NAD27 projection conversion library from USGS
!************************************************************************    
    use geo_def, only: projection
    use prec_def
    implicit none
    !Input/Output
    integer,     intent(in) :: ncalc
    real(ikind), intent(inout) :: xpts(ncalc),ypts(ncalc)
    type(projection) :: projfrom,projto
    !Internal Variables
    integer :: i
    integer :: INSYS,JNZONE,INUNIT,INSPH,IPR,JPR
    integer :: IOSYS,IOZONE,IOUNIT,IOSPH,IFLG
    real(8) :: CRDIN(2),CRDIO(2),TPARIN(15),TPARIO(15)
    
    !Input Coordinates
    INSYS = projfrom%iHorizCoordSystem !Coordinate System
    JNZONE = projfrom%iHorizZone       !Zone
    INUNIT = projfrom%iHorizUnits      !Units
    INSPH = 0                          !Sphere 

    !Output Coordinates
    IOSYS = projto%iHorizCoordSystem !Coordinate System
    IOZONE = projto%iHorizZone       !Zone
    IOUNIT = projto%iHorizUnits      !Units
    IOSPH = 0                        !Sphere 
    
    !Printing
    IPR = 0                            !Print error messages 0-yes, 1-no
    JPR = 0                            !Print parameters 0-yes, 1-no
    
    do i=1,ncalc
      CRDIN(1)=xpts(i)
      CRDIN(2)=ypts(i)
      call GTPZ0(CRDIN,INSYS,JNZONE,TPARIN,INUNIT,INSPH,IPR,JPR,&
                      CRDIO,IOSYS,IOZONE,TPARIO,IOUNIT,IOSPH,IFLG)
      xpts(i)=CRDIO(1)
      ypts(i)=CRDIO(2)
    enddo
    
    return
    end subroutine reproject27
    
!*************************************************************************
    subroutine reproject83(projfrom,projto,ncalc,xpts,ypts)
! Inteface to SPCS83 and UTM2GEO packages for NAD83 projection conversions
!*************************************************************************
    use geo_def, only: projection,aHorizCoordSystem,avg_lat
    use SPCS83_COMBINED
    use UTM2GEO
    use diag_lib
    use prec_def
    implicit none
    !Input/Output
    integer,     intent(in) :: ncalc
    real(ikind), intent(inout) :: xpts(ncalc),ypts(ncalc)
    type(projection) :: projfrom,projto    
    
    !Internal variables
    integer :: grid_zone(2)
    integer :: ICODE(3),I,WGS_84_DATUM
    real(8) :: RLON,RLAT,NORTH,EAST !Must be double
    real(8) :: lambda0
    integer :: ZCODE
    real(ikind), parameter :: rad2deg = 57.295779513082321
    logical,parameter :: EWFLAG = .FALSE.     
    
    if(projfrom%iHorizCoordSystem==0)then !From Geographic      
      if(projto%iHorizCoordSystem==1)then !To UTM
        ZCODE = projto%iHorizZone
        WGS_84_DATUM = 3
        RLON = -177.0d0 + (zcode-1)*6.0d0  !Get the central merdian for this zone
        RLAT = avg_lat !sum(adc_lat(:))/nnodes      !Set the origin latitude to zero  (could use the average over the domain)
        call get_grid_zone(RLON,RLAT,ICODE,lambda0)  !Used to get lambda0
        ICODE(1) = projto%iHorizZone 
        ICODE(2) = 0 ! This should always be zero   
        do i=1,ncalc
          RLON=xpts(i)
          RLAT=ypts(i)
          call ll2utm(RLON,RLAT,EAST,NORTH,ICODE(1:2),lambda0,WGS_84_DATUM)  !Convert to UTM
          xpts(i)=EAST
          ypts(i)=NORTH
        enddo
      elseif(projto%iHorizCoordSystem==2)then !To State Plane
        ZCODE = projto%iHorizZone
        ICODE(1) = ZCODE
        ICODE(2) = 0 ! These should always be zero
        ICODE(3) = 0 ! These should always be zero  
        do i=1,ncalc
          RLON=abs(xpts(i))
          RLAT=ypts(i)   
          CALL GEO2STPL(RLON,RLAT,ICODE,EWFLAG,NORTH,EAST) 
          xpts(i)=EAST
          ypts(i)=NORTH
        enddo
      else
        goto 200
      endif
    elseif(projfrom%iHorizCoordSystem==1)then !from UTM
      if(projto%iHorizCoordSystem==0)then !to Geographic
        call diag_print_warning("Potential Issues with UTM to LL mapping.", "- Use State Plane for grid, if needed.")
        ZCODE = projfrom%iHorizZone
        WGS_84_DATUM = 3
        RLON = -177.0d0 + (zcode-1)*6.0d0  !Get the central merdian for this zone
        RLAT = avg_lat !sum(adc_lat(:))/nnodes      !Set the origin latitude to zero  (could use the average over the domain)
        call get_grid_zone(RLON,RLAT,grid_zone,lambda0)
        ICODE(1) = projto%iHorizZone 
        ICODE(2) = 0 ! This should always be zero   
        do i=1,ncalc
          EAST =xpts(i)
          NORTH=ypts(i)
          call utm2ll(EAST,NORTH,RLON,RLAT,grid_zone,WGS_84_DATUM)
          !xpts(i)=RLON*rad2deg       !Test MEB  01/05/2024   The UTM2LL subroutine already returns in degrees.  This caused bad values.
          !ypts(i)=RLAT*rad2deg       !                       How did this ever work?
          xpts(i)=RLON
          ypts(i)=RLAT
          if(xpts(i)<180.0) xpts(i)=-xpts(i)
        enddo
      else
        goto 200
      endif
    elseif(projfrom%iHorizCoordSystem==2)then !from State Plane
      if(projto%iHorizCoordSystem==0)then !to Geographic
        ZCODE = projfrom%iHorizZone
        ICODE(1) = ZCODE
        ICODE(2) = 0 ! These should always be zero
        ICODE(3) = 0 ! These should always be zero  
        do i=1,ncalc
          EAST =xpts(i)
          NORTH=ypts(i)
          CALL STPL2GEO(NORTH,EAST,ICODE,RLAT,RLON)
          xpts(i)=RLON*rad2deg
          ypts(i)=RLAT*rad2deg
          if(xpts(i)<180.0) xpts(i)=-xpts(i)          
        enddo
      else
        goto 200  
      endif
    endif
    
    return
200 call diag_print_error('Cannot convert from ',&
      aHorizCoordSystem(projfrom%iHorizCoordSystem),' to ',&
      aHorizCoordSystem(projto%iHorizCoordSystem))
    
    end subroutine reproject83
#endif
    
!***********************************************************    
    subroutine zone2fip(azone,afip)
!***********************************************************    
    use geo_def, only: nzones,zones
    implicit none
    !Input/Output
    character(len=100),intent(in) :: azone
    character(len=100),intent(out) :: afip
    !Internal
    integer :: i,nn
    character(len=100) :: HZone
    
    !Create zones
    call zones_create
    
    HZone = azone
    call remove_underscores(HZone)  
    call uppercase(HZone)
    nn=len_trim(HZone)    
    do i=1,nzones           
      if(zones(i)%name==HZone(1:nn))then
        afip=zones(i)%fip
        exit
      endif            
    enddo        
    
    !Destroy zones
    deallocate(zones)
    nzones = 0
    
    return
    end subroutine zone2fip
    
!***********************************************************    
    subroutine zones_create
! Creates a list of UTM zones and their information
!***********************************************************
    use geo_def, only: nzones,zones
    implicit none
    
    nzones = 127    
    allocate(zones(nzones))
    
    zones(1)%name = 'ALASKA ZONE 1'
    zones(1)%fip = '5001'
    zones(1)%utm = '8 & 9'

    zones(2)%name = 'ALASKA ZONE 2'
    zones(2)%fip = '5002'
    zones(2)%utm = '7'

    zones(3)%name = 'ALASKA ZONE 3'
    zones(3)%fip = '5003'
    zones(3)%utm = '6'

    zones(4)%name = 'ALASKA ZONE 4'
    zones(4)%fip = '5004'
    zones(4)%utm = '5 & 6'

    zones(5)%name = 'ALASKA ZONE 5'
    zones(5)%fip = '5005'
    zones(5)%utm = '5'

    zones(6)%name = 'ALASKA ZONE 6'
    zones(6)%fip = '5006'
    zones(6)%utm = '4'

    zones(7)%name = 'ALASKA ZONE 7'
    zones(7)%fip = '5007'
    zones(7)%utm = '3 & 4'

    zones(8)%name = 'ALASKA ZONE 8'
    zones(8)%fip = '5008'
    zones(8)%utm = '3'

    zones(9)%name = 'ALASKA ZONE 9'
    zones(9)%fip = '5009'
    zones(9)%utm = ''

    zones(10)%name = 'ALASKA ZONE 10'
    zones(10)%fip = '5010'
    zones(10)%utm = '1 & 4'

    zones(11)%name = 'ALABAMA EAST ZONE'
    zones(11)%fip = '0101'
    zones(11)%utm = '16'

    zones(12)%name = 'ALABAMA WEST ZONE'
    zones(12)%fip = '0102'
    zones(12)%utm = '16'

    zones(13)%name = 'ARIZONA EAST ZONE'
    zones(13)%fip = '0201'
    zones(13)%utm = '12'

    zones(14)%name = 'ARIZONA CENTRAL ZONE'
    zones(14)%fip = '0202'
    zones(14)%utm = '12'

    zones(15)%name = 'ARIZONA WEST ZONE'
    zones(15)%fip = '0203'
    zones(15)%utm = '12'

    zones(16)%name = 'ARKANSAS NORTH ZONE'
    zones(16)%fip = '0301'
    zones(16)%utm = '15'

    zones(17)%name = 'ARKANSAS SOUTH ZONE'
    zones(17)%fip = '0302'
    zones(17)%utm = '15'

    zones(18)%name = 'CALIFORNIA ZONE I'
    zones(18)%fip = '0401'
    zones(18)%utm = '10'

    zones(19)%name = 'CALIFORNIA ZONE II'
    zones(19)%fip = '0402'
    zones(19)%utm = '10 & 11'

    zones(20)%name = 'CALIFORNIA ZONE III'
    zones(20)%fip = '0403'
    zones(20)%utm = '10 & 11'

    zones(21)%name = 'CALIFORNIA ZONE IV'
    zones(21)%fip = '0404'
    zones(21)%utm = '10 & 11'

    zones(22)%name = 'CALIFORNIA ZONE V'
    zones(22)%fip = '0405'
    zones(22)%utm = '10 & 11'

    zones(23)%name = 'CALIFORNIA ZONE VI'
    zones(23)%fip = '0406'
    zones(23)%utm = '11'

    zones(24)%name = 'CALIFORNIA ZONE VII'
    zones(24)%fip = '0407'
    zones(24)%utm = ''

    zones(25)%name = 'COLORADO NORTH ZONE'
    zones(25)%fip = '0501'
    zones(25)%utm = '12 & 13'

    zones(26)%name = 'COLORADO CENTRAL ZONE'
    zones(26)%fip = '0502'
    zones(26)%utm = '12 & 13'

    zones(27)%name = 'COLORADO SOUTH ZONE'
    zones(27)%fip = '0503'
    zones(27)%utm = '12 & 13'

    zones(28)%name = 'CONNECTICUT'
    zones(28)%fip = '0600'
    zones(28)%utm = '18 & 19'

    zones(29)%name = 'DELAWARE'
    zones(29)%fip = '0700'
    zones(29)%utm = '18'

    zones(30)%name = 'FLORIDA EAST ZONE'
    zones(30)%fip = '0901'
    zones(30)%utm = '17'

    zones(31)%name = 'FLORIDA WEST ZONE'
    zones(31)%fip = '0902'
    zones(31)%utm = '17'

    zones(32)%name = 'FLORIDA NORTH ZONE'
    zones(32)%fip = '0903'
    zones(32)%utm = '16 & 17'

    zones(33)%name = 'GEORGIA EAST ZONE'
    zones(33)%fip = '1001'
    zones(33)%utm = '17'

    zones(34)%name = 'GEORGIA WEST ZONE'
    zones(34)%fip = '1002'
    zones(34)%utm = '16 & 17'

    zones(35)%name = 'HAWAII ZONE 1'
    zones(35)%fip = '5101'
    zones(35)%utm = '5'

    zones(36)%name = 'HAWAII ZONE 2'
    zones(36)%fip = '5102'
    zones(36)%utm = '4'

    zones(37)%name = 'HAWAII ZONE 3'
    zones(37)%fip = '5103'
    zones(37)%utm = '4'

    zones(38)%name = 'HAWAII ZONE 4'
    zones(38)%fip = '5104'
    zones(38)%utm = '4'

    zones(39)%name = 'HAWAII ZONE 5'
    zones(39)%fip = '5105'
    zones(39)%utm = '4'

    zones(40)%name = 'IDAHO EAST ZONE'
    zones(40)%fip = '1101'
    zones(40)%utm = '12'

    zones(41)%name = 'IDAHO CENTRAL ZONE'
    zones(41)%fip = '1102'
    zones(41)%utm = '11 & 12'

    zones(42)%name = 'IDAHO WEST ZONE'
    zones(42)%fip = '1103'
    zones(42)%utm = '11'

    zones(43)%name = 'ILLINOIS EAST ZONE'
    zones(43)%fip = '1201'
    zones(43)%utm = '16'

    zones(44)%name = 'ILLINOIS WEST ZONE'
    zones(44)%fip = '1202'
    zones(44)%utm = '15 & 16'

    zones(45)%name = 'INDIANA EAST ZONE'
    zones(45)%fip = '1301'
    zones(45)%utm = '16'

    zones(46)%name = 'INDIANA WEST ZONE'
    zones(46)%fip = '1302'
    zones(46)%utm = '16'

    zones(47)%name = 'IOWA NORTH ZONE'
    zones(47)%fip = '1401'
    zones(47)%utm = '14 & 15'

    zones(48)%name = 'IOWA SOUTH ZONE'
    zones(48)%fip = '1402'
    zones(48)%utm = '15'

    zones(49)%name = 'KANSAS NORTH ZONE'
    zones(49)%fip = '1501'
    zones(49)%utm = '14 & 15'

    zones(50)%name = 'KANSAS SOUTH ZONE'
    zones(50)%fip = '1502'
    zones(50)%utm = '14 & 15'

    zones(51)%name = 'KENTUCKY NORTH ZONE'
    zones(51)%fip = '1601'
    zones(51)%utm = '16 & 17'

    zones(52)%name = 'KENTUCKY SOUTH ZONE'
    zones(52)%fip = '1602'
    zones(52)%utm = '16 & 17'

    zones(53)%name = 'LOUISIANA NORTH ZONE'
    zones(53)%fip = '1701'
    zones(53)%utm = '15'

    zones(54)%name = 'LOUISIANA SOUTH ZONE'
    zones(54)%fip = '1702'
    zones(54)%utm = '15 & 16'

    zones(55)%name = 'MAINE EAST ZONE'
    zones(55)%fip = '1801'
    zones(55)%utm = '19'

    zones(56)%name = 'MAINE WEST ZONE'
    zones(56)%fip = '1802'
    zones(56)%utm = '19'

    zones(57)%name = 'MARYLAND'
    zones(57)%fip = '1900'
    zones(57)%utm = '17 & 18'

    zones(58)%name = 'MASSACHUSETTS MAINLAND'
    zones(58)%fip = '2001'
    zones(58)%utm = '18 & 19'

    zones(59)%name = 'MASSACHUSETTS ISLAND'
    zones(59)%fip = '2002'
    zones(59)%utm = '19'

    zones(60)%name = 'MICHIGAN NORTH ZONE'
    zones(60)%fip = '2111'
    zones(60)%utm = '16'

    zones(61)%name = 'MICHIGAN CENTRAL ZONE'
    zones(61)%fip = '2112'
    zones(61)%utm = '16 & 17'

    zones(62)%name = 'MICHIGAN SOUTH ZONE'
    zones(62)%fip = '2113'
    zones(62)%utm = '16 & 17'

    zones(63)%name = 'MINNESOTA NORTH ZONE'
    zones(63)%fip = '2201'
    zones(63)%utm = '14 & 15'

    zones(64)%name = 'MINNESOTA CENTRAL ZONE'
    zones(64)%fip = '2202'
    zones(64)%utm = '14 & 15'

    zones(65)%name = 'MINNESOTA SOUTH ZONE'
    zones(65)%fip = '2203'
    zones(65)%utm = '14 & 15'

    zones(66)%name = 'MISSISSIPPI EAST ZONE'
    zones(66)%fip = '2301'
    zones(66)%utm = '16'

    zones(67)%name = 'MISSISSIPPI WEST ZONE'
    zones(67)%fip = '2302'
    zones(67)%utm = '15 & 16'

    zones(68)%name = 'MISSOURI EAST ZONE'
    zones(68)%fip = '2401'
    zones(68)%utm = '15 & 16'

    zones(69)%name = 'MISSOURI CENTRAL ZONE'
    zones(69)%fip = '2402'
    zones(69)%utm = '15'

    zones(70)%name = 'MISSOURI WEST ZONE'
    zones(70)%fip = '2403'
    zones(70)%utm = '15'

    zones(71)%name = 'MONTANA'
    zones(71)%fip = '2500'
    zones(71)%utm = '11 & 12 & 13' !NAD83 IS SINGLE ZONE

    zones(72)%name = 'MONTANA NORTH ZONE' !NAD27
    zones(72)%fip = '2501'
    zones(72)%utm = '11 & 12 & 13'

    zones(73)%name = 'MONTANA CENTRAL ZONE' !NAD27
    zones(73)%fip = '2502'
    zones(73)%utm = '11 & 12 & 13'

    zones(74)%name = 'MONTANA SOUTH ZONE' !NAD27
    zones(74)%fip = '2503'
    zones(74)%utm = '11 & 12 & 13'

    zones(75)%name = 'NEBRASKA'
    zones(75)%fip = '2600'
    zones(75)%utm = '13 & 14 & 15' !NAD83 IS SINGLE ZONE

    zones(76)%name = 'NEBRASKA NORTH ZONE' !NAD27
    zones(76)%fip = '2601'
    zones(76)%utm = '13 & 14 & 15'

    zones(77)%name = 'NEBRASKA SOUTH ZONE' !NAD27
    zones(77)%fip = '2602'
    zones(77)%utm = '13 & 14 & 15'

    zones(78)%name = 'NEVADA EAST ZONE'
    zones(78)%fip = '2701'
    zones(78)%utm = '11'

    zones(79)%name = 'NEVADA CENTRAL ZONE'
    zones(79)%fip = '2702'
    zones(79)%utm = '11'

    zones(80)%name = 'NEVADA WEST ZONE'
    zones(80)%fip = '2703'
    zones(80)%utm = '11'

    zones(81)%name = 'NEW HAMPSHIRE'
    zones(81)%fip = '2800'
    zones(81)%utm = '18 & 19'

    zones(82)%name = 'NEW JERSEY'
    zones(82)%fip = '2900'
    zones(82)%utm = '18'

    zones(83)%name = 'NEW MEXICO EAST ZONE'
    zones(83)%fip = '3001'
    zones(83)%utm = '13'

    zones(84)%name = 'NEW MEXICO CENTRAL ZONE'
    zones(84)%fip = '3002'
    zones(84)%utm = '13'

    zones(85)%name = 'NEW MEXICO WEST ZONE'
    zones(85)%fip = '3003'
    zones(85)%utm = '12 & 13'

    zones(86)%name = 'NEW YORK EAST ZONE'
    zones(86)%fip = '3101'
    zones(86)%utm = '18'

    zones(87)%name = 'NEW YORK CENTRAL ZONE'
    zones(87)%fip = '3102'
    zones(87)%utm = '18'

    zones(88)%name = 'NEW YORK WEST ZONE'
    zones(88)%fip = '3103'
    zones(88)%utm = '17 & 18'

    zones(89)%name = 'NEW YORK LONG ISLAND ZONE'
    zones(89)%fip = '3104'
    zones(89)%utm = '18'

    zones(90)%name = 'NORTH DAKOTA NORTH ZONE'
    zones(90)%fip = '3301'
    zones(90)%utm = '13 & 14'

    zones(91)%name = 'NORTH DAKOTA SOUTH ZONE'
    zones(91)%fip = '3302'
    zones(91)%utm = '13 & 14'

    zones(92)%name = 'OHIO NORTH ZONE'
    zones(92)%fip = '3401'
    zones(92)%utm = '16 & 17'

    zones(93)%name = 'OHIO SOUTH ZONE'
    zones(93)%fip = '3402'
    zones(93)%utm = '16 & 17'

    zones(94)%name = 'OKLAHOMA NORTH ZONE'
    zones(94)%fip = '3501'
    zones(94)%utm = '13 & 14 & 15'

    zones(95)%name = 'OKLAHOMA SOUTH ZONE'
    zones(95)%fip = '3502'
    zones(95)%utm = '14 & 15'

    zones(96)%name = 'OREGON NORTH ZONE'
    zones(96)%fip = '3601'
    zones(96)%utm = '10 & 11'

    zones(97)%name = 'OREGON SOUTH ZONE'
    zones(97)%fip = '3602'
    zones(97)%utm = '10 & 11'

    zones(98)%name = 'PENNSYLVANIA NORTH ZONE'
    zones(98)%fip = '3701'
    zones(98)%utm = '17 & 18'

    zones(99)%name = 'PENNSYLVANIA SOUTH ZONE'
    zones(99)%fip = '3702'
    zones(99)%utm = '17 & 18'

    zones(100)%name = 'RHODE ISLAND'
    zones(100)%fip = '3800'
    zones(100)%utm = '19'

    zones(101)%name = 'SOUTH CAROLINA'
    zones(101)%fip = '3900'
    zones(101)%utm = '17'

    zones(102)%name = 'SOUTH CAROLINA NORTH ZONE (NAD27)'
    zones(102)%fip = '3901'
    zones(102)%utm = '17'

    zones(103)%name = 'SOUTH CAROLINA SOUTH ZONE (NAD27)'
    zones(103)%fip = '3902'
    zones(103)%utm = '17'

    zones(104)%name = 'SOUTH DAKOTA NORTH ZONE'
    zones(104)%fip = '4001'
    zones(104)%utm = '13 & 14'

    zones(105)%name = 'SOUTH DAKOTA SOUTH ZONE'
    zones(105)%fip = '4002'
    zones(105)%utm = '13 & 14'

    zones(106)%name = 'TENNESSEE'
    zones(106)%fip = '4100'
    zones(106)%utm = '16 & 17'

    zones(107)%name = 'TEXAS NORTH ZONE'
    zones(107)%fip = '4201'
    zones(107)%utm = '13 & 14'

    zones(108)%name = 'TEXAS NORTH CENTRAL ZONE'
    zones(108)%fip = '4202'
    zones(108)%utm = '13 & 14 & 15'

    zones(109)%name = 'TEXAS CENTRAL ZONE'
    zones(109)%fip = '4203'
    zones(109)%utm = '13 & 14 & 15'

    zones(110)%name = 'TEXAS SOUTH CENTRAL ZONE'
    zones(110)%fip = '4204'
    zones(110)%utm = '13 & 14 & 15'

    zones(111)%name = 'UTAH NORTH ZONE'
    zones(111)%fip = '4301'
    zones(111)%utm = '12'

    zones(112)%name = 'UTAH CENTRAL ZONE'
    zones(112)%fip = '4302'
    zones(112)%utm = '12'

    zones(113)%name = 'UTAH SOUTH ZONE'
    zones(113)%fip = '4303'
    zones(113)%utm = '12'

    zones(114)%name = 'VERMONT'
    zones(114)%fip = '4400'
    zones(114)%utm = '18 & 19'

    zones(115)%name = 'VIRGINIA NORTH ZONE'
    zones(115)%fip = '4501'
    zones(115)%utm = '17 & 18'

    zones(116)%name = 'VIRGINIA SOUTH ZONE'
    zones(116)%fip = '4502'
    zones(116)%utm = '17 & 18'

    zones(117)%name = 'WASHINGTON NORTH ZONE'
    zones(117)%fip = '4601'
    zones(117)%utm = '10 & 11'

    zones(118)%name = 'WASHINGTON SOUTH ZONE'
    zones(118)%fip = '4602'
    zones(118)%utm = '10 & 11'

    zones(119)%name = 'WEST VIRGINIA NORTH ZONE'
    zones(119)%fip = '4701'
    zones(119)%utm = '17 & 18'

    zones(120)%name = 'WEST VIRGINIA SOUTH ZONE'
    zones(120)%fip = '4702'
    zones(120)%utm = '17'

    zones(121)%name = 'WISCONSIN NORTH ZONE'
    zones(121)%fip = '4801'
    zones(121)%utm = '15 & 16'

    zones(122)%name = 'WISCONSIN CENTRAL ZONE'
    zones(122)%fip = '4802'
    zones(122)%utm = '15 & 16'

    zones(123)%name = 'WISCONSIN SOUTH ZONE'
    zones(123)%fip = '4803'
    zones(123)%utm = '15 & 16'

    zones(124)%name = 'WYOMING EAST ZONE'
    zones(124)%fip = '4901'
    zones(124)%utm = '13'

    zones(125)%name = 'WYOMING EAST CENTRAL ZONE'
    zones(125)%fip = '4902'
    zones(125)%utm = '13'

    zones(126)%name = 'WYOMING WEST CENTRAL ZONE'
    zones(126)%fip = '4903'
    zones(126)%utm = '12'

    zones(127)%name = 'WYOMING WEST ZONE'
    zones(127)%fip = '4904'
    zones(127)%utm = '12' 

    return
    end subroutine zones_create
    
    
end module geo_lib
