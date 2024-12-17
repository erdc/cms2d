!==================================================================================
module interp_lib
! Interpolation Library
!
! Contains the following routines:
! ~ 1D Interpolation ~
!     interp1dlin - Simple 1D linear interpolation function
!     interp1dquad - Simple 1D quadratic interpolation Function
!  
! ~ 2D Interpolation ~
!     Element interpolation
!       interp_coef_tri - Linear interpolation coefficients for a triangle
!       interp_coef_quad - Interpolation for a quadrilateral
!     Grid-to-grid interpolations
!       interp_coef_tel2cart - Interpolation coefficients and index from a
!                           telescoping Cartesian grid (parent grid) to 
!                           nonuniform Cartesian grid (child grid)
!       interp_coef_cart2tel - Interpolation coefficients and index from a 
!                           nonuniform Cartesian grid (parent grid) to a
!                           telescoping Cartesian grid (child grid). 
!       interp_coef_tel2tel  - Interpolation coefficients and index from a
!                           telescoping Cartesian grid to another telescoping grid
!       interp_coef_tel2pts  - Interpolation coefficients and index from a
!                           telescoping Cartesian grid to scattered points
!       interp_coef_tri2tel  - Interpolation coefficients and index from an unstructured
!                            node-based triangular mesh to a telescoping grid
!       interp_coef_tri2pts  - Interpolation coefficients from an unstructured 
!                           node-based triangular mesh to scattered points
!       interp_coef_poly2pts  - Interpolation coefficients from an unstructured 
!                           centroid-based polygonal mesh to scattered points
!       interp_coef_curv2tel - Interpolation coefficients and indes from a 
!                           curvilinear grid to a telescoping grid
!       interp_coef_curv2pts - Interpolation coefficients and index from a 
!                           curvilinear grid to a scattered points
!       interp_vec_curv2tel     - Does the spatial interpolation froma curvilinear grid 
!                           to a telescoping grid
!       xtrapfunc          - Extrapolation function
!     Cell-to-face interpolations 
!       interp_coef_cell2face - Cell-to-face interpolation coefficients
!       interp_scal_cell2face   - Scalar interpolation from cell centers to cell faces
!     Cell-to-node interpolations 
!       interp_coef_cell2node_invarea - Inverse area interpolation coefficients for a 
!                         scalar from cell centroids to nodes
!       interp_coef_cell2node_invdist - Inverse distance interpolation coefficients 
!                         scalar from cell centroids to nodes
!       interp_coef_cell2node_lstsqrs - Least squares interpolation coefficients for a 
!                         scalar from cell centroids to nodes
!       interp_scal_cell2node - Interpolates a scalar (1 array) from 
!                             cell centroids to nodes
!       interp_vec_cell2node  - Interpolates a vector (2 arrays) from 
!                             cell centroids to nodes
!     Node-to-cell interpolations
!       interp_scal_node2cell - Interpolates a scalar array from nodes to cell centers
!       interp_vec_node2cell - Interpolates a vector array from nodes to cell centers
!
! written by Alex Sanchez, USACE-CHL
!==================================================================================
    implicit none
    
contains

!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!  Begin 1D interpolations
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

!******************************************************   
    function interp1dlin(xp,yp,xi) result(yi)
! Simple 1D Linear Interpolation Function
! written by Alex Sanchez, USACE-CHL
!****************************************************** 
    use prec_def
    implicit none
    real(ikind),intent(in):: xp(2),yp(2),xi
    real(ikind):: fac,yi
    
    fac=(xi-xp(1))/(xp(2)-xp(1)) 
    fac=max(min(fac,1.0),0.0) !Avoids extrapolation
    yi=(1.0-fac)*yp(1)+fac*yp(2)
    
    return
    end function interp1dlin
    
!*********************************************************   
    subroutine interp_linear_fit(m,x,xi,inc,fac1,fac2)
! Fits a 1D Linear Interpolation Function
! yi = fac1*y(inc)+fac2*y(inc+1)
!
! Author: Alex Sanchez, USACE-CHL
!********************************************************* 
    use prec_def
    implicit none
    integer,intent(in):: m
    integer,intent(inout):: inc
    integer:: i
    real(ikind),intent(in):: x(m),xi
    real(ikind):: fac1,fac2
    
    !Avoid extrapolation
    if(xi<x(1))then
      fac1 = 1.0
      fac2 = 0.0
      inc = 1
      return
    elseif(xi>x(m))then
      fac1 = 0.0
      fac2 = 1.0
      inc = m-1
      return 
    endif
    
    !Find starting location so that x(inc)<=xi<=x(inc+np)
    if(inc>=1)then !Sequential search
      i = inc
      do inc=i,m-1 !Forwards search
        if(xi<=x(inc+1)) exit  
      enddo      
      if(xi<x(inc) .or. xi>x(inc+1))then !Check
        do inc=i,1,-1 !Backwards  search
          if(x(inc)<=xi) exit  
        enddo 
        if(xi<x(inc) .or. xi>x(inc+1))then !Check
          write(*,*) 'ERROR: Problem calculating starting location in plagr_fit'
          write(*,*) '  Using Bisection search'
          read(*,*)
          stop
        endif
      endif
    endif
    
    fac2 = (xi-x(inc))/(x(inc+1)-x(inc)) 
    fac2 = max(min(fac2,1.0),0.0) !Avoids extrapolation
    fac1 = 1.0-fac2
    
    return
    end subroutine interp_linear_fit

!******************************************************   
    function interp1dquad(xp,yp,xi) result(yi)
! Simple 1D Quadratic Interpolation Function
! Fits the equation
! y = b0+b1*(x-xp(1))+b2*(x-xp(1))*(x-xp(2))
! written by Alex Sanchez, USACE-CHL
!******************************************************   
    use prec_def
    implicit none
    real(ikind),intent(in):: xp(3),yp(3),xi
    real(ikind):: b1,b2,yi
    
    if(xi<=xp(1))then
       yi=yp(1)
       return
    elseif(xi>=xp(3))then
       yi=yp(3)
       return   
    endif
    b1=(yp(2)-yp(1))/(xp(2)-xp(1))
    b2=((yp(3)-yp(2))/(xp(3)-xp(2))-b1)/(xp(3)-xp(1))
    yi=yp(1)+b1*(xi-xp(1))+b2*(xi-xp(1))*(xi-xp(2))  !b0=yp(1)
    
    return
    end function interp1dquad

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!  End 1D interpolations
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!  Begin element interpolations
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

!**************************************************************************** 
    subroutine interp_coef_tri(xi,yi,xt,yt,w)
! Interpolation coefficients for a linear triangle. 
! Solves the equation for a 2D plane defined by the three points (xt(3),yt(3))
! using Cramer's Rule and calculates the interpolation coefficients w(3) 
! for a point at a point (xi,yi).
!
! Points are normalized from 0 to 1 in order to reduce precision errors.
!
! The method is more expensive than others but has the advantage that the 
! input triangle points (x(3),y(3)) do not need to sorted in any way.
!
! written by Alex Sanchez, USACE-CHL
!****************************************************************************   
    use prec_def
    implicit none
    !Input/Output
    real(ikind),intent(in) :: xt(3),yt(3),xi,yi
    real(ikind),intent(out) :: w(3)
    !Internal Variables
    real(ikind) :: d,xn(3),yn(3),xin,yin,xmin,ymin,xrange,yrange,sumw
    
    !Normalize to reduce precision errors
    xmin = minval(xt)
    ymin = minval(yt)
    xrange = maxval(xt)-xmin
    yrange = maxval(yt)-ymin
    xn = (xt-xmin)/xrange
    yn = (yt-ymin)/yrange
    xin = (xi-xmin)/xrange
    yin = (yi-ymin)/yrange
    
    !Calculate coefficients using Cramer's Rule
    d = xn(1)*(yn(2)-yn(3)) + xn(2)*(yn(3)-yn(1)) + xn(3)*(yn(1)-yn(2))
    w(1) = (xin*(yn(2)-yn(3)) + yin*(xn(3)-xn(2)) + xn(2)*yn(3) - xn(3)*yn(2))/d
    w(2) = (xin*(yn(3)-yn(1)) + yin*(xn(1)-xn(3)) + xn(3)*yn(1) - xn(1)*yn(3))/d
    w(3) = (xin*(yn(1)-yn(2)) + yin*(xn(2)-xn(1)) + xn(1)*yn(2) - xn(2)*yn(1))/d
    
    sumw = sum(w)
    w = w/sum(w)
    if(abs(sumw-1.0)>1.0e-5)then
      write(*,*) 'ERROR: Problem calculating interpolation coefficients'
      write(*,*) 'w = ',w
      write(*,*) 'sum = w(1:3) = ',sumw
      write(*,*) 'xt(1:3) = ',xt
      write(*,*) 'yt(1:3) = ',yt
      write(*,*) 'xi = ',xi
      write(*,*) 'yi = ',yi
      write(*,*) 'd = ',d      
      write(*,*) '  Press <enter> key to continue.'
      read(*,*)
      stop
    endif

    return
    end subroutine interp_coef_tri
    
!**************************************************************************** 
    subroutine interp_coef_quad(xi,yi,xq,yq,w)
! Interpolation coefficients for a bilinear quadrilateral. 
! Solves the equation for a 2D plane defined by the three points (xq(4),yt(4))
!
! Points are normalized from 0 to 1 in order to avoid precision errors.
!
! The method is more expensive than others but has the advantage that the 
! input quadrilatera points (xq(4),yq(4)) do not need to sorted in any way.
!
! written by Alex Sanchez, USACE-CHL
!****************************************************************************   
    use prec_def
    implicit none
    !Input/Output
    real(ikind),intent(in) :: xq(4),yq(4),xi,yi
    real(ikind),intent(out) :: w(4)
    !Internal Variables
    integer :: j
    real(ikind) :: d,wsum,xn(4),yn(4),xin,yin,xmin,ymin,xrange,yrange
    
    !Normalize to avoid precision errors
    xmin = minval(xq)
    ymin = minval(yq)
    xrange = maxval(xq)-xmin
    yrange = maxval(yq)-ymin
    xn = (xq-xmin)/xrange
    yn = (yq-ymin)/yrange
    xin = (xi-xmin)/xrange
    yin = (yi-ymin)/yrange
    
    !Inverse distance
    w=0.0
    wsum=0.0
    do j=1,4
      d=sqrt((xn(j)-xin)**2+(yn(j)-yin)**2)
      w(j)=1.0/max(d,1.0e-6)
      wsum=wsum+w(j)
    enddo
    w=w/wsum !Normalize
    wsum=sum(w)
    if(abs(wsum-1.0)>1.0e-4)then
      write(*,*) 'ERROR: Problem calculating quadrilateral interpolation coefficients'
      write(*,*) 'w(1:4) = ',w
      write(*,*) 'wsum = ',wsum
      write(*,*) 'xq(1:4) = ',xq
      write(*,*) 'yq(1:4) = ',yq
      write(*,*) 'xi = ',xi
      write(*,*) 'yi = ',yi
      write(*,*) '  Press <enter> key to continue.'
      read(*,*)
      stop
    endif    

    return
    end subroutine interp_coef_quad
    
!**************************************************************************** 
    subroutine interp_coef_poly(xi,yi,ns,np,xp,yp,w)
! Interpolation coefficients for polygon. 
!
! written by Alex Sanchez, USACE-CHL
!****************************************************************************   
    use prec_def
    implicit none
    !Input/Output
    integer,    intent(in) :: ns,np
    real(ikind),intent(in) :: xi,yi,xp(ns),yp(ns)
    real(ikind),intent(out) :: w(ns)
    !Internal Variables
    integer :: j
    real(ikind) :: d,wsum,xn(np),yn(np),xin,yin,xmin,ymin,xrange,yrange
    !real(ikind) :: Rx,Ry,Exx,Exy,Eyy,Emag,delx,dely,Gx,Gy !For Least-squares only
    
    !Normalize to avoid precision errors
    xmin = minval(xp(1:np))
    ymin = minval(yp(1:np))
    xrange = maxval(xp(1:np))-xmin
    yrange = maxval(yp(1:np))-ymin
    xn(1:np) = (xp(1:np)-xmin)/xrange
    yn(1:np) = (yp(1:np)-ymin)/yrange
    xin = (xi-xmin)/xrange
    yin = (yi-ymin)/yrange
    
    !Initialize weights
    w = 0.0
    wsum = 0.0
    
    !Inverse distance    
    do j=1,np
      d = sqrt((xn(j)-xin)**2+(yn(j)-yin)**2)
      w(j) = 1.0/max(d,1.0e-6)
      wsum = wsum + w(j)
    enddo
        
    !!Least-squares  
    !Rx=0.0; Ry=0.0
    !Exx=0.0; Exy=0.0; Eyy=0.0
    !do j=1,np
    !  delx=xn(j)-xin
    !  dely=yn(j)-yin
    !  Rx=Rx+delx
    !  Ry=Ry+dely
    !  Exx=Exx+delx*delx
    !  Eyy=Eyy+dely*dely
    !  Exy=Exy+delx*dely     
    !enddo !k
    !Emag=Exx*Eyy-Exy**2
    !Gx=(Exy*Ry-Exy*Rx)/Emag
    !Gy=(Exy*Rx-Exy*Ry)/Emag
    !do j=1,np
    !  delx=xn(j)-xin
    !  dely=yn(j)-yin 
    !  w(j)=1.0+Gx*delx+Gy*dely  
    !  wsum=wsum+w(j)
    !enddo !k
    
    !Normalize weights
    w(1:np) = w(1:np)/wsum !Normalize
    wsum = sum(w(1:np))
    
    if(abs(wsum-1.0)>1.0e-4)then
      write(*,*) 'ERROR: Problem calculating polygon interpolation coefficients'
      write(*,*) 'w(1:np) = ',w
      write(*,*) 'wsum = ',wsum
      write(*,*) 'xp(1:np) = ',xp(1:np)
      write(*,*) 'yp(1:np) = ',yp(1:np)
      write(*,*) 'xi = ',xi
      write(*,*) 'yi = ',yi
      write(*,*) '  Press <enter> key to continue.'
      read(*,*)
      stop
    endif
    
    return
    end subroutine interp_coef_poly 
    
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!  End element interpolations
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!  Begin cell-to-face interpolations
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<    

!**********************************************************************    
    subroutine interp_coef_tel2cart(nc,nD,mf,x0tel,y0tel,azthtel,&
                   xtel,ytel,dxtel,dytel,c2cf,idf,ncf, &            
                   ni,nj,x0cart,y0cart,azthcart,xi,yj,&
                   xtrapdist,nbnd,bndstr,intp,cntp)
! Calculates the interpolation coefficients and index from 
! a telescoping Cartesian grid (parent grid) to 
! a nonuniform Cartesian grid (child grid)
! All the points of the nonuniform Cartesian grid are used.
! written by Alex Sanchez, USACE-CHL
!**********************************************************************    
    use bnd_def, only: bnd_type
    use const_def, only: deg2rad,small,pi
    use geo_lib, only: intriangle
    use geo_def, only: mapid
    use prec_def
    implicit none
    !Telescoping Cartesian Grid (Parent)
    integer,    intent(in) :: nc          !# of active cells
    integer,    intent(in) :: nD          !# of active cells plus dummy cells
    integer,    intent(in) :: mf          !Max # of faces (nodes) per cell
    integer,    intent(in) :: idf(mf,nD)  !Direction of cell faces
    integer,    intent(in) :: ncf(nD)     !Number of cell faces
    integer,    intent(in) :: c2cf(mf,nD) !Connectivity table
    real(ikind),intent(in) :: x0tel,y0tel,azthtel
    real(ikind),intent(in) :: xtel(nD),ytel(nD)
    real(ikind),intent(in) :: dxtel(nD),dytel(nD)
    !Nonuniform Structured Cartesian Grid (Child)
    integer,    intent(in) :: ni,nj
    real(ikind),intent(in) :: x0cart,y0cart,azthcart
    real(ikind),intent(in) :: xi(ni),yj(nj)
    integer,    intent(inout) :: intp(0:4,ni,nj)
    real(ikind),intent(inout) :: cntp(4,ni,nj)    
    !Extrapolation distance
    real(ikind),intent(in) :: xtrapdist
    !Boundaries
    integer,intent(in) :: nbnd
    type(bnd_type),intent(in) :: bndstr(nbnd)
    !Internal variables
    integer :: i,j,ii,k,im,ibnd,nck,nck2,ij,iidistmin2
    integer :: ncart,ne,nn,nne,ncne,ncnw,ncse,ncsw
    real(ikind) :: dist,distx,disty,distmin,distmin2
    real(ikind) :: coscart,sincart,costel,sintel
    real(ikind) :: x_global,y_global,fx,fy,cfx,cfy
    real(ikind) :: xtri(3),ytri(3),wcoef(3),xp,yp
    real(ikind) :: xijtel(ni,nj),yijtel(ni,nj)    
    real(ikind), parameter :: toldist = 0.01 !m
    logical :: iscart(nD)    
    
    costel = cos(azthtel*deg2rad)
    sintel = sin(azthtel*deg2rad)
    coscart = cos(azthcart*deg2rad)
    sincart = sin(azthcart*deg2rad)
    
    !Coordinates of Cartesian (uniform or nonuniform) grid
    !with respect to Telescoping grid
    do i=1,ni
      do j=1,nj
        x_global = x0cart+xi(i)*coscart - yj(j)*sincart
        y_global = y0cart+xi(i)*sincart + yj(j)*coscart
        xijtel(i,j) =  (x_global-x0tel)*costel + (y_global-y0tel)*sintel
        yijtel(i,j) = -(x_global-x0tel)*sintel + (y_global-y0tel)*costel     
      enddo
    enddo
        
    !Determin regular cartesian cells
    do i=1,nD    
      iscart(i) = .false.   
    enddo
    ncart = 0
    do i=1,nc
      ne = 0
      nn = 0
      do k=1,ncf(i)
        if(idf(k,i)==1)then !North
          nn = nn + 1
        elseif(idf(k,i)==2)then !East
          ne = ne + 1
        endif
      enddo
      if(nn==1 .and. ne==1)then
        nne = 0
        ncse = c2cf(2,i)
        if(ncse==0) cycle
        do k=1,ncf(ncse)     
          if(idf(k,ncse)==1)then !North-East
            nne = nne + 1
          endif
        enddo
        ncne = c2cf(1,ncse)
        if(nne==1 .and. ncne>0)then
          iscart(i) = .true.
          ncart = ncart + 1
        endif
      endif
    enddo
!    write(*,*) 'Number of Regular Cartesian Cells: ',ncart
!    read(*,*)
    
    !Initialize
    intp=0
    cntp=0.0
    
di: do i=1,ni
dj:   do j=1,nj
        ij = j + (i-1)*nj
        if(mod(ij,50000)==0)then
          write(*,'(A,F6.2)') '    Percent complete: ',float(ij)*100/float(ni*nj)
        endif
        
        !if(i==79 .and. j==85)then
        !  continue
        !endif
        
        !Falls on the same point
        !distmin = 1.0e20 !Initialize
        do ii=1,nc
          distx = abs(xtel(ii)-xijtel(i,j))
          disty = abs(ytel(ii)-yijtel(i,j))  
          !dist = sqrt(distx*distx+disty*disty+small)
          !if(dist<distmin)then
          !  distmin = dist
          !endif
          if(distx<=0.01*dxtel(ii) .and. disty<=0.01*dytel(ii))then   
            !!Boundaries, find nearest boundary cell
            !do ibnd=1,nbnd
            !  do im=1,bndstr(ibnd)%ncells            
            !    if(bndstr(ibnd)%cells(im)==ii)then  
            !    endif
            !  enddo
            !enddo  
            !if(ii<=nc)then !Active cell
            !  intp(0,i,j) = 1
            !else !Inactive cell
            !  intp(0,i,j) = -1
            !endif
            intp(0,i,j) = 1
            intp(1,i,j) = ii
            cntp(1,i,j) = 1.0
            cycle dj
          endif
        enddo
        
        !Bilinear interpolation
        do ii=1,nc
          if(iscart(ii))then !Only for "Cartesian" cells surrounded by non-joint cells
            ncsw = ii
            ncnw = c2cf(1,ii)
            ncse = c2cf(2,ii)
            ncne = c2cf(1,ncse)       
            if(ncnw>nc .or. ncse>nc .or. ncne>nc) cycle !Only use active cells
            if(xtel(ncsw)<=xijtel(i,j) .and. ytel(ncsw)<=yijtel(i,j) .and. &
               xtel(ncse)>=xijtel(i,j) .and. ytel(ncnw)>=yijtel(i,j))then
              intp(0,i,j) = 4
              intp(1,i,j) = ncsw
              intp(2,i,j) = ncse
              intp(3,i,j) = ncne
              intp(4,i,j) = ncnw
              fx = (xijtel(i,j)-xtel(ncsw))/(xtel(ncse)-xtel(ncsw))
              fy = (yijtel(i,j)-ytel(ncsw))/(ytel(ncnw)-ytel(ncsw))
              fx = min(max(fx,0.0),1.0)
              fy = min(max(fy,0.0),1.0)
              cfx = 1.0 - fx
              cfy = 1.0 - fy                  
              cntp(1,i,j) = cfx*cfy
              cntp(2,i,j) = fx*cfy
              cntp(3,i,j) = fx*fy
              cntp(4,i,j) = cfx*fy       
              cycle dj        
            endif
          endif
        enddo
        
        !Use triangles for linear interpolation 
        do ii=1,nc
          distx = abs(xtel(ii)-xijtel(i,j))
          disty = abs(ytel(ii)-yijtel(i,j))            
          if(2.0*distx<=dxtel(ii)+1.0e-5 .and. 2.0*disty<=dytel(ii)+1.0e-5)then               
            xp = xijtel(i,j)
            yp = yijtel(i,j)
            xtri(1) = xtel(ii)
            ytri(1) = ytel(ii)
            !Find nearest inclosing triangle
            do k=1,ncf(ii)              
              nck = c2cf(k,ii)
              if(k<ncf(ii))then
                nck2 = c2cf(k+1,ii)
              else
                nck2 = c2cf(1,ii)
              endif  
              if(nck>nc .or. nck2>nc) cycle !Do not interpolate using land cells
              xtri(2) = xtel(nck)
              ytri(2) = ytel(nck)              
              xtri(3) = xtel(nck2)
              ytri(3) = ytel(nck2)        
              if(intriangle(xp,yp,xtri,ytri))then !wave point inside flow triangle
                call interp_coef_tri(xp,yp,xtri,ytri,wcoef) !Calc interp coefficients
                intp(0,i,j) = 3
                intp(1,i,j) = ii
                intp(2,i,j) = nck
                intp(3,i,j) = nck2
                cntp(1:3,i,j) = wcoef(1:3)
                cycle dj
              endif
            enddo !k
          endif           
        enddo !ii
        
        !Nearest cells
        distmin = 1.0e20 !Initialize
        do ii=1,nc
          distx = abs(xtel(ii)-xijtel(i,j))
          disty = abs(ytel(ii)-yijtel(i,j))  
          dist = sqrt(distx*distx+disty*disty+small)
          if(dist<distmin)then
          !if(distx<=0.1*dxtel(ii) .and. disty<=0.1*dytel(ii))then
            distmin = dist
            intp(1,i,j) = ii 
          endif
        enddo
        !Boundaries, find nearest boundary cell
        do ibnd=1,nbnd
          do im=1,bndstr(ibnd)%ncells            
            if(bndstr(ibnd)%cells(im)==intp(1,i,j))then
              intp(0,i,j) = 0 !Indicates extrapolation
              cntp(1,i,j) = xtrapfunc(distmin,xtrapdist)
              cntp(2:4,i,j) = 0.0  
              cycle dj
            endif  
          enddo !im-flow-cell
        enddo !ibnd   
        distmin2 = 1.0e20 !Initialize
        iidistmin2 = 0
        do ii=nc+1,nD
          distx = abs(xtel(ii)-xijtel(i,j))
          disty = abs(ytel(ii)-yijtel(i,j))  
          dist = sqrt(distx*distx+disty*disty+small)
          if(dist<distmin2)then
            distmin2 = dist
            iidistmin2 = ii 
          endif
        enddo
        !if(distmin<distmin2)then
        !  intp(0,i,j) = 0 !Indicates extrapolation
        !  cntp(1,i,j) = xtrapfunc(distmin,xtrapdist)
        !  cntp(2:4,i,j) = 0.0
        !else
          intp(0,i,j) = -1 !Indicates inactive cell
          intp(1,i,j) = iidistmin2
          cntp(1,i,j) = xtrapfunc(distmin2,xtrapdist)        
          cntp(2:4,i,j) = 0.0    
        !endif
        
        !!Boundaries, find nearest boundary cell
        !distmin = 1.0e20 !Initialize
        !do ibnd=1,nbnd
        !  do im=1,bndstr(ibnd)%ncells
        !    ii=bndstr(ibnd)%cells(im)
        !    distx=abs(xtel(ii)-xijtel(i,j))
        !    disty=abs(ytel(ii)-yijtel(i,j))  
        !    dist=sqrt(distx**2+disty**2+small)  
        !    if(dist<distmin)then
        !      distmin=dist
        !      intp(1,i,j)=ii 
        !    endif  
        !  enddo !im-flow-cell
        !enddo !ibnd
        !!if(intp(1,i,j)<=0)then
        !!  continue
        !!endif
        !intp(0,i,j)=0 !Indicates extrapolation (values are not replaced by interpolation)
        !cntp(1,i,j)=xtrapfunc(distmin,xtrapdist)
        !cntp(2:4,i,j)=0.0  
      enddo dj !j-wave
    enddo di !i-wave
    write(*,'(A,F6.2)') '    Percent complete: ',100.0
    
    return
    end subroutine interp_coef_tel2cart
                   
!***********************************************************************
    subroutine interp_scal_tel2cart(nD,vartel,iwettel, &
                   ni,nj,intp,cntp,varcart,valdry,iextrap)
! Interpolate from a Telescoping grid to a Nonuniform Cartesian Grid
!
! iextrap = Specifies the boundary treatment
!         0 - No extrapolation (leave input value unchanged)
!         1 - Extrapolate to zero using distance function
!         2 - Extrapolate to original values using distance function
!         3 - Extrapolate out infinitely using neighest neighbor
!
! Authors: Mingliang, Z., Alex Sanchez, Weiming, Wu
!***********************************************************************
    use prec_def
    implicit none    
    !Telescoping Cartesian Grid (Parent)
    integer,    intent(in) :: nD          !# of active cells plus dummy cells
    real(ikind),intent(in) :: vartel(nD)
    integer,    intent(in) :: iwettel(nD)
    !Nonuniform Structured Cartesian Grid (Child)
    integer,    intent(in) :: ni,nj
    integer,    intent(in) :: intp(0:4,ni,nj)
    real(ikind),intent(in) :: cntp(4,ni,nj)    
    real(ikind),intent(inout) :: varcart(ni,nj)
    !Options
    real(ikind),intent(in) :: valdry
    integer,    intent(in) :: iextrap
    !Internal
    integer:: i,ii,j,k
    real(ikind) :: csum
    real(ikind), parameter :: toldry = 0.001
    
!$OMP PARALLEL DO PRIVATE(i,ii,j,k,csum)      
    do j=1,nj
      do i=1,ni
        if(intp(0,i,j)>0)then !Not extrapolation point or inactive cell
          csum = 0.0
          varcart(i,j) = 0.0
          do k=1,intp(0,i,j)
            ii = intp(k,i,j)
            if(iwettel(ii)==1)then !Wet
              varcart(i,j) = varcart(i,j) + vartel(ii)*cntp(k,i,j)
              csum = csum + cntp(k,i,j)
            endif
          enddo
          if(csum>1.0e-6)then
            varcart(i,j) = varcart(i,j)/csum
          else    
            varcart(i,j) = valdry
          endif
        elseif(intp(0,i,j)==0 .and. iextrap>0)then !intp(1,i,j)==0 -> Extrapolation point
          ii = intp(1,i,j)
          if(iwettel(ii)==1)then !Wet  
            select case(iextrap)
            case(1) !Extrapolate to zero   
              varcart(i,j) = vartel(ii)*cntp(1,i,j)  
            case(2) !Extrapolate to original
              varcart(i,j) = (1.0-cntp(1,i,j))*varcart(i,j) + vartel(ii)*cntp(1,i,j) 
            case(3) !Extrapolate out infinitly using neighest values
              varcart(i,j) = vartel(ii)
            end select  
          else
            varcart(i,j) = valdry  
          endif   
        !else !if(intp(0,i,j)<=-1 .or. iextrap==1)then !Leave unchanged
        endif
      enddo
    enddo
!$OMP END PARALLEL DO   

    return 
    end subroutine interp_scal_tel2cart  
                   
!***********************************************************************
    subroutine interp_vec_tel2cart(nD,vecxtel,vecytel,iwettel, &
                   ni,nj,intp,cntp,vecxcart,vecycart,valdry,iextrap)
! Interpolate from a Telescoping grid to a Nonuniform Cartesian Grid
!
! Input:
!   vecxcart,vecycart = Vector on Telescoping Grid
!   iwettel = Integer equal to 1 for wet cells
!   valdry = value used at dry cells
!   iextrap = Specifies the boundary treatment
!           0 - No extrapolation (leave input value unchanged)
!           1 - Extrapolate to zero using distance function
!           2 - Extrapolate to original values using distance function
!           3 - Extrapolate out infinitely using neighest neighbor
!
! Output:
!   vecxcart,vecycart = Vector on Nonuniform Cartesian Grid
!
! Author: Alex Sanchez, USACE-CHL
!***********************************************************************
    use prec_def
    implicit none    
    !Telescoping Cartesian Grid (Parent)
    integer,    intent(in) :: nD          !# of active cells plus dummy cells
    real(ikind),intent(in) :: vecxtel(nD),vecytel(nD)
    integer,    intent(in) :: iwettel(nD)
    !Nonuniform Structured Cartesian Grid (Child)
    integer,    intent(in) :: ni,nj
    integer,    intent(inout) :: intp(0:4,ni,nj)
    real(ikind),intent(inout) :: cntp(4,ni,nj)    
    real(ikind),intent(inout):: vecxcart(ni,nj),vecycart(ni,nj)
    !Options
    real(ikind),intent(in) :: valdry
    integer,    intent(in) :: iextrap
    !Internal
    integer:: i,ii,j,k
    real(ikind) :: csum
    real(ikind), parameter :: toldry = 0.001
    
!$OMP PARALLEL DO PRIVATE(i,ii,j,k,csum)      
     do j=1,nj 
       do i=1,ni
        if(intp(0,i,j)>0)then !Not extrapolation point
          csum = 0.0
          vecxcart(i,j) = 0.0
          vecycart(i,j) = 0.0
          do k=1,intp(0,i,j)
            ii = intp(k,i,j)
            if(iwettel(ii)==1)then !Wet
              vecxcart(i,j) = vecxcart(i,j) + vecxtel(ii)*cntp(k,i,j)
              vecycart(i,j) = vecycart(i,j) + vecytel(ii)*cntp(k,i,j)
              csum = csum + cntp(k,i,j)
            endif
          enddo
          if(csum>1.0e-6)then
            vecxcart(i,j) = vecxcart(i,j)/csum
            vecycart(i,j) = vecycart(i,j)/csum
          else    
            vecxcart(i,j) = valdry
            vecycart(i,j) = valdry
          endif
        elseif(intp(0,i,j)==0 .and. iextrap>0)then !intp(1,i,j)==0 -> Extrapolation point
          ii = intp(1,i,j)
          if(iwettel(ii)==1)then !Wet  
            select case(iextrap)
            case(1) !Extrapolate to zero   
              vecxcart(i,j) = vecxtel(ii)*cntp(1,i,j)  
              vecycart(i,j) = vecytel(ii)*cntp(1,i,j)  
            case(2) !Extrapolate to original
              vecxcart(i,j) = (1.0-cntp(1,i,j))*vecxcart(i,j) + vecxtel(ii)*cntp(1,i,j) 
              vecycart(i,j) = (1.0-cntp(1,i,j))*vecycart(i,j) + vecytel(ii)*cntp(1,i,j) 
            case(3) !Extrapolate out infinitly using neighest values
              vecxcart(i,j) = vecxtel(ii)
              vecycart(i,j) = vecytel(ii)
            end select  
          else
            vecxcart(i,j) = valdry  
            vecycart(i,j) = valdry
          endif   
        !else iextrap==1  - Leave unchanged
        endif
      enddo
    enddo
!$OMP END PARALLEL DO   

    return 
    end subroutine interp_vec_tel2cart                    
    
!******************************************************************************    
    subroutine interp_coef_cart2tel(ni,nj,x0cart,y0cart,azthcart,xi,yj,dxwav,dywav,&
                   nc,nD,x0tel,y0tel,azthtel,xtel,ytel,xtrapdist,intp,cntp)
! Calculates the interpolation coefficients and index from a 
! nonuniform Cartesian grid (parent grid) to a
! telescoping Cartesian grid (child grid). 
! All the points on the telescoping grid are used.
! written by Alex Sanchez, USACE-CHL
!******************************************************************************
    use const_def, only: deg2rad,small,pi
    use geo_def, only: mapid
    use prec_def
    implicit none
    !Nonuniform Structured Cartesian Grid (Parent)
    integer,    intent(in) :: ni,nj
    real(ikind),intent(in) :: x0cart,y0cart,azthcart
    real(ikind),intent(in) :: xi(ni),yj(nj)
    real(ikind),intent(in) :: dxwav(ni),dywav(nj)
    !Telescoping Cartesian Grid (Child)
    integer,    intent(in) :: nc,nD    
    real(ikind),intent(in) :: x0tel,y0tel,azthtel
    real(ikind),intent(in) :: xtel(nD),ytel(nD)
    integer,    intent(inout) :: intp(2,nD)
    real(ikind),intent(inout) :: cntp(4,nD)
    !Extrapolation distance
    real(ikind),intent(in) :: xtrapdist
    !Internal variables
    integer :: i,ii,j
    real(ikind) :: dist,x_global,y_global
    real(ikind) :: coscart,sincart,costel,sintel
    real(ikind) :: fx,cfx,fy,cfy,distmin,distx,disty
    real(ikind) :: xflwav(nD),yflwav(nD)
    real(ikind),parameter :: toldist = 0.001 ![m]
    
    costel = cos(azthtel*deg2rad)
    sintel = sin(azthtel*deg2rad)
    coscart = cos(azthcart*deg2rad)
    sincart = sin(azthcart*deg2rad)
    
    !Coordinates with respect to wave cooridinate system
    do i=1,nD
      x_global=x0tel+xtel(i)*costel-ytel(i)*sintel
      y_global=y0tel+xtel(i)*sintel+ytel(i)*costel
      xflwav(i)= (x_global-x0cart)*coscart+(y_global-y0cart)*sincart
      yflwav(i)=-(x_global-x0cart)*sincart+(y_global-y0cart)*coscart      
    enddo
    
    !Initialize variables
    intp=0
    cntp=0.0    
            
d1: do ii=1,nc
      if(mod(ii,10000)==0)then
        write(*,'(A,F6.2)') '    Percent complete: ',float(ii)*100/float(nc)
      endif
      
      !if(mapid(ii)==9040)then
      !  continue
      !endif
      
      !Search or overlapping points
      do i=1,ni
        do j=1,nj
          distx=abs(xflwav(ii)-xi(i))
          disty=abs(yflwav(ii)-yj(j))
          dist=sqrt(distx**2+disty**2+small)      
          if(dist<toldist)then !Points fall on one another
            intp(1,ii)=i
            intp(2,ii)=j
            cntp(1,ii)=1.0
            cntp(2:4,ii)=0.0
            cycle d1
          endif                                
         enddo
      enddo  
      
      !Bilinear interpolation   
      do i=1,ni-1
        do j=1,nj-1
          if((xflwav(ii)>=xi(i) .and. xflwav(ii)<=xi(i+1)) .and. &
             (yflwav(ii)>=yj(j) .and. yflwav(ii)<=yj(j+1)))then !Flow point is between 4 wave points    
            intp(1,ii)=i
            intp(2,ii)=j        
            fx=(xflwav(ii)-xi(i))/(xi(i+1)-xi(i))
            fy=(yflwav(ii)-yj(j))/(yj(j+1)-yj(j))
            fx=min(max(fx,0.0),1.0)
            fy=min(max(fy,0.0),1.0)
            cfx=1.0-fx
            cfy=1.0-fy
            cntp(1,ii)=cfx*cfy
            cntp(2,ii)=fx*cfy
            cntp(3,ii)=fx*fy
            cntp(4,ii)=cfx*fy
            cycle d1
          endif                                
         enddo
      enddo
      
      do i=1,ni
        do j=1,nj
          distx=abs(xflwav(ii)-xi(i))
          disty=abs(yflwav(ii)-yj(j))
          dist=sqrt(distx**2+disty**2+small) 
          if(2.0*distx<=dxwav(i).and.2.0*disty<=dywav(j))then !Flow point is inside wave cell     
            intp(1,ii)=i
            intp(2,ii)=j  
            cntp(1,ii)=1.0
            cntp(2:4,ii)=0.0
            cycle d1
          endif                  
        enddo
      enddo
      
      !Extrapolation
      distmin = 1.0e25
      do i=1,ni,ni-1
        do j=1,nj
          dist=sqrt((xflwav(ii)-xi(i))**2+(yflwav(ii)-yj(j))**2+small)   
          if(dist<distmin)then  
            distmin=dist
            intp(1,ii)=-i !Negative value indicates extrapolation
            intp(2,ii)=-j !Negative value indicates extrapolation
          endif
        enddo !j
      enddo !i
      do i=1,ni
        do j=1,nj,nj-1
          dist=sqrt((xflwav(ii)-xi(i))**2+(yflwav(ii)-yj(j))**2+small)   
          if(dist<distmin)then  
            distmin=dist
            intp(1,ii)=-i !Negative value indicates extrapolation
            intp(2,ii)=-j !Negative value indicates extrapolation
          endif
        enddo !j
      enddo !i    
      cntp(1,ii)=xtrapfunc(distmin,xtrapdist)
        cntp(2:4,ii)=0.0
    enddo d1
    write(*,'(A,F6.2)') '    Percent complete: ',100.0
    
    return
    end subroutine interp_coef_cart2tel
     
!******************************************************************************    
    subroutine interp_scal_cart2tel(ni,nj,varcart,iwetcart,intp,cntp,&
                   nc,nD,vartel,valdry,iextrap)
! Calculates the interpolation coefficients and index from a 
! nonuniform Cartesian grid (parent grid) to a
! telescoping Cartesian grid (child grid). 
! All the points on the telescoping grid are used.
!
! iextrap = Specifies the boundary treatment
!         0 - No extrapolation (leave input value unchanged)
!         1 - Extrapolate to zero using distance function
!         2 - Extrapolate to original values using distance function
!         3 - Extrapolate out infinitely using neighest neighbor
!
! written by Alex Sanchez, USACE-CHL
!******************************************************************************
    use prec_def
    implicit none
    !Nonuniform Structured Cartesian Grid (Parent)
    integer,    intent(in) :: ni,nj
    real(ikind),intent(in) :: varcart(ni,nj)
    integer,    intent(in) :: iwetcart(ni,nj)
    !Telescoping Cartesian Grid (Child)
    integer,    intent(in) :: nc,nD
    real(ikind),intent(out) :: vartel(nD) 
    integer,    intent(inout) :: intp(2,nD)
    real(ikind),intent(inout) :: cntp(4,nD)
    !Options
    real(ikind),intent(in) :: valdry
    integer,    intent(in) :: iextrap
    !Internal variables
    integer :: i,ii,i1,j,j1
    real(ikind) :: csum
    real(ikind),parameter :: toldry = 0.001     
    
!$OMP PARALLEL DO PRIVATE(i,ii,i1,j,j1,csum) 
    do ii=1,nc
      i=intp(1,ii)
      j=intp(2,ii)
      i1=min(i+1,ni)
      j1=min(j+1,nj)
      if(i>=1 .and. j>=1)then !No extrapolation
        vartel(ii) = 0.0
        csum = 0.0
        if(iwetcart(i ,j )==1)then !wet
          vartel(ii) = vartel(ii) + varcart(i ,j )*cntp(1,ii)
          csum = csum + cntp(1,ii)
        endif
        if(iwetcart(i1,j )==1)then !wet
          vartel(ii) = vartel(ii) + varcart(i1,j )*cntp(2,ii)
          csum = csum + cntp(2,ii)
        endif
        if(iwetcart(i1,j1)==1)then !wet
          vartel(ii) = vartel(ii) + varcart(i1,j1)*cntp(3,ii)
          csum = csum + cntp(3,ii)
        endif
        if(iwetcart(i ,j1)==1)then !wet
          vartel(ii) = vartel(ii) + varcart(i ,j1)*cntp(4,ii)
          csum = csum + cntp(4,ii)
        endif
        if(csum>1.0e-6)then
          vartel(ii) = vartel(ii)/csum
        else
          vartel(ii) = valdry
        endif
      elseif(i<1 .and. j<1 .and. iextrap>0)then !Extrapolation
        i1 = -i; j1 = -j
        if(iwetcart(i1,j1)==1)then !wet
          select case(iextrap)
          case(1) !Extrapolate to zero   
            vartel(ii) = varcart(i1,j1)*cntp(1,ii)
          case(2) !Extrapolate to original
            vartel(ii) = (1.0-cntp(1,ii))*vartel(ii) + varcart(i1,j1)*cntp(1,ii)  
          case(3) !Extrapolate out infinitly using neighest values
            vartel(ii) = varcart(i1,j1)
          end select  
        else
          vartel(ii) = valdry
        endif
      !else !Value not replaced
      endif
    enddo
!$OMP END PARALLEL DO

    return
    end subroutine interp_scal_cart2tel
    
!******************************************************************************    
    subroutine interp_vec_cart2tel(ni,nj,vecxcart,vecycart,iwetcart,intp,cntp,&
                   nc,nD,vecxtel,vecytel,valdry,iextrap)
! Calculates the interpolation coefficients and index from a 
! nonuniform Cartesian grid (parent grid) to a
! telescoping Cartesian grid (child grid). 
! All the points on the telescoping grid are used.
!
! iextrap = Specifies the boundary treatment
!         0 - No extrapolation (leave input value unchanged)
!         1 - Extrapolate to zero using distance function
!         2 - Extrapolate to original values using distance function
!         3 - Extrapolate out infinitely using neighest neighbor
!
! written by Alex Sanchez, USACE-CHL
!******************************************************************************
    use prec_def
    implicit none
    !Nonuniform Structured Cartesian Grid (Parent)
    integer,    intent(in) :: ni,nj
    real(ikind),intent(in) :: vecxcart(ni,nj),vecycart(ni,nj)
    integer,    intent(in) :: iwetcart(ni,nj)
    !Telescoping Cartesian Grid (Child)
    integer,    intent(in) :: nc,nD
    real(ikind),intent(out) :: vecxtel(nD),vecytel(nD) 
    integer,    intent(inout) :: intp(2,nD)
    real(ikind),intent(inout) :: cntp(4,nD)
    !Options
    real(ikind),intent(in) :: valdry
    integer,    intent(in) :: iextrap
    !Internal variables
    integer :: i,ii,i1,j,j1
    real(ikind) :: csum
    real(ikind),parameter :: toldry = 0.001     
    
!$OMP PARALLEL DO PRIVATE(i,ii,i1,j,j1,csum) 
    do ii=1,nc
      i=intp(1,ii)
      j=intp(2,ii)
      i1=min(i+1,ni)
      j1=min(j+1,nj)
      if(i>=1 .and. j>=1)then !No extrapolation
        vecxtel(ii) = 0.0
        vecytel(ii) = 0.0
        csum = 0.0
        if(iwetcart(i ,j )==1)then !wet
          vecxtel(ii) = vecxtel(ii) + vecxcart(i ,j )*cntp(1,ii)
          vecytel(ii) = vecytel(ii) + vecycart(i ,j )*cntp(1,ii)
          csum = csum + cntp(1,ii)
        endif        
        if(iwetcart(i1,j )==1)then !wet    
          vecxtel(ii) = vecxtel(ii) + vecxcart(i1,j )*cntp(2,ii)
          vecytel(ii) = vecytel(ii) + vecycart(i1,j )*cntp(2,ii)
          csum = csum + cntp(2,ii)
        endif
        if(iwetcart(i1,j1)==1)then !wet  
          vecxtel(ii) = vecxtel(ii) + vecxcart(i1,j1)*cntp(3,ii)
          vecytel(ii) = vecytel(ii) + vecycart(i1,j1)*cntp(3,ii)
          csum = csum + cntp(3,ii)
        endif
        if(iwetcart(i ,j1)==1)then !wet
          vecxtel(ii) = vecxtel(ii) + vecxcart(i ,j1)*cntp(4,ii)
          vecytel(ii) = vecytel(ii) + vecycart(i ,j1)*cntp(4,ii)
          csum = csum + cntp(4,ii)
        endif
        if(csum>1.0e-6)then
          vecxtel(ii) = vecxtel(ii)/csum
          vecytel(ii) = vecytel(ii)/csum
        else
          vecxtel(ii) = valdry
          vecytel(ii) = valdry
        endif
      elseif(i<1 .and. j<1 .and. iextrap>0)then !Extrapolation
        i1 = -i; j1 = -j
        if(iwetcart(i1,j1)==1)then !wet
          select case(iextrap)
          case(1) !Extrapolate to zero   
            vecxtel(ii) = vecxcart(i1,j1)*cntp(1,ii)
            vecytel(ii) = vecycart(i1,j1)*cntp(1,ii)
          case(2) !Extrapolate to original
            vecxtel(ii) = (1.0-cntp(1,ii))*vecxtel(ii) + vecxcart(i1,j1)*cntp(1,ii)  
            vecytel(ii) = (1.0-cntp(1,ii))*vecytel(ii) + vecycart(i1,j1)*cntp(1,ii)  
          case(3) !Extrapolate out infinitly using neighest values
            vecxtel(ii) = vecxcart(i1,j1)
            vecytel(ii) = vecycart(i1,j1)
          end select  
        else
          vecxtel(ii) = valdry
          vecytel(ii) = valdry
        endif
      !else !Value not replaced
      endif
    enddo
!$OMP END PARALLEL DO

    return
    end subroutine interp_vec_cart2tel    
    
!******************************************************************************    
    subroutine interp_coef_cart2pts(ni,nj,x0cart,y0cart,azthcart,xi,yj,dxwav,dywav,&
                   npts,nintp,xpts,ypts,xtrapdist,intp,cntp)
! Calculates the interpolation coefficients and index from a 
! nonuniform Cartesian grid (parent grid) to a
! telescoping Cartesian grid (child grid). 
! All the points on the telescoping grid are used.
! written by Alex Sanchez, USACE-CHL
!******************************************************************************
    use const_def, only: deg2rad,small,pi
    use prec_def
    implicit none
    !Nonuniform Structured Cartesian Grid (Parent)
    integer,    intent(in) :: ni,nj
    real(ikind),intent(in) :: x0cart,y0cart,azthcart
    real(ikind),intent(in) :: xi(ni),yj(nj)
    real(ikind),intent(in) :: dxwav(ni),dywav(nj)
    !Child Points
    integer,    intent(in) :: npts,nintp
    real(ikind),intent(in) :: xpts(npts),ypts(npts)
    integer,    intent(inout) :: intp(2,npts)   !Cell Id's to interpolation (parent) grid to interpolate from
    real(ikind),intent(inout) :: cntp(4,npts)     !Interpolation coefficients from parent to child grid
    !Extrapolation distance
    real(ikind),intent(in) :: xtrapdist
    !Internal variables
    integer :: i,ii,j
    real(ikind) :: coscart,sincart,dist
    real(ikind) :: fx,cfx,fy,cfy,distmin,distx,disty
    real(ikind) :: xptscart(npts),yptscart(npts) !Scatter points relative to Cartesian grid
    real(ikind),parameter :: toldist = 0.001 ![m]
    
    coscart = cos(azthcart*deg2rad)
    sincart = sin(azthcart*deg2rad)
    
    !Points relative position with respect to Cartesian grid
    do ii=1,nintp
      xptscart(ii) =  (xpts(ii)-x0cart)*coscart + (ypts(ii)-y0cart)*sincart
      yptscart(ii) = -(xpts(ii)-x0cart)*sincart + (ypts(ii)-y0cart)*coscart    
    enddo    
    
    !!open(224,file='Coord_flwav.dat')
    !!do ii=1,nintp
    !!  write(224,*) xptscart(ii),yptscart(ii),0.0 !,depwave(i,j)
    !!enddo
    !!close(224) 
    
    !Initialize variables
    intp = 0
    cntp = 0.0    
            
d1: do ii=1,nintp
      if(mod(ii,2000)==0)then
        write(*,'(A,F6.2)') '    Percent complete: ',float(ii)*100/float(nintp)
      endif
      
      !Search or overlapping points
      do i=1,ni
        do j=1,nj
          distx = abs(xptscart(ii)-xi(i))
          disty = abs(yptscart(ii)-yj(j))
          dist = sqrt(distx**2+disty**2+small)      
          if(dist<toldist)then !Points fall on one another
            intp(1,ii) = i
            intp(2,ii) = j
            cntp(1,ii) = 1.0
            cntp(2:4,ii) = 0.0
            cycle d1
          endif
         enddo
      enddo   
      
      !Bilinear interpolation   
      do i=1,ni-1
        do j=1,nj-1
          if((xptscart(ii)>=xi(i) .and. xptscart(ii)<=xi(i+1)) .and. &
             (yptscart(ii)>=yj(j) .and. yptscart(ii)<=yj(j+1)))then !Flow point is between 4 wave points    
            intp(1,ii) = i
            intp(2,ii) = j
            fx = (xptscart(ii)-xi(i))/(xi(i+1)-xi(i))
            fy = (yptscart(ii)-yj(j))/(yj(j+1)-yj(j))
            fx = min(max(fx,0.0),1.0)
            fy = min(max(fy,0.0),1.0)
            cfx = 1.0 - fx
            cfy = 1.0 - fy
            cntp(1,ii) = cfx*cfy
            cntp(2,ii) = fx*cfy
            cntp(3,ii) = fx*fy
            cntp(4,ii) = cfx*fy
            cycle d1
          endif
         enddo
      enddo
      
      !Falls inside a cell
      do i=1,ni
        do j=1,nj
          distx = abs(xptscart(ii)-xi(i))
          disty = abs(yptscart(ii)-yj(j))
          dist = sqrt(distx**2+disty**2+small) 
          if(2.0*distx<=dxwav(i).and.2.0*disty<=dywav(j))then !Flow point is inside wave cell     
            intp(1,ii) = i
            intp(2,ii) = j  
            cntp(1,ii) = 1.0
            cntp(2:4,ii) = 0.0
            cycle d1
          endif
        enddo
      enddo
      
      !Extrapolation
      distmin = 1.0e25
      do i=1,ni,ni-1
        do j=1,nj
          dist=sqrt((xptscart(ii)-xi(i))**2+(yptscart(ii)-yj(j))**2+small)   
          if(dist<distmin)then  
            distmin = dist
            intp(1,ii) = -i !Negative value indicates extrapolation
            intp(2,ii) = -j !Negative value indicates extrapolation
          endif
        enddo !j
      enddo !i
      do i=1,ni
        do j=1,nj,nj-1
          dist = sqrt((xptscart(ii)-xi(i))**2+(yptscart(ii)-yj(j))**2+small)   
          if(dist<distmin)then  
            distmin = dist
            intp(1,ii) = -i !Negative value indicates extrapolation
            intp(2,ii) = -j !Negative value indicates extrapolation
          endif
        enddo !j
      enddo !i    
      cntp(1,ii) = xtrapfunc(distmin,xtrapdist)
        cntp(2:4,ii) = 0.0
    enddo d1 !ii  
    
    return
    end subroutine interp_coef_cart2pts                   
    
!**********************************************************************    
    subroutine interp_coef_tel2tel(nc,nD,mf,x0tel,y0tel,azthtel, &
                  xtel,ytel,dxtel,dytel,c2cf,idf,ncf, &            
                  nc2,nD2,x0tel2,y0tel2,azthtel2, &
                  xtel2,ytel2,xtrapdist,&
                  nintrpcells,iintpcells,intp,cntp)
! Calculates the interpolation coefficients from a telescoping Cartesian
! grid to another telescoping Cartesian grid
! written by Alex Sanchez, USACE-CHL      
!**********************************************************************
    use const_def, only: deg2rad,small,pi
    use geo_lib, only: intriangle
    use prec_def
    implicit none
    !Parent/Interpolation Telescoping Cartesian grid (Parent)
    integer,       intent(in) :: nc,nD,mf
    integer,       intent(in) :: idf(mf,nD)  !Direction of cell faces
    integer,       intent(in) :: ncf(nD)     !Number of cell faces
    integer,       intent(in) :: c2cf(mf,nD) !Connectivity table
    real(ikind),   intent(in) :: x0tel,y0tel,azthtel
    real(ikind),   intent(in) :: xtel(nD),ytel(nD)
    real(ikind),   intent(in) :: dxtel(nD),dytel(nD)
    !Child/Base Telescoping Cartesian Grid (Child)
    integer,    intent(in)    :: nc2,nD2
    real(ikind),intent(in)    :: x0tel2,y0tel2,azthtel2
    real(ikind),intent(in)    :: xtel2(nD2),ytel2(nD2)
    integer,    intent(in)    :: nintrpcells             !Number of cells to interpolate
    integer,    intent(in)    :: iintpcells(nintrpcells) !Cell Id's of child grid to interpolate
    integer,    intent(inout) :: intp(0:4,nintrpcells)   !Cell Id's to interpolation (parent) grid to interpolate from
    real(ikind),intent(inout) :: cntp(4,nintrpcells)     !Interpolation coefficients from parent to child grid
    !Extrapolation distance
    real(ikind),intent(in) :: xtrapdist
    !Internal variables
    integer :: i,ii,j,k,nck,nck2
    integer :: ncart,ne,nn,nne,ncse,ncne,ncsw,ncnw
    real(ikind) :: dist,distx,disty,distmin
    real(ikind) :: costel,sintel,costel2,sintel2
    real(ikind) :: x_global,y_global,fx,fy,cfx,cfy
    real(ikind) :: xtri(3),ytri(3),wcoef(3),xp,yp
    real(ikind) :: xijtel(nc2),yijtel(nc2)    
    real(ikind),parameter :: toldist = 0.001 ![m]
    logical :: iscart(nD)
        
    costel = cos(azthtel*deg2rad)
    sintel = sin(azthtel*deg2rad)
    costel2 = cos(azthtel2*deg2rad)
    sintel2 = sin(azthtel2*deg2rad)
    
    !Coordinates of child (current) telescoping grid
    !with respect to parent telescoping grid
    do j=1,nintrpcells
      i = iintpcells(j)
      x_global  =  x0tel2 + xtel2(i)*costel2 - ytel2(i)*sintel2
      y_global  =  y0tel2 + xtel2(i)*sintel2 + ytel2(i)*costel2
      xijtel(i) =  (x_global-x0tel)*costel + (y_global-y0tel)*sintel
      yijtel(i) = -(x_global-x0tel)*sintel + (y_global-y0tel)*costel
    enddo
    
    !open(43,file='Parent.xyz')
    !do i=1,nD !Full grid including land cells
    !  write(43,*) xtel(i),ytel(i),0.0
    !enddo
    !close(43)
    !
    !open(43,file='Child.xyz')
    !do j=1,nintrpcells
    !  i=iintpcells(j)
    !  write(43,*) xijtel(i),yijtel(i),0.0
    !enddo
    !close(43)
    
    !Determin regular cartesian cells of parent telescopting grid
    do i=1,nc
      iscart(i) = .false.   
    enddo
    ncart = 0
    do i=1,nc
      ne = 0
      nn = 0
      do k=1,ncf(i)
        if(idf(k,i)==1)then !North
          nn = nn + 1
        elseif(idf(k,i)==2)then !East
          ne = ne + 1
        endif
      enddo
      if(nn==1 .and. ne==1)then
        nne = 0
        ncse = c2cf(2,i)
        if(ncse==0) cycle
        do k=1,ncf(ncse)     
          if(idf(k,ncse)==1)then !North-East
            nne = nne + 1
          endif
        enddo
        ncne = c2cf(1,ncse)
        if(nne==1 .and. ncne>0)then
          iscart(i) = .true.
          ncart = ncart + 1
        endif
      endif
    enddo
!    write(*,*) 'Number of Regular Cartesian Cells: ',ncart
!    read(*,*)
    
    !Initialize
    intp = 0
    cntp = 0.0
    
d1: do j=1,nintrpcells
      i = iintpcells(j)
      
      !Falls on the same point   
      do ii=1,nc                
        distx = abs(xtel(ii)-xijtel(i))
        disty = abs(ytel(ii)-yijtel(i))  
        dist = sqrt(distx*distx+disty*disty+small)                                                                   
        if(dist<toldist)then 
          intp(0,j) = 1
          intp(1,j) = ii
          cntp(1,j) = 1.0
          cntp(2:4,j) = 0.0
          cycle d1
        endif  
      enddo
        
      !Bilinear interpolation
      do ii=1,nc                 
        if(iscart(ii))then 
          ncsw = ii
          ncnw = c2cf(1,ii)
          ncse = c2cf(2,ii)
          ncne = c2cf(1,ncse)
          if(ncnw>nc .or. ncse>nc .or. ncne>nc) cycle
          if(xtel(ncsw)<=xijtel(i)+small .and. ytel(ncsw)<=yijtel(i)+small .and. &
             xtel(ncse)>=xijtel(i)-small .and. ytel(ncnw)>=yijtel(i)-small )then
            intp(0,j) = 4
            intp(1,j) = ncsw
            intp(2,j) = ncse
            intp(3,j) = ncne
            intp(4,j) = ncnw
            fx = (xijtel(i)-xtel(ncsw))/(xtel(ncse)-xtel(ncsw))
            fy = (yijtel(i)-ytel(ncsw))/(ytel(ncnw)-ytel(ncsw))
            cfx = 1.0 - fx
            cfy = 1.0 - fy                  
            cntp(1,j) = cfx*cfy
            cntp(2,j) = fx*cfy
            cntp(3,j) = fx*fy
            cntp(4,j) = cfx*fy       
            cycle d1
          endif    
        endif
      enddo
        
      !Use triangles for linear interpolation 
      do ii=1,nc   
        distx = abs(xtel(ii)-xijtel(i))
        disty = abs(ytel(ii)-yijtel(i))
        if(2.0*distx<=dxtel(ii)+1.0e-5 .and. 2.0*disty<=dytel(ii)+1.0e-5)then
          xp = xijtel(i)
          yp = yijtel(i)
          xtri(1) = xtel(ii)
          ytri(1) = ytel(ii)
          !Find nearest inclosing triangle
          do k=1,ncf(ii)              
            nck = c2cf(k,ii)
            if(k<ncf(ii))then
              nck2 = c2cf(k+1,ii)
            else
              nck2 = c2cf(1,ii)
            endif  
            if(nck>nc .or. nck2>nc) cycle !Do not interpolate using land cells
            xtri(2) = xtel(nck)
            ytri(2) = ytel(nck)              
            xtri(3) = xtel(nck2)
            ytri(3) = ytel(nck2)        
            if(intriangle(xp,yp,xtri,ytri))then !wave point inside flow triangle
              call interp_coef_tri(xp,yp,xtri,ytri,wcoef) !Calc interp coefficients
              intp(0,j) = 3
              intp(1,j) = ii
              intp(2,j) = nck
              intp(3,j) = nck2
              cntp(1:3,j) = wcoef(1:3)
              cntp(4,j) = 0.0
              cycle d1 !k
            endif
          enddo !k
        endif
      enddo !ii
            
      !Boundaries, find nearest cell
      distmin = 1.0e20 !Initialize
      do ii=1,nc                
        distx = abs(xtel(ii)-xijtel(i))
        disty = abs(ytel(ii)-yijtel(i))  
        dist = sqrt(distx**2+disty**2+small)
        if(dist<distmin)then
          distmin = dist
          intp(1,j) = ii
        endif  
      enddo !ii
      intp(0,j) = 0 !Indicates extrapolation
      cntp(1,j) = xtrapfunc(distmin,xtrapdist)
      cntp(2:4,j) = 0.0  
    enddo d1
    
    return
    end subroutine interp_coef_tel2tel
                  
!**********************************************************************    
    subroutine interp_coef_tel2pts(nc,nD,mf,x0tel,y0tel,azthtel, &
                  xtel,ytel,dxtel,dytel,c2cf,idf,ncf, &            
                  npts,xpts,ypts,xtrapdist,intp,cntp)
! Calculates the interpolation coefficients from a telescoping Cartesian
! grid to another telescoping Cartesian grid
! written by Alex Sanchez, USACE-CHL      
!**********************************************************************
    use const_def, only: deg2rad,small,pi
    use geo_lib, only: intriangle
    use prec_def
    implicit none
    !Parent/Interpolation Telescoping Cartesian grid (Parent)
    integer,       intent(in) :: nc,nD,mf
    integer,       intent(in) :: idf(mf,nD)  !Direction of cell faces
    integer,       intent(in) :: ncf(nD)     !# of cell faces
    integer,       intent(in) :: c2cf(mf,nD) !Connectivity table
    real(ikind),   intent(in) :: x0tel,y0tel,azthtel
    real(ikind),   intent(in) :: xtel(nD),ytel(nD)
    real(ikind),   intent(in) :: dxtel(nD),dytel(nD)
    !Child Points
    integer,       intent(in) :: npts
    real(ikind),   intent(in) :: xpts(npts),ypts(npts)
    integer,    intent(inout) :: intp(0:4,npts)   !Cell Id's to interpolation (parent) grid to interpolate from
    real(ikind),intent(inout) :: cntp(4,npts)     !Interpolation coefficients from parent to child grid
    !Extrapolation distance
    real(ikind),intent(in) :: xtrapdist
    !Internal variables
    integer :: i,ii,k,nck,nck2
    integer :: ncart,ne,nn,nne,ncse,ncne,ncsw,ncnw
    real(ikind) :: dist,distx,disty,distmin
    real(ikind) :: costel,sintel,fx,fy,cfx,cfy
    real(ikind) :: xtri(3),ytri(3),wcoef(3),xp,yp
    real(ikind) :: xijtel(npts),yijtel(npts)    
    real(ikind),parameter :: toldist = 0.001 ![m]
    logical :: iscart(nD)
        
    costel = cos(azthtel*deg2rad)
    sintel = sin(azthtel*deg2rad)
    
    !Coordinates of child (current) telescoping grid
    !with respect to parent telescoping grid
    do i=1,npts
      xijtel(i) = (xpts(i)-x0tel)*costel + (ypts(i)-y0tel)*sintel
      yijtel(i) =-(xpts(i)-x0tel)*sintel + (ypts(i)-y0tel)*costel
    enddo
    
    !open(43,file='Parent.xyz')
    !do i=1,nD !Full grid including land cells
    !  write(43,*) xtel(i),ytel(i),0.0
    !enddo
    !close(43)
    !
    !open(43,file='Child.xyz')
    !do j=1,nintrpcells
    !  i=iintpcells(j)
    !  write(43,*) xijtel(i),yijtel(i),0.0
    !enddo
    !close(43)
    
    !Determine regular cartesian cells of parent telescopting grid
    do i=1,nc
      iscart(i) = .false.   
    enddo
    ncart = 0
    do i=1,nc
      ne = 0
      nn = 0
      do k=1,ncf(i)
        if(idf(k,i)==1)then !North
          nn = nn + 1
        elseif(idf(k,i)==2)then !East
          ne = ne + 1
        endif
      enddo
      if(nn==1 .and. ne==1)then
        nne = 0
        ncse = c2cf(2,i)
        if(ncse==0) cycle
        do k=1,ncf(ncse)     
          if(idf(k,ncse)==1)then !North-East
            nne = nne + 1
          endif
        enddo
        ncne = c2cf(1,ncse)
        if(nne==1 .and. ncne>0)then
          iscart(i) = .true.
          ncart = ncart + 1
        endif
      endif
    enddo
!    write(*,*) 'Number of Regular Cartesian Cells: ',ncart
!    read(*,*)
    
    !Initialize
    intp = 0
    cntp = 0.0
    
d1: do i=1,npts
      !Falls on the same point   
      do ii=1,nc                
        distx = abs(xtel(ii)-xijtel(i))
        disty = abs(ytel(ii)-yijtel(i))  
        dist = sqrt(distx*distx+disty*disty+small)
        if(dist<toldist)then 
          intp(0,i) = 1
          intp(1,i) = ii
          cntp(1,i) = 1.0
          cntp(2:4,i) = 0.0
          cycle d1  
        endif  
      enddo
        
      !Bilinear interpolation
      do ii=1,nc
        if(iscart(ii))then 
          ncsw = ii
          ncnw = c2cf(1,ii)
          ncse = c2cf(2,ii)
          ncne = c2cf(1,ncse)
          if(ncnw>nc .or. ncse>nc .or. ncne>nc) cycle
          if(xtel(ncsw)<=xijtel(i)+small .and. ytel(ncsw)<=yijtel(i)+small .and. &
             xtel(ncse)>=xijtel(i)-small .and. ytel(ncnw)>=yijtel(i)-small )then
            intp(0,i) = 4
            intp(1,i) = ncsw
            intp(2,i) = ncse
            intp(3,i) = ncne
            intp(4,i) = ncnw
            fx = (xijtel(i)-xtel(ncsw))/(xtel(ncse)-xtel(ncsw))
            fy = (yijtel(i)-ytel(ncsw))/(ytel(ncnw)-ytel(ncsw))
            cfx = 1.0-fx
            cfy = 1.0-fy
            cntp(1,i) = cfx*cfy
            cntp(2,i) = fx*cfy
            cntp(3,i) = fx*fy
            cntp(4,i) = cfx*fy
            cycle d1
          endif
        endif
      enddo
        
      !Use triangles for linear interpolation 
      do ii=1,nc   
        distx = abs(xtel(ii)-xijtel(i))
        disty = abs(ytel(ii)-yijtel(i))            
        if(2.0*distx<=dxtel(ii)+1.0e-5 .and. 2.0*disty<=dytel(ii)+1.0e-5)then
          xp = xijtel(i)
          yp = yijtel(i)
          xtri(1) = xtel(ii)
          ytri(1) = ytel(ii)
          !Find nearest inclosing triangle
          do k=1,ncf(ii)
            nck = c2cf(k,ii)
            if(k<ncf(ii))then
              nck2 = c2cf(k+1,ii)
            else
              nck2 = c2cf(1,ii)
            endif  
            if(nck>nc .or. nck2>nc) cycle !Do not interpolate using land cells
            xtri(2) = xtel(nck)
            ytri(2) = ytel(nck)              
            xtri(3) = xtel(nck2)
            ytri(3) = ytel(nck2)        
            if(intriangle(xp,yp,xtri,ytri))then !wave point inside flow triangle
              call interp_coef_tri(xp,yp,xtri,ytri,wcoef) !Calc interp coefficients
              intp(0,i) = 3
              intp(1,i) = ii
              intp(2,i) = nck
              intp(3,i) = nck2
              cntp(1:3,i) = wcoef(1:3)
              cntp(4,i) = 0.0
              cycle d1
            endif
          enddo !k
          exit !ii
        endif           
      enddo !ii
            
      !Boundaries, find nearest cell
      distmin = 1.0e20 !Initialize
      do ii=1,nc                
        distx = abs(xtel(ii)-xijtel(i))
        disty = abs(ytel(ii)-yijtel(i))  
        dist = sqrt(distx**2+disty**2+small)
        if(dist<distmin)then
          distmin = dist
          intp(1,i) = ii
        endif  
      enddo !ii
      intp(0,i) = 0  !Indicates extrapolation
      cntp(1,i) = xtrapfunc(distmin,xtrapdist)
      cntp(2:4,i) = 0.0  
    enddo d1
    
    return
    end subroutine interp_coef_tel2pts

!*******************************************************************************    
    subroutine interp_coef_tri2tel(ne,nn,xn,yn,e2n, &
       nD,x0tel,y0tel,azthtel,xtel,ytel, &
       xtrapdist,nintrpcells,iintpcells,intp,cntp)
! Calculates the interpolation coefficients from an unstructured 
! node-based triangular mesh (such as ADCIRC) to a telescoping Cartesian mesh
! written by Alex Sanchez, USACE-CHL
!*******************************************************************************
    use const_def, only: deg2rad,small,pi
    use geo_lib, only: intriangle,triangle_area
    use prec_def
    implicit none
    !Parent Unstructured Triangular Mesh
    integer,       intent(in) :: ne            !Number of elements
    integer,       intent(in) :: nn            !Number of nodes
    integer,       intent(in) :: e2n(3,ne)     !Element to node connectivity
    real(ikind),   intent(in) :: xn(nn),yn(nn) !Node global coordinates
    !Child Telescoping Cartesian grid
    integer,       intent(in)    :: nD
    real(ikind),   intent(in)    :: x0tel,y0tel,azthtel
    real(ikind),   intent(in)    :: xtel(nD),ytel(nD)
    integer,       intent(in)    :: nintrpcells             !Number of cells to interpolate
    integer,       intent(in)    :: iintpcells(nintrpcells) !Cell Id's of child grid to interpolate
    integer,       intent(inout) :: intp(0:4,nintrpcells)
    real(ikind),   intent(inout) :: cntp(4,nintrpcells)
    !Extrapolation
    real(ikind),   intent(in)  :: xtrapdist
    !Internal variables
    integer :: i,j,k
    real(ikind) :: xtri(3),ytri(3),wcoef(3),xp,yp,costel,sintel
    real(ikind) :: dist,distmin,ae(ne)
    
    costel = cos(azthtel*deg2rad)
    sintel = sin(azthtel*deg2rad)
    
    !Initialize
    intp = 0
    cntp = 0.0
    
    !Precompute element areas
    do k=1,ne !Elements 
      xtri(1:3) = xn(e2n(1:3,k)) !Global coordinates
      ytri(1:3) = yn(e2n(1:3,k)) !Global coordinates
      ae(k) = triangle_area(xtri,ytri)
    enddo
    
d1: do j=1,nintrpcells !Cells
      i = iintpcells(j)
      xp = x0tel + xtel(i)*costel - ytel(i)*sintel  !Global coordinates
      yp = y0tel + xtel(i)*sintel + ytel(i)*costel  !Global coordinates
      
      !Search if point is in any triangles
      do k=1,ne !Elements 
        xtri(1:3) = xn(e2n(1:3,k)) !Global coordinates
        ytri(1:3) = yn(e2n(1:3,k)) !Global coordinates
        if(intriangle(xp,yp,xtri,ytri,ae(k)))then !wave point inside flow triangle
          call interp_coef_tri(xp,yp,xtri,ytri,wcoef) !Calc interp coefficients
          intp(0,j) = 3
          intp(1:3,j) = e2n(1:3,k)
          cntp(1:3,j) = wcoef(1:3)
          cntp(4,j) = 0.0
          cycle d1
        endif
      enddo !j

      !Extrapolation to nearest node
      distmin = 1.0e20
      do k=1,nn
        dist = sqrt((xp-xn(k))**2+(yp-yn(k))**2+small)
        if(dist<distmin)then
          distmin = dist
          intp(1,j) = k
        endif
      enddo
      intp(0,j) = 0 !Indicates extrapolation 
      cntp(1,j) = xtrapfunc(distmin,xtrapdist)
      cntp(2:4,j) = 0.0
    enddo d1
    
    return
    end subroutine interp_coef_tri2tel
    
!**********************************************************************    
    subroutine interp_coef_tri2pts(ne,nn,xn,yn,e2n, &            
                npts,xpts,ypts,xtrapdist,intp,cntp)
! Calculates the interpolation coefficients from an unstructured 
! node-based triangular mesh (such as ADCIRC) to scattered points
! written by Alex Sanchez, USACE-CHL
!**********************************************************************
    use const_def, only: deg2rad,small,pi
    use geo_lib, only: intriangle
    use prec_def
    implicit none
    !Parent Unstructured Triangular Mesh
    integer,       intent(in) :: ne            !# of cells (elements)
    integer,       intent(in) :: nn            !# of nodes
    integer,       intent(in) :: e2n(3,ne)     !Cells (element) to node connectivity
    real(ikind),   intent(in) :: xn(nn),yn(nn) !Node global coordinates
    !Child Points
    integer,       intent(in)    :: npts
    real(ikind),   intent(in)    :: xpts(npts),ypts(npts)
    integer,       intent(inout) :: intp(0:3,npts)
    real(ikind),   intent(inout) :: cntp(3,npts)
    !Extrapolation
    real(ikind),   intent(in)  :: xtrapdist
    !Internal variables
    integer :: i,j,k
    real(ikind) :: xtri(3),ytri(3),wcoef(3),xp,yp
    real(ikind) :: dist,distmin
    
    !Initialize
    intp = 0
    cntp = 0.0 
    
d1: do i=1,npts !Points
      xp = xpts(i) !Global coordinates
      yp = ypts(i) !Global coordinates 
      
      !Search if point is in any triangles
      do j=1,ne !Elements
        xtri(1:3) = xn(e2n(1:3,j)) !Global coordinates
        ytri(1:3) = yn(e2n(1:3,j)) !Global coordinates
        if(intriangle(xp,yp,xtri,ytri))then !wave point inside flow triangle
          call interp_coef_tri(xp,yp,xtri,ytri,wcoef) !Calc interp coefficients
          intp(0,i) = 3
          intp(1:3,i) = e2n(1:3,j)
          cntp(1:3,i) = wcoef(1:3)
          cycle d1
        endif
      enddo !j

      !Extrapolation to nearest node
      distmin = 1.0e20
      do k=1,nn
        dist = sqrt((xp-xn(k))**2+(yp-yn(k))**2+small)
        if(dist<distmin)then
          distmin = dist
          intp(1,i) = k
        endif
      enddo
      intp(0,i) = 0 !Indicates extrapolation
      cntp(1,i) = xtrapfunc(distmin,xtrapdist)
      cntp(2:3,i) = 0.0
    enddo d1
    
    return
    end subroutine interp_coef_tri2pts

!!**********************************************************************
!    subroutine interp_coef_hybridnode2pts(ne,nn,xn,yn,c2n,ncf, &
!                npts,xpts,ypts,xtrapdist,intp,cntp)
!! Calculates the interpolation coefficients from an unstructured 
!! node-based hybrid (triangles and quadrilaterals) mesh to a 
!! scattered points
!! written by Alex Sanchez, USACE-CHL
!!**********************************************************************
!    use const_def, only: deg2rad,small,pi
!    use geo_lib, only: intriangle,inquad,triangle_area,quad_area
!    use prec_def
!    implicit none
!    !Parent Unstructured Hybrid Mesh
!    integer,       intent(in) :: ne            !# of cells (elements)
!    integer,       intent(in) :: nn            !# of nodes
!    integer,       intent(in) :: e2n(4,ne)     !Element to node connectivity
!    integer,       intent(in) :: nef(ne)       !# of faces (nodes) per cell  
!    real(ikind),   intent(in) :: xn(nn),yn(nn) !Node global coordinates
!    !Child Points
!    integer,       intent(in)    :: npts
!    real(ikind),   intent(in)    :: xpts(npts),ypts(npts)
!    integer,       intent(inout) :: intp(0:4,npts)
!    real(ikind),   intent(inout) :: cntp(4,npts)
!    !Extrapolation
!    real(ikind),   intent(in)  :: xtrapdist
!    !Internal variables
!    integer :: i,j,k
!    real(ikind) :: ac(nc)
!    real(ikind) :: xtri(3),ytri(3),wtri(3)
!    real(ikind) :: xquad(4),yquad(4),wquad(3)
!    real(ikind) :: dist,distmin,xp,yp
!    
!    !Initialize
!    intp=0
!    cntp=0.0 
!    
!    !Precompute cell areas and save
!    do j=1,ne !Elements
!      if(nef(j)==3)then !Triangle
!        xtri(1:3) = xn(e2n(1:3,j)) !Global coordinates
!        ytri(1:3) = yn(e2n(1:3,j)) !Global coordinates  
!        ac(j) = triangle_area(xtri,ytri)
!      else     !Quadrilateral
!        xquad(1:4) = xn(e2n(1:4,j)) !Global coordinates
!        yquad(1:4) = yn(e2n(1:4,j)) !Global coordinates  
!        ac(j) = quad_area(xquad,yquad)
!      endif    
!    enddo
!    
!    do i=1,npts !Points
!      xp=xpts(i) !Global coordinates
!      yp=ypts(i) !Global coordinates 
!      
!      !Search if point is in any triangles
!      do j=1,nc !Cells
!        if(ncf(j)==3)then !Triangle
!          xtri(1:3) = xn(c2n(1:3,j)) !Global coordinates
!          ytri(1:3) = yn(c2n(1:3,j)) !Global coordinates
!          if(intriangle(xp,yp,xtri,ytri,ac(j)))then !wave point inside flow triangle
!            call interp_coef_tri(xp,yp,xtri,ytri,wtri) !Calc interp coefficients
!            intp(0,i)=3
!              intp(1:3,i)=c2n(1:3,j)
!            cntp(1:3,i)=wtri(1:3)
!            exit
!          endif
!        else !Quadrilateral
!          xquad(1:4) = xn(c2n(1:4,j)) !Global coordinates
!          yquad(1:4) = yn(c2n(1:4,j)) !Global coordinates
!          if(inquad(xp,yp,xquad,yquad,ac(j)))then !wave point inside flow triangle
!            !call interp_coef_quad(xp,yp,xquad,yquad,wquad) !Calc interp coefficients
!            intp(0,i)=4
!              intp(1:4,i)=c2n(1:4,j)
!            cntp(1:4,i)=wquad(1:4)
!            exit
!          endif  
!        endif
!      enddo !j
!      
!      if(intp(0,i)>0) cycle
!
!      !Extrapolation to nearest node
!      distmin = 1.0e20
!      do k=1,nn
!        dist = sqrt((xp-xn(k))**2+(yp-yn(k))**2+small)
!        if(dist<distmin)then
!          distmin=dist
!          intp(1,i)=k          
!        endif
!      enddo
!      intp(0,i)=0 !Indicates extrapolation
!      cntp(1,i)=xtrapfunc(distmin,xtrapdist)
!      cntp(2:4,i)=0.0          
!    enddo !j-tri
!    
!    return
!    end subroutine interp_coef_hybridnode2pts

!!**********************************************************************    
!    subroutine interp_coef_poly2pts(nc,nD,mf,nn,xn,yn,c2n,ncf, &            
!                npts,xpts,ypts,xtrapdist,intp,cntp)
!! Calculates the interpolation coefficients from an unstructured 
!! hybrid (triangles and quadrilaterals) mesh to a scattered points
!! written by Alex Sanchez, USACE-CHL
!!**********************************************************************
!    use const_def, only: deg2rad,small,pi
!    use geo_lib, only: poly_area,inpoly
!    use prec_def
!    implicit none
!    !Parent Unstructured Hybrid Mesh
!    integer,       intent(in) :: nc            !# of cells
!    integer,       intent(in) :: nD            !# of cells including dummy cells
!    integer,       intent(in) :: mf            !Max # of faces (nodes) per cell  
!    integer,       intent(in) :: nn            !# of nodes
!    integer,       intent(in) :: nnc(nn)       !# of cells per node
!    integer,       intent(in) :: c2n(mf,nc)    !Cell to node connectivity
!    integer,       intent(in) :: ncf(nc)       !# of faces (nodes) per cell  
!    real(ikind),   intent(in) :: xn(nn),yn(nn) !Node global coordinates
!    !Child Points
!    integer,       intent(in)    :: npts
!    real(ikind),   intent(in)    :: xpts(npts),ypts(npts)
!    integer,       intent(inout) :: intp(0:mf,npts)
!    real(ikind),   intent(inout) :: cntp(mf,npts)
!    !Extrapolation
!    real(ikind),   intent(in)  :: xtrapdist
!    !Internal variables
!    integer :: i,j,k
!    real(ikind) :: ac(nc)
!    real(ikind) :: xtri(3),ytri(3)
!    real(ikind) :: xpoly(mf),ypoly(mf),wpoly(mf)
!    real(ikind) :: dist,distmin,xp,yp
!    
!    !Initialize
!    intp=0
!    cntp=0.0 
!    
!    !Precompute cell areas and save
!    do j=1,nc !Cells
!      do k=1,ncf(j)  
!        xpoly(k) = xn(c2n(k,j)) !Global coordinates
!        ypoly(k) = yn(c2n(k,j)) !Global coordinates
!      enddo      
!      ac(j) = poly_area(mf,ncf(i),xpoly,ypoly)
!    enddo
!    
!    do i=1,npts !Points
!      xp=xpts(i) !Global coordinates
!      yp=ypts(i) !Global coordinates 
!      
!      !Search if point is in any triangles
!      do j=1,nn !Cells
!        do k=1,nnc(j)  
!          xpoly(k) = xn(c2n(k,j)) !Global coordinates
!          ypoly(k) = yn(c2n(k,j)) !Global coordinates 
!        enddo
!        if(inpoly(xp,yp,mf,ncf(j),xpoly,ypoly))then !wave point inside flow triangle
!           call interp_coef_poly(xp,yp,mf,ncf(i),xpoly,ypoly,wpoly) !Calc interp coefficients
!           intp(0,i)=ncf(j)
!             intp(1:ncf(j),i)=c2n(1:ncf(j),j)
!          cntp(1:ncf(j),i)=wpoly(1:ncf(j))          
!          exit
!        endif
!      enddo !j
!      
!      if(intp(0,i)>0) cycle
!
!      !Extrapolation to nearest node
!      distmin = 1.0e20
!      do k=1,nn
!        dist = sqrt((xp-xn(k))**2+(yp-yn(k))**2+small)
!        if(dist<distmin)then
!          distmin=dist
!          intp(1,i)=k          
!        endif
!      enddo
!      intp(0,i)=1
!      cntp(1,i)=xtrapfunc(distmin,xtrapdist)
!      cntp(2:4,i)=0.0          
!    enddo !j-tri
!    
!    return
!    end subroutine interp_coef_poly2pts
                
!**********************************************************************    
    subroutine interp_coef_poly2cart(nc,nD,mf,xc,yc,ncf,c2cf, &
    nn,mc,nnc,n2c,ni,nj,x0cart,y0cart,azthcart,xi,yj,xtrapdist,intp,cntp)
! Calculates the interpolation coefficients from an unstructured 
! polygonal mesh to a scattered points
! written by Alex Sanchez, USACE-CHL
!**********************************************************************
    use const_def, only: deg2rad,small,pi
    use geo_def, only: areamin
    use geo_lib, only: poly_area,inpoly,intriangle
    use prec_def
    implicit none
    !Parent Unstructured Hybrid Mesh
    integer,       intent(in) :: nc            !# of cells
    integer,       intent(in) :: nD            !# of cells including dummy cells
    integer,       intent(in) :: mf            !Max # of faces (nodes) per cell
    real(ikind),   intent(in) :: xc(nD),yc(nD) !Cell centroid global coordinates
    integer,       intent(in) :: ncf(nD)       !Number of cell faces
    integer,       intent(in) :: c2cf(mf,nD)
    integer,       intent(in) :: nn            !# of nodes
    integer,       intent(in) :: mc            !Max # of cells per cell node
    integer,       intent(in) :: nnc(nn)       !# of nodes per cell (nnc(:)<=mc)
    integer,       intent(in) :: n2c(mc,nn)    !Cell to node connectivity    
    !Nonuniform Structured Cartesian Grid (Child)
    integer,       intent(in) :: ni,nj
    real(ikind),   intent(in) :: x0cart,y0cart,azthcart
    real(ikind),   intent(in) :: xi(ni),yj(nj)
    integer,       intent(inout) :: intp(0:mc,ni,nj)
    real(ikind),   intent(inout) :: cntp(mc,ni,nj)
    !Extrapolation
    real(ikind),   intent(in)  :: xtrapdist
    !Internal variables
    integer :: i,ii,j,jj,k,ij,nck,nck2,npoly,idpoly(mc)
    integer :: nintpequal,nintptri,nintppoly,nextrap
    !!real(ikind) :: apoly(nn) !Area of cell-connected polygons enclosing nodes
    real(ikind) :: xpoly(mc),ypoly(mc),wpoly(mc)
    real(ikind) :: xp,yp,coscart,sincart
    real(ikind) :: xtri(3),ytri(3),wtri(3)
    real(ikind) :: dist,distmin,tol
    integer :: idintp(ni,nj)
    
    !Initialize
    intp = 0
    cntp = 0.0 
    nintpequal = 0
    nintptri = 0 
    nintppoly = 0
    nextrap = 0
    
    coscart = cos(azthcart*deg2rad)
    sincart = sin(azthcart*deg2rad)
    
    !!!Precompute area of cell-connected polygons enclosing nodes
    !!do jj=1,nn !nodes
    !!  !!if(jj==178)then
    !!  !!  continue  
    !!  !!endif  
    !!  do k=1,nnc(jj)  
    !!    xpoly(k) = xc(n2c(k,jj)) !Global coordinates
    !!    ypoly(k) = yc(n2c(k,jj)) !Global coordinates      
    !!  enddo      
    !!  apoly(jj) = poly_area(mc,nnc(jj),xpoly,ypoly)
    !!enddo
    
    tol = 0.1*sqrt(areamin)
    
    write(*,'(A)') ' Percent   Matching  Interp.    Interp.   '
    write(*,'(A)') ' Complete   Points   Triangles  Polygons  Extrap.'
    
    do i=1,ni
d2:   do j=1,nj
        idintp(i,j) = 0  
        !!if(i==7 .and. j==14)then
        !!  continue
        !!endif  
        ij = j + (i-1)*nj  
        if(mod(ij,5000)==0)then
          write(*,'(F7.2,6I10)') float(ij)*100/float(ni*nj),&
           nintpequal,nintptri,nintppoly,nextrap
        endif  
        xp = x0cart + xi(i)*coscart - yj(j)*sincart !Global coordinates
        yp = y0cart + xi(i)*sincart + yj(j)*coscart !Global coordinates 
        
        !Falls on the same point   
        do ii=1,nc
          if(abs(xc(ii)-xp)<tol .and. abs(yc(ii)-yp)<tol)then 
            intp(0,i,j)=1
            intp(1,i,j)=ii
            cntp(1,i,j)=1.0
            nintpequal = nintpequal + 1
            idintp(i,j) = 1
            cycle d2  
          endif  
        enddo
        
        !Check if in cell-based triangle composed of neighest neighbors
        do ii=1,nc
          xtri(1)=xc(ii)
          ytri(1)=yc(ii)
          do k=1,ncf(ii)
            nck=c2cf(k,ii)
            if(k<ncf(ii))then
              nck2=c2cf(k+1,ii)
            else
              nck2=c2cf(1,ii)
            endif
            if(nck>nc .or. nck2>nc) cycle !Do not interpolate using land cells
            xtri(2)=xc(nck); xtri(3)=xc(nck2)
            ytri(2)=yc(nck); ytri(3)=yc(nck2)  
            if(intriangle(xp,yp,xtri,ytri))then !wave point inside flow triangle
              call interp_coef_tri(xp,yp,xtri,ytri,wtri) !Calc interp coefficients
              intp(0,i,j)=3
              intp(1,i,j)=ii
              intp(2,i,j)=nck
              intp(3,i,j)=nck2
              cntp(1:3,i,j)=wtri(1:3)
              nintptri = nintptri + 1
              idintp(i,j) = 2
              cycle d2
            endif !intriangle
          enddo !k
        enddo
        
        !Check if point is in cell-based polygons surrounding nodes
        do jj=1,nn !Nodes
          !!if(jj==266)then
          !!  continue
          !!endif
          xpoly = 0.0; ypoly = 0.0 !Initialize
          npoly = 0
          do k=1,nnc(jj)
            nck=n2c(jj,k)
            if(nck>nc) cycle !Do not interpolate using ghost cells
            npoly = npoly + 1
            idpoly(k) = nck
            xpoly(k) = xc(nck) !Global coordinates
            ypoly(k) = yc(nck) !Global coordinates      
          enddo !k
          if(inpoly(xp,yp,mc,npoly,xpoly,ypoly))then !wave point inside flow triangle
            !!if(jj==178)then
            !!  continue  
            !!endif   
            call interp_coef_poly(xp,yp,mc,npoly,xpoly,ypoly,wpoly) !Calc interp coefficients
            intp(0,i,j)=npoly !# of points
            intp(1:npoly,i,j)=idpoly(1:npoly) !point IDs
            cntp(1:npoly,i,j)=wpoly(1:npoly)
            nintppoly = nintppoly + 1
            idintp(i,j) = 3
            cycle d2
          endif
        enddo

        !Extrapolation to nearest cell
        distmin = 1.0e20
        do k=1,nc
          dist = sqrt((xp-xc(k))**2+(yp-yc(k))**2+small)
          if(dist<distmin)then
            distmin=dist
            intp(1,i,j)=k          
          endif
        enddo !k
        intp(0,i,j)=1
        cntp(1,i,j)=xtrapfunc(distmin,xtrapdist)
        cntp(2:mc,i,j)=0.0
        nextrap = nextrap + 1
        idintp(i,j) = 4
      enddo d2 !j   
    enddo !i
    
    write(*,'(F7.2,6I10)') 100.0,&
           nintpequal,nintptri,nintppoly,nextrap
    
    open(224,file='Poly2cart.xyz')
    do i=1,ni
      do j=1,nj
        xp=x0cart+xi(i)*coscart-yj(j)*sincart !Global coordinates
        yp=y0cart+xi(i)*sincart+yj(j)*coscart !Global coordinates   
        write(224,*) xp,yp,idintp(i,j)
      enddo
    enddo
    close(224) 
    
    return
    end subroutine interp_coef_poly2cart           
                
!*************************************************************   
    subroutine interp_coef_curv2tel(ni,nj,xcurv,ycurv, &            
       nD,x0tel,y0tel,azthtel,xtel,ytel,&
       nintrpcells,iintpcells,ijntp,cntp)
! Interpolation coefficients and index from a curvilinear grid 
! to a telescoping grid.
! written by Chris Reed, URS; Alex Sanchez, USACE-CHL    
!*************************************************************
    use const_def, only: deg2rad
    use diag_lib
    use prec_def
    implicit none
    !Input/Output
    !Parent Curvilinear grid
    integer,       intent(in)    :: ni,nj
    real(ikind),   intent(in)    :: xcurv(ni,nj),ycurv(ni,nj)
    !Child Telescoping Cartesian grid
    integer,       intent(in)    :: nD
    real(ikind),   intent(in)    :: x0tel,y0tel,azthtel
    real(ikind),   intent(in)    :: xtel(nD),ytel(nD)
    integer,       intent(in)    :: nintrpcells             !Number of cells to interpolate
    integer,       intent(in)    :: iintpcells(nintrpcells) !Cell Id's of child grid to interpolate    
    integer,       intent(inout) :: ijntp(2,nD)
    real(ikind),   intent(inout) :: cntp(2,2,nD)    
    !Internal    
    integer :: i,j,k,ii,ii1,jj1,iimax,jjmax
    real(ikind) :: Ax,Ay,Bx,By,Cx,Cy,fxip,fyip
    real(ikind) :: Area_1,Area_2,Area_3,Area_4,tot_area
    real(ikind) :: sarea(ni,nj),costel,sintel
    character(len=200) :: msg
    
    costel = cos(azthtel*deg2rad)
    sintel = sin(azthtel*deg2rad)
    
    !calculate area of each Quadrilateral grid cell
    do i=1,ni-1
      do j=1,nj-1
        !Subdivide into two triangles
        Ax = xcurv(i,j);     Ay = ycurv(i,j)
        Bx = xcurv(i+1,j);   By = ycurv(i+1,j)
        Cx = xcurv(i+1,j+1); Cy = ycurv(i+1,j+1)
        Area_1 = abs(Ax*(By-Cy)+Bx*(Cy-Ay)+Cx*(Ay-By))/2.0
        Ax = xcurv(i+1,j);   Ay = ycurv(i+1,j)
        Bx = xcurv(i+1,j+1); By = ycurv(i+1,j+1)
        Cx = xcurv(i,j);     Cy = ycurv(i,j)
        Area_2 = abs(Ax*(By-Cy)+Bx*(Cy-Ay)+Cx*(Ay-By))/2.0
        sarea(i,j) = Area_1 + Area_2
      enddo
    enddo     
            
    !interpolate Telescoping Grid onto Curvilinear grid    
    ii1 = 1; iimax = ni-1
    jj1 = 1; jjmax = nj-1  
    do k=1,nintrpcells
      ii = iintpcells(k)
      Cx = x0tel + xtel(ii)*costel - ytel(ii)*sintel  !Global coordinates
      Cy = y0tel + xtel(ii)*sintel + ytel(ii)*costel  !Global coordinates
5556  continue       
      do i=1,ni-1
        do j=1,nj-1
          !calculate four areas - sub-dividing wind grid quadilateral
          Ax = xcurv(i,j);   Ay = ycurv(i,j)
          Bx = xcurv(i+1,j); By = ycurv(i+1,j)
          Area_1 = abs(Ax*(By-Cy)+Bx*(Cy-Ay)+Cx*(Ay-By))/2.0
          Ax = xcurv(i+1,j);   Ay = ycurv(i+1,j)
          Bx = xcurv(i+1,j+1); By = ycurv(i+1,j+1)
          Area_2 = abs(Ax*(By-Cy)+Bx*(Cy-Ay)+Cx*(Ay-By))/2.0 
          Ax = xcurv(i+1,j+1); Ay = ycurv(i+1,j+1)
          Bx = xcurv(i,j+1);   By = ycurv(i,j+1)
          Area_3 = abs(Ax*(By-Cy)+Bx*(Cy-Ay)+Cx*(Ay-By))/2.0
          Ax = xcurv(i,j+1);   Ay = ycurv(i,j+1)
          Bx = xcurv(i,j);     By = ycurv(i,j)
          Area_4 = abs(Ax*(By-Cy)+Bx*(Cy-Ay)+Cx*(Ay-By))/2.0
          tot_area  = area_1 + area_2 + area_3 + area_4
          !if(tot_area < 1.02*sarea(i,j))then !we have found the proper WND cell
          if(abs(tot_area-sarea(i,j))<0.001*sarea(i,j))then !we have found the proper WND cell
            ijntp(1,k) = i
            ijntp(2,k) = j
            fxip = area_1/(area_1+area_3)
            fyip = area_4/(area_2+area_4)
            cntp(1,1,k) = (1.0-fxip)*(1.0-fyip)
            cntp(2,1,k) = (1.0-fxip)*fyip
            cntp(1,2,k) = fxip*(1.0-fyip)            
            cntp(2,2,k) = fxip*fyip
            !set search local to current WND cell
            ii1 = max(i-1,1)
            iimax = min(i+1,ni-1)
            jj1 = max(j-1,1)
            jjmax = min(j+1,nj-1)
!            write(dgunit,*)'found',n,i,j,' srch = ',ii1,iimax,jj1,jjmax
!            write(dgunit,*)n,i,j,wnd_spat(n).F1,wnd_spat(n).F2
            goto 5555  !move to next CMS grid cell
          endif  
        enddo !j
      enddo !i
      
      !we are here because no WND cell has been located - increase search domain
      !but first see if we have expanded search to full grid already
      if(ii1==1 .and. iimax==ni-1 .and. jj1==1 .and. jjmax==nj-1)then
        write(msg,*) '  No wind grid cell found for CMS cell i,x,y = ',k,Cx,Cy
        call diag_print_message(msg)
        !!Nearest neighbor
        !distmin = 1.0e20
        !do i=1,ni
        !  do j=1,nj
        !    
        !  enddo
        !enddo
        goto 5555  !move to next CMS grid cell
      else !expand search area and continue
        ii1 = max(ii1-1,1)
        iimax = min(iimax+1,ni-1)
        jj1 = max(jj1-1,1)
        jjmax = min(jjmax+1,nj-1) 
        write(msg,*) '  Expanding search: ',ii1,iimax,jj1,jjmax
        call diag_print_message(msg)
        goto 5556  !start seach over - with larger domain      
      endif
5555  continue    
    enddo  !end of all CMS cells
    
    return
    end subroutine interp_coef_curv2tel

!*******************************************************************************   
    subroutine interp_coef_curv2pts(ni,nj,xcurv,ycurv,np,nintp,xp,yp,ijntp,cntp)
! Interpolates from a curvilinear grid to scattered points
! written by Chris Reed, URS; Alex Sanchez, USACE-CHL    
!*******************************************************************************
#include "CMS_cpp.h"
    use const_def, only: deg2rad
    use diag_lib
    use prec_def
    implicit none
    !Input/Output
    !Parent Curvilinear grid
    integer,       intent(in)    :: ni,nj
    real(ikind),   intent(in)    :: xcurv(ni,nj),ycurv(ni,nj)
    !Child Interpolation Points
    integer,       intent(in)    :: np,nintp
    real(ikind),   intent(in)    :: xp(np),yp(np)
    integer,       intent(inout) :: ijntp(2,np)
    real(ikind),   intent(inout) :: cntp(2,2,np)    
    !Internal Variables    
    integer :: i,j,k,ii1,jj1,iimax,jjmax
    real(ikind) :: Ax,Ay,Bx,By,Cx,Cy,fxip,fyip,distmin,dist
    real(ikind) :: Area_1,Area_2,Area_3,Area_4,tot_area
    real(ikind) :: sarea(ni,nj)
    character(len=200) :: msg
    
    !calculate area of each quadrilateral grid cell
    do i=1,ni-1
      do j=1,nj-1
        !Subdivide into two triangles
        Ax = xcurv(i,j);     Ay = ycurv(i,j)
        Bx = xcurv(i+1,j);   By = ycurv(i+1,j)
        Cx = xcurv(i+1,j+1); Cy = ycurv(i+1,j+1)
        Area_1 = abs(Ax*(By-Cy)+Bx*(Cy-Ay)+Cx*(Ay-By))/2.0
        Ax = xcurv(i+1,j);   Ay = ycurv(i+1,j)
        Bx = xcurv(i+1,j+1); By = ycurv(i+1,j+1)
        Cx = xcurv(i,j);     Cy = ycurv(i,j)
        Area_2 = abs(Ax*(By-Cy)+Bx*(Cy-Ay)+Cx*(Ay-By))/2.0
        sarea(i,j) = Area_1 + Area_2
      enddo
    enddo     
            
    !interpolate Telescoping Grid onto Curvilinear grid    
    ii1 = 1; iimax = ni-1
    jj1 = 1; jjmax = nj-1  
    do k=1,nintp
      Cx=xp(k); Cy=yp(k) !Global coordinates
5556  continue       
      do i=1,ni-1
        do j=1,nj-1
          !calculate four areas - sub-dividing wind grid quadilateral
          Ax = xcurv(i,j);     Ay = ycurv(i,j)
          Bx = xcurv(i+1,j);   By = ycurv(i+1,j)                    
          Area_1 = abs(Ax*(By-Cy)+Bx*(Cy-Ay)+Cx*(Ay-By))/2.0
          Ax = xcurv(i+1,j);   Ay = ycurv(i+1,j)
          Bx = xcurv(i+1,j+1); By = ycurv(i+1,j+1)
          Area_2 = abs(Ax*(By-Cy)+Bx*(Cy-Ay)+Cx*(Ay-By))/2.0 
          Ax = xcurv(i+1,j+1); Ay = ycurv(i+1,j+1)
          Bx = xcurv(i,j+1);   By = ycurv(i,j+1)
          Area_3 = abs(Ax*(By-Cy)+Bx*(Cy-Ay)+Cx*(Ay-By))/2.0
          Ax = xcurv(i,j+1);   Ay = ycurv(i,j+1)
          Bx = xcurv(i,j);     By = ycurv(i,j)
          Area_4 = abs(Ax*(By-Cy)+Bx*(Cy-Ay)+Cx*(Ay-By))/2.0
          tot_area  = area_1 + area_2 + area_3 + area_4
          if(tot_area < 1.001*sarea(i,j))then !we have found the proper WND cell
            ijntp(1,k) = i
            ijntp(2,k) = j
            fxip = area_1/(area_1+area_3)
            fyip = area_4/(area_2+area_4)
            cntp(1,1,k) = (1.0-fxip)*(1.0-fyip)
            cntp(2,1,k) = (1.0-fxip)*fyip
            cntp(1,2,k) = fxip*(1.0-fyip)            
            cntp(2,2,k) = fxip*fyip
            !set search local to current WND cell
            ii1 = max(i-1,1)
            iimax = min(i+1,ni-1)
            jj1 = max(j-1,1)
            jjmax = min(j+1,nj-1)
!            write(dgunit,*)'found',n,i,j,' srch = ',ii1,iimax,jj1,jjmax
!            write(dgunit,*)n,i,j,wnd_spat(n).F1,wnd_spat(n).F2
            goto 5555  !move to next CMS grid cell
          endif  
        enddo !j
      enddo !i
      
      !we are here because no WND cell has been located - increase search domain
      !but first see if we have expanded search to full grid already
      if(ii1==1 .and. iimax==ni-1 .and. jj1==1 .and. jjmax==nj-1)then
        write(msg,*) '  Cell=',k,', x=',Cx,'m, y=',Cy,'m'
        call diag_print_warning('No wind grid cell found for: ',msg,'  Using nearest neighbor') 
        !Use neighest neighbor
        distmin=1.0e20
        do i=1,ni-1
          do j=1,nj-1
            dist = sqrt((xp(k)-xcurv(i,j))**2+(yp(k)-ycurv(i,j))**2+1.0e-20)
            if(dist<distmin)then
              distmin = dist  
              ijntp(1,k) = i
              ijntp(2,k) = j
              cntp(1,1,k) = 1.0            
            endif
          enddo !j
        enddo !i
        goto 5555  !move to next CMS grid cell
      else !expand search area and continue
        ii1 = max(ii1-1,1)
        iimax = min(iimax+1,ni-1)
        jj1 = max(jj1-1,1)
        jjmax = min(jjmax+1,nj-1) 
#ifdef DEV_MODE
        write(msg,*) '  Expanding search:',ii1,iimax,jj1,jjmax
        call diag_print_message(msg)
#endif
        goto 5556  !start seach over - with larger domain      
      endif
5555  continue    
    enddo  !end of all CMS cells
    
    return
    end subroutine interp_coef_curv2pts

!**********************************************************************  
    subroutine interp_scal_curv2tel(ni,nj,nD,nc,ijntp,cntp,&
                 undef,low,high,defval,mat1,vec1)
! Spatial interpolation from curvilinear to telescoping grids
! Uses only the values which are different from the undefined value 
! and specified bounds. If not valid points are found than the default
! value is assigned.
! written by Alex Sanchez, USACE-CHL    
!**********************************************************************
    use prec_def
    implicit none
    !Input/Output
    integer,    intent(in) :: ni,nj            !Curvilinear grid dimensions
    real(ikind),intent(in) :: mat1(ni,nj)      !Curvilinear grid matrix
    integer,    intent(in) :: nD,nc            !Telescoping grid dimensions
    integer,    intent(in) :: ijntp(2,nD)      !Interpolation index
    real(ikind),intent(in) :: cntp(2,2,nD)     !Interpolation coefficients
    real(ikind),intent(in) :: undef,low,high,defval !Used to determine valid data points
    real(ikind),intent(out) :: vec1(nD)        !Telescoping grid vector
    !Internal variables
    integer :: ii,i,j,in,jm,iin1,jjm1
    real(ikind) :: csum
    
!$OMP PARALLEL DO PRIVATE(ii,i,j,in,jm,iin1,jjm1,csum)      
      do ii=1,nc
        i=ijntp(1,ii)
        j=ijntp(2,ii)      
        if(i*j==0)then
          vec1(ii)=defval
          cycle    
        endif
        !Use only values within valid range 
        vec1(ii)=0.0
        csum=0.0
        do in=1,2
          do jm=1,2
            iin1=i+in-1
            jjm1=j+jm-1  
            if(abs(mat1(iin1,jjm1)-undef)>1.0e-3 .and. &
               mat1(iin1,jjm1)>low .and. mat1(iin1,jjm1)<high)then
              vec1(ii)  = vec1(ii)  + cntp(in,jm,ii)*mat1(iin1,jjm1)
              csum = csum + cntp(in,jm,ii)
            endif
          enddo
        enddo
        if(csum>1.0e-6)then
          vec1(ii) = vec1(ii)/csum
        else
          vec1(ii) = defval
        endif
        !if(abs(mat1(i+in,j+jm)-undef)<1.0e-3)then
        !  write(*,*) 'ERROR: No defined neighboring points found'
        !  write(*,*) 'i=',i,', j=',j
        !  stop               
        !endif
      enddo           
!$OMP END PARALLEL DO
    
    return
    end subroutine interp_scal_curv2tel    
    
!**********************************************************************  
    subroutine interp_vec_curv2tel(ni,nj,nD,nc,ijntp,cntp,&
                 undef,low,high,defval,mat1,vec1,mat2,vec2)
! Spatial interpolation from curvilinear to telescoping grids
! Uses only the values which are different from the undefined value 
! and specified bounds. If not valid points are found than the default
! value is assigned.
! written by Alex Sanchez, USACE-CHL    
!**********************************************************************
    use prec_def
    implicit none
    !Input/Output
    integer,    intent(in) :: ni,nj                   !Curvilinear grid dimensions
    real(ikind),intent(in) :: mat1(ni,nj),mat2(ni,nj) !Curvilinear grid matrix
    integer,    intent(in) :: nD,nc                   !Telescoping grid dimensions
    integer,    intent(in) :: ijntp(2,nD)             !Interpolation index
    real(ikind),intent(in) :: cntp(2,2,nD)            !Interpolation coefficients
    real(ikind),intent(in) :: undef,low,high,defval   !Used to determine valid data points
    real(ikind),intent(out) :: vec1(nD),vec2(nD)      !Telescoping grid vector
    !Internal variables
    integer :: ii,i,j,in,jm,iin1,jjm1
    real(ikind) :: csum
    
!$OMP PARALLEL DO PRIVATE(ii,i,j,in,jm,iin1,jjm1,csum)      
      do ii=1,nc
        i=ijntp(1,ii)
        j=ijntp(2,ii)      
        if(i*j==0)then
          vec1(ii)=defval
          vec2(ii)=defval
          cycle    
        endif
        !Use only values within valid range 
        vec1(ii)=0.0
        vec2(ii)=0.0
        csum=0.0
        do in=1,2
          do jm=1,2
            iin1=i+in-1
            jjm1=j+jm-1
            if(abs(mat1(iin1,jjm1)-undef)>1.0e-3 .and. &
               abs(mat2(iin1,jjm1)-undef)>1.0e-3 .and. &
               mat1(iin1,jjm1)>low .and. mat1(iin1,jjm1)<high .and. &
               mat2(iin1,jjm1)>low .and. mat2(iin1,jjm1)<high)then
              vec1(ii) = vec1(ii) + cntp(in,jm,ii)*mat1(iin1,jjm1)
              vec2(ii) = vec2(ii) + cntp(in,jm,ii)*mat2(iin1,jjm1)
              csum = csum + cntp(in,jm,ii)
            endif
          enddo
        enddo
        if(csum>1.0e-6)then
          vec1(ii) = vec1(ii)/csum
          vec2(ii) = vec2(ii)/csum
        else
          vec1(ii) = defval
          vec2(ii) = defval
        endif
        !if(abs(mat1(i+in,j+jm)-undef)<1.0e-3 .or. &
        !       abs(mat2(i+in,j+jm)-undef)<1.0e-3)then
        !  write(*,*) 'ERROR: No defined neighboring points found'
        !  write(*,*) 'i=',i,', j=',j
        !  read(*,*)
        !  stop               
        !endif               
      enddo           
!$OMP END PARALLEL DO
    
    return
    end subroutine interp_vec_curv2tel
                 
!**********************************************************************  
    subroutine interp_curv2pts(ni,nj,npts,ijntp,cntp,&
                 undef,low,high,defval,mat1,vec1,mat2,vec2)
! Spatial interpolation from curvilinear to telescoping grids
! Uses only the values which are different from the undefined value 
! and specified bounds. If not valid points are found than the default
! value is assigned.
! written by Alex Sanchez, USACE-CHL    
!**********************************************************************
    use prec_def
    implicit none
    !Input/Output
    integer,    intent(in) :: ni,nj            !Curvilinear grid dimensions
    real(ikind),intent(in) :: mat1(ni,nj)      !Curvilinear grid matrix
    real(ikind),intent(in),optional :: mat2(ni,nj) !Curvilinear grid matrix
    integer,    intent(in) :: npts              !Points dimensions
    integer,    intent(in) :: ijntp(2,npts)      !Interpolation index
    real(ikind),intent(in) :: cntp(2,2,npts)     !Interpolation coefficients
    real(ikind),intent(in) :: undef,low,high,defval !Used to determine valid data points
    real(ikind),intent(out) :: vec1(npts)        !Telescoping grid vector
    real(ikind),intent(out),optional :: vec2(npts) !Telescoping grid vector
    !Internal variables
    integer :: ii,i,j,in,jm,iin1,jjm1
    real(ikind) :: csum
    
    if(present(mat2))then
      !$OMP PARALLEL DO PRIVATE(ii,i,j,in,jm,iin1,jjm1,csum)      
      do ii=1,npts
        i=ijntp(1,ii)
        j=ijntp(2,ii)      
        if(i*j==0)then
          vec1(ii)=defval
          vec2(ii)=defval
          cycle    
        endif
        !Use only values within valid range 
        vec1(ii)=0.0; vec2(ii)=0.0
        csum=0.0
        do in=1,2
          do jm=1,2
            iin1=i+in-1
            jjm1=j+jm-1  
            if(abs(mat1(iin1,jjm1)-undef)>1.0e-3 .and. &
               abs(mat2(iin1,jjm1)-undef)>1.0e-3 .and. &
               mat1(iin1,jjm1)>low .and. mat1(i+in,jjm1)<high .and. &
               mat2(iin1,jjm1)>low .and. mat2(iin1,jjm1)<high)then
              vec1(ii)  = vec1(ii)  + cntp(in,jm,ii)*mat1(iin1,jjm1)
              vec2(ii)  = vec2(ii)  + cntp(in,jm,ii)*mat2(iin1,jjm1)
              csum = csum + cntp(in,jm,ii)
            endif
          enddo
        enddo
        if(csum>1.0e-6)then
          vec1(ii) = vec1(ii)/csum
          vec2(ii) = vec2(ii)/csum
        else
          vec1(ii) = defval
          vec2(ii) = defval
        endif
        !if(abs(mat1(i+in,j+jm)-undef)<1.0e-3 .or. &
        !       abs(mat2(i+in,j+jm)-undef)<1.0e-3)then
        !  write(*,*) 'ERROR: No defined neighboring points found'
        !  write(*,*) 'i=',i,', j=',j
        !  read(*,*)
        !  stop               
        !endif               
      enddo           
      !$OMP END PARALLEL DO
    else
      !$OMP PARALLEL DO PRIVATE(ii,i,j,in,jm,iin1,jjm1,csum)      
      do ii=1,npts
        i=ijntp(1,ii)
        j=ijntp(2,ii)      
        if(i*j==0)then
          vec1(ii)=defval
          cycle    
        endif
        !Use only values within valid range 
        vec1(ii)=0.0
        csum=0.0
        do in=1,2
          do jm=1,2
            iin1=i+in-1
            jjm1=j+jm-1    
            if(abs(mat1(iin1,jjm1)-undef)>1.0e-3 .and. &
               mat1(iin1,jjm1)>low .and. mat1(iin1,jjm1)<high)then
              vec1(ii)  = vec1(ii)  + cntp(in,jm,ii)*mat1(iin1,jjm1)
              csum = csum + cntp(in,jm,ii)
            endif
          enddo
        enddo
        if(csum>1.0e-6)then
          vec1(ii) = vec1(ii)/csum
        else
          vec1(ii) = defval
        endif
        !if(abs(mat1(i+in,j+jm)-undef)<1.0e-3)then
        !  write(*,*) 'ERROR: No defined neighboring points found'
        !  write(*,*) 'i=',i,', j=',j
        !  stop               
        !endif
      enddo           
      !$OMP END PARALLEL DO
    endif
    
    return
    end subroutine interp_curv2pts                
       
!************************************************************    
    function xtrapfunc(dist,xtrapdist)
! Extrapolation function
!************************************************************    
    use const_def, only: pi
    use prec_def
    implicit none    
    real(ikind),intent(in) :: dist,xtrapdist
    real(ikind) :: xtrapfunc
    
    xtrapfunc = 0.5+0.5*cos(pi*min(dist,xtrapdist)/xtrapdist)
    
    return
    end function xtrapfunc
    
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!  End grid-to-grid interpolations
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!  Begin cell-to-face interpolations
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

!****************************************************************************
    subroutine interp_coef_cell2face
! Cell-to-face interpolation coefficient/function
! written by Alex Sanchez, USACE-CHL
!****************************************************************************
    use size_def
    use geo_def
    use interp_def
    use prec_def
    implicit none
    integer :: i,ii,j,k,nck,jcn
    real(ikind) :: dof,dif
    
    allocate(fintp(nmaxfaces,ncellsD))
    fintp=0.0
    
    !Regular Cartesian cells
    !No reconstruction needed
    do ii=1,ncellsimple
      i=idcellsimple(ii)
      fintp(1,i)=dy(i)/(dy(i)+dy(cell2cell(1,i)))
      fintp(2,i)=dx(i)/(dx(i)+dx(cell2cell(2,i)))
      fintp(3,i)=dy(i)/(dy(i)+dy(cell2cell(3,i)))
      fintp(4,i)=dx(i)/(dx(i)+dx(cell2cell(4,i)))
    enddo
    
    !Joint Cartesian cells
    !Cell-reconstruction
    do ii=1,ncelljoint
      i=idcelljoint(ii)  
      do k=1,ncface(i)
        nck=cell2cell(k,i)  
        if(idirface(k,i)==2 .or. idirface(k,i)==4) then     !east/west faces
           fintp(k,i)=dx(i)/(dx(i)+dx(nck))
         else    !north/south faces
           fintp(k,i)=dy(i)/(dy(i)+dy(nck))
         endif
      enddo
    enddo
    
    !Polygonal cells
    !Cell-reconstruction
    do i=1,ncellpoly
      do j=1,nxyface(i)
        k=kxyface(j,i)  
      !do k=1,ncface(i)
        nck=cell2cell(k,i)  
        jcn=llec2llec(k,i)
        dof=sqrt(rnx(k,i)**2+rny(k,i)**2)
        dif=sqrt(rnx(jcn,nck)**2+rny(jcn,nck)**2)
        fintp(k,i)=dof/(dof+dif)
        fintp(jcn,nck)=1.0-fintp(k,i)
      enddo
    enddo
    
    return
    end subroutine interp_coef_cell2face

!********************************************************************
    subroutine interp_scal_cell2face(phi,ibc,phik,dphix,dphiy)
! Scalar interpolation from cell centers to cell faces
! written by Alex Sanchez, USACE-CHL
!********************************************************************
    use size_def
    use geo_def
    use flow_def, only: iwet
    use comvarbl, only: skewcor
    !use der_def
    use interp_def
    use prec_def
    implicit none
    !Input/Output
    real(ikind),intent(in) :: phi(ncellsD) !Input scalar variable
    integer,    intent(in) :: ibc !Boundary condition specification, 0-None,1-Dry-zero-gradient,2-Bnd-zero-gradient
    real(ikind),intent(out) :: phik(nmaxfaces,ncellsD)
    real(ikind),intent(inout),optional :: dphix(ncellsD),dphiy(ncellsD)
    !Internal variables
    integer :: i,ii,j,k,nck,jcn
    real(ikind) :: phitemp, phiktemp
    !logical :: isnankind
    
!--- First-Second Order Interpolation ----------------------
    select case(ibc)        
    case(-1) !Modified Zero-gradient at wet/dry faces
!$OMP PARALLEL DO PRIVATE(i,j,k,nck,jcn)
      do i=1,ncells
        do j=1,nxyface(i)
          k=kxyface(j,i)
          nck=cell2cell(k,i) !Forward connectivity
          jcn=llec2llec(k,i) !Backward connectivity
          if(iwet(i)*iwet(nck)==1)then !Wet-wet boundary
            phik(k,i)=fintp(k,i)*(phi(nck)+rx(jcn,nck)*dphix(nck)+ry(jcn,nck)*dphiy(nck)) &
                     +(1.0-fintp(k,i))*(phi(i)+rx(k,i)*dphix(i)+ry(k,i)*dphiy(i))
            !     +(1.0-fintp(k,i))*(rx(k,i)*dphix(i)+ry(k,i)*dphiy(i)) !Second order cell-reconstruction
            phik(k,i)=min(max(phik(k,i),min(phi(i),phi(nck))),max(phi(i),phi(nck))) !Enforce Local Monotonicity
          elseif(iwet(i)==1)then  !Wet-dry
            phik(k,i)=phi(i)+rx(k,i)*dphix(i)+ry(k,i)*dphiy(i)
          elseif(iwet(nck)==1)then  !Dry-wet
            phik(k,i)=phi(nck)+rx(jcn,nck)*dphix(nck)+ry(jcn,nck)*dphiy(nck)
          else !Dry-Dry
            phik(k,i)=0.5*(phi(i)+phi(nck))  
          endif
          phik(jcn,nck)=phik(k,i)
        enddo !j
      enddo !i   
!$OMP END PARALLEL DO

    case(0) !No treatment
!$OMP PARALLEL DO PRIVATE(i,j,k,nck,jcn)
      do i=1,ncells
        do j=1,nxyface(i)
          k=kxyface(j,i)
          nck=cell2cell(k,i) !Forward connectivity
          jcn=llec2llec(k,i) !Backward connectivity
          phitemp=phi(nck)

          phiktemp = fintp(k,i)*phi(nck)+(1.0-fintp(k,i))*phi(i) !First order for skewed cells
          phik(k,i)= phiktemp

          phiktemp = phik(k,i)
          phik(jcn,nck) = phiktemp
        enddo !j
      enddo !i   
!$OMP END PARALLEL DO

    case(1) !Zero-gradient at wet/dry faces
!!$OMP PARALLEL DO PRIVATE(i,j,k,nck,jcn)
      do i=1,ncells
        do j=1,nxyface(i)
          k=kxyface(j,i)
          nck=cell2cell(k,i) !Forward connectivity
          jcn=llec2llec(k,i) !Backward connectivity
          if(iwet(i)*iwet(nck)==1)then !Wet-wet boundary
            phik(k,i)=fintp(k,i)*phi(nck)+(1.0-fintp(k,i))*phi(i) !First order for skewed cells
          elseif(iwet(i)==1)then  !Wet-dry
            phik(k,i)=phi(i)
            !phik(k,i)=phi(i)+rx(k,i)*dphix(i)+ry(k,i)*dphiy(i)
          elseif(iwet(nck)==1)then  !Dry-wet
            phik(k,i)=phi(nck)
            !phik(k,i)=phi(nck)+rx(jcn,nck)*dphix(nck)+ry(jcn,nck)*dphiy(nck)
          else !Dry-Dry
            phik(k,i)=0.5*(phi(i)+phi(nck))  
          endif
          phik(jcn,nck)=phik(k,i)
        enddo !j
      enddo !i   
!!$OMP END PARALLEL DO
    
    case(2) !Zero-gradient at dry boundary faces
!$OMP PARALLEL DO PRIVATE(i,j,k,nck,jcn)
      do i=1,ncells
        do j=1,nxyface(i)
          k=kxyface(j,i)
          nck=cell2cell(k,i) !Forward connectivity
          jcn=llec2llec(k,i) !Backward connectivity
          if(iwet(nck)==0 .and. nck>ncells)then !dry boundary
            phik(k,i)=phi(i)    
          else            
            phik(k,i)=fintp(k,i)*phi(nck)+(1.0-fintp(k,i))*phi(i) !First order for skewed cells
          endif
          phik(jcn,nck)=phik(k,i)
        enddo !j
      enddo !i   
!$OMP END PARALLEL DO

    case(3) !Zero at dry faces
!$OMP PARALLEL DO PRIVATE(i,j,k,nck,jcn)
      do i=1,ncells
        do j=1,nxyface(i)
          k=kxyface(j,i)
          nck=cell2cell(k,i) !Forward connectivity
          jcn=llec2llec(k,i) !Backward connectivity
          if(iwet(i)*iwet(nck)==1)then !dry boundary
            phik(k,i)=fintp(k,i)*phi(nck)+(1.0-fintp(k,i))*phi(i) !First order for skewed cells  
          else
            phik(k,i)=0.0
          endif
          phik(jcn,nck)=phik(k,i)
        enddo !j
      enddo !i   
!$OMP END PARALLEL DO   
        
    end select

!--- Skewness correction ----------------------------------------
    if(skewcor .and. present(dphiy))then 
!$OMP PARALLEL DO PRIVATE(i,ii,j,k,nck,jcn)  
      do ii=1,ncelljoint
        i=idcelljoint(ii)  
        !if(iwet(i)==0) cycle
        !if(maxval(cell2cell(1:ncface(i),i))>ncells) cycle
        !if(minval(iwet(cell2cell(1:ncface(i),i)))==0) cycle
        do j=1,nxyface(i)
          k=kxyface(j,i)
          nck=cell2cell(k,i) !Forward connectivity
          !if(iwet(nck)==0 .or. nck>ncells) cycle
          jcn=llec2llec(k,i)  !Backward connectivity
          if(abs(fnx(k,i))>0.999)then !X-direction
            phik(k,i)=phik(k,i)+fintp(k,i)*rpy(jcn,nck)*dphiy(nck) &
                           +(1.0-fintp(k,i))*rpy(k,i)*dphiy(i) !Second order correction
          else   !Y-direction
            phik(k,i)=phik(k,i)+fintp(k,i)*rpx(jcn,nck)*dphix(nck) &
                           +(1.0-fintp(k,i))*rpx(k,i)*dphix(i) !Second order correction
          endif
          phik(k,i)=min(max(phik(k,i),min(phi(i),phi(nck))),max(phi(i),phi(nck))) !Enforce Local Monotonicity
          phik(jcn,nck)=phik(k,i)
        enddo !j
      enddo !ii
!$OMP END PARALLEL DO
!$OMP PARALLEL DO PRIVATE(i,j,k,nck,jcn)  
      do i=1,ncellpoly
        !if(iwet(i)==0) cycle
        !if(maxval(cell2cell(1:ncface(i),i))>ncells) cycle
        !if(minval(iwet(cell2cell(1:ncface(i),i)))==0) cycle
        do j=1,nxyface(i)
          k=kxyface(j,i)
          nck=cell2cell(k,i) !Forward connectivity
          !if(iwet(nck)==0 .or. nck>ncells) cycle
          jcn=llec2llec(k,i)  !Backward connectivity
          if(iwet(i)*iwet(nck)==1)then !Wet-wet boundary
            phik(k,i)=phik(k,i)+fintp(k,i)*(rpx(jcn,nck)*dphix(nck)+rpy(jcn,nck)*dphiy(nck)) &
                 +(1.0-fintp(k,i))*(rpx(k,i)*dphix(i)+rpy(k,i)*dphiy(i)) !Second order cell-reconstruction
            !phik(k,i)=phik(k,i)+fintp(k,i)*(rx(jcn,nck)*dphix(nck)+ry(jcn,nck)*dphiy(nck)) &
            !     +(1.0-fintp(k,i))*(rx(k,i)*dphix(i)+ry(k,i)*dphiy(i)) !Second order face-reconstruction
            phik(k,i)=min(max(phik(k,i),min(phi(i),phi(nck))),max(phi(i),phi(nck))) !Enforce Local Monotonicity
          elseif(iwet(i)==1)then  !Wet-dry
            phik(k,i)=phik(k,i)+rx(k,i)*dphix(i)+ry(k,i)*dphiy(i) !Extrapolation
          elseif(iwet(nck)==1)then  !Dry-wet
            phik(k,i)=phik(k,i)+rx(jcn,nck)*dphix(nck)+ry(jcn,nck)*dphiy(nck) !Extrapolation
          endif  
          phik(jcn,nck)=phik(k,i)
        enddo !j
      enddo !i
!$OMP END PARALLEL DO
    endif
!!$OMP END PARALLEL
    
    return
    end subroutine interp_scal_cell2face 
        
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!  End cell-to-face interpolations
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>    

!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!  Begin cell-to-node interpolations
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
                 
!***************************************************************************    
    subroutine interp_coef_cell2node_invarea(nn,nnc,nmnc,n2c,nc,nD,area,wc2n)
! Interpolation coefficients for a scalar from cell centroids to nodes
! written by Alex Sanchez, USACE-CHL
!***************************************************************************
    use prec_def
    implicit none
    !Input/Output
    integer,    intent(in) :: nn            !Number of nodes
    integer,    intent(in) :: nnc(nn)       !Number neighboring cells for each node
    integer,    intent(in) :: nmnc          !Maximum number of cells per node
    integer,    intent(in) :: n2c(nmnc,nn)  !Forward node-to-cell connectivity
    integer,    intent(in) :: nc            !Number of active cells
    integer,    intent(in) :: nD            !Total number of cells including dummy cells    
    real(ikind),intent(in) :: area(nD)      !Cell area
    real(ikind),intent(out):: wc2n(nmnc,nn) !Weighting coefficients for interpolation
    !Internal Variables
    integer :: i,j,k
    real(ikind) :: wghtsum 
    
    do j=1,nn
      wc2n(:,j) = 0.0
      wghtsum = 0.0
      do k=1,nnc(j)
        i = n2c(k,j)
        if(i<=nc)then
          wc2n(k,j) = 1.0/area(i)
          wghtsum = wghtsum + wc2n(k,j)
        endif
      enddo
      wc2n(:,j) = wc2n(:,j)/wghtsum
    enddo !j
    
    return
    end subroutine interp_coef_cell2node_invarea

!*************************************************************************************    
    subroutine interp_coef_cell2node_invdist(nn,nnc,nmnc,n2c,xn,yn,nc,nD,xc,yc,pow,wc2n)
! Inverse distance interpolation coefficients of a scalar from cell centroids to nodes
! written by Alex Sanchez, USACE-CHL
!*************************************************************************************
    use prec_def
    implicit none
    !Input/Output
    integer,    intent(in) :: nn            !Number of nodes
    integer,    intent(in) :: nnc(nn)       !Number neighboring cells for each node
    integer,    intent(in) :: nmnc          !Maximum number of cells per node
    integer,    intent(in) :: n2c(nmnc,nn)  !Forward node-to-cell connectivity
    real(ikind),intent(in) :: xn(nn),yn(nn) !Node coordinates
    integer,    intent(in) :: nc            !Number of active cells
    integer,    intent(in) :: nD            !Total number of cells including dummy cells    
    real(ikind),intent(in) :: xc(nD),yc(nD) !Cell centroid coordinates
    real(ikind),intent(in) :: pow           !Distance power
    real(ikind),intent(out):: wc2n(nmnc,nn) !Weighting coefficients for interpolation
    !Internal Variables
    integer :: i,j,k
    real(ikind) :: wghtsum 
    
    do j=1,nn
      wc2n(:,j) = 0.0
      wghtsum = 0.0
      do k=1,nnc(j) !# of cells per node
        i = n2c(k,j) !node to cell connectivity
        if(i<=nc)then
          wc2n(k,j) = sqrt((xn(j)-xc(i))**2+(yn(j)-yc(i))**2)**(-pow)
          wghtsum = wghtsum + wc2n(k,j)
        endif
      enddo
      wc2n(:,j) = wc2n(:,j)/wghtsum !Normalize
    enddo !j
    
    return
    end subroutine interp_coef_cell2node_invdist

!*********************************************************************************    
    subroutine interp_coef_cell2node_lstsqrs(nn,nnc,nmnc,n2c,xn,yn,nc,nD,xc,yc,wc2n)
! Least squares interpolation coefficients a scalar from cell centroids to nodes
!
! Reference:
!  Rausch, R.D., Batina, J.T., and Yang, H.T.Y. 1991. Spatial adapation procedures
!    on unstructured meshes for accurate unsteady aerodynamic flow computation. 
!    NASA Langley Research Center Technical Memorandum 104039, 17 p.
!
! written by Alex Sanchez, USACE-CHL
!*********************************************************************************
    use prec_def
    implicit none
    !Input/Output
    integer,    intent(in) :: nn            !Number of nodes
    integer,    intent(in) :: nnc(nn)       !Number neighboring cells for each node
    integer,    intent(in) :: nmnc          !Maximum number of cells per node
    integer,    intent(in) :: n2c(nmnc,nn)  !Forward node-to-cell connectivity
    real(ikind),intent(in) :: xn(nn),yn(nn) !Node coordinates
    integer,    intent(in) :: nc            !Number of active cells
    integer,    intent(in) :: nD            !Total number of cells including dummy cells    
    real(ikind),intent(in) :: xc(nD),yc(nD) !Cell centroid coordinates
    real(ikind),intent(out):: wc2n(nmnc,nn) !Weighting coefficients for interpolation
    !Internal Variables
    integer :: i,j,k
    real(ikind) :: wghtsum,delx,dely
    real(ikind) :: Exx,Exy,Eyy,Emag,Rx,Ry,Gx,Gy
    
    do j=1,nn      
      Rx = 0.0; Ry = 0.0
      Exx = 0.0; Exy = 0.0; Eyy = 0.0
      do k=1,nnc(j)
        i = n2c(k,j)
        if(i<=nc)then
          delx = xn(j)-xc(i)
          dely = yn(j)-yc(i)
          Rx = Rx + delx
          Ry = Ry + dely
          Exx = Exx + delx*delx
          Eyy = Eyy + dely*dely
          Exy = Exy + delx*dely
        endif
      enddo !k
      Emag = Exx*Eyy - Exy**2
      Gx = (Exy*Ry-Exy*Rx)/Emag
      Gy = (Exy*Rx-Exy*Ry)/Emag
      wc2n(:,j) = 0.0
      wghtsum = 0.0
      do k=1,nnc(j)
        i = n2c(k,j)
        if(i<=nc)then
          delx = xn(j)-xc(i)
          dely = yn(j)-yc(i)  
          wc2n(k,j) = 1.0 + Gx*delx + Gy*dely  
          wghtsum = wghtsum + wc2n(k,j)
        endif
      enddo !k
      wc2n(:,j) = wc2n(:,j)/wghtsum
    enddo !j
    
    return
    end subroutine interp_coef_cell2node_lstsqrs           
        
!********************************************************    
    subroutine interp_scal_cell2node(var,scalout,iwritedry)
! Interpolates a scalar from cell centroids to nodes
! written by Alex Sanchez, USACE-CHL
!********************************************************
    use size_def, only: ncells,ncellsD,nnodes
    use geo_def, only: node2cell,nncell
    use flow_def, only: iwet
    use interp_def, only: wc2n
    use prec_def
    implicit none
    !Input/Output
    integer,intent(in) :: iwritedry
    real(ikind),intent(in) :: var(ncellsD)
    real(4),intent(out) :: scalout(nnodes) !Must be single
    !Internal Variables
    integer :: i,j,k
    real(ikind) :: wghtsum
    
    if(iwritedry==1)then !Write dry nodes
      do j=1,nnodes
        scalout(j) = sum(wc2n(1:nncell(j),j)*var(node2cell(1:nncell(j),j)))
      enddo !j
    else !Do not write dry nodes
      do j=1,nnodes
        scalout(j) = 0.0
        wghtsum = 0.0
        do k=1,nncell(j)
          i = node2cell(k,j)
          if(i<=ncells .and. iwet(i)==1)then
            scalout(j) = scalout(j) + wc2n(k,j)*var(i)
            wghtsum = wghtsum + wc2n(k,j)
          endif
        enddo
        if(wghtsum>1.0e-5)then
          scalout(j) = scalout(j)/wghtsum
        else
          scalout(j)=-999.0
        endif
      enddo !j
    endif
    
    return
    end subroutine interp_scal_cell2node 
    
!*************************************************************    
    subroutine interp_vec_cell2node(varx,vary,vecout,iwritedry)
! Interpolates a vector from cell centroids to nodes/vertices
! written by Alex Sanchez, USACE-CHL
!*************************************************************
    use size_def, only: ncells,ncellsD,nnodes
    use geo_def, only: node2cell,nncell,areap
    use flow_def, only: iwet
    use interp_def, only: wc2n
    use prec_def
    implicit none
    !Input/Output
    integer,    intent(in) :: iwritedry
    real(ikind),intent(in) :: varx(ncellsD),vary(ncellsD)
    real(4),intent(out):: vecout(nnodes*2)  !Must be single for XMDF libraries
    !Internal Variables
    integer :: i,j,j2,k
    real(ikind) :: wghtsum
    
    if(iwritedry==1)then !Write dry nodes
      do j=1,nnodes
        j2 = 2*j
        vecout(j2-1) = sum(wc2n(1:nncell(j),j)*varx(node2cell(1:nncell(j),j)))
        vecout(j2)   = sum(wc2n(1:nncell(j),j)*vary(node2cell(1:nncell(j),j)))      
      enddo !j
    else !Do not write dry nodes
      do j=1,nnodes
        j2 = 2*j
        vecout(j2-1:j2) = 0.0
        wghtsum = 0.0
        do k=1,nncell(j)
          i = node2cell(k,j)
          if(i<=ncells .and. iwet(i)==1)then
            vecout(j2-1) = vecout(j2-1) + wc2n(k,j)*varx(i)
            vecout(j2)   = vecout(j2)   + wc2n(k,j)*vary(i) 
            wghtsum = wghtsum + wc2n(k,j)
          endif
        enddo
        if(wghtsum>1.0e-5)then
          vecout(j2-1:j2) = vecout(j2-1:j2)/wghtsum
        else
          vecout(j2-1:j2) = -999.0
        endif
      enddo !j
    endif
    
    return
    end subroutine interp_vec_cell2node

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!  End cell-to-node interpolations
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!  Begin node-to-cell interpolations
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

!********************************************************    
    subroutine interp_scal_node2cell(vtemp,var)
! Interpolates the node values to the cell centers
! This is a temporary routine which is needed because
! SMS only saves nodal values but CMS uses center values
! In the future, SMS should write the cell-center values
! and this subroutine should be removed.
!
! written by Alex Sanchez, USACE-CHL
!********************************************************
    use size_def, only: ncells,ncellsD,ncellsfull
    use geo_def, only: ncface,cell2node
    use prec_def
    implicit none
    !Input/Output
    real(4),     intent(in) :: vtemp(ncellsfull) !Must be single
    real(ikind),intent(out):: var(ncellsD)
    !Internal Variables
    integer :: i,j,k
    
    do i=1,ncells
      var(i) = 0.0
      do k=1,ncface(i)
        j = cell2node(k,i)
        var(i) = var(i) + vtemp(j)
      enddo
      var(i) = var(i)/real(ncface(i),kind=ikind)
    enddo !i
    
    return
    end subroutine interp_scal_node2cell
    
!********************************************************    
    subroutine interp_scal_node2face(vtemp,vark)
! Interpolates the node values to the cell centers
! This is a temporary routine which is needed because
! SMS only saves nodal values but CMS uses center values
! In the future, SMS should write the cell-center values
! and this subroutine should be removed.
!
! written by Alex Sanchez, USACE-CHL
!********************************************************
    use size_def, only: ncells,ncellsD,ncellsfull,nmaxfaces
    use geo_def, only: ncface,cell2node
    use prec_def
    implicit none
    !Input/Output
    real(4),     intent(in) :: vtemp(ncellsfull) !Must be single
    real(ikind),intent(out):: vark(nmaxfaces,ncellsD)
    !Internal Variables
    integer :: i,j,j1,k
    
    do i=1,ncells
      j1 = cell2node(i,1)
      do k=1,ncface(i)
        j = j1
        if(k<ncface(i))then
          j1 = cell2node(i,k+1)
        else
          j1 = cell2node(i,1)
        endif
        vark(k,i) = 0.5*(vtemp(j)+vtemp(j1))
      enddo
    enddo !i
    
    return
    end subroutine interp_scal_node2face
    
!********************************************************    
    subroutine interp_vec_node2cell(vtemp,varx,vary)
! Interpolates the node values to the cell centers
! This is a temporary routine which is needed because
! SMS only saves nodal values but CMS uses center values
! In the future, SMS should write the cell-center values
! and this subroutine should be removed.
!
! written by Alex Sanchez, USACE-CHL
!********************************************************
    use size_def, only: ncells,ncellsD,ncellsfull
    use geo_def, only: ncface,cell2node
    use prec_def
    implicit none
    !Input/Output
    real(4),     intent(in) :: vtemp(ncellsfull*2) !Must be single
    real(ikind),intent(out):: varx(ncellsD),vary(ncellsD)
    !Internal Variables
    integer :: i,j,j2,k
    real(ikind) :: rnfinv
    !real(ikind) :: vecx(ncellsfull),vecy(ncellsfull)
        
    !vecx=vtemp(1:2:ncellsfull*2-1)
    !vecy=vtemp(2:2:ncellsfull*2)
    do i=1,ncells
      varx(i) = 0.0
      vary(i) = 0.0
      do k=1,ncface(i)
        j = cell2node(k,i)
        j2 = j*2        
        varx(i) = varx(i) + vtemp(j2-1)
        vary(i) = vary(i) + vtemp(j2)
      enddo
      rnfinv = 1.0/real(ncface(i),kind=ikind)
      varx(i) = varx(i)*rnfinv
      vary(i) = vary(i)*rnfinv
    enddo !i
    
    return
    end subroutine interp_vec_node2cell

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!  End cell-to-node interpolations
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!  Development
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

!****************************************************************************
    subroutine interp_quad_mapping(px,py,a,b)
! Computes the mapping coefficients for calculating the normalized coordinates
! of a point within a quad.
! Only needs to be calculated once per quad
!
! Author: Alex Sanchez, USACE-CHL
!***************************************************************************
    use prec_def
    implicit none
    !Input
    real(ikind), intent(in) :: px(4),py(4)
    !Output
    real(ikind), intent(out) :: a(4),b(4)
    !Internal
    real(ikind) :: AI(4,4)
    
    !Compute coefficients
    AI(1,:) = (/ 1.0, 0.0, 0.0, 0.0 /)
    AI(2,:) = (/-1.0, 1.0, 0.0, 0.0 /)
    AI(3,:) = (/-1.0, 0.0, 0.0, 1.0 /)
    AI(4,:) = (/ 1.0,-1.0, 1.0,-1.0 /)
      
    a = matmul(AI,px)
    b = matmul(AI,py)

    return
    end subroutine interp_quad_mapping

!*****************************************************************
    subroutine interp_quad_bilinear_coef(a,b,x,y,w)
! Bilinear interpolation for quadrilaterals
!
! Converts physical (x,y) to logical (xi,yi) coordinates
! quadratic equation coeffs, aa*mm^2+bb*yi+cc=0, checks if
! the interpolation coordinate is inside the quad and then
! calculates the bilinear interpolation weights
!
! Author: Alex Sanchez, USACE-CHL
!****************************************************************
    use prec_def
    implicit none
    !Input
    real(ikind),intent(in) :: a(4),b(4) !Mapping coefficients
    real(ikind),intent(in) :: x,y       !Interpolation coordinate
    !Output
    real(ikind),intent(out) :: w(4)     !Interpolation weighhts
    !Internal
    real(ikind):: xi,yi
    real(ikind) :: aa,bb,cc,d
    
    !Compute quadratic formula coefficients
    aa = a(4)*b(3) - a(3)*b(4)
    bb = a(4)*b(1) - a(1)*b(4) + a(2)*b(3) - a(3)*b(2) + x*b(4) - y*a(4)
    cc = a(2)*b(1) - a(1)*b(2) + x*b(2) - y*a(2)

    !Solve quadratiic formula yi = (-b+sqrt(b^2-4ac))/(2a)
    d = sqrt(bb*bb - 4.0_ikind*aa*cc)
    yi = (-bb+d)/(2.0_ikind*aa)

    !Compute xi
    xi = (x-a(1)-a(3)*yi)/(a(2)+a(4)*yi)
    
    !Check if inside quad
    if(xi>=0.0 .and. xi<=1.0 .and. yi>=0.0 .and. yi<=1.0)then !inside quad
        w = (/(1.0_ikind-xi)*(1.0_ikind-yi), &
            xi*(1.0_ikind-yi), &
            xi*yi, &
            (1.0_ikind-xi)*yi /)
    else
        w = 0.0
    endif
    
    return
    end subroutine interp_quad_bilinear_coef
    
end module interp_lib
