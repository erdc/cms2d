!==================================================================================
module der_lib
! Spatial derivatives (gradient and curvature) library
!
! Contains the following:
! - Initialization - 
!     der_grad_cbfd     - Computes the spatial deriviate coefficients using a second-order Finite-Difference Approximation
!     der_grad_cbgg     - Computes the spatial deriviate coefficients using the cell-based Green-Gauss method (first order for skewed cells)
!     der_grad_cbggcr   - Computes the deriviate coefficients using the cell-based Green-Gauss method with linear cell reconstruction
!     der_grad_cbwlsfs  - Weighted cell-based least-squares (face-sharing) gradient coefficients

! - Gradients - 
!     der_grad_eval - Cell-centered gradient operator for a scalar array on arbitrary structured grids at all cells
!     der_gradvec_eval - Cell-centered gradient operator for a vector array structured grids at all cells
!     gradxyvecrecon - Cell-centered gradient operator for a vector array on arbitrary structured grids onlt reconstruction cells
!     dk2d - Calculates the gradient in the cross-cell direction with the positive directions in the outward direction
!     dxy2d -    Cell-centered gradient operator in x- and y-directions at an individual cell 
!     dx2d -     Cell-centered gradient operator in x-direction at an individual cell 
!     dy2d -     Cell-centered gradient operator in y-direction at an individual cell 

! - Curvature -
!     curvxy - Curvature operator for variable scalar arrays on arbitrary structured grids
!     d2xy2d - Calculates the second derivative in two dimensions at an indivudal cell

! - Slope limiters - 
!     limitslope - Differentiable slope limiter for Cartesian grids
!     bjlim - Barth and Jespersen type slope limiter for scalar arrays
!     bjlimvec - Barth and Jespersen type slope limiter for vector arrays
!     lcdlim - Limited Central Difference slope limiter for a scalar arrays
!     lcdlimvec - Limited Central Difference slope limiter for a vector arrays
!     vanleerslopelim - Van leer (Harmonic) slope limiter
!     vanalbadaslopelim - Van Alabada slope limiter
!     musclslopelim - MUSCL (BJ) slope limiter
!     minmodslopelim - Minmod slope limiter
!     minmodfun - Minmod function for BJ and LCD limiting
!     venkatakrishnan - Venkatakrishnan function for BJ slope limiting
!
! written by Alex Sanchez, USACE-CHL
!==================================================================================
    implicit none
    
contains

!********************************************************************
    subroutine der_alloc
! Allocates the derivative operator variables
! Author: Alex Sanchez, USACE-CHL
!********************************************************************    
    use size_def
    use geo_def
    use der_def
    implicit none
    integer :: i
    
    !Gradient operators
    call der_grad_alloc(goa) !for all cells
    call der_grad_alloc(gow) !for wet cells only
    
    !Calculation lists
    allocate(ilistcw(ncellsD)) 
    nlistcw = ncells
    do i=1,ncellsD
      ilistcw(i) = i
    enddo
    
    allocate(ilistgw(ncellsD)) 
    !*Note: simple or regular cells do not have any gradient correction terms        
    if(ncelljoint>0)then !Joint cells
      nlistgw = ncelljoint
      ilistgw(1:nlistgw) = idcelljoint(1:ncelljoint)
    else !Polygonal cells
      nlistgw  = ncells
      do i=1,ncellsD
        ilistgw(i) = i
      enddo
    endif    
    
    !Slope limiters
    call der_lim_alloc
    
    return
    end subroutine der_alloc   

!********************************************************************
    subroutine der_grad_alloc(go)
! Allocates the derivative operator variables
! Author: Alex Sanchez, USACE-CHL
!********************************************************************    
    use size_def
    !use geo_def, only: idcelljoint          never used, commented out   MEB  01/26/2022
    use der_def
    implicit none
    type(der_go_type) :: go  !Gradient Operator    
    
    !--- Derivative Operators ---------------
    allocate(go%ncx(go%nd),       go%ncy(go%nd))
    allocate(go%icx(go%nsc,go%nd),go%icy(go%nsc,go%nd))
    allocate(go%wcx(go%nsc,go%nd),go%wcy(go%nsc,go%nd))      
    go%ncx=0;   go%ncy=0  
    go%icx=0;   go%icy=0
    go%wcx=0.0; go%wcy=0.0
    
    !Only for CBGGCR scheme
    if(go%ider==2)then    
      !Only allocate the variables that will be used (all in this case)
      go%nsg = nmaxfaces+1
      allocate( go%ngx(go%nd),        go%ngy(go%nd))
      allocate( go%igx(go%nsg,go%nd), go%igy(go%nsg,go%nd))
      allocate(go%wgxy(go%nsg,go%nd),go%wgyx(go%nsg,go%nd)) !Perpendicular component
      go%ngy=0;    go%ngx=0
      go%igy=0;    go%igx=0
      go%wgxy=0.0; go%wgyx=0.0
          
      if(ncellpoly>0)then    
        allocate(go%wgxx(go%nsg,go%nd),go%wgyy(go%nsg,go%nd)) !Normal component
        go%wgxx=0.0; go%wgyy=0.0
      endif
    endif
    
    return
    end subroutine der_grad_alloc
    
!********************************************************************
    subroutine der_lim_alloc
! Allocates the derivative operator variables
! Author: Alex Sanchez, USACE-CHL
!********************************************************************    
    use size_def
    use der_def
    implicit none
    
    !Slope limiters
    if(nlim>0)then
      if(ncellsimple>0)then !Cartesian grid
        !Ratio of cell sizes
        allocate(rgx(ncellsD),rgy(ncellsD))
        rgx=1.0; rgy=1.0
      endif  
    endif
    
    return
    end subroutine der_lim_alloc    

!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
! BEGIN Coefficient Calculation routines
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
!********************************************************************
    subroutine der_grad_cbfd(go,ns,nc,idlist)
! Computes the spatial deriviate coefficients using
! a second-order Finite-Difference Approximation
! For regular grids, it reduces to the Central Difference Approximation
! written by Alex Sanchez, USACE-CHL
!*******************************************************************
    use size_def
    use der_def
    use flow_def, only: iwet
    use geo_def, only: ncface,cell2cell,x,y
    !use interp_def, only: fintp                       !never used, commented out   MEB  01/26/2022
    use diag_def
    use diag_lib
    use prec_def
    implicit none
    !Input
    type(der_go_type) :: go !Gradient Operator
    integer,intent(in) :: ns,nc,idlist(ns)
    !Internal variables
    integer :: i,ii,ncn,nce,ncs,ncw,ierr
    real(ikind),dimension(0:4) :: wcxt,wcyt
    integer,    dimension(0:4) :: icxt,icyt
    real(ikind) :: sumwx,sumwy,errx,erry,d1,d2
    real(ikind),parameter :: toldsxy = 1.0e-5
        
    if(ncellsimple/=ncells)then
      call diag_print_error('Cannot use CBFD scheme for Telescoping or Unstructured Grids')  
    endif
    
    do ii=1,nc
      i = idlist(ii)
      if(iwet(i)==0 .and. go%ibc==1)then
        go%ncx(i) = 0
        go%ncy(i) = 0
        !go%wcx(:,i) = 0.0
        !go%wcy(:,i) = 0.0
        cycle
      endif
      
      icxt(0) = i
      icyt(0) = i
      icxt(1:4) = cell2cell(1:4,i)
      icyt(1:4) = cell2cell(1:4,i)
      ncn = cell2cell(1,i)
      nce = cell2cell(2,i)
      ncs = cell2cell(3,i)
      ncw = cell2cell(4,i)
      
      !X-Derivative
      wcxt(1) = 0.0
      wcxt(3) = 0.0
      if(iwet(nce)*iwet(ncw)==1 .or. go%ibc==0)then
        d2 = x(nce)-x(i)
        d1 = x(i)-x(ncw)
        wcxt(2) =  d1/(d2*(d1+d2))
        wcxt(4) = -d2/(d1*(d1+d2))          
        wcxt(0) = -wcxt(4)-wcxt(2) !(d2**2-d1**2)/(d2*d1*(d1+d2))
      elseif(iwet(nce)==1)then !go%ibc==1
        wcxt(2) = 1.0/(x(nce)-x(i))
        wcxt(4) = 0.0
        wcxt(0) = -wcxt(2)
      elseif(iwet(ncw)==1)then !go%ibc==1
        wcxt(2) = 0.0
        wcxt(4) = 1.0/(x(i)-x(ncw))
        wcxt(0) = -wcxt(4)
      else
        wcxt(2) = 0.0
        wcxt(4) = 0.0
        wcxt(0) = 0.0
      endif
        
      !Y-Derivative
      wcyt(2) = 0.0
      wcyt(4) = 0.0
      if(iwet(nce)*iwet(ncs)==1 .or. go%ibc==0)then
        d2 = y(ncn)-y(i)
        d1 = y(i)-y(ncs)
        wcyt(1) =  d1/(d2*(d1+d2))
        wcyt(3) = -d2/(d1*(d1+d2))
        wcyt(0) = -wcyt(3)-wcyt(1) !(d2**2-d1**2)/(d2*d1*(d1+d2))
      elseif(iwet(ncn)==1)then !go%ibc==1
        wcyt(1) = 1.0/(y(ncn)-y(i))
        wcyt(3) = 0.0
        wcyt(0) = -wcyt(1)
      elseif(iwet(ncs)==1)then !go%ibc==1
        wcyt(1) = 0.0
        wcyt(3) = 1.0/(y(i)-y(ncs))
        wcyt(0) = -wcyt(3)
      else
        wcyt(1) = 0.0
        wcyt(3) = 0.0
        wcyt(0) = 0.0
      endif

      !Check weights
      sumwx = sum(wcxt(0:ncface(i)))
      sumwy = sum(wcyt(0:ncface(i)))
      errx = sumwx/maxval(abs(wcxt(0:ncface(i)))) !Normalized error
      erry = sumwy/maxval(abs(wcyt(0:ncface(i)))) !Normalized error
      if(abs(errx)>1.0e-5 .or. abs(erry)>1.0e-5)then
        write(msg2,*,iostat=ierr) 'Cell: ',i
        write(msg3,*,iostat=ierr) 'sumwx,sumwy: ',sumwx,sumwy
        write(msg4,*,iostat=ierr) 'wcxt(0:ncface(i)):',wcxt(0:ncface(i))
        write(msg5,*,iostat=ierr) 'wcyt(0:ncface(i)):',wcyt(0:ncface(i))
        call diag_print_error('Problem calculating gradient coefficients using Cell-based Finite-Difference method',msg2,msg3,msg4,msg5)
      endif
      !Only keep significant weights    
      call screen_weights(4,wcxt,wcyt,icxt,icyt,&
         go%nsc,go%ncx(i),go%ncy(i),&
         go%icx(:,i),go%icy(:,i),go%wcx(:,i),go%wcy(:,i))
      sumwx = sum(go%wcx(1:go%ncx(i),i))
      sumwy = sum(go%wcy(1:go%ncy(i),i)) 
      !Distribute only to significant weights (insures a zero-gradient for a constant function)
      go%wcx(1:go%ncx(i),i) = go%wcx(1:go%ncx(i),i) - sumwx/real(go%ncx(i),kind=ikind)
      go%wcy(1:go%ncy(i),i) = go%wcy(1:go%ncy(i),i) - sumwy/real(go%ncy(i),kind=ikind)
    enddo
    
    return
    end subroutine der_grad_cbfd

!********************************************************************
    subroutine der_grad_cbgg(go,ns,nc,idlist)
! Computes the spatial deriviate coefficients using the 
! cell-based Green-Gauss method (first order for skewed cells)
! written by Alex Sanchez, USACE-CHL
!*******************************************************************
    use size_def
    use der_def
    use flow_def, only: iwet
    use geo_def, only: ncface,cell2cell,areap,dsx,dsy
    use interp_def, only: fintp
    use diag_def
    use diag_lib
    use prec_def
    implicit none
    !Input/Ouput
    type(der_go_type) :: go !Gradient Operator
    integer,intent(in) :: ns,nc,idlist(ns)
    !Internal variables
    integer :: i,ii,k,nck,ierr                            !idval(1),nw
    real(ikind),dimension(0:nmaxfaces) :: wcxt,wcyt
    integer,    dimension(0:nmaxfaces) :: icxt,icyt
    real(ikind) :: sumwx,sumwy,errx,erry,xval,yval
    real(ikind),parameter :: toldsxy = 1.0e-5
    
    do ii=1,nc
      i = idlist(ii)
      if(iwet(i)==0 .and. go%ibc==1)then
        go%ncx(i) = 0
        go%ncy(i) = 0
        !go%wcx(:,i) = 0.0
        !go%wcy(:,i) = 0.0
        cycle
      endif
      
      icxt(0) = i
      icyt(0) = i
      !temporary variables
      wcxt(:) = 0.0
      wcyt(:) = 0.0 
      !if(go%ibc==1)then !Only use wet cells
      !  if(iwet(i)==0)then  
      !    go%ncx(i) = 0
      !    go%ncy(i) = 0
      !    cycle
      !  else !cell is wet
      !    if(sum(iwet(cell2cell(1:ncface(i),i)))>0)then !Any neighboring cells are dry
      !      idval(1) = i
      !      call der_grad_cbwlsfs(go,1,1,idval) !for wet cells neighboring dry cells use CBWLSFS
      !      cycle
      !    endif  
      !  endif
      !endif
      
      !!Cell-based Green-Gauss
      !do k=1,ncface(i)
      !  wcxt(0) = wcxt(0) + dsx(k,i)*(1.0-fintp(k,i))/areap(i)
      !  wcxt(k) = dsx(k,i)*fintp(k,i)/areap(i)
      !  wcyt(0) = wcyt(0) + dsy(k,i)*(1.0-fintp(k,i))/areap(i)
      !  wcyt(k) = dsy(k,i)*fintp(k,i)/areap(i)
      !enddo
      
      !Cell-based Green-Gauss
      do k=1,ncface(i)
        nck = cell2cell(k,i)
        icxt(k) = nck
        icyt(k) = nck
        !if(nck<=ncells .and. (iwet(nck)==1 .or. go%ibc==0))then
        if(iwet(nck)==1 .or. go%ibc==0)then
          wcxt(0) = wcxt(0) + dsx(k,i)*(1.0-fintp(k,i))
          wcxt(k) = dsx(k,i)*fintp(k,i)
          wcyt(0) = wcyt(0) + dsy(k,i)*(1.0-fintp(k,i))
          wcyt(k) = dsy(k,i)*fintp(k,i)
        else
          wcxt(0) = wcxt(0) + dsx(k,i)
          wcyt(0) = wcyt(0) + dsy(k,i)
        endif
      enddo
      wcxt = wcxt/areap(i)
      wcyt = wcyt/areap(i)
      
      !Check weights
      sumwx = sum(wcxt) ; xval = maxval(abs(wcxt))
      sumwy = sum(wcyt) ; yval = maxval(abs(wcyt))
      
      if (xval /= 0.0) errx = sumwx/maxval(abs(wcxt)) !Normalized error
      if (yval /= 0.0) erry = sumwy/maxval(abs(wcyt)) !Normalized error
      
      !if(abs(errx)>1.0e-5 .or. abs(erry)>1.0e-5)then
      if(abs(sumwx)>1.0e-5 .or. abs(sumwy)>1.0e-5)then    
        write(msg2,*,iostat=ierr) i,sumwx,sumwy
        write(msg3,*,iostat=ierr) 'wcxt: ',wcxt
        write(msg4,*,iostat=ierr) 'wcyt: ',wcyt
        call diag_print_error('Problem calculating gradient coefficients using cell-based Green-Gauss method',msg2,msg3,msg4)
      endif
      !Only keep significant weights    
      call screen_weights(nmaxfaces,wcxt,wcyt,icxt,icyt,go%nsc,go%ncx(i),go%ncy(i),&
        go%icx(:,i),go%icy(:,i),go%wcx(:,i),go%wcy(:,i))
      sumwx = sum(go%wcx(1:go%ncx(i),i))
      sumwy = sum(go%wcy(1:go%ncy(i),i)) 
      !Distribute only to significant weights (insures a zero-gradient for a constant function)
      go%wcx(1:go%ncx(i),i) = go%wcx(1:go%ncx(i),i) - sumwx/real(go%ncx(i),kind=ikind)
      go%wcy(1:go%ncy(i),i) = go%wcy(1:go%ncy(i),i) - sumwy/real(go%ncy(i),kind=ikind)
    enddo
    
    return
    end subroutine der_grad_cbgg
    
!*******************************************************************************
    subroutine der_grad_cbggcr(go,ns,nc,idlist)
! Computes the deriviate coefficients using the 
! cell-based Green-Gauss method with linear cell reconstruction
! The method is second order.
! written by Alex Sanchez, USACE-CHL
!*******************************************************************************
    use size_def
    use der_def
    use geo_def
    use flow_def, only: iwet
    use interp_def, only: fintp
    use prec_def
    implicit none
    !Input/Ouput
    type(der_go_type),intent(inout) :: go !Gradient Operator
    integer,intent(in) :: ns,nc,idlist(ns)
    !Internal
    integer :: i,ii,j,k,nck,jcn
    real(ikind),dimension(0:nmaxfaces) :: wgxxt,wgyxt,wgxyt,wgyyt
    integer,    dimension(0:nmaxfaces) :: igxt,igyt
    real(ikind),parameter :: toldsxy = 0.01 !m
    real(ikind),parameter :: toldrpxy = 0.15 !m
    
    !Joint Cells
    if(ncelljoint>0)then
      do ii=1,nc
        i = idlist(ii)  
        if(iwet(i)==0 .and. go%ibc==1)then
          go%ngx(i) = 0
          go%ngy(i) = 0
          !go%wgx(:,i) = 0.0
          !go%wgy(:,i) = 0.0
          cycle
        endif
        igxt(0) = i
        igyt(0) = i
        wgxyt(:)=0.0 !temporary variable
        do j=1,nxface(i) !no repeat east and west sides
          k=kxface(j,i)
          nck=cell2cell(k,i) !Foward connectivity
          if(nck>ncells .or. abs(dy(i)/dy(nck)-1.0)<0.1) cycle !Skip dummy or ghost cells for gradient correction
          if(iwet(nck)==0 .and. go%ibc==1) cycle
          jcn=llec2llec(k,i)  !Backwards connectivity
          wgxyt(0)=wgxyt(0)+dsx(k,i)*(1.0-fintp(k,i))*rpy(k,i)/areap(i) !Note: first component is cell center  
          wgxyt(k)=dsx(k,i)*fintp(k,i)*rpy(jcn,nck)/areap(i)
          igxt(k)=nck
        enddo
        wgyxt(:)=0.0 !temporary variable    
        do j=1,nyface(i) !no repeat north and south faces
          k=kyface(j,i)
          nck=cell2cell(k,i) !Foward connectivity
          igyt(k)=nck
          if(nck>ncells .or. abs(dx(i)/dx(nck)-1.0)<0.1) cycle !Skip dummy or ghost cells for gradient correction
          if(iwet(nck)==0 .and. go%ibc==1) cycle
          jcn=llec2llec(k,i)  !Backwards connectivity
          wgyxt(0)=wgyxt(0)+dsy(k,i)*(1.0-fintp(k,i))*rpx(k,i)/areap(i) !Note: first component is cell center
          wgyxt(k)=dsy(k,i)*fintp(k,i)*rpx(jcn,nck)/areap(i)
        enddo !k
        !Only keep significant weights
        call screen_weights(nmaxfaces,wgxyt,wgyxt,igxt,igyt,go%nsg,go%ngx(i),go%ngy(i),&
                go%igx(:,i),go%igy(:,i),go%wgxy(:,i),go%wgyx(:,i))
      enddo !ii
    else !Polyhedral
      do ii=1,nc
        i = idlist(ii)
        if(iwet(i)==0 .and. go%ibc==1)then
          go%ngx(i) = 0
          go%ngy(i) = 0
          !go%wcx(:,i) = 0.0
          !go%wcy(:,i) = 0.0
          cycle
        endif
        wgxxt(:)=0.0; wgyxt(:)=0.0 !temporary variables
        wgxyt(:)=0.0; wgyyt(:)=0.0 !temporary variables
        go%ngx(i)=1
        go%ngy(i)=1
        go%igx(1,i)=i
        go%igy(1,i)=i
        do k=1,ncface(i)
          nck=cell2cell(k,i) !Foward connectivity
          if(nck>ncells) cycle !Skip dummy or ghost cells for gradient correction
          if(iwet(nck)==0 .and. go%ibc==1) cycle
          jcn=llec2llec(k,i)  !Backwards connectivity   
          go%ngx(i) = go%ngx(i) + 1
          go%ngy(i) = go%ngy(i) + 1
          !For X-Derivative
          wgxxt(0) = wgxxt(0) + dsx(k,i)*(1.0-fintp(k,i))*rpx(k,i)/areap(i) !Note: first component is cell center
          wgxxt(k) =            dsx(k,i)*fintp(k,i)*rpx(jcn,nck)/areap(i)
          wgxyt(0) = wgxyt(0) + dsx(k,i)*(1.0-fintp(k,i))*rpy(k,i)/areap(i) !Note: first component is cell center
          wgxyt(k) =            dsx(k,i)*fintp(k,i)*rpy(jcn,nck)/areap(i)
          !For Y-Derivative
          wgyxt(0) = wgyxt(0) + dsy(k,i)*(1.0-fintp(k,i))*rpx(k,i)/areap(i) !Note: first component is cell center
          wgyxt(k) =            dsy(k,i)*fintp(k,i)*rpx(jcn,nck)/areap(i)
          wgyyt(0) = wgyyt(0) + dsy(k,i)*(1.0-fintp(k,i))*rpy(k,i)/areap(i) !Note: first component is cell center
          wgyyt(k) =            dsy(k,i)*fintp(k,i)*rpy(jcn,nck)/areap(i)       
          
          !!For X-Derivative
          !wgxxt(0) = wgxxt(0) + dsx(k,i)*(1.0-fintp(k,i))*rx(k,i)/areap(i) !Note: first component is cell center
          !wgxxt(k) = wgxxt(k) + dsx(k,i)*fintp(k,i)*rx(jcn,nck)/areap(i)
          !wgxyt(0) = wgxyt(0) + dsx(k,i)*(1.0-fintp(k,i))*ry(k,i)/areap(i) !Note: first component is cell center
          !wgxyt(k) = wgxyt(k) + dsx(k,i)*fintp(k,i)*ry(jcn,nck)/areap(i)
          !!For Y-Derivative
          !wgyxt(0) = wgyxt(0) + dsy(k,i)*(1.0-fintp(k,i))*rx(k,i)/areap(i) !Note: first component is cell center
          !wgyxt(k) = wgyxt(k) + dsy(k,i)*fintp(k,i)*rx(jcn,nck)/areap(i)
          !wgyyt(0) = wgyyt(0) + dsy(k,i)*(1.0-fintp(k,i))*ry(k,i)/areap(i) !Note: first component is cell center
          !wgyyt(k) = wgyyt(k) + dsy(k,i)*fintp(k,i)*ry(jcn,nck)/areap(i)    
          
          go%igx(go%ngx(i),i) = nck
          go%igy(go%ngy(i),i) = nck
          go%wgxx(go%ngx(i),i) = wgxxt(k)
          go%wgxy(go%ngx(i),i) = wgxyt(k)
          go%wgyx(go%ngy(i),i) = wgyxt(k)
          go%wgyy(go%ngy(i),i) = wgyyt(k)
        enddo !k
        go%wgxx(1,i) = wgxxt(0)
        go%wgxy(1,i) = wgxyt(0)
        go%wgyx(1,i) = wgyxt(0)
        go%wgyy(1,i) = wgyyt(0)
      enddo !i
    endif
    
    return
    end subroutine der_grad_cbggcr
    
!*************************************************************************
    subroutine der_grad_cbwlsfs(go,ns,nc,idlist)
! Weighted Cell-based least-squares (face-sharing) gradient coefficients
! Author: Alex Sanchez, USACE-CHL
!*************************************************************************
    use size_def
    use geo_def
    use der_def
    use diag_def
    use diag_lib
    use flow_def, only: iwet
    use prec_def
    implicit none
    !Input/Ouput
    type(der_go_type) :: go !Gradient Operator
    integer,intent(in) :: ns,nc,idlist(ns)
    !Internal
    integer :: i,ii,k,kk,nck,nw,next,idcell(1)
    real(ikind),dimension(go%nsc-1) :: delx,dely
    real(ikind),dimension(go%nsc-1) :: w  !Weights
    integer,    dimension(0:go%nsc-1) :: icxt,icyt
    real(ikind),dimension(0:go%nsc-1) :: wcxt,wcyt
    real(ikind) :: Rxx,Rxy,Ryy,Rdet,Rdetinv
    real(ikind) :: sumwx,sumwy,errx,erry,valx,valy
    real(ikind),parameter :: toldist = 0.15 !m
    real(ikind),parameter :: tolcond = 0.5
    real(ikind),parameter :: tolRder = 1.0e-3
    logical :: extsten
    
    !Note: The stencils for the x- and y-directions
    ! do not have to be same.
    ! In order to reduce the stencil size for each
    ! gradient, only the cells which or offset in 
    ! coordinate of the gradient direction
    ! (e.g. for the x-direction only the left
    ! and right cells are used).
    
    nw = go%nsc-1
    next = 0
    do ii=1,nc
      i = idlist(ii)
      !if(i==83)then
      if(i==4226)then    
        continue
      endif
      if(iwet(i)==0 .and. go%ibc==1)then
        go%ncx(i) = 0
        go%ncy(i) = 0
        !go%wcx(:,i) = 0.0
        !go%wcy(:,i) = 0.0
        cycle
      endif
    
      icxt(0) = i
      icyt(0) = i
        
      !Compute moments and other variables
      Rxx = 0.0; Ryy = 0.0; Rxy = 0.0
      do k=1,ncface(i)
        nck = cell2cell(k,i)
        !if(nck>ncells) cycle
        if(iwet(nck)==1 .or. go%ibc==0)then
          delx(k) = x(nck)-x(i)
          dely(k) = y(nck)-y(i)
          w(k) = 1.0/sqrt(delx(k)**2+dely(k)**2)**npow
          delx(k) = w(k)*delx(k)
          dely(k) = w(k)*dely(k)
          Rxx = Rxx + delx(k)**2
          Ryy = Ryy + dely(k)**2
          Rxy = Rxy + delx(k)*dely(k)
        endif
      enddo
       
      extsten = .false. !Estended stencil
      
      !Check conditioning
      Rdet = Rxx*Ryy - Rxy**2
      valx = Rxx/(Rxx+Ryy)
      valy = Ryy/(Rxx+Ryy)
      !if(valx<tolcond .or. valy<tolcond .or. abs(Rdet)<tolRder)then
      !if(maxval(cell2cell(1:ncface(i),i))<=ncells)then
      if(ncellpoly>0)then
        extsten = .true.
        !Add extended stencil to improve conditioning
        next = next + 1
        do kk=1,ncnode(i)
          k = ncface(i) + kk
          nck = cell2extcell(kk,i)
          !if(nck>ncells) cycle
          if(iwet(nck)==1 .or. go%ibc==0)then
            delx(k) = x(nck)-x(i)
            dely(k) = y(nck)-y(i)
            w(k) = 1.0/sqrt(delx(k)**2+dely(k)**2)**npow
            delx(k) = w(k)*delx(k)
            dely(k) = w(k)*dely(k)
            Rxx = Rxx + delx(k)**2
            Ryy = Ryy + dely(k)**2
            Rxy = Rxy + delx(k)*dely(k)
          endif
        enddo
      endif
      
      Rdet = Rxx*Ryy - Rxy**2
      valx = Rxx/(Rxx+Ryy)
      valy = Ryy/(Rxx+Ryy)
      if(abs(Rdet)<1.0e-10)then
      !if(valx<tolcond .or. valy<tolcond .or. abs(Rdet)<tolRder)then
        if(allocated(mapid))then   
          write(msg2,*) '  Cell: ',mapid(i)
        else
          write(msg2,*) '  Cell: ',i
        endif
        call diag_print_warning('Problem in Cell-based weighted least-squares method',&
           '  Ill conditioned matrix',msg2,'  Calculating using CBGG')
        idcell(1) = i
        call der_grad_cbgg(go,1,1,idcell)
      endif

      Rdetinv = 1.0/Rdet
        
      !Compute all derivative weights
      wcxt(:)=0.0; wcyt(:)=0.0  !Initialize
      do k=1,ncface(i)
        nck=cell2cell(k,i)
        !if(nck>ncells) cycle
        if(iwet(nck)==1 .or. go%ibc==0)then
          wcxt(k) = w(k)*(Ryy*delx(k)-Rxy*dely(k))*Rdetinv
          wcxt(0) = wcxt(0)-wcxt(k)
          wcyt(k) = w(k)*(Rxx*dely(k)-Rxy*delx(k))*Rdetinv
          wcyt(0) = wcyt(0)-wcyt(k)
          icxt(k) = nck
          icyt(k) = nck
        endif
      enddo
      if(extsten)then
        do kk=1,ncnode(i)
          k = ncface(i) + kk
          nck = cell2extcell(kk,i)
          !if(nck>ncells) cycle
          if(iwet(nck)==1 .or. go%ibc==0)then
            wcxt(k) = w(k)*(Ryy*delx(k)-Rxy*dely(k))*Rdetinv
            wcxt(0) = wcxt(0)-wcxt(k)
            wcyt(k) = w(k)*(Rxx*dely(k)-Rxy*delx(k))*Rdetinv
            wcyt(0) = wcyt(0)-wcyt(k)
            icxt(k) = nck
            icyt(k) = nck
          endif
        enddo
      endif
      !Check weights
      sumwx = sum(wcxt(0:nw))
      sumwy = sum(wcyt(0:nw))
      !Normalized error
      errx = sumwx/maxval(abs(wcxt(0:nw)))
      erry = sumwy/maxval(abs(wcyt(0:nw)))
      if(abs(errx)>1.0e-5 .or. abs(erry)>1.0e-5)then
        write(msg2,*) i,sumwx,sumwy
        write(msg3,*) wcxt(0:nw),wcyt(0:nw)  
        call diag_print_error('Problem calculating gradient coefficients with Weighted Cell-based least-squares method',msg2,msg3)
      endif
      !Only keep significant weights
      call screen_weights(nw,wcxt,wcyt,icxt,icyt,go%nsc,go%ncx(i),go%ncy(i),&
        go%icx(:,i),go%icy(:,i),go%wcx(:,i),go%wcy(:,i))
      sumwx = sum(go%wcx(1:go%ncx(i),i))
      sumwy = sum(go%wcy(1:go%ncy(i),i)) 
      !Distribute only to significant weights (insures a zero-gradient for a constant function)
      go%wcx(1:go%ncx(i),i) = go%wcx(1:go%ncx(i),i) - sumwx/real(go%ncx(i),kind=ikind)
      go%wcy(1:go%ncy(i),i) = go%wcy(1:go%ncy(i),i) - sumwy/real(go%ncy(i),kind=ikind)
    enddo
    
    return
    end subroutine der_grad_cbwlsfs    

!***********************************************************************
    subroutine screen_weights(nt,wxt,wyt,ixt,iyt,ns,nx,ny,ix,iy,wx,wy)
!***********************************************************************
    !use size_def, only: ncellsD,ndmaxfaces,nmaxfaces            ! never used, commented out   MEB  01/26/2022
    !use geo_def, only: ncface,cell2cell                         ! never used, commented out   MEB  01/26/2022
    use prec_def
    implicit none
    !Input
    integer,    intent(in) :: nt
    real(ikind),intent(in),dimension(0:nt):: wxt,wyt
    integer,    intent(in),dimension(0:nt):: ixt,iyt
    !Output
    integer,    intent(in) :: ns
    integer,    intent(out):: nx,ny
    integer,    intent(out),dimension(ns):: ix,iy
    real(ikind),intent(out),dimension(ns):: wx,wy
    !Internal variables
    integer :: k
    !real(ikind), parameter :: tolw = 5.0e-5
    real(ikind), parameter :: tolw = 1.0e-8
    
    nx = 0; ny = 0
    do k=0,nt
      if(abs(wxt(k))>tolw)then
        nx = nx + 1
        ix(nx) = ixt(k)
        wx(nx) = wxt(k)
      endif
      if(abs(wyt(k))>tolw)then
        ny = ny + 1
        iy(ny) = iyt(k)
        wy(ny) = wyt(k)
      endif
    enddo   
      
    return
    end subroutine screen_weights    

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! END Coefficient Calculation routines    
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>    
    
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
! BEGIN Gradient Operator routines    
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
!*******************************************************
    subroutine der_grad_eval(go,ilim,var,dvarx,dvary)
! Derivative operator for a scalar variable
! on arbitrary structured grids
! written by Alex Sanchez, USACE-CHL
!*******************************************************
    use size_def, only: ncells,ncellsD,ncelljoint,ncellpoly
    use geo_def, only: idcelljoint                               !ncface,cell2cell   never used, commented out   MEB  01/26/2022
    use der_def    
    use prec_def
    implicit none
    !Input/Output
    type(der_go_type),intent(in) :: go !Gradient Operator
    integer, intent(in) :: ilim !Slope limiter
    real(ikind),intent(in),dimension(ncellsD) :: var !Scalar variable    
    !Output    
    real(ikind),intent(inout),dimension(ncellsD) :: dvarx,dvary
    !Internal variables
    integer :: i,ii,k,nck
    real(ikind),dimension(ncellsD) :: dvarxm,dvarym
    !logical :: isnankind

!$OMP PARALLEL
!--- Save input derivatives for deferred corrections -----------------
    if(go%ider==2)then 
      !Save previous iteration values  
!$OMP DO PRIVATE(i,ii)
      do ii=1,ncelljoint
        i=idcelljoint(ii)    
        dvarxm(i) = dvarx(i)
        dvarym(i) = dvary(i)
      enddo 
!$OMP END DO
!$OMP DO PRIVATE(i,ii)
      do i=1,ncellpoly !No index mapping needed (all cells are polygonal)
        dvarxm(i) = dvarx(i)
        dvarym(i) = dvary(i)
      enddo 
!$OMP END DO
    endif
    
!--- Gradient Calculation -----------------------------------    
!Verbose version
!$OMP DO PRIVATE(i,k)
    do i=1,ncells
      dvarx(i) = 0.0
      do k=1,go%ncx(i)
        dvarx(i) = dvarx(i) + go%wcx(k,i)*var(go%icx(k,i))
      enddo
      dvary(i) = 0.0
      do k=1,go%ncy(i)
        dvary(i) = dvary(i) + go%wcy(k,i)*var(go%icy(k,i))
      enddo
    enddo 
!$OMP END DO

!!Compact version
!!$OMP DO PRIVATE(i,k)
!!      do i=1,ncells
!!        dvarx(i)=sum(go%wcx(1:go%ncx(i),i)*var(go%icx(1:go%ncx(i),i)))
!!        dvary(i)=sum(go%wcy(1:go%ncy(i),i)*var(go%icy(1:go%ncy(i),i)))
!!      enddo 
!!$OMP END DO
     
!--- Deferred gradient corrections -------------------------
    if(go%ider==2)then
!Verbose option        
!$OMP DO PRIVATE(i,ii,k,nck)                 
      do ii=1,ncelljoint
        i=idcelljoint(ii)
        do k=1,go%ngx(i)
          dvarx(i) = dvarx(i) + go%wgxy(k,i)*dvarym(go%igx(k,i)) !Perpendicular component
        enddo
        do k=1,go%ngy(i)
          dvary(i) = dvary(i) + go%wgyx(k,i)*dvarxm(go%igy(k,i)) !Perpendicular component
        enddo
      enddo  
!$OMP END DO
!$OMP DO PRIVATE(i,k,nck)                 
      do i=1,ncellpoly !No index mapping needed (all cells are polygonal)
        do k=1,go%ngx(i)
          nck = go%igx(k,i)
          dvarx(i) = dvarx(i) + go%wgxx(k,i)*dvarxm(nck) + go%wgxy(k,i)*dvarym(nck)
        enddo
        do k=1,go%ngy(i)
          nck = go%igy(k,i)
          dvary(i) = dvary(i) + go%wgyx(k,i)*dvarxm(nck) + go%wgyy(k,i)*dvarym(nck)
        enddo
      enddo
!$OMP END DO        
        
!!Compaction option
!!$OMP DO PRIVATE(i,ii,k,nck)
!      do ii=1,ncelljoint
!        i=idcelljoint(ii)       
!        dvarx(i)=dvarx(i)+sum(go%wgxy(1:go%ngx(i),i)*dvarym(go%igx(1:go%ngx(i),i)))
!        dvary(i)=dvary(i)+sum(go%wgyx(1:go%ngy(i),i)*dvarxm(go%igy(1:go%ngy(i),i)))
!      enddo  
!!$OMP END DO
!!$OMP DO PRIVATE(i,k,nck)                 
!      do i=1,ncellpoly !No index mapping needed (all cells are polygonal)
!        dvarx(i)=dvarx(i)+sum(go%wgxx(1:go%ngx(i),i)*dvarxm(go%igx(1:go%ngx(i),i))) &
!                         +sum(go%wgxy(1:go%ngx(i),i)*dvarym(go%igx(1:go%ngx(i),i)))
!        dvary(i)=dvary(i)+sum(go%wgyx(1:go%ngy(i),i)*dvarxm(go%igy(1:go%ngy(i),i))) &
!                         +sum(go%wgyy(1:go%ngy(i),i)*dvarym(go%igy(1:go%ngy(i),i)))
!      enddo
!!$OMP END DO
    endif     
!$OMP END PARALLEL

!--- Slope Limters ----------------------------
    if(ilim>0) call slopelim(var,ilim,dvarx,dvary)

    return
    end subroutine der_grad_eval

!*******************************************************
    subroutine der_grad_eval_recon(go,ilim,var,dvarx,dvary)
! Gradient operator for variable arrays
! on arbitrary structured grids
! Gradients are only calculated at cells where a
! reconstruction (cell- or face-based) is needed
!
! written by Alex Sanchez, USACE-CHL
!*******************************************************
    use size_def, only: ncellsD,ncelljoint,ncellpoly         !ncells             never used, commented out   MEB  01/26/2022
    use geo_def, only: idcelljoint                           !ncface,cell2cell   never used, commented out   MEB  01/26/2022
    use der_def    
    use prec_def
    implicit none
    !Input/Output
    type(der_go_type),intent(in) :: go !Gradient Operator
    integer, intent(in) :: ilim !Slope limiter
    real(ikind),intent(in),dimension(ncellsD) :: var    
    real(ikind),intent(out),dimension(ncellsD) :: dvarx,dvary
    !Internal variables
    integer :: i,ii,k !,nck
    real(ikind) :: dvarxm(ncellsD),dvarym(ncellsD)
    !logical :: isnankind

    if(ncellpoly>0)then
      call der_grad_eval(go,ilim,var,dvarx,dvary) !All cells are reconstruction cells
      return
    endif
    
!$OMP PARALLEL
!!$OMP DO PRIVATE(i,ii)              
!      do ii=1,ncellsimple
!        i=idcellsimple(ii)
!        dvarx(i)=0.0; dvary(i)=0.0 !Not computed
!      enddo       
!!$OMP END DO

!--- Save previous iteration values -----------
    if(go%ider==2)then
!$OMP DO PRIVATE(i,ii)
      do ii=1,ncelljoint
        i=idcelljoint(ii)    
        dvarxm(i) = dvarx(i)
        dvarym(i) = dvary(i)
      enddo !i
!$OMP END DO
    endif
    
!$OMP DO PRIVATE(i,ii,k)
    do ii=1,ncelljoint
      i=idcelljoint(ii)
      dvarx(i) = 0.0
      do k=1,go%ncx(i)
        dvarx(i) = dvarx(i) + go%wcx(k,i)*var(go%icx(k,i))
      enddo
      dvary(i) = 0.0
      do k=1,go%ncy(i)
        dvary(i) = dvary(i) + go%wcy(k,i)*var(go%icy(k,i))
      enddo
      !dvarx(i)=sum(go%wcx(1:go%ncx(i),i)*var(go%icx(1:go%ncx(i),i)))
      !dvary(i)=sum(go%wcy(1:go%ncy(i),i)*var(go%icy(1:go%ncy(i),i)))
    enddo !i
!$OMP END DO
    
     
    !Deferred gradient corrections  
    if(go%ider==2)then
!Verbose option        
!$OMP DO PRIVATE(i,ii,k)                 
      do ii=1,ncelljoint
        i=idcelljoint(ii)
        do k=1,go%ngx(i)
          dvarx(i) = dvarx(i) + go%wgxy(k,i)*dvarym(go%igx(k,i)) !Perpendicular component
        enddo
        do k=1,go%ngy(i)
          dvary(i) = dvary(i) + go%wgyx(k,i)*dvarxm(go%igy(k,i)) !Perpendicular component
        enddo
      enddo  
!$OMP END DO
!!Compaction option
!!$OMP DO PRIVATE(i,ii,k,nck)
!      do ii=1,ncelljoint
!        i=idcelljoint(ii)       
!        dvarx(i)=dvarx(i)+sum(go%wgxy(1:go%ngx(i),i)*dvarym(go%igx(1:go%ngx(i),i)))
!        dvary(i)=dvary(i)+sum(go%wgyx(1:go%ngy(i),i)*dvarxm(go%igy(1:go%ngy(i),i)))
!      enddo  
!!$OMP END DO
    endif !nder
!$OMP END PARALLEL

    return
    end subroutine der_grad_eval_recon    

!*****************************************************************************
    subroutine der_gradvec_eval(go,ilim,u,v,dux,duy,dvx,dvy)
! Derivative operator for a scalar variable
! on arbitrary structured grids
! written by Alex Sanchez, USACE-CHL
!*****************************************************************************
    use size_def, only: ncells,ncellsD,ncelljoint,ncellpoly
    use geo_def, only: idcelljoint                          !ncface, cell2cell   never used, commented out   MEB  01/26/2022
    use der_def
    use prec_def
    implicit none
    !Input/Output
    type(der_go_type),intent(in) :: go !Gradient Operator    
    integer, intent(in) :: ilim !Slope limiter
    real(ikind),intent(in),dimension(ncellsD) :: u,v
    real(ikind),intent(out),dimension(ncellsD) :: dux,duy,dvx,dvy
    !Internal variables
    integer :: i,ii,k,nck
    real(ikind),dimension(ncellsD) :: duxm,duym,dvxm,dvym
    !logical :: isnankind

!$OMP PARALLEL
!--- Save input derivatives for deferred corrections ------------------
    if(go%ider==2)then 
      !Save previous iteration values  
!$OMP DO PRIVATE(i,ii)
      do ii=1,ncelljoint
        i=idcelljoint(ii)    
        duxm(i) = dux(i)
        duym(i) = duy(i)
        dvxm(i) = dvx(i)
        dvym(i) = dvy(i)
      enddo !ii
!$OMP END DO
!$OMP DO PRIVATE(i,ii)
      do i=1,ncellpoly !No index mapping needed (all cells are polygonal)
        duxm(i) = dux(i)
        duym(i) = duy(i)
        dvxm(i) = dvx(i)
        dvym(i) = dvy(i)
      enddo !i
!$OMP END DO
    endif

!-- Calculate derivatives -------------------------------------
!Verbose version
!$OMP DO PRIVATE(i,k,nck)
    do i=1,ncells
      dux(i) = 0.0
      dvx(i) = 0.0
      do k=1,go%ncx(i)
        nck = go%icx(k,i)
        dux(i) = dux(i) + go%wcx(k,i)*u(nck)
        dvx(i) = dvx(i) + go%wcx(k,i)*v(nck)
      enddo
      duy(i) = 0.0
      dvy(i) = 0.0
      do k=1,go%ncy(i)
        nck = go%icy(k,i)
        duy(i) = duy(i) + go%wcy(k,i)*u(nck)
        dvy(i) = dvy(i) + go%wcy(k,i)*v(nck)
      enddo
    enddo !i
!$OMP END DO

!!Compact version
!!$OMP DO PRIVATE(i)
!    do i=1,ncells
!      dux(i)=sum(go%wcx(1:go%ncx(i),i)*u(go%icx(1:go%ncx(i),i)))
!      duy(i)=sum(go%wcy(1:go%ncx(i),i)*u(go%icy(1:go%ncy(i),i)))
!      dvx(i)=sum(go%wcx(1:go%ncx(i),i)*v(go%icx(1:go%ncx(i),i)))
!      dvy(i)=sum(go%wcy(1:go%ncx(i),i)*v(go%icy(1:go%ncy(i),i)))
!    enddo !i
!!$OMP END DO
     
!--- Deferred gradient corrections -----------------------
    if(go%ider==2)then
!$OMP DO PRIVATE(i,ii,k,nck)                 
        do ii=1,ncelljoint
          i = idcelljoint(ii)    
          do k=1,go%ngx(i)
            nck = go%igx(k,i)
            dux(i) = dux(i) + go%wgxy(k,i)*duym(nck) !Perpendicular component
            dvx(i) = dvx(i) + go%wgxy(k,i)*dvym(nck) !Perpendicular component
          enddo !k
          do k=1,go%ngy(i)
            nck = go%igy(k,i)
            duy(i) = duy(i) + go%wgyx(k,i)*duxm(nck) !Perpendicular component
            dvy(i) = dvy(i) + go%wgyx(k,i)*dvxm(nck) !Perpendicular component
          enddo !k
        enddo !ii
!$OMP END DO
!$OMP DO PRIVATE(i,k,nck)
        do i=1,ncellpoly !No index mapping needed (all cells are polygonal)        
          do k=1,go%ngx(i)
            nck = go%igx(k,i)
            dux(i) = dux(i) + go%wgxx(k,i)*duxm(nck) + go%wgxy(k,i)*duym(nck)
            dvx(i) = dvx(i) + go%wgxx(k,i)*dvxm(nck) + go%wgxy(k,i)*dvym(nck)
          enddo !k
          do k=1,go%ngy(i)
            nck = go%igy(k,i)
            duy(i) = duy(i)+go%wgyx(k,i)*duxm(nck) + go%wgyy(k,i)*duy(nck)
            dvy(i) = dvy(i)+go%wgyx(k,i)*dvxm(nck) + go%wgyy(k,i)*dvy(nck)
          enddo !k
        enddo !i
!$OMP END DO
    endif     
!$OMP END PARALLEL

!--- Slope Limters ----------------------------
    if(ilim>0) call slopelimvec(u,v,ilim,dux,duy,dvx,dvy)
    
    return
    end subroutine der_gradvec_eval    

!******************************************************************
    subroutine der_gradvec_eval_recon(go,ilim,u,v,dux,duy,dvx,dvy)
! Derivative operator for a vector at reconstruction cells
! on arbitrary structured grids
! written by Alex Sanchez, USACE-CHL
!******************************************************************
    use size_def, only: ncells,ncellsD,ncelljoint,ncellpoly
    use geo_def, only: idcelljoint                              !ncface, cell2cell   never used, commented out   MEB  01/26/2022
    use der_def    
    use prec_def
    implicit none
    !Input/Output
    type(der_go_type),intent(in) :: go !Gradient Operator
    integer, intent(in) :: ilim !Slope limiter
    real(ikind),intent(in),dimension(ncellsD) :: u,v
    real(ikind),intent(out),dimension(ncellsD) :: dux,duy,dvx,dvy
    !Internal variables
    integer :: i,ii,k,nck
    real(ikind),dimension(ncellsD) :: duxm,duym,dvxm,dvym
    !logical :: isnankind

    if(ncellpoly>0)then
      call der_gradvec_eval(go,ilim,u,v,dux,duy,dvx,dvy) !All cells are reconstruction cells
      return
    endif
    
!$OMP PARALLEL
!--- Save input derivatives for deferred corrections ----------------------
    if(go%ider==2)then 
      !Save previous iteration values  
!$OMP DO PRIVATE(i,ii)
      do ii=1,ncelljoint
        i=idcelljoint(ii)    
        duxm(i) = dux(i)
        duym(i) = duy(i)
        dvxm(i) = dvx(i)
        dvym(i) = dvy(i)
      enddo !i
!$OMP END DO
    endif
    
!--- Calculate derivatives ------------------------------------
!Verbose version
!$OMP DO PRIVATE(i,k,nck)
    do i=1,ncells
      dux(i) = 0.0
      dvx(i) = 0.0
      do k=1,go%ncx(i)
        nck = go%icx(k,i)
        dux(i) = dux(i) + go%wcx(k,i)*u(nck)
        dvx(i) = dvx(i) + go%wcx(k,i)*v(nck)
      enddo
      duy(i) = 0.0
      dvy(i) = 0.0
      do k=1,go%ncy(i)
        nck = go%icy(k,i)
        duy(i) = duy(i) + go%wcy(k,i)*u(nck)
        dvy(i) = dvy(i) + go%wcy(k,i)*v(nck)
      enddo
    enddo !i
!$OMP END DO

!!Compact version
!!$OMP DO PRIVATE(i,ii)
!    do ii=1,ncelljoint
!      i=idcelljoint(ii)  
!      dux(i)=sum(go%wcx(1:go%ncx(i),i)*u(go%icx(1:go%ncx(i),i)))
!      duy(i)=sum(go%wcy(1:go%ncx(i),i)*u(go%icy(1:go%ncy(i),i)))
!      dvx(i)=sum(go%wcx(1:go%ncx(i),i)*v(go%icx(1:go%ncx(i),i)))
!      dvy(i)=sum(go%wcy(1:go%ncx(i),i)*v(go%icy(1:go%ncy(i),i)))
!    enddo !ii
!!$OMP END DO
     
!--- Deferred gradient corrections -----------------------
    if(go%ider==2)then
!$OMP DO PRIVATE(i,ii,k,nck)                 
        do ii=1,ncelljoint
          i=idcelljoint(ii)    
          do k=1,go%ngx(i)
            nck = go%igx(k,i)
            dux(i) = dux(i) + go%wgxy(k,i)*duym(nck) !Perpendicular component
            dvx(i) = dvx(i) + go%wgxy(k,i)*dvym(nck) !Perpendicular component
          enddo !k
          do k=1,go%ngy(i)
            nck = go%igy(k,i)
            duy(i) = duy(i) + go%wgyx(k,i)*duxm(nck) !Perpendicular component
            dvy(i) = dvy(i) + go%wgyx(k,i)*dvxm(nck) !Perpendicular component
          enddo !k
        enddo !ii
!$OMP END DO
!$OMP DO PRIVATE(i,k,nck)                 
        do i=1,ncellpoly !No index mapping needed (all cells are polygonal)        
          do k=1,go%ngx(i)
            nck = go%igx(k,i)
            dux(i) = dux(i) + go%wgxx(k,i)*duxm(nck) + go%wgxy(k,i)*duym(nck)
            dvx(i) = dvx(i) + go%wgxx(k,i)*dvxm(nck) + go%wgxy(k,i)*dvym(nck)
          enddo !k
          do k=1,go%ngy(i)
            nck = go%igy(k,i)
            duy(i) = duy(i) + go%wgyx(k,i)*duxm(nck) + go%wgyy(i,k)*duy(nck)
            dvy(i) = dvy(i) + go%wgyx(k,i)*dvxm(nck) + go%wgyy(i,k)*dvy(nck)
          enddo !k
        enddo !i
!$OMP END DO
    endif   
!$OMP END PARALLEL

!--- Slope Limters ----------------------------
    if(ilim>0) call slopelimvec(u,v,ilim,dux,duy,dvx,dvy)

    return
    end subroutine der_gradvec_eval_recon    
    
!**********************************************************************
    subroutine dxy2d(go,i,var,dvarx,dvary)
! Gradient Operator in x- and y-directions
!
! Revision History: 
!   11/27/2012 Alex Sanchez USACE-CHL
!     Changed to use new gradient variables for efficiency
!***********************************************************************
    use size_def
    use der_def
    use prec_def
    implicit none
    !Input
    type(der_go_type),intent(in) :: go !Gradient Operator
    integer,          intent(in) :: i  !Cell-id    
    real(ikind),      intent(in) :: var(ncellsD) !Variable used for gradient 
    !Output
    real(ikind),      intent(out):: dvarx,dvary  !Gradients in x- and y-directions
        
    dvarx = sum(go%wcx(1:go%ncx(i),i)*var(go%icx(1:go%ncx(i),i)))
    dvary = sum(go%wcy(1:go%ncy(i),i)*var(go%icy(1:go%ncy(i),i)))
    
    return
    end subroutine dxy2d
  
!**********************************************************************
    subroutine dx2d(go,i,var,dvarx) 
! Gradient Operator in x-direction
!
! Revision History: 
!   11/27/2012 Alex Sanchez USACE-CHL
!     Changed to use new gradient variables for efficiency
!***********************************************************************
    use size_def
    use der_def
    use prec_def
    implicit none
    !Input
    type(der_go_type),intent(in) :: go !Gradient Operator
    integer,    intent(in) :: i            !Cell-id
    real(ikind),intent(in) :: var(ncellsD) !Variable used for gradient 
    !Output
    real(ikind),intent(out):: dvarx        !Gradients in x-directions
        
    dvarx = sum(go%wcx(1:go%ncx(i),i)*var(go%icx(1:go%ncx(i),i)))
    
    return
    end subroutine dx2d
      
!**********************************************************************
    subroutine dy2d(go,i,var,dvary) 
! Gradient Operator in y-direction
!
! Revision History: 
!   11/27/2012 Alex Sanchez USACE-CHL
!     Changed to use new gradient variables for efficiency
!***********************************************************************
    use size_def
    use der_def
    use prec_def
    implicit none
    !Input
    type(der_go_type),intent(in) :: go !Gradient Operator
    integer,          intent(in) :: i  !Cell-id    
    real(ikind),      intent(in) :: var(ncellsD) !Variable used for gradient 
    !Output
    real(ikind),      intent(out):: dvary        !Gradients in y-direction
      
    dvary = sum(go%wcy(1:go%ncy(i),i)*var(go%icy(1:go%ncy(i),i)))
    
    return
    end subroutine dy2d
    
!******************************************************************    
    function dk2d(i,k,phi,dphix,dphiy) result(gradk)
! Calculates the gradient in the cross-cell direction
! with the positive directions in the outward direction
! written by Alex Sanchez, USACE
!******************************************************************
    use size_def
    use geo_def, only: cell2cell,llec2llec,dn,rpx,rpy          !ncface   never used, commented out   MEB  01/26/2022
    use prec_def
    implicit none
    integer :: i,k,j,n
    real(ikind):: phi(ncellsD),dphix(ncellsD),dphiy(ncellsD),gradk
    
    n=cell2cell(k,i)
    if(n<=ncells)then !Cell-reconstruction      
      j=llec2llec(k,i)
      gradk=(phi(n)+rpx(n,j)*dphix(n)+rpy(n,j)*dphiy(n) &
           -(phi(i)+rpx(k,i)*dphix(i)+rpy(k,i)*dphiy(i)))/dn(k,i)  
    else !No reconstruction
      gradk=(phi(n)-phi(i))/dn(k,i)   
    endif
    
    return
    end function dk2d   

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>    
! END Gradient operator routines
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
! BEGIN Curvature operator routines    
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<    

!*******************************************************
    subroutine curvxy(go,dvarx,dvary,d2varxx,d2varyy)
! Curvature operator for variable arrays
! on arbitrary structured grids
! written by Alex Sanchez, USACE-CHL
!*******************************************************
    use size_def, only: ncells,ncellsD
    use der_def    
    use prec_def
    implicit none
    !Input/Output
    type(der_go_type),intent(in) :: go !Gradient Operator
    real(ikind),      intent(in) :: dvarx(ncellsD),dvary(ncellsD)
    real(ikind),     intent(out) :: d2varxx(ncellsD),d2varyy(ncellsD)
    !Internal variables
    integer :: i

!$OMP PARALLEL DO PRIVATE(i)    
      do i=1,ncells
        d2varxx(i) = sum(go%wcx(1:go%ncx(i),i)*dvarx(go%icx(1:go%ncx(i),i)))
        d2varyy(i) = sum(go%wcy(1:go%ncx(i),i)*dvary(go%icy(1:go%ncy(i),i)))
      enddo 
!$OMP END PARALLEL DO
    
    return
    end subroutine curvxy      
    
!*****************************************************************    
    subroutine d2xy2d(go,i,dvarx,dvary,d2varix2,d2variy2) 
! Calculates the second derivative in two dimensions
! written by Alex Sanchez, USACE-CHL
!*****************************************************************
    use size_def, only: ncellsD
    use der_def
    use prec_def
    implicit none
    !Input
    type(der_go_type),intent(in) :: go !Gradient Operator
    integer,intent(in) :: i
    real(ikind),intent(in) :: dvarx(ncellsD),dvary(ncellsD)
    !Output
    real(ikind),intent(out) :: d2varix2,d2variy2
    
    call dx2d(go,i,dvarx,d2varix2)
    call dy2d(go,i,dvary,d2variy2)
    
    return
    end subroutine d2xy2d

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! END Curvature Operators
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
! BEGIN Slope Limiters
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

!*****************************************************************
    subroutine slopelim(var,ilim,dvarx,dvary)
! Applies a slope limiter to the scalar array var
! written by Alex Sanchez, USACE-CHL
!*****************************************************************
    use size_def
    use prec_def
    implicit none
    integer,    intent(in) :: ilim
    real(ikind),intent(in),dimension(ncellsD) :: var
    real(ikind),intent(inout),dimension(ncellsD) :: dvarx,dvary
    
    select case(ilim)
    case(1); call limitslope(minmodslopelim,var,dvarx,dvary)
    case(2); call limitslope(vanleerslopelim,var,dvarx,dvary)
    case(3); call limitslope(vanalbadaslopelim,var,dvarx,dvary)
    case(4); call limitslope(musclslopelim,var,dvarx,dvary)
    case(5); call bjlim(minmodfun,var,dvarx,dvary)
    case(6); call bjlim(venkatakrishnan,var,dvarx,dvary)
    case(7); call lcdlim(minmodfun,var,dvarx,dvary)
    case(8); call lcdlim(venkatakrishnan,var,dvarx,dvary)
    end select
    
    return
    end subroutine slopelim
    
!*****************************************************************
    subroutine slopelimvec(u,v,ilim,dux,duy,dvx,dvy)
! Applies a slope limiter to the vector array var
! written by Alex Sanchez, USACE-CHL
!*****************************************************************
    use size_def
    use prec_def
    implicit none
    integer,    intent(in) :: ilim
    real(ikind),intent(in),dimension(ncellsD) :: u,v
    real(ikind),intent(inout),dimension(ncellsD) ::dux,duy,dvx,dvy
    
    select case(ilim)
    case(1); call limitslopevec(minmodslopelim,u,v,dux,duy,dvx,dvy)
    case(2); call limitslopevec(vanleerslopelim,u,v,dux,duy,dvx,dvy)
    case(3); call limitslopevec(vanalbadaslopelim,u,v,dux,duy,dvx,dvy)
    case(4); call limitslopevec(musclslopelim,u,v,dux,duy,dvx,dvy)
    case(5); call bjlimvec(minmodfun,u,v,dux,duy,dvx,dvy)
    case(6); call bjlimvec(venkatakrishnan,u,v,dux,duy,dvx,dvy)
    end select
    
    return
    end subroutine slopelimvec
    
!*****************************************************************
    function vanleerslopelim(r) result(phi)
! Van leer (Harmonic) slope limiter
! Reference:
!   Van Leer, B. (1974). Towards the ultimate conservative 
!     difference scheme II. Monotonicity and conservation combined  
!     in a second order scheme. J. Comp. Phys., 14, p361-70.
!*****************************************************************
    use prec_def
    implicit none
    real(ikind), intent(in) :: r
    real(ikind) :: phi,rr
    
    rr = max(r,0.0)
    phi = 4.0*rr/(rr+1.0)**2
    
    end function vanleerslopelim

!*****************************************************************
    function vanalbadaslopelim(r) result(phi)
! Van Alabada slope limiter
!
! Reference:
!   Van Albada, G.D., Van Leer, B., and Roberts, W.W. (1982). 
!     A comparative study of computational methods in cosmic 
!     gas dynamics, Astron. Astrophysics, 108, p76.
!*****************************************************************
    use prec_def
    implicit none
    real(ikind), intent(in) :: r
    real(ikind) :: phi,rr
    
    rr = max(r,0.0)
    phi = 2.0*rr/(rr**2+1.0)
    
    end function vanalbadaslopelim

!*****************************************************************
    function musclslopelim(r) result(phi)
! MUSCL (BJ) slope limiter
!*****************************************************************
    use prec_def
    implicit none
    real(ikind), intent(in) :: r
    real(ikind) :: phi,b,rr
    
    rr = max(r,0.0)
    b = 4.0/(rr+1.0)
    phi = min(rr*b,b,1.0)
    
    end function musclslopelim

!*****************************************************************
    function minmodslopelim(r) result(phi)
!Minmod slope limiter
!*****************************************************************
    use prec_def
    implicit none
    real(ikind), intent(in) :: r
    real(ikind) :: phi,b,rr
    
    rr = max(r,0.0)
    b = 2.0/(rr+1.0)
    phi = min(rr*b,b)
    
    end function minmodslopelim    
    
!*****************************************************************
    function minmodfun(r) result(psi)
! Minmod function for BJ limiting
! Reference:
!   Roe, P.L. (1986). Characteristic-based schemes for the 
!     Euler equations, Annual Review in Fluid Mechanics, 18, p337.
!*****************************************************************
    use prec_def
    real(ikind),intent(in) :: r
    real(ikind) :: psi
    
    psi=max(0.0,min(1.0,r))
    
    end function minmodfun    
    
!*****************************************************************
    function venkatakrishnan(r) result(psi)
! Venkatakrishnan function for BJ slope limiting
! Reference:
!   Venkatakrishnan, V. 1993. On the accuracy of limiters and
!      convergence to steady state solutions, AIAA Paper 93-0880.
!*****************************************************************
    use prec_def
    implicit none
    real(ikind), intent(in) :: r
    real(ikind) :: psi
    real(ikind),parameter :: e2 = 0.001 !Hard-coded for simplicity
    
    psi = max(0.0,(r*r+2.0*r+e2)/(r*r+r+2+e2))
    
    end function venkatakrishnan
        
!**********************************************************************    
    subroutine bjlim(func,phi,dphix,dphiy)
! Barth and Jespersen type slope limiter for a scalar array
!
! BJ used the minmod or venkatakrishnan functions for 
! limiting the face gradients. 
! The venkatakrishnan function produced better convergence 
! since it is differenciable
!
! Reference:
!   Barth, T.J., and Jespersen, D.C., 1989. The Design and Application
!     of Upwind Schemes on Unstructured Meshes. AIAA Paper, 89-0366.
!
! written by Alex Sanchez, USACE-CHL
!**********************************************************************
    use size_def
    use geo_def
    use prec_def
    implicit none
    !Input/Output
    real(ikind),intent(in),dimension(ncellsD) :: phi
    real(ikind),intent(inout),dimension(ncellsD) :: dphix,dphiy
    !Internal Variables
    integer :: i,ii,j,k,nck
    real(ikind) :: phimax,phimin,delmax,delmin,delface,lim,r
    
    interface
      function func(r) result(psi)
        use prec_def
        real(ikind),intent(in) :: r
        real(ikind) :: psi
      end function
    endinterface

!$OMP PARALLEL
!$OMP DO PRIVATE(i,ii,j,k,nck,phimax,phimin,delmax,delmin,delface,lim,r)    
    do ii=1,ncellsimple
      i=idcellsimple(ii)
      !X-Direction
      phimax=phi(i);  phimin=phi(i)
      do j=1,nxface(i)
        k=kxface(j,i)
        nck=cell2cell(k,i)
        phimax=max(phimax,phi(nck))
        phimin=min(phimin,phi(nck))
      enddo
      delmax=phimax-phi(i) !>=0
      delmin=phimin-phi(i) !<=0
      lim=1.0
      do j=1,nxface(i)
        k=kxface(j,i)  
        delface=rx(k,i)*dphix(i)
        if(delface>1.0e-6)then
          r=delmax/delface !Double positive
          lim=min(lim,func(r)) !Limiter
        elseif(delface<-1.0e-6)then
          r=delmin/delface !Double negative
          lim=min(lim,func(r)) !Limiter
        endif        
      enddo
      dphix(i)=dphix(i)*lim
      !Y-Direction
      phimax=phi(i); phimin=phi(i)
      do j=1,nyface(i)
        k=kyface(j,i)
        nck=cell2cell(k,i)
        phimax=max(phimax,phi(nck))
        phimin=min(phimin,phi(nck))
      enddo
      delmax=phimax-phi(i) !>=0
      delmin=phimin-phi(i) !<=0
      lim=1.0
      do j=1,nyface(i)
        k=kyface(j,i)  
        delface=ry(k,i)*dphiy(i)
        if(delface>1.0e-6)then
          r=delmax/delface !Double positive
          lim=min(lim,func(r)) !Limiter
        elseif(delface<-1.0e-6)then
          r=delmin/delface !Double negative
          lim=min(lim,func(r)) !Limiter
        endif
      enddo
      dphiy(i)=dphiy(i)*lim
    enddo
!$OMP END DO NOWAIT

!$OMP DO PRIVATE(i,ii,k,nck,phimax,phimin,delmax,delmin,delface,lim,r)    
    do ii=1,ncelljoint
      i=idcelljoint(ii)
      phimax=phi(i); phimin=phi(i)
      do k=1,ncface(i)
        nck=cell2cell(k,i)
        phimax=max(phimax,phi(nck))
        phimin=min(phimin,phi(nck))
      enddo
      delmax=phimax-phi(i) !>=0
      delmin=phimin-phi(i) !<=0
      lim=1.0
      do k=1,ncface(i)
        delface=rx(k,i)*dphix(i)+ry(k,i)*dphiy(i)
        if(delface>1.0e-6)then
          r=delmax/delface !Double positive
          lim=min(lim,func(r)) !Limiter
        elseif(delface<-1.0e-6)then
          r=delmin/delface !Double negative
          lim=min(lim,func(r)) !Limiter
        endif
      enddo
      dphix(i)=dphix(i)*lim
      dphiy(i)=dphiy(i)*lim
    enddo
!$OMP END DO NOWAIT

!$OMP DO PRIVATE(i,k,nck,phimax,phimin,delmax,delmin,delface,lim,r)    
    do i=1,ncellpoly !No index mapping required (all cells polygonal)
      phimax=phi(i); phimin=phi(i)
      do k=1,ncface(i)
        nck=cell2cell(k,i)
        phimax=max(phimax,phi(nck))
        phimin=min(phimin,phi(nck))
      enddo
      delmax=phimax-phi(i) !>=0
      delmin=phimin-phi(i) !<=0
      lim=1.0
      do k=1,ncface(i)
        delface=rx(k,i)*dphix(i)+ry(k,i)*dphiy(i)
        if(delface>1.0e-6)then
          r=delmax/delface !Double positive
          lim=min(lim,func(r)) !Limiter
        elseif(delface<-1.0e-6)then
          r=delmin/delface !Double negative
          lim=min(lim,func(r)) !Limiter
        endif
      enddo
      dphix(i)=dphix(i)*lim
      dphiy(i)=dphiy(i)*lim
    enddo
!$OMP END DO
!$OMP END PARALLEL

    return    
    end subroutine bjlim

!**********************************************************************    
    subroutine bjlimvec(func,u,v,dux,duy,dvx,dvy)
! Barth and Jespersen Type Slope Limiter for vectors
!
! BJ used the minmod or venkatakrishnan functions for 
! limiting the face gradients. 
! The venkatakrishnan function produced better convergence 
! since it is differenciable
!
! Reference:
!   Barth, T.J., and Jespersen, D.C., 1989. The Design and Application
!     of Upwind Schemes on Unstructured Meshes. AIAA Paper, 89-0366.
!
! written by Alex Sanchez, USACE-CHL
!**********************************************************************
    use size_def
    use geo_def
    use prec_def
    implicit none
    !Input/Output
    real(ikind),intent(in),dimension(ncellsD) :: u,v
    real(ikind),intent(inout),dimension(ncellsD) :: dux,duy,dvx,dvy
    !Internal Variables
    integer :: i,ii,j,k,nck
    real(ikind) :: delface,r
    real(ikind) :: umax,umin,udelmax,udelmin,ulim
    real(ikind) :: vmax,vmin,vdelmax,vdelmin,vlim
    
    interface
      function func(r) result(psi)
        use prec_def
        real(ikind),intent(in) :: r
        real(ikind) :: psi
      end function
    endinterface
    
!$OMP PARALLEL
!$OMP DO PRIVATE(i,ii,j,k,nck,umax,umin,udelmax,udelmin,ulim,vmax,vmin,vdelmax,vdelmin,vlim,delface,r)
    do ii=1,ncellsimple
      i=idcellsimple(ii)
      !X-Direction
      umax=u(i); umin=u(i)
      vmax=u(i); vmin=v(i)
      do j=1,nxface(i)
        k=kxface(j,i)
        nck=cell2cell(k,i)
        umax=max(umax,u(nck)); vmax=max(vmax,v(nck)) 
        umin=min(umin,u(nck)); vmin=min(vmin,v(nck))
      enddo
      udelmax=umax-u(i); vdelmax=vmax-v(i) 
      udelmin=umin-u(i); vdelmin=vmin-v(i)
      ulim=1.0; vlim=1.0
      do j=1,nxface(i)
        k=kxface(j,i)
        !U-Component
        delface=rx(k,i)*dux(i)  
        if(delface>1.0e-6)then
          r=udelmax/delface !Double positive
          ulim=min(ulim,func(r)) !Limiter
        elseif(delface<-1.0e-6)then
          r=udelmin/delface !Double negative
          ulim=min(ulim,func(r)) !Limiter
        endif
        !V-component
        delface=rx(k,i)*dvx(i)
        if(delface>1.0e-6)then
          r=vdelmax/delface !Double positive
          vlim=min(vlim,func(r)) !Limiter
        elseif(delface<-1.0e-6)then
          r=vdelmin/delface !Double negative
          vlim=min(vlim,func(r)) !Limiter
        endif 
      enddo
      dux(i)=dux(i)*ulim !Apply limiter
      dvx(i)=dvx(i)*vlim !Apply limiter
      !Y-Direction
      umax=u(i); vmax=v(i)
      umin=u(i); vmin=v(i)
      do j=1,nyface(i)
        k=kyface(j,i)
        nck=cell2cell(k,i)
        umax=max(umax,u(nck)); vmax=max(vmax,v(nck)) 
        umin=min(umin,u(nck)); umin=min(vmin,v(nck))
      enddo
      udelmax=umax-u(i); vdelmax=vmax-v(i) 
      udelmin=umin-u(i); vdelmin=vmin-v(i)
      ulim=1.0; vlim=1.0
      do j=1,nyface(i)
        k=kyface(j,i)  
        !U-Component
        delface=ry(k,i)*duy(i)   
        if(delface>1.0e-6)then
          r=udelmax/delface !Double positive
          ulim=min(ulim,func(r)) !Limiter
        elseif(delface<-1.0e-6)then
          r=udelmin/delface !Double negative
          ulim=min(ulim,func(r)) !Limiter
        endif
        !V-component
        delface=ry(k,i)*dvy(i)
        if(delface>1.0e-6)then
          r=vdelmax/delface !Double positive
          vlim=min(vlim,func(r)) !Limiter
        elseif(delface<-1.0e-6)then
          r=vdelmin/delface !Double negative
          vlim=min(vlim,func(r)) !Limiter
        endif
      enddo
      duy(i)=duy(i)*ulim
      dvy(i)=dvy(i)*vlim
    enddo
!$OMP END DO NOWAIT

!$OMP DO PRIVATE(i,ii,k,nck,umax,umin,udelmax,udelmin,ulim,vmax,vmin,vdelmax,vdelmin,vlim,delface,r)
    do ii=1,ncelljoint
      i=idcelljoint(ii)
      umax=u(i); umin=u(i)
      vmax=u(i); vmin=v(i)
      do k=1,ncface(i)
        nck=cell2cell(k,i)
        umax=max(umax,u(nck)); vmax=max(vmax,v(nck)) 
        umin=min(umin,u(nck)); vmin=min(vmin,v(nck))
      enddo
      udelmax=umax-u(i); vdelmax=vmax-v(i) 
      udelmin=umin-u(i); vdelmin=vmin-v(i)
      ulim=1.0; vlim=1.0
      do k=1,ncface(i)
        !U-Component
        delface=rx(k,i)*dux(i)+ry(k,i)*duy(i)   
        if(delface>1.0e-6)then
          r=udelmax/delface !Double positive
          ulim=min(ulim,func(r)) !Limiter
        elseif(delface<-1.0e-6)then
          r=udelmin/delface !Double negative
          ulim=min(ulim,func(r)) !Limiter
        endif
        !V-component
        delface=rx(k,i)*dvx(i)+ry(k,i)*dvy(i)
        if(delface>1.0e-6)then
          r=vdelmax/delface !Double positive
          vlim=min(vlim,func(r)) !Limiter
        elseif(delface<-1.0e-6)then
          r=vdelmin/delface !Double negative
          vlim=min(vlim,func(r)) !Limiter
        endif
      enddo
      dux(i)=dux(i)*ulim
      duy(i)=duy(i)*ulim
      dvx(i)=dvx(i)*vlim
      dvy(i)=dvy(i)*vlim
    enddo
!$OMP END DO NOWAIT

!$OMP DO PRIVATE(i,k,nck,umax,umin,udelmax,udelmin,ulim,vmax,vmin,vdelmax,vdelmin,vlim,delface,r)
    do i=1,ncellpoly
      umax=u(i); umin=u(i)
      vmax=u(i); vmin=v(i)
      do k=1,ncface(i)
        nck=cell2cell(k,i)
        umax=max(umax,u(nck)); vmax=max(vmax,v(nck)) 
        umin=min(umin,u(nck)); vmin=min(vmin,v(nck))
      enddo
      udelmax=umax-u(i); vdelmax=vmax-v(i) 
      udelmin=umin-u(i); vdelmin=vmin-v(i)
      ulim=1.0; vlim=1.0
      do k=1,ncface(i)
        !U-Component
        delface=rx(k,i)*dux(i)+ry(k,i)*duy(i)   
        if(delface>1.0e-6)then
          r=udelmax/delface !Double positive
          ulim=min(ulim,func(r)) !Limiter
        elseif(delface<-1.0e-6)then
          r=udelmin/delface !Double negative
          ulim=min(ulim,func(r)) !Limiter
        endif
        !V-component
        delface=rx(k,i)*dvx(i)+ry(k,i)*dvy(i)
        if(delface>1.0e-6)then
          r=vdelmax/delface !Double positive
          vlim=min(vlim,func(r)) !Limiter
        elseif(delface<-1.0e-6)then
          r=vdelmin/delface !Double negative
          vlim=min(vlim,func(r)) !Limiter
        endif
      enddo
      dux(i)=dux(i)*ulim
      duy(i)=duy(i)*ulim
      dvx(i)=dvx(i)*vlim
      dvy(i)=dvy(i)*vlim
    enddo
!$OMP END DO
!$OMP END PARALLEL

    return    
    end subroutine bjlimvec

!************************************************************    
    subroutine lcdlim(func,phi,dphix,dphiy)
! Limited Central Difference Slope Limiter
!
! Uses the minmod or venkatakrishnan functions for 
! limiting the face gradients. 
! The venkatakrishnan function produced better convergence 
! since it is differenciable
!
! Reference:
!   Hubbard, M.E. 1999. Multidimensional slope limiters for 
!     MUSCL-Type Finite Volume Schemes on Unstructured Grids.
!     Journal of Computational Physics, 155, 54–74. 
!
! written by Alex Sanchez, USACE-CHL
!************************************************************
    use size_def
    use geo_def
    use comp_lib
    use prec_def
    implicit none
    !Input/Output
    real(ikind),intent(in),dimension(ncellsD) :: phi
    real(ikind),intent(inout),dimension(ncellsD) :: dphix,dphiy
    !Internal Variables
    integer :: i,ii,j,k,nck
    real(ikind) :: delmax,delmin,delface,lim,r
    
    interface
      function func(r) result(psi)
        use prec_def
        real(ikind),intent(in) :: r
        real(ikind) :: psi
      end function
    endinterface

!$OMP PARALLEL
!$OMP DO PRIVATE(i,ii,j,k,nck,delmax,delmin,delface,lim,r)
    do ii=1,ncellsimple
      i=idcellsimple(ii)
      !X-Direction
      lim=1.0
      do j=1,nxface(i)
        k=kxface(j,i)
        nck=cell2cell(k,i)  
        delmax=max(phi(nck)-phi(i),0.0)
        delmin=min(phi(nck)-phi(i),0.0) 
        delface=rx(k,i)*dphix(i)
        if(delface>delmax)then
          r=delmax/delface !Double positive
          lim=min(lim,func(r)) !Limiter
        elseif(delface<delmin)then
          r=delmin/delface !Double negative
          lim=min(lim,func(r)) !Limiter
        endif        
      enddo
      dphix(i)=dphix(i)*lim
      !Y-Direction
      lim=1.0
      do j=1,nyface(i)
        k=kyface(j,i)
        nck=cell2cell(k,i)  
        delmax=max(phi(nck)-phi(i),0.0)
        delmin=min(phi(nck)-phi(i),0.0) 
        delface=ry(k,i)*dphiy(i)
        if(delface>delmax)then
          r=delmax/delface !Double positive
          lim=min(lim,func(r)) !Limiter
        elseif(delface<delmin)then
          r=delmin/delface !Double negative
          lim=min(lim,func(r)) !Limiter
        endif
      enddo
      dphiy(i)=dphiy(i)*lim
    enddo
!$OMP END DO NOWAIT

!$OMP DO PRIVATE(i,ii,k,nck,delmax,delmin,delface,lim,r)
    do ii=1,ncelljoint
      i=idcelljoint(ii)
      lim=1.0
      do k=1,ncface(i)
        nck=cell2cell(k,i)  
        delmax=max(phi(nck)-phi(i),0.0)
        delmin=min(phi(nck)-phi(i),0.0)  
        delface=rx(k,i)*dphix(i)+ry(k,i)*dphiy(i)
        if(delface>1.0e-6)then
          r=delmax/delface !Double positive
          lim=min(lim,func(r)) !Limiter
        elseif(delface<-1.0e-6)then
          r=delmin/delface !Double negative
          lim=min(lim,func(r)) !Limiter
        endif
      enddo
      dphix(i)=dphix(i)*lim
      dphiy(i)=dphiy(i)*lim
    enddo
!$OMP END DO NOWAIT

!$OMP DO PRIVATE(i,k,nck,delmax,delmin,delface,lim,r)
    do i=1,ncellpoly !No index mapping required (all cells polygonal)
      lim=1.0
      do k=1,ncface(i)
        nck=cell2cell(k,i)  
        delmax=max(phi(nck)-phi(i),0.0)
        delmin=min(phi(nck)-phi(i),0.0)
        delface=rx(k,i)*dphix(i)+ry(k,i)*dphiy(i)
        if(delface>delmax)then
          r=delmax/delface !Double positive
          lim=min(lim,func(r)) !Limiter
        elseif(delface<delmin)then
          r=delmin/delface !Double negative
          lim=min(lim,func(r)) !Limiter
        endif
      enddo
      dphix(i)=dphix(i)*lim
      dphiy(i)=dphiy(i)*lim
    enddo
!$OMP END DO
!$OMP END PARALLEL

    return    
    end subroutine lcdlim

!************************************************************    
    subroutine lcdlimvec(func,u,v,dux,duy,dvx,dvy)
! Barth and Jespersen Type Slope Limiter for vectors
!
! Uses the minmod or venkatakrishnan functions for 
! limiting the face gradients. 
! The venkatakrishnan function produced better convergence 
! since it is differenciable
!
! Reference:
!   Hubbard, M.E. 1999. Multidimensional slope limiters for 
!     MUSCL-Type Finite Volume Schemes on Unstructured Grids.
!     Journal of Computational Physics, 155, 54–74. 
!
! written by Alex Sanchez, USACE-CHL
!************************************************************
    use size_def
    use geo_def
    use prec_def
    implicit none
    !Input/Output
    real(ikind),intent(in),dimension(ncellsD) :: u,v
    real(ikind),intent(inout),dimension(ncellsD) :: dux,duy,dvx,dvy
    !Internal Variables
    integer :: i,ii,j,k,nck
    real(ikind) :: delmax,delmin,ulim,vlim,delface,r
    
    interface
      function func(r) result(psi)
        use prec_def
        real(ikind),intent(in) :: r
        real(ikind) :: psi
      end function
    endinterface
    
!$OMP PARALLEL
!$OMP DO PRIVATE(i,ii,j,k,nck,delmax,delmin,delface,r,ulim,vlim)
    do ii=1,ncellsimple
      i=idcellsimple(ii)
      !X-Direction
      ulim=1.0; vlim=1.0
      do j=1,nxface(i)
        k=kxface(j,i)
        nck=cell2cell(k,i)
        !U-Component
        delmax=max(u(nck)-u(i),0.0)
        delmin=min(u(nck)-u(i),0.0)
        delface=rx(k,i)*dux(i)
        if(delface>delmax)then
          r=delmax/delface !Double positive
          ulim=min(ulim,func(r)) !Limiter
        elseif(delface<delmin)then
          r=delmin/delface !Double negative
          ulim=min(ulim,func(r)) !Limiter
        endif
        !V-component
        delmax=max(v(nck)-v(i),0.0)
        delmin=min(v(nck)-v(i),0.0)
        delface=rx(k,i)*dvx(i)
        if(delface>delmax)then
          r=delmax/delface !Double positive
          vlim=min(vlim,func(r)) !Limiter
        elseif(delface<delmin)then
          r=delmin/delface !Double negative
          vlim=min(vlim,func(r)) !Limiter
        endif 
      enddo
      dux(i)=dux(i)*ulim !Apply limiter
      dvx(i)=dvx(i)*vlim !Apply limiter
      !Y-Direction
      ulim=1.0; vlim=1.0
      do j=1,nyface(i)
        k=kyface(j,i)  
        nck=cell2cell(k,i)
        !U-Component
        delmax=max(u(nck)-u(i),0.0)
        delmin=min(u(nck)-u(i),0.0)
        delface=ry(k,i)*duy(i)   
        if(delface>delmax)then
          r=delmax/delface !Double positive
          ulim=min(ulim,func(r)) !Limiter
        elseif(delface<delmin)then
          r=delmin/delface !Double negative
          ulim=min(ulim,func(r)) !Limiter
        endif
        !V-component
        delmax=max(v(nck)-v(i),0.0)
        delmin=min(v(nck)-v(i),0.0)
        delface=ry(k,i)*dvy(i)
        if(delface>delmax)then
          r=delmax/delface !Double positive
          vlim=min(vlim,func(r)) !Limiter
        elseif(delface<delmin)then
          r=delmin/delface !Double negative
          vlim=min(vlim,func(r)) !Limiter
        endif
      enddo
      duy(i)=duy(i)*ulim
      dvy(i)=dvy(i)*vlim
    enddo
!$OMP END DO NOWAIT

!$OMP DO PRIVATE(i,ii,k,nck,delmax,delmin,delface,r,ulim,vlim)
    do ii=1,ncelljoint
      i=idcelljoint(ii)
      ulim=1.0; vlim=1.0
      do k=1,ncface(i)
        nck=cell2cell(k,i)
        !U-Component
        delmax=max(u(nck)-u(i),0.0)
        delmin=min(u(nck)-u(i),0.0)
        delface=rx(k,i)*dux(i)+ry(k,i)*duy(i)   
        if(delface>delmax)then
          r=delmax/delface !Double positive
          ulim=min(ulim,func(r)) !Limiter
        elseif(delface<delmin)then
          r=delmin/delface !Double negative
          ulim=min(ulim,func(r)) !Limiter
        endif
        !V-component
        delmax=max(v(nck)-v(i),0.0)
        delmin=min(v(nck)-v(i),0.0)
        delface=rx(k,i)*dvx(i)+ry(k,i)*dvy(i)
        if(delface>delmax)then
          r=delmax/delface !Double positive
          vlim=min(vlim,func(r)) !Limiter
        elseif(delface<delmin)then
          r=delmin/delface !Double negative
          vlim=min(vlim,func(r)) !Limiter
        endif
      enddo
      dux(i)=dux(i)*ulim
      duy(i)=duy(i)*ulim
      dvx(i)=dvx(i)*vlim
      dvy(i)=dvy(i)*vlim
    enddo
!$OMP END DO NOWAIT

!$OMP DO PRIVATE(i,k,nck,delmax,delmin,delface,r,ulim,vlim)
    do i=1,ncellpoly
      ulim=1.0; vlim=1.0
      do k=1,ncface(i)
        nck=cell2cell(k,i)  
        !U-Component
        delmax=max(u(nck)-u(i),0.0)
        delmin=min(u(nck)-u(i),0.0)
        delface=rx(k,i)*dux(i)+ry(k,i)*duy(i)   
        if(delface>delmax)then
          r=delmax/delface !Double positive
          ulim=min(ulim,func(r)) !Limiter
        elseif(delface<delmin)then
          r=delmin/delface !Double negative
          ulim=min(ulim,func(r)) !Limiter
        endif
        !V-component
        delmax=max(v(nck)-v(i),0.0)
        delmin=min(v(nck)-v(i),0.0)
        delface=rx(k,i)*dvx(i)+ry(k,i)*dvy(i)
        if(delface>delmax)then
          r=delmax/delface !Double positive
          vlim=min(vlim,func(r)) !Limiter
        elseif(delface<delmin)then
          r=delmin/delface !Double negative
          vlim=min(vlim,func(r)) !Limiter
        endif
      enddo
      dux(i)=dux(i)*ulim
      duy(i)=duy(i)*ulim
      dvx(i)=dvx(i)*vlim
      dvy(i)=dvy(i)*vlim
    enddo
!$OMP END DO
!$OMP END PARALLEL

    return    
    end subroutine lcdlimvec

!***********************************************************************    
    subroutine limitslope(funslopelim,phi,dphix,dphiy)
!   Limits the slopes dphix and dphiy (gradients) for variable phi
!   Written by Alex Sanchez, USACE-CHL
!***********************************************************************    
    use size_def
    use geo_def, only: idcellsimple,cell2cell,idcelljoint,ncface,idirface
    use der_def, only: rgx,rgy
    use prec_def
    implicit none
    !Input/Output
    real(ikind),intent(in),dimension(ncellsD) :: phi
    real(ikind),intent(inout),dimension(ncellsD) :: dphix,dphiy
    !Internal Variables
    integer :: i,ii,idk,k,nck,ncn,nce,ncs,ncw
    real(ikind) :: r,phinb(nmaxfaces)
    !real(ikind) :: snb(nmaxfaces),s2
   
    interface
      function funslopelim(r) result(psi)
        use prec_def
        real(ikind), intent(in) :: r
        real(ikind) :: psi
      end function
    endinterface

!$OMP PARALLEL
!--- Limiters for Regular Cells ---
!For regular cells use standard slope limiters
!$OMP DO PRIVATE(i,ii,ncn,nce,ncs,ncw,r)
    do ii=1,ncellsimple
      i=idcellsimple(ii)
      ncn=cell2cell(1,i)
      nce=cell2cell(2,i)
      ncs=cell2cell(3,i)
      ncw=cell2cell(4,i)
      if(abs(phi(i)-phi(ncw))>1.0e-10)then
        r=(phi(nce)-phi(i))/(phi(i)-phi(ncw))*rgx(i)
        dphix(i)=dphix(i)*funslopelim(r)
      endif
      if(abs(phi(i)-phi(ncs))>1.0e-10)then
        r=(phi(ncn)-phi(i))/(phi(i)-phi(ncs))*rgy(i)
        dphiy(i)=dphiy(i)*funslopelim(r)      
      endif
    enddo
!$OMP END DO NOWAIT
!$OMP DO PRIVATE(i,ii,k,idk,nck,r,phinb) !snb
    do ii=1,ncelljoint
      i=idcelljoint(ii)
      !Use most different neighbor compared to center value
      phinb=phi(i)
      do k=1,ncface(i)
        nck=cell2cell(k,i)
        idk=idirface(k,i)
        if(abs(phi(nck)-phi(i))>abs(phinb(idk)-phi(i)))then        
          phinb(idk)=phi(nck)
        endif  
      enddo
      !!Use average neighbor for each direction (allows slight overshoot)
      !phinb=0.0; snb=0.0
      !do k=1,ncface(i)
      !  nck=cell2cell(k,i)
      !  idk=idirface(k,i)
      !  phinb(idk)=phi(nck)
      !  snb(idk)=snb(idk)+1.0
      !enddo
      !phinb=phinb/snb
      if(abs(phi(i)-phinb(4))>1.0e-10)then
        r=(phinb(2)-phi(i))/(phi(i)-phinb(4))*rgx(i)
        dphix(i)=dphix(i)*funslopelim(r)
      endif
      if(abs(phi(i)-phinb(3))>1.0e-10)then
        r=(phinb(1)-phi(i))/(phi(i)-phinb(3))*rgy(i)
        dphiy(i)=dphiy(i)*funslopelim(r)  
      endif
    enddo
!$OMP END DO
!$OMP END PARALLEL

    return
    end subroutine limitslope

!***********************************************************************    
    subroutine limitslopevec(funslopelim,u,v,dux,duy,dvx,dvy)
! Limits the slopes dphix and dphiy (gradients) for variable phi
!
! Author: Alex Sanchez, USACE-CHL
!***********************************************************************    
#include "CMS_cpp.h"
    use size_def
    use geo_def, only: idcellsimple,cell2cell,idcelljoint,ncface,idirface
    use der_def, only: rgx,rgy
    use diag_def
    use diag_lib
    use prec_def
    implicit none
    !Input/Output
    real(ikind),intent(in),dimension(ncellsD) :: u,v
    real(ikind),intent(inout),dimension(ncellsD) ::dux,duy,dvx,dvy
    !Internal variables
    integer :: i,ii,idk,k,nck,ncn,nce,ncs,ncw
    real(ikind) :: r,unb(nmaxfaces),vnb(nmaxfaces)
    !real(ikind) :: snb(nmaxfaces),s2

    interface
      function funslopelim(r) result(psi)
        use prec_def
        real(ikind), intent(in) :: r
        real(ikind) :: psi
      end function
    endinterface

!$OMP PARALLEL
!--- Limiters for Regular Cells ---
!For regular cells use standard slope limiters
!$OMP DO PRIVATE(i,ii,ncn,nce,ncs,ncw,r)
    do ii=1,ncellsimple
      i=idcellsimple(ii)
      ncn=cell2cell(1,i)
      nce=cell2cell(2,i)
      ncs=cell2cell(3,i)
      ncw=cell2cell(4,i)
      if(abs(u(i)-u(ncw))>1.0e-10)then
        r=(u(nce)-u(i))/(u(i)-u(ncw))*rgx(i)
        dux(i)=dux(i)*funslopelim(r)
      endif
      if(abs(u(i)-u(ncs))>1.0e-10)then
        r=(u(ncn)-u(i))/(u(i)-u(ncs))*rgy(i)
        duy(i)=duy(i)*funslopelim(r)
      endif
      if(abs(v(i)-v(ncw))>1.0e-10)then
        r=(v(nce)-v(i))/(v(i)-v(ncw))*rgx(i)
        dvx(i)=dvx(i)*funslopelim(r) !Bug fix, replaced dux with dvx
      endif
      if(abs(v(i)-v(ncs))>1.0e-10)then   
        r=(v(ncn)-v(i))/(v(i)-v(ncs))*rgy(i)
        dvy(i)=dvy(i)*funslopelim(r) !Bug fix, replaced duy with dvy
      endif
#ifdef DIAG_MODE
      if(isnankind(dux(i)*duy(i)*dvx(i)*dvy(i)))then
        write(msg2,*) '  dux(i), duy(i) =   ',dux(i), duy(i)
        write(msg3,*) '  dvx(i), dvy(i) =   ',dvx(i), dvy(i)
        call diag_print_error('Problem applying limiter to vectors',&
            msg2,msg3)
      endif
#endif
    enddo
!$OMP END DO NOWAIT

!$OMP DO PRIVATE(i,ii,k,idk,nck,r,unb,vnb) !snb
    do ii=1,ncelljoint
      i=idcelljoint(ii)
      !Use most different neighbor compared to center value
      unb(:)=u(i); vnb(:)=v(i)
      do k=1,ncface(i)
        nck=cell2cell(k,i)
        idk=idirface(k,i)
        if(abs(u(nck)-u(i))>abs(unb(idk)-u(i)))then        
          unb(idk)=u(nck)
        endif
        if(abs(v(nck)-v(i))>abs(vnb(idk)-v(i)))then        
          vnb(idk)=v(nck)
        endif  
      enddo
      !!Use average neighbor for each direction (allows slight overshoot)
      !phinb=0.0; snb=0.0
      !do k=1,ncface(i)
      !  nck=cell2cell(k,i)
      !  idk=idirface(k,i)
      !  phinb(idk)=phi(nck)
      !  snb(idk)=snb(idk)+1.0
      !enddo
      !phinb=phinb/snb
      if(abs(u(i)-unb(4))>1.0e-10)then
        r=(unb(2)-u(i))/(u(i)-unb(4))*rgx(i)
        dux(i)=dux(i)*funslopelim(r)
      endif
      if(abs(u(i)-unb(3))>1.0e-10)then
        r=(unb(1)-u(i))/(u(i)-unb(3))*rgy(i)
        duy(i)=duy(i)*funslopelim(r)
      endif
      if(abs(v(i)-vnb(4))>1.0e-10)then
        r=(vnb(2)-v(i))/(v(i)-vnb(4))*rgx(i)
        dvx(i)=dvx(i)*funslopelim(r)
      endif
      if(abs(v(i)-vnb(3))>1.0e-10)then
        r=(vnb(1)-v(i))/(v(i)-vnb(3))*rgy(i)
        dvy(i)=dvy(i)*funslopelim(r)  
      endif
    enddo
!$OMP END DO
!$OMP END PARALLEL

    return
    end subroutine limitslopevec
        
end module der_lib
