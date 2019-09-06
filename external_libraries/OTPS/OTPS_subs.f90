!=====================================================================
module otps_subs
! AUTHORS:
!  Gary Egbert & Lana Erofeeva
!  College of Atmospheric and Oceanic Sciences
!  104 COAS Admin. Bldg.
!  Oregon State University
!  Corvallis, OR 97331-5503
!  
!  E-mail:  egbert@coas.oregonstate.edu                                      
!  Fax:     (541) 737-2064
!  Ph.:     (541) 737-2947                                        
!  http://volkov.oce.orst.edu/tides/
!
! COPYRIGHT: OREGON STATE UNIVERSITY, 2012
! (see the file COPYRIGHT for lisence agreement)
!=====================================================================

! Set of subroutines used by extract_HC and predict_tide
contains
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine rd_mod_body(modname,zuv,cind,ncon,z,n,m,mask)
      implicit none
      character*80 modname,hname,uname,gname,fname,rmCom
      character*1 zuv
      integer ncon,n,m
      complex z(ncon,n,m)
      complex, allocatable:: uv(:,:,:)
      integer mask(n,m),cind(ncon)
      integer ic,k,i,j,funit

      uname=' '
      gname=' '
      funit=111
      open(unit=funit,file=modname,status='old')
      read(funit,'(a80)')hname
      hname=rmCom(hname)
      read(funit,'(a80)',end=1)uname
      uname=rmCom(uname)
      read(funit,'(a80)',end=1)gname
      gname=rmCom(gname)
1     close(funit)
      fname=hname
      if(zuv.ne.'z')then
       fname=uname
      endif
      allocate(uv(2,n,m))
      open(unit=funit,file=trim(fname),status='old',form='unformatted') 
      mask=0
      do ic=1,ncon
       read(funit) ! pass header
       do k=1,cind(ic)
         if(zuv.eq.'z')then
          read(funit)((uv(1,i,j),i=1,n),j=1,m)
         else
          read(funit)uv
         endif
       enddo
       if(zuv.eq.'u'.or.zuv.eq.'U'.or.zuv.eq.'z')then
         z(ic,:,:)=uv(1,:,:)
       elseif(zuv.eq.'v'.or.zuv.eq.'V')then
         z(ic,:,:)=uv(2,:,:)
       endif
       rewind(funit)
       where(z(ic,:,:).ne.0)mask=1
      enddo
      close(funit)
      deallocate(uv)
      return
      endsubroutine rd_mod_body
      
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine rd_mod_body1(fname,zuv,cind,ncon,z,n,m,mask)
      implicit none
      character*80 fname
      character*1 zuv
      integer ncon,n,m
      complex z(ncon,n,m)
      complex, allocatable:: uv(:,:,:)
      integer mask(n,m),cind(ncon)
      integer ic,k,i,j,funit

      allocate(uv(2,n,m))
      funit=111
      open(unit=funit,file=trim(fname),status='old',form='unformatted',err=1) 
      mask=0
      do ic=1,ncon
       read(funit) ! pass header
       do k=1,cind(ic)
         if(zuv.eq.'z')then
          read(funit)((uv(1,i,j),i=1,n),j=1,m)
         else
          read(funit)uv
         endif
       enddo
       if(zuv.eq.'u'.or.zuv.eq.'U'.or.zuv.eq.'z')then
         z(ic,:,:)=uv(1,:,:)
       elseif(zuv.eq.'v'.or.zuv.eq.'V')then
         z(ic,:,:)=uv(2,:,:)
       endif
       rewind(funit)
       where(z(ic,:,:).ne.0)mask=1
      enddo
      close(funit)
      deallocate(uv)
      return
1     write(*,*)'File ',trim(fname),' does not EXIST'
      stop
      endsubroutine rd_mod_body1

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! define constituent indices in a model for given set of constituent
! constituent names are assumed to be all in LOWER case
      subroutine def_con_ind(c_id,ncon,c_id_mod,nc,cind)
      use otps_constit
      implicit none
      character*4 c_id(ncmx),c_id_mod(ncmx)
      integer nc,ncon,cind(ncon),ic1,ic2

      cind=0     
      do ic1=1,ncon
       do ic2=1,nc
        if(c_id(ic1).eq.c_id_mod(ic2))cind(ic1)=ic2
       enddo
      enddo

      return
      endsubroutine def_con_ind
      
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine rd_inp(modname,lltname,zuv,c_id,ncon,&
                       APRI,Geo,outname,interp)
      use otps_constit
      implicit none      
      character*80 modname,lltname,outname,tmp,rmCom
      character*1 zuv
      logical APRI,Geo,interp 
      character*4 c_id(ncmx)
      integer k,ic,ncon

      interp=.false.
      read(*,'(a80)',err=1)modname
      modname=rmCom(modname)
      if(trim(modname).eq.'')modname='model.list'
      read(*,'(a80)',err=1)lltname
      lltname=rmCom(lltname)
      read(*,'(a1)',err=1)zuv
      read(*,'(a80)',err=1)tmp
      tmp=rmCom(tmp)
      ic=0
      k=1
      do while(k.gt.0)
       ncon=ic
       k=index(tmp,',')
       ic=ic+1
       c_id(ic)=tmp(1:max(1,k-1))
       tmp=tmp(k+1:80)
      enddo
      if(trim(tmp).ne.'')ncon=ic
      c_id(ic)=tmp(1:80)
      read(*,'(a80)')tmp
      tmp=rmCom(tmp)
      APRI=.false.
      if(trim(tmp).eq.'AP')APRI=.true.
      geo=.false.
      read(*,'(a80)',err=1)tmp ! geo/oce
      tmp=rmCom(tmp)
      if(tmp(1:3).eq.'geo')geo=.true.
      read(*,*,err=1)k
      if(k.eq.1)interp=.true.
      read(*,'(a80)',err=1)outname
      outname=rmCom(outname)
      return
1     write(*,*)'Input file of WRONG FORMAT: see setup.inp'
      write(*,*)'for right format. Recovering setup.inp ... done'
      open(unit=111,file='dumsetup',status='old')
      open(unit=112,file='setup.inp',status='unknown')
3     read (111,'(a80)',end=2)tmp
      write(112,'(a)')trim(tmp)
      go to 3
2     close(111)
      close(112)      
      write(*,*)'Please edit ''setup.inp'' for your settings'
      write(*,*)'and re-run extract_HC or predict_tide'
      stop
      endsubroutine rd_inp
                       
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine rd_mod_header(modname,zuv,n,m,th_lim,ph_lim,nc,&
                              c_id,xy_ll_sub)
      use otps_constit
      implicit none
      real*4 th_lim(2),ph_lim(2)
      character*4 c_id(ncmx)
      character*80 modname
      character*80 hname,uname,gname,fname,rmCom
      character*80 xy_ll_sub
      character*1 zuv
      integer*4 nc,n,m,funit 
      
      uname=' '
      gname=' ' 
      xy_ll_sub=' '
      funit=111
      open(unit=funit,file=modname,status='old',err=1)
      read(funit,'(a80)',err=1)hname ! ALWAYS there
      hname=rmCom(hname)
      read(funit,'(a80)',end=5)uname
      uname=rmCom(uname)
      read(funit,'(a80)',end=5)gname
      read(funit,'(a80)',end=5)fname
      xy_ll_sub=rmCom(fname)
5     close(funit)
      fname=hname
      if(zuv.ne.'z')fname=uname
      if(fname.eq.' ')then
       write(*,*)'No name for UV file is given in '
       write(*,*) trim(modname)
       stop
      endif
      if((zuv.eq.'u'.or.zuv.eq.'v').and.gname.eq.' ')then
       write(*,*)'No name for bathymetry grid file is given in '
       write(*,*) trim(modname)
       stop
      endif

      open(unit=funit,file=trim(fname),status='old',form='unformatted',err=2)
      write(*,*)trim(fname)
      read(funit)n,m,nc,th_lim,ph_lim,c_id(1:nc)
      write(*,*)n,m,nc
      write(*,*)th_lim,ph_lim
      write(*,*)c_id(1:nc)
      close(funit)
      return
1     write(*,*)'File ',trim(modname),' does NOT exist or empty:'
      write(*,*)'Check line 1 of setup.inp'
      stop
2     write(*,*)'File ',trim(fname),' does NOT exist:'
      write(*,*)'Check ',trim(modname)
      stop
      endsubroutine rd_mod_header
                              
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine rd_mod_header1(fname,n,m,nc,th_lim,ph_lim,c_id)
      use otps_constit
      implicit none
      real th_lim(2),ph_lim(2)
      character*4 c_id(ncmx)
      character*80 fname
      integer nc,n,m,funit 

      funit=111
      open(unit=funit,file=trim(fname),status='old',form='unformatted',err=2) 
      read(funit)n,m,nc,th_lim,ph_lim,c_id(1:nc)
      close(funit)
      return
2     write(*,*)'File ',trim(fname),' does NOT exist:'
      stop
      endsubroutine rd_mod_header1

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      character*80 function rmCom(str) ! removes everything after '!' 
                                       ! in a string 
      character*80 str
      integer k 
      k=index(str,'!')
      rmCom=str
      if(k>0)rmCom=str(1:k-1)
      return
      endfunction rmCom

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine interp(uv,nt,n,m,mz,th_lim,ph_lim,&
                         xlat,xlon,uv1,ierr,zuv)
!     Interpolates complex array uv(nt,n,m) at point xlat,xlon
!     th_lim and ph_lim give latitude and longitude limits of grid
        implicit none
        integer ierr,n,m,mz(n,m),k,nt
        complex uv(nt,n,m)
        complex uv1(nt)
        real ww(0:1,0:1),dlat,dlon
        real th_lim(2),ph_lim(2),xlon,xlat
	    real dtheta,dphi,xlonc
        integer iw(0:1),jw(0:1)
        character*1 zuv

        ierr = 0
        dtheta = (th_lim(2)-th_lim(1))/m
        dphi = (ph_lim(2)-ph_lim(1))/n
! check on xlon longitude convention
        xlonc=xlon
! this SHOULD NOT BE DONE for streographic grids! (in km)
! I put uv1(1)=1 in calling routine to fix this
        if(uv1(1).ne.-1)then
         if(xlonc.lt.ph_lim(1))xlonc=xlonc+360.
         if(xlonc.gt.ph_lim(2))xlonc=xlonc-360.
         uv1(1)=0.
        endif
        if(xlonc.lt.ph_lim(1).or.xlonc.gt.ph_lim(2))then
          ierr=-1
          return
        endif
        if(xlat.lt.th_lim(1).or.xlat.gt.th_lim(2))then
          ierr=-1
          return
        endif
        dlat=xlat
        dlon=xlonc

        call BSI_weights(zuv,dlat,dlon,th_lim,ph_lim,&
                         dphi,dtheta,mz,n,m,ww,iw,jw)
        if(sum(ww).lt.0.01) then
          ierr = -2
          do k = 1,nt
            uv1(k) = (0.0,0.0)
          enddo
        else
         ierr = 0
         do k = 1,nt
            uv1(k) = uv(k,iw(0),jw(0))*ww(0,0)+ &
                     uv(k,iw(1),jw(0))*ww(1,0)+ &
                     uv(k,iw(0),jw(1))*ww(0,1)+ &
                     uv(k,iw(1),jw(1))*ww(1,1)
         enddo   ! k
      endif

      return
      endsubroutine interp

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine BSI_weights(zuv,theta,phi,theta_lim,phi_lim,&
                               dx,dy,mask,n,m,ww,iw,jw)
! bilinear spline interpolation weights
! for delta-forcing
    implicit none
    character*1 zuv               ! C-grid zuv: 'u','v','z' ('U','V','Z')
    real theta, phi              ! point coordinates
    real theta_lim(2),phi_lim(2)   ! grid limits
    real dx,dy                   ! step x, step y
    integer n,m,mask(n,m)          ! grid dimensions and mask
    real xi,xj,x,y
    integer i0,j0,i1,j1 
    real ww(0:1,0:1)
    integer iw(0:1),jw(0:1)
    integer ipshft,sm

    real w00,w01,w10,w11,wtot

    if(zuv.eq.'U')zuv='u'
    if(zuv.eq.'V')zuv='v'
    if(zuv.eq.'Z')zuv='z'
    if(zuv.ne.'u'.and.zuv.ne.'v'.and.zuv.ne.'z')then
      write(*,*)'Wrong call: BSI_weights, no ',zuv,' nodes'
      stop
    endif

    if(zuv.eq.'z')then
!      find grid coordinates for lower left hand z node
!      among the square (of z nodes which surround the obs point)
!
!                 |                     
!                \|/
!     -> -z---u---z---u---z--
!         |       |       |     x is the location given by theta,phi
!         v       v       v     the four surrounding elevation
!         | x     |       |     nodes are used for the bilinear
!     -> -z---u---z---u---z-    spline interpolation
!         |      /|\      |
!                 |
      xi = (phi-phi_lim(1))/dx+.5
      xj = (theta-theta_lim(1))/dy+.5
    elseif(zuv.eq.'u')then
!       interior point: measurement of current vector
!                       direction is given as unit vector (th,ph)
!                       (th,ph)=(1.,0.)->u
!                       (th,ph)=(0.,1.)->v 

!      find grid coordinates for lower left hand z node
!      among the square (of z nodes which surround the obs point)
!               
!            
!             u---z---u--
!             |       |         x is the location given by theta,phi
!             |   v---|---v     the 8 surrounding u,v 
!             |   | x |   |     nodes are used for the bilinear
!             u---z---u---z-    spline interpolation
!                 |       |     Current vector direction is given
!                 v-------v     with a unit vector (th,ph):
!                               (1.,0)->EW, (0.,1.)->NS
      xi = (phi-phi_lim(1))/dx+1.0
      xj = (theta-theta_lim(1))/dy+.5
    elseif(zuv.eq.'v')then
      xi = (phi-phi_lim(1))/dx+.5
      xj = (theta-theta_lim(1))/dy+1.0
    endif
    if(xi.lt.1.) xi = n+xi

    i0 = int(xi)
    x = xi-i0
    j0 = int(xj)
    y = xj-j0
!!!      check to see if calculated indices are in range
    if((i0.gt.n).or.(i0.lt.1).or.(j0.gt.m).or.(j0.lt.1) ) then
!            write(*,*)'Location is out of grid range:',theta,phi
!            write(*,*)'STOP in BSI_weights'
!            write(*,*)i0,j0
!            stop
      ww=0.
      return
    endif
!
!        compute weights for bilinear spline interpolation; only
!        use ocean nodes (for which mask mask is = 1)
    j1 = ipshft(j0,1,m)
    i1 = ipshft(i0,1,n)
    ww=0.
    sm=mask(i0,j0)+mask(i0,j1)+mask(i1,j0)+mask(i1,j1)
    if(sm.gt.0)then
      w00 = (1.-x)*(1.-y)*mask(i0,j0)
      w01 = (1.-x)*y*mask(i0,j1)
      w10 = x*(1.-y)*mask(i1,j0)
      w11 = x*y*mask(i1,j1)
      wtot = w00+w01+w10+w11
      if(wtot.eq.0)then
        ww=0
        return 
      endif
!
      ww(0,0) = w00/wtot
      ww(0,1) = w01/wtot
      ww(1,0) = w10/wtot
      ww(1,1) = w11/wtot
    endif
!
    iw(0)=i0
    iw(1)=i1
    jw(0)=j0
    jw(1)=j1
! 
    return
    endsubroutine BSI_weights

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    function ipshft(i,ish,n)
!      periodic shift maps i to i+ish  , mod n;
!        (always between 1 and n;  never 0 )
    integer i,ipshft,n,ish
    ipshft = mod(i+ish+n-1,n)+1
    return
    endfunction ipshft
      
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine xy_ll_N(x,y,lon,lat)
! translate x,y to lat,lon
! Uniform grid (in km) is centered in lat=90 (North Pole)
! Y oriented along lon=0
! Lana Erofeeva, May 2002      
    implicit none
    real :: y,x,lat,lon
    real, parameter :: pi=3.14159265358979

    lat=90.-sqrt(x**2+y**2)/111.7
    lon=atan2(y,x)*180./pi
    if(lon.lt.0.)lon=lon+360.
    if(x.eq.0.and.y.eq.0.)lon=0.
      
    return
    endsubroutine xy_ll_N

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine ll_xy_N(lon,lat,x,y)
! translate lat,lon to x,y (km)
! Uniform grid (in km) is centered in lat=90 (North Pole)
! Y oriented along lon=0
! Lana Erofeeva, May 2002      
    implicit none
    real y,x,lat,lon
    real, parameter :: pi=3.14159265358979
 
    x=(90.-lat)*111.7*cos(lon/180.*pi)
    y=(90.-lat)*111.7*sin(lon/180.*pi)
      
    return
    endsubroutine ll_xy_N

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine xy_ll_S(x,y,lon,lat)
! translate x,y to lat,lon Antarctic
! Uniform grid (in km) is centered in lat=-90 (South Pole)
! Y oriented along lon=0
! Lana Erofeeva, May 2002      
      implicit none
      real y,x,lat,lon
      real, parameter :: pi=3.14159265358979

      lat=-90.+sqrt(x**2+y**2)/111.7
      lon=-atan2(y,x)*180./pi+90.
      if(lon.gt.360.)lon=lon-360.
      if(lon.lt.  0.)lon=lon+360.
      
      return
      endsubroutine xy_ll_S

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine ll_xy_S(lon,lat,x,y)
! translate lat,lon to x,y (km) Antarctic
! Uniform grid (in km) is centered in lat=-90 (South Pole)
! Y oriented along lon=0
! Lana Erofeeva, May 2002      
      real :: y,x,lat,lon
      real, parameter :: pi=3.14159265358979

      x=-(90.+lat)*111.7*cos((90+lon)/180.*pi)
      y= (90.+lat)*111.7*sin((90+lon)/180.*pi)
      
      return
      endsubroutine ll_xy_S
      
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! make_a:: computes A matrix elements for one data point 
    subroutine make_a(interp,ind,nc,time,pu,pf,w,A,l_sal) 
    use otps_constit
    implicit none    
    integer, parameter:: ncon = 17 ! for case of using weights.h
    logical interp,l_sal
    real omega(ncmx),phase(ncmx)
	integer ind(nc)
	real w(ncon,8)
	real*8 pu(ncmx),pf(ncmx)
    real*8 time
    complex A(nc),c(ncmx)
    integer i,j,nc
    
! If l_sal=.true. - NO solid Earth correction is applied
! USING beta_SE coefficients
    if(.not.interp)then
      do j=1,ncmx
        omega(j)=omega_d(j)
        phase(j)=phase_mkB(j)
      enddo
  	  do j=1,nc
        i=ind(j)
        if(i.ne.0)then
 	      c(j) = cmplx( pf(i)*cos(omega(i)*time+phase(i)+pu(i)),&
                        pf(i)*sin(omega(i)*time+phase(i)+pu(i)))
        endif
	  enddo
!  remove solid Earth tide
      if(.not.l_sal)then
        do j=1,nc
          A(j)=0. 
          if(ind(j).ne.0)A(j)=c(j)*beta_SE(ind(j))
         enddo
      else
        do j=1,nc 
          A(j)=c(j)*1.
        enddo
      endif
    else
! this is the case when w from weights.h is used -> see
! comments to ../include/weights.h
!----------------------------------------------------------c        
      omega(1:ncon)=omega_d(1:ncon)
      phase(1:ncon)=phase_mkB(1:ncon)
  	  do i=1,ncon 
 	    c(i) = cmplx( pf(i)*cos(omega(i)*time+phase(i)+pu(i)),&
                      pf(i)*sin(omega(i)*time+phase(i)+pu(i)))
	  enddo
      A=cmplx(0,0)

! ind(j)=0 means the constituent is excluded
      do i=1,ncon
        do j=1,nc 
          if(ind(j).ne.0)A(j)=A(j)+c(i)*beta_SE(i)*w(i,ind(j))
        enddo
      enddo
!---------------------------------------------------------c

    endif      

    return
    endsubroutine make_a

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
    subroutine mkw(interp,ind,nc,wr)
    !include 'weights.h'
    use OTPS_weights
    real wr(17,8)
    logical interp
    integer j,nc,ind(nc)
    
    wr=w
    if(interp)return

    do j=1,nc
      if(ind(j).ne.0)wr(ind(j),:)=0.
    enddo
    
	return
    endsubroutine mkw

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	real function height(A,P,NC)
! returns height from model array of complex constituents
    implicit none
    integer i,nc
    complex A(nc),P(nc)
    real sum
    
    if(nc.eq.0)then
      height=0.
      return
    endif
    sum=0.
! height(i)=sum_of_real(A(i)*P(i))
    do i=1,nc
      sum=sum+real(P(i))*real(A(i))-aimag(P(i))*aimag(A(i))
    enddo
    height=sum
    return
    endfunction height

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! argUMENTS and ASTROL subroutines SUPPLIED by RICHARD RAY, March 1999
! attached to OTIS by Lana Erofeeva (subroutine nodal.f)
! NOTE - "no1" in constit.h corresponds to "M1" in arguments
	subroutine nodal(dtime,latitude,pu,pf)
    use otps_constit
    implicit none
    integer index(ncmx),i
    real*8 latitude,pu(ncmx),pf(ncmx)
    real*8 dtime,arg(53),f(53),u(53),pi
    
    !!!latitude = latitude !**** NOT used **********************
    
    data pi/3.14159265358979/
! index gives correspondence between constit.h and Richard's subroutines
! constit.h:       M2,S2,K1,O1,N2,P1,K2,q1,2N2,mu2,nu2,L2,t2,
!                  J1,M1(no1),OO1,rho1,Mf,Mm,SSA,M4,
!                  MS4,MN4,M6,M8,MK3,S6,2SM2,2MK3
!
	data index/30,35,19,12,27,17,37,10,25,26,28,33,34,&
                  23,14,24,11,5,3,2,45,46,44,50,0,42,51,40,0/
! f, u - same as pf, pu in old nodal.f; arg is NOT needed; 
! dtime - MJD
    call arguments(dtime,arg,f,u)
    pu=0
    pf=1
    do i=1,ncmx
! u is returned by "arguments" in degrees
      if(index(i).gt.0)then 
        pu(i)=u(index(i))*pi/180.
        pf(i)=f(index(i))
      endif
     !write(*,*)pu(i),pf(i)
    enddo
    return
    endsubroutine nodal

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine arguments( time1, arg, f, u )
!   Kernel routine for subroutine hat53.    Calculate tidal arguments.      
      implicit none
      real*8 time1, arg(*), f(*), u(*)
      real*8 shpn(4),s,h,p,omega,pp,hour,t1,t2
      real*8 tmp1,tmp2,temp1,temp2
      real*8 cosn,cos2n,sinn,sin2n,sin3n
      real*8 zero,one,two,three,four,five
      real*8 fiften,thirty,ninety
      real*8 pi, rad
      parameter       (pi=3.141592654d0, rad=pi/180.D0)
      parameter   (zero=0.d0, one=1.d0)
      parameter   (two=2.d0, three=3.d0, four=4.d0, five=5.d0)
      parameter   (fiften=15.d0, thirty=30.d0, ninety=90.d0)
      parameter   (pp=282.94) ! solar perigee at epoch 2000.
      equivalence (shpn(1),s),(shpn(2),h),(shpn(3),p),(shpn(4),omega)

!---  Determine equilibrium arguments ---------------------------
      call astrol( time1, shpn )
      hour = (time1 - int(time1))*24.d0
      t1 = fiften*hour
      t2 = thirty*hour
      arg( 1) = h - pp                                  ! Sa
      arg( 2) = two*h                                   ! Ssa
      arg( 3) = s - p                                   ! Mm
      arg( 4) = two*s - two*h                           ! MSf
      arg( 5) = two*s                                   ! Mf
      arg( 6) = three*s - p                             ! Mt
      arg( 7) = t1 - five*s + three*h + p - ninety      ! alpha1
      arg( 8) = t1 - four*s + h + two*p - ninety        ! 2Q1
      arg( 9) = t1 - four*s + three*h - ninety          ! sigma1
      arg(10) = t1 - three*s + h + p - ninety           ! q1
      arg(11) = t1 - three*s + three*h - p - ninety     ! rho1
      arg(12) = t1 - two*s + h - ninety                 ! o1
      arg(13) = t1 - two*s + three*h + ninety           ! tau1
      arg(14) = t1 - s + h + ninety                     ! M1
      arg(15) = t1 - s + three*h - p + ninety           ! chi1
      arg(16) = t1 - two*h + pp - ninety                ! pi1
      arg(17) = t1 - h - ninety                         ! p1
      arg(18) = t1 + ninety                             ! s1
      arg(19) = t1 + h + ninety                         ! k1
      arg(20) = t1 + two*h - pp + ninety                ! psi1
      arg(21) = t1 + three*h + ninety                   ! phi1
      arg(22) = t1 + s - h + p + ninety                 ! theta1
      arg(23) = t1 + s + h - p + ninety                 ! J1
      arg(24) = t1 + two*s + h + ninety                 ! OO1
      arg(25) = t2 - four*s + two*h + two*p             ! 2N2
      arg(26) = t2 - four*s + four*h                    ! mu2
      arg(27) = t2 - three*s + two*h + p                ! n2
      arg(28) = t2 - three*s + four*h - p               ! nu2
      arg(29) = t2 - two*s + h + pp                     ! M2a
      arg(30) = t2 - two*s + two*h                      ! M2
      arg(31) = t2 - two*s + three*h - pp               ! M2b
      arg(32) = t2 - s + p + 180.d0                     ! lambda2
      arg(33) = t2 - s + two*h - p + 180.d0             ! L2
      arg(34) = t2 - h + pp                             ! t2
      arg(35) = t2                                      ! S2
      arg(36) = t2 + h - pp + 180.d0                    ! R2
      arg(37) = t2 + two*h                              ! K2
      arg(38) = t2 + s + two*h - pp                     ! eta2
      arg(39) = t2 - five*s + 4.0*h + p                 ! MNS2
      arg(40) = t2 + two*s - two*h                      ! 2SM2
      arg(41) = 1.5*arg(30)                             ! M3
      arg(42) = arg(19) + arg(30)                       ! MK3
      arg(43) = three*t1                                ! S3
      arg(44) = arg(27) + arg(30)                       ! MN4
      arg(45) = two*arg(30)                             ! M4
      arg(46) = arg(30) + arg(35)                       ! MS4
      arg(47) = arg(30) + arg(37)                       ! MK4
      arg(48) = four*t1                                 ! S4
      arg(49) = five*t1                                 ! S5
      arg(50) = three*arg(30)                           ! M6
      arg(51) = three*t2                                ! S6
      arg(52) = 7.0*t1                                  ! S7
      arg(53) = four*t2                                 ! S8

!--- Determine nodal corrections f and u -----------------------------------
      sinn = sin(omega*rad)
      cosn = cos(omega*rad)
      sin2n = sin(two*omega*rad)
      cos2n = cos(two*omega*rad)
      sin3n = sin(three*omega*rad)
      f( 1) = one                                     ! Sa
      f( 2) = one                                     ! Ssa
      f( 3) = one - 0.130*cosn                        ! Mm
      f( 4) = one                                     ! MSf
      f( 5) = 1.043 + 0.414*cosn                      ! Mf
      f( 6) = sqrt((one+.203*cosn+.040*cos2n)**2 + &
                    (.203*sinn+.040*sin2n)**2)        ! Mt

      f( 7) = one                                     ! alpha1
      f( 8) = sqrt((1.+.188*cosn)**2+(.188*sinn)**2)  ! 2Q1
      f( 9) = f(8)                                    ! sigma1
      f(10) = f(8)                                    ! q1
      f(11) = f(8)                                    ! rho1
      f(12) = sqrt((1.0+0.189*cosn-0.0058*cos2n)**2 + &
                   (0.189*sinn-0.0058*sin2n)**2)      ! O1
      f(13) = one                                     ! tau1
!!!   tmp1  = 2.*cos(p*rad)+.4*cos((p-omega)*rad)
!!!   tmp2  = sin(p*rad)+.2*sin((p-omega)*rad)         ! Doodson's
      tmp1  = 1.36*cos(p*rad)+.267*cos((p-omega)*rad)  ! Ray's
      tmp2  = 0.64*sin(p*rad)+.135*sin((p-omega)*rad)
      f(14) = sqrt(tmp1**2 + tmp2**2)                 ! M1
      f(15) = sqrt((1.+.221*cosn)**2+(.221*sinn)**2)  ! chi1
      f(16) = one                                     ! pi1
      f(17) = one                                     ! P1
      f(18) = one                                     ! S1
      f(19) = sqrt((1.+.1158*cosn-.0029*cos2n)**2 + &
                   (.1554*sinn-.0029*sin2n)**2)       ! K1
      f(20) = one                                     ! psi1
      f(21) = one                                     ! phi1
      f(22) = one                                     ! theta1
      f(23) = sqrt((1.+.169*cosn)**2+(.227*sinn)**2)  ! J1
      f(24) = sqrt((1.0+0.640*cosn+0.134*cos2n)**2 + &
                   (0.640*sinn+0.134*sin2n)**2 )      ! OO1
      f(25) = sqrt((1.-.03731*cosn+.00052*cos2n)**2 + &
                   (.03731*sinn-.00052*sin2n)**2)     ! 2N2
      f(26) = f(25)                                   ! mu2
      f(27) = f(25)                                   ! N2
      f(28) = f(25)                                   ! nu2
      f(29) = one                                     ! M2a
      f(30) = f(25)                                   ! M2
      f(31) = one                                     ! M2b
      f(32) = one                                     ! lambda2
      temp1 = 1.-0.25*cos(two*p*rad) &
             -0.11*cos((two*p-omega)*rad)-0.04*cosn
      temp2 = 0.25*sin(two*p)+0.11*sin((two*p-omega)*rad) &
             + 0.04*sinn
      f(33) = sqrt(temp1**2 + temp2**2)               ! L2
      f(34) = one                                     ! t2
      f(35) = one                                     ! S2
      f(36) = one                                     ! R2
      f(37) = sqrt((1.+.2852*cosn+.0324*cos2n)**2 + &
                   (.3108*sinn+.0324*sin2n)**2)       ! K2
      f(38) = sqrt((1.+.436*cosn)**2+(.436*sinn)**2)  ! eta2
      f(39) = f(30)**2                                ! MNS2
      f(40) = f(30)                                   ! 2SM2
      f(41) = one   ! wrong                           ! M3
      f(42) = f(19)*f(30)                             ! MK3
      f(43) = one                                     ! S3
      f(44) = f(30)**2                                ! MN4
      f(45) = f(44)                                   ! M4
      f(46) = f(44)                                   ! MS4
      f(47) = f(30)*f(37)                             ! MK4
      f(48) = one                                     ! S4
      f(49) = one                                     ! S5
      f(50) = f(30)**3                                ! M6
      f(51) = one                                     ! S6
      f(52) = one                                     ! S7
      f(53) = one                                     ! S8

         u( 1) = zero                                    ! Sa
         u( 2) = zero                                    ! Ssa
         u( 3) = zero                                    ! Mm
         u( 4) = zero                                    ! MSf
         u( 5) = -23.7*sinn + 2.7*sin2n - 0.4*sin3n      ! Mf
         u( 6) = atan(-(.203*sinn+.040*sin2n)/ &
                       (one+.203*cosn+.040*cos2n))/rad   ! Mt
         u( 7) = zero                                    ! alpha1
         u( 8) = atan(.189*sinn/(1.+.189*cosn))/rad      ! 2Q1
         u( 9) = u(8)                                    ! sigma1
         u(10) = u(8)                                    ! q1
         u(11) = u(8)                                    ! rho1
         u(12) = 10.8*sinn - 1.3*sin2n + 0.2*sin3n       ! O1
         u(13) = zero                                    ! tau1
         u(14) = atan2(tmp2,tmp1)/rad                    ! M1
         u(15) = atan(-.221*sinn/(1.+.221*cosn))/rad     ! chi1
         u(16) = zero                                    ! pi1
         u(17) = zero                                    ! P1
         u(18) = zero                                    ! S1
         u(19) = atan((-.1554*sinn+.0029*sin2n)/ &
                      (1.+.1158*cosn-.0029*cos2n))/rad   ! K1
         u(20) = zero                                    ! psi1
         u(21) = zero                                    ! phi1
         u(22) = zero                                    ! theta1
         u(23) = atan(-.227*sinn/(1.+.169*cosn))/rad     ! J1
         u(24) = atan(-(.640*sinn+.134*sin2n)/ &
                     (1.+.640*cosn+.134*cos2n))/rad     ! OO1
         u(25) = atan((-.03731*sinn+.00052*sin2n)/ &
                     (1.-.03731*cosn+.00052*cos2n))/rad ! 2N2
         u(26) = u(25)                                   ! mu2
         u(27) = u(25)                                   ! N2
         u(28) = u(25)                                   ! nu2
         u(29) = zero                                    ! M2a
         u(30) = u(25)                                   ! M2
         u(31) = zero                                    ! M2b
         u(32) = zero                                    ! lambda2
         u(33) = atan(-temp2/temp1)/rad                  ! L2
         u(34) = zero                                    ! t2
         u(35) = zero                                    ! S2
         u(36) = zero                                    ! R2
         u(37) = atan(-(.3108*sinn+.0324*sin2n)/ &
                      (1.+.2852*cosn+.0324*cos2n))/rad   ! K2
         u(38) = atan(-.436*sinn/(1.+.436*cosn))/rad     ! eta2
         u(39) = u(30)*two                               ! MNS2
         u(40) = u(30)                                   ! 2SM2
         u(41) = 1.5d0*u(30)                             ! M3
         u(42) = u(30) + u(19)                           ! MK3
         u(43) = zero                                    ! S3
         u(44) = u(30)*two                               ! MN4
         u(45) = u(44)                                   ! M4
         u(46) = u(30)                                   ! MS4
         u(47) = u(30)+u(37)                             ! MK4
         u(48) = zero                                    ! S4
         u(49) = zero                                    ! S5
         u(50) = u(30)*three                             ! M6
         u(51) = zero                                    ! S6
         u(52) = zero                                    ! S7
         u(53) = zero                                    ! S8

      return
      endsubroutine arguments

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE ASTROL( time, SHPN )
!  Computes the basic astronomical mean longitudes  s, h, p, N.
!  Note N is not N', i.e. N is decreasing with time.
!  These formulae are for the period 1990 - 2010, and were derived
!  by David Cartwright (personal comm., Nov. 1990).
!  time is UTC in decimal MJD.
!  All longitudes returned in degrees.
!  R. D. Ray    Dec. 1990
!
!  Non-vectorized version.
!
!      IMPLICIT real*8 (A-H,O-Z)
      real*8 circle,shpn,t,time
      DIMENSION  SHPN(4)
      PARAMETER  (CIRCLE=360.0D0)

      T = time - 51544.4993D0
!
!     mean longitude of moon
!     ----------------------
      SHPN(1) = 218.3164D0 + 13.17639648D0 * T
!
!     mean longitude of sun
!     ---------------------
      SHPN(2) = 280.4661D0 +  0.98564736D0 * T
!
!     mean longitude of lunar perigee
!     -------------------------------
      SHPN(3) =  83.3535D0 +  0.11140353D0 * T
!
!     mean longitude of ascending lunar node
!     --------------------------------------
      SHPN(4) = 125.0445D0 -  0.05295377D0 * T

      SHPN(1) = MOD(SHPN(1),CIRCLE)
      SHPN(2) = MOD(SHPN(2),CIRCLE)
      SHPN(3) = MOD(SHPN(3),CIRCLE)
      SHPN(4) = MOD(SHPN(4),CIRCLE)

      IF (SHPN(1).LT.0.D0) SHPN(1) = SHPN(1) + CIRCLE
      IF (SHPN(2).LT.0.D0) SHPN(2) = SHPN(2) + CIRCLE
      IF (SHPN(3).LT.0.D0) SHPN(3) = SHPN(3) + CIRCLE
      IF (SHPN(4).LT.0.D0) SHPN(4) = SHPN(4) + CIRCLE
      RETURN
      ENDSUBROUTINE ASTROL

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE CALDAT(JULIAN,MM,ID,IYYY)                                                                                            
!   ROUTINE CONVERTS JULIAN DAY TO MONTH, DAY, & YEAR.                      
!   THIS CODE IS LIFTED FROM THE BOOK:                                      
!   W.H. PRESS ET AL., NUMERICAL RECIPES, CAMBRIDGE UNIV. PRESS, 1986.      
!   THE ONLY MODIFICATION IS THAT REAL ARITHMETIC IS DONE IN R*8.           
!                                                                           
!   TO CONVERT MODIFIED JULIAN DAY, CALL THIS ROUTINE WITH                 
!     JULIAN = MJD + 2400001     
      IMPLICIT NONE
      INTEGER :: JULIAN,MM,ID,IYYY                                                              
      INTEGER, PARAMETER :: IGREG=2299161 
      INTEGER :: JALPHA,JA,JB,JC,JD,JE
      
      IF (JULIAN.GE.IGREG) THEN                                            
         JALPHA=INT((DBLE(JULIAN-1867216)-0.25D0)/36524.25D0)              
         JA=JULIAN+1+JALPHA-INT(0.25D0*DBLE(JALPHA))                       
      ELSE                                                                 
         JA=JULIAN                                                         
      ENDIF                                                                
      JB=JA+1524                                                           
      JC=INT(6680.D0+(DBLE(JB-2439870)-122.1D0)/365.25D0)                  
      JD=365*JC+INT(0.25D0*DBLE(JC))                                       
      JE=INT(DBLE(JB-JD)/30.6001D0)                                        
      ID=JB-JD-INT(30.6001D0*DBLE(JE))                                     
      MM=JE-1                                                              
      IF (MM.GT.12) MM=MM-12                                               
      IYYY=JC-4715                                                         
      IF (MM.GT.2) IYYY=IYYY-1                                             
      IF (IYYY.LE.0) IYYY=IYYY-1    
      
      RETURN                                                               
      ENDSUBROUTINE CALDAT
      
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine date_mjd(mm,id,iyyy,mjd)
! Lana Erofeeva, Feb 2002
! Debugged for consistency with CALDAT for any dates mjd>=0
! converts date to mjd
! INPUT:  id - day, mm - month, iyyy - year
! OUTPUT: mjd>0 - modified julian days
! date>=11.17.1858 corresponds to mjd=0
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc       
    implicit none
    integer dpm(12),days,i,nleap,k
    integer mm,id,iyyy,mjd
    data dpm/31,28,31,30,31,30,31,31,30,31,30,31/
      
    mjd=0
    ! NO earlier dates
    if(iyyy.lt.1858)iyyy=1858
    if(iyyy.eq.1858.and.mm.gt.11)mm=11
    if(iyyy.eq.1858.and.mm.eq.11.and.id.gt.17)id=17
    !
    days=0 
    do i=1,mm-1
      days=days+dpm(i)
      if(i.eq.2.and.int(iyyy/4)*4.eq.iyyy)days=days+1         
    enddo
    days=days+id-321
    ! leap day correction
    do k=1900,iyyy,100
      if(iyyy.eq.k.and.mm.gt.2)days=days-1
    enddo
    do k=2000,iyyy,400
      if(iyyy.eq.k.and.mm.gt.2)days=days+1
    enddo
    ! EACH 4th year is leap year
    nleap=int((iyyy-1-1860)*0.25)            
    if(iyyy.gt.1860)nleap=nleap+1
    ! EXCEPT
    do k=1900,iyyy-1,100
      if(iyyy.gt.k)nleap=nleap-1
      if(iyyy.eq.k.and.mm.gt.2)days=days-1
    enddo
    ! BUT each in the row 2000:400:... IS LEAP year again
    do k=2000,iyyy-1,400
      if(iyyy.gt.k)nleap=nleap+1
      if(iyyy.eq.k.and.mm.gt.2)days=days+1
    enddo
    mjd=365*(iyyy-1858)+nleap+days   
      
    return
    endsubroutine date_mjd

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine ptide(z1,cid,ncon,ind,lat,time_mjd,ntime,interp,zpred)
    use otps_constit
    use otps_weights
    implicit none
    integer ncon,ntime,k,ind(ncon),ierr
    complex z1(ncon)
    real lat,zpred(ntime)
    real*8 time_mjd(ntime),time,dh
    character*4 cid(ncon)
    integer SecondsPerDay
    real height
    real ww(17,8)
    complex, allocatable:: A(:)
! nodal arguments shoud be real*8
    real*8 pu(ncmx),pf(ncmx),dlat
    logical interp
!
    if(interp)call mkw(interp,ind,ncon,ww)
    allocate(A(ncon))
	SecondsPerDay=86400
    dlat=lat
    ierr=0
    dh=0.0
    do k=1,ntime
	  call nodal(time_mjd(k),dlat,pu,pf)
! to use phase shifts from constit.h time should be in seconds
! relatively Jan 1 1992 (48622mjd)       
         time=(time_mjd(k)-dble(48622))*dble(SecondsPerDay)
! .true. means NO solid Earth correction will be applied in make_a
         call make_a(.false.,ind,ncon,time,pu,pf,ww,A,.true.)
         zpred(k)=height(A,z1,ncon)
         if(interp)call infer_minor(z1,cid,ncon,time_mjd(k),dh,ierr)
         if(ierr.eq.-1)then
          write(*,*)'Not enough constituents for inference of'
          write(*,*)'minor constituents: IGNORED'
          interp=.false.
         endif
         zpred(k)=zpred(k)+dh ! add minor constituents
    enddo
    deallocate(A)
    return
    endsubroutine ptide

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    character*10 function deblank(str)
    character*10 str
    integer k
    
1   k=index(str,' ')
    if(k.gt.0)then
    str(k:k)='0'
    go to 1
    endif
    deblank=str
    
    return
    endfunction deblank

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine def_cid(nc0,cid,ind)
    use otps_constit
    implicit none
    integer nc0,ic,jc,ind(nc0),k
    character*4 cid(ncmx)
!
    k=1
    do ic = 1,nc0
      do jc = 1,ncmx             
        if(cid(ic).eq.constid(jc)) go to 4
      enddo
      write(*,*) 'WARNING: Constituent ID ',cid(ic),' is not allowed'
      ind(k)=0
      k=k+1
4     continue
      ind(k)=jc
      k=k+1
    enddo
    
    return
    endsubroutine def_cid

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine infer_minor(zmaj,cid,ncon,time,dh,ierr)
! Lana Erofeeva, June 2004
! Based on Richard Ray code perth2
! Return correction for the 16 minor constituents
! listed in subroutine
!         infer minor tides
    implicit none
    integer ncon,i,j,ni
    complex zmaj(ncon)       ! HC for GIVEN consituents
    character*4 cid(ncon)    ! GIVEN constituents
    real*8 time              ! time, mjd
    real*8 dh                ! OUTPUT: correction at given time
                            ! for 16 minor constituents
    integer ierr             ! -1 if not enough constuents
                            !    for inference 
    complex*16 zmin(18)
    complex*16 z8(8)
    real*8 hour,t1,t2,shpn(4),s,h,p,omega
    real*8 sinn,cosn,sin2n,cos2n
    real*8 u(18),f(18),arg(18)

    real*8, parameter:: pi=3.141592654D0
    real*8, parameter:: rad=pi/180.D0
    real*8, parameter:: PP=282.8D0 
    EQUIVALENCE (SHPN(1),S),(SHPN(2),H),(SHPN(3),P),(SHPN(4),omega) 

    character*4 cid8(8)      ! in order to correspnd RR coefficients
    data cid8/'q1  ','o1  ','p1  ','k1  ',&
            'n2  ','m2  ','s2  ','k2  '/
! re-order to correspond to cid8
        ierr=0
        z8=0.
        ni=0
        do i=1,8
        do j=1,ncon
        if(cid(j).eq.cid8(i))then
            z8(i)=zmaj(j)
            if(i.ne.3.and.i.ne.8)ni=ni+1
        endif 
        enddo
        enddo

        if(ni.lt.6)then
        ierr=-1 ! not enough constituents for inference
        return 
        endif
!            
        zmin(1)  = 0.263 *z8(1) - 0.0252*z8(2)  !2Q1
        zmin(2) = 0.297 *z8(1) - 0.0264*z8(2)   !sigma1
        zmin(3) = 0.164 *z8(1) + 0.0048*z8(2)   !rho1 +
        zmin(4) = 0.0140*z8(2) + 0.0101*z8(4)   !M1
        zmin(5) = 0.0389*z8(2) + 0.0282*z8(4)   !M1
        zmin(6) = 0.0064*z8(2) + 0.0060*z8(4)   !chi1
        zmin(7) = 0.0030*z8(2) + 0.0171*z8(4)   !pi1
        zmin(8) =-0.0015*z8(2) + 0.0152*z8(4)   !phi1
        zmin(9) =-0.0065*z8(2) + 0.0155*z8(4)   !theta1
        zmin(10) =-0.0389*z8(2) + 0.0836*z8(4)  !J1 +
        zmin(11) =-0.0431*z8(2) + 0.0613*z8(4)  !OO1 +
        zmin(12) = 0.264 *z8(5) - 0.0253*z8(6)  !2N2 +
        zmin(13) = 0.298 *z8(5) - 0.0264*z8(6)  !mu2 +
        zmin(14) = 0.165 *z8(5) + 0.00487*z8(6) !nu2 +
        zmin(15) = 0.0040*z8(6) + 0.0074*z8(7)  !lambda2
        zmin(16) = 0.0131*z8(6) + 0.0326*z8(7)  !L2 +
        zmin(17) = 0.0033*z8(6) + 0.0082*z8(7)  !L2 +
        zmin(18) = 0.0585*z8(7)                 !t2 + 
!
        hour = (time - int(time))*24.D0
        t1 = 15.*hour
        t2 = 30.*hour
        call astrol( time, SHPN )
!
        arg(1) = t1 - 4.*S + H + 2.*P - 90.     ! 2Q1
        arg(2) = t1 - 4.*S + 3.*H - 90.         ! sigma1
        arg(3) = t1 - 3.*S + 3.*H - P - 90.     ! rho1
        arg(4) = t1 - S + H - P + 90.           ! M1
        arg(5) = t1 - S + H + P + 90.           ! M1
        arg(6) = t1 - S + 3.*H - P + 90.        ! chi1
        arg(7) = t1 - 2.*H + PP - 90.           ! pi1
        arg(8) = t1 + 3.*H + 90.                ! phi1
        arg(9) = t1 + S - H + P + 90.           ! theta1
        arg(10) = t1 + S + H - P + 90.          ! J1
        arg(11) = t1 + 2.*S + H + 90.           ! OO1
        arg(12) = t2 - 4.*S + 2.*H + 2.*P       ! 2N2
        arg(13) = t2 - 4.*S + 4.*H              ! mu2
        arg(14) = t2 - 3.*S + 4.*H - P          ! nu2
        arg(15) = t2 - S + P + 180.D0           ! lambda2
        arg(16) = t2 - S + 2.*H - P + 180.D0    ! L2
        arg(17) = t2 - S + 2.*H + P             ! L2
        arg(18) = t2 - H + PP                   ! t2
!
!     determine nodal corrections f and u
        sinn = SIN(omega*rad)
        cosn = COS(omega*rad)
        sin2n = SIN(2.*omega*rad)
        cos2n = COS(2.*omega*rad)
        F = 1.
        f(1) = sqrt((1.0 + 0.189*cosn - 0.0058*cos2n)**2 + &
                    (0.189*sinn - 0.0058*sin2n)**2)
        f(2) = f(1)
        f(3) = f(1)
        f(4) = sqrt((1.0 + 0.185*cosn)**2 + (0.185*sinn)**2)
        f(5) = sqrt((1.0 + 0.201*cosn)**2 + (0.201*sinn)**2)
        f(6) = sqrt((1.0 + 0.221*cosn)**2 + (0.221*sinn)**2)
        f(10) = sqrt((1.0 + 0.198*cosn)**2 + (0.198*sinn)**2)
        f(11) = sqrt((1.0 + 0.640*cosn + 0.134*cos2n)**2 + &
                    (0.640*sinn + 0.134*sin2n)**2 )
        f(12) = sqrt((1.0 - 0.0373*cosn)**2 + (0.0373*sinn)**2)
        f(13) = f(12)
        f(14) = f(12)
        f(16) = f(12)
        f(17) = sqrt((1.0 + 0.441*cosn)**2 + (0.441*sinn)**2)
!
        U = 0.
        u(1) = atan2(0.189*sinn - 0.0058*sin2n, &
                    1.0 + 0.189*cosn - 0.0058*sin2n)/rad
        u(2) = u(1)
        u(3) = u(1)
        u(4) = atan2( 0.185*sinn, 1.0 + 0.185*cosn)/rad
        u(5) = atan2(-0.201*sinn, 1.0 + 0.201*cosn)/rad
        u(6) = atan2(-0.221*sinn, 1.0 + 0.221*cosn)/rad
        u(10) = atan2(-0.198*sinn, 1.0 + 0.198*cosn)/rad
        u(11) = atan2(-0.640*sinn - 0.134*sin2n, &
                    1.0 + 0.640*cosn + 0.134*cos2n)/rad
        u(12) = atan2(-0.0373*sinn, 1.0 - 0.0373*cosn)/rad
        u(13) = u(12)
        u(14) = u(12)
        u(16) = u(12)
        u(17) = atan2(-0.441*sinn, 1.0 + 0.441*cosn)/rad

!     sum over all tides
!     ------------------
    dh = 0.D0
    do i=1,18
        dh = dh + dreal(zmin(i))*f(i)*cos((arg(i)+u(i))*rad)- &
                dimag(zmin(i))*f(i)*sin((arg(i)+u(i))*rad)
    enddo
    return
    endsubroutine infer_minor

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine rd_inp_JS(setup_file_unit,modname,zuv,c_id,ncon,&
                    APRI,Geo,interp)
!
    use otps_constit
    implicit none
    integer setup_file_unit
    character*80 modname,outname,tmp,rmCom
    character*80 lltname
    character*1 zuv
    logical APRI,Geo,interp 
    character*4 c_id(ncmx)
    integer k,ic,ncon
!
    interp=.false.
    read(setup_file_unit,'(a80)',err=1)modname
    modname=rmCom(modname)
    read(setup_file_unit,'(a80)',err=1)lltname
    lltname=rmCom(lltname)
    read(setup_file_unit,'(a1)',err=1)zuv
    read(setup_file_unit,'(a80)',err=1)tmp
    tmp=rmCom(tmp)
    write(*,*)tmp
    ic=0
    k=1
    do while(k.gt.0)
    ncon=ic
    k=index(tmp,',')
    ic=ic+1
    c_id(ic)=tmp(1:k-1)
    tmp=tmp(k+1:80)
    enddo
    if(trim(tmp).ne.'')ncon=ic
    c_id(ic)=tmp(1:80)
    read(setup_file_unit,'(a80)')tmp
    tmp=rmCom(tmp)
    APRI=.false.
    if(trim(tmp).eq.'AP')APRI=.true.
    geo=.false.
    read(setup_file_unit,'(a80)',err=1)tmp ! geo/oce
    tmp=rmCom(tmp)
    if(tmp(1:3).eq.'geo')geo=.true.
    read(setup_file_unit,*,err=1)k
    if(k.eq.1)interp=.true.
    read(setup_file_unit,'(a80)',err=1)outname
    outname=rmCom(outname)
    return
1   write(*,*)'Input file of WRONG FORMAT: see setup.inp'
    write(*,*)'for right format. Recovering setup.inp ... done'
    open(unit=111,file='dumsetup',status='old')
    open(unit=112,file='setup.inp',status='unknown')
3   read (111,'(a80)',end=2)tmp
    write(112,'(a)')trim(tmp)
    go to 3
2   close(111)
    close(112)      
    write(*,*)'Please edit ''setup.inp'' for your settings'
    write(*,*)'and re-run tide_js'
    stop
    endsubroutine rd_inp_JS
                       
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine mapxy(n,m,x,y,lon,lat,SLAT,SLON,HEMI)
! (X,Y)km->(lon,lat)deg
    implicit none
    integer*4 n,m
    real*8 x(n,m),y(n,m),lat(n,m),lon(n,m),SLAT,SLON,SLAT1
    real*8 CDR,E2,E,pi,RE,SGN,CM,SL
    real*8 a1,a2,a3,a4,a5
    real*8, allocatable:: RHO(:,:),T(:,:),CHI(:,:)
    character*1 HEMI
!
    allocate(RHO(n,m),T(n,m),CHI(n,m))
!
    CDR= 57.29577951
    E2 = 6.694379852e-3           ! Eccentricity squared
    E  = sqrt(E2)
    pi = 3.14159265358979
    RE = 6378.273
!                                                                         
!*************************************************************************   
    if(HEMI.eq.'S'.or.HEMI.eq.'s')then
        SGN=-1
    else
        SGN= 1
    endif
    SLAT1=abs(SLAT)
    SL  = SLAT1/CDR
    RHO = sqrt(X**2+Y**2)
    if(n.eq.1.and.m.eq.1.and.RHO(1,1).lt.0.1)then
! Don't calculate if ONE POINT is on the equator
        lat=90.*SGN
        lon=0.0
        return
    else
        CM=cos(SL)/sqrt(1.0-E2*(sin(SL)**2)) 
        T=tan((pi/4.0)-(SL/2))/((1.0-E*sin(SL))/ &
                    (1.0+E*sin(SL)))**(E/2.0)
        if(abs(SLAT1-90.).lt.1.e-5)then
            T=RHO*sqrt((1.+E)**(1+E)*(1-E)**(1-E))/2/RE
        else
            T=RHO*T/(RE*CM)
        endif
        a1 =  5*E2**2 / 24
        a2 =    E2**3 / 12
        a3 =  7*E2**2 / 48
        a4 = 29*E2**3 /240
        a5 =  7*E2**3 /120
        CHI= (pi/2)-2*atan(T);
        lat= CHI+((E2/2) + a1 + a2)*sin(2*CHI)+(a3 + a4)*sin(4*CHI)+ &
                                    a5*sin(6*CHI)
        lat=SGN*lat*CDR
        lon= -(atan2(-X,Y)*CDR)+SLON
    endif
    deallocate(RHO,T,CHI)
    return
    endsubroutine mapxy

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine mapll(n,m,lon,lat,x,y,SLAT,SLON,HEMI)
! (lon,lat)deg->(X,Y)km
    implicit none
    integer*4 n,m
    real*8 x(n,m),y(n,m),lat(n,m),lon(n,m),SLAT,SLON
    real*8 CDR,E2,E,pi,RE,SGN,RHO,SL,TC,MC
    character*1 HEMI
    real*8, allocatable:: T(:,:),lat1(:,:),lon1(:,:)
    allocate(T(n,m),lat1(n,m),lon1(n,m))
    CDR= 57.29577951
    E2 = 6.694379852e-3           ! Eccentricity squared
    E  = sqrt(E2)
    pi = 3.14159265358979
    RE = 6378.273
!                                                                         
!*************************************************************************   
    if(HEMI.eq.'S'.or.HEMI.eq.'s')then
        SGN=-1
    else
        SGN= 1
    endif
! This only works for Southern hemisphere lat/lon
!*************************************************************************
    if(abs(SLAT).eq.90)then
        RHO=2*RE/((1+E)**(1+E)*(1-E)**(1-E))**(E/2)
    else
        SL  = abs(SLAT)/CDR
        TC  = tan(pi/4-SL/2)/((1-E*sin(SL))/(1+E*sin(SL)))**(E/2)
        MC  = cos(SL)/sqrt(1-E2*(sin(SL)**2))
        RHO = RE*MC/TC;
    endif
    !lat1 = abs(lat)/CDR
        lat1= -lat/CDR
    T   = tan(pi/4-lat1/2)/((1-E*sin(lat1))/(1+E*sin(lat1)))**(E/2)
    lon1 =-(lon-SLON)/CDR
    x   =-RHO*T*sin(lon1)
    y   = RHO*T*cos(lon1)
    deallocate(T,lat1,lon1)
    return
    endsubroutine mapll

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine xy_ll_CATs(x,y,lon,lat)
    implicit none
    real :: y,x,lat,lon
    real*8 :: SLAT,SLON
    real*8, dimension(1,1) :: dx,dy,dlon,dlat
    character*1 :: HEMI
    
    SLAT=71.;SLON=-70.;HEMI='S'
    dx=x;dy=y;
    call mapxy(1,1,dx,dy,dlon,dlat,SLAT,SLON,HEMI)
    lon=dlon(1,1);lat=dlat(1,1)
    
    return
    endsubroutine xy_ll_CATs

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine ll_xy_CATs(lon,lat,x,y)
    implicit none
    real y,x,lat,lon
    real*8 :: SLAT,SLON
    real*8, dimension(1,1) :: dx,dy,dlon,dlat
    character*1 HEMI
    
    SLAT=71.;SLON=-70.;HEMI='S'
    dlon=lon;dlat=lat;
    call mapll(1,1,dlon,dlat,dx,dy,SLAT,SLON,HEMI)
    x=dx(1,1);y=dy(1,1);
    
    return
    endsubroutine ll_xy_CATs

endmodule otps_subs