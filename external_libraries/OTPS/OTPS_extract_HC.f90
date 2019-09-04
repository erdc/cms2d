!=====================================================================
module otps_extract_hc
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
    implicit none
    
contains    
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine extract_hc(modname,ncon,cname,ndat,lat,lon,amp,pha)
!   LANA, 2004 (based on old heval and uveval) 
!   reads OTIS standard binary complex model file
!   (elevations OR transports), reads a list of locations 
!   and outputs ASCII file with the complex amplitudes/amp,phases of
!   elevations/transports/currents interpolated at the locations
! 
    use otps_constit
    implicit none
    !Input/Output
    character*80,intent(in)    :: modname
    integer,     intent(inout) :: ncon
    character*6, intent(inout) :: cname(ncon)
    integer,     intent(in)    :: ndat
    real,        intent(inout) :: lat(ndat),lon(ndat)
    real,        intent(out)   :: amp(ncon,ndat),pha(ncon,ndat)
    !Internal variables
    complex, allocatable:: z(:,:,:),z1(:),dtmp(:,:)
    complex, allocatable:: zl(:,:,:),zl1(:)
    complex d1
    real, allocatable:: depth(:,:),x(:),y(:),lon0(:)
    real, allocatable:: phase(:)
    real latp,lonp,xt,yt
    real th_lim(2),ph_lim(2),lth_lim(2),lph_lim(2)
    integer, allocatable:: mask(:,:),cind(:),lcind(:),&
                            mz(:,:),mzl(:,:)
!
    character*4 :: c_id(ncmx),c_id_mod(ncmx),lc_id(ncmx)
    character*80 :: lltname,outname,gname,lname
    character*80 :: rmCom
    character*1 :: zuv,c1
    character*80 :: xy_ll_sub
    logical :: APRI,geo,interp_micon,ll_km
    integer :: nc,n,m,k,ierr,ierr1,ic,n0,m0
    integer :: ncl,nl,ml,nmod,ibl
    real,parameter :: pi = acos(-1.0)
    real,parameter :: rad2deg = 180.0/pi
!
    ll_km=.false.
    nmod=1
    ibl=0
    lname='DATA/load_file'
    !modname = 'tpxo7'
    lltname = 'lat_lon_time'    
    zuv = 'z'
    APRI = .true.
    outname = 'sample.out'
    interp_micon = .false.
      
    do ic=1,ncon
      c_id(ic)=cname(ic)
    enddo
    !call rd_inp(modname,lltname,zuv,c_id,ncon,APRI,geo,&
    !            outname,interp_micon)
    
!    if(trim(modname)=='model.list')then
!      nmod=0
!      open(unit=12,file='model.list',status='old',err=11)
!      write(*,*)'Models to extract HC from:'
!13    read(12,'(a80)',end=12)ctmp
!      write(*,*)trim(ctmp)
!      nmod=nmod+1
!      goto 13
!12    rewind(12)
!      read(12,'(a80)') modname
!    endif
    write(*,*)
    write(*,*)'Lat/Lon file:', trim(lltname)
    if(ncon>0)write(*,*)'Constituents to include: ',c_id(1:ncon)
    if(geo.and.zuv=='z')then
       write(*,*)'Extract GEOCENTRIC tide HC'
    else
       write(*,*)'Extract OCEAN tide HC'
    endif
!
    !do imod=1,nmod
    !write(*,*)
    !if(imod>1)read(12,'(a80)') modname
    call rd_mod_header(modname,zuv,n,m,th_lim,ph_lim,nc,c_id_mod,&
                          xy_ll_sub)
    write(*,*)'Model:        ',trim(modname)
    if(trim(xy_ll_sub)=='')then
      write(*,*)'Lat limits:   ',th_lim
      write(*,*)'Lon limits:   ',ph_lim
    else
      ll_km=.true.
      if(trim(xy_ll_sub)/='xy_ll_N'.and. &
         trim(xy_ll_sub)/='xy_ll_S'.and. &
         trim(xy_ll_sub)/='xy_ll_CATs')then
        write(*,*)'No converting function ', trim(xy_ll_sub),&
                  'in the OTPS'
        stop 
      endif
      if(trim(xy_ll_sub)=='xy_ll_N')then
        call xy_ll_N(ph_lim(1),th_lim(1),xt,yt)
        if(xt>180.0)xt=xt-360.0
        write(*,*)'Lat,Lon lower left corner:',yt,xt
        call xy_ll_N(ph_lim(1),th_lim(2),xt,yt)
        if(xt>180.0)xt=xt-360.0
        write(*,*)'Lat,Lon upper left corner:',yt,xt
        call xy_ll_N(ph_lim(2),th_lim(1),xt,yt)
        if(xt>180)xt=xt-360
        write(*,*)'Lat,Lon lower right corner:',yt,xt
        call xy_ll_N(ph_lim(2),th_lim(2),xt,yt)
        if(xt>180.0)xt=xt-360.0
        write(*,*)'Lat,Lon upper right corner:',yt,xt
      elseif(trim(xy_ll_sub)=='xy_ll_S')then
        call xy_ll_S(ph_lim(1),th_lim(1),xt,yt)
        if(xt>180.0)xt=xt-360.0
        write(*,*)'Lat,Lon lower left corner:',yt,xt
        call xy_ll_S(ph_lim(1),th_lim(2),xt,yt)
        if(xt>180.0)xt=xt-360.0
        write(*,*)'Lat,Lon upper left corner:',yt,xt
        call xy_ll_S(ph_lim(2),th_lim(1),xt,yt)
        if(xt>180.0)xt=xt-360.0
        write(*,*)'Lat,Lon lower right corner:',yt,xt
        call xy_ll_S(ph_lim(2),th_lim(2),xt,yt)
        if(xt>180.0)xt=xt-360.0
        write(*,*)'Lat,Lon upper right corner:',yt,xt
      elseif(trim(xy_ll_sub)=='xy_ll_CATs')then
        call xy_ll_CATs(ph_lim(1),th_lim(1),xt,yt)
        if(xt>180.0)xt=xt-360.0
        write(*,*)'Lat,Lon lower left corner:',yt,xt
        call xy_ll_CATs(ph_lim(1),th_lim(2),xt,yt)
        if(xt>180.0)xt=xt-360.0
        write(*,*)'Lat,Lon upper left corner:',yt,xt
        call xy_ll_CATs(ph_lim(2),th_lim(1),xt,yt)
        if(xt>180.0)xt=xt-360.0
        write(*,*)'Lat,Lon lower right corner:',yt,xt
        call xy_ll_CATs(ph_lim(2),th_lim(2),xt,yt)
        if(xt>180.0)xt=xt-360.0
        write(*,*)'Lat,Lon upper right corner:',yt,xt
      endif
    endif
    
    write(*,*)'Constituents: ',c_id_mod(1:nc)
    if(trim(xy_ll_sub)/='')then
      write(*,*)'Model is on uniform grid in km'
      write(*,*)'Function to convert x,y to lat,lon:',trim(xy_ll_sub)
    endif 
!
    if(zuv=='z')write(*,*)'Output elevations (m)'
    if(zuv=='U')write(*,*)'Output WE transport (m^2/s)'
    if(zuv=='u')write(*,*)'Output WE velocity  (cm/s)'
    if(zuv=='V')write(*,*)'Output SN transport (m^2/s)'
    if(zuv=='v')write(*,*)'Output SN velocity  (cm/s)'
!
    !if(zuv=='z')write(11,*)'Elevations (m)'
    !if(zuv=='U')write(11,*)'WE transport (m^2/s)'
    !if(zuv=='u')write(11,*)'WE velocity  (cm/s)'
    !if(zuv=='V')write(11,*)'SN transport (m^2/s)'
    !if(zuv=='v')write(11,*)'SN velocity  (cm/s)'
!
    if(ncon==0)then
      ibl=1
      ncon=nc
      c_id=c_id_mod
      write(*,*)'Constituents to include: ',c_id(1:ncon)
    endif
!
    allocate(cind(ncon))
    call def_con_ind(c_id,ncon,c_id_mod,nc,cind)

! Latitude/Longitude
!    ndat=0
!    open(unit=1,file=lltname,status='old',err=1)
!3   read(1,*,end=2)dum,dum
    !ndat=ndat+1
    !goto 3
!2   rewind(1)
    !allocate(lat(ndat),lon(ndat),lon0(ndat))
    if(trim(xy_ll_sub)/='')allocate(x(ndat),y(ndat))
    do k=1,ndat
      !read(1,*)lat(k),lon(k)
      if(trim(xy_ll_sub)=='xy_ll_N')then
        call ll_xy_N(lon(k),lat(k),x(k),y(k))
        !write(*,*)lat(k),lon(k),x(k),y(k)
      elseif(trim(xy_ll_sub)=='xy_ll_S')then
        call ll_xy_S(lon(k),lat(k),x(k),y(k))
      elseif(trim(xy_ll_sub)=='xy_ll_CATs')then
        call ll_xy_CATs(lon(k),lat(k),x(k),y(k))
      endif
      lon0(k)=lon(k)
      if(trim(xy_ll_sub)=='')then ! check on lon convention
        if(lon(k)>ph_lim(2))lon(k)=lon(k)-360
        if(lon(k)<ph_lim(1))lon(k)=lon(k)+360
      endif
    enddo ! k
    !close(1)

! Read model
    allocate(z(ncon,n,m),z1(ncon),mask(n,m))
    write(*,'(a)')' Reading model...'
    call rd_mod_body(modname,zuv,cind,ncon,z,n,m,mask)
    write(*,*)'done'
    if(zuv=='z'.and.geo)then
    call rd_mod_header1(lname,nl,ml,ncl,lth_lim,lph_lim,lc_id)
    allocate(lcind(ncon))
    call def_con_ind(c_id,ncon,lc_id,ncl,lcind)
    allocate(zl(ncon,nl,ml),zl1(ncon),mzl(nl,ml))
    write(*,*)
    write(*,'(a)')' Reading file to apply load corrections...'
    call rd_mod_body1(lname,'z',lcind,ncon,zl,nl,ml,mzl)
    write(*,*)'done'
    endif
    if(zuv=='u'.or.zuv=='v')then
      allocate(depth(n,m),dtmp(n,m),mz(n,m))
! currents case: need to read grid
      open(unit=1,file=modname,status='old')
      read(1,*)gname
      read(1,*)gname
      read(1,'(a80)')gname
      gname=rmCom(gname)
      close(1)
      open(unit=1,file=trim(gname),status='old',&
            form='unformatted',err=4)
      read(1)n0,m0 ! ignore the rest of record
      if(n0/=n.or.m0/=m)then
        write(*,*)'Wrong grid in ',trim(gname)
        write(*,*)'Grid size and model size are different!'
        stop
      endif
      read(1) ! pass iobc
      read(1) depth
      read(1) mz
      close(1)
      dtmp=depth
      deallocate(depth)
    endif
    
    allocate(phase(ncon))
    
!! output format
!      fmt='(f9.3,x,f9.3,x'
!      c1=','
!      do ic=1,ncon
!       if(ic==ncon)c1=')'
!       if(APRI)then
!        fmt=trim(fmt)//'f8.3,x,f8.1,x'//c1
!       else
!        fmt=trim(fmt)//'f8.3,x,f8.3,x'//c1
!       endif
!      enddo
!
    !if(APRI)then
    !  write(11,*)'  Lat     Lon       ',&
    !              (trim(c_id(ic)),'_amp  ',&
    !              trim(c_id(ic)),'_ph   ',ic=1,ncon)
    !else
    !  write(11,*)'  Lat       Lon     ',&
    !             (trim(c_id(ic)),'_Re   ',&
    !              trim(c_id(ic)),'_Im   ',ic=1,ncon)
    !endif
    c1=zuv ! since interp change zuv (U->u, V->v)
    latp=0.0
    lonp=0.0
    do k=1,ndat
      if(lat(k)/=latp.or.lon(k)/=lonp)then
        if(trim(xy_ll_sub)=='')then
         call interp(z,ncon,n,m,mask,th_lim,ph_lim,&
                    lat(k),lon(k),z1,ierr,c1)
        else
         z1(1)=-1
         call interp(z,ncon,n,m,mask,th_lim,ph_lim,&
                    y(k),x(k),z1,ierr,c1)
        endif
        amp(k,1:ncon)=z1(1:ncon)
        if(ierr==0)then
         if(zuv=='u'.or.zuv=='v')then
          if(trim(xy_ll_sub)=='')then
           call interp(dtmp,1,n,m,mz,th_lim,ph_lim,&
                       lat(k),lon(k),d1,ierr1,'z')
          else
           d1=-1
           call interp(dtmp,1,n,m,mz,th_lim,ph_lim,&
                       y(k),x(k),d1,ierr1,'z')
          endif
          z1=z1/real(d1)*100.0 ! currents cm/s
         elseif(zuv=='z'.and.geo)then       
           call interp(zl,ncon,nl,ml,mzl,lth_lim,lph_lim,&
                    lat(k),lon(k),zl1,ierr1,'z')
           z1=z1+zl1     ! apply load correction to get geocentric tide
         endif  
        endif
        if(ierr==0)then
         if(APRI)then
          phase=atan2(-imag(z1),real(z1))*rad2deg
          pha(k,1:ncon)=phase(1:ncon)          
          !write(11,fmt)lat(k),lon0(k),&
          !      (abs(z1(ic)),phase(ic),ic=1,ncon)
         else
           !write(11,fmt)lat(k),lon0(k),&
           !         (real(z1(ic)),imag(z1(ic)),ic=1,ncon)
         endif
        else
       !   write(11,'(1x,f8.3,x,f8.3,a70)')lat(k),lon(k),&
       !'************* Site is out of model grid OR land ***************'
        endif
        latp=lat(k)
        lonp=lon(k)
      endif  
    enddo !k=1,ndat

!--- Deallocate variables ----------------------------------
    deallocate(z,z1,mask,cind,lon0)
    if(zuv=='u'.or.zuv=='v')deallocate(dtmp,mz)
    if(trim(xy_ll_sub)/='')deallocate(x,y)
    if(ibl==1.and.geo)then
      ncon=0
      deallocate(zl,zl1,mzl,lcind)
    endif
    !enddo ! imod
    !close(11)
    !write(*,*)'Results are in ',trim(outname)
    if(geo.and.ibl==0)deallocate(zl,zl1,mzl,lcind)
    return
    
!--- Error messages ---------------------------------------    
1   write(*,*)'Lat Lon file ',trim(lltname),' not found'
    write(*,*)'Check setup file, line 2.'
    stop
4   write(*,*)'Grid file ',trim(gname),' not found'
    write(*,*)'Check file ',trim(modname)
    stop
11  write(*,*)'File ''model.list'' was NOT found...'
    write(*,*)'TO CREATE please do:'
    write(*,*)'ls -1 DATA/Model_*>model.list'
    stop
    endsubroutine extract_hc
      
endmodule otps_extract_hc
