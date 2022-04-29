   
!*************************************************************
    subroutine bnd_vel_flux(nbndcells,icells,kfaces)
! Calculatest the velocity at a flux boundary condition
!*************************************************************     
    use geo_def, only: dsxy,cell2cell
    use flow_def, only: visk,hk,flux,iwet,u,v,su,sv,sp,acoef
    use comp_lib, only: upwindcoef
    use prec_def, only: ikind
    
    implicit none
    !Input/Output
    integer,intent(in) :: nbndcells
    integer,intent(in) :: icells(nbndcells),kfaces(nbndcells)
    !Internal Variables
    integer :: i,j,k,nck
    real(ikind) :: ddk
    
    do j=1,nbndcells
      i=icells(j)
      k=kfaces(j)
      nck=cell2cell(k,i)
      ddk=visk(k,i)*hk(k,i)*dsxy(k,i)
      acoef(k,i)=upwindcoef(ddk,flux(k,i))
      su(i)=su(i)+acoef(k,i)*u(nck)
      sv(i)=sv(i)+acoef(k,i)*v(nck)
      sp(i)=sp(i)-acoef(k,i)
      acoef(k,i)=0.0
    enddo    

    return
    end subroutine bnd_vel_flux
    
!**************************************************************************
    subroutine bnd_vel_openwse(nbndcells,icells,kfaces)
! Calculates the velocity at an open wse boundary condition
! written by Alex Sanchez, USACE-CHL
!*************************************************************************
    use size_def, only: ncellsimple,ncelljoint,ncellpoly
    use geo_def, only: ds,dsxy,idirface,cell2cell,fnx,fny
    use flow_def, only: visk,flux,iwet,u,v,h,hk,su,sv,sp,acoef
    use comp_lib, only: upwindcoef
    use prec_def, only: ikind
    
    implicit none
    !Input/Output
    integer,intent(in) :: nbndcells
    integer,intent(in) :: icells(nbndcells),kfaces(nbndcells)
    !Internal Variables
    integer :: i,j,k,nck
    real(ikind) :: ddk
    
    interface
      function forcex(i)
        use prec_def
        integer,intent(in) :: i
        real(ikind) :: forcex
      end function
    endinterface
    
    interface
      function forcey(i)
        use prec_def
        integer,intent(in) :: i
        real(ikind) :: forcey
      end function
    endinterface
    
    do j=1,nbndcells
      i=icells(j)
      k=kfaces(j)    
      nck=cell2cell(k,i)
      ddk=visk(k,i)*hk(k,i)*dsxy(k,i)
      acoef(k,i)=upwindcoef(ddk,flux(k,i))
      if(abs(fny(k,i))>0.9999)then  !North/South Boundary
        u(nck)=iwet(i)*u(i)*h(i)/h(nck)
        v(nck)=fny(k,i)*flux(k,i)/(ds(k,i)*hk(k,i))
        sv(i)=sv(i)+acoef(k,i)*v(nck)
      elseif(abs(fnx(k,i))>0.9999)then !East/West Boundary
        u(nck)=fnx(k,i)*flux(k,i)/(ds(k,i)*hk(k,i))
        v(nck)=iwet(i)*v(i)*h(i)/h(nck)
        su(i)=su(i)+acoef(k,i)*u(nck)
      else                             !Oblique Boundary
        !u(nck)=iwet(i)*u(i)*h(i)/h(nck)
        !v(nck)=iwet(i)*v(i)*h(i)/h(nck)
        u(nck)=fnx(k,i)*flux(k,i)/(ds(k,i)*hk(k,i))
        v(nck)=fny(k,i)*flux(k,i)/(ds(k,i)*hk(k,i))
        !!v(nck)=fny(k,i)*flux(k,i)/ds(k,i)/hk(k,i)    !y-directed
        !!v(nck)=v(nck)+abs(fnx(k,i))*v(i)*h(i)/h(nck) !x-directed  
        su(i)=su(i)+acoef(k,i)*u(nck)
        sv(i)=sv(i)+acoef(k,i)*v(nck)        
      endif
      sp(i)=sp(i)-acoef(k,i)
      acoef(k,i)=0.0
    enddo
    
    return
    end subroutine bnd_vel_openwse

!**************************************************************************
    subroutine bnd_vel_openwsevel(nbndcells,icells,kfaces,ubnd,vbnd)
! Calculates the velocity at an open wse and velocity boundary condition
! written by Alex Sanchez, USACE-CHL
!*************************************************************************
    use geo_def, only: ds,dsxy,cell2cell,fnx,fny
    use flow_def, only: visk,flux,iwet,u,v,h,hk,su,sv,sp,acoef
    use comp_lib, only: upwindcoef
    use prec_def, only: ikind
    
    implicit none
    !Input/Output
    integer,intent(in) :: nbndcells
    integer,intent(in) :: icells(nbndcells),kfaces(nbndcells)
    real(ikind),intent(in) :: ubnd(nbndcells),vbnd(nbndcells)
    !Internal Variables
    integer :: i,j,k,nck
    real(ikind) :: ddk
    
    do j=1,nbndcells
      i=icells(j)
      k=kfaces(j)
      nck=cell2cell(k,i)
      if(flux(k,i)<0.0)then !Inflow
        ddk=visk(k,i)*hk(k,i)*dsxy(k,i)
        acoef(k,i)=upwindcoef(ddk,flux(k,i))
        u(nck)=iwet(i)*ubnd(j)
        v(nck)=iwet(i)*vbnd(j)
        !u(nck)=fnx(k,i)*flux(k,i)/(ds(k,i)*hk(k,i))
        !v(nck)=fny(k,i)*flux(k,i)/(ds(k,i)*hk(k,i))
        su(i)=su(i)+acoef(k,i)*u(nck)
        sv(i)=sv(i)+acoef(k,i)*v(nck)
        sp(i)=sp(i)-acoef(k,i)
      else !Outflow
        u(nck)=iwet(i)*u(i)*h(i)/h(nck)
        v(nck)=iwet(i)*v(i)*h(i)/h(nck)
        !u(nck)=fnx(k,i)*flux(k,i)/(ds(k,i)*hk(k,i))
        !v(nck)=fny(k,i)*flux(k,i)/(ds(k,i)*hk(k,i))
      endif
      acoef(k,i)=0.0
    enddo    
    
    return
    end subroutine bnd_vel_openwsevel    
    
!**************************************************************************
    subroutine bnd_vel_openwse_test(nbndcells,icells,kfaces)
! Calculates the velocity at an open wse boundary condition
! written by Alex Sanchez, USACE-CHL
!*************************************************************************
    use size_def, only: ncellsimple,ncelljoint,ncellpoly
    use geo_def, only: ds,dsxy,idirface,cell2cell,fnx,fny
    use flow_def, only: visk,flux,iwet,u,v,h,hk,su,sv,sp,acoef,uv,fc,&
        dux,duy,dvx,dvy,dpx,dpy,vis,us,vs
    use fric_def, only: cfrict
    use comp_lib, only: upwindcoef
    use der_lib, only: dx2d,dy2d
    use der_def, only: gow
    use bnd_def, only: veldamp
    use prec_def, only: ikind
    
    implicit none
    !Input/Output
    integer,intent(in) :: nbndcells
    integer,intent(in) :: icells(nbndcells),kfaces(nbndcells)
    !Internal Variables
    integer :: i,j,k,nck
    real(ikind) :: ddk,d2uiy2,d2vix2
    
    interface
      function forcex(i)
        use prec_def
        integer,intent(in) :: i
        real(ikind) :: forcex
      end function
    endinterface
    
    interface
      function forcey(i)
        use prec_def
        integer,intent(in) :: i
        real(ikind) :: forcey
      end function
    endinterface
    
    if(ncellpoly==0)then
    do j=1,nbndcells
      i=icells(j)
      k=kfaces(j)    
      nck=cell2cell(k,i)
      ddk=visk(k,i)*hk(k,i)*dsxy(k,i)
      acoef(k,i)=upwindcoef(ddk,flux(k,i))
      if(abs(fny(k,i))>0.9999)then  !North/South Boundary
        u(nck)=iwet(i)*u(i)*h(i)/h(nck) !Inward flow
        if(veldamp<1.0e-6 .or. flux(k,i)>0.0)then !Note: Always use for outward flow
          v(nck)=fny(k,i)*flux(k,i)/(ds(k,i)*hk(k,i))
        else
          call dx2d(gow,i,dvx,d2vix2)  
          !v(nck)=(forcey(i)+(-fc(i)*u(i)-dpy(i)-vis(i)*d2vix2)*h(i))/(cfrict(i)*(uv(i)+damp*v(i))+1.0e-3)
          v(nck)=(forcey(i)+(-fc(i)*u(i)-dpy(i)-vis(i)*d2vix2)*h(i))/(cfrict(i)*uv(i)+veldamp)
          v(nck)=v(nck)+vs(i) !Add stokes velocity
          !v(nck)=(damp*forcey(i)+(-fc(i)*u(i)-dpy(i)-vis(i)*d2vix2)*h(i))/(cfrict(i)*uv(i)+1.0e-6)
        endif
        sv(i)=sv(i)+acoef(k,i)*v(nck) !Inward flow
      elseif(abs(fnx(k,i))>0.9999)then !East/West Boundary
        if(veldamp<1.0e-6 .or. flux(k,i)>0.0)then  !Note: Always use for outward flow
          u(nck)=fnx(k,i)*flux(k,i)/(ds(k,i)*hk(k,i))
        else
          call dy2d(gow,i,dux,d2uiy2)
          !u(nck)=(forcex(i)+(fc(i)*v(i)-dpx(i)-vis(i)*d2uiy2)*h(i))/(cfrict(i)*(uv(i)+damp*u(i))+1.0e-3)
          u(nck)=(forcex(i)+(fc(i)*v(i)-dpx(i)-vis(i)*d2uiy2)*h(i))/(cfrict(i)*uv(i)+veldamp)
          u(nck)=u(nck)+us(i) !Add stokes velocity
          !u(nck)=(damp*forcex(i)+(fc(i)*v(i)-dpx(i)-vis(i)*d2uiy2)*h(i))/(cfrict(i)*uv(i)+1.0e-6)
        endif
        v(nck)=iwet(i)*v(i)*h(i)/h(nck)
        su(i)=su(i)+acoef(k,i)*u(nck)
      else                             !Oblique Boundary
        !u(nck)=iwet(i)*u(i)*h(i)/h(nck)
        !v(nck)=iwet(i)*v(i)*h(i)/h(nck)
        u(nck)=fnx(k,i)*flux(k,i)/(ds(k,i)*hk(k,i))
        v(nck)=fny(k,i)*flux(k,i)/(ds(k,i)*hk(k,i))
        !!v(nck)=fny(k,i)*flux(k,i)/(ds(k,i)*hk(k,i))    !y-directed
        !!v(nck)=v(nck)+abs(fnx(k,i))*v(i)*h(i)/h(nck) !x-directed  
        su(i)=su(i)+acoef(k,i)*u(nck)
        sv(i)=sv(i)+acoef(k,i)*v(nck)        
      endif
      sp(i)=sp(i)-acoef(k,i)
      acoef(k,i)=0.0
    enddo
    else
    do j=1,nbndcells
      i=icells(j)
      k=kfaces(j)    
      nck=cell2cell(k,i)
      ddk=visk(k,i)*hk(k,i)*dsxy(k,i)
      acoef(k,i)=upwindcoef(ddk,flux(k,i))
      !u(nck)=iwet(i)*u(i)*h(i)/h(nck)
      !v(nck)=iwet(i)*v(i)*h(i)/h(nck)
      u(nck)=fnx(k,i)*flux(k,i)/(ds(k,i)*hk(k,i))
      v(nck)=fny(k,i)*flux(k,i)/(ds(k,i)*hk(k,i))
      !!v(nck)=fny(k,i)*flux(k,i)/(ds(k,i)*hk(k,i))    !y-directed
      !!v(nck)=v(nck)+abs(fnx(k,i))*v(i)*h(i)/h(nck) !x-directed  
      su(i)=su(i)+acoef(k,i)*u(nck)
      sv(i)=sv(i)+acoef(k,i)*v(nck)
      sp(i)=sp(i)-acoef(k,i)
      acoef(k,i)=0.0
    enddo
    endif
    
    return
    end subroutine bnd_vel_openwse_test    