!***********************************************************************
    subroutine flow_imp_flux_skew
! Intercell flux calculation for implicit temporal scheme
! using a Rhie and Chow type momentum interpolation scheme
!
! written by Weiming Wu, NCCHE, Oct. 2008
! Modified by Alex Sanchez, USACE, 2013
!***********************************************************************
    use size_def
    use geo_def
    use der_def
    use der_lib, only: der_grad_eval,der_gradvec_eval
    use interp_def
    use flow_def
    use comp_lib
    use comvarbl, only: relax,skewcor
    use const_def, only: small
    use prec_def
    implicit none
    integer :: i,j,k,ii,nck,jcn
    real(ikind) :: fk,fp,apuvolpf,Huf,Hvf,dpf
    
    call der_gradvec_eval(gow,nlim,Hu,Hv,dHux,dHuy,dHvx,dHvy)
    call der_grad_eval(gow,nlim,apuareap,dapuareapx,dapuareapy)
    !!call der_gradvec_eval(gow,nlim,Hu,Hv,dHux,dHuy,dHvx,dHvy)
    !!call der_grad_eval(gow,nlim,apuareap,dapuareapx,dapuareapy)
!$OMP PARALLEL
    !--- Flux Calculation ------------------ 
!$OMP DO PRIVATE(i,ii,j,k,nck,fk,fp)
    do ii=1,ncellsimple
      i=idcellsimple(ii)
      do j=1,nxface(i) !east and west sides
        k=kxface(j,i)    
        nck=cell2cell(k,i) !Forward connectivity     
        if(iwet(i)*iwet(nck)==0)then !Closed Boundaries
          flux(k,i)=0.0
        elseif(nck>ncells)then !Open Boundaries
          flux(k,i)=hk(k,i)*ds(k,i)*(fnx(k,i)*Hu(i) &
              -apuareap(i)*hk(k,i)*(p(nck)-p(i))/dn(k,i))
        else
          fk=fintp(k,i); fp=1.0-fk
          flux(k,i)=hk(k,i)*ds(k,i)*(fnx(k,i)*(fk*Hu(nck)+fp*Hu(i)) &
              -(fk*apuareap(nck)+fp*apuareap(i))*hk(k,i)*(p(nck)-p(i))/dn(k,i))
        endif
        flux(llec2llec(k,i),nck)=-flux(k,i)
      enddo
      do j=1,nyface(i) !north and south faces
        k=kyface(j,i)
        nck=cell2cell(k,i) !Forward connectivity
        if(iwet(i)*iwet(nck)==0)then !Closed Boundaries
          flux(k,i)=0.0
        elseif(nck>ncells)then !Open Boundaries
          !fk=fintp(k,i); fp=1.0-fk  
          !flux(k,i)=hk(k,i)*ds(k,i)*(fk*v(nck)+fp*v(i)) !Linear interpolation
          !flux(k,i)=ds(k,i)*(fk*v(nck)*h(nck)+fp*v(i)*h(i)) !Linear interpolation  
          flux(k,i)=hk(k,i)*ds(k,i)*(fny(k,i)*Hv(i) &
             -apuareap(i)*hk(k,i)*(p(nck)-p(i))/dn(k,i))
        else
          fk=fintp(k,i); fp=1.0-fk
          flux(k,i)=hk(k,i)*ds(k,i)*(fny(k,i)*(fk*Hv(nck)+fp*Hv(i)) &
              -(fk*apuareap(nck)+fp*apuareap(i))*hk(k,i)*(p(nck)-p(i))/dn(k,i))
        endif
        flux(llec2llec(k,i),nck)=-flux(k,i)
      enddo
    enddo
!$OMP END DO     
!!$OMP END DO NOWAIT
!$OMP DO PRIVATE(i,ii,j,k,nck,jcn,fk,fp,apuvolpf,Huf,Hvf,dpf)
    do ii=1,ncelljoint
      i=idcelljoint(ii)
      do j=1,nxface(i) !east and west sides
        k=kxface(j,i)     
        nck=cell2cell(k,i) !Forward connectivity
        jcn=llec2llec(k,i) !Backward connectivity
        if(iwet(i)*iwet(nck)==0)then !Closed Boundaries
          flux(k,i)=0.0
        elseif(nck>ncells)then !Open Boundaries
          flux(k,i)=hk(k,i)*ds(k,i)*(fnx(k,i)*Hu(i) &
              -apuareap(i)*hk(k,i)*(p(nck)-p(i))/dn(k,i))
        else  !Internal cells
          fk=fintp(k,i); fp=1.0-fk
          if(abs(rpy(k,i))>1.0e-6)then !Coarse to fine
            apuvolpf=fk*apuareap(nck)+fp*(apuareap(i)+rpy(k,i)*dapuareapy(i))
            Huf=fk*Hu(nck)+fp*(Hu(i)+rpy(k,i)*dHuy(i))
            dpf=(p(nck)-p(i)-rpy(k,i)*dpy(i))/dn(k,i)
          elseif(abs(rpy(jcn,nck))>1.0e-6)then !Fine to coarse
            apuvolpf=fk*(apuareap(nck)+rpy(jcn,nck)*dapuareapy(nck))+fp*apuareap(i)  
            Huf=fk*(Hu(nck)+rpy(jcn,nck)*dHuy(nck))+fp*Hu(i)
            dpf=(p(nck)+rpy(jcn,nck)*dpy(nck)-p(i))/dn(k,i)
          else !no refinement
            apuvolpf=fk*apuareap(nck)+fp*apuareap(i)
            Huf=fk*Hu(nck)+fp*Hu(i)
            dpf=(p(nck)-p(i))/dn(k,i)
          endif
          flux(k,i)=hk(k,i)*ds(k,i)*(fnx(k,i)*Huf-apuvolpf*hk(k,i)*dpf)
        endif
        flux(jcn,nck)=-flux(k,i)
      enddo
      do j=1,nyface(i) !north and south faces
        k=kyface(j,i)
        nck=cell2cell(k,i) !Forward connectivity 
        jcn=llec2llec(k,i) !Backward connectivity
        if(iwet(i)*iwet(nck)==0)then !Closed Boundaries
          flux(k,i)=0.0
        elseif(nck>ncells)then !Open Boundaries
          flux(k,i)=hk(k,i)*ds(k,i)*(fny(k,i)*Hv(i) &
              -apuareap(i)*hk(k,i)*(p(nck)-p(i))/dn(k,i))
        else  !Internal cells
          fk=fintp(k,i); fp=1.0-fk
          if(abs(rpx(k,i))>1.0e-6)then !Coarse to fine
            apuvolpf=fk*apuareap(nck)+fp*(apuareap(i)+rpx(k,i)*dapuareapx(i))
            Hvf=fk*Hv(nck)+fp*(Hv(i)+rpx(k,i)*dHvx(i))
            dpf=(p(nck)-p(i)-rpx(k,i)*dpx(i))/dn(k,i)
          elseif(abs(rpx(jcn,nck))>1.0e-6)then !Fine to coarse
            apuvolpf=fk*(apuareap(nck)+rpx(jcn,nck)*dapuareapx(nck))+fp*apuareap(i)
            Hvf=fk*(Hv(nck)+rpx(jcn,nck)*dHvx(nck))+fp*Hv(i)
            dpf=(p(nck)+rpx(jcn,nck)*dpx(nck)-p(i))/dn(k,i)
          else !no refinement
            apuvolpf=fk*apuareap(nck)+fp*apuareap(i)
            Hvf=fk*Hv(nck)+fp*Hv(i)
            dpf=(p(nck)-p(i))/dn(k,i)
          endif
          flux(k,i)=hk(k,i)*ds(k,i)*(fny(k,i)*Hvf-apuvolpf*hk(k,i)*dpf)
        endif
        flux(jcn,nck)=-flux(k,i)
      enddo
    enddo
!$OMP END DO
!$OMP DO PRIVATE(i,j,k,nck,jcn,fk,fp,apuvolpf,Huf,Hvf,dpf)
    do i=1,ncellpoly
      do j=1,nxyface(i) !No Repeat of cell faces
        k=kxyface(j,i)     
        nck=cell2cell(k,i) !Forward connectivity 
        jcn=llec2llec(k,i) !Backward connectivity
        if(iwet(i)*iwet(nck)==0)then !Closed Boundaries
          flux(k,i)=0.0
        elseif(nck>ncells)then !Open Boundaries
          flux(k,i)=hk(k,i)*ds(k,i)*((fnx(k,i)*Hu(i)+fny(k,i)*Hv(i)) &
             -apuareap(i)*h(i)*(p(nck)-p(i))/dn(k,i))
        !!elseif(maxval(cell2cell(1:ncface(i),i))>ncells)then
        !!  fk=fintp(k,i); fp=1.0-fk
        !!  !apuvolpf=fk*apuareap(nck)+fp*apuareap(i)
        !!  apuvolpf=min(apuareap(nck),apuareap(i))
        !!  Huf=fk*Hu(nck)+fp*Hu(i)
        !!  Hvf=fk*Hv(nck)+fp*Hv(i)
        !!  flux(k,i)=hk(k,i)*ds(k,i)*((fnx(k,i)*Huf+fny(k,i)*Hvf) &
        !!     -apuvolpf*hk(k,i)*(p(nck)-p(i))/dn(k,i))  
        else
          fk=fintp(k,i); fp=1.0-fk
          !apuvolpf=fk*(apuareap(nck)+rpx(jcn,nck)*dapuareapx(nck)+rpy(jcn,nck)*dapuareapy(nck)) &
          !        +fp*(apuareap(i)+rpx(k,i)*dapuareapx(i)+rpy(k,i)*dapuareapy(i))
          !apuvolpf=min(apuareap(nck),apuareap(i)) !more stable
          apuvolpf=fk*apuareap(nck)+fp*apuareap(i)
          Huf=fk*(Hu(nck)+rpx(jcn,nck)*dHux(nck)+rpy(jcn,nck)*dHuy(nck)) &
             +fp*(Hu(i)+rpx(k,i)*dHux(i)+rpy(k,i)*dHuy(i))
          Hvf=fk*(Hv(nck)+rpx(jcn,nck)*dHvx(nck)+rpy(jcn,nck)*dHvy(nck)) &
             +fp*(Hv(i)+rpx(k,i)*dHvx(i)+rpy(k,i)*dHvy(i))
          dpf=(p(nck)+rpx(jcn,nck)*dpx(nck)+rpy(jcn,nck)*dpy(nck) &
              -p(i)-rpx(k,i)*dpx(i)-rpy(k,i)*dpy(i))/dn(k,i)
          flux(k,i)=hk(k,i)*ds(k,i)*((fnx(k,i)*Huf+fny(k,i)*Hvf)-apuvolpf*hk(k,i)*dpf)
        endif
        flux(jcn,nck)=-flux(k,i)
      enddo      
    enddo
!$OMP END DO
!$OMP END PARALLEL

    return
    end subroutine flow_imp_flux_skew
    
!***********************************************************************
    subroutine flow_imp_flux_skew_rc
! Intercell flux calculation for implicit temporal scheme
! using a Rhie and Chow type momentum interpolation scheme
!
! written by Weiming Wu, NCCHE, Oct. 2008
! Modified by Alex Sanchez, USACE, 2013
!***********************************************************************
    use size_def
    use geo_def
    use der_def
    use der_lib, only: der_grad_eval,der_gradvec_eval
    use interp_def
    use flow_def
    use comp_lib
    use comvarbl, only: relax,skewcor
    use const_def, only: small
    use prec_def
    implicit none
    integer :: i,j,k,ii,nck,jcn
    real(ikind) :: fk,fp,apuvolpf,Huf,Hvf,dpf,dpfi
    
    call der_gradvec_eval(gow,nlim,Hu,Hv,dHux,dHuy,dHvx,dHvy)
    call der_grad_eval(gow,nlim,apuareap,dapuareapx,dapuareapy)
    !!call der_gradvec_eval(gow,nlim,Hu,Hv,dHux,dHuy,dHvx,dHvy)
    !!call der_grad_eval(gow,nlim,apuareap,dapuareapx,dapuareapy)
!$OMP PARALLEL
    !--- Flux Calculation ------------------ 
!$OMP DO PRIVATE(i,ii,j,k,nck,fk,fp)
    do ii=1,ncellsimple
      i=idcellsimple(ii)
      do j=1,nxface(i) !east and west sides
        k=kxface(j,i)    
        nck=cell2cell(k,i) !Forward connectivity     
        if(iwet(i)*iwet(nck)==0)then !Closed Boundaries
          flux(k,i)=0.0
        elseif(nck>ncells)then !Open Boundaries
          flux(k,i)=hk(k,i)*ds(k,i)*(fnx(k,i)*Hu(i) &
              -apuareap(i)*hk(k,i)*(p(nck)-p(i))/dn(k,i))
        else
          fk=fintp(k,i); fp=1.0-fk
          flux(k,i)=hk(k,i)*ds(k,i)*(fnx(k,i)*(fk*Hu(nck)+fp*Hu(i)) &
              -(fk*apuareap(nck)+fp*apuareap(i))*hk(k,i)*(p(nck)-p(i))/dn(k,i))
        endif
        flux(llec2llec(k,i),nck)=-flux(k,i)
      enddo
      do j=1,nyface(i) !north and south faces
        k=kyface(j,i)
        nck=cell2cell(k,i) !Forward connectivity
        if(iwet(i)*iwet(nck)==0)then !Closed Boundaries
          flux(k,i)=0.0
        elseif(nck>ncells)then !Open Boundaries
          !fk=fintp(k,i); fp=1.0-fk  
          !flux(k,i)=hk(k,i)*ds(k,i)*(fk*v(nck)+fp*v(i)) !Linear interpolation
          !flux(k,i)=ds(k,i)*(fk*v(nck)*h(nck)+fp*v(i)*h(i)) !Linear interpolation  
          flux(k,i)=hk(k,i)*ds(k,i)*(fny(k,i)*Hv(i) &
             -apuareap(i)*hk(k,i)*(p(nck)-p(i))/dn(k,i))
        else
          fk=fintp(k,i); fp=1.0-fk
          flux(k,i)=hk(k,i)*ds(k,i)*(fny(k,i)*(fk*Hv(nck)+fp*Hv(i)) &
              -(fk*apuareap(nck)+fp*apuareap(i))*hk(k,i)*(p(nck)-p(i))/dn(k,i))
        endif
        flux(llec2llec(k,i),nck)=-flux(k,i)
      enddo
    enddo
!$OMP END DO     
!!$OMP END DO NOWAIT
!$OMP DO PRIVATE(i,ii,j,k,nck,jcn,fk,fp,apuvolpf,Huf,Hvf,dpf)
    do ii=1,ncelljoint
      i=idcelljoint(ii)
      do j=1,nxface(i) !east and west sides
        k=kxface(j,i)     
        nck=cell2cell(k,i) !Forward connectivity
        jcn=llec2llec(k,i) !Backward connectivity
        if(iwet(i)*iwet(nck)==0)then !Closed Boundaries
          flux(k,i)=0.0
        elseif(nck>ncells)then !Open Boundaries
          flux(k,i)=hk(k,i)*ds(k,i)*(fnx(k,i)*Hu(i) &
              -apuareap(i)*hk(k,i)*(p(nck)-p(i))/dn(k,i))
        else  !Internal cells
          fk=fintp(k,i); fp=1.0-fk
          if(abs(rpy(k,i))>1.0e-6)then !Coarse to fine
            apuvolpf=fk*apuareap(nck)+fp*(apuareap(i)+rpy(k,i)*dapuareapy(i))
            Huf=fk*Hu(nck)+fp*(Hu(i)+rpy(k,i)*dHuy(i))
            dpf=(p(nck)-p(i)-rpy(k,i)*dpy(i))/dn(k,i)
          elseif(abs(rpy(jcn,nck))>1.0e-6)then !Fine to coarse
            apuvolpf=fk*(apuareap(nck)+rpy(jcn,nck)*dapuareapy(nck))+fp*apuareap(i)  
            Huf=fk*(Hu(nck)+rpy(jcn,nck)*dHuy(nck))+fp*Hu(i)
            dpf=(p(nck)+rpy(jcn,nck)*dpy(nck)-p(i))/dn(k,i)
          else !no refinement
            apuvolpf=fk*apuareap(nck)+fp*apuareap(i)
            Huf=fk*Hu(nck)+fp*Hu(i)
            dpf=(p(nck)-p(i))/dn(k,i)
          endif
          flux(k,i)=hk(k,i)*ds(k,i)*(fnx(k,i)*Huf-apuvolpf*hk(k,i)*dpf)
        endif
        flux(jcn,nck)=-flux(k,i)
      enddo
      do j=1,nyface(i) !north and south faces
        k=kyface(j,i)
        nck=cell2cell(k,i) !Forward connectivity 
        jcn=llec2llec(k,i) !Backward connectivity
        if(iwet(i)*iwet(nck)==0)then !Closed Boundaries
          flux(k,i)=0.0
        elseif(nck>ncells)then !Open Boundaries
          flux(k,i)=hk(k,i)*ds(k,i)*(fny(k,i)*Hv(i) &
              -apuareap(i)*hk(k,i)*(p(nck)-p(i))/dn(k,i))
        else  !Internal cells
          fk=fintp(k,i); fp=1.0-fk
          if(abs(rpx(k,i))>1.0e-6)then !Coarse to fine
            apuvolpf=fk*apuareap(nck)+fp*(apuareap(i)+rpx(k,i)*dapuareapx(i))
            Hvf=fk*Hv(nck)+fp*(Hv(i)+rpx(k,i)*dHvx(i))
            dpf=(p(nck)-p(i)-rpx(k,i)*dpx(i))/dn(k,i)
          elseif(abs(rpx(jcn,nck))>1.0e-6)then !Fine to coarse
            apuvolpf=fk*(apuareap(nck)+rpx(jcn,nck)*dapuareapx(nck))+fp*apuareap(i)
            Hvf=fk*(Hv(nck)+rpx(jcn,nck)*dHvx(nck))+fp*Hv(i)
            dpf=(p(nck)+rpx(jcn,nck)*dpx(nck)-p(i))/dn(k,i)
          else !no refinement
            apuvolpf=fk*apuareap(nck)+fp*apuareap(i)
            Hvf=fk*Hv(nck)+fp*Hv(i)
            dpf=(p(nck)-p(i))/dn(k,i)
          endif
          flux(k,i)=hk(k,i)*ds(k,i)*(fny(k,i)*Hvf-apuvolpf*hk(k,i)*dpf)
        endif
        flux(jcn,nck)=-flux(k,i)
      enddo
    enddo
!$OMP END DO
!$OMP DO PRIVATE(i,j,k,nck,jcn,fk,fp,apuvolpf,Huf,Hvf,dpf,dpfi)
    do i=1,ncellpoly
      do j=1,nxyface(i) !No Repeat of cell faces
        k=kxyface(j,i)     
        nck=cell2cell(k,i) !Forward connectivity 
        jcn=llec2llec(k,i) !Backward connectivity
        if(iwet(i)*iwet(nck)==0)then !Closed Boundaries
          flux(k,i)=0.0
        elseif(nck>ncells)then !Open Boundaries
          flux(k,i)=hk(k,i)*ds(k,i)*((fnx(k,i)*Hu(i)+fny(k,i)*Hv(i)) &
             -apuareap(i)*h(i)*(p(nck)-p(i))/dn(k,i))
        !!elseif(maxval(cell2cell(1:ncface(i),i))>ncells)then
        !!  fk=fintp(k,i); fp=1.0-fk
        !!  !apuvolpf=fk*apuareap(nck)+fp*apuareap(i)
        !!  apuvolpf=min(apuareap(nck),apuareap(i))
        !!  Huf=fk*Hu(nck)+fp*Hu(i)
        !!  Hvf=fk*Hv(nck)+fp*Hv(i)
        !!  flux(k,i)=hk(k,i)*ds(k,i)*((fnx(k,i)*Huf+fny(k,i)*Hvf) &
        !!     -apuvolpf*hk(k,i)*(p(nck)-p(i))/dn(k,i))  
        else
          fk=fintp(k,i); fp=1.0-fk          
          !!Approach 1
          !!apuvolpf=fk*(apuareap(nck)+rpx(jcn,nck)*dapuareapx(nck)+rpy(jcn,nck)*dapuareapy(nck)) &
          !!        +fp*(apuareap(i)+rpx(k,i)*dapuareapx(i)+rpy(k,i)*dapuareapy(i))
          !!apuvolpf=min(apuareap(nck),apuareap(i)) !more stable
          !apuvolpf=fk*apuareap(nck)+fp*apuareap(i)
          !Huf=fk*(Hu(nck)+rpx(jcn,nck)*dHux(nck)+rpy(jcn,nck)*dHuy(nck)) &
          !   +fp*(Hu(i)+rpx(k,i)*dHux(i)+rpy(k,i)*dHuy(i))
          !Hvf=fk*(Hv(nck)+rpx(jcn,nck)*dHvx(nck)+rpy(jcn,nck)*dHvy(nck)) &
          !   +fp*(Hv(i)+rpx(k,i)*dHvx(i)+rpy(k,i)*dHvy(i))
          !dpf=(p(nck)+rpx(jcn,nck)*dpx(nck)+rpy(jcn,nck)*dpy(nck) &
          !    -p(i)-rpx(k,i)*dpx(i)-rpy(k,i)*dpy(i))/dn(k,i)
          !flux(k,i)=hk(k,i)*ds(k,i)*((fnx(k,i)*Huf+fny(k,i)*Hvf)-apuvolpf*hk(k,i)*dpf)
          
          !Approach 2 (original Rhie and Chow method, similar results to approach 1)
          apuvolpf=fk*(apuareap(nck)+rpx(jcn,nck)*dapuareapx(nck)+rpy(jcn,nck)*dapuareapy(nck)) &
                  +fp*(apuareap(i)+rpx(k,i)*dapuareapx(i)+rpy(k,i)*dapuareapy(i))
          !apuvolpf=min(apuareap(nck),apuareap(i)) !more stable
          !apuvolpf=fk*apuareap(nck)+fp*apuareap(i)
          Huf=fk*(u(nck)+rpx(jcn,nck)*dux(nck)+rpy(jcn,nck)*duy(nck)) &
             +fp*(u(i)+rpx(k,i)*dux(i)+rpy(k,i)*duy(i))
          Hvf=fk*(v(nck)+rpx(jcn,nck)*dvx(nck)+rpy(jcn,nck)*dvy(nck)) &
             +fp*(v(i)+rpx(k,i)*dvx(i)+rpy(k,i)*dvy(i))
          dpf=(p(nck)+rpx(jcn,nck)*dpx(nck)+rpy(jcn,nck)*dpy(nck) &
              -p(i)-rpx(k,i)*dpx(i)-rpy(k,i)*dpy(i))/dn(k,i)
          !dpfi=fk*(fnx(k,i)*dpx(nck)+fny(k,i)*dpy(nck))&
          !     +fp*(fnx(k,i)*dpx(i)+fny(k,i)*dpy(i))
          dpfi=fnx(k,i)*(fk*dpx(nck)+fp*dpx(i)) &
              +fny(k,i)*(fk*dpy(nck)+fp*dpy(i))          
          flux(k,i)=hk(k,i)*ds(k,i)*((fnx(k,i)*Huf+fny(k,i)*Hvf) &
             +apuvolpf*hk(k,i)*(dpfi-dpf))
        endif
        flux(jcn,nck)=-flux(k,i)
      enddo      
    enddo
!$OMP END DO
!$OMP END PARALLEL

    return
    end subroutine flow_imp_flux_skew_rc    
    
!***********************************************************************
    subroutine flow_imp_flux_noskew
! Intercell flux calculation for implicit temporal scheme
! using a Rhie and Chow type momentum interpolation scheme
!
! written by Weiming Wu, NCCHE, Oct. 2008
! Modified by Alex Sanchez, USACE, 2013
!***********************************************************************
    use size_def
    use geo_def
    use der_def
    use der_lib, only: der_grad_eval,der_gradvec_eval
    use interp_def
    use flow_def
    use comp_lib
    use comvarbl, only: relax,skewcor
    use const_def, only: small
    use prec_def
    implicit none
    integer :: i,j,k,nck
    real(ikind) :: fk,fp,apuvolpf,Huf,Hvf
    
!$OMP PARALLEL DO PRIVATE(i,j,k,nck,fk,fp,apuvolpf,Huf,Hvf)
      do i=1,ncells
        do j=1,nxyface(i) !east and west sides
          k=kxyface(j,i)    
          nck=cell2cell(k,i) !Forward connectivity     
          if(iwet(i)*iwet(nck)==0)then !Closed Boundaries
            flux(k,i)=0.0
          elseif(nck>ncells)then !Open Boundaries
            flux(k,i)=hk(k,i)*ds(k,i)*((fnx(k,i)*Hu(i)+fny(k,i)*Hv(i)) &
                -apuareap(i)*hk(k,i)*(p(nck)-p(i))/dn(k,i))
            !!flux(k,i)=ds(k,i)*((fnx(k,i)*Hu(i)+fny(k,i)*Hv(i)) &
            !!    -apuareap(i)*hk(k,i)*hk(k,i)*(p(nck)-p(i))/dn(k,i))
          else
            fk=fintp(k,i); fp=1.0-fk
            !Standard approach
            apuvolpf=fk*apuareap(nck)+fp*apuareap(i)
            Huf=fk*Hu(nck)+fp*Hu(i)
            Hvf=fk*Hv(nck)+fp*Hv(i)
            flux(k,i)=hk(k,i)*ds(k,i)*((fnx(k,i)*Huf+fny(k,i)*Hvf) &
               -apuvolpf*hk(k,i)*(p(nck)-p(i))/dn(k,i))
          endif
          flux(llec2llec(k,i),nck)=-flux(k,i)
        enddo
      enddo
!$OMP END PARALLEL DO

    return
    end subroutine flow_imp_flux_noskew
    
!***********************************************************************
    subroutine flow_imp_flux_noskew_rc
! Intercell flux calculation for implicit temporal scheme
! using a Rhie and Chow momentum interpolation scheme
!
! written by Weiming Wu, NCCHE, Oct. 2008
! Modified by Alex Sanchez, USACE, 2013
!***********************************************************************
    use size_def
    use geo_def
    use der_def
    use der_lib, only: der_grad_eval,der_gradvec_eval
    use interp_def
    use flow_def
    use comp_lib
    use comvarbl, only: relax,skewcor
    use const_def, only: small
    use prec_def
    implicit none
    integer :: i,j,k,nck
    real(ikind) :: fk,fp,apuvolpf,Huf,Hvf,dpf
    
!$OMP PARALLEL DO PRIVATE(i,j,k,nck,fk,fp,apuvolpf,Huf,Hvf,dpf)
      do i=1,ncells
        do j=1,nxyface(i) !east and west sides
          k=kxyface(j,i)    
          nck=cell2cell(k,i) !Forward connectivity     
          if(iwet(i)*iwet(nck)==0)then !Closed Boundaries
            flux(k,i)=0.0
          elseif(nck>ncells)then !Open Boundaries
            flux(k,i)=hk(k,i)*ds(k,i)*((fnx(k,i)*u(i)+fny(k,i)*v(i)) &
                +apuareap(i)*hk(k,i)*((fnx(k,i)*dpx(i)+fny(k,i)*dpy(i))-p(nck)-p(i))/dn(k,i))
            !!flux(k,i)=ds(k,i)*((fnx(k,i)*Hu(i)+fny(k,i)*Hv(i)) &
            !!    -apuareap(i)*hk(k,i)*hk(k,i)*(p(nck)-p(i))/dn(k,i))
          else
            fk=fintp(k,i); fp=1.0-fk
            !Approach 2 (original Rhie and Chow method, similar results to approach 1)
            apuvolpf=fk*apuareap(nck)+fp*apuareap(i)
            Huf=fk*u(nck)+fp*u(i)
            Hvf=fk*v(nck)+fp*v(i)
            dpf=fk*(fnx(k,i)*dpx(nck)+fny(k,i)*dpy(nck))&
               +fp*(fnx(k,i)*dpx(i)+fny(k,i)*dpy(i))
            flux(k,i)=hk(k,i)*ds(k,i)*((fnx(k,i)*Huf+fny(k,i)*Hvf) &
               +apuvolpf*hk(k,i)*(dpf-(p(nck)-p(i))/dn(k,i)))
          endif
          flux(llec2llec(k,i),nck)=-flux(k,i)
        enddo
      enddo
!$OMP END PARALLEL DO

    return
    end subroutine flow_imp_flux_noskew_rc    