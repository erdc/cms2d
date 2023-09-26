!==========================================================================
! Computational subroutines
!==========================================================================

!***************************************************************************
    subroutine defcorpara(phi,ssphi)
! Calculates parallel components of deferred corrections 
! The parallel component is the skewness correction.
!
! Corrections are calculated at the cell-faces in order to avoid
! a race condition when using OpenMP. The corrections are then summed up
! for each cell and added to the source term.
! This subroutine is used when the gradients for a variable are
! not available and therefore the gradients are calculated internally.
! Since the gradients at regular/simple cells are not used, they are not
! calculated at these cells in order to save time. 
!
! Version: 2.0
! Last updated: November 6, 2012
! Author: by Alex Sanchez, USACE-CHL
!***************************************************************************    
    use size_def, only: ncellsD
    use der_def,  only: gow,nlim
    use der_lib,  only: der_grad_eval_recon
    use prec_def, only: ikind
    implicit none
    !Input/Output
    real(ikind),intent(in) :: phi(ncellsD)
    real(ikind),intent(inout) :: ssphi(ncellsD)
    !Internal Variables
    real(ikind) :: dphix(ncellsD),dphiy(ncellsD)

    call der_grad_eval_recon(gow,nlim,phi,dphix,dphiy) !Only computes gradients where needed
    call defcorparagrad(dphix,dphiy,ssphi)    
    
    return
    end subroutine defcorpara
    
!***************************************************************************
    subroutine defcorparagrad(dphix,dphiy,ssphi)
! Calculates parallel components of deferred corrections 
! The parallel component is the skewness correction.
!
! Corrections are calculated at the cell-faces in order to avoid
! a race condition when using OpenMP. The corrections are then summed up
! for each cell and added to the source term.
!
! Version: 2.0
! Last updated: November 6, 2012
! Author: by Alex Sanchez, USACE-CHL
!***************************************************************************    
    use size_def, only: ncelljoint, ncellpoly, ncells, ncellsd, nmaxfaces
    use flow_def, only: iwet,acoef,flux
    use geo_def,  only: idcelljoint, ncface, nxface, kxface, cell2cell, llec2llec, rpx, rpy, nyface, kyface, nxyface, kxyface
    use prec_def, only: ikind
    
    implicit none
    !Input/Output
    real(ikind),intent(in),dimension(ncellsD) :: dphix,dphiy
    real(ikind),intent(inout),dimension(ncellsD) :: ssphi
    !Internal Variables
    integer :: i,j,k,ii,nck,jcn
    real(ikind),dimension(nmaxfaces,ncellsD) :: fskew
    real(ikind) :: rp,rn
    
!$OMP PARALLEL
!$OMP DO PRIVATE(i,ii)     
    do ii=1,ncelljoint
      i=idcelljoint(ii)
      fskew(1:ncface(i),i)=0.0
    enddo
!$OMP END DO      
!$OMP DO PRIVATE(i)   
    do i=1,ncellpoly
      fskew(1:ncface(i),i)=0.0
    enddo
!$OMP END DO
!$OMP DO PRIVATE(i,j,ii,k,nck,jcn)     
    do ii=1,ncelljoint
      i=idcelljoint(ii)
      if(iwet(i)==0) cycle
      
      !X-Direction
      do j=1,nxface(i) !No repeat cell faces
        k=kxface(j,i)
        nck=cell2cell(k,i) !Forwards connectivity
        if(nck>ncells .or. iwet(nck)==0) cycle
        jcn=llec2llec(k,i) !Backwards connectivity
        if(abs(rpy(k,i))>1.0e-6)then !coarse to fine
          fskew(k,i)=(acoef(k,i)+flux(k,i))*rpy(k,i)*dphiy(i)
          fskew(jcn,nck)=-fskew(k,i)
        elseif(abs(rpy(jcn,nck))>1.0e-6)then !fine to coarse
          fskew(k,i)=-acoef(k,i)*rpy(jcn,nck)*dphiy(nck)
          fskew(jcn,nck)=-fskew(k,i)
        endif        
      enddo
      
      !Y-direction
      do j=1,nyface(i) !No repeat cell faces
        k=kyface(j,i)
        nck=cell2cell(k,i) !Forwards connectivity
        if(nck>ncells .or. iwet(nck)==0) cycle
        jcn=llec2llec(k,i) !Backwards connectivity
        if(abs(rpx(k,i))>1.0e-6)then !coarse to fine
          fskew(k,i)=(acoef(k,i)+flux(k,i))*rpx(k,i)*dphix(i) 
          fskew(jcn,nck)=-fskew(k,i)
        elseif(abs(rpx(jcn,nck))>1.0e-6)then !fine to coarse
          fskew(k,i)=-acoef(k,i)*rpx(jcn,nck)*dphix(nck)
          fskew(jcn,nck)=-fskew(k,i)
        endif        
      enddo
    enddo
!$OMP END DO

!$OMP DO PRIVATE(i,j,k,nck,jcn,rp,rn)
    do i=1,ncellpoly
      if(iwet(i)==0) cycle
      do j=1,nxyface(i) !No repeat cell faces
        k=kxyface(j,i)            
        nck=cell2cell(k,i) !Forwards connectivity
        if(nck>ncells .or. iwet(nck)==0) cycle
        jcn=llec2llec(k,i) !Backwards connectivity
        rp=rpx(k,i)*dphix(i)+rpy(k,i)*dphiy(i)
        rn=rpx(jcn,nck)*dphix(nck)+rpy(jcn,nck)*dphiy(nck)
        fskew(k,i)=(acoef(k,i)+flux(k,i))*rp-acoef(k,i)*rn
        fskew(jcn,nck)=-fskew(k,i)
      enddo
    enddo
!$OMP END DO

!$OMP DO PRIVATE(i,ii)   
    do ii=1,ncelljoint
      i=idcelljoint(ii)
      ssphi(i)=ssphi(i)-sum(fskew(1:ncface(i),i))
    enddo
!$OMP END DO

!$OMP DO PRIVATE(i)   
    do i=1,ncellpoly
      ssphi(i)=ssphi(i)-sum(fskew(1:ncface(i),i))
    enddo
!$OMP END DO
!$OMP END PARALLEL

    return
    end subroutine defcorparagrad

!***************************************************************************
    subroutine defcorparagradvec(dudx,dudy,dvdx,dvdy,ssu,ssv)
! Calculates parallel components of deferred corrections 
! The parallel component is the skewness correction.
!
! Corrections are calculated at the cell-faces in order to avoid
! a race condition when using OpenMP. The corrections are then summed up
! for each cell and added to the source term.
!
! Version: 2.0
! Last updated: November 6, 2012
! Author: by Alex Sanchez, USACE-CHL
!***************************************************************************    
    use size_def, only: ncelljoint, ncellpoly, ncells, ncellsd, nmaxfaces
    use flow_def, only: iwet,acoef,flux
    use geo_def,  only: idcelljoint, ncface, nxface, nyface, kxface, kyface, cell2cell, llec2llec, rpx, rpy, nxyface, kxyface
    use prec_def, only: ikind
    
    implicit none
    !Input/Output
    real(ikind),intent(in) :: dudx(ncellsD),dudy(ncellsD),dvdx(ncellsD),dvdy(ncellsD)
    real(ikind),intent(inout) :: ssu(ncellsD),ssv(ncellsD)
    !Internal Variables
    integer :: i,j,k,ii,nck,jcn
    real(ikind),dimension(nmaxfaces,ncellsD) :: uskew,vskew
    real(ikind) :: val,rpu,rpv,rnu,rnv
    
!$OMP PARALLEL
!$OMP DO PRIVATE(i,ii)
    do ii=1,ncelljoint
      i=idcelljoint(ii)
      uskew(1:ncface(i),i)=0.0
      vskew(1:ncface(i),i)=0.0
    enddo
!$OMP END DO      
!$OMP DO PRIVATE(i)   
    do i=1,ncellpoly
      uskew(1:ncface(i),i)=0.0
      vskew(1:ncface(i),i)=0.0
    enddo
!$OMP END DO
!$OMP DO PRIVATE(i,j,ii,k,nck,jcn,val)
    do ii=1,ncelljoint
      i=idcelljoint(ii)
      if(iwet(i)==0) cycle
      !X-Direction
      do j=1,nxface(i) !No repeat cell faces
        k=kxface(j,i)
        nck=cell2cell(k,i) !Forwards connectivity
        if(nck>ncells .or. iwet(nck)==0) cycle
        jcn=llec2llec(k,i) !Backwards connectivity
        if(abs(rpy(k,i))>1.0e-6)then !coarse to fine
          val=(acoef(k,i)+flux(k,i))*rpy(k,i)
          uskew(k,i)=val*dudy(i); vskew(k,i)=val*dvdy(i)
          uskew(jcn,nck)=-uskew(k,i); vskew(jcn,nck)=-vskew(k,i)
        elseif(abs(rpy(jcn,nck))>1.0e-6)then !fine to coarse
          val=-acoef(k,i)*rpy(jcn,nck)
          uskew(k,i)=val*dudy(nck); vskew(k,i)=val*dvdy(nck)
          uskew(jcn,nck)=-uskew(k,i); vskew(jcn,nck)=-vskew(k,i)
        endif
      enddo
      !Y-Direction
      do j=1,nyface(i) !No repeat cell faces
        k=kyface(j,i)
        nck=cell2cell(k,i) !Forwards connectivity
        if(nck>ncells .or. iwet(nck)==0) cycle
        jcn=llec2llec(k,i) !Backwards connectivity
        if(abs(rpx(k,i))>1.0e-6)then !coarse to fine
          val=(acoef(k,i)+flux(k,i))*rpx(k,i)
          uskew(k,i)=val*dudx(i);     vskew(k,i)=val*dvdx(i)
          uskew(jcn,nck)=-uskew(k,i); vskew(jcn,nck)=-vskew(k,i)
        elseif(abs(rpx(jcn,nck))>1.0e-6)then !fine to coarse
          val=-acoef(k,i)*rpx(jcn,nck)
          uskew(k,i)=val*dudx(nck);   vskew(k,i)=val*dvdx(nck)
          uskew(jcn,nck)=-uskew(k,i); vskew(jcn,nck)=-vskew(k,i)
        endif
      enddo
    enddo
!$OMP END DO
!$OMP DO PRIVATE(i,j,k,nck,jcn,rpu,rnu,rpv,rnv)
    do i=1,ncellpoly
      if(iwet(i)==0) cycle  
      if(maxval(cell2cell(1:ncface(i),i))>ncells) cycle  
      do j=1,nxyface(i) !No repeat cell faces
        k=kxyface(j,i)
        nck=cell2cell(k,i) !Forwards connectivity
        if(nck>ncells .or. iwet(nck)==0) cycle
        jcn=llec2llec(k,i) !Backwards connectivity
        !Linear reconstruction corrections
        rpu=rpx(k,i)*dudx(i)+rpy(k,i)*dudy(i)
        rnu=rpx(jcn,nck)*dudx(nck)+rpy(jcn,nck)*dudy(nck)
        rpv=rpx(k,i)*dvdx(i)+rpy(k,i)*dvdy(i)
        rnv=rpx(jcn,nck)*dvdx(nck)+rpy(jcn,nck)*dvdy(nck)
        !Skewness corrections
        uskew(k,i)=(acoef(k,i)+flux(k,i))*rpu-acoef(k,i)*rnu
        vskew(k,i)=(acoef(k,i)+flux(k,i))*rpv-acoef(k,i)*rnv
        uskew(jcn,nck)=-uskew(k,i); vskew(jcn,nck)=-vskew(k,i)
      enddo
    enddo
!$OMP END DO
!$OMP DO PRIVATE(i,ii)   
    do ii=1,ncelljoint
      i=idcelljoint(ii)
      ssu(i)=ssu(i)-sum(uskew(1:ncface(i),i))
      ssv(i)=ssv(i)-sum(vskew(1:ncface(i),i))
    enddo
!$OMP END DO
!$OMP DO PRIVATE(i)   
    do i=1,ncellpoly
      ssu(i)=ssu(i)-sum(uskew(1:ncface(i),i))
      ssv(i)=ssv(i)-sum(vskew(1:ncface(i),i))
    enddo
!$OMP END DO
!$OMP END PARALLEL

    return
    end subroutine defcorparagradvec
    
!***********************************************************************
    subroutine defcorhlpa(phi,ssphi)
! Deferred correction (anti-diffusion) for HLPA scheme
! for nonuniform Cartesian grids
!
! Version: 3.0
! Last updated: January 15, 2014
! Author: by Alex Sanchez, USACE-CHL
!***********************************************************************
    use size_def, only: ncells, ncellsd, nmaxfaces
    use geo_def,  only: ncface, nxyface, kxyface, cell2cell, llec2llec, kkface
    use flow_def, only: iwet,flux
    use der_def,  only: nder,nlim
    use prec_def, only: ikind
    
    implicit none
    !Input/Output
    real(ikind),intent(in),dimension(ncellsD) :: phi
    real(ikind),intent(inout),dimension(ncellsD) :: ssphi
    !Internal Variables
    integer :: i,j,k,nck,jcn,nckk
    real(ikind),dimension(nmaxfaces,ncellsD) :: fnorm
    real(ikind) :: phinorm,phic,phid,phiu

!$OMP PARALLEL    
!--- Initialize --------------------
!$OMP DO PRIVATE(i)
    do i=1,ncells
      fnorm(1:ncface(i),i)=0.0
    enddo
!$OMP END DO
!$OMP DO PRIVATE(i,j,k,nck,nckk,jcn,phinorm,phic,phid,phiu)
    do i=1,ncells
      if(iwet(i)==0) cycle  
      do j=1,nxyface(i) !No repeat cell face count
        k=kxyface(j,i)   
        nck=cell2cell(k,i) !Forwards connectivity  
        if(nck>ncells .or. iwet(nck)==0) cycle
        jcn=llec2llec(k,i) !Backwards connectivity
        if(flux(k,i)>=0.0)then !Outflow
          nckk=cell2cell(kkface(k),i)  
          phic=phi(i)
          phid=phi(nck)
          phiu=phi(nckk)
        else              !Inflow
          nckk=cell2cell(k,nck)
          phic=phi(nck)
          phid=phi(i)
          phiu=phi(nckk)
        endif
        if(abs(phid-phiu)>1.0e-20)then !Avoid divide by zero
          phinorm=(phic-phiu)/(phid-phiu) !Normalized
          if(phinorm>0.0 .and. phinorm<1.0)then
            fnorm(k,i)=flux(k,i)*(phid-phic)*phinorm  !Deferred normal flux
            fnorm(jcn,nck)=-fnorm(k,i)
          endif
        endif
      enddo
    enddo
!$OMP END DO
!$OMP DO PRIVATE(i)
    do i=1,ncells
      ssphi(i)=ssphi(i)-sum(fnorm(1:ncface(i),i))
    enddo
!$OMP END DO
!$OMP END PARALLEL

    return
    end subroutine defcorhlpa
    
!***********************************************************************
    subroutine defcorhlpavec(u,v,ssu,ssv)
! Deferred correction (anti-diffusion) for HLPA scheme
! for nonuniform Cartesian grids
!
! Version: 3.0
! Last updated: January 15, 2014
! Author: by Alex Sanchez, USACE-CHL
!***********************************************************************
    use size_def, only: ncells, ncellsd, nmaxfaces
    use geo_def,  only: cell2cell, llec2llec, ncface, nxyface, kxyface, kkface
    use flow_def, only: iwet,flux
    use der_def,  only: nder,nlim
    use prec_def, only: ikind
    
    implicit none
    !Input/Output
    real(ikind),intent(in),dimension(ncellsD) :: u,v
    real(ikind),intent(inout),dimension(ncellsD) :: ssu,ssv
    !Internal Variables
    integer :: i,j,k,nck,jcn,nckk
    real(ikind),dimension(nmaxfaces,ncellsD) :: ufnorm,vfnorm
    real(ikind) :: unorm,uc,ud,uu
    real(ikind) :: vnorm,vc,vd,vu

!$OMP PARALLEL    
!--- Initialize --------------------
!$OMP DO PRIVATE(i)
    do i=1,ncells
      ufnorm(1:ncface(i),i)=0.0
      vfnorm(1:ncface(i),i)=0.0
    enddo
!$OMP END DO
!$OMP DO PRIVATE(i,j,k,nck,nckk,jcn,unorm,uc,ud,uu,vnorm,vc,vd,vu)
    do i=1,ncells
      if(iwet(i)==0) cycle  
      do j=1,nxyface(i) !No repeat cell face count
        k=kxyface(j,i)   
        nck=cell2cell(k,i) !Forwards connectivity  
        if(nck>ncells .or. iwet(nck)==0) cycle
        jcn=llec2llec(k,i) !Backwards connectivity
        if(flux(k,i)>=0.0)then !Outflow
          nckk=cell2cell(kkface(k),i)
          if(nckk>ncells .or. iwet(nckk)==0) cycle
          uc=u(i); ud=u(nck); uu=u(nckk)
          vc=v(i); vd=v(nck); vu=v(nckk)
        else              !Inflow
          nckk=cell2cell(k,nck)
          if(nckk>ncells .or. iwet(nckk)==0) cycle
          uc=u(nck); ud=u(i); uu=u(nckk)
          vc=v(nck); vd=v(i); vu=v(nckk)
        endif
        if(abs(ud-uu)>1.0e-20)then !Avoid divide by zero
          unorm=(uc-uu)/(ud-uu) !Normalized
          if(unorm>0.0 .and. unorm<1.0)then
            ufnorm(k,i)=flux(k,i)*(ud-uc)*unorm  !Deferred normal flux
            ufnorm(jcn,nck)=-ufnorm(k,i)
          endif
        endif
        if(abs(vd-vu)>1.0e-20)then !Avoid divide by zero
          vnorm=(vc-vu)/(vd-vu) !Normalized
          if(vnorm>0.0 .and. vnorm<1.0)then
            vfnorm(k,i)=flux(k,i)*(vd-vc)*vnorm  !Deferred normal flux
            vfnorm(jcn,nck)=-vfnorm(k,i)
          endif
        endif
      enddo
    enddo
!$OMP END DO
!$OMP DO PRIVATE(i)
    do i=1,ncells
      ssu(i)=ssu(i)-sum(ufnorm(1:ncface(i),i))
      ssv(i)=ssv(i)-sum(vfnorm(1:ncface(i),i))
    enddo
!$OMP END DO
!$OMP END PARALLEL

    return
    end subroutine defcorhlpavec
    
!***********************************************************************
    subroutine defcorgamma(schmdefcor,phi,ssphi)
! Deferred correction (anti-diffusion) for Gamma-family schemes
! for nonuniform Cartesian grids
!
! Version: 3.0
! Last updated: January 15, 2014
! Author: by Alex Sanchez, USACE-CHL
!***********************************************************************
    use size_def, only: ncellsd, ncells, nmaxfaces
    use geo_def,  only: cell2cell, llec2llec, ncface, nxyface, kxyface, kkface
    use flow_def, only: iwet,flux
    use der_def,  only: nder,nlim
    use prec_def, only: ikind
    use interp_def, only: fintp

    implicit none
    !Input/Output
    real(ikind),intent(in),dimension(ncellsD) :: phi
    real(ikind),intent(inout),dimension(ncellsD) :: ssphi
    !Internal Variables
    integer :: i,j,k,nck,nckk,jcn
    real(ikind),dimension(nmaxfaces,ncellsD) :: fnorm
    real(ikind) :: phinorm,phic,phid,phiu,phif
    
    interface
      function schmdefcor(phin)
        use prec_def
        implicit none
        real(ikind),intent(in) :: phin
        real(ikind) :: schmdefcor
      end function
    endinterface

!$OMP PARALLEL
!$OMP DO PRIVATE(i)
    do i=1,ncells
      fnorm(1:ncface(i),i)=0.0
    enddo
!$OMP END DO
!$OMP DO PRIVATE(i,j,k,nck,nckk,jcn,phinorm,phif,phic,phid,phiu)
      do i=1,ncells
        if(iwet(i)==0) cycle
        do j=1,nxyface(i) !No repeat cell face count
          k=kxyface(j,i)
          nck=cell2cell(k,i) !Forwards connectivity
          if(nck>ncells .or. iwet(nck)==0) cycle
          jcn=llec2llec(k,i) !Backwards connectivity
          if(flux(k,i)>=0.0)then !Outflow
            nckk=cell2cell(kkface(k),i)
            if(nckk>ncells .or. iwet(nckk)==0) cycle
            phic=phi(i); phid=phi(nck); phiu=phi(nckk)
          else              !Inflow
            nckk=cell2cell(k,nck)
            if(nckk>ncells .or. iwet(nckk)==0) cycle
            phic=phi(nck); phid=phi(i); phiu=phi(nckk)
          endif
          if(abs(phid-phiu)>1.0e-20)then !Avoid divide by zero
            phinorm=(phic-phiu)/(phid-phiu) !Normalized
            if(phinorm>0.0 .and. phinorm<1.0)then
              phif=fintp(k,i)*phi(nck)+(1.0-fintp(k,i))*phi(i) !Central Difference
              fnorm(k,i)=flux(k,i)*schmdefcor(phinorm)*(phif-phic) !Deferred normal flux
              fnorm(jcn,nck)=-fnorm(k,i)
            endif
          endif
        enddo
      enddo
!$OMP END DO
!$OMP DO PRIVATE(i)
    do i=1,ncells
      ssphi(i)=ssphi(i)-sum(fnorm(1:ncface(i),i))
    enddo
!$OMP END DO
!$OMP END PARALLEL
    
    return
    end subroutine defcorgamma
    
!***********************************************************************
    subroutine defcorgammavec(schmdefcor,u,v,ssu,ssv)
! Deferred correction (anti-diffusion) for Gamma-family schemes
! for nonuniform Cartesian grids
!
! Version: 3.0
! Last updated: January 15, 2014
! Author: by Alex Sanchez, USACE-CHL
!***********************************************************************
    use size_def, only: ncells, ncellsd, nmaxfaces
    use geo_def,  only: ncface, nxyface, kxyface, kkface, cell2cell, llec2llec
    use flow_def, only: iwet,flux
    use der_def,  only: nder,nlim
    use prec_def, only: ikind
    use interp_def, only: fintp

    implicit none
    !Input/Output
    real(ikind),intent(in),dimension(ncellsD) :: u,v
    real(ikind),intent(inout),dimension(ncellsD) :: ssu,ssv
    !Internal Variables
    integer :: i,j,k,nck,nckk,jcn
    real(ikind),dimension(nmaxfaces,ncellsD) :: ufnorm,vfnorm
    real(ikind) :: unorm,uc,ud,uu,uf
    real(ikind) :: vnorm,vc,vd,vu,vf
    
    interface
      function schmdefcor(phin)
        use prec_def
        implicit none
        real(ikind),intent(in) :: phin
        real(ikind) :: schmdefcor
      end function
    endinterface

!$OMP PARALLEL
!$OMP DO PRIVATE(i)
    do i=1,ncells
      ufnorm(1:ncface(i),i)=0.0
      vfnorm(1:ncface(i),i)=0.0
    enddo
!$OMP END DO
!$OMP DO PRIVATE(i,j,k,nck,nckk,jcn,unorm,uc,ud,uu,uf,vnorm,vc,vd,vu,vf)
      do i=1,ncells
        if(iwet(i)==0) cycle
        do j=1,nxyface(i) !No repeat cell face count
          k=kxyface(j,i)
          nck=cell2cell(k,i) !Forwards connectivity
          if(nck>ncells .or. iwet(nck)==0) cycle
          jcn=llec2llec(k,i) !Backwards connectivity
          if(flux(k,i)>=0.0)then !Outflow
            nckk=cell2cell(kkface(k),i)
            uc=u(i); ud=u(nck); uu=u(nckk)
            vc=v(i); vd=v(nck); vu=v(nckk)
          else              !Inflow
            nckk=cell2cell(k,nck)
            uc=u(nck); ud=u(i); uu=u(nckk)
            vc=v(nck); vd=v(i); vu=v(nckk)
          endif
          if(abs(ud-uu)>1.0e-20)then !Avoid divide by zero
            unorm=(uc-uu)/(ud-uu) !Normalized
            if(unorm>0.0 .and. unorm<1.0)then
              uf=fintp(k,i)*u(nck)+(1.0-fintp(k,i))*u(i) !Central Difference
              ufnorm(k,i)=flux(k,i)*schmdefcor(unorm)*(uf-uc) !Deferred normal flux
              ufnorm(jcn,nck)=-ufnorm(k,i)
            endif
          endif
          if(abs(vd-vu)>1.0e-20)then !Avoid divide by zero
            vnorm=(vc-vu)/(vd-vu) !Normalized
            if(vnorm>0.0 .and. vnorm<1.0)then
              vf=fintp(k,i)*v(nck)+(1.0-fintp(k,i))*v(i) !Central Difference
              vfnorm(k,i)=flux(k,i)*schmdefcor(vnorm)*(vf-vc) !Deferred normal flux
              vfnorm(jcn,nck)=-vfnorm(k,i)
            endif
          endif
        enddo
      enddo
!$OMP END DO
!$OMP DO PRIVATE(i)
    do i=1,ncells
      ssu(i)=ssu(i)-sum(ufnorm(1:ncface(i),i))
      ssv(i)=ssv(i)-sum(vfnorm(1:ncface(i),i))
    enddo
!$OMP END DO
!$OMP END PARALLEL
    
    return
    end subroutine defcorgammavec    
    
!***********************************************************************
    subroutine defcorhlpagrad(phi,dphix,dphiy,ssphi)
! Deferred correction (anti-diffusion) for hlpa scheme.
!
! Uses input gradients for linear reconstructions and the calculation 
! of the normalized variable at joint cells.
!
! Version: 2.0
! Last updated: November 6, 2012
! Author: by Alex Sanchez, USACE-CHL
!***********************************************************************
    use size_def, only: ncells, ncellsd, nmaxfaces, ncellsimple, ncelljoint, ncellpoly
    use geo_def,  only: cell2cell, llec2llec, ncface, idcellsimple, nxface, nyface, kxface, kyface, dnx, dny, idcelljoint, rpx, rpy, nxyface, kxyface
    use flow_def, only: iwet,flux,acoef
    use prec_def, only: ikind
    
    implicit none
    !Input/Output
    real(ikind),intent(in),dimension(ncellsD) :: phi,dphix,dphiy
    real(ikind),intent(inout),dimension(ncellsD) :: ssphi
    !Internal Variables
    integer :: ii,i,j,k,nck,jcn    
    real(ikind),dimension(nmaxfaces,ncellsD) :: fnorm
    real(ikind) :: phinorm,phic,phid,dndphin,rp,rn,phip,phin
    
!$OMP PARALLEL
!--- Initialize -----------------------
!$OMP DO PRIVATE(i)
    do i=1,ncells
      fnorm(1:ncface(i),i)=0.0
    enddo
!$OMP END DO
!--- Telescoping Grids --------------------------------------
!$OMP DO PRIVATE(i,ii,k,nck,jcn,phinorm,phid,phic,dndphin)
    do ii=1,ncellsimple !Regular/simple cells
      i=idcellsimple(ii)
      if(iwet(i)==0) cycle
      !X-Direction
      do j=1,nxface(i) !No repeat cell face count
        k=kxface(j,i)   
        nck=cell2cell(k,i) !Forwards connectivity
        if(nck>ncells .or. iwet(nck)==0) cycle
        jcn=llec2llec(k,i) !Backwards connectivity
        !Upwinding and normal gradient
        if(flux(k,i)>0.0)then !Outflow
          dndphin=dnx(k,i)*dphix(i)
          phic=phi(i); phid=phi(nck)
        else !Inflow
          dndphin=dnx(jcn,nck)*dphix(nck)
          phid=phi(i); phic=phi(nck)
        endif
        if(abs(dndphin)>1.0e-10)then
          phinorm=1.0-0.5*(phid-phic)/dndphin  !Normalized
          if(phinorm>0.0 .and. phinorm<1.0)then
            fnorm(k,i)=flux(k,i)*(phid-phic)*phinorm  !Deferred normal flux
            fnorm(jcn,nck)=-fnorm(k,i)
          endif
        endif
      enddo
      !Y-Direction
      do j=1,nyface(i) !No repeat cell face count
        k=kyface(j,i)   
        nck=cell2cell(k,i) !Forwards connectivity
        if(nck>ncells .or. iwet(nck)==0) cycle
        jcn=llec2llec(k,i) !Backwards connectivity
        !Upwinding and normal gradient
        if(flux(k,i)>0.0)then !Outflow
          dndphin=dny(k,i)*dphiy(i)
          phic=phi(i); phid=phi(nck)
        else !Inflow
          dndphin=dny(jcn,nck)*dphiy(nck)
          phid=phi(i); phic=phi(nck)
        endif
        if(abs(dndphin)>1.0e-10)then !To avoid divide by zero
          phinorm=1.0-0.5*(phid-phic)/dndphin  !Normalized
          if(phinorm>0.0 .and. phinorm<1.0)then
            fnorm(k,i)=flux(k,i)*(phid-phic)*phinorm  !Deferred normal flux
            fnorm(jcn,nck)=-fnorm(k,i)
          endif
        endif
      enddo
     enddo
!$OMP END DO
!$OMP DO PRIVATE(i,ii,j,k,nck,jcn,phinorm,phip,phin,phic,phid,dndphin,rp,rn)
    do ii=1,ncelljoint !Joint Cells
      i=idcelljoint(ii)
      if(iwet(i)==0) cycle
      !X-Direction
      do j=1,nxface(i) !no repeat east and west sides
        k=kxface(j,i)    
        nck=cell2cell(k,i) !Forwards connectivity
        if(nck>ncells .or. iwet(nck)==0) cycle
        jcn=llec2llec(k,i) !Backwards connectivity
        !Cell-reconstruction and sckewness correction
        phip=phi(i); phin=phi(nck)
        if(abs(rpy(k,i))>1.0e-6)then !coarse to fine
          rp=rpy(k,i)*dphiy(i); phip=phip+rp
          fnorm(k,i)=(acoef(k,i)+flux(k,i))*rp !Deferred parallel flux         
        elseif(abs(rpy(jcn,nck))>1.0e-6)then !fine to coarse
          rn=rpy(jcn,nck)*dphiy(nck); phin=phin+rn
          fnorm(k,i)=-acoef(k,i)*rn !Deferred parallel flux
        endif
        !Upwinding and normal gradient
        if(flux(k,i)>0.0)then !Outflow
          dndphin=dnx(k,i)*dphix(i)
          phic=phip; phid=phin
        else !Inflow
          dndphin=dnx(jcn,nck)*dphix(nck)
          phid=phip; phic=phin
        endif
        if(abs(dndphin)>1.0e-10)then
          phinorm=1.0-0.5*(phid-phic)/dndphin  !Normalized
          if(phinorm>0.0 .and. phinorm<1.0)then
            fnorm(k,i)=fnorm(k,i)+flux(k,i)*(phid-phic)*phinorm  !Deferred normal flux
          endif
        endif
        fnorm(jcn,nck)=-fnorm(k,i)
      enddo
      !Y-Direction
      do j=1,nyface(i) !no repeat north and south sides
        k=kyface(j,i)    
        nck=cell2cell(k,i) !Forwards connectivity
        if(nck>ncells .or. iwet(nck)==0) cycle
        jcn=llec2llec(k,i) !Backwards connectivity
        !Cell-reconstruction and sckewness correction
        phip=phi(i); phin=phi(nck)
        if(abs(rpx(k,i))>1.0e-6)then !coarse to fine
          rp=rpx(k,i)*dphix(i); phip=phip+rp
          fnorm(k,i)=(acoef(k,i)+flux(k,i))*rp !Deferred parallel flux
        elseif(abs(rpx(jcn,nck))>1.0e-6)then !fine to coarse
          rn=rpx(jcn,nck)*dphix(nck); phin=phin+rn
          fnorm(k,i)=-acoef(k,i)*rn !Deferred parallel flux
        endif
        !Upwinding and normal gradient
        if(flux(k,i)>0.0)then !Outflow
          dndphin=dny(k,i)*dphiy(i)
          phic=phip; phid=phin
        else !Inflow
          dndphin=dny(jcn,nck)*dphiy(nck)
          phid=phip; phic=phin
        endif
        if(abs(dndphin)>1.0e-10)then
          phinorm=1.0-0.5*(phid-phic)/dndphin  !Normalized
          if(phinorm>0.0 .and. phinorm<1.0)then
            fnorm(k,i)=fnorm(k,i)+flux(k,i)*(phid-phic)*phinorm  !Deferred normal flux
          endif
        endif
        fnorm(jcn,nck)=-fnorm(k,i)
      enddo
    enddo
!$OMP END DO

!--- Polyhedral Grids ------------------------------------------------------
!$OMP DO PRIVATE(i,j,k,nck,jcn,phinorm,phip,phin,phic,phid,dndphin,rp,rn)
    do i=1,ncellpoly !No index mapping required for polygonal cells (all cells are polygonal)
      if(iwet(i)==0) cycle
      do j=1,nxyface(i) !no repeat sides
        k=kxyface(j,i)    
        nck=cell2cell(k,i) !Forwards connectivity
        if(nck>ncells .or. iwet(nck)==0) cycle
        jcn=llec2llec(k,i) !Backwards connectivity
        !Cell reconstruction corrections
        rp=rpx(k,i)*dphix(i)+rpy(k,i)*dphiy(i)
        rn=rpx(jcn,nck)*dphix(nck)+rpy(jcn,nck)*dphiy(nck) 
        !Skewness corrections
        fnorm(k,i)=(acoef(k,i)+flux(k,i))*rp-acoef(k,i)*rn
        !Upwinding, linear reconstructions, and normal gradients
        if(flux(k,i)>0.0)then !Outflow
          dndphin=dnx(k,i)*dphix(i)+dny(k,i)*dphiy(i)
          phic=phi(i)+rp; phid=phi(nck)+rn 
        else !Inflow
          dndphin=dnx(jcn,nck)*dphix(nck)+dny(jcn,nck)*dphiy(nck)
          phid=phi(i)+rp; phic=phi(nck)+rn
        endif
        if(abs(dndphin)>1.0e-10)then
          phinorm=1.0-0.5*(phid-phic)/dndphin  !Normalized
          if(phinorm>0.0 .and. phinorm<1.0)then
            fnorm(k,i)=fnorm(k,i)+flux(k,i)*(phid-phic)*phinorm  !Deferred normal flux
          endif
        endif
        fnorm(jcn,nck)=-fnorm(k,i)
      enddo
    enddo
!$OMP END DO

!---Sum fluxes to get deferred correction --------------
!$OMP DO PRIVATE(i)
    do i=1,ncells
      ssphi(i)=ssphi(i)-sum(fnorm(1:ncface(i),i))
    enddo
!$OMP END DO    
!$OMP END PARALLEL  

    return
    end subroutine defcorhlpagrad

!***********************************************************************
    subroutine defcorgammagrad(schmdefcor,phi,dphix,dphiy,ssphi)
! Deferred correction (anti-diffusion) for Gamma-family schemes
!
! Version: 2.0
! Last updated: November 6, 2012
! Author: by Alex Sanchez, USACE-CHL
!***********************************************************************
    use size_def, only: ncells, ncellsd, nmaxfaces, ncellsimple, ncelljoint, ncellpoly
    use geo_def,  only: rpx, rpy, ncface, nxface, kxface, nyface, kyface, nxyface, kxyface, cell2cell, llec2llec, dnx, dny, idcellsimple, idcelljoint
    use flow_def, only: iwet,flux,acoef
    use prec_def, only: ikind
    use interp_def, only: fintp

    implicit none
    !Input/Output
    real(ikind),intent(in),dimension(ncellsD) :: phi,dphix,dphiy
    real(ikind),intent(inout),dimension(ncellsD) :: ssphi
    !Internal Variables
    integer :: ii,i,j,k,nck,jcn    
    real(ikind),dimension(nmaxfaces,ncellsD) :: fnorm
    real(ikind) :: phinorm,phip,phin,phic,phid,phif,dndphin,rp,rn

    interface
      function schmdefcor(phin)
        use prec_def
        implicit none
        real(ikind),intent(in) :: phin
        real(ikind) :: schmdefcor
      end function
    endinterface

!--- All other types of grids ------------------------
!$OMP PARALLEL
!$OMP DO PRIVATE(i)
    do i=1,ncells
      fnorm(1:ncface(i),i)=0.0
    enddo
!$OMP END DO
!--- Telescoping Grids --------------------------------------
!$OMP DO PRIVATE(i,ii,j,k,nck,jcn,phinorm,phic,phid,phif,dndphin)
    do ii=1,ncellsimple
      i=idcellsimple(ii)
      if(iwet(i)==0) cycle
      !X-Direction
      do j=1,nxface(i) !no repeat east and west sides
        k=kxface(j,i)    
        nck=cell2cell(k,i) !Forwards connectivity
        if(nck>ncells .or. iwet(nck)==0) cycle
        jcn=llec2llec(k,i) !Backwards connectivity
        !Upwinding and normal gradient
        if(flux(k,i)>0.0)then !Outflow
          dndphin=dnx(k,i)*dphix(i)
          phic=phi(i); phid=phi(nck)
        else !Inflow
          dndphin=dnx(jcn,nck)*dphix(nck)
          phic=phi(nck); phid=phi(i)
        endif
        if(abs(dndphin)>1.0e-10)then
          phinorm=1.0-0.5*(phid-phic)/dndphin  !Normalized
          if(phinorm>0.0 .and. phinorm<1.0)then
            phif=fintp(k,i)*phi(nck)+(1.0-fintp(k,i))*phi(i) !Central Difference
            fnorm(k,i)=flux(k,i)*schmdefcor(phinorm)*(phif-phic)   !Deferred normal flux
            fnorm(jcn,nck)=-fnorm(k,i)
          endif
        endif
      enddo
      !Y-Direction
      do j=1,nyface(i) !no repeat north and south faces
        k=kyface(j,i)
        nck=cell2cell(k,i) !Forwards connectivity
        if(nck>ncells .or. iwet(nck)==0) cycle
        jcn=llec2llec(k,i) !Backwards connectivity
        !Upwinding and normal gradient
        if(flux(k,i)>0.0)then !Outflow
          dndphin=dny(k,i)*dphiy(i)
          phic=phi(i); phid=phi(nck)
        else !Inflow
          dndphin=dny(jcn,nck)*dphiy(nck)
          phic=phi(nck); phid=phi(i)
        endif
        if(abs(dndphin)>1.0e-10)then
          phinorm=1.0-0.5*(phid-phic)/dndphin  !Normalized
          if(phinorm>0.0 .and. phinorm<1.0)then
            phif=fintp(k,i)*phi(nck)+(1.0-fintp(k,i))*phi(i) !Central Difference
            fnorm(k,i)=flux(k,i)*schmdefcor(phinorm)*(phif-phic)   !Deferred normal flux
            fnorm(jcn,nck)=-fnorm(k,i)
          endif
        endif
      enddo
    enddo
!$OMP END DO
!$OMP DO PRIVATE(i,ii,j,k,nck,jcn,phinorm,phip,phin,phic,phid,phif,dndphin,rp,rn)
    do ii=1,ncelljoint !Joint cells
      i=idcelljoint(ii)
      if(iwet(i)==0) cycle
      !X-Direction
      do j=1,nxface(i) !no repeat east and west sides
        k=kxface(j,i)    
        nck=cell2cell(k,i) !Forwards connectivity
        if(nck>ncells .or. iwet(nck)==0) cycle
        jcn=llec2llec(k,i) !Backwards connectivity
        !Linear cell reconstructions and sckewness correction
        phip=phi(i); phin=phi(nck)
        if(abs(rpy(k,i))>1.0e-6)then !coarse to fine
          rp=rpy(k,i)*dphiy(i); phip=phip+rp
          fnorm(k,i)=(acoef(k,i)+flux(k,i))*rp !Deferred parallel flux
        elseif(abs(rpy(jcn,nck))>1.0e-6)then !fine to coarse
          rn=rpy(jcn,nck)*dphiy(nck); phin=phin+rn
          fnorm(k,i)=-acoef(k,i)*rn !Deferred parallel flux
        endif
        !Upwinding and normal distance
        if(flux(k,i)>0.0)then !Outflow
          dndphin=dnx(k,i)*dphix(i)
          phic=phip; phid=phin
        else !Inflow
          dndphin=dnx(jcn,nck)*dphix(nck)
          phid=phip; phic=phin
        endif
        if(abs(dndphin)>1.0e-10)then
          phinorm=1.0-0.5*(phid-phic)/dndphin  !Normalized
          if(phinorm>0.0 .and. phinorm<1.0)then
            phif=fintp(k,i)*phin+(1.0-fintp(k,i))*phip !Central Difference
            fnorm(k,i)=fnorm(k,i)+flux(k,i)*schmdefcor(phinorm)*(phif-phic) !Deferred normal flux
          endif
        endif
        fnorm(jcn,nck)=-fnorm(k,i)
      enddo
      !Y-Direction
      do j=1,nyface(i) !no repeat north and south faces
        k=kyface(j,i)
        nck=cell2cell(k,i) !Forwards connectivity
        if(nck>ncells .or. iwet(nck)==0) cycle
        jcn=llec2llec(k,i) !Backwards connectivity
        !Linear cell reconstructions and sckewness correction
        phip=phi(i); phin=phi(nck)
        if(abs(rpx(k,i))>1.0e-6)then !coarse to fine
          rp=rpx(k,i)*dphix(i); phip=phip+rp
          fnorm(k,i)=(acoef(k,i)+flux(k,i))*rp !Deferred parallel flux
        elseif(abs(rpx(jcn,nck))>1.0e-6)then !fine to coarse
          rn=rpx(jcn,nck)*dphix(nck); phin=phin+rn
          fnorm(k,i)=-acoef(k,i)*rn !Deferred parallel flux
        endif
        !Upwinding and normal distances
        if(flux(k,i)>0.0)then !Outflow
          dndphin=dny(k,i)*dphiy(i)
          phic=phip; phid=phin
        else !Inflow
          dndphin=dny(jcn,nck)*dphiy(nck)
          phid=phip; phic=phin
        endif
        if(abs(dndphin)>1.0e-10)then
          phinorm=1.0-0.5*(phid-phic)/dndphin  !Normalized
          if(phinorm>0.0 .and. phinorm<1.0)then
            phif=fintp(k,i)*phin+(1.0-fintp(k,i))*phip !Central Difference
            fnorm(k,i)=fnorm(k,i)+flux(k,i)*schmdefcor(phinorm)*(phif-phic) !Deferred normal flux
          endif
        endif
        fnorm(jcn,nck)=-fnorm(k,i)
      enddo
    enddo
!$OMP END DO

    !Polygonal Meshes
!$OMP DO PRIVATE(i,j,k,nck,jcn,phinorm,phic,phid,phif,dndphin,rp,rn)
    do i=1,ncellpoly
      if(iwet(i)==0) cycle
      do j=1,nxyface(i) !no repeat east and west sides
        k=kxyface(j,i)    
        nck=cell2cell(k,i) !Forwards connectivity
        if(nck>ncells .or. iwet(nck)==0) cycle
        jcn=llec2llec(k,i) !Backwards connectivity
        !Cell reconstruction corrections
        rp=rpx(k,i)*dphix(i)+rpy(k,i)*dphiy(i)
        rn=rpx(jcn,nck)*dphix(nck)+rpy(jcn,nck)*dphiy(nck) 
        !Skewness corrections
        fnorm(k,i)=(acoef(k,i)+flux(k,i))*rp-acoef(k,i)*rn
        !Upwinding and linear reconstructions
        if(flux(k,i)>0.0)then !Outflow
          dndphin=dnx(k,i)*dphix(i)+dny(k,i)*dphiy(i)
          phic=phi(i)+rp; phid=phi(nck)+rn 
        else !Inflow
          dndphin=dnx(jcn,nck)*dphix(nck)+dny(jcn,nck)*dphiy(nck)
          phid=phi(i)+rp; phic=phi(nck)+rn
        endif
        if(abs(dndphin)>1.0e-10)then
          phinorm=1.0-0.5*(phid-phic)/dndphin  !Normalized
          if(phinorm>0.0 .and. phinorm<1.0)then
            phif=fintp(k,i)*phin+(1.0-fintp(k,i))*phip !Central Difference
            fnorm(k,i)=fnorm(k,i)+flux(k,i)*schmdefcor(phinorm)*(phif-phic) !Deferred normal flux  
          endif
        endif
        fnorm(jcn,nck)=-fnorm(k,i)
      enddo  
    enddo
!$OMP END DO
!$OMP DO PRIVATE(i)   
    do i=1,ncells
      ssphi(i)=ssphi(i)-sum(fnorm(1:ncface(i),i))
    enddo
!$OMP END DO    
!$OMP END PARALLEL
    
    return
    end subroutine defcorgammagrad
    
!***********************************************************************
    subroutine defcorhlpagradvecnew(u,v,dudx,dudy,dvdx,dvdy,ssu,ssv)
! Deferred correction (anti-diffusion) for HLPA scheme.
!
! Uses input gradients for linear reconstructions and the calculation 
! of the normalized variable at joint cells.
!
! Version: 2.0
! Last updated: November 6, 2012
! Author: by Alex Sanchez, USACE-CHL
!***********************************************************************
    use size_def, only: ncells, ncellsd, nmaxfaces
    use geo_def,  only: rpx, rpy, dnx, dny, ncface, kxyface, nxyface, cell2cell, llec2llec, cell2upwdcell
    use flow_def, only: iwet,flux,acoef
    use prec_def, only: ikind
    
    implicit none
    !Input/Output
    real(ikind),intent(in),dimension(ncellsD) :: u,v,dudx,dudy,dvdx,dvdy
    real(ikind),intent(inout),dimension(ncellsD) :: ssu,ssv
    !Internal Variables
    integer :: i,j,k,nck,nckk,jcn    
    real(ikind),dimension(nmaxfaces,ncellsD) :: ufnorm,vfnorm
    real(ikind) :: uvnorm,uc,vc,ud,vd,uu,vu,dndun,dndvn,rpu,rpv,rnu,rnv
    
!$OMP PARALLEL
!---- Initialize -----------------------------------------
!$OMP DO PRIVATE(i)
    do i=1,ncells
      ufnorm(1:ncface(i),i)=0.0
      vfnorm(1:ncface(i),i)=0.0
    enddo
!$OMP END DO
!--- All cells --------------------------------------------------------------------
!$OMP DO PRIVATE(i,j,k,nck,jcn,uvnorm,uc,vc,ud,vd,uu,vu,dndun,dndvn,rnu,rpu,rpv,rnv)
    do i=1,ncells
      if(iwet(i)==0) cycle
      do j=1,nxyface(i) !no repeat sides
        k=kxyface(j,i)    
        nck=cell2cell(k,i) !Forwards connectivity
        if(nck>ncells .or. iwet(nck)==0) cycle
        jcn=llec2llec(k,i) !Backwards connectivity        
        !Cell reconstruction corrections
        if(cell2upwdcell(k,i)==0)then  
          rpu=rpx(k,i)*dudx(i)+rpy(k,i)*dudy(i)
          rpv=rpx(k,i)*dvdx(i)+rpy(k,i)*dvdy(i)
          rnu=rpx(jcn,nck)*dudx(nck)+rpy(jcn,nck)*dudy(nck)
          rnv=rpx(jcn,nck)*dvdx(nck)+rpy(jcn,nck)*dvdy(nck)        
          !Skewness corrections
          ufnorm(k,i)=(acoef(k,i)+flux(k,i))*rpu-acoef(k,i)*rnu
          vfnorm(k,i)=(acoef(k,i)+flux(k,i))*rpv-acoef(k,i)*rnv
        endif
        !Upwinding
        if(flux(k,i)>0.0)then !Outflow
          if(cell2upwdcell(i,k)==0)then  
            dndun=dnx(k,i)*dudx(i)+dny(k,i)*dudy(i)
            dndvn=dnx(k,i)*dvdx(i)+dny(k,i)*dvdy(i)
            uc=u(i)+rpu; ud=u(nck)+rnu; uu=ud !Linear reconstruction
            vc=v(i)+rpv; vd=v(nck)+rnv; vu=vd !Linear reconstruction
          else
            nckk=cell2upwdcell(i,k)  
            dndun=0.0; dndvn=0.0  
            uc=u(i); ud=u(nck); uu=u(nckk)
            vc=v(i); vd=v(nck); vu=v(nckk)  
          endif
        else !Inflow
          if(cell2upwdcell(nck,jcn)==0)then  
            dndun=dnx(jcn,nck)*dudx(nck)+dny(jcn,nck)*dudy(nck)
            dndvn=dnx(jcn,nck)*dvdx(nck)+dny(jcn,nck)*dvdy(nck)
            ud=u(i)+rpu; uc=u(nck)+rnu; uu=ud !Linear reconstruction
            vd=v(i)+rpv; vc=v(nck)+rnv; vu=vd !Linear reconstruction
          else
            nckk=cell2upwdcell(nck,jcn)  
            dndun=0.0; dndvn=0.0  
            uc=u(nck); ud=u(i); uu=u(nckk)
            vc=v(nck); vd=v(i); vu=v(nckk)  
          endif
        endif
        if(abs(ud-uu)>1.0e-20)then !Avoid divide by zero
          uvnorm=(uc-uu)/(ud-uu) !Normalized
          if(uvnorm>0.0 .and. uvnorm<1.0)then
            ufnorm(k,i)=flux(k,i)*(ud-uc)*uvnorm  !Deferred normal flux
            ufnorm(jcn,nck)=-ufnorm(k,i)
          endif
        elseif(abs(dndun)>1.0e-10)then
          uvnorm=1.0-0.5*(ud-uc)/dndun  !Normalized
          if(uvnorm>0.0 .and. uvnorm<1.0)then
            ufnorm(k,i)=ufnorm(k,i)+flux(k,i)*(ud-uc)*uvnorm  !Deferred normal flux
          endif  
        endif
        if(abs(vd-vu)>1.0e-20)then !Avoid divide by zero
          uvnorm=(vc-vu)/(vd-vu) !Normalized
          if(uvnorm>0.0 .and. uvnorm<1.0)then
            vfnorm(k,i)=flux(k,i)*(vd-vc)*uvnorm  !Deferred normal flux
            vfnorm(jcn,nck)=-vfnorm(k,i)
          endif
        elseif(abs(dndvn)>1.0e-10)then
          uvnorm=1.0-0.5*(vd-vc)/dndvn  !Normalized
          if(uvnorm>0.0 .and. uvnorm<1.0)then
            vfnorm(k,i)=vfnorm(k,i)+flux(k,i)*(vd-vc)*uvnorm  !Deferred normal flux
          endif
        endif
        ufnorm(jcn,nck)=-ufnorm(k,i)
        vfnorm(jcn,nck)=-vfnorm(k,i)
      enddo
    enddo
!$OMP END DO

!---Sum fluxes to get deferred correction --------------
!$OMP DO PRIVATE(i)
    do i=1,ncells
      ssu(i)=ssu(i)-sum(ufnorm(1:ncface(i),i))
      ssv(i)=ssv(i)-sum(vfnorm(1:ncface(i),i))
    enddo
!$OMP END DO    
!$OMP END PARALLEL  

    return
    end subroutine defcorhlpagradvecnew    
    
!***********************************************************************
    subroutine defcorhlpagradvec(u,v,dudx,dudy,dvdx,dvdy,ssu,ssv)
! Deferred correction (anti-diffusion) for HLPA scheme.
!
! Uses input gradients for linear reconstructions and the calculation 
! of the normalized variable at joint cells.
!
! Version: 2.0
! Last updated: November 6, 2012
! Author: by Alex Sanchez, USACE-CHL
!***********************************************************************
    use size_def, only: ncells, ncellsd, nmaxfaces, ncellsimple, ncelljoint, ncellpoly
    use geo_def,  only: cell2cell, llec2llec, ncface, nxface, kxface, nyface, kyface, nxyface, kxyface, dnx, dny, rpx, rpy, idcelljoint, idcellsimple
    use flow_def, only: iwet,flux,acoef
    use prec_def, only: ikind
    
    implicit none
    !Input/Output
    real(ikind),intent(in),dimension(ncellsD) :: u,v,dudx,dudy,dvdx,dvdy
    real(ikind),intent(inout),dimension(ncellsD) :: ssu,ssv
    !Internal Variables
    integer :: ii,i,j,k,nck,jcn    
    real(ikind),dimension(nmaxfaces,ncellsD) :: ufnorm,vfnorm
    real(ikind) :: uvnorm,up,vp,un,vn,uc,vc,ud,vd,dndun,dndvn,rpu,rpv,rnu,rnv
    
!$OMP PARALLEL
!---- Initialize -----------------------------------------
!$OMP DO PRIVATE(i)
    do i=1,ncells
      ufnorm(1:ncface(i),i)=0.0
      vfnorm(1:ncface(i),i)=0.0
    enddo
!$OMP END DO
!--- Telescoping Grids --------------------------------------
!$OMP DO PRIVATE(i,ii,j,k,nck,jcn,uvnorm,ud,vd,uc,vc,dndun,dndvn)
    do ii=1,ncellsimple !Regular/simple cells
      i=idcellsimple(ii)
      if(iwet(i)==0) cycle
      !X-Direction
      do j=1,nxface(i) !No repeat cell face count
        k=kxface(j,i)   
        nck=cell2cell(k,i) !Forwards connectivity
        if(nck>ncells .or. iwet(nck)==0) cycle
        jcn=llec2llec(k,i) !Backwards connectivity
        !Upwinding and normal gradient
        if(flux(k,i)>0.0)then !Outflow
          dndun=dnx(k,i)*dudx(i)
          dndvn=dnx(k,i)*dvdx(i)
          uc=u(i); ud=u(nck)
          vc=v(i); vd=v(nck)
        else !Inflow
          dndun=dnx(jcn,nck)*dudx(nck)
          dndvn=dnx(jcn,nck)*dvdx(nck)
          uc=u(nck); ud=u(i)
          vc=v(nck); vd=v(i)
        endif
        if(abs(dndun)>1.0e-10)then
          uvnorm=1.0-0.5*(ud-uc)/dndun  !Normalized
          if(uvnorm>0.0 .and. uvnorm<1.0)then
            ufnorm(k,i)=flux(k,i)*(ud-uc)*uvnorm  !Deferred normal flux
            ufnorm(jcn,nck)=-ufnorm(k,i)
          endif 
        endif
        if(abs(dndvn)>1.0e-10)then
          uvnorm=1.0-0.5*(vd-vc)/dndvn  !Normalized
          if(uvnorm>0.0 .and. uvnorm<1.0)then
            vfnorm(k,i)=flux(k,i)*(vd-vc)*uvnorm  !Deferred normal flux
            vfnorm(jcn,nck)=-vfnorm(k,i)
          endif
        endif
      enddo
      !Y-Direction
      do j=1,nyface(i) !No repeat cell face count
        k=kyface(j,i)   
        nck=cell2cell(k,i) !Forwards connectivity
        if(nck>ncells .or. iwet(nck)==0) cycle
        jcn=llec2llec(k,i) !Backwards connectivity
        !Upwinding and normal gradient
        if(flux(k,i)>0.0)then !Outflow
          dndun=dny(k,i)*dudy(i)
          dndvn=dny(k,i)*dvdy(i)
          uc=u(i); ud=u(nck)
          vc=v(i); vd=v(nck)
        else !Inflow
          dndun=dny(jcn,nck)*dudy(nck)
          dndvn=dny(jcn,nck)*dvdy(nck)
          uc=u(nck); ud=u(i)
          vc=v(nck); vd=v(i)
        endif
        if(abs(dndun)>1.0e-10)then
          uvnorm=1.0-0.5*(ud-uc)/dndun  !Normalized
          if(uvnorm>0.0 .and. uvnorm<1.0)then
            ufnorm(k,i)=flux(k,i)*(ud-uc)*uvnorm  !Deferred normal flux
            ufnorm(jcn,nck)=-ufnorm(k,i)
          endif
        endif
        if(abs(dndvn)>1.0e-10)then
          uvnorm=1.0-0.5*(vd-vc)/dndvn  !Normalized
          if(uvnorm>0.0 .and. uvnorm<1.0)then
            vfnorm(k,i)=flux(k,i)*(vd-vc)*uvnorm  !Deferred normal flux
            vfnorm(jcn,nck)=-vfnorm(k,i)
          endif
        endif
      enddo
    enddo
!$OMP END DO
!$OMP DO PRIVATE(i,ii,j,k,nck,jcn,uvnorm,up,vp,un,vn,uc,vc,ud,vd,dndun,dndvn,rpu,rpv,rnu,rnv)
    do ii=1,ncelljoint !Joint Cells
      i=idcelljoint(ii)
      if(iwet(i)==0) cycle
      !X-Direction
      do j=1,nxface(i) !no repeat east and west sides
        k=kxface(j,i)    
        nck=cell2cell(k,i) !Forwards connectivity
        if(nck>ncells .or. iwet(nck)==0) cycle
        jcn=llec2llec(k,i) !Backwards connectivity
        up=u(i); vp=v(i); un=u(nck); vn=v(nck)
        if(abs(rpy(k,i))>1.0e-6)then !coarse to fine
          rpu=rpy(k,i)*dudy(i)
          rpv=rpy(k,i)*dvdy(i)
          up=up+rpu
          vp=vp+rpv
          ufnorm(k,i)=(acoef(k,i)+flux(k,i))*rpu !Deferred parallel flux
          vfnorm(k,i)=(acoef(k,i)+flux(k,i))*rpv !Deferred parallel flux
        elseif(abs(rpy(jcn,nck))>1.0e-6)then !fine to coarse
          rnu=rpy(jcn,nck)*dudy(nck)
          rnv=rpy(jcn,nck)*dvdy(nck)
          un=un+rnu
          vn=vn+rnv
          ufnorm(k,i)=-acoef(k,i)*rnu !Deferred parallel flux
          vfnorm(k,i)=-acoef(k,i)*rnv !Deferred parallel flux
        endif
        !Upwinding and normal gradient
        if(flux(k,i)>0.0)then !Outflow
          dndun=dnx(k,i)*dudx(i)
          dndvn=dnx(k,i)*dvdx(i)
          uc=up; ud=un
          vc=vp; vd=vn
        else !Inflow
          dndun=dnx(jcn,nck)*dudx(nck)
          dndvn=dnx(jcn,nck)*dvdx(nck)
          ud=up; uc=un
          vd=vp; vc=vn
        endif
        if(abs(dndun)>1.0e-10)then
          uvnorm=1.0-0.5*(ud-uc)/dndun  !Normalized
          if(uvnorm>0.0 .and. uvnorm<1.0)then
            ufnorm(k,i)=ufnorm(k,i)+flux(k,i)*(ud-uc)*uvnorm  !Deferred normal flux
          endif
        endif
        if(abs(dndvn)>1.0e-10)then
          uvnorm=1.0-0.5*(vd-vc)/dndvn  !Normalized
          if(uvnorm>0.0 .and. uvnorm<1.0)then
            vfnorm(k,i)=vfnorm(k,i)+flux(k,i)*(vd-vc)*uvnorm  !Deferred normal flux
          endif
        endif
        ufnorm(jcn,nck)=-ufnorm(k,i); vfnorm(jcn,nck)=-vfnorm(k,i)
      enddo
      !Y-Direction
      do j=1,nyface(i) !no repeat north and south sides
        k=kyface(j,i)
        nck=cell2cell(k,i) !Forwards connectivity
        if(nck>ncells .or. iwet(nck)==0) cycle
        jcn=llec2llec(k,i) !Backwards connectivity
        up=u(i); vp=v(i); un=u(nck); vn=v(nck)
        if(abs(rpx(k,i))>1.0e-6)then !coarse to fine
          rpu=rpx(k,i)*dudx(i)
          rpv=rpx(k,i)*dvdx(i)
          up=up+rpu
          vp=vp+rpv
          ufnorm(k,i)=(acoef(k,i)+flux(k,i))*rpu !Deferred parallel flux
          vfnorm(k,i)=(acoef(k,i)+flux(k,i))*rpv !Deferred parallel flux
        elseif(abs(rpx(jcn,nck))>1.0e-6)then !fine to coarse
          rnu=rpx(jcn,nck)*dudx(nck)
          rnv=rpx(jcn,nck)*dvdx(nck)
          un=un+rnu
          vn=vn+rnv
          ufnorm(k,i)=-acoef(k,i)*rnu !Deferred parallel flux
          vfnorm(k,i)=-acoef(k,i)*rnv !Deferred parallel flux
        endif
        !Upwinding and normal gradient
        if(flux(k,i)>0.0)then !Outflow
          dndun=dny(k,i)*dudy(i)
          dndvn=dny(k,i)*dvdy(i)
          uc=up; ud=un
          vc=vp; vd=vn
        else !Inflow
          dndun=dny(jcn,nck)*dudy(nck)
          dndvn=dny(jcn,nck)*dvdy(nck)
          ud=up; uc=un
          vd=vp; vc=vn 
        endif
        if(abs(dndun)>1.0e-10)then
          uvnorm=1.0-0.5*(ud-uc)/dndun  !Normalized
          if(uvnorm>0.0 .and. uvnorm<1.0)then
            ufnorm(k,i)=ufnorm(k,i)+flux(k,i)*(ud-uc)*uvnorm  !Deferred normal flux
          endif
        endif
        if(abs(dndvn)>1.0e-10)then
          uvnorm=1.0-0.5*(vd-vc)/dndvn  !Normalized
          if(uvnorm>0.0 .and. uvnorm<1.0)then
            vfnorm(k,i)=vfnorm(k,i)+flux(k,i)*(vd-vc)*uvnorm  !Deferred normal flux
          endif
        endif
        ufnorm(jcn,nck)=-ufnorm(k,i)
        vfnorm(jcn,nck)=-vfnorm(k,i)
      enddo
    enddo
!$OMP END DO

!--- Polyhedral Grids --------------------------------------
!$OMP DO PRIVATE(i,j,k,nck,jcn,uvnorm,uc,vc,ud,vd,dndun,dndvn,rnu,rpu,rpv,rnv)
    do i=1,ncellpoly !No index mapping required for polygonal cells (all cells are polygonal)
      if(iwet(i)==0) cycle
      do j=1,nxyface(i) !no repeat east and west sides
        k=kxyface(j,i)    
        nck=cell2cell(k,i) !Forwards connectivity
        if(nck>ncells .or. iwet(nck)==0) cycle
        jcn=llec2llec(k,i) !Backwards connectivity        
        !Cell reconstruction corrections
        rpu=rpx(k,i)*dudx(i)+rpy(k,i)*dudy(i)
        rpv=rpx(k,i)*dvdx(i)+rpy(k,i)*dvdy(i)
        rnu=rpx(jcn,nck)*dudx(nck)+rpy(jcn,nck)*dudy(nck)
        rnv=rpx(jcn,nck)*dvdx(nck)+rpy(jcn,nck)*dvdy(nck)
        !Skewness corrections
        ufnorm(k,i)=(acoef(k,i)+flux(k,i))*rpu-acoef(k,i)*rnu
        vfnorm(k,i)=(acoef(k,i)+flux(k,i))*rpv-acoef(k,i)*rnv
        !Upwinding
        if(flux(k,i)>0.0)then !Outflow
          dndun=dnx(k,i)*dudx(i)+dny(k,i)*dudy(i)
          dndvn=dnx(k,i)*dvdx(i)+dny(k,i)*dvdy(i)
          uc=u(i)+rpu; ud=u(nck)+rnu !Linear reconstruction
          vc=v(i)+rpv; vd=v(nck)+rnv !Linear reconstruction
        else !Inflow
          dndun=dnx(jcn,nck)*dudx(nck)+dny(jcn,nck)*dudy(nck)
          dndvn=dnx(jcn,nck)*dvdx(nck)+dny(jcn,nck)*dvdy(nck)
          ud=u(i)+rpu; uc=u(nck)+rnu !Linear reconstruction
          vd=v(i)+rpv; vc=v(nck)+rnv !Linear reconstruction
        endif
        if(abs(dndun)>1.0e-10)then
          uvnorm=1.0-0.5*(ud-uc)/dndun  !Normalized
          if(uvnorm>0.0 .and. uvnorm<1.0)then
            ufnorm(k,i)=ufnorm(k,i)+flux(k,i)*(ud-uc)*uvnorm  !Deferred normal flux
          endif
        endif
        if(abs(dndvn)>1.0e-10)then
          uvnorm=1.0-0.5*(vd-vc)/dndvn  !Normalized
          if(uvnorm>0.0 .and. uvnorm<1.0)then
            vfnorm(k,i)=vfnorm(k,i)+flux(k,i)*(vd-vc)*uvnorm  !Deferred normal flux
          endif
        endif
        ufnorm(jcn,nck)=-ufnorm(k,i)
        vfnorm(jcn,nck)=-vfnorm(k,i)
      enddo
    enddo
!$OMP END DO

!---Sum fluxes to get deferred correction --------------
!$OMP DO PRIVATE(i)
    do i=1,ncells
      ssu(i)=ssu(i)-sum(ufnorm(1:ncface(i),i))
      ssv(i)=ssv(i)-sum(vfnorm(1:ncface(i),i))
    enddo
!$OMP END DO    
!$OMP END PARALLEL  

    return
    end subroutine defcorhlpagradvec

!****************************************************************************
    subroutine defcorgammagradvec(schmdefcor,u,v,dudx,dudy,dvdx,dvdy,ssu,ssv)
! Deferred correction (anti-diffusion) for Gamma-family schemes
!
! Version: 3.0
! Last updated: January 15, 2014
! Author: by Alex Sanchez, USACE-CHL
!****************************************************************************
    use size_def, only: ncells, ncellsd, nmaxfaces, ncellsimple, ncelljoint, ncellpoly
    use geo_def,  only: ncface, nxface, kxface, nyface, kyface, nxyface, kxyface, dnx, dny, rpx, rpy, idcellsimple, idcelljoint, cell2cell, llec2llec
    use prec_def, only: ikind
    use comvarbl, only: skewcor
    use flow_def, only: iwet,flux,acoef
    use interp_def, only: fintp

    implicit none
    !Input/Output
    real(ikind),intent(in),dimension(ncellsD) :: u,v,dudx,dudy,dvdx,dvdy
    real(ikind),intent(inout),dimension(ncellsD) :: ssu,ssv
    !Internal Variables
    integer :: ii,i,j,k,nck,jcn    
    real(ikind),dimension(nmaxfaces,ncellsD) :: ufnorm,vfnorm
    real(ikind) :: unorm,up,un,uc,ud,uf,dndun,rpu,rnu
    real(ikind) :: vnorm,vp,vn,vc,vd,vf,dndvn,rpv,rnv

    interface
      function schmdefcor(phin)
        use prec_def
        implicit none
        real(ikind),intent(in) :: phin
        real(ikind) :: schmdefcor
      end function
    endinterface

!--- All other types of grids ------------------------
!$OMP PARALLEL
!$OMP DO PRIVATE(i)
    do i=1,ncells
      ufnorm(1:ncface(i),i)=0.0
      vfnorm(1:ncface(i),i)=0.0
    enddo
!$OMP END DO
!--- Telescoping Grids --------------------------------------
!$OMP DO PRIVATE(i,ii,j,k,nck,jcn,unorm,up,un,uc,ud,uf,dndun,rpu,rnu,vnorm,vp,vn,vc,vd,vf,dndvn,rpv,rnv)
    do ii=1,ncellsimple
      i=idcellsimple(ii)
      if(iwet(i)==0) cycle
      !X-Direction
      do j=1,nxface(i) !no repeat east and west sides
        k=kxface(j,i)    
        nck=cell2cell(k,i) !Forwards connectivity
        if(nck>ncells .or. iwet(nck)==0) cycle
        jcn=llec2llec(k,i) !Backwards connectivity
        !Upwinding and normal gradient
        if(flux(k,i)>0.0)then !Outflow
          dndun=dnx(k,i)*dudx(i)
          dndvn=dnx(k,i)*dvdx(i)
          uc=u(i); ud=u(nck)
          vc=v(i); vd=v(nck)
        else !Inflow
          dndun=dnx(jcn,nck)*dudx(nck)
          dndvn=dnx(jcn,nck)*dvdx(nck)
          uc=u(nck); ud=u(i)
          vc=v(nck); vd=v(i)
        endif
        if(abs(dndun)>1.0e-10)then
          unorm=1.0-0.5*(ud-uc)/dndun  !Normalized
          if(unorm>0.0 .and. unorm<1.0)then
            uf=fintp(k,i)*u(nck)+(1.0-fintp(k,i))*u(i) !Central Difference
            ufnorm(k,i)=flux(k,i)*schmdefcor(unorm)*(uf-uc)   !Deferred normal flux
            ufnorm(jcn,nck)=-ufnorm(k,i)
          endif
        endif
        if(abs(dndvn)>1.0e-10)then
          vnorm=1.0-0.5*(vd-vc)/dndvn  !Normalized
          if(vnorm>0.0 .and. vnorm<1.0)then
            vf=fintp(k,i)*v(nck)+(1.0-fintp(k,i))*v(i) !Central Difference
            vfnorm(k,i)=flux(k,i)*schmdefcor(vnorm)*(vf-vc)   !Deferred normal flux
            vfnorm(jcn,nck)=-vfnorm(k,i)
          endif
        endif
      enddo
      !Y-Direction
      do j=1,nyface(i) !no repeat north and south faces
        k=kyface(j,i)
        nck=cell2cell(k,i) !Forwards connectivity
        if(nck>ncells .or. iwet(nck)==0) cycle
        jcn=llec2llec(k,i) !Backwards connectivity
        !Upwinding and normal gradient
        if(flux(k,i)>0.0)then !Outflow
          dndun=dny(k,i)*dudy(i)
          dndvn=dny(k,i)*dvdy(i)
          uc=u(i); ud=u(nck)
          vc=v(i); vd=v(nck)
        else !Inflow
          dndun=dny(jcn,nck)*dudy(nck)
          dndvn=dny(jcn,nck)*dvdy(nck)
          uc=u(nck); ud=u(i)
          vc=v(nck); vd=v(i)
        endif
        if(abs(dndun)>1.0e-10)then
          unorm=1.0-0.5*(ud-uc)/dndun  !Normalized
          if(unorm>0.0 .and. unorm<1.0)then
            uf=fintp(k,i)*u(nck)+(1.0-fintp(k,i))*u(i) !Central Difference
            ufnorm(k,i)=flux(k,i)*schmdefcor(unorm)*(uf-uc)   !Deferred normal flux
            ufnorm(jcn,nck)=-ufnorm(k,i)
          endif
        endif
        if(abs(dndvn)>1.0e-10)then
          vnorm=1.0-0.5*(vd-vc)/dndvn  !Normalized
          if(vnorm>0.0 .and. vnorm<1.0)then
            vf=fintp(k,i)*v(nck)+(1.0-fintp(k,i))*v(i) !Central Difference
            vfnorm(k,i)=flux(k,i)*schmdefcor(vnorm)*(vf-vc)   !Deferred normal flux
            vfnorm(jcn,nck)=-vfnorm(k,i)
          endif
        endif
      enddo
    enddo
!$OMP END DO
!$OMP DO PRIVATE(i,ii,j,k,nck,jcn,unorm,up,un,uc,ud,uf,dndun,rpu,rnu,vnorm,vp,vn,vc,vd,vf,dndvn,rpv,rnv)
    do ii=1,ncelljoint !Joint cells
      i=idcelljoint(ii)
      if(iwet(i)==0) cycle
      !X-Direction
      do j=1,nxface(i) !no repeat east and west sides
        k=kxface(j,i)    
        nck=cell2cell(k,i) !Forwards connectivity
        if(nck>ncells .or. iwet(nck)==0) cycle
        jcn=llec2llec(k,i) !Backwards connectivity
        !Linear cell reconstructions and sckewness correction
        up=u(i); un=u(nck); vp=v(i); vn=v(nck)
        if(abs(rpy(k,i))>1.0e-6)then !coarse to fine
          rpu=rpy(k,i)*dudy(i)
          rpv=rpy(k,i)*dvdy(i) 
          up=up+rpu
          vp=vp+rpv
          ufnorm(k,i)=(acoef(k,i)+flux(k,i))*rpu !Deferred parallel flux
          vfnorm(k,i)=(acoef(k,i)+flux(k,i))*rpv !Deferred parallel flux
        elseif(abs(rpy(jcn,nck))>1.0e-6)then !fine to coarse
          rnu=rpy(jcn,nck)*dudy(nck)
          rnv=rpy(jcn,nck)*dvdy(nck) 
          un=un+rnu
          vn=vn+rnv
          ufnorm(k,i)=-acoef(k,i)*rnu !Deferred parallel flux
          vfnorm(k,i)=-acoef(k,i)*rnv !Deferred parallel flux
        endif
        !Upwinding and normal distance
        if(flux(k,i)>0.0)then !Outflow
          dndun=dnx(k,i)*dudx(i)
          dndvn=dnx(k,i)*dvdx(i)
          uc=up; ud=un
          vc=vp; vd=vn
        else !Inflow
          dndun=dnx(jcn,nck)*dudx(nck)
          dndvn=dnx(jcn,nck)*dvdx(nck)
          ud=up; uc=un
          vd=vp; vc=vn
        endif
        if(abs(dndun)>1.0e-10)then
          unorm=1.0-0.5*(ud-uc)/dndun  !Normalized
          if(unorm>0.0 .and. unorm<1.0)then
            uf=fintp(k,i)*un+(1.0-fintp(k,i))*up !Central Difference
            ufnorm(k,i)=ufnorm(k,i)+flux(k,i)*schmdefcor(unorm)*(uf-uc) !Deferred normal flux
          endif
        endif
        if(abs(dndvn)>1.0e-10)then
          vnorm=1.0-0.5*(vd-vc)/dndvn  !Normalized
          if(vnorm>0.0 .and. vnorm<1.0)then
            vf=fintp(k,i)*vn+(1.0-fintp(k,i))*vp !Central Difference
            vfnorm(k,i)=vfnorm(k,i)+flux(k,i)*schmdefcor(vnorm)*(vf-vc) !Deferred normal flux
          endif
        endif
        ufnorm(jcn,nck)=-ufnorm(k,i)
        vfnorm(jcn,nck)=-vfnorm(k,i)
      enddo
      !Y-Direction
      do j=1,nyface(i) !no repeat north and south faces
        k=kyface(j,i)
        nck=cell2cell(k,i) !Forwards connectivity
        if(nck>ncells .or. iwet(nck)==0) cycle
        jcn=llec2llec(k,i) !Backwards connectivity
        !Linear cell reconstructions and sckewness correction
        up=u(i); un=u(nck)
        vp=v(i); vn=v(nck)
        if(abs(rpx(k,i))>1.0e-6)then !coarse to fine
          rpu=rpx(k,i)*dudx(i)
          rpv=rpx(k,i)*dvdx(i) 
          up=up+rpu
          vp=vp+rpv
          ufnorm(k,i)=(acoef(k,i)+flux(k,i))*rpu !Deferred parallel flux
          vfnorm(k,i)=(acoef(k,i)+flux(k,i))*rpv !Deferred parallel flux
        elseif(abs(rpx(jcn,nck))>1.0e-6)then !fine to coarse
          rnu=rpx(jcn,nck)*dudx(nck)
          rnv=rpx(jcn,nck)*dvdx(nck) 
          un=un+rnu
          vn=vn+rnv
          ufnorm(k,i)=-acoef(k,i)*rnu !Deferred parallel flux
          vfnorm(k,i)=-acoef(k,i)*rnv !Deferred parallel flux
        endif
        !Upwinding and normal distances
        if(flux(k,i)>0.0)then !Outflow
          dndun=dny(k,i)*dudy(i)
          dndvn=dny(k,i)*dvdy(i)
          uc=up; ud=un
          vc=vp; vd=vn
        else !Inflow
          dndun=dny(jcn,nck)*dudy(nck)
          dndvn=dny(jcn,nck)*dvdy(nck)
          ud=up; uc=un
          vd=vp; vc=vn
        endif
        if(abs(dndun)>1.0e-10)then
          unorm=1.0-0.5*(ud-uc)/dndun  !Normalized
          if(unorm>0.0 .and. unorm<1.0)then
            uf=fintp(k,i)*un+(1.0-fintp(k,i))*up !Central Difference
            ufnorm(k,i)=ufnorm(k,i)+flux(k,i)*schmdefcor(unorm)*(uf-uc) !Deferred normal flux
          endif
        endif
        if(abs(dndvn)>1.0e-10)then
          vnorm=1.0-0.5*(vd-vc)/dndvn  !Normalized
          if(vnorm>0.0 .and. vnorm<1.0)then
            vf=fintp(k,i)*vn+(1.0-fintp(k,i))*vp !Central Difference
            vfnorm(k,i)=vfnorm(k,i)+flux(k,i)*schmdefcor(vnorm)*(vf-vc) !Deferred normal flux
          endif
        endif
        ufnorm(jcn,nck)=-ufnorm(k,i)
        vfnorm(jcn,nck)=-vfnorm(k,i)
      enddo
    enddo
!$OMP END DO

    !Polygonal Meshes
!$OMP DO PRIVATE(i,j,k,nck,jcn,unorm,uc,ud,uf,dndun,rpu,rnu,vnorm,vc,vd,vf,dndvn,rpv,rnv)
    do i=1,ncellpoly
      if(iwet(i)==0) cycle
      do j=1,nxyface(i) !no repeat east and west sides
        k=kxyface(j,i)    
        nck=cell2cell(k,i) !Forwards connectivity
        if(nck>ncells .or. iwet(nck)==0) cycle
        jcn=llec2llec(k,i) !Backwards connectivity
        !Cell reconstruction corrections
        rpu=rpx(k,i)*dudx(i)+rpy(k,i)*dudy(i)
        rpv=rpx(k,i)*dvdx(i)+rpy(k,i)*dvdy(i)
        rnu=rpx(jcn,nck)*dudx(nck)+rpy(jcn,nck)*dudy(nck) 
        rnv=rpx(jcn,nck)*dvdx(nck)+rpy(jcn,nck)*dvdy(nck) 
        up=u(i)+rpu; un=u(i)+rnu
        vp=v(i)+rpv; vn=v(i)+rnv
        !Skewness corrections
        ufnorm(k,i)=(acoef(k,i)+flux(k,i))*rpu-acoef(k,i)*rnu
        vfnorm(k,i)=(acoef(k,i)+flux(k,i))*rpv-acoef(k,i)*rnv
        !Upwinding and linear reconstructions
        if(flux(k,i)>0.0)then !Outflow
          dndun=dnx(k,i)*dudx(i)+dny(k,i)*dudy(i)
          dndvn=dnx(k,i)*dvdx(i)+dny(k,i)*dvdy(i)
          uc=up; ud=un 
          vc=vp; vd=vn 
        else !Inflow
          dndun=dnx(jcn,nck)*dudx(nck)+dny(jcn,nck)*dudy(nck)
          dndvn=dnx(jcn,nck)*dvdx(nck)+dny(jcn,nck)*dvdy(nck)
          ud=up; uc=un
          vd=vp; vc=vn
        endif
        if(abs(dndun)>1.0e-10)then
          unorm=1.0-0.5*(ud-uc)/dndun  !Normalized
          if(unorm>0.0 .and. unorm<1.0)then
            uf=fintp(k,i)*un+(1.0-fintp(k,i))*up !Central Difference
            ufnorm(k,i)=ufnorm(k,i)+flux(k,i)*schmdefcor(unorm)*(uf-uc) !Deferred normal flux  
          endif
        endif
        if(abs(dndvn)>1.0e-10)then
          vnorm=1.0-0.5*(vd-vc)/dndvn  !Normalized
          if(vnorm>0.0 .and. vnorm<1.0)then
            vf=fintp(k,i)*vn+(1.0-fintp(k,i))*vp !Central Difference
            vfnorm(k,i)=vfnorm(k,i)+flux(k,i)*schmdefcor(vnorm)*(vf-vc) !Deferred normal flux  
          endif
        endif
        ufnorm(jcn,nck)=-ufnorm(k,i)
        vfnorm(jcn,nck)=-vfnorm(k,i)
      enddo  
    enddo
!$OMP END DO
!$OMP DO PRIVATE(i)
    do i=1,ncells
      ssu(i)=ssu(i)-sum(ufnorm(1:ncface(i),i))
      ssv(i)=ssv(i)-sum(vfnorm(1:ncface(i),i))
    enddo
!$OMP END DO    
!$OMP END PARALLEL
    
    return
    end subroutine defcorgammagradvec
    
!******************************************************************************
    subroutine defcorgammagradvecold(schmdefcor,u,v,dudx,dudy,dvdx,dvdy,ssu,ssv)
! Deferred correction (anti-diffusion) for Gamma-family schemes
! for vectors
!
! Version: 2.0
! Last updated: November 6, 2012
! Author: by Alex Sanchez, USACE-CHL
!******************************************************************************
    use size_def, only: ncells, ncellsd, nmaxfaces, ncellsimple, ncelljoint, ncellpoly
    use geo_def,  only: cell2cell, llec2llec, ncface, nxface, kxface, nyface, kyface, nxyface, kxyface, dnx, dny, rpx, rpy, idcelljoint, idcellsimple
    use flow_def, only: iwet,flux,acoef
    use prec_def, only: ikind
    use interp_def, only: fintp
    
    implicit none
    !Input/Output
    real(ikind),intent(in),dimension(ncellsD) :: u,v,dudx,dudy,dvdx,dvdy
    real(ikind),intent(inout),dimension(ncellsD) :: ssu,ssv
    !Internal Variables
    integer :: ii,i,j,k,nck,jcn    
    real(ikind) :: uvnorm,up,un,uc,ud,uvf,dndun,vp,vn,vc,vd,dndvn,rpu,rpv,rnu,rnv
    real(ikind),dimension(nmaxfaces,ncellsD) :: ufnorm,vfnorm

    interface
      function schmdefcor(phin)
        use prec_def
        implicit none
        real(ikind),intent(in) :: phin
        real(ikind) :: schmdefcor
      end function
    endinterface
 
!--- All other types of grids ------------------------
!$OMP PARALLEL
!$OMP DO PRIVATE(i)
    do i=1,ncells
      ufnorm(1:ncface(i),i)=0.0
      vfnorm(1:ncface(i),i)=0.0
    enddo
!$OMP END DO
!$OMP DO PRIVATE(i,ii,j,k,nck,jcn,uvnorm,uc,vc,ud,vd,uvf,dndun,dndvn)
    do ii=1,ncellsimple !Regular/simple cells
      i=idcellsimple(ii)
      if(iwet(i)==0) cycle
      !X-Direction
      do j=1,nxface(i) !no repeat east and west sides
        k=kxface(j,i)    
        nck=cell2cell(k,i) !Forwards connectivity
        if(nck>ncells .or. iwet(nck)==0) cycle
        jcn=llec2llec(k,i) !Backwards connectivity  
        !Upwinding and normal gradients
        if(flux(k,i)>0.0)then !Outflow
          dndun=dnx(k,i)*dudx(i); dndvn=dnx(k,i)*dvdx(i)
          uc=v(i); ud=u(nck); vc=v(i); vd=v(nck)
        else !Inflow
          dndun=dnx(jcn,nck)*dudx(nck); dndvn=dnx(jcn,nck)*dvdx(nck)
          uc=u(nck); ud=u(i); vc=v(nck); vd=v(i)
        endif
        if(abs(dndun)>1.0e-10)then
          uvnorm=1.0-0.5*(ud-uc)/dndun  !Normalized
          if(uvnorm>0.0 .and. uvnorm<1.0)then
            uvf=fintp(k,i)*u(nck)+(1.0-fintp(k,i))*u(i) !Central Difference
            ufnorm(k,i)=flux(k,i)*schmdefcor(uvnorm)*(uvf-uc)   !Deferred normal flux
            ufnorm(jcn,nck)=-ufnorm(k,i)
          endif
        endif
        if(abs(dndvn)>1.0e-10)then
          uvnorm=1.0-0.5*(vd-vc)/dndvn  !Normalized
          if(uvnorm>0.0 .and. uvnorm<1.0)then
            uvf=fintp(k,i)*v(nck)+(1.0-fintp(k,i))*v(i) !Central Difference
            vfnorm(k,i)=flux(k,i)*schmdefcor(uvnorm)*(uvf-vc)   !Deferred normal flux
            vfnorm(jcn,nck)=-vfnorm(k,i)
          endif
        endif
      enddo
      !Y-Direction
      do j=1,nyface(i) !no repeat north and south faces
        k=kyface(j,i)
        nck=cell2cell(k,i) !Forwards connectivity
        if(nck>ncells .or. iwet(nck)==0) cycle
        jcn=llec2llec(k,i) !Backwards connectivity  
        !Upwinding and normal gradient
        if(flux(k,i)>0.0)then !Outflow
          dndun=dny(k,i)*dudy(i); dndvn=dny(k,i)*dvdy(i)
          uc=u(i); ud=u(nck); vc=v(i); vd=v(nck)
        else !Inflow
          dndun=dny(jcn,nck)*dudy(nck); dndvn=dny(jcn,nck)*dvdy(nck)
          uc=u(nck); ud=u(i); vc=v(nck); vd=v(i)
        endif
        if(abs(dndun)>1.0e-10)then
          uvnorm=1.0-0.5*(ud-uc)/dndun  !Normalized
          if(uvnorm>0.0 .and. uvnorm<1.0)then
            uvf=fintp(k,i)*u(nck)+(1.0-fintp(k,i))*u(i) !Central Difference
            ufnorm(k,i)=flux(k,i)*schmdefcor(uvnorm)*(uvf-uc)   !Deferred normal flux  
            ufnorm(jcn,nck)=-ufnorm(k,i)
          endif
        endif
        if(abs(dndvn)>1.0e-10)then
          uvnorm=1.0-0.5*(vd-vc)/dndvn  !Normalized
          if(uvnorm>0.0 .and. uvnorm<1.0)then
            uvf=fintp(k,i)*v(nck)+(1.0-fintp(k,i))*v(i) !Central Difference
            vfnorm(k,i)=flux(k,i)*schmdefcor(uvnorm)*(uvf-vc)   !Deferred normal flux  
            vfnorm(jcn,nck)=-vfnorm(k,i)
          endif
        endif
      enddo
    enddo
!$OMP END DO
!$OMP DO PRIVATE(i,ii,j,k,nck,jcn,uvnorm,up,vp,un,vn,uc,vc,vd,ud,uvf,dndun,dndvn,rpu,rpv,rnu,rnv)
    do ii=1,ncelljoint !Joint cells
      i=idcelljoint(ii)
      if(iwet(i)==0) cycle
      !X-Direction
      do j=1,nxface(i) !no repeat east and west sides
        k=kxface(j,i)    
        nck=cell2cell(k,i) !Forwards connectivity
        if(nck>ncells .or. iwet(nck)==0) cycle
        jcn=llec2llec(k,i) !Backwards connectivity  
        !Cell reconstruction and skewness corrections
        up=u(i); vp=v(i); un=u(nck); vn=v(nck)
        if(abs(rpy(k,i))>1.0e-6)then !coarse to fine
          rpu=rpy(k,i)*dudy(i); rpv=rpy(k,i)*dvdy(i)
          up=up+rpu; vp=vp+rpv
          ufnorm(k,i)=(acoef(k,i)+flux(k,i))*rpu !Deferred parallel flux
          vfnorm(k,i)=(acoef(k,i)+flux(k,i))*rpv !Deferred parallel flux
        elseif(abs(rpy(jcn,nck))>1.0e-6)then !fine to coarse
          rnu=rpy(jcn,nck)*dudy(nck); rnv=rpy(jcn,nck)*dvdy(nck)
          un=un+rnu; vn=vn+rnv
          ufnorm(k,i)=-acoef(k,i)*rnu !Deferred parallel flux
          vfnorm(k,i)=-acoef(k,i)*rnv !Deferred parallel flux          
        endif
        !Upwinding and normal gradient
        if(flux(k,i)>0.0)then !Outflow
          dndun=dnx(k,i)*dudx(i); dndvn=dnx(k,i)*dvdx(i)
          uc=up; ud=un; vc=vp; vd=vn
        else !Inflow
          dndun=dnx(jcn,nck)*dudx(nck); dndvn=dnx(jcn,nck)*dvdx(nck)
          ud=up; uc=un; vd=vp; vc=vn 
        endif
        if(abs(dndun)>1.0e-10)then
          uvnorm=1.0-0.5*(ud-uc)/dndun  !Normalized
          if(uvnorm>0.0 .and. uvnorm<1.0)then
            uvf=fintp(k,i)*un+(1.0-fintp(k,i))*up !Central Difference
            ufnorm(k,i)=ufnorm(k,i)+flux(k,i)*schmdefcor(uvnorm)*(uvf-uc)   !Deferred normal flux
          endif
        endif
        if(abs(dndvn)>1.0e-10)then
          uvnorm=1.0-0.5*(vd-vc)/dndvn  !Normalized
          if(uvnorm>0.0 .and. uvnorm<1.0)then
            uvf=fintp(k,i)*vn+(1.0-fintp(k,i))*vp !Central Difference
            vfnorm(k,i)=vfnorm(k,i)+flux(k,i)*schmdefcor(uvnorm)*(uvf-vc)   !Deferred normal flux            
          endif
        endif
        ufnorm(jcn,nck)=-ufnorm(k,i)
        vfnorm(jcn,nck)=-vfnorm(k,i)
      enddo
      !Y-Direction
      do j=1,nyface(i) !no repeat north and south faces
        k=kyface(j,i)
        nck=cell2cell(k,i) !Forwards connectivity
        if(nck>ncells .or. iwet(nck)==0) cycle
        jcn=llec2llec(k,i) !Backwards connectivity 
        !Cell reconstruction corrections and skewness corrections
        up=u(i); vp=v(i); un=u(nck); vn=v(nck)
        if(abs(rpx(k,i))>1.0e-6)then !coarse to fine
          rpu=rpx(k,i)*dudx(i); rpv=rpx(k,i)*dvdx(i)
          up=up+rpu; vp=vp+rpv
          ufnorm(k,i)=(acoef(k,i)+flux(k,i))*rpu !Deferred parallel flux
          vfnorm(k,i)=(acoef(k,i)+flux(k,i))*rpv !Deferred parallel flux
        elseif(abs(rpx(jcn,nck))>1.0e-6)then !fine to coarse
          rnu=rpx(jcn,nck)*dudx(nck); rnv=rpx(jcn,nck)*dvdx(nck)
          un=un+rnu; vn=vn+rnv
          ufnorm(k,i)=-acoef(k,i)*rnu !Deferred parallel flux
          vfnorm(k,i)=-acoef(k,i)*rnv !Deferred parallel flux
        endif
        !Upwinding
        if(flux(k,i)>0.0)then !Outflow
          dndun=dny(k,i)*dudy(i); dndvn=dny(k,i)*dvdy(i)
          uc=up; ud=un; vc=vp; vd=vn
        else !Inflow
          dndun=dny(jcn,nck)*dudy(nck); dndvn=dny(jcn,nck)*dvdy(nck)
          ud=up; uc=un; vd=vp; vc=vn 
        endif
        if(abs(dndun)>1.0e-10)then
          uvnorm=1.0-0.5*(ud-uc)/dndun  !Normalized
          if(uvnorm>0.0 .and. uvnorm<1.0)then
            uvf=fintp(k,i)*un+(1.0-fintp(k,i))*up !Central Difference
            ufnorm(k,i)=ufnorm(k,i)+flux(k,i)*schmdefcor(uvnorm)*(uvf-uc)   !Deferred normal flux
          endif
        endif
        if(abs(dndvn)>1.0e-10)then
          uvnorm=1.0-0.5*(vd-vc)/dndvn  !Normalized
          if(uvnorm>0.0 .and. uvnorm<1.0)then
            uvf=fintp(k,i)*vn+(1.0-fintp(k,i))*vp !Central Difference
            vfnorm(k,i)=vfnorm(k,i)+flux(k,i)*schmdefcor(uvnorm)*(uvf-vc)   !Deferred normal flux          
          endif
        endif
        ufnorm(jcn,nck)=-ufnorm(k,i)
        vfnorm(jcn,nck)=-vfnorm(k,i)
      enddo
    enddo
!$OMP END DO

    !Polygonal Meshes
!$OMP DO PRIVATE(i,j,k,nck,jcn,uvnorm,uc,vc,ud,vd,up,vp,un,vn,uvf,dndun,dndvn,rpu,rpv,rnu,rnv)
    do i=1,ncellpoly
      if(iwet(i)==0) cycle
      do j=1,nxyface(i) !no repeat east and west sides
        k=kxyface(j,i)    
        nck=cell2cell(k,i) !Forwards connectivity
        if(nck>ncells .or. iwet(nck)==0) cycle
        jcn=llec2llec(k,i) !Backwards connectivity 
        !Cell reconstruction corrections
        rpu=rpx(k,i)*dudx(i)+rpy(k,i)*dudy(i)
        rnu=rpx(jcn,nck)*dudx(nck)+rpy(jcn,nck)*dudy(nck)
        rpv=rpx(k,i)*dvdx(i)+rpy(k,i)*dvdy(i)
        rnv=rpx(jcn,nck)*dvdx(nck)+rpy(jcn,nck)*dvdy(nck)
        up=u(i)+rpu; un=u(i)+rnu
        vp=v(i)+rpv; vn=v(i)+rnv
        !Skewness correction
        ufnorm(k,i)=(acoef(k,i)+flux(k,i))*rpu-acoef(k,i)*rnu
        vfnorm(k,i)=(acoef(k,i)+flux(k,i))*rpv-acoef(k,i)*rnv
        !Upwinding, linear reconstructions, and normal distances
        if(flux(k,i)>0.0)then !Outflow
          dndun=dnx(k,i)*dudx(i)+dny(k,i)*dudy(i)
          dndvn=dnx(k,i)*dvdx(i)+dny(k,i)*dvdy(i)
          uc=up; vc=vp; ud=un; vd=vn
        else !Inflow
          dndun=dnx(jcn,nck)*dudx(nck)+dny(jcn,nck)*dudy(nck)
          dndvn=dnx(jcn,nck)*dvdx(nck)+dny(jcn,nck)*dvdy(nck)
          ud=up; vd=vp; uc=un; vc=vn
        endif
        if(abs(dndun)>1.0e-10)then
          uvnorm=1.0-0.5*(ud-uc)/dndun  !Normalized
          if(uvnorm>0.0 .and. uvnorm<1.0)then
            uvf=fintp(k,i)*un+(1.0-fintp(k,i))*up !Central Difference
            ufnorm(k,i)=ufnorm(k,i)+flux(k,i)*schmdefcor(uvnorm)*(uvf-uc)   !Deferred normal flux
          endif
        endif
        if(abs(dndvn)>1.0e-10)then
          uvnorm=1.0-0.5*(vd-vc)/dndvn  !Normalized
          if(uvnorm>0.0 .and. uvnorm<1.0)then
            uvf=fintp(k,i)*vn+(1.0-fintp(k,i))*vp !Central Difference
            vfnorm(k,i)=vfnorm(k,i)+flux(k,i)*schmdefcor(uvnorm)*(uvf-vc)   !Deferred normal flux
          endif
        endif
        ufnorm(jcn,nck)=-ufnorm(k,i)
        vfnorm(jcn,nck)=-vfnorm(k,i)
      enddo  
    enddo
!$OMP END DO
!$OMP DO PRIVATE(i)
    do i=1,ncells
      ssu(i)=ssu(i)-sum(ufnorm(1:ncface(i),i))
      ssv(i)=ssv(i)-sum(vfnorm(1:ncface(i),i))
    enddo
!$OMP END DO    
!$OMP END PARALLEL
    
    return
    end subroutine defcorgammagradvecold

!**************************************************    
    function isnankind(a) result(nan)
! Portable isnan function for arbitrary precision
!**************************************************
    use prec_def, only: ikind
    
    implicit none
    real(ikind),intent(in) :: a 
    logical :: nan
    
    if(a/=a)then 
      nan = .true. 
    else 
      nan = .false. 
    endif 
    
    end function isnankind
    
!**************************************************    
    function isnansingle(a) result(nan)
! Portable isnan function for arbitrary precision
!**************************************************
    implicit none
    real(4),intent(in) :: a 
    logical :: nan
    
    if(a/=a)then 
      nan = .true. 
    else 
      nan = .false. 
    endif 
    
    end function isnansingle
    
!**************************************************    
    function isnandouble(a) result(nan)
! Portable isnan function for arbitrary precision
!**************************************************
    implicit none
    real(8),intent(in) :: a 
    logical :: nan
    
    if(a/=a)then 
      nan = .true. 
    else 
      nan = .false. 
    endif 
    
    end function isnandouble    
        
!*****************************************************************************    
    function adbk(n,a,b) result(res)
! dot product of two vectors with loop unrolling, arbitrary precision, and
! and OpenMP parallelization (=a*b)
!
!Variables and Notations:
!  n: input integer for calculation limit
!  a: input vector
!  b: input vector
!  e: symbol for equal sign
!  d: symbol for dot product
!  k: indicates arbitrary precision
!
! written by Alex Sanchez, USACE-CHL
!*****************************************************************************
    use prec_def, only: ikind
    
    implicit none
    integer :: i,m
    integer,intent(in) :: n
    real(ikind),intent(in) :: a(*),b(*)
    real(ikind) :: res,temp
 
    temp = 0.0_ikind  !Initialize
    m = mod(n,5)
    do i=1,m
      temp = temp + a(i)*b(i)
    enddo
!$OMP PARALLEL DO PRIVATE(i) REDUCTION(+:temp)
    do i=m+1,n,5
      temp = temp + a(i  )*b(i  ) &
                  + a(i+1)*b(i+1) &
                  + a(i+2)*b(i+2) &
                  + a(i+3)*b(i+3) &
                  + a(i+4)*b(i+4)
    enddo    
!$OMP END PARALLEL DO
    res = temp

    return
    end function adbk
    
!*****************************************************************************    
    function adak(n,a) result(res)
! dot product of two vectors with loop unrolling, arbitrary precision, and
! and OpenMP parallelization
!
!Variables and Notations:
!  n: input integer for calculation limit
!  a: input vector
!  e: symbol for equal sign
!  d: symbol for dot product
!  k: indicates arbitrary precision
!
! written by Alex Sanchez, USACE-CHL
!*****************************************************************************
    use prec_def, only: ikind
    
    implicit none
    integer :: i,m
    integer,intent(in) :: n
    real(ikind),intent(in) :: a(*)
    real(ikind) :: res,temp
 
    temp = 0.0_ikind  !Initialize
    m = mod(n,7)
    do i=1,m
      temp = temp + a(i)*a(i)
    enddo
!$OMP PARALLEL DO PRIVATE(i) REDUCTION(+:temp)
    do i=m+1,n,7
      temp = temp + a(i  )*a(i  ) &
                  + a(i+1)*a(i+1) &
                  + a(i+2)*a(i+2) &
                  + a(i+3)*a(i+3) &
                  + a(i+4)*a(i+4) &
                  + a(i+5)*a(i+5) &
                  + a(i+6)*a(i+6)
    enddo    
!$OMP END PARALLEL DO
    res = temp

    return
    end function adak
        
!*****************************************************************************
    subroutine aesmbpak(n,s,b,a)
! computes a constant times a vector plus a vector (a=s*b+a) with 
! loop unrolling, arbitrary precision, and OpenMP parallelization
!
!Variables and Notations:
!  n: input integer for calculation limit
!  a: input/output vector
!  b: input vector
!  s: input scalar
!  e: symbol for equal sign
!  m: symbol for multiplication sign
!  k: indicates arbitrary precision
!
! written by Alex Sanchez, USACE-CHL
!*****************************************************************************
    use prec_def, only: ikind
    
    implicit none
    integer:: i,m
    integer,intent(in) :: n    
    real(ikind),intent(in) :: s,b(*)
    real(ikind),intent(inout) :: a(*)
  
    if(s==0.0_ikind .or. n<=0) return
    m = mod(n,5)
    do i=1,m
      a(i  ) = a(i  ) + s*b(i  )
    enddo
!$OMP PARALLEL DO PRIVATE(i)
    do i=m+1,n,5
      a(i  ) = a(i  ) + s*b(i  )
      a(i+1) = a(i+1) + s*b(i+1)
      a(i+2) = a(i+2) + s*b(i+2)
      a(i+3) = a(i+3) + s*b(i+3)
      a(i+4) = a(i+4) + s*b(i+4)
    enddo
!$OMP END PARALLEL DO

    return
    end subroutine aesmbpak
    
!*****************************************************************************
    subroutine ceambk(n,a,b,c)
! computes the element-by-element product of two vectors (c=a.*b) with 
! loop unrolling, arbitrary precision, and OpenMP parallelization
!
!Variables and Notations:
!  n: input integer for calculation limit
!  a: input vector
!  b: input/output vector
!  c: output vector
!  e: symbol for equal sign
!  m: symbol for multiplication sign
!  k: indicates arbitrary precision
!
! written by Alex Sanchez, USACE-CHL
!*****************************************************************************
    use prec_def, only: ikind
    
    implicit none
    integer:: i,m
    integer,intent(in) :: n
    real(ikind),intent(in) :: a(*),b(*)
    real(ikind),intent(inout) :: c(*)

    m = mod(n,8)
    do i=1,m
      c(i  ) = a(i  )*b(i  )
    enddo
!$OMP PARALLEL DO PRIVATE(i)    
    do i=m+1,n,8
      c(i  ) = a(i  )*b(i  )
      c(i+1) = a(i+1)*b(i+1)
      c(i+2) = a(i+2)*b(i+2)
      c(i+3) = a(i+3)*b(i+3)
      c(i+4) = a(i+4)*b(i+4)
      c(i+5) = a(i+5)*b(i+5)
      c(i+6) = a(i+6)*b(i+6)
      c(i+7) = a(i+7)*b(i+7)
    enddo
!$OMP END PARALLEL DO

    return
    end subroutine ceambk
    
!*****************************************************************************
    subroutine aeambk(n,a,b)
! Element-by-element product of two vectors (a=a.*b) with loop unrolling, 
! arbitrary precision, and OpemMP parallelization.
!
!Variables and Notations:
!  n: input integer for calculation limit
!  a: input/output vector
!  b: input vector
!   e: symbol for equal sign
!   m: symbol for multiplication sign
!   k: symbol for arbitrary precision
!
! written by Alex Sanchez, USACE-CHL
!*****************************************************************************
    use prec_def, only: ikind
    
    implicit none
    integer:: i,m
    integer,intent(in) :: n    
    real(ikind),intent(in) :: b(*)
    real(ikind),intent(inout) :: a(*)

    m = mod(n,8)
    do i=1,m
      a(i  ) = a(i  )*b(i  )
    enddo
!$OMP PARALLEL DO PRIVATE(i)    
    do i=m+1,n,8
      a(i  ) = a(i  )*b(i  )
      a(i+1) = a(i+1)*b(i+1)
      a(i+2) = a(i+2)*b(i+2)
      a(i+3) = a(i+3)*b(i+3)
      a(i+4) = a(i+4)*b(i+4)
      a(i+5) = a(i+5)*b(i+5)
      a(i+6) = a(i+6)*b(i+6)
      a(i+7) = a(i+7)*b(i+7)
    enddo
!$OMP END PARALLEL DO

    return
    end subroutine aeambk

!*****************************************************************************
    subroutine aeapbk(n,a,b)
! Element-by-element sum of two vectors (a=a+b) with loop unrolling, 
! arbitrary precision, and OpemMP parallelization.
!
!Variables and Notations:
!  n: input integer for calculation limit
!  a: input/output vector
!  b: input vector
!   e: symbol for equal sign
!   p: symbol for sumation sign
!   k: symbol for arbitrary precision
!
! written by Alex Sanchez, USACE-CHL
!*****************************************************************************
    use prec_def, only: ikind
    
    implicit none
    integer:: i,m
    integer,intent(in) :: n    
    real(ikind),intent(in) :: b(*)
    real(ikind),intent(inout) :: a(*)

    m = mod(n,8)
    do i=1,m
      a(i  ) = a(i  ) + b(i  )
    enddo
!$OMP PARALLEL DO PRIVATE(i)    
    do i=m+1,n,8
      a(i  ) = a(i  ) + b(i  )
      a(i+1) = a(i+1) + b(i+1)
      a(i+2) = a(i+2) + b(i+2)
      a(i+3) = a(i+3) + b(i+3)
      a(i+4) = a(i+4) + b(i+4)
      a(i+5) = a(i+5) + b(i+5)
      a(i+6) = a(i+6) + b(i+6)
      a(i+7) = a(i+7) + b(i+7)
    enddo
!$OMP END PARALLEL DO

    return
    end subroutine aeapbk
        
!*****************************************************************************
    subroutine aesmak(n,s,a)
! Product of a scalar time a vector (a=s*a) with loop unrolling, 
! arbitrary precision, and OpemMP parallelization.
!
!Variables and Notations:
!  a: input/output vector
!  s: input scalar
!  e: symbol for equal sign
!  m: symbol for multiplication sign
!  k: symbol for arbitrary precision
!
!Notes:
! Variable dimensions are indicated by as the variable size minus 1
! (i.e. single character = scalar, double character = vector, 
!  and triple character = matrix)
!
! written by Alex Sanchez, USACE-CHL
!*****************************************************************************
    use prec_def, only: ikind
    
    implicit none
    integer:: i,m
    integer,intent(in) :: n
    real(ikind),intent(in) :: s
    real(ikind),intent(inout) :: a(*)

    if(s==0.0_ikind)then
      a(1:n) = 0.0_ikind
      return
    endif
    m = mod(n,8)
    do i=1,m
      a(i  ) = s*a(i  )
    enddo
!$OMP PARALLEL DO PRIVATE(i)    
    do i=m+1,n,8
      a(i  ) = s*a(i  )
      a(i+1) = s*a(i+1)
      a(i+2) = s*a(i+2)
      a(i+3) = s*a(i+3)
      a(i+4) = s*a(i+4)
      a(i+5) = s*a(i+5)
      a(i+6) = s*a(i+6)
      a(i+7) = s*a(i+7)
    enddo
!$OMP END PARALLEL DO

    return
    end subroutine aesmak

!*****************************************************************************
    subroutine aesmbk(n,s,b,a)
! Product of a scalar time a vector (a=s*a) with loop unrolling, 
! arbitrary precision, and OpemMP parallelization.
!
!Variables and Notations:
!  a: output vector
!  b: input vector
!  s: input scalar
!  e: symbol for equal sign
!  m: symbol for multiplication sign
!  k: symbol for arbitrary precision
!
!Notes:
! Variable dimensions are indicated by as the variable size minus 1
! (i.e. single character = scalar, double character = vector, 
!  and triple character = matrix)
!
! written by Alex Sanchez, USACE-CHL
!*****************************************************************************
    use prec_def, only: ikind
    
    implicit none
    integer:: i,m
    integer,intent(in) :: n
    real(ikind),intent(in) :: s,b(*)
    real(ikind),intent(out) :: a(*)

    if(s==0.0_ikind)then
      a(1:n) = 0.0_ikind
      return
    endif
    m = mod(n,8)
    do i=1,m
      a(i  ) = s*b(i  )
    enddo
!$OMP PARALLEL DO PRIVATE(i)    
    do i=m+1,n,8
      a(i  ) = s*b(i  )
      a(i+1) = s*b(i+1)
      a(i+2) = s*b(i+2)
      a(i+3) = s*b(i+3)
      a(i+4) = s*b(i+4)
      a(i+5) = s*b(i+5)
      a(i+6) = s*b(i+6)
      a(i+7) = s*b(i+7)
    enddo
!$OMP END PARALLEL DO

    return
    end subroutine aesmbk
    
!*****************************************************************************
    function sumabsk(n,a) result(res)
! Sum of the absolute value of a vector (=sum(abs(a)) with loop unrolling, 
! arbitrary precision, and OpemMP parallelization.
!
!Variables and Notations:
!  a: input/output vector
!  s: input scalar
!  e: symbol for equal sign
!  m: symbol for multiplication sign
!  k: symbol for arbitrary precision
!
!Notes:
! Variable dimensions are indicated by as the variable size minus 1
! (i.e. single character = scalar, double character = vector, 
!  and triple character = matrix)
!
! written by Alex Sanchez, USACE-CHL
!*****************************************************************************
    use prec_def, only: ikind
    
    implicit none
    integer:: i,m
    integer,intent(in) :: n
    real(ikind),intent(in) :: a(*)
    real(ikind) :: temp,res
    
    temp = 0.0_ikind
    m = mod(n,7)
    do i=1,m
      temp = temp + abs(a(i))
    enddo
!$OMP PARALLEL DO PRIVATE(i) REDUCTION(+:temp)
    do i=m+1,n,7
      temp = temp + abs(a(i  )) &
                  + abs(a(i+1)) &
                  + abs(a(i+2)) &
                  + abs(a(i+3)) &
                  + abs(a(i+4)) &
                  + abs(a(i+5)) &
                  + abs(a(i+6))
    enddo    
!$OMP END PARALLEL DO
    res = temp

    return
    end function sumabsk
        
!*****************************************************************************
    subroutine aesumbmck(n,n2,b,c,a)
! Sum of the absolute value of a vector (a(1:n)=sum(b(1:n,1:n2)*c(1:n2))) 
! with loop unrolling, arbitrary precision, and OpemMP parallelization.
!
!Variables and Notations:
!  a: input/output vector
!  b: input vector
!  c: input vector
!  e: symbol for equal sign
!  m: symbol for multiplication sign
!  k: symbol for arbitrary precision
!
! written by Alex Sanchez, USACE-CHL
!*****************************************************************************
    use prec_def, only: ikind
    
    implicit none
    integer:: j
    integer,intent(in) :: n,n2
    real(ikind),intent(in) :: b(n,n2),c(n2)
    real(ikind),intent(out):: a(n)
    
    call aesmbk(n,c(1),b(:,1),a) !a=s*a
    do j=2,n2
      call aesmbpak(n,c(j),b(:,j),a) !a=s*b+a
    enddo
    
    return
    end subroutine aesumbmck
        
!***********************************************
    subroutine smooth1d(niter,m,n,x)
! Smooths a 1D scalar dataset
! niter = number of iterations    
! m = width of moving average    
!***********************************************
    use prec_def, only: ikind
    
    implicit none
    !Input/Output
    integer, intent(in) :: niter  !Number of iterations
    integer, intent(in) :: m      !Window width
    integer, intent(in) :: n      !Calculation size for x 
    real(ikind), intent(inout) :: x(*)
    !Internal variables
    integer :: i
    real(ikind) :: w(m)
    
    !Weights for moving average
    w = 1.0_ikind/real(m,kind=ikind)
    
    !Convolution
    do i=1,niter
      call conv1d(m,w,n,x)
    enddo  
    
    return
    end subroutine smooth1d
    
!***********************************************
    subroutine smoothampdir(numiter,m,n,a,d)
! Smooths 1D data with amplitude a and direction d
! niter = number of iterations    
! m = width of moving average
! a = amplitude    
! d = direction in radians
!***********************************************
    use prec_def, only: ikind
    
    implicit none
    !Input/Output
    integer, intent(in) :: numiter  !Number of iterations
    integer, intent(inout) :: m   !Width of moving average
    integer, intent(in) :: n      !Calculation size for x 
    real(ikind), intent(inout) :: a(*),d(*)
    !Internal variables
    integer :: i
    real(ikind) :: c(n),s(n)    
    
    !Calculate components
    do i=1,n
      c(i) = a(i)*cos(d(i))
      s(i) = a(i)*sin(d(i))    
    enddo
    
    !Smooth each component using a convolution
    call moving_average(numiter,m,n,c)
    call moving_average(numiter,m,n,s)
    
    !Recalculate amplitude and direction
    do i=1,n
      a(i)=sqrt(c(i)*c(i)+s(i)*s(i))
      d(i)=atan2(s(i),c(i))
    enddo

    return
    end subroutine smoothampdir

!***********************************************
    subroutine avgampdir(n,a,d,aavg,davg)
! Computes the average of directional data
! niter = number of iterations
! a = amplitude    
! d = direction [rad]
! aavg = average amplitude
! davg = average direction [rad]
!***********************************************
    use const_def, only: pi
    use prec_def,  only: ikind
    
    implicit none
    !Input/Output
    integer, intent(in) :: n      !Calculation size for x 
    real(ikind), intent(in) :: a(*),d(*)
    real(ikind), intent(out) :: aavg,davg  
    !Internal variables
    integer :: i
    real(ikind) :: cavg,savg   
    
    cavg = 0.0; savg = 0.0
    do i=1,n
      cavg = cavg + a(i)*cos(d(i))
      savg = savg + a(i)*sin(d(i))    
    enddo
    cavg = cavg/real(n,kind=ikind)
    savg = savg/real(n,kind=ikind)
    aavg=sqrt(cavg*cavg+savg*savg)
    davg=atan2(savg,cavg)
    if(davg<0.0)then
      davg = davg + 2.0*pi
    elseif(davg>2.0*pi)then
      davg = davg - 2.0*pi
    endif

    return
    end subroutine avgampdir
    
!***********************************************
    subroutine conv1d(m,w,n,x)
!Numerical computation of convolution in 1D
!    w = Weights with length m
!    x = Vector with length >= n
!    m = length(w)
!    n = length(x)
! written by Alex Sanchez, USACE-CHL
!***********************************************
    use prec_def, only: ikind
    
    implicit none
    !Input/Output
    integer, intent(in) :: m,n
    real(ikind), intent(in) :: w(m)
    real(ikind), intent(inout) :: x(*)    
    !Internal variables
    integer :: i,j,ns,nn,js,je
    real(ikind) :: wsuminv
    real(ikind), allocatable :: y(:)    
    
    nn = m+n-1
    allocate(y(nn))
    do i=1,nn
      y(i) = 0.0
      js = max(1,i+1-n)
      je = min(m,i)
      wsuminv = 1.0/sum(w(js:je))
      do j=js,je
        y(i) = y(i) + w(j)*x(i+1-j)*wsuminv
      enddo
    enddo
    
    if(mod(m,2)==1)then !odd
      ns = (m-1)/2
    else !even  
      ns = m/2
    endif

    do i=1,n
      x(i) = y(ns+i)
    enddo
    
    return
    end subroutine conv1d

!***********************************************
    subroutine moving_average(niter,m,n,x)
! Smooths a 1D scalar dataset
! niter = number of iterations    
! m = width of moving average
! written by Alex Sanchez, USACE-CHL
!***********************************************
    use prec_def, only: ikind
    
    implicit none
    !Input/Output
    integer, intent(in) :: niter  !Number of iterations
    integer, intent(inout) :: m   !Window width
    integer, intent(in) :: n      !Calculation size for x 
    real(ikind), intent(inout) :: x(*)
    !Internal variables
    integer :: i,iw,k,mph,mmh,nw
    real(ikind) :: x2(n)
    
    if(mod(m,2)==0)then
      m = m + 1
    endif
    mph = (m+1)/2
    mmh = (m-1)/2
    do k=1,niter
      do i=1,n
        if(i<=mph)then
            iw=i-1
        elseif(i>n-mph)then
            iw = n-i
        else
            iw = mmh
        endif
        nw=2*iw+1
        x2(i) = sum(x(i-iw:i+iw))/real(nw,kind=ikind)
      enddo
      x(1:n) = x2(1:n)
    enddo    
    
    return
    end subroutine moving_average
       
!**********************************************************************
    subroutine rotate_vector(nsize,ncalc,theta,vecx,vecy)
! Rotates a vector (vecx,vecy) by theta degrees clockwise
!***********************************************************************
    use const_def, only: deg2rad
    use prec_def,  only: ikind
    
    implicit none
    !Input/Output
    integer, intent(in) :: nsize,ncalc
    real(ikind),intent(in) :: theta
    real(ikind),intent(inout) :: vecx(nsize),vecy(nsize)
    !Internal Variables
    integer :: i
    real(ikind) :: costheta,sintheta,tempx,tempy
    
    if(abs(theta)<1.0e-4) return
    
    costheta=cos(theta*deg2rad)
    sintheta=sin(theta*deg2rad)
!$OMP PARALLEL DO PRIVATE(i,tempx,tempy)
    do i=1,ncalc
      tempx= vecx(i)*costheta+vecy(i)*sintheta
      tempy=-vecx(i)*sintheta+vecy(i)*costheta
      vecx(i)=tempx
      vecy(i)=tempy
    enddo
!$OMP END PARALLEL DO

    return 
    end subroutine rotate_vector
    
!**********************************************************************
    subroutine pca(nsize,ncalc,vx,vy,ap,am,s11,s22,Rxy)
! Rotates a vector (vecx,vecy) by theta degrees clockwise
!***********************************************************************
    use const_def, only: pi
    use prec_def,  only: ikind
    
    implicit none
    !Input/Output
    integer, intent(in) :: nsize,ncalc
    real(ikind),intent(in) :: vx(nsize),vy(nsize)    
    real(ikind),intent(out) :: ap,am,s11,s22,Rxy  
    !Internal Variables
    real(ikind) :: vxm,vym,sxx,syy,sxy,fac
        
    !Mean    
    vxm = sum(vx(1:ncalc))/ncalc
    vym = sum(vy(1:ncalc))/ncalc
    
    !Variance and covariance
    fac = 1.0/(ncalc-1)
    sxx = fac*sum((vx(1:ncalc)-vxm)**2)
    syy = fac*sum((vy(1:ncalc)-vym)**2)
    sxy = fac*sum((vx(1:ncalc)-vxm)*(vy(1:ncalc)-vym))
    
    !Angles
    ap = 0.5*atan2(2*sxy,(sxx-syy)) !principle
    am = ap + pi !minor

    !Correlation coefficient
    Rxy = sxy/(sxx**0.5*syy**0.5)

    !Basis vectors (variances parallel and perpendicular to principle axis)
    s11 = 0.5*((sxx+syy) + sqrt((sxx-syy)**2+4*sxy**2))
    s22 = 0.5*((sxx+syy) - sqrt((sxx-syy)**2+4*sxy**2))
    !s12 = 0 %principle covariance is always zero

    return 
    end subroutine pca    