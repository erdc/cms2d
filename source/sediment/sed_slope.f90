!========================================================================
! CMS Bed Slope Routines
!
!  Contains the following:
!    bedslope - Calculates bed-slope term for the bed-change equation
!    bedslopecor_dey - Corrects the critical shear stress/shields 
!                      parameter for the bedslope
!    bedslopecor_bailard - Bed slope correction to bed and suspended 
!                      concentrations based on the work of Bailard (1981)
!
! written by Alex Sanchez, USACE-CHL
!========================================================================    

!***********************************************************************    
    subroutine bedslope
! Calculates bed-slope term for the bed-change equation
!
! written by Alex Sanchez, USACE-CHL; Weiming Wu, NCCHE
!***********************************************************************    
    use size_def
    use geo_def, only: ncface,cell2cell,llec2llec,dsxy,areap,zb,nxyface,kxyface
    use interp_def
    use flow_def
    use bnd_def
    use sed_def
    use prec_def
    implicit none
    integer :: i,j,k,ks,ih,ibnd,nck
    real(ikind) :: asum,qbkip !,fp,

    do ks=1,nsed
      !$omp parallel 
      !Bed-load
      !$omp do private(i)
      do i=1,ncells
        qbk(i,ks) = h(i)*uv(i)*(1.0-rsk(i,ks))*Ctk(i,ks) !kg/m/s
      enddo
      !$omp end do
      
      !Bed-slope coefficients with wetting and drying
      !$omp do private(i,j,k,nck,qbkip)
      do i=1,ncells
        do j=1,nxyface(i)
          k=kxyface(j,i)  
        !do k=1,ncface(i)
          nck=cell2cell(k,i)
          !!Upslope value (more realistic)
          !if(zb(nck)>zb(i) .and. iwet(nck)==1)then
          !  qbkip=qbk(nck,ks)
          !elseif(iwet(i)==1)then
          !  qbkip=qbk(i,ks)  
          !else
          !  qbkip=0.0  
          !endif
          !Harmonic mean
          if(qbk(nck,ks)>1.0e-10 .and. qbk(i,ks)>1.0e-10)then
            qbkip=1.0/(fintp(k,i)/qbk(nck,ks)+(1.0-fintp(k,i))/qbk(i,ks))
          else
            qbkip=0.0  
          endif
          !qbkip=2.0*qbk(nck,ks)*qbk(i,ks)/(qbk(nck,ks)+qbk(i,ks))
          !!Linear interpolation
          !qbkip=fintp(k,i)*qbk(nck,ks)+(1.0-fintp(k,i))*qbk(i,ks)
          acoef(k,i)=iwet(i)*iwet(nck)*dsxy(k,i)*qbkip
          acoef(llec2llec(k,i),nck)=acoef(k,i)
        enddo !k
      enddo !i      
      !$omp end do
      !$omp end parallel
      
      !Apply at all boundary conditions
      do ibnd=1,nbndstr
        do j=1,bnd_str(ibnd)%ncells
          i=bnd_str(ibnd)%cells(j)
          k=bnd_str(ibnd)%faces(j)  
          acoef(k,i)=0.0
        enddo
      enddo          

      !Apply hardbottom
      do ih=1,nhard
        i=idhard(ih)
        if(zb(i)<=hardzb(i))then
          do k=1,ncface(i)
            nck=cell2cell(k,i)
            if(zb(nck)<=hardzb(nck) .or. zb(i)>zb(nck))then
              acoef(k,i)=0.0
              acoef(llec2llec(k,i),nck)=0.0   ! This is important.  Wu
            endif
          enddo !k    
        endif
      enddo !ih

      !Calculate bed-slope term
      !$omp parallel do private(i,k,nck,asum)
      do i=1,ncells
        Sb(i,ks)=0.0; asum=0.0
        do k=1,ncface(i)
          nck=cell2cell(k,i)
          Sb(i,ks)=Sb(i,ks)+acoef(k,i)*zb1(nck) !changed zb to zb1
          !  if (i.eq.487) write(*,*)'in sed_slope k, acoef(k,i),Sb(i,ks)',k,acoef(k,i),Sb(i,ks)        !bdj 01/2021
          !          Sb(i,ks)=Sb(i,ks)+acoef(k,i)*(zb1(nck)+dzb(nck)) !estimate zb
          asum=asum+acoef(k,i)
        enddo !k
        Sb(i,ks)=Sb(i,ks)-asum*zb1(i) !changed zb to zb1
!        Sb(i,ks)=Sb(i,ks)-asum*(zb1(i)+dzb(i)) !estimate zb
        Sb(i,ks)=dcoeff*Sb(i,ks)/areap(i)
      enddo !i   
      !$omp end parallel do
    enddo !ks
    
    return
    end subroutine bedslope    
    
!**************************************************************    
    subroutine bedslopecor_dey
! Corrects the critical shear stress/shields 
! parameter for the bedslope 
!
! Description:
!   Calculates a correction to the critical shear stress/shields
!   parameter using the method of Dey (2001) with recalibrated 
!   coefficients from Walstra et. al. (2007).
!   The correction is added and multiplied by the correction
!   varsigma(i,ks).
!
! References:
!   Dey, S. (2001). "Experimental study on incipient motion of 
!     sediment particles on generalized sloping fluvial beds". 
!     International Journal of Sediment Research, 16(3), 391-398.
!  Walstra, D.J.R., Van Rijn, L.C., Ormondt, M.V., Brière, C., 
!    and Talmon, A.M. (2007). "The effects of bed slope and 
!    wave skewness on sediment transport and morphology", 
!    Proceedings of Coastal Sediments '07, 14 pp.
!
! written by Alex Sanchez, USACE-CHL
!**************************************************************  
    use size_def, only: ncells
    use geo_def, only: dzbx,dzby
    use flow_def, only: u,v,uv
    use sed_def, only: a_repose,varsigma,nsed    
    use sed_lib, only: critslpcor_dey
    use prec_def
    implicit none
    integer :: i
    real(ikind) sl,st,slpcor
    
!$OMP PARALLEL DO PRIVATE(i,sl,st,slpcor)
    do i=1,ncells
      if(uv(i)<=1.0e-4) cycle
      sl = (dzbx(i)*u(i)+dzby(i)*v(i))/max(uv(i),1.0e-5) !Longitudinal slope
      st = (dzby(i)*u(i)-dzbx(i)*v(i))/max(uv(i),1.0e-5)  !Transverse slope
      slpcor = critslpcor_dey(sl,st,a_repose)
      varsigma(i,:)=slpcor*varsigma(i,:) !Both corrections are included in critial shear
    enddo
!$OMP END PARALLEL DO

    return
    end subroutine bedslopecor_dey
    
!***************************************************************    
    subroutine bedslopecor_bailard
! Bed slope correction to bed and suspended concentrations
! based on the work of Bailard (1981)
!
! Description:
!   The method of Bailard (1981) calculates a bed slope correction
!   of the bed and suspended load transport rates. Here the method
!   is applied to the total-load concentration capacity by
!   splitting it into its bed and suspended load components. 
!   This is the correction method proposed for Lund-CIRP transport
!   formula by Carmenen and Larson (2006) and also for the 
!   Soulsby-van Rijn transport formula by Soulsby (1997)
!
! References:
!   Bailard, J.A. (1981). "An energetic total load sediment 
!     transport model for a plane sloping beach".
!     Journal of Geophysical Research, 86(C11), 10938-10954.
!   Camenen, B., and Larson, M. (2007). "A unified sediment 
!    transport formula-tion for coastal inlet application". 
!    Technical report ERDC/CHL CR-07-1, US Army Engineer 
!    Research and Development Center, Vick-sburg, MS.
!  Soulsby, R.L. (1997). "Dynamics of marine sands". 
!    Thomas Telford, London.
!
! written by Alex Sanchez, USACE-CHL
!***************************************************************
    use size_def, only: ncells
    use geo_def, only: dzbx,dzby
    use flow_def, only: u,v,uv
    use fric_def, only: uelwc
    use sed_def, only: nsed,varsigma,CtstarP,rsk,&
       betaslope,effslope,wsfall
    use prec_def
    implicit none
    integer :: i,ks
    real(ikind) :: slp,fb,fs,CbstarP,CsstarP,val

!$OMP PARALLEL DO PRIVATE(i,ks,slp,fb,fs,CbstarP,CsstarP,val)    
    do i=1,ncells
      if(uv(i)<=1.0e-4) cycle
      !Note that it does not include a correction for the transverse slope
      slp = (dzbx(i)*u(i)+dzby(i)*v(i))/uv(i) !streamwise slope, positive upslope
      fb = 1.0-betaslope*slp !Bedload correction
      val = effslope*uelwc(i)*slp !Note: uelwc is combined wave and current velocity
      do ks=1,nsed
        if(CtstarP(i,ks)<1.0e-5) cycle
        fs = 1.0-val/wsfall(ks) !Suspended load correction
        !Split into bed and suspended load and apply corrections
        CsstarP = fs*rsk(i,ks)*CtstarP(i,ks)
        CbstarP = fb*(1.0-rsk(i,ks))*CtstarP(i,ks)
        !Update variables
        CtstarP(i,ks) = CbstarP + CsstarP
        rsk(i,ks) = CsstarP/max(CtstarP(i,ks),1.0e-6) !Fraction of suspended sediment
      enddo !ks
    enddo !i
!$OMP END PARALLEL DO

    return
    end subroutine bedslopecor_bailard    
