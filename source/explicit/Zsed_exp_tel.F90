!***********************************************************************
    subroutine sed_exp_tel
! Solves single and multi-sized sediment transport 
! written by Alex Sanchez, USACE-CHL;  Weiming Wu, NCCHE
!***********************************************************************
    !USE EXP_bndcond_def  
    USE EXP_Global_def, only: FUU,GVV,ACTIVE,ETAN
    use EXP_TELESCOPING,   only: numxfaces, numyfaces, numtbxfaces, numtbyfaces, numregxfaces, numregyfaces
    use EXP_TELESCOPING,   only: xsedtransq, xface_flux, xface_cells, xface_length, xface_wall, xface_wet
    use EXP_TELESCOPING,   only: ysedtransq, yface_flux, yface_cells, yface_length, yface_wall, yface_wet
    use EXP_TELESCOPING,   only: tbxfaces, tbyfaces, regxfaces, regyfaces, cellfaces
    use EXP_transport_def, only: adss,tsed_elapse
    use size_def, only: ncells, ncellsd
    use geo_def,  only: idmap,zb,zbk,dzbx,dzby
    use flow_def, only: vis, eta, h
    use sed_def,  only: zb1, ihidexpform, ibedslope, icapac, nsed, ctstarp, tstartmorph, scalemorph, calcmorph, ibt
    use sed_def,  only: errctk0, errpbk0, ctk, schmidt, ctkstar, alphat, wsfall, btk, do_bedslope, wavesedtrans
    use sed_def,  only: singlesize, d50, d90, dzb, do_aval, sedbalance
    use comvarbl, only: timehrs, dtime
    use der_def,  only: nder,nlim, goa            !added goa 
    use der_lib,  only: der_grad_eval             !renamed from 'gradxy'
    use prec_def, only: ikind
    use geo_def,  only: dx,dy,cell2cell  
    use diag_def, only: dgunit, dgfile
    use interp_lib, only: interp_scal_cell2face
      
    implicit none
    integer :: i,j,ks,ncs,ncw,ii,IDO,nce,ncn,id1,id2
    real(ikind) :: diffT,dxt,dyt,area,voln,fluxT
    real(ikind) :: Dflux(ncellsD)
    
741 format(' Sediment: iter      Ctk_res       errpbk')
751 format(5x,I10,1x,2(E13.4,1x))
349 format(' Sediment: iter      Ctk_res')
359 format(5x,I10,1x,E13.4,1x)
    
    !set zb1 = zb, which is needed for some implicit sub-routiensused here to work properly
    zb1 = zb
      
    !=== Hiding and exposure function ====
    select case(iHidExpForm) !varsigma(i,k)
    case(1); call HidExpEgiazaroff    !Egiazaroff (1965)
    case(2); call HidExpParker        !Parker et al. (1982) and others
    case(3); call HidExpWu            !Wu et al. (2000)
    case(4); call HidExpAshidaMichiue !Ashida and Michiue 1980
    case(5); call HidExpHayashi       !Hayashi et al. 1980
    end select
      
    !=== Incipient Motion Correction for Bedslop ======
    if(ibedslope==1) call bedslopecor_dey
    
    !=== Transport Capacity =====
    select case(icapac)  !CtstarP(i,ks),rs(i,ks)
    case(1); call sedcapac_lundcirp !Lund-CIRP
    case(2); call sedcapac_vanrijn  !Van Rijn
    case(3); call sedcapac_watanabe !Watanabe
    case(4); call sedcapac_soulsby  !Soulsby (1997)
    case(5); call wucapac           !Wu et al. (2000) (under testing)
    end select        
    
    !=== Transport Corrections for Bedslope ===
    if(ibedslope==2) call bedslopecor_bailard
    
    !=== Adjust Concentration Capacity =====
    !$OMP PARALLEL DO PRIVATE(i,ks)
    do i=1,ncells
      !fac=max(0.0,min(1.0,0.01*(zb1(i)-hardzb(i))/d50(i))) !Moved to coeffsourcesink_c 
      !CtstarP(i,:)=CtstarP(i,:)*ramp !*fac !*************
      do ks=1,nsed
        CtstarP(i,ks) = min(CtstarP(i,ks),200.0) !Limit concentration capacity ***************
      enddo
    enddo
    !$OMP END PARALLEL DO
    
    !=== Inflow Boundary Conditions ====
    !call sedbnd_eval  !CtstarP(nck,:)
    
    call Zsedbnd_eval_EXP_tel  
    
    !do i=1,ncells
    !  ctk(i,ks) = CtstarP(i,ks)
    !enddo
       
    !=== Determine whether to calculate morphology change ====
    if(timehrs<tStartMorph .or. scalemorph<1.e-6)then
      calcmorph=.false.
    else
      calcmorph=.true.  
    endif
    
    !=== Adaptation coefficient ==============
    call adaptcoef !alphat(i)           
    
    !=== Total load correction factor ========
    select case(ibt) !btk(i,ks)
    case(1); call btklogexp  
    case(2); call btklogrouse 
    end select
    
    !=== Mixing Layer Thickness ===========
    if(nsed>1) call mixing_layer !db(i,1)
      
    !=== Print header to screen ============
    !if(nsed>1  .and. sedcouple)then 
    !  write(*,741)
    !else
    !  write(*,349)
    !endif          
    
    !=== 3D dispersion terms =======================
    !if(q3d) call q3d_sed !s3dxk, s3dyk             
        
    !=== Solve each A-D equation =====================
    errCtk0 = 0.0_ikind; errpbk0 = 0.0_ikind
    !do itersed=1,maxitersed 
    !  errCtk = 0.0_ikind; errpbk = 0.0_ikind
    !  rsCtkmax = 0.0
    !  do ks=1,nsed   
    !    !=== Coefficients and solve for Ctk =======
    !    select case(ndsch) 
    !    case(1); call coeffsourcesink_c(hybridcoef,ks)
    !    case(2); call coeffsourcesink_c(powerlawcoef,ks)
    !    case(3); call coeffsourcesink_c(exponentialcoef,ks)
    !    case default; call coeffsourcesink_c(upwindcoef,ks)
    !    end select
    !
    !    if(ncellsimple==ncells)then !No gradients required
    !      select case(ndsch) !Anti-diffusion corrections
    !      case(4); call defcorhlpa(Ctk(:,ks),su)
    !      case(5); call defcorgamma(gammadefcor,Ctk(:,ks),su)
    !      case(6); call defcorgamma(cubistadefcor,Ctk(:,ks),su)
    !      case(7); call defcorgamma(alvsmartdefcor,Ctk(:,ks),su)
    !      case(8); call defcorgamma(hoabdefcor,Ctk(:,ks),su)
    !      end select
    !    else
    !      select case(ndsch) !Anti-diffusion corrections
    !      case(4); call defcorhlpagrad(Ctk(:,ks),dCtkx(:,ks),dCtky(:,ks),su)
    !      case(5); call defcorgammagrad(gammadefcor,Ctk(:,ks),dCtkx(:,ks),dCtky(:,ks),su)
    !      case(6); call defcorgammagrad(cubistadefcor,Ctk(:,ks),dCtkx(:,ks),dCtky(:,ks),su)
    !      case(7); call defcorgammagrad(alvsmartdefcor,Ctk(:,ks),dCtkx(:,ks),dCtky(:,ks),su)
    !      case(8); call defcorgammagrad(hoabdefcor,Ctk(:,ks),dCtkx(:,ks),dCtky(:,ks),su)
    !      case default; if(skewcor) call defcorparagrad(dCtkx(:,ks),dCtky(:,ks),su) !Skewness correction      
    !      end select
    !    endif
    !    call bound_c(ks)
    !    !call check_variables(4)
    !    call solve(acoef,su,sp,rsCtk,Ctk(:,ks),4)
    !    errCtk=max(rmom(4),errCtk)
    !    rsCtkmax=max(rsCtkmax,rsCtk)
    !    !=== Check for bad values =====
    !    where(Ctk(:,ks)<0.0) Ctk(:,ks)=0.0
    !    where(Ctk(:,ks)>200.0) Ctk(:,ks)=200.0
    !    !=== Compute gradients
    !    call der_grad_eval(Ctk(:,ks),1,nder,nlim,dCtkx(:,ks),dCtky(:,ks))
    !  enddo !ks


    !=== Set Time Average Flows =====================
    !do i=1,ncells
    !ADSS(i)%qx = ADSS(i)%qx/tsed_elapse
    !ADSS(i)%qy = ADSS(i)%qy/tsed_elapse
    !enddo
!$OMP PARALLEL DO 
    do i=1,numxfaces
      xSedTransQ(i) = xSedTransQ(i)/tsed_elapse
    enddo
!$OMP END PARALLEL DO

!$OMP PARALLEL DO 
    do i=1,numyfaces
      ySedTransQ(i) = ySedTransQ(i)/tsed_elapse
    enddo
!$OMP END PARALLEL DO


    !=== Solve each NET A-D equation =====================
    do ks=1,nsed   
!$OMP DO 
      do i=1,numxfaces
        if(xSedTransQ(i) .gt. 0) then
          xface_flux(i) = xSedTransQ(i)*Ctk(xface_cells(2,i),ks)*xface_length(i)
        else
          xface_flux(i) = xSedTransQ(i)*Ctk(xface_cells(1,i),ks)*xface_length(i)
        endif
      enddo
!$OMP END DO

!$OMP DO 
      do i=1,numyfaces
        if(ySedTransQ(i) .gt. 0) then
          yface_flux(i) = ySedTransQ(i)*Ctk(yface_cells(2,i),ks)*yface_length(i)
        else
          yface_flux(i) = ySedTransQ(i)*Ctk(yface_cells(1,i),ks)*yface_length(i)
        endif
      enddo
!$OMP END DO              
                            
      Dflux = 0.0
                            
!Diffusion Exchanges (telescoping cells)   !Chris Reed - 10/20/2016
      do ii=1,numTBXfaces
        i=TBXfaces(ii)
        if(.not. xface_wall(i) .and. xface_wet(i)) then
          id1 = xface_cells(1,i)
          id2 = xface_cells(2,i)      
            
          diffT = ((vis(id1)+vis(id2))/2.)/schmidt  !
          dxT = (dx(id1)+dx(id2))/2.
          area = xface_length(i)*(eta(id1)-zb(id1)+eta(id2)-zb(id2))/2.
          FluxT = -area*diffT*(Ctk(id1,ks)-Ctk(id2,ks))/dxT            
          Dflux(id1) = Dflux(id1) + FluxT
          Dflux(id2) = Dflux(id2) - FluxT           
        endif
      enddo

      do ii=1,numTBYfaces
        i=TBYfaces(ii)
        if(.not. yface_wall(i) .and. yface_wet(i)) then
          id1 = yface_cells(1,i)
          id2 = yface_cells(2,i)
       
          diffT = ((vis(id1)+vis(id2))/2.)/schmidt  !
          dyT = (dx(id1)+dx(id2))/2.
          area = yface_length(i)*(eta(id1)-zb(id1)+eta(id2)-zb(id2))/2.
          FluxT = -area*diffT*(Ctk(id1,ks)-Ctk(id2,ks))/dyT            
          Dflux(id1) = Dflux(id1) + FluxT
          Dflux(id2) = Dflux(id2) - FluxT               
        endif
      enddo
      
!Diffusion Exchanges (regular cells)   !Chris Reed - 10/20/2016
!$omp do private(i,id1,id2,DiffT,dxT,area,FluxT)
      do ii=1,numREGXfaces
        i=REGXfaces(ii) 
         if(.not. xface_wall(i) .and. xface_wet(i)) then
          id1 = xface_cells(1,i)
          id2 = xface_cells(2,i)
       
          diffT = ((vis(id1)+vis(id2))/2.)/schmidt  !
          dxT = (dx(id1)+dx(id2))/2.
          area = xface_length(i)*(eta(id1)-zb(id1)+eta(id2)-zb(id2))/2.
          FluxT = -area*diffT*(Ctk(id1,ks)-Ctk(id2,ks))/dxT            
          Dflux(id1) = Dflux(id1) + FluxT
          Dflux(id2) = Dflux(id2) - FluxT
        endif
      enddo
!$omp end do 

!$omp do private(i,id1,id2,DiffT,dyT,area,FluxT)
      do ii=1,numREGYfaces
        i=REGYfaces(ii) 
        if(.not. yface_wall(i) .and. yface_wet(i)) then
          id1 = yface_cells(1,i)
          id2 = yface_cells(2,i)     
       
          diffT = ((vis(id1)+vis(id2))/2.)/schmidt  !
          dyT = (dx(id1)+dx(id2))/2.
          area = yface_length(i)*(eta(id1)-zb(id1)+eta(id2)-zb(id2))/2.
          FluxT = -area*diffT*(Ctk(id1,ks)-Ctk(id2,ks))/dyT            
          Dflux(id1) = Dflux(id1) + FluxT
          Dflux(id2) = Dflux(id2) - FluxT           
        endif
      enddo
!$omp end do      

!advection + source/sink + update
!$OMP DO PRIVATE (NCE,NCN,VOLN)
      do i=1,ncells
        Ctkstar(i,ks)=CtstarP(i,ks) !*pbk(i,ks,1)       
        ncn = cell2cell(1,i)
        nce = cell2cell(2,i)
        if(active(i,3)) then
          voln = (-zb(i) + etan(i))*dx(i)*dy(i)        
          ADSS(i)%concn = Ctk(i,ks)*ADSS(i)%vol + Dflux(i) &         !Chris Reed - 10/20/2016
            +( xface_flux(cellfaces(8,i))-xface_flux(cellfaces(3,i))  &
            +  xface_flux(cellfaces(7,i))-xface_flux(cellfaces(4,i))  &
            +  yface_flux(cellfaces(6,i))-yface_flux(cellfaces(1,i))  &
            +  yface_flux(cellfaces(5,i))-yface_flux(cellfaces(2,i))  &        
            +  alphat(i)*wsfall(ks)*(Ctkstar(i,ks)-Ctk(i,ks))*dx(i)*dy(i)  ) &
              *tsed_elapse*Btk(i,ks) 
          ADSS(i)%concn = ADSS(i)%concn/voln
          ADSS(i)%vol = voln
        endif
      enddo
!$OMP END DO          

      do i=1,ncells
        if(active(i,3)) ctk(i,ks) = adss(i)%concn
      enddo
    enddo !ks 
         
    !ii=cell2cell(440,2)   
    !write(*,*)ctk(439,1),ctk(440,1),ctkstar(440,1) 
    !write(*,*)Fuu(439),Fuu(440),Fuu(ii)
           
                 
    !=== Bedslope term =================
    if(do_bedslope) call bedslope !Sb(i,ks)
         
    !=== Wave-induced sed transport ============
    if(wavesedtrans) call sed_wave !Sb(i,ks)+Sw(i) 
        
    !=== Bed Change and Sorting ==============
    if(calcmorph)then                   
      if(singlesize)then !no sorting needed
        call bedchange !dzb(i)
      else      
        call bedchangesort !db(i,1),dzb(i),dzbk(i,ks),pbk(i,ks,1)
      endif
      call struct_dzb !dzb(i) and optionally dzbk(i,ks)
      call check_hardbottom
    endif                  
       
    !!=== Print to screen ==============
    !if(mod(itersed,10)==0 .or. itersed==1)then
    !  if(nsed>1 .and. sedcouple)then
    !    write(*,751) itersed,errCtk,errpbk
    !  else
    !    write(*,359) itersed,errCtk
    !  endif
    !endif
       
    !!=== Check for convergence ============   
    !if(itersed>=maxitersed/2 .or. errCtk<1.e-12)then
    !  if(singlesize)then
    !    if(errCtk<=tolCtk) exit   !Error, concentrations
    !    chgCtk=abs(errCtk-errCtk0)  
    !    if(chgCtk<=1.0e-8) exit   !Absolute error change, concentrations
    !    chgCtk=chgCtk/max(errCtk0,tolCtk)
    !    if(chgCtk<=0.01) exit     !Relative error change, concentrations
    !    errCtk0=errCtk
    !  else   
    !    convCtk=(errCtk<=tolCtk)  !Error, concentrations
    !    convpbk=(errpbk<=tolpbk)  !Error, bed composition
    !    if(convCtk .and. convpbk) exit   
    !    chgCtk=abs(errCtk-errCtk0)  !Absolute error change, concentrations
    !    chgpbk=abs(errpbk-errpbk0)  !Absolute error change, bed composition
    !    convCtk=(convCtk .or. chgCtk<=1.e-7)
    !    convpbk=(convpbk .or. chgpbk<=1.0e-6)
    !    if(convCtk .and. convpbk) exit
    !    chgCtk=chgCtk/max(errCtk0,tolCtk) !Relative error change, concentrations
    !    chgpbk=chgpbk/max(errpbk0,tolpbk) !Relative error change, bed composition
    !    convCtk=(convCtk .or. chgCtk<=0.01)
    !    convpbk=(convpbk .or. chgpbk<=0.01)
    !    if(convCtk .and. convpbk) exit
    !    errCtk0=errCtk
    !    errpbk0=errpbk
    !  endif
    !endif
    
    !enddo !iteration loop     
    !rmom(4) = errCtk  
    
    !=== Bed sorting and gradation =============
    if(.not.singlesize .and. calcmorph)then
      call bedgrad !pbk(i,ks,2:nlay)        
    endif  
      
    !=== Write final sediment error estimates ====   
    !open(dgunit,file=dgfile,access='append')                 
    !if(nsed>1 .and. sedcouple)then                  
    !  !write(dgunit,741)
    !  !write(dgunit,751) itersed,errCtk,errpbk  
    !  if(mod(itersed,10)/=0)then        
    !    if(itersed<maxitersed-1)then
    !      !!write(*,751) itersed,errCtk,errpbk  
    !    endif
    !  endif        
    !else
    !  write(dgunit,349)
    !  write(dgunit,359) itersed,errCtk
    !  if(mod(itersed,10)/=0)then        
    !    if(itersed<maxitersed-1)then
    !      !write(*,359) itersed,errCtk
    !    endif 
    !  endif                     
    !endif  
    !close(dgunit)  
    
    !=== Percentiles ================
    if(.not.singlesize)then
      call sedpercentile(50,d50)
      call sedpercentile(90,d90)       
    endif
 
    !==== Update bed elevation ====================
    do i=1,ncells         
      if(abs(dzb(i))>0.5*h(i))then     
        open(dgunit,file=dgfile,access='append')        
        write(dgunit,*)
        write(dgunit,*) 'WARNING: Large bed change'
        close(dgunit) 
        write(*,*)
        write(*,*) 'WARNING: Large bed change'
        call print_sedvar(i,1)
      endif
      !zb(i)=zb1(i)+dzb(i)   !*****************
      zb(i)=zb(i)+dzb(i)*tsed_elapse/dtime   !*****************      
    enddo  
      
    !=== Avalanching =============
    if(do_aval) call avalanche
    
    !=== Bed-slopes ===============
    call der_grad_eval (goa,0,zb,dzbx,dzby)                          !call der_grad_eval(zb,0,nder,0,dzbx,dzby) !Bed-slope     !!Added goa - 11/3/2015 MB
    
    !=== Update bed elevation at cell faces =====
    !call interp_scal_cell2face(zb,0,zbk,dzbx,dzby)
    
    !=== Correct concentrations for depth changes =======
!!$OMP PARALLEL DO PRIVATE(i,val)
    !do i=1,ncells
    !  val=h(i)                        !save old total water depth
    !  h(i)=max(hmin,p(i)/grav-zb(i))  !New total water depth
    !  !if(h(i)>2*hmin .and. val>2*hmin .and. abs(dzb(i))<0.25*hmin)then      
    !  if(h(i)>2*hmin .and. val>2*hmin)then     
    !    Ctk(i,:) = val*Ctk(i,:)/h(i)
    !  endif  
    !enddo
!!$OMP END PARALLEL DO
    
    !=== Update sediment variables ====
    call sed_total  !cmb(i),Ct(i),Ctstar(i),rs(i),qtx(i),qty(i) 
    
    !== Time step statistics ============
    call sed_step_stat
    
    !=== Total Sediment Budget ======                                 
    if(sedbalance .and. calcmorph) call sed_balance   !do_sedbudget ==> sedbalance
      
    !=== Apply boundary condition to zb ====
    call Zbndzb_exp_tel
        
    return
    end subroutine sed_exp_tel
    
!***********************************************************************
    subroutine Zsedbnd_eval_exp_tel
! Applies sediment transport boundaries
! Blocks off dry regions which are not solved
! written by Alex Sanchez, USACE-CHL;  Weiming Wu, NCCHE
!***********************************************************************
    use EXP_Global_def,  only: ncn,nce,ncw,ncs,qx,qy,active
    use EXP_bndcond_def, only: QstringEXP
    use exp_telescoping, only: cellmap    
    use size_def,  only: ncells
    use geo_def,   only: cell2cell, ds
    use flow_def,  only: u, v, h, uv, iwet, flux
    use bnd_def,   only: nhstr, nthstr, nmhstr, nmhvstr, nqstr, h_str, th_str, mh_str, mhv_str, q_str
    use comvarbl,  only: timehrs, ramp
    use prec_def,  only: ikind
    use sed_def,   only: isedinflowbc, ctstarp, facqtotin, ctk, qtotin, rhosed, nsedflux, sed_str, pbk, pbk1,nsed, sedbnd, ctkstar
    use const_def, only: eps
    use diag_def,  only: dgunit, dgfile

    implicit none
    integer :: i,ii,j,k,ks,nck,ntimes,inc,ised,iwse
    real(ikind) :: fac,qstartot,qstarcell,qsedtot,qut,qvt

    
    if(nHstr .gt. 0) then
      do i = 1,nHstr  !for each cell string
        do j=1,H_str(i)%NCells    !for each cell in string
          ii=H_str(i)%Cells(j) 
          ncn = cellmap(1,ii)
          nce = cellmap(3,ii)
          ncs = cellmap(5,ii)    
          ncw = cellmap(7,ii)                
          quT = u(ii)*h(ii)
          qvT = v(ii)*h(ii)
            
          if(isedinflowbc==1)then
            !if(quT .gt. 0.0 .and. ncw .gt. Ncells) then  !inflow
            !  CtstarP(ii,:) = facQtotin*CtstarP(nce,:) !Capacity times loading factor  
            !else  !outflow
            !  CtstarP(ii,:) = CtstarP(ncw,:) !extraoplate from upstream value
            !endif
            !if(quT .le. 0.0 .and. nce .gt. Ncells) then
            !  CtstarP(ii,:) = facQtotin*CtstarP(ncw,:) !Capacity times loading factor   
            !else
            !  CtstarP(ii,:) = CtstarP(nce,:)   !extraoplate from upstream value  
            !endif
            !if(qvT .gt. 0.0 .and. ncs .gt. Ncells) then
            !  CtstarP(ii,:) = facQtotin*CtstarP(ncn,:) !Capacity times loading factor  
            !else
            !  CtstarP(ii,:) = CtstarP(ncs,:  )!extraoplate from upstream value
            !endif
            !if(qvT .le. 0.0 .and. ncn .gt. Ncells)then
            !  CtstarP(ii,:) = facQtotin*CtstarP(ncs,:) !Capacity times loading factor  
            !else
            !  CtstarP(ii,:) = CtstarP(ncn,:)  !extraoplate from upstream value
            !endif
             
            if(ncw .gt. Ncells) then  !potential left side boundary cell
              if(active(nce,3)) then    !it is left side boundary cell   
                if(quT .gt. 0.0) then
                  CtstarP(ii,:) = facQtotin*CtstarP(nce,:)
                elseif(quT .le. 0.0)then
                  CtstarP(ii,:) = CtstarP(nce,:) 
                endif
              else  ! it is a left side corner cell
                if(ncn .gt. Ncells) then  !upper left corner cell
                  CtstarP(ii,:)= 0.5*(CtstarP(ncs,:) + CtstarP(nce,:))
                elseif(ncs .gt. ncells) then  !lower left corner cell
                  CtstarP(ii,:)= 0.5*(CtstarP(ncn,:) + CtstarP(nce,:))
                endif
              endif
            endif
            if(nce .gt. Ncells) then  !potential right side boundary cell
              if(active(ncw,3)) then    !it is right side boundary cell   
                if(quT .lt. 0.0) then
                  CtstarP(ii,:) = facQtotin*CtstarP(ncw,:)
                elseif(quT .ge. 0.0)then
                  CtstarP(ii,:) = CtstarP(ncw,:) 
                endif
              else  ! it is a right side corner cell
                if(ncn .gt. Ncells) then  !upper right corner cell
                  CtstarP(ii,:)=0.5*(CtstarP(ncs,:) + CtstarP(ncw,:))
                elseif(ncs .gt. ncells) then  !lower right corner cell
                  CtstarP(ii,:)=0.5*(CtstarP(ncn,:) + CtstarP(ncw,:))
                endif                
              endif
            endif         
            if(ncs .gt. Ncells) then  !potential bottom boundary cell
              if(active(ncn,3)) then    !it is bottom boundary cell   
                if(qvT .gt. 0.0) then
                  CtstarP(ii,:) = facQtotin*CtstarP(ncn,:)
                elseif(qvT .le. 0.0)then
                  CtstarP(ii,:) = CtstarP(ncn,:) 
                endif
              else  ! it is a bottom corner cell
                if(ncw .gt. Ncells) then  !bottom left corner cell
                  CtstarP(ii,:)=0.5*(CtstarP(ncn,:) + CtstarP(nce,:))
                elseif(nce .gt. ncells) then  !lower right corner cell
                  CtstarP(ii,:)=0.5*(CtstarP(ncn,:) + CtstarP(ncw,:))
                endif  
              endif
            endif
            if(ncn .gt. Ncells) then  !potential bottom boundary cell
              if(active(ncs,3)) then    !it is bottom boundary cell   
                if(qvT .lt. 0.0) then
                  CtstarP(ii,:) = facQtotin*CtstarP(ncs,:)
                elseif(quT .ge. 0.0)then
                  CtstarP(ii,:) = CtstarP(ncs,:) 
                endif
              else  ! it is a top corner cell
                if(ncw .gt. Ncells) then  !top left corner cell
                  CtstarP(ii,:)=0.5*(CtstarP(ncs,:) + CtstarP(nce,:))
                elseif(nce .gt. ncells) then  !top right corner cell
                  CtstarP(ii,:)=0.5*(CtstarP(ncs,:) + CtstarP(ncw,:))
                endif 
              endif
            endif          

            !Ctk(ii,:)=pbk(ii,:,1)*CtstarP(ii,:)
            Ctk(ii,:)=CtstarP(ii,:)             
          elseif(isedinflowbc==2)then  
            if(quT .gt. 0.0 .and. ncw .gt. Ncells) CtstarP(ii,:) = Qtotin/(uv(nce)*h(nce))/rhosed  !Qtotin in kg/m/sec  
            if(quT .le. 0.0 .and. nce .gt. Ncells) CtstarP(ii,:) = Qtotin/(uv(ncw)*h(ncw))/rhosed  !Qtotin in kg/m/sec     
            if(qvT .gt. 0.0 .and. ncs .gt. Ncells) CtstarP(ii,:) = Qtotin/(uv(ncn)*h(ncn))/rhosed  !Qtotin in kg/m/sec    
            if(qvT .le. 0.0 .and. ncn .gt. Ncells) CtstarP(ii,:) = Qtotin/(uv(ncw)*h(ncw))/rhosed  !Qtotin in kg/m/sec   
            
            if(quT .gt. 0.0 .and. ncw .gt. Ncells) then  !inflow
              CtstarP(ii,:) = Qtotin/(uv(nce)*h(nce))/rhosed  !Qtotin in kg/m/sec
            else  !outflow
              CtstarP(ii,:) = CtstarP(ncw,:) !extrapolate from upstream value
            endif
            if(quT .le. 0.0 .and. nce .gt. Ncells) then
              CtstarP(ii,:) = Qtotin/(uv(ncw)*h(ncw))/rhosed  !Qtotin in kg/m/sec  
            else
              CtstarP(ii,:) = CtstarP(nce,:)   !extrapolate from upstream value  
            endif
            if(qvT .gt. 0.0 .and. ncs .gt. Ncells) then
              CtstarP(ii,:) = Qtotin/(uv(ncn)*h(ncn))/rhosed  !Qtotin in kg/m/sec  
            else
              CtstarP(ii,:) = CtstarP(ncs,:  )!extrapolate from upstream value
            endif
            if(qvT .le. 0.0 .and. ncn .gt. Ncells)then
              CtstarP(ii,:) = Qtotin/(uv(ncs)*h(ncs))/rhosed  !Qtotin in kg/m/sec
            else
              CtstarP(ii,:) = CtstarP(ncn,:)  !extrapolate from upstream value
            endif           
            
            !Ctk(ii,:)=pbk(ii,:,1)*CtstarP(ii,:)
            Ctk(ii,:)=CtstarP(ii,:)             
          endif                          
        enddo
      enddo ! end of each cell string
    endif  !nHstr strings

    if(nTHstr .gt. 0) then
      do iwse=1,nTHstr
        do j=1,TH_str(iwse)%NCells    !for each cell in string
          ii=TH_str(iwse)%Cells(j)  !Chris Reed - 10/20/2016
          ncn = cellmap(1,ii)
          nce = cellmap(3,ii)
          ncs = cellmap(5,ii)    
          ncw = cellmap(7,ii)                
          quT = u(ii)*h(ii)
          qvT = v(ii)*h(ii)
          if(isedinflowbc==1)then
            if(ncw .gt. Ncells) then  !potential left side boundary cell
              if(active(nce,3)) then    !it is left side boundary cell   
                if(quT .gt. 0.0) then
                  CtstarP(ii,:) = facQtotin*CtstarP(nce,:)
                elseif(quT .le. 0.0)then
                  CtstarP(ii,:) = CtstarP(nce,:) 
                endif
              else  ! it is a left side corner cell
                if(ncn .gt. Ncells) then  !upper left corner cell
                  CtstarP(ii,:)= 0.5*(CtstarP(ncs,:) + CtstarP(nce,:))
                elseif(ncs .gt. ncells) then  !lower left corner cell
                  CtstarP(ii,:)= 0.5*(CtstarP(ncn,:) + CtstarP(nce,:))
                endif
              endif
            endif
            if(nce .gt. Ncells) then  !potential right side boundary cell
              if(active(ncw,3)) then    !it is right side boundary cell   
                if(quT .lt. 0.0) then
                  CtstarP(ii,:) = facQtotin*CtstarP(ncw,:)
                elseif(quT .ge. 0.0)then
                  CtstarP(ii,:) = CtstarP(ncw,:) 
                endif
              else  ! it is a right side corner cell
                if(ncn .gt. Ncells) then  !upper right corner cell
                  CtstarP(ii,:)=0.5*(CtstarP(ncs,:) + CtstarP(ncw,:))
                elseif(ncs .gt. ncells) then  !lower right corner cell
                  CtstarP(ii,:)=0.5*(CtstarP(ncn,:) + CtstarP(ncw,:))
                endif                
              endif
            endif         
            if(ncs .gt. Ncells) then  !potential bottom boundary cell
              if(active(ncn,3)) then    !it is bottom boundary cell   
                if(qvT .gt. 0.0) then
                  CtstarP(ii,:) = facQtotin*CtstarP(ncn,:)
                elseif(qvT .le. 0.0)then
                  CtstarP(ii,:) = CtstarP(ncn,:) 
                endif
              else  ! it is a bottom corner cell
                if(ncw .gt. Ncells) then  !bottom left corner cell
                  CtstarP(ii,:)=0.5*(CtstarP(ncn,:) + CtstarP(nce,:))
                elseif(nce .gt. ncells) then  !lower right corner cell
                  CtstarP(ii,:)=0.5*(CtstarP(ncn,:) + CtstarP(ncw,:))
                endif  
              endif
            endif
            if(ncn .gt. Ncells) then  !potential bottom boundary cell
              if(active(ncs,3)) then    !it is bottom boundary cell   
                if(qvT .lt. 0.0) then
                  CtstarP(ii,:) = facQtotin*CtstarP(ncs,:)
                elseif(quT .ge. 0.0)then
                  CtstarP(ii,:) = CtstarP(ncs,:) 
                endif
              else  ! it is a top corner cell
                if(ncw .gt. Ncells) then  !top left corner cell
                  CtstarP(ii,:)=0.5*(CtstarP(ncs,:) + CtstarP(nce,:))
                elseif(nce .gt. ncells) then  !top right corner cell
                  CtstarP(ii,:)=0.5*(CtstarP(ncs,:) + CtstarP(ncw,:))
                endif 
              endif
            endif       
            !Ctk(ii,:)=pbk(ii,:,1)*CtstarP(ii,:)
            Ctk(ii,:)=CtstarP(ii,:)             
          elseif(isedinflowbc==2)then  
            if(quT .gt. 0.0 .and. ncw .gt. Ncells) CtstarP(ii,:) = Qtotin/(uv(nce)*h(nce))/rhosed  !Qtotin in kg/m/sec  
            if(quT .le. 0.0 .and. nce .gt. Ncells) CtstarP(ii,:) = Qtotin/(uv(ncw)*h(ncw))/rhosed  !Qtotin in kg/m/sec     
            if(qvT .gt. 0.0 .and. ncs .gt. Ncells) CtstarP(ii,:) = Qtotin/(uv(ncn)*h(ncn))/rhosed  !Qtotin in kg/m/sec    
            if(qvT .le. 0.0 .and. ncn .gt. Ncells) CtstarP(ii,:) = Qtotin/(uv(ncw)*h(ncw))/rhosed  !Qtotin in kg/m/sec   
            
            if(quT .gt. 0.0 .and. ncw .gt. Ncells) then  !inflow
              CtstarP(ii,:) = Qtotin/(uv(nce)*h(nce))/rhosed  !Qtotin in kg/m/sec
            else  !outflow
              CtstarP(ii,:) = CtstarP(ncw,:) !extrapolate from upstream value
            endif
            if(quT .le. 0.0 .and. nce .gt. Ncells) then
              CtstarP(ii,:) = Qtotin/(uv(ncw)*h(ncw))/rhosed  !Qtotin in kg/m/sec  
            else
              CtstarP(ii,:) = CtstarP(nce,:)   !extrapolate from upstream value  
            endif
            if(qvT .gt. 0.0 .and. ncs .gt. Ncells) then
              CtstarP(ii,:) = Qtotin/(uv(ncn)*h(ncn))/rhosed  !Qtotin in kg/m/sec  
            else
              CtstarP(ii,:) = CtstarP(ncs,:  )!extrapolate from upstream value
            endif
            if(qvT .le. 0.0 .and. ncn .gt. Ncells)then
              CtstarP(ii,:) = Qtotin/(uv(ncs)*h(ncs))/rhosed  !Qtotin in kg/m/sec
            else
              CtstarP(ii,:) = CtstarP(ncn,:)  !extrapolate from upstream value
            endif           
            
            !Ctk(ii,:)=pbk(ii,:,1)*CtstarP(ii,:)
            Ctk(ii,:)=CtstarP(ii,:)             
          endif                               
        enddo 
      enddo !iwse-str      
    endif  !H_tide

    if(nMHstr .gt. 0) then
      do i = 1,nMHstr  !for each cell string                        
        do j=1,MH_str(i)%NCells   !for each cell in string
          ii=MH_str(i)%Cells(j) 
          ncn = cellmap(1,ii)
          nce = cellmap(3,ii)
          ncs = cellmap(5,ii)    
          ncw = cellmap(7,ii)                
          quT = u(ii)*h(ii)
          qvT = v(ii)*h(ii)
          if(isedinflowbc==1)then
            if(ncw .gt. Ncells) then  !potential left side boundary cell
              if(active(nce,3)) then    !it is left side boundary cell   
                if(quT .gt. 0.0) then
                  CtstarP(ii,:) = facQtotin*CtstarP(nce,:)
                elseif(quT .le. 0.0)then
                  CtstarP(ii,:) = CtstarP(nce,:) 
                endif
              else  ! it is a left side corner cell
                if(ncn .gt. Ncells) then  !upper left corner cell
                  CtstarP(ii,:)= 0.5*(CtstarP(ncs,:) + CtstarP(nce,:))
                elseif(ncs .gt. ncells) then  !lower left corner cell
                  CtstarP(ii,:)= 0.5*(CtstarP(ncn,:) + CtstarP(nce,:))
                endif
              endif
            endif
            if(nce .gt. Ncells) then  !potential right side boundary cell
              if(active(ncw,3)) then    !it is right side boundary cell   
                if(quT .lt. 0.0) then
                  CtstarP(ii,:) = facQtotin*CtstarP(ncw,:)
                elseif(quT .ge. 0.0)then
                  CtstarP(ii,:) = CtstarP(ncw,:) 
                endif
              else  ! it is a right side corner cell
                if(ncn .gt. Ncells) then  !upper right corner cell
                  CtstarP(ii,:)=0.5*(CtstarP(ncs,:) + CtstarP(ncw,:))
                elseif(ncs .gt. ncells) then  !lower right corner cell
                  CtstarP(ii,:)=0.5*(CtstarP(ncn,:) + CtstarP(ncw,:))
                endif                
              endif
            endif         
            if(ncs .gt. Ncells) then  !potential bottom boundary cell
              if(active(ncn,3)) then    !it is bottom boundary cell   
                if(qvT .gt. 0.0) then
                  CtstarP(ii,:) = facQtotin*CtstarP(ncn,:)
                elseif(qvT .le. 0.0)then
                  CtstarP(ii,:) = CtstarP(ncn,:) 
                endif
              else  ! it is a bottom corner cell
                if(ncw .gt. Ncells) then  !bottom left corner cell
                  CtstarP(ii,:)=0.5*(CtstarP(ncn,:) + CtstarP(nce,:))
                elseif(nce .gt. ncells) then  !lower right corner cell
                  CtstarP(ii,:)=0.5*(CtstarP(ncn,:) + CtstarP(ncw,:))
                endif  
              endif
            endif
            if(ncn .gt. Ncells) then  !potential bottom boundary cell
              if(active(ncs,3)) then    !it is bottom boundary cell   
                if(qvT .lt. 0.0) then
                  CtstarP(ii,:) = facQtotin*CtstarP(ncs,:)
                elseif(quT .ge. 0.0)then
                  CtstarP(ii,:) = CtstarP(ncs,:) 
                endif
              else  ! it is a top corner cell
                if(ncw .gt. Ncells) then  !top left corner cell
                  CtstarP(ii,:)=0.5*(CtstarP(ncs,:) + CtstarP(nce,:))
                elseif(nce .gt. ncells) then  !top right corner cell
                  CtstarP(ii,:)=0.5*(CtstarP(ncs,:) + CtstarP(ncw,:))
                endif 
              endif
            endif       
            !Ctk(ii,:)=pbk(ii,:,1)*CtstarP(ii,:)
            Ctk(ii,:)=CtstarP(ii,:)             
          elseif(isedinflowbc==2)then  
            if(quT .gt. 0.0 .and. ncw .gt. Ncells) CtstarP(ii,:) = Qtotin/(uv(nce)*h(nce))/rhosed  !Qtotin in kg/m/sec  
            if(quT .le. 0.0 .and. nce .gt. Ncells) CtstarP(ii,:) = Qtotin/(uv(ncw)*h(ncw))/rhosed  !Qtotin in kg/m/sec     
            if(qvT .gt. 0.0 .and. ncs .gt. Ncells) CtstarP(ii,:) = Qtotin/(uv(ncn)*h(ncn))/rhosed  !Qtotin in kg/m/sec    
            if(qvT .le. 0.0 .and. ncn .gt. Ncells) CtstarP(ii,:) = Qtotin/(uv(ncw)*h(ncw))/rhosed  !Qtotin in kg/m/sec   
            
            if(quT .gt. 0.0 .and. ncw .gt. Ncells) then  !inflow
              CtstarP(ii,:) = Qtotin/(uv(nce)*h(nce))/rhosed  !Qtotin in kg/m/sec
            else  !outflow
              CtstarP(ii,:) = CtstarP(ncw,:) !extrapolate from upstream value
            endif
            if(quT .le. 0.0 .and. nce .gt. Ncells) then
              CtstarP(ii,:) = Qtotin/(uv(ncw)*h(ncw))/rhosed  !Qtotin in kg/m/sec  
            else
              CtstarP(ii,:) = CtstarP(nce,:)   !extrapolate from upstream value  
            endif
            if(qvT .gt. 0.0 .and. ncs .gt. Ncells) then
              CtstarP(ii,:) = Qtotin/(uv(ncn)*h(ncn))/rhosed  !Qtotin in kg/m/sec  
            else
              CtstarP(ii,:) = CtstarP(ncs,:  )!extrapolate from upstream value
            endif
            if(qvT .le. 0.0 .and. ncn .gt. Ncells)then
              CtstarP(ii,:) = Qtotin/(uv(ncs)*h(ncs))/rhosed  !Qtotin in kg/m/sec
            else
              CtstarP(ii,:) = CtstarP(ncn,:)  !extrapolate from upstream value
            endif           
            
            !Ctk(ii,:)=pbk(ii,:,1)*CtstarP(ii,:)
            Ctk(ii,:)=CtstarP(ii,:)             
          endif                                  
        enddo
      enddo ! end of each cell string
    endif  !H_multi

    if(nMHVstr .gt. 0) then
      do i = 1,nMHVstr  !for each cell string
        do j=1,MHV_str(i)%NCells   !for each cell in string
          ii=MHV_str(i)%Cells(j) 
          ncn = cellmap(1,ii)
          nce = cellmap(3,ii)
          ncs = cellmap(5,ii)    
          ncw = cellmap(7,ii)                
          quT = u(ii)*h(ii)
          qvT = v(ii)*h(ii)
          if(isedinflowbc==1)then
            if(ncw .gt. Ncells) then  !potential left side boundary cell
              if(active(nce,3)) then    !it is left side boundary cell   
                if(quT .gt. 0.0) then
                  CtstarP(ii,:) = facQtotin*CtstarP(nce,:)
                elseif(quT .le. 0.0)then
                  CtstarP(ii,:) = CtstarP(nce,:) 
                endif
              else  ! it is a left side corner cell
                if(ncn .gt. Ncells) then  !upper left corner cell
                  CtstarP(ii,:)= 0.5*(CtstarP(ncs,:) + CtstarP(nce,:))
                elseif(ncs .gt. ncells) then  !lower left corner cell
                  CtstarP(ii,:)= 0.5*(CtstarP(ncn,:) + CtstarP(nce,:))
                endif
              endif
            endif
            if(nce .gt. Ncells) then  !potential right side boundary cell
              if(active(ncw,3)) then    !it is right side boundary cell   
                if(quT .lt. 0.0) then
                  CtstarP(ii,:) = facQtotin*CtstarP(ncw,:)
                elseif(quT .ge. 0.0)then
                  CtstarP(ii,:) = CtstarP(ncw,:) 
                endif
              else  ! it is a right side corner cell
                if(ncn .gt. Ncells) then  !upper right corner cell
                  CtstarP(ii,:)=0.5*(CtstarP(ncs,:) + CtstarP(ncw,:))
                elseif(ncs .gt. ncells) then  !lower right corner cell
                  CtstarP(ii,:)=0.5*(CtstarP(ncn,:) + CtstarP(ncw,:))
                endif                
              endif
            endif         
            if(ncs .gt. Ncells) then  !potential bottom boundary cell
              if(active(ncn,3)) then    !it is bottom boundary cell   
                if(qvT .gt. 0.0) then
                  CtstarP(ii,:) = facQtotin*CtstarP(ncn,:)
                elseif(qvT .le. 0.0)then
                  CtstarP(ii,:) = CtstarP(ncn,:) 
                endif
              else  ! it is a bottom corner cell
                if(ncw .gt. Ncells) then  !bottom left corner cell
                  CtstarP(ii,:)=0.5*(CtstarP(ncn,:) + CtstarP(nce,:))
                elseif(nce .gt. ncells) then  !lower right corner cell
                  CtstarP(ii,:)=0.5*(CtstarP(ncn,:) + CtstarP(ncw,:))
                endif  
              endif
            endif
            if(ncn .gt. Ncells) then  !potential bottom boundary cell
              if(active(ncs,3)) then    !it is bottom boundary cell   
                if(qvT .lt. 0.0) then
                  CtstarP(ii,:) = facQtotin*CtstarP(ncs,:)
                elseif(quT .ge. 0.0)then
                  CtstarP(ii,:) = CtstarP(ncs,:) 
                endif
              else  ! it is a top corner cell
                if(ncw .gt. Ncells) then  !top left corner cell
                  CtstarP(ii,:)=0.5*(CtstarP(ncs,:) + CtstarP(nce,:))
                elseif(nce .gt. ncells) then  !top right corner cell
                  CtstarP(ii,:)=0.5*(CtstarP(ncs,:) + CtstarP(ncw,:))
                endif 
              endif
            endif       
            
            !Ctk(ii,:)=pbk(ii,:,1)*CtstarP(ii,:)
            Ctk(ii,:)=CtstarP(ii,:)             
          elseif(isedinflowbc==2)then  
            if(quT .gt. 0.0 .and. ncw .gt. Ncells) CtstarP(ii,:) = Qtotin/(uv(nce)*h(nce))/rhosed  !Qtotin in kg/m/sec  
            if(quT .le. 0.0 .and. nce .gt. Ncells) CtstarP(ii,:) = Qtotin/(uv(ncw)*h(ncw))/rhosed  !Qtotin in kg/m/sec     
            if(qvT .gt. 0.0 .and. ncs .gt. Ncells) CtstarP(ii,:) = Qtotin/(uv(ncn)*h(ncn))/rhosed  !Qtotin in kg/m/sec    
            if(qvT .le. 0.0 .and. ncn .gt. Ncells) CtstarP(ii,:) = Qtotin/(uv(ncw)*h(ncw))/rhosed  !Qtotin in kg/m/sec   
             
            if(quT .gt. 0.0 .and. ncw .gt. Ncells) then  !inflow
              CtstarP(ii,:) = Qtotin/(uv(nce)*h(nce))/rhosed  !Qtotin in kg/m/sec
            else  !outflow
              CtstarP(ii,:) = CtstarP(ncw,:) !extrapolate from upstream value
            endif
            if(quT .le. 0.0 .and. nce .gt. Ncells) then
              CtstarP(ii,:) = Qtotin/(uv(ncw)*h(ncw))/rhosed  !Qtotin in kg/m/sec  
            else
              CtstarP(ii,:) = CtstarP(nce,:)   !extrapolate from upstream value  
            endif
            if(qvT .gt. 0.0 .and. ncs .gt. Ncells) then
              CtstarP(ii,:) = Qtotin/(uv(ncn)*h(ncn))/rhosed  !Qtotin in kg/m/sec  
            else
              CtstarP(ii,:) = CtstarP(ncs,:  )!extrapolate from upstream value
            endif
            if(qvT .le. 0.0 .and. ncn .gt. Ncells)then
              CtstarP(ii,:) = Qtotin/(uv(ncs)*h(ncs))/rhosed  !Qtotin in kg/m/sec
            else
              CtstarP(ii,:) = CtstarP(ncn,:)  !extrapolate from upstream value
            endif           
            
            !Ctk(ii,:)=pbk(ii,:,1)*CtstarP(ii,:)
            Ctk(ii,:)=CtstarP(ii,:)             
          endif                                 
        enddo
      enddo ! end of each cell string
    endif  !H_multi   
 
    if(nQstr .gt. 0) then
      do i = 1,nQstr  !for each cell string
        if(QstringEXP(i)%vface) then 
          if( QstringEXP(i)%sgn .eq. 1 ) then  !south face 
            do j=1,Q_str(i)%NCells    !for each cell in string
              ii= Q_str(i)%Cells(j) 
              ncs = cellmap(5,ii)
              CtstarP(ncs,:) = facQtotin*CtstarP(ii,:)            
              !Ctk(iid,:)=pbk(iid,:,1)*CtstarP(iid,:)
              Ctk(ncs,:)=CtstarP(ncs,:) 
            enddo
          else  !north face
            do j=1,Q_str(i)%NCells    !for each cell in string
              II = Q_str(i)%Cells(j)
              ncs = cellmap(5,ii)
              CtstarP(ii,:) = facQtotin*CtstarP(ncs,:)
              !Ctk(iid,:)=pbk(iid,:,1)*CtstarP(iid,:)
              Ctk(ii,:)=CtstarP(ii,:) 
            enddo               
          endif
        else
          if( QstringEXP(i)%sgn .eq. 1 ) then  !west face 
            do j=1,Q_str(i)%NCells    !for each cell in string
              ii = Q_str(i)%Cells(j)
              ncw = cellmap(7,ii)
              CtstarP(ncw,:) = facQtotin*CtstarP(ii,:)           
              !Ctk(iid,:)=pbk(iid,:,1)*CtstarP(iid,:)
              Ctk(ncw,:)=CtstarP(ncw,:)            
            enddo
          else  !east face
            do j=1,Q_str(i)%NCells    !for each cell in string
              II = Q_str(i)%Cells(j)   
              ncw = cellmap(7,ii)          
              CtstarP(ii,:) = facQtotin*CtstarP(ncw,:)              
              !Ctk(iid,:)=pbk(iid,:,1)*CtstarP(iid,:)
              Ctk(ii,:)=CtstarP(ii,:) 
            enddo        
          endif  
        endif
      enddo ! end of NQdriver
    endif  !Q_single     
      
    goto 1000
      
  !Sediment flux boundaries, will overide river boundary conditions
    do ised=1,nsedflux
      !find out where we are in the time/value arrays
      inc = sed_str(ised)%inc
      ntimes = sed_str(ised)%ntimes
      do while ((timehrs+eps)>=sed_str(ised)%times(inc+1) .and. inc<ntimes)
        inc = inc + 1
        sed_str(ised)%inc = inc
      enddo
      fac = (timehrs-sed_str(ised)%times(inc))/ &
            (sed_str(ised)%times(inc+1)-sed_str(ised)%times(inc))     
      
      !Get sediment transport rate      
      if(sed_str(ised)%ibctype==1)then !Total sediment transport rate        
        qsedtot = (1.0-fac)*sed_str(ised)%val(inc,1) + &
                  fac*sed_str(ised)%val(inc+1,1)  !kg/sec, total per cellstring                        
        !Compute average bed composition for cell string   
        do ii=1,sed_str(ised)%ncells
          i = sed_str(ised)%cells(ii)  
          k = sed_str(ised)%faces(ii)    
          nck = cell2cell(i,k)
          pbk(nck,:,1) = 0.0
          do ks=1,nsed
            pbk(nck,ks,1) = pbk(nck,ks,1) + pbk1(i,ks)
          enddo          
        enddo    
        pbk(nck,:,1) = pbk(nck,:,1)/float(sed_str(ised)%ncells) 
        !Fractional bed composition based on average bed composition at boundary
        qsedtot = ramp*facQtotin*qsedtot !Apply ramp and loading factor
        sedbnd(ised,:) = pbk(nck,:,1)*qsedtot      
      else      !Fractional sediment transport rate                
        sedbnd(ised,:) = (1.0-fac)*sed_str(ised)%val(inc,:) + &
                         fac*sed_str(ised)%val(inc+1,:) !kg/sec, per fraction and cellstring           
        sedbnd(ised,:) = ramp*facQtotin*sedbnd(ised,:) !Apply ramp and loading factor
        
        !Sum sediment transport rate
        qsedtot = sum(sedbnd(ised,:))
        !!Compute bed composition at ghost cells           
        !if(qsedtot>1.e-5)then 
        !  do ii=1,sed_str(ised)%ncells
        !    i = sed_str(ised)%cells(ii)  
        !    k = sed_str(ised)%faces(ii)    
        !    nck = cell2cell(i,k)                     
        !    pbk(nck,:,1) = sedbnd(ised,:)/qsedtot
        !    pbk(nck,:,1) = pbk(nck,:,1)/sum(pbk(nck,:,1)) !Normalize to make sure sum is 1
        !  enddo
        !else
        !  do ii=1,sed_str(ised)%ncells
        !    i = sed_str(ised)%cells(ii)  
        !    k = sed_str(ised)%faces(ii)    
        !    nck = cell2cell(i,k)   
        !    pbk(nck,:,1) = pbk(i,:,1)
        !  enddo
        !endif  
      endif      
            
      write(*,*) 'Specified Total Inflow Sediment Transport Rate: ',qsedtot,' kg/sec'      
      !write(*,*) 'Fractional Inflow Sediment Transport Rates, mm, kg/sec'      
      !do ks=1,nsed                      
      !  write(*,*) diam(ks)*1000.0, sedbnd(ised,ks)
      !enddo !ks 
     
      qsedtot = 0.0
      do ks=1,nsed   
        qstartot = 0.0   
        do ii=1,sed_str(ised)%ncells
          i = sed_str(ised)%cells(ii)  
          k = sed_str(ised)%faces(ii) 
          nck = cell2cell(k,i)
          if(iwet(i)==0) cycle                                  
          if(flux(k,i)<0.0)then !Inflow      
            !qstarcell = h(i)*uv(i)*CtstarP(i,ks) !kg/m/sec  
            qstarcell = 1.0
            qstarcell = max(qstarcell,0.0001)                                                            
            qstartot = qstartot + ds(k,i)*qstarcell !kg/sec 
          endif 
        enddo !ii  
        
        qstartot = max(qstartot,1.0e-15)
        do ii=1,sed_str(ised)%ncells
          i = sed_str(ised)%cells(ii)  
          k = sed_str(ised)%faces(ii)    
          nck = cell2cell(k,i)  
          if(iwet(i)==0) cycle                                  
          if(flux(k,i)<0.0)then
            !qstarcell = h(i)*uv(i)*CtstarP(i,ks) !kg/m/sec  
            qstarcell = 1.0
            qstarcell = max(qstarcell,0.0001) 
            fac = qstarcell/(h(i)*uv(i)*qstartot)
            Ctkstar(nck,ks) = fac*sedbnd(ised,ks) !convert kg/sec to kg/m^3
            CtstarP(nck,ks) = Ctkstar(nck,ks)/pbk(i,ks,1) !Note i index in pbk, used in bound_c, pbk(i,ks,1)*CtstarP(nck,ks) used as boundary condition
            qsedtot = qsedtot + ds(k,i)*h(i)*uv(i)*Ctkstar(nck,ks)  !kg/sec  
          else
            Ctkstar(nck,:) = Ctkstar(i,:) 
            CtstarP(nck,:) = CtstarP(i,:)
          endif
        enddo !ii    
      enddo !ks    
      
      write(*,*) 'Calculated Total Inflow Sediment Transport Rate: ',qsedtot,' kg/sec'    
      open(dgunit,file=dgfile,access='append') 
      write(dgunit,*) 'Calculated Total Inflow Sediment Transport Rate: ',qsedtot,' kg/sec'   
      close(dgunit)
    enddo !ised

1000 continue    
    
    return
    end subroutine Zsedbnd_eval_EXP_tel
    
!***********************************************************************
    subroutine Zbndzb_exp_tel
! Applies sediment transport boundaries
! Blocks off dry regions which are not solved
! written by Alex Sanchez, USACE-CHL;  Weiming Wu, NCCHE
!***********************************************************************
    use EXP_Global_def,  only: ncn,nce,ncw,ncs,qx,qy
    use EXP_bndcond_def, only: QstringEXP
    use exp_telescoping, only: cellmap
    use prec_def,  only: ikind
    use size_def,  only: ncells
    use geo_def,   only: zb
    use flow_def,  only: u, v, h
    use bnd_def,   only: nhstr, nthstr, nmhstr, nmhvstr, mhv_str, mh_str, th_str, h_str
    use const_def, only: eps
    
    implicit none
    integer :: i,ii,j,iwse
    real(ikind) :: qut,qvt
    
    if(nHstr .gt. 0) then
      do i = 1,nHstr  !for each cell string
        do j=1,H_str(i)%NCells    !for each cell in string
          ii=H_str(i)%Cells(j) 
          ncn = cellmap(1,ii)
          nce = cellmap(3,ii)
          ncs = cellmap(5,ii)    
          ncw = cellmap(7,ii)                
          quT = u(ii)*h(ii)
          qvT = v(ii)*h(ii)
          if(quT .gt. 0.0 .and. nce .gt. Ncells) zb(ii) = zb(ncw)   
          if(quT .le. 0.0 .and. ncw .gt. Ncells) zb(ii) = zb(nce)   
          if(qvT .gt. 0.0 .and. ncn .gt. Ncells) zb(ii) = zb(ncs)  
          if(qvT .le. 0.0 .and. ncs .gt. Ncells) zb(ii) = zb(ncn)                                        
        enddo
      enddo ! end of each cell string
    endif  !nHstr strings

    if(nTHstr .gt. 0) then
      do iwse=1,nTHstr
        do j=1,TH_str(iwse)%NCells    !for each cell in string
          ii=TH_str(iwse)%Cells(j)  !added Mitch Brown - 10/20/2016 
          ncn = cellmap(1,ii)
          nce = cellmap(3,ii)
          ncs = cellmap(5,ii)    
          ncw = cellmap(7,ii)                
          quT = u(ii)*h(ii)
          qvT = v(ii)*h(ii)
          if(quT .gt. 0.0 .and. nce .gt. Ncells) zb(ii) = zb(ncw)   
          if(quT .le. 0.0 .and. ncw .gt. Ncells) zb(ii) = zb(nce)   
          if(qvT .gt. 0.0 .and. ncn .gt. Ncells) zb(ii) = zb(ncs)  
          if(qvT .le. 0.0 .and. ncs .gt. Ncells) zb(ii) = zb(ncn)           
        enddo 
      enddo !iwse-str      
    endif  !H_tide

    if(nMHstr .gt. 0) then
      do i = 1,nMHstr  !for each cell string                        
        do j=1,MH_str(i)%NCells   !for each cell in string
          ii=MH_str(i)%Cells(j) 
          ncn = cellmap(1,ii)
          nce = cellmap(3,ii)
          ncs = cellmap(5,ii)    
          ncw = cellmap(7,ii)                
          quT = u(ii)*h(ii)
          qvT = v(ii)*h(ii)
          if(quT .gt. 0.0 .and. nce .gt. Ncells) zb(ii) = zb(ncw)   
          if(quT .le. 0.0 .and. ncw .gt. Ncells) zb(ii) = zb(nce)   
          if(qvT .gt. 0.0 .and. ncn .gt. Ncells) zb(ii) = zb(ncs)  
          if(qvT .le. 0.0 .and. ncs .gt. Ncells) zb(ii) = zb(ncn)              
        enddo
      enddo ! end of each cell string
    endif  !H_multi

    if(nMHVstr .gt. 0) then
      do i = 1,nMHVstr  !for each cell string
        do j=1,MHV_str(i)%NCells   !for each cell in string
          ii=MHV_str(i)%Cells(j) 
          ncn = cellmap(1,ii)
          nce = cellmap(3,ii)
          ncs = cellmap(5,ii)    
          ncw = cellmap(7,ii)                
          quT = u(ii)*h(ii)
          qvT = v(ii)*h(ii)
          if(quT .gt. 0.0 .and. nce .gt. Ncells) zb(ii) = zb(ncw)   
          if(quT .le. 0.0 .and. ncw .gt. Ncells) zb(ii) = zb(nce)   
          if(qvT .gt. 0.0 .and. ncn .gt. Ncells) zb(ii) = zb(ncs)  
          if(qvT .le. 0.0 .and. ncs .gt. Ncells) zb(ii) = zb(ncn)             
        enddo
      enddo ! end of each cell string
    endif  !H_multi        
    
    return
    end subroutine Zbndzb_EXP_tel