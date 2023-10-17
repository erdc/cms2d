!===============================================================================
!
!     sed_imp     - Implicit solution of the sediment transport equations
!     coeffsourcesink_c - Calculates the sparse matrix coefficients for the 
!                         implicit solution scheme
!     bound_c     - Applies the sediment concentration boundary conditition
!                   for the implicit solution scheme
!
! written by Alex Sanchez, USACE-CHL
!===============================================================================
    
!***********************************************************************
    subroutine sed_imp
! Solves single and multi-sized sediment transport 
! written by Alex Sanchez, USACE-CHL;  Weiming Wu, NCCHE
!***********************************************************************
#include "CMS_cpp.h"
    use comp_lib
    use comvarbl
    use der_def, only: nder,nlim,gow,goa
    use der_lib, only: der_grad_eval
    use diag_def
    use diag_lib
    use flow_def
    use fric_def, only: bsxy
    use geo_def, only: idmap,zb,zbk,dzbx,dzby,ncface
    use interp_lib, only: interp_scal_cell2face
    use prec_def    
    use q3d_def
    use q3d_lib, only: q3d_eddyvert_mean
    use sed_def
    use size_def
    use solv_def, only: iconv
    implicit none
    integer :: i,k,ks
    real(ikind) :: val,chgpbk,chgCtk
    logical :: convCtk,convpbk
    real(ikind) :: tausx,tausy,ustars
#ifdef DIAG_MODE
    logical :: isnankind
#endif

!741 format(' Sediment: iter      Ctk_res       errpbk')
741 format(' Sediment: iter     Ctk_res       pbk_err')
751 format(5x,I9,1x,2(1pe13.4,1x))
349 format(' Sediment: iter      Ctk_res')
359 format(5x,I9,1x,1pe13.4,1x)
      
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
    
    !write(*,*)'bdj in sed_imp,ready to go to sedcapac, icapac=',icapac
    
    !=== Transport Capacity =====
    select case(icapac)  !CtstarP(i,ks),rs(i,ks)
    case(1); call sedcapac_lundcirp !Lund-CIRP
    case(2); call sedcapac_vanrijn  !Van Rijn
    case(3); call sedcapac_watanabe !Watanabe
    case(4); call sedcapac_soulsby  !Soulsby (1997)
    case(5); call wucapac           !Wu et al. (2000) (under testing)
    case(6); call sedcapac_c2shore   !bdj                                        
    end select
    
    !if(cohesivesed) call cohsedentrain !Cohesive sediment   !Wu 
    
#ifdef DIAG_MODE
    do i=1,ncells
      do ks=1,nsed
        if(isnankind(CtstarP(i,ks)))then
          call diag_print_var(i,4)
          call diag_print_error('Problem calculating transport capacity')            
        endif    
      enddo    
    enddo
#endif
    
    !=== Vertical diffusivity =======
    !Note: The vertical diffusivity is used in the total load correction factor calculation.
    if(icapac/=1)then
      do i=1,ncells  
        call surface_windwave_stress(i,tausx,tausy,ustars)
        epsvk(i,:) = q3d_eddyvert_mean(h(i),bsxy(i),ustars)/schmidt !Assume Schmidt number the same for all size classes
      enddo    
    endif
    
    !=== Horizontal Mixing Coefficient ====
!$OMP PARALLEL DO PRIVATE(i,k)
    do i=1,ncells
      do k=1,ncface(i)
        visk(k,i)=viskfl(k,i)/schmidt+cmixsed
      enddo  
    enddo
!$OMP END PARALLEL DO
    
    !=== Transport Corrections for Bedslope ===
    if(ibedslope==2) call bedslopecor_bailard
    
    !=== Adjust Concentration Capacity =====
    !Limit values
    !!$OMP PARALLEL DO PRIVATE(i,ks)
    !do ks=1,nsed
    !  do i=1,ncells
    !    CtstarP(i,ks) = min(CtstarP(i,ks),200.0) !Limit concentration capacity ***************
    !  enddo
    !enddo  
    !!$OMP END PARALLEL DO    
    !where(CtstarP>200.0) CtstarP=200.0
    
    !Hardbottom
    !do ih=1,nhard
    !  i=idhard(ih)
    !  fac=max(0.0,min(1.0,0.01*(zb1(i)-hardzb(i))/d50(i)))
    !  CtstarP(i,:)=CtstarP(i,:)*fac
      !!do ks=1,nsed
      !!  Ctkstarhard = Ctk(i,ks) + rhosed*solid*pbk(i,ks,1)*(zb(i)-hardzb(i)) &
      !!      /(alphat(i)*wsfall(ks)*max(scalemorph,1.0)*dtime) &
      !!      + Sb(i,ks)/(alphat(i)*wsfall(ks)+1.0e-6)
      !!  Ctkstarhard = max(Ctkstarhard,0.0)
      !!  CtstarP(i,ks) = min(CtstarP(i,ks),Ctkstarhard/max(pbk(i,ks,1),1.0e-6))
      !!enddo !ks
    !enddo !ih
    
    !=== Inflow Boundary Conditions ====
    call sedbnd_eval  !CtstarP(nck,:)
    
    !=== Determine whether to update morphology (bed) and bed composition ====
    if(timehrs<tStartMorph .or. scalemorph<1.0e-6)then
      calcmorph=.false.
    else
      calcmorph=.true.  
    endif    
    if(.not.calcbedcomp .or. timehrs<tStartBedComp)then
      calcbedcomp=.false.
    else
      calcbedcomp=.true.
    endif        
    
    !=== Adaptation coefficient ==============
    call adaptcoef !alphat(i)           
    
    !if(cohesivesed) call fallvel_cohsed   !By Wu
    !
    !if(cohesivesed) then   !Wu
    !   do i=1,ncells
    !      EtstarP(i,1)=min(EtstarP(i,1),wsfallcohsed(i)*200.0)
    !      CtstarP(i,1)=EtstarP(i,1)/wsfallcohsed(i)
    !   enddo
    !   if(nsed.gt.1) then
    !      do i=1,ncells
    !         EtstarP(i,2:nsed)=alphat(i)*wsfall(2:nsed)*CtstarP(i,2:nsed)
    !         EtstarP(i,2:nsed)=EtstarP(i,2:nsed)+(EtstarP(i,1)-EtstarP(i,2:nsed))*  & !affected by cohesive sed.
    !                              (max(0.0,min(1.0,(pbk1(i,1)-pcmin)/(pcmax-pcmin))))**2
    !      enddo
    !   endif
    !else
    !   do i=1,ncells
    !      EtstarP(i,:)=alphat(i)*wsfall*CtstarP(i,:)
    !   enddo
    !endif
    
    !=== Total load correction factor ========
    select case(ibt) !btk(i,ks)
    case(1); call btklogexp  
    case(2); call btklogrouse 
    end select
    
    !=== Mixing Layer Thickness ===========
    if(nsed>1) call mixing_layer !db(i,1)
      
    !=== Print header to screen ============
    if(nsed>1  .and. sedcouple)then 
      write(*,741)
    else
      write(*,349)
    endif          
    
    !=== 3D dispersion terms =======================
!#ifdef DEV_MODE
    !if(q3d) call q3d_sed !s3dxk, s3dyk             
!#endif
    
    !=== Solve each A-D equation =====================
    errCtk0 = 0.0_ikind
    errpbk0 = 0.0_ikind
    do itersed=1,maxitersed 
       !if(cohesivesed) call fallvel_cohsed   !By Wu
       errCtk = 0.0_ikind
       errpbk = 0.0_ikind
       rsCtkmax = 0.0
       do ks=1,nsed
         !=== Coefficients and solve for Ctk =======
         select case(ndsch) 
         case(2); call coeffsourcesink_c(hybridcoef,ks)
         case(3); call coeffsourcesink_c(powerlawcoef,ks)
         case(4); call coeffsourcesink_c(exponentialcoef,ks)
         case default; call coeffsourcesink_c(upwindcoef,ks)
         end select
         if(ncellsimple==ncells)then !No gradients required
           select case(ndsch) !Anti-diffusion corrections
           case(5); call defcorhlpa(Ctk(:,ks),su)
           case(6); call defcorgamma(gammadefcor,Ctk(:,ks),su)
           case(7); call defcorgamma(cubistadefcor,Ctk(:,ks),su)
           case(8); call defcorgamma(alvsmartdefcor,Ctk(:,ks),su)
           case(9); call defcorgamma(hoabdefcor,Ctk(:,ks),su)
           end select
         else
           select case(ndsch) !Anti-diffusion corrections
           case(5); call defcorhlpagrad(Ctk(:,ks),dCtkx(:,ks),dCtky(:,ks),su)
           case(6); call defcorgammagrad(gammadefcor,Ctk(:,ks),dCtkx(:,ks),dCtky(:,ks),su)
           case(7); call defcorgammagrad(cubistadefcor,Ctk(:,ks),dCtkx(:,ks),dCtky(:,ks),su)
           case(8); call defcorgammagrad(alvsmartdefcor,Ctk(:,ks),dCtkx(:,ks),dCtky(:,ks),su)
           case(9); call defcorgammagrad(hoabdefcor,Ctk(:,ks),dCtkx(:,ks),dCtky(:,ks),su)
           case default; if(skewcor) call defcorparagrad(dCtkx(:,ks),dCtky(:,ks),su) !Skewness correction      
           end select
         endif
         call bound_c(ks)
         if(debug_mode) call check_variables(4)
         call solve(acoef,su,sp,rsCtk,Ctk(:,ks),4)
         if(iconv==0)then 
           call write_debug
           stop
         endif    
         errCtk = max(rmom(4),errCtk)
         do i=1,ncells
           rsCtkmax(i) = min(rsCtkmax(i),rsCtk(i)) !calc max residuals
           Ctk(i,ks) = min(max(Ctk(i,ks),0.0),Cteqmax) !Bound values
         enddo
         !=== Compute gradients
         call der_grad_eval(gow,nlim,Ctk(:,ks),dCtkx(:,ks),dCtky(:,ks))
       enddo !ks
       
       !=== Bedslope term ====================
       if(do_bedslope) call bedslope !Sb(i,ks)
         
       !=== Wave-induced sed transport ============
       if(wavesedtrans) call sed_wave !Sb(i,ks)+Sw(i) 
        
       !=== Bed Change and Sorting ==============
       !if(calcbedcomp)then
         if(singlesize)then !no sorting needed
           call bedchange !dzb(i)
         else
           !if(consolidation)then
           !  call bedchangesort_consolid
           !else  
             call bedchangesort !db(i,1),dzb(i),dzbk(i,ks),pbk(i,ks,1)
           !endif
         endif
         call struct_dzb !dzb(i) and optionally dzbk(i,ks)
         call check_hardbottom
       !endif                  
       
       !=== Print to screen ==============
       if(mod(itersed,10)==0 .or. itersed==1)then
         if(nsed>1 .and. sedcouple)then
           write(*,751) itersed,errCtk,errpbk
         else
           write(*,359) itersed,errCtk
         endif
       endif
         
       !=== Check for convergence ============   
       if(itersed>=maxitersed/2 .or. errCtk<1.0e-13)then
       !if(itersed>=maxitersed/2)then    
        if(singlesize)then
         if(errCtk<=tolCtk) exit   !Error, concentrations
         chgCtk=abs(errCtk-errCtk0)  
         if(chgCtk<=1.0e-8) exit   !Absolute error change, concentrations
         chgCtk=chgCtk/max(errCtk0,tolCtk)
         if(chgCtk<=0.01) exit     !Relative error change, concentrations
         errCtk0=errCtk
        else   
         convCtk=(errCtk<=tolCtk)  !Error, concentrations
         convpbk=(errpbk<=tolpbk)  !Error, bed composition
         if(convCtk .and. convpbk) exit   
         !chgCtk=abs(errCtk-errCtk0)  !Absolute error change, concentrations
         !chgpbk=abs(errpbk-errpbk0)  !Absolute error change, bed composition
         !convCtk=(convCtk .or. chgCtk<=1.0e-7)
         !convpbk=(convpbk .or. chgpbk<=1.0e-6)
         if(convCtk .and. convpbk) exit
         !chgCtk=chgCtk/max(errCtk0,tolCtk) !Relative error change, concentrations
         !chgpbk=chgpbk/max(errpbk0,tolpbk) !Relative error change, bed composition
         !convCtk=(convCtk .or. chgCtk<=0.01)
         !convpbk=(convpbk .or. chgpbk<=0.01)
         !if(convCtk .and. convpbk) exit
         errCtk0=errCtk
         errpbk0=errpbk
        endif
       endif
    enddo !iteration loop
    rmom(4) = errCtk  
    
    !=== Bed sorting and gradation =============
    if(.not.singlesize .and. calcbedcomp)then
      call bedgrad !pbk(i,ks,2:nlay)
      !if(.not.consolidation) call bedgrad !pbk(i,ks,2:nlay) !Wu
    endif  
      
    !=== Write final sediment error estimates ====   
    open(dgunit,file=dgfile,access='append')                 
    if(nsed>1 .and. sedcouple)then                  
      write(dgunit,741)
      write(dgunit,751) itersed,errCtk,errpbk  
      if(mod(itersed,10)/=0)then        
        if(itersed<maxitersed-1)then
          write(*,751) itersed,errCtk,errpbk  
        endif
      endif
    else
      write(dgunit,349)
      write(dgunit,359) itersed,errCtk
      if(mod(itersed,10)/=0)then        
        if(itersed<maxitersed-1)then
          write(*,359) itersed,errCtk
        endif 
      endif
    endif  
    close(dgunit)  
    
    !=== Percentiles ================
    if(.not.singlesize)then
      call sedpercentile(50,d50)
      call sedpercentile(90,d90)       
    endif
 
    !==== Update bed elevation ====================
    !---Bed Change --------------------------
    if(calcmorph)then
      do i=1,ncells
        if(iwet(i)==1)then
          if(abs(dzb(i))>0.5*h(i) .or. abs(dzb(i))>dzbmax)then
            call diag_print_warning('Large bed change')
            call print_sedvar(i,1)
          endif
        endif
        zb(i) = zb1(i) + dzb(i)   !*****************
      enddo  
          
      !=== Avalanching ============================    
      if(do_aval) call avalanche
      
#ifdef DEV_MODE
      !=== Erosion of dry cells ===================
      if(erosdry%calc) call sed_dry
#endif

      !=== Bed-slopes =============================
      call der_grad_eval(goa,0,zb,dzbx,dzby) !Bed-slope
    
      !=== Update bed elevation at cell faces =====
      call interp_scal_cell2face(zb,0,zbk,dzbx,dzby)
    
      !=== Correct concentrations for depth changes =======
      call sed_concdepthchange
      
    else !No Bed Change
!$OMP PARALLEL
!$OMP DO PRIVATE(i)
      do i=1,ncells        
        dzb(i) = 0.0
      enddo
!$OMP END DO
      if(nsed>1)then
!$OMP DO PRIVATE(i)
        do i=1,ncells          
          dzbk(i,:) = 0.0 !Fractional bed change (for each size class)
        enddo
!$OMP END DO
      endif
!$OMP END PARALLEL
    endif
    
    !=== Update sediment variables ====
    call sed_total  !cmb(i),Ct(i),Ctstar(i),rs(i),qtx(i),qty(i) 
    
    !== Time step statistics ============
    call sed_step_stat
    
    !=== Total Sediment Budget =================                                 
    if(sedbalance .and. calcmorph) call sed_balance
      
    !=== Apply boundary condition to zb ====
    call bndzb
        
    return
    end subroutine sed_imp

!****************************************************************************
    subroutine coeffsourcesink_c(schmcoef,ks)
! Calculates the sparse matrix coefficients for the implicit solution scheme
! written by Weiming Wu, NCCHE, Oct. 2008
! modified by Alex Sanchez, USACE-CHL, Oct 2010
!****************************************************************************
#include "CMS_cpp.h"
    use size_def, only: ncells
    use geo_def,  only: areap,zb,ncface,cell2cell,dsxy
    use flow_def, only: iwet,su,h1,sp,h,h2,visk,hk,acoef,flux
    use sed_def,  only: ctk1,btk1,btk,ctk2,btk2,hardzb,dzbmax,ctk,rhosed,solid,pbk,alphat,wsfall
    use sed_def,  only: scalemorph,sb,ctstarp,ctkstar,rsk,cmixbedload,diam
    use fric_def, only: bsvel
!!#ifdef DEV_MODE
!!    use q3d_def
!!#endif
    use interp_def, only: fintp
!    use comp_lib
    use comvarbl, only: ntsch,dtime,ctsch1,ctsch2,ctsch
    use diag_lib, only: diag_print_var, diag_print_error
    use prec_def, only: ikind
    
    implicit none
    integer :: i,ks,k,nck,ih,idrym
    real(ikind) :: ddk,a0,awsvolp,fkip,rskface,dzblim
    real(ikind) :: Ctkstarhard 
    real(8)     :: TempVar             !HLI 01/13/2017
    logical :: isnankind

    interface
      function schmcoef(dk,fk) result(ak)
        use prec_def
        implicit none
        real(ikind),intent(in) :: dk,fk
        real(ikind) :: ak
      end function
    endinterface  

!.....Temporal terms
    if(ntsch==1)then
!$OMP PARALLEL DO PRIVATE(i,a0)        
      do i=1,ncells          
        a0=iwet(i)*areap(i)/dtime      
        su(i)=a0*h1(i)*Ctk1(i,ks)/btk1(i,ks)
        sp(i)=-a0*(h(i)*(1.0/btk(i,ks)-1.0)+h1(i))
      enddo    
!$OMP END PARALLEL DO        
    else
!$OMP PARALLEL DO PRIVATE(i,a0)    
      do i=1,ncells 
        a0=iwet(i)*areap(i)/dtime  
        su(i)=a0*(ctsch1*h1(i)*Ctk1(i,ks)/btk1(i,ks)-ctsch2*h2(i)*Ctk2(i,ks)/btk2(i,ks))
        sp(i)=-a0*(ctsch*h(i)*(1.0/btk(i,ks)-1.0)+ctsch1*h1(i)-ctsch2*h2(i))
      enddo
!$OMP END PARALLEL DO
    endif

    !!Hardbottom
    !do ih=1,nhard
    !  i=idhard(ih)
    !  Ctkstarhard = Ctk(i,ks) + rhosed*solid*pbk(i,ks,1)*(zb(i)-hardzb(i)) &
    !      /(alphat(i)*wsfall(ks)*max(scalemorph,1.0)*dtime) &
    !      + Sb(i,ks)/(alphat(i)*wsfall(ks)+1.0e-6)
    !  Ctkstarhard = max(Ctkstarhard,0.0)
    !  CtstarP(i,ks) = min(CtstarP(i,ks),Ctkstarhard/max(pbk(i,ks,1),1.0e-6))
    !  !Ctkstar(i,ks)=CtstarP(i,ks)*pbk(i,ks,1) 
    !  !Ctkstar(i,ks) = min(Ctkstar(i,ks),Ctkstarhard)
    !  !CtstarP(i,ks) = Ctkstar(i,ks)/max(pbk(i,ks,1),1.0e-6)
    !enddo !ih    
    
!$OMP PARALLEL DO PRIVATE(i,k,idrym,nck,ddk,awsvolp,fkip,rskface,Ctkstarhard,dzblim,Tempvar)     !MEB 01/13/2017
    do i=1,ncells
      
      !--- Limit capacity to avoid excessive erosion or eroding past the hard bottom ---
      dzblim = min(zb(i)-hardzb(i),dzbmax)
      TempVar = alphat(i)*wsfall(ks) + 1.0e-20                                         !MEB 01/13/2017
      !print*, TempVar, dtime 
      Ctkstarhard = Ctk(i,ks) + rhosed*solid*pbk(i,ks,1)*dzblim /(TempVar*max(scalemorph,1.0)*dtime) + Sb(i,ks)/TempVar   !changed scalemorph to computed val - meb 03/11/2019
      
      
      Ctkstarhard = max(Ctkstarhard,0.0)
      CtstarP(i,ks) = min(CtstarP(i,ks),Ctkstarhard/max(pbk(i,ks,1),1.0e-6))
      
      !--- Fractional Equilibrium Concentration ---------
      Ctkstar(i,ks)=CtstarP(i,ks)*pbk(i,ks,1)
      
      !--- Erosion and Deposition terms -----------------------------------------
      awsvolp=areap(i)*alphat(i)*wsfall(ks) !Note iwet is already included in alphat
      su(i)=su(i)+awsvolp*Ctkstar(i,ks)
      sp(i)=sp(i)-awsvolp
      !su(i)=su(i)+areap(i)*iwet(i)*Etkstar(i,ks)
      !if(cohesivesed.and.ks.eq.1) then       
      !  sp(i)=sp(i)-areap(i)*alphat(i)*wsfallcohsed(i)  !Note idry is already included in alphat
      !else
      !  sp(i)=sp(i)-areap(i)*alphat(i)*wsfall(ks) !Note idry is already included in alphat  
      !endif

      !--- Matrix Coefficients ----------------------------------
      do k=1,ncface(i)    
        fkip=fintp(k,i)
        rskface=fkip*rsk(cell2cell(k,i),ks)+(1.0-fkip)*rsk(i,ks)
        !ddk=rskface*visk(k,i)*hk(k,i)*dsxy(k,i)
        ddk=(rskface*visk(k,i)+(1.0-rskface)*cmixbedload*bsvel(i)*diam(ks))*hk(k,i)*dsxy(k,i) !Note: visk already includes Schmidt number and base value
        !ddk=rskface*dk(i,k)/schmidt+dsxy(k,i)*hk(k,i)*cmixsed
        acoef(k,i)=schmcoef(ddk,flux(k,i))  
!        acoef(k,i)=0.0 !bdj  !ASK BRAD                           
        if(isnankind(acoef(k,i)) .or. abs(acoef(k,i))>1.0e20)then
          call diag_print_var(i,4)
          call diag_print_error('Problem calculating transport capacity')     
        endif
      enddo      
      ! if(i.ge.1514.and.i.le.1524) then
      !    write(*,*)'bdj i, h(i), u(i), acoef(1,i), Ct(i)',i, h(i), u(i), acoef(1,i), Ct(i)
      ! endif

      !--- Wetting and drying -------------------------
      if(iwet(i)==1)then
        do k=1,ncface(i)
          nck=cell2cell(k,i)
           if(iwet(nck)==0)then    !Side is dry                 
            Ctk(nck,ks)=Ctk(i,ks)
            acoef(k,i)=0.0
          endif
        enddo
      elseif(iwet(i)==0)then
!        Ctk(i,ks)=0.0        
!        acoef(:,i)=0.0
!        sp(i)=-1.0
!        su(i)=0.0
        idrym=0
        do k=1,ncface(i)
          nck=cell2cell(k,i)
          if(nck<=ncells) idrym=idrym+iwet(nck)
        enddo             
        do k=1,ncface(i)
          nck=cell2cell(k,i)
          if(idrym>=1)then
            acoef(k,i)=iwet(nck)
          else
            acoef(k,i)=dsxy(k,i)
          endif
            if(nck>ncells) acoef(k,i)=0.0
        enddo
        su(i)=0.0
        sp(i)=0.0
      endif         
    enddo
!$OMP END PARALLEL DO

!#ifdef DEV_MODE
!    !3D Dispersion terms
!    if(q3d)then
!!$OMP PARALLEL DO PRIVATE(i,dvarx,dvary)
!      do i=1,ncells
!        call dx2d(i,s3dxk,dvarx)         
!        call dy2d(i,s3dyk,dvary) 
!        su(i)=su(i)-(dvarx+dvary)*areap(i)
!      enddo
!!$OMP END PARALLEL DO
!    endif
!#endif
    
    return
    end subroutine coeffsourcesink_c

!***********************************************************************
    subroutine bound_c(ks)
! Applies the sediment concentration boundary conditition
! for the implicit solution scheme.
! Blocks off dry regions which are not solved
! written by Alex Sanchez, USACE-CHL;  Weiming Wu, NCCHE
! Strucures added by Weiming Wu, NCCHE
!***********************************************************************
    use size_def
    use geo_def
    use flow_def, only: flux,acoef,su,sp,iwet
    use bnd_def
    use sed_def
    use const_def, only: small
    implicit none
    integer :: i,j,k,ks,nck,ibnd     

!--- All forcing boundaries -----------------------------
    do ibnd=1,nbndstr
      do j=1,bnd_str(ibnd)%ncells
        i=bnd_str(ibnd)%cells(j)
        k=bnd_str(ibnd)%faces(j)
        if(iwet(i)==0) cycle
        nck=cell2cell(k,i)
        if(flux(k,i)<0.0)then !Inflow
           Ctk(nck,ks)=pbk(i,ks,1)*CtstarP(nck,ks)
           !Ctk(nck,ks)=Ctk(i,ks)*CtstarP(nck,ks)/CtstarP(i,ks) !Better but needs testing
           su(i)=su(i)+acoef(k,i)*Ctk(nck,ks) !Flux bc
           sp(i)=sp(i)-acoef(k,i)
        else                  !Outflow         
           Ctk(nck,ks)=Ctk(i,ks)  !Copy to dummy cell       
           pbk(nck,ks,1)=pbk(i,ks,1)                             
        endif 
        acoef(k,i)=0.0      
      enddo !j
    enddo !ibnd

    call struct_Ctk(ks) !Structures

    return
    end subroutine bound_c    