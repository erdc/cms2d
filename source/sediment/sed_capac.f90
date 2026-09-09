!========================================================================
! CMS Sediment Transport Capacity/Formulas
!
! Contains the following:
!     sedcapac_soulsby  - Calculates the concentration capacity using the Soulsby-van Rijn transport equations
!     sedcapac_vanrijn  - Calculates the concentration capacity using the Van Rijn transport equations
!     sedcapac_watanabe - Calculates the concentration capacity using the Watanabe transport equation
!     sedcapac_lundcirp - Calculates the transport capacity based on the Lund-CIRP equations      
!     sedcapac_c2shore  - Calculates the transport capacity based on the C2SHORE model                           bdj 2019
!
! written by Alex Sanchez, USACE-CHL 
! written by Brad Johnson, USACE-CHL
!========================================================================

!********************************************************************
subroutine sedcapac_lundcirp
  ! Calculates the transport capacity based on the Lund-CIRP equations     
  ! written by Alex Sanchez, USACE-CHL
  !********************************************************************    
  use size_def
  use geo_def, only: mapid,dzbx,dzby
  use flow_def
  use const_def
  use wave_flowgrid_def
  use sed_def
  use sed_lib, only: critslpcor_dey
  use rol_def, only: roller
  use cms_def
  use prec_def
  implicit none
  integer :: i,ks,iripple
  real(ikind) :: tauct,tauwt,tauwmt,taucwt,taucwmt,tauctb,taucwtb,taucwmtb
  real(ikind) :: fcf,fwf,fcwf,BDpart,Ustc,Ustw,Hrms,phi,alfa,ur,gamma,um
  real(ikind) :: Qss,Qbs,Qsm,Qsr,Qbm,Qba,Qts,Uw,T,dbrk,fac, Hw             !Hw, added by Wu
  !logical :: isnankind

  iripple = 1               
  if(noptset>=3)then  !Waves
    gamma = 0.78      
    !rolfac = 2.0/rhow/sqrt(grav)
    do i=1,ncells     
      if(iwet(i)==0)then
        CtstarP(i,:)=0.0
        rsk(i,:)=1.0
        cycle 
      endif
      phi = abs(wang(i)-atan2(v(i),u(i)))     !Current-wave angle                  
      call shearlund(iripple,h(i),uv(i),Worb(i),Wper(i),phi,rhosed,rhow,D50(i),&
           tauct,tauwt,tauwmt,taucwt,taucwmt,tauctb,taucwtb,taucwmtb,fcf,fwf,fcwf)              
      if(wavesedtrans)then
        Hrms = Whgt(i)/sqrttwo  
        call brkwavratio(h(i),Hrms,Wlen(i),gamma,alfa)      
        fac=0.5+0.5*cos(pi*min(max(h(i)/max(Whgt2(i),0.01)-1.0,0.0)/2.0,1.0))
      endif
      do ks=1,nsed
        call susplund(h(i),uv(i),Whgt(i),Worb(i),rhosed,rhow,diam(ks),   &      !Whgt(i), added by Wu
             wsfall(ks),dstar(ks),tauct,tauwt,tauwmt,taucwt,taucwmt,&
             fcf,fwf,wavediss(i),taucr(ks),cak(i,ks),epsvk(i,ks),Qss,BDpart,Ustc,Ustw)                 
        call bedlund(rhosed,rhow,tauctb,taucwtb,taucwmtb,taucr(ks),Qbs)
        Qbs = scalebed*Qbs/varsigma(i,ks)       !Bed Load Capacity        
        Qss = scalesus*Qss/varsigma(i,ks)       !Suspended Load Capacity 
        Qts = Qbs + Qss                         !Total Load Capacity
        CtstarP(i,ks) = rhosed*Qts/(uv(i)*h(i)+small) !Total-load concentration capacity of total load, kg/m^3              
        CtstarP(i,ks) = min(CtstarP(i,ks),Cteqmax)
        rsk(i,ks) = max(Qss,small)/max(Qts,small) !Fraction of suspended sediment      

        !Wave induced sediment transport
        if(wavesedtrans .and. Whgt(i).gt.2*hmin .and. h(i).gt.3*hmin)then    
          call crossmean(h(i),Whgt(i),Wper(i),Wlen(i),alfa,&
               rhosed,rhow,diam(ks),wsfall(ks),            &
               fcf,taucr(ks),taucwtb,taucwmtb,cak(i,ks),epsvk(i,ks),      &
               um,ur,Qsm,Qsr,Qbm)
          call asymmetry(h(i),Hrms,Wper(i),Wlen(i),Worbrep(i),&
               rhosed,rhow,diam(ks),fwf,taucr(ks),taucwtb,taucwmtb,Qba)
          QwsP(i,ks) = rhosed*(fba*Qba+fac*(fsr*Qsr-fsm*Qsm-fbm*Qbm)) !Potential net onshore transport, kg/m/sec
        endif
      enddo
    enddo
  else         !No waves
    phi = 0.0; Uw = 0.0; T = 10.0; dbrk = 0.0
    Hw = 0.0     !added by Wu, July 2025
    do i=1,ncells  
      if(iwet(i)==0)then
        CtstarP(i,:)=0.0
        rsk(i,:)=1.0
        cycle 
      endif
      call shearlund(iripple,h(i),uv(i),Uw,T,phi,rhosed,rhow,D50(i),&
           tauct,tauwt,tauwmt,taucwt,taucwmt,tauctb,taucwtb,taucwmtb,fcf,fwf,fcwf)
      do ks=1,nsed
        call susplund(h(i),uv(i),Hw,Uw,rhosed,rhow,diam(ks),       &         !Hw, added by Wu
             wsfall(ks),dstar(ks),tauct,tauwt,tauwmt,taucwt,taucwmt,&
             fcf,fwf,dbrk,taucr(ks),cak(i,ks),epsvk(i,ks),Qss,BDpart,Ustc,Ustw)                
        call bedlund(rhosed,rhow,tauctb,taucwtb,taucwmtb,taucr(ks),Qbs)                            
        Qbs = scalebed*Qbs/varsigma(i,ks)       !Bed Load Capacity        
        Qss = scalesus*Qss/varsigma(i,ks)       !Suspended Load Capacity 
        Qts = Qbs + Qss          !Total Load Capacity     
        CtstarP(i,ks) = rhosed*Qts/max(uv(i)*h(i),small) !Total-load concentration capacity of total load, kg/m^3  
        CtstarP(i,ks) = min(CtstarP(i,ks),Cteqmax)
        rsk(i,ks) = max(Qss,small)/max(Qts,small)
      enddo
    enddo
  endif

  return
end subroutine sedcapac_lundcirp

!****************************************************************************
subroutine sedcapac_soulsby
  ! Calculates the concentration capacity using the 
  ! Soulsby-van Rijn transport equations
  !
  ! written by Alex Sanchez, USACE-CHL
  !*****************************************************************************          
  use size_def
  use geo_def, only: dzbx,dzby
  use flow_def
  use comvarbl
  use wave_flowgrid_def
  use sed_def
  use sed_lib, only: sedtrans_wavcur_soulsby,sedtrans_cur_soulsby
  use fric_def, only: cfrict 
  use cms_def
  use const_def, only: sqrttwo,small
  use prec_def
  implicit none
  integer :: i,ks
  real(ikind) :: Qbs,Qss,Qts,Cd,Urms

  if(noptset>=3)then !Waves and currents
     !$OMP PARALLEL DO PRIVATE(i,ks,Qbs,Qss,Qts,Cd,Urms)
     do i=1,ncells
        if(iwet(i)==0)then
           CtstarP(i,:)=0.0
           rsk(i,:)=1.0
           cycle 
        endif
        Cd = (0.4/(log(h(i)/min(0.0006,h(i)))-1.0))**2 !z0=0.0006 [m]
        !Cd = cfrict(i)
        Urms = Worbrep(i)/sqrttwo
        do ks=1,nsed
           !call sedtrans_wavcur_soulsby(h(i),uv(i),Cd,Urms,&
           !    diam(ks),diam(ks),d90(i),varsigma(i,ks),Qbs,Qss)
           call sedtrans_wavcur_soulsby(grav,h(i),uv(i),Cd,Urms,&
                specgrav,diam(ks),diam(ks),dstar(ks),varsigma(i,ks),Qbs,Qss) !Notes: d50=d90=dk,Uw=Worbrep
           Qbs = scalebed*Qbs    !Bed Load Capacity        
           Qss = scalesus*Qss    !Suspended Load Capacity 
           Qts = Qbs + Qss       !Total Load Capacity          
           CtstarP(i,ks) = rhosed*Qts/max(uv(i)*h(i),1.0e-5) !Total-load concentration capacity of total load, kg/m^3              
           CtstarP(i,ks) = min(CtstarP(i,ks),Cteqmax)
           rsk(i,ks) = Qss/(Qts+small) !Fraction of suspended sediment
        enddo !ks
     enddo !i
     !$OMP END PARALLEL DO
  else   !Currents only, no waves
     !$OMP PARALLEL DO PRIVATE(i,ks,Qbs,Qss,Qts)
     do i=1,ncells   
        if(iwet(i)==0)then
           CtstarP(i,:)=0.0
           rsk(i,:)=1.0
           cycle 
        endif
        do ks=1,nsed
           !call sedtrans_cur_soulsby(h(i),uv(i),&
           !  diam(ks),d90(i),dstar(ks),varsigma(i,ks),Qbs,Qss)
           call sedtrans_cur_soulsby(grav,h(i),uv(i),&
                specgrav,diam(ks),diam(ks),dstar(ks),varsigma(i,ks),Qbs,Qss) !d50=d90=dk
           Qbs = scalebed*Qbs    !Bed Load Capacity        
           Qss = scalesus*Qss    !Suspended Load Capacity 
           Qts = Qbs + Qss       !Total Load Capacity          
           CtstarP(i,ks) = rhosed*Qts/max(uv(i)*h(i),1.0e-5) !Total-load concentration capacity of total load, kg/m^3              
           CtstarP(i,ks) = min(CtstarP(i,ks),Cteqmax)
           rsk(i,ks) = Qss/(Qts+small) !Fraction of suspended sediment
        enddo !ks
     enddo !i
     !$OMP END PARALLEL DO
  endif

  return
end subroutine sedcapac_soulsby

!*****************************************************************************  
subroutine sedcapac_vanrijn
  ! Calculates the concentration capacity using the Van Rijn transport equations
  !
  ! written by Alex Sanchez, USACE-CHL
  !******************************************************************************      
  use size_def
  use flow_def
  use comvarbl
  use wave_flowgrid_def
  use sed_def
  use sed_lib, only: sedtrans_wavcur_vanrijn,sedtrans_cur_vanrijn
  use cms_def
  use const_def, only: small
  use prec_def
  implicit none
  integer :: i,ks
  real(ikind) :: Qbs,Qss,Qts

  if(noptset>=3)then
     !$OMP PARALLEL DO PRIVATE(i,ks,Qbs,Qss,Qts)        
     do i=1,ncells
        if(iwet(i)==0)then
           CtstarP(i,:)=0.0
           rsk(i,:)=1.0
           cycle 
        endif
        do ks=1,nsed
           !call sedtrans_wavcur_vanrijn(grav,h(i),uv(i),Worb(i),Wper(i),&
           !   specgrav,diam(ks),d90(i),dstar(ks),varsigma(i,ks),Qbs,Qss)
           call sedtrans_wavcur_vanrijn(grav,h(i),uv(i),Worb(i),Wper(i),&
                specgrav,diam(ks),diam(ks),dstar(ks),varsigma(i,ks),Qbs,Qss) !d50=d90=dk
           Qbs = scalebed*Qbs    !Bed Load Capacity        
           Qss = scalesus*Qss    !Suspended Load Capacity 
           Qts = Qbs + Qss       !Total Load Capacity    
           CtstarP(i,ks) = rhosed*Qts/max(uv(i)*h(i),small) !Total-load concentration capacity of total load, kg/m^3              
           CtstarP(i,ks) = min(CtstarP(i,ks),Cteqmax)
           rsk(i,ks) = Qss/(Qts+small) !Fraction of suspended sediment
        enddo !ks
     enddo !i
     !$OMP END PARALLEL DO      
  else     !No waves
     !$OMP PARALLEL DO PRIVATE(i,ks,Qbs,Qss,Qts)
     do i=1,ncells
        if(iwet(i)==0)then
           CtstarP(i,:)=0.0
           rsk(i,:)=1.0
           cycle 
        endif
        do ks=1,nsed
           !call sedtrans_cur_vanrijn(grav,h(i),uv(i),&
           !  diam(ks),d90(i),specgrav,dstar(ks),varsigma(i,ks),Qbs,Qss)
           call sedtrans_cur_vanrijn(grav,h(i),uv(i),specgrav,&
                diam(ks),diam(ks),dstar(ks),varsigma(i,ks),Qbs,Qss) !d50=d90=dk
           Qbs = scalebed*Qbs    !Bed Load Capacity        
           Qss = scalesus*Qss    !Suspended Load Capacity 
           Qts = Qbs + Qss       !Total Load Capacity
           CtstarP(i,ks) = rhosed*Qts/max(uv(i)*h(i),small) !Total-load concentration capacity of total load, kg/m^3              
           CtstarP(i,ks) = min(CtstarP(i,ks),Cteqmax)
           rsk(i,ks) = Qss/(Qts+small) !Fraction of suspended sediment          
        enddo !ks
     enddo !i
     !$OMP END PARALLEL DO
  endif

  return
end subroutine sedcapac_vanrijn

!*****************************************************************************  
subroutine sedcapac_watanabe
  ! Calculates the concentration capacity using the Watanabe transport equation
  !
  ! written by Alex Sanchez, USACE-CHL
  !*****************************************************************************          
  use size_def
  use flow_def
  use comvarbl
  use wave_flowgrid_def
  use fric_def, only: bsxy
  use sed_def
  use sed_lib, only: shearwatanabecw,sedtrans_wavcur_watanabe,&
       sedtrans_wavcur_vanrijn,sedtrans_cur_vanrijn
  use cms_def
  use const_def, only: small
  use prec_def
  implicit none
  integer :: i,ks
  real(ikind) :: phi,Qbs,Qss,Qts,taumax
  !logical :: isnankind

  if(noptset>=3)then !Waves
     !$OMP PARALLEL DO PRIVATE(i,ks,phi,Qbs,Qss,Qts,taumax)
     do i=1,ncells  
        if(iwet(i)*uv(i)<=1.0e-4)then
           CtstarP(i,:)=0.0
           rsk(i,:)=1.0
           cycle 
        endif
        phi = abs(Wang(i)-atan2(v(i),u(i)))   !Current-wave angle      
        !!        Aw = Worb(i)*Wper(i)/twopi
        !!        Rfac = max(Aw/(2.5*D50(i)),1.0) !Relative roughness
        !!        fw = exp(5.5*Rfac**(-0.2)-6.3)  !Wave friction factor
        !!        tauw = 0.5*rhow*fw*Worb(i)*Worb(i)
        !!        taumax = sqrt((bsxy(i)+tauw*cos(phi))**2+(tauw*sin(phi))**2) !Soulsby 1997          
        call shearwatanabecw(rhow,d50(i),bsxy(i),phi,Worb(i),Wper(i),taumax)
        do ks=1,nsed                   
           call sedtrans_wavcur_watanabe(uv(i),taumax,taucr(ks),varsigma(i,ks),Qts)
           !!          Qts = Awidg*max(0.0,taumax-varsigma(i,ks)*taucr(ks))*uv(i)
           Qts = rhosed*Qts
           CtstarP(i,ks) = Qts/max(uv(i)*h(i),small) !Total-load concentration capacity of total load, kg/m^3              
           !Use Van Rijn to get rsk
           call sedtrans_wavcur_vanrijn(grav,h(i),uv(i),Worb(i),Wper(i),&
                specgrav,diam(ks),diam(ks),dstar(ks),varsigma(i,ks),Qbs,Qss)
           Qts = Qbs + Qss       !Total Load Capacity   
           rsk(i,ks) = Qss/max(Qts,small) !Fraction of suspended sediment      
           CtstarP(i,ks) = (scalebed*(1.0-rsk(i,ks))+scalesus*rsk(i,ks))*CtstarP(i,ks)
           CtstarP(i,ks) = min(CtstarP(i,ks),Cteqmax)
        enddo !ks
     enddo !i
     !$OMP END PARALLEL DO      
  else  !No Waves
     !$OMP PARALLEL DO PRIVATE(i,ks,phi,Qbs,Qss,Qts,taumax)
     do i=1,ncells 
        if(iwet(i)*uv(i)<=1.0e-4)then
           CtstarP(i,:)=0.0
           rsk(i,:)=1.0
           cycle 
        endif
        do ks=1,nsed                    
           call sedtrans_wavcur_watanabe(uv(i),bsxy(i),taucr(ks),varsigma(i,ks),Qts)
           !!          Qts = Awidg*max(0.0,bsxy(i)-varsigma(i,ks)*taucr(ks))*uv(i)
           Qts = rhosed*Qts
           CtstarP(i,ks) = Qts/max(uv(i)*h(i),small) !Total-load concentration capacity of total load, kg/m^3              
           !Use Van Rijn to get rsk
           call sedtrans_cur_vanrijn(grav,h(i),uv(i),&
                specgrav,diam(ks),diam(ks),dstar(ks),varsigma(i,ks),Qbs,Qss)
           Qts = Qbs + Qss       !Total Load Capacity      
           rsk(i,ks) = Qss/max(Qts,small) !Fraction of suspended sediment        
           CtstarP(i,ks) = (scalebed*(1.0-rsk(i,ks))+scalesus*rsk(i,ks))*CtstarP(i,ks)
           CtstarP(i,ks) = min(CtstarP(i,ks),Cteqmax)
        enddo
     enddo
     !$OMP END PARALLEL DO      
  endif

  return
end subroutine sedcapac_watanabe

!********************************************************************
subroutine sedcapac_c2shore
  ! Calculates the transport capacity based on the C2SHORE model
  ! written by Brad Johnson, USACE-CHL
  !********************************************************************    
  use size_def
  use geo_def, only: mapid,dzbx,dzby,x,y,cell2cell
  use flow_def
  use const_def
  use wave_flowgrid_def
  use sed_def
  use sed_lib, only: critslpcor_dey
  use rol_def, only: roller
  use cms_def
  use prec_def
  implicit none
  integer :: i,ks,iripple
  real(ikind) :: tauct,tauwt,tauwmt,taucwt,taucwmt,tauctb,taucwtb,taucwmtb
  real(ikind) :: fcf,fwf,fcwf,BDpart,Ustc,Ustw,Hrms,sigT,phi,alfa,ur,gamma,um
  real(ikind) :: Qss,Qbs,Qsm,Qsr,Qbm,Qba,Qts,Uw,T,dbrk,fac
  real(ikind) :: CSPs,CSPb,CSDf,CSDb,CSefff,CSwf,CSsg,CSVs,qb
  !real(ikind) :: CSPs,CSPb,CSDf,CSDb,CSefff,CSeffb,CSwf,CSsg,CSVs      !CSeffb now defined in sed_def and initialized in sed_default - bdj 6/7/19

  iripple = 1               
  gamma = 0.78      
  !Parallel statements added 6/20/2018 - meb  
!!$OMP PARALLEL DO PRIVATE (i,ks,CSDb,CSefff,CSeffb,CSsg,CSwf,Hrms,sigT,CSDf,CSPs,CSVs)  ! this parr loop commented by bdj 2020-12-18 
  do i=1,ncells     
    if(iwet(i)==0)then
      CtstarP(i,:)=0.0
      rsk(i,:)=1.0
      cycle 
    endif

    do ks=1,nsed
      CSDb = max(wavediss(i),0.)
      CSefff = 2.*CSeffb
      CSsg = rhosed/1000.
      CSwf = wsfall(ks)
      Hrms = Whgt(i)/sqrt(2.)  
      sigT = (Hrms/sqrt(8.))*(Wlen(i)/Wper(i))/h(i)
      call get_CSDf(u(i),v(i),sigT,Wang(i),Wper(i),CSwf,CSDf)
      call prob_susload(u(i),v(i),sigT,Wang(i),Wper(i),CSwf,CSPs)
      call prob_bedload(sigT,Wper(i),CSsg,diam(1),u(i),v(i),CSPb)
      qb = rhosed*(CSPb*CSblp*sigT**3.)/(9.81*(CSsg-1.))
      CSVs = CSPs*(CSDf*CSefff + CSDb*CSeffb)/(9810.*(CSsg-1)*CSwf);
      CtstarP(i,ks) = rhosed*CSVs/(h(i)+small) !changing C2SHORE convention of depth integrated 
    ! volumetric concentration Vs [m] to mass concentration CStarP [ kg/m^3]              
      if(wavesedtrans)then
        ! following the previous convention, wave-related transport is BY DEFINITION in direction of wave propagation  
        ! Asymmetry related will be pos and return current related will be neg
        !QwsP(i,ks) = rhosed*(.0000001)*x(i) !Potential net onshore transport, kg/m/sec
        !QwsP(i,ks) = (1-CSslp)*sqrt(us(i)**2.+vs(i)**2.)*h(i)*CtstarP(i,ks) !
        QwsP(i,ks) = -CSslp*sqrt(us(i)**2.+vs(i)**2.)*h(i)*CtstarP(i,ks) !only return-current transport here, kg/m/sec
        QwsP(i,ks) = QwsP(i,ks) + qb ! the addition of bedload
        !QwsP(i,ks) = -CSslp*sqrt(us(i)**2.+vs(i)**2.)*h(i)*1 !commented 2021-01-08
        !QwsP(i,ks) = -rhosed*(csslp)*us(i) !Potential net onshore transport, kg/m/sec
        endif
     enddo
   enddo
!!$OMP END PARALLEL DO  

  return
end subroutine sedcapac_c2shore

subroutine prob_susload(u,v,sigT,alpha,Tp,CSwf,CSPs)
  ! calculates the probability of sediment suspention, Ps
  ! written by Brad Johnson, USACE-CHL;
  !************************************************************************      
  use prec_def
  implicit none
  integer :: i,numsteps
  real(ikind) :: u,v,sigT,alpha,Tp,CSwf,CSPs
  real(ikind) :: fw,rho,mag_r,r,f,dr,Uwc,Vwc,Ua,diss

  fw = 0.02  
  rho = 1000.;
  numsteps = 100 
  mag_r = 5.
  dr = 2.*mag_r/numsteps
  CSPs = 0.
  do i = 1,numsteps
     r = -mag_r + 2.*mag_r*(float(i)-1.)/float(numsteps-1)
     f = 1./sqrt(2.*3.14)*exp(-.5*r**2.);
     Uwc = 1.*abs(u)+1.*sigT*cos(alpha)*r;
     Vwc = 1.*abs(v)+1.*sigT*sin(alpha)*r;
     Ua = sqrt(Uwc**2.+Vwc**2.);
     diss = .5*fw*rho*Ua**3.;
     if(((diss/rho)**(0.33333)).gt.CSwf) then
        CSPs = CSPs+dr*f
     endif
  enddo
  return
end subroutine prob_susload

subroutine get_CSDf(u,v,sigT,alpha,Tp,CSwf,CSDf)
  ! calculates the energy dissipation in the BBL
  ! written by Brad Johnson, USACE-CHL;
  !************************************************************************      
  use prec_def
  implicit none
  integer :: i,numsteps
  real(ikind) :: u,v,sigT,alpha,Tp,CSwf,CSDf
  real(ikind) :: fw,rho,mag_r,r,f,dr,Uwc,Vwc,Ua,diss

  fw = 0.02  
  rho = 1000.;
  numsteps = 100 
  mag_r = 5.
  dr = 2.*mag_r/numsteps
  CSDf = 0.
  do i = 1,numsteps
     r = -mag_r + 2.*mag_r*(float(i)-1.)/float(numsteps-1)
     f = 1./sqrt(2.*3.14)*exp(-.5*r**2.);
     Uwc = 1.*abs(u)+1.*sigT*cos(alpha)*r;
     Vwc = 1.*abs(v)+1.*sigT*sin(alpha)*r;
     Ua = sqrt(Uwc**2.+Vwc**2.);
     diss = .5*fw*rho*Ua**3.;
     CSDf = CSDf + dr*diss*f;
  enddo

  return
end subroutine get_CSDf


!******************************************************************  
subroutine cohsedentrain
!   Cohesive sediment entrainment rate
!   by W. Wu, Clarkson Univ.
!******************************************************************          
  use size_def
  use geo_def
  use flow_def
  use fric_def 
  use const_def
  use wave_flowgrid_def
  use sed_def
  use cms_def
  use prec_def
  implicit none
  integer :: i,k,nck,kkdf,kkk
  real(ikind) :: aw,taubc,taubw,taub,coefksr,coefks,riplenw,riphegw,Urms,phi,fw,worbi  !worbind
  real(ikind) :: bedchangei,taubgr,taubcgr,taubwgr,coefksgr,fwgr,gammaw
  
  !Calculate crtical shear stress for erosion, tau_ce, for cohesive sediment
  if(methcoherodcr.eq.2) then   !depth-varying tau_ce, pares of tau_cr and depth below the bed
    do i=1,ncells
      bedchangei=zb0(i)-zb(i)-0.5*dzb(i)
      if(bedchangei.le.coheroddep(1)) coherodcr(i)=coherodtauce(1)
      do k=2,numbtaucedep
        if((bedchangei.gt.coheroddep(k-1)).and.(bedchangei.le.coheroddep(k))) then
          coherodcr(i)=coherodtauce(k-1)+(coherodtauce(k)-coherodtauce(k-1))*  &
                           (bedchangei-coheroddep(k-1))/(coheroddep(k)-coheroddep(k-1))
        endif
      enddo
      if(bedchangei.gt.coheroddep(numbtaucedep)) coherodcr(i)=coherodtauce(numbtaucedep)
    enddo
  elseif(methcoherodcr.eq.4) then  !function of mud dry density 
    do i=1,ncells
      coherodcr(i)=coherodcralpha*rhobedcoh(i,1)**coherodcrbeta
    enddo
  elseif(methcoherodcr.eq.5) then  !function of excess mud dry density by Nicholson and O'Conor (1986)
    do i=1,ncells
      coherodcr(i)=coherodcratrho0+coherodcrtau*(rhobedcoh(i,1)-rhobedcohercr0)**coherodcrn
    enddo
  endif
    
  !Calculate the entrainment rate for cohesive sediment, E
  if(cmswave)then  !Waves 
    do i=1,ncells
      phi=abs(wang(i)-atan2(v(i),u(i)))     !Current-wave angle              	
      !phi=wang(i)-atan2(v(i),u(i))     !Current-wave angle              	
      worbi=Worb(i)
      !worbind=1.0
      do k=1,ncface(i)
        nck=cell2cell(k,i)   !loconect(i,k)
        if(iwet(nck).eq.0) then
          !worbi=0.0; worbind=0.0;  kkdf=kkface(idirface(i,k))
          !do kkk=1,ncface(i)    !Assume boundary face does not split
          !   if(idirface(i,kkk).eq.kkdf) then
          !      !worbi=worbi+Worb(loconect(i,kkk))
          !      worbi=worbi+Worb(cell2cell(kkk,i))
          !      worbind=worbind+1.0
          !   endif
          !enddo                   
          !worbi=worbi/worbind
          worbi=worbi*0.5 
        endif   
      enddo
      Urms=worbi/1.41421356      !sqtwo  
      aw=Urms*Wper(i)/2.0/pi  !Wave excursion
      aw=max(0.00000000001, aw)

      !riplenw=aw/(1.0+0.00187*aw/d50(i)*(1.0-exp(-(0.0002*aw/d50(i))**1.5)))  !Soulsby and Whitehouse (2005)
      !  !!riphegw=0.15*(1.0-exp(-(5000.0*d50(i)/aw)**3.5))*riplenw
      !  !!coefksr=12.0*riphegw**2/riplenw  !Form roughness
      !riphegw=0.15*(1.0-exp(-(min(10.0,5000.0*d50(i)/aw))**3.5))   !Delta/L
      !coefksr=12.0*riphegw**2*riplenw  !Form roughness   

      gammaw=Urms**2/((rhosed/rhow-1.0)*grav*d50(i))
      if(gammaw.le.10.0) then      !ripple geometry, Van Rijn (1993)
        riphegw=0.22*aw
        coefksr=12.0*riphegw*0.18
      elseif((gammaw.gt.10.0).and.(gammaw.le.250.0)) then
        riphegw=0.0028*aw*(2.5-0.01*gammaw)**5
        coefksr=12.0*riphegw*0.02*(2.5-0.01*gammaw)**2.5
      elseif(gammaw.gt.250.0) then
        coefksr=0.0
      endif            
      
      coefks=d50(i)+coefksr
        !coefks=1.5*d90(i)+coefksr
        !!coefks=1.5*d90(i)+coefksr*worbind         
      coefksgr=d50(i)    !Grain roughness
        !!coefksgr=1.5*d90(i)    !Grain roughnes

      !fw=0.237*(aw/coefks)**(-0.52)   !Soulsby's wave friction coefficient   
      fw=min(0.237*(coefks/aw)**0.52, 0.15)   !Soulsby's wave friction coefficient  
      taubw=0.25*fw*Urms**2*rhow    !Bed shear stress by current    
      fwgr=min(0.237*(coefksgr/aw)**0.52, 0.15)   !Soulsby's wave friction coefficient  
      taubwgr=0.25*fwgr*Urms**2*rhow    !Bed shear stress by current    
                
      taubc=bsxy(i)    !rhow*grav*(abs(uv(i))*coefman(i))**2*h(i)**(-0.3333333)  !Bed shear stress by current	
      taubcgr=taubc*(cohmangrain/coefman(i))**1.5    ! grain bed shear stress by current	

      taub=sqrt(taubc**2+taubw**2+2.0*taubc*taubw*cos(phi))
      cohbsxy(i)=max(taub, 0.00000000001)

      taubgr=sqrt(taubcgr**2+taubwgr**2+2.0*taubcgr*taubwgr*cos(phi))
      EtstarP(i,1)=coherodm(i)*(max(0.0, taubgr/coherodcr(i)-1.0))**coherodn
        !EtstarP(i,1)=coherodm(i)*(max(0.0, cohbsxy(i)/coherodcr(i)-1.0))**coherodn
      rsk(i,1)=1.0
    enddo
  else
    do i=1,ncells
      cohbsxy(i)=max(bsxy(i), 0.00000000001) 
      taubcgr=max(bsxy(i)*(cohmangrain/coefman(i))**1.5, 0.00000000001)    !Grain bed shear stress
      EtstarP(i,1)=coherodm(i)*(max(0.0, taubcgr/coherodcr(i)-1.0))**coherodn
        !EtstarP(i,1)=coherodm(i)*(max(0.0, cohbsxy(i)/coherodcr(i)-1.0))**coherodn
      rsk(i,1)=1.0
    enddo
  endif
    
  return
endsubroutine cohsedentrain
    

