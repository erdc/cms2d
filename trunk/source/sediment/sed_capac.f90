!========================================================================
! CMS Sediment Transport Capacity/Formulas
!
! Contains the following:
!     sedcapac_soulsby - Calculates the concentration capacity using the 
!                    Soulsby-van Rijn transport equations
!     sedcapac_vanrijn - Calculates the concentration capacity using the 
!                    Van Rijn transport equations
!     sedcapac_watanabe - Calculates the concentration capacity using the 
!                    Watanabe transport equation
!
! written by Alex Sanchez, USACE-CHL
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
    real(ikind) :: Qss,Qbs,Qsm,Qsr,Qbm,Qba,Qts,Uw,T,dbrk,fac
    !logical :: isnankind

    iripple = 1               
	if(noptset>=3)then  !Waves
	  gamma = 0.78	  
!	  rolfac = 2.0/rhow/sqrt(grav)
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
!          call brwratio(h(i),Hrms,gamma,alfa)
          call brkwavratio(h(i),Hrms,Wlen(i),gamma,alfa)      
!          fac = 0.5+0.5*cos(pi*min(max(h(i)/max(Whgt2(i),0.01)-10.0,0.0)/10.0,1.0))
!          fac=0.5+0.5*cos(pi*min(max(max(h(i)/max(Whgt2(i),0.01)-1.5,0.0),1.0))
          fac=0.5+0.5*cos(pi*min(max(h(i)/max(Whgt2(i),0.01)-1.0,0.0)/2.0,1.0))
          !fac=1.0
!          if(roller) Qr = rolfac*Esr(i)/sqrt(max(h(i),hmin)) !Surface roller flux
        endif
	    do ks=1,nsed
          call susplund(h(i),uv(i),Worb(i),rhosed,rhow,diam(ks),   &
              wsfall(ks),dstar(ks),tauct,tauwt,tauwmt,taucwt,taucwmt,&
              fcf,fwf,wavediss(i),taucr(ks),cak(i,ks),epsvk(i,ks),Qss,BDpart,Ustc,Ustw)                 
          call bedlund(rhosed,rhow,tauctb,taucwtb,taucwmtb,taucr(ks),Qbs)
!!          call bedlund(rhosed,rhow,tauct,taucwt,taucwmt,taucr(ks),Qbs)
          Qbs = scalebed*Qbs/varsigma(i,ks)       !Bed Load Capacity        
          Qss = scalesus*Qss/varsigma(i,ks)       !Suspended Load Capacity 
          Qts = Qbs + Qss          !Total Load Capacity
	      CtstarP(i,ks) = rhosed*Qts/(uv(i)*h(i)+small) !Total-load concentration capacity of total load, kg/m^3              
          CtstarP(i,ks) = min(CtstarP(i,ks),Cteqmax)
!          rsk(i,ks) = Qss/(Qts+small) !Fraction of suspended sediment  
!          if(rsk(i,ks).lt.1.e-6) rsk(i,ks) = 1.0  !*************        
          rsk(i,ks) = max(Qss,small)/max(Qts,small) !Fraction of suspended sediment      
                    
          !Wave induced sediment transport
          if(wavesedtrans .and. Whgt(i).gt.2*hmin .and. h(i).gt.3*hmin)then    
!            if(.not.roller)then                    
              call crossmean(h(i),Whgt(i),Wper(i),Wlen(i),alfa,&
                 rhosed,rhow,diam(ks),wsfall(ks),            &
                 fcf,taucr(ks),taucwtb,taucwmtb,cak(i,ks),epsvk(i,ks),      &
                 um,ur,Qsm,Qsr,Qbm)
!            else   
!              call crossmean_rol(h(i),Hrms,Wper(i),Wlen(i),Qr, &
!                 rhosed,rhow,diam(ks),wsfall(ks),            &
!                 fcf,taucrks,taucwtb,taucwmtb,crcw,epscw,      &
!                 um,ur,Qsm,Qsr,Qbm)   
!            endif     
            call asymmetry(h(i),Hrms,Wper(i),Wlen(i),Worbrep(i),&
               rhosed,rhow,diam(ks),fwf,taucr(ks),taucwtb,taucwmtb,Qba)
            QwsP(i,ks) = rhosed*(fba*Qba+fac*(fsr*Qsr-fsm*Qsm-fbm*Qbm)) !Potential net onshore transport, kg/m/sec
!            QwsP(i,ks) = min(0.0001,QwsP(i,ks)) !Limit transport to unrealistic transport near wet/dry interface
!            QwsP(i,ks) = max(-0.0001,QwsP(i,ks)) !Limit transport to unrealistic transport near wet/dry interface
            !if(isnankind(QwsP(i,ks)))then
            !  write(*,*)
            !  write(*,*) 'ERROR: Qws(i,ks)=NaN'
            !  write(*,*) 'i,mapid(i),ks',i,mapid(i),ks
            !  write(*,*) 'h(i)',h(i)
            !  write(*,*) 'um,ur',um,ur
            !  write(*,*) 'cak(i,ks),epsvk(i,ks)',cak(i,ks),epsvk(i,ks)
            !  write(*,*) 'Whgt(i),Worbrep(i),Wper(i)',Whgt(i),Worbrep(i),Wper(i)
            !  write(*,*) 'Qsr,Qba,Qsm,Qbm',Qsr,Qba,Qsm,Qbm
            !  write(*,*) 'Press any key to continue.'
            !  read(*,*)
            !  stop
            !endif
          endif
        enddo
      enddo
	else         !No waves
	  phi = 0.0; Uw = 0.0; T = 10.0; dbrk = 0.0
	  do i=1,ncells  
        if(iwet(i)==0)then
          CtstarP(i,:)=0.0
          rsk(i,:)=1.0
          cycle 
        endif   
	    call shearlund(iripple,h(i),uv(i),Uw,T,phi,rhosed,rhow,D50(i),&
             tauct,tauwt,tauwmt,taucwt,taucwmt,tauctb,taucwtb,taucwmtb,fcf,fwf,fcwf)
	    do ks=1,nsed
          call susplund(h(i),uv(i),Uw,rhosed,rhow,diam(ks),       &
             wsfall(ks),dstar(ks),tauct,tauwt,tauwmt,taucwt,taucwmt,&
             fcf,fwf,dbrk,taucr(ks),cak(i,ks),epsvk(i,ks),Qss,BDpart,Ustc,Ustw)                
          call bedlund(rhosed,rhow,tauctb,taucwtb,taucwmtb,taucr(ks),Qbs)                            
          Qbs = scalebed*Qbs/varsigma(i,ks)       !Bed Load Capacity        
          Qss = scalesus*Qss/varsigma(i,ks)       !Suspended Load Capacity 
          Qts = Qbs + Qss          !Total Load Capacity     
          CtstarP(i,ks) = rhosed*Qts/max(uv(i)*h(i),small) !Total-load concentration capacity of total load, kg/m^3  
          CtstarP(i,ks) = min(CtstarP(i,ks),Cteqmax)
!          rsk(i,ks) = Qss/(Qts+small)  !max(Qts,small)  !Fraction of suspended sediment                
!          if(rsk(i,ks).lt.0.0001) rsk(i,ks) = 1.0  !*************        
          rsk(i,ks) = max(Qss,small)/max(Qts,small)
        enddo
      enddo	
	endif   

    return
    endsubroutine sedcapac_lundcirp    
    
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
    endsubroutine sedcapac_soulsby 
     
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
    endsubroutine sedcapac_vanrijn
     
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
          !if(isnankind(rsk(i,ks)))then
          !  continue
          !endif
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
    endsubroutine sedcapac_watanabe