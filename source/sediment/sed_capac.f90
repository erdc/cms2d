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
  !write(*,*)'bdj in sedcapac_lundcirp , noptset = ',noptset                                                              

  iripple = 1               
  if(noptset>=3)then  !Waves
     gamma = 0.78      
     !      rolfac = 2.0/rhow/sqrt(grav)
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
              !if (i.eq.487) write(*,*)'in sedcapac_lundcirp, now setting QwsP(487,ks)',QwsP(i,ks)                                   !added bdj 01/2021
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

!********************************************************************
subroutine sedcapac_cshore
  ! Calculates the transport capacity based on the CSHORE model
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

  !write(*,*)'bdj in sedcapac_cshore , noptset = ',noptset
  iripple = 1               
  gamma = 0.78      
  !Parallel statements added 6/20/2018 - meb  
!  !!$OMP PARALLEL DO PRIVATE (i,ks,CSDb,CSefff,CSeffb,CSsg,CSwf,Hrms,sigT,CSDf,CSPs,CSVs)  ! this parr loop commented by bdj 2020-12-18 
  do i=1,ncells     
     if(iwet(i)==0)then
        CtstarP(i,:)=0.0
        rsk(i,:)=1.0
        cycle 
     endif

     do ks=1,nsed
        CSDb = max(wavediss(i),0.)
        CSefff = 2.*CSeffb
        !if(i.lt.10) write(*,*) 'CSeffb,CSblp,CSclp = ',CSeffb,CSblp,CSslp
        !CSeffb = .005
        CSsg = rhosed/1000.
        CSwf = wsfall(ks)
        Hrms = Whgt(i)/sqrt(2.)  
        sigT = (Hrms/sqrt(8.))*(Wlen(i)/Wper(i))/h(i)
        call get_CSDf(u(i),v(i),sigT,Wang(i),Wper(i),CSwf,CSDf)
        call prob_susload(u(i),v(i),sigT,Wang(i),Wper(i),CSwf,CSPs)
        call prob_bedload(sigT,Wper(i),CSsg,diam(1),u(i),v(i),CSPb)
        qb = rhosed*(CSPb*CSblp*sigT**3.)/(9.81*(CSsg-1.))
        CSVs = CSPs*(CSDf*CSefff + CSDb*CSeffb)/(9810.*(CSsg-1)*CSwf);
        CtstarP(i,ks) = rhosed*CSVs/(h(i)+small) !changing CSHORE convention of depth integrated 
        ! if (i.eq.487)         write(*,*)'i,CSPs,CSDf,CSDb,CSVs,CSwf',i,CSPs,CSDf,CSDb,CSVs,CSwf
        ! if (i.eq.487)         write(*,*)'x(cell2cell(4,i)),x(cell2cell(2,i))',x(cell2cell(4,i)),x(cell2cell(2,i))
        ! if (i.eq.487)         write(*,*)'h(cell2cell(4,i)),h(i),h(cell2cell(2,i))',h(cell2cell(4,i)),h(i),h(cell2cell(2,i))
        ! if (i.eq.487)         write(*,*)'Hs(cell2cell(4,i)),Hs(i),Hs(cell2cell(2,i))',Whgt(cell2cell(4,i)),Whgt(i),Whgt(cell2cell(2,i))
        ! if (i.eq.487)         write(*,*)'bdj i,x,h,Hrms,CSPb,sigT,qb',i,x(i),h(i),Hrms,CSPb,sigT,qb
        ! if (i.eq.487)         write(*,*)'bdj i,CSDb,CSDf',i,CSPs,CSVs,CSDb,CSDf
        ! volumentric conentration Vs [m] to mass concentration CStarP [ kg/m^3]              
        ! write(*,*)'bdj i, h(i), Hrms, wavediss(i),Vs,CtstarP(i,ks)', i, h(i), Hrms, wavediss(i),CSVs,CtstarP(i,ks)
        ! write(2001,*) i, h(i), Hrms, CtstarP(i,ks), wavediss(i)
        if(wavesedtrans)then
          ! following the previous convention, wave-related trasport id BY DEFINITION in direction of wave propagation  
          ! Asymettry related will be pos and return current related will be neg
          !QwsP(i,ks) = rhosed*(.0000001)*x(i) !Potential net onshore transport, kg/m/sec
          QwsP(i,ks) = -CSslp*sqrt(us(i)**2.+vs(i)**2.)*h(i)*CtstarP(i,ks) !only return-current transport here, kg/m/sec
          !QwsP(i,ks) = (1-CSslp)*sqrt(us(i)**2.+vs(i)**2.)*h(i)*CtstarP(i,ks) !
          QwsP(i,ks) = QwsP(i,ks) + qb ! the addition of bedload
          !if (i.eq.487)         write(*,*)'bdj i,x,h,Hrms,CSPb,sigT,qb',i,x(i),h(i),Hrms,CSPb,sigT,qb
          !if (i.eq.487)         write(*,*)'bdj i,x,y,us,QwsP(i,ks)',i,x(i),y(i),sqrt(us(i)**2.+vs(i)**2.),QwsP(i,ks)
          !QwsP(i,ks) = -CSslp*sqrt(us(i)**2.+vs(i)**2.)*h(i)*1 !commented 2021-01-08

          !QwsP(i,ks) = -rhosed*(csslp)*us(i) !Potential net onshore transport, kg/m/sec
          ! if (i.eq.487) write(*,*)'in sedcapac_cshore, now setting QwsP(487,ks)',QwsP(i,ks)
          ! if (i.eq.487) write(*,*)'in sedcapac_cshore, u(487),us(487),vs(487)',u(487),us(487),vs(487)
        endif
     enddo
   enddo
!  !!$OMP END PARALLEL DO  

  return
endsubroutine sedcapac_cshore

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
     ! write(*,*) CSPs
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
  !write(*,*),'bdj u,v,sigT,alpha,Tp,CSwf,CSPs',u,v,sigT,alpha,Tp,CSwf,CSPs
  return
endsubroutine prob_susload

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
     ! write(*,*) CSPs
     r = -mag_r + 2.*mag_r*(float(i)-1.)/float(numsteps-1)
     f = 1./sqrt(2.*3.14)*exp(-.5*r**2.);
     Uwc = 1.*abs(u)+1.*sigT*cos(alpha)*r;
     Vwc = 1.*abs(v)+1.*sigT*sin(alpha)*r;
     Ua = sqrt(Uwc**2.+Vwc**2.);
     diss = .5*fw*rho*Ua**3.;
     CSDf = CSDf + dr*diss*f;
  enddo

  !write(*,*),'bdj u,v,sigT,alpha,Tp,CSwf,CSDf',u,v,sigT,alpha,Tp,CSwf,CSDf
  return
endsubroutine get_CSDf
