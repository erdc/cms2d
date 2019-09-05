
!******************************************************************  
    subroutine wucapac
! Calculates the concentration capacity using the Wu transport equations
! written by Weiming Wu, NCCHE,  Oct. 2011; Revised, June 2013
!******************************************************************          
    use size_def, only: ncells
    use flow_def, only: h,u,v,uv,iwet,rhow,hmin,grav
    use fric_def, only: coefman,bsxy,cfrict
    use const_def, only: small,pi,twopi
    use wave_flowgrid_def, only: Wper,Worbrep,Wang
    use sed_def
    use cms_def, only: noptset
    use prec_def
    implicit none
    integer :: i,ks
    real(ikind) :: Qbs,Qss,Qts,taubcrwu,val,cman_p,fac,taubcr,coefaa
    real(ikind) :: uw,tw,aw,wcangle,coefksg,coefksrw,coefksw,fcgr,fwgr !,taubgrc,taubgrw,taubgr
    real(ikind) :: riphegw,riplenw  !,riphegc,riplenc
    real(ikind) :: taubgronsh,taubgroffsh,fcwgr,fw,taub,taubc,taubw
    real(ikind) :: deltau,twc,twt,ac,at,phib,uc,XU,rw
    real(ikind) :: alphaonsh,alphaoffsh,phibonsh,phiboffsh 
    real(ikind) :: uwonsh,uwoffsh,uw2onsh,uw2offsh,ucw2onsh,ucw2offsh
    real(ikind) :: ucwonsh,ucwoffsh,sonsh,soffsh
    logical :: isnankind

    taubcrwu = 0.03*grav*(rhosed-rhow)
    rw = 0.0

    if(noptset>=3)then  !Waves + Current
       do i=1,ncells
          uc=uv(i)
          !wcangle=abs(wang(i)-atan2(v(i),u(i)))     !Current-wave angle                  
          wcangle=wang(i)-atan2(v(i),u(i))     !Current-wave angle                  
          if(wcangle.lt.0.0) wcangle=wcangle+twopi
          !uw = Worb(i)/1.41421356      !sqtwo  
          uw = Worbrep(i) 
          tw = Wper(i)

          coefaa=12.0   !Coefficient Ar
          aw=max(0.01, uw*tw/twopi)  !Wave excursion
         
          coefksg=3.0*d90(i)  !Grain roughness   !For well-sorted sediments
          !coefksg=1.5*d90(i)  !Grain roughness    !For graded sediments
            !!ripple height and length
          riplenw=aw/(1.0+0.00187*aw/d50(i)*(1.0-exp(-(0.0002*aw/d50(i))**1.5)))
          riphegw=0.15*(1.0-exp(-(5000.0*d50(i)/aw)**3.5))*riplenw
          coefksrw=coefaa*riphegw**2/riplenw  !Form roughness
          coefksw=coefksg+coefksrw

          !riplenc=0.245*(d50(i)*1000.0)**0.35   !Raudikivi ripples
          !riphegc=0.074*riplenc/(d50(i)*1000.0)**0.253

            !============current-related grain shear stress==================        
          cman_p = d50(i)**0.1666667/20.0   ! Grain Manning n
          fcgr=2.0*(min(cman_p/coefman(i),1.0))**1.5*cfrict(i)    !skin friction coefficient by current 
          fwgr=0.237*(aw/coefksg)**(-0.52)   !Soulsby's method for skin friction coefficient by waves         

            ! For Bed load
          XU=uc**2/(uc**2+0.5*uw**2+small)
          fcwgr=fcgr*XU+fwgr*(1.0-XU)  !Average skin friction coefficient

          deltau=1.0+8.0*rw**2
          if(rw.lt.0.001) then
             twc=0.5*tw
          else
             twc=tw/pi*acos((sqrt(deltau)-1.0)/4.0/rw)
          endif     
          twt=tw-twc
          ac=pi*twc/tw
          at=pi*twt/tw
          uw2onsh =0.5*uw**2*abs(1.0+rw**2+13.0/6.0*rw*sin(ac)/ac+sin(2.0*ac)/12.0/ac)
          uw2offsh=0.5*uw**2*abs(-1.0-rw**2+13.0/6.0*rw*sin(at)/at-sin(2.0*at)/12.0/at)
          uwonsh =sqrt(uw2onsh)
          uwoffsh=sqrt(uw2offsh)
          ucw2onsh =uc**2+uw2onsh+2.0*uc*uwonsh*cos(wcangle)
          ucw2offsh=uc**2+uw2offsh+2.0*uc*uwoffsh*cos(pi-wcangle)
          taubgronsh =0.5*fcwgr*rhow*ucw2onsh
          taubgroffsh=0.5*fcwgr*rhow*ucw2offsh
                        !Angle between Uc and Ucwm,onshore
          ucwonsh=sqrt(ucw2onsh)
          sonsh=max(0.5*(uwonsh+ucwonsh+uc),ucwonsh,uc)
          alphaonsh=2.0*asin(min(1.0,sqrt((sonsh-uc)*(sonsh-ucwonsh)/(ucwonsh+small)/(uc+small))))
                        !Angle between Uc and Ucwm,offshore
          ucwoffsh=sqrt(ucw2offsh)
          soffsh=max(0.5*(uwoffsh+ucwoffsh+uc),ucwoffsh,uc)
          alphaoffsh=2.0*asin(min(1.0,sqrt((soffsh-uc)*(soffsh-ucwoffsh)/(ucwoffsh+small)/(uc+small))))

             !Bed shear stress for suspended load
          taubc=cfrict(i)*uc**2*rhow !Bed shear stress by current
          fw=0.237*(aw/coefksw)**(-0.52)   !Soulsby's wave friction coefficient   
          taubw=0.25*fw*uw**2*rhow    !Bed shear stress by current    
          taub=sqrt(taubc**2+taubw**2+2.0*taubc*taubw*cos(wcangle))          

          do ks=1,nsed
             taubcr = taubcrwu*diam(ks)*varsigma(i,ks) !For uniform sediment sizes,   Eq. 3.45       
             fac = sqrt(s1grav*diam(ks)**3)
               !Bed Load          
             phibonsh =0.0053*(max(taubgronsh /taubcr-1.0,0.0))**2.2  !Onshore
             phiboffsh=0.0053*(max(taubgroffsh/taubcr-1.0,0.0))**2.2  !Offshore
             phibonsh =phibonsh *twc/tw
             phiboffsh=phiboffsh*twt/tw
             phib=sqrt(max(0.0,phibonsh**2+phiboffsh**2+2.0*phibonsh*phiboffsh*cos(alphaonsh+alphaoffsh)))  
             if(isnankind(phib)) phib=0.0
             Qbs=scalebed*phib*fac*iwet(i)      
               !Suspended load
             val=max((taub/taubcr-1.0),0.0)*uc/wsfall(ks)
             Qss=scalesus*0.0000262*fac*val**1.74 *iwet(i)
             Qts=Qbs+Qss                         
             CtstarP(i,ks)=rhosed*Qts/(uc*max(h(i),3.0*hmin)+small) !Total-load concentration capacity of total load, kg/m^3              
             rsk(i,ks)=Qss/(Qts+small) !Fraction of suspended sediment  
             if(CtstarP(i,ks).gt.Cteqmax) CtstarP(i,ks)=Cteqmax
          enddo 
       enddo 
    else             !Current only
       do i=1,ncells
          do ks=1,nsed
             taubcr = taubcrwu*diam(ks)*varsigma(i,ks) !For uniform sediment sizes,   Eq. 3.45       
             fac = sqrt(s1grav*diam(ks)**3)
               !Bed Load, Wu et al. 2000    
             cman_p = d50(i)**0.166667/20.0   !grain roughness   
             val = max(0.0, (cman_p/coefman(i))**1.5*bsxy(i)/taubcr-1.0)
             Qbs = scalebed*0.0053*fac*val**2.2 *iwet(i) !m^2/s          
               !Suspended Load, Wu et al. 2000    
             val = max(0.0, bsxy(i)/taubcr-1.0)*uv(i)/wsfall(ks)
             Qss = scalesus*0.0000262*fac*val**1.74 *iwet(i)  !m^2/s   
             Qts=Qbs+Qss                         
             CtstarP(i,ks)=rhosed*Qts/(uv(i)*max(h(i),3.0*hmin)+small) !Total-load concentration capacity of total load, kg/m^3              
             rsk(i,ks)=Qss/(Qts+small) !Fraction of suspended sediment  
             if(CtstarP(i,ks).gt.Cteqmax)  CtstarP(i,ks)=Cteqmax
          enddo 
       enddo 
    endif

    return
    endsubroutine

