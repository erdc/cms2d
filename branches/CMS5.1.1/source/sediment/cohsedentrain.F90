!******************************************************************  
    subroutine cohsedentrain
!   Cohesive sediment entrainment rate
!******************************************************************          
    use case_size
    use fl2d
    use comvarbl
    use wave_flgrid
    use flowglobal
    use sedmod
    use option
    use precision
    implicit none
    integer :: kkdf,kkk,nck
    real(ikind) :: aw,taubc,taubw,taub,coefksr,coefks,riplen,ripheg,Urms,phi,fw,worbind,worbi

    if(methcoherodcr.eq.1) then   !Critical shear stress for erosion by Nicholson and O'Conor (1986)
       do i=1,ncells
          coherodcr(i)=coherodcratrho0+coherodcrtau*(rhobedcoh(i,1)-rhobedcohercr0)**coherodcrn
       enddo
    endif
    
    if(cmswave)then  !Waves 

       do i=1,ncells
          phi=abs(wang(i)-atan2(v(i),u(i)))     !Current-wave angle              	
          worbi=Worb(i)
          do k=1,ncface(i)
             nck=loconect(i,k)
             if(idry(nck).eq.0) then
                !worbi=0.0; worbind=0.0;  kkdf=kkface(idirface(i,k))
                !do kkk=1,ncface(i)    !Assume boundary face does not split
                !   if(idirface(i,kkk).eq.kkdf) then
                !      worbi=worbi+Worb(loconect(i,kkk))
                !      worbind=worbind+1.0
                !   endif
                !enddo                   
                !worbi=worbi/worbind
                worbi=worbi*0.5
             endif   
          enddo
          Urms=worbi/1.41421356      !sqtwo  
          aw=Urms*Wper(i)/2.0/pi  !Wave excursion
          aw=max(0.000001, aw)

          riplen=aw/(1.0+0.00187*aw/d50(i)*(1.0-exp(-(0.0002*aw/d50(i))**1.5)))
          ripheg=0.15*(1.0-exp(-(5000.0*d50(i)/aw)**3.5))*riplen
          coefksr=12.0*ripheg**2/riplen  !Form roughness
          coefks=1.5*d90(i)+coefksr

          taubc=densit*grav*(abs(uv(i))*coefman(i))**2*h(i)**(-0.3333333)  !Bed shear stress by current	

          fw=0.237*(aw/coefks)**(-0.52)   !Soulsby's wave friction coefficient   
          taubw=0.25*fw*Urms**2*densit    !Bed shear stress by current    
          taub=sqrt(taubc**2+taubw**2+2.0*taubc*taubw*cos(phi))
          cohbsxy(i)=taub

          EtstarP(i,1)=coherodm(i)*(max(0.0, taub/coherodcr(i)-1.0))**coherodn
          rsk(i,1)=1.0
       enddo

    else

       do i=1,ncells
          cohbsxy(i)=bsxy(i)
          EtstarP(i,1)=coherodm(i)*(max(0.0, bsxy(i)/coherodcr(i)-1.0))**coherodn
          rsk(i,1)=1.0
       enddo

    endif

    return
    endsubroutine
    
