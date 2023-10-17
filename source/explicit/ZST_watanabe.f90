!***********************************************************************
      subroutine ST_watanabe()
!***********************************************************************
      use EXP_Global_def,    only: ncn, nce, thetac, etan, cdx, cdy, waves
      USE EXP_transport_def, only: xks, tcr, rhowdiv2, rhowdiv8, twopi, a0divrhowgrav, qsx, qsy
      use wave_flowgrid_def, only: whgt, wper, wlen
      use sed_def,   only: rhosed,d50 
      use flow_def,  only: rhow, grav
      use prec_def,  only: ikind
      use size_def,  only: ncells
      use const_def, only: pi
      use geo_def,   only: zb,cell2cell 
       
      implicit none
      !local vars
      integer i
      real(ikind) tdepth,uc2,fcEXP,tcEXP,uw,aw,rfac,fw,tw,tbmaxEXP,uc
      real(ikind) qtot,uvel,vvel,ARG,singled50
     
      call update_Q_forSedTrans()   !this fixes WSE BC cells lower vlocity from reducing Sed Trans when
                                    !extrap < 1.0 is used    
      
      !watanabe formula for bed load

!$omp parallel                           !NLH 07/24/2008
!$omp do private(tdepth,uvel,vvel,uc2,fcEXP,tcEXP,uw,aw,rfac,arg,fw,tw,tbmaxEXP,uc,qtot,xks,tcr,SINGLED50,ncn,nce)
      do i=1,ncells
        SINGLED50 = D50(I)
        xks = 2.5*SINGLED50                          
        tcr = (rhosed-rhow)*grav*SINGLED50*thetac
        NCN=cell2cell(1,i)
        NCE=cell2cell(2,i)
        tdepth = -zb(i) + etan(i)
        uvel = cdx(i) !((qxn(i)+qxn(NCE))/2.)/tdepth
        vvel = cdy(i) !((qyn(i)+qyn(NCN))/2.)/tdepth    
        Uc2 = uvel**2+vvel**2
        if(waves) then
          fcEXP = 0.24/(log10(12*tdepth/(xks)))**2
          tcEXP = rhowdiv8*fcEXP*uc2
          !WhgtT = min(whgt(i),0.78*(etan(i)+depth(i)))
          Uw = pi*Whgt(i)/(Wper(i)*sinh(twopi*tdepth/Wlen(i)))
          Aw = Uw*Wper(i)/(twopi)
          Rfac = max(Aw/(xks),1.0)
          arg = -6.3 + 5.5/Rfac**0.2
          fw = exp(arg)
          tw = rhowdiv2*fw*Uw**2
          tbmaxEXP = tcEXP*(1 + 1.2*(tw/(tw+tcEXP+1.e-20))**3.2)
        else
          fcEXP = 0.24/log10(12*tdepth/(xks))**2
          tbmaxEXP = rhowdiv8*fcEXP*uc2
        endif
        Uc = sqrt(Uc2) + 1.e-20
        qtot = A0divrhowgrav*max(0.0,tbmaxEXP-tcr)*Uc
        qsx(i) = (uvel/Uc)*qtot
        qsy(i) = (vvel/Uc)*qtot
      enddo
!$OMP end do
!$OMP END PARALLEL

      return
      end subroutine
