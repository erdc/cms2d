      subroutine ST_cohesive()
	use EXP_Global_def 
      USE EXP_bndcond_def
      USE EXP_transport_def 
      use sed_def
      use wave_flowgrid_def
      use flow_def
      use comvarbl
      use const_def, only: pi
      use size_def 
      use geo_def, only: zb,cell2cell
                
      implicit none
      integer i
      real(ikind) tdepth,uvel,vvel,uc2,fc_exp,tc,uw,aw,re,fw
      real(ikind) ratio
      real(ikind) funcws,ustc_exp,ustw,argg,tw,ws,dcc,dcw
          
      !cohesive sediment transport

      if(waves) then
!!$omp parallel do 
!!$omp+ private(tdepth,uvel,vvel,U2,fc_exp,tc,Aw,Re,fw,tw,funcWS)
!!$omp+ private(ws,ustc_exp,USTW,DCC,DCW,ratio)
        do i=1,ncells
          tdepth = -ZB(i) + etan(i)
          uvel = ((qxn(i)+qxn(cell2cell(2,i)))/2.)/tdepth
          vvel = ((qyn(i)+qyn(cell2cell(1,i)))/2.)/tdepth        
          Uc2 = uvel**2+vvel**2	
          fc_exp = 0.24/((log10(12*tdepth/(0.0001)))**2)
          tc = rhowdiv8*fc_exp*uc2
          argg = twopi*tdepth/Wlen(i)
          argg=min(argg,50.0)
          Uw = pi*Whgt(i)/(Wper(i)*sinh(argg))
          Aw = Uw*Wper(i)/(twopi)
          Re = Uw*Aw/viscos
          fw = 0.0521/(Re+100.0)**0.187
          tw = rhowdiv2*fw*Uw**2
          COHES(i)%tbmax = sqrt(tc**2+tw**2)    !added (i) - 05/15/09 reed
          COHES(i)%depo = 0
          COHES(i)%eros = 0	
          if(COHES(i)%tbmax.le.COHES(i)%tcrit_D) then
            funcWS =max(1.0 - abs(COHES(i)%conc-CHparms%c_peak)/CHparms%c_max,0.0)
            ws= chparms%ws_max*funcWS
            COHES(i)%depo = COHES(i)%conc*ws
          elseif(COHES(i)%tbmax.gt.COHES(i)%tcrit_E) then
            COHES(i)%eros = CHparms%E*(COHES(i)%tbmax-COHES(i)%tcrit_E)
          endif
          ustc_exp = sqrt(fc_exp/8)*Uc2
          USTW = sqrt(fw/2)*Uw
          DCC = 5.93*tdepth*ustc_exp
          DCW = 0.5*USTW*tdepth  !gamma = 0.5
          ratio = (Whgt(i)/tdepth)**3
          COHES(i)%DiffC = CHparms%Dfac*((1.0-ratio)*DCC + ratio*DCW)
        enddo
!!$omp end parallel do
      else
!!$omp parallel do 
!!$omp+ private(tdepth,uvel,vvel,U2,fc_exp,tc,Aw,Re,fw,tw,funcWS)
!!$omp+ private(ws,ustc_exp,USTW,DCC,DCW,ratio)
        do i =1,ncells
          tdepth = -ZB(i) + etan(i)
          uvel = ((qxn(i)+qxn(cell2cell(2,i)))/2.)/tdepth
          vvel = ((qyn(i)+qyn(cell2cell(1,i)))/2.)/tdepth	
          Uc2 = uvel**2+vvel**2	
          fc_exp = 0.24/((log10(12*tdepth/(0.0001)))**2)
          COHES(i)%tbmax = rhowdiv8*fc_exp*uc2
          COHES(i)%depo = 0
          COHES(i)%eros = 0	
          if(COHES(i)%tbmax.le.COHES(i)%tcrit_D) then
            funcWS =max(1.0 - abs(COHES(i)%conc-CHparms%c_peak)/CHparms%c_max,0.0)
            ws= CHparms%ws_max*funcWS
            COHES(i)%depo = COHES(i)%conc*ws
          elseif(COHES(i)%tbmax.gt.COHES(i)%tcrit_E) then
            COHES(i)%eros = CHparms%E*(COHES(i)%tbmax-COHES(i)%tcrit_E)
          endif
          ustc_exp = sqrt(fc_exp/8)*Uc2	
          COHES(i)%diffC = CHparms%Dfac*5.93*tdepth*ustc_exp	
        enddo	
!!$omp end parallel do        
      endif
      
      end subroutine
