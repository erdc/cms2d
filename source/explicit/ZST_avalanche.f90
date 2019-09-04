
      subroutine ST_avalanche()
	use EXP_Global_def 
      USE EXP_bndcond_def
      USE EXP_transport_def     
      use flow_def 
      use comvarbl     
      use sed_def
      use size_def
      use geo_def, only: dx,dy,zb,cell2cell
 
      !USE SYNOPTIC_VARS    !Alex, Sep 23, 2009
      implicit none
      !local vars 
      integer i
      real slope,q_aval    

      !avalanching - add avalanche flux to qsx or qsy 
      !rate is set so that avalanche of 5% of excess slope 
      !will be reduced in the next time step

!$OMP PARALLEL DO PRIVATE (NCE,NCN,SLOPE,Q_AVAL) 
      do i=1,ncells
        nce = cell2cell(2,i)
        ncn = cell2cell(1,i)
        slope = 2*(-zb(nce)+zb(i))/(dx(i)+dx(nce))
        if(slope.gt.a_repose) then  !transport will occurr
          q_aval = (slope-a_repose)  
          q_aval = min(q_aval,rate_avalanche)          !/(1./dx(i)+1./dx(nce)) 
          if(iwet(nce).eq. 1 .or. iwet(i) .eq. 1 ) then
            qsx(i) = qsx(i) + q_aval  !*0.01/tsed_elapse  !cr 8/25/08
          endif
        elseif(slope.lt.-a_repose) then
          q_aval = (a_repose+slope) !  /(1./dx(i)+1./dx(nce)) 
          q_aval = max(q_aval,-rate_avalanche)
          if(iwet(nce) .eq. 1 .or. iwet(i) .eq. 1 ) then
            qsx(nce) = qsx(nce) + q_aval !*0.01/tsed_elapse                !cr 8/25/08
          endif
        endif
        
        slope = 2*(-zb(ncn)+zb(i))/(dy(i)+dy(ncn))
        if(slope.gt.a_repose) then  !transport will occur
          q_aval = (slope-a_repose) !*(dy(i)+dy(ncn)))
          q_aval = min(q_aval,rate_avalanche)  !    /(1./dy(i)+1./dy(ncn)) 
          if(iwet(ncn) .eq. 1 .or. iwet(i) .eq. 1 ) then
            qsy(i) = qsy(i) + q_aval  !*0.01  !/tsed_elapse                    !cr 8/25/08
          endif
        elseif(slope.lt.-a_repose) then
          q_aval = (a_repose+slope)   !      /(1./dy(i)+1./dy(ncn)) 
          q_aval = max(q_aval,-rate_avalanche)
          if(iwet(ncn) .eq. 1 .or. iwet(i) .eq. 1 ) then
            qsy(ncn) = qsy(ncn) + q_aval  !*0.01  !/tsed_elapse                !cr 8/25/08
          endif
        endif
      enddo
!$OMP END PARALLEL DO    

      end subroutine  
