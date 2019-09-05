      subroutine update_wse_bc()
      use EXP_Global_def
      USE EXP_transport_def 
      USE EXP_bndcond_def   
      use flow_def
      use comvarbl, only: timehrs,ramp,dtj
      use bnd_def
      use sed_def
      use geo_def, only: zb
      use met_def, only: windconst,windvar      
      
      
      implicit none
      integer i,j,inc,k,iwse
      real(ikind) helev,hval1,hval2,htval1,htval2,hfac,val1, &
           val2,tval1,tval2,fac,elev,wsebnd
      
      helev = 0.0
      if(MODIFY_H) then                                                 !get to the correct time in the HMOD array.
        do while (timeHRS .gt. HMOD%dtime(HMOD%inc+1))
          HMOD%inc = HMOD%inc + 1
        enddo
        hval1 = HMOD%dvalue(HMOD%inc)
        hval2 = HMOD%dvalue(HMOD%inc+1)
        htval1 = HMOD%dtime(HMOD%inc)
        htval2 = HMOD%dtime(HMOD%inc+1)
        hfac = (timeHRS-htval1) / (htval2-htval1)                       !compute weighting factor for the additional WSE
        helev = hval1 + hfac*(hval2-hval1)                              !compute weighted value of the additional WSE
      endif

      if(nHstr .gt. 0) then
        do i = 1,nHstr  !for each cell string
          !find out where we are in the time/value arrays
          inc = H_str(i)%inc
          do while (timeHRS.gt.H_str(i)%times(inc+1))
            inc = inc+1
            H_str(i)%inc = inc
          enddo
          val1 = H_str(i)%wsecurv(inc)
          val2 = H_str(i)%wsecurv(inc+1) 
          tval1 = H_str(i)%times(inc)
          tval2 = H_str(i)%times(inc+1)    
          fac = (timeHRS-tval1) / (tval2-tval1)                         !compute the weighting factor for the WSE forcing
          elev = val1 + fac*(val2-val1)                                 !compute the weighted value of WSE
          elev = ramp*(elev + helev)                                    !apply the ramp to the elevation plus added WSE from a different signal 
          do j=1,H_str(i)%NCells    !for each cell in string
            eta(H_str(i)%Cells(j)) = elev
            if(eta(H_str(i)%Cells(j))-zb(H_str(i)%Cells(j)).lt.drydep)then
              eta(H_str(i)%Cells(j)) = &
                0.5d0*drydep+zb(H_str(i)%Cells(j))
            endif
            etan(H_str(i)%Cells(j)) = eta(H_str(i)%Cells(j))
          enddo
        enddo ! end of each cell string
      endif  !H_single


      !this has been modified to ignore wave induced phsse difference (psi)and the offset, the 
      !which allows for calc of wse once for all cells
      if(nTHstr .gt. 0) then
      do iwse=1,nTHstr
      !!do j=1,TH_STR(iwse)%ncells
        wsebnd = 0.0  !TH_STR(iwse)%wseoffset
        if(TH_STR(iwse)%istidal)then
          do k=1,TH_STR(iwse)%ntc  
            wsebnd = wsebnd + TH_STR(iwse)%f(k)*TH_STR(iwse)%amp(k) &
                *cos(TH_STR(iwse)%speed(k)*(timehrs+dtj) + TH_STR(iwse)%vu(k) &
                - TH_STR(iwse)%phase(k)) !!+ TH_STR(iwse)%psi(j,k))
          enddo
        else
          do k=1,TH_STR(iwse)%ntc  
            wsebnd = wsebnd + TH_STR(iwse)%amp(k)*cos(TH_STR(iwse)%speed(k)*timehrs &
                - TH_STR(iwse)%phase(k)) !!+ TH_STR(iwse)%psi(j,k))
          enddo
        endif

        do j=1,TH_str(iwse)%NCells    !for each cell in string
          eta(TH_str(iwse)%Cells(j)) = wsebnd
          if(eta(TH_str(iwse)%Cells(j))-  &
             zb(TH_str(iwse)%Cells(j)).lt.drydep)then  
            eta(TH_str(iwse)%Cells(j)) = 0.5d0*drydep+zb(TH_str(iwse)%Cells(j))
          endif
          etan(TH_str(iwse)%Cells(j)) = eta(TH_str(iwse)%Cells(j))
        enddo
        
      !enddo !j-cell
    enddo !iwse-str      
      endif  !H_tide

      if(nMHstr .gt. 0) then
        do i = 1,nMHstr  !for each cell string
          !find out where we are in the time/value arrays
          inc = MH_str(i)%inc
          do while (timeHRS.gt.MH_str(i)%times(inc+1))
            inc = inc+1
            MH_str(i)%inc = inc
          enddo
          tval1 = MH_str(i)%times(inc)    
          tval2 = MH_str(i)%times(inc+1)
          fac = (timeHRS-tval1) / (tval2-tval1)                         !compute the weighting factor for the WSE forcing
          do j=1,MH_str(i)%NCells   !for each cell in string
            val1 = MH_str(i)%wsedata(j,inc)    
            val2 = MH_str(i)%wsedata(j,inc+1)
            elev = val1 + fac*(val2-val1)                               !compute the weighted value of WSE
            elev = ramp*(elev + helev)                                  !apply the ramp to the elevation plus added WSE from a different signal 
            eta(MH_str(i)%Cells(j)) = elev
            if(eta(MH_str(i)%Cells(j))- zb(MH_str(i)%Cells(j)).lt.drydep)then  
              eta(MH_str(i)%Cells(j)) = 0.5d0*drydep+zb(MH_str(i)%Cells(j))
            endif
            etan(MH_str(i)%Cells(j)) = eta(MH_str(i)%Cells(j))
          enddo
        enddo ! end of each cell string
      endif  !H_multi

      if(nMHVstr .gt. 0) then
        do i = 1,nMHVstr  !for each cell string
          !find out where we are in the time/value arrays
          inc = MHV_str(i)%incwse
          do while (timeHRS.gt.MHV_str(i)%timeswse(inc+1))
            inc = inc+1
            MHV_str(i)%incwse = inc
          enddo
          tval1 = MHV_str(i)%timeswse(inc)    
          tval2 = MHV_str(i)%timeswse(inc+1)
          fac = (timeHRS-tval1) / (tval2-tval1)                         !compute the weighting factor for the WSE forcing
          do j=1,MHV_str(i)%NCells   !for each cell in string
            val1 = MHV_str(i)%wsedata(j,inc)    
            val2 = MHV_str(i)%wsedata(j,inc+1)
            elev = val1 + fac*(val2-val1)                               !compute the weighted value of WSE
            elev = ramp*(elev + helev)                                  !apply the ramp to the elevation plus added WSE from a different signal 
            eta(MHV_str(i)%Cells(j)) = elev
            if(eta(MHV_str(i)%Cells(j))- zb(MHV_str(i)%Cells(j)).lt.drydep)then
              eta(MHV_str(i)%Cells(j)) = 0.5d0*drydep+zb(MHV_str(i)%Cells(j))
            endif
            etan(MHV_str(i)%Cells(j)) = eta(MHV_str(i)%Cells(j))
          enddo
        enddo ! end of each cell string
      endif  !H_multi

      !updates of MV (i.e. HUV) bc's are done in q and u extrapolations
      if(WABC .or. windvar .or. windconst) then   !update etas to include wabc
        if(nHstr .gt. 0) then
          do i = 1,nHstr  !for each cell string
            do j=1,H_str(i)%NCells    !for each cell in string
              eta(H_str(i)%Cells(j))  = eta(H_str(i)%Cells(j))  + SWABC(i)%eta(j)
              etan(H_str(i)%Cells(j)) = etan(H_str(i)%Cells(j)) + SWABC(i)%eta(j)
            enddo
          enddo ! end of each cell string
        endif

        if(nTHstr .gt. 0) then     
          i=1
          do j=1,TH_str(i)%NCells    !for each cell in string
            eta(TH_str(i)%Cells(j)) = eta(TH_str(i)%Cells(j)) + TWABC(i)%eta(j)
            etan(TH_str(i)%Cells(j)) = etan(TH_str(i)%Cells(j)) + TWABC(i)%eta(j)
          enddo
        endif

        if(nMHstr .gt. 0) then
          do i = 1,nMHstr  !for each cell string                              !NDRIVER - changed 3/27/2009
            do j=1,MH_str(i)%NCells    !for each cell in string
              eta(MH_str(i)%Cells(j)) = eta(MH_str(i)%Cells(j))  + MWABC(i)%eta(j)
              etan(MH_str(i)%Cells(j))= etan(MH_str(i)%Cells(j)) + MWABC(i)%eta(j)
            enddo
          enddo ! end of each cell string
        endif
      endif
      
      !NEED TO ADD CASE OF nMHVstr > 0 
      
      end subroutine
