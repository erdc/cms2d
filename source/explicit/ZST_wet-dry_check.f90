!***********************************************************************
      subroutine ST_wet_dry_check()
!***********************************************************************
      use EXP_Global_def,    only: ncn, nce, ncs, ncw
      USE EXP_transport_def, only: qsx, qsy
      use flow_def, only: iwet
      use size_def, only: ncells
      use geo_def, only: cell2cell

      implicit none 
      integer i

      !check to make sure we are not violating wet/dry boundaries
!$omp parallel do private (ncn,nce,ncs,ncw) 
      do i=1,ncells
        ncn = cell2cell(1,i)
        nce = cell2cell(2,i)
        ncs = cell2cell(3,i)    
        ncw = cell2cell(4,i)    
        if(qsx(i).gt.0) then
          if(iwet(NCE) .eq. 0) qsx(i) = 0
        else
          if(iwet(NCW) .eq. 0) qsx(i) = 0
        endif
        if(qsy(i).gt.0) then
          if(iwet(NCN) .eq. 0) qsy(i) = 0
        else
          if(iwet(NCS) .eq. 0) qsy(i) = 0
        endif    
      enddo        
!$omp end parallel do

      return
      end subroutine
