      subroutine update_wetdry()
    use EXP_Global_def 
      USE EXP_bndcond_def
      USE EXP_transport_def       
      use sed_def
      use flow_def
      use comvarbl 
      use size_def    
      use geo_def, only: zb,cell2cell
      
      implicit none
      integer i
      real(ikind) depth1,depth2
            

!!! LIMIT FLOW BETWEEN CELLS THAT ARE DRY   
! GWH slight reformatting 07/20/08 ! appears this can be a parallel loop
!$omp parallel do private(depth1,depth2,NCW,NCS)
      DO I=1,NCELLS
        NCS = cell2cell(3,i)
        NCW = cell2cell(4,i)
        DEPTH1 = -zb(NCW)+ETA(NCW) 
        DEPTH2 = -zb(I)+ETA(I)
        IF(DEPTH1.LE.DRYDEP.and.DEPTH2.LE.DRYDEP)THEN
          QXN(I) = 0.D0
        ELSEIF(DEPTH1.LE.DRYDEP.and.DEPTH2.GT.DRYDEP)THEN
          IF(QXN(I) .GT. 0.d0) QXN(I) = 0.D0  !8/31/06
        ELSEIF(DEPTH1.GT.DRYDEP.and.DEPTH2.LE.DRYDEP)THEN
          IF(QXN(I) .LT. 0.D0) QXN(I) = 0.D0
        ENDIF
        DEPTH1 = -zb(NCS) + ETA(NCS)
        DEPTH2 = -zb(I) + ETA(I)
        IF(DEPTH1.LE.DRYDEP.and.DEPTH2.LE.DRYDEP)THEN
          QYN(I) = 0.D0
        ELSEIF(DEPTH1.LE.DRYDEP.and.DEPTH2.GT.DRYDEP)THEN
          IF(QYN(I) .GT. 0.D0) QYN(I) = 0.D0
        ELSEIF(DEPTH1.GT.DRYDEP.and.DEPTH2.LE.DRYDEP)THEN
          IF(QYN(I) .LT. 0.D0) QYN(I) = 0.D0
        ENDIF
      ENDDO
!$OMP end parallel do    

      end subroutine