!*************************************************************
      subroutine update_wetdry_tel()
!*************************************************************
      use EXP_Global_def,  only: drydep
      use exp_telescoping, only: xface_qn, yface_qn, numxfaces, numyfaces, xface_cells, yface_cells
      use flow_def, only: eta 
      use prec_def, only: ikind
      use geo_def,  only: zb, cell2cell
      
      implicit none
      integer i,id1,id2
      real(ikind) depth1,depth2
            
!!! LIMIT FLOW BETWEEN CELLS THAT ARE DRY   

!$omp parallel
!$omp do private(depth1,depth2,id1,id2)
      do i=1,numxfaces
    ! if( .not. xface_wall) then  !look at possible wet dry adjustments
      id1 = xface_cells(1,i)
      id2 = xface_cells(2,i)
        DEPTH1 = -zb(id2)+ETA(id2) 
        DEPTH2 = -zb(id1)+ETA(id1)
        IF(DEPTH1.LE.DRYDEP.and.DEPTH2.LE.DRYDEP)THEN
          xface_qn(i) = 0.D0
        ELSEIF(DEPTH1.LE.DRYDEP.and.DEPTH2.GT.DRYDEP)THEN
          IF(xface_qn(i) .GT. 0.d0) xface_qn(i) = 0.D0  
        ELSEIF(DEPTH1.GT.DRYDEP.and.DEPTH2.LE.DRYDEP)THEN
          IF(xface_qn(i) .LT. 0.D0) xface_qn(i) = 0.D0
        ENDIF
    ! endif
      enddo
!$OMP end do           
 
!$omp do private(depth1,depth2,id1,id2)
      do i=1,numyfaces
        !if( .not. yface(i).wall) then  !look at possible wet dry adjustments
        id1 = yface_cells(1,i)
        id2 = yface_cells(2,i)        
        DEPTH1 = -zb(id2) + ETA(id2)
        DEPTH2 = -zb(id1) + ETA(id1)
        IF(DEPTH1.LE.DRYDEP.and.DEPTH2.LE.DRYDEP)THEN
          yface_qn(i) = 0.D0
        ELSEIF(DEPTH1.LE.DRYDEP.and.DEPTH2.GT.DRYDEP)THEN
          IF(yface_qn(i) .GT. 0.D0) yface_qn(i) = 0.D0
        ELSEIF(DEPTH1.GT.DRYDEP.and.DEPTH2.LE.DRYDEP)THEN
          IF(yface_qn(i) .LT. 0.D0) yface_qn(i) = 0.D0
        ENDIF
  !    endif
      enddo
!$OMP end  do    
!$OMP end parallel 

!!$omp parallel do  
!!$omp+ private(depth1,depth2,NCW,NCS)
!      DO I=1,NCELLS
!        NCS = cell2cell(3,i)
!        NCW = cell2cell(4,i)
!        DEPTH1 = -zb(NCW)+ETA(NCW) 
!        DEPTH2 = -zb(I)+ETA(I)
!        IF(DEPTH1.LE.DRYDEP.and.DEPTH2.LE.DRYDEP)THEN
!          QXN(I) = 0.D0
!        ELSEIF(DEPTH1.LE.DRYDEP.and.DEPTH2.GT.DRYDEP)THEN
!          IF(QXN(I) .GT. 0.d0) QXN(I) = 0.D0  !8/31/06
!        ELSEIF(DEPTH1.GT.DRYDEP.and.DEPTH2.LE.DRYDEP)THEN
!          IF(QXN(I) .LT. 0.D0) QXN(I) = 0.D0
!        ENDIF
!        DEPTH1 = -zb(NCS) + ETA(NCS)
!        DEPTH2 = -zb(I) + ETA(I)
!        IF(DEPTH1.LE.DRYDEP.and.DEPTH2.LE.DRYDEP)THEN
!          QYN(I) = 0.D0
!        ELSEIF(DEPTH1.LE.DRYDEP.and.DEPTH2.GT.DRYDEP)THEN
!          IF(QYN(I) .GT. 0.D0) QYN(I) = 0.D0
!        ELSEIF(DEPTH1.GT.DRYDEP.and.DEPTH2.LE.DRYDEP)THEN
!          IF(QYN(I) .LT. 0.D0) QYN(I) = 0.D0
!        ENDIF
!      ENDDO
!!$OMP end parallel do    

      return
      end subroutine
    
!*************************************************************
      subroutine update_wetdry_tel_pre()
!*************************************************************
      use EXP_Global_def,  only: drydep
      use exp_telescoping, only: numxfaces, numyfaces, xface_wet, yface_wet, xface_cells, yface_cells, xface_qn, yface_qn
      use flow_def, only: eta
      use prec_def, only: ikind
      use geo_def,  only: zb, cell2cell
      
      implicit none
      integer i,id1,id2
      real(ikind) depth1,depth2
            
!!! LIMIT FLOW BETWEEN CELLS THAT ARE DRY   

    xface_wet = .true.

!$omp parallel
!$omp do private(depth1,depth2,id1,id2)
      do i=1,numxfaces
    ! if( .not. xface_wall) then  !look at possible wet dry adjustments
      id1 = xface_cells(1,i)
      id2 = xface_cells(2,i)
        DEPTH1 = -zb(id2)+ETA(id2) 
        DEPTH2 = -zb(id1)+ETA(id1)
        IF(DEPTH1.LE.DRYDEP.and.DEPTH2.LE.DRYDEP)THEN
          xface_qn(i) = 0.D0
          xface_wet(i) = .false.
        ENDIF
    ! endif
      enddo
!$OMP end do        

    yface_wet = .true.
 
!$omp do private(depth1,depth2,id1,id2)
      do i=1,numyfaces
        !if( .not. yface(i).wall) then  !look at possible wet dry adjustments
        id1 = yface_cells(1,i)
        id2 = yface_cells(2,i)        
        DEPTH1 = -zb(id2) + ETA(id2)
        DEPTH2 = -zb(id1) + ETA(id1)
        IF(DEPTH1.LE.DRYDEP.and.DEPTH2.LE.DRYDEP)THEN
          yface_qn(i) = 0.D0
          yface_wet(i) = .false.
        ENDIF
  !    endif
      enddo
!$OMP end  do    
!$OMP end parallel 

!!$omp parallel do  
!!$omp+ private(depth1,depth2,NCW,NCS)
!      DO I=1,NCELLS
!        NCS = cell2cell(3,i)
!        NCW = cell2cell(4,i)
!        DEPTH1 = -zb(NCW)+ETA(NCW) 
!        DEPTH2 = -zb(I)+ETA(I)
!        IF(DEPTH1.LE.DRYDEP.and.DEPTH2.LE.DRYDEP)THEN
!          QXN(I) = 0.D0
!        ELSEIF(DEPTH1.LE.DRYDEP.and.DEPTH2.GT.DRYDEP)THEN
!          IF(QXN(I) .GT. 0.d0) QXN(I) = 0.D0  !8/31/06
!        ELSEIF(DEPTH1.GT.DRYDEP.and.DEPTH2.LE.DRYDEP)THEN
!          IF(QXN(I) .LT. 0.D0) QXN(I) = 0.D0
!        ENDIF
!        DEPTH1 = -zb(NCS) + ETA(NCS)
!        DEPTH2 = -zb(I) + ETA(I)
!        IF(DEPTH1.LE.DRYDEP.and.DEPTH2.LE.DRYDEP)THEN
!          QYN(I) = 0.D0
!        ELSEIF(DEPTH1.LE.DRYDEP.and.DEPTH2.GT.DRYDEP)THEN
!          IF(QYN(I) .GT. 0.D0) QYN(I) = 0.D0
!        ELSEIF(DEPTH1.GT.DRYDEP.and.DEPTH2.LE.DRYDEP)THEN
!          IF(QYN(I) .LT. 0.D0) QYN(I) = 0.D0
!        ENDIF
!      ENDDO
!!$OMP end parallel do    

      return
      end subroutine 