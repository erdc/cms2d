!*************************************************************
      subroutine update_variables_tel()
!*************************************************************
      use EXP_Global_def,  only: etan, linktodummiestel, num_linktodummies 
      use EXP_TELESCOPING, only: numxfaces, numyfaces, xface_q, yface_q, xface_qn, yface_qn
      use flow_def, only: eta
      use size_def, only: ncells

      implicit none
      integer i,id1,id2
      
!$OMP PARALLEL 
!$omp do 
      do i=1,ncells
        eta(i) = etan(i)
      enddo
!$OMP END DO

!$OMP DO 
      do i=1,numxFaces
        !if(.not. xface_wall) 
        xface_q(i) = xface_qn(i)
      enddo
!$OMP END DO

!$OMP DO 
      do i=1,numyFaces
        !if(.not. yface(i).wall) 
        yface_q(i) = yface_qn(i)
      enddo
!$OMP END DO
!$OMP END PARALLEL

      !assign wave values to dummy cells
      !copy diff values to dummy cells
      do i=1,num_linktodummies
        id1 = linktodummiesTel(1,i) 
        id2 = linktodummiesTel(2,i)              
        eta(id2) = eta(id1)
      enddo
      
      return
      end subroutine
