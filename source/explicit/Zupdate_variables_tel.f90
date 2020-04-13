      subroutine update_variables_tel()
      use EXP_Global_def 
      USE EXP_bndcond_def
      USE EXP_transport_def 
      use sed_def
      use flow_def
      use comvarbl
      use size_def
      use EXP_TELESCOPING

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
      
      end subroutine
