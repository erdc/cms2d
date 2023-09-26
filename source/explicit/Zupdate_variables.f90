!*************************************************************
      subroutine update_variables()
!*************************************************************
      use EXP_Global_def, only: qx, qy, qxn, qyn, linktodummies, etan
      use flow_def, only: eta
      use size_def, only: ncells, ncellsd

      implicit none
      integer i,ii,jj
      
!$OMP PARALLEL DO 
      do i=1,ncells
        qx(i) = qxn(i)
        qy(i) = qyn(i)
        eta(i) = etan(i)
      enddo
!$OMP END PARALLEL DO

      !assign wave values to dummy cells
      ii=0
      do i=ncells+1,ncellsD
        ii=ii+1
        jj = linktodummies(ii)
        eta(i) = eta(jj)
      enddo
      
      return
      end subroutine
