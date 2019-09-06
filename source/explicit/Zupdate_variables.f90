      subroutine update_variables()
    use EXP_Global_def 
      USE EXP_bndcond_def
      USE EXP_transport_def 
      use sed_def
      use flow_def
      use comvarbl
      use size_def

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
      
      end subroutine
