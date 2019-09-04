      subroutine update_salinity_bc_tel()
	  use EXP_Global_def 
      USE EXP_bndcond_def
      USE EXP_transport_def 
      use bnd_def
      use sal_def
      use sed_def
      use flow_def
	  use comvarbl, only: timehrs
      use geo_Def, only: dx,dy,cell2cell

      implicit none
      integer i,j
      integer ido,ii,isal
      
      call bndsaleval
      
      do isal=1,nsalstr
      do j=1,sal_str(isal)%ncells
        i=sal_str(isal)%cells(j)
        sal(i)=sal_str(isal)%salbnd
      enddo
      enddo      

      end subroutine