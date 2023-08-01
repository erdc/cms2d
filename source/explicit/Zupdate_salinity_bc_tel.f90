!*************************************************************
      subroutine update_salinity_bc_tel()
!*************************************************************
      use sal_def,  only: nsalstr, sal_str, sal
      use comvarbl, only: timehrs

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