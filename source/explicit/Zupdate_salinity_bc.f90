      subroutine update_salinity_bc()
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
      

        
!*************************************************************
! update low boundary face salt fluxes since they may not be 
! calculated in the main flux calculation algorithm
!*************************************************************        
        
        do j=1,nQstr  
          if(QstringEXP(j)%vface) then
            IDO = Q_str(j)%NCells
            do ii=1,IDO    
              i=Q_str(j)%cells(ii)           
              if(salt(i)%qy.gt.0) then
                ncs = cell2cell(3,i)
                Gvv(i) = salt(i)%qy*sal(ncs)*dx(i)
              else
                Gvv(i) = salt(i)%qy*sal(i)*dx(i)
              endif    
            enddo
          else
            IDO = Q_str(j)%NCells
            do ii=1,IDO    
              i=Q_str(j)%cells(ii)                
              if(salt(i)%qx.gt.0) then
                ncw = cell2cell(4,i)
                Fuu(i) = salt(i)%qx*sal(ncw)*dy(i)
              else
                Fuu(i) = salt(i)%qx*sal(i)*dy(i)
              endif
            enddo
          endif
        enddo 


      end subroutine