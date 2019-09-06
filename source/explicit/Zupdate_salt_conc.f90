    subroutine update_salt_conc()
      use EXP_Global_def
      USE EXP_transport_def     
      USE EXP_bndcond_def
      use flow_def
      use comvarbl
      use geo_def, only: dx,dy,zb,cell2cell
      use sed_def
      use sal_def 
      use size_def       
      
      implicit none
      integer i
      real(ikind) difft,dxt,dyt,area
      
      
!$OMP PARALLEL 
!$OMP DO PRIVATE (NCW,NCS)
      do i=1,ncells
        ncs = cell2cell(3,i)
        ncw = cell2cell(4,i)
        if(salt(i)%qx.gt.0) then
          Fuu(i) = salt(i)%qx*sal(ncw)*dy(i)
        else
          Fuu(i) = salt(i)%qx*sal(i)*dy(i)
        endif
        if(salt(i)%qy.gt.0) then
          Gvv(i) = salt(i)%qy*sal(ncs)*dx(i)
        else
          Gvv(i) = salt(i)%qy*sal(i)*dx(i)
        endif
      enddo
!$OMP END DO
      
!$OMP DO PRIVATE (NCW,NCS,DIFFT,DXT,DYT,AREA)
      do i=1,ncells
        ncs = cell2cell(3,i)
        ncw = cell2cell(4,i)    
        if(active(i,1)) then    
          diffT = min(salt(i)%diffC,salt(ncw)%diffC) !/2.
          dxT = (dx(i)+dx(ncw))/2.
          area = dy(i)*(eta(i)-zb(i)+eta(ncw)-zb(ncw))/2.
          Fuu(i) = Fuu(i)-area*diffT*(sal(i)-sal(ncw))/dxT
        endif
        if(active(i,2)) then
          diffT = min(salt(i)%diffC,salt(ncs)%diffC) !/2.
          dyT = (dy(i)+dy(ncs))/2.
          area = dx(i)*(eta(i)-zb(i)+eta(ncs)-zb(ncs))/2.
          Gvv(i) = Gvv(i)-area*diffT*(sal(i)-sal(ncs))/dyT
        endif
      enddo
!$OMP END DO

!$OMP DO PRIVATE (NCE,NCN,VOLN)
      do i=1,ncells
        ncn = cell2cell(1,i)
        nce = cell2cell(2,i)
        if(active(i,3)) then
          sal1(i) = sal(i)*salt(i)%vol +(Fuu(i)-Fuu(nce) + Gvv(i)-Gvv(ncn))*tsalt_elapse
          voln = (-zb(i) + etan(i))*dx(i)*dy(i)
          sal1(i) = sal1(i)/voln
          salt(i)%vol = voln
        endif
      enddo
!$OMP END DO      

!update variables, set flow average to zero for next sed trans time step
!$OMP DO
      do i=1,ncells
        if(active(i,3)) sal(i) = sal1(i)
      enddo
!$OMP END DO      

!$OMP DO
      do i=1,ncellsD
        salt(i)%qx = 0
        salt(i)%qy = 0
      enddo  
!$OMP END DO      
!$OMP END PARALLEL        

      end subroutine