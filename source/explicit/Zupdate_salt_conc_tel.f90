!*************************************************************
    subroutine update_salt_conc_tel()
!*************************************************************
      use EXP_Global_def,    only: active, etan
      USE EXP_transport_def, only: voln, salt, tsalt_elapse
      use EXP_TELESCOPING,   only: numxfaces, numyfaces, xtransq, ytransq, xface_flux, yface_flux, xface_length, yface_length, xface_wall, yface_wall, cellfaces
      use EXP_TELESCOPING,   only: xface_cells, yface_cells
      use flow_def, only: eta
      use geo_def,  only: dx,dy,zb,cell2cell,mapid
      use prec_def, only: ikind
      use sal_def,  only: sal, sal1
      use size_def, only: ncells, ncellsd
      
      implicit none
      integer i,nm,np
      real(ikind) difft,dxt,dyt,area
      real vdiff,volchange
      
!$OMP PARALLEL 
!$OMP DO 
      do i=1,numxfaces
        if(xTransQ(i) .gt. 0) then
          xface_flux(i) = xTransQ(i)*sal(xface_cells(2,i))*xface_length(i)
        else
          xface_flux(i) = xTransQ(i)*sal(xface_cells(1,i))*xface_length(i)
        endif
      enddo
!$OMP END DO

!$OMP DO 
      do i=1,numyfaces
        if(yTransQ(i) .gt. 0) then
          yface_flux(i) = yTransQ(i)*sal(yface_cells(2,i))*yface_length(i)
        else
          yface_flux(i) = yTransQ(i)*sal(yface_cells(1,i))*yface_length(i)
        endif
      enddo
!$OMP END DO
      
!$OMP DO PRIVATE (Np,Nm,DIFFT,DXT,DYT,AREA)
      do i=1,numxfaces
        if( .not. xface_wall(i)) then    
          np = xface_cells(2,i)
          nm = xface_cells(1,i)          
          diffT = 0.0 !min(salt(np)%diffC,salt(nm)%diffC) !/2.
          dxT = (dx(nm)+dx(np))/2.
          dyT = xface_length(i)
          area = dyT*(eta(nm)-zb(nm)+eta(np)-zb(np))/2.
          !!Fuu(i) = Fuu(i)-area*diffT*(sal(np)-sal(nm))/dxT
        endif
      enddo
!$OMP END DO        
        
!$OMP DO PRIVATE (Np,Nm,DIFFT,DXT,DYT,AREA)
      do i=1,numyfaces
        if( .not. yface_wall(i)) then    
          np = yface_cells(2,i)
          nm = yface_cells(1,i)          
          diffT = 0.0 !min(salt(np)%diffC,salt(nm)%diffC) !/2.
          dyT = (dy(nm)+dy(np))/2.
          dxT = yface_length(i)
          area = dxT*(eta(nm)-zb(nm)+eta(np)-zb(np))/2.
          !!Gvv(i) = Gvv(i)-area*diffT*(sal(np)-sal(nm))/dyT
        endif
      enddo
!$OMP END DO

        call balancecheck_tel_sal_flux()

!$OMP DO PRIVATE (VOLN)
      do i=1,ncells
        if(active(i,3)) then
          sal1(i) = sal(i)*salt(i)%vol             &
            +( xface_flux(cellfaces(8,i))-xface_flux(cellfaces(3,i))  &
            +  xface_flux(cellfaces(7,i))-xface_flux(cellfaces(4,i))  &
            +  yface_flux(cellfaces(6,i))-yface_flux(cellfaces(1,i))  &
            +  yface_flux(cellfaces(5,i))-yface_flux(cellfaces(2,i)))*tsalt_elapse
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
!$OMP END PARALLEL  

        call balancecheck_tel_sal_mass()



!  open(unit=909,file='check.txt')
! 
!  do i=1,numyfaces
!      Gvv(i) = yTransQ(i)*yface(i).length
!  enddo
!  do i=1,numxfaces
!      Fuu(i) = xTransQ(i)*xface_length
!  enddo
!
!   do i=1,ncells
!    if(active(i,3)) then
!     volchange = 
!. +( Fuu(cellfaces(8,i))-Fuu(cellfaces(3,i))
!. +  Fuu(cellfaces(7,i))-Fuu(cellfaces(4,i))
!. +  Gvv(cellfaces(6,i))-Gvv(cellfaces(1,i))                        
!. +  Gvv(cellfaces(5,i))-Gvv(cellfaces(2,i))  
!.  )*tsalt_elapse
!      voln = (-zb(i) + etan(i))*dx(i)*dy(i)
!     Vdiff = voln-salt(i).vol
!     write(909,"(2i5,4f23.14)")i,mapid(i),salt(i).vol,voln,
!.     Vdiff,volchange
!     salt(i).vol = voln
!
!   endif
!   enddo 
!   
!   close(909)
!
!   pause

      end subroutine