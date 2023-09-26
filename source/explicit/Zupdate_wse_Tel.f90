      subroutine update_wse_tel()
      use EXP_Global_def
      USE EXP_bndcond_def
      USE EXP_transport_def 
      use EXP_Structures_def
      use flow_def
      use comvarbl
      use sed_def
      use size_def
      use geo_def, only: dx,dy,zb,cell2cell
      use exp_telescoping
      use BalanceCheck_def 
      
      implicit none
      !local vars
      integer i,j,ii,jj,id1,id2
      real totdepth
      logical :: isnankind
                          
      do ii=1,numREGCells
        i=REGCells(ii) 
        !if(i.eq.3281) then
        !  write (*,*) 'i, etan(i)',i, etan(i)
        !  write (*,*) 'xface_qn(cellfaces(7,i))',xface_qn(cellfaces(7,i))
        !  write (*,*) 'xface_qn(cellfaces(3,i))',xface_qn(cellfaces(3,i))
        !  write (*,*) 'cellfaces(7,i)',cellfaces(7,i)
        !  write (*,*) 'cellfaces(3,i)',cellfaces(3,i)
        !  write (*,*) 'dy(i)',dy(i)
        !  write (*,*) 'yface_qn(cellfaces(5,i))',yface_qn(cellfaces(5,i))
        !  write (*,*) 'yface_qn(cellfaces(1,i))',yface_qn(cellfaces(1,i))
        !  write (*,*) 'cellfaces(5,i)',cellfaces(5,i)
        !  write (*,*) 'cellfaces(1,i)',cellfaces(1,i)
        !  write (*,*) 'dx(i)',dx(i)
        !  write (*,*) 'etan(i) becomes NaN'
        !  write (*,*) ' '
        !endif
        if(active(i,3)) then
          etan(i)=(xface_qn(cellfaces(7,i))-xface_qn(cellfaces(3,i)))*dy(i)  &
                 +(yface_qn(cellfaces(5,i))-yface_qn(cellfaces(1,i)))*dx(i)
          etan(i) = eta(i) + etan(i)*dt/(dx(i)*dy(i))            
        endif
        totdepth = etan(i)-zb(i)
        if(totdepth.le.0.1d0*drydep) etan(i) = 0.1d0*drydep+zb(i)
        iwet(i) = 1
        if(totdepth.le.1.0_ikind*drydep) iwet(i) = 0  !if(totdepth.le.2.0_ikind*drydep) wet(i) = .false.
      enddo
!$OMP PARALLEL
                 
!$omp do private(totdepth,i)
      do ii=1,numTBCells
      i=TBCells(ii) 
        if(active(i,3)) then
         etan(i)= -xface_qn(cellfaces(3,i))*xface_Length(cellfaces(3,i)) &
                  -xface_qn(cellfaces(4,i))*xface_Length(cellfaces(4,i)) & 
                  +xface_qn(cellfaces(7,i))*xface_Length(cellfaces(7,i)) &
                  +xface_qn(cellfaces(8,i))*xface_Length(cellfaces(8,i)) &
                  -yface_qn(cellfaces(1,i))*yface_Length(cellfaces(1,i)) &
                  -yface_qn(cellfaces(2,i))*yface_Length(cellfaces(2,i)) &
                  +yface_qn(cellfaces(5,i))*yface_Length(cellfaces(5,i)) &
                  +yface_qn(cellfaces(6,i))*yface_Length(cellfaces(6,i)) 
         etan(i) = eta(i) + etan(i)*dt/(dx(i)*dy(i))    
        endif
        totdepth = etan(i)-zb(i)
        if(totdepth.le.0.1d0*drydep) etan(i) = 0.1d0*drydep+zb(i)
        iwet(i) = 1
        if(totdepth.le.1.0_ikind*drydep) iwet(i) = 0  !if(totdepth.le.2.0_ikind*drydep) wet(i) = .false.
      enddo
!$omp end do

!$OMP END PARALLEL

!do this for rubble mound structures:
      if(structures) then
        if(SRM_on) then   !rouble mounds
          do j=1,SRM%ncells
            I = SRM%cells(j)
         if(active(i,3)) then
         etan(i)= -xface_qn(cellfaces(3,i))*xface_Length(cellfaces(3,i)) &
                  -xface_qn(cellfaces(4,i))*xface_Length(cellfaces(4,i)) &
                  +xface_qn(cellfaces(7,i))*xface_Length(cellfaces(7,i)) &
                  +xface_qn(cellfaces(8,i))*xface_Length(cellfaces(8,i)) &
                  -yface_qn(cellfaces(1,i))*yface_Length(cellfaces(1,i)) &
                  -yface_qn(cellfaces(2,i))*yface_Length(cellfaces(2,i)) &
                  +yface_qn(cellfaces(5,i))*yface_Length(cellfaces(5,i)) &
                  +yface_qn(cellfaces(6,i))*yface_Length(cellfaces(6,i))
         etan(i) = (1./SRM%por(j))*etan(i)*dt/(dx(i)*dy(i))              
        endif
        totdepth = etan(i)-zb(i)
        if(totdepth.le.0.1d0*drydep) etan(i) = 0.1d0*drydep +zb(i)
        iwet(i) = 1
        if(totdepth.le.1.0_ikind*drydep) iwet(i) = 0  !if(totdepth.le.2.0_ikind*drydep) wet(i) = .false.
      enddo
      endif !rm on
      endif !rm structure

      !copy diff values to dummy cells
         do i=1,num_linktodummies
         id1 = linktodummiesTel(1,i) 
         id2 = linktodummiesTel(2,i)              
         etan(id2) = etan(id1)
         enddo
      
      end subroutine