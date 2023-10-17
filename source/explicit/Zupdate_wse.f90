!*************************************************************
      subroutine update_wse()
!*************************************************************
      use EXP_Global_def,     only: etan, active, nce, ncw, ncs, ncn, qxn, qyn, dt, drydep, linktodummies
      use EXP_Structures_def, only: structures, srm, srm_on
      use flow_def, only: eta, iwet
      use prec_def, only: ikind
      use size_def, only: ncells, ncellsd
      use geo_def, only: dx,dy,zb,cell2cell
      
      implicit none
      !local vars
      integer i,j,ii,jj
      real totdepth

!$omp parallel do private(nce, ncn, totdepth)
      do i=1,ncells
        if(active(i,3)) then
          nce = cell2cell(2,i)
          ncn = cell2cell(1,i)
          etan(i) = eta(i) + ((qxn(i)-qxn(nce))*dy(i)+  &
               (qyn(i)-qyn(ncn))*dx(i))*dt/(dx(i)*dy(i))
        endif
        totdepth = etan(i)-zb(i)
        if(totdepth.le.0.1d0*drydep) etan(i) = 0.1d0*drydep+zb(i)
        iwet(i) = 1
        if(totdepth.le.1.0_ikind*drydep) iwet(i) = 0  !if(totdepth.le.2.0_ikind*drydep) wet(i) = .false.
      enddo
!$omp end parallel do


!do this for rubble mound structures:
      if(structures) then
        if(SRM_on) then   !rouble mounds
          do j=1,SRM%ncells
            I = SRM%cells(j)
            if(active(i,3)) then
              nce = cell2cell(2,i)
              ncn = cell2cell(1,i)           
              etan(i) = eta(i) + (1./SRM%por(j))*((qxn(i)-qxn(nce))*dy(i)+  &
                       (qyn(i)-qyn(ncn))*dx(i))*dt/(dx(i)*dy(i))     
            endif
            totdepth = etan(i)-zb(i)
            if(totdepth.le.0.1d0*drydep) etan(i) = 0.1d0*drydep +zb(i)
            iwet(i) = 1
            if(totdepth.le.1.0_ikind*drydep) iwet(i) = 0  !if(totdepth.le.2.0_ikind*drydep) wet(i) = .false.
          enddo
        endif !rm on
      endif !rm structure

      !copy diff values to dummy cells
      ii=0
      do i=ncells+1,ncellsD
        ii=ii+1
        jj = linktodummies(ii)
        etan(i) = etan(jj)
      enddo
      
      return
      end subroutine