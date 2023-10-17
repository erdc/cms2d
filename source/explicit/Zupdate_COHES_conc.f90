!*******************************************************************************
      subroutine update_COHES_conc
!*******************************************************************************
      use EXP_Global_def,  only: etan, ncw, nce, ncn, ncs, fuu, gvv, active
      USE EXP_bndcond_def, only: qstringexp
      USE EXP_transport_def, only: cohes, tsed_elapse, cohes_flow_bc, voln, chparms
      use bnd_def,  only: nqstr, q_str
      use flow_def, only: eta
      use prec_def, only: ikind
      use size_def, only: ncells, ncellsd
      use geo_def,  only: dx, dy, zb, cell2cell    

      implicit none
      !local vars
      integer i,j,ii
      real(ikind) difft,dxt,dyt,area
      
!$OMP PARALLEL 
       !update C based on erosion and deposition
!$OMP DO
      do i=1,ncells
        COHES(i)%conc=COHES(i)%conc-(COHES(i)%depo-COHES(i)%eros) * tsed_elapse/(-zb(i) + etan(i))
        !this needs to be done beofre conc fluc bc since time-ave flows used
        !to calculate fluxes
        COHES(i)%qx = COHES(i)%qx/tsed_elapse
        COHES(i)%qy = COHES(i)%qy/tsed_elapse
      enddo
!$OMP END DO
        
!$OMP SINGLE
      if(cohes_flow_bc) then
      !also need to update cells with flow bc on north and west faces
      !since there IDs are > ncells  
        do j = 1,nQstr  !for each cell string
          if(QstringEXP(j)%vface) then  !N or S face
            if(Q_str(j)%cells(1).gt.ncells) then !north face and need to update
              do i=1,Q_str(j)%NCells    
                ii=Q_str(j)%cells(i)
                COHES(ii)%qy = COHES(ii)%qy/tsed_elapse
              enddo
            endif
          else  !E or W face
            if(Q_str(j)%cells(1).gt.ncells)    then !east face and need to update
              do i=1,Q_str(j)%NCells
                ii=Q_str(j)%cells(i)     
                COHES(ii)%qx = COHES(ii)%qx/tsed_elapse
              enddo
            endif
          endif
        enddo ! end of NQdriver        
        !update conc bc 
        !include "update_COHES_flow_bc.fi"
        call update_COHES_flow_bc()
      endif  !end cohes_flow_bc
!$OMP END SINGLE

!$OMP DO PRIVATE (NCW,NCS)
      do i=1,ncells
        ncs = cell2cell(3,i)
        ncw = cell2cell(4,i)
        if(COHES(i)%qx.gt.0) then
          Fuu(i) = COHES(i)%qx*COHES(ncw)%Conc*dy(i)
        else
          Fuu(i) = COHES(i)%qx*COHES(i)%Conc*dy(i)
        endif
        if(COHES(i)%qy.gt.0) then
          Gvv(i) = COHES(i)%qy*COHES(ncs)%Conc*dx(i)
        else
          Gvv(i) = COHES(i)%qy*COHES(i)%Conc*dx(i)
        endif
      enddo
!$OMP END DO      
      
!$OMP DO PRIVATE (NCW,NCS,DIFFT,DXT,DYT,AREA)
      do i=1,ncells
        ncs = cell2cell(3,i)
        ncw = cell2cell(4,i)    
        if(active(i,1)) then    
          diffT = (COHES(i)%diffC+COHES(ncw)%diffC)/2.
          dxT = (dx(i)+dx(ncw))/2.
          area = dy(i)*(eta(i)-zb(i)+eta(ncw)-zb(ncw))/2.
          Fuu(i) = Fuu(i)-area*diffT*(COHES(i)%conc-COHES(ncw)%conc)/dxT
        endif
        if(active(i,2)) then
          diffT = (COHES(i)%diffC+COHES(ncs)%diffC)/2.
          dyT = (dy(i)+dy(ncs))/2.
          area = dx(i)*(eta(i)-zb(i)+eta(ncs)-zb(ncs))/2.
          Gvv(i) = Gvv(i)-area*diffT*(COHES(i)%conc-COHES(ncs)%conc)/dyT
        endif
      enddo
!$OMP END DO

!$OMP DO PRIVATE (NCE,NCN,VOLN)
      do i=1,ncells
        nce = cell2cell(2,i)
        ncn = cell2cell(1,i)
        if(active(i,3)) then
          COHES(i)%concn = COHES(i)%conc*COHES(i)%vol + (Fuu(i)-Fuu(nce) + Gvv(i)-Gvv(ncn))*tsed_elapse
          voln = (-zb(i) + etan(i))*dx(i)*dy(i)
          COHES(i)%concn = COHES(i)%concn/voln
          COHES(i)%vol = voln
        endif
      enddo
!$OMP END DO   

!$OMP DO
      do i=1,ncellsD
        COHES(i)%qx = 0
        COHES(i)%qy = 0   
      enddo
!$OMP END DO   
    
      !update variables, set flow average to zero for next sed trans time step
!$OMP DO
      do i=1,ncells
        if(active(i,3)) then
          COHES(i)%conc = COHES(i)%concn
        endif
      enddo
!$OMP END DO      

      !update BCs at WSE forcing cells
!$OMP DO PRIVATE (NCN,NCE,NCS,NCW)
      do i=1,ncells
        if(.not. active(i,3)) then
          COHES(i)%conc = CHparms%wse_bc    
        endif
      enddo
!$OMP END DO      
!$OMP END PARALLEL

      end subroutine