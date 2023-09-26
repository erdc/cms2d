!*******************************************************************************
      subroutine update_ADSS_conc()
!*******************************************************************************
      use EXP_Global_def,    only: etan, ncw, ncs, nce, ncn, fuu, gvv, active
      USE EXP_transport_def, only: tsed_elapse, voln, adss
      use size_def, only: ncells
      use geo_def,  only: zb, cell2cell, dx, dy
      use flow_def, only: eta
      use prec_def, only: ikind

      !USE SYNOPTIC_VARS    !Alex, Sep 23, 2009
      
      implicit none
      !local vars
      integer i
      real(ikind) difft,dxt,dyt,area
      
!$OMP PARALLEL 
      !update C based on erosion and deposition
!$OMP DO      
      do i=1,ncells
        ADSS(i)%conc=ADSS(i)%conc-(ADSS(i)%depo-ADSS(i)%eros) * tsed_elapse/(-zb(i) + etan(i))
        ADSS(i)%qx = ADSS(i)%qx/tsed_elapse
        ADSS(i)%qy = ADSS(i)%qy/tsed_elapse
      enddo
!$OMP END DO

      CALL ADEQ_CONC_BC()

!$OMP DO PRIVATE (NCW,NCS)
      do i=1,ncells
        ncs = cell2cell(3,i)
        ncw = cell2cell(4,i)
        if(ADSS(i)%qx.gt.0) then
          Fuu(i) = ADSS(i)%qx*ADSS(ncw)%Conc*dy(i)
        else
          Fuu(i) = ADSS(i)%qx*ADSS(i)%Conc*dy(i)
        endif
        if(ADSS(i)%qy.gt.0) then
          Gvv(i) = ADSS(i)%qy*ADSs(ncs)%Conc*dx(i)
        else
          Gvv(i) = ADSS(i)%qy*ADSS(i)%Conc*dx(i)
        endif
      enddo
!$OMP END DO      
      
!$OMP DO PRIVATE (NCW,NCS,DIFFT,DXT,DYT,AREA)
      do i=1,ncells
        ncs = cell2cell(3,i)
        ncw = cell2cell(4,i)    
        if(active(i,1)) then    
          diffT = (ADSS(i)%diffC+ADSS(ncw)%diffC)/2.
          dxT = (dx(i)+dx(ncw))/2.
          area = dy(i)*(eta(i)-zb(i)+eta(ncw)-zb(ncw))/2.
          Fuu(i) = Fuu(i) - area*diffT*(ADSS(i)%conc-ADSS(ncw)%conc)/dxT
        endif
        if(active(i,2)) then
          diffT = (ADSS(i)%diffC+ADSS(ncs)%diffC)/2.
          dyT = (dy(i)+dy(ncs))/2.
          area = dx(i)*(eta(i)-zb(i)+eta(ncs)-zb(ncs))/2.
          Gvv(i) = Gvv(i) - area*diffT*(ADSS(i)%conc-ADSS(ncs)%conc)/dyT
        endif
      enddo
!$OMP END DO      

!$OMP DO PRIVATE (NCE,NCN,VOLN)
      do i=1,ncells
        ncn = cell2cell(1,i)
        nce = cell2cell(2,i)
        if(active(i,3)) then
          ADSS(i)%concn = ADSS(i)%conc*ADSS(i)%vol + (Fuu(i)-Fuu(nce) + Gvv(i)-Gvv(ncn))*tsed_elapse
          voln = (-zb(i) + etan(i))*dx(i)*dy(i)
          ADSS(i)%concn = ADSS(i)%concn/voln
          ADSS(i)%vol = voln
        endif
      enddo
!$OMP END DO      
      
      !update variables, set flow average to zero for next sed trans time step
!$OMP DO
      do i=1,ncells
        if(active(i,3)) then
          ADSS(i)%conc = ADSS(i)%concn
        endif
        ADSS(i)%qx = 0
        ADSS(i)%qy = 0
      enddo
!$OMP END DO

!$OMP END PARALLEL      

      end subroutine
      
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! sets concentrations at flow and wse boundaries      
      
!*******************************************************************************
      subroutine ADeq_conc_BC()
!*******************************************************************************
      USE EXP_transport_def, only: adss
      use EXP_Global_def,    only: ncn, ncs, ncw, nce, qx, qy, gvv, fuu
      USE EXP_bndcond_def,   only: qstringexp
      use size_def, only: ncells
      use geo_def,  only: cell2cell, dx, dy
      use bnd_def,  only: nqstr, nhstr, q_str, h_str

      implicit none
      !local vars
      integer i,j,ii
      real quT,qvT

       !set concentrations at river boundaries (extrapolate form interior)
      if(nQstr .gt. 0) then
        do i = 1,nQstr  !for each cell string
          if(QstringEXP(i)%vface) then 
              
           if( QstringEXP(i)%sgn .eq. 1 ) then  !south face 
            do j=1,Q_str(i)%NCells    !for each cell in string
              ii= Q_str(i)%Cells(j)
              ncs = cell2cell(3,ii)
              ADSS(ncs)%conc = ADSS(ii)%conc
            enddo
           else  !north face
            do j=1,Q_str(i)%NCells    !for each cell in string
              II = Q_str(i)%Cells(j)
              ncs = cell2cell(3,ii)
              ADSS(ii)%conc = ADSS(ncs)%conc
            enddo               
          endif
              
          else
              
          if( QstringEXP(i)%sgn .eq. 1 ) then  !west face 
            do j=1,Q_str(i)%NCells    !for each cell in string
              ii = Q_str(i)%Cells(j)
              ncw = cell2cell(4,ii)
              ADSS(ncw)%conc = ADSS(ii)%conc         
            enddo
          else  !east face
            do j=1,Q_str(i)%NCells    !for each cell in string
              II = Q_str(i)%Cells(j)   
              ncw = cell2cell(4,ii)
              ADSS(ii)%conc = ADSS(ncw)%conc
            enddo        
          endif  
            
         endif
        enddo ! end of NQdriver
      endif  !Q_single      
      
      
      
! calculate fluxes at river BC cells, since they may not be included in the 
! calculation in the update_adeq_conc subroutine
        do j=1,NQstr  
          if(QstringEXP(j)%vface) then
            do ii=1,Q_str(j)%NCells    !calculate GVV
              i=Q_str(j)%cells(ii)
              if(ADSS(i)%qy.gt.0) then
                ncs = cell2cell(3,i)
                Gvv(i) = ADSS(i)%qy* ADSS(ncs)%Conc*dx(i)
              else
                Gvv(i) = ADSS(i)%qy* ADSS(i)%conc*dx(i)
              endif    
            enddo
          else 
            do ii=1,Q_str(j)%NCells    !calculate FUU
              i=Q_str(j)%cells(ii)
              if(ADSS(i)%qx.gt.0) then
                ncw = cell2cell(4,i)
                Fuu(i) = ADSS(i)%qx*ADSS(ncw)%Conc*dy(i)
              else
                Fuu(i) = ADSS(i)%qx*ADSS(i)%Conc*dy(i)
              endif
            enddo
          endif
        enddo ! end of NQstr    

       if(nHstr .gt. 0) then
        do i = 1,nHstr  !for each cell string
          do j=1,H_str(i)%NCells    !for each cell in string
            ii=H_str(i)%Cells(j) 
            ncn = cell2cell(1,ii)
            nce = cell2cell(2,ii)
            ncs = cell2cell(3,ii)    
            ncw = cell2cell(4,ii)                
            quT = qx(ii)+qx(nce)
            qvT = qy(ii)+qy(ncn)

            if(quT .gt. 0.0 .and. ncw .gt. Ncells) then  !inflow
                ADSS(ii)%conc = ADSS(nce)%conc 
            else  !outflow
                !ADSS(ii)%conc = ADSS(nce)%conc               
            endif
            if(quT .le. 0.0 .and. nce .gt. Ncells) then
                ADSS(ii)%conc = ADSS(ncw)%conc  
            else
                !ADSS(ii)%conc = ADSS(ncw)%conc                
            endif
            if(qvT .gt. 0.0 .and. ncs .gt. Ncells) then
                ADSS(ii)%conc = ADSS(ncn)%conc  
            else
                !ADSS(ii)%conc = ADSS(ncn)%conc               
            endif
            if(qvT .le. 0.0 .and. ncn .gt. Ncells)then
                ADSS(ii)%conc = ADSS(ncs)%conc  
            else
                !ADSS(ii)%conc = ADSS(ncs)%conc                
            endif
           
          enddo
        enddo ! end of each cell string
      endif  !nHstr strings       

      end subroutine adeq_conc_bc