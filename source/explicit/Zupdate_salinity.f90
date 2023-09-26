!*************************************************************
      subroutine update_salinity()
!*************************************************************
      use EXP_Global_def,    only: qx, qxn, qy, qyn, dt, dtsalt
      USE EXP_bndcond_def,   only: qstringexp
      USE EXP_transport_def, only: salt, tsalt_elapse
      use bnd_def,  only: nqstr, q_str
      use flow_def, only: iwet, vis
      use size_def, only: ncells, ncellsd

      implicit none
      integer i,j,ii
          

!$OMP PARALLEL DO 
        do i = 1,ncells
          salt(i)%qx = salt(i)%qx + qxn(i)*dt
          salt(i)%qy = salt(i)%qy + qyn(i)*dt
        enddo
!$OMP END PARALLEL DO

        if(nQstr .gt. 0) then
          !also need to update cells with flow bc on north and west faces
          !since there IDs are > ncells  
          do j = 1,nQstr  !for each cell string
            if(QstringEXP(j)%vface) then  !N or S face
              if(Q_str(j)%cells(1).gt.ncells) then !north face and need to update
                do i=1,Q_str(j)%NCells    
                  ii=Q_str(j)%cells(i)
                  salt(ii)%qy = salt(ii)%qy + qy(ii)*dt
                enddo
              endif
            else  !E or W face
              if(Q_str(j)%cells(1).gt.ncells)    then !east face and need to update
                do i=1,Q_str(j)%NCells
                  ii=Q_str(j)%cells(i)     
                  salt(ii)%qx = salt(ii)%qx + qx(ii)*dt
                enddo
              endif
            endif
          enddo ! end of NQdriver        
        endif  !end if(Q_single)     
        
        
!!******************************************************
!!  calculate salt update for transport time step       
!!*******************************************************           
      
        tsalt_elapse = tsalt_elapse + dt
        
        if(tsalt_elapse .ge. dtsalt) then 
          
          do i=1,ncellsD
            salt(i)%diffC = 0.0
            if(iwet(i) .eq. 1) salt(i)%diffC = vis(i)    
          enddo    
      
          !this needs to be done before conc fluc bc since time-ave flows used
          !to calculate fluxes
!$OMP PARALLEL DO 
          do i=1,ncells
            salt(i)%qx = salt(i)%qx/tsalt_elapse
            salt(i)%qy = salt(i)%qy/tsalt_elapse
          enddo
!$OMP END PARALLEL DO
      
          if(nQstr .gt. 0) then
            !also need to update cells with flow bc on north and west faces
            !since there IDs are > ncells  
            do j = 1,nQstr  !for each cell string
              if(QstringEXP(j)%vface) then  !N or S face
                if(Q_str(j)%cells(1).gt.ncells) then !north face and need to update
                  do i=1,Q_str(j)%NCells    
                    ii=Q_str(j)%cells(i)
                    salt(ii)%qy = salt(ii)%qy/tsalt_elapse
                  enddo
                endif
              else  !E or W face
                if(Q_str(j)%cells(1).gt.ncells) then !east face and need to update
                  do i=1,Q_str(j)%NCells
                    ii=Q_str(j)%cells(i)     
                    salt(ii)%qx = salt(ii)%qx/tsalt_elapse
                  enddo
                endif
              endif
            enddo ! end of NQdriver    
          endif  !end if(Q_single)
      
          call update_salinity_bc()
    
          call update_salt_conc()
          
          tsalt_elapse = 0.0
          
        endif !end "if tsed > dtsed"
      
      end subroutine