!*************************************************************
      subroutine update_salinity_tel()
!*************************************************************
      use EXP_Global_def,    only: dt, dtsalt
      USE EXP_transport_def, only: tsalt_elapse, salt
      use EXP_TELESCOPING,   only: numxfaces, numyfaces, xtransq, ytransq, xface_qn, yface_qn
      use flow_def, only: iwet, vis
      use size_def, only: ncellsd

      implicit none
      integer i,j,ii
          
!$OMP PARALLEL DO 
        do i = 1,numxfaces
          xTransQ(i) = xTransQ(i) + xface_qn(i)*dt
        enddo
!$OMP END PARALLEL DO

!$OMP PARALLEL DO 
        do i = 1,numyfaces
          yTransQ(i) = yTransQ(i) + yface_qn(i)*dt
        enddo
!$OMP END PARALLEL DO
        
        
!!******************************************************
!!  calculate salt update for transport time step       
!!*******************************************************           
      
        tsalt_elapse = tsalt_elapse + dt
        
        if(tsalt_elapse .ge. dtsalt) then 
          
          do i=1,ncellsD
            salt(i)%diffC = 0.0
            if(iwet(i) .eq. 1) salt(i)%diffC = vis(i)    
          enddo    
      
          !this needs to be done before conc fluc bc 
          !since time-ave flows used to calculate fluxes
!$OMP PARALLEL DO 
          do i=1,numxfaces
            xTransQ(i) = xTransQ(i)/tsalt_elapse
          enddo
!$OMP END PARALLEL DO

!$OMP PARALLEL DO 
          do i=1,numyfaces
            yTransQ(i) = yTransQ(i)/tsalt_elapse
          enddo
!$OMP END PARALLEL DO
      
          call update_salinity_bc_tel()
    
          call update_salt_conc_tel()
          
          xTransQ=0.0
          yTransQ=0.0
          tsalt_elapse = 0.0
          
        endif !end "if tsed > dtsed"
        
      
      end subroutine