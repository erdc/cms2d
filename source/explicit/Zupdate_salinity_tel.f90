      subroutine update_salinity_tel()
    use EXP_Global_def 
      USE EXP_bndcond_def
      USE EXP_transport_def       
      use bnd_def
      use sed_def
      use flow_def
      use comvarbl      
      use sal_def
      use size_def
      use EXP_TELESCOPING

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