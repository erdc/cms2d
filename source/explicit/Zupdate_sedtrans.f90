      subroutine update_sedtrans()
      use EXP_Global_def 
      USE EXP_bndcond_def
      USE EXP_transport_def      
      use bnd_def
      use sed_def
      use flow_def
      use comvarbl 
      use size_def    

      implicit none
      !local vars
      integer :: i,j,ii

      if(sedtransEXP) then
        !if AD, then accumulate flow to get averages when updated SS conc    
        if(adeq) then                                                   
!$omp parallel do                
          do i = 1,ncells
            ADSS(i)%qx = ADSS(i)%qx + qxn(i)*dt
            ADSS(i)%qy = ADSS(i)%qy + qyn(i)*dt
          enddo
!$omp end parallel do            
        endif    
        
        !if COHES, then accumulate flow to get averages when updated SS conc    
        if(cohesive) then                                               
!$omp parallel do                
          do i = 1,ncells
            COHES(i)%qx = COHES(i)%qx + qxn(i)*dt
            COHES(i)%qy = COHES(i)%qy + qyn(i)*dt
          enddo
!$omp end parallel do            
        endif    
        
        if(cohes_flow_bc) then
          !also need to update cells with flow bc on north and west faces
          !since there IDs are > ncells  
          do j = 1,NQstr  !for each cell string
            if(QstringEXP(j)%vface) then  !N or S face
              if(Q_Str(j)%cells(1).gt.ncells) then      !north face and need to update
                do i=1,Q_Str(j)%NCells    
                  ii=Q_Str(j)%cells(i)
                  COHES(ii)%qy = COHES(ii)%qy + qyn(ii)*dt
                enddo
              endif
            else  !E or W face
              if(Q_Str(j)%cells(1).gt.ncells) then      !east face and need to update
                do i=1,Q_Str(j)%NCells
                  ii=Q_Str(j)%cells(i)     
                  COHES(ii)%qx = COHES(ii)%qx + qxn(ii)*dt
                enddo
              endif
            endif
          enddo ! end of NQdriver        
        endif  !end cohes_flow_bc    
      
        tsed_elapse = tsed_elapse + dt
        
        if(tsed_elapse .ge. dtsed ) then
          if(watanabe) then
            call ST_watanabe()
            call ST_wet_dry_check()
            if (do_aval)  call ST_avalanche()    
            call ST_hardbottom()
            call TL_prep_for_output 
          elseif(lundcirp) then
            call ST_lundcirp()
            call ST_wet_dry_check()
            if (do_aval) call ST_avalanche()    
            call ST_hardbottom() 
            call TL_prep_for_output            
          elseif(adeq) then
            call ST_adeq()
            call ST_wet_dry_check()
            if (do_aval) call ST_avalanche()
            call ST_hardbottom_AD()     
            call update_ADSS_conc()  
            call ADSS_prep_for_output          
          elseif(cohesive) then
            call ST_cohesive()
            call ST_hardbottom_COHES()     
            call update_COHES_conc()        
          endif
          !calculate bed chagnes
          call ST_bed_changes_3()
          tsed_elapse = 0.0    
        endif

        tmorph_elapse = tmorph_elapse + dt
        if(tmorph_elapse .ge. dtmorph ) then
          call ST_morphology()
          tmorph_elapse = 0.0
        endif
      endif
      
      end subroutine update_sedtrans
      
      !*************************************************************
      !Maps explicit variable to implicit variable for output for
      ! Total Load sediment transport routines
      !  Chris Reed   6/28/2013
      !*************************************************************      
      subroutine ADSS_prep_for_output
      use EXP_Global_def 
      USE EXP_bndcond_def
      USE EXP_transport_def      
      use bnd_def
      use sed_def
      use flow_def
      use const_def, only: small
      use size_def 
      use out_def, only: write_fracsusp
      use geo_def, only: cell2cell
      
      implicit none
      integer i
      real(ikind) Q,qsusX,QsusY,qbedX,qbedY,Qsus,Qtot
      
      !map to output variable for concentration
      ct = ADSS%conc*rhosed
      
      !map total transport and rs (if requested)
      if(write_fracsusp)then
        do i=1,ncells
          nce = cell2cell(2,i)
          Q = (qx(i)+qx(nce))/2.
          qsusX = ct(i)*Q 
          qbedX = qsx(i)*rhosed
          ncn = cell2cell(1,i)
          Q = (qy(i)+qy(ncn))/2.
          qsusY = ct(i)*Q 
          qbedY = qsy(i)*rhosed
          qtx(i) = qsusX + qbedX
          qty(i) = qsusY + qbedY
          Qsus = sqrt(qsusX**2+qsusY**2) 
          Qtot = sqrt(qtx(i)**2+qty(i)**2)
          rs(i) = Qsus / max(Qtot,small)
        enddo
      else  ! only mapp total load
        do i=1,ncells
          nce = cell2cell(2,i)
          Q = (qx(i)+qx(nce))/2.
          qtx(i) = ct(i)*Q + qsx(i)*rhosed
          Q = (qy(i)+qy(nce))/2.      
          qty(i) = ct(i)*Q + qsx(i)*rhosed
        enddo
      endif
      
      end subroutine ADSS_prep_for_output
      
      !*************************************************************
      !Maps explicit variable to implicit variable for output for
      ! Total Load sediment transport routines
      !  Chris Reed   6/28/2013
      !*************************************************************
      subroutine TL_prep_for_output
      USE EXP_transport_def      
      use sed_def

      implicit none
      
      qtx = qsx
      qty = qsy
      
      end subroutine TL_prep_for_output
      