      subroutine ST_bed_changes_3()
	use EXP_Global_def 
      USE EXP_bndcond_def
      USE EXP_transport_def      
      use sed_def
      use flow_def
      use comvarbl  
      use size_def 
      use geo_def, only: dx,dy,zb,cell2cell  
      use bnd_Def, only: nQstr, Q_str
      implicit none
      !local vars
      real slp,slopet,vol
      integer i,j,ii

!$OMP PARALLEL 
!$OMP DO PRIVATE (NCN,NCE,NCS,NCW,SLP,SLOPET,VOL)      
      ! update bed laod rates for bed slope effects
      do i=1,ncells
        ncn = cell2cell(1,i)
        nce = cell2cell(2,i)
        ncs = cell2cell(3,i)	
        ncw = cell2cell(4,i)	
        if(active(i,3)) THEN
        
        if(qsx(i).gt.0) then
          if(nce.le.ncells) then
            SLP = 2*(-zb(nce)+zb(i))/ (dx(nce)+dx(i))	
            SLOPET = 1 + SLP*SLPFAC          
            qsx(i) = qsx(i)*SLOPET 
          endif
        else
          if(ncw.le.ncells) then
            SLP = 2*(-zb(ncw)+zb(i))/(dx(ncw)+dx(i))
            SLOPET = 1 + SLP*SLPFAC  	
           qsx(i) = qsx(i) *SLOPET 
          endif
        endif
        if(qsy(i).gt.0) then
          if(ncn.le.ncells) then
            SLP = 2*(-zb(ncn)+zb(i))/(dy(ncn)+dy(i))
            SLOPET = 1 + SLP*SLPFAC
            qsy(i) = qsy(i)*SLOPET 
          endif
        else
          if(ncs.le.ncells) then
            SLP = 2*(-zb(ncs)+zb(i))/(dy(ncs)+dy(i))
            SLOPET = 1 + SLP*SLPFAC
            qsy(i) = qsy(i)*SLOPET 
          endif
        endif

        ENDIF
              
      enddo
!$OMP END DO      
!$OMP END PARALLEL 

   

! update bed due to bed load trasnport
!$OMP PARALLEL 
!$OMP DO PRIVATE (NCN,NCE,NCS,NCW,SLP,SLOPET,VOL)      
      do i=1,ncells
        ncn = cell2cell(1,i)
        nce = cell2cell(2,i)
        ncs = cell2cell(3,i)	
        ncw = cell2cell(4,i)	
        if(active(i,3)) THEN  ! not a wse cell
        
        if(qsx(i).gt.0) then
          if(nce.le.ncells) then       
            vol = qsx(i)*tsed_elapse
            bed(i) =  bed(i) - vol/dx(i)
            bed(nce)= bed(nce) + vol/dx(nce)
          endif
          !if(.not. active(i,1)) bed(i) = 0.0          
        else
          if(ncw.le.ncells) then
            vol = qsx(i)*tsed_elapse
            bed(i) =  bed(i) + vol/dx(i)
            bed(ncw)= bed(ncw)-vol/dx(ncw)
          endif
         !if(.not. active(cell2cell(i,2),1)) bed(i) = 0.0 
        endif
        if(qsy(i).gt.0) then
          if(ncn.le.ncells) then
            vol = qsy(i)*tsed_elapse
            bed(i) = bed(i) - vol/dy(i)
            bed(ncn)= bed(ncn) + vol/dy(ncn)
          endif
        else
          if(ncs.le.ncells) then
            vol = qsy(i)*tsed_elapse
            bed(i) =  bed(i) + vol/dy(i)
            bed(ncs)= bed(ncs) -vol/dy(ncs)
          endif
        endif
        
        ELSE  !it is a wse cell - only update wse cell bed if trans is out of domain

        if(qsx(i).gt.0) then
          if(nce.le.ncells.and.ncw.le.ncells) then    !interior cell    
            vol = qsx(i)*tsed_elapse
            bed(i) =  bed(i) - vol/dx(i)
            bed(nce) = bed(nce) + vol/dx(nce)
           elseif(nce.gt.ncells.and.ncw.le.ncells) then      !out of the domain   
            vol = qsx(i)*tsed_elapse
            bed(i) =  bed(i) - vol/dx(i)          
          endif
        else
          if(ncw.le.ncells.and.nce.le.ncells) then 	!interior cell  
            vol = qsx(i)*tsed_elapse
            bed(i) =  bed(i) + vol/dx(i)
            bed(ncw) =  bed(ncw) - vol/dx(ncw)  
            elseif(ncw.gt.ncells.and.nce.le.ncells) then  !out of the domain 
            vol = qsx(i)*tsed_elapse
            bed(i) =  bed(i) + vol/dx(i)                     
          endif
        endif
        if(qsy(i).gt.0) then
          if(ncn.le.ncells.and.ncs.le.ncells) then
            vol = qsy(i)*tsed_elapse
            bed(i) = bed(i) - vol/dy(i)
            bed(ncn) = bed(ncn) + vol/dy(ncn)            
          elseif(ncn.gt.ncells.and.ncs.le.ncells) then
            vol = qsy(i)*tsed_elapse
            bed(i) = bed(i) - vol/dy(i) 
          endif        
        else
          if(ncs.le.ncells.and.ncn.le.ncells) then
            vol = qsy(i)*tsed_elapse
            bed(i) =  bed(i) + vol/dy(i)
            bed(ncs) =  bed(ncs) - vol/dy(ncs)  
            elseif(ncs.gt.ncells.and.ncn.le.ncells) then
            vol = qsy(i)*tsed_elapse
            bed(i) =  bed(i) + vol/dy(i)                    
          endif
        endif
        
        ENDIF
                      
      enddo
!$OMP END DO      
!$OMP END PARALLEL 

      !no bed change at cells with flow into the grid
      if(nQstr .gt. 0) then
        do i = 1,nQstr  !for each cell string
          if(QstringEXP(i)%vface) then 
              
           if( QstringEXP(i)%sgn .eq. 1 ) then  !south face 
            do j=1,Q_str(i)%NCells	!for each cell in string
              ii= Q_str(i)%Cells(j)
              if(qy(ii) .gt. 0.0)  bed(ii) = 0.0
            enddo
           else  !north face
            do j=1,Q_str(i)%NCells	!for each cell in string
              II = Q_str(i)%Cells(j)
              ncs = cell2cell(3,Q_str(i)%Cells(j))
              if(qy(ii) .lt. 0.0) bed(ncs) = 0.0
            enddo               
          endif
              
          else
              
          if( QstringEXP(i)%sgn .eq. 1 ) then  !west face 
            do j=1,Q_str(i)%NCells	!for each cell in string
              ii = Q_str(i)%Cells(j)
              if(qx(ii) .gt. 0) bed(ii) = 0.0            
            enddo
          else  !east face
            do j=1,Q_str(i)%NCells	!for each cell in string
              II = Q_str(i)%Cells(j)   
              ncw = cell2cell(4,Q_str(i)%Cells(j))             
              if(qx(ii) .lt. 0) bed(ncw) = 0.0 
            enddo        
          endif  
            
         endif
        enddo ! end of NQdriver
      endif  !Q_single     

! update bed due to erosion/deposition for adeq and cohesive schemes
      if(adeq) then
!$OMP PARALLEL DO 
        do i=1,ncells
          bed(i) = bed(i) + (ADSS(i)%depo - ADSS(i)%eros)*tsed_elapse
        enddo
!$OMP END PARALLEL DO
      endif
      if(cohesive) then
!$OMP PARALLEL DO 
        do i=1,ncells
          bed(i) = bed(i) + (COHES(i)%depo - COHES(i)%eros)*tsed_elapse
        enddo
!$OMP END PARALLEL DO
      endif
      
      end subroutine  
 
 
 