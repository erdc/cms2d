!***********************************************************************
      subroutine ST_morphology()  
!***********************************************************************
      use EXP_Global_def,     only: etan, drydep, linktodummies, adeq
      USE EXP_transport_def,  only: bed, adss, cohesive, cohes
      use EXP_Structures_def, only: structures, srm_on, srm
      use size_def, only: ncells, ncellsd
      use geo_def,  only: zb, dx, dy
      use sed_def,  only: scalemorph,poros
   
      implicit none
      !local vars   
      integer i,j,ii,jj

!if rubble mound then bed needs to be scaled up for rock porosity
      if(structures) then
        if(SRM_on) then
          do j=1,SRM%ncells
            I = SRM%cells(j)
            bed(i) = bed(i)/SRM%por(j)
          enddo          
        endif
      endif      
      
!$OMP PARALLEL DO     
      do i=1,ncells
        zb(i) = zb(i) + bed(i)*scalemorph*POROS                !Changing all /POROS to *POROS - 8/28/2008 meb     
      enddo
!$OMP END PARALLEL DO

!reset bed to zero and adjust for wetting and drying
!$OMP PARALLEL DO 
      do i=1,ncells
        bed(i) = 0.0d0
        if((etan(i)-zb(i)).le.0.1d0*drydep) etan(i) = 0.1*drydep +zb(i)
      enddo    
!$OMP END PARALLEL DO      

! set dummy cell depths to interior cell depths
      ii=0
      do i=ncells+1,ncellsD
        ii=ii+1
        jj = linktodummies(ii)
        zb(i) = zb(jj)
      enddo

! adjsut volumes for adeq and cohesive routines      
      if(adeq) then                                                     
!$OMP PARALLEL DO 
        do i=1,ncells                                            
          adss(i)%vol = (-zb(i)+etan(i))*dx(i)*dy(i)
        enddo
!$OMP END PARALLEL DO        
      endif
      
      if(cohesive) then                                                 
!$OMP PARALLEL DO 
        do i=1,ncells
          COHES(i)%vol = (-zb(i)+etan(i))*dx(i)*dy(i)
        enddo
!$OMP END PARALLEL DO        
      endif     
        
      return
      end subroutine ST_morphology
