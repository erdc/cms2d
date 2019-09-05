!***********************************************************************
    subroutine flow_eddyvis
!   update eddy viscosity
!   by Weiming Wu, Oct. 2008
!   modified by Alex Sanchez, USACE 2011
!***********************************************************************
    use size_def
    use geo_def
    use flow_def
    use fric_def, only: cfrict,bsvel
    use struct_def
    use comvarbl
    use wave_flowgrid_def
    use cms_def
    use const_def, only: cappa,small
    use interp_def, only: fintp
    use interp_lib, only: interp_scal_cell2face
    use prec_def
    implicit none
    integer :: i,k,kk,kkdf,nck
    real(ikind) :: ss12,stens,cidrym,vissum

!$OMP PARALLEL    
    selectcase(mturbul)
    case(0)
!$OMP DO PRIVATE(i)          
      do i=1,ncells
        vis(i)=cviscon
      enddo
!$OMP END DO

    case(1) ! Falconer method (1980)
      if(noptset<3)then   
!$OMP DO PRIVATE(i)          
        do i=1,ncells
          vis(i)=cviscon+cvisbot*cfrict(i)*uv(i)*h(i)    !cvisbot=0.578
          vis(i)=min(vis(i),cvismax)
        enddo
!$OMP END DO      
      else   
!$OMP DO PRIVATE(i)
        do i=1,ncells
          vis(i)=cviscon+cvisbot*cfrict(i)*uv(i)*h(i)    !cvisbot=0.578
          vis(i)=vis(i)+cviswav*Whgt(i)*Worb(i) !Bottom friction contribution
          vis(i)=vis(i)+cviswavbrk*h(i)*(max(wavediss(i),1.0e-15)/rhow)**0.33333 !Breaking contribution  
          vis(i)=min(vis(i),cvismax)
        enddo
!$OMP END DO
      endif           
          
    case(2) ! Depth-averaged parobolic model
      if(noptset<3)then   
!$OMP DO PRIVATE(i)        
        do i=1,ncells
          vis(i)=cviscon+cvisbot*bsvel(i)*h(i)
          vis(i)=min(vis(i),cvismax)
        enddo
!$OMP END DO
      else   
!$OMP DO PRIVATE(i)            
        do i=1,ncells
          vis(i)=cviscon+cvisbot*bsvel(i)*h(i)  !Bottom current shear contribution
          vis(i)=vis(i)+cviswav*Whgt(i)*Worb(i) !Bottom wave shear contribution
          vis(i)=vis(i)+cviswavbrk*h(i)*(max(wavediss(i),1.0e-15)/rhow)**0.33333 !Breaking contribution  
          vis(i)=min(vis(i),cvismax)
        enddo
!$OMP END DO     
      endif
          
    case(3) !Wu sub-grid model
      if(noptset<3)then
!$OMP DO PRIVATE(i,ss12,stens)
        do i=1,ncells
          if(iwet(i)==1)then
            ss12=0.5*(duy(i)+dvx(i))            
            stens=2.0*(dux(i)*dux(i)+2.0*ss12*ss12+dvy(i)*dvy(i))
            vis(i)=cviscon+sqrt(cvishor2areaavg2*stens+(cvisbot*bsvel(i)*h(i))**2)
            vis(i)=min(vis(i),cvismax)
          else
            vis(i)=1.0e-6
          endif               
        enddo 
!$OMP END DO
      else
!$OMP DO PRIVATE(i,ss12,stens)
        do i=1,ncells
          if(iwet(i)==1)then
            ss12=0.5*(duy(i)+dvx(i))            
            stens=2.0*(dux(i)*dux(i)+2.0*ss12*ss12+dvy(i)*dvy(i))
            vis(i)=cviscon+sqrt(cvishor2areaavg2*stens+(cvisbot*bsvel(i)*h(i))**2)
            vis(i)=vis(i)+cviswav*Whgt(i)*Worb(i) !Bottom friction contribution
            vis(i)=vis(i)+cviswavbrk*h(i)*(max(wavediss(i),1.0e-15)/rhow)**0.33333 !Breaking contribution  
            vis(i)=min(vis(i),cvismax)
          else
            vis(i)=1.0e-6  
          endif               
        enddo 
!$OMP END DO   
      endif
       
    case(4) ! for mixing length model    
      if(noptset<3)then   
!$OMP DO PRIVATE(i,ss12,stens)                  
        do i=1,ncells
          if(iwet(i)==1)then
            ss12=0.5*(duy(i)+dvx(i))                
            stens=2.0*(dux(i)*dux(i)+2.0*ss12*ss12+dvy(i)*dvy(i))
            vis(i)=cviscon+sqrt((cappa*min(diswall(i),cvishor*h(i)))**4*stens  &
                          +(cvisbot*bsvel(i)*h(i))**2)      
            vis(i)=min(vis(i),cvismax)     
          else
            vis(i)=1.0e-6                  
          endif               
        enddo
!$OMP END DO          
      else
!$OMP DO PRIVATE(i,ss12,stens)                
        do i=1,ncells
          if(iwet(i)==1)then
            ss12=0.5*(duy(i)+dvx(i))                
            stens=2.0*(dux(i)*dux(i)+2.0*ss12*ss12+dvy(i)*dvy(i))
            vis(i)=cviscon+sqrt((cappa*min(diswall(i),cvishor*h(i)))**4*stens  &
                           +(cvisbot*bsvel(i)*h(i))**2)   
            vis(i)=vis(i)+cviswav*Whgt(i)*Worb(i) !Bottom wave shear contribution
            vis(i)=vis(i)+cviswavbrk*h(i)*(max(wavediss(i),1.0e-15)/rhow)**0.33333 !Breaking contribution  
            vis(i)=min(vis(i),cvismax)
          else
            vis(i)=1.0e-6  
          endif               
        enddo
!$OMP END DO            
      endif

    case(5) !Subgrid turbulence model - Alex
      if(noptset<3)then  
!$OMP DO PRIVATE(i,stens)        
        do i=1,ncells
          if(iwet(i)==1)then
            stens=sqrt(2.0*(dux(i)*dux(i)+dvy(i)*dvy(i))+(duy(i)+dvx(i))**2)
            vis(i)=cviscon+cvisbot*bsvel(i)*h(i)+cvishor2areaavg*stens
            vis(i)=min(vis(i),cvismax)
          else
            vis(i)=1.0e-6  
          endif               
        enddo 
!$OMP END DO
      else
!$OMP DO PRIVATE(i,stens)        
        do i=1,ncells
          if(iwet(i)==1)then
            stens=sqrt(2.0*(dux(i)*dux(i)+dvy(i)*dvy(i))+(duy(i)+dvx(i))**2)
            vis(i)=cviscon+cvisbot*bsvel(i)*h(i)+cvishor2areaavg*stens !Bottom and horizontal current shear constributions
            vis(i)=vis(i)+cviswav*Whgt(i)*Worb(i) !Bottom wave shear contribution
            vis(i)=vis(i)+cviswavbrk*h(i)*(max(wavediss(i),1.0e-15)/rhow)**0.33333 !Breaking contribution  
            vis(i)=min(vis(i),cvismax)
          else
            vis(i)=1.0e-6  
          endif
        enddo 
!$OMP END DO
      endif
    endselect

!--- Copy to ghost cells at all forcing BC's -----------
    call bndcopy2ghost(vis)
    
!--- Boundaries Near Dry Cells --------------
    if(ncellpoly==0)then !Cartesian mesh
!$OMP DO PRIVATE(i,k,kk,kkdf,nck,cidrym,vissum)    
      do i=1,ncells
        if(iwet(i)==1)then
          do k=1,ncface(i)
            nck=cell2cell(k,i)
            if(iwet(nck)==0 .or. nck>ncells)then
              kkdf=kkface(idirface(k,i))       
              cidrym=0.0; vissum=0.0; 
              do kk=1,ncface(i)
                if(idirface(kk,i)==kkdf)then   
                  cidrym=cidrym+1.0
                  vissum=vissum+vis(cell2cell(kk,i))-cviscon
                endif
              enddo
              vis(i)=cviscon+vissum/(cidrym+small)  
              vis(nck)=vis(i)                
            endif
          enddo        
        endif
      enddo
!$OMP END DO
    else !Polyhedral mesh
!$OMP DO PRIVATE(i,k,nck)    
      do i=1,ncells
        if(iwet(i)==1)then
          do k=1,ncface(i)
            nck=cell2cell(k,i)
            if(iwet(nck)==0 .or. nck>ncells)then
              vis(nck)=vis(i)                
            endif
          enddo        
        endif
      enddo
!$OMP END DO
    endif
    
!--- Corner wet nodes -----------------------
!$OMP DO PRIVATE(i,k,nck,cidrym,vissum)   
    do i=1,ncells
      if(icorner(i)==1)then
        cidrym=0.0; vissum=0.0
        do k=1,ncface(i)
          nck=cell2cell(k,i)
          cidrym=cidrym+iwet(nck)
          vissum=vissum+(vis(nck)-cviscon)*iwet(nck)
        enddo             
        vis(i)=cviscon+vissum/(cidrym+small)
      endif
    enddo
!$OMP END DO    
!$OMP END PARALLEL 
    
    call interp_scal_cell2face(vis,3,visk) !Interpolate viscosity to cell face, zero at wet/dry boundaries
    !call interp_scal_cell2face(vis,2,visk) !Interpolate viscosity to cell face, zero-gradient at wet/dry boundaries
    call struct_vis !Structures

    return
    endsubroutine flow_eddyvis

!***************************************************************
    subroutine distance2wall(diswall)
! Calculates the minimum distance to a wall from the cell center
! Author: Weiming Wu, Clarkson University
! Modified by Alex Sanchez, USACE-CHL
!***************************************************************
    use size_def, only: ncells,ncellsD
    use geo_def, only: x,y,areap
    use flow_def, only: iwet
    use prec_def
    implicit none
    !Input/Output
    real(ikind),intent(out) :: diswall(ncellsD)
    !Internal Variables
    integer :: i,ii
    real(ikind) :: distmin,distemp
    
    do i=1,ncells
      distmin=1.0e6 !Initialize
!For some reason the OMP section below is not behaving properly for the sudden expansion
!test case. In some instances it runs fine but in others it stalls out.
!!$OMP PARALLEL DO PRIVATE(ii,distemp), REDUCTION(MIN:distmin)
      do ii=1,ncellsD
         if(iwet(ii)==0 .and. i/=ii)then
           !distemp=sqrt( (x(i)-x(ii))**2+(y(i)-y(ii))**2 )-0.5*min(dx(ii),dy(ii))  
           distemp=sqrt((x(i)-x(ii))**2+(y(i)-y(ii))**2)-0.5*sqrt(areap(ii)) !Changed for unstructured cells
           distmin=min(distemp,distmin)
         endif
       enddo
!!$OMP END PARALLEL DO
     diswall(i)=distmin
    enddo
    
    return
    endsubroutine distance2wall
          