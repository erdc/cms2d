!***********************************************************************
     subroutine diswallEXP
!***********************************************************************
     use EXP_Global_def, only: active     
     use size_def, only: ncells,ncellsD
     use prec_def, only: ikind
     use flow_def, only: diswall, iwet
     use geo_def,  only: dx,dy,x,y
     
     implicit none
     integer i,ii
     real(ikind) distemp

     do i=1,ncells
       diswall(i)=1000000000.0
       do ii=1,ncells
         if(iwet(ii).eq.0.and.i.ne.ii) then
           distemp=sqrt( (x(i)-x(ii))**2+(y(i)-y(ii))**2 )-0.5*min(dx(ii),dy(ii))
           if(diswall(i).gt.distemp) diswall(i)=distemp
         endif
       enddo
            
       do ii=ncells,ncellsD
         if(.not. active(i,3) .and. i.ne.ii) then
           distemp=sqrt( (x(i)-x(ii))**2+(y(i)-y(ii))**2 )-0.5*min(dx(ii),dy(ii))
           if(diswall(i).gt.distemp) diswall(i)=distemp
         endif
       enddo           
     enddo
     
     return
     end subroutine diswallEXP