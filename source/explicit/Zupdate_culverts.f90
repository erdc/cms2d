!*******************************************************************************
      subroutine update_culverts()
!*******************************************************************************
      use EXP_Global_def,     only: etan, dt, adeq
      USE EXP_transport_def,  only: salt, voln, adss
      use EXP_Structures_def, only: structures, cul_on, cul
      use sal_def,   only: saltrans, sal, sal1
      use flow_def,  only: eta, iwet
      use prec_def,  only: ikind
      use geo_def,   only: dx,dy,zb    
      use const_def, only: pi 

      implicit none
      integer j,id1,id2
      real(ikind) xh1,xh2,dh,have,hgt,area,per,rh,flow,arg

      if(structures) then
      if(CUL_on) then
      
      do j=1,cul%num
        ID1 = cul%cells1(j)
        ID2 = cul%cells2(j) 
        xH1 = max(eta(ID1),cul%invert(j))
        xH2 = max(eta(ID2),cul%invert(j))
        DH = xH1 - xH2     
        Have = (xH1+xH2)/2.0
        hgt = cul%invert(j) + cul%dia(j)
        if(Have .le. cul%invert(j)) then  !average elev below invert - no flow
          cul%flow(j) = 0.0  
        else     !average elev above invert - calculate flow
          if(have .gt. hgt) then  
            Area = pi*(cul%dia(j)**2)/4.0
            Per = pi*cul%dia(j)
            RH = Area/Per
            Flow = Area*sqrt(abs(DH)/( 0.9/9.8 + 6.25*(cul%fric(j)**2)*cul%length(j)/(RH**1.33) ) )
            cul%flow(j) = sign(flow,DH)
          else
            arg = (Have - cul%invert(j))/cul%dia(j)        
            Area = pi*(cul%dia(j)**2)/4.0
            Area = Area*(-1.15*arg**3 + 1.72*arg**2 + 0.43*arg)
            Per = pi*cul%dia(j)
            Per = Per*(2.15*arg**3 - 3.30*arg**2 + 2.15*arg)
            RH = Area/Per
            Flow = Area*sqrt(abs(DH)/( 0.9/9.8 + 6.25*(cul%fric(j)**2)*cul%length(j)/(RH**1.33) ) )
            cul%flow(j) = sign(flow,DH)              
          endif
        endif
      enddo 
      
      !wet dry check 
      do j=1,cul%num
        ID1 = cul%cells1(j)
        ID2 = cul%cells2(j)      
        if(cul%flow(j).gt.0 .and. iwet(ID1) .eq. 0 ) cul%flow(j) = 0.0
        if(cul%flow(j).lt.0 .and. iwet(ID2) .eq. 0) cul%flow(j) = 0.0      
      enddo      
      
      !update etan
      do j=1,cul%num
        ID1 = cul%cells1(j)
        ID2 = cul%cells2(j)      
        !if(cul%flow(j).gt.0) then
        area = dx(id1)*dy(id1)
        etan(ID1) = etan(ID1) - cul%flow(j)*dt/area
        area = dx(id2)*dy(id2)
        etan(ID2) = etan(ID2) + cul%flow(j)*dt/area     
        !endif
      enddo
      
      if(saltrans) then
        do j=1,cul%num
          ID1 = cul%cells1(j)       
          if(cul%flow(j).ge.0) cul%saltF(j) = cul%flow(j)*sal(id1)
          ID2 = cul%cells2(j)       
          if(cul%flow(j).lt.0) cul%saltF(j) = cul%flow(j)*sal(id2)            
        enddo
       
        do j=1,cul%num
          ID1 = cul%cells1(j)       
          sal1(id1) = sal(id1)*salt(id1)%vol
          ID2 = cul%cells2(j)       
          sal1(id2) = sal(id2)*salt(id2)%vol           
        enddo            
      
        do j=1,cul%num
          ID1 = cul%cells1(j)
          ID2 = cul%cells2(j)       
          sal1(id1) = sal1(id1) - cul%saltF(j)*dt
          sal1(id2) = sal1(id2) + cul%saltF(j)*dt     
        enddo
       
        do j=1,cul%numui
          ID1 = cul%cellUI(j)
          voln = (-zb(id1) + etan(id1))*dx(id1)*dy(id1)
          sal1(id1) = sal1(id1)/voln
          salt(id1)%vol = voln     
        enddo 
       
        do j=1,cul%numui
          ID1 = cul%cellUI(j)       
          sal(id1) = sal1(id1)
        enddo 
      endif !end saltrans
       
      if(adeq) then  !(using saltF to track transport of adeq)
        do j=1,cul%num
          ID1 = cul%cells1(j)       
          if(cul%flow(j).ge.0) cul%saltF(j) = cul%flow(j)*adss(id1)%conc
          ID2 = cul%cells2(j)       
          if(cul%flow(j).lt.0) cul%saltF(j) = cul%flow(j)*adss(id2)%conc            
        enddo
       
        do j=1,cul%num
          ID1 = cul%cells1(j)       
          adss(id1)%concn = adss(id1)%conc*adss(id1)%vol
          ID2 = cul%cells2(j)       
          adss(id2)%concn = adss(id2)%conc*adss(id2)%vol           
        enddo            
      
        do j=1,cul%num
          ID1 = cul%cells1(j)
          ID2 = cul%cells2(j)       
          adss(id1)%concn = adss(id1)%concn - cul%saltF(j)*dt
          adss(id2)%concn = adss(id2)%concn + cul%saltF(j)*dt     
        enddo
       
        do j=1,cul%numui
          ID1 = cul%cellUI(j)
          voln = (-zb(id1) + etan(id1))*dx(id1)*dy(id1)
          adss(id1)%concn = adss(id1)%concn/voln
          adss(id1)%vol = voln     
        enddo 
       
        do j=1,cul%numui
          ID1 = cul%cellUI(j)       
          adss(id1)%conc = adss(id1)%concn
        enddo 
      endif !end adeq  (Ad scheme suspended sediment)       
 
      do j=1,cul%numui
        ID1 = cul%cellUI(j)
        eta(id1) = etan(id1)
      enddo

      endif !cul%on
      endif !structures
      
      end subroutine