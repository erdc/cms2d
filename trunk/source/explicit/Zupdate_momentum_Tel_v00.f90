     subroutine update_momentum_tel()
      use EXP_Global_def
      USE EXP_bndcond_def	
      use EXP_Structures_def
	  use flow_def 
	  use geo_def
      use sed_def
      use const_def, only: pi,deg2rad    
      use met_def, only: tauwindx,tauwindy,pressatm
      use wave_flowgrid_def, only: wavestrx,wavestry     
      use fric_def, only: cfrict,uelwc,coefman
      USE EXP_transport_def 
      use NupdateMod
      use sal_def  
      use exp_telescoping
      
      implicit none
      integer i,j,id,id1,id2,id3,id4,ii,advf3,advf4
      real(ikind) dyt,dxt,difft,deptht,hplus,HmnEXP,hgt,coriolis
      real(ikind) spd,deltax,deltay,tauwind
      real(ikind) detadx,detady,cd,depth,FluxEXP
      real(ikind) advdif_i,advdif_c
      real(ikind) Nmann

!!$omp parallel 
        
      ! update density differential due to salinity  
      !if(saltrans .and. saltsimD) then
!!$omp do       
     !  do i=1,ncells
         !RHOPrim(i) = 1.0 + 0.000808*salt(i).Conc
     !  enddo
!!$omp end do       
    !   endif   
      
      
!!$omp do private(i,id1,id2,hplus,HmnEXP,hgt,Coriolis,spd,cd,deltax,detadx,tauwind,depth,advdif_I,advdif_C,advf3,advf4)
      do ii=1,numTBXfaces
        i=TBXfaces(ii)
        if(.not. xface_wall(i)) then
          id1 = xface_cells(1,i)
          id2 = xface_cells(2,i)
          id3 = xface_cells(3,i)
          id4 = xface_cells(4,i)
          
          DELTAX = (DX(Id1)+DX(id2))/2.0_ikind
          if(xface_basic_orientation(i)) then
            HPLUS  = (-zb(id1)) +eta(id1)  
            HmnEXP = (-zb(id2)) +eta(id2)*xface_gcoef(2,i) + eta(id3)*xface_gcoef(3,i) + eta(id4)*xface_gcoef(4,i)
          else
            HPLUS  = (-zb(id1)) +eta(id1)*xface_gcoef(1,i) + eta(id3)*xface_gcoef(3,i) + eta(id4)*xface_gcoef(4,i)   
            HmnEXP = (-zb(id2)) +eta(id2)                         
          endif
          DETADX = (( (RHOPrim(id1)*HPLUS +pressatm(id1))*HPLUS - &
                   (RHOPrim(id2)*HmnEXP +pressatm(id2))*HmnEXP)- &
                   ((RHOPrim(id1)+RHOPrim(id2))/2)* &
                   ((pressatm(id1)+pressatm(id2))/2. + HPLUS + HmnEXP)* &
                   ( (-zb(Id1))- (-zb(id2))))/DELTAX   
             
          ! inline advection and diffusion 
          if(xface_advF(2,i) .gt. 0) then
            ADVDIF_I = (xface_advdif_i(1,i)+xface_advdif_i(2,i))/2 
          else
            ADVDIF_I = xface_advdif_i(1,i)
          endif
          advf3 =  xface_advF(3,i)
          advf4 =  xface_advF(4,i)         
          if(advf4 .gt. 0) then
            ADVDIF_I = ADVDIF_I - (xface_advdif_i(1,advf3) + xface_advdif_i(1,advf4))/2.
          else
            ADVDIF_I = ADVDIF_I - xface_advdif_i(xface_side(ii),advf3)
          endif
             
          
          ADVDIF_I = ADVDIF_I/DELTAX 
         
          !cross term advectionand diffusion       
          ADVDIF_C = xface_advdif_C(i) - (xface_advdif_C(xface_CadvF(6,i)) +xface_advdif_C(xface_CadvF(4,i)))/2.0 
          ADVDIF_C = ADVDIF_C/xface_length(i)
    
          tauwind  = 0.5*(tauwindx(id1)+tauwindx(id2))          
          CORIOLIS = 0.0  !fcoriolis*QYc(i)*0.0
          
          !HGT = (HPLUS+HmnEXP)/2.0_ikind 
          !spd = (uelwc(id1)+uelwc(id2))/2.0_ikind
          !cd =  (cfrict(id1)+Cfrict(id2))/2.0_ikind          
         
          !HGT = (HPLUS*dx(id2)+HmnEXP*dx(id1))/(dx(id1)+dx(id2)) 
          !spd = (uelwc(id1)*dx(id2)+uelwc(id2)*dx(id1))/(dx(id1)+dx(id2)) 
          !cd =  (cfrict(id1)*dx(id2)+Cfrict(id2)*dx(id1))/(dx(id1)+dx(id2)) 
          
          !!Option A
          HGT = (HPLUS+HmnEXP)/2.0_ikind 
          spd = (uelwc(id1)*xface_gcoef(1,i)+uelwc(id2)*xface_gcoef(2,i)+uelwc(id3)*xface_gcoef(3,i)+uelwc(id4)*xface_gcoef(4,i))/2.0_ikind        
          cd =  (cfrict(id1)*xface_gcoef(1,i)+Cfrict(id2)*xface_gcoef(2,i)+cfrict(id3)*xface_gcoef(3,i)+cfrict(id4)*xface_gcoef(4,i))/2.0_ikind           

          !!Option B
          !HGT = (HPLUS+HmnEXP)/2.0_ikind 
          !spd = abs(xface_vel(i))
          !Nmann = (coefman(id1)*dy(id2)+coefman(id2)*dy(id1))/(dY(id1)+dy(ID2))
          !cd = grav*Nmann*Nmann/HGT**0.333      
              
          xface_QN(i) = (xface_Q(i) + DT*((-GRAV)*DETADX/2.0_ikind &
                     +  0.5*(wavestrx(id1) + wavestrx(id2)) - ADVDIF_I - ADVDIF_C &
                     + tauwind + CORIOLIS) )/(1+dt*Cd*spd/HGT)
                 
          !xface_QN(i) = xface_Q(i) + DT*(- ADVDIF_I - ADVDIF_C )     
        
        endif  ! end "not a wall"
      enddo
!!$omp end do    


!!$omp do private(i,id1,id2,hplus,HmnEXP,hgt,Coriolis,spd,cd,deltay,detady,tauwind,depth,advdif_i,advdif_c,advf3)
      do ii=1,numTBYfaces
        i=TBYfaces(ii)
        if(.not. yface_wall(i)) then
          id1 = yface_cells(1,i)
          id2 = yface_cells(2,i)
          id3 = yface_cells(3,i)
          id4 = yface_cells(4,i)         
 
          DELTAY = (DY(Id1)+DY(id2))/2.0_ikind
          if(yface_basic_orientation(i)) then
            HPLUS  = (-zb(id1)) +eta(id1)  
            HmnEXP = (-zb(id2)) +eta(id2)*yface_gcoef(2,i) + eta(id3)*yface_gcoef(3,i) + eta(id4)*yface_gcoef(4,i)
          else
            HPLUS  = (-zb(id1)) +eta(id1)*yface_gcoef(1,i) + eta(id3)*yface_gcoef(3,i) + eta(id4)*yface_gcoef(4,i)   
            HmnEXP = (-zb(id2)) +eta(id2)
          endif
          
          DETADY = (( (RHOPrim(id1)*HPLUS+pressatm(id1))*HPLUS-  &
                   (RHOPrim(id2)*HmnEXP+pressatm(id2))*HmnEXP)-  & 
                   ((RHOPrim(id1)+RHOPrim(id2))/2)* &
                   ((pressatm(id1)+pressatm(id2))/2. + HPLUS + HmnEXP)* &
                   ( (-zb(Id1)) - (-zb(id2))))/DELTAY 
         
          !inline advection and diffusion 
          if(yface_advF(2,i) .gt. 0) then
            ADVDIF_I = (yface_advdif_i(1,i)+yface_advdif_i(2,i))/2 
          else
            ADVDIF_I = yface_advdif_i(1,i)
          endif
          advf3 = yface_advF(3,i)
          advf4 =  yface_advF(4,i)         
          if(advf4 .gt. 0) then
            ADVDIF_I = ADVDIF_I - (yface_advdif_i(1,advf3) + yface_advdif_i(1,advf4))/2.
          else
            ADVDIF_I = ADVDIF_I - yface_advdif_i(yface_side(ii),advf3)
          endif
         
          ADVDIF_I = ADVDIF_I/DELTAY 

          !cross term advection and diffusion       
          ADVDIF_C = yface_advdif_C(i) - (yface_advdif_C(yface_CadvF(6,i)) + yface_advdif_C(yface_CadvF(4,i)))/2.0
          ADVDIF_C = ADVDIF_C/yface_length(i)        
         
          tauwind = 0.5*(tauwindx(id1)+tauwindx(id2))          
          
          CORIOLIS = 0.0  !-fcoriolis*QXc(i)*0.0
          
          !HGT = (HPLUS+HmnEXP)/2.0_ikind 
          !spd = (uelwc(id1)+uelwc(id2))/2.0_ikind
          !cd =  (cfrict(id1)+Cfrict(id2))/2.0_ikind
          
          !HGT = (HPLUS*DY(id2)+HmnEXP*DY(id1))/(dY(id1)+dy(ID2))  
          !spd = (uelwc(id1)*dy(id2)+uelwc(id2)*dy(id1))/(dY(id1)+dy(ID2))  
          !cd =  (cfrict(id1)*dy(id2)+Cfrict(id2)*dy(id1))/(dY(id1)+dy(ID2))  
          
          !Option A
          HGT = (HPLUS+HmnEXP)/2.0_ikind 
          spd = (uelwc(id1)*yface_gcoef(1,i)+uelwc(id2)*yface_gcoef(2,i)+uelwc(id3)*yface_gcoef(3,i)+uelwc(id4)*yface_gcoef(4,i))/2.0_ikind 
          cd =  (cfrict(id1)*yface_gcoef(1,i)+Cfrict(id2)*yface_gcoef(2,i)+cfrict(id3)*yface_gcoef(3,i)+cfrict(id4)*yface_gcoef(4,i))/2.0_ikind  
          
          hgt=max(0.000001,hgt)
          spd=max(0.000001,spd)
          cd =max(0.000001,cd)
          
          !!Option B
          !HGT = (HPLUS+HmnEXP)/2.0_ikind 
          !spd = abs(yface_vel(i))
          !Nmann = (coefman(id1)*dy(id2)+coefman(id2)*dy(id1))/(dY(id1)+dy(ID2))
          !cd = grav*Nmann*Nmann/HGT**0.333         
          
          yface_QN(i) = (yface_Q(i) + DT*((-GRAV)*DETADY/2.0_ikind &
                      + 0.5*(wavestry(id1) + wavestry(id2)) - ADVDIF_I - ADVDIF_C &
                      + tauwind + CORIOLIS) )/(1+dt*Cd*spd/HGT)
           
        endif  ! end "not a wall"
      enddo
!!$omp end do    
     

!!$OMP END PARALLEL


! !if there are structures then modify flows
! if(structures) then
!   if(SRM_on) then   !rouble mounds
!     do j=1,srmU.ncells
!       I = srmU.cells(j)
!       if(active(i,1)) then
!         ncw = cell2cell(4,i)
!         DELTAX = SRMU.L(j)                      
!         HPLUS = -zb(i)+eta(i)  
!         HmnEXP = -zb(ncw)+eta(ncw)    
!         DETADX = (( (RHOPrim(i)*HPLUS+pressatm(i))*HPLUS
!.              - (RHOPrim(ncw)*HmnEXP+pressatm(ncw))*HmnEXP)-   
!.                ((RHOPrim(i)+RHOPrim(ncw))/2)*
!.                ((pressatm(i)+pressatm(ncw))/2. + HPLUS + HmnEXP)*
!.                          ( -zb(I)+zb(NCw)))/DELTAX      
!         HGT = (HPLUS+HmnEXP)/2.0_ikind         
!         CORIOLIS = fcoriolis*QY(I)*cos(azimuth_fl*deg2rad)  
!         spd = (uv(i)+uv(ncw))/2.0_ikind
!         QXN(I) = (QX(I)+DT*(-GRAV*DETADX/2.0_ikind 
!.         + CORIOLIS) )/
!.                    (1+dt*srmu.a(j)*grav+dt*srmu.b(j)*grav*spd )   
!       endif     
!     enddo
!
!     do j=1,srmV%ncells
!       I = srmV%cells(j)
!       if(active(i,2)) then
!         ncs = cell2cell(3,i)
!         DELTAY = srmV%L(j)
!         HPLUS = -zb(i)+eta(i)  
!         HmnEXP = -zb(ncs)+eta(ncs)
!         DETADY = (( (RHOPrim(i)*HPLUS+pressatm(i))*HPLUS
!.                - (RHOPrim(ncs)*HmnEXP+pressatm(ncs))*HmnEXP)-   
!.                 ((RHOPrim(i)+RHOPrim(ncs))/2)*
!.                ((pressatm(i)+pressatm(ncs))/2. + HPLUS + HmnEXP)*
!.                  ( -zb(I)+zb(NCS)))/DELTAY     
!         HGT = (HPLUS+HmnEXP)/2.0_ikind
!         spd = (uv(i)+uv(ncs))/2.0_ikind
!         CORIOLIS = -fcoriolis*QX(I)*cos(azimuth_fl*deg2rad) 
!         QYN(I) = (QY(I)+DT*(-GRAV*DETADY/2.0_ikind 
!.        + CORIOLIS))/
!.                   (1+dt*srmV%a(j)*grav+dt*srmV%b(j)*grav*spd )   
!       endif
!     enddo 
!   endif
! endif
   
      end subroutine
