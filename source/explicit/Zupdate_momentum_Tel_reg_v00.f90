    subroutine update_momentum_tel_reg()
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
      integer i,j,id,id1,id2,id3,ii
      real(ikind) dyt,dxt,difft,deptht,hplus,HmnEXP,hgt,coriolis
      real(ikind) spd,deltax,deltay,tauwind
      real(ikind) detadx,detady,cd,depth,FluxEXP
      real(ikind) advdif_i,advdif_c
      real(ikind) Nmann

      
!!$omp parallel 
        
      ! update density differential due to salinity  
      !  if(saltrans .and. saltsimD) then
!!$omp do       
      !    do i=1,ncells
             !RHOPrim(i) = 1.0 + 0.000808*salt(i).Conc
      !    enddo
!!$omp end do       
      !  endif   
!!$omp end parallel
      
!There was an OpenMP issue here with the following 2 DO loops.  I made each one its own PARALLEL DO instead of putting individual DOs inside one PARALLEL block.  MEB  01/16/19
     
!$omp parallel do private(i,id1,id2,deltax,hplus,HmnEXP,detadx,advdif_I,advdif_C,tauwind,Coriolis,hgt,spd,cd,Qyc)
      do ii=1,numREGXfaces
        i=REGXfaces(ii) 
                  
        if(.not. xface_wall(i)) then
          id1 = xface_cells(1,i)
          id2 = xface_cells(2,i)
      
!         DELTAX = (dx(id1)+dx(id2))/2.     
!         depth = (-zb(id1)+eta(id1) -zb(id2)+eta(id2))/2.           
!         HPLUS = eta(id1)  + pressatm(id1) + RHOPrim(id1)
!         HmnEXP = eta(id2)  + pressatm(id2) + RHOPrim(id2)
!         DETADX = depth*(HPLUS-HmnEXP)/DELTAX 
          
          DELTAX = (DX(Id1)+DX(id2))/2.0_ikind
          HPLUS  = (-zb(id1)) +eta(id1)  
          HmnEXP = (-zb(id2)) +eta(id2)                         
          DETADX = (( (RHOPrim(id1)*HPLUS +pressatm(id1))*HPLUS - &
                   (RHOPrim(id2)*HmnEXP +pressatm(id2))*HmnEXP)- &
                   ((RHOPrim(id1)+RHOPrim(id2))/2)* &
                   ((pressatm(id1)+pressatm(id2))/2. + HPLUS + HmnEXP)* &
                   ( (-zb(Id1)) - (-zb(id2))))/DELTAX   
           
          !inline advection and diffusion 
          ADVDIF_I = xface_advdif_i(1,i) - xface_advdif_i(1,xface_advF(3,i))           
          ADVDIF_I = ADVDIF_I/DELTAX 
          
          !cross term advectionand diffusion       
          ADVDIF_C = xface_advdif_C(i) - xface_advdif_C(xface_CadvF(6,i))
          ADVDIF_C = ADVDIF_C/xface_length(i)

          tauwind = 0.5*(tauwindx(id1)+tauwindx(id2))          
          
          CORIOLIS = 0.0 !fcoriolis*QYc(i)*0.0
          
          !option A
          HGT = (HPLUS+HmnEXP)/2.0_ikind 
          spd = (uelwc(id1)+uelwc(id2))/2.0_ikind
          cd =  (cfrict(id1)+Cfrict(id2))/2.0_ikind
          
          !option B
          !HGT = (HPLUS+HmnEXP)/2.0_ikind 
          !spd = abs(xface_vel(i))         
          !Nmann = (coefman(id1)*dy(id2)+coefman(id2)*dy(id1))/(dY(id1)+dy(ID2))
          !cd = grav*Nmann*Nmann/HGT**0.333    
          
          hgt=max(0.000001,hgt)
          spd=max(0.000001,spd)
          cd =max(0.000001,cd)
          
          xface_QN(i) = (xface_Q(i) + DT*( (-GRAV) *DETADX/2.0_ikind &
                        + 0.5*(wavestrx(id1) + wavestrx(id2)) - ADVDIF_I - ADVDIF_C &
                        + tauwind + CORIOLIS) )/(1+dt*Cd*spd/HGT)
          
          !xface_QN(i) = xface_Q(i) + DT*(- ADVDIF_I - ADVDIF_C )             
        endif  ! end "not a wall"
      enddo
!$omp end parallel do    
     
!$omp parallel do private(i,id1,id2,deltay,hplus,HmnEXP,detady,advdif_i,advdif_c,tauwind,Coriolis,hgt,spd,cd,QXc)
      do ii=1,numREGYfaces
        i=REGYfaces(ii) 
          
        if(.not. yface_wall(i)) then
          id1 = yface_cells(1,i)
          id2 = yface_cells(2,i)
      
!         DELTAY =(dy(id1)+dy(id2))/2.
!         depth = (-zb(id1)+eta(id1) -zb(id2)+eta(id2))/2.  
!         HPLUS = eta(id1) + pressatm(id1) + RHOPrim(id1)
!         HmnEXP = eta(id2)  + pressatm(id2) + RHOPrim(id2)
!         DETADY = depth*(HPLUS-HmnEXP)/DELTAY

          DELTAY = (DY(Id1)+DY(id2))/2.0_ikind
          HPLUS  = (-zb(id1))+eta(id1)  
          HmnEXP = (-zb(id2))+eta(id2)

          DETADY = (( (RHOPrim(id1)*HPLUS +pressatm(id1))*HPLUS - &
                   (RHOPrim(id2)*HmnEXP +pressatm(id2))*HmnEXP)- & 
                   ((RHOPrim(id1)+RHOPrim(id2))/2)* &
                   ((pressatm(id1)+pressatm(id2))/2. + HPLUS + HmnEXP)* &
                   ( (-zb(Id1))- (-zb(id2))))/DELTAY  
         
          !inline advection and diffusion 
          ADVDIF_I = yface_advdif_i(1,i) - yface_advdif_i(1,yface_advF(3,i))       
          ADVDIF_I = ADVDIF_I/DELTAY 
  
          !cross term advection and diffusion       
          ADVDIF_C = yface_advdif_C(i)- yface_advdif_C(yface_CadvF(6,i))
          ADVDIF_C = ADVDIF_C/yface_length(i)        
         
          tauwind = 0.5*(tauwindx(id1)+tauwindx(id2))          
          
          CORIOLIS = 0.0  !-fcoriolis*QXc(i)*0.0
          
          !!Option A
          HGT = (HPLUS+HmnEXP)/2.0_ikind 
          spd = (uelwc(id1)+uelwc(id2))/2.0_ikind
          cd =  (cfrict(id1)+Cfrict(id2))/2.0_ikind
          
          !!option B
          !HGT = (HPLUS+HmnEXP)/2.0_ikind 
          !spd = abs(yface_vel(i))          
          !Nmann = (coefman(id1)*dy(id2)+coefman(id2)*dy(id1))/(dY(id1)+dy(ID2))
          !cd = grav*Nmann*Nmann/HGT**0.333          

          yface_QN(i) = (yface_Q(i) + DT*((-GRAV)*DETADY/2.0_ikind &
                      + 0.5*(wavestry(id1) + wavestry(id2)) - ADVDIF_I - ADVDIF_C &
                      + tauwind + CORIOLIS) )/(1+dt*Cd*spd/HGT)

          !if(i.eq.cellfaces(5,6805)) write(*,*)"- - ",i,ADVDIF_I,ADVDIF_C
          !if(i.eq.cellfaces(5,6806)) write(*,*)"- - ",i,ADVDIF_I,ADVDIF_C
          !if(i.eq.cellfaces(5,4020)) write(*,*)"- - ",i,ADVDIF_I,ADVDIF_C         
          !if(i.eq.cellfaces(5,3973)) write(*,*)"- - ",i,ADVDIF_I,ADVDIF_C           
          !
          ! yface_QN(i)= yface_Q(i) + DT*(- ADVDIF_I - ADVDIF_C )                
           
        endif  ! end "not a wall"
      enddo
!$omp end parallel do    

! !if there are structures then modify flows
! if(structures) then
!   if(SRM_on) then   !rouble mounds
!     do j=1,SRMU%ncells
!       I = SRMU%cells(j)
!       if(active(i,1)) then
!         ncw = cell2cell(4,i)
!         DELTAX = SRMU%L(j)                      
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
!.                    (1+dt*SRMU%a(j)*grav+dt*SRMU%b(j)*grav*spd )   
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

      return
      end subroutine
