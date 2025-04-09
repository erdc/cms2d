!*************************************************************
      subroutine update_momentum()
!*************************************************************
      use EXP_Global_def,     only: mixing,ncn,nce,ncs,ncw,ncsw,fuu,ue,gvv,ve,fuv,gvu,advect,qx
      use EXP_Global_def,     only: num_fg_m_cells,fg_m_cells,num_fg_a_cells,fg_a_cells,qy,advectx,advecty,active,rhoprim,dt,qyn,qxn
      use EXP_Structures_def, only: structures,srm_on,srmu,srmv
      use wave_flowgrid_def,  only: wavestrx,wavestry     
      use flow_def,  only: vis,eta,u,v,fcoriolis,grav,uv
      use geo_def,   only: cell2cell,zb,dx,dy,azimuth_fl
      use comvarbl,  only: timesecs
      use const_def, only: pi,deg2rad    
      use met_def,   only: tauwindx,tauwindy,pressatm
      use fric_def,  only: cfrict,uelwc
      use prec_def,  only: ikind
      use size_def,  only: ncells
      
      implicit none
      integer i,j,id
      real(ikind) dyt,dxt,difft,deptht,FluxEXP,hplus,HmnEXP,hgt,coriolis
      real(ikind) spd,deltax,deltay,tauwind
      real(ikind) detadx,detady,cd
      
      !changed order of MIXING and ADVECT 8/20/08 cr
      !calculate mixing terms
      
      if(MIXING) then
!$omp parallel do private(ncn,nce,ncs,ncw,ncsw,dyT,dxT,DiffT,DepthT)
        do i=1,ncells
          ncn = cell2cell(1,i)
          nce = cell2cell(2,i)
          ncs = cell2cell(3,i)
          ncw = cell2cell(4,i)
          ncsw = cell2cell(4,ncs)
          Fuu(i)=- vis(i)*( (-zb(i)) +eta(i))*(uE(nce)-uE(i))/dx(i)
          Gvv(i)=- vis(i)*( (-zb(i)) +eta(i))*(vE(ncn)-vE(i))/dy(i)
          dyT = (dy(i)+dy(ncs))/2.
          dxT = (dx(i)+dx(ncw))/2.
          if(ncsw.gt.0) then
            DiffT = (vis(i)+vis(ncs)+vis(ncw)+vis(ncsw))/4.
            DepthT = ( (-zb(i)) -zb(ncs)-zb(ncw)-zb(ncsw) + eta(i)+eta(ncs)+eta(ncw)+eta(ncsw))/4.
          else
            DiffT = (vis(i)+vis(ncs)+vis(ncw))/3.
            DepthT = ( (-zb(i)) -zb(ncs)-zb(ncw) + eta(i)+eta(ncs)+eta(ncw))/3.
          endif
        
          Fuv(i) = -DiffT*depthT*(uE(i)-uE(ncs))/dyT
          Gvu(i) = -DiffT*depthT*(vE(i)-vE(ncw))/dxT
        enddo
!$omp end parallel do

        !eliminate cross-terms adjacent to boundaries
        do i=1,num_FG_M_cells
          id = FG_M_cells(i)
          Fuv(id) = 0
          Gvu(id) = 0
        enddo 
      ELSE
        Fuu=0
        Gvv=0
        Fuv=0
        Gvu=0
      ENDIF !end diffusion
    
!$omp parallel 
      !calculate advective terms
      if(ADVECT) then
!$OMP do private(nce, FluxEXP, ncn)                     !NLH 07/18/08
        do i=1,ncells
          nce = cell2cell(2,i)
          FluxEXP = (QX(i)+QX(nce))/2.
          if(FluxEXP.gt.0) then
            Fuu(i) = Fuu(i) + FluxEXP*UE(i)
          else
            Fuu(i) = Fuu(i) + FluxEXP*UE(nce)
          endif
          ncn = cell2cell(1,i)
          FluxEXP = (QY(i)+QY(ncn))/2.
          if(FluxEXP.gt.0) then
            Gvv(i) = Gvv(i) + FluxEXP*VE(i)
          else
            Gvv(i) = Gvv(i) + FluxEXP*VE(ncn)
          endif     
          
        enddo
!$omp end do       

!$omp do private(ncs, ncw, ncsw, FluxEXP)      !NLH 07/18/08
        do i=1,ncells
          ncs = cell2cell(3,i)
          ncw = cell2cell(4,i)
          ncsw = cell2cell(4,ncs)
          FluxEXP = (QY(i)+QY(ncw))/2.
          if(FluxEXP.gt.0) then
            Fuv(i) = Fuv(i) + FluxEXP*UE(ncs)
          else
            Fuv(i) = Fuv(i) + FluxEXP*UE(i)
          endif
          FluxEXP = (QX(i)+QX(ncs))/2.
          if(FluxEXP.gt.0) then
            Gvu(i) = Gvu(i) + FluxEXP*VE(ncw)
          else
            Gvu(i) = Gvu(i) + FluxEXP*VE(i)
          endif 
            
        enddo
!$omp end do

!$OMP single
        !added 8/20/08 cr
        !get Fuv and Gvu for any special dummy cells
        do i=1,num_FG_A_cells
          FluxEXP = (QY(FG_A_cells(i,1))+QY(FG_A_cells(i,4)))/2.
          if(FluxEXP.gt.0) then
            Fuv(FG_A_cells(i,1)) = FluxEXP*U(FG_A_cells(i,3))
          else
            Fuv(FG_A_cells(i,1)) = FluxEXP*U(FG_A_cells(i,1))
          endif    
          FluxEXP = (QX(FG_A_cells(i,2))+QX(FG_A_cells(i,3)))/2.
          if(FluxEXP.gt.0) then
            Gvu(FG_A_cells(i,2)) = FluxEXP*V(FG_A_cells(i,4))
          else
            Gvu(FG_A_cells(i,2)) = FluxEXP*V(FG_A_cells(i,1))
          endif    
        enddo    
!$OMP END single      
      endif !ADVECT
!$OMP END PARALLEL

      if(mixing .or. advect) then
!$omp parallel do private(DXT, DYT)              !NLH 07/18/08
        do I=1,ncells
          DXT = (dx(i)+dx(cell2cell(4,i)))/2.0
          ADVECTX(i) = (Fuu(i)-Fuu(cell2cell(4,i)))/DXT + (Fuv(cell2cell(1,i))-Fuv(i))/DY(i)
          DYT = (DY(i)+DY(cell2cell(3,i)))/2.0
          ADVECTY(i) = (Gvv(i)-Gvv(cell2cell(3,i)))/DYT + (Gvu(cell2cell(2,i))-Gvu(i))/DX(i)
        enddo
!$omp end parallel do
      endif
             
      ! update density differential due to salinity  
      !CR - 11/20/2009 modified to be able to turn this off.
!      if(saltrans .and. saltsimD) then
!!$omp parallel do       
!       do i=1,ncells
!         RHOPrim(i) = 0.000808*salt(i).Conc
!       enddo
!!$omp end parallel do       
!      endif   
    
      if (any(eta.gt.1.0e+3)) call diag_print_error('ETA values greater than 100 are evident.  Reduce the timestep.')

!$omp parallel do private(ncw,deltax,hplus,hmnexp,detadx,tauwind,coriolis,hgt,spd,cd)
      do i=1,ncells
        if(active(i,1)) then
          ncw = cell2cell(4,i)
          DELTAX = (DX(I)+DX(NCW))/2.0_ikind
          HPLUS = -zb(i)+eta(i)  
          HmnEXP = -zb(ncw)+eta(ncw)                         
          DETADX = (( (RHOPrim(i)*HPLUS +pressatm(i))*HPLUS -           &
                   (RHOPrim(ncw)*HmnEXP +pressatm(ncw))*HmnEXP) -       &
                   ((RHOPrim(i)+RHOPrim(ncw))/2) *                      &
                   ((pressatm(i)+pressatm(ncw))/2. + HPLUS + HmnEXP) *  &
                   ( (-zb(I))- (-zb(NCw))))/DELTAX   
          tauwind = 0.5*(tauwindx(i)+tauwindx(ncw))          
          CORIOLIS = fcoriolis*QY(I)*cos(azimuth_fl*deg2rad)          
          HGT = (HPLUS+HmnEXP)/2.0_ikind 
          spd = (uelwc(i)+uelwc(ncw))/2.0_ikind
          cd =  (cfrict(i)+Cfrict(ncw))/2.0_ikind
          QXN(I) = (QX(I)+DT*((-GRAV)*DETADX/2.0_ikind + wavestrx(i) &
                 - ADVECTX(I) + tauwind + CORIOLIS) )/(1+dt*Cd*spd/HGT)
        endif
        
        if(active(i,2)) then
          ncs = cell2cell(3,i)
          DELTAY = (DY(I)+DY(NCS))/2.0_ikind
          HPLUS  = (-zb(i)) +eta(i)  
          HmnEXP = (-zb(ncs)) +eta(ncs)
          DETADY = (( (RHOPrim(i)*HPLUS+pressatm(i))*HPLUS -           &
                   (RHOPrim(ncs)*HmnEXP+pressatm(ncs))*HmnEXP) -       &
                   ((RHOPrim(i)+RHOPrim(ncs))/2) *                     &
                   ((pressatm(i)+pressatm(ncs))/2. + HPLUS + HmnEXP)*  &
                   ( (-zb(I))- (-zb(NCS))))/DELTAY    
          CORIOLIS = (-fcoriolis) *QX(I)*cos(azimuth_fl*deg2rad)  
          tauwind =0.5*(tauwindy(i)+tauwindy(ncs))
          HGT = (HPLUS+HmnEXP)/2.0_ikind                           
          spd = (uelwc(i)+uelwc(ncs))/2.0_ikind
          cd =  (cfrict(i)+Cfrict(ncs))/2.0_ikind
          QYN(I) = (QY(I)+DT*(-GRAV*DETADY/2.0_ikind + wavestry(i) &
                 - ADVECTY(i) + tauwind + CORIOLIS))/(1+DT*Cd*spd/HGT)
        endif
      enddo
!$omp end parallel do   

      !if there are structures then modify flows
      if(structures) then
        if(SRM_on) then   !rouble mounds
          do j=1,SRMU%ncells
            I = SRMU%cells(j)
            if(active(i,1)) then
              ncw = cell2cell(4,i)
              DELTAX = SRMU%L(j)                      
              HPLUS = -zb(i)+eta(i)  
              HmnEXP = -zb(ncw)+eta(ncw)    
              DETADX = (( (RHOPrim(i)*HPLUS+pressatm(i))*HPLUS -        &
                    (RHOPrim(ncw)*HmnEXP+pressatm(ncw))*HmnEXP) -       &
                    ((RHOPrim(i)+RHOPrim(ncw))/2) *                     &
                    ((pressatm(i)+pressatm(ncw))/2. + HPLUS + HmnEXP) * &
                    (-zb(I)+zb(NCw)))/DELTAX      
              HGT = (HPLUS+HmnEXP)/2.0_ikind         
              CORIOLIS = fcoriolis*QY(I)*cos(azimuth_fl*deg2rad)  
              spd = (uv(i)+uv(ncw))/2.0_ikind
              QXN(I) = (QX(I)+DT*((-GRAV)*DETADX/2.0_ikind + CORIOLIS) )/ &
                       (1+dt*SRMU%a(j)*grav+dt*SRMU%b(j)*grav*spd)
            endif     
          enddo

          do j=1,srmV%ncells
            I = srmV%cells(j)
            if(active(i,2)) then
              ncs = cell2cell(3,i)
              DELTAY = srmV%L(j)
              HPLUS = -zb(i)+eta(i)  
              HmnEXP = -zb(ncs)+eta(ncs)
              DETADY = (( (RHOPrim(i)*HPLUS+pressatm(i))*HPLUS -        &
                    (RHOPrim(ncs)*HmnEXP+pressatm(ncs))*HmnEXP) -       &
                    ((RHOPrim(i)+RHOPrim(ncs))/2)*                      &
                    ((pressatm(i)+pressatm(ncs))/2. + HPLUS + HmnEXP) * &
                    (-zb(I)+zb(NCS)))/DELTAY     
              HGT = (HPLUS+HmnEXP)/2.0_ikind
              spd = (uv(i)+uv(ncs))/2.0_ikind
              CORIOLIS = (-fcoriolis) *QX(I)*cos(azimuth_fl*deg2rad) 
              QYN(I) = (QY(I)+DT*((-GRAV)*DETADY/2.0_ikind + CORIOLIS))/ &
                       (1+dt*srmV%a(j)*grav+dt*srmV%b(j)*grav*spd )   
            endif
          enddo 
        endif
      endif
   
      end subroutine
