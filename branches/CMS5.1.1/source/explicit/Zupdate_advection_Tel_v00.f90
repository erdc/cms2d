      subroutine update_advection_tel()
      use EXP_Global_def, only: imix, iadv
	  use flow_def, only: eta, vis
	  use geo_def,  only: dx, dy, zb 
      use prec_def, only: ikind
      use EXP_TELESCOPING, only: numtbxfaces,tbxfaces,xface_advf,xface_cells,xface_q,xface_vel,xface_advdif_i,xface_wall,xface_cadvf,xface_length,xface_advdif_C
      USE EXP_TELESCOPING, only: numtbyfaces,tbyfaces,yface_advf,yface_cells,yface_q,yface_vel,yface_advdif_i,yface_wall,yface_cadvf,yface_length,yface_advdif_C
      !use const_def, only: pi,deg2rad    
      !use met_def, only: tauwindx,tauwindy,pressatm
      !use wave_flowgrid_def, only: wavestrx,wavestry     
      !use fric_def, only: cfrict,uelwc
      !USE EXP_transport_def 
      !use NupdateMod
      !use sal_def  
      
      implicit none
      integer i,j,id,id1,id2,id3,ii,istop   !Chris Reed - 10/20/2016
      integer advf1,advf2,cadvf1,cadvf2,cadvf5,c1,c2,c15,c25,cadvf3,cadvf7,cadvf8
      real(ikind) dyt,dxt,FluxEXP
      real(ikind) vadvectp1,vadvectp2,vadvectm3
      real(ikind) uadvectp1,uadvectp2,uadvectm3
      real(ikind) diff_P1,diff_p2,diff_P,diff_T,depth_T

!!$omp parallel             !Remove comments when debugging has completed    
!!$omp do private(i,id1,id2,advf1,advf2,FluxEXP,UADVECTP1,UADVECTP2,DIFF_T,DIFF_P,DYT,UADVECTM3,DIFF_P1,DIFF_P2,depth_t,cadvf1,cadvf2,cadvf7,cadvf8,cadvf5,cadvf3,c1,c2,c15,c25)
      do ii=1,numTBXfaces
        i=TBXfaces(ii) 
          
        if(xface_advF(1,i) .gt. 0) then  !Filter out faces on east side of grid
          id1 = xface_cells(1,i)
          id2 = xface_cells(2,i) 
          
          !inline advection   
          advf1 = xface_advF(1,i)
          advf2 = xface_advF(2,i)
          FluxEXP = (xface_q(i) + xface_Q(advf1))/2.0
          if(FluxEXP.gt.0) then
            UADVECTP1 = FluxEXP*xface_vel(i)
          else
            UADVECTP1 = FluxEXP*xface_vel(advf1)
          endif
          if(advf2 .gt. 0) then
            FluxEXP = (xface_q(i) + xface_q(advf2))/2.0
            if(FluxEXP.gt.0) then
              UADVECTP2 = FluxEXP*xface_vel(i)
            else
              UADVECTP2 = FluxEXP*xface_vel(advf2)
            endif
          else
            UADVECTP2 = 0
          endif
        
          ! if(i .eq. 1 .or. i .eq. 2) write(*,*)i,advf1,FluxEXP,xface_vel(i),xface_vel(ADVF1)        

          !inline diffusion
          DIFF_P2 = 0.0
          DIFF_P1 = 0.0
          if(imix .eq. 1) then
            ID1 = xface_cells(1,i)
            Diff_P1 = vis(ID1)*(eta(ID1)-zb(ID1))*(xface_vel(i)-xface_vel(advf1))/dx(ID1)         
            if(advf2 .gt. 0) then            
              Diff_P2 = vis(ID1)*(eta(ID1)-zb(ID1))*(xface_vel(i)-xface_vel(advf2))/dx(ID1)                 
            else
              Diff_P2 = 0.0
            endif
          endif
        
          xface_advdif_I(1,i) = IADV*UADVECTP1 + DIff_P1
          xface_advdif_I(2,i) = IADV*UADVECTP2 + Diff_P2
        endif  !east boundary faces filter
        
        if(.not. xface_wall(i)) then
          !cross term advection 
          cadvf1 = xface_CadvF(1,i)
          cadvf2 = xface_CadvF(2,i)
          cadvf5 = xface_CadvF(5,i)
          cadvf3 = xface_CadvF(3,i)  
          cadvf7 = xface_CadvF(7,i)
          cadvf8 = xface_CadvF(8,i) 
          FLuxExp = (yface_q(CadvF1)*yface_length(CadvF1) + yface_q(cadvf2)*yface_length(cadvf2)+yface_q(CadvF7)*yface_length(CadvF7) + yface_q(cadvf8)*yface_length(cadvf8))/ &       
                    (yface_length(CadvF1) + yface_length(cadvf2)+yface_length(CadvF7) + yface_length(cadvf8)) 
          if(FluxEXP.lt.0) then
            UADVECTM3 = FluxEXP*(xface_vel(cadvf5)+xface_vel(cadvf3))/2.
          else
            UADVECTM3 = FluxEXP*xface_vel(i)
          endif 

          !cross-term diffusion
          DIFF_P = 0
          if(imix .eq. 1) then
            c1 = xface_cells(1,i)
            c2 = xface_cells(2,i)
            c15 = xface_cells(1,cadvf5)
            c25 = xface_cells(2,cadvf5)
            DIFF_T = ( vis(c1) + vis(c2) + vis(c15) + vis(c25) )/4.0
            DEPTH_T = ( eta(c1)+eta(c2) +  eta(c15)+eta(c25)  )/4.0 -(zb(c1)+zb(c2) + zb(c15)+zb(c25) )/4.0  
            DYT = (xface_length(i)+xface_length(cadvf5))/2.0
            DIFF_P = DIFF_T*DEPTH_T*(xface_vel(cadvf5) -xface_vel(i)) / DYT
          endif
        
!Chris Reed - 10/20/2016 - testing for issues
        istop =0
        if(c15 .eq. 0) then
            write(*,*)'c15 ',i
            istop=1
            endif
        if(c25 .eq. 0) then
            write(*,*)'c25 ',i
            istop = 1
            endif
         if(cadvf1  .eq. 0) then
            write(*,*)'cadvf1  ',i
            istop=1
            endif
        if(cadvf2  .eq. 0) then
            write(*,*)'cadvf2  ',i
            istop = 1
            endif       
          if(cadvf5  .eq. 0) then
            write(*,*)'cadvf5  ',i
            istop=1
            endif
        if(cadvf3  .eq. 0) then
            write(*,*)'cadvf3  ',i
            istop = 1
            endif        
          if(cadvf7  .eq. 0) then
            write(*,*)'cadvf7  ',i
            istop=1
            endif
        if(cadvf8  .eq. 0) then
            write(*,*)'cadvf8  ',i
            istop = 1
        endif 
        
        if(istop .eq. 1) stop
        
        
        
          xface_advdif_C(i) = IADV*UADVECTM3 - DIFF_P
         
          !Qyc(i) = (Qy(xface_CadvF(1,i))+Qy(xface_CadvF(2,i))+Qy(xface_CadvF(3,i))+Qy(xface_CadvF(4,i)))/4.0         
        endif
      enddo
!!$omp end do  

!!$omp do private(i,id1,id2,vadvectp1,vadvectp2,fluxEXP,diff_t,vadvectm3,diff_p1,diff_p2,diff_p,dxt,depth_t,advf1,advf2,cadvf1,cadvf2,cadvf7,cadvf8,cadvf5,cadvf3,c1,c2,c15,c25)
      do ii=1,numTBYfaces
        i=TBYfaces(ii)
        if(yface_advF(1,i) .gt. 0) then  !Filter out faces on North boundary
          id1 = yface_cells(1,i)
          id2 = yface_cells(2,i)
         
          !inline advection
          advf1 = yface_advF(1,i)
          advf2 = yface_advF(2,i)
          FluxEXP = (yface_q(i) + yface_q(advf1))/2.0
          if(FluxEXP.gt.0) then
            VADVECTP1 = FluxEXP*yface_vel(i)
          else
            VADVECTP1 = FluxEXP*yface_vel(advf1)
          endif
          if( yface_advF(2,i) .gt. 0) then
            FluxEXP = (yface_q(i) + yface_q(advf2))/2.0
            if(FluxEXP.gt.0) then
              VADVECTP2 = FluxEXP*yface_vel(i)
            else
              VADVECTP2 = FluxEXP*yface_vel(advf2)
            endif
          else
            VADVECTP2 = 0.0   
          endif               
          
          !inline diffusion
          DIFF_P1 = 0.0
          DIFF_P2 = 0.0
          if(imix .eq. 1) then
            ID1 = yface_cells(1,i)
            Diff_P1 = vis(ID1)*(eta(ID1)-zb(ID1))*(yface_vel(i)-yface_vel(advf1))/dy(ID1)        
            if(yface_advF(2,i) .gt. 0) then            
              Diff_P2 = vis(ID1)*(eta(ID1)-zb(ID1))*(yface_vel(i)-yface_vel(advf2))/dy(ID1)               
            else
              Diff_P2 = 0.0 
            endif
          endif
         
          yface_advdif_I(1,i) = IADV*VADVECTP1 + DIff_P1
          yface_advdif_I(2,i) = IADV*VADVECTP2 + DIff_P2
         
          !if(i.eq.cellfaces(5,6805)) write(*,*)i,yface_advdif_I(1,i),yface_advdif_I(2,i)
          !if(i.eq.cellfaces(5,6806)) write(*,*)i,yface_advdif_I(1,i),yface_advdif_I(2,i)
          !if(i.eq.cellfaces(5,4020)) write(*,*)i,yface_advdif_I(1,i),yface_advdif_I(2,i)          
          !if(i.eq.cellfaces(5,3973)) write(*,*)i,yface_advdif_I(1,i),yface_advdif_I(2,i)           
        endif  !north boundary faces filter

        if(.not. yface_wall(i)) then      
          !advection cross-term
          cadvf1 = yface_CadvF(1,i)
          cadvf2 = yface_CadvF(2,i)
          cadvf5 = yface_CadvF(5,i)
          cadvf3 = yface_CadvF(3,i)   
          cadvf7 = yface_CadvF(7,i)
          cadvf8 = yface_CadvF(8,i)          
          FLuxExp = (xface_q(cadvf1)*xface_length(cadvf1) + xface_q(cadvf2)*xface_length(cadvf2)+xface_q(cadvf7)*xface_length(cadvf7) + xface_q(cadvf8)*xface_length(cadvf8))/ &       
                    (xface_length(cadvf1) + xface_length(cadvf2)+xface_length(cadvf7) + xface_length(cadvf8)) 
          if(FluxEXP.lt.0) then
            VADVECTM3 = FluxEXP*(yface_vel(cadvf5)+yface_vel(cadvf3))/2.
          else
            VADVECTM3= FluxEXP*yface_vel(i)
          endif         
          !cross-term diffusion
          DIFF_P = 0.0
          if(imix .eq. 1) then
            c1 = yface_cells(1,i)
            c2 = yface_cells(2,i)
            c15 = yface_cells(1,cadvf5)
            c25 = yface_cells(2,cadvf5)
            DIFF_T = ( vis(c1) + vis(c2) +  vis(c15) + vis(c25) )/4.0
            DEPTH_T = ( eta(c1)+eta(c2) +  eta(c15)+eta(c25)  )/4.0  -(zb(c1)+zb(c2) +  zb(c15)+zb(c25) )/4.0  
            DXT = (yface_length(i)+yface_length(cadvf5))/2.0
            DIFF_P = DIFF_T*DEPTH_T*(yface_vel(cadvf5) -yface_vel(i)) / DXT
          endif
        
          yface_advdif_C(i) = IADV*VADVECTM3 - DIFF_P
          !Qxc(i) = (Qx(yface_CadvF(1,i))+Qx(yface_CadvF(2,i))+Qx(yface_CadvF(3,i))+Qx(yface_CadvF(4,i)))/4.0           
        endif
      enddo
!!$omp end do    
!!$OMP END PARALLEL

    return
    end subroutine
    
    
