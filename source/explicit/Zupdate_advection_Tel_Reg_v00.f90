!*******************************************************************************
     subroutine update_advection_tel_reg()
!*******************************************************************************
     use EXP_Global_def,    only: imix, iadv
     use EXP_TELESCOPING,   only: numregxfaces, regxfaces, xface_cells, xface_advf, xface_q, xface_vel, xface_advdif_i, xface_wall, xface_cadvf, xface_advdif_c
     use EXP_TELESCOPING,   only: numregyfaces, regyfaces, yface_cells, yface_advf, yface_q, yface_vel, yface_advdif_i, yface_wall, yface_cadvf, yface_advdif_c
     use EXP_TELESCOPING,   only: xface_length, yface_length, numspecx, numspecy, specx, specy
     use wave_flowgrid_def, only: wavestrx,wavestry     
     use flow_def,  only: vis, eta
     use geo_def,   only: zb, dx, dy
     use const_def, only: pi, deg2rad    
     use met_def,   only: tauwindx, tauwindy, pressatm
     use fric_def,  only: cfrict, uelwc
     use prec_def,  only: ikind
      
     implicit none
     integer i,j,id,id1,id2,id3,ii
     real(ikind) dyt,dxt,FluxEXP
     real(ikind) vadvectp1,vadvectp2,vadvectm3
     real(ikind) uadvectp1,uadvectp2,uadvectm3
     real(ikind) diff_P1,diff_p2,diff_P,diff_T,depth_T
      
     integer ADVF1,CADV5,C1,C2,C15,C25

!$omp parallel       
!$omp do private(i,id1,id2,advf1,uadvectp1,fluxEXP,diff_t,uadvectm3,diff_p1,dyt,depth_t,diff_p,c15,c25,cadv5)

     do ii=1,numREGXfaces
       i=REGXfaces(ii)  
       id1 = xface_cells(1,i)
       id2 = xface_cells(2,i) 
       if(xface_advF(1,i) .gt. 0) then  !Filter out faces on east side of grid
         ADVF1 = xface_advF(1,i)
          
         !inline advection   
         FluxEXP = (xface_q(i) + xface_q(ADVF1))/2.0
         if(FluxEXP.gt.0) then
           UADVECTP1 = FluxEXP*xface_vel(i)
         else
           UADVECTP1 = FluxEXP*xface_vel(ADVF1)
         endif
          
         !inline diffusion
         DIFF_p1 = 0.0
         if(imix .eq. 1) then
           ID1 = xface_cells(1,i)
           Diff_P1 = vis(ID1)*(eta(ID1)-zb(ID1))*(xface_vel(i)-xface_vel(ADVF1))/dx(ID1)         
         endif
        
         xface_advdif_I(1,i) = IADV*UADVECTP1 + DIff_P1
         xface_advdif_I(2,i) = 0
       endif  !east boundary faces filter
        
       if(.not. xface_wall(i)) then
         CADV5 =  xface_CadvF(5,i)
      
         !cross term advection 
         FLuxExp = (yface_q(xface_CadvF(1,i)) + yface_q(xface_CadvF(2,i)))/2     
         if(FluxEXP.lt.0) then
           UADVECTM3 = FluxEXP*xface_vel(CADV5)
         else
           UADVECTM3 = FluxEXP*xface_vel(i)
         endif 
 
         !C1 = xface_cells(1,i)
         !C2 = xface_cells(2,i)
         C15 = xface_cells(1,xface_CadvF(5,i))
         C25 = xface_cells(2,xface_CadvF(5,i))
 
         !cross-term diffusion
         diff_p = 0.0
         if(imix .eq. 1) then
           DIFF_T = ( vis(id1) + vis(id2) + vis(C15) + vis(C25) )/4.0
           DEPTH_T = ( eta(id1)+eta(id2) + eta(C15)+eta(C25) )/4.0 &
                    -(zb(id1)+zb(id2) + zb(C15)+zb(C25) )/4.0  
           DYT = (xface_length(i)+xface_length(CADV5))/2.0
           DIFF_P = DIFF_T*DEPTH_T*(xface_vel(CADV5) -xface_vel(i)) / DYT
         endif
        
         xface_advdif_C(i) = IADV*UADVECTM3 - DIFF_P
       endif
     enddo
!$omp end do  

     !i=6805
     !write(*,"(9i6)")i,cellfaces(:,6805)
     !i=6806          
     !write(*,"(9i6)")i,cellfaces(:,6806)
     !i=4020         
     !write(*,"(9i6)")i,cellfaces(:,4020)
     !i=3973         
     !write(*,"(9i6)")i,cellfaces(:,3973)          

!$omp do private(i,id1,id2,vadvectp1,fluxEXP,diff_t,vadvectm3,diff_p1,dxt,depth_t,diff_p,c15,c25,cadv5)
     do ii=1,numREGYfaces
       i=REGYfaces(ii) 
       id1 = yface_cells(1,i)
       id2 = yface_cells(2,i)
      
       if(yface_advF(1,i) .gt. 0) then  !Filter out faces on North boundary
         ADVF1 = yface_advF(1,i)  
      
         !inline diffusion     
         FluxEXP = (yface_q(i) + yface_Q(ADVF1))/2.0
         if(FluxEXP.gt.0) then
           VADVECTP1 = FluxEXP*yface_vel(i)
         else
           VADVECTP1 = FluxEXP*yface_vel(ADVF1)
         endif

        !inline diffusion
         diff_p1 = 0.0
         if(imix .eq. 1) then
           ID1 = yface_cells(1,i)
           Diff_P1 = vis(ID1)*(eta(ID1)-zb(ID1))*(yface_vel(i)-yface_vel(ADVF1))/dy(ID1)        
         endif
         
         yface_advdif_I(1,i) = IADV*VADVECTP1 + DIff_P1
         yface_advdif_I(2,i) = 0
         
       endif  !north boundary faces filter

       if(.not. yface_wall(i)) then   
         CADV5 =  yface_CadvF(5,i)          
      
         !advection cross-term
         FLuxExp = (xface_q(yface_CadvF(1,i)) + xface_q(yface_CadvF(2,i)))/2      
         if(FluxEXP.lt.0) then
           VADVECTM3 = FluxEXP*yface_vel(CADV5)
         else
           VADVECTM3= FluxEXP*yface_vel(i)
         endif     

         !C1 = Yface_cells(1,i)
         !C2 = Yface_cells(2,i)
         C15 = Yface_cells(1,yface_CadvF(5,i))
         C25 = Yface_cells(2,yface_CadvF(5,i))          
           
         !cross-term diffusion
         diff_P = 0.0
         if(imix .eq. 1) then
           DIFF_T = ( vis(id1) + vis(id2) +  vis(C15) + vis(c25) )/4.0
           DEPTH_T = ( eta(id1)+eta(id2) + eta(C15)+eta(C25)  )/4.0 &
                    -(zb(id1)+zb(id2) +  zb(C15)+zb(C25) )/4.0  
           DXT = (yface_length(i)+yface_length(CADV5))/2.0
           DIFF_P = DIFF_T*DEPTH_T*(yface_vel(CADV5) -yface_vel(i)) / DXT
         endif
        
         yface_advdif_C(i) = IADV*VADVECTM3 - DIFF_P
     
       endif
     enddo
!$omp end do  

!$OMP END PARALLEL

     do ii=1,numspecX
       i=specX(ii)
       !inline advection   
       advf1 = xface_advF(1,i)
       FluxEXP = (xface_q(i) + xface_Q(advf1))/2.0
       if(FluxEXP.gt.0) then
         UADVECTP1 = FluxEXP*xface_vel(i)
       else
         UADVECTP1 = FluxEXP*xface_vel(advf1)
       endif        
     enddo
         
     do ii=1,numspecY
       i=specY(ii)        
       !inline advection
       advf1 = yface_advF(1,i)
       FluxEXP = (yface_q(i) + yface_q(advf1))/2.0
       if(FluxEXP.gt.0) then
         VADVECTP1 = FluxEXP*yface_vel(i)
       else
         VADVECTP1 = FluxEXP*yface_vel(advf1)
       endif
     enddo

     return
     end subroutine