      subroutine update_q_and_u_extrapolations()
      use EXP_Global_def 
      USE EXP_transport_def
      USE EXP_bndcond_def	      
      use flow_def
      use comvarbl, only: timehrs
      use bnd_def
      use sed_def

      implicit none
      integer :: i
            
      do i=1,num_ext_N
        if(qy(ext_N(i,2)).gt.0) then
          qy(ext_N(i,1)) = qy(ext_N(i,2))*fac_DW
          vE(ext_N(i,1)) = vE(ext_N(i,2))*fac_DW
          qyn(ext_N(i,1)) = qy(ext_N(i,1))
        else
          qy(ext_N(i,1)) = qy(ext_N(i,2))*fac_UW
          vE(ext_N(i,1)) = vE(ext_N(i,2))*fac_UW
          qyn(ext_N(i,1)) = qy(ext_N(i,1))
        endif
      enddo

      do i=1,num_ext_S
        if(qy(ext_S(i,2)).gt.0) then
          qy(ext_S(i,1)) = qy(ext_S(i,2))*fac_UW
          vE(ext_S(i,1)) = vE(ext_S(i,2))*fac_UW
          qyn(ext_S(i,1)) = qy(ext_S(i,1))
        else
          qy(ext_S(i,1)) = qy(ext_S(i,2))*fac_DW
          vE(ext_S(i,1)) = vE(ext_S(i,2))*fac_DW
          qyn(ext_S(i,1)) = qy(ext_S(i,1))
        endif
      enddo

      do i=1,num_ext_E
        if(qx(ext_E(i,2)).gt.0) then
          qx(ext_E(i,1)) = qx(ext_E(i,2))*fac_DW
          uE(ext_E(i,1)) = uE(ext_E(i,2))*fac_DW
          qxn(ext_E(i,1)) = qx(ext_E(i,1))
        else
          qx(ext_E(i,1)) = qx(ext_E(i,2))*fac_UW
          uE(ext_E(i,1)) = uE(ext_E(i,2))*fac_UW
          qxn(ext_E(i,1)) = qx(ext_E(i,1))
        endif
      enddo

      do i=1,num_ext_W
        if(qx(ext_W(i,2)).gt.0) then
          qx(ext_W(i,1)) = qx(ext_W(i,2))*fac_UW
          uE(ext_W(i,1)) = uE(ext_W(i,2))*fac_UW
          qxn(ext_W(i,1)) = qx(ext_W(i,1))
        else
          qx(ext_W(i,1)) = qx(ext_W(i,2))*fac_DW
          uE(ext_W(i,1)) = uE(ext_W(i,2))*fac_DW
          qxn(ext_W(i,1)) = qx(ext_W(i,1))
        endif
      enddo

      return
      end subroutine