      subroutine update_q_and_u_extrapolations_Tel()
    use EXP_Global_def 
    USE EXP_transport_def
      USE EXP_bndcond_def          
    use flow_def
    use comvarbl, only: timehrs
      use bnd_def
      use sed_def
      use exp_telescoping

      implicit none
      integer :: i
              
      do i=1,num_ext_N
        if(yface_q(ext_N(i,2)) .gt. 0) then
          yface_q(ext_N(i,1)) = yface_q(ext_N(i,2))*fac_DW
          yface_vel(ext_N(i,1)) = yface_vel(ext_N(i,2))*fac_DW
          yface_qn(ext_N(i,1)) = yface_q(ext_N(i,1))
        else
          yface_q(ext_N(i,1)) = yface_q(ext_N(i,2))*fac_UW
          yface_vel(ext_N(i,1)) = yface_vel(ext_N(i,2))*fac_UW
          yface_qn(ext_N(i,1)) = yface_q(ext_N(i,1))
        endif
      enddo
      do i=1,num_ext_s
        if(yface_q(ext_s(i,2)) .gt. 0) then
          yface_q(ext_s(i,1)) = yface_q(ext_s(i,2))*fac_DW
          yface_vel(ext_s(i,1)) = yface_vel(ext_s(i,2))*fac_DW
          yface_qn(ext_s(i,1)) = yface_q(ext_s(i,1))
        else
          yface_q(ext_s(i,1)) = yface_q(ext_s(i,2))*fac_UW
          yface_vel(ext_s(i,1)) = yface_vel(ext_s(i,2))*fac_UW
          yface_qn(ext_s(i,1)) = yface_q(ext_s(i,1))
        endif
      enddo

    do i=1,num_ext_W
       if(xface_q(ext_W(i,2)) .gt. 0) then
        xface_q(ext_W(i,1)) = xface_q(ext_W(i,2))*fac_DW
        xface_vel(ext_W(i,1)) = xface_vel(ext_W(i,2))*fac_DW
        xface_qn(ext_W(i,1)) = xface_q(ext_W(i,1))
        else        
        xface_q(ext_W(i,1)) = xface_q(ext_W(i,2))*fac_UW
        xface_vel(ext_W(i,1)) = xface_vel(ext_W(i,2))*fac_UW
        xface_qn(ext_W(i,1)) = xface_q(ext_W(i,1))
      endif
    enddo

    do i=1,num_ext_E
      if(xface_q(ext_E(i,2)) .gt. 0) then
        xface_q(ext_E(i,1)) = xface_q(ext_E(i,2))*fac_DW
        xface_vel(ext_E(i,1)) = xface_vel(ext_E(i,2))*fac_DW
        xface_qn(ext_E(i,1)) = xface_q(ext_E(i,1))
        else
        xface_q(ext_E(i,1)) = xface_q(ext_E(i,2))*fac_UW
        xface_vel(ext_E(i,1)) = xface_vel(ext_E(i,2))*fac_UW
        xface_qn(ext_E(i,1)) = xface_q(ext_E(i,1))
      endif
      enddo

      return
      end subroutine