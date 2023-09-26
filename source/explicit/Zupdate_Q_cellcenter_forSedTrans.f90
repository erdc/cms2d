!*************************************************************
      subroutine update_Q_forSedTrans()
!*************************************************************
      !this is done to prevent upstream exptrapolation of velocity (and Q) 
      !using values less than 1.0 cuasing reduced sed trans on WSE boundary cells
      use EXP_Global_def,  only: ncn, ncs, ncw, nce, qx, qy, cdx, cdy, num_ext_n, num_ext_s, num_ext_w, num_ext_e, etan, qxn, qyn
      USE EXP_bndcond_def, only: EXT_N,EXT_S,EXT_E,EXT_W
      use prec_def, only: ikind
      use size_def, only: ncells
      use geo_def,  only: zb, cell2cell
          
      implicit none
      integer i
      real(ikind) dep_lc

!$omp parallel do private(Dep_lc,ncn,nce)     
      do i = 1,ncells
      DEP_lc = -zb(i) + etan(i)
      ncn = cell2cell(1,i)
      nce = cell2cell(2,i)      
      cdx(i) = ((qxn(i)+qxn(nce))/2. + 1.e-10)/dep_lc
      cdy(i) = ((qyn(i)+qyn(ncn))/2.)/dep_lc 
      enddo   
!$omp end parallel do       
              
      do i=1,num_ext_N
        Dep_lc = -zb(ext_N(i,2)) + etan(ext_N(i,2))
        cdy(ext_N(i,2)) = qy(ext_N(i,2))/Dep_lc
      enddo
  
      do i=1,num_ext_S
        Dep_lc = -zb(ext_S(i,1)) + etan(ext_S(i,1))
        cdy(ext_S(i,1)) = qy(ext_S(i,2))/Dep_lc    
      enddo
  
      do i=1,num_ext_E
        Dep_lc = -zb(ext_E(i,2)) + etan(ext_E(i,2))
        cdx(ext_E(i,2)) = (qx(ext_E(i,2)) + 1.e-10)/Dep_lc
      enddo
  
      do i=1,num_ext_W
        Dep_lc = -zb(ext_W(i,1)) + etan(ext_W(i,1))
        cdx(ext_W(i,1)) = (qx(ext_W(i,2)) + 1.e-10)/Dep_lc    
      enddo
    
      end subroutine