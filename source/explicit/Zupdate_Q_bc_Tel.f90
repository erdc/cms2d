!*************************************************************
      subroutine update_Q_bc_tel
!*************************************************************
      USE EXP_bndcond_def, only: qstringexp
      use EXP_TELESCOPING, only: xface_q, xface_qn, yface_q, yface_qn
      use comvarbl, only: timehrs,ramp
      use bnd_def,  only: nqstr, q_str
      use prec_def, only: ikind
      use geo_def,  only: dx,dy
      use size_def, only: ncells,ncellsD
      
      implicit none     
      integer i,j,inc,ido,iid
      real(ikind) val1,val2,tval1,tval2,fac,flow

      if(nQstr .gt. 0) then
        do i = 1,nQstr  !for each cell string
          !find out where we are in the time/value arrays
          inc = Q_str(i)%inc
          do while (timeHRS.gt.Q_str(i)%times(inc+1))
            inc = inc+1
            Q_str(i)%inc = inc
          enddo
          val1 = Q_str(i)%qcurv(inc)    
          val2 = Q_str(i)%qcurv(inc+1)    
          tval1 = Q_str(i)%times(inc)    
          tval2 = Q_str(i)%times(inc+1)
          fac = (timeHRS-tval1) / (tval2-tval1)
          flow = (val1 + fac*(val2-val1))  !/Q_str(i)%NCells
          if(QstringEXP(i)%vface) then
            IDO = Q_str(i)%NCells
            do j=1,IDO    !for each cell in string
              IID = Q_str(i)%Cells(j)    
              yface_q(IID) = QstringEXP(i)%sgn*ramp*flow  !/dx(IID)
              yface_qn(IID) = yface_q(IID)
              !write(*,*)'Y: ',flow,qy(IID),dx(IID),IID
            enddo
          else
            IDO = Q_str(i)%NCells
            do j=1,IDO    !for each cell in string
              IID = Q_str(i)%Cells(j)    
              xface_q(IID) = QstringEXP(i)%sgn*ramp*flow  !/dy(IID)
               xface_qn(IID) =  xface_q(IID)
               !write(*,*)'X: ',flow,qx(IID),dy(IID),IID           
            enddo
          endif
        enddo ! end of NQdriver
      endif  !Q_single

      end subroutine