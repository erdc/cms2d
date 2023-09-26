!***********************************************************************
    subroutine derivativeEXP()
! calculates derivatives for Explicit Solver
!***********************************************************************
    use size_def, only: ncells
    use flow_def, only: dux, duy, dvx, dvy, u, v
    use geo_def,  only: dx,dy,cell2cell
    
    implicit none
    integer ncn,nce,ncw,ncs,i
    real dxT,dyT

    do i=1,ncells
      ncn=cell2cell(1,i)
      nce=cell2cell(2,i)
      ncs=cell2cell(3,i)
      ncw=cell2cell(4,i) 
      dux(i) = (u(nce)-u(i))/dx(i)
      dvy(i) = (v(ncn)-v(i))/dy(i)
      dyT = 0.25*dy(ncn)+0.5*dy(i)+0.25*dy(ncs)
      duy(i) = (u(ncn)-u(ncs))/dyT
      dxT = 0.25*dx(nce)+0.5*dx(i)+0.25*dx(ncw)    
      dvx(i) = (v(nce)-v(ncw))/dxT         
    enddo
    
    return
    end subroutine
        
!***********************************************************************
    subroutine UpdateDummyCells(var)
! calculates derivatives for Explicit Solver
!***********************************************************************
    use size_def, only: ncells,ncellsD
    use prec_def, only: ikind
    use EXP_Global_def, only: linktodummies
    
    implicit none
    real(ikind) var(ncellsD)
    integer ii,i,jj
    
    ii=0
    do i=ncells+1,ncellsD
      ii=ii+1
      jj = linktodummies(ii)
      var(i) = var(jj)
    enddo
    
    return
    end subroutine