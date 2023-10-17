!===========================================================
!  CMS Interpolation routines
!
! Contains the following:
!   interp_init - Initializes interpolation variables
!
! written by Alex Sanchez, USACE-CHL
!===========================================================

!***********************************************************
    subroutine interp_init
! Interpolation functions at cell faces   
!*********************************************************** 
    use size_def
    use geo_def
    use interp_def
    use interp_lib
    use prec_def
    implicit none
    real(ikind) :: pow
    
    !Cell-to-face Interpolation
    call interp_coef_cell2face
    
    !Cell-to-node Interpolation
    if(ncellpoly>0)then
      allocate(wc2n(nmaxcells,nnodes))
      wc2n=0.0
      nnintp = 1 !************ HARD CODED for now *****************
      select case(nnintp)
      case(1) !Inverse area
        call interp_coef_cell2node_invarea(nnodes,nncell,nmaxcells,&
           node2cell,ncells,ncellsD,areap,wc2n)
      case(2) !Inverse distance (Shepard interpolation)
        pow = 1.0 !Hard-coded for now  
        call interp_coef_cell2node_invdist(nnodes,nncell,nmaxcells,&
          node2cell,xn,yn,ncells,ncellsD,x,y,pow,wc2n)
      case(3) !Least-squares
        call interp_coef_cell2node_lstsqrs(nnodes,nncell,nmaxcells,&
          node2cell,xn,yn,ncells,ncellsD,x,y,wc2n) !Note: xc=x, yc=y for unstructured meshes
      end select
    endif
    
    return
    end subroutine interp_init

!From Chris' files - 10/26/2015
!*************************************************    
    subroutine interp2facegrad(phi,dphix,dphiy,phik)
! Interpolates pressure variables including 
! pressure correction and total water depth
! to the cell faces    
! written by Alex Sanchez, USACE
!*************************************************
    use size_def
    use geo_def, only: ncface,nxyface,kxyface,cell2cell,llec2llec,&
       idirface,dsx,dsy,dx,dy,areap,rpx,rpy,kkface
    use interp_def, only: fintp
    use flow_def, only: iwet,iextrap
    use prec_def
    implicit none
    integer :: i,j,k,nck
    real(ikind) :: fk
    real(ikind) :: phi(ncellsD),phik(ncellsD,nmaxfaces),dphix(ncellsD),dphiy(ncellsD)

!$OMP PARALLEL DO PRIVATE(i,j,k,nck,fk)  
    do i=1,ncells
      do j=1,nxyface(i)
        k=kxyface(i,j)        
!      do k=1,ncface(i)      
        nck=cell2cell(i,k)
        if(iwet(nck)==1)then
          fk=fintp(i,k)
          phik(i,k)=fk*phi(nck)+(1.0-fk)*phi(i) 
        else
          phik(i,k)=phi(i)
        endif
        phik(nck,llec2llec(i,k))=phik(i,k)
      enddo
    enddo     
!$OMP END PARALLEL DO

!$OMP PARALLEL DO PRIVATE(i,k,nck)  
    do i=1,ncells    
      dphix(i)=0.0; dphiy(i)=0.0
      do k=1,ncface(i)
        dphix(i)=dphix(i)+dsx(i,k)*phik(i,k)
        dphiy(i)=dphiy(i)+dsy(i,k)*phik(i,k)      
      enddo
      dphix(i)=dphix(i)/areap(i)
      dphiy(i)=dphiy(i)/areap(i)
    enddo
!$OMP END PARALLEL DO   
    
    return
    end subroutine interp2facegrad

