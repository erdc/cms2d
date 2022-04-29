!*******************************************************************
    subroutine iccgstab5(nmaxiter,acoef,ap,b,phi)
! Incomplete Cholesky pre-conditioned Conjugate Gradient 
! Stabilized solver for non-symmetric matrices
!
! Author: Alex Sanchez, USACE-CHL
! Last modified: 11/26/2012
!*******************************************************************
    use size_def, only: ncells,ncellsD,nmaxfaces
    use geo_def, only: cell2cell
    use prec_def
    implicit none
    !Input/Output
    integer,intent(in) :: nmaxiter
    real(ikind),intent(in) :: acoef(nmaxfaces,ncellsD),ap(ncells),b(ncellsD)
    real(ikind),intent(inout) :: phi(ncellsD)
    !Internal Variables
    integer :: i,niter
    real(ikind) :: res0,resn,resm,resmax,bet,beto,alf,gam,om,ukreso,vres,vv
    real(ikind), dimension(ncellsD) :: res,d,zk,pk,uk,vk,reso

    resmax = 0.001
    
    !Initial Residual Vector
    res0=0.0
    !$omp parallel do private(i) reduction(+:res0)
    do i=1,ncells
      res(i)=b(i)+sum(acoef(1:4,i)*phi(cell2cell(1:4,i)))-ap(i)*phi(i) !Initial Residual Vector
      res0=res0+abs(res(i))
    enddo
    !$omp end parallel do

    !Preconditioning Matrix Diagonal
    d = 0.0
    do i=1,ncells
      d(i)=1.0/(ap(i)-d(cell2cell(4,i))*acoef(4,i)**2-d(cell2cell(3,i))*acoef(3,i)**2)
      d(i)=d(i)+1.0e-20
    enddo

    !Initialize working arrays
    !$omp parallel do private(i)
    do i=1,ncellsD
      reso(i)=res(i)
      pk(i)=0.0
      uk(i)=0.0
      zk(i)=0.0
      vk(i)=0.0  
    enddo
    !$omp end parallel do
    alf=1.0
    beto=1.0
    gam=1.0

!.....START INNER ITERATIONS
    do niter=1,nmaxiter
      !Beta and Omega  
      bet=0.0
      !$omp parallel do private(i) reduction(+:bet)
      do i=1,ncells
        bet=bet+res(i)*reso(i)
      enddo
      !$omp end parallel do
      om=bet*gam/(alf*beto+1.0e-20)
      beto=bet

      !Search vector
      !$omp parallel do private(i)
      do i=1,ncells
        pk(i)=res(i)+om*(pk(i)-alf*uk(i))
      enddo
      !$omp end parallel do
      
      !Foward substitution
      do i=1,ncells
        zk(i)=(pk(i)+acoef(4,i)*zk(cell2cell(4,i))+acoef(3,i)*zk(cell2cell(3,i)))*d(i)
      enddo
      !$omp parallel do private(i)
      do i=1,ncells
        zk(i)=zk(i)/(d(i)+1.0e-20)
      enddo
      !$omp end parallel do
      
      !Backward Substitution
      do i=ncells,1,-1
        zk(i)=(zk(i)+acoef(2,i)*zk(cell2cell(2,i))+acoef(1,i)*zk(cell2cell(1,i)))*d(i)
      enddo

      !Caculate uk = A.pk 
      !$omp parallel do private(i)
      do i=1,ncells
        uk(i)=ap(i)*zk(i)-sum(acoef(1:4,i)*zk(cell2cell(1:4,i)))  
      enddo
      !$omp end parallel do
      
      !Scalar product and gamma
      ukreso=0.0
      !$omp parallel do private(i) reduction(+:ukreso)
      do i=1,ncells
        ukreso=ukreso+uk(i)*reso(i)
      enddo
      !$omp end parallel do
      gam=bet/(ukreso+1.0e-20)

      !Update variable and calculate residual
      !$omp parallel do private(i)
      do i=1,ncells
        phi(i)=phi(i)+gam*zk(i)
        res(i)=res(i)-gam*uk(i)  
      enddo
      !$omp end parallel do
      
      !Solve for y (zk) 
      !Foward substitution
      do i=1,ncells
        zk(i)=(res(i)+acoef(4,i)*zk(cell2cell(4,i))+acoef(3,i)*zk(cell2cell(3,i)))*d(i)
      enddo
      !$omp parallel do private(i)
      do i=1,ncells
        zk(i)=zk(i)/(d(i)+1.0e-20)
      enddo
      !$omp end parallel do

      !Backward substitution
      do i=ncells,1,-1
        zk(i)=(zk(i)+acoef(2,i)*zk(cell2cell(2,i))+acoef(1,i)*zk(cell2cell(1,i)))*d(i)
      enddo

      !Calculate v = A.Y (vk = A.zk)
      !$omp parallel do private(i)
      do i=1,ncells
        vk(i)=ap(i)*zk(i)-sum(acoef(1:4,i)*zk(cell2cell(1:4,i)))  
      enddo
      !$omp end parallel do
      
      !Alpha
      vres=0.0
      vv=0.0
      !$omp parallel do private(i) reduction(+:vres) reduction(+:vv)
      do i=1,ncells
        vres=vres+vk(i)*res(i)
        vv=vv+vk(i)*vk(i)
      enddo
      !$omp end parallel do
      alf=vres/(vv+1.0e-20)

      !Update variable and residual vectors
      resn=0.0
      !$omp parallel do private(i) reduction(+:resn)
      do i=1,ncells
        phi(i)=phi(i)+alf*zk(i)
        res(i)=res(i)-alf*vk(i)
        resn=resn+abs(res(i))
      enddo
      !$omp end parallel do
      
      !Check convergence
      resm=resn/max(res0,1.0e-10)
      if(resm<resmax) return
    enddo
    
    !write(*,*) 'niter=',niter,',resn= ',resn,' resm= ',resm !'res0= ',res0,

    return
    end subroutine iccgstab5