!***********************************************************************
    subroutine pbicgstab(nmaxiter,acoef,ap,su,phi)
!   this is bicgstab solver       
!***********************************************************************      
    use size_def, only: ncells,ncellsD,nmaxfaces
    use geo_def, only: cell2cell,ncface
    use flow_def, only: iwet
    use const_def, only: small
    use solv_def 
    use prec_def
    !Intput/Output
    integer,intent(in) :: nmaxiter
    real(ikind),intent(in) :: acoef(nmaxfaces,ncellsD),ap(ncells),su(ncellsD)
    real(ikind),intent(inout) :: phi(ncellsD)
    !Internal Variables
    integer :: ierr,i,iter      
    real(ikind) :: pk(ncells),r(ncells),uk(ncells),vk(ncells),r0(ncells),zzk(ncellsD)
    real(ikind) :: omg,bet,bet0,sum1,sum2,ukr0,res,res0,term,alf,gam,adbk

    res0=0
!$OMP PARALLEL DO PRIVATE(i,term) REDUCTION(+:res0)
    do i=1,ncells
      term=sum(acoef(1:ncface(i),i)*phi(cell2cell(1:ncface(i),i)))+su(i)-ap(i)*phi(i)
      r(i)=term+small
      r0(i)=r(i)
      res0=res0+abs(r0(i))
      pk(i)=0.0
      uk(i)=0.0
      vk(i)=0.0
    enddo
!$OMP END PARALLEL DO
!-------
    call ilu_cr(acoef,ap)
!    call ilu0(ncells)
    call ilutp(ncells,ierr)
!------
    zzk=0.0
    alf=1.0; bet0=1.0; gam=1.0

    !..ITERATION LOOPS BEGIN
    do iter=1,nmaxiter
       !bet=sum(r(1:ncells)*r0(1:ncells))
       bet=adbk(ncells,r,r0)
       omg=bet*gam/(alf*bet0+small)
       bet0=bet

      !..CALCULATE pk
!$OMP PARALLEL DO PRIVATE(i)          
      do i=1,ncells
        pk(i)=r(i)+omg*(pk(i)-alf*uk(i))
      enddo
!$OMP END PARALLEL DO        

      call lusol(ncells, pk(1:ncells), zzk(1:ncells))

      ukr0=0.0
!$OMP PARALLEL DO PRIVATE(i) REDUCTION(+:ukr0)
      do i=1,ncells
         uk(i)=ap(i)*zzk(i)-sum(acoef(1:ncface(i),i)*zzk(cell2cell(1:ncface(i),i)))
         ukr0=ukr0+uk(i)*r0(i)
      enddo
!$OMP END PARALLEL DO           
      gam=bet/ukr0

      !..UPDATE (phi) AND CALCULATE R 
!$OMP PARALLEL DO PRIVATE(i)     
      do i=1,ncells
        phi(i)=phi(i)+gam*zzk(i)
        r(i)=r(i)-gam*uk(i)
      enddo
!$OMP END PARALLEL DO        

      call lusol(ncells, r(1:ncells), zzk(1:ncells))        

      !..CALCULATE vk and ALPHA (alf)
      sum1=0.0; sum2=0.0
!$OMP PARALLEL DO PRIVATE(i) REDUCTION(+:sum1), REDUCTION(+:sum2) 
      do i=1,ncells
        vk(i)=ap(i)*zzk(i)-sum(acoef(1:ncface(i),i)*zzk(cell2cell(1:ncface(i),i))) 
        sum1=sum1+vk(i)*r(i)
        sum2=sum2+vk(i)*vk(i)  
      enddo     
!$OMP END PARALLEL DO        
      alf=sum1/max(sum2,small)

      !..UPDATE VARIABLE (phi) AND RESIDUAL (r) VECTORS
      res=0.0
!$OMP PARALLEL DO PRIVATE(i) REDUCTION(+:res)     
      do i=1,ncells
        phi(i)=phi(i)+alf*zzk(i)
        r(i)=r(i)-alf*vk(i)
        res=res+abs(r(i))
      enddo
!$OMP END PARALLEL DO        
      if(res/res0.lt.0.001) return 
    enddo
    
    return
    end subroutine pbicgstab
    
!*******************************************************************
    subroutine iccg5(nmaxiter,acoef,ap,b,phi)
! Incomplete Cholesky preconditioned Conjugate Gradient 
! solver for 5-point stencil symmetric matrices 
!
! Desciption:
!   Incomplete Cholesky pre-conditioned Conjugate Gradient 
!   solver for 5-point stencil symmetric matrices.
!   The method is substantially faster than the SIP method
!   but is only applicable to symmetric matrices such as
!   the pressure correction equation. The solver has been
!   parallelized using OpenMP except for two loops.
!
! Input:
!   nmaxiter - maximum number of inner loop iterations
!   acoef - coefficients for north, east, south, and west neighbors
!   ap - coefficient for diagonal
!   b - R.H.S
!   phi - initial value for solution vector
!
! Output:
!   phi - final value for solution vector
!
! References:
!   Ferziger, J.H., and Peric, M. (2002). Computational Methods for 
!      Fluid Dynamics, Springer, 3rd Edition, 423 p. ISBN 3-540-42074-6
!
! Author: Alex Sanchez, USACE-CHL
! Last modified: 11/26/2012
!*******************************************************************
#include "CMS_cpp.h"
    use size_def, only: ncells,ncellsD,nmaxfaces
    use geo_def, only: cell2cell
    use solv_def, only: iconv
#ifdef DIAG_MODE
    use diag_lib
#endif
    use prec_def
    implicit none
    !Input/Output
    integer,intent(in) :: nmaxiter
    real(ikind),intent(in) :: acoef(nmaxfaces,ncellsD),ap(ncells),b(ncellsD)
    real(ikind),intent(inout) :: phi(ncellsD)
    !Internal Variables
    integer :: i,niter,nitermin
    real(ikind) :: res0,resn,resm,resmax,sk,s0,s2,beta,alpha,dr
    real(ikind), dimension(ncellsD) :: res,d,zk,pk
#ifdef DIAG_MODE
    logical :: isnankind
#endif

    nitermin = nmaxiter/2
    resmax = 0.001
    
    !Initialize Ghost/Dummy Cells
    !$omp parallel do private(i)
    do i=1,ncellsD
      d(i)=0.0; zk(i)=0.0; pk(i)=0.0
    enddo
    !$omp end parallel do
    
    !Initial Residual Vector
    res0=0.0
    !$omp parallel do private(i) reduction(+:res0)
    do i=1,ncells
      res(i)=b(i)+sum(acoef(1:4,i)*phi(cell2cell(1:4,i)))-ap(i)*phi(i) !Initial Residual Vector
#ifdef DIAG_MODE
      if(isnankind(res(i)) .or. abs(res(i))>1.0e6)then
        call diag_print_error('Problem in ICCG solver: Residual equal to NaN')
      endif
#endif      
      res0=res0+abs(res(i))
    enddo
    !$omp end parallel do

    !Preconditioning Matrix Diagonal
    do i=1,ncells
      !d(i)=1.0/(ap(i)-d(cell2cell(4,i))*acoef(4,i)**2-d(cell2cell(3,i))*acoef(3,i)**2)
      dr=ap(i)-d(cell2cell(4,i))*acoef(4,i)**2-d(cell2cell(3,i))*acoef(3,i)**2
      d(i)=1.0/max(dr,1.0e-3)
      !d(i)=max(d(i),1.0e-10)
#ifdef DIAG_MODE
      if(isnankind(d(i)))then
        call diag_print_error('Problem in ICCG solver: Preconditioner diagonal equal to NaN')
      endif
#endif
    enddo

    s0=1.0e20 !Initialize
    do niter=1,nmaxiter !Inner iteration loop
      !Forward substitution
      do i=1,ncells
        zk(i)=(res(i)+acoef(4,i)*zk(cell2cell(4,i))+acoef(3,i)*zk(cell2cell(3,i)))*d(i)
#ifdef DIAG_MODE
        if(isnankind(zk(i)) .or. zk(i)>1.0e10)then
          call diag_print_error('Problem in ICCG solver: NaN found in forward substitution')
        endif
#endif        
      enddo
      !$omp parallel do private(i)
      do i=1,ncells
        zk(i)=zk(i)/d(i)
#ifdef DIAG_MODE
        if(isnankind(zk(i)) .or. zk(i)>1.0e10)then
          call diag_print_error('Problem in ICCG solver: Problem calculating zk(i)')
        endif
#endif
      enddo
      !$omp end parallel do

      !Backward Substitution
      do i=ncells,1,-1
        zk(i)=(zk(i)+acoef(2,i)*zk(cell2cell(2,i))+acoef(1,i)*zk(cell2cell(1,i)))*d(i)
        !zk(i)=max(min(zk(i),1.0e6),-1.0e6)
#ifdef DIAG_MODE
        if(isnankind(zk(i)) .or. zk(i)>1.0e30)then
          call diag_print_error('Problem in ICCG solver: NaN found in backward substitution')
        endif
#endif
      enddo
      
      sk=0.0; s2=0.0; resn=0.0  !Initialize before parallel region
      !$omp parallel do private(i) reduction(+:sk)
      do i=1,ncells
        sk=sk+res(i)*zk(i) !Auxiliary vector
      enddo
      !$omp end parallel do
      
      beta=sk/max(s0,1.0e-10)
      
      !$omp parallel do private(i)
      do i=1,ncells
        pk(i)=zk(i)+beta*pk(i) !New search vector
#ifdef DIAG_MODE
        if(isnankind(pk(i)))then
          call diag_print_error('Problem in ICCG solver: Search vector equal to NaN')
        endif
#endif
      enddo
      !$omp end parallel do
      
      !$omp parallel do private(i) reduction(+:s2)
      do i=1,ncells 
        zk(i)=ap(i)*pk(i)-sum(acoef(1:4,i)*pk(cell2cell(1:4,i)))
#ifdef DIAG_MODE
        if(isnankind(zk(i)))then
          call diag_print_error('Problem in ICCG solver: zk equal to NaN')
        endif
#endif
        s2=s2+pk(i)*zk(i)
      enddo
      !$omp end parallel do
      
      alpha=sk/max(s2,1.0e-10)
      
      !$omp parallel do private(i) reduction(+:resn)
      do i=1,ncells
        phi(i)=phi(i)+alpha*pk(i) !Update variable
#ifdef DIAG_MODE
        if(isnankind(phi(i)))then
          call diag_print_error('Problem in ICCG solver: Variable update equal to NaN')
        endif
#endif
        res(i)=res(i)-alpha*zk(i) !New residual
        resn=resn+abs(res(i))
      enddo
      !$omp end parallel do
      
      s0=sk

      !Check convergence
      resm=resn/max(res0,1.0e-10)
      !write(*,*) 'niter =',niter,',  resm = ',resm
      if(niter>nitermin .and. resm<resmax)then
        !write(*,*) 'Exit ICCG: ',niter,resm  
        return 
      elseif(resm>10.0)then
        iconv=0
        return
      endif
    enddo
      
    return
    end subroutine iccg5    

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
    
!***********************************************************************       
    subroutine pgmres(nmaxiter,acoef,ap,su,phi)   
!***********************************************************************
    use size_def, only: ncells,ncellsD,nmaxfaces
    use geo_def, only: cell2cell,ncface
    use const_def, only: small
    use solv_def 
    use prec_def
    integer,intent(in) :: nmaxiter
    real(ikind),intent(in) :: acoef(nmaxfaces,ncellsD),ap(ncells),su(ncellsD)
    real(ikind),intent(inout) :: phi(ncellsD)
    !Internal Variables
    integer :: im,ierr,itr,i,ii,j,i1,k,kk,k1
    integer, parameter :: kmax=3       
    real(ikind) :: vv(ncells,kmax+1),rhs(ncellsD),hh(kmax+1,kmax),c(kmax),s(kmax)
    real(ikind) :: eps1,t,rs1(kmax+1),term,eps11,gam,ro,adak,adbk
        
    im=kmax
    eps1=0.001
!$OMP PARALLEL DO PRIVATE(ii)
    do ii=1,ncellsD
      rhs(ii)=0.0
    enddo  
!$OMP END PARALLEL DO    

    !call ilu_cr(acoef,ap)
    call coef2csr(acoef,ap)
    call ilutp(ncells,ierr)
!    call ilu0(ncells)
    
    do itr=1,nmaxiter
      !--- Initial residual vector --------------------------  
      !call ceaambk(ncells,aa_matrix,ja,ia,phi,vv(:,1))
      !vv(1:ncells,1) = su(1:ncells) - vv(1:ncells,1)
      !ro = sqrt(adak(ncells,vv(1:ncells,1)))
      
      ro=0.0
!!$OMP PARALLEL DO PRIVATE(ii,term) REDUCTION(+:ro)   
      do ii=1,ncells
        term=sum(acoef(1:ncface(ii),ii)*phi(cell2cell(1:ncface(ii),ii)))
        vv(ii,1)=term+su(ii)-ap(ii)*phi(ii)
!        if (isnan(vv(ii,1))) then
!          stop '"vv(ii,1)" is a NaN'
!        endif
        ro=ro+vv(ii,1)*vv(ii,1)
      enddo 
!!$OMP END PARALLEL DO 
      ro=sqrt(ro+small)     
       
      if(itr.eq.1) eps11=eps1*ro
      t=1.0/(ro+small)
      call aesmak(ncells,t,vv(1:ncells,1))
      
      !---Initialize first term of RHS of Hessenberg system
      rs1(1)=ro
!       rs1(2:kmax+1)=0.0
!       hh(2:kmax+1,1:kmax+1)=0.0

      i=0
 4    i=i+1
      i1=i+1

      call lusol(ncells,vv(1:ncells,i),rhs(1:ncells))       
!$OMP PARALLEL DO PRIVATE(ii,term)     
      do ii=1,ncells
        term=sum(acoef(1:ncface(ii),ii)*rhs(cell2cell(1:ncface(ii),ii)))
        vv(ii,i1)=ap(ii)*rhs(ii)-term
      enddo
!$OMP END PARALLEL DO
       !call ceaambk(ncells,aa_matrix,ja,ia,rhs,vv(:,i1))

      !---Modified Gram-Schmidt ---------------------------------
      do j=1,i
        hh(j,i)=adbk(ncells,vv(1:ncells,j),vv(1:ncells,i1)) !hh(j,i)=sum(vv(1:ncells,j)*vv(1:ncells,i1))
        call aesmbpak(ncells,-hh(j,i),vv(1:ncells,j),vv(1:ncells,i1)) !vv(1:ncells,i1)=-hh(j,i)*vv(1:ncells,j)+vv(1:ncells,i1)
      enddo
      t=sqrt(adak(ncells,vv(1:ncells,i1))+small)
      hh(i1,i)=t
      t=1.0/max(t,small)
      call aesmak(ncells,t,vv(1:ncells,i1))
      
      !---perfrom previous transformations  on i-th column of h
      if(i>1)then
        do k=2,i
          k1=k-1
          t=hh(k1,i)
          hh(k1,i)=c(k1)*t+s(k1)*hh(k,i)
          hh(k,i)=-s(k1)*t+c(k1)*hh(k,i)
        enddo
      endif
      gam=sqrt(hh(i,i)**2+hh(i1,i)**2)

      !---if gamma is zero then any small value will do...
      if(gam==0.0) gam=small

      !----get next plane rotation
      c(i)=hh(i,i)/gam
      s(i)=hh(i1,i)/gam
      rs1(i1)=-s(i)*rs1(i)
      rs1(i)=c(i)*rs1(i)

      !----detrermine residual norm and test for convergence
      hh(i,i)=c(i)*hh(i,i)+s(i)*hh(i1,i)
      ro=abs(rs1(i1))
      if(i.lt.im.and.(ro>eps11)) goto 4

      !---now compute solution. first solve upper triangular system.
      rs1(i)=rs1(i)/hh(i,i)
      do ii=2,i
        k=i-ii+1
        k1=k+1
        t=rs1(k)
        do j=k1,i
           t=t-hh(k,j)*rs1(j)
        enddo
        rs1(k)=t/hh(k,k)
      enddo

       !----form linear combination of v(*,i)'s to get solution
!$OMP PARALLEL DO PRIVATE(kk) 
       do kk=1,ncells
          rhs(kk)=sum(vv(kk,1:i)*rs1(1:i))
       enddo
!$OMP END PARALLEL DO 

       !---call preconditioner.
       call lusol(ncells,rhs(1:ncells),rhs(1:ncells))
       call aeapbk(ncells,phi(1:ncells),rhs(1:ncells))
       if (ro.le.eps11) return

!       !---else compute residual vector and continue..
!       do j=1,i
!          jj=i1-j+1
!          rs1(jj-1)=-s(jj-1)*rs1(jj)
!          rs1(jj)=c(jj-1)*rs1(jj)
!       enddo
!       do j=1,i1
!          t=rs1(j)
!          if(j==1) t=t-1.0
!          vv(1,j)=vv(1,j)+t*vv(1,j)
!       enddo

    enddo

    return
    end subroutine pgmres

!*****************************************************************************
    subroutine sip5(nmaxiter,acoef,ap,b,phi)
! Strongly Implicit Procedure (SIP) solver for a 5 cell stencil.
!
! Desciption:
!   ILU solver after Stone (1968). Partially parallelized with OpenMP.
!   Only for structured grids.
!   The index for should be such that (phi(ij),ij=1,ncells)
!   follows the same ordering as (((phi(i,j),i=1,nx),j=1,ny)) would
!   using a structured data structure.
!   The ordering used here allows compuational saving since inactive
!   regions are skipped and only one do loop is necessary for each section.
!
! Input:
!   nmaxiter - maximum number of inner loop iterations
!   acoef - coefficients for north, east, south, and west neighbors
!   ap - coefficient for diagonal
!   b - R.H.S
!   phi - initial value for solution vector
!
! Output:
!   phi - final value for solution vector
! 
! References:
!   Stone, H.L. (1968) Iterative Solution of Implicit Approximations of
!      Multidimensional Partial Differential Equations,
!      SIAM Journal of Numerical Analysis, 5(3), 530–538. doi:10.1137/0705044
!   Ferziger, J.H., and Peric, M. (2002). Computational Methods for 
!      Fluid Dynamics, Springer, 3rd Edition, 423 p. ISBN 3-540-42074-6
!
! Author: Alex Sanchez, USACE-CHL
! Last modified: 01/24/2013
!*****************************************************************************
#include "CMS_cpp.h"
    use size_def, only: ncells,ncellsD,nmaxfaces
    use geo_def, only: cell2cell
    use solv_def, only: iconv
    use prec_def
#ifdef DIAG_MODE
    use diag_def
    use diag_lib
#endif   
    implicit none
    !Input/Output
    integer,intent(in) :: nmaxiter
    real(ikind),intent(in) :: acoef(nmaxfaces,ncellsD),ap(ncells),b(ncellsD)
    real(ikind),intent(inout) :: phi(ncellsD)
    !Internal Variables
    integer :: i,niter,nitermin,ncs,ncw
    real(ikind) :: alpha,resmax,ressum,resn,res0,rsm,term1,term2,realncells
    real(ikind) :: ressumt
    real(ikind), dimension(ncellsD) :: Lw,Ls,Lpr,r,Un,Ue,res
#ifdef DIAG_MODE
    logical :: isnankind
#endif   

    !Using alpha=0.0 produces the standard ILU solver
    !For optimum value of alpha SIP is about 6 times faster than ILu
    !The optimum value usually lies between 0.92 - 0.96 Ferziger and Peric (2002).
    !For alpha values larger than the optimum value, the method does not converge.
    !A safe value for alpha is 0.92 which provides about 5 times 
    !the speed of the standard ILU.
    alpha = 0.92
    resmax = 0.001
    nitermin = nmaxiter/2
    realncells = real(ncells,kind=ikind)
      
    !Initialize vectors
!$OMP PARALLEL DO PRIVATE(i) IF(ncellsD-ncells>1000)
    do i=ncells+1,ncellsD
      Lw(i) = 0.0; Ls(i) = 0.0; Lpr(i) = 0.0
      Un(i) = 0.0; Ue(i) = 0.0; r(i) = 0.0
    enddo
!$OMP END PARALLEL DO
    
    !=== Factorise the equations =======================================
    do i=1,ncells
      ncs = cell2cell(3,i)
      ncw = cell2cell(4,i)
      Lw(i) = -acoef(3,i)/(1.0+alpha*(Un(ncs)))
      Ls(i) = -acoef(4,i)/(1.0+alpha*(Ue(ncw)))
      term1 = alpha*Lw(i)*Un(ncs)
      term2 = alpha*Ls(i)*Ue(ncw)
      Lpr(i) = 1.0/(ap(i)+term1+term2-(Ls(i)*Un(ncw)+Lw(i)*Ue(ncs)))
      Un(i) = -(acoef(2,i)+term1)*Lpr(i)
      Ue(i) = -(acoef(1,i)+term2)*Lpr(i)
#ifdef DIAG_MODE
      if(abs(Lw(i))>1.0e20 .or. abs(Ls(i))>1.0e20 .or. abs(Lpr(i))>1.0e20 .or. &
         abs(Un(i))>1.0e20 .or. abs(Ue(i))>1.0e20)then
        write(msg2,*) 'i = ',i
        write(msg2,*) 'Lpr(i) = ',Lpr(i)
        write(msg3,*) 'Lw(i) = ',Lw(i),', Ls(i) = ',Ls(i)
        write(msg4,*) 'Un(i) = ',Un(i),', Ue(i) = ',Ue(i)
        call diag_print_warning('Problem in SIP factorization',msg2,msg3,msg4)  
         endif         
#endif
    enddo

    !=== Inner iteration loop ================================
    do niter=1,nmaxiter
      !--- Residual ------------------------------------------------------
      ressum = 0.0; ressumt = 0.0
!$OMP PARALLEL FIRSTPRIVATE(ressumt)
!$OMP DO PRIVATE(i)
      do i=1,ncells
        res(i) = b(i) + sum(acoef(1:4,i)*phi(cell2cell(1:4,i))) - ap(i)*phi(i) 
#ifdef DIAG_MODE
        if(isnankind(res(i)) .or. abs(res(i))>1.0e20)then
          write(msg,*)  'niter = ',niter,', resn = ',resn,', res0 = ',res0
          write(msg2,*) 'i = ',i,', b(i) = ',b(i),', ap(i) = ',ap(i),', phi(i) = ',phi(i)
          write(msg3,*) 'acoef(1:4,i) = ',acoef(1:4,i)
          write(msg4,*) 'phi(cell2cell(1:4,i))) = ',phi(cell2cell(1:4,i))
          call diag_print_error('Problem calculating SIP residuals',&
            msg2,'res(i) = NaN',msg3,msg4)
        endif
#endif
        ressumt = ressumt + res(i)**2  
      enddo
!$OMP END DO
!$OMP CRITICAL
      ressum = ressum + ressumt
!$OMP END CRITICAL
!$OMP END PARALLEL    
      resn = sqrt(ressum/realncells)
      if(niter==1) res0=resn
      !--- Forward Sweep -------------------------------------------------
      do i=1,ncells
        r(i) = Lpr(i)*(res(i)-Lw(i)*r(cell2cell(3,i))-Ls(i)*r(cell2cell(4,i)))
#ifdef DIAG_MODE
        if(isnankind(r(i)) .or. abs(r(i))>1.0e20)then
          write(msg,*)  'niter = ',niter,', resn = ',resn,', res0 = ',res0
          write(msg2,*) 'i = ',i,'r(i) = ',r(i)
          write(msg3,*) 'Lpr(i) = ',Lpr(i),', Lw(i) = ',Lw(i),', Ls(i) = ',Ls(i)
          write(msg4,*) 'r(cell2cell(3,i)) = ',r(cell2cell(3,i))
          write(msg5,*) 'r(cell2cell(4,i)) = ',r(cell2cell(4,i))
          call diag_print_error('Problem in SIP foward sweep',&
            msg,msg2,msg3,msg4,msg5)
        endif
#endif         
      enddo
      !--- Backwards Sweep ---------------------------------------
      do i=ncells,1,-1
        r(i) = r(i) - Ue(i)*r(cell2cell(1,i)) - Un(i)*r(cell2cell(2,i))
#ifdef DIAG_MODE
        if(isnankind(r(i)) .or. abs(r(i))>1.0e20)then
          write(msg,*)  'niter = ',niter,', resn = ',resn,', res0 = ',res0  
          write(msg2,*) 'i = ',i,', Ue(i) = ',Ue(i),', Un(i) = ',Un(i)
          write(msg3,*) 'r(cell2cell(1,i)) = ',r(cell2cell(1,i))
          write(msg4,*) 'r(cell2cell(2,i)) = ',r(cell2cell(2,i))
          call diag_print_error('Problem in SIP backward sweep',&
            msg,msg2,' r(i) = NaN',msg3,msg4)
        endif
#endif
      enddo
      !--- Correct Solution ------------------
!$OMP PARALLEL DO PRIVATE(i)
      do i=1,ncells
        phi(i) = phi(i) + r(i)
      enddo
!$OMP END PARALLEL DO  
      !Print residual, and check for convergance
      rsm = resn/max(res0,1.0e-6)
      !write(*,*) 'SIP: ',niter,rsm  
      if(niter>nitermin .and. rsm<resmax)then
        !write(*,*) 'Exit SIP: ',niter,rsm  
        return
      elseif(rsm>10.0)then
        iconv = 0
        return  
      endif
    enddo
    
    return
    end subroutine sip5

!********************************************************************************
    subroutine gauss_seidel(nmaxiter,acoef,ap,su,phi)
! Gauss-Seidel solver
!
! Description:
!   Standard Gauss-Seidel solver with OpenMP parallelization.
!  
! written by Weiming Wu, NCCHE, Oct. 2008
! modified by Alex Sanchez, USACE-CHL
!*********************************************************************************
    use size_def, only: ncells,ncellsD,nmaxfaces
    use geo_def, only: cell2cell,ncface
    use prec_def
    implicit none
    !Input/Output
    integer,intent(in) :: nmaxiter
    real(ikind),intent(in) :: acoef(nmaxfaces,ncellsD),ap(ncells),su(ncellsD)
    real(ikind),intent(inout) :: phi(ncellsD)
    !Internal Variables
    integer :: i,iter
    real(ikind) :: term

!   inner iteration loop
!$OMP PARALLEL
    do iter=1,nmaxiter
!$OMP DO PRIVATE(i,term)
       do i=1,ncells
          term=sum(acoef(1:ncface(i),i)*phi(cell2cell(1:ncface(i),i))) !Must be separate from line below
          phi(i)=(term+su(i))/ap(i)
       enddo
!$OMP END DO
    enddo
!$OMP END PARALLEL

    return
    end subroutine gauss_seidel

!***********************************************************************
    subroutine gauss_seidel_SOR(nmaxiter,relaxsor,acoef,ap,su,phi)
! Gauss-Seidel solver with Successive-Over-Relaxation
!
! Description:
!   Standard Gauss-Seidel solver with OpenMP parallelization.

!  written by Weiming Wu, NCCHE, Oct. 2008
!  modified by Alex Sanchez, USACE-CHL
!***********************************************************************
    use size_def, only: ncells,ncellsD,nmaxfaces
    use geo_def, only: cell2cell,ncface
    use prec_def
    implicit none
    !Input/Output
    integer,intent(in) :: nmaxiter
    real(ikind),intent(in) :: relaxsor
    real(ikind),intent(in) :: acoef(nmaxfaces,ncellsD),ap(ncells),su(ncellsD)
    real(ikind),intent(inout) :: phi(ncellsD)
    !Internal Variables    
    integer :: iter,i
    real(ikind) :: comrelaxsor,term
    
    comrelaxsor=1.0-relaxsor
!inner iteration loop
!$OMP PARALLEL
    do iter=1,nmaxiter
!$OMP DO PRIVATE(i,term)
       do i=1,ncells
          term=sum(acoef(1:ncface(i),i)*phi(cell2cell(1:ncface(i),i))) !Must be separate from line below
          phi(i)=relaxsor*(term+su(i))/ap(i)+comrelaxsor*phi(i)
      enddo
!$OMP END DO
    enddo    
!$OMP END PARALLEL

    return
    end subroutine gauss_seidel_SOR

!********************************************************************************
    subroutine solv_gs(nmaxiter,acoef,ap,su,phi)
! Gauss-Seidel solver 
! written by Weiming Wu, NCCHE, Oct. 2008
! modified by Alex Sanchez, USACE-CHL
!*********************************************************************************
    use size_def, only: ncells,ncellsD,nmaxfaces
    use geo_def, only: cell2cell,ncface
    use prec_def
    implicit none
    !Input/Output
    integer,intent(in) :: nmaxiter
    real(ikind),intent(in) :: acoef(nmaxfaces,ncellsD),ap(ncells),su(ncellsD)
    real(ikind),intent(inout) :: phi(ncellsD)
    !Internal Variables
    integer :: i,ii,iter,m
    real(ikind) :: term

!   inner iteration loop
    m=mod(ncells,5)
!$OMP PARALLEL
    do iter=1,nmaxiter
!$OMP MASTER        
       do i=1,m
          term=sum(acoef(1:ncface(i),i)*phi(cell2cell(1:ncface(i),i))) !Must be separate from line below
          phi(i)=(term+su(i))/ap(i)
       enddo
!$OMP END MASTER       
!$OMP DO PRIVATE(i,ii,term)
       do i=m+1,ncells,5
          term=sum(acoef(1:ncface(i),i)*phi(cell2cell(1:ncface(i),i))) !Must be separate from line below
          phi(i)=(term+su(i))/ap(i)
          ii=i+1
          term=sum(acoef(1:ncface(ii),ii)*phi(cell2cell(1:ncface(ii),ii))) !Must be separate from line below
          phi(ii)=(term+su(ii))/ap(ii)
          ii=i+2
          term=sum(acoef(1:ncface(ii),ii)*phi(cell2cell(1:ncface(ii),ii))) !Must be separate from line below
          phi(ii)=(term+su(ii))/ap(ii)
          ii=i+3
          term=sum(acoef(1:ncface(ii),ii)*phi(cell2cell(1:ncface(ii),ii))) !Must be separate from line below
          phi(ii)=(term+su(ii))/ap(ii)
          ii=i+4
          term=sum(acoef(1:ncface(ii),ii)*phi(cell2cell(1:ncface(ii),ii))) !Must be separate from line below
          phi(ii)=(term+su(ii))/ap(ii)
      enddo
!$OMP END DO
    enddo
!$OMP END PARALLEL

    return
    end subroutine solv_gs
    
!***********************************************************************
    subroutine solv_sor(nmaxiter,relaxsor,acoef,ap,su,phi)
! Gauss-Seidel solver with Successive-Over-Relaxation
!  written by Weiming Wu, NCCHE, Oct. 2008
!  modified by Alex Sanchez, USACE-CHL
!***********************************************************************
    use size_def, only: ncells,ncellsD,nmaxfaces
    use geo_def, only: cell2cell,ncface
    use prec_def
    implicit none
    !Input/Output
    integer,intent(in) :: nmaxiter
    real(ikind),intent(in) :: relaxsor
    real(ikind),intent(in) :: acoef(nmaxfaces,ncellsD),ap(ncells),su(ncellsD)
    real(ikind),intent(inout) :: phi(ncellsD)
    !Internal Variables    
    integer :: iter,i,ii,m
    real(ikind) :: comrelaxsor,term
    
    comrelaxsor=1.0-relaxsor
    m=mod(ncells,5)
!inner iteration loop
!$OMP PARALLEL 
    do iter=1,nmaxiter
!$OMP MASTER        
       do i=1,m
          term=sum(acoef(1:ncface(i),i)*phi(cell2cell(1:ncface(i),i))) !Must be separate from line below
          phi(i)=relaxsor*(term+su(i))/ap(i)+comrelaxsor*phi(i)
       enddo
!$OMP END MASTER
!$OMP DO PRIVATE(i,ii,term)
       do i=m+1,ncells,5
          term=sum(acoef(1:ncface(i),i)*phi(cell2cell(1:ncface(i),i))) !Must be separate from line below
          phi(i)=relaxsor*(term+su(i))/ap(i)+comrelaxsor*phi(i)
          ii=i+1
          term=sum(acoef(1:ncface(ii),ii)*phi(cell2cell(1:ncface(ii),ii))) !Must be separate from line below
          phi(ii)=relaxsor*(term+su(ii))/ap(ii)+comrelaxsor*phi(ii)
          ii=i+2
          term=sum(acoef(1:ncface(ii),ii)*phi(cell2cell(1:ncface(ii),ii))) !Must be separate from line below
          phi(ii)=relaxsor*(term+su(ii))/ap(ii)+comrelaxsor*phi(ii)
          ii=i+3
          term=sum(acoef(1:ncface(ii),ii)*phi(cell2cell(1:ncface(ii),ii))) !Must be separate from line below
          phi(ii)=relaxsor*(term+su(ii))/ap(ii)+comrelaxsor*phi(ii)
          ii=i+4
          term=sum(acoef(1:ncface(ii),ii)*phi(cell2cell(1:ncface(ii),ii))) !Must be separate from line below
          phi(ii)=relaxsor*(term+su(ii))/ap(ii)+comrelaxsor*phi(ii)
      enddo
!$OMP END DO
    enddo
!$OMP END PARALLEL

    return
    end subroutine solv_sor

!***********************************************************************
    subroutine solv_sgs(nmaxiter,acoef,ap,su,phi)
! Symmetric Gauss-Seidel solver
! written by Alex Sanchez, USACE-CHL
!***********************************************************************
    use size_def, only: ncells,ncellsD,nmaxfaces
    use geo_def, only: cell2cell,ncface
    use prec_def
    implicit none
    !Input/Output
    integer,intent(in) :: nmaxiter
    real(ikind),intent(in) :: acoef(nmaxfaces,ncellsD),ap(ncells),su(ncellsD)
    real(ikind),intent(inout) :: phi(ncellsD)
    !Internal Variables
    integer :: iter,i
    real(ikind) :: term
    
!inner iteration loop
!$OMP PARALLEL
    do iter=1,nmaxiter/2
!$OMP DO PRIVATE(i,term)
       do i=1,ncells
          term=sum(acoef(1:ncface(i),i)*phi(cell2cell(1:ncface(i),i))) !Must be separate from line below
          phi(i)=(term+su(i))/ap(i)
       enddo
!$OMP END DO
!$OMP DO PRIVATE(i,term)
       do i=ncells,1,-1
          term=sum(acoef(1:ncface(i),i)*phi(cell2cell(1:ncface(i),i))) !Must be separate from line below
          phi(i)=(term+su(i))/ap(i)
       enddo
!$OMP END DO
    enddo
!$OMP END PARALLEL

    return
    end subroutine solv_sgs 
    
!***********************************************************************
    subroutine solv_ssor(nmaxiter,relaxsor,acoef,ap,su,phi)
! Symmetric Gauss-Seidel solver with Successive-Over-Relaxation
! written by Alex Sanchez, USACE-CHL
!***********************************************************************
    use size_def, only: ncells,ncellsD,nmaxfaces
    use geo_def, only: cell2cell,ncface
    use prec_def
    implicit none
    !Input/Output
    integer,intent(in) :: nmaxiter
    real(ikind),intent(in) :: relaxsor
    real(ikind),intent(in) :: acoef(nmaxfaces,ncellsD),ap(ncells),su(ncellsD)
    real(ikind),intent(inout) :: phi(ncellsD)
    !Internal Variables
    integer :: iter,i
    real(ikind) :: comrelaxsor,term
    
    comrelaxsor=1.0-relaxsor
    
!inner iteration loop
!$OMP PARALLEL
    do iter=1,nmaxiter/2
!$OMP DO PRIVATE(i,term)
       do i=1,ncells
          term=sum(acoef(1:ncface(i),i)*phi(cell2cell(1:ncface(i),i))) !Must be separate from line below
          phi(i)=relaxsor*(term+su(i))/ap(i)+comrelaxsor*phi(i)
      enddo
!$OMP END DO
!$OMP DO PRIVATE(i,term)
       do i=ncells,1,-1
          term=sum(acoef(1:ncface(i),i)*phi(cell2cell(1:ncface(i),i))) !Must be separate from line below
          phi(i)=relaxsor*(term+su(i))/ap(i)+comrelaxsor*phi(i)
       enddo
!$OMP END DO
    enddo
!$OMP END PARALLEL

    return
    end subroutine solv_ssor
    
!***********************************************************************
    subroutine solv_ssorac(nmaxiter,acoef,ap,su,phi)
! Symmetric Gauss-Seidel solver with Successive-Over-Relaxation 
! and ACceleration
! written by Alex Sanchez, USACE-CHL
!***********************************************************************
    use size_def, only: ncells,ncellsD,nmaxfaces
    use geo_def, only: cell2cell,ncface
    use comvarbl,only: relaxsor
    use prec_def
    implicit none
    !Input/Output
    integer,intent(in) :: nmaxiter
    real(ikind),intent(in) :: acoef(nmaxfaces,ncellsD),ap(ncells),su(ncellsD)
    real(ikind),intent(inout) :: phi(ncellsD)
    !Internal Variables
    integer :: iter,i
    real(ikind) :: comrelaxsor,term,relaxsor2
    
!inner iteration loop
!$OMP PARALLEL SHARED(ncface,acoef,ap,su,phi,relaxsor2,comrelaxsor)
    do iter=1,nmaxiter/2        
       relaxsor2 = 1.0 + 0.8*tanh(0.3*real(iter-1,kind=ikind))
       comrelaxsor = 1.0 - relaxsor2
!$OMP DO PRIVATE(i,term)
       do i=1,ncells
          term=sum(acoef(1:ncface(i),i)*phi(cell2cell(1:ncface(i),i))) !Must be separate from line below
          phi(i)=relaxsor2*(term+su(i))/ap(i)+comrelaxsor*phi(i)
      enddo
!$OMP END DO
!$OMP DO PRIVATE(i,term)
       do i=ncells,1,-1
          term=sum(acoef(1:ncface(i),i)*phi(cell2cell(1:ncface(i),i))) !Must be separate from line below
          phi(i)=relaxsor2*(term+su(i))/ap(i)+comrelaxsor*phi(i)
      enddo
!$OMP END DO
    enddo
!$OMP END PARALLEL
    return
    end subroutine solv_ssorac

!*****************************************************************************
    subroutine ilu0(n)
!*****************************************************************************    
    use solv_def 
    use size_def
    use prec_def
    implicit none
    integer :: n,ju0,i,ii,js,j,jcol,jf,jm,jrow,jj,jw
    integer :: iw(2*n)
    real(ikind) :: tl

    ju0 = n+2
    jlu(1) = ju0
    
!$OMP PARALLEL DO PRIVATE(i)     
    do i=1,n
        iw(i)=0
    enddo
!$OMP END PARALLEL DO 

    !--main loop
    do ii=1,n
        js = ju0     

        ! generating row number ii of L and U.
        do j=ia(ii),ia(ii+1)-1
            jcol = ja(j)
            if (jcol .eq. ii) then
                alu(ii) = aa_matrix(j)    
                iw(jcol) = ii
                ju(ii)  = ju0             
            else
                alu(ju0) = aa_matrix(j)   
                jlu(ju0) = ja(j)
                iw(jcol) = ju0
                ju0 = ju0+1
            endif
        enddo
        jlu(ii+1) = ju0
        jf = ju0-1
        jm = ju(ii)-1

        !--exit if diagonal element is reached.
        do j=js, jm    
            jrow = jlu(j)
            tl = alu(j)*alu(jrow)
            alu(j) = tl

            !--perform  linear combination
            do jj = ju(jrow), jlu(jrow+1)-1     
                jw = iw(jlu(jj))
                if (jw .ne. 0) alu(jw) = alu(jw) - tl*alu(jj)
            enddo
        enddo
        alu(ii) = 1.0/alu(ii)

        !--reset pointer iw to zero
        iw(ii) = 0
        do i = js, jf
            iw(jlu(i)) = 0
        enddo
    enddo

    return
    end subroutine ilu0
    
!*********************************************************************** 
    subroutine ilutp(n,ierr)
!*********************************************************************** 
    use solv_def 
    use size_def
    use comvarbl
    use prec_def
    use diag_lib, only: diag_print_error
    use diag_def, only: msg, msg2, msg3
    implicit none
    integer n,jw(2*n),iwk,iperm(2*n),ierr
    real(ikind) ::  w(n+1),s,tmp,tnorm,xmax,xmax0,fact,abs,t
    integer k,i,j,jrow,ju0,ii,j1,j2,jpos,len,imax,lenu,lenl,jj,mbloc,icut

    !if(nsolv==4) then
    !    lfil=8
    !else
    !    lfil=35
    !endif
    !droptol=1.0e-3
    !permtol=0.5
    
    iwk=80*no_zero
    mbloc=n
     
!----------------------------------------------------------------------- 
!     initialize ju0 (points to next element to be added to alu,jlu) and pointer array.
!-----------------------------------------------------------------------
    ju0 = n+2
    jlu(1) = ju0

    !--integer double pointer array.
!$OMP PARALLEL DO PRIVATE(j)     
    do j=1, n
      jw(n+j)  = 0
      iperm(j) = j
      iperm(n+j) = j
    enddo
!$OMP END PARALLEL DO     
!-----------------------------------------------------------------------
!     beginning of main loop.
!-----------------------------------------------------------------------
    do 500 ii = 1, n
      j1 = ia(ii)
      j2 = ia(ii+1) - 1
      tnorm = 0.0
      do k=j1,j2
        tnorm = tnorm+abs(aa_matrix(k))
      enddo
      tnorm = tnorm/(j2-j1+1)

      !-- unpack L-part and U-part of row of A in arrays  w  --
      lenu = 1
      lenl = 0
      jw(ii) = ii
      w(ii) = 0.0
      jw(n+ii) = ii

      do j = j1, j2
        k = iperm(n+ja(j))
        t = aa_matrix(j)
        if (k .lt. ii) then
          lenl = lenl+1
          jw(lenl) = k
          w(lenl) = t
          jw(n+k) = lenl
        else if (k .eq. ii) then
          w(ii) = t
        else
          lenu = lenu+1
          jpos = ii+lenu-1 
          jw(jpos) = k
          w(jpos) = t
          jw(n+k) = jpos
        endif
      enddo
      jj = 0
      len = 0 

      !-- eliminate previous rows
150   jj = jj+1
      if (jj .gt. lenl) goto 160
!-----------------------------------------------------------------------
!     in order to do the elimination in the correct order we must select
!     the smallest column index among jw(k), k=jj+1, ..., lenl.
!-----------------------------------------------------------------------
      jrow = jw(jj)
      k = jj

      !-- determine smallest column index
      do j=jj+1,lenl
        if (jw(j) .lt. jrow) then
          jrow = jw(j)
          k = j
        endif
      enddo

      if (k .ne. jj) then
        !-- exchange in jw
        j = jw(jj)
        jw(jj) = jw(k)
        jw(k) = j
        !-- exchange in jr
        jw(n+jrow) = jj
        jw(n+j) = k
        !-- exchange in w
        s = w(jj)
        w(jj) = w(k)
        w(k) = s
      endif

      !-- zero out element in row by resetting jw(n+jrow) to zero.     
      jw(n+jrow) = 0

      !--get the multiplier for row to be eliminated: jrow
      fact = w(jj)*alu(jrow)

      !-- drop term if small     
      if (abs(fact) .le. droptol) goto 150

      !-- combine current row and row jrow
      do k = ju(jrow), jlu(jrow+1)-1
        s = fact*alu(k)
        !-- new column number
        j = iperm(n+jlu(k))
        jpos = jw(n+j)
        if (j .ge. ii) then
          !-- dealing with upper part.
          if (jpos .eq. 0) then
            !-- this is a fill-in element
            lenu = lenu+1
            i = ii+lenu-1 
            jw(i) = j
            jw(n+j) = i 
            w(i) = - s
          else
            !--  no fill-in element --
            w(jpos) = w(jpos) - s
          endif
        else
          !-- dealing with lower part.
          if (jpos .eq. 0) then
            !--  this is a fill-in element
            lenl = lenl+1
            jw(lenl) = j
            jw(n+j) = lenl
            w(lenl) = - s
          else
            !-- this is not a fill-in element
            w(jpos) = w(jpos) - s
          endif
        endif
      enddo
     
      !-- store this pivot element -- (from left to right -- no danger of overlap with the working elements in L (pivots). 
      len = len+1 
      w(len) = fact
      jw(len)  = jrow
      goto 150

160   continue

      !--reset double-pointer to zero (U-part)     
      do k=1, lenu
        jw(n+jw(ii+k-1)) = 0
      enddo

      !--update L-matrix
      lenl = len 
      len = min0(lenl,lfil)
     
      !--sort by quick-split
      call qsplit (w,jw,lenl,len)

      !--     store L-part -- in original coordinates ..
      do k=1, len
        alu(ju0) =  w(k)  
        jlu(ju0) = iperm(jw(k))
        ju0 = ju0+1
      enddo

      !--  save pointer to beginning of row ii of U
      ju(ii) = ju0

      !-- update U-matrix -- first apply dropping strategy
      len = 0
      do k=1, lenu-1
        if (abs(w(ii+k)) .gt. droptol*tnorm) then 
          len = len+1
          w(ii+len) = w(ii+k) 
          jw(ii+len) = jw(ii+k) 
        endif
      enddo
      lenu = len+1
      len = min0(lenu,lfil)
      call qsplit (w(ii+1), jw(ii+1), lenu-1,len)

      !-- determine next pivot -- 
      imax = ii
      xmax = abs(w(imax))
      xmax0 = xmax
      icut = ii - 1 + mbloc - mod(ii-1,mbloc)
      do k=ii+1,ii+len-1
        t = abs(w(k))
        if (t .gt. xmax .and. t*permtol .gt. xmax0 .and.jw(k) .le. icut) then
          imax = k
          xmax = t
        endif
      enddo

      !--exchange w's
      tmp = w(ii)
      w(ii) = w(imax)
      w(imax) = tmp

      !--update iperm and reverse iperm
      j = jw(imax)
      i = iperm(ii)
      iperm(ii) = iperm(j)
      iperm(j) = i

      !--reverse iperm
      iperm(n+iperm(ii)) = ii
      iperm(n+iperm(j)) = j
 
      !--copy U-part in original coordinates     
      do k=ii+1,ii+len-1 
        jlu(ju0) = iperm(jw(k))
        alu(ju0) = w(k)
        ju0 = ju0+1
      enddo

      !-- store inverse of diagonal element of u
      if (w(ii) .eq. 0.0) w(ii) = (1.0e-4 + droptol)*tnorm
      alu(ii) = 1.0d0/ w(ii) 

      !-- update pointer to beginning of next row of U.
      jlu(ii+1) = ju0
      !--end main loop
500 continue

    !-- permute all column indices of LU ...
    do k=jlu(1),jlu(n+1)-1
      jlu(k) = iperm(n+jlu(k))
    enddo

    !--     ...and of A
    do k=ia(1), ia(n+1)-1
      ja(k)=iperm(n+ja(k))
    enddo

    ierr = 0

    return
    end subroutine ilutp
    
!***********************************************************************
     subroutine indexm 
!     index for matrix A
!***********************************************************************
    use size_def
    use geo_def, only: ncface,cell2cell
    use struct_def
    use comvarbl
    use solv_def 
    implicit none
    integer :: i,ntt2,k,kk,ntt1,ntt3,nterm,nterm1,nterm2

    nterm=0
    ia(1)=1

    do i=1,ncells
      ntt2=0
      nterm2=0    
      do k=1,ncface(i)
        if(cell2cell(k,i).lt.i)then
          nterm1=ncells
          do kk=1,ncface(i)
            if(cell2cell(kk,i).lt.i .and. cell2cell(kk,i).lt.nterm1 &
                 .and. cell2cell(kk,i).gt.nterm2)then
              nterm1=cell2cell(kk,i) 
                 ntt1=kk                    
            endif
          enddo
          nterm=nterm+1
          ntt2=ntt2+1
          nterm2=cell2cell(ntt1,i)
          ja(nterm)=cell2cell(ntt1,i)
          nposition1(i,ntt2)=ntt1
        endif
      enddo
      nlow(i)=ntt2
        
!--------
      nterm=nterm+1
      ja(nterm)=i
!--------
      ntt3=0
      nterm2=0
      do k=1,ncface(i)
        if(cell2cell(k,i).gt.i .and. cell2cell(k,i)<=ncells)then
          nterm1=ncells
          do kk=1,ncface(i)
            if(cell2cell(kk,i).gt.i .and. cell2cell(kk,i)<=ncells .and.  &
               cell2cell(kk,i).lt.nterm1 .and. cell2cell(kk,i).gt.nterm2)then
              nterm1=cell2cell(kk,i) 
              ntt1=kk 
            endif
          enddo
          nterm=nterm+1
          ntt3=ntt3+1
          nterm2=cell2cell(ntt1,i)
          ja(nterm)=cell2cell(ntt1,i)
          nposition2(i,ntt3)=ntt1
        endif
      enddo
      nup(i)=ntt3
      ia(i+1)=nterm+1
    enddo
    no_zero=nterm

    call diagonal_pointer_cr !iua
    call matrix_pointer_cr !iaa

    return
    end subroutine indexm

!**********************************************************************
    subroutine diagonal_pointer_cr 
!     DIAGONAL_POINTER_CR finds diagonal entries in a sparse compressed row matrix.
!*********************************************************************** 
    use size_def
    use solv_def
    implicit none
    integer :: i,k

    iua(1:ncells)=-1
    do i=1,ncells
      do k=ia(i),ia(i+1)-1
        if (ja(k)== i) then
            iua(i)=k
        end if
      enddo
    enddo
    
    return
    end subroutine

!*********************************************************************** 
    subroutine matrix_pointer_cr
!*********************************************************************** 
    use size_def, only: ncells
    use solv_def
    use prec_def
    implicit none
    !Internal Variables
    integer :: i,k,nterm1

    !----Form A matrix using a sparse compressed row      
    nterm1=0 
    do i=1,ncells
      do k=1,nlow(i)
        nterm1=nterm1+1  
        iaa(i,nposition1(i,k))=nterm1
        !!aa_matrix(nterm1)=-acoef(nposition1(i,k),i)
      enddo
      nterm1=nterm1+1
      iaa(i,0)=nterm1
      !!aa_matrix(nterm1)=ap(i)
      do k=1,nup(i)
        nterm1=nterm1+1  
        iaa(i,nposition2(i,k))=nterm1
        !!aa_matrix(nterm1)=-acoef(nposition2(i,k),i)
      enddo
    enddo

    return
    end subroutine    

!*********************************************************************** 
    subroutine coef2csr(acoef,ap)
! Form A matrix using a sparse compressed row      
!*********************************************************************** 
    use size_def
    use solv_def
    use prec_def
    implicit none
    !Input/Output
    real(ikind),intent(in) :: acoef(nmaxfaces,ncellsD),ap(ncells)
    !Internal Variables
    integer :: i,k,kk

    !----Form A matrix using a sparse compressed row      
!$OMP PARALLEL DO PRIVATE(i,k,kk)  
    do i=1,ncells
      do k=1,nlow(i)
        kk=nposition1(i,k)
        aa_matrix(iaa(i,kk))=-acoef(kk,i)
      enddo
      aa_matrix(iaa(i,0))=ap(i)
      do k=1,nup(i)
        kk=nposition2(i,k)
        aa_matrix(iaa(i,kk))=-acoef(kk,i)
      enddo
    enddo
!$OMP END PARALLEL DO

    return
    end subroutine    
    
!*********************************************************************** 
    subroutine ilu_cr(acoef,ap)
!     ILU_CR computes the incomplete LU factorization of a matrix.
!*********************************************************************** 
    use size_def
    use solv_def
    use prec_def
    implicit none
    !Input/Output
    real(ikind),intent(in) :: acoef(nmaxfaces,ncellsD),ap(ncells)
    !Internal Variables
    integer :: i,k,nterm1

    !----Form A matrix using a sparse compressed row      
    nterm1=0 
    do i=1,ncells
      do k=1,nlow(i)
        nterm1=nterm1+1  
        aa_matrix(nterm1)=-acoef(nposition1(i,k),i)
      enddo
      nterm1=nterm1+1
      aa_matrix(nterm1)=ap(i)
      do k=1,nup(i)
        nterm1=nterm1+1  
        aa_matrix(nterm1)=-acoef(nposition2(i,k),i)
      enddo
    enddo

    return
    end subroutine
      
!*****************************************************************************
    subroutine qsplit(a,ind,n,ncut)
!*****************************************************************************
!     does a quick-sort split of a real array.
!-----------------------------------------------------------------------
    use prec_def
    implicit none
    integer :: n,mid,j
    integer :: ind(n),ncut,itmp,first,last
    real(ikind) :: a(n),tmp, abskey

    first = 1
    last = n
    if (ncut .lt. first .or. ncut .gt. last) return

1      mid = first
    abskey = abs(a(mid))
    do 2 j=first+1, last
        if (abs(a(j)) .gt. abskey) then
            mid = mid+1
            !--interchange
            tmp = a(mid)
            itmp = ind(mid)
            a(mid) = a(j)
            ind(mid) = ind(j)
            a(j)  = tmp
            ind(j) = itmp
        endif
2      continue

    !--interchange
    tmp = a(mid)
    a(mid) = a(first)
    a(first)  = tmp

    itmp = ind(mid)
    ind(mid) = ind(first)
    ind(first) = itmp

    !--test for while loop
    if (mid .eq. ncut) return
    if (mid .gt. ncut) then
        last = mid-1
    else
        first = mid+1
    endif
    goto 1

    end subroutine qsplit

!*****************************************************************************
    subroutine lusol(n, y, x)
! This routine solves the system (LU) x = y,        
!*****************************************************************************
    use solv_def
    use size_def       
    use prec_def
    implicit none
    integer n,i,k         
    real(ikind) x(n), y(n)
    
    !-forward solve
    do i=1,n
      x(i) = y(i)
      do k=jlu(i),ju(i)-1                        
        x(i)=x(i)-alu(k)*x(jlu(k))            
      enddo
    enddo

    !-backward solve.
    do i = n, 1, -1
      do k=ju(i),jlu(i+1)-1                      
        x(i)=x(i)-alu(k)*x(jlu(k))
      enddo
      x(i) = alu(i)*x(i)
    enddo
    
    return
    end subroutine lusol    

!*****************************************************************************
    subroutine ceaambk(n,aa,ja,ia,b,c)
! Multiplication of sparse matrix in CSR format times a vector (c = aa*b)
! with loop unrolling, arbitrary precision, and OpemMP parallelization.
!
!Variables and Notations:
!  n: order of the matrix
!  aa(*): input sparse matrix in CSR format
!  ia(n+1): input index for s
!  ja(*): = input index for s
!  b(*): = input vector to be multiplied
!  c(*): = output vector
!  e: symbol for equal sign
!  m: symbol for multiplication
!  8: indicates double precision
!
! written by Alex Sanchez, USACE-CHL
!*****************************************************************************
    use prec_def
    implicit none
    integer:: i,k1,k2    
    integer,intent(in):: n,ia(n+1),ja(*)
    real(ikind),intent(in) :: aa(*),b(*)
    real(ikind),intent(out):: c(*)

!$OMP PARALLEL DO PRIVATE(i,k1,k2)
    do i=1,n 
      k1=ia(i)
      k2=ia(i+1)-1
      c(i)=sum(aa(k1:k2)*b(ja(k1:k2)))
      !c(i)=dot_product(aa(k1:k2),b(ja(k1:k2))
    enddo
!$OMP END PARALLEL DO 

    return
    end subroutine ceaambk
    
