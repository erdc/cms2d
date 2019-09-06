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
    endsubroutine iccg5    
