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
    endsubroutine sip5
