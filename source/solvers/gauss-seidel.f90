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
