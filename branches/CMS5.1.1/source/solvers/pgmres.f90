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
!$OMP PARALLEL DO PRIVATE(ii,term) REDUCTION(+:ro)   
      do ii=1,ncells
        term=sum(acoef(1:ncface(ii),ii)*phi(cell2cell(1:ncface(ii),ii)))
        vv(ii,1)=term+su(ii)-ap(ii)*phi(ii)
        ro=ro+vv(ii,1)*vv(ii,1)
      enddo 
!$OMP END PARALLEL DO 
      ro=sqrt(ro+small)     
       
      if(itr.eq.1) eps11=eps1*ro
      t=1.0/(ro+small)
      call aesmak(ncells,t,vv(1:ncells,1))
      
      !---Initialize first term of RHS of Hessenberg system
      rs1(1)=ro
!	   rs1(2:kmax+1)=0.0
!	   hh(2:kmax+1,1:kmax+1)=0.0

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
    endsubroutine pgmres
