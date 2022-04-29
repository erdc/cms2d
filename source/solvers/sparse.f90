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
    