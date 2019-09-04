!==============================================    
module LU_lib
! Lower-Upper decomposition by Crout's algorithm
! with partial pivoting. 
!==============================================
    implicit none
    
contains

!***************************************************
    subroutine LU_decomp(n,A,ind,ierr)
!***************************************************    
    use prec_def
    implicit none
    !Input/Output
    integer,    intent(in) :: n
    real(ikind),intent(inout) :: A(n,n)
    integer,    intent(out) :: ind(n)
    integer,    intent(out) :: ierr
    !Internal Variables
    integer :: i,j,k,imax
    real(ikind), parameter :: small=1.0e-8
    real(ikind) :: vv(n),maxval,val,sumval

    ierr=0
    do i=1,n
      maxval=0.0_ikind
      do j=1,n
        if(abs(A(i,j))>maxval) maxval=abs(A(i,j))
      enddo !j
      if(maxval<small) then
        ierr = -1
        return
      endif
      vv(i) = 1.0_ikind/maxval !Scaling
    enddo !i

    do j=1,n
      do i=1,j-1
        sumval = A(i,j)
        do k=1,i-1
          sumval = sumval - A(i,k)*A(k,j) 
        enddo !k
        A(i,j) = sumval
      enddo !i
      maxval = 0.0_ikind
      do i=j,n
        sumval = A(i,j)
        do k=1,j-1
          sumval = sumval - A(i,k)*A(k,j) 
        enddo !k
        A(i,j) = sumval
        val = vv(i)*abs(sumval)
        if(val>=maxval)then
          imax = i
          maxval = val
        endif
      enddo !i   
      if(j/=imax)then
        do k=1,n
          val = A(imax,k)
          A(imax,k) = A(j,k)
          A(j,k) = val
        enddo !k
        vv(imax) = vv(j)
      endif
      ind(j) = imax
      if(abs(A(j,j))<small) A(j,j) = small
      if(j/=n)then
        val = 1.0_ikind/A(j,j)
        do i=j+1,n
          A(i,j) = A(i,j)*val
        enddo !i
      endif 
    enddo !j

    return
    endsubroutine LU_decomp
    
!***************************************************    
     subroutine LU_subs(n,A,ind,b)
!***************************************************    
    use prec_def
    implicit none
    !Input/Output
    integer,    intent(in) :: n
    real(ikind),intent(in) :: A(n,n)
    integer,    intent(in) :: ind(n)
    real(ikind),intent(inout) :: b(n)
    !Internal Variables
    integer :: i,ii,j,jj
    real(ikind) :: sumval
    
    ii = 0    
    do i=1,n
      jj = ind(i)
      sumval = b(jj)
      b(jj) = b(i)
      if(ii/=0)then
        do j=ii,i-1
          sumval = sumval - A(i,j)*b(j)
        enddo !j
      elseif(sumval/=0.0)then
        ii = i
      endif
      b(i) = sumval !Intermediate vector
    enddo !i

    do i=n,1,-1
      sumval = b(i)
      if(i<n)then
        do j=i+1,n
          sumval = sumval - A(i,j)*b(j)
        enddo !j
      endif
      b(i) = sumval/A(i,i)
    enddo !i

    return
    endsubroutine LU_subs
     
endmodule LU_lib     
    