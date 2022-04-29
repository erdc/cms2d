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