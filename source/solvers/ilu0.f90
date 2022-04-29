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