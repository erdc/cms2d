!========================================================================
! CMS Hiding and Exposure Routines
!
! Contains the following
!   HidExpEgiazaroff - Calculates hiding and exposure coefficient
!                      based on Egiazaroff (1965)
!   HidExpAshidaMichiue - Calculates hiding and exposure coefficient 
!                         based on shida and Michiue (1971)
!   HidExpHayashi  - Calculates hiding and exposure coefficient
!                    based on Ashida and Michiue (1971)
!   HidExpParker - Calculates hiding and exposure coefficient based on 
!                  Parker et al. (1982) and others
!   HidExpWu - Calculates hiding and exposure coefficient 
!              based on Wu et al. (2000) 
!
! written by Alex Sanchez, USACE-CHL
!========================================================================
    
!************************************************************************
    subroutine HidExpEgiazaroff
! Calculates hiding and exposure coefficient based on Egiazaroff (1965)
!
! Reference:
!   Egiazaroff, I.V. 1965. Calculation of nonuniform sediment concentration, 
!     Journal of Hydraulics Division, ASCE, 91(H14), 225-247. 
!
! written by Alex Sanchez, USACE-CHL
!************************************************************************
    use size_def, only: ncells
    use sed_def, only: nsed,diam,pbk,varsigma,varsigmamin,varsigmamax
    use sed_lib, only: hidexp_egiazaroff
    use prec_def
    implicit none
    integer :: i,k
    real(ikind) :: dm
    
!$OMP PARALLEL DO PRIVATE(i,k,dm)
    do i=1,ncells
      dm = sum(pbk(i,:,1)*diam(:)) !Mean diameter
      do k=1,nsed          
        varsigma(i,k) = hidexp_egiazaroff(diam(k),dm) !Egiazaroff (1965)
        varsigma(i,k) = max(min(varsigma(i,k),varsigmamax),varsigmamin) !Limit values
      enddo !k
    enddo !i  
!$OMP END PARALLEL DO

    return
    endsubroutine hidexpEgiazaroff
    
!************************************************************************
    subroutine HidExpAshidaMichiue
! Calculates hiding and exposure coefficient based on 
! Ashida and Michiue (1971)
!
! Reference:
!   Ashida, K., and Michiue, M. 1971. An investigation of river bed 
!     degradation downstream of a dam. Proceedings of the 14th Congress 
!     of IAHR, Paris, France, 3, 247-256. 
!
! written by Alex Sanchez, USACE-CHL
!************************************************************************
    use size_def, only: ncells
    use sed_def, only: nsed,diam,pbk,varsigma,varsigmamin,varsigmamax
    use sed_lib, only: hidexp_ashida_michiue
    use prec_def
    implicit none
    integer :: i,k
    real(ikind) :: dm
    
!$OMP PARALLEL DO PRIVATE(i,k,dm)
    do i=1,ncells
      dm = sum(pbk(i,:,1)*diam(:)) !Mean diameter
      do k=1,nsed          
        varsigma(i,k) = hidexp_ashida_michiue(diam(k),dm)
        varsigma(i,k) = max(min(varsigma(i,k),varsigmamax),varsigmamin) !Limit values
      enddo !k
    enddo !i  
!$OMP END PARALLEL DO

    return
    endsubroutine HidExpAshidaMichiue    

!************************************************************************
    subroutine HidExpHayashi
! Calculates hiding and exposure coefficient based on 
! Hayashi et al. (1980)
!
! Reference:
!   Hayashi, T.S., Ozaki, and Ichibashi, T. 1981. Study on bed load 
!     transport of sediment mixture, Proceedings of the 24th Japanese 
!     Conference on Hydraulics, Japan.
! 
! written by Alex Sanchez, USACE-CHL
!************************************************************************
    use size_def, only: ncells
    use sed_def, only: nsed,diam,pbk,varsigma,varsigmamin,varsigmamax
    use sed_lib, only: hidexp_hayashi
    use prec_def
    implicit none
    integer :: i,k
    real(ikind) :: dm
    
!$OMP PARALLEL DO PRIVATE(i,k,dm)
    do i=1,ncells
      dm = sum(pbk(i,:,1)*diam(:)) !Mean diameter
      do k=1,nsed          
        varsigma(i,k) = hidexp_hayashi(diam(k),dm)
        varsigma(i,k) = max(min(varsigma(i,k),varsigmamax),varsigmamin) !Limit values
      enddo !k
    enddo !i  
!$OMP END PARALLEL DO

    return
    endsubroutine HidExpHayashi      
    
!************************************************************************
    subroutine HidExpParker
! Calculates hiding and exposure coefficient based on 
! Parker et al. (1982) and others
! written by Alex Sanchez, USACE-CHL
!************************************************************************
    use size_def, only: ncells
    use sed_def, only: nsed,diam,d50,mhe,varsigma,varsigmamin,varsigmamax
    use sed_lib, only: hidexp_parker
    use prec_def
    implicit none
    integer :: i,k
    
!$OMP PARALLEL DO PRIVATE(i,k) COLLAPSE(2)
    do i=1,ncells
      do k=1,nsed
        varsigma(i,k) = hidexp_parker(diam(k),d50(i),mhe)
        varsigma(i,k) = max(min(varsigma(i,k),varsigmamax),varsigmamin) !Limit values
      enddo !k 
    enddo !i
!$OMP END PARALLEL DO

    return     
    endsubroutine hidexpparker    
    
!************************************************************************
    subroutine HidExpWu
! Calculates hiding and exposure coefficient based on Wu et al. (2000)
! written by Alex Sanchez, USACE-CHL
!************************************************************************
    use size_def, only: ncells
    use sed_def, only: nsed,diam,pbk,mhe,varsigma,varsigmamin,varsigmamax
    use sed_lib, only: hidexp_wu
    use prec_def
    implicit none
    integer :: i,k
    real(ikind) :: pbm(nsed)

!$OMP PARALLEL DO PRIVATE(i,k)    
    do i=1,ncells
      pbm = pbk(i,:,1)  
      do k=1,nsed
        varsigma(i,k) = hidexp_wu(nsed,diam,k,pbm,mhe)
        varsigma(i,k) = max(min(varsigma(i,k),varsigmamax),varsigmamin) !Limit values
      enddo !k
    enddo !i
!$OMP END PARALLEL DO

    return
    endsubroutine HidExpWu
