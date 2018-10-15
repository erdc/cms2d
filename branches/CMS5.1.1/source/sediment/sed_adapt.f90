!========================================================================
! CMS Sediment Adaptation Coefficient
!
! Contains the following
!     adaptcoef   - Computes the total-load adpatation coefficient
!     adaptlenbed - Computes the bed-load adpatation length
!     adaptlensus - Computes the suspended-load adapatation length
!========================================================================
    
!************************************************************************
    subroutine adaptcoef
! Calculates adaptation length time sediment fall velocity    
! written by Alex Sanchez, USACE-CHL;  Weiming Wu, NCCHE      
!************************************************************************      
    use size_def
    use flow_def
    use comvarbl
    use sed_def
    use sed_lib, only: fallvel_soulsby,fallvel_wu_wang
    use prec_def
    implicit none
    integer :: i
    real(ikind) :: wf(ncellsD),dstar50,shape50
    
!    kt=max(nsed/2,1) !Transport grain size ********
    
    if(singlesize)then
      wf=wsfall(1)
    else
      if(iws==1)then !Soulsby   
        do i=1,ncells  
          dstar50 = d50(i)*d2dstar
          wf(i) = fallvel_soulsby(viscos,d50(i),dstar50,0.0_ikind)
        enddo
      else !Wu & Wang
        do i=1,ncells    
          dstar50 = d50(i)*d2dstar
          shape50 = sum(pbk(i,1:nsed,1)*coreyshape(1:nsed)) !Approximate
          wf(i) = fallvel_wu_wang(viscos,shape50,d50(i),dstar50) !shape could be interpolated from size class shape factors
        enddo
      endif
    endif
    
    selectcase(iadapttot)
    case(1) !Constant length
      !$omp parallel do private(i)
      do i=1,ncells
        alphat(i) = uv(i)*h(i)/(Ltot*wf(i))*iwet(i)
      enddo
      !$omp end parallel do 
      
    case(2) !Constant time 
      !$omp parallel do private(i)
      do i=1,ncells          
        alphat(i) = h(i)/(Ttot*wf(i))*iwet(i)
      enddo
      !$omp end parallel do 
      
    case(3) !Spatially variable
      !$omp parallel do private(i)
      do i=1,ncells
        alphat(i) = uv(i)*h(i)/(vLtot(i)*wf(i))*iwet(i)
      enddo
      !$omp end parallel do 
      
    case(4) !Erosion and deposition lengths
      !$omp parallel do private(i)
      do i=1,ncells
        if(dzb(i)>0.0)then !Deposition
          alphat(i) = uv(i)*h(i)/(Ldeptot*wf(i))*iwet(i)
        else  !Erosion
          alphat(i) = uv(i)*h(i)/(Lerotot*wf(i))*iwet(i)
        endif
      enddo
      !$omp end parallel do
               
    case(5) !Max of Bed and Suspended Lengths      
      call adaptlenbed     !Bed adaptation length   
      call adaptlensus(wf) !Suspended adaptation length 
      !$omp parallel do private(i)
      do i=1,ncells
        vLtot(i) = max(vLbed(i),vLsus(i))
        alphat(i) = uv(i)*h(i)/(vLtot(i)*wf(i))*iwet(i)
      enddo
      !$omp end parallel do
      
    case(6) !Weighted average of Bed and Suspended Lengths                      
      call adaptlenbed     !Bed adaptation length   
      call adaptlensus(wf) !Suspended adaptation length  
      !$omp parallel do private(i)
      do i=1,ncells
        vLtot(i) = (1.0-rs(i))*vLbed(i) + rs(i)*vLsus(i)
        alphat(i) = uv(i)*h(i)/(vLtot(i)*wf(i))*iwet(i)
      enddo
      !$omp end parallel do
      
    endselect
    
    return
    endsubroutine adaptcoef  

!***************************************************        
    subroutine adaptlenbed
! Calculates the bed load adaptation length
! written by Alex Sanchez, USACE-CHL
!***************************************************   
    use size_def
    use flow_def, only: uv,h,hdry
    use fric_def, only: bsvel
    use sed_def, only: iadaptbed,vlbed,Tbed,Lbed,Ldepbed,Lerobed,dzb,fbed
    use prec_def
    implicit none
    integer :: i
    real(ikind) :: Lbedmin
    
    Lbedmin = 5.0*hdry
    selectcase(iadaptbed)
    case(2) !Adapatation time
      !$omp parallel do private(i)
      do i=1,ncells
        vLbed(i) = h(i)/Tbed
        vLbed(i) = max(vLbed(i),Lbedmin)
      enddo  
      !$omp end parallel do 
      
    case(4) !Erosion and Deposition lengths
      !$omp parallel do private(i)    
      do i=1,ncells
        if(dzb(i)>0.0)then !Deposition
          vLbed(i) = Ldepbed
        else !Erosion
          vLbed(i) = Lerobed
        endif
        vLbed(i) = max(vLbed(i),Lbedmin)
      enddo  
      !$omp end parallel do 
        
    case(5) !Depth-dependant
      !$omp parallel do private(i)    
      do i=1,ncells  
        vLbed(i) = h(i)*fbed
        vLbed(i) = max(vLbed(i),Lbedmin)
      enddo  
      !$omp end parallel do   
        
    endselect
    
    return
    endsubroutine adaptlenbed
    
!***************************************************        
    subroutine adaptlensus(wf)
! Calculates the suspended load adpatation length
! written by Alex Sanchez, USACE-CHL
!***************************************************    
    use size_def
    use flow_def, only: uv,h,hdry
    use fric_def, only: coefman,bsvel,z0
    use sed_def, only: iadaptsus,vLsus,Tsus,Lsus,alphasus,&
        Ldepsus,Lerosus,dzb,d90,nsed
    use sed_lib, only: adaptsusp_armanini_disilvio,adaptsusp_lin,adaptsusp_gallappatti
    use flow_def, only: sqrtgrav
    use prec_def
    implicit none
    !Input/Output
    real(ikind),intent(in) :: wf(ncellsD)
    !Internal Variables
    integer :: i,kt
    real(ikind) :: alphas,Ts,vLsusmin
    
    kt = max(nsed/2,1)
    vLsusmin = 5.0*hdry
    selectcase(iadaptsus)              
      case(2) !Constant Time
        !$omp parallel do private(i)  
        do i=1,ncells  
          vLsus(i) = h(i)/Tsus
          vLsus(i) = max(vLsus(i),vLsusmin)
        enddo  
        !$omp end parallel do
        
      case(3) !Constant Alpha
        !$omp parallel do private(i)  
        do i=1,ncells  
          vLsus(i) = uv(i)*h(i)/(alphasus*wf(i))
          vLsus(i) = max(vLsus(i),vLsusmin)
        enddo  
        !$omp end parallel do
        
      case(4) !Erosion/deposition
        !$omp parallel do private(i)  
        do i=1,ncells
          if(dzb(i)>0.0)then 
            vLsus(i) = Ldepsus !Deposition
          else
            vLsus(i) = Lerosus !Erosion
          endif
          vLsus(i) = max(vLsus(i),vLsusmin)
        enddo  
        !$omp end parallel do
          
      case(5) !Armanini and di Silvio (1986)  
        !$omp parallel do private(i)  
        do i=1,ncells
          alphas = adaptsusp_armanini_disilvio(sqrtgrav,h(i),d90(i),wf(i),bsvel(i))
!!          alphas = armanini_disilvio_test(h(i),z0(i),wf(i),bsvel(i))
          vLsus(i) = uv(i)*h(i)/(alphas*wf(i))
          vLsus(i) = max(vLsus(i),vLsusmin)
        enddo  
        !$omp end parallel do
        
      case(6) !Lin (1984)
        !$omp parallel do private(i)    
        do i=1,ncells
          alphas = adaptsusp_lin(wf(i),bsvel(i)) 
          vLsus(i) = uv(i)*h(i)/(alphas*wf(i))
          vLsus(i) = max(vLsus(i),vLsusmin)
        enddo      
        !$omp end parallel do
        
      case(7) !Gallappatti (1983)
        !$omp parallel do private(i)    
        do i=1,ncells 
          Ts = adaptsusp_gallappatti(h(i),uv(i),wf(i),bsvel(i))
          vLsus(i) = uv(i)*Ts
          vLsus(i) = max(vLsus(i),vLsusmin)
        enddo          
        !$omp end parallel do
        
      case(8) !Xbeach
        !$omp parallel do private(i)    
        do i=1,ncells 
          Ts = max(0.05*h(i)/wf(i),0.2)
          vLsus(i) = uv(i)*Ts
          vLsus(i) = max(vLsus(i),vLsusmin)
        enddo   
        !$omp end parallel do 
        
    endselect
    
    return
    endsubroutine adaptlensus
