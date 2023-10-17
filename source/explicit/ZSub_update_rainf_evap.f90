!***********************************************************************
      subroutine rain_evap()
!***********************************************************************
      use EXP_Global_def,    only: etan,rain,evap,RF_frac_ro,rain_time,adeq,rf_unit
      use EXP_transport_def, only: adss,salt 
      use flow_def, only: iwet
      use size_def, only: ncells,ncellsD 
      use geo_def,  only: dx,dy      
      use sal_def,  only: saltrans, sal
       
      implicit none 
      real amount,Wamount,Damount,smass 
      integer i  
      
      !RAINFALL  
      if(rain_time.gt.3600) then
        rain_time = 0
        read(RF_unit,*)rain,evap
        Wamount = rain-evap
        Damount = rain*RF_frac_ro  
            
!$omp parallel do  private(amount)       
        do i=1,ncells
          amount=Wamount
          if(iwet(i) .eq. 0)amount=Damount
          etan(i) = etan(i) + amount
        enddo 
!$omp end parallel do
           
        if(saltrans) then
!$omp parallel do private(amount,smass)           
          do i=1,ncells
            amount=Wamount
            if(iwet(i) .eq. 0)amount=Damount
            smass = salt(i)%vol*sal(i)
            salt(i)%vol = salt(i)%vol + amount*dx(i)*dy(i)
            sal(i) = smass/salt(i)%vol
          enddo  
!$omp end parallel do                  
        endif
          
        if(adeq) then
!$omp parallel do private(amount,smass)          
          do i=1,ncells
            amount=Wamount
            if(iwet(i) .eq. 0)amount=Damount
            smass = adss(i)%vol*adss(i)%conc
            adss(i)%vol = adss(i)%vol + amount*dx(i)*dy(i)
            adss(i)%conc = smass/adss(i)%vol
          enddo 
!$omp end parallel do                   
        endif          
      endif

      return
      end subroutine
