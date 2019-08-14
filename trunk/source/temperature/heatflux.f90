!================================================================= 
    subroutine heatfluxsrc
!   With shortwave radiation data
!================================================================= 
       use size_def, only: ncells
       use geo_def, only: areap
       use met_def, only: uwind,vwind,wndx,wndy,windvar
       use flow_def, only: iwet,su,sp,rhow
       use heat_def
       use comvarbl, only: ctime
       implicit none
       integer :: i
       real :: clc,emair,halwdwn,halw,vpair,evapwind,hsensib,vpw,ctermh,hevap,hasw
       
       hasw=solarms   !Using measured shortwave radiation data
      
       !Long wave atmospheric radiation.
       clc=1.0+0.17*cloud**2
       emair=0.938*(airtmp+273.16)**2/100000.0
       halwdwn=0.97*emair*((airtmp+273.16)/100.0)**4*clc*5.669
       halw=halwdwn-0.96*2.7316**4*5.669   !in W/m2

       vpair=6.108*exp(17.27*dewpt/(dewpt+273.3))    !VPair in mb and T in oC      

       do i=1,ncells
 !        write(40,*)ctime,wndx,wndy
         if(windvar) then
            velwind=sqrt(uwind(i)**2+vwind(i)**2)
         else
            velwind=sqrt(wndx*wndx+wndy*wndy)
 !           write(40,*)ctime,wndx,wndy,velwind
         endif
 !         velwind=10.0
 !         velwind=2.0
          evapwind=6.9+0.345*velwind**2       ! Edinger et al. 1974   !wind speed at 7 m above 

          hsensib=0.62*evapwind*airtmp        ! in W/m2   Linearization

          vpw=6.108*exp(17.27*heat(i)/(heat(i)+273.3)) !VPW in mb and T in oC
          hevap=evapwind*(vpair-vpw)          ! in W/m2

          !Compute net heat flux from all sources.
          ctermh=1.0/4186.0/rhow*iwet(i)*areap(i)      !Cp=4186.J/Kg.oC          
          sp(i)=sp(i)-0.62*evapwind*ctermh                  &
                     -0.96*5.669*0.00000001*ctermh*         &
                      (heat(i)**3+4.0*heat(i)**2*273.16     &
                     +6.0*heat(i)*273.16**2+4.0*273.16**3 )
                   ! -0.96*5.669*0.01*ctermh*    &
                   !     ((0.01*heat(i))**3+4.0*(0.01*heat(i))**2*2.7316   &
                   !      +6.0*0.01*heat(i)*2.7316**2+4.0*2.7316**3 )
          suheat0(i)=suheat0(i)+(hasw+halw+hevap+hsensib)*ctermh
       enddo

       return
    end subroutine heatfluxsrc

!================================================================= 
    subroutine read_heat_solar_xmdf
!================================================================= 
#include "CMS_cpp.h"    
    use heat_def
    use comvarbl,    only: mpfile
#ifdef XMDF_IO    
    use in_xmdf_lib, only: read_dataseth5
#else
    use diag_lib, only: diag_print_error
#endif      
    implicit none
    integer :: ierr,i
    real(ikind),pointer:: vtemp(:)
      
#ifdef XMDF_IO    
    do i=1,5
      select case (i)
      case (1)
        call read_dataseth5(mpfile,heatpath,"AirTemp",numhtflux,vtemp)
        allocate(clouddat(numhtflux),dewptdat(numhtflux),airtmpdat(numhtflux),dathtflux(numhtflux),solarmsdat(numhtflux))
        airtmpdat=vtemp
      case (2)
        call read_dataseth5(mpfile,heatpath,"CloudCover",numhtflux,vtemp)
        clouddat=vtemp
      case (3)
        call read_dataseth5(mpfile,heatpath,"DewPoint",numhtflux,vtemp)
        dewptdat=vtemp
      case(4)
        call read_dataseth5(mpfile,heatpath,"SolarRadiation",numhtflux,vtemp)
        solarmsdat=vtemp
      case(5)
        call read_dataseth5(mpfile,heatpath,"Times",numhtflux,vtemp)
        dathtflux=vtemp
      end select
    enddo
#else
    call diag_print_error('Temperature not available without XMDF at present time')
#endif
    
    deallocate(vtemp)

    return
    end subroutine read_heat_solar_xmdf

!================================================================= 
    subroutine read_heat_solar       
!================================================================= 
       use heat_def
       implicit none
       integer :: iht
       
       !open(unit=41,file=file41)

       read(41,*) numhtflux
       write(*,*) 'Temperature num=',numhtflux

       allocate ( clouddat(numhtflux),dewptdat(numhtflux) )
       allocate ( airtmpdat(numhtflux),dathtflux(numhtflux) )
       allocate ( solarmsdat(numhtflux) )

       do iht=1,numhtflux
          read(41,*)  dathtflux(iht),airtmpdat(iht),dewptdat(iht),   &
                        clouddat(iht),solarmsdat(iht)
       enddo     
  
       close(41)

       return
    end subroutine read_heat_solar
    
!================================================================= 
    subroutine heatdatainterpol
!================================================================= 
       use comvarbl, only: ctime
       use heat_def
       implicit none
       integer :: iht
       real :: theat
       
       theat=ctime/3600.0

       if(theat.le.dathtflux(1)) then
            airtmp =airtmpdat(1)
            dewpt =dewptdat(1)
            cloud =clouddat(1)
            solarms=solarmsdat(1)
       endif
       if(theat.ge.dathtflux(numhtflux)) then
            airtmp =airtmpdat(numhtflux)
            dewpt =dewptdat(numhtflux)
            cloud =clouddat(numhtflux)
            solarms=solarmsdat(numhtflux)
       endif
       do iht=1,numhtflux-1
          if(theat.gt.dathtflux(iht).and.   &
                        theat.le.dathtflux(iht+1)) then 
             airtmp=airtmpdat(iht)+(airtmpdat(iht+1)-airtmpdat(iht))   &
                                   *(theat-dathtflux(iht))              &
                                   /(dathtflux(iht+1)-dathtflux(iht))
             dewpt=dewptdat(iht)+(dewptdat(iht+1)-dewptdat(iht))       &
                                    *(theat-dathtflux(iht))             &
                                    /(dathtflux(iht+1)-dathtflux(iht))
             cloud=clouddat(iht)+(clouddat(iht+1)-clouddat(iht))       &
                                   *(theat-dathtflux(iht))              &
                                   /(dathtflux(iht+1)-dathtflux(iht))
             solarms=solarmsdat(iht)+(solarmsdat(iht+1)-solarmsdat(iht))  &
                                  *(theat-dathtflux(iht))                 &
                                  /(dathtflux(iht+1)-dathtflux(iht))
          endif
       enddo       
        !write(*,*) ctime,airtmp,dewpt, cloud, solarms         
    return
    end subroutine heatdatainterpol
    