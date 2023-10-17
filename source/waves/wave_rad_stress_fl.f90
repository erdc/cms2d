!*************************************************************************
    subroutine wave_rad_stress_fl
! Calculates the wave radiation stresses on the flow grid    
!*************************************************************************
    use size_def, only: ncells,ncellsD
    use wave_flowgrid_def, only: wavstrx,wavstry,whgt,wper,wlen,wunitx,wunity
    use flow_def, only: u,v,h,rhow,grav
    use const_def, only: twopi
    use der_lib, only: dx2d,dy2d
    use der_def, only: gow
    use prec_def
    implicit none
    integer :: i
    real(ikind) :: Sxx(ncellsD),Sxy(ncellsD),Syy(ncellsD)
    real(ikind) :: c,cn,hwkw,E,dvarx,dvary
    logical :: isnankind
    
    do i=1,ncells
      c=wlen(i)/wper(i)
      hwkw=twopi/wlen(i)*h(i)
      hwkw=min(hwkw,10.0)
      cn=0.5+hwkw/sinh(2.0*hwkw)
      E=0.0625*grav*whgt(i)*whgt(i) !Energy divided by density
      !Radiation stresses        
      Sxx(i)=E*(cn*(1.0+wunitx(i)*wunitx(i))-0.5)
      Syy(i)=E*(cn*(1.0+wunity(i)*wunity(i))-0.5)
      Sxy(i)=E*cn*wunity(i)*wunity(i)
    enddo
    
    do i=1,ncells
      call dx2d(gow,i,Sxx,dvarx)
      call dy2d(gow,i,Sxy,dvary)
      wavstrx(i) = - dvarx - dvary
      call dx2d(gow,i,Sxy,dvarx)         
      call dy2d(gow,i,Syy,dvary)
      wavstry(i) = - dvarx - dvary
      if(isnankind(wavstrx(i)) .or. isnankind(wavstry(i)) .or. &
          abs(wavstrx(i))>10.0 .or. abs(wavstry(i))>10.0)then
        write(*,*) 'ERROR: Found NaN'
        write(*,*) 'i:          ',i
        write(*,*) 'h(i):       ',h(i) 
        write(*,*) 'wavstrx(i): ',wavstrx(i)
        write(*,*) 'wavstry(i): ',wavstry(i)
        write(*,*) 'whgt(i):    ',whgt(i)
        write(*,*) 'wunitx(i):  ',wunitx(i)
        write(*,*) 'wunity(i):  ',wunity(i)
        write(*,*) 'wlen(i):    ',wlen(i)
        write(*,*) 'wper(i):    ',wper(i)
        write(*,*) 'Sxx(i):     ',Sxx(i)
        write(*,*) 'Sxy(i):     ',Sxy(i)
        write(*,*) 'Syy(i):     ',Syy(i)
        write(*,*) '  Press <enter> key to continue'
        read(*,*)
      endif
    enddo
    
    return
    end subroutine wave_rad_stress_fl