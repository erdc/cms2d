!*************************************************************************
    subroutine wave_rad_stress_fl
! Calculates the wave radiation stresses on the flow grid    
!*************************************************************************
    use prec_def
    use size_def,  only: ncells,ncellsD
    use flow_def,  only: u,v,h,rhow,grav
    use const_def, only: twopi
    use der_lib,   only: dx2d,dy2d
    use der_def,   only: gow
    use diag_def,  only: msg,msg2,msg3,msg4,msg5,msg6,msg7,msg8,msg9,msg10,msg11,msg12,msg13
    use diag_lib,  only: diag_print_error
    use wave_flowgrid_def, only: wavstrx,wavstry,whgt,wper,wlen,wunitx,wunity
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
      if(isnankind(wavstrx(i)) .or. isnankind(wavstry(i)) .or. abs(wavstrx(i))>10.0 .or. abs(wavstry(i))>10.0)then
        write(msg,*)  'ERROR: Found NaN'
        write(msg2,*)  'i:          ',i
        write(msg3,*) 'h(i):       ',h(i) 
        write(msg4,*) 'wavstrx(i): ',wavstrx(i)
        write(msg5,*) 'wavstry(i): ',wavstry(i)
        write(msg6,*) 'whgt(i):    ',whgt(i)
        write(msg7,*) 'wunitx(i):  ',wunitx(i)
        write(msg8,*) 'wunity(i):  ',wunity(i)
        write(msg9,*) 'wlen(i):    ',wlen(i)
        write(msg10,*) 'wper(i):    ',wper(i)
        write(msg11,*) 'Sxx(i):     ',Sxx(i)
        write(msg12,*) 'Sxy(i):     ',Sxy(i)
        write(msg13,*) 'Syy(i):     ',Syy(i)
        call diag_print_error(msg, msg2, msg3, msg4, msg5, msg6, msg7, msg8, msg9, msg10, msg11, msg12, msg13)
      endif
    enddo
    
    return
    end subroutine wave_rad_stress_fl