!*************************************************************************
    subroutine sed_dry
! Erosion of dry cells 
!
! Description:
!   If a wet cell has erosion and is next to a dry cell, than a part of 
!   the erosion volume is transferred to the neighboring dry cell. 
!
! Author: Alex Sanchez, USACE-CHL
!         Mitchell Brown, USACE-CHL
!*************************************************************************
#include "CMS_cpp.h"
    use sed_def
    use size_def
    use flow_def
    use geo_def
    use diag_def
    use diag_lib
    use prec_def
    use const_def, only: pi
    implicit none
    integer :: i,k,nck,ndry,nwet,inddry(nmaxfaces),ndryerode
    integer :: isdryerode(ncellsD) 
    real(ikind) :: volnck,volsum,volorig,volnew,volbal,dzbdrymin
    real(ikind) :: wdry(nmaxfaces),wdrysum,slope
        
#ifdef DIAG_MODE
    do i=1,ncells
      if(iwet(i)==0 .and. abs(dzb(i))>1.0e-13)then
          call diag_print_error('Bed change at dry cell')
      endif
    enddo
#endif

    nwet = 0
    isdryerode = 0
    do i=1,ncells
      if(dzb(i)>-1.0e-15 .or. iwet(i)==0) cycle !Skip all wet or depositional cases

      !Erosion of wet cell
      ndry = 0
      inddry = 0
      wdry = 0.0
      do k=1,ncface(i)
        nck = cell2cell(k,i)
        slope = (zb(nck)-zb(i))/dc(k,i) !Always positive
        if(iwet(nck)==0 .and. nck<=ncells .and. slope>0.05)then    
          ndry = ndry + 1
          inddry(ndry) = nck
          !wdry(ndry) = slope-0.05
          !wdry(ndry) = 2.0*min(max(wdry(ndry),0.0),0.5)
          wdry(ndry) = (slope-erosdry%slopemin)/(erosdry%slopemax-erosdry%slopemin)
          wdry(ndry) = max(min(wdry(ndry),1.0),0.0)
        endif
      enddo
      if(ndry==0) cycle !Skip if no neighboring dry cells
      
      nwet = nwet + 1
      !Move erosion to neighboring dry cells
      volorig = dzb(i)*areap(i)
      !dzb(i) = dzb(i)*(1.0-erosdry%fac)
      wdry(1:ndry) = erosdry%fac*wdry(1:ndry)/real(ndry,kind=ikind)
      wdrysum = sum(wdry(1:ndry))
      dzb(i) = dzb(i)*(1.0-wdrysum)
      zb(i) = zb1(i) + dzb(i)
      volnew = dzb(i)*areap(i) !New erosional volume
      !volnck = erosdry%fac*volorig/real(ndry,kind=ikind) !Distribute volume evenly among dry cells
      volsum = 0.0
      do k=1,ndry
        nck = inddry(k)
        isdryerode(nck) = 1
        volnck = wdry(k)*volorig
        dzb(nck) = volnck/areap(nck) !dzb(nck) = 0 before this step since the cell is dry
        volsum = volsum + volnck
        zb(nck) = zb1(nck) + dzb(nck)
      enddo
!#ifdef DIAG_MODE
      volbal = volnew+volsum-volorig
      if(abs(volbal)>1.0e-6*abs(volorig))then
        write(msg2,*) ' i = ',i
        write(msg3,*) 'volorig,volnew,volsum,volbal = ',volorig,volnew,volsum,volbal
        call diag_print_warning('Problem distributing erosion at dry cells',msg2,msg3)
      !!else
      !!  if(allocated(mapid))then
      !!    write(msg2,*) '  id = ',mapid(i)
      !!  else
      !!    write(msg2,*) '  id = ',i 
      !!  endif
      !!  write(msg3,*) '  ndry = ',ndry
      !!  write(msg4,*) '  (inddry(k),k=1,ndry) = ',(inddry(k),k=1,ndry)
      !!  write(msg5,*) '  volorig = ',volorig
      !!  write(msg6,*) '  volnck = ',volnck
      !!  write(msg7,*) '  volnew = ',volnew
      !!  call diag_print_message('Erosion at dry cell',msg2,msg3,msg4,msg5,msg6,msg7)
      endif
!#endif
    enddo

    if(nwet>0)then
      ndryerode = 0
      dzbdrymin = 0.0
      do i=1,ncells
        if(isdryerode(i)==1)then
          ndryerode = ndryerode + 1 
          if(dzbdrymin>dzb(i))then
            dzbdrymin = dzb(i)  
          endif
        endif
      enddo
      write(msg, '(A)')           '  Erosion of dry cells:'
      write(msg2,'(A,I6)')        '    Number of wet cells with reduced erosion: ',nwet
      write(msg3,'(A,I6)')        '    Number of dry cells with erosion:         ',ndryerode
      write(msg4,'(A,1pe13.2,A)') '    Maximum erosion of dry cell:             ',dzbdrymin,' m'
      call diag_print_message(msg2,msg3,msg4)
    endif
    
    return
    end subroutine sed_dry
