!==============================================================================
module diag_lib
!==============================================================================    
    implicit none
    
contains    

!******************************************************************         
    subroutine diag_print_var(i,neq)
! Prints the model variables to screen and debug file
!******************************************************************       
    use geo_def
    use flow_def
    use fric_def
    use diag_def
    use comvarbl, only: niter
    use wave_flowgrid_def
    use sed_def
    use cms_def, only: noptset
    use const_def, only: deg2rad
    use prec_def
    implicit none
    integer :: i,j,k,neq,ks,junit(2)

    junit(1)=6
    junit(2)=dgunit
        
    open(dgunit,file=dgfile,access='append') 
    
    do j=1,2
      write(junit(j),*)
      selectcase(neq)
        case(1)
          write(junit(j),*) 'ERROR in solving the pressure equation'
        case(2)
          write(junit(j),*) 'ERROR in solving the U-momentum equation'
        case(3)
          write(junit(j),*) 'ERROR in solving the V-momentum equation'
        case(4)
          write(junit(j),*) 'ERROR in solving the sediment transport equation'
        case(5)
          write(junit(j),*) 'ERROR in solving the salinity transport equation'      
        case(6)
          write(junit(j),*) 'ERROR in solving the temperature transfer equation'
      endselect
      write(junit(j),*) 'Flow Iteration:',niter
      if(allocated(mapid))then
        write(junit(j),*) 'Cell ID:',mapid(i)
      else
        write(junit(j),*) 'Cell ID:',i  
      endif
      if(igridtype==0)then
        write(junit(j),*)'Global Coordinates: ',xc(i),yc(i)     
        write(junit(j),*)     ' i: ',icol(i),', j: ',irow(i)
      else
        write(junit(j),*)'Global Coordinates: ',x(i),y(i)  
      endif
    enddo    
    
    do j=1,2
      write(junit(j),*) 'iwet(i),iwet1(i)'  
      write(junit(j),*)  iwet(i),iwet1(i)  
      write(junit(j),*) 'u(i),v(i),uv(i),u1(i),v1(i)'
      write(junit(j),*)  u(i),v(i),uv(i),u1(i),v1(i)     
      write(junit(j),*) 'zb(i),p(i)*gravinv,pp(i),h(i),h1(i)'
      write(junit(j),*)  zb(i),p(i)*gravinv,pp(i),h(i),h1(i)
      write(junit(j),*) 'cell2cell(1:ncface(i),i)'
      write(junit(j),*)  cell2cell(1:ncface(i),i)
      write(junit(j),*) 'iwet(cell2cell(1:ncface(i),i))'
      write(junit(j),*)  iwet(cell2cell(1:ncface(i),i))
      write(junit(j),*) 'pp(cell2cell(1:ncface(i),i)'
      write(junit(j),*)  pp(cell2cell(1:ncface(i),i))
      write(junit(j),*) 'h(cell2cell(1:ncface(i),i)'
      write(junit(j),*)  h(cell2cell(1:ncface(i),i))
      write(junit(j),*) 'p(cell2cell(1:ncface(i),i)'
      write(junit(j),*)  p(cell2cell(1:ncface(i),i))
      write(junit(j),*) '(flux(k,i),k=1,ncface(i))'
      write(junit(j),*)  (flux(k,i),k=1,ncface(i))
      write(junit(j),*) '(acoef(k,i),k=1,ncface(i))'
      write(junit(j),*)  (acoef(k,i),k=1,ncface(i))
      write(junit(j),*) '(hk(k,i),k=1,ncface(i))'
      write(junit(j),*)  (hk(k,i),k=1,ncface(i))
      write(junit(j),*) '(dsxy(k,i),k=1,ncface(i))'
      write(junit(j),*)  (dsxy(k,i),k=1,ncface(i))
      write(junit(j),*) '(visk(k,i),k=1,ncface(i))'
      write(junit(j),*)  (visk(k,i),k=1,ncface(i))
      write(junit(j),*) 'su(i),sp(i)'
      write(junit(j),*)  su(i),sp(i)                 
      write(junit(j),*) 'cfrict(i)'   
      write(junit(j),*)  cfrict(i)  
      write(junit(j),*) 'vis(i),uelwc(i)'  
      write(junit(j),*)  vis(i),uelwc(i)      
      write(junit(j),*) 'dpx(i),dpy(i)'
      write(junit(j),*)  dpx(i),dpy(i)
      write(junit(j),*) 'dppx(i),dppy(i)'
      write(junit(j),*)  dppx(i),dppy(i)
      write(junit(j),*) 'dux(i),duy(i),dvx(i),dvy(i)'
      write(junit(j),*)  dux(i),duy(i),dvx(i),dvy(i)
    enddo
    
    if(neq==4)then !Sediment      
      do j=1,2
        write(junit(j),*) 'Sediment Iteration:',niter
        if(singlesize)then
          write(junit(j),*) 'alphat(i),dzb(i)'  
          write(junit(j),*)  alphat(i),dzb(i)   
          write(junit(j),*) 'wsfall(1),varsigma(i,1),Sb(i,1)'
          write(junit(j),*)  wsfall(1),varsigma(i,1),Sb(i,1)
          write(junit(j),*) 'btk(i,1),btk1(i,1)'
          write(junit(j),*)  btk(i,1),btk1(i,1)
          write(junit(j),*) 'CtstarP(i,1),Ctkstar(i,1)'
          write(junit(j),*)  CtstarP(i,1),Ctkstar(i,1)
          write(junit(j),*) 'Ctk(i,1),Ctk1(i,1)'
          write(junit(j),*)  Ctk(i,1),Ctk1(i,1)
        else
          write(junit(j),*) 'alphat(i)'  
          write(junit(j),*)  alphat(i)
          write(junit(j),*) '(wsfall(ks),ks=1,nsed)'  
          write(junit(j),*)  (wsfall(ks),ks=1,nsed)
          write(junit(j),*) '(Sb(i,ks),ks=1,nsed)'
          write(junit(j),*)  (Sb(i,ks),ks=1,nsed)
          write(junit(j),*) '(varsigma(i,ks),ks=1,nsed)'
          write(junit(j),*)  (varsigma(i,ks),ks=1,nsed)
          write(junit(j),*) '(CtstarP(i,ks),ks=1,nsed)'
          write(junit(j),*)  (CtstarP(i,ks),ks=1,nsed)
          write(junit(j),*) '(btk(i,ks),ks=1,nsed)'
          write(junit(j),*)  (btk(i,ks),ks=1,nsed)
          write(junit(j),*) '(btk1(i,ks),ks=1,nsed)'
          write(junit(j),*)  (btk1(i,ks),ks=1,nsed)
          write(junit(j),*) '(Ctk(i,ks),ks=1,nsed)'
          write(junit(j),*)  (Ctk(i,ks),ks=1,nsed)
          write(junit(j),*) '(Ctk1(i,ks),ks=1,nsed)'
          write(junit(j),*)  (Ctk1(i,ks),ks=1,nsed)
          write(junit(j),*) 'dzb(i)'
          write(junit(j),*)  dzb(i)   
          write(junit(j),*) '(dzbk(i,ks),ks=1,nsed)'
          write(junit(j),*)  (dzbk(i,ks),ks=1,nsed)
          write(junit(j),*) '(pbk(i,ks,1),ks=1,nsed)'
          write(junit(j),*)  (pbk(i,ks,1),ks=1,nsed)
          write(junit(j),*) '(pbk(i,ks,2),ks=1,nsed)'
          write(junit(j),*)  (pbk(i,ks,2),ks=1,nsed)
          write(junit(j),*) '(pbk1(i,ks),ks=1,nsed)'
          write(junit(j),*)  (pbk1(i,ks),ks=1,nsed)
          write(junit(j),*) 'db(i,1),db1(i,1)'
          write(junit(j),*)  db(i,1),db1(i,1)
          write(junit(j),*) 'db(i,2),db1(i,2)'
          write(junit(j),*)  db(i,2),db1(i,2)
          write(junit(j),*) '(rsk(i,ks),ks=1,nsed)'
          write(junit(j),*)  (rsk(i,ks),ks=1,nsed)  
        endif
      enddo
    endif
    
    if(noptset>=3)then
      do j=1,2  
        write(junit(j),*) 'Whgt(i),Wper(i),Wlen(i)'  
        write(junit(j),*)  Whgt(i),Wper(i),Wlen(i)
        write(junit(j),*) 'wavestrx(i),wavestry(i)'
        write(junit(j),*)  wavestrx(i),wavestry(i)
        write(junit(j),*) 'wunitx(i),wunity(i)'
        write(junit(j),*)  wunitx(i),wunity(i)
        write(junit(j),*) 'Worb(i),Worbrep(i)'
        write(junit(j),*)  Worb(i),Worbrep(i)
        write(junit(j),*) 'us(i),vs(i)'
        write(junit(j),*)  us(i),vs(i)
      enddo          
    endif

    close(dgunit)
    
    return
    endsubroutine diag_print_var 

!******************************************************************************    
    subroutine diag_print_error(msg1,msg2,msg3,msg4,msg5,msg6,msg7,msg8,msg9)
! Prints a error message(s) to the screen and diagnostic file    
! written by Alex Sanchez, USACE-CHL       
!******************************************************************************
    use diag_def, only: dgfile,dgunit
    implicit none
    !Input/Output
    character(len=*),intent(in) :: msg1
    character(len=*),intent(in),optional :: msg2,msg3,msg4,msg5,msg6,msg7,msg8,msg9
    !Internal variables
    integer :: iunit(2),i
    logical :: dgopen
    
    iunit = (/6, dgunit/)  !Note: 6 is the screen
111 format(1x,A)
222 format(1x,A,A)    
    inquire(file=dgfile,opened=dgopen)
    if(.not.dgopen) open(dgunit,file=dgfile,access='append') 
    do i=1,2
      write(iunit(i),*)
      write(iunit(i),222) 'ERROR: ',trim(msg1)
      if(present(msg2)) write(iunit(i),111) trim(msg2)
      if(present(msg3)) write(iunit(i),111) trim(msg3)
      if(present(msg4)) write(iunit(i),111) trim(msg4)
      if(present(msg5)) write(iunit(i),111) trim(msg5)
      if(present(msg6)) write(iunit(i),111) trim(msg6)
      if(present(msg7)) write(iunit(i),111) trim(msg7)
      if(present(msg8)) write(iunit(i),111) trim(msg8)
      if(present(msg9)) write(iunit(i),111) trim(msg9)
    enddo
    close(dgunit)  
    write(*,111) 'Press any key to continue.'
    read(*,*)
    stop
    
    return
    endsubroutine diag_print_error
    
!******************************************************************************    
    subroutine diag_print_warning(msg1,msg2,msg3,msg4,msg5,msg6,msg7,msg8)
! Prints a warning message(s) to the screen and diagnostic file    
! written by Alex Sanchez, USACE-CHL    
!******************************************************************************
    use diag_def, only: dgfile,dgunit
    implicit none
    !Input/Output
    character(len=*),intent(in) :: msg1
    character(len=*),intent(in),optional :: msg2,msg3,msg4,msg5,msg6,msg7,msg8
    !Internal variables
    integer :: iunit(2),i
    logical :: dgopen
    
    iunit = (/6, dgunit/) !Note: 6 is the screen

111 format(1x,A)
222 format(1x,A,A)
    inquire(file=dgfile,opened=dgopen)
    if(.not.dgopen) open(dgunit,file=dgfile,access='append') 
    do i=1,2
      write(iunit(i),*)
      write(iunit(i),222) 'WARNING: ',trim(msg1)
      if(present(msg2)) write(iunit(i),111) trim(msg2)
      if(present(msg3)) write(iunit(i),111) trim(msg3)
      if(present(msg4)) write(iunit(i),111) trim(msg4)
      if(present(msg5)) write(iunit(i),111) trim(msg5)
      if(present(msg6)) write(iunit(i),111) trim(msg6)
      if(present(msg7)) write(iunit(i),111) trim(msg7)
      if(present(msg8)) write(iunit(i),111) trim(msg8)
      !write(iunit(i),*)
    enddo
    if(.not.dgopen) close(dgunit)
    
    return
    endsubroutine diag_print_warning   
    
!******************************************************************************    
    subroutine diag_print_message(msg1,msg2,msg3,msg4,msg5,msg6,msg7,msg8)
! Prints a message(s) to the screen and diagnostic file    
! written by Alex Sanchez, USACE-CHL    
!******************************************************************************
    use diag_def, only: dgfile,dgunit
    implicit none
    !Input/Output
    character(len=*),intent(in) :: msg1
    character(len=*),intent(in),optional :: msg2,msg3,msg4,msg5,msg6,msg7,msg8
    !Internal variables
    integer :: iunit(2),i
    logical :: dgopen
    
    iunit = (/6, dgunit/)
111 format(1x,A) 
    inquire(file=dgfile,opened=dgopen)
    if(.not.dgopen) open(dgunit,file=dgfile,access='append') 
    do i=1,2
      write(iunit(i),111) trim(msg1)
      if(present(msg2)) write(iunit(i),111) trim(msg2)
      if(present(msg3)) write(iunit(i),111) trim(msg3)
      if(present(msg4)) write(iunit(i),111) trim(msg4)
      if(present(msg5)) write(iunit(i),111) trim(msg5)
      if(present(msg6)) write(iunit(i),111) trim(msg6)
      if(present(msg7)) write(iunit(i),111) trim(msg7)
      if(present(msg8)) write(iunit(i),111) trim(msg8)
    enddo
    if(.not.dgopen) close(dgunit)  
    
    return
    endsubroutine diag_print_message      
    
endmodule diag_lib
    