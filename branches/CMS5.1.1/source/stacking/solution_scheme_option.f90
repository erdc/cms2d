!***********************************************************************
    subroutine solution_scheme_option()
!     by Chris Reed, R&R, June, 2012
!***********************************************************************
    use prec_def
    !use size_def_3D   !For 3D
!    use flow_def
!    use bnd2d
!    use qwsbnd2d
!    use met_def
!    use bicgsilu
!    use comvarbl
!	use cms_def
!    use fl_wavegrid
!	use wave_flgrid
    !use wavestress3D     !For 3D
    use comvarbl, only: ctlfile,nfsch   !STACK:
    use geo_def, only: telfile,igridtype
!    use bnd_def
!    use global_coordinate, only: azimuth_fl
!    use sedmod
!    use sal_def
!    use fric_def
!    use statsmod
!    use vegmod
!    use hotmod
!    use dermod
!    use outmod
!!    use q3dmod
    !!use prec_def
!    implicit none
!    integer :: i,ii,j,k,nck,ierr
!    real(ikind) :: fx,fkip
!    logical :: ok,foundcard
    
    implicit none
    integer :: k,ierr
    character*37 :: cardname
    logical :: foundcard,foundfile
    
    
!--- Set defaults -------------------------
    call geo_default   !Geospatial
    call flow_default  !Hydrodynamics
!    call fric_default  !Frictions
!    call bnd_default   !Boundary conditions
!    call hot_default   !Hot start
!    call sed_default   !Sediment transport 
!    call sal_default   !Salinity
!    call wind_default  !Wind
!    call veg_default   !Vegetation
!    call rol_default   !Roller
!    call out_default   !Output
!    call stat_default  !Statistics
!!    call q3d_default   !Quasi-3D 

!--- Read Card File -------------------------------------------------   
    !call fileparts(ctlfile,apath,aname,aext) 
    !astring = trim(aname) // '.' // trim(aext)
    !write(*,*) 'Reading CMS-Flow Card File: ',trim(astring)
    open(77,file=ctlfile,status='unknown')
	do k=1,1000
	  read(77,*,iostat=ierr) cardname
      if(ierr/=0) exit
	  if(cardname(1:14)=='END_PARAMETERS') exit
	  if(cardname(1:1)=='!' .or. cardname(1:1)=='#') cycle

 !!   open(77,file=ctlfile,status='unknown')
	!!do k=1,1000
	!!  read(77,*,err=171,end=171) cardname
	!!  if(cardname(1:14).eq.'END_PARAMETERS') exit
	!!  if(cardname(1:1).eq.'!' .or. cardname(1:1).eq.'#') cycle
      call geo_cards(cardname,foundcard);    if(foundcard) cycle     
	  !call ignr_cards(cardname,foundcard); if(foundcard) cycle
	  call flow_cards(cardname,foundcard);   if(foundcard) cycle 
	  !call bnd_cards(cardname,foundcard);  if(foundcard) cycle
	  !call fric_cards(cardname,foundcard); if(foundcard) cycle
	  !call hot_cards(cardname,foundcard);  if(foundcard) cycle
	  !call sed_cards(cardname,foundcard);  if(foundcard) cycle
	  !call sal_cards(cardname,foundcard);  if(foundcard) cycle
	  !call wind_cards(cardname,foundcard); if(foundcard) cycle  
	  !call veg_cards(cardname,foundcard);  if(foundcard) cycle
	  !call rol_cards(cardname,foundcard);  if(foundcard) cycle
	  !call out_cards(cardname,foundcard);  if(foundcard) cycle
!!	  call q3d_cards(cardname,foundcard);  if(foundcard) cycle
	  !call stat_cards(cardname,foundcard)
!	  if(.not.foundcard)then
!	    write(*,*) cardname,' not found'
!	    open(dgunit,file=dgfile,access='append') 
!	    write(dgunit,*) cardname,' not found'
!	    close(dgunit)
!	  endif
	enddo
171 close(77)   
    

    inquire(file=telfile,exist=foundfile)
    if(foundfile) igridtype = 1  
    
!    write(*,*)
!    write(*,*)"From Solution_Scheme_option"
!    write(*,*)'nfsch =',nfsch
!    write(*,*)'igridtype = ',igridtype
!    write(*,*)
    

    return
    endsubroutine    
