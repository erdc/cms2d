!***************************************************************************   
!***************************************************************************   
    subroutine struct_default
! Sets the default parameters for structure variables
! written by Alex Sanchez, USACE-CHL; Weiming Wu, NCCHE
!***************************************************************************
    use struct_def
    implicit none
    
    numweir = 0        ! Weirs
    numculvert = 0     ! Culverts
    numtidegate = 0    ! Tidal Gate
    nrubmoundcells = 0 ! Rubble Mounds

    return
    end subroutine struct_default
    
!***************************************************************************   
    subroutine struct_cards(cardname,foundcard)
! Sets the default parameters for flow variables
! written by Alex Sanchez, USACE-CHL; Weiming Wu, NCCHE
!***************************************************************************        
    use struct_def
    use flow_def, only: viscos,grav
    use const_def, only: pi
    use prec_def
    implicit none    
    integer :: itg,idtg,maxtidegate,maxtidegateopt  !Wu
    integer :: iwr,idwr,maxweir,icv,irm,idrm,maxrubmound    !Wu
    integer :: ierr
    character(len=37) :: cardname
    logical :: foundcard    
    
    foundcard = .true.
    select case (cardname) 
      
!=== Tidal Gates ====================================            
      case('TIDE_GATE')         !Wu
        backspace(77)
        read(77,*) cardname,numtidegate    !No. of tide gates
        if(numtidegate>0)then
          allocate(ntidegate(0:numtidegate),orienttidegate(numtidegate),coeftidegate(numtidegate,2),  &
                   openhgttidegate(numtidegate),elevtidegate(numtidegate),methtidegate(numtidegate),  &
                   orienttgbay(numtidegate),orienttgsea(numtidegate),Qtottidegate(numtidegate))       
          backspace(77)
          read(77,*) cardname, numtidegate,(ntidegate(itg),itg=1,numtidegate)              
          ntidegate(0) = 0
          maxtidegate = ntidegate(numtidegate)
          allocate(idtidegate(maxtidegate),opentidegate(maxtidegate),coeftglateral(maxtidegate), &
                   qtidegate(maxtidegate),dqtgdzdown(maxtidegate),dqtgdzup(maxtidegate))
          backspace(77)
          read(77,*) cardname, numtidegate,(ntidegate(itg),itg=1,numtidegate),                         &
                     (idtidegate(idtg),idtg=1,maxtidegate),(coeftglateral(idtg),idtg=1,maxtidegate),   &
                     (orienttidegate(itg),coeftidegate(itg,1),coeftidegate(itg,2),openhgttidegate(itg),&
                     elevtidegate(itg),methtidegate(itg),itg=1,numtidegate)                 
          dqtgdzdown=0.0;  dqtgdzup=0.0   
          do itg=1,numtidegate
             orienttgsea(itg) = orienttidegate(itg)   !Define direction of sea side at local coordinate
             if(orienttidegate(itg)==1) then
                orienttgbay(itg) = 3
             elseif(orienttidegate(itg)==2) then
                orienttgbay(itg) = 4
             elseif(orienttidegate(itg)==3) then
                orienttgbay(itg) = 1
             elseif(orienttidegate(itg)==4) then
                orienttgbay(itg) = 2
             endif
          enddo
        endif

      case('TIDE_GATE_OPERATION_SCHEDULE')         !Wu     !Tide gates operation schedule:
        allocate(mtidegateopt(numtidegate),ntidegateopt(0:numtidegate))
        backspace(77)
        read(77,*) cardname,(mtidegateopt(itg),itg=1,numtidegate),(ntidegateopt(itg),itg=1,numtidegate)              
        ntidegateopt(0) = 0
        maxtidegateopt = ntidegateopt(numtidegate)
        allocate(tidegateopt(maxtidegateopt))
        backspace(77)
        read(77,*) cardname,(mtidegateopt(itg),itg=1,numtidegate),(ntidegateopt(itg),itg=1,numtidegate),  &     
                            (tidegateopt(idtg),idtg=1,maxtidegateopt)        

      case('TIDE_GATE_BEGIN')
        call tide_gate_block
        
      case('TIDEGATE_STRUCT_BEGIN')
        call new_tg_block

!=== Weirs ====================================            
      case('WEIR')         !Wu
        backspace(77)
        read(77,*) cardname,numweir    !No. of weirs
        if(numweir>0)then
          allocate(nweir(0:numweir),orientweir(numweir),iweirtype(numweir),methweir(numweir),  &
                   coefweir(numweir,2),elevweir(numweir),orientweirbay(numweir),               &
                   orientweirsea(numweir),Qtotweir(numweir) )
          backspace(77)
          read(77,*) cardname, numweir,(nweir(iwr),iwr=1,numweir)              
          nweir(0) = 0
          maxweir = nweir(numweir)
          allocate(idweir(maxweir),qweir(maxweir),coefweirlateral(maxweir),dqweirdzdown(maxweir),dqweirdzup(maxweir)) 
          backspace(77)
          read(77,*) cardname, numweir,(nweir(iwr),iwr=1,numweir),(idweir(idwr),idwr=1,maxweir),  &
                     (coefweirlateral(idwr),idwr=1,maxweir),                                      &
                     (orientweir(iwr),iweirtype(iwr),coefweir(iwr,1),coefweir(iwr,2),             &
                      elevweir(iwr),methweir(iwr),iwr=1,numweir)              
          dqweirdzdown = 0.0;  dqweirdzup = 0.0   
          do iwr=1,numweir
             orientweirsea(iwr) = orientweir(iwr)   !Define direction of sea side at local coordinate
             if(orientweir(iwr)==1) then
                orientweirbay(iwr) = 3
             elseif(orientweir(iwr)==2) then
                orientweirbay(iwr) = 4
             elseif(orientweir(iwr)==3) then
                orientweirbay(iwr) = 1
             elseif(orientweir(iwr)==4) then
                orientweirbay(iwr) = 2
             endif
          enddo
        endif

      case('WEIR_BEGIN')         !Ver 5.0 and before
        call weir_block
        
      case('WEIR_STRUCT_BEGIN')  !Ver 5.1 and later
        call new_weir_block

!=== Culverts ====================================            
      case('CULVERT')         !Wu
        backspace(77)
        read(77,*) cardname,numculvert    !No. of culverts
        if(numculvert>0)then
          allocate(idculvert(numculvert,2),iculverttype(numculvert),iculvertflap(numculvert),        &
                   culvertrad(numculvert),culvertwidth(numculvert),culvertheight(numculvert),        &
                   culvertelevbay(numculvert),culvertelevsea(numculvert),culvertlength(numculvert),  &
                   cvheadlossbayentr(numculvert),cvheadlossbayexit(numculvert),cvheadlossseaentr(numculvert),  &
                   cvheadlossseaexit(numculvert),culvertfrict(numculvert),culvertmann(numculvert),   &
                   qculvert(numculvert),dqcvdzdown(numculvert),dqcvdzup(numculvert),                 &
                   uvculvert(numculvert),angleculvertbay(numculvert),angleculvertsea(numculvert) )   
          backspace(77)
          read(77,*) cardname,numculvert,(idculvert(icv,1),idculvert(icv,2),                    &
                     iculverttype(icv),iculvertflap(icv),culvertwidth(icv),culvertheight(icv),  &
                     culvertelevbay(icv),culvertelevsea(icv),culvertlength(icv),                &
                     cvheadlossbayentr(icv),cvheadlossbayexit(icv),cvheadlossseaentr(icv),      &
                     cvheadlossseaexit(icv),culvertfrict(icv),culvertmann(icv),                 &
                     angleculvertbay(icv),angleculvertsea(icv),icv=1,numculvert) 
          do icv=1,numculvert
            culvertrad(icv)=0.5*culvertwidth(icv)
            angleculvertbay(icv)=angleculvertbay(icv)*pi/180.0
            angleculvertsea(icv)=angleculvertsea(icv)*pi/180.0
          enddo
          dqcvdzdown=0.0;  dqcvdzup=0.0 
        endif
      
      case('CULVERT_BEGIN')
        call culvert_block
        
      case('CULVERT_STRUCT_BEGIN')  !Ver 5.3.12 and later (SMS 13.4 and later)
        call new_culvert_block
    
!=== Rubble Mound Structures ====================================            
      rmblock=0
      case('RUBBLE_MOUND')         !Wu
        backspace(77)
        read(77,*) cardname,nrubmoundcells    !No. of rubble mound structures
        if(nrubmoundcells>0)then
          allocate(nrubmound(0:nrubmoundcells),rubmounddia(nrubmoundcells),rubmoundporo(nrubmoundcells),   & 
                   rubmounda(nrubmoundcells),rubmoundb(nrubmoundcells),methrubmoundab(nrubmoundcells))   
          backspace(77)
          read(77,*) cardname,nrubmoundcells,(nrubmound(irm),irm=1,nrubmoundcells)              
          nrubmound(0) = 0
          maxrubmound = nrubmound(nrubmoundcells)
          backspace(77)
          read(77,*) cardname,nrubmoundcells,(nrubmound(irm),irm=1,nrubmoundcells),  &
                     (idrubmound(idrm),idrm=1,maxrubmound),                    &
                     (rubmounddia(irm),rubmoundporo(irm),methrubmoundab(irm),irm=1,nrubmoundcells)   
          do irm=1,nrubmoundcells                   !coefficients a and b 
            if(methrubmoundab(irm)==1) then     !Sidiropoulou et al. (2007)
              rubmounda(irm)=0.00333*rubmounddia(irm)**(-1.5)  *rubmoundporo(irm)**0.06         
              rubmoundb(irm)=0.194  *rubmounddia(irm)**(-1.265)*rubmoundporo(irm)**(-1.14)
            elseif(methrubmoundab(irm)==2) then  !Kadlec and Knight (1996)
              rubmounda(irm)=255.0*viscos*(1.0-rubmoundporo(irm))/grav/rubmoundporo(irm)**3.7/rubmounddia(irm)**2         
              rubmoundb(irm)=2.0*(1.0-rubmoundporo(irm))/grav/rubmoundporo(irm)**3/rubmounddia(irm)
            elseif(methrubmoundab(irm)==3) then   !Ward (1964)
              rubmounda(irm)=360.0*viscos/grav/rubmounddia(irm)**2         
              rubmoundb(irm)=10.44/grav/rubmounddia(irm)
            endif
          enddo
        endif
        
      case('RUBBLE_MOUND_ID_DATASET')
        backspace(77)
        read(77,*) cardname, arubmoundfile, arubmoundpath
          
      case('RUBBLE_MOUND_BEGIN')
        call read_rubble_mound(cardname,foundcard)        

      case default
        foundcard = .false.
        
    end select

    return
    end subroutine struct_cards

!****************************************************    
    subroutine struct_init
!****************************************************
#include "CMS_cpp.h"
    use geo_def,  only: zb,mapid   !(hli, 11/19/13)  
    use size_def, only: ncellsfull  !(hli, 11/19/13) 
    use flow_def, only: grav  !(hli, 11/19/13) 
    use struct_def
    use flow_def, only: viscos
    use comvarbl, only: SMS_ver
#ifdef XMDF_IO
    use in_xmdf_lib, only: readscalh5
#endif
    use in_lib, only: readscalTxt
    implicit none
    integer :: i,i1,itg,im,iwr,icv,irm,kk,ierr   !hli(11/19/13)
    character(len=10) :: aext    
    integer :: j,idum
    logical :: exists

    !Allocate for full datasets and initialize
    allocate(permeability(ncellsfull))
    allocate(rockdiam(ncellsfull))
    allocate(structporo(ncellsfull))
    allocate(structbaseD(ncellsfull))
    allocate(structmeth(ncellsfull))
    permeability=0.0; rockdiam=0.0; structporo=0.0; structbaseD=0.0; structmeth=0

    inquire(file=arubmoundfile, exist=exists)
    if (exists) then   !First check if an ID file was specified. 
      if (SMS_ver >= 13.3) then
        call test_new_struct_init_v2
      else
        call new_struct_init
      endif  
    else                              !If not this is the old way.  
      nrubmoundcells = 0
      do i=1,irubmound
        nrubmoundcells = nrubmoundcells + RM_struct(i)%ncells 
      enddo

      allocate(nrubmound(0:nrubmoundcells),rubmounddia(nrubmoundcells),rubmoundporo(nrubmoundcells), & 
        rubmounda(nrubmoundcells),rubmoundb(nrubmoundcells),methrubmoundab(nrubmoundcells), &
        rubmoundbotm(nrubmoundcells),idrubmound(nrubmoundcells))  

      do irm=1,nrubmoundcells
        nrubmound(irm)=irm
      enddo

      nrubmound(0) = 0
      kk=0
      do irm=1,irubmound
        do j=1,RM_struct(irm)%ncells
          kk=kk+1
          idrubmound(kk)     = RM_struct(irm)%cells(j)         !don't convert, already the SMS Cell ID number
          rubmounddia(kk)    = RM_struct(irm)%rockdia_const
          rubmoundporo(kk)   = RM_struct(irm)%structporo_const
          rubmoundbotm(kk)   = RM_struct(irm)%structbaseD_const
          methrubmoundab(kk) = RM_struct(irm)%rubmoundmeth
        enddo
      enddo
    endif  
    
    if (rmblock == 1 .and. nrubmoundcells > 0) then 
      do irm=1,nrubmoundcells     
        select case (methrubmoundab(irm))
        case (1)             !Sidiropoulou et al. (2007)
          rubmounda(irm)=0.00333*rubmounddia(irm)**(-1.5) *rubmoundporo(irm)**0.06         
          rubmoundb(irm)=0.194  *rubmounddia(irm)**(-1.265)*rubmoundporo(irm)**(-1.14)
        case (2)             !Kadlec and Knight (1996)
          rubmounda(irm)=255.0*viscos*(1.0-rubmoundporo(irm))/grav/rubmoundporo(irm)**3.7/rubmounddia(irm)**2 
          rubmoundb(irm)=2.0*(1.0-rubmoundporo(irm))/grav/rubmoundporo(irm)**3/rubmounddia(irm)
        case (3)             !Ward (1964)
          rubmounda(irm)=360.0*viscos/grav/rubmounddia(irm)**2         
          rubmoundb(irm)=10.44/grav/rubmounddia(irm)
        end select
      enddo
    endif
    
!   (hli, 12/10/12)
!--- Map structure ids from full to active grid -----
    !Tidal Gates
    if(itidegate .gt. 0) call tg_init   !If new multi block tidegate called, then initialize
    if(numtidegate>0) call map_cell_full2active(ntidegate(numtidegate),idtidegate)

    !Weirs    
    if(iweir .gt. 0) call weir_init     !If new multi block weir called, then initialize
    if(numweir>0) call map_cell_full2active(nweir(numweir),idweir)
    
    !Culverts
    if(iculvert .gt. 0) call culvert_init
    if(numculvert>0)then
      call map_cell_full2active(numculvert,idculvert(:,1))
      call map_cell_full2active(numculvert,idculvert(:,2))    
    endif
    
    !Rubble Mounds   
    if(nrubmoundcells>0) call map_cell_full2active(nrubmound(nrubmoundcells),idrubmound)

!--- Modify Bathymetry at structures --------------------    
    !Tidal Gate
    do itg=1,numtidegate
      do im=ntidegate(itg-1)+1,ntidegate(itg)
        i=idtidegate(im)
        zb(i)=max(zb(i),elevtidegate(itg))  
      enddo
    enddo

    !Weir
    do iwr=1,numweir
      do im=nweir(iwr-1)+1,nweir(iwr)
        i=idweir(im)
        zb(i)=max(zb(i),elevweir(iwr))  
      enddo
    enddo

    !Culvert
    do icv=1,numculvert
      i =idculvert(icv,1)
      i1=idculvert(icv,2)
      zb(i)=min(zb(i),culvertelevbay(icv))  
      zb(i1)=min(zb(i1),culvertelevsea(icv))  
    enddo

!---- Rubble Mound Structure ------------
    do irm=1,nrubmoundcells
      do im=nrubmound(irm-1)+1,nrubmound(irm)
        i=idrubmound(im)
        zb(i)=-rubmoundbotm(irm)   ! minus the depth below reference datum, i.e. still water level 
      enddo
    enddo
    
    return
    end subroutine struct_init

    
!*****************************************************************************************
    subroutine new_struct_init
!*****************************************************************************************
#include "CMS_cpp.h"
    use geo_def,  only: mapid
    use size_def, only: ncellsfull
    use diag_lib, only: diag_print_error
    use in_lib,   only: readscalTxt
    use struct_def
#ifdef XMDF_IO
    use in_xmdf_lib, only: readscalh5
#endif
    implicit none
    
    logical :: use_dataset, exists
    character(len=10)  :: aext
    character(len=100) :: afile
    integer :: ierr,i,j,irm,max,num,kk, ival

    !Read permeable cell info from dataset
    call fileext(trim(arubmoundfile),aext)      
    select case (aext)
    case('h5')
#ifdef XMDF_IO
      call readscalh5(arubmoundfile,arubmoundpath,permeability,ierr)
#endif
    case('txt')
      call readscalTxt(arubmoundfile,permeability,ierr)
    end select
      
    !Count total number of permeable cells
    nrubmoundcells=0
    if(ierr<0) call diag_print_error("Cannot open Rubble Mound ID dataset")
    do i=1,ncellsfull
      if(permeability(i).gt.0.0) then
        nrubmoundcells=nrubmoundcells+1
      endif
    enddo

    !Allocate rubble mound arrays to correct size
    if(nrubmoundcells.gt.0)then
      allocate(nrubmound(0:nrubmoundcells))
      allocate(rubmounda(nrubmoundcells),rubmoundb(nrubmoundcells),methrubmoundab(nrubmoundcells))
      allocate(rubmoundbotm(nrubmoundcells),rubmounddia(nrubmoundcells),rubmoundporo(nrubmoundcells))
      allocate(idrubmound(nrubmoundcells)) 
      nrubmound = 0; rubmounda = 0.0; rubmoundb = 0.0; methrubmoundab = 0
      rubmoundbotm = 0.0; rubmounddia = 0.0; rubmoundporo = 0.0; idrubmound = 0 
    endif
        
    !Assign each rubble mound cell an index number
    do irm=1,nrubmoundcells
      nrubmound(irm) = irm  
    enddo

    !Allocate each rubble mound structures 'cells' array.
    max = maxval(permeability)
    do i = 1, max
      num = count(permeability==float(i))
      if (.not.allocated(RM_struct(i)%cells)) allocate(RM_struct(i)%cells(num))
      RM_struct(i)%ncells = num
    enddo
    
    !Assign cells ids to array.
    nrubmound(0) = 0
    do i=1,ncellsfull
      if(permeability(i).gt.0.0) then
        ival = int(permeability(i))
        RM_struct(ival)%inc = RM_struct(ival)%inc + 1
        RM_struct(ival)%cells(RM_struct(ival)%inc) = mapid(i)
      endif
    enddo
    !Will assign cells to idrubmound later

!Read rock diameter from dataset
    !if any of the rubble mound blocks have a dataset for rock diameter, 
    !  read it in. Go through and reset values to constants later.
    use_dataset = .false.
    do i=1,irubmound
      afile = RM_struct(i)%rockdia_dset(1)
      inquire(file=trim(afile), exist=exists)
      if (exists) use_dataset = .true.
    enddo
    if (use_dataset) then 
      call fileext(trim(arockdiamfile),aext)      
      select case (aext)
      case('h5')
#ifdef XMDF_IO
        call readscalh5(arockdiamfile,arockdiampath,rockdiam,ierr)
#endif
      case('txt')
        call readscalTxt(arockdiamfile,rockdiam,ierr)
      end select
          
      kk=0
      do i=1,ncellsfull
        if(rockdiam(i).gt.0.0) then
          kk=kk+1
          rubmounddia(kk)=rockdiam(i)
        endif
      enddo
    endif
    
!Read porosity from dataset
    !if any of the rubble mound blocks have a dataset for porosity, 
    !  read it in. Go through and reset values to constants later.
    use_dataset = .false.
    do i=1,irubmound
      afile = RM_struct(i)%structporo_dset(1)
      inquire(file=afile, exist=exists)
      if (exists) use_dataset = .true.
    enddo
    if (use_dataset) then
      call fileext(trim(astructporofile),aext)      
      select case (aext)
      case('h5')
#ifdef XMDF_IO
        call readscalh5(astructporofile,astructporopath,structporo,ierr)
#endif
      case('txt')
        call readscalTxt(astructporofile,structporo,ierr)
      end select
      
      kk=0
      do i=1,ncellsfull
        if(structporo(i).gt.0.0) then
          kk=kk+1
          rubmoundporo(kk)=structporo(i)
        endif
      enddo
    endif

!Read Base depth from dataset
    !if any of the rubble mound blocks have a dataset for base depth, 
    !  read it in. Go through and reset values to constants later.
    use_dataset = .false.
    do i=1,irubmound
      afile = RM_struct(i)%structbaseD_dset(1)
      inquire(file=afile, exist=exists)
      if (exists) use_dataset = .true.
    enddo
    if (use_dataset) then
      call fileext(trim(astructbaseDfile),aext)      
      select case (aext)
      case('h5')
#ifdef XMDF_IO
        call readscalh5(astructbaseDfile,astructbaseDpath,structbaseD,ierr)
#endif
      case('txt')
        call readscalTxt(astructbaseDfile,structbaseD,ierr)
      end select
         
      kk=0
      do i=1,ncellsfull
        if(structbaseD(i).gt.0.0) then
          kk=kk+1
          rubmoundbotm(kk)=structbaseD(i)
        endif
      enddo
    endif

!Read rubble mound method from dataset
    !if any of the rubble mound blocks have a dataset for method, 
    !  read it in. Go through and reset values to constants later.
    use_dataset = .false.
    do i=1,irubmound
      afile = RM_struct(i)%structmeth_dset(1)
      inquire(file=afile, exist=exists)
      if (exists) use_dataset = .true.
    enddo
    if (use_dataset) then
      call fileext(trim(astructmethfile),aext)      
      select case (aext)
      case('h5')
#ifdef XMDF_IO
        call readscalh5(astructmethfile,astructmethpath,structmeth,ierr)
#endif
      case('txt')
        call readscalTxt(astructmethfile,structmeth,ierr)
      end select
          
      kk=0
      do i=1,ncellsfull
        if(structmeth(i).gt.0.0) then
          kk=kk+1
          methrubmoundab(kk)=structmeth(i)
        endif
      enddo
    endif
    
!Now go through and set constants (only if present) as well as the name and cell id
    kk=0
    do irm=1,irubmound
      do j=1,RM_struct(irm)%ncells
        kk=kk+1
        idrubmound(kk)   = RM_struct(irm)%cells(j) !don't convert, already the SMS Cell ID number
        if(RM_struct(irm)%rubmoundmeth .gt. 0)        methrubmoundab(kk) = RM_struct(irm)%rubmoundmeth
        if(RM_struct(irm)%rockdia_const .gt. 0.0)     rubmounddia(kk)    = RM_struct(irm)%rockdia_const
        if(RM_struct(irm)%structporo_const .gt. 0.0)  rubmoundporo(kk)   = RM_struct(irm)%structporo_const
        if(RM_struct(irm)%structbaseD_const .gt. 0.0) rubmoundbotm(kk)   = RM_struct(irm)%structbaseD_const
      enddo
    enddo
    
    return
    end subroutine new_struct_init
    
!*****************************************************************************************
    subroutine test_new_struct_init_v2
#include "CMS_cpp.h"
    use geo_def,  only: mapid
    use size_def, only: ncellsfull
    use diag_lib, only: diag_print_error
    use in_lib,   only: readscalTxt
    use struct_def
#ifdef XMDF_IO
    use in_xmdf_lib, only: readscalh5
#endif
    implicit none
    
    logical :: use_dataset, exists
    character(len=10)  :: aext
    character(len=100) :: afile
    integer :: ierr,i,j,irm,max,num,kk, ival
    integer, allocatable :: temparray(:)

    !Read all datasets in from disk
    call fileext(trim(arubmoundfile),aext)
    select case(aext)
    case('h5')
#ifdef XMDF_IO
      if(arubmoundpath(1:5)/='Empty') call readscalh5(arubmoundfile,arubmoundpath,permeability,ierr)
      if(arockdiampath(1:5)/='Empty') call readscalh5(arubmoundfile,arockdiampath,rockdiam,ierr)
      if(astructporopath(1:5)/='Empty') call readscalh5(arubmoundfile,astructporopath,structporo,ierr)
      if(astructbaseDpath(1:5)/='Empty') call readscalh5(arubmoundfile,astructbaseDpath,structbaseD,ierr)
#endif
    case('txt')
      call readscalTxt(arubmoundfile,permeability,ierr)
      call readscalTxt(arockdiamfile,rockdiam,ierr)
      call readscalTxt(astructporofile,structporo,ierr)
      call readscalTxt(astructbaseDfile,structbaseD,ierr)
    end select

    allocate(temparray(ncellsfull))
    temparray=0
    
    nrubmoundcells = 0
    if(ierr < 0) call diag_print_error("Cannot open one of the Rubble Mound datasets")
    do i=1,ncellsfull
      if(permeability(i).gt.0.0) then
        nrubmoundcells = nrubmoundcells+1
        temparray(nrubmoundcells) = mapid(i)
      endif
    enddo
  
    !Allocate rubble mound arrays to correct size
    if(nrubmoundcells.gt.0)then
      allocate(nrubmound(0:nrubmoundcells))
      allocate(rubmounda(nrubmoundcells),rubmoundb(nrubmoundcells),methrubmoundab(nrubmoundcells))
      allocate(rubmoundbotm(nrubmoundcells),rubmounddia(nrubmoundcells),rubmoundporo(nrubmoundcells))
      allocate(idrubmound(nrubmoundcells)) 
      nrubmound = 0; rubmounda = 0.0; rubmoundb = 0.0; methrubmoundab = 0
      rubmoundbotm = 0.0; rubmounddia = 0.0; rubmoundporo = 0.0; idrubmound = 0 
    endif

    !Assign each rubble mound cell an index number and id value
    do irm=1,nrubmoundcells
      idrubmound(irm) = temparray(irm)
      nrubmound(irm) = irm
    enddo
    deallocate(temparray)
  
    !Allocate each rubble mound structures 'cells' array.
    max = maxval(permeability)
    do i = 1, max
      num = count(permeability==float(i))
      if (.not.allocated(RM_struct(i)%cells)) allocate(RM_struct(i)%cells(num))
      RM_struct(i)%ncells = num
    enddo
  
    kk=0
    do i=1,ncellsfull
      if(permeability(i).gt.0.0) then
        kk = kk + 1
        if(allocated(rockdiam))    rubmounddia(kk)    = rockdiam(i)
        if(allocated(structporo))  rubmoundporo(kk)   = structporo(i)
        if(allocated(structbaseD)) rubmoundbotm(kk)   = structbaseD(i)
        if(allocated(structmeth))  methrubmoundab(kk) = structmeth(i)
    
        !Assign Cell ids
        ival = int(permeability(i))
        RM_struct(ival)%inc = RM_struct(ival)%inc + 1
        RM_struct(ival)%cells(RM_struct(ival)%inc) = mapid(i)
      endif
    enddo
    
    kk=0
    do irm=1,irubmound
      do j=1,RM_struct(irm)%ncells
        kk=kk+1
        if(RM_struct(irm)%rubmoundmeth .gt. 0)        methrubmoundab(kk) = RM_struct(irm)%rubmoundmeth
        if(RM_struct(irm)%rockdia_const .gt. 0.0)     rubmounddia(kk)    = RM_struct(irm)%rockdia_const
        if(RM_struct(irm)%structporo_const .gt. 0.0)  rubmoundporo(kk)   = RM_struct(irm)%structporo_const
        if(RM_struct(irm)%structbaseD_const .gt. 0.0) rubmoundbotm(kk)   = RM_struct(irm)%structbaseD_const
      enddo
    enddo    
    
    return
    end subroutine test_new_struct_init_v2
    
    
!*****************************************************************************************
    subroutine struct_print()
! Prints the structures setup to the Screen and Diagnostic file
! Author: Alex Sanchez, USACE-CHL
!*****************************************************************************************
    use const_def
    use struct_def
    use diag_def,   only: dgfile,dgunit
    use prec_def
    use tool_def, only: vstrlz, adjustc
    
    implicit none
    integer :: i,j,icv,iwr,itg,iunit(2),ierr,ival
 
111 format(' ',A,T40,A)
112 format(' ',A,I0)
241 format(' ',A,T40,I0)
242 format(' ',A,T40,I0,A)
243 format(' ',A,T40,F5.2)
244 format(' ',A,T40,A,A)
245 format(' ',A,T40,A,T65,A)
    
    iunit = (/6, dgunit/)    
    open(dgunit,file=dgfile,access='append') 
    do i=1,2
      write(iunit(i),*)  
      if(numtidegate+numweir+numculvert+nrubmoundcells > 0)then
        write(iunit(i),111) '  Structures:'
      else
        write(iunit(i),111) '  Structures:','OFF'  
        cycle
      endif
      !---- Rubble Mound Structures -------
      !Added more information to diagnostic summary - MEB  09/25/23
      if (irubmound > 0) then
        write(iunit(i),241)   '    Number of Rubble Mounds:',irubmound 
        !if (use_dset) write(iunit(i),111)   '    Rubble Mound Dataset File:',trim(rm_struct(j)%rockdia_dset(1))
        do j=1,irubmound
          write(iunit(i),112) '      Rubble Mound Structure: ',j
          if (rm_struct(j)%name .ne. '') write(iunit(i),111) '        Name:',trim(rm_struct(j)%name) 
          write(iunit(i),241) '        Number of cells:',rm_struct(j)%ncells
          if (rm_struct(j)%structporo_const  .gt. 0.0) write(iunit(i),244) '        Porosity Constant:',      trim(vstrlz(rm_struct(j)%structporo_const,'(f0.2)'))
          if (rm_struct(j)%structbaseD_const .gt. 0.0) write(iunit(i),244) '        Base Depth Constant:',    trim(vstrlz(rm_struct(j)%structbaseD_const,'(f0.2)')),' m'
          if (rm_struct(j)%rockdia_const     .gt. 0.0) write(iunit(i),244) '        Rock Diameter Constant:', trim(vstrlz(rm_struct(j)%rockdia_const,'(f0.2)')),' m'
          if (rm_struct(j)%structporo_const  .eq. 0.0) write(iunit(i),245) '        Porosity Dataset:',       trim(rm_struct(j)%structporo_dset(2)), trim(rm_struct(j)%structporo_dset(1))
          if (rm_struct(j)%structbaseD_const .eq. 0.0) write(iunit(i),245) '        Base Depth Dataset:',     trim(rm_struct(j)%structbaseD_dset(2)),trim(rm_struct(j)%structbaseD_dset(1))
          if (rm_struct(j)%rockdia_const     .eq. 0.0) write(iunit(i),245) '        Rock Diameter Dataset:',  trim(rm_struct(j)%rockdia_dset(2)),    trim(rm_struct(j)%rockdia_dset(1))
          if (rm_struct(j)%rubmoundmeth      .eq. 0) then
            write(iunit(i),245) '        Method Dataset:', trim(rm_struct(j)%structmeth_dset(2)), trim(rm_struct(j)%structmeth_dset(1))
          else
            ival = rm_struct(j)%rubmoundmeth
            if (ival .eq. 1) write(iunit(i),111) '        Method:','Sidiropoulou et al. (2007)'
            if (ival .eq. 2) write(iunit(i),111) '        Method:','Kadlec and Knight (1996)'
            if (ival .eq. 3) write(iunit(i),111) '        Method:','Ward (1964)'
          endif
        enddo
      endif
      
!---- Tidal Gate --------------------
        ! '       |     Coefficients     |      Sea      |    Opening    |    Bottom    |               |  Schedule   |'
        ! '    ID |  Distr.   Bay   Sea  |  Orientation  |      Hgt      |     Elev     |     Method    |    Type     |'
655   format(I7,' |  ',F6.2,x,2(F6.2),'  |    ',A5,'     |   ',F6.2,'    |   ',F6.2,'   |   Approach ',I1,'  |',A12,'|')
      if(numtidegate>0)then
        write(iunit(i),241) '    Number of Tidal Gates:',numtidegate
        write(iunit(i),*)   ''
        write(iunit(i),111) '       |      Coefficients     |     Sea      |   Opening   |   Bottom   |               |  Schedule  |'
        write(iunit(i),111) '    ID |   Distr.   Bay   Sea  | Orientation  |     Hgt     |    Elev    |     Method    |    Type    |'
        write(iunit(i),111) '    ---------------------------------------------------------------------------------------------------'
        do itg=1,numtidegate
          write(iunit(i),655,iostat=ierr) itg, coeftglateral(itg),coeftidegate(itg,1), coeftidegate(itg,2),&
            trim(aorienttgsea(itg)),openhgttidegate(itg),elevtidegate(itg),methtidegate(itg),adjustc(aschedtype(itg))
        enddo
      endif
      continue

            
!---- Weir --------------------
        ! '       |             |    Coefficients     |     Sea      |    Crest    |'   
        ! '    ID |  Crest Type | Distr.   Bay   Sea  | Orientation  |  Elevation  |   Method'
656   format(I7,' |   ',A5,'    | ',F6.2,x,2(F6.2),'  |    ',A5,'    |  ',F6.2,'   | Approach ',I1)
      if(numweir>0)then
        write(iunit(i),241) '    Number of Weirs:',numweir
        write(iunit(i),*)   ''
        write(iunit(i),111) '       |            |     Coefficients     |     Sea     |   Crest   |'
        write(iunit(i),111) '    ID | Crest Type |  Distr.   Bay   Sea  | Orientation | Elevation |   Method   |'
        write(iunit(i),111) '    -------------------------------------------------------------------------------'
        do iwr=1,numweir
          write(iunit(i),656,iostat=ierr) iwr,aweirtype(iwr),      &
            coefweirlateral(iwr),coefweir(iwr,1),coefweir(iwr,2),  &
            aorientweirsea(iwr),elevweir(iwr),methweir(iwr)
        enddo
      endif
      
      ! '       |         |            |           |           |     Elevation     |      Exit       |      Entry      |     Friction       |  Outflow Angles ' 
      ! '    ID |  Type   |  Rad/Width |  Height   |  Length   |     Bay   Sea     |    Bay   Sea    |    Bay   Sea    |  Darcy  Mannings   |   Bay      Sea' 
657 format(I7,' |  ',A3,' |  'F6.2,'   | ',F5.2,'  | ',F5.2,'  | ',F5.2,x,F5.2,  ' | ',F4.2,x,F4.2,' | ',F4.2,x,F4.2,' | ',F5.2,3x,F5.2,'   | ',F6.2,2x,F6.2)
      !---- Culvert --------------------
      if(numculvert>0)then
        write(iunit(i),241) '    Number of Culverts:',numculvert
        write(iunit(i),111) '                                                          |   Loss Coefficients   |' 
        write(iunit(i),111) '       |      |           |        |        |  Elevation  |   Exit    |   Entry   |    Friction     | Outflow Angles |' 
        write(iunit(i),111) '    ID | Type | Rad/Width | Height | Length |  Bay   Sea  | Bay   Sea | Bay   Sea | Darcy  Mannings |   Bay    Sea   |' 
        write(iunit(i),111) '    ------------------------------------------------------------------------------------------------------------------'
        do icv=1,numculvert
          write(iunit(i),657,iostat=ierr) icv,aculverttype(icv),      &
             culvertwidth(icv),culvertheight(icv),culvertlength(icv), &
             culvertelevbay(icv),culvertelevsea(icv),                 &
             cvheadlossbayexit(icv),cvheadlossbayentr(icv),           &
             cvheadlossseaexit(icv),cvheadlossseaentr(icv),           &
             culvertfrict(icv),culvertmann(icv),                      &
             angleculvertbay(icv)*rad2deg,angleculvertsea(icv)*rad2deg
        enddo
        continue
      endif
    enddo
    
    close(dgunit)
    
    return
    end subroutine struct_print
    
!***********************************************************************
    subroutine struct_wetdry
! Judge Dry or Wet   
! by Weiming Wu, NCCHE, Oct. 2008
!***********************************************************************
    use geo_def, only: ncface,cell2cell,idirface,zb
    use flow_def, only: iwet,p,hmin,grav
    use comvarbl, only: timehrs,ctime
    use const_def, only: pi
    use struct_def
    use prec_def
    implicit none
    integer :: i,k,im,itg,ktg,iwr,nck
    real(ikind) :: ctimetg,wslbay,snodebay,wslsea,snodesea
    real(ikind) :: tgateoptstr,tgateoptperiod,tgateopendura
    
!--- Tide Gate --------------------------------------------------------
    if(numtidegate>0) opentidegate=0            !=1, open; =0, close
    do itg=1,numtidegate
      if(mtidegateopt(itg)==1)then              !1 for constant time period
        tgateoptstr   =tidegateopt(ntidegateopt(itg-1)+1)
        tgateoptperiod=tidegateopt(ntidegateopt(itg-1)+2)
        tgateopendura =tidegateopt(ntidegateopt(itg-1)+3)
        ctimetg=(timehrs-tgateoptstr)/tgateoptperiod
        ctimetg=(ctimetg-floor(ctimetg))*tgateoptperiod !Alex
        if(ctimetg<=tgateopendura) then
          do im=ntidegate(itg-1)+1,ntidegate(itg)
            opentidegate(im)=1
          enddo
        endif
      elseif(mtidegateopt(itg)==2)then          !2 for specified schedule
        do ktg=1,(ntidegateopt(itg)-ntidegateopt(itg-1))/2
          ctimetg=ctime/3600.0
          tgateoptstr  =tidegateopt(ntidegateopt(itg-1)+2*ktg-1)
          tgateopendura=tidegateopt(ntidegateopt(itg-1)+2*ktg  )
          if(ctimetg>=tgateoptstr .and. ctimetg<=(tgateoptstr+tgateopendura)) then
             do im=ntidegate(itg-1)+1,ntidegate(itg)
               opentidegate(im)=1
             enddo
             exit
          endif
        enddo
      elseif(mtidegateopt(itg)==3)then          !3 for flap gate
        wslbay=0.0
        snodebay=0.0
        wslsea=0.0
        snodesea=0.0
        do im=ntidegate(itg-1)+1,ntidegate(itg)
          i=idtidegate(im)
          do k=1,ncface(i)
            nck=cell2cell(k,i)
            if(idirface(k,i)==orienttgbay(itg)) then
              wslbay=wslbay+p(nck)/grav
              snodebay=snodebay+1.0
            endif
            if(idirface(k,i)==orienttgsea(itg)) then
              wslsea=wslsea+p(nck)/grav
              snodesea=snodesea+1.0
            endif
          enddo
        enddo
        wslbay=wslbay/snodebay
        wslsea=wslsea/snodesea
        if(wslbay>wslsea) then
          do im=ntidegate(itg-1)+1,ntidegate(itg)
            opentidegate(im)=1
          enddo
        endif
      elseif(mtidegateopt(itg)==4)then          !4 for no control 
        do im=ntidegate(itg-1)+1,ntidegate(itg)
          opentidegate(im)=1
        enddo
      endif
    enddo

    do itg=1,numtidegate
      do im=ntidegate(itg-1)+1,ntidegate(itg)
        i=idtidegate(im)
        if(opentidegate(im)==0) then
          iwet(i)=0  
        else
          iwet(i)=0
          if(methtidegate(itg)==2) then
            do k=1,ncface(i)
              nck=cell2cell(k,i)
              if(iwet(nck)==1.and.p(nck)/grav>zb(i)+hmin) iwet(i)=1
            enddo
          endif
        endif
      enddo
    enddo

!--- Weir -----------------------------------------------------------------
    do iwr=1,numweir
      do im=nweir(iwr-1)+1,nweir(iwr)
        i=idweir(im)
        iwet(i)=0
        if(methweir(iwr)==2)then 
          do k=1,ncface(i)
            nck=cell2cell(k,i)
            if(iwet(nck)==1.and.p(nck)/grav>zb(i)+hmin) iwet(i)=1
          enddo
        endif
      enddo
    enddo
    
    return
    end subroutine struct_wetdry
    
!***************************************************************************
    subroutine struct_uv
! Calculates and adds source/sink terms for u and v momentum equations    
! written by Weiming Wu and modified by Alex Sanchez, USACE-CHL
!***************************************************************************    
    use size_def
    use geo_def, only: ncface,cell2cell,llec2llec,signface,idirface,ds,areap
    use flow_def, only: iwet,h,h1,h2,p,dpx,dpy,u1,u2,v1,v2,uv,acoef,su,sv,spu,spv,spuv0,ssu0,ssv0,grav
    use struct_def
    use comvarbl, only: dtime,ntsch,ctsch1,ctsch2,relax
    use const_def, only: small,pi
    implicit none
    integer :: i,ii,k,im,imm,nck
    integer :: itg,iwr,irm,icv,i1,isRMS
    real(ikind) :: tgflux,wslbay,wslsea,snodesea,snodebay  
    real(ikind) :: val,fac1,fac2
    
!--- Tidal Gates ---------------------------------------------
    do itg=1,numtidegate
      if(methtidegate(itg)==1) then
        do im=ntidegate(itg-1)+1,ntidegate(itg)
          i=idtidegate(im)
          if(opentidegate(im)==1) then   !Gate open
            do k=1,ncface(i)
              nck=cell2cell(k,i)
              if(mod(idirface(k,i),2)==0) then
                if(idirface(k,i)==orienttgbay(itg)) then
                  tgflux=max(0.0,-qtidegate(im))
                elseif(idirface(k,i)==orienttgsea(itg)) then
                  tgflux=max(0.0,qtidegate(im))
                endif
                su(nck)=su(nck)+signface(idirface(k,i))*tgflux*ds(k,i)*tgflux/h(nck)
                spu(nck)=spu(nck)-abs(qtidegate(im))*ds(k,i)
                spv(nck)=spv(nck)-acoef(llec2llec(k,i),nck)
              else
                if(idirface(k,i)==orienttgbay(itg)) then
                  tgflux=max(0.0,-qtidegate(im))
                elseif(idirface(k,i)==orienttgsea(itg)) then
                  tgflux=max(0.0,qtidegate(im))
                endif
                sv(nck)=sv(nck)+signface(idirface(k,i))*tgflux*ds(k,i)*tgflux/h(nck)
                spv(nck)=spv(nck)-abs(qtidegate(im))*ds(k,i)  
                spu(nck)=spu(nck)-acoef(llec2llec(k,i),nck)
              endif
              acoef(llec2llec(k,i),nck)=0.0
            enddo
          endif
        enddo
      elseif(methtidegate(itg)==2)then !Bottom friction formulation
        do im=ntidegate(itg-1)+1,ntidegate(itg)
          i=idtidegate(im)
          if(opentidegate(im)==1) then   !Gate open
            wslbay=0.0
            snodebay=0.0
            wslsea=0.0
            snodesea=0.0
            do k=1,ncface(i)
              nck=cell2cell(k,i)
              if(idirface(k,i)==orienttgbay(itg)) then
                wslbay=wslbay+p(nck)/grav
                snodebay=snodebay+1.0
              endif
              if(idirface(k,i)==orienttgsea(itg)) then
                wslsea=wslsea+p(nck)/grav
                snodesea=snodesea+1.0
              endif
            enddo
            wslbay=wslbay/snodebay
            wslsea=wslsea/snodesea
            if(wslbay>=wslsea)then
              val=grav*(coeftidegate(itg,1)*coeftglateral(im))**2/h(i)**0.33333*uv(i)*areap(i)*iwet(i)
              spu(i)=spu(i)-val
              spv(i)=spv(i)-val
              do k=1,ncface(i)
                nck=cell2cell(k,i)
                if(idirface(k,i)==orienttgsea(itg).or.idirface(k,i)==orienttgbay(itg))then
                  val=grav*(coeftidegate(itg,1)*coeftglateral(im))**2/h(nck)**0.33333  &
                                *uv(nck)*areap(nck)*iwet(nck)
                  spu(nck)=spu(nck)-val
                  spv(nck)=spv(nck)-val
                endif  
              enddo
            else
              val=grav*(coeftidegate(itg,2)*coeftglateral(im))**2/h(i)**0.33333*uv(i)*areap(i)*iwet(i)
              spu(i)=spu(i)-val
              spv(i)=spv(i)-val
              do k=1,ncface(i)
                nck=cell2cell(k,i)
                if(idirface(k,i)==orienttgsea(itg).or.idirface(k,i)==orienttgbay(itg))then
                  val=grav*(coeftidegate(itg,2)*coeftglateral(im))**2/h(nck)**0.33333  &
                               *uv(nck)*areap(nck)*iwet(nck)  
                  spu(nck)=spu(nck)-val
                  spv(nck)=spv(nck)-val
                endif  
              enddo
            endif
          endif
        enddo
      endif
    enddo

!--- Weirs -------------------------------
    do iwr=1,numweir
      if(methweir(iwr)==1) then
        do im=nweir(iwr-1)+1,nweir(iwr)
          i=idweir(im)
          do k=1,ncface(i)
            nck=cell2cell(k,i)
            if(mod(idirface(k,i),2)==0) then
              if(idirface(k,i)==orientweirbay(iwr)) then
                tgflux=max(0.0,-qweir(im))
              elseif(idirface(k,i)==orientweirsea(iwr)) then
                tgflux=max(0.0,qweir(im))
              endif
              su(nck)=su(nck)+signface(idirface(k,i))*tgflux*ds(k,i)*tgflux/h(nck)
              spu(nck)=spu(nck)-abs(qweir(im))*ds(k,i)
              spv(nck)=spv(nck)-acoef(llec2llec(k,i),nck)
            else
              if(idirface(k,i)==orientweirbay(iwr)) then
                tgflux=max(0.0,-qweir(im))
              elseif(idirface(k,i)==orientweirsea(iwr)) then
                tgflux=max(0.0,qweir(im))
              endif
              su(nck)=su(nck)+signface(idirface(k,i))*tgflux*ds(k,i)*tgflux/h(nck)
              spv(nck)=spv(nck)-abs(qweir(im))*ds(k,i)  
              spu(nck)=spu(nck)-acoef(llec2llec(k,i),nck)
            endif
            acoef(llec2llec(k,i),nck)=0.0
          enddo
        enddo
      elseif(methweir(iwr)==2)then !Bottom friction formulation
        do im=nweir(iwr-1)+1,nweir(iwr)
          i=idweir(im)
          wslbay=0.0
          snodebay=0.0
          wslsea=0.0
          snodesea=0.0
          do k=1,ncface(i)
            nck=cell2cell(k,i)
            if(idirface(k,i)==orientweirbay(iwr)) then
              wslbay=wslbay+p(nck)/grav
              snodebay=snodebay+1.0
            endif
            if(idirface(k,i)==orientweirsea(iwr)) then
              wslsea=wslsea+p(nck)/grav
              snodesea=snodesea+1.0
            endif
          enddo
          wslbay=wslbay/snodebay
          wslsea=wslsea/snodesea
          if(wslbay>=wslsea) then
            val=grav*(coefweir(iwr,1)*coefweirlateral(im))**2/h(i)**0.33333*uv(i)*areap(i)*iwet(i)  
            spu(i)=spu(i)-val
            spv(i)=spv(i)-val
            do k=1,ncface(i)
              nck=cell2cell(k,i)
              if(idirface(k,i)==orientweirbay(iwr).or.idirface(k,i)==orientweirsea(iwr))then
                 val=grav*(coefweir(iwr,1)*coefweirlateral(im))**2/h(nck)**0.33333  &
                                *uv(nck)*areap(nck)*iwet(nck) 
                 spu(nck)=spu(nck)-val
                 spv(nck)=spv(nck)-val
              endif   
            enddo
          else
            val=grav*(coefweir(iwr,2)*coefweirlateral(im))**2/h(i)**0.33333*uv(i)*areap(i)*iwet(i)  
            spu(i)=spu(i)-val
            spv(i)=spv(i)-val
            do k=1,ncface(i)
              nck=cell2cell(k,i)
              if(idirface(k,i)==orientweirbay(iwr).or.idirface(k,i)==orientweirsea(iwr))then
                val=grav*(coefweir(iwr,2)*coefweirlateral(im))**2/h(nck)**0.33333  &
                               *uv(nck)*areap(nck)*iwet(nck)  
                spu(nck)=spu(nck)-val
                spv(nck)=spv(nck)-val
              endif  
            enddo
          endif
        enddo
      endif
    enddo

!--- Culverts --------------------------------
    do icv=1,numculvert
      i =idculvert(icv,1)
      i1=idculvert(icv,2)
      if(qculvert(icv)<0.0) then
        if(abs(angleculvertbay(icv)-pi)<=pi) then  !Skip if the angle is beyond 0 and 2*pi.
          su(i)=su(i)+abs(qculvert(icv)*uvculvert(icv))*cos(angleculvertbay(icv))
          sv(i)=sv(i)+abs(qculvert(icv)*uvculvert(icv))*sin(angleculvertbay(icv))
          spu(i)=spu(i)-abs(qculvert(icv))
          spv(i)=spv(i)-abs(qculvert(icv))
        endif
      else
        if(abs(angleculvertsea(icv)-pi)<=pi) then
          su(i1)=su(i1)+abs(qculvert(icv)*uvculvert(icv))*cos(angleculvertsea(icv))
          sv(i1)=sv(i1)+abs(qculvert(icv)*uvculvert(icv))*sin(angleculvertsea(icv))
          spu(i1)=spu(i1)-abs(qculvert(icv))
          spv(i1)=spv(i1)-abs(qculvert(icv))
        endif
      endif
    enddo

!--- Rubble Mound Structure --------------
    do irm=1,nrubmoundcells
      do im=nrubmound(irm-1)+1,nrubmound(irm)
        i=idrubmound(im)
        if(ntsch==1)then
          fac1=areap(i)*h1(i)/dtime
          spuv0(i)=-fac1
          ssu0(i)=fac1*u1(i)
          ssv0(i)=fac1*v1(i)
        else
          val=areap(i)/dtime
          fac1=val*ctsch1*h1(i)
          fac2=val*ctsch2*h2(i)
          spuv0(i)=-fac1+fac2
          ssu0(i)=fac1*u1(i)-fac2*u2(i)
          ssv0(i)=fac1*v1(i)-fac2*v2(i)
        endif
        su(i)=ssu0(i)-dpx(i)*h(i)*areap(i)
        sv(i)=ssv0(i)-dpy(i)*h(i)*areap(i)
        val=(rubmounda(irm)+rubmoundb(irm)*uv(i))*grav*h(i)*areap(i)*iwet(i)
        spu(i)=spu(i)-val
        spv(i)=spv(i)-val
        do k=1,ncface(i)
          nck=cell2cell(k,i)
          if(nck>ncells) cycle
          isRMS=0
          do imm=nrubmound(irm-1)+1,nrubmound(irm) !Consider internal nodes (not RMS) adjacent to i
            ii=idrubmound(imm)
            if(nck==ii) isRMS=1
          enddo  
          if(isRMS==0)then
            val=0.5*(rubmounda(irm)+rubmoundb(irm)*uv(nck)) *grav*h(nck)*areap(nck)*iwet(nck)
            spu(nck)=spu(nck)-val
            spv(nck)=spv(nck)-val
          endif                               
        enddo       
      enddo
    enddo
    
    return
    end subroutine struct_uv

!**********************************************************************************    
    subroutine struct_pp
! Calculates and adds source/sink terms for the pressure correction equation
! written by Weiming Wu and modified by Alex Sanchez, USACE-CHL
!**********************************************************************************
    use geo_def, only: ncface,cell2cell,idirface,ds,areap
    use flow_def, only: iwet,su,sp,p,p1,p2,pp,grav
    use struct_def
    use comvarbl, only: ntsch,dtime,ntsch,ctsch,ctsch1,ctsch2
    use prec_def
    implicit none
    integer :: i,i1,itg,iwr,icv,k,nck,im,irm
    real(ikind) :: timetmp1,gravdtime
    
    !--- Tidal Gate -------------------------------------------------------------
    do itg=1,numtidegate
      if(methtidegate(itg)==1) then
        do im=ntidegate(itg-1)+1,ntidegate(itg)
          i=idtidegate(im)
          if(opentidegate(im)==1) then   !Gate open
            do k=1,ncface(i)
              nck=cell2cell(k,i)
              if(idirface(k,i)==orienttgbay(itg)) su(nck)=su(nck)-qtidegate(im)*ds(k,i)
              if(idirface(k,i)==orienttgsea(itg)) su(nck)=su(nck)+qtidegate(im)*ds(k,i)
              if(idirface(k,i)==orienttgbay(itg)) then
                if(qtidegate(im)<0.0) then
                  sp(nck)=sp(nck)-abs(dqtgdzdown(im))*ds(k,i)/grav  !Grav appears due to dq/dpp=(dg/dzb)/grav
                else
                  sp(nck)=sp(nck)-abs(dqtgdzup(im))*ds(k,i)/grav
                endif
              endif
              if(idirface(k,i)==orienttgsea(itg)) then
                if(qtidegate(im)>0.0) then
                  sp(nck)=sp(nck)-abs(dqtgdzdown(im))*ds(k,i)/grav
                else
                  sp(nck)=sp(nck)-abs(dqtgdzup(im))*ds(k,i)/grav
                endif
              endif
            enddo
          endif
        enddo
      endif
    enddo

!--- Weir -----------------------------------------
    do iwr=1,numweir
      if(methweir(iwr)==1) then
        do im=nweir(iwr-1)+1,nweir(iwr)
          i=idweir(im)
          do k=1,ncface(i)
            nck=cell2cell(k,i)
            if(idirface(k,i)==orientweirbay(iwr)) su(nck)=su(nck)-qweir(im)*ds(k,i)
            if(idirface(k,i)==orientweirsea(iwr)) su(nck)=su(nck)+qweir(im)*ds(k,i)
            if(idirface(k,i)==orientweirbay(iwr)) then
              if(qweir(im)<0.0) then
                sp(nck)=sp(nck)-abs(dqweirdzdown(im))*ds(k,i)/grav
              else
                sp(nck)=sp(nck)-abs(dqweirdzup(im))*ds(k,i)/grav
              endif
            endif
            if(idirface(k,i)==orientweirsea(iwr)) then
              if(qweir(im)>0.0) then
                sp(nck)=sp(nck)-abs(dqweirdzdown(im))*ds(k,i)/grav
              else
                sp(nck)=sp(nck)-abs(dqweirdzup(im))*ds(k,i)/grav
              endif
            endif
          enddo
        enddo
      endif
    enddo

!--- Culvert --------------------------------
    do icv=1,numculvert
      i =idculvert(icv,1)
      i1=idculvert(icv,2)
      su(i)=su(i)-qculvert(icv)
      su(i1)=su(i1)+qculvert(icv)
      sp(i)=sp(i)-dqcvdzup(icv)/grav
      sp(i1)=sp(i1)+dqcvdzdown(icv)/grav
      su(i)=su(i)-dqcvdzdown(icv)/grav*pp(i1)
      su(i1)=su(i1)+dqcvdzup(icv)/grav*pp(i)
    enddo
    
!--- Rubble Mound Structure --------------
    gravdtime=grav*dtime
    do irm=1,nrubmoundcells
      i=idrubmound(irm)
      timetmp1=iwet(i)*areap(i)/gravdtime*(1.0-rubmoundporo(irm))
      if(ntsch==1)then
        su(i)=su(i)+(p(i)-p1(i))*timetmp1
      else
        su(i)=su(i)+(ctsch*p(i)-ctsch1*p1(i)+ctsch2*p2(i))*timetmp1
      endif  
      sp(i)=sp(i)+timetmp1
    enddo 
    
    return
    end subroutine struct_pp
    
!*******************************************************************************************    
    subroutine struct_flux
! Flux boundary condition at structures
! written by Weiming Wu and modified by Alex Sanchez, USACE-CHL
!*******************************************************************************************
    use geo_def, only: ncface,cell2cell,idirface,zb
    use flow_def, only: p,iwet,hmin,grav,sqrtgrav
    use struct_def
    use const_def, only: small,pi
    use prec_def
    use comvarbl,only: ctime   !(hli,12/13/13)
    implicit none
    integer :: i,k,im,nck,itg,iwr,icv,ibay,isea,idrycheck,idrybaysea
    !Weir
    real(ikind) :: wslbay,wslsea,snodesea,snodebay,alphasubmerg,dalphadzdown,dalphadzup
    real(ikind) :: weirdepth,weirdepth1,weirdepthd,weirtemp0,weirtemp1
    real(ikind) :: expweirdepthd1
    !Tidal Gate
    real(ikind) :: tidegatetemp,tidegatetemp1
    !Culvert
    real(ikind) :: cvperimeter,cvarea,cvwidth,cvdepth,cvdepthrev,cvtheta,dcvwidthdz
    real(ikind) :: dcvareadz,cvhalfwidth,culverttemp,culverttemp1,cvoutletref,cvhydrradius
    real(ikind) :: qculvertcritcl,qculvertnormaltemp,qculvertnormaltemp1,qculvertnormal
    real(ikind) :: culvertheadloss
    real(ikind) :: qculvert0(numculvert),dqculvertupdown
    
!--- Tidal Gate -----------------------
    do itg=1,numtidegate
      if(methtidegate(itg)==1) then
        do im=ntidegate(itg-1)+1,ntidegate(itg)
          i=idtidegate(im)
          if(opentidegate(im)==1) then   !Gate open
            wslbay=0.0
            snodebay=0.0
            wslsea=0.0
            snodesea=0.0
            idrycheck=1
            do k=1,ncface(i)
              nck=cell2cell(k,i)
              if(idirface(k,i)==orienttgbay(itg)) then
                wslbay=wslbay+p(nck)/grav
                snodebay=snodebay+1.0
                if(iwet(nck)==0) idrycheck=0
              endif
              if(idirface(k,i)==orienttgsea(itg)) then
                wslsea=wslsea+p(nck)/grav
                snodesea=snodesea+1.0
                if(iwet(nck)==0) idrycheck=0
              endif
            enddo
            wslbay=wslbay/snodebay
            wslsea=wslsea/snodesea
            tidegatetemp=openhgttidegate(itg)*coeftglateral(im)*sqrt(2.0*grav)*idrycheck
            if(wslbay>=wslsea) then
              tidegatetemp=tidegatetemp*coeftidegate(itg,1)*min(1.0,max(0.0,wslbay-zb(i))/openhgttidegate(itg))  &  !Depth Limit
                          *min(1.0,max(0.0,wslbay-zb(i))/max(small,wslbay-elevtidegate(itg)))  !Zb limit
              tidegatetemp1=sqrt(wslbay-wslsea)
              qtidegate(im)=tidegatetemp*tidegatetemp1
              dqtgdzdown(im)=-tidegatetemp*0.5/(tidegatetemp1+small)
              dqtgdzup(im)=tidegatetemp*0.5/(tidegatetemp1+small)
            elseif(wslsea>wslbay) then
              tidegatetemp=tidegatetemp*coeftidegate(itg,2)*min(1.0,max(0.0,wslsea-zb(i))/openhgttidegate(itg))  &
                          *min(1.0,max(0.0,wslsea-zb(i))/max(small,wslsea-elevtidegate(itg)))
              tidegatetemp1=sqrt(wslsea-wslbay)
              qtidegate(im)=-tidegatetemp*tidegatetemp1
              dqtgdzdown(im)=tidegatetemp*0.5/(tidegatetemp1+small)
              dqtgdzup(im)=-tidegatetemp*0.5/(tidegatetemp1+small)
            endif
          else
            qtidegate(im)=0.0
            dqtgdzdown(im)=0.0
            dqtgdzup(im)=0.0
          endif
        enddo
      endif
    enddo
    
!!$OMP SECTION
!--- Weir --------------------
    do iwr=1,numweir
      if(methweir(iwr)==1) then
        do im=nweir(iwr-1)+1,nweir(iwr)
          i=idweir(im)
          wslbay=0.0
          snodebay=0.0
          wslsea=0.0
          snodesea=0.0
          idrycheck=1
          do k=1,ncface(i)
            nck=cell2cell(k,i)
            if(idirface(k,i)==orientweirbay(iwr)) then
              wslbay=wslbay+p(nck)/grav
              snodebay=snodebay+1.0
              if(iwet(nck)==0) idrycheck=0
            endif
            if(idirface(k,i)==orientweirsea(iwr)) then
              wslsea=wslsea+p(nck)/grav
              snodesea=snodesea+1.0
              if(iwet(nck)==0) idrycheck=0
            endif
          enddo
          wslbay=wslbay/snodebay
          wslsea=wslsea/snodesea
          if(wslbay>=wslsea) then
            weirdepth=max(0.0,wslbay-zb(i))
            weirdepth1=weirdepth+small
            weirdepthd=wslbay-wslsea
            expweirdepthd1=exp(-8.5*weirdepthd/weirdepth1)
            alphasubmerg=1.0-expweirdepthd1
            dalphadzdown=-expweirdepthd1*8.5/weirdepth1
            dalphadzup=expweirdepthd1*8.5*(weirdepth1-weirdepthd)/weirdepth1**2
            weirtemp0=coefweir(iwr,1)*coefweirlateral(im)*sqrtgrav*weirdepth**0.5*idrycheck
            weirtemp1=weirtemp0*weirdepth
            qweir(im)=alphasubmerg*weirtemp1
            dqweirdzdown(im)=dalphadzdown*weirtemp1
            dqweirdzup(im)=dalphadzup*weirtemp1+alphasubmerg*1.5*weirtemp0
          elseif(wslsea>wslbay) then
            weirdepth=max(0.0,wslsea-zb(i))
            weirdepth1=weirdepth+small
            weirdepthd=wslsea-wslbay
            expweirdepthd1=exp(-8.5*weirdepthd/weirdepth1)
            alphasubmerg=1.0-expweirdepthd1
            dalphadzdown=-expweirdepthd1*8.5/weirdepth1
            dalphadzup=expweirdepthd1*8.5*(weirdepth1-weirdepthd)/weirdepth1**2
            weirtemp0=coefweir(iwr,2)*coefweirlateral(im)*sqrtgrav*weirdepth**0.5*idrycheck
            weirtemp1=weirtemp0*weirdepth
            qweir(im)=-alphasubmerg*weirtemp1
            dqweirdzdown(im)=-dalphadzdown*weirtemp1
            dqweirdzup(im)=-dalphadzup*weirtemp1-alphasubmerg*1.5*weirtemp0
          endif
        enddo
      endif
    enddo
    
!!$OMP SECTION
!--- Culvert -----------------------
    do icv=1,numculvert
      qculvert0(icv)=qculvert(icv)
      ibay=idculvert(icv,1)    !Bay node -Upstream
      isea=idculvert(icv,2)    !Sea node -Downstream
      wslbay=p(ibay)/grav
      wslsea=p(isea)/grav
      if(wslbay>wslsea+small) then  
        if(wslbay<=max(culvertelevbay(icv),culvertelevsea(icv))+hmin) then
          qculvert(icv)=0.0
          dqcvdzdown(icv)=0.0
          dqcvdzup(icv)=0.0
          uvculvert(icv)=0.0
          cycle
        endif
        if(wslsea>culvertelevsea(icv)+culvertheight(icv)) then  !Exit and inlet submerged
          if(iculverttype(icv)==1) then        !Circular culvert
            cvarea=pi*culvertrad(icv)**2
            cvhydrradius=0.5*culvertrad(icv)
          elseif(iculverttype(icv)==2) then    !Box culvert
            cvarea=culvertwidth(icv)*culvertheight(icv)
            cvhydrradius=cvarea/2.0/(culvertwidth(icv)+culvertheight(icv))
          endif
          culvertheadloss=cvheadlossbayentr(icv)+cvheadlossseaexit(icv)  &
                         +culvertfrict(icv)*culvertlength(icv)/4.0/cvhydrradius
          culverttemp=cvarea*sqrt(2.0*grav/culvertheadloss)
          culverttemp1=sqrt(wslbay-wslsea)
          qculvert(icv)=culverttemp*culverttemp1
          dqculvertupdown=culverttemp*0.5/(culverttemp1+small)
          dqcvdzdown(icv)=-dqculvertupdown
          dqcvdzup(icv)=dqculvertupdown
        else   
          if(wslbay>culvertelevbay(icv)+1.5*culvertheight(icv)) then  !Exit unsubmerged but inlet submerged
            if(iculverttype(icv)==1) then        !Circular culvert
              cvarea=pi*culvertrad(icv)**2
              cvhydrradius=0.5*culvertrad(icv)
            elseif(iculverttype(icv)==2) then    !Box culvert
              cvarea=culvertwidth(icv)*culvertheight(icv)
              cvhydrradius=cvarea/2.0/(culvertwidth(icv)+culvertheight(icv))
            endif
            culvertheadloss=cvheadlossbayentr(icv)+1.0  &
                           +culvertfrict(icv)*culvertlength(icv)/4.0/cvhydrradius
            culverttemp=cvarea*sqrt(2.0*grav/culvertheadloss)
            cvoutletref=culvertelevsea(icv)+culvertheight(icv)
            culverttemp1=sqrt(wslbay-cvoutletref)
            qculvert(icv)=culverttemp*culverttemp1
            dqcvdzdown(icv)=0.0
            dqcvdzup(icv)=culverttemp*0.5/(culverttemp1+small)
          else                 !Exit and inlet unsubmerged 
            cvdepth=max(0.6666*(wslbay-culvertelevbay(icv)),(wslsea-culvertelevsea(icv)))
            if(iculverttype(icv)==1) then        !Circular culvert
              if(cvdepth<=culvertrad(icv)) then
                cvdepthrev=culvertrad(icv)-cvdepth
                cvtheta=acos(cvdepthrev/culvertrad(icv))   
                cvhalfwidth=sqrt(culvertrad(icv)**2-cvdepthrev**2)
                cvarea=cvtheta*culvertrad(icv)**2-cvdepthrev*cvhalfwidth
                cvwidth=2.0*cvhalfwidth+small
                dcvwidthdz=2.0*cvdepthrev/(cvhalfwidth+small)
              else  !if(cvdepth>culvertrad(icv)) then
                cvdepthrev=cvdepth-culvertrad(icv)
                cvtheta=pi-acos(cvdepthrev/culvertrad(icv))
                cvhalfwidth=sqrt(culvertrad(icv)**2-cvdepthrev**2)
                cvarea=cvtheta*culvertrad(icv)**2+cvdepthrev*cvhalfwidth
                cvwidth=2.0*cvhalfwidth+small
                dcvwidthdz=-2.0*cvdepthrev/(cvhalfwidth+small)
              endif 
              cvperimeter=2.0*cvtheta*culvertrad(icv)
            elseif(iculverttype(icv)==2) then        !Box culvert
              cvarea=culvertwidth(icv)*cvdepth
              cvwidth=culvertwidth(icv)
              dcvwidthdz=0.0
              cvperimeter=culvertwidth(icv)+2.0*cvdepth
            endif
            dcvareadz=cvwidth   !dA/dzb
            qculvertcritcl=cvarea*sqrt(grav*cvarea/cvwidth)   !Critical flow
            cvhydrradius=cvarea/(cvperimeter+small)
            culvertheadloss=cvheadlossbayentr(icv)+cvheadlossseaexit(icv)  &
                           +2.0*grav*culvertmann(icv)**2*culvertlength(icv)/cvhydrradius**1.3333
            qculvertnormaltemp=cvarea*sqrt(2.0*grav/culvertheadloss)
            qculvertnormaltemp1=sqrt(wslbay-wslsea)
            qculvertnormal=qculvertnormaltemp*qculvertnormaltemp1     !Normal flow
            if(wslsea>culvertelevsea(icv).and.qculvertnormal<=qculvertcritcl) then
              qculvert(icv)=qculvertnormal
              dqculvertupdown=qculvertnormaltemp*0.5/(qculvertnormaltemp1+small)
              dqcvdzdown(icv)=-dqculvertupdown
              dqcvdzup(icv)=dqculvertupdown
            else                                        !Critical flow
              qculvert(icv)=qculvertcritcl
              dqcvdzdown(icv)=0.0
              dqcvdzup(icv)=0.5*sqrt(grav*cvarea/cvwidth**3)*(3.0*cvwidth*dcvareadz-cvarea*dcvwidthdz)
            endif
          endif
        endif
        uvculvert(icv)=qculvert(icv)/(cvarea+small)
      elseif(wslsea>wslbay+small) then  
        if(iculvertflap(icv)==1) then  ! with flap gate
          qculvert(icv)=0.0
          dqcvdzdown(icv)=0.0
          dqcvdzup(icv)=0.0
          uvculvert(icv)=0.0
          cycle
        endif
        if(wslsea<=max(culvertelevbay(icv),culvertelevsea(icv))+hmin) then
          qculvert(icv)=0.0
          dqcvdzdown(icv)=0.0
          dqcvdzup(icv)=0.0
          uvculvert(icv)=0.0
          cycle
        endif
        if(wslbay>culvertelevbay(icv)+culvertheight(icv)) then  !Exit and inlet submerged
          if(iculverttype(icv)==1) then        !Circular culvert
             cvarea=pi*culvertrad(icv)**2
             cvhydrradius=0.5*culvertrad(icv)
          elseif(iculverttype(icv)==2) then    !Box culvert
             cvarea=culvertwidth(icv)*culvertheight(icv)
             cvhydrradius=cvarea/2.0/(culvertwidth(icv)+culvertheight(icv))
          endif
          culvertheadloss=cvheadlossseaentr(icv)+cvheadlossbayexit(icv)  &
                         +culvertfrict(icv)*culvertlength(icv)/4.0/cvhydrradius
          culverttemp=cvarea*sqrt(2.0*grav/culvertheadloss)
          culverttemp1=sqrt(wslsea-wslbay)
          qculvert(icv)=-culverttemp*culverttemp1
          dqculvertupdown=culverttemp*0.5/(culverttemp1+small)
          dqcvdzdown(icv)=-dqculvertupdown
          dqcvdzup(icv)=dqculvertupdown
        else   
          if(wslsea>culvertelevsea(icv)+1.5*culvertheight(icv)) then  !Exit unsubmerged but inlet submerged
            if(iculverttype(icv)==1) then        !Circular culvert
              cvarea=pi*culvertrad(icv)**2
              cvhydrradius=0.5*culvertrad(icv)
            elseif(iculverttype(icv)==2) then    !Box culvert
              cvarea=culvertwidth(icv)*culvertheight(icv)
              cvhydrradius=cvarea/2.0/(culvertwidth(icv)+culvertheight(icv))
            endif
            culvertheadloss=cvheadlossseaentr(icv)+1.0  &
                           +culvertfrict(icv)*culvertlength(icv)/4.0/cvhydrradius
            culverttemp=cvarea*sqrt(2.0*grav/culvertheadloss)
            cvoutletref=culvertelevbay(icv)+culvertheight(icv)
            culverttemp1=sqrt(wslsea-cvoutletref)
            qculvert(icv)=-culverttemp*culverttemp1
            dqcvdzdown(icv)=-culverttemp*0.5/(culverttemp1+small)
            dqcvdzup(icv)=0.0
          else                 !Exit and inlet unsubmerged 
            cvdepth=max(0.6666*(wslsea-culvertelevsea(icv)),(wslbay-culvertelevbay(icv)))
            if(iculverttype(icv)==1) then        !Circular culvert
              if(cvdepth<=culvertrad(icv)) then
                cvdepthrev=culvertrad(icv)-cvdepth
                cvtheta=acos(cvdepthrev/culvertrad(icv))   
                cvhalfwidth=sqrt(culvertrad(icv)**2-cvdepthrev**2)
                cvarea=cvtheta*culvertrad(icv)**2-cvdepthrev*cvhalfwidth
                cvwidth=2.0*cvhalfwidth+small
                dcvwidthdz=2.0*cvdepthrev/(cvhalfwidth+small)
              else  !if(cvdepth>culvertrad(icv)) then
                cvdepthrev=cvdepth-culvertrad(icv)
                cvtheta=pi-acos(cvdepthrev/culvertrad(icv))
                cvhalfwidth=sqrt(culvertrad(icv)**2-cvdepthrev**2)
                cvarea=cvtheta*culvertrad(icv)**2+cvdepthrev*cvhalfwidth
                cvwidth=2.0*cvhalfwidth+small
                dcvwidthdz=-2.0*cvdepthrev/(cvhalfwidth+small)
              endif 
              cvperimeter=2.0*cvtheta*culvertrad(icv)
            elseif(iculverttype(icv)==2) then        !Box culvert
              cvarea=culvertwidth(icv)*cvdepth
              cvwidth=culvertwidth(icv)
              dcvwidthdz=0.0
              cvperimeter=culvertwidth(icv)+2.0*cvdepth
            endif
            dcvareadz=cvwidth   !dA/dzb
            qculvertcritcl=cvarea*sqrt(grav*cvarea/cvwidth)   !Critical flow
            cvhydrradius=cvarea/(cvperimeter+small)
            culvertheadloss=cvheadlossseaentr(icv)+cvheadlossbayexit(icv)  &
                           +2.0*grav*culvertmann(icv)**2*culvertlength(icv)/cvhydrradius**1.3333
            qculvertnormaltemp=cvarea*sqrt(2.0*grav/culvertheadloss)
            qculvertnormaltemp1=sqrt(wslsea-wslbay)
            qculvertnormal=qculvertnormaltemp*qculvertnormaltemp1     !Normal flow
            if(wslbay>culvertelevbay(icv).and.qculvertnormal<=qculvertcritcl) then
              qculvert(icv)=-qculvertnormal
              dqculvertupdown=qculvertnormaltemp*0.5/(qculvertnormaltemp1+small)
              dqcvdzdown(icv)=-dqculvertupdown
              dqcvdzup(icv)=dqculvertupdown
            else                                        !Critical flow
              qculvert(icv)=-qculvertcritcl
              dqcvdzdown(icv)=-0.5*sqrt(grav*cvarea/cvwidth**3)*(3.0*cvwidth*dcvareadz-cvarea*dcvwidthdz)
              dqcvdzup(icv)=0.0
            endif
          endif
        endif
        uvculvert(icv)=qculvert(icv)/(cvarea+small)
      else   
        qculvert(icv)=0.0
        dqcvdzdown(icv)=0.0
        dqcvdzup(icv)=0.0
        uvculvert(icv)=0.0
      endif
      qculvert(icv)=0.5*(qculvert0(icv)+qculvert(icv))   !Under-relaxation
      idrybaysea=iwet(ibay)*iwet(isea)
      qculvert(icv)=qculvert(icv)*idrybaysea
      dqcvdzdown(icv)=dqcvdzdown(icv)*idrybaysea
      dqcvdzup(icv)=dqcvdzup(icv)*idrybaysea
    enddo
!!$OMP END SECTIONS
        
    return
    end subroutine struct_flux
    
!***********************************************************************
    subroutine fluxthrustructure
!   calculate flow discharge through tide gates and weirs
!   by Weiming Wu, NCCHE,Dec. 2011
!***********************************************************************
    use size_def
    use geo_def
    use flow_def
    use struct_def
    use comvarbl,only: ctime
    use diag_def
    use diag_lib
    use prec_def
    use tool_def, only: vstrlz       
    implicit none
    integer :: i,k,im,itg,iwr,icv,ierr, counter

    !Tidal Gate
    do itg=1,numtidegate
      if(methtidegate(itg)==1) then
        Qtottidegate(itg)=0.0
        do im=ntidegate(itg-1)+1,ntidegate(itg)
          i=idtidegate(im)
          if(opentidegate(im)==1) then   !Gate open
             do k=1,ncface(i)
               if(idirface(k,i)==orienttgbay(itg)) &
                 Qtottidegate(itg)=Qtottidegate(itg)+qtidegate(im)*ds(k,i)
             enddo
          endif
        enddo
      elseif(methtidegate(itg)==2) then
        Qtottidegate(itg)=0.0
        do im=ntidegate(itg-1)+1,ntidegate(itg)
          i=idtidegate(im)
          if(opentidegate(im)==1) then   !Gate open
            do k=1,ncface(i)
              if(idirface(k,i)==orienttgbay(itg)) then  !Positive from bay to sea
                if(orienttgbay(itg)<=2) then
                  Qtottidegate(itg)=Qtottidegate(itg)-signface(idirface(k,i))*flux(k,i)
                else
                  Qtottidegate(itg)=Qtottidegate(itg)+signface(idirface(k,i))*flux(k,i)
                endif
              endif
            enddo
          endif
        enddo
      endif
    enddo

    !Weir
    do iwr=1,numweir
      if(methweir(iwr)==1) then
        Qtotweir(iwr)=0.0
        do im=nweir(iwr-1)+1,nweir(iwr)
          i=idweir(im)
          do k=1,ncface(i)
            if(idirface(k,i)==orientweirbay(iwr)) &
               Qtotweir(iwr)=Qtotweir(iwr)+qweir(im)*ds(k,i)
          enddo
        enddo
      elseif(methweir(iwr)==2) then
        Qtotweir(iwr)=0.0
        do im=nweir(iwr-1)+1,nweir(iwr)
          i=idweir(im)
          do k=1,ncface(i)
            if(idirface(k,i)==orientweirbay(iwr)) then
              if(orientweirbay(iwr)<=2) then     !Positive from bay to sea
                Qtotweir(iwr)=Qtotweir(iwr)-signface(idirface(k,i))*flux(k,i)
              else
                Qtotweir(iwr)=Qtotweir(iwr)+signface(idirface(k,i))*flux(k,i)
              endif
            endif
          enddo
        enddo
      endif
    enddo

    !Structures Screen Output
461 format('   Culvert ',I0,' Flow, ',A,A)
462 format('   Weir ',I0,' Flow, ',A,A)
463 format('   Tide Gate ',I0,' Flow, ',A,A)
464 format('   Tide Gate ',I0,' CLOSED')
    
    if(numculvert > 0 .or. numweir > 0 .or. numtidegate > 0) then
      if (numtidegate > 0) counter = count(mtidegateopt .eq. 4)                  !If there are tidegates, find how many are UNCONTROLLED
      if (counter > 0 .and. numtidegate > counter) call diag_print_message('')   !Print a blank line, except if all tidegates are UNCONTROLLED
      do icv=1,numculvert
        write(msg,461,iostat=ierr) icv, trim(vstrlz(qculvert(icv),'(f0.5)')),' m^3/s'
        call diag_print_message(msg)
      enddo
      do iwr=1,numweir
        write(msg,462,iostat=ierr) iwr, trim(vstrlz(qtotweir(iwr),'(f0.5)')),' m^3/s'
        call diag_print_message(msg)
      enddo
      do itg=1,numtidegate
        if (mtidegateopt(itg) .ne. 4) then                                       !Test if this tidegate is UNCONTROLLED and if it is, then do not print anything out. 
          if (qtottidegate(itg) .ne. 0.0) then
            write(msg,463,iostat=ierr) icv, trim(vstrlz(qtottidegate(itg),'(f0.5)')),' m^3/s'
          else
            write(msg,464,iostat=ierr) icv
          endif
          call diag_print_message(msg)
        endif
      enddo
      call diag_print_message('')
    endif
      
    return
    end subroutine fluxthrustructure

!**************************************************    
    subroutine struct_pond(iwet,ipond)
!**************************************************
    use size_def, only: ncellsD
    use geo_def, only: ncface,cell2cell
    use struct_def
    implicit none
    !Input/Output
    integer,intent(in) :: iwet(ncellsD)
    integer,intent(inout) :: ipond(ncellsD)
    !Internal Variables
    integer :: i,k,im,nck,itg,iwr,icv,i1
    
    !Tidal Gate 
    do itg=1,numtidegate
      do im=ntidegate(itg-1)+1,ntidegate(itg)
        i=idtidegate(im)
        ipond(i)=iwet(i)
        do k=1,ncface(i)
          nck=cell2cell(k,i)
          if(iwet(nck)==1) ipond(nck)=iwet(nck)
        enddo
      enddo
    enddo

    !Weir
    do iwr=1,numweir
      do im=nweir(iwr-1)+1,nweir(iwr)
        i=idweir(im)
        ipond(i)=iwet(i)
        do k=1,ncface(i)
          nck=cell2cell(k,i)
          if(iwet(nck)==1) ipond(nck)=iwet(nck)
        enddo
      enddo
    enddo

    !Culvert
    do icv=1,numculvert
      i=idculvert(icv,1)
      i1=idculvert(icv,2)
      ipond(i)=iwet(i)
      ipond(i1)=iwet(i1)
    enddo
    
    return
    end subroutine struct_pond

!**************************************************    
    subroutine struct_velbnd
!**************************************************    
    use size_def, only: ncellsD,nmaxfaces
    use geo_def, only: ncface,cell2cell,llec2llec,&
         idirface
    use flow_def, only: u,v
    use struct_def
    implicit none
    integer :: i,k,nck,itg,im,iwr
    
!--- Tidal Gate ---------------------------------------------------
    do itg=1,numtidegate
      if(methtidegate(itg)==1)then
        do im=ntidegate(itg-1)+1,ntidegate(itg)
          i=idtidegate(im)
          if(opentidegate(im)==1)then    !Gate open
            do k=1,ncface(i)
              nck=cell2cell(k,i)
              if(mod(idirface(k,i),2)==0)then
                if(idirface(k,i)==orienttgbay(itg).or.idirface(k,i)==orienttgsea(itg)) &
                   u(nck)=u(nck)*2.0     !The structure is treated as dry node, cdflux=0
              else
                if(idirface(k,i)==orienttgbay(itg).or.idirface(k,i)==orienttgsea(itg)) &
                   v(nck)=v(nck)*2.0     !The structure is treated as dry node, cdflux=0
              endif  
            enddo
          endif
        enddo
      endif
    enddo
    
!--- Weir ----------------------------
    do iwr=1,numweir
      if(methweir(iwr)==1) then
        do im=nweir(iwr-1)+1,nweir(iwr)
          i=idweir(im)
          do k=1,ncface(i)
             nck=cell2cell(k,i)
             if(mod(idirface(k,i),2)==0) then
                if(idirface(k,i)==orientweirbay(iwr).or.idirface(k,i)==orientweirsea(iwr)) &
                   u(nck)=u(nck)*2.0     !The structure is treated as dry node, cdflux=0
             else
                if(idirface(k,i)==orientweirbay(iwr).or.idirface(k,i)==orientweirsea(iwr)) &
                   v(nck)=v(nck)*2.0     !The structure is treated as dry node, cdflux=0
             endif  
          enddo
        enddo
      endif
    enddo
    
    return
    end subroutine struct_velbnd
    
!**************************************************    
    subroutine struct_vis
!**************************************************    
    use size_def, only: ncellsD,nmaxfaces
    use geo_def, only: ncface,cell2cell,llec2llec
    use flow_def, only: visk
    use struct_def
    implicit none
    integer :: i,k,nck,jcn,irm
    
!--- Rubble Mound Structures ----------------------  
    do irm=1,nrubmoundcells
      i=idrubmound(irm)
      do k=1,ncface(i)
        visk(k,i) = 0.0
        nck=cell2cell(k,i)
        jcn=llec2llec(k,i)
        visk(jcn,nck)=0.0
      enddo
    enddo   

    return
    end subroutine struct_vis
    
!**************************************************    
    subroutine struct_dzb
!**************************************************
    use struct_def
    use sed_def, only: nsed,dzb,dzbk
    implicit none
    integer :: i,itg,im,iwr,irm
    
    !Tidal Gate
    do itg=1,numtidegate   !Wu, Feb 28, 2010  
      do im=ntidegate(itg-1)+1,ntidegate(itg)
        i=idtidegate(im)
        dzb(i)=0.0
        if(nsed>1) dzbk(i,:)=0.0
      enddo
    enddo
          
    !Weir
    do iwr=1,numweir   !Wu, Sept, 2011  
      do im=nweir(iwr-1)+1,nweir(iwr)
        i=idweir(im)
        dzb(i)=0.0
        if(nsed>1) dzbk(i,:)=0.0
      enddo
    enddo

    !Rubble Mound Structure, !Wu, Oct. 2011    
    do irm=1,nrubmoundcells
      i=idrubmound(irm)
      dzb(i)=dzb(i)/rubmoundporo(irm)
      if(nsed>1) dzbk(i,:)=dzbk(i,:)/rubmoundporo(irm)
    enddo
    
    return
    end subroutine struct_dzb
    
!***********************************************************************
    subroutine struct_Ctk(ks)
! Applies sediment transport boundaries
! Blocks off dry regions which are not solved
! written by Alex Sanchez, USACE-CHL;  Weiming Wu, NCCHE
! Strucures added by Weiming Wu, NCCHE
!***********************************************************************
    use geo_def, only: ncface,cell2cell,llec2llec,idirface,ds
    use flow_def, only: acoef,su,sp
    use struct_def
    use sed_def, only: Ctk,rsk
    use const_def, only: small
    use prec_def
    implicit none
    integer :: i,i1,k,ks,nck,im,itg,iwr,icv
    real(ikind) :: Ctkbay,Ctksea,rtkbay,rtksea,snodebay,snodesea,tgflux
    
    !--- Tide Gate ---------------------------------------------------------------
    do itg=1,numtidegate   !Wu, Feb 28, 2010  
      if(methtidegate(itg)==1)then
        do im=ntidegate(itg-1)+1,ntidegate(itg)
          i=idtidegate(im)
          if(opentidegate(im)==1)then    !Gate open
            Ctkbay=0.0
            rtkbay=0.0
            snodebay=0.0
            Ctksea=0.0
            rtksea=0.0
            snodesea=0.0
            do k=1,ncface(i)
              nck=cell2cell(k,i)
              if(idirface(k,i)==orienttgbay(itg))then
                Ctkbay=Ctkbay+Ctk(nck,ks)
                rtkbay=rtkbay+Ctk(nck,ks)*rsk(nck,ks)  
                snodebay=snodebay+1.0
              endif
              if(idirface(k,i)==orienttgsea(itg))then
                Ctksea=Ctksea+Ctk(nck,ks)
                rtksea=rtksea+Ctk(nck,ks)*rsk(nck,ks)
                snodesea=snodesea+1.0
              endif
            enddo
            rtkbay=rtkbay/(Ctkbay+small)
            rtksea=rtksea/(Ctksea+small)  !Note the sequence of these lines
            Ctkbay=Ctkbay/snodebay
            Ctksea=Ctksea/snodesea
            do k=1,ncface(i)
              nck=cell2cell(k,i)
              if(idirface(k,i)==orienttgbay(itg))then
                tgflux=max(0.0,-qtidegate(im))
                su(nck)=su(nck)+tgflux*ds(k,i)*Ctksea*rtksea    !Only suspended load passes
                sp(nck)=sp(nck)-(tgflux+max(0.0,qtidegate(im))*rtkbay)*ds(k,i)
                acoef(llec2llec(k,i),nck)=0.0
              elseif(idirface(k,i)==orienttgsea(itg))then
                tgflux=max(0.0,qtidegate(im))
                su(nck)=su(nck)+tgflux*ds(k,i)*Ctkbay*rtkbay
                sp(nck)=sp(nck)-(tgflux+max(0.0,-qtidegate(im))*rtksea)*ds(k,i)
                acoef(llec2llec(k,i),nck)=0.0
              endif
            enddo
          endif
        enddo
      endif
    enddo

!--- Weir --------------------------------------------------
    do iwr=1,numweir   !Wu, Sept. 15, 2011  
      if(methweir(iwr)==1)then
        do im=nweir(iwr-1)+1,nweir(iwr)
          i=idweir(im)
          Ctkbay=0.0
          rtkbay=0.0
          snodebay=0.0
          Ctksea=0.0
          rtksea=0.0
          snodesea=0.0
          do k=1,ncface(i)
            nck=cell2cell(k,i)
            if(idirface(k,i)==orientweirbay(iwr))then
              Ctkbay=Ctkbay+Ctk(nck,ks) 
              rtkbay=rtkbay+Ctk(nck,ks)*rsk(nck,ks) 
              snodebay=snodebay+1.0
            endif
            if(idirface(k,i)==orientweirsea(iwr))then
              Ctksea=Ctksea+Ctk(nck,ks)*rsk(nck,ks)
              rtksea=rtksea+Ctk(nck,ks)*rsk(nck,ks)
              snodesea=snodesea+1.0
            endif
          enddo
          rtkbay=rtkbay/(Ctkbay+small)
          rtksea=rtksea/(Ctksea+small)  !Note the sequence of these lines
          Ctkbay=Ctkbay/snodebay
          Ctksea=Ctksea/snodesea
          do k=1,ncface(i)
            nck=cell2cell(k,i)
            if(idirface(k,i)==orientweirbay(iwr))then
              tgflux=max(0.0,-qweir(im))
              su(nck)=su(nck)+tgflux*ds(k,i)*Ctksea*rtksea   !Only suspended load passes
              sp(nck)=sp(nck)-(tgflux+max(0.0,qweir(im))*rtkbay)*ds(k,i)
              acoef(llec2llec(k,i),nck)=0.0
            elseif(idirface(k,i)==orientweirsea(iwr))then
              tgflux=max(0.0,qweir(im))
              su(nck)=su(nck)+tgflux*ds(k,i)*Ctkbay*rtkbay
              sp(nck)=sp(nck)-(tgflux+max(0.0,-qweir(im))*rtksea)*ds(k,i)
              acoef(llec2llec(k,i),nck)=0.0
            endif
          enddo
        enddo
      endif
    enddo
      
!--- Culvert ------------------------------------------------
    do icv=1,numculvert
      i=idculvert(icv,1)
      i1=idculvert(icv,2)
      if(qculvert(icv)>0.0)then
        sp(i)=sp(i)-qculvert(icv)*rsk(i,ks)
        su(i1)=su(i1)+qculvert(icv)*Ctk(i,ks)*rsk(i,ks)
        sp(i1)=sp(i1)-qculvert(icv)
      else   
        su(i)=su(i)-qculvert(icv)*Ctk(i1,ks)*rsk(i1,ks)
        sp(i)=sp(i)+qculvert(icv)
        sp(i1)=sp(i1)+qculvert(icv)*rsk(i1,ks)
      endif
    enddo
    
    return
    end subroutine struct_Ctk
    
!***********************************************************************
    subroutine struct_hardzb
!***********************************************************************
    use struct_def
    use sed_def, only: hardzb
    implicit none
    integer :: i,itg,im
    
    !Tidal Gate
    do itg=1,numtidegate       !Wu, Mar 1, 2010
       do im=ntidegate(itg-1)+1,ntidegate(itg)
          i=idtidegate(im)
          hardzb(i)=max(hardzb(i),elevtidegate(itg))
       enddo
    enddo
    
    return
    end subroutine struct_hardzb

!*****************************************************************    
    subroutine struct_sal
! Calculates source and sink terms for salinity at structures
! written by Weiming Wu, NCCHE
! modified by Alex Sanchez, USACE-CHL
!*****************************************************************    
    use size_def
    use geo_def
    use flow_def, only: acoef,sp
    use struct_def
    use sal_def
    use prec_def
    implicit none
    !Internal Variables
    integer :: i,k,nck,itg,iwr,icv,im,i1
    real(ikind) :: salbay,salsea,snodebay,snodesea,tgflux
    
    !--- Tide Gates -------------------------------------
    do itg=1,numtidegate     
      if(methtidegate(itg)==1) then
        do im=ntidegate(itg-1)+1,ntidegate(itg)
          i=idtidegate(im)
          if(opentidegate(im)==1) then   !Gate open
            salbay=0.0
            snodebay=0.0
            salsea=0.0
            snodesea=0.0
            do k=1,ncface(i)
               nck=cell2cell(k,i)
               if(idirface(k,i)==orienttgbay(itg)) then
                  salbay=salbay+sal(nck)
                  snodebay=snodebay+1.0
               endif
               if(idirface(k,i)==orienttgsea(itg)) then
                  salsea=salsea+sal(nck)
                  snodesea=snodesea+1.0
               endif
            enddo
            salbay=salbay/snodebay
            salsea=salsea/snodesea
            do k=1,ncface(i)
              nck=cell2cell(k,i)
              if(idirface(k,i)==orienttgbay(itg)) then
                 tgflux=max(0.0,-qtidegate(im))
                 susal0(nck)=susal0(nck)+tgflux*ds(k,i)*salsea
                 sp(nck)=sp(nck)-abs(qtidegate(im))*ds(k,i)
                 acoef(llec2llec(k,i),nck)=0.0
              elseif(idirface(k,i)==orienttgsea(itg)) then
                 tgflux=max(0.0,qtidegate(im))
                 susal0(nck)=susal0(nck)+tgflux*ds(k,i)*salbay
                 sp(nck)=sp(nck)-abs(qtidegate(im))*ds(k,i)
                 acoef(llec2llec(k,i),nck)=0.0
              endif
            enddo
          endif
        enddo
      endif
    enddo
       
!--- Weir --------------------------------------
    do iwr=1,numweir     
      if(methweir(iwr)==1) then
        do im=nweir(iwr-1)+1,nweir(iwr)
          i=idweir(im)
          salbay=0.0
          snodebay=0.0
          salsea=0.0
          snodesea=0.0
          do k=1,ncface(i)
            nck=cell2cell(k,i)
            if(idirface(k,i)==orientweirbay(iwr)) then
              salbay=salbay+sal(nck)
              snodebay=snodebay+1.0
            endif
            if(idirface(k,i)==orientweirsea(iwr)) then
              salsea=salsea+sal(nck)
              snodesea=snodesea+1.0
            endif
          enddo
          salbay=salbay/snodebay
          salsea=salsea/snodesea
          do k=1,ncface(i)
            nck=cell2cell(k,i)
            if(idirface(k,i)==orientweirbay(iwr)) then
              tgflux=max(0.0,-qweir(im))
              susal0(nck)=susal0(nck)+tgflux*ds(k,i)*salsea
              sp(nck)=sp(nck)-abs(qweir(im))*ds(k,i)
              acoef(llec2llec(k,i),nck)=0.0
            elseif(idirface(k,i)==orientweirsea(iwr)) then
              tgflux=max(0.0,qweir(im))
              susal0(nck)=susal0(nck)+tgflux*ds(k,i)*salbay
              sp(nck)=sp(nck)-abs(qweir(im))*ds(k,i)
              acoef(llec2llec(k,i),nck)=0.0
            endif
          enddo
        enddo
      endif
    enddo
            
!--- Culvert ------------------------------------------
    do icv=1,numculvert
      i=idculvert(icv,1)
      i1=idculvert(icv,2)
      if(qculvert(icv)>0.0) then
        sp(i)=sp(i)-qculvert(icv)
        susal0(i1)=susal0(i1)+qculvert(icv)*sal(i)
        sp(i1)=sp(i1)-qculvert(icv)
      else   
        susal0(i)=susal0(i)-qculvert(icv)*sal(i1)
        sp(i)=sp(i)+qculvert(icv)
        sp(i1)=sp(i1)+qculvert(icv)
      endif
    enddo
    
    return
    end subroutine struct_sal
