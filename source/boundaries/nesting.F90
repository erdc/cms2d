!====================================================================
!  Nesting CMS routines
!
! Contains:
!   nestwse_init - Initializes a nested WSE BC
!   nestwsevel_init - Initializes a nested WSE and velocity BC
!   tidalwse_init - Initializes a nested tidal database WSE BC
!   tidalwsevel_init - Initializes a nested tidal database WSE
!                      and velocity block
!   error_invalid_dataset - Error message for dataset read
!
!
! written by Alex Sanchez, USACE-CHL
!====================================================================

!*********************************************************************
    subroutine nestwse_init
! Initialized a Nested Water Level Boundary Condition
! written by Alex Sanchez, USACE-CHL
!*********************************************************************
       use bnd_def, only: NH_str, nNHstr, NH_type
       use geo_lib, only: proj_default
       implicit none
       integer :: i
       type(NH_type), allocatable :: NH_temp(:)

       nNHstr = nNHstr + 1
       if (nNHstr == 1) then
          allocate (NH_str(nNHstr))
       else
          allocate (NH_temp(nNHstr - 1))
          do i = 1, nNHstr - 1
             NH_temp(i) = NH_str(i)
          end do
          deallocate (NH_str)
          allocate (NH_str(nNHstr))
          do i = 1, nNHstr - 1
             NH_str(i) = NH_temp(i)
          end do
          deallocate (NH_temp)
       end if

       !Set defaults
       NH_str(nNHstr)%wseoffset = 0.0
       NH_str(nNHstr)%ncells = 0
       NH_str(nNHstr)%ntiwse = 1   !Interpolation order 1-linear, 2-quadratic
       NH_str(nNHstr)%wseout = .false.
       NH_str(nNHstr)%wsefile = ''

       return
    end subroutine nestwse_init

!*********************************************************************
    subroutine nestwsevel_init
! Initialized a Nested Water Level and Velocity Boundary Condition
! written by Alex Sanchez, USACE-CHL
!*********************************************************************
       use bnd_def, only: NHV_str, nNHVstr, NHV_type
       use geo_lib, only: proj_default
       implicit none
       integer :: i
       type(NHV_type), allocatable :: NHV_temp(:)

       nNHVstr = nNHVstr + 1
       if (nNHVstr == 1) then
          allocate (NHV_str(nNHVstr))
       else
          allocate (NHV_temp(nNHVstr - 1))
          do i = 1, nNHVstr - 1
             NHV_temp(i) = NHV_str(i)
          end do
          deallocate (NHV_str)
          allocate (NHV_str(nNHVstr))
          do i = 1, nNHVstr - 1
             NHV_str(i) = NHV_temp(i)
          end do
          deallocate (NHV_temp)
       end if

       !Set defaults
       NHV_str(nNHVstr)%wseoffset = 0.0
       NHV_str(nNHVstr)%ncells = 0
       NHV_str(nNHVstr)%ntiwse = 1   !Interpolation order 1-linear, 2-quadratic
       NHV_str(nNHVstr)%ntivel = 1   !Interpolation order 1-linear, 2-quadratic
       NHV_str(nNHVstr)%wseout = .false.
       NHV_str(nNHVstr)%velout = .false.
       NHV_str(nNHVstr)%wsefile = ''
       NHV_str(nNHVstr)%velfile = ''

       return
    end subroutine nestwsevel_init

!*********************************************************************
    subroutine parsim_init
! Initialized a Nested Water Level and Velocity Boundary Condition
! written by Alex Sanchez, USACE-CHL
!*********************************************************************
       use bnd_def, only: ParSim, nParSim, ParSim_type
       use geo_lib, only: proj_default
       implicit none
       integer :: i
       type(ParSim_type), allocatable :: ParSim_temp(:)

       nParSim = nParSim + 1
       if (nParSim == 1) then
          allocate (ParSim(nParSim))
       else
          allocate (ParSim_temp(nParSim - 1))
          do i = 1, nParSim - 1
             ParSim_temp(i) = ParSim(i)
          end do
          deallocate (ParSim)
          allocate (ParSim(nParSim))
          do i = 1, nParSim - 1
             ParSim(i) = ParSim_temp(i)
          end do
          deallocate (ParSim_temp)
       end if

       !Set defaults
       ParSim(nParSim)%timestarthr = -999.0 !Undefined
       ParSim(nParSim)%tjuldaypar = -999.0  !Undefined
       ParSim(nParSim)%ctlfilepar = ''
       ParSim(nParSim)%grdfilepar = ''
       ParSim(nParSim)%wsefilepar = ''
       ParSim(nParSim)%wsepathpar = ''
       ParSim(nParSim)%velfilepar = ''
       ParSim(nParSim)%velpathpar = ''
       ParSim(nParSim)%incwsepar = 0
       ParSim(nParSim)%incvelpar = 0
       ParSim(nParSim)%ntiwsepar = 1   !Interpolation order 1-linear, 2-quadratic
       ParSim(nParSim)%ntivelpar = 1   !Interpolation order 1-linear, 2-quadratic
       ParSim(nParSim)%velpar = .false. !Velocities on parent grid
       call proj_default(ParSim(nParSim)%projpar)

       return
    end subroutine parsim_init

!*********************************************************************
    subroutine tidalwse_init
! Initialized a nested tidal database wse boundary condition
! written by Alex Sanchez, USACE-CHL
!*********************************************************************
       use bnd_def, only: NTH_str, nNTHstr, NTH_type
       use geo_lib, only: proj_default
       implicit none
       integer :: i, ntc, nbndcells
       type(NTH_type), allocatable :: NTH_temp(:)

       nNTHstr = nNTHstr + 1
       if (nNTHstr == 1) then
          allocate (NTH_str(nNTHstr))
       else
          allocate (NTH_temp(nNTHstr - 1))
          do i = 1, nNTHstr - 1
             ntc = NTH_str(i)%ntc
             nbndcells = NTH_str(i)%ncells
             allocate (NTH_temp(i)%amp(nbndcells, ntc))
             allocate (NTH_temp(i)%phase(nbndcells, ntc))
             allocate (NTH_temp(i)%f(ntc))
             allocate (NTH_temp(i)%vu(ntc))
             allocate (NTH_temp(i)%speed(ntc))
             allocate (NTH_temp(i)%name(ntc))
             NTH_temp(i) = NTH_str(i)
          end do
          deallocate (NTH_str)
          allocate (NTH_str(nNTHstr))
          do i = 1, nNTHstr - 1
             ntc = NTH_temp(i)%ntc
             nbndcells = NTH_temp(i)%ncells
             allocate (NTH_str(i)%amp(nbndcells, ntc))
             allocate (NTH_str(i)%phase(nbndcells, ntc))
             allocate (NTH_str(i)%f(ntc))
             allocate (NTH_str(i)%vu(ntc))
             allocate (NTH_str(i)%speed(ntc))
             allocate (NTH_str(i)%name(ntc))
             NTH_str(i) = NTH_temp(i)
          end do
          deallocate (NTH_temp)
       end if

       !Set defaults
       NTH_str(nNTHstr)%wseadjust = .true.
       NTH_str(nNTHstr)%wseoffset = 0.0
       NTH_str(nNTHstr)%nssi = 0
       NTH_str(nNTHstr)%nssw = 0
       NTH_str(nNTHstr)%wseout = .false.
       call proj_default(NTH_str(nNTHstr)%projtdb)

       return
    end subroutine tidalwse_init

!*********************************************************************
    subroutine tidalwsevel_init
! Reads a nested tidal database wse and velocity block
! written by Alex Sanchez, USACE-CHL
!*********************************************************************
       use bnd_def, only: NTHV_str, nNTHVstr, NTHV_type
       use geo_lib, only: proj_default
       implicit none
       integer :: i, ntc, nbndcells
       type(NTHV_type), allocatable :: NTHV_temp(:)

       nNTHVstr = nNTHVstr + 1
       if (nNTHVstr == 1) then
          allocate (NTHV_str(nNTHVstr))
       else
          allocate (NTHV_temp(nNTHVstr - 1))
          do i = 1, nNTHVstr - 1
             ntc = NTHV_str(i)%ntc
             nbndcells = NTHV_str(i)%ncells
             allocate (NTHV_temp(i)%amp(nbndcells, ntc))
             allocate (NTHV_temp(i)%phase(nbndcells, ntc))
             allocate (NTHV_temp(i)%f(ntc))
             allocate (NTHV_temp(i)%vu(ntc))
             allocate (NTHV_temp(i)%speed(ntc))
             allocate (NTHV_temp(i)%name(ntc))
             NTHV_temp(i) = NTHV_str(i)
          end do
          deallocate (NTHV_str)
          allocate (NTHV_str(nNTHVstr))
          do i = 1, nNTHVstr - 1
             ntc = NTHV_temp(i)%ntc
             nbndcells = NTHV_temp(i)%ncells
             allocate (NTHV_str(i)%amp(nbndcells, ntc))
             allocate (NTHV_str(i)%phase(nbndcells, ntc))
             allocate (NTHV_str(i)%f(ntc))
             allocate (NTHV_str(i)%vu(ntc))
             allocate (NTHV_str(i)%speed(ntc))
             allocate (NTHV_str(i)%name(ntc))
             NTHV_str(i) = NTHV_temp(i)
          end do
          deallocate (NTHV_temp)
       end if

       !Set defaults
       NTHV_str(nNTHVstr)%wseoffset = 0.0
       NTHV_str(nNTHVstr)%nssi = 0
       NTHV_str(nNTHVstr)%nssw = 0
       NTHV_str(nNTHVstr)%wseout = .false.
       NTHV_str(nNTHVstr)%velout = .false.
       call proj_default(NTHV_str(nNTHVstr)%projtdb)

       return
    end subroutine tidalwsevel_init

!***********************************************************************
    subroutine error_invalid_dataset(afile, apath, adataset)
!***********************************************************************
       use diag_lib
       implicit none
       character(len=*) :: afile, apath, adataset

       call diag_print_error('Dataset not found: '//trim(adataset), '   File: '//trim(afile), '   Path: '//trim(apath))

       return
    end subroutine error_invalid_dataset

