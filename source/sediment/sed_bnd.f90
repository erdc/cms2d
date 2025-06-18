!========================================================================
! CMS Sediment Transport Boundary Conditions
!
! Contains the following:
!     sedbnd_init - Initializes the sediment transport boundary conditions
!     sedbnd_eval - Evaluates the sediment transport boundary conditions
!     bndzb       - Applies the bed elevation boundary condition for the
!                   implicit solution scheme
!
! written by Alex Sanchez, USACE-CHL
!========================================================================

!***********************************************************************
    subroutine sedbnd_init
!***********************************************************************
#include "CMS_cpp.h"
       use size_def
       use comvarbl, only: mpfile, flowpath
       use geo_def, only: ncface, cell2cell
       use bnd_def, only: nHstr, H_str, nTHstr, TH_str, nMHstr, MH_str, nQstr, Q_str
       use const_def, only: READONLY
       use sed_def
       use diag_lib
#ifdef XMDF_IO
       use xmdf
#endif
       use prec_def
       implicit none
       integer :: i, ks, ierr, nscell, nstimes
       integer(XID) PID, GID
       real(4), allocatable :: ftemp(:)
       character(len=200) :: astr, apath, aname
       character(len=10) :: aext
       type(sed_driver), allocatable :: sed_temp(:)
       logical :: ok
       character(len=300) :: msg

#ifdef XMDF_IO
       !====== Sediment Fluxes ==============================
       if (nsedflux > 0) then
          if (len_trim(sedfluxfile) == 0) then
             sedfluxfile = mpfile
          end if
          call fileparts(sedfluxfile, apath, aname, aext)
          if (len_trim(apath) == 0) then
             sedfluxfile = trim(flowpath)//sedfluxfile
          end if
          inquire (file=sedfluxfile, exist=ok)
          if (.not. ok) then
             write (msg, *) 'Could not find sediment boundary flux file: ', sedfluxfile
             call diag_print_error(msg)
          end if

          call XF_OPEN_FILE(sedfluxfile, READONLY, PID, ierr)
          do i = 1, nsedflux
             call XF_OPEN_GROUP(PID, sedfluxpath(i), GID, ierr)
             call XF_GET_PROPERTY_NUMBER(GID, 'Sed_Times', nstimes, ierr) !Determine number of time steps
             if (ierr < 0) cycle
             allocate (sed_str(nsedflux))
             sed_str(nsedflux)%ntimes = nstimes
             allocate (sed_str(nsedflux)%val(nstimes, nsed))
             allocate (sed_str(nsedflux)%times(nstimes))
             if (allocated(ftemp)) deallocate (ftemp)
             allocate (ftemp(nstimes))
             call XF_READ_PROPERTY_FLOAT(GID, 'Sed_Times', nstimes, ftemp(1), ierr)
             sed_str(nsedflux)%times = ftemp
             call XF_READ_PROPERTY_FLOAT(GID, 'Sediment', nstimes, ftemp(1), ierr)
             if (ierr == 0) then
                sed_str(nsedflux)%ibctype = 1 !Single source rate
                sed_str(nsedflux)%val(:, 1) = ftemp
             else
                sed_str(nsedflux)%ibctype = 2 !Fractional source rates
                do ks = 1, nsed
                   write (astr, 852) ks
                   call XF_READ_PROPERTY_FLOAT(GID, trim(astr), nstimes, ftemp(1), ierr)
                   if (ierr < 0) then
                      call diag_print_error('Invalid input format for sediment boundary condition')
                   end if
                   sed_str(nsedflux)%val(:, ks) = ftemp
                end do
                deallocate (ftemp)
             end if
             call XF_GET_PROPERTY_NUMBER(GID, 'Cells', nscell, ierr) !Determine number of cells in string
             sed_str(nsedflux)%ncells = nscell
             allocate (sed_str(nsedflux)%cells(nscell))
             call XF_READ_PROPERTY_INT(GID, 'Cells', nscell, sed_str(i)%cells(1), ierr)
          end do
          call XF_CLOSE_GROUP(GID, ierr)
          call XF_CLOSE_FILE(PID, ierr)
       end if

852    format('Sediment_', I0)
!=== Check hydro cell strings for sediment flux information ====
!--- River BC --------------------------------------------------
       do i = 1, nQstr
          call XF_OPEN_FILE(Q_str(i)%bidfile, READONLY, PID, ierr)
          call XF_OPEN_GROUP(PID, Q_str(i)%bidpath, GID, ierr)
          call XF_GET_PROPERTY_NUMBER(GID, 'Sed_Times', nstimes, ierr)
          if (ierr < 0) cycle
          nsedflux = nsedflux + 1
          !Increase sed_str size
          if (allocated(sed_str)) then
             allocate (sed_temp(nsedflux - 1))
             sed_temp = sed_str
             deallocate (sed_str)
             allocate (sed_str(nsedflux))
             sed_str(1:nsedflux - 1) = sed_temp
             deallocate (sed_temp)
          else
             allocate (sed_str(nsedflux))
          end if
          sed_str(nsedflux)%ntimes = nstimes
          allocate (sed_str(nsedflux)%val(nstimes, nsed))
          allocate (sed_str(nsedflux)%times(nstimes))
          if (allocated(ftemp)) deallocate (ftemp)
          allocate (ftemp(nstimes))
          call XF_READ_PROPERTY_FLOAT(GID, 'Sed_Times', &
                                      nstimes, sed_str(nsedflux)%times(1), ierr)
          call XF_READ_PROPERTY_FLOAT(GID, 'Sediment', nstimes, ftemp(1), ierr)
          if (ierr == 0) then
             sed_str(nsedflux)%ibctype = 1 !Single transport rate
             sed_str(nsedflux)%val(:, 1) = ftemp
          else
             sed_str(nsedflux)%ibctype = 2 !Fractional transport rates
             do ks = 1, nsed
                write (astr, 852) ks
                call XF_READ_PROPERTY_FLOAT(GID, trim(astr), nstimes, ftemp(1), ierr)
                if (ierr < 0) then
                   call diag_print_error('Invalid input format for sediment boundary condition')
                else
                   sed_str(nsedflux)%val(:, ks) = ftemp
                end if
             end do
             deallocate (ftemp)
          end if
          call XF_CLOSE_GROUP(GID, ierr)
          call XF_CLOSE_FILE(PID, ierr)

          !Cell id's and faces
          nscell = Q_str(i)%ncells
          sed_str(nsedflux)%ncells = nscell
          allocate (sed_str(nsedflux)%cells(nscell))
          allocate (sed_str(nsedflux)%faces(nscell))
          sed_str(nsedflux)%cells = Q_str(i)%cells
          sed_str(nsedflux)%faces = Q_str(i)%faces
       end do

!--- Tidal/Harmonic Water Level BC -----------------------------
       do i = 1, nTHstr
          call XF_OPEN_FILE(TH_str(i)%bidfile, READONLY, PID, ierr)
          call XF_OPEN_GROUP(PID, TH_str(i)%bidpath, GID, ierr)
          call XF_GET_PROPERTY_NUMBER(GID, 'Sed_Times', nstimes, ierr)
          if (ierr < 0) cycle
          nsedflux = nsedflux + 1
          !Increase sed_str size
          allocate (sed_temp(nsedflux - 1))
          sed_temp = sed_str
          deallocate (sed_str)
          allocate (sed_str(nsedflux))
          sed_str(1:nsedflux - 1) = sed_temp
          deallocate (sed_temp)
          sed_str(nsedflux)%ntimes = nstimes
          allocate (sed_str(nsedflux)%val(nstimes, nsed))
          allocate (sed_str(nsedflux)%times(nstimes))
          if (allocated(ftemp)) deallocate (ftemp)
          allocate (ftemp(nstimes))
          call XF_READ_PROPERTY_FLOAT(GID, 'Sed_Times', nstimes, ftemp(1), ierr)
          sed_str(nsedflux)%times = ftemp
          call XF_READ_PROPERTY_FLOAT(GID, 'Sediment', nstimes, ftemp(1), ierr)
          if (ierr == 0) then
             sed_str(nsedflux)%ibctype = 1 !Single transport rate
             sed_str(nsedflux)%val(:, 1) = ftemp
          else
             sed_str(nsedflux)%ibctype = 2 !Fractional transport rates
             do ks = 1, nsed
                write (astr, 852) ks
                call XF_READ_PROPERTY_FLOAT(GID, trim(astr), nstimes, ftemp(1), ierr)
                if (ierr < 0) then
                   call diag_print_error('Invalid input format for sediment boundary condition')
                end if
                sed_str(nsedflux)%val(:, ks) = ftemp
             end do
             deallocate (ftemp)
          end if
          call XF_CLOSE_GROUP(GID, ierr)
          call XF_CLOSE_FILE(PID, ierr)

          !Cell id's and faces
          nscell = TH_str(i)%ncells
          sed_str(nsedflux)%ncells = nscell
          allocate (sed_str(nsedflux)%cells(nscell))
          allocate (sed_str(nsedflux)%faces(nscell))
          sed_str(nsedflux)%cells = TH_str(i)%cells
          sed_str(nsedflux)%faces = TH_str(i)%faces
       end do

!--- Single Water Level BC ------------------------------------------------------
       do i = 1, nHstr
          call XF_OPEN_FILE(H_str(i)%bidfile, READONLY, PID, ierr)
          call XF_OPEN_GROUP(PID, H_str(i)%bidpath, GID, ierr)
          call XF_GET_PROPERTY_NUMBER(GID, 'Sed_Times', nstimes, ierr) !Determine number of time steps
          if (ierr < 0) cycle
          nsedflux = nsedflux + 1
          allocate (sed_str(nsedflux))
          sed_str(nsedflux)%ntimes = nstimes
          allocate (sed_str(nsedflux)%val(nstimes, nsed))
          allocate (sed_str(nsedflux)%times(nstimes))
          if (allocated(ftemp)) deallocate (ftemp)
          allocate (ftemp(nstimes))
          call XF_READ_PROPERTY_FLOAT(GID, 'Sed_Times', nstimes, ftemp(1), ierr)
          sed_str(nsedflux)%times = ftemp
          call XF_READ_PROPERTY_FLOAT(GID, 'Sediment', nstimes, ftemp(1), ierr)
          if (ierr == 0) then
             sed_str(nsedflux)%ibctype = 1 !Single transport rate
             sed_str(nsedflux)%val(:, 1) = ftemp
          else
             sed_str(nsedflux)%ibctype = 2 !Fractional transport rates
             do ks = 1, nsed
                write (astr, 852) ks
                call XF_READ_PROPERTY_FLOAT(GID, trim(astr), nstimes, ftemp(1), ierr)
                if (ierr < 0) then
                   call diag_print_error('Invalid input format for sediment boundary condition')
                end if
                sed_str(nsedflux)%val(:, ks) = ftemp
             end do
             deallocate (ftemp)
          end if
          call XF_CLOSE_GROUP(GID, ierr)
          call XF_CLOSE_FILE(PID, ierr)

          !Cell id's and faces
          nscell = H_str(i)%ncells
          sed_str(nsedflux)%ncells = nscell
          allocate (sed_str(nsedflux)%cells(nscell))
          allocate (sed_str(nsedflux)%faces(nscell))
          sed_str(nsedflux)%cells = H_str(i)%cells
          sed_str(nsedflux)%faces = H_str(i)%faces
       end do

!--- Multiple Water Level BC -------------------------
       do i = 1, nMHstr
          call XF_OPEN_FILE(MH_str(i)%bidfile, READONLY, PID, ierr)
          call XF_OPEN_GROUP(PID, MH_str(i)%bidpath, GID, ierr)
          call XF_GET_PROPERTY_NUMBER(GID, 'Sed_Times', nstimes, ierr)
          if (ierr < 0) cycle
          nsedflux = nsedflux + 1
          !Increase sed_str size
          allocate (sed_temp(nsedflux - 1))
          sed_temp = sed_str
          deallocate (sed_str)
          allocate (sed_str(nsedflux))
          sed_str(1:nsedflux - 1) = sed_temp
          deallocate (sed_temp)
          sed_str(nsedflux)%ntimes = nstimes
          allocate (sed_str(nsedflux)%val(nstimes, nsed))
          allocate (sed_str(nsedflux)%times(nstimes))
          if (allocated(ftemp)) deallocate (ftemp)
          allocate (ftemp(nstimes))
          call XF_READ_PROPERTY_FLOAT(GID, 'Sed_Times', nstimes, ftemp(1), ierr)
          sed_str(nsedflux)%times = ftemp
          call XF_READ_PROPERTY_FLOAT(GID, 'Sediment', nstimes, ftemp(1), ierr)
          if (ierr == 0) then
             sed_str(nsedflux)%ibctype = 1 !Single transport rate
             sed_str(nsedflux)%val(:, 1) = ftemp
          else
             sed_str(nsedflux)%ibctype = 2 !Fractional transport rates
             do ks = 1, nsed
                write (astr, 852) ks
                call XF_READ_PROPERTY_FLOAT(GID, trim(astr), nstimes, ftemp(1), ierr)
                if (ierr < 0) then
                   call diag_print_error('Invalid input format for sediment boundary condition')
                end if
                sed_str(nsedflux)%val(:, ks) = ftemp
             end do
             deallocate (ftemp)
          end if
          call XF_CLOSE_GROUP(GID, ierr)
          call XF_CLOSE_FILE(PID, ierr)

          !Cell id's and faces
          nscell = MH_str(i)%ncells
          sed_str(nsedflux)%ncells = nscell
          allocate (sed_str(nsedflux)%cells(nscell))
          allocate (sed_str(nsedflux)%faces(nscell))
          sed_str(nsedflux)%cells = MH_str(i)%cells
          sed_str(nsedflux)%faces = MH_str(i)%faces
       end do

       !===== Sediment sources =========================
       if (nsedsource > 0) then
          if (len_trim(sedsourcefile) == 0) then
             sedsourcefile = mpfile
          end if
          call fileparts(sedsourcefile, apath, aname, aext)
          if (len_trim(apath) == 0) then
             sedsourcefile = trim(flowpath)//sedsourcefile
          end if
          inquire (file=sedsourcefile, exist=ok)
          if (.not. ok) then
             write (msg, *) 'Could not find sediment source file: ', sedsourcefile
             call diag_print_error(msg)
          end if
          call XF_OPEN_FILE(sedsourcefile, READONLY, PID, ierr)
          do i = 1, nsedsource
             call XF_OPEN_GROUP(PID, sedsourcepath(i), GID, ierr)
             call XF_GET_PROPERTY_NUMBER(GID, 'Sed_Times', nstimes, ierr) !Determine number of time steps
             if (ierr < 0) cycle
             allocate (sed_str(nsedsource))
             sed_str(nsedflux)%ntimes = nstimes
             allocate (sed_str(nsedsource)%val(nstimes, nsed))
             allocate (sed_str(nsedsource)%times(nstimes))
             if (allocated(ftemp)) deallocate (ftemp)
             allocate (ftemp(nstimes))
             call XF_READ_PROPERTY_FLOAT(GID, 'Sed_Times', nstimes, ftemp(1), ierr)
             sed_str(nsedsource)%times = ftemp
             call XF_READ_PROPERTY_FLOAT(GID, 'Sediment', nstimes, ftemp(1), ierr)
             if (ierr == 0) then
                sed_str(nsedsource)%ibctype = 3 !Single source rate
                sed_str(nsedsource)%val(:, 1) = ftemp
             else
                sed_str(nsedsource)%ibctype = 4 !Fractional source rates
                do ks = 1, nsed
                   write (astr, 852) ks
                   call XF_READ_PROPERTY_FLOAT(GID, trim(astr), nstimes, ftemp(1), ierr)
                   if (ierr < 0) then
                      call diag_print_error('Invalid input format for sediment source')
                   end if
                   sed_str(nsedsource)%val(:, ks) = ftemp
                end do
                deallocate (ftemp)
             end if
             call XF_GET_PROPERTY_NUMBER(GID, 'Cells', nscell, ierr) !Determine number of cells in string
             sed_str(nsedsource)%ncells = nscell
             allocate (sed_str(nsedsource)%cells(nscell))
             call XF_READ_PROPERTY_INT(GID, 'Cells', nscell, sed_str(i)%cells(1), ierr)
!        do j=1,nscell
!          sed_str(nsedflux)%cells(j) = H_str(i)%cells(j)
!        enddo
          end do
          call XF_CLOSE_GROUP(GID, ierr)
          call XF_CLOSE_FILE(PID, ierr)
       end if

       nsedbc = nsedflux + nsedsource
       allocate (sedbnd(nsedbc, nsed))

#endif

       return
    end subroutine sedbnd_init

!***********************************************************************
    subroutine sedbnd_eval
! Applies sediment transport boundaries
! Blocks off dry regions which are not solved
! written by Alex Sanchez, USACE-CHL;  Weiming Wu, NCCHE
!***********************************************************************
       use size_def
       use geo_def
       use flow_def
       use struct_def
       use bnd_def
       use comvarbl
       use diag_def
       use diag_lib
       use wave_flowgrid_def
       use sed_def
       use const_def, only: eps
       use prec_def
       implicit none
       integer :: i, ii, j, k, ks, nck, ntimes, inc, ised, ibnd
       real(ikind) :: fac, qstartot, qstarcell, qsedtot

!--- All forcing boundaries -----------------------------
       do ibnd = 1, nbndstr
          do j = 1, bnd_str(ibnd)%ncells
             i = bnd_str(ibnd)%cells(j)
             k = bnd_str(ibnd)%faces(j)
             nck = cell2cell(k, i)
             if (flux(k, i) < 0.0) then !Inflow
                if (isedinflowbc == 1) then
                   CtstarP(nck, :) = facQtotin*CtstarP(i, :) !Capacity times loading factor
                elseif (isedinflowbc == 2) then
                   CtstarP(nck, :) = Qtotin/(uv(i)*h(i))/rhosed  !Qtotin in kg/m/sec
                end if
!        else
!          pbk(nck,:,:) = pbk(i,:,:)
!          pbk1(nck,:) = pbk1(i,:)
             end if
          end do !j
       end do !ibnd

       !Sediment flux boundaries, will overide river boundary conditions
       do ised = 1, nsedflux
          !find out where we are in the time/value arrays
          inc = sed_str(ised)%inc
          ntimes = sed_str(ised)%ntimes
          do while ((timehrs + eps) >= sed_str(ised)%times(inc + 1) .and. inc < ntimes)
             inc = inc + 1
             sed_str(ised)%inc = inc
          end do
          fac = (timehrs - sed_str(ised)%times(inc))/ &
                (sed_str(ised)%times(inc + 1) - sed_str(ised)%times(inc))

          !Get sediment transport rate
          if (sed_str(ised)%ibctype == 1) then !Total sediment transport rate
             qsedtot = (1.0 - fac)*sed_str(ised)%val(inc, 1) + &
                       fac*sed_str(ised)%val(inc + 1, 1)  !kg/sec, total per cellstring
             !Compute average bed composition for cell string
             do ii = 1, sed_str(ised)%ncells
                i = sed_str(ised)%cells(ii)
                k = sed_str(ised)%faces(ii)
                nck = cell2cell(k, i)
                pbk(nck, :, 1) = 0.0
                do ks = 1, nsed
                   pbk(nck, ks, 1) = pbk(nck, ks, 1) + pbk1(i, ks)
                end do
             end do
             pbk(nck, :, 1) = pbk(nck, :, 1)/float(sed_str(ised)%ncells)
             !Fractional bed composition based on average bed composition at boundary
             qsedtot = ramp*facQtotin*qsedtot !Apply ramp and loading factor
             sedbnd(ised, :) = pbk(nck, :, 1)*qsedtot
          else      !Fractional sediment transport rate
             sedbnd(ised, :) = (1.0 - fac)*sed_str(ised)%val(inc, :) + &
                               fac*sed_str(ised)%val(inc + 1, :) !kg/sec, per fraction and cellstring
             sedbnd(ised, :) = ramp*facQtotin*sedbnd(ised, :) !Apply ramp and loading factor
             !Sum sediment transport rate
             qsedtot = sum(sedbnd(ised, :))
!!        !Compute bed composition at ghost cells
!!        if(qsedtot>1.e-5)then
!!          do ii=1,sed_str(ised)%ncells
!!            i = sed_str(ised)%cells(ii)
!!            k = sed_str(ised)%faces(ii)
!!            nck = cell2cell(k,i)
!!            pbk(nck,:,1) = sedbnd(ised,:)/qsedtot
!!            pbk(nck,:,1) = pbk(nck,:,1)/sum(pbk(nck,:,1)) !Normalize to make sure sum is 1
!!          enddo
!!        else
!!          do ii=1,sed_str(ised)%ncells
!!            i = sed_str(ised)%cells(ii)
!!            k = sed_str(ised)%faces(ii)
!!            nck = cell2cell(k,i)
!!            pbk(nck,:,1) = pbk(i,:,1)
!!          enddo
!!        endif
          end if

          write (msg, *) 'Specified Total Inflow Sediment Transport Rate: ', qsedtot, ' kg/sec'
          call diag_print_message(msg)
!      write(*,*) 'Fractional Inflow Sediment Transport Rates, mm, kg/sec'
!      do ks=1,nsed
!        write(*,*) diam(ks)*1000.0, sedbnd(ised,ks)
!      enddo !ks

          qsedtot = 0.0
          do ks = 1, nsed
             qstartot = 0.0
             do ii = 1, sed_str(ised)%ncells
                i = sed_str(ised)%cells(ii)
                k = sed_str(ised)%faces(ii)
                nck = cell2cell(k, i)
                if (iwet(i) == 0) cycle
                if (flux(k, i) < 0.0) then !Inflow
!            qstarcell = h(i)*uv(i)*CtstarP(i,ks) !kg/m/sec
                   qstarcell = 1.0
                   qstarcell = max(qstarcell, 0.0001)
                   qstartot = qstartot + ds(k, i)*qstarcell !kg/sec
                end if
             end do !ii
             qstartot = max(qstartot, 1.0e-15)
             do ii = 1, sed_str(ised)%ncells
                i = sed_str(ised)%cells(ii)
                k = sed_str(ised)%faces(ii)
                nck = cell2cell(k, i)
                if (iwet(i) == 0) cycle
                if (flux(k, i) < 0.0) then
!            qstarcell = h(i)*uv(i)*CtstarP(i,ks) !kg/m/sec
                   qstarcell = 1.0
                   qstarcell = max(qstarcell, 0.0001)
                   fac = qstarcell/(h(i)*uv(i)*qstartot)
                   Ctkstar(nck, ks) = fac*sedbnd(ised, ks) !convert kg/sec to kg/m^3
                   !HLI 01/13/2017
                   CtstarP(nck, ks) = Ctkstar(nck, ks)/max(pbk(i, ks, 1), 1.0e-20) !Note i index in pbk, used in bound_c, pbk(i,ks,1)*CtstarP(nck,ks) used as boundary condition
                   qsedtot = qsedtot + ds(k, i)*h(i)*uv(i)*Ctkstar(nck, ks)  !kg/sec
                else
                   Ctkstar(nck, :) = Ctkstar(i, :)
                   CtstarP(nck, :) = CtstarP(i, :)
                end if
             end do !ii
          end do !ks
          write (msg, *) 'Calculated Total Inflow Sediment Transport Rate: ', qsedtot, ' kg/sec'
          call diag_print_message(msg)
       end do !ised

       return
    end subroutine sedbnd_eval

!***********************************************************************
    subroutine bndzb
! Update zb at all the boundaries
! written by Weiming Wu, NCCHE, Aug. 2009
! modified by Alex Sanchez, USACE-CHL
!***********************************************************************
       use size_def
       use geo_def, only: zb, cell2cell
       use bnd_def, only: nbndstr, bnd_str
       use sed_def, only: dzb
       use prec_def
       implicit none
       integer :: i, k, j, ibnd, nck

!--- All forcing boundaries ----------------
       do ibnd = 1, nbndstr
          do j = 1, bnd_str(ibnd)%ncells
             i = bnd_str(ibnd)%cells(j)
             k = bnd_str(ibnd)%faces(j)
             nck = cell2cell(k, i)
             dzb(nck) = dzb(i)
!        zb(nck)=zb(nck)+dzb(nck)
             zb(nck) = zb(i) !Alex, bug fix, June 8, 2010
!        pbk(nck,:,:) = pbk(i,:,:)
!        pbk1(nck,:) = pbk1(i,:)
          end do
       end do

       return
    end subroutine bndzb
