!*********************************************************************************************
    subroutine par_wse_eval(wsefilepar, wsepathpar, nptspar, t2hrs, timehrst, nti, inc, timewsehrs, wsepar)
! Reads the parent solution file and performs spatial and temporal interpolations.
! The wse is first spatially interpolated and then interpolated in time.
! A description of the input/output variables is provided below
! written by Alex Sanchez, USACE-CHL
!*********************************************************************************************
#include "CMS_cpp.h"
       use comvarbl, only: timehrs
       use flow_def, only: iwet
       use in_lib, only: readstep63
       use nest_lib
       use diag_lib
       use prec_def
       implicit none
       !Input/Output
       character(len=*), intent(in)    :: wsefilepar          !Name and path of parent water level solution file
       character(len=*), intent(in)    :: wsepathpar          !Name and path of parent water level solution path
       integer, intent(in)    :: nptspar             !# of solution points of the parent grid
       real(ikind), intent(in)    :: t2hrs               !Conversion factor from the parent solution time to hours
       real(ikind), intent(in)    :: timehrst            !Parent solution starting time with respect to the child grid [hrs]
       integer, intent(in)    :: nti                 !Order of interpolation
       integer, intent(inout) :: inc                 !Time increment of the parent water level solution
       real(ikind), intent(inout) :: timewsehrs(nti + 1)   !Current velocitiy solution times used for temporal interpolation
       real(ikind), intent(inout) :: wsepar(nptspar, nti + 1) !Water levels at times 1 and 2 used for interpolation
       !Internal variables
       integer :: k, ierr, kth, nti1
       !real(ikind) :: wse(nptspar) !,fac
       character :: aext*10

       nti1 = nti + 1

       !Read dataset if necessary
       if (inc == 0) then
          call fileext(wsefilepar, aext)
          if (aext(1:2) == 'h5') then !XMDF file
#ifdef XMDF_IO
             do k = 1, nti1
                inc = inc + 1
                call readparscalsteph5(wsefilepar, wsepathpar, nptspar, inc, timewsehrs(k), wsepar(:, k), ierr)
                if (ierr < 0) call error_invalid_dataset(wsefilepar, wsepathpar, 'Water Elevation')
             end do
#else
             call diag_print_error('Cannot read parent grid water level from *.h5 file without XMDF libraries')
#endif
          elseif (aext(1:2) == '63') then !ADCIRC ASCII file
             open (63, file=wsefilepar, iostat=ierr)
             read (63, *)
             read (63, *)
             do k = 1, nti1
                inc = inc + 1
                call readstep63(nptspar, timewsehrs(k), wsepar(:, k), ierr)
             end do
          else
             call diag_print_error('Format not understood for file: ', wsefilepar)
          end if
          timewsehrs(:) = timewsehrs(:)*t2hrs + timehrst
       end if

       !Read new time(s) step
       kth = max(nti, 2) !Better than using the last value. Note nti is at least 1
       ierr = 0
       do while (timehrs >= timewsehrs(kth) .and. ierr >= 0)
          inc = inc + 1
          do k = 1, nti
             timewsehrs(k) = timewsehrs(k + 1)
             wsepar(:, k) = wsepar(:, k + 1)
          end do
          call fileext(wsefilepar, aext)
          if (aext(1:2) == 'h5') then !XMDF file
#ifdef XMDF_IO
             call readparscalsteph5(wsefilepar, wsepathpar, nptspar, inc, timewsehrs(nti1), wsepar(:, nti1), ierr)
             if (ierr < 0) call error_invalid_dataset(wsefilepar, wsepathpar, 'Water Elevation')
#else
             call diag_print_error('Cannot read parent grid water level from *.h5 file without XMDF libraries')
#endif
          else
             call readstep63(nptspar, timewsehrs(nti1), wsepar(:, nti1), ierr)
          end if
          timewsehrs(nti1) = timewsehrs(nti1)*t2hrs + timehrst
       end do

       return
    end subroutine par_wse_eval

!*******************************************************************************************************
    subroutine par_vel_eval(velfilepar, velpathpar, nptspar, t2hrs, timehrst, nti, inc, timevelhrs, upar, vpar)
! Reads the parent solution file and performs spatial and temporal interpolations.
! The wse is first spatially interpolated and then interpolated in time.
! A description of the input/output variables is provided below
! written by Alex Sanchez, USACE-CHL
!*******************************************************************************************************
#include "CMS_cpp.h"
       use comvarbl, only: timehrs
       use flow_def, only: iwet
       use in_lib, only: readstep64
       use nest_lib
       use diag_lib
       use prec_def
       implicit none
       !Input/Output
       character(len=*), intent(in)    :: velfilepar          !File name of parent current velocity solution file
       character(len=*), intent(in)    :: velpathpar          !File path of parent current velocity solution file
       integer, intent(in)    :: nptspar             !# of solution points of the parent grid
       real(ikind), intent(in)    :: t2hrs               !Conversion factor from the parent solution time to hours
       real(ikind), intent(in)    :: timehrst            !Parent solution starting time with respect to the child grid
       integer, intent(in)    :: nti                 !Order of interpolation
       integer, intent(inout) :: inc                 !Time increment of the parent current velocity solution
       real(ikind), intent(inout) :: timevelhrs(nti + 1)   !Current velocitiy solution times used for temporal interpolation
       real(ikind), intent(inout) :: upar(nptspar, nti + 1) !U-Velocity data on parent grid
       real(ikind), intent(inout) :: vpar(nptspar, nti + 1) !V-Velocity data on parent grid
       !Internal variables
       integer :: j, k, ierr, kth, nti1
       !real(ikind) :: ux(nptspar) !U-Velocity data on parent grid
       !real(ikind) :: vy(nptspar) !V-Velocity data on parent grid
       character :: aext*10

       nti1 = nti + 1

       !Read dataset if necessary
       if (inc == 0) then
          call fileext(velfilepar, aext)
          if (aext(1:2) == 'h5') then !XMDF file
#ifdef XMDF_IO
             do k = 1, nti1
                inc = inc + 1
                call readparvecsteph5(velfilepar, velpathpar, inc, nptspar, timevelhrs(k), upar(:, k), vpar(:, k), ierr)
                if (ierr < 0) call error_invalid_dataset(velfilepar, velpathpar, 'Current Velocity')
             end do
#else
             call diag_print_error('Cannot read parent grid velocity from *.h5 file without XMDF libraries')
#endif
          elseif (aext(1:2) == '64') then !ADCIRC ASCII file
             open (64, file=velfilepar)
             read (64, *)
             read (64, *)
             do k = 1, nti1
                inc = inc + 1
                call readstep64(nptspar, timevelhrs(k), upar(:, k), vpar(:, k), ierr)
             end do
             !elseif(aext(1:3)=='tsd')then  !TSD file (rows=[time_step],columns=[time,(u(i),v(i),i=1,nbndcells)])
             !call read_tsd(...)
          else
             call diag_print_error('Format not understood for file: ', velfilepar)
          end if
          timevelhrs(:) = timevelhrs(:)*t2hrs + timehrst
       end if

       !Read new time(s) step
       kth = max(nti, 2)
       ierr = 0
       do while (timehrs >= timevelhrs(kth) .and. ierr >= 0)
          inc = inc + 1
          do k = 1, nti
             timevelhrs(k) = timevelhrs(k + 1)
             do j = 1, nptspar
                upar(j, k) = upar(j, k + 1)
                vpar(j, k) = vpar(j, k + 1)
             end do
          end do
          call fileext(velfilepar, aext)
          if (aext(1:2) == 'h5') then !XMDF file
#ifdef XMDF_IO
             call readparvecsteph5(velfilepar, velpathpar, inc, nptspar, &
                                   timevelhrs(nti1), upar(:, nti1), vpar(:, nti1), ierr)
             if (ierr < 0) then
                call error_invalid_dataset(velfilepar, velpathpar, 'Current Velocity')
             end if
#else
             call diag_print_error('Cannot read parent grid velocity from *.h5 file without XMDF libraries')
#endif
          else
             call readstep64(nptspar, timevelhrs(nti1), upar(:, nti1), vpar(:, nti1), ierr)
          end if
          timevelhrs(nti1) = timevelhrs(nti1)*t2hrs + timehrst
       end do

       return
    end subroutine par_vel_eval
