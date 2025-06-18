!=======================================================================
!  Derivative CMS subroutines
!=======================================================================

!***********************************************************************
    subroutine der_default
! Sets the Spatial Derivative default values
! Author: Alex Sanchez, USACE-CHL
!***********************************************************************
       use der_def, only: nder, npow, nlim

       implicit none

       !Derivatives
       nder = 2 !Spatial derivative scheme
       npow = 1 !Inverse distance power, 1 seems to be more less sensitive to ill conditioning

       !Slope limiter
       nlim = 0            !Slope limiter

       return
    end subroutine der_default

!***********************************************************************
    subroutine der_cards(cardname, foundcard)
! Spatial Derivative cards
! Author: Alex Sanchez, USACE-CHL
!***********************************************************************
       use der_def, only: nder, npow, nlim, alim
       use diag_lib, only: diag_print_error

       implicit none
       integer :: i
       character(len=37) :: cardname, cdum
       logical :: foundcard

       foundcard = .true.
       select case (cardname)
       case ('SPATIAL_DERIVATIVE_METHOD', 'SPATIAL_DERIVATIVE_SCHEME', &
             'SPATIAL_DERIVATIVE', 'SPATIAL_DERIVATIVES')
          backspace (77)
          read (77, *) cardname, cdum
          select case (cdum)
          case ('CBGG', 'CELL-BASED_GREEN-GAUSS', 'CELL_BASED_GREEN_GAUSS')
             nder = 1  !Cell-based Green-Gauss (1st order for skewed cells)
          case ('CBGGCR', 'CELL-BASED_GREEN-GAUSS_CELL-RECON')
             nder = 2  !Cell-based Green-Gauss with cell-recontruction (second order)
!      case('CBGGFR','CELL-BASED_GREEN-GAUSS_FACE-RECON')
!        nder = 3  !Cell-based Green-Gauss with face-recontruction (second order)
!      case('CBGGLR','CELL-BASED_GREEN-GAUSS_LINEAR')
!        nder = 4  !Cell-based Green-Gauss with linear interpolation (second order)
!      case('NBGG','NODE-BASED_GREEN-GAUSS')
!        nder = 5  !Nodal-based Green-Gauss (has an extended stensil, first order)
          case ('CBWLSFS', 'CELL-BASED_WEIGHTED-LEAST-SQUARES')
             nder = 6  !Compact stenil cell-based least-squares (uses face-sharing neighboring cells) (first order)
!      case('CBWLSNS','EXTENDED_CELL-BASED_LEAST-SQUARES')
!        nder = 7  !Extended stenil cell-based least-squares (uses node-sharing neighboring cells) (first order, but more stable)
          case ('CBFD', 'CELL-BASED_FINITE-DIFFERENCE')
             nder = 8
          case default
             call diag_print_error('Invalid spatial derivative method or scheme')
          end select

       case ('INVERSE_DISTANCE_POWER')
          backspace (77)
          read (77, *) cardname, npow
          npow = min(max(npow, 0), 2)

       case ('SLOPE_LIMITER')
          backspace (77)
          read (77, *) cardname, cdum
          do i = 0, size(alim)
             if (cdum == alim(i)) then
                nlim = i
                exit
             end if
          end do

       case default
          foundcard = .false.
       end select

       return
    end subroutine der_cards

!***********************************************************************
    subroutine der_init
! Initializes derivatives module
! Author: Alex Sanchez, USACE-CHL
!***********************************************************************
       use size_def, only: ncells, ncellsd, ncelljoint, ncellsimple, ncellpoly, nmaxfaces, ndmaxfaces, nmaxextcells
       use flow_def, only: dHuxm, dHuym, dHvxm, dHvym, dapuareaxm, dapuareaym
       use comvarbl, only: skewcor
       use diag_lib, only: diag_print_error
       use geo_def, only: idcellsimple, idcelljoint, idirface, cell2cell, x, y, ncface
       use der_def, only: goa, gow, nder, nlistcw, ilistcw, nlistgw, ilistgw, der_go_type, nlim, rgx, rgy
       use der_lib, only: der_grad_cbfd, der_grad_cbwlsfs, der_alloc, der_grad_cbgg, der_grad_cbggcr

       implicit none
       integer :: i, ii, k, ncn, nce, ncs, ncw, nsc
       integer :: nlistca, nlistga
       integer :: ilistca(ncellsD), ilistga(ncellsD)

       !Boundary condition
       goa%ibc = 0 !No treatment
       gow%ibc = 1 !Only wet cells

       !Derivative scheme (use same scheme for now)
       if ((nder == 2 .and. ncelljoint == 0 .and. ncellpoly == 0) .or. &
           (nder == 2 .and. .not. skewcor)) then
          nder = 1 !No gradient corrections necessary for regular Cartesian grids
       end if
       goa%ider = nder
       gow%ider = nder

       !Size of operator
       goa%nd = ncellsD
       gow%nd = ncellsD

       !Maximum size of stencil for each direction
       if (ncellsimple == ncells) then !Simple Cartesian
          nsc = ndmaxfaces + 1
       elseif (ncelljoint > 0) then !Telescoping
          nsc = nmaxfaces + 1
       else !Polygonal
          nsc = nmaxfaces + nmaxextcells + 1
       end if
       goa%nsc = nsc
       gow%nsc = nsc

       !--- Allocate Gradient Operator -------
       call der_alloc

       !--- Calculation list ---------------------------
       nlistca = nlistcw
       ilistca = ilistcw
       nlistga = nlistgw
       ilistga = ilistgw

       !--- Initialize Gradient Operator -------------------------------
       select case (nder)
       case (1) !Cell-basd-Green-Gauss (first order)
          call der_grad_cbgg(goa, ncellsD, nlistca, ilistca) !Initialize all cells
          call der_grad_cbgg(gow, ncellsD, nlistcw, ilistcw) !Initialize wet cells
       case (2) !Cell-basd-Green-Gauss with cell-reconstruction
          call der_grad_cbgg(goa, ncellsD, nlistca, ilistca) !Initialize all cells
          call der_grad_cbgg(gow, ncellsD, nlistcw, ilistcw) !Initialize wet cells
          call der_grad_cbggcr(goa, ncellsD, nlistca, ilistca) !Initialize all cells
          call der_grad_cbggcr(gow, ncellsD, nlistcw, ilistcw) !Initialize wet cells
          !case(3) !Cell-basd-Green-Gauss with face-reconstruction
          !  call cbggfr_init
          !case(4) !Cell-basd-Green-Gauss with linear interpolation
          !  call cbgglin_init
          !case(6) !Node-based-Green-Gauss (second order)
          !  call nbgg_init
       case (6) !Weighted Least-squares method (compact stensil)
          call der_grad_cbwlsfs(goa, ncellsD, nlistca, ilistca) !Initialize all cells
          call der_grad_cbwlsfs(gow, ncellsD, nlistcw, ilistcw) !Initialize wet cells
          !case(7) !Least-squares method (extended stensil)
          !  call cbwlsns_init !CBWLSNS
       case (8) !Cell-based Finite-Difference (second order)
          call der_grad_cbfd(goa, ncellsD, nlistca, ilistca) !Initialize all cells
          call der_grad_cbfd(gow, ncellsD, nlistcw, ilistcw) !Initialize wet cells

       case default
          call diag_print_error('Choice must be CBGG, CBGGCR, CBWLSFS, or CBFS')

       end select

       call der_update

       !Slope limiters
       if (nlim == 0 .or. (ncellsimple + ncelljoint) == 0) return

       !Ratio of cell sizes
       do ii = 1, ncellsimple
          i = idcellsimple(ii)
          ncn = cell2cell(1, i)
          nce = cell2cell(2, i)
          ncs = cell2cell(3, i)
          ncw = cell2cell(4, i)
          rgx(i) = (x(i) - x(ncw))/(x(nce) - x(i))
          rgy(i) = (y(i) - y(ncs))/(y(ncn) - y(i))
       end do

       do ii = 1, ncelljoint
          i = idcelljoint(ii)
          do k = 1, ncface(i)
             select case (idirface(k, i))
             case (1); ncn = cell2cell(k, i)
             case (2); nce = cell2cell(k, i)
             case (3); ncs = cell2cell(k, i)
             case (4); ncw = cell2cell(k, i)
             end select
          end do
          rgx(i) = (x(i) - x(ncw))/(x(nce) - x(i))
          rgy(i) = (y(i) - y(ncs))/(y(ncn) - y(i))
       end do

       return
    end subroutine der_init

!***********************************************************************
    subroutine der_update()
! Determines what cells need to be reevaluated based on the wetting and
! drying state
!***********************************************************************
#include "CMS_cpp.h"
       use size_def, only: ncells, ncellsD, ncelljoint, ncellpoly
       use geo_def, only: idcelljoint
       use flow_def, only: iwet, iwet1
       use geo_def, only: cell2cell, ncface
       use der_def, only: nlistcw, ilistcw, nlistgw, ilistgw, nder, gow
       use der_lib, only: der_grad_cbgg, der_grad_cbggcr, der_grad_cbwlsfs, der_grad_cbfd

       implicit none
       !Internal variables
       integer :: i, ii, k
       integer :: ichg(ncellsD)

       !Find changing states and flag stencils
       ichg = 0 !Initialize
       do i = 1, ncells
          if (iwet(i) /= iwet1(i)) then
             ichg(i) = 1
             do k = 1, ncface(i)
                ichg(cell2cell(k, i)) = 1
             end do
          end if
       end do

       !Cell-centroid list
       nlistcw = 0
       do i = 1, ncells
          !if(ichg(i)==1 .and. iwet(i)==1)then !if dry then no need to calculate
          if (ichg(i) == 1) then !if dry then no need to calculate
             nlistcw = nlistcw + 1
             ilistcw(nlistcw) = i
          end if
       end do

       !Gradient list (only for CBGGCR method)
       nlistgw = 0
       if (nder == 2) then
          !Joint Cells
          do ii = 1, ncelljoint
             i = idcelljoint(ii)
             !if(ichg(i)==1 .and. iwet(i)==1)then !if dry then no need to calculate
             if (ichg(i) == 1) then !if dry then no need to calculate
                nlistgw = nlistgw + 1
                ilistgw(nlistgw) = i
             end if
          end do

          !Polygonal cells
          if (ncellpoly > 0) then
             nlistgw = nlistcw
             ilistgw = ilistcw
          end if
          !do i=1,ncellpoly
          !  !if(ichg(i)==1 .and. iwet(i)==1)then !if dry then no need to calculate
          !  if(ichg(i)==1)then !if dry then no need to calculate
          !    nlistgw = nlistgw + 1
          !    ilistgw(nlistgw) = i
          !  endif
          !enddo
       end if

#ifdef DIAG_MODE
       write (*, *) 'nlistcw: ', nlistcw, ', nlistgw: ', nlistgw
#endif

       if (nlistcw == 0 .and. nlistgw == 0) return !Early exit if no changes are needed

       select case (nder)
       case (1) !Cell-basd-Green-Gauss (first order)
          call der_grad_cbgg(gow, ncellsD, nlistcw, ilistcw) !Output
       case (2) !Cell-basd-Green-Gauss with cell-recontruction
          call der_grad_cbgg(gow, ncellsD, nlistcw, ilistcw) !Output
          call der_grad_cbggcr(gow, ncellsD, nlistcw, ilistcw) !Output
          !case(3) !Cell-basd-Green-Gauss with face-recontruction
          !  call cbggfr_init
          !case(4) !Cell-basd-Green-Gauss with linear interpolation
          !  call cbgglin_init
          !case(6) !Node-based-Green-Gauss (second order)
          !  call nbgg_init
       case (6) !Weighted Least-squares method (compact stensil)
          call der_grad_cbwlsfs(gow, ncellsD, nlistcw, ilistcw) !Output
          !case(7) !Least-squares method (extended stensil)
          !  call cbwlsns_init !CBWLSNS
       case (8) !Cell-based Finite-Difference (second order)
          call der_grad_cbfd(gow, ncellsD, nlistcw, ilistcw) !Output
       end select

       return
    end subroutine der_update
