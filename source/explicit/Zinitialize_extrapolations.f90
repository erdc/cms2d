!********************************************************************************    
      subroutine initialize_extrapolations_CWR()
!********************************************************************************    
#include "CMS_cpp.h"    
      use EXP_Global_def,  only: num_ext_w, num_ext_e, num_ext_n, num_ext_s
      USE EXP_bndcond_def, only: ext_w, ext_e, ext_n, ext_s
      use bnd_def,  only: nhstr, nmhstr, nthstr, h_str, mh_str, th_str
      use size_def, only: ncells
      use geo_def, only: cell2cell

      implicit none
      integer i,j,id

      num_ext_W =0
      num_ext_E = 0
      num_ext_N = 0
      num_ext_S = 0

      if(nHstr .gt. 0 .or. nMHstr .gt. 0 .or. nTHstr .gt. 0) then
        !get count of each type of extrapolation
        if(nHstr .gt. 0) then
          do i=1,nHstr
            do j=1,H_Str(i)%ncells
              id = H_Str(i)%cells(j)
              if(cell2cell(1,id).gt.ncells) num_ext_N = num_ext_N + 1
              if(cell2cell(2,id).gt.ncells) num_ext_E = num_ext_E + 1
              if(cell2cell(3,id).gt.ncells) then
                if(cell2cell(1,id).gt.0) num_ext_S = num_ext_S + 1
              endif
              if(cell2cell(4,id).gt.ncells) then
                if(cell2cell(2,id).gt.0) num_ext_W = num_ext_W + 1
              endif
            enddo
          enddo
        endif
        if(nTHstr .gt. 0) then
          i=1
          do j=1,TH_str(i)%ncells
            id = TH_str(i)%cells(j)
            if(cell2cell(1,id).gt.ncells) num_ext_N = num_ext_N + 1
            if(cell2cell(2,id).gt.ncells) num_ext_E = num_ext_E + 1
            if(cell2cell(3,id).gt.ncells) then
              if(cell2cell(1,id).gt.0) num_ext_S = num_ext_S + 1
            endif
            if(cell2cell(4,id).gt.ncells) then
              if(cell2cell(2,id).gt.0) num_ext_W = num_ext_W + 1
            endif
          enddo
        endif    
        if(nMHstr .gt. 0) then
          do i=1,nMHstr
            do j=1,MH_Str(i)%ncells
              id = MH_Str(i)%cells(j)
              if(cell2cell(1,id).gt.ncells) num_ext_N = num_ext_N + 1
              if(cell2cell(2,id).gt.ncells) num_ext_E = num_ext_E + 1
              if(cell2cell(3,id).gt.ncells) then
                if(cell2cell(1,id).gt.0) num_ext_S = num_ext_S + 1
              endif
              if(cell2cell(4,id).gt.ncells) then
                if(cell2cell(2,id).gt.0) num_ext_W = num_ext_W + 1
              endif
            enddo
          enddo
        endif

        allocate(ext_w(num_ext_W,2),ext_n(num_ext_N,2),ext_s(num_ext_S,2),ext_e(num_ext_E,2))

        num_ext_W =0
        num_ext_E = 0
        num_ext_N = 0
        num_ext_S = 0

        if(nHstr .gt. 0) then
          do i=1,nHstr
            do j=1,H_Str(i)%ncells
              id = H_Str(i)%cells(j)
              if(cell2cell(1,id).gt.ncells) then
                num_ext_N = num_ext_N + 1
                ext_n(num_ext_N,1) = cell2cell(1,id)
                ext_N(num_ext_N,2) = id
              endif
              if(cell2cell(2,id).gt.ncells) then
                num_ext_E = num_ext_E + 1
                ext_E(num_ext_E,1) = cell2cell(2,id)
                ext_E(num_ext_E,2) = id
              endif
              if(cell2cell(3,id).gt.ncells) then
                if(cell2cell(1,id).gt.0) then
                  num_ext_S = num_ext_S + 1
                  ext_S(num_ext_S,1) = id
                  ext_S(num_ext_S,2) = cell2cell(1,id)
                endif
              endif
              if(cell2cell(4,id).gt.ncells) then
                if(cell2cell(2,id).gt.0) then
                  num_ext_W = num_ext_W + 1
                  ext_W(num_ext_W,1) = id
                  ext_W(num_ext_W,2) = cell2cell(2,id)
                endif
              endif
            enddo
          enddo
        endif

        if(nTHstr .gt. 0) then
          i=1
          do j=1,TH_Str(i)%ncells
            id = TH_Str(i)%cells(j)
            if(cell2cell(1,id).gt.ncells) then
              num_ext_N = num_ext_N + 1
              ext_n(num_ext_N,1) = cell2cell(1,id)
              ext_N(num_ext_N,2) = id
            endif
            if(cell2cell(2,id).gt.ncells) then
              num_ext_E = num_ext_E + 1
              ext_E(num_ext_E,1) = cell2cell(2,id)
              ext_E(num_ext_E,2) = id
            endif
            if(cell2cell(3,id).gt.ncells) then
              if(cell2cell(1,id).gt.0) then
                num_ext_S = num_ext_S + 1
                ext_S(num_ext_S,1) = id
                ext_S(num_ext_S,2) = cell2cell(1,id)
              endif
            endif
            if(cell2cell(4,id).gt.ncells) then
              if(cell2cell(2,id).gt.0) then
                num_ext_W = num_ext_W + 1
                ext_W(num_ext_W,1) = id
                ext_W(num_ext_W,2) = cell2cell(2,id)
              endif
            endif
          enddo
        endif

        if(nMHstr .gt. 0) then
          do i=1,nMHstr
            do j=1,MH_Str(i)%ncells
              id = MH_Str(i)%cells(j)
              if(cell2cell(1,id).gt.ncells) then
                num_ext_N = num_ext_N + 1
                ext_n(num_ext_N,1) = cell2cell(1,id)
                ext_N(num_ext_N,2) = id
              endif
              if(cell2cell(2,id).gt.ncells) then
                num_ext_E = num_ext_E + 1
                ext_E(num_ext_E,1) = cell2cell(2,id)
                ext_E(num_ext_E,2) = id
              endif
              if(cell2cell(3,id).gt.ncells) then
                if(cell2cell(1,id).gt.0) then
                  num_ext_S = num_ext_S + 1
                  ext_S(num_ext_S,1) = id
                  ext_S(num_ext_S,2) = cell2cell(1,id)
                endif
              endif
              if(cell2cell(4,id).gt.ncells) then
                if(cell2cell(2,id).gt.0) then
                  num_ext_W = num_ext_W + 1
                  ext_W(num_ext_W,1) = id
                  ext_W(num_ext_W,2) = cell2cell(2,id)
                endif
              endif
            enddo
          enddo
        endif
                
        
#ifdef DEBUG
        open(unit=505,file='check_extrapolations.csv')
        write(505,*)fac_UW,fac_DW
        if(num_ext_N.gt.0) then 
          write(505,*)num_ext_n      
          do i=1,num_ext_n
            write(505,*)'N ',i,ext_n(i,1),ext_N(i,2)
          enddo
        endif
        if(num_ext_S.gt.0) then 
          write(505,*)num_ext_S      
          do i=1,num_ext_S
            write(505,*)'S ',i,ext_S(i,1),ext_S(i,2)
          enddo
        endif    
        if(num_ext_E.gt.0) then 
          write(505,*)num_ext_E      
          do i=1,num_ext_E
            write(505,*)'E ',i,ext_E(i,1),ext_E(i,2)
          enddo
        endif      
        if(num_ext_W.gt.0) then 
          write(505,*)num_ext_W      
          do i=1,num_ext_W
            write(505,*)'W ',i,ext_W(i,1),ext_W(i,2)
          enddo
        endif    
        close(505)
        write(DGUNIT,*)'finished initializing bc extrapolation for u,v'
#endif
      endif !H_single or H_multi or H_tide
      
      return
      end subroutine
