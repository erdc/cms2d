      subroutine update_uv_from_qxqy()
	use EXP_Global_def
      USE EXP_bndcond_def
      USE EXP_transport_def 
      use sed_def
      use flow_def
      use comvarbl
      use size_def
      use geo_def, only: zb,cell2cell
      
      implicit none   
      integer i,id,istop
      real(ikind) waterdepth  
       
      istop=0
!$omp parallel do private(id, waterdepth)  !NLH 07/14/2008
      do i=1,ncells
        id = cell2cell(4,i)
        waterdepth = ( -zb(i)+eta(i) -zb(id)+eta(id))/2.0
        uE(i) = qx(i)/waterdepth
        id = cell2cell(3,i)
        waterdepth = ( -zb(i)+eta(i) -zb(id)+eta(id))/2.0
        vE(i) = qy(i)/waterdepth
      enddo
!$omp end parallel do

      do i=ncells+1,ncellsD
        id = cell2cell(4,i)
        if(id.gt.0) then
          waterdepth = ( -zb(i)+eta(i) -zb(id)+eta(id))/2.0
          uE(i) = qx(i)/waterdepth
        endif
        id = cell2cell(3,i)
        if(id.gt.0) then
          waterdepth = ( -zb(i)+eta(i) -zb(id)+eta(id))/2.0
          vE(i) = qy(i)/waterdepth
        endif
      enddo

      end subroutine