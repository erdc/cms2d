      subroutine update_uv_from_qxqy_Tel()
	use EXP_Global_def
      USE EXP_bndcond_def
      USE EXP_transport_def 
      use sed_def
      use flow_def
      use comvarbl
      use size_def
      use geo_def, only: zb,cell2cell
      use exp_telescoping
      
      implicit none   
      integer i,id,istop,id1,id2
      real(ikind) waterdepth  
       
      istop=0
      
      
!$omp parallel

!$omp do private(id1, id2, waterdepth)  !NLH 07/14/2008
      do i=1,numXfaces
        id1 = xface_cells(1,i)
        id2 = xface_cells(2,i)
        waterdepth = ( -zb(id1)+eta(id1) -zb(id2)+eta(id2))/2.0
        xface_vel(i) = xface_q(i)/waterdepth
      enddo
!$omp end do

!$omp do private(id1, id2, waterdepth)       
      do i=1,numYfaces
        id1 = yFace_cells(1,i)
        id2 = yface_cells(2,i)
        waterdepth = ( -zb(id1)+eta(id1) -zb(id2)+eta(id2))/2.0
        yface_vel(i) = yface_q(i)/waterdepth
      enddo
!$omp end do  

!$omp end parallel

      !do i=1,ncells
      !  id = cell2cell(4,i)
      !  waterdepth = ( -zb(i)+eta(i) -zb(id)+eta(id))/2.0
      !  uE(i) = qx(i)/waterdepth
      !  id = cell2cell(3,i)
      !  waterdepth = ( -zb(i)+eta(i) -zb(id)+eta(id))/2.0
      !  vE(i) = qy(i)/waterdepth
      !enddo


      !do i=ncells+1,ncellsD
      !  id = cell2cell(4,i)
      !  if(id.gt.0) then
      !    waterdepth = ( -zb(i)+eta(i) -zb(id)+eta(id))/2.0
      !    uE(i) = qx(i)/waterdepth
      !  endif
      !  id = cell2cell(3,i)
      !  if(id.gt.0) then
      !    waterdepth = ( -zb(i)+eta(i) -zb(id)+eta(id))/2.0
      !    vE(i) = qy(i)/waterdepth
      !  endif
      !enddo

      end subroutine