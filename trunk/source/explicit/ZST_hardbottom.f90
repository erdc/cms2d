      subroutine ST_hardbottom()
	  use EXP_Global_def 
      USE EXP_bndcond_def
      USE EXP_transport_def
      use flow_def
      use comvarbl
      use sed_def, only: poros,hardbed,scalemorph,nhard,idhard
      use geo_def, only: dx,dy,zb
        
      implicit none
      integer i,id
      real tdepth   

      !check for hardbottom conditions and set qs's to zero if hardbottom is encountered
      if(hardbottom) then
        do i=1,nhard
          id = idhard(i)
	      tdepth = -zb(id)- bed(id)*scalemorph*POROS          &
            + abs(qsx(id))*tsed_elapse*scalemorph*POROS/dx(i) &
            + abs(qsy(id))*tsed_elapse*scalemorph*POROS/dy(i)
          if(tdepth.gt.hardbed(i)) then
            qsx(id) = 0
            qsy(id) = 0
          endif
        enddo
      endif	
      
      end subroutine
