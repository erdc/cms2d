!***********************************************************************
    subroutine ST_hardbottom_AD
!***********************************************************************
    use EXP_Global_def,    only: hardbottom
    USE EXP_transport_def, only: bed, qsy, qsx, adss 
    use sed_def, only: poros,hardbed,scalemorph,nhard,idhard
    use geo_def, only: zb     
      
    implicit none
    !local vars
    integer i,id
    real tdepth 

    !check for hardbottom conditions and set qs's and eros to zero if hardbottom is encountered
    if(hardbottom) then
      do i=1,nhard
        id = idhard(i)
        tdepth = -zb(id)- bed(id)*scalemorph*POROS  
        if(tdepth.gt.hardbed(i)) then
          qsx(id) = 0
          qsy(id) = 0
          ADSS(id)%eros = 0
        endif
      enddo
    endif    
    
    return
    end subroutine
