!***********************************************************************
    subroutine pbicgstab(nmaxiter,acoef,ap,su,phi)
!   this is bicgstab solver       
!***********************************************************************      
    use size_def, only: ncells,ncellsD,nmaxfaces
    use geo_def, only: cell2cell,ncface
    use flow_def, only: iwet
    use const_def, only: small
	use solv_def 
    use prec_def
    !Intput/Output
    integer,intent(in) :: nmaxiter
    real(ikind),intent(in) :: acoef(nmaxfaces,ncellsD),ap(ncells),su(ncellsD)
    real(ikind),intent(inout) :: phi(ncellsD)
    !Internal Variables
    integer :: ierr,i,iter      
    real(ikind) :: pk(ncells),r(ncells),uk(ncells),vk(ncells),r0(ncells),zzk(ncellsD)
    real(ikind) :: omg,bet,bet0,sum1,sum2,ukr0,res,res0,term,alf,gam,adbk

	res0=0
!$OMP PARALLEL DO PRIVATE(i,term) REDUCTION(+:res0)
    do i=1,ncells
      term=sum(acoef(1:ncface(i),i)*phi(cell2cell(1:ncface(i),i)))+su(i)-ap(i)*phi(i)
	  r(i)=term+small
      r0(i)=r(i)
      res0=res0+abs(r0(i))
      pk(i)=0.0
      uk(i)=0.0
      vk(i)=0.0
    enddo
!$OMP END PARALLEL DO
!-------
	call ilu_cr(acoef,ap)
!    call ilu0(ncells)
    call ilutp(ncells,ierr)
!------
	zzk=0.0
    alf=1.0; bet0=1.0; gam=1.0

    !..ITERATION LOOPS BEGIN
    do iter=1,nmaxiter
       !bet=sum(r(1:ncells)*r0(1:ncells))
       bet=adbk(ncells,r,r0)
	   omg=bet*gam/(alf*bet0+small)
	   bet0=bet

      !..CALCULATE pk
!$OMP PARALLEL DO PRIVATE(i)          
      do i=1,ncells
        pk(i)=r(i)+omg*(pk(i)-alf*uk(i))
      enddo
!$OMP END PARALLEL DO        

      call lusol(ncells, pk(1:ncells), zzk(1:ncells))

      ukr0=0.0
!$OMP PARALLEL DO PRIVATE(i) REDUCTION(+:ukr0)
      do i=1,ncells
	     uk(i)=ap(i)*zzk(i)-sum(acoef(1:ncface(i),i)*zzk(cell2cell(1:ncface(i),i)))
         ukr0=ukr0+uk(i)*r0(i)
      enddo
!$OMP END PARALLEL DO           
      gam=bet/ukr0

      !..UPDATE (phi) AND CALCULATE R 
!$OMP PARALLEL DO PRIVATE(i)     
      do i=1,ncells
        phi(i)=phi(i)+gam*zzk(i)
	    r(i)=r(i)-gam*uk(i)
      enddo
!$OMP END PARALLEL DO        

      call lusol(ncells, r(1:ncells), zzk(1:ncells))	    

      !..CALCULATE vk and ALPHA (alf)
      sum1=0.0; sum2=0.0
!$OMP PARALLEL DO PRIVATE(i) REDUCTION(+:sum1), REDUCTION(+:sum2) 
      do i=1,ncells
        vk(i)=ap(i)*zzk(i)-sum(acoef(1:ncface(i),i)*zzk(cell2cell(1:ncface(i),i))) 
		sum1=sum1+vk(i)*r(i)
		sum2=sum2+vk(i)*vk(i)  
      enddo	 
!$OMP END PARALLEL DO        
      alf=sum1/max(sum2,small)

      !..UPDATE VARIABLE (phi) AND RESIDUAL (r) VECTORS
	  res=0.0
!$OMP PARALLEL DO PRIVATE(i) REDUCTION(+:res)     
      do i=1,ncells
        phi(i)=phi(i)+alf*zzk(i)
	    r(i)=r(i)-alf*vk(i)
        res=res+abs(r(i))
      enddo
!$OMP END PARALLEL DO        
      if(res/res0.lt.0.001) return 
    enddo
    
    return
    endsubroutine pbicgstab