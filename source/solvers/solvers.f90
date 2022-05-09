!***********************************************************************
    subroutine solv_default
!***********************************************************************
    use comvarbl, only: nsolv
    use solv_def
    implicit none

    nsolv = 4      !Solver, 1-Gauss-Seidel, 2-SOR, 3-BiCGSTAB  4-GMRES   
    asolv(0) = 'ADI'  !Still under testing
    asolv(1) = 'GAUSS-SEIDEL'
    asolv(2) = 'GAUSS-SEIDEL-SOR'
    asolv(3) = 'BICGSTAB'
    asolv(4) = 'GMRES'
    asolv(5) = 'HYBRID'
    asolv(6) = 'SIP'
    asolv(7) = 'ICCG'
    asolv(8) = 'ICCGSTAB' !Still under testing
    
    !Implicit Under-relaxation
    ppsolv%relax = 1.0   !SIMPLEC does not need under-relaxation
    velsolv%relax = 0.8
    salsolv%relax = 0.8
    heatsolv%relax = 0.8
    Ctksolv%relax = 0.8
    
    droptol = 1.0e-3
    permtol = 0.5
    lfil = 8
    
    iconv = 1
    
    return
    end subroutine solv_default
    
!*************************************************************   
    subroutine solv_block(kunit,xsolv)
! Reads wind data from Model Parameters file
! written by Alex Sanchez, USACE-CHL
!*************************************************************
    use solv_def, only: asolv,solver_options
    use prec_def
    implicit none
    !Input/Output
    integer,intent(in) :: kunit
    type(solver_options) :: xsolv
    !Internal Variables
    integer :: i,k,ierr
    character(len=37) :: cardname,cdum
    logical :: foundcard
    
d1: do k=1,10
      foundcard = .true.
      read(kunit,*,iostat=ierr) cardname
      if(ierr/=0) exit
      if(cardname(1:1)=='!' .or. cardname(1:1)=='#') cycle
      select case(cardname)    
      case('MATRIX_SOLVER','SOLVER_TYPE')
        backspace(77)
        read(77,*) cardname, cdum
        do i=0,size(asolv)
          if(cdum==asolv(i))then
            xsolv%isolv = i
            exit
          endif
        enddo         
        
      case('RELAXATION')
        backspace(77) 
        read(77,*) cardname, xsolv%relax
          
      case('SUCCESSIVE_OVER_RELAXATION_CONSTANT','SOR_CONSTANT','SOR') 
        backspace(77) 
        read(77,*) cardname, xsolv%relaxsor          
        
      case('ITERATIONS') 
        backspace(77) 
        read(77,*) cardname, xsolv%nswp
        xsolv%nswp0=xsolv%nswp
    
      case('SOLVER_END','END') 
        exit d1
      
      case default
          foundcard = .false.  
                          
      end select
    enddo d1
    
    return
    end subroutine solv_block
    
!***********************************************************************
    subroutine solv_init
!***********************************************************************    
    use size_def
    use comvarbl, only: nsolv
    use solv_def
    implicit none
    integer :: nt,nf

!    !Settings for solver 
!    if(ncelljoint>0 .and. (nsolv==6 .or. nsolv==7))then
!      nsolv = 4
!    endif
!    select case(nsolv)
!    case(0) !ADI
!      if(nswp0(1)==0) nswp(1)=40        
!      if(nswp0(2)==0) nswp(2)=10
!      if(nswp0(3)==0) nswp(3)=10
!      nswp(4)=10
!      nswp(5)=10    
!    case(1) !Gauss-Seidel
!      if(nswp0(1)==0) nswp(1)=100
!      if(nswp0(2)==0) nswp(2)=30
!      if(nswp0(3)==0) nswp(3)=30
!      nswp(4)=15
!      nswp(5)=15
!      if(maxit0==0) maxit=30+ncells/5000
!    case(2) !Gauss-Seidel-SOR
!      if(nswp0(1)==0) nswp(1)=80        
!      if(nswp0(2)==0) nswp(2)=20
!      if(nswp0(3)==0) nswp(3)=20
!      nswp(4)=15
!      nswp(5)=15        
!!      nswp(1)=100
!!      nswp(2)=15
!!      nswp(3)=15
!!      nswp(4)=15
!!      nswp(5)=10
!      if(maxit0==0) maxit=30+ncells/5000
!    case(3) !BiCGSTAB
!      if(nswp0(1)==0) nswp(1)=15
!      if(nswp0(2)==0) nswp(2)=5
!      if(nswp0(3)==0) nswp(3)=5
!      nswp(4)=5
!      nswp(5)=5
!      if(maxit0==0) maxit=20+ncells/10000
!    case(4) !GMRES
!      if(nswp0(1)==0) nswp(1)=20
!      if(nswp0(2)==0) nswp(2)=5
!      if(nswp0(3)==0) nswp(3)=5
!      nswp(4)=3
!!      nswp(4)=5
!      nswp(5)=3
!      if(maxit0==0) maxit=20+ncells/10000 
!    case(5) !Hybrid
!      if(nswp0(1)==0) nswp(1)=20
!      if(nswp0(2)==0) nswp(2)=5
!      if(nswp0(3)==0) nswp(3)=5
!      nswp(4)=3
!!      nswp(4)=5
!      nswp(5)=3
!     if(maxit0==0) maxit=20+ncells/10000  
!    case(6) !SIP
!      if(nswp0(1)==0) nswp(1)=40
!      if(nswp0(2)==0) nswp(2)=7
!      if(nswp0(3)==0) nswp(3)=7
!      nswp(4)=3
!      nswp(5)=3
!      if(maxit0==0) maxit=20+ncells/10000
!    case(7) !ICCG
!      if(nswp0(1)==0) nswp(1)=30
!      if(nswp0(2)==0) nswp(2)=7
!      if(nswp0(3)==0) nswp(3)=7
!      nswp(4)=3
!      nswp(5)=3
!      if(maxit0==0) maxit=20+ncells/10000
!    !case(8) !CGSTAB !Still under testing
!    !  if(nswp0(1)==0) nswp(1)=30
!    !  if(nswp0(2)==0) nswp(2)=7
!    !  if(nswp0(3)==0) nswp(3)=7
!    !  nswp(4)=3
!    !  nswp(5)=3
!    !  if(maxit0==0) maxit=20+ncells/10000      
!    end select      
!    nswp0=nswp      
!    maxit0=maxit      
!    rmom0=1.e10
!    rmom=0.0   
    
    if(nsolv<3 .or. nsolv>5) return
    
    nt = 40
    nf = nmaxfaces+1
    allocate(aa_matrix(nf*ncellsD))
    allocate(alu(nt*ncellsD))
    aa_matrix = 0.0; alu=0.0
    allocate(ia(ncells+1),ja(nt*ncellsD),iua(ncells))
    allocate(nposition1(nf*ncellsD,nmaxfaces),nposition2(nf*ncellsD,nmaxfaces))
    allocate(iaa(nf*ncellsD,0:nmaxfaces))
    allocate(nup(ncells),nlow(ncells),jlu(nt*ncellsD),ju(ncells)) 
    
    call indexm
    
    if(nsolv==4)then
       lfil=8
    else
       lfil=35
    endif
    
    return
    end subroutine solv_init    
    
!*******************************************************************************
    subroutine solve(acoef,ss,sp,res,phi,n)
!   call appropriate solver
!   variable index n: 1(pp), 2(u), 3(v)
!   by Weiming Wu, NCCHE, Oct. 2008
!*******************************************************************************
    use size_def, only: ncells,ncellsD,nmaxfaces
    use geo_def, only: ncface,cell2cell
    use flow_def, only: iwet
    use solv_def
    use comvarbl, only: relax,rmom,nsolv,nswp,nswp0,relaxsor
    use diag_def, only: debug_mode
    use diag_lib
    use prec_def
    implicit none
    !Input/Output
    integer,intent(in) :: n
    real(ikind),intent(in) :: acoef(nmaxfaces,ncellsD),sp(ncellsD)
    real(ikind),intent(inout) :: ss(ncellsD),phi(ncellsD)
    real(ikind),intent(out) :: res(ncellsD) !Normalized residuals
    !Internal Variables
    integer :: i
    real(ikind) :: ap(ncells),rinv,relaxm,adak
    logical :: isnankind

    rinv=1.0/relax(n)
    relaxm=1.0-relax(n)

!$OMP PARALLEL DO PRIVATE(i)
    do i=1,ncells
      ap(i)=(sum(acoef(1:ncface(i),i))-sp(i))*rinv  !Under-relaxation
      ss(i)=ss(i)+relaxm*phi(i)*ap(i)            !Under-relaxation
      res(i)=iwet(i)*(sum(acoef(1:ncface(i),i)*phi(cell2cell(1:ncface(i),i))) &
            +ss(i)-ap(i)*phi(i))/ap(i)
    enddo
!$OMP END PARALLEL DO
    rmom(n)=sqrt(adak(ncells,res)/real(ncells,kind=ikind))   !RMS of normalized residuals   
    
    select case(nsolv)
!!      case(0)
!!        call tdma2d(phi,ap,n)
      case(1)
        call gauss_seidel(nswp(n),acoef,ap,ss,phi)
      case(2)
        if(n==1 .and. rmom(1)<1.0e-3)then
          !call gauss_seidel_SOR(nswp(n),relaxsor,acoef,ap,ss,phi)
          call solv_ssor(nswp(n),relaxsor,acoef,ap,ss,phi)
        else
          !call gauss_seidel(nswp(n),acoef,ap,ss,phi)
          call solv_sgs(nswp(n),acoef,ap,ss,phi)
        endif
      case(3)
        call pbicgstab(nswp(n),acoef,ap,ss,phi)
      case(4)
        call pgmres(nswp(n),acoef,ap,ss,phi)
      case(5) !Hybrid pbicgstab-Gauss-Seidel 
        !if(niter<5 .and. n==1 .and. rmom(n)>1.0e-3)then
        if(rmom(1)>1.0e-4 .and. n==1)then
!!          call pgmres(nswp(n),acoef,ap,su,phi)
          call pbicgstab(nswp(n),acoef,ap,ss,phi)
        else 
          nswp(1)=nswp0(1)*4
          nswp(2:3)=nswp0(2:3)*4
          nswp(4:5)=nswp0(2:3)*3
          if(n==1 .and. rmom(n)<1.0e-3)then
            call gauss_seidel_SOR(nswp(n),relaxsor,acoef,ap,ss,phi)
          else
            call gauss_seidel(nswp(n),acoef,ap,ss,phi)
          endif 
          nswp(1:3)=nswp0(1:3)
        endif  
      case(6) !SIP
        call sip5(nswp(n),acoef,ap,ss,phi)
      case(7) !ICCG
        if(n==1)then
          call iccg5(nswp(n),acoef,ap,ss,phi) !Only for symmetric matrices
        else  
          call sip5(nswp(n),acoef,ap,ss,phi)
        endif  
      case(8) !ICCGSTAB !Still under testing
        if(n==1)then
          call iccg5(nswp(n),acoef,ap,ss,phi) !Only for symmetric matrices
          !call iccgstab5(nswp(n),acoef,ap,ss,phi) !Still under testing
        else  
          call iccgstab5(nswp(n),acoef,ap,ss,phi) !Still under testing
          !call sip5(nswp(n),acoef,ap,ss,phi)
        endif
    end select

    iconv = 1
    if(debug_mode)then
      do i=1,ncells
        if(isnankind(phi(i)))then         
          write(*,*) 'isnan(phi(i))'
          call diag_print_var(i,n)
          iconv = 0
          exit
        endif
      enddo
    endif
    
    return
    end subroutine solve
    