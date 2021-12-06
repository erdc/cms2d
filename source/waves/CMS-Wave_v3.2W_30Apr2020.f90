Subroutine CMS_Wave_inline !(noptset,nsteer)     !Wu
!
! Note: Presently used for Inline Steering, MEB 10/4/2018
!       Adding '_inline' to all subroutines (and references), MEB 10/9/2018
!

!***************************************************************C
!                      CMS-Wave (ver.3.2)                       C
!   (Wave Action Balance Equation with Diffraction effect)      C
!     SPECTRAL WAVE TRANSFORMATION MODEL IN SHALLOW WATER       C
!                INCLUDING WAVE DIFFRACTION                     C
!             in wave-current co-existing field                 C
!    PROGRAMED BY H.MASE and T.TAKAYAMA, JULY 1998, FEB 2001    C
!    FINISHED BY H.MASE AT UNIVERSITY OF LIVERPOOL, 30/12/'02   C
!    UPDATED BY LIHWA LIN, USACE ERDC, 30 April, 2020           C
!***************************************************************C
!---------------------------------------------------------------
!  Basic equation to be solved is wave action balance equation
!  with independent variables of x, y, q (angle).  
!  Intrinsic frequency is treated as a dependent variable.  
!  Energy dissipation term is formulated by using a modified
!  Miche's breaker index in order to include current effects.
!---------------------------------------------------------------
#include "CMS_cpp.h"
      use GLOBAL_INLINE
      use cms_def, only: noptset,nsteer,dtsteer,wavsimfile,wavepath
      use hot_def, only: coldstart
      use comvarbl, only: ctime
      use diag_lib, only: diag_print_message
      use diag_def, only: dgunit
!$include "omp_lib.h"
!      dimension dep0(ipmx,jpmx)
!      dimension fsp(npf),xc(mpd),yc(mpd),wc(mpd),wcd(mpd)
!      dimension refltx(ipmx,jpmx),reflty(ipmx,jpmx)
       double precision edate, jdate, iwind_date, icur_date, ieta_date
!
!      REAL, ALLOCATABLE :: dep0(:,:),dvarxx(:),dvaryy(:),   &   !Wu/Zhang
!      fsp(:),xc(:),yc(:),wc(:),wcd(:),hs13(:),              &   !Wu/Zhang
!      refltx(:,:),reflty(:,:)                                   !Wu/Zhang
      common /depreflx/dep0(ipmx,jpmx),fsp(npf),xc(komx),yc(komx),    &
                       wc(komx),wcd(komx),hs13(komx),                 &
                       refltx(ipmx,jpmx),reflty(ipmx,jpmx),&       !Wu
                       exx(igpx,jgpx),eyy(igpx,jgpx)
      common /dxxdyy/dvarxx(ipmx),dvaryy(jpmx)                   !Wu
      common /fl2wav/depin(ipmx,jpmx),etain(ipmx,jpmx), &        !Alex
             uin(ipmx,jpmx),vin(ipmx,jpmx)                       !Alex
      common /rsa/ sxx(ipmx,jpmx),sxy(ipmx,jpmx),syy(ipmx,jpmx)
      common /rsb/ wxrs(ipmx,jpmx),wyrs(ipmx,jpmx)
      common /rsc/ sxxx(ipmx,jpmx),sxyx(ipmx,jpmx)
      common /rsd/ sxyy(ipmx,jpmx),syyy(ipmx,jpmx)
      common /rsm/ cosa(mpd),sina(mpd),d1(ipmx,jpmx)
      common /rsw/ wk2(jpmx,npf,mpd),cgp(ipmx,jpmx)
      common /rsz/ igetfile20,nship,disx(ipmx),disy(jpmx)
      common /ship/ShipL(30),ShipB(30),ShipD(30),ShipS(30)
      common /struc0/ijstruc1,ijstruc2,ijstruc3,ijstruc4,ismall
      common /struc1/istruc1(komx),jstruc1(komx),dstruc1(komx)
      common /struc2/istruc2(NOMX),jstruc2(NOMX),dstruc2(NOMX)
      common /struc3/istruc3(komx),jstruc3(komx),dstruc3(komx) &
                    ,k3(komx),dstruc33(komx)
      common /struc4/istruc4(komx),jstruc4(komx),dstruc4(komx) &
                    ,kstruc4(komx),k4(komx),dstruc44(komx)
      common /swell/sw13(igpx,jgpx),sa13(igpx,jgpx) &
                   ,tw13(igpx,jgpx),ta13(igpx,jgpx) &
                   ,dw13(igpx,jgpx),da13(igpx,jgpx)
      COMMON /VPAI/PAI2,PAI,HPAI,RAD,akap,imod,iprp,island,imd,iprpp     &
                   ,nonln,igrav,isolv,ixmdf,iproc,imud,iwnd,depmin0
      COMMON /DATA/NF,MD,IMAX,JMAX,IGMX,JGMX,JCPB,JCPE,JCB,JCE,NFF,MDD
      COMMON /DATB/ICK3,ICK4,ICHOICE,HS0,wd,ws,ph0,wdd(mpd),aslop(ipmx)
      COMMON /DATC/TP,PRD,IBND,VL,TIDE,KPEAK,IBACK,NBLOCK,IWIND,WSMAG
      COMMON /DATD/DX,DY,DXX,DMESH,DTH,kdate,idate,depmax,depmax0
      COMMON /DATG/IBK,DBAR,WL0,WCC(JGPX),HSB(JGPX),HSG(JGPX),DCD(JGPX)
      COMMON /DVAR/DVARX(IGPX),DVARY(JGPX),ETA(IGPX,JGPX),HSK(JGPX)
      COMMON /IJBP/IJB(IGPX,JGPX),DEP(IPMX,JPMX),DEPS(IPMX),DBIG(IPMX)
      COMMON /REFA/KRMX,KR(2,6*IPMX),RK(6*IPMX),yangl(6*IPMX)
      COMMON /REFB/KRF,KCR(2,6*IPMX),RKR(6*IPMX),xangl(6*IPMX)
      COMMON /BREK/DEPM(JGPX),DMNJ(JGPX),SLF(JGPX),wlmn(JGPX),cmn(jgpx)  &
                   ,sigm(jgpx),IWVBK,ibk3
      COMMON /WAVI/H13(IGPX,JGPX),T13(IGPX,JGPX),DMN(IGPX,JGPX)
      COMMON /WAVR/H13R(IGPX,JGPX),T13R(IGPX,JGPX),DMNR(IGPX,JGPX)
      COMMON /WAVS/H13S(IGPX,JGPX),IBR(IGPX,JGPX),DISS(IGPX,JGPX)
      COMMON /WAVF/H13F(IGPX,JGPX),DMNF(IGPX,JGPX),T13F(IGPX,JGPX),T13MIN
      COMMON /OUTP/KOUT,IJSP(2,KOMX),IRS,IBREAK,ICUR,IWET,INST
      COMMON /OUTN/nest,inest(komx),jnest(komx),ix1(igpx),ix2(igpx)
      COMMON /FDS/FCN(NPF),DCM(MPD),DSFD(NPF,MPD)
      COMMON /DFCN/DF(NPF),FFCN(NPF),DFINP(NPF),PL1E(NPF),PL1T(NPF)
      common /uvp/u(ipmx,jpmx),v(ipmx,jpmx)
      common /uv1/u1(ipmx,jpmx),v1(ipmx,jpmx)
      common /uvwind/u10(ipmx,jpmx),v10(ipmx,jpmx)
      common /fric/bfric(ipmx,jpmx),amud(ipmx,jpmx)
      common /FileNames/ OptsFile, DepFile, CurrFile, EngInFile,    &
                         WaveFile, ObsFile, EngOutFile, NestFile,   &
                         BreakFile, RadsFile, StrucFile, SurgeFile, &
                         MudFile, FricFile, FrflFile, BrflFile,     &
                         SpecFile, WindFile, XMDFFile, SetupFile,   &  !Mitch 3/22/2017
                         SeaFile, SwellFile, ShipFile                  !Mitch 3/22/2017
      common /origin/x0,y0,azimuth,isteer,iidate,sinaz,cosaz
      common /comsav/depmin,g,iview,iplane,iwave,cflat,irunup,iroll
      common /wavegrid/ ni,nj                                          !Wu
      common /wavenum/ itms,ibf,iark,iarkr,bf,ark,arkr                 !Wu
      common /variables/hsmin,hs13n,inest1,nestin,nestin1,nestin2,azimnest      !Wu/Zhan/Alex
      common /logicalva/getfile4,getfile5,getfile7,getfile8            !Wu/Zhang 8/1/2009
!
      logical getfile,getfile1,getfile2,getfile3,getfile4,getfile5
      logical getfile6,getfile7,getfile8,getfile9,getfile10
      logical getfile11,getfile12,getfile15,getfile16,getfile17
      logical getfile18,getfile19,getfile20,getfile21
      CHARACTER*180 WaveFile,ObsFile,EngOutFile,BreakFile,RadsFile
      character*180 text
      character*30 text1
      REAL ITER_START, ITER_END
      INTEGER IOS  ! (IOS < 0) == END OF FILE; (IOS > 0) == IO ERROR
      INTEGER II,JJ,KK, INV
      integer numthreads,NTHR
!
! ... Input file variables
      CHARACTER*180  SimFile, OptsFile, DepFile, CurrFile, EngInFile
! ... Output/Input variable
      CHARACTER*180 NestFile, StrucFile, SurgeFile
      CHARACTER*180 MudFile,  FricFile,  FrflFile,  BrflFile, WindFile
      CHARACTER*180 SpecFile, XMDFFile,  SetupFile, ShipFile
      CHARACTER*180 SeaFile,  SwellFile                                 !Mitch 03/22/2017
      CHARACTER*80 :: cardname                                          !Mitch 10/18/2021
      
      logical :: foundfile, foundcard
! ... Output file variables
      INTEGER      :: iunit(2)     
!
      iunit = (/6,dgunit/)

      if(noptset.eq.3.and.nsteer.gt.1) goto 70                              !Wu/Zhang
  
!      ALLOCATE (fsp(npf),dep0(ipmx,jpmx),dvarxx(ipmx),dvaryy(jpmx))        !Wu/Zhang
!      ALLOCATE (xc(komx),yc(komx),wc(komx),wcd(komx),hs13(komx))           !Wu/Zhang
!      ALLOCATE (refltx(ipmx,jpmx),reflty(ipmx,jpmx))                       !Wu/Zhang

      rk =0.5
      rkr=0.3
      eta=0.0
      g=9.806
      PAI=3.141592654
      HPAI=PAI/2.0
      PAI2=2.0*PAI
      RAD=PAI/180.0
      kdate=0
      iidate=0
!
      if(noptset.ne.3)then

#ifdef MERGED_CODE
        call diag_print_message(' ')
        call diag_print_message('**********************************************')
        call diag_print_message('CMS-Wave V-3.2 Inline, last update 19 Oct 2021')
        call diag_print_message('**********************************************')
        call diag_print_message('  Point of Contact:')
        call diag_print_message('  Lihwa Lin, USACE ERDC')
        call diag_print_message('  mail to: Lihwa.Lin@usace.army.mil')
        call diag_print_message('**********************************************')
        call diag_print_message(' ')
#else
        print*,'CMS-Wave V-3.2 Inline, last update 19 Oct 2021'
        print*,'**********************************************'
        print*,'  Point of Contact:'
        print*,'  Lihwa Lin, USACE ERDC'
        print*,'  mail to: Lihwa.Lin@usace.army.mil'
        print*,'**********************************************'
        print*,' '
#endif
      endif

!Read filenames from .sim file and get correct paths set right.
      !call GetWaveFilenames_inline (SimFile)  !MEB 10/18/21 Removing this routine.  All it did was to make the next assignment to SimFile and then Call STWfiles_inline which is called here anyway.
      SimFile = trim(wavepath) // WavSimFile
      CALL STWfiles_inline (SimFile)

      inquire(file=OptsFile,exist=foundfile) 
      if(.not.foundfile)then
        write(*,*) 'ERROR: Could not find file: ',trim(OptsFile)
      endif
      open (11, file = OptsFile,  status = 'old')
      inquire(file=EngInFile,exist=foundfile) 
      if(.not.foundfile)then
        write(*,*) 'ERROR: Could not find file: ',trim(EngInFile)
      endif
      open (8,  file = EngInFile, status = 'old')
      inquire(file=DepFile,exist=foundfile) 
      if(.not.foundfile)then
        write(*,*) 'ERROR: Could not find file: ',trim(DepFile)
      endif
      open (15, file = DepFile,   status = 'old')
      open (23, file = StrucFile, status = 'unknown')
      read(23,'(a180)',end=339,err=339) text
      backspace(23)
      go to 340
339   continue
      close(23)
340   inquire(file='struct.dat',exist=getfile)
      if(getfile) then
        write(*,*) ' *** struct.dat FILE FOUND ***'
        write(*,*) '     Read struct.dat file'
        write(*,*) ' '
        if(noptset.eq.3)then
          write(9,*) ' *** struct.dat FILE FOUND ***'
          write(9,*) '     Read struct.dat file'
          write(9,*) ' '
        endif
        open(unit=23,file='struct.dat',status='old')
      end if
!
      iwnd=0
      iwind=0
      inquire(file='wind.dat',exist=getfile10)
      if(getfile10.and.iwnd.eq.0) then
        write(*,*) ' *** wind FILE FOUND ***'
        write(*,*) '      Read wind.dat'
        write(*,*) ' '
        if(noptset.eq.3)then
          write(9,*) ' *** wind FILE FOUND ***'
          write(9,*) '     Read wind file'
          write(9,*) ' '
        endif       
        open(unit=27,file='wind.dat',status='old')
        read(27,*,end=181,err=181) niwindc, njwindc
        iwind=1
        go to 180
      end if
181   continue
!
      inquire(file=WindFile,exist=getfile10)
      if(getfile10.and.iwnd.eq.0) then
        write(*,*) ' *** wind FILE FOUND ***'
        write(*,*) '     Read wind file'
        write(*,*) ' '
        if(noptset.eq.3)then
          write(9,*) ' *** wind FILE FOUND ***'
          write(9,*) '     Read wind file'
          write(9,*) ' '
        endif       
        open(unit=27,file= WindFile,status='old')
        read(27,*,end=180,err=180) niwindc, njwindc
        iwind=1
      end if
180   continue      
      if(iwind.eq.0) iwnd=1
!
      iwind1=0
      inquire(file='wind1.dat',exist=getfile16)
      if(getfile16.and.iwnd.eq.0) then
        open(unit=7,file='wind1.dat',status='old')
        read(7,*) iwind1
        close(7)
        do l=1,179
          if(SimFile(l:l+1).eq.'  ') exit
        end do
        if(l.lt.14) l=14
        if (SimFile(l-13:l-5).eq.'\1swsteer') then
          iwind1=0
          write(*,*) '*** Wind Steering No.1 Cycle ***'
        end if
      end if
!        
      if(getfile10.and.iwnd.eq.0) then
        open(unit=7,file='wind1.dat',status='unknown')
        write(7,*) iwind1+1
        write(*,*) ' '
        write(*,*) ' *** Please check Cycle Number ***'
        write(*,*) '     You are in Cycle',iwind1+1
        write(*,*) '     Check Cycle Number is correct not!!'
        write(*,*) ' *** Modify wind1.dat for correct cycle?'
        write(*,*) ' '
        close(7)
      end if 
!
!Set default values for parameters
      iprpp=1
      icur=0
      ibreak=0
      irs=0
      kout=0
      ibnd=0
      iwet=0
      ibf=0
      iark=0
      iarkr=0
      akap=1.0
!     bf=0.001
      bf=0.0
      ark=0.5
      arkr=0.3
      iwvbk=0
      iback=0
      nonln=0
      igrav=0
      irunup=0
      imud=0
      iwnd=0
      isolv=0
      ixmdf=0
      iproc=1
      iview=0
      iroll=0
      
!Read parameters from OptsFile (.std)
      itxt=0
      read(11,'(a150)') text
      read(text,*) cardname, itxt      !This line should have the card 'CMS_WAVE_STD' followed by a 1 or 2.  If not, the original format will be read.
      
      select case (cardname)
      case('CMS_WAVE_STD')
        if(itxt.eq.2) then             !New format effort to be more readable - card first, then value(s)  MEB 10/18/21
          call diag_print_message('Reading Options File with Card Format Version 2',' ')
          do  
            read(11,*,iostat=ierr) cardname
            if(ierr .ne. 0) exit
            if(cardname(1:14)=='END_PARAMETERS') exit
            if(cardname(1:1)=='!' .or. cardname(1:1)=='#') cycle
            call cmswave_cards(cardname,foundcard)
            if(foundcard) then
              cycle
            else
              write(*,*) 'Card in CMS-Wave options file not recognized: ',trim(cardname)
            endif
          enddo
        else                           !Previous format - value(s), then card  MEB 10/18/21
          call diag_print_message('Reading Options File with Card Format Version 1',' ')
          backspace(11)
          do 
            read(11,'(a80)',iostat=ierr) text
            if (ierr .ne. 0) exit
            do 138 k=1,75
              kk3=k+3
              kk4=k+4
              kk5=k+5
              if(text(k:kk4).eq.'!free') then
                goto 139
              else if(text(k:kk4).eq.'!iprp') then
                read(text,*)iprpp
                goto 139
              else if(text(k:kk4).eq.'!icur') then
                read(text,*)icur
                goto 139
              else if(text(k:kk4).eq.'!ibre') then
                read(text,*)ibreak
                goto 139
              else if(text(k:kk3).eq.'!irs') then
                read(text,*)irs
                goto 139
              else if(text(k:kk4).eq.'!kout') then
                read(text,*)kout
                if(iabs(kout).ge.1) then
                  read (11,*) (ijsp(1,nn),ijsp(2,nn),nn=1,iabs(kout))
                end if
                goto 139
              else if(text(k:kk5).eq.'!inest') then
                read(text,*)nest
                if(nest.ge.1) then
                  read (11,*) (inest(nn),jnest(nn),nn=1,nest)
                end if
                goto 139
              else if(text(k:kk4).eq.'!ibnd') then
                read(text,*)ibnd
                goto 139
              else if(text(k:kk4).eq.'!iwet') then
                read(text,*)iwet
                goto 139
              else if(text(k:kk3).eq.'!ibf') then
                read(text,*)ibf
                goto 139
              else if(text(k:kk4).eq.'!iark'.and.text(kk5:kk5).ne.'r')then
                read(text,*)iark
                goto 139
              else if(text(k:kk5).eq.'!iarkr') then
                read(text,*)iarkr
                goto 139
              else if(text(k:kk4).eq.'!akap') then
                read(text,*)akap
                goto 139
              else if(text(k:k+2).eq.'!bf') then
                read(text,*)bf
                goto 139
              else if(text(k:kk3).eq.'!ark'.and.text(kk4:kk4).ne.'r') then
                read(text,*)ark
                goto 139
              else if(text(k:kk4).eq.'!arkr') then
                read(text,*)arkr
                goto 139
              else if(text(k:kk4).eq.'!iwvb') then
                read(text,*)iwvbk
                goto 139
              else if(text(k:kk5).eq.'!nonln') then
                read(text,*)nonln
                goto 139
              else if(text(k:kk5).eq.'!igrav') then
                read(text,*)igrav
                goto 139
              else if(text(k:kk5).eq.'!irunu') then
                read(text,*)irunup
                goto 139
              else if(text(k:kk4).eq.'!imud') then
                read(text,*)imud
                goto 139
              else if(text(k:kk4).eq.'!iwnd') then
                read(text,*)iwnd
                goto 139
              else if(text(k:kk5).eq.'!isolv') then
                read(text,*)isolv
                goto 139
              else if(text(k:kk5).eq.'!ixmdf') then
                read(text,*)ixmdf
                goto 139
              else if(text(k:kk5).eq.'!iproc') then
                read(text,*)iproc
                goto 139
              else if(text(k:kk5).eq.'!iview') then
                read(text,*)iview
                goto 139
              else if(text(k:kk5).eq.'!iroll') then
                read(text,*)iroll
                goto 139
              end if
138         continue
139         continue
          enddo
        endif  
      case default                     !Very first format - all parameters specified on one line
        read(text,*,iostat=ierr)iprpp, icur, ibreak, irs, kout, ibnd,  &
          iwet,ibf,iark,iarkr,akap,bf,ark,arkr,iwvbk,nonln,igrav,irunup,  &
          imud,iwnd,isolv,ixmdf,iproc,iview,iroll
      end select

!
      iibreak=0
      inquire(file='std.dat',exist=getfile17)
      if(getfile17) then
        open(unit=20,file='std.dat',status='old')
        read(20,'(a150)') text
        close(20)
        read(text,*,end=133,err=133)icc, icc, iibreak, icc, icc1, icc,  &
          iwet,ibf,iark,iarkr,akap,bf,ark,arkr,iwvbk,nonln,igrav,irunup,  &
          imud,iwnd,isolv,ixmdf,iproc,iview,iroll
        kout=icc1
      end if
133   continue
!
      if(ibreak.lt.iibreak) then
        write(*,*) '*** Reset dissipation index to',iibreak,' ***'
        ibreak=iibreak
      end if
      if(akap.gt.4.) then
        write(*,*) '*** diffraction coef akap > 4, reset to 4 ***'
        if(noptset.eq.3)then
          write(9,*) '*** diffraction coef akap > 4, reset to 4 ***'
        endif
        akap=4.
      end if
!     if(akap.lt.1.) then
!       write(*,*) '*** diffraction coef akap < 1, reset to 1 ***'
!       akap=1.
!     end if
      if(bf.gt.1.) then
        write(*,*) '*** bed friction index bf > 1, reset to 1 ***'
        if(noptset.eq.3)then
          write(9,*) '*** bed friction index bf > 1, reset to 1 ***'
        endif
        bf=1.
      end if
!     if(ibf.eq.0) bf=0.001
      if(ibf.eq.0) bf=0.
!     if(bf.lt..001) bf=.001
      if(ark.gt.1.) then
        write(*,*) '*** norm reflect index ark > 1, reset to 1 ***'
        if(noptset.eq.3)then
          write(9,*) '*** norm reflect index ark > 1, reset to 1 ***'
        endif
        ark=1.
      end if
      if(arkr.gt.1.) then
        write(*,*) '*** sea reflect index arkr > 1, reset to 1 ***'
        if(noptset.eq.3)then
          write(9,*) '*** sea reflect index arkr > 1, reset to 1 ***'
        endif
        arkr=1.
      end if
!     if(iview.ge.1.and.iarkr.ne.0) then
!       write(*,*) '*** backward reflection > 0, reset to 0 ***'
!       iarkr=0
!     end if
! standard output files
!
      do l=1,179
        if(SimFile(l:l+1).eq.'  ') exit
      end do
      if(l.lt.12) l=12
      if (SimFile(l-11:l-5).eq.'swsteer') ixmdf=0
!
      if(iabs(kout).ge.1.and.ixmdf.ne.2) then
        if(kout.ge.1) open (10, file = EngOutFile, status = 'unknown')
        open (12, file = ObsFile,status='unknown',access='append')
      end if
!
      dindex=0.
! --- dindex=0 to set dx=dy=dmesh
! --- dindex=dy, so dx and dy can be same or different
! --- dindex=999, given variable dx(i) and dy(j) arrarys
!     at the end of the dep file
!
      read(15,'(a150)') text
      read(text,*,end=330,err=330) ni,nj,dmesh,dindex
330   continue
!
      if (ni.gt.ipmx .or. nj.gt.jpmx) then
        print*,'****************************************************'
        print*,'This version of CMS-Wave is statically dimensioned'
        print*,'and has a limit of NI = 2500, NJ = 2500'
        print*,'in the near future a dynamically dimensioned version'
        print*,'will become available. '
        print*,'****************************************************'
        stop
      endif
!
!  read/write initial depths dep(i, j):
      depmax=0.
      depmax0=0.
      depmin=10000.
      depmin0=10000.
      dstruc33=0.
      tide0=0.
      do j = nj, 1, -1
        read (15, *) (dep0(i, j), i = 1, ni)
        if(depmax0.lt.dep0(1,j)) depmax0=dep0(1,j)
        if(depmin0.gt.dep0(1,j)) depmin0=dep0(1,j)
        do i=1,ni
          if(dep0(i,j).lt.-20.) dep0(i,j)=-20.
          if(iwet.eq.1.and.dep0(i,j).lt..01) dep0(i,j)=-10.
          if(depmax.lt.dep0(i,j)) depmax=dep0(i,j)
          if(depmin.gt.dep0(i,j)) depmin=dep0(i,j)
        end do
      enddo
      hsmin=0.001
      if(depmax.gt.3.) hsmin=.05
      hs13n=depmax0*0.78
      hs13nn=hs13n
!     if(depmax.lt.20.and.nonln.eq.1) then
!     if(depmax.ne.depmin) then
!     nonln=0
!     write(*,*) '*** Max depth < 20 m'
!     write(*,*) '*** No need for wave-wave interaction calc'
!     write(*,*) '*** Reset nonln = 0'
!     end if
!     end if
!
      do i=1,ni
        do j=1,nj
          depin(i,j) = dep0(i,j)  !Alex, save depths *********************
        enddo
      enddo
!
      igetfile4=0
!
      if(noptset.ne.3)then
        inquire(file=SurgeFile,exist=getfile4)
        if(getfile4) then
          open (21,file=SurgeFile,status='old')
          read (21,'(a30)',end=119,err=119) text1
          read (21,'(a150)',end=119,err=119) text
          READ(text,*) ieta_date
!
          igetfile4=1
          read (21,*,end=119,err=119) ((eta(i,j),i=1,ni),j=nj,1,-1)
!
          do i=1,ni
            do j=2,nj-1
              if(abs(eta(i,j)).ge.900.) then
                if(abs(eta(i,j+1)).lt.900..and.abs(eta(i,j-1)).lt.900.) then
                  eta(i,j)=(eta(i,j+1)+eta(i,j-1))/2.
                end if
              end if
            end do
            if(dep0(i,2).eq.dep0(i,1)) eta(i,1)=eta(i,2)
            if(dep0(i,nj-1).eq.dep0(i,nj)) eta(i,nj)=eta(i,nj-1)
          end do
!
          do j=1,nj
            do i=1,ni
              if(abs(eta(i,j)).ge.900.) eta(i,j)=0.
              dep0(i,j)=dep0(i,j)+eta(i,j)
!             if(dep0(i,j).lt.-0.1) dep0(i,j)=-0.1
            end do
          end do
!
119       continue
        end if
      endif 
!      
      do i=1,ni
        do j=1,nj
          etain(i,j) = eta(i,j)  !Alex, save eta *********************
        enddo
      enddo      
!      
      dvarxx=dmesh
      dvaryy=dmesh
!
      if(dindex.ne.0..and.dindex.ne.999.) then
        dvaryy=dindex
        write(*,*) ' '
        write(*,*) 'Constant dx=',dmesh
        write(*,*) 'Constant dy=',dindex
        if(noptset.eq.3)then
          write(9,*) ' '
          write(9,*) 'Constant dx=',dmesh
          write(9,*) 'Constant dy=',dindex
        endif
      end if
!
      if(dindex.eq.999.) then
        write(*,*) ' '
        read(15,*,end=332,err=332) (dvarxx(i),i=1,ni)
        read(15,*,err=332) (dvaryy(j),j=1,nj)
        write(*,*) ' *** Variable dx and dy arrays ***'
        if(noptset.eq.3)then
          write(9,*) ' '
        write(*,*) '*** Variable dx & dy arrays ***'
        write(*,*) ' '
        endif
        go to 338
332     continue
        dvarxx=dmesh
        dvaryy=dmesh
        write(*,*) '*** Incomplete dx dy arrays ***'
        write(*,*) '  So switch to square cell grid'
      end if
338   continue
      close(15)
!
      dvarxxt=0.
      do i=1,ni
        disx(i)=dvarxxt
        dvarxxt=dvarxxt+dvarxx(i)
      end do
      disx(i)=dvarxxt
!
      dvaryyt=0.
      do j=1,nj
        disy(j)=dvaryyt
        dvaryyt=dvaryyt+dvaryy(j)
      end do
      disy(j)=dvaryyt
!
      ijstruc1=0
      ijstruc3=0
      ijstruc4=0
      read(23,*,IOSTAT=IOS) instruc
      DO ijk=1,instruc
        jstruc=0
        kstruc=0
        dummy=-1000.
        read(23,'(a180)',IOSTAT=IOS) text
        IF (IOS.NE.0) EXIT
        read(text,*,IOSTAT=IOS) istruc,jstruc,kstruc,dummy
!       IF (IOS.NE.0) CYCLE
        if(jstruc.eq.0) CYCLE
        if(kstruc.le.1) then
          if(dummy.lt.-999.) dummy=-10.
          dep0(istruc,jstruc)=dummy
          ijstruc1=ijstruc1+1
          istruc1(ijstruc1)=istruc
          jstruc1(ijstruc1)=jstruc
          dstruc1(ijstruc1)=dummy
        end if
!
        if(kstruc.eq.10) then
          dep0(istruc,jstruc)=-10.
        end if
!
        if(kstruc.eq.3.or.kstruc.ge.6) then
          if(dummy.ge..05) then
            if(kstruc.ge.6) then
              cc=dep0(istruc,jstruc)
              if(dummy.gt.dep0(istruc,jstruc)) dep0(istruc,jstruc)=dummy
            end if
            if(dep0(istruc,jstruc).ge.dummy) then
              ijstruc3=ijstruc3+1
              if(kstruc.ge.6) dstruc33(ijstruc3)=cc
              istruc3(ijstruc3)=istruc
              jstruc3(ijstruc3)=jstruc
              k3(ijstruc3)=kstruc
              dstruc3(ijstruc3)=dummy
!             dstruc3 is the floating breakwater draft if k3=3
!             where dstruc3 is the draft
!             dstruc3 is a flexible porous breakwater  if k3=6
!             dstruc3 is a normal permeable breakwater if k3=7
!             where dstruc3 is nominal depth without breakwater
            end if
            if(kstruc.ge.6) dep0(istruc,jstruc)=dep0(istruc,jstruc)+eta(istruc,jstruc)
          end if
        end if
!
        if(kstruc.eq.4.or.kstruc.eq.5) then
          ijstruc4=ijstruc4+1
          istruc4(ijstruc4)=istruc
          jstruc4(ijstruc4)=jstruc
          k4(ijstruc4)=kstruc
          if(dummy.lt.-999.) dummy=-dep0(istruc,jstruc)+eta(istruc,jstruc)
          dstruc44(ijstruc4)=dummy
!         dstruc44 is the breakwater elevation whatever given in the *.struct
          dep0(istruc,jstruc)=eta(istruc,jstruc)-dummy
        end if
      ENDDO
      rewind(23)
!
      ijstruc2=0
      if(irunup.ge.1.and.depmin.lt.-.01) then
        do j=1,nj
          do i=1,ni
            if(dep0(i,j).lt..05.and.dep0(i,j).ge.-3.) then
              ijstruc2=ijstruc2+1
              istruc2(ijstruc2)=i
              jstruc2(ijstruc2)=j
              dstruc2(ijstruc2)=-dep0(i,j)
              if(ijstruc2.ge.149999) exit
            end if
          end do
        end do
        if(ijstruc2.ge.149999) then
          print*,' '
          print*,'***************************************************'
          print*,'************ WARNING  WARNING  WARNING ************'
          print*,'***************************************************'
          print*,'This version of CMS-Wave is statically dimensioned.'
          print*,'It has a limit of total wave run-up cells  = 150000.'
          print*,'The total wave runup cells found is over the limit.'
          print*,'This Run-up setting covers only partial model area.'
          print*,'Suggest switch to non-automatic mode, use *.struct.'
          print*,'***************************************************'
          print*,' '
        end if
        write(*,*) '*** Automatic Runup trigger ***'
      else
        read(23,*,IOSTAT=IOS) instruc
        DO ijk=1,instruc
          jstruc=0
          kstruc=0
          dummy=-1000.
          read(23,'(a180)',IOSTAT=IOS) text
          IF (IOS.NE.0) EXIT
          read(text,*,IOSTAT=IOS) istruc,jstruc,kstruc,dummy
!         IF (IOS.NE.0) CYCLE
          if(jstruc.eq.0) CYCLE
          if(kstruc.eq.2) then
            ijstruc2=ijstruc2+1
            if(ijstruc2.gt.149999) then
              print*,' '
              print*,'***************************************************'
              print*,'This version of CMS-Wave is statically dimensioned.'
              print*,'It has a limit of total wave run-up cells  = 150000.'
              print*,'The input total wave run-up cells is over the limit'
              print*,'Please revise the data in *.struct file and re-run.'
              print*,'***************************************************'
              stop
            end if
            istruc2(ijstruc2)=istruc
            jstruc2(ijstruc2)=jstruc
            if(dummy.ge.0.) then
              dstruc2(ijstruc2)=dummy-eta(istruc,jstruc)
              dep0(istruc,jstruc)=-dstruc2(ijstruc2)
            else
              dstruc2(ijstruc2)=-dep0(istruc,jstruc)
            end if
          end if
        END DO      
        close(23)
      end if

      if(ijstruc1.ge.1) then
        write(*,*) ' '
        write(*,*) '*** Land/reef/trench feature ***'
        write(*,*) '   (need feature depth info,'
        write(*,*) '    can be negative for land)'
        write(*,*) 'Total Struc 1 Cell(s)=',ijstruc1
        if(noptset.eq.3)then
          write(9,*) ' '
          write(9,*) '*** Land/reef/trench feature ***'
          write(9,*) '   (need feature depth info,'
          write(9,*) '    can be negative for land)'
          write(9,*) 'Total Struc 1 Cell(s)=',ijstruc1
        endif  
      end if
      if(ijstruc2.ge.1) then
        write(*,*) ' '
        write(*,*) '*** Wave run-up calculation ***'
        write(*,*) '    (input max beach height,'
        write(*,*) '     or top structure elev.,'
        write(*,*) '     no effect if elev < 0 m)'
        write(*,*) 'Total Struc 2 Cell(s)=  ',ijstruc2
        if(noptset.eq.3)then
          write(9,*) ' '
          write(9,*) '*** Wave run-up calculation ***'
          write(9,*) '    (input max beach height,'
          write(9,*) '     or top structure elev.,'
          write(9,*) '     no effect if elev < 0 m)'
          write(9,*) 'Total Struc 2 Cell(s)=',ijstruc2
        endif 
        irs=2
      end if

      if(ijstruc3.ge.1) then
        ic3=0
        ic6=0
        ic7=0
        do k=1,ijstruc3
          if(k3(k).eq.3) ic3=ic3+1
          if(k3(k).eq.6) ic6=ic6+1
          if(k3(k).eq.7) ic7=ic7+1
        end do

        if(ic3.ge.1) then
          write(*,*) ' '
          write(*,*) '*** Floating Breakwater ***'
          write(*,*) '    (need average draft'
          write(*,*) '     < depth & > 0.05 m)'
          write(*,*) 'Total Struc 3 Cell(s)=',ic3
          if(noptset.eq.3)then
            write(9,*) ' '
            write(9,*) '*** Floating Breakwater ***'
            write(9,*) '    (need average draft'
            write(9,*) '     < depth & > 0.05 m)'
            write(9,*) 'Total Struc 3 Cell(s)=',ic3
          endif
        end if
        if(ic6.ge.1) then
          write(*,*) ' '
          write(*,*) '*** Flexible Porous Breakwater ***'
          write(*,*) '    (need without-breakwater'
          write(*,*) '     depth & > 0.05 m)'
          write(*,*) 'Total Struc 6 Cell(s)=',ic6      
          if(noptset.eq.3)then
            write(9,*) ' '
            write(9,*) '*** Flexible Porous Breakwater ***'
            write(9,*) '    (need without-breakwater'
            write(9,*) '     depth & > 0.05 m)'
            write(9,*) 'Total Struc 6 Cell(s)=',ic6
          endif       
        end if
        if(ic7.ge.1) then
          write(*,*) ' '
          write(*,*) '*** Normal Permeable Breakwater ***'
          write(*,*) '    (need without-breakwater'
          write(*,*) '     depth & > 0.05 m)'
          write(*,*) 'Total Struc 7 Cell(s)=',ic7
          if(noptset.eq.3)then
            write(9,*) ' '
            write(9,*) '*** Normal Permeable Breakwater ***'
            write(9,*) '    (need without-breakwater'
            write(9,*) '     depth & > 0.05 m)'
            write(9,*) 'Total Struc 7 Cell(s)=',ic7
          endif
        end if
      end if
      if(ijstruc4.ge.1) then
        write(*,*) ' '
        write(*,*) '*** Bottom-mound Breakwater ***'
        write(*,*) '    (input crest elevation,'
        write(*,*) '     can be submerged < 0 m)' 
        write(*,*) 'Total Struc 4/5 Cell(s)=',ijstruc4
        if(noptset.eq.3)then
          write(9,*) ' '
          write(9,*) '*** Bottom-mound Breakwater ***'
          write(9,*) '    (input crest elevation,'
          write(9,*) '     can be submerged < 0 m)'
          write(9,*) 'Total Struc 4/5 Cell(s)=',ijstruc4
          write(*,*) ' '
        endif
      end if
!
      if(iark.eq.0) then
        ark = 0.0
      elseif(iark.eq.2) then
        inquire(file=FrflFile,exist=getfile1)
        if(getfile1) then
          write(*,*) ' *** FrflFile Found ***'
          write(*,*) '     Read forward reflection coef file'
          write(*,*) ' '
          if(noptset.eq.3)then
            write(9,*) ' *** FrflFile Found ***'
            write(9,*) '     Read forward reflection coef file'
            write(9,*) ' '
          endif
          open(unit=18,file=FrflFile,status='old')
          read(18,*) kbi,kbj,qmesh
          if(kbi.ne.ni.or.kbj.ne.nj.or.dmesh.ne.qmesh) then
            write(*,*) 'Wrong forward reflection field file'
            close(18)
            stop
          end if
          do j=nj,1,-1
            read(18,*) (reflty(i,j),i=1,ni)
          enddo
          close(18)
          iark=2
        else
          iark=1
        end if
      endif
!
      if(iarkr.eq.0) then
        arkr = 0.0
      elseif(iarkr.eq.2) then
        inquire(file=BrflFile,exist=getfile2)
        if(getfile2) then
          write(*,*) ' *** BrflFile FOUND ***'
          write(*,*) '     Read backward reflection coef file'
          write(*,*) ' '
          if(noptset.eq.3)then
            write(9,*) ' *** BrflFile FOUND ***'
            write(9,*) '     Read backward reflection coef file'
          endif
          open(unit=18,file=BrflFile,status='old')
          read(18,*) kbi,kbj,qmesh
          if(kbi.ne.ni.or.kbj.ne.nj.or.dmesh.ne.qmesh) then
            write(*,*) 'Wrong backward reflection field file'
            if(noptset.eq.3)then
              write(9,*) 'Wrong backward reflection field file'
            endif
            close(18)
            stop
          end if
          do j=nj,1,-1
            read(18,*) (refltx(i,j),i=1,ni)
          enddo
          close(18)
          iarkr=2
        else
          iarkr=1
        end if
      endif
!
      inquire(file=MudFile,exist=getfile7)
      if(getfile7.and.imud.le.0) then
        write(*,*) ' '
        write(*,*) ' *** MudFile Found ***'
        write(*,*) ' '
        if(noptset.eq.3)then
          write(9,*) ' '
          write(9,*) ' *** MudFile Found ***'
          write(9,*) ' '
        endif                   
        open(unit=29,file=MudFile,status='old')
        read(29,*,end=190,err=190) kbi,kbj,qmesh
        if(kbi.ne.ni.or.kbj.ne.nj.or.dmesh.ne.qmesh) then
          write(*,*) 'Wrong mud field file'
          if(noptset.eq.3)then
            write(9,*) 'Wrong mud field file'
          endif
          close(29)
          stop
        end if
        close(29)
        go to 191
      end if
  190 imud=1
      close(29)
191   continue
!
      if(ibf.eq.2 .or. ibf.eq.4) then          !Only try to read the files if a variable type has been chosen, otherwise keep moving.  MEB  11/15/2021
        inquire(file=FricFile,exist=getfile8)
        if(getfile8) then
          write(*,*) ' '
          write(*,*) ' *** Friction File FOUND ***'
          write(*,*) '     Read friction coef file'
          write(*,*) ' '
          open(unit=28,file=FricFile,status='old')
          read(28,*) kbi,kbj
          if(kbi.ne.ni.or.kbj.ne.nj) then
            write(*,*) 'Wrong friction field file'
            close(28)
            stop
          end if
          close(28)
        else
          write(*,*) '' 
          write(*,*) ' *** Friction File expected but NOT FOUND ***'
          stop
        end if
      end if
!
      if(itxt.eq.0) then
        !  read in special output points
        if (iabs(kout) .ge. 1) then
          do nn = 1, iabs(kout)
            read (11, *) ijsp(1,nn), ijsp(2,nn)
          enddo
        endif

        !  read in nesting output points
        nest=0
        read (11,*,end=109) nest
        do nn = 1, nest
          read (11, *) inest(nn), jnest(nn)
        enddo
      endif
!
      if(nest.ne.0) open (13, file = NestFile, status = 'unknown')
!     azimuth read from .sim file in Subroutine STWfiles

109   continue
!
      if(noptset.eq.3) then
        irs = 1
        if(ijstruc2.ge.1) irs=2
!       ixmdf=1
      endif
!
      if(ixmdf.eq.0) then
        open (66, file = WaveFile,   status = 'unknown')
!       write(66,*) IGMX,JGMX,dmesh
        write (66, *) ni, nj, dmesh
      end if
!
! ***
      if(ixmdf.ge.1) go to 131
! ***
      if(ibreak.gt.0)then
        open (17, file = BreakFile, status = 'unknown')
        write (17, *) ni, nj, dmesh
      endif
!
      inquire(file=TotalFile,exist=getfile11)                          !Mitch 03/22/2017
      inquire(file=SwellFile,exist=getfile12)                          !Mitch 03/22/2017
      if(getfile11) iwet=-2
      if(getfile12) iwet=-3
! --- iwet=-1 for swell & sea (2 files)
!     iwet=-2 for total (1 file)
!     iwet=-3 for swell, sea, and total (3 files)
!
      isteer=0
      do l=1,179
        if(SimFile(l:l+1).eq.'  ') exit
      end do
      if(l.lt.12) l=12
      if (SimFile(l-11:l-5).eq.'swsteer') then
        isteer=1
        do ll=l-12,l-18,-1
          if(SimFile(ll:ll).eq.'\') exit
        end do
        if(ll.le.0) ll=0
        read(SimFile(ll+1:l-12),*) iidate
      end if
!
      if(iwet.lt.0) then
        if (isteer.eq.1) then
          if(SimFile(l-13:l-12).eq.'\1') then
            if(iwet.ne.-1) then
              open(95,file='totalFile',status='unknown')
              write (95, *) ni, nj, dmesh
            end if
            if(iwet.ne.-2) then
              open(96,file=SwellFile,status='unknown')                 !Mitch 03/22/2017
              open(97,file=SeaFile,status='unknown')                   !Mitch 03/22/2017
              write (96, *) ni, nj, dmesh
              write (97, *) ni, nj, dmesh
            end if
          else
            if(iwet.ne.-1) then
              open(95,file='totalFile',status='unknown',access='append')
            end if
            if(iwet.ne.-2) then
              open(96,file=SwellFile,status='unknown',access='append') !Mitch 03/22/2017
              open(97,file=SeaFile,status='unknown',access='append')   !Mitch 03/22/2017
            end if
          end if
        else
          if(iwet.ne.-1) then
            open(95,file='totalFile',status='unknown')
            write (95, *) ni, nj, dmesh
          end if
          if(iwet.ne.-2) then
            open(96,file=SwellFile,status='unknown')                     !Mitch 03/22/2017
            open(97,file=SeaFile,status='unknown')                       !Mitch 03/22/2017
            write (96, *) ni, nj, dmesh
            write (97, *) ni, nj, dmesh
          end if
        end if
      end if
!
      if(kout.ge.0) then
        if(irs.ge.1.and.ixmdf.eq.0)then
          open (18, file = RadsFile, status = 'unknown')
          write (18, *) ni, nj, dmesh
          if(irs.ge.2) then
            if(iwet.ne.-2) then
              open(98,file='setup.wav',status='unknown')
              write(98,*) ni,nj,dmesh
            else
              open(98,file='setup.wav',status='unknown',access='append')
              if(isteer.eq.0) then
                write(98,*) ni,nj,dmesh
              else
                if(SimFile(l-13:l-12).eq.'\1') then
                  close(98)
                  open(98,file=SetupFile,status='unknown')
                  write(98,*) ni,nj,dmesh
                end if
              end if
            end if
          end if
        endif
      end if
!
! ***
131   continue
! ***
      if(iprpp.eq.2.and.iview.eq.2) then
        write(*,*) '*** Conflict: iprp(=iview)= 2, reset to 0 ***'
        iprpp=0
      end if
!
      if(iprpp.eq.2.and.iview.eq.1) then
        write(*,*) '*** Given: iview=1 but iprp=2, reset to 0 ***'
        iview=0
      end if
!
      iwave=0
      nestin1=1
      nestin2=1
      getfile18=.false.       !Added to initialize a value for this variable and get past a usage below.  MEB  11/17/2021
      if(iview.ge.1)then
        ! *** keep above line as new compiler will mess up iview without it
        if(iview.eq.1) go to 280
        inquire(file=SpecFile,exist=getfile18)
        if(getfile18) then
          write(*,*) ' *** 2nd Wave FILE FOUND ***'
          write(*,*) '     Read 2nd wave file'
          write(*,*) ' '
          iwave=1
          open(unit=24,file=SpecFile,status='old')
          read(24,'(a80)',end=280) text
          nestin2=1
          read(text,*,end=281,err=281) nff,mdd,nestin2
281       continue
          write(*,*) '2nd Wave File header =',nff,mdd,nestin2
          read(24,*,end=280,err=280) (cc,nn=1,nff)
        end if
280     continue
!
        iwave1=0
        inquire(file='wave1.dat',exist=getfile19)
        if(getfile19) then
          open(unit=26,file='wave1.dat',status='old')
          read(26,*) iwave1
          close(26)
          do l=1,179
            if(SimFile(l:l+1).eq.'  ') exit
          end do
          if(l.lt.14) l=14
          if(SimFile(l-13:l-5).eq.'\1swsteer') then
            iwave1=0
            write(*,*) '*** 2nd Wave Steering No.1 Cycle ***'
          end if
        end if
!
        !You can only use this after checking the existence of a file. The inquiry is skipped over in some cases.  Added an initialization statement above to correct.
        if(getfile18) then
          do l=1,179
            if(SimFile(l:l+1).eq.'  ') exit
          end do
          if(l.lt.14) l=14
          if (SimFile(l-13:l-5).eq.'\2swsteer') iwave1=1
          if (SimFile(l-11:l-5).ne.'swsteer') go to 350
          do i=1,iwave1
            do j=1,nestin2
              READ(24,*,err=410,end=410) icc
              READ(24,*,err=410,end=410) ((cc,MM=1,MDD),NN=1,NFF)
            end do
          end do
!
          read(24,*,err=334,end=334) icc
          backspace(24)
          go to 335
334       continue
          rewind(24)
          read(24,*) nff,mdd
          read(24,*) (cc,nn=1,nff)
          iwave1=0
335       continue
!
          open(unit=25,file='wave1.dat',status='unknown')
          write(25,*) iwave1+1
          write(*,*) ' '
          write(*,*) ' *** Please check Cycle Number ***'
          write(*,*) '     You are in Cycle',iwave1+1
          write(*,*) '     Check Cycle Number is correct not!!'
          write(*,*) ' *** Modify wave1.dat for correct cycle?'
          write(*,*) ' '
          close(25)
        end if
350     continue
      end if
!
      if(iprpp.eq.1.and.iwave.eq.0) iview=0
!
      if(ibnd.eq.0) then
        READ(8,*) NFF,MDD
        if(NFF.GT.NPF)THEN
          write(*,*) 'ERROR: Number of input frequencies ',NFF
          write(*,*) '  is larger than then maximum value of ',NPF
          write(9,*) 'ERROR: Number of input frequencies ',NFF
          write(9,*) '  is larger than then maximum value of ',NPF
          write(*,*) '  Press any key to continue'
          read(*,*)
          stop
        endif
      else
        inquire(file='nest.dat',exist=getfile5)
        if(getfile5) then
          close(8)
          open(unit=8,file='nest.dat',status='old')
          write(*,*) ' '
          write(*,*) ' *** nest.dat FILE FOUND ***'
          write(*,*) '     Read new nested spectral input file'
          write(*,*) ' '
        end if
!
        inest1=0
        inquire(file='nest1.dat',exist=getfile6)
        if(getfile6) then
          open(unit=7,file='nest1.dat',status='old')
          read(7,*) inest1
          close(7)
          do l=1,179
            if(SimFile(l:l+1).eq.'  ') exit
          end do
          if(l.lt.14) l=14
          if (SimFile(l-13:l-5).eq.'\1swsteer') then
            inest1=0
            write(*,*) '*** Start Steering No.1 Cycle ***'
          end if
        end if
!
        if(getfile5) then
          open(unit=7,file='nest1.dat',status='unknown')
          write(7,*) inest1+1
          write(*,*) ' '
          write(*,*) ' *** Please check Cycle Number ***'
          write(*,*) '     You are in Cycle',inest1+1
          write(*,*) '     Check Cycle Number is correct not!!'
          write(*,*) ' *** Modify nest1.dat for correct cycle?'
          write(*,*) ' '
          close(7)
        end if
!
        read(8,*) NFF,MDD,NESTIN1,azimnest
        if(azimuth.ge.330..and.abs(azimnest).le.10.) then
          azimnest=azimnest+360.
        end if
        if(abs(azimuth).le.10.and.azimnest.ge.330.) then
          azimnest=azimnest-360.
        end if
      end if
!
      READ(8,*) (FFCN(NN),NN=1,NFF)
!
      if(kout.ge.1.and.ixmdf.ne.2) then
        if(iarkr.eq.0) then
          write (10, *) nff, mdd
        else
          write (10, *) nff, mdd*2
        end if
        write (10, 9015) (ffcn(k), k = 1, nff)
9015    format (5(E15.7))
      end if
!
!Hot start option
      n=0  !Initialize N=0, use N later to write out correct ITMS
      if(noptset.eq.3 .and. nsteer.eq.1 .and. .not.coldstart)then !Alex
        n=int(ctime/dtsteer)
        write(*,*)
        write(9,*) 
        write(*,*) 'Skipping ',n,' spectra'        
        write(9,*) 'Skipping ',n,' spectra'
!
        if(ibnd.eq.0) then
          do i=1,n
            read(8,*,end=410) eDate,ws,wd,fp,Tide
            idate = int(mod(edate,100000.))
            kdate = int(edate/100000.)
            if(edate.lt.99999999.) then
              kdate=0
              idate=int(edate)
            end if
            read(8,*,end=410) ((DSFD(NN,MM),MM=1,MDD),NN=1,NFF)
          enddo
        else
          if(getfile5) then
            do i=1,n
              do j=1,nestin
                read(8,'(a150)',err=1333,end=1333) test
                READ(8,*,err=1333,end=1333) ((DSFD(NN,MM),MM=1,MDD),NN=1,NF)
              end do
            end do
          else
            do i=1,n
              do j=1,nestin
                read(8,*)
                READ(8,*,end=410) ((DSFD(NN,MM),MM=1,MDD),NN=1,NF)
              end do
            end do
          end if
          go to 1343
!
1333      continue
          keepi2=i-2
          write(*,*) ' '
          write(*,*) ' *** Reset Input Spectrum to No.', keepi2+1
          write(*,*) ' '
          if(keepi2.lt.0) keepi2=0
          rewind(8)
          READ(8,*)
          READ(8,*) (FFCN(NN),NN=1,NFF)
          do i=1,keepi2
            do j=1,nestin
              read(8,'(a150)',err=1333,end=1333) test
              READ(8,*,err=1333,end=1333) ((DSFD(NN,MM),MM=1,MDD),NN=1,NF)
            end do
          end do
!
1343      continue
          if(getfile5) then
            read(8,'(a150)',err=1333,end=1333) text
          else
            read(8,'(a150)',end=420) text
          end if
          hs13(1)=0.
          read(text,*,end=1421,err=1421) eDate,ws,wd,fp,Tide,xc(1),yc(1),hs13(1)
          if(abs(hs13(1)).gt.900.) hs13(1)=0.
          idate = int(mod(edate,100000.))
          kdate = int(edate/100000.)
          if(edate.lt.99999999.) then
            kdate=0
            idate=int(edate)
          end if
1421      continue
!
        end if
!
        if(iview.eq.2) then
          if(getfile18) then
            do i=1,n
              do j=1,nestin2
                READ(24,*,err=410,end=410) icc
                READ(24,*,err=410,end=410) ((cc,MM=1,MDD),NN=1,NFF)
              end do
            end do
!
            read(24,*,err=1334,end=1334) icc
            backspace(24)
            go to 1335
1334        continue
            rewind(24)
            read(24,*) nff,mdd
            read(24,*) (cc,nn=1,nff)
1335        continue
          end if
        endif
      endif
!
      inst=0
      if(nest.ge.1) then
        write (13, *) nff, mdd, nest, azimuth
        write (13, 9015) (ffcn(k), k = 1, nff)
        if(iview.ge.1) then
          dummy=azimuth+180.
          if(dummy.ge.360.) dummy=dummy-360.
          inquire(file='wav.spc',exist=getfile9)
          if(getfile9) then
            do l=1,179
              if(SimFile(l:l+1).eq.'  ') exit
            end do
            if(l.lt.14) l=14
            if (SimFile(l-13:l-5).eq.'\1swsteer') then
              open(30,file='wav.spc',status='unknown')
              write (30, *) nff, mdd, nest, dummy
              write (30, 9015) (ffcn(k), k = 1, nff)
            else
              if(SimFile(l-11:l-5).ne.'swsteer') then
                open(30,file='wav.spc',status='unknown')
                write (30, *) nff, mdd, nest, dummy
                write (30, 9015) (ffcn(k), k = 1, nff)
              else
                open(30,file='wav.spc',status='unknown',access='append')
              end if
            end if
          else
            open(30,file='wav.spc',status='unknown')
            write (30, *) nff, mdd, nest, dummy
            write (30, 9015) (ffcn(k), k = 1, nff)
          end if
        end if
        do l=1,179
          if(SimFile(l:l+1).eq.'  ') exit
        end do
        if(l.lt.14) l=14
        if(SimFile(l-13:l-5).eq.'\1swsteer') then
          inst=1
          open(14,file='nst.dat',status='unknown')
          write (14, *) nff, mdd, nest, azimuth
          write (14, 9015) (ffcn(k), k = 1, nff)
          close(14)
        endif
        inquire(file='nst.dat',exist=getfile9)
        if(getfile9) then
          inst=1
          open(14,file='nst.dat',status='unknown',access='append')
        end if
      end if
!
      itms = 0 + n  !N is number of times that were skipped during hot start or 0 if cold start.
      imod = 0
      igetfile20 = 0
!
      inquire(file=ShipFile,exist=getfile20)
      if(getfile20) then
        write(*,*) ' '
        write(*,*) ' *** ShipFile Found ***'
        write(*,*) ' '
        igetfile20=1
        open(unit=39,file=ShipFile,status='old')
        read(39,*) nship
        if(nship.eq.0) then
          write(*,*) 'Wrong Shiptrack file'
!         call PressReturn('ERROR:')
          close(39)
          stop
        end if
        read(39,*) (shipL(n),shipB(n),shipD(n),shipS(n),n=1,nship)
      end if
!
      if(iprpp.eq.1.or.iprpp.eq.-2) then
        if(iwind.ge.1) then
          iwind=0
          close(27)
        end if
      end if
      if(iwnd.eq.1) then
        iwind=0
        close(27)
      end if
!
      iwbk=iwvbk
70    continue
      iwvbk=iwbk
      hs13n=hs13nn
!
      if(iproc.eq.0) iproc=1
      if(isolv.ne.0) iproc=1
      if(iproc.ge.2) isolv=0
!
!$    NTHR = omp_get_max_threads()
      NTHR = min(NTHR,ni/300)
!
      IF (IPROC .GT. 1 .and. NTHR.gt.1)THEN
        IF (IPROC .GT. NTHR) THEN
          IPR = NTHR
        ELSE
          IPR = IPROC
        ENDIF
        IPROC = IPR
      ENDIF
!$    if(IPROC.ge.1) call omp_set_num_threads(IPROC) !Alex, set in CMS-Flow (initia.f90)
!
      write(*,*) 'IPROC=', iproc
!
      if(noptset.eq.3 .and. irs.eq.0) irs = 1   !Alex, always output radiation stress if in steering
!
!Reworked this section to always output to screen and diagnostic file.  MEB 10/19/2021
      do i=1,2
        write(iunit(i),*) 'iprp    icur    ibk     irs    kout     ibnd     iwet    ibf     iark    iarkr'
        write(iunit(i),'(3(i3,5x),i3,3x,i5,5x,5(i3,5x))') &
                  iprpp,  icur,   ibreak,  irs,   kout,  &
                  ibnd,   iwet,    ibf,    iark,   iarkr
        write(iunit(i),*)
        write(iunit(i),*) 'akap    bf      ark     arkr   iwvbk    nonln    igrav   irunup  imud    iwnd '
        write(iunit(i),'(4(f7.4,1x),6(i3,5x))') akap, bf, ark, arkr, iwvbk, nonln,igrav,irunup,imud,iwnd
        write(iunit(i),*)
        write(iunit(i),*) 'isolv   ixmdf   iproc   iview   iroll'
        write(iunit(i),'(6(i3,5x))') isolv,ixmdf,iproc,iview,iroll
        write(iunit(i),*)
        if (gamma_bj78 .gt. -1) then
          write(iunit(i),'(A,F6.2)') ' User defined gamma value for Battjes-Janssen 1978: ', gamma_bj78
        endif
        write(iunit(i),*) ' '
      enddo
!
      if(noptset.eq.3 .and. irs.eq.0) irs = 1   !Alex, always output radiation stress if in steering
      if(noptset.eq.3 .and. ijstruct2.ge.1) irs = 2
!      
      itms = itms + 1
      if(getfile5) then
        if(itms.eq.2) then
          write(*,*) ' '
          write(*,*) 'inside time loop, itms = ', itms
          write(*,*) ' '
!         call PressReturn('SUCCESS')
          if(isteer.eq.1) stop
        end if
      end if
      iplane=0
!
770   continue
      if(iview.ge.1) iplane=iplane+1
      if(iplane.eq.1) imod = 0
      ibk3=0
      nestin=nestin1
      if(iplane.eq.2) nestin=nestin2
!
      aslop=1.
      ijb=1
      ibr=0
      cgp=0.0
      wdd=0.0
      fsp=0.0
      sxx=0.0
      syy=0.0
      sxy=0.0
      sxxx=0.0
      sxyx=0.0
      sxyy=0.0
      syyy=0.0
      amud=0.0
      diss=0.0
      dsfd=0.0
      deps=999.
      kstruc4=0
      iprp=iprpp
      wcc=1.
      ismall=10000
      imd=1
      irs0=0
      ichoice=0
!
      if(iplane.eq.0) then
        h13r=0.0
        t13r=1.0
        dmnr=0.0
        h13f=0.0
        t13f=1.0
        dmnf=0.0
      end if
!
      if(iplane.le.1) then
        u1 =0.0
        v1 =0.0
        u10=0.0
        v10=0.0
        sw13=0.0
        sa13=0.0
        tw13=0.0
        ta13=0.0
        dw13=0.0
        da13=0.0
        wxrs=0.0
        wyrs=0.0
      end if
!
      write(*,*) ' '
      print *, 'inside time loop, itms = ', itms
      write(*,*) ' '
!      
      if(noptset.eq.3) write(9,*) 'inside time loop, itms = ', itms
!
      nf=nff
      md=mdd
      fcn=ffcn
      md2=(md+1)/2
!      if(iwind.ge.2) iwind=1
      if(iplane.eq.2) then
        wd=wd/rad
        if(wd.lt.0.) wd=wd+360.
        wd=wd-180.
        if(iwave.eq.1) then
          read(24,'(a150)',end=424) text
          READ(text,*,err=423,end=424) JDATE,ws1,wd1,fp1,Tide1,xc(1),yc(1),hs13(1)
          if(abs(hs13(1)).gt.900.) hs13(1)=0.
423       continue
          if(abs(ws1).lt..001) wd1=0.
          ws=ws1
          wd=wd1
          if(abs(fp1).gt..001) fp=fp1
          if(iprp.eq.1.or.iprp.eq.-2) ws=0.
          write(*,*) 'Wave Spc Input Index =',JDate
          go to 424
        else
          hs0=0.
          sum=0.
          if(abs(wd).ge.120..and.iwind.eq.0) then
            h13=0.
            t13=0.
            dmn=0.
            imod=2
            irs0=1
!
            do j = 1, nj
              do i = 1, ni
                d1(i, j) = dep0(ni-i+1, nj-j+1) + Tide
              enddo
              dvary(j)=dvaryy(nj-j+1)
            enddo
!
            do i=1,ni
              dvarx(i)=dvarxx(ni-i+1)
            end do
!
            go to 80
          end if
          go to 425
        end if
      end if
!
      if(noptset.ne.3)then
        if(itms.ge.2) then
          if(getfile4) then !Alex
            read(21,'(a150)',end=219,err=219) text
            READ(text,*) ieta_date
            open (15,file=DepFile,status='old')
            read (15,*)
            do j=nj,1,-1
              read (15,*) (dep0(i,j),i=1,ni)
              read (21,*,end=219,err=219) (eta(i,j),i=1,ni)
              do i=1,ni
                if(dep0(i,j).lt.-20.) dep0(i,j)=-20.
                if(iwet.eq.1.and.dep0(i,j).lt..01) dep0(i,j)=-10.
              end do
            end do
            close(15)
!
            do i=1,ni
              do j=2,nj-1
                if(abs(eta(i,j)).ge.900.) then
                  if(abs(eta(i,j+1)).lt.900..and.abs(eta(i,j-1)).lt.900.) then
                    eta(i,j)=(eta(i,j+1)+eta(i,j-1))/2.
                  end if
                end if
              end do
              if(dep0(i,2).eq.dep0(i,1)) eta(i,1)=eta(i,2)
              if(dep0(i,nj-1).eq.dep0(i,nj)) eta(i,nj)=eta(i,nj-1)
            end do
!
            do j=1,nj
              do i=1,ni
                if(eta(i,j).lt.-900.0) eta(i,j)=0.0
                dep0(i,j)=dep0(i,j)+eta(i,j)
!               if(dep0(i,j).lt.-0.1) dep0(i,j)=-0.1
              enddo
            enddo
!          
219         continue
          else
            eta=0.0
!           do j=1,nj
!             do i=1,ni
!               if(dep0(i,j).lt.-0.1) dep0(i,j)=-0.1
!             enddo
!           enddo
          end if
        end if   
      else
!
        !Save depths and wse before updating the depth to include eta
        do i=1,ni
          do j=1,nj
            dep0(i,j)=depin(i,j)          
            eta(i,j)=etain(i,j)
            if(abs(eta(i,j)).ge.900.) eta(i,j)=0.0
            dep0(i,j)=dep0(i,j)+eta(i,j)
          enddo
        enddo
      endif
!!!
      if(ibnd.eq.0) then
        READ(8,*,end=420) eDate,ws,wd,fp,Tide
        idate = int(mod(edate,100000.))
        kdate = int(edate/100000.)
        if(edate.lt.99999999.) then
          kdate=0
          idate=int(edate)
        end if
      else
        if(getfile5) then
          do i=1,inest1
            do j=1,nestin
              read(8,'(a150)',err=333,end=333) test
              READ(8,*,err=333,end=333) ((DSFD(NN,MM),MM=1,MDD),NN=1,NF)
            end do
          end do
        else
          do i=1,inest1
            do j=1,nestin
              read(8,*)
              READ(8,*,end=410) ((DSFD(NN,MM),MM=1,MDD),NN=1,NF)
            end do
          end do
        end if
        go to 343
!
333     continue
        keepi2=i-2
        write(*,*) ' '
        write(*,*) ' *** Reset Input Spectrum to No.', keepi2+1
        write(*,*) ' '
        if(keepi2.lt.0) keepi2=0
        rewind(8)
        READ(8,*)
        READ(8,*) (FFCN(NN),NN=1,NFF)
        do i=1,keepi2
          do j=1,nestin
            read(8,'(a150)',err=333,end=333) test
            READ(8,*,err=333,end=333) ((DSFD(NN,MM),MM=1,MDD),NN=1,NF)
          end do
        end do
!
343     continue
        if(getfile5) then
          read(8,'(a150)',err=333,end=333) text
        else
          read(8,'(a150)',end=420) text
        end if
        hs13(1)=0.
        read(text,*,end=421,err=421) &
        eDate,ws,wd,fp,Tide,xc(1),yc(1),hs13(1)
        if(abs(hs13(1)).gt.900.) hs13(1)=0.
        idate = int(mod(edate,100000.))
        kdate = int(edate/100000.)
        if(edate.lt.99999999.) then
          kdate=0
          idate=int(edate)
        end if
421     continue
      end if
!    
      if(iprp.eq.1.or.iprp.eq.-2) ws=0.
      if(kdate.gt.0) then
        if(idate.le.9999) then
          idate1=mod(kdate,10)*100000+idate
          kdate1=kdate/10
          write(*,'("Wave Input Index =",i7,i6)') kdate1,idate1
        else
          write(*,'("Wave Input Index =",i8,i5)') kdate,idate
        end if
      else
        if(idate.lt.itms) idate=itms
        if(idate.lt.iidate) idate=iidate
        write(*,'("Wave Input Index =",i8)') idate
      end if
      write(*,*) ' '
!
      if(imod.eq.1.or.imod.eq.3) then
        nc=ni
        ni=nj
        nj=nc
      end if
      imod=0
!
424   continue
!
      iengspc=8
      if(iplane.eq.2) iengspc=24
      READ(iengspc,*,end=410) ((DSFD(NN,MM),MM=1,MDD),NN=1,NF)
!
      DTH = 180./float(MDD+1)
      if(iview.ge.1) then
        DTH = 180./float(MDD)
        dsfd=dsfd*float(mdd)/float(mdd+1)
      end if
      DTH = DTH*RAD
! 
      if(wd.gt.270.) wd=wd-360.
      if(wd.lt.-270.) wd=wd+360.
!
      DO NN=2,NF-1
        DF(NN)=0.5*(FFCN(NN+1)-FFCN(NN-1))
      ENDDO
!     DF(1)=0.5*(FFCN(2)+FFCN(1))
      DF(1)=DF(2)
      DF(NF)=DF(NF-1)
      DFINP=DF
!
      sum=0.
      do mm=1,mdd
        do nn=1,nf
          sum=sum+dsfd(nn,mm)*df(nn)
        ENDDO
      ENDDO
      hs0=4.*sqrt(sum*dth)

425   continue

      if(hs0.lt..01.and.fp.gt.fcn(1)) fp=0.04
      if(ibnd.ge.1) go to 435
      if(hs0.lt..02.and.ws.ge..1) then
        nf1=1
        do nn=nf,1,-1
          if(fcn(nn).ge..1) nf1=nn
        end do
        if(nf1.ne.1) then
          do nn=nf1,nf
            fcn(nn-nf1+1)=fcn(nn)
          end do
        end if
        nf=nf-nf1+1
        if(fcn(1).ne..1) then
          do nn=nf,1,-1
            fcn(nn+1)=fcn(nn)
          enddo
          fcn(1)=.1
          nf=nf+1
        end if
!
        DO NN=2,NF-1
          DF(NN)=0.5*(FCN(NN+1)-FCN(NN-1))
        ENDDO
        DF(1)=DF(2)
        DF(NF)=DF(NF-1)
        DFINP=DF
      end if
!
435   continue
!
      if(iview.eq.1.and.iplane.eq.1) then
        if(ws.lt.2..or.cos(wd*rad).gt..3) iplane=0
      end if
!
      if(hs0.lt..01.and.abs(ws-.5).lt.0.49) ws=1.
      if(iview.ge.1.and.abs(ws-.5).lt.0.49) ws=1.
      if(iprp.le.-1.and.mdd.eq.35) then
        write(*,*) '********************** '
        write(*,*) '* FAST MODE SELECTED * '
        write(*,*) '********************** '
!
        md=5
        if(ws.ge..1) md=7
        immd=mdd/md
        do mm=1,mdd
          newmm=(mm-1)/immd+1
          do nn=1,nf
            if(mod(mm-1,immd).eq.0) then
              dsfd(nn,newmm)=dsfd(nn,mm)
            else
              dsfd(nn,newmm)=dsfd(nn,newmm)+dsfd(nn,mm)
            end if
          end do
        end do
        dsfd=dsfd*float(md+1)/float(mdd+1)
        DTH = 180./float(MD+1)
        if(iview.ge.1) then
          DTH = 180./float(MD)
          dsfd=dsfd*float(md)/float(md+1)
        end if
        DTH = DTH*RAD
        md2=(md+1)/2
        ichoice=2
        write(*,*) 'New Total Direction Bins =',md
      end if
!
      cc=HPAI
      if(iview.ge.1) cc=cc+DTH/2.
      DO MM=1,MD
!       *** should above be mm=imd,md but no big deal (get what needed)
        DCM(MM)=DTH*MM-cc
        cosa(mm)=cos(dcm(mm))
        sina(mm)=sin(dcm(mm))
      ENDDO
!
      sum=0.
      do mm=imd,md
        do nn=1,nf
          sum=sum+dsfd(nn,mm)
        ENDDO
      ENDDO
!
      if(sum.lt..0001) then
        nnf=nf
        ph0=g/ws*.9
        if(ph0.gt..36) ph0=.36
        do nn=nnf,1,-1
          if(fcn(nn).gt.ph0) nnf=nn
        end do
        dsfd(nnf,(imd+md)/2)=0.0001
      end if
!
      if(ibnd.ge.1) then
        wc(1)=0.
        wcx=0.
        wcy=0.
        do mm=imd,md
          do nn=1,nf
            wc(1)=wc(1)+dsfd(nn,mm)
            wcx=wcx+dsfd(nn,mm)*cosa(mm)
            wcy=wcy+dsfd(nn,mm)*sina(mm)
          end do
        end do
        wdx=wcx
        wdy=wcy
        wcd(1)=atan2(wcy,wcx)
        do kn=2,nestin
          read(iengspc,'(a150)',end=420) text
          hs13(kn)=0.
          read(text,*,end=422,err=422) eDate,ws1,wd1,fp1,Tide1,xc(kn),yc(kn),hs13(kn)
          if(abs(hs13(kn)).gt.900.) hs13(kn)=0.
          idate = int(mod(edate,100000.))
          kdate = int(edate/100000.)
          if(edate.lt.99999999.) then
            kdate=0
            idate=int(edate)
          end if
          if(kdate.eq.0) then
            if(idate.lt.itms) idate=itms
            if(idate.lt.iidate) idate=iidate
          end if
!
422       continue
          if(ws1.gt.ws) ws=ws1
          if(fp1.lt.fp) fp=fp1
          if(iplane.eq.2) fp=fp1
          if(tide1.gt.tide) tide=tide1
          read(iengspc,*,end=410) ((wk2(1,nn,mm),mm=1,mdd),nn=1,nf)
!
          if(hs0.lt..001) wk2(1,nf,md2)=1.0
!
          if(ichoice.eq.2) then
            do mm=1,mdd
              newmm=(mm-1)/immd+1
              do nn=1,nf
                if(mod(mm-1,immd).eq.0) then
                  wk2(1,nn,newmm)=wk2(1,nn,mm)
                else
                  wk2(1,nn,newmm)=wk2(1,nn,newmm)+wk2(1,nn,mm)
                end if
              end do
            end do
            wk2=wk2*float(md+1)/float(mdd+1)
          end if
!
          wc(kn)=0.
          wcx=0.
          wcy=0.
          do mm=imd,md
            do nn=1,nf
              dsfd(nn,mm)=dsfd(nn,mm)+wk2(1,nn,mm)
              wc(kn)=wc(kn)+wk2(1,nn,mm)
              wcx=wcx+wk2(1,nn,mm)*cosa(mm)
              wcy=wcy+wk2(1,nn,mm)*sina(mm)
            end do
          end do
          wdx=wdx+wcx
          wdy=wdy+wcy
          wcd(kn)=atan2(wcy,wcx)
        end do
        wwd=atan2(wdy,wdx)
        do mm=imd,md
          do nn=1,nf
            dsfd(nn,mm)=dsfd(nn,mm)/float(nestin)
          end do
        end do
        sum=0.
        summ=0.
        do kn=1,nestin
          sum=sum+wc(kn)
          summ=summ+hs13(kn)**2
          wcd(kn)=(wcd(kn)-wwd)/dth
        end do
        sum=sum/float(nestin)
        hs13n=sqrt(summ/float(nestin))
        do kn=1,nestin
          wc(kn)=wc(kn)/sum
        end do
!
        ikn=nestin+1
        xc(ikn)=x0
        yc(ikn)=y0
        cc1=(xc(ikn)-xc(1))**2+(yc(ikn)-yc(1))**2
        wc(ikn)=wc(1)
        wcd(ikn)=wcd(1)
        do kn=2,nestin
          cc=(xc(ikn)-xc(kn))**2+(yc(ikn)-yc(kn))**2
          if(cc.lt.cc1) then
            wc(ikn)=wc(kn)
            wcd(ikn)=wcd(kn)
          end if
        end do
!
        ydis=0.
        do j=1,nj
          ydis=ydis+dvaryy(j)
        ENDDO  !613 continue
        ikn=nestin+2
        xc(ikn)=x0-ydis*sinaz
        yc(ikn)=y0+ydis*cosaz
        cc1=(xc(ikn)-xc(1))**2+(yc(ikn)-yc(1))**2
        wc(ikn)=wc(1)
        wcd(ikn)=wcd(1)
        do kn=2,nestin
          cc=(xc(ikn)-xc(kn))**2+(yc(ikn)-yc(kn))**2
          if(cc.lt.cc1) then
            wc(ikn)=wc(kn)
            wcd(ikn)=wcd(kn)
            cc1=cc
          end if
        end do
!
        ydis=dvaryy(1)/2.
        if(iplane.eq.2) ydis=dvaryy(nj)/2.
        x0p=x0+dvarxxt*cosaz-dvaryyt*sinaz
        y0p=y0+dvarxxt*sinaz+dvaryyt*cosaz
        do j=1,nj
          if(iplane.le.1) then
            if(j.ge.2) ydis=ydis+dvaryy(j-1)
            xcc=x0-ydis*sinaz
            ycc=y0+ydis*cosaz
          else
            if(j.ge.2) ydis=ydis+dvaryy(nj+2-j)
            xcc=x0p+ydis*sinaz
            ycc=y0p-ydis*cosaz
          end if
!
          sum=0.
          sum1=0.
          sum2=0.
!
          do kn=1,nestin+2
            cc=(xcc-xc(kn))**2+(ycc-yc(kn))**2
            if(cc.lt.25.) then
              wcc(j)=wc(kn)
              dcd(j)=wcd(kn)
!             CYCLE
              go to 612
            end if
            sum=sum+1./cc
            sum1=sum1+wc(kn)/cc
            sum2=sum2+wcd(kn)/cc
          end do
          wcc(j)=sum1/sum
          dcd(j)=sum2/sum
  612     continue
!         write(*,*) 'j,xcc,ycc,wcc=',j,xcc,ycc,wcc(j),dcd(j)
        end do
!
        cc=azimnest-azimuth
        wd=wd+cc
        if(cc.gt.360.) cc=cc-360.
        if(cc.lt.-360.) cc=cc+360.
        azim=cc*float(mdd+1)/180.
        if(ichoice.eq.2) azim=cc*float(md+1)/180.
        if(abs(azim).lt..01) go to 21
        write(*,*) ' '
        write(*,*) 'Child Grid total rotation bin(s) =', nint(azim)
        write(*,*) ' '
        iazim=int(azim)
        DO nn=1,nf  !do 28
          if(iabs(iazim).lt.1) go to 24
          wk2=0.
          do mm=imd,md  !do 23
            m1=mm+iazim
            if(m1.lt.imd) m1=imd
            if(m1.gt.md) m1=md
            wk2(1,nn,m1)=wk2(1,nn,m1)+dsfd(nn,mm)
          ENDDO  !23 continue
          do mm=imd,md  !do 25
            dsfd(nn,mm)=wk2(1,nn,mm)
          ENDDO  !25 continue
   24     continue
          dcdmm=azim-float(iazim)
          if(abs(dcdmm).lt..01) go to 30
          if(dcdmm.gt.0.) then
            c1=1.-dcdmm
            do mm=imd+1,md  ! do 37
              wk2(1,nn,mm)=dsfd(nn,mm)*c1+dsfd(nn,mm-1)*dcdmm
            ENDDO  !37 continue
            wk2(1,nn,imd)=dsfd(nn,imd)*c1
          else
            dcdmm=-dcdmm
            c1=1.-dcdmm
            do mm=imd,md-1  !do 27
              wk2(1,nn,mm)=dsfd(nn,mm)*c1+dsfd(nn,mm+1)*dcdmm
            ENDDO  ! 27 continue
            wk2(1,nn,md)=dsfd(nn,md)*c1
          end if
          do mm=imd,md  !do 29
            dsfd(nn,mm)=wk2(1,nn,mm)
          ENDDO  ! 29 continue
   30     continue
        ENDDO  ! 28 continue
   21   continue
      end if
!!!!
      big=0.
      ibig=1
      do mm=imd,md
        sum=0.
        do nn=1,nf
          sum=sum+dsfd(nn,mm)
        end do
        if(sum.gt.big) then
          big=sum
          ibig=mm
        end if
      end do
      prd=dcm(ibig)
      dbar=prd
!     CLOSE(8)
!
! --- Reset starting frequency and frequency range so
!     skip those high frequency bins with negligible
!     energy and also those low frequency with very
!     little energy
!
      DO NN=2,NF-1
        DF(NN)=0.5*(FCN(NN+1)-FCN(NN-1))
      ENDDO
!     DF(1)=0.5*(FCN(2)+FCN(1))
      DF(1)=DF(2)
      DF(NF)=DF(NF-1)
      DFINP=DF
!
      sum=0.
      do mm=imd,md
        do nn=1,nf
          sum=sum+dsfd(nn,mm)*df(nn)
        ENDDO
      ENDDO
      hs0=4.*sqrt(sum*dth)
!
!     write(*,*) 'iview,iplane=',iview,iplane
!     write(*,*) ' hs0,  hs13n=',hs0,hs13n
!
      if(iview.eq.1.and.iplane.eq.2) hs13n=hs13nn
!
      if(hs13n.lt..0001) hs13n=hs0
      if(ibnd.ge.1.or.hs0.gt.hs13n) then
        cc=(hs13n/hs0)**2
        hs0=hs13n
        do mm=imd,md
          do nn=1,nf
            dsfd(nn,mm)=dsfd(nn,mm)*cc
          end do
        end do
      end if
!
      if(iplane.eq.2) go to 426
!
      write(*,*) 'Original 4*sqrt(E) :',hs0,'m'
      if(noptset.eq.3) write(9,*) 'Original 4*sqrt(E) :',hs0,'m'
      write(*,*) 'isteer,iidate =',isteer,iidate
      hs1=hs0
!
      wsmag = ws
!
      write(*,*) 'Wsmag =', wsmag
!
      if(iwind.ge.1) ws=1.
      if(ws.lt..1.and.isteer.eq.0) then
        if(hs0.lt.hsmin) then
          igmx=ni
          jgmx=nj
          d1=dep0+Tide
          H13=0.01
          T13=1.
          DMN=0.
          SOP=0.
          H13R=0.
          T13R=1.
          DMNR=0.

          j1=ni/2
        
 9999     format ('Column ', i0)
 9020     format('Average wave height = ',f10.4)
          do iii=1,ni
            write(*,9999) iii
            write(*,9020) h13(iii,j1)
          end do
!
          go to 90
        end if
      end if
!
426   continue
      if(iplane.eq.1) go to 427
      if(iplane.eq.2) then
        write(*,*) ' '
        write(*,*) '*** Set Local Upper-Right Corner as Origin'
        imod=2
        ! It is the easiest one: no need to change ni & nj numbers
        ! define d1(1,1) is dep0(ni,nj), d1(2,1) is dep0(ni-1,nj),
        ! d1(3,1) is dep0(ni-2,nj), d1(1,2) is dep0(ni,j1-1), etc.
        hs1=hs0
        go to 427
      end if
!
      if(ibnd.ge.1.and.hs0.lt..01) hs0=.02
      if(ws.ge..1.and.hs0.lt..01) iprp=2
!     if(ws.ge..1.and.hs0.lt..01) then
      if(ws.ge..1.and.iprp.eq.2) then
        if(wd.gt.180.) wd=wd-360.
        if(wd.lt.-180.) wd=wd+360.
!
        depimod0=0.
        conimod0=0.
        do j=1,nj
          if(dep0(1,j).gt.0.) then
            depimod0=depimod0+dep0(1,j)
            conimod0=conimod0+1.
          end if
        ENDDO
        if(conimod0.ge.1.) then
          depimod0=depimod0/conimod0
        else
          depimod0=10.
        end if
!
        if(abs(wd).le.45.) then
          write(*,*) '*** Set Local Lower-Left Corner as Origin'
          imod=0
        else if(abs(wd-90.).lt.45.) then
          write(*,*) '*** Set Local Lower-Right Corner as Origin'
          imod=1
        ! So, d1(1,1) is dep0(ni,1), d1(2,1) is dep0(ni,2),
        !     d1(3,1) is dep0(ni,3), ..  d1(1,2) is dep0(ni-1,1),
        !     d1(1,ni) is dep0(1,1), ..  d1(nj,ni) is dep0(1,nj), etc.
          wd=wd-90.
!
          depimod1=0.
          conimod1=0.
          do i=1,ni
            if(dep0(i,1).gt.0.) then
              depimod1=depimod1+dep0(i,1)
              conimod1=conimod1+1.
            end if
          ENDDO
          if(conimod1.ge.1.) then
            depimod1=depimod1/conimod1
          else
            depimod1=10.
          end if
!
          cc=sqrt(depimod1/depimod0)
          if(cc.gt.1.) cc=1.
!
          do nn=1,nf
            do mm=1,md2-1
              md2mm=md2+mm
              dsfd(nn,mm)=dsfd(nn,md2mm)*cc
              dsfd(nn,md2mm)=0.
            ENDDO
            if(md2*2.eq.md) then
              dsfd(nn,md2)=dsfd(nn,md)*cc
              dsfd(nn,md)=0.
            else
              dsfd(nn,md2)=0.
            end if
          ENDDO
!
          sum=0.
          do mm=imd,md
            do nn=1,nf
              sum=sum+dsfd(nn,mm)*df(nn)
            ENDDO
          ENDDO
          hs0=4.*sqrt(sum*dth)
          write(*,*) 'Original 4*sqrt(E) =',hs0,'m'
          hs1=hs0
        else if(abs(wd+90.).lt.45.) then
          write(*,*) '*** Set Local Upper-Left Corner as Origin'
          imod=3
          ! So, d1(1,1) is dep0(1,nj), d1(2,1) is dep0(1,nj-1),
          !     d1(3,1) is dep0(1,nj-2), ..  d1(1,2) is dep0(2,nj),
          !     d1(1,ni) is dep0(ni,nj), ..d1(nj,ni) is dep0(ni,1), etc.
          wd=wd+90.
!
          depimod3=0.
          conimod3=0.
          do i=1,ni
            if(dep0(i,nj).gt.0.) then
            depimod3=depimod3+dep0(i,nj)
            conimod3=conimod3+1.
            end if
          ENDDO
          if(conimod3.ge.1.) then
            depimod3=depimod3/conimod3
          else
            depimod3=10.
          end if
!
          cc=sqrt(depimod3/depimod0)
          if(cc.gt.1.) cc=1.
!
           do nn=1,nf
             do mm=imd,md2-1
               md2mm=md2+mm
               dsfd(nn,md2mm)=dsfd(nn,mm)*cc
               dsfd(nn,mm)=0.
             ENDDO
             if(md2*2.eq.md) then
               dsfd(nn,md)=dsfd(nn,md2)*cc
               dsfd(nn,md2)=0.
             else
               dsfd(nn,md)=0.
             end if
           ENDDO
!  
          sum=0.
          do mm=imd,md
            do nn=1,nf
              sum=sum+dsfd(nn,mm)*df(nn)
            ENDDO
          ENDDO
          hs0=4.*sqrt(sum*dth)
          write(*,*) 'Original 4*sqrt(E) =',hs0,'m'
          hs1=hs0
        else
          write(*,*) '*** Set Local Upper-Right Corner as Origin'
          imod=2
          ! It is the easiest one: no need to change ni & nj numbers
          ! define d1(1,1) is dep0(ni,nj), d1(2,1) is dep0(ni-1,nj),
          ! d1(3,1) is dep0(ni-2,nj), d1(1,2) is dep0(ni,j1-1), etc.
          if(wd.lt.0.) wd=wd+360.
          wd=wd-180.
          dsfd=0.
!
          hs0=0.
          write(*,*) 'Original 4*sqrt(E) =',hs0,'m'
          hs1=hs0
        end if
        ! --- done reset wd for different imod index
      end if
      write(*,*) ' '
!     if(hs0.gt.2.) then
!       if(ws.ge.20..and.abs(wd).gt.40.) then
!         iprp=3
!         write(*,*) 'Reset iprp = 3 (no lateral source/sink)'
!       end if
!     end if
!
427   continue
!
! determine water depth for this event (for different imod)
!
      if(imod.eq.0) then
        do j = 1, nj
          do i = 1, ni
            d1(i, j) = dep0(i, j) + Tide
            ijb(i,j) = 1
          enddo
          dvary(j)=dvaryy(j)
        enddo
        do i=1,ni
         dvarx(i)=dvarxx(i)
        end do
      end if
!
      if(imod.eq.2) then
        do j = 1, nj
          do i = 1, ni
            d1(i, j) = dep0(ni-i+1, nj-j+1) + Tide
            ijb(i,j) = 1
          enddo
          dvary(j)=dvaryy(nj-j+1)
        enddo
        do i=1,ni
          dvarx(i)=dvarxx(ni-i+1)
        end do
      end if
!
      if(imod.eq.1) then
        nc=ni
        ni=nj
        nj=nc
        do j = 1, nj
          do i = 1, ni
            d1(i, j) = dep0(nj-j+1, i) + Tide
            ijb(i,j) = 1
          enddo
          dvary(j)=dvarxx(nj-j+1)
        enddo
        do i=1,ni
          dvarx(i)=dvaryy(i)
        end do
      end if
!
      if(imod.eq.3) then
        nc=ni
        ni=nj
        nj=nc
        do j = 1, nj
          do i = 1, ni
            d1(i, j) = dep0(j, ni-i+1) + Tide
            ijb(i,j) = 1
          enddo
          dvary(j)=dvarxx(j)
        enddo
        do i=1,ni
          dvarx(i)=dvaryy(ni-i+1)
        end do
      end if
!
      do j=1,nj
        do i=1,ni
          if(d1(i,j).le.-.01) d1(i,j)=-2.0002
!         modified 23 April 2009 Lihwa
! --- backward reflection did not work if set d1(i,j)=0 here
        end do
      end do
!
      if(hs0.lt..005.and.ws.ge..1) then
        do j=1,nj
          d1(1,j)=-100.
        end do
        prd=wd*rad
        dbar=prd
      end if
!
      do j=1,nj
        if(d1(1,j).lt..01) d1(1,j)=.01
      end do
!
      do j=2,nj-1
        if(d1(ni,j).lt..01) then
          if(d1(ni,j+1).ge..0.and.d1(ni,j-1).ge..0) then
            d1(ni,j)=(d1(ni,j+1)+d1(ni,j-1))/2.
            if(d1(ni,j).lt..01) d1(ni,j)=.01
          end if
        end if
      end do
      if(d1(ni,nj).gt..0.and.d1(ni,nj-1).le..0) d1(ni,nj)=0.
      if(d1(ni,1).gt..0.and.d1(ni,2).le..0) d1(ni,1)=0.
!
      do i=2,ni
        if(d1(i,1).ge..01.and.d1(i,2).lt..01) d1(i,1)=0.
        if(d1(i,nj).ge..01.and.d1(i,nj-1).lt..01) d1(i,nj)=0.
      end do
!
      do j=2,nj-1
        do i=2,ni
          if(d1(i,j).ge..0.and.d1(i,j).lt..01) d1(i,j)=0.
        end do
      end do
!
      do j=2,nj-1
        do i=2,ni
          if(d1(i,j).ge..01.and.d1(i,j+1).le..0.and.d1(i,j-1).le..0) &
          d1(i,j)=0.
        end do
      end do
!
      write(*,*) ' '
      if(ibf.eq.0) goto 402
      if(getfile8) then
        open(unit=28,file=FricFile,status='old')
        read(28,*) kbi,kbj
          if(imod.eq.0) then
            do j=nj,1,-1
              read(28,*,err=401,end=401) (bfric(i,j),i=1,ni)
            enddo
          elseif(imod.eq.2) then
            do j=1,nj
              read(28,*,err=401,end=401) (bfric(i,j),i=ni,1,-1)
            enddo
          elseif(imod.eq.1) then
            do i=ni,1,-1
              read(28,*,err=401,end=401) (bfric(i,j),j=nj,1,-1)
            enddo
          else
            do i=1,ni
              read(28,*,err=401,end=401) (bfric(i,j),j=1,nj)
            enddo
          end if
        close(28)
        go to 400
      else
        go to 402
      end if
401   continue
      write(*,*) ' *** Error reading friction field data ***'
402   continue
      write(*,*) 'Use friction coef=', bf
      write(*,*) ' '
!
      if(bf*hs0.gt..05.and.ibf.eq.1) ibf=3
!
      do j=1,nj
        do i=1,ni
          bfric(i,j)=bf
        end do
      end do
!
400   continue
!
      bfricmax=0.
      do i=1,ni
        do j=1,nj
          if(bfric(i,j).lt..001) bfric(i,j)=.001
          if(ibf.ge.3) then
            cc=d1(i,j)
            if(cc.lt..3) cc=.3
            bfric(i,j)=bfric(i,j)**2/cc**.333333
          end if
          if(bfricmax.lt.bfric(i,j)) bfricmax=bfric(i,j)
        end do
      end do
      
!!!!! Mud
      if(getfile7.and.imud.le.0) then
        open(unit=29,file=MudFile,status='old')
        read(29,*) kbi,kbj
        if(imod.eq.0) then
          do j=nj,1,-1
            read(29,*) (amud(i,j),i=1,ni)
          enddo
        elseif(imod.eq.2) then
          do j=1,nj
            read(29,*) (amud(i,j),i=ni,1,-1)
          enddo
        elseif(imod.eq.1) then
          do i=ni,1,-1
            read(29,*) (amud(i,j),j=nj,1,-1)
          enddo
        else
          do i=1,ni
            read(29,*) (amud(i,j),j=1,nj)
          enddo
        end if
        close(29)
      end if
      
!!!!!Currents
      
      if(noptset.ne.3)then
        if (icur .ge. 1) then
          if(iplane.eq.2) go to 428
          if(itms.eq.1.or.icur.ge.2) then
            open (16, file = CurrFile, status = 'old')
            read(16,*) nic, njc
          end if
          if(imod.eq.0.or.imod.eq.2) then
            if ((nic .ne. ni) .or. (njc .ne. nj)) then
              print *,' Error: current field size does not match depth grid size'
              stop
            endif
          else
            if ((njc .ne. ni) .or. (nic .ne. nj)) then
              print *,' Error: current field size does not match depth grid size'
              stop
            endif
          endif
          
          ! read constant current field
          read(16,'(a150)') text
          READ(text,*) icur_date
          write(*,*) 'Current Index =',icur_date
          write(*,*) ' '
          if(imod.eq.0) then
            do j = nj, 1, -1
              read (16, *) (u1(i, j), v1(i, j), i = 1, ni)
            enddo
          end if
          if(imod.eq.2) then
            do j = 1, nj
              read (16, *) (u1(i, j), v1(i, j), i = ni, 1, -1)
            enddo
          end if
          if(imod.eq.1) then
            do i = ni, 1, -1
              read (16, *) (v1(i, j), u1(i, j), j = nj, 1, -1)
            enddo
          end if
          if(imod.eq.3) then
            do i = 1, ni
              read (16, *) (v1(i, j), u1(i, j), j = 1, nj)
            enddo
          end if

          if(igrav.ge.1) then
            u1=u1*1.2
            v1=v1*1.2
          end if

          read(16,'(a150)',end=161,err=165) text
          READ(text,*) icur_date
          go to 162
165       write(*,*) 'Current Input Index Error'
          stop
161       close(16)
          icur=3
          go to 163
162       continue
          backspace(unit=16)
163       continue

428       continue
!
          if(imod.eq.2) then
            u1=-u1
            v1=-v1
          end if
          if(imod.eq.1) v1=-v1
          if(imod.eq.3) u1=-u1
!
          do i=1,ni
            if(abs(d1(i,1)-d1(i,2)).lt..01) then
              if(abs(v1(i,1)-v1(i,2)).gt..01) v1(i,1)=v1(i,2)
            end if
            if(abs(d1(i,nj)-d1(i,nj-1)).lt..01) then
              if(abs(v1(i,nj)-v1(i,nj-1)).gt..01) v1(i,nj)=v1(i,nj-1)
            end if
          end do
!
          if(hs0.lt..002) go to 303
! ----- wave asymmetry effect
          if(imod.eq.0.and.bfricmax.ge..01) then
            do i=1,ni
              do j=1,nj
                if(abs(u1(i,j)).lt.abs(v1(i,j))) then
                  u1(i,j)=-abs(v1(i,j))
                end if
              end do
            end do
          end if
303       continue
        endif
      else
        if(noptset.eq.3)then !Alex
          if(imod.eq.0)then
            do i=1,ni
              do j=1,nj
                u1(i,j) = uin(i,j)
                v1(i,j) = vin(i,j)
                if(abs(u1(i,j)).ge.900.) u1(i,j)=0.0
                if(abs(v1(i,j)).ge.900.) v1(i,j)=0.0
              enddo
            enddo
          elseif(imod.eq.1)then
            do i=1,ni
              do j=1,nj
                u1(i,j) = uin(nj-j+1,i)
                v1(i,j) = vin(nj-j+1,i)
                if(abs(u1(i,j)).ge.900.) u1(i,j)=0.0
                if(abs(v1(i,j)).ge.900.) v1(i,j)=0.0
              enddo
            enddo
          elseif(imod.eq.2)then
            do i=1,ni
              do j=1,nj
                u1(i,j) = uin(ni-i+1,nj-j+1)
                v1(i,j) = vin(ni-i+1,nj-j+1)
                if(abs(u1(i,j)).ge.900.) u1(i,j)=0.0
                if(abs(v1(i,j)).ge.900.) v1(i,j)=0.0
              enddo
            enddo
          elseif(imod.eq.3)then
            do i=1,ni
              do j=1,nj
                u1(i,j) = uin(j,ni-i+1)
                v1(i,j) = vin(j,ni-i+1)
                if(abs(u1(i,j)).ge.900.) u1(i,j)=0.0
                if(abs(v1(i,j)).ge.900.) v1(i,j)=0.0
              enddo
            enddo
          endif
        endif 
        if(imod.eq.2)then
          u1=-u1
          v1=-v1
        end if
      endif

!!!!! Wind
      
      if(iwind.ge.1) then
        if(iplane.eq.2) go to 429
        do k=1,iwind1
          read(27,'(a150)',err=258,end=258) text
          READ(text,*,err=258,end=258) iwind_date
          go to 259
258       continue
          rewind(unit=27)
          read(27,'(a150)') text1
          read(27,'(a150)') text
          READ(text,*) iwind_date
259       continue
          write(*,*) 'Wind Index =',iwind_date
          if(mod(imod,2).eq.0) then
            do j=1,nj
              read(27,*) (u10(i,j),v10(i,j), i=1,ni)
            end do
          else
            do i=1,ni
              read(27,*) (u10(i,j),v10(i,j), j=1,nj)
            end do
          end if
        end do
        read(27,'(a150)',err=268,end=268) text
        READ(text,*,err=268,end=268) iwind_date
        go to 269
268     continue
        rewind(unit=27)
        read(27,'(a150)') text1
        read(27,'(a150)') text
        READ(text,*) iwind_date
269     continue
        write(*,*) 'Wind Index =',iwind_date
        write(*,*) ' '
!
        if(imod.eq.0) then
          do j = nj, 1, -1
            read (27, *) (u10(i, j), v10(i, j), i = 1, ni)
          enddo
        end if
!
        if(imod.eq.2) then
          do j = 1, nj
            read (27, *) (u10(i, j), v10(i, j), i = ni, 1, -1)
          enddo
        end if
        if(imod.eq.1) then
          do i = ni, 1, -1
            read (27, *) (v10(i, j), u10(i, j), j = nj, 1, -1)
          enddo
        end if
        if(imod.eq.3) then
          do i = 1, ni
            read (27, *) (v10(i, j), u10(i, j), j = 1, nj)
          enddo
        end if
!
429     continue
!
        if(imod.eq.2) then
          u10=-u10
          v10=-v10
        end if
        if(imod.eq.1) v10=-v10
        if(imod.eq.3) u10=-u10
!
        ws=0.
        do j=1,nj
          do i=1,ni
            cc2=sqrt(u10(i,j)**2+v10(i,j)**2)
              if(cc2.gt.cc1) then
                cc1=cc2
                wd=atan2(v10(i,j),u10(i,j)+.0001)/rad
              end if
            ws=ws+cc2
          end do
        end do
        ws=ws/float(ni*nj*2)

        if(hs0.lt..005.and.ws.ge..1) then
          do j=1,nj
            d1(1,j)=-100.
          end do
          prd=wd*rad
          dbar=prd
        end if

        write(*,*) 'Max internal-grid wind (m/sec, deg)=',cc1,wd
        write(*,*) ' '
      end if

! --- Note: STWAVE depth (dep0 and d1)is the center of
!           a CMS-Wave cell - so CMS-Wave dep is one more
!           row and one more column than STWAVE.
!
      imax=ni+1
      jmax=nj+1
!
      do k=1,ijstruc2
        if(imod.eq.0) then
          i=istruc2(k)
        end if
        if(imod.eq.2) then
          i=imax-istruc2(k)
        end if
        if(imod.eq.1) then
          i=jstruc2(k)
        end if
        if(imod.eq.3) then
          i=imax-jstruc2(k)
        end if
        if(i.lt.ismall) ismall=i
      end do

      do k=1,ijstruc3
        if(igetfile4.eq.0) then
          if(k3(k).ge.6) dstruc33(k)=dstruc33(k)-tide0+tide
        end if

        if(imod.eq.0) then
          i=istruc3(k)
        end if
        if(imod.eq.2) then
          i=imax-istruc3(k)
        end if
        if(imod.eq.1) then
          i=jstruc3(k)
        end if
        if(imod.eq.3) then
          i=imax-jstruc3(k)
        end if
        if(i.lt.ismall) ismall=i
      end do

      tide0=tide

      do k=1,ijstruc4
        if(imod.eq.0) then
          i=istruc4(k)
          j=jstruc4(k)
          ii=istruc4(k)+1
          if(ii.ge.ni) ii=ni
          jj=jstruc4(k)
        end if
        if(imod.eq.2) then
          i=imax-istruc4(k)
          j=jmax-jstruc4(k)
          ii=istruc4(k)-1
          if(ii.lt.1) ii=1
          jj=jstruc4(k)
        end if
        if(imod.eq.1) then
          i=jstruc4(k)
          j=jmax-istruc4(k)
          ii=istruc4(k)
          jj=jstruc4(k)+1
          if(jj.gt.ni) jj=ni
        end if
        if(imod.eq.3) then
          i=imax-jstruc4(k)
          j=istruc4(k)
          ii=istruc4(k)
          jj=jstruc4(k)-1
          if(jj.lt.1) jj=1
        end if
        dstruc4(k)=tide+eta(istruc4(k),jstruc4(k))-dstruc44(k)
        cc=(eta(istruc4(k),jstruc4(k))-eta(ii,jj))/dvarx(i)
        if(cc.gt.0.0) kstruc4(k)=1
        d1(i,j)=dstruc4(k)
!
        isame=0
        do k2=1,ijstruc2
          if(imod.eq.0) then
            i2=istruc2(k2)
            j2=jstruc2(k2)
          end if
          if(imod.eq.2) then
            i2=imax-istruc2(k2)
            j2=jmax-jstruc2(k2)
          end if
          if(imod.eq.1) then
            i2=jstruc2(k2)
            j2=jmax-istruc2(k2)
          end if
          if(imod.eq.3) then
            i2=imax-jstruc2(k2)
            j2=istruc2(k2)
          end if
          if(i2.eq.i.and.j2.eq.j) then
            isame=1
            exit
          end if
        end do
!
! --- the following cannot be true if it is also a runup cell
        if(isame.eq.0.and.d1(i,j).lt..01) d1(i,j)=.01
        if(i.lt.ismall) ismall=i
      end do
      ismall=ismall-1

!!!!!!
      
      ix1=0
      ix2=0
      do 290 i=1,ni-1
        dsum=0.
        dsum1=0.
        do j=1,nj
          dsum=dsum+d1(i,j)+d1(i+1,j)
          if(i+2.le.ni) dsum1=dsum1+d1(i,j)+d1(i+1,j)+d1(i+2,j)
        end do
        dsum=dsum/float(nj)
        dsum1=dsum1/float(nj)
        do j=1,nj
          if(i+1.ge.ismall) then
            ix1(i)=1
            ix2(i)=1
            go to 291
          end if
          if(abs(d1(i,j)+d1(i+1,j)-dsum).gt..00001) then
            ix1(i)=1
            go to 291
          end if
        end do

        do j=1,nj
          if(i+2.le.ni) then
            if(abs(d1(i,j)+d1(i+1,j)+d1(i+2,j)-dsum1).gt..00001) then
              ix2(i)=1
              go to 290
            end if
          end if
        end do
290   continue
291   continue
      do i=2,ni
        if(ix1(i-1).eq.1) ix1(i)=1
        if(ix2(i-1).eq.1) ix2(i)=1
      end do

      dep(1,1)=d1(1,1)
      dep(1,jmax)=d1(1,nj)
      keep=0
      do j=2,nj
        dep(1,j)=(d1(1,j-1)+d1(1,j))/2.
        if(dep(1,j).lt..01) dep(1,j)=.01
        if(dep(1,j).gt..01) keep=1
      ENDDO

      do i=2,ni
        dep(i,1)=(d1(i-1,1)+d1(i,1))/2.
        dep(i,jmax)=(d1(i-1,nj)+d1(i,nj))/2.
        do j=2,nj
          if(d1(i,j).lt.deps(i)) deps(i)=d1(i,j)
          dep(i,j)=(d1(i-1,j-1)+d1(i,j-1)+d1(i-1,j)+d1(i,j))/4.
          if(dep(i,j).lt..01) dep(i,j)=.01
        ENDDO
        if(d1(i,1).lt.deps(i)) deps(i)=d1(i,1)
      ENDDO
      deps(1)=deps(2)
!
      if(keep.eq.0) then
        dep(2,1)=d1(2,1)
        dep(2,jmax)=d1(2,nj)
        do j=2,nj
          dep(2,j)=(d1(2,j-1)+d1(2,j))/2.
          if(dep(2,j).lt..01) dep(2,j)=.01
        ENDDO
      end if
!
      dep(imax,1)=d1(ni,1)
      dep(imax,jmax)=d1(ni,nj)
      do j=2,nj
        dep(imax,j)=(d1(ni,j-1)+d1(ni,j))/2.
        if(dep(imax,j).lt..01) dep(imax,j)=.01
      ENDDO
!
      nblock=0
      inquire(file='block.dat',exist=getfile3)
      if(getfile3) then
        nblock=1
        write(*,*) ' *** block.dat FILE FOUND ***'
        write(*,*) '     Read block cell file'
        write(*,*) ' '
        open(unit=19,file='block.dat',status='old')
        read(19,*) kblock
        do k=1,kblock
          read(19,*) iblock,jblock
          if(iblock.ne.1) d1(iblock,jblock)=0.
        ENDDO
        close(19)
      end if
!
      do i=1,ni
        ijb(i,1)=2
        ijb(i,nj)=2
      enddo
!
      ni1=ni-1
      nj1=nj-1
      do j=1,nj
        do i=1,ni
          if(d1(i,j).le.0.) then
            ijb(i,j)=0
          end if
        enddo
      enddo
!
!
      do j=2,nj1
        do i=2,ni
          j1=j+1
          ib=i-1
          jb=j-1
          if(d1(i,j).gt.0..and.d1(i,jb).le.0.) ijb(i,j)=4
          if(d1(i,j).gt.0..and.d1(i,j1).le.0.) ijb(i,j)=5
          if(d1(i,j).gt.0..and.d1(ib,j).le.0.) then
            ijb(i,j)=8
            if(d1(i,jb).le.0.) ijb(i,j)=6
            if(d1(i,j1).le.0.) ijb(i,j)=7
          end if
        enddo
      enddo
!
      do i=2,ni
        ib=i-1
        if(d1(i,1).gt.0..and.d1(ib,1).le.0.) ijb(i,1)=9
        if(d1(i,nj).gt.0..and.d1(ib,nj).le.0.) ijb(i,nj)=9
      enddo

      if(nblock.eq.1) then
        open(unit=19,file='block.dat',status='old')
        read(19,*) kblock
        do k=1,kblock
          read(19,*) iblock,jblock
          ijb(iblock,jblock)=-1
          if(ijb(iblock,jblock+1).eq.4) ijb(iblock,jblock+1)=2
          if(ijb(iblock,jblock-1).eq.5) ijb(iblock,jblock-1)=2
        enddo
        close(19)
      end if
!
      ikap=1
      if(depmin0.gt.1.) then
        do i=2,ni
          do j=2,nj1
            if(iabs(ijb(i,j)).ne.1.and.d1(i,j).lt.5.) go to 299
          end do
          ikap=i
        end do
299     continue
      end if

!!!!! Forward reflection      
      
      krmx=0
      do j=2,nj1
        do i=2,ni1
          if(ijb(i,j).eq.4.or.ijb(i,j).eq.5) then
            i1=i+1
            j1=j+1
            ib=i-1
            jb=j-1
            krmx=krmx+1
            if(krmx.gt.6*ipmx) go to 5000                      !Changed to 6*ipmx to duplicate the inline limits - MEB 02/02/2020
            kr(1,krmx)=i
            kr(2,krmx)=j
            yangl(krmx)=0.
            if(ijb(ib,jb).ge.4.and.ijb(i1,jb).lt.4) yangl(krmx)=45.
            if(ijb(ib,j1).ge.4.and.ijb(i1,j1).lt.4) yangl(krmx)=-45.
            if(ijb(i1,jb).ge.4.and.ijb(ib,jb).lt.4) yangl(krmx)=-45.
            if(ijb(i1,j1).ge.4.and.ijb(ib,j1).lt.4) yangl(krmx)=45.
            if(ijb(ib,j).ge.4.and.ijb(i1,j1).ge.4) yangl(krmx)=30.
            if(ijb(ib,j).ge.4.and.ijb(i1,jb).ge.4) yangl(krmx)=-30.
            if(ijb(i1,j).ge.4.and.ijb(ib,j1).ge.4) yangl(krmx)=-30.
            if(ijb(i1,j).ge.4.and.ijb(ib,jb).ge.4) yangl(krmx)=30.
            if(iark.eq.0 .or. iark.eq.1) then
              rk(krmx)=ark
            else
              rk(krmx)=reflty(i,j)
              if(rk(krmx).gt.1.) rk(krmx)=ark
            end if
          end if
        enddo
      enddo
5000  continue

      if(iark.ge.1) then
        write(*,*) 'Total forward reflection cells  =',krmx
        if(krmx.gt.6*ipmx) then                                    !Changed to 6*ipmx to duplicate the inline limits - MEB 02/02/2020
          write(*,*) 'Total forward reflection cells >',6*ipmx
          write(*,*) 'Need to increase total reflection cells'
          write(*,*) 'or turn off reflection calc -- Run Terminated'
          stop
        end if
        write(*,*) ' '
      end if

!!!!! Backward reflection 
      
      krf=0
      do j=1,nj
        do i=2,ni1
          i1=i+1
          if(ijb(i,j).ne.0.and.ijb(i1,j).eq.0) then
            j1=j+1
            if(j1.gt.nj) j1=nj
            ib=i-1
            jb=j-1
            if(jb.lt.1) jb=1
            krf=krf+1
            if(krf.gt.6*ipmx) go to 5001                        !Changed to 6*ipmx to duplicate the inline limits - MEB 02/02/2020
            kcr(1,krf)=i
            kcr(2,krf)=j
            xangl(krf)=0.
            if(ijb(i1,j1).eq.0.and.ijb(i1,jb).ne.0) xangl(krf)=30.
            if(ijb(i1,jb).eq.0.and.ijb(i1,j1).ne.0) xangl(krf)=-30.
            if(ijb(i,j1).eq.0.and.ijb(i,jb).ne.0) then
              xangl(krf)=45.
              if(ijb(i1,jb).eq.0) xangl(krf)=30.
            end if
            if(ijb(i,jb).eq.0.and.ijb(i,j1).ne.0) then
              xangl(krf)=-45.
              if(ijb(i1,j1).eq.0) xangl(krf)=-30.
            end if
            if(iarkr.eq.0 .or. iarkr.eq.1) then
              rkr(krf)=arkr
            else
              rkr(krf)=refltx(i,j)
              if(rkr(krf).gt.1.) rkr(krf)=arkr
            end if
            if(j.eq.1) then
              xangl(krf)=0.
              if(ijb(i,j1) .eq.0) xangl(krf)=-30.
              if(ijb(ib,j1).eq.0) xangl(krf)=-45.
            end if
            if(j.eq.nj) then
              xangl(krf)=0.
              if(ijb(i,jb) .eq.0) xangl(krf)=30.
              if(ijb(ib,jb).eq.0) xangl(krf)=45.
            end if
          end if
        enddo
      enddo
5001  continue
!
      if(iarkr.ge.1) then
        write(*,*) 'Total backward reflection cells =',krf
        if(krf.gt.6*ipmx) then                                     !Changed to 6*ipmx to duplicate the inline limits - MEB 02/02/2020
          write(*,*) 'Total backward reflection cells >',6*ipmx
          write(*,*) 'Need to increase total reflection cells'
          write(*,*) 'allowed -- Stop run'
!         call PressReturn('ERROR:')
          stop
        end if
        write(*,*) ' '
      end if
      
      if(noptset.ne.3)then !Alex, Non inline model
        if(icur.ge.1) then
          u(1,1)=u1(1,1)
          u(1,jmax)=u1(1,nj)
          v(1,1)=v1(1,1)
          v(1,jmax)=v1(1,nj)
          do j=2,nj
            u(1,j)=(u1(1,j-1)+u1(1,j))/2.
            v(1,j)=(v1(1,j-1)+v1(1,j))/2.
          enddo

          do i=2,ni
            u(i,1)=(u1(i-1,1)+u1(i,1))/2.
            u(i,jmax)=(u1(i-1,nj)+u1(i,nj))/2.
            v(i,1)=(v1(i-1,1)+v1(i,1))/2.
            v(i,jmax)=(v1(i-1,nj)+v1(i,nj))/2.
            do j=2,nj
              u(i,j)=(u1(i-1,j-1)+u1(i,j-1)+u1(i-1,j)+u1(i,j))/4.
              v(i,j)=(v1(i-1,j-1)+v1(i,j-1)+v1(i-1,j)+v1(i,j))/4.
            enddo
          enddo

          u(imax,1)=u1(ni,1)
          u(imax,jmax)=u1(ni,nj)
          v(imax,1)=v1(ni,1)
          v(imax,jmax)=v1(ni,nj)
          do j=2,nj
            u(imax,j)=(u1(ni,j-1)+u1(ni,j))/2.
            v(imax,j)=(v1(ni,j-1)+v1(ni,j))/2.
          enddo
        end if
      endif  

      if(hs0.lt..0001) then
        mwd=md2
        if(ws.ge..1) then
          mwd=(hpai+wd*rad)/dth
          if(mwd.gt.md) mwd=md
          if(mwd.lt.imd) mwd=imd
        end if
        dsfd(nf,mwd)=.0001
        fp=fcn(nf)
      end if
!
      big=0.
      ibig=nf
      do nn=1,nf
        sum=0.
        do mm=imd,md
          sum=sum+dsfd(nn,mm)
        enddo
        fsp(nn)=sum
        if(sum.gt.big) then
          ibig=nn
          big=sum
        end if
      enddo
!
      write(*,*) 'Total frequency bin =',nf
      write(*,*) 'exist frequency bin =',(fcn(nn),nn=1,nf)
      write(*,*) 'Normalized 1-D spec =',(fsp(nn)/big,nn=1,nf)
      write(*,*) 'Peak Index=',ibig,fp,fcn(ibig)
      
      if(noptset.eq.3 .and. nsteer.eq.1)then
        write(9,*) 'Total frequency bin =',nf
        write(9,*) 'exist frequency bin =',(fcn(nn),nn=1,nf)
        write(9,*) 'Normalized 1-D spec =',(fsp(nn)/big,nn=1,nf)
        write(9,*) 'Peak Index=',ibig,fp,fcn(ibig)
      endif        

	  if(fp.eq.0.) fp=fcn(ibig)

      if(ws.gt.10.) iwvbk=2
      if(iwind.ge.1) iwvbk=2
!
!     Switch to BJ breaking for larger wind or read wind field input ***
!
      wd=wd*rad
      if(ws.ge..1) then
        a1=0.0002
        cc=a1*g*8./pai
        ccc=a1*g*1.333/pai
        ph0=g/ws*.9
        if(ph0.gt.pai2) ph0=pai2
        cc1=wd+HPAI
        if(iview.ge.1) cc1=cc1+DTH/2.
        do mm=imd,md
          dd=(DTH*MM-cc1)/2.
          wdd(mm)=ccc*cos(dd)**4
        enddo

        if(hs0.lt..05) then
          ph0hz=ph0/pai2*3.
          if(ph0hz.gt..33) ph0hz=.33
          ! --- limit ph0hz to .33 avoid initial large height for small wind.
          do nn=1,nf
            om=pai2*fcn(nn)
            om51=om**5*exp(0.74*min((ph0hz/fcn(nn))**4,10.))  &
               /dth/(fcn(nf)-fcn(nf-1))
            do mm=imd,md
              dd=-hpai-dth/2.+dth*mm-wd/2.
              if(abs(dd).lt.hpai) then
                dsfd(nn,mm)=wdd(mm)/om51*g*2000.
              end if
! --- wind wave calibration use above line and two ph0 calculation.
!             if(fcn(nn).lt..2) dsfd(nn,mm)=   &
!               dsfd(nn,mm)/exp((.2-fcn(nn))*100.)
! --- limit fcn to .2 avoid large period in upwind lake area.
            end do
          end do
        end if
!
        go to 661
      end if
!
      sum=0.
      do nn=1,ibig-2
!       if(ibig-nn.le.2) go to 61
        if(fsp(nn)/big.gt..05) exit  !go to 61
        cc=df(nn)/df(ibig)
        sum=sum+fsp(nn)*cc
        do mm=imd,md
          dsfd(ibig,mm)=dsfd(ibig,mm)+dsfd(nn,mm)*cc
        enddo
      enddo
61    ib=nn
      big=big+sum
      if(ib.eq.ibig) ib=ibig-1
      if(ib.le.0) ib=1
!
      sum=0.
      do nn=nf,ibig+2,-1
!       if(nn-ibig.le.2) go to 63
        if(fsp(nn)/big.gt..05) exit !go to 63
        cc=df(nn)/df(ibig)
        sum=sum+fsp(nn)*cc
        do mm=imd,md
          dsfd(ibig,mm)=dsfd(ibig,mm)+dsfd(nn,mm)*cc
        enddo
      enddo
63    ie=nn
      big=big+sum
      if(ie.eq.ibig) ie=ibig+1
      if(ie.gt.nf) ie=nf
      fsp(ibig)=big
!
      ib0=ib-1

      if(ie-ib0.lt.nf) then
        write(*,*) ' '
        write(*,*) 'Renew frequency bin =',ie-ib0
        write(*,*) 'Write frequency bin =',(fcn(nn),nn=ib,ie)
        write(*,*) 'Normalized 1-D spec =',(fsp(nn)/big,nn=ib,ie)
        write(*,*) 'New Peak Index=',ibig-ib0,fcn(ibig)
        if(ibnd.ge.1) fp=fcn(ibig)
        sum=0.
        do mm=imd,md
          do nn=ib,ie
            sum=sum+dsfd(nn,mm)*df(nn)
          enddo
        enddo
        hss=4.*sqrt(sum*dth)
        write(*,*) 'Truncated 4*sqrt(E) =',hss,'m'
      else
        write(*,*) 'NO CHANGE OF INPUT SPECTRUM'
        go to 661
      end if

      nf=ie-ib0
      ibig=ibig-ib0

      if(ib.eq.1) go to 661
      do nn=ib,ie
        ii=nn-ib0
        fcn(ii)=fcn(nn)
        df(ii)=df(nn)
        fsp(ii)=fsp(nn)
        do mm=imd,md
          dsfd(ii,mm)=dsfd(nn,mm)
        enddo
      enddo
  661 continue
!
      if(nf.ge.18.and.df(2).le.02) then
        ib=mod(ibig,3)
        if(ib.eq.0) then
          ib=3
          do mm=imd,md
            dsfd(2,mm)=dsfd(1,mm)+dsfd(2,mm)
          enddo
          fsp(2)=fsp(1)+fsp(2)
        end if
!
        ii=0
        do nn=ib,nf,3
          n1=nn+1
          if(nn.eq.1) then
            ii=ii+1
            fsp(ii)=(fsp(nn)+fsp(n1))/3.
            do mm=imd,md
              dsfd(ii,mm)=(dsfd(nn,mm)+dsfd(n1,mm))/3.
            enddo
            cycle  !go to 88
          end if
!
          c2=(fcn(n1)-fcn(nn-1))/2.
          if(nn-2.eq.0) then
            c1=c2
          else
            c1=(fcn(nn)-fcn(nn-2))/2.
          end if
!
          if(nn.eq.nf) then
            ii=ii+1
            fsp(ii)=(fsp(nf-1)+fsp(nf))/3.
            fcn(ii)=fcn(nn)
            do mm=imd,md
              dsfd(ii,mm)=(dsfd(nf-1,mm)+dsfd(nf,mm))/3.
            enddo
            cycle  !go to 88
          end if
!
          if(n1+1.gt.nf) then
            c3=c2
          else
            c3=(fcn(n1+1)-fcn(nn))/2.
          end if
!
          cc=c1+c2+c3
! 
          ii=ii+1
          fcn(ii)=fcn(nn)
! 
          fsp(ii)=(fsp(nn-1)*c1+fsp(nn)*c2+fsp(n1)*c3)/cc
!         write(*,*) 'ii,fsp=',ii,fsp(ii),fsp(nn-1),fsp(nn),fsp(n1)
!         write(*,*) 'ii,cc =',ii,cc,c1,c2,c3
          if(nn.eq.ibig) then
            big=fsp(ii)
            ibig=ii
          end if
!
          do mm=imd,md
            dsfd(ii,mm)=(dsfd(nn-1,mm)*c1+dsfd(nn,mm)*c2+dsfd(n1,mm)*c3)/cc
          enddo
        enddo
!
        if(n1.eq.nf-1) then
          fsp(ii)=fsp(ii)+fsp(nf)*c3/cc
          do mm=imd,md
            dsfd(ii,mm)=dsfd(ii,mm)+dsfd(nf,mm)*c3/cc
          enddo
        end if
!
        nf=ii
        if(ibig.gt.nf) ibig=nf
        write(*,*) ' '
        write(*,*) 'New group total frequency bin =',nf
        write(*,*) 'New group fixed frequency bin =',(fcn(nn),nn=1,nf)
        write(*,*) 'Normalized 1-D spec =',(fsp(nn)/big,nn=1,nf)
        write(*,*) 'New Peak Index=',ibig,fcn(ibig)
        if(ibnd.ge.1) fp=fcn(ibig)
      end if
!
      if(ws .ge. 0.1) then
        do mm=imd,md
          dsfd(nf+1,mm)=dsfd(nf,mm)/2.
        enddo
        fcn(nf+1)=2.*fcn(nf)-fcn(nf-1)
        nf=nf+1
        fsp(nf)=fsp(nf-1)/2.
!
        if(ws .ge. 25.) then
          f0=2.*fcn(1)-fcn(2)
          if(f0.gt..001.and.fcn(1).gt..1) then
            do nn=nf,1,-1
              do mm=imd,md
                dsfd(nn+1,mm)=dsfd(nn,mm)
              enddo
              fcn(nn+1)=fcn(nn)
              fsp(nn+1)=fsp(nn)
            enddo
            fcn(1)=f0
            do mm=imd,md
              dsfd(1,mm)=0.
            enddo
            nf=nf+1
            ibig=ibig+1
          end if
        end if
!
        write(*,*) ' '
        write(*,*) 'Final total frequency bin =',nf
        write(*,*) 'Final fixed frequency bin =',(fcn(nn),nn=1,nf)
        write(*,*) 'Normalized 1-D spec =',(fsp(nn)/big,nn=1,nf)
        write(*,*) 'New Peak Index=',ibig,fcn(ibig)
        if(ibnd.ge.1) fp=fcn(ibig)
      end if
!
! -- Convert STWAVE spectrum (density) to total energy
!    in each cell
!
! ** frequency interval
      fspsum=fsp(1)
      DO NN=2,NF-1
        fspsum=fspsum+fsp(nn)
        DF(NN)=0.5*(FCN(NN+1)-FCN(NN-1))
      enddo
      fspsum=fspsum+fsp(nf)
      DF(1)=DF(2)
      DF(NF)=DF(NF-1)

      sum=0.
      sumx=0.
      sumy=0.
      do mm=imd,md
        do nn=1,nf
          dsfd(nn,mm)=dsfd(nn,mm)*dth*df(nn)
          sum=sum+dsfd(nn,mm)
          sumx=sumx+dsfd(nn,mm)*cosa(mm)
          sumy=sumy+dsfd(nn,mm)*sina(mm)
        enddo
      enddo
      wdir=atan2(sumy,sumx)

      hs0=4.*sqrt(sum)
      write(*,*) ' '
      write(*,*) 'Check 4*sqrt(E) =',hs0,'m'
      if(hs1.gt..01.and.ibnd.eq.0) then
        if(ws.ge.5. .and. iwnd.le.1) then
          if(ws.ge.50.) ws=50.
          hsmax=ws**1.5/sqrt(fp*g)/17.8*cos(wd-wdir)**4*abs(sin(wd))         ! ----- 17.8 = (4.3)^1.5*2 where 2 is a correction factor (2.5 in reality?)
          if(hsmax.ge.hs0*2.) hsmax=hs0*2.
          if(hsmax.gt.hs0) then
            write(*,*) 'Max Windwave(m) =',hsmax,'m'
            write(*,*) '*** Inflat incident wave ***'
            dsfd=dsfd*(hsmax/hs0)**2
            hs0=hsmax
            d1max=d1(1,1)
            do jj=1,nj
              if(d1(1,jj).gt.d1max) d1max=d1(1,jj)
            end do
            do jj=1,nj
              wcc(jj)=wcc(jj)*(abs(d1(1,jj))/d1max)**0.25
            end do
          end if
        end if
      end if
!
      flat=1.
      if(ws.lt..01.and.fspsum.gt..01) then
        flat=1.+exp(-10.*(1.-fsp(ibig)/fspsum))
        write(*,*)
        write(*,*) 'Effect transition between Hs (H_mo) & Hrms (H_mono)'
        write(*,*) 'Inflation of wave energy (in calculation) by', flat
        write(*,*)
      end if
      cflat=sqrt(16./flat)
!
      dsfd=dsfd*flat
!
      HSG=HS0
!
      write(*,*) 'Wsmag =', wsmag
!
      WRITE(*,*) ' '
      if(iplane.le.1.or.iwave.eq.1) then
        WRITE(*,*) '*** Completion of Reading Spectrum ***'
        WRITE(*,*) ' '
      endif
!
      eta=0.
      TP=1/fp
      TS=TP/1.05
      write(*,*) 'TS(=TP/1.05)',TS
      write(*,*) ' '
      WL0=1.56*TS*TS
!
      if(iplane.eq.2.and.hs0.lt..001) then
        if(abs(wd).ge.2.1.or.ws.lt..1) then
          h13=0.
          t13=0.
          dmn=0.
          go to 80
        end if
      end if
!
      if(hs0.ge..001.and.imod.eq.0) then
! ----- wave asymmetry effect
        if(icur.ge.1.and.bfricmax.lt..005) then
          hs0fp2=hs0*fp**2
          do i=1,ni-1
            do j=1,nj
              slope=(d1(i,j)-d1(i+1,j))/dvarx(i)
              if(slope.lt..02) slope=.02
              if(hs0fp2/slope**2.lt.20.) then
                if(abs(u1(i,j)).lt.abs(v1(i,j)))then
                  u1(i,j)=-abs(v1(i,j))
                  u(i,j)=-abs(v(i,j))
                end if
              end if
            end do
          end do
        end if
      endif
!
      small=1.
      kpeak=1
      do nn=1,nf
        sdiff=abs(fcn(nn)-fp)
        if(sdiff.lt.small) then
          small=sdiff
          kpeak=nn
        end if
      enddo
!---------------------------------------------------------------
!  input data for each region 
!    6.-10. are data in subroutine dspec in subroutine input1
!   11. yes or no of output, detail region and wave reflection
!   12. mesh number
!   13. depth (just declare), 14. depth (from file)
!   15. cell  (just declare), 16. cell  (from file)
!   16. current data from file
!   20.-22. if ick3=1, positions of reflection in y direction
!   23.-25. if ick4=1, positions of reflection in x direction
!---------------------------------------------------------------
!
      IGMX=IMAX-1
      JGMX=JMAX-1
!
      ICK3=0
      ICK4=0
      if(iark.ge.1) ICK3=1
      if(iarkr.ge.1) ICK4=1
!
      DO JJ=1,JGMX
        DMNJ(JJ)=DBAR
      enddo
!
!-----------------------------------------------------
!  initial integrated energy density for each region
!-----------------------------------------------------
      CALL INITL_inline
!---------------------------------------------------
!  calculation of wave action balance equation
!  (main of calculation including subroutine calc) 
!---------------------------------------------------
      CALL SETIN_inline(IKAP)
!--------------------------------
!  output of calculated results
!--------------------------------
80    continue
      if (ibreak.ge.2) call dfds_inline
      if (irs.ge.1.and.irs0.eq.0) call rstress_inline
90    continue

      
      CALL OUTFILE_inline

!----------------------
!  end of calculation
!----------------------
 1000 CONTINUE
!
      if(iplane.eq.1) then
        if(iwet.eq.-1.or.iwet.eq.-3) then
          exx=sw13
          eyy=sa13
          do j=1,jgmx
            do i=1,igmx
              sw13(i,j)=exx(imax-i,jmax-j)
              sa13(i,j)=eyy(imax-i,jmax-j)
            end do
          end do
          exx=tw13
          eyy=ta13
          do j=1,jgmx
            do i=1,igmx
              tw13(i,j)=exx(imax-i,jmax-j)
              ta13(i,j)=eyy(imax-i,jmax-j)
            end do
          end do
          exx=dw13
          eyy=da13
          do j=1,jgmx
            do i=1,igmx
              dw13(i,j)=exx(imax-i,jmax-j)
              da13(i,j)=eyy(imax-i,jmax-j)
            end do
          end do
        end if
        go to 770
      end if 
!
#ifdef XMDF_IO
      if (ixmdf.ge.1) then
        CALL XMDFOUT_inline
      endif
#endif

      if(noptset.eq.3) goto 4100   !Wu/Zhang/Alex
      go to 70

410   continue
      print *, 'end of file encountered reading boundary spectrum'

420   continue
      if(iwave.eq.1) close(24)
!     call PressReturn('SUCCESS')

!      DEALLOCATE (dep0)                 !Wu/Zhang
!      DEALLOCATE (fsp,xc,yc,wc,wcd)     !Wu
!      DEALLOCATE (refltx,reflty)        !Wu/Zhang

4100  continue           !Wu/Zhang 

contains
#ifdef XMDF_IO
!********************************************************************************
      subroutine xmdfout_inline
!********************************************************************************      
      use comvarbl, only: reftime
      use xmdf
      USE GLOBAL_INLINE, ONLY: NPF,MPD,MPD2,IPMX,JPMX,KOMX,NOMX,IGPX,JGPX
      use cms_def, only: noptset,noptzb,nsteer,dtsteer           !Alex
      use comvarbl, only: ctime
      REAL, ALLOCATABLE :: height(:),period(:),dir(:),brkdiss(:),radstr(:)  !Alex 
      REAL, ALLOCATABLE :: wave(:),depth(:),surge(:),currents(:)            !Alex      
      REAL, ALLOCATABLE :: surgew(:)            !Alex      
    
      common /fl2wav/depin(ipmx,jpmx),etain(ipmx,jpmx), &        !Alex
             uin(ipmx,jpmx),vin(ipmx,jpmx)                       !Alex
      COMMON /VPAI/PAI2,PAI,HPAI,RAD,akap,imod,iprp,island,imd,iprpp     &
                   ,nonln,igrav,isolv,ixmdf,iproc,imud,iwnd,depmin0
      COMMON /DATA/NF,MD,IMAX,JMAX,IGMX,JGMX,JCPB,JCPE,JCB,JCE,NFF,MDD
      COMMON /DATB/ICK3,ICK4,ICHOICE,HS0,wd,ws,ph0,wdd(mpd),aslop(ipmx)
      COMMON /DATC/TP,PRD,IBND,VL,TIDE,KPEAK,IBACK,NBLOCK,IWIND,WSMAG
      COMMON /DATD/DX,DY,DXX,DMESH,DTH,kdate,idate,depmax,depmax0
      COMMON /DVAR/DVARX(IGPX),DVARY(JGPX),ETA(IGPX,JGPX),HSK(JGPX)
      COMMON /SPECB/SOP(KOMX,NPF,MPD2),SR(2*IPMX,NPF,MPD)
      COMMON /OUTP/KOUT,IJSP(2,KOMX),IRS,IBREAK,ICUR,IWET,INST
      COMMON /OUTN/nest,inest(komx),jnest(komx),ix1(igpx),ix2(igpx)
      COMMON /IJBP/IJB(IGPX,JGPX),DEP(IPMX,JPMX),DEPS(IPMX),DBIG(IPMX)
      COMMON /WAVI/H13(IGPX,JGPX),T13(IGPX,JGPX),DMN(IGPX,JGPX)
      COMMON /WAVR/H13R(IGPX,JGPX),T13R(IGPX,JGPX),DMNR(IGPX,JGPX)
      COMMON /WAVS/H13S(IGPX,JGPX),IBR(IGPX,JGPX),DISS(IGPX,JGPX)
      COMMON /FDS/FCN(NPF),DCM(MPD),DSFD(NPF,MPD)
      COMMON /DFCN/DF(NPF),FFCN(NPF),DFINP(NPF),PL1E(NPF),PL1T(NPF)
      common /rsb/ wxrs(ipmx,jpmx),wyrs(ipmx,jpmx)
      common /rsm/ cosa(mpd),sina(mpd),d1(ipmx,jpmx)
      common /uv1/u1(ipmx,jpmx),v1(ipmx,jpmx)
      common /struc0/ijstruc1,ijstruc2,ijstruc3,ijstruc4,ismall
      common /struc1/istruc1(komx),jstruc1(komx),dstruc1(komx)
      common /struc2/istruc2(NOMX),jstruc2(NOMX),dstruc2(NOMX)
      common /struc3/istruc3(komx),jstruc3(komx),dstruc3(komx) &
                    ,k3(komx),dstruc33(komx)
      common /struc4/istruc4(komx),jstruc4(komx),dstruc4(komx) &
                    ,kstruc4(komx),k4(komx),dstruc44(komx)
      common /wavenum/ itms,ibf,iark,iarkr,bf,ark,arkr       !Wu
      common /origin/x0,y0,azimuth,isteer,iidate,sinaz,cosaz
      common /FileNames/ OptsFile, DepFile, CurrFile, EngInFile,    &
                         WaveFile, ObsFile, EngOutFile, NestFile,   &
                         BreakFile, RadsFile, StrucFile, SurgeFile, &
                         MudFile, FricFile, FrflFile, BrflFile,     &
                         SpecFile, WindFile, XMDFFile, SetupFile,   &  !Mitch 3/22/2017
                         SeaFile, SwellFile, ShipFile                  !Mitch 3/22/2017
      logical*2 :: FOUND=.FALSE.
      CHARACTER*180 WaveFile,ObsFile,EngOutFile,BreakFile,RadsFile
      CHARACTER*180 OptsFile,DepFile,CurrFile,EngInFile
      CHARACTER*180 NestFile,StrucFile,SurgeFile
      CHARACTER*180 MudFile,FricFile,FrflFile,BrflFile,WindFile
      CHARACTER*180 SpecFile, XMDFFile,SetupFile,ShipFile
      CHARACTER*180 SeaFile, SwellFile                                 !Mitch 3/22/2017
      CHARACTER*20  PREFIX
      REAL*8 TIME2
      INTEGER ERROR, PID, DID, DGID, COUNT
      real :: cosaz,sinaz
      
      ALLOCATE (height(IGMX*JGMX),period(IGMX*JGMX),dir(IGMX*JGMX))
      ALLOCATE (brkdiss(IGMX*JGMX),radstr(2*IGMX*JGMX),wave(2*IGMX*JGMX))
      ALLOCATE (depth(IGMX*JGMX),surge(IGMX*JGMX),currents(2*IGMX*JGMX))
      ALLOCATE (surgew(IGMX*JGMX))

      cosaz = cos(AZIMUTH*rad)
      sinaz = sin(AZIMUTH*rad)
      
      K=1; L=1   
      if (imod .eq. 0 .or. imod .eq. 2) then                  
        DO J=1,jgmx         
          DO I=1,igmx 
            depth(K)=depin(i,j)
            surge(K)=etain(i,j)
            currents(L)  =uin(i,j)*cosaz-vin(i,j)*sinaz
            currents(L+1)=uin(i,j)*sinaz+vin(i,j)*cosaz
            K=K+1
            L=L+2
          ENDDO              
        ENDDO
      else
        DO J=1,igmx         
          DO I=1,jgmx 
            depth(K)=depin(i,j)
            surge(K)=etain(i,j)
            currents(L)  =uin(i,j)*cosaz-vin(i,j)*sinaz
            currents(L+1)=uin(i,j)*sinaz+vin(i,j)*cosaz
            K=K+1
            L=L+2
          ENDDO              
        ENDDO
      endif
      
      NIJ=IGMX*JGMX
      K=1
      L=1
      cosazz=cosaz
      sinazz=sinaz
      if (mod(imod,2).eq.0) then
        if (imod .eq. 0) then
          ibeg = 1
          iend = igmx
          iinc = 1
          jbeg = 1
          jend = jgmx
          jinc = 1
        else
          ibeg = igmx
          iend = 1
          iinc = -1
          jbeg = jgmx
          jend = 1
          jinc = -1
          cosazz=-cosaz
          sinazz=-sinaz
        endif                              
        DO J=jbeg,jend,jinc 
          DO I=ibeg,iend,iinc
            HEIGHT(K) = H13S(I,J)
            PERIOD(K) = T13(I,J)
            DIR(K) = DMN(I,J)+AZIMUTH
            IF (DIR(K).LT.0) DIR(K)=DIR(K)+360.
            WAVE(L)=HEIGHT(K)*COS(DIR(K)*RAD)
            WAVE(L+1)=HEIGHT(K)*SIN(DIR(K)*RAD)
            IF (irs .gt. 0) THEN
              radstr(L)  =WXRS(I,J)*cosazz-WYRS(I,J)*sinazz
              radstr(L+1)=WXRS(I,J)*sinazz+WYRS(I,J)*cosazz
            ENDIF
            IF (IBREAK.EQ.1) THEN
              BRKDISS(K)=IBR(I,J)
            ELSEIF (IBREAK.EQ.2) THEN
              BRKDISS(K)=DISS(I,J)
            ENDIF
            K=K+1
            L=L+2
          ENDDO              
        ENDDO
      else
        if (imod .eq. 1) then
          ibeg = 1
          iend = igmx
          iinc = 1
          jbeg = jgmx
          jend = 1
          jinc = -1
          cosazz=-sinaz
          sinazz=cosaz
        else
          ibeg = igmx
          iend = 1
          iinc = -1
          jbeg = 1
          jend = jgmx
          jinc = 1
          cosazz=sinaz
          sinazz=-cosaz
        endif

        DO I=ibeg, iend, iinc
          DO J=jbeg, jend, jinc
            surgew(K)=eta(i,j)
            HEIGHT(K) = H13S(I,J)
            PERIOD(K) = T13(I,J)
            DIR(K) = DMN(I,J)+AZIMUTH
            IF (DIR(K).LT.0) DIR(K)=DIR(K)+360.
            WAVE(L)=HEIGHT(K)*COS(DIR(K)*RAD)
            WAVE(L+1)=HEIGHT(K)*SIN(DIR(K)*RAD)
            IF (irs .gt. 0) THEN
              radstr(L)  =WXRS(I,J)*cosazz-WYRS(I,J)*sinazz
              radstr(L+1)=WXRS(I,J)*sinazz+WYRS(I,J)*cosazz
            ENDIF
            IF (IBREAK.EQ.1) THEN
              BRKDISS(K)=IBR(I,J)
            ELSEIF (IBREAK.EQ.2) THEN
              BRKDISS(K)=DISS(I,J)
            ENDIF
            depth(K)=d1(i,j)
            K=K+1
            L=L+2
          ENDDO              
        ENDDO
      endif
!      
      if(noptset.ne.3)then !Alex, not inline steering
        TIME2=float(IDATE)
      else !if(nsteer.eq.1)then   !Alex, inline steering
        TIME2=float(nsteer-1)*dtsteer/3600.0 !Time in hours
!        TIME2=ctime/3600.0 !Time in hours
!      else
!        TIME2=(ctime+dtsteer)/3600.0 !Time in hours 
      endif    
!
      CALL XF_INITIALIZE(ERROR)
      IF (ITMS==1) THEN
        INQUIRE(FILE=XMDFFile,EXIST=FOUND)
        IF (FOUND) THEN
          OPEN(1000,FILE=XMDFFile)
          CLOSE(1000,STATUS='DELETE',ERR=956)
        ENDIF
      ENDIF

      CALL XF_OPEN_FILE (TRIM(XMDFFile),READWRITE,PID,ERROR)
      IF (ERROR.LT.0) THEN
        CALL XF_CREATE_FILE (TRIM(XMDFFile),READWRITE,PID,ERROR)
        WRITE(*,*) ' '
        IF (ERROR.LT.0) THEN
          WRITE(*,*) 'ERROR CREATING XMDF FILE ',TRIM(XMDFFile)
          STOP
         ELSE
          WRITE(*,*) 'BINARY FILE CREATED: ',TRIM(XMDFFile)
        ENDIF
      ENDIF

      CALL XF_OPEN_GROUP (PID,'Dataset',DGID,ERROR)
      IF (ERROR.LT.0) THEN 
        CALL XF_CREATE_GENERIC_GROUP (PID, 'Dataset', DGID, ERROR)
        IF (ERROR < 0) THEN
          WRITE(*,*) 'COULD NOT CREATE DATASET - '//'Dataset'
          STOP
        ENDIF
      ENDIF

     ! WRITE TO WAVE DATASET
      PREFIX='Wave'
      CALL XF_OPEN_GROUP (DGID, TRIM(PREFIX), DID, ERROR)
      IF (ERROR.LT.0) THEN
        CALL XF_CREATE_VECTOR_DATASET(DGID,TRIM(PREFIX),'none',  &
               TS_HOURS,0, DID,ERROR)
        if (error.gt.0) THEN
          WRITE(*,*) 'CREATED DATASET - '//TRIM(PREFIX)
          ELSE
            WRITE(*,*) 'COULD NOT CREATE DATASET - '//TRIM(PREFIX)
            STOP
        ENDIF
! removed XF_VECTORS_IN_LOCAL_COORDS (DID,ERROR) for proper vector direction
        CALL XF_VECTOR_2D_DATA_LOCS(DID,GRID_LOC_CENTER,GRID_LOC_CENTER,ERROR)
      ENDIF
      IF(noptset.eq.3) CALL XF_DATASET_REFTIME(DID,REFTIME,ERROR)
      CALL XF_WRITE_VECTOR_TIMESTEP(DID,TIME2,NIJ,2,WAVE,ERROR)
      CALL XF_CLOSE_GROUP (DID, ERROR)

      ! WRITE TO RAD STRESS DATASET
      IF (IRS .GT. 0) THEN
        PREFIX='RadStress'
        CALL XF_OPEN_GROUP (DGID, TRIM(PREFIX), DID, ERROR)
          IF (ERROR.LT.0) THEN
            CALL XF_CREATE_VECTOR_DATASET(DGID,TRIM(PREFIX),'none', &
                   TS_HOURS,0, DID,ERROR)
            IF (error.gt.0) THEN
              WRITE(*,*) 'CREATED DATASET - '//TRIM(PREFIX)
        	ELSE
      	  WRITE(*,*) 'COULD NOT CREATE DATASET - '//TRIM(PREFIX)
      	  STOP
            ENDIF
! removed XF_VECTORS_IN_LOCAL_COORDS (DID,ERROR) for proper vector direction
          CALL XF_VECTOR_2D_DATA_LOCS(DID,GRID_LOC_CENTER,GRID_LOC_CENTER,ERROR)
        ENDIF
        IF(noptset.eq.3) CALL XF_DATASET_REFTIME(DID,REFTIME,ERROR)
        CALL XF_WRITE_VECTOR_TIMESTEP(DID,TIME2,NIJ,2,radstr,ERROR)
        CALL XF_CLOSE_GROUP (DID, ERROR)
      ENDIF

      ! WRITE TO HEIGHT DATASET
      PREFIX='Height'
      CALL XF_OPEN_GROUP (DGID, TRIM(PREFIX), DID, ERROR)
      IF (ERROR.LT.0) THEN
        CALL XF_CREATE_SCALAR_DATASET(DGID,TRIM(PREFIX),'none',  &
               TS_HOURS,0, DID,ERROR)
        if (error.gt.0) THEN
            WRITE(*,*) 'CREATED DATASET - '//TRIM(PREFIX)
      	ELSE
      	  WRITE(*,*) 'COULD NOT CREATE DATASET - '//TRIM(PREFIX)
      	  STOP
          ENDIF
          CALL XF_SCALAR_DATA_LOCATION (DID,GRID_LOC_CENTER,ERROR)
        ENDIF
        IF(noptset.eq.3) CALL XF_DATASET_REFTIME(DID,REFTIME,ERROR)
        CALL XF_WRITE_SCALAR_TIMESTEP(DID,TIME2,NIJ,HEIGHT, ERROR)
        CALL XF_CLOSE_GROUP (DID, ERROR)

      ! WRITE TO PERIOD DATASET
      PREFIX='Period'
      CALL XF_OPEN_GROUP (DGID, TRIM(PREFIX), DID, ERROR)
      IF (ERROR.LT.0) THEN
        CALL XF_CREATE_SCALAR_DATASET(DGID,TRIM(PREFIX),'none',  &
               TS_HOURS,0, DID,ERROR)
        if (error.gt.0) THEN
          WRITE(*,*) 'CREATED DATASET - '//TRIM(PREFIX)
        ELSE
            WRITE(*,*) 'COULD NOT CREATE DATASET - '//TRIM(PREFIX)
            STOP
        ENDIF
        CALL XF_SCALAR_DATA_LOCATION (DID,GRID_LOC_CENTER,ERROR)
      ENDIF
      IF(noptset.eq.3) CALL XF_DATASET_REFTIME(DID,REFTIME,ERROR)
      CALL XF_WRITE_SCALAR_TIMESTEP(DID,TIME2,NIJ,PERIOD, ERROR)
      CALL XF_CLOSE_GROUP (DID, ERROR)

      ! WRITE TO DIRECTION DATASET
      PREFIX='Direction'
      CALL XF_OPEN_GROUP (DGID, TRIM(PREFIX), DID, ERROR)
      IF (ERROR.LT.0) THEN
        CALL XF_CREATE_SCALAR_DATASET(DGID,TRIM(PREFIX),'none',  &
               TS_HOURS,0, DID,ERROR)
        if (error.gt.0) THEN
          WRITE(*,*) 'CREATED DATASET - '//TRIM(PREFIX)
      	ELSE
      	  WRITE(*,*) 'COULD NOT CREATE DATASET - '//TRIM(PREFIX)
      	  STOP
          ENDIF
          CALL XF_SCALAR_DATA_LOCATION (DID,GRID_LOC_CENTER,ERROR)
        ENDIF
        IF(noptset.eq.3) CALL XF_DATASET_REFTIME(DID,REFTIME,ERROR)
        CALL XF_WRITE_SCALAR_TIMESTEP(DID,TIME2,NIJ,DIR, ERROR)
        CALL XF_CLOSE_GROUP (DID, ERROR)
      
        ! WRITE TO /DISSIPATION DATASET
        IF (IBREAK .GT. 0) THEN
          PREFIX='Breaking-Dissipation'
          CALL XF_OPEN_GROUP (DGID, TRIM(PREFIX), DID, ERROR)
          IF (ERROR.LT.0) THEN
            CALL XF_CREATE_SCALAR_DATASET(DGID,TRIM(PREFIX),'none',  &
                   TS_HOURS,0, DID,ERROR)
            if (error.gt.0) THEN
              WRITE(*,*) 'CREATED DATASET - '//TRIM(PREFIX)
        	ELSE
        	  WRITE(*,*) 'COULD NOT CREATE DATASET - '//TRIM(PREFIX)
      	      STOP
            ENDIF
            CALL XF_SCALAR_DATA_LOCATION (DID,GRID_LOC_CENTER,ERROR)
          ENDIF
          IF(noptset.eq.3) CALL XF_DATASET_REFTIME(DID,REFTIME,ERROR)
          CALL XF_WRITE_SCALAR_TIMESTEP(DID,TIME2,NIJ,BRKDISS, ERROR)
          CALL XF_CLOSE_GROUP (DID, ERROR)
        ENDIF
      
      ! WRITE TO /DISSIPATION DATASET
      IF (IBREAK .GT. 0) THEN
        PREFIX='Breaking_Dissipation'
        CALL XF_OPEN_GROUP (DGID, TRIM(PREFIX), DID, ERROR)
        IF (ERROR.LT.0) THEN
          CALL XF_CREATE_SCALAR_DATASET(DGID,TRIM(PREFIX),'none',  &
                 TS_HOURS,0, DID,ERROR)
          if (error.gt.0) THEN
            WRITE(*,*) 'CREATED DATASET - '//TRIM(PREFIX)
            ELSE
            WRITE(*,*) 'COULD NOT CREATE DATASET - '//TRIM(PREFIX)
              STOP
          ENDIF
          CALL XF_SCALAR_DATA_LOCATION (DID,GRID_LOC_CENTER,ERROR)
        ENDIF
        IF(noptset.eq.3) CALL XF_DATASET_REFTIME(DID,REFTIME,ERROR)
        CALL XF_WRITE_SCALAR_TIMESTEP(DID,TIME2,NIJ,BRKDISS, ERROR)
        CALL XF_CLOSE_GROUP (DID, ERROR)
      ENDIF
        
      ! WRITE NEW WATER DEPTH IF IN STEERING
      IF(noptset.eq.3)THEN
        if(noptzb>0)then
          PREFIX='Depth' !Alex
          CALL XF_OPEN_GROUP (DGID, TRIM(PREFIX), DID, ERROR)
          IF (ERROR.LT.0) THEN
            CALL XF_CREATE_SCALAR_DATASET(DGID,TRIM(PREFIX),'none',  &
                   TS_HOURS,0, DID,ERROR)
            if (error.gt.0) THEN
              WRITE(*,*) 'CREATED DATASET - '//TRIM(PREFIX)
              ELSE
              WRITE(*,*) 'COULD NOT CREATE DATASET - '//TRIM(PREFIX)
                STOP
            ENDIF
            CALL XF_SCALAR_DATA_LOCATION (DID,GRID_LOC_CENTER,ERROR)
          ENDIF
          IF(noptset.eq.3) CALL XF_DATASET_REFTIME(DID,REFTIME,ERROR)
          CALL XF_WRITE_SCALAR_TIMESTEP(DID,TIME2,NIJ,DEPTH, ERROR)
          CALL XF_CLOSE_GROUP (DID, ERROR)
        endif
          
        PREFIX='Water_Level' !Alex
        CALL XF_OPEN_GROUP (DGID, TRIM(PREFIX), DID, ERROR)
        IF (ERROR.LT.0) THEN
          CALL XF_CREATE_SCALAR_DATASET(DGID,TRIM(PREFIX),'none',  &
                 TS_HOURS,0, DID,ERROR)
          if (error.gt.0) THEN
            WRITE(*,*) 'CREATED DATASET - '//TRIM(PREFIX)
            ELSE
            WRITE(*,*) 'COULD NOT CREATE DATASET - '//TRIM(PREFIX)
              STOP
          ENDIF
          CALL XF_SCALAR_DATA_LOCATION (DID,GRID_LOC_CENTER,ERROR)
        ENDIF
        IF(noptset.eq.3) CALL XF_DATASET_REFTIME(DID,REFTIME,ERROR)
        CALL XF_WRITE_SCALAR_TIMESTEP(DID,TIME2,NIJ,SURGE, ERROR)
        CALL XF_CLOSE_GROUP (DID, ERROR)
          
        PREFIX='Currents'
        CALL XF_OPEN_GROUP (DGID, TRIM(PREFIX), DID, ERROR)
        IF (ERROR.LT.0) THEN
          CALL XF_CREATE_VECTOR_DATASET(DGID,TRIM(PREFIX),'none', &
                 TS_HOURS,0, DID,ERROR)
          if (error.gt.0) THEN
            WRITE(*,*) 'CREATED DATASET - '//TRIM(PREFIX)
          ELSE
              WRITE(*,*) 'COULD NOT CREATE DATASET - '//TRIM(PREFIX)
              STOP
          ENDIF
! removed XF_VECTORS_IN_LOCAL_COORDS (DID,ERROR) for proper vector direction
          CALL XF_VECTOR_2D_DATA_LOCS(DID,GRID_LOC_CENTER,GRID_LOC_FACE_J,ERROR)
        ENDIF
        IF(noptset.eq.3) CALL XF_DATASET_REFTIME(DID,REFTIME,ERROR)
        CALL XF_WRITE_VECTOR_TIMESTEP(DID,TIME2,NIJ,2,Currents,ERROR)
        CALL XF_CLOSE_GROUP (DID, ERROR)    
      ENDIF

      CALL XF_CLOSE_GROUP (DGID, ERROR)
      CALL XF_CLOSE_FILE (PID,ERROR)
      return

956   write(*,'(A,A50)') 'ERROR: Could not access: ',XMDFFile
      write(*,'(A)') 'Press <RETURN> to continue'
      read(*,*)
      return 
      end subroutine 
#endif

      END  !SUBROUTINE OR PROGRAM, Alex
!
!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!  INITIAL WAVE CONDITION AT UPWAVE BOUNDARY
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE INITL_inline
!--------------------------------------------------
!  setup of mesh size, angular-frequency spectrum
!--------------------------------------------------
      USE GLOBAL_INLINE, ONLY: NPF,MPD,IPMX,JPMX,IGPX,JGPX
      COMMON /DATA/NF,MD,IMAX,JMAX,IGMX,JGMX,JCPB,JCPE,JCB,JCE,NFF,MDD
      COMMON /DATB/ICK3,ICK4,ICHOICE,HS0,wd,ws,ph0,wdd(mpd),aslop(ipmx)
      COMMON /DATC/TP,PRD,IBND,VL,TIDE,KPEAK,IBACK,NBLOCK,IWIND,WSMAG
      COMMON /DATD/DX,DY,DXX,DMESH,DTH,kdate,idate,depmax,depmax0
      COMMON /DATG/IBK,DBAR,WL0,WCC(JGPX),HSB(JGPX),HSG(JGPX),DCD(JGPX)
      COMMON /DVAR/DVARX(IGPX),DVARY(JGPX),ETA(IGPX,JGPX),HSK(JGPX)
      COMMON /SPECA/SCP(JGPX,MPD),SI(JGPX,NPF,MPD)
      COMMON /SPECD/SIMAX(JGPX,NPF,MPD),SJJ(JGPX)
      COMMON /FDS/FCN(NPF),DCM(MPD),DSFD(NPF,MPD)
      COMMON /BREK/DEPM(JGPX),DMNJ(JGPX),SLF(JGPX),wlmn(JGPX),cmn(jgpx)  &
                   ,sigm(jgpx),IWVBK,ibk3
      COMMON /VPAI/PAI2,PAI,HPAI,RAD,akap,imod,iprp,island,imd,iprpp     &
                   ,nonln,igrav,isolv,ixmdf,iproc,imud,iwnd,depmin0
      common /comsav/depmin,g,iview,iplane,iwave,cflat,irunup,iroll

! ** mesh size
      DX=DMESH
      DY=DX
      dxx=dvarx(1)

! -- initial input spectrum, height, direction at all j points
      DO 20 MM=imd,MD
      DO 20 NN=1,NF
      DO 20 JJ=1,JGMX
        !- si(jj,nn,mm) is initial input spectrum for calculation in first region
        SI(JJ,NN,MM)=DSFD(NN,MM)*wcc(jj)
        simax(jj,nn,mm)=si(jj,nn,mm)
   20 CONTINUE

      IF(IBND.EQ.2) THEN
        do 22 jj=1,jgmx    
          if(abs(dcd(jj)).lt..01) go to 22
          mmdcd=int(dcd(jj))
          do 28 nn=1,nf
            if(iabs(mmdcd).lt.1) go to 24
            slf=0.
            do mm=imd,md           !23
              m1=mm+mmdcd
              if(m1.lt.imd) m1=imd
              if(m1.gt.md) m1=md
              slf(m1)=slf(m1)+si(jj,nn,mm)
            enddo                  !23
      
            do mm=imd,md           !25
              si(jj,nn,mm)=slf(mm)
            enddo                  !25
   24       continue
            dcdmm=dcd(jj)-float(mmdcd)
            if(abs(dcdmm).lt..01) go to 30
            if(dcdmm.gt.0.) then
              c1=1.-dcdmm
              do mm=imd+1,md       !26
                slf(mm)=si(jj,nn,mm)*c1+si(jj,nn,mm-1)*dcdmm
              enddo                !26
              slf(imd)=si(jj,nn,imd)*c1
            else
              dcdmm=-dcdmm
              c1=1.-dcdmm
              do mm=imd,md-1       !27
                slf(mm)=si(jj,nn,mm)*c1+si(jj,nn,mm+1)*dcdmm
              enddo                !27
              slf(md)=si(jj,nn,md)*c1
            end if
            do mm=imd,md           !29
              si(jj,nn,mm)=slf(mm)
            enddo                  !29
   30       continue
            do mm=imd,md           !31
              simax(jj,nn,mm)=si(jj,nn,mm)
            enddo                  !31
   28     continue
   22   continue
      END IF

      DO JJ=1,JGMX                 !21
        HSB(JJ)=HS0
        DMNJ(JJ)=PRD
        SJJ(JJ)=(HS0/cflat)**2
      ENDDO                        !21

      RETURN
      END


!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!  wave spectrum calculation in forward (incident waves) 
!    and backward (reflected waves) direction
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE SETIN_inline(IKAP)
!-------------------------------------------------------
!  in this subroutine, calc(ii) is the main subroutine 
!    which calculates wave action balance equation
!-------------------------------------------------------
      USE GLOBAL_INLINE, ONLY: NPF,MPD,MPD2,IPMX,JPMX,KOMX,NOMX,IGPX,JGPX
      use wave_lib, only: wave_Hmax
      COMMON /VPAI/PAI2,PAI,HPAI,RAD,akap,imod,iprp,island,imd,iprpp     &
                   ,nonln,igrav,isolv,ixmdf,iproc,imud,iwnd,depmin0
      COMMON /DATA/NF,MD,IMAX,JMAX,IGMX,JGMX,JCPB,JCPE,JCB,JCE,NFF,MDD
      COMMON /DATB/ICK3,ICK4,ICHOICE,HS0,wd,ws,ph0,wdd(mpd),aslop(ipmx)
      COMMON /DATC/TP,PRD,IBND,VL,TIDE,KPEAK,IBACK,NBLOCK,IWIND,WSMAG
      COMMON /DATD/DX,DY,DXX,DMESH,DTH,kdate,idate,depmax,depmax0
      COMMON /DATG/IBK,DBAR,WL0,WCC(JGPX),HSB(JGPX),HSG(JGPX),DCD(JGPX)
      COMMON /DVAR/DVARX(IGPX),DVARY(JGPX),ETA(IGPX,JGPX),HSK(JGPX)
      COMMON /SPECA/SCP(JGPX,MPD),SI(JGPX,NPF,MPD)
      COMMON /SPECB/SOP(KOMX,NPF,MPD2),SR(2*IPMX,NPF,MPD)
      COMMON /SPECC/SJ(JGPX),SJF(JGPX,NPF),FJF(JGPX)
      COMMON /SPECD/SIMAX(JGPX,NPF,MPD),SJJ(JGPX)
      COMMON /OUTP/KOUT,IJSP(2,KOMX),IRS,IBREAK,ICUR,IWET,INST
      COMMON /OUTN/nest,inest(komx),jnest(komx),ix1(igpx),ix2(igpx)
      COMMON /REFA/KRMX,KR(2,6*IPMX),RK(6*IPMX),yangl(6*IPMX)
      COMMON /REFB/KRF,KCR(2,6*IPMX),RKR(6*IPMX),xangl(6*IPMX)
      COMMON /IJBP/IJB(IGPX,JGPX),DEP(IPMX,JPMX),DEPS(IPMX),DBIG(IPMX)
      common /uvp/u(ipmx,jpmx),v(ipmx,jpmx)
      COMMON /WAVI/H13(IGPX,JGPX),T13(IGPX,JGPX),DMN(IGPX,JGPX)
      COMMON /WAVR/H13R(IGPX,JGPX),T13R(IGPX,JGPX),DMNR(IGPX,JGPX)
      COMMON /WAVS/H13S(IGPX,JGPX),IBR(IGPX,JGPX),DISS(IGPX,JGPX)
      COMMON /WAVF/H13F(IGPX,JGPX),DMNF(IGPX,JGPX),T13F(IGPX,JGPX),T13MIN
      COMMON /BREK/DEPM(JGPX),DMNJ(JGPX),SLF(JGPX),wlmn(JGPX),cmn(jgpx)  &
                   ,sigm(jgpx),IWVBK,ibk3
      COMMON /FDS/FCN(NPF),DCM(MPD),DSFD(NPF,MPD)
      COMMON /DFCN/DF(NPF),FFCN(NPF),DFINP(NPF),PL1E(NPF),PL1T(NPF)
      common /rsa/ sxx(ipmx,jpmx),sxy(ipmx,jpmx),syy(ipmx,jpmx)
      common /rsm/ cosa(mpd),sina(mpd),d1(ipmx,jpmx)
      common /rsw/ wk2(jpmx,npf,mpd),cgp(ipmx,jpmx)
      common /uv1/u1(ipmx,jpmx),v1(ipmx,jpmx)
      common /struc0/ijstruc1,ijstruc2,ijstruc3,ijstruc4,ismall
      common /struc2/istruc2(NOMX),jstruc2(NOMX),dstruc2(NOMX)
      common /struc3/istruc3(komx),jstruc3(komx),dstruc3(komx) &
                    ,k3(komx),dstruc33(komx)
      common /struc4/istruc4(komx),jstruc4(komx),dstruc4(komx) &
                    ,kstruc4(komx),k4(komx),dstruc44(komx)
      common /swell/sw13(igpx,jgpx),sa13(igpx,jgpx) &
                   ,tw13(igpx,jgpx),ta13(igpx,jgpx) &
                   ,dw13(igpx,jgpx),da13(igpx,jgpx)
      common /comsav/depmin,g,iview,iplane,iwave,cflat,irunup,iroll
!
      INTEGER, ALLOCATABLE :: KRC(:),JR(:)
      REAL, ALLOCATABLE :: RKK(:),angl(:),hsgg(:),ijbp(:,:),   &
                           h13a(:),h13b(:),dd11(:)
      INTEGER II,JJ,KK,iddd
!
!-----------------------------------------------------
!  dep(i,j) is temporally defined water depth for
!    forward and backward calculations, respectively
!-----------------------------------------------------
!  replace of variables
!-----------------------------------------------------
      ALLOCATE (KRC(2*IPMX),JR(2*IPMX),RKK(2*IPMX),angl(2*IPMX))
      ALLOCATE (hsgg(jgpx),h13a(jgpx),h13b(jgpx),dd11(jgpx))
      ALLOCATE (IJBP(IGPX,JGPX))
      IJBP=IJB
      t13min=1./fcn(nf)
      H13=0.
      T13=t13min
      DMN=0.
      ETA=0.
      SOP=0.
      HSK=HS0
      H13A=HS0
      H13B=HS0

! -- reflection in y direction
  102 IF(ICK4.NE.1) GOTO 104

! -- reflection in x direction
      KRC=0
      JR=0
      RKK=0
  104 CONTINUE

!=================================================
!  start of forward marching calculation
!    ibk=0 means that wave breaking is accounted
!=================================================
      IBK=0
      II=1
      leap=1
      island=0

! -- calc(ii) is main subroutine which returns calculated values
!------------------------------------------------------------------------
!  by CALC, si(jj,nn,mm),sj(jj) are obtained
!------------------------------------------------------------------------
   60 continue

      CALL CALC_inline(II,IKAP)
!------------------------------------------------------------------------
!  by CALC, si(jj,nn,mm),sj(jj) are obtained
!------------------------------------------------------------------------
!
      if(iplane.le.1) go to 411
!
      if(ii.ge.2.and.ws.ge.3.) then
        cc=1.02+abs(sin(wd))*ws*0.001
        do 401 jj=2,jgmx
          if(sin(dbig(jj)).le.0.) goto 401
          if(sj(jj).lt..001) go to 401
          if(d1(ii,jj-1).lt.1.) go to 401
          if(sj(jj-1).lt..01) go to 401
!         cc1=d1(ii,jj)-d1(ii,jj-1)
!         if(cc1.gt.10.) go to 401
          if(sj(jj)/sj(jj-1).gt.cc) then
            sj(jj)=sj(jj-1)*cc
            do l=imd,md
             do k=1,nf
               si(jj,k,l)=si(jj-1,k,l)*cc
             end do
            end do
          end if
  401   continue
!
        do 408 jj=jgmx-1,1,-1
          if(sin(dbig(jj)).ge.0.) goto 408
          if(sj(jj).lt..001) go to 408
          if(d1(ii,jj+1).lt.1.) go to 408
          if(sj(jj+1).lt..01) go to 408
!         cc1=d1(ii,jj)-d1(ii,jj+1)
!         if(cc1.gt.10.) go to 408
          if(sj(jj)/sj(jj+1).gt.cc) then
            sj(jj)=sj(jj+1)*cc
            do l=imd,md
             do k=1,nf
               si(jj,k,l)=si(jj+1,k,l)*cc
             end do
            end do
          end if
  408   continue
!
        do 412 jj=1,jgmx
          if(d1(ii-1,jj).lt..01) then
            if(d1(ii,jj).lt..1) goto 412
            cc=sqrt(sj(jj))/d1(ii,jj)
            if(cc.gt..1) then
              sj(jj)=.01*d1(ii,jj)**2
              do l=imd,md
               do k=1,nf
                 si(jj,k,l)=si(jj+1,k,l)/cc**2*.01
               end do
              end do
            end if
          end if
  412   continue
      end if
  411 continue
!
      if(leap.eq.1) then
        if(ws.ge..1) then
          if(ii.eq.1) then
            if(ibnd.eq.0) then
              do 405 jj=1,jgmx
                if(d1(ii,jj).lt..01) go to 405
                SS=SJ(JJ)
                if(SS.lt.1.0E-15) SS=1.0E-15
                hsi=cflat*sqrt(SS)
                cc=hsi/hsg(jj)
                if(cc.lt.1.2) go to 405
                sj(jj)=(hsg(jj)/cflat)**2
                do l=imd,md
                  do k=1,nf
                    si(jj,k,l)=si(jj,k,l)/cc**2
                  end do
                end do
  405         continue
            end if
          else
            do 505 jj=1,jgmx
              if(d1(ii,jj).lt..01) go to 505
              SS=SJ(JJ)
              if(SS.lt.1.0E-5) go to 505
              hsi=cflat*sqrt(SS)
              if(h13(ii-1,jj)/hsi.gt..85) go to 505
              if(d1(ii-1,jj).lt..01) then
                h0=0.5
              else
                h0=h13(ii-1,jj)
                if(h0.lt..012.and.d1(ii,jj).ge..1) h0=.012
                if(h0.lt..012) go to 505
              end if
              cc=hsi/h0
              if(cc.lt.1.5) go to 505
!             sj(jj)=h0**2/16.*2.
              sj(jj)=(h0/cflat)**2*2.
              do l=imd,md
                do k=1,nf
                  si(jj,k,l)=si(jj,k,l)/cc**2*2.
                end do
              end do
  505       continue
          end if
        else
          if(ii.ge.2) then
            do 605 jj=1,jgmx
              if(d1(ii,jj).lt..01) go to 605
              SS=SJ(JJ)
              hsi=cflat*sqrt(SS)
              if(hsi.lt..5) go to 605
              h0=h13(ii-1,jj)
              if(h0.lt..5) go to 605
              cc=hsi/h0
              if(cc.lt.1.15) go to 605
!             sj(jj)=h0**2/16.*1.32
              sj(jj)=(h0/cflat)**2*1.32
              do l=imd,md
                do k=1,nf
                  si(jj,k,l)=si(jj,k,l)/cc**2*1.32
                end do
              end do
  605       continue
          end if
        end if
      end if

      if(nblock.eq.1) then
        do 509 jj=5,jgmx-4
          if(ijb(ii,jj).eq.2) then
            if(ijb(ii,jj+1).eq.0) then
              do l=imd,md
                do k=1,nf
                  si(jj,k,l)=si(jj-1,k,l)
                end do
              end do
              sj(jj)=sj(jj-1)
            end if
            if(ijb(ii,jj-1).eq.0) then
              do l=imd,md
                do k=1,nf
                  si(jj,k,l)=si(jj+1,k,l)
                end do
              end do
              sj(jj)=sj(jj+1)
            end if
          end if
509     continue
      end if

      sum1=0.
      sum2=0.
      do 47 mm=imd,md
        do 47 nn=1,nf
          sum1=sum1+wk2(1,nn,mm)
          sum2=sum2+wk2(2,nn,mm)
47    continue
      sum1=sum1/float(nf*(md-imd+1))
      sum2=sum2/float(nf*(md-imd+1))

!     itrack=0
      if(abs(sum1/sum2-1.).gt..05) then
!       itrack=1
        do 49 mm=imd,md
          do 49 nn=1,nf
            wk2(1,nn,mm)=wk2(2,nn,mm)
49      continue
        do l=imd,md
          do k=1,nf
            si(1,k,l)=si(2,k,l)
          end do
        end do
        sj(1)=sj(2)
      end if

!     if(igrav.eq.3) then
!       if(iview.eq.1.and.iplane.eq.2) then
!     DO 404 JJ=1,JGMX
! --- Nonlinear wave-wave interactions for same direction
!     if(d1(ii,jj).lt..01) go to 404
!     do 407 mm=imd,md
!     do 406 nn=nf,3,-1
!     if(si(jj,nn-1,mm).gt.si(jj,nn,mm)*10.) goto 407
!     cz=.2/tanh(.1*wk2(jj,nn,mm)*d1(ii,jj))
!     cc=si(jj,nn,mm)*tanh(cz*fcn(nn-1)**2)
!     si(jj,nn-1,mm)=si(jj,nn-1,mm)+cc*0.5
!     si(jj,nn-2,mm)=si(jj,nn-2,mm)+cc*0.5
!     si(jj,nn,mm)=si(jj,nn,mm)-cc
!406   continue
!407   continue
!404   continue
!       end if
!       end if

      if(leap.gt.1) then
        cc=float(leap)
        DO 28 MM=imd,MD
        DO 28 NN=1,NF
        DO 28 JJ=1,JGMX
          c1=(si(jj,nn,mm)-simax(jj,nn,mm))/cc
          si(jj,nn,mm)=simax(jj,nn,mm)
          simax(jj,nn,mm)=c1
   28   CONTINUE
        DO 29 JJ=1,JGMX
          c1=(sj(jj)-sjj(jj))/cc
          sj(jj)=sjj(jj)
          sjj(jj)=c1
   29   CONTINUE
      end if

      do 339 iii=ii-leap+1,ii
        if(leap.gt.1) then
          DO 22 MM=imd,MD
          DO 22 NN=1,NF
          DO 22 JJ=1,JGMX
            si(jj,nn,mm)=si(jj,nn,mm)+simax(jj,nn,mm)
   22     CONTINUE
          DO 23 JJ=1,JGMX
            sj(jj)=sj(jj)+sjj(jj)
   23     CONTINUE
        end if

        do 712 k=1,ijstruc3
          if(imod.eq.0) then
            IF(istruc3(k).NE.III) GOTO 712
            JOP=jstruc3(k)
            c3=dvarx(iii)
          end if
          if(imod.eq.2) then
            IF(istruc3(k).NE.IGMX-III+1) GOTO 712
            JOP=JGMX-jstruc3(k)+1
            c3=dvarx(iii)
          end if
          if(imod.eq.1) then
            IF(jstruc3(k).NE.III) GOTO 712
            JOP=JGMX-istruc3(k)+1
            c3=dvary(jop)
          end if
          if(imod.eq.3) then
            IF(jstruc3(k).NE.IGMX-III+1) GOTO 712
            JOP=istruc3(k)
            c3=dvary(jop)
          end if
          if(k3(k).ge.6) then
            cc=dstruc33(k)
            if(cc.lt.0.) cc=0.
            c1=0.64*(h13(iii-1,jop)/c3)**0.31+0.4*cc/(h13(iii-1,jop)+.01)
!           tc=2.44*sqrt(c3/h13(iii-1,jop))
            c2=cflat*sqrt(sj(jop))/(h13(iii-1,jop)+.01)
!           write(*,*) iii,jop,c1,c2
            if(c1.gt..99) c1=.99
            if(c1.lt..01) c1=.01
            tc=1/c1**2
            if(k3(k).eq.6.and.c2.gt..9) tc=1.25*c2**2
            if(k3(k).eq.7.and.c2.gt.1.01) tc=1.5*c2**2
          else
            ibk3=1
            ipeak=int((dmnj(jop)+hpai-dth)/dth+1.5)
            if (ipeak.lt.imd) ipeak=imd
            if (ipeak.gt.md) ipeak=md
            fp=1./TP
            if(iii.gt.1) fp=1./t13(iii-1,jop)
            small=1.
            kkpeak=kpeak
            do nn=1,nf
              cc=abs(fp-fcn(nn))
              if(cc.lt.small) then
                kkpeak=nn
                small=cc
              end if
            end do
            cc=abs(wk2(jop,kkpeak,ipeak))
            c1=cc*d1(iii,jop)
            c2=cc*abs(d1(iii,jop)-dstruc3(k))
            tc=1.+(cc*c3*sinh(min(c1,10.))/cosh(min(c2,10.)))**2
!           tc=(2.*c1+sinh(min(2*c1,10.)))/(2.*c2+sinh(min(2*c2,10.)))
!           write(*,*) 'iii=',iii,jop,cc,tc
          end if
          sj(jop)=sj(jop)/tc
          DO 713 MM=imd,MD
          DO 713 NN=1,NF
            SI(JOP,NN,MM)=SI(JOP,NN,MM)/tc
  713     CONTINUE
  712   continue

! --* output of spectrum at prescribed positions
        DO 20 K=1,KOUT
          if(imod.eq.0) then
            IF(IJSP(1,K).NE.III) GOTO 20
            JOP=IJSP(2,K)
          end if
          if(imod.eq.2) then
            IF(IJSP(1,K).NE.IGMX-III+1) GOTO 20
            JOP=JGMX-IJSP(2,K)+1
          end if
          if(imod.eq.1) then
            IF(IJSP(2,K).NE.III) GOTO 20
            JOP=JGMX-IJSP(1,K)+1
          end if
          if(imod.eq.3) then
            IF(IJSP(2,K).NE.IGMX-III+1) GOTO 20
            JOP=IJSP(1,K)
          end if
          DO 30 MM=imd,MD
          DO 30 NN=1,NF
            SOP(K,NN,MM)=SI(JOP,NN,MM)/DF(NN)/DTH
   30     CONTINUE
   20   CONTINUE
! --* no output of spectrum, just stores wave quantities

! --* output of spectrum at prescribed positions
        DO 920 K=1,nest
          if(imod.eq.0) then
            IF(inest(K).NE.III) GOTO 920
            JOP=jnest(K)
          end if
          if(imod.eq.2) then
            IF(inest(K).NE.IGMX-III+1) GOTO 920
            JOP=JGMX-jnest(K)+1
          end if
          if(imod.eq.1) then
            IF(jnest(K).NE.III) GOTO 920
            JOP=JGMX-inest(K)+1
          end if
          if(imod.eq.3) then
            IF(jnest(K).NE.IGMX-III+1) GOTO 920
            JOP=inest(K)
          end if
          DO 930 MM=imd,MD
          DO 930 NN=1,NF
            SOP(K+iabs(KOUT),NN,MM)=SI(JOP,NN,MM)/DF(NN)/DTH
  930     CONTINUE
  920   CONTINUE
! --* no output of spectrum, just stores wave quantities

        write(*,9999) iii
 9999   format ('Column ', i0)
        if(mod(iii,10).eq.0) write(9,9999) iii
!
        j1=jgmx/2
        ss=sj(1)
        if(ss.lt.1.0E-15) ss=1.0E-15
        hj1=0.
        h13ave=0.
        h13add=0.
        DO 40 JJ=1,JGMX
! -- total energy at j position
          if(d1(iii,jj).lt..01) go to 40
          SS=SJ(JJ)
!         IF(SS.LT.1.0E-15) GOTO 40
          if(ss.lt.1.0E-15) ss=1.0E-15
          HSI=cflat*SQRT(SS)
          cc1=.3+.2*tanh(hsi)
          if(iwvbk.ge.5) cc1=2.
          cc2=.3+cc1
! --- cc1=.31 & cc2=.64 in 1Jan2011 Version
          edplim=cc2*d1(iii,jj)
          edplim1=edplim
          if(abs(dmnj(jj)).gt.hpai) dmnj(jj)=wd
          if(iwvbk.eq.6) go to 490
          i1=iii+1
          if(iwvbk.le.4) then
            dep11=dep(iii,jj)
            dep12=dep(i1,jj)
            dep21=dep(iii,jj+1)
            dep22=dep(i1,jj+1)
            slx=(dep11+dep21-dep12-dep22)/dvarx(iii)/2.
            sly=(dep11+dep12-dep21-dep22)/dvary(jj)/2.
            slj=slx*cos(dmnj(jj))+sly*sin(dmnj(jj))
            if(slj.ge..04) edplim=d1(iii,jj)*(.6+.5*tanh((slj/.05)**2))
          end if
!
          if(iii.ge.2) then 
            if (ibr(iii-1,jj).eq.1) edplim=cc2*d1(iii,jj)
          end if
          edplim1=edplim
          ipeak=int((dmnj(jj)+hpai-dth)/dth+1.5)
          if (ipeak.lt.imd) ipeak=imd
          if (ipeak.gt.md) ipeak=md
          fp=1./TP
          if(iii.gt.1) fp=1./t13(iii-1,jj)
          small=1.
          kkpeak=kpeak
          do nn=1,nf
            cc=abs(fp-fcn(nn))
            if(cc.lt.small) then
              kkpeak=nn
              small=cc
            end if
          end do
          wkpeak=abs(wk2(jj,kkpeak,ipeak))
          cgp(iii,jj)=wkpeak
          if(wkpeak.ne.0.) then
            cc=cc1
            if(d1(iii,jj).lt.hsi) cc=cc1*hsi/d1(iii,jj)
            estlim=cc*pai/wkpeak*tanh(wkpeak*d1(iii,jj))
            if(d1(iii,jj).lt.hs0) cc=cc1*(d1(iii,jj)+hs0)/hs0/2.
            estlim1=cc*pai/wkpeak*tanh(wkpeak*d1(iii,jj))
            edplim=min(edplim,estlim)
            if(irs.le.1) edplim1=min(edplim,estlim1)
          end if
!
          if(hsi.gt.edplim) then
            if(min(hsi,hs0/3.).gt.d1(iii,jj)*1.3*tanh(tp/8.)) ibr(iii,jj)=1
          end if
!
          if(hsi.gt.edplim1) ibr(iii,jj)=1
!
          cc=min(edplim,edplim1)
!
          if(hsi.gt.cc) then
            if(d1(iii,jj).ge..01) then
              do l=imd,md
                do k=1,nf
                  si(jj,k,l)=si(jj,k,l)*(cc/hsi)**2
                end do
              end do
            end if
            hsi=cc
          end if
!
490       continue
!
          if(ibr(iii,jj).eq.1) go to 410
          m2=iii+3
          if(m2.gt.igmx) m2=igmx
          n1=jj-3
          n2=jj+3
          if(n1.lt.1) n1=1
          if(n2.gt.jgmx) n2=jgmx
          do 409 nn=n1,n2
          do 409 mm=iii,m2
            if(d1(mm,nn).lt..01.and.ijb(mm,nn).ne.-1) then
              ic1=(mm-iii)**2+(nn-jj)**2
              ic2=((mm-iii)*dvarx(iii))**2+((nn-jj)*dvary(jj))**2
              if(ic1.le.2.or.ic2.le.10.*hsi**2) then
!     above line edited 25 May 2009 - Lihwa
                ibr(iii,jj)=1
                go to 410
              end if
            end if
409       continue
410       continue
!
! -- calculated significant wave height
          H13(III,JJ)=HSI
          if(hsi.gt.hj1) then
            j1=jj
            hj1=hsi
          end if 
          if(hsi.gt..005) then
            h13ave=h13ave+hsi
            h13add=h13add+1.
          end if
! -- hsb(jj) is used in breaker index
          HSB(JJ)=HSI
   40   CONTINUE
!
        do 501 jj=1,jgmx
          fjf(jj)=fcn(nf)
          big=0.
          nbig=nf
          sum2=0.
          do 502 nn=nf,1,-1
            sjf(jj,nn)=0.
            do 506 mm=imd,md
              sjf(jj,nn)=sjf(jj,nn)+si(jj,nn,mm)
  506       continue
            if(fcn(nn).lt..04) go to 511
            if(sjf(jj,nn).gt.big) then
              big=sjf(jj,nn)
              fjf(jj)=fcn(nn)
              nbig=nn
            end if
  511       continue
            if(irs.ge.2)then
              sum2=sum2+fcn(nn)*sjf(jj,nn)**2/df(nn)
            end if
  502     continue
!
          if(iwet.lt.0) then
            if(fjf(jj).gt..09) then
              do nn13=1,nf
                nbig=nn13
                if(fcn(nn13).gt..09) exit
              end do
              if(nbig.eq.1) nbig=0
            else
              if(nbig+1.lt.nf) then
                if(fcn(nbig+1).le..11) nbig=nbig+1
              end if
            end if
            sum=0.
            ssum=0.
            do nn13=1,nbig
              sum=sum+sjf(jj,nn13)
              ssum=ssum+sjf(jj,nn13)*fcn(nn13)
            end do
            sum1=0.
            ssum1=0.
            do nn13=nbig+1,nf
              sum1=sum1+sjf(jj,nn13)
              ssum1=ssum1+sjf(jj,nn13)*fcn(nn13)
            end do
            cc1=sum*cflat**2
            cc2=sum1*cflat**2

            cc=h13(iii,jj)**2/(cc1+cc2+.0001)
            if(cc.lt..99) then
              cc1=cc1*cc
              cc2=cc2*cc
            end if

            cc3=sw13(iii,jj)**2
            cc4=sa13(iii,jj)**2
            ctw13=10.
            cta13=1.
            if(sum.gt..0001) ctw13=sum/ssum
            if(sum1.gt..0001) cta13=sum1/ssum1
            if(iplane.le.1) then
              tw13(iii,jj)=ctw13
              ta13(iii,jj)=cta13
            else
              tw13(iii,jj)=(ctw13*cc1+tw13(iii,jj)*cc3)/(cc1+cc3+0.01)
              ta13(iii,jj)=(cta13*cc2+ta13(iii,jj)*cc4)/(cc2+cc4+0.01)
            end if
            if(iplane.le.1) then
              sw13(iii,jj)=sqrt(cc1)
              sa13(iii,jj)=sqrt(cc2)
            else
              sw13(iii,jj)=sqrt(cc1+cc3)
              sa13(iii,jj)=sqrt(cc2+cc4)
            end if
!           if(jj.eq.j1) write(*,*) 'sw13,sa13=',sw13(iii,jj),sa13(iii,jj)
!           if(jj.eq.j1) write(*,*) 'tw13,ta13=',tw13(iii,jj),ta13(iii,jj)
          end if
!
          if(irs.ge.2) then
            IF(IBR(III,JJ).EQ.0.and.h13(iii,jj).gt..001) then
              qp=2.*sum2/(h13(iii,jj)/cflat)**4
              hsk(jj)=(h13(iii,jj)/sqrt(qp)+hsk(jj))/2.
              if(abs(h13a(jj)-h13(iii,jj)).lt.h13a(jj)*.05) then
                h13a(jj)=(h13(iii,jj)+h13a(jj))/2.
              end if
              h13b(jj)=h13(iii,jj)*.6666+h13b(jj)*.3333
            END IF
          end if
          t13(iii,jj)=1./fjf(jj)
!
          if(ibk3.eq.0) then
            if(hs0.lt..3.and.ws.gt..1) then
              cc=(sqrt(h13(iii,jj)/.04)+t13min)/2.
              if(cc.lt.1.) cc=1.
              if(t13(iii,jj).gt.cc) t13(iii,jj)=cc
              fjf(jj)=1./t13(iii,jj)
            end if
          end if
!
          sumx=0.
          sumy=0.
          sbig=0.
          dbig(jj)=wd
          do 507 mm=imd,md
            scp(jj,mm)=0.
            do 503 nn=1,nf
              scp(jj,mm)=scp(jj,mm)+si(jj,nn,mm)
  503       continue
            if(scp(jj,mm).gt.sbig) then
              sbig=scp(jj,mm)
              dbig(jj)=dcm(mm)
            end if
            sumx=sumx+scp(jj,mm)*cosa(mm)
            sumy=sumy+scp(jj,mm)*sina(mm)
  507     continue
          dmnj(jj)=atan2(sumy,sumx)
          dmn(iii,jj)=dmnj(jj)/rad
!
          if(iwet.lt.0) then
            sumx=0.
            sumy=0.
            do 1507 mm=imd,md
              scp(jj,mm)=0.
              do 1503 nn=1,nbig
                scp(jj,mm)=scp(jj,mm)+si(jj,nn,mm)
 1503         continue
              sumx=sumx+scp(jj,mm)*cosa(mm)
              sumy=sumy+scp(jj,mm)*sina(mm)
 1507       continue
            cdw13=atan2(sumy,sumx)
            sumx=0.
            sumy=0.
            do 2507 mm=imd,md
              scp(jj,mm)=0.
              do 2503 nn=nbig+1,nf
                scp(jj,mm)=scp(jj,mm)+si(jj,nn,mm)
 2503         continue
              sumx=sumx+scp(jj,mm)*cosa(mm)
              sumy=sumy+scp(jj,mm)*sina(mm)
 2507       continue
            cda13=atan2(sumy,sumx)
!
            if(iplane.le.1) then
              dw13(iii,jj)=cdw13/rad
              da13(iii,jj)=cda13/rad
            else
              cx=cc1*cos(cdw13)+cc3*cos(dw13(iii,jj)*rad+pai)
              cy=cc1*sin(cdw13)+cc3*sin(dw13(iii,jj)*rad+pai)
              dw13(iii,jj)=atan2(cy,cx)/rad
              cx=cc2*cos(cda13)+cc4*cos(da13(iii,jj)*rad+pai)
              cy=cc2*sin(cda13)+cc4*sin(da13(iii,jj)*rad+pai)
              da13(iii,jj)=atan2(cy,cx)/rad
            end if
          end if
!
          if(ICK4.eq.1) then
! -- jbe is reflection points in x direction?
            DO 80 K=1,KRF
              IF(KCR(1,K).NE.III) GOTO 80
              JK=KCR(2,K)
              IF(JK.EQ.JJ) go to 70
   80       CONTINUE
            GOTO 501
            
   70       CONTINUE

            DO 609 NN=1,NF
            DO 609 MM=imd,MD
!             SR(K,NN,MM)=SI(JK,NN,MM)*cos(dmnj(jj)-xangl(k)*rad)**2
              SR(K,NN,MM)=SI(JK,NN,MM)
  609       CONTINUE
          end if
!
  501   continue
!
        if(ws.lt..1) go to 403
        DO 400 jj=1,jgmx
          hsgg(jj)=hsg(jj)
          if(d1(iii,jj).gt.2.5) go to 400
          d11=d1(iii,jj)
          if(d11.lt..01) d11=.01
          cc=5./(d11+2.5)
          if(jj.eq.1) then
            hsgg(jj)=hsg(2)/cc
            go to 400
          end if
          if(jj.eq.jgmx) then
            hsgg(jj)=hsg(jgmx-1)/cc
            go to 400
          end if
          if(d1(iii,jj-1).lt..5.and.d1(iii,jj+1).lt..5) then
            hsgg(jj)=0.
            go to 400
          end if
          if(d1(iii,jj-1).lt.2.5.and.d1(iii,jj+1).lt.2.5) then
            hsgg(jj)=(hsg(jj)+hsg(jj-1)+hsg(jj+1))/3./cc
            go to 400
          end if
          if(d1(iii,jj+1).lt.2.5) hsgg(jj)=(hsg(jj)+hsg(jj+1))/2./cc
          if(d1(iii,jj-1).lt.2.5) hsgg(jj)=(hsg(jj)+hsg(jj-1))/2./cc
  400   CONTINUE
!
        do 402 jc=2,jgmx-1
          hsg(jc)=(hsgg(jc)+hsgg(jc-1)+hsgg(jc+1))/3.
  402   continue
        hsg(1)=(hsgg(1)+hsgg(2))/2.
        hsg(jgmx)=(hsgg(jgmx)+hsgg(jgmx-1))/2.
  403 continue
!
      if(h13add.gt..5) h13ave=h13ave/h13add
      if(iplane.le.1) then
        write(*,9020) h13ave
9020    format('Average wave height = ',2f10.4)
        if(mod(iii,10).eq.0) write(9,9020) h13ave
      else
        write(*,9021) h13ave
9021    format('Seaward wave height = ',2f10.4)
        if(mod(iii,10).eq.0) write(9,9021) h13ave
      end if

      i3=iii+1

      if(irs.ge.1) then
        call sxycalc_inline(iii)
        if(irs.ge.2) then
          if(iii.gt.1.and.iii.lt.igmx) then
            deplow=.05
            if(depmax.gt.1.) deplow=.5
            i1=iii-1
            hsgg=0.
            do jj=2,jgmx-1
              if(d1(iii,jj).gt..001) then
                c1=sxx(iii,jj)-sxx(i1,jj)
                c2=(sxy(iii,jj+1)-sxy(iii,jj-1))*dvarx(iii)/dvary(jj)
                cc=d1(i1,jj)-dvarx(iii)*.05
                if(cc.lt.deplow) cc=deplow
                c3=(c1+c2/2.)/max(d1(iii,jj),cc)
                if(c3.lt.-.001) c3=-.001
                if(c3.gt..001) c3=.001
                eta(iii,jj)=eta(i1,jj)-c3
                if(ibr(iii,jj).eq.1)then
                eta(iii,jj)=eta(i1,jj)+(d1(i1,jj)-d1(iii,jj))/5.
                if(eta(iii,jj).lt..0) eta(iii,jj)=0.
                end if
            if(d1(iii,jj).gt.2..and.eta(iii,jj).gt..01) eta(iii,jj)=.01
              end if
            end do
            do jj=2,jgmx-1
              if(d1(iii,jj).gt..001) then
                if(d1(iii,jj-1).le..001) eta(iii,jj)=eta(iii,jj+1)
                if(d1(iii,jj+1).le..001) eta(iii,jj)=eta(iii,jj-1)
                if(d1(i3,jj).le..001) eta(iii,jj)=  &
                 ((h13(i1,jj)+h13(iii,jj))*.085+eta(i1,jj))/2.
              end if
            end do
            do jj=2,jgmx-1
              if(d1(iii,jj).gt..001) hsgg(jj)=  &
               (eta(iii,jj-1)+eta(iii,jj+1)+eta(iii,jj))/3.
            end do
            do jj=2,jgmx-1
              if(d1(iii,jj).gt..001)eta(iii,jj)=(hsgg(jj)+eta(i1,jj))/2.
            end do
            if(d1(iii,1).gt..001) eta(iii,1)=eta(iii,2)
            if(d1(iii,jgmx).gt..001) eta(iii,jgmx)=eta(iii,jgmx-1)
            iddd=0

            do jj=1,jgmx
                if(d1(i3,jj).gt..05.and.eta(iii,jj).gt.0.01) then
                eta(iii,jj)=eta(iii,jj)/2.+0.005
                end if
              if(d1(iii,jj).gt..001.and.d1(i3,jj).le.h13(iii,jj)/2.)then
                do 811 k=1,ijstruc2
                  if(imod.eq.0) then
                    IF(istruc2(k).NE.i3) GOTO 811
                    JOP=jstruc2(k)
                  end if
                  if(imod.eq.2) then
                    IF(istruc2(k).NE.IGMX-III) GOTO 811
                    JOP=JGMX-jstruc2(k)+1
                  end if
                  if(imod.eq.1) then
                    IF(jstruc2(k).NE.i3) GOTO 811
                    JOP=JGMX-istruc2(k)+1
                  end if
                  if(imod.eq.3) then
                    IF(jstruc2(k).NE.IGMX-III) GOTO 811
                    JOP=istruc2(k)
                  end if
                  if(jj.ne.jop) goto 811
                  cc1=(d1(i1,jj)+dstruc2(k)-tide-eta(i1,jj))/5.
                  cc2=(d1(i1,jj)-d1(iii,jj))/dvarx(iii)
                  if(d1(i3,jj).le..001.and.aslop(jj).eq.1.) then
                    if(cc2.lt.0.25.and.cc2.gt.0.0) then
                    aslop(jj)=.299+cc2**2.8*34
                    else
                    aslop(jj)=.99999
                    end if
                  end if
!             if(iabs(iii-190).le.1.and.jj.eq.143) then
!             write(*,*) 'ok3',iii,eta(iii,jj),eta(iii,jj+1),eta(iii,jj-1)
!             write(*,*) 'ok4',iii,eta(i1,jj),dstruc2(k),cc1
!             end if
                  eta(iii,jj)=eta(i1,jj)+cc1
!
             if(d1(i1,jj).lt..01) go to 813
             cc=h13(iii,jj)/h13(i1,jj)
             if(cc.lt.1.05) go to 813
             h13(iii,jj)=h13(iii,jj)*1.1/cc
             eta(iii,jj)=eta(iii,jj)*1.1/cc
             sj(jj)=sj(jj)/cc**2
             do l=imd,md
               do kk=1,nf
               si(jj,kk,l)=si(jj,kk,l)/cc**2
             end do
             end do
  813        continue
!
         iddd=1
         if(d1(iii,jj).ge.d1(i1,jj)*2.) eta(iii,jj)=eta(i1,jj)-cc1
                  dd=2.3*hsk(jj)
                  if(tp.gt.3.) dd=dd*(1.+1.5*tanh(tp-3.))
         if(dd.gt.h13b(jj)*1.5) dd=h13b(jj)*1.5
         if(dd.lt.h13a(jj)) dd=h13a(jj)
                  if(eta(iii,jj).gt.dd*aslop(jj)) then
                  eta(iii,jj)=(eta(iii,jj)+eta(i1,jj))/2.
!                 if(cc1/dvarx(iii).lt..25) eta(iii,jj)=eta(i1,jj)/2.
                  goto 811
                  end if
                  ! --- The 2% upswash is assumed identical to wave setup,
                  !     so the total wave runup is twice the wave setup.
                  !     write(*,*) iii,jj,k,dstruc2(k),eta(iii,jj)
                  !     write(*,*) iii,jj,h13(iii,jj),sxx(iii,jj),sxy(iii,jj)
                  cc=tide+2.*eta(iii,jj)-dstruc2(k)
                  if(cc.le..001) go to 811
                  ijb(i3,jj)=1
                  cc0=d1(i3,jj)
                  if(eta(iii,jj).gt.dstruc2(k)-tide) then
                  d1(i3,jj)=eta(iii,jj)-dstruc2(k)+tide
                  if(d1(i3,jj).lt..01) d1(i3,jj)=.01
                  if(jj-1.ge.1) then
                  dep(i3,jj)=(d1(i3,jj)+d1(iii,jj)+d1(iii,jj-1))/3.
                  else
                  dep(i3,jj)=(d1(i3,jj)+d1(iii,jj)+d1(iii,jj))/3.
                  end if
                  if(jj+1.le.jgmx) then
                  dep(i3,jj+1)=(d1(i3,jj)+d1(iii,jj)+d1(iii,jj+1))/3.
                  else
                  dep(i3,jj+1)=(d1(i3,jj)+d1(iii,jj)+d1(iii,jj))/3.
                  end if
                  else
                  d1(i3,jj)=(cc+d1(iii,jj))/2.
                  if(d1(i3,jj).lt..01) d1(i3,jj)=.01
                  if(jj-1.ge.1) then
                  dep(i3,jj)=(d1(i3,jj)+d1(iii,jj)+d1(iii,jj-1))/3.
                  else
                  dep(i3,jj)=(d1(i3,jj)+d1(iii,jj)+d1(iii,jj))/3.
                  end if
                  if(jj+1.le.jgmx) then
                  dep(i3,jj+1)=(d1(i3,jj)+d1(iii,jj)+d1(iii,jj+1))/3.
                  else
                  dep(i3,jj+1)=(d1(i3,jj)+d1(iii,jj)+d1(iii,jj))/3.
                  end if
                  go to 811
                  end if

        do iadd=i3+1,i3+30
          if(iadd.gt.igmx) go to 811
          if(cc0.eq.-2.0002) go to 811
          if(d1(iadd,jj).eq.cc0) then
          ijb(iadd,jj)=1
          if(iadd.eq.i3+1) eta(iii,jj)=eta(iii,jj)/2.
          d1(iadd,jj)=d1(i3,jj)
          dep(iadd,jj)=d1(i3,jj)
          dep(iadd,jj+1)=d1(i3,jj)
          else
          go to 811
          end if
        end do
  811           continue
              end if
            end do
          end if

           if(iddd.eq.1) then
            do jj=2,jgmx-1
            if(ijb(i3,jj).ge.4) then
            if(ijb(i3,jj+1).ge.1.and.ijb(i3,jj-1).ge.1.) then
            ijb(i3,jj)=1
            end if
            end if
            hsgg(jj)=(eta(iii,jj-1)+eta(iii,jj+1)+eta(iii,jj))/3.
            if(d1(iii,jj).gt..001) then
            if(d1(iii,jj-1).le..001.and.d1(iii,jj+1).le..001)  &
            hsgg(jj)=eta(iii,jj)
      if(d1(iii,jj-1).le..001) hsgg(jj)=(eta(iii,jj)+eta(iii,jj+1))/2.
!
      if(d1(iii,jj+1).le..001) hsgg(jj)=(eta(iii,jj)+eta(iii,jj-1))/2.
!
            end if
            end do
            hsgg(1)=hsgg(2)
            hsgg(jgmx)=hsgg(jgmx-1)
            do jj=2,jgmx-1
              dd11(jj)=(d1(i3,jj+1)+d1(i3,jj-1))/2.
              if(d1(i3,jj+1).lt..001) dd11(jj)=d1(i3,jj-1)/2.
              if(d1(i3,jj-1).lt..001) dd11(jj)=d1(i3,jj+1)/2.
            end do
            dd11(1)=dd11(2)
            dd11(jgmx)=dd11(jgmx-1)
            do jj=1,jgmx
              if(d1(iii,jj).gt..001)then
              eta(iii,jj)=hsgg(jj)
              else
              if(jj.eq.1) then
              eta(iii,1)=(hsgg(1)+hsgg(2))/2.
              elseif(jj.eq.jgmx) then
              eta(iii,jgmx)=(hsgg(jgmx)+hsgg(jgmx-1))/2.
              else
              eta(iii,jj)=(hsgg(jj+1)+hsgg(jj)+hsgg(jj-1))/3.
              end if
              end if

!             do 815 k=1,ijstruc2
!                 if(imod.eq.0) then
!                   IF(istruc2(k).NE.iii) GOTO 815
!                   JOP=jstruc2(k)
!                 end if
!                 if(imod.eq.2) then
!                   IF(istruc2(k).NE.IGMX-I1) GOTO 815
!                   JOP=JGMX-jstruc2(k)+1
!                 end if
!                 if(imod.eq.1) then
!                   IF(jstruc2(k).NE.iii) GOTO 815
!                   JOP=JGMX-istruc2(k)+1
!                 end if
!                 if(imod.eq.3) then
!                   IF(jstruc2(k).NE.IGMX-I1) GOTO 815
!                   JOP=istruc2(k)
!                 end if
!                 if(jj.ne.jop) goto 815
!                 cc=dstruc2(k)-tide
!             if(max(2.*eta(iii,jj),eta(iii,jj)+h13(iii,jj)/2.).lt.cc) then
!             eta(iii,jj)=0.
!             d1(iii,jj)=0.
!             end if
! 815         continue
              jja1=jj+1
              jjc1=jj-1
              if(jj.eq.1) jjc1=2
              if(jj.eq.jgmx) jja1=jgmx-1
              if(d1(i3,jja1).gt..001.or.d1(i3,jjc1).gt..001) then
              if(d1(i3,jj).lt.dd11(jj)) then
              d1(i3,jj)=dd11(jj)
              dep(i3,jj)=(dd11(jj)+d1(iii,jj))/2.
              dep(i3,jj+1)=dep(i3,jj)
              end if
              end if
            end do
            end if

          if(iii.eq.igmx) then
            do jj=1,jgmx
              if(d1(igmx,jj).gt.0001) eta(igmx,jj)=eta(i1,jj)
              if(d1(1,jj).gt..0001) eta(1,jj)=eta(2,jj)
            end do
          end if
        end if
      end if

      do 812 k=1,ijstruc4
        if(imod.eq.0) then
          IF(istruc4(k).NE.i3) GOTO 812
          JOP=jstruc4(k)
        end if
        if(imod.eq.2) then
          IF(istruc4(k).NE.IGMX-III) GOTO 812
          JOP=JGMX-jstruc4(k)+1
        end if
        if(imod.eq.1) then
          IF(jstruc4(k).NE.i3) GOTO 812
          JOP=JGMX-istruc4(k)+1
        end if
        if(imod.eq.3) then
          IF(jstruc4(k).NE.IGMX-III) GOTO 812
          JOP=istruc4(k)
        end if

          if(d1(i3,jop).eq..01) then
          cc=dstruc4(k)+eta(iii,jop)
        if(k4(k).eq.4)d1(i3,jop)=cc+h13(iii,jop)*.89
        if(k4(k).eq.5)d1(i3,jop)=cc+h13(iii,jop)*.73
        if(d1(i3,jop).lt..01) then
        d1(i3,jop)=.01
        c1=abs(cc)/(h13(iii,jop)+.01)
        if(k4(k).eq.4.and.c1.gt..5) d1(i3,jop)=h13(iii,jop)/5./c1
        if(k4(k).eq.5.and.c1.gt..5) d1(i3,jop)=h13(iii,jop)/20./c1
        end if
          else
          cc=1.
          c1=h13(iii,jop)/2.
          if(d1(i3,jop).gt.c1) cc=c1/d1(i3,jop)
        d1(i3,jop)=d1(i3,jop)+c1*cc
!       the following line may not be the best
!       if(d1(i3+1,jop).le.d1(i3,jop))d1(i3+1,jop)=(d1(i3,jop)+d1(i3+1,jop))/2.
        if(d1(i3+1,jop).le.d1(i3,jop).and.kstruc4(k).eq.0) &
        d1(i3+1,jop)=d1(i3,jop)
          end if

!       if(dep(i3,jop).lt.d1(i3,jop)) then
        dep(i3,jop)=(d1(i3,jop)+d1(i3+1,jop))/2.
        if(dep(i3+1,jop).lt.dep(i3,jop)) dep(i3+1,jop)=dep(i3,jop)
!       end if
!       if(dep(i3,jop+1).lt.d1(i3,jop)) then
!       dep(i3,jop+1)=(d1(i3,jop)+d1(i3+1,jop))/2.
        dep(i3,jop+1)=dep(i3,jop)
!       if(dep(i3+1,jop+1).lt.dep(i3,jop+1)) dep(i3+1,jop+1)=dep(i3,jop+1)
        if(dep(i3+1,jop+1).lt.dep(i3,jop)) dep(i3+1,jop+1)=dep(i3,jop)
!       end if

        do iadd=i3+1,i3+30
          if(iadd.gt.igmx) go to 812
          if(d1(iadd,jop).eq..01) then
          if(k4(k).eq.4)d1(iadd,jop)=h13(iii,jop)*.38
          if(k4(k).eq.5)d1(iadd,jop)=h13(iii,jop)*.27
          if(d1(iadd,jop).gt.d1(i3,jop)) d1(iadd,jop)=d1(i3,jop)
          else
            go to 812
          end if
        end do
  812 continue

      if(ibreak.eq.2) then
        call dissip1_inline(iii)
      end if

      if(ibreak.eq.3) then
        call dissip_inline(iii)
      end if
339   continue

!     if(leap.gt.1) THEN 
        simax=si
        sjj=sj
!     end if
      if(h13(ii,j1).lt.1.0E-15.and.ii.ge.igmx/2) goto 50
      IF(II.GE.IGMX) GOTO 50
! -- repetition of calculation for next ii location

      leap=1
      dxx=0.
      if(island.eq.1) go to 337

      do 336 jj=2,jgmx-1
      if(ijb(ii,jj).eq.-1) go to 336
      if(d1(ii,jj).lt..01.or.ibr(ii,jj).eq.1) then
      island=1
      go to 337
      end if
336   continue

      do 333 ic=ii+2,igmx
      if(ic.ge.ismall) go to 331
      do 335 jj=2,jgmx-1
      cc=d1(ic,jj)
      if(ijb(ic,jj).eq.-1) cc=10.
      if(cc.lt.25.) go to 331
335   continue
      leap=leap+1
333   continue
331   continue
      if(leap.gt.10.and.ws.ge..1) leap=10
      if(ws.ge.10.) then
      if(leap.gt.3) leap=3
      if(deps(ii).lt.30.) leap=2
      if(deps(ii).lt.20.) leap=1
      if(ii.lt.20) leap=1
      end if

337   continue
      if(ii+leap.gt.igmx) leap=igmx-ii

      do kk=ii,ii+leap-1
      dxx=dxx+dvarx(kk)
      if(dxx.gt.250..and.leap.ge.2) then
      dxx=dxx-dvarx(kk)
      leap=kk-ii
      if(leap.lt.1) then
      leap=1
      dxx=dxx+dvarx(kk)
      end if
      goto 601
      end if
      end do
601   continue

      II=II+leap

      if(leap.gt.1) then
      do 508 jj=1,jgmx
      t13(ii-1,jj)=t13(ii-leap,jj)
508   continue
      end if

      GOTO 60

   50 IF(ICK4.NE.1.or.KRF.EQ.0) GOTO 1000
!  ** check whether reflection cal is needed or not

!==================================================
!  preparation of backward marching
!    invert of i number, change of cell conditions
!    invert of current directions
!    ibk=1 means no wave breaking considred
!==================================================
      IBK=1
      IBACK=1

!     if(ick4.eq.1) then
      write(*,*) ' '
      write(*,*) '*** Backward Reflection Calculation ***'
      write(*,*) ' '
!     end if

!     open(unit=129,file='depp.dat',status='unknown')
!    do 229 j=jmax,1,-1
!     write(129,126) (dep(i,j),i=1,imax)
!229  continue
!     close(unit=129)

! -- II is original grid cell number
!    INV is new grid cell number in reflection direction
! -- iteration for ii
      ni2=igmx/2
      do ii=1,ni2
      inv=imax-ii
      cc=dvarx(inv)
      dvarx(inv)=dvarx(ii)
      dvarx(ii)=cc
      end do

      DO 110 II=1,IGMX
      INV=IMAX-II
      I1=II+1

! -- iteration for jj
      DO 110 JJ=1,JMAX
      
! -- for all points (KRF) defined as reflection points,
!    check their positions of (KI,KJ).  If (KI,KJ) is
!    equal to present position (II,JJ), IJTR is set 1.
! -- check of reflection points
      IJTR=0
      DO 310 K=1,KRF
      KI=KCR(1,K)
      IF(KI.NE.II) GOTO 310
      KJ=KCR(2,K)
      IF(KJ.NE.JJ) GOTO 310
      IJTR=1
  310 CONTINUE

! -- ijtr is set to 1 if reflection points are there
! -- quantities are defined on new grid cell number INV
      if(ii.gt.ni2) go to 311
      cc=DEP(INV,JJ)
      DEP(INV,JJ)=DEP(II,JJ)
      DEP(II,JJ)=cc

      cc=D1(INV,JJ)
      D1(INV,JJ)=D1(II,JJ)
      D1(II,JJ)=cc

      cc=u(inv,jj)
      u(inv,jj)=-u(ii,jj)
      u(ii,jj)=-cc

      cc=v(inv,jj)
      v(inv,jj)=v(ii,jj)
      v(ii,jj)=cc

      cc=u1(inv,jj)
      u1(inv,jj)=-u1(ii,jj)
      u1(ii,jj)=-cc

      cc=v1(inv,jj)
      v1(inv,jj)=v1(ii,jj)
      v1(ii,jj)=cc
311   continue
!     IF(II.GE.IMAX) GOTO 110
!     IF(JJ.GE.JMAX) GOTO 110
      IJBV=IJBP(II,JJ)
      IJBV1=IJBP(I1,JJ)
      IJB1=IJBV+1
      if(ijb1.eq.0) ijb1=1
      GOTO (1,2,3,4,5,6,7,8,9,10),IJB1

!  ** which means 0->0
    1 IJB(INV,JJ)=IJBV
      GOTO 110

!  ** which means 1->1, (1,0): 1->8, (1,1)~(1,5): 1->1
    2 IF(IJTR.EQ.1) GOTO 1
      IF(IJBV1.LE.0) GOTO 21
      IF(IJBV1.LT.6) GOTO 1
   21 IJB(INV,JJ)=8
      GOTO 110

!  ** which means 2->3 (energy flow out)
    3 IJB(INV,JJ)=3
      GOTO 110

!  ** which means 3->3, (3,0): 3->9, (3,1)~(3,5): 3->3
    4 IF(IJTR.EQ.1) GOTO 1
      IF(IJBV1.LE.0) GOTO 41
      IF(IJBV1.LT.6) GOTO 1
   41 IJB(INV,JJ)=9
      GOTO 110

!  ** which means 4->4, (4,0): 4->6, (4,1)~(4,5): 4->4
    5 IF(IJTR.EQ.1) GOTO 1
      IF(IJBV1.LE.0) GOTO 51
      IF(IJBV1.LT.6) GOTO 1
   51 IJB(INV,JJ)=6
      GOTO 110

!  ** which means 5->5, (5,0): 5->7, (5,1)~(5,5): 5->5
    6 IF(IJTR.EQ.1) GOTO 1
      IF(IJBV1.LE.0) GOTO 61
      IF(IJBV1.LT.6) GOTO 1
   61 IJB(INV,JJ)=7
      GOTO 110

!  ** which means 6->4
    7 IJB(INV,JJ)=4
      GOTO 110

!  ** which means 7->5
    8 IJB(INV,JJ)=5
      GOTO 110

!  ** which means 8->1
    9 IJB(INV,JJ)=1
      GOTO 110

!  ** which means 9->3
   10 IJB(INV,JJ)=3

  110 CONTINUE
! -- end of iteration for jj
! -- end of iteration for ii
!
      do 223 i=1,igmx
      if(ijb(i,1).eq.3) ijb(i,1)=2
      if(ijb(i,jgmx).eq.3) ijb(i,jgmx)=2
 223  continue

!     open(unit=22,file='ijbv.dat',status='unknown')
!    do 222 j=jgmx,1,-1
!     write(22,26) (ijb(i,j),i=igmx,1,-1)
!26   format((200i1))
!222  continue
!     close(unit=22)

      do j=1,jmax
      dep(imax,j)=dep(igmx,j)
      end do
!     open(unit=122,file='depv.dat',status='unknown')
!    do 225 j=jmax,1,-1
!     write(122,126) (dep(i,j),i=imax,1,-1)
!126  format((200f5.1))
!225  continue
!     close(unit=122)
!
!     open(unit=222,file='d1v.dat',status='unknown')
!    do 227 j=jgmx,1,-1
!     write(222,226) (d1(i,j),i=igmx,1,-1)
!226  format((200f5.1))
!227  continue
!     close(unit=222)

      IF(ICK3.NE.1) GOTO 112

! -- reflection in y direction
! -- i grid number and structure angle are changed
      DO 111 K=1,KRMX
      KR(1,K)=IMAX-KR(1,K)
      yangl(k)=-yangl(k)
  111 CONTINUE
  112 CONTINUE

!===========================================================
!  start of backward marching
!    ii: new coordinate
!   inv: original coordinate
!  coordinates are changed compared to previous definition
!===========================================================
      II=1
      si=0.
      IRC=0
! -- ick4=0 means no reflection in x direction (when backward marching)
      ICK4=0
      wskeep=ws
      ws=0.
!
  150 INV=IMAX-II
      if (mod(ii-1,25).eq.0) write(*,9998) ii
 9998 format ('Calc. Backward Reflection at ',i0)

! -- check whether reflection points in x direction are there or not
!     for new number II (that is, original number INV),
!     if it is reflection position,
! -- position JR(KK) and refl coefficient RKK(KK) are stored.
      KK=0
      DO 120 K=1,KRF
      IF(KCR(1,K).NE.INV) GOTO 120
      KK=KK+1
      KRC(KK)=K
      JR(KK)=KCR(2,K)
      RKK(KK)=RKR(K)
      angl(kk)=xangl(k)
  120 CONTINUE
      IF(KK.NE.0) GOTO 130
      IF(IRC.NE.0) GOTO 140
! -- when no reflection points
      II=II+1
      GOTO 150

! -- when reflection points exist, input wave spectrum is given
  130 continue
      DO 160 JJ=1,JGMX
! -- jj is reflecton point?
      DO 170 K=1,KK
      IF(JR(K).EQ.JJ) GOTO 180
  170 CONTINUE
      GOTO 160
! -- setup of input energy at reflection points
  180 continue
!     sj(jj)=0.
      DO 210 NN=1,NF
      alpha=angl(k)*rad
      ishift=int(2.0*alpha/dth)
      amari=2.0*alpha-float(ishift)*dth
      if(amari.ge.0.5) ishift=ishift+1  
      DO 210 MM=imd,MD
      mmd=mm-ishift
      if(mmd.ge.imd.and.mmd.le.md) then
      SI(JJ,NN,MMd)=SR(KRC(K),NN,MM)*RKK(K)**2
!     sj(jj)=sj(jj)+si(jj,nn,mmd)
      end if
  210 CONTINUE
! -- end do for number of frequency
  160 CONTINUE
! -- end do for j grid points

  140 IRC=1

      CALL CALC_inline(II,IKAP)
!------------------------------------------------------------------------
!  by CALC, si(jj,nn,mm),sj(jj) are obtained
!------------------------------------------------------------------------

!     if(itrack.eq.1) then
!     if(ii.eq.53.and.jj.le.2) then
!     write(*,*) 'jj',((wk2(jj,nn,mm),nn=1,nf),mm=imd,md)
!     write(*,*) 'jj',((si(jj,nn,mm),nn=1,nf),mm=imd,md)
!     end if
!     sj(1)=sj(2)
!     end if

      DO 220 JJ=1,JGMX
      cc=(D1(II,JJ)**2-H13(II,JJ)**2)/cflat**2
!
      if(cc.lt.0.) cc=0.
      if(SJ(JJ).GT.cc) THEN
         if(sj(jj).gt.1.e-15) then
         do l=imd,md
         do k=1,nf
         si(jj,k,l)=si(jj,k,l)*cc/SJ(JJ)
         end do
         end do
         end if
!
! -- Here, the reflected wave energy is si, directly
!    obtained from the forward propagation (so, it
!    can already break); however, the saved reflected
!    wave sj may not have the same energy because the
!    total energy can violate the depth limitation.
!    The question is: do we use the proper si to
!    begin the reflection?
!
        SJ(JJ)=cc
      end if
!
      SS=SJ(JJ)
!
      IF(SS.LT.1.0E-15) GOTO 220
      H13R(INV,JJ)=cflat*SQRT(SS)
      T13R(INV,JJ)=1./FJF(JJ)
!
      sumx=0.
      sumy=0.
      do 5077 mm=imd,md
      scp(jj,mm)=0.
      do 5033 nn=1,nf
      scp(jj,mm)=scp(jj,mm)+si(jj,nn,mm)
 5033 continue
      sumx=sumx+scp(jj,mm)*cosa(mm)
      sumy=sumy+scp(jj,mm)*sina(mm)
 5077 continue
      DEG=atan2(sumy,sumx)/rad
!
      PD=180.0
      IF(DEG.LT.0) PD=-180.0
      DMNR(INV,JJ)=PD-DEG
  220 CONTINUE
!
      do 233 k=1,kout
      if(IJSP(1,k).eq.inv) then
!
         jop=ijsp(2,k)
         do mm=imd,md
         mm1=md+mm-imd+1
         do nn=1,nf
         sop(k,nn,mm1)=si(jop,nn,mm)/df(nn)/dth
         end do
         end do
!
      end if
  233 continue
!
      if(irs.ge.1) call sxycalc_inline(inv)
!
      IF(II.GE.IGMX) GOTO 250
      II=II+1
      GOTO 150
! -- repetition of backward calculation for next ii location
!
  250 ICK4=1
      IBACK=0
!
      ws=wskeep
!
      do ii=1,ni2
       inv=imax-ii
       cc=dvarx(inv)
       dvarx(inv)=dvarx(ii)
       dvarx(ii)=cc
       do JJ=1,JGMX
        cc=D1(INV,JJ)
        D1(INV,JJ)=D1(II,JJ)
        D1(II,JJ)=cc
       end do
      end do
!
 1000 DEALLOCATE (KRC,JR,RKK,angl,hsgg,h13a,h13b)
      RETURN     
      END

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!  main of calculation by using sub. sgma, veloc, setab, gsm
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE CALC_inline(II,IKAP)
!------------------------------------------------------------
!  independent variables are x,y,q(angle)
!  intrinsic frequency is obtained form dispersion relation
!------------------------------------------------------------
      USE GLOBAL_INLINE, ONLY: NPF,MPD,IPMX,JPMX,KOMX,NOMX,IGPX,JGPX
      COMMON /VPAI/PAI2,PAI,HPAI,RAD,akap,imod,iprp,island,imd,iprpp     &
                   ,nonln,igrav,isolv,ixmdf,iproc,imud,iwnd,depmin0
      COMMON /DATA/NF,MD,IMAX,JMAX,IGMX,JGMX,JCPB,JCPE,JCB,JCE,NFF,MDD
      COMMON /DATC/TP,PRD,IBND,VL,TIDE,KPEAK,IBACK,NBLOCK,IWIND,WSMAG
      COMMON /DATB/ICK3,ICK4,ICHOICE,HS0,wd,ws,ph0,wdd(mpd),aslop(ipmx)
      COMMON /DATD/DX,DY,DXX,DMESH,DTH,kdate,idate,depmax,depmax0
      COMMON /FDS/FCN(NPF),DCM(MPD),DSFD(NPF,MPD)
      COMMON /IJBP/IJB(IGPX,JGPX),DEP(IPMX,JPMX),DEPS(IPMX),DBIG(IPMX)
      COMMON /SPECA/SCP(JGPX,MPD),SI(JGPX,NPF,MPD)
      COMMON /SPECC/SJ(JGPX),SJF(JGPX,NPF),FJF(JGPX)
      COMMON /WAVI/H13(IGPX,JGPX),T13(IGPX,JGPX),DMN(IGPX,JGPX)
      COMMON /DFCN/DF(NPF),FFCN(NPF),DFINP(NPF),PL1E(NPF),PL1T(NPF)
      common /brek/depm(jgpx),dmnj(jgpx),slf(jgpx),wlmn(jgpx),cmn(jgpx)  &
                   ,sigm(jgpx),IWVBK,ibk3
      common /sgma01/sgma0(jgpx,mpd),sgma1(jgpx,mpd)
      common /ccg/cwk(jpmx,2,mpd),cgk(jpmx,2,mpd)
      common /rsm/ cosa(mpd),sina(mpd),d1(ipmx,jpmx)
!      dimension scpp(jgpx,mpd),ssmax(jgpx)
      REAL, ALLOCATABLE :: scpp(:,:),ssmax(:)
      ALLOCATE (scpp(jgpx,mpd),ssmax(jgpx))

! -- small lower limit value
      eliml=1.0e-10

! -- initialization
      SJ=0.

! -- IFC shows a frequency components **
      IFC=1
  120 FC=FCN(IFC)
      DFC=DF(IFC)

! -- initialization
      SCP=0.
      sgma0=-1.
      sgma1=-1.

! -- find out of start and end of j (jb and je) **
!
      do k=1,jgmx
      if(ijb(ii,k).ge.8) then
      if(ijb(ii-1,k).eq.2) ijb(ii,k)=2
      if(ijb(ii-1,k).eq.1) ijb(ii,k)=1
      end if
      end do
!
      if(ijb(ii,1).eq.1) ijb(ii,1)=2
      if(ijb(ii,jgmx).eq.1) ijb(ii,jgmx)=2
      if(ijb(ii,1).eq.0) then
      do k=2,jgmx
      if(ijb(ii,k).ne.0) then
      ijb(ii,k)=4
      go to 845
      end if
      end do
      end if
845   continue

      if(ijb(ii,jgmx).eq.0) then
      do k=jgmx-1,1,-1
      if(ijb(ii,k).ne.0) then
      ijb(ii,k)=5
      go to 846
      end if
      end do
      end if
846   continue

      JCP=0
      JJ=1
   30 JJB=IJB(II,JJ)
      IF(JJB.GT.0) GOTO 20
   50 JJ=JJ+1
      IF(JJ.GT.JGMX) GOTO 90
      GOTO 30
   20 IF(JCP.NE.0) GOTO 40
      JB=JJ
      JCP=1
      GOTO 50
   40 IF(JJB.EQ.1) GOTO 50
      IF(JJB.EQ.8) GOTO 50
      JE=JJ

!--------------------------------------------------
!  intrinsic angular wave frequency 
!    sgma0(ii-1,j), sgma1(ii,j)
!  if there is no solution of dispersion relation 
!    and Cg+U<0, sigma=-1.0 is set
!  sigma is used to define wave action,
!    and to check whether wave exists or not
!--------------------------------------------------
      call sgma_inline(ii,jb,je,ifc)
!      write(*,6666) fc,ii
! 6666 format(' Freq.',f10.3,/,'  sub.sgma at I='i3,' is OK')

! -- input wave action
      DO 42 MM=imd,MD
      DO 42 JJSI=JB,JE
      if(sgma0(jjsi,mm).lt.eliml) then
      scp(jjsi,mm)=0.0001
      else

!- si(jj,ifc,mm) is input spectrum
      scp(jjsi,mm)=si(jjsi,ifc,mm)/sgma0(jjsi,mm)
!- si(jj,ifc,mm) is input spectrum

      end if
   42 CONTINUE

      if(iback.eq.1.or.akap.eq.0.) go to 488
      aakap=10.
      if(ii.gt.1) then
      amd2=float(mdd)/2.
      if(ichoice.eq.2) amd2=float(md)/2.
      do 43 jjsi=jb,je
      if(jjsi.eq.1.or.jjsi.eq.jgmx) go to 43
      if(ijb(ii,jjsi).eq.1.and.ijb(ii,jjsi-1).ge.6) then
      large=imd
      slarge=0.
      do 44 mm=imd,md
      if(scp(jjsi,mm).gt.slarge) then
      slarge=scp(jjsi,mm)
      large=mm
      end if
   44 continue
      if(large.eq.imd) goto 48
      do 45 mm=imd,large-1
      alarge=float(large-1)
      cc=alarge*(2.*alarge+1.)*(alarge+1.)
      scp(jjsi,mm)=max(slarge*float(large-mm)/cc*15.*aakap,scp(jjsi,mm))
!
   45 continue
      md4=min0(md,md-(md-large)/(2+nint(aakap**3)/9))
      do 455 mm=large,md4
      alarge=float(md4-large+1)
      cc=alarge*(2.*alarge+1.)*(alarge+1.)
      scp(jjsi,mm)=max(slarge*float(mm-large+1)/cc*100.  &
                  /(1.+abs(float(large)-amd2)),scp(jjsi,mm))
  455 continue
      end if
   48 continue
      if(ijb(ii,jjsi).eq.1.and.ijb(ii,jjsi+1).ge.6) then
      large=imd
      slarge=0.
      do 46 mm=imd,md
        if(scp(jjsi,mm).gt.slarge) then
          slarge=scp(jjsi,mm)
          large=mm
        end if
   46 continue
      if(large.eq.md) goto 49
      do 47 mm=large+1,md
      alarge=float(md-large)
      cc=alarge*(2.*alarge+1.)*(alarge+1.)
      scp(jjsi,mm)=max(slarge*float(mm-large)/cc*15.*aakap,scp(jjsi,mm))
!
   47 continue
      md1=max0(imd,imd+(large-imd)/(2+nint(aakap**3)/9))
      do 477 mm=md1,large
      alarge=float(large-md1+1)
      cc=alarge*(2.*alarge+1.)*(alarge+1.)
      scp(jjsi,mm)=max(slarge*float(large-mm+1)/cc*100.  &
                  /(1.+abs(float(large)-amd2)),scp(jjsi,mm))
  477 continue
      end if
   49 continue
   43 continue
      end if
  488 continue
!
      if(ii.eq.1) goto 90
!
!       if(ws.ge..1) then
        do jjsi=jb,je
        ssmax(jjsi)=1.0e-6
        do mm=imd,md
        ssmax(jjsi)=ssmax(jjsi)+scp(jjsi,mm)
        scpp(jjsi,mm)=scp(jjsi,mm)
        end do
        end do
!       end if

! -- c and cg in current field are calculated
      CALL VELOC_inline(II,JB,JE,FC)
! -- cases of sg, cw, ak=-1 are contained, GSM matrix is modified after
!      write(*,6667) ii
! 6667 format('  sub.veloc at I='i3,' is OK')

! -- matrix to be solved is made
! -- matrix size of nmx=md*jmsh is defined in subroutine setab 
      CALL SETAB_inline(II,JB,JE,NMX,FC,DFC,IFC,IKAP)
!      write(*,6668) ii
! 6668 format('  sub.setab at I='i3,' is OK')

! -- solution of matrix
      MARK=0
      if(isolv.eq.0) CALL GSR_inline(II,JB,JE,NMX,MARK)
      if(isolv.eq.1) CALL ADI_inline(II,JB,JE,NMX,MARK)
      if(isolv.eq.2) CALL GSM_inline(II,JB,JE,NMX,MARK)

!     if(mark.ge.1) then
      if(ws.ge..1.or.mark.ge.1) then
        do 1101 jjsi=jb,je
        if(ssmax(jjsi).lt.1.0e-5) goto 1101
        sum=0.
        do mm=imd,md
        sum=sum+scp(jjsi,mm)
        end do
        cc=sum/ssmax(jjsi)
        cc1=1.+0.5*tanh(d1(ii,jjsi)/5.)
        if(cc.gt.cc1) then
        do mm=imd,md
        scp(jjsi,mm)=scpp(jjsi,mm)*cc1
        end do
        elseif(cc.lt.0) then
        do mm=imd,md
        scp(jjsi,mm)=scpp(jjsi,mm)
        end do
        end if
 1101   continue
      end if
!
      scpp=scp
      do jjsi=jb+1,je-1
        c1=d1(ii,jjsi-1)
        c2=d1(ii,jjsi)+.0001
        c3=d1(ii,jjsi+1)
        cc=c1+c2+c3
        do mm=imd,md
        scp(jjsi,mm)=(c1*scpp(jjsi-1,mm)+c2*scpp(jjsi,mm)+c3*scpp(jjsi+1,mm))/cc
        end do
      end do
!
      if(je-jb.ge.1) then
        do mm=imd,md
        c1=d1(ii,jb)
        c2=d1(ii,jb+1)+.0001
        scp(jb,mm)=(c1*scpp(jb,mm)+c2*scpp(jb+1,mm))/(c1+c2)
        c1=d1(ii,je)
        c2=d1(ii,je-1)+.0001
        scp(je,mm)=(c1*scpp(je,mm)+c2*scpp(je-1,mm))/(c1+c2)
        end do
      end if
!
!      write(*,6669) ii
! 6669 format('  sub.gsm at I='i3,' is OK')

      IF(JJ.GE.JGMX) GOTO 90
! -- next find out of jb and je
      JCP=0
      GOTO 50

!==========================================================
!  after calculation for all jj at ii, results are stored
!==========================================================
!----------------------------------------------
!  check of reflection points in -x direction
!----------------------------------------------
   90 DO 60 JBE=1,JGMX

! -- transform of wave action to spectral energy
      do 81 mm=imd,md
      scp(jbe,mm)=scp(jbe,mm)*sgma1(jbe,mm)
      if(sgma1(jbe,mm).lt.0.0) scp(jbe,mm)=0.0
   81 continue

! -- mb is begining point and me is ending point
      mb=0
      me=0
      do 82 mm=imd,md-1
      if(sgma1(jbe,mm)*sgma1(jbe,mm+1).gt.0.0) go to 82
      if(sgma1(jbe,mm)*sgma1(jbe,mm+1).lt.0.0.and.sgma1(jbe,mm+1).lt.0.0) mb=mm+1 
      if(sgma1(jbe,mm)*sgma1(jbe,mm+1).lt.0.0.and.sgma1(jbe,mm+1).gt.0.0) me=mm
   82 continue

      if(mb.eq.0.and.me.eq.0) then
      go to 700

      else if(mb.eq.0.and.me.ne.0) then
      mb=1
      scp(jbe,mb)=0.0
!     dscp=(scp(jbe,me+1)-0.0)/float(me)
      dscp=scp(jbe,me+1)/float(me)
      do 83 mm=mb+1,me
      scp(jbe,mm)=scp(jbe,mm-1)+dscp
   83 continue

      else if(mb.ne.0.and.me.eq.0) then
      me=md
!     dscp=(0.0-scp(jbe,mb-1))/float(me-(mb-1))
      dscp=scp(jbe,mb-1)/float(mb-me-1)
      do 84 mm=mb,me
      scp(jbe,mm)=scp(jbe,mm-1)+dscp
   84 continue

      else if(mb.ne.0.and.me.ne.0) then
!     dscp=(scp(jbe,me+1)-scp(jbe,mb-1))/float((me+1)-(mb-1))
      dscp=(scp(jbe,me+1)-scp(jbe,mb-1))/float(me-mb+2)
      do 85 mm=mb,me
      scp(jbe,mm)=scp(jbe,mm-1)+dscp
   85 continue

      else
      go to 700
      end if

! -- calculated results
  700 continue
      DO 609 MM=imd,MD
      DC=DCM(MM)
      if(scp(jbe,mm).lt.0.) scp(jbe,mm)=0.
      sjm=scp(jbe,mm)
      SJ(JBE)=SJ(JBE)+SJM
      SI(JBE,IFC,MM)=SJM
  609 CONTINUE
   60 CONTINUE
! -- next frequency component
      IF(IFC.GE.NF) GOTO 110
      IFC=IFC+1
      GOTO 120

!--------------------------------------------------------
!  at this stage si(jbe,nn,mm) is obtained
!  revison of spectral energy at no wave solutioin 
!  for frequency components
!    1. integration of si(jbe,nn,mm) with respect to mm
!    2. spectral form of f^-5 is assumed
!    3. additional energy part is added to sj, sjf, sjd
!    4. note --> si(jbe,nn,mm) is integreted energy
!--------------------------------------------------------

  110 DEALLOCATE (scpp,ssmax)
      RETURN
      END


!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!   subroutine for intrinsic angular frequency
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      subroutine sgma_inline(ii,jb,je,nn)
!--------------------------------------
!  sgma0(j,m) is at ii-1
!  sgma1(j,m) is at ii
!    for determination of wave action
!    and for check of wave existing
!--------------------------------------
      USE GLOBAL_INLINE, ONLY: NPF,MPD,IPMX,JPMX,igpx,jgpx
      common /vpai/pai2,pai,hpai,rad,akap,imod,iprp,island,imd,iprpp     &
                   ,nonln,igrav,isolv,ixmdf,iproc,imud,iwnd,depmin0
      common /data/nf,md,imax,jmax,igmx,jgmx,jcpb,jcpe,jcb,JCE,NFF,MDD
      common /ijbp/ijb(igpx,jgpx),dep(ipmx,jpmx),deps(ipmx),DBIG(IPMX)
      common /rsm/ cosa(mpd),sina(mpd),d1(ipmx,jpmx)
      common /uvp/u(ipmx,jpmx),v(ipmx,jpmx)
      common /uv1/u1(ipmx,jpmx),v1(ipmx,jpmx)
      common /sgma01/sgma0(jgpx,mpd),sgma1(jgpx,mpd)
      common /fds/fcn(npf),dcm(mpd),dsfd(npf,mpd)
      common /rsw/ wk2(jpmx,npf,mpd),cgp(ipmx,jpmx)
      om=pai2*fcn(nn)
!
! -- The following ii and ii-1 routines are better than original
!    because d1, u1, v1 are already cell center values.
!
! -- at ii location
      do 1 jj=jb,je
      do 1 mm=imd,md
      depm=d1(ii,jj)
      if(depm.lt..01) depm=.01
      um=u1(ii,jj)
      vm=v1(ii,jj)
      call wccg_inline(dcm(mm),depm,um,vm,om,cw,cg,sig,akk)
      sgma1(jj,mm)=sig
      wk2(jj,nn,mm)=akk
   1  continue
!
! -- at (ii-1) location
      if(ii.eq.1) then
      do 2 mm=imd,md
      do 2 jj=jb,je
      sgma0(jj,mm)=sgma1(jj,mm)
    2 continue
      else
      i1=ii-1
      do 20 jj=jb,je
      depm=d1(i1,jj)
      if(depm.lt..01) depm=.01
      um=u1(i1,jj)
      vm=v1(i1,jj)
      do 3 mm=imd,md
      call wccg_inline(dcm(mm),depm,um,vm,om,cw,cg,sig,akk)
      sgma0(jj,mm)=sig
    3 continue
   20 continue
      end if

      return
      end subroutine


!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!  subroutine for C and Cg in wave-current field
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      subroutine veloc_inline(ii,jb,je,fc)
!---------------------------------------------------
!  variables needed for determination of
!   wave breaking dissipation are also calculated
!---------------------------------------------------
      USE GLOBAL_INLINE, ONLY: NPF,MPD,IPMX,JPMX,igpx,jgpx
      common /data/nf,md,imax,jmax,igmx,jgmx,jcpb,jcpe,jcb,JCE,NFF,MDD
      COMMON /DATB/ICK3,ICK4,ICHOICE,HS0,wd,ws,ph0,wdd(mpd),aslop(ipmx)
      COMMON /DATC/TP,PRD,IBND,VL,TIDE,KPEAK,IBACK,NBLOCK,IWIND,WSMAG
      COMMON /DVAR/DVARX(IGPX),DVARY(JGPX),ETA(IGPX,JGPX),HSK(JGPX)
      common /fds/fcn(npf),dcm(mpd),dsfd(npf,mpd)
      common /vpai/pai2,pai,hpai,rad,akap,imod,iprp,island,imd,iprpp     &
                   ,nonln,igrav,isolv,ixmdf,iproc,imud,iwnd,depmin0
      common /DATD/DX,DY,DXX,DMESH,DTH,kdate,idate,depmax,depmax0
      common /ijbp/ijb(igpx,jgpx),dep(ipmx,jpmx),deps(ipmx),DBIG(IPMX)
      common /rsm/ cosa(mpd),sina(mpd),d1(ipmx,jpmx)
      common /ccg/cwk(jpmx,2,mpd),cgk(jpmx,2,mpd)
      COMMON /BREK/DEPM(JGPX),DMNJ(JGPX),SLF(JGPX),wlmn(JGPX),cmn(jgpx)  &
                   ,sigm(jgpx),IWVBK,ibk3
      common /uvp/u(ipmx,jpmx),v(ipmx,jpmx)
      common /uv1/u1(ipmx,jpmx),v1(ipmx,jpmx)
      common /wavi/h13(igpx,jgpx),t13(igpx,jgpx),dmn(igpx,jgpx)

      g=9.806
      jje=je-jb+1
! -- om is angular frequency of component wave
      om=pai2*fc

!----------------------------------------------------
!    11 -> lower-left
!    12 -> lower-right
!    21 -> upper-left
!    22 -> upper-right
!     m -> center
!    for example of cwk(j,1,m) and cwk(j,2,m),
!   '1' is left hand, '2' right hand side position
!----------------------------------------------------
      i1=ii+1
      dep11=dep(ii,jb)
      if(dep11.lt..01) dep11=.01
      u11=u(ii,jb)
      v11=v(ii,jb)
      do 100 mm=imd,md
      call wccg_inline(dcm(mm),dep11,u11,v11,om,cw,cg,sig,akk)
      cwk(1,1,mm)=cw
      cgk(1,1,mm)=cg
  100 continue
      dep12=dep(i1,jb)
      if(dep12.lt..01) dep12=.01
      u12=u(i1,jb)
      v12=v(i1,jb)
      do 110 mm=imd,md
      call wccg_inline(dcm(mm),dep12,u12,v12,om,cw,cg,sig,akk)
      cwk(1,2,mm)=cw
      cgk(1,2,mm)=cg
  110 continue

! -- repetition for j
      do 10 jj=1,jje
      jjb=jj+jb
      jj1=jj+1
      dep21=dep(ii,jjb)
      if(dep21.lt..01) dep21=.01
      u21=u(ii,jjb)
      v21=v(ii,jjb)
      do 120 mm=imd,md
      call wccg_inline(dcm(mm),dep21,u21,v21,om,cw,cg,sig,akk)
      cwk(jj1,1,mm)=cw
      cgk(jj1,1,mm)=cg
  120 continue
      dep22=dep(i1,jjb)
      if(dep22.lt..01) dep22=.01
      u22=u(i1,jjb)
      v22=v(i1,jjb)
      do 130 mm=imd,md
      call wccg_inline(dcm(mm),dep22,u22,v22,om,cw,cg,sig,akk)
      cwk(jj1,2,mm)=cw
      cgk(jj1,2,mm)=cg
  130 continue

!------------------------------------------------------------
!  mean values bellow are used in wave breaking dissipation
!    omt3 is angular frequency corresponding to T3
!    dmnj(jjb) is mean direction at (ii-1) location
!------------------------------------------------------------
      jjb1=jjb-1
      t3=tp
      if(ii.ge.2) t3=t13(ii-1,jjb1)
      if(ii.eq.1.and.depmax0.lt..01) t3=1.
      omt3=pai2/t3
      dep1=dep11+dep12
      dep2=dep21+dep22
      slx=(dep11+dep21-dep12-dep22)/(dvarx(ii)+dvarx(i1))
      sly=(dep1-dep2)/(dvary(jjb)+dvary(jjb1))
      depb=(dep1+dep2)/4.0
      if(ick3.eq.1) then
      if(ijb(ii,jjb1).ge.4) depb=d1(ii,jjb1)
      end if
      if(ick4.eq.1) then
      if(d1(ii,jjb1).gt..01.and.d1(i1,jjb1).le.0.)depb=d1(ii,jjb1)
! --* above depb needs to be smoothed in general except for
!     backward reflection location to work
      end if
      um=u1(ii,jjb1)
      vm=v1(ii,jjb1)
      if(depb.lt..01) depb=.01
      call wccg_inline(dmnj(jjb1),depb,um,vm,omt3,cw,cg,sig,akk)
      sigm(jj)=sig
      cmn(jj)=cw
      wlmn(jj)=pai2/akk
      depm(jj)=depb
      if(ijb(ii,jjb1).gt.5) goto 20
      dmean=dmnj(jjb1)
      slf(jj)=slx*cos(dmean)+sly*sin(dmean)
      goto 30
   20 slf(jj)=0.
! -- above mean values are for use of wave breaking dissipation
!
   30 dep11=dep21
      dep12=dep22
   10 continue
! -- end of repetition for j

      return
      end
!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!  subroutine for dispersion relation in current
!      (om-kcos(q)*u-ksin(q)*v)^2=g*k*tanh(k*d)
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      subroutine wccg_inline(q,d,u,v,om,cw,cg,sg,akk)
!-----------------------------------------------------
!  when no solution and cg+u<0, sg=-1,....., are set
!-----------------------------------------------------
      real k,kn
      g=9.806
!     pi=3.1415927
      eps=1.0e-05
      icon=0
!
      an=atan2(v,u+0.000001)
      cqan=cos(q-an)
      ww=sqrt(u**2+v**2)
!
! ** first guess
      x=om*sqrt(d/g)
      y=x*x/(1.-exp(-x**2.5))**.4
      k=y/d
      cusv=cos(q)*u+sin(q)*v
      if(ww.lt..01) go to 42
      if(abs(cqan).lt..5) go to 42
      rkd=min(y,10.)
! ** iteration
      kn=k
      do 30 i=1,40
      c1=k*cqan*ww-om
      f=c1**2-g*k*tanh(rkd)
      fp=2*c1*cusv-g*tanh(rkd)-g*k*d/cosh(rkd)**2
      if(abs(fp).lt.1.e-10) go to 31
      kn=k-f/fp
      if(kn.eq.0.0) go to 40
      if((abs(kn-k)/kn).lt.eps) go to 40
      k=kn
! ** skip of large wavenumber
      vlimit=k*d
      if(abs(vlimit).ge.30.or.vlimit.le..0) go to 31
! ** skip of large wavenumber
   30 continue
   31 icon=1

   40 continue
! ** skip of large wavenumber even convergence
!
      k=kn
   42 continue
      if(k.lt..00001) icon=1
      if(icon.eq.0) then
        rkd=min(k*d,10.)
        sg=sqrt(g*k*tanh(rkd))
        cw=sg/k
        cg=cw*(0.5+k*d/sinh(min(2.*rkd,10.)))
        akk=k
! ** when no solution
          if(cg+cusv.le.0.0) then
          akk=y/d
          sg=om
          cw=sg/akk
          cg=cw*(0.5+y/sinh(min(2.0*y,10.)))
          end if
      else if(icon.eq.1) then
          akk=y/d
          sg=om
          cw=sg/akk
          cg=cw*(0.5+y/sinh(min(2.0*y,10.)))
      end if
!
      return
      end
!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!   setting of matrix components                     
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE SETAB_inline(II,JB,JE,NMX,FC,DFC,IFC,IKAP)
!  akap= coefficient of diffraction term
      USE GLOBAL_INLINE, ONLY: NPF,MPD,IPMX,JPMX,KOMX,NOMX,IGPX,JGPX,MPMX
      COMMON /VPAI/PAI2,PAI,HPAI,RAD,akap,imod,iprp,island,imd,iprpp     &
                   ,nonln,igrav,isolv,ixmdf,iproc,imud,iwnd,depmin0
      COMMON /DATA/NF,MD,IMAX,JMAX,IGMX,JGMX,JCPB,JCPE,JCB,JCE,NFF,MDD
      COMMON /DATB/ICK3,ICK4,ICHOICE,HS0,wd,ws,ph0,wdd(mpd),aslop(ipmx)
      COMMON /DATC/TP,PRD,IBND,VL,TIDE,KPEAK,IBACK,NBLOCK,IWIND,WSMAG
      COMMON /DATD/DX,DY,DXX,DMESH,DTH,kdate,idate,depmax,depmax0
      COMMON /DATG/IBK,DBAR,WL0,WCC(JGPX),HSB(JGPX),HSG(JGPX),DCD(JGPX)
      COMMON /DVAR/DVARX(IGPX),DVARY(JGPX),ETA(IGPX,JGPX),HSK(JGPX)
      COMMON /IJBP/IJB(IGPX,JGPX),DEP(IPMX,JPMX),DEPS(IPMX),DBIG(IPMX)
      COMMON /CCG/CWK(JPMX,2,MPD),CGK(JPMX,2,MPD)
      COMMON /SPECA/SCP(JGPX,MPD),SI(JGPX,NPF,MPD)
      COMMON /AB/AA(5,MPMX),IA(5,MPMX),B(MPMX),X(MPMX)
      COMMON /FDS/FCN(NPF),DCM(MPD),DSFD(NPF,MPD)
      COMMON /DFCN/DF(NPF),FFCN(NPF),DFINP(NPF),PL1E(NPF),PL1T(NPF)
      COMMON /BREK/DEPM(JGPX),DMNJ(JGPX),SLF(JGPX),wlmn(JGPX),cmn(jgpx)  &
                   ,sigm(jgpx),IWVBK,ibk3
      common /rsm/ cosa(mpd),sina(mpd),d1(ipmx,jpmx)
      common /uvp/u(ipmx,jpmx),v(ipmx,jpmx)
      common /uv1/u1(ipmx,jpmx),v1(ipmx,jpmx)
      common /uvwind/u10(ipmx,jpmx),v10(ipmx,jpmx)
      common /fric/bfric(ipmx,jpmx),amud(ipmx,jpmx)
      common /sgma01/sgma0(jgpx,mpd),sgma1(jgpx,mpd)
      COMMON /OUTN/nest,inest(komx),jnest(komx),ix1(igpx),ix2(igpx)
      COMMON /WAVI/H13(IGPX,JGPX),T13(IGPX,JGPX),DMN(IGPX,JGPX)
      COMMON /WAVR/H13R(IGPX,JGPX),T13R(IGPX,JGPX),DMNR(IGPX,JGPX)
      common /comsav/depmin,g,iview,iplane,iwave,cflat,irunup,iroll
!
      g2=g/2.
      om=pai2*fc
      om2=om**2
      om4=4.*om2
      om8=om**8/(pai2*dth)**2/dfc**3/g*2.
      om5=om**5/pai2**2/dfc**4/4./g*2.
      om2g=om2/sqrt(g)/2.
      om35=om**3.5
      wss=ws*tanh(float(ii)*50./dvarx(1))
!     ifast=0
!     if(iprpp.ne.-1.and.dxx.gt.300.) ifast=1
      winp1=0.00002*om2/g/dth/.03
      deca1=0.01*om35/g/dth/dfc
      if(md.eq.7.and.iprp.eq.2) then
        winp1=0.000027*om2/g/dth/.03
        deca1=0.012*om35/g/dth/dfc
      end if
!
!     if(iwind.ge.1) then
!       winp1=0.00002*om2/g/dth/.03
!       deca1=0.01*om35/g/dth/dfc
!     end if
      if(dvarx(ii).gt.300.) then
        winp1=winp1*sqrt(300./dvarx(ii))
        deca1=deca1*(300./dvarx(ii))**.6
      end if
!
      if(igrav.eq.3) then
        winp1=winp1*2.
      end if
!
!     winp1=0.000021*om2/g/dth/.03
      winp2=winp1
!     winp1=0.00003*om2/g/dth/.03
!     deca1=0.023*om35/g/dth/dfc
      aakap=akap
      if(ii.le.ikap.and.akap.ge.1.) aakap=1.
      afact=1.
      if(ichoice.eq.2) afact=(float(mdd+1)/float(md+1))**1.4
      if(ix1(ii).eq.0.and.aakap.gt.1.) aakap=1.
      winln=10.
      if(ws.ge..1) then
!       winln=5.
        if(hs0.lt..5) winln=1.5
        ph00=ph0*3.
        if(ph00.gt.1.25) ph00=1.25
! --- the above ph00 keep long waves in the upwind lake area.
        om51=om**5*exp(0.74*min((ph00/om)**4,10.))/dth/.03
      end if
!
      gamma1=0.0001*aakap
      gamma2=0.0001*aakap
!
      JMESH=JE-JB+1
      NMX=(MD-imd+1)*JMESH
      NN=0
      DO 10 N=1,NMX
        DO 20 M=1,5
          AA(M,N)=0.
          IA(M,N)=0
   20   CONTINUE
        B(N)=0.
   10 CONTINUE
      jgrav=0
      if(igrav.ge.1) then
        do jm=jb,je
          if(ijb(ii,jm).eq.4.or.ijb(ii,jm).eq.5) then
            jgrav=jgrav+1
          end if
        end do
      end if

! -- iteration for J mesh
      DO 30 JM=1,JMESH
        JBM=JB+JM
        JJ=JBM-1
        t3=tp
        if(ii.ge.2) t3=t13(ii-1,jj)
        if(ii.eq.1.and.depmax0.lt..01) t3=1.
        om1=pai2/t3/om
!
        if(iwind.ge.1) then
          ws=sqrt(v10(ii,jj)**2+u10(ii,jj)**2)
          wd=atan2(v10(ii,jj),u10(ii,jj))
          if(ws.lt.wsmag) ws=wsmag
        end if
!
        if(ws.gt..1) then
          om1=(pai2/(t3+1.))/om
          if(nonln.eq.1) winp1=winp2*1.5
        end if
!
        nln=0
        if(nonln.eq.1) then
          cc=om1
          if(om1.gt.winln) cc=winln
          cc1=(cosh(min(9.424,om**2*abs(d1(ii,jj))/g))/6200.)
          cc1=min(tanh((d1(ii,jj)/100.)**2),cc1)
          cc2=cc**12*cc1**2
          om52=om5*cc2
          om81=om8*cc2
          if(om1.lt.winln) nln=1
        end if
!
        hsij=hs0
        if(ii.ge.2) hsij=h13(ii-1,jj)
        cc=hsij/(abs(d1(ii,jj))+.01)
        if(cc.gt.1.) cc=1.
        if(imud.ge.0) then
          aamud=amud(ii,jj)*cc
        else
          aamud=sqrt(2.*amud(ii,jj)/om)*cc
        end if
!
!  following gam is changed in coefficiets A1, A2, A3      
        gam=2./om/dvary(jj)**2*gamma1
!  above gam should be changed due to spatial variation of intrin-
!  sic aungular frequency due to spatial change of current field
!
        JBV=IJB(II,JJ)
        J1=JM+1
        II1=II+1
        hx1=(dep(ii,jj)+dep(ii,jbm))/2.
        hx2=(dep(ii1,jj)+dep(ii1,jbm))/2.
        hy1=(dep(ii,jj)+dep(ii1,jj))/2.
        hy2=(dep(ii,jbm)+dep(ii1,jbm))/2.
        ht=d1(ii,jj)
! --  do not use the smooth for above ht, can cause asymmetric
        if(ht.lt..01) ht=.01
        ux1=(u(ii,jj)+u(ii,jbm))/2.
        ux2=(u(ii1,jj)+u(ii1,jbm))/2.
        uy1=(u(ii,jj)+u(ii1,jj))/2.
        uy2=(u(ii,jbm)+u(ii1,jbm))/2.
        vx1=(v(ii,jj)+v(ii,jbm))/2.
        vx2=(v(ii1,jj)+v(ii1,jbm))/2.
        vy1=(v(ii,jj)+v(ii1,jj))/2.
        vy2=(v(ii,jbm)+v(ii1,jbm))/2.
        um=u1(ii,jj)
        vm=v1(ii,jj)
        dhx=hx2-hx1
        dhy=hy2-hy1
        dux=ux2-ux1
        dvx=vx2-vx1
        duy=uy2-uy1
        dvy=vy2-vy1
        JJM=JM
!
! -- setup of energy dissipation term
! -- JJM= local #, JJ= GLOBAL_INLINE #
        cab=0.
        if (IWVBK.eq.1) then
          CALL WVBRK1_inline(CAB,JBV,II,JJM,JJ,fc)
        elseif (IWVBK.eq.2) then
          CALL WVBRK2_inline(CAB,JBV,II,JJM,JJ,om)
        elseif (IWVBK.eq.3) then
          CALL WVBRK3_inline(CAB,JBV,II,JJM,JJ)
        elseif (IWVBK.eq.4) then
          CALL WVBRK4_inline(CAB,JBV,II,JJM,JJ)
	    elseif (IWVBK.eq.0) then
          CALL WVBRK_inline(CAB,JBV,II,JJM,JJ,um,vm,fc)
        end if
!
        thc=float(imd-1)*dth-hpai
        if(iview.ge.1) thc=thc-dth/2.
! -- iteration for k direction
        DO 60 KK=imd,MD
          thc=thc+dth
!
          CWX1=(CWK(J1,1,kk)+CWK(JM,1,kk))/2.0
          CWX2=(CWK(J1,2,kk)+CWK(JM,2,kk))/2.0
          CWY1=(CWK(JM,1,kk)+CWK(JM,2,kk))/2.0
          CWY2=(CWK(J1,1,kk)+CWK(J1,2,kk))/2.0
          CGX1=(CGK(J1,1,kk)+CGK(JM,1,kk))/2.0
          CGX2=(CGK(J1,2,kk)+CGK(JM,2,kk))/2.0
          CGY1=(CGK(JM,1,kk)+CGK(JM,2,kk))/2.0
          CGY2=(CGK(J1,1,kk)+CGK(J1,2,kk))/2.0
          ccg1=cwy1*cgy1
          ccg2=cwy2*cgy2
          ccgm=0.5*(ccg1+ccg2)
          ccgn=0.5*(cwx1*cgx1+cwx2+cgx2)
          cgu1=(cgx1*cosa(kk)+ux1)/dxx
          cgu2=(cgx2*cosa(kk)+ux2)/dxx
          cgv1=(cgy1*sina(kk)+vy1)/dvary(jj)
          cgv2=(cgy2*sina(kk)+vy2)/dvary(jj)
          ccgt=(ccgm/dxx**2+ccgn/dvary(jj)**2)/2./om*gamma1
          if(cgu1.le.0.0) then
            uuu1=0.0
          else
            uuu1=cgu1
          end if
          if(cgu2.le.0.0) then
            uuu2=0.0
          else
            uuu2=cgu2
          end if
          vvv1=cgv1
          vvv2=cgv2
!
          q1=dcm(kk)-dth/2.0
          call wccg_inline(q1,ht,um,vm,om,cw1,cg1,sg1,akk1)
          rkd1=akk1*ht
          sgh1=g2*akk1**2/cosh(min(10.,rkd1))**2/sg1
          q2=dcm(kk)+dth/2.0
          call wccg_inline(q2,ht,um,vm,om,cw2,cg2,sg2,akk2)
          rkd2=akk2*ht
          sgh2=g2*akk2**2/cosh(min(10.,rkd2))**2/sg2
          rkd=(rkd1+rkd2)/2.
          akk=(akk1+akk2)/2.
          cgg=(cg1+cg2)/2.
          agg=akk*cgg/dth**2*gamma2
          ann=(cg1+cg2)/(cw1+cw2)
          aan=(1.+(2.*ann-1.)**2*cosh(min(10.,2.*rkd)))/2./ann**2-1.
          bbn=aan/ann/om
!
          cc=1./dth
          if(igrav.ge.1) then
            if(igrav.eq.2.or.jgrav.ge.2) then
              cc=tanh(5.*fc*hsij**2)*tanh(sqrt(g*hsij)/(ws+.01))*cc
            end if
          end if
          sin2q1=sin(2.*q1)/2.
          sin2q2=sin(2.*q2)/2.
          vt1=sin2q1*dux/dxx+sin(q1)**2*dvx/dxx-(cos(q1)**2*duy+sin2q1*dvy) &
              /dvary(jj)+(sgh1*dhx/dxx*sin(q1)-sgh1*dhy/dvary(jj)*cos(q1))/akk1
          vt1=vt1*cc
          vt2=sin2q2*dux/dxx+sin(q2)**2*dvx/dxx-(cos(q2)**2*duy+sin2q2*dvy) &
              /dvary(jj)+(sgh2*dhx/dxx*sin(q2)-sgh2*dhy/dvary(jj)*cos(q2))/akk2
          vt2=vt2*cc
!
          NN=NN+1
!
! -- coef.of A1
          ccg=gam*ccgm
!              *(2.0*pai2*fc)/abs(sgma1(jb+jm-1,kk))
! -- absolute w was good, so above statement was omitted.
!
          winp=0.
          sinp=0.
          fd1=1.
!
!         wc=ws*1.05
          fd1=(1.+ws*.3)/sqrt(cgy1**2+ws**2+2.*cgy1*ws*cos(wd-thc)+.1)
!         fd1=(1.+ws*.8)/sqrt(cgy1**2+ws**2+2.*cgy1*ws*cos(wd-thc)+.1)
          if(fd1.gt.2.) fd1=2.
!
          ac=0.
          if(scp(jj,kk).gt..0) ac=2.*sqrt(scp(jj,kk)*om)
!         fd2=1./tanh(om*(abs(d1(ii,jj))+.01)/cwy1)
          fd2=1./tanh(rkd)
          deca=deca1*(ac/cwy1)**1.5*cgy1*fd1*fd2
          if(iback.eq.0) deca=deca/(1.+5.*tanh(hsg(jj)))
!
          f1=1.
          if(ws.ge..1.and.iback.eq.0) then
            cc2=cos((wd-thc)/2.)
!           f1=wss*cc2-cgy1
            f1=wss*cos(wd-thc)-cgy1
!           if(d1(ii,jj).lt.3.) f1=f1*tanh(d1(ii,jj)/3.)
            if(d1(ii,jj).lt.5.) f1=f1*(2.-tanh(d1(ii,jj)/2.))
!           if(ii.ge.2) then
!             if(d1(ii-1,jj).gt..1) then
!               f1=f1*(1.+cos(dmn(ii-1,jj)*rad-wd)**2)/2.
!             end if
!           end if
!
            if(f1.lt.0.) f1=0.
            if(nonln.eq.1.and.f1.lt.1.) f1=1.
            f2=1.
            if(cgy1.lt.ws) f2=ws/cgy1
            f3=log10(f2)
            winp=winp1*f1*f3/f2
            cc1=1.-tanh(d1(ii,jj)/5.)
            cc3=1.
            if(ii.ge.2) then
              cc3=cos((dbig(jj)-wd)/2.)
!             if(d1(ii,jj).eq.d1(ii-1,jj)) then
              winp=winp/(.8+3.*tanh(hsg(jj))+.5*cc1)*(cc2*cc3)**10
!             else
!             winp=winp/(.8+3.*tanh(hsg(jj))+.5*cc1)*(cc2*cc3)**50
!             end if
            end if
!           write(69,'(2i5,5f13.6)') ii,jj,fc,winp,f1,f2,f3
!
!           if(iwind.ge.1) then
!             ph00=33./ws
!             if(ph00.gt.1.25) ph00=1.25
!             om51=om**5*exp(0.74*min((ph00/om)**4,10.))/dth/.03
!             wdd4=0.005*cos(thc-wd)**40
!             0.005=0.0002*g*8/3.14159
!             sinp=wdd4*f1/f2/om51*afact
!           else
            sinp=wdd(kk)*f1/f2/om51*afact
!           end if
!
            if(ii.ge.2) then
              if (d1(ii-1,jj).lt..1) sinp=0.
            endif
!           sinp=sinp*(cc2*cc3)**4
            if(scp(jj,kk).lt.sinp) scp(jj,kk)=sinp
!           if(ifast.eq.1) then
!             if(deca.gt.winp) deca=winp
!           end if
          end if
!
          amp=.01
          if(ii.ge.2) amp=h13(ii-1,jj)/sqrt(ht)
          cc=bfric(ii,jj)
          if(cc.gt..08.and.amp.gt..5) then
            if(tp.ge.8.) then
              cc=.08+(cc-.08)/4.
            else
              cc=.08+(cc-.08)/2.
            end if
          end if
          fric=om2g*cc*amp/sinh(min(5.,rkd))**2
!
          winp=winp*sqrt(afact)
          if(iprp.eq.3.and.jbv.eq.2) then
            winp=0.
            deca=0.
          end if
!         if(hs0.lt..3.and.ws.lt..01) deca=0.
          if(ws.lt..01.and.iback.eq.0) deca=0.
!
          if(imud.ge.0) then
            a1=uuu2+max(cab,deca)+ccg-winp+fric+om4/cwy1**2*aamud
!     amud is the viscosity that has the unit of m*m/s
!     (for the sea water at 35 F, amud=0.0000018 m*m/s)
          else
            cc=rkd1+rkd2+sinh(min(5.,rkd1+rkd2))
            a1=uuu2+max(cab,deca)+ccg-winp+fric+aamud*cgg*akk**2/cc
          end if
!
          if(nln.eq.1) then
            cc=scp(jj,kk)+1.e-10
            if(rkd.lt.pai) then
              c1=scp(jj,kk)
              if(kk.ne.imd) c1=scp(jj,kk-1)
              c2=scp(jj,kk)
              if(kk.ne.md) c2=scp(jj,kk+1)
              a1=a1-bbn*(c1**3+c2**3-2.*scp(jj,kk)**3)*akk**3*om81*ann**4/cc
            end if
! --- B= k^3 sig^8 n^4 [(fp/f)^4 N]^3/g = k^3 sig^5 n^4 [(fp/f)^4 F]^3/(g pai2^2)
            c1=0.
            c2=0.
            if(ifc.ne.1) c1=si(jj,ifc-1,kk)
            if(ifc.ne.nf) c2=si(jj,ifc+1,kk)
            a1=a1-ann*(3.*c2**3-2.*si(jj,ifc,kk)**3-c1**3)*akk**3*om52*ann**4/cc
          end if
!
          if(vvv2.ge.0.0) a1=a1+vvv2
          if(vvv1.lt.0.0) a1=a1-vvv1
          if(vt1.ge.0.0.and.vt2.ge.0.0) a1=a1+vt2
          if(vt1.lt.0.0.and.vt2.ge.0.0) a1=a1+vt2-vt1
          if(vt1.lt.0.0.and.vt2.lt.0.0) a1=a1-vt1
          if(jj.ne.jb.and.jj.ne.je) a1=a1+ccgt
          a1=a1+2.*agg
          aa(1,nn)=a1
!         if(sgma1(jj,kk).lt.0.0) aa(1,nn)=1.0
          ia(1,nn)=nn

! -- coef.of A2
          a2=gam*(-ccg1+0.5*ccgm)
!            *(2.0*pai2*fc)/abs(sgma1(jb+jm-1,kk))
! -- absolute w was good, so above statement was omitted.
          if(vvv1.ge.0.0) a2=a2-vvv1
          aa(2,nn)=a2
!         if(sgma1(jj,kk).lt.0.0) aa(2,nn)=0.0
          ia(2,nn)=nn-(md-imd+1)

! -- coef.of A3
          a3=gam*(-ccg2+0.5*ccgm)
!            *(2.0*pai2*fc)/abs(sgma1(jb+jm-1,kk))
! -- absolute w was good, so above statement was omitted.
          if(vvv2.lt.0.0) a3=a3+vvv2
          aa(3,nn)=a3
!         if(sgma1(jj,kk).lt.0.0) aa(3,nn)=0.0
          ia(3,nn)=nn+md-imd+1

! -- coef.of A4
          a4=0.0
          if(vt1.ge.0.0) a4=-vt1
          a4=a4-agg
          if(mod(nn-1,md-imd+1).eq.0) then
            aa(4,nn)=0.0
            ia(4,nn)=0
          else
            aa(4,nn)=a4
!           if(sgma1(jj,kk).lt.0.0) aa(4,nn)=0.0
            ia(4,nn)=nn-1
          end if

! -- coef.of A5
          a5=0.0
          if(vt2.lt.0.0) a5=vt2
          a5=a5-agg
          if(mod(nn,md-imd+1).eq.0) then
            aa(5,nn)=0.0
            ia(5,nn)=0
          else
            aa(5,nn)=a5
!           if(sgma1(jj,kk).lt.0.0) aa(5,nn)=0.0
            ia(5,nn)=nn+1
          end if

! -- coef.of B(N)
          if(jbv.eq.6.or.jbv.eq.7.or.jbv.eq.8.or.jbv.eq.9) then
            b(nn)=0.0
          else
!           beta=1.
!           scpc=0.
!           if(jj.ne.jb.and.jj.ne.je) then
!             do ikk=imd,md
!               scpc=scpc+scp(jj,ikk)/cosh(abs(kk-ikk)*dth*beta)**2
!             end do
!             scpc=scpc*dth*beta/2.
!           end if
            b(nn)=uuu1*scp(jj,kk)
            if(jj.ne.jb.and.jj.ne.je) b(nn)=b(nn)+ccgt*scp(jj,kk)
!           if(jj.ne.jb.and.jj.ne.je) b(nn)=b(nn)+akk*cgg*scpc
!           if(sgma1(jj,kk).lt.0.0) b(nn)=0.0
          end if
!
! -- boundary conditions at j=1
          if(jm.eq.1.and.jbv.eq.2) then
            a1=a1+a2
            aa(1,nn)=a1
!           if(sgma1(jj,kk).lt.0.0) aa(1,nn)=1.0
            aa(2,nn)=0.0
            ia(2,nn)=0
          else if(jm.eq.1.and.jbv.eq.3) then
            aa(2,nn)=0.0
            ia(2,nn)=0 
          else if(jm.eq.1.and.jbv.eq.9) then
            aa(2,nn)=0.0
            ia(2,nn)=0
            b(nn)=0.0
          else if(jm.eq.1.and.jbv.eq.4) then
            vkr=0.0
            agl=0.0
!** ckr_inline(ii,jj),akr_inline(ii,jj) are functions **
            if(ick3.eq.1) then
              vkr=ckr_inline(ii,jj)
              agl=akr_inline(ii,jj)
            end if
!** ckr_inline(ii,jj),akr_inline(ii,jj) are functions **
            a2=vkr**2*a2
!************************************
            if(thc.lt.agl*rad) a2=0.0
!************************************
            na2=nn+md-imd+2-2*kk + int(2.*agl*rad/dth)
            if(na2.eq.ia(1,nn)) then
              aa(1,nn)=a1+a2
!             if(sgma1(jj,kk).lt.0.0) aa(1,nn)=1.0
            else if(na2.eq.ia(3,nn)) then
              aa(3,nn)=a3+a2
!             if(sgma1(jj,kk).lt.0.0) aa(3,nn)=0.0
            else if(na2.eq.ia(4,nn)) then
              aa(4,nn)=a4+a2
!             if(sgma1(jj,kk).lt.0.0) aa(4,nn)=0.0
            else if(na2.eq.ia(5,nn)) then
              aa(5,nn)=a5+a2
!             if(sgma1(jj,kk).lt.0.0) aa(5,nn)=0.0
            else
              aa(2,nn)=a2
!             if(sgma1(jj,kk).lt.0.0) aa(2,nn)=0.0
              ia(2,nn)=na2
            end if
          else if(jm.eq.1.and.jbv.eq.6) then
            vkr=0.0
            agl=0.0
!** ckr_inline(ii,jj),akr_inline(ii,jj) are functions **
            if(ick3.eq.1) then
              vkr=ckr_inline(ii,jj)
              agl=akr_inline(ii,jj)
            end if
!** ckr_inline(ii,jj),akr_inline(ii,jj) are functions **
            a2=vkr**2*a2
!***********************************
            if(thc.lt.agl*rad) a2=0.0
!***********************************
            na2=nn+md-imd+2-2*kk + int(2.*agl*rad/dth)
            b(nn)=0.0
            if(na2.eq.ia(1,nn)) then
              aa(1,nn)=a1+a2
!             if(sgma1(jj,kk).lt.0.0) aa(1,nn)=1.0
            else if(na2.eq.ia(3,nn)) then
              aa(3,nn)=a3+a2
!             if(sgma1(jj,kk).lt.0.0) aa(3,nn)=0.0
            else if(na2.eq.ia(4,nn)) then
              aa(4,nn)=a4+a2
!             if(sgma1(jj,kk).lt.0.0) aa(4,nn)=0.0
            else if(na2.eq.ia(5,nn)) then
              aa(5,nn)=a5+a2
!             if(sgma1(jj,kk).lt.0.0) aa(5,nn)=0.0
            else
              aa(2,nn)=a2
!             if(sgma1(jj,kk).lt.0.0) aa(2,nn)=0.0
              ia(2,nn)=na2
            end if
          end if
!
! -- boundary conditions at j=end
          if(jm.eq.jmesh.and.jbv.eq.2) then
            a1=a1+a3
            aa(1,nn)=a1
!           if(sgma1(jj,kk).lt.0.0) aa(1,nn)=1.0
            aa(3,nn)=0.0
            ia(3,nn)=0
          else if(jm.eq.jmesh.and.jbv.eq.3) then
            aa(3,nn)=0.0
            ia(3,nn)=0
          else if(jm.eq.jmesh.and.jbv.eq.9) then
            aa(3,nn)=0.0
            ia(3,nn)=0
            b(nn)=0.0
          else if(jm.eq.jmesh.and.jbv.eq.5) then
            vkr=0.0
            agl=0.0
!** ckr_inline(ii,jj),akr_inline(ii,jj) are functions **
            if(ick3.eq.1) then
              vkr=ckr_inline(ii,jj)
              agl=akr_inline(ii,jj)
            end if
!** ckr_inline(ii,jj),akr_inline(ii,jj) are functions **
            a3=vkr**2*a3
!**********************************
            if(thc.gt.agl*rad) a3=0.0
!**********************************
            na3=nn+md-imd+2-2*kk + int(2.*agl*rad/dth)
            if(na3.eq.ia(1,nn)) then
              aa(1,nn)=a1+a3
!             if(sgma1(jj,kk).lt.0.0) aa(1,nn)=1.0
            else if(na3.eq.ia(2,nn)) then
              aa(2,nn)=a2+a3
!             if(sgma1(jj,kk).lt.0.0) aa(2,nn)=0.0
            else if(na3.eq.ia(4,nn)) then
              aa(4,nn)=a4+a3
!             if(sgma1(jj,kk).lt.0.0) aa(4,nn)=0.0
            else if(na3.eq.ia(5,nn)) then
              aa(5,nn)=a5+a3
!             if(sgma1(jj,kk).lt.0.0) aa(5,nn)=0.0
            else
              aa(3,nn)=a3
!             if(sgma1(jj,kk).lt.0.0) aa(3,nn)=0.0
              ia(3,nn)=na3
            end if
          else if(jm.eq.jmesh.and.jbv.eq.7) then
            vkr=0.0
            agl=0.0
!** ckr_inline(ii,jj),akr_inline(ii,jj) are functions **
            if(ick3.eq.1) then
              vkr=ckr_inline(ii,jj)
              agl=akr_inline(ii,jj)
            end if
!** ckr_inline(ii,jj),akr_inline(ii,jj) are functions **
            a3=vkr**2*a3
! --********************************
            if(thc.gt.agl*rad) a3=0.0
! --********************************
            na3=nn+md-imd+2-2*kk + int(2.*agl*rad/dth)
            b(nn)=0.0
            if(na3.eq.ia(1,nn)) then
              aa(1,nn)=a1+a3
!             if(sgma1(jj,kk).lt.0.0) aa(1,nn)=1.0
            else if(na3.eq.ia(2,nn)) then
 
              aa(2,nn)=a2+a3
!             if(sgma1(jj,kk).lt.0.0) aa(2,nn)=0.0
            else if(na3.eq.ia(4,nn)) then
              aa(4,nn)=a4+a3
!             if(sgma1(jj,kk).lt.0.0) aa(4,nn)=0.0
            else if(na3.eq.ia(5,nn)) then
              aa(5,nn)=a5+a3
!             if(sgma1(jj,kk).lt.0.0) aa(5,nn)=0.0
            else
              aa(3,nn)=a3
!             if(sgma1(jj,kk).lt.0.0) aa(3,nn)=0.0
              ia(3,nn)=na3
            end if
          end if
!
! -- iteration for k
   60   CONTINUE
! -- iteration for j
   30 CONTINUE

      RETURN
      END SUBROUTINE SETAB_inline

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!  subroutine to obtain wave breaking energy dissipation term
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE WVBRK2_inline(CAB,JBV,II,JM,JJ,om)
!---------------------------------------------------------------
!  energy dissipation term is Battjes and Janssen's(1978)
!    bore model with a breaker parameter based on Miche's criterion of
!    gama=0.73
!  parameters used in this subroutine are calculated 
!    in subroutine veloc
!---------------------------------------------------------------
      USE GLOBAL_inline, ONLY: NPF,MPD,IPMX,JPMX,IGPX,JGPX
      USE global_inline, ONLY: gamma_bj78                    !added MEB 10/19/2021
      COMMON /DATD/DX,DY,DXX,DMESH,DTH,kdate,idate,depmax,depmax0
      COMMON /DATG/IBK,DBAR,WL0,WCC(JGPX),HSB(JGPX),HSG(JGPX),DCD(JGPX)
      COMMON /BREK/DEPM(JGPX),DMNJ(JGPX),SLF(JGPX),wlmn(JGPX),cmn(jgpx)  &
                   ,sigm(jgpx),IWVBK,ibk3
      COMMON /VPAI/PAI2,PAI,HPAI,RAD,akap,imod,iprp,island,imd,iprpp     &
                   ,nonln,igrav,isolv,ixmdf,iproc,imud,iwnd,depmin0

! --- ibk=1 and jbv>=6 means no wave breaking
      IF(IBK.EQ.1) GOTO 10
      IF(JBV.GE.6) GOTO 10

	  slj=slf(jm)
      if(slj.ge.0.04) go to 10
      
      if(gamma_bj78 .ne. -1) then
        gama = gamma_bj78          !allow use of user-specified value if requested  MEB 10/19/2021
      else
        gama=0.73                  !otherwise use default value.
      endif
    
!     alfabj=1.0/1.414
	  alfabj=0.707
!
      H13=HSB(JJ)
!     hrms=h13/sqrt(2.0)
	  hrms=h13/1.414
	  dep=depm(jm)
	  wnk=pai2/wlmn(jm)
	  sig=sigm(jm)
	
      hmax=gama*dep 

	  if(wnk.lt.0.0) goto 10

	  if(hmax.gt.0.0.and.hrms.gt.0.0) then
	    B=hrms/hmax
	  else
	    B=0.0
	  endif
	  if(B.le.0.5) then
	    Q0=0.0
	  elseif(B.lt.1.0) then
	    Q0=(2.0*B-1.0)**2
	  endif
	  if(B.le.0.2) then
	    Qb=0.0
      elseif(B.lt.1.0) then
	    B2=B*B
	    EQ0=exp((Q0-1.0)/B2)
	    Qb=Q0-B2*(Q0-EQ0)/(B2-EQ0)
	  else
	    Qb=1.0
	  endif

!     Dab=0.25*alfabj*Qb*(sig/pai2)*hmax**2
	  Dab=0.04*alfabj*Qb*sig*hmax**2
!     cab=Dab/(sig*(1.0/8.0)*hrms**2)
	  cab=8.0*dab/sig/hrms**2*om
	  goto 20
 
   10 CAB=0.0

   20 continue
      RETURN
      END SUBROUTINE WVBRK2_inline

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!  subroutine to obtain wave breaking energy dissipation term
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE WVBRK4_inline(CAB,JBV,II,JM,JJ)
!---------------------------------------------------------------
!  energy dissipation term is Battjes and Janssen's(2007)
!    also based on Alsina and Baldock (2007)
!---------------------------------------------------------------
      USE GLOBAL_inline, ONLY: NPF,MPD,IPMX,JPMX,IGPX,JGPX
      COMMON /DATD/DX,DY,DXX,DMESH,DTH,kdate,idate,depmax,depmax0
      common /DATG/IBK,DBAR,WL0,WCC(JGPX),HSB(JGPX),HSG(JGPX),DCD(JGPX)
      COMMON /BREK/DEPM(JGPX),DMNJ(JGPX),SLF(JGPX),wlmn(JGPX),cmn(jgpx)  &
                   ,sigm(jgpx),IWVBK,ibk3
      COMMON /VPAI/PAI2,PAI,HPAI,RAD,akap,imod,iprp,island,imd,iprpp     &
                   ,nonln,igrav,isolv,ixmdf,iproc,imud,iwnd,depmin0
      common /WAVS/H13S(IGPX,JGPX),IBR(IGPX,JGPX),DISS(IGPX,JGPX)
    
      cab=0.0
    
! --- ibk=1 and jbv>=6 means no wave breaking
      if(IBK.EQ.1.or.JBV.GE.6) return
		
      wnk=pai2/wlmn(jm)
	
      i=max(ii-1,1)
      j=max(jj-1,1)
      j2=min(jj+1,JGPX)
      ijbr=ibr(i,jj)+ibr(i,j)+ibr(i,j2)

      y2=0.76*wnk*depm(jm)+0.29  !Grasmeijer
      y=max(0.64,y2) 
	
      Hmax=0.88/wnk*tanh(y*wnk*d/0.88)       
      Hmax=max(Hmax,0.64*depm(jm))
      if(Hsb(jj).le.0.72*Hmax.and.ijbr.eq.0) return  !Alex, for random waves

      fp=sigm(jm)/pai2 !Frequency, Hz
      Hrms=Hsb(jj)/1.414

      R=Hmax/Hrms
      Qb=1.0+0.7523*(R**3+1.5*R)/exp(R**2)-erf_inline(R)
!     0.7523=4/3/sqrt(pi)
      Db=-2.60127*fp/depm(jm)*Qb*Hrms**3
!     2.60127=3/16*sqrt(pi)*9.8, Note density not included
      cab=-0.8154944*Db/Hrms**2
!     Dissipation coefficient, 0.8154944=8.0/9.8

      RETURN
      END SUBROUTINE WVBRK4_inline

!*****************************************     
      real function erf_inline(x)
!Calculates the error function
!used in wave breaking formulation of Alsina and Baldock (2007)
!*****************************************
      tol=0.0001
      fac = 1.1284*x
!     1.1284=2/sqrt(pi)
      eps = tol/fac
!
        if (x.gt.3.0) then
        erf_inline = 1.0
        return	
        endif
      E=1.0
      S=1.0
      do k=1,40
      dk = dble(k)
      E = -((2.0*dk-1.0)*x**2.0)/((2.0*dk+1.0)*k)*E
      if (abs(E).le.eps) then
      exit
      endif
      S = S + E
      enddo
      erf_inline = fac*S
      if (erf_inline.gt.1.0) erf_inline = 1.0
      return
      end function erf_inline

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!  subroutine to obtain wave breaking energy dissipation term
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE WVBRK_inline(CAB,JBV,II,JM,JJ,um,vm,fc)
!---------------------------------------------------------------
!  energy dissipation term is the extended (by Sakai et al.)
!    Goda's breaker index in order to include current effects
!    the Rayleigh distribution assumption for wave height is used
!  parameters used in this subroutine are calculated 
!    in subroutine veloc
!---------------------------------------------------------------
      USE GLOBAL_inline, ONLY: A,NPF,MPD,IPMX,JPMX,IGPX,JGPX
      COMMON /DATC/TP,PRD,IBND,VL,TIDE,KPEAK,IBACK,NBLOCK,IWIND,WSMAG
      COMMON /DATD/DX,DY,DXX,DMESH,DTH,kdate,idate,depmax,depmax0
      COMMON /DATG/IBK,DBAR,WL0,WCC(JGPX),HSB(JGPX),HSG(JGPX),DCD(JGPX)
      COMMON /BREK/DEPM(JGPX),DMNJ(JGPX),SLF(JGPX),wlmn(JGPX),cmn(jgpx)  &
                   ,sigm(jgpx),IWVBK,ibk3
      COMMON /VPAI/PAI2,PAI,HPAI,RAD,akap,imod,iprp,island,imd,iprpp     &
                   ,nonln,igrav,isolv,ixmdf,iproc,imud,iwnd,depmin0

! -- ibk=1 and jbv>=6 means no wave breaking
!     CAM=1.0
      ced=1.0
      slj=slf(jm)
!
      if(slj.ge.0.04) go to 10
      if(slj.le.0.00) go to 10
!
      ALF=1.6
      g=9.806
!     dmean=dmnj(jj)
!     ujm=um*cos(dmean)+vm*sin(dmean)
!     if(ujm.lt.0.0) then
!     qjm=(-ujm)*depm(jm)
!     qstar=qjm/(g*g*TP**3)
!     ed=qstar*slj**(0.25)/(depm(jm)/wL0)
      if(um.lt.0.) then
        qstar=sqrt(um**2+vm**2)/g**2/TP**3
        ed=qstar*slj**.25*WL0
	    if    (ed.lt.0.0005) then
          ced=1.0
        elseif(ed.gt.0.0024) then
          ced=0.506
	    else
          ced=1.13-260.0*ed
        endif
      else
        ced=1.0
      endif
!
      IF(IBK.EQ.1) GOTO 10
      IF(JBV.GE.6) GOTO 10
!
!     SL43=SLJ**(4.0/3.0)
      SL43=SLJ**1.333333
      TANB=1.0+15.0*SL43
!     EXD=EXP(-1.5*PAI*DEPM(JM)*TANB/WL0)
      EXD=EXP(-4.712389*DEPM(JM)*TANB/WL0)
      HB=A*WL0*(1.0-EXD)*ced               
!     DLHB=A*1.5*PAI*DXX*TANB*SLJ*EXD*ced/2.0
      DLHB=A*0.75*PAI*DXX*TANB*SLJ*EXD*ced
      H13=HSB(JJ)
      X1=(ALF*(HB+DLHB)/H13)**2
      X2=(ALF*(HB-DLHB)/H13)**2
      IF(X1.GT.25.) GOTO 30
!     PX1=PAI*X1/4.0
      PX1=.7853982*X1
      PR1=(1.0+PX1)/EXP(PX1)
      GOTO 40
   30 PR1=0.0
   40 IF(X2.GT.25.) GOTO 50
!     PX2=PAI*X2/4.0
      PX2=.7853982*X2
      PR2=(1.0+PX2)/EXP(PX2)
      GOTO 60
   50 PR2=0.0
!  60 CAB=CAM*(PR2-PR1)/(1.0-PR1)*CMN(JM)/DXX
   60 CAB=(PR2-PR1)/(1.0-PR1)*CMN(JM)/DXX
!  60 CAB=(PR2-PR1)/(1.-PR1)*fc
      GOTO 20

   10 CAB=0.0

   20 continue
      RETURN
      END SUBROUTINE WVBRK_inline

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!  subroutine to obtain wave breaking energy dissipation term
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE WVBRK3_inline(CAB,JBV,II,JM,JJ)
!---------------------------------------------------------------
!  energy dissipation term is Chawla and Kirby (2002)
!    breaker model in order to include both current and depth limited
!    wave breaking
!    following Thornton and Guza's approach (1983)
!  parameters used in this subroutine are calculated in subroutine velo!  
!---------------------------------------------------------------
      USE GLOBAL_inline, ONLY: NPF,MPD,IPMX,JPMX,IGPX,JGPX
      COMMON /DATD/DX,DY,DXX,DMESH,DTH,kdate,idate,depmax,depmax0
      COMMON /DATG/IBK,DBAR,WL0,WCC(JGPX),HSB(JGPX),HSG(JGPX),DCD(JGPX)
      COMMON /BREK/DEPM(JGPX),DMNJ(JGPX),SLF(JGPX),wlmn(JGPX),cmn(jgpx)  &
                   ,sigm(jgpx),IWVBK,ibk3
      COMMON /VPAI/PAI2,PAI,HPAI,RAD,akap,imod,iprp,island,imd,iprpp     &
                   ,nonln,igrav,isolv,ixmdf,iproc,imud,iwnd,depmin0

! -- ibk=1 and jbv>=6 means no wave breaking
      IF(IBK.EQ.1) GOTO 10
      IF(JBV.GE.6) GOTO 10

	  slj=slf(jm)
      if(slj.ge.0.04) go to 10
      g=9.806
      gama=0.6
      beta=0.4
      H13=HSB(JJ)
!     hrms=h13/sqrt(2.0)
	  hrms=h13/1.414
	  dep=depm(jm)
	  wnk=pai2/wlmn(jm)
	  sig=sigm(jm)

	  if(wnk.lt.0.0) goto 10

	  tkh=wnk/tanh(wnk*dep)
	  ab=1.0+(tkh*hrms/gama)**2
!     dab1=1.0-1.0/ab**(5.0/2.0)
	  dab1=1.0-1.0/ab**2.5
	  dab2=(tkh/gama)**2
	  dab3=sqrt(g*tkh)
!     dab3=sqrt(g*tkh)*g*wnk
! --- is it g*wnk missing in the equation?
!     dab=3.0*beta*wnk*dab3*dab2*dab1*hrms**5/(32.0*sqrt(pai))
	  dab=0.0529*beta*wnk*dab3*dab2*dab1*hrms**5
!     cab=dab/(sig*(1.0/8.0)*hrms**2)
	  cab=8.0*dab/sig/hrms**2
	  goto 20
 
   10 CAB=0.0

20    continue
      
      RETURN
      END SUBROUTINE WVBRK3_inline

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!  subroutine to obtain wave breaking energy dissipation term
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE WVBRK1_inline(CAB,JBV,II,JM,JJ,fc)
!---------------------------------------------------------------
!  energy dissipation term is formulated by using a modified
!    Miche's breaker index in order to include current effects
!  parameters used in this subroutine are calculated 
!    in subroutine veloc
!---------------------------------------------------------------
      USE GLOBAL_inline, ONLY: NPF,MPD,IPMX,JPMX,IGPX,JGPX
      COMMON /DATD/DX,DY,DXX,DMESH,DTH,kdate,idate,depmax,depmax0
      COMMON /DATG/IBK,DBAR,WL0,WCC(JGPX),HSB(JGPX),HSG(JGPX),DCD(JGPX)
      COMMON /BREK/DEPM(JGPX),DMNJ(JGPX),SLF(JGPX),wlmn(JGPX),cmn(jgpx)  &
                   ,sigm(jgpx),IWVBK,ibk3
      COMMON /VPAI/PAI2,PAI,HPAI,RAD,akap,imod,iprp,island,imd,iprpp     &
                   ,nonln,igrav,isolv,ixmdf,iproc,imud,iwnd,depmin0

! -- skip of large wavenumber
      if(wlmn(jm).le.0.0) go to 10
! -- skip of large wavenumber

! -- ibk=1 and jbv>=6 means no wave breaking
      IF(IBK.EQ.1) GOTO 10
      IF(JBV.GE.6) GOTO 10
!
      gmh=5.
      ALF=1.6
      SLJ=SLF(JM)
      IF(SLJ.LE.0) GOTO 10
      if(slj.ge.0.04) go to 10
      if(slj.le.0.1) gm=0.8+5.0*slj
      if(slj.gt.0.1) gm=1.3
!     gmh=min(gm*pai2*depm(jm)/(0.88*wlmn(jm)),10.)
      gmh=min(7.14*gm*depm(jm)/wlmn(jm),10.)
      HB=0.14*wlmn(jm)*tanh(gmh)
!     DLHB=0.14*slj*dxx*(gm/0.88)*pai/cosh(gmh)**2
      DLHB=0.5*slj*dxx*gm/cosh(gmh)**2
      H13=HSB(JJ)
      X1=(ALF*(HB+DLHB)/H13)**2
      X2=(ALF*(HB-DLHB)/H13)**2
      IF(X1.GT.25.) GOTO 30
!     PX1=PAI*X1/4.0
      PX1=.7853982*X1
      PR1=(1.0+PX1)/EXP(PX1)
      GOTO 40
   30 PR1=0.0
   40 IF(X2.GT.25.) GOTO 50
!     PX2=PAI*X2/4.0
      PX2=.7853982*X2
      PR2=(1.0+PX2)/EXP(PX2)
      GOTO 60
   50 PR2=0.0
!  60 CAB=beta*(PR2-PR1)/(1.0-PR1)*CMN(JM)/DXX
   60 CAB=(PR2-PR1)/(1.0-PR1)*CMN(JM)/DXX
!  60 CAB=(PR2-PR1)/(1.-PR1)*fc
      GOTO 20

   10 CAB=0.0

   20 continue
      RETURN
      END


!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!  FINDING REFLECTION COEFF.
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      FUNCTION CKR_inline(II,JJ)
!  ****   FINDING REFLECTION COEFF. CALLED BY SETAB
      USE GLOBAL_INLINE, ONLY: KOMX,NOMX,IPMX,JPMX
      COMMON /REFA/KRMX,KR(2,6*IPMX),RK(6*IPMX),yangl(6*IPMX)
      K=1
   20 IF(KR(1,K).EQ.II) GOTO 10
   40 K=K+1
      IF(K.LE.KRMX) GOTO 20
      CKR_inline=0
      GOTO 30
   10 IF(KR(2,K).NE.JJ) GOTO 40
      CKR_inline=RK(K)
   30 RETURN
      END FUNCTION


!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!  FINDING structure angle
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      FUNCTION aKR_inline(II,JJ)
!  ****   FINDING structure angle CALLED BY SETAB
      USE GLOBAL_INLINE, ONLY: KOMX,NOMX,IPMX,JPMX
      COMMON /REFA/KRMX,KR(2,6*IPMX),RK(6*IPMX),yangl(6*IPMX)
      K=1
   20 IF(KR(1,K).EQ.II) GOTO 10
   40 K=K+1
      IF(K.LE.KRMX) GOTO 20
      aKR_inline=0
      GOTO 30
   10 IF(KR(2,K).NE.JJ) GOTO 40
      aKR_inline=yangl(K)
   30 RETURN
      END FUNCTION


!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!  CALCULATION OF Matrix BY GAUS-SEIDEL METHOD
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE GSM_inline(II,JB,JE,NMX,MARK)
      USE GLOBAL_INLINE, ONLY: NPF,MPD,IPMX,JPMX,IGPX,JGPX,MPMX
      COMMON /DATA/NF,MD,IMAX,JMAX,IGMX,JGMX,JCPB,JCPE,JCB,JCE,NFF,MDD
      COMMON /DATB/ICK3,ICK4,ICHOICE,HS0,wd,ws,ph0,wdd(mpd),aslop(ipmx)
      COMMON /DATD/DX,DY,DXX,DMESH,DTH,kdate,idate,depmax,depmax0
      COMMON /SPECA/SCP(JGPX,MPD),SI(JGPX,NPF,MPD)
      COMMON /AB/AA(5,MPMX),IA(5,MPMX),B(MPMX),X(MPMX)
      COMMON /VPAI/PAI2,PAI,HPAI,RAD,akap,imod,iprp,island,imd,iprpp     &
                   ,nonln,igrav,isolv,ixmdf,iproc,imud,iwnd,depmin0
!     dimension x0(mpmx)
      REAL,ALLOCATABLE :: x0(:)
      ALLOCATE (x0(mpmx))
!
!     ifast=0
!     if(iprpp.ne.-1.and.dxx.gt.300.) ifast=1
!
      LIM=1000
      DLTA=0.001
      if(dmesh.le..5) then
        LIM=10000
        DLTA=.000001
      end if
      ws1000=ws*1000.
!
      JUDG=0
      XXMAX=0
      DO 20 I=1,NMX
        CA=AA(1,I)
        IF(ABS(CA).LT.1.0E-4) GOTO 21
        XX=B(I)/CA
        GOTO 22
   21   XX=0
   22   IF(XX.GT.1.0E-20) JUDG=1
        X(I)=XX
        x0(i)=xx
        ABXX=ABS(XX)
        IF(XXMAX.LT.ABXX) XXMAX=ABXX
   20 CONTINUE
      IF(JUDG.EQ.0) GOTO 11
      XLIM=XXMAX/1000.0
      ICC=0
   70 ICC=ICC+1
      XMAX=0.0
      DO 10 I=1,NMX
        IF(ABS(AA(1,I)).LT.1.0E-4) GOTO 10
        IP=0
        S=0.
        DO 30 J=1,5
          VAA=AA(J,I)
          IF(ABS(VAA).LT.1.0E-4) GOTO 30
          K=IA(J,I)
          IF(J.EQ.1) GOTO 40
          IF(I.EQ.K) GOTO 30
          IF (K.GT.0) THEN
            if(x(k).gt.1.0E-18) then
              if(abs(s).gt.100.) s=0.
              S=S+VAA*X(K)
            endif
          ENDIF
          GOTO 30
   40       IF(I.NE.K) THEN
            MARK=1
            GOTO 50
            END IF
          CA=VAA
          IP=1
   30   CONTINUE
          IF(IP.EQ.0) THEN
          MARK=1
          GOTO 50
          END IF
        PX=X(I)
        XX=(B(I)-S)/CA
!       if(ws.ge..1) then
!         if(ifast.eq.1) then
!           if(abs(xx-px).gt.x0(i)*ws1000) xx=px
!         end if
!       end if
        X(I)=XX
        IF(ABS(PX).LT.XLIM) GOTO 10
        PPX=(XX-PX)/PX
        PPX=ABS(PPX)
        IF(XMAX.LT.PPX) XMAX=PPX
   10 CONTINUE
        IF(ICC.GT.LIM) THEN
        MARK=2
        GOTO 60
        END IF
      IF(XMAX.GT.DLTA) GOTO 70
   11 MX=0
      DO 80 J=JB,JE
        DO 80 M=imd,MD
          MX=MX+1
          if(x(mx).gt.5000.) x(mx)=0.
          SCP(J,M)=X(MX)
   80 CONTINUE
      GOTO 90
   50 write(*,*) 'GSM TYPE 1:',II
      go to 90
   60 write(*,*) 'GSM TYPE 2:',II,ICC,XMAX
!
   90 DEALLOCATE (x0)
   
      RETURN
      END

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!  subroutine to output the results 
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE OUTFILE_inline
      USE GLOBAL_INLINE, ONLY: NPF,MPD,MPD2,IPMX,JPMX,KOMX,NOMX,IGPX,JGPX
      double precision edate
      REAL, ALLOCATABLE :: SPCT(:),SPCT1(:),dep0(:,:),dsss(:,:),tmp(:,:)
      REAL, ALLOCATABLE :: gxx1(:), gyy1(:)
!     REAL, ALLOCATABLE :: gxr(:,:),gyr(:,:)
      
      !Change the following variables to allocatable and then allocate later   MEB  11/16/21
      real, allocatable :: sp1(:,:),aship(:,:),eship(:,:)
      real, allocatable :: xship1(:),yship1(:),xship2(:),yship2(:)
      !Dimension sp1(npf,mpd),aship(410,410),eship(410,410)
      !Dimension xship1(mpd),yship1(mpd),xship2(mpd),yship2(mpd)
      
      
      COMMON /VPAI/PAI2,PAI,HPAI,RAD,akap,imod,iprp,island,imd,iprpp     &
                   ,nonln,igrav,isolv,ixmdf,iproc,imud,iwnd,depmin0
      COMMON /DATA/NF,MD,IMAX,JMAX,IGMX,JGMX,JCPB,JCPE,JCB,JCE,NFF,MDD
      COMMON /DATB/ICK3,ICK4,ICHOICE,HS0,wd,ws,ph0,wdd(mpd),aslop(ipmx)
      COMMON /DATC/TP,PRD,IBND,VL,TIDE,KPEAK,IBACK,NBLOCK,IWIND,WSMAG
      COMMON /DATD/DX,DY,DXX,DMESH,DTH,kdate,idate,depmax,depmax0
      COMMON /DVAR/DVARX(IGPX),DVARY(JGPX),ETA(IGPX,JGPX),HSK(JGPX)
      COMMON /SPECB/SOP(KOMX,NPF,MPD2),SR(2*IPMX,NPF,MPD)
      COMMON /OUTP/KOUT,IJSP(2,KOMX),IRS,IBREAK,ICUR,IWET,INST
      COMMON /OUTN/nest,inest(komx),jnest(komx),ix1(igpx),ix2(igpx)
      COMMON /IJBP/IJB(IGPX,JGPX),DEP(IPMX,JPMX),DEPS(IPMX),DBIG(IPMX)
      COMMON /WAVI/H13(IGPX,JGPX),T13(IGPX,JGPX),DMN(IGPX,JGPX)
      COMMON /WAVR/H13R(IGPX,JGPX),T13R(IGPX,JGPX),DMNR(IGPX,JGPX)
      COMMON /WAVS/H13S(IGPX,JGPX),IBR(IGPX,JGPX),DISS(IGPX,JGPX)
      COMMON /WAVF/H13F(IGPX,JGPX),DMNF(IGPX,JGPX),T13F(IGPX,JGPX),T13MIN
      COMMON /FDS/FCN(NPF),DCM(MPD),DSFD(NPF,MPD)
      COMMON /DFCN/DF(NPF),FFCN(NPF),DFINP(NPF),PL1E(NPF),PL1T(NPF)
      common /dxxdyy/dvarxx(ipmx),dvaryy(jpmx)                   !Wu
      common /rsb/ wxrs(ipmx,jpmx),wyrs(ipmx,jpmx)
      common /rsm/ cosa(mpd),sina(mpd),d1(ipmx,jpmx)
      common /rsz/ igetfile20,nship,disx(ipmx),disy(jpmx)
      common /ship/ShipL(30),ShipB(30),ShipD(30),ShipS(30)
      common /uv1/u1(ipmx,jpmx),v1(ipmx,jpmx)
      common /struc0/ijstruc1,ijstruc2,ijstruc3,ijstruc4,ismall
      common /struc1/istruc1(komx),jstruc1(komx),dstruc1(komx)
      common /struc2/istruc2(NOMX),jstruc2(NOMX),dstruc2(NOMX)
      common /struc3/istruc3(komx),jstruc3(komx),dstruc3(komx) &
                    ,k3(komx),dstruc33(komx)
      common /struc4/istruc4(komx),jstruc4(komx),dstruc4(komx) &
                    ,kstruc4(komx),k4(komx),dstruc44(komx)
      common /swell/sw13(igpx,jgpx),sa13(igpx,jgpx) &
                   ,tw13(igpx,jgpx),ta13(igpx,jgpx) &
                   ,dw13(igpx,jgpx),da13(igpx,jgpx)
      common /origin/x0,y0,azimuth,isteer,iidate,sinaz,cosaz
      common /FileNames/ OptsFile, DepFile, CurrFile, EngInFile,    &
                         WaveFile, ObsFile, EngOutFile, NestFile,   &
                         BreakFile, RadsFile, StrucFile, SurgeFile, &
                         MudFile, FricFile, FrflFile, BrflFile,     &
                         SpecFile, WindFile, XMDFFile, SetupFile,   &  !Mitch 3/22/2017
                         SeaFile, SwellFile, ShipFile                  !Mitch 3/22/2017
      common /comsav/depmin,g,iview,iplane,iwave,cflat,irunup,iroll                         
      logical getfile22
      CHARACTER*180 WaveFile,ObsFile,EngOutFile,BreakFile,RadsFile
      CHARACTER*180 OptsFile,DepFile,CurrFile,EngInFile
      CHARACTER*180 NestFile,StrucFile,SurgeFile
      CHARACTER*180 MudFile,FricFile,FrflFile,BrflFile,WindFile
      CHARACTER*180 SpecFile, XMDFFile,SetupFile,ShipFile
      CHARACTER*180 SeaFile, SwellFile, TotalFile                      !Mitch 3/22/2017
      ALLOCATE (SPCT(MPD2),SPCT1(MPD2),dep0(ipmx,jpmx),dsss(ipmx,jpmx))
      ALLOCATE (gxx1(ipmx),gyy1(jpmx),tmp(ipmx,jpmx))
!     ALLOCATE (gxr(ipmx,jpmx),gyr(ipmx,jpmx))
!
      g=9.806
      ammd=float(mdd-1)/float(md-1)
!
      if(ick4.eq.0) then
        mddd=mdd
      else
        mddd=mdd*2
      end if
        
      allocate(sp1(npf,mddd))         !This variable normally is of dimension: npf, mpd.  That is too small in some cases with backward reflection.
      allocate(aship(410,410),eship(410,410))
      allocate(xship1(mpd),yship1(mpd),xship2(mpd),yship2(mpd))
        
!
      wdrad=wd
      if(iprp.ne.1) wdrad=wd/rad
      IF(nest.GE.1) THEN
        DO 920 KK=1,nest
          kn=kk+iabs(kout)
          i=inest(kk)
          j=jnest(kk)
!
          xdis=0.
          ydis=0.
          do ki=1,i-1
            xdis=xdis+dvarxx(ki)
          end do
          xdis=xdis+dvarxx(i)/2.
          do kj=1,j-1
            ydis=ydis+dvaryy(kj)
          end do
          ydis=ydis+dvaryy(j)/2.
          xc=xdis*cosaz-ydis*sinaz+x0
          yc=xdis*sinaz+ydis*cosaz+y0
!
          if(imod.eq.0) then
            fpp=1./t13(i,j)
            h13ij=h13(i,j)
          end if
          if(imod.eq.1) then
            fpp=1./t13(j,jgmx-i+1)
            h13ij=h13(j,jgmx-i+1)
          end if
          if(imod.eq.2) then
            fpp=1./t13(igmx-i+1,jgmx-j+1)
            h13ij=h13(igmx-i+1,jgmx-j+1)
          end if
          if(imod.eq.3) then
            fpp=1./t13(igmx-j+1,i)
            h13ij=h13(igmx-j+1,i)
          end if
!
9013  format(i8,i5,2f8.2,f8.4,f8.2,2f15.2,f8.3)
9014  format(5x,i8,2f8.2,f8.4,f8.2,2f15.2,f8.3)
9023  format(i7,i6,2f8.2,f8.4,f8.2,2f15.2,f8.3)
          if(kdate.gt.0) then
            if(idate.le.9999) then
              idate1=mod(kdate,10)*100000+idate
              kdate1=kdate/10
              if(iplane.le.1) then
                WRITE(13,9023) kdate1,idate1,ws,wdrad,fpp,Tide,xc,yc,h13ij
                if(inst.eq.1)WRITE(14,9023) kdate1,idate1,ws,wdrad,fpp,Tide,xc,yc,h13ij
              else
                WRITE(30,9023) kdate1,idate1,ws,wdrad,fpp,Tide,xc,yc,h13ij
              end if
            else
              if(iplane.le.1) then
                WRITE(13,9013) kdate,idate,ws,wdrad,fpp,Tide,xc,yc,h13ij
                if(inst.eq.1)WRITE(14,9013) kdate,idate,ws,wdrad,fpp,Tide,xc,yc,h13ij
              else
                WRITE(30,9013) kdate,idate,ws,wdrad,fpp,Tide,xc,yc,h13ij
              end if
            end if
          else
            if(iplane.le.1) then
              WRITE(13,9014) idate,ws,wdrad,fpp,Tide,xc,yc,h13ij
              if(inst.eq.1)WRITE(14,9014) idate,ws,wdrad,fpp,Tide,xc,yc,h13ij
            else
              WRITE(30,9014) idate,ws,wdrad,fpp,Tide,xc,yc,h13ij
            end if
          end if
!
          n00=0
          DO 935 n0=1,NFF
          if(ffcn(n0).lt.fcn(1)) then
!           if(iplane.le.1) then
!           write(13,510) (0.,mm=1,mdd)
!     if(inst.eq.1) write(14,510) (0.,mm=1,mdd)
!           else
!                   write(30,510) (0.,mm=1,mdd)
!
!          end if
          n00=n00+1
          do mm=1,mddd           !changed this from mdd.   Size of mddd depends on if backward reflection is enabled.   MEB  11/16/2021
             sp1(n00,mm)=0.
          enddo
          go to 935
          end if
!
          if(ffcn(n0).eq.fcn(1)) then
            IF(ICHOICE.LE.1) THEN
!             if(iplane.le.1) then
!               WRITE(13,510) (SOP(Kn,1,MM),MM=1,MDD)
!               if(inst.eq.1) WRITE(14,510) (SOP(Kn,1,MM),MM=1,MDD)
!             else
!               WRITE(30,510) (SOP(Kn,1,MM),MM=1,MDD)
!             end if
              n00=n00+1
              do mm=1,mdd
                sp1(n00,mm)=sop(kk,1,mm)
              enddo
            ELSE
              c0=float(mdd+1)/float(md+1)
              DO mm=1,mdd
                m1=int(float(mm-1)/ammd)+1
                m2=m1+1
                if(m2.gt.md) m2=md
                c3=(float(m1)*ammd-float(mm-1))/ammd
                c4=c0*(1.-c3)
                c3=c0*c3
                SPCT(MM)=SOP(Kn,1,M1)*c3+SOP(Kn,1,M2)*c4
              END DO
!             if(iplane.le.1) then
!               WRITE(13,510) (SPCT(MM),MM=1,MDD)
!               if(inst.eq.1) WRITE(14,510) (SPCT(MM),MM=1,MDD)
!             else
!               WRITE(30,510) (SPCT(MM),MM=1,MDD)
!             end if
              n00=n00+1
              do mm=1,mdd
                sp1(n00,mm)=spct(mm)
              enddo
            END IF
            go to 935
          end if
!
          if(ffcn(n0).gt.fcn(nf)) then
          c1=(ffcn(nff)-ffcn(n0))/(ffcn(nff)-fcn(nf))
          IF(ICHOICE.LE.1) THEN
!             if(iplane.le.1) then
!             WRITE(13,510) (SOP(Kn,NF,MM)*c1,MM=1,MDD)
!         if(inst.eq.1) WRITE(14,510) (SOP(Kn,NF,MM)*c1,MM=1,MDD)
!             else
!                       WRITE(30,510) (SOP(Kn,NF,MM)*c1,MM=1,MDD)
!             end if
              n00=n00+1
              do mm=1,mdd
                sp1(n00,mm)=sop(kk,nf,mm)*c1
              enddo
          ELSE
            c0=float(mdd+1)/float(md+1)
            DO mm=1,mdd
            m1=int(float(mm-1)/ammd)+1
            m2=m1+1
            if(m2.gt.md) m2=md
            c3=(float(m1)*ammd-float(mm-1))/ammd
            c4=c0*(1.-c3)
            c3=c0*c3
            SPCT(MM)=SOP(Kn,NF,M1)*c3+SOP(Kn,NF,M2)*c4
            END DO
!             if(iplane.le.1) then
!             WRITE(13,510) (SPCT(MM)*c1,MM=1,MDD)
!         if(inst.eq.1) WRITE(14,510) (SPCT(MM)*c1,MM=1,MDD)
!             else
!                       WRITE(30,510) (SPCT(MM)*c1,MM=1,MDD)
!             end if
              n00=n00+1
              do mm=1,mdd
                sp1(n00,mm)=spct(mm)*c1
              enddo
          END IF
          go to 935
          end if
!
            DO 930 NN=1,NF-1
            nn1=nn+1
            if(ffcn(n0).gt.fcn(nn).and.ffcn(n0).le.fcn(nn1)) then
            c0=fcn(nn1)-fcn(nn)
            c2=(ffcn(n0)-fcn(nn))/c0
            c1=(fcn(nn1)-ffcn(n0))/c0
            IF(ICHOICE.LE.1) THEN
!             if(iplane.le.1) then
!             WRITE(13,510) (SOP(Kn,NN,MM)*c1+SOP(Kn,NN1,MM)*c2,MM=1,MDD)
!     if(inst.eq.1)WRITE(14,510)(SOP(Kn,NN,MM)*c1+SOP(Kn,NN1,MM)*c2,MM=1,MDD)
!             else
!                  WRITE(30,510)(SOP(Kn,NN,MM)*c1+SOP(Kn,NN1,MM)*c2,MM=1,MDD)
!             end if
              n00=n00+1
              do mm=1,mdd
                sp1(n00,mm)=sop(kk,nn,mm)*c1+sop(kk,nn1,mm)*c2
              enddo
            ELSE
            c0=float(mdd+1)/float(md+1)
            DO mm=1,mdd
            m1=int(float(mm-1)/ammd)+1
            m2=m1+1
            if(m2.gt.md) m2=md
            c3=(float(m1)*ammd-float(mm-1))/ammd
            c4=c0*(1.-c3)
            c3=c0*c3
            SPCT(MM)=SOP(Kn,NN,M1)*c3+SOP(Kn,NN,M2)*c4
            SPCT1(MM)=SOP(Kn,NN1,M1)*c3+SOP(Kn,NN1,M2)*c4
            END DO
!             if(iplane.le.1) then
!             WRITE(13,510) (SPCT(MM)*c1+SPCT1(MM)*c2,MM=1,MDD)
!     if(inst.eq.1) WRITE(14,510)(SPCT(MM)*c1+SPCT1(MM)*c2,MM=1,MDD)
!             else
!                   WRITE(30,510)(SPCT(MM)*c1+SPCT1(MM)*c2,MM=1,MDD)
!             end if
              n00=n00+1
              do mm=1,mdd
                sp1(n00,mm)=spct(mm)*c1+spct1(mm)*c2
              enddo
            END IF
            go to 935
            end if
  930       CONTINUE
!
  935     CONTINUE
!
          sum=0.
          sum1=0.
!         sum2=0.
          do mm=1,mdd
          sum1=sum1+sp1(nff,mm)
          enddo
!
          do nn=nff,2,-1
            sum0=0.
            do mm=1,mdd
            sum0=sum0+sp1(nn-1,mm)
            enddo
            if(sum0.eq.0..and.sum1.gt.0.) then
              do mm=1,mdd
              sp1(nn-1,mm)=sp1(nn,mm)*.3
              sp1(nn,mm)=sp1(nn,mm)*.7
              enddo
              sum1=sum1*.7
              sum0=sum1*.3
            endif
            cc=sum1*dfinp(nn)
            sum=sum+cc
            sum1=sum0
!           sum2=sum2+cc*ffcn(nn)
          end do
            cc=sum1*dfinp(1)
            sum=sum+cc
!           sum2=sum2+cc*ffcn(1)
!         t13fij=sum/sum2
          h13fij=4.*sqrt(sum*pai/float(1+mdd))

!     write(*,*) 'ok13,h13ij,h13fij,t13fij,nff,mdd=',h13ij,h13fij,t13fij,nff,mdd
!     write(*,*) 'dfinp(1:5)=',dfinp(1),dfinp(2),dfinp(3),dfinp(4),dfinp(5)
!
          c1=1.
          if(h13fij.gt..001) then
          c1=(h13ij/h13fij)**2
          end if
!
        do nn=1,nff
            if(iplane.le.1) then
            WRITE(13,510) (sp1(nn,mm)*c1,mm=1,mdd)
          if(inst.eq.1) WRITE(14,510) (sp1(nn,mm)*c1,mm=1,mdd)
            else
                        WRITE(30,510) (sp1(nn,mm)*c1,mm=1,mdd)
            end if
        enddo
!
  920   CONTINUE
      END IF
!
! -- total significant wave height
!
      DO 81 J=1,JGMX
        DO 81 I=1,IGMX
          if(ick4.eq.1) then
            cx=h13(i,j)*cos(dmn(i,j)*rad)+h13r(i,j)*cos(dmnr(i,j)*rad)
            cy=h13(i,j)*sin(dmn(i,j)*rad)+h13r(i,j)*sin(dmnr(i,j)*rad)
            dmn(i,j)=atan2(cy,cx)/rad
          end if
!
            H13S(I,J)=H13(I,J)+H13R(I,J)
!
          if(iplane.eq.1) then
            H13f(I,J)=H13S(I,J)
            dmnf(I,J)=dmn(I,J)*rad
            t13f(I,J)=t13(I,J)
          elseif(iplane.eq.2) then
            imax1=imax-i
            jmax1=jmax-j
            cc1=h13s(i,j)**2
            cc2=h13f(imax1,jmax1)**2
            H13S(I,j)=SQRT(cc1+cc2)
            cx=cc1*cos(dmn(i,j)*rad)+cc2*cos(dmnf(imax1,jmax1)+pai)
            cy=cc1*sin(dmn(i,j)*rad)+cc2*sin(dmnf(imax1,jmax1)+pai)
            dmn(i,j)=atan2(cy,cx)/rad
            T13(i,j)=(T13(i,j)*cc1+T13f(imax1,jmax1)*cc2)/(cc1+cc2+1.e-10)
          end if
!
!      if(iwind.eq.2.and.d1(i,j).eq..9995) then
!      h13s(i,j)=0.
!      t13(i,j)=1.
!      end if
!
!         IF(H13S(I,J).LT..01) T13(I,J)=1.
!         if(depmax.lt.1.) go to 85
!         IF(H13S(I,J).LT..1.AND.T13(I,J).LT.1.) T13(I,J)=1.
!  85   continue
!     write(20,*) i,j,h13s(i,j)
   81 CONTINUE
!
!     if(iplane.eq.2) idate=idate+1
!
      if(nblock.eq.1) then
        open(unit=19,file='block.dat',status='old')
        read(19,*) kblock
        do 309 k=1,kblock
          read(19,*) iblock,jblock
          if(iblock.eq.1) H13S(IBLOCK,JBLOCK)=0.
309     continue
        close(19)
      end if
!
!     Special output for imod = 0 and simgrid.dat exists
!
      if(imod.eq.0) then
!
!     Use ksim=0 below without simgrid.dat for the entire grid 
      ksim=-1
!
      do ii=1,igmx
        if(ii.eq.1) then
          gxx1(1)=dvarx(1)/2.
        else
          gxx1(ii)=gxx1(ii-1)+(dvarx(ii)+dvarx(ii-1))/2.
        end if
      end do
!
      do jj=1,jgmx
        if(jj.eq.1) then
          gyy1(1)=dvary(1)/2.
        else
          gyy1(jj)=gyy1(jj-1)+(dvary(jj)+dvary(jj-1))/2.
        end if
      end do
!
      inquire(file='simgrid.dat',exist=getfile22)
      if(getfile22) then
        open(unit=35,file='simgrid.dat',status='old')
          read(35,*,end=190,err=190) x00,y00,az0
          read(35,*,end=190,err=190) isim,jsim,smesh
          ksim=1
  190   continue
        close(35)
      end if
!
      if(iplane.eq.1) go to 430
      if(ksim.eq.1) then
      open(25,file='wave.pts',status='unknown')
      if(kdate.gt.0) then
        if(idate.le.9999) then
          idate1=mod(kdate,10)*100000+idate
          kdate1=kdate/10
          write(25,'(i7,i6)') kdate1,idate1
        else
          write(25,'(i8,i5)') kdate,idate
        end if
      else
      write(25,'(i8)') idate
      end if
!
      cosaz0=cos(az0*rad)
      sinaz0=sin(az0*rad)
      az1=az0-azimuth
      cosaz1=cos(az1*rad)
      sinaz1=sin(az1*rad)
!
        do is=1,isim
        gsi=smesh*(float(is)-0.5)
          do js=1,jsim
          gsj=smesh*(float(js)-0.5)
          gxs=gsi*cosaz0-gsj*sinaz0+x00
          gys=gsi*sinaz0+gsj*cosaz0+y00
!     gxs,gys are true coordinates
!
          gxs1=(gsi+x00-x0)*cosaz1-(gsj+y00-y0)*sinaz1
          gys1=(gsi+x00-x0)*sinaz1+(gsj+y00-y0)*cosaz1
!     gxs1,gys1 are projection to local (gxx1,gyy1) plane
!
            do jj=1,jgmx-1
            do ii=1,igmx-1
            ii1=ii+1
            jj1=jj+1
            cc1=(gxs1-gxx1(ii))*(gxs1-gxx1(ii1))
            cc2=(gys1-gyy1(jj))*(gys1-gyy1(jj1))
            if(cc1.le..0.and.cc2.le..0) goto 112
            end do
            end do
!
  112     continue
          a1=abs((gxs1-gxx1(ii))*(gys1-gyy1(jj)))
          a2=abs((gxs1-gxx1(ii1))*(gys1-gyy1(jj)))
          a3=abs((gxs1-gxx1(ii))*(gys1-gyy1(jj1)))
          a4=abs((gxs1-gxx1(ii1))*(gys1-gyy1(jj1)))
          hh13=0.
          tt13=1.
          dmn13=0.
          if(h13s(ii,jj).lt..001.and.t13(ii,jj).lt.1.01) a4=0.
          if(h13s(ii,jj1).lt..001.and.t13(ii,jj1).lt.1.01) a2=0.
          if(h13s(ii1,jj).lt..001.and.t13(ii1,jj).lt.1.01) a3=0.
          if(h13s(ii1,jj1).lt..001.and.t13(ii1,jj1).lt.1.01) a1=0.
          a0=a1+a2+a3+a4
          if(a0.gt..01) then
          hh13=(h13s(ii,jj)*a4+h13s(ii,jj1)*a2    &
          +h13s(ii1,jj)*a3+h13s(ii1,jj1)*a1)/a0
          tt13=(t13(ii,jj)*a4+t13(ii,jj1)*a2    &
          +t13(ii1,jj)*a3+t13(ii1,jj1)*a1)/a0
          uu13=h13s(ii,jj)*cos((dmn(ii,jj)-az1)*rad)*a4   &
          +h13s(ii,jj1)*cos((dmn(ii,jj1)-az1)*rad)*a2    &
          +h13s(ii1,jj)*cos((dmn(ii1,jj)-az1)*rad)*a3    &
          +h13s(ii1,jj1)*cos((dmn(ii1,jj1)-az1)*rad)*a1
          vv13=h13s(ii,jj)*sin((dmn(ii,jj)-az1)*rad)*a4   &
          +h13s(ii,jj1)*sin((dmn(ii,jj1)-az1)*rad)*a2    &
          +h13s(ii1,jj)*sin((dmn(ii1,jj)-az1)*rad)*a3    &
          +h13s(ii1,jj+1)*sin((dmn(ii1,jj+1)-az1)*rad)*a1
          dmn13=atan2(vv13,uu13)/rad
          end if
          write(25,'(2f12.3,f7.2,2f7.1)') gxs,gys,hh13,tt13,dmn13
!
          end do
        end do
      end if
!
!
      if(ksim.eq.0) then
      open(25,file='wave.pts',status='unknown')
      if(kdate.gt.0) then
        if(idate.le.9999) then
          idate1=mod(kdate,10)*100000+idate
          kdate1=kdate/10
          write(25,'(i7,i6)') kdate1,idate1
        else
          write(25,'(i8,i5)') kdate,idate
        end if
      else
      write(25,'(i8)') idate
      end if
!
      do jj=1,jgmx
      do ii=1,igmx
      gxrr=gxx1(ii)*cosaz-gyy1(jj)*sinaz+x0
      gyrr=gxx1(ii)*sinaz+gyy1(jj)*cosaz+y0
      write(25,'(2f12.3,f7.2,2f7.1)') &
      gxrr,gyrr,h13s(ii,jj),t13(ii,jj),dmn(ii,jj)
      end do
      end do
!
      end if
!
!
!     if(ksim.eq.1) then
!
!     gxx=0.
!     do ii=1,igmx
!     if(ii.eq.1) then
!     gxx=gxx+dvarx(ii)/2.
!     else
!     gxx=gxx+(dvarx(ii)+dvarx(ii-1))/2.
!     end if
!     gyy=0.
!     do jj=1,jgmx
!     if(jj.eq.1) then
!     gyy=gyy+dvary(jj)/2.
!     else
!     gyy=gyy+(dvary(jj)+dvary(jj-1))/2.
!     end if
!     gxr(ii,jj)=gxx*cosaz-gyy*sinaz+x0-x00
!     gyr(ii,jj)=gxx*sinaz+gyy*cosaz+y0-y00
!     end do
!     end do
!
!     open(25,file='wave.pts',status='unknown')
!       if(kdate.le.9999) then
!       idate1=mod(kdate,10)*100000+idate
!       kdate1=kdate/10
!         write(25,'(i7,i6)') kdate1,idate1
!       else
!         write(25,'(i8,i5)') kdate,idate
!       end if
!     cosaz0=cos(az0*rad)
!     sinaz0=sin(az0*rad)
!       do is=1,isim
!       gsi=smesh*(float(is)-0.5)
!         do js=1,jsim
!         gsj=smesh*(float(js)-0.5)
!         gxs=gsi*cosaz0-gsj*sinaz0
!         gys=gsi*sinaz0+gsj*cosaz0
!
!         do jj=1,jgmx-1
!         do ii=1,igmx-1
!         ii1=ii+1
!         jj1=jj+1
!         a0=(dvarx(ii)+dvarx(ii1))*(dvary(jj)+dvary(jj1))/4.
!         a1=abs(gxs*(gyr(ii1,jj)-gyr(ii,jj))+gys*(gxr(ii,jj)  &
!         -gxr(ii1,jj))+gxr(ii1,jj)*gyr(ii,jj)-gxr(ii,jj)*gyr(ii1,jj))/2.
!         a2=abs(gxs*(gyr(ii,jj1)-gyr(ii,jj))+gys*(gxr(ii,jj)  &
!         -gxr(ii,jj1))+gxr(ii,jj1)*gyr(ii,jj)-gxr(ii,jj)*gyr(ii,jj1))/2.
!         a3=abs(gxs*(gyr(ii,jj1)-gyr(ii1,jj1))+gys*(gxr(ii1,jj1)  &
!         -gxr(ii,jj1))+gxr(ii,jj1)*gyr(ii1,jj1)-gxr(ii1,jj1)*gyr(ii,jj1))/2.
!         a4=abs(gxs*(gyr(ii1,jj)-gyr(ii1,jj1))+gys*(gxr(ii1,jj1)  &
!         -gxr(ii1,jj))+gxr(ii1,jj)*gyr(ii1,jj1)-gxr(ii1,jj1)*gyr(ii1,jj))/2.
!         if(abs(a0-a1-a2-a3-a4).lt.smesh/5.) goto 111
!         end do
!         end do
!
! 111   continue
!       hsmall=0.
!       tsmall=1.
!       weight=0.
!       if(h13s(ii,jj).gt..001.and.t13(ii,jj).gt.1.) then
!       cc3=1./(sqrt((gxs-gxr(ii,jj))**2+(gys-gyr(ii,jj))**2)+.1)
!       weight=weight+cc3
!       hsmall=hsmall+h13s(ii,jj)*cc3
!       tsmall=tsmall+t13(ii,jj)*cc3
!       end if
!       if(h13s(ii1,jj).gt..001.and.t13(ii1,jj).gt..001) then
!       cc3=1./(sqrt((gxs-gxr(ii1,jj))**2+(gys-gyr(ii1,jj))**2)+.1)
!       weight=weight+cc3
!       hsmall=hsmall+h13s(ii1,jj)*cc3
!       tsmall=tsmall+t13(ii1,jj)*cc3
!       end if
!       if(h13s(ii,jj1).gt..001.and.t13(ii,jj1).gt..001) then
!       cc3=1./(sqrt((gxs-gxr(ii,jj1))**2+(gys-gyr(ii,jj1))**2)+.1)
!       weight=weight+cc3
!       hsmall=hsmall+h13s(ii,jj1)*cc3
!       tsmall=tsmall+t13(ii,jj1)*cc3
!       end if
!       if(h13s(ii1,jj1).gt..001.and.t13(ii1,jj1).gt..001) then
!       cc3=1./(sqrt((gxs-gxr(ii1,jj1))**2+(gys-gyr(ii1,jj1))**2)+.1)
!       weight=weight+cc3
!       hsmall=hsmall+h13s(ii1,jj1)*cc3
!       tsmall=tsmall+t13(ii1,jj1)*cc3
!       end if
!       if(weight.gt.0.) then
!       hsmall=hsmall/weight
!       tsmall=tsmall/weight
!       end if
!
!       write(25,'(2f12.3,f7.2,2f7.1)') gxs+x00,gys+y00,hsmall,tsmall
!
!         end do
!       end do
!     end if
!
!     if(ksim.eq.0) then
!     open(25,file='wave.pts',status='unknown')
!       if(idate.le.9999) then
!         idate1=mod(kdate,10)*100000+idate
!         kdate1=kdate/10
!         write(25,'(i7,i6)') kdate1,idate1
!       else
!         write(25,'(i8,i5)') kdate,idate
!       end if
!     gxx=0.
!     do ii=1,igmx
!     if(ii.eq.1) then
!     gxx=gxx+dvarx(ii)/2.
!     else
!     gxx=gxx+(dvarx(ii)+dvarx(ii-1))/2.
!     end if
!     gyy=0.
!     do jj=1,jgmx
!     if(jj.eq.1) then
!     gyy=gyy+dvary(jj)/2.
!     else
!     gyy=gyy+(dvary(jj)+dvary(jj-1))/2.
!     end if
!     gxrr=gxx*cosaz-gyy*sinaz+x0
!     gyrr=gxx*sinaz+gyy*cosaz+y0
!     write(25,'(2f12.3,f7.2,2f7.1)') &
!     gxrr,gyrr,h13s(ii,jj),t13(ii,jj),dmn(ii,jj)
!     end do
!     end do
!     end if
!
      end if
!
!     End special output for imod = 0 and simgrid.dat exists
!
430   continue
!
      if(icur.ge.1) go to 291
      if(wd.ne.0.) go to 291
      if(ick4.ne.0) go to 291
      if(nblock.ne.0) go to 291
      if(ijstruc1.ne.0) go to 291
!     if(ijstruc2.ne.0) go to 291
!     if(ijstruc3.ne.0) go to 291
!     if(ijstruc4.ne.0) go to 291
!
      do 290 i=1,igmx-1
        nj2=jgmx/2
        do j=1,jgmx
          if(ix1(i).eq.1) go to 291
        end do
        h13ij=h13s(i,nj2)
        do j=1,jgmx
          h13s(i,j)=h13ij
        end do
        t13ij=t13(i,nj2)
        do j=1,jgmx
          t13(i,j)=t13ij
        end do
        dmnij=dmn(i,nj2)
        do j=1,jgmx
          dmn(i,j)=dmnij
        end do
        if(ibreak.eq.1) then
          ibrij=ibr(i,nj2)
          do j=1,jgmx
            ibr(i,j)=ibrij
          end do
        end if
        if(ibreak.ge.2) then
          dissij=diss(i,nj2)
          do j=1,jgmx
            diss(i,j)=dissij
          end do
        end if
        if(irs.ge.1) then
          wxrsij=wxrs(i,nj2)
          wyrsij=wyrs(i,nj2)
          if(i+2.le.igmx) then
            if(ix2(i).eq.1) go to 293
          end if
          do j=1,jgmx
            wxrs(i,j)=wxrsij
            wyrs(i,j)=wyrsij
          end do
293       continue
        end if
290   continue
291   continue
!
!     if(iplane.eq.1) RETURN
      if(iplane.eq.1) go to 431
!
      if (ibreak.ge.2.and.ixmdf.eq.0) then
      if(kdate.gt.0) then
        if(idate.le.9999) then
          idate1=mod(kdate,10)*100000+idate
          kdate1=kdate/10
          write(17,'(i7,i6)') kdate1,idate1
        else
          write(17,'(i8,i5)') kdate,idate
        end if
      else
      write(17,'(i8)') idate
      end if
        if(imod.eq.0) then
          do j=jgmx,1,-1
            write(17,*) (diss(i,j),i=1,igmx)
          end do
        end if
        if(imod.eq.1) then
          do i=igmx,1,-1
            write(17,*) (diss(i,j),j=jgmx,1,-1)
          end do
        end if
        if(imod.eq.2) then
          do j=1,jgmx
            write(17,*) (diss(i,j),i=igmx,1,-1)
          end do
        end if
        if(imod.eq.3) then
          do i=1,igmx
            write(17,*) (diss(i,j),j=1,jgmx)
          end do
        end if
      end if
!
431   continue
      dsss=0.
!
       if(irs.ge.2) then
        if(kout.ge.0) then
        if(iplane.ne.1.and.ixmdf.eq.0) then
          if(isteer.eq.1.and.idate.lt.10000) then
            write(98,*) iidate
          else
           if(kdate.gt.0) then
            if(idate.le.9999) then
              idate1=mod(kdate,10)*100000+idate
              kdate1=kdate/10
              write(98,'(i7,i6)') kdate1,idate1
            else
              write(98,'(i8,i5)') kdate,idate
            end if
           else
             write(98,'(i8)') idate
           end if
          end if
        end if
        end if

!
        if(iplane.ge.1) then
          open(unit=68,file='setup.dat',status='unknown')
          if(iplane.eq.1) then
            if(isteer.eq.1.and.idate.lt.10000) then
              write(68,*) iidate
            else
              write(68,*) Idate
            end if
          else
            read(68,*) edate
            idate = int(mod(edate,100000.))
            kdate = int(edate/100000.)
          end if
        end if
!
        open (15, file = DepFile,   status = 'old')
          read (15, *) ni,nj
          do j = nj, 1, -1
          read (15, *) (dep0(i, j), i = 1, ni)
          end do
        close(15)
        do k=1,ijstruc1
        dep0(istruc1(k),jstruc1(k))=dstruc1(k)
        end do
        do k=1,ijstruc2
        dep0(istruc2(k),jstruc2(k))=-dstruc2(k)
        end do
        do k=1,ijstruc4
        dep0(istruc4(k),jstruc4(k))=-dstruc44(k)
        end do
!
        do 53 j=1,jgmx
        do 53 i=1,igmx
!
        if(imod.eq.0) then
        dep00=dep0(i,j)
        elseif(imod.eq.1)then
        dep00=dep0(ni-j+1,i)
        elseif(imod.eq.2)then
        dep00=dep0(ni-i+1,nj-j+1)
        else
        dep00=dep0(j,nj-i+1)
        end if
        if(dep00.lt.-20.) dep00=-20.
!
        cc=eta(i,j)
        if(eta(i,j).eq.0..and.h13s(i,j).gt.0.001) then
        if(d1(i,j).lt.hs0) then
        eta(i,j)=h13s(i,j)/2.
        d1(i,j)=eta(i,j)
        end if
        end if
        dsss(i,j)=eta(i,j)+h13s(i,j)/2.
        if(eta(i,j).eq.0..and.h13s(i,j).lt.hs0/4.) dsss(i,j)=0.
!
        if(irs.ge.2.and.cc.gt..01) then
        if(i.eq.1.or.i.eq.igmx) goto 55
!
          eta(i,j)=cc
          if(j.eq.1) then
          if(ijb(i+1,j).eq.2) go to 55
          if(eta(i+1,j).ge..01.and.eta(i-1,j).ge..01.and.    &
          eta(i,j+1).ge..01) goto 55
          cc=2.5*eta(i,j)
          if(cc.gt.dsss(i,j)) dsss(i,j)=cc
          go to 55
          end if
!
          if(j.eq.jgmx) then
          if(ijb(i+1,j).eq.2) go to 55
          if(eta(i+1,j).ge..01.and.eta(i-1,j).ge..01.and.    &
          eta(i,j-1).ge..01) goto 55
          cc=2.5*eta(i,j)
          if(cc.gt.dsss(i,j)) dsss(i,j)=cc
          go to 55
          end if
!
!       if(ijb(i+1,j).eq.1) go to 55
!       if(eta(i+1,j).ge..01.and.eta(i-1,j).ge..01.and.      &
!       eta(i,j+1).ge..01.and.eta(i,j-1).ge..01) goto 55
        cc=2.5*eta(i,j)
!
! ----- instead of 2*eta(i,j) for 2% exceedence runup,
!       final best result is 2.5*eta(i,j) and
!       reported in selhts.out and setup.wav as the
!       maximum water level, dsss(i,j).
!
        if(h13s(i,j).le..00001) then
        dsss(i,j)=0.
        go to 55
        end if
!
        if(cc.gt.dsss(i,j)) dsss(i,j)=cc
        end if
55      continue
        if(eta(i,j)+tide+dep00.lt.0..and.d1(i,j).gt.0.) then
        eta(i,j)=d1(i,j)-dep00-tide
!       dsss(i,j)=eta(i,j)+h13s(i,j)/2.
        end if
!
        if(eta(i,j).gt.45.) eta(i,j)=45.0
        if(eta(i,j).le.-5.) eta(i,j)=-5.0
!
        if(dsss(i,j).lt.1.e-7.or.irunup.eq.2) go to 53
        eta(i,j)=eta(i,j)+tide
        dsss(i,j)=dsss(i,j)+tide
53      continue
      else
        if(irunup.eq.2) then
          dsss=eta+h13s/2.
        else
          dsss=eta+tide+h13s/2.
        end if
      end if
!
      if(irs.ge.2) then
      if(iplane.eq.1) then
        DO 521 i=igmx,1,-1
        write(68,67) (eta(i,j),j=jgmx,1,-1)
  521   CONTINUE
        DO 522 i=igmx,1,-1
        write(68,97) (dsss(i,j),j=jgmx,1,-1)
  522   CONTINUE
        DO 523 i=igmx,1,-1
        write(68,699) (float(ibr(i,j)),j=jgmx,1,-1)
  523   CONTINUE
        close(68)
      end if
!
      if(iplane.eq.2) then
        DO 527 i=1,igmx
        read(68,*) (tmp(i,j),j=1,jgmx)
        DO 527 j=1,jgmx
        if(tmp(i,j).gt.eta(i,j)) eta(i,j)=tmp(i,j)
  527   CONTINUE
        DO 528 i=1,igmx
        read(68,*) (tmp(i,j),j=1,jgmx)
        DO 528 j=1,jgmx
        if(tmp(i,j).gt.dsss(i,j)) dsss(i,j)=tmp(i,j)
  528   CONTINUE
        DO 529 i=1,igmx
        read(68,*) (tmp(i,j),j=1,jgmx)
        DO 529 j=1,jgmx
        if(nint(tmp(i,j)).gt.ibr(i,j)) ibr(i,j)=1
  529   CONTINUE
        close(68)
      end if
      end if
!
! *** Very last time check iplane = 1 ***
!
      if(iplane.eq.1) then
!
! *** if iplane = 1 output save location spectral data
!
      IF(iabs(KOUT).GE.1.and.ixmdf.ne.2) THEN
        DO 220 KK=1,iabs(KOUT)
          i=ijsp(1,kk)
          j=ijsp(2,kk)
          if(kout.le.0) go to 238
          if(kdate.gt.0) then
            if(idate.le.9999) then
              idate1=mod(kdate,10)*100000+idate
              kdate1=kdate/10
              write(10,422) kdate1,idate1,i,j,imod*90,h13f(i,j)
            else
              WRITE(10,420) kdate,idate,i,j,imod*90,h13f(i,j)
            end if
          else
          WRITE(10,421) idate,i,j,imod*90,h13f(i,j)
          end if
!
          n00=0
          DO 235 N0=1,NFF
          if(ffcn(n0).lt.fcn(1)) then
!         write(10,510) (0.,mm=1,mddd)
          n00=n00+1
          do mm=1,mddd
             sp1(n00,mm)=0.
          enddo
          go to 235
          end if
!
          if(ffcn(n0).eq.fcn(1)) then
          IF(ICHOICE.LE.1) THEN
!           WRITE(10,510) (SOP(KK,1,MM),MM=1,MDDD)
              n00=n00+1
              do mm=1,mddd
                sp1(n00,mm)=sop(kk,1,mm)
              enddo
          ELSE
            c0=float(mdd+1)/float(md+1)
            DO mm=1,mddd
            m1=int(float(mm-1)/ammd)+1
            m2=m1+1
            if(m2.gt.md) m2=md
            c3=(float(m1)*ammd-float(mm-1))/ammd
            c4=c0*(1.-c3)
            c3=c0*c3
            SPCT(MM)=SOP(KK,1,M1)*c3+SOP(KK,1,M2)*c4
            END DO
!           WRITE(10,510) (SPCT(MM),MM=1,MDDD)
              n00=n00+1
              do mm=1,mddd
                sp1(n00,mm)=spct(mm)
              enddo
          END IF
          go to 235
          end if
!
          if(ffcn(n0).gt.fcn(nf)) then
          c1=(ffcn(nff)-ffcn(n0))/(ffcn(nff)-fcn(nf))
          IF(ICHOICE.LE.1) THEN
!           WRITE(10,510) (SOP(KK,NF,MM)*c1,MM=1,MDDD)
              n00=n00+1
              do mm=1,mddd
                sp1(n00,mm)=sop(kk,nf,mm)*c1
              enddo
          ELSE
            c0=float(mdd+1)/float(md+1)
            DO mm=1,mddd
            m1=int(float(mm-1)/ammd)+1
            m2=m1+1
            if(m2.gt.md) m2=md
            c3=(float(m1)*ammd-float(mm-1))/ammd
            c4=c0*(1.-c3)
            c3=c0*c3
            SPCT(MM)=SOP(KK,NF,M1)*c3+SOP(KK,NF,M2)*c4
            END DO
!           WRITE(10,510) (SPCT(MM)*c1,MM=1,MDDD)
              n00=n00+1
              do mm=1,mddd
                sp1(n00,mm)=spct(mm)*c1
              enddo
          END IF
          go to 235
          end if
!
            DO 230 NN=1,NF-1
            nn1=nn+1
            if(ffcn(n0).gt.fcn(nn).and.ffcn(n0).le.fcn(nn1)) then
            c0=fcn(nn1)-fcn(nn)
            c2=(ffcn(n0)-fcn(nn))/c0
            c1=(fcn(nn1)-ffcn(n0))/c0
            IF(ICHOICE.LE.1) THEN
!           WRITE(10,510) (SOP(KK,NN,MM)*c1+SOP(KK,NN1,MM)*c2,MM=1,MDDD)
              n00=n00+1
              do mm=1,mddd
                sp1(n00,mm)=sop(kk,nn,mm)*c1+sop(kk,nn1,mm)*c2
              enddo
            ELSE
            c0=float(mdd+1)/float(md+1)
            DO mm=1,mddd
            m1=int(float(mm-1)/ammd)+1
            m2=m1+1
            if(m2.gt.md) m2=md
            c3=(float(m1)*ammd-float(mm-1))/ammd
            c4=c0*(1.-c3)
            c3=c0*c3
            SPCT(MM)=SOP(KK,NN,M1)*c3+SOP(KK,NN,M2)*c4
            SPCT1(MM)=SOP(KK,NN1,M1)*c3+SOP(KK,NN1,M2)*c4
            END DO
!           WRITE(10,510) (SPCT(MM)*c1+SPCT1(MM)*c2,MM=1,MDDD)
              n00=n00+1
              do mm=1,mddd
                sp1(n00,mm)=spct(mm)*c1+spct1(mm)*c2
              enddo
            END IF
            go to 235
            end if
  230       CONTINUE
!
  235     CONTINUE
!
          sum=0.
          sum1=0.
          sum2=0.0000001
          do mm=1,mddd
          sum1=sum1+sp1(nff,mm)
          enddo
!
          do nn=nff,2,-1
            sum0=0.
            do mm=1,mddd
            sum0=sum0+sp1(nn-1,mm)
            enddo
            if(sum0.eq.0..and.sum1.gt.0.) then
              do mm=1,mddd
              sp1(nn-1,mm)=sp1(nn,mm)*.3
              sp1(nn,mm)=sp1(nn,mm)*.7
              enddo
              sum1=sum1*.7
              sum0=sum1*.3
            endif
            cc=sum1*dfinp(nn)
            sum=sum+cc
            sum1=sum0
            sum2=sum2+cc*ffcn(nn)
          end do
            cc=sum1*dfinp(1)
            sum=sum+cc
            sum2=sum2+cc*ffcn(1)
          t13fij=sum/sum2
          h13fij=4.*sqrt(sum*pai/float(1+mddd))
          PL1E(kk)=sum
          PL1T(kk)=t13fij

!      write(*,*) 'ok12,h13f,h13fij,t13fij,nff,mddd=',h13f(i,j),h13fij,t13fij,nff,mddd
!      write(*,*) 'dfinp(1:5)=',dfinp(1),dfinp(2),dfinp(3),dfinp(4),dfinp(5)

          c1=1.
          if(h13fij.gt..001) then
          c1=(h13f(i,j)/h13fij)**2
          end if
!
          do nn=1,nff
            WRITE(10,510) (sp1(nn,mm)*c1,mm=1,mddd)
          enddo
!
  238     continue
!
  220   CONTINUE
      END IF
!
!
        if(igetfile20.eq.1) then
!
        read(39,*) tship1,((xship1(n),yship1(n)),n=1,nship)
        read(39,*) tship2,((xship2(n),yship2(n)),n=1,nship)
        backspace(39)
!
!         write(*,*) tship1,((xship1(n),yship1(n)),n=1,nship)
          write(*,*) tship2,((xship2(n),yship2(n)),n=1,nship)
!
          do 299 n=1,nship

          if(xship1(n).lt.0.001.or.xship2(n).lt..001) go to 299
          xship0= (xship2(n)-x0)*cosaz+(yship2(n)-y0)*sinaz
          yship0=-(xship2(n)-x0)*sinaz+(yship2(n)-y0)*cosaz
!
          nix=1
        do i=1,igmx
          if(xship0.le.disx(i)) then
            nix=i
            exit
          end if
        end do
!
          niy=1
        do j=1,jgmx
          if(yship0.le.disy(j)) then
            niy=j
            exit
          end if
        end do
!
          c4=yship1(n)-yship2(n)
          c3=xship1(n)-xship2(n)
!
          aship0=atan2(c4,c3)-azimuth*rad
          vship0=sqrt(c3**2+c4**2)/(tship2-tship1)/3600.
          tship0=pai*vship0/g*(1.+tanh(5./vship0))
!         hship00=min(g*(tship0/12.)**2,sqrt(1.125/tship0))
!
          shipR=sqrt(ShipB(n)*ShipD(n)*g/ShipL(n))
          ship20LD=20.*(shipL(n)+shipD(n))
!
          cb=0.2+0.5*tanh(ShipR/vship0)
          alpha=2.35*(1.-cb)
          dship0=d1(nix,niy)
!
          Fstar=vship0/sqrt(g*ShipL(n))*exp(alpha*ShipD(n)/dship0)
          beta=1+8.*tanh(0.45*(ShipL(n)/ShipS(n)-2.))**3
      cship0=vship0*ShipL(n)*ShipD(n)/dship0/g/tship0/5.
      hship0=cship0*beta*(Fstar-0.1)**2*(2.*ShipL(n)/ShipB(n))**0.333
!
!         write(*,*) n,xship0,yship0,aship0,dship0,hship0,hship00
! *** xship0, yship0, aship0 are local disx, disy, and tail angle
!
          nix=7000
          do nix=1,igmx
      if(xship0.ge.disx(nix).and.xship0.lt.disx(nix+1)) exit
          end do
!
          niy=7000
          do niy=1,jgmx
      if(yship0.ge.disy(niy).and.yship0.lt.disy(niy+1)) exit
          end do
      if(nix.gt.igmx.or.niy.gt.jgmx) go to 299
!
          xship0=disx(nix)+dvarxx(nix)/2.
          yship0=disy(niy)+dvaryy(niy)/2.
!
!     write(*,*) 'nix,niy,h13=',nix,niy,h13f(nix,niy),hship0,tship0,vship0
          c3=h13f(nix,niy)*cos(dmnf(nix,niy))+hship0*cos(aship0+pai)
          c4=h13f(nix,niy)*sin(dmnf(nix,niy))+hship0*sin(aship0+pai)
          dmnf(nix,niy)=atan2(c4,c3)
      t13f(nix,niy)=(t13f(nix,niy)*h13f(nix,niy)+tship0*hship0)/(h13f(nix,niy)+hship0)
          h13f(nix,niy)=h13f(nix,niy)+hship0
!
      Eship=0.
      Aship=0.
!
          do 391 kni=-200,200
          do 391 knj=-200,200
!
      if(kni.eq.0.and.knj.eq.0) go to 391
          kii=kni+nix
          kjj=knj+niy
      if(kii.gt.igmx) go to 391
      if(kjj.gt.jgmx) go to 391
      if(kii.le.0)    go to 391
      if(kjj.le.0)    go to 391
!     write(*,*) 'kii,kjj,d1=',kii,kjj,d1(kii,kjj)
      if(d1(kii,kjj).lt.0.1) go to 391
!
          xshipc=(disx(kii)+dvarxx(kii)/2.)-xship0
          yshipc=(disy(kjj)+dvaryy(kjj)/2.)-yship0
          zshipc=sqrt(xshipc**2+yshipc**2)
!
      if(zshipc.gt.1000.) go to 391
!
          ashipc=atan2(yshipc,xshipc)-aship0
          if(ashipc.gt.pai)  ashipc=ashipc-pai2
          if(ashipc.gt.pai)  ashipc=ashipc-pai2
          if(ashipc.lt.-pai) ashipc=ashipc+pai2
          if(ashipc.lt.-pai) ashipc=ashipc+pai2
!
          inear=0
          if(zshipc.le.50.) inear=1
!
          if(inear.eq.1.and.abs(ashipc).ge.hpai/.8) go to 391
          if(inear.eq.0.and.abs(ashipc).gt.hpai/2.) go to 391
!
          c2=hship0*abs(sin(4.0*ashipc))*ship20LD/(ship20LD+zshipc)
          c2=c2*tanh(10.-zshipc/100.)
          c5=aship0-2.*ashipc+pai
          kni201=kni+201
          knj201=knj+201
          Eship(kni201,knj201)=c2**2
          Aship(kni201,knj201)=c5
          c3=h13f(kii,kjj)*cos(dmnf(kii,kjj))+c2*cos(c5)
          c4=h13f(kii,kjj)*sin(dmnf(kii,kjj))+c2*sin(c5)
          dmnf(kii,kjj)=atan2(c4,c3)
      t13f(kii,kjj)=(t13f(kii,kjj)*h13f(kii,kjj)+tship0*c2)/(h13f(kii,kjj)+c2)
          h13f(kii,kjj)=h13f(kii,kjj)+c2
!
  391     continue
!
          Eship=Eship*g/32.
!
          do 392 kni201=2,400
          do 392 knj201=2,400
!
          kii=kni201+nix-201
          kjj=knj201+niy-201
!
          if(kii.ge.igmx) go to 392
          if(kjj.ge.jgmx) go to 392
          if(kii.le.1)    go to 392
          if(kjj.le.1)    go to 392
        if(d1(kii,kjj).lt.0.1) go to 392
!
          dx2=dvarxx(kii)+(dvarxx(kii-1)+dvarxx(kii+1))/2.
          c2=Eship(kni201+1,knj201)*cos(Aship(kni201+1,knj201))**2
          c3=Eship(kni201-1,knj201)*cos(Aship(kni201-1,knj201))**2
          dy2=dvaryy(kjj)+(dvaryy(kjj-1)+dvaryy(kjj+1))/2.
          c4=Eship(kni201,knj201+1)*sin(Aship(kni201,knj201+1)*2.)
          c5=Eship(kni201,knj201-1)*sin(Aship(kni201,knj201-1)*2.)
          wxrs(kii,kjj)=wxrs(kii,kjj)+(c3-c2)/dx2+(c5-c4)/dy2
          c2=Eship(kni201,knj201+1)*sin(Aship(kni201,knj201+1))**2
          c3=Eship(kni201,knj201-1)*sin(Aship(kni201,knj201-1))**2
          c4=Eship(kni201+1,knj201)*sin(Aship(kni201+1,knj201)*2.)
          c5=Eship(kni201-1,knj201)*sin(Aship(kni201-1,knj201)*2.)
          wyrs(kii,kjj)=wyrs(kii,kjj)+(c3-c2)/dy2+(c5-c4)/dx2
!
  392     continue
!
  299     continue
!
        end if
!
!***** There is need to reset back reflection variables here
!
        h13r = 0.
        t13r = t13min
        dmnr = 0.
        RETURN
!
      endif
!
! *** Return here if iplane = 1 ***
!
      if(imod.eq.0) then
      DO 152 J=1,JGMX
      if(iview.eq.0) then
      do i=1,igmx
      if(dmn(i,j).lt.-99.9) dmn(i,j)=-99.9 
      if(dw13(i,j).lt.-99.9) dw13(i,j)=-99.9
      if(da13(i,j).lt.-99.9) da13(i,j)=-99.9
      end do
      end if
152   CONTINUE
      end if
!
      if(imod.eq.2) then
      DO 252 J=1,JGMX
      DO 252 I=1,IGMX
      DMN(I,J)=DMN(I,J)+180.
      if(iwet.eq.-1.or.iwet.eq.-3) then
      dw13(i,j)=dw13(i,j)+180.
      da13(i,j)=da13(i,j)+180.
      end if
      if(iview.eq.0) then
      if(dmn(i,j).lt.-99.9) dmn(i,j)=-99.9
      if(dw13(i,j).lt.-99.9) dw13(i,j)=-99.9
      if(da13(i,j).lt.-99.9) da13(i,j)=-99.9
      end if
252   CONTINUE
      end if
!
      if(imod.eq.1.or.imod.eq.3) then
      cc=90.*float(imod)
      DO 253 J=1,JGMX
      DO 253 I=1,IGMX
      DMN(I,J)=DMN(I,J)+cc
      if(iwet.eq.-1.or.iwet.eq.-3) then
      dw13(i,j)=dw13(i,j)+cc
      da13(i,j)=da13(i,j)+cc
      end if
      if(iview.eq.0) then
      if(dmn(i,j).lt.-99.9) dmn(i,j)=-99.9
      if(dw13(i,j).lt.-99.9) dw13(i,j)=-99.9
      if(da13(i,j).lt.-99.9) da13(i,j)=-99.9
      end if
253   CONTINUE
      end if
!
      if(ixmdf.ge.1) go to 531
! ***
      if (ibreak.eq.1) then
      if(kdate.gt.0) then
        if(idate.le.9999) then
          idate1=mod(kdate,10)*100000+idate
          kdate1=kdate/10
          write(17,'(i7,i6)') kdate1,idate1
        else
          write(17,'(i8,i5)') kdate,idate
        end if
      else
      write(17,'(i8)') idate
      end if
        if(imod.eq.0) then
          do j=jgmx,1,-1
            write(17,'(16i5)') (ibr(i,j),i=1,igmx)
          end do
        end if
        if(imod.eq.1) then
          do i=igmx,1,-1
            write(17,'(16i5)') (ibr(i,j),j=jgmx,1,-1)
          end do
        end if
        if(imod.eq.2) then
          do j=1,jgmx
            write(17,'(16i5)') (ibr(i,j),i=igmx,1,-1)
          end do
        end if
        if(imod.eq.3) then
          do i=1,igmx
            write(17,'(16i5)') (ibr(i,j),j=1,jgmx)
          end do
        end if
      end if
!
      if(kdate.gt.0) then
        if(idate.le.9999) then
          idate1=mod(kdate,10)*100000+idate
          kdate1=kdate/10
          write(66,'(i7,i6)') kdate1,idate1
        else
          write(66,'(i8,i5)') kdate,idate
        end if
      else
      write(66,'(i8)') idate
      end if
      if(iwet.lt.0) then
      if(isteer.eq.1.and.idate.lt.10000) then
        if(iwet.ne.-1) write(95,*) iidate
        if(iwet.ne.-2) then
        write(96,*) iidate
        write(97,*) iidate
        end if
      else
       if(kdate.eq.0) then
        if(iwet.ne.-1) write(95,'(1x,i8)') idate
        if(iwet.ne.-2) then
        write(96,'(1x,i8)') idate
        write(97,'(1x,i8)') idate
        end if
       else
        if(iwet.ne.-1) then
          if(idate.le.9999) then
            idate1=mod(kdate,10)*100000+idate
            kdate1=kdate/10
            write(95,'(i7,i6)') kdate1,idate1
          else
            write(95,'(i8,i5)') kdate,idate
          end if
        end if
        if(iwet.ne.-2) then
          if(idate.le.9999) then
            idate1=mod(kdate,10)*100000+idate
            kdate1=kdate/10
            write(96,'(i7,i6)') kdate1,idate1
            write(97,'(i7,i6)') kdate1,idate1
          else
            write(96,'(i8,i5)') kdate,idate
            write(97,'(i8,i5)') kdate,idate
          end if
        end if
       end if
      end if
      end if
!

!     write(*,*) 'imod =', imod

      if(imod.eq.0) then
       DO 50 J=JGMX,1,-1
        if(depmax.ge.1..and.irs.le.1) then
         write(66,70) (H13S(I,J),I=1,IGMX)
        else
         write(66,67) (H13S(I,J),I=1,IGMX)
        end if
         if(iwet.lt.0) then
         if(iwet.ne.-1) write(95,67) (H13S(i,j),i=1,igmx)
         if(iwet.ne.-2) then
         write(96,67)(sw13(i,j),i=1,igmx)
         write(97,67)(sa13(i,j),i=1,igmx)
         end if
         end if
      if(irs.ge.2.and.kout.ge.0) write(98,67) (eta(i,j),i=1,igmx)
50    CONTINUE
67    format(16f7.3)
70    format(16f6.2)
      DO 51 J=JGMX,1,-1
        if(depmax.ge.1.) then
         write(66,68) (T13(I,J),I=1,IGMX)
        else
         write(66,67) (T13(I,J),I=1,IGMX)
        end if
         if(iwet.lt.0) then
         if(iwet.ne.-1) write(95,68) (T13(i,j),i=1,igmx)
         if(iwet.ne.-2) then
         write(96,68)(tw13(i,j),i=1,igmx)
         write(97,68)(ta13(i,j),i=1,igmx)
         end if
         end if
      if(irs.ge.2.and.kout.ge.0) write(98,97) (dsss(i,j),i=1,igmx)
97    format(16f9.4)
68    format(16f5.1)
51    CONTINUE
      DO 52 J=JGMX,1,-1
      write(66,699)  (DMN(I,J),I=1,IGMX)
        if(iwet.lt.0) then
        if(iwet.ne.-1) write(95,69) (dmn(i,j),i=1,igmx)
        if(iwet.ne.-2) then
        write(96,69)(dw13(i,j),i=1,igmx)
        write(97,69)(da13(i,j),i=1,igmx)
        end if
        end if
      if(irs.ge.2.and.kout.ge.0)   &
      write(98,699) (dmn(i,j)+sign(90.,eta(i,j))-90.,i=1,igmx)
69    format(16f6.1)
699   format(16f6.0)
52    CONTINUE
      end if
!
      if(imod.eq.2) then
       DO 2050 J=1,JGMX
        if(depmax.ge.1..and.irs.le.1) then
         write(66,70)  (H13S(I,J),I=IGMX,1,-1)
        else
         write(66,67) (H13S(I,J),I=IGMX,1,-1)
        end if
         if(iwet.lt.0) then
         if(iwet.ne.-1) write(95,67) (H13S(i,j),i=igmx,1,-1)
         if(iwet.ne.-2) then
         write(96,67)(sw13(i,j),i=igmx,1,-1)
         write(97,67)(sa13(i,j),i=igmx,1,-1)
         end if
         end if
        if(irs.ge.2.and.kout.ge.0) write(98,67)(eta(i,j),i=igmx,1,-1)
2050  CONTINUE
       DO 2051 J=1,JGMX
        if(depmax.ge.1.) then
         write(66,68) (T13(I,J),I=IGMX,1,-1)
        else
         write(66,67) (T13(I,J),I=IGMX,1,-1)
        end if
         if(iwet.lt.0) then
         if(iwet.ne.-1) write(95,68) (T13(i,j),i=igmx,1,-1)
         if(iwet.ne.-2) then
         write(96,68)(tw13(i,j),i=igmx,1,-1)
         write(97,68)(ta13(i,j),i=igmx,1,-1)
         end if
         end if
        if(irs.ge.2.and.kout.ge.0) write(98,97)(dsss(i,j),i=igmx,1,-1)
2051  CONTINUE
      DO 2052 J=1,JGMX
      write(66,69)  (DMN(I,J),I=IGMX,1,-1)
        if(iwet.lt.0) then
        if(iwet.ne.-1) write(95,69) (dmn(i,j),i=igmx,1,-1)
        if(iwet.ne.-2) then
        write(96,69)(dw13(i,j),i=igmx,1,-1)
        write(97,69)(da13(i,j),i=igmx,1,-1)
        end if
        end if
       if(irs.ge.2.and.kout.ge.0)  &
       write(98,699) (dmn(i,j)+sign(90.,eta(i,j))-90.,i=igmx,1,-1)
2052   CONTINUE
      end if
!
      if(imod.eq.1) then
      DO 1050 I=IGMX,1,-1
      if(depmax.ge.1..and.irs.le.1) then
      write(66,70)  (H13S(I,J),J=JGMX,1,-1)
      else
      write(66,67) (H13S(I,J),J=JGMX,1,-1)
      end if
        if(iwet.lt.0) then
        if(iwet.ne.-1) write(95,67) (H13S(i,j),j=jgmx,1,-1)
        if(iwet.ne.-2) then
        write(96,67)(sw13(i,j),j=jgmx,1,-1)
        write(97,67)(sa13(i,j),j=jgmx,1,-1)
        end if
        end if
      if(irs.ge.2.and.kout.ge.0) write(98,67)(eta(i,j),j=jgmx,1,-1)
1050  CONTINUE
      DO 1051 I=IGMX,1,-1
      if(depmax.ge.1.) then
      write(66,68) (T13(I,J),J=JGMX,1,-1)
      else
      write(66,67) (T13(I,J),J=JGMX,1,-1)
      end if
        if(iwet.lt.0) then
        if(iwet.ne.-1) write(95,68) (T13(i,j),j=jgmx,1,-1)
        if(iwet.ne.-2) then
        write(96,68)(tw13(i,j),j=jgmx,1,-1)
        write(97,68)(ta13(i,j),j=jgmx,1,-1)
        end if
        end if
      if(irs.ge.2.and.kout.ge.0) write(98,97)(dsss(i,j),j=jgmx,1,-1)
1051  CONTINUE
      DO 1052 I=IGMX,1,-1
      write(66,69)  (DMN(I,J),J=JGMX,1,-1)
        if(iwet.lt.0) then
        if(iwet.ne.-1) write(95,69) (dmn(i,j),j=jgmx,1,-1)
        if(iwet.ne.-2) then
        write(96,69)(dw13(i,j),j=jgmx,1,-1)
        write(97,69)(da13(i,j),j=jgmx,1,-1)
        end if
        end if
      if(irs.ge.2.and.kout.ge.0)  &
      write(98,699) (dmn(i,j)+sign(90.,eta(i,j))-90.,j=jgmx,1,-1)
1052  CONTINUE
      end if
!
      if(imod.eq.3) then
      DO 3050 I=1,IGMX
      if(depmax.ge.1..and.irs.le.1) then
      write(66,70)  (H13S(I,J),J=1,JGMX)
      else
      write(66,67) (H13S(I,J),J=1,JGMX)
      end if
        if(iwet.lt.0) then
        if(iwet.ne.-1) write(95,67) (H13S(i,j),j=1,jgmx)
        if(iwet.ne.-2) then
        write(96,67)(sw13(i,j),j=1,jgmx)
        write(97,67)(sa13(i,j),j=1,jgmx)
        end if
        end if
      if(irs.ge.2.and.kout.ge.0) write(98,67)(eta(i,j),j=1,jgmx)
3050  CONTINUE
      DO 3051 I=1,IGMX
      if(depmax.ge.1.) then
      write(66,68) (T13(I,J),J=1,JGMX)
      else
      write(66,67) (T13(I,J),J=1,JGMX)
      end if
        if(iwet.lt.0) then
        if(iwet.ne.-1) write(95,68) (T13(i,j),j=1,jgmx)
        if(iwet.ne.-2) then
        write(96,68)(tw13(i,j),j=1,jgmx)
        write(97,68)(ta13(i,j),j=1,jgmx)
        end if
        end if
      if(irs.ge.2.and.kout.ge.0) write(98,97)(dsss(i,j),j=1,jgmx)
3051  CONTINUE
      DO 3052 I=1,IGMX
      write(66,69)  (DMN(I,J),J=1,JGMX)
        if(iwet.lt.0) then
        if(iwet.ne.-1) write(95,69) (DMN(i,j),j=1,jgmx)
        if(iwet.ne.-2) then
        write(96,69)(dw13(i,j),j=1,jgmx)
        write(97,69)(da13(i,j),j=1,jgmx)
        end if
        end if
      if(irs.ge.2.and.kout.ge.0)  &
      write(98,699) (dmn(i,j)+sign(90.,eta(i,j))-90.,j=1,jgmx)
3052  CONTINUE
      end if
! ***
531   CONTINUE
! ***
!
      if (kout.ge.0) then
      if (irs.ge.1.and.ixmdf.eq.0) then
      if(kdate.gt.0) then
        if(idate.le.9999) then
          idate1=mod(kdate,10)*100000+idate
          kdate1=kdate/10
          write(18,'(i7,i6)') kdate1,idate1
        else
          write(18,'(i8,i5)') kdate,idate
        end if
      else
      write(18,'(i8)') idate
      end if
          if(imod.eq.0) then
          do j=jgmx,1,-1
             write(18,*) (wxrs(i,j),wyrs(i,j),i=1,igmx)
          enddo
          end if
          if(imod.eq.1) then
          do i=igmx,1,-1
             write(18,*) (-wyrs(i,j),wxrs(i,j),j=jgmx,1,-1)
          enddo
          end if
          if(imod.eq.2) then
          do j=1,jgmx
             write(18,*) (-wxrs(i,j),-wyrs(i,j),i=igmx,1,-1)
          enddo
          end if
          if(imod.eq.3) then
          do i=1,igmx
             write(18,*) (wyrs(i,j),-wxrs(i,j),j=1,jgmx)
          enddo
          end if
      end if
      end if
!
      if(iabs(kout).GE.1.and.ixmdf.ne.2) then
!
        do 120 kk=1,iabs(kout)
          i=ijsp(1,kk)
          j=ijsp(2,kk)
!
          if(imod.eq.0) then
          ii=i
          jj=j
          end if
          if(imod.eq.1) then
          ii=j
          jj=jgmx-i+1
          end if
          if(imod.eq.2) then
          ii=igmx-i+1
          jj=jgmx-j+1
          end if
          if(imod.eq.3) then
          ii=igmx-j+1
          jj=i
          end if
!
          if(kout.le.0) go to 138
!
          h13cij=sqrt(abs(h13s(ii,jj)**2-h13f(i,j)**2))
!
          if(kdate.gt.0) then
            if(idate.le.9999) then
              idate1=mod(kdate,10)*100000+idate
              kdate1=kdate/10
              write(10,422) kdate1,idate1,i,j,imod*90,h13cij
  422     FORMAT(i7,i6,3i5,f8.3)
            else
              WRITE(10,420) kdate,idate,i,j,imod*90,h13cij
  420     FORMAT(i8,i5,3i5,f8.3)
            end if
          else
          WRITE(10,421) idate,i,j,imod*90,h13cij
  421     FORMAT(5x,i8,3i5,f8.3)
          end if
!
          n00=0
          DO 135 N0=1,NFF
          if(ffcn(n0).lt.fcn(1)) then
!         write(10,510) (0.,mm=1,mddd)
          n00=n00+1
          do mm=1,mddd
             sp1(n00,mm)=0.
          enddo
          go to 135
          end if
!
          if(ffcn(n0).eq.fcn(1)) then
          IF(ICHOICE.LE.1) THEN
!           WRITE(10,510) (SOP(KK,1,MM),MM=1,MDDD)
              n00=n00+1
              do mm=1,mddd
                sp1(n00,mm)=sop(kk,1,mm)
              enddo
          ELSE
            c0=float(mdd+1)/float(md+1)
            DO mm=1,mddd
            m1=int(float(mm-1)/ammd)+1
            m2=m1+1
            if(m2.gt.md) m2=md
            c3=(float(m1)*ammd-float(mm-1))/ammd
            c4=c0*(1.-c3)
            c3=c0*c3
            SPCT(MM)=SOP(KK,1,M1)*c3+SOP(KK,1,M2)*c4
            END DO
!           WRITE(10,510) (SPCT(MM),MM=1,MDDD)
              n00=n00+1
              do mm=1,mddd
                sp1(n00,mm)=spct(mm)
              enddo
          END IF
          go to 135
          end if
!
          if(ffcn(n0).gt.fcn(nf)) then
          c1=(ffcn(nff)-ffcn(n0))/(ffcn(nff)-fcn(nf))
          IF(ICHOICE.LE.1) THEN
!           WRITE(10,510) (SOP(KK,NF,MM)*c1,MM=1,MDDD)
              n00=n00+1
              do mm=1,mddd
                sp1(n00,mm)=sop(kk,nf,mm)*c1
              enddo
          ELSE
            c0=float(mdd+1)/float(md+1)
            DO mm=1,mddd
            m1=int(float(mm-1)/ammd)+1
            m2=m1+1
            if(m2.gt.md) m2=md
            c3=(float(m1)*ammd-float(mm-1))/ammd
            c4=c0*(1.-c3)
            c3=c0*c3
            SPCT(MM)=SOP(KK,NF,M1)*c3+SOP(KK,NF,M2)*c4
            END DO
!           WRITE(10,510) (SPCT(MM)*c1,MM=1,MDDD)
              n00=n00+1
              do mm=1,mddd
                sp1(n00,mm)=spct(mm)*c1
              enddo
          END IF
          go to 135
          end if
!
            DO 130 NN=1,NF-1
            nn1=nn+1
            if(ffcn(n0).gt.fcn(nn).and.ffcn(n0).le.fcn(nn1)) then
            c0=fcn(nn1)-fcn(nn)
            c2=(ffcn(n0)-fcn(nn))/c0
            c1=(fcn(nn1)-ffcn(n0))/c0
            IF(ICHOICE.LE.1) THEN
!           WRITE(10,510) (SOP(KK,NN,MM)*c1+SOP(KK,NN1,MM)*c2,MM=1,MDDD)
              n00=n00+1
              do mm=1,mddd
                sp1(n00,mm)=sop(kk,nn,mm)*c1+sop(kk,nn1,mm)*c2
              enddo
            ELSE
            c0=float(mdd+1)/float(md+1)
            DO mm=1,mddd
            m1=int(float(mm-1)/ammd)+1
            m2=m1+1
            if(m2.gt.md) m2=md
            c3=(float(m1)*ammd-float(mm-1))/ammd
            c4=c0*(1.-c3)
            c3=c0*c3
            SPCT(MM)=SOP(KK,NN,M1)*c3+SOP(KK,NN,M2)*c4
            SPCT1(MM)=SOP(KK,NN1,M1)*c3+SOP(KK,NN1,M2)*c4
            END DO
!           WRITE(10,510) (SPCT(MM)*c1+SPCT1(MM)*c2,MM=1,MDDD)
              n00=n00+1
              do mm=1,mddd
                sp1(n00,mm)=spct(mm)*c1+spct1(mm)*c2
              enddo
            END IF
            go to 135
            end if
  130       CONTINUE
!
  135     CONTINUE
!
          sum=0.
          sum1=0.
          sum2=0.0000001
          do mm=1,mddd
          sum1=sum1+sp1(nff,mm)
          enddo
!
          do nn=nff,2,-1
            sum0=0.
            do mm=1,mddd
            sum0=sum0+sp1(nn-1,mm)
            enddo
            if(sum0.eq.0..and.sum1.gt.0.) then
              do mm=1,mddd
              sp1(nn-1,mm)=sp1(nn,mm)*.3
              sp1(nn,mm)=sp1(nn,mm)*.7
              enddo
              sum1=sum1*.7
              sum0=sum1*.3
            endif
            cc=sum1*dfinp(nn)
            sum=sum+cc
            sum1=sum0
            sum2=sum2+cc*ffcn(nn)
          end do
            cc=sum1*dfinp(1)
            sum=sum+cc
            sum2=sum2+cc*ffcn(1)
          t13fij=sum/sum2
          h13fij=4.*sqrt(sum*pai/float(1+mddd))

!      write(*,*) 'ok11,h13cij,h13fij,t13fij=',h13cij,h13fij,t13fij
!      write(*,*) 'dfinp(1:5)=',dfinp(1),dfinp(2),dfinp(3)
!
          if(iview.eq.1.and.igrav.eq.3) then
          t13(ii,jj)=(PL1T(kk)*PL1E(kk)+sum*t13fij)/(PL1E(kk)+sum)
!      write(*,*) 'ok10,h13cij,h13fij,t13fij=',h13cij,h13fij,t13(ii,jj)
          end if
!
          c1=1.
          if(h13fij.gt..001) then
          c1=(h13cij/h13fij)**2
          end if
!
          do nn=1,nff
            WRITE(10,510) (sp1(nn,mm)*c1,mm=1,mddd)
          enddo
!
! 510 FORMAT(8F10.4)
  510       format(17(1x,f8.3))
  138     continue
!
          ll=0
          if(icur.eq.0) then
          ll=1
          if(ijstruc2.ge.1) then
          do k=1,ijstruc2
          if(istruc2(k).eq.i.and.jstruc2(k).eq.j) then
          ll=0
          go to 601
          end if
          end do
          end if
          if(ijstruc4.ge.1) then
          do k=1,ijstruc4
          if(istruc4(k).eq.i.and.jstruc4(k).eq.j) then
          ll=0
          go to 601
          end if
          end do
          end if
  601     continue
          end if
!
          cc1=0.
!         cc3=0.
          if(d1(ii,jj).gt.0.001.and.ll.eq.0) then
          cc2=0.
          i14=ii-14
          if(i14.lt.1) i14=1
          do ij=i14,ii
          if(d1(ij,jj).lt..01) go to 901
          if(icur.ge.1) then
          cc1=cc1+sqrt(u1(ij,jj)**2+v1(ij,jj)**2)*d1(ij,jj)
          else
          if(ibr(ij,jj).eq.0) then
          cc1=cc1+0.18*g**0.55*h13s(ij,jj)**1.45*t13(ij,jj)**0.1
!         cc3=cc3+0.3*sqrt(g)*h13s(ij,jj)**1.5
          else
          cc1=cc1+2.2*g**0.55*h13s(ij,jj)**1.45*t13(ij,jj)**0.1
!         cc3=cc3+2.*sqrt(g)*h13s(ij,jj)**1.5
          end if
          end if
          cc2=cc2+1.
  901     continue
          end do
          cc1=cc1/cc2
!         cc3=cc3/cc2
          end if

            ccc=u1(ii,jj)
            vv1=v1(ii,jj)
!
           cosaz0=cos(azimuth*rad)
           sinaz0=sin(azimuth*rad)
!
           uu1=ccc*cosaz0-vv1*sinaz0
           vv1=ccc*sinaz0+vv1*cosaz0
!
           if(igrav.ge.1) then
            uu1=uu1/1.2
            vv1=vv1/1.2
           end if
!
          if(kdate.gt.0) then
            if(idate.le.9999) then
            idate1=mod(kdate,10)*100000+idate
            kdate1=kdate/10
          write(12,907)kdate1,idate1,i,j,h13s(ii,jj),t13(ii,jj),dmn(ii,jj), &
          sw13(ii,jj),tw13(ii,jj),dw13(ii,jj), &
          sa13(ii,jj),ta13(ii,jj),da13(ii,jj), &
          ibr(ii,jj),dsss(ii,jj),cc1,uu1,vv1
  907     format(i7,i6, 2(1x,i4),f7.3, f6.2, f7.1, &
          2(1x,f7.2, f6.2, f7.1),i6, f7.3,5f7.2)
            else
          write(12,905)kdate,idate,i,j,h13s(ii,jj),t13(ii,jj),dmn(ii,jj), &
          sw13(ii,jj),tw13(ii,jj),dw13(ii,jj), &
          sa13(ii,jj),ta13(ii,jj),da13(ii,jj), &
          ibr(ii,jj),dsss(ii,jj),cc1,uu1,vv1
  905     format(i8,i5, 2(1x,i4),f7.3, f6.2, f7.1, &
          2(1x,f7.2, f6.2, f7.1),i6, f7.3,5f7.2)
            end if
          else
          write(12,906)idate,i,j,h13s(ii,jj),t13(ii,jj),dmn(ii,jj), &
          sw13(ii,jj),tw13(ii,jj),dw13(ii,jj), &
          sa13(ii,jj),ta13(ii,jj),da13(ii,jj), &
          ibr(ii,jj),dsss(ii,jj),cc1,uu1,vv1
  906     format(5x,i8, 2(1x,i4),f7.3, f6.2, f7.1, &
          2(1x,f7.2, f6.2, f7.1),i6, f7.3,5f7.2)
          end if
!
  120   CONTINUE
      end if

!
      DEALLOCATE (SPCT,SPCT1,dep0,dsss,gxx1,gyy1,tmp)
      RETURN
      END SUBROUTINE
!
!********************************************************************************
!      LOGICAL FUNCTION ExistFile_inline (FileName)
!---- ExistFile -------------------------------------------------------S
!  PURPOSE - Determines whether or not a file exists.
!  ARGUMENT DEFINITIONS
!     FileName - Name of the file.
!  ARGUMENT DECLARATIONS
!      CHARACTER FileName*(*)
!----------------------------------------------------------------------F
!
!      INQUIRE (FILE = FileName, EXIST = ExistFile_inline)
!      RETURN
!      END FUNCTION
!
!********************************************************************************
      SUBROUTINE STWfiles_inline (SimFile)
!---STWfiles---------------------------------------------------------S
! PURPOSE Reads the names of the STWAVE data files from the sim file
! ARGUMENT DEFINITIONS
!     SimFile  - Name of the simulation file.
! ARGUMENT DECLARATIONS
      CHARACTER     SimFile*180
      CHARACTER     Name*180
!----------------------------------------------------------------------F
! VARIABLE DECLARATIONS
      LOGICAL       ExistFile
      CHARACTER     FileName*180
      CHARACTER     ThePath*180
      CHARACTER     CARD*6
!     Input file variables
      CHARACTER*180  OptsFile, DepFile, CurrFile, EngInFile, StrucFile
      CHARACTER*10  text1
      CHARACTER*80  text2
!     Output/Input file variable
      CHARACTER *180 NestFile, SurgeFile
!     Output file variables
      CHARACTER*180 WaveFile,ObsFile,EngOutFile,BreakFile,RadsFile 
      CHARACTER*180 MudFile, FricFile, FrflFile, BrflFile, WindFile
      CHARACTER*180 SpecFile, XMDFFile,SetupFile,ShipFile
      CHARACTER*180 SeaFile, SwellFile, TotalFile                      !Mitch 3/22/2017
      common /FileNames/ OptsFile, DepFile, CurrFile, EngInFile,    &
                         WaveFile, ObsFile, EngOutFile, NestFile,   &
                         BreakFile, RadsFile, StrucFile, SurgeFile, &
                         MudFile, FricFile, FrflFile, BrflFile,     &
                         SpecFile, WindFile, XMDFFile, SetupFile,   &  !Mitch 3/22/2017
                         SeaFile, SwellFile, ShipFile                  !Mitch 3/22/2017
      COMMON /VPAI/PAI2,PAI,HPAI,RAD,akap,imod,iprp,island,imd,iprpp     &
                   ,nonln,igrav,isolv,ixmdf,iproc,imud,iwnd,depmin0
      common /origin/x0,y0,azimuth,isteer,iidate,sinaz,cosaz

! Fill in names with defaults in case card does not exist
      OptsFile   = 'options.std'
      DepFile    = 'dep.in'
      CurrFile   = 'current.in'
      EngInFile  = 'spec.in'
      EngOutFile = 'spec.out'
      WaveFile   = 'wavfld'
      ObsFile    = 'selhts.out'
      NestFile   = 'nest.out'                                           
      BreakFile  = 'break' 
      RadsFile   = 'radstress'
      StrucFile  = 'struct.dat'
      SurgeFile  = 'eta.in'
      MudFile    = 'mud.dat'
      FricFile   = 'friction.dat'
      FrflFile   = 'forward.dat'
      BrflFile   = 'backward.dat'
      WindFile   = 'wind.dat'
      SpecFile   = 'wave.spc'
      XMDFFile   = 'cms-wave_sol.h5'
      SHIPFILE   = 'shipwakeinfo.txt'
      SetupFile  = 'setup.wav'
      SeaFile    = 'sea.wav'
      SwellFile  = 'swell.wav'
!
      IF (SimFile(1:4) .NE. 'none') THEN
        OPEN (41, FILE = SimFile, STATUS = 'OLD')
! Make sure the file is an STWAVE sim file
        READ (41, '(A10,A80)',END = 200) text1, text2
        READ (text1,'(A6)') CARD
        READ (text2,*,end=200,err=200) x0,y0,azimuth
!       READ (41, '(A6,2x,3f15.5)', END = 200) CARD, x0, y0, azimuth
        cosaz=cos(azimuth*rad)
        sinaz=sin(azimuth*rad)
        IF (CARD(1:6).EQ.'STWAVE'.OR.CARD(1:5).EQ.'WABED'.OR.  &
        CARD(1:3).EQ.'CMS') THEN
10        READ (41, '(2A)', END = 100) CARD, FileName
          FileName = ADJUSTL(FileName)
          IF (CARD(1:4).EQ.'OPTS') THEN
             OptsFile   = FileName
          ELSEIF (CARD(1:3).EQ.'DEP') THEN
            DepFile    = FileName
          ELSEIF (CARD(1:4).EQ.'CURR') THEN
            CurrFile   = FileName
          ELSEIF (CARD(1:5).EQ.'SPEC ') THEN
            EngInFile  = FileName
          ELSEIF (CARD(1:4).EQ.'WAVE') THEN
            WaveFile   = FileName
          ELSEIF (CARD(1:4).EQ.'OBSE') THEN
            EngOutFile = FileName
          ELSEIF (CARD(1:5).EQ.'SELHT') THEN
            ObsFile = FileName
          ELSEIF (CARD(1:4).EQ.'NEST') THEN                              
            NestFile   = FileName
          ELSEIF (CARD(1:5).EQ.'BREAK') THEN
            BreakFile = FileName
          ELSEIF (CARD(1:4).EQ.'RADS') THEN
            RadsFile = FileName
          ELSEIF (CARD(1:5).EQ.'STRUC') THEN
            StrucFile = FileName
          ELSEIF (CARD(1:3).EQ.'SUR'.OR.CARD(1:3).EQ.'ETA') THEN
            SurgeFile = FileName
          ELSEIF (CARD(1:3).EQ.'MUD') THEN
            MudFile   = FileName
          ELSEIF (CARD(1:4).EQ.'FRIC') THEN
            FricFile  = FileName
          ELSEIF (CARD(1:4).EQ.'FREF') THEN
            FrflFile  = FileName
          ELSEIF (CARD(1:4).EQ.'BREF') THEN
            BrflFile  = FileName
          ELSEIF (CARD(1:4).EQ.'WIND') THEN
            WindFile  = FileName
          ELSEIF (CARD(1:5).EQ.'SPEC2') THEN
            SpecFile  = FileName
          ELSEIF (CARD(1:5).EQ.'SHIPW') THEN
                  ShipFile  = FileName
          ELSEIF (CARD(1:4).EQ.'XMDF') THEN
            XMDFFile  = FileName
          ENDIF
          GO TO 10
100     ENDIF
        ILOC=INDEX(SimFile,'/',BACK=.TRUE.)
        ThePath=NAME(1:ILOC)
        ILOC=INDEX(OptsFile,'/',BACK=.TRUE.)
        OptsFile=trim(ThePath)//OptsFile(ILOC+1:)//'                   '
        ILOC=INDEX(DepFile,'/',BACK=.TRUE.)
        DepFile=trim(ThePath)//DepFile(ILOC+1:)//'                     '
        ILOC=INDEX(CurrFile,'/',BACK=.TRUE.)
        CurrFile=trim(ThePath)//CurrFile(ILOC+1:)//'                   '
        ILOC=INDEX(EngInFile,'/',BACK=.TRUE.)
        EngInFile=trim(ThePath)//EngInFile(ILOC+1:)//'                 '
        ILOC=INDEX(WaveFile,'/',BACK=.TRUE.)
        WaveFile=trim(ThePath)//WaveFile(ILOC+1:)//'                   '
        ILOC=INDEX(EngOutFile,'/',BACK=.TRUE.)
        EngOutFile=trim(ThePath)//EngOutFile(ILOC+1:)//'               '
        ILOC=INDEX(NestFile,'/',BACK=.TRUE.)
        NestFile=trim(ThePath)//NestFile(ILOC+1:)//'                   '
        ILOC=INDEX(BreakFile,'/',BACK=.TRUE.)
        BreakFile=trim(ThePath)//BreakFile(ILOC+1:)//'                 '
        ILOC=INDEX(RadsFile,'/',BACK=.TRUE.)
        RadsFile=trim(ThePath)//RadsFile(ILOC+1:)//'                   '
        ILOC=INDEX(StrucFile,'/',BACK=.TRUE.)
        StrucFile=trim(ThePath)//StrucFile(ILOC+1:)//'                 '
        ILOC=INDEX(MudFile,'/',BACK=.TRUE.)
        MudFile=trim(ThePath)//MudFile(ILOC+1:)//'                     '
        ILOC=INDEX(FricFile,'/',BACK=.TRUE.)
        FricFile=trim(ThePath)//FricFile(ILOC+1:)//'                   '
        ILOC=INDEX(FrflFile,'/',BACK=.TRUE.)
        FrflFile=trim(ThePath)//FrflFile(ILOC+1:)//'                   '
        ILOC=INDEX(BrflFile,'/',BACK=.TRUE.)
        BrflFile=trim(ThePath)//BrflFile(ILOC+1:)//'                   '
        ILOC=INDEX(WindFile,'/',BACK=.TRUE.)
        WindFile=trim(ThePath)//WindFile(ILOC+1:)//'                   '
        ILOC=INDEX(SpecFile,'/',BACK=.TRUE.)
        SpecFile=trim(ThePath)//SpecFile(ILOC+1:)//'                   '
        ILOC=INDEX(ObsFile,'/',BACK=.TRUE.)
        ObsFile=trim(ThePath)//ObsFile(ILOC+1:)//'                     '
        ILOC=INDEX(ShipFile,'/',BACK=.TRUE.)
        ShipFile=trim(ThePath)//ShipFile(ILOC+1:)//'                   '
        ILOC=INDEX(XMDFFile,'/',BACK=.TRUE.)
        XMDFFile=trim(ThePath)//XMDFFile(ILOC+1:)//'                   '
        ILOC=INDEX(SurgeFile,'/',BACK=.TRUE.)
        SurgeFile=trim(ThePath)//SurgeFile(ILOC+1:)  &
                              //'                   '
200     CLOSE (41)
      ENDIF
      END SUBROUTINE
!
!********************************************************************************
      SUBROUTINE GetPath_inline (NAME, PATH)
!--GetPath -----------------------------------------------------
! PURPOSE - get the path from the name - used to get proj directory
! ARGUMENT DECLARATION
      CHARACTER*180 NAME,PATH
!---------------------------------------------------------------
      INTEGER COUNT
!
      COUNT = 180
      PATH = NAME

! PC
 200   IF (COUNT .GT. 0) THEN
         IF (PATH(COUNT:COUNT) .NE. '\') THEN
           PATH(COUNT:COUNT) = ' '
           COUNT = COUNT - 1
           GO TO 200
        ENDIF
      ENDIF

      RETURN
      END SUBROUTINE
!
!********************************************************************************
      SUBROUTINE StripPath_inline(FNAME)
!--StripPath----------------------------------------------------
! PURPOSE-Strip the path from the file name used
!---------------------------------------------------------------
      CHARACTER FNAME*180

      ISLASH = 0
      ILAST = 0
      ISTART = 0
      DO I = 1,179

! PC
        IF(FNAME(I:I).EQ.'\') ISLASH = I
! UNIX - LINUX
!       IF(FNAME(I:I).EQ.'/') ISLASH = I
! END
        IF(FNAME(I:I).NE.' ' .AND. ISTART.EQ.0    &
           .AND.  FNAME(I:I).NE.CHAR(9) )ISTART = I
      END DO
      DO I = 180,1,-1
        IF(FNAME(I:I).NE.' ' .AND. ILAST.EQ.0)ILAST = I
      END DO
      IF (ISLASH.GT.0) THEN
          FNAME = FNAME(ISLASH+1:ILAST)
      ELSE
          FNAME = FNAME(ISTART:ILAST)
      END IF

      RETURN
      END SUBROUTINE
!
!********************************************************************************
      SUBROUTINE SetPath_inline(PATH,FNAME)
!--SetPath------------------------------------------------------
! PURPOSE This subroutine takes a filename and a path and puts
!         the path at the beginning of the filename.
!---------------------------------------------------------------
      CHARACTER*180 PATH
      CHARACTER*180 FNAME
      CHARACTER*180 NEWNAME
      INTEGER COUNT,I
!
      COUNT = 180
888   IF (COUNT .GT. 0) then
        if (PATH(COUNT:COUNT) .EQ. ' ') THEN
          COUNT = COUNT - 1
          GO TO 888
        endif
      ENDIF
      COUNT = COUNT + 1
      I=1
      NEWNAME = PATH
998   IF (FNAME(I:I) .NE. ' ') THEN
        GO TO 999
      ELSE
        I = I + 1
      ENDIF
      GO TO 998
999   IF (COUNT .LT. 181 .AND. I .LT. 181) THEN
!       IF (FNAME(I:I) .NE. ' ') THEN
          NEWNAME(COUNT:COUNT) = FNAME(I:I)
          COUNT = COUNT + 1
!       ENDIF
        I = I + 1
        GO TO 999
      ENDIF
      FNAME = NEWNAME
      END SUBROUTINE
!

!********************************************************************************
      subroutine dissip1_inline(i)
      USE GLOBAL_INLINE, ONLY: NPF,MPD,IPMX,JPMX,IGPX,JGPX
      common /rsm/ cosa(mpd),sina(mpd),d1(ipmx,jpmx)
      common /rsw/ wk2(jpmx,npf,mpd),cgp(ipmx,jpmx)
      COMMON /DATA/NF,MD,IMAX,JMAX,IGMX,JGMX,JCPB,JCPE,JCB,JCE,NFF,MDD
      COMMON /DATD/DX,DY,DXX,DMESH,DTH,kdate,idate,depmax,depmax0
      COMMON /SPECA/SCP(JGPX,MPD),SI(JGPX,NPF,MPD)
      COMMON /FDS/FCN(NPF),DCM(MPD),DSFD(NPF,MPD)
      COMMON /rse/ex(ipmx,jpmx),ey(ipmx,jpmx)
      COMMON /VPAI/PAI2,PAI,HPAI,RAD,akap,imod,iprp,island,imd,iprpp     &
                   ,nonln,igrav,isolv,ixmdf,iproc,imud,iwnd,depmin0
      common /comsav/depmin,g,iview,iplane,iwave,cflat,irunup,iroll
      REAL, ALLOCATABLE :: exx(:),eyy(:)
      ALLOCATE (exx(ipmx),eyy(jpmx))

      do j=2,jgmx-1
        sumxx=.0
        sumyy=.0
        if(d1(i,j).gt..001) then
            d1ij=d1(i,j)
            cc=(d1(i,j)+d1(i,j-1)+d1(i,j+1))/3.
            if(cc.gt..001) d1ij=cc
          do l=imd,md
            do k=1,nf
              if (wk2(j, k, l) .gt. 1.0e-10) then
                tkd=min(2.*wk2(j,k,l)*d1ij,10.)
                fns=.5*(1.+2.*wk2(j,k,l)*d1ij/sinh(tkd))
                cc=si(j,k,l)*fns*pai2*fcn(k)/wk2(j,k,l)
                sumxx=sumxx+cc*cosa(l)
                sumyy=sumyy+cc*sina(l)
              end if
            end do
          end do
        end if
        if(iplane.le.1) then
          ex(i,j)=sumxx
          ey(i,j)=sumyy
        else
          ex(i,j)=ex(i,j)+sumxx
          ey(i,j)=ey(i,j)+sumyy
        end if
        exx(j)=ex(i,j)
        eyy(j)=ey(i,j)
      end do

      exx(1)=exx(2)
      eyy(1)=eyy(2)
      exx(jgmx)=exx(jgmx-1)
      eyy(jgmx)=eyy(jgmx-1)

      do j=2,jgmx-1
      ex(i,j)=(exx(j-1)+exx(j)+exx(j+1))/3.
      ey(i,j)=(eyy(j-1)+eyy(j)+eyy(j+1))/3.
      end do

      ex(i,1)=ex(i,2)
      ey(i,1)=ey(i,2)
      ex(i,jgmx)=ex(i,jgmx-1)
      ey(i,jgmx)=ey(i,jgmx-1)

      DEALLOCATE (exx,eyy)
      return
      end subroutine

!********************************************************************************
      subroutine dissip_inline(i)
      USE GLOBAL_INLINE, ONLY: NPF,MPD,IPMX,JPMX,IGPX,JGPX
      common /rsm/ cosa(mpd),sina(mpd),d1(ipmx,jpmx)
      common /rsw/ wk2(jpmx,npf,mpd),cgp(ipmx,jpmx)
      COMMON /DATA/NF,MD,IMAX,JMAX,IGMX,JGMX,JCPB,JCPE,JCB,JCE,NFF,MDD
      COMMON /DATC/TP,PRD,IBND,VL,TIDE,KPEAK,IBACK,NBLOCK,IWIND,WSMAG
      COMMON /DATD/DX,DY,DXX,DMESH,DTH,kdate,idate,depmax,depmax0
      COMMON /SPECA/SCP(JGPX,MPD),SI(JGPX,NPF,MPD)
      COMMON /rse/ex(ipmx,jpmx),ey(ipmx,jpmx)
      COMMON /VPAI/PAI2,PAI,HPAI,RAD,akap,imod,iprp,island,imd,iprpp     &
                   ,nonln,igrav,isolv,ixmdf,iproc,imud,iwnd,depmin0
      COMMON /BREK/DEPM(JGPX),DMNJ(JGPX),SLF(JGPX),wlmn(JGPX),cmn(jgpx)  &
                   ,sigm(jgpx),IWVBK,ibk3
      COMMON /WAVI/H13(IGPX,JGPX),T13(IGPX,JGPX),DMN(IGPX,JGPX)
      common /comsav/depmin,g,iview,iplane,iwave,cflat,irunup,iroll
      REAL, ALLOCATABLE :: exx(:),eyy(:)
      ALLOCATE (exx(ipmx),eyy(jpmx))

      do j=2,jgmx-1
        cgpeak=0.
         if (d1(i,j).gt..001) then
            d1ij=d1(i,j)
            cc=(d1(i,j)+d1(i,j-1)+d1(i,j+1))/3.
            if(cc.gt..001) d1ij=cc
            ipeak=int((dmnj(j)+hpai-dth)/dth+1.5)
            if (ipeak.lt.imd) ipeak=imd
            if (ipeak.gt.md) ipeak=md
            wkpeak=wk2(j,kpeak,ipeak)
            tkd=min(2.*wkpeak*d1ij,10.)
      cgpeak=h13(i,j)**2*hpai/t13(i,j)/wkpeak*(.125+.25*wkpeak*d1ij/sinh(tkd))
          endif
          if(iplane.le.1) then
            ex(i,j)=cgpeak*cos(dmnj(j))
            ey(i,j)=cgpeak*sin(dmnj(j))
          else
            ex(i,j)=ex(i,j)+cgpeak*cos(dmnj(j))
            ey(i,j)=ey(i,j)+cgpeak*sin(dmnj(j))
          end if
        exx(j)=ex(i,j)
        eyy(j)=ey(i,j)
       enddo

      exx(1)=exx(2)
      eyy(1)=eyy(2)
      exx(jgmx)=exx(jgmx-1)
      eyy(jgmx)=eyy(jgmx-1)

      do j=2,jgmx-1
      ex(i,j)=(exx(j-1)+exx(j)+exx(j+1))/3.
      ey(i,j)=(eyy(j-1)+eyy(j)+eyy(j+1))/3.
      end do

      ex(i,1)=ex(i,2)
      ey(i,1)=ey(i,2)
      ex(i,jgmx)=ex(i,jgmx-1)
      ey(i,jgmx)=ey(i,jgmx-1)

      DEALLOCATE (exx,eyy)
      RETURN
      END SUBROUTINE

!********************************************************************************
      subroutine dfds_inline
      USE GLOBAL_INLINE, ONLY: NPF,MPD,IPMX,JPMX,IGPX,JGPX
      COMMON /DATA/NF,MD,IMAX,JMAX,IGMX,JGMX,JCPB,JCPE,JCB,JCE,NFF,MDD
      COMMON /DATD/DX,DY,DXX,DMESH,DTH,kdate,idate,depmax,depmax0
      COMMON /DVAR/DVARX(IGPX),DVARY(JGPX),ETA(IGPX,JGPX),HSK(JGPX)
      COMMON /FDS/FCN(NPF),DCM(MPD),DSFD(NPF,MPD)
      COMMON /DFCN/DF(NPF),FFCN(NPF),DFINP(NPF),PL1E(NPF),PL1T(NPF)
      COMMON /WAVS/H13S(IGPX,JGPX),IBR(IGPX,JGPX),DISS(IGPX,JGPX)
      COMMON /rse/ex(ipmx,jpmx),ey(ipmx,jpmx)
      common /comsav/depmin,g,iview,iplane,iwave,cflat,irunup,iroll
!     dimension exx(ipmx,jpmx),eyy(ipmx,jpmx)
      REAL, ALLOCATABLE :: exx(:,:),eyy(:,:)
      ALLOCATE (exx(ipmx,jpmx),eyy(ipmx,jpmx))

      ni1=igmx-1
      nj1=jgmx-1
         
!     X-derivatives on x-boundaries:
      dx2=dvarx(1)+dvarx(2)
      dx2ni=dvarx(igmx)+dvarx(ni1)
      do j=1,jgmx
        exx(1,j)=(-3.*ex(1,j)+4.*ex(2,j)-ex(3,j))/dx2
        exx(igmx,j)=(3.*ex(igmx,j)-4.*ex(ni1,j)+ex(igmx-2,j))/dx2ni
      end do

!     Y-derivatives on y-boundaries:
      dy2=dvary(1)+dvary(2)
      dy2nj=dvary(jgmx)+dvary(nj1)
      do i=1,igmx
        eyy(i,1)=(-3.*ey(i,1)+4.*ey(i,2)-ey(i,3))/dy2
        eyy(i,jgmx)=(3.*ey(i,jgmx)-4.*ey(i,nj1)+ey(i,jgmx-2))/dy2nj
      end do

!     X-derivatives on internal grid pts & y-boundaries:
      do j=1,jgmx
        do i=2,igmx-1
          dx2=dvarx(i)+(dvarx(i-1)+dvarx(i+1))/2.
          exx(i,j)=(ex(i+1,j)-ex(i-1,j))/dx2
        end do
      end do

!     Y-derivatives on internal grid pts & x-boundaries:
      do i=1,igmx
       do j=2,nj1
         dy2=dvary(j)+(dvary(j-1)+dvary(j+1))/2.
         eyy(i,j)=(ey(i,j+1)-ey(i,j-1))/dy2
         diss(i,j)=(exx(i,j-1)+eyy(i,j-1)+exx(i,j)+eyy(i,j)+exx(i,j+1)+eyy(i,j+1))/3.
         if(isnan(diss(i,j))) diss(i,j)=0.
         if(diss(i,j).gt.0.) diss(i,j)=0.
	 end do
      end do

      do i=1,igmx
        diss(i,1) =diss(i,2)
        diss(i,jgmx)=diss(i,nj1)
      end do

      diss=diss*g

      if(iplane.eq.1) then
      exx=-ex
      eyy=-ey
      do j=1,jgmx
         do i=1,igmx
         ex(i,j)=exx(imax-i,jmax-j)
         ey(i,j)=eyy(imax-i,jmax-j)
         end do
      end do
      end if

      DEALLOCATE (exx,eyy)
      return
      end subroutine

!********************************************************************************
      subroutine sxycalc_inline(i)
      USE GLOBAL_INLINE, ONLY: NPF,MPD,IPMX,JPMX,IGPX,JGPX
      common /rsa/ sxx(ipmx,jpmx),sxy(ipmx,jpmx),syy(ipmx,jpmx)
      common /rsb/ wxrs(ipmx,jpmx),wyrs(ipmx,jpmx)
      common /rsc/ sxxx(ipmx,jpmx),sxyx(ipmx,jpmx)
      common /rsd/ sxyy(ipmx,jpmx),syyy(ipmx,jpmx)
      common /rsm/ cosa(mpd),sina(mpd),d1(ipmx,jpmx)
      common /rsw/ wk2(jpmx,npf,mpd),cgp(ipmx,jpmx)
      COMMON /DATA/NF,MD,IMAX,JMAX,IGMX,JGMX,JCPB,JCPE,JCB,JCE,NFF,MDD
      COMMON /DATC/TP,PRD,IBND,VL,TIDE,KPEAK,IBACK,NBLOCK,IWIND,WSMAG
      COMMON /DATD/DX,DY,DXX,DMESH,DTH,kdate,idate,depmax,depmax0
      COMMON /SPECA/SCP(JGPX,MPD),SI(JGPX,NPF,MPD)
      COMMON /VPAI/PAI2,PAI,HPAI,RAD,akap,imod,iprp,island,imd,iprpp     &
                   ,nonln,igrav,isolv,ixmdf,iproc,imud,iwnd,depmin0
      g=9.806
      nj1=jgmx-1

      do j=2,nj1
         sumxx=.0
         sumyy=.0
         sumxy=.0
         if(d1(i,j).gt..001) then
            d1ij=d1(i,j)
            cc=(d1(i,j)+d1(i,j-1)+d1(i,j+1))/3.
            if(cc.gt..01) d1ij=cc
            do l=imd,md
              do k=1,nf
                if (wk2(j, k, l) .gt. 1.0e-10) then
                  tkd=min(2.*wk2(j,k,l)*d1ij,10.)
                  fns=.5*(1.+2.*wk2(j,k,l)*d1ij/sinh(tkd))
                  sumxx=sumxx+si(j,k,l)*(fns*(1.+cosa(l)**2)-.5)
                  sumyy=sumyy+si(j,k,l)*(fns*(1.+sina(l)**2)-.5)
                  sumxy=sumxy+si(j,k,l)*fns*sina(l)*cosa(l)
                end if
              end do
            end do
         end if
           if(iback.eq.0) then
              sxx(i,j)=sumxx*g
              syy(i,j)=sumyy*g
              sxy(i,j)=sumxy*g
           else
             sxx(i,j)=sxx(i,j)+sumxx*g
             syy(i,j)=syy(i,j)+sumyy*g
             sxy(i,j)=sxy(i,j)-sumxy*g
           end if
      end do

      sxx(i,1)=sxx(i,2)
      syy(i,1)=syy(i,2)
      sxy(i,1)=sxy(i,2)
      sxx(i,jgmx)=sxx(i,nj1)
      syy(i,jgmx)=syy(i,nj1)
      sxy(i,jgmx)=sxy(i,nj1)

!     do j=1,jgmx
!     write(20,*) j,d1(i,j),sxx(i,j),sxy(i,j),syy(i,j)
!     end do

      return
      end subroutine

!********************************************************************************
      subroutine rstress_inline
      USE GLOBAL_INLINE, ONLY: NPF,MPD,IPMX,JPMX,IGPX,JGPX,KOMX
      common /rsa/ sxx(ipmx,jpmx),sxy(ipmx,jpmx),syy(ipmx,jpmx)
      common /rsb/ wxrs(ipmx,jpmx),wyrs(ipmx,jpmx)
      common /rsc/ sxxx(ipmx,jpmx),sxyx(ipmx,jpmx)
      common /rsd/ sxyy(ipmx,jpmx),syyy(ipmx,jpmx)
      common /rsm/ cosa(mpd),sina(mpd),d1(ipmx,jpmx)
      common /rsw/ wk2(jpmx,npf,mpd),cgp(ipmx,jpmx)
      COMMON /DATA/NF,MD,IMAX,JMAX,IGMX,JGMX,JCPB,JCPE,JCB,JCE,NFF,MDD
      COMMON /DATB/ICK3,ICK4,ICHOICE,HS0,wd,ws,ph0,wdd(mpd),aslop(ipmx)
      COMMON /DATC/TP,PRD,IBND,VL,TIDE,KPEAK,IBACK,NBLOCK,IWIND,WSMAG
      COMMON /DATD/DX,DY,DXX,DMESH,DTH,kdate,idate,depmax,depmax0
      COMMON /DVAR/DVARX(IGPX),DVARY(JGPX),ETA(IGPX,JGPX),HSK(JGPX)
      COMMON /WAVI/H13(IGPX,JGPX),T13(IGPX,JGPX),DMN(IGPX,JGPX)
      COMMON /WAVS/H13S(IGPX,JGPX),IBR(IGPX,JGPX),DISS(IGPX,JGPX)
      COMMON /VPAI/PAI2,PAI,HPAI,RAD,akap,imod,iprp,island,imd,iprpp     &
                   ,nonln,igrav,isolv,ixmdf,iproc,imud,iwnd,depmin0
      COMMON /OUTP/KOUT,IJSP(2,KOMX),IRS,IBREAK,ICUR,IWET,INST
      common /comsav/depmin,g,iview,iplane,iwave,cflat,irunup,iroll
      REAL, ALLOCATABLE :: exx(:,:),eyy(:,:),sr(:,:)
      ALLOCATE (exx(ipmx,jpmx),eyy(ipmx,jpmx),sr(ipmx,jpmx))
!
      br = 0.025
      ceff = 0.25
!
!     Intensity of roller is 0 (none) to 4 (strongest)
!     br is the wave-front slope (assumed to b 0.1 or less)
!     ceff is the fraction of the broken wave energy into roller
!     Dally and Brown put 1 but Tajima used 0.5
!     For iroll=4, br = 0.1 & ceff = 1
!
      Sr=0.0
!
      ni1=igmx-1
      nj1=jgmx-1
!
      if(ibreak.le.1) go to 111
      if(iroll.eq.0) go to 111
      br=br*float(iroll)
      ceff=ceff*float(iroll)
!
      do j=1,jgmx
        do i=1,igmx
          wkpeak=cgp(i,j)
          cgp(i,j)=pai2/(t13(i,j)+1.e-10)/(wkpeak+1.e-10)
          c1=sqrt(g*max(d1(i,j),.01))
          if(cgp(i,j).gt.c1) cgp(i,j)=c1
!          ak=sinh(min(2.*wkpeak*d1(i,j),10.))
!          cgp(i,j)=pai/t13(i,j)/wkpeak*(1.+2.*wkpeak*d1(i,j)/ak)
        enddo
      enddo
!
      do i=1,ni1
        i1=i+1
        do j=2,nj1
          dy2=dvary(j)+(dvary(j-1)+dvary(j+1))/2.
          cyij1=cgp(i,j+1)*sin(dmn(i,j+1)*rad)
          cyijp1=cgp(i,j-1)*sin(dmn(i,j-1)*rad)
          DcySrDy=(cyij1*Sr(i,j+1)-cyijp1*Sr(i,j-1))/dy2
          if(d1(i,j).lt..01) cycle
          if(h13(i,j).lt..01) cycle
          Dr=g*br*Sr(i,j)/cgp(i,j) !Stive and De Vriend
          Dw=ceff*diss(i1,j)
          Wij=-DcySrDy-Dr+Dw
          cxij=cgp(i,j)*cos(dmn(i,j)*rad)
          cxi1j=cgp(i1,j)*cos(dmn(i1,j)*rad)
          if(cxi1j.lt..5) cxi1j=.5
          if(dvarx(i).gt.10.) then
            if(ibr(i1,j).eq.0.and.Wij.lt.-1.e-7) Wij=-1.e-7
            if(ibr(i1,j-1)*ibr(i1,j+1).eq.0.and.Wij.lt.-1.e-7) Wij=-1.e-7
          end if
          Sr(i1,j)=(cxij*Sr(i,j)+dvarx(i)*Wij)/cxi1j
!         if(Sr(i1,j).gt.0.) Sr(i1,j)=0.
        enddo
        Sr(i1,1)=sr(i1,2)
        Sr(i1,jgmx)=sr(i1,nj1)
      enddo
!
      do i=1,igmx
        do j=1,jgmx
          sxx(i,j)=sxx(i,j)-Sr(i,j)*cos(dmn(i,j)*rad)**2
          sxy(i,j)=sxy(i,j)-Sr(i,j)*sin(2.*dmn(i,j)*rad)
          syy(i,j)=syy(i,j)-Sr(i,j)*sin(dmn(i,j)*rad)**2
        enddo
      enddo
!
  111 continue
!
!     X-derivatives on x-boundaries:
      dx2=dvarx(1)+dvarx(2)
      dx2ni=dvarx(igmx)+dvarx(ni1)
      do j=1,jgmx
        sxxx(1,j)=(-3.*sxx(1,j)+4.*sxx(2,j)-sxx(3,j))/dx2
        sxyx(1,j)=(-3.*sxy(1,j)+4.*sxy(2,j)-sxy(3,j))/dx2
        sxxx(igmx,j)=(3.*sxx(igmx,j)-4.*sxx(ni1,j)+sxx(igmx-2,j))/dx2ni
        sxyx(igmx,j)=(3.*sxy(igmx,j)-4.*sxy(ni1,j)+sxy(igmx-2,j))/dx2ni
      end do

!     Y-derivatives on y-boundaries:
      do i=1,igmx
        dy2=dvary(1)+dvary(2)
        dy2nj=dvary(jgmx)+dvary(nj1)
        sxyy(i,1)=(-3.*sxy(i,1)+4.*sxy(i,2)-sxy(i,3))/dy2
        syyy(i,1)=(-3.*syy(i,1)+4.*syy(i,2)-syy(i,3))/dy2
        sxyy(i,jgmx)=(3.*sxy(i,jgmx)-4.*sxy(i,nj1)+sxy(i,jgmx-2))/dy2nj
        syyy(i,jgmx)=(3.*syy(i,jgmx)-4.*syy(i,nj1)+syy(i,jgmx-2))/dy2nj
      end do

!     X-derivatives on internal grid pts & y-boundaries:
      do j=1,jgmx
        do i=2,igmx-1
          dx2=dvarx(i)+(dvarx(i-1)+dvarx(i+1))/2.
          sxxx(i,j)=(sxx(i+1,j)-sxx(i-1,j))/dx2
          sxyx(i,j)=(sxy(i+1,j)-sxy(i-1,j))/dx2
        end do
      end do

!     Y-derivatives on internal grid pts & x-boundaries:
      do j=2,jgmx-1
        do i=1,igmx
          dy2=dvary(j)+(dvary(j-1)+dvary(j+1))/2.
          sxyy(i,j)=(sxy(i,j+1)-sxy(i,j-1))/dy2
          syyy(i,j)=(syy(i,j+1)-syy(i,j-1))/dy2
        end do
      end do
!
      do i=1,igmx
        sxxx(i,1)=sxxx(i,2)
        sxyy(i,1)=sxyy(i,2)
        syyy(i,1)=syyy(i,2)
        sxyx(i,1)=sxyx(i,2)
        sxxx(i,jgmx)=sxxx(i,nj1)
        sxyy(i,jgmx)=sxyy(i,nj1)
        syyy(i,jgmx)=syyy(i,nj1)
        sxyx(i,jgmx)=sxyx(i,nj1)
      end do
!
      if(iplane.eq.2) then
        exx=wxrs
        eyy=wyrs
      end if
!
!     Radiation stress. Set stress values to zero on dry cells:
      do j=1,jgmx
        do i=1,igmx
          if (d1(i,j).lt..01) then
            wxrs(i,j)=.0
            wyrs(i,j)=.0
          else
            wxrs(i,j)=-sxxx(i,j)-sxyy(i,j)
            wyrs(i,j)=-sxyx(i,j)-syyy(i,j)
!
!     Rest stress values to a threshold on next-to-land cells:
!     modified by lihwa 27aug03
!     Check four neighbour cells (right, left, upper, lower)
            iadd1=i+1
            isub1=i-1
            jadd1=j+1
            jsub1=j-1
            if(iadd1.gt.igmx) iadd1=igmx
            if(isub1.lt.1)  isub1=1
            if(jadd1.gt.jgmx) jadd1=jgmx
            if(jsub1.lt.1)  jsub1=1
            cc=0.
            if(d1(iadd1,j).le..0) cc=cc+1.
            if(d1(isub1,j).le..0) cc=cc+1.
            if(d1(i,jadd1).le..0) cc=cc+1.
            if(d1(i,jsub1).le..0) cc=cc+1.
               if(cc.gt..5) then
               ss=sqrt(wxrs(i,j)**2+wyrs(i,j)**2+1.e-10)
               sdis=sqrt(dvarx(i)**2+dvary(j)**2)
               ssc=0.0005*sdis/cc/ss
               if(ssc.lt.1.) then
               wxrs(i,j)=wxrs(i,j)*ssc
               wyrs(i,j)=wyrs(i,j)*ssc
               end if
               go to 301
               end if
!     Check four diagonal neighbour cells
               if(d1(iadd1,jadd1).le..0) cc=cc+1.
               if(d1(iadd1,jsub1).le..0) cc=cc+1.
               if(d1(isub1,jadd1).le..0) cc=cc+1.
               if(d1(isub1,jsub1).le..0) cc=cc+1.
               if(cc.gt..5) then
               ss=sqrt(wxrs(i,j)**2+wyrs(i,j)**2+1.e-10)
               sdis=sqrt(dvarx(i)**2+dvary(j)**2)
               ssc=0.0005*sdis/cc/ss
               if(ssc.lt.1.) then
               wxrs(i,j)=wxrs(i,j)*ssc
               wyrs(i,j)=wyrs(i,j)*ssc
               end if 
               end if
  301          continue
!     end of modification

               if(ibr(i,j).eq.0) then
               cc1=d1(i,j)-hs0
               if(cc1.lt.0.) cc1=0.
               cc=exp(min(cc1,10.))
               wxrs(i,j)=wxrs(i,j)/cc
               wyrs(i,j)=wyrs(i,j)/cc
               end if
            endif
         enddo
      enddo
!
      if(iplane.eq.1) then
        exx=-wxrs
        eyy=-wyrs
        do j=1,jgmx
          do i=1,igmx
            wxrs(i,j)=exx(imax-i,jmax-j)
            wyrs(i,j)=eyy(imax-i,jmax-j)
          end do
        end do
      end if
!
      if(iplane.eq.2) then
        wxrs=wxrs+exx
        wyrs=wyrs+eyy
      end if
!
      DEALLOCATE (exx,eyy,sr)
      return
      end subroutine

!********************************************************************************
      SUBROUTINE RUNNING_TIME_inline(TIME_BEGIN,TIME_END)

      REAL TIME_BEGIN,TIME_END,RTIME
      INTEGER HH,MM
	REAL SS

      RTIME=TIME_END-TIME_BEGIN
      HH=RTIME/3600.
      MM=(RTIME-HH*3600.)/60.
      SS=RTIME-HH*3600-MM*60

      WRITE(*,200) HH,MM,SS
 200  FORMAT(' CPU-TIME=',1X,I4,'h',I2,'m',F10.7,'s')
      
      RETURN
      END SUBROUTINE

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!  CALCULATION OF Matrix BY GAUS-SEIDEL METHOD---revised version
!  By Zhang & Wu, NCCHE, Oct. 1, 2009
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE GSR_inline(II,JB,JE,NMX,MARK)
      USE GLOBAL_INLINE, ONLY: NPF,MPD,IPMX,JPMX,IGPX,JGPX,MPMX
      COMMON /DATA/NF,MD,IMAX,JMAX,IGMX,JGMX,JCPB,JCPE,JCB,JCE,NFF,MDD
      COMMON /DATB/ICK3,ICK4,ICHOICE,HS0,wd,ws,ph0,wdd(mpd),aslop(ipmx)
      COMMON /DATD/DX,DY,DXX,DMESH,DTH,kdate,idate,depmax,depmax0
      COMMON /SPECA/SCP(JGPX,MPD),SI(JGPX,NPF,MPD)
      COMMON /AB/AA(5,MPMX),IA(5,MPMX),B(MPMX),X(MPMX)
      COMMON /VPAI/PAI2,PAI,HPAI,RAD,akap,imod,iprp,island,imd,iprpp     &
                   ,nonln,igrav,isolv,ixmdf,iproc,imud,iwnd,depmin0
!     dimension x0(mpmx)
      REAL,ALLOCATABLE :: X0(:)
      ALLOCATE (X0(mpmx))
!
!     ifast=0
!     if(iprpp.ne.-1.and.dxx.gt.300.) ifast=1
!
      LIM=1000
      DLTA=0.001
      if(dmesh.le..5) then
        LIM=10000
        DLTA=.000001
      end if
      ws1000=ws*1000.
!
      JUDG=0
      XXMAX=0
      DO I=1,NMX            
         IF(ABS(AA(1,I)).LT.1.0E-4) THEN
           X(I)=0.0
         ELSE
           X(I)=B(I)/AA(1,I)
         ENDIF
         IF(X(I).GT.1.0E-20) JUDG=1
         X0(I)=X(I)
         ABXX=ABS(X(I))
         IF(XXMAX.LT.ABXX) XXMAX=ABXX
        ENDDO
!
      IF(JUDG.EQ.0) GOTO 11
      XLIM=XXMAX/1000.0
      ICC=0
      XMAX=100.0
!
!----Iteration loop
        IO=0
      DO WHILE (XMAX.GT.DLTA.AND.ICC.LT.LIM)
        ICC=ICC+1
        XMAX=0.0
!$omp parallel
!$omp do private (S, I, J, PX, PPX) REDUCTION(MAX:XMAX)
        DO I=1,NMX
         IF(ABS(AA(1,I)).GE.1.0E-4) THEN
           IP=0
           S=0.
           IF(I.NE.IA(1,I)) THEN
           IO=1
           ELSE
             IP=1
           ENDIF
!
           DO J=2,5
              IF(ABS(AA(J,I)).GE.1.0E-4) THEN
                IF(I.NE.IA(J,I).AND.IA(J,I).GT.0) THEN
                  IF(x(IA(J,I)).gt.1.0E-18) then
                    IF(abs(s).gt.100.) s=0.
                    S=S+AA(J,I)*X(IA(J,I))
                  ENDIF
                ENDIF
              ENDIF
           ENDDO
!           IF(IP.EQ.0) GOTO 50
           PX=X(I)
           if(abs(AA(1,I)).LT.1.0E-4) then
           XX=0.0
           else
           XX=(B(I)-S)/AA(1,I)
           end if
!       if(ws.ge..1) then
!         if(ifast.eq.1) then
!           if(abs(xx-px).gt.x0(i)*ws1000) xx=px
!         end if
!       end if
           X(I)=XX
           IF(ABS(PX).GE.XLIM) THEN
             PPX=(XX-PX)/PX
             PPX=ABS(PPX)
             IF(XMAX.LT.PPX) XMAX=PPX
           ENDIF
        ENDIF
     ENDDO
!$omp end do
!$omp end parallel
     IF(IO.EQ.1) THEN
     MARK=1
     GOTO 50
     END IF
   ENDDO
     IF(ICC.GE.LIM) THEN
     MARK=2
     GOTO 60
     END IF
!
   11 MX=0
      DO J=JB,JE
        DO M=Imd,MD
          MX=MX+1
          IF(X(MX).gt.5000.) X(MX)=0.
            SCP(J,M)=X(MX)
        ENDDO
      ENDDO
      GOTO 90
!
   50 continue !write(*,*) 'GSR TYPE 1:',II
      GOTO 90
   60 continue !write(*,*) 'GSR TYPE 2:',II,ICC,XMAX
   90 DEALLOCATE (x0)
!
      RETURN
      END SUBROUTINE
!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!  CALCULATION OF Matrix BY TDMA and Line by line
!  By Zhang & Wu, NCCHE, Oct. 1, 2009
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&      
      SUBROUTINE ADI_inline(II,JB,JE,NMX,MARK)
      USE GLOBAL_INLINE, ONLY: NPF,MPD,IPMX,JPMX,IGPX,JGPX,MPMX
      COMMON /DATA/NF,MD,IMAX,JMAX,IGMX,JGMX,JCPB,JCPE,JCB,JCE,NFF,MDD
      COMMON /DATB/ICK3,ICK4,ICHOICE,HS0,wd,ws,ph0,wdd(mpd),aslop(ipmx)
      COMMON /DATD/DX,DY,DXX,DMESH,DTH,kdate,idate,depmax,depmax0
      COMMON /SPECA/SCP(JGPX,MPD),SI(JGPX,NPF,MPD)
      COMMON /AB/AA(5,MPMX),IA(5,MPMX),B(MPMX),X(MPMX)
      COMMON /VPAI/PAI2,PAI,HPAI,RAD,akap,imod,iprp,island,imd,iprpp     &
                   ,nonln,igrav,isolv,ixmdf,iproc,imud,iwnd,depmin0
      REAL, ALLOCATABLE :: AN(:,:), AS(:,:), AE(:,:), AW(:,:), AP(:,:)
      REAL, ALLOCATABLE :: SU(:,:), F(:,:), F0(:,:)
      REAL, ALLOCATABLE :: AA1(:), BB1(:), CC1(:), DD1(:)
      INTEGER,ALLOCATABLE :: JJ(:,:), MM(:,:), NE(:,:), NW(:,:)      
      ALLOCATE (AN(JPMX,MPD),AS(JPMX,MPD),AE(JPMX,MPD),AW(JPMX,MPD)      &
               ,AP(JPMX,MPD),SU(JPMX,MPD),F(0:JPMX,0:MPD),F0(0:JPMX,0:MPD))
      ALLOCATE (JJ(JPMX,MPD),MM(JPMX,MPD),NE(JPMX,MPD),NW(JPMX,MPD))                                                                 
      ALLOCATE (AA1(0:JPMX),BB1(JPMX),CC1(0:JPMX),DD1(JPMX))

      LIM=500
      DLTA=0.001
      IF(dmesh.LE..5) THEN
        LIM=5000
        DLTA=.000001
      ENDIF

      !--- coefficient matrix
      F=0.0
      MX=0
      JUDG=0
      DO J=JB,JE
        DO M=Imd,MD
          MX=MX+1          
          AN(J,M)=-AA(5,MX)
          AS(J,M)=-AA(4,MX)
          AE(J,M)=-AA(3,MX)
          AW(J,M)=-AA(2,MX)
          AP(J,M)=AA(1,MX)
          SU(J,M)=B(MX)
          IF(AP(J,M).GT.1.0E-4) F(J,M)=B(MX)/AP(J,M)
          IF(F(J,M).GT.1.0E-20) JUDG=1 
           
          !----find reflection node                 
          NE(J,M)=IA(3,MX)
          NW(J,M)=IA(2,MX)
          JJ(J,M)=J
          MM(J,M)=M                         
          IF(J.EQ.JE.AND.NE(J,M).GT.0) THEN
            JJ(J,M)=NE(J,M)/(MD-Imd+1)
            MM(J,M)=NE(J,M)-JJ(J,M)*(MD-Imd+1)  
            IF(MM(J,M).EQ.0)THEN
               MM(J,M)=MD-Imd+1
            ELSE
               JJ(J,M)=JJ(J,M)+1
            ENDIF           
          ENDIF 
          IF(J.EQ.JB.AND.NW(J,M).GT.0) THEN
            JJ(J,M)=NW(J,M)/(MD-Imd+1)
            MM(J,M)=NW(J,M)-JJ(J,M)*(MD-Imd+1)  
            IF(MM(J,M).EQ.0)THEN
               MM(J,M)=MD-Imd+1
            ELSE
               JJ(J,M)=JJ(J,M)+1
            ENDIF 
          ENDIF                    
        ENDDO
      ENDDO
      IF(JUDG.EQ.0) GOTO 11

      !----Iteration loop 
      DO Iter=1,LIM
        XMAX=0.0
         
        !C-- LINE-RELAXATION USING TDMA ( AT M DIRECTION )
        AA1(Imd-1)=0.
        DO J=JB,JE
           CC1(Imd-1)=0.0    
          DO M=Imd,MD
             F0(J,M)=F(J,M)                       
             AA1(M)=AN(J,M)
             BB1(M)=AS(J,M)
             IF(J.EQ.JE.AND.NE(J,M).GT.0) F(J+1,M)=F(JJ(J,M),MM(J,M))
             IF(J.EQ.JB.AND.NW(J,M).GT.0) F(J-1,M)=F(JJ(J,M),MM(J,M))
             CC1(M)=AW(J,M)*F(J-1,M)+AE(J,M)*F(J+1,M)+SU(J,M)          
             DD1(M)=AP(J,M)
             TERM1=DD1(M)-BB1(M)*AA1(M-1)+1.E-10
             TERM =1./TERM1
             AA1(M)=AA1(M)*TERM
             CC1(M)=(CC1(M)+BB1(M)*CC1(M-1))*TERM
          ENDDO
          DO M=MD,Imd,-1
            F(J,M)=AA1(M)*F(J,M+1)+CC1(M)
          ENDDO
        ENDDO

        !----  ( AT J DIRECTION )
        AA1(JB-1)=0.
        DO M=Imd,MD           
          CC1(JB-1)=0.0     
          IF(NW(JB,M).GT.0) CC1(JB-1)=F(JJ(JB,M),MM(JB,M))                
          DO J=JB,JE
             AA1(J)=AE(J,M)
             BB1(J)=AW(J,M)
             CC1(J)=AS(J,M)*F(J,M-1)+AN(J,M)*F(J,M+1)+SU(J,M)
             DD1(J)=AP(J,M)
             TERM1=DD1(J)-BB1(J)*AA1(J-1)+1.E-10
             TERM =1./TERM1
             AA1(J)=AA1(J)*TERM
             CC1(J)=(CC1(J)+BB1(J)*CC1(J-1))*TERM
          ENDDO
          DO J=JE,JB,-1
             IF(J.EQ.JE.AND.NE(J,M).GT.0) F(J+1,M)=F(JJ(J,M),MM(J,M))
             F(J,M)=AA1(J)*F(J+1,M)+CC1(J)
             IF(abs(F(J,M)).GT.5000.) F(J,M)=0.

             !.....CHECK CONVERGENCE OF INNER ITERATIONS
             PPX=(F(J,M)-F0(J,M))/(F0(J,M)+1.0E-20)
             PPX=ABS(PPX)
             IF(XMAX.LT.PPX) XMAX=PPX
          ENDDO
        ENDDO
        IF(XMAX.LT.DLTA) EXIT
      ENDDO
      IF(Iter.EQ.LIM) THEN
      MARK=2
      GOTO 60
      END IF

 11   CONTINUE
      DO M=IMD,MD
         DO J=JB,JE
            IF(F(J,M).gt.5000.) F(J,M)=0.
            SCP(J,M)=F(J,M)
         ENDDO
      ENDDO

      GOTO 90
 60   WRITE(*,*) 'ADI TYPE 2',II,Iter,XMAX
 90   DEALLOCATE (AN,AS,AE,AW,AP,SU,F,F0,AA1,BB1,CC1,DD1,JJ,MM,NE,NW)  

      RETURN
      END SUBROUTINE
    
!***********************************************************************   
    subroutine cmswave_cards(aCard,foundcard)
! Reads cards from the CMS-Wave options file (*.std)
! written by Mitchell Brown, USACE-CHL  10/18/2021
!***********************************************************************
    use diag_lib,      only: diag_print_warning, diag_print_error
    use global_inline, only: gamma_bj78, KOMX, JGPX, IGPX
    
    implicit integer (i-n), real (a-h,o-z)
    
    COMMON /VPAI/PAI2,PAI,HPAI,RAD,akap,imod,iprp,island,imd,iprpp,     &
                 nonln,igrav,isolv,ixmdf,iproc,imud,iwnd,depmin0
    COMMON /OUTP/KOUT,IJSP(2,KOMX),IRS,IBREAK,ICUR,IWET,INST
    COMMON /DATC/TP,PRD,IBND,VL,TIDE,KPEAK,IBACK,NBLOCK,IWIND,WSMAG
    COMMON /wavenum/ itms,ibf,iark,iarkr,bf,ark,arkr                 !Wu
    COMMON /BREK/DEPM(JGPX),DMNJ(JGPX),SLF(JGPX),wlmn(JGPX),cmn(JGPX),  &
                 sigm(JGPX),IWVBK,ibk3    
    COMMON /comsav/depmin,g,iview,iplane,iwave,cflat,irunup,iroll
    COMMON /OUTN/nest,inest(komx),jnest(komx),ix1(igpx),ix2(igpx)
    
    character(len=*),intent(inout) :: aCard
    logical,intent(out)            :: foundcard
    
    character*80  :: aVal
    character*200 :: text
    
    foundcard = .true.
    
    select case (aCard)
    case('WV_PROPAGATION_TYPE') 
      backspace(11)
      read(11,*) aCard, aVal
      if     (aVal == 'WIND_AND_SPECTRA') then        ! 0 - wave generation and spectra (use wind if provided)
        iprpp = 0
      elseif (aVal == 'SPECTRA_ONLY') then            ! 1 - propagation with spectra only (neglect wind input)
        iprpp = 1
      elseif (aVal == 'FAST-MODE') then               !-1 - fast-mode (wave generation and spectra)
        iprpp = -1    
      elseif (aVal == 'FAST-MODE_SPECTRA_ONLY') then  !-2 - fast-mode w/spectra only
        iprpp = -2
      else
        call diag_print_error ('Bad selection for WV_PROPAGATION_TYPE')
      endif
    case('WV_PLANE_DEFINITION')
      backspace(11)
      read(11,*) aCard, aVal
      if     (aVal == 'HALF-PLANE') then
        iview = 0
      elseif (aVal == 'FULL-PLANE') then
        iview = 1
      elseif (aVal == 'FULL-PLANE_WITH_REVERSE_SPECTRA') then
        iview = 2
      else
        call diag_print_error ('Bad selection for WV_PLANE_DEFINITION')
      endif
    case('WV_MATRIX_SOLVER')
      backspace(11)
      read(11,*) aCard, aVal
      if     (aVal == 'GAUSS-SEIDEL') then
        isolv = 0
      elseif (aVal == 'ADI') then
        isolv = 1
      else
        call diag_print_error ('Bad selection for WV_MATRIX_SOLVER')
      endif
    case('WV_CURRENT_TYPE')       
      backspace(11)
      read(11,*) aCard, aVal
      if      (aVal == 'OFF') then                 !0 - no current
        icur = 0
      else if (aVal == 'MULTIPLE') then            !1 - multiple currents in seq. order
        icur = 1
      else if (aVal == 'SINGLE') then              !2 - read only first current record
        icur = 2
      else
        call diag_print_error ('Bad selection for WV_CURRENT_TYPE')
      endif
    case('WV_BOUNDARY_NESTING')   
      backspace(11)
      read(11,*) aCard, aVal
      if      (aVal == 'OFF') then                 !0 - no nesting
        ibnd = 0
      else if (aVal == 'LINEAR') then              !1 - linear interp for 'n' points
        ibnd = 1
      else if (aVal == 'MORPHIC') then             !2 - morphic interp for 'n' points
        ibnd = 2
      else
        call diag_print_error ('Bad selection for WV_BOUNDARY_NESTING')
      endif
    case('WV_BOTTOM_FRICTION')    
      backspace(11)
      read(11,*) aCard, aVal                       
      if     (aVal == 'OFF') then                        !0 - no bottom friction
        ibf = 0
      elseif (aVal == 'CONSTANT_DARCY_WEISBACH') then    !1 - constant Darcy-Weisbach coefficient (=bf)
        ibf = 1
      elseif (aVal == 'VARIABLE_DARCY_WEISBACH') then    !2 - variable Darcy-Weisbach coefficient (friction.dat)
        ibf = 2          
      elseif (aVal == 'CONSTANT_MANNINGS') then          !3 - constant Manning coefficient (=bf)
        ibf = 3
      elseif (aVal == 'VARIABLE_MANNINGS') then          !4 - variable Manning coefficient (friction.dat)
        ibf = 4
      else
        call diag_print_error ('Bad selection for WV_BOTTOM_FRICTION')
      endif
    case('WV_FWD_REFLECTION')     
      backspace(11)
      read(11,*) aCard, aVal
      if     (aVal == 'OFF') then      !0 - no forward reflection
        iark = 0                                   
      elseif (aVal == 'CONSTANT') then !1 - Use a constant coefficient with 'WV_FWD_REFLECTION_COEFF' card.
        iark = 1                                   
      elseif (aVal == 'VARIABLE') then !2 - Use a variable coefficient with values stored in '<project>.fref' file
        iark = 2                                  
      else
        call diag_print_error ('Bad selection for WV_FWD_REFLECTION')
      endif
    case('WV_BWD_REFLECTION')     
      backspace(11)
      read(11,*) aCard, aVal
      if     (aVal == 'OFF') then      !0 - no backward reflection
        iarkr = 0                                   
      elseif (aVal == 'CONSTANT') then !1 - Use a constant coefficient with 'WV_BWD_REFLECTION_COEFF' card.
        iarkr = 1                                   
      elseif (aVal == 'VARIABLE') then !2 - Use a variable coefficient with values stored in '<project>.bref' file
        iarkr = 2                                  
      else
        call diag_print_error ('Bad selection for WV_BWD_REFLECTION')
      endif
    case('WV_WETTING_DRYING')     
      backspace(11)
      read(11,*) aCard, aVal
      if     (aVal == 'ON') then                       ! 0 - normal wetting/drying
        iwet = 0                                   
      elseif (aVal == 'OFF') then                      ! 1 - no wetting/drying
        iwet = 1                                   
      elseif (aVal == 'ON_WITH_SEA-SWELL_FILES') then  !-1 - normal wet/dry, output swell and local sea files
        iwet = -1                                  
      else
        call diag_print_error ('Bad selection for WV_WETTING_DRYING')
      endif
    case('WV_DIFFRACTION_INTENSITY')
      backspace(11)
      read(11,*) aCard, akap
    case('WV_NUM_THREADS')
      backspace(11)
      read(11,*) aCard, iproc
      
            
    case('WV_NESTING_CELLS') 
      backspace(11)
      read(11,*) aCard, nest                       !0 - no nest cells
      if(nest.ge.1) then                           !n - list of nesting cells to read in  
        backspace(11)
        read (11,'(A)') text
        read (text,*) aCard, nest, (inest(nn),jnest(nn),nn=1,nest)   !This reads card, total, and all values from one line
      end if
    case('WV_OBSERVATION_CELLS') 
      backspace(11)
      read(11,*) aCard, kout                       !0 - no obs output
      if(kout.ge.1) then                           !n - output of spectra and parameters at 'n' selected cells
        backspace(11)
        read (11,'(A)') text 
        read (text,*) aCard, kout, (ijsp(1,nn),ijsp(2,nn),nn=1,kout)  !This reads card, total, and all values from one line
      end if

    
    case('WV_BOTTOM_FRICTION_COEFF')
      backspace(11)
      read(11,*) aCard, bf
    case('WV_FWD_REFLECTION_COEFF')
      backspace(11)
      read(11,*) aCard, ark                        !limit 0.0 <= ark < 1.0
    case('WV_BWD_REFLECTION_COEFF')
      backspace(11)
      read(11,*) aCard, arkr                       !limit 0.0 <= arkr < 1.0
      
      
    case('WV_BREAKING_FORMULA')
      backspace(11)
      read(11,*) aCard, aVal
      if     (aVal == 'EXT_GODA') then             !0 - Extended Goda
        iwvbk = 0
      elseif (aVal == 'EXT_MICHE') then            !1 - Extended Miche
        iwvbk = 1
      elseif (aVal == 'BATTJES_JANSSEN_78') then   !2 - Battjes and Janssen (1978)
        iwvbk = 2
      elseif (aVal == 'CHAWLA_KIRBY') then         !3 - Chawla and Kirby
        iwvbk = 3  
      elseif (aVal == 'BATTJES_JANSSEN_07') then   !4 - Battjes and Janssen (2007)
        iwvbk = 4
      else
        call diag_print_error ('Bad selection for WV_BREAKING_FORMULA')
      endif
    !Added to have user-specified Gamma for certain breaking formulae
    case('WV_SET_GAMMA_BJ78')
      backspace(11)
      read(11,*) aCard, gamma_bj78     !initial restriction: 0.4 <= gamma_bj78 <= 0.8
      if ((gamma_bj78 .lt. 0.4) .or. (gamma_bj78 .gt. 0.8)) then
        call diag_print_warning('WV_SET_GAMMA_BJ78 - Outside normal range of 0.4 <= [value] <= 0.8')
      else if ((gamma_bj78 .le. 0.0) .or. (gamma_bj78 .ge. 1.0)) then
        call diag_print_error  ('WV_SET_GAMMA_BJ78 - Invalid value')
      endif
      
    case('WV_ROLLER_EFFECT')
      backspace(11)
      read(11,*) aCard, aVal
      if     (aVal == 'OFF') then
        iroll = 0
      elseif (aVal == '25_PERCENT') then
        iroll = 1
      elseif (aVal == '50_PERCENT') then
        iroll = 2
      elseif (aVal == '75_PERCENT') then
        iroll = 3
      elseif (aVal == '100_PERCENT') then
        iroll = 4
      else
        call diag_print_error ('Bad selection for WV_ROLLER_EFFECT')
      endif
      

    case('WV_ENABLE_INFRAGRAVITY')
      backspace(11)
      read(11,*) aCard, aVal
      if     (aVal == 'OFF') then
        igrav = 0
      elseif (aVal == 'ON') then
        igrav = 1
      else
        call diag_print_error ('Bad selection for WV_ENABLE_INFRAGRAVITY')
      endif
    case('WV_ENABLE_MUDDY_BED')  
      backspace(11)
      read(11,*) aCard, aVal
      if     (aVal == 'OFF') then
        imud = 1
      elseif (aVal == 'ON') then
        imud = 0                 !opposite to what seems natural
      else
        call diag_print_error ('Bad selection for WV_ENABLE_MUDDY_BED')
      endif
    case('WV_ENABLE_NONLINEAR_WAVES')
      backspace(11)
      read(11,*) aCard, aVal
      if     (aVal == 'OFF') then
        nonln = 0
      elseif (aVal == 'ON') then
        nonln = 1                  
      else
        call diag_print_error ('Bad selection for WV_ENABLE_NONLINEAR_WAVES')
      endif
    !case('WV_ENABLE_ROLLER')
    !  backspace(11)
    !  read(11,*) aCard, aVal
    !  if     (aVal == 'OFF') then
    !    iroll = 0
    !  elseif (aVal == 'ON') then
    !    iroll = 1
    !  else
    !    call diag_print_error ('Bad selection for WV_ENABLE_ROLLER')
    !  endif
    case('WV_ENABLE_RUNUP')
      backspace(11)
      read(11,*) aCard, aVal
      if     (aVal == 'OFF') then
        irunup = 0
      elseif (aVal == 'REL_TO_ABS_DATUM') then
        irunup = 1
      elseif (aVal == 'RED_TO_UPDATED_MWL') then
        irunup = 2
      else
        call diag_print_error ('Bad selection for WV_ENABLE_RUNUP')
      endif
    case('WV_ENABLE_WIND')
      backspace(11)
      read(11,*) aCard, aVal
      if     (aVal == 'OFF') then
        iwnd = 1
      elseif (aVal == 'ON') then
        iwnd = 0                     !opposite to what seems natural
      else
        call diag_print_error ('Bad selection for WV_ENABLE_WIND')
      endif

    
    case('WV_BREAKING_OUTPUT')    
      backspace(11)
      read(11,*) aCard, aVal
      if     (aVal == 'OFF') then
        ibreak = 0
      elseif (aVal == 'BREAKING_INDICES') then
        ibreak = 1
      elseif (aVal == 'DISSIPATION_VALUES') then
        ibreak = 2   
      else
        call diag_print_error ('Bad selection for WV_BREAKING_OUTPUT')
      endif
    case('WV_RAD_STRESS_OUTPUT')   
      backspace(11)
      read(11,*) aCard, aVal
      if     (aVal == 'OFF') then                  !0 - no rad output 
        irs = 0
      elseif (aVal == 'RAD_STRESS_FILE') then      !1 - output wave radiation stresses only (*.rad)
        irs = 1
      elseif (aVal == 'STRESS_SETUP_FILES') then   !2 - output wave radiation stresses and setup/maximum water level
        irs = 2      
      else
        call diag_print_error ('Bad selection for WV_RAD_STRESS_OUTPUT')
      endif
    case('WV_OUTPUT_XMDF')
      backspace(11)
      read(11,*) aCard, aVal
      if     (aVal == 'OFF') then
        ixmdf = 0                                  !0 - ASCII OUTPUT ONLY
      elseif (aVal == 'ON') then       
        ixmdf = 1                                  !1 - XMDF OUTPUT ONLY
      else
        call diag_print_error ('Bad selection for WV_OUTPUT_XMDF')
      endif
      
        
    case default
      foundcard = .false.

    end select
    
    return
    end subroutine cmswave_cards
