C============================================================================
      module nos_tide
C  NOS Tidal Module
C  Modified from NOS tidal prediction software
C
C Technical Contact(s):   Name:  Aijun Zhang             Org:  NOS/CO-OPS
C                         Phone: 301-713-2890x127      E-Mail: aijun.zhang@noaa.gov
C
C Abstract:
C         This program is modified from pred.f of Chris Zervas so that 
C		it can make prediction of multiple years. 
C		Change call CONCTJ and CONJTC to call julian.
C  		Also call equarg.f to calculate XODE and VPU, instead of reading 
C		from data file 'yr'.
C
C
C Usage:   nos_pred "$BEGINDATE" "$ENDDATE" $KINDAT $DELT $CONV $XMAJOR $filein $fileout
C
C	   Called by 
C
C Input Parameters:  
C		BEGINDATE="200801011230"
C         	ENDDATE=  "200812311230"
C         	KINDAT=1, for current prediction; =2 for water level prediction
C         	DELT is time interval of output time series in hours
C         	CONV: Units convertion of predicted variable
C         	XMAJOR is principle current direction in degrees
C         	filein is input file name which includes tide constituents
C         	fileout is output name which contains predicted water level or current time series      
C============================================================================
      implicit none      
      contains
C***********************************************************************      
      subroutine nos_tide_init(iyr,cname,fn,eqarg)
C***********************************************************************      
      use prec_def
      implicit none
      !Input/Output
      integer,intent(in) :: iyr
      character(len=*) :: cname(37)
      real(ikind),intent(out) :: fn(37),eqarg(37)
      !Internal variables
      integer :: i
      
CCCCCCCCCCCCCC  used by equarg.f
      real spd(180),fff(180),vau(180),VPU(180),XODE(180),a(180)
      character*10 labl(180),ALIST(37)
      DATA (ALIST(I),I=1,37) /'M(2)      ','S(2)      ','N(2)      ',
     1                        'K(1)      ','M(4)      ','O(1)      ',
     2                        'M(6)      ','MK(3)     ','S(4)      ',
     3                        'MN(4)     ','NU(2)     ','S(6)      ',
     4                        'MU(2)     ','2N(2)     ','OO(1)     ',
     5                        'LAMBDA(2) ','S(1)      ','M(1)      ',
     6                        'J(1)      ','MM        ','SSA       ',
     7                        'SA        ','MSF       ','MF        ',
     8                        'RHO(1)    ','Q(1)      ','T(2)      ',
     9                        'R(2)      ','2Q(1)     ','P(1)      ',
     1                        '2SM(2)    ','M(3)      ','L(2)      ',
     2                        '2MK(3)   ','K(2)      ','M(8)      ',
     3                        'MS(4)     '/
      DATA (A(I), I=1,37)/ 28.9841042,  30.0000000,  28.4397295,
     115.0410686,57.9682084,13.9430356,86.9523127,44.0251729,
     260.0000000,57.4238337,28.5125831,90.0000000,27.9682084,
     327.8953548,16.1391017,29.4556253,15.0000000,14.4966939,
     415.5854433, 0.5443747, 0.0821373, 0.0410686, 1.0158958,
     5 1.0980331,13.4715145,13.3986609,29.9589333,30.0410667,
     612.8542862,14.9589314,31.0158958,43.4761563,29.5284789,
     742.9271398,30.0821373, 115.9364169,58.9841042/

CCCCCCCCCCCCCC


C node factor and equilibrium arguments for the middle of each year (day 183 or 184) are used in 
C the same calender year regardless of the length of time series. This is consistent with CO-OPS
C tidal prediction programs
      call equarg(37,iyr,1,1,365,labl(1),fff(1),vau(1))

      do i=1,37
        cname(i) = labl(i)
        fn(i) = fff(i)
        eqarg(i) = vau(i)
      enddo
      
      return
      endsubroutine
      
***************************************************************************      
      subroutine equarg(nsped,IYR,ICM,ICD,length,label,fff,vau)
c
c     program to calculate equilibrium arguments and node factors
c     by e.e.long      (slight revision by b.b.parker)
c       major revisions by len hickman 6/11/86 to make this a subroutine
c       of lsqha.  v is calculated for the beginning of the series.
c       u and f are adjusted to the midpoint of the series (regardless
c       of length).
*     More revisions by Geoff French 9/1/86
***************************************************************************
      implicit none
      integer nsped,IYR,ICM,ICD,length,iyer,j,isub,inum
      integer juday,iz,mt,ipick
      real jday,jday0,jday1,jbase_date,JULIAN,yearb,monthb,hourb
      real speed,spd,fff,vau,dayb,daym,grbs,grms,xyer,tm,gonl,tml
      real e,f,vuu,vou,daybb,daymm,grbss,grmss,cxx,oex,con,u,q,ui
      real s,xl,pm,pl,sl,ps,plm,skyn,vi,v,xi,vpp
      real vp,p,aul,aum,cra,cqa
      real aw,ai,ae,ae1,asp
      character*10 lname,label
      dimension label(180),vau(180),fff(180)
      dimension spd(180),vuu(180),vou(180)
      common/locat/tm,gonl
      common/costx/cxx(30),oex(5)
      common/fad/ipick
      common/vee/tml,con,u,q,ui
      common/boxa/s,xl,pm,pl,sl,ps,plm,skyn,vi,v,xi,vpp
      common/boxb/vp,p,aul,aum,cra,cqa
      common/boxs/aw,ai,ae,ae1,asp

*     ******************************************************************
*     *                                                                *
*     *   nsped = number of constituents to be calculated              *
*     *                                                                *
*     *                                                                *
*     *     ICM  = month of first data point                        *
*     *     ICD    = day  of first data point                         *
*     *     IYR   = year of first data point                         *
*     *     length  = length of series to be analyzed (in days)        *
*     *                                                                *
*     ******************************************************************
      yearb=IYR
      monthb=1.
      dayb=1.
      hourb=0.
      jbase_date=JULIAN(yearb,monthb,dayb,hourb)

      iyer =IYR
      dayb = 0.0
      daym = 183.0
      grbs = 0.0
      grms = 12.0

      if (mod(iyer,4) .eq. 0) then
        daym = 184.0
        grms = 0.0
      end if

      xyer = iyer

*     compute basic astronimical constants for the time period wanted
      call astro(xyer,dayb,daym,grbs,grms)

      tm = 0.0
      gonl = 0.0
      tml = 0.0

*     look up constituent parameters by matching the name
      do 100 j = 1,nsped
      lname = label(j)
      speed = 0.0
      call name (speed,lname,isub,inum,2)
      spd(j) = speed
      call vanduf (speed,e,f,2)
      fff(j) = f
      vou(j) = e
  100 continue

!      month = ICM
!      iday = ICD
!60      call CONCTJ(juday,month,iday,iyer)
60    yearb=IYR
      monthb=ICM
      dayb=ICD
      hourb=0.
      jday=JULIAN(yearb,monthb,dayb,hourb)
      juday=jday-jbase_date+1
      do 200 j = 1,nsped
      vuu(j) = vou(j) + spd(j)*(float(juday-1)*24.00)
  200 continue
      call twopi(vuu,nsped)
      daybb = juday
      daymm = juday
      grbss = 0.0
      grmss = 0.0
      daymm = juday + length/2

      call astro(xyer,daybb,daymm,grbss,grmss)

      do 277 iz = 1,nsped
      speed = spd(iz)
      call vanduf(speed,e,f,2)
      vau(iz) = e
  277 fff(iz) = f

*     round node to 3 decimals, vo+u to 1 place
      do 222 mt = 1,nsped
      fff(mt) = anint(fff(mt)*1000.) * 0.001
      vau(mt) = anint(vau(mt)*10.) * 0.1
  222 continue

      return
      endsubroutine equarg
      
************************************************************************
      subroutine astro(xyer,dayb,daym,grbs,grms)
************************************************************************      
      implicit none
      integer nyear,noe,ii,iii
      real xyer,dayb,daym,grbs,grms
      real pinv,xcen,xsx,xpx,xhx,xp1x,xnx,oex,t
      real xw,xi,aw,ai,ae,ae1,asp
      real cxx,doby,amit,ami,cplex,dicf,xcet,cdif,farm,farx
      real v,vi,vib,vb,vp,vpp,xib,vpb,vppb,xx,ang,anb,atb
      real cxsb,cxpb,cxhb,cxp1b,cig,cvx,cex,pvc,pvcp,an,at
      real pgx,xpg,raxe,raxn,rxx,rax,qxx,spxx,cqa
      real cra,um2,ul2,zes,zec,crav
      real pm,pl,sl,ps,plm,skyn,p,aul,aum,u,ui,q
      real tm,gonl,tml,con
      real s,xl
      common/locat/tm,gonl
      common/costx/cxx(30),oex(5)
      common/boxa/s,xl,pm,pl,sl,ps,plm,skyn,vi,v,xi,vpp
      common/boxb/vp,p,aul,aum,cra,cqa
      common/vee/tml,con,u,q,ui
      common/boxs/aw,ai,ae,ae1,asp
      common/boxxs/vib,vb,xib,vpb,vppb,cxsb,cxpb,cxhb,cxp1b

      pinv = 57.29578
      nyear = ifix(xyer)
      call orbit(xcen,xsx,xpx,xhx,xp1x,xnx,oex,t,xyer,5)
      xw = 23.4522944 - .0130125*t - .00000164*t**2 + .000000503*t**3
      xi = 5.14537628
      aw = xw*0.0174533
      ai = xi*0.0174533
      ae = 0.0548997
      ae1 = 0.01675104 - 0.0000418*t - .000000126*t**2
      asp = .46022931
      do 30 noe = 1,30
   30 cxx(noe) = 0.0
      if(dayb.gt.0.0) dayb = dayb - 1.0
      if(daym.gt.0.0) daym = daym - 1.0
      doby = 0.0
      amit = 0.0
      ami = xyer - xcen
      cplex = xcen/400.0 + 0.0001
      dicf = cplex - aint(cplex)
      if(ami.eq.0.0) go to 32
      xcet = xcen + 1.0
      cdif = xyer - xcet
      doby = cdif/4.0 + 0.0001
      amit = aint(doby)
      if(dicf.lt.0.001) amit = amit + 1.0
   32 farm = 0.25*ami
      farx = farm - aint(farm)
      cxx(1) = xsx + 129.384820*ami + 13.1763968*(dayb + amit) + 0.54901
     16532*grbs
      cxx(2) = xpx + 40.6624658*ami + 0.111404016*(dayb + amit) + 0.0046
     141834*grbs
      cxx(3) = xhx - 0.238724988*ami + 0.985647329*(dayb + amit) + 0.041
     1068639*grbs
      cxx(4) = xp1x + 0.01717836*ami + 0.000047064*(dayb + amit) + 0.000
     1001961*grbs
      cxx(5) = xpx + 40.6624658*ami + 0.111404016*(daym + amit) + 0.0046
     141834*grms
      cxx(6) = xnx - 19.3281858*ami - 0.052953934*(daym + amit) - 0.0022
     106414*grms
   40 cxx(7) = xpx + 40.6624658*ami + 0.111404016*(daym + amit) + 0.0046
     141834*grbs
      cxx(8) = xnx - 19.328185764*ami - 0.0529539336*(dayb + amit) - 0.0
     106414*grbs
      call twopi(cxx, 8)
   41 do 100 ii = 1,8
  100 cxx(ii) = float(ifix(cxx(ii)*100.0 + 0.5))*0.01
      ang = cxx(8)
      call table6(vib,vb,xib,vpb,vppb,xx,xx,xx,xx,xx,ang,anb,atb)
      cxx(26) = vib
      cxx(27) = vb
      cxx(28) = xib
      cxx(29) = vpb
      cxx(30) = vppb
      cxsb = cxx(1)
      cxpb = cxx(2)
      cxhb = cxx(3)
      cxp1b= cxx(4)
      ang = cxx(6)
      call table6(vi,v,xi,vp,vpp,cig,cvx,cex,pvc,pvcp,ang,an,at)
      cxx(9 ) = vi
      cxx(10) = v
      cxx(11) = xi
      cxx(12) = vp
      cxx(13) = vpp
  230 do 333 ii = 9,13
  333 cxx(ii) = float(ifix(cxx(ii)*100.0 + 0.5))*0.01
      pgx = cxx(5) - cxx(11)
      pgx = float(ifix(pgx*100.0 + 0.5))*0.01
      call twopi(pgx, 1)
      xpg = pgx*0.0174533
      cxx(14) = pgx
      raxe = sin(2.0*xpg)
      raxn = (cos(0.5*at)**2/(6.0*sin(0.5*at)**2)) - cos(2.0*xpg)
      rxx = 0.0
      if(raxe.eq.0.0.or.raxn.eq.0.0) go to 232
      rax = raxe/raxn
      if(rax.gt.3450.0) go to 232
        rxx   = atan(rax )*pinv
      cxx(22) = rxx
  232 cra = sqrt(1.0 - 12.0*(sin(0.5*at)**2/cos(0.5*at)**2)*cos(2.0*xpg)
     1 + 36.0*(sin(0.5*at)**4/cos(0.5*at)**4))
      um2 = 2.0*(cxx(11) - cxx(10))
      cxx(21) = um2
      cxx(24) = cra
      ul2 = um2 - rxx
  404 ul2 = ul2 + 180.0
  405 cxx(15) = ul2
      zes = (5.0*cos(aw) - 1.0)*sin(xpg)
      zec = (7.0*cos(aw) + 1.0)*cos(xpg)
      call fitan(zes,zec,qxx,spxx,2)
      cxx(23) = qxx
      crav = 0.5*um2 + qxx + 090.0
      cxx(16) = crav
      cqa = sqrt(0.25 + 1.5*((cos(aw)/cos(0.5*aw)**2)*cos(2.0*xpg)) + 2.
     125*(cos(aw)**2/cos(0.5*aw)**4))
      cxx(25) = cqa
      do 444 iii = 14,23
  444 cxx(iii) = float(ifix(cxx(iii)*100.0 + 0.5 ))*0.01
      pm   = cxx(1)
      pl   = cxx(2)
      sl   = cxx(3)
      ps   = cxx(4)
      plm  = cxx(5)
      skyn = cxx(6)
      vi   = cxx(9)
      v    = cxx(10)
      xi   = cxx(11)
      vp   = cxx(12)
      vpp  = cxx(13)
      p    = cxx(14)
      aul  = cxx(15)
      aum  = cxx(16)
      cra = cxx(24)
      cqa = cxx(25)
      u = v*0.0174533
      q = p*0.0174533
      ui = vi*0.0174533
      return
      endsubroutine astro
      
************************************************************************
      subroutine vanduf(speed,e,f,itype)
*     Order of constituents is same as in NAMES common
************************************************************************
      implicit none
      integer j,itype,ipick,ipik,ipikk,ipikkk,ipiikk
      real speed,e,f
      real tm,gonl
      real tml,con,u,q,ui
      real s,xl,pm,pl,sl,ps,plm,skyn,vi,v,xi,vpp
      real vp,p,aul,aum,cra,cqa
      real aw,ai,ae,ae1,asp
      real spd
      integer ms
      dimension spd(180),ms(180)
      common/locat/tm,gonl
      common/fad/ipick
      common/vee/tml,con,u,q,ui
      common/boxa/s,xl,pm,pl,sl,ps,plm,skyn,vi,v,xi,vpp
      common/boxb/vp,p,aul,aum,cra,cqa
      common/boxs/aw,ai,ae,ae1,asp
      common /speeds/spd
      common /mmss/ms

  500 format('0*** (v + u) not computed for constituent of speed',
     1 f12.7,'  ****',' Execution terminated')
      con = sl + tml
      do 600 j = 1,175
      ipick = j
      if(speed.eq.spd(j)) go to 611
  600 continue
      ipick = 176
  611 if(ipick.gt.164) go to 620
      if(ipick.gt.124) go to 618
      if(ipick.gt.88) go to 616
      if(ipick.gt.38) go to 614
      go to (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23
     1,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38),ipick
  614 ipik = ipick - 38
      go to (39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58
     1,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80
     2,81,82,83,84,85,86,87,88),ipik
  616 ipikk = ipick - 88
      go to (89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,10
     16,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,
     2123,124),ipikk
  618 ipikkk = ipick - 124
      go to (125,126,127,128,129,130,131,132,133,134,135,136,137,138,139
     1,140,141,142,143,144,145,146,147,148,149,150,151,152,153,154,155,
     2156,157,158,159,160,161,162,163,164),ipikkk
  620 ipiikk = ipick - 164
      go to (165,166,167,168,169,170,171,172,173,174,175,488),ipiikk

    1 e = 2.0*(con - pm + xi - v)
  201 f   = ((cos(0.5*aw)**4*cos(0.5*ai)**4)/(cos(0.5*ui)**4))
      go to 888

    2 e =  2.0*(con + xi - v) - 3.0*pm + pl
      go to 201

    3 e = 2.0*tml
  203 f = 1.0

      go to 888
    4 e = con - v - 2.0*(pm - xi) - 90.0
      go to 218

    5 e = con - vp + 90.0
  205 f = 1.0/sqrt(0.8965*sin(2.0*ui)**2 + 0.6001*sin(2.0*ui)*cos(u) + 0
     1.1006)
      go to 888

    6 e = 2.0*con - vpp
  206 f = 1.0/sqrt(19.0444*sin(ui)**4 + 2.7702*sin(ui)**2*cos(2.0*u) + 0
     1.0981)
      go to 888

    7 e = 2.0*con - pm - pl + aul
  207 f = ((cos(0.5*aw)**4*cos(0.5*ai)**4)/(cos(0.5*ui)**4))*(1.0/cra)
      go to 888

    8 e = 2.0*(con + xi - v + pl) - 4.0*pm
      go to 201

    9 e = sl - ps + 180.0 + 2.0*tml
      go to 203

   10 e = 2.0*tml - (sl - ps)
      go to 203

   11 e = 2.0*(con + xi - v - sl) - pm + pl + 180.0
      go to 201

   12 e = 2.0*(con + xi - v + sl) - 4.0*pm
      go to 201

   13 e = 2.0*(con + xi - v + sl) - 3.0*pm - pl
      go to 201

   14 e = con + pm - pl - v + 90.0
  214 f = (sin(2.0*aw)*(1.0 - 1.5*sin(ai)**2))/sin(2.0*ui)
      go to 888

   15 e = con - pm + aum
      f =  ((sin(aw)*cos(0.5*aw)**2*cos(0.5*ai)**4)/(sin(ui)*cos(0.5*ui)
     1**2))*(1.0/cqa)
      go to 888
   16 e = con - v + 2.0*(pm - xi)+ 90.0
  216 f = (sin(aw)*sin(.5*aw)**2*cos(.5*ai)**4)/(sin(ui)*sin(.5*ui)**2)
      go to 888
   17 e = tml + 270.0 - sl
      go to 203
   18 e = con - v - 3.0*pm + 2.0*xi + pl - 90.0
  218 f = (sin(aw)*cos(.5*aw)**2*cos(.5*ai)**4)/(sin(ui)*cos(.5*ui)**2)
      go to 888
   19 e = con - v - 4.0*pm + 2.0*xi + 2.0*pl - 90.0
      go to 218
   20 e = con - v - 3.0*pm + 2.0*xi - pl + 2.0*sl - 90.0
      go to 218
   21 e = 4.0*(con - pm + xi - v)
  221 f   = ((cos(0.5*aw)**4*cos(0.5*ai)**4)/(cos(0.5*ui)**4))**2
      go to 888
   22 e = 6.0*(con - pm + xi - v)
  222 f   = ((cos(0.5*aw)**4*cos(0.5*ai)**4)/(cos(0.5*ui)**4))**3
      go to 888
   23 e = 8.0*(con - pm + xi - v)
  223 f   = ((cos(0.5*aw)**4*cos(0.5*ai)**4)/(cos(0.5*ui)**4))**4
      go to 888
   24 e = 4.0*tml
      go to 203
   25 e = 6.0*tml
      go to 203
   26 e = 3.0*(con - pm + xi - v) + 180.0
  226 f = (cos(0.5*aw )**6*cos(0.5*ai)**6)/cos(0.5*ui)**6
      go to 888
   27 e = tml + 180.0
      go to 203
   28 e = 2.0*(con - pm + xi - v) + ( con - vp + 90.0)
  228 f = ((cos(0.5*aw)**4*cos(0.5*ai)**4)/(cos(0.5*ui)**4))**1*(1./sqrt
     1(0.8965*sin(2.0*ui)**2 + 0.6001*sin(2.0*ui)*cos(u) + 0.1006))
      go to 888
   29 e = 4.0*(con - pm + xi - v) - (con - vp + 90.0)
  229 f = ((cos(0.5*aw)**4*cos(0.5*ai)**4)/(cos(0.5*ui)**4))**2*(1./sqrt
     1(0.8965*sin(2.0*ui)**2 + 0.6001*sin(2.0*ui)*cos(u) + 0.1006))
      go to 888
   30 e = 4.0*(con + xi - v) + pl - 5.0*pm
      go to 221
   31 e = 2.0*(con - pm + xi - v) + 2.0*tml
      go to 201
   32 e = 4.0*tml - 2.0*(con - pm + xi - v)
      go to 201
   33 e = 2.0*(pm - xi)
  233 f = (sin(aw)**2*cos(0.5*ai)**4)/sin(ui)**2
      go to 888
   34 e = 2.0*tml - 2.0*(con - pm + xi - v)
      go to 201
   35 e = pm - pl
  235 f = ((2./3.- sin(aw)**2)*(1.- 1.5*sin(ai)**2))/(2./3.- sin(ui)**2)
      go to 888
   36 e = sl
      go to 203
   37 e = 2.0*sl
      go to 203
   38 e = 2.0*(con - v - 2.0*(pm - xi) - 90.0) - (tml+270.0 - sl)
  238 f = ((sin(aw)*cos(0.5*aw)**2*cos(0.5*ai)**4)/(sin(ui)*cos(0.5*ui)*
     1*2))**2
      go to 888
   39 e = 2.0*(con - pm + xi - v) - (tml + 270.0 - sl)
      go to 201
   40 e = con + 2.0*sl - v - pm - pl + 90.0
      go to 258
   41 e = 2.0*(tml + 270.0 - sl) - (con - v - 2.0*(pm - xi) - 90.0)
      go to 218
   42 e = 2.0*tml - (con - v - 2.0*(pm - xi) - 90.0)
      go to 218
   43 e = 2.0*tml + pm - pl
      go to 221
   44 e = 4.0*(con + xi - v) - 5.0*pm - 2.0*tml  + pl
      go to 221
   45 e = con - v - 2.0*(pm - xi) + tml - sl + 180.0
      go to 218
   46 e = 4.0*con - 2.0*(pm - xi + v + tml) - vpp
      go to 259
   47 e = 4.0*(con + xi - v) - 6.0*pm + 2.0*(pl - tml)
      go to 221
   48 e = 6.0*(con - pm) + 4.0*(xi - v - tml) + aul
      go to 261
   49 e = 6.0*con - 5.0*pm + 4.0*(xi - v - tml) - pl + aul
      go to 261
   50 e = 2.0*(tml + pm - xi + v) - vpp
      go to 259
   51 e = 4.0*(xi - pm - v) + 2.0*(tml + vpp)
  251 f = ((cos(.5*aw)**4*cos(.5*ai)**4)/(cos(.5*ui)**4))**2*(1./sqrt(19
     1.0444*sin(ui)**4 + 2.7702*sin(ui)**2*cos(2.*u) + .0981)**2)
      go to 888
   52 e = 6.0*con - 4.0*tml - 3.0*pm + 2.0*(xi - v) - pl - vpp + aul
  252 cra = sqrt(1.0 - 12.0*(sin(.5*ui)**2/cos(.5*ui)**2)*cos(2.0*q) + 3
     16.0*(sin(.5*ui)**4/cos(.5*ui)**4))
      f = ((cos(.5*aw)**4*cos(.5*ai)**4)/(cos(.5*ui)**4))**2*(1./sqrt(19
     1.0444*sin(ui)**4 + 2.7702*sin(ui)**2*cos(2.*u) + .0981))*(1./cra)
      go to 888
   53 e = 6.0*con - 4.0*tml + 2.0*(xi - v - pm - vpp)
  253 f = ((cos(.5*aw)**4*cos(.5*ai)**4)/(cos(.5*ui)**4))**1*(1./sqrt(19
     1.0444*sin(ui)**4 + 2.7702*sin(ui)**2*cos(2.*u) + .0981)**2)
      go to 888
   54 e = 4.0*tml - 2.0*con + pl - pm + vpp
  254 f = ((cos(.5*aw)**4*cos(.5*ai)**4)/(cos(.5*ui)**4))**2*(1./sqrt(19
     1.0444*sin(ui)**4 + 2.7702*sin(ui)**2*cos(2.*u) + .0981))
      go to 888
   55 e = 4.0*con - 2.0*(tml + vpp) + pm - pl
      go to 251
   56 e = 2.0*tml + con - v - 2.0*(pm - xi) - 90.0
      go to 218
   57 e = 2.0*tml + con - vp + 90.0
      go to 205
   58 e = 3.0*(con - v) + 4.0*xi - 5.0*pm + pl - 90.0
  258 f = ((cos(.5*aw)**4*cos(.5*ai)**4)/(cos(.5*ui)**4))**1*((sin(aw)*c
     1os(.5*aw)**2*cos(.5*ai)**4)/(sin(ui)*cos(.5*ui)**2))
      go to 888
   59 e = 4.0*con - 2.0*(pm - xi + v) - vpp
  259 f = ((cos(.5*aw)**4*cos(.5*ai)**4)/(cos(.5*ui)**4))**1*(1./sqrt(19
     1.0444*sin(ui)**4 + 2.7702*sin(ui)**2*cos(2.*u) + .0981))
      go to 888
   60 e = 2.0*(tml + con + xi - v) - 3.0*pm + pl
      go to 201
   61 e = 6.0*con - 5.0*pm + 4.0*(xi - v) - 2.0*tml - pl + aul
  261 f = ((cos(.5*aw)**4*cos(.5*ai)**4)/(cos(.5*ui)**4))**3*(1./sqrt(1.
     10 - 12.0*(sin(.5*ui)**2/cos(.5*ui)**2)*cos(2.0*q) + 36.0*(sin(.5*u
     2i)**4/cos(.5*ui)**4)))
      go to 888
   62 e = 6.0*(con - pm + xi - v) - 2.0*tml
      go to 222
   63 e = 4.0*con - 3.0*pm + 2.0*(xi - v) - pl + aul
  263 f = ((cos(.5*aw)**4*cos(.5*ai)**4)/(cos(.5*ui)**4))**2*(1./sqrt(1.
     10 - 12.0*(sin(.5*ui)**2/cos(.5*ui)**2)*cos(2.0*q) + 36.0*(sin(.5*u
     2i)**4/cos(.5*ui)**4)))
      go to 888
   64 e = 4.0*(con + xi - v) - 6.0*pm + 2.0*pl
      go to 221
   65 e = 2.0*(con + tml) - pm - pl + aul
      go to 207
   66 e = 5.0*(con - v) - 7.0*pm + 6.0*xi + pl - 90.0
  266 f = ((cos(.5*aw)**4*cos(.5*ai)**4)/(cos(.5*ui)**4))**2*((sin(aw)*c
     1os(.5*aw)**2*cos(.5*ai)**4)/(sin(ui)*cos(.5*ui)**2))
      go to 888
   67 e = 5.0*(con - v) + 6.0*(xi - pm) - 90.0
      go to 266
   68 e = 5.0*con + 4.0*(xi - pm - v) - vp + 90.0
      go to 229
   69 e = 3.0*con + 2.0*(xi - pm - v + tml) - vp + 90.0
      go to 228
   70 e = 5.0*con + 2.0*(xi - pm - v) - (vp + vpp) + 90.0
  270 f = ((cos(.5*aw)**4*cos(.5*ai)**4)/(cos(.5*ui)**4))**1*(1./sqrt(19
     1.0444*sin(ui)**4 + 2.7702*sin(ui)**2*cos(2.*u) + .0981))*(1./sqrt(
     2.8965*sin(2.*ui)**2 + .6001*sin(2.*ui)*cos(u) + 0.1006))
      go to 888
   71 e = 4.0*(con - pm + xi - v) + tml - sl + 270.0
      go to 221
   72 e = 6.0*(con - pm + xi - v) - tml + sl - 270.0
      go to 222
   73 e = 5.0*(con - pm) + 4.0*(xi - v) + pl - vp + 90.0
      go to 229
   74 e = 2.0*(con - pm + xi - v) + 4.0*tml
      go to 201
   75 e = 6.0*(con + xi - v) - 7.0*pm + pl
      go to 222
   76 e = 4.0*(con + xi - v) - 5.0*pm + 2.0*tml + pl
      go to 221
   77 e = 4.0*(con - pm + xi - v) + 2.0*tml
      go to 221
   78 e = 8.0*con - 9.0*pm + 6.0*(xi - v) - 2.0*tml + pl + aul
  278 f = ((cos(.5*aw)**4*cos(.5*ai)**4)/(cos(.5*ui)**4))**4*(1./sqrt(1.
     10 - 12.0*(sin(.5*ui)**2/cos(.5*ui)**2)*cos(2.0*q) + 36.0*(sin(.5*u
     2i)**4/cos(.5*ui)**4)))
      go to 888
   79 e = 6.0*(con + xi - v) - 8.0*pm + 2.0*pl
      go to 222
   80 e = 4.0*con - 3.0*pm + 2.0*(xi - v + tml) - pl + aul
      go to 263
   81 e = 6.0*con - 5.0*pm + 4.0*(xi - v) - pl + aul
      go to 261
   82 e = 4.0*con + 2.0*(xi - pm - v + tml) - vpp
      go to 259
   83 e = 8.0*(con - pm) + 6.0*(xi - v) - 2.0*tml + aul
      go to 278
   84 e = 8.0*con - 7.0*pm + 6.0*(xi - v) - 2.0*tml - pl + aul
      go to 278
   85 e = 6.0*con + 4.0*(xi - pm - v) - vpp
      go to 254
   86 e = 7.0*(con - v) - 9.0*pm + 8.0*xi + pl - 90.0
  286 f = ((cos(.5*aw)**4*cos(.5*ai)**4)/(cos(.5*ui)**4))**3*((sin(aw)*c
     1os(.5*aw)**2*cos(.5*ai)**4)/(sin(ui)*cos(.5*ui)**2))
      go to 888
   87 e = 7.0*con + 6.0*(xi - v) - 8.0*pm + 2.0*pl - vp + 90.0
  287 f = ((cos(0.5*aw)**4*cos(0.5*ai)**4)/(cos(0.5*ui)**4))**3*(1./sqrt
     1(0.8965*sin(2.0*ui)**2 + 0.6001*sin(2.0*ui)*cos(u) + 0.1006))
      go to 888
   88 e = 5.0*(con - v) + 6.0*(xi - pm) + 2.0*tml - 90.0
      go to 266
   89 e = 5.0*con + 4.0*(xi - pm) - 3.0*v + 2.0*tml - vpp - 90.0
  289 f = ((cos(.5*aw)**4*cos(.5*ai)**4)/(cos(.5*ui)**4))**1*(1./sqrt(19
     1.0444*sin(ui)**4 + 2.7702*sin(ui)**2*cos(2.*u) + .0981))*( (sin(aw
     2)*cos(.5*aw)**2*cos(.5*ai)**4)/(sin(ui)*cos(.5*ui)**2))
      go to 888
   90 e = 6.0*con - 7.0*pm + 6.0*(xi - v) + 2.0*tml + pl
      go to 222
   91 e = 6.0*(con - pm + xi - v) + 2.0*tml
      go to 222
   92 e = 4.0*(con - pm + xi - v) + 4.0*tml
      go to 221
   93 e = 8.0*(con + xi - v) - 10.0*pm + 2.0*pl
      go to 223
   94 e = 8.0*(con + xi - v) - 9.0*pm + pl
      go to 223
   95 e = 6.0*con - 5.0*pm + 4.0*(xi - v) + 2.0*tml - pl + aul
      go to 261
   96 e = 10.0*con - 9.0*pm + 8.0*(xi - v) - 2.0*tml - pl + aul
  296 f = ((cos(.5*aw)**4*cos(.5*ai)**4)/(cos(.5*ui)**4))**5*(1./sqrt(1.
     10 - 12.0*(sin(.5*ui)**2/cos(.5*ui)**2)*cos(2.0*q) + 36.0*(sin(.5*u
     2i)**4/cos(.5*ui)**4)))
      go to 888
   97 e = 8.0*con - 7.0*pm + 6.0*(xi - v) - pl + aul
      go to 278
   98 e = 8.0*con + 6.0*(xi - pm - v) - vpp
  298 f = ((cos(.5*aw)**4*cos(.5*ai)**4)/(cos(.5*ui)**4))**3*(1./sqrt(19
     1.0444*sin(ui)**4 + 2.7702*sin(ui)**2*cos(2.*u) + .0981))
      go to 888
   99 e = 6.0*con + 4.0*(xi - pm - v) + 2.0*tml - vpp
      go to 254
  100 e = 9.0*con - 10.0*pm + 8.0*(xi - v) + 2.0*pl - vp + 90.0
  300 f = ((cos(0.5*aw)**4*cos(0.5*ai)**4)/(cos(0.5*ui)**4))**4*(1./sqrt
     1(0.8965*sin(2.0*ui)**2 + 0.6001*sin(2.0*ui)*cos(u) + 0.1006))
      go to 888
  101 e = 9.0*(con - pm) + 8.0*(xi - v) + pl - vp + 90.0
      go to 300
  102 e = 9.0*con + 8.0*(xi - pm - v) - vp + 90.0
      go to 300
  103 e = 7.0*con + 6.0*(xi - pm - v) + 2.0*tml - vp + 90.0
      go to 287
  104 e = 10.0*(con + xi - v) - 11.0*pm + pl
  304 f = ((cos(.5*aw)**4*cos(.5*ai)**4)/(cos(.5*ui)**4))**5
      go to 888
  105 e = 10.0*(con - pm + xi - v)
      go to 304
  106 e = 8.0*(con + xi - v) - 9.0*pm + pl + 2.0*tml
      go to 223
  107 e = 8.0*(con - pm + xi - v) + 2.0*tml
      go to 223
  108 e = 8.0*con - 7.0*pm + 6.0*(xi - v) - pl + aul
      go to 278
  109 e = 6.0*(con - pm + xi - v) + 4.0*tml
      go to 222
  110 e = 9.0*con + 8.0*(xi - pm - v) + 2.0*tml - vp + 90.0
      go to 300
  111 e = 10.0*(con + xi - v) - 11.0*pm + pl
      go to 304
  112 e = 10.0*(con - pm + xi - v) + 2.0*tml
      go to 304
  113 e = 10.0*con - 9.0*pm + 8.0*(xi - v) + 2.0*tml - pl + aul
      go to 296
  114 e = 8.0*(con - pm + xi - v) + 4.0*tml
      go to 223
  115 e = tml - 2.0*sl + ps - 90.0
      go to 203
  116 e = con + sl - ps + 90.0
      go to 203
  117 e = con + 2.0*sl + 90.0
      go to 203
  118 e = tml + pm - sl + pl - v + 90.0
      go to 258
  119 e = 2.0*con + pm - v - vp - pl + 180.0
  319 f =  ((sin(2.*aw)*(1.0 - 1.5*sin(ai)**2))/sin(2.0*ui))*(1./sqrt(.8
     1965*sin(2.*ui)**2 + .6011*sin(2.*ui)*cos(u) + 0.1006))
      go to 888
  120 e = 2.0*(con - v) - 5.0*pm + 4.0*xi + pl + 180.0
      go to 238
  121 e = 3.0*(con - v) - 4.0*(pm - xi) - 90.0
      go to 258
  122 e = 2.0*(con + tml) - vpp
      go to 206
  123 e = con + 2.0*(pm - xi - vp) + v + 270.0
  323 f = ((sin(aw)*cos(.5*aw)**2*cos(.5*ai)**4)/(sin(ui)*cos(.5*ui)**2)
     1)**1*(1.0/sqrt(0.8965*sin(2.0*ui)**2 + 0.6001*sin(2.0*ui)*cos(u) +
     2 0.1006))**2
      go to 888
  124 e = con - 2.0*v - 4.0*(pm - xi) + vp - 270.0
  324 f = ((sin(aw)*cos(.5*aw)**2*cos(.5*ai)**4)/(sin(ui)*cos(.5*ui)**2)
     1)**2*(1.0/sqrt(0.8965*sin(2.0*ui)**2 + 0.6001*sin(2.0*ui)*cos(u) +
     2 0.1006))
      go to 888
  125 e = 6.0*(con - pm) + 4.0*(xi - v - tml) + 2.0*pl - vpp
  325 go to 254
  126 e = 6.0*con - 5.0*pm + 4.0*(xi - v - tml) + pl - vpp
  326 go to 254
  127 e = 6.0*con - 4.0*tml - 3.0*pm + 2.0*(xi - v - vpp) + pl
  327 go to 253
  128 e = 6.0*con - 5.0*pm + 4.0*(xi - v) - 2.0*tml + pl - vpp
  328 go to 254
  129 e = 4.0*con - 3.0*pm + 2.0*(xi - v) + pl - vpp
  329 go to 259
  130 e = 8.0*con - 9.0*pm + 6.0*(xi - v) + 3.0*pl - 2.0*tml - vpp
  330 go to 298
  131 e = 8.0*(con - pm) + 6.0*(xi - v) + 2.0*(pl - tml) - vpp
  331 go to 298
  132 e = 8.0*con - 7.0*pm + 6.0*(xi - v) - 2.0*tml + pl - vpp
  332 go to 298
  133 e = 6.0*con - 5.0*pm + 4.0*(xi - v) + pl - vpp
  333 go to 254
  134 e = 4.0*con - 3.0*pm + 2.0*(xi - v + tml) + pl - vpp
  334 go to 259
  135 e = 10.0*con - 9.0*pm + 8.0*(xi - v) - 2.0*tml + pl - vpp
  335 f = ((cos(.5*aw)**4*cos(.5*ai)**4)/(cos(.5*ui)**4))**4*(1./sqrt(19
     1.0444*sin(ui)**4 + 2.7702*sin(ui)**2*cos(2.*u) + .0981))
      go to 888
  136 e = 8.0*con - 7.0*pm + 6.0*(xi - v) + pl - vpp
  336 go to 298
  137 e = 6.0*con - 5.0*pm + 4.0*(xi - v) - 2.0*tml + pl - vpp
  337 go to 254
  138 e = 8.0*con - 7.0*pm + 6.0*(xi - v) + 2.0*tml + pl - vpp
  338 go to 298
  139 e = 10.0*con - 9.0*pm + 8.0*(xi - v) + 2.0*tml + pl - vpp
  339 go to 335
  140 e = 4.0*(con + xi - v) - 6.0*pm + 2.0*pl + sl - tml - 270.0
  340 go to 221
  141 e = 3.0*sl - 2.0*vp - tml - 90.0
  341 f = 1.0/(.8965*sin(2.0*ui)**2 + .6001*sin(2.0*ui)*cos(u) + .1006)
      go to 888
  142 e = tml + vp - 3.0*sl + 90.0
  342 go to 205
  143 e = 2.0*tml - vp
  343 go to 205
  144 e = 2.0*(tml - sl) + vpp
  344 go to 206
  145 e = 2.0*(tml - vpp) + 4.0*sl
  345 f = (1.0/sqrt(19.0444*sin(ui)**4 + 2.7702*sin(ui)**2*cos(2.0*u) +
     10.0981))**2
      go to 888
  146 e = 2.0*tml - 2.0*(sl - ps)
  346 go to 203
  147 e = 4.0*tml + ps - sl
  347 go to 203
  148 e = 4.0*tml - 2.0*sl + vpp
  348 go to 206
  149 e = 4.0*tml + 6.0*sl - 3.0*vpp
  349 f = (1.0/sqrt(19.0444*sin(ui)**4 + 2.7702*sin(ui)**2*cos(2.0*u) +
     10.0981))**3
      go to 888
  150 e = 4.0*tml - 3.0*(sl - ps)
  350 go to 203
  151 e = con - v - 2.0*(pm - xi) + tml + 90.0
  351 go to 218
  152 e = 2.0*(pm - xi) + v + tml - con - 90.0
  352 go to 218
  153 e = pm - xi - 90.0
  353 f = 0.31920/(sin(ui) - sin(ui)**3)
      go to 888
  154 e = pm - sl
  354 go to 235
  155 e = sl - ps
  355 go to 203
  156 e = 3.0*sl
  356 f = 1.0
      go to 888
  157 e = 4.0*sl
  357 go to 356
  158 e = 6.0*sl
  358 go to 356
  159 e = 8.0*sl
  359 go to 356
  160 e = 10.0*sl
  360 go to 356
  161 e = 12.0*sl
  361 go to 356
  162 e = 24.0*sl
  362 go to 356
  163 e = 3.0*tml + 180.0
  363 go to 356
  164 e = 5.0*tml + 180.0
  364 go to 356
  165 e = 2.0*(con - v) - 4.0*(pm - xi) - 180.0
  365 f =((sin(aw)*cos(.5*aw)**2*cos(.5*ai)**4)/(sin(ui)*cos(.5*ui)**2))
     1**2
      go to 888
  166 e = tml + con - vp + 270.0
  366 go to 205
  167 e = 3.0*(con - pm) + 2.0*(xi - v) + pl - vp + 90.0
  367 f = ((cos(0.5*aw)**4*cos(0.5*ai)**4)/cos(0.5*ui)**4 )*(1.0/sqrt(0.
     18965*sin(2.0*ui)**2 + 0.6001*sin(2.0*ui)*cos(u) + 0.1006))
      go to 888
  168 e = 3.0*tml + 270.0 - sl
  368 go to 203
  169 e = 3.0*con - vpp - vp + 90.0
  369 f = 1.0/(sqrt(19.0444*sin(ui)**4 + 2.7702*sin(ui)**2*cos(2.0*u) +
     10.0981)*sqrt(0.8965*sin(2.0*ui)**2 + 0.6001*sin(2.0*ui)*cos(u) +
     20.1006))
      go to 888
  170 e = 4.0*(con - v) + 6.0*xi - 7.0*pm + pl - 180.0
  370 f = ((cos(0.5*aw)**4*cos(0.5*ai)**4)/(cos(0.5*ui)**4))*((sin(aw)*c
     1os(0.5*aw)**2*cos(.5*ai)**4)/(sin(ui)*cos(.5*ui)**2))**2
      go to 888
  171 e = 4.0*(con - v) - 6.0*(pm - xi) - 180.0
  371 go to 370
  172 e = 2.0*(tml + con - v) - 4.0*(pm - xi) - 180.0
  372 go to 365
  173 e = 7.0*tml + 180.0
  373 go to 356
  174 e = 8.0*tml
  374 go to 356
  175 e = 10.0*tml
  375 go to 356
  488 print 500, speed
      stop
  888 go to (400,390),itype
  390 e = e + float(ms(ipick))*gonl - spd(ipick)*(tm/15.0)
  400 call twopi(e,1)
      f = 1.0/f
      return
      endsubroutine vanduf
      
************************************************************************
      subroutine name(spdd,itag,isub,inum,icode)
c     this subroutine identifies the constituent by its speed,
*     name label or constituent number,
c     and makes it available for labeling.
*     ICODE = 1, by speed
*     ICODE = 2, by label
*     ICODE = 3, by number
c     it also determine the subscript of the constituent
c        order of constituent speeds***  m(2),n(2),s(2),o(1),k(1),k(2)
c      l(2),2n(2)r(2),t(2),lambda(2),mu(2),nu(2),j(1),m(1),oo(1),p(1)
c      q(1),2q(1),rho(1),m(4),m(6),m(8),s(4),s(6),m(3),s(1),mk(3),2mk(3)
c      mn(4),ms(4),2sm(2),mf,msf,mm,sa,ssa
      implicit none
      integer :: isub,inum,icode,i,j,ip
      real spd,spdd
      common /mmss/ip
      common /speeds/spd
      common /names/lable
      character*10 lable(180)*10,itag
      dimension spd(180),ip(180)

    1 format(10x,'Constituent of speed ',f12.7,' not in list.')
    2 format(10x,'Constituent ', a10,' not in the list.')
    3 format(10x,'Constituent no.',i4,' not in list.')

      if (icode.eq.1) then
*       search by speed
        do 100 j = 1,175
        if(spdd.ne.spd(j)) go to 100
        itag = lable(j)
        isub = ip(j)
        inum = j
        return
100     continue
        print 1,spdd
      else if (icode.eq.2) then
*       search by name
        do 200 i = 1,175
        if(itag.ne.lable(i)) go to 200
        spdd = spd(i)
        isub = ip(i)
        inum = i
        return
200     continue
        print 2, itag
      else if (icode.eq.3) then
*       search by number
!        if (i.gt.0.and.i.le.175) then
!          itag = lable(k)
!          spdd = spd(k)
!          isub = ip(k)
!          return
!        end if
        write(*,*)'get wrong icode, icode must be 1 or 2'
        stop
        print 3, inum
      end if

      stop '**** Execution terminated in NAME (illegal icode) ****'

      endsubroutine name
      
************************************************************************
      subroutine orbit(xcen,xsx,xpx,xhx,xp1x,xnx,oex,t,xyer,nnn)
************************************************************************      
      implicit none
      integer nnn,i,jk
      real xcen,xsx,xpx,xhx,xp1x,xnx,oex(nnn),t,xyer
      real s,p,xh,p1,xn,xcan,yr,gat,gp,col
      
      s = 13.1763968
      p = 0.1114040
      xh = 0.9856473
      p1 = 0.0000471
      xn = -.0529539
      xcan = xyer*0.01 + 0.001
      xcen = aint(xcan)*100.0
      t = -3.0
      yr = 2.5
      gat = 1600.0
      do 10 jk = 1,30
      gp = gat/400.0 + 0.00001
      col = gp - aint(gp)
      if(col.lt.0.010) go to 11
      if(gat.eq.xcen) go to 12
      yr = yr - 1.0
      go to 9
   11 if( gat.eq.xcen) go to 12
    9 gat = gat + 100.0
   10 continue
   12 t = (gat - 1900.0)*0.01
      oex(1) = 270.437422 + 307.892*t + 0.002525*t**2 + .00000189*t**3 +
     1 yr*s
      oex(2) = 334.328019 + 109.032206*t - 0.01034444*t**2 - .0000125*t*
     1*3 + yr*p
      oex(3) = 279.696678 + 0.768925*t + .0003205*t**2 + yr*xh
      oex(4) = 281.220833 + 1.719175*t + 0.0004528*t**2 + .00000333*t**3
     1 + yr*p1
      oex(5) = 259.182533 - 134.142397*t + .00210556*t**2 + .00000222*t*
     1*3 + yr*xn
      call twopi(oex,5)
      do 100 i = 1,5
  100 oex(i) = float(ifix(oex(i)*100.0 + 0.5))*0.01
      xsx = oex(1)
      xpx = oex(2)
      xhx = oex(3)
      xp1x = oex(4)
      xnx = oex(5)
      
      return
      endsubroutine orbit
      
************************************************************************
      subroutine table6(vi,v,xi,vp,vpp,cig,cvx,cex,pvc,pvcp,ang,an,at)
************************************************************************     
      implicit none
      real vi,v,xi,vp,vpp,cig,cvx,cex,pvc,pvcp,ang,an,at      
      real ax,eye
      real aw,ai,ae,ae1,asp
      real c9,c10,c11,vxx,vxxe,vxxn,term,exx,ezz,exez
      real a22,b22,vpx,vpxe,vpxn
      real a47,b47,vpye,vpyn,vpy
      common/boxs/aw,ai,ae,ae1,asp
      
      v = 0.0
      xi = 0.0
      vp = 0.0
      vpp = 0.0
      an = ang*0.0174533
      ax = ang
      eye = cos(ai)*cos(aw) - sin(ai)*sin(aw)*cos(an)
      c9 = acos(eye)*57.2957795
      vi = float(ifix(c9*100.0 + 0.5))*0.01
      cig = vi*0.0174533
      at = cig
      if(cig.eq.0.0) go to 230
      if(ax.eq.0.0.or.ax.eq.180.0) go to 230
      vxxe = sin(ai)*sin(an)
      vxxn = cos(ai)*sin(aw) + sin(ai)*cos(aw)*cos(an)
      if(vxxe.eq.0.0.or.vxxn.eq.0.0) go to 201
      vxx = vxxe/vxxn
      c10 = atan(vxx)*57.2957795
      v = float(ifix(c10*100.0 + 0.5))*0.01
      if(ax.gt.180.0.and.v.gt.0.0) v = -1.0*v
  201 cvx = v*0.0174533
      term = sin(ai)*(cos(aw)/sin(aw))
      exx = term*(sin(an)/cos(an)) + (cos(ai) - 1.0)*sin(an)
      if(exx.eq.0.0) go to 202
      ezz = term + cos(ai)*cos(an) + (sin(an)**2/cos(an))
      if(ezz.eq.0.0) go to 202
      exez = exx/ezz
      if(exez.gt.3450.0) go to 202
      c11 = atan(exez)*57.2957795
      xi = float(ifix(c11*100.0 + 0.5))*0.01
      if(ax.gt.180.0.and.xi.gt.0.0) xi = -1.0*xi
  202 cex = xi*0.0174533
      a22 = (0.5 + 0.75*ae**2)*sin(2.0*cig)
      b22 = (0.5 + 0.75*ae1**2)*sin(2.0*aw)*asp
      vpxe = a22*sin(cvx)
      vpxn = a22*cos(cvx) + b22
      if(vpxe.eq.0.0.or.vpxn.eq.0.0) go to 203
      vpx = vpxe/vpxn
      if(vpx.gt.3450.0) go to 203
      vp = atan(vpx)*57.2957795
      if(ax.gt.180.0.and.vp.gt.0.0) vp = -1.0*vp
  203 pvc = vp*0.0174533
      a47 = (0.5 + 0.75*ae**2)*sin(cig)**2
      b47 = (0.5 + 0.75*ae1**2)*asp*sin(aw)**2
      vpye = a47*sin(2.0*cvx)
      vpyn = a47*cos(2.0*cvx) + b47
      if(vpye.eq.0.0.or.vpyn.eq.0.0) go to 204
      vpy = vpye/vpyn
      if(vpy.gt.3450.0) go to 204
      vpp = atan(vpy)*57.2957795
      if(ax.gt.180.0.and.vpp.gt.0.0) vpp = -1.0*vpp
  204 pvcp = vpp*0.0174533
  230 return
      endsubroutine table6
      
C*************************************************************************
      subroutine gregorian(jday,yr,month,day,hour)
C AJ 09/26/08     include './library/gregorian.f'      
C*************************************************************************      
      implicit none
      real jday,yr,month,day,hour,dayoweek,week
      real a,b,c,d,e

      a = aint(jday+.5)
      b = a+1537
      c = aint((b-122.1)/365.25)
      d = aint(365.25*c)
      e = aint((b-d)/30.6001)
      day = b-d-aint(30.6001*e)+mod(jday+0.5d00,1.0d00)
      month = e-1-12*aint(e/14)
      yr = c-4715-aint((7+month)/10)
      hour = mod(jday+0.5d00,1.0d00)*24.
      day = aint(day)

!      dayoweek = mod(aint(jday+0.5),7)
!      week = aint((jday-2444244.5)/7)
!  2000 leap year crashed.
      if(month.eq.2. .and. day .eq. 31. ) then
       yr=yr-1.
       day = 29.
      endif
      
      endsubroutine gregorian
      
C************************************************************************* 
      function julian(yr,monthb,dayb,hourb)
C AJ 09/26/08     include './library/julian.f'      
C*************************************************************************       
      implicit none
      real yr,monthb,dayb,hourb
      real Yb,mb,julian
      
      Yb=yr
      if (yr.lt.100. .and. yr.gt.50.) then
        Yb=yr+1900.
      elseif (yr.lt.100. .and. yr.le.50.) then
        Yb=yr+2000.
      endif

      if (monthb.le. 2.) then
        Yb=yb-1
        mb=monthb+12.
      else
        Yb=Yb
        mb=monthb
      endif

      JULIAN=aint(365.25*Yb)+aint(30.6001*(mb+1))
     1    +dayb+hourb/24.+1720981.50
!      print *,'test1'  
!      print *,aint(365.25*Y)
!      print *,aint(30.6001*(m+1))
!      print *,day
!      print *,hour/24.
!      print *,julian
!      print *,'test2'
C  1980 epic....
c      JULIAN = aint(365.25*(Y-1980.)) + aint(30.6001*(m+1))
c     &        + day + hour/24.
c      write(*,'(''Julian ='',f12.4,4f8.2)') julian, yr,month,day,hour
c      END FUNCTION JULIAN
       return
       endfunction julian
      
      endmodule nos_tide