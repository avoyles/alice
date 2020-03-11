C
c alice2017.f  changes to 14 June 2017
c   test moving evap calc of barf earlier in sub 3/6/2016
c this is a version in development;please contact hmblann@gmail.com with problems
c a collaborative effort of M.Blann,J.Bisplinghoff(Bonn Uni.),W.Scobel(Uni.Hamburg),
c A.Konobeyev(IfK,Karlsruhe),W.B.Wilson and S.G.Mashnik(LANL),K.K.Gudima(IAP,Moldava)
c A.V. Ignatyuk,V.P. Luneev, Yu. Shubin(IPPE,Obninsk).We acknowledge
c contributions as well from D. Madland, Marie Giaccri, and others who have
c generously contributed to the development of this code over the past half century.
c This version differs from earlier versions in allowing isotopic targets to be used
c It contains cluster exit channels d,t,3He,7Be on demand,in addition to default n,p,4He.
c An ENDF output for 1,2,3 n,p,alpha out reactions is an option.The logic used could be
c extended to include other clusters and to higher multiplicities if required.
c Earlier versions permitted setting of coincident detector 'gates' to predict
c spectra such as might be measured in an experiment with multiple detectors.Parts of
c this coding is still in place, as are arrays of recoil sdcs and ddcs.
c For further info on thoese options, contact hmblann@gmail.com
c
c Updated Fermi decay subroutines, and ranf()==>rndm(-1.)) March 14,2010 by K.K. Gudima
c
c compilation notes:
c
c it is recommended to compile and run with the bounds check option initially:
c
c   f95 -C -o hms.x hms.f or f95 -fbounds-check -o hms.x hms.f
c  or on a 64 bit computer:
c  f95 alice2016.f -fdefault-real-8 -fdefault-double-8 -fbounds-check -o alice2016.x
c
cc After verifying that there are no arrays out of bounds for several runs,
c much time in execution can be saved by recompiling without bounds check
c and with a low level of optimization. In the past, using a high level of
c optimization has led to execution with erroneous results. It is therefore
c strongly recommended that if optimization is used, several sets of results
c for the identical input file be compared with a non-optimized compile version!
c attempts have been made to assure that this code will compile on
c different FORTRAN compilers....Good luck!
c
c Marie Giaccri table of GDR photonuclear parms from S.Goriely,M.Pearson,F.Tondeur;
c file is 'gdr-parameters-theor'.
c
c  fission width ;
c these will be Sierk folded Yukawa+exp model results,except for
c actinides which are read in from table/file tmadland,based on
c private communication,nov.1,2002 from D.Madland.
c
c May 2002 add fermi breakup routines from K.K.Gudima to handle decay of
c fragments  lighter than 'aferm',a parameter presently an input option
c with default of A=12.
c  dec 2000 begin putting hi angle/exciton coupling into HI pe routines
 
c  RECLSV:
c      stores recoil spectra vs. A,energy,Z,theta,x,y,z
c
c  RECLVEL:
c      stores recoil velocity distributions vs. emission no.,vs. x,y,z,v-resultant
c
c  RECOIL:
c      does the kinematics to boost emitted particles into lab frame,and to
c      adjust recoil momenta each time a nucleon is emitted(or alpha)
c      It is assumed that the recoil initially has the momentum of the
c      incident nucleon
c
c   RECOILWR:
c       writes the various recoil arrays stored on energy,mass,atomic no.
c       and lab angle
c
c   ANGWRT:
c        writes the cm and lab frame double differential cross sections
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c****************************************************************************
c
c*****************************************************************************
c
c       some discussion of code physics may be found in phys. rev.28,
c       1475(1983),and in llnl report no. ucid 19614(1982),and in
c       references therein.
c       A manual is also available valid up to 2008
c
c **********************************************************************
c
c
c
      program alice14
      common/kkg2/fractot(99),cs90
      common/par3/eq,sigml(999),acrs(999)
      common/gamcl/igamcl
      common/atsv/atsave
      common/fissngl/dissa(200),dissu(200),dissz(200)
       common/rcsrec/rcsorig(100)
      doubleprecision dissa,dissu,dissz
      common/libst/a53t(260)
      common/fms/sdissa,sdissz
      doubleprecision sdissa,sdissz
      common/emin/eminf(106,160)
      doubleprecision dafr1,dafr2,dzfr1,dzfr2,dufr1,dufr2
      common/azfrag/dafr1,dafr2,dzfr1,dzfr2,dufr1,dufr2
      common/krld/powkr(9999),gamkr(9999),powgkr(9999)
      REAL*8 powkr,gamkr,powgkr
      common/pmasb/pmass(100,150,30)
      common/fizzc/popfiz(13,20,150),pofizbig(13,20,150,30)
      common/gams/gamgam,csgam,csgamm(100)
      REAL*8 gamgam,csgamm,csgam
      common/outlim/limout,lim2
       common/fizz/fizbar(106,260),fizcs(100)
      common/fermst/aferm
      common/ismrdat1/englev(1000,3),spingnd(106,260)
      common/ismrdat2/ididl(106,260),nolevs(106,260)
      common/ismrdat/isoid(106,260),levnm(1000),bevj(1000,3)
      common/rotor/irotor,iidummy,arote(15,25,181)
       common/scrn/att(10),cldd(10),bexp(24,15),reacss(100),edelt
       common/scrn1/ingo,neng,nengg,ireact
        common/mccom/ncount,ipot,rum(15,24),ppmc(15,24,999),cf
      common/giani/pre, sigt(1100),sigam(999),sigpre(999),gsp(999)
     2,eqgam(999)
      common/rcc/rcsi(10,100)
      common/qq/gow(1100)
      common/memo/scrs(999,8)
      REAL*8 gow
      REAL*8 q,gamft
      common/m/ q(999),sp(999),sifis(999),gamft(999),t(999),pairx(8)
     1,scale(24,15),rd(4)
      common/temscale/a2h,a3h,a3he,a4he,a4hep
      common/temscalf/b2h,b3h,b3he,b4he,b4hep
c     common/exclus/ gib(2,999),pr(3),ge(2100),zz(3),
c    1delta(3),ts(3),pairx(8),asum(10),rate(2,999)
c     common/m/ q(999),sp(999),sifis(999),gamft(999),t(999)
c    1,scale(24,15),rd(4)
      common/cs/crsum
      common/lab10/pow(4,9999),gam(9999)
      common/lab12/pof(4,9999)
      REAL*8 pof
        REAL*8 pow,gam
      common/lab11/powaz(15,24,2200)
      common/scat/jz,ja, probxx(2,2),q1(8,1100),fcs,eloss,fcspec
      common/isomer1/xgs(690),xiso(690),x2iso(690),isonum,isoint,iso2num
      common/isomer2/zs(690),as(690),eiso(690),eiso2(690)
      common/isomer3/iziso(690),iaiso(690),csgs(690),csis(690)
     1,csgse(690,100),csise(690,100),csis2e(690,100),rm1(690,100),
     2rm2(690,100)
      common/error/rerror(15,24)
      common/evapmc/jcount,ipremc,eminn(15,24)
      common/eder/ed2,v,vhalf,vsq
      common/incr/ed,tened,efinal
      common/sf6/aoz,zoz,en(8,999)
      common/shft/mc,mp,inver,ike,ipch,kq,qval,ap,at,zp,zt,cld,barfac
      common/ug/iz ,ia
      common/nhy2/gdo,bisp
      common/pl2/sux(24,13,100)
      common/s1/jfrac,jupper
      common/sf/m3,kplt
      common/sft9/k5,jmax,pld,c(8),delt,alt
      common/par2/cncss
      common/pl8/jamma(15,24),nulim(15,24)
      common/pl1/ecm(100)
      common/nam/k2
      common/fis/ifis
      common/levop/ldopt,igam
      common/hjk/jang,iq,rcss
      common/ss/sor,rr
      common/nr34/nr3,nr4,ke5,i3d,irtst,i3t2,ijkl,lq,tem(36)
      common/parfs/mz,k6,delrr(999)
      common/shft4/k3,jcal
      common/shft5/fiss
      common/shft2/fs(25),dsp(25),brr(25),der(25),er(25)
     1,cx(25),bilsum,xmiss,sumiz
      REAL*8 fiss
      common/pl4/crs(24,13)
      common/nhy/ij,jl,ji,jq,td,ex1,ex2,tmx,av,gav,cost,b(3)
      common/paro/pq,cross
      common/hyb2/pp(15,24,999)
      common/sft5/exc(18,27),xmax
      common/pl3/jx,nepr,eb(100),rcsp(100),zee,amass,plex
      common/lym1/be(18,27,8),pair(18,27),xmas(18,27),delshl(18,27),amas
     1(18,27)
c     common/bym1/symb(18,27),symbp(18,27),symbs(18,27)
      common/pl7/mas
      common/lab3/sig(8,999),slav(8,999)
      common/pl5/na,nz,title(20)
      common/scr/k9
      common/iso/qpn(3),qpnc
      common/tst/test
      common/send/irfr,iadst,dlt
      common/gam/io,la,ppgam(999),gmspec(999),beta,ams
      common/eqin/eind(100)
      common/ist/isot,iikl,knz,idelm,fracts,fract,ccrs(32,13,100)
     1,engbuf(100),exd,een(8,999,100)
      common/cshel/shel(18,27,2)
      common/gr/rint(1100)
      common/deform/beta2(106,260)
      REAL*8 rint
      REAL*8 rtem
      common/mbc/sigg(21,4,70)
      common/rc/rzero
      common/equ/equiv
      common/nuopt/kinem,iem,kb
      common/bound/neut2,iprot,idnuc(100),idm(100),einc(100),
     1thett(100),phit(100),ppp
      common/spins/rote(15,24,181),popej(15,24,180,100),fcx(180),ai,ecmp
      common/lscpl/bx(180),sppj
      common/sping/pfpej(600,180,100),pfpj(600,180)
      common/gate1/spind(5,180),popj(15,24,180),egatl(5),egatu(5),npart(
     15),igate
       common/spins3/ajmax
       common/spins4/eprot,ealpha,eoff(8)
      common/tagspin/trgspin
      common/tgspin/levno
      common/estar/estarcnt,e6,e7,e8
      common/sigsav/sigdr(100),sigqd(100)
      common/frags/xmasf(160,106),bq(106,160,4),pf(106,160)
     1,ppf(106,160,100),ppfi(106,160),ppfz(106),ppfa(260)
     2,ppfaz(106,260)
      common/fissiso/csgsef(690,100),csisef(690,100),csis2ef(690,100)
      common/frags2/rotef(260,181),arotef(260,181)
      common/new/fpyismr(3,100)
      common/update/month,iday,iyear,nen
      common/clustr/m3set
      common/endvparm/iend,lll(6,6,6),ll(6,6)
      common/astp/tim(1)
      common/astp3/niso,no,iadlt,isotag,massno(12),mass1(12),ncsave
      common/astp2/fabund(10),frstabun(10)
      common/isotg/iztag,imtag,tix,abundc,unct
       REAL*8 ucla,pt,gam11
      REAL*8 powaz
c
      ingo=0
      month=02
      iday=01
      iyear=2016
c enter month day year for id purposes of code version used
c xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
cc
      open (7,file='plot',status='unknown')
      open (8,file='engang',status='unknown')
      open (9,file='sigs',status='unknown')
      open (10,file='labout',status='unknown')
crc
      open (11,file='labangd',status='unknown')
      open (12,file='labengd',status='unknown')
      open (13,file='recvel',status='unknown')
      open (14,file='recoils',status='unknown')
      open (15,file='insave',status='unknown')
c     open (16,file='fout',status='unknown')
      open (17,file='outpt',status='unknown')
      open (18,file='angeng',status='unknown')
      open (19,file='mcspec',status='unknown')
      open (20,file='labspec',status='unknown')
      open (21,file='fizzers',status='unknown')
      open (22,file='fizzersj',status='unknown')
      open (23,file='fizisomr',status='unknown')
      open (24,file='fizisrat',status='unknown')
      open (25,file='fizylds',status='unknown')
      open (26,file='fizybya',status='unknown')
      open (27,file='ejplanes',status='unknown')
c     open (28,file='billiso',status='old')
      open (28,file='phtlib',status='old')
      open (29,file='fragss',status='unknown')
      open (30,file='tmadland',status='old')
      open (31, file ='MASS.TBL', status='old')
      open (32, file ='SHELL.TBL', status='old')
      open (33,file='finp',status='unknown')
      open (34,file='fout',status='unknown')
      open (35,file='pmtable',status='old')
      open (37,file='fizej',status='unknown')
      open (38,file='gdrparms',status='old')
      open (39,file='sdcsfcm',status='unknown')
      open (40,file='ddcsfcm',status='unknown')
      open (41,file='sdcsflb',status='unknown')
      open (42,file='ddcsflb',status='unknown')
      open (43,file='ddcstotlb',status='unknown')
      open (44,file='sdcstotlb',status='unknown')
      open (45,file='sdcscm',status='unknown')
      open (46,file='abund',status='old')
      open (50,file='excdst',status='unknown')
      open (51,file='excdst1',status='unknown')
      open (52,file='exciton',status='unknown')
      open (53,file='ratiox',status='unknown')
      open (54,file='gammaspec',status='unknown')
      open (55,file='sdcscl',status='unknown')
      open (56,file='sdlbcl',status='unknown')
      open (57,file='temp',status='unknown')
      open (58,file='endf',status='unknown')
c
c-----------------------------------------------------------------------
      
       write(*,*)
       write(*,*)'          Welcome to Alice"s Wonderland  2014'
       write(*,*)
       write(*,*)'               version Sept.2014'
       write(*,*)
       write(*,*)' J.Bisplinghoff, W.Scobel,W.B.Wilson,S.G.Mashnik'
       write(*,*)' Yu. Konobeev,A.V.Ignatyuk,Y.Shubin,V.P.Luneev,'
       write(*,*)'             K.K.Gudima,M.Mebel,M.Blann'
       write(*,*)
       write(*,*)' Follow the screen directions to avoid a bad hare day'
       write(*,*)
       write(*,*)' Unless otherwise stated,all formats are free format'
       write(*,*)
       write(*,*)'  To read a standard file "INSAVE" enter "0" now,'
       write(*,*)
       write(*,*)'  or for input via SCREEN PROMPT,enter "1" now'
       read(*,*)inopt
      if(inopt.eq.0)write(*,*)' program will read input from file "INSAV
     1E"'
      write(*,*)' All spectra in c.m. system are given with respect to'
      write(*,*)' channel angle for each emission,not referred back to' 
      write(*,*)'                 beam direction'
       write(*,*)
      write(*,*)'   Questions or problems? email  hmblann@gmail.com '
       write(*,*)
c
      call headers(1)
      call zero(1)
      call zero(6)
c xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      call findismr
c xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
c below,put experimental fission barriers in for some heavy elements
c at present there are 66 entries in file 30,tmadland
c     do 12010 i=1,66
c     read(30,12011)izia,bfsr
c     izz=float(izia)/1000.
c     iaa=izia-1000*izz
c      fizbar(izz,iaa)=bfsr
c12010 continue
c      rewind(28)
c store table of a**5/3
      do i=1,260
      a=float(i)
      a53t(i)=a**(5./3.)
      enddo
c
c following 11 parameters fixed here rather than read from input file Oct 31-01
      aferm=12.
      nepr=1
c      isoint=1
       iparm=0
       jcal=1
       pld=9.
       m3=7
       igate=0
       td=1.
       iend=0
       
c---------------------call screenin input interrogation(screenin) or 
c---------------------existing input file 'insave' via 'stdinpt'---
         if(inopt.eq.0)call stdinpt
       if(inopt.eq.1)call screenin
c-------------------------------------------------------------------
       if(kinem.ne.0)kinem=1
c-------------------------------------
c -------------------------------
       call evapconst
c-------------------------------------------------
c identify vectors for ENDF output if selected---
        if(iend.gt.0)call endvsel
c-------------------------------------------------
       atsave=at
       noiso=0
       niso=1
       izt=zt
c---------------------------------s1--------------------------------
                    if(at.eq.1000.)then
         call isotagr(1)
      write(*,*)
      write(*,*)' An isotopic target is assumed:'
      write(*,*)' At# Mass# % abundance % uncertainty'
c below we find the stable isotopes of the target and % abundances for nat. tag
       do ii=1,287
          read(46,3333)iztag,imtag,tix,abundc,unct
                              if(iztag.eq.izt)then
           write(*,3333)iztag,imtag,tix,abundc,unct
            call isotagr(2)
            noiso=noiso+1
            frstabun(noiso)=abundc/100.
             mass1(noiso)=imtag
             niso=noiso
                                   endif
        enddo
                                   endif
              write(*,*)
c---------------------------------s1-------------------------------
c re-order isotopes of target heaviest to lightest

                     if(niso.gt.1)then
        do no=1,niso
            massno(no)=mass1(niso+1-no)
            fabund(no)=frstabun(niso+1-no)
        enddo  
         at=massno(1)
          atsave=massno(1)
                                  endif
               if(niso.lt.1)niso=1
               ksot=niso
c-------------------------------s1----------------------------------
c this ends preliminary housekeeping for natural isotopic targets
         att(1)=at
        lim2=0
        if(limout.eq.2.or.limout.eq.0)lim2=1
         if(limout.eq.2)limout=1
       if(ldopt.eq.-1)then
       call ldmagic
       endif
       
c aug 06 add zero to sig arrays
        do ii=1,8
         do kk=1,999
          sig(ii,kk)=0.
           enddo
            enddo
       continue
c ------------------------------------------
        ncsave=ncount
       isotag=0
       iadlt=0
       if(niso.gt.1)isotag=1
c++++++++++++++++++++++++s2+++++++++++++++
c    7/7/14
         call prelim
c+++++++++++++++++++++++++++++++++++++++s2+++++++++++
c now translate input switches
       if(neng.gt.0)then
       do 951 ien=2,neng
  951  eind(ien)=eind(ien-1)+edelt
       endif
       nengg=abs(neng)
          aum=eind(nengg)/ed
        if(aum.gt.990.)then
       nr=aum/990.
       ed=ed*(1.+float(nr))
      write(*,*)'maximum energy/bin width problem;reset width to ',ed,
     1'MeV'
      write(7,*)'maximum energy/bin width problem;reset width to ',ed,
     1'MeV'
       endif
      iem=40./ed
c
      call lymasf
      call gem
      call sigf
c
c  0000000000000000000000xxxxxxxxxxxxxxxxxxxxxxxxxxxxx000000000000000000
        iso=1
       if(jcal.eq.0)ifis=0
       if(jcal.eq.10)ifis=0
       mc=10
       mp=1
       if(cld.eq.0.)cld=1.02
c nov 7,2002 made cld default 1.02 rather than 1.00
       if(barfac.eq.0.)barfac=1.
       if(ed.eq.0.)ed=0.5
       ike=4
       if(pld.eq.0.)pld=9.
       if(m3.eq.0)m3=3
       gdo=0.
       tmx=0.
       if(td.eq.0.)td=1.
cbound
       if(na.eq.0)na=22
       if(nz.eq.0)nz=13
       neut2=at+ap-zp-zt-1
       if(na.ge.neut2)na=neut2-1
        iprot=zp+zt-1
       if(nz.ge.iprot)nz=iprot-1
      sor=0.7*sqrt(100.*(ap+at))
      kplt=1
      igam=1
      if(ike.gt.4)ike=ike-5
      pre=0.
      if(at.ne.0.)iecoun=0
      if(ed.eq.0.)ed=0.5
      ams=ap+at
       iparm=0
      iparm=mc+mp+na+nz
      if(iparm.gt.0)go to 80001
      mc=10
      mp=3
      inver=2
      ike=4
       if(limout.eq.0)then
      if(at.ne.0.)write(17,80002)
      if(ipch.eq.2)write(17,9877)
      if(ipch.eq.1)write(17,9878)
       endif
80001 continue
      iikl=0
      ams=ap+at
      pq=0.
      k5=0
      k9=0
      if(inver.eq.2)k9=1
      plex=nz
      k6=0
      gim=0.
      k2=1
      call plt
      nepr=0
      k3=0
c     if(ipch.eq.0)go to 18
c     do 16 iz=1,nz
c     read(16,9876)(bexp(ia,iz),ia=1,na)
c  16 write(17,9876)(bexp(ia,iz),ia=1,na)
c
c  18 read(16,7005) title
   18 continue
      izp=zp
      iza=ap
c---set proj spin 'sppj'-------------------
      if(ap.eq.0.)sppj=1.0
      if(ap.eq.1.)sppj=0.5
      if(ap.eq.2.)sppj=1.0
      if(ap.eq.3.)sppj=0.5
      if(ap.eq.4.)sppj=0.
      if(ap.gt.4.)sppj=spingnd(izp,iza)
c--------------------------------------
c
       call spincpl
       if(limout.eq.0)call tagout
       write(17,*)zp,ap,zt,at,isot,levno,' ZP,AP,ZT,AT,ISOT,LEVNO'
crc
c
      if(barfac.eq.0.)barfac=1.
c------------assure a,z indices within array limits-----------
      if(nz.gt.13)write(17,7010)
      if(na.gt.22)write(17,7015)
      if(nz.gt.13)nz=13
      if(na.gt.22)na=22
c-------------------------------------------------------------
      if(m3.gt.4)m3=7
      if(m3.eq.0)m3=3
      if(pld.eq.0.)pld=9.
      amass=ap+atsave
      zee=zp+zt
      aneut=amass-zee
      rr=amass/pld
      tened=10.*ed
      if(aneut.lt.24)na=aneut-2
      if(zee.lt.15)nz=zee-2
c --set jmax, number of inverse x-sections to calc each particle--
      jmax=float(90)/ed
      if(jmax.gt.999)jmax=999
c----------------------------------------------------------------
      write(17,7466)ldopt
      write(*,*)ldopt
c
      if(na.gt.0.and.nz.gt.0)go to 7025
      nz=5
      plex=nz
      na=10
 7025 if(cld.le.0.)cld=1.02
c 7 nov 2002 made cld default 1.02
      continue
c
      if(ldopt.eq.1)mp=1
c
c     print out input options from sr optout,then get liquid drop mass values from lymass
      call optout
      write(*,99)sppj
 99   format(' Projectile spin=',1f5.2)
      amass=atsave+ap
      call lymass(zee,amass,nz,na,mc,mp,ap,at,zp,zt,qval,ldopt)
c     get separation energies and level density ground state shifts
      call zero(4)
      ia=1
      iz=1
c     load level densities as per option selected or default Fermi gas
c     put fermi gas level densities in,to be overwritten if other options selected
      if(ldopt.eq.0)call fermld
c temp to call ldcalc
      if(ldopt.eq.1)call ldcalc
c 1= Kataria- Ramamurthy
      if(ldopt.eq.2)call obninsk
c 2= Obninsk/ Ignatiuk et al
      if(ldopt.eq.3)call ldload(ldopt)
c 3= MB Chadwick for elements Z<10,N<10
      if(ldopt.eq.4)call ldload(ldopt)
c 4= Gilbert - Cameron
cc    initialize buffers relevant to gamma emission,additional gamma coding -------------
      pre=1.
      call gama
c-------------sub4--------------------------------
      pre=0.
      kx=150./ed
      if(kx.gt.999)kx=999
      rint(1)=0.
      do 6144 it=2,kx
      do 6143 if=2,it
      jf=if-1
      eg=float(jf)*ed
      ib=(float(it)*ed-eg)*10.
      if(ib.lt.1)ib=1
      if(ib.gt.999)go to 6143
      rtem=sigt(jf)*gow(ib)*0.5
      rint(it)=rint(it)+rtem
 
 6143 continue
c   note that sigt(i) is e*e*sig(e) from subroutine gamma
 6144 continue
      do 6145 it=2,kx
      rint(it)=rint(it)*ed
 6145 continue
c-------------sub4---------------------------------------
cc    this is end of gamma emission related coding-for now  ----------------------
c     integer variables for compound nucleus mass and atomic number
      nia=amass
      iip=zee
c
c     get inverse cross sections, see call sigi also
c--------------------------s6-------------------------------------
      if(inver.eq.0)go to 35
      if(inver.eq.1)call rdinv(m3,jmax)
      if(inver.eq.2)call sigi
      if(inver.ge.2)go to 44
      continue
c
      go to 40
c
   35 call over
c
   40  k3=3
c
      call shaft
   44 continue
      ldex=-1
c MAR 970414
**
** add code to determine the minimum energy and paricle type (n, p, or alpha)
** necessary for particle emission for all nuclei.
      ix = 40./ed + 1.
      iproton = 0
      ideut=0
      itrit=0
      ihe3=0
      i7be=0
      do 41 i = 1,ix
         if(sig(2,i) .gt. 0.)go to 41
         iproton = i
 41   continue
      ialpha = 0
      do 42 i = 1,ix
         if(sig(3,i) .gt. 0.)go to 42
         ialpha = i
 42   continue
      ideut = 0
      do 842 i = 1,ix
         if(sig(4,i) .gt. 0.)go to 842
         ideut = i
  842   continue
      itrit = 0
      do 843 i = 1,ix
         if(sig(5,i) .gt. 0.)go to 843
         itrit = i
 843  continue
      i3he = 0
      do 844 i = 1,ix
         if(sig(6,i) .gt. 0.)go to 844
         i3he = i
 844   continue
      i7be = 0
      do 845 i = 1,ix
         if(sig(7,i) .gt. 0.)go to 845
         i7be = i
 845   continue
 
      eprot  =float (iproton) *ed
      ealpha = float(ialpha)*ed
c  1/23 try mod eprot,ealpha
c      eprot=eprot+ed
c      ealpha=ealpha+ed
c--1/23
      edeut=float(ideut)*ed+ed
      etrit=float(itrit)*ed+ed
      e3he=float(i3he)*ed+ed
      e7be=float(i7be)*ed+ed
      eoff(1)=0.
      eoff(2)=eprot
      eoff(3)=ealpha
      eoff(4)=edeut
      eoff(5)=etrit
      eoff(6)=e3he
      eoff(7)=e7be
c eoff array to adjust emission momenta to above barrier values in excl1,evap
      do 43 iz = 1,nz+2
         do 43 ia = 1,na
            beffp = be(iz,ia,2) + eprot+ed/2.
            beffa = be(iz,ia,3) + ealpha+ed/2.
           
            eminn(iz,ia) = min(be(iz,ia,1),beffp,beffa)
            eminn(iz,ia)=eminn(iz,ia)+ed+ed/2.
c           eminn(iz,ia)=eminn(iz,ia)+ed
 43   continue
      ratem=amass**0.3333
      do 8215 iz=10,106
      if(iz.le.30)then
      nlo=iz/2
      nhi=1.5*float(iz)
      endif
      if(iz.le.30)go to 25
      if(iz.le.50)then
      nlo=iz
      nhi=2*iz
      endif
      if(iz.le.50)go to 25
      if(iz.le.80)then
      nlo=1.2*float(iz)
      nhi=1.7*float(iz)
      endif
      if(iz.le.80)go to 25
       if(iz.le.106)then
       nlo=1.37*float(iz)
       nhi=1.83*float(iz)
       endif
   25 continue
       if(nlo.lt.10)nlo=10
       if(nhi.gt.160)nhi=160
      do 8215 n=nlo,nhi
      ia=n+iz
      z=iz
      un=n
      a=ia
      afac=ratem/a**0.333
            beffp = bq(iz,n,2) + afac*eprot*float(iz)/zee
            beffa = bq(iz,n,3) + afac*ealpha*float(iz)/zee
            eminf(iz,n) = min(bq(iz,n,1),beffp,beffa)
c           eminf(iz,n)=eminf(iz,n)+ed/2.
ctemp put ed in place ed/2
            eminf(iz,n)=eminf(iz,n)+ed
 8215 continue
c------------------------------s6------------------------------
c july 2001 add loop to store rotational energies from sierk routine
c
      if(irotor.eq.0)go to 10041
      do 10040 jz=1,13
        iiz=zee+1-jz
        if(iiz.lt.20)go to 10040
      do 10039 ja=1,22
        iai=amass-ja-jz+2
       allim=0.
      do 10038 il=1,90
       jallim=allim
       if(il.gt.5.and.il.ge.jallim)go to 10038
        jl=il-1
       call barfit1(iiz,iai,il,ibfis,irot,ilmax)
c      call barfit(iiz,iai,il,ibfis,irot,ilmax)
       allim=float(ilmax)/10000.
       j2p1=2*jl+1
       j2p2=j2p1+1
        if(j2p1.lt.1)j2p1=1
        if(j2p2.lt.1)j2p2=1
c above changes rotational energy arrays to 2J+1 index
        rote(jz,ja,j2p1)=float(irot)/10000.
        rote(jz,ja,j2p2)=float(irot)/10000.
       if(il.eq.1)rote(jz,ja,j2p1)=0.
       if(il.eq.1)rote(jz,ja,j2p2)=0.
       if(irotor.gt.0)arote(jz,ja,j2p1)=rote(jz,ja,j2p1)
       if(irotor.gt.0)arote(jz,ja,j2p2)=rote(jz,ja,j2p2)
       if(irotor.eq.0)arote(jz,ja,j2p1)=0.
       if(irotor.eq.0)arote(jz,ja,j2p2)=0.
10038  continue
10039  continue
10040  continue
c make rotational energy arrays for fission frags dec'03
       do iz=20,90
       z=float(iz)
       a=2.4*z
       iiz=z
       iai=a
       arotef(iai,1)=0.
       arotef(iai,2)=0.
       allim=0.
       do il=2,88
       jl=il-1
       call barfit1(iiz,iai,il,ibfis,irot,ilmax)
c      call barfit(iiz,iai,il,ibfis,irot,ilmax)
       allim=float(ilmax)/10000.
       j2p1=2*jl+1
       j2p2=j2p1+1
c above changes rotational energy arrays to 2J+1 index
        rotef(iai,j2p1)=float(irot)/10000.
        rotef(iai,j2p2)=float(irot)/10000.
       if(irotor.gt.0)arotef(iai,j2p1)=rotef(iai,j2p1)
       if(irotor.gt.0)arotef(iai,j2p2)=rotef(iai,j2p2)
       if(irotor.eq.0)arotef(iai,j2p1)=0.
       if(irotor.eq.0)arotef(iai,j2p2)=0.
       enddo
       enddo
       
10041  continue
       continue
       continue
c MAR 970414
 
c       xxxxxxxxxxxxxxxxx00000000000000000000000000000xxxxxxxxxxxxxxxxxx
      do 12010 i=1,66
      read(30,12011)izia,bfsr
      izz=float(izia)/1000.
      iaa=izia-1000*izz
       fizbar(izz,iaa)=bfsr
12010 continue
       rewind(28)
      call zero(8)
c ))))))))))))))))))))))))))) energy loop )))))))))))))))))))))))))
17453 dlt=ed
      do 1051 nen=1,nengg
      nenged=(eind(nen)+ed)/ed
      if(nenged.gt.990)then
      write(*,*)' excitation energy marginal: re-run with larger energy',
     1,'mesh size'
c23456789012345678901234567890123456789012345678901234567890123456789012
      if(eind(nen).gt.250.)write(*,*)' Excitation exceeds intended',
     1'upper range; use results with care!'
      go to 1051
      endif
c      call zero(5)
       call zero(7)
c
       gamgam=0.
       csgam=0.
c
c
      eq=eind(nen)
      estarcnt=0.
      e6=0.
      e7=0.
      e8=0.
      write(*,98)eq
  98  format(' now doing energy =',1f6.2)
      rcss=0.
      if(ireact.ge.0)rcss=reacss(nen)
      if(ireact.lt.0)rcss=-1.
c     av=0.
      if(ncount.gt.0)write(17,7789)ncount
      if(ncount.gt.0)write(*,7789)ncount
      if(ipot.gt.0)write(17,7780)
      if(rzero.eq.0.)rzero=1.5
      if(ncount.eq.0)go to 17453
      do 17152 ijk=1,15
      do 17252 jjk=1,24
      rum(ijk,jjk)=0.
      do 17353 kjk=1,999
      ppmc(ijk,jjk,kjk)=0.
      pp(ijk,jjk,kjk)=0.
17353 continue
17252 continue
17152 continue
      ldex=ldex+1
c
c  set more default parameters for case that only a,z of target/projectile and e incident input
      if(iparm.eq.0)jcal=1
      if(iparm.eq.0)td=1.
      if(iparm.eq.0.and.ap.eq.1.)tmx=1.
c
c for case of isotopic weighting,check to be sure that the incident energies agree each isotope
c     dif=0.
      iikl=iikl+1
      engbuf(iikl)=eq
c
      if(jcal.eq.0.and.ifis.eq.0)write(17,6781)
      if(jcal.eq.0.and.ifis.gt.0.and.ifis.ne.5)write(17,6782)
      if(jcal.eq.0.and.ifis.eq.5)write(17,6783)
c
c
c if incident energy negative,branch to graphics,zero read next problem line one,else continue
c
c     calc max. excitation each nuclide, evaluate and print options
      jdelt=1
      iecoun=iecoun+1
c     ecoun(iecoun)=eq
      if (jang.gt.200) jdelt=jang-200
      if (jang.gt.200) jang=1
      if (jang.gt.100) jdelt=jang-100
      if (jang.gt.100) jang=0
      jsw=0
      if (jcal.ge.9) jsw=1
      if (jcal.ge.9) jcal=jcal-10
c
      if(td.gt.0..and.jcal.gt.1)write(17,82)
      if(td.gt.0..and.jcal.gt.1)jcal=0
c
      if (td.eq.0..or.jcal.ne.-1) go to 81
      write (17,82)
      jcal=0
   81 do 55 m=1,999
      delrr(m)=0.
      t(m)=0.
   55 sigml(m)=0.
      do 60 n=1,999
      do 60 m=1,8
   60 en(m,n)=0.
      if(jcal.ne.1)k6=1
c  62 continue
       isotag=0
       iadlt=0
       no=1
       if(niso.gt.1)isotag=1
        ncount=ncsave
       call zero(5)
c----xxxxxxxxxxxxxx--------------xxxxxxxxxxx isotope =loop xxxxxxxxxxxxx
       do 10004 no=1,niso
        if(isotag.eq.0)then
          iadlt=0
           go to 17005
                       else
       iadlt=massno(1)-massno(no)
       ncount=float(ncsave)*fabund(no)
                       endif
17005       continue
       jz=1
       ja=1+iadlt
c here iadlt id's isotopic tag index
       if(niso.gt.1)at=massno(no)
       amass=atsave+ap-float(iadlt)
       zee=zp+zt
       na=amass-zee-2
       if(na.gt.22)na=22-iadlt
       nz=zee-2
       if(nz.gt.13)nz=13
       lz=zee
       la=amass
c jan 02 correct to target mass/z from compound
       kz=zt
       ka=at
c correct tag mass for this isotope
       trgspin=spingnd(kz,ka)
c following is to allow isomer targets
       qdelt=0.
       if(levno.gt.1)then
       n=isoid(kz,ka)
       trgspin=bevj(n,levno)
       qdelt=englev(n,levno)
       write(*,*)
       write(*,*)'levno=',levno,' target spin reset for isomer=',trgspin
       write(*,*)
       write(*,*)' q value incremented by',qdelt,' MeV'
       write(7,*)'levnno=',i2,'  target spin reset for isomer=',trgspin
       write(7,*)
       write(7,*)' q value incremented by',qdelt,' MeV'
       write(7,*)
                    endif
       if(limout.eq.0)then
       write(27,*)trgspin,'   TARGET SPIN '
                      endif
       write(*,*)
c
c here be sure that a,z are properly defined
      if(jcal.ne.1)call fisrot(a,z,an,al,delr,delsp,ero,barfac)
      if(jcal.eq.0.and.ifis.eq.5)go to 7722
      go to 7822
7722  continue
c
c calculate center of mass incident energy
 7822     pq=eq*at/amass
      if(no.gt.1)pq=eq*float(massno(no))/(massno(no)+ap)
      bar=300.
      afis=0.
      if(no.eq.1)nepr=nepr+1
      qval=0.
      if(noiso.gt.1)then
       ax=massno(no)
       amasx=massno(no)+ap
        else
         ax=at
         amasx=at+ap
                   endif
      call bind(zee,amasx,ap,ax,zp,zt,qval)
c               endif
      xmax=pq+qval+qdelt
      if(xmax.gt.300.)write(*,*)' Excitation>300MeV exceeds intended ran
     1ge-proceed skeptically'
      write(17,*)'-----------------------------------------------'
      write(17,*)
      if(nen.eq.1)then
      write(17,97)qval
 97   format(' Qval used from experimental masses (MeV) =',1f6.2)
       endif
      write(17,*)
      write(17,96)pq,xmax
 96   format('Ecm (MeV) =',1f6.2,' Initial excitation (MeV) = ',1f6.2)
      write(17,*)'-----------------------------------------------'
      write(*,*)'-----------------------------------------------'
      write(*,*)
      if(nen.eq.1)then
      write(17,97)qval
       endif
      write(*,*)
      write(17,96)pq,xmax
      write(*,*)'------------------------------------------------'

c correct for isomer excitation
      if(ap.ne.0.)go to 7153
c calculate recoil energy if the projectile is a photon
      erec=0.00053*eq*eq/at
      qval=0.
      eq=eq-erec
      equiv=erec*at
      pq=equiv
      xmax=eq+qval
 
7153  if(xmax.le.0.)go to 10004
      if(xmax.gt.999.)write(17,56)
      if(xmax.gt.300.)write(*,*)' Excitation>300MeV exceeds intended ran
     1ge-proceed with skepticism'
      xtest=xmax/ed
      if(xtest.gt.999.)write(17,5699)
      k3=2
c
      call shaft
c
      if(nepr.lt.101)ecm(nepr) = xmax
      if(nepr.lt.101)eb(nepr)=eind(nepr)
      delr=0.
      delsp=0.
      if(xmax.gt.999.)go to 10004
c
c     entrance channel cross sections to be read in as transmission coefficients
      if(rcss.ge.0..and.no.le.1)go to 75
      jkl=-rcss
      rcss=0.
c
c     read(16,65)(t(l),l=1,jkl)
c
      con=660.*amass/(pq*ap*at)
      do 70 l=1,999
      tl=l-1
      sigml(l)=(tl+tl+1.)*t(l)*con
      rcss=rcss+sigml(l)
   70 continue
c      if(limout.eq.0)then
   75 if(rcss.ne.0..and.limout.eq.0)write(17,305) rcss
c     endif
c calculate photonuclear reaction cross section from photin subroutine
      temp=0.
 
      if(ap.eq.0.)then
      call photin
      write(*,*)'no,at,rcss',no,at,rcss
      cncss=rcss
      rcsp(nepr)=rcss
      coul=0.
      temp=1.
      if(xmax.lt.eminn(1,1+iadlt))csgamm(nepr)=rcss
       idxx=1+iadlt
      rcsp(nepr)=rcss
                endif
      if(temp.eq.1.)go to 95
      if(zp.eq.1..and.ireact.eq.1)go to 90 
      if(zp.eq.0..and.ireact.eq.1)go to 90 
      if(rcss.le.0.)go to 85
      cncss = rcss
      if(ireact.eq.1)cncss=reacss(nen)
      go to 101
   85 if(zp.lt.2.)go to 90
c
c   calculate reaction cross section with parabolic model for z g.e.2
      call parap(mx,jcal)
c
      go to 95
c call optical model for n,p,or d reaction cross sections-
   90 i3=m3
      ax=atsave
      if(no.gt.1)ax=massno(no)
      if(zp.eq.0.)m3=1
      if(zp.eq.1..and.ap.eq.1.)m3=2
      if(zp.eq.1..and.ap.eq.2.)m3=4
      if(zp.eq.1..and.ap.eq.3.)m3=5
      coul=0.45*zt/(ax**0.3333)
       if(ireact.eq.1)go to 9101
      if(zp.eq.0.)coul=0.
      coultst=0.6*coul
      if(pq.le.coultst.and.zp.ge.1)then
       cross=0.
       go to 9101
                                endif
c
      if(ireact.eq.0) call over
c
 9101     m3=i3
      cncss=cross
      if(ireact.eq.1)cncss=reacss(nen)
       reacss(nen)=cncss
       rcss=cncss
      if(cncss.eq.0.)go to 1051
   95 continue
       if(zp.ge.2.)reacss(nen)=cncss
       write(*,94)nen,reacss(nen)
   94 format(' Energy index', i3,'  Reaction cross section=',1f9.4)      
c put in small rcss for proton below cb
       if(cross.eq.0.)cross=0.000000001
c here call routine to get partial cross sections for userprovided rx cs
 101  continue
c oct 01 try to get l values from user input reaction cs
      jlim=1
      jj=1
      icall=0
      if(ireact.eq.1)call crslav(nen,icall,jlim)
      continue
      if(inopt.eq.0)ratio=1.
      do 105 i=1,100
      jj=i
      if(sigml(i).lt..001)go to 110
  105 continue
  110 continue
      jlim=jj
      jj=jlim
      continue
      write(17,*)jj,jlim,'jj,jlim'
      if(isot.gt.1)rcsi(ksot,iikl)=cncss
      if(jj.le.0)jj=1
      rcss=cncss
      if (tmx.eq.1.) tmx=jj
      jjj=jj
      if(jfrac.gt.0)go to 125
      ji=1
      go to 130
  125 ji=jfrac + 1
  130 if(jupper.ne.0) jj = jupper + 1
      if (jcal.gt.0) go to 131
      if(jcal.eq.0.and.ifis.eq.5)go to 131
      if (jfrac.le.0.and.jupper.le.0) go to 131
      buff=0.
      do 132 lin=1,999
      if (lin.lt.ji) sigml(lin)=0.
      if (lin.gt.jj) sigml(lin)=0.
  132 buff=buff+sigml(lin)
      rcss=buff
      cncss=buff
  131 if(jcal.le.0.and.ifis.ne.5)go to 140
      ji=1
      jj=1
  140 if(nepr.lt.101)rcsp(nepr) = rcss
c
       ai=jlim
c
       if(ap.eq.0.)cncss=rcss
      do 1140 ii=1,100
       fcx(ii)=sigml(ii)/cncss
      if(fcx(ii).le.0.001)go to 1140
      jlim=ii
1140  continue
      ecmp=eq*at/(at+ap)
c define center of mass energy for spin transfer calcs
 
      crsum=0.
      bisp=0.
      totfis=0.
      fiss=0.
      sumiz =0.
c     begin calc on partial cross section on j
c     j is the angular momentum index;jj is the maximum and ji the minimum
c     if jcal is 1,a single calc is done for total reaction cross section,
c     but if jcal is 0,it is done partial wave by partial wave,summing
c     over partial reaction cross sections on index jl
c     111111111111111111111111111111111111111111111111111111111111
      if(jj.gt.100)write(17,9610)
      if(jj.gt.100)jj=100
      do 610 jl=ji,jj,jdelt
      crsum=crsum+sigml(jl)
      ajl=jl
      test=5*(jl/5)-ajl
      if(jl.le.2) test=0.
      if(no.le.1)then
      do 145 mm=1,999
      ppgam(mm)=0.
      do 145 mk=1,24
      do 145 ml=1,15
  145 pp(ml,mk,mm)=0.
                endif
      lmx=1
      lmn=1
      al=jl-1
      jll=al
      if(jcal.lt.0.and.td.gt.0.)go to 165
      if(jcal.ge.0)go to 165
      if(acrs(jl).le.0.)go to 610
      pp(1,1,1)=acrs(jl)
      go to 170
  165 pp(1,1,1)=sigml(jl)
  170 nia=amass
      iip=zee
      if(jcal.ge.1)pp(1,1,1)=rcss
      if(jcal.eq.0.and.ifis.eq.5)pp(1,1,1)=rcss
      if(jcal.eq.0.and.ifis.ne.5.and.sigml(jl).le.0.)go to 610
      test2=pp(1,1,1)/rcsp(nepr)
      if(test2.le.0.001)test=0.
      if(test.ne.0.)go to 171
c
      write(17,640)
      write(17,645)iip,nia
      write (17,685) rcsp(nepr)
      write(17,635)xmax,jll
      write(17,701)pp(1,1,1)
c
  171 iz=1
      if(inver.eq.2)call sigi
      swc=0.
      max=float(20)/ed
c   set up n,p binding energies for precompound
      b(1)=be(1,1+iadlt,1)
      b(2)=be(1,1+iadlt,2)
c
c  read in inverse reaction cross sections for light elements from pre-calcd array
      izea=zp+zt
c MAR 970414
      jcount=1
      ipremc=1
      pp(1,1,1)=0.
c
c  call precompound routine
      call hybrid
      if(xmax.lt.eminn(1,1+iadlt))then
       rcsp(nepr)=rcss
                                  endif
 
  610 continue
      continue
      if(no.lt.niso)go to 10004
      if(ike.eq.0)go to 9602
      k3=4
      call shaft
9602   continue
      do i = 1,22
         cx(i) = 0.
      enddo
      ndmx = 40./ed
      do 9611 iz = 1,nz
         do 9612 ia = 1,na
            sa = 0.
            do 9613 iexc = 1,ndmx
               sa = sa + pp(iz,ia,iexc)
 9613       continue
            cx(ia) = sa
 9612    continue
c                                                    endif
         mas=nia-iz+1
         if(noiso.gt.1)mas=atsave+ap-iz+1
         zoz = zee - iz + 1
         k3 = 1
         call shaft
 9611 continue
c                                                    endif
c
** write out array of particle emmissions
       do 4700 jj=1,15
          do 4690 jja=1,24
             rum(jj,jja)=0.
 4690     continue
 4700  continue
       ixmx=(xmax+ed )/ed
       gsum=0.
       do 7000 jj=1,15
          do 6500 jja=1,24
             do 6400 i=1,ixmx
                rum(jj,jja)=rum(jj,jja)+pp(jj,jja,i)
 6400        continue
             rerror(jj,jja)=sqrt(rum(jj,jja)*fcs)
             gsum=gsum+rum(jj,jja)
 6500     continue
 7000  continue
       pp(1,1,ixmx)=0.
       write(17,*)
       write(17,7008)gsum
       write(17,204)
       write(17,204)
       write(17,7009)
       write(17,204)
 
       do 7101 i =1,15
          do j=24,1,-1
             if(rum(i,j).ne.0.0) then
                k = j
                goto 6999
             end if
          end do
          goto 7101
 
 6999     l=i-1
          write(17,201)l,(rum(i,j),j=1,k)
          write(17,203)(rerror(i,j),j=1,k)
          write(17,204)
 
 7101  continue
 
       write(17,7110)
      if(ike.eq.4)k3=4
      if(ike.eq.4)call shaft
      write(17,670)fiss
c is kplt defined?have moved 152 lines to new sr
       write(23,*)
      write(23,87654) eind(nen),rcsp(nen),fizcs(nen)
       write(23,*)
      call ismrout
      if(kplt.gt.0)k2=2
      if(kplt.gt.0)call plt
      if(lim2.eq.1)call fiznorm
       write(34,87654)eind(nen),rcsp(nen),fizcs(nen)
      write(29,*)' A Distributions for fixed Z'
       do 9620 iz=15,70
       itest=0
       if(ppfz(iz).eq.0.)go to 9620
       do 9621 in=5,120
       if(ppfi(iz,in).eq.0.)go to 9621
       if(itest.gt.0)go to 9621
       itest=1
       z=iz
       a=iz+in
       i=a
       i1=i+1
       i2=i1+1
       i3=i2+1
       i4=i3+1
       i5=i4+1
       i6=i5+1
        i7=i6+1
        i8=i7+1
         i9=i8+1
       write(29,691)iz,i,i1,i2,i3,i4,i5,i6,i7,i8,i9
       write(29,692)(ppfi(iz,n),n=in,in+9)
9621  continue
9620  continue
      write(29,*)' Z Distributions for fixed A'
      do 9630 ia=40,180
      izl=float(ia)*0.4-8
      izh=float(ia)*0.6
      if(izh.le.izl)izh=izl+5
      if(izh.gt.70)izh=70
         itest=0
      do 9631 iz=izl,izh
      in=ia-iz
      if(ppfaz(iz,ia).eq.0.)go to 9631
        if(itest.gt.0)go to 9631
       itest=1
      lz=iz
      iz2=iz+1
       iz3=iz2+1
       iz4=iz3+1
       iz5=iz4+1
       iz6=iz5+1
       iz7=iz6+1
       iz8=iz7+1
      iz9=iz8+1
      iz10=iz9+1
      im=in+1
      write(29,693)ia,iz,iz2,iz3,iz4,iz5,iz6,iz7,iz8,iz9,iz10
      write(29,694)(ppfaz(mz,ia),mz=iz,iz+9)
 9631  continue
 9630 continue
10004 continue
c1051 continue
cxxxxxxxxxxxxxxxxxxxxxx00000000000000000000000000xxxxxxxxxxxxxxxxxxx
c     k2=3
c     call plt
c1151 continue
c0000000000000000000000xxxxxxxxxxxxxxxxxxxxxxxxxx00000000000000000000
c1251 continue
       inopt=0
         atem=ap+atsave
         ztem=zp+zt
        write(34,100)atem,ztem,xmax
          sdissa=0.
          sdissz=0.
c1152 continue
             do  i=1,190
         sdissa=sdissa+dissa(i)
        if(i.le.100) sdissz=sdissz+dissz(i)
           if(dissa(i).gt.0.0)  then
       dissu(i)=dissu(i)/dissa(i)
c	    dissa(i)=dissa(i)*wrun
                               endif
       if(i.le.100.and.(dissa(i).gt.0..or.dissz(i).gt.0.)) then 
           write(34,981) i,dissa(i),dissu(i),i,dissz(i)
           continue
                          elseif(i.gt.100.and.dissa(i).gt.0.) then
           write(34,981) i,dissa(i),dissu(i)
           continue
       endif 
       enddo
c
       write(34,*)'sum A fission=',sdissa,'sum Z fission =',sdissz
 1051    continue
      k2=3
      call plt
      call finish
      write(*,*)
      write(*,*)'  DAS IST ALICE   DAS IST ALICE  DAS IST ALICE  DAS
     1IST ALICE  DAS IST ALICE'
      write(*,*)

  100 format(/' fragments mass,', 
     & ' excitation energy and charge distributions from fission'/
     & '     of nucleus: A0=',F4.0,', Z0=',F4.0,', E*=',F6.1,' MeV'/
     & 3x,'  A ','f(A1)+f(A2)','  <u(A)/A> ',3x,' Z ','f(Z1)+f(Z2)') 
c10009 format('isotopically weighted cross sections and spectra follow')
  693 format(' A = ',1x,i3,1x,' Z = ',10(1x,i3,6x))
  691 format(' Z=',1i3,'A= ',10(7x,1i3))
  692 format(11x,10(1x,1pe9.2))
87654 format('Beam energy (MeV)=',1pe11.4,';  Reaction CX (mb)=',e11.4,
     1 ';  Fission CX (mb)=',e11.4)                                   
  204 format('  ')
  201 format(' Z-',i2,3x,10(1x,1pe9.3))
  203 format('+ - sd =',10(1x,1pe9.3))
 7110 format(' error is +/- one sigma statistical error')
c7788 format(' first emission gamma eq cs =',e10.3,'mb')
 7008 format(' sum over all events = ',f14.5)
 7009 format(' Z indx N-no.= N       N-1     N-2       N-3       N-4
     1    N-5       N-6       N-7       N-8       N-9       N-10    '
     2,'     N-11')
 9610 format(//,20x,'partial wave upper limit reduced to 100 h-bar ')
  305 format (36x,'compound nucleus cross section provided by user ='
     1  ,f8.2,' mb')
c  65 format(10f5.3)
 5699 format(' excitation/ed exceeds 999 dimension limit
     1next incident energy will be attempted')
56    format(' compound nucleus excitation exceeds 999 mev li
     1mit. next bombarding energy will be attempted.')
c7725 format(' mass no.=',f5.0,'at.no.=',f5.0,'fiss barrier=',f8.2)
   82 format (1h ,80('*')/' precompound calculation incompa',
     1 'tible with 1 mev bin rotation grid, jcal reset to 0.'/
     2 1h ,80('*'))
c10006 format('  energies of different isotopes do not match;abort')
 6783 format(' fission barriers from mebel subroutine')
 6781 format(//,10x,'fission barriers and rotating ground state energies
     1 calculated via rotating finite range model of sierk'//)
 6782 format(//,10x,'fission barriers and rotating ground state energies
     1 calculated from rotating liquid drop model'//)
c8788 format(1i6,1i1)
 7789 format(' monte carlo precompound option selected;no. events = '
     1 ,i10)
      if(ipot.gt.0)write(17,7780)
 7780 format(' monte carlo partial densities calculated for infinite hol
     1e depth')
12011 format(1i6,1f5.2)
 9877 format(' fission barriers user supplied and scaled pro
     1portionally with rldm vs. j')
 9878 format(' fission barriers user supplied and are not j dependent')
c77777format(//,' weighting option selected;abundance this isotope = ',
c    1f10.5)
80002 format(//,10x,  ' parameters selected internally under default opt
     1ion: mc=10,mp=3,inver=2,ed=1.,ike=4')
c6001 format(20a4)
c6003 format(5e12.4)
 3333 format(2x,i2,1x,1i3,1a4,1x,1f10.6,2x,1f8.6)
c3334 format(2x,i2,1x,1a4,1x,1i3,1f8.6)
c3334 format('At#=',1i3,'Element=',1a4,'Mass no.=',1i3,'Abundance=',1e9.
c    16)
  694 format(11x,10(1x,1pe9.2))
  981 format(2x,I4,2(1PE11.4),2x,I4,1PE11.4)
c 625 format(6e10.3)
c 630 format(2f5.1,3i1,i2,f3.0,i2,2i5,4f5.1,f9.0,i1,2f5.1,i5,f5.0)
  635 format(31x,'excitation energy of compound nucleus =',f6.1,'mev  j=
     1',i4/)
  640 format(2h  )
c23456789012345678901234567890123456789012345678901234567890123456789012
  645 format(1h1,33x,'compound nucleus atomic number = ',i3,' mass num',
     1'ber= ',i3/)
c 650 format (32x,'compound nucleus (equilibrium) cross section =',
c    1 e10.3/)
c 655 format(6f5.1,i1,f4.1,i1,i4,3i5,i1,i4,f4.1,i1,i5,f5.1,i1,i4)
c 660 format(50x,'evaporation code alice/livermore/91',//30x,' this code
c    1 is  recursive version 3/17/98')
  670 format(31h total fission cross section = ,e10.3)
c 675 format(48x,'gamf/gamtot = ',e10.3)
  685 format (41x,'total reaction cross section = ',f8.1/)
c 690 format(' das ist alice das ist alice das ist alice das ist alice')
  701 format (44x,'reaction cross section',f8.1/)
c9876 format(11f5.1)
 7010 format(' nz exceeds dimensioned limit. default to nz=9')
 7015 format(' na exceeds dimensioned limit. default to 22')
c7005 format(20a4)
c7465 format(/34x,'qval =',f5.1,'  ap =',f5.1,'  at =',f5.1,'  zp =',f5.
c    11,2x,'zt =',f5.1/1x,' cld =',f5.1,'  na =',i3,'  nz =',i3,'  mc =
c    2',i3,'  mp =',i3,'  inver =',i3,'  ike =',i3,'   ipch =',i3
c    3 ,'   pld =',f5.1,'  kplt=',i3,'  m3=',i3)
 7466 format(40x,'ldopt = ',i1)
c7475 format(1h1,///20a4///)
          stop
c      go to 1
      end
c-------------------------------------------------------------------
      subroutine screenin
c      this is screen input routine for code alice 
      common/gamcl/igamcl
      common/endvparm/iend,lll(6,6,6),ll(6,6)
       common/hidd/hiddcs
       common/rcsrec/rcsorig(100)
      common/outlim/limout,lim2
      common/rotor/irotor,iidummy,arote(15,25,181)
      common/ismrdat/isoid(106,260),levnm(1000),bevj(1000,3)
      common/ismrdat1/englev(1000,3),spingnd(106,260)
      common/ismrdat2/ididl(106,260),nolevs(106,260)
      common/isomer1/xgs(690),xiso(690),x2iso(690),isonum,isoint,iso2num
      common/isomer2/zs(690),as(690),eiso(690),eiso2(690)
      common/isomer3/iziso(690),iaiso(690),csgs(690),csis(690)
     1,csgse(690,100),csise(690,100),csis2e(690,100),rm1(690,100),
     2rm2(690,100)
      common/ijcpl/ajtarg,ajin
      common/tagspin/trgspin
      common/tgspin/levno
      common/gate1/spind(5,180),popj(15,24,180),egatl(5),egatu(5),npart(
     15),igate
      common/eferm/efermt
        common/mccom/ncount,ipot,rum(15,24),ppmc(15,24,999),cf
      common/nuopt/kinem,iem,kb
       common/scrn/att(10),cldd(10),bexp(24,15),reacss(100),edelt
       common/scrn1/ingo,neng,nengg,ireact
      common/ist/isot,iikl,knz,idelm,fracts,fract,ccrs(32,13,100)
     1,engbuf(100),exd,een(8,999,100)
      common/shft/mc,mp,inver,ike,ipch,kq,qval,ap,at,zp,zt,cld,barfac
      common/sf/m3,kplt
      common/levop/ldopt,igam
      common/shft4/k3,jcal
      common/shft5/fiss
      common/shft2/fs(25),dsp(25),brr(25),der(25),er(25)
     1,cx(25),bilsum,xmiss,sumiz
      REAL*8 fiss
      common/nhy/ij,jl,ji,jj,td,ex1,ex2,tmx,av,gav,cost,b(3)
      common/pl5/na,nz,title(20)
      common/incr/ed,tened,efinal
      common/bound/neut2,iprot,idnuc(100),idm(100),einc(100),
     1thett(100),phit(100),ppp
      common/sft9/k5,jmax,pld,c(8),delt,alt
      common/eqin/eind(100)
c     common/tang3/efermp,nexc,eincc
      common/tang3/efermp,eincc,nexc
      common/pl3/jx,nepr,eb(100),rcsp(100),zee,amass,plex
      common/clustr/m3set
      common/fermst/aferm
       inchek=0
       write(*,*)
       write(*,*)'  This code version requires library files phtlib,'
      write(*,*)'   tmadland,MASS.TBL,SHELL.TBL,pmtable,gdrparms for' 
       write(*,*)'  execution'
       write(*,*)
    1  write(*,*)'  Enter a title line to identify this problem '
       read(*,9)title
       write(15,9)title
        iend=0
        mckeva=4
        kinem=1
        inver=0
        ireact=0
        hiddcs=1
        levno=0
        ldopt=-1
        limout=0
       isot=1
       pld=9.
       m3set=3
   9   format(20a4)
       write(*,*)
       write(*,*)'    Enter 4 values on the next line :'
       write(*,*)
       write(*,*)' ZP   the projectile atomic no.'
       write(*,*)
       write(*,*)' AP    "      "      mass no.'
       write(*,*)
       write(*,*)' ZT   the target atomic no.'
       write(*,*)
       write(*,*)' AT    "   "  mass no.(use 1000 for isotopic targets)'
       write(*,*)
       write(*,*)' note that zp=ap=0 is interpreted as a photon'
       write(*,*)
c      write(*,*)'LEVNO the number of target levels,1 for ground only,2'
c      write(*,*)'      to calculate m1 as target,3 for m2 as target'
       write(*,*)
       read(*,*)zp,ap,zt,at
       write(15,*)zp,ap,zt,at
       if(ap.gt.at)write(*,*)'Inverse kinematics (AP>AT);precompound
     1 error will result; rerun with AP<AT'
       if(isot.eq.1)att(1)=at
       write(*,*)
       write(*,*)' M3SET enter 3 to use n,p,4He as ejectiles,or 7 to '
       write(*,*)'       also use d,t,3He and 7Be'
       write(*,*)
        read(*,*)m3set
        write(15,*)m3set
       write(*,*)
       write(*,*)'       make two entries on the next line: '
       write(*,*)
       write(*,*)'  ED   the ejectile channel energy histogram width'
       write(*,*)'       (this should be an integer increment of 0.1 '
       write(*,*)'       MeV;e.g. 0.1,0.5,etc. If you enter 0, default'
       write(*,*)'       will be 0.5 MeV'
       write(*,*)
       write(*,*)'  NCOUNT  Number of cascades to run for MC option'
       write(*,*)'          if "0" is entered,default is 100000 events'
       write(*,*)
       read(*,*)ed,ncount
       if(ncount.eq.0)ncount=100000
       if(ed.eq.0.)ed=0.5
       write(15,*)ed,ncount
       write(*,*)' The next entries control projectile energies:there'
       write(*,*)' are TWO possible modes;either a list of energies to'
       write(*,*)' be used,or a loop of "NENG" energies,for which '
       write(*,*)' the first energy "EINC" is given(lab frame),followed'
       write(*,*)' by the increment in energy, EDELT.The choice is made'
       write(*,*)' by entering'
       write(*,*)'                                    '
      write(*,*)' "NENG"  as a negative number (list of energies),or as'
      write(*,*)'         a positive number(loop on EDELT).'
       write(*,*)'                        '
       write(*,*)'        In either case,there is a limit of 100'
       write(*,*)'        energies per problem'
       write(*,*)
      write(*,*)'         Enter  NENG  now'
       read(*,*)neng
       write(15,*)neng
       nengg=abs(neng)
c 7/17/2012 put test to stop entry of over 100 incident energies
                       if(nengg.gt.100)then
       write(*,*)'Number of energies requested exceeds 100 limit'
       write(*,*)' please resubmit as two or more runs of 100 or fewer
     1energies per run'
       go to 1
                          endif
       if(neng.lt.0)then
       write(*,*)'         Enter a list of ',nengg,' incident'
       write(*,*)'         energies EINC '
       endif
       if(neng.gt.0)write(*,*)'  EINC   Enter the first energy EINC '
       if(neng.gt.1)write(*,*)'  EDELT  and the loop increment,now'
       if(neng.lt.0)read(*,*)(eind(i),i=1,nengg)
       if(neng.lt.0)write(15,*)(eind(i),i=1,nengg)
       if(neng.gt.1)read(*,*)ein,edelt
       if(neng.eq.1)read(*,*)ein
       if(neng.gt.1)write(15,*)ein,edelt
       if(neng.eq.1)write(15,*)ein
       if(neng.gt.0)eind(1)=ein
      write(*,*)' To begin execution of this problem,       ENTER "777"'
      write(*,*)' To re-enter data,                         ENTER "111"'
      write(*,*)' To view options to change default parms,  ENTER "333"'
      write(*,*)'Note: if interested in sub actinide fission and you '
      write(*,*)' get none, re-run, use 333 option and force ldopt=1,'
      write(*,*)' the Fermi gas level density option'
      read(*,*) istart
      write(15,*)istart
      if(istart.eq.777.and.ap.gt.1)kinem=0
      if(istart.eq.111)go to 1
      if(istart.eq.333)go to 152
      if(istart.eq.777)return
      write(*,*)' You must enter either 777,333,or 111 to continue'
      read(*,*)istart
      write(15,*)istart
      if(istart.eq.777.and.ap.gt.1)kinem=0
      if(istart.eq.777)return
      if(istart.eq.111)go to 1
      if(istart.eq.333)go to 152
 152  continue
      write(*,*)
      write(*,*)' Default options may be changed by inputting a number'
      write(*,*)' which is a sum of "PARM" integers as follows:'
      write(*,*)
      write(*,*)' To specify mass no. at or below which Fermi'
      write(*,*)' breakup is implemented (default=12)          PARM=512'
      write(*,*)
      write(*,*)' for ENDF output                              PARM=256'
       write(*,*)
c     write(*,*)' for emission of d,t,3He and 7Be in addition ' 
c     write(*,*)' to default n,p,4He (precompound+evaporation),PARM=128'
       write(*,*)
      write(*,*)' for user provided reaction cross sections,   PARM= 64'
       write(*,*)
      write(*,*)' for inverse cross section options,           PARM= 32'
       write(*,*)
      write(*,*)' ALICE selects level densities as Fermi gas or'
      write(*,*)' Kataria-Ramamurthy based on shell proximity;to'
      write(*,*)' force a choice, including Obninsk densities, PARM= 16'
       write(*,*)
      write(*,*)' for an isomeric,rather than gstate target    PARM=  8'
       write(*,*)
      write(*,*)' to use the s-wave angular momentum approx    PARM=  4'
       write(*,*)
      if(ap.gt.1.)then
      write(*,*)' for non-default heavy ion projectile options '
      write(*,*)' (Fermi coupling,Fermi energy values),        PARM=  2'
       write(*,*)' coupling, Fermi energies, PARM= 2'
      write(*,*)' note precompound clusters not available with'
      write(*,*)' PARM=2 options (presently)'
      endif
       write(*,*)
      write(*,*)' To limit output to yields /suppress spectra, PARM=  1'
       write(*,*)
      write(*,*)' Enter the sum of all PARM values (as integer)for '
      write(*,*)' options you wish to select ' 
       write(*,*)
      read(*,*)ipar
      write(15,*)ipar 
       aferm=12.
       if(ipar.ge.512)then
       write(*,*)'type the mass number at which Fermi breakup begins'
       read(*,*)aferm
       ipar=ipar-512
       endif
c      if(ipar.ge.256)iend=1
       if(ipar.ge.256)then
       write(*,*)' Enter "1" for CM ENDF output,"2" for Lab system'
       read(*,*)iend
       write(15,*)iend
       if(ipar.ge.256.and.iend.ne.1.and.iend.ne.2)then
       write(*,*)' You MUST enter a 1 or 2 to direct format of iend out
     1put; enter now!'
       read(*,*)iend
       write(15,*)iend
        endif
        endif

       if(ipar.ge.256)ipar=ipar-256
c      if(ipar.ge.128)m3set=7
       if(ipar.ge.128)ipar=ipar-128
       if(m3set.eq.7.and.ap.eq.0.)then
         write(*,*)' You have selected cluster emission for a photonucle
     1ar reaction.'
      write(*,*)' To treat the incident photon as a nucleon for cluster
     1excitation, '
        write(*,*)' enter "0",'
        write(*,*)' or, to use only the secondary nucleon cascade,'
        write(*,*)' enter "1"' 
         read(*,*)igamcl
         write(15,*)igamcl
                                  endif
c set m3set=7 to evaporate d,t,3He,7Be
       if(ipar.ge.64)ireact=1
       if(ireact.gt.0.and.neng.lt.0)then
       write(*,*)'    Now enter ',nengg,' reaction cross sections in '
       write(*,*)'    the same order as the incident energies'
       endif

       if(ireact.gt.0.and.neng.ge.0)then
      write(*,*)' Enter ',neng,'reaction cross sections,one for each'
      write(*,*)' incident energy,starting with EINC'
      endif
      if(ireact.gt.0)then
      write(*,*)' Enter', nengg,' values of reaction cross sections NOW'
      read(*,*)(reacss(i),i=1,nengg)
      write(15,*)(reacss(i),i=1,nengg)
      endif
       if(ipar.ge.64)ipar=ipar-64
       if(ipar.ge.32)then
       write(*,*)
       write(*,*)'  INVER  inverse reaction cross section option'
       write(*,*)'         "0" selects the nuclear optical model;'
       write(*,*)'         "2" selects the classical sharp cutoff '
       write(*,*)'         (faster but less accurate)'
       write(*,*)
       write(*,*)'           Enter a value now  '
       read(*,*)inver
       write(15,*)inver
       write(*,*)
        ipar=ipar-32
        endif
        if(ipar.ge.16)then
        write(*,*)' Default is to select Fermi gas or Kataria level'
        write(*,*)' densities internally. Other options are as follows:'
        write(*,*) 
       write(*,*)'      The next line has two entries:'
       write(*,*)
       write(*,*)'  LDOPT     level density option:  '
       write(*,*)'         0 Fermi gas,backshifted pairing energies'
       write(*,*)'         1 Kataria-Ramamurthy'
       write(*,*)'         2 Obninsk (note fission not working with thi
     1s option)'
        ams=ap+at
        zms=zp+zt
        iin=ams-zms+.0001
        iiz=zms+0.0001
       write(*,*)'         3 Gilbert -Cameron, not yet tested'
         if(iin.lt.10.and.iiz.lt.10)then
       write(*,*)'         4 Mark Chadwick densities for light (Z<10 nuc
     1lei'
       endif
       write(*,*)
       write(*,*)'  PLD      enter FG level density parameter here,or'
       write(*,*)'           "0" for default "9",i.e. a=A/9'
       write(*,*)
       write(*,*)'           Enter these two values now'
       read(*,*)ldopt,pld
       write(15,*)ldopt,pld
       if(pld.eq.0.)pld=9.
        ipar=ipar-16
        endif
        if(ipar.ge.8)then
       write(*,*)'  Isomers may be used as targets; for ground state'
       write(*,*)'  target enter "0", for m1 "1", for m2 "2" '
       read(*,*)levno
       write(17,*)levno
       levno=levno+1
       ipar=ipar-8
       endif
       if(ipar.ge.4)then
       write(*,*)'  enter  "0" for no rotational energy correction, or'
       write(*,*)'  enter  "1" for the "s-wave approximation" '
       read(*,*) irotor
       write(15,*)irotor
       ipar=ipar-4
       endif
       if(ipar.ge.2.and.ap.gt.1.)then
       write(*,*)'  For heavy ion reactions (mass projectile gt 1.),'
       write(*,*)'  actual coupling of Fermi energies may be selected'
       write(*,*)'  rather than a faster algorithm,and the default'
       write(*,*)'  Fermi energies of target and projectile may be'
       write(*,*)'  changed from the 35 MeV values. Enter the following'
       write(*,*)'  parameters, or enter zeroes for default:'
       write(*,*)
       write(*,*)' ICOUPLE  as "1" to use the algorithm for DDCS,or"0"'
       write(*,*)'          for actual random coupling of Fermi momenta
     1'
       write(*,*)' Enter ICOUPLE now'
                 read(*,*)icouple
                 write(15,*)icouple
                  hiddcs=icouple
       write(*,*)
       write(*,*)' EFERMT,EFERMP target and projectile Fermi energies'
       write(*,*)'               Enter these values now, or "0" for '
       write(*,*)'               default of 35 MeV'
       read(*,*)efermt,efermp
       write(15,*)efermt,efermp
       ipar=ipar-2
       endif
       if(levno.eq.0)levno=1
       if(ipar.eq.1)then
         kinem=0
       write(*,*)
       write(*,*)'            enter limout as "1" '
       write(*,*)'            to limit output to yield files only,'
       write(*,*)' LIMOUT     enter "0" for full output, which will'
       write(*,*)'            include particle emission spectra'
       write(*,*)
       read(*,*)limout
       write(15,*)limout
       endif
       return

c      write(*,*)'  IGATE    Enter "0" for no J gates,or number (1-5) if
c    1 desired'
c      read(*,*)igate
c       write(15,*)igate
        go to 12
       if(igate.gt.0)then
c23456789012345678901234567890123456789012345678901234567890123456789012
      write(*,*)' EGATL,EGATU,PARTID   enter the lower and upper gate en'
     1','ergies and particle ID'
       write(*,*)'                       (1=neutron,2=proton) for each
     1gate,one set per line'
        do 800 i=1,igate
        read(*,*)egatl(i),egatu(i),npart(i)
  800   write(15,*)egatl(i),egatu(i),npart(i)
        endif
  12   continue
       write(*,*)
       write(*,*)'  ***************************************************'
       write(*,*)
       write(*,*)'     A summary of this input file follows:'
       write(*,*)
       write(*,*)'  ***************************************************'
       write(*,9)title
       write(*,*)
       write(7,*)
       write(*,*)zp,ap,zt,at,isot,levno,' ZP,AP,ZT,AT,ISOT,LEVNO'
       write(7,*)zp,ap,zt,at,isot,levno,' ZP,AP,ZT,AT,ISOT,LEVNO'
       isot=0
       read(*,*)zp,ap,zt,at
       write(15,*)zp,ap,zt,at,levno
       isot=1
       levno=1
       if(isot.eq.1)att(1)=at
c      if(isot.gt.1)write(*,*)
c     if(isot.gt.1) write(*,*)'Enter the mass numbers of the isotopes '
c     if(isot.gt.1)write(*,*)'not entered on the previous line'
c    1t entered on previous line'
c      if(isot.gt.1)read(*,*)(att(i),i=2,isot)
c      if(isot.gt.1)write(15,*)(att(i),i=2,isot)
c      if(isot.gt.1)write(*,*)
c      if(isot.gt.1)write(*,*)'Enter the PERCENT ABUNDANCE of all '
c      if(isot.gt.1)write(*,*)'target isotopes in the order entered on '
c      if(isot.gt.1)write(*,*)'the previous two lines'
c      if(isot.gt.1)read(*,*)(cldd(i),i=1,isot)
c      if(isot.gt.1)write(15,*)(cldd(i),i=1,isot)
c       if(isot.gt.1)trgspin=0.
       write(*,*)
       write(*,*)'            enter limout as "1" '
       write(*,*)'            to limit output to yield files only,'
       write(*,*)' LIMOUT     enter "0" for full output, which will'
       write(*,*)'            include particle emission spectra'
       write(*,*)
       read(*,*)limout
       write(15,*)limout
       write(*,*)
       write(*,*)
       write(*,*)
       write(*,*)'            The next input selects a group of '
       write(*,*)'            default parameters as follows:'
       write(*,*)
       write(*,*)'        level densities-Fermi gas model,backshifted'
       write(*,*)'        level density parameter-a=A/9'
       write(*,*)'        ejectiles treated,neutron,proton,alpha'
       write(*,*)'        reaction cross sections internally generated'
       write(*,*)'        Weisskopf-Ewing model, fission'
       write(*,*)
       write(*,*)'        To accept these defaults,enter "0" now;to'
       write(*,*)'        reject,enter "1" '
       write(*,*)' IPARM'
       write(*,*)
       read(*,*)iparm
       write(15,*)iparm
       lim2=0
       if(limout.eq.2)then
       lim2=1
       limout=1
       endif
       if(iparm.eq.0)then
       jcal=1
       ireact=0
       pld=9.
       m3set=3
       cld=1.02
       endif
       write(*,*)'            The next line has one parm to determine'
       write(*,*)'            the type of calculation:'
       write(*,*)
       write(*,*)' JCAL    "0" Bohr-Wheeler fission/evaporation (use 10'
       write(*,*)'             to calc with bf=0. if J.gt.Jcrit)'
       write(*,*)'             iers for each nuclide'
       write(*,*)
       write(*,*)
       write(*,*)'         "1" Weisskopf-Ewing evaporation,no fission'
       write(*,*)
       write(*,*)'         "2" s-wave approximation,liquid drop moment '
       write(*,*)'             of inertia'
       write(*,*)
       write(*,*)'         "3" rigid body moment of inertia'
       write(*,*)
       write(*,*)'         Enter JCAL now '
       read(*,*)jcal
       if(jcal.eq.0)then
       write(*,*)'Choose from the following fission options by entering'
       write(*,*)'JCAL once more according to the following menu:'
       write(*,*)
       write(*,*)' JCAL    "0" Bohr-Wheeler fission/evaporation (use 10'
       write(*,*)'             to calc with bf=0. if J.gt.Jcrit)'
       write(*,*)'             iers for each nuclide'
       write(*,*)
       write(*,*)'        "-1" if fission calc with user supplied barr-'
       write(*,*)'             iers for each nuclide'
       write(*,*)
       write(*,*)'        "-2" user supplied barriers scaled with J as '
       write(*,*)'             for Rotating LDM'
       write(*,*)'         Enter JCAL now '
       read(*,*)jcal
       endif
       if(jcal.lt.0)then
       write(*,*)' As selected,fission barriers will now be entered '
       write(*,*)' for each isotope'
       write(*,*)' Enter the number of mass nos. for Z compound for'
       write(*,*)' which you will enter fission barriers (NA),and the'
       write(*,*)' number of Z for which barriers will be needed (NZ)'
       read(*,*)na,nz
       write(15,*)na,nz
       write(*,*)
       write(*,*)' now enter the barriers for all isotopes of CN,CN-1'
       write(*,*)' until NA values have been entered;then repeat for '
       write(*,*)' Z-1,Acn-1,Acn-2,etc,until all A and Z are read'
         do 11 iz=1,nz
         read(*,*)(bexp(ia,iz),ia=1,na)
       write(15,*)(bexp(ia,iz),ia=1,na)
 11    continue
       endif
       if(jcal.eq.0)then
       write(*,*)'         Two parameters for fission should be entered'
       write(*,*)'         on the next line:'
       write(*,*)'  BARFAC the factor which multiplies the finite range'
       write(*,*)'         rotating liquid drop fission barrier,and'
       write(*,*)'  CLD    the ratio of single particle level densities'
 
       write(*,*)'         at saddle point to ground state'
       write(*,*)'    Both these values may be entered as 1. as default'
       endif
      if(jcal.eq.0)read(*,*)barfac,cld
      if(jcal.eq.0)write(15,*)barfac,cld
       write(*,*)'      The next line has two entries:'
       write(*,*)
       write(*,*)
       write(*,*)'  LDOPT     level density option:  '
       write(*,*)'         0 Fermi gas,backshifted pairing energies'
       write(*,*)'         1 Kataria-Ramamurthy'
       write(*,*)'         4 Fermi gas,recalculated for each mass'
       write(*,*)'         6 uses calculated liquid drop masses with '
      write(*,*)'            neither shell nor pairing terms'
       write(*,*)
       write(*,*)'  PLD      enter FG level density parameter here,or'
       write(*,*)'           "0" for default "9",i.e. a=A/9'
       write(*,*)
       write(*,*)'           Enter these two values now'
       read(*,*)ldopt,pld
       write(15,*)ldopt,pld
       if(pld.eq.0.)pld=9.
 
       return
       end
      subroutine stdinpt
c      this is  standard input routine for code alice 8/98
      common/gamcl/igamcl
      common/endvparm/iend,lll(6,6,6),ll(6,6)
       common/hidd/hiddcs
       common/rcsrec/rcsorig(100)
      common/outlim/limout,lim2
      common/rotor/irotor,iidummy,arote(15,25,181)
      common/ismrdat/isoid(106,260),levnm(1000),bevj(1000,3)
      common/ismrdat1/englev(1000,3),spingnd(106,260)
      common/ismrdat2/ididl(106,260),nolevs(106,260)
      common/isomer1/xgs(690),xiso(690),x2iso(690),isonum,isoint,iso2num
      common/isomer2/zs(690),as(690),eiso(690),eiso2(690)
      common/isomer3/iziso(690),iaiso(690),csgs(690),csis(690)
     1,csgse(690,100),csise(690,100),csis2e(690,100),rm1(690,100),
     2rm2(690,100)
      common/ijcpl/ajtarg,ajin
      common/tagspin/trgspin
      common/tgspin/levno
      common/gate1/spind(5,180),popj(15,24,180),egatl(5),egatu(5),npart(
     15),igate
      common/eferm/efermt
        common/mccom/ncount,ipot,rum(15,24),ppmc(15,24,999),cf
      common/nuopt/kinem,iem,kb
       common/scrn/att(10),cldd(10),bexp(24,15),reacss(100),edelt
       common/scrn1/ingo,neng,nengg,ireact
      common/ist/isot,iikl,knz,idelm,fracts,fract,ccrs(32,13,100)
     1,engbuf(100),exd,een(8,999,100)
      common/shft/mc,mp,inver,ike,ipch,kq,qval,ap,at,zp,zt,cld,barfac
      common/sf/m3,kplt
      common/levop/ldopt,igam
      common/shft4/k3,jcal
      common/shft5/fiss
      common/shft2/fs(25),dsp(25),brr(25),der(25),er(25)
     1,cx(25),bilsum,xmiss,sumiz
      REAL*8 fiss
      common/nhy/ij,jl,ji,jj,td,ex1,ex2,tmx,av,gav,cost,b(3)
      common/pl5/na,nz,title(20)
      common/incr/ed,tened,efinal
      common/bound/neut2,iprot,idnuc(100),idm(100),einc(100),
     1thett(100),phit(100),ppp
      common/sft9/k5,jmax,pld,c(8),delt,alt
      common/eqin/eind(100)
c     common/tang3/efermp,nexc,eincc
      common/tang3/efermp,eincc,nexc
      common/pl3/jx,nepr,eb(100),rcsp(100),zee,amass,plex
      common/clustr/m3set
      common/fermst/aferm

       inchek=0
   1   read(15,9)title
   9   format(20a4)
        mckeva=4
        kinem=1
        inver=0
        ireact=0
        hiddcs=1
        levno=1
        ldopt=-1
        limout=0
       isot=1
         pld=9.
      m3set=3
        aferm=12.
       read(15,*)zp,ap,zt,at
       if(isot.eq.1)att(1)=at
       read(15,*)m3set
       read(15,*)ed,ncount
       if(ncount.eq.0)ncount=100000
       if(ed.eq.0.)ed=0.5
       read(15,*)neng
       nengg=abs(neng)
c 7/17/2012 trap for too many incident energies
      if(nengg.gt.100)then
       write(*,*)'Number of incident energies requested exceeds 100 dime
     1nsion limit'
      write(*,*)'Please enter as more than one run,each with 100 or fewe
     1r energies'
        endif
       if(neng.lt.0)read(15,*)(eind(i),i=1,nengg)
       if(neng.gt.1)read(15,*)ein,edelt
       if(neng.eq.1)read(15,*)ein
       if(neng.gt.0)eind(1)=ein
      read(15,*) istart
      if(istart.eq.777.and.ap.gt.1)kinem=0
      if(istart.eq.111)go to 1
      if(istart.eq.333)go to 152
      if(istart.eq.777)go to 150
      read(15,*)istart
      if(istart.eq.777.and.ap.gt.1)kinem=0
      if(istart.eq.777)return
      if(istart.eq.111)go to 1
      if(istart.eq.333)go to 152
 150     continue
      if(istart.eq.777)return
 152  continue
      read(15,*)ipar
         iend=0
       aferm=12.
       if(ipar.ge.512)then
       read(*,*)aferm
       ipar=ipar-512
       endif
       if(ipar.ge.256)read(15,*)iend
        if(ipar.ge.256.and.iend.ne.1.and.iend.ne.2)read(15,*)iend
       if(ipar.ge.256)ipar=ipar-256
c      if(ipar.ge.128)m3set=7
       if(ipar.ge.128)ipar=ipar-128
       if(ipar.ge.64)ireact=1

      if(ireact.gt.0.)then
      read(15,*)(reacss(i),i=1,nengg)
      endif
       if(ipar.ge.64)ipar=ipar-64
       if(ipar.ge.32)then
       read(15,*)inver
        ipar=ipar-32
        endif
        if(ipar.ge.16)then
       read(15,*)ldopt,pld
       if(pld.eq.0.)pld=9.
        ipar=ipar-16
        endif
        if(ipar.ge.8)then
       read(15,*)levno
       levno=levno+1
       ipar=ipar-8
       endif
       if(ipar.ge.4)then
       read(15,*) irotor
       ipar=ipar-4
       endif
       if(ipar.ge.2.and.ap.gt.1.)then
                 read(15,*)icouple
                  hiddcs=icouple
       read(15,*)efermt,efermp
       ipar=ipar-2
       endif
       if(levno.eq.0)levno=1
       if(ipar.eq.1)then
        kinem=0
       read(15,*)limout
       endif
       return
       write(*,*)
       write(*,*)'  ***************************************************'
       write(*,*)
       write(*,*)'     A summary of this input file follows:'
       write(*,*)
       write(*,*)'  ***************************************************'
       write(*,9)title
       write(7,9)title
       write(*,*)
       write(7,*)
       write(*,*)zp,ap,zt,at,isot,levno,' ZP,AP,ZT,AT,ISOT,LEVNO'
       write(7,*)zp,ap,zt,at,isot,levno,' ZP,AP,ZT,AT,ISOT,LEVNO'
       write(*,*)
       write(7,*)
       if(isot.gt.1)write(*,*)(att(i),i=2,isot),' Other isotope masses'
       write(*,*)
       if(isot.gt.1)write(*,*)(cldd(i),i=1,isot),' Isotopic abundances'
       write(*,*)
       write(*,*)jcal,' JCAL'
       write(*,*)
       if(jcal.lt.0)write(*,*)na,nz,' NA,NZ'
       write(*,*)
       if(jcal.lt.0)then
       do 110 iz=1,nz
  110  write(*,*)(bexp(ia,iz),ia=1,na),' BEXP(ia,iz)'
       endif
       write(*,*)
       if(jcal.eq.0)write(*,*)barfac,cld,' BARFAC,CLD'
       write(*,*)
       write(*,*)ed,inver,ldopt,pld,' ED,INVER,LDOPT,PLD'
       write(*,*)
       write(*,*)kinem,ncount,' KINEM,NCOUNT'
       write(*,*)
       write(*,*)
       write(*,*)neng,' NENG'
       write(*,*)
      if(neng.lt.0)write(*,*)(eind(i),i=1,nengg),' EINC'
       write(*,*)
       if(neng.gt.0)write(*,*)ein,edelt,' EINC,EDELT'
       write(*,*)
       write(*,*)ireact,m3,' IREACT,M3'
       if(ireact.gt.0)then
       write(*,*)(reacss(i),i=1,nengg),'Reaction cross sections'
       endif
       return
       end
c 00000000000000000000000000000000000000000000000000000000000
c--------------------------------------------------------------------
      subroutine finish
c--------------------------------------------------------------------
      common/endvparm/iend,lll(6,6,6),ll(6,6)
      common/shft/mc,mp,inver,ike,ipch,kq,qval,ap,at,zp,zt,cld,barfac
      common/clustr/m3set
      common/outlim/limout,lim2
      common/nuopt/kinem,iem,kb
       common/scrn1/ingo,neng,nengg,ireact
      write(*,*)
      write(*,*)' Calculation is completed;the following output files',
     1' have been created:'
      write(*,*)
      write(*,*)' "INSAVE " input file created from screen input'
      write(*,*)
      write(*,*)' "OUTPT"   many things-read to check-may reduce in futu
     1re'
      write(*,*)
      write(*,*)' "PLOT  "  contains excitation function data for plots'
      write(*,*)
c     write(*,*)' "SIGS  "  summary of product yields vs. incident energ
c    1y'
c     write(*,*)
      write(*,*)' "SDCSCM"  single differential cross sections in CM fra
     1me (mb/MeV)'
      write(*,*)
      write(*,*)' "ENGANG"  energy columns versus angle DDCS in CM frame
     1(mb/MeV/sr)'
      write(*,*)
      write(*,*)' "ANGENG"  angle columns versus energy in CM frame(mb/
     1MeV/sr)'
      write(*,*)
      write(*,*)' "SDCSTOTLB" total single differential spectra,evap.'
      write(*,*)'             plus fission n,p,alpha,LAB frame'
      if(kinem.gt.0.and.limout.lt.1)then
      write(*,*)
      write(*,*)' "LABOUT"  single differential cross sections in LAB fr
     1ame (mb/MeV)'
      write(*,*)
      write(*,*)
      write(*,*)' "LABANGD" DDCS in LAB frame angle columns vs. energy'
      write(*,*)'           (mb/MeV/sr)'
      write(*,*)
      write(*,*)' "LABENGD" DDCS in LAB frame,energy columns vs. angle'
      write(*,*)'           (mb/MeV/sr)'
      write(*,*)
                   if(m3set.eq.7)then
      write(*,*)' "SDCSCL  " single differential cross section clusters'
      write(*,*)'            from evaporation process,CM frame(d,t,3He,7
     1Be)'
      write(*,*)
      write(*,*)' "SDLBCL  " single differential cross section clusters'
      write(*,*)'            from evaporation process,LAB frame(d,t,3He,
     17Be)'
      write(*,*)
                 endif
c     write(*,*)
      write(*,*)' "GAMMASPEC" gamma ray spectra summed over nuclides'
c     write(*,*)
          if(iend.eq.1)then
      write(*,*)
      write(*,*)' "ENDF"    single/ddcs CM for exclusive reactions'
      write(*,*)'           presently fully implemented only for '
      write(*,*)'           nucleon in reactions'
            endif
          if(iend.eq.2)then
      write(*,*)
      write(*,*)' "ENDF"    single/ddcs LAB for exclusive reactions'
      write(*,*)'           presently fully implemented only for '
      write(*,*)'           nucleon in reactions'
      write(*,*)
           endif
      if(ap.gt.1. )then
      write(*,*)
      write(*,*)' "EXCDST"     exciton distribution for Fermi coupling'
      write(*,*)
      write(*,*)' "EXCDST1  " exciton distribution from algorithm'
      write(*,*)
      write(*,*)' "EXCITON "  exciton disribution from formula'
      write(*,*)
      write(*,*)' "RATIOX  "  ratio of coupled exciton distribution to'
      write(*,*)'             exciton formula'
      write(*,*)
            endif
      write(*,*)
      write(*,*)' "SDCSFLB"  single diff. fission fragment emitted n,p,a
     1lpha spectra'
      write(*,*)'            in LAB frame (mb/Mev*sr)'
      write(*,*)
      write(*,*)' "SDCSFCM"  single differential particle spectra in'
      write(*,*)'            fission fragment frame (not true cm)'
      write(*,*)'             plus fission n,p,alpha,LAB frame'
      write(*,*)
      write(*,*)' "DDCSFLB"  double diff. spectra of n,p,alpha emitted'
      write(*,*)'            from fission fragments in LAB frame,A vs E'
      write(*,*)
      write(*,*)' "DDCSFCM"  double diff. fission part. spectra'
      write(*,*)'            in ff frame,energy vs. angle'
      write(*,*)
      write(*,*)' "DDCSTOTLB" total double differential spectra for'
      write(*,*)'             fission+evaporation n,p,alpha,LAB frame'
      write(*,*)
      if(kinem.gt.0)then
c     write(*,*)' "RECVEL " recoil velocity distributions'
c     write(*,*)
c     write(*,*)' "RECOILS" recoil energy,mass,charge,angular distributi
c    1ons'
c     write(*,*)
c     write(*,*)' "EJPLANES"spin dependent arrays,gates,isomer data'
      endif
      write(*,*)' "FOUT"    fission product yields-A,Excitation,Z'
      write(*,*)
      write(*,*)' "FIZISOMR" fission product yields'
      write(*,*)
      write(*,*)' "FIZISRAT" fission isomer ratios '
      write(*,*)
      write(*,*)' "FIZYBYA"  fission mass yields  '
      write(*,*)
c     write(*,*)' "FIZZERS "unity normalized output of popfiz(Z,A,E/2)'
c     write(*,*)'           file for fission mass use,see subroutine'
c     write(*,*)'           fiznorm for format details'
c     write(*,*)
c     write(*,*)' "FIZZERSJ"unity normalized output of pofizbig(Z,A,E/2
c    1,J)'
c     write(*,*)
c     write(*,*)' "FIZEJ"   fission E-J planes'      
      write(*,*)' "FRAGSS"  fission fragment distributions:mass vs Z,'  
      write(*,*)'           and Z yields for fixed A'
      write(*,*)
      write(*,*)' "FIZYLDS" fission fragment yields with normalization'
      write(*,*)
                          endif
      write(*,*)'            All spectra in CM system are given with'
      write(*,*)'            respect to channel angle for each emission'
      write(*,*)'            (not referred back to beam direction)' 
      return
      end
c----------------------------------------------------------------------------------
c00000000000000000000000000000000000000000000000000000000
c 7/7/14 prelim to calc initial arrays-----
      subroutine prelim
      common/ismrdat1/englev(1000,3),spingnd(106,260)
      common/ismrdat2/ididl(106,260),nolevs(106,260)
      common/ismrdat/isoid(106,260),levnm(1000),bevj(1000,3)
      common/isomer1/xgs(690),xiso(690),x2iso(690),isonum,isoint,iso2num
      common/isomer2/zs(690),as(690),eiso(690),eiso2(690)
c++++++++++++++++++++++++s2+++++++++++++++
       do 4000 iz=5,100
       do 4000 in=2,160
       lz=iz
       la=in+iz
       ma=iz+in
       neut=la-lz
       if(neut.lt.1)go to 4000
       bn=float(neut)
       bz=float(lz)
       neut2=neut/2
       jz2=lz/2
       endff=bn/2.-neut2
       ezdff=bz/2.-jz2
       non=endff*2
       nez=2.*ezdff
        itot=non+nez
       if(itot.eq.1)nolevs(iz,ma)=1
       n=isoid(lz,la)
       if(levnm(n).lt.3)go to 4000
         zs(n)=lz
         as(n)=la
         xgs(n)=bevj(n,1)
         xiso(n)=bevj(n,2)
         eiso(n)=englev(n,2)
          x2iso(n)=bevj(n,3)
          eiso2(n)=englev(n,3)
 4000  continue
       do 4501 iz=5,100
       do 4501 in=2,160
       lz=iz
       la=in+iz
       n=isoid(lz,la)
       if(levnm(n).ne.2)go to 4501
         zs(n)=lz
         as(n)=la
         xgs(n)=bevj(n,1)
         xiso(n)=bevj(n,2)
         eiso(n)=englev(n,2)
 4501  continue
        return
         end
        subroutine tagout
      common/ist/isot,iikl,knz,idelm,fracts,fract,ccrs(32,13,100)
     1,engbuf(100),exd,een(8,999,100)
      common/shft/mc,mp,inver,ike,ipch,kq,qval,ap,at,zp,zt,cld,barfac
      common/pl5/na,nz,title(20)
      common/tgspin/levno
 7475 format(1h1,///20a4///)
      write (17,7475)title
crC
       write(8,7475)title
       write(9,7475)title
       write(10,7475)title
       write(11,7475)title
       write(58,7475)title
       write(54,7475)title
      write(12,7475)title
c     write(13,7475)title
c     write(14,7475)title
       write(7,*)zp,ap,zt,at,isot,levno,' ZP,AP,ZT,AT,ISOT,LEVNO'
       write(8,*)zp,ap,zt,at,isot,levno,' ZP,AP,ZT,AT,ISOT,LEVNO'
       write(9,*)zp,ap,zt,at,isot,levno,' ZP,AP,ZT,AT,ISOT,LEVNO'
       write(10,*)zp,ap,zt,at,isot,levno,' ZP,AP,ZT,AT,ISOT,LEVNO'
       write(11,*)zp,ap,zt,at,isot,levno,' ZP,AP,ZT,AT,ISOT,LEVNO'
       write(58,*)'      ZP            AP            ZT            AT  
     1             ISOT            LEVNO'
       write(58,*)zp,ap,zt,at,isot,levno
       write(12,*)zp,ap,zt,at,isot,levno,' ZP,AP,ZT,AT,ISOT,LEVNO'
c      write(13,*)zp,ap,zt,at,isot,levno,' ZP,AP,ZT,AT,ISOT,LEVNO'
c      write(14,*)zp,ap,zt,at,isot,levno,' ZP,AP,ZT,AT,ISOT,LEVNO'
       write(18,*)zp,ap,zt,at,isot,levno,' ZP,AP,ZT,AT,ISOT,LEVNO'
       write(19,*)zp,ap,zt,at,isot,levno,' ZP,AP,ZT,AT,ISOT,LEVNO'
       write(20,*)zp,ap,zt,at,isot,levno,' ZP,AP,ZT,AT,ISOT,LEVNO'
       write(21,*)zp,ap,zt,at,isot,levno,' ZP,AP,ZT,AT,ISOT,LEVNO'
       write(22,*)zp,ap,zt,at,isot,levno,' ZP,AP,ZT,AT,ISOT,LEVNO'
       write(23,*)zp,ap,zt,at,isot,levno,' ZP,AP,ZT,AT,ISOT,LEVNO'
       write(24,*)zp,ap,zt,at,isot,levno,' ZP,AP,ZT,AT,ISOT,LEVNO'
       write(25,*)zp,ap,zt,at,isot,levno,' ZP,AP,ZT,AT,ISOT,LEVNO'
       write(26,*)zp,ap,zt,at,isot,levno,' ZP,AP,ZT,AT,ISOT,LEVNO'
       write(27,*)zp,ap,zt,at,isot,levno,' ZP,AP,ZT,AT,ISOT,LEVNO'
       write(29,*)zp,ap,zt,at,isot,levno,' ZP,AP,ZT,AT,ISOT,LEVNO'
       write(37,*)zp,ap,zt,at,isot,levno,' ZP,AP,ZT,AT,ISOT,LEVNO'
       return
       end
      subroutine isotagr(isotg)
      common/isotg/iztag,imtag,tix,abundc,unct
       if(isotg.eq.2)then
           write(7,3333)iztag,imtag,tix,abundc,unct
           write(17,3333)iztag,imtag,tix,abundc,unct
           write(58,3333)iztag,imtag,tix,abundc,unct
                     endif
       if(isotg.eq.1)then
      write(7,*)
      write(7,*)' An isotopic target is assumed:'
       write(7,*)' At# Mass# % abundance % uncertainty'
      write(17,*)
      write(17,*)' An isotopic target is assumed:'
       write(17,*)' At# Mass# % abundance % uncertainty'
      write(58,*)
      write(58,*)' An isotopic target is assumed:'
       write(58,*)' At# Mass# % abundance % uncertainty'
                      endif
3333  format(2x,i2,1x,1i3,1a4,1x,1f10.6,2x,1f8.6)
      return
      end
c----------------------------------------------------------
crc
      subroutine headers(ikl)
      common/update/month,iday,iyear,nen
      common/outlim/limout,lim2
      common/eqin/eind(100)
c23456789012345678901234567890123456789012345678901234567890123456789012
c     write(13,*) ' code hms.f last modified ',month,iday,iyear
c     write(14,*) ' code hms.f last modified ',month,iday,iyear
       if(ikl.eq.6)then
      write(17,*) ' code hms.f last modified ',month,iday,iyear
      write(21,*) ' code hms.f last modified ',month,iday,iyear
      write(22,*) ' code hms.f last modified ',month,iday,iyear
      write(23,*) ' code hms.f last modified ',month,iday,iyear
      write(24,*) ' code hms.f last modified ',month,iday,iyear
      write(25,*) ' code hms.f last modified ',month,iday,iyear
      write(26,*) ' code hms.f last modified ',month,iday,iyear
      write(27,*) ' code hms.f last modified ',month,iday,iyear
      write(29,*) ' code hms.f last modified ',month,iday,iyear
      write(37,*) ' code hms.f last modified ',month,iday,iyear
      if(limout.lt.1)then
      write(8,*) ' code hms.f last modified ',month,iday,iyear
      write(10,*) ' code hms.f last modified ',month,iday,iyear
      write(11,*) ' code hms.f last modified ',month,iday,iyear
      write(58,*) ' code hms.f last modified ',month,iday,iyear
      write(12,*) ' code hms.f last modified ',month,iday,iyear
      write(18,*) ' code hms.f last modified ',month,iday,iyear
      write(19,*) ' code hms.f last modified ',month,iday,iyear
      write(20,*) ' code hms.f last modified ',month,iday,iyear
      write(39,*) ' code hms.f last modified ',month,iday,iyear
      write(40,*) ' code hms.f last modified ',month,iday,iyear
      write(41,*) ' code hms.f last modified ',month,iday,iyear
      write(42,*) ' code hms.f last modified ',month,iday,iyear
      write(43,*) ' code hms.f last modified ',month,iday,iyear
      write(44,*) ' code hms.f last modified ',month,iday,iyear
      write(45,*) ' code hms.f last modified ',month,iday,iyear
      write(50,*) ' code hms.f last modified ',month,iday,iyear
      write(51,*) ' code hms.f last modified ',month,iday,iyear
      write(52,*) ' code hms.f last modified ',month,iday,iyear
      write(53,*) ' code hms.f last modified ',month,iday,iyear
      write(54,*) ' code hms.f last modified ',month,iday,iyear
      write(54,*) ' CAUTION_ AS OF 11 OCT"07,THIS GAMMA RAY OUTPUT HAS'
      write(54,*)' NOT BEEN BENCHMARKED AND MAY BE NONSENSE!'
      write(55,*) ' code a104.f last modified ',month,iday,iyear
      endif
      write(7,*)
      write(7,*)' Problems  with code should be communicated to'
      write(7,*)' M.Blann c/o hmblann@gmail.com'
      write(7,*)' Those wishing to receive a copy of this code '
      write(7,*)' should request same from M.Blann or others'
      write(7,*)
                       endif
c
              if(ikl.eq.7)then
c
             if(limout.lt.1)then
       write( 8,*)
       write( 8,*)'Beam energy (MeV) ',eind(nen)
       write( 8,*)
       write(10,*)
       write(10,*)'Beam energy (MeV) ',eind(nen)
       write(10,*)
       write(11,*)
       write(11,*)'Beam energy (MeV) ',eind(nen)
       write(11,*)
       write(58,*)
       write(58,*)'Beam energy (MeV) ',eind(nen)
       write(58,*)
       write(12,*)
       write(12,*)'Beam energy (MeV) ',eind(nen)
       write(12,*)
c      write(13,*)
       write(16,*)
       write(16,*)'Beam energy (MeV) ',eind(nen)
       write(16,*)
c      write(13,*)'Beam energy (MeV) ',eind(nen)
c      write(13,*)
c      write(14,*)
c      write(14,*)'Beam energy (MeV) ',eind(nen)
c      write(14,*)
       write(18,*)
       write(18,*)'Beam energy (MeV) ',eind(nen)
       write(18,*)
       write(19,*)
       write(19,*)'Beam energy (MeV) ',eind(nen)
       write(19,*)
       write(20,*)
       write(20,*)'Beam energy (MeV) ',eind(nen)
       write(20,*)
c      write(21,*)
c      write(21,*)'Beam energy (MeV) ',eind(nen)
c      write(21,*)
c      write(22,*)
c      write(22,*)'Beam energy (MeV) ',eind(nen)
c      write(22,*)
       write(17,*)
       write(17,*)'Beam energy (MeV) ',eind(nen)
       write(17,*)
       write(23,*)
       write(23,*)'Beam energy (MeV) ',eind(nen)
       write(23,*)
       write(24,*)
       write(24,*)'Beam energy (MeV) ',eind(nen)
       write(24,*)
       write(25,*)
       write(25,*)'Beam energy (MeV) ',eind(nen)
       write(25,*)
       write(26,*)
       write(26,*)'Beam energy (MeV) ',eind(nen)
       write(26,*)
c      write(27,*)
c      write(27,*)'Beam energy (MeV) ',eind(nen)
c      write(27,*)
       write(29,*)
       write(29,*)'Beam energy (MeV) ',eind(nen)
       write(29,*)
c      write(37,*)
c      write(37,*)'Beam energy (MeV) ',eind(nen)
c      write(37,*)
       write(39,*)
       write(39,*)'Beam energy (MeV) ',eind(nen)
       write(39,*)
       write(40,*)
       write(40,*)'Beam energy (MeV) ',eind(nen)
       write(40,*)
       write(41,*)
       write(41,*)'Beam energy (MeV) ',eind(nen)
       write(41,*)
       write(42,*)
       write(42,*)'Beam energy (MeV) ',eind(nen)
       write(42,*)
       write(43,*)
       write(43,*)'Beam energy (MeV) ',eind(nen)
       write(43,*)
       write(44,*)
       write(44,*)'Beam energy (MeV) ',eind(nen)
       write(44,*)
       write(45,*)
       write(45,*)'Beam energy (MeV) ',eind(nen)
       write(45,*)
       write(54,*)
       write(54,*)'Beam energy (MeV) ',eind(nen)
       write(54,*)
       write(55,*)
       write(55,*)'Beam energy (MeV) ',eind(nen)
       write(55,*)
                     endif
                      endif
       return
        end
      subroutine zero(ikl)
      common/astp3/niso,no,iadlt,isotag,massno(12),mass1(12),ncsave
      common/endvarray/endvcm1(36,999,8),endvcm2(6,36,999,8),
     1endvcm3(10,36,999,8),cmen1(8,999),cmen2(6,8,999),cmen3(10,8,999)
      common/endvparm/iend,lll(6,6,6),ll(6,6)
      common/memof/sfcrs(999,8),dfsgp(36,999,4),efn(4,999)
      common/fisn/dfslb(36,999,4),cflab(999,4),dfetl(36,999,4)
      common/fissngl/dissa(200),dissu(200),dissz(200)
       common/rcsrec/rcsorig(100)
      doubleprecision dissa,dissu,dissz
      common/libst/a53t(260)
      common/fms/sdissa,sdissz
      doubleprecision sdissa,sdissz
      common/emin/eminf(106,160)
      doubleprecision dafr1,dafr2,dzfr1,dzfr2,dufr1,dufr2
      common/azfrag/dafr1,dafr2,dzfr1,dzfr2,dufr1,dufr2
      common/krld/powkr(9999),gamkr(9999),powgkr(9999)
      REAL*8 powkr,gamkr,powgkr
      common/pmasb/pmass(100,150,30)
      common/fizzc/popfiz(13,20,150),pofizbig(13,20,150,30)
      common/gams/gamgam,csgam,csgamm(100)
      REAL*8 gamgam,csgamm,csgam
      common/outlim/limout,lim2
       common/fizz/fizbar(106,260),fizcs(100)
      common/fermst/aferm
      common/ismrdat/isoid(106,260),levnm(1000),bevj(1000,3)
      common/ismrdat1/englev(1000,3),spingnd(106,260)
      common/ismrdat2/ididl(106,260),nolevs(106,260)
      common/rotor/irotor,iidummy,arote(15,25,181)
       common/scrn/att(10),cldd(10),bexp(24,15),reacss(100),edelt
       common/scrn1/ingo,neng,nengg,ireact
        common/mccom/ncount,ipot,rum(15,24),ppmc(15,24,999),cf
      common/giani/pre, sigt(1100),sigam(999),sigpre(999),gsp(999)
     2,eqgam(999)
      common/rcc/rcsi(10,100)
      common/qq/gow(1100)
      REAL*8 gow
      REAL*8 q,gamft
      common/m/ q(999),sp(999),sifis(999),gamft(999),t(999),pairx(8)
     1,scale(24,15),rd(4)
      common/cs/crsum
      common/par3/eq,sigml(999),acrs(999)
      common/lab10/pow(4,9999),gam(9999)
        REAL*8 pow,gam
      common/lab12/pof(4,9999)
      REAL*8 pof
      common/lab11/powaz(15,24,2200)
        REAL*8 powaz
      common/scat/jz,ja, probxx(2,2),q1(8,1100),fcs,eloss,fcspec
      common/isomer1/xgs(690),xiso(690),x2iso(690),isonum,isoint,iso2num
      common/isomer2/zs(690),as(690),eiso(690),eiso2(690)
      common/isomer3/iziso(690),iaiso(690),csgs(690),csis(690)
     1,csgse(690,100),csise(690,100),csis2e(690,100),rm1(690,100),
     2rm2(690,100)
      common/error/rerror(15,24)
      common/evapmc/jcount,ipremc,eminn(15,24)
      common/sf6/aoz,zoz,en(8,999)
      common/shft/mc,mp,inver,ike,ipch,kq,qval,ap,at,zp,zt,cld,barfac
      common/ug/iz ,ia
      common/nhy2/gdo,bisp
      common/pl2/sux(24,13,100)
      common/s1/jfrac,jupper
      common/sf/m3,kplt
      common/sft9/k5,jmax,pld,c(8),delt,alt
      common/par2/cncss
      common/pl8/jamma(15,24),nulim(15,24)
      common/pl1/ecm(100)
      common/nam/k2
      common/fis/ifis
      common/levop/ldopt,igam
      common/hjk/jang,iq,rcss
      common/ss/sor,rr
      common/nr34/nr3,nr4,ke5,i3d,irtst,i3t2,ijkl,lq,tem(36)
      common/parfs/mz,k6,delrr(999)
      common/shft4/k3,jcal
      common/shft5/fiss
      common/shft2/fs(25),dsp(25),brr(25),der(25),er(25)
     1,cx(25),bilsum,xmiss,sumiz
      REAL*8 fiss
      common/pl4/crs(24,13)
      common/nhy/ij,jl,ji,jq,td,ex1,ex2,tmx,av,gav,cost,b(3)
      common/paro/pq,cross
      common/hyb2/pp(15,24,999)
      common/sft5/exc(18,27),xmax
      common/pl3/jx,nepr,eb(100),rcsp(100),zee,amass,plex
      common/lym1/be(18,27,8),pair(18,27),xmas(18,27),delshl(18,27),amas
     1(18,27)
c     common/bym1/symb(18,27),symbp(18,27),symbs(18,27)
      common/pl7/mas
      common/lab3/sig(8,999),slav(8,999)
      common/pl5/na,nz,title(20)
      common/scr/k9
      common/iso/qpn(3),qpnc
      common/incr/ed,tened,efinal
      common/tst/test
      common/send/irfr,iadst,dlt
      common/gam/io,la,ppgam(999),gmspec(999),beta,ams
      common/eqin/eind(100)
      common/ist/isot,iikl,knz,idelm,fracts,fract,ccrs(32,13,100)
     1,engbuf(100),exd,een(8,999,100)
      common/cshel/shel(18,27,2)
      common/gr/rint(1100)
      common/deform/beta2(106,260)
      REAL*8 rint
      common/mbc/sigg(21,4,70)
      common/rc/rzero
      common/equ/equiv
      common/nuopt/kinem,iem,kb
      common/bound/neut2,iprot,idnuc(100),idm(100),einc(100),
     1thett(100),phit(100),ppp
      common/spins/rote(15,24,181),popej(15,24,180,100),fcx(180),ai,ecmp
      common/sping/pfpej(600,180,100),pfpj(600,180)
      common/gate1/spind(5,180),popj(15,24,180),egatl(5),egatu(5),npart(
     15),igate
       common/spins3/ajmax
       common/spins4/eprot,ealpha,eoff(8)
      common/tagspin/trgspin
      common/tgspin/levno
      common/estar/estarcnt,e6,e7,e8
      common/sigsav/sigdr(100),sigqd(100)
      common/frags/xmasf(160,106),bq(106,160,4),pf(106,160)
     1,pfp(106,160,100),ppfi(106,160),ppfz(106),ppfa(260)
     2,ppfaz(106,260)
      common/fissiso/csgsef(690,100),csisef(690,100),csis2ef(690,100)
      common/frags2/rotef(260,181),arotef(260,181)
      common/new/fpyismr(3,100)
      common/update/month,iday,iyear,nen
      common/sfmp1/sfmp(100,300,5)
      common/memo/scrs(999,8)
c---------------------------------
       if(ikl.eq.1) then
      do  i=1,200
      dissa(i)=0.
      dissu(i)=0.
      if(i.le.100) dissz(i)=0.
      enddo
      do i=1,3
      do j=1,100
       fpyismr(i,j)=0.
      enddo
      enddo
      do i=1,690
      do j=1,100
            csisef(i,j)=0.
             csis2ef(i,j)=0.
              csgsef(i,j)=0.
      enddo
       enddo
      do 2008 i=1,10
      do 2007 j=1,100
      rcsi(i,j)=0.
      fizcs(j)=0.
      csgamm(j)=0.
      sigdr(j)=0.
      sigqd(j)=0.
2007  continue
2008  continue
      do 2009 i=1,106
      do 2010 j=1,260
      beta2(i,j)=0.
      fizbar(i,j)=90.
      isoid(i,j)=1
      nolevs(i,j)=0
c nolevs will be 0 for o-o or e-e nuclides,1 for o-e
2010  continue
2009  continue
      do 2011 i=1,1000
      levnm(i)=1
      bevj(i,1)=0.
      bevj(i,2)=0.
      bevj(i,3)=0.
      englev(i,1)=0.
      englev(i,2)=0.
      englev(i,3)=0.
2011  continue
      do 10001 i=1,9
      do 10002 j=1,32
      do 10002 k=1,100
      engbuf(k)=0.
      ccrs(j,i,k)=0.
10002 continue
10001 continue
      do 10003 i=1,8
      do 10004 j=1,999
      do 10004 k=1,100
      een(i,j,k)=0.
10004 continue
10003 continue
                do 10005 i=1,100
                do 10005 j=1,300
                do 10005 k=1,5
10005           sfmp(i,j,k)=-1000.
                                endif
c
        if(ikl.eq.4)then
c
      do 7040 kk=1,13
      do 7040 jj=1,24
 7040 exc(kk,jj)=0.
c     do 7150 i=1,999
c     eqgam(i)=0.
c     sigt(i)=0.
c     gmspec(i)=0.
c     sigam(i)=0.
c     sigpre(i)=0.
c     rint(i)=0.
c7150 continue
      do 7151 i=1000,1100
      rint(i)=0.
 7151 sigt(i)=0.
          endif
       if(ikl.eq.5)then
c
       do 8800 jz=1,15
       do 8800 ja=1,24
       do 8800 jl=1,180
       popj(jz,ja,jl)=0.
c      rote(jz,ja,jl)=0.
c      rote(jz,ja,jl+90)=0.
c      arote(jz,ja,jl)=0.
c      arote(jz,ja,jl+90)=0.
       do 8800  ie=1,100
       popej(jz,ja,jl,ie)=0.
 8800  continue
       do 8810 izz=1,13
       do 8810 iaa=1,20
       do 8810 ie=1,150
       popfiz(izz,iaa,ie)=0.
       do 8810 jf=1,30
       pofizbig(izz,iaa,ie,jf)=0.
 8810 continue
      do 8811 i=1,100
      do 8811 j=1,150
      do 8811 k=1,30
      pmass(i,j,k)=0.
 8811 continue
      do lz=1,106
      ppfz(lz)=0.00
      do ln=1,160
      la=lz+ln
      if(la.le.260)ppfaz(lz,la)=0.
      ppfa(ln)=0.
      ppfa(ln+100)=0.
      ppfi(lz,ln)=0.
      do ien=1,100
      pfp(lz,ln,ien)=0.
      enddo
      enddo
      enddo
      do n=1,600
      do j=1,180
      pfpj(n,j)=0.
      do ie=1,100
      pfpej(n,j,ie)=0.
      enddo
      enddo
      enddo
      do i=1,4
      do j=1,999
      cflab(j,i)=0.
      efn(i,j)=0.
      do k=1,36
      dfslb(k,j,i)=0.
      dfetl(k,j,i)=0.
      enddo
      enddo
      enddo
      do 25 i=1,8
      do 25 ke=1,999
      sfcrs(ke,i)=0.
   25 scrs(ke,i)=0.
          endif
c
                          if(ikl.eq.8.and.no.eq.1)then
      do 10008 i=1,100
      do 10007 j=1,690
      csise(j,i)=0.
      csis2e(j,i)=0.
      csgse(j,i)=0.
      rm1(j,i)=0.
      rm2(j,i)=0.
10007 continue
10008 continue
                                                 endif
c   here initialize exclusive spectra arrays each incident energy
                          if(ikl.eq.9.and.no.eq.1)then
      do j=1,999
       do k=1,8
        cmen1(k,j)=0.
      do i=1,36
         endvcm1(i,j,k)=0.
          do l=1,6
           cmen2(l,k,j)=0.
           endvcm2(l,i,j,k)=0.
            enddo
             enddo
              enddo
               enddo
        do j=1,999
         do k=1,8
          do l=1,10
          cmen3(l,k,j)=0.
       do i=1,36
            endvcm3(l,i,j,k)=0.
             enddo
              enddo
               enddo
                enddo
                                           endif
c
      return
      end
c
      subroutine ismrout
c
      common/eder/ed2,v,vhalf,vsq
      common/fissiso/csgsef(690,100),csisef(690,100),csis2ef(690,100)
      common/outlim/limout,lim2
      common/sping/pfpej(600,180,100),pfpj(600,180)
       common/fizz/fizbar(106,260),fizcs(100)
      common/fermst/aferm
      common/ismrdat/isoid(106,260),levnm(1000),bevj(1000,3)
      common/ismrdat1/englev(1000,3),spingnd(106,260)
      common/ismrdat2/ididl(106,260),nolevs(106,260)
      common/rotor/irotor,iidummy,arote(15,25,181)
       common/scrn/att(10),cldd(10),bexp(24,15),reacss(100),edelt
       common/scrn1/ingo,neng,nengg,ireact
        common/mccom/ncount,ipot,rum(15,24),ppmc(15,24,999),cf
      common/giani/pre, sigt(1100),sigam(999),sigpre(999),gsp(999)
     2,eqgam(999)
      common/rcc/rcsi(10,100)
      common/qq/gow(1100)
      REAL*8 gow
      REAL*8 q,gamft
c     common/exclus/ gib(2,999),pr(3),ge(2100),zz(3),
c    1delta(3),ts(3),pairx(8),asum(10),rate(2,999)
c     common/m/ q(999),sp(999),sifis(999),gamft(999),t(999)
c    1,scale(24,15),rd(4)
      common/m/ q(999),sp(999),sifis(999),gamft(999),t(999),pairx(8)
     1,scale(24,15),rd(4)
      common/cs/crsum
      common/par3/eq,sigml(999),acrs(999)
      common/lab10/pow(4,9999),gam(9999)
        REAL*8 pow,gam
      common/lab11/powaz(15,24,2200)
      common/scat/jz,ja, probxx(2,2),q1(8,1100),fcs,eloss,fcspec
      common/isomer1/xgs(690),xiso(690),x2iso(690),isonum,isoint,iso2num
      common/isomer2/zs(690),as(690),eiso(690),eiso2(690)
      common/isomer3/iziso(690),iaiso(690),csgs(690),csis(690)
     1,csgse(690,100),csise(690,100),csis2e(690,100),rm1(690,100),
     2rm2(690,100)
      common/error/rerror(15,24)
      common/evapmc/jcount,ipremc,eminn(15,24)
      common/sf6/aoz,zoz,en(8,999)
      common/shft/mc,mp,inver,ike,ipch,kq,qval,ap,at,zp,zt,cld,barfac
      common/ug/iz ,ia
      common/nhy2/gdo,bisp
      common/pl2/sux(24,13,100)
      common/s1/jfrac,jupper
      common/sf/m3,kplt
      common/sft9/k5,jmax,pld,c(8),delt,alt
      common/par2/cncss
      common/pl8/jamma(15,24),nulim(15,24)
      common/pl1/ecm(100)
      common/nam/k2
      common/fis/ifis
      common/levop/ldopt,igam
      common/hjk/jang,iq,rcss
      common/ss/sor,rr
      common/nr34/nr3,nr4,ke5,i3d,irtst,i3t2,ijkl,lq,tem(36)
      common/parfs/mz,k6,delrr(999)
      common/pl4/crs(24,13)
      common/shft4/k3,jcal
      common/shft5/fiss
      common/shft2/fs(25),dsp(25),brr(25),der(25),er(25)
     1,cx(25),bilsum,xmiss,sumiz
      REAL*8 fiss
      common/nhy/ij,jl,ji,jq,td,ex1,ex2,tmx,av,gav,cost,b(3)
      common/paro/pq,cross
      common/hyb2/pp(15,24,999)
      common/sft5/exc(18,27),xmax
      common/pl3/jx,nepr,eb(100),rcsp(100),zee,amass,plex
      common/lym1/be(18,27,8),pair(18,27),xmas(18,27),delshl(18,27),amas
     1(18,27)
c     common/bym1/symb(18,27),symbp(18,27),symbs(18,27)
      common/pl7/mas
      common/lab3/sig(8,999),slav(8,999)
      common/pl5/na,nz,title(20)
      common/scr/k9
      common/iso/qpn(3),qpnc
      common/incr/ed,tened,efinal
      common/tst/test
      common/send/irfr,iadst,dlt
      common/gam/io,la,ppgam(999),gmspec(999),beta,ams
      common/eqin/eind(100)
      common/ist/isot,iikl,knz,idelm,fracts,fract,ccrs(32,13,100)
     1,engbuf(100),exd,een(8,999,100)
      common/cshel/shel(18,27,2)
      common/new/fpyismr(3,100)
      common/gr/rint(1100)
      REAL*8 rint
      common/mbc/sigg(21,4,70)
      common/rc/rzero
      common/equ/equiv
      common/nuopt/kinem,iem,kb
      common/bound/neut2,iprot,idnuc(100),idm(100),einc(100),
     1thett(100),phit(100),ppp
      common/spins/rote(15,24,181),popej(15,24,180,100),fcx(180),ai,ecmp
      common/gate1/spind(5,180),popj(15,24,180),egatl(5),egatu(5),npart(
     15),igate
       common/spins3/ajmax
       common/spins4/eprot,ealpha,eoff(8)
      common/tagspin/trgspin
      common/tgspin/levno
      common/estar/estarcnt,e6,e7,e8
      common/sigsav/sigdr(100),sigqd(100)
      common/frags/xmasf(160,106),bq(106,160,4),pf(106,160)
     1,pfp(106,160,100),ppfi(106,160),ppfz(106),ppfa(260)
     2,ppfaz(106,260)
      common/frags2/rotef(260,181),arotef(260,181)
      REAL*8 powaz

      REAL*8 totcsf,totcsfg,totcsf1,totcsf2,totfpy,totfpyg,totfpy1      
      REAL*8 totfpy2, totf,toty                                       
c initialize totals                                                    
      idif=0
      totcsf=0.                                                         
      totcsfg=0.                                                     
      totcsf1=0.                                                      
      totcsf2=0.                                                       
      totfpy=0.                                                         
      totfpyg=0.                                                      
      totfpy1=0.                                                       
      totfpy2=0.                                                        
c
c  0000000000000000000000000000000000000000000000000000000000000000
      if(igate.eq.0)go to 88895
       if(limout.eq.0)then
      write(27,*)'  gated first nucleon out J distributions'
      write(27,*)'  gate energies,upper,lower:'
      write(27,*)' j ',egatl(1),egatl(2),egatl(3),egatl(4),egatl(5)
      write(27,*)' j ',egatu(1),egatu(2),egatu(3),egatu(4),egatu(5),'OM
     1 partial cs'
      write(27,*)'n=1,p=2',npart(1),npart(2),npart(3),npart(4),npart(5)
      write(27,*)
      endif
      do 88894 j=1,90
      sumj=spind(1,j)+spind(2,j)+spind(3,j)+spind(4,j)+spind(5,j)
      if(sumj.eq.0.)go to 88894
       k=j-1
       if(limout.eq.0)then
      write(27,88896) k,(spind(i,j),i=1,5),sigml(j)
       endif
88896 format(1x,i4,1x,6(1x,f8.2))
88894 continue
88895 continue
       if(limout.eq.0)then
      write(27,*)'E vs. J array'
      endif
      acn=at+ap
      zcn=zt+zp
c88796 continue
      do 88890 kz=1,13
      do 88889 ka=1,22
      ikz=zcn+1.-float(kz)
      ika=acn+2.-float(kz)-float(ka)
      if(ikz.le.0.or.ika.le.0)go to 88889
      isno=isoid(ikz,ika)
         nolev=levnm(isno)
c6/16/04 ------------
      if(isno.le.1)go to 88889
        if(nolev.gt.2)call iso2(isno,kz,ka)
        if(nolev.eq.2)call isomers(kz,ka,isno)
      iiz=zee+1.-float(kz)
      iia=amass+2.-float(kz)-float(ka)
       if(limout.eq.0)then
      write(27,*)'Z,A of isomer =',iiz,iia
c change isoj to isno here
      write(27,*)'gs spin = ',xgs(isno),' isomer spin =',xiso(isno),' is
     1omer energy = ',eiso(isno)
c6/16/04
      if(levno.gt.2)write(27,*)' isomer2 spin = ',x2iso(isno),' iso
     1mer2 energy = ',eiso2(isno)
      write(27,*)' ground state cs= ',csgse(isno,nepr),' isomer cs= ',cs
     1ise(isno,nepr)
c6/16/04 change if for m2
      if(levno.gt.2)write(27,*)' isomer2 cs = ',csis2e(isno,nepr)
        endif
       idif=nolevs(iiz,iia)
       
        qum=0.
        xum=0.
       do 88891 j=1,90
       do 88891 i=1,100
      ex=float(i)*ed-ed/2.
      if(ex.ge.rote(kz,ka,j))then
      xum=xum+popej(kz,ka,j,i)
      else
      qum=qum+popej(kz,ka,j,i)
      endif
 
 
       tot=qum+xum
88891 continue
      if(tot.eq.0.)go to 78870
       if(limout.eq.0)then
      write(27,*)
      write(27,*)'iz=',kz
      write(27,*)'ia=',ka
      write(27,*)
      if(idif.eq.0)write(27,1347)(rote(kz,ka,j),j=1,19,2)
      if(idif.ne.0)write(27,1347)(rote(kz,ka,j),j=2,20,2)
      write(27,*)
      if(idif.eq.0)write(27,1348)((k-1),k=1,10)
      if(idif.ne.0)write(27,1349) 0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.
     15
      write(27,*)
      endif
      do 88888 i=1,30
      e=float(i)*ed-ed/2.
       if(limout.eq.0)then
      if(idif.eq.0)write(27,1346)e,(popej(kz,ka,j,i),j=1,19,2)
      if(idif.ne.0)write(27,1346)e,(popej(kz,ka,j,i),j=2,20,2)
      endif
88888 continue
1346  format('E',f6.2,10(1x,f6.2))
1347  format('erot(j)',2x,f6.2,9(1x,f6.2))
1348  format(' J =',4x,10(2x,i3,2x))
1349  format(' J =',4x,10(2x,f4.1,1x))
       if(limout.eq.0)then
      write(27,13461)
13461 format(' SUM')
      if(idif.eq.0)write(27,13460)(popj(kz,ka,j),j=1,19,2)
      if(idif.ne.0)write(27,13460)(popj(kz,ka,j),j=2,20,2)
      write(27,*)'sum over E and J above yrast =',xum,' mb'
      write(27,*)'sum over E and J below yrast =',qum,' mb'
      write(27,*)
      write(27,*)'iz=',kz
      write(27,*)'ia=',ka
      write(27,*)
      if(idif.eq.0)write(27,1347)(rote(kz,ka,j),j=21,39,2)
      if(idif.ne.0)write(27,1347)(rote(kz,ka,j),j=22,40,2)
      write(27,*)
      if(idif.eq.0)write(27,1348)(k,k=10,19)
      if(idif.ne.0)write(27,1349) 10.5,11.5,12.5,13.5,14.5,15.5,16.5,17.
     15,18.5,19.5
      write(27,*)
      endif
      do 78888 i=1,30
      e=float(i)*ed-ed/2.
       if(limout.eq.0)then
      if(idif.eq.0)write(27,1346)e,(popej(kz,ka,j,i),j=21,39,2)
      if(idif.ne.0)write(27,1346)e,(popej(kz,ka,j,i),j=22,40,2)
       endif
78888 continue
       if(limout.eq.0)then
      if(idif.eq.0)write(27,13460)(popj(kz,ka,j),j=21,39,2)
      if(idif.ne.0)write(27,13460)(popj(kz,ka,j),j=22,40,2)
       endif
78870 continue
 
88889 continue
13460 format(' OVER E:',1f6.3,9(1x,f6.3))
88890 continue
   
 
      write(24,*)' Ebeam  Z    A    GS/tot    m1/tot    m2/tot   Jgs  J1
     1  J2   Egs       Em1       Em2'
98872 format(//,17x,'  Production Cross Section  ',                     
     1'           Isomer Spin, Excitation                Fission-Produ',
     2 'ct Yield',/,17x,28('-'),2x,43('-'),2x,28('-'),/,                
     3' Ebeam   Z    A    cs-g      cs-m1     cs-m2   J-gs J-m1',      
     1' J-m2   E-g       E-m1      E-m2     FPY-g     FPY-m1    FPY-m2',
     2'    TotFCS    TotFPY'/,' =====  ===  === ========  ======== ',   
     3' ========  ==== ==== ==== ========  ========  ========  =======',
     4'=  ========  ========  ========  ======== ')                     
98873 format(' =====           ========  ========  ======== '           
     1,46x,'========  ========  ========  ========  ========',/,        
     2  f6.1,' **Total**',3(1pe9.2,1x),45x,5(1pe9.2,1x),//)             
98876 format(f6.1,2i5,3(1pe9.2,1x),    3(0pf5.1),   8(1pe9.2,1x))       
98875 format(f6.1,2i5,2(1pe9.2,1x),10x,2(0pf5.1),5x,2(1pe9.2,1x),10x,   
     1 2(e9.2,1x),10x,e9.2,1x,e9.2)                                     
98874 format(f6.1,2i5, (1pe9.2,1x),35x,              (1pe9.2,1x),20x,   
     1 e9.2,21x,e9.2,1x,e9.2)                                           
      write(23,98872)                                                  
      toty=0.
      totfpy=0.
      do 98877 jz=10,70
      do 98888 jn=10,140
      toty=0.
       jzt=jz/2
       jzt2=2*jz
       if(jn.lt.jzt)go to 98888
       if(jn.gt.jzt2)go to 98888
      ja=jz+jn
      n=isoid(jz,ja)
      if(n.lt.1)n=1
      toty=0.                                                           
      if(n.le.1) then                                                   

        test=ppfaz(jz,ja)
        if(fizcs(nepr).le.0.)go to 98888
        if(test.le.0.) go to 98888                                      
        fpyismr(1,nepr)=ppfaz(jz,ja)/fizcs(nepr)                        
        toty=fpyismr(1,nepr)                                            
        totfpy=totfpy+fpyismr(1,nepr)                                   
        totfpyg=totfpyg+fpyismr(1,nepr)                                 
        totf=ppfaz(jz,ja)                                               
        totcsf=totcsf+ppfaz(jz,ja)                                      
        totcsfg=totcsfg+ppfaz(jz,ja)                                   
        write(23,98874) eind(nepr),jz,ja,ppfaz(jz,ja),                  
     1  englev(n,1),fpyismr(1,nepr),totf,toty                           
         go to 98888                                                    
        endif                                                           
      if(levnm(n).eq.3)then
c ground + 2  isomers produced                                         
      call iso2f(n,jz,ja)
     
      tesx=csgsef(n,nepr)+csisef(n,nepr)+csis2ef(n,nepr)
      if(tesx.gt.0.)then
       rgs=csgsef(n,nepr)/tesx
       rx1=csisef(n,nepr)/tesx
       rx2=csis2ef(n,nepr)/tesx
      
      if(fizcs(nepr).le.0.)write(*,*)' fizcs(nepr)=0.',nepr
      if(fizcs(nepr).le.0.)go to 5454
      fpyismr(1,nepr)= csgsef(n,nepr)/fizcs(nepr)                       
      fpyismr(2,nepr)= csisef(n,nepr)/fizcs(nepr)                       
      fpyismr(3,nepr)=csis2ef(n,nepr)/fizcs(nepr)                      
5454   continue
      toty=fpyismr(1,nepr)+fpyismr(2,nepr)+fpyismr(3,nepr)              
      totfpy=totfpy+toty                                                
      totfpyg=totfpyg+fpyismr(1,nepr)                                   
      totfpy1=totfpy1+fpyismr(2,nepr)                                   
      totfpy2=totfpy2+fpyismr(3,nepr)                                   
      totf=test                                                         
      totcsf=totcsf+totf                                                
      totcsfg=totcsfg+csgsef(n,nepr)                                    
      totcsf1=totcsf1+csisef(n,nepr)                                    
      totcsf2=totcsf2+csis2ef(n,nepr)                                   
c      write(23,*)levnm(n),n,'lev,n'
      write(23,98876)eind(nepr),jz,ja,csgsef(n,nepr),csisef(n,nepr),
     1csis2ef(n,nepr),bevj(n,1),bevj(n,2),bevj(n,3),englev(n,1),
     2englev(n,2),englev(n,3),fpyismr(1,nepr),fpyismr(2,nepr),          
     3fpyismr(3,nepr),totf,toty                                         
      write(24,98896)eind(nepr),jz,ja,rgs,rx1,
     1rx2,bevj(n,1),bevj(n,2),bevj(n,3),englev(n,1),
     2englev(n,2),englev(n,3)
      endif
      go to 98888
c new 9/15
      endif
      if(levnm(n).eq.2)then
c ground + 1 isomer produced                                            
      call isomersf(jz,ja,n )
      test=csgsef(n,nepr)+csisef(n,nepr)
      if(test.gt.0..and.fizcs(nepr).gt.0.)then
      fpyismr(1,nepr)= csgsef(n,nepr)/fizcs(nepr)                      
      fpyismr(2,nepr)= csisef(n,nepr)/fizcs(nepr)                      
      toty=fpyismr(1,nepr)+fpyismr(2,nepr)                              
      totfpy=totfpy+toty                                                
      totfpyg=totfpyg+fpyismr(1,nepr)                                   
      totfpy1=totfpy1+fpyismr(2,nepr)                                   
      rgs=csgsef(n,nepr)/(csgsef(n,nepr)+csisef(n,nepr))
      rx1=csisef(n,nepr)/(csgsef(n,nepr)+csisef(n,nepr))
      totf=test                                                         
      totcsf=totcsf+totf                                                
      totcsfg=totcsfg+csgsef(n,nepr)                                    
      totcsf1=totcsf1+csisef(n,nepr)                                    
      write(23,98875)eind(nepr),jz,ja,csgsef(n,nepr),csisef(n,nepr),   
     1bevj(n,1),bevj(n,2),englev(n,1),englev(n,2),fpyismr(1,nepr),      
c    2fpyismr(2,nepr),nepr,fizcs(nepr)                                 
     2fpyismr(2,nepr),totf,toty                                         
      write(24,98895)eind(nepr),jz,ja,rgs,rx1,
     1bevj(n,1),bevj(n,2),englev(n,1),englev(n,2)
       endif
       endif
98888 continue
98877 continue
      write(23,98873) eind(nepr),totcsfg,totcsf1,totcsf2,totfpyg,      
     1 totfpy1,totfpy2,totcsf,totfpy                                   
      write(25,98901)                                                  
      totfcx=totcsf/2.                                                 
98901 format(' Ebeam Nuclide    FP CX    FP Yield',/,' ===== ======= ',
     1'=========== ===========')                                        
98902 format(f6.1,i8,2(1pe12.4))                                        
98903 format(' Ebeam  A  Cross Sectn Fission Yld',/,                    
     1       ' ===== === =========== ===========')                      
98904 format(f6.1,i4,2(1pe12.4))                                        
98905 format(' ===== === =========== ===========',/,f6.1,' Tot',       
     1 2(1pe12.4),//)
      write(26,98903)                                                  
      do 98909 ja=20,200                                               
      ya=0.                                                           
      do 98908 jz=10,70                                               
      jn=ja-jz
      if(jn.lt.5)go to 98908
      n=isoid(jz,ja)                                                  
      cxval=ppfaz(jz,ja)                                              
      if(n.gt.1) cxval=csgsef(n,nepr)                                 
      if(totfcx.gt.0.)y=cxval/totfcx                                  
      if(totfcx.eq.0.)y= 0.                                           
      ya=ya+y                                                         
      toty=toty+y                                                      
      name=10000*ja + 10*jz                                            
      if(y.gt.0.) write(25,98902) eind(nepr),name,cxval,y              
      if(n.gt.1.and.totfcx.gt.0.) then                                 
        y=csisef(n,nepr)/totfcx                                        
        ya=ya+y                                                        
        toty=toty+y                                                    
        name=name+1                                                    
        if(y.gt.0.) write(25,98902) eind(nepr),name,csisef(n,nepr),y   
        endif                                                          
      if(levnm(n).eq.3.and.totfcx.gt.0.) then                          
        y=csis2ef(n,nepr)/totfcx                                       
        ya=ya+y                                                        
        toty=toty+y                                                    
        name=name+1                                                    
        if(y.gt.0.) write(25,98902) eind(nepr),name,csis2ef(n,nepr),y  
        endif                                                          
98908 continue                                                         
98909 if(ya.gt.0.) write(26,98904) eind(nepr),ja,ppfa(ja),ya           
      write(25,98910) eind(nepr),totcsf,toty                           
98910 format(' ===== ======= =========== ===========',/,f6.1,'**Total', 
     1 2(1pe12.4),//)                                               
      write(26,98905) eind(nepr),totcsf,toty                         

98896 format(f6.1,2i5,3(1pe9.2,1x),3(0p,f5.1),3(1pe9.2,1x))
98895 format(f6.1,2i5,2(1pe9.2,1x),10x,2(0p,f5.1),5x,2(1pe9.2,1x),10x)
      return
      end
      subroutine crslav(nen,icall,jlim)
      common/outlim/limout,lim2
      common/shft/mc,mp,inver,ike,ipch,kq,qval,ap,at,zp,zt,cld,barfac
      common/eqin/eind(100)
       common/scrn/att(10),cldd(10),bexp(24,15),reacss(100),edelt
      common/par3/eq,sigml(999),acrs(999)
c     common/avergl/ind,cbb,etem,i,ie,sigtemp
c     common/avergl/cbb,ei,sigtemp,ind,i,k
      common/avergl/cbb,etem,sigtemp,ind,i,ie
      common/lab3/sig(8,999),slav(8,999)
      if(icall.eq.0)then
      e=eind(nen)*at/(at+ap)
      coul=1.4393*zt*zp
      u=at*ap/(at+ap)
      ra=(1.18*at**(1./3.))*(1.-1./(1.18*at**2./3.))
      cb=coul/ra
      eef=eind(nen)*(at/(at+ap))-cb
      if(eef.le.cb)sigml(1)=reacss(nen)
      endif
       if(icall.eq.1)then
      cb=cbb
        app=1.
      if(i.eq.3)app=4.
      if(i.eq.4)app=2.
      u=app*(at+ap-app)/(at+ap)
      e=etem
      eef=e-cb
      atx=at+ap-app
      ra=(1.18*atx**(1./3.))*(1.-1./(1.18*atx**2./3.))
      endif
      jlim=1
      ai=1
      jj=1
      if(e.lt.cb)go to 110
      alambar=4.6/sqrt(u*eef)
      pilbarsq=31.4159*alambar**2
      ra=ra+alambar
      rfin=0.
      sigtot=0.
      jlim=1
      ai=0.
      suml=0.
      sum2l=0.
      do 100 k=1,90
      l=k
      rfin=rfin+alambar
      expval=-(ra-rfin)/0.55
       if(expval.gt.80.)then
       tlferm=0.
                        else
      tlferm=1./(1.+exp(expval))
                        endif
      if(tlferm.lt.0.0001)go to 99100
      suml=suml+tlferm*float(k)
      sum2l=sum2l+(2.*float(k)-1.)*tlferm*float(k)
      jlim=l-1
      ai=jlim
      if(icall.eq.0)sigml(l)=pilbarsq*(2.*float(l-1)+1.)*tlferm
      if(icall.eq.0)then
      sigtot=sigtot+sigml(l)
      endif
  100 continue
 
99100 continue
      if(icall.eq.0)go to 99102
      if(suml.gt.0.)then
      avel=sum2l/suml-1.
      if(avel.lt.0.)avel=0.
      else
      avel=0.
      endif
      slav(i,ie)=avel
      if(icall.eq.1)return
99102     if(sigtot.gt.0.)ratio=reacss(nen)/sigtot
      if(sigtot.eq.0.)ratio=1.
      jj=1
      do 105 i=1,100
      if(sigml(i).lt..001)go to 105
      sigml(i)=sigml(i)*ratio
      jj=i
  105 continue
  110 continue
      jlim=jj
      return
      end
 
      subroutine crslav1(nen,icall,jlim,j,aeff)
      common/outlim/limout,lim2
      common/shft/mc,mp,inver,ike,ipch,kq,qval,ap,at,zp,zt,cld,barfac
      common/eqin/eind(100)
       common/scrn/att(10),cldd(10),bexp(24,15),reacss(100),edelt
      common/par3/eq,sigml(999),acrs(999)
c     common/avergl/ind,cbb,etem,i,ie,sigtemp
       common/avergl/cbb,etem,sigtemp,ind,i,ie
c     common/avergl/cbb,ei,sigtemp,ind,i,k
      common/lab3/sig(8,999),slav(8,999)
      common/lab3f/slavf(100,4,999)
      icall=1
       iz=j
      cb=cbb
        app=1.
      if(i.eq.3)app=4.
      if(i.eq.4)app=2.
      u=app*(aeff-app)/(aeff)
      e=etem
      eef=e-cb
      atx=aeff-app
      ra=(1.18*atx**(1./3.))*(1.-1./(1.18*atx**2./3.))
c     endif
      jlim=1
      ai=1
      jj=1
      if(e.lt.cb)go to 110
      alambar=4.6/sqrt(u*eef)
      pilbarsq=31.4159*alambar**2
      ra=ra+alambar
      rfin=0.
      sigtot=0.
      jlim=1
      ai=0.
      suml=0.
      sum2l=0.
      do 100 k=1,90
      l=k
      rfin=rfin+alambar
      diff=(rfin-ra)/0.55
      if(diff.gt.80.)diff=80.
      tlferm=1./(1.+exp(diff))
      if(tlferm.lt.0.0001)go to 99100
      suml=suml+tlferm*float(k)
      sum2l=sum2l+(2.*float(k)-1.)*tlferm*float(k)
      jlim=l-1
      ai=jlim
      if(icall.eq.0)sigml(l)=pilbarsq*(2.*float(l-1)+1.)*tlferm
      if(icall.eq.0)then
      sigtot=sigtot+sigml(l)
      endif
  100 continue
 
99100 continue
      if(icall.eq.0)go to 99102
      if(suml.gt.0.)then
      avel=sum2l/suml-1.
      if(avel.lt.0.)avel=0.
      else
      avel=0.
      endif
      slavf(j,i,ie)=avel
      if(icall.eq.1)return
99102     if(sigtot.gt.0.)ratio=reacss(nen)/sigtot
  110 continue
      return
      end
 
      subroutine iso2f(n,jz,ja)
c integrates E-J plane to get ground state and isomer yields for 2 isomer case
      common/eder/ed2,v,vhalf,vsq
      common/outlim/limout,lim2
      common/ismrdat/isoid(106,260),levnm(1000),bevj(1000,3)
      common/ismrdat1/englev(1000,3),spingnd(106,260)
      common/ismrdat2/ididl(106,260),nolevs(106,260)
      common/isomer1/xgs(690),xiso(690),x2iso(690),isonum,isoint,iso2num
      common/isomer2/zs(690),as(690),eiso(690),eiso2(690)
      common/isomer3/iziso(690),iaiso(690),csgs(690),csis(690)
     1,csgse(690,100),csise(690,100),csis2e(690,100),rm1(690,100),
     2rm2(690,100)
      common/fissiso/csgsef(690,100),csisef(690,100),csis2ef(690,100)
      common/isomer4/odev(13,22)
      common/pl3/jx,nepr,eb(100),rcsp(100),zee,amass,plex
      common/spins/rote(15,24,181),popej(15,24,180,100),fcx(180),ai,ecmp
      common/sping/pfpej(600,180,100),pfpj(600,180)
      common/frags2/rotef(260,181),arotef(260,181)
      common/gate1/spind(5,180),popj(15,24,180),egatl(5),egatu(5),npart(
     15),igate
      common/incr/ed,tened,efinal
      n=isoid(jz,ja)
      isno=n
      ize=iziso(n)
      iae=iaiso(n)
      ize=jz
      iae=ja
      if(n.le.1)return
      idif=nolevs(jz,ja)
      levno=levnm(n)
13460 format(' OVER E:',1f6.3,9(1x,f6.3))
      
      ix1=englev(n,2)/ed+1.
      ix2=englev(n,3)/ed+1.
       xj0=bevj(n,1)
       xj1=bevj(n,2)
       xj2=bevj(n,3)
c define energy indices for isomer levels for integration limits
c define average spins between levels expressed as 2J+1
      j0j1=(2.*xj0+2.*xj1+2.)/2.
      j1j2=(2.*xj1+2.*xj2+2.)/2.
      j0j2=(2.*xj0+2.*xj2+2.)/2.
c now divide possible spin relationships
      csg=0.
      cs1=0.
      cs2=0.
      if(xj0.lt.xj1.and.xj1.lt.xj2)go to 100
      if(xj0.lt.xj2.and.xj2.lt.xj1)go to 200
      if(xj2.lt.xj1.and.xj1.lt.xj0)go to 300
      if(xj1.lt.xj2.and.xj2.lt.xj0)go to 400
      if(xj1.lt.xj0.and.xj0.lt.xj2)go to 500
      if(xj2.lt.xj0.and.xj0.lt.xj1)go to 600
  100 if( ix1.eq.1)go to 11
c integrate ground state below excitation m1
      do 10 ie=1,ix1-1
       do 9 j=1,90
        csg=csg+pfpej(n,j,ie)
    9 continue
   10 continue
   11 continue
c ground state below half spin diff to m1
       temp=0.
       do 8 ie=ix1,100
        do 7 j=1,j0j1-1
        csg=csg+pfpej(n,j,ie)
    7 continue
      temp=temp+pfpej(n,j0j1,ie)
    8 continue
      temp=temp/2.
      csg=csg+temp
c above adds half border cs tp csg;will add to cs1 later
c integrate for m1 below m2 energy
       do 6 ie=ix1,ix2-1
        do 5 j=j0j1+1,90
         cs1=cs1+pfpej(n,j,ie)
    5 continue
    6 continue
c integrate m1 contribution from cascade down
       cs1=cs1+temp
       temp=0.
       do 4 ie=ix2,100
        do 3 j=j0j1+1,j1j2-1
         cs1=cs1+pfpej(n,j,ie)
    3 continue
      temp=temp+pfpej(n,j1j2,ie)
    4 continue
      temp=temp/2.
      cs1=cs1+temp
c integrate m2
c integrate ground state below excitation m1
       do 13 ie=ix2,100
        do 12 j=j1j2+1,90
         cs2=cs2+pfpej(n,j,ie)
   12 continue
   13 continue
       cs2=cs2+temp
       go to 1000
  200  continue
       if(ix1.eq.1)go to 211
       do 210 ie=1,ix1-1
        do 209 j=1,90
         csg=csg+pfpej(n,j,ie)
  209 continue
  210 continue
  211 continue
c ground state below half spin diff to m2
      temp=0.
       do 208 ie=ix1,ix2-1
        do 207 j=1,j1j2-1
         csg=csg+pfpej(n,j,ie)
  207 continue
      temp=temp+pfpej(n,j1j2,ie)
  208 continue
      temp=temp/2.
      csg=csg+temp
      cs1=cs1+temp
      tem02=0.
      do 220 ie=ix2,100
       do 219 j=1,j0j2-1
         csg=csg+pfpej(n,j,ie)
  219 continue
      tem02=tem02+pfpej(n,j0j2,ie)
  220 continue
      tem02=tem02/2.
      csg=csg+tem02
      cs2=cs2+tem02
c integrate for m1 below m2 energy
       do 206 ie=ix1,ix2-1
        do 205 j=j0j1+1,90
         cs1=cs1+pfpej(n,j,ie)
  205 continue
  206 continue
c integrate m2 contribution from cascade down
      tem21=0.
       do 204 ie=ix2,100
        do 203 j=j0j2+1,j1j2-1
         cs2=cs2+pfpej(n,j,ie)
  203 continue
         tem21=tem21+pfpej(n,j1j2,ie)
  204 continue
      tem21=tem21/2.
      cs2=cs2+tem21
      cs1=cs1+tem21
 
c integrate m2
c integrate ground state below excitation m1
       do 213 ie=ix2,100
        do 212 j=j1j2+1,90
         cs1=cs1+pfpej(n,j,ie)
  212 continue
  213 continue
      go to 1000
  300  continue
       if( ix1.eq.1)go to 311
       do 310 ie=1,ix1-1
        do 309 j=1,90
         csg=csg+pfpej(n,j,ie)
  309 continue
  310 continue
  311 continue
c ground state below half spin diff to m2
       do 308 ie=ix1,100
        do 307 j=j0j1+1,90
         csg=csg+pfpej(n,j,ie)
  307 continue
  308 continue
      tem01=0.
      do 320 ie=ix1,ix2-1
       do 319 j=1,j0j1-1
         cs1=cs1+pfpej(n,j,ie)
  319 continue
         tem01=tem01+pfpej(n,j0j1,ie)
  320 continue
      tem01=tem01/2.
      cs1=cs1+tem01
      csg=csg+tem01
c integrate for m1
       tem01=0.
       do 306 ie=ix2,100
        do 305 j=j1j2+1,j0j1-1
         cs1=cs1+pfpej(n,j,ie)
  305 continue
         tem01=tem01+pfpej(n,j0j1,ie)
  306 continue
      tem01=tem01/2.
      cs1=cs1+tem01
      csg=csg+tem01
c integrate m2 contribution from cascade down
      tem21=0.
       do 304 ie=ix2,100
        do 303 j=1,j1j2-1
         cs2=cs2+pfpej(n,j,ie)
  303 continue
         tem21=tem21+pfpej(n,j1j2,ie)
  304 continue
      tem21=tem21/2.
      cs2=cs2+tem21
      cs1=cs1+tem21
      go to 1000
  400  continue
       if(ix1.eq.1)go to 411
       do 410 ie=1,ix1-1
        do 409 j=1,90
         csg=csg+pfpej(n,j,ie)
  409 continue
  410 continue
  411 continue
c ground state below half spin diff to m2
       tem02=0.
       do 408 ie=ix2,100
        do 407 j=j0j2+1,90
         csg=csg+pfpej(n,j,ie)
  407 continue
      tem02=tem02+pfpej(n,j0j2,ie)
  408 continue
      tem02=tem02/2.
      csg=csg+tem02
      cs2=cs2+tem02
      tem01=0.
      do 420 ie=ix1,ix2-1
       do 419 j=j0j1+1,90
         csg=csg+pfpej(n,j,ie)
  419 continue
      tem01=tem01+pfpej(n,j0j1,ie)
  420 continue
      tem01=tem01/2.
       csg=csg+tem01
       cs1=cs1+tem01
c integrate for m1 below m2 energy
       do 406 ie=ix1,ix2-1
        do 405 j=1,j0j1-1
         cs1=cs1+pfpej(n,j,ie)
  405 continue
  406 continue
c integrate m1 contribution from cascade down
      tem21=0.
       do 404 ie=ix2,100
        do 403 j=1,j1j2-1
         cs1=cs1+pfpej(n,j,ie)
  403 continue
      tem21=tem21+pfpej(n,j1j2,ie)
  404 continue
 
      tem21=tem21/2.
      cs1=cs1+tem21
      cs2=cs2+tem21
c integrate m2
c integrate ground state below excitation m1
       do 413 ie=ix2,100
        do 412 j=j1j2+1,j0j2-1
         cs2=cs2+pfpej(n,j,ie)
  412 continue
  413 continue
       go to 1000
  500  continue
       if(ix1.eq.1)go to 511
       do 510 ie=1,ix1-1
        do 509 j=1,90
         csg=csg+pfpej(n,j,ie)
  509 continue
  510 continue
  511 continue
c ground state below half spin diff to m2
      tem01=0.
       do 508 ie=ix1,ix2-1
        do 507 j=j0j1+1,90
         csg=csg+pfpej(n,j,ie)
  507 continue
      tem01=tem01+pfpej(n,j0j1,ie)
  508 continue
      tem01=tem01/2.
      csg=csg+tem01
      cs1=cs1+tem01
c temp stop here
      tem01=0.
      tem02=0.
      do 520 ie=ix2,100
       do 519 j=j0j1+1,j0j2-1
         csg=csg+pfpej(n,j,ie)
  519 continue
      tem01=tem01+pfpej(n,j0j1,ie)
      tem02=tem02+pfpej(n,j0j2,ie)
  520 continue
      tem01=tem01/2.
      tem02=tem02/2.
      csg=csg+tem01
      cs2=cs2+tem02
      csg=csg+tem02
      cs1=cs1+tem01
c integrate for m1 below m2 energy
       do 506 ie=ix1,100
        do 505 j=1,j0j1-1
         cs1=cs1+pfpej(n,j,ie)
  505 continue
  506 continue
c integrate m2 contribution from cascade down
       do 504 ie=ix2,100
        do 503 j=j0j2+1,90
         cs2=cs2+pfpej(n,j,ie)
  503 continue
  504 continue
 
       go to 1000
  600  continue
       if(ix1.eq.1)go to 611
       do 610 ie=1,ix1-1
        do 609 j=1,90
         csg=csg+pfpej(n,j,ie)
  609 continue
  610 continue
  611 continue
c ground state below half spin diff to m2
       tem01=0.
       do 608 ie=ix1,ix2-1
        do 607 j=1,j0j1-1
         csg=csg+pfpej(n,j,ie)
  607 continue
      tem01=tem01+pfpej(n,j0j1,ie)
  608 continue
      tem01=tem01/2.
      csg=csg+tem01
      cs1=cs1+tem01
       tem01=0.
       tem02=0.
      do 620 ie=ix2,100
       do 619 j=j0j2+1,j0j1-1
         csg=csg+pfpej(n,j,ie)
  619 continue
      tem01=tem01+pfpej(n,j0j1,ie)
      tem02=tem02+pfpej(n,j0j2,ie)
  620 continue
      tem01=tem01/2.
      tem02=tem02/2.
         csg=csg+tem01+tem02
       cs1=cs1+tem01
       cs2=cs2+tem02
       do 606 ie=ix1,100
        do 605 j=j0j1+1,90
         cs1=cs1+pfpej(n,j,ie)
  605 continue
  606 continue
c integrate m2 contribution from cascade down
       do 604 ie=ix2,100
        do 603 j=1,j0j2-1
         cs2=cs2+pfpej(n,j,ie)
  603 continue
  604  continue
       continue
 1000  continue
       csgsef(n,nepr)=csg
       csisef(n,nepr)=cs1
       csis2ef(n,nepr)=cs2
      write(37,*)'ground state cs=',csg,'isomer1=',cs1,
     1'isomer2=',cs2
      write(37,*)'E(meV) J=   1-10'
c801   format('E=',1f6.3,10(1x,f6.3))
      write(37,*)'Z,A of isomer =',jz,ja
c change isoj to isno here
      write(37,*)'gs spin = ',xgs(isno),' isomer spin =',xiso(isno),' is
     1omer energy = ',eiso(isno)
      if(levno.gt.2)write(37,*)' isomer2 spin = ',x2iso(isno),' iso
     1mer2 energy = ',eiso2(isno)
      write(37,*)' ground state cs= ',csgse(isno,nepr),' isomer cs= ',cs
     1ise(isno,nepr)
      if(levno.gt.2)write(37,*)' isomer2 cs = ',csis2e(isno,nepr)
       
        qum=0.
        xum=0.
       do 88891 j=1,99
       do 88891 i=1,100
      ex=float(i)*ed-ed/2.
      if(ex.ge.rotef(iae,j))then
      xum=xum+pfpej(n,j,i)
      else
      qum=qum+pfpej(n,j,i)
      endif
 
 
       tot=qum+xum
88891 continue
      if(tot.eq.0.)go to 78870
c *******
c *********
      if(idif.eq.0)write(37,1347)(rotef(iae,j),j=1,19,2)
      if(idif.ne.0)write(37,1347)(rotef(iae,j),j=2,20,2)
      write(37,*)
      if(idif.eq.0)write(37,1348)((k-1),k=1,10)
      if(idif.ne.0)write(37,1349) 0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.
     15
      write(37,*)
      do 800 i=1,50
      sm1=0.
      do 906 ij=1,20
      sm1=sm1+pfpej(n,ij,i)
906   continue
      if(sm1.gt.0.)then
      e=float(i)*ed-ed/2.
      if(idif.eq.0)write(37,1346)e,(pfpej(n,j,i),j=1,19,2)
      if(idif.ne.0)write(37,1346)e,(pfpej(n,j,i),j=2,20,2)
c88888 continue
1346  format('E',f6.2,10(1x,f6.3))
1347  format('erot(j)',2x,f6.3,9(1x,f6.3))
1348  format(' J =',4x,10(2x,i3,2x))
1349  format(' J =',4x,10(2x,f4.1,1x))
       endif
 800   continue
      write(37,13461)
13461 format(' SUM')
      if(idif.eq.0)write(27,13460)(pfpj(n,j),j=1,19,2)
      if(idif.ne.0)write(27,13460)(pfpj(n,j),j=2,20,2)
      write(37,*)'sum over E and J above yrast =',xum,' mb'
      write(37,*)'sum over E and J below yrast =',qum,' mb'
c ****************end new algo
      if(idif.eq.0)write(37,1347)(rotef(iae,j),j=21,39,2)
      if(idif.ne.0)write(37,1347)(rotef(iae,j),j=22,40,2)
      write(37,*)
      if(idif.eq.0)write(37,1348)(k,k=10,19)
      if(idif.ne.0)write(37,1349) 10.5,11.5,12.5,13.5,14.5,15.5,16.5,17.
     15,18.5,19.5
      write(37,*)
      do 805 i=1,50
       sm2=0.
      do 907 ij=21,40
      sm2=sm2+pfpej(n,ij,i)
907   continue
       if(sm2.gt.0.)then
      e=float(i)*ed-ed/2.
      if(idif.eq.0)write(37,1346)e,(pfpej(n,j,i),j=21,39,2)
      if(idif.ne.0)write(37,1346)e,(pfpej(n,j,i),j=22,40,2)
      endif
 805   continue
      if(idif.eq.0)write(37,13460)(pfpj(n,j),j=21,39,2)
      if(idif.ne.0)write(37,13460)(pfpj(n,j),j=22,40,2)
78870 continue
       return
       end
      subroutine isomersf(iqz,iqa,n)
c integrates E-J plane to get ground state and isomer yields
      common/frags2/rotef(260,181),arotef(260,181)
      common/outlim/limout,lim2
      common/fissiso/csgsef(690,100),csisef(690,100),csis2ef(690,100)
      common/ismrdat/isoid(106,260),levnm(1000),bevj(1000,3)
      common/ismrdat1/englev(1000,3),spingnd(106,260)
      common/ismrdat2/ididl(106,260),nolevs(106,260)
      common/isomer1/xgs(690),xiso(690),x2iso(690),isonum,isoint,iso2num
      common/isomer2/zs(690),as(690),eiso(690),eiso2(690)
      common/isomer3/iziso(690),iaiso(690),csgs(690),csis(690)
     1,csgse(690,100),csise(690,100),csis2e(690,100),rm1(690,100),
     2rm2(690,100)
      common/isomer4/odev(13,22)
      common/pl3/jx,nepr,eb(100),rcsp(100),zee,amass,plex
      common/spins/rote(15,24,181),popej(15,24,180,100),fcx(180),ai,ecmp
      common/sping/pfpej(600,180,100),pfpj(600,180)
      common/gate1/spind(5,180),popj(15,24,180),egatl(5),egatu(5),npart(
     15),igate
      common/incr/ed,tened,efinal
      idif=nolevs(iqz,iqa)
      isoj=n
      isno=isoj
      if(isoj.le.1)return
      gsspin=bevj(n,1)
      spiniso=bevj(n,2)
      eisomr=englev(n,2)
c now change criteria to branch at average of gs and isomer spins
          spinav=(gsspin+spiniso)/2.
          ispinav=(2.*spinav+1.1)
         sumgs=0.
         sumiso=0.
         indxis=(eisomr/ed+1.)
           ispiniso=(2.*spiniso+1.1)
      if(spiniso.lt.gsspin)go to 400
         if(indxis.lt.2)go to 101
           do 100 ie=1,indxis-1
           do 99 ij=1,120
         sumgs=sumgs+pfpej(n,ij,ie)
  99       continue
 100       continue
      if(spiniso.lt.gsspin)go to 400
c to this point,sumgs gives integral over energies below isomer
 101       continue
c integrate cross section below isomer excitation into g.s. xs
c next sum cs of spin lt  average spin into gs
          if(indxis.lt.1)indxis=1
                 sumav=0.
          do 200 ie=indxis,100
          do 199 ij=1,ispinav-1
         sumgs=sumgs+pfpej(n,ij,ie)
 199       continue
         sumav=sumav+pfpej(n,ispinav,ie)
 200       continue
             sumav=sumav/2.
             sumgs=sumgs+sumav
c now have integrated ground state;next integrate high spin isomer
           do 300 ie=indxis,100
           do 299 ij=ispinav+1,120
         sumiso=sumiso+pfpej(n,ij,ie)
 299       continue
 300       continue
           sumiso=sumiso+sumav
      csgsef(isoj,nepr)=sumgs
      csisef(isoj,nepr)=sumiso
      go to 800

 400  continue
c do sums for case isomer spin less than ground state spin
      temp=0.
      do 500 iex=1,100
      do 499 jin=1,120
      temp=temp+pfpej(n,jin,iex)
 499  continue
 500  continue
      igspin=(2.*gsspin+1.1)
      sumav=0.
      sumiso=0.
      do 600 ie=indxis,100
      do 599 ij=1,ispinav-1
      sumiso=sumiso+pfpej(n,ij,ie)
 599  continue
           sumav=sumav+pfpej(n,ispinav,ie)
 600  continue
       sumiso=sumiso+sumav/2.
      sumgs=temp-sumiso
      csgsef(n,nepr)=sumgs
      csisef(n,nepr)=sumiso
      if(sumgs.eq.0..and.sumiso.eq.0.)go to 911
c *******
c *********
1346  format('E',f6.3,10(1x,f6.3))
1347  format('erot(j)',2x,f6.3,9(1x,f6.3))
1348  format(' J =',4x,10(2x,i3,2x))
1349  format(' J =',4x,10(2x,f4.1,1x))
 800   continue
           write(37,*)'isomer yields for A=',iqa,' Z= ',iqz
      write(37,*)' Ground state yield = ',sumgs ,' mb'
      write(37,*)' Isomer yield = ',sumiso,' mb'
c     endif
      csgsef(isoj,nepr)=sumgs
      csisef(isoj,nepr)=sumiso
c801   format('E=',1f6.3,10(1x,f6.3))
c 802   continue
      write(37,*)'Z,A of isomer =',iqz,iqa
c change isoj to isno here
      write(37,*)'gs spin = ',xgs(isno),' isomer spin =',xiso(isno),' is
     1omer energy = ',eiso(isno)
c6/16/04
c ****************end new algo
      if(idif.eq.0)write(37,1347)(rotef(iqa,j),j=1,19,2)
      if(idif.ne.0)write(37,1347)(rotef(iqa,j),j=2,20,2)
      write(37,*)
      if(idif.eq.0)write(37,1348)((k-1),k=1,10)
      if(idif.ne.0)write(37,1349) 0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.
     15
      write(37,*)
      do 807 i=1,50
      sm1=0.
      do 906 ij=1,20
      sm1=sm1+pfpej(n,ij,i)
906   continue
      if(sm1.gt.0.)then
      e=float(i)*ed-ed/2.
      if(idif.eq.0)write(37,1346)e,(pfpej(n,j,i),j=1,19,2)
      if(idif.ne.0)write(37,1346)e,(pfpej(n,j,i),j=2,20,2)
c88888 continue
       endif
 807   continue
      write(37,13461)
13461 format(' SUM')
      if(idif.eq.0)write(27,13460)(pfpj(n,j),j=1,19,2)
      if(idif.ne.0)write(27,13460)(pfpj(n,j),j=2,20,2)
c     write(37,*)'sum over E and J above yrast =',xum,' mb'
c     write(37,*)'sum over E and J below yrast =',qum,' mb'
c ****************end new algo
      if(idif.eq.0)write(37,1347)(rotef(iqa,j),j=21,39,2)
      if(idif.ne.0)write(37,1347)(rotef(iqa,j),j=22,40,2)
      write(37,*)
      if(idif.eq.0)write(37,1348)(k,k=10,19)
      if(idif.ne.0)write(37,1349) 10.5,11.5,12.5,13.5,14.5,15.5,16.5,17.
     15,18.5,19.5
      write(37,*)
      do 805 i=1,50
       sm2=0.
      do 907 ij=21,40
      sm2=sm2+pfpej(n,ij,i)
907   continue
       if(sm2.gt.0.)then
      e=float(i)*ed-ed/2.
      if(idif.eq.0)write(37,1346)e,(pfpej(n,j,i),j=21,39,2)
      if(idif.ne.0)write(37,1346)e,(pfpej(n,j,i),j=22,40,2)
c78888 continue
      endif
 805   continue
13460 format(' OVER E:',1f6.3,9(1x,f6.3))
      if(idif.eq.0)write(37,13460)(pfpj(n,j),j=21,39,2)
      if(idif.ne.0)write(37,13460)(pfpj(n,j),j=22,40,2)
c78870 continue
911   continue
      return
      end
      subroutine iso2(isno,jz,ja)
      common/outlim/limout,lim2
c integrates E-J plane to get ground state and isomer yields for 2 isomer case
      common/isomer1/xgs(690),xiso(690),x2iso(690),isonum,isoint,iso2num
      common/isomer2/zs(690),as(690),eiso(690),eiso2(690)
      common/isomer3/iziso(690),iaiso(690),csgs(690),csis(690)
     1,csgse(690,100),csise(690,100),csis2e(690,100),rm1(690,100),
     2rm2(690,100)
      common/isomer4/odev(13,22)
      common/pl3/jx,nepr,eb(100),rcsp(100),zee,amass,plex
      common/spins/rote(15,24,181),popej(15,24,180,100),fcx(180),ai,ecmp
      common/sping/pfpej(600,180,100),pfpj(600,180)
      common/gate1/spind(5,180),popj(15,24,180),egatl(5),egatu(5),npart(
     15),igate
      common/incr/ed,tened,efinal
      xj0=xgs(isno)
      xj1=xiso(isno)
      xj2=x2iso(isno)
c define energy indices for isomer levels for integration limits
      ix1=(eiso(isno)/ed +1.)
      ix2=(eiso2(isno)/ed+1.)
c define average spins between levels expressed as 2J+1
      j0j1=(2.*xj0+2.*xj1+2.)/2.
      j1j2=(2.*xj1+2.*xj2+2.)/2.
      j0j2=(2.*xj0+2.*xj2+2.)/2.
c now divide possible spin relationships
      csg=0.
      cs1=0.
      cs2=0.
      if(xj0.lt.xj1.and.xj1.lt.xj2)go to 100
      if(xj0.lt.xj2.and.xj2.lt.xj1)go to 200
      if(xj2.lt.xj1.and.xj1.lt.xj0)go to 300
      if(xj1.lt.xj2.and.xj2.lt.xj0)go to 400
      if(xj1.lt.xj0.and.xj0.lt.xj2)go to 500
      if(xj2.lt.xj0.and.xj0.lt.xj1)go to 600
  100 if( ix1.eq.1)go to 11
c integrate ground state below excitation m1
      do 10 ie=1,ix1-1
       do 9 j=1,90
        csg=csg+popej(jz,ja,j,ie)
    9 continue
   10 continue
   11 continue
c ground state below half spin diff to m1
       temp=0.
       do 8 ie=ix1,100
        do 7 j=1,j0j1-1
         csg=csg+popej(jz,ja,j,ie)
    7 continue
      temp=temp+popej(jz,ja,j0j1,ie)
    8 continue
      temp=temp/2.
      csg=csg+temp
c above adds half border cs tp csg;will add to cs1 later
c integrate for m1 below m2 energy
       do 6 ie=ix1,ix2-1
        do 5 j=j0j1+1,90
         cs1=cs1+popej(jz,ja,j,ie)
    5 continue
    6 continue
c integrate m1 contribution from cascade down
       cs1=cs1+temp
       temp=0.
       do 4 ie=ix2,100
        do 3 j=j0j1+1,j1j2-1
         cs1=cs1+popej(jz,ja,j,ie)
    3 continue
      temp=temp+popej(jz,ja,j1j2,ie)
    4 continue
      temp=temp/2.
      cs1=cs1+temp
c integrate m2
c integrate ground state below excitation m1
       do 13 ie=ix2,100
        do 12 j=j1j2+1,90
         cs2=cs2+popej(jz,ja,j,ie)
   12 continue
   13 continue
       cs2=cs2+temp
       go to 1000
  200  continue
       if(ix1.eq.1)go to 211
       do 210 ie=1,ix1-1
        do 209 j=1,90
        csg=csg+popej(jz,ja,j,ie)
  209 continue
  210 continue
  211 continue
c ground state below half spin diff to m2
      temp=0.
       do 208 ie=ix1,ix2-1
        do 207 j=1,j1j2-1
         csg=csg+popej(jz,ja,j,ie)
  207 continue
      temp=temp+popej(jz,ja,j1j2,ie)
  208 continue
      temp=temp/2.
      csg=csg+temp
      cs1=cs1+temp
      tem02=0.
      do 220 ie=ix2,100
       do 219 j=1,j0j2-1
         csg=csg+popej(jz,ja,j,ie)
  219 continue
      tem02=tem02+popej(jz,ja,j0j2,ie)
  220 continue
      tem02=tem02/2.
      csg=csg+tem02
      cs2=cs2+tem02
c integrate for m1 below m2 energy
       do 206 ie=ix1,ix2-1
        do 205 j=j0j1+1,90
         cs1=cs1+popej(jz,ja,j,ie)
  205 continue
  206 continue
c integrate m2 contribution from cascade down
      tem21=0.
       do 204 ie=ix2,100
        do 203 j=j0j2+1,j1j2-1
         cs2=cs2+popej(jz,ja,j,ie)
  203 continue
       tem21=tem21+popej(jz,ja,j1j2,ie)
  204 continue
      tem21=tem21/2.
      cs2=cs2+tem21
      cs1=cs1+tem21
 
c integrate m2
c integrate ground state below excitation m1
       do 213 ie=ix2,100
        do 212 j=j1j2+1,90
         cs1=cs1+popej(jz,ja,j,ie)
  212 continue
  213 continue
      go to 1000
  300  continue
       if( ix1.eq.1)go to 311
       do 310 ie=1,ix1-1
        do 309 j=1,90
        csg=csg+popej(jz,ja,j,ie)
  309 continue
  310 continue
  311 continue
c ground state below half spin diff to m2
       do 308 ie=ix1,100
        do 307 j=j0j1+1,90
         csg=csg+popej(jz,ja,j,ie)
  307 continue
  308 continue
      tem01=0.
      do 320 ie=ix1,ix2-1
       do 319 j=1,j0j1-1
         cs1=cs1+popej(jz,ja,j,ie)
  319 continue
      tem01=tem01+popej(jz,ja,j0j1,ie)
  320 continue
      tem01=tem01/2.
      cs1=cs1+tem01
      csg=csg+tem01
c integrate for m1
       tem01=0.
       do 306 ie=ix2,100
        do 305 j=j1j2+1,j0j1-1
         cs1=cs1+popej(jz,ja,j,ie)
  305 continue
      tem01=tem01+popej(jz,ja,j0j1,ie)
  306 continue
      tem01=tem01/2.
      cs1=cs1+tem01
      csg=csg+tem01
c integrate m2 contribution from cascade down
      tem21=0.
       do 304 ie=ix2,100
        do 303 j=1,j1j2-1
         cs2=cs2+popej(jz,ja,j,ie)
  303 continue
      tem21=tem21+popej(jz,ja,j1j2,ie)
  304 continue
      tem21=tem21/2.
      cs2=cs2+tem21
      cs1=cs1+tem21
      go to 1000
  400  continue
       if(ix1.eq.1)go to 411
       do 410 ie=1,ix1-1
        do 409 j=1,90
        csg=csg+popej(jz,ja,j,ie)
  409 continue
  410 continue
  411 continue
c ground state below half spin diff to m2
       tem02=0.
       do 408 ie=ix2,100
        do 407 j=j0j2+1,90
         csg=csg+popej(jz,ja,j,ie)
  407 continue
      tem02=tem02+popej(jz,ja,j0j2,ie)
  408 continue
      tem02=tem02/2.
      csg=csg+tem02
      cs2=cs2+tem02
      tem01=0.
      do 420 ie=ix1,ix2-1
       do 419 j=j0j1+1,90
         csg=csg+popej(jz,ja,j,ie)
  419 continue
      tem01=tem01+popej(jz,ja,j0j1,ie)
  420 continue
      tem01=tem01/2.
       csg=csg+tem01
       cs1=cs1+tem01
c integrate for m1 below m2 energy
       do 406 ie=ix1,ix2-1
        do 405 j=1,j0j1-1
         cs1=cs1+popej(jz,ja,j,ie)
  405 continue
  406 continue
c integrate m1 contribution from cascade down
      tem21=0.
       do 404 ie=ix2,100
        do 403 j=1,j1j2-1
         cs1=cs1+popej(jz,ja,j,ie)
  403 continue
      tem21=tem21+popej(jz,ja,j1j2,ie)
  404 continue
 
      tem21=tem21/2.
      cs1=cs1+tem21
      cs2=cs2+tem21
c integrate m2
c integrate ground state below excitation m1
       do 413 ie=ix2,100
        do 412 j=j1j2+1,j0j2-1
         cs2=cs2+popej(jz,ja,j,ie)
  412 continue
  413 continue
       go to 1000
  500  continue
       if(ix1.eq.1)go to 511
       do 510 ie=1,ix1-1
        do 509 j=1,90
        csg=csg+popej(jz,ja,j,ie)
  509 continue
  510 continue
  511 continue
c ground state below half spin diff to m2
      tem01=0.
       do 508 ie=ix1,ix2-1
        do 507 j=j0j1+1,90
         csg=csg+popej(jz,ja,j,ie)
  507 continue
      tem01=tem01+popej(jz,ja,j0j1,ie)
  508 continue
      tem01=tem01/2.
      csg=csg+tem01
      cs1=cs1+tem01
c temp stop here
      tem01=0.
      tem02=0.
      do 520 ie=ix2,100
       do 519 j=j0j1+1,j0j2-1
         csg=csg+popej(jz,ja,j,ie)
  519 continue
      tem01=tem01+popej(jz,ja,j0j1,ie)
      tem02=tem02+popej(jz,ja,j0j2,ie)
  520 continue
      tem01=tem01/2.
      tem02=tem02/2.
      csg=csg+tem01+tem02
      cs2=cs2+tem02
      cs1=cs1+tem01
c integrate for m1 below m2 energy
       do 506 ie=ix1,100
        do 505 j=1,j0j1-1
         cs1=cs1+popej(jz,ja,j,ie)
  505 continue
  506 continue
c integrate m2 contribution from cascade down
       do 504 ie=ix2,100
        do 503 j=j0j2+1,90
         cs2=cs2+popej(jz,ja,j,ie)
  503 continue
  504 continue
 
       go to 1000
  600  continue
       if(ix1.eq.1)go to 611
       do 610 ie=1,ix1-1
        do 609 j=1,90
        csg=csg+popej(jz,ja,j,ie)
  609 continue
  610 continue
  611 continue
c ground state below half spin diff to m2
       tem01=0.
       do 608 ie=ix1,ix2-1
        do 607 j=1,j0j1-1
         csg=csg+popej(jz,ja,j,ie)
  607 continue
      tem01=tem01+popej(jz,ja,j0j1,ie)
  608 continue
      tem01=tem01/2.
      csg=csg+tem01
      cs1=cs1+tem01
       tem01=0.
       tem02=0.
      do 620 ie=ix2,100
       do 619 j=j0j2+1,j0j1-1
         csg=csg+popej(jz,ja,j,ie)
  619 continue
      tem01=tem01+popej(jz,ja,j0j1,ie)
      tem02=tem02+popej(jz,ja,j0j2,ie)
  620 continue
      tem01=tem01/2.
      tem02=tem02/2.
         csg=csg+tem01+tem02
       cs1=cs1+tem01
       cs2=cs2+tem02
       do 606 ie=ix1,100
        do 605 j=j0j1+1,90
         cs1=cs1+popej(jz,ja,j,ie)
  605 continue
  606 continue
c integrate m2 contribution from cascade down
       do 604 ie=ix2,100
        do 603 j=1,j0j2-1
         cs2=cs2+popej(jz,ja,j,ie)
  603 continue
  604 continue
 1000  continue
       csgse(isno,nepr)=csg
       csise(isno,nepr)=cs1
       csis2e(isno,nepr)=cs2
       return
       end
      subroutine isomers(iqz,iqa,isoj)
c integrates E-J plane to get ground state and isomer yields
      common/ismrdat2/ididl(106,260),nolevs(106,260)
      common/outlim/limout,lim2
      common/isomer1/xgs(690),xiso(690),x2iso(690),isonum,isoint,iso2num
      common/isomer2/zs(690),as(690),eiso(690),eiso2(690)
      common/isomer3/iziso(690),iaiso(690),csgs(690),csis(690)
     1,csgse(690,100),csise(690,100),csis2e(690,100),rm1(690,100),
     2rm2(690,100)
      common/isomer4/odev(13,22)
      common/pl3/jx,nepr,eb(100),rcsp(100),zee,amass,plex
      common/spins/rote(15,24,181),popej(15,24,180,100),fcx(180),ai,ecmp
      common/sping/pfpej(600,180,100),pfpj(600,180)
      common/gate1/spind(5,180),popj(15,24,180),egatl(5),egatu(5),npart(
     15),igate
      common/incr/ed,tened,efinal
      ikz=zee+1.-float(iqz)
      ika=amass+2.-float(iqz)-float(iqa)
      idif=nolevs(ikz,ika)
      gsspin=xgs(isoj)
      spiniso=xiso(isoj)
      eisomr=eiso(isoj)
c now change criteria to branch at average of gs and isomer spins
          spinav=(gsspin+spiniso)/2.
          ispinav=(2.*spinav+1.)
         sumgs=0.
         sumiso=0.
         indxis=(eisomr/ed+1.)
           ispiniso=(2.*spiniso+1.)
      if(spiniso.lt.gsspin)go to 400
         if(indxis.lt.2)go to 101
           do 100 ie=1,indxis-1
           do 99 ij=1,90
               sumgs=sumgs+popej(iqz,iqa,ij,ie)
  99       continue
 100       continue
 101       continue
c integrate cross section below isomer excitation into g.s. xs
c next sum cs of spin lt  average spin into gs
          if(indxis.lt.1)indxis=1
                 sumav=0.
          do 200 ie=indxis,100
          do 199 ij=1,ispinav-1
              sumgs=sumgs+popej(iqz,iqa,ij,ie)
 199       continue
           sumav=sumav+popej(iqz,iqa,ispinav,ie)
 200       continue
             sumav=sumav/2.
             sumgs=sumgs+sumav
c now have integrated ground state;next integrate high spin isomer
           do 300 ie=indxis,100
           do 299 ij=ispinav+1,90
           sumiso=sumiso+popej(iqz,iqa,ij,ie)
 299       continue
 300       continue
           sumiso=sumiso+sumav
           write(27,*)'isomer yields for A=',iqa,' Z= ',iqz
      write(27,*)' Ground state yield = ',sumgs ,' mb'
      write(27,*)' Isomer yield = ',sumiso,' mb'
      csgse(isoj,nepr)=sumgs
      csise(isoj,nepr)=sumiso
      return
 400  continue
c do sums for case isomer spin less than ground state spin
      temp=0.
      do 500 iex=1,100
      do 499 jin=1,90
      temp=temp+popej(iqz,iqa,jin,iex)
 499  continue
 500  continue
      igspin=(2.*gsspin+1.)
      sumav=0.
      do 600 ie=indxis,100
      do 599 ij=1,ispinav-1
      sumiso=sumiso+popej(iqz,iqa,ij,ie)
 599  continue
           sumav=sumav+popej(iqz,iqa,ispinav,ie)
 600  continue
       sumiso=sumiso+sumav/2.
      sumgs=temp-sumiso
       if(limout.eq.0)then
           write(27,*)'isomer yields for A=',iqa,' Z= ',iqz
      write(27,*)' Ground state yield = ',sumgs ,' mb'
      write(27,*)' Isomer yield = ',sumiso,' mb'
      endif
      csgse(isoj,nepr)=sumgs
      csise(isoj,nepr)=sumiso
      return
      end
c----------------------------------------------------------------------------------
      subroutine ldmagic
c     routine to select best level density model based on proximity to closed shells
c     if user does not select, i.e. default routine
      common/shft/mc,mp,inver,ike,ipch,kq,qval,ap,at,zp,zt,cld,barfac
      common/levop/ldopt,igam
      ldopt=0
      zee=zp+zt
      n=ap+at-zee
      kz=zee
      if((kz-28).le.2.and.(kz-28).gt.-2)ldopt=1
      if((n-28).le.2.and.(n-28).gt.-2)ldopt=1
      if((kz-20).le.2.and.(kz-20).gt.-2)ldopt=1
      if((n-20).le.2.and.(n-20).gt.-2)ldopt=1
      if((kz-8).le.2.and.(kz-8).gt.-2)ldopt=1
      if((n-8).le.2.and.(n-8).gt.-2)ldopt=1
      
      if((kz-50).le.2.and.(kz-50).gt.-2)ldopt=1
      if((kz-82).le.2.and.(kz-82).gt.-2)ldopt=1
      if((n -50).le.2.and.(n-50).gt.-2)ldopt=1
      if((n -82).le.3.and.(n-82).gt.-2)ldopt=1
      if((n -126).le.3.and.(n-126).gt.-2)ldopt=1
      write(*,*)
      write(*,*)' ldopt set internally =',ldopt
      write(7,*)' ldopt set internally =',ldopt
      if(ldopt.eq.0)write(*,*)'Fermi gas level densities'
      if(ldopt.eq.0)write(7,*)'Fermi gas level densities'
      if(ldopt.eq.1)write(*,*)'Kataria-Ramamurthy shell dependent ld'
      if(ldopt.eq.1)write(7,*)'Kataria-Ramamurthy shell dependent ld'
      write(7,*)
      write(*,*)
      return
      end
c 00000000000000000000000000000000000000000000000000000000000
      subroutine findismr
c sr to find and hold isomer data from file billiso
      common/outlim/limout,lim2
      common/ismrdat/isoid(106,260),levnm(1000),bevj(1000,3)
      common/ismrdat1/englev(1000,3),spingnd(106,260)
      common/ismrdat2/ididl(106,260),nolevs(106,260)
      common/pl3/jx,nepr,eb(100),rcsp(100),zee,amass,plex
      common/pl5/na,nz,title(20)
      character*1 plmn
      character*1 plmn2
      character*1 plmn3
      n=1
      levnm(1)=1
c     nlines=2640
      nlines=5275
      do 100 i=1,nlines
      read(28,1000)jnum,nlev,e1,s1,plmn,e2,s2,plmn2,e3,s3,plmn3,id
      jz=jnum/1000
      ja=jnum-1000*jz
      spingnd(jz,ja)=s1
      ididl(jz,ja)=id
      if(nlev.eq.1)isoid(jz,ja)=1
      if(jnum.eq.105260)return
      if(nlev.le.1)go to 100
      n=n+1
      isoid(jz,ja)=n
      levnm(n)=nlev
      bevj(n,1)=s1
      bevj(n,2)=s2
      bevj(n,3)=s3
      englev(n,1)=0.
      englev(n,2)=e2
      englev(n,3)=e3
  100 continue
      write(*,*)'number of isomers found =',n
 1000 format(1x,1i6,2x,1i1,1e12.5,1x,1f4.1,1a1,1e12.5,1x,1f4.1,1a1,
     11e12.5,1x,1f4.1,1a1,1i1)
      rewind(28)
      return
      end
 
      subroutine fizstor(emax,fis1,spinr)
c sub to store fission contributions according to z,a,E)
      REAL*8 fis1
      common/pl3/jx,nepr,eb(100),rcsp(100),zee,amass,plex
      common/outlim/limout,lim2
      common/fizzc/popfiz(13,20,150),pofizbig(13,20,150,30)
      common/scat/jz,ja,probxx(2,2),q1(8,1100),fcs,eloss,fcspec
c EMAX is excitation of fissioning nucleus,JZ is a Z index,JA an A index
c spinr is spin of fissioning nucleus in integer unitswhen used
      jf=abs(spinr+0.5)
      if(jf.lt.1)jf=1
      if(jf.gt.30)jf=30
      iu=(emax/2.+1.)
      if(iu.gt.150)write(*,*)'fizzer excitation exceeds 300 MeV'
      if(iu.gt.150)iu=150
      if(iu.lt.1)iu=1
      if(jz.gt.13)jz=13
      if(ja.gt.20)ja=20
c here store sum of cross sections in mb versus z,a indices
c and excitation of fissioning nucleus in 2 MeV wide bin widths
c (to reduce array size).Note that we also have angular momentum
c of fissioning nucleus,which could be used as additional index,
c of interest if we are to calculate fission isomers.This has not
c yet been done due to array size considerations.
      popfiz(jz,ja,iu)=popfiz(jz,ja,iu)+fis1
      pofizbig(jz,ja,iu,jf)=pofizbig(jz,ja,iu,jf)+fis1
      return
      end

c      subroutine fizfrag(emax,spinr)
c     common/pl3/jx,nepr,eb(100),rcsp(100),zee,amass,plex
c     common/fizzc/popfiz(13,20,150),pofizbig(13,20,150,30)
c      common/spinc/psp(90,90,181)
c      common/pmasb/pmass(100,150,30)
c here emax is the excitation of a primary fission
c fragment;afrag is the primary mass number,and spinr
c is the initial angular momentum of the fissioning nucleus.
c loop thru array pofizbig(z,a,E,j) for fissioning nuclides,call
c Mebel routine to get fragments and internal excitation for fission
c of each nucleus,then use algorithm to estimate angular momenta of
c fission fragments.
c the simple algorithm divides fissioning nucleus I between frags prop.
c to moments of inertia of light and heavy.An algo based on sf of Cf252
c fragment spins-mass and energy dependent-is used to estimate spins for
c frags from saddle point modes;the two frag spin contributions are then
c assumed to couple with 2I+1 final spin weighting probability.
c The experimental Cf252 data used for algo was presented inPhys.Rev.C65,
c @2002,'Role of bending mode....',by T.M.Shneidman et al.
c begin loops:
c      do 100 jz=1,13
c      z=zee+1-float(jz)
c      do 200 ja=1,20
c      a=amass+2.-float(jz)-float(ja)
c23456789012345678901234567890123456789012345678901234567890123456789012
c      do 300 iu=1,150
c      e=2.*float(iu)-1.
c      do 400 jf=1,30
c      if(pofizbig(jz,ja,iu,jf).eq.0.)go to 400
c      crmini=pofizbig(jz,ja,iu,jf)/1000.
c       do 450 nev=1,1000
c      call mebel
c     call mebel(z,a,e,a1,z1,a2,z2,ex1,ex2,afin1,afin2,zfin1,zfin2)
c we sent the Z,A,E of fissioning nucleus,receive primary fragment
c Af,Zf,and Ef,and final fragment masses and charges from Mebel routines
c      a153=a1**(5/3)
c      a253=a2**(5/3)
c      al1=spinr*a153/(a153+a253)
c      al2=spinr-al1
c above apportion fissioning nucleus spin between the two fragments
c       jl=al1+0.5
c       jh=al2+0.5
c      afrag=a1
c      afragf=afin1
c      zfin=zfin1
c      do 500 ifrag=1,2
c  loop over light frag(1) and heavy(2)
c      if(ifrag.eq.2)then
c       afrag=a2
c       afragf=afin2
c       zfin=zfin2
c       jl=jh
c       endif
c      lf=0.0352*afrag+0.00070*afrag*emax
c       if(jl.eq.lf)then
c       i1=jl
c       i2=lf
c       else
c       i2=min(jl,lf)
c       i1=max(jl,lf)
c       endif
c 510   continue
c the position of i11 and i22 are likely interchanged here
c       i11=abs(i1-i2)+1
c       i22=i2+i1+1
c       x=rndm(-1.)
c          itran=0
c       do 600 kf=i11,i22
c       if(itran.eq.1)go to 600
c       if(x.le.psp(i11,i22,kf))then
c        lfrag=kf-1
c        itran=1
c        endif
c 600  continue
c 601  continue
c      in=afragf-zfin
c      iz=zfin
c      pmass(iz,in,kf)=pmass(iz,in,kf)+crmini
c store on neutron no. rather than mass no. to reduce storage needs
c 500  continue
c 450  continue
c 400  continue
c 300  continue
c 200  continue
c 100  continue
c      return
c      end
c
      subroutine spincpl
      common/spinc/psp(90,90,181)
c this sub calcs a table psp() to random couple I+l
      do 1 i1=1,90
      do 2 i2=1,90
      do 3 i3=1,181
      psp(i1,i2,i3)=0.
    3 continue
    2 continue
    1 continue
      do 10 i1=1,90
      l1=i1-1
      do 20 i2=1,i1
      l2=i2-1
      tot=(2*l1+1)*(2*l2+1)
      prob=0.
      do 30 i4=l1-l2+1,l1+l2+1
      i3=i4-1
      prob=((i3+i3+1))/tot+prob
c because i3 is one greater as index than vector,above reduced by 2
      psp(i1,i2,i4)=prob
   30 continue
   20 continue
   10 continue
      return
      end
     
      subroutine tfis(a,z,u,px,py,pz)
      implicit doubleprecision(a-h,o-z)
c     CHARACTER*60 FINP,FOUT
      common/tot/totke
      common/fis1c/xfis1
      common/azfrag/afr1,afr2,zfr1,zfr2,ufr1,ufr2
      common /exims/ wapsm(0:150,0:250)
      common /fiss4/ shell(150,250)
      common/rand2/iseed
      common/snglpre/lufr1,lufr2,lafr1,lafr2,lzfr1,lzfr2
      common/snglpri/ia1,ia2,iz1,iz2
      common/fissngl/dissa(200),dissu(200),dissz(200)
      dimension pcn(3),pfr1(3),pfr2(3)
c
       pcn(1)=px
       pcn(2)=py
       pcn(3)=pz
       fiscr=xfis1
         f=fiscr
C...Read random number seed and set up
      iseed=100000
      ran=ran1(iseed)
        call fisgem(u,a,z,pcn,
     1  afr1,zfr1,afr2,zfr2,ufr1,ufr2,pfr1,pfr2,nofis)
        if(nofis.eq.0)  then
              continue
        else 
          ia1=idint(afr1+0.1)
          ia2=idint(afr2+0.1)
           ua1=ufr1/afr1
           ua2=f*ufr2/afr2
            if(ia1.ge.1.and.ia1.le.200)  dissa(ia1)=dissa(ia1)+f
            if(ia2.ge.1.and.ia2.le.200)  dissa(ia2)=dissa(ia2)+f
            if(ia1.ge.1.and.ia1.le.200)  dissu(ia1)=dissu(ia1)+ua1
            if(ia2.ge.1.and.ia2.le.200)  dissu(ia2)=dissu(ia2)+ua2
          iz1=idint(zfr1+0.1)
          iz2=idint(zfr2+0.1)
          if(iz1.ge.1.and.iz1.le.100)  dissz(iz1)=dissz(iz1)+f
        if(iz2.ge.1.and.iz2.le.100)  dissz(iz2)=dissz(iz2)+f
        endif
      return
      end
C----------------------------------------------------------------------
C----------------------------------------------------------------------
      subroutine fisgem(u,a,z,pnuc,afis1,zfis1,afis2,zfis2,
     &                             ufis1,ufis2,pf1,pf2,ifiss)
C/////////////////////////////////////////////////////////////////////
C****************This routine is originally in LAHET code**************
C  FIS
c      Pick post fission parameters such as mass, charge, kinetic energy
c      and excitation energy 
C=====================================================================
C <variables> 
C     a    :   mass of nucleus before fission                     (IN)
C     z    :   charge  of nucleus before fission                  (IN)
C     u    :   excitation energy of nucleus before fission        (IN)
C  pnuc    :   momentum of nucleus before fission                 (IN)
C ifiss    : 1=fission occur, 0=no fission                       (OUT)
C  zfis1,2 :  charge of the fission fragment                     (OUT)
C  afis1,2 :  mass of the fission fragment                       (OUT)
C  ufis1,2 :  excitation energy of the fission fragment          (OUT)
c  pf1,2   :  recoil momentum of the fission fragment            (OUT)

C/////////////////////////////////////////////////////////////////////

      implicit doubleprecision(a-h,o-z)                                 
      parameter (amu=931.494d0)
      parameter (pi=3.1415926535898d0)
      parameter (inn=180, izz=98)                                       
      common/tot/totke
      common /cook1/ sz(izz), sn(inn), pz(izz), pn(inn)                 
      common/barfis/darf
      common/rand2/iseed
cc                                               
      dimension evodba(4)
      dimension pnuc(3),pf1(3),pf2(3),vres(3)
c     following data for assymetric vs symetric picking                 
      data evodba/18.8d0,18.1d0,18.1d0,18.5d0/
      data sigmaa, aamean /6.5d0,140.d0/                  
c     asym gauss width,1.0/level density parameter and high mass mean.  
      data afact, zfact, enmass /8.071323d0,0.782354d0,939.56563d0/     
c     mass diff from 1 amu for neutron,diff p & n masses and n mass     
      data af1, af2, af3 /0.2185024d0,16.70314d0,321.175d0/
c      data ifis/1/
      data ifis/0/
c          ifis=0, Furihata; ifis=1, Atchison
c***********************************************************************
c                                                                       
c       pick the masses                                                 
c                                                                       
c***********************************************************************
c     write(*,*)'from fisgem routine,a,z,u from tfis',a,z,u
c     write(*,*)'from fisgem pnuc=',pnuc
      jz=int(z)
      ja=int(a)
      jn=ja-jz
      se=energy(jz,ja-1)+energy(0,1)-energy(jz,ja)
      se=se+pz(jz)+pn(jn) 
Calculate fission barrier
      x=z*z/a
      ef=x*(af1*x-af2)+af3+se                                           
      ef=darf
      if (ef.gt.u) then    !Excited energy is below fission barrier
        ifiss=0
        return
      endif
      nck2=0  
      ifiss=1                                                          
   10 continue                                                          
      nck2=nck2+1                                                       
      if (nck2.gt.10) then                                              
c***********************************************************************
c     fission failure                                                   
c***********************************************************************
        ifiss=0                                                 
        return                                                          
      endif                                                             
      temp=z*z/a                                                        
      if (temp.le.35.d0) then                                           
        if (jz.le.88) go to 20                                          
      elseif (u.le.62.d0) then                                          
c***********************************************************************
c   high z fission mass distribution. competition for sym vs assym      
c  simple symmetric to assymetric data fit                              
c***********************************************************************
        arg=-0.36d0*u                                                   
        arg=4.87d+03*dexp(arg)                                         
        proba=arg/(1.d0+arg)                                            
        if (ran1(iseed).le.proba) then                                 
c***********************************************************************
c      assymetric fission                                               
c***********************************************************************
          a1=gaussn(aamean,sigmaa)                                      
          go to 30                                                      
        endif                                                           
      endif                                                             
c***********************************************************************
c                                                                       
c     find assymetric barrier for width computation                     
c     assymetric barrier from seaborg and vandenbosch                   
c     phys.rev. 88,507 (1952) phys.rev. 110,507 (1958)                  
c     find if ee eo oe or oo nucleus                                    
c     na=1 odd-odd,2 even-odd,3 odd-even,4 even-even                    
c                                                                       
c***********************************************************************
      if(a.gt.240) then        ! furi May 4, 1999
       in=ja-jz
       na=1                                                             
       if (jz.eq.2*(jz/2)) na=na+1                                      
       if (in.eq.2*(in/2)) na=na+2                                      
       temp=z*z/a                                                       
       ef=evodba(na)-0.36d0*temp                                        
       ef=darf
c dec'03 replace fission barriers with alice generated values
      endif                    ! furi May 4, 1999
   20 continue                                                          
C  furi ... May 4, 1999
      if(ifis.ne.0) then
       upr=u-ef                                                         
       if (upr.gt.1.d+02) upr=1.d+02                                    
       sigmas=0.425d0*(upr*(1.d0-0.005d0*upr)+9.35d0)                   
      else 
       xx=z*z/a                                                        
       bf=efms(z,a)
c replace with alice barrier
       bf=darf
       if(bf.lt.0)bf=ef
       upr=min(u-bf,450.0d0)
       sigmas=0.122d0*xx**2.d0-7.77*xx+134.0d0 +3.32d-2*upr
      endif
c***********************************************************************
c     sigmas is symmetrin fission mass width.taken from systematic of   
c     neuzil & fairhall phys.rev. 129,2705,(1963)                       
c                                                                       
c      low z fission is always symmetric                                
c                                                                       
c      high z fission is sometimes.loop back here from 1 loop if        
c      symmetric fission predicted.                                     
c***********************************************************************
C  furi ... May 4, 1999
      amean=0.5d0*a                                                     
      a1=gaussn(amean,sigmas)   
   30 continue                                                          
c***********************************************************************
c      1 loop for assymmetrin fission returns to here.                  
c***********************************************************************
      afis1=dint(a1)                                                  
c***********************************************************************
c    check for low final a                                              
c***********************************************************************
      if (afis1.lt.5.d0) afis1=5.d0                              
      if (a-afis1.lt.5.d0) afis1=a-5.d0
      afis2=a-afis1                                                 
c***********************************************************************
c        pick the charge                                                
c***********************************************************************
      z1=65.5d0*afis1/(131.d0+afis1**0.666667d0)                    
      z2=65.5d0*afis2/(131.d0+afis2**0.666667d0)                    
      z1=z1+.5d0*(z-z1-z2)                                              
c***********************************************************************
c       we use constant charge density with a 2 unit gaussian smearing  
c***********************************************************************
c   sigma = 0.75;   z1=gaussn(z1,dph)
c                                by furi 18/DEC/1997
c***********************************************************************
      if(ifis.ne.0) then
       sigz=2.0d0
      else
       sigz=0.75d0
      endif
      z1=gaussn(z1,sigz)                                                
      zfis1=dint(z1)                                                  
      zfis2=z-zfis1                                                 
c***********************************************************************
c   check for reasonable z a combinations.                              
c***********************************************************************
      if (zfis1.ge.afis1) go to 10                                  
      if (zfis2.ge.afis2) go to 10                                  
      if (zfis1.lt.1.d0) go to 10                                      
      if (zfis2.lt.1.d0) go to 10                                      
c***********************************************************************
c      compute binding energy and actual masses of fragments            
c***********************************************************************
      be0=afact*a-zfact*z-energy(jz,ja)                                 
      rm0=enmass*a-zfact*z-be0                                          
      iaf=int(afis1)
      izf=int(zfis1)
      be1=afact*afis1-zfact*zfis1-energy(izf,iaf)           
      rm1=enmass*afis1-zfact*zfis1-be1                              
      iaf=int(afis2)
      izf=int(zfis2)
      be2=afact*afis2-zfact*zfis2-energy(izf,iaf)
      rm2=enmass*afis2-zfact*zfis2-be2                              
c***********************************************************************
c      pick recoil kinetic energy.use systematic of ...........         
c      unik et.al. proc.3rd iaea symp.on phy.& chem. fision,rochester.vo
c***********************************************************************
      if(ifis.ne.0)then
       totkm=0.13323d0*z*z/a**0.33333333d0-11.4d0                       
      else
C... 4 May 1999 by furi
       x=z*z/a**0.33333d0
       if(x.le.900.d0) then
        totkm=0.131d0*x
       else if(x.le.1800.d0) then
        totkm=0.104d0*x+24.3d0
       else
        write(*,*)'Error in totk in subroutine fiss'
       endif
      endif
c***********************************************************************
c      use a width of 15% value at half height.                         
c***********************************************************************
      if(ifis.ne.0)then
       sigmak=0.084d0*totkm                                             
C... 4 May 1999 by furi
      else
       if(x.lt.1000)then
        sigmak=86.5d0
       else if (x.lt.1800.d0)then
        sigmak=5.70d-4*(x-1000.d0)**2.d0+86.5d0
       else
        write(*,*)'Error in sigmak in subroutine fiss'
       endif
      sigmak=dsqrt(sigmak)
      endif
c***********************************************************************
c   check event is energetically possible                               
c***********************************************************************
      temp2=u+be1+be2-be0                                               
      nck=0                                                             
   40 continue                                                          
      totke=gaussn(totkm,sigmak)                                        
      if (nck.gt.10) go to 10                                           
      nck=nck+1                                                         
      if (totke.gt.temp2) go to 40                                      
c***********************************************************************
c      pick excitation from equidistribution of original plus energy bal
c***********************************************************************
      temp=(temp2-totke)/a                                              
      ufis1=afis1*temp                                              
      ufis2=afis2*temp                                              
cccc                                                                    
c     find total masses, including excitation energies, at evap time    
cccc                                                                    
      amcf=rm0+u                                                        
      amc1=rm1+ufis1                                                  
      amc2=rm2+ufis2                                                  
      amdiff=amcf-amc1-amc2                                             
      if (amdiff.lt.0.d0) go to 40                                      
c***********************************************************************
c     amdiff= ekin should be satisfied                                  
c***********************************************************************
c     note that the coding  following is 'dead code' relativistic ff kinematics
c     into cm  and  lab for frags 1  and 2; this is now done in sub fizkinem
c     non-relativistically
c------------------------------------------------------------------------
      return
c... Velocity of pre-fission nucleus 
      pres=dsqrt(pnuc(1)**2+pnuc(2)**2+pnuc(3)**2)
      eres=dsqrt(pres**2+(a*amu)**2)      
      vres(1)=pnuc(1)/eres
      vres(2)=pnuc(2)/eres
      vres(3)=pnuc(3)/eres 
c... Momentum of CM system
c        9/4/2014 c this stmt     return
        return
      redm=amu*(afis1*afis2)/a
      pcm=dsqrt(totke**2+2.d0*redm*totke)
      th  = dacos( 2.0d0 * ran1(iseed) - 1.0d0 )
      ph  = 2.0d0 * pi  * ran1(iseed)

      pcmx = pcm* dsin(th) * dcos(ph) 
      pcmy = pcm* dsin(th) * dsin(ph)
      pcmz = pcm* dcos(th)

C... Velocity of fission fragment 1 in CM system
      v1x=pcmx/dsqrt(pcm**2 + (amu*afis1)**2)
      v1y=pcmy/dsqrt(pcm**2 + (amu*afis1)**2)
      v1z=pcmz/dsqrt(pcm**2 + (amu*afis1)**2)

C... Velocity of fission fragment 2 in CM system
      v2x=-pcmx/dsqrt(pcm**2 + (amu*afis2)**2)
      v2y=-pcmy/dsqrt(pcm**2 + (amu*afis2)**2)
      v2z=-pcmz/dsqrt(pcm**2 + (amu*afis2)**2)

C... Boost to Lab system
      v1x=v1x+vres(1)
      v1y=v1y+vres(2)
      v1z=v1z+vres(3)

      v2x=v2x+vres(1)
      v2y=v2y+vres(2)
      v2z=v2z+vres(3)
      
C... Kinetic Energy and velocity of fission fragment 1 in Lab system
      gam=1.d0/dsqrt(1.d0-v1x**2-v1y**2-v1z**2)
      er1= afis1*amu *(gam -1.d0)
      pm1=dsqrt(er1*(er1+2.*afis1*amu))
      pf1(1)=pm1*v1x/dsqrt(v1x**2+v1y**2+v1z**2)
      pf1(2)=pm1*v1y/dsqrt(v1x**2+v1y**2+v1z**2)
      pf1(3)=pm1*v1z/dsqrt(v1x**2+v1y**2+v1z**2)

C... Kinetic Energy and velocity of fission fragment 2 in Lab system
      gam=1.d0/dsqrt(1.d0-v2x**2-v2y**2-v2z**2)
      er2= afis2*amu *(gam -1.d0)
      pm2=sqrt(er2*(er2+2.*afis2*amu))
      pf2(1)=pm2*v2x/dsqrt(v2x**2+v2y**2+v2z**2)
      pf2(2)=pm2*v2y/dsqrt(v2x**2+v2y**2+v2z**2)
      pf2(3)=pm2*v2z/dsqrt(v2x**2+v2y**2+v2z**2)

      return                                                            
c                                                                       
c  50 format ('---> fission failed: ja=',i5,'  jz=',i5/'-       u=',1pe1
c    10.3,'      erec=',e10.3/'    zfis1=',e10.3,'     afis1=',e10.3/' 
c    2   zfis2=',e10.3,'     afis2=',e10.3)                             
      end                                                               
C----------------------------------------------------------------------
C----------------------------------------------------------------------
C----------------------------------------------------------------------
      function energy (iz,ia)                                          
C/////////////////////////////////////////////////////////////////////
C****************This routine is originally in the LAHET code********
C  ENERGY
c    Calculate excess mass
C=====================================================================
C <variables> 
C     ia   :   mass of nucleus                      (IN)
C     iz   :   charge  of nucleus                   (IN)
C  energy  :   Excess mass [MeV]                   (OUT)
C/////////////////////////////////////////////////////////////////////
      implicit doubleprecision(a-h,o-z)                                 
      common /exims/ wapsm(0:150,0:250)
      parameter (inn=180, izz=98)                                       
      common /cook1/ sz(izz), sn(inn), pz(izz), pn(inn)                
      cam(a,z,b)=8.071323d0*a-0.782354d0*z-17.0354d0*a*(1.d0-1.84619d0
     1 *(a-2.d0*z)**2/a**2)+25.8357d0*b**2*(1.d0-1.71219d0*(a-2.d0*z)
     2 **2/a**2)*(1.d0-0.62025d0/b**2)**2+0.779d0*z*(z-1.d0)*
     3 (1.d0-1.5849d0/b**2+1.2273d0/a+1.5772d0/(a*b))/b-
     4 0.4323d0*dexp(1.3333333333d0*dlog(z))*(1.d0-0.57811d0/b-
     5 0.14518d0/b**2+0.49597d0/a)/b 

      energy=1.d10
      in=ia-iz
      if (iz.ge.0.and.iz.le.150.and.ia-iz.ge.0.and.ia-iz.le.250) then   
       energy=wapsm(iz,ia-iz)                                           
      if(energy.eq.0.d0) goto 10
      return
      endif

Cameron's mass formula
 10   if(iz.eq.6.and.ia.eq.12) return
      if(iz.eq.0) return
      n=ia-iz                                                          
      a=dble(ia)
      z=dble(iz)
      a3=dexp(0.3333333333d0*dlog(a))                                  
      energy=cam(a,z,a3)+sz(iz)+sn(n)                               
      return                                                            
      end                                                               
C----------------------------------------------------------------------
CC----------------------------------------------------------------------
      function efms(z,a)
C/////////////////////////////////////////////////////////////////////
C  EFMS
C  Fission barrier given by Myer & Swaiteski (PRC60,014606,1999)
C=====================================================================
C <variables> 
C     a   :   the mass of a fissioning nucleus      (IN)
C     z   :   the charge of a fissioning nucleus    (IN)
C   efms  :   fission barrier  [MeV]                (OUT)
C/////////////////////////////////////////////////////////////////////

      implicit doubleprecision(a-h,o-z)                                 
      common /fiss4/ shell(150,250)
      parameter(x0=48.5428d0, x1=34.15d0)
      f1(t)=1.99749d-4*(x0-t)**3
      f2(t)=5.95553d-1-0.124136*(t-x1)

      efms=0.d0
      iz=int(z)
      ia=int(a)

      c=1.9+(z-80)/75
      ai=1.-2*(z/a)
      xx=1-c*ai**2
      ss=a**.66667*xx
      x=z**2/a/xx
      if(ia-iz.le.0.or.ia-iz.gt.250.or.iz.gt.150.or.iz.lt.1) then
       sh=0.d0
      else
       sh=shell(iz,ia-iz)
      endif
      if(x.ge.x1.and.x.le.x0) then
       efms=ss*f1(x)-sh
      else if(x.ge.20.and.x.lt.x0)then
       efms=ss*f2(x)-sh
      else
       efms=-1.0
      endif

      return
      end
C----------------------------------------------------------------------
C----------------------------------------------------------------------
C----------------------------------------------------------------------
       function gaussn (xmean,sd)                                       
C/////////////////////////////////////////////////////////////////////
C****************This routine is originally from the LAHET code********
C  GAUSSN
c      Gaussian randum number gemerator
c  compute random gaussian number for given                           
c  mean and s.d.                                                      
c  uses mean of sum of 12 uniform r.n"s                               
c                                                                     
C/////////////////////////////////////////////////////////////////////
      implicit doubleprecision(a-h,o-z)                                 
      common/rand2/iseed                               

      a=0.d0
      do 10 n=1,12
        a=a+ran1(iseed)
 10   continue
      gaussn=(a-6.d0)*sd+xmean                                          
      return                                                            
      end                                                               
C----------------------------------------------------------------------
C ----------------------------------------------
C ----------------------------------------------
      function ran1(iseed)
C/////////////////////////////////////////////////////////////////////
C  RAN1
C    Random number generator
C=====================================================================
C <variables> 
C  iseed   :  Seed of random number.  <0 to initialize the sequence
C/////////////////////////////////////////////////////////////////////

      implicit doubleprecision(a-h,o-z)                                 
      common/rand/rn(97),ix1,ix2,ix3
      parameter (m1=259200,ia1=7141,ic1=54733,rm1=1./m1)
      parameter (m2=134456,ia2=8121,ic2=28411,rm2=1./m2)
      parameter (m3=243000,ia3=4561,ic3=51349)
      data iff /0/

      if(iseed.lt.0.or.iff.eq.0)then
       iff=1
       ix1=mod(ic1+iseed,m1)
       ix1=mod(ia1*ix1+ic1,m1)
       ix2=mod(ix1,m2)
       ix1=mod(ia1*ix1+ic1,m1)
       ix3=mod(ix1,m3)
       do 10 j=1,97
         ix1=mod(ia1*ix1+ic1,m1)
         ix2=mod(ia2*ix2+ic2,m2)
         rn(j)=(dble(ix1)+dble(ix2)*rm2)*rm1
 10    continue
       iseed=1
      endif
      ix1=mod(ia1*ix1+ic1,m1)
      ix2=mod(ia2*ix2+ic2,m2)
      ix3=mod(ia3*ix3+ic3,m3)
      j=1+(97*ix3)/m3
      ran1=rn(j)
      rn(j)=(dble(ix1)+dble(ix2)*rm2)*rm1
      return
      end
C ----------------------------------------------
C ----------------------------------------------
      block data bdcook                                                
C/////////////////////////////////////////////////////////////////////
C****************This routine is originally in LAHET code**************
C  BD01
c      DATA of paring energy, shell correction
c       
C=====================================================================
C <variables> 
c    sz  :  shell correction for each Z
c    sn  :  shell correction for each N
c    pz  :  paring energy for each Z
c    pn  :  paring energy for each N
C  
C/////////////////////////////////////////////////////////////////////

      implicit doubleprecision(a-h,o-z)                                 

      parameter (dp0=0.d0, dp1=1.d0, dp2=2.d0, dp3=3.d0, dp4=4.d0, dph=.
     1 5d0, dp5=5.d0, dp10=1.d1, dpth=dp1/dp3, dppi=3.1415926535898d0,  
     2 dp2th=dp2/dp3)                                                   
c                                                 
      parameter (inn=180, izz=98)                                       
      common /cook1/ sz(izz), sn(inn), pz(izz), pn(inn)               
c                                                                       
c     data tables of cook et. al.  aaec/tm392, supplimented by g+c.     
c                                                                       
      data sz /8*dp0,-.11d0,-.81d0,-2.91d0,-4.17d0,-5.72d0,-7.8d0,-8.97d
     1 0,-9.7d0,-10.1d0,-10.7d0,-11.38d0,-12.07d0,-12.55d0,-13.24d0,-13.
     2 93d0,-14.71d0,-15.53d0,-16.37d0,-17.36d0,-18.6d0,-18.7d0,-18.01d0
     3 ,-17.87d0,-17.08d0,-16.6d0,-16.75d0,-16.5d0,-16.35d0,-16.22d0,-16
     4 .41d0,-16.89d0,-16.43d0,-16.68d0,-16.73d0,-17.45d0,-17.29d0,-17.4
     5 4d0,-17.82d0,-18.62d0,-18.27d0,-19.39d0,-19.91d0,-19.14d0,-18.26d
     6 0,-17.4d0,-16.42d0,-15.77d0,-14.37d0,-13.91d0,-13.1d0,-13.11d0,-1
     7 1.43d0,-10.89d0,-10.75d0,-10.62d0,-10.41d0,-10.21d0,-9.85d0,-9.47
     8 d0,-9.03d0,-8.61d0,-8.13d0,-7.46d0,-7.48d0,-7.2d0,-7.13d0,-7.06d0
     9 ,-6.78d0,-6.64d0,-6.64d0,-7.68d0,-7.89d0,-8.41d0,-8.49d0,-7.88d0,
     $ -6.3d0,-5.47d0,-4.78d0,-4.37d0,-4.17d0,-4.13d0,-4.32d0,-4.55d0,-5
     $ .04d0,-5.28d0,-6.06d0,-6.28d0,-6.87d0,-7.20d0,-7.74d0/           
      data (sn(i),i=1,110) /8*dp0,10.3d0,5.66d0,6.8d0,7.53d0,7.55d0,7.21
     1 d0,7.44d0,8.07d0,8.94d0,9.81d0,10.6d0,11.39d0,12.54d0,13.68d0,14.
     2 34d0,14.19d0,13.83d0,13.5d0,13.d0,12.13d0,12.6d0,13.26d0,14.13d0,
     3 14.92d0,15.52d0,16.38d0,17.16d0,17.55d0,18.03d0,17.59d0,19.03d0,1
     4 8.71d0,18.8d0,18.99d0,18.46d0,18.25d0,17.76d0,17.38d0,16.72d0,15.
     5 62d0,14.38d0,12.88d0,13.23d0,13.81d0,14.9d0,14.86d0,15.76d0,16.2d
     6 0,17.62d0,17.73d0,18.16d0,18.67d0,19.69d0,19.51d0,20.17d0,19.48d0
     7 ,19.98d0,19.83d0,20.2d0,19.72d0,19.87d0,19.24d0,18.44d0,17.61d0,1
     8 7.1d0,16.16d0,15.9d0,15.33d0,14.76d0,13.54d0,12.63d0,10.65d0,10.1
     9 d0,8.89d0,10.25d0,9.79d0,11.39d0,11.72d0,12.43d0,12.96d0,13.43d0,
     $ 13.37d0,12.96d0,12.11d0,11.92d0,11.d0,10.8d0,10.42d0,10.39d0,9.69
     $ d0,9.27d0,8.93d0,8.57d0,8.02d0,7.59d0,7.33d0,7.23d0,7.05d0,7.42d0
     $ ,6.75d0,6.6d0,6.38d0/                                            
      data (sn(i),i=111,150) /6.36d0,6.49d0,6.25d0,5.85d0,5.48d0,4.53d0,
     1 4.3d0,3.39d0,2.35d0,1.66d0,.81d0,0.46d0,-.96d0,-1.69d0,-2.53d0,-3
     2 .16d0,-1.87d0,-.41d0,.71d0,1.66d0,2.62d0,3.22d0,3.76d0,4.1d0,4.46
     3 d0,4.83d0,5.09d0,5.18d0,5.17d0,5.1d0,5.01d0,4.97d0,5.09d0,5.03d0,
     4 4.93d0,5.28d0,5.49d0,5.50d0,5.37d0,5.30d0/                       
      data pz /dp0,5.44d0,dp0,2.76d0,dp0,3.34d0,0.,2.7d0,dp0,1.90d0,dp0,
     1 2.12d0,dp0,2.09d0,dp0,1.62d0,dp0,1.62d0,dp0,1.83d0,dp0,1.73d0,dp0
     2 ,1.35d0,dp0,1.54d0,dp0,1.28d0,0.26d0,0.88d0,0.19d0,1.35d0,-.05d0,
     3 1.52d0,-.09d0,1.17d0,.04d0,1.24d0,0.29d0,1.09d0,.26d0,1.17d0,.23d
     4 0,1.15d0,-.08d0,1.35d0,0.34d0,1.05d0,.28d0,1.27d0,dp0,1.05d0,dp0,
     5 1.d0,.09d0,1.2d0,.2d0,1.4d0,.93d0,1.d0,-.2d0,1.19d0,.09d0,.97d0  
     6 ,dp0,.92d0,.11d0,.68d0,.05d0,.68d0,-.22d0,.79d0,.09d0,.69d0,.01d0
     7 ,.72d0,dp0,.4d0,.16d0,.73d0,dp0,.46d0,.17d0,.89d0,dp0,.79d0,dp0,.
     8 89d0,dp0,.81d0,-.06d0,.69d0,-.2d0,.71d0,-.12d0,.72d0,dp0,.77d0/  
      data (pn(i),i=1,125) /dp0,5.98d0,dp0,2.77d0,dp0,3.16d0,dp0,3.01d0 
     1 ,dp0,1.68d0,0.,1.73d0,0.,2.17d0,dp0,1.67d0,dp0,1.86d0,dp0,2.04d0 
     2 ,dp0,1.64d0,dp0,1.44d0,dp0,1.54d0,dp0,1.3d0,dp0,1.27d0,dp0,1.29d0
     3 ,.08d0,1.41d0,-.08d0,1.5d0,-.05d0,2.24d0,-.47d0,1.43d0,-.15d0,1.4
     4 4d0,.06d0,1.56d0,.25d0,1.57d0,-.16d0,1.46d0,dp0,.93d0,.01d0,.62d0
     5 ,-.5d0,1.42d0,.13d0,1.52d0,-.65d0,.8d0,-.08d0,1.29d0,-.47d0,1.25d
     6 0,-.44d0,.97d0,.08d0,1.65d0,-.11d0,1.26d0,-.46d0,1.06d0,0.22d0,1.
     7 55d0,-.07d0,1.37d0,0.1d0,1.2d0,-.27d0,.92d0,-.35d0,1.19d0,dp0,1.0
     8 5d0,-.25d0,1.61d0,-.21d0,.9d0,-.21d0,.74d0,-.38d0,.72d0,-.34d0,.9
     9 2d0,-.26d0,.94d0,.01d0,.65d0,-.36d0,.83d0,.11d0,.67d0,.05d0,1.d0,
     $ .51d0,1.04d0,.33d0,.68d0,-.27d0,.81d0,.09d0,.75d0,.17d0,.86d0,.14
     $ d0,1.1d0,-.22d0,.84d0,-.47d0,.48d0,.02d0,.88d0,.24d0,.52d0,.27d0,
     $ .41d0,-.05d0/                                                    
      data (pn(i),i=126,150) /0.38d0,.15d0,.67d0,dp0,.61d0,dp0,.78d0,dp0
     1 ,.67d0,dp0,.67d0,dp0,.79d0,dp0,.6d0,.04d0,.64d0,-.06d0,.45d0,.05d
     2 0,.26d0,-.22d0,.39d0,dp0,.39d0/                                  
c                                                                       
      end                                                               

C ----------------------------------------------
C ----------------------------------------------




       subroutine gem
c begin to modify mashnik/gudima code for use as subs in alice
      implicit doubleprecision(a-h,o-z)                                 
      common /exims/ wapsm(0:150,0:250)
      common /fiss4/ shell(150,250)
      common/rand2/iseed
      common/fissngl/dissa(200),dissu(200),dissz(200)
      do iz=0,150
        do in=0,250
          wapsm(iz,in)=0.d0
        enddo
      enddo
      do 1 j=1,10000
        read(31,*,end=99)iz,ia,wapsm(iz,ia-iz)
        
        if(ia-iz.gt.250.or.iz.gt.150) then
         write(*,*)'ERROR: Mass table reading error in subroutine SETUP'
         write(*,*)' Z or A-Z in mass.tbl is larger than 150 or 250', 
     &        iz,ia-iz
         stop
c fill buffer for ff mass excesses to have lysekil/wapstra
        endif
 1    continue
      write(*,*)' WARNING: The last record read in subroutine SETUP is '
     &     ,iz,ia,wapsm(ia,ia-iz)
      write(*,*)'mass.tbl should be less than 10000 records'
 99   close(31)

C  Read shell effect data for fission barrier calculations
      do i=1,150
        do j=1,250
          shell(i,j)=0.
        enddo
      enddo
c
      do i=1,10000
        read(32,*,end=98)j,k,sh
        if(k-j.lt.1) then
         write(*,*)'Shell table reading error : iz, ia= ',j,k
         stop
        else
         shell(j,k-j)=sh
        endif
      enddo
 98   close(32)
      return
      end
      subroutine fizfr(a1,a2,ux1,ux2,spinr,lfr1,lfr2,zf1,zf2)
      common/pl3/jx,nepr,eb(100),rcsp(100),zee,amass,plex
      common/fizzc/popfiz(13,20,150),pofizbig(13,20,150,30)
      common/spinc/psp(90,90,181)
      common/pmasb/pmass(100,150,30)
      common/libst/a53t(260)
c here emax is the excitation of a primary fission
c fragment;afrag is the primary mass number,and spinr
c is the initial angular momentum of the fissioning nucleus.
c loop thru array pofizbig(z,a,E,j) for fissioning nuclides,call
c Mebel routine to get fragments and internal excitation for fission
c of each nucleus,then use algorithm to estimate angular momenta of
c fission fragments.
c the simple algorithm divides fissioning nucleus I between frags prop.
c to moments of inertia of light and heavy.An algo based on sf of Cf252
c fragment spins-mass and energy dependent-is used to estimate spins for
c frags from saddle point modes;the two frag spin contributions are then
c assumed to couple with 2I+1 final spin weighting probability.
c The experimental Cf252 data used for algo was presented inPhys.Rev.C65,
c @2002,'Role of bending mode....',by T.M.Shneidman et al.
c begin loops:
       j=a1
c     a153=a53t(j)
       k=a2
c     a253=a53t(k)
c     al1=spinr*a153/(a153+a253)
c     al1=spinr*a1/(a1+a2)
c     al2=abs(spinr-al1)
c above apportion fissioning nucleus spin between the two fragments
c try equal spin split
       al1=spinr/2.
       al2=spinr/2.
       jl=al1
       jh=al2
      afrag=a1
      zfrag=zf1
        umax=ux1
      do 500 ifrag=1,2
c  loop over light frag(1) and heavy(2)
      if(ifrag.eq.2)then
       afrag=a2
       zfrag=zf2
       jl=jh
        umax=ux2
       endif
       alf=0.0352*afrag+0.00070*afrag*umax
       lf=alf
c here modify for z50,n82 and 49-51,81-83 products-arbitrary
       fragn=afrag-zfrag
       frnt=abs(fragn-82.)
       frzt=abs(zfrag-50.)
       if(frnt.le.2..and.frzt.le.2.)lf=alf*0.9
       if(frnt.le.1..and.frzt.le.1.)lf=alf*0.8
       if(afrag.eq.132..and.zfrag.eq.50.)lf=alf*0.7
c here modify to remove energy dependence of doubly magic 50,126?
       i2=min(jl,lf)
       i1=max(jl,lf)
       if(i2.lt.1)i2=1
       if(i1.lt.1)i1=1
       continue
       i11=i1+1
       i22=i2+1
       if(i22.lt.1)i22=1
       if(i11.lt.1)i11=1
       x=rndm(-1.)
          itran=0
       imin=i1-i2+1
       if(imin.lt.1)imin=1
       imax=i1+i2+1
       do 600 kf=imin,imax
       if(itran.eq.1)go to 600
       lfrag=i22
       if(x.le.psp(i11,i22,kf))then
        lfrag=kf-1
        itran=1
        endif
 600  continue
        if(ifrag.eq.1)then
       lfr1=lfrag
        else
       lfr2=lfrag
      endif
      continue
 500  continue
      return
      end
c
      subroutine findbeta2s
      common/deform/beta2(106,260)
      do i=1,8000
      read(35,2)iaz,d1,d2,d3,d4,d5,b2,d6,d7,d8,d9,d10,d11,d12,d13,d
      
      iz=iaz/1000
      ia=iaz-1000*iz
      if(iz.gt.106.or.ia.gt.260)go to 3
      beta2(iz,ia)=b2
   3  continue
      enddo
   2  format(1i7,15f7.0)
      return
      end
c
      subroutine fizkinem(tote,af1,af2,ecmf,amasf)
      common/endvparm/iend,lll(6,6,6),ll(6,6)
      common/history/ide(60),idn(60),idith(60),iddum(60),the(60),adelj
     1(60),spinfh(60),enghis(60)
      common/out1/nout
      common/histindx/ihist
      common/nuopt/kinem,iem,kb
      common/evapcon/zof(8),aof(8)
      common/evapconi/lzl(8),laa(8),lnn(8)
      common/eder/ed2,v,vhalf,vsq
      common/fizk/csd,boff
      common/ifizk/jk,ik,nbr,nfr
      common/labcon/prx,pry,prz,pcm,rmas,clab(999,8),dslb(36,999,8)
      common/labcon1/prxf1,pryf1,przf1,prxf2,pryf2,przf2
      common/labcon3/pcmx1,pcmy1,pcmz1,pcmx2,pcmy2,pcmz2
      common/labcon4/v1x,v1y,v1z,v2x,v2y,v2z
      common/incr/ed,tened,efinal
      common/memof/sfcrs(999,8),dfsgp(36,999,4),efn(4,999)
      common/fisn/dfslb(36,999,4),cflab(999,4),dfetl(36,999,4)
      common/tot/totke
      doubleprecision totke
c in this routine we do not continually rotate angles back to
c a lab frame due to isotropicity of emissions.Only final ejectile
c directions are calculated in lab frame - could be done for fission frags.
        aoff=boff
cxxxxxxxxxxxxxxxx
c aoff is mass of light ejectile,csds is the fission cross section
        csds=csd/ed
c//////////////////////////////////////////////////////////////////////
        if(nbr.eq.1)then
c apportion total fission KE between fragments 1 and 2
      afis=af1+af2
      ef1=tote*af2/(af1+af2)
      ef2=tote-ef1
c prx,y,z are the momentum of recoil nucleus which will undergo fission
      vrecx=prx/afis
      vrecy=pry/afis
      vrecz=prz/afis
c above are velocities of recoiling nucleus which fissioned
       thet=acos(2.*rndm(-1.)-1.)
       phi=2.*3.14159*rndm(-1.)
c calc theta and phi for first fragment wrt fissioning /recoil nucleus
        pcmf1=sqrt(2.*ef1*af1)
        pcmf2=-pcmf1
c pcmf1 should =-pcmf2
        pcmx1=pcmf1*sin(thet)*cos(phi)
        pcmy1=pcmf1*sin(thet)*sin(phi)
        pcmz1=pcmf1*cos(thet)
c for fission fragment 1:
c convert to lab velocities by adding fissioning nucleus velocities
        v1x=vrecx+pcmx1/af1
        v1y=vrecy+pcmy1/af1
        v1z=vrecz+pcmz1/af1
c next do secondfragment
       pcmx2=-pcmx1
       pcmy2=-pcmy1
       pcmz2=-pcmz1
c for fission fragment-
c convert to lab velocities by adding fissioning nucleus velocities
        v2x=vrecx+pcmx2/af2
        v2y=vrecy+pcmy2/af2
        v2z=vrecz+pcmz2/af2
c these velocities would be used to calculate lab system fission fragment
c energies
      endif
c//////////////////////////////////////////////////////////////////////
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      if(nbr.eq.2)then
c
c here we arrive having emitted a particle of mass aoff with c.m. energy
c ecmf from a frag having final mass amasf and momentum pcm 
c
c select theta and phi for the particle emitted from the fission fragment
       thet=acos(2.*rndm(-1.)-1.)
       phi=2.*3.14159*rndm(-1.)
      pcmej=sqrt(2.*aoff*ecmf)
      pcmejx=pcmej*sin(thet)*cos(phi)
      pcmejy=pcmej*sin(thet)*sin(phi)
      pcmejz=pcmej*cos(thet)
      vfx=pcmejx/aoff
      vfy=pcmejy/aoff
      vfz=pcmejz/aoff
c now add the velocity of the fissioning nucleus
c
c------------------------------------------------------------------------------
       if(nfr.eq.1)then
c  do for either frag 1 or frag 2
       vf1x=vfx+v1x
       vf1y=vfy+v1y
       vf1z=vfz+v1z
c calculate lab ke of emitted particle and lab angle
      ek=0.5*aoff*(vf1x**2+vf1y**2+vf1z**2)
      vxy=sqrt(vf1x**2+vf1y**2)
      thetlab=atan2(vxy,vf1z)      
c now decrease momentum of emitting fragment
c now adjust the fission frag velocity after emission
       pcmx1=pcmx1-pcmejx
       pcmy1=pcmy1-pcmejy
       pcmz1=pcmz1-pcmejz
       endif
c---------------------------------------------------------------------------------
c=========================================================================
       if(nfr.eq.2)then
       vf2x=vfx+v2x
       vf2y=vfy+v2y
       vf2z=vfz+v2z
c calculate lab ke of emitted particle and lab angle
      ek=0.5*aoff*(vf2x**2+vf2y**2+vf2z**2)
      vxy=sqrt(vf2x**2+vf2y**2)
      thetlab=atan2(vxy,vf2z)      
c now decrease momentum of emitting fragment
       pcmx2=pcmx2-pcmejx
       pcmy2=pcmy2-pcmejy
       pcmz2=pcmz2-pcmejz
            endif
c=========================================================================
c do emitted energy in ff frame jk
c here check that element of cross section added is fission cs for this frag.
        ik=(ecmf+ed)/ed
        if(ik.lt.1)ik=1
c------------------------------------------
              if(ik.gt.999.or.jk.gt.4)then
       write(*,*)'ik,jk',ik,jk
        return
          endif
c------------------------------------------
       sfcrs(ik,jk)=sfcrs(ik,jk)+csds
       ith=thet*11.459155+1.
       if(ith.ge.37)ith=36
       dfsgp(ith,ik,jk)=dfsgp(ith,ik,jk)+csds
c now do lab system store
c      calc new angle lab thetlab
       if(iend.eq.1)call endf(thet,ecmf,ith,ik,jk )
       ithl=thetlab*11.459155+1.
        if(ithl.gt.36)ithl=36
       mk=(ek+ed)/ed
       if(mk.lt.1)mk=1
c---------------------------------------------
           if(mk.gt.999)then
           write(*,*)'mk=',mk
           return
           endif
c---------------------------------------------
       cflab(mk,jk)=cflab(mk,jk)+csds
       dfslb(ithl,mk,jk)=dfslb(ithl,mk,jk)+csds
       if(iend.eq.2)call endf(thetlab,ek,ithl,mk,jk )
         endif
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
      return
      end
c
      subroutine evap1(amasf,zfis,ufis,fjx,tote)
c nov'03 make this sr from evap to de-excite fission fragments to ground
c or state which cannot emit nucleons/alphas
c we will identify fragments having isomers as those for which 'isoid(Z,A)>1;
c this array will give a unique number 'n' to the isomer,which will relate to
c the A,Z thru arrays AS(n),ZS(n).The array of excitation vs.  spin(as 2J+1)
c will be pfpej(n,2J+1,E)
 
      doubleprecision dafr1,dafr2,dzfr1,dzfr2,dufr1,dufr2
      common/PRODA/PROD(7,20)
      common/history/ide(60),idn(60),idith(60),iddum(60),the(60),adelj
     1(60),spinfh(60),enghis(60)
      common/out1/nout
      common/clustr/m3set
      common/histindx/ihist
      common/endvparm/iend,lll(6,6,6),ll(6,6)
      common/tot/totke
       doubleprecision totke
      common/gam/io,la,ppgam(999),gmspec(999),beta,ams
      common/eder/ed2,v,vhalf,vsq
      common/evapcon/zof(8),aof(8)
      common/evapconi/lzl(8),laa(8),lnn(8)
      common/labcon3/pcmx1,pcmy1,pcmz1,pcmx2,pcmy2,pcmz2
      common/labcon4/v1x,v1y,v1z,v2x,v2y,v2z
      common/fisn/dfslb(36,999,4),cflab(999,4),dfetl(36,999,4)
      common/azfrag/dafr1,dafr2,dzfr1,dzfr2,dufr1,dufr2
      common/lab3f/slavf(100,4,999)
      common/fizinv/sif(100,4,999)
      common/emin/eminf(106,160)
      common/frags/xmasf(160,106),bq(106,160,4),pf(106,160)
     1,pfp(106,160,100),ppfi(106,160),ppfz(106),ppfa(260)
     2,ppfaz(106,260)
      common/fizstuf/pc(3)
      REAL*8 pc
      common/ldf/gamtem(9999)
      REAL*8 gamtem
      common/krld/powkr(9999),gamkr(9999),powgkr(9999)
      REAL*8 powkr,gamkr,powgkr
      common/shell/shelmol
      common/snglpre/lufr1,lufr2,lafr1,lafr2,lzfr1,lzfr2
      common/gr/rint(1100)
      REAL*8 rint
      common/gams/gamgam,csgam,csgamm(100)
      REAL*8 gam1,gamgam,gamfis,csgamm,gfis,fis1,csgam
       common/fizz/fizbar(106,260),fizcs(100)
      common/outlim/limout,lim2
      common/shft4/k3,jcal
      common/shft5/fiss
      common/shft2/fs(25),dsp(25),brr(25),der(25),er(25)
     1,cx(25),bilsum,xmiss,sumiz
      REAL*8 fiss
      common/fermst/aferm
      common/lab11/powaz(15,24,2200)
      common/levop/ldopt,igam
      common/ismrdat1/englev(1000,3),spingnd(106,260)
      common/ismrdat/isoid(106,260),levnm(1000),bevj(1000,3)
      common/ismrdat2/ididl(106,260),nolevs(106,260)
      common/isomer1/xgs(690),xiso(690),x2iso(690),isonum,isoint,iso2num
      common/isomer2/zs(690),as(690),eiso(690),eiso2(690)
      common/isomer3/iziso(690),iaiso(690),csgs(690),csis(690)
     1,csgse(690,100),csise(690,100),csis2e(690,100),rm1(690,100),
     2rm2(690,100)
      common/ijcpl/ajtarg,ajin
      common/lin/lspin
      common/spins2/spinu,spinf,prz1
      common/lscpl/bx(180),sppj
      common/tgspin/levno
      common/rotor/irotor,iidummy,arote(15,25,181)
       common/spins4/eprot,ealpha,eoff(8)
       common/spins3/ajmax
      common/sping/pfpej(600,180,100),pfpj(600,180)
      common/spins/rote(15,24,181),popej(15,24,180,100),fcx(180),ai,ecmp
      common/gate1/spind(5,180),popj(15,24,180),egatl(5),egatu(5),npart(
     15),igate
      common/nuopt/kinem,iem,kb
      common/distb/dsgp(36,999,8)
      common/hjk/jang,iq,rcss
      common/bound/neut2,iprot,idnuc(100),idm(100),einc(100),
     1thett(100),phit(100),ppp
      common/memo/scrs(999,8)
      common/memof/sfcrs(999,8),dfsgp(36,999,4),efn(4,999)
      common/pl3/jx,nepr,eb(100),rcsp(100),zee,amass,plex
      common/angprm/aprp(500),eav3(1000),
     1eav2(1000),pot2x,thet,phi,snt(360),cst(360),snp(720),csp(720)
     2,pot
      common/labcon/prx,pry,prz,pcm,rmas,clab(999,8),dslb(36,999,8)
      common/ilcom/nu,ndum,am(8)
      common/evapmc/jcount,ipremc,eminn(15,24)
      common/trans/itran
      common/lab10/pow(4,9999),gam(9999)
      REAL*8 pow,gam,powaz
        common/mccom/ncount,ipot,rum(15,24),ppmc(15,24,999),cf
      common/scat/mnz,mnja,probxx(2,2),q1(8,1100),fcs,eloss,fcspec
      common/exclbe/idz(4),ida(4),beff(20,20,2),pairmc(18,27)
      common/incr/ed,tened,efinal
      common/hyb2/pp(15,24,999)
      common/lym1/be(18,27,8),pair(18,27),xmas(18,27),delshl(18,27),amas
     1(18,27)
c     common/bym1/symb(18,27),symbp(18,27),symbs(18,27)
      common/sf6/aoz,zoz,en(8,999)
      common/lab3/sig(8,999),slav(8,999)
      common/te/pairx(8)
      common/mcevap/gm(8),gs(8,999)
      common/mcevpp/gn(8),gq(8,999)
      common/perp/spinprp
      common/mmtcom/mmt,mzee,iares,izres,light
      common/estar/estarcnt,e6,e7,e8
      common/lab12/pof(4,9999)
      REAL*8 pof
      common/frags2/rotef(260,181),arotef(260,181)
      common/fizk/csd,boff
      common/ifizk/jk,ik,nbr,nfr
      REAL*8 gm, gs,gg,ucla,g12,x
      REAL*8 gn,gq,g13
c we call for light,then heavy frag;each is for fiscrs from evap
      elost=0.
c above definition just for ff de-excite
 
c set jz,ja=1 for fragment starting de-excitation
       jz=zfis
       ja=amasf
       jn=ja-jz
c temp fix index 0 for jz
        jjj=jz*ja*jn
        if(jjj.lt.1)return
      ispinr=fjx
c define fission frag indices
      jrot=2*ispinr+1
      emax=ufis-elost 
       spinr=fjx
      if(emax.le.eminf(jz,jn))go to 1001
c here emax is reduced by rotational energy,so later fission barrier
c correction(rather than saddle point energy) is appropriate
** start loop on evap, initialize variables
 10       ucla = 0.
c define amasf,zfis as initial fragment A,Z-be sure to transfer as such
       iares=amasf
       izres=zfis
       ja=iares
        jz=izres
       if(iares.lt.aferm)then
       call gudima(emax,kst)
       call fermfrag(kst)
       light=1
       endif
       if(light.eq.1)return
      continue
          do 20 i = 1,8
             gn(i) = 0.
             pairx(i) = 0.
             do 20 j = 1,999
                gq(i,j) = 0.
 20       continue
** store pair corrections in pairx() array
          jn=ja-jz
          if(jz.lt.3)jz=3
          if(jn.lt.3)jn=3
          pairx(1) = pf(jz,jn-1)
          pairx(2) = pf(jz-1,jn)
          pairx(3) = pf(jz-2,jn-2)
** determine normalization parameter ucla and store
         gfis=0.
          do 90 k = 1,3
          gg=0.
             umax = emax - bq(jz,jn,k)
             if(umax.le.0.0) go to 81
             ib=(10.0*umax)  - pairx(k)-tened
             if(ib.lt.1)go to 90
             in = (umax+ed)/ed
             if(in.lt.1)go to 90
             if(in.gt.iem) in=iem
             do 80  ie = 1,in
           if(ib.lt.1)go to 80
                 gg= gg + pof(1,ib)*sif(jz-lzl(k),k,ie)
                gq(k,ie) = gg
                ib = ib-tened
 80          continue
 81          continue
             gn(k) = gg
 90       continue
          continue
       ucla=0.
       do 2345 l=1,3
         if(gn(l).gt.0.)then
          ucla=ucla+gn(l)
           go to 2345
            else
              gn(l)=0.
               endif
2345   continue           
       do 2347 l=1,3
        do 2346 k2=1,999
          if(gq(l,k2).gt.0.)go to 2346
            gq(l,k2)=0.
2346   continue
2347   continue
          ixc=(emax+ed)/ed
c defer half spin couple until ready to store ff
          fspec=csd/ed
** pick random number (0. to 1.)and normalize to total cross section
          if(ucla.le.0) go to 500
        g12=gn(1)+gn(2)
        g13=g12+gn(3)
       if(g13.le.0.)go to 1001
          x = rndm(-1.)*ucla
         jk=1
       if(x.ge.gn(1))jk=2
       if(x.ge.g12)jk=3
** it's a neutron,or proton,or alpha
          continue
           ik=1
         x=rndm(-1.)*gn(jk)
          do 100 ie = 1,iem
             if(gq(jk,ie) .ge. x) go to 101
             ik = ie+1
 100      continue
 101      continue
c  redefine channel energy,then define ecm
          ech=rndm(-1.)*ed+float(ik-1)*ed
       amasf=amasf-aof(jk)
       zfis=zfis-zof(jk)
       jz=zfis
       jn=amasf-zfis
       ja=amasf
         if(jn.lt.1)return
c here we will call the fizkinem subroutine to boost into lab energies
          nbr=2
         if(limout.lt.1)then
          boff=aof(jk)
          ecm=ech*amasf/(amasf+aof(jk))
          pcm=sqrt(2.*aof(jk)*(abs(ecm)))
          qcm=pcm
ct        pcmf=pcm
          ecmf=ecm
c------------------------------------------------
          call fizkinem(tote,af1,af2,ecmf,amasf)
c------------------------------------------------
          endif
            if(jz.gt.100.or.jn.gt.160)write(*,*)'jz,jn,k',jz,jn,k
          u = emax -ech-  bq(jz,jn,jk)
          emax = u
c
       dl=slavf(jz,jk,ik)
        if(dl.gt.50)write(*,*)'dl',dl
         branch=0.5
         if(rndm(-1.).le.branch)then
         spinr=abs(spinr-dl)
         else
          spinr=spinr+dl
         endif
          if(u .gt. eminf(jz,jn))go to 10 ! enough energy to evap more
c change spin index to 2j+1
 1001       continue
  500    continue
         ispinr=abs(spinr+0.5)
c round spn to nearest integer value as ispinr
c round to nearest half integer value if odd/even nolevs=1 nuclide 
          jres=(2*ispinr+1)
          spin=ispinr
      if(nolevs(jz,ja).eq.1)then
      x=rndm(-1.)
      if(x.lt.0.5)then
       spin=abs(spin-0.5)
       else
       spin=abs(spin+0.5)
       endif
       jres=(2.*spin+1.)
       endif
        if(jres.gt.180)jres=180
          if(jres.lt.1)jres=1
c 4/03 at this point,for storeable ff,spinr is integer or half
c integer spin in decimal form,ispinr is integer part of spin,
c and jres is an integer value of 2J+1 for the integer or half
c integer spin,where nolevs(z,a)=1 is oe/eo nucleus with half
c integer spin
** end loop on evap
           ixc=(u+ed)/ed
           if(ixc.gt.100)ixc=100
           if(ixc.lt.1)ixc=1
           n=isoid(jz,ja)
           if(n.gt.1)pfpej(n,jres,ixc)=pfpej(n,jres,ixc)+csd
           pfp(jz,jn,ixc)=pfp(jz,jn,ixc)+csd
           ppfi(jz,jn)=ppfi(jz,jn)+csd
           ppfaz(jz,ja)=ppfaz(jz,ja)+csd
           ppfz(jz)=ppfz(jz)+csd
           ppfa(ja)=ppfa(ja)+csd
            ppgam(ixc)=ppgam(ixc)+csd
       gam1=0.
       return
       end
c---sub for constants
      subroutine evapconst
      common/eder/ed2,v,vhalf,vsq
      common/incr/ed,tened,efinal
      common/nuopt/kinem,iem,kb
      common/evapcon/zof(8),aof(8)
      common/evapconi/lzl(8),laa(8),lnn(8)
      ed2=ed/2.
      tened=10.*ed
      iem=40./ed
      zof(1)=0.
       zof(2)=1.
        zof(3)=2.
         zof(4)=1.
          zof(5)=1.
           zof(6)=2.
            zof(7)=4.
      aof(1)=1.
       aof(2)=1.
        aof(3)=4.
         aof(4)=2.
          aof(5)=3.
           aof(6)=3.
            aof(7)=7.
      lzl(1)=0
      lzl(2)=1
      lzl(3)=2
      lzl(4)=1
      lzl(5)=1
      lzl(6)=2
      lzl(7)=4
      laa(1)=1
      laa(2)=1
      laa(3)=4
      laa(4)=2
      laa(5)=3
      laa(6)=3
      laa(7)=7
      lnn(1)=1
      lnn(2)=0
      lnn(3)=2
      lnn(4)=1
      lnn(5)=2
      lnn(6)=1
      lnn(7)=3
      return
      end
      subroutine evap
 
** Written by M.A.Ross
** MAR 970414
c this sr is to calculate actinide fission with fission width part of
c total and selected by random number
** This subroutine calculates evaporation of excited nuclei via Weiskopf-Ewing
** using Monte Carlo techniques.The differential cross sections for n,p,4He
** evaporation are calculated as a function of energy  and are stored and totaled
** (a running sum of the cross sections is also stored).
** A random number is chosen and a particle type for evaporation
** is found (n, p, or alpha).  The energy of the evaporated particle is chosen
** knowing the relative position of the random number wrt the calculated cross
** section for that cross section (i.e, by using a look up table). The process
** is repeated until the excitation energy of all residual nuclei is below
** the threshold for paricle emission.
**
      doubleprecision dafr1,dafr2,dzfr1,dzfr2,dufr1,dufr2,xfis1,demax
      doubleprecision ca,cz,darf,px,py,pz
      doubleprecision totke
      common/misscs/ppstor
      common/PRODA/PROD(7,20)
      common/clustr/m3set
      common/history/ide(60),idn(60),idith(60),iddum(60),the(60),adelj
     1(60),spinfh(60),enghis(60)
      common/out1/nout
      common/atsv/atsave
      common/histindx/ihist
      common/endvparm/iend,lll(6,6,6),ll(6,6)
      common/gam/io,la,ppgam(999),gmspec(999),beta,ams
      common/eder/ed2,v,vhalf,vsq
      common/tr/xtran
      common/evapcon/zof(8),aof(8)
      common/evapconi/lzl(8),laa(8),lnn(8)
      common/nuopt/kinem,iem,kb
      common/fisn/dfslb(36,999,4),cflab(999,4),dfetl(36,999,4)
      common/ifizk/jk,ik,nbr,nfr
      common/fizk/csd,boff
      common/tot/totke
      common/labcon3/pcmx1,pcmy1,pcmz1,pcmx2,pcmy2,pcmz2
      common/labcon4/v1x,v1y,v1z,v2x,v2y,v2z
      common/azfrag/dafr1,dafr2,dzfr1,dzfr2,dufr1,dufr2
      common/snglpri/ia1,ia2,iz1,iz2
      common/fis1c/xfis1
      common/barfis/darf
      common/fizinv/sif(100,4,999)
      common/fizstuf/pc(3)
      REAL*8 pc
      common/ldf/gamtem(9999)
      REAL*8 gamtem
      common/krld/powkr(9999),gamkr(9999),powgkr(9999)
      REAL*8 powkr,gamkr,powgkr
      common/snglpre/lufr1,lufr2,lafr1,lafr2,lzfr1,lzfr2
      common/shell/shelmol
      common/gr/rint(1100)
      REAL*8 rint
      common/gams/gamgam,csgam,csgamm(100)
      REAL*8 gam1,gamgam,gamfis,csgamm,gfis,rgam,fis1,csgam
       common/fizz/fizbar(106,260),fizcs(100)
      common/outlim/limout,lim2
      common/shft4/k3,jcal
      common/shft5/fiss
      common/shft2/fs(25),dsp(25),brr(25),der(25),er(25)
     1,cx(25),bilsum,xmiss,sumiz
      REAL*8 fiss
      common/lab11/powaz(15,24,2200)
      common/levop/ldopt,igam
      common/ismrdat2/ididl(106,260),nolevs(106,260)
      common/isomer1/xgs(690),xiso(690),x2iso(690),isonum,isoint,iso2num
      common/isomer2/zs(690),as(690),eiso(690),eiso2(690)
      common/isomer3/iziso(690),iaiso(690),csgs(690),csis(690)
     1,csgse(690,100),csise(690,100),csis2e(690,100),rm1(690,100),
     2rm2(690,100)
      common/ijcpl/ajtarg,ajin
      common/lin/lspin
      common/tagspin/trgspin
      common/lscpl/bx(180),sppj
      common/tgspin/levno
      common/rotor/irotor,iidummy,arote(15,25,181)
       common/spins4/eprot,ealpha,eoff(8)
      common/spins2/spinu,spinf,prz1
       common/spins3/ajmax
      common/sping/pfpej(600,180,100),pfpj(600,180)
      common/spins/rote(15,24,181),popej(15,24,180,100),fcx(180),ai,ecmp
      common/gate1/spind(5,180),popj(15,24,180),egatl(5),egatu(5),npart(
     15),igate
      common/distb/dsgp(36,999,8)
      common/hjk/jang,iq,rcss
      common/bound/neut2,iprot,idnuc(100),idm(100),einc(100),
     1thett(100),phit(100),ppp
      common/memo/scrs(999,8)
      common/pl3/jx,nepr,eb(100),rcsp(100),zee,amass,plex
      common/angprm/aprp(500),eav3(1000),
     1eav2(1000),pot2x,thet,phi,snt(360),cst(360),snp(720),csp(720)
     2,pot
      common/labcon/prx,pry,prz,pcm,rmas,clab(999,8),dslb(36,999,8)
      common/ilcom/nu,ndum,am(8)
      common/evapmc/jcount,ipremc,eminn(15,24)
      common/trans/itran
      common/lab10/pow(4,9999),gam(9999)
      REAL*8 pow,gam,powaz
        common/mccom/ncount,ipot,rum(15,24),ppmc(15,24,999),cf
      common/scat/jz,ja,probxx(2,2),q1(8,1100),fcs,eloss,fcspec
      common/exclbe/idz(4),ida(4),beff(20,20,2),pairmc(18,27)
      common/incr/ed,tened,efinal
      common/hyb2/pp(15,24,999)
      common/frags/xmasf(160,106),bq(106,160,4),pf(106,160)
     1,ppf(106,160,100),ppfi(106,160),ppfz(106),ppfa(260)
     2,ppfaz(106,260)
      common/sft5/exc(18,27),xmxx
      common/lym1/be(18,27,8),pair(18,27),xmas(18,27),delshl(18,27),amas
     1(18,27)
c     common/bym1/symb(18,27),symbp(18,27),symbs(18,27)
      common/sf6/aoz,zoz,en(8,999)
      common/lab3/sig(8,999),slav(8,999)
      common/te/pairx(8)
      common/mcevap/gm(8),gs(8,999)
      common/perp/spinprp
      common/shft/mc,mp,inver,ike,ipch,kq,qval,ap,at,zp,zt,cld,barfac
      common/mmtcom/mmt,mzee,iares,izres,light
      common/estar/estarcnt,e6,e7,e8
      REAL*8 gm, gs,gg,ucla,g12,x,denom
      REAL*8 g13,gamclust,gsum,crsclst
      REAL*8 denom1
      common/fermst/aferm
c>>>>>>>> the parameter m3cl will cause cluster emission
c only for the first equilibrated nuclide if =0, or for all
c evaporation cascade nuclei for >0 values. This could be made
c an input option, but first it should be tried as 0 and 1 to 
c see if it makes any difference in cluster emission, which should
c decrease rapidly with decreasing excitation
      m3cl=0
c if m3cl=1, clusters will be emitted throughout the evap cascade, else
c only for first equilibrated nucleus
c<<<<<<<<<<<<<<<<<<
      iem=50./ed
c>>>>>>>>> we may wish iem=50./ed to get full 7Be  contribution;
c >>>>>>>>> the value 40/ed truncates a little of the spectrum
c        write(*,*)'in evap'
      numx=0
       aferm=12.           
      light=0
      cs=fcs
      fiz=0.
      csgam=0.
      gfis=0.
      kmem=0
      jfis1=0
       jfis2=0
      m3num=m3set
c >>>>>>>>>>>> here set the number of particles emitted to 3 or 7;
c>>>>>>>>>> m3set is sent over from input routines 
      xmax=xtran
      emax = xmax - eloss       ! excitation energy of nucleus
c >>>>>>>>>>>>>>>> set spin composite for rotational E correction
           izz=zee+1.-float(jz)
           iaa=atsave+ap+2.-float(ja)-float(jz)
      spintem=spinu-spinf
      spinr=spintem
         ispinr=abs(spinr+0.5)
             spin=ispinr
          if(nolevs(izz,iaa).eq.1)then
          x=rndm(-1.)
          if(x.lt.0.5)then
           spin=abs(spin-0.5)
           else
           spin=abs(spin+0.5)
          endif
           endif
      jrot=2.*spin+1.
      if(jrot.gt.180)jrot=180
      if(jrot.lt.1)jrot=1
c correct spins so 2J+1 index is 1,3,5,7...for even A nuclides,
c and 2,4,6,8 for odd A nuclides
c >>>>>>>>>>>>>>>> set spin composite for rotational E correction
      corj=arote(jz,ja,jrot)
      emax=xmax-eloss -corj
c 3/16 move fermi breakup call to here
c mmt def needs to be adjusted for iso tags
       iares=mmt-ja-jz
       izres=mzee-jz
c>> do Fermi statistics breakup for light nuclei, set far A=7 by parm 'aferm'
      if(iares.lt.aferm)then
       call gudima(emax,kst)
       call fermfrag(kst)
       light=1
        return
       endif
c<<<<<<<<<<<<<<<<<< end Fermi statistics
c we do not allow fission for A<100
c 3/6/2016 move fiss barrier calc to here
c 3/7/15 define barf before redefinition
      barf=300.
       if(amass.lt.100)go to 1000
c above to bypass fission calc if A<100 since it will be negligible
         call seltzer(izres,iares,ispinr,barf,delr)
c      write(*,*)'l6231 return seltzer barf=',barf,'ispinr',ispinr
      darf=barf
       if(barf.le.0.)then
           barf=0.
           darf=0.
            endif
 1000   continue
c-----1/23-----------------------
c below,store element if can neither emit nor fiss
          if(emax.le.eminn(jz,ja).and.emax.le.barf)then
         inx=(emax+ed)/ed
      pp(jz,ja,inx)=pp(jz,ja,inx)+cs
        return
          endif
C
c 6/11 2017 now add case of emax.lt.eminn but gt. barf- the fission barrier
c           wherein the cs element will undergo fission,add test branch ifs
c           for case of this branch positive
          ifs=0
         if(emax.le.eminn(jz,ja).and.emax.gt.fisbar)then
          ifs=1
           csd=cs
          go to 97641
           endif
c <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< end rotational E
c here emax is reduced by rotational energy,so later fission barrier
c correction(rather than saddle point energy) is appropriate
c  july 2001 correct excitation for rotational energy
c >>>>>>>>>>>>>>>>>>> begin sequence on particle emission/cluster/fission 
 10       ucla = 0.
          denom=0.
       if(ja.le.neut2.and.jz.le.iprot)go to 9005
       ppp=ppp+fcs
       if(ja.gt.22.or.jz.gt.13)return
c      return
c >>>>>>>>> if array will be out of bounds Z,A, store cs element & return
9005  continue
c >>>>>>>>initialize particle emission relative widths
          do 20 i = 1,8
             gm(i) = 0.
             pairx(i) = 0.
             do 20 j = 1,999
             gs(i,j) = 0.
 20       continue
          pairgam=0.
          pairfis=0.
        if(jz.gt.11)jz=11
        if(ja.gt.21)ja=21
c <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
** store pair corrections in pairx() array
         call pairset(pairfis,pairgam)
       if(ldopt.ge.1.and.emax.le.210.)go to 9010
c >>>>>>>>>>>> branch to 9010 for other than Fermi gas densities
** determine normalization parameter ucla and store
c >>>>>>>>>>>> calc partial widths Fermi gas here, non-FG at 9010
          do 90 k = 1,m3num
          gd=0.
                 i=k
                 if(i.gt.4)i=3
             umax = emax - be(jz,ja,k)
             ib=(10.0*umax)  - pairx(k)-tened
             if(ib.lt.1)go to 90
             in =(umax+ed)/ed
             if(in.lt.1)go to 90
             if(in.gt.iem) in=iem
             do 80  ie = 1,in
             if(ib.lt.1)go to 81
                 gd= gd + pow(i,ib)*sig(k,ie)
                 gs(k,ie) = gd
                 ib = ib-tened
 80          continue
 81          continue
                 gm(k) = gd
 90          continue
c----------------loop on 90 does n,p,4he widths for FG only----

             go to 9020
9010         continue
c <<<<<<<<<<< end FG widths <<<<<<<<<<<<<<<<<<<<<<
c >>> do partial widths for KR or Obninsk level densities 
          do 901 k=1,m3num
c neut id
          gg=0.
          jzz=jz+lzl(k)
          jaa=ja+lnn(k)
          if(jaa.gt.23)jaa=23
          if(jzz.gt.13)jzz=13
            umax = emax - be(jz,ja,k)
             if(umax.le.0.0) go to 901
              ib=(10.0*umax)  - tened
         if(ib.lt.1)go to 901
c here for KR or Obn level densities, we do not add a pairing correction
             in = umax/ed
             if(in.lt.1)go to 901
             if(in.gt.iem) in=iem
             do 801 ie = 1,in
           if(ib.lt.1)go to 811
                 gg= gg + powaz(jzz,jaa,ib)*sig(k,ie)
                gs(k,ie) = gg
                ib = ib-tened
 801         continue
 811         continue
             gm(k) = gg
 901      continue
c>>>>>>>>>>>> end Obninsk/KR level density partial widths 
9020      continue
c<<<<<<<<<<<<< add call to fission management/width routine
c we do not allow fission for A<100
       gfis=0.
       fis1=0.
         umx=emax-barf
          ib=10.0*(umx)-pairfis
          if(umx.le.0..or.ib.lt.1.or.barf.gt.35.)then
      gfis=0.
      else
c      if(ldopt.le.1) gfis=gam(ib)
c june 2017 do not allow fission except for FG level density
       if(ldopt.eq.0)gfis=gam(ib)
c.........................................................
c      write(*,*)'ldopt,gfis,gam(ib)',ldopt,gfis,gam(ib)
       endif
c     write(*,*)' jz=',jz,'ja=',ja,'emax=',emax,'umx=',umx,'J',ispinr
c      write(*,99)barf,gfis,gm(1),gm(2)
   99 format('in evap barf,gfis,gm(1),gm(2)',4e10.2)
         if(gfis.eq.0.)fis1=0.
c >>>>>>>>> end fission width calculation>>>>>>>>>>>>>>>>>>
        do 9021 jq=1,m3num
         if(gm(jq).gt.0.)go to 9021
          gm(jq)=0.
 9021 continue
c <<<<<<  loops to zero gs elements which have not been initialized
        do 9022 jq=1,m3num
         do 9023 iq=1,999
          if(gs(jq,iq).gt.0.)go to 9023
           gs(jq,iq)=0.
 9023 continue
 9022 continue
c <<<<<<<<<<<<<<<<<<<<<<<<< end loops to remove NaN's
      ucla=0.
      denom=0.
          do jq=1,3
           ucla=ucla+gm(jq)
            enddo
c >>>>> sum partial widths for n,p,4He emission
       if(ucla.gt.0.)then
        denom=ucla
         else
          denom=0.
           endif
         rgam=0. 
        g12=gm(1)+gm(2)
        g13=g12+gm(3)
c below g1x sums for pure MC calc; we are doing hybrid version
c wherein low probability events are treated every event with
c diminished cross section.
c       g14=g13+gm(4)
c       g15=g14+gm(5)
c       g16=g15+gm(6)
c       g17=g16+gm(7)
c       g18=g17+rgam
c       g15=g14+gfis
c >>>>>>>>>> above calcs particle emission sums on partial widths
c calculate partial fission width,then fission cross section fis1
c also calculate photon width and cross section
          ixc=(emax+ed-pairgam/10.)/ed
          imc=(emax+ed)/ed
       if(ixc.gt.999)ixc=999
       if(imc.gt.999)imc=999
          if(ixc.lt.1)ixc=1
          if(imc.lt.1)imc=1
          gam1=0.
          fis1=0.
          rgam=0.
c         denom=0.
            if(jz.eq.1.and.ja.eq.1)rgam=rint(ixc)
c calculate partial gamma width >>>>>>>>>>>>>>>>
c for first emitting nucleus
c         denom=ucla
          denom=denom+rgam+gfis
          denom1=denom
          if(denom.le.0.)go to 500
          if(m3num.le.3)go to 9999
c calculate cluster width fraction
c calculate cluster cross sections
          gamclust=0.
          crsclst=0.
          do kj=4,7
           gamclust=gamclust+gm(kj)
            enddo
            denom1=denom
           if(gamclust.gt.0.)then
             denom1=denom+gamclust
                              endif
c now calculate the fraction of element of cs this event going into cluster
c decay, which is a non-MC sampling to enhance events of rare events by
c using a fraction of event cs every pass of event loop, not just when MC
c selected.
           if(denom1.gt.0.)then
              crsclst=(gamclust*cs/denom1)
               else
              crsclst=0.
                endif
            cs=cs-crsclst
c here decrease event cs since a small piece will be used as cluster emission
c every event, rather than the entire element cs when randomly selected.
            if(cs.lt.0.)cs=0.
            if(cs.le.0.)go to 500
            if(crsclst.lt.0.)crsclst=0.
            if(crsclst.le.0.)go to 9999
         emaxx=emax
          jzsave=jz
          jasave=ja
          if(m3num.gt.3)call clustevap(crsclst,emaxx)
          ja=jasave
          jz=jzsave
c
 9999    continue
         if(denom1.gt.0.)then
           gamgam=rgam/denom1
           gamfis=gfis/denom1
                         else
             gamgam=0.
             gamfis=0.
              rgam=0.
              gfis=0.
                         endif
          if(cs.ge.0.)then
            gam1=gamgam*cs
             else
              gam1=0.
                endif
           if(jz.eq.1.and.ja.eq.1.and.denom.gt.0.)then
c store radiative emission cs for first nuclide
c in evaporation chain only if it is composite nucleus
           pp(jz,ja,imc)=pp(jz,ja,imc)+gam1
           itran=imc
          izz=zee+1.-float(jz)
          iaa=amass+2.-float(ja)-float(jz)
c >>>>>>>>>>>>> calculate spin of radiative capture nucleus- is it worth the time?
           ispinr=abs(spinr+0.5)
           spin=ispinr
          if(nolevs(izz,iaa).eq.1)then
          x=rndm(-1.)
          if(x.lt.0.5)then
           spin=abs(spin-0.5)
           else
           spin=abs(spin+0.5)
          endif
          endif
c change spin index to 2j+1
           jres=((2.*spin)+1.)
           if(jres.gt.180)jres=180
           if(jres.lt.1)jres=1
           if(ixc.gt.100)ixc=100
           if(imc.gt.100)imc=100
          if(jz.lt.14.and.ja.lt.25)then
          popej(jz,ja,jres,imc)=popej(jz,ja,jres,imc)+gam1
          popj(jz,ja,jres)=popj(jz,ja,jres)+gam1
          endif
          endif
c <<<<<<<<<<<<< end store radiative capture element residual nucleus
          if(cs.le.0.)go to 500
c
c >>>>>>>>>calculate fission cross section for non-actinides
          fis1=gamfis*cs
c >>>> for non-actinide fission, fis1 cross section will be emitted
c every event, whereas for actinides, cs will go into fission when
c a MC choice determines that this event will be a fission event
c below we reduce cs element for A<220 fission and for cluster events, which
c are apportioned and treated every event as fraction width/cross section
c rather than MC selected for whole element of cs.
c <<<<<<<<<<<<<<<<<<<<<<<
c
        if(gam1.gt.0.)cs=cs-gam1
             if((atsave+ap).lt.220.)then
               if(fis1.gt.0.)cs=cs-fis1
                if(cs.ge.0.)go to 19762
                  cs=0.
19762   continue
          csd=0.
          if(fis1.gt.0.)csd=fis1 
                go to 9764
                 endif
             if((atsave+ap).ge.220.)csd=cs
c set csd,the element for fission,to the whole value
c for actinide fission it is only the fission fraction in evap
c for A lt.220,i.e. non monte carlo method for very small cs
9762   x=rndm(-1.)*denom1
          if(x.le.g13.and.(atsave+at).ge.220.)go to 9765
c  it is a prticle emission if x.le.g13, the particle emission width
9764   continue
c---------------------------------
97641  continue
         fizcs(nepr)=csd+fizcs(nepr)
          fspec=cs/ed
c---------------------------------
c make running sum fiss of fission cross section over all events
          fiz=fiz+csd
c call subroutine to store fission event on A,Z,excitation (J?)
        a=atsave+ap+2.-float(ja)-float(jz)
        z=zee+1.-float(jz)
        pc(1)=0.
        pc(2)=0.
        pc(3)=0.
        ca=a
        cz=z
        xfis1=csd
        fis1=csd
        demax=emax
        px=0.
        py=0.
        pz=0.
cxxxxxxxxxxxxxxxxxxxxxx
c       go to 9765
cxxxxxxxxxxxxxxxxxxx
        if(csd.gt.0.)then
        call tfis(ca,cz,demax,px,py,pz)
cxxxxxxxxxxxxxxxxxxxxxx
c       go to 9765
cxxxxxxxxxxxxxxxxxxx
       if(dafr1.eq.0..or.dafr2.eq.0.)write(23,*)'zero evap'
       if(dafr1.eq.0..or.dafr2.eq.0.)go to 9765
        uf1=dufr1
        af1=dafr1
        zf1=dzfr1
          uf2=dufr2
          af2=dafr2
          zf2=dzfr2
          tote=totke
        call fizfr(af1,af2,uf1,uf2,spinr,jfis1,jfis2,zf1,zf2)
c call fizfr to apportion spinr  cn spin between frags,couple it
c with fission frag spin algorithm to give primary frag spins to
c be modified in evap1 due to de-excitation of primary to final frags
c frag spins are returned as jfis1 and jfis2
c spins are returned as jfis1,jfis2;float these as fj1 and fj2
c call to de-excite ff1 to particle bound
         fj1=float(jfis1)
         fj2=float(jfis2)
         nbr=1
         nfr=0
         boff=1.
        call fizkinem(tote,af1,af2,ecmf,amasf)
c>>>>>>>>>>>>>> now we de-excite the 2 fission frags by n,p,4He evaporation
        nfr=1
        call evap1(af1,zf1,uf1,fj1,tote)
c call to de-excite ff2 to particle bound
        nfr=2
        call evap1(af2,zf2,uf2,fj2,tote)
c done with fission fragments <<<<<<<<<<<<<<<<
       endif
      if((atsave+ap).lt.220.)go to 9765
c  >>>>>>>>> for non-actinides, we continue evaporation cascade;
c for actinides, if this nuclide fizzed, we return for next event
        return
 9765   continue 
        if(ifs.eq.1)return
c    above added 6/12/2017;if event is fission only,return(case E*<eminn and>Bf)
      fspec=cs/ed
          if(cs.eq.0.)go to 500
**
** pick random number (0. to 1.)and normalize to total cross section
          if(ucla.le.0) go to 500
c kk is particle type,n,p,alpha=1,2,3
          x=rndm(-1.)*ucla
          kk=1
          if(x.ge.gm(1))kk=2
          if(x.gt.g12)kk=3
          rmas=atsave+ap+2.-ja-jz-laa(kk)
         x=rndm(-1.)*gm(kk)
          ik = 1
          do 100 ie = 1,iem
             if(gs(kk,ie) .ge. x) go to 101
             ik = ie+1
 100      continue
 101      continue
c  redefine channel energy,then define ecm
          ech=rndm(-1.)*ed+float(ik-1)*ed
          u = emax -ech-  be(jz,ja,kk)
          iu = (u+ed)/ed 
         if(iu.lt.1)iu=1
          ecm=ech*rmas/(rmas+laa(kk))
          pcm=sqrt(2.*float(laa(kk))*ecm)
          qcm=pcm
c 9/20/2014 change index ik to represent cm rather than channel energy
           ik=(ecm+ed)/ed
          x=rndm(-1.)
c
          thet=acos(1.-2.*x)
          phi=rndm(-1.)*6.28318
c
          nu=kk
c------------------------------------------
          call recoil(1)
c------------------------------------------
c  >>>>>>>>>>>>>>>>>>> store the cm emission spectra from emitting nucleus
       ith=thet*11.459155+1.
       if(ith.eq.37)ith=36
       scrs(ik,kk)=scrs(ik,kk)+fspec
       dsgp(ith,ik,kk)=dsgp(ith,ik,kk)+fspec
c <<<<<<<<<<<<<< store ddcs in cm <<<<<<<<<<<<<<<<<<<<<<<<
c try to make delta j dependent on level density
       if(iend.eq.1)call endf(thet,ech,ith,ik,kk )
          ja = ja + lnn(kk)
          jz=jz+lzl(kk)
          emax=u
         iares=amass+2.-float(ja)-float(jz)
           if(iares.lt.aferm)then
              izres=zee+1.-float(jz)
                call gudima(u,kst)
                 call fermfrag(kst)
                  light=1
                  return
                   endif
c >>>>>>> reset nucleus excitation after emission, see if another 
c particle may be emitted - should add test to see if it can fiss
c but not emit particles, then make it fiss, else go to 10 to restart
c process
 
c  now calculate nucleus angular momentum >>>>>>>>>>>>>>>>
        kkg=kk
          if(kkg.gt.4)kkg=3
c--------------------------------
       dl=slav(kk,ik)
       dlo=abs(spinr-dl)
       ldlo=2.*(dlo)+1.
       dlhi=abs(spinr+dl)
       ldhi=2.*(dlhi)+1.
        if(ldhi.gt.180)ldhi=180
        if(ldlo.gt.180)ldlo=180
        if(ldhi.lt.1)ldhi=1
        if(ldlo.lt.1)ldlo=1
       uup=abs(u-rote(jz,ja,ldhi))
       ulo=abs(u-rote(jz,ja,ldlo))
          iup =uup*10.+1.
          iul =ulo*10.+1.
         branch=pow(kkg,iul)/(pow(kkg,iul)+pow(kkg,iup))
         if(rndm(-1.).le.branch)then
         spinr=abs(spinr-dl)
         else
          spinr=spinr+dl
         endif
           jspinr=spinr
       if(jz.gt.13)jz=13
       if(ja.gt.24)ja=24
      m3num=3
      if(m3cl.gt.0)m3num=m3set
c reset m3num if cluster emission desired throughout evap cascade
c-------------------------------------------------------------------
          if(u .gt. eminn(jz,ja))go to 10 ! enough energy to evap more
c-------------------------------------------------------------------
c since we are in actinide region, we must check to see that this
c residual nucleus not only has not enough energy to evap more, but
c also that it is below the fission barrier
      if(amass.ge.220.)then
       iares=mmt-ja-jz
        izres=mzee-jz
         call seltzer(izres,iares,jspinr,barf,delr)
          if(u.le.eminn(jz,ja).and.u.ge.barf)then
           csd=cs
            go to 9764
             endif
              endif
c the transfer to 9764 will put this cross section into fission
         continue
c >>>>>>>> there was not sufficient energy to evap more, (or to fiss)so 
c  store cs into residual nucleus >>>>>>>>>>>>>>
c--------------------------------------------------------
          iu=(u+ed)/ed
            if(iu.lt.1)iu=1
          pp(jz,ja,iu) = pp(jz,ja,iu) + cs
c------------------------------------------------------
         ppgam(iu)=ppgam(iu)+cs
c aug 07 put in array to be used for equilibrium gamma cascade
         ispinr=abs(spinr+0.5)
          spinr=ispinr
c here must be sure indexed n,p,alpha definition gives right izz,iaa
      izz=zee+1.-float(jz)
      iaa=atsave+ap+2.-float(ja)-float(jz)
c get correct integer or non- integer spin
      if(nolevs(izz,iaa).eq.1)then
      x=rndm(-1.)
      if(x.lt.0.5)then
       spin=abs(spinr-0.5)
       else
       spin=abs(spinr+0.5)
       endif
       endif
c change spin index to 2j+1
          jres=(2.*spin)+1.
          if(jres.gt.180)jres=180
          if(jres.lt.1)jres=1
          if(iu.le.100)then
          popej(jz,ja,jres,iu)=popej(jz,ja,jres,iu)+cs
          popj(jz,ja,jres)=popj(jz,ja,jres)+cs
           endif
        continue
** end loop on evap
 500   continue
       continue
c23456789012345678901234567890123456789012345678901234567890123456789012
       gam1=0.
cxxxxxxxxxxxxxxxxx
       return
       end
      subroutine clustevap(crsclst,emaxx)
 
      doubleprecision dafr1,dafr2,dzfr1,dzfr2,dufr1,dufr2,xfis1,demax
      doubleprecision ca,cz,darf,px,py,pz
      doubleprecision totke
      common/PRODA/PROD(7,20)
      common/clustr/m3set
      common/atsv/atsave
      common/shft/mc,mp,inver,ike,ipch,kq,qval,ap,at,zp,zt,cld,barfac
      common/endvparm/iend,lll(6,6,6),ll(6,6)
      common/history/ide(60),idn(60),idith(60),iddum(60),the(60),adelj
     1(60),spinfh(60),enghis(60)
      common/out1/nout
      common/histindx/ihist
      common/eder/ed2,v,vhalf,vsq
      common/tr/xtran
      common/evapcon/zof(8),aof(8)
      common/evapconi/lzl(8),laa(8),lnn(8)
      common/nuopt/kinem,iem,kb
      common/fisn/dfslb(36,999,4),cflab(999,4),dfetl(36,999,4)
c     common/ifizk/jk,ik,nbr,nfr
c     common/fizk/csd,boff
      common/tot/totke
      common/labcon3/pcmx1,pcmy1,pcmz1,pcmx2,pcmy2,pcmz2
      common/labcon4/v1x,v1y,v1z,v2x,v2y,v2z
      common/azfrag/dafr1,dafr2,dzfr1,dzfr2,dufr1,dufr2
      common/snglpri/ia1,ia2,iz1,iz2
      common/fis1c/xfis1
      common/barfis/darf
      common/fizinv/sif(100,4,999)
      common/fizstuf/pc(3)
      REAL*8 pc
      common/ldf/gamtem(9999)
      REAL*8 gamtem
      common/krld/powkr(9999),gamkr(9999),powgkr(9999)
      REAL*8 powkr,gamkr,powgkr
      common/snglpre/lufr1,lufr2,lafr1,lafr2,lzfr1,lzfr2
      common/shell/shelmol
      common/gr/rint(1100)
      REAL*8 rint
      common/gams/gamgam,csgam,csgamm(100)
      REAL*8 gam1,gamgam,gamfis,csgamm,gfis,rgam,fis1,csgam
       common/fizz/fizbar(106,260),fizcs(100)
      common/outlim/limout,lim2
      common/shft4/k3,jcal
      common/shft5/fiss
      common/shft2/fs(25),dsp(25),brr(25),der(25),er(25)
     1,cx(25),bilsum,xmiss,sumiz
      REAL*8 fiss
      common/fermst/aferm
      common/lab11/powaz(15,24,2200)
      common/levop/ldopt,igam
      common/ismrdat2/ididl(106,260),nolevs(106,260)
      common/isomer1/xgs(690),xiso(690),x2iso(690),isonum,isoint,iso2num
      common/isomer2/zs(690),as(690),eiso(690),eiso2(690)
      common/isomer3/iziso(690),iaiso(690),csgs(690),csis(690)
     1,csgse(690,100),csise(690,100),csis2e(690,100),rm1(690,100),
     2rm2(690,100)
      common/lin/lspin
      common/ijcpl/ajtarg,ajin
      common/lscpl/bx(180),sppj
      common/tagspin/trgspin
      common/tgspin/levno
      common/spins2/spinu,spinf,prz1
      common/rotor/irotor,iidummy,arote(15,25,181)
       common/spins4/eprot,ealpha,eoff(8)
       common/spins3/ajmax
      common/sping/pfpej(600,180,100),pfpj(600,180)
      common/spins/rote(15,24,181),popej(15,24,180,100),fcx(180),ai,ecmp
      common/gate1/spind(5,180),popj(15,24,180),egatl(5),egatu(5),npart(
     15),igate
      common/distb/dsgp(36,999,8)
      common/hjk/jang,iq,rcss
      common/bound/neut2,iprot,idnuc(100),idm(100),einc(100),
     1thett(100),phit(100),ppp
      common/memo/scrs(999,8)
      common/pl3/jx,nepr,eb(100),rcsp(100),zee,amass,plex
      common/angprm/aprp(500),eav3(1000),
     1eav2(1000),pot2x,thet,phi,snt(360),cst(360),snp(720),csp(720)
     2,pot
      common/labcon/prx,pry,prz,pcm,rmas,clab(999,8),dslb(36,999,8)
      common/ilcom/nu,ndum,am(8)
      common/evapmc/jcount,ipremc,eminn(15,24)
      common/trans/itran
      common/lab10/pow(4,9999),gam(9999)
      REAL*8 pow,gam,powaz
        common/mccom/ncount,ipot,rum(15,24),ppmc(15,24,999),cf
      common/scat/jz,ja,probxx(2,2),q1(8,1100),fcs,eloss,fcspec
      common/exclbe/idz(4),ida(4),beff(20,20,2),pairmc(18,27)
      common/incr/ed,tened,efinal
      common/hyb2/pp(15,24,999)
      common/frags/xmasf(160,106),bq(106,160,4),pf(106,160)
     1,ppf(106,160,100),ppfi(106,160),ppfz(106),ppfa(260)
     2,ppfaz(106,260)
 
      common/sft5/exc(18,27),xmxx
      common/lym1/be(18,27,8),pair(18,27),xmas(18,27),delshl(18,27),amas
     1(18,27)
c     common/bym1/symb(18,27),symbp(18,27),symbs(18,27)
      common/sf6/aoz,zoz,en(8,999)
      common/lab3/sig(8,999),slav(8,999)
      common/te/pairx(8)
      common/mcevap/qm(8),qs(8,999)
      common/perp/spinprp
      common/mmtcom/mmt,mzee,iares,izres,light
      common/estar/estarcnt,e6,e7,e8
      dimension gm(8),gs(8,999)
      REAL*8 gm, gs
      REAL*8 qm, qs,gg,ucla,g12,x,denom,cscl,gmcl,gscl
      REAL*8 g13,gamclust,gsum,crsclst,cq,fspec
c transfer cluster partial widths to this  evap routine
c as gm(ii),ii=4-7, and gs(ii,e), the running sum for each cluster
c set the particle widths to new arrays so that changed values will
c not be returned to evap or evapt
c
c>>>>>>>>>> set particle emission widths into new arrays so as not to
c change values in calling routine evap/evapt
      jres=0
      do 8003 ii=4,7
       gm(ii)=qm(ii)
       do 8002 jj=1,999
        gs(ii,jj)=qs(ii,jj)
 8002 continue
 8003 continue
c <<<<<<<<<<<<<<<< end save arrays
      jasave=ja
      jzsave=jz
c >>>>>>>>>>>>>>> save jz,ja to return original values to calling routine
      iex=60./ed
       gsum=0.
      do ii=4,7
       if(gm(ii).gt.0.)then
      gsum=gsum+gm(ii)
      endif
      enddo
       if(gsum.le.0.)return
c here start loop to evaporate clusters 4-7, for first nuclide
c reaching the evap routine; after the cluster, each nuclide
c will de-excite to ground by n,p,4He emission. No fission
c in this loop- it is handled in calling sub.
      do 10000 ii=4,7
       xmax=emaxx
      if(gm(ii).eq.0.)go to 10000
       if(emaxx.le.0.)go to 10000
c reset excitation and ja,jz for each cluster emission to value sent
c over from the precompound routine to the evaporation routine
       ja=jasave
       jz=jzsave
c loop over each cluster ii, with emission cross section cq
       cq=0.
      if(crsclst.gt.0.)cq=(gm(ii)/gsum)*crsclst
      fspec=cq/ed
         rmas=atsave+ap+2.-ja-jz-laa(ii)
         x=rndm(-1.)*gm(ii)
          ut = (xmax- be(jz,ja,ii))
           iem=(ut+ed)/ed
           if(iem.gt.iex)iem=iex
            if(iem.lt.1)go to 10000
          ik=1
          do 900 ie = 1,iem
             if(gs(ii,ie) .ge. x) go to 981
             ik = ie+1
 900      continue
 981      continue
c  redefine channel energy,then define ecm
          ech=rndm(-1.)*ed+float(ik-1)*ed
          u = xmax -ech-  be(jz,ja,ii)
          iu = (u+ed)/ed 
         if(iu.lt.1)iu=1
          ecm=ech*rmas/(rmas+laa(ii))
c 9/20/2014 change ik index to represent cm rather than channel energy
            ik=(ecm+ed)/ed
c         pcm=sqrt(2.*float(laa(ii))*(abs(ecm)))
c         qcm=pcm
          x=rndm(-1.)
 
c         thet=acos(1.-2.*x)
          phi=rndm(-1.)*6.28318
c
c         nu=ii
c here check that recoil sr does not change values for main
c evaporation sequence
c         call recoil
       ith=thet*11.459155+1.
       if(ith.eq.37)ith=36
c-------------------------------------------------
       scrs(ik,ii)=scrs(ik,ii)+fspec
       dsgp(ith,ik,ii)=dsgp(ith,ik,ii)+fspec
c-------------------------------------------------
c try to make delta j dependent on level density
c      endif
c here we have stored sdcs and ddcs for cluster, next must
c de-excite nucleus to ground
       if(iend.eq.1)call endf(thet,ech,ith,ik,nu)
          ja = ja + lnn(ii)
          jz=jz+lzl(ii)
          emax = u
          if(emax.le.0.)go to 10000
c set ja,jz, excitation emax of residual after cluster
c        if(u .le. eminn(jz,ja))write(*,*)'l6772clust return'
c----------------------------------------------
         if(u .gt. eminn(jz,ja))go to 10 ! enough energy to evap more
c move transfer to 10 until after delta 'l' calc -dec 03
           izz=zee+1.-jz
           iaa=amass+2.-ja-jz
          ucla = 0.
       iares=mmt-ja-jz
       izres=mzee-jz
       ztem=izres
       if(ja.gt.21.or.jz.gt.11)then
       ja=jasave
       jz=jzsave
       return
       endif
       if(ja.le.neut2.and.jz.le.iprot)go to 9005
c here the element cs is the cross section for each cluster
c to be emitted, and therefore for each particle in the cascade
c if we are about to exceed the Z,A array dimension, we keep
c track of the cross section (ppp) and exit
       ppp=ppp+cq
       ja=jasave
       jz=jzsave
       return
9005  continue
c here have done cluster decay, next at 10 start n,p,4He decay
c cascade until no more particles can be emitted
   10  continue
          do 20 i = 1,3
             gm(i) = 0.
             pairx(i) = 0.
             do 20 j = 1,999
                gs(i,j) = 0.
 20       continue
 
        if(jz.gt.11)jz=11
        if(ja.gt.21)ja=21
** store pair corrections in pairx() array
         call pairset(pairfis,pairgam)
       if(ldopt.ge.1.and.emax.le.210.)go to 9010
** determine normalization parameter ucla and store
          do 90 k = 1,3
c first pass m3num=7; may make it 3 for further passes
          gg=0.
             umax = emax - be(jz,ja,k)
             if(umax.le.0.)go to 90
c here calculate contribution to width of each particle type 
             ib=(10.0*umax)  - pairx(k)-tened
             if(ib.lt.1)go to 90
             in = (umax+ed)/ed
             if(in.lt.1)go to 90
             if(in.gt.iex) in=iex
             do 80  ie = 1,in
           if(ib.lt.1)go to 81
                 gg= gg + pow(k,ib)*sig(k,ie)
                gs(k,ie) = gg
                ib = ib-tened
 80          continue
 81          continue
             gm(k) = gg
 90       continue
          go to 9020
9010      continue
 
          do 901 k=1,3
c neut id
          gg=0.
          jzz=jz+lzl(k)
          jaa=ja+lnn(k)
          if(jaa.gt.23)jaa=23
          if(jzz.gt.13)jzz=13
             umax = emax - be(jz,ja,k)
           if(umax.lt.0.0) go to 811
             ib=(10.0*umax) -tened
             if(ib.lt.1)go to 901
             in = (umax-pairx(k)-tened)/ed
             if(in.lt.1)go to 901
             if(in.gt.iex)in=iex
             do 801 ie = 1,in
           if(ib.lt.1)go to 811
                 gg= gg + powaz(jzz,jaa,ib)*sig(k,ie)
                gs(k,ie) = gg
                ib = ib-tened
 801         continue
 811         continue
             gm(k) = gg
 901      continue
 
9020      continue
          ucla=0.
          do 9021 jk=1,3
          if(gm(jk).gt.0.)go to 9021
           gm(jk)=0.
 9021      continue
           do 9022 jk=1,3
          do 9022 jq=1,999
          if(gs(jk,jq).gt.0.)go to 9022
          gs(jk,jq)=0.
9022      continue
          ucla=gm(1)+gm(2)+gm(3)
          if(ucla.eq.0.)go to 9000
        g12=gm(1)+gm(2)
        g13=g12+gm(3)
           ispinr=abs(spinr+0.5)
           spin=ispinr
c 8/5/04
          izz=zee+1.-float(jz)
            iaa=atsave+ap+2.-float(ja)-float(jz)
          if(nolevs(izz,iaa).eq.1)then
          x=rndm(-1.)
          if(x.lt.0.5)then
           spin=abs(spin-0.5)
           else
           spin=abs(spin+0.5)
          endif
          endif
c  change spin index to 2j+1
           jres=(2.*spin)+1.
           if(jres.gt.180)jres=180
          if(cs.lt.0.)cs=0.
        a=atsave+ap+2.-float(ja)-float(jz)
        z=zee+1.-float(jz)
        ztem=z
** pick random number (0. to 1.)and normalize to total cross section
          if(ucla.gt.0.)then
        x=rndm(-1.)*ucla
            else
             ucla=0.
              endif
c kk is particle type,n,p,alpha=1,2,3
         if(ucla.eq.0.)go to 500
          kk=1
          if(x.ge.gm(1))kk=2
          if(x.gt.g12)kk=3
          rmas=atsave+ap+2.-ja-jz-laa(kk)
         x=rndm(-1.)*gm(kk)
          ik = 1 
          do 100 ie = 1,iem
             ik=ie
             if(gs(kk,ie) .ge. x) go to 101
 100      continue
 101      continue
c  redefine channel energy,then define ecm
          ech=rndm(-1.)*ed+(float(ik)-1)*ed
          u = emax -ech-  be(jz,ja,kk)
          iu = (u+ed)/ed 
         if(iu.lt.1)iu=1
          ecm=ech*rmas/(rmas+laa(kk))
c 9/20/2014 change ik reference from channel to cm energy
            ik=(ecm+ed)/ed
          x=rndm(-1.)
c
          thet=acos(1.-2.*x)
          phi=rndm(-1.)*6.28318
c
          nu=kk
c         call recoil
       ith=thet*11.459155+1.
       if(ith.eq.37)ith=36
c---------------------------------------------
       scrs(ik,kk)=scrs(ik,kk)+fspec
       dsgp(ith,ik,kk)=dsgp(ith,ik,kk)+fspec
c---------------------------------------------
       if(iend.eq.1)call endf(thet,ech,ith,ik,kk)
          ja = ja + lnn(kk)
          jz=jz+lzl(kk)
          emax = u
          if(u .gt. eminn(jz,ja))go to 10 ! enough energy to evap more
 9000    continue
c------------------------------------------------
          pp(jz,ja,iu) = pp(jz,ja,iu) + cq
c-----------------------------------------------
         ispinr=abs(spinr+0.5)
          spinr=ispinr
c here must be sure indexed n,p,alpha definition gives right izz,iaa
      izz=zee+1.-jz
      iaa=atsave+ap+2.-ja-jz
      if(nolevs(izz,iaa).eq.1)then
      x=rndm(-1.)
      if(x.lt.0.5)then
       spin=abs(spinr-0.5)
       else
       spin=abs(spinr+0.5)
       endif
       endif
c change spin index to 2j+1
          jres=(2.*spin)+1.
          if(jres.gt.180)jres=180
          if(jres.lt.1)jres=1
      if(iu.le.100)then
          popej(jz,ja,jres,iu)=popej(jz,ja,jres,iu)+cq
          popj(jz,ja,jres)=popj(jz,ja,jres)+cq
           endif
        continue
** end loop on evap
 500   continue
cxxxxxxxxxxxxxxxxx
10000 continue
      ja=jasave
      jz=jzsave
c end loop on cluster type
       return
       end
c
      subroutine pairset(pairfis,pairgam)
      common/te/pairx(8)
      common/lym1/be(18,27,8),pair(18,27),xmas(18,27),delshl(18,27),amas
     1(18,27)
      common/scat/jz,ja, probxx(2,2),q1(8,1100),fcs,eloss,fcspec
      common/sft5/exc(18,27),xmxx
          pairx(1) = pair(jz,ja+1)
          pairx(2) = pair(jz+1,ja)
          pairx(3) = pair(jz+2,ja+2)
          pairx(4) = pair(jz+1,ja+1)
          pairx(5) = pair(jz+1,ja+2)
          pairx(6) = pair(jz+2,ja+1)
          pairx(7) = pair(jz+4,ja+3)
          pairfis=pair(jz,ja)
          pairgam=pair(jz,ja)
       return
       end
c
       subroutine seltzer(iz,ia,il,fisbar,delr)
      common/shell/shelmol
      common/pl3/jx,nepr,eb(100),rcsp(100),zee,amass,plex
      common/outlim/limout,lim2
       common/fizz/fizbar(106,260),fizcs(100)
c subroutine to add MC fission) to agam.f code,begsn Nov. 5,2002
c known barriers will be entered on Z and neutron number in array fizbar(z,a)
c these will be scaled with angular momentum according to Sierk or RLD mass
c formula type routines.If there are no known barriers for relevant nuclide,
c the angular momentum dependent barriers of Sierk will be used (barfit1),
c or for Z.lt.20,the old RLD results will be used.
c
c        write(*,*)'in selt fizber(iz,ia)=',fizbar(iz,ia)
        if(fizbar(iz,ia).eq.90.)go to 10
c transfer to 10 for mass formula calculation,i.e. no entry in exptl table
c now for barriers from the table,we will scale with angular momentum
       call barfit1(iz,ia,0,ibfis,irot,ilmax)
c      call barfit(iz,ia,0.,ibfis,irot,ilmax)
       bar0=float(ibfis)/10000.
       delr0=float(irot)/10000.
       call barfit1(iz,ia,il,ibfis,irot,ilmax)
c      call barfit(iz,ia,il,ibfis,irot,ilmax)
       barl=float(ibfis)/10000.
       delr=float(il)/10000.
c      delr=float(irot)/10000.
c       write(*,*)'seltzer bar0,bar1',bar0,bar1
c now scale
      fisbar=fizbar(iz,ia)*barl/bar0
       if(iz.lt.20)fisbar=fizbar(iz,ia)
      return
c now do case of exptl barrier not in table
   10  continue
       if(iz.ge.20)then
       call barfit1(iz,ia,il,ibfis,irot,ilmax)
c      call barfit(iz,ia,il,ibfis,irot,ilmax)
       fisbar=float(ibfis)/10000.
       delr=float(irot)/10000.
          endif
        return
c11       in=ia-iz
c      call mlnix(iz,in,1)
c      if(iz.lt.88)shelmol=shelmol*0.7
c      fisbar=fisbar-shelmol
c      else
c      a=ia
c      z=iz
c      an=a-z
c      al=il
c      call fisrot(a,z,an,al,delr,delsp,ero,barfac)
c     fisbar=delsp-delr
c     endif
c     write(*,*)'in selt,fisbar=',fisbar
c     return
      end
c
c
       subroutine isoevent(efinal,jl)
c integrates event mode to get ground state and isomer yields
      common/outlim/limout,lim2
      common/ismrdat/isoid(106,260),levnm(1000),bevj(1000,3)
      common/scat/jz,ja, probxx(2,2),q1(8,1100),fcs,eloss,fcspec
      common/isomer1/xgs(690),xiso(690),x2iso(690),isonum,isoint,iso2num
      common/isomer2/zs(690),as(690),eiso(690),eiso2(690)
      common/isomer3/iziso(690),iaiso(690),csgs(690),csis(690)
     1,csgse(690,100),csise(690,100),csis2e(690,100),rm1(690,100),
     2rm2(690,100)
c now branch if an isomer exists
        isno=isoid(jz,ja)
c july 04 use isoid as ident
        spinav=(xgs(isno)+xiso(isno))/2.
        ispinav=(2.*spinav+1.)
        if(xgs(isno).gt.xiso(isno))go to 9004
c branch for gs spin gt isomer spin
        if(efinal.lt.eiso(isno).or.jl.le.ispinav)then
        csgs(isno)=csgs(isno)+fcs
        else
        csis(isno)=csis(isno)+fcs
        endif
        return
 9004  continue
c branch here if gs spin is greater than isomer spin
        if(efinal.lt.eiso(isno).or.jl.gt.ispinav)then
        csgs(isno)=csgs(isno)+fcs
        else
        csis(isno)=csis(isno)+fcs
        endif
       continue
c
       return
       end
c
      subroutine rdinv(m3,jmax)
      common/outlim/limout,lim2
      common/lab3/sig(8,999),slav(8,999)
      write (6,680)
      do 30 k = 1,m3
      read(16,625)(sig(k,je),je=1,jmax)
      go to (895,805,815,825),k
  895 write(17,800)
  800 format(/,38x,'user provided neutron inverse cross sections'/)
      go to 835
  805 write(17,810)
  810 format(//38x,'user provided proton inverse cross sections'/)
      go to 835
  815 write(17,820)
  820 format(//39x,'user provided alpha inverse cross sections'/)
      go to 835
  825 write(17,830)
  830 format(//37x,'user provided deuteron inverse cross sections'/)
  835 write(17,840) jmax
  840 format('0inverse reaction cross sect. for e = 1 to',i3,' mev')
      write(17,845) (sig(k,je),je=1,jmax)
  845 format ( 1h ,10f7.0)
  680 format(//33x,'inverse reaction cross sections were provided by use
     1r'/)
  625 format(6e10.3)
  30  continue
      return
      end
      SUBROUTINE FISBAR(A,Z,BF)
c============================================
c	FISSION BARRIER BY M.MEBEL  MAY/96
c============================================
        COMMON /BLOCX/  IPARQ
        COMMON /FUSR/   BARR,SMASS,SHLL
        IPARQ=0
         NZ=INT(Z)
          NA=INT(A)
        CALL LYMASM(NZ,NA,CMASS,CBARR,NOBARR)
       BF=CBARR
          IF(NOBARR.EQ.1) BF=BARR
c-------
      RETURN
      END
      SUBROUTINE LYMASM(IZ,IA,CMASS,CBARR,NOBARR)
      COMMON /BLOCX/ IPARQ
      COMMON /FFSS/  ENEX
      COMMON /FUSR/   BARR,SMASS,SHLL
 
      DIMENSION EM(10),EMP(10),XK(10),Y(2),F(2)
      DATA ZVT/1.6666666666/, ZT/.3333333333/,
     1 ZTT/.6666666667/, SR5/2.2360679775/
      DATA EM/0.00,2.00,8.00,14.00,28.00,50.00,82.00,126.00,
     1        184.00,258.00/
      DATA CAY1/0.0/, CAY2/0.0/, CAY3/2.0/, CAY4/11.0/, CAY5/8.07144/,
     * CAY6/7.28899/,
     *            D/.444/, C/5.8/, SMALC/.325/
C
C     DMASS = REMAINDER AFTER CM - NO SHELL EFFECTS SUBTRACTED
C     SHLL = CALCULATED SHELL EFFECT
C     DIFMAS= DMASS - SHLL
C
C------------------------------
      A1 = 15.4941
C..... IPARQ=0  PARAMETRS MYERS-SWIATECKI
      IF(IPARQ.NE.0) GO TO 121
      A2 = 17.9439
      A3 = 0.7053
      GAMMA = 1.7826
      GO TO 126
 121  CONTINUE
C..... IPARQ=1  PARAMETRS KRAPPE-NIX
      IF(IPARQ.NE.1) GO TO 123
      A2 = 24.70
      A3 = 0.74476032
      GAMMA = 4.0
      GO TO 126
 123  CONTINUE
C..... IPARQ=2  PARAMETRS PAULI-LEDERGERBER
      IF(IPARQ.NE.2) GO TO 124
      A2 = 19.008
      A3 = 0.720
      GAMMA = 2.840
      GO TO 126
 124  CONTINUE
C..... IPARQ=3  BF(T)   PARAMETRS MYERS-SWIATECKI
      ALEVEL=0.1
      AMPAR=ALEVEL*FLOAT(IA)
      TSQ=ENEX/AMPAR
      A2 = 17.9439*(1.-0.0063157*TSQ)
      A3 = 0.7053*(1.-0.001*TSQ)
      GAMMA = 1.7826
 126  CONTINUE
C------------------------------
      IF(IZ.NE.0) GO TO 15
      CMASS=0.0
      RETURN
   15 NOBARR=0
      DO 1 I=1,10
      EMP(I)=EM(I)**ZVT
    1 CONTINUE
      DO 2 I=1,9
      XK(I)=.600*(EMP(I+1)-EMP(I)) /(EM(I+1)-EM(I))
    2 CONTINUE
C
C     FOR DEFINITIONS OF CAY1 AND RZ,SEE UCRL-11980
C
      CAY1=3.28637900*A3**3
      RZ=.86398700/A3
      Z= FLOAT(IZ)
      ZSQ=Z**2
      N=IA-IZ
      UN= FLOAT(N)
      A= FLOAT(IA)
      A3RT=A**ZT
      A3RT2=A3RT**2
      A2RT= SQRT(A)
      SYM=((UN-Z)/A)**2
      ACOR=1.00-GAMMA*SYM
      PARMAS=CAY5*UN+CAY6*Z
      VOLNUC=-A1*ACOR*A
      SUFNUC=A2*ACOR*A3RT2
      COULMB=A3*ZSQ/A3RT
      FUZSUR=-CAY1*ZSQ/A
      ODDEV=-(1.00+2.00*(N/2)-UN+2.00*(IZ/2)-Z)/A2RT*CAY4
      WTERM=-CAY2*A3RT2* EXP(-CAY3*SYM)
      WOTNUC=PARMAS+COULMB+FUZSUR+ODDEV+WTERM
      SMASS=WOTNUC+VOLNUC+SUFNUC
      SPW=SUFNUC+WTERM
      C2=SPW/A3RT2
      X=.5*COULMB/SPW
      IF(X.GE.1.00) GO TO 4
C------------------------------
      ARQ=X
      IF(IPARQ.EQ.0) BARR=SUFNUC*XI(X)
      IF(IPARQ.EQ.2) BARR=SUFNUC*XI(X)
      IF(IPARQ.EQ.3) BARR=SUFNUC*XI(X)
      IF(IPARQ.EQ.1) BARR=SUFNUC*XIMOD(X)
      IF(IPARQ.EQ.0) F1Q=XI(X)
      IF(IPARQ.EQ.1) F1MQ=XIMOD(X)
C------------------------------
      GO TO 6
    4 BARR=0.0
    6 Y(1)=UN
      Y(2)=Z
      DO 31 J=1,2
      DO 32 I=1,9
      IF(Y(J).le.EM(I+1))go to 3 
   32 CONTINUE
      PRINT 332,J
  332 FORMAT('1FAILURE IN LYMASS - Y(',I1,') EXCEEDS LAST MAGIC NO.')
      STOP
    3 F(J)=XK(I)*(Y(J)-EM(I))-.600*(Y(J)**ZVT-EMP(I))
   31 CONTINUE
      S=(2.00/A)**ZTT*(F(1)+F(2))-SMALC*A3RT
      C2D2=C2*D**2
      EE=(C2D2+C2D2)*(1.00-X)
      FF=.425917710*C2D2*D*(1.00+X+X)/A3RT
      SSHELL=C*S
      V=SSHELL/EE
      EPS=1.500*FF/EE
      IF(EE*(1.00-3.00*V).LE.0.00) GO TO 51
      QCALC=0.00
      THETA=0.00
      SHLL=SSHELL
      GO TO 52
C
C       ESTIMATE THETA
C
   51 TO=1.00
C
C       ITERATE TO FIND EQUILIBRIUM THETA
C
  101 DO 725 IPQ=1,10
      TO2=TO**2
C----------------------------------------
C     IF (TO2.GT.170.) PRINT 500, IZ,IA
C 500 FORMAT(1X,'LYMASM',2X,2I5)
C----------------------------------------
C     EXMT2= EXP(-TO2)
      EXMT2= 1.E-20
      IF((ABS(TO2)).LT.30.) EXMT2= EXP(-TO2)
C
      T=TO-(1.00-EPS*TO-V*(3.00-TO2-TO2)*EXMT2) /
     1(-EPS+V*TO*(10.00-4.00*TO2)*EXMT2)
      IF(T.LE.0.00) GO TO 728
      IF( ABS(T-TO) .LT.1.E-4) GO TO 732
      TO=T
  725 CONTINUE
      GO TO 729
  732 T2=T**2
C     EXT2= EXP(-T2)
      EXT2= 1.E-20
      IF((ABS(T2)).LT.30.) EXT2= EXP(-T2)
C
      TEST=EE*(1.00-EPS*(T+T)-V*((4.00*T2-12.00)*T2+3.00)* EXT2)
      IF(TEST.GT.0.00) GO TO 81
  728 TO=.100
      DO 100 I=1,20
      TO2=TO**2
      GL=EE*(1.00-EPS*TO-V*(3.00-TO2-TO2)* EXP(-TO2))
      IF(GL.GT.0.00) GO TO 101
  100 CONTINUE
  729 CMASS=SMASS
      CBARR=0.00
      NOBARR=1
      RETURN
   81 THETA=T
      ALPHA0=D*SR5/A3RT
      ALPHA=ALPHA0*THETA
      SIGMA=ALPHA*(1.00+ALPHA/14.00)
      EXS= EXP(SIGMA+SIGMA)- EXP(-SIGMA)
      QCALC=4.E-3*Z*(RZ*A3RT)**2*EXS
      T2=T**2
      SHLL=T2*(EE-FF*T) + SSHELL*(1.00-T2-T2)* EXP(-T2)
   52 CMASS=SMASS+SHLL
      CBARR=BARR-SHLL
      RETURN
      END
      FUNCTION XI(Z)
C
C     6-POINT LAGRANGE INTERPOLATION
C
      DIMENSION Y(51)
      DATA Y/.25900,.255200,.250700,.245100,.2400,.23400,.228500,
     1      .22200,.21600,.2100,.20300,.196800,.1900,.18300,.175800,
     2      .1692400,.1620300,.1547800,.147500,.1401900,.1328400,
     3      .1254500,.1180100,.1105200,.1029600,.0953500,.0876800,
     4      .0799900,.0722900,.064600,.0569500,.0493700,.0419300,
     5      .0347600,.0281100,.0223600,.0176200,.0137300,.0105600,
     6      .0079800,.0059100,.0042500,.0029600,.0019700,.0012300,
     7      7.1E-4,3.6E-4,1.5E-4,4.E-5,1.E-5,0.00/
C
C     THE X VALUES ARE EVENLY SPACED - X = 0(.02)1
C
      ZBH=Z*50.00
      M=IFIX(ZBH)
      DEL=ZBH- FLOAT(M)
      M=M+1
      IF(M.LE.51) GO TO 105
      M=51
  100 XI=Y(M)
      RETURN
  105 IF (DEL.LT.1.E-4) GO TO 100
      IF(M.GE.3) GO TO 110
      DEL=DEL- FLOAT(3-M)
      M=3
      GO TO 115
  110 IF(M.LE.48) GO TO 115
      DEL=DEL+ FLOAT(M-48)
      M=48
  115 DM3=DEL-3.00
      PROD=DM3*DEL
      W6=1.00/(1.2E2*DM3)
      DM2=DM3+1.00
      PROD=DM2*PROD
       W5=-1.00/(24.00*DM2)
      DM1=DM2+1.00
      PROD=DM1*PROD
      W4=1.00/(12.00*DM1)
      DP1=DM1+2.00
      PROD=DP1*PROD
      W2=1.00/(24.00*DP1)
      DP2=DP1+1.00
      PROD=DP2*PROD
      W1=-1.00/(1.2E2*DP2)
      W3=-1.00/(12.00*DEL)
      XI=PROD*(W1*Y(M-2)+W2*Y(M-1)+W3*Y(M)+W4*Y(M+1)+W5*Y(M+2)
     1 +W6*Y(M+3))
      RETURN
      END
      FUNCTION XIMOD(Z)
C
C     6-POINT LAGRANGE INTERPOLATION
C      IN MODIFIED LIQUID-DROP FORMULA
C         ( KRAPPE [ NIX --  IAEA-SM-174/12 )
C
      DIMENSION Y(51)
      DATA  Y/
     1    0.12200, 0.12100, 0.11980, 0.11830, 0.11690, 0.11520, 0.1133,
     2    0.11130, 0.10900, 0.10670, 0.10420, 0.10150, 0.09850, 0.09540,
     3    0.09180, 0.08780, 0.08350, 0.07900, 0.07460, 0.06960, 0.06470,
     4    0.05960, 0.05420, 0.04880, 0.04350, 0.03880, 0.03400, 0.02920,
     5    0.02460, 0.02020, 0.01580, 0.01220, 0.00900, 0.00660, 0.00490,
     6    0.00360, 0.00280, 0.00220, 0.00180, 0.00140, 0.00100, 0.00090,
     7    0.00060, 0.00040, 0.00020, 0.00010, 0.00000, 0.00000, 0.00000,
     8    0.00000, 0.00000/
C
C     THE X VALUES ARE EVENLY SPACED - X = 0(.02)1
C
      ZBH=Z*50.00
      M=IFIX(ZBH)
      DEL=ZBH- FLOAT(M)
      M=M+1
      IF(M.LE.51) GO TO 105
      M=51
  100 XIMOD=Y(M)
      RETURN
  105 IF (DEL.LT.1.E-4) GO TO 100
      IF(M.GE.3) GO TO 110
      DEL=DEL- FLOAT(3-M)
      M=3
      GO TO 115
  110 IF(M.LE.48) GO TO 115
      DEL=DEL+ FLOAT(M-48)
      M=48
  115 DM3=DEL-3.00
      PROD=DM3*DEL
      W6=1.00/(1.2E2*DM3)
      DM2=DM3+1.00
      PROD=DM2*PROD
       W5=-1.00/(24.00*DM2)
      DM1=DM2+1.00
      PROD=DM1*PROD
      W4=1.00/(12.00*DM1)
      DP1=DM1+2.00
      PROD=DP1*PROD
      W2=1.00/(24.00*DP1)
      DP2=DP1+1.00
      PROD=DP2*PROD
      W1=-1.00/(1.2E2*DP2)
      W3=-1.00/(12.00*DEL)
      XIMOD=PROD*(W1*Y(M-2)+W2*Y(M-1)+W3*Y(M)+W4*Y(M+1)+W5*Y(M+2)
     1 +W6*Y(M+3))
      RETURN
      END
      subroutine optout
      common/outlim/limout,lim2
      common/shft/mc,mp,inver,ike,ipch,kq,qval,ap,at,zp,zt,cld,barfac
      common/ug/iz ,ia
      common/nhy2/gdo,bisp
      common/pl2/sux(24,13,100)
      common/s1/jfrac,jupper
      common/sf/m3,kplt
      common/sft9/k5,jmax,pld,c(8),delt,alt
      common/par2/cncss
      common/pl8/jamma(15,24),nulim(15,24)
      common/pl1/ecm(100)
      common/nam/k2
      common/hjk/jang,iq,rcss
      common/ss/sor,rr
      common/nr34/nr3,nr4,ke5,i3d,irtst,i3t2,ijkl,lq,tem(36)
      common/parfs/mz,k6,delrr(999)
      common/shft4/k3,jcal
      common/shft5/fiss
      common/shft2/fs(25),dsp(25),brr(25),der(25),er(25)
     1,cx(25),bilsum,xmiss,sumiz
      REAL*8 fiss
      common/pl4/crs(24,13)
      common/nhy/ij,jl,ji,jj,td,ex1,ex2,tmx,av,gav,cost,b(3)
      common/paro/pq,cross
      common/hyb2/pp(15,24,999)
      common/pl3/jx,nepr,eb(100),rcsp(100),zee,amass,plex
      common/sft5/exc(18,27),xmxx
      common/lym1/be(18,27,8),pair(18,27),xmas(18,27),delshl(18,27),amas
     1(18,27)
c     common/bym1/symb(18,27),symbp(18,27),symbs(18,27)
      common/pl7/mas
      common/lab3/sig(8,999),slav(8,999)
      common/pl5/na,nz,title(20)
      common/scr/k9
      common/fis/ifis
       common/scrn/att(10),cldd(10),bexp(24,15),reacss(100),edelt
      common/iso/qpn(3),qpnc
      common/incr/ed,tened,efinal
      common/levop/ldopt,igam
      common/tst/test
      common/send/irfr,iadst,dlt
      common/gam/io,la,ppgam(999),gmspec(999),beta,ams
      common/ist/isot,iikl,knz,idelm,fracts,fract,ccrs(32,13,100)
     1,engbuf(100),exd,een(8,999,100)
c
 1000 format(' qval=0. ; q value will be calculated internally')
 1001 format(' projectile mass no = ',f5.1,
     1'  projectile atomic no =',f5.1,/,
     2' target mass no =      ',f5.1,
     3' target atomic no =     ',f5.1/)
 1002 format(' cld =',f5.1,'   ratio of nucleon to fission level density
     3 parameters')
 1003 format(' cld =',f5.3,'   isotopic abundance of this isotope')
 1004 format(' na = ',i5,'   calculate up to na-1 neutrons out')
 1005 format(' nz = ',i5,'   calculate up to nz-1 charges out')
 1030 format(//,' *** input options selected ***',//)
 1031 format(' ifis =   0   rotating finite range barriers of A.J. Sierk
     1 used')
 1032 format(' ifis =',i5,'  rotating liquid drop barriers of C-P-S used
     1')
 1033 format(' barfac=',f5.1,' multiplier of fission barrier')
c1006 format(' mc=0   myers swiatecki lysekil (msl) masses incl.',
c    1' shell corr.')
c    2'       and/or ldgs at least partly provided by user (see below)')
c1008 format(' mc = ',i5,'   experimental masses used where available (s
c    2ee comment cards in main for description of mc, mp options)')
c1009 format(' mp =     0   no pairing term in masses')
c1010 format(' mp =     1   pairing term in masses, ldgs calculated from
c    1 msl formula and applied backshifted')
c1011 format(' mp =     2    pairing term in masses, ldgs calculated fro
c    1m msl formula with shell correction and backshifted')
c1012 format(' mp =     3   normal pairing shift, zero for odd-even nucl
c    1ei, delta added to excitation for odd-odd nuclei, etc.')
c1013 format(' inver =  0   optical model used to calculate inverse cros
c    1s sections')
c1014 format(' inver =  1   user supplied inverse cross sections')
c1015 format(' inver =  2   sharp cutoff method used to calculate',
c    1'inverse cross sections')
c1016 format(' pld =',f5.1,'   level density parameter = acn /',f4.1)
c1017 format(' ike =    0   no particle spectra will be printed')
c1018 format(' ike =    1   equilibrium spectra for each nuclide will be
c    1 printed')
c1019 format(' ike =    2   only precompound spectra printed')
c1020 format(' ike =    3   equilibrium spectra for each nuclide and pre
c    1compound spectra will be printed')
c1021 format(' ike =    4   precompound spectra will be printed as well
c    1as the sum',/,'(over all emitting nuclides and all partial waves)
c    2of precompound plus equilibrium spectra')
c1022 format(' m3 =     1   neutron emission only')
c1023 format(' m3 =     2   neutron and proton emission')
c1024 format(' m3 =     3   neutron, proton and alpha emission')
c1025 format(' m3 =     4   neutron, proton, alpha and deuteron emission
c    3')
c1026 format(////////////)
c1027 format(' ldopt =  0   fermi gas level density')
c1028 format(' ldopt =  1   kataria ramamurthy level density option sele
c    1cted - mp set to 1, backshifted pairing required')
c1029 format(' ldopt =  2   ignatyuk level density')
c10291 format(' ldopt = 3,gilbert-cameron level densities per p.sarkar')
c10292 format(' ldopt = 5,light element densities m.chadwick')
c2029 format(' ldopt = 4 fermi gas level density computed each a,z')
c1034 format(' ed = ',f5.2,'   energy bin mesh size in mev')
c1035 format(' ipch =   1  fission barriers read in after card #1 and ar
c    1e independent of ang. momentum')
c1036 format(' ipch =   2  fission barriers read in after card #1 and ar
c    1e scaled as rotating liquid drop model barriers')
c1037 format(' kplt =   1   excitation functions will be plotted on stan
c    1dard output if last energy input line is followed by -1., col 1-5'
c    1)
c1038 format('              gamma spectra will be printed')
c
c1039 format(' inver =  3   inverse cross sections from Chadwick file')
      write(17,1030)
      if(qval.eq.0.)write(17,1000)
      write(17,1001)ap,zp,at,zt
      if(isot.eq.1)write(17,1002)cld
      if(ifis.eq.0)write(17,1031)
      if(ifis.gt.0)write(17,1032)ifis
      write(17,1033)barfac
      if(isot.gt.1)write(17,1003)fract
      write(17,1004)na
      write(17,1005)nz
c     if(mc.eq.0)write(17,1006)
cc     if(mc.eq.1)write(17,1007)
cc     if(mc.ge.10)write(17,1008)mc
c     if(mp.eq.0)write(17,1009)
c     if(mp.eq.1)write(17,1010)
c     if(mp.eq.2)write(17,1011)
c     if(mp.eq.3)write(17,1012)
c     if(ldopt.eq.0)write(17,1027)
c     if(ldopt.eq.1)write(17,1028)
c     if(ldopt.eq.2)write(17,1029)
c     if(ldopt.eq.4)write(17,2029)
c     if(ldopt.eq.3)write(17,10291)
c     if(ldopt.eq.5)write(17,10292)
c     if(inver.eq.0)write(17,1013)
c     if(inver.eq.1)write(17,1014)
c     if(inver.eq.2)write(17,1015)
c     if(inver.eq.3)write(17,1039)
c     write(17,1034)ed
c     if(ike.eq.0)write(17,1017)
c     if(ike.eq.1)write(17,1018)
c     if(ike.eq.2)write(17,1019)
c     if(ike.eq.3)write(17,1020)
c     if(ike.eq.4)write(17,1021)
c     if(igam.eq.1)write(17,1038)
c     if(ipch.eq.1)write(17,1035)
c     if(ipch.eq.2)write(17,1036)
c     write(17,1016)pld,pld
c     if(kplt.eq.1)write(17,1037)
c     if(m3.eq.1)write(17,1022)
c     if(m3.eq.2)write(17,1023)
c     if(m3.eq.3)write(17,1024)
c     if(m3.eq.4)write(17,1025)
c
      return
      end
 
      subroutine fermld
      common/levop/ldopt,igam
      common/outlim/limout,lim2
      common/sft5/exc(18,27),xmxx
      common/lym1/be(18,27,8),pair(18,27),xmas(18,27),delshl(18,27),amas
     1(18,27)
      common/sf/m3,kplt
      common/giani/pre, sigt(1100),sigam(999),sigpre(999),gsp(999)
     2,eqgam(999)
      common/gam/io,la,ppgam(999),gmspec(999),beta,ams
      common/gr/rint(1100)
      REAL*8 rint
      common/qq/gow(1100)
      REAL*8 gow
c     common/exclus/ gib(2,999),pr(3),ge(2100),zz(3),
c    1delta(3),ts(3),pairx(8),asum(10),rate(2,999)
c     common/q/ q(999),sp(999),sifis(999),gamft(999),t(999)
c    1,scale(24,15),rd(4)
      common/q/ q(999),sp(999),sifis(999),gamft(999),t(999),pairx(8)
     1,scale(24,15),rd(4)
      common/cs/crsum
      common/lab10/pow(4,9999),gam(9999)
      REAL*8 pow,gam
      common/shft/mc,mp,inver,ike,ipch,kq,qval,ap,at,zp,zt,cld,barfac
      common/ug/iz ,ia
      common/sft9/k5,jmax,pld,c(8),delt,alt
      common/incr/ed,tened,efinal
      common/lab12/pof(4,9999)
      REAL*8 pof
      common/ss/sor,rr
      common/ldf/gamtem(9999)
      REAL*8 gamtem
      REAL*8 sq,sr,sf
      REAL*8 temp,ex,temf
      REAL*8 cq,cq1,temp1

      amass=ap+at+2.-float(ia)-float(iz)
c     sor=0.7*sqrt((ap+at)*100.)
c   june 14 2017 temp put sor=0. to verify no change...
      sor=0.
c.....................................................
       if((ap+at).gt.50.)sor=sor/2.
      cld=1.06
      rf=amass*cld/pld
      rdf=amass/(2.*pld)
      rd(1)=(amass-1.)/pld
      rd(2)=(amass-1.)/pld
      rd(3)=(amass-4.)/pld
      rd(4)=(amass-2.)/pld
      rdd=rd(1)*0.5
      rdc=1.000
      if(amass.lt.200.)rdc=(200./amass)**0.333333
      do 7050 l=1,4
      do 7050 ib=1,9999
      bj=ib
      bj=bj/10.-.05
      sq=sqrt(rd(l)*bj)
      sq=sq+sq-sor
        if(l.eq.1)then
      sf=sqrt(rf*bj)
      sf=sf+sf-sor
c add the sor correction in line above on6/8/2005
      pof(l,ib)=(1./(1.+bj**1.25))*dexp(sf)
      endif
 7050 pow(l,ib)=(1./(1.+bj**1.25))*dexp(sq)
      do 7051 ib=1,999
      bj=ib
      bj=bj/10.-.05
      sq=sqrt(rdd*bj)
      sq=sq+sq-sor
      gow(ib)=(1./(1.+bj**1.25))*dexp(sq)
 7051 continue
      const=33.0
      sfr=sor
       rdc=1.
      rf=1.28*amass/pld
      do 7055 ib=1,9999
      if(ib.gt.99)rf=1.26*amass/pld
      if(ib.gt.140)rf=1.23*amass/pld
      if(ib.gt.210)rf=1.21*amass/pld
      if(ib.gt.280)rf=1.18*amass/pld
      if(ib.gt.380)rf=1.16*amass/pld
      if(ib.gt.480)rf=1.12*amass/pld
      if(ib.gt.580)rf=1.09*amass/pld
      if(ib.gt.660)rf=1.08*amass/pld
       if(ib.gt.920)rf=1.07*amass/pld
      if(ib.gt.1200)rf=1.06*amass/pld
     
      bj=ib
      bj=bj/10.-.05
      sr=sqrt(rf*bj*rdc)
      sr=sr+sr-sfr
      gam(ib)=const*(1./(1.+bj**1.25))*dexp(sr)
 7055 gamtem(ib)=gam(ib)
c
      nia=amass
      ex=0.5
c
      do 7047 l=1,3
      temp=sqrt(ex/rd(l))
      temf=sqrt(ex/rdf)     
      ten=10.*ex+0.5
      it=ten
      if=it-1
      d=(1./temp)*dexp((ex)/temp)
      cq=pow(l,it)/d
      cq=cq/temp
      f=(1./temf)*dexp((ex)/temf)
      cf=pof(l,it)/f
      cf=cf/temf
      do 7048 ib=1,if
      e=0.10*float(ib)-0.05
      pow(l,ib)=cq*dexp((e)/temp)
      if(l.eq.1)then
      pof(l,ib)=cf*dexp((e)/temf)
      endif
7048  continue
7047  continue
      rf=1.11*amass/pld
      temp1=sqrt(ex/rf)
      d1=(1./temp1)*dexp((ex)/temp1)
      temp=sqrt(ex/rdd)
      ten=10.*ex+0.5
      it=ten
      if=it-1
      d=(1./temp)*dexp((ex)/temp)
      cq=gow(it)/d
      cq=cq/temp
      cq1=gam(it)/d1
      cq1=cq1/temp1
      do 27048 ib=1,if
      e=0.10*float(ib)-0.05
      gow(ib)=cq*dexp((e)/temp)
      gam(ib)=cq1*dexp((e)/temp1)
      gamtem(ib)=gam(ib)
c     gow(ib)=pow(1,ib)
27048  continue
      continue
      do 37047 ib=2,9999
      gam(ib)=gam(ib-1)+gam(ib)
37047 continue
      do 37048 ib=1,9999
      gam(ib)=gam(ib)*0.1
37048 continue
      return
      end
      subroutine gama
      common/outlim/limout,lim2
      common/shft/mc,mp,inver,ike,ipch,kq,qval,ap,at,zp,zt,cld,barfac
      common/giani/pre, sigt(1100),sigam(999),sigpre(999),gsp(999)
     1,eqgam(999)
      common/incr/ed,tened,efinal
      common/qq/gow(1100)
      REAL*8 gow
      common/gam/io,la,ppgam(999),gmspec(999),beta,ams
      common/ug/iz,ia
      common/pl5/na,nz,title(20)
      common/hjk/jang,iq,rcss
      common/sft5/exc(18,27),xmax
      common/jpn/sgjpn(999)
 
      data alev0/6.5/
      alev=6.5
      an=at-zt
      beta=0.
      max=30./ed
      if(pre.eq.1.)max=999
      if(pre.eq.0.)go to 3000
      mx=max
      do 1 i=1,999
      gsp(i)=0.
   1  continue
c     define lorentzian constants and squares of values
c     e0=43.4*ams**(-0.215)
c     e1=e0*(1.-beta/3.)**2
c     e1sq=e1*e1
c     sig1=0.0145*ams/e1
c     gam1=0.232*e1
c     gam1sq=gam1*gam1
c     e2=e0*(1.-0.16*beta)
c     e2sq=e2*e2
c     sig2=0.0235*ams/e2
c     gam2=0.275*e2
c     gam2sq=gam2*gam2
c     set up spectral shape on maximum excitation energy
c     below particle binding energy
       er=77.9*(1./ams**0.33333)*(1.-exp(-ams/238.))+(34.5/ams**
     10.166666)*exp(-ams/238.)
      do 2 ieps=1,max
c     eps is gamma ray energy in mev
      eps=float(ieps)*ed
      egam=eps
      epsq=eps*eps
c     sigt(ieps)=sig1*epsq*gam1sq/((epsq-e1sq)**2+epsq*gam1sq)
c    1+sig2*epsq*gam2sq/((epsq-e2sq)**2+epsq*gam2sq)
c     sigt(ieps)=sigt(ieps)*epsq
c begin gdr calculation using systematics from japan
c      gamtot=0.026*eps**1.91
c      gamtot=0.026*er
c   we use the berman-fultz expression for the single lorentzian
c   resonance energy(mev),and the form of peak cross section due to
c   berman and fultz,but with a width factor of o.0.275*e-resonance,a
c   value taken from the second width in reffo's two sum formulation.
c   this is a purely empirical,ad-hoc choice.
       gamtot=0.275*er
      tem=epsq*gamtot*gamtot
      sigr=1.2*60.*(an*zt/at)*(2./(3.14159*gamtot))
      sgjpn(ieps)=tem*sigr/(epsq*gamtot*gamtot+(epsq-er*er)**2)
c the above defines the single lorentzian giant dipole cross section(millibarns)
c
c     sigt(ieps)=epsq*sgjpn(ieps)*0.001
c     begin quasideuteron contribution to reaction cross section
c     using routine from mb chadwick based on chadwick,prc44,814(1991)
      aa=at
      an=at-zt
      az=zt
c   define the free deuteron photodisintegration cross section in mb
      if(egam.lt.2.224)qd=0.
      if(egam.lt.2.224)go to 20
      sigdeut=61.2*((egam-2.224)**1.5)/(egam**3)
c     define pauli blocking function
      f=8.3714e-02-9.8343e-3*egam
     1+4.1222e-04*egam**2-3.4762e-6*egam**3
     2+9.3537e-9*egam**4
c the above expression is valid only between 20 and 140 mev.for
c other energies chadwick matches on a levinger type exponential
c damping:
      if(egam.lt.20.)f=exp(-73.3/egam)
      if(egam.gt.140)f=exp(-24.2348/egam)
c now calculate the nuclear photoabsorption x/s via the qd mechanism
      qdxs=(alev*an*az/aa)*sigdeut*f
c
c calculate e gamma squared times cross section in barns
   20 if(egam.lt.2.224)qdxs=0.
      sigt(ieps)=0.001*epsq*sgjpn(ieps)+qdxs*epsq*0.001
    2 continue
c     do 140 ie=1,100
c     write(17,141)ie,sgjpn(ie),sigt(ie)
c 140 continue
c 141 format(' jpn gdr energy =',i5,' cross section mb=',1e10.3,
c    1'  reffo value= ',1e10.3)
      if(pre.eq.1.)return
 3000 continue
      mx=30./ed+0.01
      if(mx.gt.999)mx=999
      do 4 je=1,mx
      ie=mx-je+1
      if(ie.lt.1)go to 4
      if(ppgam(ie).eq.0.)go to 4
      e=float(ie     )*ed
c     e is energy in mev of a bin which will emit a gamma cascade
      smgm=0.
c     do normalization summation for emission of all gamma energies
      map=ie-1
      if(map.lt.1)go to 4
      do 3 ieps=1,map
      eps=float(ieps)*ed
      u=e-eps
      bi=10.*u+.5
      ib=bi
      if(ib.le.0)go to 3
      if(ib.gt.999)ib=999
      sigam(ieps)=sigt(ieps)*gow(ib)
      smgm=smgm+sigam(ieps)
    3 continue
      do 5 ieps=1,map
      spec=ppgam(ie     )*sigam(ieps)/smgm
      gmspec(ieps)=gmspec(ieps)+spec
      if((ie-ieps).le.0)go to 5
      ppgam(ie-ieps)=ppgam(ie-ieps)+spec
    5 continue
    4 continue
      bnzer=ppgam(1)
      sgm=0.
      do 14 i=1,mx
      sgm=sgm+gmspec(i)
   14 continue
      sgm=sgm+bnzer
      write(54,15)sgm,bnzer
   15 format(/'  summed gamma ray cross section = ',e10.3,
     1//,'  gamma spectrum summed over all emitting nuclides',e10.3)
       write(54,*)' Egamma(MeV)  gamma cs mb/MeV  gamma cs mb/MeV.sr'
      do 13 i=1,max
      egam=ed*float(i)
      gmspec(i)=gmspec(i)/ed
      gsp(i)=gmspec(i)/12.6
      write(54,*)egam,gmspec(i),gsp(i)
   13 continue
      mmm=(xmax+ed)/ed
      write(54,301)
      if(mmm.gt.999)mmm=999
      do 200 jx=1,mmm
      eqgam(jx)=eqgam(jx)/ed
      sigpre(jx)=sigpre(jx)+gmspec(jx)+eqgam(jx)
      gsp(jx)=sigpre(jx)/12.6
      egam=ed*float(jx)
      write(54,300)egam,gmspec(jx),eqgam(jx),sigpre(jx),gsp(jx)
  200 continue
  300 format(4x,f6.2,3x,4(2x,e10.2,2x))
  301 format(///'e gamma(mev) bound equilib unbound equil sum pre+equil
     1sum/4pi(mb/mev.sr)'//)
c 600 continue
c 800 continue
      return
      end
      subroutine photin
       common/scrn/att(10),cldd(10),bexp(24,15),reacss(100),edelt
      common/pl3/jx,nepr,eb(100),rcsp(100),zee,amass,plex
      common/outlim/limout,lim2
      common/shft/mc,mp,inver,ike,ipch,kq,qval,ap,at,zp,zt,cld,barfac
      common/par3/eq,sigml(999),acrs(999)
      common/incr/ed,tened,efinal
      common/hjk/jang,iq,rcss
      common/sft5/exc(18,27),xmax
      common/sigsav/sigdr(100),sigqd(100)
      common/deform/beta2(106,260)
 
      data alev0/6.5/
      alev=6.5
      if(alev.lt.0.0001)alev=alev0
      ams=ap+at
      zee=zp+zt
      iz=zee
      ia=ams
      beta=beta2(iz,ia)
      an=at-zt
      isw=0
      do 130 i=1,5986
      if(isw.eq.1)go to 130
      read(38,900)izz,iaa,el,eta,e1,w1,e2,w2
      if(izz.eq.iz.and.iaa.eq.ia)isw=1
 130  continue
      rewind 38
 900  format(2i4,1x,a2,1f7.3,4f7.2)
c add Marie constant of 0.17 to e1
      e1=e1+0.17
      eps=eq
      egam=eq
c     eps is gamma ray energy in mev,as is egam
c     calculate giant dipole contribution
      epsq=eps*eps
c   we use the berman-fultz expression for the single lorentzian
c   resonance energy(mev),and the form of peak cross section due to
c   berman and fultz,but with a width factor of 0.275*e-resonance,a
c   value taken from the second width in reffo's two sum formulation.
c   this is a purely empirical,ad-hoc choice.
c 4/04 add some non-degenerate giant dipole numbers
c     e0=77.9*(1./ams**0.3333)*(1.-exp(-ams/238.))+(34.5/ams**0.1666
c    1*exp(-ams/238.))
c      e1=e0*(1.+beta*0.9)
c      sig1=(72.*an*zt/ at)*0.35 

c      gam1=(0.275-0.23*beta)*e1
c      e2=e0*(1.-beta/3.)
c      sig2=(72.*an*zt/at)*0.65
c      gam2=0.275*e2
c      er=77.9*(1./ams**0.33333)*(1.-exp(-ams/238.))+(34.5/ams**
c    10.166666)*exp(-ams/238.)
c      gamtot=0.275*er
c     tem=epsq*gamtot*gamtot
c     sigr=1.2*60.*(an*zt/at)*(2./(3.14159*gamtot))
      sigr1=1.2*60.*(an*zt/at)*(2./3.14159)/(3.*w1)
      sigr2=2.*sigr1*w1/w2
      t1=1.+(eps**2-e1**2)**2/(eps**2*w1**2)
      t2=1.+(eps**2-e2**2)**2/(eps**2*w2**2)
      siggdr=sigr1/t1+sigr2/t2
c above from Marie Giacci
c     siggdr=tem*sigr/(epsq*gamtot*gamtot+(epsq-er*er)**2)
c now try to remember form for double resonance
c      siggdr2=epsq*gam1**2*sig1/((epsq-e1**2)**2+(epsq*gam1**2))
c    1  +epsq*gam2**2*sig2/((epsq-e2**2)**2+(epsq*gam2*gam2 ))
c     write(*,*)'ams,zee,eps,s0,s2',ams,zee,eps,siggdr,siggdr2
c     write(*,*)'sig1,sig2,sigr',sig1,sig2,sigr
c     write(17,1000)er,siggdr
c1000  format('calculated resonance energy = ',f7.3,' MeV, GDR cross
c     1section contribution = ',1f8.3,' mb')
c the above defines the single lorentzian giant dipole cross section(millibarns)
c     begin quasideuteron contribution to reaction cross section
c     using routine from mb chadwick based on chadwick,prc44,814(1991)
      aa=at
      an=at-zt
      az=zt
c   define the free deuteron photodisintegration cross section in mb
      if(egam.lt.2.224)qd=0.
      if(egam.lt.2.224)go to 20
      sigdeut=61.2*((egam-2.224)**1.5)/(egam**3)
c     define pauli blocking function
      f=8.3714e-02-9.8343e-3*egam
     1+4.1222e-04*egam**2-3.4762e-6*egam**3
     2+9.3537e-9*egam**4
c the above expression is valid only between 20 and 140 mev.for
c other energies chadwick matches on a levinger type exponential
c damping:
      if(egam.lt.20.)f=exp(-73.3/egam)
      if(egam.gt.140)f=exp(-24.2348/egam)
c now calculate the nuclear photoabsorption x/s via the qd mechanism
      qdxs=(alev*an*az/aa)*sigdeut*f
c
   20 continue
 
      write(17,103)egam
      write(*,103)egam
  103 format(' photonuclear cross section for incident energy',
     11f10.3,' MeV')
      rcss=qdxs+siggdr
      sigdr(nepr)=siggdr
      sigqd(nepr)=qdxs
      reacss(nepr)=rcss
      write(*,*)'rcss for incident photon =',rcss,' mb'
      write(17,100)rcss,siggdr,qdxs
  100 format(' for projectile mass zero a reaction c.s. of ',f10.5,
     1' mb was calculated  giant dipole =',f8.3,' mb and quasi-deuteron
     2=',f10.5,' mb')
      write(17,101)
  101 format(' for photon in the entrance channel,precompound is
     1treated according to Phys.Rev.28C,2286(1983)')
      return
      end
 
      subroutine shaft
      common/astp3/niso,no,iadlt,isotag,massno(12),mass1(12),ncsave
      common/outlim/limout,lim2
      dimension ms(22),asum(10)
        common/mccom/ncount,ipot,rum(15,24),ppmc(15,24,999),cf
      common/fis/ifis
      common/rcc/rcsi(10,100)
      common/ks/ksot,mtar(10),iete(10)
      common/renorm/ssm1,ssm2,ssm3,ssm4
      common/cs/crsum
      common/memo/scrs(999,8)
      common/pl7/mas
      common/s1/jfrac,jupper
      common/pl8/jamma(15,24),nulim(15,24)
      common/ug/iz,ia
      common/ss/sor,rr
      common/shft/mc,mp,inver,ike,ipch,kq,qval,ap,at,zp,zt,cld,barfac
      common/sft5/exc(18,27),xmax
      common/lym1/be(18,27,8),pair(18,27),xmas(18,27),delshl(18,27),amas
     1(18,27)
c     common/bym1/symb(18,27),symbp(18,27),symbs(18,27)
      common/pl3/jx,nepr,eb(100),rcsp(100),zee,amass,plex
      common/par3/eq,sigml(999),acrs(999)
      common/pl5/na,nz,title(20)
      common/sf6/aoz,zoz,en(8,999)
      common/lab3/sig(8,999),slav(8,999)
      common/nhy/ij,jl,ji,jj,td,ex1,ex2,tmx,av,gav,cost,b(3)
      common/nhy2/gdo,bisp
      common/sft9/k5,jmax,pld,c(8),delt,alt
      common/sf/m3,kplt
      common/nr34/nr3,nr4,ke5,i3d,irtst,i3t2,ijkl,lq,tem(36)
      common/pl4/crs(24,13)
      common/hjk/jang,iq,rcss
      common/incr/ed,tened,efinal
      common/send/irfr,iadst,dlt
      common/shft4/k3,jcal
      common/shft5/fiss
      common/shft2/fs(25),dsp(25),brr(25),der(25),er(25)
     1,cx(25),bilsum,xmiss,sumiz
      REAL*8 fiss
      common/scr/k9
      common/tst/test
      common/dist/n,ni,ke,irefr,xz(3),pod,pcs,ek,xmix,cfrac
     1,dtau,sign(36),tmix,bji(10),fin,dtup
      common/ist/isot,iikl,knz,idelm,fracts,fract,ccrs(32,13,100)
     1,engbuf(100),exd,een(8,999,100)
      common/nuopt/kinem,iem,kb
      common/distb/dsgp(36,999,8)
      if(k3.eq.1)go to 60
      if(k3.eq.2)go to 95
      if(k3.eq.3)go to 340
      if(k3.eq.4)go to 430
      if(k3.eq.5)go to 531
   60 k3=0
      ms(1)=mas
      do 65 in=2,22
   65 ms(in)=ms(1)-in+1
      do 70 im=1,na
      sumiz=sumiz+cx(im)
   70 crs(im,iz)=crs(im,iz)+cx(im)
      iiz=zoz
      if(test.ne.0.)go to 76
c
      su=0.
      do 3075 jj=1,na
      su=su+crs(jj,iz)
 3075 continue
      if(su.eq.0.)go to 1176
      write(17,480)iiz,iiz,iiz,iiz,iiz,iiz,iiz,iiz,iiz,iiz,iiz
      write(17,485)(ms(in),in=1,na)
      write(17,490)(crs(n,iz),n=1,na)
   76 if(jcal.gt.0)go to 1176
      if(test.ne.0.) go to 90
      write(17,491)(cx(i),i=1,11)
      write(17,505)( fs(ia),ia=1,11)
      write(17,515)(dsp(ia),ia=1,11)
      write(17,510)(brr(ia),ia=1,11)
      write(17,520)(der(ia),ia=1,11)
      write(17,525)( er(ia),ia=1,11)
c
1176    if(na.le.11)go to 90
c
      su=0.
      do 3076 jj=1,na
      su=su+crs(jj,iz)
 3076 continue
      if(su.eq.0.)go to 180
      write(17,480)iiz,iiz,iiz,iiz,iiz,iiz,iiz,iiz,iiz,iiz,iiz
      write(17,485)(ms(in),in=12,na)
      write(17,490)(crs(n,iz),n=12,na)
        if(jcal.gt.0)go to 90
      continue
      write(17,491)(cx(i),i=12,na)
      write(17,505)( fs(ia),ia=12,na)
      write(17,515)(dsp(ia),ia=12,na)
      write(17,510)(brr(ia),ia=12,na)
      write(17,520)(der(ia),ia=12,na)
      write(17,525)( er(ia),ia=12,na)
   90 xmiss=rcss-sumiz
      if(jcal.lt.1)xmiss=xmiss-fiss
c   add neutron width calculation
c
      do 175 n=2,na
      barnun=barnun+crs(n,iz)*float(n-1)
      brnsqn=brnsqn+crs(n,iz)*float(n-1)**2
      barnud=barnud+crs(n,iz)
  175 continue
      if(barnud.le.0.)barnu=0.
      if(barnud.le.0.)brnsq=0.
      if(barnud.le.0.)w=0.
      if(barnud.le.0.)go to 180
      barnu=barnun/barnud
      brnsq=brnsqn/barnud
      wsq=brnsq-barnu**2
      if(wsq.lt.0.)write(17,176)
      if(wsq.lt.0.)go to 180
      w=sqrt(wsq)
  180 continue
  176 format('  neutron width negative')
  177 format('  neutron multiplicity =',e10.3,'  width= ',e10.3)
      if(test.eq.0..and.su.gt.0.)write(17,177)barnu,w
      if(jcal.le.0.and.ifis.ne.5)xmiss=crsum-fiss-sumiz
      if(test.ne.0.) go to 91
      if(su.gt.0.)write(17,177)barnu,w
      profis=fiss/rcss
      if(test.eq.0.)write(17,530)fiss ,sumiz,xmiss
   91 return
   95 nz=nz+1
      na=na+1
      barnun=0.
      barnud=0.
      brnsqn=0.
      if(no.le.1)then
      do 100 iz=1,13
      do 100 ia=1,22
  100 crs(ia,iz)=0.
                 endif
      exc(1,1)=xmax
      do 115 iz=2,nz
  115 exc(iz,1)=exc(iz-1,1)-be(iz-1,1,2)
      nz=nz-1
      nna=na-1
      do 135 iz=1,nz
      do 135 ia=1,nna
  135 exc(iz,ia+1)=exc(iz,ia)-be(iz,ia,1)
      nz=nz+1
      do 150 ia=1,nna
  150 exc(nz,ia)=exc(nz-1,ia)-be(nz-1,ia,2)
      do 155 ndl=3,9,2
      if (nz.lt.ndl.or.na.lt.ndl) go to 255
      do 155 iz=ndl,nz
      do 155 ia=ndl,na
  155 exc(iz,ia)=exc(iz,ia)+28.3
  255 nz=nz-1
      na=na-1
      write (17,305) eq,at,amass,qval
      write(17,265) jcal,jfrac,jupper,jang
      if(iadst.eq.0)write(17,8857)
 8857 format(//35x,' iadst=0,no angular distributions will be calculated'//)
     1'//)
c1000 format(//,' iadst = 1  angular distributions calculated for neutro
c    1ns using nucleon-nucleon kinematics (Phys Rev C 30, 1493 (1984))')
c1001 format(//,' iadst = 2  angular distributions calculated for proton
c    1s using nucleon-nucleon kinematics (Phys Rev C 30, 1493 (1984))')
c1002 format(//,' iadst = 3  angular distributions calculated for neutro
c    1ns using Kalbach systematics (Phys. Rev. C 37, 2350 (1988)')
c     if(iadst.eq.0)goto 8850
c9857  format(//41x,'i3d folding=',i2,5x,' direction cosines for 3d folding'//)
c    1ng'//)
c9858  format(//41x,'i3d= ',i2,5x,'  use two dimensional folding'//)
c8851  format(//41x,' irfr=',i2,5x,'  no refraction'//)
c8852  format(//41x,' irfr=',i2,5x,'  snell law entrance refraction only'//)
c     1//)
c8853  format(//41x,' irfr=',i2,5x,'  snell refraction in,heisenberg out'//)
c    1//)
c8854  format(//41x,' irfr=',i2,5x,'  heisenberg refraction in and out'//)
c    1)
c
c8850  continue
       xmat=xmax/ed
      if(xmat.gt.999.)write(17,260)
      if(xmat.gt.999.)return
  260 format(25x,'excitation energy exceeds dimensioned limit.  calculat
     1ion at this energy terminated.')
      if(jcal.le.0)write(17,325)
      if(jcal.eq.1)write(17,330)
  265 format(/////32x,'  jcal=',i5,'  start j=',i5,' stop j =',i5,'  jan
     1g =',i5/)
      if(jcal.eq.2)write(17,295)
      if(jcal.eq.3)write(17,300)
      if(jang.gt.0)write(17,315)
      if(td.gt.0..and.tmx.eq.0.)write(17,320)
      if(td.gt.0..and.tmx.gt.0.)write(17,335)
c     if(td.gt.0..and.tmx.eq.0.)write(7,320)
c     if(td.gt.0..and.tmx.gt.0.)write(7,335)
  295 format(' s wave approximation with rotating liquid drop moment of
     1inertia at equilibrium deformation ')
  300 format(' s wave approximation with rigid rotor moment of inertia')
  305 format (28x,'laboratory bombardment energy = ',f6.1,'mev',3x,'targ
     1et mass = ',f6.1,//22x,'compound nucleus mass = ',f6.1,2x,'compoun
     2d nucleus formation q value = ',f6.1/)
  315 format(41h option particles remove angular momentum/)
  320 format(41x,' hybrid calculation has been selected '/)
  325 format(41x,' calculation includes fission competition'/)
  330 format(41x,' standard weisskopf-ewing option selected'/)
  335 format(41x,' geometry dependent hybrid model selected'/)
      write (6,236)
  236 format (1h ,1h )
      k3=0
      return
  340 do 345 lz=1,10
      do 345 ia=1,22
      nulim(lz,ia)=1
  345 jamma(lz,ia)=1
      delt = 0.0
      alt = 0.0
      do 360 l=1,jmax
      do 360 k = 1,3
        if(k.eq.1)go to 360
        if(k.eq.3)go to 355
      if (sig(k,l).le.0.0) delt = delt + 1.0
      go to 360
  355 if (sig(k,l).le.0.0) alt =alt + 1.0
  360 continue
      if(k9.eq.1)go to 371
      if(inver.eq.3)go to 371
      do 370 k=1,m3
      do 370 je=jmax,999
  370 sig(k,je)=sig(k,jmax)
  371 continue
        if(m3.lt.3)alt=48./ed
      jmax=999
      pelt=delt*ed-ed/2.
      aat=alt*ed-ed/2.
      if(test.eq.0.)write(17,500)delt,alt,cld
      do 400 lz=1,nz
      do 400 ia=1,na
      if(be(lz,ia,1).ge.(be(lz,ia,2)+pelt))go to 380 
      if(be(lz,ia,1).le.(be(lz,ia,3)+aat))go to 390
  380 if((be(lz,ia,2)+pelt).le.(be(lz,ia,3)+aat))go to 395 
      jamma(lz,ia)=(be(lz,ia,3)+aat)/ed
      go to 400
  390 jamma(lz,ia)=(be(lz,ia,1))/ed
      go to 400
  395 jamma(lz,ia)=(be(lz,ia,2)+pelt)/ed
  400 jamma(lz,ia)=jamma(lz,ia)+1.
      if(limout.eq.0)then
      do 415 lz=1,nz
      do 415 ia=1,na
      if (jamma(lz,ia).gt.0)go to 415 
      write(17,410)lz,ia
  410 format(' cross section for nuclide of iz = ',i2,' and ia = ' ,i2,'
     1 will be meaningless due to negative binding energy')
      jamma(lz,ia)=1
  415 continue
      endif
      if(k9.eq.1)k3=0
      if(k9.eq.1)return
      if(inver.eq.3)return
c c(x) constants 1=n,2=p,3=alpha,4=deut,5=triton,6=3He,7=7Be 
c the c() multiplier is (2s+1)*reduced mass of ejectile, divided
c by two for historical reasons.
c------------------------------
      c(1)=(amass-1.)/amass
      c(2)=c(1)
      c(3)=2. *(amass-4.)/amass
      c(4)=3.*(amass-2.)/amass
      c(5)=3.*(amass-3.)/amass
      c(6)=c(5)
      c(7)=14.*(amass-7.)/amass
c-----------------------------------
c re-store inverse cross sections as sig.inv*(2s+1)*reduced mass*energy
      do 425 l =1,999
      bl=float(l)*ed-ed/2.
      do 425 k=1,m3
  425 sig(k,l)=c(k)*sig(k,l)*bl*ed
      k3=0
      return
 430  continue
      if(ike.eq.4.and.ia.eq.na.and.iz.eq.nz)write(17,435)
      if(ike.eq.4.and.(ia.ne.na.or.iz.ne.nz)) go to 461
c
  435 format(1h /' particle spectra summed over all emitting ',
     1 'nuclei (and partial waves if applicable) with precompound ',
     2 'included.')
      write(17,440)aoz,zoz
  440 format(1h /' kinetic energy spectra from a = ',f5.1,'  z = ',
     1 f5.1,' (note that actual channel energy ',
     2 'equals ke index )')
      write(17,455)
      write(7,455)
      ssm1=0.
      ssm2=0.
      ssm3=0.
      ssm4=0.
      vr=sqrt(2.*eq*ap/((ap+at)**2))
      con=2.*at/(ap+at)
      do 450 kee=1,999
      ssm1=ssm1+en(1,kee)*ed
      ssm2=ssm2+en(2,kee)*ed
      ssm3=ssm3+en(3,kee)*ed
c     ssm4=ssm4+en(4,kee)*ed
      ake=float(kee)*ed-ed/2.
c 2/20/2001 put in approx lab energy values for zero&90 deg. emission
      cms=en(1,kee)+en(2,kee)+en(3,kee)
      if(cms.eq.0.)go to 450
      el0=0.5*(vr+sqrt(ake*con))**2
      el90=ake*con/2.
      write(7,460)ake,en(1,kee),en(2,kee),en(3,kee)
      write(17,860)ake,el0,el90,en(1,kee),en(2,kee),en(3,kee)
  450 continue
      write(17,4601)ssm1,ssm2,ssm3,ssm4
 4601 format('sum cross sections =     ',1e11.4,'  ',1e11.4,'  ',1e11.4
     1,'  ',1e11.4)
c     if(ike.ne.4.or.ia.ne.na.or.iz.ne.nz) go to 461
c 456 continue
      do 4513 i=1,2
      nmx=xmax/ed
      if(nmx.gt.999)nmx=999
      write(17,4598)ncount
4598  format(' ncount= ',i9)
      if(iadst.eq.0.and.ncount.eq.0)go to 451
      do 4511 ke=1,nmx
      if(ike.eq.4)then
      do4510 it=1,36
4510  dsgp(it,ke,i)=dsgp(it,ke,i)+(en(i,ke)-scrs(ke,i))/12.566
      else
      do4512 itt=1,36
4512  dsgp(itt,ke,i)=dsgp(itt,ke,i)+en(i,ke)/12.566
      endif
c
c
4511  continue
      confac=3.14*5./180.
      if(iadst.eq.4)write(17,1003)
      write(17,406)
      if(i.eq.1)write(17,8417)
      if(i.eq.2)write(17,8418)
      write(17,6562)
6562  format('   ')
8417  format(' for neutrons ')
8418  format(' for protons ')
  406 format('   total angular distribution  compound plus precompound')
      nend=xmax/(10.*dlt)+0.99
c
      do 519 j=1,nend
      bji(1)=10.*ed*float(j-1)+ed/2.
      do 8406 ii=2,10
8406  bji(ii)=bji(ii-1)+ed
c
      write(17,401)(bji(ii),ii=1,10)
      nii=j*10-9
      nif=nii+9
c
      do 6510 id=1,10
6510  asum(id)=0.
c
c
      do 509 ith=1,36
c
      iii=0
      do 6511 kef=nii,nif
      iii=iii+1
6511  asum(iii)=asum(iii)+sign(ith)*dsgp(ith,kef,i)*confac
c
      angle=5.*ith-2.5
      as=asum(iii)
      if(as.gt.0.)write(17,402)angle,(dsgp(ith,kef,i),kef=nii,nif)
  509  continue
c
      write(17,6562)
      write(17,6512)(asum(ii),ii=1,10)
      write(17,6562)
6512  format(' integral=',10(e10.3))
c
  519 continue
  451 continue
 1003 format(//' iadst = 4  angular distributions calculated for proton
     1s using Kalbach systematics (Phys. Rev. C 37, 2350 (1988)')
4513  continue
4620  format(' sum neutron=',e10.3,'   sum proton =',e10.3,
     1'  sum alpha=',e10.3,'  sum deut=',e10.3)
  455 format('  ke  ds/de(mb/mev)=        NEUTRONS    PROTONS     ALPHAS',
     1      deuterons')
  460 format(' ',f8.3,'                ',e11.4,'  ',e11.4,'  ',e11.4,' '
     1,e10.3)
  860 format(' ',3(1x,f8.3),'                ',e11.4,'  ',e11.4,'  '
     1,e11.4,' ',e10.3)
  461 k3=0
      son=0.
      sop=0.
      soa=0.
      sod=0.
      do 457 kee=1,999
      son=son+en(1,kee)
      sop=sop+en(2,kee)
      soa=soa+en(3,kee)
457   continue
      son=son*ed
      sop=sop*ed
      soa=soa*ed
      sod=sod*ed
      write(17,4620)son,sop,soa
      return
  531  continue
      iikl=iikl-1
      write(17,550)
  550 format('  isotopically weighted results ')
2551  format(1e10.3,10x,11(1e10.3))
      write(17,2553)
2553  format('  reaction cross sections for different isotopic targets')
      write(17,2550)(mtar(ksor),ksor=1,ksot)
      write(17,2552)(iete(lsot),lsot=1,ksot)
2552  format('  MeV      At.No.=',10(2x,i5,3x))
      do 2554 kij=1,iikl
      write(17,2551)engbuf(kij),(rcsi(lsot,kij),lsot=1,ksot)
2554  continue
2550  format('  E-lab    Mass no.',10(1x,i5,4x))
      nz=knz
      ed=exd
      do 543 lz=1,nz
      ilz=zee+1.-float(lz)
      do 549 kn=1,22
      ms(kn)=amass+2.-float(kn)-lz
549   continue
      write(17,480)ilz,ilz,ilz,ilz,ilz,ilz,ilz,ilz,ilz,ilz,ilz
      write(17,485)(ms(in),in=1,11)
      write(17,4855)
4855  format('  e(lab)  ')
      do 5421 kij=1,iikl
      write(17,990) engbuf(kij),(ccrs(n,lz,kij),n=1,11)
5421  continue
      lpa=ia
      if(lpa.le.11)go to 543
      write(17,480)ilz,ilz,ilz,ilz,ilz,ilz,ilz,ilz,ilz,ilz,ilz
      write(17,485)(ms(in),in=12,22)
      write(17,4855)
      do 5431 kij=1,iikl
      write(17,990)engbuf(kij),(ccrs(n,lz,kij),n=12,22)
5431  continue
  543 continue
      if(ike.eq.0)return
      do 650 kij=1,iikl
      write(17,5501)engbuf(kij)
5501  format('  isotopically weighted results ',
     1f10.3,' mev(lab)')
      write(17,455)
      write(7,455)
      do 650 kee=1,999
      ake=float(kee)*ed-ed/2.
      cms=een(1,kee,kij)+een(2,kee,kij)+een(3,kee,kij)
      if(cms.eq.0.)go to 650
      write(17,460)ake,een(1,kee,kij),een(2,kee,kij),een(3,kee,kij)
  650 continue
c 541 continue
      k3=0
      return
  401 format(' angle/deg. ke=',f6.2,9(4x,f6.2))
cneadb402 format(x,f5.1,5x,10(e10.3))
  402 format(1x,f5.1,5x,10(e10.3))
  480 format(15h0atomic number ,i3,10(7x,i3))
  485 format(15h mass number   ,i3,10(7x,i3))
  990 format(12(1x,e9.3))
  490 format(' er xsect s(j) ',1e10.3,10(1e10.3))
  491 format(15h er xsect   j  ,1e10.3,10(   e10.3))
  500 format (/34x,'delt = ',f10.1,6x,'alt = ',f10.1,6x,'cld =',f10.5/)
  505 format(15h fiss x sect   ,11e10.3)
  510 format(15h fiss barrier  ,11e10.3)
  515 format(15h del saddle pt ,11e10.3)
  520 format(15h del rotation  ,11e10.3)
  525 format(15h rotational en ,11e10.3/)
  530 format(' fission cross section to here=',e10.3,2x,' er xsect to he
     1re  =',e10.3,2x,' missing cross section  =',e10.3)
      end
      subroutine angout(i)
      common/outlim/limout,lim2
      common/distb/dsgp(36,999,8)
      common/renorm/ssm1,ssm2,ssm3,ssm4
      common/shft/mc,mp,inver,ike,ipch,kq,qval,ap,at,zp,zt,cld,barfac
      common/incr/ed,tened,efinal
      common/b2/nmx,iqx,b(4)
      common/sf6/aoz,zoz,en(8,999)
      common/sft5/exc(18,27),xmxx
      common/lym1/be(18,27,8),pair(18,27),xmas(18,27),delshl(18,27),amas
     1(18,27)
c     common/bym1/symb(18,27),symbp(18,27),symbs(18,27)
      common/lbcm2/indx(999)
      common/labcm/encm(999,36),thetcm(999,36),xjac(999,36)
     1,enlab(999),thetl(36),rcoil(999,36),thel(36),si(36)
      common/dist/n,ni,ke,irefr,xz(3),pod,pcs,ek,xmix,cfrac
     1,dtau,sign(36),tmix,bji(10),fin,dtup
      common/memo/scrs(999,8)
      common/par3/ep,sigml(999),acrs(999)
      common/sums/sun(999)
      common/labcm2/ercoil(1100)
      common/cslabor/cslab(4,999)
 
      if(i.eq.1)write(17,5)
      if(i.eq.2)write(17,6)
      if(i.eq.3)write(17,600)
      if(i.eq.4)write(17,601)
 600  format('  laboratory angular distributions for alphas using
     1kalbach sytematics,evap only treated')
 601  format('  laboratory angular distributions for deuterons using
     1kalbach systematics,evap plus precompound treated')
   5  format('  laboratory angular distributions for neutrons using
     1kalbach systematics  ')
   6  format('   laboratory angular distributions for protons using',
     1'kalbach systematics')
      do 7 ii=1,999
      sun(ii)=0.
   7  continue
      confac=ed*3.14159*5.0/180.
      sumnp=0.
      ian=float(nmx)/10.+0.99
      if(ian.gt.30)ian=30
      if(ian.eq.0)ian=1
      sumk=0.
      sumf=0.
      do 2 j=1,36
      dg=5.*float(j)-2.5
      arg=dg*3.14159/180.
    2 sign(j)=2.*3.14159*sin(arg)
      do 50 ie=1,nmx
      do 49 ith=1,36
      sun(ie)=sun(ie)+dsgp(ith,ie,i)*sign(ith)*confac
   49 continue
      sumf=sumf+sun(ie)
   50 continue
      do 710 ie=1,nmx
  710 sun(ie)=0.
      if(sumf.eq.0.)sumf=0.1
      if(i.eq.1)fact=ssm1/sumf
      if(i.eq.2)fact=ssm2/sumf
      if(i.eq.3)fact=ssm3/sumf
      if(i.eq.4)fact=ssm4/sumf
      write(17,5001)fact
 5001 format('  renormalization factor for ddcs = ',f8.3)
      do 700 ie=1,nmx
      do 690 ith=1,36
      dsgp(ith,ie,i)=dsgp(ith,ie,i)*fact
      sun(ie)=sun(ie)+dsgp(ith,ie,i)*sign(ith)*confac
  690 continue
      engy=float(ie)*ed-ed/2.
      sumk=sun(ie)*engy+sumk
      if(sun(ie).eq.0.)sun(ie)=1.*10.**(-20)
      cslab(i,ie)=sun(ie)
      sumnp=sumnp+sun(ie)
  700 continue
      if(i.eq.1)write(7,92)ep
      if(i.eq.2)write(7,93)ep
      if(i.eq.3)write(7,95)ep
      if(i.eq.4)write(7,96)ep
      eavg=sumk/sumnp
      gdkarm=0.001*sumnp*eavg*9.65/at
      if(i.eq.1)write(17,960)eavg,gdkarm,sumnp
      if(i.eq.2)write(17,961)eavg,gdkarm,sumnp
      if(i.eq.3)write(17,962)eavg,gdkarm,sumnp
      if(i.eq.4)write(17,963)eavg,gdkarm,sumnp
      if(i.eq.1)write(7,960)eavg,gdkarm,sumnp
      if(i.eq.2)write(7,961)eavg,gdkarm,sumnp
      if(i.eq.3)write(7,962)eavg,gdkarm,sumnp
      if(i.eq.4)write(7,963)eavg,gdkarm,sumnp
 960  format(' average neutron energy =',f7.2,'mev karma =',e10.3,
     1'*10**-15 grays*m**2,total neutron cross section=',e10.3,'mb')
 961  format(' average proton  energy =',f7.2,'mev karma =',e10.3,
     1'*10**-15 grays*m**2,total proton emission cross section=',e10.3,'mb')
     2mb')
 962  format(' average alpha   energy =',f7.2,'mev karma =',e10.3,
     1'*10**-15 grays*m**2,total alpha emission cross section=',e10.3,'mb')
     2b')
 963  format(' average deuteron energy =',f7.2,'mev karma =',e10.3,
     1'*10**-15 grays*m**2,total deuteron emission cross section=',e10.3
     2,'mb')
   95 format('alpha lab energy,sdcs for incident energy',f6.1)
   96 format('deuteron lab energy,sdcs for incident energy',f6.1)
   92 format('neutron lab energy,sdcs for incident energy',f6.1)
   93 format('proton lab energy,sdcs for incident energy ',f6.1)
      do 100 k=1,ian
      j1=10*(k-1)+1
      j2=j1+9
      sumen=0.
      do 790 is=j1,j2
      sumen=sumen+enlab(is)
790   continue
      if(sumen.gt.0.)write(17,8)(enlab(j),j=j1,j2)
      write(17,9)
      write(17,10)
      do 90 ith=1,36
      sumth=0.
      do 791 it=j1,j2
      sumth=sumth+dsgp(ith,it,i)
 791  continue
      if(sumth.gt.0.)write(17,11)thel(ith),(dsgp(ith,j,i),j=j1,j2)
   90 continue
      write(17,12)(sun(j),j=j1,j2)
      do 94 j=j1,j2
      sume=0.
      do 694 ii=j1,j2
      sume=sume+enlab(ii)
 694  continue
 
      if(sume.gt.0.)write(7,91)enlab(j),sun(j)
   94 continue
   91 format(2e10.3)
  100 continue
      if(limout.eq.0)then
      if(i.eq.1)write(17,13)sumnp
      if(i.eq.2)write(17,14)sumnp
      if(i.eq.3)write(17,603)sumnp
      if(i.eq.4)write(17,604)sumnp
      if(i.eq.2)write(17,500)
      if(i.eq.2)write(17,503)
      endif
      do 505 jj=1,50
      xe=(float(jj)-1.)*.2
      j1=10*jj-9
c     if(i.eq.2) write(17,504)xe,(ercoil(ll),ll=j1,j1+9)
  505 continue
  503 format('       .010      .030    .050    .070    .090    .11')
  500 format('recoil spectra follow in 20 kev increments')
   13 format('total neutron emission cross section = ',1e10.3)
   14 format('total proton emission cross section = ',1e10.3)
 603  format('total alpha emission cross section = ',1e10.3)
 604  format('total deuteron emission cross section =',1e10.3)
    8 format(' elab(mev) ',10(4x,f7.2))
    9 format('  angle')
   10 format(' deg lab')
   11 format(2x,f6.2,3x,10(1x,e10.3))
   12 format(' sum over angle',10(1x,e10.3))
      return
      end
c
      subroutine reclsv
      common/outlim/limout,lim2
      common/labcon/prx,pry,prz,pcm,rmas,clab(999,8),dslb(36,999,8)
      common/scat/jz,ja, probxx(2,2),q1(8,1100),fcs,eloss,fcspec
      common/pl3/jx,nepr,eb(100),rcsp(100),zee,amass,plex
      common/incr/ed,tened,efinal
      common/angprm/aprp(500),eav3(1000),
     1eav2(1000),pot2x,thet1,phi1,snt(360),cst(360),snp(720),csp(720)
     2,pot
      common/rec1/ermxc,thetrc
      common/recoils/rcletza(100,13,24,18),rcleza(100,13,24),
     1rclea(100,36),rcleat(100,36,18),rphir(36)
     2,rclsp(200),rclz(200),rcly(200),rclx(200),escale
     3,rcle(200),rclt(36),rclet(200,36)
c
c
        return
       arecl=amass+2.-float(ja)-float(jz)
       arec=2.*arecl
       renx=prx**2/arec
       reny=pry**2/arec
       renz=prz**2/arec
       ercl=renx+reny+renz
c      prsq=prx**2+pry**2+prz**2
c      ercl=prsq/(2.*arecl)
c calculate recoil energies total,x,y,z
       iex=renx*escale+1.
       iey=reny*escale+1.
       iez=renz*escale+1.
       ier=ercl*escale+1.
       if(iex.gt.200)iex=200
       if(iey.gt.200)iey=200
       if(iez.gt.200)iez=200
       if(ier.gt.200)ier=200
       rclsp(ier)=rclsp(ier)+fcs
       rclz(iez)=rclz(iez)+fcs
       rclx(iex)=rclx(iex)+fcs
       rcly(iey)=rcly(iey)+fcs
c      iercl=100.*ercl +1.
c      if(iercl.gt.999)iercl=999
c      rclsp(iercl)=rclsp(iercl)+fcs
c     vx=prx/arecl
c     vy=pry/arecl
c     vz=prz/arecl
c     vrec=sqrt(vx**2+vy**2+vz**2)
      prxy=sqrt(prx**2+pry**2)
      thetar=atan2(prxy,prz)
c     if(prx.eq.0.)prx=0.0001
c     arg=(pry/prx)
c     phir=atan(arg)
      phir=atan2(pry,prx)
      if(phir.lt.0.)phir=phir+6.2832
c     write(14,10020)phir,prx,pry,prz
c10020 format('phir=',f10.2,'prx=',f10.3,'pry=',f10.3,'prz=',f10.3)
      if(prx.lt.0..and.pry.lt.0.)phir=3.14159+phir
      if(prx.gt.0..and.pry.lt.0.)phir=3.14159+phir
c     write(14,10721)phir,prx,pry,prz
c     ipr=1.+phir*18./3.1416
c     phb(ipr)=phb(ipr)+1.
c10721 format(' second phir=',f10.2,'prx=',f10.3,'pry=',f10.3,'prz=',f10.
c     13)
c now have energy and angles of recoils
c make 4 flavors of recoil arrays
c first define indices
       ithr=1.+thetar*thetrc
       it=1.+thetar*thetrc*2.
c     if(limout.eq.0)then
c      if(ithr.gt.18)write(14,9001)thetar,thetrc,ithr
c      if(it.gt.36)write(14,9001)thetar,thetrc,ithr
c      if(ithr.lt.0.)write(14,9001)thetar,thetrc,ithr
c      endif
c9001   format('thetar,thetrc,ithr',2f10.3,i5)
       ier=1.+ercl*ermxc
       if(ier.gt.100)ier=100
c      if(phir.lt.0.)phir=6.2831+phir
c     write(14,10022)phir,prx,pry,prz
c10022 format('   third phir=',f10.2,'prx=',f10.3,'pry=',f10.3,'prz=',f10
c     1.3)
c      iphii=1.+phir*thetrc
       iia=ja+jz-1
       if(iia.gt.24)iia=24
c     if(limout.eq.0)then
c      if(iia.le.1)write(14,9000)ja,jz,ercl
c      endif
c9000   format('ja,jz,ercl',2i5,f6.2)
c above define temp phi index to check isotropicity in rphir array
       if(jz.gt.13)jz=13
       item=0
c put in item to stop warning
       if(item.eq.0)go to 10000
c      go to 10000
c above temporar
       rcletza(ier,jz,ja,ithr)=rcletza(ier,jz,ja,ithr)+fcs
       rcleza(ier,jz,ja)=rcleza(ier,jz,ja)+fcs
       rcleat(ier,iia,ithr)=rcleat(ier,iia,ithr)+fcs
       rclea(ier,iia)=rclea(ier,iia)+fcs
c      rphir(iphii)=rphir(iphii)+fcs
       rclet(ier,ithr)=rclet(ier,ithr)+fcs
       rcle(ier)=rcle(ier)+fcs
       rclt(ithr)=rclt(ithr)+fcs
10000  continue
      return
      end
      subroutine bngwrt
      common/endvarray/endvcm1(36,999,8),endvcm2(6,36,999,8),
     1endvcm3(10,36,999,8),cmen1(8,999),cmen2(6,8,999),cmen3(10,8,999)
      common/endvparm/iend,lll(6,6,6),ll(6,6)
      common/tgspin/levno
      common/ist/isot,iikl,knz,idelm,fracts,fract,ccrs(32,13,100)
     1,engbuf(100),exd,een(8,999,100)
        common/mccom/ncount,ipot,rum(15,24),ppmc(15,24,999),cf
      common/memof/sfcrs(999,8),dfsgp(36,999,4),efn(4,999)
      common/fisn/dfslb(36,999,4),cflab(999,4),dfetl(36,999,4)
      common/outlim/limout,lim2
      common/sft5/exc(18,27),xmax
      common/dist/n,ni,ke,irefr,xz(3),pod,pcs,ek,xmix,cfrac
     1,dtau,sign(36),tmix,bji(10),fin,dtup
      common/labcon/prx,pry,prz,pcm,rmas,clab(999,8),dslb(36,999,8)
      common/unit8/angles(9)
      common/mbcfudg/cfg(999)
      common/shft/mc,mp,inver,ike,ipch,kq,qval,ap,at,zp,zt,cld,barfac
      common/pl3/jx,nepr,eb(100),rcsp(100),zee,amass,plex
      common/incr/ed,tened,efinal
      common/send/irfr,iadst,dlt
      common/eder/ed2,v,vhalf,vsq
      common/distb/dsgp(36,999,8)
      common/out1/nout
      common/fis5/afmm(36),afm(36)
      common/out4/asm(36),asmm(36)
      common/exclus/ gib(8,999),pr(3),ge(2100),zz(3),
     1delta(3),ts(3),pairx(8),asum(10),rate(8,999)
      common/bill/sumethl(999,8),sumthel(36,8),sumethc(999,8),sumthec(36,8)
     1,8)
      common/bilg/sufethl(999,8),sufthel(36,8),sufethc(999,8),sufthec(36,8)
     1,3)
      confac=3.14*5./180.
      if(limout.eq.0)then
       write(42,*)zp,ap,zt,at,isot,levno,' ZP,AP,ZT,AT,ISOT,LEVNO'
       write(40,*)zp,ap,zt,at,isot,levno,' ZP,AP,ZT,AT,ISOT,LEVNO'
       write(11,*)zp,ap,zt,at,isot,levno,' ZP,AP,ZT,AT,ISOT,LEVNO'
       write(18,*)zp,ap,zt,at,isot,levno,' ZP,AP,ZT,AT,ISOT,LEVNO'
       write(43,*)zp,ap,zt,at,isot,levno,' ZP,AP,ZT,AT,ISOT,LEVNO'
       write(44,*)zp,ap,zt,at,isot,levno,' ZP,AP,ZT,AT,ISOT,LEVNO'
       write(45,*)zp,ap,zt,at,isot,levno,' ZP,AP,ZT,AT,ISOT,LEVNO'
       write(55,*)zp,ap,zt,at,isot,levno,' ZP,AP,ZT,AT,ISOT,LEVNO'
       write(54,*)zp,ap,zt,at,isot,levno,' ZP,AP,ZT,AT,ISOT,LEVNO'
      write(11,6562)
6562  format('   ')
      write(11,406)
      write(44,*)' Total SDCS fission plus evaporation'
      write(43,*)' Total DDCS fission plus evaporation'
  406 format('    angular distribution   precompound')
       call headers(7)
      endif
c
      write(17,406)
       do 7000 kk=1,3
        do 7000 i=1,999
      sumethl(i,kk)=0.
      sumethc(i,kk)=0.
      sufethl(i,kk)=0.
      sufethc(i,kk)=0.
 7000 continue
       do 7001 kk=1,3
        do 7001 i=1,36
      sumthec(i,kk)=0.
      sumthel(i,kk)=0.
      sufthec(i,kk)=0.
      sufthel(i,kk)=0.
 7001 continue
      dlt=ed
      nend=xmax/(10.*dlt)
      nen=(xmax+ed)/ed
      do 7003 j=1,3
      do 7003 ith=1,36
      do 7003 ie=1,nen
      sumethc(ie,j)=sumethc(ie,j)+dsgp(ith,ie,j)
      sumethl(ie,j)=sumethl(ie,j)+dslb(ith,ie,j)
      sumthec(ith,j)=sumthec(ith,j)+dsgp(ith,ie,j)
      sumthel(ith,j)=sumthel(ith,j)+dslb(ith,ie,j)
      sufethc(ie,j)=sufethc(ie,j)+dfsgp(ith,ie,j)
      sufethl(ie,j)=sufethl(ie,j)+dfslb(ith,ie,j)
      sufthec(ith,j)=sufthec(ith,j)+dfsgp(ith,ie,j)
      sufthel(ith,j)=sufthel(ith,j)+dfslb(ith,ie,j)
 7003 continue
c
      do 519 i=1,3
c
       write(11,*)
       write(58,*)
       write(17,*)
       write(40,*)
       write(42,*)
       write(43,*)
      if(i.eq.1)then
      write(17,99885)
      write(11,99885)
      write(40,99885)
      write(42,99885)
      write(43,99885)
c     write(44,99885)
      endif
      if(i.eq.2)then
      write(17,99886)
      write(11,99886)
      write(40,99886)
      write(42,99886)
      write(43,99886)
c     write(44,99886)
       write(11,*)
       write(17,*)
       write(40,*)
       write(42,*)
       write(43,*)
      endif
      if(i.eq.3)then
      write(17,49886)
      write(11,49886)
      write(40,49886)
      write(42,49886)
      write(43,49886)
      endif
      write(18,6562)
      write(11,*)
      write(40,*)
      write(17,*)' table entries are in mb/MeV*sr'
      write(18,*)' table entries are in mb/MeV*sr'
      write(11,*)' table entries are in mb/MeV*sr'
      write(11,*)
      write(42,*)' table entries are in mb/MeV*sr'
      write(40,*)
      write(40,*)' table entries are in mb/MeV*sr'
      write(43,*)' table entries are in mb/MeV*sr'
      write(40,*)
      write(18,*)
      write(11,*)' integrals at base of columns are mb/MeV'
      write(40,*)' integrals at base of columns are mb/MeV'
      write(42,*)' integrals at base of columns are mb/MeV'
      write(43,*)' integrals at base of columns are mb/MeV'
      write(18,*)' integrals at base of columns are mb/MeV'
      write(11,*)
      write(40,*)
      write(18,*)
crc
99885 format(' monte carlo NEUTRON angular distribution')
99886 format(' monte carlo PROTON angular distribution')
49886 format(' monte carlo ALPHA angular distribution')
  401 format('    angle    channel energy (MeV)')
c8401  format(' deg.',2x,f7.2,3x,9(1x,f8.2))
8401  format('    deg.',2x,10(1x,f8.2))
c 401 format(' angle/deg. ke=',f7.2,9(2x,f8.2))
c
      write(17,401)
       if(limout.gt.0)go to 9999
      do 519 j=1,nend
c
      bji(1)=10.*dlt*float(j-1)+ed/2.
      do 8406 ii=2,10
8406  bji(ii)=bji(ii-1)+dlt
c     write(17,401)
c     if(limout.lt.1)then
      write(18,401)
      write(11,4010)
      write(42,4010)
      write(40,401)
      write(43,4010)
4010  format('   angle,deg.  laboratory energy,MeV')
      write(18,8401)(bji(ii),ii=1,10)
      write(11,8401)(bji(ii),ii=1,10)
      write(42,8401)(bji(ii),ii=1,10)
      write(40,8401)(bji(ii),ii=1,10)
      write(43,8401)(bji(ii),ii=1,10)
c     endif
      write(17,8401)(bji(ii),ii=1,10)
      nii=j*10-9
      nif=nii+9
c
      do 6510 id=1,36
crc
      asm(id)=0.
      afm(id)=0.
crc
      afmm(id)=0.
6510  asmm(id)=0.
c
c
      do 509 ith=1,36
c
      iii=ith
c     iii=0
      do 6511 kef=nii,nif
c     iii=iii+1
      asmm(iii)=asmm(iii)+dsgp(ith,kef,i)
      afmm(iii)=afmm(iii)+dfsgp(ith,kef,i)
crc
      asm(iii)=asm(iii)+dslb(ith,kef,i)
      afm(iii)=afm(iii)+dfslb(ith,kef,i)
      dslb(ith,kef,i)=dslb(ith,kef,i)/(sign(ith)*confac)
      dfslb(ith,kef,i)=dfslb(ith,kef,i)/(sign(ith)*confac)
      dfetl(ith,kef,i)=dfslb(ith,kef,i)+dslb(ith,kef,i)
c above add sum of evap+fission emittees in lab frame
crc
      dfsgp(ith,kef,i)=dfsgp(ith,kef,i)/(sign(ith)*confac)
6511  dsgp(ith,kef,i)=dsgp(ith,kef,i)/(sign(ith)*confac)
c check the constants here!
c
      angle=5.*ith-2.5
      as=asmm(iii)
      af=afmm(iii)
      if(as.gt.0.)write(18,402)angle,(dsgp(ith,kef,i),kef=nii,nif)
      if(af.gt.0.)write(40,402)angle,(dfsgp(ith,kef,i),kef=nii,nif)
crc
       tq=asm(iii)
       tf=afm(iii)
       if(tf.gt.0.)write(42,402)angle,(dfslb(ith,kef,i),kef=nii,nif)
       if(tq.gt.0.)write(11,402)angle,(dslb(ith,kef,i),kef=nii,nif)
      if(tq.gt.0.)write(43,402)angle,(dfetl(ith,kef,i),kef=nii,nif)
crc
  509  continue
c
      write(18,*)
      write(11,*)
      write(18,59503)(sumethc(kef,i),kef=nii,nif)
      write(42,59503)(sufethl(kef,i),kef=nii,nif)
      write(40,59503)(sufethc(kef,i),kef=nii,nif)
      write(11,59503)(sumethl(kef,i),kef=nii,nif)
      write(11,6562)
      write(18,6562)
59503 format(' integral ',10(1pe9.2))
crc
crc
  402 format(1x,f8.2,1x,10(1pe9.2))
c6512  format(' integral= ',10(1pe9.2))
c
  519 continue
c add file 8 output for energy distributions
      jeng=(xmax+ed)/ed+1.
      do 99499 i=1,3
      write(8,*)
      write(40,*)
      write(42,*)
      write(12,*)
      if(i.eq.1)write(8,99498)
      if(i.eq.2)write(8,99497)
      if(i.eq.3)write(8,49497)
      if(i.eq.1)write(40,99498)
      if(i.eq.2)write(40,99497)
      if(i.eq.3)write(40,49497)
crc
      if(i.eq.1)write(12,99498)
      if(i.eq.2)write(12,99497)
      if(i.eq.3)write(12,49497)
      if(i.eq.1)write(42,99498)
      if(i.eq.2)write(42,99497)
      if(i.eq.3)write(42,49497)
      write(8,*)
      write(40,*)
      write(42,*)
      write(12,*)
 
      write(40,99506)
      write(42,99506)
      write(42,*)' table entries are in mb/MeV*sr'
      write(40,*)' table entries are in mb/MeV*sr'
      write(40,*)
      write(40,*)' the integrals at bottom of columns are the sums of all
     1l cross sections multiplied by '
      write(40,*)' 2pi*sin(theta)*(5./180.)*pi,as bins are 5 degrees wid
     1e'
      write(40,*)
      write(42,*)' the integrals at bottom of columns are the sums of al
     1l cross sections multiplied by '
      write(42,*)' 2pi*sin(theta)*(5./180.)*pi,as bins are 5 degrees wid
     1e'
       write(12,*)zp,ap,zt,at,isot,levno,' ZP,AP,ZT,AT,ISOT,LEVNO'
       write(8,*)zp,ap,zt,at,isot,levno,' ZP,AP,ZT,AT,ISOT,LEVNO'
      write(8,99506)
      write(12,99506)
      write(12,*)' table entries are in mb/MeV*sr'
      write(8,*)' table entries are in mb/MeV*sr'
      write(8,*)
      write(8,*)' the integrals at bottom of columns are the sums of all
     1cross sections multiplied by '
      write(8,*)' 2pi*sin(theta)*(5./180.)*pi,as bins are 5 degrees wide'
     1'
      write(8,*)
      write(12,*)' the integrals at bottom of columns are the sums of all
     1l cross sections multiplied by '
      write(12,*)' 2pi*sin(theta)*(5./180.)*pi,as bins are 5 degrees wide'
     1e'
      write(12,*)
crc
99498  format(' NEUTRON out energy distributions,energy vs. angle deg')
99497 format('PROTON out energy distributions,energy vs. angle deg')
49497 format('ALPHA out energy distributions,energy vs. angle deg')
      if(jeng.gt.999)jeng=999
      do 99500 index=1,4
      ji=1+(index-1)*9
      jf=ji+8
      do 99501 ian=1,9
      angles(ian)=(float(ian)*5.-2.5)+45.*(float(index)-1.)
99501 continue
      write(40,99506)
      write(42,99506)
      write(8,99506)
      write(12,99506)
99506 format('   ')
      write(40,99505)(angles(ian),ian=1,9)
      write(42,99505)(angles(ian),ian=1,9)
      write(8,99505)(angles(ian),ian=1,9)
      write(12,99505)(angles(ian),ian=1,9)
99505 format(' angle: ',9(3x,f7.2))
99503 format(1x,f7.2,1x,9(1x,1pe9.2))
      write(40,99504)
      write(42,99504)
      write(8,99504)
      write(12,99504)
99504 format(' energy ')
      do 99502 ienn=1,jeng
      engt=float(ienn)*ed-ed/2.
      st=0.
      sf=0.
      do 89504 jn=ji,jf
      st=st+dsgp(jn,ienn,i)
      sf=sf+dfsgp(jn,ienn,i)
89504 continue
      if(st.gt.0.)write(8,99503)engt,(dsgp(jan,ienn,i),jan=ji,jf)
      if(sf.gt.0.)write(40,99503)engt,(dfsgp(jan,ienn,i),jan=ji,jf)
crc
      st=0.
      sf=0.
      do 89505 jn=ji,jf
      st=st+dslb(jn,ienn,i)
      sf=sf+dfslb(jn,ienn,i)
89505 continue
      if(st.gt.0.)write(12,99503)engt,(dslb(jan,ienn,i),jan=ji,jf)
      if(sf.gt.0.)write(42,99503)engt,(dfslb(jan,ienn,i),jan=ji,jf)
crc
99502 continue
      write(8,*)
      write(12,*)
      write(40,*)
      write(42,*)
      write(8,69504)(sumthec(ith,i),ith=ji,jf)
      write(12,69504)(sumthel(ith,i),ith=ji,jf)
      write(40,69504)(sumthec(ith,i),ith=ji,jf)
      write(42,69504)(sumthel(ith,i),ith=ji,jf)
69504 format(' integral',9(1x,1pe9.2))
99500 continue
99499 continue
c4620  format(' sum neutron=',e10.3,'   sum proton =',e10.3,
c     1'  sum alpha=',e10.3,'  sum deut=',e10.3)
c 455 format('  ke  ds/de(mb/mev)=        NEUTRONS    PROTONS     ALPHAS',
c    1      deuterons')
c 460 format(' ',f8.3,'                ',e11.4,'  ',e11.4,'  ',e11.4,' '
c    1,e10.3)
9999  continue
      return
      end
c-------------------------------------------------------------------------
       subroutine endvout
      common/endvparm/iend,lll(6,6,6),ll(6,6)
      common/nuopt/kinem,iem,kb
      common/endvarray/endvcm1(36,999,8),endvcm2(6,36,999,8),
     1endvcm3(10,36,999,8),cmen1(8,999),cmen2(6,8,999),cmen3(10,8,999)
      common/unit8/angles(9)
      common/incr/ed,tened,efinal
      common/eder/ed2,v,vhalf,vsq
      common/endvo/sendv1(8),sendv2(6,8),sendv3(10,8),sumthec(36,8),
     1sumend1(999),sumend2(6,999),sumend3(10,999)
        if(iend.eq.1)write(58,*)' Results are in CM coordinates, channe
     1l energies'
      if(iend.eq.2)write(58,*)' Results are in LAB coordinates'
      if(kinem.gt.0)then
       do 7001 kk=1,3
        do 7001 i=1,36
      sumthec(i,kk)=0.
 7001 continue
      do 7003 j=1,3
      do 7003 ith=1,36
      do 7003 ie=1,999
      sumthec(ith,j)=sumthec(ith,j)+endvcm1(ith,ie,j)
c for one and only one particle out, sum each type over energy each angle
 7003 continue
       endif
      write(58,*)'  ENDF exclusive Spectra follow:'
      write(58,*)
       do k=1,3
      sendv1(k)=0.
           enddo
       do ie=1,999
       sumend1(ie)=0.
           enddo
       do k=1,3
         do ie=1,999
      sendv1(k)=sendv1(k)+cmen1(k,ie)
          enddo
           enddo
c get total emission cs for each particle type
       write(58,*)'  Exclusive neutron only proton only alpha only sdcs'
       write(58,*)
       write(58,*)'  Energy (MeV)  1 Neutron only 1 Proton only 1 Alpha 
     1only'
        write(58,*)
         do ie=1,999
          ee=float(ie)*ed-ed/2.
       sum3=cmen1(1,ie)+cmen1(2,ie)+cmen1(3,ie)
      if(sum3.gt.0.)write(58,100)ee,cmen1(1,ie),cmen1(2,ie),cmen1(3,ie)
         enddo  
c make zero suppress on sdcs output
        write(58,101)sendv1(1),sendv1(2),sendv1(3) 
  100 format(2x,1f6.2,9x,1e9.2,6x,1e9.2,6x,1e9.2)
  101 format('  sum= ',9x,1e9.2,6x,1e9.2,6x,1e9.2)    
c above should print one particle exclusive spectra, energy spectra only
c next write ddcs for single exclusive particle emission
99498  format(' NEUTRON out angular distributions,energy vs. angle deg')
99497 format('PROTON out angular distributions,energy vs. angle deg')
49497 format('ALPHA out angular distributions,energy vs. angle deg')
c begin loop to print 1 particle ddcs-n,p,4He
      if(kinem.gt.0)then
       do 99500 k=1,3
        write(58,*)
        write(58,*)' EXCLUSIVE one particle only DDCS:'
        write(58,*)
        if(k.eq.1)write(58,99498)
        if(k.eq.2)write(58,99497)
        if(k.eq.3)write(58,49497)
      i=k
      jeng=999
      do 99500 index=1,4
      ji=1+(index-1)*9
      jf=ji+8
      do 99501 ian=1,9
      angles(ian)=(float(ian)*5.-2.5)+45.*(float(index)-1.)
99501 continue
      write(58,*)'  '
      write(58,99505)(angles(ian),ian=1,9)
99505 format(' angle: ',9(3x,f7.2))
99503 format(1x,f7.2,1x,9(1x,1pe9.2))
      write(58,*)' energy'
      do 99502 ienn=1,999
      engt=float(ienn)*ed-ed/2.
      st=0.
      do 89504 jn=ji,jf
      st=st+endvcm1(jn,ienn,i)
89504 continue
      if(st.gt.0.)write(58,99503)engt,(endvcm1(jan,ienn,k),jan=ji,jf)
crc
crc
99502 continue
      write(58,69504)(sumthec(ith,k),ith=ji,jf)
69504 format(' integral',9(1x,1pe9.2))
      write(58,*)' the integrals at bottom of columns are the sums of all
     1l cross sections multiplied by '
      write(58,*)' 2pi*sin(theta)*(5./180.)*pi,as bins are 5 degrees wide'
     1e'
99500 continue
      endif
      write(58,*)
      write(58,*)' END ONE PARTICLE OUT EXCLUSIVE SDCS,DDCS; BEGIN 2'
      write(58,*)
      write(58,*)' BEGIN TWO NEUTRON, TWO PROTON, TWO ALPHA SPECTRA' 
      do l=1,6
       do k=1,3
      sendv2(l,k)=0.
           enddo
           enddo
       do ie=1,999
       sumend1(ie)=0.
           enddo
      do l=1,6
       do k=1,3
         do ie=1,999
      sendv2(l,k)=sendv2(l,k)+cmen2(l,k,ie)
          enddo
           enddo
            enddo
c  take sums over energy each particle emission on same coincidence page
       write(58,*)
       write(58,*)'  Energy (MeV)  2 Neutron only 2 Proton only 2 Alpha 
     1only sdcs'
c      write(58,*)'  Energy (MeV)  1 Neutron (and) 1 Proton (only)'  
        write(58,*)
         do ie=1,999
          sumend=cmen2(1,1,ie)+cmen2(3,2,ie)+cmen2(5,3,ie)
          ee=float(ie)*ed-ed/2.
      if(sumend.gt.0.)write(58,100)ee,cmen2(1,1,ie),cmen2(3,2,ie),cmen2(
     15,3,ie)
         enddo  
        write(58,101)sendv2(1,1),sendv2(3,2),sendv2(5,3) 
c  above nn,pp, 4He4He spectra,sdcs. Next do sdcs for np,p4He,n4He sdcs
c after which begin ddcs output.
        
c above should print two particle exclusive spectra, energy spectra only
      do l=2,6,2
       write(58,*)
      if(l.eq.2) write(58,*)'  Energy (MeV)  1 Neutron (and) 1 Proton (o  
     1nly)'
      if(l.eq.4) write(58,*)'  Energy (MeV)  1 Proton (and) 1  Alpha(onl 
     1y)'
      if(l.eq.6) write(58,*)'  Energy (MeV)  1 Neutron (and) 1  Alpha(on 
     1ly)'
         if(l.eq.2)then
          m=1
          n=2
           endif
         if(l.eq.4)then
          m=2
          n=3
           endif
         if(l.eq.6)then
          m=1
          n=3
           endif
       do ie=1,999
          ee=float(ie)*ed-ed/2.
         sumend=cmen2(l,m,ie)+cmen2(l,n,ie)
          if(sumend.gt.0.)write(58,*)ee,cmen2(l,m,ie),cmen2(l,n,ie)
        enddo   
       write(58,*)
        write(58,103)sendv2(l,m),sendv2(l,n) 
       write(58,*)
         enddo
  103 format('  sum= ',9x,1e9.2,6x,1e9.2)    
c next write two particle exclusive ddcs.
c     write(58,*)'  ENDF Exclusive Spectra follow'
c       write(58,101)sendv1(1),sendv1(2),sendv1(3) 
        
c above should print one particle exclusive spectra, energy spectra only
c next write ddcs for single particle emission
98498  format(' NEUTRON out angular distributions,energy vs. angle deg')
98497 format(' PROTON out angular distributions,energy vs. angle deg')
48497 format(' ALPHA out angular distributions,energy vs. angle deg')
       write(58,*)
         if(kinem.gt.0)then
       do 3001 kk=1,3
        do 3001 ithi=1,36
      sumthec(ithi,kk)=0.
 3001 continue
       do l=1,5,2
       if(l.eq.1)write(58,*)'  N-N (two neutron only) DDCS'
       if(l.eq.3)write(58,*)'  P-P (two proton only) DDCS'
       if(l.eq.5)write(58,*)'  4He-4He (two alpha only) DDCS'
       write(58,*)
         if(l.eq.1)j=1
         if(l.eq.3)j=2
         if(l.eq.5)j=3
      do 1003 ith=1,36
      do 1003 ie=1,999
      sumthec(ith,j)=sumthec(ith,j)+endvcm2(l,ith,ie,j)
 1003 continue
      jeng=999
      do 39500 index=1,4
      ji=1+(index-1)*9
      jf=ji+8
      do 39501 ian=1,9
      angles(ian)=(float(ian)*5.-2.5)+45.*(float(index)-1.)
39501 continue
      write(58,*)'  '
      write(58,99505)(angles(ian),ian=1,9)
      write(58,*)' energy'
      do 39502 ienn=1,999
      engt=float(ienn)*ed-ed/2.
      st=0.
      do 29504 jn=ji,jf
      st=st+endvcm2(l,jn,ienn,j)
29504 continue
      if(st.gt.0.)write(58,99503)engt,(endvcm2(l,jan,ienn,j),jan=ji,jf)
crc
crc
39502 continue
       write(58,*)
      write(58,69504)(sumthec(ith,j),ith=ji,jf)
       write(58,*)
39500 continue
       enddo
       do 59503 l=2,6,2
         if(l.eq.2)then
          m=1
          n=2
        write(58,*)' Neutron in coinc with Proton DDCS:'
        write(58,*)
           endif
         if(l.eq.4)then
        write(58,*)' Proton in coinc with Alpha DDCS:'
        write(58,*)
          m=2
          n=3
           endif
         if(l.eq.6)then
          m=1
          n=3
        write(58,*)' Neutron in coinc with Alpha DDCS:'
        write(58,*)
           endif
         do 59499 lsw=1,2
         if(lsw.eq.1)i=m
         if(lsw.eq.2)i=n
          if(lsw.eq.1.and.l.eq.2)write(58,*)' Neutron DDCS (p follows)'
          if(lsw.eq.2.and.l.eq.2)write(58,*)' Proton DDCS(coinc with n)'
          if(lsw.eq.1.and.l.eq.4)write(58,*)' Proton DDCS(alpha follows
     1)' 
          if(lsw.eq.2.and.l.eq.4)write(58,*)' Alpha DDCS(coinc with p)'
          if(lsw.eq.1.and.l.eq.6)write(58,*)' Neutron DDCS(4He follows)'
          if(lsw.eq.2.and.l.eq.6)write(58,*)' Alpha DDCS(coinc with n)'
 
       k=i
       write(58,*)
        if(k.eq.1)write(58,98498)
        if(k.eq.2)write(58,98497)
        if(k.eq.3)write(58,48497)
       write(58,*)
       do 5001 kk=1,3
        do 5001 ithi=1,36
      sumthec(ithi,kk)=0.
 5001 continue
      do 5003 ith=1,36
      do 5003 ie=1,999
      sumthec(ith,i)=sumthec(ith,i)+endvcm2(l,ith,ie,i)
 5003 continue
      jeng=999
      do 59500 index=1,4
      ji=1+(index-1)*9
      jf=ji+8
      do 59501 ian=1,9
      angles(ian)=(float(ian)*5.-2.5)+45.*(float(index)-1.)
59501 continue
      write(58,*)'  '
      write(58,99505)(angles(ian),ian=1,9)
      write(58,*)' energy'
      do 59502 ienn=1,999
      engt=float(ienn)*ed-ed/2.
      st=0.
      do 59504 jn=ji,jf
      st=st+endvcm2(l,jn,ienn,i)
59504 continue
      if(st.gt.0.)write(58,99503)engt,(endvcm2(l,jan,ienn,i),jan=ji,jf)
crc
crc
59502 continue
       write(58,*)
      write(58,69504)(sumthec(ith,i),ith=ji,jf)
       write(58,*)
59500 continue
59499 continue
59503 continue
       write(58,*)
      write(58,*)' the integrals at bottom of columns are the sums of all
     1l cross sections multiplied by '
      write(58,*)' 2pi*sin(theta)*(5./180.)*pi,as bins are 5 degrees wide'
     1e'
        endif
c next write out 3 particle exclusive spectra n,p,4He
c----------+++++++++++++++-----------------3 particle
      write(58,*)
      write(58,*)' BEGIN 3 NEUTRON, 3 PROTON, 3 ALPHA ONLY SPECTRA ' 
      write(58,*)
      do l=1,10
      do k=1,3
      sendv3(l,k)=0.
           enddo
           enddo
       do ie=1,999
       sumend1(ie)=0.
           enddo
      do l=1,10
       do k=1,3
         do ie=1,999
      sendv3(l,k)=sendv3(l,k)+cmen3(l,k,ie)
c     sumend2(l,k,ie)=cmen2(l,k,ie)+sumend2(l,k,ie)
          enddo
           enddo
            enddo
       write(58,*)
       write(58,*)' Exclusive neutron only proton only alpha only sdcs'
       write(58,*)
       write(58,*)'  Energy (MeV)  3 Neutron only 3 Proton only 3 Alpha 
     1only sdcs'
c      write(58,*)'  Energy (MeV)  1 Neutron (and) 1 Proton (only)'  
        write(58,*)
         do ie=1,999
          sumend=cmen3(1,1,ie)+cmen3(2,2,ie)+cmen3(3,3,ie)
          ee=float(ie)*ed-ed/2.
      if(sumend.gt.0.)write(58,100)ee,cmen3(1,1,ie),cmen3(2,2,ie),cmen2(
     13,3,ie)
         enddo  
       write(58,*)
        write(58,101)sendv3(1,1),sendv3(2,2),sendv3(3,3) 
       write(58,*)
c  above nnn,ppp,4He4He4He spectra,sdcs. Next do sdcs for np,p4He,n4He sdcs
c after which begin ddcs output.
        
c above should print three particle exclusive spectra, energy spectra only
      do l=4,9
         if(l.eq.4)then
           m=1
           n=2
            endif
         if(l.eq.5)then
           m=1
           n=3
            endif
         if(l.eq.6)then
           m=1
           n=2
            endif
         if(l.eq.7)then
           m=2
           n=3
            endif
         if(l.eq.8)then
           m=1
           n=3
            endif
         if(l.eq.9)then
           m=2
           n=3
            endif
       
       write(58,*)
      if(l.eq.4) write(58,*)'  Energy (MeV)  2 Neutron (and) 1 Proton (o  
     1nly)'
      if(l.eq.5) write(58,*)'  Energy (MeV)  2 Neutron (and) 1 Alpha (on 
     1ly)'
      if(l.eq.6) write(58,*)'  Energy (MeV)  1 Neutron (and) 2 Proton (o  
     1nly)'
      if(l.eq.7) write(58,*)'  Energy (MeV)  2 Proton (and) 1 Alpha (onl
     1y)'
      if(l.eq.8) write(58,*)'  Energy (MeV)  1 Neutron (and) 2 Alpha (on 
     1ly)'
      if(l.eq.9) write(58,*)'  Energy (MeV)  1 Proton (and)  2 Alpha(onl
     1y' 
       write(58,*)
         do ie=1,999
          sumend=cmen3(l,m,ie)+cmen3(l,n,ie)
          ee=float(ie)*ed-ed/2.
      if(sumend.gt.0.)write(58,100)ee,cmen3(l,m,ie),cmen3(l,n,ie)
         enddo  
       write(58,*)
        write(58,201)sendv3(l,m),sendv3(l,n) 
       write(58,*)
         enddo  
  201 format('  sum= ',9x,1e9.2,6x,1e9.2)    
c have printed three particle out except for n,p,4He sequence
       l=10
       write(58,*)
      if(l.eq.10) write(58,*)'  Energy (MeV) 1 Neutron(and)1 Proton(and
     1) 1 Alpha
     1y' 
       write(58,*)
         do ie=1,999
          sumend=cmen3(10,1,ie)+cmen3(10,2,ie)+cmen3(10,3,ie)
          ee=float(ie)*ed-ed/2.
      if(sumend.gt.0.)write(58,100)ee,cmen3(10,1,ie),cmen3(10,2,ie),
     1cmen3(10,3,ie)
         enddo  
       write(58,*)
        write(58,101)sendv3(10,1),sendv3(10,2),sendv3(10,3)
       write(58,*)
c this should be all sdcs for 3 particles in coincidence; next do ddcs
         if(kinem.gt.0)then
      do l=4,9
        do kk=1,36
         do jj=1,3
           sumthec(kk,jj)=0.
            enddo
             enddo
      do 3003 j=1,3
      do 3003 ith=1,36
      do 3003 ie=1,999
      sumthec(ith,j)=sumthec(ith,j)+endvcm3(l,ith,ie,j)
 3003 continue
         if(l.eq.4)then
           m=1
           n=2
            endif
         if(l.eq.5)then
           m=1
           n=3
            endif
         if(l.eq.6)then
           m=1
           n=2
            endif
         if(l.eq.7)then
           m=2
           n=3
            endif
         if(l.eq.8)then
           m=1
           n=3
            endif
         if(l.eq.9)then
           m=2
           n=3
            endif
         do i=1,2
           if(i.eq.1)k=m
           if(i.eq.2)k=n
       
       write(58,*)
      if(l.eq.4) write(58,*)'  DDCS  2 Neutron (and) 1 Proton (only)'  
      if(l.eq.5) write(58,*)'  DDCS  2 Neutron (and) 1 Alpha (only) '
      if(l.eq.6) write(58,*)'  DDCS  1 Neutron (and) 2 Proton (only)'  
      if(l.eq.7) write(58,*)'  DDCS  2 Proton (and) 1 Alpha (only)'
      if(l.eq.8) write(58,*)'  DDCS  1 Neutron (and) 2 Alpha (only)' 
      if(l.eq.9) write(58,*)'  DDCS  1 Proton (and)  2 Alpha(only)'
       write(58,*)
      if(k.eq.1)write(58,*)' Neutron DDCS'
      if(k.eq.2)write(58,*)' Proton DDCS'
      if(k.eq.3)write(58,*)' Alpha DDCS'
       write(58,*)
      jeng=999
      do 69500 index=1,4
      ji=1+(index-1)*9
      jf=ji+8
      do 69501 ian=1,9
      angles(ian)=(float(ian)*5.-2.5)+45.*(float(index)-1.)
69501 continue
      write(58,*)'  '
      write(58,99505)(angles(ian),ian=1,9)
      write(58,*)' energy'
      do 69502 ienn=1,999
      engt=float(ienn)*ed-ed/2.
      st=0.
      sf=0.
      do 39504 jn=ji,jf
      st=st+endvcm3(l,jn,ienn,k)
39504 continue
      if(st.gt.0.)write(58,99503)engt,(endvcm3(l,jan,ienn,k),jan=ji,jf)
crc
crc
69502 continue
       write(58,*)
      write(58,69504)(sumthec(ith,k),ith=ji,jf)
       write(58,*)
      write(58,*)' the integrals at bottom of columns are the sums of all
     1l cross sections multiplied by '
      write(58,*)' 2pi*sin(theta)*(5./180.)*pi,as bins are 5 degrees wide'
     1e'
69500 continue
       enddo
        enddo
c have completed all 3 particle ddcs except n-p-4He, which we do next, l=10
       write(58,*)
      write(58,*)' 3 particle only emission DDCS-neutron,proton,alpha'
       write(58,*)
        do kk=1,36
         do jj=1,3
           sumthec(kk,jj)=0.
            enddo
             enddo
      do 2003 j=1,3
      do 2003 ith=1,36
      do 2003 ie=1,999
      sumthec(ith,j)=sumthec(ith,j)+endvcm3(10,ith,ie,j)
 2003 continue
        do k=1,3
       write(58,*)
       if(k.eq.1)write(58,*)' Neutron (coinc with p and 4He) DDCS'
       if(k.eq.2)write(58,*)' Proton (coinc with n and 4He) DDCS'
       if(k.eq.3)write(58,*)' Alpha (coinc with n and p) DDCS'
       write(58,*)
      jeng=999
      do 49500 index=1,4
      ji=1+(index-1)*9
      jf=ji+8
      do 49501 ian=1,9
      angles(ian)=(float(ian)*5.-2.5)+45.*(float(index)-1.)
49501 continue
      write(58,*)'  '
      write(58,99505)(angles(ian),ian=1,9)
      write(58,*)' energy'
      do 49502 ienn=1,999
      engt=float(ienn)*ed-ed/2.
      st=0.
      sf=0.
      do 49504 jn=ji,jf
      st=st+endvcm3(l,jn,ienn,k)
49504 continue
      if(st.gt.0.)write(58,99503)engt,(endvcm3(l,jan,ienn,k),jan=ji,jf)
crc
crc
49502 continue
       write(58,*)
      write(58,69504)(sumthec(ith,k),ith=ji,jf)
       write(58,*)
      write(58,*)' the integrals at bottom of columns are the sums of all
     1l cross sections multiplied by '
      write(58,*)' 2pi*sin(theta)*(5./180.)*pi,as bins are 5 degrees wide'
     1e'
49500 continue
       enddo
c4620  format(' sum neutron=',e10.3,'   sum proton =',e10.3,
c    1'  sum alpha=',e10.3,'  sum deut=',e10.3)
c  455 format('  ke  ds/de(mb/mev)=        NEUTRONS    PROTONS     ALPHAS',
c     1      deuterons')
c 460 format(' ',f8.3,'                ',e11.4,'  ',e11.4,'  ',e11.4,' '
c    1,e10.3)
c9999  continue
       endif
c above ends ddcs output conditional on 'kinem'
      return
      end
      subroutine angwrt
        common/mccom/ncount,ipot,rum(15,24),ppmc(15,24,999),cf
      common/memof/sfcrs(999,8),dfsgp(36,999,4),efn(4,999)
      common/fisn/dfslb(36,999,4),cflab(999,4),dfetl(36,999,4)
      common/outlim/limout,lim2
      common/sft5/exc(18,27),xmax
      common/dist/n,ni,ke,irefr,xz(3),pod,pcs,ek,xmix,cfrac
     1,dtau,sign(36),tmix,bji(10),fin,dtup
      common/labcon/prx,pry,prz,pcm,rmas,clab(999,8),dslb(36,999,8)
      common/unit8/angles(9)
      common/mbcfudg/cfg(999)
      common/shft/mc,mp,inver,ike,ipch,kq,qval,ap,at,zp,zt,cld,barfac
      common/pl3/jx,nepr,eb(100),rcsp(100),zee,amass,plex
      common/incr/ed,tened,efinal
      common/send/irfr,iadst,dlt
      common/eder/ed2,v,vhalf,vsq
      common/distb/dsgp(36,999,8)
      common/out1/nout
      common/out4/asm(36),asmm(36)
      common/exclus/ gib(8,999),pr(3),ge(2100),zz(3),
     1delta(3),ts(3),pairx(8),asum(10),rate(8,999)
      common/bill/sumethl(999,8),sumthel(36,8),sumethc(999,8),sumthec(36,8)
     1,8)
      common/bilg/sufethl(999,8),sufthel(36,8),sufethc(999,8),sufthec(36,8)
     1,3)
      common/ist/isot,iikl,knz,idelm,fracts,fract,ccrs(32,13,100)
     1,engbuf(100),exd,een(8,999,100)
      common/tgspin/levno
c
      confac=3.14*5./180.
c     write(17,406)
      if(limout.eq.1)return
      write(11,6562)
      write(58,6562)
      write(42,6562)
6562  format('   ')
      write(11,406)
  406 format('    angular distribution   precompound')
      call headers(7)
c
       do 7000 kk=1,3
        do 7000 i=1,999
      sumethl(i,kk)=0.
      sumethc(i,kk)=0.
      sufethl(i,kk)=0.
      sufethc(i,kk)=0.
 7000 continue
       do 7001 kk=1,3
        do 7001 i=1,36
      sumthec(i,kk)=0.
      sumthel(i,kk)=0.
      sufthec(i,kk)=0.
      sufthel(i,kk)=0.
 7001 continue
      dlt=ed
      nend=xmax/(10.*dlt)
      nen=(xmax+ed)/ed
      do 7003 j=1,3
      do 7003 ith=1,36
      do 7003 ie=1,nen
      sumethc(ie,j)=sumethc(ie,j)+dsgp(ith,ie,j)
      sumethl(ie,j)=sumethl(ie,j)+dslb(ith,ie,j)
      sumthec(ith,j)=sumthec(ith,j)+dsgp(ith,ie,j)
      sumthel(ith,j)=sumthel(ith,j)+dslb(ith,ie,j)
      sufethc(ie,j)=sufethc(ie,j)+dsgp(ith,ie,j)
      sufethl(ie,j)=sufethl(ie,j)+dslb(ith,ie,j)
      sufthec(ith,j)=sufthec(ith,j)+dsgp(ith,ie,j)
      sufthel(ith,j)=sufthel(ith,j)+dslb(ith,ie,j)
 7003 continue
      do 519 i=1,3
      if(i.eq.1)write(17,99885)
      if(i.eq.2)write(17,99886)
      if(i.eq.3)write(17,49886)
crc
      if(i.eq.1)write(11,99885)
      if(i.eq.2)write(11,99886)
      if(i.eq.3)write(11,49886)
      if(i.eq.1)write(58,99885)
      if(i.eq.2)write(58,99886)
      if(i.eq.3)write(58,49886)
      if(i.eq.1)write(18,99885)
      if(i.eq.2)write(18,99886)
      if(i.eq.3)write(18,49886)
      if(i.eq.1)write(40,99885)
      if(i.eq.2)write(40,99886)
      if(i.eq.3)write(40,49886)
      if(i.eq.1)write(42,99885)
      if(i.eq.2)write(42,99886)
      if(i.eq.3)write(42,49886)
       write(18,*)zp,ap,zt,at,isot,levno,' ZP,AP,ZT,AT,ISOT,LEVNO'
       write(11,*)zp,ap,zt,at,isot,levno,' ZP,AP,ZT,AT,ISOT,LEVNO'
      write(11,6562)
       write(58,*)zp,ap,zt,at,isot,levno,' ZP,AP,ZT,AT,ISOT,LEVNO'
      write(58,6562)
       write(42,*)zp,ap,zt,at,isot,levno,' ZP,AP,ZT,AT,ISOT,LEVNO'
       write(40,*)zp,ap,zt,at,isot,levno,' ZP,AP,ZT,AT,ISOT,LEVNO'
      write(42,6562)
      write(40,6562)
      write(18,6562)
      write(17,*)' table entries are in mb/MeV*sr'
      write(11,*)' table entries are in mb/MeV*sr'
      write(58,*)' table entries are in mb/MeV*sr'
      write(40,*)' table entries are in mb/MeV*sr'
      write(42,*)' table entries are in mb/MeV*sr'
      write(18,*)' table entries are in mb/MeV*sr'
      write(11,*)
      write(58,*)
      write(42,*)
      write(18,*)
      write(11,*)' integrals at base of columns are mb/MeV'
      write(58,*)' integrals at base of columns are mb/MeV'
      write(42,*)' integrals at base of columns are mb/MeV'
      write(40,*)' integrals at base of columns are mb/MeV'
      write(18,*)' integrals at base of columns are mb/MeV'
      write(42,*)
      write(18,*)
crc
69503 format(9x,9(1x,1Pe9.2))
99885 format(' monte carlo neutron angular distribution')
99886 format(' monte carlo proton angular distribution')
49886 format(' monte carlo alpha angular distribution')
  401 format('   angle      channel energy (MeV)')
c8401  format('   deg.',1x,f7.2,9(2x,f7.2))
8401  format('   deg.',1x,10(1x,f8.2))
c 401 format(' angle/deg. ke=',f7.2,9(2x,f7.2))
      do 519 j=1,nend
      bji(1)=10.*dlt*float(j-1)+ed/2.
      do 8406 ii=2,10
8406  bji(ii)=bji(ii-1)+dlt
      write(17,401)
      write(18,401)
      write(11,401)
      write(58,401)
      write(42,401)
      write(17,8401)(bji(ii),ii=1,10)
      write(18,8401)(bji(ii),ii=1,10)
      write(58,8401)(bji(ii),ii=1,10)
      write(42,8401)(bji(ii),ii=1,10)
      nii=j*10-9
      nif=nii+9
c
      do 6510 id=1,36
crc
      asm(id)=0.
crc
6510  asmm(id)=0.
c
c
      do 509 ith=1,36
c
      iii=ith
      do 6511 kef=nii,nif
      asmm(iii)=asmm(iii)+dsgp(ith,kef,i)
crc
      asm(iii)=asm(iii)+dslb(ith,kef,i)
      dslb(ith,kef,i)=dslb(ith,kef,i)/(sign(ith)*confac)
      if(kef.le.500)then
      dfslb(ith,kef,i)=dfslb(ith,kef,i)/(sign(ith)*confac)
      dfsgp(ith,kef,i)=dfsgp(ith,kef,i)/(sign(ith)*confac)
      endif
crc
6511  dsgp(ith,kef,i)=dsgp(ith,kef,i)/(sign(ith)*confac)
c check the constants here!
csm
      angle=5.*ith-2.5
      as=asmm(iii)
      if(as.gt.0.)write(17,402)angle,(dsgp(ith,kef,i),kef=nii,nif)
      if(as.gt.0.)write(18,402)angle,(dsgp(ith,kef,i),kef=nii,nif)
crc
       tq=asm(iii)
       if(tq.gt.0.)write(11,402)angle,(dslb(ith,kef,i),kef=nii,nif)
       if(tq.gt.0.)write(58,402)angle,(dslb(ith,kef,i),kef=nii,nif)
       if(tq.gt.0.)write(42,402)angle,(dfslb(ith,kef,i),kef=nii,nif)
crc
  509  continue
c
      write(11,6562)
      write(58,6562)
      write(42,6562)
      write(18,6562)
      write(17,6562)
      write(17,*)
      write(18,*)
      write(11,*)
      write(58,*)
      write(42,*)
      write(17,59503)(sumethc(kef,i),kef=nii,nif)
      write(18,59503)(sumethc(kef,i),kef=nii,nif)
      write(11,59503)(sumethl(kef,i),kef=nii,nif)
      write(58,59503)(sumethl(kef,i),kef=nii,nif)
      write(11,6562)
      write(58,6562)
       write(42,59503)(sufethl(kef,i),kef=nii,nif)
      write(42,6562)
      write(17,6562)
      write(18,6562)
59503 format(' integral ',10(1pe9.2))
crc
crc
  402 format(1x,f7.2,1x,10(1pe9.2))
c6512  format(' integral=',10(1p,e9.2))
c
  519 continue
c add file 8 output for energy distributions
      jeng=(xmax+ed)/ed
      do 99499 i=1,3
      if(i.eq.1)write(18,99498)
      if(i.eq.2)write(18,99497)
      if(i.eq.3)write(18,49497)
crc
      if(i.eq.1)write(12,99498)
      if(i.eq.2)write(12,99497)
      if(i.eq.3)write(12,49497)
       write(12,*)zp,ap,zt,at,isot,levno,' ZP,AP,ZT,AT,ISOT,LEVNO'
       write(8,*)zp,ap,zt,at,isot,levno,' ZP,AP,ZT,AT,ISOT,LEVNO'
      write(12,*)
      write(8,*)
      write(12,*)' table entries are in mb/MeV*sr'
      write(40,*)' table entries are in mb/MeV*sr'
      write(42,*)' table entries are in mb/MeV*sr'
      write(8,*)' table entries are in mb/MeV*sr'
      write(12,*)
      write(8,*)
      write(40,*)' the integrals at bottom of columns are the sums of all
     1cross sections multiplied by '
      write(42,*)' the integrals at bottom of columns are the sums of all
     1cross sections multiplied by '
      write(8,*)' the integrals at bottom of columns are the sums of all
     1cross sections multiplied by '
      write(8,*)' 2pi*sin(theta)*(5./180.)*pi,as bins are 5 degrees wide'
     1'
      write(8,*)
      write(12,*)' the integrals at bottom of columns are the sums of all
     1l cross sections multiplied by '
      write(12,*)' 2pi*sin(theta)*(5./180.)*pi,as bins are 5 degrees wide'
     1e'
      write(8,*)
      write(12,*)
crc
99498  format(' neutron out energy distributions,energy vs. angle deg')
99497 format('proton out energy distributions,energy vs. angle deg')
49497 format('alpha out energy distributions,energy vs. angle deg')
      if(jeng.gt.999)jeng=999
      do 99500 index=1,4
      ji=1+(index-1)*9
      jf=ji+8
      do 99501 ian=1,9
      angles(ian)=(float(ian)*5.-2.5)+45.*(float(index)-1.)
99501 continue
      write(8,99506)
      write(12,99506)
99506 format('   ')
      write(8,99505)(angles(ian),ian=1,9)
      write(42,99505)(angles(ian),ian=1,9)
      write(40,99505)(angles(ian),ian=1,9)
      write(12,99505)(angles(ian),ian=1,9)
99505 format('    angle:',9(1x,f7.2,2x))
99503 format(1x,f7.2,1x,9(1x,1pe9.2))
      write(40,99504)
      write(42,99504)
      write(8,99504)
      write(12,99504)
99504 format('   energy ')
      do 99502 ienn=1,jeng
      engt=(float(ienn))*ed-ed/2.
      st=0.
      do 89504 jn=ji,jf
      st=st+dsgp(jn,ienn,i)
89504 continue
      if(st.gt.0.)write(8,99503)engt,(dsgp(jan,ienn,i),jan=ji,jf)
      if(st.gt.0.)write(40,99503)engt,(dfsgp(jan,ienn,i),jan=ji,jf)
      if(st.gt.0.)write(42,99503)engt,(dfslb(jan,ienn,i),jan=ji,jf)
crc
      st=0.
      do 89505 jn=ji,jf
      st=st+dslb(jn,ienn,i)
89505 continue
      if(st.gt.0.)write(12,99503)engt,(dslb(jan,ienn,i),jan=ji,jf)
crc
99502 continue
99500 continue
      write(8,69503)(sumthec(ith,i),ith=ji,jf)
      write(12,69503)(sumthel(ith,i),ith=ji,jf)
99499 continue
c9999  continue
c4620  format(' sum neutron=',e10.3,'   sum proton =',e10.3,
c     1'  sum alpha=',e10.3,'  sum deut=',e10.3)
c  455 format('  ke  ds/de(mb/mev)=        NEUTRONS    PROTONS     ALPHAS',
c     1      deuterons')
c  460 format(' ',f8.3,'                ',e11.4,'  ',e11.4,'  ',e11.4,' '
c     1,e10.3)
      return
      end
c
C
      subroutine recoil(irecoi)
c  parameter irecoi is added, kkg  04/13/11
c this sr will give lab system energies and angles for emitted particles
      common/debug/alphan
      common/clkin/elab
      common/endvparm/iend,lll(6,6,6),ll(6,6)
      common/history/ide(60),idn(60),idith(60),iddum(60),the(60),adelj
     1(60),spinfh(60),enghis(60)
      common/out1/nout
      common/histindx/ihist
      common/outlim/limout,lim2
 
      common/dec/reccon
      common/recar/thetar,er
      common/pl3/jx,nepr,eb(100),rcsp(100),zee,amass,plex
      common/labcon/prx,pry,prz,pcm,rmas,clab(999,8),dslb(36,999,8)
      common/ilcom/nu,ndum,am(8)
      common/incr/ed,tened,efinal
crc
      common/out2/vrl(200,25),vrlz(200,25),vrlx(200,25),vrly(200,25)
      common/out3/vscale
crc
      common/tem/philabb(36),phirec(36),phb(36),phin(36)
      common/scat/jz,ja, probxx(2,2),q1(8,1100),fcs,eloss,fcspec
cfeb98
      common/star/pcmx,pcmy,pcmz
      common/nuopt/kinem,iem,kb
      common/angprm/aprp(500),eav3(1000),
     1eav2(1000),pot2x,thet1,phi,snt(360),cst(360),snp(720),csp(720)
     2,pot
      common/shft/mc,mp,inver,ike,ipch,kq,qval,ap,at,zp,zt,cld,barfac
c
calculate momentum vectors for emitted particle
c     if(rmas.eq.0.)rmas=amass-1.
c     if(nu.eq.0)nu=1
c the above two statements should be removed
c
c  feb 98 stationary version we leave old version,but return after calc
c  leaving much dead code
c
c this sub starts with z(beam direction) prz momentum vector for recoil= beam momentum
c with pry,prx recoil momenta=0.;these are then changed with each n,p,alpha emission
c and are used to boost each particle emission energy into the lab frame;the recoil
c angles are recomputed at each stage as well,in that their prz,pry,prx components
c are recalculated.nout is the number of particle(i.e. first emission,second,third etc.)
c rmas is the mass number of the residual nucleus left after this particle has been
c emitted.
      if(ap.gt.1.)then
      if(reccon.gt.0.)go to 1000
      endif
c     nout=nout+1
c     if(nout.gt.25)nout=25
c     rmas=rmas+am(nu)
c correct rmas to initial nucleus
       phi=rndm(-1.)*6.28318
      ith=1.+thet1/0.0087267
      iph=1.+phi/0.0087267
c     ip=1.+phi*18./3.1416
c thet,phi are the emission angles wrt z beam axis in recoil frame
c     phin(ip)=phin(ip)+1.
c9999  format(' phi=',f10.3,'iph=',i5,'csp=',f10.3,'snp=',f10.3)
       if(kinem.eq.1)go to 500
       pcmx=pcm*snt(ith)*csp(iph)+pcmx
       pcmy=pcm*snt(ith)*snp(iph)+pcmy
       pcmz=pcm*cst(ith)+pcmz
       return
  500 continue
c add 11/18/2008 def pcm as ejectile c.m. momentum
c     pcm=sqrt(2.*am(nu)*elab)
      pcmx=pcm*snt(ith)*csp(iph)
      pcmy=pcm*snt(ith)*snp(iph)
      pcmz=pcm*cst(ith)
      vcm=pcm/am(nu)
      vcmx=pcmx/am(nu)
      vcmy=pcmy/am(nu)
      vcmz=pcmz/am(nu)
c velocities = p/a(ejectile)
c calculate recoil boosted velocities of ejectiles in lab frame
      vlabx=vcmx+prx/(rmas+am(nu))
      vlaby=vcmy+pry/(rmas+am(nu))
      vlabz=vcmz+prz/(rmas+am(nu))
          plax=vlabx*am(nu)
          play=vlaby*am(nu)
          plaz=vlabz*am(nu)
           prx=prx-plax
           pry=pry-play
           prz=prz-plaz
      vsqlab=vlabx**2+vlaby**2+vlabz**2
      erecoil=(prz**2+pry**2+prx**2)/(2.*rmas)
      elabch=0.5*vsqlab*am(nu)
       elab=elabch
c8/23c     elab=elabch*rmas/(rmas+am(nu))
      argy=sqrt(vlabx**2+vlaby**2)
      thlab=atan2(argy,vlabz)
c  70 thlab=atan(sqrt(vlabx**2+vlaby**2)/vlabz)
c     if(thlab.lt.0.)thlab=thlab+3.14159
      if(vlabx.eq.0.)vlabx=0.00001
c     arg=(vlaby/vlabx)
      philab=atan2(vlaby,vlabx)
      if(philab.lt.0.0) philab = 6.28318531 + philab
c       if(vlabx.lt.0..and.vlaby.lt.0.)philab=3.14159+philab
c       if(vlabx.lt.0..and.vlaby.gt.0.)philab=3.14159+philab
c       if(vlabx.gt.0..and.vlaby.lt.0.)philab=6.28318+philab
c      iphlab=philab*18./3.1416+1.
c      if(iphlab.lt.1) iphlab = 1
c      philabb(iphlab)=philabb(iphlab)+1.
c the above give us theta and phi for light ejectiles in lab frame,phi
c commented out until needed
      iel=(elab+ed)/ed
       if(iel.lt.1)iel=1
c in this routine,rmas should be emittting nucleus mass including emitted particle
      clab(iel,nu)=clab(iel,nu)+fcspec
c store lab single diff cs
      ith=thlab*11.459155+1.
       if(ith.eq.37)ith=36
      dslb(ith,iel,nu)=dslb(ith,iel,nu)+fcspec
c store lab ddcs
       if(iend.eq.2)call endf(thlab,elab,ith,iel,nu )
      return
crc
c above decrement recoil px,y,z for emitted particle momenta
c calculate recoil velocity components,then energy and thetar,phir
1000  continue
      fmas=amass+2.-ja-jz
      if(fmas.lt.4.)fmas=4.
      vx=prx/fmas
      vy=pry/fmas
      vz=prz/fmas
      vrecsq=(vx**2+vy**2+vz**2)
      vrec=sqrt(vrecsq)
      er=0.5*vrecsq*fmas
      ix=vx*vscale+50.
      iy=vy*vscale+50.
      iz=vz*vscale+50.
      it=vrec*vscale+50
      if(iz.gt.200)iz=200
      if(ix.gt.200)ix=200
      if(iy.gt.200)iy=200
      if(it.gt.200)it=200
      if(iz.lt.1)iz=1
      if(ix.lt.1)ix=1
      if(iy.lt.1)iy=1
      if(it.lt.1)it=1
      if(nout.ge.1)then
      vrl(it,nout)=vrl(it,nout)+fcspec
      vrlz(iz,nout)=vrlz(iz,nout)+fcspec
      vrly(iy,nout)=vrly(iy,nout)+fcspec
      vrlx(ix,nout)=vrlx(ix,nout)+fcspec
       endif
      vxy=sqrt(vx**2+vy**2)
      thetar=atan2(vxy,vz)
      phir=atan2(vy,vx)
      if(phir.lt.0.)phir=phir+6.2832
c      if(vx.lt.0..and.vy.lt.0.)phir=phir+3.14159
c      if(vx.gt.0..and.vy.lt.0.)phir=6.28318-phir
c      if(vx.lt.0..and.vy.gt.0.)phir=3.14159-phir
       ipr=1.+phir*18./3.1416
       if(ipr.gt.36)ipr=36
      phirec(ipr)=phirec(ipr)+1.
crc
c scale the xyz,total velocities to go to twice the initial max recoil
c velocity,beginning at index 100 as zero velocity,so that negative
c velocities may be stored.vscale=4.242*sqrt(ep(incident proj energy))/amass
      return
      end
c
      subroutine recoilfm(kst)  !c to do kinematics on fragments from fermi statistics
      common/PRODA/PROD(7,20)
      common/outlim/limout,lim2
      common/fermsum/rfsum,zsum
      common/fermst/aferm
      common/mmtcom/mmt,mzee,ifragmas,ifragzee,light
c     COMMON /PRODA/ PROD(7,20)
      common/recar/thetar,er
      common/pl3/jx,nepr,eb(100),rcsp(100),zee,amass,plex
      common/labcon/prx,pry,prz,pcm,rmas,clab(999,8),dslb(36,999,8)
      common/ilcom/nu,ndum,am(8)
      common/angprm/aprp(500),eav3(1000),
     1eav2(1000),pot2x,thet,phi,snt(360),cst(360),snp(720),csp(720)
     2,pot
      common/incr/ed,tened,efinal
crc
      common/out1/nout
      common/out2/vrl(200,25),vrlz(200,25),vrlx(200,25),vrly(200,25)
      common/out3/vscale
crc
      common/tem/philabb(36),phirec(36),phb(36),phin(36)
      common/scat/jz,ja, probxx(2,2),q1(8,1100),fcs,eloss,fcspec
cfeb98
      common/star/pcmx,pcmy,pcmz
      common/nuopt/kinem,iem,kb
c
calculate momentum vectors for emitted particle
c
c
c this sub starts with z(beam direction) prz momentum vector for recoil= beam momentum
c with pry,prx recoil momenta=0.;these are then changed with each n,p,alpha emission
c and are used to boost each particle emission energy into the lab frame;the recoil
c angles are recomputed at each stage as well,in that their prz,pry,prx components
c are recalculated.nout is the number of particle(i.e. first emission,second,third etc.)
c rmas is the mass number of the residual nucleus left after this particle has been
c emitted.
c prx,pry,prz are recoil momenta of exploding fermi sphere;they need not be changed
c with each emission as per successive cascade emission,as it is a single event
      nout=nout+1
c calc boost for all fragments using velocity of exploding nucleus
      do 80 nst=1,kst
       rf=prod(1,nst)
c     rmas=rmas+am(nu)
c correct rmas to initial nucleus
      ith=1.+thet/0.0087267
      iph=1.+phi/0.0087267
      ip=1.+phi*18./3.1416
c thet,phi are the emission angles wrt z beam axis in recoil frame
c may'02 we will get these angles from px,y,z from prod( array,
c or just use those vectors directly
      phin(ip)=phin(ip)+1.
      pcmx=pcm*snt(ith)*csp(iph)
      pcmy=pcm*snt(ith)*snp(iph)
      pcmz=pcm*cst(ith)
      vcm=pcm/rf
      vcmx=pcmx/rf
      vcmy=pcmy/rf
      vcmz=pcmz/rf
c     vcm=pcm/am(nuf
c     vcmx=pcmx/am(nu)
c     vcmy=pcmy/am(nu)
c     vcmz=pcmz/am(nu)
c velocities = p/a(ejectile)
c calculate recoil boosted velocities of ejectiles in lab frame
      vlabx=vcmx+prx/rfsum
      vlaby=vcmy+pry/rfsum
      vlabz=vcmz+prz/rfsum
      vsqlab=vlabx**2+vlaby**2+vlabz**2
      elab=0.5*vsqlab*rf
c     elab=0.5*m*v**2
      argy=sqrt(vlabx**2+vlaby**2)
      thlab=atan2(argy,vlabz)
c  70 thlab=atan(sqrt(vlabx**2+vlaby**2)/vlabz)
c     if(thlab.lt.0.)thlab=thlab+3.14159
c     if(vlabx.eq.0.)vlabx=0.00001
c     arg=(vlaby/vlabx)
      philab=atan2(vlaby,vlabx)
      if(philab.lt.0.0) philab = 6.28318531 + philab
       iphlab=philab*18./3.1416+1.
       if(iphlab.lt.1) iphlab = 1
       philabb(iphlab)=philabb(iphlab)+1.
c the above give us theta and phi for  ejectiles in lab frame,phi
c commented out until needed
      iel=(elab+ed)/ed
c in this routine,rmas(rfsum) should be emittting nucleus mass including emitted particle
      clab(iel,nu)=clab(iel,nu)+fcspec
c store lab single diff cs
      ith=thlab*11.459155+1.
       if(ith.eq.37)ith=36
      dslb(ith,iel,nu)=dslb(ith,iel,nu)+fcspec
c store lab ddcs
c     prx=prx-pcmx
c     pry=pry-pcmy
c     prz=prz-pcmz
crc
c calculate recoil velocity components,then energy and thetar,phir
c     fmas=rmas-am(nu)
c     vx=prx/fmas
c     vy=pry/fmas
c     vz=prz/fmas
      fmas=rf
       vx=vlabx
       vy=vlaby
       vz=vlabz
      vrec=sqrt(vx**2+vy**2+vz**2)
      er=0.5*vrec**2*fmas
      ix=vx*vscale+50.
      iy=vy*vscale+50.
      iz=vz*vscale+50.
      it=vrec*vscale+50
      if(iz.gt.200)iz=200
      if(ix.gt.200)ix=200
      if(iy.gt.200)iy=200
      if(it.gt.200)it=200
      if(iz.lt.1)iz=1
      if(ix.lt.1)ix=1
      if(iy.lt.1)iy=1
      if(it.lt.1)it=1
      vrl(it,nout)=vrl(it,nout)+fcspec
      vrlz(iz,nout)=vrlz(iz,nout)+fcspec
      vrly(iy,nout)=vrly(iy,nout)+fcspec
      vrlx(ix,nout)=vrlx(ix,nout)+fcspec
      vxy=sqrt(vx**2+vy**2)
      thetar=atan2(vxy,vz)
      phir=atan2(vy,vx)
      if(phir.lt.0.)phir=phir+6.2832
       ipr=1.+phir*18./3.1416
       if(ipr.gt.36)ipr=36
      phirec(ipr)=phirec(ipr)+1.
  80   continue
crc
c scale the xyz,total velocities to go to twice the initial max recoil
c velocity,beginning at index 100 as zero velocity,so that negative
c velocities may be stored.vscale=4.242*sqrt(ep(incident proj energy))/amass
      return
      end
c
      subroutine initial
c     this routine just initializes mc pre-eq related arrays
      common/isomer3/iziso(690),iaiso(690),csgs(690),csis(690)
     1,csgse(690,100),csise(690,100),csis2e(690,100),rm1(690,100),
     2rm2(690,100)
      common/lin/lspin
      common/lscpl/bx(180),sppj
      common/ijcpl/ajtarg,ajin
      common/astp3/niso,no,iadlt,isotag,massno(12),mass1(12),ncsave
      common/outlim/limout,lim2
      common/spins/rote(15,24,181),popej(15,24,180,100),fcx(180),ai,ecmp
      common/sping/pfpej(600,180,100),pfpj(600,180)
      common/gate1/spind(5,180),popj(15,24,180),egatl(5),egatu(5),npart(
     15),igate
      common/ilcom/nu,ndum,am(8)
      common/exclbe/idz(4),ida(4),beff(20,20,2),pairmc(18,27)
      common/sft5/exc(18,27),xmax
      common/lym1/be(18,27,8),pair(18,27),xmas(18,27),delshl(18,27),amas
     1(18,27)
c     common/bym1/symb(18,27),symbp(18,27),symbs(18,27)
      common/incr/ed,tened,efinal
      common/hybsub/iav,idum,rs(8,1100),c,pox,r1
      common/scat/jz,ja, probxx(2,2),q1(8,1100),fcs,eloss,fcspec
      common/exclus/ gib(8,999),pr(3),ge(2100),zz(3),
     1delta(3),ts(3),pairx(8),asum(10),rate(8,999)
      common/out2/vrl(200,25),vrlz(200,25),vrlx(200,25),vrly(200,25)
      common/distb/dsgp(36,999,8)
      common/labcon/prx,pry,prz,pcm,rmas,clab(999,8),dslb(36,999,8)
      common/recoils/rcletza(100,13,24,18),rcleza(100,13,24),
     1rclea(100,36),rcleat(100,36,18),rphir(36)
     2,rclsp(200),rclz(200),rcly(200),rclx(200),escale
     3,rcle(200),rclt(36),rclet(200,36)
c
c
c 3/25/2011 zero pairmc array hmb
      do i=1,18
        do j=1,27
           pairmc(i,j)=0.
             enddo
               enddo
       do 9002 i=1,690
       csgs(i)=0.
       csis(i)=0.
9002   continue
       am(1)=1.
       am(2)=1.
       am(3)=4.
       am(4)=2.
       am(5)=3.
       am(6)=3.
       am(7)=7.
c  calc emission/rescatter branch vs energy
c
       je1=be(1,1+iadlt,1)/ed
       je2=be(1,1+iadlt,2)/ed
       je3=be(1,1+iadlt,3)/ed
       je4=be(1,1+iadlt,4)/ed
       je5=be(1,1+iadlt,5)/ed
       je6=be(1,1+iadlt,6)/ed
       imax=xmax/ed+30.
       if(imax.gt.999)imax=999
c
       do 910 je=1,imax
       l1=je+je1
       l2=je+je2
       l3=je+je3
       l4=je+je4
       l5=je+je5
       l6=je+je6
       if(l1.gt.0)then
       q1(1,je)=rate(1,je)/(rate(1,je)+rs(1,l1))
       else
        q1(1,je)=0.
        endif
       if(l2.gt.0)then
       q1(2,je)=rate(2,je)/(rate(2,je)+rs(2,l2))
        else
         q1(2,je)=0.
          endif
       if(l3.gt.0)then
       q1(3,je)=rate(3,je)/(rate(3,je)+rs(3,l3))
        else
         q1(3,je)=0.
          endif
       if(l4.gt.0)then
       q1(4,je)=rate(4,je)/(rate(4,je)+rs(4,l4))
         else
           q1(4,je)=0.
             endif
        if(l5.gt.0)then
       
       q1(5,je)=rate(5,je)/(rate(5,je)+rs(5,l5))
         else
          q1(5,je)=0.
            endif
        if(l6.gt.0)then
       q1(6,je)=rate(6,je)/(rate(6,je)+rs(6,l6))
          else
            q1(6,je)=0.
             endif
  910  continue
c
c calculate effective energies
c effective means we will not rescatter a nucleon if it is not in
c excess of the binding energy by a constant,presently 5 MeV
       if(no.le.1)then
       do 1800 kz=1,13
       do 1800 ka=1,20
       do 1800 ip=1,2
       beff(kz,ka,ip)=be(kz,ka,ip)
c      beff(kz,ka,ip)=be(kz,ka,ip)+1.
 1800  continue
       do 5 ieng=1,999
       do 5 ipart=1,8
       clab(ieng,ipart)=0.
       do 4 itht=1,36
       dslb(itht,ieng,ipart)=0.
       dsgp(itht,ieng,ipart)=0.
    4  continue
    5  continue
cc add initialization loops on recoil arrays
       do 7771 ie=1,100
       do 7772 ia=1,26
       rclea(ie,ia)=0.
       do 7773 it=1,18
       rcleat(ie,ia,it)=0.
7773   continue
7772   continue
       do 7770 jz=1,13
       do 7769 ja=1,16
       rcleza(ie,jz,ja)=0.
       do 7767 it=1,18
       rcletza(ie,jz,ja,it)=0.
7767   continue
7769   continue
7770   continue
7771   continue
 
crc
       do 993 ie=1,200
       rclsp(ie)=0.
       rclz(ie)=0.
       rclx(ie)=0.
       rcly(ie)=0.
       do 994 jj=1,25
       vrl(ie,jj)=0.
       vrlx(ie,jj)=0.
       vrly(ie,jj)=0.
       vrlz(ie,jj)=0.
 994   continue
 993   continue
              endif
       return
       end
        subroutine precomgam
      common/astp3/niso,no,iadlt,isotag,massno(12),mass1(12),ncsave
      common/nhy/ij,jl,ji,jj,bd,bx1,bx2,cmx,cv,bav,bost,b(3)
      common/sft5/exc(18,27),xmax
      common/incr/ed,tened,efinal
      common/pregam/cff(2,9999),kfn,kfp
      common/hybsub/iav,idum,rs(8,1100),c,pox,r1
      common/exclus/ gib(8,999),pr(3),ge(2100),zz(3),
     1delta(3),ts(3),pairx(8),asum(10),rate(8,999)
      common/shft/mc,mp,inver,ike,ipch,kq,qval,ap,at,zp,zt,cld,barfac
      common/giani/pre, sigt(1100),sigam(999),sigpre(999),gsp(999)
     1,eqgam(999)
      rnex=(zp+zt)/(ap+at)
      rpex=1.-rnex
      rnex=rnex*rnex
      rpex=rpex*rpex
      rnex=rnex*0.315*10.**18
      rpex=rpex*0.315*10.**18
c       write(*,*)'kfn,kfp',kfn,kfp
           if(no.eq.1)then
        do ii=1,999
         sigpre(ii)=0.
          enddo
                      endif
coc convert hybrid photon cascade to MC physics; use cff(e)
c   as cross sect element at excitation 'e' 
      mu=(xmax+ed)/ed
      do 5558 jnu=1,2
        ked=b(jnu)/ed
c loop on excitation energy
      do 5557 jie=1,mu
       le=jie
        ke=jie-ked
        if(ke.lt.1)go to 5557
      if(jnu.eq.1)then 
       rfact=rnex/((rs(jnu,le)+rate(jnu,ke)))*ge(le)
        else
       rfact=rpex/((rs(jnu,le)+rate(jnu,ke)))*ge(le)
        endif
c     do 5556 kee=1,jie
c       ken=mu+1-kee 
        ken=jie
c      multiply gamma constant by 0.25
      sigpre(ken)=cff(jnu,jie)*sigt(ken)*(0.125)*rfact+sigpre(ken)
      continue
 5557 continue
 5558 continue
       return
       end
c---------------------------------------------------------------------- 
      subroutine hybrid
      common/clustr/m3set
      common/lscpl/bx(180),sppj
      common/exci/excdist(700),delexc,counter,ctloop,excdist1(700)
      common/atsv/atsave
      common/astp3/niso,no,iadlt,isotag,massno(12),mass1(12),ncsave
      common/pregam/cff(2,9999),kfn,kfp
      common/tempo/sp3of(999),sp3(999)
      common/ist/isot,iikl,knz,idelm,fracts,fract,ccrs(32,13,100)
     1,engbuf(100),exd,een(8,999,100)
      common/tgspin/levno
       common/fizz/fizbar(106,260),fizcs(100)
      common/fisn/dfslb(36,999,4),cflab(999,4),dfetl(36,999,4)
      common/memof/sfcrs(999,8),dfsgp(36,999,4),efn(4,999)
      common/outlim/limout,lim2
      common/fisn2/tcflab(3),tefn(3),tsfcrs(3)
c
      common/hjk/jang,iq,rcss
      common/distb/dsgp(36,999,8)
      common/rcl/amis
      common/labcon/prx,pry,prz,pcm,rmas,clab(999,8),dslb(36,999,8)
      common/fis/ifis
      common/exclbe/idz(4),ida(4),beff(20,20,2),pairmc(18,27)
      common/scat/jz,ja, probxx(2,2),q1(8,1100),fcs,eloss,fcspec
        common/mccom/ncount,ipot,rum(15,24),ppmc(15,24,999),cf
      common/mul/cor,cres
      common/rc/rzero
      common/giani/pre, sigt(1100),sigam(999),sigpre(999),gsp(999)
     1,eqgam(999)
      common/sms/smscr(999),dep(3)
      common/exclus/ gib(8,999),pr(3),ge(2100),zz(3),
     1delta(3),ts(3),pairx(8),asum(10),rate(8,999)
      common/shft4/k3,jcal
      common/shft5/fiss
      common/shft2/fs(25),dsp(25),brr(25),der(25),er(25)
     1,cx(25),bilsum,xmiss,sumiz
      REAL*8 fiss
      common/shft/mc,mp,inver,ike,ipch,kq,qval,ap,at,zp,zt,cld,barfac
      common/nhy2/bdo,bisp
      common/hybsub/iav,idum,rs(8,1100),c,pox,r1
      common/hyb2/pp(15,24,999)
      common/parfs/moz,k6,delrr(999)
      common/sft9/k5,jmax,pld,d(8),delt,alt
      common/par2/cncs
      common/nhy/ij,jl,ji,jj,bd,bx1,bx2,cmx,cv,bav,bost,b(3)
      common/lab3/sig(8,999),slav(8,999)
      common/pl3/jx,nepr,eb(100),rcsp(100),z,a,plex
      common/par3/ep,sigml(999),acrs(999)
      common/sf6/aoz,zoz,en(8,999)
      common/iso/qpn(3),qpnc
      common/incr/ed,tened,efinal
      common/send/irfr,iadst,dlt
      common/memo/scrs(999,8)
      common/sft5/exc(18,27),xmax
      common/lym1/be(18,27,8),pair(18,27),xmas(18,27),delshl(18,27),amas
     1(18,27)
c     common/bym1/symb(18,27),symbp(18,27),symbs(18,27)
      common/dist/n,ni,ke,irefr,xz(3),pod,pcs,ek,xmix,cfrac
     1,dtau,sign(36),tmix,bji(10),fin,dtup
      common/r34/nr3,nr4,ke5,i3d,irtst,i3t2,ijkl,lq,tem(36)
      common/tst/test
      common/nuopt/kinem,iem,kb
      common/eder/ed2,v,vhalf,vsq
      common/angprm/aprp(500),eav3(1000),
     1eav2(1000),pot2x,thet1,phi1,snx(360),csx(360),sny(720),csy(720)
     2,pot
      ed2=ed/2.
        do i=1,9999
         do j=1,2
          cff(j,i)=0.
           enddo
             enddo
c-----------
       if(no.eq.1)then
      do 7152 i=1,999
      eqgam(i)=0.
      sigam(i)=0.
      sigpre(i)=0.
 7152  continue
         endif
c----------
c  zero ang dist buffer
c----------------
      if(no.eq.1)then
      do 35 jm=1,4
      do 35 ii=1,36
      do 35 kl=1,999
      dfsgp(ii,kl,jm)=0.
   35 dsgp(ii,kl,jm)=0.
                endif
c---------------
c==========
      do j=1,999
      sp3(j)=0.
      sp3of(j)=0.
      enddo
c=========
      bost=0.
      ano=ap+atsave-zp-zt-float(iadlt)
      zee=zp+zt
      z=zee
      a=ap+at-float(iadlt)
      amass=ap+atsave
      if(test.eq.0.)write(17,9898)jcal
      tmix=187./((ap+at)*xmax)
      delta(3)=0.
      rnex=(zp+zt)/(ap+at)
      rpex=1.-rnex
      rnex=rnex*rnex
      rpex=rpex*rpex
      rnex=rnex*0.315*10.**18
      rpex=rpex*0.315*10.**18
c   the above numbers are to modify gamma decay rates for polarization
      pairx(3)=0.
      zz(3)=0.
      ts (3)=0.
      qpn(3)=0.
      pre=1.
      call gama
      pre=0.
      b(3)=0.
      gdo=bdo
      td=bd
      ex1=bx1
      ex2=bx2
      tmx=cmx
      gav=bav
      cost=bost
      cnsig=cncs
      bmx=0.
c     av=cv
c dec98 av=1
      av=1.
      xz(1)=b(1)
      xz(2)=b(2)
c
c    calculate sin table for normalization in 5 deg steps
c
      do 2 i=1,36
      dg=5.*float(i)-2.5
      arg=dg*3.14159/180.
    2 sign(i)=2.*3.14159*sin(arg)
c
      if(jcal.eq.0.and.ifis.ne.5)tmx=1.
      if(jcal.eq.0.and.ifis.eq.5)tmx=0.
c
c   tcg(isoc,isor,i) are clebsch-gordan coefficients for outgoing
c   particles where indices are:isoc =initial state,1=tlower,2=
c   tupper;isor=final state isospin indices,and i=particle type
c   index,1=neutron,2=proton.
c
c   define tcg for n or p or 3he projectiles
c   first neutrons
c
c   zero buffers for four isospin/ke sets
c initialize monte carlo spectra
c---------------------------------
       if(no.le.1)then
       do 66644 ize=1,15
       do 66645 iae=1,24
       rum(ize,iae)=0.
       pairmc(ize,iae)=pair(ize,iae)/10.
       do 66647 ie=1,999
       ppmc(ize,iae,ie)=0.
66647  continue
66645  continue
66644  continue
                 endif
c--------------------------------
c
c      add default exciton parameters
c
      if(ex1.gt.0..or.ex2.gt.0.)go to 1000
      if(ap.eq.1..and.zp.eq.0.)go to 1002
      if(ap.eq.1..and.zp.eq.1.)go to 1003
 1002 exnut=3.*zt+2.*(at-zt)
      ex1=2.*exnut/(exnut+3.*zt)
      ex2=2.-ex1
cccc  temp mod for photonuclear rx
      if(ap.eq.0.) ex1=1.
      if(ap.eq.0.) ex2=1.
      if(ap.eq.0.) td=3.
      if(ap.eq.0.)gdo=0.
      if(ap.eq.0.)bdo=0.
cccc  end temp mod
c
c   define pairing correction per grimes for neutron in/n,p
c
      go to 1004
 1003 exprot=3.*(at-zt)+2.*zt
      ex2=2.*exprot/(exprot+3.*(at-zt))
      ex1=2.-ex2
c
c  define pairing for protons in,n/p out
c
 1004 if(ap.eq.1.)td=3
      cv=1.
      gav=1.00
      if(tmx.eq.0.)cost=1.
      if(tmx.gt.0.)gdo=1.
      if(test.eq.0.)write(17,1001)
      ij=0
      if(ap.eq.1.)go to 1000
      ex1=ap-zp
      ex2=zp
      td=ap+1.
      if(ap.eq.4.)td=ap
1000  continue
      if(tmx.eq.0.)cost=1.
c above dec 98
      cv=1.
      ityp=0
      if(tmx.ne.0..and.gdo.eq.0.)ityp=1
      if(tmx.ne.0..and.gdo.gt.0.)ityp=2
      if(jcal.le.0.and.ifis.ne.0)tcs=sigml(jl)
      if(test.eq.0..and.tmx.eq.0.)write(17,3200)
      if(test.eq.0..and.tmx.ge.1.)write(17,3350)
      if(test.eq.0.)write(17,2700)
      cccc=cost+1.
      if(test.eq.0.)write(17,2750)cccc
      if(test.eq.0..and.cv.eq.0.)write(17,2800)
      if(test.eq.0..and.cv.ge.1.)write(17,2801)
      if(test.eq.0..and.cv.eq.0..and.xmax.gt.55)write(17,2851)
 2851 format(' warning=the imaginary optical potential parameter set may
     1 be invalid above 55 mev channel energy.')
      delcn=pair(1,1)/10.
      pairx(1)=pair(1,2+iadlt)/10.
      pairx(2)=pair(2,1+iadlt)/10.
c
      cor=cncs
      if(ap.eq.0.)tmx=0.
      if (tmx.gt.0..and.jcal.le.0.and.ifis.ne.5) cor=sigml(jl)
      c=cost+1.
c change mfp to saturation value for A>4 projectiles
c     if(ap.gt.4.)c=cost
      ni=td
      do 5 i=1,7
      do 5 ie=1,999
    5 gib(i,ie)=sig(i,ie)/(ed*(ed*float(ie)-ed2))
c
c     calculate efermi averaged over fermi density distribution
c     replace old fermi radius with droplet model value
c
      ecm=ep*at/(at+ap)
      if(ap.gt.0.)delr=4.6/sqrt(ap*ep)
      if(ap.eq.0.)delr=4.6/sqrt(ep)
      cr=1.18*a**.333333
      crr=1.-1./cr**2
      cd=crr*cr+delr
      rmxx=cd-delr
       rrp=2.75+cd
      r2p=1.+exp((rrp-cd)/.55)
      c3=exp((-cd)/.55)+1.
      cave=(1./(1.+(.55/rrp) *(alog(c3)-alog(r2p))))**.666666
      apt=40./cave
      pot=apt
      fin=0.5*29.44/(rmxx*sqrt(ecm))
      fin=fin*180./3.14159
c
c    calculate inverse c.s.*kinetic energy/single particle density
c
      g=a/14.
c
      e=xmax
      cv=1.
c change above to insure nn scattering rates
      av=cv
      iav=av
      m=999
      at=a-ap
      l1=1
      l2=1
c  new relativistic values
      bmas=939.516
      tmas=2.*bmas
      volume=4.*3.14159*at*(rzero**3)/3.
      gfact=8.*(10.**(-6))*3.14*(volume/3.)*(1.602/(6.625*3.0))**3
      do 16 ik=1,2100
      t1=float(ik)*ed-ed
      t2=t1+ed
      ge(ik)=gfact*((t2*(t2+tmas))**1.5-(t1*(t1+tmas))**1.5)
   16   continue
      sm=0.
      do 160 i=1,350
      sm=sm+ge(i)
      kfn=i
      if(sm.gt.ano)go to 161
  160  continue
  161 sm=0.
      do 1700  i=1,350
      sm=sm+ge(i)
      kfp=i
      if(sm.gt.zee)go to 171
 1700  continue
 171  efermn=ed*float(kfn)-ed2
      efermp=ed*float(kfp)-ed2
c     wfact=8.*3.14*(1.602**3)*10.**(-18.)/(3.*(6.6256**3)*9.)
c     wfact=wfact*10.**(34.)
      wfact=1.3*(10.**14)
      do 17 i=1,999
      j=(efermn+b(1))/ed+i
      k=(efermp+b(2))/ed+i
      l=(efermp+be(1,1+iadlt,3))/ed+i
      m=(efermp+be(1,1+iadlt,4))/ed+i
      mm=(efermp+be(1,1+iadlt,5))/ed+i
      mpq=(efermp+be(1,1+iadlt,6))/ed+i
      mpr=(efermp+be(1,1+iadlt,7))/ed+i
      t=float(i)*ed-ed2
      t2=t+ed2
      t1=t2-ed
      beta=sqrt(1.-(939.516/(939.516+t))**2)
      dt=wfact*beta*((t2*(t2+tmas))**1.5-(t1*(t1+tmas))**1.5)
      rate(1,i)=dt*gib(1,i)/(ge(j))
      rate(2,i)=dt*gib(2,i)/(ge(k))
      rate(3,i)=dt*gib(3,i)/(ge(l))
      rate(4,i)=dt*gib(4,i)/(ge(m))
      rate(5,i)=dt*gib(5,i)/(ge(mm))
      rate(6,i)=dt*gib(6,i)/(ge(mpq))
      rate(7,i)=dt*gib(7,i)/(ge(mpr))
c add tentative rates for clusters until more exact are installed
c index3 = alpha
c index4=deuts, change rates according to masses only
c index 5= tritons
c index 6= 3He
c index 7=7Li
   17 continue
c
c   zero depletion  numbers
c
      csn=0.
      csp=0.
      csnn=0.
      cspp=0.
      cspn=0.
c
      nmax=sqrt(1.5*g*e)-4.
      if(nmax.lt.ni)nmax=ni
      smn=0.
      smp=0.
      l2=tmx
      if(l2.le.0)l2=1
      if(jcal.lt.1.and.ifis.ne.5)l1=jl
      if(jcal.lt.1.and.ifis.ne.5)l2=jl
      cres=0.
c
c    begin gdh calculation over partial waves l1 to l2 *****************
c
      l2=l1
      do 300 l=l1,l2
      dep(1)=1.
      dep(2)=1.
      tl=l-1
      if(tmx.le.0.) tl=0.
      cf=sigml(l)
      if(tmx.le.0..and.jcal.ge.0) cf=cncs
      if(cncs.eq.0.)cncs=.000001
      cfrac=cf/cncs
      if(cncs.eq..000001)cfrac=0.
      rmx=rmxx
      continue
c
c     calculate radius of l-th partial wave
c
      ires=0
c
      aq=ap
      if(ap.eq.0.)aq=1.
      r=4.6*tl/sqrt(aq*ep )
      r2=4.6*(tl+1.)/sqrt(aq*ep)
      r1=r
c
9999  continue
c
c    call optical model transition rates
c
      c=cccc
c     if(iav.le.0)call mfp
      xmix=xmax
c
c     if(gav)55,55,50
c
c     calculate local density averages for geometry dependent options
c
      if(r.le.rrp)go to 60
      c1=exp((r-cd)/.55)+1.
      c2=exp((r2-cd)/.55)+1.
      d2=1./c2
      d1=1./c1
      c3=2./(d1+d2)
      dave=1./c3
      c=c3*cccc
      go to 65
   60 c1=exp((r-cd)/.55)+1.
      dave=(1.+(.55/(rrp-r))*(alog(c1)-alog(r2p)))
      c1=1./dave
      c=c1*cccc
   65  pot=40.*dave**.6666666
      if(ires.eq.0)ptemp=pot
      if(ires.eq.0)csav=c
      if(ires.eq.0)rsav=r
      if(pot.le.1.)pot=1.
      pod2=(pot/30.)*(xmax-ecm)+pot
      pod=pot*1.33
c
      r2=r
      r1=r
      if(ires.eq.1)go to 9999
c  estimate refraction when coulomb trajectory bend is present
c
      continue
      n=ni
c
      pot=ptemp
      r=rsav
      c=csav
c
c     call nucleon-nucleon transition rates
c
      if(iav.ge.1)call nucmfp
c
      ert=0.
      srt=0.
      nx=(xmax+ed)/ed
      e1=ex1
      e2=ex2
      tt=td
      pss=0.
      ss=0.
      if(gdo.eq.1..and.tmx.ne.0.)nmax=ni
c
c   add loop on compound nucleus isospin
c
      nisoc=1
      if(ij.eq.1)nisoc=2
c
c     loop on exciton number *********************************************
c
      thetrc=5.72958
      pcon=11.459
      pot=35.
      pol=pot
      pot2x=pot*2.
      v=pot
      vhalf=v/2.
      vsq=v*v
        call trigstuf
      kinem=1
      idz(1)=0
       idz(2)=1
      ida(1)=1
       ida(2)=0
        call initial
      aferm=12.
      fcs=cf/float(ncount)
      fcspec=fcs/ed
c------------------------------------------------------------
      if(ncount.gt.0.and.kinem.eq.0.and.ap.ge.1.)call excl0
c call excl0 for nucleons or HI projectiles, no ddcs
c call excl4 for heavy ion projectiles with ddcs
      if(ncount.gt.0.and.kinem.gt.0.and.ap.gt.1.)call excl4
c call excl1 for nucleon projectiles, ddcs
      if(ncount.gt.0.and.ap.eq.1..and.kinem.eq.1)call excl1
      continue
c call excl3 if projectile is a photon
      if(ncount.gt.0.and.ap.eq.0.)call excl3
c------------------------------------------------------------        
       spsum=0.
       spofsum=0.
c9955 format(1i3,8(1x,e9.3))
  300 continue
c
c    end loop on orbital angular momentum for gdh option
c
      e=xmax
      iemclo=(e+ed)/ed
      if(ncount.gt.0.and.ap.eq.1.)then
      nx=iemclo
      vr=sqrt(2.*ep*ap/((ap+at)**2))
      endif
      if(ncount.gt.0.and.ap.eq.0.)then
      vr=0.03256*ep/amass
      nx=iemclo
      endif
      con=2.*(ap+at-1.)/(ap+at)
      do 305 ke=1,nx
      smn=smn+scrs(ke,1)*ed
  305 smp=smp+scrs(ke,2)*ed
      if(l1.eq.l2.and.tmx.ne.0.) write(17,310) tl
      if(test.ne.0.)go to 335
      write(17,440)smn,smp
c
      if(ike.lt.1)go to 335
c
       if(niso.gt.1.and.no.lt.niso)go to  13350
c///////////////////////
       if(limout.eq.0)then
c-----------------------+++++++++++-------
       if(ncount.gt.0)then
       write(10,*)zp,ap,zt,at,isot,levno,' ZP,AP,ZT,AT,ISOT,LEVNO'
       write(20,*)zp,ap,zt,at,isot,levno,' ZP,AP,ZT,AT,ISOT,LEVNO'
       write(45,*)zp,ap,zt,at,isot,levno,' ZP,AP,ZT,AT,ISOT,LEVNO'
c----------------------------------------------------------
      if(m3set.gt.3)then
       write(55,*)zp,ap,zt,at,isot,levno,' ZP,AP,ZT,AT,ISOT,LEVNO'
       write(54,*)zp,ap,zt,at,isot,levno,' ZP,AP,ZT,AT,ISOT,LEVNO'
       write(56,*)zp,ap,zt,at,isot,levno,' ZP,AP,ZT,AT,ISOT,LEVNO'
       endif
c----------------------------------------------------------
c-------------------------+++++++++++++++
      if(ncount.eq.0)write(17,328)
      if(ncount.gt.0)write(17,329)eb(nepr)
      if(ncount.gt.0)write(10,329)eb(nepr)
      if(ncount.gt.0)write(20,329)eb(nepr)
      if(ncount.gt.0)write(45,329)eb(nepr)
c-------------------------------------------
      if(m3set.gt.3)then
      if(ncount.gt.0)write(55,829)eb(nepr)
      if(ncount.gt.0)write(54,829)eb(nepr)
      if(ncount.gt.0)write(56,829)eb(nepr)
      endif
c------------------------------------------
      endif
c  10/03/2014 add endif here
      endif
c------------+++++++++++++++++---------------
      sum1=0.
      sum2=0.
      sum3=0.
      sum4=0.
      sum5=0.
      sum6=0.
      sum7=0.
      sum8=0.
      sum1l=0.
      sum2l=0.
      sum3l=0.
c aug. 17,06
      sum4l=0.
      sum5l=0.
      sum6l=0.
      sum7l=0.
c end aug 06
      smul1=0.
      smul2=0.
      smul3=0.
      smull1=0.
      smull2=0.
      smull3=0.
c
      do 339 kj=1,999
      sum1=sum1+scrs(kj,1)
      sum2=sum2+scrs(kj,2)
      sum3=sum3+scrs(kj,3)
      sum4=sum4+scrs(kj,4)
      sum5=sum5+scrs(kj,5)
      sum6=sum6+scrs(kj,6)
      sum7=sum7+scrs(kj,7)
      sum1l=sum1l+clab(kj,1)
      sum2l=sum2l+clab(kj,2)
      sum3l=sum3l+clab(kj,3)
      sum4l=sum4l+clab(kj,4)
      sum5l=sum5l+clab(kj,5)
      sum6l=sum6l+clab(kj,6)
      sum7l=sum7l+clab(kj,7)
      wrtest=0.
c-------------------------------------
      do ijk=1,7
      wrtest=wrtest+clab(kj,ijk)
      enddo
c-------------------------------------
      if(wrtest.le.0.)go to 339
      eg=float(kj)*ed-ed2
      ake=eg
      if(ncount.eq.0)go to 337
      sderr1=sqrt(scrs(kj,1)*fcspec)
      slerr1=sqrt(clab(kj,1)*fcspec)
      sderr2=sqrt(scrs(kj,2)*fcspec)
      slerr2=sqrt(clab(kj,2)*fcspec)
      sderr3=sqrt(scrs(kj,3)*fcspec)
      slerr3=sqrt(clab(kj,3)*fcspec)
      write(17,327)eg,scrs(kj,1),sderr1,scrs(kj,2),sderr2
     1,scrs(kj,3),sderr3
      write(45,327)eg,scrs(kj,1),sderr1,scrs(kj,2),sderr2
     1,scrs(kj,3),sderr3
c==========================
      if(m3set.gt.3)then
      test=scrs(kj,4)+scrs(kj,5)+scrs(kj,6)+scrs(kj,7)
c---------------
       if(test.gt.0.)then
      write(55,527)eg,scrs(kj,4),scrs(kj,5),scrs(kj,6),scrs(kj,7)
       endif
c---------------
        endif
c==========================
c temp change above write
      write(10,327)eg,clab(kj,1),slerr1,clab(kj,2),slerr2
     1,clab(kj,3),slerr3
      write(56,527)eg,clab(kj,4),clab(kj,5),clab(kj,6),clab(kj,7)
c     endif
      go to 339
 337  continue
      write(17,325)eg,scrs(kj,1),scrs(kj,2)
 339  continue
      write(44,*)' TOTAL SDCS in LAB  for EMITTED n,p,alpha'
      write(45,*)' TOTAL SDCS in CM  for EMITTED  n,p,4He'
c///////////////////////////////////
       if(m3set.gt.3)then
      write(55,*)' TOTAL SDCS in CM  for EMITTED d,t,3He,7Be'
      write(56,*)' TOTAL SDCS in LAB for EMITTED d,t,3He,7Be'
      endif
      write(41,*)' SDCS in LAB from FISSION for EMITTED n,p,alpha'
      write(39,*)' SDCS in CM  from FISSION for EMITTED n,p,alpha '
      write(39,627)
c     write(41,626)
      write(41,627)
      write(44,627)
      tefn(1)=0.
      tefn(2)=0.
      tefn(3)=0.
      tcflab(1)=0.
      tcflab(2)=0.
      tcflab(3)=0.
      tsfcrs(1)=0.
      tsfcrs(2)=0.
      tsfcrs(3)=0.
c-------------------------------------------------
      do iek=1,990
      sum=cflab(iek,1)+cflab(iek,2)+cflab(iek,3)
      efn(1,iek)=cflab(iek,1)+clab(iek,1)
      efn(2,iek)=cflab(iek,2)+clab(iek,2)
      efn(3,iek)=cflab(iek,3)+clab(iek,3)
      tefn(1)=tefn(1)+efn(1,iek)
      tefn(2)=tefn(2)+efn(2,iek)
      tefn(3)=tefn(3)+efn(3,iek)
      tcflab(1)=tcflab(1)+cflab(iek,1)
      tcflab(2)=tcflab(2)+cflab(iek,2)
      tcflab(3)=tcflab(3)+cflab(iek,3)
      tsfcrs(1)=tsfcrs(1)+sfcrs(iek,1)
      tsfcrs(2)=tsfcrs(2)+sfcrs(iek,2)
      tsfcrs(3)=tsfcrs(3)+sfcrs(iek,3)
c--------------------------------
      if(sum.gt.0.)then
      engy=float(iek)*ed-ed/2.
      write(39,628)engy,sfcrs(iek,1),sfcrs(iek,2),sfcrs(iek,3)
      write(41,628)engy,cflab(iek,1),cflab(iek,2),cflab(iek,3)
      write(44,628)engy,efn(1,iek),efn(2,iek),efn(3,iek)
      endif
c-------------------------------
      enddo
c----------------------------------------------------
      tcflab(1)=tcflab(1)*ed
      tcflab(2)=tcflab(2)*ed
      tcflab(3)=tcflab(3)*ed
      tsfcrs(1)=tsfcrs(1)*ed
      tsfcrs(2)=tsfcrs(2)*ed
      tsfcrs(3)=tsfcrs(3)*ed
          tefn(1)=tefn(1)*ed
          tefn(2)=tefn(2)*ed
          tefn(3)=tefn(3)*ed
       write(39,*)
       write(41,*)
       write(44,*)
       write(39,34556)tsfcrs(1),tsfcrs(2),tsfcrs(3)
       write(41,34556)tcflab(1),tcflab(2),tcflab(3)
       write(44,34556)tefn(1),tefn(2),tefn(3)
       fm1=tsfcrs(1)/fizcs(nepr)
       fm2=tsfcrs(2)/fizcs(nepr)
       fm3=tsfcrs(3)/fizcs(nepr)
       fc1=tcflab(1)/fizcs(nepr)
       fc2=tcflab(2)/fizcs(nepr)
       fc3=tcflab(3)/fizcs(nepr)
       write(39,*)
       write(39,*)'fission neutron multiplicity = ',fm1
       write(39,*)'fission  proton multiplicity = ',fm2
       write(39,*)'fission  alpha  multiplicity = ',fm3
       write(41,*)
       write(41,*)'fission neutron multiplicity = ',fm1
       write(41,*)'fission  proton multiplicity = ',fm2
       write(41,*)'fission  alpha  multiplicity = ',fm3
       write(39,*)
       write(41,*)
      sum1=sum1*ed
      sum2=sum2*ed
      sum3=sum3*ed
      sum4=sum4*ed
      sum5=sum5*ed
      sum6=sum6*ed
      sum7=sum7*ed
      sum1l=sum1l*ed
      sum2l=sum2l*ed
      sum3l=sum3l*ed
      smul1=sum1/cncs
      smul2=sum2/cncs
      smul3=sum3/cncs
      smull1=sum1l/cncs
      smull2=sum2l/cncs
      smull3=sum3l/cncs
c------------------------------------
       if(limout.eq.0)then
      write(17,44556)sum1,sum2,sum3
      write(45,44556)sum1,sum2,sum3
c000000000000000000
      if(m3set.gt.3)then
      write(55,44559)sum4,sum5,sum6,sum7
      endif
c000000000000000000
      write(17,44557)smul1,smul2,smul3
      write(10,44557)sum1l,sum2l,sum3l
      write(56,44559)sum4l,sum5l,sum6l,sum7l
      write(10,44557)smull1,smull2,smull3
      endif
c----------------------------------
c      endif
c
      continue
c
  335 if(iadst.eq.0.and.ncount.eq.0)go to 5200
      confac=3.14*5./180.
      nend=xmax/(10.*dlt)+0.99
c
      do 522 i=1,2
      do 520 j=1,nend
      bji(1)=dlt*float(10*j-9)-ed/2.
      do 8406 ii=2,10
8406  bji(ii)=bji(ii-1)+dlt
c
      continue
c
c------------------------------------
       if(limout.eq.0)then
      write(17,401)(bji(ii),ii=1,10)
      endif
c-----------------------------------
      nii=j*10-9
      nif=nii+9
c
      do 6510 id=1,10
6510  asum(id)=0.
c
c
      do 510 ith=1,36
c
      iii=0
 
      do 6511 kef=nii,nif
      iii=iii+1
6511  asum(iii)=asum(iii)+sign(ith)*dsgp(ith,kef,i)*confac
c
      as=asum(iii)
      angle=5.*ith-2.5
c----------------------------
       if(limout.eq.0)then
      if(as.gt.0.)write(17,402)angle,(dsgp(ith,kef,i),kef=nii,nif)
      endif
c---------------------------
  510 continue
c
c
  520 continue
  522 continue
5200  continue
c
      if(ike.ne.4)go to 3511
      do 351 i=1,999
      en(1,i)=en(1,i)+scrs(i,1)
c     en(4,i)=en(4,i)+scrs(i,4)
 351  en(2,i)=en(2,i)+scrs(i,2)
c
 3511 continue
c
c   set applicable reaction cs for type of calculation
c
13350 continue
      if(jcal.lt.0.)rcs=tcs
      if(jcal.ge.0)rcs=cnsig
      if(jcal.gt.0.and.ncount.gt.0)rcs=rum(1,1)
      if(jcal.gt.0.and.ncount.gt.0)cncs=rum(1,1)
      if(jcal.eq.0.and.ifis.eq.5.and.ncount.gt.0)rcs=rum(1,1)
      if(jcal.eq.0.and.ifis.eq.5.and.ncount.gt.0)cncs=rum(1,1)
      write(17,9899)cncs
      write(17,9898)jcal
c
c6666  format(5f5.2)
cc6667  format('  rzero value =',f5.2)
 328  format(' k.e. (MeV) N-cs mb/MeV P-cs mb/MeV ')
  329  format(' @',f8.2,' MeV:',/,                                      
     1       ' k.e. (MeV) N-cs mb/MeV one sigma   P-cs mb/MeV one sigma 
     2  Alpha cs mb/MeV one sigma')                                     
  829  format(' @',f8.2,' MeV:',/,                                     
     1       ' k.e. (MeV) D-cs mb/MeV T-cs mb/MeV 3He cs mb/MeV 7Be cs 
     2mb/MeV')
 627  format(' energy  neuts(mb/MeV) prots(mb/MeV) alphas(mb/MeV)')
c927  format(' energy  deuts(mb/MeV) triton(mb/MeV)    3He(mb/MeV),
c    1  7Be (mb/MeV')
 527  format(1x,f6.2,4(3x,1pe10.1))
 628  format(3x,1f7.2,3(1x,1pe10.2))
44556 format('  neut sum= ',1pe9.2,5x,'prot sum= ',1pe9.2,3x,
     1' alpha sum= ',1pe9.2)
44559 format('2H sum=',1pe9.2,1x,'3H sum=',1pe9.2,1x,
     1'3He sum=',1pe9.2,1x,'7Be sum=',1pe9.2)
34556 format('  SUM(mb) ',3(2x,1pe9.2))

44557 format(' neut mlt=',1pe9.2,' prot mlt=',1pe9.2,
     1' alpha mlt=',1pe9.2)
c9340  format('  kj =',i5,'recoil energy(MeV)=',1pe9.3,
c     1' recoil cs =',1pe9.3)
  401 format(' angle/deg. ke=',f6.2,9(4x,f6.2))
  402 format(1x,f5.1,5x,10(e10.3))
9899  format(' cncs fron hybrid mc =',e10.2)
9898  format(50x,' jcal=',i5)
c13101 format(//' csn=',e10.3,'csp=',e10.3,'csnn=',e10.3,'cspp=',e10.3,
c     1'cspn=',e10.3)
 2801 format(/41x,'transition rates are from nucleon-nucleon scattering')
     1)
  440 format(/' total precompound neutron cross section = ',e10.3,
     1 ', total precompound proton cross section = ',e10.3)
 1001 format(//30x,' precompound parameters have been selected internally'/)
     1y'/)
 2700 format(/20x,' the following parameter values were selected for the prec
     1 precompound calculation')
 2750 format(/40x,' mean free path multiplier is = ',f10.3,/)
 2800 format(' intranuclear transition rates were calculated using the i
     1maginary optical potential')
c2900 format(' td= ',f5.2,' ex1= ',f5.2,' ex2= ',f5.2,' tmx= ',f5.2,' av
c    1= ',f5.2,' gav= ',f5.2,' cost= ',f5.2,' gdo= ',f5.2,' ij= ',i5)
 3200 format(/41x,' hybrid calculation has been selected ')
 3350 format(/41x,' geometry dependent hybrid model selected')
  310 format (/' precompound cross sections for partial wave l=',
     1 f3.0,' only')
cneadb  320 format(x,f7.2,x,'ds/de neutrons',5(1x,e9.3),2x,5(1x,e9.3))
  325 format(1x,f9.3,2x,1e10.3,2x,1e10.3)
  327 format(1x,f9.3,6(2x,1pe10.3))
c---------------temp endif-----------------
c       endif
       return
      end
 
      subroutine partid(ek,ipart,n1,n11,nhol)
      common/outlim/limout,lim2
c  subroutine to give nucleon/hole identities (n,n,n-1..)
c  when one nucleon scatters in fermi sea
c  returns 'n1' as first exciton of 2p1h,'nn1' as second,and nhol as hole
      common/edep/pxx(2,1100)
      common/incr/ed,tened,efinal
      common/excid/npartt(2)
      if(ipart.lt.1)ipart=1
      n1=ipart
      npart=npartt(ipart)
      ixmx=(ek+ed)/ed
      if(ixmx.lt.1)ixmx=1
      x=rndm(-1.)
      if(x.gt.pxx(ipart,ixmx))go to 2
c set case of like-like scatter
      n11=ipart
      nhol=ipart
      return
    2 n11=npartt(ipart)
      nhol=n11
c now for this unlike scatter case,randomize which nuc we consider first
      x=rndm(-1.)
      if(x.gt.0.5)return
      n1=n11
      n11=ipart
      return
      end
      subroutine clustin1(tph,pph)
c this is routine to provide input to excl for projectiles
c with mass greater than 1
c it will receive ap,zp,the projectile mass and charge
      common/tr/xtran
      common/eferm/efermt
      common/sq/sqx,efds(40,11),pfermx(40,11)
      common/arrays/eout(300),efd(40,11),pferm(40,11)
     1,tsetlb(181,300),edist(300),eric(300),ie(40),itl(40),
     2pfpe(40,11),efinc(40,11),eferm2(40,11),phirad(100),peincsq,peinc
      common/tang3/efermp,eincc,nexc
      common/incr/ed,tened,efinal
      common/exci2/excdis2t(700)
      common/exci/excdist(700),delexc,counter,ctloop,excdist1(700)
      common/clust1/denomc,connum,ptproj
      common/labcon/prx,pry,prz,pcm,rmas,clab(999,8),dslb(36,999,8)
      common/himbc/th2,ph2,przz,px,py,pz
      common/himbc2/clthet(100),clphi(100)
      common/outlim/limout,lim2
      common/eder/ed2,v,vhalf,vsq
      common/nhy/ij,jl,ji,jq,bd,bx1,bx2,cmx,cv,bav,bost,b(3)
      common/bound/neut2,iprot,idnuc(100),idm(100),einc(100),
     1thett(100),phit(100),ppp
      common/shft/mc,mp,inver,ike,ipch,kq,qval,ap,at,zp,zt,cld,barfac
      common/sft5/exc(18,27),xmxx
      common/scat/jz,ja, probxx(2,2),q1(8,1100),fcs,eloss,fcspec
      common/exclbe/idz(4),ida(4),beff(20,20,2),pairmc(18,27)
      common/memo/scrs(999,8)
c oct 05-should add selection of angles for the HI excitons
c based on mbc model.
c     do 2 i=1,40
c     clthet(i)=0.
c     clphi(i)=0.
c     px=0.
c     py=0.
c     pz=0.
c here start px,y,z=0, then sum for nucleons out, put in initial
c values to get last nucleon parms
c set initial momentum to get final exciton by conservation
c     e=xmax
      e=xtran
      app=ap
      zpp=zp
      i=0
      ia=1
      iz=1
c case of some hole excitons as well as particle excitons
      n0=bd
      chol=bd-app
c xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      do 700 n=n0,3,-1
c xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      fhol=chol/bd
      icycle=0
      i=i+1
c branch if hole is chosen
      root=1./(float(n)-1.)
      if(rndm(-1.).gt.fhol)go to 650
      n1=3
      chol=chol-1.
 199     eps=e*(1.-rndm(-1.)**root)
      if(eps.gt.v)go to 199
      go to 201
 650     n1=2
      fp=zpp/app
      if(rndm(-1.).gt.fp)n1=1
 299  continue
      mz=1+idz(n1)
      ma=1+ida(n1)
      icycle=icycle+1
      eps=e-(rndm(-1.)**root)*(e)
      if(eps.gt.e.and.icycle.lt.2)go to 299
      if(eps.gt.e)eps=e
      if(eps.lt.0.)eps=0.
c might instead recycle on above condition
      if(n1.eq.2)zpp=zpp-1.
 201     e=e-eps
      einc(i)=eps
      idnuc(i)=n1
      app=app-1.
c     call recoilhi
c in recoilhi keep track of px,py,pz of emitted nucleons(whether or not
c actually emitted, to get andle of last exciton by momentum balance
cc error- the i and nc indices must be identical
c     clthet(i)=th2
c     clphi(i)=ph2
c the above stores are done in rohi
c the theta,phi of this emitted nucleon are returned, but not rotated
c here 'e' is the remaining excitation, and app the cluster size remaining
 700  continue
      i=i+1
      n1=1
      n11=1
      if(zpp.eq.2.)then
      n1=2
      n11=2
      endif
      if(zpp.eq.1.)then
      n1=1
      n11=2
      endif
      if(app.eq.1)then
c programmed for at most one hole remaining
      if(zpp.eq.1)then
      n1=2
      n11=3
      endif
      if(zpp.eq.0.)then
      n1=1
      n11=3
      endif
      endif
      if(app.eq.0.)then
      n1=3
      n11=3
      endif
      eps=rndm(-1.)*e
      einc(i)=eps
      e=e-eps
      idnuc(i)=n1
      app=app-1.
      ep2=e
      i=i+1
      einc(i)=ep2
      idnuc(i)=n11
      e=e-ep2
      l=bd
      do 9000 j=1,l
      jj=(rndm(-1.)*bd+0.9999999)
      jjj=jj
      k=(rndm(-1.)*10.+0.999999)
      if(jj.lt.1)jj=1
      if(k.lt.1)k=1
 
          arg=(einc(j)-eincc-efd(jj,k))/(pfermx(jj,k))
       if(abs(arg).gt.1.)then
cxxxxxxxxxxxxx

        k=10
        jj=jjj+1+(l-jjj-1)*rndm(-1.)
          arg=(einc(j)-eincc-efd(jj,k))/(pfermx(jj,k))
       if(abs(arg).gt.1.)then
        jj=l
        k=10
          arg=(einc(j)-eincc-efd(jj,k))/(pfermx(jj,k))
       endif
       endif
        if(abs(arg).gt.1.)arg=arg/abs(arg)
       arg2=(sqx+pferm(jj,k)*arg)/sqrt(2.*(einc(j)+efermt))
       if(abs(arg2).gt.1.)then
       jj=l
       k=10
          arg2=(sqx+pferm(jj,k)*arg)/sqrt(2.*(einc(j)+efermt))
       endif
       if(abs(arg2).gt.1.)arg2=arg2/abs(arg2)
      thett(j)=acos(arg2)
c     costhet=(einc(j)-connum)/denomc
c     pfinal=sqrt(2.*(einc(j)+efermt))
c     thett(j)=acos(sqx+ptproj*costhet/pfinal)
      phit(j)= rndm(-1.)*6.28318
 9000     continue
        iap=ap
       do 9065 ll=1,iap
       ix=(einc(ll)+ed)/ed
       if(ix.lt.1)ix=1
       excdist1(ix)=excdist1(ix)+delexc
9065   continue
c      e=0.
c seems that e should be zero here- might check
      
      return
      end
 
      subroutine clusterin(tph,pph)
c this is routine to provide input to excl for projectiles
c with mass greater than 1
c it will receive ap,zp,the projectile mass and charge
      common/exci/excdist(700),delexc,counter,ctloop,excdist1(700)
      common/nuopt/kinem,iem,kb
      common/tr/xtran
      common/eferm/efermt
      common/incr/ed,tened,efinal
      common/labcon/prx,pry,prz,pcm,rmas,clab(999,8),dslb(36,999,8)
      common/himbc/th2,ph2,przz,px,py,pz
      common/himbc2/clthet(100),clphi(100)
      common/outlim/limout,lim2
      common/eder/ed2,v,vhalf,vsq
      common/nhy/ij,jl,ji,jq,bd,bx1,bx2,cmx,cv,bav,bost,b(3)
      common/bound/neut2,iprot,idnuc(100),idm(100),einc(100),
     1thett(100),phit(100),ppp
      common/shft/mc,mp,inver,ike,ipch,kq,qval,ap,at,zp,zt,cld,barfac
      common/sft5/exc(18,27),xmxx
      common/scat/jz,ja, probxx(2,2),q1(8,1100),fcs,eloss,fcspec
      common/exclbe/idz(4),ida(4),beff(20,20,2),pairmc(18,27)
      common/memo/scrs(999,8)
c oct 05-should add selection of angles for the HI excitons
c based on mbc model.
         tph=0.
         pph=0.
         th2=0.
         ph2=0.
      do 2 i=1,40
      clthet(i)=0.
      clphi(i)=0.
 2    continue
      px=0.
      py=0.
      pz=0.
c here start px,y,z=0, then sum for nucleons out, put in initial
c values to get last nucleon parms
c set initial momentum to get final exciton by conservation
c     e=xmax
       e=xtran
      app=ap
      zpp=zp
      i=0
      ia=1
      iz=1
c case of some hole excitons as well as particle excitons
      if(ap.gt.1.)bd=ap
      n0=bd
      chol=bd-app
c xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      do 700 n=n0,3,-1
c xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      fhol=chol/bd
      icycle=0
      i=i+1
c branch if hole is chosen
      root=1./(float(n)-1.)
      if(rndm(-1.).gt.fhol)go to 650
      n1=3
      chol=chol-1.
c 199     x=rndm(-1.)
c     r=log(x)/(n-1)
c     eps=e-exp(r)*e
 199     eps=e*(1.-rndm(-1.)**root)
      if(eps.gt.v)go to 199
      go to 201
 650     n1=2
      fp=zpp/app
      if(rndm(-1.).gt.fp)n1=1
 299  continue
c299     x=rndm(-1.)
c     r=log(x)/(n-1)
      mz=1+idz(n1)
      ma=1+ida(n1)
      icycle=icycle+1
      eps=e-(rndm(-1.)**root)*(e)
      if(eps.gt.e.and.icycle.lt.2)go to 299
      if(eps.gt.e)eps=e
      if(eps.lt.0.)eps=0.
c might instead recycle on above condition
      if(n1.eq.2)zpp=zpp-1.
c     app=app-1.
 201     e=e-eps
      einc(i)=eps
      idnuc(i)=n1
c oct 05 add call to calc angles via rohingl
      pcm=sqrt(2.*eps)
      if(kinem.eq.1)call rohi(pdum,e,eps,tph,pph,elab,app)
      app=app-1.
c     call recoilhi
c in recoilhi keep track of px,py,pz of emitted nucleons(whether or not
c actually emitted, to get andle of last exciton by momentum balance
cc error- the i and nc indices must be identical
c     clthet(i)=th2
c     clphi(i)=ph2
c the above stores are done in rohi
c the theta,phi of this emitted nucleon are returned, but not rotated
c here 'e' is the remaining excitation, and app the cluster size remaining
 700  continue
      i=i+1
      n1=1
      n11=1
      if(zpp.eq.2.)then
      n1=2
      n11=2
      endif
      if(zpp.eq.1.)then
      n1=1
      n11=2
      endif
      if(app.eq.1)then
c programmed for at most one hole remaining
 
      if(zpp.eq.1)then
      n1=2
      n11=3
      endif
      if(zpp.eq.0.)then
      n1=1
      n11=3
      endif
      endif
      if(app.eq.0.)then
      n1=3
      n11=3
      endif
      eps=rndm(-1.)*e
      einc(i)=eps
      e=e-eps
      if(e.lt.0..or.eps.lt.0.)then
      endif
      idnuc(i)=n1
c     app=app-1.
      pcm=sqrt(2.*eps)
      if(kinem.eq.1)call rohi(pdum,e,eps,tph,pph,elab,app)
      app=app-1.
c     call recoilhi
c     clthet(i)=th2
c     clphi(i)=ph2
      ep2=e
      i=i+1
      einc(i)=ep2
      idnuc(i)=n11
c here for last nucleon probably must get last angle by p cons.
      pz=przz-pz
      py=-py
      px=-px
      argy=sqrt(px**2+py**2)
      th2=atan2(argy,pz)
      ph2=atan2(py,px)
c for last nucleon, the angles should be in recoil-Zbeam frame
      if(ph2.lt.0.0) ph2 = 6.28318531 + ph2
      clthet(i)=th2
      clphi(i)=ph2
      e=e-ep2
       e=0.
c seems that e should be zero here- might check
      return
      end
       subroutine hiin
c start dec 2000 to sub this for clusterin
c dec 2000 this is edited version of a17.f to give only energy spectra
c sept 2000 modify to give choice of equidistant
c projectile exciton spacing,or Fermi gas as input choice
c sept 12,2000 modify to cycle n-th exciton to be chosen
c by difference
c use equidistant ef for projectile,select energy
c in center of interval,or center of fermi gas interval
c this version goes to next event when total energy
c exceeds available energy
c this version will select n excitons from unrestricted
c energy coupling,where energy is determined by choosing
c a random angle for exciton having momentum calc'd from
c Eincident per nucleon plus Efermi-target,coupling with
C Fermi momentum of projectile
c version to insist on combos where only last exciton is
c selected to conserve energy
c 5/22/00 modify to include target well depth adding to einc
c per nucleon as per af.f version
      common/incr/ed,tened,efinal
      common/tr/xtran
      common/outlim/limout,lim2
      common/arrays/eout(300),efd(40,11),pferm(40,11)
     1,tsetlb(181,300),edist(300),eric(300),ie(40),itl(40),
     2pfpe(40,11),efinc(40,11),eferm2(40,11),phirad(100),peincsq,peinc
       common/tang2/anexc(40),anexci(40),eferm(40,11),psq(40,11),
     1emit(300),angl(181),thetlb(37,300),pdist(300),pd(300)
      common/shft/mc,mp,inver,ike,ipch,kq,qval,ap,at,zp,zt,cld,barfac
      common/tang3/efermp,eincc,nexc
      common/sft5/exc(18,27),emax
      common/angprm/aprp(500),eav3(1000),
     1eav2(1000),pot2x,thet1,phi1,snt(360),cst(360),snp(720),csp(720)
     2,pot
      common/eferm/efermt
      common/bound/neut2,iprot,idnuc(100),idm(100),eink(100),
     1thett(100),phit(100),ppp
      common/stort/storth(51),stor(51)
      common/sq/sqx,efds(40,11),pfermx(40,11)
      common/exci/excdist(700),delexc,counter,ctloop,excdist1(700)
      xmax=emax
   8  continue
        tout=0.
       lim=ap
      do 191 ike=1,lim
      emit(ike)=0.
      angl(ike)=0.
      phirad(ike)=0.
c
 191  continue
      do 60 n=1,nexc-1
      if(tout.gt.xmax)go to 8
      ji=n
      k=1.+9.*rndm(-1.)
c      jx=n
c     ji=n+jd
c ji is selected to id the exciton cycling whichis the first efermi
c     if(ji.gt.nexc)ji=ji-nexc
c   eincc is energy per nucleon of incident ion,actually xmax/ap
c   thet is angle of fermi energy wrt incident ion,in radians
c   costhet is cos( thet)
  12   continue
       costhet=(1.-2.*rndm(-1.))
c substitute pfermx= pferm*sqrt(2.*(efermt+eincc)), and efds() for efd+sqx
       out=efds(ji,k)+costhet*pfermx(ji,k)
       if(out.gt.emax.or.out.lt.0.)go to 12
c out is exciton energy in lab,i.e above efermi where binding
c is still being assumed to be zero
c if the emission energy exceeds the total initial energy,
c select another random angle theta
      tout=out+tout
      if(tout.gt.xmax)go to 8
c if selected energy exceeds emax,start new event
c     thlab=acos((peinc+pe*costhet)/sqrt(poutsq))
c     costhl=(peinc+pferm(ji,k)*costhet)/out
      etot=out+efermt
      ptot=sqrt(2.*etot)
      costhl=(costhet*pferm(ji,k)+sqx)/ptot
       if(abs(costhl).gt.1.)costhl=costhl/(abs(costhl))
      thlab=acos(costhl)
c now put these two results into arrays to store for output
c save lab angle and energy indices of n-1 excitons
cc      ithl=const*thlab+1.
c we will store angles on cosines to save time,then when all good
c events are collected,we will do arc cos functions on the final
c costhet arrays-this should speed the routine
      angl(ji)=costhl
      eink(ji)=out
  60  continue
      if(tout.gt.xmax)go to 8
c now do last exciton by energy conservation
      jlast=ji+1
c     if(jlast.gt.nexc)jlast=jlast-nexc
      aout=xmax-tout
      k=1.+ 9.*rndm(-1.)
c     costh =(aout-eincc-efd(jlast,k))/(sqx*pferm(jlast,k))
      costh =(aout-efds(jlast,k))/pfermx(jlast,k)
      if(abs(costh).gt.1.)go to 8
      costhl=(sqx+pferm(jlast,k)*costh)/sqrt(2.*(aout+efermt))
      angl(jlast)=costhl
      eink(jlast)=aout
c now add nucleon id loops,i.e. neuts or prots
      app=ap
      zpp=zp
       do 65 inex=1,nexc
       thett(inex)=acos(angl(inex))
       phit(inex)=rndm(-1.)*6.2832
       ith=thett(inex)*11.459155+1.
       if(ith.eq.37)ith=36
       storth(ith)=storth(ith)+1.
       fp=zpp/app
       if(rndm(-1.).gt.fp)go to 7
       idnuc(inex)=2
       zpp=zpp-1.
       app=app-1.
       go to 18
   7   idnuc(inex)=1
       app=app-1.
  18   continue
       phirad(inex)=rndm(-1.)*6.2831853
  65   continue
        iap=ap
       counter=counter+1.
       do 9065 ll=1,iap
       ix=(eink(ll))/ed+1.
       excdist(ix)=excdist(ix)+delexc
9065   continue
       return
       end
       subroutine ijlcpl
      common/spins2/spinu,spinf,prz1
      common/tagspin/trgspin
      common/lscpl/bx(180),sppj
      common/ijcpl/ajtarg,ajin
       common/ijcouple/p(90),b(90)
      common/kkg2/fractot(99),cs90
      common/spinc/psp(90,90,181)
c---------------------------------------------------
c        do i=1,90
c         p(i)=0.
c         b(i)=0.
c            enddo
c----------------------------------------------------
c      aji=abs(trgspin-spinu)
c      ajf=trgspin+spinu
c       ji=2.*aji+1.
c       jf=2.*ajf+1.
c=====================================================
c------------------------------------
c                    sum=(float(ji+jf)/2.)*(float(jf-ji)/2.+1.)
c         p(ji)=float(ji)/sum
c-----------------------------------
c        do j=ji+2,jf,2
c         p(j)=float(j)/sum+p(j-2)
c           enddo

c       do j=ji,jf,2
c        spinu=float(j-1)/2.
c         if(x.le.p(j))go to 1010
c          enddo
c1010  continue
       aji=abs(trgspin-spinu)
       ajf=trgspin+spinu
        ji=2.*aji+1.
        jf=2.*ajf+1.
c 11 june 2017 fix for jf.gt.90?
        if(jf.gt.90)jf=90
        if(ji.gt.90)ji=90
cd 10/25/2014          x=ranf(-1.)
          x=rndm(-1.)
        do j=ji,jf,2
         spinu=float(j-1)/2.
          if(x.le.psp(ji,jf,j))go to 1020
            enddo
 1020 continue
c-----------------------------------
c-----------------------------------
c set sumij array to be running sum of 2j+1 probs for coupling I-j to
c I+j where I is nuclear spin and j is ls coupled incoming channel spin
c           enddo
       return
       end
c-----------------------------------------------------------
      subroutine excl0
c calls exciton sampling,clusterin s.r. for exciton energies
      common/clustr/m3set
      common/lscpl/bx(180),sppj
      common/kkg2/fractot(99),cs90
      common/zamax/kzmax,kamax
      common/hjk/jang,iq,rcss
      common/tr/xtran
      common/atsv/atsave
      common/astp3/niso,no,iadlt,isotag,massno(12),mass1(12),ncsave
      common/endvarray/endvcm1(36,999,8),endvcm2(6,36,999,8),
     1endvcm3(10,36,999,8),cmen1(8,999),cmen2(6,8,999),cmen3(10,8,999)
      common/endvparm/iend,lll(6,6,6),ll(6,6)
      common/isomer3/iziso(690),iaiso(690),csgs(690),csis(690)
     1,csgse(690,100),csise(690,100),csis2e(690,100),rm1(690,100),
     2rm2(690,100)
      common/pregam/cff(2,9999),kfn,kfp
      common/temscale/a2h,a3h,a3he,a4he,a4hep
      common/clustr1/eahe4m,eadm,eatm,eahe3m
      common/shft/mc,mp,inver,ike,ipch,kq,qval,ap,at,zp,zt,cld,barfac
      common/history/ide(60),idn(60),idith(60),iddum(60),the(60),adelj
     1(60),spinfh(60),enghis(60)
      common/out1/nout
      common/histindx/ihist
      common/giani/pre, sigt(1100),sigam(999),sigpre(999),gsp(999)
     1,eqgam(999)
      common/par2/cncss
      common/outlim/limout,lim2
      common/lout/spincon
      common/iihist/kp
      common/levop/ldopt,igam
      common/perp/spinprp
      common/arrays/eout(300),efd(40,11),pferm(40,11)
     1,tsetlb(181,300),edist(300),eric(300),ie(40),itl(40),
     2pfpe(40,11),efinc(40,11),eferm2(40,11),phirad(100),peincsq,peinc
       common/tang2/anexc(40),anexci(40),eferm(40,11),psq(40,11),
     1emit(300),angl(181),thetlb(37,300),pdist(300),pd(300)
      common/tang3/efermp,eincc,nexc
      common/nuopt/kinem,iem,kb
      common/excid/npartt(2)
      common/tem/philabb(36),phirec(36),phb(36),phin(36)
      common/bound/neut2,iprot,idnuc(100),idm(100),einc(100),
     1thett(100),phit(100),ppp
      common/hyb2/pp(15,24,999)
      common/nhy/ij,jl,ji,jq,bd,bx1,bx2,cmx,cv,bav,bost,b(3)
      common/evapmc/jcount,ipremc,eminn(15,24)
      common/ilcom/nu,ndum,am(8)
      common/rcl/amis
      common/labcon/prx,pry,prz,pcm,rmas,clab(999,8),dslb(36,999,8)
      common/exclbe/idz(4),ida(4),beff(20,20,2),pairmc(18,27)
      common/scat/jz,ja, probxx(2,2),q1(8,1100),fcs,eloss,fcspec
        common/mccom/ncount,ipot,rum(15,24),ppmc(15,24,999),cf
      common/exclus/ gib(8,999),pr(3),ge(2100),zz(3),
     1delta(3),ts(3),pairx(8),asum(10),rate(8,999)
      common/pl3/jx,nepr,eb(100),rcsp(100),zee,amass,plex
      common/par3/ep,sigml(999),acrs(999)
      common/memo/scrs(999,8)
      common/incr/ed,tened,efinal
      common/sft5/exc(18,27),xmxx
      common/lym1/be(18,27,8),pair(18,27),xmas(18,27),delshl(18,27),amas
     1(18,27)
c     common/bym1/symb(18,27),symbp(18,27),symbs(18,27)
      common/dist/n,ni,ke,irefr,xz(3),pod,pcs,ek,xmix,cfrac
     1,dtau,sign(36),tmix,bji(10),fin,dtup
      common/error/rerror(15,24)
      common/eder/ed2,v,vhalf,vsq
      common/angprm/aprp(500),eav3(1000),
     1eav2(1000),pot2x,thet1,phi1,snt(360),cst(360),snp(720),csp(720)
     2,pot
crc
      common/out2/vrl(200,25),vrlz(200,25),vrlx(200,25),vrly(200,25)
      common/out3/vscale
      common/out4/asm(36),asmm(36)
crc
      common/distb/dsgp(36,999,8)
      common/recoils/rcletza(100,13,24,18),rcleza(100,13,24),
     1rclea(100,36),rcleat(100,36,18),rphir(36)
     2,rclsp(200),rclz(200),rcly(200),rclx(200),escale
     3,rcle(200),rclt(36),rclet(200,36)
      common/rec1/ermxc,thetrc
      common/eferm/efermt
      common/delts/del(81),delt(81)
      common/mmtcom/mmt,mzee,iares,izres,light
      common/fermst/aferm
      write(*,*)' in sr excl0-no ddcs'
      kzmax=zt+zp
      kamax=at+ap-zt-zp
c
c if ipot.gt.0 use infinite hole depth,else finite as default
c next divide number of mc events by partial rx cross sect for gdh
c model
c or should use total for hybrid model
crc
c recoil info added to this version 7/14/97;this will be known as
c 'bastille day alice.we define the following indices for various
c recoil arrays:
c energy index,1-100,ermax=1.25*e-lab
c projectile*mass(proj)/mass(cn)
c so index to store in array = erecoil*100./ermax+1.
c for theta index,use theta-recoil*18./pi+1. to give 10 degree
c intervals
c constants for particle id in exciton scatter routine
c      read(*,*)efermp
c add temp efermp def
c set z,a max index values for light elements
c     m3set=7
       ed2=ed/2.
      kzmax=zt+zp
      kamax=at+ap-zt-zp
      kinem=0
      xmax=xmxx
      xtran=xmax
      ja=1+iadlt
      jz=1
       a2h=0.20
        a3h=0.05
         a3he=0.08
          a4he=0.030*(198./amass)**0.5
c note that in this version the E dependence of a4he is reset in clustpe at each energy
      if(no.eq.1)then
        do i=1,9999
         do j=1,2
          cff(j,i)=0.
           enddo
            enddo
                 endif
       amast=atsave+ap
c      aferm=7.
c      aferm=12.           !  KKG  04/07/10
       if(efermp.eq.0.)efermp=30.
       if(efermt.eq.0.)efermt=30.
       npartt(1)=2
       npartt(2)=1
       mmt=amast+2.
       mzee=zp+zt+1.
       kp=2
       write(17,1112)
1112   format(' monte carlo precompound with pairing correction')
       if(iend.eq.1.and.no.eq.1)call zero(9)
crc
       call initial
        if(no.eq.1)then
          do i=1,690
           csgs(i)=0.
           csis(i)=0.
          enddo
                   endif
c initialize arrays
c define input/starting values for hi e,theta arrays
       eink=xmax/ap
      nexc=ap
      ef=efermp
      efermt=30.
c calculate mid point efermi with fermi gas spacing-needs correction
        nproj=ap
          if(no.eq.1)then
       do710 l=1,nproj
        anexc(l)=(float(l)/ap)**(2./3.)*efermp
  710  continue
       del(1)=anexc(1)
       do 720 l=2,nproj
       del(l)=anexc(l)-anexc(l-1)
       delt(l)=del(l)/10.
       do 720 k=1,11
       efd(l,k)=anexc(l)+float(k-1)*delt(l)
       pferm(l,k)=sqrt(2.*efd(l,k))
       eferm(l,k)=2.*efd(l,k)
  720  continue
       delt(1)=del(1)/10.
       do 721 k=1,11
       efd(1,k)=anexc(1)+float(k-1)*delt(1)
       pferm(1,k)=sqrt(2.*efd(1,k))
       eferm(1,k)=2.*efd(1,k)
  721  continue
        endif
c   here add target fermi energy-temporarily ignore qval
       do 9 je=1,300
       eout(je)=0.
       edist(je)=0.
       eric(je)=0.
c c    pdist(je)=0.
c c    pd(je)=0.
       do 19 jan=1,37
       thetlb(jan,je)=0.
   19  continue
    9  continue
c      ef=30
c define a fermi energy of 30 MeV
       peincsq=2.*(eink+efermt)
       peinc=sqrt(peincsq)
c april change to smooth as above
       do 300 jj=1,nexc
       do 300 k=1,11
       psq(jj,k)=2.*pferm(jj,k)+peinc
       efinc(jj,k)=efd(jj,k)+eink
       pfpe(jj,k)=pferm(jj,k)*peinc
       eferm2(jj,k)=eferm(jj,k)+peincsq
  300  continue
c begin changes for closed form solutions
c
       pot=35.
       pol=pot
       pot2x=pot*2.
       v=pot
       vhalf=v/2.
       vsq=v*v
cxxx
c  define z,a index offset increments for n(=1) and p(n=2)
c  emission
       idz(1)=0
       idz(2)=1
       ida(1)=1
       ida(2)=0
c
       thett(1)=0.
       phit(1)=0.
       einc(1)=xmax
       if(ap.eq.0.)einc(1)=ep
       idnuc(1)=2
       if(zp.eq.0..and.ap.eq.1.)idnuc(1)=1
c
c
       ppp=0.
       do 10103 i=2,100
       einc(i)=0.
       idnuc(i)=0.
10103  continue
c---------------------------
         iexx=(xmax+ed)/ed
         nu=2
         if(zp.eq.0.)nu=1
         cff(nu,iexx)=cff(nu,iexx)+rcss
c--------------------
       if(ap.le.1..and.bd.eq.0.)bd=1.
       if(bd.le.ap.and.ap.gt.1.)bd=ap
       if(ap.gt.1.)bd=ap
       jc=bd
       jj=1
c
c begin mc loop on number of events ncount
c
c renormalize to use only up to 90-h-bar
       cs90=0.
       do 80001 l=1,90
       cs90=cs90+sigml(l)
80001  continue
       fcs=cs90/float(ncsave)
       fcspec=fcs/ed
c--------------------------------------------------
         call spinch
c spinch sums over incident cs by partial wave l
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
       do 10000 icnt=1,ncsave
        call xjlcpl(0)
c  xjlcpl will select entrance ch l value, couple it with proj spin
c for spin 0,1/2,1----need to complete for other projectile spins
       nnn=1
       if(ap.le.1.)nnn=2
c10201  continue
       do 10100 i=nnn,jj+10
       einc(i)=0.
       idnuc(i)=0.
10100  continue
crc
       jj=1
       ii=0
       nout=0
c dec 98 for clusters
       eloss=0.
       jz=1
       ja=1+iadlt
       if(ap.gt.0.)go to 91736
       n1=1
       n11=2
       nhol=1
       if(rndm(-1.).gt.0.5)then
       n1=2
       n11=1
       endif
       idnuc(1)=n1
       idnuc(2)=n11
       if(rndm(-1.).gt.0.5)nhol=2
91736      if(jc.eq.1)go to 91732
       call clusterin(tph,pph)
c this should return,for one event,the initial exciton energies(einc()),and angles
c in radians(angl()-theta,plus exciton id idnuc())
       do 91731 jexc=1,jc
       n1=idnuc(jexc)
       ecall=einc(jexc)
       if(ecall.le.0.1)go to 91731
        jsav=ja
        izsav=jz
       call clustpe(ecall,en1,ipart,epap)
       if(n1.eq.3)go to 91733
       iares=mmt-ja-jz
       izres=mzee-jz
         if(jz.gt.15.or.ja.gt.24)go to 91731
         if(izres.lt.1.or.iares.lt.1.)go to 10200
       en1lab=einc(jexc)-be(jz,ja,n1 )
       if(en1lab.lt.0.)go to 91731
       ix=(en1lab+ed)/ed
       if(rndm(-1.).gt.q1(n1,ix))go to 91731
c  9/20/2014 convert channel energy ejectile to cm:
        ecmejct=(float(iares-1)/float(iares))*en1lab
         ix=(ecmejct+ed)/ed
c--------------------------------------------------
       scrs(ix,n1)=scrs(ix,n1)+fcspec
c--------------------------------------------------
         thet1=0.
         ith=0
       if(iend.eq.1)call endf(thet1,en1lab,ith,ix,n1)
       eloss=eloss+einc(jexc)
       einc(jexc)=0.
       ja=ja+ida(n1)
       jz=jz+idz(n1)
       go to 91731
c check if nucleons from cluster are emitted prior to rescattering
91733  continue
       nhol=1
       ehol=einc(jexc)
       if(rndm(-1.).gt.0.5)nhol=2
       call holscat(nhol,ehol,0.,0.)
       einc(jexc)=0.
91731   continue
91732   continue
       jj=jc+2
crc
c
10200  continue
ccxxxxxxxx looks as tho here, ii will be set to 1
       ii=ii+1
 
       ein=einc(ii)
       if(ein.le.0.)go to 5000
       ipart=idnuc(ii)
c nucleon of energy 'ein' above Fermi energy makes a 2p1h state with
c exciton 'n1' of energy 'en1' and 'n11' of energy 'en11',with a hole
c nhol of energy ehol=ein-en1-en11.If these excitons are energetic
c enough,they are stored to rescatter in this 'loop' if they are not
c emitted.
       call partid(ein,ipart,n1,n11,nhol)
c  partid returns the particle1,2,and hole identities for the
c  incident
c  nucleon at energy ein( as n1,n11,nhol)
c  begin additions for daughter pairing shift,and calc en1
c
       kz=jz+idz(n1)
       ka=ja+ida(n1)
       iares=mmt-ja-jz
       izres=mzee-jz
       icyle=0
        if(kz.gt.13)kz=13
        if(ka.gt.22)ka=22
        if(kz.gt.kzmax .or. ka.gt.kamax)then
          eff=0.
           en1lab=0.
           go to 2001
            end if
c
        eff=ein -pairmc(kz,ka)
       
 635   x=rndm(-1.)
       if(eff.le.v)then
        en1=eff*(1.-sqrt(1.-x))
       go to 40
       end if
       br=(eff-v)/(eff-vhalf)
       if(x.gt.br)then
       en1=eff-sqrt(v*(2.*eff-v)*(1.-x))
       else
       en1=x*(eff-vhalf)
       end if
c
c now en1 is the energy of the 1p of the 2p1h config following rescatter
c of
c the first 1p from a 2p1h config with energy ein  for the primary 1p
c
         ien=(en1+ed)/ed
          if(ien.lt.1)ien=1
           cff(n1,ien)=cff(n1,ien)+fcs
   40  continue
c
       icycle=icycle+1
       if(en1.gt.ein.and.icycle.eq.1)go to 635
       if(en1.gt.ein.and.icycle.eq.2)en1=ein-0.1
       en1lab=en1-be(jz,ja,n1 )
       e1p1h=ein -en1
c
c 2000 is to treat 1p1h left when this n1 cannot be emitted or
c      rescatters
c
           if(en1lab.lt.0.)go to 2000
       ix=(en1lab+ed)/ed
       if(rndm(-1.).gt.q1(n1,ix))go to 2001
c 9/20/2014 change from channel energy to cm for ejectile
        ecmejct=en1lab*(float(iares-1)/float(iares))
           ix=(ecmejct+ed)/ed
c
c here we see if n1 at en1 is emitted or rescatters(go to 2001
c rescat)
c if branch negative,particle rescatters so we go to 2000 for 1p1h
c calc
c-------------------------------------------------------------
       scrs(ix,n1)=scrs(ix,n1)+fcspec
c------------------------------------------------------------
         thet1=0.
         ith=0
       if(iend.eq.1)call endf(thet1,en1lab,ith,ix,n1)
c add recoil/lab conv
       eloss=eloss+en1
       ja=ja+ida(n1)
       jz=jz+idz(n1)
       go to 2000
c this is nucleon with en1
c
c now store parameters of nucleon which rescatters rather than exit
 2001  if(en1lab.le.1.)go to 2000
c only consider rescatter if nucleon exceeds binding energy by a
c constant-here it is 5 MeV
       jj=jj+1
       einc(jj)=en1
       idnuc(jj)=n1
 2000  continue
c here we do 1p1h ;we reset to the angles of the 1p1h in z-frame
c we have saved angle of 1p1h in z dir as thetph,phiph
c first,pick energy of particle n11
       kz=jz+idz(n11)
       ka=ja+ida(n11)
        if(kz.gt.kzmax .or. ka.gt.kamax)then
          eff=0.
           en1lab=0.
           go to 5000
            end if
       u1=ein -en1
       if(u1.le.be(jz,ja,n11))go to 5000
       x=rndm(-1.)
       if(u1.le.v)then
       en11=u1*x
       else
        en11=v*x+u1-v
        end if
c
c 750  continue
       if(en11.gt.u1)en11=u1-0.1
       en11lab=en11-be(jz,ja,n11)
       ehol=u1-en11
            ken=(en11+ed)/ed
            if(ken.lt.1)ken=1
            cff(n11,ken)=cff(n11,ken)+fcs
       if(ehol.gt.beff(jz,ja,nhol))call holscat(nhol,ehol,thol,
     1phol)
       if(en11lab.lt.0.01)go to 5000
c now is this 1p from the second 1p1h going to emit or scatter?
       ix=(en11lab+ed)/ed
       if(rndm(-1.).gt.q1(n11,ix))go to 5001
       ja=ja+ida(n11)
       jz=jz+idz(n11)
c 9/20/2014 change channel energy to cm
       atem=amass+2.-float(ja)-float(jz)
        en11cm=en11lab*(atem-1.)/atem
         ix=(en11cm+ed)/ed
c----------------------------------------------------------
       scrs(ix,n11)=scrs(ix,n11)+fcspec
c----------------------------------------------------------
         thet1=0.
         ith=0
       if(iend.eq.1)call endf(thet1,en1lab,ith,ix,n1)
       eloss=eloss+en11
       go to 5000
 5001  if(en11lab.le..01)go to 5000
       jj=jj+1
       einc(jj)=en11
       idnuc(jj)=n11
 5000  continue
c now test for hot nucleons,and return to top if found
       if(ii.lt.jj)go to 10200
 
c  ---------------------------------------
c now do storage
c  ----------------------------------
       efinal=xmax-eloss
       if(efinal.lt.0.)efinal=0.
       ixm=(efinal+ed)/ed
c MAR 970414
c add call to evap subroutine for MC evaporation
c skip ppmc storage but increment pp array if eloss is large enough
c to prevent evaporation
       if(jz.le.iprot.and.ja.le.neut2)go to 10005
        ppp=ppp+fcs
       go to 10000
10005  continue
       iares=mmt-ja-jz
       izres=mzee-jz
       enn=efinal
       if(efinal.lt.eminn(jz,ja)) then
       if(iares.lt.aferm)then
       call gudima(enn,kst)
       call fermfrag(kst)
       light=1
        go to 937
       endif
c----------------------------------------------------
         pp(jz,ja,ixm)=pp(jz,ja,ixm)+fcs
c----------------------------------------------------
       endif
       xtran=xmax
  937  continue
       if(efinal.ge.eminn(jz,ja))call evap
 
            if(iend.eq.0)go to 10000
          if(nout.gt.0.and.nout.le.3)call endvsort1
10000  continue
crc
       if(no.eq.niso)then
        call sortmc
       if(iend.gt.0)call endvout
                     endif
       return
       end
c-----------------------------------------------------------------
      subroutine excl1
      common/clustr/m3set
      common/kkg2/fractot(99),cs90
      common/clkin/elab
      common/dec/reccon
      common/exci/excdist(700),delexc,counter,ctloop,excdist1(700)
      common/atsv/atsave
      common/shft/mc,mp,inver,ike,ipch,kq,qval,ap,at,zp,zt,cld,barfac
      common/astp3/niso,no,iadlt,isotag,massno(12),mass1(12),ncsave
      common/endvarray/endvcm1(36,999,8),endvcm2(6,36,999,8),
     1endvcm3(10,36,999,8),cmen1(8,999),cmen2(6,8,999),cmen3(10,8,999)
      common/endvparm/iend,lll(6,6,6),ll(6,6)
      common/hjk/jang,iq,rcss
      common/pregam/cff(2,9999),kfn,kfp
      common/giani/pre, sigt(1100),sigam(999),sigpre(999),gsp(999)
     1,eqgam(999)
      common/hybsub/iav,idum,rs(8,1100),c,pox,r1
      common/temscalf/b2h,b3h,b3he,b4he,b4hep
      common/temscale/a2h,a3h,a3he,a4he,a4hep
      common/tr/xtran
      common/fisn/dfslb(36,999,4),cflab(999,4),dfetl(36,999,4)
      common/ismrdat/isoid(106,260),levnm(1000),bevj(1000,3)
      common/gams/gamgam,csgam,csgamm(100)
      REAL*8 gamgam,csgamm,csgam
      common/trans/itran
      common/outlim/limout,lim2
      common/iihist/kp
      common/histindx/ihist
      common/levop/ldopt,igam
      common/ismrdat2/ididl(106,260),nolevs(106,260)
      common/perp/spinprp
      common/isomer1/xgs(690),xiso(690),x2iso(690),isonum,isoint,iso2num
      common/isomer2/zs(690),as(690),eiso(690),eiso2(690)
      common/isomer3/iziso(690),iaiso(690),csgs(690),csis(690)
     1,csgse(690,100),csise(690,100),csis2e(690,100),rm1(690,100),
     2rm2(690,100)
      common/lin/lspin
      common/lscpl/bx(180),sppj
      common/ijcpl/ajtarg,ajin
      common/tagspin/trgspin
      common/spins2/spinu,spinf,prz1
      common/tgspin/levno
      common/rotor/irotor,iidummy,arote(15,25,181)
      common/history/ide(60),idn(60),idith(60),iddum(60),the(60),adelj
     1(60),spinfh(60),enghis(60)
       common/spins4/eprot,ealpha,eoff(8)
      common/spins/rote(15,24,181),popej(15,24,180,100),fcx(180),ai,ecmp
      common/sping/pfpej(600,180,100),pfpj(600,180)
      common/gate1/spind(5,180),popj(15,24,180),egatl(5),egatu(5),npart(
     15),igate
      common/nuopt/kinem,iem,kb
      common/excid/npartt(2)
      common/tem/philabb(36),phirec(36),phb(36),phin(36)
      common/bound/neut2,iprot,idnuc(100),idm(100),einc(100),
     1thett(100),phit(100),ppp
      common/hyb2/pp(15,24,999)
      common/evapmc/jcount,ipremc,eminn(15,24)
      common/ilcom/nu,ndum,am(8)
      common/rcl/amis
      common/labcon/prx,pry,prz,pcm,rmas,clab(999,8),dslb(36,999,8)
      common/exclbe/idz(4),ida(4),beff(20,20,2),pairmc(18,27)
      common/scat/jz,ja, probxx(2,2),q1(8,1100),fcs,eloss,fcspec
        common/mccom/ncount,ipot,rum(15,24),ppmc(15,24,999),cf
      common/exclus/ gib(8,999),pr(3),ge(2100),zz(3),
     1delta(3),ts(3),pairx(8),asum(10),rate(8,999)
      common/pl3/jx,nepr,eb(100),rcsp(100),zee,amass,plex
      common/par3/ep,sigml(999),acrs(999)
      common/memo/scrs(999,8)
      common/incr/ed,tened,efinal
      common/sft5/exc(18,27),xmxx
      common/lym1/be(18,27,8),pair(18,27),xmas(18,27),delshl(18,27),amas
     1(18,27)
c     common/bym1/symb(18,27),symbp(18,27),symbs(18,27)
      common/distb/dsgp(36,999,8)
      common/dist/n,ni,ke,irefr,xz(3),pod,pcs,ek,xmix,cfrac
     1,dtau,sign(36),tmix,bji(10),fin,dtup
      common/error/rerror(15,24)
      common/eder/ed2,v,vhalf,vsq
      common/angprm/aprp(500),eav3(1000),
     1eav2(1000),pot2x,thet1,phi1,snt(360),cst(360),snp(720),csp(720)
     2,pot
crc
      common/out1/nout
      common/out2/vrl(200,25),vrlz(200,25),vrlx(200,25),vrly(200,25)
      common/out3/vscale
      common/out4/asm(36),asmm(36)
crc
      common/recoils/rcletza(100,13,24,18),rcleza(100,13,24),
     1rclea(100,36),rcleat(100,36,18),rphir(36)
     2,rclsp(200),rclz(200),rcly(200),rclx(200),escale
     3,rcle(200),rclt(36),rclet(200,36)
      common/rec1/ermxc,thetrc
      common/mmtcom/mmt,mzee,iares,izres,light
      common/tempo/sp3of(999),sp3(999)
      common/clustr1/eahe4m,eadm,eatm,eahe3m
      common/nuc/cnu
      common/fermst/aferm
      common/deljay/delj
      common/misscs/ppstor
c
c if ipot.gt.0 use infinite hole depth,elso finite as default
c next divide number of mc events by partial rx cross sect for gdh
c model
c or should use total for hybrid model
c
crc
c recoil info added to this version 7/14/97;this will be known as
c 'bastille day alice.we define the following indices for various
c recoil arrays:
c energy index,1-100,ermax=1.25*e-lab
c projectile*mass(proj)/mass(cn)
c so index to store in array = erecoil*100./ermax+1.
c for theta index,use theta-recoil*18./pi+1. to give 10 degree
c intervals
c temporarily specify 'aferm',the mass no. at which fermi stats are used
c a2= precompound deuteron scale factor; a4= alpha scale factor;
c a3= precompound triton scale factor; a5= precompound 3He scale factor
c to use read for these constants, the values must be commented out here,
c and the read statement in MAIN must be activated
      ed2=ed/2.
      a2h=0.20
       a3h=0.05
        a3he=0.08
         a4he=0.030*(198./amass)**0.5
c note that in this version the E dependence of a4he is reset in clustpe at each energy
c primary or secondary.
c---------------------------------------------
       write(*,998)sppj
 998   format(' projectile spin = ',1f5.2)
       call lscupl
c get fraction of l-0.5 couple vs l+/- couple
c--------------------------------------------
      delexc=1./float(ncount)
c delexc=element of reaction cross section for each MC event
          if(no.eq.1)then
         do i=1,9999
         do j=1,2
         cff(j,i)=0.
          enddo
          enddo
                    endif
c set aferm=1. to turn off Fermi breakup
       aferm=12.           !  KKG  04/07/10
       pcon=36./3.1416
       amast=atsave+ap
       amis=amast+2.
       ztem=zee
       spinprp=0.
       ajtarg=trgspin
       thetrc=18./3.1416
       epap=sqrt(ep*ap)
       ermax=4.00*ep/(atsave+ap)
       ermxc=100./ermax
       escale=amast*50./ep
       vscale=50.*amast/(sqrt(2.*ap*ep))
       xmax=xmxx
       write(*,999)xmax
 999   format(' xmax = ',f6.2)
c-----------------m-
        jfis=0
        kp=2
        jfis=0
c the above are index constants for recoil arrays
crc
c
c for recoil
        przz=sqrt(2.*ep)
        prz1=przz
        prx=0.
        pry=0.
        prz=prz1
              prxi=0.
              pryi=0.
              przi=przz
c add def prz 11-19-2009
c constants for particle id in exciton scatter routine
       npartt(1)=2
       npartt(2)=1
cc
      ja=1+iadlt
      jz=1
       write(17,1112)
       mmt=amast+2.
       mzee=zp+zt+1.
       fcs=cf/float(ncsave)
       fcspec=fcs/ed
c
       pot=35.
       pol=pot
       pot2x=pot*2.
       v=pot
       vhalf=v/2.
       vsq=v*v
c initialize arrays
c------------------------
       call initial
       call trigstuf
c------------------------
c  define z,a index offset increments for n(=1) and p(n=2)
c emission
       idz(1)=0
       idz(2)=1
       ida(1)=1
       ida(2)=0
c
       thett(1)=0.
       phit(1)=0.
       idnuc(1)=2
       if(zp.eq.0..and.ap.eq.1.)idnuc(1)=1
c
       ppp=0.
         iexx=(xmax+ed)/ed
        nu=2
        if (zp.eq.0.)nu=1
        cff(nu,iexx)=cff(nu,iexx)+rcss
c
c july 2001 begin loop on composite entrance ch I
       call spinch
       jlim=ai+1.
c--------------------------------------------------------
       do 10000 icnt=1,ncount     
c-------------------------------------------------------
       light=0
       spinprp=0.
c---------------------------------------------------
         call xjlcpl(jzl)
c         xjlcpl returns random l-inc coupled with projectile spin as
c         'spinu';jzl is l-inc before spin couple
         call ijlcpl
c--------------------------------------------------
c spinu, incident wave,returned from xjlcpl is l.s wave coupled to I-trg---------------
          js=2.*spinu+1.
          if(js.lt.1)js=1
         xmax=xmxx-arote(1,1+iadlt,js)
c-------------------------------------------------------
        jj=1
c-----------------------------
       call exclzro
c-----------------------------
       continue
       einc(1)=xmax
       idnuc(1)=2
       if(zp.eq.0..and.ap.eq.1.)idnuc(1)=1
crc
       ii=0
       jj=1
c for recoil
c here set initial event recoil nucleus momenta x,y,z
       prz=przz
       prx=0.
       pry=0.
       eloss=0.
       jz=1
       ja=1+iadlt
c
       ii=ii+1
       ein=einc(ii)
       if(ein.le.ed/2.)go to 5000
c      if(ein.le.0.01)go to 5000
       thet1=thett(ii)
       phi1=phit(ii)
       ipart=idnuc(ii)
       ko=0
c 9/2/06 here start cluster pe by calculating pe cluster width,
c selecting whether nucleon scatter or cluster for this event
       rmas=amass+2.-float(ja)-float(jz)
       if(rmas.le.4.)go to 10000
       if(rmas.ge.8.)then
       if(m3set.gt.0)call clustpe(ein,en1,ipart,epap)
         endif
c nucleon of energy 'ein' above Fermi energy makes a 2p1h state with
c exciton 'n1' of energy 'en1' and 'n11' of energy 'en11',with a hole
c nhol of energy ehol=ein-en1-en11.If these excitons are energetic
c enough,they are stored to rescatter in this 'loop' if they are not
c emitted.
       if(ein.le.ed/2.)go to 5000
c         if(ein.le.0.01)go to 5000
       call partid(ein,ipart,n1,n11,nhol)
c  partid returns the particle1,2,and hole identities for the
c incident
c  nucleon at energy ein( as n1,n11,nhol)
c  begin additions for daughter pairing shift,and calc en1
c
       kz=jz+idz(n1)
       ka=ja+ida(n1)
       iares=mmt-ja-jz
       izres=mzee-jz

       if(light.gt.0)go to 5000
       if(kz.gt.13)kz=13
       if(ka.gt.22)ka=22
       icycle=0
       
        eff=ein -pairmc(kz,ka)
       x=rndm(-1.)
       if(eff.le.v)then
        en1=eff*(1.-sqrt(1.-x))
       go to 40
       end if
       br=(eff-v)/(eff-vhalf)
       if(x.gt.br)then
       en1=eff-sqrt(v*(2.*eff-v)*(1.-x))
       else
       en1=x*(eff-vhalf)
       end if
c now en1 is the energy of the 1p of the 2p1h config following rescatter
c of
c the first 1p from a 2p1h config with energy ein  for the primary 1p
c load fcs into cff() array as exciton population to use for precomp gammas
       if(en1.gt.ein)en1=ein-0.01
        ien=(en1+ed)/ed
        if(ien.lt.1)ien=1
        cff(n1,ien)=cff(n1,ien)+fcs
   40  continue
c
       en1lab=en1-be(jz,ja,n1 )
c
c here incident nucleon scatters to make 2p1h, and we consider
c 1p of 2p1h result.
        call rotangl(p3,ein ,en1,thetph,phiph,en1lab)
c here thetph and phiph are angles of ph in original z direction,
c and thet1,phi1 have been returned in common as the 1p directions
c following the scattering creating the 2p1h state.
c
c  this rot gets us the scattered nucleon angle in z dir even if
c  nucleon
c cannot be emitted-so we can get the angle of the 1p1h remaining
c this returns new thet1,phi1 for n1 with energy en1,and angles for
c the
c particle-hole pair thetph,phiph from the incident nucleon with ein
c energy
c 2000 is to treat 1p1h left when this n1 cannot be emitted or
crescatters
c
       if(en1lab.le.0.00)go to 2000
       ix=(en1lab+ed)/ed
       if(rndm(-1.).gt.q1(n1,ix))go to 2001
c here we see if n1 at en1 is emitted or rescatters(go to 2001
c rescat)
c if branch negative,particle rescatters so we go to 2000 for 1p1h
c calc
c set theta/5 deg bin index to store  scattered cross section
       ith=thet1*11.459155+1.
       if(ith.eq.37)ith=36
c add recoil/lab conv
       eloss=eloss+en1
       ja=ja+ida(n1)
       jz=jz+idz(n1)
       rmas=amass+2.-float(ja)-float(jz)
c      if(rmas.lt.1.)go to 10000
       acm=en1lab*rmas/(rmas+1.)
         ix=(acm+ed)/ed
c   above change index from channel energy to particle cm energy  9/20/2014
c-----------------------------------------------
       scrs(ix,n1)=scrs(ix,n1)+fcspec
       dsgp(ith,ix,n1)=dsgp(ith,ix,n1)+fcspec
c------------------------------------------------
       if(iend.eq.1)call endf(thet1,en1lab,ith,ix,n1)
c--------------------------
       pcm=sqrt(2.*acm)
       qcm=sqrt(2.*en1lab)
c add spin change algos july 2001
       delj=(float((jzl-1)))*(qcm/prz1)*cos(thet1)
c add perp 1/20/02
       delprp=delj*sin(thet1)*sin(phi1)/cos(thet1)
       spinprp=spinprp+delprp
       adelj(ii)=delj
       spinf=spinf+delj
       nu=n1
         ievap=0
c-------------------------------
       call recoil(12)
c------------------------------
c we have thet1,phi1 of emitted particle in z-dir; we need to get
c a lab angle and adjust the recoil vectors for the emission, get
c a lab angle. this is done in recoil sr.
c here thet1,phi1 are angles in z dir of nucleon before rescattering  in
c rescat
c this is nucleon with en1
c now store parameters of nucleon which rescatters rather than exit
 2001  if(en1lab.le.0.00)go to 2000
c only consider rescatter if nucleon exceeds binding energy by a
c constant-here it is 1 MeV
       jj=jj+1
       if(jj.gt.99)jj=99
       einc(jj)=en1
       thett(jj)=thet1
       phit(jj)=phi1
       idnuc(jj)=n1
 2000  continue
       thet1=thetph
       phi1=phiph
c here we do 1p1h ;we reset to the angles of the 1p1h in z-frame
c we have saved angle of 1p1h in z dir as thetph,phiph
c first,pick energy of particle n11
c u1 is excitation shared by 1p+1h
       u1=ein -en1
c--------------------------------------------
      if(u1.le.be(jz,ja,n11))go to 5000
c      if(u1.le.eminn(jz,ja))go to 5000
c-----------------------------------------
       kz=jz+idz(n11)
       ka=ja+ida(n11)
       iares=mmt-ja-jz
       izres=mzee-jz
       if(kz.gt.13)kz=13
       if(ka.gt.22)ka=22
       x=rndm(-1.)
       if(u1.le.v)then
       en11=u1*x
       else
        en11=v*x+u1-v
        end if
      continue
c
       continue
       en11lab=en11-be(jz,ja,n11)
       ehol=u1-en11
c store 2nd nucleon for gamma cascade
         ken=(en11+ed)/ed
         if(ken.lt.1)ken=1
         cff(n11,ken)=cff(n11,ken)+fcs
c here we must compute angle of particle from 1p1h in z dir,then of the
c hole
c i.e. pick angle of 1p from 1p1h,leaving a 1h residue
c call rotang2 to get the 1p angles,then do hole value after return
       call rotang2(p3,u1,en11,thol,phol,en11lab)
       if(ehol.gt.beff(jz,ja,nhol))call holscat1(nhol,ehol,thol,
     1phol)
        
       if(light.gt.0)go to 5000
       if(en11lab.le.0.0)go to 5000
c now is this 1p from the second 1p1h going to emit or scatter?
       ix=(en11lab+ed)/ed
c 7/3/14
         if(rndm(-1.).gt.q1(n11,ix))go to 5001
       ja=ja+ida(n11)
       jz=jz+idz(n11)
       ith=thet1*11.459155+1.
       if(ith.eq.37)ith=36
c add recoil/lab conv----------------
c 9/20/2014 change channel energy to cm
       nu=n11
       rmas=amass+2.-float(ja)-float(jz)
        mzee=zt+zp+1.-float(jz)
       ecm=en11lab*rmas/(rmas+1.)
         ix=(ecm+ed)/ed
c------------------------------------
       pcm=sqrt(2.*ecm)
       qcm=sqrt(2.*en11lab)
c-------------------------------------------------
       scrs(ix,n11)=scrs(ix,n11)+fcspec
       dsgp(ith,ix,n11)=dsgp(ith,ix,n11)+fcspec
c--------------------------------------------------
       eloss=eloss+en11
       if(iend.eq.1)call endf(thet1,en11lab,ith,ix,n11)
c----------------------------
c end spin transfer this emission
c add spin change algos july 2001
       delj=(float((jzl-1)))*(qcm/prz1)*cos(thet1)
       delprp=delj*sin(thet1)*sin(phi1)/cos(thet1)
       spinprp=spinprp+delprp
       adelj(ii)=delj
c      xq=cos(thet1)
c      if(nout.ne.1)go to 72001
c      if(igate.eq.0)go to 72001
c      jfin=(abs(spinu-delj)+1.)
c      do 72002 jgate=1,igate
c      if(n11.ne.npart(jgate))go to 72002
c      if(en11lab.lt.egatl(jgate).or.en11lab.gt.egatu(jgate))go to 72003
c      spind(jgate,jfin)=spind(jgate,jfin)+fcs
c      go to 72001
c72003  continue
c72002  continue
c72001  continue
       spinf=spinf+delj
c end spin transfer this emission
       nu=n11
           efin=xmax-eloss
        call recoil(13)
       go to 5000
c here,if it is to rescat,thet1,phi1 are its initial direction
c in the zyx frame
 5001  if(en11lab.le.0.01)go to 5000
       jj=jj+1
       if(jj.gt.99)jj=99
       einc(jj)=en11
       thett(jj)=thet1
       phit(jj)=phi1
       idnuc(jj)=n11
c----------------------------------
 5000  continue
c----------------------------------
c here calc residual exc, and if.lt.be lowest,store residual cs!
c now test for hot nucleons,and return to top if found
       rmt=amass+2.-float(ja)-float(jz)
       spinsv=spinf
       ispin=(abs(2.*(spinu-spinf))+1.)
       eftem=xmax-eloss
         if(ispin.lt.1)ispin=1
c      if(arote(jz,ja,ispin).gt.eftem.and.mmt.gt.20)go to 10205
c temp c 8/12
c      if(igate.eq.0)go to 92001
c      do 92002 jgate=1,igate
c      if(n11.ne.npart(jgate))go to 92002
c      if(enlab.lt.egatl(jgate).or.enlab.gt.egatu(jgate))go to 92003
c      jfin=(abs(spinu-adelj(ii))+1.)
c      spind(jgate,jfin)=spind(jgate,jfin)+fcs
c note if gates reinstated, the event # ii needs to be reinstated
c      go to 92001
c end spin transfer this emission
       continue
c if the above event is below yrast line after all pe emission,restart the event
       continue
c  ---------------------------------------
c now do storage
c  ----------------------------------
       efinal=xmax-eloss
       ixm=(efinal+ed)/ed
          if(ixm.lt.1)ixm=1
c      spinf=spinsv
c MAR 970414
c add call to evap subroutine for MC evaporation
c skip ppmc storage but increment pp array if eloss is large enough
c to prevent evaporation
cbound
       if(jz.gt.13)jz=13
       if(ja.gt.23)ja=23
       if(jz.le.iprot.and.ja.le.neut2)go to 10005
        if(light.gt.0)go to 10005
        ppp=ppp+fcs
10005  continue
       if(light.eq.0)then
c===============================================================
         ppstor=0.
       if(efinal.le.eminn(jz,ja)) then
c-------------------------------------------------------
         pp(jz,ja,ixm)=pp(jz,ja,ixm)+fcs
c-------------------------------------------------------
         ppstor=1.
        spinxt=spinu-spinf
        spinr=sqrt(spinxt**2+spinprp**2)
         ispinr=abs(spinr+0.5)
          spin=ispinr
      
      ijz=zee+1.-float(jz)
      ija=amass+2.-float(jz)-float(ja)
      ija=atsave+ap+2.-float(jz)-float(ja)
c---------------------------------
      if(nolevs(ijz,ija).eq.1)then
      x=rndm(-1.)
      if(x.lt.0.5)then
       spin=abs(spin-0.5)
       else
       spin=abs(spin+0.5)
       endif
       endif
c-------------------------------
          jl=(2.*spin)+1.
        if(jl.lt.1)jl=abs(jl)
        if(jl.gt.180)jl=180
c-----------------------------------
        if(ixm.le.100)then
        popej(jz,ja,jl,ixm)=popej(jz,ja,jl,ixm)+fcs
        popj(jz,ja,jl)=popj(jz,ja,jl)+fcs
        endif
c----------------------------------
c now branch if an isomer exists
c      if(isno.gt.0.and.isoint.eq.0)call isoevent(efinal,jl)
       endif
c============================================================
 
        xtran=xmax
c---------------------------------------------------------
       if(efinal.gt.eminn(jz,ja))call evap
        
c---------------------------------------------------------
       endif
c call ENDF sort/store endvsort
        if(iend.eq.0)go to 7001
       if (nout.gt.0.and.nout.le.3)call endvsort       
 7001  continue
 
c add recoil energy storage
       call reclsv
10000  continue
       if(itran.lt.1)itran=1
       csgamm(nepr)=pp(1,1,itran)
c--------------------------------------------
      if(no.eq.niso)then
c-------------------------------------------
      call sortmc
      call bngwrt
      if(iend.gt.0)call endvout
                    endif
      call precomgam
      pre=0.
      call gama
c-------------------------------------------
      pre=1.
1112   format(' monte carlo precompound with pairing correction')
      return
      end
c------------------------------------------------------------- 
      subroutine excl3
c this sr is called only for incident photons
      common/clustr/m3set
      common/lscpl/bx(180),sppj
      common/kkg2/fractot(99),cs90
      common/gamcl/igamcl
      common/temscale/a2h,a3h,a3he,a4he,a4hep
      common/clustr1/eahe4m,eadm,eatm,eahe3m
      common/tr/xtran
      common/atsv/atsave
      common/shft/mc,mp,inver,ike,ipch,kq,qval,ap,at,zp,zt,cld,barfac
      common/astp3/niso,no,iadlt,isotag,massno(12),mass1(12),ncsave
      common/endvparm/iend,lll(6,6,6),ll(6,6)
      common/history/ide(60),idn(60),idith(60),iddum(60),the(60),adelj
     1(60),spinfh(60),enghis(60)
      common/out1/nout
      common/histindx/ihist
      common/giani/pre, sigt(1100),sigam(999),sigpre(999),gsp(999)
     1,eqgam(999)
      common/gams/gamgam,csgam,csgamm(100)
      REAL*8 gamgam,csgamm,csgam
      common/trans/itran
      common/ismrdat/isoid(106,260),levnm(1000),bevj(1000,3)
      common/outlim/limout,lim2
      common/iihist/kp
      common/fc/fc(3),jlim
      common/lout/spincon
      common/lab3/sig(8,999),slav(8,999)
      common/phot/bspin(3)
      common/hjk/jang,iq,rcss
      common/levop/ldopt,igam
      common/ismrdat2/ididl(106,260),nolevs(106,260)
      common/isomer1/xgs(690),xiso(690),x2iso(690),isonum,isoint,iso2num
      common/isomer2/zs(690),as(690),eiso(690),eiso2(690)
      common/isomer3/iziso(690),iaiso(690),csgs(690),csis(690)
     1,csgse(690,100),csise(690,100),csis2e(690,100),rm1(690,100),
     2rm2(690,100)
      common/perp/spinprp
      common/lin/lspin
      common/ijcpl/ajtarg,ajin
      common/tagspin/trgspin
      common/tgspin/levno
      common/spins2/spinu,spinf,prz1
      common/rotor/irotor,iidummy,arote(15,25,181)
       common/spins4/eprot,ealpha,eoff(8)
      common/spins/rote(15,24,181),popej(15,24,180,100),fcx(180),ai,ecmp
      common/sping/pfpej(600,180,100),pfpj(600,180)
      common/gate1/spind(5,180),popj(15,24,180),egatl(5),egatu(5),npart(
     15),igate
      common/nuopt/kinem,iem,kb
      common/excid/npartt(2)
      common/tem/philabb(36),phirec(36),phb(36),phin(36)
      common/bound/neut2,iprot,idnuc(100),idm(100),einc(100),
     1thett(100),phit(100),ppp
      common/hyb2/pp(15,24,999)
      common/evapmc/jcount,ipremc,eminn(15,24)
      common/ilcom/nu,ndum,am(8)
      common/rcl/amis
      common/labcon/prx,pry,prz,pcm,rmas,clab(999,8),dslb(36,999,8)
      common/exclbe/idz(4),ida(4),beff(20,20,2),pairmc(18,27)
      common/scat/jz,ja, probxx(2,2),q1(8,1100),fcs,eloss,fcspec
        common/mccom/ncount,ipot,rum(15,24),ppmc(15,24,999),cf
      common/exclus/ gib(8,999),pr(3),ge(2100),zz(3),
     1delta(3),ts(3),pairx(8),asum(10),rate(8,999)
      common/pl3/jx,nepr,eb(100),rcsp(100),zee,amass,plex
      common/par3/ep,sigml(999),acrs(999)
      common/memo/scrs(999,8)
      common/incr/ed,tened,efinal
      common/sft5/exc(18,27),xmxx
      common/lym1/be(18,27,8),pair(18,27),xmas(18,27),delshl(18,27),amas
     1(18,27)
c     common/bym1/symb(18,27),symbp(18,27),symbs(18,27)
      common/distb/dsgp(36,999,8)
      common/dist/n,ni,ke,irefr,xz(3),pod,pcs,ek,xmix,cfrac
     1,dtau,sign(36),tmix,bji(10),fin,dtup
      common/error/rerror(15,24)
      common/eder/ed2,v,vhalf,vsq
      common/angprm/aprp(500),eav3(1000),
     1eav2(1000),pot2x,thet1,phi1,snt(360),cst(360),snp(720),csp(720)
     2,pot
crc
      common/out2/vrl(200,25),vrlz(200,25),vrlx(200,25),vrly(200,25)
      common/out3/vscale
      common/out4/asm(36),asmm(36)
crc
      common/recoils/rcletza(100,13,24,18),rcleza(100,13,24),
     1rclea(100,36),rcleat(100,36,18),rphir(36)
     2,rclsp(200),rclz(200),rcly(200),rclx(200),escale
     3,rcle(200),rclt(36),rclet(200,36)
      common/rec1/ermxc,thetrc
      common/mmtcom/mmt,mzee,iares,izres,light
      common/fermst/aferm
      common/pregam/cff(2,9999),kfn,kfp
      common/testsp/numspn(30)
       do jjj=1,30
        numspn(jjj)=0
         enddo
        ed2=ed/2.
c
c if ipot.gt.0 use infinite hole depth,elso finite as default
c next divide number of mc events by partial rx cross sect for gdh
c model
c or should use total for hybrid model
c
crc
c recoil info added to this version 7/14/97;this will be known as
c 'bastille day alice.we define the following indices for various
c recoil arrays:
c energy index,1-100,ermax=1.25*e-lab
c projectile*mass(proj)/mass(cn)
c so index to store in array = erecoil*100./ermax+1.
c for theta index,use theta-recoil*18./pi+1. to give 10 degree
c intervals
c temporarily specify 'aferm',the mass no. at which fermi stats are used
c sept 02 make l-ejectile=0.23*sqrt ep,where r-ave=1.5A(1/3)fm/1.414
       write(17,1112)
       xmax=xmxx
       amast=atsave
       ja=1+iadlt
       jz=1
       a2h=0.20
       a3h=0.05
       a3he=0.08
       a4he=0.030*(198./amass)**0.5
        sppj=1.
        if(no.eq.1)then
        do i=1,9999 
        do j=1,2
        cff(j,i)=0.
        enddo
        enddo
                 endif
       spincon=0.33*amast**(1./3.)
       aferm=12.           !  KKG  04/07/10
       spinprp=0.
       ajtarg=trgspin
        call photspin
       pcon=36./3.1416
       thetrc=18./3.1416
       ermax=4.00*ep/(atsave+ap)
       ermxc=100./ermax
       escale=amass*50./ep
       amis=amass+2.
       vscale=0.04242*sqrt(ep)/amass
c      vscale=50.*amass/(sqrt(2.*ap*ep))
c the above are index constants for recoil arrays
crc
c
       erecoil=0.00053*ep*ep/amast
       przz=sqrt(erecoil*amast*2.)
       epap=przz/1.414
       prz1=przz
       prx=0.
       pry=0.
c constants for particle id in exciton scatter routine
       npartt(1)=2
       npartt(2)=1
cc
c change to number=ncount for photons
       number=ncount
        if(iend.gt.0.and.no.eq.1)call zero(9)
       call initial
c initialize arrays
c identify nuclides having isomers
       do 9002 i=1,690
       csgs(i)=0.
       csis(i)=0.
9002   continue
c identify by z,a index those nuclides having an isomer
       mmt=amast+2.
       mzee=zp+zt+1.
       numb=ncsave
       if(niso.gt.1)numb=cfrac*float(ncsave)
       cf=rcss
c change 1/20/08 numb to ncsave      fcs=cf /(float(numb))
       fcs=cf /(float(ncsave))
 
       fcspec=fcs/ed
c begin changes for closed form solutions
c
       pot=35.
       pol=pot
       pot2x=pot*2.
       v=pot
       vhalf=v/2.
       vsq=v*v
       call trigstuf
cxxx
c  define z,a index offset increments for n(=1) and p(n=2)
c emission
       idz(1)=0
       idz(2)=1
       ida(1)=1
       ida(2)=0
c
       thett(1)=0.
       phit(1)=0.
       einc(1)=xmax
c add 7 lines june 03
       if(xmax.le.eminn(1,1) )then
        ixc=(xmax+ed)/ed
         if(ixc.lt.1)ixc=1
        pp(1,1,ixc)=rcss
       csgamm(nepr)=pp(1,1,ixc)
       return
                              endif
       idnuc(1)=2
       if(zp.eq.0..and.ap.eq.1.)idnuc(1)=1
c
c
       ppp=0.
       do 10103 i=2,100
       einc(i)=0.
       thett(i)=0.
       phit(i)=0.
       idnuc(i)=0.
10103  continue
c
c
       jz=1
       ja=1+iadlt
cxx
        kp=0
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
       do 10000 icnt=1,number
        call photonj
c photonj gets target spin/photon spin coupled spinu result
        item=2.*spinu
        numspn(item)=numspn(item)+1
           
       spinf=0.
       spinprp=0.
       light=0
       kp=1
        jj=1
10205  continue
c set switch 'light' to tell which light residues will undergo breakup
c with light=1 meaning breakup algo has been applied,so no residue to store
       do 10202 inum=1,60
       ide(inum)=0
       idn(inum)=0
       idith(inum)=0
       the(inum)=0.
       adelj(inum)=0.
       spinfh(inum)=0.
       enghis(inum)=0.
10202  continue
cxx
       delj=0.
       ihist=0
       light=0
       nout=0
c10201  continue
       do 10100 i=2,jj+1
       einc(i)=0.
       thett(i)=0.
       phit(i)=0.
       idnuc(i)=0.
10100  continue
crc
       ii=0
       jj=1
       nout=0
crc
c for recoil
       prz=przz
       prx=0.
       pry=0.
       eloss=0.
       jz=1
       ja=1+iadlt
c the above  calls  a random no. from 0. to 1.
c
c
10200  continue
       ii=ii+1
       ein=einc(ii)
       thet1=thett(ii)
       phi1=phit(ii)
       ipart=idnuc(ii)
       ko=0
       n1=1
       n11=2
       if(rndm(-1.).gt.0.5)then
       n1=2
       n11=1
       endif
       nhol=1
       if(rndm(-1.).gt.0.5)nhol=2
c 1/20/08 add call to cluster emission precompound assuming initial
c photon acts as a nucleon- else we later will do only for excited
c nucleons.
       rmas=amass+2.-float(ja)-float(jz)
c
       if(m3set.gt.0.and.ii.ge.igamcl)call clustpe(ein,en1,ipart,epap)
c
c nucleon of energy 'ein' above Fermi energy makes a 2p1h state with
c exciton 'n1' of energy 'en1' and 'n11' of energy 'en11',with a hole
c nhol of energy ehol=ein-en1-en11.If these excitons are energetic
c enough,they are stored to rescatter in this 'loop' if they are not
c emitted.
c      call partid(ein,ipart,n1,n11,nhol)
c replace partid with definitions for photon in
c  begin additions for daughter pairing shift,and calc en1
c
          if(ein.lt.0.01)go to 6000
       kz=jz+idz(n1)
       ka=ja+ida(n1)
       iares=mmt-ja-jz
       izres=mzee-jz
       if(kz.gt.13)kz=13
       if(ka.gt.22)ka=22
       icycle=0
        eff=ein -pairmc(kz,ka)
       x=rndm(-1.)
       if(eff.le.v)then
        en1=eff*(1.-sqrt(1.-x))
       go to 40
       end if
       br=(eff-v)/(eff-vhalf)
       if(x.gt.br)then
       en1=eff-sqrt(v*(2.*eff-v)*(1.-x))
       else
       en1=x*(eff-vhalf)
       end if
c
c now en1 is the energy of the 1p of the 2p1h config following rescatter
c of
c the first 1p from a 2p1h config with energy ein  for the primary 1p
c
       ien=(en1+ed)/ed
        if(ien.lt.1)ien=1
         cff(n1,ien)=cff(n1,ien)+fcs
   40  continue
c
c 740  continue
       if(en1.gt.ein)en1=ein-0.01
       en1lab=en1-be(jz,ja,n1 )
c      e1p1h=ein -en1
c
        call rotangl(p3,ein ,en1,thetph,phiph,en1lab)
c
c  this rot gets us the scattered nucleon angle in z dir even if
c  nucleon
c cannot be emitted-so we can get the angle of the 1p1h remaining
c this returns new thet1,phi1 for n1 with energy en1,and angles for
c the
c particle-hole pair thetph,phiph from the incident nucleon with ein
c energy
c 2000 is to treat 1p1h left when this n1 cannot be emitted or
crescatters
c
       if(en1lab.le.0.)go to 2000
       ix=(en1lab+ed)/ed
       if(ix.lt.1)ix=1
       if(rndm(-1.).gt.q1(n1,ix))go to 2001
c
c here we see if n1 at en1 is emitted or rescatters(go to 2001
c rescat)
c if branch negative,particle rescatters so we go to 2000 for 1p1h
c calc
c set theta/5 deg bin index to store  scattered cross section
       ith=thet1*11.459155+1.
       if(ith.eq.37)ith=36
c add recoil/lab conv
       eloss=eloss+en1
       ja=ja+ida(n1)
       jz=jz+idz(n1)
       rmas=amass+2.-float(ja)-float(jz)
       acm=en1lab*rmas/(rmas+1.)
       pcm=sqrt(2.*acm)
       qcm=sqrt(2.*en1lab)
c  9/20/2014 change  ejectile channel energy to CM
         ix=(acm+ed)/ed
c------------------------------------------------------
       scrs(ix,n1)=scrs(ix,n1)+fcspec
       dsgp(ith,ix,n1)=dsgp(ith,ix,n1)+fcspec
c------------------------------------------------------
       if(iend.eq.1)call endf(thet1,en1lab,ith,ix,n1)
c add spin change algos july 2001
       tspin=spincon*qcm
       delj=tspin*cos(thet1)
c add perp 1/20/02
       adelj(ii)=delj
       spinf=spinf+delj
       nu=n1
       call recoil(6)
       go to 2000
c here thet1,phi1 are angles in z dir of nucleon before rescattering  in
c rescat
c this is nucleon with en1
c
c now store parameters of nucleon which rescatters rather than exit
 2001  if(en1lab.le.0.01)go to 2000
c only consider rescatter if nucleon exceeds binding energy by a
c constant-here it is 5 MeV
       jj=jj+1
       if(jj.gt.99)jj=99
       einc(jj)=en1
       thett(jj)=thet1
       phit(jj)=phi1
       idnuc(jj)=n1
 2000  continue
       thet1=thetph
       phi1=phiph
c here we do 1p1h ;we reset to the angles of the 1p1h in z-frame
c we have saved angle of 1p1h in z dir as thetph,phiph
c first,pick energy of particle n11
       u1=ein -en1
       if(u1.le.be(jz,ja,n11))go to 5000
       kz=jz+idz(n11)
       ka=ja+ida(n11)
       iares=mmt-ja-jz
       izres=mzee-jz
       if(kz.gt.15)kz=15
       if(ka.gt.24)ka=24
       x=rndm(-1.)
       if(u1.le.v)then
       en11=u1*x
       else
        en11=v*x+u1-v
        end if
c
       if(en11.gt.u1)en11=u1-0.01
       en11lab=en11-be(jz,ja,n11)
       ehol=u1-en11
         ken=(en11+ed)/ed
         if(ken.lt.1)ken=1
         cff(n11,ken)=cff(n11,ken)+fcs
c here we must compute angle of particle from 1p1h in z dir,then of the
c hole
c i.e. pick angle of 1p from 1p1h,leaving a 1h residue
c call rotang2 to get the 1p angles,then do hole value after return
       call rotang2(p3,u1,en11,thol,phol,en11lab)
       if(ehol.gt.beff(jz,ja,nhol))call holscat1(nhol,ehol,thol,
     1phol)
       if(light.gt.0)go to 5000
       if(en11lab.lt.0.01)go to 5000
c now is this 1p from the second 1p1h going to emit or scatter?
       ix=(en11lab+ed)/ed
       if(rndm(-1.).gt.q1(n11,ix))go to 5001
       ja=ja+ida(n11)
       jz=jz+idz(n11)
       ith=thet1*11.459155+1.
       if(ith.eq.37)ith=36
c add recoil/lab conv
       nu=n11
       rmas=amass+2.-float(ja)-float(jz)
       ecm=en11lab*rmas/(rmas+1.)
        ix=(ecm+ed)/ed
       pcm=sqrt(2.*ecm)
       qcm=sqrt(2.*en11lab)
c----------------------------------------------------
       scrs(ix,n11)=scrs(ix,n11)+fcspec
       dsgp(ith,ix,n11)=dsgp(ith,ix,n11)+fcspec
c-----------------------------------------------------
       eloss=eloss+en11
       if(iend.eq.1)call endf(thet1,en11lab,ith,ix,n11)
c end spin transfer this emission
c add spin change algos july 2001
       tspin=spincon*qcm
       delj=tspin*cos(thet1)
       adelj(ii)=delj
c      xq=cos(thet1)
c      if(nout.ne.1)go to 72001
c      if(igate.eq.0)go to 72001
c      jfin=(abs(spinu-delj)+1.)
c      do 72002 jgate=1,igate
c      if(n11.ne.npart(jgate))go to 72002
c      if(en11lab.lt.egatl(jgate).or.en11lab.gt.egatu(jgate))go to 72003
c      spind(jgate,jfin)=spind(jgate,jfin)+fcs
c      go to 72001
c72003  continue
c72002  continue
c72001  continue
       spinf=spinf+delj
c end spin transfer this emission
       nu=n11
       call recoil(7)
       go to 5000
c here,if it is to rescat,thet1,phi1 are its initial direction
c in the zyx frame
 5001  if(en11lab.le..01)go to 5000
       jj=jj+1
       einc(jj)=en11
       thett(jj)=thet1
       phit(jj)=phi1
       idnuc(jj)=n11
c
c now test for hot nucleons,and return to top if found
 5000  continue
       kp=kp+1
       if(ii.lt.jj.and.light.eq.0.and.rmas.gt.5.)go to 10200
 
       spinsv=spinf
       ispin=(abs(2.*(spinu-spinf))+1.)
       if(ispin.gt.180)ispin=180
       if(ispin.lt.1)ispin=1
       eftem=xmax-eloss
       if(arote(jz,ja,ispin).gt.eftem.and.mmt.gt.20)go to 10205
c      do 6001 ii=1,ihist
c      ix=ide(ii)
c      n11=idn(ii)
c      ith=idith(ii)
c      nout=iddum(ii)
c      enlab=enghis(ii)
c      spinf=spinfh(ii)
c      scrs(ix,n11)=scrs(ix,n11)+fcspec
c      dsgp(ith,ix,n11)=dsgp(ith,ix,n11)+fcspec
c       endif
c7001  continue
c      if(igate.eq.0)go to 92001
c      if(nout.ne.1)go to 92001
c      do 92002 jgate=1,igate
c      if(n11.ne.npart(jgate))go to 92002
c      if(enlab.lt.egatl(jgate).or.enlab.gt.egatu(jgate))go to 92003
c      jfin=(abs(spinu-adelj(ii))+1.)
c      spind(jgate,jfin)=spind(jgate,jfin)+fcs
c      go to 92001
c92003  continue
c92002  continue
c92001  continue
c end spin transfer this emission
c96002  continue
c      write(*,*)'after hist loop'
c if the above event is below yrast line after all pe emission,restart the event
 6000  continue
c  ---------------------------------------
c now do storage
c  ----------------------------------
       efinal=xmax-eloss
       ixm=(efinal+ed)/ed
          if(ixm.lt.1)ixm=1
c add call to evap subroutine for MC evaporation
c skip ppmc storage but increment pp array if eloss is large enough
c to prevent evaporation
       if(jz.gt.13)jz=13
       if(ja.gt.23)ja=23
       if(jz.le.iprot.and.ja.le.neut2)go to 10005
        if(light.gt.0)go to 10005
        ppp=ppp+fcs
       go to 10000
10005  continue
       if(light.eq.0)then
       if(efinal.lt.eminn(jz,ja)) then
c----------------------------------------------------------
         pp(jz,ja,ixm)=pp(jz,ja,ixm)+fcs
c----------------------------------------------------------
        spinxt=spinu-spinf
        spinr=sqrt(spinxt**2+spinprp**2)
         ispinr=abs(spinr+0.5)
          spin=ispinr
      ijz=zee+1.-float(jz)
      ija=atsave+2.-float(jz)-float(ja)
      if(nolevs(ijz,ija).eq.1)then
      x=rndm(-1.)
      if(x.lt.0.5)then
       spin=abs(spin-0.5)
       else
       spin=abs(spin+0.5)
                  endif
                              endif
          jl=(2.*spin)+1.
        if(jl.gt.180)jl=180
        if(jl.lt.1)jl=abs(jl)
        if(ixm.le.100)then
        popej(jz,ja,jl,ixm)=popej(jz,ja,jl,ixm)+fcs
        popj(jz,ja,jl)=popj(jz,ja,jl)+fcs
                       endif
                                                     endif
       xtran=xmax
       if(efinal.ge.eminn(jz,ja))call evap
                     endif
 
         if(iend.eq.0)go to 10002
       if(nout.gt.0.and.nout.le.3)call endvsort
c add recoil energy storage
10002  continue
       call reclsv
10000  continue
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c10001  continue
       if(itran.eq.0)itran=1
       csgamm(nepr)=pp(1,1,itran)
       if(no.eq.niso)then
      call sortmc
      call bngwrt
         if(iend.gt.0)call endvout
                     endif
      call precomgam
      pre=0.
      call gama
      pre=1.
1112   format(' monte carlo precompound with pairing correction')
      return
      end
c
      subroutine excl4
c modified starting dec7,2000 to use random coupling angles of HI excitons
c calls exciton sampling,clusterin s.r. for exciton energies
      common/kkg2/fractot(99),cs90
      common/clustr1/eahe4m,eadm,eatm,eahe3m
       common/spins4/eprot,ealpha,eoff(8)
      common/temscalf/b2h,b3h,b3he,b4he,b4hep
      common/temscale/a2h,a3h,a3he,a4he,a4hep
      common/clustr/m3set
      common/distb/dsgp(36,999,8)
      common/atsv/atsave
      common/shft/mc,mp,inver,ike,ipch,kq,qval,ap,at,zp,zt,cld,barfac
      common/astp3/niso,no,iadlt,isotag,massno(12),mass1(12),ncsave
      common/endvparm/iend,lll(6,6,6),ll(6,6)
      common/history/ide(60),idn(60),idith(60),iddum(60),the(60),adelj
     1(60),spinfh(60),enghis(60)
      common/out1/nout
      common/histindx/ihist
      common/giani/pre, sigt(1100),sigam(999),sigpre(999),gsp(999)
     1,eqgam(999)
      common/dec/reccon
      common/tr/xtran
      common/trans/itran
      common/himbc2/clthet(100),clphi(100)
      common/himbc1/threc,phrec
      common/himbc/thx,phx,przz,px,py,pz
      common/spins2/spinu,spinf,prz1
      common/rotor/irotor,iidummy,arote(15,25,181)
      common/ismrdat2/ididl(106,260),nolevs(106,260)
      common/gams/gamgam,csgam,csgamm(100)
      REAL*8 gamgam,csgamm,csgam
      common/isomer3/iziso(690),iaiso(690),csgs(690),csis(690)
     1,csgse(690,100),csise(690,100),csis2e(690,100),rm1(690,100),
     2rm2(690,100)
      common/lin/lspin
      common/lscpl/bx(180),sppj
      common/ijcpl/ajtarg,ajin
      common/tagspin/trgspin
      common/mmtcom/mmt,mzee,iares,izres,light
      common/par2/cncss
      common/outlim/limout,lim2
      common/lout/spincon
      common/iihist/kp
      common/levop/ldopt,igam
      common/perp/spinprp
      common/sq/sqx,efds(40,11),pfermx(40,11)
      common/spins/rote(15,24,181),popej(15,24,180,100),fcx(180),ai,ecmp
      common/sping/pfpej(600,180,100),pfpj(600,180)
      common/gate1/spind(5,180),popj(15,24,180),egatl(5),egatu(5),npart(
     15),igate
      common/arrays/eout(300),efd(40,11),pferm(40,11)
     1,tsetlb(181,300),edist(300),eric(300),ie(40),itl(40),
     2pfpe(40,11),efinc(40,11),eferm2(40,11),phirad(100),peincsq,peinc
       common/tang2/anexc(40),anexci(40),eferm(40,11),psq(40,11),
     1emit(300),angl(181),thetlb(37,300),pdist(300),pd(300)
      common/tang3/efermp,eincc,nexc
      common/nuopt/kinem,iem,kb
      common/excid/npartt(2)
      common/tem/philabb(36),phirec(36),phb(36),phin(36)
      common/bound/neut2,iprot,idnuc(100),idm(100),einc(100),
     1thett(100),phit(100),ppp
      common/hyb2/pp(15,24,999)
      common/nhy/ij,jl,ji,jq,bd,bx1,bx2,cmx,cv,bav,bost,b(3)
      common/evapmc/jcount,ipremc,eminn(15,24)
      common/ilcom/nu,ndum,am(8)
      common/rcl/amis
      common/labcon/prx,pry,prz,pcm,rmas,clab(999,8),dslb(36,999,8)
      common/exclbe/idz(4),ida(4),beff(20,20,2),pairmc(18,27)
      common/scat/jz,ja, probxx(2,2),q1(8,1100),fcs,eloss,fcspec
        common/mccom/ncount,ipot,rum(15,24),ppmc(15,24,999),cf
      common/exclus/ gib(8,999),pr(3),ge(2100),zz(3),
     1delta(3),ts(3),pairx(8),asum(10),rate(8,999)
      common/pl3/jx,nepr,eb(100),rcsp(100),zee,amass,plex
      common/par3/ep,sigml(999),acrs(999)
      common/memo/scrs(999,8)
      common/incr/ed,tened,efinal
      common/sft5/exc(18,27),xmxx
      common/lym1/be(18,27,8),pair(18,27),xmas(18,27),delshl(18,27),amas
     1(18,27)
c     common/bym1/symb(18,27),symbp(18,27),symbs(18,27)
      common/dist/n,ni,ke,irefr,xz(3),pod,pcs,ek,xmix,cfrac
     1,dtau,sign(36),tmix,bji(10),fin,dtup
      common/error/rerror(15,24)
      common/eder/ed2,v,vhalf,vsq
      common/angprm/aprp(500),eav3(1000),
     1eav2(1000),pot2x,thet1,phi1,snt(360),cst(360),snp(720),csp(720)
     2,pot
crc
      common/out2/vrl(200,25),vrlz(200,25),vrlx(200,25),vrly(200,25)
      common/out3/vscale
        common/out4/asm(36),asmm(36)
crc
      common/recoils/rcletza(100,13,24,18),rcleza(100,13,24),
     1rclea(100,36),rcleat(100,36,18),rphir(36)
     2,rclsp(200),rclz(200),rcly(200),rclx(200),escale
     3,rcle(200),rclt(36),rclet(200,36)
      common/rec1/ermxc,thetrc
      common/eferm/efermt
      common/delts/del(81),delt(81)
      common/exci/excdist(700),delexc,counter,ctloop,excdist1(700)
      common/clust1/denomc,connum,ptproj
      common/exci2/excdis2t(700)
      common/hidd/hiddcs
      common/fermst/aferm
c
c if ipot.gt.0 use infinite hole depth,elso finite as default
c next divide number of mc events by partial rx cross sect for gdh
c model
c or should use total for hybrid model
c recoil info added to this version 7/14/97;this will be known as
c 'bastille day alice.we define the following indices for various
c recoil arrays:
c energy index,1-100,ermax=1.25*e-lab
c projectile*mass(proj)/mass(cn)
c so index to store in array = erecoil*100./ermax+1.
c for theta index,use theta-recoil*18./pi+1. to give 10 degree
c intervals
c constants for particle id in exciton scatter routine
c above is temp array to check method oj channel j calc
        a2h=0.20
         a3h=0.05
          a3he=0.08
           a4he=0.030*(198./amass)**0.5
        if(efermt.eq.0.)efermt=30.
        if(efermp.eq.0.)efermp=30.

       mmt=amass+2.
       mzee=zp+zt+1.
       npartt(1)=2
       npartt(2)=1
       kp=2
       ajtarg=trgspin
       thetrc=18./3.1416
       epap=sqrt(ep*ap)
       ermax=4.00*ep/(at+ap)
       ermxc=100./ermax
       escale=amass*50./ep
       vscale=50.*amass/(sqrt(2.*ap*ep))
       amis=amass+2.
c the above are index constants for recoil arrays
       przz=sqrt(2.*ep*ap)
       prz1=przz
       prx=0.
       pry=0.
       prz=przz
c+++++++++++++++++++++++++++++++++++++
       reccon=0.
       counter=0.
       cntloop=0.
       spinprp=0.
       do 9002 i=1,690
       excdist(i)=0.
       excdist1(i)=0.
       excdis2t(i)=0.
       csgs(i)=0.
       csis(i)=0.
 9002  continue
       do 9 je=1,300
       eout(je)=0.
       edist(je)=0.
       eric(je)=0.
       do 19 jan=1,37
       thetlb(jan,je)=0.
   19  continue
    9  continue
      do 10103 i=2,100
        einc(i)=0.
        idnuc(i)=0
10103  continue
c+++++++++++++++++++++++++++++++++++++
crc
       write(17,1112)
c---------------------------------------------
          if(iend.gt.0.and.no.eq.1)call zero(9)
c-------------------------------------------
       call initial
       call trigstuf
c-------------------------------------------
c initialize arrays
      nexc=ap
      ef=efermp
       eincc=xmxx/ap
       sqx=sqrt(2.*(efermt+eincc))
       pfermpmx=sqrt(2.*efermp)
c calculate mid point efermi with fermi gas spacing-needs correction
        nproj=ap
       do710 l=1,nproj
        anexc(l)=(float(l)/ap)**(2./3.)*efermp
  710  continue
       del(1)=anexc(1)
       anexc(nproj+1)=anexc(nproj)+del(1)
       do 720 l=2,nproj
       del(l)=anexc(l)-anexc(l-1)
       delt(l)=del(l)/10.
       do 720 k=1,10
       efd(l,k)=anexc(l-1)+float(k-1)*delt(l)
       efds(l,k)=efd(l,k)+eincc
       pferm(l,k)=sqrt(2.*efd(l,k))
       pfermx(l,k)=pferm(l,k)*sqx
       eferm(l,k)=2.*efd(l,k)
  720  continue
       delt(1)=del(1)/10.
       do 721 k=1,10
       efd(1,k)=float(k)*delt(1)
       pferm(1,k)=sqrt(2.*efd(1,k))
       pfermx(1,k)=pferm(1,k)*sqx
       efds(1,k)=efd(1,k)+eincc
       eferm(1,k)=2.*efd(1,k)
  721  continue
c+++++++++++++++++++++++++++++++++++
c define a fermi energy of 30 MeV
       peincsq=2.*(eincc+efermt)
       peinc=sqrt(peincsq)
c april change to smooth as above
       etotin=eincc+efermt
       denomc=2.*sqrt(efermp*etotin)
       connum=efermp+eincc
       ptproj=sqrt(2.*efermp)
c begin changes for closed form solutions
c
       pot=30.
       pol=pot
       pot2x=pot*2.
       v=pot
       vhalf=v/2.
       vsq=v*v
cxxx
c  define z,a index offset increments for n(=1) and p(n=2)
c  emission
       idz(1)=0
       idz(2)=1
       ida(1)=1
       ida(2)=0
c   6  continue
       ppp=0.
       bd=ap
       jc=bd
       jj=1
       eincc=xmxx/bd
c----------------
         cs90=0.
         do l=1,90
         cs90=cs90+sigml(l)
         enddo
        fcs=cs90/float(ncsave)
c       delexc=1./float(ncount)
        delexc=1./float(ncsave)
        fcspec=fcs/ed
c----------------
       call spinch
c-----------------------------------------------------
       do 10000 icnt=1,ncount
c---------------------------------------------------------
        call xjlcpl(jzl)
        call ijlcpl
       
c--------------------------------------------------------

       jz=1
       ja=1
       reccon=0.
c jzl is a.m. of composite before decay
       ajcpl=2.*spinu+1.
         j=jzl+1
         if(j.lt.1)j=1
        xmax=xmxx-arote(1,1,j)
       eincc=xmax/bd
c
       light=0
       ihist=0
c
       nnn=1
c10217  continue
c set switch 'light' to tell which light residues will undergo breakup
c with light=1 meaning breakup algo has been applied,so no residue to store
c now initialize arrays on parameters for this event
c///////////////////////////////////////////////////////////////
       do 10100 i=1,100
       einc(i)=0.
       thett(i)=0.
       phit(i)=0.
       idnuc(i)=0.
10100  continue
c below we zero history arrays for event- emitted nucleon pptys
       do 10202 inum=1,60
       ide(inum)=0
       idn(inum)=0
       idith(inum)=0
       the(inum)=0.
       adelj(inum)=0.
       spinfh(inum)=0.
       enghis(inum)=0.
10202  continue
cxx
c/////////////////////////////////////////////////////////////////
crc
       jj=1
       ii=1
       spinprp=0.
       spinf=0.
       delj=0.
       nout=0
c dec 98 for clusters
       prz=przz
       prx=0.
       pry=0.
       eloss=0.
       jz=1
       ja=1
c here define a cluster momentum wrt bottom of composite pot
       pmbcz=(bd)*sqrt(2.*(xmax/bd+pot))
       jbd=ap
c set recoil theta,phi to initial zero-beam direction-values
       threc=0.
       phrec=0.
       do 91735 ibd=1,jbd
       clthet(ibd)=0.
       clphi(ibd)=0.
91735  continue
       tph=0.
       pph=0.
         phi1=0.
c above define initial recoil directions
c define arrays of theta and phi for cluster nucleons from mbc before rotation
           if(jc.eq.1)go to 91732
c 2/21/01 replace with call to clusterin      call hiin
          xtran=xmax
        if(hiddcs.eq.1)then
c------------------------------------------------------------
       call clustin1(tph,pph)
       else
       call hiin
c----------------------------------------------------------
       endif
c this should return,for one event,the initial exciton energies(einc()),and angles
c for hiin, and the energies only for clustin1
c in radians(angl()-theta,plus exciton id idnuc())
       ihist=0
c-----------------------------------------------++++++++++++++++++++
       do 91731 jexc=1,jc
       n1=idnuc(jexc)
c--------------------------
       ecall=einc(jexc)
       if(n1.eq.3)go to 91733
       if(einc(jexc).lt.0.1)go to 91731
c----------------------------------------
       ipart=idnuc(ii)
       phi1=phit(ii)
       thet1=thett(ii)
        if(m3set.gt.0) call clustpe(ecall,en1,ipart,epap)
c----------------------------------------
       en1lab=einc(jexc)-be(jz,ja,n1 )
       if(en1lab.lt.0.)go to 91731
       ix=(en1lab+ed)/ed
       
       if(rndm(-1.).gt.q1(n1,ix))go to 91731
        phi1=phit(jexc)
        thet1=thett(jexc)
       ith=thet1*11.459155+1.
       if(ith.eq.37)ith=36
c add recoil/lab conv
       ja=ja+ida(n1)
       jz=jz+idz(n1)
       rmas=amass+2.-float(ja)-float(jz)
       acm=en1lab*rmas/(rmas+1.)
       pcm=sqrt(2.*en1lab)
c 9/20/2014 change energy index ix from channel to cm
        ix=(acm+ed)/ed
c----------------------------------------------------
         scrs(ix,n1)=scrs(ix,n1)+fcspec
         dsgp(ith,ix,n1)=dsgp(ith,ix,n1)+fcspec
c----------------------------------------------------
c may comment stores above to put back on history loop
c note that nout is also incremented in 'recoil'
c      iddum(ihist)=nout
c add spin change algos july 2001
       delj=(float((jzl-1)))*(pcm/prz1)*cos(thet1)
c add perp 1/20/02
       delprp=delj*sin(thet1)*sin(phi1)/cos(thet1)
       spinprp=spinprp+delprp
       spinf=spinf+delj
c note that lab sdcs,ddcs are stored/calc'd in recoil, could instead
c be put in history file herein to be used to check gated values in lab
       if(iend.eq.1)call endf(thet1,en1lab,ith,ix,n1)
       
       nu=n1
       call recoil(2)
cxxxx
       eloss=eloss+einc(jexc)
       einc(jexc)=0.
       go to 91731
c check if nucleons from cluster are emitted prior to rescattering
91733  continue
       nhol=1
       ehol=einc(jexc)
       if(rndm(-1.).gt.0.5)nhol=2
       call holscat(nhol,ehol,0.,0.)
       einc(jexc)=0.
c for now we don't keep track of nucleon emission from hole scattering
91731   continue
c----------------------------------------------------++++++++++++
c here is end of cascade of primary HI nucleons; next follow
c re-scattered HI nucleons.
91732   continue
       jj=jc+2
c
10200  continue
c---------return point hot nucleons----------------------
       ii=ii+1
       ein=einc(ii)
       if(ein.le.0.1)go to 5000
       thet1=thett(ii)
       phi1=phit(ii)
       ipart=idnuc(ii)
c----------------------------------
      if(m3set.gt.0)call clustpe(ein,en1,ipart,epap)
c----------------------------------------------
c nucleon of energy 'ein' above Fermi energy makes a 2p1h state with
c exciton 'n1' of energy 'en1' and 'n11' of energy 'en11',with a hole
c nhol of energy ehol=ein-en1-en11.If these excitons are energetic
c enough,they are stored to rescatter in this 'loop' if they are not
c emitted.
       call partid(ein,ipart,n1,n11,nhol)
c  partid returns the particle1,2,and hole identities for the
c incident
c  nucleon at energy ein( as n1,n11,nhol)
c  begin additions for daughter pairing shift,and calc en1
c
       kz=jz+idz(n1)
       ka=ja+ida(n1)
       iares=mmt-ja-jz
       izres=mzee-jz
       if(kz.gt.13)kz=13
       if(ka.gt.22)ka=22
       icycle=0
        eff=ein -pairmc(kz,ka)
 635   x=rndm(-1.)
       if(eff.le.v)then
        en1=eff*(1.-sqrt(1.-x))
       go to 40
       end if
       br=(eff-v)/(eff-vhalf)
       if(x.gt.br)then
       en1=eff-sqrt(v*(2.*eff-v)*(1.-x))
       else
       en1=x*(eff-vhalf)
       end if
c
c now en1 is the energy of the 1p of the 2p1h config following rescatter
c of
c the first 1p from a 2p1h config with energy ein  for the primary 1p
c
   40  continue
c
c 740  continue
       icycle=icycle+1
       if(en1.gt.ein.and.icycle.eq.1)go to 635
       if(en1.gt.ein.and.icycle.eq.2)en1=ein-0.01
       en1lab=en1-be(jz,ja,n1 )
       e1p1h=ein -en1
c
        call rotangl(p3,ein ,en1,thetph,phiph,en1lab)
c
c  this rot gets us the scattered nucleon angle in z dir even if
c  nucleon
c cannot be emitted-so we can get the angle of the 1p1h remaining
c this returns new thet1,phi1 for n1 with energy en1,and angles for
c the
c particle-hole pair thetph,phiph from the incident nucleon with ein
c energy
c 2000 is to treat 1p1h left when this n1 cannot be emitted or
crescatters
c
       if(en1lab.le.0.)go to 2000
       ix=(en1lab+ed)/ed
       if(rndm(-1.).gt.q1(n1,ix))go to 2001
c       emix=emix+1.
c
c here we see if n1 at en1 is emitted or rescatters(go to 2001
c rescat)
c if branch negative,particle rescatters so we go to 2000 for 1p1h
c calc
c set theta/5 deg bin index to store  scattered cross section
       ith=thet1*11.459155+1.
       if(ith.eq.37)ith=36
c add recoil/lab conv
       eloss=eloss+en1
       ja=ja+ida(n1)
       jz=jz+idz(n1)
       rmas=amass+2.-float(ja)-float(jz)
       acm=en1lab*rmas/(rmas+1.)
       pcm=sqrt(2.*en1lab)
       qcm=sqrt(2.*en1lab)
c   9/20/2014 change energy index ix from channel to cm
         ix=(acm+ed)/ed
c-------------------------------------------------
       scrs(ix,n1)=scrs(ix,n1)+fcspec
       dsgp(ith,ix,n1)=dsgp(ith,ix,n1)+fcspec
c------------------------------------------------
       if(iend.eq.1)call endf(thet1,en1lab,ith,ix,n1)
c add spin change algos july 2001
       delj=(float((jzl-1)))*(pcm/prz1)*cos(thet1)
c add perp 1/20/02
       delprp=delj*sin(thet1)*sin(phi1)/cos(thet1)
       spinprp=spinprp+delprp
       adelj(ii)=delj
       spinf=spinf+delj
       nu=n1
cxxxx
c note that lab sdcs,ddcs are stored/calc'd in recoil, could instead
c be put in history file herein to be used to check gated values in lab
       call recoil(3)
       go to 2000
c here thet1,phi1 are angles in z dir of nucleon before rescattering  in
c rescat
c this is nucleon with en1
c
c now store parameters of nucleon which rescatters rather than exit
 2001  if(en1lab.le.0.)go to 2000
c only consider rescatter if nucleon exceeds binding energy by a
c constant-here it is 5 MeV
       jj=jj+1
       if(jj.gt.99)jj=99
       einc(jj)=en1
       thett(jj)=thet1
       phit(jj)=phi1
       idnuc(jj)=n1
 2000  continue
       thet1=thetph
       phi1=phiph
c here we do 1p1h ;we reset to the angles of the 1p1h in z-frame
c we have saved angle of 1p1h in z dir as thetph,phiph
c first,pick energy of particle n11
       u1=ein -en1
       if(u1.le.be(jz,ja,n11))go to 5000
       kz=jz+idz(n11)
       ka=ja+ida(n11)
       iares=mmt-ja-jz
       izres=mzee-jz
       if(kz.gt.13)kz=13
       if(ka.gt.22)ka=22
       x=rndm(-1.)
       if(u1.le.v)then
       en11=u1*x
       else
        en11=v*x+u1-v
        end if
c
       if(en11.gt.u1)en11=u1-0.01
       en11lab=en11-be(jz,ja,n11)
       ehol=u1-en11
c here we must compute angle of particle from 1p1h in z dir,then of the
c hole
c i.e. pick angle of 1p from 1p1h,leaving a 1h residue
c call rotang2 to get the 1p angles,then do hole value after return
       call rotang2(p3,u1,en11,thol,phol,en11lab)
       if(ehol.gt.beff(jz,ja,nhol))call holscat1(nhol,ehol,thol,
     1phol)
       if(light.gt.0)go to 5000
       if(en11lab.lt.0.01)go to 5000
c now is this 1p from the second 1p1h going to emit or scatter?
       ix=(en11lab+ed)/ed
       if(rndm(-1.).gt.q1(n11,ix))go to 5001
       ja=ja+ida(n11)
       jz=jz+idz(n11)
       ith=thet1*11.459155+1.
       if(ith.ge.37)ith=36
c change channel energy to cm  for ix 9/20/2014
        rmas=amass+2.-float(ja)-float(jz)
         acm=en11lab*rmas/(rmas+1.)
           ix=(acm+ed)/ed
c----------------------------------------------------------
       scrs(ix,n11)=scrs(ix,n11)+fcspec
       dsgp(ith,ix,n11)=dsgp(ith,ix,n11)+fcspec
c-----------------------------------------------------------
       eloss=eloss+en11
       pcm=sqrt(2.*en11lab)
       delj=(float((jzl-1)))*(pcm/prz1)*cos(thet1)
       delprp=delj*sin(thet1)*sin(phi1)/cos(thet1)
       spinprp=spinprp+delprp
       if(iend.eq.1)call endf(thet1,en11lab,ith,ix,n11)
c add recoil/lab conv
       nu=n11
       pcm=sqrt(2.*en11lab)
       qcm=sqrt(2.*en11lab)
c add spin change algos july 2001
       delj=(float((jzl-1)))*(pcm/prz1)*cos(thet1)
       delprp=delj*sin(thet1)*sin(phi1)/cos(thet1)
       spinprp=spinprp+delprp
c72003  continue
c72002  continue
c72001  continue
       spinf=spinf+delj
c end spin transfer this emission
       nu=n11
c note that lab sdcs,ddcs are stored/calc'd in recoil, could instead
c be put in history file herein to be used to check gated values in lab
        call recoil(4)
       go to 5000
c here,if it is to rescat,thet1,phi1 are its initial direction
c in the zyx frame
 5001  if(en11lab.le..01)go to 5000
       jj=jj+1
       if(jj.gt.99)jj=99
       einc(jj)=en11
       thett(jj)=thet1
       phit(jj)=phi1
       idnuc(jj)=n11
 5000  continue
c now test for hot nucleons,and return to top if found
c      if(ii.lt.jj.and.light.eq.0)go to 10200
       if(ii.lt.jj.and.rmas.gt.5.)go to 10200
 
       spinsv=spinf
       ispin=(abs(2.*(spinu-spinf))+1.)
       eftem=xmax-eloss
c      if(arote(jz,ja,ispin).gt.eftem.and.mmt.gt.20)go to 10205
c      do 6001 ii=1,ihist
c again we note that lab spectra have been stored in sub recoil,
c and should instead be put into history arrays to gate on
c      ix=ide(ii)
c      n11=idn(ii)
c      ith=idith(ii)
c      nout=iddum(ii)
c      enlab=enghis(ii)
c      spinf=spinfh(ii)
c      scrs(ix,n11)=scrs(ix,n11)+fcspec
c      dsgp(ith,ix,n11)=dsgp(ith,ix,n11)+fcspec
       if(igate.eq.0)go to 92001
c      if(nout.ne.1)go to 92001
c      do 92002 jgate=1,igate
c      if(n11.ne.npart(jgate))go to 92002
c      if(enlab.lt.egatl(jgate).or.enlab.gt.egatu(jgate))go to 92003
c      jfin=(abs(spinu-adelj(ii))+1.)
c      spind(jgate,jfin)=spind(jgate,jfin)+fcs
c      go to 92001
c92003  continue
c92002  continue
92001  continue
c end spin transfer this emission
c      if(ii.lt.jj)go to 10200
c6001  continue 
c96002  continue
c6000  continue
c now do storage
c  ---------------------------------------
c  ----------------------------------
       efinal=xmax-eloss
       if(efinal.lt.0.)efinal=0.
       ixm=(efinal+ed)/ed
c      spinf=spinsv
c MAR 970414
c add call to evap subroutine for MC evaporation
c skip ppmc storage but increment pp array if eloss is large enough
c to prevent evaporation
cbound
       if(jz.gt.13)jz=13
       if(jz.le.iprot.and.ja.le.neut2)go to 10005
        if(light.gt.0)go to 10005
        ppp=ppp+fcs
       go to 10000
10005  continue
       if(light.eq.0)                               then
       if(efinal.lt.eminn(jz,ja))          then
c--------------------------------------------------------
         pp(jz,ja,ixm)=pp(jz,ja,ixm)+fcs
c------------------------------------------------------
        spinxt=spinu-spinf
        spinr=sqrt(spinxt**2+spinprp**2)
         ispinr=abs(spinr+0.5)
          spin=ispinr
        if(spin.gt.90)spin=90
      
      ijz=zee+1.-float(jz)
      ija=amass+2.-float(jz)-float(ja)
      if(nolevs(ijz,ija).eq.1)then
      x=rndm(-1.)
      if(x.lt.0.5)then
       spin=abs(spin-0.5)
       else
       spin=abs(spin+0.5)
                  endif
                              endif
          jl=(2.*spin)+1.
c 01/25/02
        if(jl.gt.180)jl=180
        if(jl.lt.1)jl=abs(jl)
        if(ixm.le.100)then
        if(jl.eq.0)jl=1
        popej(jz,ja,jl,ixm)=popej(jz,ja,jl,ixm)+fcs
        popj(jz,ja,jl)=popj(jz,ja,jl)+fcs
                       endif
                                           endif
       xtran=xmax
       if(efinal.ge.eminn(jz,ja))call evap
       reccon=1.
                                                    endif
         if(iend.eq.0)go to 10002
       if(nout.gt.0.and.nout.le.3)call endvsort
 
c add recoil energy storage
10002  continue
       call recoil(5)
       call reclsv
10000  continue
      write(50,*)' ein    excdist(inx) hiin-coupling'
      write(51,*)' ein    excdist1(inx) clusterin1'
      write(52,*)' ein    exciton(inx) exciton '
      sumexd=0.
      sumexd1=0.
      sumexciton=0.
      constant=(xmax/ap)**(ap-1.)
      do 10009 inx=1,690
      exciton=0.
      u=xmax-ed*float(inx)
      if(u.lt.0.)go to 10008
        u=u/ap
      exciton=((u )**(ap-2.))*ed/constant
      exciton=exciton*(ap-1.)
      sumexciton=sumexciton+exciton
       ratiox=excdist(inx)/exciton
10008 continue
      if(excdist(inx).eq.0..and.excdist1(inx).eq.0.)go to 10009
      sumexd=sumexd+excdist(inx)
      sumexd1=sumexd1+excdist1(inx)
      ein=float(inx)*ed-ed/2.
      exdser1=sqrt(excdist1(inx)*delexc)
      exdser=sqrt(excdist(inx)*delexc)
      raterr=(exdser/excdist(inx))*ratiox
      write(50,*)ein,excdist(inx),exdser
      write(51,*)ein,excdist1(inx),exdser1
      write(52,*)ein,exciton
      write(53,*)ein,ratiox,raterr
10009 continue
      write(50,*)' sum hiin    ',sumexd
      write(51,*)' sum clust  ',sumexd1
      write(52,*)' sum exciton   ',sumexciton
         if(no.eq.niso)then
      call sortmc
      call bngwrt
        if(iend.gt.0)call endvout
                        endif
1112   format(' monte carlo precompound with pairing correction')
      return
      end
c-----------------------------------------------------------
       subroutine lscupl
      common/outlim/limout,lim2
      common/lscpl/bx(180),sppj
c we will assume proj. spins gt.0.5 couple as +/- only as for nucleons
      bx(1)=1.
      do 10  i=2,40
      al=float(i)-sppj
      sum=(2.*al+1.)+2.*(al+1.)+1.
      bx(i)=(2.*al+1.)/sum
  10  continue
      return
      end
       subroutine lscoupl(delj)
      common/outlim/limout,lim2
      common/lscpl/bx(180),sppj
c to couple outgoing wave with spin half
      itest=0
      idelj=(abs(delj)+sppj)
c  integerize delj
      if(idelj.eq.0)then
      if(delj.ge.0.)delj=sppj
      if(delj.lt.0.)delj=-sppj
      itest=1
      endif
      if(itest.eq.1)return
      if(delj.gt.0.)then
      delj=float(idelj)+sppj
      if(rndm(-1.).lt.bx(idelj))delj=delj-sppj
       endif
      if(delj.lt.0.)then
      delj=-float(idelj)-sppj
      if(rndm(-1.).lt.bx(idelj))delj=delj+sppj
      endif
      return
      end
c
      subroutine clustpe(ein,en1,ipart,epap)
      common/lscpl/bx(180),sppj
      common/clkin/elab
      common/endvparm/iend,lll(6,6,6),ll(6,6)
      common/astp3/niso,no,iadlt,isotag,massno(12),mass1(12),ncsave
      common/out1/nout
      common/histindx/ihist
      common/history/ide(60),idn(60),idith(60),iddum(60),the(60),adelj
     1(60),spinfh(60),enghis(60)
      common/bound/neut2,iprot,idnuc(100),idm(100),einc(100),
     1thett(100),phit(100),ppp
      common/nuc/cnu
      common/shft/mc,mp,inver,ike,ipch,kq,qval,ap,at,zp,zt,cld,barfac
      common/par3/ep,sigml(999),acrs(999)
      common/spins2/spinu,spinf,prz1
      common/lin/lspin
      common/temscalf/b2h,b3h,b3he,b4he,b4hep
      common/temscale/a2h,a3h,a3he,a4he,a4hep
      common/eder/ed2,v,vhalf,vsq
      common/hybsub/iav,idum,rs(8,1100),c,pox,r1
      common/pl3/jx,nepr,eb(100),rcsp(100),zee,amass,plex
      common/exclbe/idz(4),ida(4),beff(20,20,2),pairmc(18,27)
       common/spins4/eprot,ealpha,eoff(8)
      common/scat/jz,ja, probxx(2,2),q1(8,1100),fcs,eloss,fcspec
      common/exclus/ gib(8,999),pr(3),ge(2100),zz(3),
     1delta(3),ts(3),pairx(8),asum(10),rate(8,999)
      common/memo/scrs(999,8)
      common/incr/ed,tened,efinal
      common/sft5/exc(18,27),xmxx
      common/lym1/be(18,27,8),pair(18,27),xmas(18,27),delshl(18,27),amas
     1(18,27)
c     common/bym1/symb(18,27),symbp(18,27),symbs(18,27)
      common/clustr/m3set
      common/tempo/sp3of(999),sp3(999)
      common/clustr1/eahe4m,eadm,eatm,eahe3m
      common/distb/dsgp(36,999,8)
      common/ilcom/nu,ndum,am(8)
      common/debug/alphan
      common/labcon/prx,pry,prz,pcm,rmas,clab(999,8),dslb(36,999,8)
      common/angprm/aprp(500),eav3(1000),
     1eav2(1000),pot2x,thet1,phi,snt(360),cst(360),snp(720),csp(720)
     2,pot
c in this version we assume that the alpha pick-up mechanism
c leaves holes which will equilibrate; for knock-out, we return
c the particle-hole pair to continue the precompound process.
       alphan=0.
      jzl=lspin-1
c
      eahe4=ein-be(jz,ja,3)-eoff(3)
      ead=ein-be(jz,ja,4)-eoff(4)
      eat=ein-be(jz,ja,5)-eoff(5)
      eahe3=ein-be(jz,ja,6)-eoff(6)
      eahe4t=eahe4+eoff(3)
      eadt=ead+eoff(4)
      eatt=eat+eoff(5)
      eahe3t=eahe3+eoff(6)
c     
       if(eahe4.lt.0.)eahe4=0.
       if(ead.lt.0.)ead=0.
       if(eat.lt.0.)eat=0.
       if(eahe3.lt.0.)eahe3=0.
      if(eahe4.le.0..and.ead.le.0.)return
c calculate total emission widths for 4He(walpha), 2H(wdeut),3H(wtrit)
c make initial alpha excitation inversely energy dependent
      b4he=a4he*(200./ein)
c     if(ein.lt.40.)b4he=a4he*5.
      if(ein.lt.30.)b4he=a4he*6.666
      walpha=b4he*(eahe4/ein)**4
c make unphysical, data fitting second softer alpha emission with walphap
      walphap=2.*b4he*(eahe4/ein)**8
      wdeut=a2h*(ead/ein)
      wtrit=a3h*(eat/ein)**2
      whe3=a3he*(eahe3/ein)**2
      vv=20.
      vvhalf=10.
c above hole depth limits for high energy alpha component
      br=rndm(-1.)*(1.+walpha+wdeut+wtrit+whe3+walphap)
      if(br.le.1.)return
         ialph=0
c here we will eventually add widths for other clusters to 1.+...,
c and if br.lt.1.  'return'; it is a nucleon-nucleon scatter event
      if((br-1.).le.walpha)go to 3332
c it is an alpha excitation 
      bp=br-1.-walpha
      if(bp.le.wdeut)go to 3330
c it is a  deuteron excitation
      bp=bp-wdeut
      if(bp.le.wtrit)go to 3333
c it is a triton excitation
      bp=bp-wtrit
      if(bp.le.whe3)go to 3334
       bp=bp-whe3
c it is a  3He excitation
         ialph=0
      if(bp.le.walphap)ialph=1
      if(ialph.eq.1)go to 3332
       if(bp.gt.walphap)return
c it is harder component of alpha spectrum
        return
 3330  continue
c deuterons
c      w=40.
c       w2=w/2.
c      icycle=0
      nu=4
      kz=jz+1
      ka=ja+1
      tmas=amass+2.-float(kz)-float(ka)
c mb 4-20-2011
      if(tmas.lt.0.1)return
       if(kz.gt.13)kz=13
       if(ka.gt.22)ka=22
        eff=ein -pairmc(kz,ka)
       x=rndm(-1.)
        en1=x*eff
       if(en1.gt.ein)en1=ein-0.01
       en1lab=en1-be(jz,ja,4 )
       elab=en1lab
         if(elab.lt.0.1)return
         if(en1lab.le.0.)return
       pcm=sqrt(2.*en1lab*2.)
       iex=(en1lab+ed)/ed
        iey=(en1+ed)/ed
       if(iex.lt.1)return
       if(iey.lt.1)return
c iey 10/27/2010
       q=rate(4,iex)/(rate(4,iex)+rs(4,iey))
       if(rndm(-1.).gt.q)return
         thetcl=(1.-en1lab/eadt)*3.14159/2.
          if(thetcl.lt.0.)thetcl=0.
          if(thetcl.gt.1.57)thetcl=1.57
         thet1=thetcl
       ith=thet1*11.459155+1.
c 9/20/2014 for deuts,change echannel to ecm:
       ecmd=en1lab*tmas/(tmas+2.)
        iex=(ecmd+ed)/ed
c------------------------------------------------------
       dsgp(ith,iex,4)=dsgp(ith,iex,4)+fcspec
       scrs(iex,4)=scrs(iex,4)+fcspec
c------------------------------------------------------
       call recoil(8)
c
c at this point one should convert Ecm of the  d cluster to lab energy and store
c  this done in recoil sr
       if(iend.eq.1)call endf(thet1,en1lab,ith,iex,4 )
       eloss=eloss+en1
        ein=ein-en1
        ja=ja+1
        jz=jz+1
c add angular distribution algorithm, decrement spin for cluster emission
c        thetcl=(1.-en1lab/eadm)*3.14159/2.
         spinu=spinu-cos(thetcl)*(sqrt(2.*en1lab*2.)/epap)*float(jzl)
         return
 3332 continue
c transfer for 4He emission
      tmas=amass+2.-4.-float(ja)-float(jz)
      if(tmas.lt.0.1)return
c mb 4-20-2011
      nu=3
      x=rndm(-1.)
       ealp=ein-eahe4*(x **0.2500)
       if(ialph.eq.1) ealp=ein-eahe4*(x **(1./8.))
c 40   continue
       ealab=ealp-be(jz,ja,3)
        elab=ealab
         if(elab.lt.0.1)return
        pcm=sqrt(2.*4.*ealab)
       iex=(ealab+ed)/ed
       iey=(ealp+ed)/ed
      if(iex.lt.1)return
c add 9/14/2010
      if(iey.lt.1)return
       q=rate(3,iex)/(rate(3,iex)+rs(3,iey))
        x=rndm(-1.)
        if(x.le.q)go to 401
c here if alpha is not emitted, let it rescatter
       x=rndm(-1.)
       eff=ealp
       vl=40.
       if(eff.le.vl)then
        en1=eff*(1.-sqrt(1.-x))
       go to 400
       end if
       br=(eff-vl)/(eff-20.)
       if(x.gt.br)then
       en1=eff-sqrt(vl*(2.*eff-vl)*(1.-x))
       else
       en1=x*(eff-20.)
       end if
  400  continue
       ealp=en1
       ealab=ealp-be(jz,ja,3)
        elab=ealab
         if(elab.lt.0.01)return
        pcm=sqrt(2.*4.*elab)
       iex=(ealab+ed)/ed
       iey=(ealp+ed)/ed
      if(iex.lt.1)return
c add 9/14/2010
      if(iey.lt.1)return
       q=rate(3,iex)/(rate(3,iex)+rs(3,iey))
        x=rndm(-1.)
         if(x.gt.q)return
 401      kz=jz+2
      ka=ja+2
       if(kz.gt.13)kz=13
       if(ka.gt.22)ka=22
         thetcl=(1.-ealab/eahe4t)*3.14159/2.
          if(thetcl.lt.0.)thetcl=0.
          if(thetcl.gt.1.57)thetcl=1.57
         thet1=thetcl
       ith=thet1*11.459155+1.
       if(ith.lt.1)ith=1
c  9/20/2014 change channel energy index to cm index:
         ecma=ealab*tmas/(tmas+4.)
           iex=(ecma+ed)/ed
c-----------------------------------------------------
       scrs(iex,3)=scrs(iex,3)+fcspec
       dsgp(ith,iex,3)=dsgp(ith,iex,3)+fcspec
c-----------------------------------------------------
       alphan=1.
c
       call recoil(9)
c

c add lab conversion/storage alphas 9 nov 2009
       if(iend.eq.1)call endf(thet1,ealab,ith,iex,3)
       eloss=eloss+ealp
        ein=ein-ealp
        ja=ja+2
        jz=jz+2
c        thetcl=(1.-ealab/eahe4m)*3.14159/2.
         spinu=spinu-cos(thetcl)*(sqrt(2.*ealab*4.)/epap)*float(jzl)

        return
 3333  continue
c transfer for triton emission
      nu=5
      kz=jz+1
      ka=ja+2
      

      tmas=amass+2.-float(kz)-float(ka)
       if(tmas.lt.1.)return
c  mb 4/20/2011
       if(kz.gt.13)kz=13
       if(ka.gt.22)ka=22
      eff=ein-pairmc(kz,ka)
      x=rndm(-1.)
      et=eff-eat*(x **0.5)
        elab=et-be(jz,ja,5)
         if(elab.le.0.0)return
        pcm=sqrt(2.*3.*elab)
       iex=(elab+ed)/ed
       iey=(et+ed)/ed
c  9/20/2014 change iex index from channel to cm energy
        tiex=elab*tmas/(tmas+3.)
          iex=(tiex+ed)/ed
      if(iex.lt.1)return
c add 9/14/2010
      if(iey.lt.1)return
       q=rate(5,iex)/(rate(5,iex)+rs(5,iey))
       if(rndm(-1.).gt.q)return
         thetcl=(1.-et/eatt )*3.14159/2.
           if(thetcl.lt.0.)thetcl=0.
           if(thetcl.gt.1.57)thetcl=1.57
         thet1=thetcl
       ith=thet1*11.459155+1.
c------------------------------------------------------
       scrs(iex,5)=scrs(iex,5)+fcspec
       dsgp(ith,iex,5)=dsgp(ith,iex,5)+fcspec
c-----------------------------------------------------
c
        call recoil(10)
c
c add lab conversion tritons
c        enghis(ihist)=ealab
       if(iend.eq.1)call endf(thet1,elab,ith,iex,5 )
       eloss=eloss+et
        ein=ein-et
        ja=ja+2
        jz=jz+1
      tmas=amass+2.-float(jz)-float(ja)
       if(tmas.lt.1.)return
c  mb 4/20/2011
c        thetcl=(1.-elab/eatm)*3.14159/2.
         spinu=spinu-cos(thetcl)*(sqrt(2.*elab*3.)/epap)*float(jzl)
        return
 3334 continue
c transfer for 3He emission
      nu=6
      kz=jz+2
      ka=ja+1
      tmas=amass+2.-float(kz)-float(ka)
       if(tmas.lt.1.)return
c  mb 4/20/2011
       if(kz.gt.13)kz=13
       if(ka.gt.22)ka=22
        eff=ein-pairmc(kz,ka)
      x=rndm(-1.)
      ehe3=eff-eahe3*(x **0.5)
       elab=ehe3-be(jz,ja,6)
         if(elab.le.0.0)return
        pcm=sqrt(2.*3.*elab)
       iex=(elab+ed)/ed
       iey=(ehe3+ed)/ed
      if(iex.lt.1)return
c add 9/14/2010
      if(iey.lt.1)return
       q=rate(6,iex)/(rate(6,iex)+rs(6,iey))
       if(rndm(-1.).gt.q)return
         thetcl=(1.-elab/eahe3t)*3.14159/2.
          if(thetcl.lt.0.)thetcl=0.
          if(thetcl.gt.1.57)thetcl=1.57
         thet1=thetcl
       ith=thet1*11.459155+1.
c  9/20/2014 change channel energy to cm
          tecm=elab*(tmas)/(tmas+3.)
            iex=(tecm+ed)/ed
c---------------------------------------------------
       scrs(iex,6)=scrs(iex,6)+fcspec
       dsgp(ith,iex,6)=dsgp(ith,iex,6)+fcspec
c----------------------------------------------------
c  kkg  04/13/11
       call recoil(11)
c
c add lab conversion 3He emission
       if(iend.eq.1)call endf(thet1,elab,ith,iex,5 )
       eloss=eloss+ehe3
        ein=ein-ehe3
        ja=ja+1
        jz=jz+2
c        thetcl=(1.-elab/eahe3m)*3.14159/2.
         spinu=spinu-cos(thetcl)*(sqrt(2.*elab*3.)/epap)*float(jzl)
      return
      end
c---------------------------------------------
c sub to initialize arrays on event
      subroutine exclzro
      common/history/ide(60),idn(60),idith(60),iddum(60),the(60),adelj
     1(60),spinfh(60),enghis(60)
      common/histindx/ihist
      common/mmtcom/mmt,mzee,iares,izres,light
      common/spins2/spinu,spinf,prz1
      common/deljay/delj
      common/out1/nout
      common/bound/neut2,iprot,idnuc(100),idm(100),einc(100),
     1thett(100),phit(100),ppp
       do 10202 inum=1,60
       ide(inum)=0
       idn(inum)=0
       idith(inum)=0
       the(inum)=0.
       adelj(inum)=0.
       spinfh(inum)=0.
       enghis(inum)=0.
10202  continue
cxx
       ihist=0
       light=0
       spinf=0.
       delj=0.
       nout=0
       do 10100 i=1,99
       einc(i)=0.
       thett(i)=0.
       phit(i)=0.
       idnuc(i)=0.
10100  continue
      return
      end
c-------------------------------------------------------------------
c---------------add sr xjlcpl torandom select a.m.-------7/25/2014
      subroutine xjlcpl(jzl)
      common/spins2/spinu,spinf,prz1
      common/tagspin/trgspin
      common/lscpl/bx(180),sppj
      common/ijcpl/ajtarg,ajin
      common/kkg2/fractot(99),cs90
c--------------------random select incident a.m 'l'-----------------
c sppj is projectile spin
       jzl=0
       spinu=0.
       x=rndm(-1.)
c select an incoming l value-partial orbital a.m.=spinu
        do 987 lr=1,90
          jzl=lr-1
         if(x.gt.fractot(lr))go to 987
           spinu=float(jzl)
            if(spinu.gt.89.)spinu=89.
               go to 988
  987 continue
  988 continue
c----jzl is actual incident l, lr is index (1 larger than l)-------
       ajin=float(jzl)
         ajcpl=ajin
       if(sppj.eq.0.)go to 10205
           
c----- for a.m gt.0 must couple proj. spin if gt.0 with l---
c------below only set for spin 0.5 of projectile----------
       if(jzl.eq.0.and.sppj.eq.0.5)ajin=0.5
       if(jzl.eq.0.and.sppj.eq.1.0)ajin=1.0
c above 2 lines couple 0 -l incident with spin 1/2 or 1 projectile
       if(jzl.gt.0)then
        ajin=abs(float(jzl))-sppj
         x=rndm(-1.)
          if(x.gt.bx(jzl))then
          ajin=abs(float(jzl)+sppj)
          endif
         endif
10205      spinu=ajin
      return
      end
c------------------------------------------------------------------------
      subroutine spinch
      common/par3/ep,sigml(999),acrs(999)
      common/kkg2/fractot(99),cs90
c sub to make running sum partial projectile cs to select incident ch. l
              cs90=0.
           do  lr=1,90
           fractot(lr)=0.
           cs90=cs90+sigml(lr)
           enddo
           fractot(1)=sigml(1)/cs90
             do lr=2,90
             fractot(lr)=fractot(lr-1)+sigml(lr)/cs90
              enddo
       return
       end
c
      subroutine litup(li,excit)
c      formal parameter li was added, KKG, 04/08/10 
c written May 02;the cm particle energies/angles in this crude light particle
c decay routine have not yet been recoil coupled to give proper cm-lab conversion
      common/outlim/limout,lim2
      common/levop/ldopt,igam
      common/ismrdat2/ididl(106,260),nolevs(106,260)
      common/perp/spinprp
      common/isomer1/xgs(690),xiso(690),x2iso(690),isonum,isoint,iso2num
      common/isomer2/zs(690),as(690),eiso(690),eiso2(690)
      common/isomer3/iziso(690),iaiso(690),csgs(690),csis(690)
     1,csgse(690,100),csise(690,100),csis2e(690,100),rm1(690,100),
     2rm2(690,100)
      common/lin/lspin
      common/lscpl/bx(180),sppj
      common/ijcpl/ajtarg,ajin
      common/tagspin/trgspin
      common/tgspin/levno
      common/spins2/spinu,spinf,prz1
      common/rotor/irotor,iidummy,arote(15,25,181)
      common/history/ide(60),idn(60),idith(60),iddum(60),the(60),adelj
     1(60),spinfh(60),enghis(60)
       common/spins4/eprot,ealpha,eoff(8)
      common/spins/rote(15,24,181),popej(15,24,180,100),fcx(180),ai,ecmp
      common/sping/pfpej(600,180,100),pfpj(600,180)
      common/gate1/spind(5,180),popj(15,24,180),egatl(5),egatu(5),npart(
     15),igate
      common/nuopt/kinem,iem,kb
      common/excid/npartt(2)
      common/tem/philabb(36),phirec(36),phb(36),phin(36)
      common/bound/neut2,iprot,idnuc(100),idm(100),einc(100),
     1thett(100),phit(100),ppp
      common/hyb2/pp(15,24,999)
      common/evapmc/jcount,ipremc,eminn(15,24)
      common/ilcom/nu,ndum,am(8)
      common/rcl/amis
      common/labcon/prx,pry,prz,pcm,rmas,clab(999,8),dslb(36,999,8)
      common/exclbe/idz(4),ida(4),beff(20,20,2),pairmc(18,27)
      common/scat/jz,ja, probxx(2,2),q1(8,1100),fcs,eloss,fcspec
        common/mccom/ncount,ipot,rum(15,24),ppmc(15,24,999),cf
      common/exclus/ gib(8,999),pr(3),ge(2100),zz(3),
     1delta(3),ts(3),pairx(8),asum(10),rate(8,999)
      common/shft/mc,mp,inver,ike,ipch,kq,qval,ap,at,zp,zt,cld,barfac
      common/pl3/jx,nepr,eb(100),rcsp(100),zee,amass,plex
      common/par3/ep,sigml(999),acrs(999)
      common/memo/scrs(999,8)
      common/incr/ed,tened,efinal
      common/sft5/exc(18,27),xmxx
      common/lym1/be(18,27,8),pair(18,27),xmas(18,27),delshl(18,27),amas
     1(18,27)
c     common/bym1/symb(18,27),symbp(18,27),symbs(18,27)
      common/distb/dsgp(36,999,8)
      common/dist/n,ni,ke,irefr,xz(3),pod,pcs,ek,xmix,cfrac
     1,dtau,sign(36),tmix,bji(10),fin,dtup
      common/error/rerror(15,24)
      common/eder/ed2,v,vhalf,vsq
      common/angprm/aprp(500),eav3(1000),
     1eav2(1000),pot2x,thet1,phi1,snt(360),cst(360),snp(720),csp(720)
     2,pot
crc
      common/out1/nout
      common/out2/vrl(200,25),vrlz(200,25),vrlx(200,25),vrly(200,25)
      common/out3/vscale
      common/out4/asm(36),asmm(36)
crc
      common/recoils/rcletza(100,13,24,18),rcleza(100,13,24),
     1rclea(100,36),rcleat(100,36,18),rphir(36)
     2,rclsp(200),rclz(200),rcly(200),rclx(200),escale
     3,rcle(200),rclt(36),rclet(200,36)
      common/rec1/ermxc,thetrc
      common/mmtcom/mmt,mzee,iares,izres,light
      light=1
c light=1 indicates that this sub has been called
      if(izres.gt.2)go to 200
c here treat helium isotopes
      if(izres.eq.1)go to 150
c for now,we ignore 7He and 8He
c  He6
      if(iares.eq.6)then
      ealph=(excit-1.)/3.
      e2nut=2.*(excit-1.)/3.
       if(e2nut.lt.0.)return
      ix=(e2nut*0.5+ed)/ed
      ith=rndm(-1.)*35.9999+1.
c set theta as random angle for present in light breakup
       scrs(ix,1)=scrs(ix,1)+fcspec*2.
       dsgp(ith,ix,1)=dsgp(ith,ix,1)+fcspec*2.
       ix=(ealph+ed)/ed
      ith=rndm(-1.)*35.9999+1.
       scrs(ix,3)=scrs(ix,3)+fcspec
       dsgp(ith,ix,3)=dsgp(ith,ix,3)+fcspec
      return
      endif
c  He5
      if(iares.eq.5)then
      ealph=0.2*(excit+0.9)
      eneut=0.8*(excit+0.9)
      ix=(eneut+ed)/ed
      ith=rndm(-1.)*35.9999+1.
c set theta as random angle for present in light breakup
       scrs(ix,1)=scrs(ix,1)+fcspec
       dsgp(ith,ix,1)=dsgp(ith,ix,1)+fcspec
       ix=(ealph+ed)/ed
      ith=rndm(-1.)*35.9999+1.
       scrs(ix,3)=scrs(ix,3)+fcspec
       dsgp(ith,ix,3)=dsgp(ith,ix,3)+fcspec
      return
      endif
c  He4
c temporarily,for excited alpha,we assume all gamma decay
      if(iares.eq.4)then
      ix=1
      ith=1
       scrs(ix,3)=scrs(ix,3)+fcspec
       dsgp(ith,ix,3)=dsgp(ith,ix,3)+fcspec
      return
      endif
c  He3
      if(iares.eq.3)then
      edeut=(excit-5.5)/3.
      eprot=(excit-5.5)*2./3.
      if(excit.lt.5.5)return
      ix=(edeut+ed)/ed
c     iy=(eprot-5.5)/ed
      iy=(eprot+ed)/ed
      ith=rndm(-1.)*35.9999+1.
       scrs(ix,4)=scrs(ix,4)+fcspec
       dsgp(ith,ix,4)=dsgp(ith,ix,4)+fcspec
 
      ith=rndm(-1.)*35.9999+1.
       scrs(iy,2)=scrs(iy,2)+fcspec
       dsgp(ith,iy,2)=dsgp(ith,iy,2)+fcspec
       return
       endif
       return
c above return for case z=2,but isotope not treated in 'if' tests
  150 continue
c here treat isotopes of hydrogen
c  H4
      if(iares.eq.4)then
      eneut=0.75*(excit+2.82)
      etrit=0.25*(excit+2.82)
c     ix=(etrit+ed)/ed
      iy=(eneut+ed)/ed
c     ith=rndm(-1.)*35.9999+1.
 
      ith=rndm(-1.)*35.9999+1.
       scrs(iy,1)=scrs(iy,1)+fcspec
       dsgp(ith,iy,1)=dsgp(ith,iy,1)+fcspec
       return
       endif
c  H3
      if(iares.eq.3)then
       if(excit.le.6.26)return
      eneut=0.67*(excit-6.26)
      edeut=0.33*(excit-6.26)
      ix=(edeut+ed)/ed
      iy=(eneut+ed)/ed
      ith=rndm(-1.)*35.9999+1.
       scrs(ix,4)=scrs(ix,4)+fcspec
       dsgp(ith,ix,4)=dsgp(ith,ix,4)+fcspec
 
      ith=rndm(-1.)*35.9999+1.
       scrs(iy,1)=scrs(iy,1)+fcspec
       dsgp(ith,iy,1)=dsgp(ith,iy,1)+fcspec
       return
       endif
c  H2
      if(iares.eq.2)then
       if(excit.le.2.22)return
      eneut=0.50*(excit-2.22)
      eprot=0.50*(excit-2.22)
      ix=(eprot+ed)/ed
      iy=(eneut+ed)/ed
      ith=rndm(-1.)*35.9999+1.
       scrs(ix,2)=scrs(ix,2)+fcspec
       dsgp(ith,ix,2)=dsgp(ith,ix,2)+fcspec
 
      ith=rndm(-1.)*35.9999+1.
       scrs(iy,1)=scrs(iy,1)+fcspec
       dsgp(ith,iy,1)=dsgp(ith,iy,1)+fcspec
       return
       endif
 
  200 continue
      return
      end
c
c   **********************************************************************
c
      subroutine tfer10m
      common/mmtcom/mmt,mzee,ifragmas,ifragzee,light
      COMMON /PRODA/ PROD(7,20)
c
c
      OPEN(16,FILE='testf10m.out',STATUS='UNKNOWN')
c
      ifragmas=12
      ifragzee=6
      fragexc=100.
      nrun=1000
      do  n=1,nrun
      call  gudima(fragexc,KST)
        write(16,*) 'n,KST=',n,KST
        do  k=1,kst
       write(16,101) (prod(j,k),j=1,7)
  101     format(7(1PE11.4))
        enddo
      enddo
      stop
       end
c	  		 
c FERMI STATISTICS ROUTINE FROM K.K.GUDIMA
      subroutine gudima(fragexc,KST)
      common/outlim/limout,lim2
c
c This routine with associated routines was provided by Prof. K.K.
c Gudima;some discussion of physics is to be found in report
c LA-UR-01-6804 by K.K.Gudima,S.G.Mashnik,and A.J.Sierk
c updated by K.K. Gudima March 14, 2010
c
c      test
c     produced fragments are in PROD:
c          i                          PROD(i,k)
c     --------------------------------------------------------
c          1                         mass number
c          2                         charge
c          3                          Px (MeV/c)
c          4                          Py (MeV/c) Momentum
c          5                          Pz (MeV/c)
c          6                          Ek (MeV), kinetic energy
c          7                          Mk (MeV), mass
c     --------------------------------------------------------
c     KST - number(multiplicity) of produced  fragments
c     --------------------------------------------------------
      common/mmtcom/mmt,mzee,ifragmas,ifragzee,light
      common/fermsum/asum,zsum
      COMMON /PRODA/ PROD(7,20)
C
        real*8 a,z,af,zf,u,uf,pn,pf
      DIMENSION PN(3),psum(3),pf(3),prod0(7,20)
      save ifirstf
      data ifirstf/1/,thsnd/1.0d3/
c
c     OPEN(16,FILE='testf.out',STATUS='UNKNOWN')
c
cx    call RDMINI         !random number initialization
      if(ifirstf.eq.1)  then
        call GITAB          !Fermi Decay initialization
        ifirstf=0
      endif
c
 
      A=12.             ! Mass number of excited nucleus
      Z=6.              ! Charge of excited nucleus
      PN(1)=0.          ! GeV/c
      PN(2)=0.          ! Gev/c  Momentum of excited nucleus
c     PN(3)=0.1         ! GeV/c
      PN(3)=0.0         ! GeV/c
c     U=0.100           ! Exitation energy in GeV
      U=fragexc/1000.
c convert alice MeV excitations to GeV
c
      NW=1
c
c     NRUN=100          ! Number of simulations
c     do  n=1,NRUN
c     A=ifragmas
c     Z=ifragzee
      a=ifragmas
      z=ifragzee
        asum=0
        zsum=0
        do  k=1,3
          psum(k)=0.
        enddo
        KST0=0
        KST1=0
c       CALL  RASTAR(A,Z,U,PN,KST1,NW)
         call rastar(a,z,u,pn,kst1,nw)
c       mv = KST1
        mv = kst1 
        do  m = 1,mv
          do  k = 1,7
          prod0(k,m) = prod(k,m)
          enddo
        enddo
   10 continue
      ntry = 0 
      m1 = mv
        do  m = 1,mv
        af = prod0(1,m)
        zf = prod0(2,m)
        iaf = nint(af)
        izf = nint(zf)
c All allowed intermediate states should be included below. More are
c included here than are presently allowed.       
        if ((iaf.eq.5 .and. (izf.eq.2.or.izf.eq.3)) .or.
     &      (iaf.eq.6 .and. izf.ge.4) .or.  
     &      (iaf.eq.7 .and. izf.ge.2) .or.          
     &      (iaf.eq.8 .and. (izf.eq.4.or.izf.eq.6)) .or.   
     &      (iaf.eq.9 .and. izf.eq.5).or.                ! KKG 12/17/07
     &      (iaf.eq.9 .and. izf.eq.6))     then
          pf(1) = prod0(3,m)/thsnd  
          pf(2) = prod0(4,m)/thsnd  
          pf(3) = prod0(5,m)/thsnd
          uf = 1.0d-6
          ksf = 0  
cc        write (16, *) ' Fragment m, af, zf = ', m, af, zf
c nov 8,2010 add return on af or zf lt.1
          if(af.lt.1..or.zf.lt.1.)return
          call rastar (af, zf, uf, pf, ksf, 1)
          if (ksf.gt.1)   then      ! kkg 21.02.09
            ntry = 1
            do mf = 1,ksf
cc          write (16, *) ' Decay mf, af, zf = ', 
cc   &                            mf, prod(1,mf), prod(2,mf) 
              do  k = 1,7
                if (mf.eq.1) then
                  prod0(k,m) = prod(k,mf)
                else 
                  if (k.eq.1) m1 = m1 + 1
                  prod0(k,m1) = prod(k,mf)
                endif
              enddo
            enddo
          endif 
        endif
        enddo
      mv = m1
      if (ntry.gt.0) go to 10
      KST=mv                   ! number of produced fragments
        do m = 1,mv
          do k = 1,7
          prod(k,m) = prod0(k,m)
          enddo
        enddo 
c    
        do  k=1,mv
          asum=asum+PROD(1,k)
          zsum=zsum+PROD(2,k)
          psum(1)=psum(1)+PROD(3,k)/1000.
          psum(2)=psum(2)+PROD(4,k)/1000.
          psum(3)=psum(3)+PROD(5,k)/1000.
          if(limout.eq.0)then
            write(29,100) k,(PROD(l,k),l=1,7)
  100       format(1x,I3,2F5.0,1x,5(1PE11.4))
          endif
        enddo
c       write(* ,101) asum,zsum,psum,A,Z,PN
       if(limout.eq.0)then
        write(29,101) asum,zsum,psum,A,Z,PN
  101   format(4x,77('-')/' sum',2F5.0,1x,3(1PE11.4)/
     &         ' ini',2(0PF5.0),1x,3(1PE11.4)/4x,77('-')/)
       endif
c     enddo
c     stop
      return
      end
 
C**********************************************************************
      subroutine rastar (ap, zp, up, pn, nst, nw)

c ======================================================================
c
c    Last change: 13-Aug-2003 BY NVMokhov
c    Edited by A. J. Sierk, LANL T-16, September, 2003.
c
c ======================================================================

      implicit real*8 (a-h, o-z), integer (i-n)
      real*4 yz

c ======================================================================

      common /ipaz/  iaz(8000) 
      common /proda/ prod(7,20)
      common /spdm/  dm(12,11), ms(12,11)
      common /wnq/   wnq(2000)
cc    common /yield/ ysp(12,11)

      dimension pn(3), vn(3), pnc(3), pnl(3), pn0(3,20)
      dimension ipa(20), ipz(20), amk(20), nc(20)

      data zro, thsnd /0.d0, 1.d3/

c ======================================================================

      izp = nint(zp) + 1
      inp = nint(ap - zp) + 1
c  kkg  08/11/04
      if(izp.lt.1.or.izp.gt.11.or.inp.lt.1.or.inp.gt.12)  then
      cmp = 0.938d0*ap
      else
        cmp = 0.9314943d0*ap + 0.001d0*dm(inp,izp) - 0.000511d0*zp
      endif
      pq = pn(1)**2 + pn(2)**2 + pn(3)**2
      eap = sqrt(pq + cmp**2)
c  kkg  31.01.09
        eap = eap + up
c
        do k = 1,3
        vn(k) = pn(k)/eap
        end do
      if (ap.ge.1.9d0) then
        call razval (ap, zp, up, nf, npt)
        if (nf.ge.1) then
          ia = nint(ap)
          iz = nint(zp)
          jaz = 100*ia + iz
          daz = wept (ia, iz) + up
            do l = 1,nw
            x = rndm(-1.)
            sr = zro
              do n = 1,nf
              nb = n
              sr = sr + wnq(n)
              if (br.lt.sr) go to 10
              if (x.lt.sr) go to 10
              end do
   10       nc(l) = nb
            end do
          nv = nw
          k = 0
          msa = 0
          ng = 0
            do n = 1,npt
            k = k + 1
            ipa(k) = iaz(n)/100
            ipz(k) = iaz(n) - 100*ipa(k)
            msa = msa + iaz(n)
            if (msa.ge.jaz) then
              ng = ng + 1
cc  Unused!  AJS
cc            wq = wnq(ng)
cc              do i = 1,k
cc              jp = ipz(i) + 1
cc              in = ipa(i) - jp + 2
cc              ysp(in,jp) = ysp(in,jp) + wq
cc              end do
   20         if (nv.ge.1) then
                  do l = 1,nv
                  lv = l
                  if (nc(l).eq.ng) go to 30
                  end do
                go to 40
   30           nv = nv - 1
                if (nv.ge.1) then
                  do  l = lv,nv
                  nc(l) = nc(l+1)
                  end do
                endif
                tn = daz
                  do i = 1,k
                  mv = nst + i
                  if (mv.gt.20) write (16, 100) mv
                  amk(i) = wept (ipa(i), ipz(i))
                  tn = tn - amk(i)
                  end do
                call disimp (k, amk, pn0, tn)
                  do i = 1,k
                    do j = 1,3
                    pnc(j) = pn0(j,i)
                    end do
                  en = sqrt(pnc(1)**2 + pnc(2)**2 + pnc(3)**2 + 
     &                      amk(i)**2)
                  call clpv (pnc, vn, pnl, en)
                  call pint (pnl, ct, st, cf, sf, tk, amk(i))
                  mv = nst + i
                  prod(1,mv) = dble(ipa(i))
                  prod(2,mv) = dble(ipz(i))
                  prod(3,mv) = pnl(1)*thsnd
                  prod(4,mv) = pnl(2)*thsnd
                  prod(5,mv) = pnl(3)*thsnd
                  prod(6,mv) = tk*thsnd
                  prod(7,mv) = amk(i)*thsnd      
                  end do
                nst = mv
                go to 20
              endif
   40         k = 0
              msa = 0
            endif
            end do
          return
        endif
      endif
      m1 = nst + 1
      nst = m1
      if (ap.gt.0.1d0) then
c 8 nov 2010
        if(izp.lt.1)then
         izp=1
         zp=1.
         endif
          if(inp.lt.1)inp=1
c end 8 nov
        cmn = 0.9314943d0*ap + 0.001d0*dm(inp,izp) - 0.000511d0*zp
        call pint (pn, ct, st, cf, sf, tmn, cmn)
        prod(1,m1) = ap
        prod(2,m1) = zp
        prod(3,m1) = pn(1)*thsnd
        prod(4,m1) = pn(2)*thsnd
        prod(5,m1) = pn(3)*thsnd
        prod(6,m1) = tmn*thsnd
        prod(7,m1) = cmn*thsnd      
      endif
      return

c ======================================================================

  100 format (5x,'mv > 20 in RASTAR; mv = ',i5)
cc200 format (2x,3(e10.3,1x),2i5)
cc201 format (2x,33('-')/2x,3(e10.3,1x)/)

c ======================================================================
      end


      subroutine disimp (k, amk, pn0, tn)

c ======================================================================
c
c    Modified to incorporate ISOTR subroutine inline.
c
c    Last change: 13-Aug-2003 BY NVMokhov
c    Modified by A. J. Sierk, LANL T-16, October, 2003.
c    Corrected if(rndm(-1.) -> rnd=rndm(-1.), if(rnd...  by KKG&NVM 08/14/09
c
c ======================================================================

       implicit real*8 (a-h, o-z), integer (i-n)
c ======================================================================

      dimension anl(3), pnc(3), vrs(3), pn0(3,20), amk(20)

      common /pi2pi/ pi, twpi

      data zro, one, two /0.d0, 1.d0, 2.d0/

c ======================================================================

      smk = zro
        do i = 1,k
        smk = smk + amk(i)
        end do
        do i = 1,3
        vrs(i) = zro
        end do
      tkn = tn
      kl = k - 1
        do l = 1,kl
        lk = k - l + 1
        amp = amk(lk)
        smk = smk - amp
        amr = smk
        tpr = tkn
        if (lk.ge.3) then
   10     csi = rndm(-1.)
          kk = lk - 1
          if (kk.ge.3) then
            pex = one/(1.5d0*dble(kk) - two)
            csi = csi**pex
          endif
          fk = two*sqrt(csi*(one - csi))
          rnd=rndm(-1.)
          if (rnd.gt.fk) go to 10
          rnksi = csi
          tkm = tkn*rnksi
          tpr = tkn - tkm
        endif
        pmc = sqrt(two*((amp*amr)/(amp + amr))*tpr)
        ct = one - two*rndm(-1.)
        st = sqrt(abs(one - ct**2))
        anl(3) = ct
        fi = twpi*rndm(-1.)
        anl(1) = cos(fi)*st
        anl(2) = Sin(fi)*st
          do i = 1,3
          pnc(i) = pmc*anl(i)
          pn0(i,lk) = vrs(i)*amp + pnc(i)
          vrs(i) = vrs(i) - pnc(i)/amr
          end do
        tkn = tkm
        end do
        do i = 1,3
        pn0(i,1) = vrs(i)*amr
        end do
      return

c ======================================================================
      end
      subroutine pint (p, ct, st, cf, sf, t, cm)

c ======================================================================
c
c   Defining angles from the 3-momentum components
c
c    Last change: 13-Aug-2003 BY NVMokhov
c    Edited by A. J. Sierk, LANL T-16, September, 2003.
c
c ======================================================================

      implicit real*8 (a-h, o-z)

c ======================================================================

      dimension p(3)

      data zro, one /0.d0, 1.d0/

c ======================================================================

      pz = p(3)**2
      pq = p(1)**2 + p(2)**2 + pz
      ctq = pz/pq
      t = sqrt(pq + cm**2) - cm
      if (ctq.ge.one) then
        st = zro
        ct = one
        sf = zro
        cf = one
      else
        ct = sqrt(ctq)
        if (p(3).le.zro) ct = -ct
        st = sqrt(one - ctq)
        pmh = st*sqrt(pq)
        cf = p(1)/pmh
        sf = p(2)/pmh
      endif
      return

c ======================================================================
      end
      subroutine clpv  (p1, v, p2, e1)

c ======================================================================
c
c   Transforms momentum p1 in frame moving with velocity v with respect
c   to the lab frame into the lab frame momentum p2.
c
c    Last change: 13-Aug-2003 BY NVMokhov
c    Edited by A. J. Sierk, LANL T-16, September, 2003.
c
c ======================================================================

      implicit real*8 (a-h, o-z), integer (i-n)

c ======================================================================

      dimension p1(3), v(3), p2(3)

      data zro, one /0.d0, 1.d0/

c ======================================================================

      spv = zro
      v2 = zro
        do i = 1,3
        v2 = v2 + v(i)**2
        spv = spv + p1(i)*v(i)
        end do
      gam = one/sqrt(one - v2)
      temp = gam*(spv*gam/(gam + one) + e1)
        do k = 1,3
        p2(k) = p1(k) + v(k)*temp
        end do
      return

c ======================================================================
      end

       FUNCTION wechan (ja, k, ntv, daz, tq)
c     DOUBLE PRECISION FUNCTION wechan (ja, k, ntv, daz, tq)

c ======================================================================
c
c   Last change: 13-Aug-2003 BY NVMokhov
c   Edited by A. J. Sierk, LANL T-16, September, 2003.
c   Modified 03 May 2006 by R E Prael to make the Coulomb barrier slightly
c      penetrable (purely an artificial procedure that reduces the channel
c      probability by 10.d-06 below the Coulomb energy.) Necessary for the 
c      breakup of Be8 at low excitation.
c
c   Last change: 14-Aug-2009 BY NVMokhov
c ======================================================================

      implicit real*8 (a-h, o-z), integer (i-n)

c ======================================================================

      common /gaf/    gaf(20)
      common /kappa4/ akp
      common /spdm/   dm(12,11), ms(12,11)
      common /vak/    vak

      dimension ntv(20)

      data zro, one /0.d0, 1.d0/
      save ifirst,exp1
c ======================================================================

      exp1 = exp(one)
      bnq = daz
      rmq = -0.1d0
      spm = one
      akm = one
      akv = dble(ja)
      if (akp.gt.zro) akv = akp*dble(k)
      akv = akv*vak
        do i = 1,k
        ia = ntv(i)/100
        iz = ntv(i) - ia*100
        akm = akm*dble(ia)
        spm = spm*dble(ms(ia-iz+1,iz+1))
        bnq = bnq - wept (ia, iz)
        end do
      tn = bnq
      if (tn.gt.zro) then
        akm = akm/dble(ja)
        akm = akm*sqrt(akm)*spm
        if (k.gt.2) then
          ki = k - 1
          rpm = one
          vmk = one
cc        teq = 2.71828182846d0*tn/(1.5d0*dble(k) - 2.5d0)
cc        teq = exp1*tn/(1.5d0*dble(k) - 2.5d0)
          teq = exp1*tn/(1.5d0*dble(ki) - one)
          vtk = akv*teq*sqrt(teq)
            do i = 1,ki
            mrs = 1
            ik = i + 1
            vmk = vmk*vtk
              do j = ik,k
              if (ntv(i).eq.ntv(j)) mrs = mrs + 1
              end do
            rpm = rpm*dble(mrs)
            end do
          gam = gaf(k)
          rmq = vmk*gam*akm/(teq*rpm)
        else
          rmq = 1.1283792d0*akv*akm*sqrt(tn)
          if (ntv(1).eq.ntv(2)) rmq = 0.5d0*rmq
        endif
        if (tn.lt.tq) rmq = rmq*1.d-06
      endif
      wechan = rmq
      return

c ======================================================================
      end

      subroutine gitab 

c ======================================================================
c
c    Last change: 13-Aug-2003 BY NVMokhov
c    Modified by A. J. Sierk, LANL T-16, September-November, 2003.
c
c ======================================================================

c     implicit real*8 (a-h, o-z), integer (i-n)

c ======================================================================

      common /constan8/ thrd, twthrd
      common /gaf/     gaf(20)
      common /kappa4/  akp
      common /mfp/     mf, mp
      common /pi2pi/      pi, twpi
      common /spdm/    dm(12,11), ms(12,11)
      common /vak/     vak
c
      common /ato3rd8/  a3rd(300)

c     dimension ipa(12), jpz(12), dmp(12), msp(12), amp(12)
      dimension  ms1(12,11), md(12,11)

      data  ms1/      !2*S + 1  with some empirical multipliers for
c                               common isotopes.
     &     0,  2,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 
     &     0,  3,  2,  0,  0,  0,  0,  0,  0,  0,  0,  0, 
     &     0,  2,  1,  4,  1,  0,  1,  0,  0,  0,  0,  0, 
     &     0,  0,  0,  4,  6,  5,  4,  0,  4,  0,  0,  0, 
     &     0,  0,  0,  6,  1,  4,  7,  2,  1,  0,  1,  0, 
     &     0,  0,  0,  5,  4, 14, 16, 18,  4,  5,  0,  0, 
     &     0,  0,  0,  4,  1,  4, 11, 14, 13,  8,  1,  0, 
     &     0,  0,  0,  0,  0,  3,  2, 18, 10,  5,  2,  3, 
     &     0,  0,  0,  0,  0,  4,  1, 10, 17,  6,  1,  6, 
     &     0,  0,  0,  0,  0,  0,  0,  1,  6,  3,  2,  5, 
     &     0,  0,  0,  0,  0,  0,  0,  2,  1,  2,  1,  4/

      data ((md(i,j),i=1,12), j=1,9)/      !Mass excesses in keV
     &   99000,  8071, 99999, 99999, 99999, 99999, 99999, 99999, 99999,
     &   99999, 99999, 99999, 
     &   99000, 13136, 14950, 25900, 36800, 41900, 99999, 99999, 99999,
     &   99999, 99999, 99999, 
     &   99000, 14931,  2425, 11390, 17594, 26110, 31598, 40820, 48810,
     &   65000, 75240, 89260, 
     &   99000, 25300, 11680, 14086, 14908, 20946, 24954, 33050, 40800,
     &   52940, 61570, 72280, 
     &   99000, 38000, 18375, 15769,  4942, 11348, 12607, 20174, 25080,
     &   33700, 39900, 51210, 
     &   99000, 99999, 27870, 22921, 12416, 12051,  8668, 13369, 16562,
     &   23660, 28970, 37080,
     &   99000, 99999, 35090, 28914, 15699, 10651,   000,  3125,  3020,
     &    9873, 13694, 21040,
     &   99000, 99999, 99999, 39700, 25300, 17338,  5345,  2863,   101, 
     &    5683,  7870, 13120, 
     &   99000, 99999, 99999, 42700, 32050, 23111,  8007,  2855,  4737,
     &     809,   782,  3334/

      data ((md(i,j), i=1,12),j=10,11)/
     &   99000, 99999, 99999, 99999, 39700, 33600, 16800, 10680,  1952,
     &     873,  1487,    17, 
     &   99000, 99999, 99999, 99999, 49400, 36400, 23990, 16490,  5307, 
     &    1751,  7042,  5732/

      data one, two /1.d0, 2.d0/

c ======================================================================
c
      pi = 3.141592653589793d0
      twpi = two*pi
      thrd = 1.d0/3.d0
      twthrd = 2.d0/3.d0
      do  k = 1,300        
          a3rd(k) = dble(k)**thrd
      enddo 
c
      mf = 1999
      mp = 8000
      akp = -one
c   NOTE:  proton number is j - 1
c   NOTE:  neutron number is i - 1
        do i = 1,12
          do j = 1,11
          ms(i,j) = ms1(i,j)
          dm(i,j) = 0.001d0*dble(md(i,j))
          end do
        end do
      ms(1,2) = 2
      dm(1,2) = 7.289d0
      dm(9,9) = -dm(9,9)
      dm(10,9) = -dm(10,9)
      dm(11,9) = -dm(11,9)
      dm(11,10) = -dm(11,10)
      dm(12,10) = -dm(12,10)
      dm(11,11) = -dm(11,11)
      dm(12,11) = -dm(12,11)
      vl = 1.4d0*sqrt(0.94d0)*5.07d0
      vol = vl*vl
      vol = vol*vl
      vak = vol*thrd*sqrt(two/pi)
      gaf(1) = one
      gaf(2) = one
        do k = 3,20
        qki = 1.5d0*dble(k-1) - one
        qk = one/qki
        gq = one + qk*(one + qk*(one - qk*139.d0/180.d0)/24.d0)/12.d0
        gaf(k) = sqrt(qk*0.1591549d0)/gq
        end do
c   Unused quantities; AJS 11/4/03.
c       do j = 1,11
c         do i = 1,12
c         jpz(i) = j - 1
c         ipa(i) = j + i - 2
c         dmp(i) = dm(i,j)
c         msp(i) = ms(i,j)
c         amp(i) = 0.9314943d0*dble(ipa(i)) + 0.001d0*dmp(i) -
c    &             0.000511d0*dble(jpz(i))
c         end do
c       end do
c
      return

c ======================================================================
      end

       FUNCTION wept (ia, iz)
c     DOUBLE PRECISION FUNCTION wept (ia, iz)

c ======================================================================
c
c    Last change: 13-Aug-2003 BY NVMokhov
c    Edited by A. J. Sierk, LANL T-16, September, 2003.
c    Changed by KKG 25 Apr 2006  
c    Modified 03 May 2006 by R E Prael to set mass of a "proton star",
c       thus allowing breakup under all conditions.
c
c ======================================================================

      implicit real*8 (a-h, o-z), integer (i-n)

c ======================================================================

      common /spdm/ dm(12,11), ms(12,11)

c ======================================================================

      j = iz
      i = ia - j
      if(i.eq.0.and.j.gt.1) then      ! ia=iz
c    Modified 03 May 2006 by R E Prael (IA*proton mass)      
        wept=dble(ia)*(0.9309833d0 + 0.001d0*dm(1,2)) 
      elseif (i.lt.0 .or. i.gt.11 .or. j.le.0 .or. j.gt.10)  then
c Changed by KKG 25 Apr 2006    
ccc, SGM, 06/08/06 write (*, *) ' In WEPT : ia, iz = ', ia, iz ! 05/31/02 KKG
        wept = 0.938d0*dble(ia)
      else
        dmp = dm(i+1,j+1)
        if (dmp.le.90.d0) then
          wept = 0.9314943d0*dble(ia) + 0.001d0*dmp - 0.000511d0*dble(j)
        else
          wept = 0.938d0*dble(ia)
        endif
      endif
      return

c ======================================================================
      end

      subroutine crack (k, ia, iz, np)

c ======================================================================
c
c    Last change: 13-Aug-2003 BY NVMokhov
c    Edited by A. J. Sierk, LANL T-16, September, 2003.
c
c ======================================================================

c     implicit real*8 (a-h, o-z), integer (i-n)

c ======================================================================

      common /mpraz/ mpa(200), mpz(200), mra(200), mrz(200)
      common /ipaz/  iaz(8000)

      dimension maz(20)

c ======================================================================

      kp = k
      nw = 0
      call divaz (ia, iz, ia, iz, kp, nr)
      if (nr.lt.1) return
      nk = np - kp
        do n = 1,nr
        nk = nk + kp
        iaz(nk+1) = 100*mpa(n) + mpz(n)
        iaz(nk+2) = 100*mra(n) + mrz(n)
        nw = nw + 1
        end do
      nn = kp*nw + np
      if (kp.lt.3) then
        np = nn
      else
        nq = nw
        nw = 0
        nc = 0
        mk = kp - 1
          do m = 2,mk
          km = kp - m + 1
          nk = np - kp
          mi = m - 1
   10     continue
          nk = nk + kp
          nm = nk + m
   20     continue
            do i = 1,mi
            maz(i) = iaz(nk+i)
            end do
          la = iaz(nm-1)/100
          lz = iaz(nm-1) - la*100
          ka = iaz(nm)/100
          kz = iaz(nm) - ka*100
          call divaz (la, lz, ka, kz, km, nr)
          if (nr.lt.1) go to 30
          iaz(nm) = 100*mpa(1) + mpz(1)
          iaz(nm+1) = 100*mra(1) + mrz(1)
          nw = nw + 1
          if (nr.ge.2) then
            do n = 2,nr
              do i = 1,mi
              iaz(nn+i) = maz(i)
              end do
            nm = nn + m
            iaz(nm) = 100*mpa(n) + mpz(n)
            iaz(nm+1) = 100*mra(n) + mrz(n)
            nn = nn + kp
            nw = nw + 1
            end do
          endif
          nc = nc + 1
          if (nc.lt.nq) go to 10
          go to 40
   30     nq = nq - 1
          nn = nn - kp
          if (nq.lt.1) then
            np = nn
            return
          endif
          n1 = nk + 1
          n2 = nn
            do ni = n1,n2
            iaz(ni) = iaz(ni+kp)
            end do
          if (nc.lt.nq) go to 20
   40     nq = nw
          nc = 0
          nw = 0
          end do
        np = nn
      endif
      return

c ======================================================================
      end

       subroutine divaz (la, lz, ia, iz, k, nr)

c ======================================================================
c
c     Last change: 13-Aug-2003 BY NVMokhov
c   Edited by A. J. Sierk, LANL T-16, September, 2003.
c   map replaced with mak by KKG&NVM 08/14/09
c
c ======================================================================

c     implicit real*8 (a-h, o-z), integer (i-n)

c ======================================================================

      common /mpraz/ mpa(200), mpz(200), mra(200), mrz(200)
      common /spdm/  dm(12,11), ms(12,11)

c ======================================================================

      npr = 200
      m = 0
      nmx = la
      nmk = ia - k + 1
      nmx = min (nmx, nmk)
      nmn = (ia - 1)/k + 1
      nmn = min(nmx, nmn)
      nm = nmx - nmn + 1
        do n = 1,nm
        mak = nmx - n + 1
        mar = ia - mak
        j1 = iz - mar
        j1 = max (j1,0)
        jz = 0
        if (mak.eq.la) jz = lz
        j1 = max(j1, jz)
        j2 = iz
        if (k.lt.3 .and. mak.eq.mar) j2 = iz/2
        j2 = min (mak, j2)
        j1 = j1 + 1
        j2 = j2 + 1
        if (j2.ge.j1) then
          do j = j1,j2
          mzp = j - 1
          mzr = iz - mzp
          jp = mzp + 1
          ip = mak - mzp + 1
c kkg 11/02/04
          if(ip.le.0.or.jp.le.0.or.ip.gt.12.or.jp.gt.11) go to 10  
          if (ms(ip,jp).ge.1) then
            if (k.le.2) then
              jp = mzr + 1
              ip = mar - mzr + 1
c kkg 11/02/04
c             if (ms(ip,jp).lt.1) go to 10
              if(ip.le.0.or.jp.le.0.or.ip.gt.12.or.jp.gt.11) go to 10  
            endif
            m = m + 1
            if (m.gt.200) then
              write (16, 100) npr, n
              return
            endif
            mpa(m) = mak
            mpz(m) = mzp
            mra(m) = mar
            mrz(m) = mzr
          endif
   10     continue
          end do
        endif
        end do
      nr = m
      return

c ======================================================================

  100 format (10x,'Array size is exceeded in DIVAZ; ',i3,'Particles ',
     &         3x,'n = ',i3)

c ======================================================================
      end

      subroutine razval (an, zn, up, nf, np)

c ======================================================================
c
c    Last change: 13-Aug-2003 BY NVMokhov
c    Edited by A. J. Sierk, LANL T-16, September, 2003.
c    Modified by R E Prael, 03 May 2006.
c
c ======================================================================
c 10-29-14 remove c from line below as possible real8 errr

c     implicit real*8 (a-h, o-z), integer (i-n)
      implicit real*8 (a-h, o-z), integer (i-n)

c ======================================================================

      common /ipaz/ iaz(8000)
      common /wnq/  wnq(2000)
      common /mfp/  mf, mp

      dimension ntv(20)

      data zro /0.d0/

c ======================================================================

      ia = nint(an)
      mpt = mp - ia
      iz = nint(zn)
      daz = up + wept (ia, iz)
      nf = 0
      np = 0
      srq = zro
      ka = ia - 1
      if (ka.gt.2) then
        do k = 2,ka
        kp = k
        nt = np
        nk = nt - kp
        call crack (kp, ia, iz, np)
        if (np.gt.nt) then
   10     continue
          nk = nk + kp
          if (nk.gt.mpt) then
            write (16, 100) ia, iz, kp, nf, nt, nk, np
            return
          endif
   20       do i = 1,k
            ntv(i) = iaz(nk+i)
            end do
          tq = tcul (kp, ntv)
          rmq = wechan (ia, kp, ntv, daz, tq)
          if (rmq.gt.zro) then        ! REP 03 May 2006         
            nf = nf + 1
            if (nf.gt.mf) then
              write (16, 100) ia, iz, kp, nf, nt, nk, np
              return
            endif
            srq = srq + rmq
            wnq(nf) = rmq
            if ((np - nk).gt.kp) go to 10
          else
            np = np - kp
            ni = nk + 1
            if (np.ge.ni) then
                do n = ni,np
                iaz(n) = iaz(n+kp)
                end do
              go to 20
            else
              if ((np.le.nt) .and. (kp.gt.2) .and. srq.gt.zro) then ! REP 03 May 2006       
                if (nf.ge.1) then
                  do n = 1,nf
                  wnq(n) = wnq(n)/srq
                  end do
                endif
                return
              endif
            endif
          endif
        endif
        end do
      endif
        do i = 1,ia
        jz = 1
        if (i.gt.iz) jz = 0
        ntv(i) = 100 + jz
        end do
      tq = tcul (ia, ntv)
      rmq = wechan (ia, ia, ntv, daz, tq)
      if (rmq.ge.zro) then
        nf = nf + 1
        wnq(nf) = rmq
        srq = srq + rmq
          do i = 1,ia
          iaz(np+i) = ntv(i)
          end do
        np = np + ia
      endif
      if (nf.ge.1) then
        do n = 1,nf
        wnq(n) = wnq(n)/srq
        end do
      endif
      return

c ======================================================================

  100 format (10x,'In RAZVAL, ia, iz, kp, nf, nt, nk, np = ',7(i4))
c 200 format (3x,'In RAZVAL An = ',1pe10.3,2x,', Zn = ',e10.3,2x,
c    &           ', and Up = ',e10.3)

c ======================================================================
      end

      FUNCTION tcul (k, ntv)
c     DOUBLE PRECISION FUNCTION tcul (k, ntv)

c ======================================================================
c
c     Last change: 13-Aug-2003 BY NVMokhov
c     Edited by A. J. Sierk, LANL T-16, September, 2003.
c
c ======================================================================

      implicit real*8 (a-h, o-z), integer (i-n)

c ======================================================================

      common /ato3rd8/  a3rd(300)
      common /constan8/ thrd, twthrd

      dimension ntv(20), ia(20), iz(20)

      data fkap /1.d0/

      data zro, one /0.d0, 1.d0/

c ======================================================================

      tcul = zro
      a = zro
      z = zro
      coef = (0.6d0*1.44d0/1.3d0)*(one/(one + fkap)**thrd) 
        do  i = 1,k
        ia(i) = ntv(i)/100 
        iz(i) = ntv(i) - ia(i)*100
        a = a + dble(ia(i))
        z = z + dble(iz(i))
        enddo
      ec = zro    
      iaa = nint(a)
        do  i = 1,k
        zi = dble(iz(i)) 
        ec = ec + zi*zi/a3rd(ia(i))
        enddo
      ec = z*z/a3rd(iaa) - ec
      tcul = 0.001d0*coef*ec
      tcul = max (zro, tcul)
cc
c     write (*, *) ' k, a, z, tcul = ', k, a, z, tcul 
cc
      return

c ======================================================================
      end
* INITIALIZATION:
* At initialization stage, IJKLIN,NTOTIN,NTOT2N must have the above default
* values 54217137, 0, 0
* You can read SEED card as described above to start with another seed,
* e.g., 64217136, 0 0
* If so, the following call must be made in your initialization routine:
*
*      IF(IJKLIN.NE.54217137.OR.NTOTIN.NE.0.OR.NTOT2N.NE.0) THEN
*         CALL RM48IN(IJKLIN,NTOTIN,NTOT2N)
*      END IF
* CALL of this DOUBLE PRECISION FUNCTION with REAL argument (= -1.)
*     x=rndm(-1.)
********************************************
      FUNCTION RNDM (RDUMMY)
c     IMPLICIT REAL*8 (A-H,O-Z)
c     REAL*8 RNDNUM,RDUMMY
c 1/28 make rndnum real8
      REAL*8 RNDNUM,RDUMMY
      DIMENSION RNDNUM (2)
c 1/28 temp c  out real8
      CALL RM48 ( RNDNUM, 1 )
      RNDM = RNDNUM (1)
      RETURN
      END
c ******************************************
      SUBROUTINE RDMINI
      ISEED1=54217137
      ISEED2=0
      ISEED3=0
      CALL RM48IN (ISEED1,ISEED2,ISEED3)
      RETURN
      END
c ******************************************
      SUBROUTINE RDMIN (ISEED1,ISEED2,ISEED3)
      REAL*8 RNDNUM
c 2/5 place real8 rndnum above as test
      write(*,*) ' ISEED1,2,3=',ISEED1,ISEED2,ISEED3
      if(ISEED1.lt.0.or.ISEED1.gt.900000000)  ISEED1=1234567
      if(ISEED2.lt.0.or.ISEED2.gt.999999999)  ISEED1=0
      if(ISEED3.lt.0.or.ISEED3.gt.900000000)  ISEED1=0
      CALL RM48IN (ISEED1,ISEED2,ISEED3)
      write(*,*) ' ISEED1,2,3=',ISEED1,ISEED2,ISEED3
      RETURN
c ******************************************
      ENTRY RDMOUT(ISEED1,ISEED2,ISEED3)
      CALL RM48UT (ISEED1,ISEED2,ISEED3)
      RETURN
      END
*************************************************************
      SUBROUTINE RM48(RVEC,LENV)
C     Double-precision version of
C Universal random number generator proposed by Marsaglia and Zaman
C in report FSU-SCRI-87-50
C        based on RANMAR, modified by F. James, to generate vectors
C        of pseudorandom numbers RVEC of length LENV, where the numbers
C        in RVEC are numbers with at least 48-bit mantissas.
*                                                                      *
*     REVISION: 21-JAN-1997    BY NVM                                  *
*                                                                      *
C   Input and output entry points: RM48IN, RM48UT.
C!!! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C!!!  Calling sequences for RM48:                                    ++
C!!!      CALL RM48 (RVEC, LEN)     returns a vector RVEC of LEN     ++
C!!!                   64-bit random floating point numbers between  ++
C!!!                   zero and one.                                 ++
C!!!      CALL RM48IN(I1,N1,N2)   initializes the generator from one ++
C!!!                   64-bit integer I1, and number counts N1,N2    ++
C!!!                  (for initializing, set N1=N2=0, but to restart ++
C!!!                    a previously generated sequence, use values  ++
C!!!                    output by RM48UT)                            ++
C!!!      CALL RM48UT(I1,N1,N2)   outputs the value of the original  ++
C!!!                  seed and the two number counts, to be used     ++
C!!!                  for restarting by initializing to I1 and       ++
C!!!                  skipping N2*100000000+N1 numbers.              ++
C!!! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

C for 32-bit machines, use IMPLICIT DOUBLE PRECISION

c     IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c 2/5 add imp real8
      implicit real*8 (a-h,o-z)
      DIMENSION RVEC(*)
      COMMON/R48ST1/U(97),C,I97,J97
      PARAMETER (MODCNS=1000000000)
      SAVE CD, CM, TWOM24,  ZERO, ONE, NTOT, NTOT2, IJKL
      DATA NTOT,NTOT2,IJKL/-1,0,0/

      DATA LUNOUT/16/
C
      IF (NTOT .GE. 0)  GO TO 50
C
C        Default initialization. User has called RM48 without RM48IN.
      IJKL = 54217137
      NTOT = 0
      NTOT2 = 0
      KALLED = 0
      GO TO 1
C
      ENTRY      RM48IN(IJKLIN, NTOTIN,NTOT2N)
C         Initializing routine for RM48, may be called before
C         generating pseudorandom numbers with RM48.   The input
C         values should be in the ranges:  0<=IJKLIN<=900 OOO OOO
C                                          0<=NTOTIN<=999 999 999
C                                          0<=NTOT2N<<999 999 999!
C To get the standard values in Marsaglia's paper, IJKLIN=54217137
C                                            NTOTIN,NTOT2N=0
      IJKL = IJKLIN
      NTOT = MAX(NTOTIN,0)
      NTOT2= MAX(NTOT2N,0)
      KALLED = 1
C          always come here to initialize
    1 CONTINUE
      IJ = IJKL/30082
      KL = IJKL - 30082*IJ
      I = MOD(IJ/177, 177) + 2
      J = MOD(IJ, 177)     + 2
      K = MOD(KL/169, 178) + 1
      L = MOD(KL, 169)
c     WRITE(LUNOUT,'(A,I10,2X,2I10)')
c    & ' RM48 INITIALIZED:',IJKL,NTOT,NTOT2
CCC      PRINT '(A,4I10)', '   I,J,K,L= ',I,J,K,L
      ONE = 1.D+00
      HALF = 0.5D+00
      ZERO = 0.D+00
      DO 2 II= 1, 97
      S = 0.D+00
      T = HALF
      DO 3 JJ= 1, 48
         M = MOD(MOD(I*J,179)*K, 179)
         I = J
         J = K
         K = M
         L = MOD(53*L+1, 169)
         IF (MOD(L*M,64) .GE. 32)  S = S+T
    3    T = HALF*T
    2 U(II) = S
      TWOM24 = ONE
      DO 4 I24= 1, 24
    4 TWOM24 = HALF*TWOM24
      C  =   362436.D+00*TWOM24
      CD =  7654321.D+00*TWOM24
      CM = 16777213.D+00*TWOM24
      I97 = 97
      J97 = 33
C       Complete initialization by skipping
C            (NTOT2*MODCNS + NTOT) random numbers
      DO 45 LOOP2= 1, NTOT2+1
      NOW = MODCNS
      IF (LOOP2 .EQ. NTOT2+1)  NOW=NTOT
      IF (NOW .GT. 0)  THEN
c     WRITE(LUNOUT,'(A,I15)') ' RM48IN SKIPPING OVER ',NOW
          DO 40 IDUM = 1, NTOT
          UNI = U(I97)-U(J97)
          IF (UNI .LT. ZERO)  UNI=UNI+ONE
          U(I97) = UNI
          I97 = I97-1
          IF (I97 .EQ. 0)  I97=97
          J97 = J97-1
          IF (J97 .EQ. 0)  J97=97
          C = C - CD
          IF (C .LT. ZERO)  C=C+CM
   40     CONTINUE
      ENDIF
   45 CONTINUE
      IF (KALLED .EQ. 1)  RETURN
C
C          Normal entry to generate LENV random numbers
   50 CONTINUE
      DO 100 IVEC= 1, LENV
      UNI = U(I97)-U(J97)
      IF (UNI .LT. ZERO)  UNI=UNI+ONE
      U(I97) = UNI
      I97 = I97-1
      IF (I97 .EQ. 0)  I97=97
      J97 = J97-1
      IF (J97 .EQ. 0)  J97=97
      C = C - CD
      IF (C .LT. ZERO)  C=C+CM
      UNI = UNI-C
      IF (UNI .LT. ZERO) UNI=UNI+ONE
      RVEC(IVEC) = UNI
CC             Replace exact zeros by uniform distr. *2**-24
C         IF (UNI .EQ. 0.)  THEN
C         ZUNI = TWOM24*U(2)
CC             An exact zero here is very unlikely, but let's be safe.
C         IF (ZUNI .EQ. 0.) ZUNI= TWOM24*TWOM24
C         RVEC(IVEC) = ZUNI
C         ENDIF
  100 CONTINUE
      NTOT = NTOT + LENV
         IF (NTOT .GE. MODCNS)  THEN
         NTOT2 = NTOT2 + 1
         NTOT = NTOT - MODCNS
         ENDIF
      RETURN
C           Entry to output current status
      ENTRY RM48UT(IJKLUT,NTOTUT,NTOT2T)
      IJKLUT = IJKL
      NTOTUT = NTOT
      NTOT2T = NTOT2
CNVM
c     WRITE(LUNOUT,'(//A,I10,2X,2I10)')
c    & ' RM48UT OUTPUT:',IJKL,NTOT,NTOT2
CNVM
      RETURN
      END
c
c ***************************************************
c
      subroutine fermfrag(kst)
      common/distb/dsgp(36,999,8)
      common/outlim/limout,lim2
      COMMON /PRODA/prod (7,20)
      common/levop/ldopt,igam
      common/ismrdat2/ididl(106,260),nolevs(106,260)
      common/perp/spinprp
      common/isomer1/xgs(690),xiso(690),x2iso(690),isonum,isoint,iso2num
      common/isomer2/zs(690),as(690),eiso(690),eiso2(690)
      common/isomer3/iziso(690),iaiso(690),csgs(690),csis(690)
     1,csgse(690,100),csise(690,100),csis2e(690,100),rm1(690,100),
     2rm2(690,100)
      common/lin/lspin
      common/lscpl/bx(180),sppj
      common/ijcpl/ajtarg,ajin
      common/tagspin/trgspin
      common/tgspin/levno
      common/spins2/spinu,spinf,prz1
      common/rotor/irotor,iidummy,arote(15,25,181)
      common/history/ide(60),idn(60),idith(60),iddum(60),the(60),adelj
     1(60),spinfh(60),enghis(60)
       common/spins4/eprot,ealpha,eoff(8)
      common/spins/rote(15,24,181),popej(15,24,180,100),fcx(180),ai,ecmp
      common/sping/pfpej(600,180,100),pfpj(600,180)
      common/gate1/spind(5,180),popj(15,24,180),egatl(5),egatu(5),npart(
     15),igate
      common/nuopt/kinem,iem,kb
      common/excid/npartt(2)
      common/tem/philabb(36),phirec(36),phb(36),phin(36)
      common/bound/neut2,iprot,idnuc(100),idm(100),einc(100),
     1thett(100),phit(100),ppp
      common/hyb2/pp(15,24,999)
      common/evapmc/jcount,ipremc,eminn(15,24)
      common/ilcom/nu,ndum,am(8)
      common/rcl/amis
      common/labcon/prx,pry,prz,pcm,rmas,clab(999,8),dslb(36,999,8)
      common/exclbe/idz(4),ida(4),beff(20,20,2),pairmc(18,27)
      common/scat/jz,ja, probxx(2,2),q1(8,1100),fcs,eloss,fcspec
        common/mccom/ncount,ipot,rum(15,24),ppmc(15,24,999),cf
      common/exclus/ gib(8,999),pr(3),ge(2100),zz(3),
     1delta(3),ts(3),pairx(8),asum(10),rate(8,999)
      common/shft/mc,mp,inver,ike,ipch,kq,qval,ap,at,zp,zt,cld,barfac
      common/pl3/jx,nepr,eb(100),rcsp(100),zee,amass,plex
      common/par3/ep,sigml(999),acrs(999)
      common/memo/scrs(999,8)
      common/incr/ed,tened,efinal
      common/sft5/exc(18,27),xmxx
      common/lym1/be(18,27,8),pair(18,27),xmas(18,27),delshl(18,27),amas
     1(18,27)
c     common/bym1/symb(18,27),symbp(18,27),symbs(18,27)
      common/dist/n,ni,ke,irefr,xz(3),pod,pcs,ek,xmix,cfrac
     1,dtau,sign(36),tmix,bji(10),fin,dtup
      common/error/rerror(15,24)
      common/eder/ed2,v,vhalf,vsq
      common/angprm/aprp(500),eav3(1000),
     1eav2(1000),pot2x,thet1,phi1,snt(360),cst(360),snp(720),csp(720)
     2,pot
crc
      common/out1/nout
      common/out2/vrl(200,25),vrlz(200,25),vrlx(200,25),vrly(200,25)
      common/out3/vscale
      common/out4/asm(36),asmm(36)
crc
      common/recoils/rcletza(100,13,24,18),rcleza(100,13,24),
     1rclea(100,36),rcleat(100,36,18),rphir(36)
     2,rclsp(200),rclz(200),rcly(200),rclx(200),escale
     3,rcle(200),rclt(36),rclet(200,36)
      common/rec1/ermxc,thetrc
      common/mmtcom/mmt,mzee,iares,izres,light
      common/fermst/aferm
c
c loop over all fragments numbering 'kst'
       do 1000 nfr=1,kst
       izz=prod(2,nfr)
       iaa=prod(1,nfr)
       if(izz.gt.5)go to 600
        if(izz.eq.5)go to 500
        if(izz.eq.4)go to 400
        if(izz.eq.3)go to 300
        if(izz.eq.2)go to 200
        if(izz.eq.1)go to 100
        if(izz.eq.0)go to 50
   50  continue
c store neutron
       n1=1
       ix=(prod(6,nfr)+ed)/ed
c prod(6,nfr) is the cm kinetic energy of frag nfr
       ith=thet1*11.459155+1.
       if(ith.eq.37)ith=36
       scrs(ix,n1)=scrs(ix,n1)+fcspec
       dsgp(ith,ix,n1)=dsgp(ith,ix,n1)+fcspec
       go to 1000
c do z=1 products
  100  if(iaa.eq.1)then
       n1=2
       else
       n1=4
c assign all H isotopes A gt 1. as deuterons for present
       endif
       ix=(prod(6,nfr)+ed)/ed
       ith=thet1*11.459155+1.
       if(ith.eq.37)ith=36
       scrs(ix,n1)=scrs(ix,n1)+fcspec
       dsgp(ith,ix,n1)=dsgp(ith,ix,n1)+fcspec
       go to 1000
  200  continue
c treat Z=2 fragments
       n1=3
       if(iaa.lt.5)then
       ix=(prod(6,nfr)+ed)/ed
       ith=thet1*11.459155+1.
       if(ith.eq.37)ith=36
       scrs(ix,n1)=scrs(ix,n1)+fcspec
       dsgp(ith,ix,n1)=dsgp(ith,ix,n1)+fcspec
       endif
c assume that He5 gives He4+neut,energies according to mass
c not applying n-alpha momentum cons. to this point
       if(iaa.eq.5)then
       ixa=(0.9+prod(6,nfr)+ed)*0.8/ed
       if(ixa.lt.1)ixa=1
       nn=1
       ixn=(0.2*(0.9+prod(6,nfr)+ed))/ed
       if(ixn.lt.1)ixn=1
       ith=thet1*11.459155+1.
       if(ith.eq.37)ith=36
       scrs(ixa,n1)=scrs(ixa,n1)+fcspec
       dsgp(ith,ixa,n1)=dsgp(ith,ixa,n1)+fcspec
       scrs(ixn,nn)=scrs(ixn,nn)+fcspec
       dsgp(ith,ixn,nn)=dsgp(ith,ixn,nn)+fcspec
       endif
c next do He6,
       if(iaa.gt.5)then
        kan=ja+iaa-izz
        kzn=jz+izz
       if(kan.lt.23.and.kzn.lt.14)then
c      pp(kan,kzn,1)=pp(kan,kzn,1)+fcs
       pp(kzn,kan,1)=pp(kzn,kan,1)+fcs
       endif
       endif
       go to 1000
  300  continue
c  do Li isotopes,beginning with Li5
       if(iaa.eq.5)then
       ixa=(2.0+prod(6,nfr)+ed)*0.80 /ed
       if(ixa.lt.1)ixa=1
       n1=3
       nn=2
c id proton as second particle,alpha as first from Li5 decay
       ixn=(0.200*(2.0+prod(6,nfr)+ed))/ed
       if(ixn.lt.1)ixn=1
       ith=thet1*11.459155+1.
       if(ith.eq.37)ith=36
       scrs(ixa,n1)=scrs(ixa,n1)+fcspec
       dsgp(ith,ixa,n1)=dsgp(ith,ixa,n1)+fcspec
c     proton
       scrs(ixn,nn)=scrs(ixn,nn)+fcspec
       dsgp(ith,ixn,nn)=dsgp(ith,ixn,nn)+fcspec
       endif
c do other Li isotopes
       if(iaa.ne.5)then
        kan=ja+iaa-izz
        kzn=jz+izz
       if(kan.lt.23.and.kzn.lt.14)then
c      pp(kan,kzn,1)=pp(kan,kzn,1)+fcs
       pp(kzn,kan,1)=pp(kzn,kan,1)+fcs
       endif
       endif
       go to 1000
 400   continue
c treat Berylium isotopes(Z=4)
c we will decay Be6 to  alpha plus 2 prots
c we will decay Be8 to 2 alphas
       if(iaa.ne.6.and.iaa.ne.8)then
        kan=ja+iaa-izz
        kzn=jz+izz
       if(kan.lt.23.and.kzn.lt.14)then
c      pp(kan,kzn,1)=pp(kan,kzn,1)+fcs
       pp(kzn,kan,1)=pp(kzn,kan,1)+fcs
       endif
       endif
       if(iaa.eq.8)then
       n1=3
       ix=0.5*((prod(6,nfr)+ed)/ed)
       if(ix.lt.1)ix=1
       ith=thet1*11.459155+1.
       if(ith.eq.37)ith=36
       scrs(ix,n1)=scrs(ix,n1)+fcspec*2.
       dsgp(ith,ix,n1)=dsgp(ith,ix,n1)+fcspec*2.
       endif
       if(iaa.eq.6)then
c assume Be6 gives 2 prots +alpha
       n1=3
       ix=0.6666*(9.27+(prod(6,nfr)+ed)/ed)
       if(ix.lt.1)ix=1
       ith=thet1*11.459155+1.
       if(ith.eq.37)ith=36
       scrs(ix,n1)=scrs(ix,n1)+fcspec
       dsgp(ith,ix,n1)=dsgp(ith,ix,n1)+fcspec
       iy=0.1666*(9.27+(prod(6,nfr)+ed)/ed)
       if(iy.lt.1)iy=1
       nn=2
       scrs(iy,nn)=scrs(iy,nn)+fcspec*2.
       dsgp(ith,iy,nn)=dsgp(ith,iy,nn)+fcspec*2.
       endif
       go to 1000
 500   continue
c here we treat Boron isotopes(Z=5)
c we will decay B9 to a proton plus 2 alphas
       if(iaa.ne.9)then
        kan=ja+iaa-izz
        kzn=jz+izz
       if(kan.lt.23.and.kzn.lt.14)then
c      pp(kan,kzn,1)=pp(kan,kzn,1)+fcs
       pp(kzn,kan,1)=pp(kzn,kan,1)+fcs
       endif
       endif
       if(iaa.eq.9)then
       n1=3
       ix=0.88888*(0.2 +(prod(6,nfr)+ed)/ed)/2.
 
       if(ix.lt.1)ix=1
       ith=thet1*11.459155+1.
       if(ith.eq.37)ith=36
       scrs(ix,n1)=scrs(ix,n1)+fcspec*2.
       dsgp(ith,ix,n1)=dsgp(ith,ix,n1)+fcspec*2.
       iy=0.1111*(0.20+(prod(6,nfr)+ed)/ed)
       if(iy.lt.1)iy=1
       nn=2
       scrs(iy,nn)=scrs(iy,nn)+fcspec
       dsgp(ith,iy,nn)=dsgp(ith,iy,nn)+fcspec
       endif
       go to 1000
 
  600  continue
c catch all other Z here
        kan=ja+iaa-izz
        kzn=jz+izz
       if(kan.lt.23.and.kzn.lt.14)then
       pp(kzn,kan,1)=pp(kzn,kan,1)+fcs
       endif
1000   continue
       return
       end
      subroutine photonj
      common/ijcpl/ajtarg,ajin
      common/tagspin/trgspin
      common/spins2/spinu,spinf,prz1
      common/phot/bspin(3)
      common/fc/fc(3),jlim
c------------------------------------
       if(jlim.eq.1)then
        spinu=1.
         return
          endif
c----------------------------------------
        if(jlim.eq.2)then
         spinu=1.5
          x=rndm(-1.)
          if(x.lt.fc(1))then
           spinu=0.5
            return
             endif
             endif
c--------------------------------------------
         if(jlim.eq.3)then
          spinu=ajtarg
           x=rndm(-1.)
            if(x.lt.fc(1))then
            spinu=abs(ajtarg-1.)
             return
              endif
               if(x.lt.(fc(1)+fc(2)))then
                spinu=ajtarg+1.
                 endif
                  endif
       return
        end

c xxxo
      subroutine photspin
c sub to couple entry channel spin randomly for excl3 photon precompd
      common/ijcpl/ajtarg,ajin
      common/tagspin/trgspin
      common/spins2/spinu,spinf,prz1
      common/phot/bspin(3)
      common/fc/fc(3),jlim
       ajtarg=trgspin
c for photons incident,set initial spin parameters
c use 'jlim' to flag photon coupling options vs tag spin
       if(ajtarg.eq.0.)then
       jlim=1
       bspin(1)=1.
       bspin(2)=0.
       fc(1)=1.0
       fc(2)=0.
       endif
       if(ajtarg.eq.0.5)then
       jlim=2
       bspin(1)=0.5
       bspin(2)=1.5
       fc(1)=2./6.
       fc(2)=4./6.
       endif
c      write(*,*)'ajtarg',ajtarg
       if(ajtarg.gt.0.5)then
       jlim=3
       bspin(1)=abs(ajtarg-1.)
       bspin(2)=ajtarg+1.
       bspin(3)=ajtarg
       den=2.*bspin(1)+2.*bspin(2)+3.+2.*bspin(3)
       fc(1)=(2.*bspin(1)+1.)/den
       fc(2)=(2.*bspin(2)+1.)/den
       fc(3)=(2.*bspin(3)+1.)/den
        endif
      return
       end
c
c---------------------------------------------------------------------------

c---------------------------------------------------------------
      subroutine endvsel
      common/endvparm/iend,lll(6,6,6),ll(6,6)
c this sub assigns a different 'page' to each type of exclusive
c reaction regardless of the sequence of emission, i.e. whether
c it is (x,np) or (x,pn), etc. goes on the same page.
      ll(1,1)=1
      ll(1,2)=2
      ll(2,1)=2
      ll(2,2)=3
      ll(2,3)=4
      ll(3,2)=4
      ll(3,3)=5
      ll(1,3)=6
      ll(3,1)=6
c           
      lll(1,1,1)=1      
      lll(2,2,2)=2
      lll(3,3,3)=3
      lll(1,1,2)=4      
      lll(1,2,1)=4      
      lll(2,1,1)=4      
      lll(1,1,3)=5      
      lll(1,3,1)=5      
      lll(3,1,1)=5      
      lll(1,2,2)=6      
      lll(2,1,2)=6      
      lll(2,2,1)=6      
      lll(2,2,3)=7      
      lll(2,3,2)=7      
      lll(3,2,2)=7      
      lll(1,3,3)=8      
      lll(3,1,3)=8      
      lll(3,3,1)=8      
      lll(3,3,2)=9      
      lll(3,2,3)=9      
      lll(2,3,3)=9      
      lll(1,2,3)=10      
      lll(1,3,2)=10     
      lll(3,1,2)=10     
      lll(3,2,1)=10     
      lll(2,3,1)=10     
      lll(2,1,3)=10     
      return
      end
c------------------------------------------------------------------
      subroutine endf(thet1,en1lab,ith,ix,n1)
      common/history/ide(60),idn(60),idith(60),iddum(60),the(60),adelj
     1(60),spinfh(60),enghis(60)
      common/endvarray/endvcm1(36,999,8),endvcm2(6,36,999,8),
     1endvcm3(10,36,999,8),cmen1(8,999),cmen2(6,8,999),cmen3(10,8,999)
      common/endvparm/iend,lll(6,6,6),ll(6,6)
      common/histindx/ihist
      common/out1/nout
       ihist=ihist+1
       nout=nout+1
       the(ihist)=thet1
       idith(ihist)=ith
       ide(ihist)=ix
       idn(ihist)=n1
       enghis(ihist)=en1lab
      return
      end
c----------------------------------------------------------------
c sub endvsort1 should sort 1,2,3 particle exclusive events for sdcs (only n,p,4He initially)
      subroutine endvsort1
      common/nuopt/kinem,iem,kb
      common/scat/jz,ja, probxx(2,2),q1(8,1100),fcs,eloss,fcspec
      common/histindx/ihist
      common/out1/nout
      common/history/ide(60),idn(60),idith(60),iddum(60),the(60),adelj
     1(60),spinfh(60),enghis(60)
      common/endvarray/endvcm1(36,999,8),endvcm2(6,36,999,8),
     1endvcm3(10,36,999,8),cmen1(8,999),cmen2(6,8,999),cmen3(10,8,999)
      common/endvparm/iend,lll(6,6,6),ll(6,6)
c     nout=iddum(ihist)
        if(nout.eq.3)go to 200
        if(nout.eq.2)go to 100
        ix=ide(1)
        n11=idn(1)
c       enlab=enghis(1)
c         if(kinem.eq.1)then
c       ith=idith(1)
c        endvcm1(ith,ix,n11)=endvcm1(ith,ix,n11)+fcspec
c         endif
         cmen1(n11,ix)=cmen1(n11,ix)+fcspec
c here ok part id up to 8 for 1 particle out only
        return
  100    continue
         n11=idn(1)
         n12=idn(2)
          if(n11.gt.3.or.n12.gt.3)return
c do not store other than n,p,4He
         i=ll(n11,n12)
          do ih=1,2
           ix=ide(ih)
           n11=idn(ih)
c          enlab=enghis(ih)
c         if(kinem.eq.1)then
c          ith=idith(ih)
c         endvcm2(i,ith,ix,n11)=endvcm2(i,ith,ix,n11)+fcspec
c          endif
         cmen2(i,n11,ix)=cmen2(i,n11,ix)+fcspec
            enddo
            return
  200     continue
         n11=idn(1)
         n12=idn(2)
         n13=idn(3)
          if (n11.gt.3.or.n12.gt.3.or.n13.gt.3)return
         i=lll(n11,n12,n13)
          do ihist=1,3
           ix=ide(ihist)
           n11=idn(ihist)
c         if(kinem.eq.1)then
c          ith=idith(ihist)
c          enlab=enghis(ihist)
c         endvcm3(i,ith,ix,n11)=endvcm3(i,ith,ix,n11)+fcspec
c           endif
         cmen3(i,n11,ix)=cmen3(i,n11,ix)+fcspec
            enddo
       return
       end
c sub endvsort should sort 1,2,3 particle exclusive events for ddcs (only n,p,4He initially)
      subroutine endvsort
      common/nuopt/kinem,iem,kb
      common/scat/jz,ja, probxx(2,2),q1(8,1100),fcs,eloss,fcspec
      common/histindx/ihist
      common/out1/nout
      common/history/ide(60),idn(60),idith(60),iddum(60),the(60),adelj
     1(60),spinfh(60),enghis(60)
      common/endvarray/endvcm1(36,999,8),endvcm2(6,36,999,8),
     1endvcm3(10,36,999,8),cmen1(8,999),cmen2(6,8,999),cmen3(10,8,999)
      common/endvparm/iend,lll(6,6,6),ll(6,6)
c     nout=iddum(ihist)
        if(nout.eq.3)go to 200
        if(nout.eq.2)go to 100
        ix=ide(1)
        n11=idn(1)
c       enlab=enghis(1)
          if(kinem.eq.1)then
        ith=idith(1)
         endvcm1(ith,ix,n11)=endvcm1(ith,ix,n11)+fcspec
          endif
         cmen1(n11,ix)=cmen1(n11,ix)+fcspec
c here ok part id up to 8 for 1 particle out only
        return
  100    continue
         n11=idn(1)
         n12=idn(2)
          if(n11.gt.3.or.n12.gt.3)return
c do not store other than n,p,4He
         i=ll(n11,n12)
          do ih=1,2
           ix=ide(ih)
           n11=idn(ih)
c          enlab=enghis(ih)
          if(kinem.eq.1)then
           ith=idith(ih)
          endvcm2(i,ith,ix,n11)=endvcm2(i,ith,ix,n11)+fcspec
           endif
         cmen2(i,n11,ix)=cmen2(i,n11,ix)+fcspec
            enddo
            return
  200     continue
         n11=idn(1)
         n12=idn(2)
         n13=idn(3)
          if (n11.gt.3.or.n12.gt.3.or.n13.gt.3)return
         i=lll(n11,n12,n13)
          do ihist=1,3
           ix=ide(ihist)
           n11=idn(ihist)
          if(kinem.eq.1)then
           ith=idith(ihist)
c          enlab=enghis(ihist)
          endvcm3(i,ith,ix,n11)=endvcm3(i,ith,ix,n11)+fcspec
            endif
         cmen3(i,n11,ix)=cmen3(i,n11,ix)+fcspec
            enddo
       return
       end
c
       subroutine reclvel
      common/outlim/limout,lim2
      common/par3/ep,sigml(999),acrs(999)
      common/out2/vrl(200,25),vrlz(200,25),vrlx(200,25),vrly(200,25)
      common/out3/vscale
      common/pl3/jx,nepr,eb(100),rcsp(100),zee,amass,plex
      return
crc
       vcomp=1.414*sqrt(ep)/amass
       if(limout.eq.0)then
       write(13,98790)vcomp
98790  format(' composite nucleus velocity =',1pe9.2)
       write(13,98789)
       endif
98789  format(' velocity units=sqrt(2M*E)/M,M=amu,E=MeV')
       do 98785 j=1,10
       vsum=0.
       vsumx=0.
       vsumy=0.
       vsumz=0.
       sumr=0.
       sumx=0.
       sumy=0.
       sumz=0.
       if(limout.eq.0)then
       write(13,98786)j
       write(13,98788)
      endif
98788  format('       velocity  rec vel   z comp     y comp    x comp.')
98786  format(' emission number = ',i3)
       do 98782 i=1,200
       vr=((float(i)-0.5)-50.)/vscale
       vx=vr
       vy=vr
       vz=vr
       sumr=sumr+vr*vrl(i,j)
       sumx=sumx+vx*vrlx(i,j)
       sumy=sumy+vy*vrly(i,j)
       sumz=sumz+vz*vrlz(i,j)
       vsum=vsum+vrl(i,j)
       vsumx=vsumx+vrlx(i,j)
       vsumy=vsumy+vrly(i,j)
       vsumz=vsumz+vrlz(i,j)
       if(vrlx(i,j).eq.0..and.vrly(i,j).eq.0.)go to 98782
       if(limout.eq.0)then
       write(13,98787)vr,vrl(i,j),vrlz(i,j),vrly(i,j),vrlx(i,j)
       endif
98787  format(5x,5(1x,1pe9.2))
98782  continue
       if(limout.eq.0)then
       write(13,98784)
       endif
       avrec=0.
       avrecx=0.
       avrecy=0.
       avrecz=0.
       if(vsum.eq.0.)go to 98780
       avrec=sumr/vsum
98780  if(vsumx.eq.0.)go to 98779
       avrecx=sumx/vsumx
98779  if(vsumy.eq.0.)go to 98778
       avrecy=sumy/vsumy
98778  if(vsumz.eq.0.)go to 98781
       avrecz=sumz/vsumz
98781  continue
       stest=avrec+avrecz+avrecy+avrecx
       if(stest.eq.0.)go to 98785
       if(limout.eq.0)then
       write(13,98783)j,avrec,avrecz,avrecy,avrecx
       write(13,98793)j,vsum,vsumz,vsumy,vsumx
       endif
98793  format(' emission no.= ',i5,'sum vel cs = ',1pe9.2,'sumz=',
     1 1pe9.2,'sumy=',1pe9.2,'sumx=',1pe9.2)
98783  format(4x,i3,7x,4(4x,1pe9.2))
98784  format('emission no. avg recoil v   avg z vel    avg y vel    av
     1g x vel')
98785  continue
      return
      end
      subroutine recoilwr(ermxc)
      common/outlim/limout,lim2
c this subroutine is to write recoil distributions
      common/recar/thetar,er
      common/tem/philabb(36),phirec(36),phb(36),phin(36)
      common/pl3/jx,nepr,eb(100),rcsp(100),zee,amass,plex
      common/recoils/rcletza(100,13,24,18),rcleza(100,13,24),
     1rclea(100,36),rcleat(100,36,18),rphir(36)
     2,rclsp(200),rclz(200),rcly(200),rclx(200),escale
     3,rcle(200),rclt(36),rclet(200,36)
      common/ot/amm(26)
      common/vicdat/sumat(18,26)
      common/prnt/th(37)
cc here return to suppress recoil writing for any limout option
      return
      do 400 i=1,18
      do 401 j=1,26
      sumat(i,j)=0.
 401  continue
 400  continue
c
c print phi tables just to check symmetry
c      if(limout.eq.0)then
c     write(14,1007)
c      write(14,1004)
c     write(14,1007)
c     endif
c     do 105 i=1,36
c1004  format(' phi distribution follows')
c      if(limout.eq.0)then
c     write(14,1005)i,rphir(i),philabb(i),phirec(i),phb(i),phin(i)
c     endif
      continue
c print recoils vs. theta
c      if(limout.eq.0)then
c     write(14,3008)
c     endif
c     do 3007 k=1,4
c     n1=9*(k-1)+1
c     n2=n1+9
c     if(n2.gt.36)n2=36
c     do 3009 l=1,9
c3009  th(l)=-2.5+5.*float(n1)
c      if(limout.eq.0)then
c     write(14,3005)(th(i),i=n1,n2)
c     write(14,3006)(rclt(j),j=n1,n2)
c     endif
c3007  continue
c3005  format('thet= ',9(3x,f4.1,2x))
c3006  format(7x,9(1x,e8.1))
c3008  format('recoil angular distribution')
c2000  format(i5,10(1x,f5.0))
c1005   format(i6,5(1x,1pe9.2))
c      if(limout.eq.0)then
c     write(14,1007)
c     write(14,1008)
c     write(14,1007)
c     endif
c1008  format(' Energy distributions vs. mass number')
c     do 135 i=1,9
c 135 amm(i)=amass+1.-float(i)
c      if(limout.eq.0)then
c     write(14,1007)
c     write(14,1009)(amm(i),i=1,9)
c     write(14,1007)
c     endif
c1009  format(' Erec MeV  A = ',9(f4.0,6x))
c     sumt=0.
c     do 111 ie=1,100
c     erc=(float(ie)-1.)/ermxc
c     sum=0.
c     do 112 i=1,9
c112  sum=sum+rclea(ie,i)
c     sumt=sumt+sum
c     if(sum.eq.0.)go to 111
c      if(limout.eq.0)then
c     write(14,1010)erc,(rclea(ie,i),i=1,9)
c     endif
      continue
c      if(limout.eq.0)then
c     write(14,1007)
c     endif
c     do 136 i=10,18
c136  amm(i)=amass+1.-float(i)
c      if(limout.eq.0)then
c     write(14,1009)(amm(i),i=10,18)
1011  format(' Erec MeV   A = ',8(f5.0,5x))
1010  format(10(1x,1pe9.2))
c     write(14,1007)
c     endif
c     do 116 ie=1,100
c     erc=(float(ie)-1.)/ermxc
c     sum=0.
c     do 115 i=10,18
c115  sum=sum+rclea(ie,i)
c     sumt=sumt+sum
c     if(sum.eq.0.)go to 116
c      if(limout.eq.0)then
c     write(14,1010)erc,(rclea(ie,i),i=10,18)
c     endif
      continue
      write(14,1007)
      do 137 i=19,26
 137  amm(i)=amass+1.-float(i)
       if(limout.eq.0)then
      write(14,1007)
      write(14,1011)(amm(i),i=19,26)
      write(14,1007)
      endif
      do 117 ie=1,100
      erc=(float(ie)-1.)/ermxc
      sum=0.
      do 118 i=19,26
 118  sum=sum+rclea(ie,i)
      sumt=sumt+sum
      if(sum.eq.0.)go to 117
       if(limout.eq.0)then
      write(14,1010)erc,(rclea(ie,i),i=19,26)
      endif
 117  continue
       if(limout.eq.0)then
      write(14,1007)
      write(14,2100)sumt
2100  format(' Cross section for all recoils vs A =',1pe9.2)
      write(14,1007)
      write(14,1012)
1012  format('  Recoil energy spectra vs. Z and A of recoils')
      write(14,1007)
      endif
      sumcsaze=0.
      sumeaze=0.
      do 123 ijz=1,10
      ized=1.+zee -float(ijz)
       if(limout.eq.0)then
      write(14,1013)ized
      endif
1013  format(' Spectra for Atomic number = ',i4)
       do 138 i=1,8
 138   amm(i)=amass-float(ijz)+2.-float(i)
       if(limout.eq.0)then
      write(14,1007)
      write(14,1014)((amm(i)),i=1,8)
      write(14,1007)
      endif
1014  format(' Erec MeV  A = ',8(f5.0,5x))
      do 119 ie=1,100
      erc=(float(ie)-1.)/ermxc
      sum=0.
      do 120 ia=1,8
      sum=sum+rcleza(ie,ijz,ia)
      sumcsaze=sumcsaze+rcleza(ie,ijz,ia)
      sumeaze=sumeaze+erc*rcleza(ie,ijz,ia)
  120 continue
      if(sum.eq.0.)go to 119
       if(limout.eq.0)then
      write(14,1015)erc,(rcleza(ie,ijz,ia),ia=1,8)
       endif
1015  format(2x,1pe9.2,3x,1pe9.2,7(1pe9.2,1x))
  119 continue
       do 139 i=9,18
 139      amm(i)=amass-float(ijz)+2.-float(i)
       if(limout.eq.0)then
      write(14,1007)
      write(14,1014)(  amm(i),i=9,16)
      write(14,1007)
       endif
      do 121 ie=1,100
      erc=(float(ie)-1.)/ermxc
      sum=0.
      do 122 ia=9,16
      sum=sum+rcleza(ie,ijz,ia)
      sumcsaze=sumcsaze+rcleza(ie,ijz,ia)
      sumeaze=sumeaze+erc*rcleza(ie,ijz,ia)
  122 continue
      if(sum.eq.0.)go to 121
       if(limout.eq.0)then
      write(14,1015)erc,(rcleza(ie,ijz,ia),ia=9,16)
      endif
  121 continue
  123 continue
      avgen=sumeaze/sumcsaze
       if(limout.eq.0)then
      write(14,1007)
      write(14,2028)avgen,sumcsaze
      write(14,1007)
      write(14,1021)
      write(14,1007)
      write(14,1025)
      write(14,1007)
      endif
2028  format(' Average recoil energy =',1pe9.2,
     1' Summed recoil cross section = ',1pe9.2)
1025  format(' E recoil MeV  E total MeV E projected E x(transverse) E y
     1 ')
1021  format('  Recoil energy spectra:total,projected(Ez),transverse(Ex,
     1Ey)')
      sumr=0.
      sumrt=0.
      sumz=0.
      sumzt=0.
      sumx=0.
      sumxt=0.
      sumy=0.
      sumyt=0.
      do 140 ie=1,200
      en=(float(ie)-ed/2.)/escale
      sumr=sumr+rclsp(ie)
      sumz=sumz+rclz(ie)
      sumx=sumx+rclx(ie)
      sumy=sumy +rcly(ie)
      sumrt=sumrt+en*rclsp(ie)
      sumzt=sumzt+en*rclz(ie)
      sumxt=sumxt+en*rclx(ie)
      sumyt=sumyt+en*rcly(ie)
      temp=rclsp(ie)+rclz(ie)+rclx(ie)+rcly(ie)
      if(temp.eq.0.)go to 140
       if(limout.eq.0)then
      write(14,1022)en,rclsp(ie),rclz(ie),rclx(ie),rcly(ie)
      endif
 140  continue
      write(14,1007)
1022  format(2x,1pe9.2,4x,4(1x,1pe9.2,1x))
      recav=sumrt/sumr
      recavz=sumzt/sumz
      recavx=sumxt/sumx
      recavy=sumyt/sumy
      tkerma=recav*sumr/1000.
       if(limout.eq.0)then
      write(14,1007)
      write(14,1023)recav,recavz,recavx,recavy
      write(14,1007)
      write(14,1024)sumr,sumz,sumx,sumy
      write(14,1007)
      write(14,2024)sumr,tkerma
 2024 format(' Total cross section =',1pe9.2,
     1' Total recoil Kerma (MeV*barns) =',1pe9.2)
 1023 format(' Ave E rec =',1pe9.2,'Ave Z proj =',1pe9.2,'Ave X proj ='
     1,1pe9.2,' Ave Y proj = ',1pe9.2)
 1024 format(' sum of events recoils,total,z,x,y =',4(1x,1pe9.2))
c  write table of recoil energy vs angle in 10 deg angular bins
1007  format('   ')
      write(14,1007)
      write(14,1001)
      write(14,1007)
1001  format('  table of recoil energy distributions vs angle follows')
      write(14,1007)
      endif
      do 104 ia=1,26
      suma=0.
 
      do 106 ie=1,100
      do 107 ith=1,9
      sumat(ith,ia)=sumat(ith,ia)+rcleat(ie,ia,ith)
      suma=suma+rcleat(ie,ia,ith)
 107  continue
 106  continue
      sumb=0.
 
      do 109 ie=1,100
      do 108 ith=10,18
      sumb=sumb+rcleat(ie,ia,ith)
      sumat(ith,ia)=sumat(ith,ia)+rcleat(ie,ia,ith)
 108  continue
 109  continue
      if(suma.eq.0..and.sumb.eq.0.)go to 104
      amas=amass+1.-float(ia)
       if(limout.eq.0)then
      write(14,1007)
      write(14,1002)amas
      write(14,1007)
      endif
      if(suma.eq.0.)go to 110
1002  format(' recoil E vs theta for mass no. = ',f6.0)
       if(limout.eq.0)then
      write(14,1007)
      write(14,1003)
      write(14,1007)
1003  format(' Erec MeV    0-10  deg 10-20 deg 20-30 deg 30-40 deg ',
     1'40-50 deg 50-60 deg 60-70 deg 70-80 deg 80-90 deg ')
1000  format(1x,1pe9.2,1x,9(1x,1pe9.2))
      endif
      do 100 ie=1,100
      erc=(float(ie)-1.)/ermxc
      sum=0.
      do 101 i=1,9
  101 sum=sum+rcleat(ie,ia,i)
      if(sum.eq.0.)go to 100
       if(limout.eq.0)then
      write(14,1000)erc,(rcleat(ie,ia,ii),ii=1,9)
      endif
  100 continue
  110 continue
      if(sumb.eq.0.)go to 104
       if(limout.eq.0)then
      write(14,1007)
      write(14,1006)
      write(14,1007)
      endif
1006  format(' Erec MeV    90-100deg 100-110dg 110-120dg 120-130dg ',
     1'130-140dg 140-150dg 150-160dg 160-170dg 170-180dg')
      do 102 ie=1,100
      erc=(float(ie)-1.)/ermxc
      sum=0.
      do 103 i=10,18
  103 sum=sum+rcleat(ie,ia,i)
      if(sum.eq.0.)go to 102
       if(limout.eq.0)then
      write(14,1000)erc,(rcleat(ie,ia,i),i=10,18)
      endif
  102 continue
  104 continue
3001  format(' Arec Mass   0-10  deg 10-20 deg 20-30 deg 30-40 deg ',
     1'40-50 deg 50-60 deg 60-70 deg 70-80 deg 80-90 deg ')
3002  format(' Arec Mass   90-100deg 100-110dg 110-120dg 120-130dg ',
     1'130-140dg 140-150dg 150-160dg 160-170dg 170-180dg')
       if(limout.eq.0)then
      write(14,1007)
      write(14,3001)
      write(14,1007)
      endif
      do 299 ia=1,26
       amas=amass+1.-float(ia)
       if(limout.eq.0)then
      write(14,3000)amas,(sumat(ith,ia),ith=1,9)
      endif
 299     continue
       if(limout.eq.0)then
      write(14,1007)
      write(14,3002)
      write(14,1007)
      endif
       do 300 ia=1,26
       amas=amass+1.  -float(ia)
       if(limout.eq.0)then
      write(14,3000)amas,(sumat(ith,ia),ith=10,18)
      endif
  300 continue
3000  format(f5.0,5x,9(1x,1pe9.2))
      return
       end
c add subroutine to calc tables for angular distributions
      subroutine trigstuf
      common/incr/ed,tened,efinal
      common/outlim/limout,lim2
      common/mbcfudg/cfg(999)
      common/eder/ed2,v,vhalf,vsq
c     common/hybsub/iav,idum,rs(2,1100),c,pot,r1
      common/angprm/aprp(500),eav3(1000),
     1eav2(1000),pot2x,thet1,phi1,snt(360),cst(360),snp(720),csp(720)
     2,pot
       do 2 jth=1,360
       arg=(3.14159/360.)*(float(jth)-1.)+3.14159/720.
       snt(jth)=sin(arg)
       cst(jth)=cos(arg)
    2  continue
       do 3 jth=1,720
       arg=(3.14159/360.)*(float(jth)-1.)+3.14159/720.
       snp(jth)=sin(arg)
       csp(jth)=cos(arg)
    3  continue
c
c next store tables which will be used to calculate 'a' parameter of ang dist kernel
c
       ipot=pot
c next calc eav for 2/3 exciton residues where v=pot=efermi
       do 22003 ien=1,ipot-1
c the next 2 lines are true for an equidistant model
       eav3(ien)=pot
       eav2(ien)=pot-float(ien)
cmbc replace eav3 with a Fermi gas result, based on sqrt(e) dependence:
       t1mc=(2./5.)*(pot**2.5 - ((pot-float(ien))**2.5))
       t2mc=(1./3.)*float(ien)*( pot**1.5 - ((pot-float(ien))**1.5) )
       t3mc=(2./3.)*( pot**1.5 -  ((pot-float(ien))**1.5) )
       eav3(ien)=(t1mc+t2mc)/t3mc
cmbc Note that the FG result for eav2 = same as equidistant, so leave unchanged.
       if(eav2(ien).lt.1.)eav2(ien)=0.5
22003  continue
       do 22004 ien=ipot,999
       eav3(ien)=(float(ien)+pot)/2.
       eav2(ien)=0.500
cmbc this corresponds to regions where the bottom of the well is felt.
c Instead, use a Fermi gas result:
       eav3(ien)=(3./5.)*pot + (float(ien)/2.)
22004  continue
c now we will modify the eav values such that aa(e-incid)*aa(e-exit)/eav(E)=a
c i.e. multiply by 2*n/3,where n is the final exciton no.
       do 22005 ien=1,999
       eav2(ien)=eav2(ien)*0.666666
       eav3(ien)=eav3(ien)*1.333333
       cfg(ien)=sqrt(float(ien))/9.3
       if(ien.gt.85)cfg(ien)=1.
22005 continue
c the first theta,t1,for scattering will be given by
c   theta1=arc cos((1./a)*ln(e**(-a)-ranf*(e**a-e**(-a))
c
c  we will next store aparm(1000),and aprp(1000),as(e**a-e**(-a)),and e**(-a),on index mod 10
c with the above definitions,and iena=10.*a,the value of theta is given by
c theta=arc cos(((1./a)ln(aprp(iena)-ranf*aparm(iena))
co this will be in radians,so divide by 0.0872665 after adding 0.087266462,
c  and the index number for the trig function will be given,such that an
c  index no. 1 means an angle falling 0 to 5 deg,etc.
       return
       end
c
      subroutine pickang2(einv,eoutv,arg,th2,elab)
      common/outlim/limout,lim2
      common/pmom/bio
      common/mbcfudg/cfg(999)
      common/angprm/aprp(500),eav3(1000),
     1eav2(1000),pot2x,th,ph,snt(360),cst(360),snp(720),csp(720)
     2,pot
       eres=einv-eoutv
c eres is energy of p+h remaining after rescattered first exciton has the energy en1
c selected -where ein  is the energy of the first scattered exciton
       ieres=eres+1.
       if(ieres.lt.1)ieres=1
c calculate the chadwick 'a' parameter for second scattered particle
c      apr=2.*sqrt(pio)/eav3(ieres)
       ieout=elab
       if(ieout.lt.1)ieout=1
       apr=cfg(ieout)*2.*sqrt(bio)/eav2(ieres)
c      apr=cfg(ieout)*2.*sqrt(einv*eoutv)/eav2(ieres)
       if(apr.eq.0.)then
         th2=0.
         arg=1.
         go to 78
        endif
       xx=rndm(-1.)
       if(xx.eq.1.)xx=0.99999
       if(apr.lt.30.)go to 76
       arg=(apr+log(1.-xx))/apr
       go to 77
  76   xxx=exp(apr)
       arg=log((1.-xx)*xxx+xx/xxx)/apr
   77    if(abs(arg).gt.1.)arg=arg/abs(arg)
       th2=acos(arg)
  78   continue
       return
       end
c
c write pickang for nucleons of incident HI
      subroutine pickanghi(einv,eoutv,arg,thx,elab,app)
      common/himbc1/threc,phrec
      common/himbc/th2,ph2,przz,px,py,pz
      common/outlim/limout,lim2
      common/pmom/bio
      common/mbcfudg/cfg(999)
      common/angprm/aprp(500),eav3(1000),
     1eav2(1000),pot2x,th,ph,snt(360),cst(360),snp(720),csp(720)
     2,pot
       eres=einv-eoutv
c eres is energy of p+h remaining after rescattered first exciton has the energy en1
c selected -where ein  is the energy of the first scattered exciton
       ieres=eres+1.
       if(ieres.lt.1)ieres=1
c calculate the chadwick 'a' parameter for second scattered particle
c      apr=2.*sqrt(einv*eoutv)/eav3(ieres)
       ieout=elab
       if(ieout.lt.1)ieout=1
c      apr=cfg(ieout)*2.*sqrt(einv*eoutv)/eav3(ieres)
c      apr=cfg(ieout)*2.*sqrt(bio)/eav3(ieres)
c      apr=cfg(ieout)*3.*(bio)/eres
c here check where bio is defined
c      if(apr.eq.0.)then
c        th2=0.
c        arg=1.
c        go to 78
c       endif
c below we use new def of bio from rohi
       apr=cfg(ieout)*bio
       xx=rndm(-1.)
       if(xx.eq.1.)xx=0.99999
       if(apr.lt.30.)go to 76
       arg=(apr+log(1.-xx))/apr
       go to 77
  76   xxx=exp(apr)
       arg=log((1.-xx)*xxx+xx/xxx)/apr
 77    if(abs(arg).gt.1.)arg=arg/abs(arg)
       th2=acos(arg)
       ph2=rndm(-1.)*6.283185
       continue
c here we mbc select the angle theta of the outward pointed nucleon
c from the cluster of previous mass app.
c should this be defined as th,ph or th2,ph2?
c now we need to define the angles of the recoil, then rotate back into Z 
c direction, then if necessary, rotate the th2,ph2 into Z direction.
       if(limout.eq.0)then
       if(einv.gt.198)write(19,100)einv,eoutv,eres,th2
      endif
 100   format(' from pi0-einv,eoutv,eres,th2 ',4f10.3)
c      do 900 ien=1,100,5
c900      write(19,5555)cfg(ien),eav3(ien),ien
c5555  format(' cfg,eav3,ien ',2f10.3,i5)
       return
       end
c
      subroutine pickang(einv,eoutv,arg,th2,elab)
      common/outlim/limout,lim2
      common/pmom/bio
      common/mbcfudg/cfg(999)
      common/angprm/aprp(500),eav3(1000),
     1eav2(1000),pot2x,th,ph,snt(360),cst(360),snp(720),csp(720)
     2,pot
       eres=einv-eoutv
c eres is energy of p+h remaining after rescattered first exciton has the energy en1
c selected -where ein  is the energy of the first scattered exciton
       ieres=eres+1.
       if(ieres.lt.1)ieres=1
c calculate the chadwick 'a' parameter for second scattered particle
c      apr=2.*sqrt(einv*eoutv)/eav3(ieres)
       ieout=elab
       if(ieout.lt.1)ieout=1
c      apr=cfg(ieout)*2.*sqrt(einv*eoutv)/eav3(ieres)
       apr=cfg(ieout)*2.*sqrt(bio)/eav3(ieres)
c      if(apr.eq.0.)then
c        th2=0.
c        arg=1.
c        go to 78
c       endif
       xx=rndm(-1.)
       if(xx.eq.1.)xx=0.99999
       if(apr.lt.30.)go to 76
       arg=(apr+log(1.-xx))/apr
       go to 77
  76   xxx=exp(apr)
       arg=log((1.-xx)*xxx+xx/xxx)/apr
 77    if(abs(arg).gt.1.)arg=arg/abs(arg)
       th2=acos(arg)
       continue
       if(limout.eq.0)then
       if(einv.gt.198)write(19,100)einv,eoutv,eres,th2
      endif
 100   format(' from pi0-einv,eoutv,eres,th2 ',4f10.3)
c      do 900 ien=1,100,5
c900      write(19,5555)cfg(ien),eav3(ien),ien
c5555  format(' cfg,eav3,ien ',2f10.3,i5)
       return
       end
c
c add subroutine to do angle rotation
       subroutine rotate
c define function for x,y,z after scattering
c      gq(ith1,ith2,iph2)=cst(ith1)*cst(ith2)-snt(ith1)*snt(ith2)*csp(iph2)
c    1h2)
c
c      aq(ith1,ith2,iph1,iph2)=cst(ith1)*snt(ith2)*csp(iph1)*csp(iph2)
c    1 -snt(ith2)*snp(iph1)*snp(iph2)+csp(iph1)*cst(ith2)*snt(ith1)
c
c      bq(ith1,ith2,iph1,iph2)=csp(iph2)*snp(iph1)*cst(ith1)*snt(ith2)
c    1 +csp(iph1)*snp(iph2)*snt(ith2)+snt(ith1)*snp(iph1)*cst(ith2)
      common/angprm/aprp(500),eav3(1000),
     1eav2(1000),pot2x,th,ph,snt(360),cst(360),snp(720),csp(720)
     2,pot
       common/indices/ith1,ith2,iph1,iph2
       common/rotout/zed,aqn,bqn
       zed=cst(ith1)*cst(ith2)-snt(ith1)*snt(ith2)*csp(iph2)
 
c
       aqn=cst(ith1)*snt(ith2)*csp(iph1)*csp(iph2)
     1 -snt(ith2)*snp(iph1)*snp(iph2)+csp(iph1)*cst(ith2)*snt(ith1)
c
       bqn=csp(iph2)*snp(iph1)*cst(ith1)*snt(ith2)
     1 +csp(iph1)*snp(iph2)*snt(ith2)+snt(ith1)*snp(iph1)*cst(ith2)
       return
       end
c
c
       subroutine rotang3(einc,eout,thnu1,phnu1,elab)
      common/pmom/bio
      common/iihist/ii
      common/shft/mc,mp,inver,ike,ipch,kq,qval,ap,at,zp,zt,cld,barfac
      common/par3/ep,sigml(999),acrs(999)
       common/indices/ith1,ith2,iph1,iph2
       common/rotout/zed,aqn,bqn
      common/phitst/phian(36),phian2(36),phian3(36),phian4(36),pcon
     1,phian5(36)
c
c
      common/mbcfudg/cfg(999)
      common/angprm/aprp(500),eav3(1000),
     1eav2(1000),pot2x,th,ph,snt(360),cst(360),snp(720),csp(720)
     2,pot
c
       einv=einc+pot
       eoutv=eout+pot
       bio=einv*eoutv
c substitute photon momenta for first scattering for photon induced rx
c
       call pickang(einv,eoutv,arg,th2,elab)
c      call pickanghi(einv,eoutv,arg,th2,elab,app)
c old hi version has app in call to hi pickang
c but app may not be defined in this routine-
c
       ph2=rndm(-1.)*6.283185
       ith2=1.+th2/0.0087267
       ith1=1.+thnu1/0.0087267
       iph2=1.+ph2/0.0087267
       iph1=1.+phnu1/0.0087267
c now calculate the new theta and phi angles related to original axis system
       if(ith1.lt.1)then
       ith1=1
       endif
       if(ith2.lt.1)then
       ith2=1
       endif
       if(ith1.lt.1)then
       ith1=1
       endif
       if(ith2.lt.1)then
       ith2=1
       endif
       call rotate
       if(abs(zed).gt.1.)zed=zed/abs(zed)
       thnu1=acos(zed)
       phnu1=atan2(bqn,aqn)
        if(phnu1.lt.0.)phnu1=6.28318+phnu1
      continue
      return
      end
c
       subroutine rohi(pout,einc,eout,tph,pph,elab,app)
      common/himbc2/clthet(100),clphi(100)
      common/himbc1/threc,phrec
      common/himbc/th2,ph2,przz,px,py,pz
      common/hik/yph,xph,zph
      common/pmom/bio
      common/shft/mc,mp,inver,ike,ipch,kq,qval,ap,at,zp,zt,cld,barfac
      common/iihist/ii
      common/par3/ep,sigml(999),acrs(999)
       common/indices/ith1,ith2,iph1,iph2
       common/rotout/zed,aqn,bqn
      common/phitst/phian(36),phian2(36),phian3(36),phian4(36),pcon
     1,phian5(36)
      common/mbcfudg/cfg(999)
      common/angprm/aprp(500),eav3(1000),
     1eav2(1000),pot2x,thet1,phi1,snt(360),cst(360),snp(720),csp(720)
     2,pot
c
c app is the emitting projectile mass no.
       pcm=sqrt(2.*eout)
       einv=(einc+eout)+pot*app
c einv= total excitation from well bottom for app excitons
c the value passed as einc is the total excitation decremented by eout this event
c einv= Ecm+Q+A*V, where V is Efermi at top of nucleon fill
c and eoutv is e-nucleon out+V
       eoutv=eout+pot
       if(app.eq.1.)then
        th2=rndm(-1.)*3.14159
        ph2=rndm(-1.)*6.283184
       go to 3000
       endif
c define 'bio as the 'an' parameter of mbc,which for a HI will be
c an=3*sqrt((einv*A*eoutv))/((A-1)*m*rho*eav)
c above in hi sub
c the above expression is still  missing the rho in denom.
c define eres as the final proj fragment excitation measured from well bottom
c in the definition of 'bio' the mbc a-sub-n parm,the eav is replaced by eres
c due to cancellation of the 'n' or app-1. term in the denominator
       eres=einc+(app-1.)*pot
       bio=3.*sqrt(app*einv*eoutv)/(eres)
       call pickanghi(einv,eoutv,arg,thq,elab,app)
c pickang sends back the nucleon angle according to mbc
c as th2; also ph2
c600   ph2=rndm(-1.)*6.283185
c     arg=cos(th2)
c set arrays of unrotated theta,phi for xth exciton
c now conserve momentum in cluster system via px,y,z?
c initially we will wait-
c     return
c now do math to get recoil directions
 3000   continue
        pkin=sqrt(einv*2.*app)
        pkout=sqrt(eoutv*2.)
        pkrecsq=pkin*pkin+pkout*pkout-2.*pkin*pkout*cos(th2)
        threc=asin((pkout/sqrt(pkrecsq))*sin(th2))
        if(ph2.le.3.14159)then
         phrec=ph2+3.14159
         else
         phrec=ph2-3.14159
         endif
         if(phrec.lt.0.)phrec=phrec+6.28318
       write(*,*)'pkout,pkrecsq,th2,einv,eoutv',pkout,pkrecsq,th2,einv,
     1eoutv
         write(*,*)'app,threc,phrec before rotating',app,threc,phrec
c these recoil directions are un-rotated to beam direction, as are
c the ejectile directions. Unless this is the first emission, the rotat-
c ions should be performed, and angles wrt beam direction returned for use
c      if(pout.gt.0.)go to 1719
c     tph=0.
c
c     pph=0.
c     go to 1720
c1719     ptem=sqrt(2.*eoutv)
c     yph=ptem*sin(th2)
c     xph=sqrt(2.*einv*app)-ptem*arg
c     arg=cos(th2)
c     tph=atan2(yph,xph)
c  above to calc  ph angle by mom. cons in 2p1h collision
c1720      pph=ph2+3.141592654
c      if(pph.gt.6.283185307)pph=pph-6.283185
c now have conserved momenta to calc tph,pph the angles of
c the particle-hole residue
c check if above subtract is pi or 2pi
c       thet1=threc
c      tph=threc
c      pph=phrec
c       phi1=phrec
       ith2=1.+th2/0.0087267
       ith1=1.+thet1/0.0087267
       iph2=1.+ph2/0.0087267
       iph1=1.+phi1/0.0087267
c stopped here tempZZ
c now calculate the new theta and phi angles related to original axis system
c for the scattered particle of the 2p1h,next do it for the ph pair
       if(ith1.lt.1)then
       ith1=1
       endif
       if(ith2.lt.1)then
       ith2=1
       endif
       write(*,*)' before rotate,th2,ph2,thet1,phi1,app',th2,ph2,thet1,
     1phi1,app
       call rotate
       if(abs(zed).gt.1.)zed=zed/abs(zed)
       thet1=acos(zed)
       i=(ap)-(app)+1.
       write(*,*)'i',i
       clthet(i)=thet1
       phi1=atan2(bqn,aqn)
       if(phi1.lt.0.)phi1=6.28318+phi1
       clphi(i)=phi1


c here, clthet and clphi should be emitted  particle angle rotated directions
       write(*,*)'  after rotate,th2,ph2,thet1,phi1,app',clthet(i),clph
     1i(i),thet1,phi1,app
c here add sum over emitted momenta to get final exciton by conservation
      ith=1.+thet1/0.0087267
      iph=1.+phi1/0.0087267
c thet,phi are the emission angles wrt z beam axis in recoil frame
c9999  format(' phi=',f10.3,'iph=',i5,'csp=',f10.3,'snp=',f10.3)
      continue
c      write(*,*)'pcm in recoil',pcm
      px=px+pcm*snt(ith)*csp(iph)
      py=py+pcm*snt(ith)*snp(iph)
      pz=pz+pcm*cst(ith)
      write(*,*)'px,py,pz in rohi',px,py,pz
      continue
c now thet1 and phi1 represent the scattered nucleon angles
c rotated back to the original beam direction for those cases
c where it started with non-zero initial values.returned in common
c add rot for p+h
       write(*,*)'  before rotate,tph,pph,,app',tph,pph,app
      ith2=tph/0.0087267+1.
      iph2=pph/0.0087267+1.
       if(ith1.lt.1)then
       ith1=1
       endif
       if(ith2.lt.1)then
       ith2=1
       endif
      call rotate
       if(abs(zed).gt.1.)zed=zed/abs(zed)
       tph=acos(zed)
       pph=atan2(bqn,aqn)
        if(pph.lt.0.)pph=6.28318+pph
       write(*,*)'  after rotate,tph,pph,,app',tph,pph,app
       continue
c now tph and pph represent the scattered residual cluster angles
c rotated back to the original beam direction for those cases
c where it started with non-zero initial values.returned in call
c argument
      return
      end
c
c
       subroutine rotangl(pout,einc,eout,tph,pph,elab)
      common/hik/yph,xph,zph
      common/pmom/bio
      common/shft/mc,mp,inver,ike,ipch,kq,qval,ap,at,zp,zt,cld,barfac
      common/iihist/ii
      common/par3/ep,sigml(999),acrs(999)
       common/indices/ith1,ith2,iph1,iph2
       common/rotout/zed,aqn,bqn
      common/phitst/phian(36),phian2(36),phian3(36),phian4(36),pcon
     1,phian5(36)
      common/mbcfudg/cfg(999)
      common/angprm/aprp(500),eav3(1000),
     1eav2(1000),pot2x,thet1,phi1,snt(360),cst(360),snp(720),csp(720)
     2,pot
c
       einv=einc+pot
       eoutv=eout+pot
       bio=einv*eoutv
c      bio=sqrt(app*einv*eoutv)
c above in hi sub
c change for photons 9/02
       if(ap.eq.0..and.ii.eq.1)then
c      bio=0.
       bio=2.*0.00053*einc*einc*eoutv
       einvg=bio/eoutv
       einv=einc
       endif
c substitute photon momenta for first scattering for photon induced rx
       call pickang(einv,eoutv,arg,th2,elab)
c      call pickanghi(einv,eoutv,arg,th2,elab,app)
c pickang sends back the nucleon angle according to mbc
c as th2; select phi as ph2
       ph2=rndm(-1.)*6.283185
      arg=cos(th2)
      if(ap.eq.0..and.ii.eq.1)then
      pout=abs(2.*(einvg+eoutv)-4.*arg*sqrt(einvg*eoutv))
      else
      pout=abs(2.*(einv+eoutv)-4.*arg*sqrt(einv*eoutv))
      endif
       if(pout.gt.0.)go to 1719
      tph=0.
 
      pph=0.
      go to 1720
 1719     ptem=sqrt(2.*eoutv)
      yph=ptem*sin(th2)
      xph=sqrt(2.*einv)-ptem*arg
c     arg=cos(th2)
      tph=atan2(yph,xph)
c  above to calc  ph angle by mom. cons in 2p1h collision
 1720      pph=ph2+3.141592654
       if(pph.gt.6.283185307)pph=pph-6.283185
c now have conserved momenta to calc tph,pph the angles of
c the particle-hole residue
c check if above subtract is pi or 2pi
c     write(*,*)'exl1 thet1,phi1,tph,pph before rot',thet1,phi1,tph,pph
c     write(*,*)'th2,ph2,elab',th2,ph2,elab
       ith2=1.+th2/0.0087267
       ith1=1.+thet1/0.0087267
       iph2=1.+ph2/0.0087267
       iph1=1.+phi1/0.0087267
c now calculate the new theta and phi angles related to original axis system
c for the scattered particle of the 2p1h,next do it for the ph pair
       if(ith1.lt.1)then
       ith1=1
       endif
       if(ith2.lt.1)then
       ith2=1
       endif
       call rotate
       if(abs(zed).gt.1.)zed=zed/abs(zed)
       thet1=acos(zed)
       phi1=atan2(bqn,aqn)
       if(phi1.lt.0.)phi1=6.28318+phi1
c     write(*,*)'exl1 thet1,phi1,tph,pph after rot',thet1,phi1,tph,pph
c now thet1 and phi1 represent the scattered nuc3yyleon angles
c rotated back to the original beam direction for those cases
c where it started with non-zero initial values.returned in common
c add rot for p+h
      ith2=tph/0.0087267+1.
      iph2=pph/0.0087267+1.
       if(ith1.lt.1)then
       ith1=1
       endif
       if(ith2.lt.1)then
       ith2=1
       endif
c     write(*,*)' do ph part'
c     write(*,*)'exl1 thet1,phi1,tph,pph before rot',thet1,phi1,tph,pph
      call rotate
       if(abs(zed).gt.1.)zed=zed/abs(zed)
       tph=acos(zed)
       pph=atan2(bqn,aqn)
        if(pph.lt.0.)pph=6.28318+pph
c     write(*,*)' do ph part'
c     write(*,*)'exl1 thet1,phi1,tph,pph after rot',thet1,phi1,tph,pph
c now tph and pph represent the scattered nucleon angles
c rotated back to the original beam direction for those cases
c where it started with non-zero initial values.returned in call
c argument
      return
      end
c
c
c add subroutine to do angle rotation
c this is for the 1p1h config giving 1p at angle th3,ph3 wrt 1p1h directions
c tph1,pph1 from calling routine,entered here as thet1,phi1,returned here the
c final angles in original coordinate system as thnu1,phnu1
c
       subroutine rotang2(pin,einc,eout,thol,phol,elab)
      common/labcon/prx,pry,prz,pcm,rmas,clab(999,8),dslb(36,999,8)
      common/pmom/bio
      common/shft/mc,mp,inver,ike,ipch,kq,qval,ap,at,zp,zt,cld,barfac
      common/iihist/ii
      common/par3/ep,sigml(999),acrs(999)
       common/indices/ith1,ith2,iph1,iph2
       common/rotout/zed,aqn,bqn
      common/phitst/phian(36),phian2(36),phian3(36),phian4(36),pcon
     1,phian5(36)
c
      common/mbcfudg/cfg(999)
      common/angprm/aprp(500),eav3(1000),
     1eav2(1000),pot2x,thet1,phi1,snt(360),cst(360),snp(720),csp(720)
     2,pot
       einv=einc+pot
       eoutv=eout+pot
       bio=einv*eoutv
c
       call pickang2(einv,eoutv,arg,th2,elab)
c
           ph2=rndm(-1.)*6.283185
       ith2=1.+th2/0.0087267
       ith1=1.+thet1/0.0087267
       iph2=1.+ph2/0.0087267
       iph1=1.+phi1/0.0087267
      phol=ph2+3.141592654
      ptem=sqrt(2.*eoutv)
      xmhol=sqrt(2.*einv)-arg*ptem
      ymhol=sqrt(2.*eoutv)*sin(th2)
      thol=atan2(ymhol,xmhol)
       if(phol.gt.6.283185307)phol=phol-6.283185307
       if(ith1.lt.1)then
       ith1=1
       endif
       if(ith2.lt.1)then
       ith2=1
       endif
       call rotate
       if(abs(zed).gt.1.)zed=zed/abs(zed)
       thet1=acos(zed)
       phi1=atan2(bqn,aqn)
c the above are rotated particle directionns;must rot hole next
       if(phi1.lt.0.)phi1=6.28318+phi1
       continue
 
c
c now calculate the new theta and phi angles related to original axis system
       ith2=thol/0.0087267+1.
       iph2=phol/0.0087267+1.
       if(ith1.lt.1)then
       ith1=1
       endif
       if(ith2.lt.1)then
       ith2=1
       endif
       call rotate
       if(abs(zed).gt.1.)zed=zed/abs(zed)
        thol=acos(zed)
         phol=atan2(bqn,aqn)
        if(phol.lt.0.)phol=6.28318+phol
      return
      end
c
      subroutine holscat1(nhol,ehol,thethol,phihol)
      common/endvparm/iend,lll(6,6,6),ll(6,6)
      common/lout/spincon
      common/fermst/aferm
      common/lab3/sig(8,999),slav(8,999)
      common/out1/nout
      common/histindx/ihist
      common/history/ide(60),idn(60),idith(60),iddum(60),the(60),adelj
     1(60),spinfh(60),enghis(60)
      common/perp/spinprp
       common/spins4/eprot,ealpha,eoff(8)
      common/spins2/spinu,spinf,prz1
      common/lin/lspin
      common/ijcpl/ajtarg,ajin
      common/lscpl/bx(180),sppj
      common/tagspin/trgspin
      common/tgspin/levno
      common/nuopt/kinem,iem,kb
      common/labcon/prx,pry,prz,pcm,rmas,clab(999,8),dslb(36,999,8)
      common/ilcom/nu,ndum,am(8)
      common/distb/dsgp(36,999,8)
      common/edep/pxx(2,1100)
      common/exclbe/idz(4),ida(4),beff(20,20,2),pairmc(18,27)
      common/scat/jz,ja, probxx(2,2),q1(8,1100),fcs,eloss,fcspec
      common/memo/scrs(999,8)
      common/sft5/exc(18,27),xmxx
      common/lym1/be(18,27,8),pair(18,27),xmas(18,27),delshl(18,27),amas
     1(18,27)
c     common/bym1/symb(18,27),symbp(18,27),symbs(18,27)
      common/pl3/jx,nepr,eb(100),rcsp(100),zee,amass,plex
      common/incr/ed,tened,efinal
      common/eder/ed2,v,vhalf,vsq
      common/angprm/aprp(500),eav3(1000),
     1eav2(1000),pot2x,th1,ph1,snt(360),cst(360),snp(720),csp(720)
     2,pot
      common/mmtcom/mmt,mzee,iares,izres,light
        ed2=ed/2.
       iehol=(ehol+ed)/ed
       if(iehol.eq.0)iehol=1
       nscat=nhol
       x=rndm(-1.)
       if(x.gt.pxx(nhol,iehol).and.nhol.eq.1)nscat=2
       if(x.gt.pxx(nhol,iehol).and.nhol.eq.2)nscat=1
        kz=jz+idz(nscat)
       ka=ja+ida(nscat)
       iares=mmt-ja-jz
       izres=mzee-jz
c      if(iares.lt.aferm)then
c      write(*,*)'fermi stats called'
c      call gudima(ein,kst)
c      call fermfrag(kst)
c      light=1
c      return
c      endif
c      if(izres.lt.3)call litup(light,ehol)
c      if(izres.lt.6.and.iares.lt.9)call litup(light,ehol)
c      if(light.gt.0)return
       if(kz.gt.13)kz=13
       if(ka.gt.22)ka=22
       x=rndm(-1.)
       eff=ehol-pairmc(kz,ka)
       if(eff.gt.v)eff=v
       enscat=eff*(1.-sqrt(1.-x))
 
       if(enscat.gt.ehol)enscat=ehol
       elabh=enscat-be(jz,ja,nscat)
       if(elabh.le.ed)go to 100
       ieh=(elabh+ed)/ed
       if(rndm(-1.).gt.q1(nscat,ieh))return
       ja=ja+ida(nscat)
       jz=jz+idz(nscat)
       call rotang3(ehol,enscat,thethol,phihol,elabh)
 
       ith=thethol*11.459155+1.
       if(ith.eq.37)ith=36
       scrs(ieh,nscat)=scrs(ieh,nscat)+fcspec
       dsgp(ith,ieh,nscat)=dsgp(ith,ieh,nscat)+fcspec
c sept 02 put this emission into history file
       if(iend.eq.1)call endf(thethol,elabh,ith,ieh,nscat )
       eloss=eloss+enscat
        kz=jz+idz(nhol)
       ka=ja+ida(nhol)
       if(kz.gt.12)kz=12
       if(ka.gt.22)ka=22
c       ek1=enscat
       nu=nscat
       rmas=amass+2.-float(ka)-float(kz)
       ecm=elabh*rmas/(rmas+1.)
       pcm=sqrt(2.*elabh)
       qcm=sqrt(2.*(ecm-eoff(nscat)))
       tspin=spincon*qcm
       delj=tspin*cos(thethol)
       delprp=delj*sin(thethol)*sin(phihol)/cos(thethol)
       spinprp=spinprp+delprp
       spinf=spinf+delj
       if(kinem.eq.1)pcm=sqrt(2.*ecm)
       nu=nscat
       if(tspin.gt.30.)write(*,*)tspin,'tspin'
c  kkg  04/13/11
      if(rmas.lt.0.1)  then 
      write(99,*) 'HOLSCAS1(irecoi=9):rmas,amass,ka,kz=',rmas,amass,ka,
     1kz 
      return 
      endif       
       	if(rmas.ge.1.)call recoil(14)
c
c
  100  continue
       return
       end
      subroutine holscat(nhol,ehol,thethol,phihol)
      common/zamax/kzmax,kamax
      common/out1/nout
      common/histindx/ihist
      common/history/ide(60),idn(60),idith(60),iddum(60),the(60),adelj
     1(60),spinfh(60),enghis(60)
      common/endvparm/iend,lll(6,6,6),ll(6,6)
      common/lout/spincon
      common/perp/spinprp
      common/ilcom/nu,ndum,am(8)
      common/edep/pxx(2,1100)
      common/exclbe/idz(4),ida(4),beff(20,20,2),pairmc(18,27)
      common/scat/jz,ja, probxx(2,2),q1(8,1100),fcs,eloss,fcspec
      common/memo/scrs(999,8)
      common/sft5/exc(18,27),xmxx
      common/lym1/be(18,27,8),pair(18,27),xmas(18,27),delshl(18,27),amas
     1(18,27)
c     common/bym1/symb(18,27),symbp(18,27),symbs(18,27)
      common/pl3/jx,nepr,eb(100),rcsp(100),zee,amass,plex
      common/incr/ed,tened,efinal
      common/eder/ed2,v,vhalf,vsq
      common/shft/mc,mp,inver,ike,ipch,kq,qval,ap,at,zp,zt,cld,barfac
      common/angprm/aprp(500),eav3(1000),
     1eav2(1000),pot2x,th1,ph1,snt(360),cst(360),snp(720),csp(720)
     2,pot
c
      kzmax=zt+zp
      kamax=at+ap-zt-zp
       iehol=(ehol+ed)/ed
       if(iehol.eq.0)iehol=1
       nscat=nhol
       x=rndm(-1.)
       if(x.gt.pxx(nhol,iehol).and.nhol.eq.1)nscat=2
       if(x.gt.pxx(nhol,iehol).and.nhol.eq.2)nscat=1
        kz=jz+idz(nscat)
       ka=ja+ida(nscat)
        if(kz.gt.13)kz=13
        if(ka.gt.22)ka=22
c 11/8/2010
        if(kz.gt.kzmax .or. ka.gt.kamax)then
          eff=0.
           en1lab=0.
           go to 100
            end if
       x=rndm(-1.)
       eff=ehol-pairmc(kz,ka)
       if(eff.gt.v)eff=v
       enscat=eff*(1.-sqrt(1.-x))
       if(enscat.gt.ehol)enscat=ehol
       elabh=enscat-be(jz,ja,nscat)
       if(elabh.le.ed/2.)go to 100
       ieh=(elabh+ed)/ed
       if(rndm(-1.).gt.q1(nscat,ieh))return
       ja=ja+ida(nscat)
       jz=jz+idz(nscat)
       scrs(ieh,nscat)=scrs(ieh,nscat)+fcspec
       if(iend.eq.1)then
        theta=rndm(-1.)*3.14159
        ith=theta*57.296/36.+1.
        if(ith.gt.36)ith=36
       if(iend.eq.1)call endf(theta,elabh,ith,ieh,nscat )
       endif
       eloss=eloss+enscat
        kz=jz+idz(nhol)
       ka=ja+ida(nhol)
  100  continue
       return
       end
c
      subroutine sortmc
      common/outlim/limout,lim2
      common/nuopt/kinem,iem,kb
      common/isomer1/xgs(690),xiso(690),x2iso(690),isonum,isoint,iso2num
      common/isomer2/zs(690),as(690),eiso(690),eiso2(690)
      common/isomer3/iziso(690),iaiso(690),csgs(690),csis(690)
     1,csgse(690,100),csise(690,100),csis2e(690,100),rm1(690,100),
     2rm2(690,100)
      common/tem/philabb(36),phirec(36),phb(36),phin(36)
      common/hyb2/pp(15,24,999)
      common/labcon/prx,pry,prz,pcm,rmas,clab(999,8),dslb(36,999,8)
      common/unit8/angles(9)
      common/mbcfudg/cfg(999)
        common/mccom/ncount,ipot,rum(15,24),ppmc(15,24,999),cf
      common/exclus/ gib(8,999),pr(3),ge(2100),zz(3),
     1delta(3),ts(3),pairx(8),asum(10),rate(8,999)
      common/shft/mc,mp,inver,ike,ipch,kq,qval,ap,at,zp,zt,cld,barfac
      common/sft5/exc(18,27),xmax
      common/pl3/jx,nepr,eb(100),rcsp(100),zee,amass,plex
      common/par3/ep,sigml(999),acrs(999)
      common/sf6/aoz,zoz,en(8,999)
      common/memo/scrs(999,8)
      common/incr/ed,tened,efinal
      common/send/irfr,iadst,dlt
      common/dist/n,ni,ke,irefr,xz(3),pod,pcs,ek,xmix,cfrac
     1,dtau,sign(36),tmix,bji(10),fin,dtup
      common/error/rerror(15,24)
      common/eder/ed2,v,vhalf,vsq
      common/out1/nout
      common/out2/vrl(200,25),vrlz(200,25),vrlx(200,25),vrly(200,25)
      common/out3/vscale
      common/out4/asm(36),asmm(36)
crc
      common/distb/dsgp(36,999,8)
      common/recoils/rcletza(100,13,24,18),rcleza(100,13,24),
     1rclea(100,36),rcleat(100,36,18),rphir(36)
     2,rclsp(200),rclz(200),rcly(200),rclx(200),escale
     3,rcle(200),rclt(36),rclet(200,36)
      common/rec1/ermxc,thetrc
      common/scat/jz,ja, probxx(2,2),q1(8,1100),fcs,eloss,fcspec
crc
       do 700 jj=1,15
       do 690 jja=1,24
       rum(jj,jja)=0.
  690  continue
  700  continue
       ixmx=(xmax+ed )/ed
       gsum=0.
       do 7000 jj=1,15
       do 6500 jja=1,24
       do 6400 i=1,ixmx
       rum(jj,jja)=rum(jj,jja)+pp(jj,jja,i)
6400   continue
       rerror(jj,jja)=sqrt(rum(jj,jja)*fcs)
       gsum=gsum+rum(jj,jja)
6500   continue
7000   continue
       write(17,7008)gsum
       pp(1,1,ixmx)=0.
       write(17,204)
       write(17,7011)
7011   format(' equilibrated cascade yields returned for weisskopf
     1 decay')
       write(17,204)
7008   format(' sum over all events = ',f10.5)
       write(17,7009)
       write(17,204)
7009   format(' Z indx N-no.= N       N-1     N-2       N-3       N-4
 
     1    N-5       N-6       N-7       N-8       N-9       N-10    '
 
     2,'     N-11')
       do 7001 i =1,15
 
        do j=24,1,-1
         if(rum(i,j).ne.0.0) then
          k = j
          goto 6999
         end if
        end do
        goto 7001
 
 6999      l=i-1
        write(17,201)l,(rum(i,j),j=1,k)
       write(17,203)(rerror(i,j),j=1,k)
       write(17,204)
  204 format('  ')
 
 7001   continue
 
       write(17,7010)
7010   format(' error is +/- one sigma statistical error')
 201   format(' Z-',i2,3x,10(1x,1pe9.3))
 203   format('+ - sd =',10(1x,1pe9.3))
c 402   format(1x,f5.1,1x,10(1p,e8.1))
c dec98 defer write of ddcs until equil cs added
       call reclvel
c now add outputs on various recoil energy arrays.
c make subroutine recloilwr for this
       call recoilwr(ermxc)
crc
       return
       end
c
c
c*s
c   NAME OF ROUTINE: ranf
c
c   CALLING ROUTINE: most source and transport routines
c
c   AUTHOR:  D. A. Resler, LLNL, September 29, 1994
c
c   PURPOSE:
c
c     This routine produces a new random number every time it is
called
c
c     entry point rnset:
c
      REAL*8 function ranf()
c
c   NAME OF ROUTINE: ranf
c
c   CALLING ROUTINE: most source and transport routines
c
c   AUTHOR:  D. A. Resler, LLNL, September 29, 1994
c
c   PURPOSE:
c
c     This routine produces a new random number every time it is called
c
c     entry point rnset:
c
c       Provides different seeds to the random number generator.
c
c   REFERENCES:
c
c      George Marsaglia and Arif Zaman,
c      "Some portable very-long-period random number generators"
c      Computers in Physics,
c      Vol. 8, No. 1, Jan/Feb 1994, p. 117.
c
c   RESTRICTIONS:
c
c   SUBROUTINES AND FUNCTIONS CALLED:
c
c   GLOSSARY:
c
c     rnset:
c       indx = the index specifying the set of seeds to use
c              currently indx must be between 0 and 200 inclusive
c
c       rnset = 0d0  -->  seed set OK
c               1d0  -->  invalid value for 'indx'
c
c
c   CHANGES:
c
c 09/29/94 chs Initiated in PEREGRINE
c 10/25/94 dar Add entry point to set seeds
c 01/31/95 llc Use IMPLICIT NONE instead of IMPLICIT UNDEFINED (A-Z);
c              later is not supported by DEC Alpha OSF/1 compiler.
c 11/03/95 smh Added return to ranset (for clarity),
c              Added entry point ranf4 , return sp number uniform on ( 0,+1)
c              Added entry point ranfu , return    number uniform on (-1,+1)
c              Added entry point ranf4u, return sp number uniform on (-1,+1)
c              Note that normalization factors return numbers on `open` intervals
c              Note that physics requires open interval
c 12/15/95 smh Switched to independent functions and a block data sub.
c
c*f
 
      common /ran_com/ i , j , k , n
      integer          i , j , k , n
 
      integer mzran
 
c ......................................................................
 
      mzran = i - k
 
      IF (mzran.lt.0) mzran = mzran + 2147483579
 
      i = j
      j = k
      k = mzran
      n = 69069 * n + 1013904243
 
c ... REAL*8 on open interval (0,1)
 
      ranf = 0.5d0 + 0.232830643d-9 * (mzran + n)
 
      return
 
      end
 
cccccccccccccccccccccccccccccccccccccccccc
 
      REAL*8 function ranfu()
 
      common /ran_com/ i , j , k , n
      integer          i , j , k , n
 
      integer mzran
 
c ......................................................................
 
      mzran = i - k
 
      if (mzran.lt.0) mzran = mzran + 2147483579
 
      i = j
      j = k
      k = mzran
      n = 69069 * n + 1013904243
 
c ... REAL*8 on open interval (0,1)
 
      ranfu = 0.465661287d-9 * (mzran + n)
 
      return
 
      end
 
cccccccccccccccccccccccccccccccccccccccccc
 
      REAL function ranf4u()
 
      common /ran_com/ i , j , k , n
      integer          i , j , k , n
 
      integer mzran
 
c ......................................................................
 
      mzran = i - k
 
      if (mzran.lt.0) mzran = mzran + 2147483579
 
      i = j
      j = k
      k = mzran
      n = 69069 * n + 1013904243
 
c ... REAL on open interval (0,1)
 
      ranf4u = 0.4656612e-9 * (mzran + n)
 
      return
 
      end
 
cccccccccccccccccccccccccccccccccccccccccc
 
      REAL*8 function rnset(indx)
 
      common /ran_com/ i , j , k , n
      integer          i , j , k , n
 
      common /rans_com/ is, js, ks, ns
      integer is(0:200), js(0:200), ks(0:200), ns(0:200)
 
      integer indx
 
c ......................................................................
 
      if (indx.lt.0.or.indx.gt.200) then
 
        rnset = 1.0d0
 
      else
 
c ..... set the seeds i, j, k, n
 
        i = is(indx)
        j = js(indx)
        k = ks(indx)
        n = ns(indx)
 
        rnset = 0.0d0
 
      endif
 
      return
 
      end
 
cccccccccccccccccccccccccccccccccccccccccc
 
      subroutine rnget(iseed, jseed, kseed, nseed)
 
      common /ran_com/ i , j , k , n
      integer          i , j , k , n
 
      integer iseed, jseed, kseed, nseed
 
c ......................................................................
 
c ... get the seeds i, j, k, n
 
      iseed = i
      jseed = j
      kseed = k
      nseed = n
 
      return
 
      end
 
cccccccccccccccccccccccccccccccccccccccccc
 
      block data rancom
 
      common /ran_com/ i , j , k , n
      integer          i , j , k , n
 
      common /rans_com/ is, js, ks, ns
      integer is(0:200), js(0:200), ks(0:200), ns(0:200)
 
c ... The seeds were generated by calling a test version of 'ranf'
c     2 trillion times.  At every 10 billion calls, the seed information
c     was recorded in a file for inclusion here.
 
      data i,j,k,n/521288629,362436069,16163801,1131199299/
 
      data   is(  0),    js(  0),    ks(  0),    ns(  0)
     >  /  521288629,  362436069,   16163801, 1131199299/
      data   is(  1),    js(  1),    ks(  1),    ns(  1)
     >  / 1013417550, 2021587047,  793237685, -308031677/
      data   is(  2),    js(  2),    ks(  2),    ns(  2)
     >  /  743256993,  362093083,  339630377,-2028281021/
      data   is(  3),    js(  3),    ks(  3),    ns(  3)
     >  / 2049403850, 1278762488, 2071188061,  265418563/
      data   is(  4),    js(  4),    ks(  4),    ns(  4)
     >  /  696800831, 1146290771, 1501321453,-2016867517/
      data   is(  5),    js(  5),    ks(  5),    ns(  5)
     >  / 1420161510,  978067693,  569455188, -285204669/
      data   is(  6),    js(  6),    ks(  6),    ns(  6)
     >  / 1696213944, 1323660753,  992991869, 1165439811/
      data   is(  7),    js(  7),    ks(  7),    ns(  7)
     >  / 2095931264, 1377640598,  548487367,-1959901373/
      data   is(  8),    js(  8),    ks(  8),    ns(  8)
     >  /  819840680, 1627780128, 1722969211,-1071293629/
      data   is(  9),    js(  9),    ks(  9),    ns(  9)
     >  / 1388408751, 1461698873,  993331411, -463704253/
      data   is( 10),    js( 10),    ks( 10),    ns( 10)
     >  /  915073478, 2122472120, 1590153721, -137133245/
      data   is( 11),    js( 11),    ks( 11),    ns( 11)
     >  / 1728824372, 1484230701, 1225420891,  -91580605/
      data   is( 12),    js( 12),    ks( 12),    ns( 12)
     >  / 1932693086, 1152830081,  122608524, -327046333/
      data   is( 13),    js( 13),    ks( 13),    ns( 13)
     >  / 1111388652,  192759078, 1298310681, -843530429/
      data   is( 14),    js( 14),    ks( 14),    ns( 14)
     >  /  269917017,  286668054, 1254138045,-1641032893/
      data   is( 15),    js( 15),    ks( 15),    ns( 15)
     >  / 1766961713,  926119556, 1867800318, 1575413571/
      data   is( 16),    js( 16),    ks( 16),    ns( 16)
     >  /   13021851,  415274785,  855468610,  215874371/
      data   is( 17),    js( 17),    ks( 17),    ns( 17)
     >  / 1216791231, 1457651852,  799260263,-1424683197/
      data   is( 18),    js( 18),    ks( 18),    ns( 18)
     >  / 1254200136, 1251144008, 1247686438,  948708163/
      data   is( 19),    js( 19),    ks( 19),    ns( 19)
     >  /  866821275, 1485996330, 2003827084,-1253886141/
      data   is( 20),    js( 20),    ks( 20),    ns( 20)
     >  / 1856807102,  297236216, 1589356795,  557468483/
      data   is( 21),    js( 21),    ks( 21),    ns( 21)
     >  / 1574235846,  438824661, 1508522588, 2087804739/
      data   is( 22),    js( 22),    ks( 22),    ns( 22)
     >  / 1708591214,  928587522, 1625435021, -957844669/
      data   is( 23),    js( 23),    ks( 23),    ns( 23)
     >  /  748103857, 1401908328,  714636665,   10454851/
      data   is( 24),    js( 24),    ks( 24),    ns( 24)
     >  / 1473780165, 2087909014, 1366364078,  697736003/
      data   is( 25),    js( 25),    ks( 25),    ns( 25)
     >  / 2036218191,  542848634, 1126727648, 1103998787/
      data   is( 26),    js( 26),    ks( 26),    ns( 26)
     >  / 2085840008, 1722332183,  229208295, 1229243203/
      data   is( 27),    js( 27),    ks( 27),    ns( 27)
     >  /  368392468, 1821288836, 1225876574, 1073469251/
      data   is( 28),    js( 28),    ks( 28),    ns( 28)
     >  / 1200591972,  828065665,  387107710,  636676931/
      data   is( 29),    js( 29),    ks( 29),    ns( 29)
     >  /   42389948, 2038629135, 1922438844,  -81133757/
      data   is( 30),    js( 30),    ks( 30),    ns( 30)
     >  /  858155273,  455694602, 2078953128,-1079962813/
      data   is( 31),    js( 31),    ks( 31),    ns( 31)
     >  / 1567487196, 1171424283, 1063625771, 1935157059/
      data   is( 32),    js( 32),    ks( 32),    ns( 32)
     >  /  343029734, 1717280897,  697755387,  374291267/
      data   is( 33),    js( 33),    ks( 33),    ns( 33)
     >  / 1588700217, 1094676910,  613708470,-1467592893/
      data   is( 34),    js( 34),    ks( 34),    ns( 34)
     >  / 1831456196, 1139379080,   42765459,  704471875/
      data   is( 35),    js( 35),    ks( 35),    ns( 35)
     >  /  527243397, 1242990401, 1742567071,-1699449021/
      data   is( 36),    js( 36),    ks( 36),    ns( 36)
     >  / 2014510677,  331693581,   44538930,  -89420989/
      data   is( 37),    js( 37),    ks( 37),    ns( 37)
     >  /  950881674,  914431491,  414494926, 1239588675/
      data   is( 38),    js( 38),    ks( 38),    ns( 38)
     >  / 1679762001, 1279550335,  207467822,-2007387325/
      data   is( 39),    js( 39),    ks( 39),    ns( 39)
     >  /  322161359, 1071627473,  386266493,-1240414397/
      data   is( 40),    js( 40),    ks( 40),    ns( 40)
     >  /  452128837,  973586925,  101927879, -754459837/
      data   is( 41),    js( 41),    ks( 41),    ns( 41)
     >  /  840767649, 1538223507,  666817672, -549523645/
      data   is( 42),    js( 42),    ks( 42),    ns( 42)
     >  / 1469750540,  381270084, 2092014849, -625605821/
      data   is( 43),    js( 43),    ks( 43),    ns( 43)
     >  / 1945068262, 1819911439,  222156715, -982706365/
      data   is( 44),    js( 44),    ks( 44),    ns( 44)
     >  /  124972818, 1625607724,  389916186,-1620825277/
      data   is( 45),    js( 45),    ks( 45),    ns( 45)
     >  /  382205210, 1858695634,  680639953, 1755004739/
      data   is( 46),    js( 46),    ks( 46),    ns( 46)
     >  /  773251578,  569290382,  682797090,  554849091/
      data   is( 47),    js( 47),    ks( 47),    ns( 47)
     >  /  744042674, 1603033885,  575162700, -926324925/
      data   is( 48),    js( 48),    ks( 48),    ns( 48)
     >  /  868575903,  293478224,  133687189, 1606449987/
      data   is( 49),    js( 49),    ks( 49),    ns( 49)
     >  / 2017687168,  867922419,  260002558, -436760765/
      data   is( 50),    js( 50),    ks( 50),    ns( 50)
     >  / 1771509433, 1098924612,  179705323, 1533977411/
      data   is( 51),    js( 51),    ks( 51),    ns( 51)
     >  /  300765823,  954800766, 1300094664,-1071270077/
      data   is( 52),    js( 52),    ks( 52),    ns( 52)
     >  /  420104577,  922261311,  257130586,  337431363/
      data   is( 53),    js( 53),    ks( 53),    ns( 53)
     >  /  657282786,  292983851, 1321708006, 1465114435/
      data   is( 54),    js( 54),    ks( 54),    ns( 54)
     >  /  735242062, 1374009347,  114715275,-1983188157/
      data   is( 55),    js( 55),    ks( 55),    ns( 55)
     >  / 1896578982, 1645241098, 1134131111,-1417541821/
      data   is( 56),    js( 56),    ks( 56),    ns( 56)
     >  /  679687938,  983370975,  782953306,-1132913853/
      data   is( 57),    js( 57),    ks( 57),    ns( 57)
     >  / 1405118890, 1696417556,  543549592,-1129304253/
      data   is( 58),    js( 58),    ks( 58),    ns( 58)
     >  / 1840739087,  966234354, 1962685837,-1406713021/
      data   is( 59),    js( 59),    ks( 59),    ns( 59)
     >  /  634026637,  951934551, 1720678670,-1965140157/
      data   is( 60),    js( 60),    ks( 60),    ns( 60)
     >  / 1350601724, 1125406639,  165185997, 1490381635/
      data   is( 61),    js( 61),    ks( 61),    ns( 61)
     >  /  463603434, 1951530915,  113203881,  369917763/
      data   is( 62),    js( 62),    ks( 62),    ns( 62)
     >  /  179522345,  522484601,  126299736,-1031564477/
      data   is( 63),    js( 63),    ks( 63),    ns( 63)
     >  / 1291383330,  786861304, 1611444374, 1580902211/
      data   is( 64),    js( 64),    ks( 64),    ns( 64)
     >  / 1491455160, 1951921038,  198136012, -382616765/
      data   is( 65),    js( 65),    ks( 65),    ns( 65)
     >  / 1726674313,   72257224,  306998123, 1667813187/
      data   is( 66),    js( 66),    ks( 66),    ns( 66)
     >  /  503437758, 1477539377,  799241873, -857742525/
      data   is( 67),    js( 67),    ks( 67),    ns( 67)
     >  /  517717256, 1192183549,  796958292,  630650691/
      data   is( 68),    js( 68),    ks( 68),    ns( 68)
     >  /  346161708,  664229615, 1331536856, 1838025539/
      data   is( 69),    js( 69),    ks( 69),    ns( 69)
     >  / 1009230873, 1999412387,  893103246,-1530585277/
      data   is( 70),    js( 70),    ks( 70),    ns( 70)
     >  / 2094226656, 1728554803, 1163542708, -885247165/
      data   is( 71),    js( 71),    ks( 71),    ns( 71)
     >  / 1226321866,  958734844, 1310546555, -520927421/
      data   is( 72),    js( 72),    ks( 72),    ns( 72)
     >  /  630392817,  853533098, 1063448412, -437626045/
      data   is( 73),    js( 73),    ks( 73),    ns( 73)
     >  /  934799378,  899852415,  807383913, -635343037/
      data   is( 74),    js( 74),    ks( 74),    ns( 74)
     >  /  134490037, 1709931467,   16791313,-1114078397/
      data   is( 75),    js( 75),    ks( 75),    ns( 75)
     >  / 2099097490,   34006906, 1519338347,-1873832125/
      data   is( 76),    js( 76),    ks( 76),    ns( 76)
     >  / 1195140651, 1264339087, 1816228238, 1380363075/
      data   is( 77),    js( 77),    ks( 77),    ns( 77)
     >  / 1917857987,  632056288, 1618471385,   58572611/
      data   is( 78),    js( 78),    ks( 78),    ns( 78)
     >  / 1944941534, 2087394370, 1463139173,-1544236221/
      data   is( 79),    js( 79),    ks( 79),    ns( 79)
     >  /  151608279, 1225988314,  643000034,  866903875/
      data   is( 80),    js( 80),    ks( 80),    ns( 80)
     >  / 1424256193, 1578912300,  335986147,-1297941693/
      data   is( 81),    js( 81),    ks( 81),    ns( 81)
     >  / 1241417776, 1485324956, 1107423106,  551161667/
      data   is( 82),    js( 82),    ks( 82),    ns( 82)
     >  / 2137821368,  576353659, 1880537674, 2119246659/
      data   is( 83),    js( 83),    ks( 83),    ns( 83)
     >  /  904595386,  908973260,  550148161, -888654013/
      data   is( 84),    js( 84),    ks( 84),    ns( 84)
     >  /  238162371, 1111888832, 1041062287,  117394243/
      data   is( 85),    js( 85),    ks( 85),    ns( 85)
     >  / 1406447457,  638264260, 1700136421,  842424131/
      data   is( 86),    js( 86),    ks( 86),    ns( 86)
     >  / 1289653512,  960770339,  690376939, 1286435651/
      data   is( 87),    js( 87),    ks( 87),    ns( 87)
     >  /  699188683, 2079430653, 1553015440, 1449428803/
      data   is( 88),    js( 88),    ks( 88),    ns( 88)
     >  / 1034091741, 1087681807,  834856837, 1331403587/
      data   is( 89),    js( 89),    ks( 89),    ns( 89)
     >  / 1472564639, 1079527482,  206795572,  932360003/
      data   is( 90),    js( 90),    ks( 90),    ns( 90)
     >  / 1935638959,  979979222,  334871275,  252298051/
      data   is( 91),    js( 91),    ks( 91),    ns( 91)
     >  / 1289670849, 2138224005, 1323212690, -708782269/
      data   is( 92),    js( 92),    ks( 92),    ns( 92)
     >  / 1961524036,  199219564,  568432161,-1950880957/
      data   is( 93),    js( 93),    ks( 93),    ns( 93)
     >  /  324326725, 2136714753,  491750884,  820969283/
      data   is( 94),    js( 94),    ks( 94),    ns( 94)
     >  /  787059924, 1332899615,  895340626, -983166141/
      data   is( 95),    js( 95),    ks( 95),    ns( 95)
     >  / 1690806733,  440486884,  273618507, 1226647363/
      data   is( 96),    js( 96),    ks( 96),    ns( 96)
     >  / 2083018261, 1459237202,  721080211,-1139524797/
      data   is( 97),    js( 97),    ks( 97),    ns( 97)
     >  / 1892512131, 1491767918,  779812658,  508251971/
      data   is( 98),    js( 98),    ks( 98),    ns( 98)
     >  /  953198476,  991072060, 2130685652, 1875010371/
      data   is( 99),    js( 99),    ks( 99),    ns( 99)
     >  /  276189351, 2078036710, 1227930095,-1334216893/
      data   is(100),    js(100),    ks(100),    ns(100)
     >  / 2009757041, 1555488289, 1725690622, -529495229/
      data   is(101),    js(101),    ks(101),    ns(101)
     >  / 1014564040, 1896201581, 2090485145,   -5791933/
      data   is(102),    js(102),    ks(102),    ns(102)
     >  / 1691526742,  527135590,  262660994,  236892995/
      data   is(103),    js(103),    ks(103),    ns(103)
     >  /  546881082,  281524268, 1850533317,  198559555/
      data   is(104),    js(104),    ks(104),    ns(104)
     >  /   84548922,  144208096, 1386439813, -120792253/
      data   is(105),    js(105),    ks(105),    ns(105)
     >  / 1153077296,  551354209, 1392554801, -721162429/
      data   is(106),    js(106),    ks(106),    ns(106)
     >  /  357621270, 1091169784, 1494378352,-1602550973/
      data   is(107),    js(107),    ks(107),    ns(107)
     >  / 1015851432,  179401555, 1071853537, 1530009411/
      data   is(108),    js(108),    ks(108),    ns(108)
     >  /  726696319, 2047479224, 1380101434,   86584131/
      data   is(109),    js(109),    ks(109),    ns(109)
     >  /  129328670, 1716500020, 1297444061,-1637859517/
      data   is(110),    js(110),    ks(110),    ns(110)
     >  / 1260139593,  983290131,  152792373,  651645763/
      data   is(111),    js(111),    ks(111),    ns(111)
     >  /  392112457,  346051619,  999943340,-1634834621/
      data   is(112),    js(112),    ks(112),    ns(112)
     >  /  123679258,  362122696, 1746731684,   92633923/
      data   is(113),    js(113),    ks(113),    ns(113)
     >  /  882651809,  524071327,  560814956, 1539084099/
      data   is(114),    js(114),    ks(114),    ns(114)
     >  / 1357437319, 1023637502,  381128610,-1590451389/
      data   is(115),    js(115),    ks(115),    ns(115)
     >  / 2063506503, 1765623623, 1718838704, -706037949/
      data   is(116),    js(116),    ks(116),    ns(116)
     >  / 1090899958,  544579423, 1702385605, -102642877/
      data   is(117),    js(117),    ks(117),    ns(117)
     >  /  274905376, 1559375026, 1702475062,  219733827/
      data   is(118),    js(118),    ks(118),    ns(118)
     >  /  267310444, 1975961023, 1944634884,  261092163/
      data   is(119),    js(119),    ks(119),    ns(119)
     >  / 1305308308, 1167660155, 1195518880,   21432131/
      data   is(120),    js(120),    ks(120),    ns(120)
     >  / 2146161879, 1192995066,   12956994, -499246269/
      data   is(121),    js(121),    ks(121),    ns(121)
     >  / 1809437720, 1719963222,  779421674,-1300943037/
      data   is(122),    js(122),    ks(122),    ns(122)
     >  /  967194072, 1041128112,  897109411, 1911309123/
      data   is(123),    js(123),    ks(123),    ns(123)
     >  / 1663628879, 1351226721,  691981404,  547575619/
      data   is(124),    js(124),    ks(124),    ns(124)
     >  /  665474344, 1098337700, 1191519172,-1097176253/
      data   is(125),    js(125),    ks(125),    ns(125)
     >  / 2145980375, 1995308825, 1195497901, 1272020803/
      data   is(126),    js(126),    ks(126),    ns(126)
     >  /  739473203,  201370171,  448346140, -934767805/
      data   is(127),    js(127),    ks(127),    ns(127)
     >  / 1052261745,  261227100,  917465567,  872392515/
      data   is(128),    js(128),    ks(128),    ns(128)
     >  /  669351975, 1995310693, 1958736665,-1896432829/
      data   is(129),    js(129),    ks(129),    ns(129)
     >  /  761375636,  918014304,  422761964, -651309245/
      data   is(130),    js(130),    ks(130),    ns(130)
     >  /  783769283, 1694674718, 1428289344,  312795971/
      data   is(131),    js(131),    ks(131),    ns(131)
     >  / 1841657625, 1562041223,   29205595,  995882819/
      data   is(132),    js(132),    ks(132),    ns(132)
     >  /   19734051, 1295139672, 1170776359, 1397951299/
      data   is(133),    js(133),    ks(133),    ns(133)
     >  /  944303723,  130899444,  417717848, 1519001411/
      data   is(134),    js(134),    ks(134),    ns(134)
     >  / 1108130801, 1446138813,  472339326, 1359033155/
      data   is(135),    js(135),    ks(135),    ns(135)
     >  / 1440463259, 1895759879, 1190063174,  918046531/
      data   is(136),    js(136),    ks(136),    ns(136)
     >  / 1995946762,  435003548,  840120935,  196041539/
      data   is(137),    js(137),    ks(137),    ns(137)
     >  /  564551504,  741552439,   53498545, -806981821/
      data   is(138),    js(138),    ks(138),    ns(138)
     >  /   21935489, 1277688289,  398952870,-2091023549/
      data   is(139),    js(139),    ks(139),    ns(139)
     >  /  733093335,  331568157, 1042185968,  638883651/
      data   is(140),    js(140),    ks(140),    ns(140)
     >  /  143557321,  193796695,  117369009,-1207194813/
      data   is(141),    js(141),    ks(141),    ns(141)
     >  /  727557789,  652610570, 2077920346,  960675651/
      data   is(142),    js(142),    ks(142),    ns(142)
     >  /  828345684, 1057716030, 1379723676,-1447439549/
      data   is(143),    js(143),    ks(143),    ns(143)
     >  /  416752718,  161127010,  630304153,  158394179/
      data   is(144),    js(144),    ks(144),    ns(144)
     >  /   15730625,  893179360,  716484453, 1483209539/
      data   is(145),    js(145),    ks(145),    ns(145)
     >  / 1088077911,  422216780, 1719749836,-1767960765/
      data   is(146),    js(146),    ks(146),    ns(146)
     >  /  497451102, 2045536655,  225106351,-1005182141/
      data   is(147),    js(147),    ks(147),    ns(147)
     >  / 2018608972, 1197104204, 1801814030, -523421885/
      data   is(148),    js(148),    ks(148),    ns(148)
     >  /  306370044, 1474602286, 1764588819, -322679997/
      data   is(149),    js(149),    ks(149),    ns(149)
     >  / 1087380334, 1852148039,  154518015, -402956477/
      data   is(150),    js(150),    ks(150),    ns(150)
     >  / 2120459640, 1709027598, 2061856491, -764251325/
      data   is(151),    js(151),    ks(151),    ns(151)
     >  / 1831038150,  953110656,  294187074,-1406564541/
      data   is(152),    js(152),    ks(152),    ns(152)
     >  /  624987994, 1099601476,  532831116, 1965071171/
      data   is(153),    js(153),    ks(153),    ns(153)
     >  /  235920569, 1787552080, 1515395455,  760721219/
      data   is(154),    js(154),    ks(154),    ns(154)
     >  /  513737855, 1875791744,  571664731, -724647101/
      data   is(155),    js(155),    ks(155),    ns(155)
     >  /  744499401,  707560137, 1305708936, 1803933507/
      data   is(156),    js(156),    ks(156),    ns(156)
     >  /  638119166,  215275139,  862675565, -243471549/
      data   is(157),    js(157),    ks(157),    ns(157)
     >  /  985583543, 1769013121, 1141540686, 1723072323/
      data   is(158),    js(158),    ks(158),    ns(158)
     >  / 1662873212,  484373291, 1282462630, -886369469/
      data   is(159),    js(159),    ks(159),    ns(159)
     >  /  333913009, 1198987070,   62340298,  518137667/
      data   is(160),    js(160),    ks(160),    ns(160)
     >  /  380202123, 1819273579, 1967279667, 1641626435/
      data   is(161),    js(161),    ks(161),    ns(161)
     >  / 2146866740, 1603615015,  217996661,-1810870461/
      data   is(162),    js(162),    ks(162),    ns(162)
     >  /  718152277, 1490562980, 1019756025,-1249418429/
      data   is(163),    js(163),    ks(163),    ns(163)
     >  /  151670360,  520405074, 1481233944, -968984765/
      data   is(164),    js(164),    ks(164),    ns(164)
     >  / 1741862501,  687410730,  401810445, -969569469/
      data   is(165),    js(165),    ks(165),    ns(165)
     >  / 1557915042, 1610656902, 1170113139,-1251172541/
      data   is(166),    js(166),    ks(166),    ns(166)
     >  /  524778872, 1902501047, 1734332922,-1813793981/
      data   is(167),    js(167),    ks(167),    ns(167)
     >  /  774017850,  862321102,  353951141, 1637533507/
      data   is(168),    js(168),    ks(168),    ns(168)
     >  /  857970272,  215417267, 1553575508,  512875331/
      data   is(169),    js(169),    ks(169),    ns(169)
     >  / 1591097042,  612865728,  533555028, -892801213/
      data   is(170),    js(170),    ks(170),    ns(170)
     >  / 1167151568,  164790388, 1727527791, 1715471171/
      data   is(171),    js(171),    ks(171),    ns(171)
     >  /  553586074,  955635867,  789360780, -252242109/
      data   is(172),    js(172),    ks(172),    ns(172)
     >  /  253527164,  161595022, 1039955998, 1793993539/
      data   is(173),    js(173),    ks(173),    ns(173)
     >  / 1232709612,  530402430,  658913930, -735756477/
      data   is(174),    js(174),    ks(174),    ns(174)
     >  /  997706914, 1678873932, 1300069408,  748442435/
      data   is(175),    js(175),    ks(175),    ns(175)
     >  / 2026101010, 1500636615, 1042135968, 1951622979/
      data   is(176),    js(176),    ks(176),    ns(176)
     >  /  510373388,  345031046,  186909249,-1421182141/
      data   is(177),    js(177),    ks(177),    ns(177)
     >  /  939754062,  643658332, 1774147506, -780038333/
      data   is(178),    js(178),    ks(178),    ns(178)
     >  / 1598779809, 1162539976,  512364295, -419912893/
      data   is(179),    js(179),    ks(179),    ns(179)
     >  / 1358681365, 1820296401,  541401358, -340805821/
      data   is(180),    js(180),    ks(180),    ns(180)
     >  / 1475687713,  498817651,  684112869, -542717117/
      data   is(181),    js(181),    ks(181),    ns(181)
     >  / 1700007714,   70147687, 1575423863,-1025646781/
      data   is(182),    js(182),    ks(182),    ns(182)
     >  / 2023387475, 1250654222, 1591924661,-1789594813/
      data   is(183),    js(183),    ks(183),    ns(183)
     >  / 1108304872,  690557157, 1590592588, 1460406083/
      data   is(184),    js(184),    ks(184),    ns(184)
     >  /  944052160, 1857099499,  805939629,  134421315/
      data   is(185),    js(185),    ks(185),    ns(185)
     >  /  376379172,   89257943,  144968020,-1472581821/
      data   is(186),    js(186),    ks(186),    ns(186)
     >  /  146244701, 1967420015,  656522023,  934363971/
      data   is(187),    js(187),    ks(187),    ns(187)
     >  / 1323993072, 1250563403, 1778411148,-1234675901/
      data   is(188),    js(188),    ks(188),    ns(188)
     >  / 1572129676, 1213252917, 1655247997,  610233155/
      data   is(189),    js(189),    ks(189),    ns(189)
     >  /  817711696, 2137497831,  324637395,-2120843453/
      data   is(190),    js(190),    ks(190),    ns(190)
     >  / 1631797516,  867054282,  563735660, -837971133/
      data   is(191),    js(191),    ks(191),    ns(191)
     >  / 1825518888, 1475905597,   48796258,  163882819/
      data   is(192),    js(192),    ks(192),    ns(192)
     >  /  612671762, 1301109284, 1090264286,  884718403/
      data   is(193),    js(193),    ks(193),    ns(193)
     >  / 1860484777, 1216359013,  240192098, 1324535619/
      data   is(194),    js(194),    ks(194),    ns(194)
     >  / 1311552069,  920174858, 1583250034, 1483334467/
      data   is(195),    js(195),    ks(195),    ns(195)
     >  /  784327023, 1645274234,  549332955, 1361114947/
      data   is(196),    js(196),    ks(196),    ns(196)
     >  /   81634451, 2115485774, 1583439436,  957877059/
      data   is(197),    js(197),    ks(197),    ns(197)
     >  / 1203552899, 1051458461, 1034370209,  273620803/
      data   is(198),    js(198),    ks(198),    ns(198)
     >  /  399477757, 1087867726,  463957724, -691653821/
      data   is(199),    js(199),    ks(199),    ns(199)
     >  /  148045743, 2130714138, 1128779874,-1937946813/
      data   is(200),    js(200),    ks(200),    ns(200)
     >  /   72477969, 2064488191, 1911492604,  829709123/
 
      end
c-------------------------------------------------------------------------
c
c
c
      subroutine adisout
      common/outlim/limout,lim2
      dimension asum(10)
      common/sft5/exc(18,27),xmax
      common/send/irfr,iadst,dlt
      common/dist/n,ni,ke,irefr,xz(3),pod,pcs,ek,xmix,cfrac
     1,dtau,sign(36),tmix,bji(10),fin,dtup
      common/incr/ed,tened,efinal
      common/nuopt/kinem,iem,kb
      common/distb/dsgp(36,999,8)
      ed2=ed/2.
      confac=3.14*5./180.
      nend=xmax/(10.*dlt)+0.99
      do 520 l=1,2
c
      do 520 j=1,nend
      bji(1)=dlt*float(10*j-9)-ed/2.
      do 8406 ii=2,10
8406  bji(ii)=bji(ii-1)+dlt
c
      continue
c
      write(17,401)(bji(ii),ii=1,10)
      nii=j*10-9
      nif=nii+9
c
  401 format(' angle/deg. ke=',f6.2,9(4x,f6.2))
  402 format(1x,f5.1,5x,10(e10.3))
      do 6510 id=1,10
6510  asum(id)=0.
c
c
      do 510 i=1,2
      do 510 ith=1,36
c
      iii=0
      do 6511 kef=nii,nif
      iii=iii+1
6511  asum(iii)=asum(iii)+sign(ith)*dsgp(ith,kef,i)*confac
c
      angle=5.*ith-2.5
      write(17,402)angle,(dsgp(ith,kef,i),kef=nii,nif)
  510 continue
c
c
  520 continue
      return
      end
c
c
      subroutine pion
c
c pion from NN scattering versus E from Ver West and Arndt,Phys Rev.C25,1979
c     (1982)
c
      common/outlim/limout,lim2
      common/hybsub/iav,idum,rs(8,1100),c,pox,r1
      common/incr/ed,tened,efinal
      common/pies/sppplus(1100),sppzero(1100),snnmins(1100),snnzero(1100
     1),spnzero(1100),spnplus(1100),spnmins(1100),pxxpi(1100),pxypi(1100
     2),fppplus(1100),fnnmins(1100),fpnplus(1100),fpnzero(1100)
     3,ppinel(1100),pninel(1100)
      common/angprm/aprp(500),eav3(1000),
     1eav2(1000),pot2x,thet1,phi1,snt(360),cst(360),snp(720),csp(720)
     2,pot
      do 50 ie=1,1100
      ppinel(ie)=0.
      pninel(ie)=0.
      sppplus(ie)=0.
      sppzero(ie)=0.
      snnmins(ie)=0.
      snnzero(ie)=0.
      spnzero(ie)=0.
      spnplus(ie)=0.
      spnmins(ie)=0.
      pxxpi(ie)=0.
      pxypi(ie)=0.
      fppplus(ie)=0.
      fnnmins(ie)=0.
      fpnplus(ie)=0.
      fpnzero(ie)=0.
   50 continue
c
      imin=300./ed
      do 100 ie=imin,1100
      e1=float(ie)*ed+pot
      bta=sqrt(1.-(939.516/(939.516+e1))**2)
      aa=(10.63/bta**2-29.92/bta+42.9)
      bb=(34.1/bta**2-82.2/bta+82.2)
c     write(17,5)tl
      g0=134.3
      em0=1203.
      empi=138.0
      emn=938.9
      emd=2.*emn
      alpha=6.03
      beta=1.700
      gsq=g0*g0
      s=4.*emn*emn+2.*emn*e1
      psq=s/4.-emn*emn
      spin=(sqrt(s)-emn)**2
      prsqs=(s-(emd-empi)**2)*(s-(emd+empi)**2)/(4.*s)
      prs=sqrt(prsqs)
      s0=(emn+em0)**2
      p0sq=s0/4.-emn**2
      p0=sqrt(p0sq)
      sigd=(1./psq*2.)*alpha*((prs/p0)**beta)*(em0*em0*gsq)/((spin-em0*em0)
     1m0)**2+em0*em0*gsq)
      sigd=sigd*1.22e+6
c   divide by 4
      sigd=sigd/4.
c begin 3 body states
c
      alpha=15.28
      beta=0.
      em0c=1220.
      em0=1245.
      gam=137.4
      g0=120.
      zplus=(sqrt(s)-emn-em0c)*(2./g0)
      zminus=(emn+empi-em0c)*(2./g0)
      emav=em0c+(1./(atan(zplus)-atan(zminus)))*
     1(g0/4.)*log((1.+zplus**2)/(1.+zminus**2))
      sstar=emav**2
      prsq=(s-(emn-emav)**2)*(s-(emn+emav)**2)/(4.*s)
      pr=sqrt(prsq)
      qsqstar=(sstar-(emn-empi)**2)*(sstar-(emn+empi)**2)/(4.*sstar)
      qsqm0=(em0**2-(emn-empi)**2)*(em0**2-(emn+empi)**2)/(4.*em0**2)
      q0=sqrt(qsqm0)
      q=sqrt(qsqstar)
c
      sig10=(0.35e+6)*(1./(2.*psq))*alpha*(pr/p0)**beta
      sig10=sig10*(em0**2*gam**2*(q/q0)**3)/((sstar-em0**2)**2+em0**2*
     1gam**2)
      sig10=sig10*3.50
      alpha=3.772
      beta=1.262
      em0c=1220.
      em0=1188.
      gam=99.02
      g0=120.
      zplus=(sqrt(s)-emn-em0c)*(2./g0)
      zminus=(emn+empi-em0c)*(2./g0)
      emav=em0c+(1./(atan(zplus)-atan(zminus)))*
     1(g0/4.)*log((1.+zplus**2)/(1.+zminus**2))
      sstar=emav**2
      prsq=(s-(emn-emav)**2)*(s-(emn+emav)**2)/(4.*s)
      pr=sqrt(prsq)
      qsqstar=(sstar-(emn-empi)**2)*(sstar-(emn+empi)**2)/(4.*sstar)
      qsqm0=(em0**2-(emn-empi)**2)*(em0**2-(emn+empi)**2)/(4.*em0**2)
      q0=sqrt(qsqm0)
      q=sqrt(qsqstar)
c
      sig11=(0.35e+6)*(1./(2.*psq))*alpha*(pr/p0)**beta
      sig11=sig11*(em0**2*gam**2*(q/q0)**3)/((sstar-em0**2)**2+em0**2*
     1gam**2)
      sig11=sig11*3.66
c
c begin sig01
c
      alpha=146.3
      beta=0.
c     em0c=1220.
      em0c=1430.
      g0=200.
      em0=1472.
      gam=26.49
c     g0=120.
      zplus=(sqrt(s)-emn-em0c)*(2./g0)
      zminus=(emn+empi-em0c)*(2./g0)
      emav=em0c+(1./(atan(zplus)-atan(zminus)))*
     1(g0/4.)*log((1.+zplus**2)/(1.+zminus**2))
      sstar=emav**2
      prsq=(s-(emn-emav)**2)*(s-(emn+emav)**2)/(4.*s)
      pr=sqrt(prsq)
      qsqstar=(sstar-(emn-empi)**2)*(sstar-(emn+empi)**2)/(4.*sstar)
      qsqm0=(em0**2-(emn-empi)**2)*(em0**2-(emn+empi)**2)/(4.*em0**2)
      q0=sqrt(qsqm0)
      q=sqrt(qsqstar)
c
      sig01=(0.35e+6)*(1./(2.*psq))*alpha*(pr/p0)**beta
      sig01=sig01*(em0**2*gam**2*(q/q0)**3)/((sstar-em0**2)**2+em0**2*
     1gam**2)
      sig01=sig01*3.66
      continue
      ppinel(ie)=sigd+sig10+2.*sig11
      pninel(ie)=0.5*(sigd+sig10+2.*sig11+3.*sig01)
c     write(17,101)
c     write(17,102)e1,ppel,ppinel,pnel,pninel,sig10,sig11,sig01,sigd
c101  format(' E inc   Spp elas  Spp inel Spn elas    Spn inel     ',
c    1' S10      S11       S01         S10D')
c102  format(1x,f7.1,8(1x,1pe9.3))
      sppplus(ie)=sigd+sig10+sig11
      sppzero(ie)=sig11
      snnmins(ie)=sigd+sig11+sig10
      snnzero(ie)=sig11
      spnzero(ie)=0.5*(sigd+sig10+sig01)
      spnplus(ie)=0.5*(sig11+sig01)
      spnmins(ie)=spnplus(ie)
c     pxxpi(ie)=ppinel/ppel
c     pxypi(ie)=pninel/pnel
c     sum=sppplus(ie)+sppzero(ie)
c     if(sum.gt.0.)fppplus(ie)=sppplus(ie)/sum
c     if(sum.gt.0.)fnnmins(ie)=snnmins(ie)/sum
c     sum=spnzero(ie)+spnplus(ie)+spnmins(ie)
c     if(sum.eq.0.)go to 100
c     fpnplus(ie)=spnplus(ie)/sum
c     fpnzero(ie)=spnzero(ie)/sum+fpnplus(ie)
c
c     fpnmins(ie)=spnmins(ie)/sum
c     write(17,200)
c     write(17,201)e1,pxxpi(ie),pxypi(ie),fppplus(ie),fnnmins(ie),
c    1fpnplus(ie),fpnzero(ie)
c200  format(' E inc   pxxpi     pxypi     fppplus   fnnmins   ',
c    1'fpnplus   fpnzero ')
c201  format(1x,f7.1,6(1x,1pe9.3))
c
c    pxxpi gives probability a XX collision will give a pion
c    pxypi gives probability a XY collision(n,p) will give a pion,else elastic
c    fppplus gives probability that a pion formed in a pp collision will be pi+
c    then 1.-fppplus is probability the pion will be pi0
c    fnnmins is probability that a pion from nn collision is a pi-,else it is pi0
c    fpnplus is probability that a pn collision making a pion gives a pi+
c    fpnzero is probability a pn collision makes a pi0 added to the pi+ probability,
c    so that a number gt. fpnplus and less than fpnzero means pi0 production
 100  continue
      return
      end
 
      subroutine nucmfp
      common/outlim/limout,lim2
      common/eder/ed2,v,vhalf,vsq
      common/edep/pxx(2,1100)
      common/pl3/lz,nepr,eb(100),rcsp(100),z,a,plex
      common/hybsub/iav,idum,rs(8,1100),c,pox,r1
      common/incr/ed,tened,efinal
      common/rc/rzero
      common/angprm/aprp(500),eav3(1000),
     1eav2(1000),pot2x,thet1,phi1,snt(360),cst(360),snp(720),csp(720)
     2,pot
      common/pies/sppplus(1100),sppzero(1100),snnmins(1100),snnzero(1100
     1),spnzero(1100),spnplus(1100),spnmins(1100),pxxpi(1100),pxypi(1100
     2),fppplus(1100),fnnmins(1100),fpnplus(1100),fpnzero(1100)
     3,ppinel(1100),pninel(1100)
c
c     call pion
c
      w0=10.+0.345*(a-2.*z)
      e0=0.228*a
      do 2 ii=1,2
      do 3 kk=1,1100
      pxx(ii,kk)=0.
      rs(ii,kk)=0.
      rs(ii+2,kk)=0.
      rs(ii+4,kk)=0.
   3  continue
   2  continue
      an=a-z
      den=1./c
      dc=312.5*den
      volume=4.*(3.14/3.)*rzero**3*a
      dcn=den*3.*(10.**22)/volume
c
         kee=50./ed-ed/2.
      do 20 ke=1,1100
c
      ek=float(ke)*ed-ed/2.
      e1=pot+ek
      bta=sqrt(1.-(939.516/(939.516+e1))**2)
      aa=((10.63/bta**2-29.92/bta+42.9)+ppinel(ke))/a
      bb=((34.1/bta**2-82.2/bta+82.2)+pninel(ke))/a
      sig1=an*aa+z*bb
      sig2=z*aa+an*bb
      re=pot/e1
       if(re.gt.0.5)go to 10
      pa=1.-1.4*re
      go to 15
   10 pa=1.-1.4*re+.4*re*(2.-1./re)**2.5
   15 rs(1,ke)=dcn*sig1*pa*bta*a
      pxx(1,ke)=an*aa/sig1
      pxx(2,ke)=z*aa/sig2
      en=float(ke)*ed-ed/2.
      wa=(en/e0)*0.25*w0
      rs(3,ke)=wa*3.33*(10.**21)
      rs(4,ke)=wa*1.65*(10.**21)
      rs(5,ke)=wa*3.33*(10.**21)
      rs(6,ke)=wa*3.33*(10.**21)
   20 rs(2,ke)=dcn*sig2*pa*bta*a
      return
      end
      subroutine parap(mx,jcal)
c     quick calculation of nuclear cross sections via optical model
c     with parabolic barrier - see thomas, phys rev 116,703 (1959) -
c     program below by w. g. winn
      common/outlim/limout,lim2
      dimension cncs(999)
      common/parfs/lz,k6,delrr(999)
      common/par2/cncss
      common/par3/eq,sigml(999),acrs(999)
      common/shft/mc,mp,inver,ike,ipch,kq,qval,ap,at,zp,zt,cld,barfac
 
      en=eq
      do 10 i=1,999
      cncs(i)=0.
      acrs(i)=0.
   10 sigml(i)=0.
      ccs=0.
      lmax=200
      lcm=0
      write(17,20)
   20 format(' parabolic barrier approximation for reaction cross sectio
     1ns')
c     lmax=lmax+1
      coul=(1.4393)*zt*zp
       if(ap.le.4)go to 30
      continue
      ar=at**0.333333+ap**0.333333
      ra=1.22*ar
      f=67.0
      go to 35
   30 ar=at**.33333
       ar=at**.333333
      ra=1.17*ar
      f=1100.0
   35 continue
      d=0.574
      u=at*ap/(at+ap)
      h2=41.814
      acn=ap+at
      crotl=h2/(1.16*acn**1.666666)
      cr2=d*coul/f
      cr3p=h2*d/(f*u)
      if(lcm.gt.0)go to 45
      e=en*at/(at+ap)
      el=en
      go to 50
   45 e=en
      el=en*(at+ap)/at
   50  continue
      write(17,55)
   55 format(3h0  )
      write(17,60)
      write(17,65)
   60 format('projectile  target   projectile  target    energy    e
     1nergy ')
   65 format('   mass       mass      charge    charge     lab
     1 cm   ')
      write(17,70)
   70 format(3h   )
      write(17,75) ap,at,zp,zt,el,e
   75 format(4x,f7.3,5(3x,f7.3))
      write(17,55)
      write(17,80)
   80 format(' angular transmission barrier bar-radius  barrier    ba
     1rrier   relative curvature  probable sigma(l) rotor-e effrot effsi
     2gma   ')
      write(17,85)
   85 format('momenta coefficients  radii   constant    (mev)    cur
     1vature  sigma(l)   sign     l-energy  (mb)     (mev)   (mev)  (mb)
     2      ')
      write(17,70)
      write(17,90)
   90 format('    l         tl        rb        r0       vlb
     1wl       rsigl    sd2vlb     sigxl2   sigl     erotl  erotli sumdx
     2l   ')
      write(17,70)
      pwave=sqrt(2.*931.16*u*e)
      wavel=19.732/pwave
      arwave=3141.59*(wavel**2)
      sumal2=0.0
      sumtl=0.0
      sumdxl=0.0
      erotli=0.0
      do 200 k=1,lmax
      l=k-1
      r=ra
      al=float(l)
      cr3=cr3p*al*(al+1)
   95 rb=ra-d*alog(cr2/r**2+cr3/r**3)
      deltr=abs(rb-r)
      if(deltr.lt.0.001)go to 105
      r=rb
      go to 95
  105 continue
      r0=rb/ar
      vlb=coul/rb+(h2/2./u)*(al*(al+1)/rb**2)-f*exp(-(rb-ra)/d)
      d2vlb=2.*coul/rb**3+(h2/u)*3.*al*(al+1)/rb**4
     1 -(f/d**2)*exp(-(rb-ra)/d)
       if(d2vlb.eq.0.)go to 115
       if(d2vlb.gt.0.)go to 120
      sd2vlb=-1.0
      go to 125
  115 sd2vlb=0.0
      go to 125
  120 sd2vlb=1.0
  125 continue
      wl2=abs((h2/u)*d2vlb)
      wl=sqrt(wl2)
      tl=1./(1.+exp(2.*3.14159*(vlb-e)/wl))
      rsigl=tl*(2.*al+1.)
      sigl=arwave*rsigl
      sigml(k)=sigl
      if(jcal.le.0.or.jcal.eq.2)go to 130
      go to 135
  130 erotl=delrr(k)
      go to 140
  135 erotl=crotl*al*(al+1.)
      delrr(k)=erotl
  140 sigxl2=al*(al+1.)*rsigl
       if((erotl-0.5).gt.erotli)go to 155
      sumdxl=sumdxl+sigl
      if(tl.lt.0.0001)go to  155
      continue
      go to 175
  155 if(l.le.0)go to 170 
      continue
      write(17,165) erotli,sumdxl
      lx=erotl+.5
      if(tl.lt..0001)lx=lx+1
      acrs(k-1)=sumdxl
      mx=lx
  165 format(3h+  ,107x,f5.1,1x,f7.2)
  170 continue
      erotli=erotli+1.0
      sumdxl=sigl
  175 continue
      write(17,180) l,tl,rb,r0,vlb,wl,rsigl,sd2vlb,sigxl2,sigl,erotl
  180 format(3x,i3,4x,f10.8,2x,f6.3,4x,f6.3,4x,f6.1,4x,f6.3,4x,f6.2,4x,f
     16.3,1x,f10.2,2x,f7.3,3x,f7.3)
      if(l.gt.0.)go to 195 
      if(tl.ge.0.0001)go to 195
      continue
      write(17,165) erotli,sumdxl
  195 continue
      sumtl=sumtl+rsigl
      sumal2=sumal2+sigxl2
       tlx=tl-.0001
      if(tlx.lt.0.)go to 205
c     if(tl-.0001)205,200,200
  200 continue
  205 sigma=3141.59*(wavel**2)*sumtl
      write(17,210) sigma
  210 format(25h0cross section (mbarn) = ,f10.3)
      avel2=sumal2/sumtl
      avroe=avel2*crotl
      rmsl2=sqrt(avel2)
      write(17,215) rmsl2,avel2,avroe
  215 format(23h0root mean square l =  ,f6.2,10x,16h rmsl squared = ,f10
     1.3,10x,33h average rotation energy (mev) = ,f6.3)
      do 230 i=1,999
      ccs=ccs+sigml(i)
      cncs(i)=ccs
      if(sigml(i).gt.0.)go to 230
      jj=i-1
      go to 235
  230 jj=i
  235 continue
      write(17,240)cncs(jj)
      write(17,245)(cncs(i),i=1,jj)
      cncss=cncs(jj)
  240 format (38x,'total reaction cross section =',e10.3,' mb,'//
     1 '         cumulative partial reaction cross sections'/)
  245 format(9x,10e10.3)
      return
      end
      subroutine over
c     calculates the transmission coefficients
      common/eder/ed2,b,vhalf,vsq
      common/outlim/limout,lim2
      common/pl3/lz,nepr,eb(100),rcsp(100),zee,amass,plex
      common/sf/m3,kplt
      common/sft9/k5,jmax,pld,c(8),delt,alt
      common/paro/pq,cross
      common/lab3/sig(8,999),slav(8,999)
      common/par3/dq,tl(999),acrs(999)
      common /lab2/t(4,81),v(80),v1(3)
      common/incr/ed,tened,efinal
c  particle type is index 'k'; 
c  1= neutron
c  2=proton
c  3=alpha
c  4=deuteron
c  july 06 add 5=triton and
c  6= 3He
c  7= 7Be
c optical model parms from deuts were used modifeied only for mass and spin
c
c     m3=7
       ed2=ed/2.
      eq=pq
      jmin=1
      m2=1
      am=amass**.333333
c     ii=0.14*zee/ed
c     iik=0.28*zee/ed
c 7/17/2012--------------
       ii=0.07*zee/ed
        iik=0.14*zee/ed
c -------------
      ik=ii/2
      if(eq.gt.0.)m2=m3
      if(eq.gt.0.)jmin=jmax
      if(eq.gt.0.)go to 5
      delt=0.0
      alt=0.0
    5 do 150 k=m2,m3
      tsum=0.
      go to(10,15,20,25,26,26,26),k
   10 xjp=0.50
      xmp=1.00
      xmt=amass-1.0
      zp=0.0
      zt=zee
      p=48.00
      rv=1.322-.00076*amass+4.0*amass*amass*1.0e-06
     1-8.0*(amass**3)*1.0e-09
      av=.660
      w=9.00
      rw=1.266-.00037*amass+2.0*amass*amass*1.0e-06
     1-4.0*(amass**3)*1.0e-09
      rv=1.151+1.77*(amass-2.*zee)/(amass**1.3333)
      rw=rv
      aw=0.64
      vso=7.00
      rso=rw
      aso=aw
      rclmb=.001
      s=1.00
      xlm=80.0
      go to 30
   15 xjp=0.50
      rso=1.25
      xmp=1.00
      xmt=amass-1.0
      zp=1.0
      zt=zee-1.0
      vso=7.50
      aso=.51
      rclmb=1.25
      s=0.00
      xlm=80.0
      symm=(amass-2.*zee)/amass
      delpd=27.*(amass-2.*zee)/amass+.4*zee/am
      rv=1.2
      rw=1.55
      aw=.5
      av=.6
      p=60.
      w=5
      go to 30
   20 xjp=0.0
      xmp=4.0
      xmt=amass-4.0
      zp=2.00
      zt=zee-2.0
      p=50.2
c     rv=1.2+1.5/am
      rv=1.3+1.5/am
      av=.564
      w=12.30
      rw=rv
      aw=.564
      rso=1.
      vso=.001
      aso=1.0
      rclmb=1.30
      s=0.0
      xlm=80.0
      go to 30
   25 xjp=1.00
      xmp=2.00
      xmt=amass-2.
      zp=1.
      zt=zee-1.
      p=79.+2.*zee/am
      rv=1.15
      av=.81
      w=16.
      rw=1.01+1.26/am
      aw=.68
      rso=.98
      vso=5.6
      aso=1.00
      rclmb=1.30
      s=1.00
      xlm=80.
      go to 30
   26 if(k.eq.5)then
      xjp=0.5
      xmp=3.00
      xmt=amass-3.
      zp=1.
      zt=zee-1.
c     go to 27
      endif
      if(k.eq.5)go to 27
      if(k.eq.6)then
      xjp=0.5
      xmp=3.00
      xmt=amass-3.
      zp=2.
      zt=zee-2.
      p=50.2
      rv=1.2+1.5/am
      av=0.564
      w=12.3
      rw=rv
      aw=0.564
      rso=1.
      vso=0.001
      aso=1.0
      xlm=80.
      go to 30
c     go to 27
      endif
c     if(k.eq.6)go to 27
      if(k.eq.7)then
      xjp=1.5
      xmp=7.00
      xmt=amass-7.
      zp=4.
      zt=zee-4.
c     go to 27
      endif
      
   27   p=79.+2.*zee/am
      rv=1.15
      av=.81
      w=16.
      rw=1.01+1.26/am
      aw=.68
      rso=.98
      vso=5.6
      if(k.eq.6)then
       rv=1.2+1.5/am
       rw=rv
       w=12.3
       aw=0.564
       endif
c     aw=.68
c     rso=.98
c     vso=5.6
      aso=1.00
      rclmb=1.30
      s=1.00
      xlm=80.
   30 xf=xmt**.3333
      xm=xmt/(xmt + xmp)
      xjt=0.0
      lmx=xlm
      do 90 je=jmin,jmax
       engy=float(je)*ed-ed2
      if(eq.gt.0.)engy=eq
      e=engy
      if(k.eq.2.and.je.le.ik)go to 65
      if(k.eq.3.and.je.le.ii)go to 65
      if(k.eq.4.and.je.le.ik)go to 65
      if(k.eq.5.and.je.le.ik)go to 65
      if(k.eq.6.and.je.le.ii)go to 65
      if(k.eq.7.and.je.le.iik)go to 65
      v(1)=rv*xf
      v(2)=av
      v(3)=rw*xf
      v(4)=aw
      v(5)=s
      v(6)=vso
      v(9)=rclmb*xf
      v(10)=rso*xf
      v(11)=aso
      v1(1)=xjp
      v1(2)=+0.0
      w1=0.04783258*xmp*xm
      stplth=0.1
      v1(3)=xlm
      v(7)=p
      v(8)=w
      h4=0.04783*xmp*xm*e
      h1=sqrt(h4)
      z=(0.03478*zp*zt*xmp*xm)/h1
      c6=31.42/h4
      do 40 j=1,4
      do 40 i=1,81
   40 t(j,i)=0.0
      call tlj(h4,n1,j1,h1,z,w1,stplth)
      tsum=0.0
      tlsum=0.0
      tlsum2=0.
c feb 2002 add tlsum to get average tl
      intrpo=2.*xjp
      if(intrpo.lt.1)go to 55
      do 50 i=1,lmx
      tk=i-1
      tl(i) = c6* (t(1,i)*tk +t(2,i)*(tk+1.))
      tlsum=tlsum+tl(i)*(tk+1.)
       tlsum2=tlsum2+c6*(t(1,i)+t(2,i))*(tk+1.)
   50 tsum=tsum+tl(i)
      go to 65
   55 do 60 i=1,lmx
      tlsum=tlsum+t(1,i)*c6*(2*i+1)*(i)
      tlsum2=tlsum2+t(1,i)*c6*(i)
   60 tsum = tsum+t(1,i)*c6*(2*i-1)
   65 continue
      if(tsum.gt.0.01) go to 85
      go to(85,70,75,80),k
   70 delt=delt+1.0
      tsum =0.0
      go to 85
   75 alt=alt+1.0
      tsum =0.0
      go to 85
   80 tsum=0.
   85 if(eq.eq.0.)then
      sig(k,je)=tsum
      if(tsum.gt.0.)then
      slav(k,je)=tlsum/tlsum2-1.
      if(slav(k,je).lt.0.)slav(k,je)=0.
      else
      slav(k,je)=0.
      endif
      endif
      if(eq.gt.0.)cross=tsum
      if(eq.gt.0.)jj=i
   90 continue
      ajmax=float(jmax)*ed-ed2
      if(eq.gt.0.)go to 160
      go to (95,105,115,125,915,815,715),k
   95 write(17,100)
  100 format(/,38x,'optical model neutron inverse cross sections'/)
      go to 135
  105 write(17,110)
  110 format(//38x,'optical model proton inverse cross sections'/)
      go to 135
  115 write(17,120)
  120 format(//39x,'optical model alpha inverse cross sections'/)
      go to 135
  125 write(17,130)
  130 format(//37x,'optical model deuteron inverse cross sections'/)
      go to 135
  915 write(17,920)
  920 format(//37x,'optical model triton inverse cross sections'/)
      go to 135
  815 write(17,820)
  820 format(//37x,'optical model 3He inverse cross sections'/)
      go to 135
  715 write(17,720)
  720 format(//37x,'optical model 7Be inverse cross sections'/)
  135 write(17,140)ajmax
  140 format(' inverse reaction cross sect. for e=1 to',f8.3,'mev')
      write(17,145) (sig(k,je),je=1,jmax)
      write(17,145) (slav(k,je),je=1,jmax)
  145 format ( 1h ,10f9.3)
  150 continue
      return
  160 continue
      write(17,165)cross
  165 format(' reaction cross section computed by optical model subrouti
     1ne',f10.2,' mb')
      eq=0.
      write(17,170)
  170 format(' partial reaction cross sections,l=0 to lmax ')
      write(17,145)(tl(i),i=1,80)
      return
      end
      subroutine sigf
      common/eder/ed2,v,vhalf,vsq
      common/fizinv/sif(100,4,999)
      common/outlim/limout,lim2
      common/fis/ifis
      common/nhy/ij,jl,ji,jj,td,ex1,ex2,tmx,av,gav,cost,b(3)
      common/sft9/k5,jmax,pld,c(8),delt,alt
      common/ug/iz,ia
      common/pl3/lz,nepr,eb(100),rcsp(100),zee,amass,plex
      common/lab3/sig(8,999),slav(8,999)
      common/shft4/k3,jcal
      common/shft5/fiss
      common/shft2/fs(25),dsp(25),brr(25),der(25),er(25)
     1,cx(25),bilsum,xmiss,sumiz
      REAL*8 fiss
      common/incr/ed,tened,efinal
c     common/avergl/ind,cbb,ei,i,ik,sigtemp
c      common/avergl/cbb,etem,sigtemp,ind,i,ie
c      common/avergl/cbb,etem,sigtemp,ind,i,ie
      common/avergl/cbb,ei,sigtemp,ind,i,ik
      common/lab3f/slavf(100,4,999)
 
      ed2=ed/2.
      icall=1
      do 6 j=1,100
      do 6 k=1,4
      do 6 i=1,999
      slavf(j,k,i)=0.
  6   sif(j,k,i)=0.
c     if(jz.eq.0)jz=1
c1
      do jz=5,100
c
      j=jz
      aeff=float(jz)*2.4
      zeff=float(jz)
      rp=1.21*((aeff-1.)**.3333+1.)
      ra=1.21*((aeff-4.)**.3333+1.587)
      rd=1.21*((aeff-2.)**.3333+1.260)
      conrp=31.42*rp*rp
      conra=31.42*ra*ra
      conrd=31.42*rd*rd
      vp=(zeff-1.)*1.15/(rp+1.6)
      va=(zeff-2.)*2.64/(ra+1.6)
      vd=(zeff-1.)*1.32/(rd+1.6)
      rmp=1.-1./aeff
      rma=2.-8./aeff
      rmd=3.-6./aeff
c
      do 10 ik=1,999
c
      ei=float(ik)*ed-ed2
      sif(jz,1,ik)=31.42*(rp+3.4/sqrt(ei+.5))**2
      cbb=0.
      i=1
 
      if(ei.lt.300)call crslav1(nen,icall,jlim,j,aeff)
      if(vd.ge.ei) go to 11
       i=4
      cbb=vd
      if(ei.lt.300)call crslav1(nen,icall,jlim,j,aeff)
      sif(jz,4,ik) = conrd*(1.-vd/ei)
   11 if(vp.ge.ei)go to 10
       i=2
       cbb=vp
      if(ei.lt.300)call crslav1(nen,icall,jlim,j,aeff)
      sif(jz,2,ik)=conrp*(1.-vp/ei)
      if(va.ge.ei) go to 10
       i=3
       cbb=va
      if(ei.lt.300)call crslav1(nen,icall,jlim,j,aeff)
      sif(jz,3,ik)=conra*(1.-va/ei)
10    continue
      if (jl.ne.ji) go to 160
      ajmax=100.*ed-ed2
      if(jz.gt.1)go to 160
c     if(limout.eq.0)then
c     do 150 k=1,4
c     go to (95,105,115,125),k
c  95 write(17,100)
c 100 format(/,38x,'sharp cutoff neutron inverse cross sections'/)
c     go to 135
c 105 write(17,110)
c 110 format(//38x,'sharp cutoff proton inverse cross sections'/)
c     go to 135
c 115 write(17,120)
c 120 format(//39x,'sharp cutoff alpha inverse cross sections'/)
c     go to 135
c 125 write(17,130)
c 130 format(//37x,'sharp cutoff deuteron inverse cross sections'/)
c 135 write(17,140)ajmax
c 140 format('0inverse reaction cross sect. for e = 1 to',f8.3,' mev')
c     write(17,145) (sig(k,je),je=1,100)
c     write(17,146) (slav(k,je),je=1,100)
c 145 format ( 1h ,10f7.0)
c 146 format ( 1h ,10f7.3)
c 150 continue
c     endif
  160 do 200 k=1,999
      ei=(float(k)*ed-ed2)*ed
      sif(jz,1,k)=sif(jz,1,k)*ei*rmp
      sif(jz,2,k)=sif(jz,2,k)*ei*rmp
      sif(jz,3,k)=sif(jz,3,k)*ei*rma
      sif(jz,4,k)=sif(jz,4,k)*ei*rmd
200   continue
      enddo
      k3=3
c     call shaft
      return
      end
      subroutine sigi
      common/eder/ed2,v,vhalf,vsq
      common/outlim/limout,lim2
      common/fis/ifis
      common/nhy/ij,jl,ji,jj,td,ex1,ex2,tmx,av,gav,cost,b(3)
      common/sft9/k5,jmax,pld,c(8),delt,alt
      common/ug/iz,ia
      common/pl3/lz,nepr,eb(100),rcsp(100),zee,amass,plex
      common/lab3/sig(8,999),slav(8,999)
      common/shft4/k3,jcal
      common/shft5/fiss
      common/shft2/fs(25),dsp(25),brr(25),der(25),er(25)
     1,cx(25),bilsum,xmiss,sumiz
      REAL*8 fiss
      common/incr/ed,tened,efinal
c     common/avergl/ind,cbb,ei,i,ik,sigtemp
c     common/avergl/cbb,ei,sigtemp,ind,i,k
       common/avergl/cbb,ei,sigtemp,ind,i,ik
      ed2=ed/2.
 
      icall=1
      do 6 i=1,999
      do 6 k=1,8
      slav(k,i)=0.
  6   sig(k,i)=0.
      if(iz.eq.0)iz=1
      az=iz
      aeff=amass+2.-az-az
      zeff=zee-az+1.
      rp=1.21*((aeff-1.)**.3333+1.)
      ra=1.21*((aeff-4.)**.3333+1.587)
      rd=1.21*((aeff-2.)**.3333+1.260)
      conrp=31.42*rp*rp
      conra=31.42*ra*ra
      conrd=31.42*rd*rd
      vp=(zeff-1.)*1.15/(rp+1.6)
      va=(zeff-2.)*2.64/(ra+1.6)
      vd=(zeff-1.)*1.32/(rd+1.6)
      rmp=1.-1./aeff
      rma=2.-8./aeff
      rmd=3.-6./aeff
c
      do 10 ik=1,999
c
      ei=float(ik)*ed-ed2
      sig(1,ik)=31.42*(rp+3.4/sqrt(ei+.5))**2
      cbb=0.
      i=1
 
      if(ei.lt.300)call crslav(nen,icall,jlim)
      if(vd.ge.ei) go to 11
       i=4
      cbb=vd
      if(ei.lt.300)call crslav(nen,icall,jlim)
      sig(4,ik) = conrd*(1.-vd/ei)
   11 if(vp.ge.ei)go to 10
       i=2
       cbb=vp
      if(ei.lt.300)call crslav(nen,icall,jlim)
      sig(2,ik)=conrp*(1.-vp/ei)
      if(va.ge.ei) go to 10
       i=3
       cbb=va
      if(ei.lt.300)call crslav(nen,icall,jlim)
      sig(3,ik)=conra*(1.-va/ei)
10    continue
c have algos only for n,p,d,4He- must add for t,3He,7Be
      if (jl.ne.ji) go to 160
      ajmax=100.*ed-ed2
      if(iz.gt.1)go to 160
      do 150 k=1,4
      go to (95,105,115,125),k
   95 write(17,100)
  100 format(/,38x,'sharp cutoff neutron inverse cross sections'/)
      go to 135
  105 write(17,110)
  110 format(//38x,'sharp cutoff proton inverse cross sections'/)
      go to 135
  115 write(17,120)
  120 format(//39x,'sharp cutoff alpha inverse cross sections'/)
      go to 135
  125 write(17,130)
  130 format(//37x,'sharp cutoff deuteron inverse cross sections'/)
  135 write(17,140)ajmax
  140 format('0inverse reaction cross sect. for e = 1 to',f8.3,' mev')
      write(17,145) (sig(k,je),je=1,100)
      write(17,146) (slav(k,je),je=1,100)
  145 format ( 1h ,10f7.0)
  146 format ( 1h ,10f7.3)
  150 continue
  160 do 200 k=1,999
      ei=(float(k)*ed-ed2)*ed
      sig(1,k)=sig(1,k)*ei*rmp
      sig(2,k)=sig(2,k)*ei*rmp
      sig(3,k)=sig(3,k)*ei*rma
      sig(4,k)=sig(4,k)*ei*rmd
200   continue
      k3=3
      call shaft
      return
      end
      subroutine lymasf
      common/shft/mc,mp,inver,ike,ipch,kq,qval,ap,at,zp,zt,cld,barfac
      common/outlim/limout,lim2
      dimension em(10),xk(10),y(2),f(2),           emp(10)
      common/frags/xmasf(160,106),bq(106,160,4),pf(106,160)
     1,ppf(106,160,100),ppfi(106,160),ppfz(106),ppfa(260)
     2,ppfaz(106,260)
      common/sft5/exc(18,27),xmxx
      common/lym1/be(18,27,8),pair(18,27),xmas(18,27),delshl(18,27),amas
     1(18,27)
      common/bym1/symb(18,27),symbp(18,27),symbs(18,27)
      common/fis/ifis
      common/sf/m3,kplt
      common/cshel/shel(18,27,2)
c     data blank,tab,for,rinp,abe/4h    ,4h tab,4h msl,4h inp,4h abe/
c 7/17/2012 change data blank--------------------
      character*4 blank,tab,for,rinp,abe
      character*4 symbpr,symbs,symb,symbp
      character*4 aaa,aab
      data blank,tab,for,rinp,abe/'    ',' tab',' msl',' inp',' abe'/
c-----------------------------------------------
      do kz=1,106
      do kn=1,160
      xmasf(kn,kz)=0.
      pf(kz,kn)=0.
      do ik=1,100
      ppf(kz,kn,ik)=0.
      enddo
      do in=1,4
      bq(kz,kn,in)=0.
      enddo
      enddo
      enddo
      ldopt=0
      if(ldopt.eq.6)mc=1
      if(ldopt.eq.6)mp=0
      del=0.
      if(mp.eq.0.)del=2.
      if(mp.eq.3.)del=1.
      do 1 i=1,106
      do 2 j=1,160
      xmasf(j,i)=0.
   2  continue
   1  continue
c
c      above statements define pairing treatment
c
      ibind=0
      if(ldopt.eq.6)ibind=1
      em(1)=0.0
      em(2)=2.0
      em(3)=8.0
      em(4)=14.0
      em(5)=28.0
      em(6)=50.0
      em(7)=82.0
      em(8)=126.0
      em(9)=184.0
      em(10)=258.0
      cay1=1.15303
      cay3=200.0
      cay4=11.0
      cay5=8.07144
      cay6=7.28899
      gamma=1.7826
      a1=15.4941
      a2=17.9439
      a3=0.7053
      d=0.444
      c=5.8
      smalc=0.325
      pval=0.
      do 15 i=1,10
      emp(i)=em(i)**(5.0/3.0)
   15 continue
      do 20 i=1,9
      xk(i)=0.6*(emp(i+1)-emp(i))/(em(i+1)-em(i))
   20 continue
      rz=.863987/a3
      l=0
      z=1.0
      do 215 iz=10,106
      if(iz.le.30)then
      nlo=iz/2
      nhi=1.5*float(iz)
      endif
      if(iz.le.30)go to 25
      if(iz.le.50)then
      nlo=iz
      nhi=2*iz
      endif
      if(iz.le.50)go to 25
      if(iz.le.80)then
      nlo=1.2*float(iz)
      nhi=1.7*float(iz)
      endif
      if(iz.le.80)go to 25
       if(iz.le.106)then
       nlo=1.37*float(iz)
       nhi=1.83*float(iz)
       endif
   25 continue
       if(nlo.lt.10)nlo=10
       if(nhi.gt.160)nhi=160
      do 215 n=nlo,nhi
      ia=n+iz
      z=iz
      un=n
      a=ia
      a3rt=a**(1.0/3.0)
      a2rt=sqrt(a)
      a3rt2=a3rt**2.0
      zsq=z**2.0
      sym=((un-z)/a)**2
      acor=1.0-gamma*sym
      parmas=cay5*un+cay6*z
      volnuc=-1.0*a1*acor*a
      sufnuc=a2*acor*a3rt2
      coulmb=a3*zsq/a3rt
      fuzsur=-1.0*cay1*zsq/a
      oddev=-1.0*(1.0+2.0*(n/2)-un+2.*(iz/2)-z)/sqrt(a)*cay4
      pairr=-oddev
      symbpr=for
      if(mp.eq.0)oddev=0.
      if(mp.eq.0)go to 11
      pairr=(2.0*(n/2)-un+2.0*(iz/2)-z+del)/sqrt(a)*cay4
      symbpr = for
   11 continue
      if(sym.gt.0.4)wterm=0.
      wterm=0.
      wotnuc=parmas+coulmb+fuzsur+oddev+wterm
      smass=wotnuc+volnuc+sufnuc
      continue
      c2=(sufnuc+wterm)/(a**(2.0/3.0))
      x=coulmb/(2.0*(sufnuc+wterm))
      barr=0.0
      y(1)=un
      y(2)=z
      do 165 j=1,2
      do 150 i=1,9
      if (y(j).le.em(i+1))go to 160 
  150 continue
      stop
  160 f(j)=xk(i)*(y(j)-em(i))-.6*(y(j)**(5./3.)-emp(i))
  165 continue
      s=(2.0/a)**(2.0/3.0)*(f(1)+f(2))-smalc*a**(1./3.)
      ee=2.*c2*d**2*(1.0-x)
      ff=.42591771*c2*d**3*(1.+2.*x)/a3rt
      sshell=c*s
      v=sshell/ee
      eps=1.5*ff/ee
      if(ee*(1.-3.*v).le.0.0) go to 170
      qcalc=0.0
      theta=0.0
      shll=sshell
      go to 210
  170 to=1.0
  175 do 180 ipq=1,10
      t=to-(1.-eps*to-v*(3.-2.*to**2)*exp(-to**2))/(-eps+v*(10.*to-4.
     1 *to**3)*exp(-to**2))
      if (t.le.0.0) go to 190
      if (abs(t-to) .lt.0.0001) go to 185
      to=t
  180 continue
      go to 200
  185 if (2.*ee*(1.-2.*eps*t-v*(3.-12.*t**2+4.*t**4)*exp(-t**2))
     1 .gt.0.0) go to 205
  190 do 195 i=1,20
      to=float(i)/10.
      gl=ee*(1.-eps*to-v*(3.-2.*to**2)*exp(-to**2))
      if (gl.ge.0.0) go to 175
  195 continue
  200 continue
      go to 215
  205 theta=t
      alpha0=d*sqrt(5.)/a**(1./3.)
      alpha=alpha0*theta
      sigma=alpha*(1.+alpha/14.)
      qcalc=.004*z*(rz*a3rt)**2*(exp(2.*sigma)-exp(-sigma))
      sheldef=sshell*(1.-2.*t**2)*exp(-t**2)
      shll=ee*t**2-ff*t**3+sshell*(1.-2.*t**2)*exp(-t**2)
  210 if(mc.ne.1.or.mp.ne.0) go to 211
c     pair(jz,ja)=pair(jz,ja)-shll
c      shll = 0.
c     if(ldopt.ne.6)shll=0.
  211 cmass=smass+shll
c     xms(ja,jz)=cmass
c     amas(jz,ja)=smass
c     xmas(jz,ja)=cmass
c     symbs(jz,ja)=for
c     delshl(jz,ja)=shll
c     shel(jz,ja,2)=sheldef
      xmasf(n,iz)=cmass
c
c     if(mp.eq.2) pair(jz,ja)=pair(jz,ja)-shll
c     xq(jz)=cmass
  215 continue
      continue
      call bindf
      return
      end
      subroutine lymass(zee,amass,nz,na,mc,mp,ap,at,zp,zt,qval,ldopt)
      common/outlim/limout,lim2
      dimension em(10),xk(10),y(2),f(2),xms(24,15),           emp(10)
      dimension xq(30)
      common/sft5/exc(18,27),xmxx
      common/lym1/be(18,27,8),pair(18,27),xmas(18,27),delshl(18,27),amas
     1(18,27)
      common/bym1/symb(18,27),symbp(18,27),symbs(18,27)
      common/fis/ifis
      common/sf/m3,kplt
      common/cshel/shel(18,27,2)
c     data blank,tab,for,rinp,abe/4h    ,4h tab,4h msl,4h inp,4h abe/
c 7/17/2012 change data blank--------------------
      character*4 blank,tab,for,rinp,abe
      character*4 symbpr,symbs,symb,symbp
      character*4 aaa,aab
      data blank,tab,for,rinp,abe/'    ',' tab',' msl',' inp',' abe'/
c-----------------------------------------------
      aset=0.
      zset=0.
      na=na+1
      nz=nz+1
      if(na.gt.21)then
      na=21
      aset=1.
                  endif
      if(nz.gt.11)then
      nz=11
      zset=1.
                  endif
      if(ldopt.eq.6)mc=1
      if(ldopt.eq.6)mp=0
c     data blank,tab,for,rinp,abe/4h    ,4h tab,4h msl,4h inp,4h abe/
c 12/23/2010 set  outer be values to large values
      do 1234 i=1,18
      do 1234 k=1,27
      do 1234 j=1,8
1234  be(i,k,j)=2000.
      do 1 i=1,18
      do 1 k=1,27
c 3/15/2011 enlarge loop indices to18,27
      symbp(i,k)=blank
      symb(i,k)=blank
      symbs(i,k)=blank
      xmas(i,k)=0.
      amas(i,k)=0.
      delshl(i,k)=0.
      pair(i,k)=0.
      if(i.gt.13)go to 1
      do 10088 l=1,8
      be(i,k,l)=0.
10088 continue
    1 continue
      del=0.
      if(mp.eq.0.)del=2.
      if(mp.eq.3.)del=1.
c
c      above statements define pairing treatment
c
      ibind=0
      if(ldopt.eq.6)ibind=1
      if(mc.lt.10) go to 6
      ibind = 1
      mc=mc-10
    6 em(1)=0.0
      em(2)=2.0
      em(3)=8.0
      em(4)=14.0
      em(5)=28.0
      em(6)=50.0
      em(7)=82.0
      em(8)=126.0
      em(9)=184.0
      em(10)=258.0
      cay1=1.15303
      cay3=200.0
      cay4=11.0
      cay5=8.07144
      cay6=7.28899
      gamma=1.7826
      a1=15.4941
      a2=17.9439
      a3=0.7053
      d=0.444
      c=5.8
      smalc=0.325
      pval=0.
      do 15 i=1,10
      emp(i)=em(i)**(5.0/3.0)
   15 continue
      do 20 i=1,9
      xk(i)=0.6*(emp(i+1)-emp(i))/(em(i+1)-em(i))
   20 continue
      rz=.863987/a3
      l=0
      z=1.0
   25 kz=zee
      ka=amass
      if(qval.eq.0..and.pval.eq.0.)go to 30
      go to 35
   30 nnz=3
      nna=1
      go to 90
   35 nnz=nz+2
      nna=na+2
      if(nna.gt.24)nna=24
      if(nnz.gt.15)nnz=15
 
      write (17,36)
   36 format (1h1,53x,'mass options'/)
c add new ldopt=6 option
      if(ldopt.eq.6)write(17,9986)
9986  format(' liquid drop shell plus pairing correction substituted'
     1 ,'for pairing correction')
      if (mc.eq.1.and.mp.ne.0) mc=0
      if(ibind.eq.1)write(17,86)
   86 format(20x,' experimental masses are used where tabul'
     1,'ated; liquid drop values otherwise.')
      if(mc.eq.0) write(17,70)
      if(mc.eq.1) write(17,65)
      if (mc.eq.2) write (6,70)
      if(mp.eq.0)write(17,75)
       if(mp.eq.1)write(17,80)
      if(mp.eq.3)write(17,81)
   81 format(' normal pairing shift with odd-even reference point')
   65 format(40x,'liquid drop without shell correction term')
   70 format(41x,'liquid drop with shell correction term')
   75 format(21x,'without pairing, i.e. level density pairing ',
     1 'shift absorbed in binding energies')
   80 format(10x,'with pairing, level density pairing shift',
     1 ' calc. from msl formula and applied in backshifted',
     2 ' convention')
      if (mp.eq.2) write(17,85)
   85 format (20x,'msl shell correction term included in level ',
     1 'density ground state shift')
   90 do 215 jz=1,nnz
      do 215 ja=1,nna
      if(qval.eq.0..and.pval.eq.0..and.ap.eq.0.)go to 215
      if(qval.eq.0..and.pval.eq.0.)go to 95
      go to 115
   95 if(jz.lt.2)go to 100
      if(jz.eq.2)go to 105
      go to 110
  100 z=zee
      a=amass
      ia=a
      iz=z
      n=ia-iz
      un=amass-zee
      go to 125
  105 z=zt
      a=at
      ia=a
      iz=z
      n=ia-iz
      un=at-zt
      go to 125
  110 z=zp
      a=ap
      iz=z
      ia=a
      n=ia-iz
      un=ap-zp
      go to 125
  115 ia=ka+2-ja-jz
      iz=kz+1-jz
      n=ia-iz
      z=iz
      un=n
      a=ia
  125 a3rt=a**(1.0/3.0)
      a2rt=sqrt(a)
      a3rt2=a3rt**2.0
      zsq=z**2.0
      sym=((un-z)/a)**2
      acor=1.0-gamma*sym
      parmas=cay5*un+cay6*z
      volnuc=-1.0*a1*acor*a
      sufnuc=a2*acor*a3rt2
      coulmb=a3*zsq/a3rt
      fuzsur=-1.0*cay1*zsq/a
      oddev=-1.0*(1.0+2.0*(n/2)-un+2.*(iz/2)-z)/sqrt(a)*cay4
      pair(jz,ja)=-oddev
      symbp(jz,ja)=for
      if(mp.eq.0)oddev=0.
      if(mp.eq.0)go to 11
      pair(jz,ja)=(2.0*(n/2)-un+2.0*(iz/2)-z+del)/sqrt(a)*cay4
      symbp(jz,ja) = for
   11 continue
      if(sym.gt.0.4)wterm=0.
      wterm=0.
      wotnuc=parmas+coulmb+fuzsur+oddev+wterm
      smass=wotnuc+volnuc+sufnuc
      xms(ja,jz)=smass
      xq(jz)=smass
      continue
      c2=(sufnuc+wterm)/(a**(2.0/3.0))
      x=coulmb/(2.0*(sufnuc+wterm))
      barr=0.0
      y(1)=un
      y(2)=z
      do 165 j=1,2
      do 150 i=1,9
      if (y(j).le.em(i+1))go to 160 
  150 continue
      stop
  160 f(j)=xk(i)*(y(j)-em(i))-.6*(y(j)**(5./3.)-emp(i))
  165 continue
      s=(2.0/a)**(2.0/3.0)*(f(1)+f(2))-smalc*a**(1./3.)
      ee=2.*c2*d**2*(1.0-x)
      ff=.42591771*c2*d**3*(1.+2.*x)/a3rt
      sshell=c*s
      v=sshell/ee
      eps=1.5*ff/ee
      if(ee*(1.-3.*v).le.0.0) go to 170
      qcalc=0.0
      theta=0.0
      shll=sshell
      go to 210
  170 to=1.0
  175 do 180 ipq=1,10
      t=to-(1.-eps*to-v*(3.-2.*to**2)*exp(-to**2))/(-eps+v*(10.*to-4.
     1 *to**3)*exp(-to**2))
      if (t.le.0.0) go to 190
      if (abs(t-to) .lt.0.0001) go to 185
      to=t
  180 continue
      go to 200
  185 if (2.*ee*(1.-2.*eps*t-v*(3.-12.*t**2+4.*t**4)*exp(-t**2))
     1 .gt.0.0) go to 205
  190 do 195 i=1,20
      to=float(i)/10.
      gl=ee*(1.-eps*to-v*(3.-2.*to**2)*exp(-to**2))
      if (gl.ge.0.0) go to 175
  195 continue
  200 continue
      go to 215
  205 theta=t
      alpha0=d*sqrt(5.)/a**(1./3.)
      alpha=alpha0*theta
      sigma=alpha*(1.+alpha/14.)
      qcalc=.004*z*(rz*a3rt)**2*(exp(2.*sigma)-exp(-sigma))
      sheldef=sshell*(1.-2.*t**2)*exp(-t**2)
      shll=ee*t**2-ff*t**3+sshell*(1.-2.*t**2)*exp(-t**2)
  210 if(mc.ne.1.or.mp.ne.0) go to 211
      pair(jz,ja)=pair(jz,ja)-shll
c      shll = 0.
      if(ldopt.ne.6)shll=0.
  211 cmass=smass+shll
      xms(ja,jz)=cmass
      amas(jz,ja)=smass
       xmas(jz,ja)=cmass
      symbs(jz,ja)=for
      delshl(jz,ja)=shll
c     shel(jz,ja,2)=sheldef
c
      if(mp.eq.2) pair(jz,ja)=pair(jz,ja)-shll
      xq(jz)=cmass
  215 continue
       qtest=0.
       if(qval.eq.0..and.pval.eq.0.)qtest=1.
      if(qval.eq.0..and.pval.eq.0..and.ap.eq.0.)go to 231
      if(qval.eq.0..and.pval.eq.0.)go to 220
      go to 240
  220 if(zp.gt.20.)go to 230
      continue
      ize=zp+.1
      ine=ap-zp+.1
      call mass(ize,ine,xq(3),er,iret)
  230 pval=xq(2)+xq(3)-xq(1)
  231 continue
c sept'06 add temp
      if(qtest.gt.0.)qval=pval
c 7/17/2012 clarify experimental Q value used
      write(17,*)'-------------------------------------'
      write(17,*)' Experimental Qvalue will override M&S if available'
      write(17,*)'-------------------------------------'
c gloris comment
      if(ibind.eq.0.and.qval.eq.0.)qval=pval
      pval=0.0000001
      go to 25
       if(na.gt.21)na=21
       if(nz.gt.11)nz=11
  240 do 245 jz=1,nz
      do 245 ja=1,na
      symb(jz,ja) = for
      be(jz,ja,1)=8.07+xms(ja+1,jz)-xms(ja,jz)
      be(jz,ja,2)=7.29+xms(ja,jz+1)-xms(ja,jz)
      be(jz,ja,4)=13.1290+xms(ja+1,jz+1)-xms(ja,jz)
      be(jz,ja,3)=2.42+xms(ja+2,jz+2)-xms(ja,jz)
      be(jz,ja,5)=14.95+xms(ja+2,jz+1)-xms(ja,jz)
      be(jz,ja,6)=14.93+xms(ja+1,jz+2)-xms(ja,jz)
      be(jz,ja,7)=15.77+xms(ja+3,jz+4)-xms(ja,jz)
c     be(jz,ja,7)=15.77+xms(ja+3,jz+3)-xms(ja,jz)
c     be(jz,ja,8)=15.77+xms(ja+4,jz+3)-xms(ja,jz)
c here use index 7 for 6Li,8 for 7 Li- but these need check for consistent use
c 12/23/2010 put be=500. if product will have fewer than 2 neutrons
c or if proton number lt.1
      itestz=kz+1-jz
      itesta=ka+2-ja-jz
      itestn=itesta-itestz
      if(itestz.lt.8) be(jz,ja,7)=500.
      if(itestz.lt.8) be(jz,ja,8)=500.
      if(itestz.le.3) be(jz,ja,3)=500.
      if(itestz.le.3) be(jz,ja,6)=500.
      if(itestz.lt.2) be(jz,ja,4)=500.
      if(itestz.lt.2) be(jz,ja,5)=500.
      if(itestz.lt.2) be(jz,ja,2)=500.
      if(itestn.lt.6) be(jz,ja,7)=500.
      if(itestn.lt.8) be(jz,ja,8)=500.
      if(itestn.lt.3) be(jz,ja,3)=500.
      if(itestn.lt.2) be(jz,ja,6)=500.
      if(itestn.lt.3) be(jz,ja,5)=500.
      if(itestn.lt.2) be(jz,ja,4)=500.
      if(itestn.lt.2) be(jz,ja,1)=500.
      if(itestn.lt.2) be(jz,ja,2)=500.
       
      
      do kkk=1,7
      if(be(jz,ja,kkk).lt.-10.)be(jz,ja,kkk)=-10.
      enddo
  245 continue
      if(ldopt.ne.6)go to 8211
      amx=pair(1,1)
      do 9211 iz=1,nz+2
      do 9211 ia=1,na+2
      if(amx.lt.pair(iz,ia))amx=pair(iz,ia)
9211  continue
      do 9212 iz=1,nz+2
      do 9212 ia=1,na+2
      pair(iz,ia)=pair(iz,ia)-amx
9212  continue
8211     if(ldopt.eq.6)mc=10
      if(ldopt.eq.6)mp=1
      if(ibind.eq.1) call binden(zee,amass,nz,na,ap,at,zp,zt,
     1   qval,m3,mp,mc,ldopt)
      if(qval.eq.0.)qval=pval
      if(mc.ne.2) go to 260
      do 255 iz=1,nz
      do 255 ia=1,na
      read(16,270)be1,be2,be3,be4,pdel
      if(be1.eq.0.) go to 400
      be(iz,ia,1)=be1
      be(iz,ia,2)=be2
      be(iz,ia,3)=be3
      be(iz,ia,4)=be4
      symb(iz,ia)=rinp
  400 if(pdel.eq.0.) go to 255
      pair(iz,ia)=pdel
      symbp(iz,ia)=rinp
  255 continue
      if(limout.eq.0)then
      write (6,275)
      endif
  260 continue
      write(17,280)
      do 265 iz=1,nz
      mz=zee+1-iz
      do 265 ia=1,na
      ma=amass+2-ia-iz
      write(17,285)mz,ma,(be(iz,ia,k),k=1,4),symb(iz,ia),
     1 pair(iz,ia),symbp(iz,ia),delshl(iz,ia),symbs(iz,ia)
      if(mp.le.0)pair(iz,ia)=0.
      pair(iz,ia)=10.*pair(iz,ia)
  265 continue
      write(17,281)
      do 865 iz=1,nz
      mz=zee+1-iz
      do 865 ia=1,na
      ma=amass+2-ia-iz
      write(17,286)mz,ma,(be(iz,ia,k),k=5,7),symb(iz,ia),
     1 pair(iz,ia),symbp(iz,ia),delshl(iz,ia),symbs(iz,ia)
      if(mp.le.0)pair(iz,ia)=0.
  865 continue
      continue
      write (6,268)
      aaa=tab
      aab=abe
  268 format (1h1)
      if(aset.eq.0.)na=na-1
      if(zset.eq.0.)nz=nz-1
  270 format(5f10.5)
  275 format (20x,'some binding energies or level density ground ',
     1 'state shifts provided by user')
  280 format (/30x,'binding energies and level density ground ',
     1 'state shifts used'//22x,' iz = z-index of nucleus,',
     2 ' ia = a-index of nucleus in program isotope table'/32x,
     3 'msl = calculated by myers swiatecki lysekill mass formula'/
     4 32x,'tab = taken from 1971 mass table'/32x,'inp = provided',
     5 ' by user'/32x,'abe = absorbed in binding energy'//
     6 10x,'  iz  ia     neutron      proton      alpha     ',
     7 'deuteron    source    gs shift  source shell corr   source'/)
  281 format (/30x,'binding energies and level density ground ',
     1 'state shifts used'//22x,' iz = z-index of nucleus,',
     2 ' ia = a-index of nucleus in program isotope table'/32x,
     3 'msl = calculated by myers swiatecki lysekill mass formula'/
     4 32x,'tab = taken from 1971 mass table'/32x,'inp = provided',
     5 ' by user'/32x,'abe = absorbed in binding energy'//
     6 10x,'  iz  ia     triton     3He         7Be        ',
     7 ' source    gs shift  source shell corr   source'/)
  285 format(10x,2i4,2x,4(f10.5,2x),3x,a4,2(3x,f10.5,3x,a4))
  286 format(10x,2i4,2x,3(f10.5,2x),3x,a4,2(3x,f10.5,3x,a4))
      return
      end
c
      subroutine bindf
      common/frags/xmasf(160,106),bq(106,160,4),pf(106,160)
     1,ppf(106,160,100),ppfi(106,160),ppfz(106),ppfa(260)
     2,ppfaz(106,260)
c xmasf for fission calcs,dimensioned on N number(to 160) and Z (to 106)
      common/outlim/limout,lim2
      common/sft5/exc(18,27),xmxx
      common/lym1/be(18,27,8),pair(18,27),xmas(18,27),delshl(18,27),amas
     1(18,27)
      common/bym1/symb(18,27),symbp(18,27),symbs(18,27)
      dimension exces(4),res(4),err(4)
c     data blank,tab,for,rinp,abe/4h    ,4h tab,4h msl,4h inp,4h abe/
c 7/17/2012 change data blank--------------------
      character*4 blank,tab,for,rinp,abe
      character*4 symbpr,symbs,symb,symbp
      character*4 aaa,aab
      data blank,tab,for,rinp,abe/'    ',' tab',' msl',' inp',' abe'/
c-----------------------------------------------
      mp=1
      del=0.
      if(mp.eq.0.)del=2.
      if(mp.eq.3.)del=1.
      call mass(0,1,exces(1),err(1),iret)
      call mass(1,0,exces(2),err(2),iret)
      call mass(2,2,exces(3),err(3),iret)
      call mass(1,1,exces(4),err(4),iret)
      do 215 iz=10,106
      if(iz.le.30)then
      nlo=iz/2
      nhi=1.5*float(iz)
      endif
      if(iz.le.30)go to 25
      if(iz.le.50)then
      nlo=iz
      nhi=2*iz
      endif
      if(iz.le.50)go to 25
      if(iz.le.80)then
      nlo=1.2*float(iz)
      nhi=1.7*float(iz)
      endif
      if(iz.le.80)go to 25
       if(iz.le.106)then
       nlo=1.37*float(iz)
       nhi=1.83*float(iz)
       endif
   25 continue
       if(nlo.lt.10)nlo=10
       if(nhi.gt.160)nhi=160
      do 215 n=nlo,nhi
      ia=n+iz
c 115 ia=ka+2-ja-jz
c     iz=kz+1-jz
c     n=ia-iz
      z=iz
      un=n
      a=ia
       a=z+un
      pf(iz,n)=10.*(2.0*(n/2)-un+2.0*(iz/2)-z+del)*11./sqrt(a)
      imass=iz+n
      izz=iz
      inn=imass-izz
      izz1=izz-1
      inn1=inn-1
      inn2=inn-2
      izz2=izz-2
      iret=0
      call mass(izz,inn,xmasq,err(1),iret)
      call mass(izz,inn1,res(1),err(1),iret)
      call mass(izz1,inn,res(2),err(2),iret)
      call mass(izz2,inn2,res(3),err(3),iret)
      call mass(izz1,inn1,res(4),err(4),iret)
      xmasf(n,iz)=xmasq
      if(iret.gt.0) go to 54
      if(mp.ne.0)go to 55
      xmasq=xmasq+pf(iz,n)/10.
      res(1)=res(1)+pf(iz,n-1)/10.
      res(2)=res(2)+pf(iz-1,n)/10.
      res(3)=res(3)+pf(iz-2,n-2)/10.
      res(4)=res(4)+pf(iz-1,n-1)/10.
   54 continue
   55 do 60 i=1,3
60    bq(iz,n,i)=-xmasq+exces(i)+res(i)
      continue
 215  continue
      return
      end
      subroutine bind(zee,amass,ap,at,zp,zt,qval)
      common/outlim/limout,lim2
      common/sft5/exc(18,27),xmxx
      common/lym1/be(18,27,8),pair(18,27),xmas(18,27),delshl(18,27),amas
     1(18,27)
      common/bym1/symb(18,27),symbp(18,27),symbs(18,27)
      dimension exces(8),res(8),err(8)
c     data blank,tab,for,rinp,abe/4h    ,4h tab,4h msl,4h inp,4h abe/
c 7/17/2012 change data blank--------------------
      character*4 blank,tab,for,rinp,abe
      character*4 symbpr,symbs,symb,symbp
      character*4 aaa,aab
      data blank,tab,for,rinp,abe/'    ',' tab',' msl',' inp',' abe'/
c-----------------------------------------------
       qval=0.
      izee=zee+.1
      imass=amass+.1
      iap=ap+.1
      iat=at+.1
      izp=zp+.1
      izt=zt+.1
      inp=iap-izp
      inx=iat-izt
      iret=0
      n=imass-izee
      if(ap.eq.0..and.qval.eq.0.)go to 1
      call mass(izp,inp,exces(1),err(1),iret)
      call mass(izt,inx,exces(2),err(2),iret)
      call mass(izee,n,exces(3),err(3),iret)
    1 continue
      if(qval.eq.0..and.ap.eq.0.)qval=0.000001
      if(qval.eq.0.000001)go to 1001
      qval=-exces(3)+exces(2)+exces(1)
 1001 continue
       return
       end
      subroutine binden(zee,amass,nz,na,ap,at,zp,zt,qval,m3,mp,mc,ldopt)
      common/outlim/limout,lim2
      common/sft5/exc(18,27),xmxx
      common/lym1/be(18,27,8),pair(18,27),xmas(18,27),delshl(18,27),amas
     1(18,27)
      common/bym1/symb(18,27),symbp(18,27),symbs(18,27)
      dimension exces(8),res(8),err(8)
c     data blank,tab,for,rinp,abe/4h    ,4h tab,4h msl,4h inp,4h abe/
c 7/17/2012 change data blank--------------------
      character*4 blank,tab,for,rinp,abe
      character*4 symbpr,symbs,symb,symbp
      character*4 aaa,aab
      data blank,tab,for,rinp,abe/'    ',' tab',' msl',' inp',' abe'/
c-----------------------------------------------
      izee=zee+.1
      imass=amass+.1
      iap=ap+.1
      iat=at+.1
      izp=zp+.1
      izt=zt+.1
      inp=iap-izp
      inx=iat-izt
      iret=0
      n=imass-izee
      if(ap.eq.0..and.qval.eq.0.)go to 1
      call mass(izp,inp,exces(1),err(1),iret)
      call mass(izt,inx,exces(2),err(2),iret)
      call mass(izee,n,exces(3),err(3),iret)
 1    if(iret.gt.0) go to 10
      if(qval.ne.0.)go to 10
      if(qval.eq.0..and.ap.eq.0.)qval=0.000001
      if(qval.eq.0.000001)go to 1001
      qval=-exces(3)+exces(2)+exces(1)
 1001 continue
    4 format(//,20x,'********************************************'//)
      write(17,4)
      write(17,5)qval
      write(17,4)
      write(*,4)
      write(*,5)qval
      write(*,4)
5     format(20x,' qval used (calculated from mass table) = ',f5.1)
10    iret=0
      call mass(0,1,exces(1),err(1),iret)
      call mass(1,0,exces(2),err(2),iret)
      call mass(2,2,exces(3),err(3),iret)
      call mass(1,1,exces(4),err(4),iret)
      call mass(1,2,exces(5),err(5),iret)
      call mass(2,1,exces(6),err(6),iret)
      call mass(4,3,exces(7),err(7),iret)
      if(iret.gt.0) return
      write(17,6)
6     format(30x,'binding energies calculated from mass table',
     1 ' where possible.')
      if(mp.eq.0)write(17,7)
7     format(35x,'pairing removed from masses in table')
      if (mp.eq.0.and.mc.eq.1) write(17,8)
    8 format (35x,'shell correction removed from masses in table')
      if(na.gt.21)na=21
      if(nz.gt.11)nz=11
      do 65 iz=1,nz
      do 65 ia=1,na
      izz=izee-iz+1
      inn=imass-izee-ia+1
      izz1=izz-1
      inn1=inn-1
      inn2=inn-2
      izz2=izz-2
      izz4=izz-4
      inn3=inn-3
      iret=0
      call mass(izz,inn,xmass,err(1),iret)
      call mass(izz,inn1,res(1),err(1),iret)
c     if(m3.ge.2)call mass(izz1,inn,res(2),err(2),iret)
c     if(m3.ge.3)call mass(izz2,inn2,res(3),err(3),iret)
c     if(m3.eq.4)call mass(izz1,inn1,res(4),err(4),iret)
      call mass(izz1,inn,res(2),err(2),iret)
      call mass(izz2,inn2,res(3),err(3),iret)
      call mass(izz1,inn1,res(4),err(4),iret)
      call mass(izz1,inn2,res(5),err(5),iret)
      call mass(izz2,inn1,res(6),err(6),iret)
      call mass(izz4,inn3,res(7),err(7),iret)
      xmas(iz,ia)=xmass
      if(iret.gt.0) go to 65
      if(mp.ne.0)go to 55
      xmass=xmass+pair(iz,ia)
      res(1)=res(1)+pair(iz,ia+1)
      res(2)=res(2)+pair(iz+1,ia)
      res(3)=res(3)+pair(iz+2,ia+2)
      res(4)=res(4)+pair(iz+1,ia+1)
      res(5)=res(5)+pair(iz+1,ia+2)
      res(6)=res(6)+pair(iz+2,ia+1)
      res(7)=res(7)+pair(iz+4,ia+3)
   55 do 60 i=1,7
60    be(iz,ia,i)=-xmass+exces(i)+res(i)
cccc      delshl(iz,ia)=xmas(iz,ia)-amas(iz,ia)
      symb(iz,ia)=tab
65    continue
      aaa=abe
      aab=for
      return
      end
      subroutine fisrot(a,z,an,al,delr,delsp,ero,barfac)
      common/outlim/limout,lim2
      dimension x1b(6,11),x2b(6,11),x3b(10,20),x1h(6,11),x2h(6,11),x3h(1
     10,20)
      common/pl3/lz,nepr,eb(100),rcsp(100),zee,amass,plex
      common/parfs/mz,k6,delrr(999)
      data x1b/.28,.243,.221,.208,.195,.18,.211,.186,.17,.1506,.136,.12,
     1.152,.131,.1155,.096,.0795,.0625,.09725,.0795,.065,.0506,.0375,.02
     153,.05771,.0455,.03414,.0235,.014,.0065,.03325,.0235,.0153,.0081,.
     1001,.0,.01625,.009,.0032,.0,.0,.0,.0071,.0,.0,.0,.0,.0,.0,.0,.0,.0
     1,.0,.0,.0,.0,.0,.0,.0,.0,0.,0.,0.,0.,0.,0./
      data x1h/.0,.0,.0,.0,.0,.0,-.0057,-.0058,-.006,-.0061,-.0062,-.006
     13,-.0193,-.0203,-.0211,-.022,-.023,-.0245,-.0402,-.0427,-.0456,-.0
     1497,-.054,-.0616,-.0755,-.0812,-.0899,-.0988,-.109,-.12,-.1273,-.1
     1356,-.147,-.1592,-.1745,-.1897,-.1755,-.1986,-.2128,-.2296,-.251,-
     1.26,-.255,-.271,-.291,-.301,-.327,-.335,-.354,-.36,-.365,-.372,-.4
     103,-.42,-.35,-.35,-.35,-.35,-.35,-.35,-.35,-.35,-.35,-.35,-.35,
     1-.35/
      data x2b/.18,.1695,.1515,.133,.1155,.0949,.1495,.1363,.1165,.099,.
     10815,.0594,.12,.1032,.0864,.0678,.0469,.028,.09,.0725,.0556,.037,.
     1019,.0057,.0625,.045,.0304,.016,.005,0.,.0406,.0264,.0151,.0052,0.
     1,0.,.0253,.0144,.0027,0.,0.,0.,.0141,.006,0.,0.,0.,0.,.0065,.0008,
     10.,0.,0.,0.,.002,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0./
      data x2h/0.,0.,0.,0.,0.,0.,-.0018,-.0019,-.00215,-.0024,-.0025,-.0
     103,-.0063,-.00705,-.0076,-.0083,-.0091,-.0095,-.015,-.0158,-.0166,
     1-.0192,-.0217,-.025,-.0245,-.0254,-.029,-.0351,-.0478,-.0613,-.038
     1 7,-.0438,-.0532,-.0622,-.0845,-.0962,-.0616,-.0717,-.0821,-.0972,
     1-.1123,-.1274,-.0793,-.1014,-.1138,-.1262,-.1394,-.1526,-.12,-.134
     1,-.1503,-.1666,-.1829,-.1992,-.1528,-.171,-.1907,-.2104,-.2301,-.2
     1498,-.23,-.23,-.23,-.23,-.23,-.23/
      data x3h/0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,-.00012,-.00014,-.00016,-.0
     10018,-.0002,-.00024,-.00029,-.00036,-.00065,-.00089,-.00047,-.0005
     1,-.00058,-.00065,-.00074,-.00085,-.00101,-.00124,-.00138,-.00178,-
     1.001,-.00105,-.00124,-.00138,-.00156,-.00179,-.00275,-.00292,-.003
     1,-.003,-.00176,-.0019,-.00211,-.00235,-.00263,-.00298,-.00449,-.00
     153,-.0053,-.0053,-.003,-.00308,-.00318,-.00352,-.00392,-.00417,-.0
     1062,-.0062,-.0062,-.0062,-.00374,-.0041,-.00444,-.00488,-.00521,-.
     100545,-.0066,-.0066,-.0066,-.0066,-.0053,-.0055,-.00585,-.0064,-.0
     10695,-.007,-.007,-.007,-.007,-.007,-.00632,-.007,-.00742,-.00792,-
     1.00856,-.009,-.009,-.009,-.009,-.009,-.0079,-.0085,-.01022,-.0119,
     1-.012,-.012,-.012,-.012,-.012,-.012,-.00944,-.0102,-.0142,-.0182,-
     1.019,-.019,-.019,-.019,-.019,-.019,-.0112,-.0133,-.0182,-.0238,-.0
     124,-.024,-.024,-.024,-.024,-.024,-.01303,-.0178,-.0226,-.0274,-.02
     18,-.028,-.028,-.028,-.028,-.028,-.0165,-.0254,-.0343,-.0343,-.034,
     1-.034,-.034,-.034,-.034,-.034,-.0203,-.033,-.04,-.04,-.04,-.04,-.0
     14,-.04,-.04,-.04,-.025,-.0406,-.046,-.047,-.047,-.047,-.047,-.047,
     1-.047,-.047,-.03036,-.0482,-.048,-.048,-.048,-.048,-.048,-.048,-.0
     148,-.048,-.0363,-.0558,-.056,-.056,-.056,-.056,-.056,-.056,-.056,-
     1.056,-.04234,-.0634,-.064,-.064,-.064,-.064,-.064,-.064,-.064,-.06
     14,-.064,-.064,-.064,-.064,-.064,-.064,-.064,-.064,-.064,-.064/
      data x3b/.0949,.0755,.0564,.0382,.0223,.0121,.00588,.00242,.00069,
     1.0001,.0873,.0684,.049,.0306,.0162,.0074,.00267,.00055,0.,0.,.0801
     1,.061,.0418,.0235,.0108,.00373,.00071,0.,0.,0.,.073,.054,.035,.017
     18,.0062,.00125,0.,0.,0.,0.,.0661,.047,.0284,.012,.0025,0.,0.,0.,0.
     1,0.,.0594,.0404,.022,.0065,0.,0.,0.,0.,0.,0.,.0528,.034,.0159,.002
     1,0.,0.,0.,0.,0.,0.,.0465,.0277,.01,0.,0.,0.,0.,0.,0.,0.,.0401,.021
     17,.0044,0.,0.,0.,0.,0.,0.,0.,.0339,.0158,.00024,0.,0.,0.,0.,0.,0.,
     10.,.028,.0106,0.,0.,0.,0.,0.,0.,0.,0.,.0219,.0064,0.,0.,0.,0.,0.,0
     1.,0.,0.,.0164,.0025,0.,0.,0.,0.,0.,0.,0.,0.,.0122,0.,0.,0.,0.,0.,0
     1.,0.,0.,0.,.0085,0.,0.,0.,0.,0.,0.,0.,0.,0.,.0057,0.,0.,0.,0.,0.,0
     1.,0.,0.,0.,.0035,0.,0.,0.,0.,0.,0.,0.,0.,0.,.0016,0.,0.,0.,0.,0.,0
     1.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0
     1.,0./
      if(k6.gt.0)a=amass
      if(k6.gt.0)z=zee
      if(k6.gt.0)an=a-z
      if(barfac.eq.0.)barfac=1.
      paren=1.-1.7826*((an-z)/a)**2
      eso=17.9439*paren*a**.666666
      x=0.019655*z*(z/a)/paren
      jm=1
      if(k6.gt.0)jm=300
      do 100 j=1,jm
      if(k6.gt.0)al=j-1
      ero=34.548*al*al/a**1.666666
      y=1.9254*al*al/(paren*a**2.3333333)
      ix=20.*x+1.
      cx=ix
      bx=20.*x+1.
      dx=bx-cx
      if(x.gt.0.25)go to 30
      by=10.*y+1.
      if(by.gt.9.)by=9.
      by=9.
      if(by.gt.1.)go to 25
      by=1.
   25 iy=by
      cy=iy
      dy=by-cy
      h1=(x1h(ix+1,iy)-x1h(ix,iy))*dx+x1h(ix,iy)
      h2=(x1h(ix+1,iy+1)-x1h(ix,iy+1))*dx+x1h(ix,iy+1)
      hf=(h2-h1)*dy+h1
      b2=(x1b(ix+1,iy+1 )-x1b(ix,iy+1))*dx+x1b(ix,iy+1)
      b1=(x1b(ix+1,iy)-x1b(ix,iy))*dx+x1b(ix,iy)
      bf=(b2-b1)*dy+b1
      go to 95
   30 if(x.gt.0.5)go to 60
      by=20.*y+1.
      if(by.gt.10.)by=10.
      if(by.le.1.)by=1.
      ix=ix-5
      iy=by
      cy=iy
      dy=by-cy
      h1=(x2h(ix+1,iy)-x2h(ix,iy))*dx+x2h(ix,iy)
      h2=(x2h(ix+1,iy+1)-x2h(ix,iy+1))*dx+x2h(ix,iy+1)
      hf=(h2-h1)*dy+h1
      b1=(x2b(ix+1,iy)-x2b(ix,iy))*dx+x2b(ix,iy)
      b2=(x2b(ix+1,iy+1 )-x2b(ix,iy+1))*dx+x2b(ix,iy+1)
      bf=(b2-b1)*dy+b1
      go to 95
   60 if(x.gt.0.95)x=0.95
      ix=20.*x+1.
      ix=ix-10
      by=100.*y+1.
      if(by.gt.19.)by=19.
      continue
      if(by.le.1.)by=1.
      continue
      iy=by
      cy=iy
      dy=by-cy
      h1=(x3h(ix+1,iy)-x3h(ix,iy))*dx+x3h(ix,iy)
      h2=(x3h(ix+1,iy+1)-x3h(ix,iy+1))*dx+x3h(ix,iy+1)
      hf=(h2-h1)*dy+h1
      b1=(x3b(ix+1,iy)-x3b(ix,iy))*dx+x3b(ix,iy)
      b2=(x3b(ix+1,iy+1 )-x3b(ix,iy+1))*dx+x3b(ix,iy+1)
      bf=(b2-b1)*dy+b1
   95 delr=ero+hf*eso
  100 if(k6.gt.0)delrr(j)=delr
      k6=0
      bf=bf*barfac
      delsp=delr+bf*eso
      return
      end
      subroutine mass(mz,n,exces,errr,iret)
      common/outlim/limout,lim2
c     using index(z*1000+n) finds mass excess(mev) and its error(kev)
c     from 1983 wapstra adjustments as presented in the ensdf
c     data base and rendered here by w.p. trower 10feb89.
      dimension itable(2999),amass(2999),err(2999),
     * i1 (100),i2 (100),i3 (100),i4 (100),i5 (100),i6 (100),i7 (100),
     * i8 (100),i9 (100),i10(100),i11(100),i12(100),i13(100),i14(100),
     * i15(100),i16(100),i17(100),i18(100),i19(100),i20(100),i21(100),
     * i22(100),i23(100),
     * a1 (100),a2 (100),a3 (100),a4 (100),a5 (100),a6 (100),a7 (100),
     * a8 (100),a9 (100),a10(100),a11(100),a12(100),a13(100),a14(100),
     * a15(100),a16(100),a17(100),a18(100),a19(100),a20(100),a21(100),
     * a22(100),a23(100),
     * e1 (100),e2 (100),e3 (100),e4 (100),e5 (100),e6 (100),e7 (100),
     * e8 (100),e9 (100),e10(100),e11(100),e12(100),e13(100),e14(100),
     * e15(100),e16(100),e17(100),e18(100),e19(100),e20(100),e21(100),
     * e22(100),e23(100)
      equivalence           (i1 (1),itable(   1)),(i2 (1),itable( 101)),
     *(i3 (1),itable(201 )),(i4 (1),itable( 301)),(i5 (1),itable( 401)),
     *(i6 (1),itable( 501)),(i7 (1),itable( 601)),(i8 (1),itable( 701)),
     *(i9 (1),itable( 801)),(i10(1),itable( 901)),(i11(1),itable(1001)),
     *(i12(1),itable(1101)),(i13(1),itable(1201)),(i14(1),itable(1301)),
     *(i15(1),itable(1401)),(i16(1),itable(1501)),(i17(1),itable(1601)),
     *(i18(1),itable(1701)),(i19(1),itable(1801)),(i20(1),itable(1901)),
     *(i21(1),itable(2001)),(i22(1),itable(2101)),(i23(1),itable(2201))
      equivalence          (a1 (1),amass(   1)),(a2 (1),amass( 101)),
     *(a3 (1),amass( 201)),(a4 (1),amass( 301)),(a5 (1),amass( 401)),
     *(a6 (1),amass( 501)),(a7 (1),amass( 601)),(a8 (1),amass( 701)),
     *(a9 (1),amass( 801)),(a10(1),amass( 901)),(a11(1),amass(1001)),
     *(a12(1),amass(1101)),(a13(1),amass(1201)),(a14(1),amass(1301)),
     *(a15(1),amass(1401)),(a16(1),amass(1501)),(a17(1),amass(1601)),
     *(a18(1),amass(1701)),(a19(1),amass(1801)),(a20(1),amass(1901)),
     *(a21(1),amass(2001)),(a22(1),amass(2101)),(a23(1),amass(2201))
      equivalence        (e1 (1),err(   1)),(e2 (1),err( 101)),
     *(e3 (1),err( 201)),(e4 (1),err( 301)),(e5 (1),err( 401)),
     *(e6 (1),err( 501)),(e7 (1),err( 601)),(e8 (1),err( 701)),
     *(e9 (1),err( 801)),(e10(1),err( 901)),(e11(1),err(1001)),
     *(e12(1),err(1101)),(e13(1),err(1201)),(e14(1),err(1301)),
     *(e15(1),err(1401)),(e16(1),err(1501)),(e17(1),err(1601)),
     *(e18(1),err(1701)),(e19(1),err(1801)),(e20(1),err(1901)),
     *(e21(1),err(2001)),(e22(1),err(2101)),(e23(1),err(2201))
      data nt/2213/
      data i1/
     *      -999,         0,       999,         1,         4,      1000,
     *      1001,      1002,      1003,      2001,      2002,      2003,
     *      2004,      2005,      2006,      2007,      3001,      3002,
     *      3003,      3004,      3005,      3006,      3007,      3008,
     *      4002,      4003,      4004,      4005,      4006,      4007,
     *      4008,      4009,      4010,      5002,      5003,      5004,
     *      5005,      5006,      5007,      5008,      5009,      5010,
     *      5011,      5012,      6002,      6003,      6004,      6005,
     *      6006,      6007,      6008,      6009,      6010,      6011,
     *      6012,      6013,      6014,      7003,      7004,      7005,
     *      7006,      7007,      7008,      7009,      7010,      7011,
     *      7012,      7013,      7014,      8004,      8005,      8006,
     *      8007,      8008,      8009,      8010,      8011,      8012,
     *      8013,      8014,      8015,      9005,      9006,      9007,
     *      9008,      9009,      9010,      9011,      9012,      9013,
     *      9014,      9015,      9016,     10006,     10007,     10008,
     *     10009,     10010,     10011,     10012/
      data i2/
     *     10013,     10014,     10015,     10016,     10017,     11007,
     *     11008,     11009,     11010,     11011,     11012,     11013,
     *     11014,     11015,     11016,     11017,     11018,     11019,
     *     11020,     11021,     11022,     11023,     12008,     12009,
     *     12010,     12011,     12012,     12013,     12014,     12015,
     *     12016,     12017,     12018,     12019,     12020,     12021,
     *     12022,     12023,     13009,     13010,     13011,     13012,
     *     13013,     13014,     13015,     13016,     13017,     13018,
     *     13019,     13020,     13021,     13022,     13023,     14010,
     *     14011,     14012,     14013,     14014,     14015,     14016,
     *     14017,     14018,     14019,     14020,     14021,     14022,
     *     14023,     14024,     15011,     15012,     15013,     15014,
     *     15015,     15016,     15017,     15018,     15019,     15020,
     *     15021,     15022,     15023,     15024,     15025,     16012,
     *     16013,     16014,     16015,     16016,     16017,     16018,
     *     16019,     16020,     16021,     16022,     16023,     16024,
     *     16025,     16026,     17013,     17014/
      data i3/
     *     17015,     17016,     17017,     17018,     17019,     17020,
     *     17021,     17022,     17023,     17024,     17025,     17026,
     *     17027,     18014,     18015,     18016,     18017,     18018,
     *     18019,     18020,     18021,     18022,     18023,     18024,
     *     18025,     18026,     18027,     18028,     19015,     19016,
     *     19017,     19018,     19019,     19020,     19021,     19022,
     *     19023,     19024,     19025,     19026,     19027,     19028,
     *     19029,     19030,     20016,     20017,     20018,     20019,
     *     20020,     20021,     20022,     20023,     20024,     20025,
     *     20026,     20027,     20028,     20029,     20030,     20031,
     *     21017,     21018,     21019,     21020,     21021,     21022,
     *     21023,     21024,     21025,     21026,     21027,     21028,
     *     21029,     21030,     21031,     22018,     22019,     22020,
     *     22021,     22022,     22023,     22024,     22025,     22026,
     *     22027,     22028,     22029,     22030,     22031,     22032,
     *     23019,     23020,     23021,     23022,     23023,     23024,
     *     23025,     23026,     23027,     23028/
      data i4/
     *     23029,     23030,     23031,     23032,     23033,     24020,
     *     24021,     24022,     24023,     24024,     24025,     24026,
     *     24027,     24028,     24029,     24030,     24031,     24032,
     *     24033,     24034,     25021,     25022,     25023,     25024,
     *     25025,     25026,     25027,     25028,     25029,     25030,
     *     25031,     25032,     25033,     25034,     25035,     26023,
     *     26024,     26025,     26026,     26027,     26028,     26029,
     *     26030,     26031,     26032,     26033,     26034,     26035,
     *     26036,     26037,     27024,     27025,     27026,     27027,
     *     27028,     27029,     27030,     27031,     27032,     27033,
     *     27034,     27035,     27036,     27037,     27038,     28025,
     *     28026,     28027,     28028,     28029,     28030,     28031,
     *     28032,     28033,     28034,     28035,     28036,     28037,
     *     28038,     28039,     28040,     28041,     29026,     29027,
     *     29028,     29029,     29030,     29031,     29032,     29033,
     *     29034,     29035,     29036,     29037,     29038,     29039,
     *     29040,     29041,     29042,     30027/
      data i5/
     *     30028,     30029,     30030,     30031,     30032,     30033,
     *     30034,     30035,     30036,     30037,     30038,     30039,
     *     30040,     30041,     30042,     30043,     30044,     30045,
     *     30046,     30047,     30048,     31030,     31031,     31032,
     *     31033,     31034,     31035,     31036,     31037,     31038,
     *     31039,     31040,     31041,     31042,     31043,     31044,
     *     31045,     31046,     31047,     31048,     31049,     31050,
     *     32031,     32032,     32033,     32034,     32035,     32036,
     *     32037,     32038,     32039,     32040,     32041,     32042,
     *     32043,     32044,     32045,     32046,     32047,     32048,
     *     32049,     32050,     32051,     33032,     33033,     33034,
     *     33035,     33036,     33037,     33038,     33039,     33040,
     *     33041,     33042,     33043,     33044,     33045,     33046,
     *     33047,     33048,     33049,     33050,     33051,     33052,
     *     34033,     34034,     34035,     34036,     34037,     34038,
     *     34039,     34040,     34041,     34042,     34043,     34044,
     *     34045,     34046,     34047,     34048/
      data i6/
     *     34049,     34050,     34051,     34052,     34053,     35034,
     *     35035,     35036,     35037,     35038,     35039,     35040,
     *     35041,     35042,     35043,     35044,     35045,     35046,
     *     35047,     35048,     35049,     35050,     35051,     35052,
     *     35053,     35054,     35055,     36035,     36036,     36037,
     *     36038,     36039,     36040,     36041,     36042,     36043,
     *     36044,     36045,     36046,     36047,     36048,     36049,
     *     36050,     36051,     36052,     36053,     36054,     36055,
     *     36056,     36057,     37036,     37037,     37038,     37039,
     *     37040,     37041,     37042,     37043,     37044,     37045,
     *     37046,     37047,     37048,     37049,     37050,     37051,
     *     37052,     37053,     37054,     37055,     37056,     37057,
     *     37058,     37059,     37060,     37061,     37062,     38039,
     *     38040,     38041,     38042,     38043,     38044,     38045,
     *     38046,     38047,     38048,     38049,     38050,     38051,
     *     38052,     38053,     38054,     38055,     38056,     38057,
     *     38058,     38059,     38060,     38061/
      data i7/
     *     38062,     39040,     39041,     39042,     39043,     39044,
     *     39045,     39046,     39047,     39048,     39049,     39050,
     *     39051,     39052,     39053,     39054,     39055,     39056,
     *     39057,     39058,     39059,     39060,     39061,     39062,
     *     40041,     40042,     40043,     40044,     40045,     40046,
     *     40047,     40048,     40049,     40050,     40051,     40052,
     *     40053,     40054,     40055,     40056,     40057,     40058,
     *     40059,     40060,     40061,     40062,     40063,     41044,
     *     41045,     41046,     41047,     41048,     41049,     41050,
     *     41051,     41052,     41053,     41054,     41055,     41056,
     *     41057,     41058,     41059,     41060,     41061,     41062,
     *     41063,     41064,     42045,     42046,     42047,     42048,
     *     42049,     42050,     42051,     42052,     42053,     42054,
     *     42055,     42056,     42057,     42058,     42059,     42060,
     *     42061,     42062,     42063,     42064,     42065,     43046,
     *     43047,     43048,     43049,     43050,     43051,     43052,
     *     43053,     43054,     43055,     43056/
      data i8/
     *     43057,     43058,     43059,     43060,     43061,     43062,
     *     43063,     43064,     43065,     43066,     44047,     44048,
     *     44049,     44050,     44051,     44052,     44053,     44054,
     *     44055,     44056,     44057,     44058,     44059,     44060,
     *     44061,     44062,     44063,     44064,     44065,     44066,
     *     44067,     45048,     45049,     45050,     45051,     45052,
     *     45053,     45054,     45055,     45056,     45057,     45058,
     *     45059,     45060,     45061,     45062,     45063,     45064,
     *     45065,     45066,     45067,     45068,     46049,     46050,
     *     46051,     46052,     46053,     46054,     46055,     46056,
     *     46057,     46058,     46059,     46060,     46061,     46062,
     *     46063,     46064,     46065,     46066,     46067,     46068,
     *     46069,     46070,     47050,     47051,     47052,     47053,
     *     47054,     47055,     47056,     47057,     47058,     47059,
     *     47060,     47061,     47062,     47063,     47064,     47065,
     *     47066,     47067,     47068,     47069,     47070,     47071,
     *     47072,     47073,     47074,     48051/
      data i9/
     *     48052,     48053,     48054,     48055,     48056,     48057,
     *     48058,     48059,     48060,     48061,     48062,     48063,
     *     48064,     48065,     48066,     48067,     48068,     48069,
     *     48070,     48071,     48072,     48073,     48074,     48075,
     *     49052,     49053,     49054,     49055,     49056,     49057,
     *     49058,     49059,     49060,     49061,     49062,     49063,
     *     49064,     49065,     49066,     49067,     49068,     49069,
     *     49070,     49071,     49072,     49073,     49074,     49075,
     *     49076,     49077,     49078,     49079,     49080,     49081,
     *     49082,     50052,     50053,     50054,     50055,     50056,
     *     50057,     50058,     50059,     50060,     50061,     50062,
     *     50063,     50064,     50065,     50066,     50067,     50068,
     *     50069,     50070,     50071,     50072,     50073,     50074,
     *     50075,     50076,     50077,     50078,     50079,     50080,
     *     50081,     50082,     50083,     51053,     51054,     51055,
     *     51056,     51057,     51058,     51059,     51060,     51061,
     *     51062,     51063,     51064,     51065/
      data i10/
     *     51066,     51067,     51068,     51069,     51070,     51071,
     *     51072,     51073,     51074,     51075,     51076,     51077,
     *     51078,     51079,     51080,     51081,     51082,     51083,
     *     51084,     52054,     52055,     52056,     52057,     52058,
     *     52059,     52060,     52061,     52062,     52063,     52064,
     *     52065,     52066,     52067,     52068,     52069,     52070,
     *     52071,     52072,     52073,     52074,     52075,     52076,
     *     52077,     52078,     52079,     52080,     52081,     52082,
     *     52083,     52084,     52085,     53055,     53056,     53057,
     *     53058,     53059,     53060,     53061,     53062,     53063,
     *     53064,     53065,     53066,     53067,     53068,     53069,
     *     53070,     53071,     53072,     53073,     53074,     53075,
     *     53076,     53077,     53078,     53079,     53080,     53081,
     *     53082,     53083,     53084,     53085,     53086,     54056,
     *     54057,     54058,     54059,     54060,     54061,     54062,
     *     54063,     54064,     54065,     54066,     54067,     54068,
     *     54069,     54070,     54071,     54072/
      data i11/
     *     54073,     54074,     54075,     54076,     54077,     54078,
     *     54079,     54080,     54081,     54082,     54083,     54084,
     *     54085,     54086,     54087,     54088,     55058,     55059,
     *     55060,     55061,     55062,     55063,     55064,     55065,
     *     55066,     55067,     55068,     55069,     55070,     55071,
     *     55072,     55073,     55074,     55075,     55076,     55077,
     *     55078,     55079,     55080,     55081,     55082,     55083,
     *     55084,     55085,     55086,     55087,     55088,     55089,
     *     55090,     55091,     55092,     55093,     56061,     56062,
     *     56063,     56064,     56065,     56066,     56067,     56068,
     *     56069,     56070,     56071,     56072,     56073,     56074,
     *     56075,     56076,     56077,     56078,     56079,     56080,
     *     56081,     56082,     56083,     56084,     56085,     56086,
     *     56087,     56088,     56089,     56090,     56091,     56092,
     *     56093,     57070,     57071,     57072,     57073,     57074,
     *     57075,     57076,     57077,     57078,     57079,     57080,
     *     57081,     57082,     57083,     57084/
      data i12/
     *     57085,     57086,     57087,     57088,     57089,     57090,
     *     57091,     57092,     57093,     58072,     58073,     58074,
     *     58075,     58076,     58077,     58078,     58079,     58080,
     *     58081,     58082,     58083,     58084,     58085,     58086,
     *     58087,     58088,     58089,     58090,     58091,     58092,
     *     58093,     59073,     59074,     59075,     59076,     59077,
     *     59078,     59079,     59080,     59081,     59082,     59083,
     *     59084,     59085,     59086,     59087,     59088,     59089,
     *     59090,     59091,     59092,     59093,     60074,     60075,
     *     60076,     60077,     60078,     60079,     60080,     60081,
     *     60082,     60083,     60084,     60085,     60086,     60087,
     *     60088,     60089,     60090,     60091,     60092,     60093,
     *     60094,     61075,     61076,     61077,     61078,     61079,
     *     61080,     61081,     61082,     61083,     61084,     61085,
     *     61086,     61087,     61088,     61089,     61090,     61091,
     *     61092,     61093,     61094,     61095,     62076,     62077,
     *     62078,     62079,     62080,     62081/
      data i13/
     *     62082,     62083,     62084,     62085,     62086,     62087,
     *     62088,     62089,     62090,     62091,     62092,     62093,
     *     62094,     62095,     62096,     63077,     63078,     63079,
     *     63080,     63081,     63082,     63083,     63084,     63085,
     *     63086,     63087,     63088,     63089,     63090,     63091,
     *     63092,     63093,     63094,     63095,     63096,     63097,
     *     64078,     64079,     64080,     64081,     64082,     64083,
     *     64084,     64085,     64086,     64087,     64088,     64089,
     *     64090,     64091,     64092,     64093,     64094,     64095,
     *     64096,     64097,     64098,     65079,     65080,     65081,
     *     65082,     65083,     65084,     65085,     65086,     65087,
     *     65088,     65089,     65090,     65091,     65092,     65093,
     *     65094,     65095,     65096,     65097,     65098,     65099,
     *     66080,     66081,     66082,     66083,     66084,     66085,
     *     66086,     66087,     66088,     66089,     66090,     66091,
     *     66092,     66093,     66094,     66095,     66096,     66097,
     *     66098,     66099,     66100,     66101/
      data i14/
     *     67081,     67082,     67083,     67084,     67085,     67086,
     *     67087,     67088,     67089,     67090,     67091,     67092,
     *     67093,     67094,     67095,     67096,     67097,     67098,
     *     67099,     67100,     67101,     67102,     67103,     68082,
     *     68083,     68084,     68085,     68086,     68087,     68088,
     *     68089,     68090,     68091,     68092,     68093,     68094,
     *     68095,     68096,     68097,     68098,     68099,     68100,
     *     68101,     68102,     68103,     68104,     68105,     69082,
     *     69083,     69084,     69085,     69086,     69087,     69088,
     *     69089,     69090,     69091,     69092,     69093,     69094,
     *     69095,     69096,     69097,     69098,     69099,     69100,
     *     69101,     69102,     69103,     69104,     69105,     69106,
     *     69107,     70082,     70083,     70084,     70085,     70086,
     *     70087,     70088,     70089,     70090,     70091,     70092,
     *     70093,     70094,     70095,     70096,     70097,     70098,
     *     70099,     70100,     70101,     70102,     70103,     70104,
     *     70105,     70106,     70107,     70108/
      data i15/
     *     71082,     71083,     71084,     71085,     71086,     71087,
     *     71088,     71089,     71090,     71091,     71092,     71093,
     *     71094,     71095,     71096,     71097,     71098,     71099,
     *     71100,     71101,     71102,     71103,     71104,     71105,
     *     71106,     71107,     71108,     71109,     72082,     72083,
     *     72084,     72085,     72086,     72087,     72088,     72089,
     *     72090,     72091,     72092,     72093,     72094,     72095,
     *     72096,     72097,     72098,     72099,     72100,     72101,
     *     72102,     72103,     72104,     72105,     72106,     72107,
     *     72108,     72109,     72110,     72111,     72112,     73083,
     *     73084,     73085,     73086,     73087,     73088,     73089,
     *     73090,     73091,     73092,     73093,     73094,     73095,
     *     73096,     73097,     73098,     73099,     73100,     73101,
     *     73102,     73103,     73104,     73105,     73106,     73107,
     *     73108,     73109,     73110,     73111,     73112,     73113,
     *     74084,     74085,     74086,     74087,     74088,     74089,
     *     74090,     74091,     74092,     74093/
      data i16/
     *     74094,     74095,     74096,     74097,     74098,     74099,
     *     74100,     74101,     74102,     74103,     74104,     74105,
     *     74106,     74107,     74108,     74109,     74110,     74111,
     *     74112,     74113,     74114,     74115,     74116,     75086,
     *     75087,     75088,     75089,     75090,     75091,     75092,
     *     75093,     75094,     75095,     75096,     75097,     75098,
     *     75099,     75100,     75101,     75102,     75103,     75104,
     *     75105,     75106,     75107,     75108,     75109,     75110,
     *     75111,     75112,     75113,     75114,     75115,     75116,
     *     75117,     76087,     76088,     76089,     76090,     76091,
     *     76092,     76093,     76094,     76095,     76096,     76097,
     *     76098,     76099,     76100,     76101,     76102,     76103,
     *     76104,     76105,     76106,     76107,     76108,     76109,
     *     76110,     76111,     76112,     76113,     76114,     76115,
     *     76116,     76117,     76118,     76119,     76120,     77089,
     *     77090,     77091,     77092,     77093,     77094,     77095,
     *     77096,     77097,     77098,     77099/
      data i17/
     *     77100,     77101,     77102,     77103,     77104,     77105,
     *     77106,     77107,     77108,     77109,     77110,     77111,
     *     77112,     77113,     77114,     77115,     77116,     77117,
     *     77118,     77119,     77120,     77121,     78090,     78091,
     *     78092,     78093,     78094,     78095,     78096,     78097,
     *     78098,     78099,     78100,     78101,     78102,     78103,
     *     78104,     78105,     78106,     78107,     78108,     78109,
     *     78110,     78111,     78112,     78113,     78114,     78115,
     *     78116,     78117,     78118,     78119,     78120,     78121,
     *     78122,     78123,     79094,     79095,     79096,     79097,
     *     79098,     79099,     79100,     79101,     79102,     79103,
     *     79104,     79105,     79106,     79107,     79108,     79109,
     *     79110,     79111,     79112,     79113,     79114,     79115,
     *     79116,     79117,     79118,     79119,     79120,     79121,
     *     79122,     79123,     79124,     79125,     80095,     80096,
     *     80097,     80098,     80099,     80100,     80101,     80102,
     *     80103,     80104,     80105,     80106/
      data i18/
     *     80107,     80108,     80109,     80110,     80111,     80112,
     *     80113,     80114,     80115,     80116,     80117,     80118,
     *     80119,     80120,     80121,     80122,     80123,     80124,
     *     80125,     80126,     80127,     81098,     81099,     81100,
     *     81101,     81102,     81103,     81104,     81105,     81106,
     *     81107,     81108,     81109,     81110,     81111,     81112,
     *     81113,     81114,     81115,     81116,     81117,     81118,
     *     81119,     81120,     81121,     81122,     81123,     81124,
     *     81125,     81126,     81127,     81128,     81129,     82101,
     *     82102,     82103,     82104,     82105,     82106,     82107,
     *     82108,     82109,     82110,     82111,     82112,     82113,
     *     82114,     82115,     82116,     82117,     82118,     82119,
     *     82120,     82121,     82122,     82123,     82124,     82125,
     *     82126,     82127,     82128,     82129,     82130,     82131,
     *     82132,     83105,     83106,     83107,     83108,     83109,
     *     83110,     83111,     83112,     83113,     83114,     83115,
     *     83116,     83117,     83118,     83119/
      data i19/
     *     83120,     83121,     83122,     83123,     83124,     83125,
     *     83126,     83127,     83128,     83129,     83130,     83131,
     *     83132,     83133,     84108,     84109,     84110,     84111,
     *     84112,     84113,     84114,     84115,     84116,     84117,
     *     84118,     84119,     84120,     84121,     84122,     84123,
     *     84124,     84125,     84126,     84127,     84128,     84129,
     *     84130,     84131,     84132,     84133,     84134,     85109,
     *     85110,     85111,     85112,     85113,     85114,     85115,
     *     85116,     85117,     85118,     85119,     85120,     85121,
     *     85122,     85123,     85124,     85125,     85126,     85127,
     *     85128,     85129,     85130,     85131,     85132,     85133,
     *     85134,     85135,     86113,     86114,     86115,     86116,
     *     86117,     86118,     86119,     86120,     86121,     86122,
     *     86123,     86124,     86125,     86126,     86127,     86128,
     *     86129,     86130,     86131,     86132,     86133,     86134,
     *     86135,     86136,     87114,     87115,     87116,     87117,
     *     87118,     87119,     87120,     87121/
      data i20/
     *     87122,     87123,     87124,     87125,     87126,     87127,
     *     87128,     87129,     87130,     87131,     87132,     87133,
     *     87134,     87135,     87136,     87137,     87138,     87139,
     *     87140,     87141,     88116,     88117,     88118,     88119,
     *     88120,     88121,     88122,     88123,     88124,     88125,
     *     88126,     88127,     88128,     88129,     88130,     88131,
     *     88132,     88133,     88134,     88135,     88136,     88137,
     *     88138,     88139,     88140,     88141,     88142,     89120,
     *     89121,     89122,     89123,     89124,     89125,     89126,
     *     89127,     89128,     89129,     89130,     89131,     89132,
     *     89133,     89134,     89135,     89136,     89137,     89138,
     *     89139,     89140,     89141,     89142,     89143,     90122,
     *     90123,     90124,     90125,     90126,     90127,     90128,
     *     90129,     90130,     90131,     90132,     90133,     90134,
     *     90135,     90136,     90137,     90138,     90139,     90140,
     *     90141,     90142,     90143,     90144,     90145,     91124,
     *     91125,     91126,     91127,     91128/
      data i21/
     *     91129,     91130,     91131,     91132,     91133,     91134,
     *     91135,     91136,     91137,     91138,     91139,     91140,
     *     91141,     91142,     91143,     91144,     91145,     91146,
     *     91147,     92134,     92135,     92136,     92137,     92138,
     *     92139,     92140,     92141,     92142,     92143,     92144,
     *     92145,     92146,     92147,     92148,     93136,     93137,
     *     93138,     93139,     93140,     93141,     93142,     93143,
     *     93144,     93145,     93146,     93147,     93148,     93149,
     *     93150,     94137,     94138,     94139,     94140,     94141,
     *     94142,     94143,     94144,     94145,     94146,     94147,
     *     94148,     94149,     94150,     94151,     94152,     95138,
     *     95139,     95140,     95141,     95142,     95143,     95144,
     *     95145,     95146,     95147,     95148,     95149,     95150,
     *     95151,     95152,     95153,     96139,     96140,     96141,
     *     96142,     96143,     96144,     96145,     96146,     96147,
     *     96148,     96149,     96150,     96151,     96152,     96153,
     *     96154,     96155,     97140,     97141/
      data i22/
     *     97142,     97143,     97144,     97145,     97146,     97147,
     *     97148,     97149,     97150,     97151,     97152,     97153,
     *     97154,     97155,     98141,     98142,     98143,     98144,
     *     98145,     98146,     98147,     98148,     98149,     98150,
     *     98151,     98152,     98153,     98154,     98155,     98156,
     *     99142,     99143,     99144,     99145,     99146,     99147,
     *     99148,     99149,     99150,     99151,     99152,     99153,
     *     99154,     99155,     99156,     99157,    100143,    100144,
     *    100145,    100146,    100147,    100148,    100149,    100150,
     *    100151,    100152,    100153,    100154,    100155,    100156,
     *    100157,    101146,    101147,    101148,    101149,    101150,
     *    101151,    101152,    101153,    101154,    101155,    101156,
     *    101157,    102149,    102150,    102151,    102152,    102153,
     *    102154,    102155,    102156,    102157,    103150,    103151,
     *    103152,    103153,    103154,    103155,    103156,    103157,
     *    104151,    104152,    104153,    104154,    104155,    104156,
     *    104157,    105152,    105153,    105154/
      data i23/
     *    105155,    105156,    105157,    106154,    106155,    106156,
     *    106157,    107155,    107156,    107157,    108156,    108157,
     *    109157,         0,         0,         0,         0,         0,
     *         0,         0,         0,         0,         0,         0,
     *         0,         0,         0,         0,         0,         0,
     *         0,         0,         0,         0,         0,         0,
     *         0,         0,         0,         0,         0,         0,
     *         0,         0,         0,         0,         0,         0,
     *         0,         0,         0,         0,         0,         0,
     *         0,         0,         0,         0,         0,         0,
     *         0,         0,         0,         0,         0,         0,
     *         0,         0,         0,         0,         0,         0,
     *         0,         0,         0,         0,         0,         0,
     *         0,         0,         0,         0,         0,         0,
     *         0,         0,         0,         0,         0,         0,
     *         0,         0,         0,         0,         0,         0,
     *         0,         0,         0,         0/
      data a1/
     * 139.05655,  134.9734, 140.07855,  8.071369, 32.285476,  7.28903 ,
     * 13.12896, 14.94991 , 25.84    , 14.93132 ,  2.42492 , 11.39    ,
     * 17.5923  , 26.11    , 31.598   , 40.81    , 25.12    , 11.68    ,
     * 14.0856  , 14.9068  , 20.9454  , 24.9539  , 33.83    , 40.9     ,
     * 18.374   , 15.7687  ,  4.9390 , 11.3414  , 12.607   , 20.174   ,
     * 25.077   , 34.95    , 41.02    , 27.87    , 22.9203  , 12.4158  ,
     * 12.05078 ,  8.668   , 13.3695  , 16.5623  , 23.664   , 28.97    ,
     * 37.64    , 44.01    , 35.095   , 28.9132  , 15.7017  , 10.6501  ,
     *  0.      ,  3.125025,  3.01991 ,  9.8732  , 13.694   , 21.03    ,
     * 24.89    , 32.76    , 38.03    , 39.7     , 24.91    , 17.3381  ,
     *  5.34552 ,  2.863436,  0.1015  ,  5.6821  ,  7.871   , 13.117   ,
     * 15.873   , 22.1     , 26.05    , 32.06    , 23.111   ,  8.00656 ,
     *  2.8555  , -4.73703 , -0.8093  , -0.7822  ,  3.3322  ,  3.7963  ,
     *  8.13    ,  9.44    , 17.46    , 33.61    , 16.77    , 10.68    ,
     *  1.95154 ,  0.8732  , -1.4874  , -0.01733 , -0.048   ,  2.83    ,
     *  3.35    ,  8.75    , 12.54    , 23.989   , 16.48    ,  5.319   ,
     *  1.751   , -7.0462  , -5.7354  , -8.0266/
      data a2/
     * -5.1555  , -5.95    , -2.16    ,  0.44    ,  6.75    , 25.32    ,
     * 12.929   ,  6.841   , -2.1886  , -5.1846  , -9.5314  , -8.4195  ,
     * -9.359   , -6.906   , -5.65    , -1.14    ,  2.64    ,  8.2     ,
     * 11.81    , 16.53    , 21.45    , 26.64    , 17.572   , 10.914   ,
     * -0.3966  , -5.4731  ,-13.9331  ,-13.1925  ,-16.214   ,-14.5862  ,
     *-15.0188  ,-10.728   , -9.1     , -3.79    , -1.77    ,  3.93    ,
     *  6.94    , 13.56    , 18.04    ,  6.767   , -0.055   , -8.9154  ,
     *-12.2099  ,-17.1968  ,-16.8506  ,-18.215   ,-15.89    ,-15.09    ,
     *-11.18    , -9.27    , -4.36    , -1.44    ,  3.91    , 10.755   ,
     *  3.827   , -7.144   ,-12.3853  ,-21.4924  ,-21.895   ,-24.4332  ,
     *-22.9502  ,-24.0808  ,-20.57    ,-19.86    ,-14.54    ,-12.76    ,
     * -7.      , -4.66    , 11.26    , -0.75    , -7.161   ,-16.9505  ,
     *-20.2074  ,-24.4407  ,-24.3058  ,-26.338   ,-24.5579  ,-24.94    ,
     *-20.89    ,-19.1     ,-14.66    ,-12.3     , -7.62    ,  4.13    ,
     * -3.16    ,-14.063   ,-19.0452  ,-26.01618 ,-26.58651 ,-29.93225 ,
     *-28.84689 ,-30.66444 ,-26.89659 ,-26.862   ,-23.      ,-22.52    ,
     *-17.87    ,-16.42    ,  4.84    , -7.07/
      data a3/
     *-13.33    ,-21.0039  ,-24.43992 ,-29.01372 ,-29.52215 ,-31.76175 ,
     *-29.79823 ,-29.804   ,-27.54    ,-27.4     ,-24.42    ,-23.13    ,
     *-20.01    , -2.18    , -9.38    ,-18.379   ,-23.0488  ,-30.23139 ,
     *-30.9479  ,-34.7147  ,-33.242   ,-35.0396  ,-33.0674  ,-34.42    ,
     *-31.98    ,-32.262   ,-29.72    ,-29.72    , -1.48    ,-11.168   ,
     *-17.426   ,-24.7994  ,-28.8017  ,-33.8066  ,-33.5348  ,-35.5597  ,
     *-35.0227  ,-36.592   ,-35.81    ,-36.611   ,-35.42    ,-35.698   ,
     *-32.124   ,-30.79    , -6.44    ,-13.16    ,-22.06    ,-27.276   ,
     *-34.8469  ,-35.1383  ,-38.5478  ,-38.4094  ,-41.4699  ,-40.8134  ,
     *-43.138   ,-42.343   ,-44.216   ,-41.291   ,-39.571   ,-35.94    ,
     * -4.46    ,-14.18    ,-20.527   ,-28.6434  ,-32.1241  ,-36.1887  ,
     *-37.8154  ,-41.0699  ,-41.7592  ,-44.3313  ,-44.493   ,-46.555   ,
     *-44.538   ,-43.22    ,-40.04    , -9.064   ,-15.7     ,-25.122   ,
     *-29.321   ,-37.5491  ,-39.0073  ,-44.1257  ,-44.9319  ,-48.4871  ,
     *-48.5581  ,-51.4262  ,-49.7272  ,-49.464   ,-46.83    ,-45.43    ,
     * -8.22    ,-17.92    ,-23.8     ,-31.875   ,-37.0753  ,-42.0048  ,
     *-44.472   ,-47.9563  ,-49.2197  ,-52.1997/
      data a4/
     *-51.4396  ,-51.847   ,-49.889   ,-49.15    ,-46.11    ,-13.22    ,
     *-19.46    ,-29.472   ,-34.554   ,-42.818   ,-45.3289  ,-50.2579  ,
     *-51.4483  ,-55.4152  ,-55.2834  ,-56.931   ,-55.1059  ,-55.291   ,
     *-52.69    ,-52.05    ,-12.47    ,-22.65    ,-29.22    ,-37.611   ,
     *-42.626   ,-48.2397  ,-50.7034  ,-54.6872  ,-55.554   ,-57.7092  ,
     *-56.9084  ,-57.488   ,-55.83    ,-55.477   ,-52.9     ,-24.47    ,
     *-34.47    ,-40.218   ,-48.331   ,-50.9436  ,-56.2508  ,-57.4774  ,
     *-60.6041  ,-60.1789  ,-62.1522  ,-60.6619  ,-61.407   ,-58.919   ,
     *-58.896   ,-55.19    ,-27.42    ,-34.3     ,-42.64    ,-48.0093  ,
     *-54.026   ,-56.038   ,-59.3425  ,-59.8443  ,-62.2265  ,-61.6471  ,
     *-62.8971  ,-61.424   ,-61.839   ,-59.791   ,-59.16    ,-29.41    ,
     *-39.21    ,-45.33    ,-53.902   ,-56.0774  ,-60.2251  ,-61.1536  ,
     *-64.4707  ,-64.2196  ,-66.7457  ,-65.5128  ,-67.098   ,-65.1248  ,
     *-66.028   ,-63.742   ,-63.482   ,-60.46    ,-31.63    ,-38.5     ,
     *-47.38    ,-51.6624  ,-56.353   ,-58.3438  ,-61.9812  ,-62.797   ,
     *-65.5787  ,-65.4235  ,-67.261   ,-66.2556  ,-67.303   ,-65.56    ,
     *-65.741   ,-62.982   ,-62.82    ,-32.61/
      data a5/
     *-42.21    ,-47.26    ,-54.185   ,-56.343   ,-61.17    ,-62.2116  ,
     *-66.0017  ,-65.9102  ,-68.8987  ,-67.8793  ,-70.0061  ,-68.417   ,
     *-69.56    ,-67.322   ,-68.134   ,-65.41    ,-65.707   ,-62.7     ,
     *-62.46    ,-58.91    ,-57.96    ,-47.54    ,-51.999   ,-56.69    ,
     *-58.837   ,-62.6544  ,-63.724   ,-66.8784  ,-67.085   ,-69.3225  ,
     *-68.905   ,-70.1415  ,-68.5912  ,-69.705   ,-68.06    ,-68.466   ,
     *-66.44    ,-66.41    ,-63.56    ,-62.76    ,-59.38    ,-57.99    ,
     *-47.39    ,-54.43    ,-56.41    ,-61.622   ,-62.656   ,-66.978   ,
     *-67.097   ,-70.5615  ,-69.9059  ,-72.5836  ,-71.2947  ,-73.4234  ,
     *-71.8575  ,-73.2145  ,-71.2154  ,-71.863   ,-69.53    ,-69.38    ,
     *-66.31    ,-65.38    ,-61.24    ,-47.31    ,-52.07    ,-56.65    ,
     *-58.88    ,-63.08    ,-64.34    ,-67.893   ,-68.228   ,-70.955   ,
     *-70.8611  ,-73.0351  ,-72.2907  ,-73.9185  ,-72.816   ,-73.639   ,
     *-72.165   ,-72.535   ,-70.078   ,-69.88    ,-66.08    ,-63.51    ,
     *-46.86    ,-54.08    ,-56.29    ,-61.59    ,-63.09    ,-67.897   ,
     *-68.215   ,-72.215   ,-72.1713  ,-75.2541  ,-74.6016  ,-77.0281  ,
     *-75.9191  ,-77.7621  ,-76.3918  ,-77.5961/
      data a6/
     *-75.343   ,-75.952   ,-72.42    ,-70.54    ,-66.71    ,-46.79    ,
     *-51.19    ,-56.59    ,-59.03    ,-63.64    ,-65.3     ,-69.161   ,
     *-70.302   ,-73.237   ,-73.454   ,-76.07    ,-75.891   ,-77.977   ,
     *-77.499   ,-79.011   ,-77.778   ,-78.607   ,-75.64    ,-73.88    ,
     *-70.72    ,-68.42    ,-64.26    ,-46.49    ,-53.97    ,-56.89    ,
     *-62.14    ,-64.246   ,-68.969   ,-70.227   ,-74.151   ,-74.442   ,
     *-77.892   ,-77.697   ,-80.591   ,-79.983   ,-82.431   ,-81.477   ,
     *-83.262   ,-80.706   ,-79.687   ,-76.72    ,-74.96    ,-71.37    ,
     *-68.68    ,-64.15    ,-46.59    ,-51.75    ,-57.28    ,-60.58    ,
     *-64.95    ,-66.98    ,-70.837   ,-72.173   ,-75.461   ,-76.202   ,
     *-79.044   ,-79.746   ,-82.1644  ,-82.7434  ,-84.5923  ,-82.6     ,
     *-81.713   ,-79.353   ,-77.794   ,-74.836   ,-72.679   ,-68.529   ,
     *-65.808   ,-61.14    ,-58.28    ,-54.06    ,-50.86    ,-57.89    ,
     *-63.65    ,-65.34    ,-70.19    ,-71.47    ,-75.997   ,-76.788   ,
     *-80.64    ,-81.099   ,-84.5177  ,-84.8746  ,-87.9162  ,-86.21    ,
     *-85.942   ,-83.661   ,-82.956   ,-80.121   ,-78.836   ,-75.09    ,
     *-72.89    ,-68.8     ,-66.49    ,-62.18/
      data a7/
     *-60.02    ,-58.24    ,-61.19    ,-65.95    ,-68.18    ,-72.38    ,
     *-74.23    ,-77.839   ,-79.278   ,-83.0136  ,-84.294   ,-87.7022  ,
     *-86.4881  ,-86.347   ,-84.844   ,-84.235   ,-82.348   ,-81.214   ,
     *-78.3     ,-76.27    ,-72.37    ,-70.13    ,-66.72    ,-64.38    ,
     *-58.79    ,-64.18    ,-66.36    ,-71.43    ,-73.15    ,-77.98    ,
     *-79.348   ,-83.626   ,-84.869   ,-88.7697  ,-87.8928  ,-88.4566  ,
     *-87.1193  ,-87.2679  ,-85.6593  ,-85.442   ,-82.95    ,-81.288   ,
     *-77.74    ,-76.62    ,-73.1     ,-71.76    ,-67.61    ,-66.74    ,
     *-69.58    ,-74.18    ,-76.43    ,-80.622   ,-82.659   ,-86.638   ,
     *-86.4507  ,-87.2098  ,-86.3679  ,-86.7836  ,-85.605   ,-85.6083  ,
     *-83.528   ,-82.327   ,-79.95    ,-78.88    ,-76.35    ,-75.11    ,
     *-71.78    ,-70.14    ,-67.44    ,-72.83    ,-75.004   ,-80.172   ,
     *-82.2     ,-86.808   ,-86.804   ,-88.4132  ,-87.7092  ,-88.7921  ,
     *-87.5421  ,-88.1132  ,-85.9674  ,-86.186   ,-83.513   ,-83.559   ,
     *-80.61    ,-80.48    ,-77.14    ,-76.43    ,-72.51    ,-68.      ,
     *-70.97    ,-75.98    ,-78.938   ,-83.606   ,-84.157   ,-86.018   ,
     *-85.819   ,-87.222   ,-86.429   ,-87.3244/
      data a8/
     *-86.0174  ,-86.325   ,-84.573   ,-84.606   ,-82.48    ,-82.14    ,
     *-79.63    ,-78.96    ,-75.99    ,-74.91    ,-68.4     ,-74.41    ,
     *-77.27    ,-82.567   ,-83.449   ,-86.071   ,-86.111   ,-88.225   ,
     *-87.618   ,-89.2199  ,-87.9506  ,-89.0995  ,-87.2606  ,-88.098   ,
     *-85.937   ,-86.33    ,-83.71    ,-83.7     ,-80.81    ,-80.34    ,
     *-76.92    ,-69.11    ,-72.97    ,-78.34    ,-79.63    ,-82.6     ,
     *-83.168   ,-85.519   ,-85.59    ,-87.413   ,-86.803   ,-88.027   ,
     *-86.954   ,-87.853   ,-86.37    ,-86.862   ,-85.09    ,-85.014   ,
     *-82.94    ,-82.32    ,-79.73    ,-78.84    ,-70.15    ,-76.37    ,
     *-77.8     ,-81.299   ,-82.192   ,-85.207   ,-85.431   ,-87.902   ,
     *-87.455   ,-89.397   ,-88.419   ,-89.91    ,-88.374   ,-89.522   ,
     *-87.604   ,-88.337   ,-86.02    ,-86.329   ,-83.68    ,-83.54    ,
     *-80.49    ,-80.11    ,-70.9     ,-73.07    ,-76.76    ,-78.12    ,
     *-81.22    ,-82.02    ,-84.78    ,-85.118   ,-87.076   ,-86.944   ,
     *-88.407   ,-87.605   ,-88.72    ,-87.458   ,-88.218   ,-86.623   ,
     *-87.041   ,-84.99    ,-84.95    ,-82.72    ,-82.25    ,-79.58    ,
     *-78.59    ,-75.77    ,-74.55    ,-69.99/
      data a9/
     *-74.31    ,-75.69    ,-79.7     ,-80.62    ,-83.974   ,-84.339   ,
     *-87.132   ,-86.99    ,-89.26    ,-88.536   ,-90.351   ,-89.255   ,
     *-90.5823  ,-89.0511  ,-90.0229  ,-88.0924  ,-88.721   ,-86.417   ,
     *-86.71    ,-83.94    ,-83.973   ,-80.95    ,-80.58    ,-77.32    ,
     *-68.41    ,-70.42    ,-74.42    ,-75.97    ,-79.589   ,-80.59    ,
     *-83.57    ,-84.135   ,-86.505   ,-86.41    ,-88.391   ,-87.994   ,
     *-89.367   ,-88.57    ,-89.534   ,-88.247   ,-88.943   ,-87.45    ,
     *-87.73    ,-85.8     ,-85.84    ,-83.58    ,-83.42    ,-81.06    ,
     *-80.42    ,-77.81    ,-77.01    ,-74.      ,-73.03    ,-69.99    ,
     *-68.55    ,-64.8     ,-66.92    ,-71.47    ,-73.27    ,-77.29    ,
     *-78.37    ,-82.09    ,-82.63    ,-85.83    ,-85.939   ,-88.654   ,
     *-88.328   ,-90.557   ,-90.032   ,-91.523   ,-90.3964  ,-91.6516  ,
     *-90.0657  ,-91.1017  ,-89.2018  ,-89.9454  ,-87.8202  ,-88.2371  ,
     *-85.8982  ,-86.021   ,-83.504   ,-83.31    ,-80.63    ,-80.19    ,
     *-77.37    ,-76.61    ,-69.99    ,-59.27    ,-64.09    ,-66.49    ,
     *-70.67    ,-72.53    ,-76.25    ,-77.53    ,-80.84    ,-81.589   ,
     *-84.421   ,-84.67    ,-87.002   ,-86.816/
      data a10/
     *-88.641   ,-87.995   ,-89.472   ,-88.421   ,-89.5907  ,-88.3259  ,
     *-89.2229  ,-87.619   ,-88.258   ,-86.4     ,-86.705   ,-84.6     ,
     *-84.631   ,-82.36    ,-82.02    ,-79.73    ,-79.04    ,-74.      ,
     *-70.31    ,-58.05    ,-60.51    ,-65.62    ,-67.65    ,-72.14    ,
     *-73.47    ,-77.3     ,-78.32    ,-81.76    ,-82.25    ,-85.28    ,
     *-85.11    ,-87.647   ,-87.178   ,-89.38    ,-88.542   ,-90.309   ,
     *-89.1717  ,-90.5252  ,-89.025   ,-90.0673  ,-88.286   ,-88.993   ,
     *-87.008   ,-87.348   ,-85.206   ,-85.217   ,-82.99    ,-82.41    ,
     *-77.85    ,-74.41    ,-69.48    ,-52.55    ,-57.76    ,-60.49    ,
     *-64.97    ,-67.12    ,-71.12    ,-72.86    ,-76.3     ,-77.52    ,
     *-80.61    ,-81.25    ,-83.81    ,-83.98    ,-86.263   ,-86.075   ,
     *-87.939   ,-87.368   ,-88.8468  ,-87.912   ,-88.984   ,-87.738   ,
     *-88.506   ,-86.897   ,-87.455   ,-85.71    ,-85.91    ,-83.97    ,
     *-83.813   ,-79.51    ,-76.5     ,-72.31    ,-68.88    ,-51.75    ,
     *-54.38    ,-59.88    ,-62.13    ,-66.91    ,-68.67    ,-73.02    ,
     *-74.29    ,-78.05    ,-78.82    ,-82.03    ,-82.49    ,-85.54    ,
     *-85.261   ,-87.6596  ,-87.1915  ,-89.162/
      data a11/
     *-88.323   ,-89.8608  ,-88.6974  ,-89.8811  ,-88.426   ,-89.29    ,
     *-87.665   ,-88.125   ,-86.509   ,-86.431   ,-82.385   ,-80.13    ,
     *-75.7     ,-73.02    ,-68.36    ,-65.55    ,-51.61    ,-54.71    ,
     *-59.55    ,-62.3     ,-66.23    ,-68.24    ,-72.2     ,-73.77    ,
     *-77.09    ,-78.16    ,-81.05    ,-81.72    ,-84.092   ,-84.334   ,
     *-86.231   ,-85.926   ,-87.536   ,-86.859   ,-88.079   ,-87.16    ,
     *-88.093   ,-86.913   ,-87.668   ,-86.361   ,-86.561   ,-82.9     ,
     *-80.715   ,-77.076   ,-74.515   ,-70.59    ,-67.79    ,-63.41    ,
     *-60.24    ,-55.69    ,-52.38    ,-47.59    ,-56.93    ,-61.95    ,
     *-63.96    ,-68.47    ,-70.14    ,-74.36    ,-75.26    ,-78.82    ,
     *-79.51    ,-82.66    ,-82.78    ,-85.478   ,-85.1     ,-87.299   ,
     *-86.721   ,-88.453   ,-87.572   ,-88.972   ,-87.873   ,-88.909   ,
     *-87.736   ,-88.276   ,-84.928   ,-83.294   ,-79.771   ,-77.91    ,
     *-74.07    ,-71.87    ,-68.04    ,-65.1     ,-61.26    ,-58.51    ,
     *-53.39    ,-77.98    ,-78.88    ,-81.38    ,-81.6     ,-83.76    ,
     *-83.74    ,-85.57    ,-85.27    ,-86.673   ,-86.04    ,-87.13    ,
     *-86.531   ,-87.238   ,-84.328   ,-83./
      data a12/
     *-80.025   ,-78.32    ,-74.85    ,-72.99    ,-69.37    ,-66.97    ,
     *-63.91    ,-61.19    ,-57.89    ,-79.4     ,-79.86    ,-82.44    ,
     *-82.57    ,-84.87    ,-84.657   ,-86.5     ,-85.91    ,-87.575   ,
     *-86.973   ,-88.089   ,-85.446   ,-84.542   ,-81.615   ,-80.442   ,
     *-77.1     ,-75.76    ,-72.16    ,-70.41    ,-67.29    ,-65.51    ,
     *-62.26    ,-75.34    ,-78.07    ,-78.77    ,-80.91    ,-81.38    ,
     *-83.2     ,-83.138   ,-84.844   ,-84.701   ,-86.027   ,-83.799   ,
     *-83.077   ,-80.76    ,-79.636   ,-76.78    ,-75.47    ,-72.46    ,
     *-70.988   ,-68.59    ,-67.16    ,-64.56    ,-75.57    ,-76.21    ,
     *-79.17    ,-79.4     ,-82.14    ,-82.04    ,-84.481   ,-84.213   ,
     *-85.959   ,-84.012   ,-83.757   ,-81.441   ,-80.935   ,-78.156   ,
     *-77.418   ,-74.385   ,-73.694   ,-70.957   ,-70.16    ,-67.37    ,
     *-65.77    ,-71.28    ,-74.1     ,-75.14    ,-77.52    ,-78.41    ,
     *-80.48    ,-81.07    ,-82.969   ,-81.424   ,-81.28    ,-79.45    ,
     *-79.052   ,-76.874   ,-76.074   ,-73.607   ,-73.4     ,-71.27    ,
     *-70.669   ,-68.47    ,-67.1     ,-64.48    ,-71.34    ,-72.09    ,
     *-75.41    ,-75.942   ,-78.986   ,-79.526/
      data a13/
     *-81.974   ,-80.66    ,-80.992   ,-79.276   ,-79.346   ,-77.147   ,
     *-77.061   ,-74.587   ,-74.773   ,-72.569   ,-72.466   ,-70.202   ,
     *-69.38    ,-66.87    ,-65.2     ,-67.21    ,-69.98    ,-71.59    ,
     *-74.38    ,-75.645   ,-77.998   ,-77.114   ,-77.555   ,-76.266   ,
     *-76.452   ,-74.798   ,-74.663   ,-72.897   ,-73.379   ,-71.749   ,
     *-71.829   ,-70.094   ,-69.473   ,-67.25    ,-66.059   ,-63.45    ,
     *-67.19    ,-68.48    ,-71.94    ,-72.95    ,-76.1     ,-75.505   ,
     *-76.278   ,-75.131   ,-75.766   ,-74.198   ,-74.719   ,-72.895   ,
     *-73.718   ,-72.082   ,-72.547   ,-70.835   ,-70.702   ,-68.573   ,
     *-67.954   ,-65.518   ,-64.26    ,-62.94    ,-66.2     ,-67.86    ,
     *-70.96    ,-70.67    ,-71.495   ,-71.102   ,-71.632   ,-70.869   ,
     *-71.316   ,-70.16    ,-71.26    ,-70.103   ,-70.773   ,-69.48    ,
     *-69.544   ,-67.848   ,-67.473   ,-65.66    ,-64.69    ,-62.12    ,
     *-62.86    ,-64.57    ,-67.98    ,-67.89    ,-69.325   ,-68.902   ,
     *-70.127   ,-69.146   ,-70.395   ,-69.166   ,-70.536   ,-69.434   ,
     *-70.419   ,-69.178   ,-69.683   ,-68.065   ,-68.19    ,-66.39    ,
     *-65.977   ,-63.622   ,-62.594   ,-59.94/
      data a14/
     *-58.37    ,-61.85    ,-62.22    ,-63.803   ,-63.74    ,-65.023   ,
     *-64.637   ,-66.064   ,-65.54    ,-66.89    ,-66.2     ,-67.342   ,
     *-66.397   ,-67.208   ,-66.051   ,-66.387   ,-64.939   ,-64.908   ,
     *-63.081   ,-62.292   ,-60.28    ,-58.806   ,-56.25    ,-58.02    ,
     *-58.5     ,-60.62    ,-60.67    ,-62.623   ,-62.36    ,-64.      ,
     *-63.42    ,-65.2     ,-64.573   ,-66.063   ,-65.209   ,-66.347   ,
     *-65.177   ,-65.952   ,-64.531   ,-64.935   ,-63.299   ,-62.999   ,
     *-60.931   ,-60.118   ,-57.728   ,-56.493   ,-53.77    ,-51.      ,
     *-51.74    ,-54.18    ,-54.63    ,-56.81    ,-56.97    ,-58.79    ,
     *-58.7     ,-60.57    ,-60.46    ,-62.01    ,-61.56    ,-62.738   ,
     *-61.99    ,-62.939   ,-61.888   ,-62.552   ,-61.321   ,-61.282   ,
     *-59.804   ,-59.219   ,-57.383   ,-56.267   ,-53.86    ,-52.3     ,
     *-49.6     ,-46.32    ,-47.27    ,-50.12    ,-50.74    ,-53.38    ,
     *-53.62    ,-56.023   ,-55.93    ,-58.06    ,-57.81    ,-59.75    ,
     *-59.37    ,-60.99    ,-60.176   ,-61.595   ,-60.598   ,-61.578   ,
     *-60.374   ,-60.772   ,-59.315   ,-59.264   ,-57.56    ,-56.953   ,
     *-54.704   ,-53.502   ,-50.997   ,-49.706/
      data a15/
     *-38.38    ,-39.63    ,-42.77    ,-43.72    ,-46.63    ,-47.24    ,
     *-49.85    ,-50.26    ,-52.51    ,-52.66    ,-54.74    ,-54.69    ,
     *-56.38    ,-56.12    ,-57.47    ,-57.11    ,-58.081   ,-57.332   ,
     *-57.836   ,-56.743   ,-56.887   ,-55.577   ,-55.173   ,-53.3946  ,
     *-52.3954  ,-50.336   ,-49.13    ,-46.69    ,-32.82    ,-34.44    ,
     *-37.86    ,-38.96    ,-42.3     ,-43.09    ,-46.06    ,-46.48    ,
     *-49.179   ,-49.39    ,-51.79    ,-51.65    ,-53.79    ,-53.47    ,
     *-55.21    ,-54.73    ,-56.13    ,-55.44    ,-56.39    ,-55.29    ,
     *-55.849   ,-54.486   ,-54.581   ,-52.8932  ,-52.4474  ,-50.476   ,
     *-49.7928  ,-47.4171  ,-46.063   ,-43.29    ,-41.5     ,-25.92    ,
     *-29.58    ,-31.      ,-34.6     ,-35.74    ,-38.92    ,-39.91    ,
     *-42.69    ,-43.44    ,-45.88    ,-46.31    ,-48.37    ,-48.61    ,
     *-50.28    ,-50.33    ,-51.54    ,-51.47    ,-52.49    ,-51.85    ,
     *-52.49    ,-51.48    ,-51.735   ,-50.54    ,-50.366   ,-48.94    ,
     *-48.445   ,-46.437   ,-45.299   ,-42.844   ,-41.403   ,-38.62    ,
     *-23.78    ,-25.55    ,-29.36    ,-30.61    ,-34.2     ,-35.15    ,
     *-38.36    ,-39.02    ,-41.899   ,-42.37/
      data a16/
     *-44.91    ,-45.02    ,-47.24    ,-47.24    ,-48.97    ,-48.71    ,
     *-50.15    ,-49.59    ,-50.68    ,-49.73    ,-50.45    ,-49.307   ,
     *-49.648   ,-48.259   ,-48.25    ,-46.37    ,-45.71    ,-43.393   ,
     *-42.517   ,-39.912   ,-38.676   ,-35.49    ,-34.27    ,-20.71    ,
     *-22.3     ,-26.11    ,-27.39    ,-30.84    ,-31.91    ,-34.85    ,
     *-35.82    ,-38.45    ,-39.09    ,-41.34    ,-41.68    ,-43.56    ,
     *-43.61    ,-45.24    ,-45.18    ,-46.23    ,-45.79    ,-46.62    ,
     *-45.85    ,-46.56    ,-45.45    ,-45.814   ,-44.218   ,-43.826   ,
     *-41.933   ,-41.224   ,-39.025   ,-37.987   ,-35.54    ,-34.361   ,
     *-31.79    ,-16.45    ,-20.46    ,-21.87    ,-25.64    ,-26.74    ,
     *-30.11    ,-30.88    ,-33.934   ,-34.57    ,-37.26    ,-37.54    ,
     *-39.95    ,-40.07    ,-42.03    ,-41.93    ,-43.55    ,-43.01    ,
     *-44.35    ,-43.53    ,-44.6     ,-43.51    ,-44.257   ,-42.811   ,
     *-43.007   ,-41.227   ,-41.145   ,-38.995   ,-38.717   ,-36.403   ,
     *-35.893   ,-33.406   ,-32.441   ,-29.7     ,-28.3     ,-13.17    ,
     *-17.14    ,-18.56    ,-22.14    ,-23.32    ,-26.36    ,-27.43    ,
     *-30.22    ,-31.06    ,-33.4     ,-34.02/
      data a17/
     *-36.      ,-36.29    ,-38.02    ,-37.95    ,-39.46    ,-39.15    ,
     *-40.32    ,-39.54    ,-40.31    ,-39.176   ,-39.73    ,-38.35    ,
     *-38.46    ,-36.72    ,-36.716   ,-34.857   ,-34.543   ,-32.538   ,
     *-31.702   ,-29.46    ,-28.29    ,-25.93    ,-11.15    ,-12.6     ,
     *-16.51    ,-17.71    ,-21.22    ,-22.1     ,-25.326   ,-25.96    ,
     *-28.94    ,-29.47    ,-31.95    ,-32.35    ,-34.35    ,-34.38    ,
     *-36.18    ,-35.74    ,-37.33    ,-36.61    ,-37.85    ,-36.83    ,
     *-37.832   ,-36.499   ,-37.338   ,-35.71    ,-36.311   ,-34.487   ,
     *-34.787   ,-32.821   ,-32.671   ,-30.446   ,-29.93    ,-27.43    ,
     *-26.625   ,-23.75    ,-12.84    ,-14.21    ,-17.34    ,-18.57    ,
     *-21.5     ,-22.58    ,-24.99    ,-25.8     ,-27.83    ,-28.33    ,
     *-30.13    ,-30.24    ,-31.85    ,-31.58    ,-33.11    ,-32.53    ,
     *-33.8     ,-32.896   ,-33.88    ,-32.796   ,-33.49    ,-32.278   ,
     *-32.591   ,-31.165   ,-31.165   ,-29.606   ,-29.119   ,-27.32    ,
     *-26.413   ,-24.37    ,-23.153   ,-20.22    , -8.27    ,-11.89    ,
     *-12.94    ,-16.323   ,-17.11    ,-20.26    ,-20.76    ,-23.52    ,
     *-23.89    ,-26.26    ,-26.17    ,-28.55/
      data a18/
     *-28.17    ,-30.2     ,-29.6     ,-31.3     ,-30.54    ,-32.      ,
     *-31.15    ,-32.238   ,-31.07    ,-31.851   ,-30.565   ,-30.979   ,
     *-29.571   ,-29.529   ,-27.687   ,-27.37    ,-25.292   ,-24.716   ,
     *-22.312   ,-20.969   ,-16.27    , -7.92    , -9.35    ,-12.32    ,
     *-13.55    ,-16.17    ,-17.07    ,-19.4     ,-20.02    ,-22.13    ,
     *-22.47    ,-24.4     ,-24.7     ,-26.24    ,-25.97    ,-27.46    ,
     *-27.09    ,-28.29    ,-27.52    ,-28.42    ,-27.52    ,-28.07    ,
     *-27.075   ,-27.205   ,-26.003   ,-25.784   ,-24.369   ,-23.846   ,
     *-22.278   ,-21.048   ,-16.778   ,-13.662   , -9.263   , -7.73    ,
     *-11.06    ,-11.65    ,-14.63    ,-15.06    ,-17.72    ,-17.9     ,
     *-20.43    ,-20.34    ,-22.55    ,-22.25    ,-24.24    ,-23.77    ,
     *-25.45    ,-24.82    ,-26.12    ,-25.27    ,-26.27    ,-25.3     ,
     *-25.957   ,-24.811   ,-25.132   ,-23.792   ,-23.809   ,-22.476   ,
     *-21.772   ,-17.638   ,-14.752   ,-10.493   , -7.572   , -3.25    ,
     * -0.1879  , -7.38    , -9.71    ,-10.73    ,-12.94    ,-13.6     ,
     *-15.66    ,-16.26    ,-17.98    ,-17.99    ,-19.67    ,-19.57    ,
     *-20.94    ,-20.41    ,-21.49    ,-20.81/
      data a19/
     *-21.59    ,-20.74    ,-21.085   ,-20.048   ,-20.078   ,-18.894   ,
     *-18.282   ,-14.815   ,-11.872   , -8.146   , -5.254   , -1.219   ,
     *  1.71    ,  5.96    , -7.98    , -8.37    ,-11.01    ,-11.17    ,
     *-13.47    ,-13.41    ,-15.5     ,-15.27    ,-17.04    ,-16.59    ,
     *-17.99    ,-17.35    ,-18.36    ,-17.56    ,-18.206   ,-17.168   ,
     *-17.492   ,-16.391   ,-15.977   ,-12.457   ,-10.394   , -6.676   ,
     * -4.494   , -0.542   ,  1.759   ,  5.83    ,  8.3517  , -0.81    ,
     * -3.11    , -3.97    , -6.14    , -6.94    , -8.77    , -8.97    ,
     *-10.77    ,-10.79    ,-12.31    ,-11.92    ,-13.05    ,-12.5     ,
     *-13.29    ,-12.56    ,-12.902   ,-11.992   ,-11.672   , -8.64    ,
     * -6.603   , -3.403   , -1.269   ,  2.226   ,  4.373   ,  8.089   ,
     * 10.52    , 14.29    , -1.61    , -4.      , -4.12    , -6.31    ,
     * -6.22    , -8.07    , -7.78    , -9.18    , -8.67    , -9.68    ,
     * -8.97    , -9.623   , -8.779   , -8.682   , -5.723   , -4.342   ,
     * -1.192   ,  0.232   ,  3.634   ,  5.198   ,  8.829   , 10.589   ,
     * 14.41    , 16.367   ,  3.83    ,  3.08    ,  0.93    ,  0.63    ,
     * -1.29    , -1.44    , -2.98    , -2.72/
      data a20/
     * -3.84    , -3.41    , -4.2     , -3.61    , -3.573   , -0.98    ,
     *  0.289   ,  2.96    ,  4.293   ,  7.036   ,  8.609   , 11.451   ,
     * 13.255   , 16.36    , 18.381   , 21.63    , 23.84    , 27.2     ,
     * 29.59    , 33.14    ,  6.02    ,  5.8     ,  3.54    ,  3.48    ,
     *  1.63    ,  1.79    ,  0.4     ,  0.8     , -0.22    ,  0.31    ,
     *  0.074   ,  2.51    ,  3.269   ,  5.863   ,  6.63    ,  9.365   ,
     * 10.25    , 12.938   , 14.301   , 17.234   , 18.803   , 21.987   ,
     * 23.6626  , 27.1726  , 28.936   , 32.48    , 34.46    ,  8.87    ,
     *  8.59    ,  7.07    ,  7.23    ,  6.09    ,  6.37    ,  5.97    ,
     *  8.06    ,  8.684   , 10.82    , 11.54    , 13.73    , 14.5     ,
     * 16.603   , 17.818   , 20.2     , 21.615   , 24.298   , 25.849   ,
     * 28.89    , 30.72    , 33.76    , 35.91    , 39.24    , 12.      ,
     * 12.06    , 10.65    , 10.89    , 10.27    , 12.16    , 12.346   ,
     * 14.45    , 14.646   , 16.916   , 17.183   , 19.243   , 19.98    ,
     * 22.283   , 23.18    , 25.805   , 26.748   , 29.58    , 30.8586  ,
     * 33.8121  , 35.4444  , 38.7294  , 40.607   , 44.25    , 17.66    ,
     * 17.66    , 17.      , 18.59    , 18.49/
      data a21/
     * 20.18    , 20.31    , 21.94    , 22.31    , 23.78    , 24.31    ,
     * 26.015   , 26.825   , 28.852   , 29.876   , 32.162   , 33.4222  ,
     * 35.923   , 37.4859  , 40.337   , 42.32    , 45.54    , 47.64    ,
     * 51.27    , 27.17    , 28.87    , 29.208   , 31.181   , 31.598   ,
     * 33.78    , 34.586   , 36.914   , 38.142   , 40.9155  , 42.4417  ,
     * 45.3872  , 47.306   , 50.5709  , 52.711   , 33.74    , 35.22    ,
     * 35.62    , 37.28    , 38.      , 39.95    , 41.0387  , 43.37    ,
     * 44.8683  , 47.4516  , 49.3069  , 52.21    , 54.26    , 57.41    ,
     * 59.921   , 38.39    , 38.349   , 40.02    , 40.333   , 42.16    ,
     * 42.879   , 45.086   , 46.1602  , 48.5849  , 50.1224  , 52.952   ,
     * 54.7139  , 57.751   , 59.801   , 63.174   , 65.365   , 43.17    ,
     * 44.34    , 44.64    , 46.      , 46.63    , 48.42    , 49.385   ,
     * 51.491   , 52.9312  , 55.4632  , 57.171   , 59.8763  , 61.893   ,
     * 64.991   , 67.23    , 70.59    , 48.02    , 47.87    , 49.15    ,
     * 49.39    , 51.09    , 51.701   , 53.696   , 54.8007  , 57.1773  ,
     * 58.449   , 60.998   , 62.6133  , 65.528   , 67.388   , 70.746   ,
     * 72.985   , 76.65    , 53.19    , 54.17/
      data a22/
     * 54.27    , 55.6     , 56.1     , 57.7     , 58.682   , 60.69    ,
     * 61.812   , 64.01    , 65.485   , 68.099   , 69.8436  , 72.948   ,
     * 75.23    , 78.53    , 58.24    , 58.02    , 59.17    , 59.33    ,
     * 60.91    , 61.459   , 63.377   , 64.0873  , 66.15    , 67.239   ,
     * 69.7179  , 71.167   , 74.128   , 76.03    , 79.296   , 81.337   ,
     * 63.82    , 64.69    , 64.72    , 65.96    , 66.38    , 67.93    ,
     * 68.55    , 70.27    , 71.11    , 73.27    , 74.507   , 77.263   ,
     * 79.0077  , 81.99    , 84.09    , 87.16    , 69.36    , 69.04    ,
     * 70.1     , 70.13    , 71.54    , 71.885   , 73.5     , 74.063   ,
     * 76.      , 76.817   , 79.339   , 80.897   , 83.787   , 85.482   ,
     * 88.585   , 76.04    , 77.08    , 77.26    , 78.6     , 79.02    ,
     * 80.54    , 81.24    , 83.49    , 84.842   , 87.522   , 89.03    ,
     * 91.82    , 82.78    , 82.856   , 84.33    , 84.723   , 86.87    ,
     * 87.796   , 90.22    , 91.42    , 94.018   , 88.64    , 89.73    ,
     * 90.05    , 91.74    , 92.67    , 94.75    , 95.85    , 98.1     ,
     * 94.31    , 94.28    , 95.89    , 96.35    , 98.3     , 99.02    ,
     *101.24    ,100.39    ,101.55    ,102.07/
      data a23/
     *103.44    ,104.16    ,105.97    ,106.91    ,108.22    ,108.47    ,
     *110.12    ,114.51    ,114.8     ,115.96    ,120.13    ,121.24    ,
     *128.21    ,  0.      ,  0.      ,  0.      ,  0.      ,  0.      ,
     *  0.      ,  0.      ,  0.      ,  0.      ,  0.      ,  0.      ,
     *  0.      ,  0.      ,  0.      ,  0.      ,  0.      ,  0.      ,
     *  0.      ,  0.      ,  0.      ,  0.      ,  0.      ,  0.      ,
     *  0.      ,  0.      ,  0.      ,  0.      ,  0.      ,  0.      ,
     *  0.      ,  0.      ,  0.      ,  0.      ,  0.      ,  0.      ,
     *  0.      ,  0.      ,  0.      ,  0.      ,  0.      ,  0.      ,
     *  0.      ,  0.      ,  0.      ,  0.      ,  0.      ,  0.      ,
     *  0.      ,  0.      ,  0.      ,  0.      ,  0.      ,  0.      ,
     *  0.      ,  0.      ,  0.      ,  0.      ,  0.      ,  0.      ,
     *  0.      ,  0.      ,  0.      ,  0.      ,  0.      ,  0.      ,
     *  0.      ,  0.      ,  0.      ,  0.      ,  0.      ,  0.      ,
     *  0.      ,  0.      ,  0.      ,  0.      ,  0.      ,  0.      ,
     *  0.      ,  0.      ,  0.      ,  0.      ,  0.      ,  0.      ,
     *  0.      ,  0.      ,  0.      ,  0./
      data e1/
     *     0.33 ,     2.5  ,     0.33 ,     0.013,     0.052,     0.011,
     *     0.022,     0.03 ,   380.   ,     0.03 ,     0.05 ,    50.   ,
     *     1.   ,    30.   ,     7.   ,   120.   ,   300.   ,    50.   ,
     *     0.7  ,     0.8  ,     0.8  ,     2.   ,   250.   ,   110.   ,
     *     5.   ,     0.8  ,     0.11 ,     0.4  ,     0.4  ,     6.   ,
     *    15.   ,   500.   ,   300.   ,    70.   ,     1.2  ,     1.   ,
     *     0.3  ,     0.4  ,     1.3  ,     1.1  ,    21.   ,    22.   ,
     *   500.   ,   700.   ,    24.   ,     2.2  ,     0.5  ,     0.9  ,
     *     0.   ,     0.016,     0.025,     0.8  ,     4.   ,    40.   ,
     *   130.   ,   240.   ,   910.   ,   400.   ,   150.   ,     1.   ,
     *     0.27 ,     0.024,     0.04 ,     2.3  ,    15.   ,    20.   ,
     *    19.   ,   400.   ,   600.   ,    40.   ,    10.   ,     0.08 ,
     *     0.5  ,     0.05 ,     0.4  ,     0.9  ,     2.8  ,     2.1  ,
     *    50.   ,    90.   ,   620.   ,   400.   ,   130.   ,     8.   ,
     *     0.24 ,     0.7  ,     0.14 ,     0.2  ,     8.   ,    30.   ,
     *   170.   ,   400.   ,   500.   ,    20.   ,    50.   ,     5.   ,
     *     0.6,     2.   ,     2.   ,     1.6/
      data e2/
     *     2.6  ,    10.   ,   100.   ,    70.   ,   400.   ,   400.   ,
     *    12.   ,     7.   ,     2.1  ,     1.7  ,     0.9  ,     0.9  ,
     *     6.   ,    16.   ,    40.   ,   140.   ,   150.   ,   250.   ,
     *   580.   ,   740.   ,  1140.   ,  3570.   ,    27.   ,    16.   ,
     *     1.5  ,     1.5  ,     0.7  ,     0.7  ,     0.8  ,     0.9  ,
     *     2.1  ,    30.   ,   210.   ,   410.   ,  1580.   ,   500.   ,
     *   760.   ,   650.   ,    90.   ,    25.   ,     4.   ,     1.1  ,
     *     0.8  ,     0.7  ,     0.7  ,     5.   ,    40.   ,    70.   ,
     *   300.   ,   400.   ,   500.   ,   510.   ,   600.   ,    19.   ,
     *    10.   ,     3.   ,     0.7  ,     0.6  ,     0.6  ,     0.6  ,
     *     0.7  ,     2.4  ,    50.   ,   300.   ,   410.   ,   400.   ,
     *   500.   ,   600.   ,   300.   ,    40.   ,     4.   ,     1.8  ,
     *     1.7  ,     0.6  ,     0.6  ,     1.5  ,     1.2  ,    80.   ,
     *    60.   ,   300.   ,   400.   ,   400.   ,   500.   ,   160.   ,
     *    50.   ,     3.   ,     1.5  ,     0.24 ,     0.21 ,     0.21 ,
     *     0.21 ,     0.25 ,     0.26 ,    12.   ,   200.   ,    40.   ,
     *   300.   ,   400.   ,   300.   ,    50./
      data e3/
     *     8.   ,     0.7  ,     0.27 ,     0.07 ,     0.07 ,     0.1  ,
     *     0.14 ,    19.   ,   500.   ,   160.   ,   200.   ,    60.   ,
     *   300.   ,    50.   ,    30.   ,     3.   ,     1.3  ,     0.27 ,
     *     0.6  ,     0.8  ,     5.   ,     1.3  ,     1.4  ,    40.   ,
     *    70.   ,    20.   ,    60.   ,    40.   ,   300.   ,    20.   ,
     *     8.   ,     1.4  ,     0.9  ,     1.1  ,     1.1  ,     1.1  ,
     *     1.5  ,    10.   ,    40.   ,    11.   ,    16.   ,     8.   ,
     *    24.   ,   400.   ,    40.   ,    22.   ,     5.   ,     2.1  ,
     *     1.2  ,     1.2  ,     1.2  ,     1.2  ,     1.3  ,     1.3  ,
     *     4.   ,     4.   ,     4.   ,     4.   ,     8.   ,    50.   ,
     *   300.   ,   200.   ,     4.   ,     1.5  ,     1.3  ,     2.3  ,
     *     2.2  ,     1.3  ,     1.3  ,     2.2  ,     5.   ,     4.   ,
     *    16.   ,    20.   ,   300.   ,    11.   ,    40.   ,     6.   ,
     *     7.   ,     1.4  ,     2.8  ,     1.3  ,     1.1  ,     1.1  ,
     *     1.1  ,     1.1  ,     1.4  ,     7.   ,   100.   ,   300.   ,
     *   300.   ,   200.   ,   100.   ,    17.   ,     1.4  ,     1.3  ,
     *     3.   ,     1.3  ,     1.5  ,     1.5/
      data e4/
     *     1.8  ,     3.   ,    15.   ,   100.   ,   300.   ,   180.   ,
     *   150.   ,    20.   ,    14.   ,     7.   ,     2.7  ,     1.6  ,
     *     1.6  ,     1.6  ,     1.6  ,     1.6  ,     1.6  ,    10.   ,
     *   200.   ,   300.   ,   400.   ,   200.   ,   100.   ,    24.   ,
     *     1.6  ,     1.5  ,     2.4  ,     1.6  ,     1.8  ,     1.5  ,
     *     1.5  ,     3.   ,    30.   ,    29.   ,   100.   ,   160.   ,
     *    60.   ,    15.   ,    12.   ,     2.2  ,     1.4  ,     1.4  ,
     *     1.5  ,     1.5  ,     1.5  ,     1.5  ,     4.   ,    20.   ,
     *    15.   ,    60.   ,   200.   ,   100.   ,    18.   ,     1.5  ,
     *     1.5  ,     2.5  ,     1.5  ,     1.8  ,     1.5  ,     1.5  ,
     *     1.7  ,    19.   ,    20.   ,    20.   ,    50.   ,   180.   ,
     *    50.   ,    11.   ,    11.   ,     3.   ,     1.5  ,     1.5  ,
     *     1.5  ,     1.5  ,     1.5  ,     1.5  ,     1.6  ,     1.6  ,
     *    16.   ,    19.   ,    18.   ,   200.   ,   300.   ,   200.   ,
     *   300.   ,     2.5  ,     2.1  ,     2.6  ,     1.8  ,     5.   ,
     *     1.5  ,     1.5  ,     1.8  ,     1.9  ,     8.   ,    50.   ,
     *     8.   ,    20.   ,   300.   ,   130./
      data e5/
     *   100.   ,    40.   ,    11.   ,    16.   ,    10.   ,     2.2  ,
     *     1.7  ,     1.8  ,     1.6  ,     1.6  ,     1.7  ,     1.8  ,
     *     3.   ,    11.   ,     6.   ,    40.   ,    19.   ,   150.   ,
     *   190.   ,   280.   ,   280.   ,   400.   ,    28.   ,   100.   ,
     *     4.   ,     2.   ,     3.   ,     1.7  ,     2.1  ,     2.9  ,
     *     3.   ,     2.3  ,     2.5  ,     6.   ,    70.   ,     7.   ,
     *   150.   ,   200.   ,   200.   ,   130.   ,   300.   ,   190.   ,
     *   320.   ,   250.   ,   100.   ,    13.   ,     5.   ,     6.   ,
     *     4.   ,     1.5  ,     1.8  ,     1.5  ,     1.5  ,     1.4  ,
     *     2.1  ,     1.6  ,     1.9  ,     4.   ,   100.   ,    30.   ,
     *   120.   ,   140.   ,   400.   ,   410.   ,    50.   ,   100.   ,
     *   100.   ,    30.   ,    20.   ,     4.   ,     4.   ,     4.   ,
     *     2.2  ,     1.6  ,     1.8  ,     2.2  ,    10.   ,     6.   ,
     *    24.   ,     6.   ,    25.   ,   220.   ,   300.   ,   400.   ,
     *   400.   ,   300.   ,    40.   ,   200.   ,   200.   ,    12.   ,
     *    11.   ,     1.5  ,     1.5  ,     1.5  ,     1.5  ,     1.5  ,
     *     1.6  ,     1.8  ,     1.9  ,     2.1/
      data e6/
     *     4.   ,    15.   ,   100.   ,   130.   ,   400.   ,   500.   ,
     *   360.   ,   300.   ,   200.   ,   250.   ,    15.   ,    20.   ,
     *    11.   ,     3.   ,     4.   ,     2.4  ,     2.4  ,     5.   ,
     *     5.   ,     5.   ,    26.   ,    19.   ,    60.   ,   120.   ,
     *   130.   ,   400.   ,   400.   ,   420.   ,   240.   ,   140.   ,
     *    60.   ,    20.   ,    12.   ,    29.   ,     8.   ,     6.   ,
     *     8.   ,     6.   ,     5.   ,     3.   ,     3.   ,     3.   ,
     *     5.   ,     5.   ,    14.   ,    50.   ,    30.   ,    80.   ,
     *    80.   ,   120.   ,   620.   ,   470.   ,   110.   ,    70.   ,
     *    40.   ,    30.   ,    27.   ,    19.   ,    23.   ,    17.   ,
     *    21.   ,     4.   ,     3.   ,     2.9  ,     2.9  ,     5.   ,
     *     8.   ,    13.   ,     8.   ,    14.   ,    15.   ,    26.   ,
     *    28.   ,    30.   ,    50.   ,    60.   ,   120.   ,   150.   ,
     *   300.   ,   200.   ,    30.   ,    40.   ,     9.   ,    21.   ,
     *     4.   ,     4.   ,     2.6  ,     2.6  ,     2.6  ,     4.   ,
     *     3.   ,     5.   ,    12.   ,    12.   ,     7.   ,    60.   ,
     *    40.   ,    70.   ,    80.   ,   130./
      data e7/
     *   610.   ,   500.   ,   400.   ,    70.   ,    90.   ,    60.   ,
     *   170.   ,    11.   ,    20.   ,     2.9  ,     3.   ,     2.8  ,
     *     2.8  ,     3.   ,    10.   ,    11.   ,     6.   ,     6.   ,
     *    40.   ,    60.   ,    60.   ,   100.   ,   410.   ,   510.   ,
     *   300.   ,   510.   ,   100.   ,   350.   ,   100.   ,   200.   ,
     *     8.   ,    10.   ,     4.   ,     2.4  ,     2.4  ,     2.4  ,
     *     2.4  ,     2.6  ,     2.6  ,     3.   ,     3.   ,    20.   ,
     *    70.   ,    80.   ,   140.   ,   400.   ,   540.   ,   400.   ,
     *   300.   ,    60.   ,   200.   ,    19.   ,     5.   ,     3.   ,
     *     2.8  ,     2.5  ,     2.5  ,     2.1  ,     4.   ,     2.7  ,
     *     6.   ,    13.   ,    30.   ,    70.   ,    70.   ,   230.   ,
     *   410.   ,   500.   ,   310.   ,   300.   ,    16.   ,     6.   ,
     *    13.   ,     4.   ,     4.   ,     2.4  ,     2.1  ,     2.   ,
     *     2.   ,     2.   ,     2.1  ,     6.   ,     6.   ,    21.   ,
     *   200.   ,   310.   ,   300.   ,   400.   ,   580.   ,   400.   ,
     *   300.   ,   200.   ,    26.   ,     4.   ,     5.   ,     6.   ,
     *     6.   ,     5.   ,     4.   ,     2.1/
      data e8/
     *     2.3  ,    24.   ,    10.   ,    11.   ,    70.   ,   200.   ,
     *   300.   ,   300.   ,   400.   ,   500.   ,   400.   ,   300.   ,
     *    90.   ,    13.   ,    12.   ,     8.   ,     8.   ,     6.   ,
     *     2.2  ,     2.2  ,     2.2  ,     2.3  ,     2.4  ,     5.   ,
     *     5.   ,    10.   ,   300.   ,   610.   ,   300.   ,   300.   ,
     *   400.   ,   400.   ,   300.   ,   150.   ,    13.   ,    30.   ,
     *    12.   ,    10.   ,    20.   ,    17.   ,    21.   ,     3.   ,
     *     3.   ,     5.   ,    10.   ,    18.   ,    50.   ,    22.   ,
     *   100.   ,   210.   ,   300.   ,   400.   ,   400.   ,   300.   ,
     *   300.   ,    21.   ,    16.   ,    12.   ,    18.   ,     4.   ,
     *     4.   ,     5.   ,     5.   ,     5.   ,     6.   ,     4.   ,
     *     4.   ,    18.   ,    50.   ,    25.   ,   150.   ,   150.   ,
     *   300.   ,   160.   ,   400.   ,   300.   ,   150.   ,    90.   ,
     *   120.   ,   110.   ,    50.   ,     7.   ,     9.   ,     6.   ,
     *     5.   ,     5.   ,     3.   ,     3.   ,     4.   ,    28.   ,
     *    20.   ,   110.   ,    70.   ,   120.   ,    50.   ,   100.   ,
     *    70.   ,   100.   ,   190.   ,   500./
      data e9/
     *   300.   ,   180.   ,   200.   ,    18.   ,    11.   ,    10.   ,
     *     6.   ,     7.   ,     6.   ,     4.   ,     3.   ,     3.   ,
     *     2.9  ,     2.9  ,     2.9  ,     2.9  ,     3.   ,    13.   ,
     *    20.   ,    60.   ,    19.   ,   150.   ,   210.   ,   300.   ,
     *   540.   ,   420.   ,   300.   ,   200.   ,    18.   ,    30.   ,
     *     8.   ,    11.   ,     7.   ,    30.   ,     6.   ,     5.   ,
     *     4.   ,     4.   ,     4.   ,     4.   ,     5.   ,   300.   ,
     *    18.   ,   170.   ,    28.   ,    50.   ,    30.   ,    50.   ,
     *    80.   ,    80.   ,    70.   ,   170.   ,   170.   ,   320.   ,
     *   240.   ,   630.   ,   500.   ,   400.   ,    80.   ,   200.   ,
     *   100.   ,    30.   ,    10.   ,    16.   ,     7.   ,     5.   ,
     *     4.   ,     3.   ,     3.   ,     3.   ,     2.9  ,     2.9  ,
     *     2.9  ,     2.7  ,     2.7  ,     2.7  ,     2.7  ,     1.6  ,
     *     2.1  ,    11.   ,    25.   ,    40.   ,   120.   ,   120.   ,
     *   140.   ,    80.   ,   210.   ,   630.   ,   500.   ,   500.   ,
     *   400.   ,   300.   ,    19.   ,   200.   ,   200.   ,    25.   ,
     *    24.   ,   200.   ,    20.   ,     6./
      data e10/
     *     9.   ,     4.   ,     8.   ,     8.   ,     2.7  ,     2.7  ,
     *     2.3  ,     2.3  ,     3.   ,    30.   ,     6.   ,    40.   ,
     *    22.   ,   100.   ,    70.   ,    80.   ,   210.   ,   150.   ,
     *   200.   ,   630.   ,   500.   ,   400.   ,    70.   ,   200.   ,
     *    70.   ,   160.   ,   200.   ,   200.   ,   200.   ,   100.   ,
     *    30.   ,    23.   ,     8.   ,    19.   ,    25.   ,     3.   ,
     *     2.   ,     1.7  ,     2.4  ,     2.4  ,     4.   ,     3.   ,
     *     4.   ,     4.   ,     4.   ,    21.   ,    80.   ,   110.   ,
     *   100.   ,   100.   ,   300.   ,   660.   ,   540.   ,   500.   ,
     *   400.   ,   300.   ,    50.   ,   280.   ,   200.   ,   160.   ,
     *   200.   ,   200.   ,   100.   ,   200.   ,    19.   ,     6.   ,
     *     5.   ,     4.   ,     2.6  ,     5.   ,     4.   ,     4.   ,
     *     4.   ,    10.   ,     4.   ,    20.   ,    30.   ,    60.   ,
     *    22.   ,    50.   ,    80.   ,    80.   ,   120.   ,   630.   ,
     *   510.   ,   400.   ,    80.   ,   290.   ,   280.   ,   260.   ,
     *   320.   ,   280.   ,   160.   ,   280.   ,    60.   ,   410.   ,
     *    16.   ,     2.1  ,     2.1  ,     7./
      data e11/
     *     6.   ,     1.6  ,     1.9  ,     1.5  ,     4.   ,     4.   ,
     *     7.   ,     7.   ,    11.   ,     7.   ,     7.   ,    50.   ,
     *    60.   ,    60.   ,    90.   ,   100.   ,   410.   ,   510.   ,
     *   500.   ,   310.   ,   180.   ,   130.   ,   100.   ,    80.   ,
     *    60.   ,    70.   ,    40.   ,    40.   ,    17.   ,    24.   ,
     *     8.   ,     6.   ,    14.   ,     8.   ,     6.   ,    12.   ,
     *     6.   ,     6.   ,     8.   ,     6.   ,     6.   ,    23.   ,
     *     7.   ,    16.   ,    22.   ,    30.   ,    40.   ,    60.   ,
     *    70.   ,   100.   ,   130.   ,   370.   ,   560.   ,   490.   ,
     *   410.   ,   410.   ,   350.   ,   310.   ,   300.   ,   300.   ,
     *   250.   ,   200.   ,   100.   ,    18.   ,    14.   ,     8.   ,
     *     8.   ,     8.   ,     7.   ,     6.   ,     6.   ,     6.   ,
     *     6.   ,     6.   ,     6.   ,    10.   ,    24.   ,    40.   ,
     *    50.   ,    80.   ,   130.   ,   110.   ,   510.   ,   620.   ,
     *   580.   ,   220.   ,   400.   ,    50.   ,   200.   ,   100.   ,
     *    50.   ,   100.   ,    30.   ,    12.   ,    70.   ,    50.   ,
     *     5.   ,     4.   ,     4.   ,    25./
      data e12/
     *     7.   ,    50.   ,   110.   ,   100.   ,    90.   ,   300.   ,
     *   330.   ,   500.   ,   630.   ,   300.   ,   220.   ,   210.   ,
     *   200.   ,   110.   ,    18.   ,    50.   ,    50.   ,    11.   ,
     *     7.   ,     4.   ,     4.   ,     4.   ,     4.   ,     5.   ,
     *    50.   ,    80.   ,    80.   ,   140.   ,   300.   ,   220.   ,
     *   500.   ,   300.   ,   220.   ,   200.   ,   200.   ,    50.   ,
     *    50.   ,    15.   ,     7.   ,     7.   ,     4.   ,     4.   ,
     *     4.   ,     5.   ,     8.   ,    70.   ,    40.   ,   100.   ,
     *    11.   ,   200.   ,   300.   ,   300.   ,   320.   ,   360.   ,
     *    60.   ,    60.   ,   200.   ,    50.   ,    19.   ,     9.   ,
     *     3.   ,     3.   ,     3.   ,     3.   ,     3.   ,     3.   ,
     *     4.   ,     4.   ,     4.   ,     4.   ,    30.   ,   200.   ,
     *   300.   ,   400.   ,   310.   ,   280.   ,    60.   ,    40.   ,
     *    30.   ,    50.   ,     4.   ,     4.   ,     4.   ,     7.   ,
     *     3.   ,    15.   ,     5.   ,    20.   ,     6.   ,    80.   ,
     *    16.   ,   100.   ,   200.   ,   300.   ,   450.   ,   160.   ,
     *   300.   ,    13.   ,    16.   ,    11./
      data e13/
     *     4.   ,     4.   ,     6.   ,     3.   ,     3.   ,     3.   ,
     *     3.   ,     3.   ,     3.   ,     3.   ,     3.   ,     3.   ,
     *    14.   ,   200.   ,   200.   ,   400.   ,    40.   ,   100.   ,
     *    40.   ,    21.   ,     5.   ,    10.   ,     4.   ,    21.   ,
     *     6.   ,     8.   ,     3.   ,     3.   ,     4.   ,     4.   ,
     *     4.   ,    10.   ,     7.   ,    80.   ,     9.   ,   200.   ,
     *   400.   ,   300.   ,   300.   ,    50.   ,    11.   ,    23.   ,
     *     4.   ,     5.   ,     8.   ,     4.   ,     3.   ,     3.   ,
     *     3.   ,     3.   ,     3.   ,     3.   ,     3.   ,     4.   ,
     *     3.   ,     4.   ,   120.   ,   400.   ,   300.   ,   150.   ,
     *    29.   ,    40.   ,     6.   ,    11.   ,     5.   ,    15.   ,
     *     5.   ,    50.   ,    13.   ,     5.   ,     4.   ,     4.   ,
     *     4.   ,     4.   ,     4.   ,    60.   ,    50.   ,   150.   ,
     *   250.   ,    70.   ,    60.   ,   200.   ,    11.   ,    23.   ,
     *     6.   ,     5.   ,     9.   ,    12.   ,     7.   ,     7.   ,
     *     5.   ,     4.   ,     4.   ,     3.   ,     3.   ,     3.   ,
     *     3.   ,     3.   ,     3.   ,    60./
      data e14/
     *   300.   ,   210.   ,   200.   ,    29.   ,    60.   ,    10.   ,
     *    12.   ,    24.   ,   200.   ,    50.   ,    30.   ,     7.   ,
     *    15.   ,     4.   ,     5.   ,     4.   ,     4.   ,     3.   ,
     *     3.   ,     6.   ,   100.   ,    20.   ,    50.   ,   280.   ,
     *   300.   ,    60.   ,   210.   ,    12.   ,    60.   ,   200.   ,
     *    90.   ,   110.   ,    12.   ,    28.   ,    11.   ,     4.   ,
     *     6.   ,     4.   ,     4.   ,     3.   ,     3.   ,     3.   ,
     *     3.   ,     4.   ,     4.   ,     5.   ,   200.   ,   430.   ,
     *   410.   ,   210.   ,   210.   ,    30.   ,    80.   ,   210.   ,
     *   200.   ,   200.   ,   100.   ,   100.   ,    60.   ,     6.   ,
     *    20.   ,     4.   ,    12.   ,     4.   ,     4.   ,     3.   ,
     *     3.   ,     3.   ,     6.   ,     5.   ,    50.   ,    50.   ,
     *   200.   ,   500.   ,   450.   ,   280.   ,   310.   ,    60.   ,
     *   210.   ,    16.   ,   210.   ,   200.   ,   220.   ,   200.   ,
     *   100.   ,   100.   ,    20.   ,     8.   ,     5.   ,     5.   ,
     *     5.   ,     3.   ,     3.   ,     3.   ,     3.   ,     3.   ,
     *     3.   ,     3.   ,     4.   ,    11./
      data e15/
     *   540.   ,   540.   ,   430.   ,   410.   ,   220.   ,    70.   ,
     *    60.   ,   320.   ,   220.   ,   210.   ,   300.   ,   200.   ,
     *   200.   ,   160.   ,   100.   ,    80.   ,     6.   ,    20.   ,
     *     4.   ,     4.   ,     4.   ,     4.   ,     3.   ,     3.   ,
     *     3.   ,    25.   ,    40.   ,    70.   ,   580.   ,   590.   ,
     *   500.   ,   450.   ,   280.   ,   310.   ,    60.   ,   220.   ,
     *    17.   ,   420.   ,   230.   ,   400.   ,   300.   ,   220.   ,
     *   130.   ,   100.   ,   200.   ,   200.   ,    50.   ,   100.   ,
     *     4.   ,     4.   ,     3.   ,     2.8  ,     2.8  ,     2.7  ,
     *     2.8  ,     2.8  ,     7.   ,    30.   ,    60.   ,   650.   ,
     *   550.   ,   540.   ,   430.   ,   410.   ,   220.   ,   210.   ,
     *   210.   ,   380.   ,   280.   ,   340.   ,   320.   ,   310.   ,
     *   200.   ,   200.   ,   200.   ,   190.   ,   220.   ,   100.   ,
     *   100.   ,   100.   ,     6.   ,   100.   ,     6.   ,     3.   ,
     *     3.   ,     3.   ,     3.   ,    26.   ,    14.   ,    60.   ,
     *   580.   ,   590.   ,   500.   ,   460.   ,   280.   ,   310.   ,
     *    60.   ,   220.   ,    18.   ,   370./
      data e16/
     *   240.   ,   350.   ,   300.   ,   280.   ,   270.   ,   220.   ,
     *   300.   ,   200.   ,   200.   ,   100.   ,   100.   ,    16.   ,
     *     5.   ,     7.   ,     3.   ,     3.   ,     3.   ,     3.   ,
     *     3.   ,     3.   ,     4.   ,   200.   ,   210.   ,   550.   ,
     *   550.   ,   430.   ,   410.   ,   230.   ,   220.   ,   210.   ,
     *   380.   ,   350.   ,   400.   ,   370.   ,   370.   ,   280.   ,
     *   350.   ,   360.   ,   280.   ,   200.   ,   210.   ,    50.   ,
     *    30.   ,   100.   ,   100.   ,     9.   ,     5.   ,     3.   ,
     *     3.   ,     3.   ,     3.   ,    10.   ,   200.   ,    11.   ,
     *   200.   ,   590.   ,   500.   ,   460.   ,   280.   ,   320.   ,
     *    60.   ,   230.   ,    18.   ,   380.   ,   240.   ,   350.   ,
     *   300.   ,   310.   ,   280.   ,   320.   ,   220.   ,   200.   ,
     *   200.   ,   220.   ,   100.   ,   100.   ,     3.   ,     3.   ,
     *     3.   ,     3.   ,     3.   ,     3.   ,     3.   ,     4.   ,
     *     4.   ,     4.   ,     5.   ,   500.   ,    40.   ,   550.   ,
     *   440.   ,   420.   ,   240.   ,   220.   ,   220.   ,   380.   ,
     *   350.   ,   400.   ,   380.   ,   370./
      data e17/
     *   290.   ,   290.   ,   310.   ,   300.   ,   320.   ,   230.   ,
     *   140.   ,   250.   ,   200.   ,    20.   ,   100.   ,    11.   ,
     *    14.   ,   200.   ,     4.   ,     5.   ,     4.   ,     4.   ,
     *    13.   ,    60.   ,    21.   ,   200.   ,   500.   ,   460.   ,
     *   280.   ,   320.   ,    60.   ,   240.   ,    19.   ,   380.   ,
     *   240.   ,   350.   ,   300.   ,   300.   ,   280.   ,   320.   ,
     *   230.   ,   210.   ,   200.   ,   200.   ,   110.   ,   200.   ,
     *     7.   ,    12.   ,     7.   ,     8.   ,     5.   ,     4.   ,
     *     4.   ,     4.   ,     4.   ,     4.   ,     6.   ,    19.   ,
     *    21.   ,    50.   ,   240.   ,   230.   ,   230.   ,   390.   ,
     *   350.   ,   400.   ,   380.   ,   420.   ,   290.   ,   290.   ,
     *   310.   ,   310.   ,   320.   ,   300.   ,   150.   ,   300.   ,
     *   200.   ,    16.   ,    50.   ,    17.   ,   100.   ,    16.   ,
     *     4.   ,     5.   ,     4.   ,     4.   ,     4.   ,    50.   ,
     *    16.   ,   200.   ,    16.   ,   300.   ,   330.   ,    60.   ,
     *   240.   ,    19.   ,   380.   ,   240.   ,   350.   ,   300.   ,
     *   300.   ,   280.   ,   320.   ,   230./
      data e18/
     *   210.   ,   200.   ,   280.   ,   100.   ,    70.   ,   200.   ,
     *   100.   ,    25.   ,    50.   ,     5.   ,     5.   ,     4.   ,
     *     4.   ,     4.   ,     4.   ,     4.   ,     5.   ,     4.   ,
     *     7.   ,    21.   ,   150.   ,   380.   ,   490.   ,   460.   ,
     *   500.   ,   430.   ,   430.   ,   350.   ,   290.   ,   370.   ,
     *   370.   ,   360.   ,   320.   ,   310.   ,   710.   ,   210.   ,
     *   190.   ,   150.   ,   140.   ,    60.   ,    80.   ,   220.   ,
     *     9.   ,    16.   ,    17.   ,     4.   ,     4.   ,     4.   ,
     *     4.   ,     6.   ,     5.   ,    12.   ,    13.   ,   380.   ,
     *   240.   ,   350.   ,   300.   ,   300.   ,   280.   ,   320.   ,
     *   230.   ,   220.   ,   200.   ,   320.   ,   300.   ,   300.   ,
     *   200.   ,   210.   ,   130.   ,    80.   ,   100.   ,    40.   ,
     *    11.   ,     9.   ,     4.   ,     4.   ,     4.   ,     4.   ,
     *     4.   ,     4.   ,     4.   ,     4.   ,     6.   ,   100.   ,
     *     2.8  ,   440.   ,   350.   ,   310.   ,   370.   ,   380.   ,
     *   370.   ,   330.   ,   320.   ,   700.   ,   180.   ,   160.   ,
     *   120.   ,   100.   ,    60.   ,    70./
      data e19/
     *    40.   ,    40.   ,     8.   ,    11.   ,     5.   ,     5.   ,
     *     4.   ,     4.   ,     6.   ,     5.   ,    11.   ,    12.   ,
     *   100.   ,   100.   ,   280.   ,   330.   ,   230.   ,   220.   ,
     *   200.   ,   320.   ,   300.   ,   300.   ,   200.   ,   210.   ,
     *   130.   ,    80.   ,   100.   ,    30.   ,    11.   ,     9.   ,
     *     5.   ,     5.   ,     4.   ,     4.   ,     4.   ,     5.   ,
     *     4.   ,     4.   ,     6.   ,   100.   ,     2.8  ,   300.   ,
     *   420.   ,   380.   ,   380.   ,   340.   ,   330.   ,   700.   ,
     *   180.   ,   160.   ,   120.   ,    90.   ,    60.   ,    70.   ,
     *    40.   ,    40.   ,     8.   ,    11.   ,     5.   ,     5.   ,
     *    13.   ,     6.   ,     7.   ,     6.   ,    11.   ,    13.   ,
     *    80.   ,   100.   ,   230.   ,   200.   ,   330.   ,   300.   ,
     *   300.   ,   200.   ,   220.   ,   130.   ,    80.   ,   100.   ,
     *    30.   ,    11.   ,    10.   ,     6.   ,     9.   ,    11.   ,
     *     9.   ,    11.   ,     6.   ,     5.   ,     4.   ,     6.   ,
     *   100.   ,     2.8  ,   380.   ,   350.   ,   330.   ,   690.   ,
     *   170.   ,   160.   ,   110.   ,    70./
      data e20/
     *    60.   ,    50.   ,    40.   ,    40.   ,     9.   ,    12.   ,
     *    10.   ,    13.   ,    15.   ,     6.   ,     8.   ,     7.   ,
     *    11.   ,    50.   ,     4.   ,    50.   ,    90.   ,   140.   ,
     *    90.   ,   200.   ,   280.   ,   340.   ,   300.   ,   310.   ,
     *   200.   ,   220.   ,   130.   ,    80.   ,   100.   ,    30.   ,
     *    12.   ,    11.   ,    10.   ,    12.   ,    14.   ,    14.   ,
     *    15.   ,     8.   ,     6.   ,     4.   ,     6.   ,     3.   ,
     *     2.8  ,     2.8  ,     4.   ,   160.   ,   280.   ,   180.   ,
     *   170.   ,   120.   ,    90.   ,    80.   ,    70.   ,    60.   ,
     *    40.   ,    13.   ,    50.   ,    50.   ,    50.   ,    50.   ,
     *     7.   ,     8.   ,     7.   ,    11.   ,     3.   ,     3.   ,
     *     4.   ,   150.   ,   200.   ,   100.   ,   200.   ,   200.   ,
     *   230.   ,   130.   ,   100.   ,   100.   ,    40.   ,    16.   ,
     *    50.   ,    23.   ,    13.   ,    16.   ,    17.   ,    18.   ,
     *     9.   ,     6.   ,     4.   ,     6.   ,     3.   ,     2.4  ,
     *     2.4  ,     2.1  ,     2.1  ,     6.   ,    50.   ,   130.   ,
     *   100.   ,    90.   ,    90.   ,   210./
      data e21/
     *   210.   ,   200.   ,    70.   ,    70.   ,    70.   ,    70.   ,
     *    12.   ,    10.   ,     7.   ,    12.   ,     3.   ,     2.8  ,
     *    11.   ,     2.4  ,     8.   ,   100.   ,   200.   ,    50.   ,
     *   300.   ,    30.   ,   100.   ,    21.   ,    10.   ,     6.   ,
     *    50.   ,     6.   ,     3.   ,     2.2  ,     2.2  ,     2.1  ,
     *     2.2  ,     2.1  ,     2.1  ,     5.   ,    90.   ,    50.   ,
     *    50.   ,   100.   ,   100.   ,     9.   ,     2.4  ,    50.   ,
     *     2.1  ,     2.1  ,     3.   ,    60.   ,   100.   ,   200.   ,
     *    11.   ,   140.   ,    23.   ,    50.   ,     8.   ,    60.   ,
     *     6.   ,     6.   ,     2.2  ,     2.2  ,     2.1  ,     2.1  ,
     *     2.1  ,     3.   ,     5.   ,    14.   ,    19.   ,   220.   ,
     *   210.   ,   210.   ,   140.   ,   140.   ,    50.   ,     5.   ,
     *    12.   ,     2.1  ,     2.2  ,     3.   ,     2.3  ,     3.   ,
     *    22.   ,   100.   ,   200.   ,   250.   ,   200.   ,   210.   ,
     *    30.   ,   100.   ,     6.   ,     6.   ,     2.2  ,     2.4  ,
     *     2.1  ,     2.2  ,     2.8  ,     5.   ,     5.   ,     5.   ,
     *    11.   ,   200.   ,   300.   ,   290./
      data e22/
     *   290.   ,   200.   ,   200.   ,   200.   ,     6.   ,    50.   ,
     *     5.   ,   100.   ,     6.   ,    21.   ,     3.   ,     7.   ,
     *   200.   ,   200.   ,   320.   ,   200.   ,   290.   ,    30.   ,
     *   140.   ,     6.   ,     6.   ,     2.5  ,   100.   ,     6.   ,
     *     2.3  ,     2.9  ,     5.   ,     5.   ,     7.   ,    12.   ,
     *   360.   ,   350.   ,   290.   ,   200.   ,   210.   ,   220.   ,
     *    50.   ,   110.   ,    50.   ,   100.   ,     7.   ,    21.   ,
     *     3.   ,     8.   ,   200.   ,   220.   ,   320.   ,   280.   ,
     *   350.   ,    40.   ,   170.   ,    21.   ,   100.   ,    20.   ,
     *   100.   ,    21.   ,     5.   ,     6.   ,     5.   ,     7.   ,
     *     7.   ,   350.   ,   280.   ,   290.   ,   300.   ,   210.   ,
     *   230.   ,   210.   ,   140.   ,     8.   ,    29.   ,   200.   ,
     *   200.   ,   180.   ,    26.   ,   220.   ,    26.   ,   100.   ,
     *    26.   ,    30.   ,   200.   ,    11.   ,   290.   ,   360.   ,
     *   290.   ,   240.   ,   210.   ,   140.   ,    50.   ,    60.   ,
     *   270.   ,   200.   ,   230.   ,   200.   ,   110.   ,   200.   ,
     *   200.   ,   300.   ,   410.   ,   350./
      data e23/
     *   240.   ,   220.   ,   150.   ,   450.   ,   460.   ,   360.   ,
     *   120.   ,   420.   ,   610.   ,   560.   ,   920.   ,   920.   ,
     *   420.   ,     0.   ,     0.   ,     0.   ,     0.   ,     0.   ,
     *     0.   ,     0.   ,     0.   ,     0.   ,     0.   ,     0.   ,
     *     0.   ,     0.   ,     0.   ,     0.   ,     0.   ,     0.   ,
     *     0.   ,     0.   ,     0.   ,     0.   ,     0.   ,     0.   ,
     *     0.   ,     0.   ,     0.   ,     0.   ,     0.   ,     0.   ,
     *     0.   ,     0.   ,     0.   ,     0.   ,     0.   ,     0.   ,
     *     0.   ,     0.   ,     0.   ,     0.   ,     0.   ,     0.   ,
     *     0.   ,     0.   ,     0.   ,     0.   ,     0.   ,     0.   ,
     *     0.   ,     0.   ,     0.   ,     0.   ,     0.   ,     0.   ,
     *     0.   ,     0.   ,     0.   ,     0.   ,     0.   ,     0.   ,
     *     0.   ,     0.   ,     0.   ,     0.   ,     0.   ,     0.   ,
     *     0.   ,     0.   ,     0.   ,     0.   ,     0.   ,     0.   ,
     *     0.   ,     0.   ,     0.   ,     0.   ,     0.   ,     0.   ,
     *     0.   ,     0.   ,     0.   ,     0.   ,     0.   ,     0.   ,
     *     0.   ,     0.   ,     0.   ,     0./
      nzno=mz*1000+n
      nl=((mz/10)*2-4)*100
      if(nl.le.0) nl=1
      do 11 j=nl,nt
      if(itable(j).eq.nzno) go to 12
   11 continue
      iret=iret+1
      return
   12 exces=amass(j)
      errr=err(j)
      return
      end
      subroutine tlj(h4,n1,j1,h1,z,w1,stplth)
      REAL*8    sg(2),s(16,4),f(101,2),g(101,3),r(2,3),h(502),a(502),b(5
     102),q(502),d(999),e(999),y(999),u(999),u1(7),y1(7)
      common /lab2/t(4,81),v(80),v1(3)
cneadb
      REAL*8 t1,t2,t3,t4,t5,t6,t7,t8,t9
cneadb
      b1 = 0.
      i1=2.0*v1(1) + 2
      i3=i1+2
      i4=i1*(i1-2)
      j1=i1-1
      b7=v(1)
      v(1)=abs(v(1))
      rb=v(1)+9.*v(2)
      rb=amax1(rb,v(10)+9.0*v(11))
      t1=v(3)+9.*v(4)
       if(rb.lt.t1)rb=t1
      continue
      irb=10.*rb
      rb=irb
      rb=rb/10.
      m4=rb/stplth
      if(m4.gt.299)m4=299
      continue
      d1=rb/float(m4)
      m4=m4+1
      m1=m4+1
      m3=m4+2
      n1=(v1(3)+1.0)
      v1(3)=float(n1-1)
      z2=z*z
      dt=2.*d1
      dtt=d1*d1/12.
      w2=w1
      w1=w1*dtt
      dt3=-2.*dtt
      h2=1.+dtt*h4
      hy=dt3*h1*z
      s(16,1)=1.
      s(16,2)=1.
      s(16,3)=0.
      s(16,4)=0.
      t1=1.
      t2=2.
      t3=0.
      do 30 i=1,15
      j=16-i
      t7=t1*z/t2
      t6=-t3*(t3+1.)+z2
      t4=t6/t2
      t5=(2.+t6)/t2
      s(j,1)=t7*s(j+1,1)-t4*s(j+1,3)
      s(j,2)=t7*s(j+1,2)-t5*s(j+1,4)
      s(j,3)=t7*s(j+1,3)+t4*s(j+1,1)
      s(j,4)=t7*s(j+1,4)+t5*s(j+1,2)
      t1=t1+2.
      t2=t2+2.
   30 t3=t3+1.
      tm=rb
      if(z.eq.0.)go to 50
      do 45 i=1,4
cneadbt1=(exp((13.816+alog(abs(s(1,i))))/15.0))/h1
      t1=(exp((13.816+log(abs(s(1,i))))/15.0))/h1
        if((t1-tm).gt.0.)tm=t1
      continue
   45 continue
   50 m2=(tm-rb)/d1 + 1
      tm=rb+d1*float(m2)
      t1=z2+16.
cneadbsg(1)=-z+z*(alog(t1))/2.0+3.5*atan(z/4.0)-atan(z)-atan(z/2.0)-atan
      sg(1)=-z+z*(log(t1))/2.0+3.5*atan(z/4.0)-atan(z)-atan(z/2.0)-atan
     1(z/3.)-z*(1.+(z2-48.)/(30.*t1*t1)+(z2*z2-160.*z2+1280.)/(105.*t1*
     2t1*t1*t1))/(12.*t1)
      sg(2)=sg(1)-1.5707963+atan(z)
      t1=tm+dt
      do 65 i=1,2
      t1=t1-d1
      t2=t1*h1
cneadbt3=t2-z*alog(2.0*t2)
      t3=t2-z*log(2.0*t2)
      do 60 j=1,2
      t7=0.
      t8=0.
      do 55 k=1,15
      t7=(t7+s(k,j))/t2
   55 t8=(t8+s(k,j+2))/t2
      t7=t7+1.
      t4=t3+sg(j)
      t5=cos(t4)
      t6=sin(t4)
   60 g(j,i)=t7*t5-t8*t6
      r(1,i)=h2+hy/t1
   65 r(2,i)=r(1,i)+dt3/(t1*t1)
      do 70 i=1,m2
      t1=t1-d1
      r(1,3)=h2+hy/t1
      r(2,3)=r(1,3)+dt3/(t1*t1)
      do 70 j=1,2
      g(j,3)=((12.-10.*r(j,2))*g(j,2)-r(j,1)*g(j,1))/r(j,3)
      do 70 k=1,2
      r(j,k)=r(j,k+1)
   70 g(j,k)=g(j,k+1)
      t1=t1+d1
      do 130 i=1,2
      t2=t1*h1
      t3=1./t2
      t1=t1-d1
      t7=0.
      n4=0
      n2=(t2/1.4142)*sqrt(25.-2.*z*t3+10.*sqrt((z*t3-.5)**2+6.))
      if((n2-n1-8).lt.0.)n2=n1+8
      continue
   80 if ((n2-500).gt.0.)go to 230 
      t6=float(n2)
      n3=n2+1
      t5=t7
      h(n3+1)=0.
      h(n3)=1.0e-20
      a(n3)=sqrt(z2+(t6+1.0)**2)/(t6+1.0)
      do 105 k=1,n2
      m=n3-k
      if((n4-m).ge.0.)go to 95
      a(m)=sqrt(z2+t6*t6)/t6
      b(m)=(2.*t6+1.)*(z/(t6*(t6+1.))+t3)
      t6=t6-1.
   95 h(m)=(b(m)*h(m+1)-a(m+1)*h(m+2))/a(m)
      if((abs(h(m+2))-(10.0**30)).le.0.)go to 105 
      h(m)=h(m)/(10.0**25)
      h(m+1)=h(m+1)/(10.0**25)
  105 continue
      n4=n2
      n2=n2+10
      t7=h(2)/h(1)
      if((abs((t5-t7)/t7)-0.0001).le.0.)go to 115 
      if(n2.le.500) go to 80
      n2 = 500
  115 t5=1./(a(1)*(h(1)*g(2,i)-h(2)*g(1,i)))
      f(1,i)=t5*h(1)
      f(2,i)=t5*h(2)
      do 130 k=3,n1
      if((abs(h(k)/h(k-1))-10.0**15).le.0.)go to 125 
      t5=t5/(10.0**25)
  125 f(k,i)=t5*h(k)
  130 g(k,i)=(b(k-2)*g(k-1,i)-a(k-2)*g(k-2,i))/a(k-1)
      q(1)=0.
      q(2)=0.
      h(1)=0.
      h(2)=0.
      a(1)=0.
      a(2)=1.0e-20
      b(1)=0.
      b(2)=1.0e-20
      hz3=h2+1.5*hy/v(9)
      hz=-hy/(2.*(v(9)**3))
      t101=1.0/exp(v(10)/v(11))
      t1=1.0/exp(v(1)/v(2))
      t102=exp(d1/v(11))
      t2=exp(d1/v(2))
       if(b7.ge.0.)go to 140
      a2=v(5)*v(5)/16.
      v(5)=1.
      a1=4.*w2*a2
      b1=exp(4.0*a2*h4)
  140 v5=v(5)
      if(v5.ge.0.)go to 150
      v(5)=-v5
      v4=v(4)
      v(4)=0.69*v(4)
  150 t9=v(8)*(1.-v(5))
      t3=1.0/exp(v(3)/v(4))
      t4=exp(d1/v(4))
      t6=0.
      do 200 k=1,m4
      t6=t6+d1
      y(k)=t6*t6
      t1=t1*t2
      t101=t101*t102
      t3=t3*t4
      t5=v(7)/(1.+t1)
      t7=1./(1.+t3)
      if(v5.ge.0.)go to 170
      t8=((t6-v(3))/v4)**2
      if((t8-10.).lt.0.)go to 165
      t8=0.
      go to 185
  165 t8=v(8)*exp(-t8)
      go to 185
  170 t8=4.*t3*t7*t7*v(8)
      if(b7.ge.0.)go to 185
      p5=a2/(t5*t5+t8*t8)
      p8=t1/(1.+t1)
      p9=t3*t7
      p6=-t5*p8/v(2)
      p7=t8*(1.-2.*p9)/v(4)
      p8=p6*(1.-2.*p8)/v(2)
      p9=t8*(1.-6.*p9*(1.-p9))/(v(4)*v(4))
      u2=p5*((t5*p6+t8*p7)*2./t6+t5*p8+t8*p9)
      y2=p5*((t5*p7-t8*p6)*2./t6+t5*p9-t8*p8)
      y1(1)=t8/(b1*b1+2.*b1*a1*t5)
      u1(1)=(t5+t8*a1*y1(1))/(b1+a1*t5)
      do 180 j=1,6
      p5=a1*y1(j)-y2
      p6=sin(p5)
      p5=cos(p5)
      p7=1.0/(b1*exp(a1*u1(j)-u2))
      u1(j+1)=(t5*p5+t8*p6)*p7
  180 y1(j+1)=(t8*p5-t5*p6)*p7
      t5=u1(7)-((u1(7)-u1(6))**2)/(u1(7)-2.*u1(6)+u1(5))
      t8=y1(7)-((y1(7)-y1(6))**2)/(y1(7)-2.*y1(6)+y1(5))
  185 h(k+2)=w1*(t9*t7+v(5)*t8)
      u(k)=2.*w1*v(6)*t101/(v(11)*t6*((1.+t101)**2))
      if((t6-v(9)).ge.0.)go to 195
      e(k)=w1*t5+hz3+hz*y(k)
      go to 200
  195 e(k)=w1*t5+h2+hy/t6
  200 continue
      t8=0.
      do 225 i=1,n1
      i2=2*i
      t1=dtt*t8*(t8+1.)
      t8=t8+1.
      do 205 k=1,m4
  205 d(k)=e(k)-t1/y(k)
      do 225 j=1,j1
      l=i2-i3+2*j
      if((iabs(i2-i1)-l).gt.0.)then
      t(j,i)=0.
       else
      t9=(float(l*(l+2)-i2*(i2-2)-i4))/4.
      do 220 k=1,m4
      q(k+2)=d(k)+t9*u(k)
      t3=12.-10.*q(k+1)
      t4=10.*h(k+1)
      t1=t3*a(k+1)-q(k)*a(k)+t4*b(k+1)+h(k)*b(k)
      t2=t3*b(k+1)-q(k)*b(k)-t4*a(k+1)-h(k)*a(k)
      t3=q(k+2)**2+h(k+2)**2
      a(k+2)=(q(k+2)*t1+h(k+2)*t2)/t3
  220 b(k+2)=(q(k+2)*t2-h(k+2)*t1)/t3
      t3=a(m3)**2+b(m3)**2
      t1=(a(m1)*a(m3)+b(m1)*b(m3))/t3
      t2=(a(m3)*b(m1)-a(m1)*b(m3))/t3
      t5=f(i,2)-f(i,1)*t1
      t6=f(i,1)*t2
      t3=t5-g(i,1)*t2
      t4=g(i,2)-g(i,1)*t1+t6
      t7 = (t3/t4)*t3+ t4
      t1=(t3*t5+t4*t6)/t7
      t2=(t4*t5-t3*t6)/t7
      t7 = t4
      t1 = t1/t7
      if(t1.lt..001)t1=0.
      t2 = t2/t7
      t90=abs(t2)
      if(t90.lt..001)t2=0.
      t(j,i)=4.*(t1-t1**2-t2**2)
      if(t(j,i).lt.0.)t(j,i)=0.
                               endif
  225 continue
  230 continue
      return
      end
c
      subroutine fiznorm
      common/pl3/lz,nepr,eb(100),rcsp(100),zee,amass,plex
      common/fizzc/popfiz(13,20,150),pofizbig(13,20,150,30)
c this subroutine is to convert the file of fizzer populations
c from mb per unit to probability
       cstot=0.
       cjtot=0.
      do 5 iz=1,13
      do 4 ia=1,20
      do 3 ie=1,150
      cstot=cstot+popfiz(iz,ia,ie)
      do 2 jf=1,30
      cjtot=cjtot+pofizbig(iz,ia,ie,jf)
   2  continue
   3  continue
   4  continue
   5  continue
      if(cjtot.le.0..or.cstot.le.0.)go to 108
      do 8 iz=1,13
      do 7 ia=1,20
      do 6 ie=1,150
      popfiz(iz,ia,ie)=popfiz(iz,ia,ie)/cstot
      do 12 jf=1,30
      pofizbig(iz,ia,ie,jf)=pofizbig(iz,ia,ie,jf)/cjtot
  12  continue
   6  continue
   7  continue
   8  continue
108   continue
      write(21,*)' Z A E/2   Nf(ZAF) unity normalized'
      write(22,*)' Z  A  J E/2   Nfj(ZAJE) unity normalized'
      do 11 iz=1,13
      ize=zee+1-iz
      do 10 ia=1,20
      iae=amass+2-ia-iz
      do 9 ne=1,141,10
      ss=0.
      do 99 ie=ne,ne+9
      ss=ss+popfiz(iz,ia,ie)
 99   continue
       if(ss.gt.0.)then
      write(21,20)ize,iae,ne,(popfiz(iz,ia,ie),ie=ne,ne+9)
       endif
      do 13 jf=1,30
      s=0.
      do 22 ie=ne,ne+9
22    s=s+pofizbig(iz,ia,ie,jf)
c     write(22,*)'jfiss= ',jf
      if(s.gt.0.)then
      write(22,23)ize,iae,jf,ne,(pofizbig(iz,ia,ie,jf),ie=ne,ne+9)
      endif
  13  continue
   9  continue
  10  continue
  11  continue
  20  format(2i3,1x,i3,10(1x,1pe8.2))
  23  format(4i3,10(1x,1pe8.2))
      return
      end
c
       subroutine plt
      common/astp3/niso,no,iadlt,isotag,massno(12),mass1(12),ncsave
      common/atsv/atsave
      dimension ma(22)
       common/scrn/att(10),cldd(10),bexp(24,15),reacss(100),edelt
       common/rcsrec/rcsorig(100)
       common/scrn1/ingo,neng,nengg,ireact
      common/shft/mc,mp,inver,ike,ipch,kq,qval,ap,at,zp,zt,cld,barfac
      common/fizzc/popfiz(13,20,150),pofizbig(13,20,150,30)
      common/gams/gamgam,csgam,csgamm(100)
      REAL*8 gamgam,csgamm,csgam
      common/outlim/limout,lim2
       common/fizz/fizbar(106,260),fizcs(100)
      common/eqin/eind(100)
      common/shft4/k3,jcal
      common/shft5/fiss
      common/shft2/fs(25),dsp(25),brr(25),der(25),er(25)
     1,cx(25),bilsum,xmiss,sumiz
      REAL*8 fiss
      common/sf/m3,kplt
      common/ismrdat2/ididl(106,260),nolevs(106,260)
      common/pl1/ecm(100)
      common/pl2/sux(24,13,100)
      common/nam/k2
      common/pl3/lz,nepr,eb(100),rcsp(100),zee,amass,plex
      common/ug/iz,ia
      common/pl5/na,nz,title(20)
      common/pl4/crs(24,13)
      common/bound/neut2,iprot,idnuc(100),idm(100),einc(100),
     1thett(100),phit(100),ppp
      common/isomer1/xgs(690),xiso(690),x2iso(690),isonum,isoint,iso2num
      common/isomer2/zs(690),as(690),eiso(690),eiso2(690)
      common/isomer3/iziso(690),iaiso(690),csgs(690),csis(690)
     1,csgse(690,100),csise(690,100),csis2e(690,100),rm1(690,100),
     2rm2(690,100)
 
      common/sigsav/sigdr(100),sigqd(100)
      common/ismrdat/isoid(106,260),levnm(1000),bevj(1000,3)
C
      common/testout/sumcs(50)
c     switch=1.
       go to(5,15,25),k2
    5 continue
c       if(no.eq.1)then
        if(no.eq.0)then
      do 10 iz=1,13
      do 10 ia=1,24
      do 10 nep = 1,100
   10 sux(ia,iz,nep) = 0.0
                   endif
      return
   15 do 20 im=1,na
      do 20 iz=1,nz
      sux(im,iz,nepr)=crs(im,iz)
   20 continue
      return
   25 nzl=1
      write(17,901)
      write(17,900)ppp
 900  format(' out of bounds cross section(mb)= ',1pe9.2)
 901  format(' ')
      nzpl=nz
      bl=0.
   30 if(nzl.gt.plex)return
c 10/2/2014 k2 cond
      mz=zee-bl
c     ma(1)=amass-bl
      ma(1)=atsave+ap-bl
      do 35 ik=2,22
   35 ma(ik)=ma(ik-1)-1
      write(17,80)title
      write(17,70)
      write(17,60)mz,mz,mz,mz,mz,mz,mz,mz,mz,mz,mz
      write(17,65)(ma(i),i=1,11)
c----------------------------
      do 55 n=1,nengg
c     if(rcsp(n).gt.0.)write(17,75)ecm(n),eb(n),rcsp(n),(sux(ia,nzl,n),
c    1ia=1,11)
   55 continue
c----------------------------
      if(na.lt.12)go to 555
      write(17,80)title
      write(17,70)
      write(17,60)mz,mz,mz,mz,mz,mz,mz,mz,mz,mz,mz
      write(17,65)(ma(i),i=12,22)
c     do 551 n=1,nengg
c     if(rcsp(n).gt.0.)write(17,75)ecm(n),eb(n),rcsp(n),(sux(ia,nzl,n),
c    1ia=12,22)
  551 continue
  555 continue
c 6/96 add device 7 output sum
c
c feb26'03 add call to sub here to normalize fissioning nuclei to
c unity from mb
c--------------------------------------------------------------------
c     if(switch.ne.1.)go to 654
c try nzl fix
        if(k2.le.2)return
      masss=atsave+ap
      mz=zee
      write(7,901)
c     write(7,900)ppp
      write(7,901)
      write(7,80)title
c23456789012345678901234567890123456789012345678901234567890123456789012
      write(7,*)' Ebeam  Rx.cs used  Rx. cs algo fiss cs (mb) capture cs
     1 (mb)'
      do 5678 nen=1,nengg
      i=nen
      if(sux(1,1,i).eq.0..and.csgamm(i).gt.0.)sux(1,1,i)=csgamm(i)
      write(7,5679)eind(nen),reacss(nen),rcsorig(nen),fizcs(nen),
     1sux(1,1,nen)
 5678 continue
 5679 format(1x,1f6.2,2x,1pe9.3,3(4x,1pe8.2))
 
c-------------------------------------------------
      if(ap.eq.0.)then
      write(7,8766)
8765  format(6(1x,1pe10.3))
      do 8768 i=1,nengg
8768  write(7,8765)eind(i),rcsp(i),sigdr(i),sigqd(i),fizcs(i),
     1sux(1,1,i)
      else
      write(7,6766)
      endif
c-----------------------------------------------
8766  format(' Ebeam (MeV) Rx cs (mb) GDR (mb)   QD (mb)  fiss cs (mb) capture cs (mb)')
     1apture cs (mb)')
6766  format(' Ebeam (MeV) Rx cs (mb) fiss cs (mb) capture cs (mb)')
      write(7,*)
      write(*,*)
      write(7,*)' Ebeam  Z   A   Total    grnd st  isomer1 isomer2  isomer
     1er1  isomer2  didl nolevs'
      write(7,*)' (MeV)           cs.      cs.      cs.     cs.     to tot
     1otal to total'
c--------------------------------------
C
       do nen=1,nengg
         sumcs(nen)=0.
           enddo
c--------------------------------------
      do 650 iz=1,nz
      mzf=mz+1-iz
c-----------------------------------
      do 648 ia=1,na
      masout=masss+2-ia-iz
      n=isoid(mzf,masout)
      isno=n
      nolev=levnm(n)
      continue
c----------------------------------
      do 646 nen=1,nengg
      if(isno.le.1)go to 9646
      if(sux(ia,iz,nen).gt.0.)then
      rm1(isno,nen)=csise(isno,nen)/sux(ia,iz,nen)
      rm2(isno,nen)=csis2e(isno,nen)/sux(ia,iz,nen)
      else
      rm1(isno,nen)=0.
      rm2(isno,nen)=0.
      endif
c--------------------------------------
C
c--------------------------------------
        nolev=levnm(isno)
      if(nolev.gt.2.and.sux(ia,iz,nen).gt.0.)then
      write(7,6721)eb(nen),mzf,masout,sux(ia,iz,nen),csgse(isno
     1,nen),csise(isno,nen),csis2e(isno,nen),rm1(isno,nen),rm2(isno,nen)
     2,ididl(mzf,masout),nolev
      endif
      if(sux(ia,iz,nen).gt.0..and.nolev.eq.2)then
      write(7,6521)eb(nen),mzf,masout,sux(ia,iz,nen),csgse(isno
     1,nen),csise(isno,nen),rm1(isno,nen),ididl(mzf,masout),nolev
       go to 6461
      endif
9646  continue
       if(nolev.eq.1)then
      if(sux(ia,iz,nen).gt.0.)write(7,652)eb(nen),mzf,masout,sux(ia,
     1iz,nen),nolev
       endif
c        endif
6461  continue
       sumcs(nen)=sumcs(nen)+sux(ia,iz,nen)
c        if(sux(ia,iz,nen).gt.0.)then
c      write(7,*)'nen=',nen
c     write(7,*)'z=',mzf,'a=',masout,'nen=',nen,sux(ia,iz,nen),'ia=',ia
c         endif
  646 continue
  648 continue
  650 continue
C
C
         do ien=1,nengg
       write(7,*)' Ebeam =',eb(ien),'sumcs=',sumcs(ien)
           enddo
       continue
  654 continue
c--------------------------------------------------------------
      switch=0.
  652 format(1f6.2,1x,1i3,1x,1i3,1x,1pe8.2,52x,1i1)
 6721 format(1f6.2,1x,1i3,1x,1i3,6(1x,1pe8.2),3x,1i1,3x,1i1)
 6521 format(1f6.2,1x,1i3,1x,1i3,3(1x,1pe8.2),10x,(1pe8.2),12x,1i1,3x
     1,1i1)
      nzl1 = nzl -1
      nzl = nzl + 1
      bl=bl+1.
c     go to 30
   60 format('  exc  elab  rcs ',11(3x,i3,4x))
   65 format(' (mev) (mev) (mb)',11(3x,i3,4x))
   70 format(' excitation function data ')
   75 format(1x,f5.1,1x,f5.1,1x,f6.1,11(1x,f8.3,1x))
   80 format(1h1,20a4)
      return
      end
 
      subroutine barfit1(iz,ia,il,ibfis,irot,ilmax)

c ======================================================================
c
c    This program is for testing the subroutines momfit and barfit to see
c    that barriers and moments of inertia are correctly calculated.
c
c        Written by A. J. Sierk  LANL   T-9   February, 1984.
c
c ======================================================================

      implicit real*8 (a-h, o-z)
      implicit integer (i-n)

      real*8 elmax
c ======================================================================
c
c  tape 6 is standard output on a SUN workstation.
c
c  1  write (6, 100)
c
c  tape 5 is standard input on a SUN workstation.
c
c     read (5, *, err=1) iswitch
c  2  write (6, 120)
c     read (5, *, err=2) iz, ia
c     if (iz.eq.-2) go to 99
        el=dble(il)
      zn = dble (iz)
      an = dble (ia)
      call elmaxc (zn, an, elmax)
c     write(*,*)'barfit el=',el
       ilmax=10000.*elmax
c     open (7, file='bar.out')
c     write (7, 130) iz, ia, elmax
c     write (6, 130) iz, ia, elmax
c     if (iswitch.eq.1) go to 2
c  5  if (iswitch.eq.2) write (6, 140)
c     read (5, *, err=5) el0, eldel, elmax
c     if (eldel.eq.0.d0) go to 2
      nmax = 0
c         eldel=1.
c         if(el.gt.elmax)el=elmax
c     write (6, 150)
c     write (7, 150)
c     if (eldel.ne.0.d0) nmax = int((elmax - el0)/eldel) + 1
c       do i = 1,nmax
c       el = eldel*(i-1) + el0
        if (el.gt.elmax) el = elmax
c       call momfit (iz, ia, el, aimin, aimid, aimax, elmax)
        call barfit (iz, ia, el, bfis, egs, elmax)
         irot=10000.*egs
          ibfis=10000.*bfis
c       write (6, 160) el, bfis, egs, aimin, aimid, aimax
c       write (7, 160) el, bfis, egs, aimin, aimid, aimax
        if (el.eq.elmax)return 
c       end do 
c     call momfit (iz, ia, 1000.d0, aimin, aimid, aimax, elmax)
c     write (7, 170) elmax, aimin, aimid, aimax
c     write (6, 170) elmax, aimin, aimid, aimax
c  20 go to 2
   99 return

c ======================================================================

  100 format (/2x,'Enter 1 for finding only Lmax; '/2x,'Enter 2 ',
     & 'for table(s) of GS energies, barriers and moments vs. L:')
  120 format (2x,'Enter iz and ia; iz = -2 to exit program:')
  130 format (/5x,'Z = ',i3,', A = ',i3,', Lmax = ',f6.1)
  140 format (2x,'Enter the range of L values: L0, delta L and Lmax:')
  150 format (/7x,'L',8x,'Bf',6x,'Egs',8x,' Imin',4x,' Imid',4x,
     &        ' Imax'/)
  160 format (4x,f6.2,2f9.2,4x,3f9.3)
  170 format (/6x,'Lmax, Moments of inertia at Lmax'//5x,f5.1,2x,3f8.3)

c ======================================================================
      end

      subroutine barfit (iz, ia, el, bfis, egs, elmax)
c     subroutine barfit (iz, ia, il, bfis, egs, elmax)

c ======================================================================
c
c    This subroutine returns the barrier height bfis, the ground-state
c    energy egs, in MeV, and the angular momentum at which the fission
c    barrier disappears,  Lmax,  in units of h-bar,
c    when called with integer arguments iz, the atomic number,
c    ia, the atomic mass number, and el, the angular momentum in units
c    of h-bar, (Planck's constant divided by 2*pi).
c
c         The fission barrier for  el = 0  is calculated from a 7th order
c    fit in two variables to 638 calculated fission barriers for z values
c    from 20 to 110.  These  638 barriers are fit with an rms deviation of
c    0.10 MeV by this 49-parameter function.
c    If  barfit  is called with (iz,ia) values outside the range of the fit
c    the barrier height is set to 0.0, and a message is
c    printed on the default output file.
c
c         For el values not equal to zero, the values of
c    L at which the barrier is  80%  and  20%  of the L=0 value are
c    respectively fit to 20-parameter functions of  Z  and  A, over a more
c    restricted range of  A  values, than is the case for  L = 0.
c    The value of L where the barrier disappears, Lmax, for 61 nuclei,
c    is fit to a 35-parameter function of Z and A,  with the same range of
c    Z  and  A  values as  l-80  and  l-20.
c         Once again, if an  (iz,ia) pair is outside of the range of
c    validity of the fit, the barrier value is set to  0.0  and a message
c    is printed.  These three values  (Bfis(L=0),L-80, and L-20) and the
c    constraints  of  Bfis = 0 and  d(Bfis)/dL = 0 at L = Lmax and L = 0
c    lead to a fifth-order fit to Bfis(L) for L> L-20.  The first three
c    constraints lead to a third-order fit for the region L < L-20.
c
c         The ground-state energies are calculated from a 175-parameter
c    fit in Z, A, and L to 329 ground-state energies for 36 different
c    Z  and  A  values.
c    (The range of Z and A is the same as for L-80, L-20, and L-max)
c
c         The calculated barriers from which the fits were
c    made were calculated in 1983-1985 by A. J. Sierk of Los Alamos
c    National Laboratory   Group T-9, using  Yukawa-plus-exponential double
c    folded nuclear energy, exact Couloub diffuseness corrections,
c    and diffuse-matter moments of inertia. The parameters of the model
c    are those derived by Moller and Nix in 1979:
c    r-0 = 1.16 fm, as = 21.13 MeV, kappa-s = 2.3  a = 0.68 fm.
c    The diffuseness of the matter and charge distributions used
c    corresponds to a surface diffuseness parameter (defined by Myers)
c    of 0.99 fm.  The calculated barriers for L = 0 are
c    accurate to a little less than 0.1 MeV;  the output from this
c    subroutinle
c    is a little less accurate.  Worst errors may be as large
c    as 0.5 MeV; characteristic uncertainty is in the range of 0.1-0.2
c    MeV.   The values of egs are generally approximated to within
c    about 0.1-0.2 MeV;  the largest deviation is about 0.5 MeV,
c    near L-I for light nuclei.
c         The rms deviation of Lmax from the 61 input values is 0.31
c    h-bar.  The approximate value is nearly always within
c    0.5 h-bar of the calculated one.
c
c    Below is a table of test values to check implementation
c    of the program.
c    Z, A,  L    Egnd st  Fiss Bar      Moments of inertia     Lmax
c
c   28, 58, 0    0.00     33.14        0.816 3.603 3.603      46.1
c         ,25   21.61     19.50        0.778 3.662 3.662      46.1
c         ,40   50.12      2.97        0.723 3.647 3.648      46.1
c         ,45.8 59.14      0.00        0.746 3.160 3.160      46.1
c   65,153, 0    0.00     28.88        0.621 3.698 3.698      82.3
c         ,50   19.07     16.16        0.615 3.638 3.638      82.3
c         ,80   45.37      0.23        0.618 2.736 2.760      82.3
c         ,82.1 47.04      0.00        0.682 2.231 2.276      82.3
c   93,229, 0    0.00      3.76        0.715 1.747 1.747      68.1
c         ,45    8.22      1.26        0.766 1.577 1.577      68.1
c         ,68.1 17.96      0.00        1.053 1.053 1.236      68.1
c
c    written by A. J. Sierk,  LANL  T-9
c    Version 1.0   February, 1984
c    Version 1.1   January, 1985  Improved coefficients in egs and Lmax
c    Version 1.2   September, 1985  Improved Lmax, egs coefficients
c    Version 1.21  June, 1986   minor changes made
c    Version 1.3   February, 1996  Moved elmax to ELMAXC subroutine
c                  Improved Lmax coefficients.
c
c        Copyright, 1996,  The Regents of the University of California.
c        This software was produced under a U. S. Government contract
c        (W-7405-ENG-36) by the Los Alamos National Laboratory, which is
c        operated by the University of California for the U. S. Department
c        of Energy.  The U. S. Government is licensed to use, reproduce,
c        and distribute this software.  Permission is granted to the public
c        to copy and use this software without charge, provided that this
c        notice and any statement of authorship are reproduced on all
c        copies.  Neither the Government nor the University makes any
c        warranty, expressed or implied, or assumes any liability
c        or responsibility for the use of this software.
c
c ======================================================================

      implicit integer (i-n)
      implicit real*8 (a-h, o-z)

      dimension elzcof(7,7), elmcof(5,4), emncof(5,4), pa(10), pz(10), 
c     dimension elzcof(7,7), elmcof(5,4), emncof(5,4), pa(7), pz(7), 
     &          pl(10)
c     dimension emxcof(7,5),egs4(5,7),egs1(5,7),egs2(5,7),egs5(5,7),
c    & egs3(5,7), egscof(5,7,5)
      dimension egs4(5,7), egs1(5,7), egs2(5,7), egs5(5,7), egs3(5,7),
     &          egscof(5,7,5)

      equivalence (egs1,egscof), (egs2,egscof(1,1,2)),
     &            (egs3,egscof(1,1,3)), (egs4,egscof(1,1,4)),
     &            (egs5,egscof(1,1,5))

      data emncof
     &/-9.01100d+2,-1.40818d+3, 2.77000d+3,-7.06695d+2, 8.89867d+2,
     &  1.35355d+4,-2.03847d+4, 1.09384d+4,-4.86297d+3,-6.18603d+2,
     & -3.26367d+3, 1.62447d+3, 1.36856d+3, 1.31731d+3, 1.53372d+2,
     &  7.48863d+3,-1.21581d+4, 5.50281d+3,-1.33630d+3, 5.05367d-2/
      data elmcof
     & /1.84542d+3,-5.64002d+3, 5.66730d+3,-3.15150d+3, 9.54160d+2,
     & -2.24577d+3, 8.56133d+3,-9.67348d+3, 5.81744d+3,-1.86997d+3,
     &  2.79772d+3,-8.73073d+3, 9.19706d+3,-4.91900d+3, 1.37283d+3,
     & -3.01866d+1, 1.41161d+3,-2.85919d+3, 2.13016d+3,-6.49072d+2/
      data elzcof
     & /5.11819909d+5,-1.30303186d+6, 1.90119870d+6,-1.20628242d+6,
     &  5.68208488d+5, 5.48346483d+4,-2.45883052d+4,
     & -1.13269453d+6, 2.97764590d+6,-4.54326326d+6, 3.00464870d+6,
     & -1.44989274d+6,-1.02026610d+5, 6.27959815d+4,
     &  1.37543304d+6,-3.65808988d+6, 5.47798999d+6,-3.78109283d+6,
     &  1.84131765d+6, 1.53669695d+4,-6.96817834d+4,
     & -8.56559835d+5, 2.48872266d+6,-4.07349128d+6, 3.12835899d+6,
     & -1.62394090d+6, 1.19797378d+5, 4.25737058d+4,
     &  3.28723311d+5,-1.09892175d+6, 2.03997269d+6,-1.77185718d+6,
     &  9.96051545d+5,-1.53305699d+5,-1.12982954d+4,
     &  4.15850238d+4, 7.29653408d+4,-4.93776346d+5, 6.01254680d+5,
     & -4.01308292d+5, 9.65968391d+4,-3.49596027d+3,
     & -1.82751044d+5, 3.91386300d+5,-3.03639248d+5, 1.15782417d+5,
     & -4.24399280d+3,-6.11477247d+3, 3.66982647d+2/
      data egs1 /
     & -1.781665232d6,-2.849020290d6, 9.546305856d5, 2.453904278d5,
     &  3.656148926d5,
     &  4.358113622d6, 6.960182192d6,-2.381941132d6,-6.262569370d5,
     & -9.026606463d5,
     & -4.804291019d6,-7.666333374d6, 2.699742775d6, 7.415602390d5,
     &  1.006008724d6,
     &  3.505397297d6, 5.586825123d6,-2.024820713d6,-5.818008462d5,
     & -7.353683218d5,
     & -1.740990985d6,-2.759325148d6, 1.036253535d6, 3.035749715d5,
     &  3.606919356d5,
     &  5.492532874d5, 8.598827288d5,-3.399809581d5,-9.852362945d4,
     & -1.108872347d5,
     & -9.229576432d4,-1.431344258d5, 5.896521547d4, 1.772385043d4,
     &  1.845424227d4/
      data egs2 /
     &  4.679351387d6, 7.707630513d6,-2.718115276d6,-9.845252314d5,
     & -1.107173456d6,
     & -1.137635233d7,-1.870617878d7, 6.669154225d6, 2.413451470d6,
     &  2.691480439d6,
     &  1.237627138d7, 2.030222826d7,-7.334289876d6,-2.656357635d6,
     & -2.912593917d6,
     & -8.854155353d6,-1.446966194d7, 5.295832834d6, 1.909275233d6,
     &  2.048899787d6,
     &  4.290642787d6, 6.951223648d6,-2.601557110d6,-9.129731614d5,
     & -9.627344865d5,
     & -1.314924218d6,-2.095971932d6, 8.193066795d5, 2.716279969d5,
     &  2.823297853d5,
     &  2.131536582d5, 3.342907992d5,-1.365390745d5,-4.417841315d4,
     & -4.427025540d4/
      data egs3 /
     & -3.600471364d6,-5.805932202d6, 1.773029253d6, 4.064280430d5,
     &  7.419581557d5,
     &  8.829126250d6, 1.422377198d7,-4.473342834d6,-1.073350611d6,
     & -1.845960521d6,
     & -9.781712604d6,-1.575666314d7, 5.161226883d6, 1.341287330d6,
     &  2.083994843d6,
     &  7.182555931d6, 1.156915972d7,-3.941330542d6,-1.108259560d6,
     & -1.543982755d6,
     & -3.579820035d6,-5.740079339d6, 2.041827680d6, 5.981648181d5,
     &  7.629263278d5,
     &  1.122573403d6, 1.777161418d6,-6.714631146d5,-1.952833263d5,
     & -2.328129775d5,
     & -1.839672155d5,-2.871137706d5, 1.153532734d5, 3.423868607d4,
     &  3.738902942d4/
      data egs4 /
     &  2.421750735d6, 4.107929841d6,-1.302310290d6,-5.267906237d5,
     & -6.197966854d5,
     & -5.883394376d6,-9.964568970d6, 3.198405768d6, 1.293156541d6,
     &  1.506909314d6,
     &  6.387411818d6, 1.079547152d7,-3.517981421d6,-1.424705631d6,
     & -1.629099740d6,
     & -4.550695232d6,-7.665548805d6, 2.530844204d6, 1.021187317d6,
     &  1.141553709d6,
     &  2.182540324d6, 3.646532772d6,-1.228378318d6,-4.813626449d5,
     & -5.299974544d5,
     & -6.518758807d5,-1.070414288d6, 3.772592079d5, 1.372024952d5,
     &  1.505359294d5,
     &  9.952777968d4, 1.594230613d5,-6.029082719d4,-2.023689807d4,
     & -2.176008230d4/
      data egs5 /
     & -4.902668827d5,-8.089034293d5, 1.282510910d5,-1.704435174d4,
     &  8.876109934d4,
     &  1.231673941d6, 2.035989814d6,-3.727491110d5, 4.071377327d3,
     & -2.375344759d5,
     & -1.429330809d6,-2.376692769d6, 5.216954243d5, 7.268703575d4,
     &  3.008350125d5,
     &  1.114306796d6, 1.868800148d6,-4.718718351d5,-1.215904582d5,
     & -2.510379590d5,
     & -5.873353309d5,-9.903614817d5, 2.742543392d5, 9.055579135d4,
     &  1.364869036d5,
     &  1.895325584d5, 3.184776808d5,-9.500485442d4,-3.406036086d4,
     & -4.380685984d4,
     & -2.969272274d4,-4.916872669d4, 1.596305804d4, 5.741228836d3,
     &  6.669912421d3/

c ======================================================================
c
c    The program starts here

c---1/30/15 el=float(il)
c          el=dble(il)
      if (iz.lt.19 .or. iz.gt.111) go to 900
      if (iz.gt.102 .and. el.gt.0.d0) go to 910
      z = dble(iz)
      a = dble(ia)
      amin = 1.2d0*z + 0.01d0*z*z
      amax = 5.8d0*z - 0.024d0*z*z
      if (a.lt.amin .or. a.gt.amax) go to 920
      aa = dble(2.5d-3*a)
      zz = dble(1.d-2*z)
      bfis0 = 0.d0
      call lpoly (zz, 7, pz)
      call lpoly (aa, 7, pa)
        do i = 1,7
          do j = 1,7
          bfis0 = bfis0 + elzcof(j,i)*pz(j)*pa(i)
          end do
        end do
      egs = 0.d0
      bfis = bfis0
      amin2 = 1.4d0*z + 0.009d0*z*z
      amax2 = 20.d0 + 3.0d0*z
      if ((a.lt.amin2-5.d0 .or. a.gt.amax2+10.d0) .and. el.gt.0.d0)
     &      go to 930
      el80 = 0.d0
      el20 = 0.d0
      elmax = 0.d0
        do i = 1,4
          do j = 1,5
          el80 = el80 + elmcof(j,i)*pz(j)*pa(i)
          el20 = el20 + emncof(j,i)*pz(j)*pa(i)
          end do
        end do
      call elmaxc (z, a, elmax)
      if (el.lt.0.5d0) return
      x = el20/elmax
      y = el80/elmax
        if (el.le.el20) then
        q = 0.2d0/(el20**2*el80**2*(el20 - el80))
        qa =  q*(4.d0*el80**3 - el20**3)
        qb = -q*(4.d0*el80**2 - el20**2)
        bfis = bfis*(1.d0 + qa*el**2 + qb*el**3)
      else
        aj = (-20.d0*x**5 + 25.d0*x**4 - 4.d0)*(y - 1.d0)**2*y*y
        ak = (-20.d0*y**5 + 25.d0*y**4 - 1.d0)*(x - 1.d0)**2*x*x
        q = 0.2d0/((y - x)*((1.d0 - x)*(1.d0 - y)*x*y)**2)
        qa =  q*(aj*y - ak*x)
        qb = -q*(aj*(2.d0*y + 1.d0) - ak*(2.d0*x + 1.d0))
        z = el/elmax
        a1 = 4.d0*z**5 - 5.d0*z**4 + 1.d0
        a2 = qa*(2.d0*z + 1.d0)
        bfis = bfis*(a1 + (z - 1.d0)*(a2 + qb*z)*z*z*(z - 1.d0))
      endif
      if (bfis.le.0.0d0) bfis = 0.0d0
      if (el.gt.elmax) bfis = 0.0d0
c
c    Now calculate rotating ground-state energy
c
      if (el.gt.elmax .and. el.ne.1000.d0) return
      ell = el/elmax
      if (el.eq.1000.d0) ell = 1.d0
      call lpoly (ell, 9, pl)
        do k = 1,5
          do l = 1,7
            do m = 1,5
            egs = egs + egscof(m,l,k)*pz(l)*pa(k)*pl(2*m-1)
            end do
          end do
        end do
      egs = max(0.0d0, egs)
      return
  899 bfis = 0.0d0
      egs = 0.0d0
      elmax = 0.0d0
      return
  900 write(7,1000)
      go to 899
  910 write(7,1010)
      go to 899
  920 write(7,1020) ia, iz
      go to 899
  930 write(7,1030) ia, iz, el
      go to 899

c ======================================================================

 1000 format(/10x,'* * * * Barfit called with  Z  less than 19 or ',
     & ' greater than 111.'/15x,'bfis is set to 0.0. * * * *')
 1010 format(/10x,'* * * * Barfit called with  Z  greater than 102',
     & ' and  L  not equal to zero.'/15x,'bfis is set to 0.0. * * * *')
 1020 format(/10x,'* * * * Barfit called with  A =',i3,', outside ',
     & /10x,'the allowed values for Z = ',i3,'* * * *')
 1030 format(/10x,'* * * * Barfit called with  A  =',i3,', outside',
     & /10x,' the allowed values for Z = ',i3,'for nonzero  L =',
     & f6.1,' * * * *')

c ======================================================================
      end

      subroutine momfit (iz, ia, el, saimin, saimid, saimx, elmax)

c ======================================================================
c
c    This subroutine returns the three principal-axis moments of inertia,
c    and the angular momentum at which the fission
c    barrier disappears,  Lmax,  in units of h-bar,
c    when called with integer arguments iz, the atomic number,
c    ia, the atomic mass number, and el, the angular momentum in units
c    of h-bar, (Planck's constant divided by 2*pi).
c
c         The moments of inertia at L = 0, 70% of Lmax, 95% of Lmax,
c    and at Lmax are fit to Nz x Na parameter functions, where Nz = 6,
c    Na = 5 for major moments of inertia, and Na = 4 for minor moments
c    of inertia.  For a given value of Z and A, the moments at these 
c    values of L are evaluated from the fitting functions.  Then a f
c    function of L which passes through these values is calculated.
c    The value of the maximum moment of inertia is typically 
c    approximated to within about 1%, while the minor moment is usually
c    somewhat better approximated.  A separate 4x4 function is used to
c    approximate the maximum moment for Z>=80, since the shape at Lmax
c    is axially symmetric and oblate in this region.   The moments of
c    inertia are given in units of 2/5 M-zero (R-zero)**2, the moment
c    of inertia of a rigidly rotating, sharp-surfaced sphere.  
c    R-zero = 1.16*A**(1/3) fermis,
c    and M-zero *c-squared is 931.5016*A - .511004*Z  MeV.
c    If  momfit  is called with (iz,ia) values outside the range 
c    of the fit the moments are set to 0.0, and a message is
c    printed on the default output file.
c
c         The moments of inertia from which the fits were made were
c    calculated in 1983-1985 by A. J. Sierk of Los Alamos National
c    Laboratory,  Group T-9, using  Yukawa-plus-exponential double
c    folded nuclear energy, exact Couloub diffuseness corrections,
c    and diffuse-matter moments of inertia. The parameters of the model
c    are those derived by Moller and Nix in 1979:
c    r-0 = 1.16 fm, as = 21.13 MeV, kappa-s = 2.3  a = 0.68 fm.
c    The diffuseness of the matter and charge distributions used
c    corresponds to a surface diffuseness parameter (defined by Myers)
c    of 0.99 fm.
c
c    Below is a table of test values to check implementation
c    of the program.
c    Z, A,  L    Egnd st  Fiss Bar      Moments of inertia     Lmax
c
c   28, 58, 0    0.00     33.14        0.816 3.603 3.608      46.1
c         ,25   21.36     19.50        0.778 3.662 3.662      46.1
c         ,40   49.66      2.97        0.724 3.648 3.650      46.1
c         ,46.1 59.14      0.00        0.746 3.160 3.160      46.1
c   65,153, 0    0.00     28.88        0.621 3.698 3.698      82.3
c         ,50   19.00     16.16        0.615 3.639 3.639      82.3
c         ,80   45.24      0.26        0.616 2.765 2.788      82.3
c         ,82.3 47.04      0.00        0.682 2.231 2.276      82.3
c   93,229, 0    0.00      3.76        0.715 1.747 1.747      68.1
c         ,45    8.21      1.26        0.765 1.578 1.578      68.1
c         ,68.1 17.96      0.00        1.053 1.053 1.236      68.1
c
c    written by A. J. Sierk,  LANL  T-9
c    Version 1.0   October, 1984
c    Version 2.0   October, 1985.
c    Version 2.1   June, 1986   minor changes made
c
c        Copyright, 1986,  The Regents of the University of California.
c        This software was produced under a U. S. Government contract
c        (W-7405-ENG-36) by the Los Alamos National Laboratory, which is
c        operated by the University of California for the U. S. Department
c        of Energy.  The U. S. Government is licensed to use, reproduce,
c        and distribute this software.  Permission is granted to the public
c        to copy and use this software without charge, provided that this
c        notice and any statement of authorship are reproduced on all
c        copies.  Neither the Government nor the University makes any
c        warranty, expressed or implied, or assumes any liability
c        or responsibility for the use of this software.
c
c ======================================================================
c
      implicit integer (i-n)
      implicit real*8 (a-h, o-z)

      dimension pa(10), pz(10), bi95c(6,4)
c     dimension emxcof(7,5), aizroc(6,5), ai70c(6,5), ai95c(6,5)
      dimension aizroc(6,5), ai70c(6,5), ai95c(6,5)
      dimension ai952c(6,5), aimax2c(6,5), aimaxc(6,5), aimax3c(4,4)
      dimension aimax4c(4,4), bimaxc(6,4), bizroc(6,4), bi70c(6,4)

      data aizroc /
     &  2.34441624d4,-5.88023986d4, 6.37939552d4,-4.79085272d4,
     &  2.27517867d4,-5.35372280d3,
     & -4.19782127d4, 1.09187735d5,-1.24597673d5, 9.93997182d4,
     & -4.95141312d4, 1.19847414d4,
     &  4.18237803d4,-1.05557152d5, 1.16142947d5,-9.00443421d4,
     &  4.48976290d4,-1.10161792d4,
     & -8.27172333d3, 2.49194412d4,-3.39090117d4, 3.33727886d4,
     & -1.98040399d4, 5.37766241d3,
     &  5.79695749d2,-1.61762346d3, 2.14044262d3,-3.55379785d3,
     &  3.25502799d3,-1.15583400d3/
      data ai70c /
     &  3.11420101d4,-7.54335155d4, 7.74456473d4,-4.79993065d4,
     &  2.23439118d4,-4.81961155d3,
     & -7.24025043d4, 1.72276697d5,-1.72027101d5, 1.03891065d5,
     & -4.83180786d4, 1.08040504d4,
     &  7.14932917d4,-1.72792523d5, 1.75814382d5,-1.07245918d5,
     &  4.86163223d4,-1.10623761d4,
     & -2.87206866d4, 6.76667976d4,-6.50167483d4, 3.67161268d4,
     & -1.74755753d4, 4.67495427d3,
     &  1.67914908d4,-3.97304542d4, 3.81446552d4,-2.04628156d4,
     &  7.20091899d3,-1.49978283d3/
      data ai95c /
     & -6.17201449d5, 1.45561724d6,-1.47514522d6, 9.37798508d5,
     & -3.74435017d5, 7.81254880d4,
     &  1.24304280d6,-2.94179116d6, 3.00170753d6,-1.92737183d6,
     &  7.79238772d5,-1.64803784d5,
     & -1.49648799d6, 3.52658199d6,-3.56784327d6, 2.26413602d6,
     & -9.02243251d5, 1.88619658d5,
     &  7.27293223d5,-1.72140677d6, 1.75634889d6,-1.12885888d6,
     &  4.57150814d5,-9.74833991d4,
     & -3.75965723d5, 8.83032946d5,-8.87134867d5, 5.58350462d5,
     & -2.20433857d5, 4.62178756d4/
      data aimaxc /
     & -1.07989556d6, 2.54617598d6,-2.56762409d6, 1.62814115d6,
     & -6.39575059d5, 1.34017942d5,
     &  2.17095357d6,-5.13081589d6, 5.19610055d6,-3.31651644d6,
     &  1.31229476d6,-2.77511450d5,
     & -2.66020302d6, 6.26593165d6,-6.31060776d6, 3.99082969d6,
     & -1.56447660d6, 3.25613262d5,
     &  1.29464191d6,-3.05746938d6, 3.09487138d6,-1.97160118d6,
     &  7.79696064d5,-1.63704652d5,
     & -7.13073644d5, 1.67482279d6,-1.67984330d6, 1.05446783d6,
     & -4.10928559d5, 8.43774143d4/
      data ai952c /
     & -7.37473153d5, 1.73682827d6,-1.75850175d6, 1.11320647d6,
     & -4.41842735d5, 9.02463457d4,
     &  1.49541980d6,-3.53222507d6, 3.59762757d6,-2.29652257d6,
     &  9.21077757d5,-1.90079527d5,
     & -1.80243593d6, 4.24319661d6,-4.29072662d6, 2.71416936d6,
     & -1.07624953d6, 2.20863711d5,
     &  8.86920591d5,-2.09589683d6, 2.13507675d6,-1.36546686d6,
     &  5.48868536d5,-1.14532906d5,
     & -4.62131503d5, 1.08555722d6,-1.09187524d6, 6.87308217d5,
     & -2.70986162d5, 5.61637883d4/
      data aimax2c /
     & -1.16343311d6, 2.74470544d6,-2.77664273d6, 1.76933559d6,
     & -7.02900226d5, 1.49345081d5,
     &  2.36929777d6,-5.60655122d6, 5.7043177d6,-3.66528765d6,
     &  1.47006527d6,-3.15794626d5,
     & -2.82646077d6, 6.66086824d6,-6.72677653d6, 4.27484625d6,
     & -1.69427298d6, 3.58429081d5,
     &  1.39112772d6,-3.29007553d6, 3.34544584d6,-2.14723142d6,
     &  8.61118401d5,-1.84500129d5,
     & -7.21329917d5, 1.69371794d6,-1.69979786d6, 1.07037781d6,
     & -4.20662028d5, 8.80728361d4/
      data aimax3c /
     & -2.88270282d3, 5.30111305d3,-3.07626751d3, 6.56709396d2,
     &  5.84303930d3,-1.07450449d4, 6.24110631d3,-1.33480875d3,
     & -4.20629939d3, 7.74058373d3,-4.50256063d3, 9.65788439d2,
     &  1.23820134d3,-2.28228958d3, 1.33181316d3,-2.87363568d2/
      data aimax4c /
     & -3.34060345d3, 6.26384099d3,-3.77635848d3, 8.57180868d2,
     &  6.76377873d3,-1.26776571d4, 7.64206952d3,-1.73406840d3,
     & -4.74821371d3, 8.89857519d3,-5.36266252d3, 1.21614216d3,
     &  1.46369384d3,-2.74251101d3, 1.65205435d3,-3.74262365d2/
      data bizroc /
     &  5.88982505d2,-1.35630904d3, 1.32932125d3,-7.78518395d2,
     &  2.73122883d2,-3.49600841d1,
     & -9.67701343d2, 2.24594418d3,-2.24303790d3, 1.35440047d3,
     & -4.96538939d2, 6.66791793d1,
     &  1.17090267d3,-2.71181535d3, 2.67008958d3,-1.58801770d3,
     &  5.66896359d2,-8.21530057d1,
     & -3.83031864d2, 9.05191483d2,-9.30560410d2, 5.96618532d2,
     & -2.34403480d2, 3.97909172d1/
      data bi70c /
     &  2.32414810d3,-5.42381778d3, 5.40202710d3,-3.26923144d3,
     &  1.18318943d3,-1.93186467d2,
     & -4.38084778d3, 1.03523570d4,-1.05573803d4, 6.59901160d3,
     & -2.47601209d3, 4.19497260d2,
     &  4.35377377d3,-1.01728647d4, 1.01311246d4,-6.14038462d3,
     &  2.21957562d3,-3.62854365d2,
     & -1.84533539d3, 4.41613298d3,-4.59403284d3, 2.95951225d3,
     & -1.14630148d3, 2.02702459d2/
      data bi95c /
     &  1.55359266d3,-3.58209715d3, 3.50693744d3,-2.03992913d3,
     &  7.05498010d2,-1.49075519d2,
     & -2.86876240d3, 6.77107086d3,-6.90300614d3, 4.20246063d3,
     & -1.50290693d3, 3.13662258d2,
     &  2.60138185d3,-5.95414919d3, 5.70261588d3,-3.17188958d3,
     &  9.89207911d2,-1.76320647d2,
     & -1.75198402d3, 4.16635208d3,-4.25212424d3, 2.59953301d3,
     & -9.09813362d2, 1.51070448d2/
      data bimaxc /
     &  4.17708254d3,-8.59358778d3, 6.46392215d3,-8.84972189d2,
     & -1.59735594d3, 1.39662071d3,
     & -1.56318394d4, 3.54574417d4,-3.35945173d4, 1.65495998d4,
     & -3.32021998d3,-1.46150905d3,
     &  1.41292811d4,-3.11818487d4, 2.77454429d4,-1.19628827d4,
     &  1.28008968d3, 1.66111636d3,
     & -1.92878152d4, 4.56505796d4,-4.66413277d4, 2.89229633d4,
     & -1.07284346d4, 1.50513815d3/

      data pi /3.141592653589793d0/

c ======================================================================
c
c    The program starts here
c
      if (iz.lt.19 .or. iz.gt.101) go to 900
      z = dble(iz)
      a = dble(ia)
      amin = 1.2d0*z + 0.01d0*z*z
      amax = 5.8d0*z - 0.024d0*z*z
      if (a.lt.amin .or. a.gt.amax) go to 910
      aa = a/400.
      zz = z/100.
      amin2 = 1.4d0*z + 0.009d0*z*z
      amax2 = 20.d0 + 3.0d0*z
      if ((a.lt.amin2-5.d0 .or. a.gt.amax2+10.d0) .and. el.gt.0.d0) 
     &     go to 920
      call lpoly (aa, 6, pa)
      call lpoly (zz, 7, pz)
      call elmaxc (z, a, elmax)
      ell = el/elmax
      if (el.eq.1000.d0) ell = 1.0d0
      aizro = 0.d0
      ai70 = 0.d0
      aimax = 0.d0
      ai95 = 0.d0
      aimin = 0.0
      bizro = 0.d0
      bi70 = 0.d0
      bimax = 0.d0
      bi95 = 0.d0
      aimax2 = 0.d0
      ai952 = 0.d0
c
c    Now calculate rotating moments of inertia
c
      if (el.gt.elmax .and. el.lt.1000.d0) return
        do l=1,6
          do k = 1,5
          aizro = aizro + aizroc(l,k)*pz(l)*pa(k)
          ai70 = ai70 + ai70c(l,k)*pz(l)*pa(k)
          ai95 = ai95 + ai95c(l,k)*pz(l)*pa(k)
          aimax = aimax + aimaxc(l,k)*pz(l)*pa(k)
          ai952 = ai952 + ai952c(l,k)*pz(l)*pa(k)
          aimax2 = aimax2 + aimax2c(l,k)*pz(l)*pa(k)
          end do
          do k = 1,4
          bizro = bizro + bizroc(l,k)*pz(l)*pa(k)
          bi70 = bi70 + bi70c(l,k)*pz(l)*pa(k)
          bi95 = bi95 + bi95c(l,k)*pz(l)*pa(k)
          bimax = bimax + bimaxc(l,k)*pz(l)*pa(k)
          end do
        end do
      ff1 = 1.0d0
      ff2 = 0.0d0
      fg1 = 1.0d0
      fg2 = 0.0d0
      if (iz.gt.70) then
c       aimaxh = 0.0
c       aimidh = 0.0
        aimaxh = 0.d0
        aimidh = 0.d0
          do l = 1,4
            do k = 1,4
            aimaxh = aimaxh + aimax3c(l,k)*pz(l)*pa(k)
            aimidh = aimidh + aimax4c(l,k)*pz(l)*pa(k)
            end do
          end do
        if (iz.gt.80) ff1 = 0.0d0
        if (iz.ge.80) fg1 = 0.0d0
        if (bimax.gt.0.95d0) fg1 = 0.0d0
        if (aimaxh.gt.aimax) ff1 = 0.0d0
        ff2 = 1.0d0 - ff1
        fg2 = 1.0d0 - fg1
        aimax = aimax*ff1 + ff2*aimaxh
        aimax2 = aimax2*ff1 + ff2*aimidh
      endif
      bimax = bimax*fg1 + aimidh*fg2
      saizro = max(aizro, 0.d0)
      sai70 = max(ai70, 0.d0)
      sai95 = max(ai95, 0.d0)
      saimax = max(aimax, 0.d0)
      sai952 = max(ai952, 0.d0)
      simax2 = max(aimax2, 0.d0)
      sbimax = max(bimax, 0.d0)
      sbi70 = max(bi70, 0.d0)
      sbi95 = max(bi95, 0.d0)
      sbizro = max(bizro, 0.d0)
      q1 = -3.148849569d0
      q2 =  4.465058752d0
      q3 = -1.316209183d0
      q4 =  2.26129233d0
      q5 = -4.94743352d0
      q6 =  2.68614119d0
      gam = - 20.d0*log(dabs(saizro-sai95)/dabs(saizro-saimax))
      aa = q1*saizro + q2*sai70 + q3*sai95
      bb = q4*saizro + q5*sai70 + q6*sai95
      gam2 = - 20.d0*log(dabs(saizro-sai952)/dabs(saizro-simax2))
      aa2 = q1*saizro + q2*sai70 + q3*sai952
      bb2 = q4*saizro + q5*sai70 + q6*sai952
      aa3 = q1*sbizro + q2*sbi70 + q3*sbi95
      bb3 = q4*sbizro + q5*sbi70 + q6*sbi95
      gam3 = 60.d0
      alpha = pi*(ell - 0.7d0)
      beta = 5.d0*pi*(ell - 0.9d0)
      silt = saizro + aa*ell**2 + bb*ell**4
      sjlt = sbizro + aa3*ell**2 + bb3*ell**4
      silt2 = saizro + aa2*ell**2 + bb2*ell**4
      sigt = saizro + (saimax-saizro)*exp(gam*(ell-1.0d0))
      sjgt = sbi95 + (sbimax-sbi95)*exp(gam3*(ell-1.0d0))
      sigt2 = saizro + (simax2-saizro)*exp(gam2*(ell-1.0d0))
      f1 = silt*(cos(alpha))**2 + sigt*(sin(alpha))**2
      f2 = silt*(cos(beta))**2 + sigt*(sin(beta))**2
      f1m = silt2*(cos(alpha))**2 + sigt2*(sin(alpha))**2
      f2m = silt2*(cos(beta))**2 + sigt2*(sin(beta))**2
      f3 = sjlt*(cos(alpha))**2 + sjgt*(sin(alpha))**2
      f4 = sjlt*(cos(beta))**2 + sjgt*(sin(beta))**2
      if (ell.le.0.7d0) then
c
c    ell is less than 0.7,  use ilt
c
        saimin = sjlt
        saimx = silt
        saimid = silt2
        if (ff2.gt.0.01d0 .and. fg2.gt.0.01d0) go to 90
      elseif (ell.le.0.95d0) then
c
c    ell is greater than 0.7, less than 0.95  use first l. c.
c
        saimx = f1
        saimin = f3
        saimid = f1m
        if (ff2.gt.0.01d0 .and. fg2.gt.0.01d0) go to 90
      else
c
c    ell is greater than 0.95,  use 2nd l. c.
c
        saimx = f2
        saimin = f4
        saimid = f2m
        if (ff2.gt.0.01d0 .and. fg2.gt.0.01d0) go to 90
      endif
      go to 100
c
c    For nuclei with z gt 80 use 4th order function with separate fit to
c    imax and imax2
c
   90 q1 = 4.001600640d0
      q2 = 0.960784314d0
      q3 = 2.040816327d0
      aa3 =  q1*sai70 - q2*saimax - (1.d0+q3)*saizro
      bb3 = -q1*sai70 + (1.d0+q2)*saimax + q3*saizro
      aa4 =  q1*sai70 - q2*simax2 - (1.d0+q3)*saizro
      bb4 = -q1*sai70 + (1.d0+q2)*simax2 + q3*saizro
      saimx = saizro + aa3*ell**2 + bb3*ell**4
      saimid = saizro + aa4*ell**2 + bb4*ell**4
  100 if (saimid.gt.saimx) saimid = saimx
      saimin = amax1(saimin, 0.d0)
      return

  900 write(7,1000)
  905 aimax = 0.0d0
      saimx = 0.0d0
      saimin = 0.0d0
      aimin = 0.0d0
      saimid = 0.0d0
      aimid = 0.0d0
      elmax = 0.0d0
      return
  910 write(7,1010) ia, iz
      go to 905
  920 write(7,1020) ia, iz, el
      go to 905

c ======================================================================

 1000 format(/10x,'* * * * Momfit called with  Z  less than 19 or ',
     & ' greater than 101.'/18x,'moments are set to 0.0. * * * *')
 1010 format(/10x,'* * * * Momfit called with  A =',i3,', outside ',
     & 'the allowed values for Z = ',i3,' * * * *')
 1020 format(/10x,'* * * * Momfit called with  A  =',i3,', outside',
     & /10x,'the allowed values for Z = ',i3,'for nonzero  L =',f6.1,
     & ' * * * *')

c ======================================================================
      end

      subroutine lpoly (x, n, pl)

c ======================================================================
c
c    This subroutine calculates the ordinary Legendre Polynomials of
c    order 0 to n-1 of argument  x  and stores them in the vector
c    pl.  They are calculated by recursion relation from the first two
c    polynomials.
c
c    written by A. J. Sierk   LANL  T-9  February,1984
c
c    NOTE:  pl and x must be real*8 !
c
      implicit integer (i-n)
      real*8 pl, x

      dimension pl(10)

c ======================================================================

      pl(1) = 1.d0
      pl(2) = x
      do i = 3,n
        pl(i) = (dble(2*i-3)*x*pl(i-1) - dble(i-2)*pl(i-2))/dble(i-1)
        end do
      return

c ======================================================================
      end

      subroutine elmaxc (zn, an, elmax)

c ======================================================================

      real*8 an, elmax, zn

c ======================================================================
c
c   This subroutine evaluates a 2-dimensional fitting function
c   to arrive at an approximation to Lmax(Z,A), the value of
c   angular momentum in units of h-bar at which the fission
c   saddle point and the fission barrier vanishes.  Lmax is the
c   maximum amount of angular momentum which the nucleus can
c   sustain without fissioning.
c
c   The 2-dimensional fit was made to 192 calculated values of
c   Lmax for values of Z from 10 to 110, with A values spanning
c   a range larger than any nuclei which might be formed in
c   a reaction.
c
c   The 192 calculated values of Lmax were calculated in the
c   three-quadratic-surface axially symmetric parametrization,
c   or in the triaxial Legendre Polynomial parametrization
c   using Yukawa-plus-exponential nuclear energy, diffuse
c   Coulomb energy, and diffuse-matter moments of inertia by
c   A. J. Sierk of Los Alamos National Laboratory, Group T-2.
c
c   This subroutine written February, 1994.
c   Updated constants and extended Z and A range of fit, May, 1996.
c
c ======================================================================

c mb remove duplicate elmax below 1/28/15
      real*8 b(40), pla(5), plz(8), xa, xz
c     real*8 b(40), elmax, pla(5), plz(8), xa, xz

      integer ia, ib, iz, npa, npz

      data npa, npz /5,8/

      data b
     & /-3.30787006E+05, 8.34310803E+05,-4.11955431E+05, 9.82177374E+05,
     &   3.40653902E+05, 8.39479538E+05,-2.13869250E+06, 1.11540268E+06,
     &  -2.52645986E+06,-7.72016242E+05,-1.00538896E+06, 2.59234090E+06,
     &  -1.49269622E+06, 3.06066515E+06, 7.02514732E+05, 7.97067135E+05,
     &  -2.08607738E+06, 1.33444088E+06,-2.52998737E+06,-3.45973915E+05,
     &  -4.15169782E+05, 1.10638461E+06,-7.71715934E+05, 1.45998634E+06,
     &   6.29116062E+04, 1.07407080E+05,-3.04281993E+05, 2.18876998E+05,
     &  -5.37873263E+05, 2.14561453E+04, 1.64345870E+03, 6.52424621E+03,
     &   3.41057615E+03, 1.03029197E+05,-1.17188782E+04,-5.44933809E+03,
     &   1.23028722E+04,-1.25511319E+04,-5.30918539E+03, 1.17305352E+03/

c ======================================================================

      xz = zn/100.d0
      xa = an/320.d0
      elmax = 0.d0
      call lpoly2 (xz, plz, npz)
      call lpoly2 (xa, pla, npa)
        do iz = 1,npz
          do ia = 1,npa
          ib = ia + npa*(iz-1)
          elmax = elmax + b(ib)*plz(iz)*pla(ia)
          end do
        end do

      return

c ======================================================================
c
      end

      subroutine lpoly2 (x, ple, n)

c ======================================================================

      real*8 x, ple(n)

      integer i, n

c ======================================================================

      ple(1) = 1.d0
      ple(2) = x
      ple(3) = 0.5d0*(3.d0*x*x - 1.d0)
      ple(4) = 0.5d0*x*(5.d0*x*x - 3.d0)
      if (n.ge.5) then
        do i = 5,n
        ple(i) = (dble(2*i-3)*x*ple(i-1) - dble(i-2)*ple(i-2))/dble(i-1)
        end do
      endif
      return

c ======================================================================
      end
c ======================================================================
      subroutine ldload(ldopt)
      common/atsv/atsave
      common/lab11/powaz(15,24,2200)
      REAL*8 pow,gam,powaz
      common/ss/sor,rr
      common/lab10/pow(4,9999),gam(9999)
cneadb activate the following statemet if running on 32-bit word machines
cneadb
      common/ug/iz ,ia
      common/sf/m3,kplt
      common/sft5/exc(18,27),xmxx
      common/lym1/be(18,27,8),pair(18,27),xmas(18,27),delshl(18,27),amas
     1(18,27)
c     common/bym1/symb(18,27),symbp(18,27),symbs(18,27)
      common/pl3/lz,nepr,eb(100),rcsp(100),zee,amass,plex
      common/shft/mc,mp,inver,ike,ipch,kq,qval,ap,at,zp,zt,cld,barfac
      common/lab12/pof(4,9999)
      REAL*8 pof
 
      go to 600
      do 30 jk=1,4
      do 30 ie=1,9999
      pow(jk,ie)=0.
      pof(jk,ie)=0.
   30 continue
       m5=m3
       if(m5.gt.4)m5=4
c     do 10 jk=1,m5
      if(jk.eq.1) then
      n=0
      m=1
      else if(jk.eq.2) then
      n=1
      m=0
      else if(jk.eq.3) then
      n=2
      m=2
      else if(jk.eq.4) then
      n=1
      m=1
      end if
      itz=iz+n
      ita=ia+m
 600  continue     
         jk=1
c jk is now a dummy if using powaz array
      do 10 iz=1,15
       do 9 ia=1,24
c     excit=exc(iz,ia)
      a=atsave+ap+ 2. -float(ia)-float(iz)
      iia=a+0.001
      iiz=zee+1.001-float(iz)
      zz=iiz
      iin=iia-iiz
       itest=iin*iiz
        if(itest.lt.4)go to 9
        if(amass.lt.4..or.iiz.lt.3.or.iin.lt.1)go to 9
      excit=exc(iz,ia)
      epair=pair(iz,ia)/10.
      if(ldopt.eq.3)call gc(a,zz,excit,jk)
      if(ldopt.eq.4)call mmbc(iiz,iia,jk)
    9 continue
   10 continue
      return
      end
 
c*****************************************************************************
      subroutine mmbc(iiz,iia,jk)
c*****************************************************************************
c Calculates TOTAL [constant-temp+Fermi-gas] level density for light nuclei
c (from Gilbert-Cameron eqns.) M.B. Chadwick, Oct. 29th. 1993.
c See Refs. Young, Arthur and Chadwick, LANL Document LA-12343 (1993)
c and Gilbert and Cameron, Can. J. Phys. 43, 1446 (1965)
c Constant-temperature and Fermi-gas regions match continuously in both
c the level density and it's first derivative.
c*****************************************************************************
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c INPUT & OUTPUT quantities:-
c iiz=nucleus Z, iia=nucleus A.
c calcs. array rhotot= the TOTAL (not 'observable") level density
c The following parameters are defined in the DATA blocks:-
c  id     =nuclide indicator (1000*Z+A)
c  a      ="a" parameter (in Fermi-gas expression)
c  pair   =pairing parameter (in Fermi-gas expression)
c  e0gc   = E-0 parameter (in constant-T expression)
c  tgc    =temparature param. (in constant-T expression)
c  ematgc =matching energy between constant-T and Fermi-gas regions.
c also, sigcon =0.146= spin cut-off expression (needed for constant-T expr.)
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c **NUCEI NOT INCLUDED**:- anything with N>9,Z>9, and in addition
c nuclei with id=1000Z+A=1,2,3,4,5,6,7,8,9,1001,1002,1003,2002,2003,
c 3003,9018 not included (none should be important in n,p+16O,14N,12C
c reactions)
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c DATA statements built from file LEVDEN, with mkdata.f, see
c Chadwick's LLNL pd4/codes/gcplt/levdens-codes/alice-gc area.
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      dimension a1(0:10,0:10),pair1(0:10,0:10),tgc1(0:10,0:10),
     & e0gc1(0:10,0:10),ematgc1(0:10,0:10)
      common/lab10/pow(4,9999),gam(9999)
      common/lab11/powaz(15,24,2200)
      REAL*8 pow,gam,powaz
c arrays dimensioned for N and Z from 0 to 9.
      iin=iia-iiz
      pi=3.1415927
      amass=float(iia)
 
c DATA blocks of level density parameters follow. Chadwick's Oct.93------------
c work, matching P.G. Young's extension of Cook systematics to experimental
c low-lying cumulative number of levels, within a Gilbert-Cameron approach.
      data ((a1(jjn,jjz),jjz=0,9),jjn=0,9)
     &/0.000,   0.000,   0.000,   0.000,   0.555,   0.685,   0.811,
     &0.934,   1.054,   1.170,   0.000,   0.000,   0.000,   0.570,
     &0.704,   0.834,   0.961,   1.085,   1.205,   1.321,   0.000,
     &0.000,   0.600,   0.741,   0.879,   1.013,   1.144,   1.271,
     &1.395,   1.516,   0.000,   0.629,   0.777,   0.922,   1.065,
     &1.203,   1.337,   1.469,   1.596,   1.721,   0.000,   0.815,
     &0.967,   1.116,   1.262,   1.404,   1.543,   1.678,   1.810,
     &1.938,   0.000,   1.012,   1.168,   1.321,   1.471,   1.817,
     &1.760,   1.899,   2.034,   2.166,   0.000,   1.220,   1.380,
     &1.537,   1.691,   1.841,   1.988,   2.131,   2.270,   2.406,
     &0.000,   1.440,   1.604,   1.765,   1.923,   2.077,   2.227,
     &1.988,   2.130,   2.658,   0.000,   1.671,   1.839,   2.004,
     &2.166,   2.324,   2.478,   2.100,   2.777,   2.921,   0.000,
     &1.912,   2.085,   2.253,   2.419,   2.581,   2.739,   2.894,
     &3.045,   0.000/
 
      data ((pair1(jjn,jjz),jjz=0,9),jjn=0,9)
     &/0.000,   0.000,   0.000,   0.000,   4.500,   0.000,   3.950,
     &0.000,   3.400,   0.000,   0.000,   0.000,   0.000,   0.000,
     &4.500,   0.000,   3.950,   0.000,   3.400,   0.000,   0.000,
     &0.000,  10.300,   5.250,   9.750,   5.250,   9.200,   5.250,
     &8.650,   5.250,   0.000,   0.000,   5.050,   0.000,   4.500,
     &0.000,   3.950,   0.000,   3.400,   0.000,   0.000,   4.700,
     &9.750,   4.700,   9.200,   4.700,   8.650,   4.700,   8.100,
     &4.700,   0.000,   0.000,   5.050,   0.000,   4.500,   0.000,
     &3.950,   0.000,   3.400,   0.000,   0.000,   4.150,   9.200,
     &4.150,   8.650,   4.150,   8.100,   4.150,   7.550,   4.150,
     &0.000,   0.000,   5.050,   0.000,   4.500,   0.000,   3.950,
     &0.000,   3.400,   0.000,   0.000,   3.580,   8.630,   3.580,
     &8.080,   3.580,   7.530,   3.580,   6.980,   3.580,   0.000,
     &0.000,   5.050,   0.000,   4.500,   0.000,   3.950,   0.000,
     &3.400,   0.000/
 
      data ((tgc1(jjn,jjz),jjz=0,9),jjn=0,9)
     & /0.000,   0.000,   0.000,   0.000, 226.085, 149.446,   7.408,
     & 81.365,   6.055,  52.549,   0.000,   0.000,   0.000,  14.829,
     &141.530, 101.472,   6.609,  60.852,   5.561,  41.547,   0.000,
     &  0.000,  19.192,   8.265,  10.482,   6.877,   8.059,   6.124,
     &  7.250,   5.645,   0.000,  12.737,   9.466,   9.384,   6.229,
     &  6.128,  13.200,  33.875,   4.681,  25.007,   0.000,   7.489,
     & 44.685,   7.152,   5.817,   6.996,   5.789,   5.081,   6.294,
     &  4.774,   0.000,  69.698,   6.169,   4.834,   5.653,   3.962,
     &  4.687,   4.059,  18.184,  16.137,   0.000,   5.626,   7.294,
     &  4.737,   5.544,   4.051,   3.606,   4.021,   4.490,   3.543,
     &  0.000,  35.195,   5.224,  23.822,   5.000,   3.513,   3.974,
     &  3.320,   3.180,   2.322,   0.000,   4.486,   6.301,   4.104,
     &  4.552,   4.194,   3.371,   3.698,   2.287,   2.941,   0.000,
     & 20.441,   4.616,  14.986,   4.131,   3.341,   3.097,   2.748,
     &  2.537,   0.000/
 
      data ((e0gc1(jjn,jjz),jjz=0,9),jjn=0,9)
     &  /0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,
     &   0.000,   0.000,   0.000,   0.000,   0.000,   0.000, -32.436,
     &   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,   0.000,
     &   0.000, -21.677,   0.000, -17.574,   0.000,   0.000,   0.000,
     &   0.000,   0.000,   0.000, -23.941,  -4.083, -20.640,  -0.546,
     &  -7.481, -21.296,   0.000,   0.000,   0.000,   0.000,   0.000,
     &-125.600,  -7.271,   2.207, -15.786,  -4.259,   0.000,   0.000,
     &   0.000,   0.000,   0.000,   0.000,  -2.720,  -5.801,  -4.103,
     &  -4.060,  -5.307, -37.625,   0.000,   0.000,   0.000,   0.000,
     &  -1.176,  -5.526,  -0.598,   4.649,   0.627,  -6.213,  -2.233,
     &   0.000,   0.000,   0.000,   0.000,  -9.024,  -3.647,   0.284,
     &  -1.863,   1.369,  -0.614,   0.000,   0.000,   0.000,   0.000,
     &  -4.951,  -7.948,   1.795,  -1.125,   6.244,  -1.724,   0.000,
     &   0.000,   0.000,   0.000,   0.000,  -6.324,  -1.534,  -3.789,
     &   0.345,   0.000/
 
      data ((ematgc1(jjn,jjz),jjz=0,9),jjn=0,9)
     &  /0.000,   0.000,   0.000,   0.000,   8.656,   3.386,  14.199,
     &   2.508,   8.917,   2.023,   0.000,   0.000,   0.000,  74.147,
     &   7.796,   2.797,  10.082,   2.174,   7.670,   1.803,   0.000,
     &   0.000,  15.545,  14.580,  70.836,   9.963,  12.308,   8.486,
     &  11.055,   7.733,   0.000,  57.461,  10.714,  49.038,  22.189,
     &  23.130,   5.998,   1.632,   6.104,   1.408,   0.000,  13.381,
     &  12.254,  36.741,  30.957,  49.991,  40.640,   6.949,   9.846,
     &   6.531,   0.000,   2.324,   8.932,   4.393,  31.931,  14.132,
     &  26.351,  16.913,   4.606,   1.139,   0.000,   8.170,  11.637,
     &  21.527,  41.920,  19.862,  20.840,   5.914,  38.374,  22.166,
     &   0.000,   1.663,   7.438,   1.375,  35.777,  12.951,   5.595,
     &   9.283,  13.037,   2.488,   0.000,   6.154,  10.335,   5.524,
     &  37.743,  30.399,  24.025,  19.236,   9.105,  18.730,   0.000,
     &   1.277,   6.708,   1.099,   5.896,  17.336,  19.550,  12.223,
     &  14.027,   0.000/
c END OF DATA----------------------------------------------------
 
c define level density parameters:-
      id=1000*iiz+(iin+iiz)
      if(id.le.1003)return
      if(id.eq.2002.or.id.eq.2003)return
      if(id.eq.3003)return
      if(id.ge.9018)return
c define id and return if illegal entry
        iia=iin+iiz
        a=a1(iin,iiz)
        pair=pair1(iin,iiz)
        tgc=tgc1(iin,iiz)
        e0gc=e0gc1(iin,iiz)
        ematgc=ematgc1(iin,iiz)
        sigcon=0.146
c this is the Facchini constant in the spin cut-off.
c     id=1000*iiz+(iin+iiz)
c     write(17,8)id,a,pair,tgc,e0gc,ematgc,xnlgcl,ecgcl,
c    & xnlgc,ecgc,sigcon
c8    format(i4,5f8.3,f4.0,f7.3,f4.0,f7.3,f6.3)
 
c If nucleus not included in this subroutine, STOP:-
      if(iiz.ge.10.or.iin.ge.10.or.a.lt.0.01)return
c     if(iiz.ge.10.or.iin.ge.10.or.a.lt.0.01)go to 999
 
c Calculate total level density and put in array rhotot:-
      do 200 ne=1,2200
         excit=0.1*float(ne-1) + 0.05
c        energy grid 0.05 MeV to 299.95 MeV step 0.1 MeV
         ures=excit-pair
         if(ures.le.0)ures=0.
         sig=sqrt(sigcon*(amass**(2./3.))*sqrt(a*ures))
c if excit is in the constant-temp region, using this spin cut-off can
c lead to problems - it is possible to have a zero total, but non-zero
c observable, level density. To avoid this, I assume an e-dependence
c of excit**0.25 (not U**0.25)
c of the spin cut-off in the constant-T region.
         sigmat=sqrt( sigcon*(amass**(2./3.))*sqrt(a*(ematgc-pair)) )
         if(excit.le.ematgc)sig=sigmat*((excit/ematgc)**0.25)
 
         if(excit.le.ematgc)call contemp(excit,e0gc,tgc,rho)
         if(excit.gt.ematgc)call fgas(excit,a,pair,sig,rho)
c Now we convert the "observable" level density to the "total" level density:
c        pow(jk,ne)=rho*sqrt(2.*pi)*sig
         powaz(iiz,iia,ne)=rho*sqrt(2.*pi)*sig
 200  continue
 
      return
      continue
      write(17,*)'Level density for nucleus',id,' not defined in '
      write(17,*)'Light-nucleus Gilbert-Cameron Subroutine RHOGC'
c      stop
      return
      end
 
      subroutine fgas(excit,a,pair,sig,rho)
c determines Fermi-gas observable level density (Gilbert-Cameron)
      u=excit-pair
      if(u.le.0.)then
         rho=0.
         go to 70
      endif
      con=12.*sqrt(2.)
      rho=exp(2.*sqrt(a*u))
      rho=rho/(con*sig*u*((a*u)**0.25))
 70   continue
      return
      end
 
      subroutine contemp(excit,e0gc,tgc,rho)
c determines constant-temperature observable level density (Gilbert-Cameron)
      rho=(1./tgc)*exp((excit-e0gc)/tgc)
      return
      end
 
      subroutine ldcalc4
c this is for gamma emission densities
c 5/8/03 make this do a KR level density with no shell effect for o-o nuclei
c then make gamma level densities;probable error using ald/2.
      common/atsv/atsave
      common/levop/ldopt,igam
      common/qq/gow(1100)
      REAL*8 gow
      common/shft/mc,mp,inver,ike,ipch,kq,qval,ap,at,zp,zt,cld,barfac
      common/lab11/powaz(15,24,2200)
c
      common/krld/powkr(9999),gamkr(9999),powgkr(9999)
      REAL*8 powkr,gamkr,powgkr
      common/pl5/na,nz,title(20)
      common/lab10/pow(4,9999),gam(9999)
cneadb activate the following statement if running on 32-bit word machine
      REAL*8 pow,gam,powaz
       REAL*8 dens,cnorm,arg,ssor
cenadb
      common/ss/sor,rr
      common/pl3/lz,nepr,eb(100),rcsp(100),zee,amass,plex
      common/sft5/exc(18,27),xmxx
      common/lym1/be(18,27,8),pair(18,27),xmas(18,27),delshl(18,27),amas
     1(18,27)
c     common/bym1/symb(18,27),symbp(18,27),symbs(18,27)
 
c we need to provide path to epair in loop
      integer hilim
      pi=3.14159
      pisq=9.8696
      srtpi=1.7725/12.
      tpi=6.283
      con=1./(12.*sqrt(2.))
      gf=6./pisq
c     if(na.gt.22)na=22
c     if(nz.gt.13)nz=13
       iz=1
       ia=1
c     do 200 iz=1,nz+2
c     do 200 ia=1,na+2
      a=atsave+ap+2.-iz-ia
      a=a-1.
c substitute acn-1 for loop over masses 5/8/03
      a3=a**0.33333
      a23=a3*a3
c     bofn=be(iz,ia,1)
      bofn=be(1,1,1)
      bofp=be(1,1,2)
c     bofp=be(iz,ia,2)
c     bofn=0.
c     bofp=0.
c     if(bofn.eq.0..or.bofp.eq.0.) then
c     aldm=a/9.
c     else
          if(bofn.lt.1.)bofn=1.
          if(bofp.lt.1.)bofp=1.
      aldm=0.11035*a-0.25132*a23+1.4*a23/bofn + 0.4*a23/bofp
c     end if
      aldm=aldm/2.
c     epair=pair(iz,ia)/10.
c     epair=-2.*11./sqrt(a)
      epair=0.
      w=0.185*a3
      piw=pi*w
      em2bar=0.240*a23
      conmul=pi/(12.*sqrt(12.*em2bar))
      umin=2.5+150./a
      delsh=0.
c put in zero shell correction
c     ud=umin-(epair+delshl(iz,ia))-0.1
c     ud=umin-(epair+delsh)-0.1
      ud=umin
      if(ud.lt.0.1)umin=-(epair+delsh)+0.1
      if(umin.gt.2.)umin=2.
      if(umin.lt.1.0)umin=1.0
      ilow=10.*umin
cneadbif(excit.gt.50.)excit=50.    letter from the author 24-March-1992
      hilim=1100
c jan 02 set for 50 mev max
      if(hilim.lt.1)go to 900
      u=110.-epair+0.050
      dels=delshl(iz,ia)
c     dels=delsh
      if(hilim.lt.ilow)hilim=ilow
cneadb                            letter from the author 24-March-1992
      if((hilim+ilow-1).gt.19999)hilim=hilim-ilow
cneadb
      do 850 iie=ilow,hilim
      ie=hilim+ilow-iie
c     u=u-0.100
      u=float(ie)*0.1-0.05
      if(u.lt.0.05)go to 850
      tp=sqrt(u/aldm)
   10 t=tp
      t1=t
      t2=t+0.005
      u1=aldm*t1*t1-dels
      u2=aldm*t2*t2 -dels
      xx1=piw*t1
      if(xx1.gt.15.)go to 700
      xx2=piw*t2
      u1=u1+dels*(xx1*xx1*cosh(xx1)/(sinh(xx1))**2)
      u2=u2+dels*(xx2*xx2*cosh(xx2)/(sinh(xx2))**2)
  700 continue
      tp=t1+(u-u1)*((t2-t1)/(u2-u1))
      if(abs(t-tp).gt.0.005)go to 10
      t=tp
      s=2.*aldm*t
      xxx=piw*t
      if(xxx.gt.15.)go to 99
      s=s+(dels/t)*(xxx*xxx*cosh(xxx)/((sinh(xxx))**2)-xxx/sinh(xxx))
   99 aeff=s/2./t
      ssor=s-sor
      arg=dexp(ssor)
      dens=arg*conmul/((aeff**0.75)*sqrt(t)*(u**1.25+t))
c     dens=dexp(s-sor)*conmul/((aeff**0.75)*sqrt(t)*(u**1.25+t))
      gow(ie)=dens
c     powaz(iz,ia,ie)=dens
  850 continue
      cnorm=dens/exp(u/t)
c  now add constant temperature part of level density
      do 860 iie=1,ilow
      ie=ilow+1-iie
      u=float(ie)*0.1-0.05
      gow(ie)=cnorm*exp(u/t)
      if(gow(ie).le.0.)write(*,*)'iie,ie,gow(ie),cnorm,conmul',
     1iie,ie,gow(ie),cnorm,conmul
c     powaz(iz,ia,ie)=cnorm*exp(u/t)
c     u=u-0.1
  860 continue
  900 continue
      continue
      continue
      continue
      continue
 
      continue
      return
      end
 
      subroutine ldcalc3
c 5/8/03 make this do a KR level density with no shell effect for o-o nuclei
c then make fission width from this
      common/ldf/gamtem(9999)
      REAL*8 gamtem
      common/atsv/atsave
      common/levop/ldopt,igam
      common/qq/gow(1100)
      REAL*8 gow
      common/shft/mc,mp,inver,ike,ipch,kq,qval,ap,at,zp,zt,cld,barfac
      common/lab11/powaz(15,24,2200)
c
      common/krld/powkr(9999),gamkr(9999),powgkr(9999)
      REAL*8 powkr,gamkr,powgkr
      common/pl5/na,nz,title(20)
      common/lab10/pow(4,9999),gam(9999)
cneadb activate the following statement if running on 32-bit word machine
      REAL*8 pow,gam,powaz
       REAL*8 dens,cnorm,arg,ssor
cenadb
      common/ss/sor,rr
      common/pl3/lz,nepr,eb(100),rcsp(100),zee,amasq,plex
      common/sft5/exc(18,27),xmxx
      common/lym1/be(18,27,8),pair(18,27),xmas(18,27),delshl(18,27),amas
     1(18,27)
c     common/bym1/symb(18,27),symbp(18,27),symbs(18,27)
 
c we need to provide path to epair in loop
      integer hilim
c tepm cld wire
c     cl=1.06
      pi=3.14159
      pisq=9.8696
      srtpi=1.7725/12.
      tpi=6.283
      con=1./(12.*sqrt(2.))
      gf=6./pisq
       iz=1
       ia=1
      a=atsave+ap+2.-iz-ia
      a=a-1.
c substitute acn-1 for loop over masses 5/8/03
      a3=a**0.33333
      a23=a3*a3
      bofn=0.
      bofp=0.
      aldm=(a/9.)
      bldm=(a/9.)
      epair=0.
      w=0.185*a3
      piw=pi*w
      em2bar=0.240*a23
      conmul=pi/(12.*sqrt(12.*em2bar))
      umin=2.5+150./a
      delsh=0.
      ud=umin
      if(umin.lt.1.0)umin=1.0
      ilow=10.*umin
cneadbif(excit.gt.50.)excit=50.    letter from the author 24-March-1992
c     if(excit.gt.60.)excit=60.
c     hilim=10.*excit
      hilim=600
c jan 02 set for 50 mev max
c     if(hilim.lt.1)go to 900
c     u=excit-epair+0.050
c     u=60.-epair+0.050
c     dels=delshl(iz,ia)
      dels=delsh
c     if(hilim.lt.ilow)hilim=ilow
cneadb                            letter from the author 24-March-1992
c     if((hilim+ilow-1).gt.999)hilim=hilim-ilow
cneadb
        rf=1.02
      do 850 iie=ilow,hilim
      ib=iie
      ie=iie
c here add multipliers for aldm as in fermld
c       rf=1.02
c     if(ib.gt.70)rf=1.04
c     if(ib.gt.140)rf=1.03
c     if(ib.gt.210)rf=1.02
c     if(ib.gt.280)rf=1.05
c     if(ib.gt.360)rf=1.04
c     if(ib.gt.420)rf=1.03
c     if(ib.gt.480)rf=1.02
c     if(ib.gt.540)rf=1.00
c     u=u-0.100
      u=float(ib)*0.1-0.05
      aldm=bldm*rf
c above puts fudge factors into the fission level density parameter
      if(u.lt.0.1)go to 850
      tp=sqrt(u/aldm)
       if(ib.eq.ilow)tsave=tp
   10 t=tp
      t1=t
      t2=t+0.005
      u1=aldm*t1*t1-dels
      u2=aldm*t2*t2 -dels
      xx1=piw*t1
      if(xx1.gt.15.)go to 700
      xx2=piw*t2
      u1=u1+dels*(xx1*xx1*cosh(xx1)/(sinh(xx1))**2)
      u2=u2+dels*(xx2*xx2*cosh(xx2)/(sinh(xx2))**2)
  700 continue
      tp=t1+(u-u1)*((t2-t1)/(u2-u1))
      if(abs(t-tp).gt.0.005)go to 10
      t=tp
      s=2.*aldm*t
      xxx=piw*t
      if(xxx.gt.15.)go to 99
      s=s+(dels/t)*(xxx*xxx*cosh(xxx)/((sinh(xxx))**2)-xxx/sinh(xxx))
   99 aeff=s/2./t
      ssor=s-sor
      arg=dexp(ssor)
      dens=arg*conmul/((aeff**0.75)*sqrt(t)*(u**1.25+t))
c     dens=dexp(s-sor)*conmul/((aeff**0.75)*sqrt(t)*(u**1.25+t))
      powgkr(ie)=dens
c     powaz(iz,ia,ie)=dens
  850 continue
      cnorm=powgkr(ilow)/exp(float(ilow)*0.1/tsave)
c  now add constant temperature part of level density
      u=float(ilow)*0.1
      do 860 iie=1,ilow
      ie=ilow+1-iie
      powgkr(ie)=cnorm*exp(u/t)
      u=u-0.1
  860 continue
      continue
      continue
      continue
      continue
      do 950 ie=2,601
      powgkr(ie)=powgkr(ie-1)+powgkr(ie)
 950  continue
         do 1953 ie=601,2200
1953   powgkr(ie)=powgkr(600)*gam(ie)/gam(600)
      continue
      continue
      do 1952 ie=1,2200
      gam(ie)=powgkr(ie)*0.1*33.
1952  continue
      return
      end
c
      subroutine ldcalc2
c 5/8/03 make this do a KR level density with no shell effect for o-o nuclei
      common/levop/ldopt,igam
      common/lab11/powaz(15,24,2200)
      common/qq/gow(1100)
      REAL*8 gow
c
      common/krld/powkr(9999),gamkr(9999),powgkr(9999)
      REAL*8 powkr,gamkr,powgkr
      common/pl5/na,nz,title(20)
      common/lab10/pow(4,9999),gam(9999)
cneadb activate the following statement if running on 32-bit word machine
      REAL*8 pow,gam,powaz
       REAL*8 dens,cnorm,arg,ssor
cenadb
      common/ss/sor,rr
      common/pl3/lz,nepr,eb(100),rcsp(100),zee,amass,plex
      common/sft5/exc(18,27),xmxx
      common/lym1/be(18,27,8),pair(18,27),xmas(18,27),delshl(18,27),amas
     1(18,27)
c     common/bym1/symb(18,27),symbp(18,27),symbs(18,27)
      common/atsv/atsave
      common/shft/mc,mp,inver,ike,ipch,kq,qval,ap,at,zp,zt,cld,barfac
 
c we need to provide path to epair in loop
      integer hilim
      pi=3.14159
      pisq=9.8696
      srtpi=1.7725/12.
      tpi=6.283
      con=1./(12.*sqrt(2.))
      gf=6./pisq
       iz=1
       ia=1
      a=atsave+ap+2.-iz-ia
      a=a-1.
c substitute acn-1 for loop over masses 5/8/03
      a3=a**0.33333
      a23=a3*a3
      bofn=0.
      bofp=0.
      aldm=a/9.
      epair=0.
      w=0.185*a3
      piw=pi*w
      em2bar=0.240*a23
      conmul=pi/(12.*sqrt(12.*em2bar))
      umin=2.5+150./a
      delsh=0.
      ud=umin
      if(ud.lt.0.1)umin=-(epair+delsh)+0.1
      if(umin.gt.2.)umin=2.
      if(umin.lt.1.0)umin=1.0
      ilow=10.*umin
      hilim=9999
      if(hilim.lt.1)go to 900
      u=1000.-epair+0.050
      dels=delsh
      if(hilim.lt.ilow)hilim=ilow
cneadb                            letter from the author 24-March-1992
      if((hilim+ilow-1).gt.9999)hilim=hilim-ilow
cneadb
      do 850 iie=ilow,hilim
      ie=hilim+ilow-iie
      u=float(ie)*0.100-.05
      if(u.lt.0.05)go to 850
      tp=sqrt(u/aldm)
   10 t=tp
      t1=t
      t2=t+0.005
      u1=aldm*t1*t1-dels
      u2=aldm*t2*t2 -dels
      xx1=piw*t1
      if(xx1.gt.15.)go to 700
      xx2=piw*t2
c     u1=u1+dels*(xx1*xx1*cosh(xx1)/(sinh(xx1))**2)
c     u2=u2+dels*(xx2*xx2*cosh(xx2)/(sinh(xx2))**2)
  700 continue
      tp=t1+(u-u1)*((t2-t1)/(u2-u1))
      if(abs(t-tp).gt.0.005)go to 10
      t=tp
      s=2.*aldm*t
      xxx=piw*t
      if(xxx.gt.15.)go to 99
c     s=s+(dels/t)*(xxx*xxx*cosh(xxx)/((sinh(xxx))**2)-xxx/sinh(xxx))
   99 aeff=s/2./t
      ssor=s-sor
      arg=dexp(ssor)
      dens=arg*conmul/((aeff**0.75)*sqrt(t)*(u**1.25+t))
c     dens=exp(s-sor)*conmul/((aeff**0.75)*sqrt(t)*(u**1.25+t))
      powkr(ie)=dens
c     powaz(iz,ia,ie)=dens
  850 continue
      cnorm=dens/exp(u/t)
c  now add constant temperature part of level density
      do 860 iie=1,ilow
      ie=ilow+1-iie
      u=float(ie)*0.1-0.05
      powkr(ie)=cnorm*exp(u/t)
  860 continue
  900 continue
      continue
      continue
      return
      end
 
      subroutine ldcalc
      common/atsv/atsave
      common/levop/ldopt,igam
      common/qq/gow(1100)
      REAL*8 gow
      common/lab11/powaz(15,24,2200)
      common/krld/powkr(9999),gamkr(9999),powgkr(9999)
      REAL*8 powkr,gamkr,powgkr,ssor,arg
c
      common/pl5/na,nz,title(20)
      common/lab10/pow(4,9999),gam(9999)
cneadb activate the following statement if running on 32-bit word machine
      REAL*8 pow,gam,powaz
       REAL*8 dens,cnorm
cenadb
      common/ss/sor,rr
      common/pl3/lz,nepr,eb(100),rcsp(100),zee,amass,plex
      common/sft5/exc(18,27),xmxx
      common/lym1/be(18,27,8),pair(18,27),xmas(18,27),delshl(18,27),amas
     1(18,27)
c     common/bym1/symb(18,27),symbp(18,27),symbs(18,27)
      common/shft/mc,mp,inver,ike,ipch,kq,qval,ap,at,zp,zt,cld,barfac
c we need to provide path to epair in loop
      integer hilim
      sor= 7.*sqrt(amass)
      if(amass.gt.50.)sor=sor/2.
      pi=3.14159
      pisq=9.8696
      srtpi=1.7725/12.
      tpi=6.283
      con=1./(12.*sqrt(2.))
      gf=6./pisq
      if(na.gt.22)na=22
      if(nz.gt.13)nz=13
      do 200 iz=1,nz+2
      do 200 ia=1,na+2
      a=atsave+ap+2.-iz-ia
      a3=a**0.33333
      a23=a3*a3
      bofn=be(iz,ia,1)
      bofp=be(iz,ia,2)
c     if(bofn.eq.0..or.bofp.eq.0.) then
c     aldm=a/9.
c     else
      if(bofn.lt.1.)bofn=1.
      if(bofp.lt.1.)bofp=1.
      aldm=0.11035*a-0.25132*a23+1.4*a23/bofn + 0.4*a23/bofp
c     end if
c     epair=pair(iz,ia)/10.
      epair=0.
      w=0.185*a3
      piw=pi*w
      em2bar=0.240*a23
      conmul=pi/(12.*sqrt(12.*em2bar))
      umin=2.5+150./a
      ud=umin-(epair+delshl(iz,ia))-0.1
      if(ud.lt.0.1)umin=-(epair+delshl(iz,ia))+0.1
      if(umin.gt.2.)umin=2.
      if(umin.le.1.0)umin=1.0
      ilow=10.*umin
      if(excit.gt.60.)excit=60.
      hilim=600
 
c jan 02 set for 50 mev max
      if(hilim.lt.1)go to 900
      u=60.-epair+0.050
      dels=delshl(iz,ia)
      if(hilim.lt.ilow)hilim=ilow
cneadb                            letter from the author 24-March-1992
      if((hilim+ilow-1).gt.699)hilim=hilim-ilow
cneadb
      do 850 iie=ilow,hilim
      ie=hilim+ilow-iie
c     u=u-0.100
      u=float(ie)*0.1-0.05
      if(u.lt.0.1)go to 850
      tp=sqrt(u/aldm)
   10 t=tp
      t1=t
      t2=t+0.005
      u1=aldm*t1*t1-dels
      u2=aldm*t2*t2 -dels
      xx1=piw*t1
      if(xx1.gt.15.)go to 700
      xx2=piw*t2
      if(xx2.gt.15.)go to 700
      u1=u1+dels*(xx1*xx1*cosh(xx1)/(sinh(xx1))**2)
      u2=u2+dels*(xx2*xx2*cosh(xx2)/(sinh(xx2))**2)
  700 continue
      tp=t1+(u-u1)*((t2-t1)/(u2-u1))
      if(abs(t-tp).gt.0.005)go to 10
      t=tp
      s=2.*aldm*t
      xxx=piw*t
      if(xxx.gt.15.)go to 99
      s=s+(dels/t)*(xxx*xxx*cosh(xxx)/((sinh(xxx))**2)-xxx/sinh(xxx))
   99 aeff=s/2./t
      ssor=s-sor
      arg=dexp(ssor)
      dens=arg*conmul/((aeff**0.75)*sqrt(t)*(u**1.25+t))
      powaz(iz,ia,ie)=dens
  850 continue
      cnorm=dens/exp(u/t)
c  now add constant temperature part of level density
      do 860 iie=1,ilow
      ie=ilow+1-iie
      u=float(ie)*0.1-.05
      powaz(iz,ia,ie)=cnorm*exp(u/t)
  860 continue
  900 continue
  200 continue
      call ldcalc2
      call ldcalc3
      call ldcalc4
      do 901 iz=1,nz+2
      do 901 ia=1,na+2
      r=powaz(iz,ia,600)/pow(1,600)
      do 901 ie=601,2200
      powaz(iz,ia,ie)=r*pow(1,ie)
 901  continue
c here extend kr table to 200 MeV using kr form for zero shell correction
      continue
      continue
      return
      end
c the ib indices to be used in alice must now be correct thermo values
c  **************************************        08.10.90
* 27.07.2006 
*
*ak july06, original version: A.V.Ignatyuk 08.10.90  ippe
*
*
*       
      subroutine obninsk
      common/pl5/na,nz,title(20)
      common/lab10/pow(4,9999),gam(9999)
cneadb activate the following statement if running on 32-bit word machine
      REAL*8 pow,gam,powaz
      common/qq/gow(1100)
      REAL*8 gow,bb
      REAL*8 ania,sq11,flou, par,res, sq12
      common/krld/powkr(9999),gamkr(9999),powgkr(9999)
      REAL*8 powkr,gamkr,powgkr
      common/lab11/powaz(15,24,2200)
      common/sft5/exc(18,27),xmxx
      common/lym1/be(18,27,8),pair(18,27),xmas(18,27),delshl(18,27),amas
     1(18,27)
c     common/bym1/symb(18,27),symbp(18,27),symbs(18,27)
      common/shft/mc,mp,inver,ike,ipch,kq,qval,ap,at,zp,zt,cld,barfac
      common/sfmp1/sfmp(100,300,5)
      REAL*8 sqcld
          dimension par(10),res(7)
          dimension am(5)
*
* ibsf=0 old data 
*     =1 new data (ripl-2) and "as" from new systematics
*     =2 new data, including "as" (not recommended at high energy)
*
       gamkr(1)=0.
c      gam(1)=0.
      ibsf=1
      maxobn=2200
*
*
         write(*,*)'sfmp=',sfmp(26,56,1)
      if(sfmp(26,56,1).lt.0.0) call readsfm
*
      if(na.gt.22)na=22
      if(nz.gt.13)nz=13
      nia=ifix(ap+at+0.0001e0)
      ania=dfloat(nia)
      iip=ifix(zp+zt+0.0001e0)
       iin=nia-iip
      am(1)=28.
      am(2)=50.
      am(3)=82.
      am(4)=126.
      am(5)=184.
* scaling factor for level desity 
      sq11=-dsqrt(ania*100./9.d+0)
      flou=dexp(sq11)
*
*
* main cycle for iz,ia
      do 9000 iz=1,nz+2
      do 9000 ia=1,na+2
c------------------------------
       iza=iz*ia
      kz=iz
      ja=ia
*   mother nucleus
         iaa=nia-kz-ja+2
         izz=iip-kz+1
         inn=iaa-izz
       if(inn.lt.1.or.izz.lt.3)go to 9000
      a=float(iaa)
      z=float(izz)
c
       x=z
       inx=0
       ind=0
   8   iz1=0
       do 6 i=1,5
       x1=abs(x-am(i))
       if(x1.eq.0.)iz1=3
       if(x1.eq.1.)iz1=2
   6   continue
       ind=ind+iz1
       inx=inx+1
       x=a-z
       if(inx.lt.2)goto 8
       akp=1.d0
       if(ind.eq.6)akp=0.4
       if(ind.eq.5)akp=0.6
       if(ind.eq.4)akp=0.8
c
c------------------------------------------------------------
c   for double magic nuclei zm&nm del0= 0.4*del00
c   for nonmagic nuclei
c        nm+1 or nm-1 and zm  del0=0.6*del00
c        nm+1 or nm-1 and zm+1 orzm-1  del0=0.8*del00
c============================================================
      par(1)=z
       par(2)=a
*
                       par(3)=0.073*a+0.115*a**0.66666666
        if(ibsf.eq.1)  par(3)=0.118*a-0.172*a**0.66666666
         par(4)=.4/a**.333333333
          par(5)=delshl(kz,ja)
           par(6)=akp*12.0/sqrt(a)
            par(7)=.4*akp
             par(8)=30./a**.66666666
               par(9)=0. ! def(ja,kz)
               if(izz.gt.100.or.iaa.gt.300) goto 10 
c for double magic nucleus 208-pb-82, position of first 2+ level
              if(z.eq.82..and.a.eq.208.)par(8)=4.1
c
c redefinition of systematics parameters if they are available in tapesf.dat
c as
        if(ibsf.eq.0) goto 10
        if(ibsf.eq.1) goto 2
        if(sfmp(izz,iaa,1).gt.-999.) par(3)=sfmp(izz,iaa,1)
c shift
    2   if(sfmp(izz,iaa,2).gt.-999.) par(7)=sfmp(izz,iaa,2)
c e2+
        if(sfmp(izz,iaa,3).gt.-999.) par(8)=sfmp(izz,iaa,3)
c dw
        if(sfmp(izz,iaa,4).gt.-999.) par(5)=sfmp(izz,iaa,4)
c beta
        if(sfmp(izz,iaa,5).gt.-999.) par(9)=sfmp(izz,iaa,5)
c
c
   10   continue
c       gam(1)=0.d+00
c
c
c
      do 8000 ib=1,maxobn
c energy bin for level density calculation = 0.1 MeV
      bj=float(ib)
      bj=bj/10.-.05
      res(1)=bj
co
      ipc=1
      call bcs(par,res,ipc)
c
      sq12=res(2)
c ---------------------------------  flou is a scaling factor
c       if(l.ne.7)                              then
      if(iza.eq.1)then
              bb=dexp(sq12)*flou
                   pow(1,ib)=bb
                   pow(2,ib)=bb
                   pow(3,ib)=bb
                   pow(4,ib)=bb
                   powkr(ib)=bb
                   powgkr(ib)=bb
c                                               else
c gamma channel 
        if(ib.le.1100) gow(ib)=bb

c fission channel 
c       if(ib.gt.1)then
c                  gam(ib)=33.d+00*dexp(sq12)*flou+gam(ib-1)
c                  gamkr(ib)=33.d+00*dexp(sq12)*flou+gamkr(ib-1)
c                                               endif
                                                endif

c------------------------------------------------21-mar-1992
                 powaz(iz,ia,ib)=dexp(sq12)*flou
 8000 continue
c
c
c
 9000 continue
c end obninsk ld calc for residual nuclei- next do for fissioning nucleus
c end level density calc; below is for fission widths
      continue
         go to 9501
      nia=ifix(ap+at+0.0001e0)
         a=nia
       gamkr(1)=0.
       gam(1)=0.
      ibsf=1
      maxobn=2200
*
      rewind(80)
*
      if(sfmp(26,56,1).lt.0.0) call readsfm
      nia=ifix(ap+at+0.0001e0)
      ania=dfloat(nia)
      iip=ifix(zp+zt+0.0001e0)
       sqcld=sqrt(cld)
       iz=1
       ia=1
c------------------------------
      kz=iz
      ja=ia
*   mother nucleus
         iaa=nia-kz-ja+2
         izz=iip-kz+1
      a=float(iaa)
      z=float(izz)
c
       x=z
       inx=0
       ind=0
       iz1=0
       do 16 i=1,5
       x1=abs(x-am(i))
       if(x1.eq.0.)iz1=3
       if(x1.eq.1.)iz1=2
  16   continue
       ind=ind+iz1
       inx=inx+1
       x=a-z
       if(inx.lt.2)inx=2
       akp=1.d0
       if(ind.eq.6)akp=0.4
       if(ind.eq.5)akp=0.6
       if(ind.eq.4)akp=0.8
c
c------------------------------------------------------------
c   for double magic nuclei zm&nm del0= 0.4*del00
c   for nonmagic nuclei
c        nm+1 or nm-1 and zm  del0=0.6*del00
c        nm+1 or nm-1 and zm+1 orzm-1  del0=0.8*del00
c============================================================
      par(1)=z
       par(2)=a
*
                       par(3)=0.073*a+0.115*a**0.66666666
        if(ibsf.eq.1)  par(3)=0.118*a-0.172*a**0.66666666
         par(4)=.4/a**.333333333
           par(5)=0.
           par(6)=0.
c set shell /pairing corrections to 0. for fission widths
            par(7)=.4*akp
             par(8)=30./a**.66666666
               par(9)=0. ! def(ja,kz)
               if(izz.gt.100.or.iaa.gt.300) goto 101 
c for double magic nucleus 208-pb-82, position of first 2+ level
c             if(z.eq.82..and.a.eq.208.)par(8)=4.1
c
c redefinition of systematics parameters if they are available in tapesf.dat
c as
        if(sfmp(izz,iaa,1).gt.-999.) par(3)=sfmp(izz,iaa,1)
          ibsf=1
        if(ibsf.eq.0) goto 101
        if(ibsf.eq.1) goto 21
c       if(sfmp(izz,iaa,1).gt.-999.) par(3)=sfmp(izz,iaa,1)
c shift
   21  if(sfmp(izz,iaa,2).gt.-999.) par(7)=sfmp(izz,iaa,2)
c e2+
        if(sfmp(izz,iaa,3).gt.-999.) par(8)=sfmp(izz,iaa,3)
c dw
c       if(sfmp(izz,iaa,4).gt.-999.) par(5)=sfmp(izz,iaa,4)
c beta
        if(sfmp(izz,iaa,5).gt.-999.) par(9)=sfmp(izz,iaa,5)
c
c
  101  continue
c       gam(1)=0.d+00
c
c
c
c
c     if(na.gt.22)na=22
c     if(nz.gt.13)nz=13
c      ipc=0
c       call bcs(par,res,ipc)
c     sq12=res(2)*sqcld
c     gam(1)=dexp(sq12)*flou
c     gamkr(1)=gam(1)
      do 8001 ib=1,maxobn
c energy bin for level density calculation = 0.1 MeV
      bj=float(ib)
      bj=bj/10.-.05
      res(1)=bj
c
       ipc=1
      call bcs(par,res,ipc)
c
      sq12=res(2)*sqcld
      gam(1)=dexp(sq12)*flou
c     gamkr(1)=33.0*dexp(sq12)*flou
      gamkr(1)=gam(1)
c here try to modify so that the fission ld may have modified 'a' parm
c ---------------------------------  flou is a scaling factor
c       if(l.ne.7)                              then
c     if(iza.eq.1)then
              bb=dexp(sq12)*flou
c                  pow(1,ib)=bb
c                  pow(2,ib)=bb
c                  pow(3,ib)=bb
c                  pow(4,ib)=bb
c                  powkr(ib)=bb
c                  powgkr(ib)=bb
c                                               else
c gamma channel 
c       if(ib.le.1100) gow(ib)=bb

c fission channel 
        if(ib.gt.1)then
                   gam(ib)=dexp(sq12)*flou+gam(ib-1)
c we have put the 0.1 meshsize into the 33.0 factor
c                  gamkr(ib)=33.d+00*dexp(sq12)*flou+gamkr(ib-1)
                   gamkr(ib)=gam(ib)
c                                               endif
                                                endif

c------------------------------------------------21-mar-1992
c                powaz(iz,ia,ib)=dexp(sq12)*flou
 8001 continue
c
         do ib=1,maxobn
           gam(ib)=gam(ib)*3.3
           gamkr(ib)=gam(ib)
            enddo
c
c
      continue
 9501 continue
c             write(*,*)' obninsk level densities'
c     write(*,*)' ib gam(ib) pow(1,1,ib  pow(1,2,ib   pow(1,3,ib'
c        do j=1,200,5
c      write(*,*)j, gam(j), powaz(1,1,j), powaz(1,2,j),powaz(1,3,j)
c       enddo
c========redefine pow(*,1) and gow
c            do i1=1,4
c            pow(i1,1)=0.0
c            enddo
c               gow(1)=0.0
      return
      end
c
      subroutine bcs(par,res,ipc)
c
c*****phenomen.description bcs+coll******
c***************************************************************
c  input data:
c      par(1)= z
c      par(2)= a
c      par(3)= a_asymptotic, old systematics: 0.073*a +  0.115*a**(2/3)
c      par(4)= gamma  =  0.4/a**(1/3)
c      par(5)= shell correction
c      par(6)= delta = 12.0/a**(1/2)
c      par(7)= dshift
c      par(8)= position of first 2+ level
c      par(9)=deformation=beta from lysekil formula
c      res(1)=  excitation energy u in mev
c 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  output results:
c      res(2)=  logarithm nat( state density )
c      res(3)=  logarithm nat( level density )
c      res(4)=    sigma**2  /spin cutoff parameter/
c      res(5)=   temperature in mev
c 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit REAL*8 (a-h,o-z)
      dimension par(9),res(7)
      fu(x,dx,gx)=1.+(1.-dexp(-gx*x))*dx/x
      tm2=.24*par(2)**.6666666
      z=par(1)
      a=par(2)
      om2=par(8)
      if(om2.le.0.) om2=30./a**.66666666
      om3=50./a**.66666666
      bet=par(9)
      if(par(9).lt.0.15)bet=0.d0
      u=res(1)
      if(u.gt.0.)go to 11
   10 do 12 i=1,6
   12 res(i+1)=0.d0
      return
   11 ap=par(3)
      dl0=par(6)
      gam=par(4)
c***********************   14-mar-1992 *******************
      dw=par(5)
      nnz=idint(z)
      nz=(nnz/2)*2
      iz=0
      if(nnz.gt.nz) iz=1
      nnnn=idint(a-z)
      nn=((nnnn)/2)*2
c
      in=0
      if(nnnn.gt.nn) in=1
      del=dl0*dfloat(iz+in)+par(7)
      u=u+del
       if(u.le.0.)go to 10
      tkp=0.567*dabs(dl0)
      tkp2=tkp**2
      ax=ap*(1.+dw*gam)
      do 8 lx=1,5
      x=ax*tkp2
      akp=ap*fu(x,dw,gam)
      if(dabs(akp-ax)/akp.le.0.001d0) go to 9
      ax=akp
    8 continue
    9 ecoh=0.152*akp*dl0**2
      ukp=akp*tkp2+ecoh
      skp=2.*akp*tkp
      dtkp=45.84*akp**3*tkp**5
      fkp=0.608*akp*tm2*(1.-0.6667*bet)
      fkp1=0.608*akp*tm2*(1.+0.333*bet)
      if(u.ge.ukp) go to 30
      if(ukp.ne.0.) fi2=1-u/ukp
*mar02
      if(fi2.le.0.d+0) goto 30
      if(fi2.ge.1.d+0) fi2=0.9999999d+0
*
      fi=dsqrt(fi2)
      t=2.*tkp*fi/dlog((fi+1.)/(1.-fi))
      t1=t
      sk=skp*tkp/t*(1-fi2)
      det=dtkp*(1.-fi2)*(1.+fi2)**2
c      det=dtkp*(1.-fi2)**3*(1.+fi2)**2
      fm=fkp*tkp/t*(1.-fi2)
      fm1=fkp1*(1.+2.*tkp/t*(1.-fi2))/3.
      go to 40
   30 ux=u-ecoh
      as=ap*fu(ux,dw,gam)
      if(as.le.0.d0) goto 10
      if(ux.le.0.d0) goto 10
*
      t=dsqrt(ux/as)
      t1=t
      sk=2.*as*t
      det=45.84*as**3*t1**5
      fm=0.608*as*tm2*(1.-.6667*bet)
      fm1=0.608*as*tm2*(1.+0.333333*bet)
   40 sig=fm*t1
      sig1=fm1*t1
      if(dabs(bet).gt.0.d0) sig=fm1**.6666666*fm**.3333333*t1
      res(2)=sk-dlog(dsqrt(det))
      if(u.lt.1.)then
         endif
      res(3)=sk-dlog(dsqrt(6.283*sig*det))
      res(4)=sig
      res(5)=t1
      cga=.0075*a**.33333333
      if(bet.eq.0.d0)go to 51
      q3=1.d0
      om=1.d0
      call qvibr(t1,om,cga,3,q2)
      call qrot(a,bet,sig1,u,qr)
      go to 52
   51 call qvibr(t1,om2,cga,5,q2)
      call qvibr(t1,om3,cga,7,q3)
c      call qrot(a,bet,sig1,u,qr)
      qr=1.d0
   52 res(6)=q2*q3
      res(7)=qr
      col=dlog(q2*q3*qr)
      res(2)=res(2)+col
      res(3)=res(3)+col
      return
      end
c
      subroutine qvibr(t,om,cga,lam,q)
c***** qvibr including damping ***
      implicit REAL*8 (a-h,o-z)
      q=1.d0
      if(t.lt.0.05) goto 10
      gam=cga*(om**2+(2.*3.141593*t)**2)/2.
      fn=dexp(-gam/om)/(dexp(om/t)-1.)
      if(fn.lt.0.) goto 10
      u=lam*om*fn
      s=lam*((1.+fn)*dlog(1.+fn)-fn*dlog(fn))
      q=dexp(s-u/t)
   10 return
      end
c
      subroutine qrot(a,bet,sig,u,qr)
c***** qrot including damping ***
      implicit REAL*8 (a-h,o-z)
      fr(u)=1./(1.+dexp((u-ucr)/dcr))
      ucr=120.*bet*bet*a**.33333333
      dcr=1400.*bet*bet/a**.66666666
      if(bet.le.0.)then
      qr=1.d0
      else
      qr=fr(u)*(sig-1.)+1.
      return
                    endif
      end
*dec05
      subroutine readsfm
      common/sfmp1/sfmp(100,300,5)
      open(80,file='tapesf.dat')
                do 1000 i1=1,100
                do 1000 i2=1,300
                do 1000 i3=1,5
1000            sfmp(i1,i2,i3)=-1000.
c skip comments
      do 1 iii=1,23
 1     read(80,*,err=9999,end=9999)
c        
      do 2000 iii=1,1111111
c                                     as        shift     e2+
      read(80,4,end=2001) iz,name,ia,par3,cas, par7,csh, par8,ce2,
     #                     dw,beta
4     format(i3,1x,a2,i4,    f8.3,a1,  f8.3,a1,  f8.3, a1, 2f8.3)
             if(iz.gt.100.or.ia.gt.300) then
                                        print *,'error. see readsfm'
                                        stop
                                        endif
c 
      sfmp(iz,ia,1)=par3   ! as
      sfmp(iz,ia,2)=par7   ! shift
      sfmp(iz,ia,3)=par8   ! e2+
      sfmp(iz,ia,4)=dw     ! 
      sfmp(iz,ia,5)=beta   ! 
2000  continue
c
2001  close(80)
      return
c
c
9999  write(*,*)' error. file tapesf.dat with superf. model parameters',
     *' is absent or prepared incorrectly.'
      stop
      end
      subroutine gc(a,zz,excit,k)
      common/lab11/powaz(15,24,2200)
      REAL*8 pow,gam,powaz
      common/lab10/pow(4,9999),gam(9999)
cneadb activate the following statement if running on 32-bit word machine
      dimension sz(100),sn(150),pz(100),pn(150)
c
c     data tables of Cook et. al. [At.Data & Nucl. data tables 39,281,(1988)]
c
      data (sz(i1),i1=1,85)/
     1 -.19,-.38,-.57,-.75,-.94,-1.13,-1.32,-1.51,-1.70,
     x -2.10,-2.91,-4.17,-5.72,-7.8,-8.97,-9.7,-10.1,-10.7,-11.38
     1 ,-12.07,-12.55,-13.24,-13.93,-14.71,-15.53,-16.37,-17.36,-18.6,
     2 -18.7,-18.01,-17.87,-17.08,-16.6,-16.75,-16.5,-16.35,-16.22,
     3 -16.41,-16.89,-16.43,-16.68,-16.73,-17.45,-17.29,-17.44,-17.82,
     4 -18.62,-18.27,-19.39,-19.91,-19.14,-18.26,-17.4,-16.42,-15.77,
     5 -14.37,-13.91,-13.1,-13.11,-11.43,-10.89,-10.75,-10.62,-10.41,
     6 -10.21,-9.85,-9.47,-9.03,-8.61,-8.13,-7.46,-7.48,-7.2,-7.13,-7.06
     7 ,-6.78,-6.64,-6.64,-7.68,-7.89,-8.41,-8.49,-7.88,-6.3,-5.47/
      data (sz(il),il=86,100)/-4.78,
     1  -4.37,-4.17,-4.13,-4.32,-4.55,-5.04,-5.28,-6.06,-6.28,-6.87,
     2  -7.20,-7.74,2*0./
      data (sn(i1),i1=1,97)/
     1 0.62,1.24,1.85,2.47,3.09,3.71,4.33,4.95,5.56,6.18,
     x 6.8,7.53,7.55,7.21,7.44,8.07,8.94,9.81,10.6,11.39,
     1 12.54,13.68,14.34,14.19,13.83,13.5,13.,12.13,12.6,13.26,14.13,
     2 14.92,15.52,16.38,17.16,17.55,18.03,17.59,19.03,18.71,18.8,18.99,
     3 18.46,18.25,17.76,17.38,16.72,15.62,14.38,12.88,13.23,13.81,14.9,
     4 14.86,15.76,16.2,17.62,17.73,18.16,18.67,19.69,19.51,20.17,19.48,
     5 19.98,19.83,20.2,19.72,19.87,19.24,18.44,17.61,17.1,16.16,15.9,
     6 15.33,14.76,13.54,12.63,10.65,10.1,8.89,10.25,9.79,11.39,11.72,
     7 12.43,12.96,13.43,13.37,12.96,12.11,11.92,11.,10.8,10.42,10.39/
      data  (sn(il),il=98,150)/
     1  9.69,9.27,8.93,8.57,8.02,7.59,7.33,7.23,7.05,7.42,6.75,6.6,6.38,
     2  6.36,6.49,6.25,5.85,5.48,4.53,4.3,3.39,2.35,1.66,.81,
     3  0.46,-.96,-1.69,-2.53,-3.16,-1.87,-.41,.71,1.66,2.62,3.22,3.76,
     4  4.1,4.46,4.83,5.09,5.18,5.17,5.1,5.01,4.97,5.09,5.03,4.93,5.28,
     5  5.49,5.50,5.37,5.30/
      data pz/0.,5.05,0.,4.50,0.,3.95,0.,3.40,0.,2.90,
     x 0.,2.46,0.,2.09,0.,1.62,0.,1.62,0.,1.83,0.,1.73,
     1 0.,1.35,0.,1.54,0.,1.28,0.26,0.88,0.19,1.35,-.05,1.52,-.09,1.17,
     2 .04,1.24,0.29,1.09,.26,1.17,.23,1.15,-.08,1.35,0.34,1.05,.28,1.27
     3 ,0.,1.05,0.,1.,.09,1.2,.2,1.4,.93,1.,-.2,1.19,.09,.97,0.,.92,.11,
     4 .68,.05,.68,-.22,.79,.09,.69,.01,.72,0.,.4,.16,.73,0.,.46,.17,
     5 .89,0.,.79,0.,.89,0.,.81,-.06,.69,-.2,.71,-.12,.72,0.,.77,2*0./
      data (pn(i1),i1=1,125)/0.,5.25,0.,4.70,0.,4.15,0.,3.58,0.,3.05,
     x 0.,2.67,0.,1.8,0.,1.67,0.,1.86,0.,2.04,0.,1.64,0.,1.44,
     1 0.,1.54,0.,1.3,0.,1.27,0.,1.29,.08,1.41,-.08,1.5,-.05,2.24,-.47,
     2 1.43,-.15,1.44,.06,1.56,.25,1.57,-.16,1.46,0.,.93,.01,.62,-.5,
     3 1.42,.13,1.52,-.65,.8,-.08,1.29,-.47,1.25,-.44,.97,.08,1.65,-.11,
     4 1.26,-.46,1.06,0.22,1.55,-.07,1.37,0.1,1.2,-.27,.92,-.35,1.19,0.,
     5 1.05,-.25,1.61,-.21,.9,-.21,.74,-.38,.72,-.34,.92,-.26,.94,.01,
     6 .65,-.36,.83,.11,.67,.05,1.,.51,1.04,.33,.68,-.27,.81,.09,.75,
     7 .17,.86,.14,1.1,-.22,.84,-.47,.48,.02,.88,.24,.52,.27,.41,-.05/
      data  (pn(il),il=126,150)/
     x  .38,.15,.67,0.,.61,0.,.78,0.,.67,0.,.67,0.,.79,
     1 0.,.6,.04,.64,-.06,.45,.05,.26,-.22,.39,0.0,.39/
c
c     Spin cut off constant by Gilbert & Cameron Can.J.Phys.43,1446,(1965)
c
c     Spin cut off constant by Reffo & Herman Nuovo Cimento 34,261,(1988)
c
      reffo = 0.146
      decon=0.142
      xa = a
      ia = xa
      iz = zz
      n = ia-iz
      if(iz.eq.0)iz=1
      if(n.eq.0)n=1
      epair = pz(iz) + pn(n)
      ss = sz(iz) + sn(n)
      if(ss.ge.-1.0.and.ss.le.1.0)decon=0.12
      ag = xa*(0.00917*ss + decon)
      if(ag.le.0.)ag=xa*(0.1375-0.0000836*xa)
      cons=5.0571*sqrt(reffo)*(xa**0.3333)
      em = epair + 2.25/ag +0.1
      ict = 0
      de = 10.0
      do 100 ir=1,4
      de = de*0.1
      do 200 j=1,200
      u = em - epair
      temp = 1./(sqrt(ag/u) - 1.5/u)
      if(temp.le.0)go to 38
      ict=ict+1
      del2 = em+temp*(alog(cons*sqrt(ag*u**3)/temp)-2.*sqrt(ag*u))
      if(ict.eq.1) sign0 = sign(1.,del2)
      sign2 = sign(1.,del2)
      if(sign2.ne.sign0) go to 45
      del1 = del2
 38   em= em + de
 200   continue
 45   em = em-de
 100   continue
      dxx = abs(del1-del2)
      if(dxx.gt.1.0e-30)go to 250
      emat = 2.25/ag + epair + 0.1
      write(17,11)ia,iz
 11   format(/,'xxxx Constant temp. & Fermi gas level densities do not
     1match for ia =',i3,' iz =',i3,' xxxx',/)
      go to 251
 250  emat = em + de*(del1/(del1-del2))
 251  u = emat - epair
      temp=1./(sqrt(ag/u)-1.5/u)
      ilow = 10.*(emat+0.05)
c     jilim = 10.*excit
      ex=10.*excit
      iex=ex
      bix=iex
      dex=ex-bix
c     u=excit+0.1-epair
      u=300.1+dex/10.-epair
      do 850 iie = ilow,2200
      ie = 2200+ilow-iie
      u=u-0.1
      if(u.le.0.)go to 850
c 2/7/2016 e 5 lines lower undefined;set=u but needs checking
       e=u
      aurt=sqrt(ag*u)
      din=cons*u*aurt
      xnum=exp(2.*aurt)
c     pow(k,ie)=xnum/din
      powaz(iz,ia,ie)=exp(e/temp)/temp
 850  continue
      do 860 iie =1,ilow
      ie=ilow+1-iie
      e=float(ie)*0.1
      u=e - epair
c     pow(k,ie) = exp(e/temp)/temp
      powaz(iz,ia,ie)=exp(e/temp)/temp
 860  continue
      return
      end
 
c     function molnix (iz, in, itype)
       subroutine mlnix (iz, in, itype)
 
c ======================================================================
c
c   This function extracts either the ground-state microscopic
c   correction, for itype = 1; the experimental (or theoretical, if
c   Mex is not measured) mass excess for itype = 2; or the energy
c   shift from the Moller, Nix & Kratz calculated pairing gaps, for
c   itype = 3  from the tables contained in the block data subroutine
c   MOLLNIX. The inputs are the proton and neutron numbers, and the
c   value of itype.
c
c   Called by: DELTAM PRECOF SHELL
c
c   Calls:  MNMACRO
c
c   Written by A. J. Sierk  LANL  T-2  May, 1996.
c   Modifications (AJS); July, October, 1997.
c   Corrected bug (AJS); September, 1998
c
c ======================================================================
 
      integer iz, in, nmin, nmax, nmina, nmaxa, in1, iz1, iz2, iz3, iz4
      integer ina, itype, nma, nmi, nmx
 
      REAL molnix, emx1, emx2, emx3, emx4, emxa, emic1, emic2, emic3,
     &     emic4, emica, epair1, epair2, epair3, epair4, epaira
      REAL eprh1, eprh2, eprl1, eprl2, pair
      REAL a, cevap, cfis, un, z
 
      logical nok
 
      common/shell/shelmol
      common /pairc/  cevap, cfis
      common /emcnix/ emic1(30,27), emic2(35,13), emic3(40,17),
     &                emic4(45,36)
      common /emxnix/ emx1(30,27), emx2(35,13), emx3(40,17), emx4(45,36)
      common /eprnix/ epair1(30,27), epair2(35,13), epair3(40,17),
     &                epair4(45,36)
      common /mxnix2/ emxa(15,7), nmina(7), nmaxa(7)
      common /mcnix2/ emica(15,7), epaira(15,7)
      common /nsmnix/ nmin(93), nmax(93)
 
c ======================================================================
 
      un = REAL(in)
      z = REAL(iz)
      a = z + un
      iz1 = iz - 7
      iz2 = iz - 34
      iz3 = iz - 47
      iz4 = iz - 64
      nok = .true.
      if (iz.lt.8) then
        nok = in.ge.nmina(iz) .and. in.le.nmaxa(iz)
        ina = in - nmina(iz) + 1
        nmi = nmina(iz)
        nma = nmaxa(iz)
        nmx = nmaxa(iz) - nmina(iz) + 1
      else
        nok = in.ge.nmin(iz1) .and. in.le.nmax(iz1)
        in1 = in - nmin(iz1) + 1
        nmi = nmin(iz1)
        nma = nmax(iz1)
        nmx = nmax(iz1) - nmin(iz1) + 1
      endif
      molnix = 0.0
      shelmol=molnix
      if (iz.gt.100) then
        write (16, 1200) iz, in
        return
      endif
 
      if (nok) then
        if (iz.lt.8) then
          if (itype.eq.1) then
            molnix = emica(ina,iz)
          elseif (itype.eq.2) then
            molnix = emxa(ina,iz)
          elseif (itype.eq.3) then
            molnix = epaira(ina,iz)
          endif
        elseif (iz.lt.35 .and. iz.ge.8) then
          if (itype.eq.1) then
            molnix = emic1(in1,iz1)
          elseif (itype.eq.2) then
            molnix = emx1(in1,iz1)
          elseif (itype.eq.3) then
            molnix = epair1(in1,iz1)
          endif
        elseif (iz.lt.48 .and. iz.ge.35) then
          if (itype.eq.1) then
            molnix = emic2(in1,iz2)
          elseif (itype.eq.2) then
            molnix = emx2(in1,iz2)
          elseif (itype.eq.3) then
            molnix = epair2(in1,iz2)
          endif
        elseif (iz.lt.65 .and. iz.ge.48) then
          if (itype.eq.1) then
            molnix = emic3(in1,iz3)
          elseif (itype.eq.2) then
            molnix = emx3(in1,iz3)
          elseif (itype.eq.3) then
            molnix = epair3(in1,iz3)
          endif
        elseif (iz.lt.101 .and. iz.ge.65) then
          if (itype.eq.1) then
            molnix = emic4(in1,iz4)
          elseif (itype.eq.2) then
            molnix = emx4(in1,iz4)
          elseif (itype.eq.3) then
            molnix = epair4(in1,iz4)
          endif
        endif
      shelmol=molnix
      else
c  (nok is false) IN is outside the table for this value of IZ:
 
        if (itype.eq.2) then
c         molnix = mnmacro(iz, in)
         return
        elseif (itype.eq.1) then
c  For nuclei outside the table, set the microscopic correction to that
c  for the first or last nuclide (for this value of iz) in the table.
          if (iz.lt.8) then
            if (in.lt.nmina(iz) .and. in.gt.0) molnix = emica(ina, iz)
            if (in.gt.nmaxa(iz)) molnix = emica(nmx, iz)
          elseif (iz.ge.8 .and. iz.lt.35) then
            if (in.lt.nmin(iz1) .and. in.gt.0) molnix = emic1(1, iz1)
            if (in.gt.nmax(iz1)) molnix = emic1(nmx, iz1)
          elseif (iz.ge.35 .and. iz.lt.48) then
            if (in.lt.nmin(iz1) .and. in.gt.0) molnix = emic2(1, iz2)
            if (in.gt.nmax(iz1)) molnix = emic2(nmx, iz2)
          elseif (iz.ge.48 .and. iz.lt.65) then
            if (in.lt.nmin(iz1) .and. in.gt.0) molnix = emic3(1, iz3)
            if (in.gt.nmax(iz1)) molnix = emic3(nmx, iz3)
          elseif (iz.ge.65 .and. iz.lt.101) then
            if (in.lt.nmin(iz1) .and. in.gt.0) molnix = emic4(1, iz4)
            if (in.gt.nmax(iz1)) molnix = emic4(nmx, iz4)
          endif
        elseif (itype.eq.3) then
          if (iz.lt.8) then
            eprl1 = epaira(1, iz)
            eprl2 = epaira(2, iz)
            eprh1 = epaira(nmx-1, iz)
            eprh2 = epaira(nmx, iz)
          elseif (iz.lt.35 .and. iz.ge.8) then
            eprl1 = epair1(1, iz1)
            eprl2 = epair1(2, iz1)
            eprh1 = epair1(nmx-1, iz1)
            eprh2 = epair1(nmx, iz1)
          elseif (iz.lt.48 .and. iz.ge.35) then
            eprl1 = epair2(1, iz2)
            eprl2 = epair2(2, iz2)
            eprh1 = epair2(nmx-1, iz2)
            eprh2 = epair2(nmx, iz2)
          elseif (iz.lt.65 .and. iz.ge.48) then
            eprl1 = epair3(1, iz3)
            eprl2 = epair3(2, iz3)
            eprh1 = epair3(nmx-1, iz3)
            eprh2 = epair3(nmx, iz3)
          elseif (iz.lt.101 .and. iz.ge.65) then
            eprl1 = epair4(1, iz4)
            eprl2 = epair4(2, iz4)
            eprh1 = epair4(nmx-1, iz4)
            eprh2 = epair4(nmx, iz4)
          endif
          if (in.lt.nmi .and. in.gt.0) then
            if (in.eq.nmi-1 .or. in.eq.nmi-3 .or. in.eq.nmi-5) then
              molnix = eprl2
            elseif (in.eq.nmi-2.or.in.eq.nmi-4 .or. in.eq.nmi-6) then
              molnix = eprl1
            else
              pair = (1.0 - z + 2.0*REAL(iz/2)) +
     &               (1.0 - un + 2.0*REAL(in/2))
              molnix = cevap*pair/sqrt(a)
            endif
          elseif (in.gt.nma) then
            if (in.eq.nma+1 .or. in.eq.nma+3 .or. in.eq.nma+5) then
              molnix = eprh1
            elseif (in.eq.nma+2.or.in.eq.nma+4 .or. in.eq.nma+6) then
              molnix = eprh2
            else
              pair = (1.0 - z + 2.0*REAL(iz/2)) +
     &               (1.0 - un + 2.0*REAL(in/2))
              molnix = cevap*pair/sqrt(a)
            endif
          endif
        endif
      endif
      return
 
c ======================================================================
 
 1200 format (1x,'MOLNIX called with Z,N = ',i3,',',i4,'; outside the ',
     &       'range of the table in MOLLNIX.')
 
c ======================================================================
      end
 
      block data mollnix
 
c ======================================================================
c
c   This block data routine enters the microscopic correction to the
c   macroscopic mass of the spherical macroscopic ground state and the
c   experimental mass excess, (or the calculated one for nuclei without
c   a measured mass excess) according to the Finite-Range Liquid Drop
c   Model of nuclear masses. [P. Moller, J. R. Nix, W. D. Myers, and
c   W. J. Swiatecki, Atomic Data Nucl. Data Tables, 59, 185 (1995)]
c   The less accurate FRLDM is used (instead of the FRDM), because that
c   is the model used to calculate fission barrier heights in the Sierk
c   BARFIT routine which is used in CEM97a.
c
c   The routine also enters the pairing gap (energy shift) to be used
c   in nuclear level density formulas, derived from
c   the tables of Moller, Nix, & Kratz (1996). For odd-Z, even-N nuclei,
c   the energy shift is delta-n from the M, N, & K table;
c   for odd-N, even-Z, delta-p; and for even-even nuclei, the correction
c   is delta-p + delta-n.
c
c   For each value of Z, the values of N from Nmin(Z) to Nmax(Z)
c   are included in the table.
c   The Z index of the arrays nmin, nmax, emic1, epair1, and emx1
c   is equal to Z - 7, since the tables only contain entries for
c   Z,N > 7. For emic2, emx2, & epair2, the index is Z - 34;
c   For emic3, emx3, & epair3, the index is Z - 47; and
c   for emicr, emxr, & epair4, the index is Z - 70.
c   The N index corresponds to N - Nmin(Z) + 1.
c   The table contains Z values from 8 to 100, and N values from the
c   most neutron-deficient isotope in the original table to
c   approximately 5-10  neutrons above the most neutron rich isotope
c   which exists in nature (or has a long half life, for elements with no
c   stable isotopes).  This asymmetry acknowledges the fact that the
c   isotopes encountered in the decay of nuclei will start near the
c   line of beta stability, and move to the neutron deficient region,
c   as evaporation and pre-equilbrium decay proceed.
c
c   Tables added for 1 <= Z <= 7, 11/97
c
c   Written by A. J. Sierk,  LANL  T-2,  April-May, 1996.
c   Modified by AJS, November, 1997.
c   Some numbers updated to 2000 mass tables (Nuclear Wallet Cards),
c   AJS, March, 2002.
c
c ======================================================================
 
      integer nmin, nmina, nmax, nmaxa
 
      REAL a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13,
     &     a14, a15, a16, a17, a18, a19, a20, a21, a22, a23, a24, a25
      REAL b1, b2, b3, b4, b5, b6, b7, b8, b9, b10, b11, b12, b13, b14,
     &     b15, b16, b17, b18, b19, b20, b21, b22, b23, b24, b25, b26,
     &     b27
      REAL c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13,
     &     c14, c15, c16, c17, c18, c19, c20, c21, c22, c23, c24, c25
      REAL emic1, emic2, emic3, emic4, emica, emx1, emx2, emx3, emx4,
     &     emxa, epair1, epair2, epair3, epair4, epaira
 
      dimension a1(30,9), a2(30,6), a3(30,6), a4(30,6), a5(35,4),
     &          a6(35,4), a7(35,4), a8(35,1), a9(40,3), a10(40,4),
     &          a11(40,4), a12(40,4), a13(40,2), a14(45,3), a15(45,3),
     &          a16(45,3), a17(45,3), a18(45,3), a19(45,3), a20(45,3),
     &          a21(45,3), a22(45,3), a23(45,3), a24(45,3), a25(45,3)
      dimension b1(30,9), b2(30,6), b3(30,6), b4(30,5), b5(30,1),
     &          b6(35,3), b7(35,4), b8(35,4), b9(35,2), b10(40,3),
     &          b11(40,3), b12(40,3), b13(40,3), b14(40,3), b15(40,2),
     &          b16(45,3), b17(45,3), b18(45,3), b19(45,3), b20(45,3),
     &          b21(45,3), b22(45,3), b23(45,3), b24(45,3), b25(45,3),
     &          b26(45,3), b27(45,3)
      dimension c1(30,10), c2(30,6), c3(30,6), c4(30,5), c5(35,4),
     &          c6(35,4), c7(35,4), c8(35,1), c9(40,4), c10(40,4),
     &          c11(40,4), c12(40,4), c13(40,1), c14(45,3), c15(45,3),
     &          c16(45,3), c17(45,3), c18(45,3), c19(45,3), c20(45,3),
     &          c21(45,3), c22(45,3), c23(45,3), c24(45,3), c25(45,3)
 
      equivalence (emic1(1, 1),  a1(1,1)), (emic1(1,10),  a2(1,1)),
     &            (emic1(1,16),  a3(1,1)), (emic1(1,22),  a4(1,1))
      equivalence (emic2(1, 1),  a5(1,1)), (emic2(1, 5),  a6(1,1)),
     &            (emic2(1, 9),  a7(1,1)), (emic2(1,13),  a8(1,1))
      equivalence (emic3(1, 1),  a9(1,1)), (emic3(1, 4), a10(1,1)),
     &            (emic3(1, 8), a11(1,1)), (emic3(1,12), a12(1,1)),
     &            (emic3(1,16), a13(1,1))
      equivalence (emic4(1, 1), a14(1,1)), (emic4(1, 4), a15(1,1)),
     &            (emic4(1, 7), a16(1,1)), (emic4(1,10), a17(1,1)),
     &            (emic4(1,13), a18(1,1)), (emic4(1,16), a19(1,1)),
     &            (emic4(1,19), a20(1,1)), (emic4(1,22), a21(1,1)),
     &            (emic4(1,25), a22(1,1)), (emic4(1,28), a23(1,1)),
     &            (emic4(1,31), a24(1,1)), (emic4(1,34), a25(1,1))
 
      equivalence (emx1(1, 1),  b1(1,1)), (emx1(1,10),  b2(1,1)),
     &            (emx1(1,16),  b3(1,1)), (emx1(1,22),  b4(1,1)),
     &            (emx1(1,27),  b5(1,1))
      equivalence (emx2(1, 1),  b6(1,1)), (emx2(1, 4),  b7(1,1)),
     &            (emx2(1, 8),  b8(1,1)), (emx2(1,12),  b9(1,1))
      equivalence (emx3(1, 1), b10(1,1)), (emx3(1, 4), b11(1,1)),
     &            (emx3(1, 7), b12(1,1)), (emx3(1,10), b13(1,1)),
     &            (emx3(1,13), b14(1,1)), (emx3(1,16), b15(1,1))
      equivalence (emx4(1, 1), b16(1,1)), (emx4(1, 4), b17(1,1)),
     &            (emx4(1, 7), b18(1,1)), (emx4(1,10), b19(1,1)),
     &            (emx4(1,13), b20(1,1)), (emx4(1,16), b21(1,1)),
     &            (emx4(1,19), b22(1,1)), (emx4(1,22), b23(1,1)),
     &            (emx4(1,25), b24(1,1)), (emx4(1,28), b25(1,1)),
     &            (emx4(1,31), b26(1,1)), (emx4(1,34), b27(1,1))
 
      equivalence (epair1(1, 1),  c1(1,1)), (epair1(1,11),  c2(1,1)),
     &            (epair1(1,17),  c3(1,1)), (epair1(1,23),  c4(1,1))
      equivalence (epair2(1, 1),  c5(1,1)), (epair2(1, 5),  c6(1,1)),
     &            (epair2(1, 9),  c7(1,1)), (epair2(1,13),  c8(1,1))
      equivalence (epair3(1, 1),  c9(1,1)), (epair3(1, 5), c10(1,1)),
     &            (epair3(1, 9), c11(1,1)), (epair3(1,13), c12(1,1)),
     &            (epair3(1,17), c13(1,1))
      equivalence (epair4(1, 1), c14(1,1)), (epair4(1, 4), c15(1,1)),
     &            (epair4(1, 7), c16(1,1)), (epair4(1,10), c17(1,1)),
     &            (epair4(1,13), c18(1,1)), (epair4(1,16), c19(1,1)),
     &            (epair4(1,19), c20(1,1)), (epair4(1,22), c21(1,1)),
     &            (epair4(1,25), c22(1,1)), (epair4(1,28), c23(1,1)),
     &            (epair4(1,31), c24(1,1)), (epair4(1,34), c25(1,1))
 
      common /nsmnix/ nmin(93), nmax(93)
      common /emcnix/ emic1(30,27), emic2(35,13), emic3(40,17),
     &                emic4(45,36)
      common /emxnix/ emx1(30,27), emx2(35,13), emx3(40,17), emx4(45,36)
      common /eprnix/ epair1(30,27), epair2(35,13), epair3(40,17),
     &                epair4(45,36)
      common /mxnix2/ emxa(15,7), nmina(7), nmaxa(7)
      common /mcnix2/ emica(15,7), epaira(15,7)
 
      data nmina /0,  1,  1,  2,  2,  2,  4/
 
      data nmaxa /5,  8,  8, 10, 12, 14, 15/
 
      data nmin /                             4,  6,  6,
     &            8,  8,  8,  8,  8,  8,  8,  9, 10, 10,
     &           11, 12, 13, 14, 15, 16, 17, 18, 19, 21,
     &           22, 23, 24, 25, 26, 27, 29, 30, 31, 32,
     &           33, 35, 36, 37, 38, 40, 41, 42, 43, 44,
     &           46, 47, 48, 49, 51, 52, 53, 55, 56, 58,
     &           59, 61, 62, 64, 65, 67, 69, 70, 72, 73,
     &           75, 77, 78, 80, 81, 83, 85, 87, 88, 90,
     &           92, 93, 95, 97, 99,100,102,104,106,108,
     &          109,111,113,115,117,119,121,123,125,126/
 
      data nmax /                            16, 18, 20,
     &           22, 22, 24, 24, 26, 27, 28, 29, 31, 32,
     &           33, 34, 38, 39, 40, 41, 42, 44, 44, 50,
     &           51, 52, 53, 54, 58, 59, 62, 62, 63, 64,
     &           66, 68, 69, 70, 72, 74, 75, 81, 82, 83,
     &           85, 86, 87, 88, 90, 91, 92, 94, 95, 97,
     &           98,100,101,103,109,111,113,114,116,117,
     &          119,121,122,124,125,127,129,131,132,134,
     &          136,137,139,141,143,144,146,148,150,152,
     &          153,155,157,159,159,160,160,162,164,165/
 
c   Microscopic corrections to the ground state mass; Consists of
c   shell, pairing, and macroscopic deformation energy. (The deformation
c   of the ground state shape is due to microscopic corrections to the
c   macroscopic potential energy surface.)
c   Non-zero values only for N >= Nmin(Z) and N <= Nmax(Z).
c
c Z = 1 to 7   Difference between measured mass excess and calculated
c              spherical macroscopic values.
      data emica /
     & -40.60,-10.69, -7.41, -9.42, -8.11,-17.42, 9*0.0,
     &  -7.08, -1.73, -4.16, -3.45, -7.22, -9.00,-12.89,-13.46, 7*0.0,
     &  -9.98, -4.05,  0.66,  2.66,  1.08, -1.07, -2.23, -2.89, 7*0.0,
     &  -3.54,  2.84,  2.36,  2.58,  1.84,  0.91,  1.37,  1.35, -0.12,
     &  6*0.0,
     &  -7.51,  1.23,  2.52,  4.01,  2.86,  2.72,  2.38,  2.50,  2.30,
     &   1.85,  1.44, 4*0.0,
     &  -9.16, -0.71,  2.23,  3.28,  1.27,  0.64,  0.40,  1.12,  2.41,
     &   1.86,  1.37, -0.03, -1.07, 2*0.0,
     &   1.35,  3.16,  0.89,  0.48,  0.69,  3.13,  3.57,  3.59,  2.61,
     &   1.64, -0.20, -1.53, 3*0.0/
 
c Z = 8 to 16
      data a1
     & /1.28, 3.03, 1.06, 1.09, 2.05, 4.51, 3.20, 3.82, 2.70, 2.02,
     &  1.14, 0.66, 0.62, 17*0.0,
     &  3.04, 3.14, 3.55, 5.36, 4.94, 6.04, 4.75, 4.51, 3.31, 2.72,
     &  2.80, 2.38, 2.41, 17*0.0,
     &  2.53, 4.00, 2.65, 5.12, 5.16, 5.39, 4.62, 4.83, 3.78, 3.41,
     &  3.46, 3.21, 3.11, 2.68, 2.63, 15*0.0,
     &  2.70, 4.97, 4.68, 5.16, 4.33, 5.64, 4.56, 4.69, 4.14, 3.96,
     &  4.47, 3.99, 3.80, 3.64, 2.98, 15*0.0,
     &  1.90, 4.04, 3.99, 4.23, 3.65, 4.76, 3.73, 4.89, 3.92, 4.73,
     &  4.20, 4.09, 3.78, 3.97, 3.26, 15*0.0,
     &  1.53, 3.81, 3.97, 4.99, 4.43, 5.08, 3.66, 4.41, 4.57, 4.50,
     &  4.15, 4.00, 3.67, 3.67, 3.40, 3.34, 2.84, 13*0.0,
     &  0.27, 2.46, 2.83, 3.76, 2.78, 3.11, 2.42, 2.95, 3.84, 4.21,
     &  3.71, 3.58, 3.09, 3.36, 3.11, 3.23, 2.59, 13*0.0,
     & -0.11, 1.61, 2.96, 3.76, 4.38, 3.67, 2.43, 3.97, 4.10, 4.52,
     &  4.05, 3.85, 3.33, 3.53, 3.28, 3.55, 2.87, 2.74, 2.18,11*0.0,
     &  0.01, 1.82, 2.90, 3.05, 2.93, 3.46, 3.69, 4.15, 4.43, 4.53,
     &  4.12, 3.90, 3.37, 3.63, 3.47, 3.31, 2.85, 2.71, 2.16, 1.84,
     & 10*0.0/
 
c Z = 17 to 22
      data a2
     &/-0.40, 1.09, 2.33, 3.70, 3.86, 3.25, 3.41, 3.90, 4.29, 3.64,
     &  3.63, 3.88, 3.36, 3.45, 3.35, 3.62, 3.18, 3.12, 2.65, 2.26,
     &  1.05, 9*0.0,
     &  0.86, 1.18, 2.63, 2.60, 2.84, 2.29, 3.49, 3.92, 3.37, 3.55,
     &  3.29, 2.75, 2.93, 3.47, 3.79, 3.47, 3.17, 2.41, 1.66, 1.14,
     &  1.13, 9*0.0,
     &  1.39, 2.60, 2.65, 2.64, 2.27, 2.66, 3.13, 3.32, 3.14, 2.78,
     &  2.35, 2.48, 2.75, 3.21, 2.89, 2.49, 1.75, 1.33, 0.44, 0.49,
     &  1.08, 1.71, 8*0.0,
     &  1.53, 2.55, 2.40, 2.32, 1.79, 2.25, 2.69, 2.91, 2.64, 2.29,
     &  2.25, 2.53, 2.51, 2.99, 2.61, 2.33, 1.47, 0.66, 0.11, 0.09,
     &  0.85, 1.49, 2.12, 7*0.0,
     &  1.17, 2.44, 2.26, 1.95, 2.44, 2.95, 3.01, 2.91, 2.46, 2.57,
     &  2.61, 2.76, 3.08, 2.66, 2.27, 1.31, 0.94, 0.20, 0.32, 0.85,
     &  1.66, 2.60, 3.28, 7*0.0,
     &  1.25, 1.89, 1.77, 2.15, 2.60, 2.81, 3.26, 2.88, 2.75, 2.95,
     &  3.03, 3.41, 3.01, 2.67, 1.68, 0.62, 0.23, 0.25, 1.02, 1.72,
     &  2.37, 2.84, 3.15, 7*0.0/
 
c Z = 23 to 28
      data a3
     & /1.42, 1.56, 2.16, 2.37, 3.01, 3.69, 3.32, 3.21, 3.38, 3.45,
     &  3.85, 3.35, 1.93, 1.98, 0.81, 0.53, 0.54, 1.35, 2.17, 2.58,
     &  3.29, 3.54, 4.08, 3.92, 4.46, 4.08, 4*0.0,
     &  1.04, 1.56, 2.04, 2.38, 3.39, 3.00, 2.93, 3.10, 3.21, 3.55,
     &  3.09, 1.45, 1.75, 0.36, 0.12, 0.06, 0.79, 1.61, 2.03, 2.73,
     &  2.97, 3.40, 3.32, 4.52, 4.34, 4.54, 4*0.0,
     &  1.24, 1.48, 2.34, 2.16, 2.54, 2.60, 2.69, 2.82, 2.18, 1.37,
     &  1.19, 0.40, 0.03,-0.28,-0.20, 0.72, 1.53, 2.01, 2.62, 2.88,
     &  3.24, 3.44, 4.15, 3.91, 4.20, 3.94, 4*0.0,
     &  0.84, 1.77, 2.35, 1.93, 1.83, 2.03, 2.06, 2.42, 1.94, 0.66,
     &  0.52,-0.73,-1.02,-0.98, 0.05, 0.75, 1.27, 2.01, 2.17, 2.50,
     &  2.71, 3.13, 2.72, 2.60, 2.24, 2.82, 4*0.0,
     &  1.46, 1.09, 1.28, 0.68, 0.99, 0.67, 0.97, 0.35, 0.12,-0.86,
     & -1.18,-2.18,-1.54,-0.89, 0.12, 0.79, 1.42, 1.79, 2.39, 2.44,
     &  2.90, 2.56, 2.59, 2.23, 2.46, 1.78, 4*0.0,
     &  1.00, 0.44, 0.32, 0.50, 0.47, 0.85, 0.28,-0.10,-1.15,-2.18,
     & -2.74,-2.68,-1.59,-0.68,-0.17, 0.61, 1.11, 1.59, 1.64, 1.77,
     &  1.55, 1.24, 1.02, 1.69, 0.80, 0.99, 0.55, 3*0.0/
 
c Z = 29 to 34
      data a4
     &/ 0.61, 0.11, 0.78, 0.34, 0.99, 0.39, 0.38,-0.62,-0.94,-2.31,
     & -1.20,-0.64, 0.36, 0.90, 1.55, 1.88, 2.43, 2.35, 2.71, 2.37,
     &  2.62, 2.30, 2.53, 1.92, 1.78, 1.04, 4*0.0,
     &  1.04, 1.19, 1.90, 1.13, 1.14, 0.24,-0.29,-1.01,-0.31, 0.38,
     &  1.20, 1.62, 2.30, 2.59, 2.99, 2.98, 3.19, 3.02, 3.37, 2.94,
     &  3.25, 2.59, 2.44, 1.35, 1.01, 0.11,-0.42,-1.39,-1.98,-2.83,
     &  1.67, 2.30, 1.73, 1.77, 1.14, 0.61,-0.20, 0.69, 1.25, 1.96,
     &  2.33, 2.79, 3.02, 3.49, 3.74, 3.94, 3.84, 3.98, 3.77, 3.99,
     &  3.38, 2.87, 2.13, 1.76, 0.89, 0.30,-0.60,-1.35,-1.86,-1.54,
     &  2.49, 1.96, 2.01, 1.49, 1.08, 0.57, 1.14, 1.84, 2.39, 2.72,
     &  3.11, 3.29, 3.88, 3.80, 4.26, 4.21, 4.35, 4.15, 4.26, 3.91,
     &  3.76, 2.56, 2.19, 1.27, 0.71,-0.12,-0.60,-1.10,-0.75,-0.10,
     &  2.39, 2.41, 1.95, 1.75, 1.28, 1.95, 2.37, 2.86, 3.08, 3.42,
     &  3.52, 3.91, 3.87, 4.38, 4.43, 4.53, 4.44, 4.46, 4.16, 4.12,
     &  3.05, 2.64, 1.69, 1.08, 0.41,-0.39,-0.80,-0.43, 0.25, 0.97,
     &  2.38, 1.97, 1.76, 1.63, 2.12, 2.57, 3.02, 3.22, 3.47, 3.52,
     &  3.87, 3.86, 4.34, 4.44, 4.80, 4.49, 4.54, 4.16, 4.15, 3.30,
     &  2.90, 1.93, 1.30, 0.47, 0.02,-0.42, 0.09, 0.58, 1.37, 1.90/
 
c Z = 35 to 38
      data a5
     &/ 2.42, 2.79, 2.17, 2.70, 3.04, 3.51, 3.68, 3.80, 3.67, 3.93,
     &  3.90, 4.25, 4.58, 4.53, 4.73, 4.65, 4.35, 4.28, 4.07, 3.33,
     &  2.30, 1.60, 0.69, 0.35,-0.33, 0.23, 0.85, 1.57, 2.39, 2.87,
     &  3.22, 3.60, 3.87, 2*0.0,
     &  2.69, 2.18, 2.60, 3.00, 3.53, 3.67, 3.80, 3.66, 3.90, 3.84,
     &  4.17, 3.84, 4.13, 4.07, 4.94, 4.47, 4.36, 4.39, 3.77, 2.74,
     &  1.67, 0.96, 0.39,-0.40, 0.26, 0.93, 1.77, 2.55, 3.05, 3.36,
     &  3.76, 3.91, 3.83, 2*0.0,
     &  3.21, 3.41, 4.08, 4.21, 4.17, 4.06, 4.35, 3.97, 3.83, 3.47,
     &  3.43, 3.63, 3.90, 4.71, 4.60, 4.56, 4.02, 2.94, 2.14, 1.13,
     &  0.45,-0.35, 0.31, 1.11, 2.03, 2.80, 3.18, 3.48, 3.75, 3.87,
     &  3.67, 3.68, 3.51, 3.59, 0.0 ,
     &  3.00, 3.68, 3.83, 4.27, 4.03, 4.05, 3.55, 3.34, 3.03, 3.02,
     &  3.09, 3.44, 5.45, 5.22, 4.27, 3.76, 2.70, 1.95, 0.80, 0.05,
     & -0.97,-0.07, 0.82, 1.81, 2.77, 3.37, 3.37, 3.45, 3.63, 3.49,
     &  3.49, 3.37, 3.38, 2*0.0/
 
c Z = 39 to 42
      data a6
     &/ 3.63, 3.88, 4.11, 4.07, 4.00, 3.54, 3.27, 2.84, 2.83, 2.80,
     &  3.19, 3.53, 5.26, 4.30, 3.80, 2.73, 1.95, 0.78, 0.08,-1.19,
     & -0.11, 0.85, 1.89, 2.71, 3.46, 3.66, 3.76, 3.60, 3.47, 3.38,
     &  3.29, 3.32, 3.28, 2*0.0,
     &  3.78, 4.15, 4.00, 3.94, 3.60, 3.31, 2.91, 2.80, 2.79, 3.14,
     &  5.21, 4.98, 4.00, 3.50, 2.42, 1.65, 0.43,-0.38,-1.63,-0.47,
     &  0.46, 1.53, 2.55, 3.20, 3.57, 3.80, 3.70, 3.62, 3.53, 3.43,
     &  3.41, 3.42, 3.52, 2*0.0,
     &  4.06, 4.39, 4.36, 3.97, 3.74, 3.20, 3.06, 3.05, 3.35, 3.93,
     &  5.03, 4.04, 3.52, 2.46, 1.67, 0.46,-0.34,-1.44,-0.57, 0.44,
     &  1.51, 2.46, 3.00, 3.41, 3.81, 3.79, 3.84, 3.66, 3.56, 3.52,
     &  3.53, 3.66, 3.58, 3.64,  0.0,
     &  4.05, 4.41, 4.81, 3.31, 3.26, 4.17, 4.09, 4.84, 4.61, 3.59,
     &  3.09, 1.99, 1.21,-0.04,-0.82,-2.12,-1.07,-0.12, 0.97, 1.80,
     &  2.50, 3.04, 3.51, 3.71, 3.92, 3.96, 3.86, 3.87, 3.81, 3.92,
     &  3.86, 3.92, 3.70, 3.55,  0.0/
 
c Z = 43 to 46
      data a7
     &/ 4.15, 4.42, 4.23, 4.33, 3.80, 3.68, 4.52, 4.27, 3.23, 2.76,
     &  1.62, 0.84,-0.39,-1.16,-2.45,-1.37,-0.40, 0.70, 1.28, 2.08,
     &  2.68, 3.19, 3.57, 3.80, 3.81, 3.89, 3.94, 3.98, 4.01, 3.97,
     &  4.79, 4.57, 4.27, 3.99,  0.0,
     &  4.05, 3.86, 5.24, 4.74, 4.64, 3.77, 3.52, 2.49, 1.98, 0.85,
     &  0.08,-1.20,-1.95,-3.24,-2.08,-1.11, 0.01, 0.59, 1.39, 2.04,
     &  2.58, 3.03, 3.37, 3.57, 3.80, 3.93, 4.00, 4.09, 4.04, 4.36,
     &  4.12, 3.86, 3.57, 3.31,  0.0,
     &  3.38, 4.33, 3.88, 3.89, 3.00, 2.78, 1.71, 1.24, 0.09,-0.60,
     & -1.86,-2.59,-3.86,-2.82,-1.93,-0.99,-0.05, 0.71, 1.41, 1.98,
     &  2.47, 3.00, 3.40, 3.48, 3.78, 3.74, 4.38, 4.26, 4.00, 3.76,
     &  3.48, 3.22, 2.96, 2.74, 2.45,
     &  3.12, 3.13, 2.26, 2.04, 0.96, 0.48,-0.70,-1.44,-2.87,-3.48,
     & -4.87,-3.69,-2.86,-1.55,-0.80,-0.04, 0.65, 1.32, 1.85, 2.41,
     &  2.83, 3.11, 3.38, 3.40, 3.48, 3.87, 3.59, 3.36, 3.08, 2.87,
     &  2.55, 2.40, 2.06, 1.61, 1.07/
 
c Z = 47
      data a8
     &/ 2.03, 1.12, 0.91,-0.17,-0.58,-1.78,-2.52,-3.82,-4.56,-5.67,
     & -4.81,-3.84,-2.72,-1.81,-1.12,-0.28, 0.51, 1.09, 1.69, 2.22,
     &  2.64, 2.94, 3.14, 3.22, 3.46, 3.24, 2.99, 2.70, 2.53, 2.20,
     &  2.09, 1.62, 1.23, 0.63, 0.01/
 
c Z = 48 to 50
      data a9
     & /0.19,-0.04,-1.13,-1.61,-2.81,-3.54,-5.00,-5.64,-7.05,-5.83,
     & -5.10,-3.67,-2.56,-1.86,-0.95,-0.35, 0.32, 0.87, 1.37, 1.77,
     &  2.14, 2.45, 2.55, 2.66, 2.56, 2.55, 2.37, 2.25, 1.88, 1.80,
     &  0.84, 0.35,-0.44,-1.01,-2.42,-3.48,-4.52,-5.74,-7.13,-8.82,
     & -1.60,-2.75,-3.23,-4.44,-5.19,-6.53,-7.33,-8.69,-7.49,-6.49,
     & -5.25,-4.12,-3.15,-2.11,-1.42,-0.72,-0.13, 0.43, 1.01, 1.35,
     &  1.74, 1.84, 2.05, 1.98, 2.09, 1.74, 1.61, 1.07, 0.73, 0.00,
     & -0.46,-1.44,-2.08,-3.33,-4.23,-5.68,-6.83,-8.30,-9.59,-10.58,
     & -4.15,-4.62,-5.85,-6.61,-8.04,-8.63,-10.02,-9.02,-8.09,-6.67,
     & -5.49,-4.36,-3.34,-2.25,-1.51,-0.73,-0.41, 0.13, 0.37, 0.89,
     &  0.81, 1.03, 0.94, 0.96, 0.75, 0.70, 0.19,-0.28,-0.99,-1.53,
     & -2.51,-3.23,-4.36,-5.30,-6.58,-7.67,-9.14,-10.43,-11.55,-10.79/
 
c Z = 51 to 54
      data a10
     &/-4.11,-4.82,-6.10,-6.83,-8.22,-7.10,-6.18,-4.93,-3.82,-2.88,
     & -1.93,-1.10,-0.41, 0.17, 0.73, 1.23, 1.54, 1.86, 1.84, 1.93,
     &  1.71, 1.67, 1.37, 1.23, 0.77, 0.52,-0.15,-0.60,-1.56,-2.13,
     & -3.42,-4.26,-5.87,-6.85,-8.43,-9.60,-10.64,-10.03,-8.77,-7.64,
     & -3.73,-5.26,-5.64,-7.11,-5.89,-5.03,-3.72,-3.00,-1.93,-1.17,
     & -0.40, 0.26, 0.85, 1.35, 1.68, 2.06, 2.10, 2.30, 2.67, 2.47,
     &  2.40, 2.10, 1.91, 1.57, 1.32, 0.68, 0.26,-0.60,-1.19,-2.46,
     & -3.52,-4.74,-5.79,-7.33,-8.53,-9.57,-8.90,-7.79,-6.64,-5.39,
     & -3.73,-4.31,-5.30,-4.59,-3.89,-2.67,-1.91,-1.15,-0.43, 0.19,
     &  0.80, 1.22, 1.73, 1.88, 2.28, 2.36, 2.57, 2.58, 2.62, 2.52,
     &  2.62, 2.45, 2.16, 1.95, 1.47, 1.05, 0.37,-0.30,-1.27,-2.11,
     & -3.43,-4.55,-6.14,-7.27,-8.39,-7.74,-6.61,-4.60,-3.69,-2.79,
     & -3.15,-4.22,-3.39,-2.72,-1.83,-1.15,-0.47, 0.13, 0.77, 1.19,
     &  1.55, 1.90, 2.08, 2.39, 2.56, 2.71, 2.75, 2.77, 2.71, 2.61,
     &  2.44, 2.37, 2.08, 1.88, 1.53, 1.03, 0.49,-0.29,-1.09,-2.37,
     & -3.59,-4.92,-6.10,-7.20,-6.62,-5.35,-3.86,-2.86,-2.16,-1.50/
 
c Z = 55 to 58
      data a11
     &/-2.31,-1.80,-1.11,-0.52,-0.29, 0.16, 0.82, 1.36, 1.68, 1.80,
     &  2.01, 2.30, 2.46, 2.66, 2.67, 2.69, 2.59, 2.57, 2.37, 2.47,
     &  2.21, 2.07, 1.77, 1.40, 0.91, 0.38,-0.25,-1.26,-2.26,-3.90,
     & -5.10,-6.18,-5.67,-4.31,-2.84,-2.12,-1.69,-1.07,-0.67,-0.06,
     & -1.04,-0.38,-0.13, 0.25, 0.93, 1.45, 1.65, 1.80, 1.78, 1.96,
     &  2.23, 2.40, 2.60, 2.66, 2.66, 2.57, 2.53, 2.37, 2.36, 2.21,
     &  2.22, 1.97, 1.78, 1.37, 0.95, 0.40,-0.52,-1.42,-3.01,-4.18,
     & -5.29,-4.66,-3.53,-2.05,-1.45,-1.12,-0.66,-0.24, 0.34, 0.58,
     &  0.03, 0.29, 0.83, 0.91, 1.21, 1.41, 1.39, 1.32, 1.50, 1.74,
     &  2.01, 2.26, 2.32, 2.37, 2.34, 2.36, 2.25, 2.34, 2.24, 2.35,
     &  2.18, 2.09, 1.75, 1.38, 0.87, 0.19,-0.63,-2.22,-3.36,-4.50,
     & -3.76,-2.73,-1.22,-0.76,-0.58,-0.26, 0.24, 0.70, 0.60, 0.75,
     &  1.00, 1.29, 1.36, 1.26, 1.15, 1.04, 1.14, 1.35, 1.56, 1.81,
     &  1.92, 2.05, 2.05, 2.14, 2.09, 2.25, 2.22, 2.34, 2.37, 2.25,
     &  1.97, 1.74, 1.30, 0.71,-0.08,-1.57,-2.72,-3.86,-3.11,-2.07,
     & -0.71,-0.34,-0.13, 0.22, 0.63, 0.88, 0.64, 0.70, 0.44, 0.64/
 
c Z = 59 to 62
      data a12
     &/ 0.97, 0.92, 0.68, 0.52, 0.35, 0.35, 0.55, 0.79, 1.10, 1.26,
     &  1.45, 1.51, 1.63, 1.67, 1.85, 1.95, 2.19, 2.20, 2.26, 2.19,
     &  2.08, 1.64, 1.16, 0.41,-0.98,-2.07,-3.26,-2.51,-1.47,-0.26,
     &  0.19, 0.34, 0.77, 0.86, 0.79, 0.51, 0.57, 0.22, 0.40, 0.21,
     &  0.40, 0.19,-0.05,-0.04, 0.17, 0.42, 0.72, 0.94, 1.14, 1.22,
     &  1.38, 1.41, 1.62, 1.75, 2.03, 2.31, 2.42, 2.30, 2.22, 1.91,
     &  1.52, 0.74,-0.58,-1.73,-2.88,-2.14,-1.04, 0.17, 0.66, 0.75,
     &  1.01, 0.95, 0.78, 0.45, 0.46, 0.16, 0.21, 0.03, 0.20, 0.00,
     & -0.14,-0.43,-0.47,-0.27,-0.03, 0.28, 0.47, 0.72, 0.83, 1.02,
     &  1.05, 1.29, 1.35, 1.70, 1.95, 2.28, 2.16, 2.23, 2.12, 1.78,
     &  1.10,-0.15,-1.24,-2.45,-1.70,-0.61, 0.56, 0.92, 0.95, 1.09,
     &  1.01, 0.80, 0.47, 0.41, 0.05, 0.08,-0.11,-0.03,-0.27,-0.15,
     & -0.47,-0.28,-0.06, 0.33, 0.49, 0.70, 0.85, 0.97, 1.05, 1.15,
     &  1.32, 1.65, 1.86, 2.15, 2.09, 2.19, 2.13, 1.83, 1.18, 0.00,
     & -1.09,-2.29,-1.51,-0.41, 0.72, 1.18, 1.28, 1.41, 1.36, 1.06,
     &  0.72, 0.64, 0.29, 0.26, 0.03, 0.12,-0.19,-0.11,-0.41,-0.31/
 
c Z = 63 to 64
      data a13
     &/-0.28, 0.04, 0.26, 0.36, 0.52, 0.61, 0.67, 0.67, 0.79, 0.98,
     &  1.32, 1.58, 1.93, 1.91, 1.93, 2.03, 1.88, 1.28, 0.21,-0.77,
     & -2.00,-1.13,-0.14, 0.89, 1.22, 1.37, 1.46, 1.36, 1.15, 0.84,
     &  0.70, 0.42, 0.40, 0.10, 0.13,-0.20,-0.20,-0.48,-0.44,-0.71,
     &  0.49, 0.52, 0.61, 0.65, 0.59, 0.56, 0.63, 0.82, 1.16, 1.42,
     &  1.56, 1.72, 1.70, 2.55, 1.69, 1.09, 0.20,-0.90,-2.06,-1.17,
     & -0.16, 0.88, 1.34, 1.54, 1.64, 1.65, 1.42, 1.16, 1.04, 0.77,
     &  0.72, 0.41, 0.38, 0.05, 0.00,-0.32,-0.33,-0.64,-0.50,-0.71/
 
c Z = 65 to 67
      data a14
     & /0.29, 0.31, 0.26, 0.18, 0.08, 0.11, 0.29, 0.68, 0.96, 1.30,
     &  1.25, 1.35, 2.27, 1.56, 1.00, 0.20,-0.65,-1.94,-0.88, 0.01,
     &  0.93, 1.30, 1.46, 1.57, 1.60, 1.39, 1.16, 1.07, 0.82, 0.73,
     &  0.48, 0.36, 0.05,-0.09,-0.41,-0.43,-0.75,-0.71,-0.92,-0.67,
     & -0.72,-0.48,-0.51,-0.18,-0.18,
     &  0.05,-0.09,-0.18,-0.21,-0.05, 0.32, 0.54, 0.80, 0.94, 1.03,
     &  2.01, 1.24, 0.67,-0.01,-1.07,-2.17,-1.32,-0.31, 0.77, 1.27,
     &  1.54, 1.66, 1.73, 1.60, 1.41, 1.31, 1.11, 1.00, 0.73, 0.56,
     &  0.25, 0.06,-0.29,-0.36,-0.70,-0.70,-1.00,-0.76,-0.88,-0.66,
     & -0.75,-0.52,-0.52,-0.19,-0.15,
     & -0.83,-0.90,-0.73,-0.37,-0.10, 0.33, 0.28, 0.47, 1.64, 0.96,
     &  0.42,-0.23,-0.92,-2.06,-1.07,-0.23, 0.69, 1.16, 1.47, 1.55,
     &  1.61, 1.53, 1.38, 1.34, 1.11, 1.00, 0.70, 0.51, 0.18,-0.04,
     & -0.44,-0.51,-0.90,-1.00,-1.27,-1.09,-1.27,-1.15,-1.30,-1.05,
     & -1.08,-0.71,-0.67,-0.28,-0.18/
 
c Z = 68 to 70
      data a15
     &/-1.16,-1.05,-0.62,-0.41,-0.01, 0.08, 0.31, 1.35, 0.65, 0.12,
     & -0.52,-1.24,-2.34,-1.35,-0.51, 0.36, 1.00, 1.43, 1.63, 1.77,
     &  1.80, 1.66, 1.56, 1.38, 1.25, 0.96, 0.76, 0.37, 0.14,-0.28,
     & -0.43,-0.86,-0.91,-1.26,-1.14,-1.37,-1.27,-1.51,-1.34,-1.38,
     & -1.02,-1.03,-0.54,-0.73,-0.45,
     & -1.33,-1.05,-0.59,-0.66,-0.41, 0.80, 0.15,-0.40,-1.01,-1.43,
     & -2.43,-1.48,-0.60, 0.16, 0.80, 1.20, 1.42, 1.61, 1.67, 1.68,
     &  1.56, 1.37, 1.27, 0.93, 0.76, 0.35, 0.09,-0.37,-0.51,-0.98,
     & -1.08,-1.43,-1.37,-1.62,-1.61,-1.96,-1.77,-1.86,-1.56,-1.57,
     & -1.17,-1.23,-0.93,-1.09,-0.84,
     & -1.31,-0.87,-0.88, 0.92, 0.42,-0.29,-0.81,-1.48,-1.87,-2.92,
     & -1.92,-1.04, 0.11, 0.56, 0.98, 1.38, 1.64, 1.75, 1.84, 1.83,
     &  1.63, 1.49, 1.22, 0.99, 0.65, 0.37,-0.06,-0.26,-0.68,-0.82,
     & -1.19,-1.15,-1.55,-1.55,-1.91,-1.84,-1.96,-1.71,-1.77,-1.48,
     & -1.43,-1.24,-1.48,-1.15,-1.36/
 
c Z = 71 to 73
      data a16
     &/-1.27, 0.29,-0.14,-0.84,-1.31,-1.91,-2.29,-3.30,-2.32,-1.42,
     & -0.36, 0.15, 0.59, 1.00, 1.30, 1.49, 1.63, 1.72, 1.64, 1.43,
     &  1.19, 1.06, 0.75, 0.46, 0.09,-0.12,-0.55,-0.65,-1.10,-1.10,
     & -1.49,-1.63,-2.02,-2.00,-2.13,-1.95,-2.07,-1.71,-2.00,-1.65,
     & -1.90,-1.70,-1.93,-1.69,-1.66,
     & -0.38,-1.18,-1.67,-2.35,-2.83,-3.98,-2.89,-2.02,-0.97,-0.11,
     &  0.27, 0.73, 1.12, 1.37, 1.56, 1.72, 1.73, 1.68, 1.46, 1.28,
     &  1.02, 0.78, 0.44, 0.24,-0.12,-0.31,-0.73,-0.79,-1.22,-1.36,
     & -1.79,-1.74,-1.89,-1.80,-1.93,-1.71,-1.73,-1.73,-2.07,-1.94,
     & -2.12,-2.15,-2.12,-2.45,-3.03,
     & -1.86,-2.38,-3.07,-3.50,-4.67,-3.54,-2.72,-1.53,-0.80,-0.32,
     &  0.22, 0.59, 0.91, 1.15, 1.36, 1.43, 1.44, 1.36, 1.17, 0.96,
     &  0.74, 0.46, 0.28,-0.02,-0.22,-0.65,-0.73,-1.15,-1.41,-1.80,
     & -1.84,-2.07,-1.90,-2.19,-2.03,-2.25,-2.31,-2.48,-2.44,-2.66,
     & -2.82,-2.57,-2.93,-3.55,-3.97/
 
c Z = 74 to 76
      data a17
     &/-3.67,-4.22,-5.24,-4.26,-3.44,-2.09,-1.28,-0.66,-0.12, 0.29,
     &  0.66, 0.94, 1.19, 1.31, 1.38, 1.38, 1.28, 1.14, 1.00, 0.77,
     &  0.62, 0.38, 0.14,-0.20,-0.40,-0.79,-1.11,-1.50,-1.55,-1.82,
     & -1.79,-2.01,-2.04,-2.32,-2.31,-2.61,-2.57,-3.21,-3.29,-3.01,
     & -3.40,-4.07,-4.60,-5.45,-5.98,
     & -5.10,-6.17,-5.10,-4.31,-3.06,-2.10,-1.51,-0.85,-0.35, 0.07,
     &  0.42, 0.79, 0.91, 1.06, 1.10, 1.05, 0.93, 0.88, 0.62, 0.49,
     &  0.19, 0.01,-0.33,-0.58,-0.93,-1.23,-1.66,-1.84,-2.06,-2.02,
     & -2.27,-2.23,-2.57,-2.55,-2.94,-3.22,-3.85,-3.94,-4.13,-4.05,
     & -4.69,-5.32,-6.15,-6.68,-7.61,
     & -6.10,-5.30,-3.86,-2.81,-1.92,-1.22,-0.69,-0.21, 0.14, 0.50,
     &  0.81, 0.96, 1.01, 1.07, 1.02, 0.92, 0.93, 0.76, 0.55, 0.34,
     &  0.09,-0.26,-0.54,-0.88,-1.25,-1.51,-1.71,-1.73,-1.98,-1.90,
     & -2.24,-2.30,-2.89,-3.29,-3.94,-4.18,-4.61,-4.60,-5.40,-6.08,
     & -6.95,-7.52,-8.63,-8.82,-8.02/
 
c Z = 77 to 79
      data a18
     &/-5.21,-4.04,-2.96,-2.07,-1.37,-0.77,-0.36, 0.12, 0.45, 0.75,
     &  0.96, 1.17, 1.19, 1.21, 0.85, 1.02, 0.64, 0.35, 0.16,-0.17,
     & -0.29,-0.69,-1.00,-1.37,-1.56,-1.39,-1.55,-1.65,-2.16,-2.44,
     & -3.02,-3.51,-4.10,-4.57,-5.07,-5.32,-6.42,-6.98,-7.86,-8.45,
     & -9.50,-9.69,-8.74,-7.86,-6.54,
     & -3.54,-2.62,-1.67,-0.89,-0.52, 0.07, 0.47, 0.80, 0.93, 1.16,
     &  1.12, 1.25, 1.22, 1.08, 1.03, 0.92, 0.66, 0.51, 0.21,-0.05,
     & -0.34,-0.47,-0.76,-0.15,-0.50,-0.95,-1.42,-1.98,-2.58,-3.28,
     & -4.00,-4.75,-5.36,-6.05,-7.11,-7.73,-8.60,-9.17,-10.25,-10.40,
     & -9.58,-8.54,-7.19,-6.02,-4.91,
     & -4.03,-3.03,-2.18,-1.36,-0.71,-0.06, 0.41, 0.86, 1.17, 1.45,
     &  1.58, 1.67, 1.64, 1.58, 1.42, 1.25, 0.98, 0.75, 0.46, 0.18,
     & -0.15,-0.39,-0.79,-1.17,-1.65,-2.14,-2.76,-3.33,-4.05,-4.75,
     & -5.53,-6.26,-6.90,-7.98,-8.59,-9.44,-10.06,-11.14,-11.39,-10.41,
     & -9.47,-8.24,-6.93,-5.77,-4.65/
 
c Z = 80 to 82
      data a19
     &/-2.79,-1.91,-1.23,-0.56,-0.04, 0.43, 0.77, 1.11, 1.30, 1.44,
     &  1.47, 1.45, 1.31, 1.18, 0.90, 0.66, 0.42, 0.11,-0.24,-0.55,
     & -0.99,-1.41,-1.94,-2.43,-3.04,-3.67,-4.49,-5.09,-5.96,-6.72,
     & -7.49,-8.34,-9.07,-9.97,-10.65,-11.82,-12.00,-11.13,-10.13,-8.91,
     & -7.60,-6.44,-5.31,-4.28,-3.37,
     & -1.99,-1.20,-0.64,-0.04, 0.32, 0.69, 0.87, 1.12, 1.11, 1.22,
     &  1.06, 1.04, 0.80, 0.60, 0.31,-0.01,-0.41,-0.80,-1.32,-1.79,
     & -2.43,-2.91,-3.72,-4.38,-5.13,-5.87,-6.64,-7.49,-8.27,-9.24,
     & -9.97,-10.82,-11.57,-12.58,-12.82,-11.89,-10.89,-9.59,-8.40,
     & -7.34,-6.16,-5.15,-4.12,-3.12,-2.33,
     & -1.50,-0.96,-0.44,-0.07, 0.34, 0.48, 0.73, 0.74, 0.84, 0.80,
     &  0.80, 0.61, 0.56, 0.21,-0.02,-0.44,-0.70,-1.33,-1.69,-2.43,
     & -2.88,-3.68,-4.27,-5.08,-5.75,-6.59,-7.40,-8.22,-9.14,-10.02,
     & -10.99,-11.82,-12.68,-12.84,-11.97,-10.92,-9.65,-8.50,-7.24,
     & -6.20,-5.15,-4.28,-3.46,-2.68,-1.98/
 
c Z = 83 to 85
      data a20
     &/ 0.62, 1.01, 1.39, 1.83, 1.84, 2.02, 2.12, 1.96, 1.91, 1.69,
     &  1.55, 1.23, 0.96, 0.53, 0.17,-0.40,-0.82,-1.49,-1.98,-2.77,
     & -3.42,-4.24,-4.94,-5.75,-6.57,-7.39,-8.28,-9.08,-9.98,-10.77,
     & -11.70,-11.95,-11.07,-9.99,-8.76,-7.64,-6.34,-5.41,-4.36,-3.55,
     & -2.83,-1.98,-1.30,-0.69,-0.06,
     &  1.33, 1.32, 1.17, 1.12, 0.95, 1.03, 0.98, 1.25, 1.35, 1.42,
     &  1.40, 0.88, 0.62, 0.45,-0.01,-0.50,-1.11,-1.73,-2.31,-3.15,
     & -3.80,-4.64,-5.43,-6.26,-7.13,-8.07,-9.01,-9.80,-10.71,-10.84,
     & -9.99,-8.97,-7.71,-6.57,-5.34,-4.39,-3.69,-2.91,-1.96,-1.13,
     & -0.54, 0.03, 0.20, 0.56, 0.55,
     &  0.90, 0.70, 0.57, 0.65, 0.71, 0.86, 1.04, 1.43, 1.16, 0.94,
     &  0.67, 0.46, 0.28,-0.03,-0.36,-0.87,-1.50,-2.15,-2.88,-3.59,
     & -4.56,-5.07,-6.22,-6.94,-7.77,-8.52,-9.44,-9.58,-8.69,-7.76,
     & -6.52,-5.37,-4.25,-3.60,-3.01,-2.18,-1.82,-1.16,-0.55, 0.02,
     &  0.10, 0.45, 0.41, 0.71, 0.64/
 
c Z = 86 to 88
      data a21
     &/ 0.46, 0.32, 0.43, 0.45, 0.67, 0.95, 1.21, 1.34, 1.36, 1.10,
     &  0.87, 0.60, 0.42, 0.06,-0.22,-0.69,-1.45,-2.13,-2.90,-3.72,
     & -4.26,-5.20,-5.83,-6.66,-7.45,-8.35,-8.46,-7.56,-6.58,-5.30,
     & -4.30,-3.51,-2.71,-2.35,-1.72,-1.46,-0.65,-0.18, 0.32, 0.33,
     &  0.60, 0.50, 0.74, 0.61, 0.75,
     &  0.02, 0.07, 0.27, 0.36, 0.65, 1.04, 1.05, 1.37, 1.16, 0.87,
     &  0.63, 0.35, 0.03,-0.37,-0.83,-1.38,-2.28,-3.06,-3.67,-4.32,
     & -4.90,-5.65,-6.39,-7.27,-7.33,-6.44,-5.55,-4.29,-3.23,-3.00,
     & -2.51,-2.22,-1.74,-1.26,-0.69,-0.37,-0.02, 0.05, 0.35, 0.28,
     &  0.47, 0.31, 0.44, 0.21, 0.39,
     &  0.19, 0.45, 0.77, 1.04, 1.18, 1.14, 1.57, 1.30, 1.16, 0.83,
     &  0.61, 0.21,-0.26,-0.81,-1.62,-2.32,-2.91,-3.51,-3.99,-4.74,
     & -5.44,-6.24,-6.36,-5.53,-4.59,-3.39,-2.34,-2.55,-1.95,-1.86,
     & -1.44,-0.85,-0.31,-0.35, 0.03, 0.07, 0.36, 0.21, 0.31, 0.13,
     &  0.25, 0.03, 0.23,-0.11, 0.06/
 
c Z = 89 to 91
      data a22
     &/ 0.39, 0.54, 0.72, 0.87, 1.30, 1.64, 1.49, 1.20, 0.95, 0.56,
     &  0.12,-0.40,-0.97,-1.55,-2.35,-2.85,-3.21,-3.86,-4.50,-5.32,
     & -5.48,-4.60,-3.67,-2.71,-1.60,-2.20,-1.79,-1.75,-1.10,-0.79,
     & -0.51,-0.56,-0.24,-0.25, 0.01,-0.21,-0.12,-0.28,-0.17,-0.37,
     & -0.17,-0.48,-0.30,-0.55,-0.35,
     &  0.76, 0.85, 1.11, 1.55, 1.87, 1.62, 1.42, 1.04, 0.61, 0.07,
     & -0.56,-1.15,-1.72,-2.24,-2.55,-3.22,-3.75,-4.56,-4.66,-3.84,
     & -2.92,-2.04,-1.14,-1.65,-1.27,-1.22,-0.38,-0.49, 0.07,-0.18,
     & -0.04,-0.28,-0.17,-0.41,-0.34,-0.53,-0.38,-0.58,-0.38,-0.71,
     & -0.52,-0.75,-0.54,-0.89,-0.60,
     &  0.51, 0.80, 1.31, 1.52, 1.43, 1.43, 1.07, 0.89, 0.34,-0.23,
     & -0.77,-1.31,-1.67,-1.96,-2.55,-3.05,-3.79,-3.91,-3.03,-2.14,
     & -0.83, 0.16,-1.41,-0.93,-0.75,-0.39,-0.17,-0.25,-0.51,-0.38,
     & -0.65,-0.56,-0.80,-0.72,-0.97,-0.82,-1.01,-0.78,-1.14,-0.94,
     & -1.15,-0.93,-1.08,-1.01,-0.94/
 
c Z = 92 to 94
      data a23
     &/ 1.24, 1.44, 1.41, 1.53, 1.16, 1.05, 0.73, 0.14,-0.39,-0.90,
     & -1.32,-1.49,-2.07,-2.52,-3.26,-3.33,-2.39,-1.61,-1.23,-0.46,
     & -0.84,-0.37,-0.13, 0.19,-0.16,-0.21,-0.54,-0.45,-0.77,-0.74,
     & -1.04,-1.00,-1.22,-1.06,-1.25,-1.01,-1.36,-1.16,-1.37,-1.16,
     & -1.32,-1.22,-1.18,-0.54,-0.23,
     &  1.00, 1.14, 1.07, 0.95, 0.70, 0.38,-0.14,-0.58,-0.99,-1.09,
     & -1.59,-1.92,-2.52,-2.68,-1.71,-0.98, 0.31,-0.66,-0.38, 0.14,
     &  0.06, 0.00,-0.36,-0.44,-0.78,-0.75,-1.12,-1.15,-1.45,-1.46,
     & -1.64,-1.52,-1.69,-1.62,-1.84,-1.63,-1.84,-1.62,-1.78,-1.70,
     & -1.62,-1.00,-0.84,-0.30,-0.31,
     &  0.91, 0.91, 0.60, 0.48, 0.07,-0.35,-0.73,-0.81,-1.30,-1.63,
     & -2.28,-2.34,-1.41,-0.64,-0.61,-0.19,-0.06, 0.36, 0.35, 0.23,
     & -0.12,-0.25,-0.59,-0.62,-1.00,-1.07,-1.41,-1.44,-1.68,-1.56,
     & -1.90,-1.72,-1.95,-1.77,-2.02,-1.85,-2.25,-1.97,-1.91,-1.26,
     & -1.09,-0.59,-0.62,-0.34,-0.47/
 
c Z = 95 to 97
      data a24
     &/ 0.43, 0.25, 0.25,-0.13,-0.42,-0.52,-0.94,-1.18,-1.70,-1.88,
     & -0.79,-0.19, 1.06, 0.26, 0.43, 0.55, 0.23, 0.07,-0.28,-0.41,
     & -0.75,-0.86,-1.25,-1.37,-1.72,-1.79,-2.05,-1.96,-2.31,-2.17,
     & -2.41,-2.24,-2.50,-2.32,-2.73,-2.43,-2.37,-1.70,-1.51,-1.00,
     & -1.03,-0.75,-0.91, 2*0.0,
     &  0.05,-0.02,-0.32,-0.44,-0.84,-1.09,-1.67,-1.76,-0.78,-0.35,
     &  1.19, 0.41, 0.70, 0.76, 0.48, 0.31,-0.02,-0.17,-0.52,-0.63,
     & -1.04,-1.19,-1.53,-1.64,-1.92,-2.02,-2.29,-2.19,-2.63,-2.53,
     & -2.84,-2.74,-2.97,-2.77,-2.69,-2.04,-1.86,-1.47,-1.33,-1.07,
     & -1.22,-1.07, 3*0.0,
     & -0.20,-0.28,-0.67,-0.84,-1.28,-1.48,-0.38,-0.25, 1.34, 0.82,
     &  0.60, 0.47, 0.30, 0.12,-0.20,-0.34,-0.70,-0.85,-1.25,-1.43,
     & -1.75,-1.91,-2.37,-2.35,-2.65,-2.77,-3.09,-3.02,-3.34,-3.27,
     & -3.53,-3.33,-3.26,-2.84,-2.41,-2.03,-1.85,-1.58,-1.72,-1.57,
     & 5*0.0/
 
c Z = 98 to 100
      data a25
     &/-0.74,-0.94,-1.37,-1.52,-0.54,-0.24, 1.28, 0.83, 0.67, 0.61,
     &  0.38, 0.26, 0.03,-0.14,-0.48,-0.64,-1.04,-1.24,-1.56,-1.76,
     & -2.24,-2.27,-2.77,-2.77,-3.10,-3.07,-3.42,-3.39,-3.68,-3.65,
     & -3.46,-3.04,-2.62,-2.24,-2.27,-1.85,-2.00,-1.86,-2.07,-1.95,
     & 5*0.0,
     & -1.23,-1.34,-0.37,-0.27, 1.36, 0.49, 0.34, 0.29, 0.09, 0.04,
     & -0.24,-0.39,-0.71,-0.91,-1.29,-1.55,-1.89,-2.14,-2.50,-2.54,
     & -3.10,-3.11,-3.47,-3.47,-3.85,-3.86,-4.34,-4.16,-3.96,-3.56,
     & -3.14,-2.76,-2.80,-2.50,-2.52,-2.39,-2.64,-2.51,-2.52,-2.03,
     & 5*0.0,
     & -1.62,-0.66,-0.31, 1.13, 0.46, 0.29, 0.27, 0.05, 0.01,-0.14,
     & -0.28,-0.56,-0.78,-1.14,-1.32,-1.67,-1.92,-2.29,-2.36,-2.95,
     & -2.99,-3.38,-3.41,-3.81,-3.87,-4.39,-4.26,-4.06,-3.70,-3.56,
     & -2.95,-3.01,-2.76,-2.96,-2.76,-3.00,-2.90,-2.92,-2.43,-2.15,
     & 5*0.0/
 
c ======================================================================
c
c   Mass excess table; Experimental value if measured; Moller, Nix,
c   Myers, & Swiatecki value if unmeasured. Non-zero values only
c   for N >= Nmin(Z) and N <= Nmax(Z).
c
c Z = 1 to 7   Measured values only.
      data emxa /
     &   7.29, 13.14, 14.95, 25.90, 36.80, 41.90, 9*0.0,
     &  14.93,  2.43, 11.39, 17.59, 26.11, 31.60, 40.82, 48.81, 7*0.0,
     &  25.30, 11.68, 14.09, 14.91, 20.95, 24.95, 33.05, 40.80, 7*0.0,
     &  18.38, 15.77,  4.94, 11.35, 12.61, 20.17, 25.08, 35.16, 39.90,
     &  6*0.0,
     &  27.87, 22.92, 12.42, 12.05,  8.67, 13.37, 16.56, 23.66, 28.97,
     &  37.08, 43.70, 4*0.0,
     &  35.09, 28.91, 15.70, 10.65,  0.00,  3.13,  3.02,  9.87, 13.69,
     &  21.04, 24.92, 32.80, 37.60, 2*0.0,
     &  25.30, 17.34,  5.35,  2.86,  0.10,  5.68,  7.87, 13.12, 15.86,
     &  21.77, 25.23, 32.10, 3*0.0/
 
c Z = 8 to 16  Measured values, or Moller-Nix calculations.
      data b1
     & / 32.06, 23.11,  8.01,  2.86, -4.74, -0.81, -0.78,  3.33,  3.80,
     &    8.06,  9.28, 14.65, 19.00, 17*0.0,
     &   16.80, 10.68,  1.95,  0.87, -1.49, -0.02, -0.05,  2.79,  3.33,
     &    7.54, 11.27, 18.3 , 25.0 , 17*0.0,
     &   23.99, 16.49,  5.32,  1.75, -7.04, -5.73, -8.03, -5.16, -5.95,
     &   -2.06,  0.44,  6.96, 11.12, 18.00, 22.20, 15*0.0,
     &   12.93,  6.84, -2.19, -5.18, -9.53, -8.42, -9.36, -6.90, -5.60,
     &   -1.14,  2.66,  8.29, 12.7 , 18.3 , 26.  , 15*0.0,
     &   17.57, 10.91, -0.40, -5.47,-13.93,-13.19,-16.21,-14.59,-15.02,
     &  -10.66, -9.07, -3.35, -0.82,  5.09,  8.5 , 15*0.0,
     &   25.70, 17.60,  6.77, -0.05, -8.91,-12.21,-17.20,-16.85,-18.22,
     &  -15.87,-14.97,-11.08, -8.59, -3.21, -0.32, 5.9 ,  9.6 , 13*0.0,
     &   31.89, 23.11, 10.75,  3.83, -7.14,-12.39,-21.49,-21.90,-24.43,
     &  -22.95,-24.08,-20.49,-19.96,-14.36,-12.48, -6.5 ,-3.7 , 13*0.0,
     &   43.99, 32.58, 19.72, 10.97, -0.75, -7.16,-16.95,-20.20,-24.44,
     &  -24.31,-26.34,-24.56,-24.86,-20.25,-18.99,-14.5 ,-12.6 , -8.3 ,
     &   -4.8 , 11*0.0,
     &   53.68, 41.85, 26.63, 16.71,  4.1 , -3.16,-14.06,-19.05,-26.02,
     &  -26.59,-29.93,-28.85,-30.66,-26.90,-26.86,-23.16,-22.8 ,-18.6 ,
     &  -17.2 ,-12.5 , 10*0.0/
 
c Z = 17 to 22
      data b2
     & / 67.57, 53.42, 37.86, 27.22, 14.11,  4.52, -7.06,-13.33,-21.00,
     &  -24.44,-29.01,-29.52,-31.76,-29.80,-29.80,-27.56,-27.40,-24.69,
     &  -24.0 ,-20.0 ,-18.9 , 9*0.0,
     &   64.17, 45.74, 34.72, 19.52, 10.32, -2.18, -9.38,-18.38,-23.05,
     &  -30.23,-30.95,-34.72,-33.24,-35.04,-33.07,-34.42,-31.97,-32.26,
     &  -29.72,-29.72,-25.90, 9*0.0,
     &   59.57, 46.47, 30.91, 19.69,  6.74, -1.42,-11.17,-17.42,-24.80,
     &  -28.80,-33.81,-33.53,-35.56,-35.02,-36.59,-35.81,-36.61,-35.42,
     &  -35.69,-32.12,-30.32,-25.4 , 8*0.0,
     &   70.57, 56.86, 39.33, 27.62, 12.79,  4.44, -6.44,-13.16,-22.06,
     &  -27.28,-34.85,-35.14,-38.55,-38.41,-41.47,-40.81,-43.13,-42.34,
     &  -44.21,-41.29,-39.57,-35.90,-32.5 , 7*0.0,
     &   68.99, 52.47, 38.98, 23.97, 13.85,  2.31, -5.87,-14.17,-20.53,
     &  -28.64,-32.12,-36.19,-37.82,-41.07,-41.76,-44.33,-44.49,-46.55,
     &  -44.54,-43.22,-40.5 ,-38.0 ,-34.0 , 7*0.0,
     &   61.75, 48.69, 32.23, 21.61,  8.43,  0.01, -8.9 ,-15.71,-25.12,
     &  -29.32,-37.55,-39.01,-44.13,-44.93,-48.49,-48.56,-51.43,-49.73,
     &  -49.46,-46.8 ,-45.6 ,-41.7 ,-39.1 , 7*0.0/
 
c Z = 23 to 28
      data b3
     & / 61.32, 44.76, 32.81, 19.01,  9.51, -1.44, -9.82,-19.46,-23.72,
     &  -31.87,-37.07,-42.00,-44.47,-47.96,-49.22,-52.20,-51.44,-51.85,
     &  -49.89,-49.1 ,-46.2 ,-44.3 ,-40.3 ,-37.9 ,-33.1 ,-30.79, 4*0.0,
     &   54.42, 42.03, 27.00, 16.83,  4.74, -4.01,-13.55,-19.41,-29.47,
     &  -34.55,-42.82,-45.33,-50.26,-51.45,-55.41,-55.28,-56.93,-55.10,
     &  -55.29,-52.39,-51.9 ,-47.8 ,-46.8 ,-42.8 ,-41.2 ,-35.49, 4*0.0,
     &   54.53, 38.91, 27.82, 14.21,  4.81, -6.42,-14.17,-23.41,-29.20,
     &  -37.61,-42.62,-48.24,-50.70,-54.69,-55.55,-57.71,-56.91,-57.49,
     &  -55.90,-55.47,-52.8 ,-51.6 ,-48.5 ,-46.8 ,-43.1 ,-40.9 , 4*0.0,
     &   48.27, 36.93, 22.64, 12.13, -0.64, -8.60,-18.12,-24.58,-34.47,
     &  -40.22,-48.33,-50.94,-56.25,-57.47,-60.60,-60.18,-62.15,-60.66,
     &  -61.41,-58.92,-58.90,-55.5 ,-54.9 ,-51.3 ,-50.3 ,-46.6 , 4*0.0,
     &   49.24, 33.70, 22.44,  8.85, -0.32,-11.67,-19.06,-28.93,-34.35,
     &  -42.64,-48.01,-54.02,-56.04,-59.34,-59.84,-62.22,-61.65,-62.90,
     &  -61.43,-61.84,-59.79,-59.16,-56.1 ,-55.3 ,-51.8 ,-51.0 , 4*0.0,
     &   43.50, 31.18, 16.75,  7.14, -5.21,-12.81,-22.91,-29.38,-39.21,
     &  -45.33,-53.90,-56.08,-60.22,-61.15,-64.47,-64.22,-66.74,-65.51,
     &  -67.10,-65.12,-66.03,-63.74,-63.49,-60.4 ,-59.5 ,-55.9 ,-54.7 ,
     &   3*0.0/
 
c Z = 29 to 33
      data b4
     & / 43.85, 28.74, 18.35,  5.30, -3.28,-14.68,-22.28,-32.43,-38.55,
     &  -47.31,-51.66,-56.35,-58.34,-61.98,-62.80,-65.58,-65.42,-67.26,
     &  -66.25,-67.30,-65.54,-65.74,-62.96,-62.76,-60.15,-59.49, 4*0.0,
     &   28.17, 14.46,  5.66, -7.14,-14.98,-25.95,-32.70,-42.29,-47.26,
     &  -54.18,-56.34,-61.17,-62.21,-66.00,-65.91,-68.90,-67.88,-70.00,
     &  -68.42,-69.56,-67.32,-68.13,-65.41,-65.71,-62.47,-62.0 ,-58.6 ,
     &  -57.2 ,-54.17,-52.16,
     &   27.10, 17.01,  4.15, -4.85,-16.09,-24.12,-34.01,-39.82,-46.31,
     &  -52.00,-56.69,-58.84,-62.65,-63.72,-66.88,-67.08,-69.32,-68.91,
     &  -70.14,-68.59,-69.70,-68.05,-68.46,-66.20,-65.87,-63.66,-62.5 ,
     &  -59.1 ,-58.0 ,-52.80,
     &   26.77, 12.77,  3.53, -8.76,-16.93,-27.77,-33.42,-41.49,-45.97,
     &  -54.43,-56.41,-61.62,-62.65,-66.98,-67.10,-70.56,-69.91,-72.59,
     &  -71.30,-73.42,-71.86,-73.21,-71.21,-71.86,-69.49,-69.45,-66.30,
     &  -65.54,-60.57,-57.79,
     &   25.36, 14.94,  2.45, -6.61,-17.54,-24.30,-32.90,-38.50,-45.97,
     &  -52.07,-56.65,-58.88,-63.08,-64.34,-67.89,-68.23,-70.96,-70.86,
     &  -73.03,-72.29,-73.92,-72.82,-73.64,-72.12,-72.53,-70.24,-69.88,
     &  -65.54,-62.88,-57.99/
 
c Z = 34
      data b5
     & / 24.54, 10.99,  1.66,-10.03,-17.21,-26.85,-32.73,-41.27,-46.11,
     &  -53.56,-56.30,-61.78,-63.37,-67.89,-68.22,-72.21,-72.17,-75.25,
     &  -74.60,-77.03,-75.92,-77.76,-76.39,-77.59,-75.34,-75.95,-72.43,
     &  -70.54,-66.58,-63.87/
 
c Z = 35 to 37
      data b6
     &  /23.63, 13.79,  1.38, -6.82,-16.81,-23.71,-32.51,-38.50,-46.35,
     &  -50.60,-56.75,-59.2 ,-63.56,-65.31,-69.14,-70.29,-73.23,-73.45,
     &  -76.07,-75.89,-77.97,-77.50,-79.01,-77.78,-78.61,-75.64,-73.85,
     &  -70.73,-68.56,-64.61,-61.55,-56.62,-52.52, 2*0.0,
     &   23.42, 10.05,  1.50, -9.48,-16.54,-26.40,-32.60,-41.45,-46.38,
     &  -54.1 ,-56.89,-62.17,-64.24,-68.98,-70.17,-74.16,-74.44,-77.89,
     &  -77.69,-80.59,-79.98,-82.43,-81.48,-83.26,-80.71,-79.69,-76.72,
     &  -74.96,-71.35,-68.83,-64.09,-61.00,-56.35, 2*0.0,
     &   13.34,  1.94, -6.01,-16.09,-23.47,-32.50,-38.35,-46.47,-51.7 ,
     &  -57.22,-60.48,-64.83,-66.93,-70.79,-72.17,-75.46,-76.19,-79.09,
     &  -79.75,-82.17,-82.75,-84.59,-82.60,-81.70,-79.35,-77.79,-74.81,
     &  -72.70,-68.53,-65.86,-61.23,-58.38,-54.27,-50.89,  0.0 /
 
c Z = 38 to 41
      data b7
     & / 10.34,  2.19, -8.88,-15.98,-26.13,-32.46,-41.65,-47.09,-54.96,
     &  -58.0 ,-63.17,-65.48,-70.30,-71.52,-76.01,-76.80,-80.64,-81.10,
     &  -84.52,-84.88,-87.92,-86.21,-85.94,-83.65,-82.92,-80.16,-78.84,
     &  -75.16,-72.98,-68.80,-66.61,-62.15,-60.22, 2*0.0,
     &   13.50,  2.32, -5.97,-16.13,-23.49,-32.85,-39.28,-47.48,-52.20,
     &  -58.4 ,-61.52,-66.01,-68.18,-72.33,-74.23,-77.85,-79.28,-83.02,
     &  -84.30,-87.70,-86.49,-86.35,-84.83,-84.24,-82.35,-81.24,-78.35,
     &  -76.26,-72.44,-70.20,-67.30,-64.91,-61.89, 2*0.0,
     &   11.22,  2.86, -8.36,-15.92,-26.11,-32.75,-41.84,-47.24,-54.85,
     &  -58.9 ,-64.18,-66.46,-70.64,-73.15,-77.81,-79.35,-83.63,-84.87,
     &  -88.77,-87.89,-88.46,-87.12,-87.27,-85.66,-85.44,-82.95,-81.27,
     &  -77.77,-76.61,-73.46,-71.74,-68.38,-66.67, 2*0.0,
     &   14.28,  3.33, -5.13,-15.56,-23.06,-32.47,-38.81,-46.60,-51.11,
     &  -59.0 ,-59.34,-67.15,-69.83,-74.18,-76.53,-80.58,-82.66,-86.64,
     &  -86.45,-87.21,-86.37,-86.78,-85.61,-85.61,-83.53,-82.33,-79.94,
     &  -78.94,-76.35,-75.32,-72.23,-70.86,-66.88,-64.55,   0.0/
 
c Z = 42 to 45
      data b8
     & /  3.56, -7.04,-14.11,-25.38,-31.81,-39.58,-45.01,-51.90,-55.82,
     &  -62.81,-67.7 ,-72.70,-75.00,-80.17,-82.21,-86.81,-86.81,-88.41,
     &  -87.71,-88.79,-87.54,-88.11,-85.97,-86.19,-83.51,-83.56,-80.85,
     &  -80.33,-77.34,-76.26,-72.93,-71.18,-67.21,-65.28,   0.0,
     &    4.19, -3.90,-14.06,-21.22,-30.60,-36.94,-43.92,-49.06,-56.59,
     &  -60.72,-67.5 ,-70.70,-75.99,-78.94,-83.60,-84.16,-86.02,-85.82,
     &  -87.22,-86.43,-87.32,-86.02,-86.34,-84.57,-84.60,-82.49,-82.29,
     &  -79.78,-79.09,-75.94,-74.61,-70.40,-68.73,-65.02,   0.0,
     &    4.94, -6.11,-12.17,-22.39,-28.89,-38.43,-44.10,-52.80,-57.13,
     &  -64.34,-68.6 ,-74.57,-77.27,-82.56,-83.45,-86.07,-86.11,-88.23,
     &  -87.62,-89.22,-87.95,-89.10,-87.26,-88.09,-85.93,-86.32,-83.91,
     &  -83.65,-80.85,-80.1 ,-76.31,-75.29,-71.70,-70.17,   0.0,
     &    5.10, -2.26,-12.61,-19.84,-29.59,-36.06,-44.97,-50.42,-58.15,
     &  -62.62,-69.26,-72.94,-78.34,-79.62,-82.58,-83.17,-85.51,-85.59,
     &  -87.41,-86.78,-88.02,-86.95,-87.85,-86.36,-86.86,-85.00,-85.01,
     &  -82.9 ,-81.71,-79.24,-78.34,-75.39,-73.98,-70.51,-68.65/
 
c Z = 46 to 47
      data b9
     & / -4.79,-12.20,-22.76,-29.40,-39.15,-45.10,-54.00,-58.67,-66.29,
     &  -70.02,-76.18,-77.80,-81.29,-82.15,-85.23,-85.43,-87.93,-87.48,
     &  -89.39,-88.41,-89.91,-88.37,-89.52,-87.60,-88.35,-86.03,-86.34,
     &  -83.69,-83.49,-80.40,-79.95,-76.37,-75.5 ,-71.64,-70.21,
     &   -2.26,-13.01,-20.47,-30.38,-37.08,-46.15,-51.95,-59.91,-64.55,
     &  -71.16,-72.9 ,-76.76,-78.15,-81.19,-82.00,-84.79,-85.11,-87.07,
     &  -86.94,-88.41,-87.60,-88.72,-87.46,-88.22,-86.62,-87.03,-84.94,
     &  -84.99,-82.56,-82.24,-79.64,-78.56,-75.77,-74.5 ,-70.94/
 
c Z = 48 to 50
      data b10
     & / -5.11,-12.74,-23.48,-30.41,-40.27,-46.55,-55.76,-60.46,-68.14,
     &  -70.19,-74.29,-75.66,-79.42,-80.65,-83.98,-84.33,-87.13,-86.99,
     &  -89.25,-88.51,-90.35,-89.25,-90.58,-89.05,-90.02,-88.09,-88.72,
     &  -86.43,-86.71,-83.90,-83.97,-80.9 ,-80.68,-77.31,-76.71,-73.32,
     &  -72.33,-68.53,-67.27,-65.20,
     &   -3.02,-13.98,-21.70,-31.74,-38.80,-48.05,-53.99,-62.07,-64.90,
     &  -69.51,-70.5 ,-74.60,-76.08,-79.48,-80.61,-83.56,-84.11,-86.49,
     &  -86.47,-88.39,-87.99,-89.37,-88.57,-89.54,-88.25,-88.94,-87.23,
     &  -87.70,-85.73,-85.84,-83.58,-83.43,-80.88,-80.48,-77.81,-76.99,
     &  -74.36,-72.97,-69.99,-68.20,
     &   -6.26,-14.13,-24.98,-32.21,-42.32,-48.51,-57.69,-60.86,-66.29,
     &  -68.28,-71.6 ,-73.23,-77.43,-78.56,-82.01,-82.64,-85.83,-85.94,
     &  -88.66,-88.33,-90.56,-90.03,-91.52,-90.40,-91.65,-90.07,-91.10,
     &  -89.20,-89.94,-87.82,-88.24,-85.90,-86.02,-83.51,-83.34,-80.62,
     &  -80.24,-77.38,-76.62,-71.1 /
 
c Z = 51 to 53
      data b11
     & /-11.86,-19.82,-29.93,-37.01,-46.33,-50.44,-56.31,-59.21,-63.82,
     &  -66.28,-70.30,-72.17,-76.25,-77.18,-80.21,-81.60,-84.42,-84.68,
     &  -87.00,-86.82,-88.64,-88.00,-89.47,-88.42,-89.59,-88.32,-89.22,
     &  -87.62,-88.26,-86.40,-86.71,-84.61,-84.63,-82.39,-82.02,-79.92,
     &  -78.96,-74.03,-69.74,-65.13,
     &   -9.50,-20.62,-27.50,-37.66,-42.09,-49.06,-52.03,-57.76,-60.23,
     &  -65.67,-67.58,-72.28,-73.47,-77.3 ,-78.04,-81.29,-82.35,-85.32,
     &  -85.11,-87.72,-87.18,-89.40,-88.55,-90.30,-89.17,-90.52,-89.03,
     &  -90.07,-88.29,-88.99,-87.01,-87.35,-85.21,-85.21,-82.96,-82.44,
     &  -77.83,-74.42,-69.56,-66.65,
     &   -7.47,-15.30,-25.11,-30.80,-38.06,-42.13,-48.24,-51.74,-57.56,
     &  -59.75,-64.24,-66.67,-71.12,-72.8 ,-75.94,-77.54,-80.33,-81.27,
     &  -83.67,-83.78,-86.28,-86.07,-87.93,-87.36,-88.84,-87.92,-88.99,
     &  -87.74,-88.50,-86.93,-87.44,-85.70,-85.88,-83.95,-83.79,-79.50,
     &  -76.50,-72.30,-68.84,-63.35/
 
c Z = 54 to 56
      data b12
     & / -4.62,-15.26,-20.97,-28.99,-33.81,-40.99,-44.70,-50.64,-53.65,
     &  -59.92,-62.06,-66.39,-68.51,-72.63,-74.0 ,-78.  ,-78.7 ,-81.82,
     &  -82.55,-85.17,-85.25,-87.66,-87.19,-89.17,-88.33,-89.86,-88.70,
     &  -89.88,-88.42,-89.28,-87.65,-88.12,-86.44,-86.42,-82.38,-80.12,
     &  -75.65,-73.00,-68.32,-65.50,
     &   -8.73,-17.05,-22.79,-30.19,-35.33,-41.82,-45.51,-51.66,-54.24,
     &  -59.33,-62.4 ,-66.24,-68.43,-72.34,-73.90,-77.15,-78.12,-81.05,
     &  -81.74,-84.10,-84.35,-86.25,-85.93,-87.50,-86.90,-88.06,-87.16,
     &  -88.08,-86.90,-87.59,-86.34,-86.55,-82.89,-80.71,-77.06,-74.47,
     &  -70.52,-67.71,-63.32,-60.16,
     &   -7.22,-13.12,-21.58,-26.97,-34.19,-38.15,-44.56,-48.19,-54.11,
     &  -57.04,-61.98,-64.  ,-68.9 ,-70.3 ,-74.35,-75.70,-79.07,-79.5 ,
     &  -82.68,-82.8 ,-85.41,-85.07,-87.27,-86.69,-88.44,-87.56,-88.95,
     &  -87.86,-88.89,-87.73,-88.27,-84.92,-83.28,-79.73,-77.83,-73.95,
     &  -71.78,-68.05,-65.04,-61.49/
 
c Z = 57 to 59
      data b13
     &  /-1.28, -9.86,-15.78,-23.73,-28.87,-35.67,-40.14,-46.24,-49.83,
     &  -54.92,-57.77,-62.19,-64.63,-68.62,-70.55,-73.97,-75.40,-78.17,
     &  -78.8 ,-81.35,-81.65,-83.7 ,-83.73,-85.3 ,-85.24,-86.65,-86.02,
     &  -87.13,-86.53,-87.24,-84.33,-82.94,-80.04,-78.19,-74.90,-72.98,
     &  -69.16,-67.24,-63.2 ,-61.21,
     &   -6.36,-14.79,-20.55,-28.59,-33.28,-40.08,-43.87,-49.66,-52.68,
     &  -57.76,-60.26,-64.82,-66.82,-70.80,-72.29,-75.61,-76.51,-79.32,
     &  -79.73,-82.12,-82.24,-84.7 ,-84.63,-86.49,-85.90,-87.57,-86.96,
     &  -88.09,-85.44,-84.54,-81.62,-80.44,-77.10,-75.70,-72.18,-70.4 ,
     &  -66.80,-65.0 ,-62.22,-59.96,
     &   -3.52,-10.08,-18.38,-24.05,-31.28,-35.83,-41.74,-45.39,-50.52,
     &  -53.62,-58.22,-60.81,-64.86,-66.88,-70.28,-71.68,-74.46,-75.41,
     &  -77.82,-78.32,-80.92,-81.37,-83.20,-83.14,-84.83,-84.70,-86.03,
     &  -83.80,-83.08,-80.76,-79.64,-76.74,-75.47,-72.5 ,-70.99,-68.00,
     &  -66.79,-64.57,-62.41,-58.76/
 
c Z = 60 to 62
      data b14
     & / -9.78,-15.87,-24.08,-28.75,-35.30,-39.06,-44.84,-47.99,-53.21,
     &  -55.89,-60.54,-62.67,-66.67,-67.9 ,-71.49,-72.26,-75.23,-75.88,
     &  -79.16,-79.51,-81.27,-82.04,-84.48,-84.20,-85.96,-84.01,-83.76,
     &  -81.44,-80.94,-78.16,-77.42,-74.39,-73.69,-70.96,-70.16,-68.16,
     &  -66.65,-63.07,-61.07,-57.15,
     &   -4.98,-13.36,-18.98,-25.90,-30.31,-36.20,-40.01,-45.30,-48.57,
     &  -53.30,-56.05,-60.12,-62.27,-65.65,-67.06,-69.90,-71.3 ,-73.72,
     &  -74.45,-77.52,-78.39,-80.49,-81.09,-82.97,-81.43,-81.28,-79.46,
     &  -79.05,-76.88,-76.08,-73.61,-73.40,-71.27,-70.67,-68.41,-67.77,
     &  -64.74,-62.89,-59.54,-57.29,
     &  -10.13,-17.96,-22.50,-28.95,-32.90,-38.86,-42.20,-47.61,-50.42,
     &  -55.25,-57.38,-61.40,-62.94,-66.43,-67.71,-70.84,-72.07,-74.62,
     &  -75.94,-78.99,-79.53,-81.98,-80.66,-81.01,-79.28,-79.35,-77.15,
     &  -77.06,-74.59,-74.77,-72.57,-72.47,-70.20,-69.37,-66.8 ,-66.17,
     &  -62.95,-61.27,-57.70,-55.62/
 
c Z = 63 to 64
      data b15
     & / -6.55,-11.88,-18.84,-23.47,-29.60,-33.61,-39.20,-42.69,-47.61,
     &  -50.33,-54.44,-56.52,-60.05,-61.89,-65.19,-66.38,-70.4 ,-71.63,
     &  -74.36,-75.65,-78.00,-77.13,-77.55,-76.24,-76.45,-74.80,-74.66,
     &  -72.90,-73.38,-71.75,-71.83,-70.09,-69.47,-67.21,-66.06,-63.89,
     &  -62.35,-59.27,-57.32,-53.89,
     &  -10.30,-15.11,-21.93,-26.10,-32.42,-36.05,-41.62,-44.45,-49.15,
     &  -51.34,-55.67,-57.42,-61.36,-61.89,-66.13,-68.4 ,-71.36,-72.95,
     &  -76.10,-75.37,-76.28,-75.14,-75.77,-74.20,-74.72,-72.89,-73.72,
     &  -72.08,-72.55,-70.83,-70.70,-68.57,-67.95,-65.52,-64.29,-61.99,
     &  -60.60,-57.27,-55.38,-51.62/
 
c Z = 65 to 67
      data b16
     & / -4.49,-11.71,-16.59,-23.03,-27.34,-33.05,-36.49,-41.23,-44.00,
     &  -48.22,-50.76,-54.67,-55.72,-59.90,-61.92,-66.4 ,-68.0 ,-70.76,
     &  -70.59,-71.50,-71.12,-71.63,-70.73,-71.32,-70.15,-71.26,-70.10,
     &  -70.77,-69.48,-69.54,-67.85,-67.47,-65.68,-64.61,-62.1 ,-61.11,
     &  -58.29,-56.56,-53.29,-51.01,-47.25,-44.65,-40.54,-37.52,-33.08,
     &   -8.24,-15.35,-19.75,-26.13,-29.68,-35.04,-37.96,-42.85,-45.31,
     &  -49.81,-50.88,-55.69,-57.80,-62.9 ,-64.38,-67.91,-67.69,-69.32,
     &  -68.76,-70.13,-69.15,-70.40,-69.16,-70.53,-69.43,-70.42,-69.18,
     &  -69.68,-68.07,-68.19,-66.39,-65.98,-63.62,-62.59,-59.94,-59.34,
     &  -55.6 ,-54.44,-50.81,-48.71,-44.73,-42.29,-37.90,-35.04,-30.31,
     &  -10.15,-16.66,-20.80,-26.27,-29.72,-34.54,-37.76,-42.25,-43.71,
     &  -48.53,-51.18,-55.44,-57.74,-61.67,-62.06,-63.63,-63.65,-65.02,
     &  -64.65,-66.06,-66.49,-66.89,-66.18,-67.34,-66.39,-67.21,-66.05,
     &  -66.39,-64.99,-64.91,-63.08,-62.29,-60.08,-58.81,-56.25,-54.5 ,
     &  -52.08,-50.14,-46.70,-44.30,-40.42,-37.57,-33.32,-30.14,-25.52/
 
 
c Z = 68 to 70
      data b17
     & / -8.60,-12.90,-18.89,-22.50,-27.93,-31.10,-36.13,-37.81,-43.21,
     &  -45.95,-50.75,-53.61,-57.93,-58.15,-60.55,-60.46,-62.62,-62.22,
     &  -64.26,-63.42,-65.70,-64.57,-66.06,-65.20,-66.35,-65.18,-65.95,
     &  -64.53,-64.93,-63.30,-63.00,-60.93,-60.12,-57.73,-56.49,-54.44,
     &  -53.00,-49.71,-47.86,-44.06,-41.69,-37.54,-34.73,-30.46,-27.55,
     &   -9.10,-13.22,-18.68,-22.59,-27.68,-29.76,-35.21,-38.52,-43.39,
     &  -46.06,-50.80,-51.60,-54.00,-54.55,-56.64,-56.89,-58.9 ,-59.30,
     &  -60.7 ,-60.2 ,-62.04,-61.54,-62.74,-61.99,-62.94,-61.89,-62.55,
     &  -61.32,-61.28,-59.80,-59.22,-57.38,-56.26,-53.87,-52.32,-50.46,
     &  -48.66,-45.38,-43.13,-39.45,-36.78,-32.85,-29.97,-25.84,-22.70,
     &   -4.90,-10.95,-14.89,-19.02,-22.88,-28.95,-32.32,-37.81,-40.54,
     &  -45.88,-46.72,-49.63,-49.84,-53.31,-53.41,-56.02,-55.7 ,-58.07,
     &  -58.00,-59.94,-59.4 ,-61.39,-60.18,-61.59,-60.60,-61.58,-60.37,
     &  -60.77,-59.32,-59.26,-57.56,-56.95,-54.70,-53.50,-50.99,-49.70,
     &  -47.68,-45.97,-42.38,-40.30,-36.30,-34.00,-29.99,-27.25,-22.91/
 
c Z = 71 to 73
      data b18
     & / -5.00, -9.45,-13.82,-19.97,-23.84,-29.35,-32.60,-38.00,-39.40,
     &  -42.38,-43.22,-46.09,-47.34,-49.72,-50.26,-52.52,-52.92,-54.77,
     &  -55.04,-56.26,-56.11,-57.47,-57.09,-58.08,-57.31,-57.83,-56.74,
     &  -56.89,-55.58,-55.17,-53.39,-52.40,-50.36,-49.12,-46.69,-46.00,
     &  -42.96,-40.85,-37.66,-35.25,-31.71,-29.15,-25.29,-22.38,-17.98,
     &   -5.77,-12.57,-16.56,-22.69,-26.14,-32.22,-33.59,-37.14,-38.06,
     &  -41.13,-42.25,-46.06,-45.90,-49.18,-49.05,-51.42,-51.59,-53.74,
     &  -53.73,-55.58,-54.81,-56.71,-56.04,-56.39,-56.08,-55.85,-54.49,
     &  -54.58,-52.89,-52.45,-50.47,-49.79,-47.42,-46.06,-43.29,-41.50,
     &  -39.55,-37.94,-34.56,-32.51,-28.65,-26.47,-22.11,-19.93,-15.90,
     &   -2.70, -7.27,-13.50,-17.44,-23.63,-25.49,-29.17,-30.50,-33.76,
     &  -35.33,-38.30,-39.90,-42.54,-43.13,-45.53,-46.15,-48.31,-48.67,
     &  -50.59,-50.68,-52.23,-51.47,-53.10,-52.50,-53.25,-51.47,-51.73,
     &  -50.54,-50.37,-48.94,-48.44,-46.44,-45.30,-42.84,-41.40,-38.62,
     &  -38.29,-35.19,-33.28,-29.92,-27.91,-23.79,-21.68,-18.14,-15.80/
 
c Z = 74 to 76
      data b19
     & / -5.99,-10.13,-16.70,-18.81,-23.02,-24.26,-27.98,-29.48,-32.98,
     &  -34.23,-38.36,-38.36,-41.90,-41.87,-44.47,-44.83,-47.17,-47.24,
     &  -49.21,-48.99,-50.56,-49.97,-51.25,-50.38,-50.44,-49.31,-49.65,
     &  -48.26,-48.25,-46.37,-45.71,-43.39,-42.51,-39.91,-38.67,-35.48,
     &  -34.31,-32.57,-30.93,-26.82,-25.19,-21.74,-19.96,-16.41,-14.34,
     &   -0.65, -7.36, -9.90,-14.23,-16.10,-19.76,-21.80,-25.26,-26.93,
     &  -30.16,-31.97,-34.86,-35.51,-38.10,-38.91,-41.27,-41.83,-43.78,
     &  -44.08,-45.70,-45.65,-46.92,-45.78,-46.60,-45.84,-47.41,-45.45,
     &  -45.81,-44.22,-43.83,-41.93,-41.22,-39.02,-37.99,-35.58,-34.36,
     &  -32.87,-31.30,-28.11,-26.05,-23.02,-21.38,-18.25,-16.22,-12.92,
     &   -2.49, -7.34, -9.10,-13.19,-15.02,-18.96,-20.67,-24.35,-25.79,
     &  -30.12,-30.17,-33.93,-34.14,-36.90,-37.46,-39.95,-40.04,-42.19,
     &  -42.12,-43.91,-43.50,-45.04,-44.29,-44.54,-44.43,-44.26,-42.81,
     &  -43.00,-41.22,-41.14,-38.99,-38.71,-36.40,-35.89,-33.40,-32.44,
     &  -29.70,-28.30,-25.75,-24.59,-21.55,-20.00,-16.92,-14.70, -9.44/
 
c Z = 77 to 79
      data b20
     & / -0.34, -4.38, -6.54,-10.36,-12.42,-16.07,-17.95,-21.24,-23.37,
     &  -26.36,-27.00,-29.67,-30.65,-33.09,-34.05,-35.92,-36.49,-38.44,
     &  -38.45,-40.03,-39.60,-40.88,-39.97,-41.15,-39.17,-40.06,-38.33,
     &  -38.46,-36.71,-36.72,-34.84,-34.54,-32.54,-31.70,-29.46,-28.29,
     &  -26.39,-25.17,-22.56,-21.08,-18.37,-16.19,-11.21, -7.68, -2.06,
     &    1.08, -3.24, -5.12, -9.09,-11.10,-14.79,-16.31,-21.23,-21.16,
     &  -25.32,-25.42,-28.25,-28.93,-31.60,-31.92,-34.15,-34.29,-36.18,
     &  -35.98,-37.60,-37.03,-37.79,-37.22,-37.83,-36.49,-37.33,-35.70,
     &  -36.30,-34.49,-34.79,-32.82,-32.67,-30.44,-29.93,-27.42,-26.63,
     &  -23.75,-23.69,-21.05,-19.27,-14.45,-11.19, -5.58, -1.91,  3.72,
     &    5.76,  3.41, -0.56, -2.62, -6.32, -8.10,-11.52,-14.16,-17.16,
     &  -17.47,-20.36,-21.41,-24.04,-24.85,-27.20,-27.72,-29.80,-30.00,
     &  -31.71,-31.59,-32.96,-32.45,-33.53,-32.89,-33.87,-32.78,-33.42,
     &  -32.29,-32.59,-31.17,-31.16,-29.60,-29.12,-27.28,-26.41,-24.42,
     &  -23.15,-21.50,-19.86,-15.31,-12.19, -7.13, -3.36,  1.91,  5.76/
 
c Z = 80 to 82
      data b21
     & /  6.90,  4.82,  0.65, -1.18, -5.05, -6.63,-11.88,-11.52,-16.32,
     &  -15.91,-18.96,-19.79,-22.60,-23.15,-25.69,-25.97,-28.10,-28.07,
     &  -29.93,-29.55,-31.13,-30.69,-31.82,-31.08,-32.25,-31.07,-31.85,
     &  -30.57,-30.98,-29.57,-29.53,-27.69,-27.37,-25.29,-24.72,-22.31,
     &  -20.97,-16.27,-14.66, -9.63, -6.30, -1.06,  2.38,  7.73, 11.23,
     &   10.09,  7.89,  3.98,  2.04, -1.61, -3.35, -7.97, -8.19,-11.35,
     &  -12.53,-15.42,-16.33,-18.90,-19.61,-21.84,-22.29,-24.25,-24.41,
     &  -26.13,-26.02,-27.44,-27.04,-28.27,-27.73,-28.40,-27.52,-28.14,
     &  -27.07,-27.20,-26.00,-25.78,-24.37,-23.84,-22.28,-21.05,-16.77,
     &  -13.65, -9.26, -6.63, -1.91,  1.54,  6.47, 10.03, 15.21, 18.80,
     &   15.97, 11.55,  9.46,  5.33,  3.55, -0.36, -1.87, -6.82, -6.73,
     &   -9.99,-10.94,-13.94,-14.56,-17.31,-17.73,-20.17,-20.27,-22.55,
     &  -22.38,-24.41,-23.99,-25.72,-25.09,-26.50,-25.27,-26.73,-25.30,
     &  -25.96,-24.81,-25.13,-23.79,-23.81,-22.48,-21.77,-17.64,-14.75,
     &  -10.49, -7.57, -3.11, -0.19,  4.71,  7.72, 12.67, 15.85, 20.94/
 
c Z = 83 to 85
      data b22
     & / 20.60, 16.43, 14.14, 10.46,  8.23,  4.72,  2.98, -0.45, -1.92,
     &   -5.01, -6.17, -8.97, -9.89,-12.41,-13.05,-15.32,-15.68,-17.69,
     &  -17.96,-19.64,-19.54,-20.92,-20.40,-21.47,-20.81,-21.58,-20.73,
     &  -21.08,-20.05,-20.08,-18.89,-18.28,-14.81,-11.87, -8.14, -5.24,
     &   -1.22,  1.71,  6.10,  8.99, 13.44, 16.65, 21.32, 24.55, 29.41,
     &   22.29, 17.67, 15.22, 11.00,  8.95,  5.28,  3.75,  0.68, -0.32,
     &   -3.17, -3.91, -6.97, -7.58, -9.91,-10.35,-12.64,-12.87,-14.94,
     &  -14.80,-16.73,-16.33,-17.92,-17.35,-18.58,-17.55,-18.20,-17.17,
     &  -17.49,-16.39,-15.98,-12.46,-10.39, -6.68, -4.49, -0.54,  1.76,
     &    6.00,  8.35, 13.11, 15.89, 20.43, 23.22, 27.59, 30.43, 34.87,
     &   24.87, 20.43, 17.94, 14.21, 12.32,  9.08,  7.70,  5.10,  3.66,
     &    0.83, -0.23, -2.66, -3.28, -5.43, -5.84, -7.84, -8.94,-10.74,
     &  -10.76,-12.29,-11.90,-13.03,-12.49,-13.28,-12.56,-12.90,-11.99,
     &  -11.67, -8.64, -6.60, -3.40, -1.27,  2.23,  4.38,  8.08, 10.53,
     &   14.36, 16.93, 21.09, 23.85, 27.73, 30.51, 34.53, 37.53, 41.76/
 
c Z = 86 to 88
      data b23
     & / 28.33, 25.76, 21.59, 19.60, 15.96, 14.62, 11.43, 10.32,  7.29,
     &    6.17,  3.27,  2.52,  0.03, -0.46, -2.67, -2.93, -5.27, -5.39,
     &   -7.40, -7.33, -8.77, -8.67, -9.70, -8.97, -9.62, -8.78, -8.68,
     &   -5.72, -4.34, -1.19,  0.23,  3.63,  5.20,  8.83, 10.59, 14.49,
     &   16.37, 20.80, 23.07, 26.85, 29.17, 33.09, 35.63, 39.76, 42.45,
     &   31.41, 28.98, 25.26, 23.28, 20.05, 18.76, 15.65, 14.67, 11.72,
     &   10.50,  7.90,  7.05,  4.74,  4.12,  2.02,  0.65, -1.26, -1.41,
     &   -2.96, -2.70, -3.83, -3.40, -4.20, -3.60, -3.57, -0.98,  0.29,
     &    2.96,  4.29,  7.04,  8.61, 11.46, 13.27, 16.38, 18.38, 21.63,
     &   23.85, 27.30, 29.59, 33.63, 36.08, 39.80, 42.44, 46.33, 49.28,
     &   33.18, 31.31, 27.65, 26.18, 22.75, 21.35, 18.60, 17.33, 14.39,
     &   13.44, 10.78, 10.11,  7.57,  7.09,  4.55,  4.25,  2.29,  2.41,
     &    0.89,  0.80, -0.23,  0.31,  0.08,  2.51,  3.27,  5.86,  6.63,
     &    9.36, 10.25, 12.94, 14.30, 17.23, 18.80, 21.99, 23.66, 27.17,
     &   28.94, 32.66, 34.96, 38.62, 40.86, 44.73, 47.30, 51.29, 54.08/
 
c Z = 89 to 91
      data b24
     & / 37.36, 35.32, 31.87, 30.20, 27.40, 26.30, 23.29, 21.94, 19.20,
     &   18.10, 15.53, 14.66, 12.31, 11.70,  8.89,  8.62,  7.09,  7.24,
     &    6.10,  6.38,  5.97,  8.06,  8.69, 10.82, 11.54, 13.73, 14.50,
     &   16.60, 17.82, 20.20, 21.63, 24.30, 25.85, 28.89, 30.90, 34.27,
     &   35.91, 39.49, 41.69, 45.19, 47.74, 51.36, 54.14, 58.06, 61.10,
     &   39.79, 38.01, 34.58, 33.54, 30.56, 29.18, 26.06, 24.91, 21.92,
     &   20.96, 18.13, 17.46, 15.03, 14.77, 12.94, 12.85, 11.13, 10.89,
     &   10.25, 12.16, 12.35, 14.45, 14.65, 16.92, 17.18, 19.36, 19.98,
     &   22.28, 23.18, 25.80, 26.75, 29.58, 30.86, 33.81, 35.44, 38.73,
     &   40.61, 44.25, 46.67, 50.24, 52.64, 56.56, 59.22, 63.25, 66.23,
     &   47.70, 44.25, 42.83, 39.68, 38.03, 35.05, 33.49, 30.69, 29.30,
     &   26.47, 25.42, 22.98, 22.45, 20.59, 20.16, 17.68, 17.68, 17.02,
     &   18.60, 18.48, 20.32, 22.18, 21.94, 22.31, 23.78, 24.32, 26.01,
     &   26.83, 28.86, 29.89, 32.17, 33.42, 35.92, 37.49, 40.33, 42.33,
     &   45.34, 47.64, 50.77, 53.56, 57.12, 59.76, 63.61, 66.34, 70.64/
 
c Z = 92 to 94
      data b25
     & / 50.98, 47.38, 45.74, 42.44, 40.82, 37.66, 36.44, 33.17, 32.08,
     &   29.24, 28.60, 26.44, 25.97, 23.86, 23.55, 22.15, 23.83, 23.60,
     &   25.03, 25.10, 26.07, 26.15, 27.37, 27.17, 29.81, 29.21, 31.18,
     &   31.60, 33.78, 34.59, 36.92, 38.14, 40.92, 42.44, 45.39, 47.31,
     &   50.57, 52.71, 56.53, 58.78, 62.59, 64.98, 69.23, 72.39, 77.12,
     &   55.25, 51.91, 50.16, 46.94, 45.36, 42.30, 40.81, 37.98, 36.93,
     &   34.80, 33.99, 31.95, 31.37, 29.83, 31.14, 30.81, 32.74, 31.02,
     &   32.25, 32.33, 33.49, 33.28, 34.45, 33.75, 35.22, 35.62, 38.15,
     &   38.83, 39.95, 41.04, 43.37, 44.87, 47.45, 49.31, 52.32, 54.26,
     &   57.41, 59.92, 63.61, 65.94, 69.85, 72.97, 77.18, 80.46, 84.73,
     &   58.13, 54.60, 52.91, 49.62, 48.19, 44.96, 43.89, 41.36, 40.51,
     &   38.06, 37.37, 35.53, 36.75, 36.05, 36.67, 35.95, 36.97, 36.55,
     &   37.73, 37.07, 38.20, 37.83, 39.25, 38.35, 40.02, 40.34, 42.16,
     &   42.88, 45.09, 46.16, 48.58, 50.12, 52.95, 54.71, 57.75, 59.80,
     &   63.17, 65.39, 69.31, 72.08, 76.28, 79.14, 83.37, 86.24, 90.58/
 
c Z = 95 to 97
      data b26
     & / 62.57, 59.17, 57.72, 54.48, 53.08, 50.47, 49.28, 46.87, 45.92,
     &   43.89, 44.87, 43.95, 45.39, 43.40, 44.07, 43.30, 43.78, 43.03,
     &   43.77, 43.35, 44.38, 44.27, 45.53, 45.69, 47.26, 47.75, 48.42,
     &   49.38, 51.50, 52.93, 55.47, 57.17, 59.88, 61.89, 64.99, 67.32,
     &   70.81, 73.57, 77.43, 80.27, 84.14, 86.98, 90.94, 2*0.0,
     &   65.57, 62.22, 60.77, 57.72, 56.50, 53.67, 52.60, 50.26, 51.07,
     &   49.58, 51.27, 48.89, 49.64, 48.41, 48.88, 47.73, 48.45, 47.60,
     &   48.59, 48.08, 49.29, 49.03, 50.58, 49.38, 52.51, 51.70, 53.70,
     &   54.80, 57.18, 58.45, 61.00, 62.61, 65.53, 67.39, 70.75, 72.98,
     &   76.64, 79.68, 83.69, 86.15, 90.09, 92.67, 3*0.0,
     &   70.63, 67.57, 65.96, 63.15, 61.82, 59.32, 59.84, 58.01, 59.35,
     &   57.18, 57.01, 55.55, 55.73, 54.53, 54.87, 54.00, 54.58, 54.00,
     &   54.82, 54.49, 55.67, 55.64, 56.96, 57.39, 59.12, 58.68, 60.70,
     &   61.81, 64.88, 65.48, 68.11, 69.84, 72.95, 75.22, 79.25, 81.56,
     &   85.25, 87.68, 91.27, 93.82, 5*0.0/
 
c Z = 98 to 100
      data b27
     & / 73.87, 70.63, 69.25, 66.40, 66.76, 64.69, 65.91, 63.41, 63.26,
     &   61.47, 61.56, 60.01, 60.39, 59.11, 59.66, 58.68, 59.46, 58.74,
     &   59.88, 59.43, 60.68, 59.32, 62.17, 61.46, 63.38, 64.09, 66.13,
     &   67.24, 69.72, 71.17, 74.13, 76.03, 79.29, 81.34, 85.57, 87.80,
     &   91.35, 93.53, 97.25, 99.64, 5*0.0,
     &   79.07, 76.21, 76.16, 73.85, 74.77, 71.81, 71.26, 69.45, 69.17,
     &   67.65, 67.59, 66.28, 66.47, 65.41, 65.82, 65.00, 65.73, 65.20,
     &   66.19, 66.14, 67.20, 67.45, 68.97, 68.55, 70.29, 71.11, 73.95,
     &   74.51, 77.29, 79.01, 81.99, 84.08, 87.99, 90.07, 93.39, 95.53,
     &   98.85,101.22,105.01,107.96, 5*0.0,
     &   83.92, 83.80, 81.33, 82.02, 78.87, 78.27, 76.07, 75.73, 73.83,
     &   73.85, 72.17, 72.36, 70.89, 71.28, 70.15, 70.84, 69.92, 70.86,
     &   70.40, 71.39, 70.12, 72.70, 71.89, 74.53, 74.07, 75.98, 76.81,
     &   79.34, 80.90, 83.79, 85.48, 88.58, 91.08, 94.19, 96.05, 99.35,
     &  101.34,105.08,107.68,111.95, 5*0.0/
 
c ======================================================================
 
c   Energy shifts from calculated pairing gaps to be used for level
c   density:
 
c Z = 1 to 7 [Very crude approximation; 12/sqrt(A)]
      data epaira /
     &   0.00,  0.00,  6.00,  0.00,  4.90,  0.00, 9*0.0,
     &   6.93, 12.00,  5.37,  9.80,  4.54,  8.49,  4.00,  7.59, 7*0.0,
     &   0.00,  5.37,  0.00,  4.54,  0.00,  4.00,  0.00,  3.62, 7*0.0,
     &   9.80,  4.54,  8.49,  4.00,  7.59,  3.62,  6.93,  3.33,  6.41,
     &  6*0.0,
     &   4.54,  0.00,  4.00,  0.00,  3.62,  0.00,  3.33,  0.00,  3.10,
     &   0.00,  2.91, 4*0.0,
     &   8.49,  4.00,  7.59,  3.62,  6.93,  3.33,  6.41,  3.10,  6.00,
     &   2.91,  5.66,  2.75,  5.37, 2*0.0,
     &   3.62,  0.00,  3.33,  0.00,  3.10,  0.00,  2.91,  0.00,  2.75,
     &   0.00,  2.62,  0.00, 3*0.0/
 
c Z = 8 to 17
      data c1
     &    /8.00, 3.80, 7.65, 3.65, 7.18, 3.54, 6.47, 3.41, 5.96, 3.41,
     &     5.59, 3.41, 5.24, 17*0.0,
     &     0.00, 4.00, 0.00, 3.55, 0.00, 3.05, 0.00, 2.53, 0.00, 2.30,
     &     0.00, 1.92, 0.00, 17*0.0,
     &     2.80, 6.10, 2.80, 6.14, 2.84, 6.20, 2.91, 5.61, 2.88, 5.12,
     &     2.75, 4.64, 2.72, 4.56, 2.72, 15*0.0,
     &     3.50, 0.00, 3.30, 0.00, 3.04, 0.00, 2.75, 0.00, 2.43, 0.00,
     &     2.22, 0.00, 2.00, 0.00, 1.88, 15*0.0,
     &     2.30, 5.67, 2.40, 5.41, 2.52, 5.24, 2.51, 5.16, 2.47, 4.87,
     &     2.40, 4.58, 2.37, 4.29, 2.36, 15*0.0,
     &     3.30, 0.00, 2.81, 0.00, 2.63, 0.00, 2.74, 0.00, 2.16, 0.00,
     &     2.16, 0.00, 1.97, 0.00, 1.87, 0.00, 1.66, 13*0.0,
     &     5.42, 2.13, 5.00, 2.31, 4.95, 2.49, 5.25, 2.38, 4.48, 2.31,
     &     4.37, 2.26, 4.36, 2.24, 4.16, 2.23, 3.92, 13*0.0,
     &     3.41, 0.00, 2.77, 0.00, 2.44, 0.00, 2.49, 0.00, 2.14, 0.00,
     &     2.06, 0.00, 2.11, 0.00, 1.93, 0.00, 1.83, 0.00, 1.66,11*0.0,
     &     5.08, 1.71, 4.55, 2.01, 4.64, 2.11, 4.27, 1.96, 4.04, 2.11,
     &     4.11, 2.05, 4.14, 2.03, 3.96, 2.08, 3.94, 2.14, 3.81, 2.09,
     &    10*0.0,
     &     3.43, 0.00, 2.77, 0.00, 2.40, 0.00, 2.28, 0.00, 2.07, 0.00,
     &     2.17, 0.00, 2.09, 0.00, 1.95, 0.00, 1.84, 0.00, 1.71, 0.00,
     &     1.58, 9*0.0/
 
c Z = 18 to 23
      data c2
     &    /1.67, 4.46, 1.78, 4.32, 1.91, 4.31, 1.79, 3.84, 2.00, 4.01,
     &     1.99, 4.06, 1.97, 3.91, 1.96, 3.93, 1.97, 3.72, 2.05, 3.63,
     &     2.06, 9*0.0,
     &     2.77, 0.00, 2.38, 0.00, 2.27, 0.00, 2.04, 0.00, 1.99, 0.00,
     &     2.11, 0.00, 2.00, 0.00, 1.94, 0.00, 1.72, 0.00, 1.63, 0.00,
     &     1.47, 0.00, 8*0.0,
     &     4.39, 1.69, 4.14, 1.82, 4.14, 1.89, 3.94, 1.93, 3.96, 1.98,
     &     4.08, 1.96, 3.96, 2.04, 3.97, 2.05, 3.75, 2.12, 3.75, 2.14,
     &     3.71, 2.16, 3.74, 7*0.0,
     &     0.00, 2.37, 0.00, 2.26, 0.00, 2.03, 0.00, 1.98, 0.00, 2.11,
     &     0.00, 1.99, 0.00, 1.98, 0.00, 1.81, 0.00, 1.74, 0.00, 1.55,
     &     0.00, 1.56, 0.00, 7*0.0,
     &     4.24, 1.61, 4.03, 1.84, 3.92, 1.88, 3.87, 1.89, 4.00, 1.87,
     &     3.88, 1.84, 3.81, 1.83, 3.62, 1.96, 3.70, 1.98, 3.52, 2.04,
     &     3.61, 2.01, 3.57, 7*0.0,
     &     0.00, 2.27, 0.00, 2.10, 0.00, 1.98, 0.00, 2.12, 0.00, 2.04,
     &     0.00, 1.97, 0.00, 1.77, 0.00, 1.72, 0.00, 1.51, 0.00, 1.57,
     &     0.00, 1.58, 0.00, 1.58 ,0.00, 1.52, 4*0.0/
 
c Z = 24 to 29
      data c3
     &    /3.87, 1.59, 3.76, 1.68, 3.80, 1.82, 3.96, 1.81, 3.85, 1.78,
     &     3.74, 1.65, 3.53, 1.75, 3.56, 1.85, 3.45, 1.86, 3.43, 1.87,
     &     3.44, 1.89, 3.48, 1.85, 3.31, 1.85, 4*0.0,
     &     0.00, 2.16, 0.00, 2.01, 0.00, 2.15, 0.00, 2.06, 0.00, 1.91,
     &     0.00, 1.80, 0.00, 1.71, 0.00, 1.61, 0.00, 1.56, 0.00, 1.56,
     &     0.00, 1.57, 0.00, 1.49, 0.00, 1.49, 4*0.0,
     &     3.71, 1.58, 3.63, 1.64, 3.80, 1.64, 3.71, 1.61, 3.59, 1.55,
     &     3.32, 1.53, 3.28, 1.59, 3.08, 1.68, 3.25, 1.74, 3.30, 1.75,
     &     3.31, 1.71, 3.24, 1.69, 3.20, 1.70, 4*0.0,
     &     0.00, 2.07, 0.00, 2.17, 0.00, 2.10, 0.00, 1.97, 0.00, 1.77,
     &     0.00, 1.72, 0.00, 1.51, 0.00, 1.50, 0.00, 1.52, 0.00, 1.53,
     &     0.00, 1.53, 0.00, 1.51, 0.00, 1.51, 4*0.0,
     &     3.45, 1.47, 3.66, 1.51, 3.63, 1.53, 3.53, 1.54, 3.29, 1.54,
     &     3.27, 1.57, 3.07, 1.59, 3.10, 1.60, 3.17, 1.66, 3.18, 1.66,
     &     3.19, 1.67, 3.18, 1.66, 3.16, 1.65, 3.09 , 3*0.0,
     &     0.00, 2.18, 0.00, 2.12, 0.00, 1.98, 0.00, 1.78, 0.00, 1.72,
     &     0.00, 1.52, 0.00, 1.50, 0.00, 1.50, 0.00, 1.50, 0.00, 1.51,
     &     0.00, 1.51, 0.00, 1.50, 0.00, 1.45, 4*0.0/
 
c Z = 30 to 34
      data c4
     &    /1.45, 3.70, 1.54, 3.49, 1.49, 3.35, 1.49, 3.24, 1.48, 3.02,
     &     1.44, 2.94, 1.48, 3.02, 1.49, 3.01, 1.49, 2.98, 1.49, 2.96,
     &     1.46, 2.96, 1.46, 2.90, 1.51, 2.80, 1.50, 2.60, 1.42, 2.68,
     &     2.14, 0.00, 1.95, 0.00, 1.82, 0.00, 1.72, 0.00, 1.60, 0.00,
     &     1.53, 0.00, 1.52, 0.00, 1.53, 0.00, 1.51, 0.00, 1.50, 0.00,
     &     1.50, 0.00, 1.43, 0.00, 1.32, 0.00, 1.22, 0.00, 1.29, 0.00,
     &     1.54, 3.50, 1.55, 3.36, 1.52, 3.28, 1.50, 3.10, 1.48, 2.98,
     &     1.46, 2.98, 1.43, 2.99, 1.47, 2.99, 1.49, 2.99, 1.49, 2.97,
     &     1.47, 2.92, 1.46, 2.81, 1.45, 2.70, 1.43, 2.78, 1.42, 2.58,
     &     1.93, 0.00, 1.83, 0.00, 1.74, 0.00, 1.59, 0.00, 1.53, 0.00,
     &     1.55, 0.00, 1.54, 0.00, 1.53, 0.00, 1.51, 0.00, 1.46, 0.00,
     &     1.47, 0.00, 1.37, 0.00, 1.28, 0.00, 1.37, 0.00, 1.19, 0.00,
     &     1.54, 3.38, 1.55, 3.26, 1.54, 3.12, 1.52, 3.06, 1.50, 3.06,
     &     1.41, 2.96, 1.42, 2.99, 1.54, 3.00, 1.50, 2.97, 1.50, 2.97,
     &     1.47, 2.84, 1.47, 2.75, 1.46, 2.86, 1.44, 2.64, 1.45, 0.00/
 
c Z = 35 to 38
      data c5
     &    /1.82, 0.00, 1.76, 0.00, 1.61, 0.00, 1.57, 0.00, 1.57, 0.00,
     &     1.56, 0.00, 1.54, 0.00, 1.52, 0.00, 1.47, 0.00, 1.56, 0.00,
     &     1.43, 0.00, 1.32, 0.00, 1.40, 0.00, 1.21, 0.00, 1.19, 0.00,
     &     1.20, 0.00, 1.22, 2*0.0,
     &     1.52, 3.26, 1.49, 3.07, 1.45, 3.02, 1.44, 3.00, 1.42, 2.97,
     &     1.37, 2.87, 1.40, 2.81, 1.49, 2.99, 1.48, 3.13, 1.52, 3.01,
     &     1.50, 2.85, 1.48, 2.88, 1.47, 2.67, 1.48, 2.67, 1.48, 2.70,
     &     1.46, 2.59, 1.42, 2*0.0,
     &     0.00, 1.72, 0.00, 1.66, 0.00, 1.57, 0.00, 1.54, 0.00, 1.48,
     &     0.00, 1.41, 0.00, 1.51, 0.00, 1.61, 0.00, 1.50, 0.00, 1.34,
     &     0.00, 1.40, 0.00, 1.20, 0.00, 1.21, 0.00, 1.23, 0.00, 1.17,
     &     0.00, 1.12, 0.00, 1.08, 0.00,
     &     3.14, 1.42, 3.11, 1.43, 3.03, 1.38, 2.91, 1.35, 2.77, 1.26,
     &     2.67, 1.26, 3.16, 1.53, 3.15, 1.53, 3.04, 1.53, 2.85, 1.51,
     &     2.92, 1.50, 2.69, 1.50, 2.75, 1.49, 2.66, 1.40, 2.58, 1.40,
     &     2.52, 1.40, 2.49, 2*0.0/
 
c Z = 39 to 42
      data c6
     &    /0.00, 1.76, 0.00, 1.69, 0.00, 1.58, 0.00, 1.50, 0.00, 1.42,
     &     0.00, 1.36, 0.00, 1.63, 0.00, 1.52, 0.00, 1.34, 0.00, 1.41,
     &     0.00, 1.19, 0.00, 1.21, 0.00, 1.26, 0.00, 1.23, 0.00, 1.16,
     &     0.00, 1.11, 0.00, 2*0.0,
     &     3.10, 1.38, 3.03, 1.31, 2.90, 1.31, 2.80, 1.29, 2.70, 1.28,
     &     3.15, 1.50, 3.14, 1.50, 3.04, 1.51, 2.85, 1.54, 2.94, 1.52,
     &     2.72, 1.54, 2.77, 1.50, 2.68, 1.38, 2.57, 1.33, 2.51, 1.33,
     &     2.46, 1.34, 2.43, 2*0.0,
     &     0.00, 1.72, 0.00, 1.60, 0.00, 1.51, 0.00, 1.43, 0.00, 1.53,
     &     0.00, 1.65, 0.00, 1.53, 0.00, 1.35, 0.00, 1.41, 0.00, 1.18,
     &     0.00, 1.21, 0.00, 1.22, 0.00, 1.22, 0.00, 1.17, 0.00, 1.13,
     &     0.00, 1.11, 0.00, 1.08, 0.00,
     &     1.36, 2.90, 1.35, 2.82, 1.30, 2.90, 1.34, 3.11, 1.46, 3.12,
     &     1.46, 3.01, 1.47, 2.82, 1.48, 2.89, 1.50, 2.72, 1.54, 2.75,
     &     1.47, 2.66, 1.42, 2.59, 1.33, 2.50, 1.31, 2.43, 1.30, 2.42,
     &     1.30, 2.39, 1.30, 2.33,  0.0/
 
c Z = 43 to 46
      data c7
     &    /1.56, 0.00, 1.57, 0.00, 1.57, 0.00, 1.67, 0.00, 1.68, 0.00,
     &     1.56, 0.00, 1.36, 0.00, 1.43, 0.00, 1.19, 0.00, 1.19, 0.00,
     &     1.21, 0.00, 1.22, 0.00, 1.18, 0.00, 1.15, 0.00, 1.13, 0.00,
     &     1.17, 0.00, 1.12, 0.00,  0.0,
     &     1.28, 2.82, 1.40, 2.99, 1.44, 3.12, 1.44, 3.13, 1.44, 3.01,
     &     1.45, 2.83, 1.45, 2.91, 1.41, 2.60, 1.43, 2.59, 1.38, 2.58,
     &     1.37, 2.59, 1.33, 2.53, 1.32, 2.48, 1.31, 2.46, 1.30, 2.54,
     &     1.36, 2.48, 1.36, 2.45,  0.0,
     &     1.57, 0.00, 1.60, 0.00, 1.68, 0.00, 1.70, 0.00, 1.58, 0.00,
     &     1.39, 0.00, 1.45, 0.00, 1.19, 0.00, 1.19, 0.00, 1.20, 0.00,
     &     1.22, 0.00, 1.24, 0.00, 1.23, 0.00, 1.21, 0.00, 1.17, 0.00,
     &     1.12, 0.00, 1.10, 0.00, 1.09,
     &     2.87, 1.27, 2.96, 1.27, 2.97, 1.26, 2.85, 1.26, 2.66, 1.26,
     &     2.74, 1.26, 2.47, 1.26, 2.44, 1.22, 2.45, 1.25, 2.47, 1.28,
     &     2.53, 1.28, 2.54, 1.27, 2.51, 1.26, 2.46, 1.27, 2.41, 1.27,
     &     2.38, 1.27, 2.38, 1.28, 2.40/
 
c Z = 47
      data c8
     &    /0.00, 1.69, 0.00, 1.71, 0.00, 1.60, 0.00, 1.39, 0.00, 1.44,
     &     0.00, 1.19, 0.00, 1.20, 0.00, 1.20, 0.00, 1.21, 0.00, 1.23,
     &     0.00, 1.25, 0.00, 1.25, 0.00, 1.18, 0.00, 1.15, 0.00, 1.12,
     &     0.00, 1.12, 0.00, 1.14, 0.00/
 
c Z = 48 to 51
      data c9
     &    /2.70, 1.01, 2.72, 1.00, 2.60, 1.00, 2.40, 0.99, 2.46, 1.00,
     &     2.19, 0.99, 2.20, 0.99, 2.24, 1.05, 2.29, 1.13, 2.40, 1.20,
     &     2.46, 1.24, 2.50, 1.24, 2.50, 1.22, 2.38, 1.23, 2.36, 1.23,
     &     2.48, 1.27, 2.48, 1.27, 2.56, 1.30, 2.42, 1.31, 2.20, 1.31,
     &     0.00, 1.72, 0.00, 1.61, 0.00, 1.40, 0.00, 1.46, 0.00, 1.19,
     &     0.00, 1.21, 0.00, 1.22, 0.00, 1.22, 0.00, 1.22, 0.00, 1.24,
     &     0.00, 1.26, 0.00, 1.29, 0.00, 1.32, 0.00, 1.33, 0.00, 1.32,
     &     0.00, 1.30, 0.00, 1.25, 0.00, 1.14, 0.00, 0.99, 0.00, 1.15,
     &     3.09, 1.35, 2.98, 1.35, 2.79, 1.35, 2.86, 1.36, 2.60, 1.39,
     &     2.61, 1.39, 2.62, 1.31, 2.55, 1.28, 2.55, 1.29, 2.59, 1.35,
     &     2.63, 1.40, 2.68, 1.38, 2.73, 1.40, 2.76, 1.39, 2.80, 1.39,
     &     2.78, 1.38, 2.70, 1.37, 2.56, 1.37, 2.39, 1.37, 2.54, 1.38,
     &     1.61, 0.00, 1.41, 0.00, 1.46, 0.00, 1.20, 0.00, 1.19, 0.00,
     &     1.20, 0.00, 1.21, 0.00, 1.22, 0.00, 1.24, 0.00, 1.25, 0.00,
     &     1.26, 0.00, 1.27, 0.00, 1.29, 0.00, 1.30, 0.00, 1.29, 0.00,
     &     1.28, 0.00, 1.20, 0.00, 1.03, 0.00, 1.18, 0.00, 0.95, 0.00/
 
c Z = 52 to 55
      data c10
     &    /1.20, 2.61, 1.18, 2.65, 1.17, 2.37, 1.16, 2.30, 1.10, 2.30,
     &     1.10, 2.33, 1.12, 2.35, 1.09, 2.36, 1.08, 2.35, 1.10, 2.36,
     &     1.12, 2.35, 1.12, 2.38, 1.12, 2.43, 1.13, 2.40, 1.13, 2.49,
     &     1.15, 2.39, 1.15, 2.22, 1.14, 2.34, 1.14, 2.09, 1.15, 2.17,
     &     1.42, 0.00, 1.45, 0.00, 1.22, 0.00, 1.22, 0.00, 1.22, 0.00,
     &     1.22, 0.00, 1.23, 0.00, 1.22, 0.00, 1.20, 0.00, 1.17, 0.00,
     &     1.22, 0.00, 1.24, 0.00, 1.24, 0.00, 1.23, 0.00, 1.19, 0.00,
     &     1.14, 0.00, 1.07, 0.00, 1.20, 0.00, 0.96, 0.00, 0.94, 0.00,
     &     1.27, 2.72, 1.21, 2.43, 1.16, 2.36, 1.13, 2.33, 1.10, 2.34,
     &     1.13, 2.35, 1.14, 2.34, 1.16, 2.32, 1.17, 2.30, 1.18, 2.27,
     &     1.19, 2.28, 1.16, 2.29, 1.16, 2.32, 1.16, 2.36, 1.17, 2.46,
     &     1.22, 2.30, 1.21, 2.42, 1.20, 2.15, 1.16, 2.09, 1.16, 2.06,
     &     0.00, 1.27, 0.00, 1.23, 0.00, 1.22, 0.00, 1.24, 0.00, 1.20,
     &     0.00, 1.17, 0.00, 1.16, 0.00, 1.14, 0.00, 1.10, 0.00, 1.10,
     &     0.00, 1.10, 0.00, 1.11, 0.00, 1.10, 0.00, 1.11, 0.00, 1.08,
     &     0.00, 1.21, 0.00, 0.96, 0.00, 0.92, 0.00, 0.91, 0.00, 0.89/
 
c Z = 56 to 59
      data c11
     &    /2.49, 1.15, 2.36, 1.08, 2.37, 1.11, 2.37, 1.13, 2.35, 1.16,
     &     2.31, 1.16, 2.30, 1.17, 2.28, 1.17, 2.25, 1.18, 2.22, 1.19,
     &     2.23, 1.17, 2.27, 1.18, 2.28, 1.19, 2.31, 1.22, 2.35, 1.25,
     &     2.47, 1.25, 2.23, 1.18, 2.06, 1.14, 2.06, 1.15, 2.06, 1.19,
     &     0.00, 1.27, 0.00, 1.26, 0.00, 1.24, 0.00, 1.20, 0.00, 1.15,
     &     0.00, 1.13, 0.00, 1.11, 0.00, 1.07, 0.00, 1.05, 0.00, 1.05,
     &     0.00, 1.06, 0.00, 1.05, 0.00, 1.10, 0.00, 1.08, 0.00, 1.21,
     &     0.00, 0.97, 0.00, 0.92, 0.00, 0.92, 0.00, 0.90, 0.00, 0.87,
     &     1.09, 2.36, 1.05, 2.31, 1.06, 2.28, 1.06, 2.23, 1.08, 2.22,
     &     1.09, 2.21, 1.12, 2.21, 1.16, 2.23, 1.18, 2.20, 1.18, 2.23,
     &     1.19, 2.27, 1.21, 2.33, 1.26, 2.39, 1.28, 2.51, 1.28, 2.25,
     &     1.20, 2.08, 1.14, 2.06, 1.14, 2.07, 1.15, 2.04, 1.15, 2.02,
     &     1.31, 0.00, 1.26, 0.00, 1.22, 0.00, 1.17, 0.00, 1.14, 0.00,
     &     1.11, 0.00, 1.08, 0.00, 1.05, 0.00, 1.03, 0.00, 1.02, 0.00,
     &     1.06, 0.00, 1.08, 0.00, 1.09, 0.00, 1.22, 0.00, 0.97, 0.00,
     &     0.93, 0.00, 0.95, 0.00, 0.95, 0.00, 0.91, 0.00, 0.86, 0.00/
 
c Z = 60 to 63
      data c12
     &    /2.28, 1.01, 2.23, 0.99, 2.17, 1.00, 2.15, 1.02, 2.14, 1.04,
     &     2.13, 1.07, 2.14, 1.10, 2.13, 1.18, 2.22, 1.20, 2.26, 1.23,
     &     2.35, 1.30, 2.40, 1.30, 2.53, 1.30, 2.28, 1.22, 2.12, 1.15,
     &     2.11, 1.13, 2.08, 1.12, 2.04, 1.12, 1.99, 1.13, 1.99, 1.15,
     &     0.00, 1.24, 0.00, 1.18, 0.00, 1.14, 0.00, 1.12, 0.00, 1.09,
     &     0.00, 1.06, 0.00, 1.04, 0.00, 1.04, 0.00, 1.02, 0.00, 1.05,
     &     0.00, 1.09, 0.00, 1.22, 0.00, 0.98, 0.00, 0.95, 0.00, 0.97,
     &     0.00, 0.98, 0.00, 0.93, 0.00, 0.88, 0.00, 0.87, 0.00, 0.86,
     &     1.00, 2.19, 1.01, 2.14, 1.00, 2.12, 1.02, 2.12, 1.05, 2.11,
     &     1.05, 2.09, 1.06, 2.16, 1.19, 2.23, 1.23, 2.36, 1.31, 2.41,
     &     1.32, 2.54, 1.32, 2.29, 1.22, 2.13, 1.16, 2.10, 1.10, 2.06,
     &     1.07, 2.01, 1.07, 1.98, 1.08, 1.96, 1.09, 1.96, 1.11, 1.92,
     &     1.19, 0.00, 1.15, 0.00, 1.12, 0.00, 1.09, 0.00, 1.06, 0.00,
     &     1.04, 0.00, 1.03, 0.00, 1.01, 0.00, 1.04, 0.00, 1.09, 0.00,
     &     1.22, 0.00, 0.98, 0.00, 0.96, 0.00, 0.99, 0.00, 0.99, 0.00,
     &     0.96, 0.00, 0.91, 0.00, 0.89, 0.00, 0.87, 0.00, 0.83, 0.00/
 
c Z = 64
      data c13
     &    /2.17, 1.00, 2.15, 1.01, 2.10, 1.03, 2.11, 1.03, 2.07, 1.05,
     &     2.09, 1.14, 2.17, 1.30, 2.35, 1.31, 2.45, 1.35, 2.57, 1.34,
     &     2.32, 1.24, 2.13, 1.14, 2.09, 1.08, 2.04, 1.05, 2.01, 1.05,
     &     1.97, 1.05, 1.95, 1.06, 1.94, 1.07, 1.92, 1.07, 1.91, 1.09/
 
c Z = 65 to 67
      data c14
     &    /0.00, 1.14, 0.00, 1.10, 0.00, 1.08, 0.00, 1.04, 0.00, 1.02,
     &     0.00, 1.00, 0.00, 1.03, 0.00, 1.04, 0.00, 1.22, 0.00, 0.98,
     &     0.00, 0.96, 0.00, 1.00, 0.00, 1.01, 0.00, 0.98, 0.00, 0.94,
     &     0.00, 0.92, 0.00, 0.89, 0.00, 0.86, 0.00, 0.84, 0.00, 0.82,
     &     0.00, 0.80, 0.00, 0.80, 0.00,
     &     0.96, 2.08, 0.98, 2.08, 0.98, 2.01, 0.98, 1.99, 1.05, 2.09,
     &     1.25, 2.31, 1.28, 2.52, 1.37, 2.64, 1.38, 2.40, 1.25, 2.15,
     &     1.15, 2.07, 1.06, 2.05, 1.04, 2.01, 1.04, 1.99, 1.04, 1.96,
     &     1.03, 1.95, 1.03, 1.91, 1.04, 1.89, 1.05, 1.88, 1.07, 1.90,
     &     1.10, 1.91, 1.13, 1.91, 1.15,
     &     0.00, 1.10, 0.00, 1.03, 0.00, 1.01, 0.00, 1.00, 0.00, 1.03,
     &     0.00, 1.05, 0.00, 1.24, 0.00, 1.00, 0.00, 0.96, 0.00, 1.01,
     &     0.00, 1.02, 0.00, 1.01, 0.00, 0.97, 0.00, 0.95, 0.00, 0.92,
     &     0.00, 0.89, 0.00, 0.85, 0.00, 0.82, 0.00, 0.81, 0.00, 0.80,
     &     0.00, 0.78, 0.00, 0.77, 0.00/
 
c Z = 68 to 70
      data c15
     &    /2.02, 0.89, 1.94, 0.91, 1.94, 0.99, 2.03, 1.13, 2.21, 1.18,
     &     2.26, 1.31, 2.57, 1.32, 2.34, 1.26, 2.16, 1.13, 2.08, 1.05,
     &     2.04, 1.02, 2.00, 1.02, 1.99, 1.02, 1.97, 1.01, 1.93, 1.01,
     &     1.89, 1.01, 1.87, 1.02, 1.85, 1.02, 1.85, 1.04, 1.85, 1.05,
     &     1.85, 1.08, 1.86, 1.10, 1.86,
     &     1.04, 0.00, 1.01, 0.00, 1.00, 0.00, 1.04, 0.00, 1.07, 0.00,
     &     1.25, 0.00, 1.01, 0.00, 0.97, 0.00, 1.00, 0.00, 1.02, 0.00,
     &     1.03, 0.00, 1.02, 0.00, 0.99, 0.00, 0.95, 0.00, 0.90, 0.00,
     &     0.86, 0.00, 0.85, 0.00, 0.83, 0.00, 0.82, 0.00, 0.78, 0.00,
     &     0.77, 0.00, 0.75, 0.00, 0.75,
     &     0.88, 1.89, 0.90, 2.16, 1.07, 2.12, 1.12, 2.22, 1.24, 2.54,
     &     1.26, 2.30, 1.19, 2.13, 1.11, 2.08, 1.05, 2.03, 1.01, 2.02,
     &     1.00, 2.00, 1.00, 1.96, 0.98, 1.93, 0.99, 1.90, 0.99, 1.87,
     &     1.00, 1.86, 1.00, 1.84, 1.00, 1.82, 1.01, 1.80, 1.03, 1.80,
     &     1.04, 1.81, 1.04, 1.80, 1.07/
 
c Z = 71 to 73
      data c16
     &    /0.00, 1.11, 0.00, 1.05, 0.00, 1.08, 0.00, 1.24, 0.00, 1.01,
     &     0.00, 0.98, 0.00, 1.00, 0.00, 1.02, 0.00, 1.03, 0.00, 1.02,
     &     0.00, 1.00, 0.00, 0.97, 0.00, 0.94, 0.00, 0.91, 0.00, 0.87,
     &     0.00, 0.85, 0.00, 0.83, 0.00, 0.79, 0.00, 0.78, 0.00, 0.78,
     &     0.00, 0.76, 0.00, 0.77, 0.00,
     &     1.04, 2.13, 1.06, 2.15, 1.16, 2.46, 1.19, 2.22, 1.12, 2.09,
     &     1.06, 2.04, 1.03, 2.01, 0.98, 2.00, 0.97, 1.98, 0.98, 1.97,
     &     0.99, 1.94, 0.98, 1.89, 0.98, 1.85, 0.96, 1.84, 0.96, 1.81,
     &     0.96, 1.80, 0.97, 1.77, 0.98, 1.77, 0.99, 1.76, 0.97, 1.72,
     &     1.00, 1.81, 1.12, 1.88, 1.15,
     &     1.07, 0.00, 1.09, 0.00, 1.24, 0.00, 1.03, 0.00, 1.00, 0.00,
     &     1.00, 0.00, 1.02, 0.00, 1.03, 0.00, 1.03, 0.00, 1.01, 0.00,
     &     0.97, 0.00, 0.93, 0.00, 0.92, 0.00, 0.89, 0.00, 0.87, 0.00,
     &     0.84, 0.00, 0.81, 0.00, 0.79, 0.00, 0.77, 0.00, 0.75, 0.00,
     &     0.78, 0.00, 0.77, 0.00, 0.68/
 
c Z = 74 to 76
      data c17
     &    /2.12, 1.07, 2.37, 1.10, 2.15, 1.02, 2.04, 0.98, 1.97, 0.93,
     &     1.96, 0.94, 1.97, 0.93, 1.96, 0.94, 1.96, 0.97, 1.94, 0.98,
     &     1.91, 0.99, 1.86, 0.96, 1.83, 0.95, 1.82, 0.93, 1.78, 0.92,
     &     1.73, 0.92, 1.69, 0.92, 1.69, 0.92, 1.69, 0.95, 1.74, 1.07,
     &     1.85, 1.09, 1.83, 1.14, 1.92,
     &     0.00, 1.27, 0.00, 1.04, 0.00, 1.04, 0.00, 1.02, 0.00, 1.02,
     &     0.00, 1.03, 0.00, 1.05, 0.00, 1.04, 0.00, 1.03, 0.00, 0.99,
     &     0.00, 0.96, 0.00, 0.93, 0.00, 0.90, 0.00, 0.87, 0.00, 0.84,
     &     0.00, 0.81, 0.00, 0.78, 0.00, 0.76, 0.00, 0.78, 0.00, 0.78,
     &     0.00, 0.73, 0.00, 0.77, 0.00,
     &     0.99, 2.05, 0.92, 2.05, 0.89, 1.90, 0.85, 1.92, 0.90, 1.93,
     &     0.89, 1.94, 0.90, 1.94, 0.94, 1.93, 0.98, 1.91, 0.96, 1.87,
     &     0.95, 1.83, 0.94, 1.82, 0.93, 1.78, 0.92, 1.75, 0.91, 1.75,
     &     0.91, 1.73, 0.93, 1.70, 0.93, 1.72, 0.95, 1.79, 1.01, 1.77,
     &     1.05, 1.82, 1.10, 1.98, 1.09/
 
c Z = 77 to 79
      data c18
     &    /0.00, 1.09, 0.00, 1.05, 0.00, 1.05, 0.00, 1.04, 0.00, 1.05,
     &     0.00, 1.06, 0.00, 1.06, 0.00, 1.04, 0.00, 0.93, 0.00, 0.92,
     &     0.00, 0.89, 0.00, 0.88, 0.00, 0.89, 0.00, 0.89, 0.00, 0.85,
     &     0.00, 0.81, 0.00, 0.79, 0.00, 0.79, 0.00, 0.74, 0.00, 0.76,
     &     0.00, 0.88, 0.00, 0.75, 0.00,
     &     0.79, 1.87, 0.78, 1.85, 0.83, 1.88, 0.83, 1.90, 0.88, 1.92,
     &     0.91, 1.94, 1.01, 1.92, 1.01, 1.90, 0.99, 1.88, 0.98, 1.84,
     &     0.97, 1.82, 0.95, 1.90, 0.89, 1.89, 0.90, 1.86, 0.90, 1.83,
     &     0.92, 1.76, 0.92, 1.71, 0.92, 1.68, 0.95, 1.73, 0.98, 1.87,
     &     0.98, 1.73, 0.97, 1.79, 0.96,
     &     1.08, 0.00, 1.08, 0.00, 1.08, 0.00, 1.08, 0.00, 1.08, 0.00,
     &     1.08, 0.00, 1.07, 0.00, 1.06, 0.00, 1.04, 0.00, 1.04, 0.00,
     &     1.03, 0.00, 1.03, 0.00, 1.01, 0.00, 0.98, 0.00, 0.92, 0.00,
     &     0.86, 0.00, 0.79, 0.00, 0.76, 0.00, 0.79, 0.00, 0.89, 0.00,
     &     0.76, 0.00, 0.82, 0.00, 0.87/
 
c Z = 80 to 82
      data c19
     &    /1.82, 0.76, 1.87, 0.79, 1.90, 0.82, 1.90, 0.82, 1.93, 0.86,
     &     1.95, 0.87, 1.94, 0.87, 1.93, 0.87, 1.94, 0.87, 1.94, 0.88,
     &     1.95, 0.91, 1.96, 0.94, 1.93, 0.95, 1.88, 0.95, 1.82, 0.95,
     &     1.75, 0.94, 1.68, 0.92, 1.69, 0.87, 1.77, 0.87, 1.64, 0.88,
     &     1.71, 0.88, 1.76, 0.88, 1.79,
     &     1.12, 0.00, 1.10, 0.00, 1.09, 0.00, 1.09, 0.00, 1.08, 0.00,
     &     1.08, 0.00, 1.08, 0.00, 1.09, 0.00, 1.09, 0.00, 1.09, 0.00,
     &     1.09, 0.00, 1.07, 0.00, 1.04, 0.00, 0.98, 0.00, 0.90, 0.00,
     &     0.82, 0.00, 0.75, 0.00, 0.90, 0.00, 0.76, 0.00, 0.83, 0.00,
     &     0.89, 0.00, 0.93, 0.00, 0.95,
     &     1.03, 2.14, 1.02, 2.12, 1.02, 2.11, 1.02, 2.10, 1.02, 2.09,
     &     1.02, 2.10, 1.03, 2.12, 1.03, 2.12, 1.02, 2.12, 1.01, 2.11,
     &     1.01, 2.10, 1.01, 2.06, 1.00, 2.01, 1.00, 1.93, 0.99, 1.83,
     &     0.98, 1.75, 1.00, 1.91, 0.99, 1.76, 0.99, 1.83, 0.99, 1.88,
     &     1.00, 1.88, 1.00, 1.87, 0.99/
 
c Z = 83 to 85
      data c20
     &    /0.00, 1.09, 0.00, 1.08, 0.00, 1.07, 0.00, 1.07, 0.00, 1.07,
     &     0.00, 1.08, 0.00, 1.09, 0.00, 1.09, 0.00, 1.09, 0.00, 1.07,
     &     0.00, 1.04, 0.00, 0.98, 0.00, 0.92, 0.00, 0.82, 0.00, 0.76,
     &     0.00, 0.90, 0.00, 0.77, 0.00, 0.84, 0.00, 0.90, 0.00, 0.79,
     &     0.00, 0.81, 0.00, 0.79, 0.00,
     &     1.07, 2.02, 1.05, 1.96, 1.04, 1.94, 1.04, 1.93, 1.07, 1.93,
     &     1.08, 1.85, 0.87, 2.01, 0.88, 2.01, 0.88, 2.00, 0.91, 1.96,
     &     0.90, 1.91, 0.90, 1.84, 0.90, 1.74, 0.89, 1.66, 0.89, 1.80,
     &     0.89, 1.66, 0.89, 1.73, 0.89, 1.72, 0.89, 1.68, 0.88, 1.63,
     &     0.91, 1.67, 0.95, 1.68, 0.98,
     &     0.00, 0.92, 0.00, 0.91, 0.00, 0.91, 0.00, 0.98, 0.00, 0.97,
     &     0.00, 0.97, 0.00, 0.96, 0.00, 1.00, 0.00, 0.98, 0.00, 0.98,
     &     0.00, 0.92, 0.00, 0.82, 0.00, 0.77, 0.00, 0.90, 0.00, 0.77,
     &     0.00, 0.84, 0.00, 0.77, 0.00, 0.73, 0.00, 0.73, 0.00, 0.73,
     &     0.00, 0.72, 0.00, 0.70, 0.00/
 
c Z = 86 to 88
      data c21
     &    /1.92, 0.98, 1.90, 1.00, 1.85, 1.05, 1.94, 1.10, 1.84, 0.88,
     &     1.83, 0.88, 1.84, 0.89, 1.83, 0.91, 1.90, 0.93, 1.86, 0.94,
     &     1.89, 0.97, 1.82, 0.97, 1.75, 0.99, 1.89, 0.98, 1.76, 0.98,
     &     1.83, 0.93, 1.70, 0.86, 1.58, 0.86, 1.59, 0.90, 1.63, 0.93,
     &     1.67, 0.98, 1.71, 1.01, 1.71,
     &     0.91, 0.00, 0.91, 0.00, 0.92, 0.00, 0.92, 0.00, 0.95, 0.00,
     &     0.94, 0.00, 0.93, 0.00, 0.89, 0.00, 0.88, 0.00, 0.82, 0.00,
     &     0.82, 0.00, 0.78, 0.00, 0.90, 0.00, 0.78, 0.00, 0.85, 0.00,
     &     0.72, 0.00, 0.72, 0.00, 0.72, 0.00, 0.73, 0.00, 0.72, 0.00,
     &     0.70, 0.00, 0.69, 0.00, 0.67,
     &     1.85, 1.00, 1.91, 1.08, 1.94, 0.97, 1.85, 0.90, 1.85, 0.92,
     &     1.86, 0.94, 1.84, 0.96, 1.84, 0.98, 1.81, 1.01, 1.85, 1.04,
     &     1.83, 1.06, 1.96, 1.06, 1.84, 1.03, 1.85, 0.90, 1.58, 0.85,
     &     1.57, 0.86, 1.59, 0.87, 1.61, 0.90, 1.66, 0.96, 1.68, 0.98,
     &     1.69, 1.00, 1.67, 1.01, 1.68/
 
c Z = 89 to 91
 
      data c22
     &    /0.89, 0.00, 0.90, 0.00, 0.94, 0.00, 0.94, 0.00, 0.93, 0.00,
     &     0.89, 0.00, 0.84, 0.00, 0.81, 0.00, 0.79, 0.00, 0.79, 0.00,
     &     0.90, 0.00, 0.78, 0.00, 0.80, 0.00, 0.72, 0.00, 0.73, 0.00,
     &     0.73, 0.00, 0.74, 0.00, 0.73, 0.00, 0.72, 0.00, 0.69, 0.00,
     &     0.68, 0.00, 0.67, 0.00, 0.66,
     &     1.91, 1.00, 1.85, 0.90, 1.88, 0.94, 1.88, 0.96, 1.86, 0.98,
     &     1.84, 1.00, 1.83, 1.03, 1.87, 1.08, 1.89, 1.12, 2.02, 1.11,
     &     1.90, 1.06, 1.83, 0.87, 1.58, 0.83, 1.57, 0.83, 1.62, 0.87,
     &     1.63, 0.91, 1.63, 0.91, 1.65, 0.93, 1.63, 0.95, 1.64, 0.97,
     &     1.65, 0.99, 1.65, 1.01, 1.69,
     &     0.00, 0.91, 0.00, 0.84, 0.00, 0.83, 0.00, 0.89, 0.00, 0.84,
     &     0.00, 0.83, 0.00, 0.78, 0.00, 0.80, 0.00, 0.90, 0.00, 0.79,
     &     0.00, 0.85, 0.00, 0.75, 0.00, 0.75, 0.00, 0.76, 0.00, 0.74,
     &     0.00, 0.73, 0.00, 0.72, 0.00, 0.69, 0.00, 0.68, 0.00, 0.67,
     &     0.00, 0.67, 0.00, 0.67, 0.00/
 
c Z = 92 to 94
      data c23
     &    /1.00, 1.83, 1.00, 1.84, 1.02, 1.84, 0.99, 1.86, 1.02, 1.86,
     &     1.06, 1.88, 1.11, 1.94, 1.15, 2.05, 1.15, 1.94, 1.03, 1.79,
     &     0.90, 1.61, 0.83, 1.64, 0.83, 1.62, 0.83, 1.60, 0.86, 1.59,
     &     0.87, 1.60, 0.88, 1.59, 0.89, 1.59, 0.92, 1.61, 0.95, 1.63,
     &     0.97, 1.67, 0.99, 1.68, 1.01,
     &     0.00, 0.82, 0.00, 0.80, 0.00, 0.85, 0.00, 0.84, 0.00, 0.78,
     &     0.00, 0.82, 0.00, 0.90, 0.00, 0.79, 0.00, 0.75, 0.00, 0.80,
     &     0.00, 0.82, 0.00, 0.79, 0.00, 0.76, 0.00, 0.74, 0.00, 0.73,
     &     0.00, 0.71, 0.00, 0.69, 0.00, 0.68, 0.00, 0.67, 0.00, 0.68,
     &     0.00, 0.68, 0.00, 0.69, 0.00,
     &     1.00, 1.83, 1.04, 1.87, 1.05, 1.89, 1.07, 1.89, 1.12, 1.97,
     &     1.18, 2.08, 1.17, 1.97, 1.02, 1.69, 0.91, 1.65, 0.85, 1.64,
     &     0.82, 1.61, 0.82, 1.58, 0.82, 1.59, 0.83, 1.57, 0.83, 1.57,
     &     0.86, 1.57, 0.88, 1.58, 0.90, 1.60, 0.93, 1.62, 0.94, 1.64,
     &     0.97, 1.66, 0.98, 1.67, 1.01/
 
c Z = 95 to 97
      data c24
     &    /0.00, 0.80, 0.00, 0.85, 0.00, 0.78, 0.00, 0.84, 0.00, 0.89,
     &     0.00, 0.79, 0.00, 0.77, 0.00, 0.87, 0.00, 0.83, 0.00, 0.81,
     &     0.00, 0.78, 0.00, 0.76, 0.00, 0.74, 0.00, 0.72, 0.00, 0.70,
     &     0.00, 0.69, 0.00, 0.69, 0.00, 0.69, 0.00, 0.69, 0.00, 0.69,
     &     0.00, 0.68, 0.00, 2*0.0,
     &     1.07, 1.91, 1.07, 1.89, 1.14, 2.00, 1.20, 2.10, 1.19, 1.95,
     &     1.16, 1.73, 0.92, 1.72, 0.84, 1.68, 0.83, 1.65, 0.82, 1.60,
     &     0.80, 1.57, 0.80, 1.57, 0.81, 1.54, 0.81, 1.55, 0.84, 1.53,
     &     0.85, 1.56, 0.88, 1.58, 0.89, 1.58, 0.91, 1.60, 0.94, 1.62,
     &     0.97, 1.62, 3*0.0,
     &     0.00, 0.79, 0.00, 0.85, 0.00, 0.90, 0.00, 0.79, 0.00, 0.88,
     &     0.00, 0.88, 0.00, 0.86, 0.00, 0.84, 0.00, 0.81, 0.00, 0.78,
     &     0.00, 0.76, 0.00, 0.74, 0.00, 0.72, 0.00, 0.70, 0.00, 0.69,
     &     0.00, 0.69, 0.00, 0.68, 0.00, 0.68, 0.00, 0.68, 0.00, 0.65,
     &     5*0.0/
 
c Z = 98 to 100
      data c25
     &    /1.13, 2.02, 1.20, 2.11, 1.20, 1.94, 1.17, 1.69, 0.81, 1.69,
     &     0.81, 1.70, 0.83, 1.67, 0.82, 1.64, 0.81, 1.61, 0.81, 1.58,
     &     0.81, 1.56, 0.80, 1.52, 0.80, 1.51, 0.81, 1.51, 0.83, 1.54,
     &     0.85, 1.54, 0.87, 1.56, 0.89, 1.59, 0.93, 1.59, 0.96, 1.62,
     &     5*0.0,
     &     0.00, 0.90, 0.00, 0.81, 0.00, 0.87, 0.00, 0.87, 0.00, 0.89,
     &     0.00, 0.87, 0.00, 0.84, 0.00, 0.81, 0.00, 0.78, 0.00, 0.76,
     &     0.00, 0.73, 0.00, 0.71, 0.00, 0.69, 0.00, 0.70, 0.00, 0.67,
     &     0.00, 0.69, 0.00, 0.67, 0.00, 0.65, 0.00, 0.65, 0.00, 0.65,
     &     5*0.0,
     &     2.11, 1.20, 1.95, 1.17, 1.67, 0.80, 1.67, 0.80, 1.68, 0.80,
     &     1.68, 0.80, 1.65, 0.80, 1.62, 0.80, 1.59, 0.80, 1.56, 0.79,
     &     1.53, 0.79, 1.50, 0.79, 1.49, 0.79, 1.51, 0.80, 1.50, 0.83,
     &     1.53, 0.85, 1.53, 0.87, 1.53, 0.90, 1.57, 0.93, 1.58, 0.94,
     &     5*0.0/
c ======================================================================
      end
