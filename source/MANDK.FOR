c       mandk
c	single coil form version
c	compute mutual inductance between a tilted pancake spiral primary
c	coil and a solenoidal secondary coil
c	this version steps through various coil positions!!!
c
c	Mark S. Rzeszotarski, Ph.D.
c
c	16-nov-95
c       this version: 29-oct-1997
c	fixed ds*dsp scaling error on 7-jun-98 - version 2.0
c	added operating characteristics section 7/4/98 - version 2.1
c	version 2.4 - 12/18/98 - changed figure of merit, cleanups.
c	version 2.5 - 1/15-99 - changed copyright to 1998-1999
c	version 2.6 - 09/29/99 - fixed N=1 solenoid error.
c	version 2.7 - 03/24/000 - added toroid table option.
c	version 2.8 - 09/06/2000 - fixed bug in Rproximity - was 
c					overestimating this.
c	version 3.0 - 12/29/2000 - added ini file for certain defaults
c				    added Lundin solenoid formula
c	version 3.1 - 12/31/2000 - allow pos1>pos2, change step size
c				   to reduce error slightly
c	version 3.2 - 01/10/2001 - improved accuracy of flat spiral calcs
c
c
	integer*2 itor
	real*8 xa,xb,ya,yb,za,zb,xm,ym,zm
	real*8 xap,xbp,yap,ybp,zap,zbp
	real*8 xmp(8190),ymp(8190),zmp(8190)
     1,steps,stepsp,dx,dy,dz,ct(8)
     1,dxp(8190),dyp(8190),dzp(8190)
	real*8 piece,sum,r,coseps,sum1
	real*8 xst,yst,zst,xstp,ystp,zstp,radian
c
c	read .ini file
c	
	htincr=0.5
	ctmax=35.
	isw1=0
	isw2=0
	open(2,file='mandk.ini', status='old',form='formatted', err=1234)
	read (2,1,err=1234) htincr
	read (2,1,err=1234) ctmax
	read (2,2,err=1234) isw1
	read (2,2,err=1234) isw2
    1   format(f16.8)
    2   format(i5)
	close(2)
	goto 1235
 1234   continue
	write (*,1233)
 1233   format(' Error reading (or locating) the file: mandk.ini',/
     1,' It should be placed in the same folder as mandk.exe',/
     1,' Default values will be used instead.')
 1235   continue
	ct(1)=0.
	do 3 j=2,8
    3   ct(j)=ctmax*float(j-1)/7.
	radian=6.283185307179/360.0
c
c	obtain parameters for primary coil
c
   32   write (*,690)
	read (*,*) itype
  690   format(' Mutual Inductance Calculation Program V3.2',/
     1,' Copyright (c) 1998, 1999, 2000 by Mark'
     1,' S. Rzeszotarski, Ph.D.',/,/
     1,' Specify Primary Coil Parameters:'
     1,/,' Select 3=flat spiral (pancake)'
     1,', 2=inverted cone (saucer), 1=solenoid ')
c  695   format(i5)
	if (itype.lt.1.or.itype.gt.3) goto 32
	if (itype.eq.1) then
		write (*,1700)
	else
		write (*,700)
	endif
	read (*,*) diam1
	elevat=0.
	if (itype.eq.1) then
                write (*,621)
		read (*,*) elevat
		diam2=diam1
	else
		write (*,702)
		read (*,*) diam2
	endif
	write (*,703)
	read (*,*) priN
	npri=nint(priN)
        if (itype.eq.2) then
		write (*,601)
		read (*,*) elevat
	endif
	d1=diam1
	d2=diam2
	h=elevat
  601   format(' Specify outer turn height in inches')
  621   format(' Specify primary coil height in inches')
	write (*,801)
	read (*,*) wd
  801   format(' Specify primary coil wire diameter in inches')
  700   format(' Specify inner primary diameter in inches')
 1700   format(' Specify primary coil diameter in inches')
c  701   format(f16.8)
  702   format(' Specify outer primary diameter in inches')
  703   format(' Specify number of primary coil turns of wire')
c  704   format(i6)
c
c	obtain parameters for solenoidal secondary coil
c
	write (*,705)
	Read (*,*) diams
	write (*,706)
	read (*,*) height
	write (*,707)
	read (*,*) secn
	nsec=nint(secn)
	write (*,802)
	read (*,*) wds
        write (*,508)
        read (*,*) pos1
	write (*,507)
        read (*,*) pos2
  508   format(' Specify elevation in inches for first iteration '
     1,/,' A negative value means the bottom wire of the'
     1,/,' secondary is below the bottom wire of the primary')
  507   format(' Specify elevation in inches for last iteration '
     1,/,' A positive value means that the bottom wire of the'
     1,/,' secondary is above the bottom wire of the primary coil')
c	
c	
cc	write (*,651)
cc  651   format(' Do you want a toroid capacitance table (1=yes, 0=no)')
cc	read (*,*) itor
	itor=isw1
c	
c	
	write (*,599)
  599   format(' Specify the file name for the output text data.',/
     1,' Recommended naming convention:  result.dat   ')
	write (1,600)
  600   format(' Mutual Inductance Program V3.2',/
     1,' Copyright (c) 1998, 1999, 2000 by Mark S. Rzeszotarski, Ph.D.')
  802   format(' Specify secondary coil wire diameter in inches')
  705   format(' Specify diameter of solenoidal secondary coil'
     1,' in inches')
  706   format(' Specify secondary height in inches')
  707   format(' Specify number of turns of wire in secondary coil')
  708   format(' Specify height in inches above plane of primary',/
     1,' for first turn of the secondary to be positioned')
c
c
c	solenoid
	write (*,490)
  490   format(/,/,' Extensive calculations are underway...',/
     1,'  Please be patient.  This will take several minutes.'/,)
	if (itype.eq.1) then
		write (1,501)
		write (*,501)
  501   format(/,' Solenoidal Primary Coil Geometry')
	endif
c	saucer
	if (itype.eq.2) then
		write (*,502)
		write (1,502)
  502   format(/,' Inverted Cone Primary Coil Geometry')
	endif
c	flat spiral
	if (itype.eq.3) then
		write (1,503)
		write (*,503)
  503   format(/,' Flat Spiral Primary Coil Geometry')
	endif
        ipos=0
	if (pos1.gt.pos2) htincr=-htincr
        do 8000 posit=pos1,pos2,htincr
        ipos=ipos+1
        secpos=posit
	diam1=diam1*2.54
	diam2=diam2*2.54
	elevat=elevat*2.54
	pitch=height/secn
	s=wds/pitch
	diams=diams*2.54
	height=height*2.54
	secpos=secpos*2.54
c
c	initial parameters for solenoidal secondary
c
	rad1=diams/2.0
c
c	allow maximum of 1% error in solenoid wire position
        pct=0.001
	pct=0.0025
        ss=2.0*sqrt(rad1**2-(rad1*(1.0-pct))**2)
	totdis=6.283185*rad1
        steps=totdis/ss
	nstep=nint(steps)
	nst=nstep
	theta1=360.0/nstep
	nstep=nstep*secn
c
c	initial parameters for primary coil
c
	radp1=diam1/2.0
	radp2=diam2/2.0
	thick=radp2-radp1
	totdis=6.283185*radp2
c	allow maximum of 1% error in outer turn wire position
        pct=0.001
	pct=0.002
        ss=2.0*sqrt(radp2**2-(radp2*(1.0-pct))**2)
   34   CONTINUE
        stepsp=totdis/ss
	nstepp=nint(stepsp)
	theta2=360.0/nstepp
	nstepp=nstepp*priN
c	write (*,445) nstep,nstepp
c  445   format(' nstep,nstepp=',2i6)
	IF (NSTEPP.GT.8190) THEN
		SS=SS*1.05
		GOTO 34
	ENDIF
	scale9=thick/float(nstepp-1)
	scale8=elevat/float(nstepp-1)
c	=================================================================
C	MAKE A TABLE OF PRIMARY LOOP PARAMETERS
C
c
c	initialize parameters to loop around primary
c
	xstp=radp1
	ystp=0.0
	zstp=0.0
	xap=xstp
	yap=ystp
	zap=zstp
	sum1=0.
	angle2=0.
	if (posit.eq.pos1) cm=0.
c
c	inner loop computes positions around the primary coil
c
	do 1001 j=1,nstepp
	angle2=angle2+theta2
	rad=radp1+scale9*float(j-1)
	xstp=rad*cos(angle2*radian)
	ystp=rad*sin(angle2*radian)
	zstp=scale8*float(j-1)
	xbp=xstp
	ybp=ystp
	zbp=zstp
	xmp(J)=(xap+xbp)/2.0
	ymp(J)=(yap+ybp)/2.0
	zmp(J)=(zap+zbp)/2.0
	dxp(J)=xbp-xap
	dyp(J)=ybp-yap
	dzp(J)=zbp-zap
	dl=sqrt((xbp-xap)**2+(ybp-yap)**2+(zbp-zap)**2)
	if (posit.eq.pos1) cm=cm+dl
	xap=xbp
	yap=ybp
	zap=zbp
 1001   continue
c	secondary wire length calculation
	if (posit.eq.pos1) then
		cccc=6.283185*rad1
		ccc=height/secn
		cms=secn*sqrt((cccc)**2+(ccc)**2)
	endif
c
C	=====================================================================
c
c	initialize parameters to loop around primary
c
	xstp=radp1
	ystp=0.0
	zstp=0.0
	xap=xstp
	yap=ystp
	zap=zstp
	sum1=0.
	angle2=0.
c
c	initialize parameters for looping
c
	sum=0.
	xst=rad1
	yst=0.0
	zst=secpos
	xa=xst
	ya=yst
	za=zst
	angle1=0.
c
c	outer loop handles turns around the secondary
c
	ito=0
	do 2000 i=1,nstep
	angle1=angle1+theta1
	xst=rad1*cos(radian*angle1)
	yst=rad1*sin(radian*angle1)
	zst=secpos+angle1/360.0*height/secn
	xb=xst
	yb=yst
	zb=zst
	xm=0.5*(xa+xb)
	ym=0.5*(ya+yb)
	zm=0.5*(za+zb)
	dx=xb-xa
	dy=yb-ya
	dz=zb-za
c
c	initialize parameters to loop around primary
c
c	xstp=radp1
c	ystp=0.0
c	zstp=0.0
c	xap=xstp
c	yap=ystp
c	zap=zstp
	sum1=0.
c	angle2=0.
c
c	inner loop computes positions around the primary coil
c
	do 1000 j=1,nstepp
C	angle2=angle2+theta2
C	rad=radp1+scale9*float(j-1)
C	xstp=rad*cos(angle2*radian)
C	ystp=rad*sin(angle2*radian)
C	zstp=scale8*float(j-1)
C	xbp=xstp
C	ybp=ystp
C	zbp=zstp
C	xmp=(xap+xbp)/2.0
C	ymp=(yap+ybp)/2.0
C	zmp=(zap+zbp)/2.0
C	dxp=xbp-xap
C	dyp=ybp-yap
C	dzp=zbp-zap
	r=sqrt((xm-xmp(J))**2+(ym-ymp(J))**2+(zm-zmp(J))**2)
	coseps=(dx*dxp(J)+dy*dyp(J)+dz*dzp(J))
	piece=coseps/r
	sum1=sum1+piece
C	xap=xbp
C	yap=ybp
C	zap=zbp
 1000   continue
	sum=sum+sum1
	xa=xb
	ya=yb
	za=zb
 2000   continue
c	convert to microhenries
c        write (1,554) sum,ds,dsp
c  554   format(' sum,ds,dsp:',3f16.8)
c	sum=sum*ds*dsp/1000.
	sum=sum/1000.
	diam1=diam1/2.54
	diam2=diam2/2.54
	height=height/2.54
	diams=diams/2.54
	secpos=secpos/2.54
	elevat=elevat/2.54
	turns=secn
	if (ipos.eq.1) then
	if (itype.ne.3) then
	write (*,9) diam1,diam2,priN,elevat
     1,wd,diams,height,secn,wds
	write (1,9) diam1,diam2,priN,elevat
     1,wd,diams,height,secn,wds
	else
	write (*,49) diam1,diam2,priN
     1,wd,diams,height,secn,wds
	write (1,49) diam1,diam2,priN
     1,wd,diams,height,secn,wds
	endif
	endif
    9   format(' Primary coil inside diameter (inches)= '
     1,f10.3,/
     1,' Primary coil outside diameter (inches)= ',f10.3,/
     1,' Number of primary coil turns = ',f10.3,/
     1,' Last turn elevation in inches = ',f10.3,/
     1,' Wire diameter (inches)= ',f10.4,/,/
     1,' Solenoidal Secondary Coil Geometry',/
     1,' Secondary coil diameter (inches)=',f10.3,/
     1,' Secondary coil height (inches)= ',f10.3,/
     1,' Number of secondary coil turns =',f10.3,/
     1,' Secondary coil wire diameter (inches)= ',f10.4,/)
   49   format(' Primary coil inside diameter (inches)= '
     1,f10.3,/
     1,' Primary coil outside diameter (inches)= ',f10.3,/
     1,' Number of primary coil turns = ',f10.3,/
     1,' Wire diameter (inches)= ',f10.4,/,/
     1,' Solenoidal Secondary Coil Geometry',/
     1,' Secondary coil diameter (inches)=',f10.3,/
     1,' Secondary coil height (inches)= ',f10.3,/
     1,' Number of secondary coil turns =',f10.3,/
     1,' Secondary coil wire diameter (inches)= ',f10.4,/)
	if (ipos.eq.1) then
	if (itype.NE.1) call spiral(henry1,d1,d2,priN,wd)
c	if (itype.eq.1) call wheelr(henry2,diam2
c     1,wd/elevat*priN,elevat,priN)
	if (itype.eq.1) call lundin(henry2,diam2,elevat,priN)
c        write (1,555) henry2,henry4,diam1,wd,h,priN
c  555   format(' henry2,henry4,diam1,wd,h,priN',/,1x,6f16.8)
c	if (itype.eq.2) call wheelr(hnry2,(d1+d2)/2.0
c     1,wd/h*priN,h,priN)
	if (itype.eq.2) call lundin(hnry2,(d1+d2)/2.0
     1,h,priN)
	if (itype.eq.2) call malcolm(henry3,d1,d2
     1,h,priN,henry1,hnry2)
c        write (1,553) henry1,henry2,henry3
c  553   format(' henry1,henry2,henry3',3f16.8)
	if (itype.eq.1) henry=henry2
	if (itype.eq.2) henry=henry3
	if (itype.eq.3) henry=henry1
	write (*,11) henry
	write (1,11) henry
   11   format(' Calculational Results'
     1,/,' Primary coil inductance in microhenries:  Lp =',f16.2)
c	call wheelr(henrys,diams,s,height,turns)
	call lundin(henrys,diams,height,turns)
c       write (1,555) henrys,henry4,diams,wd,height,turns
	write (*,12) henrys
	write (1,12) henrys
   12   format(' Secondary coil inductance in microhenries:  Ls ='
     1,f16.2)
	fom=sqrt(henrys/henry)
	write (*,113) fom
	write (1,113) fom
  113   format(' Figure of Merit (square'
     1,' root of Ls/Lp) =',f10.2)
	feetp=(cm/2.54)/12.0
	feets=(cms/2.54)/12.0
	write (*,491) feetp,feets
	write (1,491) feetp,feets
  491   format(' Primary coil wire length in feet = 'f10.2,/
     1,' Secondary coil wire length in feet = 'f10.2)
c
c	compute wire size info for secondary
	amx=1.e30
c	find nearest wire gauge
	do 322 j=1,40
	call wire(j,awg,ohms,warea,wdiam,weight)
	if (abs(wds-wdiam).lt.amx) then
		amx=abs(wds-wdiam)
		iwg=j
	endif
  322   continue
	call wire(iwg,awg,ohms,warea,wdiam,weight)
	Rdc=ohms*feets/1000.
	write (*,323) Rdc
	write (1,323) Rdc
  323   FORMAT(' DC secondary resistance in ohms:',f10.2)
c
	call cap(capsp,height,diams)
	write (*,13) capsp
	write (1,13) capsp
   13   format(' Secondary coil distributed capacitance in'
     1,' picofarads:',f10.2,/
     1,' (Medhurst formula, assumes one end of secondary is grounded)'
     1,/,/,' Mutual Inductance Results'
     1,/,' Position is the secondary coil bottom wire'
     1,' position in inches'
     1,/,' above the bottom wire of the primary coil.'
     1,/,' A negative value means the bottom wire of the'
     1,/,' secondary is below the bottom wire of the primary.'
     1,/,' A positive value means that the bottom wire of the'
     1,/,' secondary is above the bottom wire of the primary coil.'
     1,/,' M = Mutual Inductance in microhenries'
     1,/,' K = Coefficient of Coupling: K =  M /'
     1,' square root ( Lp x Ls )',/,/
     1,8x,'   Position           M                K')
	endif
	couple=sum/sqrt(henry*henrys)
	write (*,10) secpos,sum,couple
	write (1,10) secpos,sum,couple
   10   format(1x,f16.3,f16.3,f16.4)
 8000   continue
	c1=ct(1)+capsp
	c2=ct(2)+capsp
	c3=ct(3)+capsp
	c4=ct(4)+capsp
	c5=ct(5)+capsp
	c6=ct(6)+capsp
	c7=ct(7)+capsp
	c8=ct(8)+capsp
	write (*,888)
	write (1,888)
  888   format(/,/,' Operating Characteristics',/
     1,' Ctoroid is the added capacitance of the top terminal.',/
     1,' Fres is the resonant operating frequency in kilohertz.',/
     1,' Cpri is the required primary capacitance in nanofarads.',/
     1,' Zpri is the primary surge impedance in ohms.',/
     1,' Zsec is the secondary surge impedance in ohms.',/
     1,' Qsec is the secondary coil quality factor Q.',/
     1,' Rac is the secondary coil RF resistance in ohms.')
	if (isw2.eq.1) write (*,781)
	if (isw2.eq.1) write (1,781)
  781    format(' Lambda is the secondary wire 1/4 wavelength'
     1,' in feet')
c     1,12x,'Ct=0pF  Ct=5pF Ct=10pF'
c     1,'  Ct=15pF Ct=20pF Ct=30pF Ct=40pF Ct=50pF')
	write (*,789) (ct(j),j=1,8)
	write (1,789) (ct(j),j=1,8)
  789   format(/,' Ctoroid(pF)',f6.1,7f8.1)
	f1=1./(6.283185*sqrt(henrys*c1*1.e-18))/1000.
	f2=1./(6.283185*sqrt(henrys*c2*1.e-18))/1000.
	f3=1./(6.283185*sqrt(henrys*c3*1.e-18))/1000.
	f4=1./(6.283185*sqrt(henrys*c4*1.e-18))/1000.
	f5=1./(6.283185*sqrt(henrys*c5*1.e-18))/1000.
	f6=1./(6.283185*sqrt(henrys*c6*1.e-18))/1000.
	f7=1./(6.283185*sqrt(henrys*c7*1.e-18))/1000.
	f8=1./(6.283185*sqrt(henrys*c8*1.e-18))/1000.
	write (*,889) f1,f2,f3,f4,f5,f6,f7,f8
	write (1,889) f1,f2,f3,f4,f5,f6,f7,f8
  889   format(' Fres(kHz)',f8.1,7f8.1)
	cp1=1./(henry*1.e-6)*(1./(6.283185*f1*1000.))**2*1.e6
	cp2=1./(henry*1.e-6)*(1./(6.283185*f2*1000.))**2*1.e6
	cp3=1./(henry*1.e-6)*(1./(6.283185*f3*1000.))**2*1.e6
	cp4=1./(henry*1.e-6)*(1./(6.283185*f4*1000.))**2*1.e6
	cp5=1./(henry*1.e-6)*(1./(6.283185*f5*1000.))**2*1.e6
	cp6=1./(henry*1.e-6)*(1./(6.283185*f6*1000.))**2*1.e6
	cp7=1./(henry*1.e-6)*(1./(6.283185*f7*1000.))**2*1.e6
	cp8=1./(henry*1.e-6)*(1./(6.283185*f8*1000.))**2*1.e6
	cp1=cp1*1000.
	cp2=cp2*1000.
	cp3=cp3*1000.
	cp4=cp4*1000.
	cp5=cp5*1000.
	cp6=cp6*1000.
	cp7=cp7*1000.
	cp8=cp8*1000.
	fom2=sqrt((cp2*1.e-9)/(c2*1.e-12))
	fom3=sqrt((cp3*1.e-9)/(c3*1.e-12))
	fom4=sqrt((cp4*1.e-9)/(c4*1.e-12))
	fom5=sqrt((cp5*1.e-9)/(c5*1.e-12))
	fom6=sqrt((cp6*1.e-9)/(c6*1.e-12))
	fom7=sqrt((cp7*1.e-9)/(c7*1.e-12))
	fom8=sqrt((cp8*1.e-9)/(c8*1.e-12))
	write (*,887) cp1,cp2,cp3,cp4,cp5,cp6,cp7,cp8
	write (1,887) cp1,cp2,cp3,cp4,cp5,cp6,cp7,cp8
  887   format(' Cpri (nF)',f8.2,7f8.2)
	cp1=cp1/1000.
	cp2=cp2/1000.
	cp3=cp3/1000.
	cp4=cp4/1000.
	cp5=cp5/1000.
	cp6=cp6/1000.
	cp7=cp7/1000.
	cp8=cp8/1000.
	zp1=sqrt(henry/cp1)
	zp2=sqrt(henry/cp2)
	zp3=sqrt(henry/cp3)
	zp4=sqrt(henry/cp4)
	zp5=sqrt(henry/cp5)
	zp6=sqrt(henry/cp6)
	zp7=sqrt(henry/cp7)
	zp8=sqrt(henry/cp8)
	write (*,886) zp1,zp2,zp3,zp4,zp5,zp6,zp7,zp8
	write (1,886) zp1,zp2,zp3,zp4,zp5,zp6,zp7,zp8
  886   format(' Zpri(ohms)',f7.1,7f8.1)
c
	zp1=sqrt(henrys/c1*1.e6)
	zp2=sqrt(henrys/c2*1.e6)
	zp3=sqrt(henrys/c3*1.e6)
	zp4=sqrt(henrys/c4*1.e6)
	zp5=sqrt(henrys/c5*1.e6)
	zp6=sqrt(henrys/c6*1.e6)
	zp7=sqrt(henrys/c7*1.e6)
	zp8=sqrt(henrys/c8*1.e6)
	write (*,686) zp1,zp2,zp3,zp4,zp5,zp6,zp7,zp8
	write (1,686) zp1,zp2,zp3,zp4,zp5,zp6,zp7,zp8
  686   format(' Zsec(ohms)',f8.0,7f8.0)
c
c	compute wire size info for secondary
	amx=1.e30
c	find nearest wire gauge
	do 22 j=1,40
	call wire(j,awg,ohms,warea,wdiam,weight)
	if (abs(wds-wdiam).lt.amx) then
		amx=abs(wds-wdiam)
		iwg=j
	endif
   22   continue
	call wire(iwg,awg,ohms,warea,wdiam,weight)
	Rdc=ohms*feets/1000.
	ss=height/secn
c	write (*,21) iwg,wdiam,wds
c   21   format(' iwg,wdiam,wds,ss:',i4,3f10.4)
	call getRac(Rskin,Rdc,F1*1000.,wdiam)
	call getRpx(Rprox,Rskin,height,diams,wds,ss)
c	Rac=Rskin+Rprox
	Rac = Rprox
	q1=6.283185307*f1*1000.*henrys*1.e-6/Rac
	r1=Rac
c	write (1,24) Rdc,Rskin,Rprox,Rac,F1
c	write (*,24) Rdc,Rskin,Rprox,Rac,F1
c   24   format(' Rdc, Rskin, Rprox, Rac, Fres:',5f10.3)
c	write (*,27) ohms,feets
c	write (1,27) ohms,feets
c   27   format(' ohms=',f10.2,'    feets=',f10.2)
	call getRac(Rskin,Rdc,F2*1000.,wdiam)
	call getRpx(Rprox,Rskin,height,diams,wds,ss)
c	Rac=Rskin+Rprox
	Rac = Rprox
	r2=Rac
	q2=6.283185307*f2*1000.*henrys*1.e-6/Rac
	call getRac(Rskin,Rdc,F3*1000.,wdiam)
	call getRpx(Rprox,Rskin,height,diams,wds,ss)
	Rac=Rskin+Rprox
	Rac = Rprox
	q3=6.283185307*f3*1000.*henrys*1.e-6/Rac
	r3=Rac
	call getRac(Rskin,Rdc,F4*1000.,wdiam)
	call getRpx(Rprox,Rskin,height,diams,wds,ss)
	Rac=Rskin+Rprox
	Rac = Rprox
	q4=6.283185307*f4*1000.*henrys*1.e-6/Rac
	r4=Rac
	call getRac(Rskin,Rdc,F5*1000.,wdiam)
	call getRpx(Rprox,Rskin,height,diams,wds,ss)
	Rac=Rskin+Rprox
	Rac = Rprox
	q5=6.283185307*f5*1000.*henrys*1.e-6/Rac
	r5=Rac
	call getRac(Rskin,Rdc,F6*1000.,wdiam)
	call getRpx(Rprox,Rskin,height,diams,wds,ss)
	Rac=Rskin+Rprox
	Rac = Rprox
	q6=6.283185307*f6*1000.*henrys*1.e-6/Rac
	r6=Rac
	call getRac(Rskin,Rdc,F7*1000.,wdiam)
	call getRpx(Rprox,Rskin,height,diams,wds,ss)
	Rac=Rskin+Rprox
	Rac = Rprox
	q7=6.283185307*f7*1000.*henrys*1.e-6/Rac
	r7=Rac
	call getRac(Rskin,Rdc,F8*1000.,wdiam)
	call getRpx(Rprox,Rskin,height,diams,wds,ss)
	Rac=Rskin+Rprox
	Rac = Rprox
	q8=6.283185307*f8*1000.*henrys*1.e-6/Rac
	r8=Rac
	write (*,23) q1,q2,q3,q4,q5,q6,q7,q8
	write (1,23) q1,q2,q3,q4,q5,q6,q7,q8
   23   format(' Qsec     ',8f8.1)
	write (*,28) r1,r2,r3,r4,r5,r6,r7,r8
	write (1,28) r1,r2,r3,r4,r5,r6,r7,r8
   28   format(' Rac (ohms)',f7.1,7f8.1)
	r1=300000000./(f1*1000.)*39.36/12.0/4.0
	r2=300000000./(f2*1000.)*39.36/12.0/4.0
	r3=300000000./(f3*1000.)*39.36/12.0/4.0
	r4=300000000./(f4*1000.)*39.36/12.0/4.0
	r5=300000000./(f5*1000.)*39.36/12.0/4.0
	r6=300000000./(f6*1000.)*39.36/12.0/4.0
	r7=300000000./(f7*1000.)*39.36/12.0/4.0
	r8=300000000./(f8*1000.)*39.36/12.0/4.0
	if (isw2.eq.1) write (*,29) r1,r2,r3,r4,r5,r6,r7,r8
	if (isw2.eq.1) write (1,29) r1,r2,r3,r4,r5,r6,r7,r8
   29   format(' Lambda    ',f7.1,7f8.1)
c	write (*,41) fom1,fom2,fom3,fom4,fom5,fom6,fom7,fom8
c	write (1,41) fom1,fom2,fom3,fom4,fom5,fom6,fom7,fom8
c   41   format(' F.O.M.    ',f7.1,7f8.1)
c
c	
	if (itor.eq.1) call toroid
	write (*,291)
  291   format(/,/,' Calculations completed.',/
     1,' The results text file can now be reviewed with any'
     1,' text viewer.')
	stop
	end
c
c
c	spiral
c	computes inductance of a spiral-wound coil using round wire
c	or tubing
c
c	Mark S. Rzeszotarski, Ph.D.
c
c	31 oct 95
c
c	reference: Radio Designer's Handbook pp 445-6
c
c	all units in inches please!
c	result is in microhenries
c
	subroutine spiral(henry,d1,d2,enn,wd)
c	write (*,111) d1,d2,n,wd
c  111   format(1x,' d1,d2,n,wd:',2f12.3,i6,f12.4)
c	compute mean radius
	a=(d1+d2)/4.0
c	compute width
	c=d2/2.0-d1/2.0
	if (c.ge.0.2*a) then
	sheet=a*enn/(8.0*a+11.0*c)*a*enn
	else
	sheet=0.03193*2.303*a*enn*enn*(alog10(8.0*a/c)-0.5/2.303)
	endif
c	
c	stephan's formula
c	
	sheet=0.03193*a*enn*enn*(2.303*alog10(8.0*a/c)-0.5
     1+(c*c/(96.*a*a))*(2.303*alog10(8.*a/c)+(43./12.)))
	p=c/enn
	aa=2.3*alog10(1.7*wd/p)
	bb=0.336*(1.0-((2.5/enn)+3.8/(enn**2)))
	if (enn.le.1.0) bb=0.0
	if (enn.gt.1.0.and.enn.le.2.0) bb=0.113*(enn-1.)
	if (enn.gt.2.0.and.enn.le.3.0) bb=0.050*(enn-2.)+0.115
	if (enn.gt.3.0.and.enn.lt.4.0) bb=0.030*(enn-3.)+0.165
c	henry=sheet
c	if (enn.lt.4.) return
	henry=sheet-0.0319*a*enn*(aa+bb)
ccc	write (*,222) sheet,aa,bb,henry
ccc  222   format(' Ls, A, B, L0:',4f16.8)
ccc	pause
c	
c	use sheet formula (Stephan's equation) without corrections
	henry=sheet
	return
	end
c
c
C	WHEELER
C	COMPUTES INDUCTANCE USING WHEELERS FORMULA
	subroutine wheelr(result,diam,s,height,turns)
	radius=diam/2.0
	pi=3.14159265
	CALL GETk(ak,diam/height)
	call getab(a,b,turns,s)
	als=0.10028*(radius*radius)
	als=als*(turns*turns)
	als=als*ak/height
	result=als*(1.0-(height*(a+b))/(pi*radius*ak*turns))
	return
	END
C
	SUBROUTINE GETk(ak,doverl)
	DIMENSION R(16),hh(16)
	data r/100.,50.,20.,10.,5.,2.,1.,.75,.5,.4,.3,.2,.1,.05
     1,.02,.01/
        data hh/.03496,.061098,.123615,.203324
     1,.319825,.52551,.688423,.747762,.818136,.849853,.883803
     1,.920093,.958807,.979092,.991562,.995769/
	if (doverl.lt.r(16)) goto 100
	if (doverl.gt.r(1)) goto 200
c
c	interpolate between two points
c
	do 10 j=1,16
	if (doverl.eq.r(j)) then
		ak=hh(j)
		return
	endif
   10   continue
	do 15 j=1,15
	if (doverl.le.r(j).and.doverl.gt.r(j+1)) then
		ak=hh(j)*((doverl-r(j+1))/(r(j)-r(j+1)))
     1+hh(j+1)*((r(j)-doverl)/(r(j)-r(j+1)))
		return
	endif
   15   continue
c
  100   continue
	ak=1.0
	return
c
  200   continue
        ak=0.0350*(1.0-(0.01-doverl)/0.01)
c	added next line on 9/29/99
	ak=0.03496
	return
	END
c
c
	subroutine getab(a,b,turns,s)
	DIMENSION R(13),hh(13)
	data r/1000.,500.,200.,100.,50.,20.,10.,7.,5.,4.,3.,2.,1./
	data hh/0.34,0.335,.333,.329,.318,.296,.267,.245,.218
     1,.195,.165,.113,0.0/
c	assume close wound
	a=2.3*alog10(1.7*s)
c	next line added on 9/29/99
	b=0.0
	if (turns.lt.r(13)) goto 100
	if (turns.gt.r(1)) goto 200
        if (turns.ge.4.0) then
                A=2.3*alog10(1.7*s)
                B=0.336*(1.0-(2.5/turns+3.8/(turns)**2))
                return
        endif
c
c	interpolate between two points
c
	do 10 j=1,13
	if (turns.eq.r(j)) then
		b=hh(j)
		return
	endif
   10   continue
	do 15 j=1,12
	if (turns.le.r(j).and.turns.gt.r(j+1)) then
		b=hh(j)*((turns-r(j+1))/(r(j)-r(j+1)))
     1+hh(j+1)*((r(j)-turns)/(r(j)-r(j+1)))
		return
	endif
   15   continue
c
c	very short height/diam ratio
c
  100   continue
	b=0.0
	return
c
c	very high height/diam ratio
c
  200   continue
	b=0.34
	return
	END
C	MEDHURST
C	COMPUTES DISTRIBUTED CAPACITANCE OF A TESLA COIL USING
C	MEDHURST'S FORMULA
	subroutine cap(c,height,diam)
	CALL GETmH(H,HEIGHT,DIAM)
	c=H*DIAM*2.54
c	write (*,10) c,height,diam
c   10   format(' C, H, D:',3f16.8)
	return
	END
C
	SUBROUTINE GETmH(H,HEIGHT,DIAMET)
	DIMENSION R(33),hh(33)
	data r/50.,40.,30.,25.,20.,15.,10.,9.,8.,7.,6.,5.,4.5,4.,3.5
     1,3.,2.5,2.,1.5,1.,.9,.8,.7,.6,.5,.45,.4,.35,.3,.25,.2,.15,.1/
	data hh/5.8,4.6,3.4,2.9,2.36,1.86,1.32,1.33,1.12,1.01,.92
     1,.81,.77,.72,.67,.61,.56,.5,.47,3*.46,.47,.48,.5,.52,.54
     1,.57,.6,.64,.7,.79,.96/
	RATIO=HEIGHT/DIAMET
	if (ratio.lt.r(33)) goto 100
	if (ratio.gt.r(1)) goto 200
c
c	interpolate between two points
c
	do 10 j=1,33
	if (ratio.eq.r(j)) then
		h=hh(j)
		return
	endif
   10   continue
	do 15 j=1,32
	if (ratio.le.r(j).and.ratio.gt.r(j+1)) then
		h=hh(j)*((ratio-r(j+1))/(r(j)-r(j+1)))
     1+hh(j+1)*((r(j)-ratio)/(r(j)-r(j+1)))
		return
	endif
   15   continue
c
c	very short height/diam ratio
c
  100   continue
	h=1.0
d	write (*,11)
d  11   format(' No data available for this height/diameter ratio')
	return
c
c	very high height/diam ratio
c
  200   continue
	h=ratio/10.*1.2
d	write (*,11)
	return
	END
c
c
	subroutine malcolm(henryp,diam1,diam2,elevat,priN,h1,h2)
	r=(diam1+diam2)/4.0
	c=(diam2-diam1)/2.0
	h=elevat
	henryp=priN**2*r**2/(8.0*r+10.0*h+11.0*c)
	r1=diam1/2.0
	r2=diam2/2.0
	thick=r2-r1
	angle=atan(elevat/thick)
	f1=1.0-angle/1.5708
	f2=1.0-f1
	henryp=h1*f1+h2*f2
	return
	end
c
c
	subroutine wire(iawg,awg,ohms,warea,wdiam,weight)
c	wire
c	generates wire table
c
c	mark s rzeszotarski, phd
c
c
c	write (1,5)
c    5   format(' A.W.G.  Ohms/1000 ft   Area (mils)'
c     1,'   Diameter (mils)   Pounds/1000 ft')
	wg=float(iawg)
	awg=wg
	rft=0.0001/(0.793**wg)
	rmft=1000.*rft
	ohms=rmft
	wx=106000./(1.261**wg)
	warea=wx
	wd=0.325/(1.123**wg)
	wdiam=wd
	wt=wx*0.00327
	weight=wt
c	write (1,20) iawg,rmft,wx,wd,wt
c   20   format(1x,i4,3x,f10.2,5x,f10.2,4x,f10.4,9x,f10.3)
	return
	end
c
c	subroutine getRac(Rac,Rdc,Freq)
c	x=sqrt(Freq/Rdc)
c	if (x.lt.25.) then
c		Rac=1.0*Rdc
c		return
c	endif
c	if (x.ge.180.) then
c		ratio=1.25+(x-100.0)/(180.0-100.0)*0.75
c		Rac=Rdc*ratio
cc		write (*,10)
c		return
c	endif
cc   10   format(1x,'Error in computation')
c	if (x.ge.25.0.and.x.lt.50.0) then
c		ratio=1.00+(x-25.0)/(50.0-25.0)*0.02
c		Rac=Rdc*ratio
c		return
c	endif
c	if (x.ge.50.0.and.x.lt.75.0) then
c		ratio=1.02+(x-50.0)/(75.0-50.0)*0.08
c		Rac=Rdc*ratio
c		return
c	endif
c	if (x.ge.75.0.and.x.lt.100.0) then
c		ratio=1.10+(x-75.0)/(100.0-75.0)*0.15
c		Rac=Rdc*ratio
c		return
c	endif
c	if (x.ge.100.0.and.x.lt.180.0) then
c		ratio=1.25+(x-100.0)/(180.0-100.0)*0.75
c		Rac=Rdc*ratio
c		return
c	endif
c	return
c	end
c
c
	subroutine getRpx(Rprox,Rskin,height,diams,wdiam,ss)
	call proxim(Rprx,Height,Diams,wdiam,ss)
	Rprox=rprx*Rskin
	return
	end
c
c	proxim
c	computes proximity effect
c
c	Rprox is proximity effect resistance ratio
c	Height of solenoid in inches
c	Diam of solenoid in inches
c	wdiam is diameter of wire in inches
c	wspace is wire turn-to-turn spacing in inches
c
	subroutine proxim(Rprox,Height,Diam,wdiam,wspace)
	dimension table(11,10),htod(11),ds(10)
	data ds/1.,.9,.8,.7,.6,.5,.4,.3,.2,.1/
	data htod/0.0,0.2,0.4,0.6,0.8,1.0,2.0,4.0,6.0,8.0,10.0/
	data table/5.31
     1,5.45,5.65,5.80,5.80,5.55,4.10,3.54,3.31,3.20,3.23
     1,3.73,3.84,3.99,4.11,4.17,4.10,3.30,3.05,2.92,2.90,2.93
     1,2.74,2.83,2.97,3.10,3.20,3.17,2.74,2.60,2.60,2.62,2.65
     1,2.12,2.20,2.28,2.38,2.44,2.47,2.32,2.27,2.29,2.34,2.37
     1,1.74,1.77,1.83,1.89,1.92,1.94,1.98,2.01,2.03,2.08,2.10
     1,1.44,1.48,1.54,1.60,1.64,1.67,1.74,1.78,1.80,1.81,1.83
     1,1.26,1.29,1.33,1.38,1.42,1.45,1.50,1.54,1.56,1.57,1.58
     1,1.16,1.19,1.21,1.22,1.23,1.24,1.28,1.32,1.34,1.34,1.35
     1,1.07,1.08,1.08,1.10,1.10,1.10,1.13,1.15,1.16,1.16,1.17
     1,1.02,1.02,1.03,1.03,1.03,1.03,1.04,1.04,1.04,1.04,1.04/
	dovers=wdiam/wspace
	if (dovers.lt.0.1) then
		Rprox=1.0
		return
	endif
	if (dovers.gt.1.0) then
c		write (*,10)
		dovers=1.0
c   10   format(' Error - wire diameter greater than wire spacing'
c     1,/,' Setting wire spacing equal to wire diameter')
	endif
	hoverd=height/diam
cc	if (hoverd.gt.10.0.or.hoverd.lt.0.0) then
cc		write (*,11)
cc		stop
cc   11   format(' Error - H/D ratio exceeds 10:1')
cc	endif
	if (hoverd.gt.10.) hoverd=10.
	do 12 j=10,1,-1
	if (hoverd.ge.htod(j)) then
		loc1=j
		loc2=j+1
		h1=htod(loc1)
		h2=htod(loc2)
		p1=(hoverd-h1)/(h2-h1)
		p2=1.0-p1
		goto 13
	endif
   12   continue
   13   continue
	do 20 j=2,10
	if (dovers.ge.ds(j)) then
		iloc1=j
		iloc2=j-1
		ds1=ds(iloc1)
		ds2=ds(iloc2)
		d1=(dovers-ds2)/(ds1-ds2)
		d2=1.0-d1
		goto 21
	endif
   20   continue
   21   continue
	a=table(loc1,iloc2)*p2*d2
	b=table(loc1,iloc1)*p2*d1
	c=table(loc2,iloc1)*p1*d1
	d=table(loc2,iloc2)*p1*d2
	Rprox=a+b+c+d
c	write (*,22) hoverd,dovers,p1,p2,d1,d2,loc1,loc2
c     1,iloc1,iloc2,rprox,a,b,c,d
c   22   format(1x,'H/D=',f7.2,'  D/S=',f7.2,' p1,p2,d1,d2:',4f6.2
c     1,/,'  loc1,loc2,iloc1,iloc2',4i4,'   Rprox=',f7.2
c     1,/,1x,4f7.2)
c	pause
	return
	end
c
c
c	toroid
c	makes a table of toroid capacitance values
c	
c	reference: Australian Journal of Physics, V7 1953, pp 347-350
c	
c	
	subroutine toroid
	real*8 value(12),c(12)
	data value/1.,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12./
c	open(unit=1,file='toroid.dat')
	write (1,10) (value(j),j=1,12)
   10   format(/,/,/
     1,' Table of Toroid Capacitance in picofarads (pF)',/
     1,' Reference: Australian Journal of Physics, V7 1953,'
     1,' pp 347-350',/
     1,' 1st Column Is Toroid Outside Diameter In Inches',/,/
     1,24x,' Torus Circle Diameter in Inches',/,1x,6x,12f6.0)
	do 100 k=4,100,2
	diam=float(k)
	do 20 j=1,12
	d1=float(j)*2.54
	d2=diam*2.54
	ratio=d1/d2
	if (ratio.lt.0.25) then
	c(j)=1.8*(d2-d1)/(alog((8.0*(d2-d1)/d1)))
	else
	c(j)=0.37*d2+0.23*d1
	endif
   20   continue
	write (1,30) diam,(c(j),j=1,12)
   30   format(1x,f6.0,12f6.2)
  100   continue
	return
	end
	subroutine getRac(Rac,Rdc,Freq,wdiam)
c	
c	Rdc is the dc resistance of the wire (required input)
c	Freq is the frequency in hertz (required input)
c	wdiam is the required wire diameter in inches
c	xx is the parameter x in Grover
c	Rprime is the ratio of Rac/Rdc solved for from x
c
	real*4 xx(95),Rprime(95)
	data xx/0.,.5,.6,.7,.8,.9,1.,1.1,1.2,1.3,1.4,1.5,1.6,1.7
     1,1.8,1.9,2.0,2.2,2.4,2.6,2.8,3.0,3.2,3.4,3.6,3.8,4.0,4.2
     1,4.4,4.6,4.8,5.0,5.2,5.4,5.6,5.8,6.0,6.2,6.4,6.6,6.8,7.0
     1,7.2,7.4,7.6,7.8,8.0,8.2,8.4,8.6,8.8,9.0,9.2,9.4,9.6,9.8
     1,10.0,10.5,11.,1.5,12.,12.5,13.,13.5,14.,14.5,15.,16.,17.
     1,18.,19.,20.,21.,22.,23.,24.,25.,26.,28.,30.,32.,34.,36.
     1,38.,40.,42.,44.,46.,48.,50.,60.,70.,80.,90.,100./
	data rprime/1.0,1.0003,1.0007,1.0012,1.0021,1.0034,1.005
     1,1.008,1.011,1.015,1.02,1.026,1.033,1.042,1.052,1.064,1.078
     1,1.111,1.152,1.201,1.256,1.318,1.385,1.456,1.529,1.603,1.678
     1,1.752,1.826,1.899,1.971,2.043,2.114,2.184,2.254,2.324,2.394
     1,2.463,2.533,2.603,2.673,2.743,2.813,2.884,2.954,3.024,3.094
     1,3.165,3.235,3.306,3.376,3.446,3.517,3.587,3.658,3.728,3.799
     1,3.975,4.151,4.327,4.504,4.680,4.856,5.033,5.209,5.386,5.562
     1,5.915,6.268,6.621,6.974,7.328,7.681,8.034,8.387,8.741,9.094
     1,9.447,10.15,10.86,11.57,12.27,12.98,13.69,14.4,15.1,15.81
     1,16.52,17.22,17.93,21.47,25.,28.54,32.07,35.61/
	wdcm=wdiam*2.54
	radius=wdcm/2.0
	x=0.2142*radius*sqrt(Freq)
	if (x.le.0.0) goto 100
	if (x.gt.100.) goto 200
c
c	interpolate between two points
c
c	do 10 j=1,95
c		if (x.eq.xx(j)) then
c		Rac=Rdc*Rprime(j)
c	endif
c	return
c   10   continue
	do 15 j=1,94
	if (x.ge.xx(j).and.x.lt.xx(j+1)) then
	fractn=(x-xx(j))/(xx(j+1)-xx(j))
	r=rprime(j)*(1.0-fractn)+rprime(j+1)*(fractn)
	Rac=Rdc*r
c	write (2,11) x,freq,wdiam,r,rdc,rac
c   11   format(' x,freq,wdiam,r,rdc,rac',6f10.4)
	return
	endif
   15   continue
c
c	no skin effect
c
  100   continue
	Rac=Rdc
	return
c
c	very high skin effect - use approximation above top of table
c
  200   continue
	r=x/(2.0*sqrt(2.0))
	Rac=Rdc*r
	return
	end
c	
c	
c	Lundin single layer circular coil inductance formula
c	Proc. IEEE V73 No 9 Sep 1985, 1428-1429.
c	
c	
	subroutine lundin(result,diam,height,turns)
c	units in meters
	pi=3.141592653589
	a=diam/2.*0.0254
	b=height*0.0254
	u=4.0*pi*1.e-7
	if ((2.0*a).le.b) then
		result=u*turns**2*pi*a**2/b*(f1((4.0*a**2)/(b**2))
     1-(8.0*a)/(3.0*pi*b))
	else
		result=u*turns**2*a*((alog(8.0*a/b)-0.5)
     1*f1((b**2)/(4.0*a**2))+f2((b**2)/(4.0*a**2)))
	endif
	result=result*1.e6
	return
	end
c	
	function f1(x)
	if (x.lt.0.0.or.x.gt.1.0) write (*,10) x
   10   format(' Illegal value for x:',f16.8)
	f1=(1.0+0.383901*x+0.017108*x*x)/(1.0+0.258952*x)
	return
	end
c	
	function f2(x)
	if (x.lt.0.0.or.x.gt.1.0) write (*,10) x
   10   format(' Illegal value for x:',f16.8)
	f2=(0.093842*x+0.002029*x*x-0.000801*x*x*x)
	return
	end
