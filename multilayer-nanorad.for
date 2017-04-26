	subroutine multilayer_nanorad(w,kx,TeTE,TeTM)

*----Declaration of variables
*    nm=maximum number of layers
	parameter (nm=100)
*    nl=number of layers-1; nll is the layer where the field is observed; ss is the emitting layer
	integer nl,nll,ss
*    A, B, C, and D are the T-matrix coefficients in each layer (for both polarizations)
	complex*16 ATE(0:nm),BTE(0:nm),CTE(0:nm),DTE(0:nm)
	complex*16 ATM(0:nm),BTM(0:nm),CTM(0:nm),DTM(0:nm)
*    rfTE,and rfTM are the Fresnel reflections coefficients (TE and TM polarizations) at each interface
	complex*16 rfTE(0:nm-1,0:nm-1),rfTM(0:nm-1,0:nm-1)
*    tfTE,and tfTM are the Fresnel transmission coefficients (TE and TM polarizations) at each interface
	complex*16 tfTE(0:nm-1,0:nm-1),tfTM(0:nm-1,0:nm-1)
*    kz(0:nm-1) is the z-component of the wavevector in each layer. kv(0:nm-1) is the wavector in a layer. 
	complex*16 kz(0:nm-1),kv(0:nm-1)
*    diel(0:nm-1) is the dielectric constant in each layer
	complex*16 diel(0:nm-1)
*    dt(0,nm) is the thickness of a layer; zint is the location of an interface (see figure in dissertation)  
	double precision dt(0:nm),zint(0:nm)
	double precision kx,k0,w,c0,pi
*    im is the complex constant
	complex*16 im
*    Transmitted energy in TM and TE polarizations
	double precision TeTE,TeTM
*    gEij and gHij are the Weyl components of the Green's tensor (i.e., DGF) - We define here only
*    the components that are useful for calculation of the Poynting vector in the z-direction
	complex*16 gExx,gExz,gEyy,gHxy,gHyx,gHyz ,gEzx,gEzz,gHzy
*    funcTE and funcTM are used to calculate the integral over the volume of the emitter
	double precision funcTE(3),funcTM(3)
*    zs is the location of the source point; zsa and zsb are used for the integration over the
*    volume of the emitter (Simpson's method)
	double precision zs,zsa,zsb
*    N_sl is the number of layers in the spatial discretization(sub-layers); del_z is the volume of the sublayers
	integer N_sl
	double precision del_z
*    Sm is the scattering matrix, expjm is the positive exponential term, Qc is the z-component of the wavevector,
*    and Um is the matrix U defined as the product of the propagation matrix and impendance of an interface
	complex*16 Um(0:nm,2,2),Sm(0:nm,0:nm,2,2),expjm(0:nm),Qc(0:nm)
*    Forward and backward propagating source in the emitting layer
	complex*16 S_for,S_bac  

*    Variables to calculate dielectric constants of metals and polar crystals
	double precision einf,wlo,wto,gam
	complex*16 num_diel,den_diel
	double precision wp,re_diel,im_diel
*    Reflectivity of the film: for analytical expressions of Drevillon
	complex*16 RTE,RTM
*    Film reflection coefficients for analytical expression of LDOS within the cavity formed by 2 films
	complex*16 RTM1,RTM3 ! films 1 and 3
	complex*16 RTE1,RTE3 ! films 1 and 3
*    Variables for different tests
	complex*16 ATE2,TTM3,termTM
	
**********************INPUT OF THE OPTICAL PROPERTIES OF EACH LAYERS AND THEIR THICKNESSES*****************************
	pi=4.d0*atan(1.d0)
	c0=2.998d8 ! Speed of light in vacuum
	k0=w/c0 ! Magnitude of wavevector in vaccum
	im=dcmplx(0.d0,1.d0) ! complex constant

*----Number of layers
	nl=5 ! Number of layers - 1 (since half-space is denoted 0)
	nll=3 ! Layer where we want to calculate the flux
	ss=1 ! emitting layer

*----Parameters for spatial integration
	N_sl=10
	del_z=1.d-9

*----Optical properties in each layer
*    Medium 0 is vacuum
	diel(0)=dcmplx(1.d0,0.d0)
*    Medium 1 is SiC
	einf=6.7d0
	wlo=182.53d12
	wto=149.37d12
	gam=89.66d10
	num_diel=dcmplx(((w**2)-(wlo**2)),gam*w)
	den_diel=dcmplx(((w**2)-(wto**2)),gam*w)
	diel(1)=einf*(num_diel/den_diel)
*    Media 2 and 3 are vacuum
	diel(2)=dcmplx(1.d0,0.d0)
	diel(3)=diel(2)
	diel(4)=diel(1)
	diel(5)=diel(2)

*----Thickness of each layer. Note that the thicknesses of the two half-spaces (0 and nl) should be zero
*    zint is the z-location of an interface. For example, layer 1 is bounded by zint(1) and zint(2)
	dt(0)=0.d0
	dt(1)=10.d-9
	dt(2)=10.d-9
	dt(3)=0.d-9
	dt(4)=100.d-9
	dt(5)=0.d0

*----Coordinates of the interface of a layer. For layer j, the left interface is zint(j) and the right interface
*    is zint(j+1)
	zint(1)=0.d0
	zint(2)=zint(1)+dt(1)
	zint(3)=zint(2)+dt(2)
	zint(4)=zint(3)+dt(3)
	zint(5)=zint(4)+dt(4)

********************************************END OF INPUT***************************************************************
*----Wavevector in each layer
	do i=0,nl
		kv(i)=sqrt(diel(i)*(k0**2))
	enddo

*----z-component of the wavevector in each layer
	do i=0,nl
		kz(i)=sqrt((kv(i)**2)-(kx**2))
	enddo

*----Fresnel's reflection and transmission coefficients at each interface
	do i=0,nl
		do j=0,nl
			rfTE(i,j)=(kz(i)-kz(j))/(kz(i)+kz(j))
			rfTM(i,j)=((kz(i)*diel(j))-diel(i)*kz(j))
     &				/((kz(i)*diel(j))+kz(j)*diel(i))

			tfTE(i,j)=2.d0*kz(i)/(kz(i)+kz(j))
			tfTM(i,j)=2.d0*sqrt(diel(i))*sqrt(diel(j))*kz(i)
     &/(diel(i)*kz(j)+diel(j)*kz(i))
		enddo
	enddo


*----Loop for the volume integration
	TeTE=0.d0
	TeTM=0.d0

	do 1 i=1,N_sl
		zsa=((i-1)*del_z)+zint(ss)
		zsb=(i*del_z)+zint(ss)

		do 2 k=1,3
			if (k.eq.1) then
				zs=zsa
			elseif (k.eq.2) then
				zs=(zsa+zsb)/2.d0
			elseif (k.eq.3) then
				zs=zsb
			endif		

*----Fields transmission and reflection coefficients in TE-MODES (modified T-matrix approach)------------------*
*    This gives the coefficients A, B, C, and D in each layer.
*    The solution is based on the modified T-matrix approach, also called "Scattering matrix" method.

!    Calculation of S-matrix relative to layer 0. In that case, the matrix S(0,0)=I, where I is the identity matrix.
	Sm(0,0,1,1)=dcmplx(1.d0,0.d0)
	Sm(0,0,1,2)=dcmplx(0.d0,0.d0)
	Sm(0,0,2,1)=dcmplx(0.d0,0.d0)
	Sm(0,0,2,2)=dcmplx(1.d0,0.d0)
	expjm(0)=dcmplx(1.d0,0.d0)
	Qc(0)=kz(0)

	do j=1,ss ! Calculation of S-matrix S(0,1) to S(0,ss)
		expjm(j)=cdexp(im*kz(j)*dt(j))
		Qc(j)=kz(j)
		Um(j,2,1)=0.5d0*expjm(j-1)*(2.d0*rfTE(j-1,j)/tfTE(j-1,j))
		Um(j,2,2)=0.5d0*expjm(j-1)*(2.d0/tfTE(j-1,j))
		Sm(0,j,1,1)=Sm(0,j-1,1,1)*(2.d0*expjm(j-1)
     &	*(1.d0/((2.d0/tfTE(j-1,j))-(Sm(0,j-1,1,2)
     &	*expjm(j-1)*expjm(j-1)*(2.d0*rfTE(j-1,j)/tfTE(j-1,j))))))
		Sm(0,j,1,2)=((Sm(0,j-1,1,2)*expjm(j-1)
     &	*expjm(j-1)*(2.d0/tfTE(j-1,j)))-(2.d0*rfTE(j-1,j)
     &    /tfTE(j-1,j)))
     &	/((2.d0/tfTE(j-1,j))-(Sm(0,j-1,1,2)*expjm(j-1)*expjm(j-1)
     &	*(2.d0*rfTE(j-1,j)/tfTE(j-1,j))))
		Sm(0,j,2,1)=Sm(0,j-1,2,2)*Um(j,2,1)*Sm(0,j,1,1)+Sm(0,j-1,2,1)
		Sm(0,j,2,2)=Sm(0,j-1,2,2)*Um(j,2,1)*Sm(0,j,1,2)
     &	+(Sm(0,j-1,2,2)*Um(j,2,2))
	enddo

!    Calculation of S-matrix relative to layer ss. In that case, the matrix S(ss,ss)=I, where I is the identity matrix.
	Sm(ss,ss,1,1)=dcmplx(1.d0,0.d0)
	Sm(ss,ss,1,2)=dcmplx(0.d0,0.d0)
	Sm(ss,ss,2,1)=dcmplx(0.d0,0.d0)
	Sm(ss,ss,2,2)=dcmplx(1.d0,0.d0)
	expjm(ss)=cdexp(im*kz(ss)*dt(ss))
	Qc(ss)=kz(ss)

	do j=ss+1,nl ! Calculation of S-matrix S(ss,ss+1) to S(ss,nl)
		expjm(j)=cdexp(im*kz(j)*dt(j))
		Qc(j)=kz(j)
		Um(j,2,1)=0.5d0*expjm(j-1)*(2.d0*rfTE(j-1,j)/tfTE(j-1,j))
		Um(j,2,2)=0.5d0*expjm(j-1)*(2.d0/tfTE(j-1,j))
		Sm(ss,j,1,1)=Sm(ss,j-1,1,1)*(2.d0*expjm(j-1)
     &	*(1.d0/((2.d0/tfTE(j-1,j))-(Sm(ss,j-1,1,2)
     &	*expjm(j-1)*expjm(j-1)*(2.d0*rfTE(j-1,j)/tfTE(j-1,j))))))
		Sm(ss,j,1,2)=((Sm(ss,j-1,1,2)*expjm(j-1)
     &	*expjm(j-1)*(2.d0/tfTE(j-1,j)))-(2.d0*rfTE(j-1,j)
     &    /tfTE(j-1,j)))
     &	/((2.d0/tfTE(j-1,j))-(Sm(ss,j-1,1,2)*expjm(j-1)*expjm(j-1)
     &	*(2.d0*rfTE(j-1,j)/tfTE(j-1,j))))
		Sm(ss,j,2,1)=Sm(ss,j-1,2,2)*Um(j,2,1)*Sm(ss,j,1,1)
     &    +Sm(ss,j-1,2,1)
		Sm(ss,j,2,2)=Sm(ss,j-1,2,2)*Um(j,2,1)*Sm(ss,j,1,2)
     &	+(Sm(ss,j-1,2,2)*Um(j,2,2))
	enddo

!----Calculation of ATE & BTE

!    Calculation of the emission term on the boundary zint(ss) (i.e., left boundary).
!    For A & B, the source is propagating in the forward direction (i.e., z-positive direction).
!    This source is included in the relation between layer ss and nl (A(ss)+source).
	S_for=1.d0*cdexp(im*kz(ss)*(zint(ss)-zs))

!    Values of some coefficients
	ATE(0)=dcmplx(0.d0,0.d0)
	BTE(nl)=dcmplx(0.d0,0.d0)
	BTE(ss)=(Sm(ss,nl,2,1)*S_for)
     &/(1.d0-(Sm(ss,nl,2,1)*Sm(0,ss,1,2)))
	BTE(0)=Sm(0,ss,2,2)*BTE(ss)
	ATE(ss)=Sm(0,ss,1,2)*BTE(ss)
	ATE(nl)=Sm(ss,nl,1,1)*(ATE(ss)+S_for)

!    Recursive scheme to find all coefficients in layer 1 to ss-1
	do j=1,ss-1
		BTE(j)=BTE(0)/Sm(0,j,2,2)
		ATE(j)=Sm(0,j,1,2)*BTE(j)
	enddo

!    Recursive scheme to find all coefficients in layer ss+1 to nl
	do j=ss+1,nl-1
		BTE(j)=(BTE(ss)-(Sm(ss,j,2,1)*(ATE(ss)+S_for)))/Sm(ss,j,2,2) 
		ATE(j)=(Sm(ss,j,1,1)*(ATE(ss)+S_for))+(Sm(ss,j,1,2)*BTE(j))
	enddo

!----Calculation of CTE & DTE

!    Calculation of the emission term on the boundary zint(ss) (i.e., left boundary).
!    For C & D, the source is propagating in the backward direction (i.e., z-negative direction).
!    This source is included in the relation between layer 0 and ss (D(ss)+source).
	S_bac=1.d0*cdexp(im*kz(ss)*(zs-zint(ss)))

!    Values of some coefficients
	CTE(0)=dcmplx(0.d0,0.d0)
	DTE(nl)=dcmplx(0.d0,0.d0)
	CTE(ss)=(Sm(0,ss,1,2)*S_bac)/(1.d0-(Sm(0,ss,1,2)*Sm(ss,nl,2,1)))
	CTE(nl)=Sm(ss,nl,1,1)*CTE(ss)
	DTE(ss)=Sm(ss,nl,2,1)*CTE(ss)
	DTE(0)=Sm(0,ss,2,2)*(DTE(ss)+S_bac)

!    Recursive scheme to find all coefficients in layer 1 to ss-1
	do j=1,ss-1
		DTE(j)=DTE(0)/Sm(0,j,2,2)
		CTE(j)=Sm(0,j,1,2)*DTE(j)
	enddo

!    Recursive scheme to find all coefficients in layer ss+1 to nl
	do j=ss+1,nl-1
		DTE(j)=(DTE(ss)-(Sm(ss,j,2,1)*CTE(ss)))/Sm(ss,j,2,2)
		CTE(j)=(Sm(ss,j,1,1)*CTE(ss))+(Sm(ss,j,1,2)*DTE(j))
	enddo

*---------------------------------------End of calculations for TE-MODES------------------------------*


*----Fields transmission and reflection coefficients in TM-MODES (modified T-matrix approach)------------------*
*    This gives the coefficients A, B, C, and D in each layer.
*    The solution is based on the modified T-matrix approach, also called "Scattering matrix" method.

!    Calculation of S-matrix relative to layer 0. In that case, the matrix S(0,0)=I, where I is the identity matrix.
	Sm(0,0,1,1)=dcmplx(1.d0,0.d0)
	Sm(0,0,1,2)=dcmplx(0.d0,0.d0)
	Sm(0,0,2,1)=dcmplx(0.d0,0.d0)
	Sm(0,0,2,2)=dcmplx(1.d0,0.d0)
	expjm(0)=dcmplx(1.d0,0.d0)
	Qc(0)=kz(0)

	do j=1,ss ! Calculation of S-matrix S(0,1) to S(0,ss)
		expjm(j)=cdexp(im*kz(j)*dt(j))
		Qc(j)=kz(j)
		Um(j,2,1)=0.5d0*expjm(j-1)*(2.d0*rfTM(j-1,j)/tfTM(j-1,j))
		Um(j,2,2)=0.5d0*expjm(j-1)*(2.d0/tfTM(j-1,j))
		Sm(0,j,1,1)=Sm(0,j-1,1,1)*(2.d0*expjm(j-1)
     &	*(1.d0/((2.d0/tfTM(j-1,j))-(Sm(0,j-1,1,2)
     &	*expjm(j-1)*expjm(j-1)*(2.d0*rfTM(j-1,j)/tfTM(j-1,j))))))
		Sm(0,j,1,2)=((Sm(0,j-1,1,2)*expjm(j-1)
     &	*expjm(j-1)*(2.d0/tfTM(j-1,j)))-(2.d0*rfTM(j-1,j)
     &    /tfTM(j-1,j)))
     &	/((2.d0/tfTM(j-1,j))-(Sm(0,j-1,1,2)*expjm(j-1)*expjm(j-1)
     &	*(2.d0*rfTM(j-1,j)/tfTM(j-1,j))))
		Sm(0,j,2,1)=Sm(0,j-1,2,2)*Um(j,2,1)*Sm(0,j,1,1)+Sm(0,j-1,2,1)
		Sm(0,j,2,2)=Sm(0,j-1,2,2)*Um(j,2,1)*Sm(0,j,1,2)
     &	+(Sm(0,j-1,2,2)*Um(j,2,2))
	enddo

!    Calculation of S-matrix relative to layer ss. In that case, the matrix S(ss,ss)=I, where I is the identity matrix.
	Sm(ss,ss,1,1)=dcmplx(1.d0,0.d0)
	Sm(ss,ss,1,2)=dcmplx(0.d0,0.d0)
	Sm(ss,ss,2,1)=dcmplx(0.d0,0.d0)
	Sm(ss,ss,2,2)=dcmplx(1.d0,0.d0)
	expjm(ss)=cdexp(im*kz(ss)*dt(ss))
	Qc(ss)=kz(ss)

	do j=ss+1,nl ! Calculation of S-matrix S(ss,ss+1) to S(ss,nl)
		expjm(j)=cdexp(im*kz(j)*dt(j))
		Qc(j)=kz(j)
		Um(j,2,1)=0.5d0*expjm(j-1)*(2.d0*rfTM(j-1,j)/tfTM(j-1,j))
		Um(j,2,2)=0.5d0*expjm(j-1)*(2.d0/tfTM(j-1,j))
		Sm(ss,j,1,1)=Sm(ss,j-1,1,1)*(2.d0*expjm(j-1)
     &	*(1.d0/((2.d0/tfTM(j-1,j))-(Sm(ss,j-1,1,2)
     &	*expjm(j-1)*expjm(j-1)*(2.d0*rfTM(j-1,j)/tfTM(j-1,j))))))
		Sm(ss,j,1,2)=((Sm(ss,j-1,1,2)*expjm(j-1)
     &	*expjm(j-1)*(2.d0/tfTM(j-1,j)))-(2.d0*rfTM(j-1,j)/tfTM(j-1,j)))
     &	/((2.d0/tfTM(j-1,j))-(Sm(ss,j-1,1,2)*expjm(j-1)*expjm(j-1)
     &	*(2.d0*rfTM(j-1,j)/tfTM(j-1,j))))
		Sm(ss,j,2,1)=Sm(ss,j-1,2,2)*Um(j,2,1)*Sm(ss,j,1,1)+Sm(ss,j-1,2,1)
		Sm(ss,j,2,2)=Sm(ss,j-1,2,2)*Um(j,2,1)*Sm(ss,j,1,2)
     &	+(Sm(ss,j-1,2,2)*Um(j,2,2))
	enddo

!----Calculation of ATM & BTM

!    Calculation of the emission term on the boundary zint(ss) (i.e., left boundary).
!    For A & B, the source is propagating in the forward direction (i.e., z-positive direction).
!    This source is included in the relation between layer ss and nl (A(ss)+source).
	S_for=1.d0*cdexp(im*kz(ss)*(zint(ss)-zs))

!    Values of some coefficients
	ATM(0)=dcmplx(0.d0,0.d0)
	BTM(nl)=dcmplx(0.d0,0.d0)
	BTM(ss)=(Sm(ss,nl,2,1)*S_for)
     &/(1.d0-(Sm(ss,nl,2,1)*Sm(0,ss,1,2)))
	BTM(0)=Sm(0,ss,2,2)*BTM(ss)
	ATM(ss)=Sm(0,ss,1,2)*BTM(ss)
	ATM(nl)=Sm(ss,nl,1,1)*(ATM(ss)+S_for)

!    Recursive scheme to find all coefficients in layer 1 to ss-1
	do j=1,ss-1
		BTM(j)=BTM(0)/Sm(0,j,2,2)
		ATM(j)=Sm(0,j,1,2)*BTM(j)
	enddo

!    Recursive scheme to find all coefficients in layer ss+1 to nl
	do j=ss+1,nl-1
		BTM(j)=(BTM(ss)-(Sm(ss,j,2,1)*(ATM(ss)+S_for)))/Sm(ss,j,2,2) 
		ATM(j)=(Sm(ss,j,1,1)*(ATM(ss)+S_for))+(Sm(ss,j,1,2)*BTM(j))
	enddo

!----Calculation of CTM & DTM

!    Calculation of the emission term on the boundary zint(ss) (i.e., left boundary).
!    For C & D, the source is propagating in the backward direction (i.e., z-negative direction).
!    This source is included in the relation between layer 0 and ss (D(ss)+source).
	S_bac=1.d0*cdexp(im*kz(ss)*(zs-zint(ss)))

!    Values of some coefficients
	CTM(0)=dcmplx(0.d0,0.d0)
	DTM(nl)=dcmplx(0.d0,0.d0)
	CTM(ss)=(Sm(0,ss,1,2)*S_bac)/(1.d0-(Sm(0,ss,1,2)*Sm(ss,nl,2,1)))
	CTM(nl)=Sm(ss,nl,1,1)*CTM(ss)
	DTM(ss)=Sm(ss,nl,2,1)*CTM(ss)
	DTM(0)=Sm(0,ss,2,2)*(DTM(ss)+S_bac)

!    Recursive scheme to find all coefficients in layer 1 to ss-1
	do j=1,ss-1
		DTM(j)=DTM(0)/Sm(0,j,2,2)
		CTM(j)=Sm(0,j,1,2)*DTM(j)
	enddo

!    Recursive scheme to find all coefficients in layer ss+1 to nl
	do j=ss+1,nl-1
		DTM(j)=(DTM(ss)-(Sm(ss,j,2,1)*CTM(ss)))/Sm(ss,j,2,2)
		CTM(j)=(Sm(ss,j,1,1)*CTM(ss))+(Sm(ss,j,1,2)*DTM(j))
	enddo

*---------------------------------------End of calculations for TM-MODES------------------------------*

*----Calculation of Weyl components of DGF
	gExx=(im*kz(nll))/(2.d0*kv(ss)*kv(nll))
     &	*(ATM(nll)-BTM(nll)-CTM(nll)+DTM(nll))
	gExz=(im*kz(nll)*kx)/(2.d0*kz(ss)*kv(ss)*kv(nll))
     &	*(-ATM(nll)+BTM(nll)-CTM(nll)+DTM(nll))
	gEyy=im/(2.d0*kz(ss))
     &	*(ATE(nll)+BTE(nll)+CTE(nll)+DTE(nll))
	gEzx=(im*kx)/(2.d0*kv(ss)*kv(nll))
     &	*(-ATM(nll)-BTM(nll)+CTM(nll)+DTM(nll))
	gEzz=(im*kx*kx)/(2.d0*kz(ss)*kv(ss)*kv(nll))
     &	*(ATM(nll)+BTM(nll)+CTM(nll)+DTM(nll))
	gHxy=(kz(nll))/(2.d0*kz(ss))
     &	*(ATE(nll)-BTE(nll)+CTE(nll)-DTE(nll))
	gHyx=(kv(nll))/(2.d0*kv(ss))
     &	*(-ATM(nll)-BTM(nll)+CTM(nll)+DTM(nll))
	gHyz=(kv(nll)*kx)/(2.d0*kv(ss)*kz(ss))
     &	*(ATM(nll)+BTM(nll)+CTM(nll)+DTM(nll))
	gHzy=(kv(nll)*kx)/(2.d0*kv(nll)*kz(ss))
     &	*(-ATE(nll)-BTE(nll)-CTE(nll)-DTE(nll))

*----Formulae of ED taken from Zhang's book and paper (consistent with definition
*    given in book chapter of Joulain)
*    u=4*((e0/4)<E**2>+(u0/4)<H**2>)
	funcTE(k)=
     &(((w**3.d0)*dimag(diel(ss)))/(2.d0*(pi**2.d0)*(c0**4.d0)))
     &*((gEyy*dconjg(gEyy))
     &+((1.d0/(k0**2.d0))*((gHxy*dconjg(gHxy))+(gHzy*dconjg(gHzy)))))

	funcTM(k)=
     &(((w**3.d0)*dimag(diel(ss)))/(2.d0*(pi**2.d0)*(c0**4.d0)))
     &*((gExx*dconjg(gExx))+(gExz*dconjg(gExz))
     &+(gEzx*dconjg(gEzx))+(gEzz*dconjg(gEzz))
     &+((1.d0/(k0**2.d0))*((gHyx*dconjg(gHyx))+(gHyz*dconjg(gHyz)))))

2	continue

*----Calculation of energy transmitted: Integration by Simpson's method
	TeTE=TeTE
     &	+((zsb-zsa)/6.d0)*(funcTE(1)+4.d0*funcTE(2)+funcTE(3))
	TeTM=TeTM
     &	+((zsb-zsa)/6.d0)*(funcTM(1)+4.d0*funcTM(2)+funcTM(3))

1	continue
	write(6,*) 'k_rho=',kx	
	write(6,*) '1-Num-TE=',TeTE	

*----Analytical expression of Drevillon et al. for a single emitting film
!MF	RTE=(rfTE(0,1)+rfTE(1,2)*cdexp(2.d0*im*kz(1)*(zint(2)-zint(1))))
!MF     &/(1.d0+rfTE(0,1)*rfTE(1,2)*cdexp(2.d0*im*kz(1)*(zint(2)-zint(1))))
!MF	RTM=(rfTM(0,1)+rfTM(1,2)*cdexp(2.d0*im*kz(1)*(zint(2)-zint(1))))
!MF     &/(1.d0+rfTM(0,1)*rfTM(1,2)*cdexp(2.d0*im*kz(1)*(zint(2)-zint(1))))

!MF	TeTE=((kx**2.d0)*dimag(RTE)*dexp(-2.d0*dimag(kz(2))*dt(2)))
!MF     &/(2.d0*(pi**2.d0)*w*cdabs(kz(2)))
!MF	TeTM=((kx**2.d0)*dimag(RTM)*dexp(-2.d0*dimag(kz(2))*dt(2)))
!MF     &/(2.d0*(pi**2.d0)*w*cdabs(kz(2)))




!----Analytical formulation when looking at distance z above film 1 in cavity (see p. 47 of notes)
!    The equation below is correct, and give the exact same result than when solving the problem
!    with the pure numerical approach.
	RTE1=(rfTE(0,1)+rfTE(1,2)*cdexp(2.d0*im*kz(1)*dt(1)))
     &/(1.d0+rfTE(0,1)*rfTE(1,2)*cdexp(2.d0*im*kz(1)*dt(1)))

	RTE3=(rfTE(3,4)+rfTE(4,5)*cdexp(2.d0*im*kz(4)*dt(4)))
     &/(1.d0+rfTE(3,4)*rfTE(4,5)*cdexp(2.d0*im*kz(4)*dt(4)))

	termTE=(cdabs(1.d0+RTE3*dexp(-2.d0*dimag(kz(2))*dt(3)))**2.d0)
     &-(2.d0*(cdabs(kz(2))**2.d0)
     &*dreal(RTE3)*(dexp(-2.d0*dimag(kz(2))*dt(3))))/(kx**2.d0)

	TeTE=(((1.d0/(2.d0*(pi**2.d0)*w))*(kx**2.d0))*(1.d0/cdabs(kz(2)))
     &*dimag(RTE1)*termTE*dexp(-2.d0*dimag(kz(2))*dt(2)))
     &/(cdabs(1.d0-RTE1*RTE3*cdexp(2.d0*im*kz(2)*(dt(2)+dt(3))))**2.d0)

	 
	write(6,*) '2-Ana-TE=',TeTE

	!stop

	end subroutine




