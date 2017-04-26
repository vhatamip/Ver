*****The objective of this program is to calculate the energy density (ED) and local density 
*****of electromagnetic states (LDOS) above a film emitter.

*****We assume plane waves of the form exp[i(kx-wt)]. Therefore, exp(ikx) is a wave traveling in the 
*****z-positive direction, and exp(-ikx) is a wave traveling in the z-negative direction.

*****Mathieu Francoeur
*****Radiative Transfer Laboratory
*****Department of Mechanical Engineering
*****University of Kentucky, Lexington, KY 40506

*****Last update: November 2008

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    VALIDATION: The code has been validated using the following paper:									 !
!    A. Drevillon, P. Ben-Abdallah, K. Joulain, and G. Domingues, "Local density of states of			     !
!    non-propagative electromagnetic field at the surface polaritons frequency of thin metallic films",	 !
!    Sumitted to JQSRT, 2008. Results of Fig. 2.															 !
!																										 !			
!    Other validation: Biehs et al., "Thermal radiation and near-field energy density of thin metallic     !
!    films", The European Physical Journal B 55, 237-251, 2007. Results of Fig. 8.                         !
!    We look at the LDOS in TM-polarization of Bismuth film (zobs=10 nm above the film in vacuum)			 !
!    No able to validate these results: the Optical properties of Bi seem to be the problem		         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	program bulk LDOS ED_main
*----Declaration of variables
*    jm is a variable used for integration with Simpson's method; mm is the maximum number of frequencies
	parameter (jm=3,mm=200000)
*    useful constants
	double precision c0,pi,h,kb
*    w is the frequency; t0 is the temperature of the emitter
	double precision w,t0
*    theta1 is the mean energy of a Planck oscillator in thermal equilibrium
	double precision theta1 
*    ldos and ed are the monochromatic local density of electromagnetic states and
*    energy density (prop=propagating, evan=evanescent, tot=prop+evan)
	double precision ldos_prop(0:mm),ldos_evan(0:mm),ldos_tot(0:mm)
	double precision ed_prop(0:mm),ed_evan(0:mm),ed_tot(0:mm)
	double precision ed_propTM(0:mm),ed_propTE(0:mm),ed_evanTM(0:mm),
     &				 ed_evanTE(0:mm)
	double precision ldos_propTM(0:mm),ldos_propTE(0:mm),
     &				 ldos_evanTM(0:mm),ldos_evanTE(0:mm)
*    kx is the x-component of the wavevector (same in all layers); kxa and kxb are defined for calculation of the integrals
*    using Simpson's method	
	double precision kx,kxa,kxb,k_rho	
*    dfunct, dfunctev are used to perform the integration over kx with Simpson's method
	double precision dfunct(jm),dfunctev(jm),dfunctTM(jm),dfunctTE(jm)
*    General transmission coefficients of energy (TE and TM modes) calculated with the subroutine (T-matrix method)
	double precision TeTE,TeTM
*    weight is the interval of frequencies for which the integration is performed
	double precision weight
*    sum are used to perform integration
	double precision sum1,sum2,sum3
	double precision sumTM1,sumTM2,sumTM3,sumTE1,sumTE2,sumTE3
*    variables for convergence on kx
	double precision conv1,p_sum

*******************************************************INPUT BY THE USER*****************************************
	t0=300.d0 ! temperature of the bulk
	c0=2.998d8 ! speed of light in vacuum
	pi=4.d0*atan(1.d0)
	h=1.055d-34 ! circular Planck's constant
	kb=1.381d-23 ! Boltzmann constant

*----Initialization of the sum
	sum1=0.d0
	sum2=0.d0
	sum3=0.d0
	sumTM1=0.d0
	sumTM2=0.d0
	sumTM3=0.d0
	sumTE1=0.d0
	sumTE2=0.d0
	sumTE3=0.d0

*----Initialization of the coefficients for transmission of energy
	TeTE=0.d0
	TeTM=0.d0

*----Open files to print monochromatic fluxes
!	open (unit=100,file='Mono_ed_prop.txt',status='unknown')
!	open (unit=200,file='Mono_ed_evan.txt',status='unknown')
!	open (unit=300,file='Mono_ed_tot.txt',status='unknown')
!	open (unit=400,file='Mono_ldos_prop.txt',status='unknown')
!	open (unit=500,file='Mono_ldos_evan.txt',status='unknown')
!	open (unit=600,file='Mono_ldos_tot.txt',status='unknown')

!	open (unit=700,file='Mono_ed_propTM.txt',status='unknown')
!	open (unit=800,file='Mono_ed_propTE.txt',status='unknown')
!	open (unit=900,file='Mono_ed_evanTM.txt',status='unknown')
!	open (unit=1000,file='Mono_ed_evanTE.txt',status='unknown')

!	open (unit=1100,file='Mono_ldos_propTM.txt',status='unknown')
!	open (unit=1200,file='Mono_ldos_propTE.txt',status='unknown')
	open (unit=1300,file='Mono_ldos_evanTM.txt',status='unknown')
!	open (unit=1400,file='Mono_ldos_evanTE.txt',status='unknown')

*----Loop for all frequencies.
cccccStudy of surface phonon-polaritons coupling in SiC
c	weight=0.5d12
c	do 1 m=1,501
c	w=50.d12+((m-1)*weight)

	weight=0.1d12
	do 1 m=1,401
c	weight=0.5d12
c	do 1 m=1,81
	w=150.d12+((m-1)*weight)
      k_rho=w/c0
cccccFrequencies for validation of Drevillon et al. (Al film)
c	weight=4.d12
c	do 1 m=1,120
c	w=17.47d15*1.d-2*m

	theta1=h*w/(dexp((h*w)/(kb*t0))-1.d0) ! mean energy of Planck's oscillator

*******************************PROPAGATING COMPONENT OF THE RADIATIVE HEAT FLUX*******************************
*----Loop for all wavevectors (integration over kx)
	do 2 j=1,1999
!c	do 2 j=1,1
		kxa=(j-1.d0)*((w/c0)/2000.d0)
		kxb=j*((w/c0)/2000.d0)
          
*    Loop to perform the integration by Simpson's method
		do 3 l=1,3
			if (l.eq.1) then
				kx=kxa
			elseif (l.eq.2) then
				kx=(kxa+kxb)/2.d0
			elseif (l.eq.3) then
				kx=kxb
			endif

*----Calculation of the general transmission coefficient. Output = TeTE and TeTM 
			call multilayer_nanorad(w,kx,TeTE,TeTM) 

*----Integration over all kx (define a function that we evaluate) by Simpson method
			dfunct(l)=(TeTE+TeTM)*kx
			dfunctTM(l)=TeTM*kx
			dfunctTE(l)=TeTE*kx
3		continue

		sum1=sum1+((kxb-kxa)/6.d0)
     &	*(dfunct(1)+4.d0*dfunct(2)+dfunct(3))
		sumTM1=sumTM1+((kxb-kxa)/6.d0)
     &	*(dfunctTM(1)+4.d0*dfunctTM(2)+dfunctTM(3))
		sumTE1=sumTE1+((kxb-kxa)/6.d0)
     &	*(dfunctTE(1)+4.d0*dfunctTE(2)+dfunctTE(3))

2	continue

*******************************END OF PROPAGATING COMPONENT OF THE RADIATIVE HEAT FLUX*******************************


*******************************EVANESCENT COMPONENT OF THE RADIATIVE HEAT FLUX - LOOP 1******************************

*----Loop for all wavevectors (integration over kx)
c	do 4 j=102,600
c		kxa=((j-1)*(w/c0))/(100.d0)
c		kxb=(j*(w/c0))/(100.d0)
c	do 4 j=30002,180000
c		kxa=((j-1)*(w/c0))/(30000.d0)
c		kxb=(j*(w/c0))/(30000.d0)
	do 4 j=40002,240000
		kxa=((j-1)*(w/c0))/(40000.d0)
		kxb=(j*(w/c0))/(40000.d0)

*    Loop to perform the integration by Simpson's method
		do 5 l=1,3
			if (l.eq.1) then
				kx=kxa
			elseif (l.eq.2) then
				kx=(kxa+kxb)/2.d0
			elseif (l.eq.3) then
				kx=kxb
			endif

*----Calculation of the general transmission coefficient of energy. Output = TeTE and TeTM 
			call multilayer_nanorad(w,kx,TeTE,TeTM)

*----Integration over all kx (define a function that we evaluate) by Simpson method

			dfunctev(l)=(TeTE+TeTM)*kx
			dfunctTM(l)=TeTM*kx
			dfunctTE(l)=TeTE*kx
 
5		continue

		sum2=sum2+((kxb-kxa)/6.d0)*(dfunctev(1)+4.d0*dfunctev(2)
     &	+dfunctev(3))
		sumTM2=sumTM2+((kxb-kxa)/6.d0)
     &	*(dfunctTM(1)+4.d0*dfunctTM(2)+dfunctTM(3))
		sumTE2=sumTE2+((kxb-kxa)/6.d0)
     &	*(dfunctTE(1)+4.d0*dfunctTE(2)+dfunctTE(3))

4	continue

*******************************END OF EVANESCENT COMPONENT OF THE RADIATIVE HEAT FLUX - LOOP 1***************************


*******************************EVANESCENT COMPONENT OF THE RADIATIVE HEAT FLUX - LOOP 2**********************************

*----Loop for all wavevectors-second loop (integration over kx)
	j=2
	p_sum=0.d0
	conv1=1.d0

	do 6 while (conv1.gt.1.d-8) ! definition of the convergence criterion
c	do 6 j=2,8000
c		kxa=((6.d0*w)/c0)+((j-1)*(w/c0))/(200.0d0)
c		kxb=((6.d0*w)/c0)+(j*(w/c0))/(200.0d0)
		kxa=((6.d0*w)/c0)+((j-1)*(w/c0))/(1000.0d0)
		kxb=((6.d0*w)/c0)+(j*(w/c0))/(1000.0d0)

		j=j+1

*    Loop to perform the integration by Simpson's method
		do 7 l=1,3
			if (l.eq.1) then
				kx=kxa
			elseif (l.eq.2) then
				kx=(kxa+kxb)/2.d0
			elseif (l.eq.3) then
				kx=kxb
			endif

*----Calculation of the general transmission coefficient. Output = TeTE and TeTM 
			call multilayer_nanorad(w,kx,TeTE,TeTM)

*----Integration over all kx (define a function that we evaluate) by Simpson method
			dfunctev(l)=(TeTE+TeTM)*kx
			dfunctTM(l)=TeTM*kx
			dfunctTE(l)=TeTE*kx
  
7	continue

		sum3=sum3+((kxb-kxa)/6.d0)*(dfunctev(1)+4.d0*dfunctev(2)
     &	+dfunctev(3))
		sumTM3=sumTM3+((kxb-kxa)/6.d0)
     &	*(dfunctTM(1)+4.d0*dfunctTM(2)+dfunctTM(3))
		sumTE3=sumTE3+((kxb-kxa)/6.d0)
     &	*(dfunctTE(1)+4.d0*dfunctTE(2)+dfunctTE(3))

		conv1=(dabs(sum3-p_sum))/sum3
		p_sum=sum3 ! assignation of the last sum to p_sum

6	continue

*******************************END EVANESCENT COMPONENT OF THE RADIATIVE HEAT FLUX - LOOP 2*******************************


*********************************************CALCULATION OF ED AND LDOS***************************************************

*----Calculation of the monochromatic energy density and ldos. Note that LDOS is simply the enrgy density divided
*    by the mean energy of a Planck's oscillator.
	ed_prop(m)=theta1*sum1
	ed_evan(m)=theta1*(sum2+sum3)
	ed_tot(m)=theta1*(sum1+sum2+sum3)

	ldos_prop(m)=sum1
	ldos_evan(m)=sum2+sum3
	ldos_tot(m)=sum1+sum2+sum3

!	write(6,*) 'ldos_evan = ',w,j,ldos_evan(m)
!	write(100,*) w,ed_prop(m)
!	write(200,*) w,ed_evan(m)
!	write(300,*) w,ed_tot(m)
!	write(400,*) w,ldos_prop(m)
!	write(500,*) w/17.47d15,ldos_evan(m)
!	write(500,*) w,ldos_evan(m)
!	write(600,*) w,ldos_tot(m)

	ed_propTM(m)=theta1*sumTM1
	ed_propTE(m)=theta1*sumTE1
	ed_evanTM(m)=theta1*(sumTM2+sumTM3)
	ed_evanTE(m)=theta1*(sumTE2+sumTE3)

	ldos_propTM(m)=sumTM1
	ldos_propTE(m)=sumTE1
	ldos_evanTM(m)=sumTM2+sumTM3
	ldos_evanTE(m)=sumTE2+sumTE3

!	write(700,*) w,ed_propTM(m)
!	write(800,*) w,ed_propTE(m)
!	write(900,*) w,ed_evanTM(m)
!	write(1000,*) w,ed_evanTE(m)

	write(6,*) 'ldos_evan TM = ',w/1.d12,j,ldos_evanTM(m)

!	write(1100,*) w,ldos_propTM(m)
!	write(1200,*) w,ldos_propTE(m)
	write(1300,*) w,ldos_evanTM(m)
!	write(1400,*) w,ldos_evanTE(m)

*----sum1,sum2,sum3 are re-initialized at 0 for calculations at next frequency
	sum1=0.d0
	sum2=0.d0
	sum3=0.d0
	sumTM1=0.d0
	sumTM2=0.d0
	sumTM3=0.d0
	sumTE1=0.d0
	sumTE2=0.d0
	sumTE3=0.d0

1	continue ! End of loop for all frequencies

!	close(100)
!	close(200)
!	close(300)
!	close(400)
!	close(500)
!	close(600)
!	close(700)
!	close(800)
!	close(900)
!	close(1000)
!	close(1100)
!	close(1200)
	close(1300)
!	close(1400)

******************************************END OF CALCULATION OF ED AND LDOS***************************************

	end	program







