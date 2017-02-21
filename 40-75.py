import matplotlib.pyplot as plt
from numpy import *
from scipy import constants as const
import csv

psiaccuracy=1/3600/180*const.pi
theta_pstep=1/60/180*const.pi

P=0.00305
Omega=2*const.pi/P
c=const.c*100
R_L=c/Omega
R_p=1e6
alpha=40/180*const.pi
zeta_view=75/180*const.pi
ratio_core=0.75
ratio_annule=0.5

phase_gamma=list()
intensity_gamma=list()
reader=csv.reader(open('gamma.txt','r'))
for line in reader:
	c1,c2=line
	c1=float(c1)
	c2=float(c2)
	phase_gamma.append(c1*360)
	intensity_gamma.append(c2)

phase_x=list()
intensity_x=list()
reader=csv.reader(open('x.txt','r'))
for line in reader:
	c1,c2=line
	c1=float(c1)
	c2=float(c2)
	phase_x.append(c1*360)
	intensity_x.append(c2)

phase_radio=list()
intensity_radio=list()
reader=csv.reader(open('radio.txt','r'))
for line in reader:
	c1,c2=line
	c1=float(c1)
	c2=float(c2)
	phase_radio.append(c1*360)
	intensity_radio.append(c2)


def corotationToMagnetic(zeta,Phi):
	costheta=cos(zeta)*cos(alpha)+sin(zeta)*sin(alpha)*cos(Phi)
	theta=arccos(costheta)
	sinpsi=sin(Phi)/sin(theta)*sin(zeta)
	cospsi=-(cos(zeta)-cos(alpha)*costheta)/(sin(alpha)*sin(theta))
	if sinpsi>=0:
		psi=arccos(cospsi)
	else:
		psi=-arccos(cospsi)
	return(theta,psi)


def magneticToCorotation(theta,psi):
	coszeta=cos(alpha)*cos(theta)-sin(alpha)*sin(theta)*cos(psi)
	zeta=arccos(coszeta)
	cosPhi=(cos(theta)-coszeta*cos(alpha))/(sin(alpha)*sin(zeta))
	sinPhi=sin(psi)/sin(zeta)*sin(theta)
	if sinPhi>=0:
		Phi=arccos(cosPhi)
	else:
		Phi=-arccos(cosPhi)
	return(zeta,Phi)


def r(R_e,theta):
	r=R_e*sin(theta)**2
	return(r)


def nullChargeSurface(psi):
	cotalpha=1/tan(alpha)
	secpsi=1/cos(psi)
	if cos(psi)>=0:
		theta_N=1/2*arccos((sqrt(8*cotalpha**2*secpsi**2+9)-cotalpha**2*secpsi**2)/(3*cotalpha**2*secpsi**2+3))
	else:
		theta_N=1/2*arccos((sqrt(8*cotalpha**2*secpsi**2+9)-cotalpha**2*secpsi**2)/(3*cotalpha**2*secpsi**2+3))
	#choose the solution near north pole.
	return(theta_N)


def coreBoundray(psi):
	theta_mu=arctan(1/(tan(alpha)*cos(psi)))
	sintheta=sin(theta_mu)
	if cos(psi)>=0:
		cosine=(-sintheta**2+sqrt(sintheta**4-10*sintheta**2+9))/3
	else:
		cosine=(-sintheta**2-sqrt(sintheta**4-10*sintheta**2+9))/3
	theta_N=1/2*arccos(cosine)
	Theta_N=arccos(cos(alpha)*cos(theta_N)-sin(alpha)*sin(theta_N)*cos(psi))
	#print(theta_N/const.pi*180)
	R_e=R_L/(sin(theta_N)**2*sin(Theta_N))
	theta_C=arcsin(sqrt(R_p/R_e))
	return(theta_C)


def coreBoundray2(psi):
	theta_mu=const.pi+arctan(1/(tan(alpha)*cos(psi)))
	sintheta=sin(theta_mu)
	if cos(psi)>=0:
		cosine=(-sintheta**2-sqrt(sintheta**4-10*sintheta**2+9))/3
	else:
		cosine=(-sintheta**2+sqrt(sintheta**4-10*sintheta**2+9))/3
	theta_N=const.pi-1/2*arccos(cosine)
	Theta_N=arccos(cos(alpha)*cos(theta_N)-sin(alpha)*sin(theta_N)*cos(psi))
	#print(theta_N/const.pi*180)
	R_e=R_L/(sin(theta_N)**2*sin(Theta_N))
	theta_C=const.pi-arcsin(sqrt(R_p/R_e))
	return(theta_C)


def annuleBoundray(psi):
	a=4*sin(alpha)**2
	b=5*sin(2*alpha)*cos(psi)
	c=4*sin(alpha)**2-6*(sin(alpha)**2*cos(psi)**2-cos(alpha)**2)
	d=-sin(2*alpha)*cos(psi)
	root=roots([a,b,c,d])
	if cos(psi)>=0:
		cotan=max(root)
	else:
		cotan=min(root)
	theta_M=arctan(1/cotan)
	R_e=R_L*(1+cotan**2)**(3/2)/sqrt(1+cotan**2-(cos(alpha)*cotan-sin(alpha)*cos(psi))**2)
	theta_A=arcsin(sqrt(R_p/R_e))
	return(theta_A)


def annuleBoundray2(psi):
	a=4*sin(alpha)**2
	b=5*sin(2*alpha)*cos(psi)
	c=4*sin(alpha)**2-6*(sin(alpha)**2*cos(psi)**2-cos(alpha)**2)
	d=-sin(2*alpha)*cos(psi)
	root=roots([a,b,c,d])
	'''if cos(psi)>=0:
		cotan=max(root)
	else:
		cotan=min(root)'''
	cotan=root[2]
	theta_M=arctan(1/cotan)
	R_e=R_L*(1+cotan**2)**(3/2)/sqrt(1+cotan**2-(cos(alpha)*cotan-sin(alpha)*cos(psi))**2)
	theta_A=const.pi-arcsin(sqrt(R_p/R_e))
	return(theta_A)


def height(zeta,Phi):
	theta_mu,psi=corotationToMagnetic(zeta,Phi)
	sintheta=sin(theta_mu)
	if cos(theta_mu)>=0:
		cosine=(-sintheta**2+sqrt(sintheta**4-10*sintheta**2+9))/3
	else:
		cosine=(-sintheta**2-sqrt(sintheta**4-10*sintheta**2+9))/3
	theta=1/2*arccos(cosine)
	#print(Phi/const.pi*180, psi/const.pi*180, theta/const.pi*180, sqrt(sintheta**4-10*sintheta**2+9))
	theta_C=coreBoundray(psi)
	theta_A=annuleBoundray(psi)
	h_core=R_p/sin(ratio_core*theta_C)**2*sin(theta)**2
	h_annule=R_p/sin(ratio_annule*theta_A+(1-ratio_annule)*theta_C)**2*sin(theta)**2
	return(h_core,h_annule)


def height2(zeta,Phi):
	theta_mu,psi=corotationToMagnetic(zeta,Phi)
	sintheta=sin(theta_mu)
	if cos(theta_mu)>=0:
		cosine=(-sintheta**2-sqrt(sintheta**4-10*sintheta**2+9))/3
	else:
		cosine=(-sintheta**2+sqrt(sintheta**4-10*sintheta**2+9))/3
	theta=const.pi-1/2*arccos(cosine)
	#print(Phi/const.pi*180, psi/const.pi*180, theta/const.pi*180, sqrt(sintheta**4-10*sintheta**2+9))
	theta_C=coreBoundray2(psi)
	theta_A=annuleBoundray2(psi)
	h_core=R_p/sin(const.pi-ratio_core*(const.pi-theta_C))**2*sin(theta)**2
	h_annule=R_p/sin(ratio_annule*theta_A+(1-ratio_annule)*theta_C)**2*sin(theta)**2
	return(h_core,h_annule)


def correction(theta_p, Phi0):
	theta_mu,psi2=corotationToMagnetic(zeta_view,Phi0)
	sintheta_mu=sin(theta_mu)
	if cos(theta_mu)>=0:
		cosine=(-sintheta_mu**2+sqrt(sintheta_mu**4-10*sintheta_mu**2+9))/3
	else:
		cosine=(-sintheta_mu**2-sqrt(sintheta_mu**4-10*sintheta_mu**2+9))/3
	theta=1/2*arccos(cosine)
	r=R_p*(sin(theta)/sin(theta_p))**2
	zeta,none=magneticToCorotation(theta, psi2)
	while True:
		psi1=psi2
		squre=1-(r*sin(zeta)/R_L)**2
		#print('squre=',squre)
		if squre<0:
			return(psi2,0,0)
		dPhi_abe=arctan(R_L*sqrt(squre)/(r*sin(zeta)))
		tantheta_mu=3*sin(2*theta)/(1+3*cos(2*theta))
		if theta<=const.pi/2 and tantheta_mu>=0:
			thata_mu=arctan(tantheta_mu)
		elif theta_mu>const.pi/2 and tantheta_mu<0:
			theta_mu=arctan(tantheta_mu)+2*const.pi
		else:
			thata_mu=arctan(tantheta_mu)+const.pi
		sigma=theta_mu-theta
		dPhi_ret=r*cos(sigma)/R_L
		dPhi_msb=1.2*(r/R_L)**3*sin(alpha)**2
		#print('dPhi_abe=',dPhi_abe,'dPhi_ret=',dPhi_ret,'dPhi_msb=',dPhi_msb)
		Phi=Phi0-(dPhi_abe)#TODO:add dPhi_ret and dPhi_msb
		theta_mu,psi2=corotationToMagnetic(zeta_view,Phi)
		sintheta_mu=sin(theta_mu)
		if cos(theta_mu)>=0:
			cosine=(-sintheta_mu**2+sqrt(sintheta_mu**4-10*sintheta_mu**2+9))/3
		else:
			cosine=(-sintheta_mu**2-sqrt(sintheta_mu**4-10*sintheta_mu**2+9))/3
		theta=1/2*arccos(cosine)
		r=R_p*(sin(theta)/sin(theta_p))**2
		zeta,none=magneticToCorotation(theta, psi2)
		if abs(psi2-psi1)<psiaccuracy:
			break

	rho=r*sin(zeta)
	return(psi2, r, rho)


def area(Phi1,Phi2):
	bordertheta=list()
	borderpsi=list()
	for Phi in linspace(Phi1,Phi2,100):
		print(Phi)
		theta_p=0
		allowtheta=list()
		allowpsi=list()
		while True:
			theta_p+=theta_pstep
			#print(theta_p)
			psi, r, rho=correction(theta_p, Phi)
			#print(psi)
			theta_A=annuleBoundray(psi)
			if rho==0:
				continue
			elif theta_p>theta_A:
				break
			else:
				allowtheta.append(theta_p/const.pi*180)
				allowpsi.append(psi)

		if Phi==Phi1:
			bordertheta=allowtheta
			borderpsi=allowpsi
		elif Phi==Phi2:
			bordertheta+=allowtheta[::-1]
			borderpsi+=allowpsi[::-1]
		else:
			bordertheta.append(allowtheta[-1])
			borderpsi.append(allowpsi[-1])
			bordertheta.insert(0,allowtheta[0])
			borderpsi.insert(0,allowpsi[0])


		

		'''theta_A=annuleBoundray(psi)
		if theta>theta_A:
			theta2=theta_A
		else:
			theta2=theta
		Theta,none=magneticToCorotation(theta,psi)
		theta1=arcsin(sqrt(R_p/R_L*sin(Theta))*sin(theta))
		borderpsi.insert(0,psi)
		bordertheta.insert(0,theta1/const.pi*180)
		borderpsi.append(psi)
		bordertheta.append(theta2/const.pi*180)'''
	return(borderpsi,bordertheta)


def area2(Phi1,Phi2):
	bordertheta=list()
	borderpsi=list()
	for Phi in linspace(Phi1,Phi2,100):
		theta_mu,psi=corotationToMagnetic(zeta_view,Phi)
		sintheta=sin(theta_mu)
		if cos(theta_mu)>=0:
			cosine=(-sintheta**2-sqrt(sintheta**4-10*sintheta**2+9))/3
		else:
			cosine=(-sintheta**2+sqrt(sintheta**4-10*sintheta**2+9))/3
		theta=const.pi-1/2*arccos(cosine)
		theta_A=annuleBoundray2(psi)
		if theta<theta_A:
			theta2=theta_A
		else:
			theta2=theta
		
		Theta,none=magneticToCorotation(theta,psi)
		theta1=const.pi-arcsin(sqrt(R_p/R_L*sin(Theta))*sin(theta))
		if theta1>theta2:
			borderpsi.insert(0,psi)
			bordertheta.insert(0,-theta1/const.pi*180)
			borderpsi.append(psi)
			bordertheta.append(-theta2/const.pi*180)
	return(borderpsi,bordertheta)


h_annule=list()
h_core=list()
h_core2=list()
h_annule2=list()
Philist=linspace(-const.pi, const.pi, 200)
psilist=linspace(-const.pi, const.pi, 200)
for Phi in Philist:
	h1,h2=height(zeta_view,Phi)
	h_core.append(h1/R_L)
	h_annule.append(h2/R_L)
	h3,h4=height2(zeta_view,Phi)
	h_core2.append(h3/R_L)
	h_annule2.append(h4/R_L)
	none,psi=corotationToMagnetic(zeta_view,Phi)

theta_core=list()
theta_annule=list()
for psi in psilist:
	t1=coreBoundray(psi)
	theta_core.append(t1/const.pi*180)
	t2=annuleBoundray(psi)
	theta_annule.append(t2/const.pi*180)
dash_core=list()
dash_annule=list()
for i in range(0,len(theta_core)):
	dash_core.append(ratio_core*theta_core[i])
	dash_annule.append(ratio_annule*theta_annule[i]+(1-ratio_annule)*theta_core[i])

theta_core2=list()
theta_annule2=list()
for psi in psilist:
	t12=coreBoundray2(psi)
	theta_core2.append(-t12/const.pi*180)
	t22=annuleBoundray2(psi)
	theta_annule2.append(-t22/const.pi*180)
dash_core2=list()
dash_annule2=list()
for i in range(0,len(theta_core)):
	dash_core2.append(ratio_core*(theta_core2[i]+180)-180)
	dash_annule2.append(ratio_annule*theta_annule2[i]+(1-ratio_annule)*theta_core2[i])


plt.xlabel("$\Phi(^\circ)$",fontsize=15)
plt.xlim(-180,180)
plt.ylabel("$\\frac{h}{R_L}$",fontsize=15)
plt.semilogy(Philist/const.pi*180, h_core, 'r',markersize=20,linewidth=1.5,label="N core(0.75)")
plt.semilogy(Philist/const.pi*180, h_annule,'g', markersize=20,linewidth=1.5,label="N annule(0.5)")
plt.semilogy(Philist/const.pi*180, h_core2, 'r--',markersize=20,linewidth=1.5,label="S core(0.75)")
plt.semilogy(Philist/const.pi*180, h_annule2,'g--', markersize=20,linewidth=1.5,label="S annule(0.5)")
plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=1,
           ncol=2, mode="expand", borderaxespad=0.)
plt.savefig('general.pdf')
plt.close('all')


borderpsi_gamma1,bordertheta_gamma1=area(-0.3983*2*const.pi,-0.0699*2*const.pi)
borderpsi_gamma2,bordertheta_gamma2=area(0.0711*2*const.pi, 0.2598*2*const.pi)
plt.figure(figsize=(10,15))
a1 = plt.subplot(321,projection='polar')
a1.plot(psilist, theta_core,'k', markersize=20,linewidth=1.5,label="critical field line")
a1.plot(psilist, dash_core,'r', markersize=20,linewidth=1.5)
a1.plot(psilist, theta_annule,'k', markersize=20,linewidth=1.5,label="last open filed line")
a1.plot(psilist, dash_annule,'g', markersize=20,linewidth=1.5)
a1.fill(borderpsi_gamma1,bordertheta_gamma1,'b')
a1.fill(borderpsi_gamma2,bordertheta_gamma2,'b')


borderpsi2_gamma1,bordertheta2_gamma1=area2(-0.3983*2*const.pi,-0.0699*2*const.pi)
borderpsi2_gamma2,bordertheta2_gamma2=area2(0.0711*2*const.pi, 0.2598*2*const.pi)
a2 = plt.subplot(322, projection='polar')
a2.set_rlim(-180,-164)
a2.plot(psilist, theta_core2,'k', markersize=20,linewidth=1.5,label="critical field line")
a2.plot(psilist, dash_core2,'r--', markersize=20,linewidth=1.5)
a2.plot(psilist, theta_annule2,'k', markersize=20,linewidth=1.5,label="last open filed line")
a2.plot(psilist, dash_annule2,'g--', markersize=20,linewidth=1.5)
a2.fill(borderpsi2_gamma1,bordertheta2_gamma1,'b')
a2.fill(borderpsi2_gamma2,bordertheta2_gamma2,'b')

h_annule=list()
h_core=list()
h_core2=list()
h_annule2=list()
Philist=linspace(-0.3983*2*const.pi,-0.0699*2*const.pi, 100)
for Phi in Philist:
	h1,h2=height(zeta_view,Phi)
	h_core.append(h1/R_L)
	h_annule.append(h2/R_L)
	h3,h4=height2(zeta_view,Phi)
	h_core2.append(h3/R_L)
	h_annule2.append(h4/R_L)
	none,psi=corotationToMagnetic(zeta_view,Phi)
b=plt.subplot(312)
b.set_xlabel("$\Phi(^\circ)$",fontsize=15)
b.set_xlim(-180,180)
b.set_ylabel("$\\frac{h}{R_L}$",fontsize=15)
b.semilogy(Philist/const.pi*180, h_core, 'r',markersize=20,linewidth=1.5,label="N core(0.75)")
b.semilogy(Philist/const.pi*180, h_annule,'g', markersize=20,linewidth=1.5,label="N annule(0.5)")
b.semilogy(Philist/const.pi*180, h_core2, 'r--',markersize=20,linewidth=1.5,label="S core(0.75)")
b.semilogy(Philist/const.pi*180, h_annule2,'g--', markersize=20,linewidth=1.5,label="S annule(0.5)")
h_annule=list()
h_core=list()
h_core2=list()
h_annule2=list()
Philist=linspace(0.0711*2*const.pi, 0.2598*2*const.pi, 100)
for Phi in Philist:
	h1,h2=height(zeta_view,Phi)
	h_core.append(h1/R_L)
	h_annule.append(h2/R_L)
	h3,h4=height2(zeta_view,Phi)
	h_core2.append(h3/R_L)
	h_annule2.append(h4/R_L)
	none,psi=corotationToMagnetic(zeta_view,Phi)
b.semilogy(Philist/const.pi*180, h_core, 'r',markersize=20,linewidth=1.5)
b.semilogy(Philist/const.pi*180, h_annule,'g', markersize=20,linewidth=1.5)
b.semilogy(Philist/const.pi*180, h_core2, 'r--',markersize=20,linewidth=1.5)
b.semilogy(Philist/const.pi*180, h_annule2,'g--', markersize=20,linewidth=1.5)
b.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=1,
           ncol=2, mode="expand", borderaxespad=0.)
for Phi in Philist:
	h_c,h_a=height(zeta_view,Phi)
	h_c/=R_L
	h_a/=R_L
c=plt.subplot(313)
c.set_xlim(-180,180)
c.set_ylabel("$\gamma$-ray flux",fontsize=15)
c.plot(phase_gamma,intensity_gamma,'b',linewidth=1.5)
plt.savefig("B1821gamma.pdf")
plt.close('all')


borderpsi_x1,bordertheta_x1=area(-0.2941*2*const.pi,-0.2132*2*const.pi)
borderpsi_x2,bordertheta_x2=area(0.2696*2*const.pi,0.3407*2*const.pi)
plt.figure(figsize=(10,10))
a1 = plt.subplot(321,projection='polar')
a1.plot(psilist, theta_core,'k', markersize=20,linewidth=1.5,label="critical field line")
a1.plot(psilist, dash_core,'r', markersize=20,linewidth=1.5)
a1.plot(psilist, theta_annule,'k', markersize=20,linewidth=1.5,label="last open filed line")
a1.plot(psilist, dash_annule,'g', markersize=20,linewidth=1.5)
a1.fill(borderpsi_x1,bordertheta_x1,'b')
a1.fill(borderpsi_x2,bordertheta_x2,'b')

borderpsi2_x1,bordertheta2_x1=area2(-0.2941*2*const.pi,-0.2132*2*const.pi)
borderpsi2_x2,bordertheta2_x2=area2(0.2696*2*const.pi,0.3407*2*const.pi)
a2 = plt.subplot(322, projection='polar')
a2.set_rlim(-180,-164)
a2.plot(psilist, theta_core2,'k', markersize=20,linewidth=1.5,label="critical field line")
a2.plot(psilist, dash_core2,'r--', markersize=20,linewidth=1.5)
a2.plot(psilist, theta_annule2,'k', markersize=20,linewidth=1.5,label="last open filed line")
a2.plot(psilist, dash_annule2,'g--', markersize=20,linewidth=1.5)
a2.fill(borderpsi2_x1,bordertheta2_x1,'b')
a2.fill(borderpsi2_x2,bordertheta2_x2,'b')

h_annule=list()
h_core=list()
h_core2=list()
h_annule2=list()
Philist=linspace(-0.2941*2*const.pi,-0.2132*2*const.pi, 100)
for Phi in Philist:
	h1,h2=height(zeta_view,Phi)
	h_core.append(h1/R_L)
	h_annule.append(h2/R_L)
	h3,h4=height2(zeta_view,Phi)
	h_core2.append(h3/R_L)
	h_annule2.append(h4/R_L)
	none,psi=corotationToMagnetic(zeta_view,Phi)
b=plt.subplot(312)
b.set_xlabel("$\Phi(^\circ)$",fontsize=15)
b.set_xlim(-180,180)
b.set_ylabel("$\\frac{h}{R_L}$",fontsize=15)
b.semilogy(Philist/const.pi*180, h_core, 'r',markersize=20,linewidth=1.5,label="N core(0.75)")
b.semilogy(Philist/const.pi*180, h_annule,'g', markersize=20,linewidth=1.5,label="N annule(0.5)")
b.semilogy(Philist/const.pi*180, h_core2, 'r--',markersize=20,linewidth=1.5,label="S core(0.75)")
b.semilogy(Philist/const.pi*180, h_annule2,'g--', markersize=20,linewidth=1.5,label="S annule(0.5)")
h_annule=list()
h_core=list()
h_core2=list()
h_annule2=list()
Philist=linspace(0.2696*2*const.pi,0.3407*2*const.pi, 100)
for Phi in Philist:
	h1,h2=height(zeta_view,Phi)
	h_core.append(h1/R_L)
	h_annule.append(h2/R_L)
	h3,h4=height2(zeta_view,Phi)
	h_core2.append(h3/R_L)
	h_annule2.append(h4/R_L)
	none,psi=corotationToMagnetic(zeta_view,Phi)
b.semilogy(Philist/const.pi*180, h_core, 'r',markersize=20,linewidth=1.5)
b.semilogy(Philist/const.pi*180, h_annule,'g', markersize=20,linewidth=1.5)
b.semilogy(Philist/const.pi*180, h_core2, 'r--',markersize=20,linewidth=1.5)
b.semilogy(Philist/const.pi*180, h_annule2,'g--', markersize=20,linewidth=1.5)
b.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=1,
           ncol=2, mode="expand", borderaxespad=0.)
for Phi in Philist:
	h_c,h_a=height(zeta_view,Phi)
	h_c/=R_L
	h_a/=R_L
c=plt.subplot(313)
c.set_xlim(-180,180)
c.set_ylabel("X-ray flux",fontsize=15)
c.plot(phase_x,intensity_x,'b',linewidth=1.5)
plt.savefig("B1821x.pdf")
plt.close('all')


borderpsi_radio1,bordertheta_radio1=area(-0.2970*2*const.pi,-0.2170*2*const.pi)
borderpsi_radio2,bordertheta_radio2=area(0.0196*2*const.pi,0.0912*2*const.pi)
borderpsi_radio3,bordertheta_radio3=area(0.1851*2*const.pi,0.2749*2*const.pi)
plt.figure(figsize=(10,10))
a1 = plt.subplot(321,projection='polar')
a1.plot(psilist, theta_core,'k', markersize=20,linewidth=1.5,label="critical field line")
a1.plot(psilist, dash_core,'r', markersize=20,linewidth=1.5)
a1.plot(psilist, theta_annule,'k', markersize=20,linewidth=1.5,label="last open filed line")
a1.plot(psilist, dash_annule,'g', markersize=20,linewidth=1.5)
a1.fill(borderpsi_radio1,bordertheta_radio1,'b')
a1.fill(borderpsi_radio2,bordertheta_radio2,'b')
a1.fill(borderpsi_radio3,bordertheta_radio3,'b')

borderpsi2_radio1,bordertheta2_radio1=area2(-0.2970*2*const.pi,-0.2170*2*const.pi)
borderpsi2_radio2,bordertheta2_radio2=area2(0.0196*2*const.pi,0.0912*2*const.pi)
borderpsi2_radio3,bordertheta2_radio3=area2(0.1851*2*const.pi,0.2749*2*const.pi)
a2 = plt.subplot(322, projection='polar')
a2.set_rlim(-180,-164)
a2.plot(psilist, theta_core2,'k', markersize=20,linewidth=1.5,label="critical field line")
a2.plot(psilist, dash_core2,'r--', markersize=20,linewidth=1.5)
a2.plot(psilist, theta_annule2,'k', markersize=20,linewidth=1.5,label="last open filed line")
a2.plot(psilist, dash_annule2,'g--', markersize=20,linewidth=1.5)
a2.fill(borderpsi2_radio1,bordertheta2_radio1,'b')
a2.fill(borderpsi2_radio2,bordertheta2_radio2,'b')
a2.fill(borderpsi2_radio3,bordertheta2_radio3,'b')

h_annule=list()
h_core=list()
h_core2=list()
h_annule2=list()
Philist=linspace(-0.2970*2*const.pi,-0.2170*2*const.pi, 100)
for Phi in Philist:
	h1,h2=height(zeta_view,Phi)
	h_core.append(h1/R_L)
	h_annule.append(h2/R_L)
	h3,h4=height2(zeta_view,Phi)
	h_core2.append(h3/R_L)
	h_annule2.append(h4/R_L)
	none,psi=corotationToMagnetic(zeta_view,Phi)
b=plt.subplot(312)
b.set_xlabel("$\Phi(^\circ)$",fontsize=15)
b.set_xlim(-180,180)
b.set_ylabel("$\\frac{h}{R_L}$",fontsize=15)
b.semilogy(Philist/const.pi*180, h_core, 'r',markersize=20,linewidth=1.5,label="N core(0.75)")
b.semilogy(Philist/const.pi*180, h_annule,'g', markersize=20,linewidth=1.5,label="N annule(0.5)")
b.semilogy(Philist/const.pi*180, h_core2, 'r--',markersize=20,linewidth=1.5,label="S core(0.75)")
b.semilogy(Philist/const.pi*180, h_annule2,'g--', markersize=20,linewidth=1.5,label="S annule(0.5)")
h_annule=list()
h_core=list()
h_core2=list()
h_annule2=list()
Philist=linspace(0.0196*2*const.pi,0.0912*2*const.pi, 100)
for Phi in Philist:
	h1,h2=height(zeta_view,Phi)
	h_core.append(h1/R_L)
	h_annule.append(h2/R_L)
	h3,h4=height2(zeta_view,Phi)
	h_core2.append(h3/R_L)
	h_annule2.append(h4/R_L)
	none,psi=corotationToMagnetic(zeta_view,Phi)
b.semilogy(Philist/const.pi*180, h_core, 'r',markersize=20,linewidth=1.5)
b.semilogy(Philist/const.pi*180, h_annule,'g', markersize=20,linewidth=1.5)
b.semilogy(Philist/const.pi*180, h_core2, 'r--',markersize=20,linewidth=1.5)
b.semilogy(Philist/const.pi*180, h_annule2,'g--', markersize=20,linewidth=1.5)
h_annule=list()
h_core=list()
h_core2=list()
h_annule2=list()
Philist=linspace(0.1851*2*const.pi,0.2749*2*const.pi, 100)
for Phi in Philist:
	h1,h2=height(zeta_view,Phi)
	h_core.append(h1/R_L)
	h_annule.append(h2/R_L)
	h3,h4=height2(zeta_view,Phi)
	h_core2.append(h3/R_L)
	h_annule2.append(h4/R_L)
	none,psi=corotationToMagnetic(zeta_view,Phi)
b.semilogy(Philist/const.pi*180, h_core, 'r',markersize=20,linewidth=1.5)
b.semilogy(Philist/const.pi*180, h_annule,'g', markersize=20,linewidth=1.5)
b.semilogy(Philist/const.pi*180, h_core2, 'r--',markersize=20,linewidth=1.5)
b.semilogy(Philist/const.pi*180, h_annule2,'g--', markersize=20,linewidth=1.5)
b.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=1,
           ncol=2, mode="expand", borderaxespad=0.)
for Phi in Philist:
	h_c,h_a=height(zeta_view,Phi)
	h_c/=R_L
	h_a/=R_L
c=plt.subplot(313)
c.set_xlim(-180,180)
c.set_ylabel("radio flux",fontsize=15)
c.plot(phase_radio,intensity_radio,'b',linewidth=1.5)
plt.savefig("B1821radio.pdf")
plt.close('all')