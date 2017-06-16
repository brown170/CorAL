// <<BEGIN-copyright>>
// 
//                 The GNU General Public License (GPL) Version 2, June 1991
// 
// Copyright (c) 2013, Lawrence Livermore National Security, LLC. Produced at the Lawrence 
// Livermore National Laboratory. Written by Ron Soltz (soltz1@llnl.gov), David A. Brown 
// (dbrown@bnl.gov) and Scott Pratt (pratts@pa.msu.edu).
// 
// CODE-CODE-643336 All rights reserved. 
// 
// This file is part of CorAL, Version: 1.17.
// 
// Please see the file LICENSE.TXT in the main directory of this source code distribution.
// 
// This program is free software; you can redistribute it and/or modify it under the terms of 
// the GNU General Public License (as published by the Free Software Foundation) version 2, 
// dated June 1991.
// 
// This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
// without even the IMPLIED WARRANTY OF MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
// See the terms and conditions of the GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License along with this program; 
// if not, write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, 
// MA 02111-1307 USA
// 
// <<END-copyright>>
#ifndef __CWAVEFUNCTION_WF_Xipi_phaseshift_CC
#define __CWAVEFUNCTION_WF_Xipi_phaseshift_CC
#include "wavefunction.h"

// See wavefunction.cc for parent class

CWaveFunction_Xipi_phaseshift::CWaveFunction_Xipi_phaseshift(string  parsfilename){
	int iq,ichannel;
	double q;

	ParsInit(parsfilename);

	m1=1321.3;
	m2=139.58;
	IDENTICAL=0;
	q1q2=-1;
	if(COULOMB==0) q1q2=0;
	nchannels=1;

	ellmax=1;
	InitArrays();

	ell[0]=1;

	InitWaves();

	channelweight[0]=(4.0/3.0);

	ichannel=0;
	for(iq=0;iq<nqmax;iq++){
		q=qarray[iq];
		get_phaseshifts_Xistar(q,delta[0][iq],ddeltadq[0][iq]);
		if (q1q2!=0)
			CoulWave::phaseshift_CoulombCorrect(ell[ichannel],q,eta[iq],
			delta[ichannel][iq],ddeltadq[ichannel][iq]);
		printf("ichannel=%d, q=%g, E=%g, delta=%g, ddelta/dq=%g\n",
			ichannel,q,sqrt(q*q+m1*m1)+sqrt(q*q+m2*m2),
			delta[ichannel][iq]*180.0/PI,ddeltadq[ichannel][iq]*180.0/PI);
		Wepsilon[ichannel][iq]=ddeltadq[ichannel][iq]
			-GetIW(ell[ichannel],epsilon,q,q1q2,eta[iq],delta[ichannel][iq])
			+GetIW(ell[ichannel],epsilon,q,q1q2,eta[iq],0.0);
		Wepsilon[ichannel][iq]=3.0*Wepsilon[ichannel][iq]
			/(4.0*PI*pow(epsilon,3));
	}
}

double CWaveFunction_Xipi_phaseshift::CalcPsiSquared(int iq,double r,double ctheta){
	double psisquared,x,dpsi2,q;
	complex<double> psi,hstar1,psi0,delpsi;
	int ichannel;

	if(iq>=nqmax){
		psisquared=1.0;
	}
	else{
		q=qarray[iq];
		x=q*r/HBARC;
		psi0=planewave[iq]->planewave(r,ctheta);

		if(r<epsilon){
			psisquared=real(psi0*conj(psi0));
			for(ichannel=0;ichannel<nchannels;ichannel++){
				dpsi2=channelweight[ichannel]*2.0*PI*Wepsilon[ichannel][iq]
					*pow(HBARC,3)/(q*q);
				psisquared+=dpsi2;
			}
		}
		else{
			hstar1=partwave[1][iq]->GetPhiIncoming(r)/x;
			delpsi=0.5*hstar1*(Misc::ceiphi(-2.0*delta[0][iq])-1.0)
				*ci*(3.0)*ctheta;
			// m_s equals 1/2, ell=1
			psi=psi0+(4.0/9.0)*delpsi;
			psisquared=real(psi*conj(psi));

			//This accounts for the other three channels, where the incoming
			// wave was either in the Xi0-pi0 channel or if m_ell=1, m_s=-1/2.
			//First, the channel where m_ell=0,m_s=1/2 but the Xi0-pi0
			psisquared+=(8.0/81.0)*real(delpsi*conj(delpsi));
			//Now, the other two channels where m_ell=1, m_s=-1/2
			delpsi=0.5*hstar1*(Misc::ceiphi(-2.0*delta[0][iq])-1.0)
				*ci*(3.0)*sqrt((1.0-ctheta*ctheta)/2.0);
			psisquared+=(12.0/81.0)*real(delpsi*conj(delpsi));

		}
	}

	psisquared*=RelativisticCorrection(r,iq);
	return psisquared;

}

void CWaveFunction_Xipi_phaseshift::get_phaseshifts_Xistar(double q,double &delta,double &ddeltadq){
	double tandelta,Epi,EXi,Lambda;
	const double MXi=1321.3,Lambda0=9.1,MXistar=1530.0;
	const double q0=sqrt(Misc::triangle(MXistar,MXi,MPI));
	EXi=sqrt(MXi*MXi+q*q);
	Epi=sqrt(MPI*MPI+q*q);
	Lambda=Lambda0*pow(q/q0,3);
	tandelta=0.5*Lambda/(MXistar-EXi-Epi);
	delta=atan(tandelta);
	if(delta<0.0) delta+=PI;
	ddeltadq=pow(cos(delta),2)*tandelta*((3.0/q)+(q/(MXistar-EXi-Epi))
		*((1.0/EXi)+(1.0/Epi)));
}

#endif
