/*
 SuperCollider real time audio synthesis system
 Copyright (c) 2002 James McCartney. All rights reserved.
 http://www.audiosynth.com
 
 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 2 of the License, or
 (at your option) any later version.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301  USA
 */

// NLUGens by yota morimoto (http://yota.tehis.net/)

#include "SC_PlugIn.h"

#define LATTICE 10
#define GENEBIT 16

static InterfaceTable *ft;

struct Logist0 : public Unit {
	double x;
	float counter;
};
struct Logist1 : public Unit {
	double x, xm, frac;
	float counter;
};
struct Logist3 : public Unit {
	double x, xm1, xm2, xm3, frac, c0, c1, c2, c3;
	float counter;
};

struct Nagumo : public Unit {
	double u, v;
};

struct FIS : public Unit {
};

struct CML0 : public Unit {
	double x[LATTICE];
	float counter;
};
struct CML1 : public Unit {
	double xm, frac;
	double x[LATTICE];
	float counter;
};
struct CML3 : public Unit {
	double xm1, xm2, xm3, frac, c0, c1, c2, c3;
	double x[LATTICE];
	float counter;
};

struct GCM0 : public Unit {
	double x[LATTICE];
	float counter;
};
struct GCM1 : public Unit {
	double xm, frac;
	double x[LATTICE];
	float counter;
};
struct GCM3 : public Unit {
	double xm1, xm2, xm3, frac, c0, c1, c2, c3;
	double x[LATTICE];
	float counter;
};

struct HCM0 : public Unit {
	unsigned short x[GENEBIT];
	float counter;
};
struct HCM1 : public Unit {
	double xm, frac;
	unsigned short x[GENEBIT];
	float counter;
};
struct HCM3 : public Unit {
	double xm1, xm2, xm3, frac, c0, c1, c2, c3;
	unsigned short x[GENEBIT];
	float counter;
};

struct TLogist : public Logist0 {
	double trig;
};

extern "C" {
	void Logist0_next(Logist0 *unit, int inNumSamples);
	void Logist0_Ctor(Logist0 *unit);
	void Logist1_next(Logist1 *unit, int inNumSamples);
	void Logist1_Ctor(Logist1 *unit);
	void Logist3_next(Logist3 *unit, int inNumSamples);
	void Logist3_Ctor(Logist3 *unit);
	void Nagumo_next(Nagumo *unit, int inNumSamples);
	void Nagumo_Ctor(Nagumo *unit);
	void FIS_next(FIS *unit, int inNumSamples);
	void FIS_Ctor(FIS *unit);
	void CML0_next(CML0 *unit, int inNumSamples);
	void CML0_Ctor(CML0 *unit);
	void CML1_next(CML1 *unit, int inNumSamples);
	void CML1_Ctor(CML1 *unit);
	void CML3_next(CML3 *unit, int inNumSamples);
	void CML3_Ctor(CML3 *unit);
	void GCM0_next(GCM0 *unit, int inNumSamples);
	void GCM0_Ctor(GCM0 *unit);
	void GCM1_next(GCM1 *unit, int inNumSamples);
	void GCM1_Ctor(GCM1 *unit);
	void GCM3_next(GCM3 *unit, int inNumSamples);
	void GCM3_Ctor(GCM3 *unit);
	void HCM0_next(HCM0 *unit, int inNumSamples);
	void HCM0_Ctor(HCM0 *unit);
	void HCM1_next(HCM1 *unit, int inNumSamples);
	void HCM1_Ctor(HCM1 *unit);
	void HCM3_next(HCM3 *unit, int inNumSamples);
	void HCM3_Ctor(HCM3 *unit);
	void TLogist_next(TLogist *unit, int inNumSamples);
	void TLogist_Ctor(TLogist *unit);
}

// calc 3rd order interpolation coefs from four points
static inline void ipol3Coef(double xm3, double xm2, double xm1, double x, double &c0, double &c1, double &c2, double &c3)
{
	c0 = xm2;
	c1 = 0.5l * (xm1 - xm3);
	c2 = xm3 - (2.5l * xm2) + xm1 + xm1 - 0.5l * x;
	c3 = 0.5l * (x - xm3) + 1.5l * (xm2 - xm1);
}
// do 3rd order interpolation using coefs
static inline double ipol3(double frac, double c0, double c1, double c2, double c3)
{
	return ((c3 * frac + c2) * frac + c1) * frac + c0;
}

inline double logist(double r, double x);
inline double logist(double r, double x)
{
	return 1.l - r * x * x;
}

void Logist0_next(Logist0 *unit, int inNumSamples)
{
	float *out = ZOUT(0);
	float freq = ZIN0(0);
	double r = ZIN0(1);
	
	double x = unit->x;
	float counter = unit->counter;
	
	float spc;
	if(freq < SAMPLERATE)
		spc = SAMPLERATE / sc_max(freq, 0.001f);
	else spc = 1.f;
	
	LOOP(inNumSamples,
		 if(counter >= spc){
			 counter -= spc;
			 x = logist(r, x);
		 }
		 counter++;		
		 ZXP(out) = x;
		 )
	unit->x = x;
	unit->counter = counter;
}

void Logist0_Ctor(Logist0 *unit)
{
	SETCALC(Logist0_next);
	unit->x = IN0(2);
	unit->counter = 0.f;
	Logist0_next(unit, 1);
}

void Logist1_next(Logist1 *unit, int inNumSamples)
{
	float *out = ZOUT(0);
	float freq = ZIN0(0);
	double r = ZIN0(1);

	double x = unit->x;
	double xm = unit->xm;
	double frac = unit->frac;
	float counter = unit->counter;

	float spc;
	float slope;
	if(freq < SAMPLERATE){
		spc = SAMPLERATE / sc_max(freq, 0.001f);
		slope = 1.f / spc;
	} 
	else spc = slope = 1.f;

	double dx = x - xm;
	
	LOOP(inNumSamples,
		if(counter >= spc){
			counter -= spc;
			frac = 0.l;
			xm = x;
			x = logist(r, x);
			dx = x - xm;
		}
		counter++;		
		ZXP(out) = xm + dx * frac;
		frac += slope;
	)
	unit->x = x;
	unit->xm = xm;
	unit->frac = frac;
	unit->counter = counter;
}

void Logist1_Ctor(Logist1 *unit)
{
	SETCALC(Logist1_next);
	unit->x = unit->xm = IN0(2);
	unit->counter = 0.f;
	unit->frac = 0.l;
	Logist1_next(unit, 1);
}

void Logist3_next(Logist3 *unit, int inNumSamples)
{
	float *out = ZOUT(0);
	float freq = ZIN0(0);
	double r = ZIN0(1);

	double x = unit->x;
	double xm1 = unit->xm1;
	double xm2 = unit->xm2;
	double xm3 = unit->xm3;
	double c0 = unit->c0;
	double c1 = unit->c1;
	double c2 = unit->c2;
	double c3 = unit->c3;
		
	float counter = unit->counter;
	double frac = unit->frac;

	float spc;
	double slope;
	if(freq < SAMPLERATE){
		spc = SAMPLERATE / sc_max(freq, 0.001f);
		slope = 1.l / spc;
	} 
	else spc = slope = 1.l;

	LOOP(inNumSamples,
		if(counter >= spc){
			counter -= spc;
			frac = 0.l;
			xm3 = xm2;
			xm2 = xm1;
			xm1 = x;
			x = logist(r, x);
			ipol3Coef(xm3, xm2, xm1, x, c0, c1, c2, c3);
		}
		counter++;		
		ZXP(out) = ipol3(frac, c0, c1, c2, c3);
		frac += slope;
	)
	unit->x = x;
	unit->counter = counter;
	unit->frac = frac;
	unit->xm1 = xm1;
	unit->xm2 = xm2;
	unit->xm3 = xm3;
	unit->c0 = c0;
	unit->c1 = c1;
	unit->c2 = c2;
	unit->c3 = c3;
}

void Logist3_Ctor(Logist3 *unit)
{
	SETCALC(Logist3_next);
	unit->x = unit->xm1 = unit->xm2 = unit->xm3 = IN0(2);
	unit->c0 = unit->c1 = unit->c2 = unit->c3 = 0.l;
	unit->counter = 0.f;
	unit->frac = 0.l;
	Logist3_next(unit, 1);
}

void Nagumo_next(Nagumo *unit, int inNumSamples)
{
	float *out = ZOUT(0);

	double uh = ZIN0(0);
	double vh = ZIN0(1);
	float *pulse = ZIN(2);
		
	double u = unit->u;
	double v = unit->v;
	
	LOOP(inNumSamples,
		float zPulse = ZXP(pulse);
		u += uh * (10.l * (- v + u - 0.3333333l*u*u*u + zPulse));
		v += vh * (u - 0.8l * v + 0.7l);
		ZXP(out) = u * 0.3;
	)
	unit->u = u;
	unit->v = v;
}

void Nagumo_Ctor(Nagumo *unit)
{
	SETCALC(Nagumo_next);
	unit->u = 0.1l;
	unit->v = 0.l;
	Nagumo_next(unit, 1);
}

void FIS_next(FIS *unit, int inNumSamples)
{
	float *out = ZOUT(0);
	float *r = ZIN(0);
	float *x = ZIN(1);
	int n = ZIN0(2);

	LOOP(inNumSamples,
		double zx = ZXP(x);
		double zr = ZXP(r);
		for(int i=0; i<n; i++)
			zx = sin(zr * zx);
		ZXP(out) = zx;
	)	
}

void FIS_Ctor(FIS *unit)
{
	SETCALC(FIS_next);
	FIS_next(unit, 1);
}

void CML0_next(CML0 *unit, int inNumSamples)
{
	float *out = ZOUT(0);
	float freq = ZIN0(0);
	double r = ZIN0(1);
	double g = ZIN0(2);
	double x[LATTICE];

	for (int i=0; i<LATTICE; i++) x[i] = unit->x[i];	
	float counter = unit->counter;

	float spc;
	float slope;
	if(freq < SAMPLERATE){
		spc = SAMPLERATE / sc_max(freq, 0.001f);
		slope = 1.f / spc;
	} 
	else spc = slope = 1.f;

	LOOP(inNumSamples,
		if(counter >= spc){
			counter -= spc;
			for (int i=1; i<LATTICE-1; i++)
				// x[i] = x[i] rather than new = old
				x[i] = (1.l - g) * logist(r, x[i]) + 0.5 * g * (logist(r, x[i+1]) + logist(r, x[i-1]));//no wrapping
		}
		counter++;
		ZXP(out) = x[5];
	)
	unit->counter = counter;
	for (int i=0; i<LATTICE; i++) unit->x[i] = x[i];
}

void CML0_Ctor(CML0 *unit)
{
	SETCALC(CML0_next);
	for (int i=0; i<LATTICE; i++) unit->x[i] = IN0(3);
	unit->counter = 0.f;
	CML0_next(unit, 1);
}

void CML1_next(CML1 *unit, int inNumSamples)
{
	float *out = ZOUT(0);
	float freq = ZIN0(0);
	double r = ZIN0(1);
	double g = ZIN0(2);
	double x[LATTICE];

	for (int i=0; i<LATTICE; i++) x[i] = unit->x[i];
	double xm = unit->xm;
	double frac = unit->frac;
	float counter = unit->counter;
		
	float spc;
	float slope;
	if(freq < SAMPLERATE){
		spc = SAMPLERATE / sc_max(freq, 0.001f);
		slope = 1.f / spc;
	} 
	else spc = slope = 1.f;
	
	double dx = x[5] - xm;
	
	LOOP(inNumSamples,
		 if(counter >= spc){
			 counter -= spc;
			 frac = 0.l;
			 xm = x[5];
			 for (int i=1; i<LATTICE-1; i++)
				 x[i] = (1.l - g) * logist(r, x[i]) + 0.5 * g * (logist(r, x[i+1]) + logist(r, x[i-1]));
			 dx = x[5] - xm;
		 }
		 counter++;		
		 ZXP(out) = xm + dx * frac;
		 frac += slope;
	)
	for (int i=0; i<LATTICE; i++) unit->x[i] = x[i];
	unit->xm = xm;
	unit->frac = frac;
	unit->counter = counter;
}

void CML1_Ctor(CML1 *unit)
{
	SETCALC(CML1_next);
	for (int i=0; i<LATTICE; i++) unit->x[i] = IN0(3);
	unit->counter = 0.f;
	unit->frac = 0.l;
	CML1_next(unit, 1);
}

void CML3_next(CML3 *unit, int inNumSamples)
{
	float *out = ZOUT(0);
	float freq = ZIN0(0);
	double r = ZIN0(1);
	double g = ZIN0(2);
	double x[LATTICE];

	for (int i=0; i<LATTICE; i++) x[i] = unit->x[i];	
	double xm1 = unit->xm1;
	double xm2 = unit->xm2;
	double xm3 = unit->xm3;
	double c0 = unit->c0;
	double c1 = unit->c1;
	double c2 = unit->c2;
	double c3 = unit->c3;
	double frac = unit->frac;		
	float counter = unit->counter;

	float spc;
	double slope;
	if(freq < SAMPLERATE){
		spc = SAMPLERATE / sc_max(freq, 0.001f);
		slope = 1.l / spc;
	} 
	else spc = slope = 1.l;

	LOOP(inNumSamples,
		if(counter >= spc){
			counter -= spc;
			frac = 0.l;
			xm3 = xm2;
			xm2 = xm1;
			xm1 = x[5];
			for (int i=1; i<LATTICE-1; i++)
				x[i] = (1.l - g) * logist(r, x[i]) + 0.5 * g * (logist(r, x[i+1]) + logist(r, x[i-1]));
			ipol3Coef(xm3, xm2, xm1, x[5], c0, c1, c2, c3);
		}
		counter++;
		ZXP(out) = ipol3(frac, c0, c1, c2, c3);
		frac += slope;
	)
	unit->counter = counter;
	unit->frac = frac;
	unit->xm1 = xm1;
	unit->xm2 = xm2;
	unit->xm3 = xm3;
	unit->c0 = c0;
	unit->c1 = c1;
	unit->c2 = c2;
	unit->c3 = c3;
	for (int i=0; i<LATTICE; i++) unit->x[i] = x[i];
}

void CML3_Ctor(CML3 *unit)
{
	SETCALC(CML3_next);
	for (int i=0; i<LATTICE; i++) unit->x[i] = IN0(3);
	unit->xm1 = unit->xm2 = unit->xm3 = unit->x[0];
	unit->c0 = unit->c1 = unit->c2 = unit->c3 = 0.l;
	unit->counter = 0.f;
	unit->frac = 0.l;
	CML3_next(unit, 1);
}

void GCM0_next(GCM0 *unit, int inNumSamples)
{
	float *out = ZOUT(0);
	float freq = ZIN0(0);
	double r = ZIN0(1);
	double g = ZIN0(2);
	double x[LATTICE];
	double reciprocal = 1.l / LATTICE;

	for (int i=0; i<LATTICE; i++) x[i] = unit->x[i];
	float counter = unit->counter;

	float spc;
	if(freq < SAMPLERATE)
		spc = SAMPLERATE / sc_max(freq, 0.001f);
	else spc = 1.f;
	
	double sum = 0;
	for (int i=0; i<LATTICE; i++) sum += logist(r, x[i]);//in theory should be in LOOP
	
	LOOP(inNumSamples,
		if(counter >= spc){
			counter -= spc;
			for (int i=0; i<LATTICE; i++)
				// x[i] = x[i] rather than new = old
				x[i] = (1.l - g) * logist(r, x[i]) + g * reciprocal * sum;
		}
		counter++;
		ZXP(out) = x[5];
	)
	unit->counter = counter;
	for (int i=0; i<LATTICE; i++) unit->x[i] = x[i];
}

void GCM0_Ctor(GCM0 *unit)
{
	SETCALC(GCM0_next);
	for (int i=0; i<LATTICE; i++) unit->x[i] = IN0(3);
	unit->counter = 0.f;
	GCM0_next(unit, 1);
}

void GCM1_next(GCM1 *unit, int inNumSamples)
{
	float *out = ZOUT(0);
	float freq = ZIN0(0);
	double r = ZIN0(1);
	double g = ZIN0(2);
	double x[LATTICE];
	double reciprocal = 1.l / LATTICE;
	
	for (int i=0; i<LATTICE; i++) x[i] = unit->x[i];
	double xm = unit->xm;
	double frac = unit->frac;
	float counter = unit->counter;
	
	float spc;
	float slope;
	if(freq < SAMPLERATE){
		spc = SAMPLERATE / sc_max(freq, 0.001f);
		slope = 1.f / spc;
	} 
	else spc = slope = 1.f;
	
	double dx = x[5] - xm;
	
	double sum = 0;
	for(int i=0; i<LATTICE; i++) sum += logist(r, x[i]);
	
	LOOP(inNumSamples,
		 if(counter >= spc){
			 counter -= spc;
			 frac = 0.l;
			 xm = x[5];
			 for (int i=0; i<LATTICE; i++) x[i] = (1.l - g) * logist(r, x[i]) + g * reciprocal * sum;
			 dx = x[5] - xm;
		 }
		 counter++;
		 ZXP(out) = xm + dx * frac;
		 frac += slope;
	)
	for(int i=0; i<LATTICE; i++) unit->x[i] = x[i];
	unit->xm = xm;
	unit->frac = frac;
	unit->counter = counter;
}

void GCM1_Ctor(GCM1 *unit)
{
	SETCALC(GCM1_next);
	for(int i=0; i<LATTICE; i++) unit->x[i] = IN0(3);
	unit->counter = 0.f;
	unit->frac = 0.l;
	GCM1_next(unit, 1);
}

void GCM3_next(GCM3 *unit, int inNumSamples)
{
	float *out = ZOUT(0);
	float freq = ZIN0(0);
	double r = ZIN0(1);
	double g = ZIN0(2);
	double x[LATTICE];
	double reciprocal = 1.l / LATTICE;
	
	for (int i=0; i<LATTICE; i++) x[i] = unit->x[i];	
	double xm1 = unit->xm1;
	double xm2 = unit->xm2;
	double xm3 = unit->xm3;
	double c0 = unit->c0;
	double c1 = unit->c1;
	double c2 = unit->c2;
	double c3 = unit->c3;
	double frac = unit->frac;		
	float counter = unit->counter;
	
	float spc;
	double slope;
	if(freq < SAMPLERATE){
		spc = SAMPLERATE / sc_max(freq, 0.001f);
		slope = 1.l / spc;
	} 
	else spc = slope = 1.l;

	double sum = 0;
	for (int i=0; i<LATTICE; i++) sum += logist(r, x[i]);
	
	LOOP(inNumSamples,
		 if(counter >= spc){
			 counter -= spc;
			 frac = 0.l;
			 xm3 = xm2;
			 xm2 = xm1;
			 xm1 = x[5];
			 for (int i=0; i<LATTICE; i++) x[i] = (1.l - g) * logist(r, x[i]) + g * reciprocal * sum;
 			 ipol3Coef(xm3, xm2, xm1, x[5], c0, c1, c2, c3);
		 }
		 counter++;
		 ZXP(out) = ipol3(frac, c0, c1, c2, c3);
		 frac += slope;
	)
	unit->counter = counter;
	unit->frac = frac;
	unit->xm1 = xm1;
	unit->xm2 = xm2;
	unit->xm3 = xm3;
	unit->c0 = c0;
	unit->c1 = c1;
	unit->c2 = c2;
	unit->c3 = c3;
	for (int i=0; i<LATTICE; i++) unit->x[i] = x[i];
}

void GCM3_Ctor(GCM3 *unit)
{
	SETCALC(GCM3_next);
	for (int i=0; i<LATTICE; i++) unit->x[i] = IN0(3);
	unit->xm1 = unit->xm2 = unit->xm3 = unit->x[0];
	unit->c0 = unit->c1 = unit->c2 = unit->c3 = 0.l;
	unit->counter = 0.f;
	unit->frac = 0.l;
	GCM3_next(unit, 1);
}

inline unsigned flip(unsigned x, unsigned bit);
inline unsigned flip(unsigned x, unsigned bit)
{
	return x ^ (1UL << bit);
}
inline double i2f(unsigned short s);
inline double i2f(unsigned short s)
{
	return s / 32768.l - 1.l;
}
inline unsigned short f2i(double f);
inline unsigned short f2i(double f)
{
	return (unsigned short)(f * 32768.l + 32768.l);
}
void HCM0_next(HCM0 *unit, int inNumSamples)
{
	float *out = ZOUT(0);
	float freq = ZIN0(0);
	double r = ZIN0(1);
	double g = ZIN0(2);
	unsigned short x[GENEBIT];
	double reciprocal = 1.l / GENEBIT;
	
	for (int i=0; i<GENEBIT; i++) x[i] = unit->x[i];
	float counter = unit->counter;
	
	float spc;
	if(freq < SAMPLERATE)
		spc = SAMPLERATE / sc_max(freq, 0.001f);
	else spc = 1.f;
	
	double sum = 0;
	double tmp;
	
	LOOP(inNumSamples,
		 if(counter >= spc){
			 counter -= spc;
			 for (int i=0; i<GENEBIT; i++){
				 for (int j=0; j<GENEBIT; j++) sum += logist(r, i2f(flip(x[j], i)));
				 tmp = (1.l - g) * logist(r, i2f(x[i])) + g * reciprocal * sum;
				 x[i] = f2i(tmp);
			 }
		 }
		 counter++;
		 ZXP(out) = i2f(x[4]);
	)
	for(int i=0; i<GENEBIT; i++) unit->x[i] = x[i];
	unit->counter = counter;
}

void HCM0_Ctor(HCM0 *unit)
{
	SETCALC(HCM0_next);
	for (int i=0; i<GENEBIT; i++) unit->x[i] = 1;
	unit->counter = 0.f;
	HCM0_next(unit, 1);
}

void HCM1_next(HCM1 *unit, int inNumSamples)
{
	float *out = ZOUT(0);
	float freq = ZIN0(0);
	double r = ZIN0(1);
	double g = ZIN0(2);
	unsigned short x[GENEBIT];
	double reciprocal = 1.l / GENEBIT;
	
	for (int i=0; i<GENEBIT; i++) 
		x[i] = unit->x[i];
	double xm = unit->xm;
	double frac = unit->frac;
	float counter = unit->counter;
	
	float spc;
	float slope;
	if(freq < SAMPLERATE){
		spc = SAMPLERATE / sc_max(freq, 0.001f);
		slope = 1.f / spc;
	} 
	else spc = slope = 1.f;
	
	double dx = i2f(x[5]) - xm;

	double sum = 0;
	double tmp;
	
	LOOP(inNumSamples,
		 if(counter >= spc){
			 counter -= spc;
			 frac = 0.l;
			 xm = i2f(x[5]);
			 for (int i=0; i<GENEBIT; i++){
				 for (int j=0; j<GENEBIT; j++)
					 sum += logist(r, i2f(flip(x[j], i)));
				 tmp = (1.l - g) * logist(r, i2f(x[i])) + g * reciprocal * sum;
				 x[i] = f2i(tmp);
			 }
			 dx = i2f(x[5]) - xm;
		 }
		 counter++;
		 ZXP(out) = xm + dx * frac;
		 frac += slope;
	)
	for(int i=0; i<GENEBIT; i++) unit->x[i] = x[i];
	unit->xm = xm;
	unit->frac = frac;
	unit->counter = counter;
}

void HCM1_Ctor(HCM1 *unit)
{
	SETCALC(HCM1_next);
	for (int i=0; i<GENEBIT; i++) unit->x[i] = 1;
	unit->counter = 0.f;
	unit->frac = 0.l;
	HCM1_next(unit, 1);
}

void HCM3_next(HCM3 *unit, int inNumSamples)
{
	float *out = ZOUT(0);
	float freq = ZIN0(0);
	double r = ZIN0(1);
	double g = ZIN0(2);
	unsigned short x[GENEBIT];
	double reciprocal = 1.l / GENEBIT;

	float counter = unit->counter;
	double frac = unit->frac;
	double xm1 = unit->xm1;
	double xm2 = unit->xm2;
	double xm3 = unit->xm3;
	double c0 = unit->c0;
	double c1 = unit->c1;
	double c2 = unit->c2;
	double c3 = unit->c3;

	for (int i=0; i<GENEBIT; i++) 
		x[i] = unit->x[i];
	
	float spc;
	double slope;
	if(freq < SAMPLERATE){
		spc = SAMPLERATE / sc_max(freq, 0.001f);
		slope = 1.l / spc;
	} 
	else spc = slope = 1.l;
	
	double sum = 0;
	double tmp;
	
	LOOP(inNumSamples,
		 if(counter >= spc){
			 counter -= spc;
			 frac = 0.l;
			 xm3 = xm2;
			 xm2 = xm1;
			 xm1 = i2f(x[4]);
			 for (int i=0; i<GENEBIT; i++){
				 for (int j=0; j<GENEBIT; j++)
					 sum += logist(r, i2f(flip(x[j], i)));
				 tmp = (1.l - g) * logist(r, i2f(x[i])) + g * reciprocal * sum;
				 x[i] = f2i(tmp);
			 }
			 ipol3Coef(xm3, xm2, xm1, i2f(x[4]), c0, c1, c2, c3);
		 }
		 counter++;
		 ZXP(out) = ipol3(frac, c0, c1, c2, c3);
		 frac += slope;
	)
	unit->counter = counter;
	unit->frac = frac;
	unit->xm1 = xm1;
	unit->xm2 = xm2;
	unit->xm3 = xm3;
	unit->c0 = c0;
	unit->c1 = c1;
	unit->c2 = c2;
	unit->c3 = c3;
	for(int i=0; i<GENEBIT; i++) unit->x[i] = x[i];
}

void HCM3_Ctor(HCM3 *unit)
{
	SETCALC(HCM3_next);
	for (int i=0; i<GENEBIT; i++) unit->x[i] = 1;
	unit->xm1 = unit->xm2 = unit->xm3 = unit->x[4];
	unit->c0 = unit->c1 = unit->c2 = unit->c3 = 0.l;
	unit->counter = 0.f;
	unit->frac = 0.l;
	HCM3_next(unit, 1);
}

void TLogist_next(TLogist *unit, int inNumSamples)
{	
	float trig = ZIN0(2);
	if (trig > 0.f && unit->trig <= 0.f) {
		double r = ZIN0(0);
		ZOUT0(0) = unit->x = 1.f - r * unit->x * unit->x;
	} else {
		ZOUT0(0) = unit->x;
	}
	unit->trig = trig;
}

void TLogist_Ctor(TLogist *unit)
{	
	double r = ZIN0(0);
	ZOUT0(0) = unit->x = ZIN0(1);
	SETCALC(TLogist_next);
	unit->trig = ZIN0(2);
}

PluginLoad(NL)
{
	ft = inTable;
	DefineSimpleUnit(Logist0);
	DefineSimpleUnit(Logist1);
	DefineSimpleUnit(Logist3);
	DefineSimpleUnit(Nagumo);
	DefineSimpleUnit(FIS);
	DefineSimpleUnit(CML0);
	DefineSimpleUnit(CML1);
	DefineSimpleUnit(CML3);
	DefineSimpleUnit(GCM0);
	DefineSimpleUnit(GCM1);
	DefineSimpleUnit(GCM3);
	DefineSimpleUnit(HCM0);
	DefineSimpleUnit(HCM1);
	DefineSimpleUnit(HCM3);
	DefineSimpleUnit(TLogist);
}
