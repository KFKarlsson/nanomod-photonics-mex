//Copyright (C) 2016-2018 Fredrik Karlsson (karlsson.fredrik@gmail.com)
//All rights reserved
//
//Modified Drude-code 20160502
//Modified live play 20170326
//Modified poynting box data 20170627
//Modified added active sides of poynting box (six numbers 0 or 1) after the PYBN paramaters

#include "mex.h"  /* Always include this */
#include <stdio.h>
#include <math.h>

#define SRCN	11 //Number or source parameters: i, j, k, field type:(0, 1, 2, 3, 4, or 5), point type: (0 single or 1 cluster), nx, ny, nz, phase-x, -y, -z
#define PRBN	5 //Number or probe parameters: i, j, k and field type:(0, 1, 2, 3, 4, or 5), tfsf type: (0:DF, 1:SF, 2:TF)
#define PYBN	11 //Number or poynting box parameters: box center position i, j, k, boxSizei, boxSizej, boxSizek, pybData size, startIteration, tfsf type: (0:DF, 1:SF, 2:TF), interleave, cnt
#define DFTN	10//Number or DFT parameters: frequency dimension (1:y, 2:x, 3:z) plane index field type sizex sizey, startIteration, tfsf type: (0:DF, 1:SF, 2:TF), interleave, cnt
#define EX		0 //Source/probe/dft types
#define EY		1
#define EZ		2
#define HX		3
#define HY		4
#define HZ		5
#define PLW		6
#define EDP		10	//Electric dipole
#define MDP		11	//Magentic dipole

#define TE		1
#define TM		2

#define SINGLEPOINT		0
#define CLUSTERPOINT	1


#define _ijkEx	(i)+(j)*Sy+(k)*Sy*sx	//Ex big in Y&Z  //i&k>0 
#define _ijkEy	(i)+(j)*sy+(k)*sy*Sx	//Ey big in X&Z  //j&k>0 
#define _ijkEz	(i)+(j)*Sy+(k)*Sy*Sx	//Ez big in X&Y	 //i&j>0 
#define _ijkHx	(i)+(j)*sy+(k)*sy*Sx
#define _ijkHy	(i)+(j)*Sy+(k)*Sy*sx
#define _ijkHz	(i)+(j)*sy+(k)*sy*sx

#define _dEx_dz	(i)+(j)*Sy+(k+1)*Sy*sx
#define _dEy_dz	(i)+(j)*sy+(k+1)*sy*Sx
#define _dEz_dy	(i+1)+(j)*Sy+(k)*Sy*Sx
#define _dEx_dy	(i+1)+(j)*Sy+(k)*Sy*sx
#define _dEy_dx	(i)+(j+1)*sy+(k)*sy*Sx
#define _dEz_dx	(i)+(j+1)*Sy+(k)*Sy*Sx

#define _dHx_dz	(i)+(j)*sy+(k-1)*sy*Sx
#define _dHy_dz	(i)+(j)*Sy+(k-1)*Sy*sx
#define _dHx_dy	(i-1)+(j)*sy+(k)*sy*Sx
#define _dHz_dy	(i-1)+(j)*sy+(k)*sy*sx
#define _dHz_dx	(i)+(j-1)*sy+(k)*sy*sx
#define _dHy_dx	(i)+(j-1)*Sy+(k)*Sy*sx


#define iMin 0
#define iMax 1
#define jMin 2
#define jMax 3
#define kMin 4
#define kMax 5

#define DF 0
#define SF 1
#define TF 2

#define CW 0
#define PULSE 1

#define SYM_FULL 0

#define BC_PML 0
#define BC_CYC 1

#define BCi 0		//Boundary type, PML or CYCLIC
#define BCj 1
#define BCk 2
#define SYMi 3		//Type EVEN or ODD
#define SYMj 4
#define SYMk 5
#define SYMy 6		//Domain size, even or odd
#define SYMx 7
#define SYMz 8

#define pi 3.141592653589793

int max ( int a, int b ) { return a > b ? a : b; }
int min ( int a, int b ) { return a > b ? b : a; }
double sign ( double a ) { return a >= 0 ? 1.0 : -1.0; }

struct SplitField {
	double *x;
	double *y;
	double *z;
	double *xy;
	double *xz;
	double *yx;
	double *yz;
	double *zx;
	double *zy;
};

struct Map1D {
	int	src_type;
	double amplitude;
	double src_f;
	double src_df;
	int MX;
	int MY;
	double nx;
	double ny;
	double dr;
	double signX;
	double signY;	
	double nr;
	double nrc;
	double nrc0;
	double* map0;
	double* map1;
	double* map2;
	double* mapX;
	double* mapY;
	int nMin;
	int nMax;
	int N;
};

double sourceSignal(double src_f, double src_df, double src_amplitude, double src_phase, double src_type, double x, double vt, double nr) {
	double signal;
	double sigma = 0.4/(nr*src_df);
	double x0 = 6*sigma;
	double kVec = 2*pi*src_f*nr;
	double X = (x0+x) - vt;
	if (src_type == PULSE) signal = src_amplitude*sin(kVec*X+src_phase)*exp(-(X*X) / (2*sigma*sigma));
	else if (src_type == CW) signal = src_amplitude*sin(kVec*X+src_phase)/pow(1 + exp(X/sigma), 2);
	return signal;
}


void addValueToFieldAtClusterPoint(int i, int j, int k, double *F, int fieldType, int sx, int sy, int sz, double value) {
    int Sx = sx+1;
    int Sy = sy+1;
    
	switch (fieldType) {  //Works only in interior, perhaps not at symmetry boundary.
		case EX:
			F[(i)+(j)*Sy+(k)*Sy*sx] += value/4;
			F[(i+1)+(j)*Sy+(k)*Sy*sx] += value/4;
			F[(i)+(j)*Sy+(k+1)*Sy*sx] += value/4;
			F[(i+1)+(j)*Sy+(k+1)*Sy*sx] += value/4;
			break;
		case EY:
			F[(i)+(j)*sy+(k)*sy*Sx] += value/4;
			F[(i)+(j+1)*sy+(k)*sy*Sx] += value/4;
			F[(i)+(j)*sy+(k+1)*sy*Sx] += value/4;
			F[(i)+(j+1)*sy+(k+1)*sy*Sx] += value/4;
			break;
		case EZ:
			F[(i)+(j)*Sy+(k)*Sy*Sx] += value/4;
			F[(i+1)+(j)*Sy+(k)*Sy*Sx] += value/4;
			F[(i)+(j+1)*Sy+(k)*Sy*Sx] += value/4;
			F[(i+1)+(j+1)*Sy+(k)*Sy*Sx] += value/4;
			break;
		case HX:
			F[(i)+(j)*sy+(k)*sy*Sx] += value/2;
			F[(i)+(j+1)*sy+(k)*sy*Sx] += value/2;
			break;
		case HY:
			F[(i)+(j)*Sy+(k)*Sy*sx] += value/2;
			F[(i+1)+(j)*Sy+(k)*Sy*sx] += value/2;
			break;
		case HZ:
			F[(i)+(j)*sy+(k)*sy*sx] += value/2;
			F[(i)+(j)*sy+(k+1)*sy*sx] += value/2;	
			break;
	}
	return;
}

void addValueToFieldAtSinglePoint(int i, int j, int k, double *F, int fieldType, int sx, int sy, int sz, double value) {
    int Sx = sx+1;
    int Sy = sy+1;
    
	switch (fieldType) {  //Works only in interior, perhaps not at symmetry boundary.
		case EX:
			F[(i)+(j)*Sy+(k)*Sy*sx] += value;
			break;
		case EY:
			F[(i)+(j)*sy+(k)*sy*Sx] += value;
			break;
		case EZ:
			F[(i)+(j)*Sy+(k)*Sy*Sx] += value;
			break;
		case HX:
			F[(i)+(j)*sy+(k)*sy*Sx] += value;
			break;
		case HY:
			F[(i)+(j)*Sy+(k)*Sy*sx] += value;
			break;
		case HZ:
			F[(i)+(j)*sy+(k)*sy*sx] += value;
			break;
	}
	return;
}

double getRawValueAtClusterPoint(int i, int j, int k, double *F, int fieldType, int sx, int sy, int sz) {
	double value;
    int Sx = sx+1;
    int Sy = sy+1;
    
	switch (fieldType) {
		case EX:
			value = (F[(i)+(j)*Sy+(k)*Sy*sx] + F[(i+1)+(j)*Sy+(k)*Sy*sx] + F[(i)+(j)*Sy+(k+1)*Sy*sx] + F[(i+1)+(j)*Sy+(k+1)*Sy*sx])/4;
			break;
		case EY:
			value = (F[(i)+(j)*sy+(k)*sy*Sx] + F[(i)+(j+1)*sy+(k)*sy*Sx] + F[(i)+(j)*sy+(k+1)*sy*Sx] + F[(i)+(j+1)*sy+(k+1)*sy*Sx])/4;
			break;
		case EZ:
			value = (F[(i)+(j)*Sy+(k)*Sy*Sx] + F[(i+1)+(j)*Sy+(k)*Sy*Sx] + F[(i)+(j+1)*Sy+(k)*Sy*Sx] + F[(i+1)+(j+1)*Sy+(k)*Sy*Sx])/4;
			break;
		case HX:
			value = (F[(i)+(j)*sy+(k)*sy*Sx] + F[(i)+(j+1)*sy+(k)*sy*Sx])/2;
			break;
		case HY:
			value = (F[(i)+(j)*Sy+(k)*Sy*sx] + F[(i+1)+(j)*Sy+(k)*Sy*sx])/2;
			break;
		case HZ:
			value = (F[(i)+(j)*sy+(k)*sy*sx] + F[(i)+(j)*sy+(k+1)*sy*sx])/2;
			break;
	}
	return value;
}

double getRawValueAtSinglePoint(int i, int j, int k, double *F, int fieldType, int sx, int sy, int sz) {
	double value = 0;
    int Sx = sx+1;
    int Sy = sy+1;

    
	switch (fieldType) {
		case EX:
			if (j<sx) value = F[(i)+(j)*Sy+(k)*Sy*sx];
			break;
		case EY:
			if (i<sy) value = F[(i)+(j)*sy+(k)*sy*Sx];
			break;
		case EZ:
			if (k<sz) value = F[(i)+(j)*Sy+(k)*Sy*Sx];
			break;
		case HX:
			if ((i<sy)&&(k<sz)) value = F[(i)+(j)*sy+(k)*sy*Sx];
			break;
		case HY:
			if ((j<sx)&&(k<sz)) value = F[(i)+(j)*Sy+(k)*Sy*sx];
			break;
		case HZ:
			if ((i<sy)&&(j<sx)) value = F[(i)+(j)*sy+(k)*sy*sx];
			break;
	}
	return value;
}


double getValueAtSinglePointModifyChecked(int i, int j, int k, double *F, int fieldType, int sx, int sy, int sz, double *tfsfBox, double *map0, double *mapX, double *mapY, double nx, double ny, double dr, int signX, int signY, int tfsf) {
	double value = 0;
    int Sx = sx+1;
    int Sy = sy+1;
	
	int TFSFy = tfsfBox[iMin];
	int TFSFY = tfsfBox[iMax];
	int TFSFx = tfsfBox[jMin];
	int TFSFX = tfsfBox[jMax];
	int TFSFz = tfsfBox[kMin];
	int TFSFZ = tfsfBox[kMax];
	int x0, y0;
    
	double mod = 0;
	double TFSF = 0;		//SF if -1, TF if 1, DF if 0
	if (tfsf==SF) TFSF = -1; else if (tfsf==TF) TFSF = 1;
	
	switch (fieldType) {
		int withinTF;
		case EX:  //TE		
			withinTF = ((i>=TFSFy)&&(i<=TFSFY)&&(j>=TFSFx)&&(j<TFSFX)&&(k>TFSFz)&&(k<=TFSFZ))?1:-1; //1 inside TF -1 outside TF
			if (!(withinTF+TFSF)) {		//sum zero if inside TF and SF is requested, or if inside SF and TF is requested 
				if (signX>0) x0 = 1+(j-TFSFx); else x0 = TFSFX-j;
				if (signY>0) y0 = 1+(i-TFSFy); else y0 = 1+TFSFY-i;
				mod = signY*TFSF*mapX[10+(int)round((x0*nx+y0*ny) / dr)];
			} 
			if (j<sx) value = F[(i)+(j)*Sy+(k)*Sy*sx] - mod;
			break;
			
		case EY:  //TE
			withinTF = ((i>=TFSFy)&&(i<TFSFY)&&(j>=TFSFx)&&(j<=TFSFX)&&(k>TFSFz)&&(k<=TFSFZ))?1:-1;
			if (!(withinTF+TFSF)) { //sum zero if inside TF and SF is requested, or if inside SF and TF is requested
				if (signX>0) x0 = 1+(j-TFSFx); else x0 = 1+TFSFX-j;
				if (signY>0) y0 = 1+(i-TFSFy); else y0 = TFSFY-i;
				mod = signX*TFSF*mapY[10+(int)round((x0*nx+y0*ny) / dr)];
			}
			if (i<sy) value = F[(i)+(j)*sy+(k)*sy*Sx] - mod;
			break;
			
		case EZ:  //TM
			withinTF = ((i>=TFSFy)&&(i<=TFSFY)&&(j>=TFSFy)&&(j<=TFSFX)&&(k>TFSFz)&&(k<=TFSFZ))?1:-1;
			if (!(withinTF+TFSF)) { //sum zero if inside TF and SF is requested, or if inside SF and TF is requested
				if (signX>0) x0 = 1+(j-TFSFx); else x0 = 1+TFSFX-j;
				if (signY>0) y0 = 1+(i-TFSFy); else y0 = 1+TFSFY-i;
				mod = TFSF*map0[10+(int)round((x0*nx+y0*ny) / dr)];
			}
			if (k<sz) value = F[(i)+(j)*Sy+(k)*Sy*Sx] + mod;
			break;
			
		case HX:  //TM
			withinTF = ((i>=TFSFy)&&(i<TFSFY)&&(j>=TFSFx)&&(j<=TFSFX)&&(k>TFSFz)&&(k<=TFSFZ))?1:-1; //1 inside TF -1 outside TF
			if (!(withinTF+TFSF)) {		//sum zero if inside TF and SF is requested, or if inside SF and TF is requested 
				if (signX>0) x0 = 1+(j-TFSFx); else x0 = 1+TFSFX-j;
				if (signY>0) y0 = 1+(i-TFSFy); else y0 = TFSFY-i;
				mod = signY*TFSF*mapX[10+(int)round((x0*nx+y0*ny) / dr)];
			} 	
			if ((i<sy)&&(k<sz)) value = F[(i)+(j)*sy+(k)*sy*Sx]+mod;
			
			break;
		case HY:  //TM
			withinTF = ((i>=TFSFy)&&(i<=TFSFY)&&(j>=TFSFx)&&(j<TFSFX)&&(k>TFSFz)&&(k<=TFSFZ))?1:-1;
			if (!(withinTF+TFSF)) { //sum zero if inside TF and SF is requested, or if inside SF and TF is requested
				if (signX>0) x0 = 1+(j-TFSFx); else x0 = TFSFX-j;
				if (signY>0) y0 = 1+(i-TFSFy); else y0 = 1+TFSFY-i;
				mod =  signX*TFSF*mapY[10+(int)round((x0*nx+y0*ny) / dr)];
			}
			if ((j<sx)&&(k<sz)) value = F[(i)+(j)*Sy+(k)*Sy*sx]+mod;
			break;
			
		case HZ:	//TE
			withinTF = ((i>=TFSFy)&&(i<TFSFY)&&(j>=TFSFy)&&(j<TFSFX)&&(k>TFSFz)&&(k<=TFSFZ))?1:-1;
			if (!(withinTF+TFSF)) { //sum zero if inside TF and SF is requested, or if inside SElinF and TF is requested
				if (signX>0) x0 = 1+(j-TFSFx); else x0 = TFSFX-j;
				if (signY>0) y0 = 1+(i-TFSFy); else y0 = TFSFY-i;
				mod = TFSF*map0[10+(int)round((x0*nx+y0*ny) / dr)];
			}
			if ((i<sy)&&(j<sx)) value = F[(i)+(j)*sy+(k)*sy*sx] + mod;
			break;
	}
	
	return value;
}


double getValueAtSinglePointModifyWithoutCheck(int i, int j, int k, double *F, int fieldType, int sx, int sy, int sz, double *tfsfBox, double *map0, double *mapX, double *mapY, double nx, double ny, double dr, int signX, int signY, int tfsf) {
	double value = 0;
    int Sx = sx+1;
    int Sy = sy+1;
	
	int TFSFy = tfsfBox[iMin];
	int TFSFY = tfsfBox[iMax];
	int TFSFx = tfsfBox[jMin];
	int TFSFX = tfsfBox[jMax];
	int x0, y0;
    
	double mod = 0;
	double TFSF = 0;		//SF if -1, TF if 1, DF if 0
	if (tfsf==SF) TFSF = -1; else if (tfsf==TF) TFSF = 1;
	
	switch (fieldType) {
		case EX:  //TE		
			if (signX>0) x0 = 1+(j-TFSFx); else x0 = TFSFX-j;
			if (signY>0) y0 = 1+(i-TFSFy); else y0 = 1+TFSFY-i;
			mod = signY*TFSF*mapX[10+(int)round((x0*nx+y0*ny) / dr)];
			if (j<sx) value = F[(i)+(j)*Sy+(k)*Sy*sx] - mod;
			break;
			
		case EY:  //TE
			if (signX>0) x0 = 1+(j-TFSFx); else x0 = 1+TFSFX-j;
			if (signY>0) y0 = 1+(i-TFSFy); else y0 = TFSFY-i;
			mod = signX*TFSF*mapY[10+(int)round((x0*nx+y0*ny) / dr)];
			if (i<sy) value = F[(i)+(j)*sy+(k)*sy*Sx] - mod;
			break;
			
		case EZ:  //TM
			if (signX>0) x0 = 1+(j-TFSFx); else x0 = 1+TFSFX-j;
			if (signY>0) y0 = 1+(i-TFSFy); else y0 = 1+TFSFY-i;
			mod = TFSF*map0[10+(int)round((x0*nx+y0*ny) / dr)];
			if (k<sz) value = F[(i)+(j)*Sy+(k)*Sy*Sx] + mod;
			break;
			
		case HX:  //TM
			if (signX>0) x0 = 1+(j-TFSFx); else x0 = 1+TFSFX-j;
			if (signY>0) y0 = 1+(i-TFSFy); else y0 = TFSFY-i;
			mod = signY*TFSF*mapX[10+(int)round((x0*nx+y0*ny) / dr)];
			if ((i<sy)&&(k<sz)) value = F[(i)+(j)*sy+(k)*sy*Sx]+mod;
			
			break;
		case HY:  //TM
			if (signX>0) x0 = 1+(j-TFSFx); else x0 = TFSFX-j;
			if (signY>0) y0 = 1+(i-TFSFy); else y0 = 1+TFSFY-i;
			mod =  signX*TFSF*mapY[10+(int)round((x0*nx+y0*ny) / dr)];
			if ((j<sx)&&(k<sz)) value = F[(i)+(j)*Sy+(k)*Sy*sx]+mod;
			break;
			
		case HZ:	//TE
			if (signX>0) x0 = 1+(j-TFSFx); else x0 = TFSFX-j;
			if (signY>0) y0 = 1+(i-TFSFy); else y0 = TFSFY-i;
			mod = TFSF*map0[10+(int)round((x0*nx+y0*ny) / dr)];
			if ((i<sy)&&(j<sx)) value = F[(i)+(j)*sy+(k)*sy*sx]+mod;
			break;
	}
	
	return value;
}

double getValueAtClusterPointModifyWithoutCheck(int i, int j, int k, double *F, int fieldType, int sx, int sy, int sz, double *tfsfBox, double *map0, double *mapX, double *mapY, double nx, double ny, double dr, int signX, int signY, int tfsf) {
	double value;
	
	value = getValueAtSinglePointModifyWithoutCheck(i, j, k, F, fieldType, sx, sy, sz, tfsfBox, map0, mapX, mapY, nx, ny, dr, signX, signY, tfsf);

	switch (fieldType) {
		case EX:
			//value = (F[(i)+(j)*Sy+(k)*Sy*sx] + F[(i+1)+(j)*Sy+(k)*Sy*sx] + F[(i)+(j)*Sy+(k+1)*Sy*sx] + F[(i+1)+(j)*Sy+(k+1)*Sy*sx])/4;
			value += getValueAtSinglePointModifyWithoutCheck(i+1, j, k  , F, fieldType, sx, sy, sz, tfsfBox, map0, mapX, mapY, nx, ny, dr, signX, signY, tfsf);
			value += getValueAtSinglePointModifyWithoutCheck(  i, j, k+1, F, fieldType, sx, sy, sz, tfsfBox, map0, mapX, mapY, nx, ny, dr, signX, signY, tfsf);
			value += getValueAtSinglePointModifyWithoutCheck(i+1, j, k+1, F, fieldType, sx, sy, sz, tfsfBox, map0, mapX, mapY, nx, ny, dr, signX, signY, tfsf);
			value /= 4;
			break;
		case EY:
			//value = (F[(i)+(j)*sy+(k)*sy*Sx] + F[(i)+(j+1)*sy+(k)*sy*Sx] + F[(i)+(j)*sy+(k+1)*sy*Sx] + F[(i)+(j+1)*sy+(k+1)*sy*Sx])/4;
			value += getValueAtSinglePointModifyWithoutCheck(i, j+1, k  , F, fieldType, sx, sy, sz, tfsfBox, map0, mapX, mapY, nx, ny, dr, signX, signY, tfsf);
			value += getValueAtSinglePointModifyWithoutCheck(i,   j, k+1, F, fieldType, sx, sy, sz, tfsfBox, map0, mapX, mapY, nx, ny, dr, signX, signY, tfsf);
			value += getValueAtSinglePointModifyWithoutCheck(i, j+1, k+1, F, fieldType, sx, sy, sz, tfsfBox, map0, mapX, mapY, nx, ny, dr, signX, signY, tfsf);
			value /= 4;
			break;
		case EZ:
			//value = (F[(i)+(j)*Sy+(k)*Sy*Sx] + F[(i+1)+(j)*Sy+(k)*Sy*Sx] + F[(i)+(j+1)*Sy+(k)*Sy*Sx] + F[(i+1)+(j+1)*Sy+(k)*Sy*Sx])/4;
			value += getValueAtSinglePointModifyWithoutCheck(i+1,   j, k, F, fieldType, sx, sy, sz, tfsfBox, map0, mapX, mapY, nx, ny, dr, signX, signY, tfsf);
			value += getValueAtSinglePointModifyWithoutCheck(  i, j+1, k, F, fieldType, sx, sy, sz, tfsfBox, map0, mapX, mapY, nx, ny, dr, signX, signY, tfsf);
			value += getValueAtSinglePointModifyWithoutCheck(i+1, j+1, k, F, fieldType, sx, sy, sz, tfsfBox, map0, mapX, mapY, nx, ny, dr, signX, signY, tfsf);
			value /= 4;
			break;
		case HX:
			//value = (F[(i)+(j)*sy+(k)*sy*Sx] + F[(i)+(j+1)*sy+(k)*sy*Sx])/2;
			value += getValueAtSinglePointModifyWithoutCheck(i, j+1, k, F, fieldType, sx, sy, sz, tfsfBox, map0, mapX, mapY, nx, ny, dr, signX, signY, tfsf);
			value /= 2;
			break;
		case HY:
			//value = (F[(i)+(j)*Sy+(k)*Sy*sx] + F[(i+1)+(j)*Sy+(k)*Sy*sx])/2;
			value += getValueAtSinglePointModifyWithoutCheck(i+1, j, k, F, fieldType, sx, sy, sz, tfsfBox, map0, mapX, mapY, nx, ny, dr, signX, signY, tfsf);
			value /= 2;
			break;
		case HZ:
			//value = (F[(i)+(j)*sy+(k)*sy*sx] + F[(i)+(j)*sy+(k+1)*sy*sx])/2;
			value += getValueAtSinglePointModifyWithoutCheck(i, j, k+1, F, fieldType, sx, sy, sz, tfsfBox, map0, mapX, mapY, nx, ny, dr, signX, signY, tfsf);
			value /= 2;
			break;
	}
	return value;
}

double getValueAtClusterPointModifyCheckAll(int i, int j, int k, double *F, int fieldType, int sx, int sy, int sz, double *tfsfBox, double *map0, double *mapX, double *mapY, double nx, double ny, double dr, int signX, int signY, int tfsf) {
	double value;
	
	value = getValueAtSinglePointModifyChecked(i, j, k, F, fieldType, sx, sy, sz, tfsfBox, map0, mapX, mapY, nx, ny, dr, signX, signY, tfsf);

	switch (fieldType) {
		case EX:
			//value = (F[(i)+(j)*Sy+(k)*Sy*sx] + F[(i+1)+(j)*Sy+(k)*Sy*sx] + F[(i)+(j)*Sy+(k+1)*Sy*sx] + F[(i+1)+(j)*Sy+(k+1)*Sy*sx])/4;
			value += getValueAtSinglePointModifyChecked(i+1, j, k  , F, fieldType, sx, sy, sz, tfsfBox, map0, mapX, mapY, nx, ny, dr, signX, signY, tfsf);
			value += getValueAtSinglePointModifyChecked(  i, j, k+1, F, fieldType, sx, sy, sz, tfsfBox, map0, mapX, mapY, nx, ny, dr, signX, signY, tfsf);
			value += getValueAtSinglePointModifyChecked(i+1, j, k+1, F, fieldType, sx, sy, sz, tfsfBox, map0, mapX, mapY, nx, ny, dr, signX, signY, tfsf);
			value /= 4;
			break;
		case EY:
			//value = (F[(i)+(j)*sy+(k)*sy*Sx] + F[(i)+(j+1)*sy+(k)*sy*Sx] + F[(i)+(j)*sy+(k+1)*sy*Sx] + F[(i)+(j+1)*sy+(k+1)*sy*Sx])/4;
			value += getValueAtSinglePointModifyChecked(i, j+1, k  , F, fieldType, sx, sy, sz, tfsfBox, map0, mapX, mapY, nx, ny, dr, signX, signY, tfsf);
			value += getValueAtSinglePointModifyChecked(i,   j, k+1, F, fieldType, sx, sy, sz, tfsfBox, map0, mapX, mapY, nx, ny, dr, signX, signY, tfsf);
			value += getValueAtSinglePointModifyChecked(i, j+1, k+1, F, fieldType, sx, sy, sz, tfsfBox, map0, mapX, mapY, nx, ny, dr, signX, signY, tfsf);
			value /= 4;
			break;
		case EZ:
			//value = (F[(i)+(j)*Sy+(k)*Sy*Sx] + F[(i+1)+(j)*Sy+(k)*Sy*Sx] + F[(i)+(j+1)*Sy+(k)*Sy*Sx] + F[(i+1)+(j+1)*Sy+(k)*Sy*Sx])/4;
			value += getValueAtSinglePointModifyChecked(i+1,   j, k, F, fieldType, sx, sy, sz, tfsfBox, map0, mapX, mapY, nx, ny, dr, signX, signY, tfsf);
			value += getValueAtSinglePointModifyChecked(  i, j+1, k, F, fieldType, sx, sy, sz, tfsfBox, map0, mapX, mapY, nx, ny, dr, signX, signY, tfsf);
			value += getValueAtSinglePointModifyChecked(i+1, j+1, k, F, fieldType, sx, sy, sz, tfsfBox, map0, mapX, mapY, nx, ny, dr, signX, signY, tfsf);
			value /= 4;
			break;
		case HX:
			//value = (F[(i)+(j)*sy+(k)*sy*Sx] + F[(i)+(j+1)*sy+(k)*sy*Sx])/2;
			value += getValueAtSinglePointModifyChecked(i, j+1, k, F, fieldType, sx, sy, sz, tfsfBox, map0, mapX, mapY, nx, ny, dr, signX, signY, tfsf);
			value /= 2;
			break;
		case HY:
			//value = (F[(i)+(j)*Sy+(k)*Sy*sx] + F[(i+1)+(j)*Sy+(k)*Sy*sx])/2;
			value += getValueAtSinglePointModifyChecked(i+1, j, k, F, fieldType, sx, sy, sz, tfsfBox, map0, mapX, mapY, nx, ny, dr, signX, signY, tfsf);
			value /= 2;
			break;
		case HZ:
			//value = (F[(i)+(j)*sy+(k)*sy*sx] + F[(i)+(j)*sy+(k+1)*sy*sx])/2;
			value += getValueAtSinglePointModifyChecked(i, j, k+1, F, fieldType, sx, sy, sz, tfsfBox, map0, mapX, mapY, nx, ny, dr, signX, signY, tfsf);
			value /= 2;
			break;
	}
	return value;
}


double getValueAtClusterPoint(int i, int j, int k, double *F, int fieldType, int sx, int sy, int sz, double *tfsfBox, double *map0, double *mapX, double *mapY, double nx, double ny, double dr, int signX, int signY, int tfsf, double withinTFflag) {
	switch ((int)withinTFflag) {
		case 2:
		return getRawValueAtClusterPoint(i, j, k, F, fieldType, sx, sy, sz);
		case 0:
		return getValueAtClusterPointModifyWithoutCheck(i, j, k, F, fieldType, sx, sy, sz, tfsfBox, map0, mapX, mapY, nx, ny, dr, signX, signY, tfsf);
		default:
		return getValueAtClusterPointModifyCheckAll(i, j, k, F, fieldType, sx, sy, sz, tfsfBox, map0, mapX, mapY, nx, ny, dr, signX, signY, tfsf);
	}
}


void updateTFSF_ETE(double *Ex, double *Ey, double *Exy, double *Eyx, int *bc, double *tfsfBox, double *map0, double nx, double ny, double dr, int signX, int signY, int sx, int sy, int sz, double dt, double nr) {
	
	int i, j, k, ijk, x0, y0;
	double val;
	
	int TFSFy = tfsfBox[iMin];
	int TFSFY = tfsfBox[iMax];
	int TFSFx = tfsfBox[jMin];
	int TFSFX = tfsfBox[jMax];
	int TFSFz = tfsfBox[kMin];
	int TFSFZ = tfsfBox[kMax];
	
    int Sx = sx+1;
	int Sy = sy+1;
	
	double A1 = dt/(nr*nr);
	
	for (k=(TFSFz+1); k<=min(TFSFZ, sz-1); k++) for (j=TFSFx; j<min(TFSFX, sx); j++) {
		i=TFSFy; ijk = _ijkEx;		
		if (signX>0) x0 = 1+(j-TFSFx); else x0 = TFSFX-j;
		if (signY>0) y0 = 0;  else y0 = 1+TFSFY-TFSFy;
		val = A1*map0[10+(int)round((x0*nx+y0*ny) / dr)];
		Exy[ijk] -= val; Ex[ijk]  -= val;
		
		if (bc[SYMi]==SYM_FULL) {
			i=TFSFY; ijk = _ijkEx;	
			if (signY>0) y0 = 1+TFSFY-TFSFy; else y0 = 0;
			val = A1*map0[10+(int)round((x0*nx+y0*ny) / dr)];
			Exy[ijk] += val; Ex[ijk] += val;
		}
	}
	
	for (k=(TFSFz+1); k<=min(TFSFZ, sz-1); k++) for (i=TFSFy; i<min(TFSFY,sy); i++) {
		j=TFSFx; ijk = _ijkEy;
		if (signY>0) y0 = 1+(i-TFSFy); else y0 = TFSFY-i;
		if (signX>0) x0 = 0;  else x0 = TFSFX-TFSFx+1;
		val = A1*map0[10+(int)round((x0*nx+y0*ny) / dr)];
		Eyx[ijk] += val; Ey[ijk] += val;

		if (bc[SYMj]==SYM_FULL) {
			j=TFSFX; ijk = _ijkEy;
			if (signX>0) x0 = TFSFX-TFSFx+1; else x0 = 0;
			val = A1*map0[10+(int)round((x0*nx+y0*ny) / dr)];
			Eyx[ijk] -= val; Ey[ijk] -= val;
		}
	}
}

void updateTFSF_ETM(double *Ex, double *Exz, double *Ey, double *Eyz, double *Ez, double *Ezx, double *Ezy, int *bc, double *tfsfBox, double *mapX, double *mapY, double nx, double ny, double dr, int signX, int signY, int sx, int sy, int sz, double dt, double nr) {

	int i, j, k, ijk, x0, y0;
	double val;
	
	int TFSFy = tfsfBox[iMin];
	int TFSFY = tfsfBox[iMax];
	int TFSFx = tfsfBox[jMin];
	int TFSFX = tfsfBox[jMax];
	int TFSFz = tfsfBox[kMin];
	int TFSFZ = tfsfBox[kMax];
	
    int Sx = sx+1;
	int Sy = sy+1;
	
	double A1 = dt/(nr*nr);
	
	
	for (j=TFSFx; j<min(TFSFX, sx); j++) for (i=TFSFy; i<=min(TFSFY, sy-1); i++) {
		k = TFSFz+1; ijk = _ijkEx;
		if (signX>0) x0 = 1+j-TFSFx; else x0 = TFSFX-j;
		if (signY>0) y0 = 1+i-TFSFy; else y0 = 1+TFSFY-i;
		val = signX*A1*mapY[10+(int)round((x0*nx+y0*ny) / dr)];
		Exz[ijk] += val; Ex[ijk] += val;
		
		if (bc[SYMk]==SYM_FULL) {
			k = TFSFZ+1; ijk = _ijkEx;
			val = signX*A1*mapY[10+(int)round((x0*nx+y0*ny) / dr)];
			Exz[ijk] -= val; Ex[ijk] -= val;
		}
	}
	
	for (j=TFSFx; j<=min(TFSFX, sx-1); j++) for (i=TFSFy; i<min(TFSFY, sy); i++) {
		k = TFSFz+1; ijk = _ijkEy;
		if (signX>0) x0 = 1+(j-TFSFx); else x0 = 1+TFSFX-j;
		if (signY>0) y0 = 1+i-TFSFy; else y0 = TFSFY-i;
		val = signY*A1*mapX[10+(int)round((x0*nx+y0*ny) / dr)];
		Eyz[ijk] -= val;  Ey[ijk] -= val;
		
		if (bc[SYMk]==SYM_FULL) {
			k = TFSFZ+1; ijk = _ijkEy;
			val = signY*A1*mapX[10+(int)round((x0*nx+y0*ny) / dr)];
			Eyz[ijk] += val; Ey[ijk] += val;
		}	
	}
							
	for (k=(TFSFz+1); k<=min(TFSFZ, sz-1); k++) for (i=TFSFy; i<=min(TFSFY, sy-1); i++) {
		j = TFSFx; ijk = _ijkEz;
		if (signY>0) y0 = 1+(i-TFSFy); else y0 = 1+TFSFY-i;
		if (signX>0) x0 = 0; else x0 = 1+TFSFX-TFSFx;
		val = signX*A1*mapY[10+(int)round((x0*nx+y0*ny) / dr)];
		Ezx[ijk] -= val; Ez[ijk] -= val;
		
		if (bc[SYMj]==SYM_FULL) {
			j=TFSFX; ijk = _ijkEz;
			if (signX>0) x0 = 1+TFSFX-TFSFx; else x0 = 0;
			val = signX*A1*mapY[10+(int)round((x0*nx+y0*ny) / dr)];
			Ezx[ijk] += val; Ez[ijk] += val;
		}
	}
	
	for (k=(TFSFz+1); k<=min(TFSFZ, sz-1); k++) for (j=TFSFx; j<=min(TFSFX, sx-1); j++) {
		i = TFSFy; ijk = _ijkEz;
		if (signX>0) x0 = 1+(j-TFSFx); else x0 = 1+TFSFX-j;
		if (signY>0) y0 = 0; else y0 = 1+(TFSFY-TFSFy);
		val = signY*A1*mapX[10+(int)round((x0*nx+y0*ny) / dr)];
		Ezy[ijk] += val; Ez[ijk] += val;
		
		if (bc[SYMi]==SYM_FULL) {
			i=TFSFY; ijk = _ijkEz;
			if (signY>0)  y0 = 1+(TFSFY-TFSFy); else y0 = 0;
			val = signY*A1*mapX[10+(int)round((x0*nx+y0*ny) / dr)];
			Ezy[ijk] -= val; Ez[ijk] -= val;
		}
	}
}

void updateTFSF_HTE(double *Hx, double *Hxz, double *Hy, double *Hyz, double *Hz, double *Hzx, double *Hzy, int *bc, double *tfsfBox, double *mapX, double *mapY, double nx, double ny, double dr, int signX, int signY, int sx, int sy, int sz, double dt) {

	int i, j, k, ijk, x0, y0;
	double val;
	
	int TFSFy = tfsfBox[iMin];
	int TFSFY = tfsfBox[iMax];
	int TFSFx = tfsfBox[jMin];
	int TFSFX = tfsfBox[jMax];
	int TFSFz = tfsfBox[kMin];
	int TFSFZ = tfsfBox[kMax];
	
    int Sx = sx+1;
	int Sy = sy+1;
	
	
	for (j=TFSFx; j<=min(TFSFX, sx); j++) for (i=TFSFy; i<min(TFSFY, sy); i++) {
		k = TFSFz; ijk = _ijkHx;
		if (signX>0) x0 = j-TFSFx; else x0 = TFSFX-j;
		if (signY>0) y0 = 1+i-TFSFy; else y0 = TFSFY-i;
		x0 += 1; 
		val = signX*dt*mapY[10+(int)round((x0*nx+y0*ny) / dr)];
		Hxz[ijk] += val; Hx[ijk] += val/2;
		
		if (bc[SYMk]==SYM_FULL) {
			k = TFSFZ; ijk = _ijkHx;
			val = signX*dt*mapY[10+(int)round((x0*nx+y0*ny) / dr)];
			Hxz[ijk] -= val; Hx[ijk] -= val/2;
		}
	}
	
	for (j=TFSFx; j<min(TFSFX, sx); j++) for (i=TFSFy; i<=min(TFSFY, sy); i++) {
		k = TFSFz; ijk = _ijkHy;
		if (signX>0) x0 = 1+j-TFSFx; else x0 = TFSFX-j;
		if (signY>0) y0 = i-TFSFy; else y0 = TFSFY-i;
		y0 += 1;
		val = signY*dt*mapX[10+(int)round((x0*nx+y0*ny) / dr)];
		Hyz[ijk] -= val;  Hy[ijk] -= val/2;
		
		if (bc[SYMk]==SYM_FULL) {
			k = TFSFZ; ijk = _ijkHy;
			val = signY*dt*mapX[10+(int)round((x0*nx+y0*ny) / dr)];
			Hyz[ijk] += val; Hy[ijk] += val/2;
		}	
	}
	
    						
    for (k=(TFSFz+1); k<=min(TFSFZ, sz); k++) for (i=TFSFy; i<min(TFSFY, sy); i++) {
		j=(TFSFx-1); ijk = _ijkHz;
		if (signY>0) y0 = 1+(i-TFSFy); else y0 = TFSFY-i;
		if (signX>0) x0 = 0; else x0 = TFSFX-TFSFx;
		x0 += 1; 
		val = signX*dt*mapY[10+(int)round((x0*nx+y0*ny) / dr)];
		Hzx[ijk] -= val; Hz[ijk] -= val/2;
		
		if (bc[SYMj]==SYM_FULL) {
			j=TFSFX; ijk = _ijkHz;
			if (signX>0) x0 = TFSFX-TFSFx; else x0 = 0;
			x0 += 1;
			val = signX*dt*mapY[10+(int)round((x0*nx+y0*ny) / dr)];
			Hzx[ijk] += val; Hz[ijk] += val/2;
		}
	}
	
	for (k=(TFSFz+1); k<=min(TFSFZ, sz); k++) for (j=TFSFx; j<min(TFSFX, sx); j++) {
		i = (TFSFy-1); ijk = _ijkHz;
		if (signX>0) x0 = 1+(j-TFSFx); else x0 = TFSFX-j;
		if (signY>0) y0 = 0; else y0 = (TFSFY-TFSFy);
		y0 += 1;
		val = signY*dt*mapX[10+(int)round((x0*nx+y0*ny) / dr)];
		Hzy[ijk] += val; Hz[ijk] += val/2;
		
		if (bc[SYMi]==SYM_FULL) {
			i=TFSFY; ijk = _ijkHz;
			if (signY>0)  y0 = (TFSFY-TFSFy); else y0 = 0;
			y0 += 1;
			val = signY*dt*mapX[10+(int)round((x0*nx+y0*ny) / dr)];
			Hzy[ijk] -= val; Hz[ijk] -= val/2;
		}
	}					
}

void updateTFSF_HTM(double *Hx, double *Hy, double *Hxy, double *Hyx, int *bc, double *tfsfBox, double *map0, double nx, double ny, double dr, int signX, int signY, int sx, int sy, int sz) {
	
	int i, j, k, ijk, x0, y0;
	double val;
	
	int TFSFy = tfsfBox[iMin];
	int TFSFY = tfsfBox[iMax];
	int TFSFx = tfsfBox[jMin];
	int TFSFX = tfsfBox[jMax];
	int TFSFz = tfsfBox[kMin];
	int TFSFZ = tfsfBox[kMax];
	
    int Sx = sx+1;
	int Sy = sy+1;

	double dt = 0.5;
	
	for (k=(TFSFz+1); k<=min(TFSFZ, sz-1); k++) for (j=TFSFx; j<=min(TFSFX, sx); j++) {
		i=TFSFy-1; ijk = _ijkHx;		
		if (signX>0) x0 = 1+(j-TFSFx); else x0 = 1+TFSFX-j;
		if (signY>0) y0 = 1;  else y0 = 1+TFSFY-TFSFy;
		val = dt*map0[10+(int)round((x0*nx+y0*ny) / dr)];
		Hxy[ijk] += val; Hx[ijk]  += val/2;
		
		if (bc[SYMi]==SYM_FULL) {
			i=TFSFY; ijk = _ijkHx;	
			if (signY>0) y0 = 1+TFSFY-TFSFy; else y0 = 1;
			val = dt*map0[10+(int)round((x0*nx+y0*ny) / dr)];
			Hxy[ijk] -= val; Hx[ijk] -= val/2;
		}
	}
	
	for (k=(TFSFz+1); k<=min(TFSFZ, sz-1); k++) for (i=TFSFy; i<=min(TFSFY, sy); i++) {
		j=TFSFx-1; ijk = _ijkHy;
		if (signY>0) y0 = 1+(i-TFSFy); else y0 = 1+TFSFY-i;
		if (signX>0) x0 = 1;  else x0 = TFSFX-TFSFx+1;
		val = dt*map0[10+(int)round((x0*nx+y0*ny) / dr)];
		Hyx[ijk] -= val; Hy[ijk] -= val/2;
	
		if (bc[SYMj]==SYM_FULL) {
			j=TFSFX; ijk = _ijkHy;
			if (signX>0) x0 = TFSFX-TFSFx+1; else x0 = 1;
			val = dt*map0[10+(int)round((x0*nx+y0*ny) / dr)];
			Hyx[ijk] += val; Hy[ijk] += val/2;
		}
	}
}

void cfdtd3d_pml(double *options, double *bcD, int sx, int sy, int sz,
				struct SplitField E,
				struct SplitField H,
				struct SplitField J,
				int N, int T0,
				double *SRCsignal, double *SRCparam, int SRCn, double *SRC, double *SRCfield,
				double *PRBparam, int PRBn, double *PRB,
				int PMLn, double *PML,
				double *PYBparam, int PYBn, double *PYBr, double *PYBi, double *PYBf, int PYBfn,
				double *DFTparam, int DFTn, double *DFTr, double *DFTi,
				unsigned char *MIndx, int MPARAMn, double *MPARAM,
				double *PLWparam, double *PLWprb, double *tfsfBox,
				double *play,
				int mapN, double *mapTE0, double *mapTE1, double *mapTE2, double *mapTEX, double *mapTEY,
				          double *mapTM0, double *mapTM1, double *mapTM2, double *mapTMX, double *mapTMY,
						  double *PYBtfsf) {
				
    //Up to 15% faster using direct pointers instead of struct 
	double *Ex  = E.x;
	double *Exy = E.xy;
	double *Exz = E.xz;
	double *Ey  = E.y;
	double *Eyx = E.yx;
	double *Eyz = E.yz;
	double *Ez  = E.z;
	double *Ezx = E.zx;
	double *Ezy = E.zy;
	 
	double *Hx  = H.x;
	double *Hxy = H.xy;
	double *Hxz = H.xz;
	double *Hy  = H.y;
	double *Hyx = H.yx;
	double *Hyz = H.yz;
	double *Hz  = H.z;
	double *Hzx = H.zx;
	double *Hzy = H.zy;
	
	double *Pxy = J.xy;
	double *Pxz = J.xz;
	double *Pyx = J.yx;
	double *Pyz = J.yz;
	double *Pzx = J.zx;
	double *Pzy = J.zy;

    double dt = options[0];
    int playFieldType = (int)options[1];
    int playDim = (int)options[2];
    int playIndx = (int)options[3];
    int playTFSF = (int)options[4];
		
    int Sx = sx+1;
    int Sy = sy+1;
    int Sz = sz+1;
    
    int ssx = sx-1;
    int ssy = sy-1;
    int ssz = sz-1;
    
    //For plasmonic materials size of Aux-fields
    int	SX = (int)options[7];
	int	SY = (int)options[6];
	int	SZ = (int)options[8];

	int	gTime = (int)options[9];	//Global time step
	
    int nCP = 6;
    
	
    int t, i, j, k, ijk;

    int bc[9];
    for (i=0; i<9; i++) {bc[i] = (int)bcD[i];}
	
	int TFSFy = tfsfBox[iMin];
	int TFSFY = tfsfBox[iMax];
	int TFSFx = tfsfBox[jMin];
	int TFSFX = tfsfBox[jMax];
	int TFSFz = tfsfBox[kMin];
	int TFSFZ = tfsfBox[kMax];
	
	//Source param field (EX, EY, EZ, HX, HY, HZ, TE, TM), pos (i, j, k), type (CW or pulse), amplitude, phase, f, df, MX, MY, nr
		
	//double amplitudeTE = 0.7071*10;
	//double amplitudeTM = 0*0.7071*10;
	
	double amplitudeTE = PLWparam[0];
	double amplitudeTM = PLWparam[1];
	double phaseTE = 0;
	double phaseTM = PLWparam[2];
	int MX = (int)PLWparam[3];
	int MY = (int)PLWparam[4];

	double nr = PLWparam[5];
	double mapSx = dt/nr;
	double mapSy = dt/nr;
	int mx = abs(MX);
	int my = abs(MY);

	int signX = 1; if (MX<0) signX = -1;
	int signY = 1; if (MY<0) signY = -1;

	double dr = 1/sqrt(mx*mx + my*my);
	double px = mx*dr;
	double py = my*dr;
	
	if (mx == 0) { px = 1;	py = 0; dr = 1; } //Enable horisontal wave

	int nMin = (max(mx, my));
	int nMax = (mapN-max(mx, my));

	double nx = mx*dr;
	double ny = my*dr;

	
	/*double lambda = 10; //In units of dx;	
	double src_f = 1/lambda;
	double src_df = 0.01;*/
	
	/*double lambdaMin = 2*30;
	double lambdaMax = 2*80;
	double src_f  = (1/lambdaMin + 1/lambdaMax)/2;
	double src_df = (1/lambdaMin - 1/lambdaMax);*/
	double src_f = SRCsignal[0];
	double src_df = SRCsignal[1];
	int src_type = (int)SRCsignal[2];
	
	double omega = 2*pi*src_f;
	double kVec0 = omega;
	
	//Numerical dispersion
	//Zero angle
	double _A = 0.5;
	double _B = 0.0;
	double _C = ((nr*nr)/(dt*dt))*sin(omega*dt/2)*sin(omega*dt/2);
	
	double kVec = 1;
	for (i=0; i<10; i++) {
		kVec -= (sin(_A*kVec)*sin(_A*kVec) + sin(_B*kVec)*sin(_B*kVec)-_C) / ( _A*sin(2*_A*kVec) + _B*sin(2*_B*kVec) );
	}
	double nrc0 = kVec/kVec0;   //Effective refractive index for zero inclination (MX = 0; MY = 0)
	
	//Actual angle
	 _A = px/2;
	 _B = py/2;
	 _C = ((nr*nr)/(dt*dt))*sin(omega*dt/2)*sin(omega*dt/2);
	
	 kVec = 1;
	for (i=0; i<10; i++) {
		kVec -= (sin(_A*kVec)*sin(_A*kVec) + sin(_B*kVec)*sin(_B*kVec)-_C) / ( _A*sin(2*_A*kVec) + _B*sin(2*_B*kVec) );
	}
	double nrc = kVec/kVec0;   //Unclear if adjustment of refractive index to compensates for numerical dispersion at wave injection is significant.
	
	
	
    //set up PML conductivities
    double ZXH[sx], ZYH[sy], ZZH[sz];
    double ZXE[sx], ZYE[sy], ZZE[sz];
    double zx, zy, zz;
    int id, jd, kd;
    int sym;
    
    if (bc[SYMj]==SYM_FULL) sym = 0; else sym = sx;
	for (j=0; j<sx; j++) {
		jd = min(j, sx-1-j + sym);  //sym removes PML at high index for symmetric geometry
		if (jd<PMLn) zx = 1-PML[jd*2]*dt; else zx = 1;
		ZXH[j] = zx;
		
		jd = min(j,sx-j + sym);		//sym removes PML at high index for symmetric geometry 
		if ((jd>0)&&(jd<PMLn)) zx = 0-PML[jd*2-1]*dt; else zx = 0;
		ZXE[j] = zx;
    }
    
    if (bc[SYMi]==SYM_FULL) sym = 0; else sym = sy;
    for (i=0; i<sy; i++) {
		id = min(i,sy-1-i + sym);	//sym removes PML at high index for symmetric geometry 
		if (id<PMLn) zy = 1-PML[id*2]*dt; else zy = 1;
		ZYH[i] = zy;
		
		id = min(i, sy-i + sym);	//sym removes PML at high index for symmetric geometry 
        if ((id>0)&&(id<PMLn)) zy = 0-PML[id*2-1]*dt; else zy = 0;
        ZYE[i] = zy;
	}

    if (bc[SYMk]==SYM_FULL) sym = 0; else sym = sz;	
	for (k=0; k<sz; k++) {
		kd = min(k,sz-1-k + sym);	//sym removes PML at high index for symmetric geometry 
		if (kd<PMLn) zz = 1-PML[kd*2]*dt; else zz = 1;
		ZZH[k] = zz;
		
		kd = min(k,sz-k + sym);		//sym removes PML at high index for symmetric geometry 
        if ((kd>0)&&(kd<PMLn)) zz = 0-PML[kd*2-1]*dt; else zz = 0;
        ZZE[k] = zz;
	}
	
	
	
	
	
	
	
	
	//Initialization (first run)
	
	if (T0==0) {  //Sets up TFSF data for PYB: 2: use raw data for entire cluster point, 1: check each component individually, 0: adjust data for entire cluster point
		int indx = 0;	
		for ( int PYBbox = 0 ; PYBbox < PYBn ; PYBbox++ ) {
			int posi =       (int)PYBparam[(PYBbox+1)*PYBN-11];
			int posj =       (int)PYBparam[(PYBbox+1)*PYBN-10];
			int posk =       (int)PYBparam[(PYBbox+1)*PYBN-9];
			int BSi 	=    (int)PYBparam[(PYBbox+1)*PYBN-8];
			int BSj 	=    (int)PYBparam[(PYBbox+1)*PYBN-7];
			int BSk 	=    (int)PYBparam[(PYBbox+1)*PYBN-6];
			int tfsf =       (int)PYBparam[(PYBbox+1)*PYBN-3];
			int side0 = (int)PYBparam[(PYBbox+1)*PYBN+0];
			int side1 = (int)PYBparam[(PYBbox+1)*PYBN+1];
			int side2 = (int)PYBparam[(PYBbox+1)*PYBN+2];
			int side3 = (int)PYBparam[(PYBbox+1)*PYBN+3];
			int side4 = (int)PYBparam[(PYBbox+1)*PYBN+4];
			int side5 = (int)PYBparam[(PYBbox+1)*PYBN+5];
			
			double TFSF = 0;		//SF if -1, TF if 1, DF if 0
			if (tfsf==SF) TFSF = -1; else if (tfsf==TF) TFSF = 1;
				
			int TFSFy = tfsfBox[iMin];
			int TFSFY = tfsfBox[iMax];
			int TFSFx = tfsfBox[jMin];
			int TFSFX = tfsfBox[jMax];
			int TFSFz = tfsfBox[kMin];
			int TFSFZ = tfsfBox[kMax];
			
			//withinTF+TFSF is zero when adjustment is required.
			
    		for (int b=-BSk; b<=BSk; b++) {
	    		for (int a=-BSj; a<=BSj; a++) {
	    			//constant y: (EzHx-ExHz)
		    		i = posi-BSi;
		    		j = posj+a;
		    		k = posk+b;
					
		    		if ((side0==1)&&(i>=0)&&(j>=0)&&(k>=0)&&(i<sy)&&(j<sx)&&(k<sz)) {
						double withinTF = 0;
						if ((i>=TFSFy)&&(i<(TFSFY-1))&&(j>=TFSFx)&&(j<(TFSFX-1))&&(k>TFSFz)&&(k<=(TFSFZ-1))) withinTF =  1;
						else if ((i<(TFSFy-1))||(i>TFSFY)||(j<(TFSFx-1))||(j>TFSFX)||(k<TFSFz)||(k>TFSFZ))   withinTF = -1;
						PYBtfsf[indx++] = (!TFSF)?2:fabs(withinTF+TFSF);
		    		} 
					
		    		i = posi+BSi;
		    		if ((side1==1)&&(i>=0)&&(j>=0)&&(k>=0)&&(i<sy)&&(j<sx)&&(k<sz)) {
						double withinTF = 0;
						if ((i>=TFSFy)&&(i<(TFSFY-1))&&(j>=TFSFx)&&(j<(TFSFX-1))&&(k>TFSFz)&&(k<=(TFSFZ-1))) withinTF =  1;
						else if ((i<(TFSFy-1))||(i>TFSFY)||(j<(TFSFx-1))||(j>TFSFX)||(k<TFSFz)||(k>TFSFZ))   withinTF = -1;
						PYBtfsf[indx++] = (!TFSF)?2:fabs(withinTF+TFSF);
		    		}
		    	}
		    }
		    for (int b=-BSk; b<=BSk; b++) {
	    		for (int a=-BSi; a<=BSi; a++) {
		    		//constant x: (EyHz-EzHy)
		    		i = posi+a;
		    		j = posj-BSj;
		    		k = posk+b;
										
		    		if ((side2==1)&&(i>=0)&&(j>=0)&&(k>=0)&&(i<sy)&&(j<sx)&&(k<sz)) {
						double withinTF = 0;
						if ((i>=TFSFy)&&(i<(TFSFY-1))&&(j>=TFSFx)&&(j<(TFSFX-1))&&(k>TFSFz)&&(k<=(TFSFZ-1))) withinTF =  1;
						else if ((i<(TFSFy-1))||(i>TFSFY)||(j<(TFSFx-1))||(j>TFSFX)||(k<TFSFz)||(k>TFSFZ))   withinTF = -1;
						PYBtfsf[indx++] = (!TFSF)?2:fabs(withinTF+TFSF);
					} 
					
		    		j = posj+BSj;
		    		if ((side3==1)&&(i>=0)&&(j>=0)&&(k>=0)&&(i<sy)&&(j<sx)&&(k<sz)) {
						double withinTF = 0;
						if ((i>=TFSFy)&&(i<(TFSFY-1))&&(j>=TFSFx)&&(j<(TFSFX-1))&&(k>TFSFz)&&(k<=(TFSFZ-1))) withinTF =  1;
						else if ((i<(TFSFy-1))||(i>TFSFY)||(j<(TFSFx-1))||(j>TFSFX)||(k<TFSFz)||(k>TFSFZ))   withinTF = -1;
						PYBtfsf[indx++] = (!TFSF)?2:fabs(withinTF+TFSF);
		    		} 
		    	}
		    }
		    
		    for (int b=-BSj; b<=BSj; b++) {
	    		for (int a=-BSi; a<=BSi; a++) {
		    		//constant z: (ExHy-EyHx)
		    		i = posi+a;
		    		j = posj+b;
		    		k = posk-BSk;		
					
		    		if ((side4==1)&&(i>=0)&&(j>=0)&&(k>=0)&&(i<sy)&&(j<sx)&&(k<sz)) {
						double withinTF = 0;
						if ((i>=TFSFy)&&(i<(TFSFY-1))&&(j>=TFSFx)&&(j<(TFSFX-1))&&(k>TFSFz)&&(k<=(TFSFZ-1))) withinTF =  1;
						else if ((i<(TFSFy-1))||(i>TFSFY)||(j<(TFSFx-1))||(j>TFSFX)||(k<TFSFz)||(k>TFSFZ))   withinTF = -1;
						PYBtfsf[indx++] = (!TFSF)?2:fabs(withinTF+TFSF);
		    		}
					
		    		k = posk+BSk;
		    		if ((side5==1)&&(i>=0)&&(j>=0)&&(k>=0)&&(i<sy)&&(j<sx)&&(k<sz)) {
						double withinTF = 0;
						if ((i>=TFSFy)&&(i<(TFSFY-1))&&(j>=TFSFx)&&(j<(TFSFX-1))&&(k>TFSFz)&&(k<=(TFSFZ-1))) withinTF =  1;
						else if ((i<(TFSFy-1))||(i>TFSFY)||(j<(TFSFx-1))||(j>TFSFX)||(k<TFSFz)||(k>TFSFZ))   withinTF = -1;
						PYBtfsf[indx++] = (!TFSF)?2:fabs(withinTF+TFSF);
		    		}
	    		}
	    	}	
			
		}
			
	}
	
	
	
	
	int m = -1; //Material at current point (default no material)
	double *A, *B, *C;
	int CPn;
	
	//main loop start
    for ( t = 0 ; t < N ; t++ ) {
		
		//update E-field
		
		if (amplitudeTM) {    //1DMAP TM
			for ( i = 0 ; i < mapN ; i++ ) {				//Update 1DMAP
				mapTM2[i] = mapTM1[i];
				mapTM1[i] = mapTM0[i];	
			}
			for ( i = nMin ; i < nMax ; i++ ) {
				mapTM0[i] = -mapTM2[i] + 2*(1-mapSx*mapSx-mapSy*mapSy)*mapTM1[i]  + (mapSx*mapSx)*(mapTM1[i+mx] + mapTM1[i-mx]) + (mapSy*mapSy)*(mapTM1[i+my] + mapTM1[i-my]);
			}
			
			//1D-MAP source
			for ( i = 0; i <nMin; i++) {
				mapTM0[(nMin-1)-i] = sourceSignal(src_f, src_df, amplitudeTM*(nrc0/nrc), phaseTM, src_type, -i*dr-0.5*(px+py), (T0+t-0.5)*dt/nr, nr);
			}
			
			for ( i = nMax; i <mapN; i++) {  //MUR Absorbing boundary
				mapTM0[i] = mapTM1[i-1] + (dt/nrc-dr)/(dt/nrc+dr)*(mapTM0[i-1] - mapTM1[i]);
			}
		}
		
		if (amplitudeTE) {	 //1DMAP TE
			for ( i = nMin ; i < mapN ; i++ ) {			//Compute E-field
				    mapTEX[i] = mapTEX[i] + (dt/(nr*nr))*(mapTE0[i-my] - mapTE0[i]);
					mapTEY[i] = mapTEY[i] + (dt/(nr*nr))*(mapTE0[i] - mapTE0[i-mx]);
			}
		}
		
        int apx = 0;	//Reset auxiliary field indices
        int apy = 0;
        int apz = 0;
        
		for (k=0; k<sz; k++) {			//3D E
	    	zz = ZZE[k];	
	        for (j=0; j<sx; j++) {
	            zx = ZXE[j];
	            for (i=0; i<sy; i++) {
	                zy = ZYE[i];
	                					
	                int m0 = MIndx[(i)+(j)*sy+(k)*sy*sx];
	                
	                if ((i>0)&&(k>0)) {
		                int mp = max(MIndx[(i-1)+(j)*sy+(k)*sy*sx], MIndx[(i)+(j)*sy+(k-1)*sy*sx]);
		                mp =  max(mp, MIndx[(i-1)+(j)*sy+(k-1)*sy*sx]);
		                if (m0>mp) mp = m0;	//mp = max value of neighbouring materials
	                    if (mp!=m) {
	                    	m = mp;
	                    	A = &MPARAM[m*MPARAMn];
							if (m>=254) A = &MPARAM[0]; //reserved material index 254:PMC, 255:PEC, set vacuum properties
	                    	CPn = (int)A[6];
	                    	B = &A[7];
	                    	C = &A[7+2*CPn];
	                    }
							
	                	ijk = _ijkEx;
	                	double  oExy = Exy[ijk];
	                	double  oExz = Exz[ijk];
	                	double 	ooExy;
	                	double 	ooExz;
	                	
	                	//Standard-update
		                Exy[ijk] = oExy*(zy + A[0]) + A[1]*(Hz[_ijkHz]-Hz[_dHz_dy]);
		                Exz[ijk] = oExz*(zz + A[0]) - A[1]*(Hy[_ijkHy]-Hy[_dHy_dz]);
						
						if (CPn>=0) {		//Drude
							Exy[ijk] -= A[2]*Pxy[apx];
							Exz[ijk] -= A[2]*Pxz[apx];
						
							if (CPn>0) {	//Critical Points
								ooExy = Pxy[apx+SX];
			                	ooExz = Pxz[apx+SX];
				                Exy[ijk] -= A[5]*ooExy;
				                Exz[ijk] -= A[5]*ooExz;
				                
				                for (int n = 0; n<CPn; n++) {
				                	int m = 2+2*n;
				                	Exy[ijk] -= ( B[0+2*n]*Pxy[apx+(m+0)*SX] + B[1+2*n]*Pxy[apx+(m+1)*SX] );
				                	Exz[ijk] -= ( B[0+2*n]*Pxz[apx+(m+0)*SX] + B[1+2*n]*Pxz[apx+(m+1)*SX] );
				                }
				 		    }
		                }

						if (m==255) { Exy[ijk]=0; Exz[ijk] = 0; }			//if PEC material						
		                Ex[ijk]  = Exy[ijk] + Exz[ijk];				
		                
		                if (CPn>=0) {		//Drude
		                	Pxy[apx] = A[3]*Pxy[apx] - A[4]*(oExy+Exy[ijk]);
		                	Pxz[apx] = A[3]*Pxz[apx] - A[4]*(oExz+Exz[ijk]);
          
			                if (CPn>0) {	//Critical Points
								double oPxy, oPxz, ooPxy, ooPxz;
								
				                for (int n = 0; n<CPn; n++) { //Compute polarization for critical point n
				                	int m = 2+2*n;
					                oPxy  = Pxy[apx+SX*(m+0)];
					                oPxz  = Pxz[apx+SX*(m+0)];
					                ooPxy = Pxy[apx+SX*(m+1)];
					                ooPxz = Pxz[apx+SX*(m+1)];
				
					                Pxy[apx+SX*(m+0)] = C[0+5*n]*oPxy + C[1+5*n]*ooPxy + C[2+5*n]*Exy[ijk] + C[3+5*n]*oExy + C[4+5*n]*ooExy;
					                Pxz[apx+SX*(m+0)] = C[0+5*n]*oPxz + C[1+5*n]*ooPxz + C[2+5*n]*Exz[ijk] + C[3+5*n]*oExz + C[4+5*n]*ooExz;
					                Pxy[apx+SX*(m+1)] = oPxy; //Store old polarization for critical point 1
					                Pxz[apx+SX*(m+1)] = oPxz;
				               }
				                
			                	Pxy[apx+SX*1] = oExy; //Store old E-field
			                	Pxz[apx+SX*1] = oExz;
				            }
				            apx++;  //Auxiliary field index (same index for all).
			            }   
	                }
	                
	                if ((j>0)&&(k>0)) {
	                	int mp = max(MIndx[(i)+(j-1)*sy+(k)*sy*sx], MIndx[(i)+(j)*sy+(k-1)*sy*sx]);
	                	mp =  max(mp, MIndx[(i)+(j-1)*sy+(k-1)*sy*sx]);
		                if (m0>mp) mp = m0;	//mp = max value of neighbouring materials
	                    if (mp!=m) {
	                    	m = mp;
	                    	A = &MPARAM[m*MPARAMn];
							if (m>=254) A = &MPARAM[0]; //reserved material index 254:PMC, 255:PEC, set vacuum properties
	                    	CPn = (int)A[6];
	                    	B = &A[7];
	                    	C = &A[7+2*CPn];
	                    }
	                    
	                	ijk = _ijkEy;
	                	double oEyz = Eyz[ijk];
	                	double oEyx = Eyx[ijk];
	                	double ooEyz;
	                	double ooEyx;
	                	
	                	//Standard-update
		                Eyz[ijk] = oEyz*(zz + A[0]) + A[1]*(Hx[_ijkHx]-Hx[_dHx_dz]);
		                Eyx[ijk] = oEyx*(zx + A[0]) - A[1]*(Hz[_ijkHz]-Hz[_dHz_dx]);
						
						if (CPn>=0) {		//Drude
							Eyz[ijk] -= A[2]*Pyz[apy];
							Eyx[ijk] -= A[2]*Pyx[apy];
												
							if (CPn>0) {	//Critical Points
								ooEyz = Pyz[apy+SY];
			                	ooEyx = Pyx[apy+SY];
				                Eyz[ijk] -= A[5]*ooEyz;
				                Eyx[ijk] -= A[5]*ooEyx;
				                
				                for (int n = 0; n<CPn; n++) {
				                	int m = 2+2*n;
				                	Eyz[ijk] -= ( B[0+2*n]*Pyz[apy+(m+0)*SY] + B[1+2*n]*Pyz[apy+(m+1)*SY] );
				                	Eyx[ijk] -= ( B[0+2*n]*Pyx[apy+(m+0)*SY] + B[1+2*n]*Pyx[apy+(m+1)*SY] );
				                }
				            }
						}
						
						if (m==255) {Eyz[ijk]=0; Eyx[ijk]=0;}			//if PEC material
		                Ey[ijk]  = Eyz[ijk] + Eyx[ijk];
		                
						
		                if (CPn>=0) {		//Drude
		                	Pyz[apy] = A[3]*Pyz[apy] - A[4]*(oEyz+Eyz[ijk]);
		                	Pyx[apy] = A[3]*Pyx[apy] - A[4]*(oEyx+Eyx[ijk]);
		                
			                if (CPn>0) {	//Critical Points
								double oPyz, oPyx, ooPyz, ooPyx;
								
				                for (int n = 0; n<CPn; n++) { //Compute polarization for critical point n
				                	int m = 2+2*n;
					                oPyz  = Pyz[apy+SY*(m+0)];
					                oPyx  = Pyx[apy+SY*(m+0)];
					                ooPyz = Pyz[apy+SY*(m+1)];
					                ooPyx = Pyx[apy+SY*(m+1)];
				
					                Pyz[apy+SY*(m+0)] = C[0+5*n]*oPyz + C[1+5*n]*ooPyz + C[2+5*n]*Eyz[ijk] + C[3+5*n]*oEyz + C[4+5*n]*ooEyz;
					                Pyx[apy+SY*(m+0)] = C[0+5*n]*oPyx + C[1+5*n]*ooPyx + C[2+5*n]*Eyx[ijk] + C[3+5*n]*oEyx + C[4+5*n]*ooEyx;
					                Pyz[apy+SY*(m+1)] = oPyz; //Store old polarization for critical point 1
					                Pyx[apy+SY*(m+1)] = oPyx;
					            }
	
				                Pyz[apy+SY] = oEyz; //Store old E-field
				                Pyx[apy+SY] = oEyx;
			            	}
			            	apy++;  //Auxiliary field index (same index for all).
			            } 
		            }
		            
		            if ((i>0)&&(j>0)) {
		            	int mp = max(MIndx[(i-1)+(j)*sy+(k)*sy*sx], MIndx[(i)+(j-1)*sy+(k)*sy*sx]);
		            	mp =  max(mp, MIndx[(i-1)+(j-1)*sy+(k)*sy*sx]);
		                if (m0>mp) mp = m0;	//mp = max value of neighbouring materials
	                    if (mp!=m) {
	                    	m = mp;
	                    	A = &MPARAM[m*MPARAMn];
							if (m>=254) A = &MPARAM[0]; //reserved material index 254:PMC, 255:PEC, set vacuum properties
	                    	CPn = (int)A[6];
	                    	B = &A[7];
	                    	C = &A[7+2*CPn];	
	                    }
	                    
		            	ijk = _ijkEz;
		            	double oEzx = Ezx[ijk];
	                	double oEzy = Ezy[ijk];
	                	double ooEzx;
	                	double ooEzy;
	                	
	                	//Standard-update
						Ezx[ijk] = oEzx*(zx + A[0]) + A[1]*(Hy[_ijkHy]-Hy[_dHy_dx]);
		                Ezy[ijk] = oEzy*(zy + A[0]) - A[1]*(Hx[_ijkHx]-Hx[_dHx_dy]);
		                
		                if (CPn>=0) {		//Drude
		                	Ezx[ijk] -= A[2]*Pzx[apz];
		                	Ezy[ijk] -= A[2]*Pzy[apz];

			                if (CPn>0) {	//Critical Points	
								ooEzx = Pzx[apz+SZ];
			                	ooEzy = Pzy[apz+SZ];
				                Ezx[ijk] -= A[5]*ooEzx;    //Sum of all critical points C[4+0*5] + C[4+1*5] + C[4+2*5]
				                Ezy[ijk] -= A[5]*ooEzy;
				                
				                for (int n = 0; n<CPn; n++) {
				                	int m = 2+2*n;
				                	Ezx[ijk] -= ( B[0+2*n]*Pzx[apz+(m+0)*SZ] + B[1+2*n]*Pzx[apz+(m+1)*SZ] );
				                	Ezy[ijk] -= ( B[0+2*n]*Pzy[apz+(m+0)*SZ] + B[1+2*n]*Pzy[apz+(m+1)*SZ] );
				                }
				            }
			            }
			            
						if (m==255) {Ezx[ijk]=0; Ezy[ijk]=0; }			//if PEC material
		                Ez[ijk]  = Ezx[ijk] + Ezy[ijk];
						
												
						
						
		                if (CPn>=0) {		//Drude
		                	Pzx[apz] = A[3]*Pzx[apz] - A[4]*(oEzx+Ezx[ijk]);
		                	Pzy[apz] = A[3]*Pzy[apz] - A[4]*(oEzy+Ezy[ijk]);
		                
			                if (CPn>0) {	//Critical Points
								double oPzx, oPzy, ooPzx, ooPzy;
								
				                for (int n = 0; n<CPn; n++) { //Compute polarization for critical points
				                	int m = 2+2*n;
					                oPzx  = Pzx[apz+SZ*(m+0)];
					                oPzy  = Pzy[apz+SZ*(m+0)];
					                ooPzx = Pzx[apz+SZ*(m+1)];
					                ooPzy = Pzy[apz+SZ*(m+1)];
				
					                Pzx[apz+SZ*(m+0)] = C[0+5*n]*oPzx + C[1+5*n]*ooPzx + C[2+5*n]*Ezx[ijk] + C[3+5*n]*oEzx + C[4+5*n]*ooEzx;
					                Pzy[apz+SZ*(m+0)] = C[0+5*n]*oPzy + C[1+5*n]*ooPzy + C[2+5*n]*Ezy[ijk] + C[3+5*n]*oEzy + C[4+5*n]*ooEzy;
					                Pzx[apz+SZ*(m+1)] = oPzx; //Store old polarization for critical point 1
					                Pzy[apz+SZ*(m+1)] = oPzy;
				                }
				                
			                	Pzx[apz+SZ] = oEzx; //Store old E-field
			                	Pzy[apz+SZ] = oEzy;
			                	
				            } 
				            apz++;  //Auxiliary field index (same index for all).
				        }
		            }
		            
	            }
	        }
	    }					
		
		
		if (amplitudeTM) updateTFSF_ETM(Ex,Exz,Ey,Eyz,Ez,Ezx,Ezy,bc,tfsfBox,mapTMX,mapTMY,nx,ny,dr,signX,signY,sx,sy,sz,dt,nr);
		if (amplitudeTE) updateTFSF_ETE(Ex, Ey, Exy, Eyx, bc, tfsfBox, mapTE0, nx, ny, dr, signX, signY, sx, sy, sz, dt, nr);
		
		
		//point dipole sources
		for ( int source = 0 ; source < SRCn ; source++ ) { 
			i = 			(int)SRCparam[(source+1)*SRCN-11];
			j = 			(int)SRCparam[(source+1)*SRCN-10];
			k = 			(int)SRCparam[(source+1)*SRCN-9];
			int srcType = 	(int)SRCparam[(source+1)*SRCN-8];
			int pntType = 	(int)SRCparam[(source+1)*SRCN-7];
			double phase =  SRCparam[(source+1)*SRCN-3];
			
			int m0 = MIndx[(i)+(j)*sy+(k)*sy*sx];
			double amplitude = 10.0;
			double nr = 1.0;  //Not very relevant here, just modifies phase of signal not spectrum.
			double signal = 0.0;
			switch (srcType) {
				case EX:
					if ((i>0)&&(k>0)) {
		                int mp = max(MIndx[(i-1)+(j)*sy+(k)*sy*sx], MIndx[(i)+(j)*sy+(k-1)*sy*sx]);
		                mp =  max(mp, MIndx[(i-1)+(j)*sy+(k-1)*sy*sx]);
						mp =  max(mp, m0); //mp = max value of neighbouring materials
		                A = &MPARAM[mp*MPARAMn]; if (mp>=254) A = &MPARAM[0];
						m=-1; //Update material parameters for next E-field update
						signal = sourceSignal(src_f, src_df, amplitude, phase, src_type, 0.0, (T0+t)*dt/nr, nr);
						SRC[t+source*N] = signal;
						if (pntType == CLUSTERPOINT) {
							addValueToFieldAtClusterPoint(i, j, k, Exy, EX, sx, sy, sz, A[1]*signal/2);
							addValueToFieldAtClusterPoint(i, j, k, Exz, EX, sx, sy, sz, A[1]*signal/2);
							addValueToFieldAtClusterPoint(i, j, k, Ex,  EX, sx, sy, sz, A[1]*signal);
						} else if (pntType == SINGLEPOINT) {
							addValueToFieldAtSinglePoint(i, j, k, Exy, EX, sx, sy, sz, A[1]*signal/2);
							addValueToFieldAtSinglePoint(i, j, k, Exz, EX, sx, sy, sz, A[1]*signal/2);
							addValueToFieldAtSinglePoint(i, j, k, Ex,  EX, sx, sy, sz, A[1]*signal);
						}
					}
					break;
				case EY:
					if ((j>0)&&(k>0)) {
	                	int mp = max(MIndx[(i)+(j-1)*sy+(k)*sy*sx], MIndx[(i)+(j)*sy+(k-1)*sy*sx]);
	                	mp =  max(mp, MIndx[(i)+(j-1)*sy+(k-1)*sy*sx]);
						mp =  max(mp, m0); //mp = max value of neighbouring materials
						A = &MPARAM[mp*MPARAMn]; if (mp>=254) A = &MPARAM[0];
						m=-1; //Update material parameters for next E-field update
						signal = sourceSignal(src_f, src_df, amplitude, phase, src_type, 0.0, (T0+t)*dt/nr, nr);
						SRC[t+source*N] = signal;
						if (pntType == CLUSTERPOINT) {
							addValueToFieldAtClusterPoint(i, j, k, Eyx, EY, sx, sy, sz, A[1]*signal/2);
							addValueToFieldAtClusterPoint(i, j, k, Eyz, EY, sx, sy, sz, A[1]*signal/2);
							addValueToFieldAtClusterPoint(i, j, k, Ey,  EY, sx, sy, sz, A[1]*signal);
						} else if (pntType == SINGLEPOINT) {
							addValueToFieldAtSinglePoint(i, j, k, Eyx, EY, sx, sy, sz, A[1]*signal/2);
							addValueToFieldAtSinglePoint(i, j, k, Eyz, EY, sx, sy, sz, A[1]*signal/2);
							addValueToFieldAtSinglePoint(i, j, k, Ey,  EY, sx, sy, sz, A[1]*signal);
						}
					}
					break;
				case EZ:
					 if ((i>0)&&(j>0)) {
		            	int mp = max(MIndx[(i-1)+(j)*sy+(k)*sy*sx], MIndx[(i)+(j-1)*sy+(k)*sy*sx]);
		            	mp =  max(mp, MIndx[(i-1)+(j-1)*sy+(k)*sy*sx]);
						mp =  max(mp, m0); //mp = max value of neighbouring materials
						A = &MPARAM[mp*MPARAMn]; if (mp>=254) A = &MPARAM[0];
						m=-1; //Update material parameters for next E-field update
						signal = sourceSignal(src_f, src_df, amplitude, phase, src_type, 0.0, (T0+t)*dt/nr, nr);
						SRC[t+source*N] = signal;
						if (pntType == CLUSTERPOINT) {
							addValueToFieldAtClusterPoint(i, j, k, Ezx, EZ, sx, sy, sz, A[1]*signal/2);
							addValueToFieldAtClusterPoint(i, j, k, Ezy, EZ, sx, sy, sz, A[1]*signal/2);
							addValueToFieldAtClusterPoint(i, j, k, Ez,  EZ, sx, sy, sz, A[1]*signal);
						} else if (pntType == SINGLEPOINT) {
							addValueToFieldAtSinglePoint(i, j, k, Ezx, EZ, sx, sy, sz, A[1]*signal/2);
							addValueToFieldAtSinglePoint(i, j, k, Ezy, EZ, sx, sy, sz, A[1]*signal/2);
							addValueToFieldAtSinglePoint(i, j, k, Ez,  EZ, sx, sy, sz, A[1]*signal);
						}
					}
					break;
				case HX:
					signal = sourceSignal(src_f, src_df, amplitude, phase, src_type, 0.0, (T0+t)*dt/nr, nr);
					SRC[t+source*N] = signal;
					if (pntType == CLUSTERPOINT) {
						addValueToFieldAtClusterPoint(i, j, k, Hxy, HX, sx, sy, sz, signal*dt/2);
						addValueToFieldAtClusterPoint(i, j, k, Hxz, HX, sx, sy, sz, signal*dt/2);
					} else if (pntType == SINGLEPOINT) {
						addValueToFieldAtSinglePoint(i, j, k, Hxy, HX, sx, sy, sz, signal*dt/2);
						addValueToFieldAtSinglePoint(i, j, k, Hxz, HX, sx, sy, sz, signal*dt/2);
					}
					break;
				case HY:
					signal = sourceSignal(src_f, src_df, amplitude, phase, src_type, 0.0, (T0+t)*dt/nr, nr);
					SRC[t+source*N] = signal;
					if (pntType == CLUSTERPOINT) {
						addValueToFieldAtClusterPoint(i, j, k, Hyx, HY, sx, sy, sz, signal*dt/2);
						addValueToFieldAtClusterPoint(i, j, k, Hyz, HY, sx, sy, sz, signal*dt/2);
					} else if (pntType == SINGLEPOINT) {
						addValueToFieldAtSinglePoint(i, j, k, Hyx, HY, sx, sy, sz, signal*dt/2);
						addValueToFieldAtSinglePoint(i, j, k, Hyz, HY, sx, sy, sz, signal*dt/2);
					}
					break;
				case HZ:
					signal = sourceSignal(src_f, src_df, amplitude, phase, src_type, 0.0, (T0+t)*dt/nr, nr);
					SRC[t+source*N] = signal;
					if (pntType == CLUSTERPOINT) {
						addValueToFieldAtClusterPoint(i, j, k, Hzx, HZ, sx, sy, sz, signal*dt/2);
						addValueToFieldAtClusterPoint(i, j, k, Hzy, HZ, sx, sy, sz, signal*dt/2);
					} else if (pntType == SINGLEPOINT) {
						addValueToFieldAtSinglePoint(i, j, k, Hzx, HZ, sx, sy, sz, signal*dt/2);
						addValueToFieldAtSinglePoint(i, j, k, Hzy, HZ, sx, sy, sz, signal*dt/2);
					}
					break;
				case PLW:
					//H1D[i] += SRC[t+source*N]*dt;
					break;
				case EDP:
					if ((i>0)&&(j>0)&&(k>0)) {
						int mp = max(MIndx[(i-1)+(j)*sy+(k)*sy*sx], MIndx[(i)+(j)*sy+(k-1)*sy*sx]);
		                mp =  max(mp, MIndx[(i-1)+(j)*sy+(k-1)*sy*sx]);
						mp =  max(mp, MIndx[(i)+(j-1)*sy+(k)*sy*sx]);
	                	mp =  max(mp, MIndx[(i)+(j-1)*sy+(k-1)*sy*sx]);
		            	mp =  max(mp, MIndx[(i-1)+(j-1)*sy+(k)*sy*sx]);
						mp =  max(mp, m0);
						
		                A = &MPARAM[mp*MPARAMn]; if (mp>=254) A = &MPARAM[0];
						m=-1; //Update material parameters for next E-field update
						signal = sourceSignal(src_f, src_df, amplitude, 0, src_type, 0.0, (T0+t)*dt/nr, nr);
						SRC[t+source*N] = signal;
						
						double nx = SRCparam[(source+1)*SRCN-6];
						double ny = SRCparam[(source+1)*SRCN-5];
						double nz = SRCparam[(source+1)*SRCN-4];
						
						double phasex = SRCparam[(source+1)*SRCN-3];
						double phasey = SRCparam[(source+1)*SRCN-2];
						double phasez = SRCparam[(source+1)*SRCN-1];
						
						signal = sourceSignal(src_f, src_df, amplitude, phasex, src_type, 0.0, (T0+t)*dt/nr, nr);
						addValueToFieldAtClusterPoint(i, j, k, Exy, EX, sx, sy, sz, nx*A[1]*signal/2);
						addValueToFieldAtClusterPoint(i, j, k, Exz, EX, sx, sy, sz, nx*A[1]*signal/2);
						addValueToFieldAtClusterPoint(i, j, k, Ex,  EX, sx, sy, sz, nx*A[1]*signal);
						
						signal = sourceSignal(src_f, src_df, amplitude, phasey, src_type, 0.0, (T0+t)*dt/nr, nr);
						addValueToFieldAtClusterPoint(i, j, k, Eyx, EY, sx, sy, sz, ny*A[1]*signal/2);
						addValueToFieldAtClusterPoint(i, j, k, Eyz, EY, sx, sy, sz, ny*A[1]*signal/2);
						addValueToFieldAtClusterPoint(i, j, k, Ey,  EY, sx, sy, sz, ny*A[1]*signal);
						
						signal = sourceSignal(src_f, src_df, amplitude, phasez, src_type, 0.0, (T0+t)*dt/nr, nr);
						addValueToFieldAtClusterPoint(i, j, k, Ezx, EZ, sx, sy, sz, nz*A[1]*signal/2);
						addValueToFieldAtClusterPoint(i, j, k, Ezy, EZ, sx, sy, sz, nz*A[1]*signal/2);
						addValueToFieldAtClusterPoint(i, j, k, Ez,  EZ, sx, sy, sz, nz*A[1]*signal);
					}
					break;
				case MDP: {
						signal = sourceSignal(src_f, src_df, amplitude, 0, src_type, 0.0, (T0+t)*dt/nr, nr);
						SRC[t+source*N] = signal;
						
						double nx = SRCparam[(source+1)*SRCN-6];
						double ny = SRCparam[(source+1)*SRCN-5];
						double nz = SRCparam[(source+1)*SRCN-4];
						
						double phasex = SRCparam[(source+1)*SRCN-3];
						double phasey = SRCparam[(source+1)*SRCN-2];
						double phasez = SRCparam[(source+1)*SRCN-1];
						
						signal = sourceSignal(src_f, src_df, amplitude, phasex, src_type, 0.0, (T0+t)*dt/nr, nr);
						addValueToFieldAtClusterPoint(i, j, k, Hxy, HX, sx, sy, sz, nx*dt*signal/2);
						addValueToFieldAtClusterPoint(i, j, k, Hxz, HX, sx, sy, sz, nx*dt*signal/2);
						
						signal = sourceSignal(src_f, src_df, amplitude, phasey, src_type, 0.0, (T0+t)*dt/nr, nr);
						addValueToFieldAtClusterPoint(i, j, k, Hyx, HY, sx, sy, sz, ny*dt*signal/2);
						addValueToFieldAtClusterPoint(i, j, k, Hyz, HY, sx, sy, sz, ny*dt*signal/2);
						
						signal = sourceSignal(src_f, src_df, amplitude, phasez, src_type, 0.0, (T0+t)*dt/nr, nr);
						addValueToFieldAtClusterPoint(i, j, k, Hzx, HZ, sx, sy, sz, nz*dt*signal/2);
						addValueToFieldAtClusterPoint(i, j, k, Hzy, HZ, sx, sy, sz, nz*dt*signal/2);
					}
					break;
					
			}
	    }
		
		//Symmetry BC at sx
		if (bc[SYMj]!=SYM_FULL) {
			int n;
			int sign = bc[SYMj];
			int symx = bc[SYMx];
			for (k=0; k<=sz; k++) {
		        j = sx;
	            for (i=0; i<=sy; i++) {
	            	ijk = _ijkEy;
	            	if (i<sy) { Eyz[ijk] = sign*Eyz[(i)+(j-symx)*sy+(k)*sy*Sx];
	            				Eyx[ijk] = sign*Eyx[(i)+(j-symx)*sy+(k)*sy*Sx];
	            				Ey[ijk]  = Eyz[ijk] + Eyx[ijk];
	            	} 
	            	
	            	ijk = _ijkEz;
	            	if (k<sz) { Ezx[ijk] = sign*Ezx[(i)+(j-symx)*Sy+(k)*Sy*Sx];
	            	            Ezy[ijk] = sign*Ezy[(i)+(j-symx)*Sy+(k)*Sy*Sx];
	            	            Ez[ijk]  = Ezx[ijk] + Ezy[ijk];
	            	}
	            }
			}
		}
		
		//Symmetry BC sy
		if (bc[SYMi]!=SYM_FULL) {
			int n; int nCP = 6;
			int sign = bc[SYMi];
			int symy = bc[SYMy];
			for (k=0; k<=sz; k++) {
	            for (j=0; j<=sx; j++) {	
	            	i = sy;
	            	ijk = _ijkEx;
	            	if (j<sx) { Exy[ijk] = sign*Exy[(i-symy)+(j)*Sy+(k)*Sy*sx];
	            				Exz[ijk] = sign*Exz[(i-symy)+(j)*Sy+(k)*Sy*sx];
	            				Ex[ijk]  = Exy[ijk] + Exz[ijk];
	            	}
	            	
	            	ijk = _ijkEz;
	            	if (k<sz) { Ezx[ijk] = sign*Ezx[(i-symy)+(j)*Sy+(k)*Sy*Sx];
	            	            Ezy[ijk] = sign*Ezy[(i-symy)+(j)*Sy+(k)*Sy*Sx];
	            	            Ez[ijk]  = Ezx[ijk] + Ezy[ijk];
	            	}
	            }
			}
		}
		
		//Symmetry BC at sz
		if (bc[SYMk]!=SYM_FULL) {	
			int n; int nCP = 6;
			int sign = bc[SYMk];
			int symz = bc[SYMz];
		   	k = sz;
		   	for (j=0; j<=sx; j++) {
	            for (i=0; i<=sy; i++) {
	            	ijk = _ijkEx;
	            	if (j<sx) { Exy[ijk] = sign*Exy[(i)+(j)*Sy+(k-symz)*Sy*sx];
	            				Exz[ijk] = sign*Exz[(i)+(j)*Sy+(k-symz)*Sy*sx];
	            				Ex[ijk]  = Exy[ijk] + Exz[ijk];
	            	}
	            	
	            	ijk = _ijkEy;
	            	if (i<sy) { Eyz[ijk] = sign*Eyz[(i)+(j)*sy+(k-symz)*sy*Sx];
	            				Eyx[ijk] = sign*Eyx[(i)+(j)*sy+(k-symz)*sy*Sx];
	            				Ey[ijk]  = Eyz[ijk] + Eyx[ijk];
	            	} 
	            }
			}
		}
		
		
		//update H-field
		
		if (amplitudeTM) {	//1D-MAP TM
			for ( i = 0 ; i < mapN ; i++ ) {			//Compute H-field
		    	mapTMX[i] = mapTMX[i] + (mapSx*nr)*(mapTM0[i] - mapTM0[i+my])*0.5; //half step forward (final step completed at end)
				mapTMY[i] = mapTMY[i] + (mapSy*nr)*(mapTM0[i+mx] - mapTM0[i])*0.5; //half step forward (final step completed at end)
			}
		}
		if (amplitudeTE) { 	//1D-MAP TE
			for ( i = 0 ; i < mapN ; i++ ) {				//Update 1DMAP  TE
				mapTE2[i] = mapTE1[i];						//mapTE2 only used temporarely
				mapTE1[i] = mapTE0[i];	
			}
			
			for ( i = nMin ; i < nMax ; i++ ) {
				mapTE0[i] = -mapTE2[i] + 2*(1-mapSx*mapSx-mapSy*mapSy)*mapTE1[i] + (mapSx*mapSx)*(mapTE1[i+mx] + mapTE1[i-mx]) + (mapSy*mapSy)*(mapTE1[i+my] + mapTE1[i-my]);
			}
			
			//1D-MAP source
			for ( i = 0; i <nMin; i++) {
				mapTE0[(nMin-1)-i] = sourceSignal(src_f, src_df, amplitudeTE*nr*(nrc/nrc0), phaseTE, src_type, -i*dr, (T0+t)*dt/nr, nr);
			}
			
			for ( i = nMax; i <mapN; i++) {  //MUR Absorbing boundary
				mapTE0[i] = mapTE1[i-1] + (dt/nrc-dr)/(dt/nrc+dr)*(mapTE0[i-1] - mapTE1[i]);
			}
			
			for ( i = 0; i < mapN; i++) {
				mapTE2[i] = mapTE0[i];	//Store full step in mapTE2
				mapTE0[i] = (mapTE0[i] + mapTE1[i])/2; //mapTE0 now holds Hz half step
			}
		}


		
	    for (k=0; k<=sz; k++) {
	    	zz = ZZH[k];
	        for (j=0; j<=sx; j++) {
	            zx = ZXH[j];
	            for (i=0; i<=sy; i++) {
	                zy = ZYH[i];
					
					if ((i<sy)&&(k<sz)) {
		                ijk = _ijkHx;						
		                Hxy[ijk] = Hxy[ijk]*zy - (Ez[_dEz_dy]-Ez[_ijkEz])*dt;
		                Hxz[ijk] = Hxz[ijk]*zz + (Ey[_dEy_dz]-Ey[_ijkEy])*dt;
						
						Hx[ijk]  = (Hx[ijk] + Hxy[ijk] + Hxz[ijk])/2;
	                }
	                
	                if ((j<sx)&&(k<sz)) {
		                ijk = _ijkHy;						
		                Hyz[ijk] = Hyz[ijk]*zz - (Ex[_dEx_dz] - Ex[_ijkEx])*dt;
		                Hyx[ijk] = Hyx[ijk]*zx + (Ez[_dEz_dx] - Ez[_ijkEz])*dt;
				
		                Hy[ijk]  = (Hy[ijk] + Hyz[ijk] + Hyx[ijk])/2;
                	}
                	
                	if ((j<sx)&&(i<sy)) {
		                ijk = _ijkHz;						
		                Hzx[ijk] = Hzx[ijk]*zx - (Ey[_dEy_dx]-Ey[_ijkEy])*dt;
		            	Hzy[ijk] = Hzy[ijk]*zy + (Ex[_dEx_dy]-Ex[_ijkEx])*dt;

						Hz[ijk]  = (Hz[ijk] + Hzx[ijk] + Hzy[ijk])/2;
		            }            
	            }
	        }
		}
		
		//H now contains half step forward, at same time step as E. The sum of the two split field components contains the H at full step forward.
		
		
		if (amplitudeTM) updateTFSF_HTM(Hx, Hy, Hxy, Hyx, bc, tfsfBox, mapTM0, nx, ny, dr, signX, signY, sx, sy, sz);
		if (amplitudeTE) updateTFSF_HTE(Hx,Hxz,Hy,Hyz,Hz,Hzx,Hzy,bc,tfsfBox,mapTEX,mapTEY,nx,ny,dr,signX,signY,sx,sy,sz,dt);
		
		 
		//point probes
		if (amplitudeTE) {
			i = nMin;		//1DMAP
			PLWprb[t+0*N] = (mapTE0[i-1]+mapTE0[i])/2;	//Hz									//TE0 is half step of Hz (TE2 holds full step)
			PLWprb[t+1*N] = sqrt(mapTEX[i]*mapTEX[i] + mapTEY[i]*mapTEY[i])*((my<mx)?sign(mapTEY[i]):sign(mapTEX[i])); //Ex+Ey	
		}
		
		if (amplitudeTM) {
			i = nMin;		//1DMAP
			PLWprb[t+2*N] = (mapTM0[i]+mapTM0[i+1])/2; //Ez
			PLWprb[t+3*N] = sqrt(mapTMX[i]*mapTMX[i] + mapTMY[i]*mapTMY[i])*((my<mx)?sign(mapTMY[i]):sign(mapTMX[i])); //Hx+Hy TMX and TMY is half step of Hx and Hy
		}
		
		for ( int probe = 0 ; probe < PRBn ; probe++ ) { 		//3D
			i = 			(int)PRBparam[(probe+1)*PRBN-5];
			j = 			(int)PRBparam[(probe+1)*PRBN-4];
			k = 			(int)PRBparam[(probe+1)*PRBN-3];
			int prbtype = 	(int)PRBparam[(probe+1)*PRBN-2];
			int tfsf =      (int)PRBparam[(probe+1)*PRBN-1];
			int withinTF = 2;
			switch (prbtype) {
				case EX:
					if (amplitudeTE) withinTF = 1;
					PRB[t+probe*N] = getValueAtClusterPoint(i, j, k, Ex, EX, sx, sy, sz, tfsfBox, mapTE0, mapTEX, mapTEY, nx, ny, dr, signX, signY, tfsf, withinTF);
					break;
				case EY:
					if (amplitudeTE) withinTF = 1;
					PRB[t+probe*N] = getValueAtClusterPoint(i, j, k, Ey, EY, sx, sy, sz, tfsfBox, mapTE0, mapTEX, mapTEY, nx, ny, dr, signX, signY, tfsf, withinTF);
					break;
				case EZ:
					if (amplitudeTM) withinTF = 1;
					PRB[t+probe*N] = getValueAtClusterPoint(i, j, k, Ez, EZ, sx, sy, sz, tfsfBox, mapTM0, mapTMX, mapTMY, nx, ny, dr, signX, signY, tfsf, withinTF);
					break;
				case HX:
					if (amplitudeTM) withinTF = 1;
					PRB[t+probe*N] = getValueAtClusterPoint(i, j, k, Hx, HX, sx, sy, sz, tfsfBox, mapTM0, mapTMX, mapTMY, nx, ny, dr, signX, signY, tfsf, withinTF);
					break;
				case HY:
					if (amplitudeTM) withinTF = 1;
					PRB[t+probe*N] = getValueAtClusterPoint(i, j, k, Hy, HY, sx, sy, sz, tfsfBox, mapTM0, mapTMX, mapTMY, nx, ny, dr, signX, signY, tfsf, withinTF);
					break;
				case HZ:
					if (amplitudeTE) withinTF = 1;
					PRB[t+probe*N] = getValueAtClusterPoint(i, j, k, Hz, HZ, sx, sy, sz, tfsfBox, mapTE0, mapTEX, mapTEY, nx, ny, dr, signX, signY, tfsf, withinTF);
					break;
			}
	    }
		
		for ( int source = 0 ; source < SRCn ; source++ ) { 
			i = 			(int)SRCparam[(source+1)*SRCN-11];
			j = 			(int)SRCparam[(source+1)*SRCN-10];
			k = 			(int)SRCparam[(source+1)*SRCN-9];
			int srcType = 	(int)SRCparam[(source+1)*SRCN-8];
			int pntType = 	(int)SRCparam[(source+1)*SRCN-7];
			switch (srcType) {
				case EX:
					if (pntType == CLUSTERPOINT) SRCfield[t+source*N] = getRawValueAtClusterPoint(i,j,k, Ex, EX, sx, sy, sz);
					else if (pntType == SINGLEPOINT) SRCfield[t+source*N] = getRawValueAtSinglePoint(i,j,k, Ex, EX, sx, sy, sz);
					break;
				case EY:
					if (pntType == CLUSTERPOINT) SRCfield[t+source*N] = getRawValueAtClusterPoint(i,j,k, Ey, EY, sx, sy, sz);
					else if (pntType == SINGLEPOINT) SRCfield[t+source*N] = getRawValueAtSinglePoint(i,j,k, Ey, EY, sx, sy, sz);
					break;
				case EZ:
					if (pntType == CLUSTERPOINT) SRCfield[t+source*N] = getRawValueAtClusterPoint(i,j,k, Ez, EZ, sx, sy, sz);
					else if (pntType == SINGLEPOINT) SRCfield[t+source*N] = getRawValueAtSinglePoint(i,j,k, Ez, EZ, sx, sy, sz);
					break;
				case HX:
					if (pntType == CLUSTERPOINT) SRCfield[t+source*N] = getRawValueAtClusterPoint(i,j,k, Hx, HX, sx, sy, sz);
					else if (pntType == SINGLEPOINT) SRCfield[t+source*N] = getRawValueAtSinglePoint(i,j,k, Hx, HX, sx, sy, sz);
					break;
				case HY:
					if (pntType == CLUSTERPOINT) SRCfield[t+source*N] = getRawValueAtClusterPoint(i,j,k, Hy, HY, sx, sy, sz);
					else if (pntType == SINGLEPOINT) SRCfield[t+source*N] = getRawValueAtSinglePoint(i,j,k, Hy, HY, sx, sy, sz);
					break;
				case HZ:
					if (pntType == CLUSTERPOINT) SRCfield[t+source*N] = getRawValueAtClusterPoint(i,j,k, Hz, HZ, sx, sy, sz);
					else if (pntType == SINGLEPOINT) SRCfield[t+source*N] = getRawValueAtSinglePoint(i,j,k, Hz, HZ, sx, sy, sz);
					break;
				case EDP: {
						double nx = SRCparam[(source+1)*SRCN-6];
						double ny = SRCparam[(source+1)*SRCN-5];
						double nz = SRCparam[(source+1)*SRCN-4];
						SRCfield[t+source*N]  = nx*getRawValueAtClusterPoint(i,j,k, Ex, EX, sx, sy, sz);
						SRCfield[t+source*N] += ny*getRawValueAtClusterPoint(i,j,k, Ey, EY, sx, sy, sz);
						SRCfield[t+source*N] += nz*getRawValueAtClusterPoint(i,j,k, Ez, EZ, sx, sy, sz);
				} break;	
				case MDP: {
						double nx = SRCparam[(source+1)*SRCN-6];
						double ny = SRCparam[(source+1)*SRCN-5];
						double nz = SRCparam[(source+1)*SRCN-4];
						SRCfield[t+source*N]  = nx*getRawValueAtClusterPoint(i,j,k, Hx, HX, sx, sy, sz);
						SRCfield[t+source*N] += ny*getRawValueAtClusterPoint(i,j,k, Hy, HY, sx, sy, sz);
						SRCfield[t+source*N] += nz*getRawValueAtClusterPoint(i,j,k, Hz, HZ, sx, sy, sz);
					} break;	
				}
	    }
		
		
		//poynting boxes E x H
		int indx = 0;
		int PYBtfsfBoxIndx = 0;
		for ( int PYBbox = 0 ; PYBbox < PYBn ; PYBbox++ ) {			
			int posi =       (int)PYBparam[(PYBbox+1)*PYBN-11]; //Box positions
			int posj =       (int)PYBparam[(PYBbox+1)*PYBN-10];
			int posk =       (int)PYBparam[(PYBbox+1)*PYBN-9];
			int BSi  =    	 (int)PYBparam[(PYBbox+1)*PYBN-8];	//Box sizes
			int BSj  =    	 (int)PYBparam[(PYBbox+1)*PYBN-7];
			int BSk  =    	 (int)PYBparam[(PYBbox+1)*PYBN-6];
			int pybSize = 	 (int)PYBparam[(PYBbox+1)*PYBN-5];
			int startT = 	 (int)PYBparam[(PYBbox+1)*PYBN-4];
			int tfsf =       (int)PYBparam[(PYBbox+1)*PYBN-3];
			int interleave = (int)PYBparam[(PYBbox+1)*PYBN-2];
			int side0 = (int)PYBparam[(PYBbox+1)*PYBN+0];
			int side1 = (int)PYBparam[(PYBbox+1)*PYBN+1];
			int side2 = (int)PYBparam[(PYBbox+1)*PYBN+2];
			int side3 = (int)PYBparam[(PYBbox+1)*PYBN+3];
			int side4 = (int)PYBparam[(PYBbox+1)*PYBN+4];
			int side5 = (int)PYBparam[(PYBbox+1)*PYBN+5];
						
			double ex, ey, ez, hx, hy, hz;
			
			
			if (((T0+t)>=startT)&&(!((T0+t)%interleave)))  {
	        	for ( int fn = 0 ; fn < PYBfn ; fn++ ) {
		        	double f = PYBf[fn];
		        	double iSin0 = sin(f*dt*(T0+t));
					double rCos0 = cos(f*dt*(T0+t));
					double iSin, rCos, F;
					double withinTF; //-1 if out side, 1 inside and 0 part inside part outside
					int PYBtfsfIndx = 0;
										
		    		for (int b=-BSk; b<=BSk; b++) {
			    		for (int a=-BSj; a<=BSj; a++) {
			    			//constant y: (EzHx-ExHz)
				    		i = posi-BSi;
				    		j = posj+a;
				    		k = posk+b;
							
							if (abs(a)<BSj) F = 1; else F = 0.7071;
							if (abs(b)==BSk) F *= 0.7071;
							iSin = F*iSin0;
							rCos = F*rCos0;
							
				    		if ((side0==1)&&(i>=0)&&(j>=0)&&(k>=0)&&(i<sy)&&(j<sx)&&(k<sz)) {
								if (!amplitudeTE) withinTF = 2; else withinTF = PYBtfsf[PYBtfsfBoxIndx + PYBtfsfIndx];
								ex = getValueAtClusterPoint(i, j, k, Ex, EX, sx, sy, sz, tfsfBox, mapTE0, mapTEX, mapTEY, nx, ny, dr, signX, signY, tfsf, withinTF);
								hz = getValueAtClusterPoint(i, j, k, Hz, HZ, sx, sy, sz, tfsfBox, mapTE0, mapTEX, mapTEY, nx, ny, dr, signX, signY, tfsf, withinTF);
								if (!amplitudeTM) withinTF = 2; else withinTF = PYBtfsf[PYBtfsfBoxIndx + PYBtfsfIndx];
								ez = getValueAtClusterPoint(i, j, k, Ez, EZ, sx, sy, sz, tfsfBox, mapTM0, mapTMX, mapTMY, nx, ny, dr, signX, signY, tfsf, withinTF);
								hx = getValueAtClusterPoint(i, j, k, Hx, HX, sx, sy, sz, tfsfBox, mapTM0, mapTMX, mapTMY, nx, ny, dr, signX, signY, tfsf, withinTF);
								PYBtfsfIndx++;
					    		PYBr[indx] -= ez*rCos; PYBi[indx++] -= ez*iSin;
					    		PYBr[indx] += hx*rCos; PYBi[indx++] += hx*iSin;
					    		PYBr[indx] -= ex*rCos; PYBi[indx++] -= ex*iSin;
					    		PYBr[indx] += hz*rCos; PYBi[indx++] += hz*iSin;
				    		} 
							
				    		i = posi+BSi;
				    		if ((side1==1)&&(i>=0)&&(j>=0)&&(k>=0)&&(i<sy)&&(j<sx)&&(k<sz)) {
								if (!amplitudeTE) withinTF = 2; else withinTF = PYBtfsf[PYBtfsfBoxIndx + PYBtfsfIndx];
								ex = getValueAtClusterPoint(i, j, k, Ex, EX, sx, sy, sz, tfsfBox, mapTE0, mapTEX, mapTEY, nx, ny, dr, signX, signY, tfsf, withinTF);
								hz = getValueAtClusterPoint(i, j, k, Hz, HZ, sx, sy, sz, tfsfBox, mapTE0, mapTEX, mapTEY, nx, ny, dr, signX, signY, tfsf, withinTF);
								if (!amplitudeTM) withinTF = 2; else withinTF = PYBtfsf[PYBtfsfBoxIndx + PYBtfsfIndx];
								ez = getValueAtClusterPoint(i, j, k, Ez, EZ, sx, sy, sz, tfsfBox, mapTM0, mapTMX, mapTMY, nx, ny, dr, signX, signY, tfsf, withinTF);
								hx = getValueAtClusterPoint(i, j, k, Hx, HX, sx, sy, sz, tfsfBox, mapTM0, mapTMX, mapTMY, nx, ny, dr, signX, signY, tfsf, withinTF);
								PYBtfsfIndx++;
					    		PYBr[indx] += ez*rCos; PYBi[indx++] += ez*iSin;
					    		PYBr[indx] += hx*rCos; PYBi[indx++] += hx*iSin;
					    		PYBr[indx] += ex*rCos; PYBi[indx++] += ex*iSin;
					    		PYBr[indx] += hz*rCos; PYBi[indx++] += hz*iSin;
				    		}
				    	}
				    }
				    for (int b=-BSk; b<=BSk; b++) {
			    		for (int a=-BSi; a<=BSi; a++) {
				    		//constant x: (EyHz-EzHy)
				    		i = posi+a;
				    		j = posj-BSj;
				    		k = posk+b;
							
							if (abs(a)<BSi) F = 1; else F = 0.7071;
							if (abs(b)==BSk) F *= 0.7071;
							iSin = F*iSin0;
							rCos = F*rCos0;
												
				    		if ((side2==1)&&(i>=0)&&(j>=0)&&(k>=0)&&(i<sy)&&(j<sx)&&(k<sz)) {
								if (!amplitudeTE) withinTF = 2; else withinTF = PYBtfsf[PYBtfsfBoxIndx + PYBtfsfIndx];
								ey = getValueAtClusterPoint(i, j, k, Ey, EY, sx, sy, sz, tfsfBox, mapTE0, mapTEX, mapTEY, nx, ny, dr, signX, signY, tfsf, withinTF);
								hz = getValueAtClusterPoint(i, j, k, Hz, HZ, sx, sy, sz, tfsfBox, mapTE0, mapTEX, mapTEY, nx, ny, dr, signX, signY, tfsf, withinTF);
								if (!amplitudeTM) withinTF = 2; else withinTF = PYBtfsf[PYBtfsfBoxIndx + PYBtfsfIndx];
								ez = getValueAtClusterPoint(i, j, k, Ez, EZ, sx, sy, sz, tfsfBox, mapTM0, mapTMX, mapTMY, nx, ny, dr, signX, signY, tfsf, withinTF);
								hy = getValueAtClusterPoint(i, j, k, Hy, HY, sx, sy, sz, tfsfBox, mapTM0, mapTMX, mapTMY, nx, ny, dr, signX, signY, tfsf, withinTF);
					    		PYBtfsfIndx++;
								PYBr[indx] -= ey*rCos; PYBi[indx++] -= ey*iSin;
					    		PYBr[indx] += hz*rCos; PYBi[indx++] += hz*iSin;
					    		PYBr[indx] -= ez*rCos; PYBi[indx++] -= ez*iSin;
					    		PYBr[indx] += hy*rCos; PYBi[indx++] += hy*iSin;
							} 
							
				    		j = posj+BSj;
				    		if ((side3==1)&&(i>=0)&&(j>=0)&&(k>=0)&&(i<sy)&&(j<sx)&&(k<sz)) {
								if (!amplitudeTE) withinTF = 2; else withinTF = PYBtfsf[PYBtfsfBoxIndx + PYBtfsfIndx];
								ey = getValueAtClusterPoint(i, j, k, Ey, EY, sx, sy, sz, tfsfBox, mapTE0, mapTEX, mapTEY, nx, ny, dr, signX, signY, tfsf, withinTF);
								hz = getValueAtClusterPoint(i, j, k, Hz, HZ, sx, sy, sz, tfsfBox, mapTE0, mapTEX, mapTEY, nx, ny, dr, signX, signY, tfsf, withinTF);
								if (!amplitudeTM) withinTF = 2; else withinTF = PYBtfsf[PYBtfsfBoxIndx + PYBtfsfIndx];
								ez = getValueAtClusterPoint(i, j, k, Ez, EZ, sx, sy, sz, tfsfBox, mapTM0, mapTMX, mapTMY, nx, ny, dr, signX, signY, tfsf, withinTF);
								hy = getValueAtClusterPoint(i, j, k, Hy, HY, sx, sy, sz, tfsfBox, mapTM0, mapTMX, mapTMY, nx, ny, dr, signX, signY, tfsf, withinTF);
								PYBtfsfIndx++;
					    		PYBr[indx] += ey*rCos; PYBi[indx++] += ey*iSin;
					    		PYBr[indx] += hz*rCos; PYBi[indx++] += hz*iSin;
					    		PYBr[indx] += ez*rCos; PYBi[indx++] += ez*iSin;
					    		PYBr[indx] += hy*rCos; PYBi[indx++] += hy*iSin;
				    		} 
				    	}
				    }
				    
				    for (int b=-BSj; b<=BSj; b++) {
			    		for (int a=-BSi; a<=BSi; a++) {
				    		//constant z: (ExHy-EyHx)
				    		i = posi+a;
				    		j = posj+b;
				    		k = posk-BSk;
							
							if (abs(a)<BSi) F = 1; else F = 0.7071;
							if (abs(b)==BSj) F *= 0.7071;
							iSin = F*iSin0;
							rCos = F*rCos0;
							
				    		if ((side4==1)&&(i>=0)&&(j>=0)&&(k>=0)&&(i<sy)&&(j<sx)&&(k<sz)) {
								if (!amplitudeTE) withinTF = 2; else withinTF = PYBtfsf[PYBtfsfBoxIndx + PYBtfsfIndx];
								ex = getValueAtClusterPoint(i, j, k, Ex, EX, sx, sy, sz, tfsfBox, mapTE0, mapTEX, mapTEY, nx, ny, dr, signX, signY, tfsf, withinTF);
								ey = getValueAtClusterPoint(i, j, k, Ey, EY, sx, sy, sz, tfsfBox, mapTE0, mapTEX, mapTEY, nx, ny, dr, signX, signY, tfsf, withinTF);
								if (!amplitudeTM) withinTF = 2; else withinTF = PYBtfsf[PYBtfsfBoxIndx + PYBtfsfIndx];
								hx = getValueAtClusterPoint(i, j, k, Hx, HX, sx, sy, sz, tfsfBox, mapTM0, mapTMX, mapTMY, nx, ny, dr, signX, signY, tfsf, withinTF);
								hy = getValueAtClusterPoint(i, j, k, Hy, HY, sx, sy, sz, tfsfBox, mapTM0, mapTMX, mapTMY, nx, ny, dr, signX, signY, tfsf, withinTF);
					    		PYBtfsfIndx++;
								PYBr[indx] -= ex*rCos; PYBi[indx++] -= ex*iSin;
					    		PYBr[indx] += hy*rCos; PYBi[indx++] += hy*iSin;
					    		PYBr[indx] -= ey*rCos; PYBi[indx++] -= ey*iSin;
					    		PYBr[indx] += hx*rCos; PYBi[indx++] += hx*iSin;
				    		}
							
				    		k = posk+BSk;
				    		if ((side5==1)&&(i>=0)&&(j>=0)&&(k>=0)&&(i<sy)&&(j<sx)&&(k<sz)) {
								if (!amplitudeTE) withinTF = 2; else withinTF = PYBtfsf[PYBtfsfBoxIndx + PYBtfsfIndx];
								ex = getValueAtClusterPoint(i, j, k, Ex, EX, sx, sy, sz, tfsfBox, mapTE0, mapTEX, mapTEY, nx, ny, dr, signX, signY, tfsf, withinTF);
								ey = getValueAtClusterPoint(i, j, k, Ey, EY, sx, sy, sz, tfsfBox, mapTE0, mapTEX, mapTEY, nx, ny, dr, signX, signY, tfsf, withinTF);
								if (!amplitudeTM) withinTF = 2; else withinTF = PYBtfsf[PYBtfsfBoxIndx + PYBtfsfIndx];
								hx = getValueAtClusterPoint(i, j, k, Hx, HX, sx, sy, sz, tfsfBox, mapTM0, mapTMX, mapTMY, nx, ny, dr, signX, signY, tfsf, withinTF);
								hy = getValueAtClusterPoint(i, j, k, Hy, HY, sx, sy, sz, tfsfBox, mapTM0, mapTMX, mapTMY, nx, ny, dr, signX, signY, tfsf, withinTF);
					    		PYBtfsfIndx++;
								PYBr[indx] += ex*rCos; PYBi[indx++] += ex*iSin;
					    		PYBr[indx] += hy*rCos; PYBi[indx++] += hy*iSin;
					    		PYBr[indx] += ey*rCos; PYBi[indx++] += ey*iSin;
					    		PYBr[indx] += hx*rCos; PYBi[indx++] += hx*iSin;
				    		}
			    		}
			    	}					
				}
				PYBparam[(PYBbox+1)*PYBN-1]++; //Count up
			}
			PYBtfsfBoxIndx += pybSize;
		}
			
		
		//DFT
		int dftDataIndex = 0;
		for ( int dft = 0 ; dft < DFTn ; dft++ ) {
			//frequency dimension plane fieldtype sizex sizey
			double f = 		DFTparam[(dft+1)*DFTN-10];
			int dim =  (int)DFTparam[(dft+1)*DFTN-9];
			int indx = (int)DFTparam[(dft+1)*DFTN-8];
			int fieldType = (int)DFTparam[(dft+1)*DFTN-7];
			int sizea = (int)DFTparam[(dft+1)*DFTN-6];
			int sizeb = (int)DFTparam[(dft+1)*DFTN-5];
			int startT = (int)DFTparam[(dft+1)*DFTN-4];
			int tfsf = (int)DFTparam[(dft+1)*DFTN-3];
			int interleave = (int)DFTparam[(dft+1)*DFTN-2];
			if (((T0+t)>=startT)&&(!((T0+t)%interleave))) {
				double fieldValue;
				double iSin = sin(f*dt*(T0+t));
				double rCos = cos(f*dt*(T0+t));
				for (int b=0; b<sizeb; b++) {
					for (int a=0; a<sizea; a++) {
						switch (dim) {
							case 1: 
								i = indx;
				    			j = a;
				    			k = b;
				    			break;
				    		case 2:
				    			i = a;
				    			j = indx;
				    			k = b;
				    			break;
				    		case 3:
				    			i = a;
				    			j = b;
				    			k = indx;
				    			break;
						}
						int withinTF = 1; //Check always - default
						switch (fieldType) {
							case EX:
								if (!amplitudeTE) withinTF = 2;
								fieldValue = getValueAtClusterPoint(i, j, k, Ex, EX, sx, sy, sz, tfsfBox, mapTE0, mapTEX, mapTEY, nx, ny, dr, signX, signY, tfsf, withinTF);
								break;
							case EY:
								if (!amplitudeTE) withinTF = 2;
								fieldValue = getValueAtClusterPoint(i, j, k, Ey, EY, sx, sy, sz, tfsfBox, mapTE0, mapTEX, mapTEY, nx, ny, dr, signX, signY, tfsf, withinTF);
								break;
							case EZ:
								if (!amplitudeTM) withinTF = 2;
								fieldValue = getValueAtClusterPoint(i, j, k, Ez, EZ, sx, sy, sz, tfsfBox, mapTM0, mapTMX, mapTMY, nx, ny, dr, signX, signY, tfsf, withinTF);
								break;
							case HX:
								if (!amplitudeTM) withinTF = 2;
								fieldValue = getValueAtClusterPoint(i, j, k, Hx, HX, sx, sy, sz, tfsfBox, mapTM0, mapTMX, mapTMY, nx, ny, dr, signX, signY, tfsf, withinTF);
								break;
							case HY:
								if (!amplitudeTM) withinTF = 2;
								fieldValue = getValueAtClusterPoint(i, j, k, Hy, HY, sx, sy, sz, tfsfBox, mapTM0, mapTMX, mapTMY, nx, ny, dr, signX, signY, tfsf, withinTF);
								break;
							case HZ:
								if (!amplitudeTE) withinTF = 2;
								fieldValue = getValueAtClusterPoint(i, j, k, Hz, HZ, sx, sy, sz, tfsfBox, mapTE0, mapTEX, mapTEY, nx, ny, dr, signX, signY, tfsf, withinTF);
								break;
						}
						DFTr[a+b*sizea+dftDataIndex] += rCos*fieldValue;
						DFTi[a+b*sizea+dftDataIndex] += iSin*fieldValue;
					}
				}
				DFTparam[(dft+1)*DFTN-1]++;  //Count up
			}
			dftDataIndex += sizea*sizeb;
		}    
    
		//Only at last iteration
		if (t==(N-1)) {							//Record data for live-view preview
			int A = 1, B = 1;
			switch (playDim) {
			case 1: 
				A = sx; //a = j;
				B = sz; //b = k;
				break;
			case 2:
				A = sy; //a = i;
				B = sz; //b = k;
				break;
			case 3:
				A = sy;	//a = i
				B = sx; //b = j;
				break;
			}
						
			for (int b=0; b<B; b++) {			//2D view
		    	for (int a=0; a<A; a++) {
					switch (playDim) {
						case 1: 
							i = playIndx;
			    			j = a;
			    			k = b;
			    			break;
			    		case 2:
			    			i = a;
			    			j = playIndx;
			    			k = b;
			    			break;
			    		case 3:
			    			i = a;
			    			j = b;
			    			k = playIndx;
			    			break;
					}
					double fieldValue;
					int withinTF = 2;  //Raw default
					switch (playFieldType) {
						case EX:
							if (amplitudeTE) withinTF = 1;
							fieldValue = getValueAtClusterPoint(i, j, k, Ex, EX, sx, sy, sz, tfsfBox, mapTE0, mapTEX, mapTEY, nx, ny, dr, signX, signY, playTFSF, withinTF);
							break;
						case EY:
							if (amplitudeTE) withinTF = 1;
							fieldValue = getValueAtClusterPoint(i, j, k, Ey, EY, sx, sy, sz, tfsfBox, mapTE0, mapTEX, mapTEY, nx, ny, dr, signX, signY, playTFSF, withinTF);
							break;
						case EZ:
							if (amplitudeTM) withinTF = 1;
							fieldValue = getValueAtClusterPoint(i, j, k, Ez, EZ, sx, sy, sz, tfsfBox, mapTM0, mapTMX, mapTMY, nx, ny, dr, signX, signY, playTFSF, withinTF);
							break;
						case HX:
							if (amplitudeTM) withinTF = 1;
							fieldValue = getValueAtClusterPoint(i, j, k, Hx, HX, sx, sy, sz, tfsfBox, mapTM0, mapTMX, mapTMY, nx, ny, dr, signX, signY, playTFSF, withinTF);
							break;
						case HY:
							if (amplitudeTM) withinTF = 1;
							fieldValue = getValueAtClusterPoint(i, j, k, Hy, HY, sx, sy, sz, tfsfBox, mapTM0, mapTMX, mapTMY, nx, ny, dr, signX, signY, playTFSF, withinTF);
							break;
						case HZ:
							if (amplitudeTE) withinTF = 1;
							fieldValue = getValueAtClusterPoint(i, j, k, Hz, HZ, sx, sy, sz, tfsfBox, mapTE0, mapTEX, mapTEY, nx, ny, dr, signX, signY, playTFSF, withinTF);
							break;
					}
		        	play[a+b*A] =  fieldValue;
		        }
		    }
		}
		
		
		//save H-field advanced a full step
		if (amplitudeTM) {	//1D-MAP TM
			for ( i = 0 ; i < mapN ; i++ ) {			//Compute H-field
		    	mapTMX[i] = mapTMX[i] + (mapSx*nr)*(mapTM0[i] - mapTM0[i+my])*0.5; //final half step forward
				mapTMY[i] = mapTMY[i] + (mapSy*nr)*(mapTM0[i+mx] - mapTM0[i])*0.5; //final half step forward
			}
		}
		
		if (amplitudeTE) {
			for ( i = 0; i < mapN; i++) {
				mapTE0[i] = mapTE2[i];	//Restore full step in mapTE0
			}
		}
			
	    for (k=0; k<=sz; k++) {
	        for (j=0; j<=sx; j++) {
	            for (i=0; i<=sy; i++) {
					if ((i<sy)&&(k<sz)) {
		                ijk = _ijkHx;						
						Hx[ijk] = Hxy[ijk] + Hxz[ijk];
	                }
	                
	                if ((j<sx)&&(k<sz)) {
		                ijk = _ijkHy;												
		                Hy[ijk]  = Hyz[ijk] + Hyx[ijk];
	            	}
	            	
	            	if ((j<sx)&&(i<sy)) {
		                ijk = _ijkHz;						
						Hz[ijk]  = Hzx[ijk] + Hzy[ijk];
		            }            
	            }
	        }
		}//H now contains full step forward, at half a time-step ahead of E. 
		
	
	gTime++; //Increment global time step
	
	} //Main loop end
	
    options[9] = gTime;	//Update global time step
	
}




void mexFunction(int nlhs, mxArray *plhs[],       /* Output variables */
                 int nrhs, const mxArray *prhs[]) /* Input variables  */
{   
    int sx, sy, sz, N, T0, SRCn, PRBn, PMLn, PYBn, PYBfn, DFTn, MPARAMn;
    double *Ex, *Exy, *Exz, *Ey, *Eyx, *Eyz, *Ez, *Ezx, *Ezy, *Hx, *Hxy, *Hxz, *Hy, *Hyx, *Hyz, *Hz, *Hzx, *Hzy, *SRC, *SRCfield, *PRB, *PML, *PLWprb;
    double *PYBr, *PYBi, *DFTr, *DFTi, *DFTparam, *PYBparam, *PYBf,  *SRCsignal, *SRCparam, *PRBparam, *MPARAM, *play;
    double *Jxy, *Jxz, *Jyx, *Jyz, *Jzx, *Jzy;
    unsigned char *MIndx;
    double  *PLWparam, *options, *bc, *tfsfBox;
    const mwSize *dims;
    struct SplitField E, H, J;
	int mapN;
	double  *mapTE0, *mapTE1, *mapTE2, *mapTEX, *mapTEY;
	double  *mapTM0, *mapTM1, *mapTM2, *mapTMX, *mapTMY;
	double  *PYBtfsf;
	
    if (nrhs>0) {
        options = mxGetPr(prhs[0]);
        tfsfBox = mxGetPr(prhs[1]);
        bc = mxGetPr(prhs[2]);
        dims = mxGetDimensions(prhs[3]);
        sy = dims[0]-1; sx = dims[1]; sz = dims[2]-1;
        E.x = mxGetPr(prhs[3]);
        E.xy= mxGetPr(prhs[4]);
        E.xz= mxGetPr(prhs[5]);
        E.y = mxGetPr(prhs[6]);
        E.yx= mxGetPr(prhs[7]);
        E.yz= mxGetPr(prhs[8]);
        E.z = mxGetPr(prhs[9]);
        E.zx= mxGetPr(prhs[10]);
        E.zy= mxGetPr(prhs[11]);
        H.x = mxGetPr(prhs[12]);
        H.xy= mxGetPr(prhs[13]);
        H.xz= mxGetPr(prhs[14]);
        H.y = mxGetPr(prhs[15]);
        H.yx= mxGetPr(prhs[16]);
        H.yz= mxGetPr(prhs[17]);
        H.z = mxGetPr(prhs[18]);
        H.zx= mxGetPr(prhs[19]);
        H.zy= mxGetPr(prhs[20]);
        N = mxGetScalar(prhs[21]);
        T0 = mxGetScalar(prhs[22]);
		SRCsignal = mxGetPr(prhs[23]);
        SRCparam = mxGetPr(prhs[24]);	//pointsources
        SRCn = mxGetN(prhs[24]); //Same
        SRC = mxGetPr(prhs[25]);
		SRCfield = mxGetPr(prhs[26]);
        PRBparam = mxGetPr(prhs[27]);	//probes
        PRBn = mxGetN(prhs[27]); //Same
        PRB = mxGetPr(prhs[28]);
        PML= mxGetPr(prhs[29]);								//PML
        PMLn = (mxGetM(prhs[29])+1)/2; //Same
        PYBparam = mxGetPr(prhs[30]);	//poynting box
        PYBn = mxGetN(prhs[30]); //Same
        PYBf = mxGetPr(prhs[31]);	//poynting box frequencies
        PYBfn= mxGetM(prhs[31]);	//same number of frequencies
        PYBr = mxGetPr(prhs[32]); //real part
        PYBi = mxGetPi(prhs[32]); //Same imaginary part
        DFTparam = mxGetPr(prhs[33]);	//poynting box
        DFTn = mxGetN(prhs[33]); //Same
        DFTr = mxGetPr(prhs[34]); //real part
        DFTi = mxGetPi(prhs[34]); //Same imaginary part
        
        J.xy= mxGetPr(prhs[35]);
        J.xz= mxGetPr(prhs[36]);
        J.yx= mxGetPr(prhs[37]);
        J.yz= mxGetPr(prhs[38]);
        J.zx= mxGetPr(prhs[39]);
        J.zy= mxGetPr(prhs[40]);
        
        MPARAM = mxGetPr(prhs[41]);
        MPARAMn = mxGetM(prhs[41]); //Same
        
        MIndx = (unsigned char*) mxGetData(prhs[42]);
        PLWparam = mxGetPr(prhs[43]);
        
        PLWprb = mxGetPr(prhs[44]);
        play = mxGetPr(prhs[45]);
		mapN = mxGetM(prhs[46])/2;
		mapTE0 = mxGetPr(prhs[46]);
		mapTE1 = mxGetPr(prhs[47]);
		mapTE2 = mxGetPr(prhs[48]);
		mapTEX = mxGetPr(prhs[49]);
		mapTEY = mxGetPr(prhs[50]);
		mapTM0 = &mapTE0[mapN];
		mapTM1 = &mapTE1[mapN];
		mapTM2 = &mapTE2[mapN];
		mapTMX = &mapTEX[mapN];
		mapTMY = &mapTEY[mapN];
		PYBtfsf = mxGetPr(prhs[51]);
		
		
        mexPrintf(".");
        //if (!((1+T0/(N))%50)) mexPrintf("\n");
        
        cfdtd3d_pml(options, bc,
        		sx, sy, sz,
				E, H, J,
				N, T0,
				SRCsignal, SRCparam, SRCn, SRC, SRCfield,
				PRBparam, PRBn, PRB,
				PMLn, PML, 
				PYBparam, PYBn, PYBr, PYBi, PYBf, PYBfn,
				DFTparam, DFTn, DFTr, DFTi,
				MIndx, MPARAMn, MPARAM,
				PLWparam, PLWprb, tfsfBox,
				play, mapN, mapTE0, mapTE1, mapTE2, mapTEX, mapTEY, mapTM0, mapTM1, mapTM2, mapTMX, mapTMY, PYBtfsf);
    }
    return;
}
