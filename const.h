#pragma once

const double PI = acos(-1);
const int nx = 250;
const double dx = (double)1 / (nx - 1);
const double cmax = 2;
const double cmin = 0.5;

const int nrt = 1;
const int nst = 1;

const double t = 3.0;
const double dt = 0.5*dx / 1.0;
const int nt = (int)(t / dt + 0.5);
const double f0 = cmin / (10 * dx);
const double t0 = 4 / f0;

const double apara = 0.015;
const int nxa = 20;



struct fwdpar {
	int *nrx;
	int *nsx;
	int nr;
	int ns;
	int nxd;
	double *src;
	int type;

	double *p;
	double *pold;
	double *pnew;
	double *d2p;
	double *gg;

	double *pwave;
	double *p0;
	double **p0x;
};
