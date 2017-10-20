#pragma once

#include "const.h"
#include "cg.h"

using namespace std;

void forward(double *c, fwdpar par) {
	if (par.type) {
		for (int i = 0; i < nt; i++) {
			for (int j = 0; j < nx; j++) {
				int ij = i*nx + j;
				double tt = (double)(i + 1)*dt - t0;
				if (j == par.ns) {
					par.src[ij] = (1 - 2 * tt*tt*PI*PI*f0*f0)*exp(-PI*PI*f0*f0*tt*tt)*dt;
				}
				else {
					par.src[ij] = 0;
				}
			}
		}
	}
	for (int i = 0; i < nx; i++) {
		par.p[i] = 0;
		par.pnew[i] = 0;
		par.pold[i] = 0;
		par.d2p[i] = 0;
		for (int j = 0; j < nt; j++) {
			par.pwave[i*nt + j] = 0;
		}
	}
	
	for (int i = 0; i < nt; i++) {
		for (int j = 1; j < nx - 1; j++) {
			par.d2p[j] = (par.p[j + 1] - 2 * par.p[j] + par.p[j - 1]) / dx / dx;
			par.d2p[j] = par.d2p[j] *c[j];
		}

		for (int j = 0; j < nx; j++) {
			par.pnew[j] = 2 * par.p[j] - par.pold[j] + par.d2p[j] * dt*dt;
		}

		for (int j = 0; j < nx; j++) {
			par.pnew[j] += par.src[i*nx + j] * dt*dt;
			par.pold[j] = par.p[j] * par.gg[j];
			par.p[j] = par.pnew[j] * par.gg[j];
		}

		par.p[0] = 0;
		par.p[nx - 1] = 0;

		for (int j = 2; j < nx - 1; j++) {
			par.pwave[j*nt + i] = par.p[j];
		}
	}
}

double *getc(double *cd,int nxd) {
	int *points = mat::createInt(nxd);
	for (int i = 0; i < nxd; i++) {
		points[i] = (int)((double)i*(nx - 1) / (nxd - 1) + 0.5);
	}

	double *c = mat::create(nx);
	for (int i = 0; i < nxd - 1; i++) {
		int dtx = points[i + 1] - points[i];
		double dty = cd[i + 1] - cd[i];
		for (int j = points[i]; j <= points[i + 1]; j++) {
			c[j] = cd[i] + dty / dtx*(j - points[i]);
		}
	}
	return c;
}
cg::fwdval gradient(double *cd, struct fwdpar par) {
	int nxd = par.nxd;
	double *c = getc(cd, nxd);
	par.type = 1;
	forward(c, par);

	// 在此测试并入par的效率
	double *ssr = mat::create(nx*nx,0);
	ssr[par.nr*nx + par.nr] = 1;
	double *ssrt = mat::transpose(ssr, nx, nx);
	double *srcx = mat::multiply(ssrt,ssr,nx,nx,nx);
	double *p = mat::copy(par.pwave, nx*nt);
	double *dp = mat::add(par.pwave, par.p0, 1, -1, nx*nt);
	double g = 0;

	srcx = mat::multiply(srcx, dp, nx, nt, nx);

	for (int i = 0; i < nt; i++) {
		for (int j = 0; j < nx; j++) {
			if (j == par.nr) {
				par.src[i*nx + j] = srcx[par.nr*nt + nt - 1 - i];
			}
			else {
				par.src[i*nx + j] = 0;
			}
		}
		double dg = par.pwave[par.nr*nt + i] - par.p0[par.nr*nt + i];
		g += dg*dg;
	}

	g *= dt;
	double *gp = mat::create(nx,0);
	par.type = 0;
	forward(c, par);
	for (int i = 0; i < nt; i++) {
		for (int j = 1; j < nx - 1; j++) {
			int ii = nt - 2 - i;
			if (ii < 0){
				ii = 0;
			}
			gp[j] += par.pwave[j*nt + ii] * (p[(j + 1)*nt + i] - 2 * p[j*nt + i] + p[(j - 1)*nt + i])*dt / dx / dx;
		}
	}

	free(dp);
	free(ssr);
	free(ssrt);
	free(p);
	free(srcx);

	double *ccd = mat::create(nx*nxd,0);
	int *points = mat::createInt(nxd);
	for (int i = 0; i < nxd; i++) {
		points[i] = (int)((double)i*(nx - 1) / (nxd - 1) + 0.5);
	}
	for (int i = 0; i < nxd - 1; i++) {
		int dtx = points[i + 1] - points[i];
		for (int j = points[i]; j <= points[i + 1]; j++) {
			ccd[j*nxd + i] = 1 - (double)(j - points[i]) / dtx;
		}
	}
	for (int i = 1; i < nxd; i++) {
		int dtx = points[i] - points[i-1];
		for (int j = points[i-1]; j <= points[i]; j++) {
			ccd[j*nxd + i] = (double)(j - points[i-1]) / dtx;
		}
	}

	double *gpd = mat::multiply(gp, ccd, 1, nxd, nx);

	return{ g,gpd};
}
cg::fwdval gradient_x(double *cd, struct fwdpar par) {
	double g = 0;
	double *gpd = mat::create(par.nxd,0);
	for (int i = 0; i < nst; i++) {
		par.ns = par.nsx[i];
		par.p0 = par.p0x[i];
		for (int j = 0; j < nrt; j++) {
			//nr个接收点只需一次正演
			par.nr = par.nrx[j];
			cg::fwdval result = gradient(cd, par);
			g += result.f;
			mat::add(gpd, result.g, gpd, par.nxd);
			free(result.g);
		}
	}
	return{ g,gpd };
}

double getg(double *p, double *p0, int nr){
	double g = 0;
	for (int i = 0; i < nt; i++){
		int ij = nr*nt + i;
		double dg = p[ij] - p0[ij];
		g += dg*dg/2;
	}
	g *= dt;
	return g;
}

cg::fwdval gradient_d(double *cd, struct fwdpar par) {
	double g = 0;
	double *gpd = mat::create(par.nxd, 0);
	double dc = 0.001;

	for (int i = 0; i < nst; i++) {
		par.ns = par.nsx[i];
		par.p0 = par.p0x[i];
		for (int j = 0; j < nrt; j++) {
			par.nr = par.nrx[j];
			forward(getc(cd,par.nxd), par);
			double gg = getg(par.pwave, par.p0, par.nr);
			g += gg;

			for (int k = 0; k < par.nxd; k++){
				cd[k] += dc;
				forward(getc(cd, par.nxd), par);
				gpd[k] += (getg(par.pwave, par.p0, par.nr)-gg)/dc;
				cd[k] -= dc;
			}
		}
	}
	return{ g, gpd };
}