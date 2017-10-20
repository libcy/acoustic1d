#include <iostream>
#include "forward.h"

using namespace std;

int main() {

	double *c0 = mat::create(nx);
	for (int i = 0; i < nx; i++) {
		c0[i] = (double)1 + sin((double)i / 10)*sin((double)i / 20) / 5;
	}

	int nrx[nrt];
	int nsx[nst];
	for (int i = 0; i < nrt; i++) {
		nrx[i] = (int)(nx*(i + 1) / (nrt + 1) + 0.5);
	}
	for (int i = 0; i < nst; i++) {
		nsx[i] = (int)(nx*(i + 1) / (nst + 1) + 0.5);
	}
	nrx[0] = (int)(nx / 3 + 0.5) - 1;
	nsx[0] = (int)(nx * 2 / 3 + 0.5) - 1;

	int ns = nsx[0];
	double *src = mat::create(nt*nx, 0);

	struct fwdpar par;
	par.nrx = nrx;
	par.nsx = nsx;
	par.nr = nrx[0];
	par.ns = nsx[0];
	par.src = mat::create(nt*nx);
	par.type = 1;

	par.p = mat::create(nx);
	par.pold = mat::create(nx);
	par.pnew = mat::create(nx);
	par.d2p = mat::create(nx);
	par.gg = mat::create(nx);
	par.pwave = mat::create(nx*nt);

	for (int i = 0; i < nx; i++) {
		if (i > nx - nxa) {
			par.gg[i] = exp(-(apara * (i - nx + nxa)) * (apara * (i - nx + nxa)));
		}
		else {
			par.gg[i] = 1;
		}
	}

	double *p0x[nst];
	for (int i = 0; i < nst; i++) {
		par.ns = par.nsx[i];
		forward(c0, par);
		p0x[i] = mat::copy(par.pwave, nx*nt);
	}
	par.p0x = p0x;

	double *cx[2];
	cx[0] = c0;

	int tryn[] = { 7, 13, 25, 49, 97 };
	cg::val iresult;
	for (int i = 0; i < 5; i++) {
		par.nxd = tryn[i];
		double *c = mat::create(par.nxd, 0.9);

		if (i == 0) {
			for (int j = 0; j < par.nxd; j++) {
				c[j] = 0.9;
			}
		}
		else {
			double *tmpc = getc(iresult.x, tryn[i - 1]);
			for (int j = 0; j < par.nxd; j++) {
				c[j] = tmpc[(int)((double)j*(nx - 1) / (par.nxd - 1) + 0.5)];
			}
			free(tmpc);
			free(iresult.x);
			mat::plot(c, par.nxd, i + 9);
		}
		iresult = cg::calc(&gradient_x, par, c, par.nxd, 2, 2, 1, 0, 0, 200, 0, 0, 0);
		mat::plot(iresult.x, par.nxd, i);
	}

	cx[1] = getc(iresult.x, par.nxd);
	mat::plot(cx, 2, nx);

	//system("pause");
	return 0;
}