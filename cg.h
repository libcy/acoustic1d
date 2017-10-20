#pragma once

#include "mat.h"

namespace cg {
	struct val {
		double *x;
		double f;
		double ginf;
		double step;
		int nstep;
		int neval;
		int stop;
	};
	struct lsval {
		double alpha;
		double fn;
		double *gn;
		int neval;
	};
	struct fwdval {
		double f;
		double *g;
	};
	struct cg::lsval ls(cg::fwdval(*func)(double *, struct fwdpar), struct fwdpar fpar, double *x, double f, double *g, double *h, int n,
		double alpha_max, int ls_method, double ls_gradient, double ls_step, int ls_neval) {
		double alpha = 0, fn = f;
		double *gn = mat::copy(g, n);
		int neval = 1;
		double dfi0 = mat::dot(h, gn, n);
		double fi0 = f, slope0 = ls_step*dfi0, slopethr = ls_gradient*dfi0;
		double b = alpha_max > 1 ? 1 : alpha_max;
		cg::fwdval result = (*func)(mat::add(x, h, 1, b, n), fpar);
		double fib = result.f;
		g = result.g;
		//printf("%e\n", fib);
		double dfib = mat::dot(g, h, n);
		double eps = pow(2, -52);
		int fiok = 0;
		double *xh = mat::create(n);
		if (ls_method == 1) {

		}
		else {
			int more = 1;
			double dfia = dfi0;
			while (more) {
				if (fib <= fi0 + slope0*b) {
					fiok = 1;
					alpha = b;
					fn = fib;
					mat::copy(g, gn, n);
					dfia = dfib;
					if (neval < ls_neval&&b < alpha_max&&dfib < slopethr) {
						if (2.5*b >= alpha_max) {
							b = alpha_max;
						}
						else {
							b *= 2;
						}
						mat::add(x, h, xh, 1, b, n);
						result = (*func)(xh, fpar);
						fib = result.f;
						free(g);
						g = result.g;
						neval++;
						dfib = mat::dot(g, h, n);
					}
					else {
						more = 0;
					}
				}
				else {
					fiok = 0;
					more = 0;
				}
			}
			if (!fiok) {
				double xfd11 = alpha;
				double xfd12 = fn;
				double xfd13 = dfia;
				double xfd21 = b;
				double xfd22 = fib;
				double xfd23 = dfib;
				double xfd31 = b;
				double xfd32 = fib;
				double xfd33 = dfib;
				while (neval < ls_neval&&xfd21 - xfd11>0 && (xfd33<slopethr || xfd32>fi0 + slope0*xfd31)) {
					double c;
					{
						double a = xfd11;
						double b = xfd21;
						double d = b - a;
						double dfia = xfd13;
						double C = xfd22 - xfd12 - d*dfia;
						if (C >= 5 * n*eps*b) {
							double A = a - 0.5*dfia*(d*d / C);
							d *= 0.1;
							if (a + d > A) {
								c = a + d;
							}
							else {
								c = A;
							}
							if (c > b - d) {
								c = b - d;
							}
						}
						else {
							c = (a + b) / 2;
						}
					}

					mat::add(x, h, xh, 1, c, n);
					cg::fwdval result = (*func)(xh, fpar);
					neval++;
					double fic = result.f;
					free(g);
					g = result.g;
					xfd31 = c;
					xfd32 = fic;
					xfd33 = mat::dot(g, h, n);

					if (fic < fi0 + slope0*c) {
						xfd11 = xfd31;
						xfd12 = xfd32;
						xfd13 = xfd33;
						alpha = c;
						fn = fic;
						mat::copy(g, gn, n);
					}
					else {
						xfd21 = xfd31;
						xfd22 = xfd32;
						xfd23 = xfd33;
					}
				}
			}
		}
		free(xh);
		free(g);
		return{ alpha,fn,gn,neval };
	}
	struct cg::val calc(cg::fwdval(*func)(double *, struct fwdpar), struct fwdpar par, double *x, int n,
		int cg_method, int ls_method, double upper_bound, double min_gradient, double min_step, int max_neval,
		double ls_gradient, double ls_step, int ls_neval) {
		x = mat::copy(x, n);
		int k = 1, kmax = max_neval;
		ls_method = 2;
		if (ls_gradient <= 0) {
			if (ls_method == 2) {
				ls_gradient = 1e-1;
			}
			else {
				ls_gradient = 1e-6;
			}
		}
		if (ls_step <= 0) {
			if (ls_method == 2) {
				ls_step = 1e-2;
			}
			else {
				ls_step = 1e-6;
			}
		}
		if (ls_neval <= 0) {
			ls_neval = 10;
		}
		cg::fwdval initval = (*func)(x, par);
		double f = initval.f;
		double *g = initval.g;
		double ng = mat::linf(g, n);
		double n2g = mat::norm(g, n);
		double gg = n2g*n2g;
		double gam = 0;
		double *h = mat::create(n, 0);
		double nh = 0;
		int neval = 1;
		printf("neval: %d\n", neval);
		int found = (ng <= min_gradient);
		double alph = 0;

		double *gpr = mat::create(n);
		double *sch = mat::create(n);
		while (!found) {
			double nhpr = nh;
			double ggpr = gg;
			mat::copy(g, gpr, n);
			k++;
			mat::add(h, g, h, gam, -1, n);
			nh = mat::norm(h, n);

			if (mat::dot(g, h, n) >= 1e-3*n2g*nh) {
				mat::copy(g, h, -1, n);
				nh = n2g;
			}
			double sc;
			if (k > 2) {
				sc = 0.9*alph*nhpr / nh;
			}
			else {
				sc = upper_bound / 32 / nh;
			}

			mat::copy(h, sch, sc, n);
			cg::lsval lval = cg::ls(func, par, x, f, g, sch, n,
				32, ls_method, ls_gradient, ls_step, kmax - neval > ls_neval ? ls_neval : kmax - neval);

			alph = sc*lval.alpha;
			f = lval.fn;
			g = lval.gn;
			int dval = lval.neval;
			ng = mat::linf(g, n);
			n2g = mat::norm(g, n);
			gg = n2g*n2g;

			if (cg_method == 1) {
				gam = gg / ggpr;
			}
			else {
				mat::add(g, gpr, gpr, 1, -1, n);
				gam = mat::dot(gpr, g, n) / ggpr;
			}

			mat::add(x, h, x, 1, alph, n);
			neval += dval;
			printf("neval: %d\n", neval);
			if (ng <= min_gradient) {
				found = 1;
			}
			else if (alph*nh <= min_step*(min_step*mat::norm(x, n))) {
				found = 2;
			}
			else if (neval >= kmax) {
				found = 3;
			}
			//printf("%e\n", alph);
		}

		free(g);
		free(h);
		free(gpr);
		free(sch);
		return{ x,f,ng,alph*nh ,k - 1,neval,found };
	}
}
