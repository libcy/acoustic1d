#pragma once

#include<stdio.h>
#include<stdlib.h>
#include<math.h>

namespace mat {
	double *create(int n) {
		return (double *)malloc(n*sizeof(double));
	}
	int *createInt(int n) {
		return (int *)malloc(n*sizeof(int));
	}
	double *create(int n, double val) {
		double *a = (double *)malloc(n*sizeof(double));
		for (int i = 0; i < n; i++) {
			a[i] = val;
		}
		return a;
	}
	void add(double *a, double *b, double *c, int n) {
		for (int i = 0; i < n; i++) {
			c[i] = a[i] + b[i];
		}
	}
	void add(double *a, double *b, double *c, double ca, double cb, int n) {
		for (int i = 0; i < n; i++) {
			c[i] = ca*a[i] + cb*b[i];
		}
	}
	double *add(double *a, double *b, int n) {
		double *c = mat::create(n);
		for (int i = 0; i < n; i++) {
			c[i] = a[i] + b[i];
		}
		return c;
	}
	double *add(double *a, double *b, double ca, double cb, int n) {
		double *c = mat::create(n);
		for (int i = 0; i < n; i++) {
			c[i] = ca*a[i] + cb*b[i];
		}
		return c;
	}
	void copy(double *a, double *b, int n) {
		for (int i = 0; i < n; i++) {
			b[i] = a[i];
		}
	}
	void copy(double *a, double *b, double ca, int n) {
		for (int i = 0; i < n; i++) {
			b[i] = ca*a[i];
		}
	}
	double *copy(double *a, int n) {
		double *b = mat::create(n);
		for (int i = 0; i < n; i++) {
			b[i] = a[i];
		}
		return b;
	}
	double *copy(double *a, double ca, int n) {
		double *b = mat::create(n);
		for (int i = 0; i < n; i++) {
			b[i] = ca*a[i];
		}
		return b;
	}
	double *transpose(double *a, int m, int n) {
		double *b = mat::create(m* n);
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				b[i*n + j] = a[j*m + i];
			}
		}
		return b;
	}
	double *multiply(double *a, double *b, int m, int n, int p) {
		double *c = mat::create(m*n, 0);
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				int ij = i * n + j;
				for (int k = 0; k < p; k++) {
					c[ij] += a[i*p + k] * b[k*n + j];
				}
			}
		}
		return c;
	}
	double dot(double *a, double *b, int n) {
		double ab = 0;
		for (int i = 0; i < n; i++) {
			ab += a[i] * b[i];
		}
		return ab;
	}
	double norm(double *a, int n) {
		double na = 0;
		for (int i = 0; i < n; i++) {
			na += a[i] * a[i];
		}
		return sqrt(na);
	}
	double linf(double *a, int n) {
		double na = 0;
		for (int i = 0; i < n; i++) {
			if (a[i] < 0) {
				if (-a[i] > na) {
					na = -a[i];
				}
			}
			else if (a[i] > na) {
				na = a[i];
			}
		}
		return na;
	}
	double max(double *a, int n) {
		double amax = a[0];
		for (int i = 1; i < n; i++) {
			if (a[i] > amax) {
				amax = a[i];
			}
		}
		return amax;
	}
	double min(double *a, int n) {
		double amin = a[0];
		for (int i = 1; i < n; i++) {
			if (a[i] < amin) {
				amin = a[i];
			}
		}
		return amin;
	}
	void plot(double *y, int n, int filename){
		double ymin = mat::min(y, n);
		double ymax = mat::max(y, n);
		double cy;
		if (ymax == ymin){
			ymin = 0;
			cy = 0.5 / ymax;
		}
		else{
			cy = 1 / (ymax - ymin);
		}

		FILE *wave = fopen("wave", "w");
		for (int i = 0; i < n; i++) {
			fprintf(wave, "%f %f\n", (double)2 * i / (n - 1), (y[i] - ymin) * cy);
		}
		fclose(wave);

		char plotwave[150];
		sprintf(plotwave, "gmt psxy wave -R0/2/0/1 -Jx3i -P -Ba0.2f0.02:\"%.1e ~ %.1e\":/a0.2f0.02:\"\":SWne -Y19 -Wblue > wave%d.ps", ymin, ymax,filename);
		system(plotwave);
	}
	void plot(double *y, int n) {
		mat::plot(y, n, 0);
		system("start wave0.ps");
	}
	void plot(double **y, int m, int n, int filename) {
		double ymin = y[0][0];
		double ymax = y[0][0];
		FILE *wave = fopen("wave", "w");
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < m; j++) {
				if (y[j][i] < ymin) {
					ymin = y[j][i];
				}
				if (y[j][i] > ymax) {
					ymax = y[j][i];
				}
			}
		}
		double cy;
		if (ymax == ymin){
			ymin = 0;
			cy = 0.5 / ymax;
		}
		else{
			cy = 1 / (ymax - ymin);
		}


		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				fprintf(wave, "%f %f\n", (double)2 * j / (n - 1), (y[i][j] - ymin) * cy);
			}
			fprintf(wave, ">\n");
		}
		fclose(wave);

		char plotwave[150];
		sprintf(plotwave, "gmt psxy wave -R0/2/0/1 -Jx3i -P -Ba0.2f0.02:\"%.1e ~ %.1e\":/a0.2f0.02:\"\":SWne -Y19 -Wblue > wave%d.ps", ymin, ymax, filename);
		system(plotwave);
	}
	void plot(double **y, int m, int n){
		plot(y, m, n, 0);
		system("start wave0.ps");
	}
}