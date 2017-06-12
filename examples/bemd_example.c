/* Copyright 2017 Perttu Luukko

 * This file is part of libeemd.

 * libeemd is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.

 * libeemd is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.

 * You should have received a copy of the GNU General Public License
 * along with libeemd.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include <gsl/gsl_math.h>
const double pi = M_PI;

#include "eemd.h"

const unsigned int num_directions = 64;
const unsigned int num_siftings = 10;
const size_t num_imfs = 4;
const char outfile[] = "bemd_example.out";

const size_t N = 512;

int main(void) {
	libeemd_error_code err;
	// Decompose a synthetic signal made out of rotating components
	double complex* inp = malloc(N*sizeof(double complex));
	for (size_t i=0; i<N; i++) {
		const double t = 2*pi*i/(double)N;
		inp[i] = cos(0.3*t)*cexp(2*I*t) + 0.3*fabs(sin(2.3*t))*cexp(17*I*t);
	}
	// Use evenly spaced angles as directions
	double* directions = malloc(num_directions*sizeof(double));
	for (size_t d=0; d<num_directions; d++) {
		directions[d] = 2*pi*(d+1)/(double)num_directions;
	}
	// Allocate memory for output data
	double complex* outp = malloc(num_imfs*N*sizeof(double complex));
	// Run BEMD
	err = bemd(inp, N, directions, num_directions, outp, num_imfs, num_siftings);
	if (err != EMD_SUCCESS) {
		emd_report_if_error(err);
		exit(1);
	}
	// Write output to file
	FILE* fp = fopen(outfile, "w");
	for (size_t j=0; j<N; j++) {
		const double complex z = inp[j];
		fprintf(fp, "%f%+fj ", creal(z), cimag(z));
	}
	fprintf(fp, "\n");
	for (size_t i=0; i<num_imfs; i++) {
		for (size_t j=0; j<N; j++) {
			const double complex z = outp[i*N+j];
			fprintf(fp, "%f%+fj ", creal(z), cimag(z));
		}
		fprintf(fp, "\n");
	}
	printf("Done! Run bemd_example_plot.py to visualize the results.\n");
	// Cleanup
	fclose(fp);
	free(inp); inp = NULL;
	free(directions); directions = NULL;
	free(outp); outp = NULL;
}
