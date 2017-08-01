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

#include "memd.h"

memd_sifting_workspace* allocate_memd_sifting_workspace(size_t N, lock* output_lock) {
	memd_sifting_workspace* w = malloc(sizeof(memd_sifting_workspace));
	w->N = N;
	w->projected_signal = malloc(N*sizeof(double));
	w->maxx = malloc(N*sizeof(double));
	w->maxy = malloc(N*sizeof(double));
	w->maxspline = malloc(N*sizeof(double));
	// Spline evaluation requires 5*m-10 doubles where m is the number of
	// extrema. The worst case scenario is that every point is an extrema, so
	// use m=N to be safe.
	const size_t spline_workspace_size = (N > 2)? 5*N-10 : 0;
	w->spline_workspace = malloc(spline_workspace_size*sizeof(double));
	w->output_lock = output_lock;
	return w;
}

void free_memd_sifting_workspace(memd_sifting_workspace* w) {
	free(w->projected_signal); w->projected_signal = NULL;
	free(w->maxx); w->maxx = NULL;
	free(w->maxy); w->maxy = NULL;
	free(w->maxspline); w->maxspline = NULL;
	free(w->spline_workspace); w->spline_workspace = NULL;
	free(w); w = NULL;
}


static libeemd_error_code _memd_sift_once(double complex* restrict x, size_t D, size_t N, double const* restrict directions, size_t num_directions, memd_sifting_workspace* w) {
	libeemd_error_code errcode = EMD_SUCCESS;
	double complex* restrict m = calloc(N, sizeof(double complex));
	double* const px = w->projected_signal;
	// TODO: handle different directions in parallel
	for (size_t direction_i=0; direction_i<num_directions; direction_i++) {
		const double phi = directions[direction_i];
		const double sin_phi = sin(phi);
		const double cos_phi = cos(phi);
		// Project signal
		for (size_t i=0; i<N; i++) {
			const double a = creal(x[i]);
			const double b = cimag(x[i]);
			px[i] = a*cos_phi + b*sin_phi;
		}
		// Find maxima
		emd_find_maxima(px, N, w->maxx, w->maxy, &(w->num_max));
		// Fit spline
		errcode = emd_evaluate_spline(w->maxx, w->maxy, w->num_max, w->maxspline, w->spline_workspace);
		if (errcode != EMD_SUCCESS) {
			return errcode;
		}
		// Add to m
		for (size_t i=0; i<N; i++) {
			m[i] += cexp(phi*I) * (w->maxspline)[i];
		}
	}
	// Scale m
	complex_array_mult(m, N, 2.0/(double)num_directions);
	// Subtract mean from input
	complex_array_sub(m, N, x);
	// Done
	free(m); m = NULL;
	return errcode;
}

libeemd_error_code memd(double const* restrict input, size_t D, size_t N,
		double const* restrict directions, size_t num_directions,
		double* restrict output, size_t M,
		unsigned int num_siftings) {
	gsl_set_error_handler_off();
	if (M == 0) {
		M = emd_num_imfs(N);
	}
	libeemd_error_code memd_err = EMD_SUCCESS;
	// Create a read-write copy of input data
	double* const x = malloc(N*Dsizeof(double));
	array_copy(input, N*D, x);
	double* const res = malloc(N*D*sizeof(double));
	// For the first iteration, the residual is the original input data
	array_copy(input, N*D, res);
	memd_sifting_workspace* w = allocate_memd_sifting_workspace(N, NULL);
	// Loop over all IMFs to be separated from input
	for (size_t imf_i=0; imf_i<M-1; imf_i++) {
		if (imf_i != 0) {
			// Except for the first iteration, restore the previous residual
			// and use it as an input
			array_copy(res, N*D, x);
		}
		// Perform siftings on x until it is an IMF
		for (unsigned int sift_counter=0; sift_counter<num_siftings; sift_counter++) {
			memd_err = _memd_sift_once(x, D, N, directions, num_directions, w);
			if (memd_err != EMD_SUCCESS) {
				return memd_err;
			}
		}
		// Subtract this IMF from the saved copy to form the residual for
		// the next round
		array_sub(x, N*D, res);
		// Write the discovered IMF to the output matrix
		array_copy(x, N*D, output+N*imf_i);
	}
	// Save final residual
	array_copy(res, N*D, output+N*(M-1));
	free_memd_sifting_workspace(w);
	free(res);
	free(x);
	return memd_err;
}
