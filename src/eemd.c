/* Copyright 2013 Perttu Luukko

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

#include "eemd.h"

// If we are using OpenMP for parallel computation, we need locks to ensure
// that the same output data is not written by several threads at the same
// time.
#ifdef _OPENMP
typedef omp_lock_t lock;
inline static void init_lock(lock* l) { omp_init_lock(l); }
inline static void destroy_lock(lock* l) { omp_destroy_lock(l); }
inline static void get_lock(lock* l) { omp_set_lock(l); }
inline static void release_lock(lock* l) { omp_unset_lock(l); }
#else
// If we don't use OpenMP, we provide a dummy lock that does nothing. This
// avoids littering the code with too many #ifdefs for _OPENMP.
typedef char lock;
inline static void init_lock(__attribute__((unused)) lock* l) {}
inline static void destroy_lock(__attribute__((unused)) lock* l) {}
inline static void get_lock(__attribute__((unused)) lock* l) {}
inline static void release_lock(__attribute__((unused)) lock* l) {}
#endif


// Helper functions for working with data arrays
inline static void array_copy(double const* restrict src, size_t n, double* restrict dest) {
	memcpy(dest, src, n*sizeof(double));
}

inline static void array_add(double const* src, size_t n, double* dest) {
	for (size_t i=0; i<n; i++)
		dest[i] += src[i];
}

inline static void array_add_to(double const* src1, double const* src2, size_t n, double* dest) {
	for (size_t i=0; i<n; i++)
		dest[i] = src1[i] + src2[i];
}

inline static void array_addmul_to(double const* src1, double const* src2, double val, size_t n, double* dest) {
	for (size_t i=0; i<n; i++)
		dest[i] = src1[i] + val*src2[i];
}

inline static void array_sub(double const* src, size_t n, double* dest) {
	for (size_t i=0; i<n; i++)
		dest[i] -= src[i];
}

inline static void array_mult(double* dest, size_t n, double val) {
	for (size_t i=0; i<n; i++)
		dest[i] *= val;
}

// In the following part the necessary workspace memory structures for several
// EMD operations are defined

// For sifting we need arrays for storing the found extrema of the signal, and memory required
// to form the spline envelopes
typedef struct {
	// Number of samples in the signal
	size_t N;
	// Found extrema
	double* restrict maxx;
	double* restrict maxy;
	double* restrict minx;
	double* restrict miny;
	// Upper and lower envelope spline values
	double* restrict maxspline;
	double* restrict minspline;
	// Extra memory required for spline evaluation
	double* restrict spline_workspace;
} sifting_workspace;

sifting_workspace* allocate_sifting_workspace(size_t N) {
	sifting_workspace* w = malloc(sizeof(sifting_workspace));
	w->N = N;
	w->maxx = malloc(N*sizeof(double));
	w->maxy = malloc(N*sizeof(double));
	w->minx = malloc(N*sizeof(double));
	w->miny = malloc(N*sizeof(double));
	w->maxspline = malloc(N*sizeof(double));
	w->minspline = malloc(N*sizeof(double));
	// Spline evaluation requires 5*m-10 doubles where m is the number of
	// extrema. The worst case scenario is that every point is an extrema, so
	// use m=N to be safe.
	const size_t spline_workspace_size = (N > 2)? 5*N-10 : 0;
	w->spline_workspace = malloc(spline_workspace_size*sizeof(double));
	return w;
}

void free_sifting_workspace(sifting_workspace* w) {
	free(w->spline_workspace); w->spline_workspace = NULL;
	free(w->minspline); w->minspline = NULL;
	free(w->maxspline); w->maxspline = NULL;
	free(w->miny); w->miny = NULL;
	free(w->minx); w->minx = NULL;
	free(w->maxy); w->maxy = NULL;
	free(w->maxx); w->maxx = NULL;
	free(w); w = NULL;
}


// For EMD we need space to do the sifting and somewhere to save the residual from the previous run.
// We also leave room for an array of locks to protect multi-threaded EMD.
typedef struct {
	size_t N;
	// Previous residual for EMD
	double* restrict res;
	// What is needed for sifting
	sifting_workspace* restrict sift_w;
	// A pointer for shared locks. These locks are used to make EMD thread-safe
	// even when several threads run EMD with the same output matrix (we'll do
	// this in EEMD).
	lock** locks;
} emd_workspace;

emd_workspace* allocate_emd_workspace(size_t N) {
	emd_workspace* w = malloc(sizeof(emd_workspace));
	w->N = N;
	w->res = malloc(N*sizeof(double));
	w->sift_w = allocate_sifting_workspace(N);
	w->locks = NULL; // The locks are assumed to be allocated and freed independently
	return w;
}

void free_emd_workspace(emd_workspace* w) {
	free_sifting_workspace(w->sift_w);
	free(w->res); w->res = NULL;
	free(w); w = NULL;
}


// EEMD needs a random number generator in addition to emd_workspace. We also need a place to store
// the member of the ensemble (input signal + realization of noise) to be worked on.
typedef struct {
	size_t N;
	// The random number generator
	gsl_rng* r;
	// The ensemble member signal
	double* restrict x;
	// What is needed for running EMD
	emd_workspace* restrict emd_w;
} eemd_workspace;

eemd_workspace* allocate_eemd_workspace(size_t N) {
	eemd_workspace* w = malloc(sizeof(eemd_workspace));
	w->N = N;
	w->r = gsl_rng_alloc(gsl_rng_mt19937);
	w->x = malloc(N*sizeof(double));
	w->emd_w = allocate_emd_workspace(N);
	return w;
}

void set_rng_seed(eemd_workspace* w, unsigned long int rng_seed) {
	gsl_rng_set(w->r, rng_seed);
}

void free_eemd_workspace(eemd_workspace* w) {
	free_emd_workspace(w->emd_w);
	free(w->x); w->x = NULL;
	gsl_rng_free(w->r); w->r = NULL;
	free(w); w = NULL;
}

// Forward declaration of a helper function used internally for making a single
// EMD run with a preallocated workspace
static libeemd_error_code _emd(double* restrict input, emd_workspace* restrict w,
		double* restrict output, size_t M,
		unsigned int S_number, unsigned int num_siftings);

// Forward declaration of a helper function for applying the sifting procedure to
// input until it is reduced to an IMF according to the stopping criteria given
// by S_number and num_siftings
static libeemd_error_code _sift(double* restrict input, sifting_workspace*
		restrict w, unsigned int S_number, unsigned int num_siftings, unsigned int*
		sift_counter);

// Forward declaration of a helper function for parameter validation shared by functions eemd and ceemdan
static inline libeemd_error_code _validate_eemd_parameters(unsigned int ensemble_size, double noise_strength, unsigned int S_number, unsigned int num_siftings);

// Main EEMD decomposition routine definition
libeemd_error_code eemd(double const* restrict input, size_t N,
		double* restrict output, size_t M,
		unsigned int ensemble_size, double noise_strength, unsigned int
		S_number, unsigned int num_siftings, unsigned long int rng_seed) {
	gsl_set_error_handler_off();
	// Validate parameters
	libeemd_error_code validation_result = _validate_eemd_parameters(ensemble_size, noise_strength, S_number, num_siftings);
	if (validation_result != EMD_SUCCESS) {
		return validation_result;
	}
	// For empty data we have nothing to do
	if (N == 0) {
		return EMD_SUCCESS;
	}
	if (M == 0) {
		M = emd_num_imfs(N);
	}
	// The noise standard deviation is noise_strength times the standard deviation of input data
	const double noise_sigma = (noise_strength != 0)? gsl_stats_sd(input, 1, N)*noise_strength : 0;
	// Initialize output data to zero
	memset(output, 0x00, M*N*sizeof(double));
	// Each thread gets a separate workspace if we are using OpenMP
	eemd_workspace** ws = NULL;
	// The locks are shared among all threads
	lock** locks;
	// Don't start unnecessary threads if the ensemble is small
	#ifdef _OPENMP
	if (omp_get_num_threads() > (int)ensemble_size) {
		omp_set_num_threads(ensemble_size);
	}
	#endif
	unsigned int ensemble_counter = 0;
	// The following section is executed in parallel
	libeemd_error_code emd_err = EMD_SUCCESS;
	#pragma omp parallel
	{
		#ifdef _OPENMP
		const int num_threads = omp_get_num_threads();
		const int thread_id = omp_get_thread_num();
		#if EEMD_DEBUG >= 1
		#pragma omp single
		fprintf(stderr, "Using %d thread(s) with OpenMP.\n", num_threads);
		#endif
		#else
		const int num_threads = 1;
		const int thread_id = 0;
		#endif
		#pragma omp single
		{
			ws = malloc(num_threads*sizeof(eemd_workspace*));
			locks = malloc(M*sizeof(lock*));
			for (size_t i=0; i<M; i++) {
				locks[i] = malloc(sizeof(lock));
				init_lock(locks[i]);
			}
		}
		// Each thread allocates its own workspace
		ws[thread_id] = allocate_eemd_workspace(N);
		eemd_workspace* w = ws[thread_id];
		// All threads share the same array of locks
		w->emd_w->locks = locks;
		// Loop over all ensemble members, dividing them among the threads
		#pragma omp for
		for (size_t en_i=0; en_i<ensemble_size; en_i++) {
			// Check if an error has occured in other threads
			#pragma omp flush(emd_err)
			if (emd_err != EMD_SUCCESS) {
				continue;
			}
			// Initialize ensemble member as input data + noise
			if (noise_strength == 0.0) {
				array_copy(input, N, w->x);
			}
			else {
				// set rng seed based on ensemble member to ensure
				// reproducibility even in a multithreaded case
				set_rng_seed(w, rng_seed+en_i);
				for (size_t i=0; i<N; i++) {
					w->x[i] = input[i] + gsl_ran_gaussian(w->r, noise_sigma);
				}
			}
			// Extract IMFs with EMD
			emd_err = _emd(w->x, w->emd_w, output, M, S_number, num_siftings);
			#pragma omp flush(emd_err)
			#pragma omp atomic
			ensemble_counter++;
			#if EEMD_DEBUG >= 1
			fprintf(stderr, "Ensemble iteration %u/%u done.\n", ensemble_counter, ensemble_size);
			#endif
		}
		// Free resources
		free_eemd_workspace(w);
		#pragma omp single
		{
			free(ws); ws = NULL;
			for (size_t i=0; i<M; i++) {
				destroy_lock(locks[i]);
				free(locks[i]);
			}
			free(locks); locks = NULL;
		}
	} // End of parallel block
	if (emd_err != EMD_SUCCESS) {
		return emd_err;
	}
	// Divide output data by the ensemble size to get the average
	if (ensemble_size != 1) {
		const double one_per_ensemble_size = 1.0/ensemble_size;
		array_mult(output, N*M, one_per_ensemble_size);
	}
	return EMD_SUCCESS;
}

// Main CEEMDAN decomposition routine definition
libeemd_error_code ceemdan(double const* restrict input, size_t N,
		double* restrict output, size_t M,
		unsigned int ensemble_size, double noise_strength, unsigned int
		S_number, unsigned int num_siftings, unsigned long int rng_seed) {
	gsl_set_error_handler_off();
	// Validate parameters
	libeemd_error_code validation_result = _validate_eemd_parameters(ensemble_size, noise_strength, S_number, num_siftings);
	if (validation_result != EMD_SUCCESS) {
		return validation_result;
	}
	// For empty data we have nothing to do
	if (N == 0) {
		return EMD_SUCCESS;
	}
	// For M == 1 the only "IMF" is the residual
	if (M == 1) {
		memcpy(output, input, N*sizeof(double));
		return EMD_SUCCESS;
	}
	if (M == 0) {
		M = emd_num_imfs(N);
	}
	const double one_per_ensemble_size = 1.0/ensemble_size;
	// Initialize output data to zero
	memset(output, 0x00, M*N*sizeof(double));
	// Each thread gets a separate workspace if we are using OpenMP
	eemd_workspace** ws = NULL;
	// All threads need to write to the same row of the output matrix
	// so we need only one shared lock
	lock* output_lock = malloc(sizeof(lock));
	init_lock(output_lock);
	// The threads also share the same precomputed noise
	double* noises = malloc(ensemble_size*N*sizeof(double));
	// Since we need to decompose this noise by EMD, we also need arrays for storing
	// the residuals
	double* noise_residuals = malloc(ensemble_size*N*sizeof(double));
	// Don't start unnecessary threads if the ensemble is small
	#ifdef _OPENMP
	if (omp_get_num_threads() > (int)ensemble_size) {
		omp_set_num_threads(ensemble_size);
	}
	#endif
	int num_threads;
	// The following section is executed in parallel
	#pragma omp parallel
	{
		#ifdef _OPENMP
		num_threads = omp_get_num_threads();
		const int thread_id = omp_get_thread_num();
		#if EEMD_DEBUG >= 1
		#pragma omp single
		fprintf(stderr, "Using %d thread(s) with OpenMP.\n", num_threads);
		#endif
		#else
		num_threads = 1;
		const int thread_id = 0;
		#endif
		#pragma omp single
		{
			ws = malloc(num_threads*sizeof(eemd_workspace*));
		}
		// Each thread allocates its own workspace
		ws[thread_id] = allocate_eemd_workspace(N);
		eemd_workspace* w = ws[thread_id];
		// Precompute and store white noise, since for each mode of the data we
		// need the same mode of the corresponding realization of noise
		#pragma omp for
		for (size_t en_i=0; en_i<ensemble_size; en_i++) {
			// set rng seed based on ensemble member to ensure
			// reproducibility even in a multithreaded case
			set_rng_seed(w, rng_seed+en_i);
			for (size_t j=0; j<N; j++) {
				noises[N*en_i+j] = gsl_ran_gaussian(w->r, 1.0);
			}
		}
	} // Return to sequental mode
	// Allocate memory for the residual shared among all threads
	double* restrict res = malloc(N*sizeof(double));
	// For the first iteration the residual is the input signal
	array_copy(input, N, res);
	// Each mode is extracted sequentially, but we use parallelization in the inner loop
	// to loop over ensemble members
	for (size_t imf_i=0; imf_i<M; imf_i++) {
		// Provide a pointer to the output vector where this IMF will be stored
		double* const imf = &output[imf_i*N];
		// Then we go parallel to compute the different ensemble members
		libeemd_error_code sift_err = EMD_SUCCESS;
		#pragma omp parallel
		{
			#ifdef _OPENMP
			const int thread_id = omp_get_thread_num();
			#else
			const int thread_id = 0;
			#endif
			eemd_workspace* w = ws[thread_id];
			unsigned int sift_counter = 0;
			#pragma omp for
			for (size_t en_i=0; en_i<ensemble_size; en_i++) {
				// Check if an error has occured in other threads
				#pragma omp flush(sift_err)
				if (sift_err != EMD_SUCCESS) {
					continue;
				}
				// Provide a pointer to the noise vector and noise residual used by
				// this ensemble member
				double* const noise = &noises[N*en_i];
				double* const noise_residual = &noise_residuals[N*en_i];
				// Initialize input signal as data + noise.
				// The noise standard deviation is noise_strength times the
				// standard deviation of input data divided by the standard
				// deviation of the noise. This is used to fix the SNR at each
				// stage.
				const double noise_sd = gsl_stats_sd(noise, 1, N);
				const double noise_sigma = (noise_sd != 0)? noise_strength*gsl_stats_sd(res, 1, N)/noise_sd : 0;
				array_addmul_to(res, noise, noise_sigma, N, w->x);
				// Sift to extract first EMD mode
				sift_err = _sift(w->x, w->emd_w->sift_w, S_number, num_siftings, &sift_counter);
				#pragma omp flush(sift_err)
				// Sum to output vector
				get_lock(output_lock);
				array_add(w->x, N, imf);
				release_lock(output_lock);
				// Extract next EMD mode of the noise. This is used as the noise for
				// the next mode extracted from the data
				if (imf_i == 0) {
					array_copy(noise, N, noise_residual);
				}
				else {
					array_copy(noise_residual, N, noise);
				}
				sift_err = _sift(noise, w->emd_w->sift_w, S_number, num_siftings, &sift_counter);
				#pragma omp flush(sift_err)
				array_sub(noise, N, noise_residual);
			}
		} // Parallel section ends
		if (sift_err != EMD_SUCCESS) {
			return sift_err;
		}
		// Divide with ensemble size to get the average
		array_mult(imf, N, one_per_ensemble_size);
		// Subtract this IMF from the previous residual to form the new one
		array_sub(imf, N, res);
	}
	// Save final residual
	get_lock(output_lock);
	array_add(res, N, output+N*(M-1));
	release_lock(output_lock);
	// Free global resources
	for (int thread_id=0; thread_id<num_threads; thread_id++) {
		free_eemd_workspace(ws[thread_id]);
	}
	free(ws); ws = NULL;
	free(res); res = NULL;
	free(noise_residuals); noise_residuals = NULL;
	free(noises); noises = NULL;
	destroy_lock(output_lock);
	free(output_lock); output_lock = NULL;
	return EMD_SUCCESS;
}

static inline libeemd_error_code _validate_eemd_parameters(unsigned int ensemble_size, double noise_strength, unsigned int S_number, unsigned int num_siftings) {
	if (ensemble_size < 1) {
		return EMD_INVALID_ENSEMBLE_SIZE;
	}
	if (noise_strength < 0) {
		return EMD_INVALID_NOISE_STRENGTH;
	}
	if (ensemble_size == 1 && noise_strength > 0) {
		return EMD_NOISE_ADDED_TO_EMD;
	}
	if (ensemble_size > 1 && noise_strength == 0) {
		return EMD_NO_NOISE_ADDED_TO_EEMD;
	}
	if (S_number == 0 && num_siftings == 0) {
		return EMD_NO_CONVERGENCE_POSSIBLE;
	}
	return EMD_SUCCESS;
}

// Helper function for applying the sifting procedure to input until it is
// reduced to an IMF according to the stopping criteria given by S_number and
// num_siftings. The required number of siftings is saved to sift_counter.
static libeemd_error_code _sift(double* restrict input, sifting_workspace*
		restrict w, unsigned int S_number, unsigned int num_siftings,
		unsigned int* sift_counter) {
	const size_t N = w->N;
	// Provide some shorthands to avoid excessive '->' operators
	double* const maxx = w->maxx;
	double* const maxy = w->maxy;
	double* const minx = w->minx;
	double* const miny = w->miny;
	// Initialize counters that keep track of the number of siftings
	// and the S number
	*sift_counter = 0;
	unsigned int S_counter = 0;
	// Numbers of extrema and zero crossings are initialized to dummy values
	size_t num_max = (size_t)(-1);
	size_t num_min = (size_t)(-1);
	size_t num_zc = (size_t)(-1);
	size_t prev_num_max = (size_t)(-1);
	size_t prev_num_min = (size_t)(-1);
	size_t prev_num_zc = (size_t)(-1);
	while (num_siftings == 0 || *sift_counter < num_siftings) {
		(*sift_counter)++;
		#if EEMD_DEBUG >= 1
		if (*sift_counter == 10000) {
			fprintf(stderr, "Something is probably wrong. Sift counter has reached 10000.\n");
		}
		#endif
		prev_num_max = num_max;
		prev_num_min = num_min;
		prev_num_zc = num_zc;
		// Find extrema and count zero crossings
		emd_find_extrema(input, N, maxx, maxy, &num_max, minx, miny, &num_min, &num_zc);
		// Check if we are finished based on the S-number criteria
		if (S_number != 0) {
			const int max_diff = (int)num_max - (int)prev_num_max;
			const int min_diff = (int)num_min - (int)prev_num_min;
			const int zc_diff = (int)num_zc - (int)prev_num_zc;
			if (abs(max_diff)+abs(min_diff)+abs(zc_diff) <= 1) {
				S_counter++;
				if (S_counter >= S_number) {
					const int num_diff = (int)num_min + (int)num_max - 4 - (int)num_zc;
					if (abs(num_diff) <= 1) {
						// Number of extrema has been stable for S_number steps
						// and the number of *interior* extrema and zero
						// crossings differ by at most one -- we are converged
						// according to the S-number criterion
						break;
					}
				}
			}
			else {
				S_counter = 0;
			}
		}
		// Fit splines, choose order of spline based on the number of extrema
		libeemd_error_code max_errcode = emd_evaluate_spline(maxx, maxy, num_max, w->maxspline, w->spline_workspace);
		if (max_errcode != EMD_SUCCESS) {
			return max_errcode;
		}
		libeemd_error_code min_errcode = emd_evaluate_spline(minx, miny, num_min, w->minspline, w->spline_workspace);
		if (min_errcode != EMD_SUCCESS) {
			return min_errcode;
		}
		// Subtract envelope mean from the data
		for (size_t i=0; i<N; i++) {
			input[i] -= 0.5*(w->maxspline[i] + w->minspline[i]);
		}
	}
	return EMD_SUCCESS;
}

// Helper function for extracting all IMFs from input using the sifting
// procedure defined by _sift. The contents of the input array are destroyed in
// the process.
static libeemd_error_code _emd(double* restrict input, emd_workspace* restrict w,
		double* restrict output, size_t M,
		unsigned int S_number, unsigned int num_siftings) {
	// Provide some shorthands to avoid excessive '->' operators
	const size_t N = w->N;
	double* const res = w->res;
	lock** locks = w->locks;
	if (M == 0) {
		M = emd_num_imfs(N);
	}
	// We need to store a copy of the original signal so that once it is
	// reduced to an IMF we have something to subtract the IMF from to form
	// the residual for the next iteration
	array_copy(input, N, res);
	// Loop over all IMFs to be separated from input
	unsigned int sift_counter;
	for (size_t imf_i=0; imf_i<M-1; imf_i++) {
		if (imf_i != 0) {
			// Except for the first iteration, restore the previous residual
			// and use it as an input
			array_copy(res, N, input);
		}
		// Perform siftings on input until it is an IMF
		libeemd_error_code sift_err = _sift(input, w->sift_w, S_number, num_siftings, &sift_counter);
		if (sift_err != EMD_SUCCESS) {
			return sift_err;
		}
		// Subtract this IMF from the saved copy to form the residual for
		// the next round
		array_sub(input, N, res);
		// Add the discovered IMF to the output matrix. Use locks to ensure
		// other threads are not writing to the same row of the output matrix
		// at the same time
		get_lock(locks[imf_i]);
		array_add(input, N, output+N*imf_i);
		release_lock(locks[imf_i]);
		#if EEMD_DEBUG >= 2
		fprintf(stderr, "IMF %zd saved after %u siftings.\n", imf_i+1, sift_counter);
		#endif
	}
	// Save final residual
	get_lock(locks[M-1]);
	array_add(res, N, output+N*(M-1));
	release_lock(locks[M-1]);
	return EMD_SUCCESS;
}

size_t emd_num_imfs(size_t N) {
	if (N == 0) {
		return 0;
	}
	if (N <= 3) {
		return 1;
	}
	return (size_t)(log2(N));
}

// Helper functions for printing what error codes mean
void emd_report_to_file_if_error(FILE* file, libeemd_error_code err) {
	if (err == EMD_SUCCESS) {
		return;
	}
	fprintf(file, "libeemd error: ");
	switch (err) {
		case EMD_INVALID_ENSEMBLE_SIZE :
			fprintf(file, "Invalid ensemble size (zero or negative)\n");
			break;
		case EMD_INVALID_NOISE_STRENGTH :
			fprintf(file, "Invalid noise strength (negative)\n");
			break;
		case EMD_NOISE_ADDED_TO_EMD :
			fprintf(file, "Positive noise strength but ensemble size is one (regular EMD)\n");
			break;
		case EMD_NO_NOISE_ADDED_TO_EEMD :
			fprintf(file, "Ensemble size is more than one (EEMD) but noise strength is zero\n");
			break;
		case EMD_NO_CONVERGENCE_POSSIBLE :
			fprintf(file, "Stopping criteria invalid: would never converge\n");
			break;
		case EMD_NOT_ENOUGH_POINTS_FOR_SPLINE :
			fprintf(file, "Spline evaluation tried with insufficient points\n");
			break;
		case EMD_INVALID_SPLINE_POINTS :
			fprintf(file, "Spline evaluation points invalid\n");
			break;
		case EMD_GSL_ERROR :
			fprintf(file, "Error reported by GSL library\n");
			break;
		default :
			fprintf(file, "Error code with unknown meaning. Please file a bug!\n");
	}
}

void emd_report_if_error(libeemd_error_code err) {
	emd_report_to_file_if_error(stderr, err);
}
