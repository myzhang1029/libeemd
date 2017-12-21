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

// Common error reporting and validation routines

#include "error.h"

libeemd_error_code validate_eemd_parameters(unsigned int ensemble_size, double noise_strength, unsigned int S_number, unsigned int num_siftings) {
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

// Helper functions for printing what error codes mean

void emd_error_string(libeemd_error_code err, char* err_string) {
	if (err == EMD_SUCCESS) {
		return;
	}
	switch (err) {
		case EMD_INVALID_ENSEMBLE_SIZE :
			strcpy(err_string, "Invalid ensemble size (zero or negative)");
			break;
		case EMD_INVALID_NOISE_STRENGTH :
			strcpy(err_string, "Invalid noise strength (negative)");
			break;
		case EMD_NOISE_ADDED_TO_EMD :
			strcpy(err_string, "Positive noise strength but ensemble size is one (regular EMD)");
			break;
		case EMD_NO_NOISE_ADDED_TO_EEMD :
			strcpy(err_string, "Ensemble size is more than one (EEMD) but noise strength is zero");
			break;
		case EMD_NO_CONVERGENCE_POSSIBLE :
			strcpy(err_string, "Stopping criteria invalid: would never converge");
			break;
		case EMD_NOT_ENOUGH_POINTS_FOR_SPLINE :
			strcpy(err_string, "Spline evaluation tried with insufficient points");
			break;
		case EMD_INVALID_SPLINE_POINTS :
			strcpy(err_string, "Spline evaluation points invalid");
			break;
		case EMD_GSL_ERROR :
			strcpy(err_string, "Error reported by GSL library");
			break;
		case EMD_NO_CONVERGENCE_IN_SIFTING :
			strcpy(err_string, "Convergence not reached after sifting 10000 times");
			break;
		default :
			strcpy(err_string, "Error code with unknown meaning. Please file a bug!");
			break;
	}
}

void emd_report_to_file_if_error(FILE* file, libeemd_error_code err) {
	if (err == EMD_SUCCESS) {
		return;
	}
	fputs("libeemd error: ", file);
	char* err_string = malloc(64 * sizeof(char));
	emd_error_string(err, err_string);
	fputs(err_string, file);
	free(err_string); err_string = NULL;
}

void emd_report_if_error(libeemd_error_code err) {
	emd_report_to_file_if_error(stderr, err);
}
