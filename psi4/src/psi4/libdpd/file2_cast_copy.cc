/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2019 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This file is part of Psi4.
 *
 * Psi4 is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, version 3.
 *
 * Psi4 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License along
 * with Psi4; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

/*! \file
    \ingroup DPD
    \brief Enter brief description of file here
*/


// This function is used for making copies of Fock mstrix and TEI in single-precision
// Used for mixed/single-precision calculation
// ZW 10/2019

#include <cstdio>
#include <cstring>
#include "dpd.h"

namespace psi {

int DPD::file2_cast_copy_dtof(dpdfile2<double> *InFile, int outfilenum, const char *label) {
    int h, i, j, row, col, my_irrep, rowtot, coltot;
    double ***matrix;
    float ***tmp_matrix;
    dpdfile2<float> OutFile;

    my_irrep = InFile->my_irrep;

    file2_init_sp(&OutFile, outfilenum, InFile->my_irrep, InFile->params->pnum, InFile->params->qnum, label);

    file2_mat_init(InFile);
    file2_mat_rd(InFile);
    file2_mat_init_sp(&OutFile);
    
    tmp_matrix = (float ***)malloc(OutFile.param->nirreps * sizeof(float **));

    for (h = 0; h < OutFile.params->nirreps; h++) {
        rowtot = OutFile.params->rowtot[h];
        coltot = OutFile.params->coltot[h ^ my_irrep]; 
        for (i=0; i < rowtot; i++){
            for (j=0; j < coltot; j++){
            tmp_matrix[h][i][j] = static_cast<float>(InFile->matrix[h][i][j]);
            }
        }
        if (rowtot && coltot)
            memcpy((void *)&(OutFile.matrix[h][0][0]), (const void *)&(tmp_matrix[h][0][0]),
                   sizeof(float) * rowtot * coltot);
    }
    
    file2_mat_wrt_sp(&OutFile);
    file2_mat_close_sp(&OutFile);
    file2_mat_close(InFile);
    file2_close_sp(&OutFile);

    free(tmp_matrix);

    return 0;
}

// This function is used for casting dpd files from single-precision to double-precision

int DPD::file2_cast_copy_ftod(dpdfile2<float> *InFile, int outfilenum, const char *label) {
    int h, i, j, row, col, my_irrep, rowtot, coltot;
    float ***matrix;
    double ***tmp_matrix;
    dpdfile2<double> OutFile;

    my_irrep = InFile->my_irrep;

    file2_init(&OutFile, outfilenum, InFile->my_irrep, InFile->params->pnum, InFile->params->qnum, label);

    file2_mat_init_sp(InFile);
    file2_mat_rd_sp(InFile);
    file2_mat_init(&OutFile);
    
    tmp_matrix = (double ***)malloc(OutFile.param->nirreps * sizeof(double **));

    for (h = 0; h < OutFile.params->nirreps; h++) {
        rowtot = OutFile.params->rowtot[h];
        coltot = OutFile.params->coltot[h ^ my_irrep]; 
        for (i=0; i < rowtot; i++){
            for (j=0; j < coltot; j++){
            tmp_matrix[h][i][j] = static_cast<double>(InFile->matrix[h][i][j]);
            }
        }
        if (rowtot && coltot)
            memcpy((void *)&(OutFile.matrix[h][0][0]), (const void *)&(tmp_matrix[h][0][0]),
                   sizeof(double) * rowtot * coltot);
    }
    
    file2_mat_wrt(&OutFile);
    file2_mat_close(&OutFile);
    file2_mat_close_sp(InFile);
    file2_close(&OutFile);

    free(tmp_matrix);

    return 0;
}


}  // namespace psi
