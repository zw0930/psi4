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
#include <cstdio>
#include "psi4/libqt/qt.h"
#include "dpd.h"

namespace psi {

/* file2_axpy(): Evaluates the standard operation a * X + Y -> Y for
 ** dpdfile2's.
 **
 ** Arguments:
 **   dpdfile2<double>*FileA: A pointer to the leftmost dpdfile2.
 **   dpdfile2<double>*FileB: A pointer to the rightmost (and target) dpdfile2.
 **   double alpha: The scalar prefactor in the multiplication.
 **   int transA: A boolean indicating that we should use the transpose of
 **               FileA
 */

int DPD::file2_axpy_sp(dpdfile2<float> *FileA, dpdfile2<float> *FileB, float alpha, int transA) {
    int h, nirreps, my_irrep;
    int row, col;

    nirreps = FileA->params->nirreps;
    my_irrep = FileA->my_irrep;

    file2_mat_init_sp(FileA);
    file2_mat_init_sp(FileB);
    file2_mat_rd_sp(FileA);
    file2_mat_rd_sp(FileB);

    for (h = 0; h < nirreps; h++) {
        if (!transA) {
            for (row = 0; row < FileA->params->rowtot[h]; row++)
                for (col = 0; col < FileA->params->coltot[h ^ my_irrep]; col++)
                    FileB->matrix[h][row][col] += alpha * FileA->matrix[h][row][col];

        } else {
            for (row = 0; row < FileB->params->rowtot[h]; row++)
                for (col = 0; col < FileB->params->coltot[h ^ my_irrep]; col++)
                    FileB->matrix[h][row][col] += alpha * FileA->matrix[h ^ my_irrep][col][row];
        }
    }

    file2_mat_wrt_sp(FileB);
    file2_mat_close_sp(FileA);
    file2_mat_close_sp(FileB);

    return 0;
}

}  // namespace psi
