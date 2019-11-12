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
#include <cmath>
#include "dpd.h"
#include "psi4/libqt/qt.h"

namespace psi {

/* dpd_buf4_axpy(): Evaluates the standard operation a * X + Y -> Y for dpd
** four-index buffers.
**
** Arguments:
**   dpdbuf4<double>*BufX: A pointer to the leftmost dpdbuf4.
**   dpdbuf4<double>*BufY: A pointer to the rightmost (and target)
**                        dpdbuf4.
**   double alpha: The scalar prefactor in the multiplication.
*/
int DPD::buf4_axpy_sp(dpdbuf4<float> *BufX, dpdbuf4<float> *BufY, float alpha) {
    int h, nirreps, my_irrep;
    int row, col, incore, n, nbuckets;
    long int length;
    long int memoryd, rows_per_bucket, rows_left;
    float *X, *Y;

    nirreps = BufX->params->nirreps;
    my_irrep = BufX->file.my_irrep;

#ifdef DPD_TIMER
    timer_on("buf4_axpy");
#endif

    for (h = 0; h < nirreps; h++) {
        memoryd = (dpd_memfree() - BufX->file.params->coltot[h ^ my_irrep]) / 2; /* use half the memory for each buf4 */
        if (BufX->params->rowtot[h] && BufX->params->coltot[h ^ my_irrep]) {
            rows_per_bucket = memoryd / BufX->params->coltot[h ^ my_irrep];

            /* enough memory for the whole matrix? */
            if (rows_per_bucket > BufX->params->rowtot[h]) rows_per_bucket = BufX->params->rowtot[h];

            if (!rows_per_bucket) dpd_error("buf4_axpy: Not enough memory for one row!", "outfile");

            nbuckets = (int)ceil(((double)BufX->params->rowtot[h]) / ((double)rows_per_bucket));

            rows_left = BufX->params->rowtot[h] % rows_per_bucket;

            incore = 1;
            if (nbuckets > 1) {
                incore = 0;
#if DPD_DEBUG
                outfile->Printf("buf4_axpy: memory information.\n");
                outfile->Printf("buf4_axpy: rowtot[%d] = %d\n", h, BufX->params->rowtot[h]);
                outfile->Printf("buf4_axpy: nbuckets = %d\n", nbuckets);
                outfile->Printf("buf4_axpy: rows_per_bucket = %d\n", rows_per_bucket);
                outfile->Printf("buf4_axpy: rows_left = %d\n", rows_left);
                outfile->Printf("buf4_axpy: out-of-core algorithm used\n");
#endif
            }
        } else
            incore = 1;

        if (incore) {
            buf4_mat_irrep_init_sp(BufX, h);
            buf4_mat_irrep_rd_sp(BufX, h);

            buf4_mat_irrep_init_sp(BufY, h);
            buf4_mat_irrep_rd_sp(BufY, h);

            length = ((long)BufX->params->rowtot[h]) * ((long)BufX->params->coltot[h ^ my_irrep]);
            if (length) {
                X = &(BufX->matrix[h][0][0]);
                Y = &(BufY->matrix[h][0][0]);
                C_SAXPY(length, alpha, X, 1, Y, 1);
            }

            buf4_mat_irrep_wrt_sp(BufY, h);

            buf4_mat_irrep_close_sp(BufX, h);
            buf4_mat_irrep_close_sp(BufY, h);
        } else {
            buf4_mat_irrep_init_block_sp(BufX, h, rows_per_bucket);
            buf4_mat_irrep_init_block_sp(BufY, h, rows_per_bucket);

            length = ((long)rows_per_bucket) * ((long)BufX->params->coltot[h ^ my_irrep]);
            X = &(BufX->matrix[h][0][0]);
            Y = &(BufY->matrix[h][0][0]);

            for (n = 0; n < (rows_left ? nbuckets - 1 : nbuckets); n++) {
                buf4_mat_irrep_rd_block_sp(BufX, h, n * rows_per_bucket, rows_per_bucket);
                buf4_mat_irrep_rd_block_sp(BufY, h, n * rows_per_bucket, rows_per_bucket);

                C_SAXPY(length, alpha, X, 1, Y, 1);

                buf4_mat_irrep_wrt_block_sp(BufY, h, n * rows_per_bucket, rows_per_bucket);
            }

            if (rows_left) {
                length = ((long)rows_left) * ((long)BufX->params->coltot[h ^ my_irrep]);

                buf4_mat_irrep_rd_block_sp(BufX, h, n * rows_per_bucket, rows_left);
                buf4_mat_irrep_rd_block_sp(BufY, h, n * rows_per_bucket, rows_left);

                C_SAXPY(length, alpha, X, 1, Y, 1);

                buf4_mat_irrep_wrt_block_sp(BufY, h, n * rows_per_bucket, rows_left);
            }

            buf4_mat_irrep_close_block_sp(BufX, h, rows_per_bucket);
            buf4_mat_irrep_close_block_sp(BufY, h, rows_per_bucket);
        }
    }

#ifdef DPD_TIMER
    timer_off("buf4_axpy");
#endif

    return 0;
}

}  // namespace psi
