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
#include "psi4/libqt/qt.h"
#include "dpd.h"

namespace psi {

double DPD::buf4_dot_sp(dpdbuf4<float> *BufA, dpdbuf4<float> *BufB) {
    int h, nirreps, n, my_irrep;
    float dot;
    int incore, nbuckets;
    long int memoryd, rows_per_bucket, rows_left;

    nirreps = BufA->params->nirreps;
    my_irrep = BufA->file.my_irrep;

    dot = 0.0;

    for (h = 0; h < nirreps; h++) {
        memoryd = dpd_memfree();

        if (BufA->params->rowtot[h] && BufA->params->coltot[h ^ my_irrep]) {
            /* Compute the memory for one row of A/B */
            if (BufA->params->coltot[h ^ my_irrep]) /* NB: we need at least one row of both A and B */
                rows_per_bucket = memoryd / (2 * BufA->params->coltot[h ^ my_irrep]);
            else
                rows_per_bucket = -1;

            if (rows_per_bucket > BufA->params->rowtot[h]) rows_per_bucket = BufA->params->rowtot[h];

            if (!rows_per_bucket) dpd_error("buf4_dot: Not enough memory for one row!", "outfile");

            nbuckets = (int)ceil((double)BufA->params->rowtot[h] / (double)rows_per_bucket);

            rows_left = BufA->params->rowtot[h] % rows_per_bucket;

            incore = 1;
            if (nbuckets > 1) incore = 0;

        } else
            incore = 1;

        if (incore) {
            buf4_mat_irrep_init_sp(BufA, h);
            buf4_mat_irrep_init_sp(BufB, h);
            buf4_mat_irrep_rd_sp(BufA, h);
            buf4_mat_irrep_rd_sp(BufB, h);

            dot += dot_block_sp(BufA->matrix[h], BufB->matrix[h], BufA->params->rowtot[h],
                             BufA->params->coltot[h ^ my_irrep], 1.0);

            buf4_mat_irrep_close_sp(BufA, h);
            buf4_mat_irrep_close_sp(BufB, h);
        } else {
            buf4_mat_irrep_init_block_sp(BufA, h, rows_per_bucket);
            buf4_mat_irrep_init_block_sp(BufB, h, rows_per_bucket);

            for (n = 0; n < (rows_left ? nbuckets - 1 : nbuckets); n++) {
                buf4_mat_irrep_rd_block_sp(BufA, h, n * rows_per_bucket, rows_per_bucket);
                buf4_mat_irrep_rd_block_sp(BufB, h, n * rows_per_bucket, rows_per_bucket);

                dot += dot_block_sp(BufA->matrix[h], BufB->matrix[h], rows_per_bucket, BufA->params->coltot[h ^ my_irrep],
                                 1.0);
            }

            if (rows_left) {
                buf4_mat_irrep_rd_block_sp(BufA, h, n * rows_per_bucket, rows_left);
                buf4_mat_irrep_rd_block_sp(BufB, h, n * rows_per_bucket, rows_left);

                dot += dot_block_sp(BufA->matrix[h], BufB->matrix[h], rows_left, BufA->params->coltot[h ^ my_irrep], 1.0);
            }

            buf4_mat_irrep_close_block_sp(BufA, h, rows_per_bucket);
            buf4_mat_irrep_close_block_sp(BufB, h, rows_per_bucket);
        }
    }

    return dot;
}

}  // namespace psi
