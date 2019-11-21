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
#include <cstring>
#include <cmath>
#include "dpd.h"

namespace psi {

/* dpd_buf4_cast_copy(): Casts the an existing four-index dpdbuf4<double>and copies it into another file.
**
** Arguments:
**   dpdbuf4<double>*InBuf: A pointer to the given dpd buffer.
**   int outfilenum: The PSI unit number for the new buffer.
**   char *label: A string labelling for this buffer.
**
** NB: The buffer and file pq/rs parameters are assumed to be
** identical for the copy, obviously.  Hence, anti flag must be off.
**
** Converted to use buf4 only rather than assumptions about file4.
** TDC, September 1999
*/

int DPD::buf4_cast_copy_dtof(dpdbuf4<double> *InBuf, int outfilenum, const char *label) {
    int h, i, j, row, col, my_irrep;
    long int rowtot, coltot;
    int nbuckets, incore, n;
    long int memoryd, rows_per_bucket, rows_left, size;
    float ***tmp_matrix;
    dpdbuf4<float> OutBuf;
    
    my_irrep = InBuf->file.my_irrep;

    buf4_init_sp(&OutBuf, outfilenum, InBuf->file.my_irrep, InBuf->params->pqnum, InBuf->params->rsnum,
              InBuf->params->pqnum, InBuf->params->rsnum, 0, label);
    
    tmp_matrix = (float ***)malloc(OutBuf.params->nirreps * sizeof(float **));

    for (h = 0; h < InBuf->params->nirreps; h++) {
        memoryd = dpd_memfree() / 2; /* use half the memory for each buf4 */

        rowtot = InBuf->params->rowtot[h];
        coltot = InBuf->params->coltot[h ^ my_irrep];

        if (rowtot && coltot) {
            rows_per_bucket = memoryd / coltot;
            /* enough memory for the whole matrix? */
            if (rows_per_bucket > rowtot) rows_per_bucket = rowtot;

            if (!rows_per_bucket) dpd_error("buf4_scmcopy: Not enough memory for one row!", "outfile");

            nbuckets = (int)ceil(((double)rowtot) / ((double)rows_per_bucket));

            rows_left = rowtot % rows_per_bucket;

            incore = 1;
            if (nbuckets > 1) {
                incore = 0;
#if DPD_DEBUG
                outfile->Printf("buf4_copy: memory information.\n");
                outfile->Printf("buf4_copy: rowtot[%d] = %d.\n", h, InBuf->params->rowtot[h]);
                outfile->Printf("buf4_copy: nbuckets = %d\n", nbuckets);
                outfile->Printf("buf4_copy: rows_per_bucket = %d\n", rows_per_bucket);
                outfile->Printf("buf4_copy: rows_left = %d\n", rows_left);
                outfile->Printf("buf4_copy: out-of-core algorithm used\n");
#endif
            }

            if (incore) {
                buf4_mat_irrep_init(InBuf, h);
                buf4_mat_irrep_rd(InBuf, h);

                buf4_mat_irrep_init_sp(&OutBuf, h);

                if (rowtot && coltot){
                    for (i=0; i < rowtot; i++){
                        for (j=0; j< coltot; j++){
                            tmp_matrix[h][i][j] = static_cast<float>(InBuf->matrix[h][i][j]);
                        }
                    }
                    memcpy((void *)&(OutBuf.matrix[h][0][0]), (const void *)&(tmp_matrix[h][0][0]),
                           sizeof(float) * rowtot * coltot);
                }
                buf4_mat_irrep_wrt_sp(&OutBuf, h);

                buf4_mat_irrep_close_sp(&OutBuf, h);
                buf4_mat_irrep_close(InBuf, h);
            } else {
                buf4_mat_irrep_init_block(InBuf, h, rows_per_bucket);
                buf4_mat_irrep_init_block_sp(&OutBuf, h, rows_per_bucket);
                
                coltot = InBuf->params->coltot[h ^ my_irrep];
                size = ((long)rows_per_bucket) * ((long)coltot);
                
                for (n = 0; n < (rows_left ? nbuckets - 1 : nbuckets); n++) {
                    buf4_mat_irrep_rd_block(InBuf, h, n * rows_per_bucket, rows_per_bucket);
                    for (i=0; i < rows_per_bucket; i++){
                        for (j=0; j< coltot; j++){
                            tmp_matrix[h][i][j] = static_cast<float>(InBuf->matrix[h][i][j]);
                        }
                    }
                    memcpy((void *)&(OutBuf.matrix[h][0][0]), (const void *)&(tmp_matrix[h][0][0]),
                           ((long)sizeof(float)) * size);

                    buf4_mat_irrep_wrt_block_sp(&OutBuf, h, n * rows_per_bucket, rows_per_bucket);
                }
                if (rows_left) {
                    size = ((long)rows_left) * ((long)coltot);

                    buf4_mat_irrep_rd_block(InBuf, h, n * rows_per_bucket, rows_left);
                    for (i=0; i < rows_left; i++){
                        for (j=0; j< coltot; j++){
                            tmp_matrix[h][i][j] = static_cast<float>(InBuf->matrix[h][i][j]);
                        }
                    }
                    memcpy((void *)&(OutBuf.matrix[h][0][0]), (const void *)&(tmp_matrix[h][0][0]),
                           ((long)sizeof(float)) * size);

                    buf4_mat_irrep_wrt_block_sp(&OutBuf, h, n * rows_per_bucket, rows_left);
                }

                buf4_mat_irrep_close_block(InBuf, h, rows_per_bucket);
                buf4_mat_irrep_close_block_sp(&OutBuf, h, rows_per_bucket);
            }
        }
    }
    free(tmp_matrix);
    buf4_close_sp(&OutBuf);

    return 0;
}

int DPD::buf4_cast_copy_ftod(dpdbuf4<float> *InBuf, int outfilenum, const char *label) {
    int h, i, j, row, col, my_irrep;
    long int rowtot, coltot;
    int nbuckets, incore, n;
    long int memoryd, rows_per_bucket, rows_left, size;
    double ***tmp_matrix;
    dpdbuf4<double> OutBuf;
    
    my_irrep = InBuf->file.my_irrep;

    buf4_init(&OutBuf, outfilenum, InBuf->file.my_irrep, InBuf->params->pqnum, InBuf->params->rsnum,
              InBuf->params->pqnum, InBuf->params->rsnum, 0, label);
    
    tmp_matrix = (double ***)malloc(OutBuf.params->nirreps * sizeof(double **));

    for (h = 0; h < InBuf->params->nirreps; h++) {
        memoryd = dpd_memfree() / 2; /* use half the memory for each buf4 */

        rowtot = InBuf->params->rowtot[h];
        coltot = InBuf->params->coltot[h ^ my_irrep];

        if (rowtot && coltot) {
            rows_per_bucket = memoryd / coltot;
            /* enough memory for the whole matrix? */
            if (rows_per_bucket > rowtot) rows_per_bucket = rowtot;

            if (!rows_per_bucket) dpd_error("buf4_scmcopy: Not enough memory for one row!", "outfile");

            nbuckets = (int)ceil(((double)rowtot) / ((double)rows_per_bucket));

            rows_left = rowtot % rows_per_bucket;

            incore = 1;
            if (nbuckets > 1) {
                incore = 0;
#if DPD_DEBUG
                outfile->Printf("buf4_copy: memory information.\n");
                outfile->Printf("buf4_copy: rowtot[%d] = %d.\n", h, InBuf->params->rowtot[h]);
                outfile->Printf("buf4_copy: nbuckets = %d\n", nbuckets);
                outfile->Printf("buf4_copy: rows_per_bucket = %d\n", rows_per_bucket);
                outfile->Printf("buf4_copy: rows_left = %d\n", rows_left);
                outfile->Printf("buf4_copy: out-of-core algorithm used\n");
#endif
            }

            if (incore) {
                buf4_mat_irrep_init_sp(InBuf, h);
                buf4_mat_irrep_rd_sp(InBuf, h);

                buf4_mat_irrep_init(&OutBuf, h);

                if (rowtot && coltot){
                    for (i=0; i < rowtot; i++){
                        for (j=0; j< coltot; j++){
                            tmp_matrix[h][i][j] = static_cast<double>(InBuf->matrix[h][i][j]);
                        }
                    }
                    memcpy((void *)&(OutBuf.matrix[h][0][0]), (const void *)&(tmp_matrix[h][0][0]),
                           sizeof(double) * rowtot * coltot);
                }
                buf4_mat_irrep_wrt(&OutBuf, h);

                buf4_mat_irrep_close(&OutBuf, h);
                buf4_mat_irrep_close_sp(InBuf, h);
            } else {
                buf4_mat_irrep_init_block_sp(InBuf, h, rows_per_bucket);
                buf4_mat_irrep_init_block(&OutBuf, h, rows_per_bucket);
                
                coltot = InBuf->params->coltot[h ^ my_irrep];
                size = ((long)rows_per_bucket) * ((long)coltot);
                
                for (n = 0; n < (rows_left ? nbuckets - 1 : nbuckets); n++) {
                    buf4_mat_irrep_rd_block_sp(InBuf, h, n * rows_per_bucket, rows_per_bucket);
                    for (i=0; i < rows_per_bucket; i++){
                        for (j=0; j< coltot; j++){
                            tmp_matrix[h][i][j] = static_cast<double>(InBuf->matrix[h][i][j]);
                        }
                    }
                    memcpy((void *)&(OutBuf.matrix[h][0][0]), (const void *)&(tmp_matrix[h][0][0]),
                           ((long)sizeof(double)) * size);

                    buf4_mat_irrep_wrt_block(&OutBuf, h, n * rows_per_bucket, rows_per_bucket);
                }
                if (rows_left) {
                    size = ((long)rows_left) * ((long)coltot);

                    buf4_mat_irrep_rd_block_sp(InBuf, h, n * rows_per_bucket, rows_left);
                    for (i=0; i < rows_left; i++){
                        for (j=0; j< coltot; j++){
                            tmp_matrix[h][i][j] = static_cast<double>(InBuf->matrix[h][i][j]);
                        }
                    }
                    memcpy((void *)&(OutBuf.matrix[h][0][0]), (const void *)&(tmp_matrix[h][0][0]),
                           ((long)sizeof(double)) * size);

                    buf4_mat_irrep_wrt_block(&OutBuf, h, n * rows_per_bucket, rows_left);
                }

                buf4_mat_irrep_close_block_sp(InBuf, h, rows_per_bucket);
                buf4_mat_irrep_close_block(&OutBuf, h, rows_per_bucket);
            }
        }
    }
    free(tmp_matrix);
    buf4_close(&OutBuf);

    return 0;
}


}  // namespace psi
