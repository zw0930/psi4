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
#include <cstdlib>
#include <cmath>
#include "psi4/libqt/qt.h"
#include "dpd.h"
#include "psi4/psi4-dec.h"
#include "psi4/libpsi4util/PsiOutStream.h"
namespace psi {
int DPD::contract222_mp(dpdfile2<float> *X, dpdfile2<float> *Y, dpdfile2<double> *Z, int target_X, int target_Y, float alpha, double beta) {
    int h, nirreps, Xtrans, Ytrans, *numlinks;
    int GX, GY, GZ;
    int Hx, Hy, Hz;
    int symlink;
    int hz, row, col;
    float **TMP;

#ifdef DPD_DEBUG
    int *xrow, *xcol, *yrow, *ycol, *zrow, *zcol;
#endif

    nirreps = X->params->nirreps;
    GX = X->my_irrep;
    GY = Y->my_irrep;
    GZ = Z->my_irrep;

    file2_mat_init_sp(X);
    file2_mat_rd_sp(X);
    file2_mat_init_sp(Y);
    file2_mat_rd_sp(Y);
    //file2_mat_init_sp(Z_tmp);
    file2_mat_init(Z);
    if (std::fabs(beta) > 0.0) file2_mat_rd(Z);
    
    
    if (target_X == 0) {
        Xtrans = 0;
        numlinks = X->params->coltot;
        symlink = GX;
    } else if (target_X == 1) {
        Xtrans = 1;
        numlinks = X->params->rowtot;
        symlink = 0;
    } else {
        outfile->Printf("Junk X index %d in contract222\n", target_X);
        exit(PSI_RETURN_FAILURE);
    }
    if (target_Y == 0)
        Ytrans = 1;
    else if (target_Y == 1)
        Ytrans = 0;
    else {
        outfile->Printf("Junk Y index %d in contract222\n", target_Y);
        exit(PSI_RETURN_FAILURE);
    }

#ifdef DPD_DEBUG
    if (Xtrans) {
        xrow = X->params->coltot;
        xcol = X->params->rowtot;
    } else {
        xrow = X->params->rowtot;
        xcol = X->params->coltot;
    }

    if (Ytrans) {
        yrow = Y->params->coltot;
        ycol = Y->params->rowtot;
    } else {
        yrow = Y->params->rowtot;
        ycol = Y->params->coltot;
    }

    zrow = Z->params->rowtot;
    zcol = Z->params->coltot;

    if ((zrow != xrow) || (zcol != ycol) || (xcol != yrow)) {
        outfile->Printf("** Alignment error in contract222 **\n");
        dpd_error("dpd_contract222", "outfile");
    }
#endif

    /* loop over row irreps of X */
    for (Hx = 0; Hx < nirreps; Hx++) {
        if ((!Xtrans) && (!Ytrans)) {
            Hy = Hx ^ GX;
            Hz = Hx;
        } else if ((!Xtrans) && (Ytrans)) {
            Hy = Hx ^ GX ^ GY;
            Hz = Hx;
        } else if ((Xtrans) && (!Ytrans)) {
            Hy = Hx;
            Hz = Hx ^ GX;
        } else /*(( Xtrans)&&( Ytrans))*/ {
            Hy = Hx ^ GY;
            Hz = Hx ^ GX;
        }

        // TMP: temporarily holds the product of X and Y in single-precision
        TMP = dpd_block_matrix_sp(Z->params->rowtot[Hz], Z->params->coltot[Hz,GZ]);
        if (Z->params->rowtot[Hz] && Z->params->coltot[Hz ^ GZ] && numlinks[Hx ^ symlink]) {
            C_SGEMM(Xtrans ? 't' : 'n', Ytrans ? 't' : 'n', Z->params->rowtot[Hz], Z->params->coltot[Hz ^ GZ],
                    numlinks[Hx ^ symlink], alpha, &(X->matrix[Hx][0][0]), X->params->coltot[Hx ^ GX],
                    &(Y->matrix[Hy][0][0]), Y->params->coltot[Hy ^ GY], 0, &(TMP[0][0]),
                    Z->params->coltot[Hz ^ GZ]);
        }
       
        // Cast and add the matrix to Z
        for (row = 0; row < Z->params->rowtot[Hz]; row++){
		for (col= 0; col < Z->params->coltot[Hz,GZ]; col++){
		    Z->matrix[Hz][row][col] = beta * Z->matrix[Hz][row][col] + static_cast<double>(TMP[row][col]);
		}
	} 
	
        free_dpd_block_sp(TMP, Z->params->rowtot[Hz], Z->params->coltot[Hz,GZ]);
        /*
    newmm(X->matrix[Hx], Xtrans, Y->matrix[Hy], Ytrans, Z->matrix[Hz],
        Z->params->rowtot[Hz], numlinks[Hx^symlink], Z->params->coltot[Hz^GZ],
        alpha, beta);
    */
    } // loop over row irreps of X

    file2_mat_wrt(Z);
    file2_mat_close_sp(X);
    file2_mat_close_sp(Y);
    file2_mat_close(Z);
    //file2_mat_close(Z_tmp);

    return 0;
}
}  // namespace psi
