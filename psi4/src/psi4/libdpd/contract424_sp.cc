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
#include "psi4/libpsi4util/PsiOutStream.h"

namespace psi {

/* dpd_contract424(): Contracts four-index and two-index quantities to
 ** give a product four-index quantity.
 **
 ** Arguments:
 **   dpdbuf4<double>*X: A pointer to the leftmost dpd four-index
 **               buffer in the product.
 **   dpdfile2<double>*Y: A pointer to the rightmost dpd two-index
 **                file in the product.
 **   dpdbuf4<double>*Z: A pointer to the dpd four-index buffer target.
 **   int sum_X: Indicates which index (values of 0, 1, 2, and 3) of X
 **              is to be summed.
 **   int sum_Y: Indicates which index (values of 0 and 1) of Y is to be summed.
 **   int Ztrans: Boolean to indicate whether the final product must
 **                be bra-ket transposed in order for its indices to
 **                match those of the target, Z.
 **   double alpha: A prefactor for the product alpha * X * Y.
 **   double beta: A prefactor for the target beta * Z.
 */
int DPD::contract424_sp(dpdbuf4<float> *X, dpdfile2<float> *Y, dpdbuf4<float> *Z, int sum_X, int sum_Y, int Ztrans, float alpha, float beta) {
    int h, nirreps, GX, GY, GZ, hxbuf, hzbuf, h0, Hx, Hy, Hz, GsX, GsZ;
    int rking = 0, symlink;
    int Xtrans, Ytrans;
    int *numlinks, *numrows, *numcols;
    int incore;
    long int core, memoryd, core_total, rowtot, coltot, maxrows;
    int xcount, zcount, scount, Ysym;
    int rowx, rowz, colx, colz;
    int pq, rs, r, s, Gr, Gs;
    dpdtrans4<float> Xt, Zt;
    float ***Xmat;
    float ***Zmat;
#ifdef DPD_DEBUG
    int *xrow, *xcol, *yrow, *ycol, *zrow, *zcol;
#endif

    nirreps = X->params->nirreps;
    GX = X->file.my_irrep;
    GY = Y->my_irrep;
    GZ = Z->file.my_irrep;

    memoryd = dpd_main.memory;
    incore = 1; /* default */

    file2_mat_init_sp(Y);
    file2_mat_rd_sp(Y);

    if (sum_Y == 0) {
        Ytrans = 0;
        numlinks = Y->params->rowtot;
        symlink = 0;
    } else if (sum_Y == 1) {
        Ytrans = 1;
        numlinks = Y->params->coltot;
        symlink = GY;
    } else {
        outfile->Printf("Junk Y index %d\n", sum_Y);
        exit(PSI_RETURN_FAILURE);
    }

    if ((sum_X == 1) || (sum_X == 2)) trans4_init_sp(&Xt, X);

    if (Ztrans) trans4_init_sp(&Zt, Z);

    /*  if(std::fabs(beta) > 0.0) dpd_buf4_scm(Z, beta); */
    buf4_scm_sp(Z, beta);

#ifdef DPD_DEBUG
    if (Ytrans) {
        yrow = Y->params->coltot;
        ycol = Y->params->rowtot;
    } else {
        yrow = Y->params->rowtot;
        ycol = Y->params->coltot;
    }
#endif

    for (hxbuf = 0; hxbuf < nirreps; hxbuf++) {
        incore = 1; /* default */

        if (sum_X < 2) {
            if (!Ztrans)
                hzbuf = hxbuf ^ GX;
            else
                hzbuf = hxbuf ^ GX ^ GZ;
        } else {
            if (!Ztrans)
                hzbuf = hxbuf;
            else
                hzbuf = hxbuf ^ GZ;
        }

        /* Compute the core requirements for the straight contraction */
        core_total = 0;
        /** X terms **/
        coltot = X->params->coltot[hxbuf ^ GX];
        if (coltot) {
            maxrows = DPD_BIGNUM / coltot;
            if (maxrows < 1) {
                outfile->Printf("\nLIBDPD Error: cannot compute even the number of rows in contract424.\n");
                dpd_error("contract424", "outfile");
            }
        } else
            maxrows = DPD_BIGNUM;
        rowtot = X->params->rowtot[hxbuf];
        for (; rowtot > maxrows; rowtot -= maxrows) {
            if (core_total > (core_total + maxrows * coltot))
                incore = 0;
            else
                core_total += maxrows * coltot;
        }
        if (core_total > (core_total + rowtot * coltot)) incore = 0;
        core_total += rowtot * coltot;

        if (sum_X == 1 || sum_X == 2) core_total *= 2; /* we need room to transpose the X buffer */

        /** Z terms **/
        coltot = Z->params->coltot[hzbuf ^ GZ];
        if (coltot) {
            maxrows = DPD_BIGNUM / coltot;
            if (maxrows < 1) {
                outfile->Printf("\nLIBDPD Error: cannot compute even the number of rows in contract424.\n");
                dpd_error("contract424", "outfile");
            }
        } else
            maxrows = DPD_BIGNUM;
        rowtot = Z->params->rowtot[hzbuf];
        for (; rowtot > maxrows; rowtot -= maxrows) {
            if (core_total > (core_total + maxrows * coltot))
                incore = 0;
            else
                core_total += maxrows * coltot;
        }
        if (core_total > (core_total + rowtot * coltot)) incore = 0;
        core_total += rowtot * coltot;

        if (core_total > memoryd) incore = 0;

        /* Force incore for all but a "normal" 424 contraction for now */
        if (Ztrans || sum_X == 0 || sum_X == 1 || sum_X == 2) incore = 1;

        if (incore) {
            /*       dpd_buf4_scm(Z, beta); */
            buf4_mat_irrep_init_sp(Z, hzbuf);
            if (std::fabs(beta) > 0.0) buf4_mat_irrep_rd_sp(Z, hzbuf);
            if (Ztrans) {
                trans4_mat_irrep_init_sp(&Zt, hzbuf);
                trans4_mat_irrep_rd_sp(&Zt, hzbuf);
                buf4_mat_irrep_close_sp(Z, hzbuf);
                trans4_mat_irrep_shift31_sp(&Zt, hzbuf);
                rking = 1;
                numrows = Zt.shift.rowtot[hzbuf];
                numcols = Zt.shift.coltot[hzbuf];
                Zmat = Zt.shift.matrix[hzbuf];
#ifdef DPD_DEBUG
                zrow = Zt.shift.rowtot[hzbuf];
                zcol = Zt.shift.coltot[hzbuf];
#endif
            } else {
                buf4_mat_irrep_shift31_sp(Z, hzbuf);
                rking = 1;
                numrows = Z->shift.rowtot[hzbuf];
                numcols = Z->shift.coltot[hzbuf];
                Zmat = Z->shift.matrix[hzbuf];
#ifdef DPD_DEBUG
                zrow = Z->shift.rowtot[hzbuf];
                zcol = Z->shift.coltot[hzbuf];
#endif
            }

            if (sum_X == 0) {
                buf4_mat_irrep_init_sp(X, hxbuf);
                buf4_mat_irrep_rd_sp(X, hxbuf);
                buf4_mat_irrep_shift13_sp(X, hxbuf);
                Xmat = X->shift.matrix[hxbuf];
                Xtrans = 1;
#ifdef DPD_DEBUG
                xrow = X->shift.coltot[hxbuf];
                xcol = X->shift.rowtot[hxbuf];
#endif
            } else if (sum_X == 1) {
                buf4_mat_irrep_init_sp(X, hxbuf);
                buf4_mat_irrep_rd_sp(X, hxbuf);
                trans4_mat_irrep_init_sp(&Xt, hxbuf);
                trans4_mat_irrep_rd_sp(&Xt, hxbuf);
                buf4_mat_irrep_close_sp(X, hxbuf);
                trans4_mat_irrep_shift31_sp(&Xt, hxbuf);
                rking = 1;
                Xmat = Xt.shift.matrix[hxbuf];
                Xtrans = 0;
#ifdef DPD_DEBUG
                xrow = Xt.shift.rowtot[hxbuf];
                xcol = Xt.shift.coltot[hxbuf];
#endif
            } else if (sum_X == 2) {
                buf4_mat_irrep_init_sp(X, hxbuf);
                buf4_mat_irrep_rd_sp(X, hxbuf);
                trans4_mat_irrep_init_sp(&Xt, hxbuf);
                trans4_mat_irrep_rd_sp(&Xt, hxbuf);
                buf4_mat_irrep_close_sp(X, hxbuf);
                trans4_mat_irrep_shift13_sp(&Xt, hxbuf);
                Xmat = Xt.shift.matrix[hxbuf];
                Xtrans = 1;
#ifdef DPD_DEBUG
                xrow = Xt.shift.coltot[hxbuf];
                xcol = Xt.shift.rowtot[hxbuf];
#endif
            } else if (sum_X == 3) {
                buf4_mat_irrep_init_sp(X, hxbuf);
                buf4_mat_irrep_rd_sp(X, hxbuf);
                buf4_mat_irrep_shift31_sp(X, hxbuf);
                rking = 1;
                Xmat = X->shift.matrix[hxbuf];
                Xtrans = 0;
#ifdef DPD_DEBUG
                xrow = X->shift.rowtot[hxbuf];
                xcol = X->shift.coltot[hxbuf];
#endif
            }

            if (rking)
                for (Hz = 0; Hz < nirreps; Hz++) {
                    if (!Xtrans && !Ytrans) {
                        Hx = Hz;
                        Hy = Hz ^ GX;
                    } else if (!Xtrans && Ytrans) {
                        Hx = Hz;
                        Hy = Hz ^ GX ^ GY;
                    } else if (Xtrans && !Ytrans) {
                        Hx = Hz ^ GX;
                        Hy = Hz ^ GX;
                    } else if (Xtrans && Ytrans) {
                        Hx = Hz ^ GX;
                        Hy = Hz ^ GX ^ GY;
                    }
#ifdef DPD_DEBUG
                    if ((xrow[Hz] != zrow[Hz]) || (ycol[Hz] != zcol[Hz]) || (xcol[Hz] != yrow[Hz])) {
                        outfile->Printf("** Alignment error in contract424 **\n");
                        outfile->Printf("** Irrep: %d; Subirrep: %d **\n", hxbuf, Hz);
                        dpd_error("dpd_contract424", "outfile");
                    }
#endif
                    /* outfile->Printf("Hx %d Hy %d Hz %d\n",Hx,Hy,Hz);
         outfile->Printf("numrows %d numlinks %d numcols %d\n",numrows[Hz],numlinks[Hy],numcols[Hz]); */
                    newmm_rking_sp(Xmat[Hx], Xtrans, Y->matrix[Hy], Ytrans, Zmat[Hz], numrows[Hz], numlinks[Hy ^ symlink],
                                numcols[Hz], alpha, 1.0);
                }
            else
                for (Hz = 0; Hz < nirreps; Hz++) {
                    if (!Xtrans && !Ytrans) {
                        Hx = Hz;
                        Hy = Hz ^ GX;
                    } else if (!Xtrans && Ytrans) {
                        Hx = Hz;
                        Hy = Hz ^ GX ^ GY;
                    } else if (Xtrans && !Ytrans) {
                        Hx = Hz ^ GX;
                        Hy = Hz ^ GX;
                    } else if (Xtrans && Ytrans) {
                        Hx = Hz ^ GX;
                        Hy = Hz ^ GX ^ GY;
                    }
#ifdef DPD_DEBUG
                    if ((xrow[Hz] != zrow[Hz]) || (ycol[Hz] != zcol[Hz]) || (xcol[Hz] != yrow[Hz])) {
                        outfile->Printf("** Alignment error in contract424 **\n");
                        outfile->Printf("** Irrep: %d; Subirrep: %d **\n", hxbuf, Hz);
                        dpd_error("dpd_contract424", "outfile");
                    }
#endif
                    if (numrows[Hz] && numcols[Hz] && numlinks[Hy ^ symlink]) {
                        if (!Xtrans && !Ytrans) {
                            C_SGEMM('n', 'n', numrows[Hz], numcols[Hz], numlinks[Hy ^ symlink], alpha,
                                    &(Xmat[Hz][0][0]), numlinks[Hy ^ symlink], &(Y->matrix[Hy][0][0]), numcols[Hz], 1.0,
                                    &(Zmat[Hz][0][0]), numcols[Hz]);
                        } else if (Xtrans && !Ytrans) {
                            C_SGEMM('t', 'n', numrows[Hz], numcols[Hz], numlinks[Hy ^ symlink], alpha,
                                    &(Xmat[Hz][0][0]), numrows[Hz], &(Y->matrix[Hy][0][0]), numcols[Hz], 1.0,
                                    &(Zmat[Hz][0][0]),  numcols[Hz]);
                        } else if (!Xtrans && Ytrans) {
                            C_SGEMM('n', 't', numrows[Hz], numcols[Hz], numlinks[Hy ^ symlink], alpha,
                                    &(Xmat[Hz][0][0]), numlinks[Hy ^ symlink], &(Y->matrix[Hy][0][0]),
                                    numlinks[Hy ^ symlink], 1.0, &(Zmat[Hz][0][0]), numcols[Hz]);
                        } else {
                            C_SGEMM('t', 't', numrows[Hz], numcols[Hz], numlinks[Hy ^ symlink], alpha,
                                    &(Xmat[Hz][0][0]), numrows[Hz], &(Y->matrix[Hy][0][0]), numlinks[Hy ^ symlink], 1.0,
                                    &(Zmat[Hz][0][0]), numcols[Hz]);
                        }
                    }
                }

            if (sum_X == 0)
                buf4_mat_irrep_close_sp(X, hxbuf);
            else if (sum_X == 1)
                trans4_mat_irrep_close_sp(&Xt, hxbuf);
            else if (sum_X == 2)
                trans4_mat_irrep_close_sp(&Xt, hxbuf);
            else if (sum_X == 3)
                buf4_mat_irrep_close_sp(X, hxbuf);

            if (Ztrans) {
                buf4_mat_irrep_init_sp(Z, hzbuf);
                trans4_mat_irrep_wrt_sp(&Zt, hzbuf);
                trans4_mat_irrep_close_sp(&Zt, hzbuf);
            }

            buf4_mat_irrep_wrt_sp(Z, hzbuf);
            buf4_mat_irrep_close_sp(Z, hzbuf);

        }      /* end if(incore) */
        else { /* out-of-core for "normal" 424 contractions */
               /* Prepare the input buffer for the X factor and the target*/
#ifdef DPD_DEBUG
            outfile->Printf("\t424 out-of-core: %d\n", hxbuf);
#endif
            buf4_mat_irrep_row_init_sp(X, hxbuf);
            buf4_mat_irrep_row_init_sp(Z, hzbuf);

            /* Loop over rows of the X factor and the target */
            for (pq = 0; pq < Z->params->rowtot[hzbuf]; pq++) {
                buf4_mat_irrep_row_zero_sp(X, hxbuf, pq);
                buf4_mat_irrep_row_rd_sp(X, hxbuf, pq);

                buf4_mat_irrep_row_zero_sp(Z, hzbuf, pq);

                if (std::fabs(beta) > 0.0) buf4_mat_irrep_row_rd_sp(Z, hzbuf, pq);

                xcount = zcount = 0;

                for (Gr = 0; Gr < nirreps; Gr++) {
                    GsX = Gr ^ hxbuf ^ GX;
                    GsZ = Gr ^ hzbuf ^ GZ;

                    rowx = X->params->rpi[Gr];
                    colx = X->params->spi[GsX];
                    rowz = Z->params->rpi[Gr];
                    colz = Z->params->spi[GsZ];

                    if (rowx && colx && colz) {
                        C_SGEMM('n', Ytrans ? 't' : 'n', rowx, colz, colx, alpha, &(X->matrix[hxbuf][0][xcount]), colx,
                                &(Y->matrix[Ytrans ? GsZ : GsX][0][0]), Ytrans ? colx : colz, 1.0,
                                &(Z->matrix[hzbuf][0][zcount]), colz);
                    }

                    xcount += rowx * colx;
                    zcount += rowz * colz;
                }

                buf4_mat_irrep_row_wrt_sp(Z, hzbuf, pq);
            }

            buf4_mat_irrep_row_close_sp(X, hxbuf);
            buf4_mat_irrep_row_close_sp(Z, hzbuf);
        }
    }

    if ((sum_X == 1) || (sum_X == 2)) trans4_close_sp(&Xt);

    if (Ztrans) trans4_close_sp(&Zt);

    file2_mat_close_sp(Y);

    return 0;
}

}  // namespace psi
