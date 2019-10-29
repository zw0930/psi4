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
    \ingroup CCENERGY
    \brief Enter brief description of file here
*/
#include <cstdio>
#include <cstdlib>
#include "psi4/libdpd/dpd.h"
#include "psi4/libqt/qt.h"
#include "MOInfo.h"
#include "Params.h"
#include "psi4/libciomr/libciomr.h"
#include "psi4/cc/ccwave.h"


// Didn't finish ROHF & UHF

namespace psi {
namespace ccenergy {

void CCEnergyWavefunction::Fae_build_sp() {
    int a, e, nirreps;
    int ef, m, f, M, A, Gm, Ga, Ge, Gf, nrows, ncols;
    float *X;
    dpdfile2<float> tIA_sp, tia_sp;
    dpdfile2<float> FME_sp, Fme_sp;
    dpdfile2<float> fAB_sp, fab_sp, fIA_sp, fia_sp;
    dpdfile2<float> FAE_sp, Fae_sp;
    dpdfile2<float> FAEt_sp, Faet_sp;
    dpdbuf4<float> F_anti_sp, F_sp, D_anti_sp, D_sp;
    dpdbuf4<float> tautIJAB_sp, tautijab_sp, tautIjAb_sp, taut;

    nirreps = moinfo_.nirreps;

    if (params_.ref == 0) { /** RHF **/
        global_dpd_->file2_initi_sp(&fAB_sp, PSIF_CC_OEI, 0, 1, 1, "fAB_sp");
        global_dpd_->file2_copy_sp(&fAB_sp, PSIF_CC_OEI, "FAE_sp");
        global_dpd_->file2_close_sp(&fAB_sp);
    } else if (params_.ref == 1) { /** ROHF **/
        global_dpd_->file2_init_sp(&fAB_sp, PSIF_CC_OEI, 0, 1, 1, "fAB_sp");
        global_dpd_->file2_copy_sp(&fAB_sp, PSIF_CC_OEI, "FAE_sp");
        global_dpd_->file2_close_sp(&fAB_sp);

        global_dpd_->file2_init_sp(&fab_sp, PSIF_CC_OEI, 0, 1, 1, "fab_sp");
        global_dpd_->file2_copy_sp(&fab_sp, PSIF_CC_OEI, "Fae_sp");
        global_dpd_->file2_close_sp(&fab_sp);
    } else if (params_.ref == 2) { /** UHF **/
        global_dpd_->file2_init_sp(&fAB_sp, PSIF_CC_OEI, 0, 1, 1, "fAB_sp");
        global_dpd_->file2_copy_sp(&fAB_sp, PSIF_CC_OEI, "FAE_sp");
        global_dpd_->file2_close_sp(&fAB_sp);

        global_dpd_->file2_init_sp(&fab_sp, PSIF_CC_OEI, 0, 3, 3, "fab_sp");
        global_dpd_->file2_copy_sp(&fab_sp, PSIF_CC_OEI, "Fae_sp");
        global_dpd_->file2_close_sp(&fab_sp);
    }

    if (params_.ref == 0) { /** RHF **/
        global_dpd_->file2_init_sp(&FAE_sp, PSIF_CC_OEI, 0, 1, 1, "FAE_sp");

        global_dpd_->file2_mat_init_sp(&FAE_sp);
        global_dpd_->file2_mat_rd_sp(&FAE_sp);

        /*
          for(int h = 0; h < moinfo.nirreps; h++) {
          for(int a = 0; a < FAE.params->rowtot[h]; a++)
          FAE.matrix[h][a][a] = 0;
          }
        */

        global_dpd_->file2_mat_wrt_sp(&FAE_sp);
        global_dpd_->file2_mat_close_sp(&FAE_sp);
        global_dpd_->file2_close_sp(&FAE_sp);
    } else if (params_.ref == 1) { /** ROHF **/
        global_dpd_->file2_init_sp(&FAE_sp, PSIF_CC_OEI, 0, 1, 1, "FAE_sp");
        global_dpd_->file2_init_sp(&Fae_sp, PSIF_CC_OEI, 0, 1, 1, "Fae_sp");

        global_dpd_->file2_mat_init_sp(&FAE_sp);
        global_dpd_->file2_mat_rd_sp(&FAE_sp);
        global_dpd_->file2_mat_init_sp(&Fae_sp);
        global_dpd_->file2_mat_rd_sp(&Fae_sp);

        for (int h = 0; h < moinfo_.nirreps; h++) {
            for (int a = 0; a < FAE_sp.params->rowtot[h]; a++)
                for (int e = 0; e < FAE_sp.params->coltot[h]; e++) FAE_sp.matrix[h][a][e] *= (1 - (a == e));

            for (int a = 0; a < Fae_sp.params->rowtot[h]; a++)
                for (int e = 0; e < Fae_sp.params->coltot[h]; e++) Fae_sp.matrix[h][a][e] *= (1 - (a == e));
        }

        global_dpd_->file2_mat_wrt_sp(&FAE_sp);
        global_dpd_->file2_mat_close_sp(&FAE_sp);
        global_dpd_->file2_mat_wrt_sp(&Fae_sp);
        global_dpd_->file2_mat_close_sp(&Fae_sp);

        global_dpd_->file2_close_sp(&FAE_sp);
        global_dpd_->file2_close_sp(&Fae_sp);
    } else if (params_.ref == 2) { /** UHF **/
        global_dpd_->file2_init_sp(&FAE_sp, PSIF_CC_OEI, 0, 1, 1, "FAE_sp");
        global_dpd_->file2_init_sp(&Fae_sp, PSIF_CC_OEI, 0, 3, 3, "Fae_sp");

        global_dpd_->file2_mat_init_sp(&FAE_sp);
        global_dpd_->file2_mat_rd_sp(&FAE_sp);
        global_dpd_->file2_mat_init_sp(&Fae_sp);
        global_dpd_->file2_mat_rd_sp(&Fae_sp);

        for (int h = 0; h < moinfo_.nirreps; h++) {
            for (int a = 0; a < FAE_sp.params->rowtot[h]; a++)
                for (int e = 0; e < FAE_sp.params->coltot[h]; e++) FAE_sp.matrix[h][a][e] *= (1 - (a == e));

            for (int a = 0; a < Fae_sp.params->rowtot[h]; a++)
                for (int e = 0; e < Fae_sp.params->coltot[h]; e++) Fae_sp.matrix[h][a][e] *= (1 - (a == e));
        }

        global_dpd_->file2_mat_wrt_sp(&FAE_sp);
        global_dpd_->file2_mat_close_sp(&FAE_sp);
        global_dpd_->file2_mat_wrt_sp(&Fae_sp);
        global_dpd_->file2_mat_close_sp(&Fae_sp);

        global_dpd_->file2_close_sp(&FAE_sp);
        global_dpd_->file2_close_sp(&Fae_sp);
    }

    if (params_.ref == 0) { /** RHF **/
        global_dpd_->file2_init_sp(&FAE_sp, PSIF_CC_OEI, 0, 1, 1, "FAE_sp");
        global_dpd_->file2_init_sp(&fIA_sp, PSIF_CC_OEI, 0, 0, 1, "fIA_sp");
        global_dpd_->file2_init_sp(&tIA_sp, PSIF_CC_OEI, 0, 0, 1, "tIA_sp");
        global_dpd_->contract222_sp(&tIA_sp, &fIA_sp, &FAE_sp, 1, 1, -0.5, 1);
        global_dpd_->file2_close_sp(&tIA_sp);
        global_dpd_->file2_close_sp(&fIA_sp);
        global_dpd_->file2_close_sp(&FAE_sp);

        /*
          dpd_file2_init(&FAE, CC_OEI, 0, 1, 1, "FAE");
          dpd_buf4_init(&F_anti, CC_FINTS, 0, 10, 5, 10, 5, 1, "F <ia|bc>");
          dpd_buf4_init(&F, CC_FINTS, 0, 10, 5, 10, 5, 0,"F <ia|bc>");
          dpd_file2_init(&tIA, CC_OEI, 0, 0, 1, "tIA");
          dpd_dot13(&tIA, &F_anti, &FAE, 0, 0, 1.0, 1.0);
          dpd_dot13(&tIA, &F, &FAE, 0, 0, 1.0, 1.0);
          dpd_file2_close(&tIA);
          dpd_buf4_close(&F_anti);
          dpd_buf4_close(&F);
          dpd_file2_close(&FAE);
        */

        /* Out-of-core algorithm for F->FAE added 3/20/05 - TDC */
        /* Fae <-- t(m,f) [2 <ma|fe> - <ma|ef>] */
        global_dpd_->file2_init_sp(&FAE_sp, PSIF_CC_OEI, 0, 1, 1, "FAE_sp");
        global_dpd_->file2_mat_init_sp(&FAE_sp);
        global_dpd_->file2_mat_rd_sp(&FAE_sp);
        global_dpd_->file2_init_sp(&tIA_sp, PSIF_CC_OEI, 0, 0, 1, "tIA_sp");
        global_dpd_->file2_mat_init_sp(&tIA_sp);
        global_dpd_->file2_mat_rd_sp(&tIA_sp);
        global_dpd_->buf4_init_sp(&F_sp, PSIF_CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc> sp");
        for (int Gma = 0; Gma < nirreps; Gma++) {
            global_dpd_->buf4_mat_irrep_row_init_sp(&F_sp, Gma);
            X = init_array(F_sp.params->coltot[Gma]);

            for (int ma = 0; ma < F_sp.params->rowtot[Gma]; ma++) {
                global_dpd_->buf4_mat_irrep_row_rd_sp(&F_sp, Gma, ma);
                m = F_sp.params->roworb[Gma][ma][0];
                a = F_sp.params->roworb[Gma][ma][1];
                Gm = F_sp.params->psym[m];
                Ga = Ge = Gm ^ Gma; /* Fae is totally symmetric */
                Gf = Gm;            /* T1 is totally symmetric */
                M = m - F_sp.params->poff[Gm];
                A = a - F_sp.params->qoff[Ga];

                zero_arr_sp(X, F_sp.params->coltot[Gma]);

                /* build spin-adapted F-integrals for current ma */
                for (int fe = 0; fe < F_sp.params->coltot[Gma]; fe++) {
                    f = F_sp.params->colorb[Gma][fe][0];
                    e = F_sp.params->colorb[Gma][fe][1];
                    ef = F_sp.params->colidx[e][f];
                    X[fe] = 2.0 * F_sp.matrix[Gma][0][fe] - F_sp.matrix[Gma][0][ef];
                }

                nrows = moinfo_.virtpi[Gf];
                ncols = moinfo_.virtpi[Ge];
                if (nrows && ncols)
                    C_SGEMV('t', nrows, ncols, 1.0, &X[F.col_offset[Gma][Gf]], ncols, tIA_sp.matrix[Gm][M], 1, 1.0,
                            FAE_sp.matrix[Ga][A], 1);
            }

            free(X);
            global_dpd_->buf4_mat_irrep_row_close_sp(&F_sp, Gma);
        }
        global_dpd_->buf4_close_sp(&F_sp);
        global_dpd_->file2_mat_close_sp(&tIA_sp);
        global_dpd_->file2_close_sp(&tIA_sp);
        global_dpd_->file2_mat_wrt_sp(&FAE_sp);
        global_dpd_->file2_mat_close_sp(&FAE_sp);
        global_dpd_->file2_close_sp(&FAE_sp);

        global_dpd_->file2_init_sp(&FAE_sp, PSIF_CC_OEI, 0, 1, 1, "FAE_sp");

        global_dpd_->buf4_init_sp(&D_sp, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D 2<ij|ab> - <ij|ba> sp");
        global_dpd_->buf4_init_sp(&tautIjAb_sp, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tautIjAb_sp");
        global_dpd_->contract442_sp(&tautIjAb_sp, &D_sp, &FAE_sp, 3, 3, -1, 1);
        global_dpd_->buf4_close_sp(&D_sp);
        global_dpd_->buf4_close_sp(&tautIjAb_sp);

        /* Build the tilde intermediates */
        global_dpd_->file2_copy_sp(&FAE_sp, PSIF_CC_OEI, "FAEt_sp");
        global_dpd_->file2_close_sp(&FAE_sp);

        global_dpd_->file2_init_sp(&FAEt_sp, PSIF_CC_OEI, 0, 1, 1, "FAEt_sp");

        global_dpd_->file2_init_sp(&tIA_sp, PSIF_CC_OEI, 0, 0, 1, "tIA_sp");
        global_dpd_->file2_init_sp(&FME_sp, PSIF_CC_OEI, 0, 0, 1, "FME_sp");
        global_dpd_->contract222_sp(&tIA_sp, &FME_sp, &FAEt_sp, 1, 1, -0.5, 1);
        global_dpd_->file2_close_sp(&tIA_sp);
        global_dpd_->file2_close_sp(&FME_sp);

        global_dpd_->file2_close_sp(&FAEt_sp);
    } else if (params_.ref == 1) { /** ROHF **/
        global_dpd_->file2_init(&FAE, PSIF_CC_OEI, 0, 1, 1, "FAE");
        global_dpd_->file2_init(&Fae, PSIF_CC_OEI, 0, 1, 1, "Fae");

        global_dpd_->file2_init(&fIA, PSIF_CC_OEI, 0, 0, 1, "fIA");
        global_dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
        global_dpd_->contract222(&tIA, &fIA, &FAE, 1, 1, -0.5, 1);
        global_dpd_->file2_close(&tIA);
        global_dpd_->file2_close(&fIA);

        global_dpd_->file2_init(&fia, PSIF_CC_OEI, 0, 0, 1, "fia");
        global_dpd_->file2_init(&tia, PSIF_CC_OEI, 0, 0, 1, "tia");
        global_dpd_->contract222(&tia, &fia, &Fae, 1, 1, -0.5, 1);
        global_dpd_->file2_close(&tia);
        global_dpd_->file2_close(&fia);

        global_dpd_->buf4_init(&F_anti, PSIF_CC_FINTS, 0, 10, 5, 10, 5, 1, "F <ia|bc>");
        global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc>");
        global_dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
        global_dpd_->file2_init(&tia, PSIF_CC_OEI, 0, 0, 1, "tia");

        global_dpd_->dot13(&tIA, &F_anti, &FAE, 0, 0, 1.0, 1.0);
        global_dpd_->dot13(&tia, &F, &FAE, 0, 0, 1.0, 1.0);

        global_dpd_->dot13(&tia, &F_anti, &Fae, 0, 0, 1.0, 1.0);
        global_dpd_->dot13(&tIA, &F, &Fae, 0, 0, 1.0, 1.0);

        global_dpd_->file2_close(&tIA);
        global_dpd_->file2_close(&tia);
        global_dpd_->buf4_close(&F_anti);
        global_dpd_->buf4_close(&F);

        global_dpd_->buf4_init(&D_anti, PSIF_CC_DINTS, 0, 2, 5, 2, 5, 0, "D <ij||ab> (i>j,ab)");

        global_dpd_->buf4_init(&tautIJAB, PSIF_CC_TAMPS, 0, 2, 5, 2, 7, 0, "tautIJAB");
        global_dpd_->buf4_init(&tautijab, PSIF_CC_TAMPS, 0, 2, 5, 2, 7, 0, "tautijab");

        global_dpd_->contract442(&tautIJAB, &D_anti, &FAE, 2, 2, -1, 1);
        global_dpd_->contract442(&tautijab, &D_anti, &Fae, 2, 2, -1, 1);

        global_dpd_->buf4_close(&D_anti);
        global_dpd_->buf4_close(&tautIJAB);
        global_dpd_->buf4_close(&tautijab);

        global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
        global_dpd_->buf4_init(&tautIjAb, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tautIjAb");

        global_dpd_->contract442(&tautIjAb, &D, &Fae, 3, 3, -1, 1);
        global_dpd_->contract442(&tautIjAb, &D, &FAE, 2, 2, -1, 1);

        global_dpd_->buf4_close(&D);
        global_dpd_->buf4_close(&tautIjAb);

        /* Build the tilde intermediates */
        global_dpd_->file2_copy(&FAE, PSIF_CC_OEI, "FAEt");
        global_dpd_->file2_copy(&Fae, PSIF_CC_OEI, "Faet");

        global_dpd_->file2_close(&FAE);
        global_dpd_->file2_close(&Fae);

        global_dpd_->file2_init(&FAEt, PSIF_CC_OEI, 0, 1, 1, "FAEt");
        global_dpd_->file2_init(&Faet, PSIF_CC_OEI, 0, 1, 1, "Faet");

        global_dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
        global_dpd_->file2_init(&FME, PSIF_CC_OEI, 0, 0, 1, "FME");
        global_dpd_->contract222(&tIA, &FME, &FAEt, 1, 1, -0.5, 1);
        global_dpd_->file2_close(&tIA);
        global_dpd_->file2_close(&FME);

        global_dpd_->file2_init(&tia, PSIF_CC_OEI, 0, 0, 1, "tia");
        global_dpd_->file2_init(&Fme, PSIF_CC_OEI, 0, 0, 1, "Fme");
        global_dpd_->contract222(&tia, &Fme, &Faet, 1, 1, -0.5, 1);
        global_dpd_->file2_close(&tia);
        global_dpd_->file2_close(&Fme);

        global_dpd_->file2_close(&FAEt);
        global_dpd_->file2_close(&Faet);
    } else if (params_.ref == 2) { /** UHF **/

        global_dpd_->file2_init(&FAE, PSIF_CC_OEI, 0, 1, 1, "FAE");
        global_dpd_->file2_init(&Fae, PSIF_CC_OEI, 0, 3, 3, "Fae");

        global_dpd_->file2_init(&fIA, PSIF_CC_OEI, 0, 0, 1, "fIA");
        global_dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
        global_dpd_->contract222(&tIA, &fIA, &FAE, 1, 1, -0.5, 1);
        global_dpd_->file2_close(&tIA);
        global_dpd_->file2_close(&fIA);

        global_dpd_->file2_init(&fia, PSIF_CC_OEI, 0, 2, 3, "fia");
        global_dpd_->file2_init(&tia, PSIF_CC_OEI, 0, 2, 3, "tia");
        global_dpd_->contract222(&tia, &fia, &Fae, 1, 1, -0.5, 1);
        global_dpd_->file2_close(&tia);
        global_dpd_->file2_close(&fia);

        global_dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
        global_dpd_->file2_init(&tia, PSIF_CC_OEI, 0, 2, 3, "tia");

        global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 20, 5, 20, 5, 1, "F <IA|BC>");
        global_dpd_->dot13(&tIA, &F, &FAE, 0, 0, 1, 1);
        global_dpd_->buf4_close(&F);

        global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 27, 29, 27, 29, 0, "F <iA|bC>");
        global_dpd_->dot13(&tia, &F, &FAE, 0, 0, 1, 1);
        global_dpd_->buf4_close(&F);

        global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 30, 15, 30, 15, 1, "F <ia|bc>");
        global_dpd_->dot13(&tia, &F, &Fae, 0, 0, 1, 1);
        global_dpd_->buf4_close(&F);

        global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 24, 28, 24, 28, 0, "F <Ia|Bc>");
        global_dpd_->dot13(&tIA, &F, &Fae, 0, 0, 1, 1);
        global_dpd_->buf4_close(&F);

        global_dpd_->file2_close(&tIA);
        global_dpd_->file2_close(&tia);

        global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 2, 5, 2, 5, 0, "D <IJ||AB> (I>J,AB)");
        global_dpd_->buf4_init(&taut, PSIF_CC_TAMPS, 0, 2, 5, 2, 7, 0, "tautIJAB");
        global_dpd_->contract442(&taut, &D, &FAE, 2, 2, -1, 1);
        global_dpd_->buf4_close(&taut);
        global_dpd_->buf4_close(&D);

        global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 22, 28, 22, 28, 0, "D <Ij|Ab>");
        global_dpd_->buf4_init(&taut, PSIF_CC_TAMPS, 0, 22, 28, 22, 28, 0, "tautIjAb");
        global_dpd_->contract442(&taut, &D, &FAE, 2, 2, -1, 1);
        global_dpd_->buf4_close(&taut);
        global_dpd_->buf4_close(&D);

        global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 12, 15, 12, 15, 0, "D <ij||ab> (i>j,ab)");
        global_dpd_->buf4_init(&taut, PSIF_CC_TAMPS, 0, 12, 15, 12, 17, 0, "tautijab");
        global_dpd_->contract442(&taut, &D, &Fae, 2, 2, -1, 1);
        global_dpd_->buf4_close(&taut);
        global_dpd_->buf4_close(&D);

        global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 22, 28, 22, 28, 0, "D <Ij|Ab>");
        global_dpd_->buf4_init(&taut, PSIF_CC_TAMPS, 0, 22, 28, 22, 28, 0, "tautIjAb");
        global_dpd_->contract442(&taut, &D, &Fae, 3, 3, -1, 1);
        global_dpd_->buf4_close(&taut);
        global_dpd_->buf4_close(&D);

        /* Build the tilde intermediates */
        global_dpd_->file2_copy(&FAE, PSIF_CC_OEI, "FAEt");
        global_dpd_->file2_copy(&Fae, PSIF_CC_OEI, "Faet");

        global_dpd_->file2_close(&FAE);
        global_dpd_->file2_close(&Fae);

        global_dpd_->file2_init(&FAEt, PSIF_CC_OEI, 0, 1, 1, "FAEt");
        global_dpd_->file2_init(&Faet, PSIF_CC_OEI, 0, 3, 3, "Faet");

        global_dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
        global_dpd_->file2_init(&FME, PSIF_CC_OEI, 0, 0, 1, "FME");
        global_dpd_->contract222(&tIA, &FME, &FAEt, 1, 1, -0.5, 1);
        global_dpd_->file2_close(&tIA);
        global_dpd_->file2_close(&FME);

        global_dpd_->file2_init(&tia, PSIF_CC_OEI, 0, 2, 3, "tia");
        global_dpd_->file2_init(&Fme, PSIF_CC_OEI, 0, 2, 3, "Fme");
        global_dpd_->contract222(&tia, &Fme, &Faet, 1, 1, -0.5, 1);
        global_dpd_->file2_close(&tia);
        global_dpd_->file2_close(&Fme);

        global_dpd_->file2_close(&FAEt);
        global_dpd_->file2_close(&Faet);
    }
}



}  // namespace ccenergy
}  // namespace psi
