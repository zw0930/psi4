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
#include "Params.h"
#include "MOInfo.h"
#include "psi4/cc/ccwave.h"

namespace psi {
namespace ccenergy {

void CCEnergyWavefunction::t1_build() {
    dpdfile2<double> newtIA, newtia, tIA, tia, fIA, fia;
    dpdfile2<double> FAE, Fae, FMI, Fmi, FME, Fme;
    dpdbuf4<double> tIJAB, tijab, tIjAb, tiJaB, T2;
    dpdbuf4<double> C, C_anti, D, F_anti, F, E_anti, E;
    int Gmi, Gm, Gi, Ga, m, a, A, nrows, ncols;

    if (params_.ref == 0) { /** RHF **/
        global_dpd_->file2_init(&fIA, PSIF_CC_OEI, 0, 0, 1, "fIA");
        global_dpd_->file2_copy(&fIA, PSIF_CC_OEI, "New tIA");
        global_dpd_->file2_close(&fIA);

        global_dpd_->file2_init(&newtIA, PSIF_CC_OEI, 0, 0, 1, "New tIA");

        global_dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");

        global_dpd_->file2_init(&FAE, PSIF_CC_OEI, 0, 1, 1, "FAE");
        global_dpd_->contract222(&tIA, &FAE, &newtIA, 0, 0, 1, 1);
        global_dpd_->file2_close(&FAE);

        global_dpd_->file2_init(&FMI, PSIF_CC_OEI, 0, 0, 0, "FMI");
        global_dpd_->contract222(&FMI, &tIA, &newtIA, 1, 1, -1, 1);
        global_dpd_->file2_close(&FMI);

        global_dpd_->file2_close(&tIA);

        global_dpd_->file2_init(&FME, PSIF_CC_OEI, 0, 0, 1, "FME");

        global_dpd_->buf4_init(&tIjAb, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "2 tIAjb - tIBja");
        global_dpd_->contract422(&tIjAb, &FME, &newtIA, 0, 0, 1, 1);
        global_dpd_->buf4_close(&tIjAb);

        global_dpd_->file2_close(&FME);

        global_dpd_->buf4_init(&C_anti, PSIF_CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia||jb>");
        global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");

        global_dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");

        global_dpd_->dot14(&tIA, &C_anti, &newtIA, 0, 1, -1, 1);
        global_dpd_->dot13(&tIA, &D, &newtIA, 0, 0, 1, 1);

        global_dpd_->file2_close(&tIA);

        global_dpd_->buf4_close(&C_anti);
        global_dpd_->buf4_close(&D);

        /*
          dpd_buf4_init(&Z, CC_TMP0, 0, 10, 0, 10, 0, 0, "Z(ma,mi)");
          dpd_buf4_init(&F, CC_FINTS, 0, 10, 5, 10, 5, 0, "F 2<ia|bc> - <ia|cb>");
          dpd_buf4_init(&tIjAb, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
          dpd_contract444(&F, &tIjAb, &Z, 0, 0, 1.0, 0.0);
          dpd_buf4_close(&tIjAb);
          dpd_buf4_close(&F);
          dpd_trace42_13(&Z, &newtIA, 1, 1.0, 1.0);
          dpd_buf4_close(&Z);
        */

        /* t(i,a) <-- (2 t(mi,ef) - t(mi,fe)) <ma|ef> */
        /* out-of-core version replacing the *stupid* code above 3/22/05, TDC */
        global_dpd_->file2_mat_init(&newtIA);
        global_dpd_->file2_mat_rd(&newtIA);
        global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "2 tIjAb - tIjBa");
        global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc>");
        for (int Gma = 0; Gma < moinfo_.nirreps; Gma++) {
            Gmi = Gma; /* T1 is totally symmetric */

            global_dpd_->buf4_mat_irrep_row_init(&F, Gma);
            global_dpd_->buf4_mat_irrep_init(&T2, Gmi);
            global_dpd_->buf4_mat_irrep_rd(&T2, Gmi);

            for (int ma = 0; ma < F.params->rowtot[Gma]; ma++) {
                global_dpd_->buf4_mat_irrep_row_rd(&F, Gma, ma);

                m = F.params->roworb[Gma][ma][0];
                a = F.params->roworb[Gma][ma][1];
                Gm = F.params->psym[m];
                Ga = F.params->qsym[a];
                Gi = Ga; /* T1 is totally symmetric */
                A = a - F.params->qoff[Ga];

                nrows = moinfo_.occpi[Gi];
                ncols = F.params->coltot[Gma];

                if (nrows && ncols && moinfo_.virtpi[Ga])
                    C_DGEMV('n', nrows, ncols, 1.0, T2.matrix[Gmi][T2.row_offset[Gmi][m]], ncols, F.matrix[Gma][0], 1,
                            1.0, &newtIA.matrix[Gi][0][A], moinfo_.virtpi[Ga]);
            }

            global_dpd_->buf4_mat_irrep_close(&T2, Gmi);
            global_dpd_->buf4_mat_irrep_row_close(&F, Gma);
        }
        global_dpd_->buf4_close(&F);
        global_dpd_->buf4_close(&T2);
        global_dpd_->file2_mat_wrt(&newtIA);
        global_dpd_->file2_mat_close(&newtIA);

        global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 11, 0, 11, 0, 0, "E 2<ai|jk> - <ai|kj>");
        global_dpd_->buf4_init(&tIjAb, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
        global_dpd_->contract442(&E, &tIjAb, &newtIA, 1, 3, -1, 1);
        global_dpd_->buf4_close(&E);
        global_dpd_->buf4_close(&tIjAb);

        if (params_.just_residuals) {
            global_dpd_->file2_close(&newtIA);
            return;
        }

        //    dpd_file2_copy(&newtIA, CC_OEI, "New tIA Increment");
        global_dpd_->file2_close(&newtIA);

        /*
            dpd_file2_init(&newtIA, CC_OEI, 0, 0, 1, "New tIA Increment");
            if(params.local && local.filter_singles) {
              local_filter_T1(&newtIA);
            }
            else {
              dpd_file2_init(&dIA, CC_OEI, 0, 0, 1, "dIA");
              dpd_file2_dirprd(&dIA, &newtIA);
              dpd_file2_close(&dIA);
            }
            dpd_file2_close(&newtIA);

            dpd_file2_init(&tIA, CC_OEI, 0, 0, 1, "tIA");
            dpd_file2_copy(&tIA, CC_OEI, "New tIA");
            dpd_file2_close(&tIA);
            dpd_file2_init(&newtIA, CC_OEI, 0, 0, 1, "New tIA");
            dpd_file2_init(&tIA, CC_OEI, 0, 0, 1, "New tIA Increment");
            dpd_file2_axpy(&tIA, &newtIA, 1, 0);
            dpd_file2_close(&tIA);
            dpd_file2_close(&newtIA);
        */
    } else if (params_.ref == 1) { /** ROHF **/

        global_dpd_->file2_init(&fIA, PSIF_CC_OEI, 0, 0, 1, "fIA");
        global_dpd_->file2_copy(&fIA, PSIF_CC_OEI, "New tIA");
        global_dpd_->file2_close(&fIA);

        global_dpd_->file2_init(&fia, PSIF_CC_OEI, 0, 0, 1, "fia");
        global_dpd_->file2_copy(&fia, PSIF_CC_OEI, "New tia");
        global_dpd_->file2_close(&fia);

        global_dpd_->file2_init(&newtIA, PSIF_CC_OEI, 0, 0, 1, "New tIA");
        global_dpd_->file2_init(&newtia, PSIF_CC_OEI, 0, 0, 1, "New tia");

        global_dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
        global_dpd_->file2_init(&tia, PSIF_CC_OEI, 0, 0, 1, "tia");

        global_dpd_->file2_init(&FAE, PSIF_CC_OEI, 0, 1, 1, "FAE");
        global_dpd_->file2_init(&Fae, PSIF_CC_OEI, 0, 1, 1, "Fae");

        global_dpd_->contract222(&tIA, &FAE, &newtIA, 0, 0, 1, 1);
        global_dpd_->contract222(&tia, &Fae, &newtia, 0, 0, 1, 1);

        global_dpd_->file2_close(&FAE);
        global_dpd_->file2_close(&Fae);

        global_dpd_->file2_close(&tIA);
        global_dpd_->file2_close(&tia);

        global_dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
        global_dpd_->file2_init(&tia, PSIF_CC_OEI, 0, 0, 1, "tia");

        global_dpd_->file2_init(&FMI, PSIF_CC_OEI, 0, 0, 0, "FMI");
        global_dpd_->file2_init(&Fmi, PSIF_CC_OEI, 0, 0, 0, "Fmi");

        global_dpd_->contract222(&FMI, &tIA, &newtIA, 1, 1, -1, 1);
        global_dpd_->contract222(&Fmi, &tia, &newtia, 1, 1, -1, 1);

        global_dpd_->file2_close(&FMI);
        global_dpd_->file2_close(&Fmi);
        global_dpd_->file2_close(&tIA);
        global_dpd_->file2_close(&tia);

        global_dpd_->file2_init(&FME, PSIF_CC_OEI, 0, 0, 1, "FME");
        global_dpd_->file2_init(&Fme, PSIF_CC_OEI, 0, 0, 1, "Fme");

        global_dpd_->buf4_init(&tIJAB, PSIF_CC_TAMPS, 0, 0, 5, 2, 7, 0, "tIJAB");
        global_dpd_->buf4_init(&tijab, PSIF_CC_TAMPS, 0, 0, 5, 2, 7, 0, "tijab");

        global_dpd_->dot13(&FME, &tIJAB, &newtIA, 0, 0, 1, 1);
        global_dpd_->dot13(&Fme, &tijab, &newtia, 0, 0, 1, 1);

        global_dpd_->buf4_close(&tIJAB);
        global_dpd_->buf4_close(&tijab);

        global_dpd_->buf4_init(&tIjAb, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");

        global_dpd_->dot24(&Fme, &tIjAb, &newtIA, 0, 0, 1, 1);
        global_dpd_->dot13(&FME, &tIjAb, &newtia, 0, 0, 1, 1);

        global_dpd_->buf4_close(&tIjAb);

        global_dpd_->file2_close(&FME);
        global_dpd_->file2_close(&Fme);

        global_dpd_->buf4_init(&C_anti, PSIF_CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia||jb>");
        global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");

        global_dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
        global_dpd_->file2_init(&tia, PSIF_CC_OEI, 0, 0, 1, "tia");

        global_dpd_->dot14(&tIA, &C_anti, &newtIA, 0, 1, -1, 1);
        global_dpd_->dot13(&tia, &D, &newtIA, 0, 0, 1, 1);

        global_dpd_->dot14(&tia, &C_anti, &newtia, 0, 1, -1, 1);
        global_dpd_->dot13(&tIA, &D, &newtia, 0, 0, 1, 1);

        global_dpd_->file2_close(&tIA);
        global_dpd_->file2_close(&tia);

        global_dpd_->buf4_close(&C_anti);
        global_dpd_->buf4_close(&D);

        global_dpd_->buf4_init(&F_anti, PSIF_CC_FINTS, 0, 10, 7, 10, 5, 1, "F <ia|bc>");
        global_dpd_->buf4_init(&tIJAB, PSIF_CC_TAMPS, 0, 0, 7, 2, 7, 0, "tIJAB");
        global_dpd_->buf4_init(&tijab, PSIF_CC_TAMPS, 0, 0, 7, 2, 7, 0, "tijab");

        global_dpd_->contract442(&tIJAB, &F_anti, &newtIA, 1, 1, 1, 1);
        global_dpd_->contract442(&tijab, &F_anti, &newtia, 1, 1, 1, 1);

        global_dpd_->buf4_close(&tIJAB);
        global_dpd_->buf4_close(&tijab);
        global_dpd_->buf4_close(&F_anti);

        global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc>");
        global_dpd_->buf4_init(&tIjAb, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
        global_dpd_->buf4_init(&tiJaB, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tiJaB");

        global_dpd_->contract442(&tiJaB, &F, &newtIA, 1, 1, 1, 1);
        global_dpd_->contract442(&tIjAb, &F, &newtia, 1, 1, 1, 1);

        global_dpd_->buf4_close(&F);
        global_dpd_->buf4_close(&tIjAb);
        global_dpd_->buf4_close(&tiJaB);

        global_dpd_->buf4_init(&E_anti, PSIF_CC_EINTS, 0, 11, 2, 11, 0, 1, "E <ai|jk>");
        global_dpd_->buf4_init(&tIJAB, PSIF_CC_TAMPS, 0, 2, 5, 2, 7, 0, "tIJAB");
        global_dpd_->buf4_init(&tijab, PSIF_CC_TAMPS, 0, 2, 5, 2, 7, 0, "tijab");

        global_dpd_->contract442(&E_anti, &tIJAB, &newtIA, 1, 3, -1, 1);
        global_dpd_->contract442(&E_anti, &tijab, &newtia, 1, 3, -1, 1);

        global_dpd_->buf4_close(&E_anti);
        global_dpd_->buf4_close(&tIJAB);
        global_dpd_->buf4_close(&tijab);

        global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 11, 0, 11, 0, 0, "E <ai|jk>");
        global_dpd_->buf4_init(&tIjAb, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
        global_dpd_->buf4_init(&tiJaB, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tiJaB");

        global_dpd_->contract442(&E, &tiJaB, &newtIA, 1, 3, -1, 1);
        global_dpd_->contract442(&E, &tIjAb, &newtia, 1, 3, -1, 1);

        global_dpd_->buf4_close(&E);
        global_dpd_->buf4_close(&tIjAb);
        global_dpd_->buf4_close(&tiJaB);

        /*
            dpd_file2_init(&dIA, CC_OEI, 0, 0, 1, "dIA");
            dpd_file2_dirprd(&dIA, &newtIA);
            dpd_file2_close(&dIA);

            dpd_file2_init(&dia, CC_OEI, 0, 0, 1, "dia");
            dpd_file2_dirprd(&dia, &newtia);
            dpd_file2_close(&dia);
        */

        global_dpd_->file2_close(&newtIA);
        global_dpd_->file2_close(&newtia);
    } else if (params_.ref == 2) { /*** UHF ***/

        global_dpd_->file2_init(&fIA, PSIF_CC_OEI, 0, 0, 1, "fIA");
        global_dpd_->file2_copy(&fIA, PSIF_CC_OEI, "New tIA");
        global_dpd_->file2_close(&fIA);

        global_dpd_->file2_init(&fia, PSIF_CC_OEI, 0, 2, 3, "fia");
        global_dpd_->file2_copy(&fia, PSIF_CC_OEI, "New tia");
        global_dpd_->file2_close(&fia);

        global_dpd_->file2_init(&newtIA, PSIF_CC_OEI, 0, 0, 1, "New tIA");
        global_dpd_->file2_init(&newtia, PSIF_CC_OEI, 0, 2, 3, "New tia");

        global_dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
        global_dpd_->file2_init(&tia, PSIF_CC_OEI, 0, 2, 3, "tia");

        global_dpd_->file2_init(&FAE, PSIF_CC_OEI, 0, 1, 1, "FAE");
        global_dpd_->contract222(&tIA, &FAE, &newtIA, 0, 0, 1, 1);
        global_dpd_->file2_close(&FAE);

        global_dpd_->file2_init(&Fae, PSIF_CC_OEI, 0, 3, 3, "Fae");
        global_dpd_->contract222(&tia, &Fae, &newtia, 0, 0, 1, 1);
        global_dpd_->file2_close(&Fae);

        global_dpd_->file2_init(&FMI, PSIF_CC_OEI, 0, 0, 0, "FMI");
        global_dpd_->contract222(&FMI, &tIA, &newtIA, 1, 1, -1, 1);
        global_dpd_->file2_close(&FMI);

        global_dpd_->file2_init(&Fmi, PSIF_CC_OEI, 0, 2, 2, "Fmi");
        global_dpd_->contract222(&Fmi, &tia, &newtia, 1, 1, -1, 1);
        global_dpd_->file2_close(&Fmi);

        global_dpd_->file2_close(&tIA);
        global_dpd_->file2_close(&tia);

        global_dpd_->file2_init(&FME, PSIF_CC_OEI, 0, 0, 1, "FME");
        global_dpd_->file2_init(&Fme, PSIF_CC_OEI, 0, 2, 3, "Fme");

        global_dpd_->buf4_init(&tIJAB, PSIF_CC_TAMPS, 0, 0, 5, 2, 7, 0, "tIJAB");
        global_dpd_->dot13(&FME, &tIJAB, &newtIA, 0, 0, 1, 1);
        global_dpd_->buf4_close(&tIJAB);

        global_dpd_->buf4_init(&tijab, PSIF_CC_TAMPS, 0, 10, 15, 12, 17, 0, "tijab");
        global_dpd_->dot13(&Fme, &tijab, &newtia, 0, 0, 1, 1);
        global_dpd_->buf4_close(&tijab);

        global_dpd_->buf4_init(&tIjAb, PSIF_CC_TAMPS, 0, 22, 28, 22, 28, 0, "tIjAb");
        global_dpd_->dot24(&Fme, &tIjAb, &newtIA, 0, 0, 1, 1);
        global_dpd_->dot13(&FME, &tIjAb, &newtia, 0, 0, 1, 1);
        global_dpd_->buf4_close(&tIjAb);

        global_dpd_->file2_close(&FME);
        global_dpd_->file2_close(&Fme);

        global_dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
        global_dpd_->file2_init(&tia, PSIF_CC_OEI, 0, 2, 3, "tia");

        global_dpd_->buf4_init(&C, PSIF_CC_CINTS, 0, 20, 20, 20, 20, 0, "C <IA||JB>");
        global_dpd_->dot14(&tIA, &C, &newtIA, 0, 1, -1, 1);
        global_dpd_->buf4_close(&C);

        global_dpd_->buf4_init(&C, PSIF_CC_CINTS, 0, 30, 30, 30, 30, 0, "C <ia||jb>");
        global_dpd_->dot14(&tia, &C, &newtia, 0, 1, -1, 1);
        global_dpd_->buf4_close(&C);

        global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 23, 29, 23, 29, 0, "D <iJ|aB>");
        global_dpd_->dot13(&tia, &D, &newtIA, 0, 0, 1, 1);
        global_dpd_->buf4_close(&D);

        global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 22, 28, 22, 28, 0, "D <Ij|Ab>");
        global_dpd_->dot13(&tIA, &D, &newtia, 0, 0, 1, 1);
        global_dpd_->buf4_close(&D);

        global_dpd_->file2_close(&tIA);
        global_dpd_->file2_close(&tia);

        global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 20, 7, 20, 5, 1, "F <IA|BC>");
        global_dpd_->buf4_init(&tIJAB, PSIF_CC_TAMPS, 0, 0, 7, 2, 7, 0, "tIJAB");
        global_dpd_->contract442(&tIJAB, &F, &newtIA, 1, 1, 1, 1);
        global_dpd_->buf4_close(&tIJAB);
        global_dpd_->buf4_close(&F);

        global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 30, 17, 30, 15, 1, "F <ia|bc>");
        global_dpd_->buf4_init(&tijab, PSIF_CC_TAMPS, 0, 10, 17, 12, 17, 0, "tijab");
        global_dpd_->contract442(&tijab, &F, &newtia, 1, 1, 1, 1);
        global_dpd_->buf4_close(&tijab);
        global_dpd_->buf4_close(&F);

        global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 28, 26, 28, 26, 0, "F <Ab|Ci>");
        global_dpd_->buf4_init(&tIjAb, PSIF_CC_TAMPS, 0, 22, 28, 22, 28, 0, "tIjAb");
        global_dpd_->contract442(&tIjAb, &F, &newtIA, 0, 2, 1, 1);
        global_dpd_->buf4_close(&tIjAb);
        global_dpd_->buf4_close(&F);

        global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 24, 28, 24, 28, 0, "F <Ia|Bc>");
        global_dpd_->buf4_init(&tIjAb, PSIF_CC_TAMPS, 0, 22, 28, 22, 28, 0, "tIjAb");
        global_dpd_->contract442(&tIjAb, &F, &newtia, 1, 1, 1, 1);
        global_dpd_->buf4_close(&tIjAb);
        global_dpd_->buf4_close(&F);

        global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 21, 2, 21, 0, 1, "E <AI|JK>");
        global_dpd_->buf4_init(&tIJAB, PSIF_CC_TAMPS, 0, 2, 5, 2, 7, 0, "tIJAB");
        global_dpd_->contract442(&E, &tIJAB, &newtIA, 1, 3, -1, 1);
        global_dpd_->buf4_close(&E);

        global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 31, 12, 31, 10, 1, "E <ai|jk>");
        global_dpd_->buf4_init(&tijab, PSIF_CC_TAMPS, 0, 12, 15, 12, 17, 0, "tijab");
        global_dpd_->contract442(&E, &tijab, &newtia, 1, 3, -1, 1);
        global_dpd_->buf4_close(&E);

        global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 22, 24, 22, 24, 0, "E <Ij|Ka>");
        global_dpd_->buf4_init(&tIjAb, PSIF_CC_TAMPS, 0, 22, 28, 22, 28, 0, "tIjAb");
        global_dpd_->contract442(&E, &tIjAb, &newtIA, 2, 2, -1, 1);
        global_dpd_->buf4_close(&E);
        global_dpd_->buf4_close(&tIjAb);

        global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 23, 27, 23, 27, 0, "E <iJ|kA>");
        global_dpd_->buf4_init(&tIjAb, PSIF_CC_TAMPS, 0, 23, 29, 23, 29, 0, "tiJaB");
        global_dpd_->contract442(&E, &tIjAb, &newtia, 2, 2, -1, 1);
        global_dpd_->buf4_close(&E);
        global_dpd_->buf4_close(&tIjAb);

        /*
            dpd_file2_init(&dIA, CC_OEI, 0, 0, 1, "dIA");
            dpd_file2_dirprd(&dIA, &newtIA);
            dpd_file2_close(&dIA);

            dpd_file2_init(&dia, CC_OEI, 0, 2, 3, "dia");
            dpd_file2_dirprd(&dia, &newtia);
            dpd_file2_close(&dia);
        */

        global_dpd_->file2_close(&newtIA);
        global_dpd_->file2_close(&newtia);
    }
}

// Mixed-precision
void CCEnergyWavefunction::t1_build_mp() {
    dpdfile2<double> newtIA, tIA, fIA, file2_tmp_double;
    dpdfile2<float> tIA_sp, fIA_sp, file2_tmp_float;
    dpdfile2<double> FAE, FMI, FME;
    dpdfile2<float> FAE_sp, FMI_sp, FME_sp;
    dpdbuf4<double> tIjAb, T2;
    dpdbuf4<float> tIjAb_sp, T2_sp;
    dpdbuf4<double> C_anti, D, F, E;
    dpdbuf4<float> C_anti_sp, D_sp, F_sp, E_sp;
    int Gmi, Gm, Gi, Ga, m, a, A, nrows, ncols;
    int row, col;
    float **TMP;

    if (params_.ref == 0) { /** RHF **/
        // **1
        global_dpd_->file2_init(&fIA, PSIF_CC_OEI, 0, 0, 1, "fIA");
        global_dpd_->file2_copy(&fIA, PSIF_CC_OEI, "New tIA");
        global_dpd_->file2_close(&fIA);
        
        // **2
        global_dpd_->file2_init(&newtIA, PSIF_CC_OEI, 0, 0, 1, "New tIA");

        global_dpd_->file2_init_sp(&tIA_sp, PSIF_CC_OEI, 0, 0, 1, "tIA_sp");
        
        global_dpd_->file2_init_sp(&FAE_sp, PSIF_CC_OEI, 0, 1, 1, "FAE_sp");
        // Contractions of tIA and FAE in single-precision (the second term in t1 eqn ) -> cast to double-precision -> add to the residual 
        global_dpd_->contract222_mp(&tIA_sp, &FAE_sp, &newtIA, 0, 0, 1, 1);
        global_dpd_->file2_close_sp(&FAE_sp);
        
        // **3
        global_dpd_->file2_init_sp(&FMI_sp, PSIF_CC_OEI, 0, 0, 0, "FMI_sp");
        global_dpd_->contract222_mp(&FMI_sp, &tIA_sp, &newtIA, 1, 1, -1, 1);
        global_dpd_->file2_close_sp(&FMI_sp);
            
        global_dpd_->file2_close_sp(&tIA_sp);

        // **4
        global_dpd_->file2_init_sp(&FME_sp, PSIF_CC_OEI, 0, 0, 1, "FME_sp");

        global_dpd_->buf4_init_sp(&tIjAb_sp, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "2 tIAjb - tIBja sp");
        global_dpd_->contract422_mp(&tIjAb_sp, &FME_sp, &newtIA, 0, 0, 1, 1);
        global_dpd_->buf4_close_sp(&tIjAb_sp);

        global_dpd_->file2_close_sp(&FME_sp);
               
        // **5
        global_dpd_->buf4_init_sp(&C_anti_sp, PSIF_CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia||jb> sp");
        global_dpd_->buf4_init_sp(&D_sp, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab> sp");

        global_dpd_->file2_init_sp(&tIA_sp, PSIF_CC_OEI, 0, 0, 1, "tIA_sp");
        global_dpd_->dot14_mp(&tIA_sp, &C_anti_sp, &newtIA, 0, 1, -1, 1);
        global_dpd_->dot13_mp(&tIA_sp, &D_sp, &newtIA, 0, 0, 1, 1);
        
        global_dpd_->file2_close_sp(&tIA_sp);

        global_dpd_->buf4_close_sp(&C_anti_sp);
        global_dpd_->buf4_close_sp(&D_sp);
        
        // **6
        /*
          dpd_buf4_init(&Z, CC_TMP0, 0, 10, 0, 10, 0, 0, "Z(ma,mi)");
          dpd_buf4_init(&F, CC_FINTS, 0, 10, 5, 10, 5, 0, "F 2<ia|bc> - <ia|cb>");
          dpd_buf4_init(&tIjAb, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
          dpd_contract444(&F, &tIjAb, &Z, 0, 0, 1.0, 0.0);
          dpd_buf4_close(&tIjAb);
          dpd_buf4_close(&F);
          dpd_trace42_13(&Z, &newtIA, 1, 1.0, 1.0);
          dpd_buf4_close(&Z);
        */

        /* t(i,a) <-- (2 t(mi,ef) - t(mi,fe)) <ma|ef> */
        /* out-of-core version replacing the *stupid* code above 3/22/05, TDC */
        
        global_dpd_->file2_mat_init(&newtIA);
        global_dpd_->file2_mat_rd(&newtIA);
        //global_dpd_->file2_init_sp(&TMP, PSIF_CC_TMP0, 0, 0, 1, "TMP newtIA");
        //global_dpd_->file2_mat_init_sp(&TMP);
        //global_dpd_->file2_mat_rd_sp(&TMP);
        global_dpd_->buf4_init_sp(&T2_sp, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "2 tIjAb - tIjBa sp");
        global_dpd_->buf4_init_sp(&F_sp, PSIF_CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc> sp");
        for (int Gma = 0; Gma < moinfo_.nirreps; Gma++) {
            Gmi = Gma; /* T1 is totally symmetric */

            global_dpd_->buf4_mat_irrep_row_init_sp(&F_sp, Gma);
            global_dpd_->buf4_mat_irrep_init_sp(&T2_sp, Gmi);
            global_dpd_->buf4_mat_irrep_rd_sp(&T2_sp, Gmi);

            for (int ma = 0; ma < F.params->rowtot[Gma]; ma++) {
                global_dpd_->buf4_mat_irrep_row_rd_sp(&F_sp, Gma, ma);

                m = F_sp.params->roworb[Gma][ma][0];
                a = F_sp.params->roworb[Gma][ma][1];
                Gm = F_sp.params->psym[m];
                Ga = F_sp.params->qsym[a];
                Gi = Ga; /* T1 is totally symmetric */
                A = a - F_sp.params->qoff[Ga];

                nrows = moinfo_.occpi[Gi];
                ncols = F_sp.params->coltot[Gma];

                if (nrows && ncols && moinfo_.virtpi[Ga]){
                    C_SGEMV('n', nrows, ncols, 1.0, T2_sp.matrix[Gmi][T2_sp.row_offset[Gmi][m]], ncols, F_sp.matrix[Gma][0], 1,
                            0.0, &(TMP[0][0]), moinfo_.virtpi[Ga]);
                    for (row = 0; row < moinfo_.virtpi[Ga]; row++){
                         newtIA.matrix[Gi][0][A+row] += static_cast<double>(TMP[0][row]);   
                    }
               }
            }

            global_dpd_->buf4_mat_irrep_close_sp(&T2_sp, Gmi);
            global_dpd_->buf4_mat_irrep_row_close_sp(&F_sp, Gma);
        }
        global_dpd_->buf4_close_sp(&F_sp);
        global_dpd_->buf4_close_sp(&T2_sp);
        //global_dpd_->file2_close_sp(&TMP);
        free(TMP);
        global_dpd_->file2_mat_wrt(&newtIA);
        global_dpd_->file2_mat_close(&newtIA);

        global_dpd_->buf4_init_sp(&E_sp, PSIF_CC_EINTS, 0, 11, 0, 11, 0, 0, "E 2<ai|jk> - <ai|kj> sp");
        global_dpd_->buf4_init_sp(&tIjAb_sp, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb_sp");
        global_dpd_->contract442_mp(&E_sp, &tIjAb_sp, &newtIA, 1, 3, -1, 1);
        global_dpd_->buf4_close_sp(&E_sp);
        global_dpd_->buf4_close_sp(&tIjAb_sp);
        
               
        if (params_.just_residuals) {
            global_dpd_->file2_close(&newtIA);
            return;
        }

        //    dpd_file2_copy(&newtIA, CC_OEI, "New tIA Increment");
        global_dpd_->file2_close(&newtIA);

        /*
            dpd_file2_init(&newtIA, CC_OEI, 0, 0, 1, "New tIA Increment");
            if(params.local && local.filter_singles) {
              local_filter_T1(&newtIA);
            }
            else {
              dpd_file2_init(&dIA, CC_OEI, 0, 0, 1, "dIA");
              dpd_file2_dirprd(&dIA, &newtIA);
              dpd_file2_close(&dIA);
            }
            dpd_file2_close(&newtIA);

            dpd_file2_init(&tIA, CC_OEI, 0, 0, 1, "tIA");
            dpd_file2_copy(&tIA, CC_OEI, "New tIA");
            dpd_file2_close(&tIA);
            dpd_file2_init(&newtIA, CC_OEI, 0, 0, 1, "New tIA");
            dpd_file2_init(&tIA, CC_OEI, 0, 0, 1, "New tIA Increment");
            dpd_file2_axpy(&tIA, &newtIA, 1, 0);
            dpd_file2_close(&tIA);
            dpd_file2_close(&newtIA);
        */
    } // RHF
} //build_t1_mp


// Single-precision
void CCEnergyWavefunction::t1_build_sp() {

    dpdfile2<float> newtIA_sp, tIA_sp, fIA_sp;
    dpdfile2<float> FAE_sp, FMI_sp, FME_sp;
    dpdbuf4<float> tIjAb_sp, T2_sp;
    dpdbuf4<float> C_anti_sp, D_sp, F_sp, E_sp;
    int Gmi, Gm, Gi, Ga, m, a, A, nrows, ncols;
   
    if (params_.ref == 0) { /** RHF **/
        global_dpd_->file2_init_sp(&fIA_sp, PSIF_CC_OEI, 0, 0, 1, "fIA_sp");
        global_dpd_->file2_copy_sp(&fIA_sp, PSIF_CC_OEI, "New tIA sp");
        global_dpd_->file2_close_sp(&fIA_sp);

        global_dpd_->file2_init_sp(&newtIA_sp, PSIF_CC_OEI, 0, 0, 1, "New tIA sp");

        global_dpd_->file2_init_sp(&tIA_sp, PSIF_CC_OEI, 0, 0, 1, "tIA_sp");

        global_dpd_->file2_init_sp(&FAE_sp, PSIF_CC_OEI, 0, 1, 1, "FAE_sp");
        global_dpd_->contract222_sp(&tIA_sp, &FAE_sp, &newtIA_sp, 0, 0, 1, 1);
        global_dpd_->file2_close_sp(&FAE_sp);

        global_dpd_->file2_init_sp(&FMI_sp, PSIF_CC_OEI, 0, 0, 0, "FMI_sp");
        global_dpd_->contract222_sp(&FMI_sp, &tIA_sp, &newtIA_sp, 1, 1, -1, 1);
        global_dpd_->file2_close_sp(&FMI_sp);

        global_dpd_->file2_close_sp(&tIA_sp);

        global_dpd_->file2_init_sp(&FME_sp, PSIF_CC_OEI, 0, 0, 1, "FME_sp");

        global_dpd_->buf4_init_sp(&tIjAb_sp, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "2 tIAjb - tIBja sp");
        global_dpd_->contract422_sp(&tIjAb_sp, &FME_sp, &newtIA_sp, 0, 0, 1, 1);
        global_dpd_->buf4_close_sp(&tIjAb_sp);

        global_dpd_->file2_close_sp(&FME_sp);

        global_dpd_->buf4_init_sp(&C_anti_sp, PSIF_CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia||jb> sp");
        global_dpd_->buf4_init_sp(&D_sp, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab> sp");

        global_dpd_->file2_init_sp(&tIA_sp, PSIF_CC_OEI, 0, 0, 1, "tIA_sp");

        global_dpd_->dot14_sp(&tIA_sp, &C_anti_sp, &newtIA_sp, 0, 1, -1, 1);
        global_dpd_->dot13_sp(&tIA_sp, &D_sp, &newtIA_sp, 0, 0, 1, 1);

        global_dpd_->file2_close_sp(&tIA_sp);

        global_dpd_->buf4_close_sp(&C_anti_sp);
        global_dpd_->buf4_close_sp(&D_sp);

        /*
          dpd_buf4_init(&Z, CC_TMP0, 0, 10, 0, 10, 0, 0, "Z(ma,mi)");
          dpd_buf4_init(&F, CC_FINTS, 0, 10, 5, 10, 5, 0, "F 2<ia|bc> - <ia|cb>");
          dpd_buf4_init(&tIjAb, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
          dpd_contract444(&F, &tIjAb, &Z, 0, 0, 1.0, 0.0);
          dpd_buf4_close(&tIjAb);
          dpd_buf4_close(&F);
          dpd_trace42_13(&Z, &newtIA, 1, 1.0, 1.0);
          dpd_buf4_close(&Z);
        */

        /* t(i,a) <-- (2 t(mi,ef) - t(mi,fe)) <ma|ef> */
        /* out-of-core version replacing the *stupid* code above 3/22/05, TDC */
        global_dpd_->file2_mat_init_sp(&newtIA_sp);
        global_dpd_->file2_mat_rd_sp(&newtIA_sp);
        global_dpd_->buf4_init_sp(&T2_sp, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "2 tIjAb - tIjBa sp");
        global_dpd_->buf4_init_sp(&F_sp, PSIF_CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc> sp");
        for (int Gma = 0; Gma < moinfo_.nirreps; Gma++) {
            Gmi = Gma; /* T1 is totally symmetric */

            global_dpd_->buf4_mat_irrep_row_init_sp(&F_sp, Gma);
            global_dpd_->buf4_mat_irrep_init_sp(&T2_sp, Gmi);
            global_dpd_->buf4_mat_irrep_rd_sp(&T2_sp, Gmi);

            for (int ma = 0; ma < F_sp.params->rowtot[Gma]; ma++) {
                global_dpd_->buf4_mat_irrep_row_rd_sp(&F_sp, Gma, ma);

                m = F_sp.params->roworb[Gma][ma][0];
                a = F_sp.params->roworb[Gma][ma][1];
                Gm = F_sp.params->psym[m];
                Ga = F_sp.params->qsym[a];
                Gi = Ga; /* T1 is totally symmetric */
                A = a - F_sp.params->qoff[Ga];

                nrows = moinfo_.occpi[Gi];
                ncols = F_sp.params->coltot[Gma];

                if (nrows && ncols && moinfo_.virtpi[Ga])
                    C_SGEMV('n', nrows, ncols, 1.0, T2_sp.matrix[Gmi][T2_sp.row_offset[Gmi][m]], ncols, F_sp.matrix[Gma][0], 1,

                            1.0, &newtIA_sp.matrix[Gi][0][A], moinfo_.virtpi[Ga]);
            }

            global_dpd_->buf4_mat_irrep_close_sp(&T2_sp, Gmi);
            global_dpd_->buf4_mat_irrep_row_close_sp(&F_sp, Gma);
        }
        global_dpd_->buf4_close_sp(&F_sp);
        global_dpd_->buf4_close_sp(&T2_sp);
        global_dpd_->file2_mat_wrt_sp(&newtIA_sp);
        global_dpd_->file2_mat_close_sp(&newtIA_sp);

        global_dpd_->buf4_init_sp(&E_sp, PSIF_CC_EINTS, 0, 11, 0, 11, 0, 0, "E 2<ai|jk> - <ai|kj> sp");
        global_dpd_->buf4_init_sp(&tIjAb_sp, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb_sp");
        global_dpd_->contract442_sp(&E_sp, &tIjAb_sp, &newtIA_sp, 1, 3, -1, 1);
        global_dpd_->buf4_close_sp(&E_sp);
        global_dpd_->buf4_close_sp(&tIjAb_sp);

        if (params_.just_residuals) {
            global_dpd_->file2_close_sp(&newtIA_sp);
            return;
        }

        //    dpd_file2_copy(&newtIA, CC_OEI, "New tIA Increment");
        global_dpd_->file2_close_sp(&newtIA_sp);

        /*
            dpd_file2_init(&newtIA, CC_OEI, 0, 0, 1, "New tIA Increment");
            if(params.local && local.filter_singles) {
              local_filter_T1(&newtIA);
            }
            else {
              dpd_file2_init(&dIA, CC_OEI, 0, 0, 1, "dIA");
              dpd_file2_dirprd(&dIA, &newtIA);
              dpd_file2_close(&dIA);
            }
            dpd_file2_close(&newtIA);

            dpd_file2_init(&tIA, CC_OEI, 0, 0, 1, "tIA");
            dpd_file2_copy(&tIA, CC_OEI, "New tIA");
            dpd_file2_close(&tIA);
            dpd_file2_init(&newtIA, CC_OEI, 0, 0, 1, "New tIA");
            dpd_file2_init(&tIA, CC_OEI, 0, 0, 1, "New tIA Increment");
            dpd_file2_axpy(&tIA, &newtIA, 1, 0);
            dpd_file2_close(&tIA);
            dpd_file2_close(&newtIA);
        */
    }
    }  // build_t1_sp

}  // namespace ccenergy
}  // namespace psi
