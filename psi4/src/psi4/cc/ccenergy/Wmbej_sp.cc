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
#include "psi4/libpsi4util/PsiOutStream.h"

namespace psi {
namespace ccenergy {

/* Wmbej_build(): Build the Wmbej intermediate.
**
** Wmbej = <mb||ej> + t(j,f) <mb||ef> - t(n,b) <mn||ej>
**         - [1/2 t(jn,fb) + t(j,f) t(n,b)] <mn||ef>
**
** Spin cases for UHF and ROHF orbitals:
** -------------------------------------
**
**
** TDC
** May 2000
*/

void CCEnergyWavefunction::Wmbej_build_sp() {
    dpdbuf4<float> WMbEj, WMbeJ, W;
    dpdbuf4<float> C, D, E, F, t2, Y;
    dpdfile2<float> tIA;
    int Ge, Gf, nrows, ncols, nlinks;

    timer_on("C->Wmbej");

    /* W(mb,je) <-- <mb||ej> */

    if (params_.ref == 0) { /** RHF **/
        global_dpd_->buf4_init_sp(&C, PSIF_CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia|jb> sp");
        global_dpd_->buf4_scmcopy_sp(&C, PSIF_CC_TMP0, "WMbeJ_sp", -1);
        global_dpd_->buf4_close_sp(&C);

        global_dpd_->buf4_init_sp(&D, PSIF_CC_DINTS, 0, 10, 11, 10, 11, 0, "D <ij|ab> (ib,aj) sp");
        global_dpd_->buf4_copy_sp(&D, PSIF_CC_TMP0, "WMbEj_sp");
        global_dpd_->buf4_close_sp(&D);
    } 
   
    timer_off("C->Wmbej");

    timer_on("F->Wmbej");

    if (params_.ref == 0) { /** RHF **/
        global_dpd_->file2_init_sp(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA_sp");

        global_dpd_->buf4_init_sp(&F, PSIF_CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc> sp");
        global_dpd_->buf4_init_sp(&WMbEj, PSIF_CC_TMP0, 0, 10, 11, 10, 11, 0, "WMbEj_sp");
        global_dpd_->contract424_sp(&F, &tIA, &WMbEj, 3, 1, 0, 1, 1); /* should be OOC-capable in libdpd */
        global_dpd_->buf4_close_sp(&WMbEj);
        global_dpd_->buf4_close_sp(&F);

        /*
        dpd_buf4_init(&F, CC_FINTS, 0, 11, 5, 11, 5, 0, "F <ai|bc>");

        dpd_buf4_init(&Z, CC_TMP0, 0, 11, 11, 11, 11, 0, "Z(bM,eJ)");
        dpd_contract424(&F, &tIA, &Z, 3, 1, 0, -1, 0);
        dpd_buf4_sort(&Z, CC_TMP0, qpsr, 10, 10, "Z(Mb,Je)");
        dpd_buf4_close(&Z);
        dpd_buf4_init(&Z, CC_TMP0, 0, 10, 10, 10, 10, 0, "Z(Mb,Je)");
        dpd_buf4_init(&WMbeJ, CC_TMP0, 0, 10, 10, 10, 10, 0, "WMbeJ");
        dpd_buf4_axpy(&Z, &WMbeJ, 1.0);
        dpd_buf4_close(&WMbeJ);
        dpd_buf4_close(&Z);

        dpd_buf4_close(&F);
        */

        /* W(Mb,Je) <-- t(J,F) <Mb|Fe> */
        /* OOC code added to replace above on 3/23/05, TDC */
        global_dpd_->buf4_init_sp(&F, PSIF_CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc> sp");
        global_dpd_->buf4_init_sp(&W, PSIF_CC_TMP0, 0, 10, 10, 10, 10, 0, "WMbeJ_sp");
        global_dpd_->file2_mat_init_sp(&tIA);
        global_dpd_->file2_mat_rd_sp(&tIA);

        for (int Gmb = 0; Gmb < moinfo_.nirreps; Gmb++) {
            global_dpd_->buf4_mat_irrep_row_init_sp(&W, Gmb);
            global_dpd_->buf4_mat_irrep_row_init_sp(&F, Gmb);

            for (int mb = 0; mb < F.params->rowtot[Gmb]; mb++) {
                global_dpd_->buf4_mat_irrep_row_rd_sp(&W, Gmb, mb);
                global_dpd_->buf4_mat_irrep_row_rd_sp(&F, Gmb, mb);

                for (int Gj = 0; Gj < moinfo_.nirreps; Gj++) {
                    Gf = Gj;       /* T1 is totally symmetric */
                    Ge = Gmb ^ Gf; /* <mb|fe> is totally symmetric */

                    nrows = moinfo_.occpi[Gj];
                    ncols = moinfo_.virtpi[Ge];
                    nlinks = moinfo_.virtpi[Gf];
                    if (nrows && ncols && nlinks)
                        C_SGEMM('n', 'n', nrows, ncols, nlinks, -1.0, tIA.matrix[Gj][0], nlinks,
                                &F.matrix[Gmb][0][F.col_offset[Gmb][Gf]], ncols, 1.0,
                                &W.matrix[Gmb][0][W.col_offset[Gmb][Gj]], ncols);
                }

                global_dpd_->buf4_mat_irrep_row_wrt_sp(&W, Gmb, mb);
            }

            global_dpd_->buf4_mat_irrep_row_close_sp(&F, Gmb);
            global_dpd_->buf4_mat_irrep_row_close_sp(&W, Gmb);
        }

        global_dpd_->file2_mat_close_sp(&tIA);
        global_dpd_->buf4_close_sp(&W);
        global_dpd_->buf4_close_sp(&F);

        global_dpd_->file2_close_sp(&tIA);
    }  
  
    timer_off("F->Wmbej");
     timer_on("E->Wmbej");

    if (params_.ref == 0) { /** RHF **/
        global_dpd_->file2_init_sp(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA_sp");

        global_dpd_->buf4_init_sp(&E, PSIF_CC_EINTS, 0, 11, 0, 11, 0, 0, "E <ai|jk> sp");
        global_dpd_->buf4_init_sp(&WMbEj, PSIF_CC_TMP0, 0, 10, 11, 10, 11, 0, "WMbEj_sp");
        global_dpd_->contract424_sp(&E, &tIA, &WMbEj, 3, 0, 1, -1, 1);
        global_dpd_->buf4_close_sp(&WMbEj);
        global_dpd_->buf4_close_sp(&E);

        global_dpd_->buf4_init_sp(&E, PSIF_CC_EINTS, 0, 0, 10, 0, 10, 0, "E <ij|ka> sp");
        global_dpd_->buf4_init_sp(&WMbeJ, PSIF_CC_TMP0, 0, 10, 10, 10, 10, 0, "WMbeJ_sp");
        global_dpd_->contract424_sp(&E, &tIA, &WMbeJ, 1, 0, 1, 1, 1);
        global_dpd_->buf4_close_sp(&WMbeJ);
        global_dpd_->buf4_close_sp(&E);

        global_dpd_->file2_close_sp(&tIA);
      }  
    timer_off("E->Wmbej");

    /* Convert to (ME,JB) for remaining terms */

    if (params_.ref == 0) { /** RHF **/

        global_dpd_->buf4_init_sp(&WMbEj, PSIF_CC_TMP0, 0, 10, 11, 10, 11, 0, "WMbEj_sp");
        global_dpd_->buf4_sort_sp(&WMbEj, PSIF_CC_HBAR, prsq, 10, 10, "WMbEj_sp");
        global_dpd_->buf4_close_sp(&WMbEj);

        global_dpd_->buf4_init_sp(&WMbeJ, PSIF_CC_TMP0, 0, 10, 10, 10, 10, 0, "WMbeJ_sp");
        global_dpd_->buf4_sort_sp(&WMbeJ, PSIF_CC_HBAR, psrq, 10, 10, "WMbeJ_sp");
        global_dpd_->buf4_close_sp(&WMbeJ);

    }
     
    // Check WMbEj
    outfile->Printf("Check WMbEj_sp: E(sort)\n");
    global_dpd_->buf4_init_sp(&WMbEj, PSIF_CC_HBAR, 0, 10, 10, 10, 10, 0, "WMbEj_sp");
    global_dpd_->buf4_print_sp(&WMbEj, "outfile", 1);
    global_dpd_->buf4_close_sp(&WMbEj);

 
    timer_on("X->Wmbej");

    if (params_.ref == 0) { /** RHF **/
        global_dpd_->file2_init_sp(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA_sp");

        /*** ABAB ***/

        global_dpd_->buf4_init_sp(&W, PSIF_CC_HBAR, 0, 10, 10, 10, 10, 0, "WMbEj_sp");
        global_dpd_->buf4_init_sp(&t2, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIAjb_sp");
        global_dpd_->buf4_init_sp(&D, PSIF_CC_DINTS, 0, 10, 10, 10, 10, 0, "D 2<ij|ab> - <ij|ba> (ia,jb) sp");
        global_dpd_->contract444_sp(&D, &t2, &W, 0, 1, 0.5, 1);
        global_dpd_->buf4_close_sp(&D);
        global_dpd_->buf4_close_sp(&t2);
        global_dpd_->buf4_init_sp(&t2, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "tjAIb_sp");
        global_dpd_->buf4_init_sp(&D, PSIF_CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij|ab> (ia,jb) sp");
        global_dpd_->contract444_sp(&D, &t2, &W, 0, 1, -0.5, 1);
        global_dpd_->buf4_close_sp(&D);
        global_dpd_->buf4_close_sp(&t2);
        global_dpd_->buf4_close_sp(&W);

        global_dpd_->buf4_init_sp(&Y, PSIF_CC_TMP0, 0, 10, 0, 10, 0, 0, "Y (ME,JN) sp");
        global_dpd_->buf4_init_sp(&D, PSIF_CC_DINTS, 0, 10, 11, 10, 11, 0, "D <ij|ab> (ia,bj) sp");
        global_dpd_->contract244_sp(&tIA, &D, &Y, 1, 2, 1, 1, 0);
        global_dpd_->buf4_close_sp(&D);
        global_dpd_->buf4_init_sp(&W, PSIF_CC_HBAR, 0, 10, 10, 10, 10, 0, "WMbEj_sp");
        global_dpd_->contract424_sp(&Y, &tIA, &W, 3, 0, 0, -1, 1);
        global_dpd_->buf4_close_sp(&W);
        global_dpd_->buf4_close_sp(&Y);

        /*** ABBA ***/

        global_dpd_->buf4_init_sp(&W, PSIF_CC_HBAR, 0, 10, 10, 10, 10, 0, "WMbeJ_sp");
        global_dpd_->buf4_init_sp(&t2, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIbjA_sp");
        global_dpd_->buf4_init_sp(&D, PSIF_CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij|ab> (ib,ja) sp");
        global_dpd_->contract444_sp(&D, &t2, &W, 0, 1, 0.5, 1);
        global_dpd_->buf4_close_sp(&D);
        global_dpd_->buf4_close_sp(&t2);
        global_dpd_->buf4_close_sp(&W);

        global_dpd_->buf4_init_sp(&Y, PSIF_CC_TMP0, 0, 10, 0, 10, 0, 0, "Y (ME,JN) sp");
        global_dpd_->buf4_init_sp(&D, PSIF_CC_DINTS, 0, 10, 11, 10, 11, 0, "D <ij|ab> (ib,aj) sp");
        global_dpd_->contract244_sp(&tIA, &D, &Y, 1, 2, 1, 1, 0);
        global_dpd_->buf4_close_sp(&D);
        global_dpd_->buf4_init_sp(&W, PSIF_CC_HBAR, 0, 10, 10, 10, 10, 0, "WMbeJ_sp");
        global_dpd_->contract424_sp(&Y, &tIA, &W, 3, 0, 0, 1, 1);
        global_dpd_->buf4_close_sp(&W);
        global_dpd_->buf4_close_sp(&Y);

        global_dpd_->file2_close_sp(&tIA);
    } 
    timer_off("X->Wmbej");
}
}  // namespace ccenergy
}  // namespace psi
