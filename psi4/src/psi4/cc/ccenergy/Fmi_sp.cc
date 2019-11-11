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
#include "Params.h"
#include "MOInfo.h"
#include "psi4/cc/ccwave.h"

namespace psi {
namespace ccenergy {

void CCEnergyWavefunction::Fmi_build_sp() {
    dpdfile2<float> FMI, FMIt, fIJ, fIA;
    dpdfile2<float> tIA, FME;
    dpdbuf4<float> E, E_anti, D;
    dpdbuf4<float> tautIjAb;

    if (params_.ref == 0) { /** RHF **/
        global_dpd_->file2_init_sp(&fIJ, PSIF_CC_OEI, 0, 0, 0, "fIJ_sp");
        global_dpd_->file2_copy_sp(&fIJ, PSIF_CC_OEI, "FMI_sp");
        global_dpd_->file2_close_sp(&fIJ);
    } 
    if (params_.ref == 0) { /** RHF **/
        global_dpd_->file2_init_sp(&FMI, PSIF_CC_OEI, 0, 0, 0, "FMI_sp");
        global_dpd_->file2_mat_init_sp(&FMI);
        global_dpd_->file2_mat_rd_sp(&FMI);

        /*
        for(int h = 0; h < moinfo.nirreps; h++) {
          for(int m = 0; m < FMI.params->rowtot[h]; m++)
            FMI.matrix[h][m][m] = 0;
        }
        */

        global_dpd_->file2_mat_wrt_sp(&FMI);
        global_dpd_->file2_mat_close_sp(&FMI);
        global_dpd_->file2_close_sp(&FMI);
    } 
    if (params_.ref == 0) { /** RHF **/
        global_dpd_->file2_init_sp(&FMI, PSIF_CC_OEI, 0, 0, 0, "FMI_sp");

        global_dpd_->file2_init_sp(&fIA, PSIF_CC_OEI, 0, 0, 1, "fIA_sp");
        global_dpd_->file2_init_sp(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA_sp");
        global_dpd_->contract222_sp(&fIA, &tIA, &FMI, 0, 0, 0.5, 1);
        global_dpd_->file2_close_sp(&tIA);
        global_dpd_->file2_close_sp(&fIA);

        global_dpd_->buf4_init_sp(&E_anti, PSIF_CC_EINTS, 0, 11, 0, 11, 0, 1, "E <ai|jk> sp");
        global_dpd_->buf4_init_sp(&E, PSIF_CC_EINTS, 0, 11, 0, 11, 0, 0, "E <ai|jk> sp");

        global_dpd_->file2_init_sp(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA_sp");

        global_dpd_->dot13_sp(&tIA, &E_anti, &FMI, 1, 1, 1.0, 1.0);
        global_dpd_->dot13_sp(&tIA, &E, &FMI, 1, 1, 1.0, 1.0);

        global_dpd_->file2_close_sp(&tIA);

        global_dpd_->buf4_close_sp(&E_anti);
        global_dpd_->buf4_close_sp(&E);

        global_dpd_->buf4_init_sp(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D 2<ij|ab> - <ij|ba> sp");
        global_dpd_->buf4_init_sp(&tautIjAb, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tautIjAb_sp");
        global_dpd_->contract442_sp(&D, &tautIjAb, &FMI, 0, 0, 1, 1);
        global_dpd_->buf4_close_sp(&tautIjAb);
        global_dpd_->buf4_close_sp(&D);

        /* Build the tilde intermediate */
        global_dpd_->file2_copy_sp(&FMI, PSIF_CC_OEI, "FMIt_sp");
        global_dpd_->file2_close_sp(&FMI);

        global_dpd_->file2_init_sp(&FMIt, PSIF_CC_OEI, 0, 0, 0, "FMIt_sp");

        global_dpd_->file2_init_sp(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA_sp");
        global_dpd_->file2_init_sp(&FME, PSIF_CC_OEI, 0, 0, 1, "FME_sp");
        global_dpd_->contract222_sp(&FME, &tIA, &FMIt, 0, 0, 0.5, 1);
        global_dpd_->file2_close_sp(&FME);
        global_dpd_->file2_close_sp(&tIA);

        global_dpd_->file2_close_sp(&FMIt);
    } 
}
}  // namespace ccenergy
}  // namespace psi
