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
#include "psi4/libciomr/libciomr.h"
#include "psi4/libdpd/dpd.h"
#include "psi4/libqt/qt.h"
#include "Params.h"
#include "psi4/cc/ccwave.h"

namespace psi {
namespace ccenergy {

void CCEnergyWavefunction::Wmnij_build_sp() {
    dpdbuf4<float> A_anti, A;
    dpdbuf4<float> WMNIJ, Wmnij, WMnIj, W;
    dpdfile2<float> tIA, tia;
    dpdbuf4<float> Eijka, Eijka_anti, Eaijk, Eaijk_anti;
    dpdbuf4<float> D_anti, D, tauIJAB, tauijab, tauIjAb;

    if (params_.ref == 0) { /** RHF **/
        global_dpd_->buf4_init_sp(&A, PSIF_CC_AINTS, 0, 0, 0, 0, 0, 0, "A <ij|kl> sp");
        global_dpd_->buf4_copy_sp(&A, PSIF_CC_HBAR, "WMnIj_sp");
        global_dpd_->buf4_close_sp(&A);
    } 
    if (params_.ref == 0) { /** RHF **/
        global_dpd_->buf4_init_sp(&WMnIj, PSIF_CC_HBAR, 0, 0, 0, 0, 0, 0, "WMnIj_sp");
        global_dpd_->file2_init_sp(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA_sp");

        global_dpd_->buf4_init_sp(&Eaijk, PSIF_CC_EINTS, 0, 11, 0, 11, 0, 0, "E <ai|jk> sp");
        global_dpd_->contract244_sp(&tIA, &Eaijk, &WMnIj, 1, 0, 1, 1, 1);
        global_dpd_->buf4_close_sp(&Eaijk);

        global_dpd_->buf4_init_sp(&Eijka, PSIF_CC_EINTS, 0, 0, 10, 0, 10, 0, "E <ij|ka> sp");
        global_dpd_->contract424_sp(&Eijka, &tIA, &WMnIj, 3, 1, 0, 1, 1);
        global_dpd_->buf4_close_sp(&Eijka);

        global_dpd_->file2_close_sp(&tIA);
        global_dpd_->buf4_close_sp(&WMnIj);
    }   
    if (params_.ref == 0) { /** RHF **/
        global_dpd_->buf4_init_sp(&WMnIj, PSIF_CC_HBAR, 0, 0, 0, 0, 0, 0, "WMnIj_sp");
        global_dpd_->buf4_init_sp(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab> sp");
        global_dpd_->buf4_init_sp(&tauIjAb, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tauIjAb_sp");
        global_dpd_->contract444_sp(&D, &tauIjAb, &WMnIj, 0, 0, 1, 1);
        global_dpd_->buf4_close_sp(&tauIjAb);
        global_dpd_->buf4_close_sp(&D);
        global_dpd_->buf4_close_sp(&WMnIj);
    }
}
}  // namespace ccenergy
}  // namespace psi
