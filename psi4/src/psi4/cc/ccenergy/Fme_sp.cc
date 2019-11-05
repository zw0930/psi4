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
#include "psi4/cc/ccwave.h"

namespace psi {
namespace ccenergy {

void CCEnergyWavefunction::Fme_build_sp() {
    dpdfile2<float> FME_sp, Fme_sp, fIA_sp, fia_sp, tIA_sp, tia_sp;
    dpdbuf4<float> D_anti_sp, D_sp;

    if (params_.ref == 0) { /** RHF **/
        global_dpd_->file2_init_sp(&fIA_sp, PSIF_CC_OEI, 0, 0, 1, "fIA_sp");
        global_dpd_->file2_copy_sp(&fIA_sp, PSIF_CC_OEI, "FME_sp");
        global_dpd_->file2_close_sp(&fIA_sp);

        global_dpd_->file2_init_sp(&FME_sp, PSIF_CC_OEI, 0, 0, 1, "FME_sp");

        global_dpd_->buf4_init_sp(&D_anti_sp, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij||ab> sp");
        global_dpd_->buf4_init_sp(&D_sp, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab> sp");
        global_dpd_->file2_init_sp(&tIA_sp, PSIF_CC_OEI, 0, 0, 1, "tIA_sp");

        global_dpd_->dot13_sp(&tIA_sp, &D_anti_sp, &FME_sp, 0, 0, 1.0, 1.0);
        global_dpd_->dot13_sp(&tIA_sp, &D_sp, &FME_sp, 0, 0, 1.0, 1.0);

        global_dpd_->file2_close_sp(&tIA_sp);
        global_dpd_->buf4_close_sp(&D_anti_sp);
        global_dpd_->buf4_close_sp(&D_sp);

        global_dpd_->file2_close_sp(&FME_sp);
    }
}
}  // namespace ccenergy
}  // namespace psi
