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
#include "psi4/cc/ccwave.h"

namespace psi {
namespace ccenergy {

void CCEnergyWavefunction::Z_build_sp() {
    dpdbuf4<float> Z_sp;
    dpdbuf4<float> tauIjAb, F, tau;

    timer_on("Z");

    if (params_.ref == 0) { /** RHF **/
        /* ZMbIj = <Mb|Ef> * tau(Ij,Ef) */
        /* OOC code added 3/23/05  -TDC */
        global_dpd_->buf4_init_sp(&Z_sp, PSIF_CC_MISC, 0, 10, 0, 10, 0, 0, "ZMbIj_sp");
        global_dpd_->buf4_init_sp(&F, PSIF_CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc> sp");
        global_dpd_->buf4_init_sp(&tau, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tauIjAb_sp");
        global_dpd_->contract444_sp(&F, &tau, &Z_sp, 0, 0, 1, 0);
        global_dpd_->buf4_close_sp(&tau);
        global_dpd_->buf4_close_sp(&F);
        global_dpd_->buf4_close_sp(&Z_sp);
       
    }  
    timer_off("Z");
}

}  // namespace ccenergy
}  // namespace psi
