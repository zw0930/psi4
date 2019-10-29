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

void CCEnergyWavefunction::ET2() {
    dpdfile2<double> tIA, tia;
    dpdbuf4<double> newtIJAB, newtijab, newtIjAb;
    dpdbuf4<double> E, t2, t2a, t2b;

    if (params_.ref == 0) { /** RHF **/
        global_dpd_->buf4_init(&newtIjAb, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "New tIjAb");

        global_dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");

        global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 11, 0, 11, 0, 0, "E <ai|jk>");
        global_dpd_->contract424(&E, &tIA, &newtIjAb, 1, 0, 0, -1, 1);
        global_dpd_->buf4_close(&E);

        global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 10, 0, 10, 0, 0, "E <ia|jk>");
        global_dpd_->contract244(&tIA, &E, &newtIjAb, 0, 0, 1, -1, 1);
        global_dpd_->buf4_close(&E);

        global_dpd_->file2_close(&tIA);

        global_dpd_->buf4_close(&newtIjAb);
    } else if (params_.ref == 1) { /** ROHF **/

        global_dpd_->buf4_init(&newtIJAB, PSIF_CC_TAMPS, 0, 2, 5, 2, 7, 0, "New tIJAB");
        global_dpd_->buf4_init(&newtijab, PSIF_CC_TAMPS, 0, 2, 5, 2, 7, 0, "New tijab");
        global_dpd_->buf4_init(&newtIjAb, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "New tIjAb");

        global_dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
        global_dpd_->file2_init(&tia, PSIF_CC_OEI, 0, 0, 1, "tia");

        /*** AA ***/
        global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 11, 2, 11, 0, 1, "E <ai|jk>");
        global_dpd_->buf4_init(&t2, PSIF_CC_TMP0, 0, 2, 5, 2, 5, 0, "T (I>J,AB)");
        global_dpd_->contract424(&E, &tIA, &t2, 1, 0, 0, -1, 0);
        global_dpd_->buf4_sort(&t2, PSIF_CC_TMP0, pqsr, 2, 5, "T (I>J,BA)");
        global_dpd_->buf4_close(&t2);
        global_dpd_->buf4_init(&t2a, PSIF_CC_TMP0, 0, 2, 5, 2, 5, 0, "T (I>J,AB)");
        global_dpd_->buf4_init(&t2b, PSIF_CC_TMP0, 0, 2, 5, 2, 5, 0, "T (I>J,BA)");
        global_dpd_->buf4_axpy(&t2b, &t2a, -1);
        global_dpd_->buf4_close(&t2b);
        global_dpd_->buf4_axpy(&t2a, &newtIJAB, 1);
        global_dpd_->buf4_close(&t2a);
        global_dpd_->buf4_close(&E);

        /*** BB ***/
        global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 11, 2, 11, 0, 1, "E <ai|jk>");
        global_dpd_->buf4_init(&t2, PSIF_CC_TMP0, 0, 2, 5, 2, 5, 0, "T (I>J,AB)");
        global_dpd_->contract424(&E, &tia, &t2, 1, 0, 0, -1, 0);
        global_dpd_->buf4_sort(&t2, PSIF_CC_TMP0, pqsr, 2, 5, "T (I>J,BA)");
        global_dpd_->buf4_close(&t2);
        global_dpd_->buf4_init(&t2a, PSIF_CC_TMP0, 0, 2, 5, 2, 5, 0, "T (I>J,AB)");
        global_dpd_->buf4_init(&t2b, PSIF_CC_TMP0, 0, 2, 5, 2, 5, 0, "T (I>J,BA)");
        global_dpd_->buf4_axpy(&t2b, &t2a, -1);
        global_dpd_->buf4_close(&t2b);
        global_dpd_->buf4_axpy(&t2a, &newtijab, 1);
        global_dpd_->buf4_close(&t2a);
        global_dpd_->buf4_close(&E);

        /*** AB ***/

        global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 11, 0, 11, 0, 0, "E <ai|jk>");
        global_dpd_->contract424(&E, &tia, &newtIjAb, 1, 0, 0, -1, 1);
        global_dpd_->buf4_close(&E);
        global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 10, 0, 10, 0, 0, "E <ia|jk>");
        global_dpd_->contract244(&tIA, &E, &newtIjAb, 0, 0, 1, -1, 1);
        global_dpd_->buf4_close(&E);

        global_dpd_->file2_close(&tIA);
        global_dpd_->file2_close(&tia);

        global_dpd_->buf4_close(&newtIJAB);
        global_dpd_->buf4_close(&newtijab);
        global_dpd_->buf4_close(&newtIjAb);

    } else if (params_.ref == 2) { /*** UHF ***/

        global_dpd_->buf4_init(&newtIJAB, PSIF_CC_TAMPS, 0, 2, 5, 2, 7, 0, "New tIJAB");
        global_dpd_->buf4_init(&newtijab, PSIF_CC_TAMPS, 0, 12, 15, 12, 17, 0, "New tijab");
        global_dpd_->buf4_init(&newtIjAb, PSIF_CC_TAMPS, 0, 22, 28, 22, 28, 0, "New tIjAb");

        global_dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
        global_dpd_->file2_init(&tia, PSIF_CC_OEI, 0, 2, 3, "tia");

        /*** AA ***/
        global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 21, 2, 21, 0, 1, "E <AI|JK>");
        global_dpd_->buf4_init(&t2, PSIF_CC_TMP0, 0, 2, 5, 2, 5, 0, "T (I>J,AB)");
        global_dpd_->contract424(&E, &tIA, &t2, 1, 0, 0, -1, 0);
        global_dpd_->buf4_sort(&t2, PSIF_CC_TMP0, pqsr, 2, 5, "T (I>J,BA)");
        global_dpd_->buf4_close(&t2);
        global_dpd_->buf4_init(&t2a, PSIF_CC_TMP0, 0, 2, 5, 2, 5, 0, "T (I>J,AB)");
        global_dpd_->buf4_init(&t2b, PSIF_CC_TMP0, 0, 2, 5, 2, 5, 0, "T (I>J,BA)");
        global_dpd_->buf4_axpy(&t2b, &t2a, -1);
        global_dpd_->buf4_close(&t2b);
        global_dpd_->buf4_axpy(&t2a, &newtIJAB, 1);
        global_dpd_->buf4_close(&t2a);
        global_dpd_->buf4_close(&E);

        /*** BB ***/
        global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 31, 12, 31, 10, 1, "E <ai|jk>");
        global_dpd_->buf4_init(&t2, PSIF_CC_TMP0, 0, 12, 15, 12, 15, 0, "T (i>j,ab)");
        global_dpd_->contract424(&E, &tia, &t2, 1, 0, 0, -1, 0);
        global_dpd_->buf4_sort(&t2, PSIF_CC_TMP0, pqsr, 12, 15, "T (i>j,ba)");
        global_dpd_->buf4_close(&t2);
        global_dpd_->buf4_init(&t2a, PSIF_CC_TMP0, 0, 12, 15, 12, 15, 0, "T (i>j,ab)");
        global_dpd_->buf4_init(&t2b, PSIF_CC_TMP0, 0, 12, 15, 12, 15, 0, "T (i>j,ba)");
        global_dpd_->buf4_axpy(&t2b, &t2a, -1);
        global_dpd_->buf4_close(&t2b);
        global_dpd_->buf4_axpy(&t2a, &newtijab, 1);
        global_dpd_->buf4_close(&t2a);
        global_dpd_->buf4_close(&E);

        /*** AB ***/

        global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 22, 26, 22, 26, 0, "E <Ij|Ak>");
        global_dpd_->contract424(&E, &tia, &newtIjAb, 3, 0, 0, -1, 1);
        global_dpd_->buf4_close(&E);
        global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 24, 22, 24, 22, 0, "E <Ia|Jk>");
        global_dpd_->contract244(&tIA, &E, &newtIjAb, 0, 0, 1, -1, 1);
        global_dpd_->buf4_close(&E);

        global_dpd_->file2_close(&tIA);
        global_dpd_->file2_close(&tia);

        global_dpd_->buf4_close(&newtIJAB);
        global_dpd_->buf4_close(&newtijab);
        global_dpd_->buf4_close(&newtIjAb);
    }
}


void CCEnergyWavefunction::ET2_mp() {
    dpdfile2<float> tIA_sp;
    dpdbuf4<double> newtIjAb, buf4_tmp_double;
    dpdbuf4<float> buf4_tmp_float;
    dpdbuf4<float> E_sp;

    if (params_.ref == 0) { /** RHF **/
        global_dpd_->buf4_init(&newtIjAb, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "New tIjAb");
        global_dpd_->buf4_init_sp(&buf4_tmp_float, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "buf4_tmp_float");
        global_dpd_->file2_init(&tIA_sp, PSIF_CC_OEI, 0, 0, 1, "tIA_sp");

        global_dpd_->buf4_init_sp(&E_sp, PSIF_CC_EINTS, 0, 11, 0, 11, 0, 0, "E <ai|jk> sp");
        global_dpd_->contract424_sp(&E_sp, &tIA_sp, &buf4_tmp_float, 1, 0, 0, -1, 0);
        global_dpd_->buf4_close_sp(&E_sp);
        global_dpd_->buf4_cast_ftod_copy(&buf4_tmp_float, PSIF_CC_TAMPS, "buf4_tmp_double");
        global_dpd_->buf4_init(&buf4_tmp_double, PSIF_CC_TMP0, 0, 0, 5, 0, 5, 0, "buf4_tmp_double");
        global_dpd_->buf4_axpy(&buf4_tmp_double, &newtIjAb, 1);
        global_dpd_->buf4_close_sp(&buf4_tmp_float);
        global_dpd_->buf4_close(&buf4_tmp_double);

        global_dpd_->buf4_init_sp(&E_sp, PSIF_CC_EINTS, 0, 10, 0, 10, 0, 0, "E <ia|jk> sp");
        global_dpd_->buf4_init_sp(&buf4_tmp_float, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "buf4_tmp_float");
        global_dpd_->contract244_sp(&tIA_sp, &E_sp, &buf4_tmp_float, 0, 0, 1, -1, 0);
        global_dpd_->buf4_cast_ftod_copy(&buf4_tmp_float, PSIF_CC_TAMPS, "buf4_tmp_double" );
        global_dpd_->buf4_init_dp(&buf4_tmp_double, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "buf4_tmp_double");
        glabal_dpd_->buf4_axpy(&buf4_tmp_double, &newtIjAb, 1);
        global_dpd_->buf4_close_sp(&E_sp);
        global_dpd_->buf4_close(&buf4_tmp_double);
        global_dpd_->buf4_close_sp(&buf4_tmp_float);

        global_dpd_->file2_close_sp(&tIA_sp);

        global_dpd_->buf4_close(&newtIjAb);
    } 
}


void CCEnergyWavefunction::ET2_sp() {
    dpdfile2<float> tIA;
    dpdbuf4<float> newtIjAb;
    dpdbuf4<float> E;

    if (params_.ref == 0) { /** RHF **/
        global_dpd_->buf4_init_sp(&newtIjAb, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "New tIjAb sp");

        global_dpd_->file2_init_sp(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA_sp");

        global_dpd_->buf4_init_sp(&E, PSIF_CC_EINTS, 0, 11, 0, 11, 0, 0, "E <ai|jk> sp");
        global_dpd_->contract424_sp(&E, &tIA, &newtIjAb, 1, 0, 0, -1, 1);
        global_dpd_->buf4_close_sp(&E);

        global_dpd_->buf4_init_sp(&E, PSIF_CC_EINTS, 0, 10, 0, 10, 0, 0, "E <ia|jk>");
        global_dpd_->contract244_sp(&tIA, &E, &newtIjAb, 0, 0, 1, -1, 1);
        global_dpd_->buf4_close_sp(&E);

        global_dpd_->file2_close_sp(&tIA);

        global_dpd_->buf4_close_sp(&newtIjAb);
    } 
}
}  // namespace ccenergy
}  // namespace psi
