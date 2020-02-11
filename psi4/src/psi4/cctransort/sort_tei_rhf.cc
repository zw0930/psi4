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

#include "psi4/libpsio/psio.hpp"
#include "psi4/libdpd/dpd.h"

namespace psi {
namespace cctransort {

void sort_tei_rhf(std::shared_ptr<PSIO> psio, int print) {
    dpdbuf4<double> K;
    dpdbuf4<float> K_test;

    psio->open(PSIF_CC_AINTS, PSIO_OPEN_OLD);
    global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, "ij", "kl", "i>=j+", "k>=l+", 0, "MO Ints (OO|OO)");
    global_dpd_->buf4_sort(&K, PSIF_CC_AINTS, prqs, "ij", "kl", "A <ij|kl>");

    //Test the buf4_copy function
    global_dpd_->buf4_close(&K);
    global_dpd_->buf4_init(&K, PSIF_CC_AINTS, 0, "ij", "kl", 0, "A <ij|kl>");    
    //global_dpd_->buf4_copy(&K, PSIF_CC_AINTS, "A <ij|kl> copy"); 
   // global_dpd_->buf4_init(&K_test, PSIF_CC_AINTS, 0, "ij", "kl", 0, "A <ij|kl> copy");

    //Cast and make a copy of buf4 in single-precision
    //if ((precision == 1)||(precision == 2)) global_dpd_->buf4_cast_copy(&K, PSIF_CC_AINTS, "ij", "kl", "A <ij|kl> sp");
    global_dpd_->buf4_cast_copy_dtof(&K, PSIF_CC_AINTS,  "A <ij|kl> sp");
    
    // Test the cast copy function
    //global_dpd_->buf4_init_sp(&K_test, PSIF_CC_AINTS, 0, "ij", "kl", 0, "A <ij|kl> sp");

    global_dpd_->buf4_close(&K);
    if (print > 6) {
        global_dpd_->buf4_init(&K, PSIF_CC_AINTS, 0, "ij", "kl", 0, "A <ij|kl>");
        global_dpd_->buf4_print(&K, "outfile", 1);
        global_dpd_->buf4_close(&K);
    }
    psio->close(PSIF_CC_AINTS, 1);

    psio->open(PSIF_CC_BINTS, PSIO_OPEN_OLD);
    global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, "ab", "cd", "a>=b+", "c>=d+", 0, "MO Ints (VV|VV)");
    global_dpd_->buf4_sort(&K, PSIF_CC_BINTS, prqs, "ab", "cd", "B <ab|cd>");
    //Cast and make a copy of buf4 in single-precision
    //if ((precision == 1)||(precision == 2)) global_dpd_->buf4_cast_copy(&K, PSIF_CC_AINTS, "ab", "cd", "B <ab|cd> sp");
    global_dpd_->buf4_close(&K);
    global_dpd_->buf4_init(&K, PSIF_CC_BINTS, 0, "ab", "cd", 0, "B <ab|cd>");
    global_dpd_->buf4_cast_copy_dtof(&K, PSIF_CC_BINTS, "B <ab|cd> sp");
    global_dpd_->buf4_close(&K);
    if (print > 6) {
        global_dpd_->buf4_init(&K, PSIF_CC_BINTS, 0, "ab", "cd", 0, "B <ab|cd>");
        global_dpd_->buf4_print(&K, "outfile", 1);
        global_dpd_->buf4_close(&K);
    }
    psio->close(PSIF_CC_BINTS, 1);

    psio->open(PSIF_CC_CINTS, PSIO_OPEN_OLD);
    global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, "ij", "ab", "i>=j+", "a>=b+", 0, "MO Ints (OO|VV)");
    global_dpd_->buf4_sort(&K, PSIF_CC_CINTS, prqs, "ia", "jb", "C <ia|jb>");
    //Cast and make a copy of buf4 in single-precision
    //if ((precision == 1)||(precision == 2)) global_dpd_->buf4_cast_copy(&K, PSIF_CC_AINTS, "ia", "jb", "C <ia|jb> sp");
    global_dpd_->buf4_close(&K);
    global_dpd_->buf4_init(&K, PSIF_CC_CINTS, 0, "ia", "jb", 0, "C <ia|jb>");
    global_dpd_->buf4_cast_copy_dtof(&K, PSIF_CC_CINTS, "C <ia|jb> sp");

    global_dpd_->buf4_close(&K);
    if (print > 6) {
        global_dpd_->buf4_init(&K, PSIF_CC_CINTS, 0, "ia", "jb", 0, "C <ia|jb>");
        global_dpd_->buf4_print(&K, "outfile", 1);
        global_dpd_->buf4_close(&K);
    }
    psio->close(PSIF_CC_CINTS, 1);

    psio->open(PSIF_CC_DINTS, PSIO_OPEN_OLD);
    global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, "ia", "jb", "ia", "jb", 0, "MO Ints (OV|OV)");
    global_dpd_->buf4_sort(&K, PSIF_CC_DINTS, prqs, "ij", "ab", "D <ij|ab>");
    //Cast and make a copy of buf4 in single-precision
    //if ((precision == 1)||(precision == 2)) global_dpd_->buf4_cast_copy(&K, PSIF_CC_AINTS, "jb", "ia", "D <jb|ia> sp");
    global_dpd_->buf4_close(&K);
    global_dpd_->buf4_init(&K, PSIF_CC_DINTS, 0, "ij", "ab", 0, "D <ij|ab>");
    global_dpd_->buf4_cast_copy_dtof(&K, PSIF_CC_DINTS, "D <ij|ab> sp");
    global_dpd_->buf4_close(&K);
    if (print > 6) {
        global_dpd_->buf4_init(&K, PSIF_CC_DINTS, 0, "ij", "ab", 0, "D <ij|ab>");
        global_dpd_->buf4_print(&K, "outfile", 1);
        global_dpd_->buf4_close(&K);
    }
    psio->close(PSIF_CC_DINTS, 1);

    psio->open(PSIF_CC_EINTS, PSIO_OPEN_OLD);
    global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, "ij", "ka", "i>=j+", "ka", 0, "MO Ints (OO|OV)");
    global_dpd_->buf4_sort(&K, PSIF_CC_EINTS, sqrp, "ai", "jk", "E <ai|jk>");
    //Cast and make a copy of buf4 in single-precision
    //if ((precision == 1)||(precision == 2)) global_dpd_->buf4_cast_copy(&K, PSIF_CC_AINTS, "ai", "jk", "E <ai|jk> sp");
    global_dpd_->buf4_close(&K);
    global_dpd_->buf4_init(&K, PSIF_CC_EINTS, 0, "ai", "jk", 0, "E <ai|jk>");
    global_dpd_->buf4_cast_copy_dtof(&K, PSIF_CC_EINTS, "E <ai|jk> sp");
    global_dpd_->buf4_close(&K);
    if (print > 6) {
        global_dpd_->buf4_init(&K, PSIF_CC_EINTS, 0, "ai", "jk", 0, "E <ai|jk>");
        global_dpd_->buf4_print(&K, "outfile", 1);
        global_dpd_->buf4_close(&K);
    }
    psio->close(PSIF_CC_EINTS, 1);

    psio->open(PSIF_CC_FINTS, PSIO_OPEN_OLD);
    global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, "ia", "bc", "ia", "b>=c+", 0, "MO Ints (OV|VV)");
    global_dpd_->buf4_sort(&K, PSIF_CC_FINTS, prqs, "ia", "bc", "F <ia|bc>");
    //Cast and make a copy of buf4 in single-precision
    //if ((precision == 1)||(precision == 2)) global_dpd_->buf4_cast_copy(&K, PSIF_CC_AINTS, "ia", "bc", "F <ia|bc> sp");
    global_dpd_->buf4_close(&K);
    global_dpd_->buf4_init(&K, PSIF_CC_FINTS, 0, "ia", "bc", 0, "F <ia|bc>");
    global_dpd_->buf4_cast_copy_dtof(&K, PSIF_CC_FINTS, "F <ia|bc> sp");

    global_dpd_->buf4_close(&K);
    if (print > 6) {
        global_dpd_->buf4_init(&K, PSIF_CC_FINTS, 0, "ia", "bc", 0, "F <ia|bc>");
        global_dpd_->buf4_print(&K, "outfile", 1);
        global_dpd_->buf4_close(&K);
    }
    global_dpd_->buf4_init(&K, PSIF_CC_FINTS, 0, "ia", "bc", 0, "F <ia|bc>");
    global_dpd_->buf4_sort(&K, PSIF_CC_FINTS, qpsr, "ai", "bc", "F <ai|bc>");
    global_dpd_->buf4_close(&K);
    global_dpd_->buf4_init(&K, PSIF_CC_FINTS, 0, "ai", "bc", 0, "F <ai|bc>");
    global_dpd_->buf4_cast_copy_dtof(&K, PSIF_CC_FINTS, "F <ai|bc> sp");
    global_dpd_->buf4_close(&K);
    psio->close(PSIF_CC_FINTS, 1);
}

}  // namespace cctransort
}  // namespace psi
