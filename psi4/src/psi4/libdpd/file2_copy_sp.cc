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
#include <cstring>
#include "dpd.h"

namespace psi {

int DPD::file2_copy_sp(dpdfile2<float> *InFile, int outfilenum, const char *label) {
    int h, row, col, my_irrep, rowtot, coltot;
    float ***matrix;
    dpdfile2<float> OutFile;

    my_irrep = InFile->my_irrep;

    file2_init_sp(&OutFile, outfilenum, InFile->my_irrep, InFile->params->pnum, InFile->params->qnum, label);

    file2_mat_init_sp(InFile);
    file2_mat_rd_sp(InFile);
    file2_mat_init_sp(&OutFile);

    for (h = 0; h < OutFile.params->nirreps; h++) {
        rowtot = OutFile.params->rowtot[h];
        coltot = OutFile.params->coltot[h ^ my_irrep];
        if (rowtot && coltot)
            memcpy((void *)&(OutFile.matrix[h][0][0]), (const void *)&(InFile->matrix[h][0][0]),
                   sizeof(float) * rowtot * coltot);
    }

    file2_mat_wrt_sp(&OutFile);
    file2_mat_close_sp(&OutFile);
    file2_mat_close_sp(InFile);
    file2_close_sp(&OutFile);

    return 0;
}

}  // namespace psi
