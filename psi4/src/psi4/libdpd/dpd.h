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
#ifndef _psi_src_lib_libdpd_dpd_h
#define _psi_src_lib_libdpd_dpd_h

#include <cstdio>
#include <string>
#include "psi4/psifiles.h"
#include "psi4/libpsio/config.h"
#include "psi4/pragma.h"
PRAGMA_WARNING_PUSH
PRAGMA_WARNING_IGNORE_DEPRECATED_DECLARATIONS
#include <memory>
PRAGMA_WARNING_POP
#include <vector>
#include "psi4/psi4-dec.h"

// Testing -TDC
#include "dpdmospace.h"

namespace psi {

#define T3_TIMER_ON (0)

#define DPD_BIGNUM 2147483647 /* the four-byte signed int limit */
/* #define ALL_BUF4_SORT_OOC */


struct dpdparams4 {
    int nirreps;   /* No. of irreps */
    int pqnum;     /* Pair number for the row indices */
    int rsnum;     /* Pair number for the column indices */
    int *rowtot;   /* Row dimension for each irrep */
    int *coltot;   /* Column dimension for each irrep */
    int **rowidx;  /* Row index lookup array */
    int **colidx;  /* Column index lookup array */
    int ***roworb; /* Row index -> orbital index lookup array */
    int ***colorb; /* Column index -> orbital index lookup array */
    int *ppi;      /* Number of p indices per irrep */
    int *qpi;      /* Number of q indices per irrep */
    int *rpi;      /* Number of r indices per irrep */
    int *spi;      /* Number of s indices per irrep */
    int *poff;     /* Orbital offset for p */
    int *qoff;     /* Orbital offset for q */
    int *roff;     /* Orbital offset for r */
    int *soff;     /* Orbital offset for s */
    int *psym;     /* Orbital symmetry for index p */
    int *qsym;     /* Orbital symmetry for index q */
    int *rsym;     /* Orbital symmetry for index r */
    int *ssym;     /* Orbital symmetry for index s */
    int perm_pq;   /* Can p and q be permuted? */
    int perm_rs;   /* Can r and s be permuted? */
    int peq;       /* Can p and q be equal? */
    int res;       /* Can r and s be equal? */
    int **start13; /* returns the starting row of irrep matrix h for orbital index p */
};

struct dpdfile4 {
    int dpdnum; /* dpd structure reference */
    char label[PSIO_KEYLEN];
    int filenum;
    int my_irrep;         /* Total irrep of this quantity */
    psio_address *lfiles; /* File address for each submatrix by ROW irrep */
    dpdparams4 *params;
    int incore;
    double ***matrix;
};

struct dpdfile4_sp {
    int dpdnum; /* dpd structure reference */
    char label[PSIO_KEYLEN];
    int filenum;
    int my_irrep;         /* Total irrep of this quantity */
    psio_address *lfiles; /* File address for each submatrix by ROW irrep */
    dpdparams4 *params;
    int incore;
    float ***matrix;
};


template <typename U>
struct dpdshift4 {
    int shift_type;
    int **rowtot;
    int **coltot;
    U ****matrix;
};

template <typename U> 
struct dpdbuf4<double>{
    int dpdnum; /* dpd structure reference */
    int anti;   /* Is this buffer antisymmetric? */
    dpdparams4 *params;
    dpdfile4 file;
    dpdfile4_sp file_sp;
    dpdshift4<U> shift;
    int **row_offset;
    int **col_offset;
    U ***matrix;
};




template <typename U>
struct dpdtrans4 {
    U ***matrix;
    dpdshift4<U> shift;
    dpdbuf4<U> buf;
};



struct dpdparams2 {
    int nirreps; /* No. of irreps */
    int pnum;
    int qnum;
    int *rowtot;  /* Row dimension for each submatrix */
    int *coltot;  /* Column dimension for each submatrix */
    int *rowidx;  /* Row index lookup array */
    int *colidx;  /* Column index lookup array */
    int **roworb; /* Row index -> orbital index lookup array */
    int **colorb; /* Column index -> orbital index lookup array */
    int *ppi;     /* Number of p indices per irrep */
    int *qpi;     /* Number of q indices per irrep */
    int *poff;    /* Orbital offset for p */
    int *qoff;    /* Orbital offset for q */
    int *psym;    /* Orbital symmetry for index p */
    int *qsym;    /* Orbital symmetry for index q */
};

template <typename U> 
struct dpdfile2<double>{
    int dpdnum; /* dpd structure reference */
    char label[PSIO_KEYLEN];
    int filenum;
    int my_irrep;
    psio_address *lfiles;
    dpdparams2 *params;
    int incore;
    U ***matrix;
};


/* DPD File4 Cache entries */
struct dpd_file4_cache_entry {
    int dpdnum;                  /* dpd structure reference */
    int filenum;                 /* libpsio unit number */
    int irrep;                   /* overall symmetry */
    int pqnum;                   /* dpd pq value */
    int rsnum;                   /* dpd rs value */
    char label[PSIO_KEYLEN];     /* libpsio TOC keyword */
    double ***matrix;            /* pointer to irrep blocks */
    int size;                    /* size of entry in double words */
    size_t access;               /* access time */
    size_t usage;                /* number of accesses */
    size_t priority;             /* priority level */
    int lock;                    /* auto-deletion allowed? */
    int clean;                   /* has this file4 changed? */
    dpd_file4_cache_entry *next; /* pointer to next cache entry */
    dpd_file4_cache_entry *last; /* pointer to previous cache entry */
};

/* DPD File2 Cache entries */
struct dpd_file2_cache_entry {
    dpd_file2_cache_entry() : next(nullptr), last(nullptr) {}
    int dpdnum;                  /* dpd structure reference */
    int filenum;                 /* libpsio unit number */
    int irrep;                   /* overall symmetry */
    int pnum;                    /* dpd p value */
    int qnum;                    /* dpd q value */
    char label[PSIO_KEYLEN];     /* libpsio TOC keyword */
    double ***matrix;            /* pointer to irrep blocks */
    int size;                    /* size of entry in double words */
    int clean;                   /* has this file2 changed? */
    dpd_file2_cache_entry *next; /* pointer to next cache entry */
    dpd_file2_cache_entry *last; /* pointer to previous cache entry */
};

/* DPD global parameter set */
struct dpd_data {
    int nirreps;
    int num_subspaces;
    int num_pairs;
    int *numorbs;
    int **orboff;
    int **pairtot;
    int **orbspi;
    int **orbsym;
    int **orbidx2;
    int ***pairidx;
    int ***orbs2;
    int ****pairorb;
    dpdparams2 **params2;
    dpdparams4 **params4;
};

struct thread_data {
    dpdbuf4<double> *CIjAb;
    dpdbuf4<double> *WAbEi;
    dpdbuf4<double> *WMbIj;
    int do_singles;
    dpdbuf4<double> *Dints;
    dpdfile2<double> *SIA;
    int do_doubles;
    dpdfile2<double> *FME;
    dpdbuf4<double> *WmAEf;
    dpdbuf4<double> *WMnIe;
    dpdbuf4<double> *SIjAb;
    int *occpi;
    int *occ_off;
    int *virtpi;
    int *vir_off;
    double omega;
    dpdfile2<double> *fIJ;
    dpdfile2<double> *fAB;
    int Gi;
    int Gj;
    int Gk;
    int first_ijk;
    int last_ijk;
    std::string outfile;
    int thr_id;
    dpdfile2<double> SIA_local;
    dpdbuf4<double> SIjAb_local;
    int newtrips;
};

struct dpd_gbl {
    long int memory;    /* Total memory requested by the user */
    long int memused;   /* Total memory used (cache + other) */
    long int memcache;  /* Total memory in cache (locked and unlocked) */
    long int memlocked; /* Total memory locked in the cache */

    // The default C'tor will zero everything out properly
    dpd_gbl()
        : file2_cache(nullptr),
          file4_cache(nullptr),
          file4_cache_most_recent(0),
          file4_cache_least_recent(1),
          file4_cache_lru_del(0),
          file4_cache_low_del(0) {}
    dpd_file2_cache_entry *file2_cache;
    dpd_file4_cache_entry *file4_cache;
    size_t file4_cache_most_recent;
    size_t file4_cache_least_recent;
    size_t file4_cache_lru_del;
    size_t file4_cache_low_del;
    int cachetype;
    int *cachefiles;
    int **cachelist;
    dpd_file4_cache_entry *file4_cache_priority;
};

/* Useful for the generalized 4-index sorting function */
enum indices {
    pqrs,
    pqsr,
    prqs,
    prsq,
    psqr,
    psrq,
    qprs,
    qpsr,
    qrps,
    qrsp,
    qspr,
    qsrp,
    rqps,
    rqsp,
    rpqs,
    rpsq,
    rsqp,
    rspq,
    sqrp,
    sqpr,
    srqp,
    srpq,
    spqr,
    sprq
};

/* Useful for the 3-index sorting function dpd_3d_sort() */
enum pattern { abc, acb, cab, cba, bca, bac };

class PSI_API DPD {
   public:
    // These used to live in the dpd_data struct
    int nirreps;
    int num_subspaces;
    int num_pairs;
    int *numorbs;
    int **orboff;
    int **pairtot;
    int **orbspi;
    int **orbsym;
    int **orbidx2;
    int ***pairidx;
    int ***orbs2;
    int ****pairorb;
    dpdparams2 **params2;
    dpdparams4 **params4;

    std::vector<DPDMOSpace> moSpaces;

    DPD(int dpd_num, int nirreps, long int memory, int cachetype, int *cachefiles, int **cachelist,
        dpd_file4_cache_entry *priority, int num_subspaces, std::vector<int *> &spaceArrays);
    DPD(int dpd_num, int nirreps, long int memory, int cachetype, int *cachefiles, int **cachelist,
        dpd_file4_cache_entry *priority, int num_subspaces, std::vector<DPDMOSpace> &moSpaces);
    DPD();

    ~DPD();

    int init(int dpd_num, int nirreps, long int memory, int cachetype, int *cachefiles, int **cachelist,
             dpd_file4_cache_entry *priority, int num_subspaces, ...);
    int init(int dpd_num, int nirreps, long int memory, int cachetype, int *cachefiles, int **cachelist,
             dpd_file4_cache_entry *priority, int num_subspaces, std::vector<int *> &spaceArrays);

    void dpd_error(const char *caller, std::string out_fname);

    double **dpd_block_matrix(size_t n, size_t m);
    void free_dpd_block(double **array, size_t n, size_t m);

    int contract222(dpdfile2<double> *X, dpdfile2<double> *Y, dpdfile2<double> *Z, int target_X, int target_Y, double alpha, double beta);
    int contract442(dpdbuf4<double> *X, dpdbuf4<double> *Y, dpdfile2<double> *Z, int target_X, int target_Y, double alpha, double beta);
    int contract422(dpdbuf4<double> *X, dpdfile2<double> *Y, dpdfile2<double> *Z, int trans_Y, int trans_Z, double alpha, double beta);
    int contract244(dpdfile2<double> *X, dpdbuf4<double> *Y, dpdbuf4<double> *Z, int sum_X, int sum_Y, int trans_Z, double alpha, double beta);
    int contract424(dpdbuf4<double> *X, dpdfile2<double> *Y, dpdbuf4<double> *Z, int sum_X, int sum_Y, int trans_Z, double alpha, double beta);
    int contract444(dpdbuf4<double> *X, dpdbuf4<double> *Y, dpdbuf4<double> *Z, int target_X, int target_Y, double alpha, double beta);
    // not used
    int contract444_df(dpdbuf4<double> *B, dpdbuf4<double> *tau_in, dpdbuf4<double> *tau_out, double alpha, double beta);

    /* Need to consolidate these routines into one general function */
    int dot23(dpdfile2<double> *T, dpdbuf4<double> *I, dpdfile2<double> *Z, int transt, int transz, double alpha, double beta);
    int dot24(dpdfile2<double> *T, dpdbuf4<double> *I, dpdfile2<double> *Z, int transt, int transz, double alpha, double beta);
    int dot13(dpdfile2<double> *T, dpdbuf4<double> *I, dpdfile2<double> *Z, int transt, int transz, double alpha, double beta);
    int dot14(dpdfile2<double> *T, dpdbuf4<double> *I, dpdfile2<double> *Z, int transt, int transz, double alpha, double beta);

    int trace42_13(dpdbuf4<double> *A, dpdfile2<double> *B, int transb, double alpha, double beta);

    int file2_init(dpdfile2<double> *File, int filenum, int irrep, int pnum, int qnum, const char *label);
    int file2_close(dpdfile2<double> *File);
    int file2_mat_init(dpdfile2<double> *File);
    int file2_mat_close(dpdfile2<double> *File);
    int file2_mat_rd(dpdfile2<double> *File);
    int file2_mat_wrt(dpdfile2<double> *File);
    int file2_mat_wrt_sp(dpdfile2<float> *File);
//**
    int file2_print(dpdfile2<double> *File, std::string out_fname);
    int file2_mat_print(dpdfile2<double> *File, std::string out_fname);
    int file2_copy(dpdfile2<double> *InFile, int outfilenum, const char *label);
    int file2_copy_sp(dpdfile2<float> *InFile, int outfilenum, const char *label);

// ** 
    int file2_cast_copy_dtof(dpdfile2<double> *InFile, int outfilenum, const char *label);
    int file2_cast_copy_ftod(dpdfile2<float> *InFile, int outfilenum, const char *label);



    int file2_dirprd(dpdfile2<double> *FileA, dpdfile2<double> *FileB);
    double file2_dot(dpdfile2<double> *FileA, dpdfile2<double> *FileB);
    int file2_scm(dpdfile2<double> *InFile, double alpha);
//**   
    double file2_dot_self(dpdfile2<double> *BufX);
    double file2_trace(dpdfile2<double> *InFile);
    int file2_axpy(dpdfile2<double> *FileA, dpdfile2<double> *FileB, double alpha, int transA);
    int file2_axpy_sp(dpdfile2<float> *FileA, dpdfile2<float> *FileB, float alpha, int transA);
    int file2_axpbycz(dpdfile2<double> *FileA, dpdfile2<double> *FileB, dpdfile2<double> *FileC, double a, double b, double c);
    int file2_axpbycz_sp(dpdfile2<float> *FileA, dpdfile2<float> *FileB, dpdfile2<float> *FileC, float a, float b, float c);


    int file4_init(dpdfile4 *File, int filenum, int irrep, int pqnum, int rsnum, const char *label);
    int file4_init_nocache(dpdfile4 *File, int filenum, int irrep, int pqnum, int rsnum, const char *label);
    int file4_close(dpdfile4 *File);
    int file4_mat_irrep_init(dpdfile4 *File, int irrep);
    int file4_mat_irrep_close(dpdfile4 *File, int irrep);
    int file4_mat_irrep_rd(dpdfile4 *File, int irrep);
    int file4_mat_irrep_wrt(dpdfile4 *File, int irrep);
    int file4_mat_irrep_row_init(dpdfile4 *File, int irrep);
    int file4_mat_irrep_row_close(dpdfile4 *File, int irrep);
    int file4_mat_irrep_row_rd(dpdfile4 *File, int irrep, int row);
    int file4_mat_irrep_row_wrt(dpdfile4 *File, int irrep, int row);
    int file4_mat_irrep_row_zero(dpdfile4 *File, int irrep, int row);
    int file4_print(dpdfile4 *File, std::string out_fname);
    int file4_mat_irrep_rd_block(dpdfile4 *File, int irrep, int start_pq, int num_pq);
    int file4_mat_irrep_wrt_block(dpdfile4 *File, int irrep, int start_pq, int num_pq);

    int buf4_init(dpdbuf4<double> *Buf, int inputfile, int irrep, int pqnum, int rsnum, int file_pqnum, int file_rsnum,
                  int anti, const char *label);
    int buf4_init(dpdbuf4<double> *Buf, int inputfile, int irrep, std::string pq, std::string rs, std::string file_pq,
                  std::string file_rs, int anti, const char *label);
    int buf4_init(dpdbuf4<double> *Buf, int inputfile, int irrep, std::string pq, std::string rs, int anti, const char *label);
    int pairnum(std::string);
//*
    double buf4_trace(dpdbuf4<double> *Buf);
    int buf4_close(dpdbuf4<double> *Buf);
    int buf4_mat_irrep_init(dpdbuf4<double> *Buf, int irrep);
    int buf4_mat_irrep_close(dpdbuf4<double> *Buf, int irrep);
    int buf4_mat_irrep_rd(dpdbuf4<double> *Buf, int irrep);
    int buf4_mat_irrep_wrt(dpdbuf4<double> *Buf, int irrep);
//**
    int buf4_print(dpdbuf4<double> *Buf, std::string out_fname, int print_data);
    int buf4_copy(dpdbuf4<double> *InBuf, int outfilenum, const char *label);

//**
    int buf4_cast_copy_dtof(dpdbuf4<double> *InBuf, int outfilenum, const char *label);
    int buf4_cast_copy_ftod(dpdbuf4<float> *InBuf, int outfilenum, const char *label);




    int buf4_sort(dpdbuf4<double> *InBuf, int outfilenum, enum indices index, int pqnum, int rsnum, const char *label);
    int buf4_sort(dpdbuf4<double> *InBuf, int outfilenum, enum indices index, std::string pq, std::string rs,
                  const char *label);
    int buf4_sort_ooc(dpdbuf4<double> *InBuf, int outfilenum, enum indices index, int pqnum, int rsnum, const char *label);
    int buf4_sort_axpy(dpdbuf4<double> *InBuf, int outfilenum, enum indices index, int pqnum, int rsnum, const char *label,
                       double alpha);
    
//**
    int buf4_axpy(dpdbuf4<double> *BufX, dpdbuf4<double> *BufY, double alpha);
    int buf4_axpbycz(dpdbuf4<double> *FileA, dpdbuf4<double> *FileB, dpdbuf4<double> *FileC, double a, double b, double c);
//*
    int buf4_dirprd(dpdbuf4<double> *BufA, dpdbuf4<double> *BufB);
    double buf4_dot(dpdbuf4<double> *BufA, dpdbuf4<double> *BufB);
    double buf4_dot_self(dpdbuf4<double> *BufX);
//**
    int buf4_scm(dpdbuf4<double> *InBuf, double alpha);
    int buf4_scmcopy(dpdbuf4<double> *InBuf, int outfilenum, const char *label, double alpha);
    int buf4_symm(dpdbuf4<double> *Buf);
    int buf4_symm2(dpdbuf4<double> *Buf1, dpdbuf4<double> *Buf2);
    int buf4_mat_irrep_shift13(dpdbuf4<double> *Buf, int irrep);
    int buf4_mat_irrep_shift31(dpdbuf4<double> *Buf, int irrep);
    int buf4_mat_irrep_row_init(dpdbuf4<double> *Buf, int irrep);
    int buf4_mat_irrep_row_close(dpdbuf4<double> *Buf, int irrep);
    int buf4_mat_irrep_row_zero(dpdbuf4<double> *Buf, int irrep, int row);
    int buf4_mat_irrep_row_rd(dpdbuf4<double> *Buf, int irrep, int pq);
    int buf4_mat_irrep_row_wrt(dpdbuf4<double> *Buf, int irrep, int pq);
//**
    int buf4_mat_irrep_init_block(dpdbuf4<double> *Buf, int irrep, int num_pq);
    int buf4_mat_irrep_close_block(dpdbuf4<double> *Buf, int irrep, int num_pq);
    int buf4_mat_irrep_rd_block(dpdbuf4<double> *Buf, int irrep, int start_pq, int num_pq);
    int buf4_mat_irrep_wrt_block(dpdbuf4<double> *Buf, int irrep, int start_pq, int num_pq);
//*
    int buf4_dump(dpdbuf4<double> *DPDBuf, struct iwlbuf *IWLBuf, int *prel, int *qrel, int *rrel, int *srel, int bk_pack,
                  int swap23);
    int trans4_init(dpdtrans4<double> *Trans, dpdbuf4<double> *Buf);
    int trans4_close(dpdtrans4<double> *Trans);
    int trans4_mat_irrep_init(dpdtrans4<double> *Trans, int irrep);
    int trans4_mat_irrep_close(dpdtrans4<double> *Trans, int irrep);
    int trans4_mat_irrep_rd(dpdtrans4<double> *Trans, int irrep);
    int trans4_mat_irrep_wrt(dpdtrans4<double> *Trans, int irrep);
//*
    int trans4_mat_irrep_shift13(dpdtrans4<double> *Trans, int irrep);
    int trans4_mat_irrep_shift31(dpdtrans4<double> *Trans, int irrep);

    int mat4_irrep_print(double **matrix, dpdparams4 *Params, int irrep, int my_irrep, std::string out_fname);
    int mat4_irrep_print_sp(float **matrix, dpdparams4 *Params, int irrep, int my_irrep, std::string out_fname);



    void file2_cache_init();
    void file2_cache_close();
    void file2_cache_print(std::string out_fname);
    dpd_file2_cache_entry *file2_cache_scan(int filenum, int irrep, int pnum, int qnum, const char *label, int dpdnum);
    dpd_file2_cache_entry *dpd_file2_cache_last();
    int file2_cache_add(dpdfile2<double> *File);
    int file2_cache_del(dpdfile2<double> *File);
    int file4_cache_del_low();
    void file2_cache_dirty(dpdfile2<double> *File);

    void file4_cache_init();
    void file4_cache_close();
    void file4_cache_print(std::string out_fname);
    void file4_cache_print_screen();
//*
    int file4_cache_get_priority(dpdfile4 *File);

    dpd_file4_cache_entry *file4_cache_scan(int filenum, int irrep, int pqnum, int rsnum, const char *label,
                                            int dpdnum);
    dpd_file4_cache_entry *file4_cache_last();
    int file4_cache_add(dpdfile4 *File, size_t priority);
    int file4_cache_del(dpdfile4 *File);
    dpd_file4_cache_entry *file4_cache_find_lru();
    int file4_cache_del_lru();
    void file4_cache_dirty(dpdfile4 *File);
    void file4_cache_lock(dpdfile4 *File);
    void file4_cache_unlock(dpdfile4 *File);

    void sort_3d(double ***Win, double ***Wout, int nirreps, int h, int *rowtot, int **rowidx, int ***roworb, int *asym,
                 int *bsym, int *aoff, int *boff, int *cpi, int *coff, int **rowidx_out, enum pattern index, int sum);

    void T3_AAA(double ***W1, int nirreps, int I, int Gi, int J, int Gj, int K, int Gk, dpdbuf4<double> *T2, dpdbuf4<double> *F,
                dpdbuf4<double> *E, dpdfile2<double> *fIJ, dpdfile2<double> *fAB, int *occpi, int *occ_off, int *virtpi, int *vir_off,
                double omega);

    void T3_AAB(double ***W1, int nirreps, int I, int Gi, int J, int Gj, int K, int Gk, dpdbuf4<double> *T2AA, dpdbuf4<double> *T2AB,
                dpdbuf4<double> *T2BA, dpdbuf4<double> *FAA, dpdbuf4<double> *FAB, dpdbuf4<double> *FBA, dpdbuf4<double> *EAA, dpdbuf4<double> *EAB, dpdbuf4<double> *EBA,
                dpdfile2<double> *fIJ, dpdfile2<double> *fij, dpdfile2<double> *fAB, dpdfile2<double> *fab, int *aoccpi, int *aocc_off, int *boccpi,
                int *bocc_off, int *avirtpi, int *avir_off, int *bvirtpi, int *bvir_off, double omega);

    void T3_RHF(double ***W1, int nirreps, int I, int Gi, int J, int Gj, int K, int Gk, dpdbuf4<double> *T2, dpdbuf4<double> *F,
                dpdbuf4<double> *E, dpdfile2<double> *fIJ, dpdfile2<double> *fAB, int *occpi, int *occ_off, int *virtpi, int *vir_off,
                double omega);

    void T3_RHF_ic(double ***W1, int nirreps, int I, int Gi, int J, int Gj, int K, int Gk, dpdbuf4<double> *T2, dpdbuf4<double> *F,
                   dpdbuf4<double> *E, dpdfile2<double> *fIJ, dpdfile2<double> *fAB, int *occpi, int *occ_off, int *virtpi, int *vir_off,
                   double omega);

    void cc3_sigma_RHF(dpdbuf4<double> *CIjAb, dpdbuf4<double> *WAbEi, dpdbuf4<double> *WMbIj, int do_singles, dpdbuf4<double> *Dints, dpdfile2<double> *SIA,
                       int do_doubles, dpdfile2<double> *FME, dpdbuf4<double> *WAmEf, dpdbuf4<double> *WMnIe, dpdbuf4<double> *SIjAb, int *occpi,
                       int *occ_off, int *virtpi, int *vir_off, double omega, std::string out_fname, int newtrips);

    void cc3_sigma_RHF_ic(dpdbuf4<double> *CIjAb, dpdbuf4<double> *WAbEi, dpdbuf4<double> *WMbIj, int do_singles, dpdbuf4<double> *Dints, dpdfile2<double> *SIA,
                          int do_doubles, dpdfile2<double> *FME, dpdbuf4<double> *WAmEf, dpdbuf4<double> *WMnIe, dpdbuf4<double> *SIjAb, int *occpi,
                          int *occ_off, int *virtpi, int *vir_off, double omega, std::string out_fname, int nthreads,
                          int newtrips);

    void cc3_sigma_UHF_AAA(dpdbuf4<double> *CMNEF, dpdbuf4<double> *WABEI, dpdbuf4<double> *WMBIJ, int do_singles, dpdbuf4<double> *Dints_anti,
                           dpdfile2<double> *SIA, int do_doubles, dpdfile2<double> *FME, dpdbuf4<double> *WMAFE, dpdbuf4<double> *WMNIE, dpdbuf4<double> *SIJAB,
                           int *aoccpi, int *aocc_off, int *avirtpi, int *avir_off, double omega,
                           std::string out_fname);

    void cc3_sigma_UHF_BBB(dpdbuf4<double> *Cmnef, dpdbuf4<double> *Wabei, dpdbuf4<double> *Wmbij, int do_singles, dpdbuf4<double> *Dijab_anti,
                           dpdfile2<double> *Sia, int do_doubles, dpdfile2<double> *Fme, dpdbuf4<double> *Wmafe, dpdbuf4<double> *Wmnie, dpdbuf4<double> *Sijab,
                           int *boccpi, int *bocc_off, int *bvirtpi, int *bvir_off, double omega,
                           std::string out_fname);

    void cc3_sigma_UHF_AAB(dpdbuf4<double> *C2AA, dpdbuf4<double> *C2AB, dpdbuf4<double> *C2BA, dpdbuf4<double> *FAA, dpdbuf4<double> *FAB, dpdbuf4<double> *FBA,
                           dpdbuf4<double> *EAA, dpdbuf4<double> *EAB, dpdbuf4<double> *EBA, int do_singles, dpdbuf4<double> *DAA, dpdbuf4<double> *DAB,
                           dpdfile2<double> *SIA, dpdfile2<double> *Sia, int do_doubles, dpdfile2<double> *FME, dpdfile2<double> *Fme, dpdbuf4<double> *WMAFE,
                           dpdbuf4<double> *WMaFe, dpdbuf4<double> *WmAfE, dpdbuf4<double> *WMNIE, dpdbuf4<double> *WMnIe, dpdbuf4<double> *WmNiE,
                           dpdbuf4<double> *SIJAB, dpdbuf4<double> *SIjAb, int *aoccpi, int *aocc_off, int *boccpi, int *bocc_off,
                           int *avirtpi, int *avir_off, int *bvirtpi, int *bvir_off, double omega,
                           std::string out_fname);

    void cc3_sigma_UHF_BBA(dpdbuf4<double> *C2BB, dpdbuf4<double> *C2AB, dpdbuf4<double> *C2BA, dpdbuf4<double> *FBB, dpdbuf4<double> *FAB, dpdbuf4<double> *FBA,
                           dpdbuf4<double> *EBB, dpdbuf4<double> *EAB, dpdbuf4<double> *EBA, int do_singles, dpdbuf4<double> *DBB, dpdbuf4<double> *DBA,
                           dpdfile2<double> *SIA, dpdfile2<double> *Sia, int do_doubles, dpdfile2<double> *FME, dpdfile2<double> *Fme, dpdbuf4<double> *Wmafe,
                           dpdbuf4<double> *WMaFe, dpdbuf4<double> *WmAfE, dpdbuf4<double> *Wmnie, dpdbuf4<double> *WMnIe, dpdbuf4<double> *WmNiE,
                           dpdbuf4<double> *Sijab, dpdbuf4<double> *SIjAb, int *aoccpi, int *aocc_off, int *boccpi, int *bocc_off,
                           int *avirtpi, int *avir_off, int *bvirtpi, int *bvir_off, double omega,
                           std::string out_fname);


// *functions for single-precision and mixed-precision
    float **dpd_block_matrix_sp(size_t n, size_t m);
    void free_dpd_block_sp(float **array, size_t n, size_t m);

    int contract222_sp(dpdfile2<float> *X, dpdfile2<float> *Y, dpdfile2<float> *Z, int target_X, int target_Y, float alpha, float beta);
    int contract442_sp(dpdbuf4<float> *X, dpdbuf4<float> *Y, dpdfile2<float> *Z, int target_X, int target_Y, float alpha, float beta);
    int contract422_sp(dpdbuf4<float>*X, dpdfile2<float> *Y, dpdfile2<float> *Z, int trans_Y, int trans_Z, float alpha, float beta);
    int contract244_sp(dpdfile2<float> *X, dpdbuf4<float> *Y, dpdbuf4<float> *Z, int sum_X, int sum_Y, int trans_Z, float alpha, float beta);
    int contract424_sp(dpdbuf4<float> *X, dpdfile2<float> *Y, dpdbuf4<float> *Z, int sum_X, int sum_Y, int trans_Z, float alpha, float beta);
    int contract444_sp(dpdbuf4<float> *X, dpdbuf4<float> *Y, dpdbuf4<float> *Z, int target_X, int target_Y, float alpha, float beta);

    int contract222_mp(dpdfile2<float> *X, dpdfile2<float> *Y, dpdfile2<double> *Z, int target_X, int target_Y, float alpha, double beta);
    int contract442_mp(dpdbuf4<float> *X, dpdbuf4<float> *Y, dpdfile2<double> *Z, int target_X, int target_Y, float alpha, double beta);
    int contract422_mp(dpdbuf4<float> *X, dpdfile2<float> *Y, dpdfile2<double> *Z, int trans_Y, int trans_Z, float alpha, double beta);
    int contract244_mp(dpdfile2<float> *X, dpdbuf4<float> *Y, dpdbuf4<double> *Z, int sum_X, int sum_Y, int trans_Z, float alpha, double beta);
    int contract424_mp(dpdbuf4<float> *X, dpdfile2<float> *Y, dpdbuf4<double> *Z, int sum_X, int sum_Y, int trans_Z, float alpha, double beta);
    int contract444_mp(dpdbuf4<float> *X, dpdbuf4<float> *Y, dpdbuf4<double> *Z, int target_X, int target_Y, float alpha, double beta);

    // not used
    // int contract444_df(dpdbuf4<double> *B, dpdbuf4<double> *tau_in, dpdbuf4<double> *tau_out, double alpha, double beta);

    /* Need to consolidate these routines into one general function */
    int dot23_sp(dpdfile2<float> *T, dpdbuf4<float> *I, dpdfile2<float> *Z, int transt, int transz, float alpha, float beta);
    int dot24_sp(dpdfile2<float> *T, dpdbuf4<float> *I, dpdfile2<float> *Z, int transt, int transz, float alpha, float beta);
    int dot13_sp(dpdfile2<float> *T, dpdbuf4<float> *I, dpdfile2<float> *Z, int transt, int transz, float alpha, float beta);
    int dot14_sp(dpdfile2<float> *T, dpdbuf4<float> *I, dpdfile2<float> *Z, int transt, int transz, float alpha, float beta);

    int dot23_mp(dpdfile2<float> *T, dpdbuf4<float> *I, dpdfile2<double> *Z, int transt, int transz, float alpha, double beta);
    int dot24_mp(dpdfile2<float> *T, dpdbuf4<float> *I, dpdfile2<double> *Z, int transt, int transz, float alpha, double beta);
    int dot13_mp(dpdfile2<float> *T, dpdbuf4<float> *I, dpdfile2<double> *Z, int transt, int transz, float alpha, double beta);
    int dot14_mp(dpdfile2<float> *T, dpdbuf4<float> *I, dpdfile2<double> *Z, int transt, int transz, float alpha, double beta);



    //int trace42_13(dpdbuf4<double> *A, dpdfile2<double> *B, int transb, double alpha, double beta);

    int file2_init_sp(dpdfile2<float> *File, int filenum, int irrep, int pnum, int qnum, const char *label);
    int file2_close_sp(dpdfile2<float> *File);
    int file2_mat_init_sp(dpdfile2<float> *File);
    int file2_mat_close_sp(dpdfile2<float> *File);
    int file2_mat_rd_sp(dpdfile2<float> *File);
   // int file2_mat_wrt(dpdfile2<U> *File);
//**
    // int file2_print(dpdfile2<float> *File, std::string out_fname);
    // int file2_mat_print(dpdfile2<float> *File, std::string out_fname);
    //int file2_copy(dpdfile2<U> *InFile, int outfilenum, const char *label);
   // int file2_dirprd(dpdfile2<U> *FileA, dpdfile2<U> *FileB);
    double file2_dot_sp(dpdfile2<float> *FileA, dpdfile2<float> *FileB);
    int file2_scm_sp(dpdfile2<float> *InFile, float alpha);
//**   
    double file2_dot_self_sp(dpdfile2<float> *BufX);
    // double file2_trace(dpdfile2<U> *InFile);
    //int file2_axpy(dpdfile2<U> *FileA, dpdfile2<U> *FileB, W alpha, int transA);
    //int file2_axpbycz(dpdfile2<U> *FileA, dpdfile2<U> *FileB, dpdfile2<U> *FileC, W a, W b, W c);

/*
    int file4_init(dpdfile4<double> *File, int filenum, int irrep, int pqnum, int rsnum, const char *label);
    int file4_init_nocache(dpdfile4<double> *File, int filenum, int irrep, int pqnum, int rsnum, const char *label);
    int file4_close(dpdfile4<double> *File);
    int file4_mat_irrep_init(dpdfile4<double> *File, int irrep);
    int file4_mat_irrep_close(dpdfile4<double> *File, int irrep);
    int file4_mat_irrep_rd(dpdfile4<double> *File, int irrep);
    int file4_mat_irrep_wrt(dpdfile4<double> *File, int irrep);
    int file4_mat_irrep_row_init(dpdfile4<double> *File, int irrep);
    int file4_mat_irrep_row_close(dpdfile4<double> *File, int irrep);
    int file4_mat_irrep_row_rd(dpdfile4<double> *File, int irrep, int row);
    int file4_mat_irrep_row_wrt(dpdfile4<double> *File, int irrep, int row);
    int file4_mat_irrep_row_zero(dpdfile4<double> *File, int irrep, int row);
    int file4_print(dpdfile4<double> *File, std::string out_fname);
    int file4_mat_irrep_rd_block(dpdfile4<double> *File, int irrep, int start_pq, int num_pq);
    int file4_mat_irrep_wrt_block(dpdfile4<double> *File, int irrep, int start_pq, int num_pq);
*/
    int buf4_init_sp(dpdbuf4<float> *Buf, int inputfile, int irrep, int pqnum, int rsnum, int file_pqnum, int file_rsnum,
                  int anti, const char *label);
    int buf4_init_sp(dpdbuf4<float> *Buf, int inputfile, int irrep, std::string pq, std::string rs, std::string file_pq,
                  std::string file_rs, int anti, const char *label);
    int buf4_init_sp(dpdbuf4<float> *Buf, int inputfile, int irrep, std::string pq, std::string rs, int anti, const char *label);
    //int pairnum(std::string);
//*
    //double buf4_trace(dpdbuf4<U2> *Buf);
    int buf4_close_sp(dpdbuf4<float> *Buf);
    int buf4_mat_irrep_init_sp(dpdbuf4<float> *Buf, int irrep);
    int buf4_mat_irrep_close_sp(dpdbuf4<float> *Buf, int irrep);
    int buf4_mat_irrep_rd_sp(dpdbuf4<float> *Buf, int irrep);
    int buf4_mat_irrep_wrt_sp(dpdbuf4<float> *Buf, int irrep);
//**
    int buf4_print_sp(dpdbuf4<float> *Buf, std::string out_fname, int print_data);
    int buf4_copy_sp(dpdbuf4<float> *InBuf, int outfilenum, const char *label);
    int buf4_sort_sp(dpdbuf4<float> *InBuf, int outfilenum, enum indices index, int pqnum, int rsnum, const char *label);
    int buf4_sort_sp(dpdbuf4<float> *InBuf, int outfilenum, enum indices index, std::string pq, std::string rs,
                  const char *label);
    //int buf4_sort_ooc_sp(dpdbuf4<float> *InBuf, int outfilenum, enum indices index, int pqnum, int rsnum, const char *label);
    int buf4_sort_axpy_sp(dpdbuf4<float> *InBuf, int outfilenum, enum indices index, int pqnum, int rsnum, const char *label,
                       float alpha);
//**
    int buf4_axpy_sp(dpdbuf4<float> *BufX, dpdbuf4<float> *BufY, float alpha);
    int buf4_axpbycz_sp(dpdbuf4<float> *FileA, dpdbuf4<float> *FileB, dpdbuf4<float> *FileC, float a, float b, float c);
//*
    int buf4_dirprd_sp(dpdbuf4<float> *BufA, dpdbuf4<float> *BufB);
    double buf4_dot_sp(dpdbuf4<float> *BufA, dpdbuf4<float> *BufB);
    double buf4_dot_self_sp(dpdbuf4<float> *BufX);
//**
    int buf4_scm_sp(dpdbuf4<float> *InBuf, float alpha);
    int buf4_scmcopy_sp(dpdbuf4<float> *InBuf, int outfilenum, const char *label, float alpha);
   // int buf4_symm(dpdbuf4<double> *Buf);
   // int buf4_symm2(dpdbuf4<double> *Buf1, dpdbuf4<double> *Buf2);
    int buf4_mat_irrep_shift13_sp(dpdbuf4<float> *Buf, int irrep);
    int buf4_mat_irrep_shift31_sp(dpdbuf4<float> *Buf, int irrep);
    int buf4_mat_irrep_row_init_sp(dpdbuf4<float> *Buf, int irrep);
    int buf4_mat_irrep_row_close_sp(dpdbuf4<float> *Buf, int irrep);
    int buf4_mat_irrep_row_zero_sp(dpdbuf4<float> *Buf, int irrep, int row);
    int buf4_mat_irrep_row_rd_sp(dpdbuf4<float> *Buf, int irrep, int pq);
    int buf4_mat_irrep_row_wrt_sp(dpdbuf4<float> *Buf, int irrep, int pq);
//**
    int buf4_mat_irrep_init_block_sp(dpdbuf4<float> *Buf, int irrep, int num_pq);
    int buf4_mat_irrep_close_block_sp(dpdbuf4<float> *Buf, int irrep, int num_pq);
    int buf4_mat_irrep_rd_block_sp(dpdbuf4<float> *Buf, int irrep, int start_pq, int num_pq);
    int buf4_mat_irrep_wrt_block_sp(dpdbuf4<float> *Buf, int irrep, int start_pq, int num_pq);
//*
    //int buf4_dump(dpdbuf4<double> *DPDBuf, struct iwlbuf *IWLBuf, int *prel, int *qrel, int *rrel, int *srel, int bk_pack,
             //     int swap23);
    int trans4_init_sp(dpdtrans4<float> *Trans, dpdbuf4<float> *Buf);
    int trans4_close_sp(dpdtrans4<float> *Trans);
    int trans4_mat_irrep_init_sp(dpdtrans4<float> *Trans, int irrep);
    int trans4_mat_irrep_close_sp(dpdtrans4<float> *Trans, int irrep);
    int trans4_mat_irrep_rd_sp(dpdtrans4<float> *Trans, int irrep);
    int trans4_mat_irrep_wrt_sp(dpdtrans4<float> *Trans, int irrep);
//*
    int trans4_mat_irrep_shift13_sp(dpdtrans4<float> *Trans, int irrep);
    int trans4_mat_irrep_shift31_sp(dpdtrans4<float> *Trans, int irrep);

    //int mat4_irrep_print(double **matrix, dpdparams4 *Params, int irrep, int my_irrep, std::string out_fname);



};  // Dpd class

/*
 * Static variables/functions to mimic the old C machinery
 */
extern dpd_gbl dpd_main;
extern PSI_API DPD *global_dpd_;
extern PSI_API int dpd_default;
extern DPD *dpd_list[2];
extern PSI_API int dpd_set_default(int dpd_num);
extern int dpd_init(int dpd_num, int nirreps, long int memory, int cachetype, int *cachefiles, int **cachelist,
                    dpd_file4_cache_entry *priority, int num_subspaces, std::vector<int *> &spaceArrays);
extern int dpd_close(int dpd_num);
extern long int PSI_API dpd_memfree();
extern void dpd_memset(long int memory);

}  // Namespace psi

#endif /* _psi_src_lib_libdpd_dpd_h */

