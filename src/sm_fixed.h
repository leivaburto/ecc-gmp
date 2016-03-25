#ifndef __SM_FIXED_H__
#define __SM_FIXED_H__
#include <string.h>
#include <stdio.h>
#include <math.h>
#include "ec.h"

#define BGMW_WINDOW 4

#define COMB_MAX 256
#define COMB_WINDOW 8 /* window width */

void sm_pre_bgmw(ecc_jcb_t (*ROP), ecc_jcb_t *, char *, ecc_curve *);
void sm_post_bgmw(ecc_jcb_t (*ROP), char *);
void sm_bgmw(ecc_jcb_t *ROP, ecc_jcb_t (*Pi), ecc_jcb_t *, char *, ecc_curve *);

void sm_pre_comb(int, int, ecc_jcb_t *, ecc_curve *, ecc_jcb_t (*Pt), ecc_jcb_t (*a));
void sm_post_comb(ecc_jcb_t (*Pt), ecc_jcb_t (*a));
void sm_comb(ecc_jcb_t *, int, int, ecc_jcb_t *, mpz_ptr, ecc_curve *, int K[], ecc_jcb_t (*Pt));

#endif
