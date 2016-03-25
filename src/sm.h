#ifndef __SM_H__
#define __SM_H__
#include <string.h>
#include "naf.h"
#include "ec.h"
#include "utils.h"

#define POS_NEG(x) ((x%2) ? (-1) : (1))

#define SW_WINDOW           5
#define SW_UPPER_BOUND_PI   (( (1 << SW_WINDOW) - POS_NEG(SW_WINDOW) )/3)*2
#define NAF_UPPER_BOUND_PI  1 << (NAF_WINDOW-1)

void ltr(ecc_jcb_t *, ecc_afn_t *, mpz_ptr, ecc_curve *);
void rtl(ecc_jcb_t *, ecc_jcb_t *, mpz_ptr, ecc_curve *);
void sm_bin_naf(ecc_jcb_t *, ecc_jcb_t *, mpz_ptr, ecc_curve *);
void sm_odd_pi(ecc_jcb_t *, ecc_jcb_t (*Pi), int, ecc_curve *);
void sm_wnaf(ecc_jcb_t *, ecc_jcb_t *, mpz_ptr, ecc_curve *);
void sm_sw(ecc_jcb_t *, ecc_jcb_t *, mpz_ptr, ecc_curve *);

#endif
