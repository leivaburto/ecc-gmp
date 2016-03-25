#ifndef __EC_H__
#define __EC_H__
#include "fp.h"

typedef struct { mpz_t p, a; int max_bits, max_str; } ecc_curve;
typedef struct { mpz_t x, y, z; } ecc_jcb_t;
typedef struct { mpz_t x, y; } ecc_afn_t;

void ecc_init(mpz_ptr, char *);
void ecc_clear(mpz_ptr);
void ecc_init_curve(ecc_curve *, char *, char *, int, int);
void ecc_clear_curve(ecc_curve *);
void ecc_init_jcb(ecc_jcb_t *);
void ecc_init_setstr_jcb(ecc_jcb_t *, char *, char *, char *);
void ecc_init_set_jcb(ecc_jcb_t *, ecc_jcb_t *);
void ecc_cp_jcb(ecc_jcb_t *, ecc_jcb_t *);
void ecc_clear_jcb(ecc_jcb_t *);
void ecc_init_afn(ecc_afn_t *);
void ecc_init_setstr_afn(ecc_afn_t *, char *, char *);
void ecc_init_set_afn(ecc_afn_t *, ecc_afn_t *);
void ecc_clear_afn(ecc_afn_t *);
void ecc_jcb_to_afn(ecc_afn_t *, ecc_jcb_t *, ecc_curve *);
void ecc_get_add_inv(ecc_jcb_t *, ecc_jcb_t *, ecc_curve *);
void ecc_doubling(ecc_jcb_t *, ecc_jcb_t *, ecc_curve *);
void ecc_add(ecc_jcb_t *, ecc_jcb_t *, ecc_jcb_t *, ecc_curve *);
void ecc_add_mix(ecc_jcb_t *, ecc_jcb_t *, ecc_afn_t *, ecc_curve *);

#endif
