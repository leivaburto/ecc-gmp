#ifndef __NAF_H__
#define __NAF_H__
#include <math.h>
#include "fp.h"

#define NAF_WINDOW 5

void naf_init_bits(mpz_ptr, int);
void naf_set_bit(mpz_ptr, int, short);
short naf_get_bit(mpz_ptr, int);
void naf_get_substr(mpz_ptr, char[], int, int);
void naf_convert(mpz_ptr, mpz_ptr, int *);
int naf_convert_inverse(char[], int);

void wnaf_init_bits(mpz_ptr, int);
void wnaf_mods(int *, mpz_ptr);
void wnaf_set_bit(mpz_ptr, int, short);
short wnaf_get_bit(mpz_ptr, int);
void wnaf_convert(mpz_ptr, mpz_ptr, int *);

#endif
