#ifndef FP_H
#define FP_H
#include <gmp.h>

#ifndef CUSTOM_CURVE
#define CUSTOM_CURVE
#define CURVE_P         "115792089210356248762697446949407573530086143415290314195533631308867097853951"
#define CURVE_A         "115792089210356248762697446949407573530086143415290314195533631308867097853948"
#define CURVE_MAX_BITS 	256
#define CURVE_MAX_STR   80
#endif

void fp_add_mpz(mpz_ptr, mpz_ptr, mpz_ptr, mpz_ptr);
void fp_sub_mpz(mpz_ptr, mpz_ptr, mpz_ptr, mpz_ptr);
void fp_mul_mpz(mpz_ptr, mpz_ptr, mpz_ptr, mpz_ptr);
void fp_pow_si(mpz_ptr, mpz_ptr, long, mpz_ptr);
void fp_mul_si(mpz_ptr, mpz_ptr, long, mpz_ptr);

#endif
