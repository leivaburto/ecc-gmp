#include "fp.h"

void fp_add_mpz(mpz_ptr rop, mpz_ptr u, mpz_ptr v, mpz_ptr p)
{
	mpz_add(rop, u, v);
	mpz_mod(rop, rop, p);
}

void fp_sub_mpz(mpz_ptr rop, mpz_ptr u, mpz_ptr v, mpz_ptr p)
{
	mpz_sub(rop, u, v);
	mpz_mod(rop, rop, p);
}

void fp_mul_mpz(mpz_ptr rop, mpz_ptr u, mpz_ptr v, mpz_ptr p)
{
	mpz_mul(rop, u, v);
	mpz_mod(rop, rop, p);
}

void fp_pow_si(mpz_ptr rop, mpz_ptr u, long v, mpz_ptr p)
{
    mpz_powm_ui(rop, u, v, p);
}

void fp_mul_si(mpz_ptr rop, mpz_ptr u, long v, mpz_ptr p)
{
	mpz_mul_si(rop, u, v);
	mpz_mod(rop, rop, p);
}
