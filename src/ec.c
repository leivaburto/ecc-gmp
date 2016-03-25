#include "ec.h"
#include <stdio.h>
#include <string.h>

void ecc_init(mpz_ptr p, char *p_str)
{
	mpz_init_set_str(p, p_str, 10);
}

void ecc_clear(mpz_ptr p)
{
	mpz_clear(p);
}

void ecc_init_curve(ecc_curve *c, char *p, char *a, int max_bits, int max_str)
{
    c->max_bits  = max_bits;
    c->max_str   = max_str;
    mpz_init_set_str(c->p, p, 10);
    mpz_init_set_str(c->a, a, 10);
}

void ecc_clear_curve(ecc_curve *c)
{
    mpz_clears(c->p, c->a, NULL);
}

void ecc_init_jcb(ecc_jcb_t *P)
{
	mpz_inits(P->x, P->y, P->z, NULL);
}

void ecc_init_setstr_jcb(ecc_jcb_t *P, char *x, char *y, char *z)
{
	mpz_inits(P->x, P->y, P->z, NULL);

	mpz_set_str(P->x, x, 10);
	mpz_set_str(P->y, y, 10);
	mpz_set_str(P->z, z, 10);
}

void ecc_init_set_jcb(ecc_jcb_t *P, ecc_jcb_t *Q)
{
	mpz_inits(P->x, P->y, P->z, NULL);

	mpz_set(P->x, Q->x);
	mpz_set(P->y, Q->y);
	mpz_set(P->z, Q->z);
}

/**/
void ecc_cp_jcb(ecc_jcb_t *ROP, ecc_jcb_t *P)
{
	mpz_set(ROP->x, P->x);
	mpz_set(ROP->y, P->y);
	mpz_set(ROP->z, P->z);
}
/**/

void ecc_clear_jcb(ecc_jcb_t *P)
{
	mpz_clears(P->x, P->y, P->z, NULL);
}

void ecc_init_afn(ecc_afn_t *P)
{
	mpz_inits(P->x, P->y, NULL);
}

void ecc_init_set_afn(ecc_afn_t *P, ecc_afn_t *Q)
{
	mpz_inits(P->x, P->y, NULL);

	mpz_set(P->x, Q->x);
	mpz_set(P->y, Q->y);
}

void ecc_init_setstr_afn(ecc_afn_t *P, char *x, char *y)
{
	mpz_inits(P->x, P->y, NULL);

	mpz_set_str(P->x, x, 10);
	mpz_set_str(P->y, y, 10);
}

void ecc_clear_afn(ecc_afn_t *P)
{
	mpz_clears(P->x, P->y, NULL);
}

void ecc_jcb_to_afn(ecc_afn_t *ROP, ecc_jcb_t *P, ecc_curve *curve)
{
	mpz_t g, s, inv, tmp;	
	mpz_inits(g, s, inv, tmp, NULL);

	mpz_gcdext(g, s, inv, curve->p, P->z); // inv = Z^1

	fp_mul_mpz(tmp, inv, inv, curve->p); // (Z^1)^2
	fp_mul_mpz(ROP->x, P->x, tmp, curve->p); // X*((Z^1)^2)

	fp_mul_mpz(tmp, tmp, inv, curve->p); // (Z^1)^3
	fp_mul_mpz(ROP->y, P->y, tmp, curve->p); // Y*((Z^1)^3)

	mpz_clears(g, s, inv, tmp, NULL);
}

void ecc_get_add_inv(ecc_jcb_t *ROP, ecc_jcb_t *P, ecc_curve *curve)
{
	/* -P = P(x, -y, z) */
	mpz_set(ROP->x, P->x);
	fp_mul_si(ROP->y, P->y, -1, curve->p);
	mpz_set(ROP->z, P->z);
}

void ecc_doubling(ecc_jcb_t *ROP, ecc_jcb_t *P, ecc_curve *curve)
{
    mpz_t alfa, alfa1, alfa2, beta, x3_2, y3_2, y2, z3_2;
    if((strcmp(mpz_get_str(NULL, 10, P->x), "0") == 0) && (strcmp(mpz_get_str(NULL, 10, P->y), "0") == 0) && (strcmp(mpz_get_str(NULL, 10, P->z), "0") == 0))
    {
        ecc_cp_jcb(ROP, P);
    }
    else
    {
        mpz_inits(alfa, alfa1, alfa2, beta, x3_2, y3_2, z3_2, y2, NULL);
        // alfa = 3*(x1)**2 + a * (z1**4)
        fp_mul_mpz(alfa1, P->x, P->x, curve->p);
        fp_mul_mpz(alfa2, P->z, P->z, curve->p);
        fp_mul_mpz(alfa2, alfa2, alfa2, curve->p);
        fp_mul_si(alfa1, alfa1, 3, curve->p);
        fp_mul_mpz(alfa2, alfa2, curve->a, curve->p);
        fp_add_mpz(alfa, alfa1, alfa2, curve->p);

        // y2 = y1**2
        fp_mul_mpz(y2, P->y, P->y, curve->p);
        // beta = 4*x1*(y1**2)
        fp_mul_mpz(beta, y2, P->x, curve->p);
        fp_mul_si(beta, beta, 4, curve->p);
        
        mpz_set(ROP->z, P->z);
        
        // z3 = y1*z1
        fp_mul_mpz(ROP->z, P->z, P->y, curve->p);
        
        // z3 = 2*y1*z1
        fp_mul_si(ROP->z, ROP->z, 2, curve->p);
        
        mpz_init_set(x3_2, beta);
        
        // x3_2 = 2*beta
        fp_mul_si(x3_2, x3_2, 2, curve->p);
        mpz_set(ROP->x, alfa);
        
        // x3 = alfa**2 - 2*beta
        fp_mul_mpz(ROP->x, ROP->x, ROP->x, curve->p);
        fp_sub_mpz(ROP->x, ROP->x, x3_2, curve->p);
        
        mpz_init_set(y3_2, P->y);
        mpz_set(ROP->y, beta);
        
        // y3 = alfa(beta -x3) 
        fp_sub_mpz(ROP->y, ROP->y, ROP->x, curve->p);
        fp_mul_mpz(ROP->y, ROP->y, alfa, curve->p);
        
        // y3_2 = 8*y1**4
        fp_mul_mpz(y3_2, y2, y2, curve->p);
        fp_mul_si(y3_2, y3_2, 8, curve->p);
        
        // y3 = alfa(beta -x3) - 8*y1**4
        fp_sub_mpz(ROP->y, ROP->y, y3_2, curve->p);
        
        mpz_clears(alfa, alfa1, alfa2, beta, x3_2, y3_2, z3_2, y2, NULL);
    }
}

void ecc_add(ecc_jcb_t *ROP, ecc_jcb_t *P, ecc_jcb_t *Q, ecc_curve *curve)
{
    mpz_t alpha, beta, beta_pow2, Z2_pow2, tmp1, tmp2, g, s;
    mpz_inits(alpha, beta, beta_pow2, Z2_pow2, tmp1, tmp2, g, s, NULL);
    
    if((mpz_cmp(P->x, Q->x) == 0) && (mpz_cmp(P->y, Q->y) == 0) && (mpz_cmp(P->z, Q->z) == 0))
    {
		ecc_doubling(ROP, P, curve);
	}
    else if((mpz_sgn(P->x) == 0) && (mpz_sgn(P->y) == 0) && (mpz_sgn(P->z) == 0))
    {
        ecc_cp_jcb(ROP, Q);
    }
    else if((mpz_sgn(Q->x) == 0) && (mpz_sgn(Q->y) == 0) && (mpz_sgn(Q->z) == 0))
    {
        ecc_cp_jcb(ROP, P);
    }
	else
	{
		/* Z2^2 */    
		fp_mul_mpz(Z2_pow2, Q->z, Q->z, curve->p); // Z1^2
		
		/* alpha = Z1^3*Y2 - Z2^3*Y1 */
		fp_pow_si(alpha, P->z, 3, curve->p);
		fp_mul_mpz(alpha, alpha, Q->y, curve->p);
		fp_pow_si(tmp1, Q->z, 3, curve->p);
		fp_mul_mpz(tmp1, tmp1, P->y, curve->p);
		fp_sub_mpz(alpha, alpha, tmp1, curve->p);

		/* beta = Z1^2*X2 - Z2^2*X1 */
		fp_pow_si(beta, P->z, 2, curve->p);
		fp_mul_mpz(beta, beta, Q->x, curve->p);
		fp_pow_si(tmp1, Q->z, 2, curve->p);
		fp_mul_mpz(tmp1, tmp1, P->x, curve->p);
		fp_sub_mpz(beta, beta, tmp1, curve->p);
		
		/* beta^2 */
		fp_pow_si(beta_pow2, beta, 2, curve->p);

		/* X3 = alpha^2 - beta^3 - 2*Z2^2*X1*beta^2 */
		fp_pow_si(tmp1, alpha, 2, curve->p);
		fp_pow_si(tmp2, beta, 3, curve->p);
		fp_sub_mpz(tmp1, tmp1, tmp2, curve->p);
		fp_mul_si(tmp2, Z2_pow2, 2, curve->p);
		fp_mul_mpz(tmp2, tmp2, P->x, curve->p);
		fp_mul_mpz(tmp2, tmp2, beta_pow2, curve->p);
		fp_sub_mpz(ROP->x, tmp1, tmp2, curve->p);

		/* Y3 = alpha*(Z2^2*X1*beta^2 - X3) - Z2^3*Y1*beta^3*/
		fp_mul_mpz(tmp1, Z2_pow2, P->x, curve->p);
		fp_mul_mpz(tmp1, tmp1, beta_pow2, curve->p);
		fp_sub_mpz(tmp1, tmp1, ROP->x, curve->p);
		fp_mul_mpz(tmp1, tmp1, alpha, curve->p);
		fp_mul_mpz(tmp2, Z2_pow2, Q->z, curve->p);
		fp_mul_mpz(tmp2, tmp2, P->y, curve->p);
		fp_mul_mpz(tmp2, tmp2, beta_pow2, curve->p);
		fp_mul_mpz(tmp2, tmp2, beta, curve->p);
		fp_sub_mpz(ROP->y, tmp1, tmp2, curve->p);
			
		/* Z3 = Z1*Z2*beta */
		fp_mul_mpz(tmp1, P->z, Q->z, curve->p);
		fp_mul_mpz(ROP->z, tmp1, beta, curve->p);
	}

    mpz_clears(alpha, beta, beta_pow2, Z2_pow2, tmp1, tmp2, g, s, NULL);
}

void ecc_add_mix(ecc_jcb_t *ROP, ecc_jcb_t *P, ecc_afn_t *Q, ecc_curve *curve)
{
	mpz_t alfa, beta, alfa_sqr, beta_sqr, beta_cub, Z1_sqr, Z1_cub, tmp;
	mpz_inits(alfa, beta, alfa_sqr, beta_sqr, beta_cub, Z1_sqr, Z1_cub, tmp, NULL);

	fp_mul_mpz(Z1_sqr, P->z, P->z, curve->p); // Z1^2
	fp_mul_mpz(Z1_cub, Z1_sqr, P->z, curve->p); // Z1^3

	// beta = (Z1^2)*X2 - X1
	fp_mul_mpz(tmp, Z1_sqr, Q->x, curve->p);
	fp_sub_mpz(beta, tmp, P->x, curve->p);

	// alfa = (Z1^3)*Y2 - Y1
	fp_mul_mpz(tmp, Z1_cub, Q->y, curve->p);
	fp_sub_mpz(alfa, tmp, P->y, curve->p);

	fp_mul_mpz(beta_sqr, beta, beta, curve->p); // beta^2
	fp_mul_mpz(beta_cub, beta_sqr, beta, curve->p); // beta^3

	fp_mul_mpz(alfa_sqr, alfa, alfa, curve->p); // alfa^2

	// X3 = alfa^2 - beta^3 - 2*X1*beta^2
	fp_sub_mpz(ROP->x, alfa_sqr, beta_cub, curve->p);
	fp_mul_si(tmp, P->x, 2, curve->p);
	fp_mul_mpz(tmp, tmp, beta_sqr, curve->p);
	fp_sub_mpz(ROP->x, ROP->x, tmp, curve->p);

	// Y3 = alfa*(X1*beta^2 - X3) - Y1*beta^3
	fp_mul_mpz(tmp, P->x, beta_sqr, curve->p);
	fp_sub_mpz(tmp, tmp, ROP->x, curve->p);
	fp_mul_mpz(tmp, alfa, tmp, curve->p);
	mpz_set(ROP->y, tmp);
	fp_mul_mpz(tmp, P->y, beta_cub, curve->p);
	fp_sub_mpz(ROP->y, ROP->y, tmp, curve->p);

	// Z3 = Z1*beta
	fp_mul_mpz(ROP->z, P->z, beta, curve->p);

	mpz_clears(alfa, beta, alfa_sqr, beta_sqr, beta_cub, Z1_sqr, Z1_cub, tmp, NULL);
}
