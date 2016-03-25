#include "sm.h"

void ltr(ecc_jcb_t *ROP, ecc_afn_t *P, mpz_ptr k, ecc_curve *curve)
{
    int i, l;
    ecc_jcb_t Q;
    l = mpz_sizeinbase(k, 2);    
    /* 
     * First bit is alway 1 (if going from left-to-right),
     * so the first iteration is Q <- P (adding z=1 coord).
     */
    ecc_init_jcb(&Q);
    mpz_set(Q.x, P->x);
    mpz_set(Q.y, P->y);
    mpz_set_str(Q.z, "1", 10); 
    for(i = l-2; i >= 0; i--)
    {
        ecc_doubling(ROP, &Q, curve);
        ecc_cp_jcb(&Q, ROP);
        if(mpz_tstbit(k, i))
        {
            ecc_add_mix(ROP, &Q, P, curve);
            ecc_cp_jcb(&Q, ROP);
        }
    }
}

void rtl(ecc_jcb_t *ROP, ecc_jcb_t *P, mpz_ptr k, ecc_curve *curve)
{
    int i, l, flag = 0;
    ecc_jcb_t Q;
    l = mpz_sizeinbase(k, 2);
    for(i = 0; i < l; i++)
    {
        if(mpz_tstbit(k, i))
        {
			if (!flag)
            {
                ecc_init_set_jcb(&Q, P);
				flag = 1;
			}
			else
            {
				ecc_add(ROP, &Q, P, curve);
                ecc_cp_jcb(&Q, ROP);
			}
        }
        ecc_doubling(ROP, P, curve);
        ecc_cp_jcb(P, ROP);
    }
    ecc_cp_jcb(ROP, &Q); 
}

void sm_bin_naf(ecc_jcb_t *ROP, ecc_jcb_t *P, mpz_ptr k, ecc_curve *curve)
{
    int i, size;
    short k_i;
    ecc_jcb_t Q, _P; /* _P es (-P) */
    mpz_t k_naf;

	mpz_init(k_naf);
    /* Left-to-right => iteration 0: Q = P */
	ecc_init_set_jcb(&Q, P);
	ecc_init_jcb(&_P);
	ecc_get_add_inv(&_P, P, curve);
	/* Obtener representación de k en NAF */
	naf_init_bits(k_naf, curve->max_bits);
	naf_convert(k_naf, k, &size);

    for(i = size-2; i >= 0; i--)
    {
        ecc_doubling(ROP, &Q, curve);
        ecc_cp_jcb(&Q, ROP);
        k_i = naf_get_bit(k_naf, i);
        if(k_i != 0)
        {
            if(k_i == 1) /* 1 || -1 */
                ecc_add(ROP, &Q, P, curve);
            else
                ecc_add(ROP, &Q, &_P, curve);
            ecc_cp_jcb(&Q, ROP);
        }
    }
    ecc_clear_jcb(&_P);
    ecc_clear_jcb(&Q);
	mpz_clear(k_naf);
}

/*
 * Computes iP for i = {1,3,5,...,limit}
 */
void sm_odd_pi(ecc_jcb_t *P, ecc_jcb_t (*Pi), int upper_bound, ecc_curve *curve)
{
    int i;
    ecc_jcb_t P2, ROP;
    ecc_init_jcb(&P2);
    ecc_init_jcb(&ROP);
    /* P */
    ecc_init_jcb(&Pi[0]);
    ecc_init_set_jcb(&Pi[1], P);
    ecc_get_add_inv(&Pi[0], &Pi[1], curve);
    /* 2P */
    ecc_doubling(&P2, P, curve);
    for(i = 3; i <= upper_bound; i = i + 2)
    {
        ecc_init_jcb(&Pi[i-1]);
        ecc_init_jcb(&Pi[i]);
        ecc_add(&Pi[i], &Pi[i-2], &P2, curve);
        ecc_get_add_inv(&Pi[i-1], &Pi[i], curve);
    }
    ecc_clear_jcb(&ROP);
    ecc_clear_jcb(&P2);
}

void sm_clear_odd_pi(ecc_jcb_t (*Pi), int upper_bound)
{
    int i;
	for(i = 0; i < upper_bound; i++)
		ecc_clear_jcb(&Pi[i]);
}


void sm_wnaf(ecc_jcb_t *ROP, ecc_jcb_t *P, mpz_ptr k, ecc_curve *curve)
{
    int i, size, limit, flag = 0;
    short k_i;
    ecc_jcb_t Q, P_i[NAF_UPPER_BOUND_PI];
    mpz_t k_naf;

	mpz_init(k_naf);
    ecc_init_setstr_jcb(&Q, "0", "0", "0");
    /* Precomputo */
    sm_odd_pi(P, P_i, NAF_UPPER_BOUND_PI, curve);
	/* Obtener representación de k en W-NAF  */
	wnaf_init_bits(k_naf, curve->max_bits);
	wnaf_convert(k_naf, k, &size);
    for(i = size-1; i >= 0; i--)
    {
        ecc_doubling(ROP, &Q, curve);
        ecc_cp_jcb(&Q, ROP);
        k_i = wnaf_get_bit(k_naf, i);
        if(k_i != 0)
        {
			if(k_i > 0)
			{
				if(!flag)
				{
                    ecc_cp_jcb(&Q, &P_i[k_i]);
					flag = 1;
				}
				else
				{				
					ecc_add(ROP, &Q, &P_i[k_i], curve);
                    ecc_cp_jcb(&Q, ROP);
				}
			}
			else
			{
				if(!flag)
				{
                    ecc_cp_jcb(&Q, &P_i[-1*k_i-1]);
					flag = 1;
				}
				else
				{
					ecc_add(ROP, &Q, &P_i[-1*k_i-1], curve);
                    ecc_cp_jcb(&Q, ROP);
				}
			}
        }
    }
    ecc_cp_jcb(ROP, &Q);

	sm_clear_odd_pi(P_i, NAF_UPPER_BOUND_PI);
	ecc_clear_jcb(&Q);
	mpz_clear(k_naf);
}

void sm_sw(ecc_jcb_t *ROP, ecc_jcb_t *P, mpz_ptr k, ecc_curve *curve)
{
    int i, j, size, subk_size, t, u, flag = 0;
    short k_i;
    char subk_naf_str[curve->max_str];
    ecc_jcb_t Q, P_i[(int)(SW_UPPER_BOUND_PI)];
    mpz_t k_naf, subk_naf, subk;

	mpz_inits(k_naf, subk_naf, subk, NULL);
	ecc_init_jcb(&Q);

    /* Precomputo */
    sm_odd_pi(P, P_i, (int)(SW_UPPER_BOUND_PI), curve);
	/* Obtener representación de k en NAF */
	naf_init_bits(k_naf, curve->max_bits);
	naf_convert(k_naf, k, &size);
    naf_init_bits(subk_naf, curve->max_bits);

    i = size-1;
    while(i >= 0)
    {
        k_i = naf_get_bit(k_naf, i);
        if(k_i == 0)
        {
            t = 1;
            u = 0;
        }
        else
        {
            t = 1;
            /* Find largest t <= w such that k_i, ..., k_i - t+1 is odd */
            for(j = 1; j < SW_WINDOW; j++)
            {
                k_i = naf_get_bit(k_naf, i-j);
                if(k_i != 0)
                    t = j+1;
            }
            naf_get_substr(k_naf, subk_naf_str, i, t);
            u = naf_convert_inverse(subk_naf_str, t);
        }
        /* Se calcula Q <- Q*2^t i.e. t doublings */
        if(flag)
        {
            for(j = 1; j <= t; j++)
            {
                ecc_doubling(ROP, &Q, curve);
                ecc_cp_jcb(&Q, ROP);
            }
        }
        /* Dependiendo de si u es mayor o menor que cero se suma P_u
         * o se resta P_-u */        
        if(u > 0)
        {
            if(!flag)
            {
                ecc_cp_jcb(&Q, &P_i[u]);
                flag = 1;
            }
            else
            {
                ecc_add(ROP, &Q, &P_i[u], curve);
                ecc_cp_jcb(&Q, ROP);
            }
        }
        else if(u < 0)
        {
            if(!flag)
            {
                ecc_cp_jcb(&Q, &P_i[-1*u-1]);
                flag = 1;
            }
            else
            {
                ecc_add(ROP, &Q, &P_i[-1*u-1], curve);
                ecc_cp_jcb(&Q, ROP);
            }
        }
        i = i-t;
    }
    ecc_cp_jcb(ROP, &Q);
    
    sm_clear_odd_pi(P_i, SW_UPPER_BOUND_PI);
    ecc_clear_jcb(&Q);
	mpz_clears(k_naf, subk_naf, subk, NULL);
}
