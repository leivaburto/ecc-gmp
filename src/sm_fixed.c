#include "sm_fixed.h"

/* Precomputation of d points, P_{i} = 2^{w*i}*P */
/* Argument k must be represented in 2^w base, where w is size of window */
/* Save result in an array of points */
void sm_pre_bgmw(ecc_jcb_t (*ROP), ecc_jcb_t *P, char *k, ecc_curve *curve)
{	
	int i, j, d;
	ecc_jcb_t Q;
	ecc_init_jcb(&Q);
	/* Size of k in base 2^w  */
	d = strlen(k);
    /* 2^w*0 => P */
    ecc_init_set_jcb(&ROP[0], P);
    /* 2^{w*i}*P*/
	for(i = 1; i < d; i++)
	{
        /* w dublings over the previous value */
        ecc_init_set_jcb(&ROP[i], &ROP[i-1]);
        for(j = 0; j < BGMW_WINDOW; j++)
        {
            ecc_doubling(&Q, &ROP[i], curve);
            ecc_cp_jcb(&ROP[i], &Q);
        }
	}
	ecc_clear_jcb(&Q);
}

void sm_post_bgmw(ecc_jcb_t (*ROP), char *k)
{
	int i, d;
	d = strlen(k);
	for(i = 0; i < d; i++)
		ecc_clear_jcb(&ROP[i]);
}

void sm_bgmw(ecc_jcb_t *ROP, ecc_jcb_t (*Pi), ecc_jcb_t *P, char *k, ecc_curve *curve)
{
	int i, j, limit, d, dif;
	short flag, flag_sum;
	ecc_jcb_t B, Q;
	ecc_init_jcb(&B);
	ecc_init_jcb(&Q);

	d = strlen(k);
	limit = 1 << BGMW_WINDOW; /* 2^w */
	flag = 0;
	for(j = limit-1; j > 0; j--)
	{
		for(i = 0; i < d; i++)
		{
            /* Deal with hexadecimal (and greater) base */
            if(isdigit(k[i])) dif = k[i]-'0'; else dif = k[i]-'a'+10; 
			if(dif == j)
			{
				if(!flag) /* First iteration */
				{
                    ecc_cp_jcb(&B, &Pi[d-i-1]);
					flag = 1;
				}
				else
				{
					ecc_add(&Q, &B, &Pi[d-i-1], curve);
					ecc_cp_jcb(&B, &Q);
				}
			}
		}
		if(flag == 1)
		{
            ecc_cp_jcb(ROP, &B);
			flag = 2;
		}
		else if(flag == 2)
		{
			ecc_add(&Q, ROP, &B, curve);
			ecc_cp_jcb(ROP, &Q);
		}
	}
	ecc_clear_jcb(&B);
	ecc_clear_jcb(&Q);
}


void sm_pre_comb(int d, int e, ecc_jcb_t *P, ecc_curve *curve, ecc_jcb_t (*Pt), ecc_jcb_t (*a))
{
    short flag = 0;
    long long i, j;
    ecc_jcb_t ROP, ROP2, zero;

    ecc_init_jcb(&ROP);
    ecc_init_jcb(&ROP2);
    ecc_init_jcb(&zero);
    ecc_init_set_jcb(&a[0], P);

    for(i = 1; i < COMB_WINDOW; i++)
    {
        ecc_cp_jcb(&ROP, P);
        for(j = 1; j <= d*i; j++)
        {
            ecc_doubling(&ROP2, &ROP, curve);
            ecc_cp_jcb(&ROP, &ROP2);
        }
        ecc_init_set_jcb(&a[i], &ROP);
    }

    for(i = 0; i < (1 << COMB_WINDOW); i++)
    {
        ecc_cp_jcb(&ROP, &zero);
        for(j = COMB_WINDOW-1; j >= 0; j--)
        {
            if(i & (1 << j))
            {
                if(!flag)
                {
                    ecc_cp_jcb(&ROP, &a[j]);
                    flag = 1;
                }
                else
                {
                    ecc_add(&ROP2, &a[j], &ROP, curve);
                    ecc_cp_jcb(&ROP, &ROP2);
                }
            }
        }

        ecc_init_set_jcb(&Pt[i], &ROP);
    }

	ecc_clear_jcb(&ROP);
    ecc_clear_jcb(&ROP2);
    ecc_clear_jcb(&zero);
}

void sm_post_comb(ecc_jcb_t (*Pt), ecc_jcb_t (*a))
{
    long long i;
    for(i = 0; i < (1 << COMB_WINDOW); i++)
        ecc_clear_jcb(&Pt[i]);
    
    for(i = 0; i < COMB_WINDOW; i++)
        ecc_clear_jcb(&a[i]);
}

void sm_comb(ecc_jcb_t *ROP, int d, int e, ecc_jcb_t *P, mpz_ptr k, ecc_curve *curve, int K[COMB_MAX], ecc_jcb_t (*Pt))
{
    long long i, j, tmp;
    ecc_jcb_t ROP2;

    ecc_init_jcb(&ROP2);

    for(i = d-1; i >= 0; i--)
    {
        tmp = 0;
        for(j = COMB_WINDOW-1; j >= 0; j--)
        {
            if(mpz_tstbit(k, d*j + i))
            {
                tmp = tmp | 1 << j; 
            }
        }
        K[d-1-i] = tmp;
    }

    ecc_cp_jcb(ROP, &Pt[K[0]]);
    for(i = 1; i < d; i++)
    {
        ecc_doubling(&ROP2, ROP, curve);
        ecc_cp_jcb(ROP, &ROP2);

        ecc_add(&ROP2, ROP, &Pt[K[i]], curve);
        ecc_cp_jcb(ROP, &ROP2);
    }

    ecc_clear_jcb(&ROP2);
}
