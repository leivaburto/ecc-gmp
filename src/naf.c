#include <stdio.h>
#include "naf.h"

void naf_init_bits(mpz_ptr k, int max_bits)
{
	mpz_init2(k, max_bits*2);
}

void naf_set_bit(mpz_ptr k_naf, int bit, short val)
{
	bit <<= 2; // 2*bit
	if(val == 1)
		mpz_setbit(k_naf, bit);
	else if(val == -1) // -1 => {11}
	{
		mpz_setbit(k_naf, bit);
		mpz_setbit(k_naf, bit+1);
	}
}

short naf_get_bit(mpz_ptr naf_k, int bit)
{
	bit <<= 2; // 2*bit
	if(mpz_tstbit(naf_k, bit))
	{
		if(mpz_tstbit(naf_k, bit+1))
			return -1;
		else
			return 1;
	}
	else
		return 0;
}

void naf_get_substr(mpz_ptr naf_k, char naf_subk[], int bgn, int end)
{
    int i, v = 0;
    for(i = 0; i < end; i++)
    {
        if(naf_get_bit(naf_k, bgn-i) == 1)
            naf_subk[i] = '1';
        else if(naf_get_bit(naf_k, bgn-i) == -1)
            naf_subk[i] = '-';
        else
            naf_subk[i] = '0';
    }
    naf_subk[end] = '\0';
}

void naf_convert(mpz_ptr k_naf, mpz_ptr k, int *size)
{
	int i = 0;
	mpz_t kmod, ki;
	mpz_inits(kmod, ki, NULL);

	while(mpz_sgn(k)) // k > 0
	{
		if(mpz_odd_p(k))
		{
			mpz_tdiv_r_ui(kmod, k, 4); // kmod = k%4
			mpz_ui_sub(ki, 2, kmod);
			mpz_sub(k, k, ki);
		}
		else
			mpz_set_ui(ki, 0);
		
		mpz_cdiv_q_2exp(k, k, 1); // k = k/(2^1) <- right shifts
		naf_set_bit(k_naf, i, mpz_get_si(ki)); // set i-th bit (actually, i and i+1)
		i++;
	}
	*size = i;

	mpz_clears(kmod, ki, NULL);
}

int naf_convert_inverse(char naf_str[], int size)
{
    int i, v = 0;
    for(i = size-1; i >= 0; i--)
    {
        if(naf_str[i] == '-')
            v = v - (1 << size-i-1);
        else if(naf_str[i] == '1')
            v = v + (1 << size-i-1);
    }
    return v;
}


void naf_print(mpz_ptr naf_k, int limit)
{
	int i;
	for(i = limit-1; i >= 0; i--)
		printf("%d", naf_get_bit(naf_k, i));
	printf("\n");
}


void wnaf_init_bits(mpz_ptr k, int max_bits)
{
	/* |W| bits para cada número...
	 * Primer bit para signo y W-1 restantes para almacenar el número */
	mpz_init2(k, max_bits*NAF_WINDOW);
}

void wnaf_mods(int *val, mpz_ptr k)
{
	short limit = 1;
	mpz_t kmod;
	mpz_init(kmod);

	limit <<= NAF_WINDOW; // 2^W
	mpz_tdiv_r_ui(kmod, k, limit); // kmod = kmod2^W

	if(mpz_cmp_ui(kmod, (limit >> 1)) >= 0) // kmod >= 2^(W-1)
		mpz_sub_ui(kmod, kmod, limit);
	
	*val = mpz_get_si(kmod);

	mpz_clear(kmod);
}

void wnaf_set_bit(mpz_ptr k_naf, int bit_pos, short val)
{
	short bit_val, bit_count;
	bit_pos <<= NAF_WINDOW; // WINDOW_SIZE*bit
	if(val < 0) // Activar primer bit para identificar valores negativos
	{
		mpz_setbit(k_naf, bit_pos);
		val = -val; // Ahora trabajar con valor positivo
	}
	// Convertir val to binary string
	bit_val = 0;
	bit_count = 1;
	while(val > 0)
	{
		bit_val = val%2;
		if(bit_val)
			mpz_setbit(k_naf, bit_pos + bit_count);
		val >>= 1; // val /= 2
		bit_count++;
	}
}

short wnaf_get_bit(mpz_ptr naf_k, int bit_pos)
{
	short bit_count, val;
	// Empezar a comprobar bits de derecha a izq (desde extremo derecho de la ventana)
	bit_pos <<= NAF_WINDOW;
	val = 0;
	for(bit_count = 1; bit_count <= NAF_WINDOW; bit_count++)
		if(mpz_tstbit(naf_k, bit_pos + bit_count))
			val += (1 << (bit_count-1));

	// Comprobar signo al extremo izquierdo de la ventana
	if(mpz_tstbit(naf_k, bit_pos))
		val *= -1;
	return val;
}

void wnaf_convert(mpz_ptr k_naf, mpz_ptr k, int *size)
{
	int i = 0, ki;

	while(mpz_sgn(k)) // k > 0
	{
		if(mpz_odd_p(k))
		{
			wnaf_mods(&ki, k); // ki = k mods 2^W
			// k = k - ki
			// GMP no tiene signed integers operations, así que que:
			// Si ki es negativo, k = k + (-1)ki
			// Si ki es positivo, k = k - ki
			if(ki > 0)
				mpz_sub_ui(k, k, (unsigned long int)ki);
			else
				mpz_add_ui(k, k, (unsigned long int)(ki*-1));
		}
		else
			ki = 0;
		
		mpz_cdiv_q_2exp(k, k, 1); // k = k/(2^1) <- TODO : CHANGE TO right shifts
		wnaf_set_bit(k_naf, i, (short)ki); // set "i-th bit"
		i++;
	}
	*size = i;
}

/*
void wnaf_print(mpz_ptr naf_k, int limit)
{
	int i;
	// Imprimir de izq a derecha
	for(i = limit-1; i >= 0; i--)
	{
		printf("%d", wnaf_get_bit(naf_k, i));
	}
	printf("\n");
}
*/
