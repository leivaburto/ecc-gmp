#include "utils.h"

int ecc_read_jcb_from_file(ecc_jcb_t *P, char file_path[], int max_str)
{
    FILE *fp;
    char x[max_str], y[max_str], z[max_str];
    
    fp = fopen(file_path, "r");
    if(fp == NULL)
        return -1;

    fscanf(fp, "%s", x);
    fscanf(fp, "%s", y);
    fscanf(fp, "%s", z);

    ecc_init_setstr_jcb(P, x, y, z);

    if(ferror(fp))
        return -1;

    fclose(fp);
    return 0;
}

int ecc_read_afn_from_file(ecc_afn_t *P, char file_path[], int max_str)
{
    FILE *fp;
    char x[max_str], y[max_str];
    
    fp = fopen(file_path, "r");
    if(fp == NULL)
        return -1;

    fscanf(fp, "%s", x);
    fscanf(fp, "%s", y);

    ecc_init_setstr_afn(P, x, y);

    if(ferror(fp))
        return -1;

    fclose(fp);
    return 0;
}

int ecc_read_k_from_file(mpz_ptr k, char file_path[], int max_str)
{
    FILE *fp;
    char k_str[max_str];
    
    fp = fopen(file_path, "r");
    if(fp == NULL)
        return -1;

    fscanf(fp, "%s", k_str);

    ecc_init(k, k_str);

    if(ferror(fp))
        return -1;

    fclose(fp);
    return 0;
}

void ecc_print_jcb(ecc_jcb_t *P)
{
	printf("(%s : %s : %s)\n", mpz_get_str(NULL, 10, P->x), mpz_get_str(NULL, 10, P->y), mpz_get_str(NULL, 10, P->z));
}

void ecc_print_afn(ecc_afn_t *P)
{
	printf("(%s : %s : 1)\n", mpz_get_str(NULL, 10, P->x), mpz_get_str(NULL, 10, P->y));
}

void diff_time(struct timespec *t_fin, struct timespec *t_ini, struct timespec *delta )
{
    if( ( (*t_fin).tv_nsec - (*t_ini).tv_nsec ) < 0 )
    {
        if( (*t_fin).tv_sec == (*t_ini).tv_sec )
        {
            (*delta).tv_sec  = 0;
            (*delta).tv_nsec = 1000000000 + (*t_fin).tv_nsec - (*t_ini).tv_nsec;
        }
        else
        {
            (*delta).tv_sec  = (*t_fin).tv_sec - (*t_ini).tv_sec - 1;
            (*delta).tv_nsec = 1000000000 + (*t_fin).tv_nsec - (*t_ini).tv_nsec;
        }
    }
    else
    {
        if( (*t_fin).tv_sec == (*t_ini).tv_sec )
        {
            (*delta).tv_sec  = 0;
            (*delta).tv_nsec = (*t_fin).tv_nsec - (*t_ini).tv_nsec;
        }
        else
        {
            (*delta).tv_sec  = (*t_fin).tv_sec - (*t_ini).tv_sec;
            (*delta).tv_nsec = (*t_fin).tv_nsec - (*t_ini).tv_nsec;
        }
    }
}
