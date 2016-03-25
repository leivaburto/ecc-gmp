#ifndef __UTILS_H__
#define __UTILS_H__
#include <stdio.h>
#include <time.h>
#include <sys/time.h>
#include "ec.h"

int ecc_read_jcb_from_file(ecc_jcb_t *, char[], int);
int ecc_read_afn_from_file(ecc_afn_t *, char[], int);
int ecc_read_k_from_file(mpz_ptr, char[], int);
void ecc_print_jcb(ecc_jcb_t *);
void ecc_print_afn(ecc_afn_t *);

#endif
