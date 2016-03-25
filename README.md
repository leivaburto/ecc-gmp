# GMP-ECC

## Introduction

A proof-of-concept library implementing scalar multiplication (or point multiplication) for elliptic curves.
The code is written in C using GMP http://gmplib.org . 

Algorithms for scalar multiplication with unknown point:

* Left-to-right
* Right-to-left
* NAF
* wNAF
* Sliding window

Algorithms for scalar multiplication with fixed point:

* Fixed-base windowing (BGMW)
* Fixed-base comb

Parameters of the elliptic curve belong to NIST P-256.

## Authors

* Israel Leiva
* Cristóbal Leiva
* Félix Pérez
