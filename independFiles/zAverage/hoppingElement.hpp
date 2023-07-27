// This header file computes the hopping elements for z-averaging
// Based on  `Vivaldo L. Campo, Jr. and Luiz N. Oliveira
// Phys. Rev. B 72, 104432(2015)
#pragma once
#include "mpreal.h" // This file is copied from https://github.com/advanpix/mpreal
#include <vector>
mpfr::mpreal        onsite(mpfr::mpreal lambda, const mpfr::mpreal &ZZ, int i);
std::vector<double> calcTn(double doubleZZ, double dlambda);
