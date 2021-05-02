/*
    Copyright (C) 2016 William Hart
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly.h"

void fmpz_mpoly_get_monomial_ui(ulong * ulong_exp, const fmpz_t expin, const mpoly_ctx_t mctx)
{
    slong i;
    fmpz_t prime;
    fmpz_t exp;

    fmpz_init(prime);
    fmpz_one(prime);

    fmpz_init(exp);
    fmpz_set(exp, expin);

    for (i = 0; i < mctx->nvars; i++) {
        fmpz_nextprime(prime, prime, 1);
        ulong_exp[i] = 0;
        while (fmpz_divisible(exp, prime)) {
            ulong_exp[i] ++;
            fmpz_divexact(exp, exp, prime);
        }
    }

    fmpz_clear(prime);
    fmpz_clear(exp);
}
