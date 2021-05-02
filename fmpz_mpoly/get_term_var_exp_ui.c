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

ulong fmpz_mpoly_get_term_var_exp_ui(const fmpz_mpoly_t A, slong i,
                                         slong var, const fmpz_mpoly_ctx_t ctx)
{
    ulong result = 0;
    fmpz_t prime;
    fmpz_t exp;

    if ((ulong) i >= (ulong) A->length)
    {
        flint_throw(FLINT_ERROR, "Index out of range in fmpz_mpoly_get_term_var_exp_ui");
    }

    fmpz_init(prime);
    fmpz_one(prime);

    fmpz_init(exp);
    fmpz_set(exp, A->new_exps + i);

    do {
        fmpz_nextprime(prime, prime, 1);
    } while (var--);

    while (fmpz_divisible(exp, prime)) {
        result ++;
        fmpz_divexact(exp, exp, prime);
    }

    fmpz_clear(prime);
    fmpz_clear(exp);

    return result;
}
