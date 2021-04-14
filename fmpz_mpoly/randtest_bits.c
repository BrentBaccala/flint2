/*
    Copyright (C) 2017, 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly.h"

void fmpz_mpoly_randtest_bits(fmpz_mpoly_t A, flint_rand_t state,
              slong length, flint_bitcnt_t coeff_bits, flint_bitcnt_t exp_bits,
                                                    const fmpz_mpoly_ctx_t ctx)
{
    slong i;

    fmpz_mpoly_zero(A, ctx);
    fmpz_mpoly_fit_length(A, length, ctx);

    for (i = 0; i < length; i++)
    {
        fmpz_randtest(A->coeffs + i, state, coeff_bits);
        fmpz_randtest_unsigned(A->new_exps + i, state, exp_bits);
    }
    A->length = length;

    fmpz_mpoly_sort_terms(A, ctx);
    fmpz_mpoly_combine_like_terms(A, ctx);
}
