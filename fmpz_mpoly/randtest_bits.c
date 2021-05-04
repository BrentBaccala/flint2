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
    ulong * ulong_exp;
    fmpz_t exp;
    TMP_INIT;

    fmpz_mpoly_zero(A, ctx);
    fmpz_mpoly_fit_length(A, length, ctx);

    TMP_START;
    ulong_exp = (ulong *) TMP_ALLOC(ctx->minfo->nvars*sizeof(ulong));
    fmpz_init(exp);

    for (i = 0; i < length; i++)
    {
        fmpz_randtest(A->coeffs + i, state, coeff_bits);

        /* Generate a random bit pattern for our exponent(s), but its
         * factorization might contain a prime that's larger than any
         * of our variable primes, so unpack it and repack it to
         * ensure that all of our monomials correspond to unique
         * exponent numbers.
         */
        fmpz_randtest_unsigned(exp, state, exp_bits);
        if (fmpz_is_zero(exp))
            fmpz_one(exp);
	fmpz_mpoly_get_monomial_ui(ulong_exp, exp, ctx->minfo);
	fmpz_mpoly_set_monomial_ui(A->new_exps + i, ulong_exp, ctx);

    }

    _fmpz_mpoly_set_length(A, length, ctx);

    TMP_END;
    fmpz_clear(exp);

    fmpz_mpoly_sort_terms(A, ctx);
    fmpz_mpoly_combine_like_terms(A, ctx);
}
