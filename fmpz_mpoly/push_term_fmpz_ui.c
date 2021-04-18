/*
    Copyright (C) 2018, 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly.h"

void _fmpz_mpoly_push_exp_ui(fmpz_mpoly_t A,
                                 const ulong * exp, const fmpz_mpoly_ctx_t ctx)
{
    slong i;
    slong old_length = A->length;
    fmpz_t new_exp;
    fmpz_t prime;
    fmpz_t primepow;

    fmpz_init(new_exp);
    fmpz_init(prime);
    fmpz_init(primepow);
    fmpz_one(prime);

    for (i = 0; i < ctx->minfo->nvars; i++) {
      fmpz_nextprime(prime, prime, 1);
      if (exp[i] > 0) {
        fmpz_pow_ui(primepow, prime, exp[i]);
        fmpz_mul(new_exp, new_exp, primepow);
      }
    }

    fmpz_mpoly_fit_length(A, old_length + 1, ctx);
    A->length = old_length + 1;
    fmpz_set(A->new_exps + A->length - 1, new_exp);

    fmpz_clear(new_exp);
    fmpz_clear(prime);
    fmpz_clear(primepow);
}


void fmpz_mpoly_push_term_fmpz_ui(fmpz_mpoly_t A,
                 const fmpz_t c, const ulong * exp, const fmpz_mpoly_ctx_t ctx)
{
    _fmpz_mpoly_push_exp_ui(A, exp, ctx);
    fmpz_set(A->coeffs + A->length - 1, c);
}

void fmpz_mpoly_push_term_ui_ui(fmpz_mpoly_t A,
                        ulong c, const ulong * exp, const fmpz_mpoly_ctx_t ctx)
{
    _fmpz_mpoly_push_exp_ui(A, exp, ctx);
    fmpz_set_ui(A->coeffs + A->length - 1, c);
}

void fmpz_mpoly_push_term_si_ui(fmpz_mpoly_t A,
                        slong c, const ulong * exp, const fmpz_mpoly_ctx_t ctx)
{
    _fmpz_mpoly_push_exp_ui(A, exp, ctx);
    fmpz_set_si(A->coeffs + A->length - 1, c);
}

