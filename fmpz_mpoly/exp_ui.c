/*
    Copyright (C) 2021 Brent Baccala

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly.h"

void _fmpz_mpoly_exp_ui(fmpz_t new_exp, const ulong * exp, const fmpz_mpoly_ctx_t ctx)
{
    slong i;
    fmpz_t prime;
    fmpz_t primepow;

    fmpz_init(prime);
    fmpz_init(primepow);
    fmpz_one(prime);
    fmpz_one(new_exp);

    for (i = 0; i < ctx->minfo->nvars; i++) {
      fmpz_nextprime(prime, prime, 1);
      if (exp[i] > 0) {
        fmpz_pow_ui(primepow, prime, exp[i]);
        fmpz_mul(new_exp, new_exp, primepow);
      }
    }

    fmpz_clear(prime);
    fmpz_clear(primepow);
}
