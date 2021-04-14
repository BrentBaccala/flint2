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

void fmpz_mpoly_gen(fmpz_mpoly_t A, slong var, const fmpz_mpoly_ctx_t ctx)
{
    if (var < 0)
      flint_throw(FLINT_ERROR, "Negative variable in fmpz_mpoly_gen");

    fmpz_mpoly_fit_length(A, WORD(1), ctx);

    fmpz_one(A->coeffs);
    fmpz_one(A->new_exps);
    do {
      fmpz_nextprime(A->new_exps, A->new_exps, 1);
    } while (var --);

    _fmpz_mpoly_set_length(A, WORD(1), ctx);
}
