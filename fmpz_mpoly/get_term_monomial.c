/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly.h"

void fmpz_mpoly_get_term_monomial(fmpz_mpoly_t M, const fmpz_mpoly_t A,
                                           slong i, const fmpz_mpoly_ctx_t ctx)
{
    if ((ulong) i >= (ulong) A->length)
    {
        flint_throw(FLINT_ERROR, "Index out of range in fmpz_mpoly_get_term_monomial");
    }

    fmpz_mpoly_fit_length(M, WORD(1), ctx);

    fmpz_set(M->new_exps + 0, A->new_exps + i);
    fmpz_one(M->coeffs + 0);
    _fmpz_mpoly_set_length(M, 1, ctx);
}
