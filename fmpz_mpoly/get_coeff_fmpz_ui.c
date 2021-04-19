/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly.h"

void fmpz_mpoly_get_coeff_fmpz_ui(fmpz_t c, const fmpz_mpoly_t A,
                                 const ulong * exp, const fmpz_mpoly_ctx_t ctx)
{
    slong index;
    fmpz_t newexp;
    int exists;

    fmpz_init(newexp);
    _fmpz_mpoly_exp_ui(newexp, exp, ctx);

    exists = _fmpz_mpoly_monomial_exists(&index, A->new_exps, newexp, A->length);

    fmpz_zero(newexp);

    if (! exists)
    {
        fmpz_zero(c);
    }
    else
    {
        FLINT_ASSERT(index < A->length);
        fmpz_set(c, A->coeffs + index);
    }
}
