/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly.h"

ulong fmpz_mpoly_get_coeff_ui_ui(const fmpz_mpoly_t A,
                                 const ulong * exp, const fmpz_mpoly_ctx_t ctx)
{
    slong index;
    fmpz_t newexp;
    int exists;

    fmpz_init(newexp);
    fmpz_mpoly_set_monomial_ui(newexp, exp, ctx);
    exists = _fmpz_mpoly_monomial_exists(&index, A->new_exps, newexp, A->length);
    fmpz_zero(newexp);

    if (! exists)
    {
        return 0;
    }
    else
    {
        FLINT_ASSERT(index < A->length);
        return fmpz_get_ui(A->coeffs + index);
    }
}
