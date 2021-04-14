/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly.h"

int fmpz_mpoly_cmp(const fmpz_mpoly_t A, const fmpz_mpoly_t B,
                                                    const fmpz_mpoly_ctx_t ctx)
{
    int result;

    if (A->length != 1 || B->length != 1
        || !fmpz_is_one(A->coeffs + 0) || !fmpz_is_one(B->coeffs + 0))
    {
        flint_throw(FLINT_ERROR, "Inputs to cmp are not both monomials");
    }

    result = fmpz_cmp(A->new_exps + 0, B->new_exps + 0);

    return (result > 0) ? 1 : ((result < 0) ? -1 : 0);
}
