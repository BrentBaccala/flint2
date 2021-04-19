/*
    Copyright (C) 2016 William Hart
    Copyright (C) 2017 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly.h"

void fmpz_mpoly_set_si(fmpz_mpoly_t A, slong c, const fmpz_mpoly_ctx_t ctx)
{
    if (c == 0)
    {
        _fmpz_mpoly_set_length(A, 0, ctx);
        return;
    }

    fmpz_mpoly_fit_length(A, 1, ctx);
    fmpz_set_si(A->coeffs + 0, c);
    fmpz_set_ui(A->new_exps + 0, 1);
    _fmpz_mpoly_set_length(A, 1, ctx);
}
