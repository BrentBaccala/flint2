/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly.h"

void fmpz_mpoly_set_term_exp_fmpz(fmpz_mpoly_t A, 
                       slong i, fmpz * const * exp, const fmpz_mpoly_ctx_t ctx)
{
    slong nvars = ctx->minfo->nvars;
    ulong j;
    ulong * newexp;
    TMP_INIT;

    if ((ulong) i >= (ulong) A->length)
    {
        flint_throw(FLINT_ERROR, "Index out of range in fmpz_mpoly_set_term_exp_fmpz");
    }

    TMP_START;
    newexp = (ulong *) TMP_ALLOC(nvars*sizeof(ulong));

    for (j = 0; j < nvars; j ++) {
        newexp[j] = fmpz_get_ui(exp[j]);
    }

    fmpz_mpoly_set_term_exp_ui(A, i, newexp, ctx);

    TMP_END;
}
