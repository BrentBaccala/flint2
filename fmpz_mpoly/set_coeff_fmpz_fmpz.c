/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly.h"

void _fmpz_mpoly_set_coeff_fmpz_fmpz(fmpz_mpoly_t poly,
                  const fmpz_t c, const fmpz * exp, const fmpz_mpoly_ctx_t ctx)
{
    slong i, nvars = ctx->minfo->nvars;
    ulong * ulong_exp;
    TMP_INIT;

    TMP_START;

    ulong_exp = (ulong *) TMP_ALLOC(nvars*sizeof(ulong));

    for (i=0; i < nvars; i ++)
    {
        ulong_exp[i] = fmpz_get_ui(exp + i);
    }

    fmpz_mpoly_set_coeff_fmpz_ui(poly, c, ulong_exp, ctx);

    TMP_END;
}


void fmpz_mpoly_set_coeff_fmpz_fmpz(fmpz_mpoly_t poly,
                const fmpz_t c, fmpz * const * exp, const fmpz_mpoly_ctx_t ctx)
{
    slong i, nvars = ctx->minfo->nvars;
    ulong * ulong_exp;
    TMP_INIT;

    TMP_START;

    ulong_exp = (ulong *) TMP_ALLOC(nvars*sizeof(ulong));

    for (i=0; i < nvars; i ++)
    {
        ulong_exp[i] = fmpz_get_ui(exp[i]);
    }

    fmpz_mpoly_set_coeff_fmpz_ui(poly, c, ulong_exp, ctx);

    TMP_END;
}
