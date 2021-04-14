/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly.h"

int fmpz_mpoly_repack_bits(fmpz_mpoly_t A, const fmpz_mpoly_t B,
                                 flint_bitcnt_t Abits, const fmpz_mpoly_ctx_t ctx)
{
    fmpz_mpoly_set(A, B, ctx);
    return 1;
}


int fmpz_mpoly_repack_bits_inplace(
    fmpz_mpoly_t A,
    flint_bitcnt_t Abits,
    const fmpz_mpoly_ctx_t ctx)
{
    return 1;
}
