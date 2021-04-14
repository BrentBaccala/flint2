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

void fmpz_mpoly_set(fmpz_mpoly_t poly1, const fmpz_mpoly_t poly2,
                                                    const fmpz_mpoly_ctx_t ctx)
{

   fmpz_mpoly_fit_length(poly1, poly2->length, ctx);

   _fmpz_vec_set(poly1->new_exps, poly2->new_exps, poly2->length);
   _fmpz_vec_set(poly1->coeffs, poly2->coeffs, poly2->length);

   _fmpz_mpoly_set_length(poly1, poly2->length, ctx);
}
