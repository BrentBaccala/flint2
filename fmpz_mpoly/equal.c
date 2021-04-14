/*
    Copyright (C) 2016 William Hart
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include <stdlib.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_mpoly.h"

int _fmpz_mpoly_equal_new(fmpz * poly1, fmpz * exps1,
                     const fmpz * poly2, const fmpz * exps2, slong n)
{
   slong i;

   if (poly1 != poly2)
   {
      for (i = 0; i < n; i++)
      {
         if (!fmpz_equal(poly1 + i, poly2 + i))
            return 0;
      }
   }

   if (exps1 != exps2)
   {
      for (i = 0; i < n; i++)
      {
         if (!fmpz_equal(exps1 + i, exps2 + i))
            return 0;
      }
   }

   return 1;
}

int fmpz_mpoly_equal(const fmpz_mpoly_t poly1, const fmpz_mpoly_t poly2,
                                                    const fmpz_mpoly_ctx_t ctx)
{
   int r;

   if (poly1 == poly2)
      return 1;

   if (poly1->length != poly2->length)
      return 0;

   r = _fmpz_mpoly_equal_new(poly1->coeffs, poly1->new_exps,
                             poly2->coeffs, poly2->new_exps, poly2->length);

   return r;
}
