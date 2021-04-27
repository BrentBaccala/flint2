/*
    Copyright (C) 2016 William Hart

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

slong _fmpz_mpoly_add(fmpz * poly1, ulong * exps1,
                 const fmpz * poly2, const ulong * exps2, slong len2,
                 const fmpz * poly3, const ulong * exps3, slong len3, slong N,
                                                        const ulong * cmpmask)
{
  return 0;
}

slong _fmpz_mpoly_add_new(fmpz * poly1, fmpz * exps1,
                  const fmpz * poly2, const fmpz * exps2, slong len2,
                  const fmpz * poly3, const fmpz * exps3, slong len3)
{
   slong i = 0, j = 0, k = 0;

   while (i < len2 && j < len3)
   {
      int cmp = fmpz_cmp(exps2 + i, exps3 + j);

      if (cmp > 0)
      {
         fmpz_set(poly1 + k, poly2 + i);
         fmpz_set(exps1 + k, exps2 + i);
         i++;
      } else if (cmp == 0)
      {
         fmpz_add(poly1 + k, poly2 + i, poly3 + j);
         if (fmpz_is_zero(poly1 + k))
            k--;
         else
            fmpz_set(exps1 + k, exps2 + i);
         i++;
         j++;
      } else
      {
         fmpz_set(poly1 + k, poly3 + j);
         fmpz_set(exps1 + k, exps3 + j);
         j++;         
      }
      k++;
   }

   while (i < len2)
   {
      fmpz_set(poly1 + k, poly2 + i);
      fmpz_set(exps1 + k, exps2 + i);
      i++;
      k++;
   }

   while (j < len3)
   {
      fmpz_set(poly1 + k, poly3 + j);
      fmpz_set(exps1 + k, exps3 + j);
      j++;
      k++;
   }

   return k;
}

void fmpz_mpoly_add(fmpz_mpoly_t poly1, const fmpz_mpoly_t poly2,
                          const fmpz_mpoly_t poly3, const fmpz_mpoly_ctx_t ctx)
{
    slong len = 0;

   if (poly2->length == 0)
   {
      fmpz_mpoly_set(poly1, poly3, ctx);
      return;
   } else if (poly3->length == 0)
   {
      fmpz_mpoly_set(poly1, poly2, ctx);
      return;
   }

   if (poly1 == poly2 || poly1 == poly3)
   {
      fmpz_mpoly_t temp;

      fmpz_mpoly_init2(temp, poly2->length + poly3->length, ctx);

      len = _fmpz_mpoly_add_new(temp->coeffs, temp->new_exps, 
                    poly2->coeffs, poly2->new_exps, poly2->length,
                    poly3->coeffs, poly3->new_exps, poly3->length);

      fmpz_mpoly_swap(temp, poly1, ctx);

      fmpz_mpoly_clear(temp, ctx);
   } else
   {
      fmpz_mpoly_fit_length(poly1, poly2->length + poly3->length, ctx);

      len = _fmpz_mpoly_add_new(poly1->coeffs, poly1->new_exps, 
                       poly2->coeffs, poly2->new_exps, poly2->length,
                       poly3->coeffs, poly3->new_exps, poly3->length);
   }
      
   _fmpz_mpoly_set_length(poly1, len, ctx);
}
