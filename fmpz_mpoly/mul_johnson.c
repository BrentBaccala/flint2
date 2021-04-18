/*
    Copyright (C) 2017 Daniel Schultz

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

#include "heap.h"

/*
   Set poly1 to poly2*poly3 using Johnson's heap method. The function
   realocates its output and returns the length of the product.
*/
slong _fmpz_mpoly_mul_johnson_new(fmpz ** poly1, fmpz ** exp1, slong * alloc,
                 const fmpz * poly2, const fmpz * exp2, slong len2,
                 const fmpz * poly3, const fmpz * exp3, slong len3)
{
   slong i, j, k;
   slong next_loc;
   slong Q_len = 0, heap_len = 2; /* heap zero index unused */
   fmpz_mpoly_heap_s * heap;
   mpoly_heap_t * chain;
   slong * Q;
   mpoly_heap_t * x;
   fmpz * p1 = *poly1;
   fmpz * e1 = *exp1;
   ulong cy;
   ulong c[3], p[2]; /* for accumulating coefficients */
   fmpz_t exp;
   fmpz * exp_list;
   slong exp_next;
   slong * hind;
   int first, small;
   TMP_INIT;

   TMP_START;

   /* whether input coeffs are small, thus output coeffs fit in three words */
   small = _fmpz_mpoly_fits_small(poly2, len2) &&
                                           _fmpz_mpoly_fits_small(poly3, len3);

   next_loc = len2 + 4;   /* something bigger than heap can ever be */
   heap = (fmpz_mpoly_heap_s *) TMP_ALLOC((len2 + 1)*sizeof(fmpz_mpoly_heap_s));
   for (i = 0; i < len2 + 1; i++)
      fmpz_init(heap[i].exp);
   /* alloc array of heap nodes which can be chained together */
   chain = (mpoly_heap_t *) TMP_ALLOC(len2*sizeof(mpoly_heap_t));
   /* space for temporary storage of pointers to heap nodes */
   Q = (slong *) TMP_ALLOC(2*len2*sizeof(slong));
   /* allocate space for exponent vectors */
   exp_list = (fmpz *) TMP_ALLOC(len2*sizeof(fmpz));
   for (i = 0; i < len2; i++)
      fmpz_init(exp_list + i);

   /* space for heap indices */
   hind = (slong *) TMP_ALLOC(len2*sizeof(slong));
   for (i = 0; i < len2; i++)
       hind[i] = 1;

   /* start with no heap nodes and no exponent vectors in use */
   exp_next = 0;

   fmpz_init(exp);

   /* put (0, 0, exp2[0] + exp3[0]) on heap */
   x = chain + 0;
   x->i = 0;
   x->j = 0;
   x->next = NULL;

   exp_next ++;

   heap[1].next = x;
   fmpz_init(heap[1].exp);

   fmpz_mul(heap[1].exp, exp2 + 0, exp3 + 0);

   hind[0] = 2*1 + 0;

   /* output poly index starts at -1, will be immediately updated to 0 */
   k = -WORD(1);

   /* while heap is nonempty */
   while (heap_len > 1)
   {
      /* get copy of exponent field of heap top */
      fmpz_set(exp, heap[1].exp);

      /* realloc output poly ready for next product term */
      k++;
      _fmpz_mpoly_fit_length_new(&p1, &e1, alloc, k + 1);

      /* whether we are on first coeff product for this output exponent */
      first = 1;

      /* set temporary coeff to zero */
      c[0] = c[1] = c[2] = 0;

      /* while heap nonempty and contains chain with current output exponent */
      while (heap_len > 1 && fmpz_equal(heap[1].exp, exp))
      {
         /* pop chain from heap and set exponent field to be reused */
         fmpz_set(exp_list + (--exp_next), heap[1].exp);

         x = _fmpz_mpoly_heap_pop(heap, &heap_len);

         /* take node out of heap and put into store */
         hind[x->i] |= WORD(1);
         Q[Q_len++] = x->i;
         Q[Q_len++] = x->j;

         /* if output coeffs will fit in three words */
         if (small)
         {
            /* compute product of input poly coeffs */
            if (first)
            {
               smul_ppmm(c[1], c[0], poly2[x->i], poly3[x->j]);
               c[2] = -(c[1] >> (FLINT_BITS - 1));

               /* set output monomial */
	       fmpz_set(e1 + k, exp);

               first = 0; 
            } else /* addmul product of input poly coeffs */
            {
               smul_ppmm(p[1], p[0], poly2[x->i], poly3[x->j]);
               add_sssaaaaaa(cy, c[1], c[0], 0, c[1], c[0], 0, p[1], p[0]);
               c[2] += (0 <= (slong) p[1]) ? cy : cy - 1;
            }
      
            /* for every node in this chain */
            while ((x = x->next) != NULL)
            {
               /* addmul product of input poly coeffs */
               smul_ppmm(p[1], p[0], poly2[x->i], poly3[x->j]);
               add_sssaaaaaa(cy, c[1], c[0], 0, c[1], c[0], 0, p[1], p[0]);
               c[2] += (0 <= (slong) p[1]) ? cy : cy - 1;

               /* take node out of heap and put into store */
               hind[x->i] |= WORD(1);
               Q[Q_len++] = x->i;
               Q[Q_len++] = x->j;
            }
         } else /* output coeffs require multiprecision */
         {
            if (first) /* compute product of input poly coeffs */
            {
               fmpz_mul(p1 + k, poly2 + x->i, poly3 + x->j);
               
               /* set output monomial */
	       fmpz_set(e1 + k, exp);

               first = 0; 
            } else
            {  /* addmul product of input poly coeffs */
               fmpz_addmul(p1 + k, poly2 + x->i, poly3 + x->j);
            }

            /* for each node in this chain */
            while ((x = x->next) != NULL)
            {
               /* addmul product of input poly coeffs */
               fmpz_addmul(p1 + k, poly2 + x->i, poly3 + x->j);

               /* take node out of heap and put into store */
               hind[x->i] |= WORD(1);
               Q[Q_len++] = x->i;
               Q[Q_len++] = x->j;
            }
         }
      }

      /* for each node temporarily stored */
      while (Q_len > 0)
      {
         /* take node from store */
         j = Q[--Q_len];
         i = Q[--Q_len];

         /* should we go right? */
         if (  (i + 1 < len2)
            && (hind[i + 1] == 2*j + 1)
            )
         {
            x = chain + i + 1;
            x->i = i + 1;
            x->j = j;
            x->next = NULL;

            hind[x->i] = 2*(x->j+1) + 0;

            fmpz_mul(exp_list + exp_next, exp2 + x->i, exp3 + x->j);

            if (!_fmpz_mpoly_heap_insert(heap, exp_list + exp_next++, x,
                                         &next_loc, &heap_len))
               exp_next--;
         }

         /* should we go up? */
         if (  (j + 1 < len3)
            && ((hind[i] & 1) == 1)
            && (  (i == 0)
               || (hind[i - 1] >  2*(j + 2) + 1)
               || (hind[i - 1] == 2*(j + 2) + 1) /* gcc should fuse */
               )
            )
         {
            x = chain + i;
            x->i = i;
            x->j = j + 1;
            x->next = NULL;

            hind[x->i] = 2*(x->j+1) + 0;

            fmpz_mul(exp_list + exp_next, exp2 + x->i, exp3 + x->j);

            if (!_fmpz_mpoly_heap_insert(heap, exp_list + exp_next++, x,
                                         &next_loc, &heap_len))
               exp_next--;
         }
      }

      /* set output poly coeff from temporary accumulation, if not multiprec */
      if (small)
         fmpz_set_signed_uiuiui(p1 + k, c[2], c[1], c[0]);

      if (fmpz_is_zero(p1 + k))
         k--;
   }

   k++;

   (*poly1) = p1;
   (*exp1) = e1;
   
   TMP_END;

   return k;
}

void fmpz_mpoly_mul_johnson(
    fmpz_mpoly_t A,
    const fmpz_mpoly_t B,
    const fmpz_mpoly_t C,
    const fmpz_mpoly_ctx_t ctx)
{
    slong Alen;

    if (B->length == 0 || C->length == 0)
    {
        fmpz_mpoly_zero(A, ctx);
        return;
    }

    /* deal with aliasing and do multiplication */
    if (A == B || A == C)
    {
        fmpz_mpoly_t T;
        fmpz_mpoly_init2(T, B->length + C->length - 1, ctx);

        /* algorithm more efficient if smaller poly first */
        if (B->length > C->length)
        {
            Alen = _fmpz_mpoly_mul_johnson_new(&T->coeffs, &T->new_exps, &T->alloc,
                                                  C->coeffs, C->new_exps, C->length,
                                                  B->coeffs, B->new_exps, B->length);
        }
        else
        {
            Alen = _fmpz_mpoly_mul_johnson_new(&T->coeffs, &T->new_exps, &T->alloc,
                                                  B->coeffs, B->new_exps, B->length,
                                                  C->coeffs, C->new_exps, C->length);
        }

        fmpz_mpoly_swap(T, A, ctx);
        fmpz_mpoly_clear(T, ctx);
    }
    else
    {
        fmpz_mpoly_fit_length(A, B->length + C->length - 1, ctx);

        /* algorithm more efficient if smaller poly first */
        if (B->length > C->length)
        {
            Alen = _fmpz_mpoly_mul_johnson_new(&A->coeffs, &A->new_exps, &A->alloc,
                                                  C->coeffs, C->new_exps, C->length,
                                                  B->coeffs, B->new_exps, B->length);
        }
        else
        {
            Alen = _fmpz_mpoly_mul_johnson_new(&A->coeffs, &A->new_exps, &A->alloc,
                                                  B->coeffs, B->new_exps, B->length,
                                                  C->coeffs, C->new_exps, C->length);
        }
    }

    _fmpz_mpoly_set_length(A, Alen, ctx);
}
