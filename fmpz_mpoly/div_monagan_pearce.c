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
#include "longlong.h"

#include "heap.h"

/*
   Set polyq to the quotient poly2 by poly3 discarding remainder (with notional
   remainder coeffs reduced modulo the leading coeff of poly3), and return
   the length of the quotient. This version of the function assumes the
   exponent vectors all fit in a single word. The exponent vectors are
   assumed to have fields with the given number of bits. Assumes input polys
   are nonzero. Implements "Polynomial division using dynamic arrays, heaps
   and packed exponents" by Michael Monagan and Roman Pearce [1], except that
   we use a heap with smallest exponent at head. Note that if a < b then
   (n - b) < (n - b) where n is the maximum value a and b can take. The word
   "maxn" is set to an exponent vector whose fields are all set to such a
   value n. This allows division from left to right with a heap with smallest
   exponent at the head. Quotient poly is written in reverse order.
   [1] http://www.cecm.sfu.ca/~rpearcea/sdmp/sdmp_paper.pdf 
*/
slong _fmpz_mpoly_div_monagan_pearce(fmpz ** polyq, fmpz ** expq,
                slong * allocq, const fmpz * poly2, const fmpz * exp2,
            slong len2, const fmpz * poly3, const fmpz * exp3, slong len3)
{
    slong i, j, q_len, s;
    slong next_loc, heap_len = 2;
    fmpz_mpoly_heap_s * heap;
    mpoly_heap_t * chain;
    slong * store, * store_base;
    mpoly_heap_t * x;
    fmpz * q_coeff = *polyq;
    fmpz * q_exp = *expq;
    slong * hind;
    fmpz_t exp;
    fmpz_t exp_r;
    fmpz_t r, acc_lg;
    fmpz_t temp1;
    ulong acc_sm[3];
    int small;
    slong bits2, bits3;
    ulong lc_norm = 0, lc_abs = 0, lc_sign = 0, lc_n = 0, lc_i = 0;
    TMP_INIT;

    TMP_START;

    fmpz_init(acc_lg);
    fmpz_init(r);
    fmpz_init(exp);
    fmpz_init(exp_r);
    fmpz_init(temp1);

    /* whether intermediate computations q - a*b will fit in three words */
    bits2 = _fmpz_vec_max_bits(poly2, len2);
    bits3 = _fmpz_vec_max_bits(poly3, len3);
    /* allow one bit for sign, one bit for subtraction */
    small = FLINT_ABS(bits2) <= (FLINT_ABS(bits3) + FLINT_BIT_COUNT(len3) +
            FLINT_BITS - 2) && FLINT_ABS(bits3) <= FLINT_BITS - 2;

    /* alloc array of heap nodes which can be chained together */
    next_loc = len3 + 4;   /* something bigger than heap can ever be */
    heap = (fmpz_mpoly_heap_s *) TMP_ALLOC((len3 + 1)*sizeof(fmpz_mpoly_heap_s));
    for (i = 0; i < len3 + 1; i++)
       fmpz_init(heap[i].exp);
    chain = (mpoly_heap_t *) TMP_ALLOC(len3*sizeof(mpoly_heap_t));
    store = store_base = (slong *) TMP_ALLOC(2*len3*sizeof(mpoly_heap_t *));

    /* space for flagged heap indicies */
    hind = (slong *) TMP_ALLOC(len3*sizeof(slong));
    for (i = 0; i < len3; i++)
        hind[i] = 1;

    q_len = WORD(0);

    /* see description of divisor heap division in paper */
    s = len3;

    /* insert (-1, 0, exp2[0]) into heap */
    x = chain + 0;
    x->i = -WORD(1);
    x->j = 0;
    x->next = NULL;
    FMPZ_POLY_HEAP_ASSIGN(heap[1], exp2 + 0, x);

    /* precompute leading cofficient info assuming "small" case */
    if (small)
    {
        lc_abs = FLINT_ABS(poly3[0]);
        lc_sign = FLINT_SIGN_EXT(poly3[0]);
        count_leading_zeros(lc_norm, lc_abs);
        lc_n = lc_abs << lc_norm;
        invert_limb(lc_i, lc_n);
    }

    while (heap_len > 1)
    {
        fmpz_set(exp, heap[1].exp);

        _fmpz_mpoly_fit_length_new(&q_coeff, &q_exp, allocq, q_len + 1);
        fmpz_cdiv_qr(q_exp + q_len, exp_r, exp, exp3 + 0);

        /* take nodes from heap with exponent matching exp */

        if (! fmpz_is_zero(exp_r))
        {
            /* optimation: coeff arithmetic not needed */

            if (__fmpz_gt(exp3 + 0, exp))
            {
                /* optimization: no more quotient terms possible */
                goto cleanup;
            }

            do
            {
                x = _fmpz_mpoly_heap_pop(heap, &heap_len);
                do
                {
                    *store++ = x->i;
                    *store++ = x->j;
                    if (x->i != -WORD(1))
                        hind[x->i] |= WORD(1);

                } while ((x = x->next) != NULL);
            } while (heap_len > 1 && __fmpz_eq(heap[1].exp, exp));
        }
        else if (small)
        {
            /* optimization: small coeff arithmetic, acc_sm used below */

            acc_sm[0] = acc_sm[1] = acc_sm[2] = 0;
            do
            {
                x = _fmpz_mpoly_heap_pop(heap, &heap_len);
                do
                {
                    *store++ = x->i;
                    *store++ = x->j;
                    if (x->i != -WORD(1))
                        hind[x->i] |= WORD(1);

                    if (x->i == -WORD(1))
                        _fmpz_mpoly_add_uiuiui_fmpz(acc_sm, poly2 + x->j);
                    else
                        _fmpz_mpoly_submul_uiuiui_fmpz(acc_sm, poly3[x->i], q_coeff[x->j]);
                } while ((x = x->next) != NULL);
            } while (heap_len > 1 && __fmpz_eq(heap[1].exp, exp));
        }
        else
        {
            /* general coeff arithmetic */

            fmpz_zero(acc_lg);  
            do
            {
                x = _fmpz_mpoly_heap_pop(heap, &heap_len);
                do
                {
                    *store++ = x->i;
                    *store++ = x->j;
                    if (x->i != -WORD(1))
                        hind[x->i] |= WORD(1);

                    if (x->i == -WORD(1))
                        fmpz_add(acc_lg, acc_lg, poly2 + x->j);
                    else
                        fmpz_submul(acc_lg, poly3 + x->i, q_coeff + x->j);
                } while ((x = x->next) != NULL);
            } while (heap_len > 1 && __fmpz_eq(heap[1].exp, exp));
        }

        /* process nodes taken from the heap */
        while (store > store_base)
        {
            j = *--store;
            i = *--store;

            if (i == -WORD(1))
            {
                /* take next dividend term */
                if (j + 1 < len2)
                {
                    x = chain + 0;
                    x->i = i;
                    x->j = j + 1;
                    x->next = NULL;
                    _fmpz_mpoly_heap_insert(heap, exp2 + x->j, x,
                                                 &next_loc, &heap_len);
                }
            } else
            {
                /* should we go right? */
                if (  (i + 1 < len3)
                   && (hind[i + 1] == 2*j + 1)
                   )
                {
                    x = chain + i + 1;
                    x->i = i + 1;
                    x->j = j;
                    x->next = NULL;
                    hind[x->i] = 2*(x->j + 1) + 0;
                    fmpz_mul(temp1, exp3 + x->i, q_exp + x->j);
                    _fmpz_mpoly_heap_insert(heap, temp1, x,
                                                 &next_loc, &heap_len);
                }
                /* should we go up? */
                if (j + 1 == q_len)
                {
                    s++;
                } else if (  ((hind[i] & 1) == 1)
                          && ((i == 1) || (hind[i - 1] >= 2*(j + 2) + 1))
                          )
                {
                    x = chain + i;
                    x->i = i;
                    x->j = j + 1;
                    x->next = NULL;
                    hind[x->i] = 2*(x->j + 1) + 0;
                    fmpz_mul(temp1, exp3 + x->i, q_exp + x->j);
                    _fmpz_mpoly_heap_insert(heap, temp1, x,
                                                 &next_loc, &heap_len);
                }
            }
        }

        /* try to divide accumulated term by leading term */
        if (! fmpz_is_zero(exp_r))
            continue;

        if (small)
        {
            ulong d0, d1, ds = acc_sm[2];

            /* d1:d0 = abs(acc_sm[1:0]) assuming ds is sign extension of acc_sm[1] */
            sub_ddmmss(d1, d0, acc_sm[1]^ds, acc_sm[0]^ds, ds, ds);
            
            if ((acc_sm[0] | acc_sm[1] | acc_sm[2]) == 0)
                continue;

            if (ds == FLINT_SIGN_EXT(acc_sm[1]) && d1 < lc_abs)
            {
                ulong qq, rr, nhi, nlo;
                FLINT_ASSERT(0 < lc_norm && lc_norm < FLINT_BITS);
                nhi = (d1 << lc_norm) | (d0 >> (FLINT_BITS - lc_norm));
                nlo = d0 << lc_norm;
                udiv_qrnnd_preinv(qq, rr, nhi, nlo, lc_n, lc_i);
                (void) rr; /* silence compiler warning */

                if (qq == 0)
                    continue;

                if (qq <= COEFF_MAX)
                {
                    _fmpz_demote(q_coeff + q_len);
                    q_coeff[q_len] = qq;
                    if (ds != lc_sign)
                        q_coeff[q_len] = -q_coeff[q_len];
                }
                else
                {
                    small = 0;
                    fmpz_set_ui(q_coeff + q_len, qq);
                    if (ds != lc_sign)
                        fmpz_neg(q_coeff + q_len, q_coeff + q_len);
                }
            }
            else
            {
                small = 0;
                fmpz_set_signed_uiuiui(acc_lg, acc_sm[2], acc_sm[1], acc_sm[0]);
                goto large_lt_divides;
            }
        }
        else
        {
            if (fmpz_is_zero(acc_lg))
                continue;

large_lt_divides:

            fmpz_fdiv_qr(q_coeff + q_len, r, acc_lg, poly3 + 0);
            if (fmpz_is_zero(q_coeff + q_len))
                continue;
        }

        /* put newly generated quotient term back into the heap if neccesary */
        if (s > 1)
        {
            i = 1;
            x = chain + i;
            x->i = i;
            x->j = q_len;
            x->next = NULL;
            hind[x->i] = 2*(x->j + 1) + 0;
            fmpz_mul(temp1, exp3 + x->i, q_exp + x->j);
            _fmpz_mpoly_heap_insert(heap, temp1, x,
                                                 &next_loc, &heap_len);
        }
        s = 1;
        q_len++;
    }


cleanup:

    fmpz_clear(acc_lg);
    fmpz_clear(r);
    fmpz_clear(exp);
    fmpz_clear(exp_r);
    fmpz_clear(temp1);

    /* If on our last pass through the loop, we had a zero coefficient, then
     * we set q_exp + q_len without incrementing q_len, so make sure it gets
     * zero'ed out, as we don't want lingering bigints.  On the other hand,
     * if we incremented q_len on the last pass, then it might have moved
     * pass the end of the array, so don't zero if q_len == *allocq.
     */
    if (q_len < *allocq)
        fmpz_zero(q_exp + q_len);

    (*polyq) = q_coeff;
    (*expq) = q_exp;

    TMP_END;

    /* return quotient poly length */
    return q_len;
}


void fmpz_mpoly_div_monagan_pearce(fmpz_mpoly_t q, const fmpz_mpoly_t poly2, 
                          const fmpz_mpoly_t poly3, const fmpz_mpoly_ctx_t ctx)
{
    slong lenq;
    fmpz_mpoly_t temp1;
    fmpz_mpoly_struct * tq;

   /* check divisor is nonzero */
   if (poly3->length == 0)
      flint_throw(FLINT_DIVZERO, "Divide by zero in fmpz_mpoly_div_monagan_pearce");

   /* dividend zero, write out quotient */
   if (poly2->length == 0)
   {
      fmpz_mpoly_zero(q, ctx);
      return;
   }

   /* check divisor leading monomial is at most that of the dividend */
   if (__fmpz_lt(poly2->new_exps + 0, poly3->new_exps + 0))
   {
      fmpz_mpoly_zero(q, ctx);
      return;
   }

   /* take care of aliasing */
   if (q == poly2 || q == poly3)
   {
      fmpz_mpoly_init2(temp1, poly2->length/poly3->length + 1, ctx);
      tq = temp1;
   } else
   {
      fmpz_mpoly_fit_length(q, poly2->length/poly3->length + 1, ctx);
      tq = q;
   }

   /* do division with remainder */
   lenq = _fmpz_mpoly_div_monagan_pearce(&tq->coeffs, &tq->new_exps,
                         &tq->alloc, poly2->coeffs, poly2->new_exps, poly2->length, 
                                     poly3->coeffs, poly3->new_exps, poly3->length);

   /* take care of aliasing */
   if (q == poly2 || q == poly3)
   {
      fmpz_mpoly_swap(temp1, q, ctx);
      fmpz_mpoly_clear(temp1, ctx);
   } 

   _fmpz_mpoly_set_length(q, lenq, ctx);
}
