/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly.h"

/*
    sort terms in [left, right) by exponent
    assuming that bits in position > pos are already sorted
    and assuming exponent vectors fit into one word
    and assuming that all bit positions that need to be sorted are in totalmask
*/
void _fmpz_mpoly_radix_sort1(fmpz_mpoly_t A, slong left, slong right,
                               flint_bitcnt_t pos, ulong cmpmask, ulong totalmask)
{
    ulong mask = UWORD(1) << pos;
    ulong cmp = cmpmask & mask;
    slong mid, cur;

    FLINT_ASSERT(left <= right);
    FLINT_ASSERT(pos < FLINT_BITS);

    /* do nothing on lists of 0 or 1 elements */
    if (left + 1 >= right)
    {
        return;
    }

    /* return if there is no information to sort on this bit */
    if ((totalmask & mask) == WORD(0))
    {
        --pos;
        if ((slong)(pos) >= 0)
        {
            _fmpz_mpoly_radix_sort1(A, left,  right, pos, cmpmask, totalmask);
        }
        return;
    }

    /* find first 'zero' */
    mid = left;
    while (mid < right && ((A->exps + 1*mid)[0] & mask) != cmp)
    {
        mid++;
    }

    /* make sure [left,mid)  doesn't match cmpmask in position pos 'one'
                 [mid,right)    does match cmpmask in position pos 'zero' */
    cur = mid;
    while (++cur < right)
    {
        if (((A->exps + 1*cur)[0] & mask) != cmp)
        {
            fmpz_swap(A->coeffs + cur, A->coeffs + mid);
            mpoly_monomial_swap(A->exps + 1*cur, A->exps + 1*mid, 1);
            mid++;
        }
    }

    --pos;
    if ((slong)(pos) >= 0)
    {
        _fmpz_mpoly_radix_sort1(A, left,  mid, pos, cmpmask, totalmask);
        _fmpz_mpoly_radix_sort1(A, mid, right, pos, cmpmask, totalmask);
    }
}


/*
    sort terms in [left, right) by exponent
    assuming that bits in position > pos are already sorted

    TODO: Stack depth is proportional to N*FLINT_BITS
            Might turn into iterative version
            Low priority
*/
void _fmpz_mpoly_radix_sort(fmpz_mpoly_t A, slong left, slong right,
                                     flint_bitcnt_t pos, slong N, ulong * cmpmask)
{
    ulong off = pos/FLINT_BITS;
    ulong bit = pos%FLINT_BITS;
    ulong mask = UWORD(1) << bit;
    ulong cmp = cmpmask[off] & mask;
    slong mid, check;

    FLINT_ASSERT(left <= right);
    FLINT_ASSERT(pos < N*FLINT_BITS);

    /* do nothing on lists of 0 or 1 elements */
    if (left + 1 >= right)
        return;

    /* find first 'zero' */
    mid = left;
    while (mid < right && ((A->exps+N*mid)[off] & mask) != cmp)
    {
        mid++;
    }

    /* make sure [left,mid)  doesn't match cmpmask in position pos 'one'
                 [mid,right)    does match cmpmask in position pos 'zero' */
    check = mid;
    while (++check < right)
    {
        if (((A->exps + N*check)[off] & mask) != cmp)
        {
            fmpz_swap(A->coeffs + check, A->coeffs + mid);
            mpoly_monomial_swap(A->exps + N*check, A->exps + N*mid, N);
            mid++;
        }
    }

    --pos;
    if ((slong)(pos) >= 0)
    {
        _fmpz_mpoly_radix_sort(A, left,  mid, pos, N, cmpmask);
        _fmpz_mpoly_radix_sort(A, mid, right, pos, N, cmpmask);
    }
}

/* include some declarations for qsort.c (a modified GNU quicksort) */

typedef int (*__compar_d_fn_t) (const void *, const void *, void *);
typedef void (*__swap_d_fn_t) (void *, void *, void *);

void
_fmpz_mpoly_quicksort (void *const pbase, size_t total_elems, size_t size,
	    __compar_d_fn_t cmp, __swap_d_fn_t swap, void *arg);

int _fmpz_mpoly_compare (const void * a, const void * b, void * A)
{
  const fmpz * aa = (fmpz *) a;
  const fmpz * bb = (fmpz *) b;

  return -fmpz_cmp(a, b);
}

void _fmpz_mpoly_swap (void * a, void * b, void * A)
{
  fmpz * aa = (fmpz *) a;
  fmpz * bb = (fmpz *) b;
  fmpz_mpoly_struct * AA = (fmpz_mpoly_struct *) A;
  fmpz * aaa = (aa - AA->new_exps + AA->coeffs);
  fmpz * bbb = (bb - AA->new_exps + AA->coeffs);
  fmpz_swap(aa, bb);
  fmpz_swap(aaa, bbb);
}

/*
    sort the terms in A by exponent
    assuming that the exponents are valid (other than being in order)
*/
void fmpz_mpoly_sort_terms(fmpz_mpoly_t A, const fmpz_mpoly_ctx_t ctx)
{
    _fmpz_mpoly_quicksort(A->new_exps, A->length, sizeof(fmpz),
			  _fmpz_mpoly_compare, _fmpz_mpoly_swap, A);
}
