/*
    Copyright (C) 2021 Brent Baccala

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly.h"

int _fmpz_mpoly_monomial_exists(slong * index, const fmpz * poly_exps,
                 const fmpz_t exp, slong len)
{
    /* copied from glibc-2.31 */
    slong __l, __u, __idx;
    int __comparison;

    __l = 0;
    __u = len;
    while (__l < __u) {
        __idx = (__l + __u) / 2;
        __comparison = fmpz_cmp (exp, poly_exps + __idx);
        if (__comparison < 0)
          __u = __idx;
        else if (__comparison > 0)
          __l = __idx + 1;
        else
          *index = __idx;
          return 1;
    }

    *index = __l;
    return 0;
}
