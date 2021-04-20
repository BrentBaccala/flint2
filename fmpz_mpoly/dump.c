/*
    Copyright (C) 2021 Brent Baccala

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly.h"

void _fmpz_mpoly_dump (fmpz_mpoly_struct * A)
{
    ulong i;

    for (i=0; i<A->length; i++) {
        fmpz_fprint(stderr, A->new_exps + i);
	fputc(' ', stderr);
        fmpz_fprint(stderr, A->coeffs + i);
	fputc('\n', stderr);
    }
}
