/*
    Copyright (C) 2017 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly.h"

void fmpz_mpoly_set_coeff_fmpz_ui(fmpz_mpoly_t poly,
                 const fmpz_t c, const ulong * exp, const fmpz_mpoly_ctx_t ctx)
{
    slong i;
    slong index;
    fmpz_t newexp;
    int exists;

    fmpz_init(newexp);
    _fmpz_mpoly_exp_ui(newexp, exp, ctx);

    exists = _fmpz_mpoly_monomial_exists(&index, poly->new_exps, newexp, poly->length);

    if (!exists)
    {
        if (!fmpz_is_zero(c)) /* make new term only if coeff is nonzero*/
        {

            fmpz_mpoly_fit_length(poly, poly->length + 1, ctx);

            for (i = poly->length; i >= index + 1; i--)
            {
                fmpz_set(poly->coeffs + i, poly->coeffs + i - 1);
                fmpz_set(poly->new_exps + i, poly->new_exps + i - 1);
            }

            fmpz_set(poly->coeffs + index, c);
            fmpz_set(poly->new_exps + index, newexp);

            poly->length++; /* safe because length is increasing */
        }
    } else if (fmpz_is_zero(c)) /* zero coeff, remove term */
    {
        for (i = index; i < poly->length - 1; i++)
        {
            fmpz_set(poly->coeffs + i, poly->coeffs + i + 1);
            fmpz_set(poly->new_exps + i, poly->new_exps + i + 1);
        }

        _fmpz_mpoly_set_length(poly, poly->length - 1, ctx);

    } else /* term with that monomial exists, coeff is nonzero */
    {
        fmpz_set(poly->coeffs + index, c);
    }

    fmpz_zero(newexp);
}
