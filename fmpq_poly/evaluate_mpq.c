/*=============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

   Copyright (C) 2010 Sebastian Pancratz

******************************************************************************/

#include <mpir.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_poly.h"
#include "fmpq_poly.h"

void 
_fmpq_poly_evaluate_mpq(mpq_t res, const fmpz * poly, const fmpz_t den, 
                        ulong len, const mpq_t a)
{
    mpq_t rden;
    mpq_init(rden);
    _fmpz_poly_evaluate_horner_mpq(res, poly, len, a);
    mpz_set_ui(mpq_numref(rden), 1);
    fmpz_get_mpz(mpq_denref(rden), den);
    mpq_mul(res, res, rden);
    mpq_clear(rden);
}

void 
fmpq_poly_evaluate_mpq(mpq_t res, const fmpq_poly_t poly, const mpq_t a)
{
    if (res == a)
    {
        mpq_t t;
        mpq_init(t);
        _fmpq_poly_evaluate_mpq(t, poly->coeffs, poly->den, poly->length, a);
        mpq_swap(res, t);
        mpq_clear(t);
    }
    else
        _fmpq_poly_evaluate_mpq(res, poly->coeffs, poly->den, poly->length, a);
}

