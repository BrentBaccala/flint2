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
    assuming that the conversion can be done, set poly1 and poly1_shift
    so that poly2 is poly1 * X^poly1_shift
    TODO: handle multiprecision exponents
*/
void fmpz_mpoly_to_fmpz_poly(fmpz_poly_t poly1, slong * poly1_shift,
               const fmpz_mpoly_t poly2, slong var, const fmpz_mpoly_ctx_t ctx)
{
    slong i, shift, off, bits, N;
    ulong k;
    slong _shift = 0, len = poly2->length;
    fmpz * coeff = poly2->coeffs;
    ulong * exp = poly2->exps;

    if (poly2->bits > FLINT_BITS)
        flint_throw(FLINT_EXPOF, "Bits too high in fmpz_mpoly_to_fmpz_poly");

    bits = poly2->bits;
    N = mpoly_words_per_exp(bits, ctx->minfo);
    mpoly_gen_offset_shift_sp(&off, &shift, var, bits, ctx->minfo);

    fmpz_poly_zero(poly1);
    if (len > 0)
    {
        ulong mask = (-UWORD(1)) >> (FLINT_BITS - bits);
        _shift = (exp[N*(len - 1)] >> shift) & mask;
        for (i = 0; i < len; i++)
        {
            k = (exp[N*i + off] >> shift) & mask;
            k -= _shift;
            FLINT_ASSERT(((slong)k) >= 0);
            fmpz_poly_set_coeff_fmpz(poly1, k, coeff + i);
        }
    }

    *poly1_shift = _shift;
}

/*
    set poly1 to poly2 * X^shift2 in the variable var
*/
void fmpz_mpoly_from_fmpz_poly(fmpz_mpoly_t poly1, const fmpz_poly_t poly2,
                           slong shift2, slong var, const fmpz_mpoly_ctx_t ctx)
{
    flint_bitcnt_t bits;
    slong N;
    slong k;
    slong p_len;
    fmpz * p_coeff;
    ulong * p_exp;
    slong p_alloc;
    ulong * one;
    TMP_INIT;

    TMP_START;

    bits = 1 + FLINT_BIT_COUNT(FLINT_MAX(WORD(1),
                                            shift2 + fmpz_poly_degree(poly2)));
    if (bits > FLINT_BITS)
        flint_throw(FLINT_EXPOF, "Exponent overflow in fmpz_mpoly_from_fmpz_poly");
    bits = mpoly_fix_bits(bits, ctx->minfo);
    
    N = mpoly_words_per_exp_sp(bits, ctx->minfo);
    one = (ulong*) TMP_ALLOC(N*sizeof(ulong));
    mpoly_gen_monomial_sp(one, var,bits, ctx->minfo);

    fmpz_mpoly_fit_bits(poly1, bits, ctx);
    poly1->bits = bits;

    p_coeff = poly1->coeffs;
    p_exp = poly1->exps;
    p_alloc = poly1->alloc;
    p_len = 0;
    for (k = fmpz_poly_degree(poly2); k >= 0; k--)
    {
        _fmpz_mpoly_fit_length(&p_coeff, &p_exp, &p_alloc, p_len + 1, N);
        mpoly_monomial_mul_ui(p_exp + N*p_len, one, N, k + shift2);
        fmpz_poly_get_coeff_fmpz(p_coeff + p_len, poly2, k);
        p_len += !fmpz_is_zero(p_coeff + p_len);
    }

    poly1->coeffs = p_coeff;
    poly1->exps = p_exp;
    poly1->alloc = p_alloc;
    _fmpz_mpoly_set_length(poly1, p_len, ctx);

    TMP_END;
}


/*
    set A(x_var^Bstride[var]) to B/xbar^Bshifts
    it is asserted that the conversion is correct
*/
void _fmpz_mpoly_to_fmpz_poly_deflate(
    fmpz_poly_t A,
    const fmpz_mpoly_t B,
    slong var,
    const ulong * Bshift,
    const ulong * Bstride,
    const fmpz_mpoly_ctx_t ctx)
{
    slong i;
    slong len = B->length;
    ulong var_shift, var_stride;

    FLINT_ASSERT(len > 0);

    fmpz_poly_zero(A);
    var_shift = Bshift[var];
    var_stride = Bstride[var];
    for (i = 0; i < len; i++)
    {
        ulong k = fmpz_mpoly_get_term_var_exp_ui(B, i, var, ctx);
        FLINT_ASSERT(k >= var_shift);
        k -= var_shift;
        if (k != 0)
        {
            k /= var_stride;
        }
        fmpz_poly_set_coeff_fmpz(A, k, B->coeffs + i);
    }
}

/*
    set A to B(x_var^Astride[var])*xbar^Ashift
    A must be packed into bits = Abits
*/
void _fmpz_mpoly_from_fmpz_poly_inflate(
    fmpz_mpoly_t A,
    const fmpz_poly_t B,
    slong var,
    const ulong * Ashift,
    const ulong * Astride,
    const fmpz_mpoly_ctx_t ctx)
{
    slong k;
    slong Alen;
    slong Bdeg = fmpz_poly_degree(B);
    fmpz_t exp;
    fmpz_t exp2;

    FLINT_ASSERT(!fmpz_poly_is_zero(B));

    fmpz_init(exp);
    fmpz_init(exp2);

    fmpz_mpoly_gen(A, var, ctx);
    _fmpz_mpoly_exp_ui(exp, Ashift, ctx);
    fmpz_pow_ui(exp2, A->new_exps + 0, Astride[var]);

    Alen = 0;
    for (k = 0; k <= Bdeg; k ++)
    {
        fmpz_mpoly_fit_length(A, Alen + 1, ctx);
        fmpz_poly_get_coeff_fmpz(A->coeffs + Alen, B, k);
        if (!fmpz_is_zero(A->coeffs + Alen))
        {
            fmpz_set(A->new_exps + Alen, exp);
            Alen++;
        }
        fmpz_mul(exp, exp, exp2);
    }

    _fmpz_mpoly_set_length(A, Alen, ctx);

    fmpz_clear(exp);
    fmpz_clear(exp2);
}
