/*
    Copyright (C) 2016 William Hart
    Copyright (C) 2017 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include <string.h>
#include "fmpz_mpoly.h"

#define ALLOC_PER_VAR ((FLINT_BITS+4)/3)

char *
_fmpz_mpoly_get_str_pretty(const fmpz * coeffs, const fmpz * exps, slong len,
                        const char ** x_in, const mpoly_ctx_t mctx)
{
    char * str, ** x = (char **) x_in, *xtmp;
    slong i, j, bound, off;
    ulong * exponents;
    int first;
    TMP_INIT;

    if (len == 0)
    {
        str = flint_malloc(2);
        str[0] = '0';
        str[1] = '\0';
        return str;
    }

    TMP_START;

    if (x == NULL)
    {
        xtmp = (char *) TMP_ALLOC(mctx->nvars * ALLOC_PER_VAR * sizeof(char));
        x = (char **) TMP_ALLOC(mctx->nvars*sizeof(char *));
        for (i = 0; i < mctx->nvars; i++)
        {
            x[i] = xtmp + i * ALLOC_PER_VAR;
            flint_sprintf(x[i], "x%wd", i + 1);
        }
    }

    bound = 1;
    for (i = 0; i < len; i++)
        bound += fmpz_sizeinbase(coeffs + i, 10) + 1;

    exponents = (ulong *) TMP_ALLOC(mctx->nvars*sizeof(ulong));

    for (i = 0; i < len; i++) {
        fmpz_mpoly_get_monomial_ui(exponents, exps + i, mctx);
        for (j = 0; j < mctx->nvars; j++) {
            while (exponents[j] > 0) {
                bound ++;
                exponents[j] /= 10;
            }
            bound += strlen(x[j]) + 3;
        }
    }

    str = flint_malloc(bound);
    off = 0;

    for (i = 0; i < len; i++)
    {
        if (fmpz_sgn(coeffs + i) > 0 && i != 0)
            str[off++] = '+';
        if (coeffs[i] == -WORD(1))
            str[off++] = '-';
        if (coeffs[i] != WORD(1) && coeffs[i] != -WORD(1))
        {
            if (!COEFF_IS_MPZ(coeffs[i]))
                off += flint_sprintf(str + off, "%wd", coeffs[i]);
            else
                off += gmp_sprintf(str + off, "%Zd", COEFF_TO_PTR(coeffs[i]));
        }

        fmpz_mpoly_get_monomial_ui(exponents, exps + i, mctx);

        first = 1;

        for (j = 0; j < mctx->nvars; j++)
        {
            if (exponents[j] > 1)
            {
                if (!first || (coeffs[i] != WORD(1) && coeffs[i] != -WORD(1)))
                    off += flint_sprintf(str + off, "*");
                off += flint_sprintf(str + off, "%s^", x[j]);

                if (!COEFF_IS_MPZ(exponents[j]))
                    off += flint_sprintf(str + off, "%wd", exponents[j]);
                else
                    off += gmp_sprintf(str + off, "%Zd", COEFF_TO_PTR(exponents[j]));

                first = 0;
            }
            else if (exponents[j] == 1)
            {
                if (!first || (coeffs[i] != WORD(1) && coeffs[i] != -WORD(1)))
                    off += flint_sprintf(str + off, "*");
                off += flint_sprintf(str + off, "%s", x[j]);
                first = 0;
            }
        }

        if (fmpz_equal_ui(exps + i, 1) && (coeffs[i] == WORD(1) || coeffs[i] == -WORD(1)))
            off += flint_sprintf(str + off, "1");
    }

    for (i = 0; i < mctx->nvars; i++)
        fmpz_clear(exponents + i);

    TMP_END;

    return str;
}

char *
fmpz_mpoly_get_str_pretty(const fmpz_mpoly_t poly, const char ** x, const fmpz_mpoly_ctx_t ctx)
{
   return _fmpz_mpoly_get_str_pretty(poly->coeffs, poly->new_exps,
                             poly->length, x, ctx->minfo);
}
