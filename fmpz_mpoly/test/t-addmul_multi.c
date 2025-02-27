/*
    Copyright (C) 2018 Daniel Schultz
    Copyright (C) 2021 Brent Baccala

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly.h"

int
main(void)
{
    int i, j, result, max_threads = 5;
    slong tmul = 1;
    FLINT_TEST_INIT(state);
#ifdef _WIN32
    tmul = 2;
#endif

    flint_printf("addmul_multi....");
    fflush(stdout);

    /* check fixed cases (same as mul test) */
    for (i = 0; i < 1 + tmul; i++)
    {
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t f, g, h, k;
        const fmpz_mpoly_struct * fptr[] = {f, g};
        slong f_length = 2;
        const char * vars[] = {"x", "y", "z", "t"};

        fmpz_mpoly_ctx_init(ctx, 4, ORD_DEGLEX);
        fmpz_mpoly_init(f, ctx);
        fmpz_mpoly_init(g, ctx);
        fmpz_mpoly_init(h, ctx);
        fmpz_mpoly_init(k, ctx);

        /* designed to trigger the dense case; reduced exponent to 10 to make it run faster */
        fmpz_mpoly_set_str_pretty(f, "((1-x)*(1+y)*(1+z))^10", vars, ctx);
        fmpz_mpoly_set_str_pretty(g, "((1+x)*(1-y)*(1-z))^10", vars, ctx);
        fmpz_mpoly_set_str_pretty(k, "((1-x^2)*(1-y^2)*(1-z^2))^10", vars, ctx);
        fmpz_mpoly_addmul_multi(h, fptr, &f_length, 1, ctx);
        if (!fmpz_mpoly_equal(h, k, ctx))
        {
            printf("FAIL\n");
            flint_printf("Check fixed case 1\n");
            flint_abort();
        }

        /* designed to trigger the array case */
        fmpz_mpoly_set_str_pretty(f, "(1+x+y+z+t)^20", vars, ctx);
        fmpz_mpoly_set_str_pretty(g, "(1-x-y-z-t)^20", vars, ctx);
        fmpz_mpoly_set_str_pretty(k, "((1+x+y+z+t)*(1-x-y-z-t))^20", vars, ctx);
        fmpz_mpoly_addmul_multi(h, fptr, &f_length, 1, ctx);
        if (!fmpz_mpoly_equal(h, k, ctx))
        {
            printf("FAIL\n");
            flint_printf("Check fixed case 2\n");
            flint_abort();
        }

        /* designed to trigger the heap case */
        fmpz_mpoly_set_str_pretty(f, "(1+x^10)^50", vars, ctx);
        fmpz_mpoly_set_str_pretty(g, "(1+y^10)^50", vars, ctx);
        fmpz_mpoly_set_str_pretty(k, "((1+x^10)*(1+y^10))^50", vars, ctx);
        fmpz_mpoly_addmul_multi(h, fptr, &f_length, 1, ctx);
        if (!fmpz_mpoly_equal(h, k, ctx))
        {
            printf("FAIL\n");
            flint_printf("Check fixed case 3\n");
            flint_abort();
        }

        fmpz_mpoly_clear(f, ctx);
        fmpz_mpoly_clear(g, ctx);
        fmpz_mpoly_clear(h, ctx);
        fmpz_mpoly_clear(k, ctx);
        fmpz_mpoly_ctx_clear(ctx);

        flint_set_num_threads(n_randint(state, max_threads) + 1);
    }

    /* Check triple multiplication f*g*h = (f*g)*h */
    for (i = 0; i < tmul * flint_test_multiplier(); i++)
    {
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_struct f[3];
        const fmpz_mpoly_struct * fptr [] = {f + 0, f + 1, f + 2};
        slong f_length = 3;
        fmpz_mpoly_t g, h;
        slong len, len1, len2;
        flint_bitcnt_t coeff_bits, exp_bits, exp_bits1, exp_bits2;

        fmpz_mpoly_ctx_init_rand(ctx, state, 20);

        fmpz_mpoly_init(f + 0, ctx);
        fmpz_mpoly_init(f + 1, ctx);
        fmpz_mpoly_init(f + 2, ctx);
        fmpz_mpoly_init(g, ctx);
        fmpz_mpoly_init(h, ctx);

        len = n_randint(state, 100);
        len1 = n_randint(state, 100);
        len2 = n_randint(state, 100);

        coeff_bits = n_randint(state, 200);

        for (j = 0; j < 2; j++)
        {
            exp_bits = n_randint(state, 100) + 2;
            exp_bits1 = n_randint(state, 100) + 2;
            exp_bits2 = n_randint(state, 100) + 2;

            fmpz_mpoly_randtest_bits(f + 0, state, len, coeff_bits, exp_bits, ctx);
            fmpz_mpoly_randtest_bits(f + 1, state, len, coeff_bits, exp_bits, ctx);
            fmpz_mpoly_randtest_bits(f + 2, state, len1, coeff_bits, exp_bits1, ctx);

            fmpz_mpoly_mul(g, f + 0, f + 1, ctx);
            fmpz_mpoly_assert_canonical(g, ctx);
            fmpz_mpoly_mul(g, g, f + 2, ctx);
            fmpz_mpoly_assert_canonical(g, ctx);

            fmpz_mpoly_addmul_multi(h, fptr, &f_length, 1, ctx);
            fmpz_mpoly_assert_canonical(h, ctx);

            result = fmpz_mpoly_equal(g, h, ctx);
            if (!result)
            {
                printf("FAIL\n");
                flint_printf("Check f*g*h = (f*g)*h\ni = %wd, j = %wd\n", i ,j);
                flint_abort();
            }
        }

        fmpz_mpoly_clear(f + 0, ctx);
        fmpz_mpoly_clear(f + 1, ctx);
        fmpz_mpoly_clear(f + 2, ctx);
        fmpz_mpoly_clear(g, ctx);
        fmpz_mpoly_clear(h, ctx);
        fmpz_mpoly_ctx_clear(ctx);
    }

    /* Check f*(g + h) = f*g + f*h */
    for (i = 0; i < tmul * flint_test_multiplier(); i++)
    {
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t f, g, h, k1, k2, t1, t2;
        const fmpz_mpoly_struct * fptr [] = {f, g, f, h};
        slong f_lengths[] = {2,2};
        slong len, len1, len2;
        flint_bitcnt_t coeff_bits, exp_bits, exp_bits1, exp_bits2;

        fmpz_mpoly_ctx_init_rand(ctx, state, 20);

        fmpz_mpoly_init(f, ctx);
        fmpz_mpoly_init(g, ctx);
        fmpz_mpoly_init(h, ctx);
        fmpz_mpoly_init(k1, ctx);
        fmpz_mpoly_init(k2, ctx);
        fmpz_mpoly_init(t1, ctx);
        fmpz_mpoly_init(t2, ctx);

        len = n_randint(state, 100);
        len1 = n_randint(state, 100);
        len2 = n_randint(state, 100);

        coeff_bits = n_randint(state, 200);

        for (j = 0; j < 2; j++)
        {
            exp_bits = n_randint(state, 100) + 2;
            exp_bits1 = n_randint(state, 100) + 2;
            exp_bits2 = n_randint(state, 100) + 2;

            fmpz_mpoly_randtest_bits(k1, state, len, coeff_bits, exp_bits, ctx);
            fmpz_mpoly_randtest_bits(k2, state, len, coeff_bits, exp_bits, ctx);
            fmpz_mpoly_randtest_bits(f, state, len1, coeff_bits, exp_bits1, ctx);
            fmpz_mpoly_randtest_bits(g, state, len2, coeff_bits, exp_bits2, ctx);
            fmpz_mpoly_randtest_bits(h, state, len2, coeff_bits, exp_bits2, ctx);

            fmpz_mpoly_add(t1, g, h, ctx);
            fmpz_mpoly_assert_canonical(t1, ctx);
            fmpz_mpoly_mul(k1, f, t1, ctx);
            fmpz_mpoly_assert_canonical(k1, ctx);

            fmpz_mpoly_addmul_multi(k2, fptr, f_lengths, 2, ctx);
            fmpz_mpoly_assert_canonical(k2, ctx);

            result = fmpz_mpoly_equal(k1, k2, ctx);
            if (!result)
            {
                printf("FAIL\n");
                flint_printf("Check f*(g + h) = f*g + f*h\ni = %wd, j = %wd\n", i ,j);
                flint_abort();
            }
        }

        fmpz_mpoly_clear(f, ctx);
        fmpz_mpoly_clear(g, ctx);
        fmpz_mpoly_clear(h, ctx);
        fmpz_mpoly_clear(k1, ctx);
        fmpz_mpoly_clear(k2, ctx);
        fmpz_mpoly_clear(t1, ctx);
        fmpz_mpoly_clear(t2, ctx);
        fmpz_mpoly_ctx_clear(ctx);
    }

    /* Check aliasing first argument */
    for (i = 0; i < tmul * flint_test_multiplier(); i++)
    {
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t f, g, h;
        const fmpz_mpoly_struct * fptr[] = {f, g};
        slong f_length = 2;
        slong len, len1, len2, nvars;
        flint_bitcnt_t coeff_bits, exp_bound, exp_bound1, exp_bound2;

        fmpz_mpoly_ctx_init_rand(ctx, state, 4);
        nvars = ctx->minfo->nvars;

        fmpz_mpoly_init(f, ctx);
        fmpz_mpoly_init(g, ctx);
        fmpz_mpoly_init(h, ctx);

        exp_bound = 3 + n_randint(state, 1 + 100/nvars/nvars);
        exp_bound1 = 3 + n_randint(state, 1 + 100/nvars/nvars);
        exp_bound2 = 3 + n_randint(state, 1 + 100/nvars/nvars);

        len = exp_bound + 1;
        len1 = exp_bound1 + 1;
        len2 = exp_bound2 + 1;
        for (j = n_randint(state, nvars) + 2; j >= 0; j--)
        {
            len *= exp_bound + 1;
            len1 *= exp_bound1 + 1;
            len2 *= exp_bound2 + 1;
            len = FLINT_MIN(len, WORD(1000));
            len1 = FLINT_MIN(len, WORD(1000));
            len2 = FLINT_MIN(len, WORD(1000));
        }

        coeff_bits = n_randint(state, 200);

        for (j = 0; j < 2; j++)
        {
            fmpz_mpoly_randtest_bound(f, state, len1, coeff_bits, exp_bound1, ctx);
            fmpz_mpoly_randtest_bound(g, state, len2, coeff_bits, exp_bound2, ctx);
            fmpz_mpoly_randtest_bound(h, state, len, coeff_bits, exp_bound, ctx);

            flint_set_num_threads(n_randint(state, max_threads) + 1);

            fmpz_mpoly_mul(h, f, g, ctx);
            fmpz_mpoly_assert_canonical(h, ctx);
            flint_set_num_threads(n_randint(state, max_threads) + 1);
            fmpz_mpoly_addmul_multi(f, fptr, &f_length, 1, ctx);
            fmpz_mpoly_assert_canonical(f, ctx);
            result = fmpz_mpoly_equal(h, f, ctx);
            if (!result)
            {
                printf("FAIL\n");
                flint_printf("Check aliasing first arg\ni = %wd, j = %wd\n", i ,j);
                flint_abort();
            }
        }

        fmpz_mpoly_clear(f, ctx);
        fmpz_mpoly_clear(g, ctx);
        fmpz_mpoly_clear(h, ctx);
        fmpz_mpoly_ctx_clear(ctx);
    }

    /* Check aliasing second argument */
    for (i = 0; i < tmul * flint_test_multiplier(); i++)
    {
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t f, g, h;
        const fmpz_mpoly_struct * fptr[] = {f, g};
        slong f_length = 2;
        slong len, len1, len2, nvars;
        flint_bitcnt_t coeff_bits, exp_bound, exp_bound1, exp_bound2;

        fmpz_mpoly_ctx_init_rand(ctx, state, 4);
        nvars = ctx->minfo->nvars;

        fmpz_mpoly_init(f, ctx);
        fmpz_mpoly_init(g, ctx);
        fmpz_mpoly_init(h, ctx);

        exp_bound = 3 + n_randint(state, 1 + 100/nvars/nvars);
        exp_bound1 = 3 + n_randint(state, 1 + 100/nvars/nvars);
        exp_bound2 = 3 + n_randint(state, 1 + 100/nvars/nvars);

        len = exp_bound + 1;
        len1 = exp_bound1 + 1;
        len2 = exp_bound2 + 1;
        for (j = n_randint(state, nvars) + 2; j >= 0; j--)
        {
            len *= exp_bound + 1;
            len1 *= exp_bound1 + 1;
            len2 *= exp_bound2 + 1;
            len = FLINT_MIN(len, WORD(1000));
            len1 = FLINT_MIN(len, WORD(1000));
            len2 = FLINT_MIN(len, WORD(1000));
        }

        coeff_bits = n_randint(state, 200);

        for (j = 0; j < 2; j++)
        {
            fmpz_mpoly_randtest_bound(f, state, len1, coeff_bits, exp_bound1, ctx);
            fmpz_mpoly_randtest_bound(g, state, len2, coeff_bits, exp_bound2, ctx);
            fmpz_mpoly_randtest_bound(h, state, len, coeff_bits, exp_bound, ctx);

            flint_set_num_threads(n_randint(state, max_threads) + 1);

            fmpz_mpoly_mul(h, f, g, ctx);
            fmpz_mpoly_assert_canonical(h, ctx);
            flint_set_num_threads(n_randint(state, max_threads) + 1);
            fmpz_mpoly_addmul_multi(g, fptr, &f_length, 1, ctx);
            fmpz_mpoly_assert_canonical(g, ctx);
            result = fmpz_mpoly_equal(h, g, ctx);
            if (!result)
            {
                printf("FAIL\n");
                flint_printf("Check aliasing second arg\ni = %wd, j = %wd\n", i ,j);
                flint_abort();
            }
        }

        fmpz_mpoly_clear(f, ctx);
        fmpz_mpoly_clear(g, ctx);
        fmpz_mpoly_clear(h, ctx);
        fmpz_mpoly_ctx_clear(ctx);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}

