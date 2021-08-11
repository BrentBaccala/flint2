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

    flint_printf("addmul_multi_threaded....");
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
        fmpz_mpoly_addmul_multi_threaded(h, fptr, &f_length, 1, ctx);
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
        fmpz_mpoly_addmul_multi_threaded(h, fptr, &f_length, 1, ctx);
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
        fmpz_mpoly_addmul_multi_threaded(h, fptr, &f_length, 1, ctx);
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

            fmpz_mpoly_addmul_multi_threaded(h, fptr, &f_length, 1, ctx);
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

            fmpz_mpoly_addmul_multi_threaded(k2, fptr, f_lengths, 2, ctx);
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

    /* Check f*g + f*g = 2*f*g */
    for (i = 0; i < tmul * flint_test_multiplier(); i++)
    {
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t f, g, k1, k2;
        const fmpz_mpoly_struct * fptr [] = {f, g, f, g};
        slong f_lengths[] = {2,2};
        slong len1, len2;
        flint_bitcnt_t coeff_bits, exp_bits1, exp_bits2;

        fmpz_mpoly_ctx_init_rand(ctx, state, 20);

        fmpz_mpoly_init(f, ctx);
        fmpz_mpoly_init(g, ctx);
        fmpz_mpoly_init(k1, ctx);
        fmpz_mpoly_init(k2, ctx);

        len1 = n_randint(state, 100);
        len2 = n_randint(state, 100);

        coeff_bits = n_randint(state, 200);

        for (j = 0; j < 2; j++)
        {
            exp_bits1 = n_randint(state, 100) + 2;
            exp_bits2 = n_randint(state, 100) + 2;

            fmpz_mpoly_randtest_bits(f, state, len1, coeff_bits, exp_bits1, ctx);
            fmpz_mpoly_randtest_bits(g, state, len2, coeff_bits, exp_bits2, ctx);

            fmpz_mpoly_mul(k1, f, g, ctx);
            fmpz_mpoly_scalar_mul_si(k1, k1, 2, ctx);
            fmpz_mpoly_assert_canonical(k1, ctx);

            fmpz_mpoly_addmul_multi_threaded(k2, fptr, f_lengths, 2, ctx);
            fmpz_mpoly_assert_canonical(k2, ctx);

            result = fmpz_mpoly_equal(k1, k2, ctx);
            if (!result)
            {
                printf("FAIL\n");
                flint_printf("Check f*g + f*g = 2*f*g\ni = %wd, j = %wd\n", i ,j);
                flint_abort();
            }
        }

        fmpz_mpoly_clear(f, ctx);
        fmpz_mpoly_clear(g, ctx);
        fmpz_mpoly_clear(k1, ctx);
        fmpz_mpoly_clear(k2, ctx);
        fmpz_mpoly_ctx_clear(ctx);
    }

    /* Check f*g + h*k + l*m + s*t */
    for (i = 0; i < tmul * flint_test_multiplier(); i++)
    {
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t f, g, h, k, l, m, s, t, r1, r2;
        const fmpz_mpoly_struct * fptr [] = {f, g, h, k, l, m, s, t};
        slong f_lengths[] = {2,2,2,2};
        slong len1, len2;
        flint_bitcnt_t coeff_bits, exp_bits1, exp_bits2;

        fmpz_mpoly_ctx_init_rand(ctx, state, 20);

        fmpz_mpoly_init(f, ctx);
        fmpz_mpoly_init(g, ctx);
        fmpz_mpoly_init(h, ctx);
        fmpz_mpoly_init(k, ctx);
        fmpz_mpoly_init(l, ctx);
        fmpz_mpoly_init(m, ctx);
        fmpz_mpoly_init(s, ctx);
        fmpz_mpoly_init(t, ctx);
        fmpz_mpoly_init(r1, ctx);
        fmpz_mpoly_init(r2, ctx);

        len1 = n_randint(state, 100);
        len2 = n_randint(state, 100);

        coeff_bits = n_randint(state, 200);

        for (j = 0; j < 2; j++)
        {
            exp_bits1 = n_randint(state, 100) + 2;
            exp_bits2 = n_randint(state, 100) + 2;

            fmpz_mpoly_randtest_bits(f, state, len1, coeff_bits, exp_bits1, ctx);
            fmpz_mpoly_randtest_bits(g, state, len2, coeff_bits, exp_bits2, ctx);
            fmpz_mpoly_randtest_bits(h, state, len1, coeff_bits, exp_bits1, ctx);
            fmpz_mpoly_randtest_bits(k, state, len2, coeff_bits, exp_bits2, ctx);
            fmpz_mpoly_randtest_bits(l, state, len1, coeff_bits, exp_bits1, ctx);
            fmpz_mpoly_randtest_bits(m, state, len1, coeff_bits, exp_bits1, ctx);
            fmpz_mpoly_randtest_bits(s, state, len2, coeff_bits, exp_bits2, ctx);
            fmpz_mpoly_randtest_bits(t, state, len2, coeff_bits, exp_bits2, ctx);

            fmpz_mpoly_mul(r1, f, g, ctx);
            fmpz_mpoly_mul(r2, h, k, ctx);
            fmpz_mpoly_add(r1, r1, r2, ctx);
            fmpz_mpoly_mul(r2, l, m, ctx);
            fmpz_mpoly_add(r1, r1, r2, ctx);
            fmpz_mpoly_mul(r2, s, t, ctx);
            fmpz_mpoly_add(r1, r1, r2, ctx);
            fmpz_mpoly_assert_canonical(r1, ctx);

            flint_set_num_threads(n_randint(state, max_threads) + 1);

            fmpz_mpoly_addmul_multi_threaded(r2, fptr, f_lengths, 4, ctx);
            fmpz_mpoly_assert_canonical(r2, ctx);

            result = fmpz_mpoly_equal(r1, r2, ctx);
            if (!result)
            {
                printf("FAIL\n");
                flint_printf("Check f*g + h*k + l*m + s*t\ni = %wd, j = %wd\n", i ,j);
                flint_abort();
            }
        }

        fmpz_mpoly_clear(f, ctx);
        fmpz_mpoly_clear(g, ctx);
        fmpz_mpoly_clear(h, ctx);
        fmpz_mpoly_clear(k, ctx);
        fmpz_mpoly_clear(l, ctx);
        fmpz_mpoly_clear(m, ctx);
        fmpz_mpoly_clear(s, ctx);
        fmpz_mpoly_clear(t, ctx);
        fmpz_mpoly_clear(r1, ctx);
        fmpz_mpoly_clear(r2, ctx);
        fmpz_mpoly_ctx_clear(ctx);
    }

    /* Check f*g + h*k + l*m + s*t with len(f) = len(g) = 1, to ensure
     * test coverage of the case where one term is very short
     */
    for (i = 0; i < tmul * flint_test_multiplier(); i++)
    {
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t f, g, h, k, l, m, s, t, r1, r2;
        const fmpz_mpoly_struct * fptr [] = {f, g, h, k, l, m, s, t};
        slong f_lengths[] = {2,2,2,2};
        slong len1, len2;
        flint_bitcnt_t coeff_bits, exp_bits1, exp_bits2;

        fmpz_mpoly_ctx_init_rand(ctx, state, 20);

        fmpz_mpoly_init(f, ctx);
        fmpz_mpoly_init(g, ctx);
        fmpz_mpoly_init(h, ctx);
        fmpz_mpoly_init(k, ctx);
        fmpz_mpoly_init(l, ctx);
        fmpz_mpoly_init(m, ctx);
        fmpz_mpoly_init(s, ctx);
        fmpz_mpoly_init(t, ctx);
        fmpz_mpoly_init(r1, ctx);
        fmpz_mpoly_init(r2, ctx);

        len1 = n_randint(state, 100);
        len2 = n_randint(state, 100);

        coeff_bits = n_randint(state, 200);

        for (j = 0; j < 2; j++)
        {
            exp_bits1 = n_randint(state, 100) + 2;
            exp_bits2 = n_randint(state, 100) + 2;

            fmpz_mpoly_randtest_bits(f, state, 1, coeff_bits, exp_bits1, ctx);
            fmpz_mpoly_randtest_bits(g, state, 1, coeff_bits, exp_bits2, ctx);
            fmpz_mpoly_randtest_bits(h, state, len1, coeff_bits, exp_bits1, ctx);
            fmpz_mpoly_randtest_bits(k, state, len2, coeff_bits, exp_bits2, ctx);
            fmpz_mpoly_randtest_bits(l, state, len1, coeff_bits, exp_bits1, ctx);
            fmpz_mpoly_randtest_bits(m, state, len1, coeff_bits, exp_bits1, ctx);
            fmpz_mpoly_randtest_bits(s, state, len2, coeff_bits, exp_bits2, ctx);
            fmpz_mpoly_randtest_bits(t, state, len2, coeff_bits, exp_bits2, ctx);

            fmpz_mpoly_mul(r1, f, g, ctx);
            fmpz_mpoly_mul(r2, h, k, ctx);
            fmpz_mpoly_add(r1, r1, r2, ctx);
            fmpz_mpoly_mul(r2, l, m, ctx);
            fmpz_mpoly_add(r1, r1, r2, ctx);
            fmpz_mpoly_mul(r2, s, t, ctx);
            fmpz_mpoly_add(r1, r1, r2, ctx);
            fmpz_mpoly_assert_canonical(r1, ctx);

            flint_set_num_threads(n_randint(state, max_threads) + 1);

            fmpz_mpoly_addmul_multi_threaded(r2, fptr, f_lengths, 4, ctx);
            fmpz_mpoly_assert_canonical(r2, ctx);

            result = fmpz_mpoly_equal(r1, r2, ctx);
            if (!result)
            {
                printf("FAIL\n");
                flint_printf("Check f*g + h*k + l*m + s*t\ni = %wd, j = %wd\n", i ,j);
                flint_abort();
            }
        }

        fmpz_mpoly_clear(f, ctx);
        fmpz_mpoly_clear(g, ctx);
        fmpz_mpoly_clear(h, ctx);
        fmpz_mpoly_clear(k, ctx);
        fmpz_mpoly_clear(l, ctx);
        fmpz_mpoly_clear(m, ctx);
        fmpz_mpoly_clear(s, ctx);
        fmpz_mpoly_clear(t, ctx);
        fmpz_mpoly_clear(r1, ctx);
        fmpz_mpoly_clear(r2, ctx);
        fmpz_mpoly_ctx_clear(ctx);
    }

    /* Check copying a single polynomial to the output */
    for (i = 0; i < tmul * flint_test_multiplier(); i++)
    {
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t f, r;
        const fmpz_mpoly_struct * fptr [] = {f};
        slong f_lengths[] = {1};
        slong len;
        flint_bitcnt_t coeff_bits, exp_bits;

        fmpz_mpoly_ctx_init_rand(ctx, state, 20);

        fmpz_mpoly_init(f, ctx);
        fmpz_mpoly_init(r, ctx);

        len = n_randint(state, 100);

        coeff_bits = n_randint(state, 200);

        for (j = 0; j < 2; j++)
        {
            exp_bits = n_randint(state, 100) + 2;

            fmpz_mpoly_randtest_bits(f, state, len, coeff_bits, exp_bits, ctx);

            flint_set_num_threads(n_randint(state, max_threads) + 1);

            fmpz_mpoly_addmul_multi_threaded(r, fptr, f_lengths, 1, ctx);
            fmpz_mpoly_assert_canonical(r, ctx);

            result = fmpz_mpoly_equal(r, f, ctx);
            if (!result)
            {
                printf("FAIL\n");
                flint_printf("Check f\ni = %wd, j = %wd\n", i ,j);
                flint_abort();
            }
        }

        fmpz_mpoly_clear(f, ctx);
        fmpz_mpoly_clear(r, ctx);
        fmpz_mpoly_ctx_clear(ctx);
    }

    /* Check f + g */
    for (i = 0; i < tmul * flint_test_multiplier(); i++)
    {
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t f, g, r1, r2;
        const fmpz_mpoly_struct * fptr [] = {f, g};
        slong f_lengths[] = {1,1};
        slong len1, len2;
        flint_bitcnt_t coeff_bits, exp_bits1, exp_bits2;

        fmpz_mpoly_ctx_init_rand(ctx, state, 20);

        fmpz_mpoly_init(f, ctx);
        fmpz_mpoly_init(g, ctx);
        fmpz_mpoly_init(r1, ctx);
        fmpz_mpoly_init(r2, ctx);

        len1 = n_randint(state, 100);
        len2 = n_randint(state, 100);

        coeff_bits = n_randint(state, 200);

        for (j = 0; j < 2; j++)
        {
            exp_bits1 = n_randint(state, 100) + 2;
            exp_bits2 = n_randint(state, 100) + 2;

            fmpz_mpoly_randtest_bits(f, state, len1, coeff_bits, exp_bits1, ctx);
            fmpz_mpoly_randtest_bits(g, state, len2, coeff_bits, exp_bits2, ctx);

            fmpz_mpoly_add(r1, f, g, ctx);
            fmpz_mpoly_assert_canonical(r1, ctx);

            flint_set_num_threads(n_randint(state, max_threads) + 1);

            fmpz_mpoly_addmul_multi_threaded(r2, fptr, f_lengths, 2, ctx);
            fmpz_mpoly_assert_canonical(r2, ctx);

            result = fmpz_mpoly_equal(r1, r2, ctx);
            if (!result)
            {
                printf("FAIL\n");
                flint_printf("Check f + g\ni = %wd, j = %wd\n", i ,j);
                flint_abort();
            }
        }

        fmpz_mpoly_clear(f, ctx);
        fmpz_mpoly_clear(g, ctx);
        fmpz_mpoly_clear(r1, ctx);
        fmpz_mpoly_clear(r2, ctx);
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
            fmpz_mpoly_addmul_multi_threaded(f, fptr, &f_length, 1, ctx);
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
            fmpz_mpoly_addmul_multi_threaded(g, fptr, &f_length, 1, ctx);
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

