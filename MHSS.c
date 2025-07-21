#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gmp.h>
#include <mpfr.h>
//#include <emscripten.h>

//EMSCRIPTEN_KEEPALIVE
int MHSS(unsigned long *input, unsigned long *output) {
    // input[0..9] = S[0..9]
    // input[10] = a, input[11] = p, input[12] = prec(bits)
    mpz_t p_big, a_big, ai_pow, a_plus_i;
    mpz_inits(p_big, a_big, ai_pow, a_plus_i, NULL);

    unsigned long a = input[10];
    unsigned long p = input[11];
    unsigned long prec = input[12]; //precision is input[12] bits
    unsigned long prec10 = (unsigned long)round((double)(prec * 0.30102999))+10UL;
    unsigned long row_l = p - a;

    mpz_set_ui(p_big, p);
    mpz_set_ui(a_big, a);

    // 精度設定
    mpfr_prec_t mp_prec = prec;     // precision
    mpfr_t *v1  = malloc(sizeof(mpfr_t) * row_l);    // value array
    mpfr_t *cfs = malloc(sizeof(mpfr_t) * row_l);   // coefficients array
    for (unsigned long i = 0; i < row_l; i++) {
        mpfr_init2(v1[i], mp_prec);
        mpfr_set_ui(v1[i], 1, MPFR_RNDN);
        mpfr_init2(cfs[i], mp_prec);
    }


    mpfr_t pss, inv_pow, tmp1, X;
    mpfr_inits2(mp_prec, pss, inv_pow, tmp1, X, (mpfr_ptr) NULL);
    mpfr_set_ui(pss, 1, MPFR_RNDN);

    // Sの長さを求める
    int s_length = 0;
    for (int i = 0; i < 10; i++) {
        printf("%lu\n",input[i]);
        if (input[i] == 0) {
            s_length++;
        } else {
            s_length = 0;
        }
    }
    s_length = 9 - s_length;
    

    for (int si = s_length; si >= 0; si--) {
        unsigned long s_end = input[si];

        mpz_pow_ui(ai_pow, p_big, s_end);
        mpfr_set_z(tmp1, ai_pow, MPFR_RNDN);
        mpfr_ui_div(inv_pow, 1, tmp1, MPFR_RNDN);
        mpfr_mul(pss, pss, inv_pow, MPFR_RNDN);

        // cfs[i] = 1 / (a + i)^s_end
        for (unsigned long i = 0; i < row_l; i++) {
            mpz_add_ui(a_plus_i, a_big, i);
            mpz_pow_ui(ai_pow, a_plus_i, s_end);
            mpfr_set_z(tmp1, ai_pow, MPFR_RNDN);
            mpfr_ui_div(cfs[i], 1, tmp1, MPFR_RNDN);
        }

        mpfr_set(X, pss, MPFR_RNDN);
        for (long i = row_l - 1; i >= 0; i--) {
            mpfr_mul(tmp1, cfs[i], v1[i], MPFR_RNDN);
            mpfr_add(X, X, tmp1, MPFR_RNDN);
            mpfr_set(v1[i], X, MPFR_RNDN);
        }
    }

    // 出力用
    mpfr_exp_t exp_ptr = 0;
    char *str_value = (char *)malloc((prec + 32) * sizeof(char));
    if (str_value == NULL) {
        fprintf(stderr, "malloc failed\n");
        return -987;
    }

    mpfr_get_str(str_value, &exp_ptr, 10, prec10, v1[0], MPFR_RNDN);
    size_t len = strlen(str_value);
    printf("%lu\n",len);
    for (size_t i = 0; i < len; ++i) {
        if (str_value[i] == '.') {
            output[i] = 10;
        } else {
            output[i] = (unsigned long)(str_value[i] - '0');
        }
    }

    // 後始末
    free(str_value);
    for (unsigned long i = 0; i < row_l; i++) {
        mpfr_clear(v1[i]);
        mpfr_clear(cfs[i]);
    }
    mpfr_clears(pss, inv_pow, tmp1, X,(mpfr_ptr) NULL);
    mpz_clears(p_big, a_big, ai_pow, a_plus_i, NULL);
    free(v1);
    free(cfs);

    return exp_ptr;
}
