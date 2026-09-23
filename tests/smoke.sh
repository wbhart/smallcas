#!/bin/sh
set -eu

actual=$(printf "%s\n" \
    "R = PolynomialRing(ZZ)" \
    "x = R('x')" \
    "f = 3*x^2 + 2*x + 1" \
    "f" \
    "f(10)" \
    "evaluate(f, 10)" \
    "evaluate_horner(f, 10)" \
    "evaluate_dc(f, 10)" \
    "g = x + 1" \
    "f(g)" \
    "compose(f, g)" \
    "compose_horner(f, g)" \
    "compose_dc(f, g)" \
    "taylor_shift(f, 2)" \
    "taylor_shift_horner(f, 2)" \
    "taylor_shift_dc(f, 2)" \
    "derivative(f)" \
    "derivative(f, 2)" \
    "discriminant(f)" \
    "is_squarefree(f)" \
    "squarefree_part(6*(x - 1)^2*(x + 2))" \
    "degree(f)" \
    "leading_coefficient(f)" \
    "constant_coefficient(f)" \
    "coeff(f, 1)" \
    "reverse(f)" \
    "truncate(f, 2)" \
    "shift_left(f, 2)" \
    "shift_right(f, 1)" \
    "height(f)" \
    "max_abs_bits(f)" \
    "inflate(x^2 + x + 1, 2)" \
    "deflation(x^4 + x^2 + 1)" \
    "deflate(x^4 + x^2 + 1, 2)" \
    "(4*x^2 + 6*x + 2)/(2*x + 2)" \
    "quo(x^2 + 1, x + 1)" \
    "divrem(x^2 + 1, x + 1)" \
    "(x^2 + 1) % (x + 1)" \
    "pseudodiv(x^2 + 1, 2*x + 1)" \
    "inv_series(1 - x, 6)" \
    "b = x^2 + x + 1" \
    "p = preinverse(b, 4)" \
    "a = b*(x^3 + 2*x + 3) + 5*x + 7" \
    "quo_preinv(a, b, p)" \
    "content(6*x^2 + 12*x + 6)" \
    "primitive_part(6*x^2 + 12*x + 6)" \
    "(6*x^2 + 12*x + 6)/6" \
    "gcd((x - 1)*(x + 2), (x - 1)*(x + 3))" \
    "gcd_pseudo((x - 1)*(x + 2), (x - 1)*(x + 3))" \
    "gcd_hgcd((x - 1)*(x + 2), (x - 1)*(x + 3))" \
    "gcd_lr((x - 1)*(x + 2), (x - 1)*(x + 3))" \
    "resultant(x - 1, x - 2)" \
    "xgcd(x + 1, x + 2)" \
    "q, r = divrem(x^2 + 1, x + 1)" \
    "q" \
    "r" \
    "h, (u, v) = xgcd(x + 1, x + 2)" \
    "h" \
    "u" \
    "v" \
    "quit" | ./smallcas | sed 's/smallcas> //g')

expected='smallcas iter30 -- ZZ and univariate ZZ polynomials -- type quit to exit
3*x^2 + 2*x + 1
321
321
321
321
3*x^2 + 8*x + 6
3*x^2 + 8*x + 6
3*x^2 + 8*x + 6
3*x^2 + 8*x + 6
3*x^2 + 14*x + 17
3*x^2 + 14*x + 17
3*x^2 + 14*x + 17
6*x + 2
6
-8
1
x^2 + x - 2
2
3
1
2
x^2 + 2*x + 3
2*x + 1
3*x^4 + 2*x^3 + x^2
3*x + 2
3
2
x^4 + x^2 + 1
2
x^2 + x + 1
2*x + 1
x - 1
(x - 1, 2)
2
(2*x - 1, 5)
x^5 + x^4 + x^3 + x^2 + x + 1
x^3 + 2*x + 3
6
x^2 + 2*x + 1
x^2 + 2*x + 1
x - 1
x - 1
x - 1
x - 1
-1
(1, (-1, 1))
x - 1
2
1
-1
1' 

if [ "$actual" != "$expected" ]; then
    printf '%s\n' "smoke test failed" "--- actual ---" "$actual" \
        "--- expected ---" "$expected" >&2
    exit 1
fi
printf '%s\n' "smoke test passed"


bench=$(printf "%s\n" \
    "@time on" \
    "a = 0" \
    "for i = 1:5; a = a + i" \
    "@time off" \
    "a" \
    "i" \
    "quit" | ./smallcas)

bench_values=$(printf '%s\n' "$bench" | sed '/^smallcas iter30 /d; /^time: /d')
bench_times=$(printf '%s\n' "$bench" |
    awk '/^time: [0-9][0-9]*\.[0-9][0-9]* s$/ { n++ } END { print n + 0 }')
if [ "$bench_values" != "15
5" ] || [ "$bench_times" -ne 2 ]; then
    printf '%s\n' "REPL loop/timing smoke test failed" "$bench" >&2
    exit 1
fi
printf '%s\n' "REPL loop/timing smoke test passed"
