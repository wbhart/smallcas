# smallcas iter54

A deliberately tiny computer algebra system in C. The current checkpoint has `ZZ` and
univariate polynomial rings over `ZZ`, backed directly by GMP integers. GNU
Bison generates the parser, Flex generates the lexer, and GNU Readline provides
editable interactive input with persistent command history.

## Algorithms document

`algorithms.tex` is a mathematical description of the nontrivial algorithms,
with proofs of correctness and complexity; `algorithms.pdf` is the compiled
version.  The first three sections cover multiplication in `ZZ[x]`: ordinary
full multiplication; low and high truncated products; and middle products,
including transposition, Toom42 and Toom63.  The next four sections cover
division: classical and divide-and-conquer division; Mulders short division and
bidirectional exact division; reciprocal/Newton division; and classical and
fast pseudo-division.  These are followed by subresultant PRS normalization, primitive/subresultant
polynomial gcd algorithms, integral subresultant XGCD certificates,
fraction-free half-GCD, and resultants via Brown's subresultant PRS.

## WSL build

On Ubuntu under WSL:

```sh
sudo apt update
sudo apt install build-essential bison flex libgmp-dev libreadline-dev
make
make check
make test
./smallcas
```

A custom GMP installation can be selected with `CPPFLAGS` and `LDFLAGS`; the
Makefile does not pin a GMP release.

`make test-core` tests polynomial arithmetic, evaluation, division, gcd/resultant and
tuple unpacking without requiring Bison, Flex or Readline.

## Tuning

Tracked fallback thresholds live in `include/tuning_defaults.h`.  `include/tuning.h`
is a machine-local overlay, is ignored by Git, and is created automatically with no
overrides if it is missing.  This lets patches add new defaults without overwriting a
machine's measured values.

`make tune` builds a separate `SC_TUNE` core in which multiplication cutoffs are
writable variables.  Each timing point alternates the two algorithms, uses paired ratios,
aims for 1% MAD, and brackets the crossover between a 5% loss and a 5% win.  NTT cutoffs
are measured as requested/shorter length per CRT prime.  MFA is timed at successive
transform depths on the minimum-size Fermat ring allowed by SSA, using a
forward-plus-inverse round trip; if no 5% MFA win is found through depth 15 it is disabled.

The tuner covers the full-product chain, the mullo classical/DC and full-FFT
crossovers, independent mulhi full-FFT crossovers, and the balanced middle-product chain
classical -> Toom42 -> Toom63 -> FFT wraparound.  CRT-NTT and SSA have separate short
product cutoffs.  The Toom63 search begins after the tuned classical/Toom42 boundary so
the generated chain remains ordered.

Division tuning covers the recursive series-quotient and inverse bases, ordinary
quotient/divrem, Mulders short division, Newton/Karp--Markstein division, bidirectional
exact division, fast pseudo-division and the pseudo-remainder-only crossover used by PRS
algorithms.  These fast division algorithms inherit the tuned low and middle products,
hence NTT/SSA, rather than needing a separate FFT division kernel.  Division searches are
deliberately bounded in degree, number of sampled sizes and timing samples.  A crossover
requires a measured 5% loss followed by two consecutive 5% wins, so an isolated early
pocket is ignored.  A missing crossover is a valid result; optional paths such as public
Mulders, Newton or the fast pseudo-division/remainder paths remain disabled.

GCD tuning additionally chooses the small-size boundary between the primitive
pseudo-Euclidean gcd and Brown's subresultant PRS.  The search is capped at length 64; if
Brown does not win in that range, primitive gcd is still used only through that bounded
small-size range, so the public gcd retains the subresultant coefficient-growth guarantee
asymptotically.  Resultant tuning similarly chooses a bounded small Sylvester/Bareiss base
case before Brown's PRS, using the Sylvester order `deg(A)+deg(B)` as its metric.  The Brown
PRS itself needs no additional crossover: each pseudo-remainder already uses the tuned
pseudo-remainder dispatcher.  XGCD has no separate top-level crossover: its subresultant
certificate algorithm already inherits tuned pseudo-division and polynomial multiplication.

Only after every tuning stage succeeds does the tuner atomically replace
`include/tuning.h`.  Subsequent source patches therefore change `tuning_defaults.h`, not
the machine-local `tuning.h`.

The FFT layer provides cache-friendly TFT/ITFT transforms.  NTT and SSA full products use
a truncated round trip away from power-of-two endpoints.  Balanced middle products have
exact CRT-NTT and SSA cyclic-wraparound kernels.  Mullo and mulhi currently have no
substantially cheaper FFT algorithm; above their separately tuned crossovers they compute
a full exact FFT product and retain the requested edge.

## Example

```text
smallcas> R = PolynomialRing(ZZ)
smallcas> x = R('x')
smallcas> f = 3*x^2 + 2*x + 1
smallcas> f
3*x^2 + 2*x + 1
```

Calling a parent on an element performs coercion. Thus `R(7)` is the constant
polynomial 7. Mixed addition such as `x + 7` coerces 7 into `R`; `3*x` uses the
coefficient-ring scalar multiplication path.

## Integer evaluation

A polynomial bound to a variable is callable on an element of its coefficient ring:

```text
smallcas> f(10)
321
```

The normal evaluator dispatches between Horner evaluation and a divide-and-conquer
algorithm. The current untuned cutoff uses Horner for fewer than 16 coefficients (and
for arguments 0, 1 or -1), otherwise divide-and-conquer. The algorithms can also be
called explicitly:

```text
smallcas> evaluate(f, 10)
321
smallcas> evaluate_horner(f, 10)
321
smallcas> evaluate_dc(f, 10)
321
```

Horner has 12 counted lines. The divide-and-conquer implementation has 27 counted
lines and combines coefficient blocks with `a`, `a^2`, `a^4`, ... in a copied
coefficient workspace, avoiding temporary polynomial allocation.

## Interactive command line

Interactive input uses GNU Readline, so command lines can be edited and prior
commands recalled with the usual Readline keys. History is stored in
`~/.smallcas_history` and written after each nonblank command. Set
`SMALLCAS_HISTORY` to override the filename.

Piped or redirected input uses `getline()` instead, prints no prompt, and does
not read or update interactive history.

## Tuple unpacking

Pair-valued results can be unpacked directly in assignments. Both bare and
parenthesized pair patterns are accepted, and patterns may be nested:

```text
smallcas> q, r = divrem(x^2 + 1, x + 1)
smallcas> q
x - 1
smallcas> r
2
smallcas> h, (u, v) = xgcd(x + 1, x + 2)
smallcas> h
1
smallcas> u
-1
smallcas> v
1
```

A destructuring pattern must match the pair structure of the value. Ordinary
single-name assignment goes through the same pattern-assignment path.

## Dispatch and implementation layers

Dispatch is deliberately separate from mathematical implementations. Dispatch
functions are unrestricted by the 59-line rule: they may validate arguments,
inspect parents and representation metadata, choose algorithms, and use helper
functions which do those things. They may not hide part of a mathematical
algorithm.

For ordinary polynomial addition the structure is now:

```text
sc_add
  -> sc_poly_add
       -> sc_zz_poly_add
            -> sc_zz_poly_add_impl
```

`sc_zz_poly_add` is unrestricted dispatch/checking. `sc_zz_poly_add_impl` is
the actual mathematical implementation and is subject to the 59-line rule.
The same naming convention is used for subtraction, negation and scalar
multiplication. The `_impl` name is preferred to `_classical` when there is no
meaningful family of competing algorithms.

Multiplication really does have alternative algorithms:

```text
sc_mul
  -> sc_poly_mul
       -> sc_zz_poly_mul
            -> sc_zz_poly_mul_constant
            -> sc_zz_poly_mul_classical
            -> sc_zz_poly_mul_karatsuba
            -> sc_zz_poly_mul_ks
            -> sc_zz_poly_mul_toom3
```

The unrestricted `sc_zz_poly_mul` selector considers operand lengths and the
maximum absolute coefficient bit-size.  The raw bit-size scan is shared with the
public `max_abs_bits(f)` operation; it only measures coefficients and performs no
multiplication.

A balanced Kronecker-substitution implementation is also available. If the
maximum input coefficient bit-size is smaller than the shorter polynomial
length (with a small minimum-length guard), the current untuned heuristic uses
KS. The dispatcher chooses a packing width

```text
k = bits(a) + bits(b) + ceil(log2(min(length(a), length(b)))) + 1.
```

This makes every product coefficient strictly smaller than `2^(k-1)` in
absolute value. `sc_zz_poly_mul_ks` evaluates both polynomials at `2^k`, uses
one GMP integer multiplication, then recovers balanced base-`2^k` digits. The
KS implementation is 28 counted lines.

Karatsuba recursively calls the normal `ZZ[x]` multiplication dispatcher, so
subproducts may independently choose their implementation, including KS.

Balanced Toom-3 is also implemented. The dispatcher currently considers it
when the shorter operand has length at least 48, the operands differ in length
by at most a factor of 3/2, and the earlier KS heuristic has not fired. The
cutoff is intentionally untuned. With `m = ceil(max(length(a), length(b))/3)`,
the operands are written as

```text
a = a0 + a1 X + a2 X^2
b = b0 + b1 X + b2 X^2,       X = x^m,
```

and evaluated at `0, 1, -1, 2, infinity`. The five products all go back
through the ordinary `ZZ[x]` multiplication dispatcher.

The Toom-3 implementation is exactly 30 counted lines. It is not made compact
by replacing arithmetic with allocation-heavy generic polynomial expressions.
Evaluation is coefficientwise, using five GMP operations per coefficient for
each input; interpolation mutates the three middle product buffers in place;
and recombination is a direct shifted accumulation of the five blocks. The
interpolation uses exact divisions by 2 and 6.

Toom-3 uses `sc_zz_poly_toom3_ws`, whose initialization and cleanup live in
`memory.c` under the memory/representation-management exemption. That code
only creates the three non-owning block views, allocates the six evaluation
buffers and result, and manages ownership. It performs no polynomial
arithmetic. All evaluation formulas, all five recursive multiplications, the
full interpolation, and the shifted recombination remain visible in
`sc_zz_poly_mul_toom3`.


## Truncated products

Iteration 11 added genuine low, high and middle-window products. Iteration 12
added transposed-Karatsuba (Toom42) for balanced middle products, and iterations
13--14 added transposed Toom-3 (Toom63), including direct residue corrections.
`sc_zz_poly_mullow(a, b, n)` returns the first `n` coefficients of `a*b`. Its
classical implementation is 16 counted lines. A recursive truncated algorithm
splits at `k = ceil(n/2)`: it forms the complete low-block product and computes
only the required low parts of the two cross products. The high-high product is
never formed because it starts at degree at least `n`. The recursive algorithm
is 24 counted lines.

`sc_zz_poly_mulhigh(a, b, n)` returns the highest `n` product coefficients,
shifted down to degree zero. Its concrete implementation reverses only the
highest input tails that can contribute to the requested window, invokes the
low-product dispatcher, and
reverses the answer. High products therefore inherit the recursive low-product
implementation while the ancillary copying remains `O(n)` in the requested
window length.

The general operation

```text
sc_zz_poly_mulmid(a, b, start, n)
```

returns coefficients `start .. start+n-1` of `a*b`, shifted down. Its direct
classical implementation visits exactly the input coefficient pairs which
contribute to that window and is 21 counted lines. If the requested window
touches the low or high end, dispatch uses `mullow` or `mulhigh`.

For an interior window, dispatch also detects when one operand has length at
most the requested output length. After slicing the other operand, this is a
balanced middle product

```text
MP_n : (2*n - 1) by n -> n.
```

For `n > 16`, balanced middle products use transposed Karatsuba (Toom42). For
even `n = 2*k`, the long operand is split into the three overlapping blocks
`f0`, `f1`, `f2`, while the short operand is split into `g0`, `g1`. The three
recursive products are

```text
P0 = MP_k(f0 + f1, g1)
P1 = MP_k(f1, g0 - g1)
P2 = MP_k(f1 + f2, g0)
```

and recombination is

```text
(P0 + P1) + x^k*(P2 - P1).
```

The even Toom42 implementation is 27 counted lines. For odd `n`, the algorithm
recurses on the even `n-1` core after shifting the long operand by one
coefficient, then adds the omitted last row and final diagonal explicitly. The
odd-size algorithm is 23 counted lines, so no zero-padding or hidden arithmetic
helper is needed. The base case remains the direct classical middle product.

Large interior windows which do not have the balanced shape still fall back to
ordinary fast multiplication followed by extraction.

The D&C series quotient requests exactly coefficients `k .. n-1` of its
correction product. That window has the balanced shape after slicing, so for
large recursive steps it now goes through Toom42 instead of forming a complete
product. This gives the division code a genuine subquadratic middle-product
primitive rather than a full multiplication with discarded coefficients.

## Polynomial division

Iteration 10 added quotient-only and divide-and-conquer division to the existing
`ZZ[x]` division operations. The public quotient-only operation is `quo(a, b)`.
It returns the quotient even when the remainder is nonzero, provided every
leading-coefficient division needed by long division is exact in `ZZ`:

```text
smallcas> quo(x^2 + 1, x + 1)
x - 1
```

The classical quotient-only implementation does not compute and discard a
remainder. Its quotient array doubles as the live high part of the remainder;
coefficients below `deg(b)` are never updated because they cannot affect any
later quotient coefficient. It is 19 counted lines.

For larger quotients the current dispatcher uses divide-and-conquer division.
If `q = length(a) - length(b) + 1`, reversal converts the quotient problem into
truncated power-series division:

```text
rev(a) = rev(b) * rev(q)  (mod x^q).
```

The recursive series quotient computes the low half of `rev(q)`, obtains only
the required middle window of its product with the relevant prefix of `rev(b)`,
corrects the remaining coefficients, and recurses on the high half. The
recursive implementation is 29 counted lines; its classical series base case
is 18 lines. Reversal itself is a general 12-line arithmetic operation rather
than an exempt helper. Large balanced corrections now dispatch through Toom42
or Toom63 automatically.

Iteration 15 adds Mulders' short division with the simple balanced split
`beta = 1/2`. For a balanced problem with divisor and quotient length `n`, let
`n1 = ceil(n/2)` and `n2 = n - n1`. The high `n1`-term quotient block is found
by one ordinary full division. The only correction product needed for the low
block is the opposite short part of `b2*q1`; this is requested from the existing
`mulmid` dispatcher, so sufficiently large calls use Toom42 or Toom63 rather
than a full multiplication.

The implementation also reuses the remainder from that first full division.
If `u = a1*x^(n1-1) = b1*q1 + r1`, the low `n1-1` coefficients of `b1*q1`
are simply `-r1`. Thus it does not multiply `b1*q1` a second time merely to
construct the correction. After forming the corrected `2*n2-1`-term dividend,
Mulders division recurses for the low quotient block and concatenates the two
blocks. The concrete balanced recurrence is 39 counted lines.

The unrestricted Mulders wrapper also handles quotient length smaller than the
divisor length by taking only the leading `2*q-1` dividend coefficients and
leading `q` divisor coefficients. Quotients longer than the divisor currently
fall back to the existing D&C algorithm. This keeps the first Mulders
implementation close to the standard short-division problem rather than adding
an unrelated block loop.

Iteration 17 removes Mulders from the automatic public dispatcher. Iteration 18
adds Newton/Karp--Markstein division for sufficiently large unit-leading divisors.
The current untuned quotient/divrem selection is:

```text
q >= 64 and lc(b) = +/-1            Newton/Karp--Markstein
q >= 32                             middle-product D&C
otherwise                           classical quotient division
```

Thus monic divisors use Newton at the larger sizes requested in iteration 18; the
implementation also handles leading coefficient `-1`. `divrem(a, b)` follows the
same selection. Mulders remains implemented and tested through
`sc_zz_poly_quo_mulders` and `sc_zz_poly_divrem_mulders`, but it is now a
reference/experimental algorithm rather than a default dispatch choice. The
internal `sc_zz_poly_divrem_full` continues to select only classical or D&C
division for the retained Mulders implementation.

For small operands the original 20-line classical long-division `divrem`
remains available. The divisor need not be monic in the classical or D&C
algorithms: division proceeds exactly when its leading coefficient divides each
current leading coefficient. Newton division is restricted to unit-leading
divisors so that the reversed divisor has an integral power-series inverse.

### Newton division and reusable preinverses

`sc_zz_poly_inv_series` computes `f^(-1) mod x^n` when the constant coefficient
of `f` is `+/-1`. The Newton step does not form complete products. If an inverse
`g` is known to precision `m`, it asks `mulmid` for the high error block in
`f*g - 1` and uses `mullow` for the correction. The Newton implementation is 20
counted lines; the classical inverse base case is 15.

For repeated division by one divisor, `preinverse(b, n)` returns the first `n`
coefficients of the inverse of the reversed divisor. The same value can be reused:

```text
p = preinverse(b, 100)
q1 = quo_preinv(a1, b, p)
q2 = quo_preinv(a2, b, p)
```

The caller must request enough preinverse precision for the quotient being
computed. This is intentionally a lightweight toy API rather than a separate
preinverse object carrying precision metadata. `inv_series(f, n)` is also exposed
at the REPL.

For a one-off division, SmallCAS does not compute a full-precision reciprocal.
`sc_zz_poly_series_quo_newton` computes an inverse only to
`m = ceil(n/2)`, obtains the first quotient block with `mullow`, computes only the
necessary middle correction with `mulmid`, and obtains the remaining block by one
more low product with the same half-precision inverse. This is the
Karp--Markstein-style final correction. The routine is 29 counted lines.

`sc_zz_poly_divrem_newton` computes only the low `length(b)-1` coefficients of
`b*q` to recover the remainder, since the high coefficients have already been
cancelled by the quotient construction. It is 17 counted lines.

For large unit-leading divisors, exact `/` goes through the Newton `divrem` path
and simply requires the resulting remainder to vanish. For other divisors, exact
`/` retains the bidirectional path once the quotient has at least 32 coefficients.
After removing a common power of `x` from dividend and divisor, the low half of the
quotient is computed low-to-high by D&C series division, while the high half is
computed independently from reversed high coefficients. The two halves therefore
never correct one another.

Because SmallCAS does not assume the division is exact merely because `/` was used,
there is one final check. The low half already certifies the bottom quotient-length
half of `b*q`, and the reversed high half certifies the top half. Exactly
`length(b)-1` product coefficients remain unchecked, so `sc_zz_poly_divexact` asks
`mulmid` for precisely that middle window and compares it with the dividend. This
keeps the verification truncated rather than forming the full product. Divisors with
zero constant coefficient are handled by stripping their common `x`-adic valuation
first. Small cases continue through ordinary `divrem`.

`%` returns the remainder. Thus:

```text
smallcas> (4*x^2 + 6*x + 2)/(2*x + 2)
2*x + 1
smallcas> divrem(x^2 + 1, x + 1)
(x - 1, 2)
smallcas> (x^2 + 1) % (x + 1)
2
```

`pseudodiv(a, b)` remains defined for every nonzero divisor and returns `(q, r)`
with

```text
lc(b)^(deg(a)-deg(b)+1) * a = q*b + r
```

when `deg(a) >= deg(b)`, and exponent zero otherwise. It guarantees
`deg(r) < deg(b)`. Its implementation remains exactly 30 counted lines.

The division implementation counts are:

```text
sc_zz_poly_reverse_impl               12
sc_zz_poly_quo_classical              19
sc_zz_poly_series_quo_classical       18
sc_zz_poly_series_quo_dc              29
sc_zz_poly_inv_series_classical       15
sc_zz_poly_inv_series_newton          20
sc_zz_poly_series_quo_newton          29
sc_zz_poly_preinverse_newton           9
sc_zz_poly_quo_preinv                 15
sc_zz_poly_quo_newton                 17
sc_zz_poly_divrem_newton              17
sc_zz_poly_quo_dc                     12
sc_zz_poly_quo_bidirectional          26
sc_zz_poly_quo_mulders_balanced       39
sc_zz_poly_divrem_classical           20
sc_zz_poly_divrem_dc                   9
sc_zz_poly_divrem_mulders             10
sc_zz_poly_pseudodiv_impl             30
```


## Division families

The D&C quotient and Mulders quotient are retained as distinct divide-and-conquer
families. The former is reversal plus recursive power-series division; its correction
is a middle-product window. Mulders is recursive short division: it finds one quotient
block by a full division and uses an opposite short product before recursing on the
other block. Since SmallCAS now has genuine Toom42/Toom63 middle products, the D&C
family is the preferred automatic path and Mulders is kept for comparison and as a
compact reference implementation. The later Hart--Novocin half-full refinement is
therefore not currently planned.


## Polynomial views

`sc_zz_poly_view` creates a non-owning `sc_value` whose coefficient pointer
refers to a contiguous part of another polynomial. The logical length is
trimmed past trailing zero coefficients, so a view is a canonical ordinary
polynomial input as far as read-only arithmetic is concerned.

Karatsuba therefore creates views for `a0`, `a1`, `b0`, and `b1` and passes
those directly through `sc_zz_poly_add` and `sc_zz_poly_mul`. There is no
special `add_coeffwise` path merely because an operand is a view.

Views must not be freed, resized, truncated or otherwise used as owning output
objects. Direct coefficient writes in the Karatsuba recombination still target
the newly allocated result, not a view.

## Helper and macro policy

Two categories of helper are exempt from the mathematical 59-line rule:

1. Helpers used by dispatch functions for validation, metadata inspection and
   algorithm selection.
2. Helpers used for memory allocation, ownership, lifetime and representation
   management.

Neither exemption may be used to hide algebraic work.

Macros may be used in mathematical code, but a mathematical macro invocation
may replace at most one substantive C statement. It may not hide a loop,
branch, sequence of arithmetic operations, or several mathematical calls.
Representation-access macros are also permitted.

For example, `src/zz_poly.c` has:

```c
#define SC_MPZ_ADDEQ(x, y) mpz_add((x), (x), (y))
#define SC_MPZ_SUBEQ(x, y) mpz_sub((x), (x), (y))
```

Each invocation is still exactly one GMP arithmetic operation and consumes one
source line. There is deliberately no macro which hides the Karatsuba middle
term or its recombination loops.

## Balanced Toom63 middle product

Iteration 13 added the direct transposed Toom-3 (Toom63) algorithm for
`n = 3*k`. The long operand is viewed as five overlapping blocks and the short
operand as three blocks. The algorithm forms five explicit evaluation pairs,
performs five recursive balanced middle products, and interpolates directly
into the three output blocks.

The evaluation and interpolation arithmetic remains in
`sc_zz_poly_mulmid_toom63`; the workspace code only creates views, allocates
temporaries, and frees them. In particular, no general evaluation helper is
used to hide the Toom63 transforms. The interpolation needs one exact
division by 3 and two exact divisions by 2.

Iteration 14 extends Toom63 to every balanced size `n >= 48`. If
`n = m + r`, where `3 | m` and `r` is 1 or 2, it computes the `m`-term core
on the shifted long-operand view and the first `m` coefficients of the short
operand. It then adds the omitted `r` short-operand rows to the first `m`
outputs and computes the final `r` output diagonals directly. The correction
is linear work and does not change the recursive complexity.

`sc_zz_poly_mulmid_toom63` requires exactly **42 counted lines**, while the
generic residue-1/residue-2 tail routine requires 27. The direct Toom63 routine
is the reason the project-wide mathematical-function limit was raised from 30
to 42 lines. Existing algorithms have not been expanded merely because more
lines are now available.

The balanced middle-product dispatcher uses classical, Toom42 and Toom63 below
the transform range. Multiples of 3 use direct Toom63 and the other two residue classes
use the diagonal/row correction. These boundaries, together with the CRT-NTT and SSA
wraparound crossovers, are generated by `make tune`.

## GCD, content and resultant

Iteration 20 starts the `ZZ[x]` gcd/resultant layer. The public operations are

```text
content(f)
primitive_part(f)
gcd(f, g)
gcd_pseudo(f, g)
resultant(f, g)
```

Polynomial division by an integer scalar is also supported through `/`, and it
remains exact: every coefficient must be divisible by the scalar.

`content(f)` is the nonnegative gcd of the coefficients, with `content(0) = 0`.
`primitive_part(f)` divides by that content and preserves the sign of `f`. The
public gcd is normalized to have positive leading coefficient, so integer
content and primitive polynomial gcd are recombined according to Gauss's lemma.

`gcd_pseudo` is retained as a compact reference algorithm. It repeatedly takes
a pseudo-remainder and then removes its content. This primitive PRS keeps
coefficients much smaller than an unnormalised pseudo-Euclidean sequence, but
it still performs coefficient gcds at every stage.

The default `gcd` uses Brown's subresultant PRS. Its shared core computes the
last nonzero PRS polynomial together with the scalar `h` recurrence. If the
last polynomial has positive degree, its primitive part gives the primitive
gcd and the resultant is zero. If the last polynomial is constant, `h` is the
primitive resultant. Contents are then restored using

```text
Res(c*f, d*g) = c^deg(g) d^deg(f) Res(f, g).
```

Thus gcd and resultant genuinely share the same subresultant computation rather
than having unrelated implementations. The main new mathematical routine sizes
are:

```text
sc_zz_poly_content_impl                    10
sc_zz_poly_primitive_part_impl             11
sc_zz_poly_scalar_divexact_impl            19
sc_zz_poly_gcd_pseudo_impl                 31
sc_zz_poly_subres_prs_last                 39
sc_zz_poly_gcd_subresultant_impl           28
sc_zz_poly_resultant_subresultant_impl     41
sc_zz_poly_xgcd_subresultant_impl           59
```

Iteration 20 also adds `xgcd(f, g)`. It returns

```text
(h, (u, v))
```

with `u*f + v*g = h`, where `h` is the last nonzero Brown subresultant. The
Bézout cofactors are propagated through the same pseudo-remainder scaling and
exact divisions as the PRS itself, so they remain in `ZZ[x]`; no rational
polynomial layer is introduced.

For positive-degree coprime inputs, `h` is a nonzero integer which divides the
resultant, but it need not equal the resultant and need not be the smallest
positive integer in `(f, g)`. For example,

```text
f = -x^2 + x + 2
g =  x^2 - x + 1
xgcd(f, g) = (-3, (-1, -1))
resultant(f, g) = 9
```

so `(-1)*f + (-1)*g = -3`. Computing the minimal positive integer in the
intersection `(f, g) cap ZZ` is deliberately deferred until SmallCAS has
integer matrix/HNF/SNF machinery.

The direct cofactor-propagating subresultant routine is 59 counted lines. This
is the reason the project-wide mathematical-function limit is raised from 42
to 59: the additional lines are the two visible Bézout recurrences and their
exact PRS normalization, not hidden helper machinery. The subresultant and
xgcd workspaces in `memory.c` only own polynomial copies, basis cofactors and
GMP temporaries and rotate their ownership; they perform no polynomial
arithmetic.

### Experimental fraction-free half-GCD

Iteration 22 adds an explicit `gcd_hgcd(f, g)` experiment. It is deliberately
not used by the normal `gcd` dispatcher yet. The implementation is a direct
pseudo-Euclidean half-GCD over `ZZ[x]`: if pseudo-division gives

```text
lc(B)^delta A = Q B + R,
```

one pseudo-Euclidean step is represented by the integral transformation
matrix

```text
[ 0             1 ]
[ lc(B)^delta  -Q ].
```

The recursive half-GCD computes transformations from the high halves, applies
them to the full inputs, performs one pseudo-division step, and recurses again
on the next high halves. Matrix entries stay in `ZZ[x]`; neither `QQ[x]` nor
modular images are used.

To keep the pseudo-division itself subquadratic in degree, iteration 22 also
adds a scaled Newton series inverse. For a series `B` with constant term `c`,
it constructs `I_m` satisfying

```text
B I_m = c^m  (mod x^m)
```

and doubles precision with

```text
I_(2m) = I_m * (2*c^m - B*I_m)  (mod x^(2m)).
```

Reversal plus this scaled inverse gives a fast pseudo-quotient, while `mullow`
produces only the coefficients needed for the pseudo-remainder.

A raw pseudo-HGCD develops large common integer factors in its transformation
matrices. After matrix products, SmallCAS therefore divides all four entries
by their common integer content. This does not change the represented degree
reduction, and substantially reduces coefficient growth. It is still not the
full Lickteig--Roy coefficient-controlled fast subresultant algorithm, so
`gcd_hgcd` remains an experimental/reference path rather than the default gcd.

The main new routine sizes are:

```text
sc_zz_poly_inv_series_scaled_newton       30
sc_zz_poly_pseudodiv_fast                 45
sc_zz_poly_mat2_mul                       26
sc_zz_poly_mat2_apply                     17
sc_zz_poly_mat2_primitive                 31
sc_zz_poly_hgcd_pseudo                    53
sc_zz_poly_gcd_hgcd_impl                  47
```

All remain below the existing 59-line limit.

## Lickteig--Roy / quotient-boot half-GCD experiment

Iteration 23 adds an explicit `gcd_lr(f, g)` path. It is kept separate from
ordinary `gcd` for now, so Brown's subresultant PRS remains the default and a
correctness oracle.

The implementation follows the quotient-boot/half-GCD architecture used in
fast Sturm--Habicht/subresultant algorithms: recursively use only the high
parts of the current Euclidean remainders to predict a block of Euclidean
steps, compose the corresponding 2 by 2 transition matrix, apply it to the
full pair, perform the crossing Euclidean step, and recurse on the next high
parts. This is the degree-divide-and-conquer part underlying the
Lickteig--Roy approach.

A naive implementation over `QQ[x]`, or a matrix representation with one
common denominator, showed exactly the coefficient-growth problem we wanted
to avoid. Iteration 23 instead represents every temporary field polynomial as

```text
primitive ZZ[x] polynomial * one rational scalar.
```

The coefficient array therefore always remains an ordinary primitive
`ZZ[x]` polynomial. Rational arithmetic is confined to one GMP `mpq_t` scale
per temporary polynomial. After every addition or multiplication, any common
integer content is removed from the polynomial immediately and absorbed into
that scale. Matrix entries use the same representation independently, so an
unrelated denominator in one entry cannot inflate the other three.

Euclidean division of two such scaled polynomials does not construct a
rational-coefficient polynomial. It pseudo-divides their primitive integer
parts using the fast scaled-Newton pseudo-division from iteration 22, then
folds the predictable power of the divisor leading coefficient into the two
scalar factors. Thus the expensive polynomial arithmetic remains the existing
`ZZ[x]` multiplication/middle-product/pseudo-division machinery.

`gcd_lr` strips contents first, runs the quotient-boot half-GCD on the
primitive parts, and returns the primitive final Euclidean remainder with the
integer gcd of the original contents restored and positive leading
coefficient.

The LR source is deliberately exempt from the 59-line mathematical-function
limit for this iteration, as agreed for this experiment. The exemption is
local to `src/zz_poly_lr.c`; the 100-character line rule still applies. The
code is nevertheless split only along genuine representation operations
(normalize/copy/add/multiply/divide, 2 by 2 matrix operations, half-GCD, and
the gcd driver), rather than into algebra-hiding micro-helpers.

This iteration implements the fast quotient-boot/HGCD path needed for gcd. It
does **not** yet expose an operation that reconstructs an arbitrary selected
Sylvester--Habicht polynomial from the quotient boot. Brown's existing
subresultant implementation remains the source of canonical subresultants,
resultants and `xgcd` certificates.

`gcd_lr` is intentionally not selected by the normal `gcd` dispatcher yet.
The randomized regression suite compares it against Brown subresultant gcd on
inputs of lengths roughly 30--70. A separate development stress run also
checked 30 constructed-gcd cases with input lengths roughly 80--160; all
agreed with Brown and peaked at about 6 MB resident memory. This is the range
where the abandoned common-denominator matrix representation had already
shown pathological growth.

## 59-line rule

Every mathematical implementation is limited to 59 counted source lines.
`make check` applies the following counting rules:

- completely blank lines do not count;
- a line containing only optional whitespace and a single `}` does not count;
- all other lines in the mathematical function count, including its signature
  and opening brace.

This permits normal vertical spacing and conventional closing-brace placement
without reducing the visible amount of mathematics.

Dispatch files and `src/memory.c` are exempt for the reasons above. The
100-character maximum source-line length remains enforced for hand-written
sources, tests, build files and documentation.

## Testing

`make test-core` compares Karatsuba, balanced Kronecker substitution, balanced
Toom-3, Toom42 and Toom63 middle products against classical arithmetic over
signed inputs and a range of unequal lengths. It also tests quotient-only
classical and D&C division through the automatic dispatcher, while exercising Mulders
explicitly as a retained reference implementation. It tests the D&C `divrem`, explicit
Mulders `divrem`, exact division, ordinary `divrem`, and pseudo-division, including
randomized reconstruction and scaled-identity tests.
It exercises the KS, Toom-3, Karatsuba, and middle-product dispatch regions,
includes recursive Toom-3, Toom42, and Toom63 cases, explicitly tests both
nonzero residue classes modulo 3 for Toom63, checks odd and even middle-product
sizes, and checks a cancellation-heavy signed case. Because Karatsuba
and Toom-3 use non-owning block views, these
tests also exercise the view-as-input convention.

The gcd/resultant tests compare primitive pseudo-PRS gcd against subresultant
gcd on constructed examples, test contents and scalar exact division, compare
the experimental pseudo-HGCD and the LR quotient-boot gcd against the
subresultant gcd, compare fast Newton pseudo-division against the classical
pseudo-division implementation, and
compare randomized subresultant resultants against an independent Sylvester
matrix determinant computed by fraction-free Bareiss elimination. They also
verify randomized `xgcd` Bézout identities, including degree swaps, and check
that a constant terminal subresultant divides the resultant for positive-degree
inputs. The fixed defective-degree example above is included explicitly. Zero
and constant edge cases are also included.

`make test` additionally runs the parser/REPL smoke test.

## Iteration 24: quarantine of the LR experiment

The temporary rational-polynomial machinery used by `gcd_lr` is now fully private to
`src/zz_poly_lr.c`.  Its `qpoly` and `qmat2` structures, GMP `mpq_t` scales, cleanup
routines, quotient-boot half-GCD, and rational-matrix arithmetic are not declared in the
public SmallCAS header.  The rest of SmallCAS sees only the ordinary `ZZ[x]` entry point
`sc_zz_poly_gcd_lr_impl` through the existing `gcd_lr` dispatch path.

This is intentional quarantine: the LR path remains an experiment and must not become an
accidental dependency of the core object model before genuine `QQ`, `QQ[x]`, and matrix
types exist.

## Iteration 25: polynomial evaluation

`ZZ[x]` now supports evaluation at integers by Horner and divide-and-conquer
evaluation. Calling a bound polynomial, for example `f(123)`, performs evaluation.
The ordinary `evaluate(f, a)` call uses unrestricted dispatch; `evaluate_horner` and
`evaluate_dc` expose the two implementations directly for testing and comparison.
The randomized core tests compare both algorithms through degree 95 over signed
coefficients and large signed integer arguments, and also test the callable-polynomial
dispatch path.


## Iteration 26: polynomial composition

`ZZ[x]` now supports matrix-free polynomial composition.  If `f` and `g` have the same
polynomial parent, `f(g)` means `f` composed with `g`; calling `f(a)` for `a` in `ZZ`
still means scalar evaluation.  `compose(f, g)` uses unrestricted dispatch, while
`compose_horner(f, g)` and `compose_dc(f, g)` expose the two implementations.

Horner composition repeatedly multiplies the accumulated polynomial by `g` and adds
the next coefficient.  Divide-and-conquer composition starts with the coefficients as
constant polynomials and combines adjacent blocks using `g`, `g^2`, `g^4`, and so on.
All polynomial products go through the normal `ZZ[x]` multiplication dispatcher.  No
matrix or Brent--Kung machinery is used.

The randomized composition tests compare Horner, divide-and-conquer, and automatic
dispatch, and independently verify `(f compose g)(a) = f(g(a))` using scalar
evaluation.

## Iteration 27: calculus and basic invariants

`ZZ[x]` now provides formal differentiation, higher derivatives, discriminant,
squarefreeness testing, and squarefree part:

```text
derivative(f)
derivative(f, n)
nth_derivative(f, n)
discriminant(f)
is_squarefree(f)
squarefree_part(f)
```

The first derivative is a direct coefficientwise implementation.  The nth derivative
uses `k! * binomial(i, k)` for the falling-factorial coefficient multiplier, avoiding
repeated construction of intermediate derivative polynomials.

The discriminant uses the existing subresultant resultant together with
`disc(f) = (-1)^(n(n-1)/2) resultant(f, f') / lc(f)`.  As in Sage, the
discriminant of a constant or zero polynomial is defined to be zero; a nonconstant
linear polynomial has discriminant one.

Squarefreeness is understood over the fraction field: scalar content is ignored, so
nonzero linear polynomials are squarefree regardless of content.  The zero polynomial
is also defined to be squarefree.  `squarefree_part(f)` computes `f / gcd(f, f')`,
thereby removing both repeated nonconstant factors and scalar content, and normalises
the result to positive leading coefficient.  The zero polynomial maps to zero.

The substantive line counts remain small: first derivative, nth derivative,
discriminant, squarefree test, and squarefree part all fit comfortably below the
59-line project limit.  `tests/calculus.c` checks fixed derivatives, higher-order edge
cases, quadratic/cubic/linear/constant discriminants, repeated factors, scalar content,
and named dispatch through `derivative(f, n)`.


## Iteration 28: polynomial structure utilities

The remaining small `ZZ[x]` structural operations are now exposed without adding new
infrastructure:

```text
degree(f)
leading_coefficient(f)
constant_coefficient(f)
coeff(f, n)
reverse(f)
reverse(f, n)
truncate(f, n)
shift_left(f, n)
shift_right(f, n)
height(f)
max_abs_bits(f)
inflate(f, n)
deflation(f)
maximal_deflation(f)
deflate(f, n)
```

`degree(0)` is `-1`.  Coefficients beyond the current degree, and the leading or
constant coefficient of the zero polynomial, are returned as zero.  `reverse(f, n)`
regards `f` as having notional length `n`, truncating or zero-padding before reversal;
`reverse(f)` uses its current length.  `truncate` is non-mutating at the language level.
The two shifts move coefficients by powers of `x`, with right shift discarding the low
coefficients.

`height(f)` is the largest absolute coefficient and `max_abs_bits(f)` is its maximum
absolute coefficient bit-size, with both zero on the zero polynomial.  Inflation sends
`f(x)` to `f(x^n)`; the `n = 0` case therefore gives the constant `f(1)`.  `deflation(f)`
(or its alias `maximal_deflation(f)`) returns the largest exponent factor by which `f`
can be deflated, with the conventions `deflation(0) = 0` and deflation one for a
nonzero constant.  `deflate(f, n)` requires `n > 0` and rejects a polynomial having a
nonzero term whose exponent is not divisible by `n`.

The implementation reuses the existing 12-line reversal routine.  Every new
mathematical routine remains comfortably below the 59-line limit, and
`tests/structural.c` checks fixed accessors, reversal/truncation/shifts, inflation and
deflation edge cases, and named REPL dispatch.  The raw coefficient-bit scan is shared
by `max_abs_bits` and the multiplication dispatcher rather than duplicated.

## Iteration 29: Taylor shift

`ZZ[x]` now supports the integer Taylor shift `f(x) -> f(x + c)`:

```text
taylor_shift(f, c)
taylor_shift_horner(f, c)
taylor_shift_dc(f, c)
taylor_shift_conv(f, c)
```

The Horner implementation is specialised for the linear factor `x + c`.  It updates
one output coefficient array in descending order, so multiplying an intermediate
polynomial by `x + c` requires no temporary polynomial allocation.  It is quadratic
in the degree but has very small overhead.

The divide-and-conquer implementation constructs the linear polynomial `x + c` and
uses the existing divide-and-conquer composition engine.  Thus all polynomial
products continue to pass through the normal multiplication dispatcher, and no matrix
or rational-polynomial infrastructure is introduced.

Automatic dispatch is deliberately conservative.  A coarse local timing sweep showed
the generic-composition backend losing to the specialised Horner loop until roughly
4096 coefficients for small nontrivial shifts, so the current untuned cutoff is 4096;
shifts by 0, 1, or -1 stay on Horner.  Both algorithms remain directly callable for
testing and future tuning.

`tests/taylor_shift.c` checks a fixed expansion, randomized agreement of Horner,
divide-and-conquer, automatic dispatch and ordinary composition with `x + c`, plus
named function dispatch.  The public `max_abs_bits` implementation now also shares its
raw coefficient scan with the multiplication selectors.

## Iteration 30: simple REPL loops and timing

The REPL now has a deliberately small counted-loop form:

```text
smallcas> a = 0
smallcas> for i = 1:5; a = a + i
smallcas> a
15
```

The bounds are ordinary `ZZ` expressions, the range is inclusive and ascending, and
the loop body is one ordinary SmallCAS command.  If the lower bound exceeds the upper
bound, the body is not run.  Body results are suppressed, while assignments and other
side effects persist; the loop variable remains bound to its last value.  This keeps
benchmark loops from printing the result on every iteration.

REPL timing is controlled by:

```text
smallcas> @time on
smallcas> for i = 1:100; f*g
time: 0.012345678 s
smallcas> @time off
```

Timing uses `CLOCK_MONOTONIC` and reports elapsed wall time for evaluation of each
top-level command, excluding printing of the resulting value.  The `@time` control
commands themselves are not timed.  A `for` timing is for the complete loop, so it is
convenient for coarse repeated-operation benchmarks.  Since SmallCAS does not yet keep
a parsed AST, the loop timing also includes reparsing the body on every iteration.

The loop and timing commands live in a thin REPL-command layer above the existing
Bison parser; the expression grammar and the mathematical dispatch layers are
unchanged.

## Iteration 33: exact-division policy in Toom interpolation

The Toom-3 full-product interpolation uses a direct exact division by 6.  Although
one may first halve the numerator and then divide exactly by 3, benchmarks with
GMP's `mpz_divexact_ui` showed that this requires an extra linear pass and is
slower than the direct division by 6.  Its separate exact division by 2 remains a
one-bit shift.

The Toom63 middle-product interpolation was audited for the analogous issue.  It
already contains only one exact division by 3 and two exact divisions by 2; the
quantity divided by 3 is not generally even.  Its two exact halvings are expressed
as one-bit shifts.

## Iteration 34: Toom-3 interpolation-denominator optimality

`algorithms.tex` now proves that every Toom-3 interpolation using four finite
integer evaluation points and infinity has interpolation denominator divisible by
6.  For the implemented points 0, 1, -1, 2 and infinity the denominator is exactly
6, so the choice is denominator-optimal in this class.  In particular, a linear
integer-point interpolation using only shifts for exact division is impossible.
## Iteration 36: high-product tail reversal

`mulhigh` now reverses only the highest coefficients that can contribute to the
requested high window.  For an `r`-coefficient high product it materializes at
most `min(r, len(a))` and `min(r, len(b))` input coefficients, rather than
reversing both complete operands.  The multiplication is still reduced to the
existing low product, so only the ancillary copy/reversal work changes, from
operand-length dependent work to `O(r)`.

`algorithms.tex` now states the high-product-by-reversal algorithm in this
tail-only form and records the corresponding `L(r) + O(r)` coefficient cost.



## Iteration 37: middle-product theory

`algorithms.tex` now contains a full section on polynomial middle products. It
defines the balanced `(2*n-1) by n -> n` operation, proves its relation to the
transpose of multiplication by the reversed short operand, derives the block
transposition identity, and gives correctness and complexity proofs for the
classical basecase, Toom42 and Toom63. The exact evaluation/interpolation
sequences used by the companion implementation are recorded. Odd Toom42 sizes
and the two nonzero residue classes modulo 3 for Toom63 are proved correct by
explicit row/diagonal corrections, with their linear overhead quantified.

## Iteration 38: classical and divide-and-conquer division theory

`algorithms.tex` now begins the division material with classical polynomial
long division, classical truncated series quotient, divide-and-conquer series
quotient using one middle product per recursion node, conversion of ordinary
quotient computation to series division by reversal, and full-product remainder
recovery. Correctness, coefficient-operation complexity, coefficient-growth
bounds and bit-complexity bounds are proved throughout.

The reversal analysis also exposed a small implementation inefficiency:
`sc_zz_poly_quo_dc` now reverses only the highest `min(divisor_length,
quotient_length)` coefficients of the divisor, since no lower divisor
coefficient can affect the quotient. This makes the ancillary reversal/copy
work linear in the quotient length even for a very long divisor.

## Iteration 39: Mulders short division and bidirectional exact division theory

`algorithms.tex` now treats the balanced Mulders short-division recursion in
full detail.  It derives the high quotient block from a synthetic full division,
proves how the first remainder repairs the coefficients omitted from that
synthetic dividend, derives the second recursive short-division problem (with
the even/odd split handled uniformly), and proves the recurrence and bit
complexity.  The reduction from a quotient of length at most the divisor length
to a balanced `(2*n-1) / n` problem is also proved.

The exact-division subsection proves that, after removing the divisor's
valuation at `x`, the low and high halves of an exact quotient can be computed
independently by forward and reversed series division.  The only unverified
equations are a contiguous central window of exactly `divisor_length - 1`
coefficients, which is checked by one middle-product window.

The theory audit exposed the same avoidable reversal work previously fixed in
`quo_dc`: `sc_zz_poly_quo_bidirectional` now reverses only the highest
`floor(quotient_length/2)` dividend coefficients and the highest
`min(divisor_length, floor(quotient_length/2))` divisor coefficients for the
high-half quotient.  A regression with a 257-coefficient divisor and a
65-coefficient quotient covers this case.

## Iteration 40: Newton reciprocal and division theory

`algorithms.tex` now covers series inversion and the Newton division family.
It derives the classical reciprocal recurrence for a unit constant term, proves
the truncated Newton reciprocal step using one middle-product error block and
one low correction product, and gives both asymptotic and smooth leading-constant
complexities.

The section then treats reusable preinverses, showing that after a reciprocal is
precomputed each further quotient is one low product.  For a one-off quotient it
derives the Karp--Markstein construction used by the companion implementation:
compute only a half-precision reciprocal, obtain the low quotient block, form the
next residual block with a middle product, and obtain the high quotient block with
a second low product.  Ordinary polynomial quotient follows by reversal.  Newton
quotient/remainder recovery is proved to require only the low `divisor_length - 1`
coefficients of the quotient-divisor product rather than a full product.  Worst-case
coefficient-growth and bit-complexity bounds are included.

The reversal audit found the same avoidable work as in `quo_dc`:
`sc_zz_poly_quo_newton` now reverses only the highest
`min(divisor_length, quotient_length)` divisor coefficients.  A regression with a
257-coefficient divisor and 65-coefficient quotient covers this case.


## Iteration 41: classical and fast pseudo-division theory

`algorithms.tex` now treats pseudo-division over `ZZ[x]`.  It proves the
classical fraction-free elimination invariant, its quadratic degree complexity,
and coefficient/bit growth.  It then derives the scaled reciprocal
`B*I_m = lc(B)^m (mod x^m)`, proves integrality and coefficient growth directly,
and proves the scaled Newton precision-doubling step.

Reversal plus the scaled reciprocal gives the fast pseudo-quotient.  The
document proves the exact factor `lc(B)^(p-d)` that must be removed when the
Newton precision `p` is the next power of two above quotient length `d`, and
shows that only a low product of length `divisor_length - 1` is needed to
recover the pseudo-remainder.  Coefficient-operation and bit-complexity bounds
are given in both two-parameter and balanced forms.

The audit also tightened `sc_zz_poly_pseudodiv_fast`: its quotient stage now
reverses only the highest `d` dividend coefficients and the highest
`min(divisor_length, p)` divisor coefficients.  A regression with a
257-coefficient divisor and a 19-coefficient pseudo-quotient checks this case.

## Iteration 42: primitive and subresultant PRS normalization theory

`algorithms.tex` now follows the pseudo-division section with a treatment of
coefficient control across repeated pseudo-divisions.  It distinguishes the
primitive PRS, which removes the full content of each pseudo-remainder, from
Brown's subresultant PRS, which removes explicitly predictable systematic
factors without coefficient-gcd computations.

The section states the Brown--Traub fundamental theorem of subresultants in the
form needed for PRSs, derives Brown's recurrence for the normalizing factors,
proves that all resulting scalar divisions are exact and that the normalized
remainders are actual subresultants, and proves a determinant/Hadamard bound of
`O(N (tau + log N))` bits for every coefficient in the subresultant PRS.
It also makes explicit that primitive normalization is coefficient-minimal
among integral scalar associates, whereas Brown's factors optimize the tradeoff
between predictable normalization and content-computation cost and may leave
accidental content.

## Iteration 43: primitive and subresultant GCD theory

`algorithms.tex` now develops polynomial gcds over `ZZ[x]` from the PRS theory
of iteration 42.  It proves Gauss's lemma, multiplicativity of content, the
content/primitive decomposition of the canonical gcd, and correctness of the
primitive pseudo-Euclidean algorithm.  It then proves that the last nonzero
Brown-normalized subresultant is a rational associate of the polynomial gcd,
so taking its primitive part and restoring the input-content gcd gives the
canonical gcd in `ZZ[x]`.

Both classical PRS-based gcd algorithms are shown to use `O(N^2)` polynomial
coefficient operations.  The normalized remainders satisfy the subresultant
determinant bound `O(N (tau + log N))` bits.  Direct bit bounds also account
for the larger raw pseudo-remainders materialized before normalization.  The
comparison makes explicit that primitive PRS pays repeated coefficient-gcd
computations, whereas the subresultant PRS replaces them by exact divisions by
known Brown factors.

The implementation audit also added an internal remainder-only pseudo-division
kernel.  Primitive gcd and the subresultant PRS no longer construct and discard
pseudo-quotients; XGCD still uses full pseudo-division because its cofactor
recurrences require the quotient.


## Iteration 44: extended GCD and subresultant Bezout certificates

`algorithms.tex` now develops the extended subresultant PRS used by `xgcd`.
It explains why `ZZ[x]` is not a Bezout domain, proves the propagated cofactor
invariant and the exact Brown normalization of the cofactor pairs, and identifies
the returned polynomial as a subresultant Bezout certificate rather than, in
general, the canonical gcd.  In the rationally coprime case the terminal
certificate is an integer in `(A,B) cap ZZ`; its relation to the resultant and
to the least positive integer in that ideal is made explicit.  Correctness,
degree bounds, coefficient growth and coefficient/bit complexity are proved.

A regression now records the basic non-Bezout example `gcd(2,x)=1` while
`xgcd(2,x)` necessarily returns a nonunit integer certificate.

## Iteration 45: fraction-free half-GCD theory

`algorithms.tex` now develops the experimental integral half-GCD path in
matrix form.  It proves the pseudo-Euclidean transformation invariant, the
high-half transfer lemma, the two recursive half-GCD calls and crossing step,
the half-degree reduction, and the matrix-degree bound needed to justify the
truncations.  It also proves that dividing a transformation matrix by the
common integer content of all four entries preserves the represented rational
transformation up to a nonzero scalar.

The coefficient-operation recurrence is analyzed separately from coefficient
growth.  Common-content normalization is not Brown/subresultant normalization,
so the document deliberately gives only a conditional bit-complexity bound in
terms of the actual matrix coefficient size rather than claiming determinant-
size intermediates.

The audit corrected the randomized HGCD test, which had accidentally been
calling the experimental `gcd_lr` path.  Direct tests now exercise
`gcd_hgcd` and check the transformation determinant, degree reduction and gcd
preservation.

## Iteration 46: resultant theory

`algorithms.tex` now treats resultants over `ZZ[x]`: the Sylvester determinant,
root-product formula, symmetry, scaling, multiplicativity, the common-root
criterion, Euclidean and pseudo-remainder recurrences, and Brown's terminal
subresultant invariant.  It proves that a positive-degree terminal
subresultant gives resultant zero, whereas for a constant terminal member the
auxiliary Brown scalar `h` is the zeroth subresultant and hence the resultant.
Contents and the swap sign are then restored for arbitrary integer inputs.

A regression records a primitive quadratic pair for which the last constant
PRS member is `5` but the resultant and terminal Brown scalar are `25`; this
protects the distinction between the last polynomial remainder and the actual
zeroth subresultant.  Coefficient growth, output size, coefficient complexity
and direct bit complexity are included.

## Iteration 47: discriminants and squarefreeness

`algorithms.tex` now derives the discriminant from the resultant, proves that
`lc(f)` divides `resultant(f, derivative(f))` universally over the integer
coefficient ring, and proves the root-difference formula, repeated-root
criterion, scaling law and product law.  The implementation convention is
made explicit: zero and nonzero constant polynomials have discriminant zero,
while a linear polynomial has discriminant one.

Squarefreeness is defined over `QQ[x]`, so nonzero integer content is ignored.
The derivative-gcd criterion is proved, as is the structure of the squarefree
part as the positive-leading primitive product of the distinct irreducible
factors.  The zero polynomial is now correctly treated as not squarefree;
`squarefree_part(0)` remains zero, and a nonzero constant has squarefree part
one.  Coefficient and bit-complexity consequences are related directly to the
resultant, gcd and exact-division bounds from the preceding sections.

## Iteration 48: powering, evaluation, composition and Taylor shift theory

`algorithms.tex` now finishes the multiplication section with binary
polynomial powering, including exact multiplication counts, degree and
coefficient-growth bounds.  The binary-power implementation skips the leading
bit, avoiding the previous trivial square of one and multiplication by the
base; multiplication tests compare it with repeated classical products.

A new section documents integer evaluation by Horner and balanced/Estrin
schemes, polynomial composition by Horner and divide-and-conquer schemes, and
Taylor shift by the specialized in-place Horner method and by balanced
composition with `x + c`.  Correctness invariants, exact elementary operation
counts, length complexity, coefficient growth and direct bit-complexity bounds
are included.  The faster matrix-based general composition algorithms are
not used here; they are left for later linear-algebra infrastructure.

## Iteration 49: formal differentiation theory and direct higher derivatives

`algorithms.tex` now includes a formal differentiation section covering both the
coefficientwise first derivative and the direct nth derivative.  The latter uses the
falling-factorial recurrence

```text
(i + 1)_k = (i)_k * (i + 1) / (i + 1 - k)
```

with exact division, instead of recomputing `binomial(i, k)` for every output
coefficient.  This gives one pass over the surviving coefficients and avoids constructing
all lower derivatives.  The section proves correctness, exact operation counts,
coefficient growth and bit-complexity bounds, and compares the direct method with repeated
first differentiation.

The implementation now follows that recurrence.  `tests/calculus.c` additionally compares
the direct nth derivative with repeated first differentiation for every order from 0 through
25 on a degree-24 polynomial.

## Iteration 54: FFT wraparound middle products

The FFT chapter now proves the cyclic-wraparound identity for a balanced size-`n`
middle product: with cyclic length at least `2*n - 1`, every aliased coefficient lands
below the requested degrees `n - 1` through `2*n - 2`.  It gives the resulting
three-transform algorithm, coefficient bound and bit complexity, and compares the
smaller complete cyclic transform with ordinary TFT multiplication.

The implementation adds independent CRT-NTT and SSA wraparound middle-product kernels.
Tests compare both exact backends with the classical middle product across transform
boundaries and multi-prime CRT cases.  The balanced dispatcher now includes these kernels;
`make tune` measures their crossovers after tuning the classical/Toom42/Toom63 chain.


## Iteration 57: production GCD tuning

Brown's subresultant PRS previously called the classical pseudo-remainder kernel
directly, bypassing the tuned pseudo-division machinery.  Pseudo-remainder now has its
own classical/fast cutoff; the fast side uses the scaled-Newton pseudo-division and
discards the quotient.  Primitive pseudo-Euclidean gcd and Brown subresultant gcd also
have a bounded small-size dispatcher crossover.  Resultants inherit the pseudo-remainder
choice because they use the same Brown PRS.

No new division mathematics was added to `algorithms.tex`: low-product remainder
recovery, Mulders division, Newton/Karp--Markstein division, bidirectional exact division
and fast pseudo-division were already documented there.  The changes in this iteration
are dispatch and tuning choices.  The experimental fraction-free and quotient-boot HGCD
paths retain their existing base cutoff and are not part of the production gcd/xgcd
dispatch.

## Iteration 58: resultant/subresultant tuning

The public resultant now has a tuned small Sylvester/Bareiss base case before
Brown's subresultant PRS.  The cutoff metric is the Sylvester order
`deg(A) + deg(B)`, and `make tune` searches only through order 32 with the same
bounded paired timing used for the other algebraic crossovers.  The tracked
fallback uses Bareiss only at small order; Brown remains the asymptotic path.

There is no separate subresultant-PRS arithmetic cutoff to tune: every Brown
pseudo-remainder already goes through the pseudo-remainder dispatcher tuned in
iteration 57, so large PRS steps automatically inherit the fast division and
multiplication chain.  Random resultant tests now compare the public dispatcher,
the Bareiss base case and Brown's PRS against the independent Sylvester
reference determinant.


## Iteration 59: convolution Taylor shift

`ZZ[x]` now has an exact one-convolution Taylor-shift implementation, exposed as
`taylor_shift_conv(f, c)`.  It uses the identity

```text
u[n - 1 - i] = a[i] * i!
v[k] = c^k / k!
w = u * v
b[j] = w[n - 1 - j] / j!
```

but performs the divisions by factorials inside a Fermat transform ring rather
than clearing denominators by `(n - 1)!`.  The transform-root condition implies
that every factorial needed is a unit modulo `2^K + 1`; `K` is also chosen large
enough that centered reconstruction recovers the integer result uniquely.  The
FFT chapter of `algorithms.tex` gives the algorithm, the factorial-unit lemma,
correctness and bit-complexity proofs.

The automatic `taylor_shift` dispatcher is intentionally unchanged in this
iteration.  Experiments show that the convolution/D&C comparison is strongly
phase-dependent when the transform length or Fermat modulus crosses a power of
two, so a single server-derived cutoff would be misleading.  The named method
allows direct benchmarking while preserving the existing default behaviour.

## Iteration 60: phase-aware Taylor-shift dispatch

The convolution Taylor shift has two visible dyadic costs: the transform length and the
Fermat modulus exponent both round upward to powers of two.  The divide-and-conquer
composition path also has strong dyadic phases, so a single unconditional crossover is a
poor model.  Automatic dispatch therefore considers convolution only when the polynomial
length lies in the central part of its dyadic band, `5/8 <= n/P <= 13/16`, where `P` is
the least power of two at least `n`, and when the exact reconstruction-bit requirement
fills at least `5/8` of the selected Fermat modulus.  Outside those phases the existing
Horner/divide-and-conquer path is retained.

`make tune` keeps this phase rule fixed and only tunes the minimum size at which it is
worth enabling it.  It benchmarks `n = 3P/4` for `P = 512, 1024, 2048, 4096`, using
32-bit coefficients and shift 10, and stops at the first 5% convolution win.  Thus the
new stage has four bounded candidate points and cannot grow into an open-ended search.
If no win is found through `n = 3072`, automatic convolution is disabled; the named
`taylor_shift_conv` operation remains available.
