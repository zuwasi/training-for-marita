# Chapter 1: Wolfram Language Overview

## Core Idea
Mathematica is a front end (notebook) talking to a kernel that holds state. Everything is an expression; assignments are rewrite rules checked in order; you can freely mix procedural, functional, and rule-based styles. For quantum mechanics, the essentials are: control *when* right-hand sides evaluate (`=` vs `:=`, `->` vs `:>`), represent states/operators as (sparse) lists, diagonalize Hermitian matrices, and state assumptions about real variables. Strip units into dimensionless numbers before numerics.

## Frameworks Introduced
- **Immediate vs delayed assignment**: use `=` to evaluate the RHS once, now; use `:=` to store the RHS and re-evaluate it on every use. The same split applies to functions and rules.
- **Memoization** (functions that remember their results): `F[x_] := F[x] = ...`. Each call adds a specific immediate definition `F[2]=128`, which is found before the general pattern. Use it for expensive or recursive functions. Put `Clear[F]` before the definition, and persist values with `DumpSave`/`Get`.
- **Patterns with conditions**: anything-goes `x_`, type-restricted `x_Integer`, conditional `x_ /; x<=3`, test `x_?PrimeQ`, optional `b_:5`.
- **Rules and replacements**: `/.` applies once, `//.` repeats until nothing changes, and `Replace[..., {1}]` targets one level.
- **Many ways to define the factorial function**: 17 definitions covering alias, `Gamma`, `Product`, rule-based recursion, `If` recursion, `#0` anonymous recursion, `Do`/`For`, `Times @@ Range[n]`, `Fold`, `Permutations`, `//.`, strings, and `D[x^n,{x,n}]`. Choose whichever paradigm fits how you think.
- **Dimensionless multipliers**: choose a consistent unit set, divide each `Quantity` by its unit, and compute with plain numbers.

## Key Concepts
- **Pattern** `x_`: a named blank that matches any expression. It is not the symbol `x`.
- **Definition order**: specific definitions (`f[0]=1`) are tried before general patterns (`f[x_]`).
- **Anonymous function**: `#^2&`. Here `#`=`#1`, `##` means all arguments, and `#0` is the function itself. Functions are first-class objects.
- **Map / Apply**: `f /@ list` works element by element; `f @@ list` replaces the head (a reduce).
- **Three-valued logic**: `x==0` with undefined `x` stays unevaluated, so the 4th slot of `If` handles "unknown".
- **FullForm**: shows the internal tree. `a/b` is `Times[a,Power[b,-1]]` but `2/3` is `Rational[2,3]`.
- **SparseArray**: stores only nonzeros, and products scale with the nonzero count.
- **Assumptions**: by default every variable is complex. Supply realness via `ComplexExpand`, `Assuming`, or `$Assumptions`.
- **Protected symbols**: built-ins need `Unprotect` before you redefine them.

## Mental Models
- Use `=` when the RHS should be frozen, or when it simplifies symbolically into something faster (`f[n_] = Product[i,{i,n}]` gives `n!`). Use `:=` whenever the RHS needs the argument's value, loops, `Range[n]`, randomness, or recursion.
- Immediate is a hard link; delayed is a symbolic link.
- Use `:>` when the rule's RHS must run per match (side effects, `Table[m-1,{m}]`).
- Use `Module[{t},...]` or `Do` for local loop variables. `For`/`While` counters leak into `Global``.
- Use `Echo[x,label]` to inspect a value mid-pipeline without changing the result. Use `Trace`/`FullForm` when a pattern fails to match.

## Anti-patterns
- **`=` with an argument-dependent RHS**: `f[n_] = Module[...For[i=1,i<=n,...]]` runs at definition time, while `n` is unbound, so you always get 1. Recursion likewise breaks.
- **Stale memoized values**: redefining `F` without `Clear[F]` keeps the old cached `F[2]`.
- **Not clearing between redefinitions**: old `f[1]=1` or alias definitions linger. Use `Clear[f]` (or `ClearAll` to remove attributes too).
- **Assuming notebook order is kernel order**: cells execute in the order you press shift-enter, so kernel state can disagree with what you see on screen.
- **Leaked globals**: a global `i=6` left over from `For` silently pollutes later symbolic code. `x=5` turns later `Trace[x-3x+1]` into a number.
- **Assuming eigenvalues are sorted**: their order is undefined, and eigenvector phases are arbitrary.
- **Matching `x_/y_` on `2/3`**: it fails, because `Rational` is atomic.
- **Units in heavy numerics**: they are slow, can need Wolfram|Alpha online, and do not work in matrix algebra.
- **Precedence surprises**: `Sqrt@2+3` gives `3+Sqrt[2]`. When unsure, use brackets.

## Code Examples
```mathematica
Clear[F];
F[x_] := F[x] = x^9          (* memoization *)
SetDirectory[NotebookDirectory[]];
DumpSave["Fdef.mx", F];      (* later: Clear[F]; redefine; Get["Fdef.mx"] *)
g[1] = g[2] = 1;
g[n_] := g[n] = g[n-1] + g[n-2];   (* memoized Fibonacci: dramatic speedup *)
Table[g[n], {n, 100}]
```
What it demonstrates: memoization through the interplay of `:=` and `=`, plus persisting cached values.

```mathematica
f[n_] := Times @@ Range[n]                          (* most elegant: says what, not how *)
f[n_] := Fold[Times, 1, Range[n]]
f = If[#1 == 1, 1, #1*#0[#1-1]]&;                   (* recursive anonymous function *)
f[n_] := First[{1,n} //. {a_,b_/;b>0} :> {b*a,b-1}] (* repeated rule replacement *)
```
What it demonstrates: four paradigms for the same function. Call `Clear[f]` between them.

```mathematica
i = 0;
{a, b, c, d, c, a, c, b} /. c :> (i++; Echo[i, "replacement "]; f)
a //. {a -> b, b -> c}                    (* c; with /. you get b *)
Replace[a + b, x_ -> x^2, {1}]            (* a^2 + b^2 *)
x^2 - 4 /. Solve[x^2 - 4 == 0, x]         (* {0, 0} *)
```
What it demonstrates: delayed rules with side effects, `/.` vs `//.`, level-specific replacement, and using `Solve` output as rules.

```mathematica
H = {{0, 0.3, I, 0}, {0.3, 1, 0, 0}, {-I, 0, 1, -0.2}, {0, 0, -0.2, 3}};
{eval, evec} = Transpose[Sort[Transpose[Eigensystem[H]]]];
{eval, evec} = Transpose[Sort[Transpose[-Eigensystem[-H, 2,
   Method -> {"Arnoldi", "Criteria" -> "RealPart", MaxIterations -> 10^6}]]]];
Conjugate[evec].Transpose[evec] // Chop // MatrixForm     (* identity: orthonormal *)
Conjugate[evec].H.Transpose[evec] // Chop // MatrixForm   (* diagonal = eval *)
```
What it demonstrates: getting ascending (eigenvalue, eigenvector) pairs, computing only the lowest k via `-H`, and checking the result.

## Reference Tables
| Immediate | Delayed | Meaning |
|---|---|---|
| `a = RandomReal[]` | `b := RandomReal[]` | RHS evaluated once vs on every use |
| `f[x_] = ...` | `f[x_] := ...` | function body evaluated at definition vs at each call |
| `x -> y` | `x :> y` | rule RHS evaluated once vs per match |

| Brackets | Use |
|---|---|
| `( )` | grouping `2*(3-7)` |
| `[ ]` | function call `Sin[0.2]` |
| `{ }` | list `{a,b,c}` |
| `[[ ]]` | indexing `v[[2]]`, `M[[1,3]]`, `v[[-1]]` (1-based) |

| Form | Syntax | Precedence |
|---|---|---|
| standard | `Sqrt[2+3]` | unambiguous |
| prefix | `Sqrt @ 2+3` gives `3+Sqrt[2]` | high (above `*`) |
| postfix | `2+3 // Sqrt` | low (below `+`); `//N`, `//MatrixForm`, `//Timing` |

| Sparse / tensor | Purpose |
|---|---|
| `SparseArray[{{1,2}->3,{2,4}->-1},{3,10}]` | build from nonzeros |
| `ArrayRules[Ms]`, `Normal[Ms]` | inspect, convert to dense |
| `ArrayReshape[v,{2,3,4}]`, `Flatten[u,{{1,2},{3,4}}]` | reshape, partial flatten |
| `ArrayFlatten[u]` | same as `Flatten[u,{{1,3},{2,4}}]` |
| `Transpose[t,{2,3,1}]`, `TensorContract[u,{1,2}]` | permute indices, contract equal-length indices |

| Function | Use when |
|---|---|
| `Eigenvalues[H]` | energies only. Output is unsorted, so wrap it in `Sort` |
| `Eigenvalues[H, k, Method->{"Arnoldi",...}]` | you need the k largest (use `-Eigenvalues[-H,k,...]` for the smallest) |
| `Eigenvectors[H]` | vectors only |
| `Eigensystem[H]` | paired values and vectors (the usual choice) |

## Key Takeaways
1. Choosing `=` or `:=` (and `->` or `:>`) is a choice about *when evaluation happens*. Most bugs come from evaluating too early.
2. Use `F[x_] := F[x] = ...` for memoization, and always `Clear[F]` first.
3. Prefer functional idioms (`/@`, `@@`, `Range`, `Fold`, `Total`) for speed and clarity. Keep variables local with `Module`.
4. When patterns fail, check `FullForm`.
5. Sort eigensystems explicitly, and use the `-H` Arnoldi trick for ground states.
6. Declare realness (`ComplexExpand`, `Assuming`, `$Assumptions`) before using `Conjugate`/`Im`/`FullSimplify`.
7. Do numerics in dimensionless multipliers of a consistent unit set.

## Connects To
- **Ch 2 (quantum mechanics)**: Hilbert-space states become vectors, operators become (sparse) matrices, `Eigensystem` solves the Schrödinger equation, and the tensor tools of 1.10.5 (`ArrayReshape`, `Flatten`, `TensorContract`) show up in 2.4.3 for composite systems.
- **Later numerical chapters**: `SparseArray`, Arnoldi lowest-eigenvalue calls, and the dimensionless-units convention from 1.12.
