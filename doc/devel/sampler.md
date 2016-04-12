---
title: Random sampling
---

The problem: probabilistically even distribution of a value (concentration) across a discrete set of elements (cells) with
varying weights (areas or volumes), where each element must have an integer count. Suppose we distrute $n$ across $N$
elements with weights $w_i$, with total weight $w=\sum w_i$. then the goal is to sample the random variables $n_i$ such that

1.  The distribution is _fair_: $E[n_i]= nw_i/w$ .

2.  The distribution is _even_: $\text{Var}(n_i)$ is minimized.

3.  The distribution is _uncorrelated_: $\text{Cov}(n_i,n_j)=\delta_{i,j}$ is minimized.


A general approach to reducing the variation while preserving the expectation is to first allocate $n_i^{(0)}=\lfloor n w_i/w \rfloor$ to
each element, and then distribute the remaining $n'=n-\sum n_i^{(0)}$ items by a weighted sample: $n_i = n_i^{(0)} + S_i$ , where
$S=(S_1,\ldots,S_N)$ is a sample of $n'$ items with inclusion probabilities $\pi_i=n w_i/w - n_i^{(0)}$.

## Distribution routines in `demo_distribute.cc` 

Three routines are currently implemented: 'steps', 'multinomial' and 'oss'.

The 'steps' implementation is included for comparison: it mimics the routine used in the Tetexact and TetOpSplit solvers
in STEPS 0.9.1.

The 'multinomial' and 'oss' solvers are of the form described above; they allocate the integer part of the proportional
count to each element, and then draw a random sample for the difference. For the 'multinomial' sampler, this is a weighted
with-replacement sample corresponding to the multinomial distribution. For the 'oss' sampler, this is a weighted
without-replacement sample drawn using ordered systematic sampling (see [@tillé06, p. 124] 
Without-replacement samplers will give the lowest variation, as $n_i$ will be constrained to be $\lfloor n w_i/w \rfloor$
or $\lceil n w_i/w \rceil$ for each $i$. The 'oss' sampler though introduces severe correlation (its chief advantage
is its ease of implementation).

## Sampling API

The implemented samplers follow an API/concept that is broadly based on that of the random distribution concept from the
C++11 standard. All samplers sample from a population (described by an iterator range); the minimum population size and
the minimum and maximum sample sizes are determined by the sampler and its parameters.

<div class="concept-table">

Any instance `s` of a sampler implementation `S` should provide the following interface to callers.

### Types

 name      | type        | description 
-----------|-------------|----------------------------------------------
`S::size_type` | integral type | represents sample and population counts
`S::real_type` | floating point type | represents probabilities and weights
`S::param_type` | | represents parameters (e.g. weights) for sampler

### Methods

In the following,

* `g` is a uniform random number generator, passed by reference,
* `p` is of type `S::param_type`,
* `b` and `e` are input iterators describing a population by the range [`b`,`e`).
    Note that depending on the sampler, these may be required to be forward or
    random access iterators.
* `o` is an output iterator. Depending on the sampler, this may be required to
    be a random access iterator.

expression | return type | description
-----------+-------------+----------------------------------------------
`s.param()`   | `S::param_type`  | retrieve sampler parameters
`s.param(p)`  | `void`  | set sampler parameters
`s.min()`     | `S::size_type` | minimum sample size
`s.max()`     | `S::size_type` | maximum sample size
`s.size()`    | `S::size_type` | minimum population size to draw from
`s.sample(b,e,o,g)` | `S::size_type` | sample from range [`b`,`e`) to `o`, returning sample count

</div>


----

## References

---
references:
- id: tillé06
  type: book
  author:
  - family: Tillé
    given: Yves
  title: Sampling Algorithms
  issued: { year: 2006 }
  ISBN: 978-0387-30814-2
  publisher: Springer
  publisher-place: NY, USA

link-citations: true

citation-style: style.csl
...






