---
title: Simulator concepts and exceptions
---

## Exceptions

### `rdmini::operation_not_supported`

Derives from `std::logic_error`

Some operations described in the API below may be optional, in that they may not be practical
to support in every implementation, and are not crucial in providing SSA functionality. Such
operations may thus instead throw an exception of this type.

### `rdmini::invalid_value`

Derives from `std::runtime_error`

Implementations are not obliged to perform range checking; should they do so, however,
they must throw this exception when a supplied value is invalid.

### `rdmini::ssa_error`

Derives from `std::runtime_error`

Represents an internal error in an SSA implementation.

## Simulator engine interface

<!-- use special table rendering from css -->
<div class="concept-table">

Any instance `s` of a simulator implementation `S` should provide the following interface to callers.

An SSA implementation may be comprise an SSA selector and SSA process system, described below.

### Types

 name      | type        | description 
-----------|-------------|--------------------------------------
`S::count_type` | integral type | represents population counts
foo | bar | baz

### Constants

 name      | type        | description 
-----------|-------------|--------------------------------------
`S::max_process_order` | unsigned integral type | maximum order of an (elementary) reaction
`S::max_participants` | unsigned integral type | maximum number of distinct populations affected by one reaction
`S::dynamic_range` | unsigned intgral type | maximum (base 2) logarithm of ratios of propensities

### Methods

Here `g` represents a uniform random generator, passed by reference.

expression | return type | description
-----------|-------------|--------------------------------------
`s.initialise(M,t0)` | | set up simulator to simulate given `rd_model` with initial sim time `t0`
`s.count(s,c)` | `count_type` | return population count of species index `s` in cell `c`
`s.set_count(s,c,k)` | | set population count of species index `s` in cell `c` to `k`
`s.advance(g)` | double | advance simulator state by minimum time step, returning new simulation time
`s.advance(t,g)` | double | advance simulator up to time `t`, returning new simulation time


## SSA selector implementations

Let `A` denote a class implementing the SSA selector API.

### Types

name       | type        | description
-----------|-------------|--------------------------------------
`A::key_type` | implementation specific | keys representing individual procesess
`A::value_type` | floating point type | represents process propensities
`A::event_type` | implementation specific | describes generated process events


### Methods

In the following let `a` be an instance of `A`, `ac` a const instance of `A`,
`g` a uniform random number generator (see section [rand.req.urng] in C++ standard),
`n` an unsigned integral value, `k` a value of type `A::key_type`, `r` a value of type
`A::value_type`, and `ev` an event of type `A::event_type` returned by the
`next` method.

expression | return type | description
-----------|-------------|--------------------------------------
`a.reset(n)` | | initialise state to represent `n` processes, with initially zero propensity
`a.update(k,r)` | | set propensity of process `k` to `r`; may throw `rdmini::invalid_value`
`ac.size()`   | unsigned integral type | total number of represented processes
`ac.propensity(k)` | `A::value_type` | *[optional]* retrieve propensity of process `k`; may throw `rdmini::invalid_value`
`ac.total_propensity()` | `A::value_type` | *[optional]* retrieve propensity of process `k`; may throw `rdmini::invalid_value`
`a.next(g)` | `A::event_type` | generate next event drawing uniformly distributed numbers from `g`
`ev.key()`  | `A::key_type` | identifier of process in event
`ev.dt()`   | floating point type | event time delta

In practice, `A::key_type` should generally be an unsigned integral type, taking
values from the range [0, `a.size()`).

## SSA process system implementation

A process system encapsulates the dependency relations between populations and
processes.

For a process system of type `Y`, instance `y`.

### Types

 name      | type        | description 
-----------|-------------|--------------------------------------
`Y::key_type` | implementation specific | keys representing individual procesess
`Y::value_type` | floating point type | represents process propensities
`Y::count_type` | integral type | represents population counts

As for an SSA selector, `Y::key_type` should likely be an unsigned integral type.

### Constants

 name      | type        | description 
-----------|-------------|--------------------------------------
`Y::max_population_index` | unsigned integral | maximum population index
`Y::max_process_order` | unsigned integral | maximum process order
`Y::max_participants` | unsigned integral | maximum disrinct populations involved in a process
`Y::max_count` | unsigned integral | maximum population count

### Methods

In the follwing,

* `k` is of type `Y::key_type`, representing a process identifier.
* `p` is an unsigned integral value referencing a population.
* `c` is of type `Y::count_type`, representing a population count.
* `b`,`e` form an iterator range of process descriptions, where for iterators `i` in [`b`,`e`),
    * `i->left()` is a collection (multiset) of population indices (reactants)
    * `i->right()` is a collection (multiset) of population indices (products)
    * `i->rate()` is the elementary process rate (convertible to `Y::value_type`.)
* `notify` is a function object with signature equivalent to `void notify(key_type)`.

expression | return type | description
-----------|-------------|--------------------------------------
`y.clear()`   | | reset state, discarding all process information
`y.defing_processes(b,e)` | | configure system with processes described by iterator interval [`b`,`e`)
`y.size()  `  | | number of processes in system
`y.zero_populations()` | | zero population counts, invalidating all propensities
`y.propensity(k)` | `Y::value_type` | calculate propensity for process `k`
`y.count(p)`  | `Y::count_type` | population count for population `p`
`y.apply(k,notify)`  | | apply process `k` to state; call `notify(j)` for each affected process.
`y.set_count(p,c,notify)`  | | set count for population `p` to `c`; call `notify(j)` for each affected process.

As for an SSA selector, `Y::key_type` should likely be an unsigned integral type.

</div>
