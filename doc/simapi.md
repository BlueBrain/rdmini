# Simulator engine interface

Any instance `s` of a simulator implementation `S` should provide the following interface to callers.

An SSA implementation may be comprise an SSA selector and SSA process system, described below.

## Types

 name | type | description 
------|------|-------------
 `S::count_type` | integral type | represents population counts

## Constants

 name | type | description
------|------|-------------
`S::max_order` | unsigned integral type | maximum order of an (elementary) reaction
`S::max_participants` | unsigned integral type | maximum number of distinct populations affected by one reaction
`S::dynamic_range` | unsigned intgral type | maximum (base 2) logarithm of ratios of propensities

## Methods

Here `g` represents a uniform random generator, passed by reference.

expression | return type | description
-----------|-------------|------------
`s.initialise(M,t0)` | | set up simulator to simulate given `rd_model` with initial sim time `t0`
`s.count(s_id,c_id)` | `count_type` | return population count of species index `s_id` in cell `c_id`
`s.set_count(s_id,c_id,count)` | | set population count of species index `s_id` in cell `c_id`
`s.advance(g)` | double | advance simulator state by minimum time step, returning new simulation time
`s.advance(t,g)` | double | advance simulator past time `t`, returning new simulation time


# SSA selector implementations

For an SSA selector of type `A`, instance `a`, and uniform integral random generator `g`.

## Types

name | type | description
-----|------|------------
`A::key_type` | implementation specific | keys representing individual procesess
`A::value_type` | floating point type | represents process propensities
`A::event_type` | implementation specific | describes generated process events


## Methods

Here `n` is an unsigned integral type, `k` is of type `A::key_type`, `r` is of type
`A::value_type`, `ev` is a generated event of type `A::event_type`.

expression | return type | description
-----------|-------------|------------
`a.reset(n)` | | initialise state to represent `n` processes, with initially zero propensity
`a.size()`   | unsigned integral type | Total number of represented processes
`a.update(k,r)` | | set propensity of process `k` to `r`
`a.next(g)` | `A::event_type` | generate next event drawing uniformly distributed numbers from `g`
`ev.key()`  | `A::key_type` | identifier of process in event
`ev.dt()`   | floating point type | event time delta

In practice, `A::key_type` should probably be an unsigned integral type, taking
values from the range [0, `a.size()`).


# SSA process system implementation

A process system encapsulates the dependency relations between populations and
processes.

For a process system of type `Y`, instance `y`.

## Types

name | type | description
-----|------|------------
`Y::key_type` | implementation specific | keys representing individual procesess
`Y::value_type` | floating point type | represents process propensities
`Y::count_type` | integral type | represents population counts

As for an SSA selector, `Y::key_type` should likely be an unsigned integral type.

## Constants

name | type | description
-----|------|------------
`Y::max_population_index` | unsigned integral | maximum population index
`Y::max_process_order` | unsigned integral | maximum process order
`Y::max_participants` | unsigned integral | maximum disrinct populations involved in a process
`Y::max_count` | unsigned integral | maximum population count

## Methods

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
-----------|-------------|------------
`y.clear()`   | | reset state, discarding all process information
`y.defing_processes(b,e)` | | configure system with processes described by iterator interval [`b`,`e`)
`y.size()  `  | | number of processes in system
`y.zero_populations()` | | zero population counts, invalidating all propensities
`y.propensity(k)` | `Y::value_type` | calculate propensity for process `k`
`y.count(p)`  | `Y::count_type` | population count for population `p`
`y.apply(k,notify)`  | | apply process `k` to state; call `notify(j)` for each affected process.
`y.set_count(p,c,notify)`  | | set count for population `p` to `c`; call `notify(j)` for each affected process.

As for an SSA selector, `Y::key_type` should likely be an unsigned integral type.

