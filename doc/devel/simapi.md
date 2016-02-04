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

An SSA simulator engine may be comprised of an SSA selector and SSA process system, as described below.

### Types

 name      | type        | description 
-----------|-------------|--------------------------------------
`S::count_type` | integral type | represents population counts


### Constants

 name      | type        | description 
-----------|-------------|--------------------------------------
`S::max_process_order` | unsigned integral type | maximum order of an (elementary) reaction
`S::max_participants` | unsigned integral type | maximum number of distinct populations affected by one reaction
`S::max_instances` | unsigned integral type | maximum number of simulation instances
`S::dynamic_range` | unsigned integral type | maximum (base 2) logarithm of ratios of propensities

### Methods

In the following,

*  `g` represents a uniform random generator, passed by reference,
*  `t` is of type `double`,
*  `s` is an unsigned integral value representing a species index,
*  `c` is an unsigned integral value representing a cell index,
*  `k` is a population count of type `S::count_type`,
*  `n` is an unsigned integral value representing a number of instances,
*  `j` is an unsigned integral value representing a particular instance.

expression | return type | description
-----------|-------------|--------------------------------------
`s.initialise(n,M,t0)` |  | set up simulator to simulate `n` instances of the given `rd_model` `M` with initial simulation time `t0`
`s.count(s,c,j)` | `count_type` | population count of species index `s` in cell `c` in instance `j`
`s.count(s,c)`   | `count_type` | equivalent to `s.count(s,c,0)`
`s.set_count(s,c,k,j)` |  | set population count of species index `s` in cell `c` to `k` in instance `j`
`s.set_count(s,c,k)`   |  | equivalent to `s.set_count(s,c,k,0)`
`s.advance(g)` | double   | *[optional]* advance simulator state by minimum time step, returning new simulation time
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
processes, and maintains the state of populations and processes across one or more
independent instances.

For a process system of type `Y`, instance `y`.

### Types

 name      | type        | description 
-----------|-------------|--------------------------------------
`Y::key_type`   | implementation specific | keys representing individual procesess
`Y::value_type` | floating point type | represents process propensities
`Y::pop_type`   | integral type | represents population indices
`Y::count_type` | integral type | represents population counts

As for an SSA selector, `Y::key_type` and `Y::pop_type` should likely be unsigned integral types.

Individual processes are described by a process description `p` which offers methods
according to the following:

 expression | type | description
------------|------|-------------------------------------
`p.left()` | collection of `Y::pop_type` |  (multi-)set of population indices consumed by process
`p.right()` | collection of `Y::pop_type` |  (multi-)set of population indices produced by process
`p.rate()` | `Y::value_type` | elementary process rate


### Constants

 name      | type        | description 
-----------|-------------|--------------------------------------
`Y::max_population_index` | unsigned integral | maximum population index
`Y::max_process_order` | unsigned integral | maximum process order
`Y::max_participants` | unsigned integral | maximum disrinct populations involved in a process
`Y::max_count` | unsigned integral | maximum population count
`Y::max_instances` | integral type | maximum number of system instances

### Methods

In the following,

* `k` is of type `Y::key_type`, representing a process identifier,
* `p` is an unsigned integral value referencing a population,
* `c` is of type `Y::count_type`, representing a population count,
* `q` is a process description,
* `b`,`e` form an iterator range of process descriptions,
* `n` is an unsigned integral value representing a number of instances,
* `j` is an unsigned integral value representing a particular instance,
* `notify` is a function or functional object with signature equivalent to `void notify(key_type)`.

expression | return type | description
-----------|-------------|--------------------------------------
`y.initialise(n)` |          | initialise state for `n` instances and 
`y.clear()` |                | remove all processes
`y.add(q)` |                 | add process with description `q`
`y.add(b,e)` |               | add processes described by iterator interval [`b`,`e`)
`y.define_processes(b,e)` |  | configure system with processes described by iterator interval [`b`,`e`)
`y.size()`  | `size_t`       | number of processes in system
`y.n_instances()` | `size_t` | number of instances
`y.reset()` |                | zero population counts across all instances
`y.propensity(j,k)` | `Y::value_type` | calculate propensity for process `k` in instance `j`
`y.count(p,j)`  | `Y::count_type` | population count for population `p` in instance `j`
`y.count(p)`  | `Y::count_type` | equivalent to `y.count(p,0)`
`y.set_count(p,c,notify,j)` | | set count for population `p` to `c` in instance `j`; call `notify(u)` for each affected process `u`
`y.set_count(p,c,notify)` |  | equivalent to `y.set_count(p,c,notify,0)`
`y.apply(k,notify,j)` |      | apply process `k` to state of instance `j`; call `notify(u)` for each affected process `u`.
`y.apply(k,notify)` |        | equivalent to `y.apply(k,notify,0)`

As for an SSA selector, `Y::key_type` should likely be an unsigned integral type.
Adding a process to a process system may or may not preserve population counts — this is a quality
of implementation issue.
</div>
