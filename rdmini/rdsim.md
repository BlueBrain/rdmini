# Simulator engine interface

Any instance `s` of a simulator implementation `S` should provide the following interface to callers.

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
`a.next(g)` | A::event_type | generate next event drawing uniformly distributed numbers from `g`
`ev.key()`  | A::key_type | identifier of process in event
`ev.dt()`   | floating point type | event time delta

In practice, `A::key_type` should probably be an unsigned integral type, taking
values from the range [0, `a.size()`).

