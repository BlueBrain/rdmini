# Simulator engine interface

Any instance `s` of a simulator implementation `S` should provide the following interface to callers.

## Types

 name | type | description 
------|------|-------------
 `count_type` | integral type | represents population counts

## Constants

 name | type | description
------|------|-------------
`max_order` | unsigned integral type | maximum order of an (elementary) reaction
`max_participants` | unsigned integral type | maximum number of distinct populations affected by one reaction
`dynamic_range` | unsigned intgral type | maximum (base 2) logarithm of ratios of propensities

## Methods

expression | return type | description
-----------|-------------|------------
`s.initialise(M,t0)` | | set up simulator to simulate given `rd_model` with initial sim time `t0`
`s.count(s_id,c_id)` | `count_type` | return population count of species index `s_id` in cell `c_id`
`s.set_count(s_id,c_id,count)` | | set population count of species index `s_id` in cell `c_id`
`s.advance()` | double | advance simulator state by minimum time step, returning new simulation time
`s.advance(t)` | double | advance simulator past time `t`, returning new simulation time

