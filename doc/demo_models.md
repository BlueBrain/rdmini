# Example reaction–diffusion models.

The `demo/models.yaml` contains a number of simple reaction system and reaction–diffusion
system models from the literature.

Note that the units of the quantities expressed in the models file are in counts, metres,
and seconds; this leads to some very large or small quantities for typical cell volumes,
such as the 25 µm^3 cells used in the Turing example below. This problem will be
addressed in a future revision that includes general unit support. For single well-mixed
volume simulations, a 'unit' volume allows rate constants to be considered to be in s^-1^.


## Schnakenberg

An example demonstrating limit cycle behaviour, first raised in [@schnakenberg79, p. 398].
The rate constants and examples come from [@erban2007, p. 30], and have been chosen to
demonstrate stochastic behaviour that differs qualitatively from that of the corresponding
deterministic system of ODEs.

**Reaction system**

|                  |                                            |
|------------------|--------------------------------------------|
|2*A* + *B* → 3*A* | &#x2001; *k*~1~=4×10^-5^ s^-1^             |
|∅ ⇄ *A*           | &#x2001; *k*~2~=40 s^-1^,  *k*~3~=10 s^-1^ |
|∅ → *B*           | &#x2001; *k*~4~=25 s^-1^                   |

**Initial conditions**

|     |    |
|-----|----|
|*A*: | 10 |
|*B*: | 10 |

## Schlögl

Single species model reaction model exhibiting bistable behaviour, as introduced by
Schlögl in [@schlögl82, p. 150] and discussed in [@erban2007, p. 29]. The system
exhibits two stable steady states with concentrations 100 and 400; stochastic
simulation demonstrates spontaneous transition between these two states.

**Reaction system**

|            |                                                     |
|------------|-----------------------------------------------------|
|2*A* ⇄ 3*A* | &#x2001; *k*~1~=0.18 s^-1^, *k*~2~=2.5×10^-4^ s^-1^ |
|∅ ⇄ *A*     | &#x2001; *k*~3~=2200 s^-1^, *k*~4~=37.5 s^-1^       |

Note that this description has rescaled the reaction constants from per-minute
quantities to per-second quantities.

**Initial conditions**

|     |  |
|-----|--|
|*A*: | 0|

## Turing

A reaction–diffusion model exhibiting spontaneous emergence of spatial inhomogeneity,
based on Turing patterns [@turing52]. The example constructed in [@erban2007, p. 32]
uses the Schnakenberg reaction system with rate constants described below, on a
1-D mesh of length 1 mm, divided into 40 equal sized compartments.

Note that the description in the paper expresses concentrations as counts per
mesh element; the model implementation has 40 cells each of volume 25 µm^3^,
and the reaction rates below have been scaled accordingly.

**Reaction system**

|                  |                                                      |
|------------------|------------------------------------------------------|
|2*A* + *B* → 3*A* | &#x2001; *k*~1~=6.25×10^-4^ µm^6^s^-1^               |
|∅ ⇄ *A*           | &#x2001; *k*~2~=0.04 µm^3^·s^-1^,  *k*~3~=0.02 s^-1^ |
|∅ → *B*           | &#x2001; *k*~4~=0.12 µm^3^·s^-1^                     |


**Diffusion coefficients**

|     |                    |
|-----|--------------------|
|*A*: | 10^-5^ mm^2^·s^-1^ |
|*B*: | 10^-3^ mm^2^·s^-1^ |



-----

## References

---
references:
- id: schnakenberg79
  type: article-journal
  author:
  - family: Schnakenberg
    given: [ J. ]
  title: Simple chemical reaction systems with limit cycle behaviour
  container-title: Journal of Theoretical Biology
  volume: 81
  issued: { year: 1979 }
  page: 389–400
  DOI: 10.1016/0022-5193(79)90042-0
- id: schlögl82
  type: article-journal
  author:
  - family: Schlögl
    given: [ F. ]
  title: Chemical reaction models for non-equilibrium phase transitions
  container-title: Zeitscrift für Physik
  volume: 253
  issue: 2
  issued: { year: 1972, month: 4 }
  page: 147–161
  DOI: 10.1007/BF01379769
- id: erban2007
  type: webpage
  author:
  - family: Erban
    given: [ Radek ]
  - family: Chapman
    given: [ Jonathan ]
  - family: Maini
    given: [ Phillip ]
  title: A practical guide to stochastic simulations of reaction–diffusion processes
  issued: { year: 2007 }
  URL: http://arxiv.org/abs/0704.1908v2
- id: turing52
  type: article-journal
  author:
  - family: Turing
    given: [ A., M. ]
  title: The chemical basis of morphogenesis
  containter-title: Philosophical Transactions B
  volume: 237
  issue: 541
  issued: { year: 1952 }
  page: 37–72
  DOI: 10.1098/rstb.1952.0012
link-citations: true
citation-style: style.csl
...
