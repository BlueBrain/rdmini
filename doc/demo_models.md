---
title: Example reaction–diffusion models
---

The `demo/models.yaml` contains a number of simple reaction system and reaction–diffusion
system models from the literature.

Note that the units of the quantities expressed in the models file are in counts, metres,
and seconds; this leads to some very large or small quantities for typical cell volumes,
such as the 25 µm^3 cells used in the Turing example below. This problem will be
addressed in a future revision that includes general unit support. For single well-mixed
volume simulations, a 'unit' volume allows rate constants to be considered to be in s<sup>-1</sup>.


## Schnakenberg

An example demonstrating limit cycle behaviour, first raised in [@schnakenberg79, p. 398].
The rate constants and examples come from [@erban2007, p. 30], and have been chosen to
demonstrate stochastic behaviour that differs qualitatively from that of the corresponding
deterministic system of ODEs.

**Reaction system**

| reaction   | rate constants |
|------------|----------------|
| 2*A* + *B* → 3*A* | *k*<sub>1</sub>=4×10<sup>-5</sup> s<sup>-1</sup>             |
| ∅ ⇄ *A*           | *k*<sub>2</sub>=40 s<sup>-1</sup>,  *k*<sub>3</sub>=10 s<sup>-1</sup> |
| ∅ → *B*           | *k*<sub>4</sub>=25 s<sup>-1</sup>                   |

**Initial conditions**

| species | concentration |
|---------|---------------|
| *A*: | 10 |
| *B*: | 10 |

## Schlögl

Single species model reaction model exhibiting bistable behaviour, as introduced by
Schlögl in [@schlögl82, p. 150] and discussed in [@erban2007, p. 29]. The system
exhibits two stable steady states with concentrations 100 and 400; stochastic
simulation demonstrates spontaneous transition between these two states.

**Reaction system**

| reaction   | rate constants |
|------------|----------------|
| 2*A* ⇄ 3*A* | *k*<sub>1</sub>=0.18 s<sup>-1</sup>, *k*<sub>2</sub>=2.5×10<sup>-4</sup> s<sup>-1</sup> |
| ∅ ⇄ *A*     | *k*<sub>3</sub>=2200 s<sup>-1</sup>, *k*<sub>4</sub>=37.5 s<sup>-1</sup>       |

Note that this description has rescaled the reaction constants from per-minute
quantities to per-second quantities.

**Initial conditions**

| species | concentration |
|---------|---------------|
| *A*: | 0 |

## Turing

A reaction–diffusion model exhibiting spontaneous emergence of spatial inhomogeneity,
based on Turing patterns [@turing52]. The example constructed in [@erban2007, p. 32]
uses the Schnakenberg reaction system with rate constants described below, on a
1-D mesh of length 1 mm, divided into 40 equal sized compartments.

Note that the description in the paper expresses concentrations as counts per
mesh element; the model implementation has 40 cells each of volume 25 µm<sup>3</sup>,
and the reaction rates below have been scaled accordingly.

**Reaction system**

| reaction   | rate constants |
|------------|----------------|
| 2*A* + *B* → 3*A* | *k*<sub>1</sub>=6.25×10<sup>-4</sup> µm<sup>6</sup>s<sup>-1</sup>               |
| ∅ ⇄ *A*           | *k*<sub>2</sub>=0.04 µm<sup>3</sup>·s<sup>-1</sup>,  *k*<sub>3</sub>=0.02 s<sup>-1</sup> |
| ∅ → *B*           | *k*<sub>4</sub>=0.12 µm<sup>3</sup>·s<sup>-1</sup>                     |


**Diffusion coefficients**

| species | diffusivity |
|---------|-------------|
| *A*: | 10<sup>-5</sup> mm<sup>2</sup>·s<sup>-1</sup> |
| *B*: | 10<sup>-3</sup> mm<sup>2</sup>·s<sup>-1</sup> |

**Initial conditions**

| species | concentration |
|---------|---------------|
| *A*: | 8 µm<sup>-3 |
| *B*: | 0.12 µm<sup>-3 |


-----

## References

---
author:
- name: Sam Yates
  affiliation: Blue Brain Project, EPFL

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
