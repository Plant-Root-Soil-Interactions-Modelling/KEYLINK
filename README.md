# KEYLINK
In the KEYLINK model, we integrate new knowledge on soil structure and its importance for soil organic matter (SOM) stabilization and hydrology, with the existing concepts on SOM pools, and elements from food web models, i.e. those from direct trophic interactions among soil organisms. 

## Using the model from Python

```python
import keylink_core as core
cfg     = core.KeylinkConfig.from_files('.')       # all KL_*.txt parameter files of the directory
model   = core.KeylinkModel(cfg)
results = model.run(cfg.gmax)                       # any gmax vector; can be called repeatedly
results.pools[:, core.SOM]                          # daily SOM, g C m-3  (BACT .. PRED, LITTER, SOM, ROOTS, CO2)
results.soil_water, results.pores, results.SOMCN    # hydrology and nitrogen diagnostics
scenario = cfg.replace(pH=5.9, tStop=365)           # variants: validated copies of the configuration
```

`KeylinkConfig` is an immutable, validated dataclass (one field per parameter, documented in
`keylink_core.py`); `KeylinkModel.run()` allocates its own output arrays, so runs never share or
mutate state, and importing `keylink_core` has no side effects. `keylink_main.py` shows the three
run modes (single, bayesian, posterior) on top of this API.

Code layout (dependencies point downwards only; checked by `tests/test_modules.py`):

| module | concern |
|---|---|
| `keylink_core.py` | `KeylinkConfig`, `KeylinkModel` (daily loop, N budget), `KeylinkResults`, plotting |
| `keylink_foodweb.py` | growth modifiers, growth / respiration / turnover, `fCompSpecies` (ODE right-hand side) |
| `keylink_structure.py` | pore size distribution, availability of pools to consumers, bioturbation |
| `keylink_hydrology.py` | water balance, evapotranspiration, van Genuchten functions, DOM adsorption |
| `keylink_plant.py` | external plant inputs (litter, roots, exudation, C to mycorrhiza) and plant N uptake |
| `keylink_params.py`, `keylink_names.py` | parameter-file schema / reader, named indices |
| `keylink_functions.py` | compatibility facade re-exporting the four equation modules |

## Documentation

* `docs/equations.md` maps every equation of Flores et al. (2021) (Eqs. 1-37) and of the
  calculation order to the function that implements it, and lists the deviations from the papers
  (and the remaining known defects) in one place.
* Every module, class and function has a docstring stating purpose, arguments and **units**;
  `tests/test_modules.py` checks that this stays true and that the equation map only refers to
  functions that exist.
* `archive/` contains superseded scripts and the functions the model no longer uses; see
  `archive/README.md`.

## Parameter files

The `KL_*.txt` files keep their format (numbers followed by a `# comment`, one line per parameter
or one column per functional group), but they are read **by name**: `keylink_params.SCHEMA` lists
for every file the expected lines in order, a keyword that must occur in the line's comment, the
unit and the valid range. A missing, reordered, mislabelled or out-of-range line raises a
`ParameterFileError` naming the file and line, instead of silently shifting the values below it.
The group header of `KL_FaunalParams.txt` / `KL_densityParams.txt` must list the groups in the
order bact, fungi, myc (EM), bacterivores, fungivores, detritivores, engineers, herbivores,
predators. To see what a file must contain:

    python -c "import keylink_params as kp; print(kp.describe('KL_initSoil'))"

Pool, flux and pore-class indices are named in `keylink_names.py` (`Pool.SOM`, `Flux.SOM_EATEN`,
`MESO`, ...); `results.pool('SOM')`, `results.flux('SOMeaten')` and `results.as_dict()` give
named access to the daily output.

## Density dependent turnover (optional extension)

`KeylinkConfig.density_dependent_turnover` (`from_files(..., density_dependent_turnover=True)`, default on). The turnover (death)
flux of each functional group is `d * B * (B/Bref)**(beta-1)` (Georgiou et al. 2017), with
`beta` and `Bref` per group in `KL_densityParams.txt` (`Bref <= 0` uses the initial biomass;
`beta = 1` or the flag set to `False` gives the first-order turnover `d * B` of Flores et al.
2021, Eq. 12). This self-limitation removes the boom-bust behaviour of the exponential growth
equation (microbial peaks far above plausible values followed by collapse) while leaving the
turnover at the observed biomass unchanged.


## Installation and tests

Dependencies: `pip install -r requirements.txt` (numpy, scipy, matplotlib; pytest is optional).

To check that everything still works after a change, run from the model directory

    python run_tests.py          # add -v to list every test, -k NAME to select tests

(equivalent: `python -m unittest discover -s tests -t .` or `pytest tests`). The suite
(`tests/`) covers the water modifier, the daily derivative functions, the carbon and nitrogen
budgets and the density-dependent turnover; it runs on a two-year climate fixture
(`tests/fixtures/climate_1999_2000.txt`) and takes well under a minute. Shared set-up lives in
`tests/helpers.py` (`make_config(ndays, ...)`, `run_model(...)`).
