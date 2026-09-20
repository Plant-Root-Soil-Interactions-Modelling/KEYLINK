# KEYLINK
In the KEYLINK model, we integrate new knowledge on soil structure and its importance for soil organic matter (SOM) stabilization and hydrology, with the existing concepts on SOM pools, and elements from food web models, i.e. those from direct trophic interactions among soil organisms. 

## Density dependent turnover (optional extension)

Flag `DensityDependentTurnover` in `keylink_core.py` (default `True`). The turnover (death)
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
`tests/helpers.py`; `keylink_core.set_run_length(ndays, climatefile)` selects run length and
climate file without editing the parameter files.
