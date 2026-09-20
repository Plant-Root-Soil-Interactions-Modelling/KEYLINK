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

Tests: `python -m unittest discover -p "test_*.py"`
