# KEYLINK equations and where they are implemented

Equation numbers refer to

* **II** — Flores O, Deckmyn G, Curiel Yuste J, Javaux M, Uvarov A, van der Linde S, De Vos B,
  Vereecken H, Jiménez J, Vinduskova O, Schnepf A (2021) *KEYLINK: towards a more integrative soil
  representation for inclusion in ecosystem scale models — II: model description, implementation
  and testing.* PeerJ 9:e10707 — the equations of the model.
* **I** — Deckmyn G et al. (2020) *KEYLINK … I. review and model concept.* PeerJ 8:e9750 — the
  concept; it contains no equations.

Units: C pools `g C m-3`, fluxes `g C m-3 day-1`, N `g N m-3`, pore volumes and water `l m-3`,
precipitation/infiltration/drainage `mm day-1`, rates `day-1`, temperature `degC`.

## Soil structure and hydrology

| Eq. | Description | Implementation |
|---|---|---|
| 1–2 | particle density `Ds`, total porosity `P%` from bulk density | **not implemented** — the textural pore volumes per class are given in `KL_initSoil.txt` (`cfg.PVstruct`); `cfg.BD` is read but unused |
| 3 | net precipitation `Pnet = P − E` | **deviation**: `KeylinkModel.run` infiltrates gross precipitation and removes the evapotranspiration at the end of the day (`keylink_hydrology.wl`) |
| 4–8 | infiltration, runoff, drainage, filling of the pore classes | `keylink_hydrology.calcPW`. `ImaxMat` = `cfg.Ksat` (the "measured or derived maximal infiltration rate"); `ImaxPor` = macropore volume + `drainmax`. Deviations: runoff also limited by the free pore space; drainage capped at `drainmax`; `ImaxMat` is not multiplied by `(1 − SAmacro)` |
| 14 | overall hydration modifier `mH2Otot` | `keylink_structure.calcAvail` |
| 15 | `mH2O` of the pore class that is partially filled | `keylink_structure.calcAvail` |
| 16 | `mH2O` of the pore class one above it | `keylink_structure.calcAvail` (applied for every class; a class that is exactly full hands all activity to the class above, which keeps the modifier continuous in water content) |
| 27 | burrow volume `PVB` of the engineers | `keylink_structure.calcPVD` |
| 28 | aggregated fraction `Ag` of the SOM | `keylink_structure.calcAg` (coefficient c = 10, hard-coded) |
| 29 | aggregation porosity `PVAg` | `keylink_structure.calcPVD` (inside Eqs. 33–34). **Note**: the coefficient is `k = 0.002 l g C-1`; the paper prints `2 l g C-1 m-3`, which would give >1000 l m-3 of extra porosity |
| 30 | total porosity | implicit in `calcPVD` (sum of the five classes) |
| 31–34 | pore volume per size class | `keylink_structure.calcPVD` |
| 35–36 | DOM production and adsorption | **not implemented** — SOM is a single pool without a DOM fraction |
| 37 | DOM adsorption coefficient `KD` | `archive/keylink_unused.py:calcKD` (not called) |
| — | potential evapotranspiration (Thornthwaite 1948) | `keylink_hydrology.PET`, heat index/exponent in `KeylinkModel.__init__` |
| — | removal of evapotranspiration from the pore classes | `keylink_hydrology.wl` |
| — | bioturbation: SOM and litter moved out of the layer | `keylink_structure.calcBioturb`, `calcLittermove` (paper: "Calculations regarding engineers") |

## Food web and C flow

| Eq. | Description | Implementation |
|---|---|---|
| 9 | growth `G` of a pool (Monod kinetics, summed over food sources) | `keylink_foodweb.calcgrowth`, summed per group in `fCompSpecies`. Deviation: the flux is capped at half of the available substrate per day |
| 10 | `dB/dt = G − R − Dt − Pd` | `keylink_foodweb.fCompSpecies` (one expression per functional group) |
| 11 | respiration `R = rB`, `r` temperature dependent | `keylink_foodweb.calcresp`, applied in `fCompSpecies`; cumulated in the CO2 pool |
| 12 | turnover `Dt = dB` | `keylink_foodweb.calcDeath` with `beta = 1`; the necromass is routed to SOM (microbes, microbivores) or litter (meso/macrofauna) in `fCompSpecies` |
| 13 | predation `Pd = Gpred (1 + ffaec)` | `fCompSpecies` (predator, fungivore and bacterivore terms) |
| 17 | temperature modifier `mT` | `keylink_foodweb.calcmodt` (Q10 = 2, hard-coded; respiration uses `q10` from the parameter file) |
| 18–19 | pH modifier for fungi / bacteria | `keylink_foodweb.calcgmaxmod` |
| 20 | `gmaxEng(pH)` of the engineers | `keylink_foodweb.calcgmaxEng` |
| 21–22 | C:N modifier `mC:N` for fungi / bacteria | `keylink_foodweb.calcgmaxmod` (also applied to SOM consumption, with `rec = 0`) |
| 23–24 | recalcitrance modifier `mrec` | `keylink_foodweb.calcgmaxmod` |
| 25 | `mtot = mT mC:N mrec mpH mH2O` | `calcgmaxmod` (C:N, rec, pH) × `modt` (T) in `fCompSpecies` × the availability `fa` of `calcAvail` (H2O) |
| 26 | effective faeces fraction `ffaecEff` | `keylink_foodweb.calcFaec`. Deviation: the returned value is the maximum of Eq. 26 and an undocumented N-demand term (`faecshort`), which with the standard parameters saturates the faeces fraction at 1 for litter feeders |
| — | density dependent turnover `Dt = d B (B/Bref)^(beta−1)` | `keylink_foodweb.calcDeath` (extension after Georgiou et al. 2017, suggested in II; switched by `cfg.density_dependent_turnover`) |
| — | plant inputs: litter fall, root growth/turnover, exudation, C to mycorrhiza | `keylink_plant`, applied in `KeylinkModel.run` and `fCompSpecies` |
| — | fragmentation of litter eaten by engineers into SOM | `KeylinkModel.run`. **Known defect**: it reads the spare column 21 of the state instead of `LIT_EATEN_ENG` (column 20), so the flux is always 0 |

## Nitrogen

Neither paper specifies the N bookkeeping; it is implemented in `KeylinkModel.run`:
biota have a fixed C:N (`cfg.cn`), litter C:N follows from mixing with the input, mineral N changes
by plant uptake (`keylink_plant.mycNtoPlant`, `plantNuptake`) and by the net mineralisation of the
bacteria, and SOM N is the residual of the total N budget
(`Ntot = Nbiota + Nmin + Nlitter + NSOM`). Only bacteria mineralise explicitly; the N released by
the other groups ends up in SOM through the residual. Diagnostics per day:
`results.Nmin`, `.Ntot`, `.SOMCN`, `.litterCN`, `.mineralisation`.

## Calculation order

The sequence a–j of II ("Calculation order") is the body of `KeylinkModel.run`:
a) `calcAg` + `calcPVD`, b) `calcPW`, c) DOM leaching — not implemented, d) `calcAvail`,
e) `calcgmaxmod`/`calcgmaxEng`/`calcmodt`/`calcresp`, f–g) `odeint(fCompSpecies)`,
h) litter C:N and the N budget, i) engineering effects (i1 via `calcPVD` of the next day,
i2 `calcBioturb`/`calcLittermove`), j) other changes of the pore size distribution — not implemented.

## Deviations from the papers, summarised

1. Eq. 3: evapotranspiration is subtracted after infiltration, not before.
2. Eqs. 5–8: additional storage limit on runoff, drainage capped at `drainmax`, `ImaxMat` not scaled
   by `(1 − SAmacro)`; `ImaxMat = Ksat` (constant) — using the unsaturated conductivity `K(theta)`
   instead prevents a dry soil from rewetting.
3. Eq. 9: consumption capped at half the available substrate per day.
4. Eq. 16: applied to every pore class, not only below the mesopores (continuity in water content).
5. Eq. 26: additional N-demand term in `calcFaec`.
6. Eq. 29: `k = 0.002 l g C-1` instead of the printed `2 l g C-1 m-3`.
7. Eq. 17: Q10 of growth hard-coded as 2 while respiration uses the parameter file value (2.5).
8. Eqs. 1–2, 35–36 and step (c)/(j) of the calculation order are not implemented; `calcKD` (Eq. 37)
   is archived.
9. Addition: density dependent turnover (`calcDeath`), proposed as an add-on in II.
10. Known defect kept for traceability: the fragmentation flux reads the wrong state column.
