# Archive — not part of the model

Kept for reference only. Nothing here is imported or run by the model or the tests; none of it
contains the corrections made to the model since the publications.

| file | what it is |
|---|---|
| `KEYLINK_model.py` | the pre-refactoring single-script version of the model (module-level state, own copy of the ODE right-hand side). Superseded by `keylink_core.py` + the equation modules. Still contains the original errors (dropped consumption terms, mycorrhizal gmax unused, fragmentation index, N bookkeeping). |
| `KEYLINK_HCsample_G.py` | script that ran a Latin hypercube / posterior sample of gmax vectors; imports a module (`functions_KEYLINK_G`) that no longer exists. Replaced by the `posterior` run mode of `keylink_main.py`. |
| `KL_FaunalParameters.txt` | older variant of `KL_FaunalParams.txt` without the gmax row; read by nothing. |
| `keylink_unused.py` | functions no longer called by the model `calcKD` (DOM adsorption, Eq. 37 — the leaching concept is not implemented) and `inputLitter`. |
