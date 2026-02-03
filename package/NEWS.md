PSSDplusplus 0.0.2
==================

Updates
-------
- Added sediment alignment support using kg-based particle doses.
- Extended Monte Carlo outputs and thresholds to include sediment results.
- Added LHS parameter toggles and sediment-aware parameter plots.
- Fixed PSSD prep for single-species datasets and added sediment-aware PNEC labels.
- Fixed sediment alignment to respect exposure.route and derive particles/kg when only mass-based sediment doses are available.
- make_all_pSSDs now skips empty environment/ERM combinations cleanly to avoid downstream errors.
- Filled missing polydisperse length bounds when only a single length is reported to stabilize food-dilution alignment.
- PSSD generators now safely handle cases where all species lack data (returns empty matrices instead of error).
- make_all_pSSDs now reports a per-combination status summary and minimum species requirement.
- pSSD plot x-axis labels now reflect particles/L vs particles/kg based on environment.
