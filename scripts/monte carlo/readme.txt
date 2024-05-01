Scripts and data in this folder are mirrored from aq_mp_tox_shiny. Ensure the following files are up to date with that repo:
 - data/input/aoc_setup.RDS (this is ToMEx 1.0 data - copied from the main folder in aq_mp_tox_shiny; generated with RDAMaker there.)
 - data/input/aoc_z_tomex2.rds (tomex 2.0 data generated from Master_data_tidying.R in aq_mp_tox_shiny)
 - scripts/monte carlo/ref data/tomex2_input.rds (tomex 2.0 BASE file generated from collated excel scripts in aq_mp_tox_shiny master_data_tidying_script)
- scripts/functions.R (copied from aq_mp_tox_shiny main folder) - contains all particle property estimation functions
- scrips/monte carlo/ref data/gape_size.csv - (copied from aq_mp_tox_shiny main folder)