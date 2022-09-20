### Nonlinear CBC sims
## General notes
- explicit_option = "rk2"
- nx = ny = 91 ; y0 = 15
- nzed = 64
- No dissipation
- Electrostatic

## List of simulations
# adiabatic sims:
- Created in earlier benchmarking experiments; the summary_pickle files don't exist
# stella_nonlinear_2species_isl
- Uses our new splitting scheme, but treats the NL term using ISL
# stella_nonlinear_2species_leapfrog_nonlinear
- Uses our new splitting scheme, treating the NL term using a Leapfrog-like approach
# stella_nonlinear_2species_master
- Uses master
# stella_nonlinear_2species_master_archer2
- Uses master (longer sim time than stella_nonlinear_2species_master)
# stella_nonlinear_2species_nisl
- Actually a "broken" simulation, doesn't use NISL
# stella_nonlinear_2species_nisl_archer2
- NISL, run on archer2
# stella_nonlinear_2species_nisl_delt_004
- NISL, run on archer2, with smaller dt, but seems to only have run for a single timestep

## Of interest to the thesis
# stella_nonlinear_2species_master_archer2
- Example of a saturating NL sim
# stella_nonlinear_2species_nisl_archer2
- Example of a NISL sim
# stella_nonlinear_2species_leapfrog_nonlinear
- Example of the leapfrog splitting
# stella_nonlinear_2species_isl
- Example of ISL
