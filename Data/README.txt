DATA - Group 4A (CL 202, IIT Gandhinagar)
Rohan Jaybhave, Janil Jain, Preet Patel, Saksham Malav, Manish Kumawat

All files are CSV or TXT. They are produced by the code in Codes/python:
    cd Codes/python
    python run_all.py        (runs export_data.py at the end and rewrites this folder)
Units: temperature in deg C unless the column says K; solubility in mg of compound per
100 g of oil; masses in mg or g as named; x = mole fraction; gamma = activity coefficient.
"inf" means no limit (molten compound fully miscible with the oil). Oils: sunflower and
coconut (corn, olive and MCT appear only in the literature validation).

INPUTS (inputs/)
  compound_properties.csv            M, Tm, dHfus, dCp, UNIFAC groups and sources (vanillin, curcumin, piperine)
  oil_fatty_acid_compositions.csv    fatty-acid wt% of each oil, pseudo-triglyceride molar mass and groups
  unifac_subgroups.csv               UNIFAC R, Q and main group of every subgroup used
  unifac_interaction_parameters.csv  UNIFAC a_nm, a_mn (K) for every main-group pair used
  literature_solubility.csv          published solubilities used for validation (pure solvents and oils)
  flavordb_records.csv               FlavorDB2 records (IIIT-Delhi): CAS, FEMA, flavour profile, XLogP
  flavordb_spice_entities.csv        FlavorDB2 entries for the spices
  spicerx_links.csv                  SpiceRx (IIIT-Delhi) spice to compound links
  lab_measurements_TEMPLATE.csv      template for the group's own measurements
  scenarios_saturation_single.csv    60 single-compound cases: 3 compounds x 2 oils x 5 temperatures
                                     (25, 60, 100, 140, 180 C) x 2 amounts (tadka dose, 1000 mg) in 15 g oil
  scenarios_saturation_mixture.csv   24 mixture cases: 2 oils x 4 temperatures (25, 60, 100, 140 C) x
                                     3 recipes (tadka V+C+P, turmeric + pepper C+P, high loading V+C+P)
  route_inputs_R1_to_R5/*.txt        the JSON input files of the full workflow (routes R1 to R5,
                                     a fallback case, a custom oil and a deliberately wrong input)

OUTPUTS (outputs/)
  01_single_compound_solubility/
    solubility_<compound>_in_<oil>_oil_25-200C.csv
        1 C steps: status (SLE / LLE-limited / miscible), x, ideal x, gamma and its
        combinatorial and residual parts, solubility, ideal solubility, saturation point in 15 g oil
  02_saturation_point_single_compound/
    saturation_single_scenarios.csv     for each scenario: saturation point in the pan, dissolved,
                                        amount that cannot mix in, saturated (True/False),
                                        saturation temperature and the printed statement
    saturation_single_statements.txt    the statements, e.g. "At 60.0 C saturation is achieved ...
                                        More mixing cannot be done; ..."
    saturation_temperature_vs_loading.csv  temperature needed to dissolve 10 to 10000 mg per 100 g oil
  03_mixtures/
    saturation_mixture_scenarios.csv    per compound and scenario: dissolved, cannot mix in,
                                        saturation point in the mixture, gamma, saturated, temperature
                                        at which the whole recipe dissolves, statement
    saturation_mixture_statements.txt   the mixture statements
    mixture_pairs_and_triple.csv        all compounds in excess: solubility in mixture / alone
    mixture_cross_effects.csv           three compounds saturated together vs alone, 25-75 C
    recipe_extraction.csv               share of the tadka recipe dissolved, 25-180 C
    batch_example_sunflower_15g.csv     grams dissolved / left for a batch in 15 g sunflower oil
  04_validation/
    validation_report.txt               UNIFAC benchmark, pure-solvent AARD, oil comparison
    validation_pure_solvents.csv        point-by-point comparison (78 points)
    validation_pure_solvents_summary.csv
    validation_oils.csv                 published vs predicted solubility in edible oils
    validation_flavordb.csv             FlavorDB2 XLogP vs UNIFAC log Kow
  05_routes_R1_to_R5/
    output_<input>.txt                  full report of python -m flavorsol predict for every input file
  06_analysis/
    thermodynamic_consistency_checks.csv  gamma >= 1, x below ideal, S(T) rising, pure limit (all pass)
    thermo_summary.csv                  ln gamma parts at infinite dilution, dH_sol, oiling-out temperature
    sensitivity_60C_sunflower.csv       change in solubility for each perturbed input
    uncertainty_check.csv               delta method vs Monte Carlo intervals
    solubility_vs_temperature.csv       5 C sweep of all compounds and oils
    summary_headline_numbers.txt        all headline numbers of the manuscript
