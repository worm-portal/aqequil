"""
3o_mine.py - Python conversion of 3o_mine.r

Functions for mining EQ3 output (.3o) files.
"""

import pandas as pd
import numpy as np
import re
import os
import time


### Helper functions

def trimspace(string):
    """Trims away leading and trailing spaces and condenses multiple spaces between words."""
    return re.sub(r'(?<=\s)\s+|^\s+|\s+$', '', string)


def isolate_block(string, begin_str, end_str):
    """Isolate a substring by trimming off the portions before and after it."""
    result = re.sub(begin_str, '', string, flags=re.DOTALL)
    result = re.sub(end_str, '', result, flags=re.DOTALL)
    return result


### Main functions

def mine_3o(this_file,
            this_pressure,
            get_aq_dist=True,
            get_mass_contribution=True,
            get_mineral_sat=True,
            get_redox=True,
            get_charge_balance=True,
            get_ion_activity_ratios=True,
            get_fugacity=True,
            get_basis_totals=True,
            get_solid_solutions=True,
            mass_contribution_other=True,
            verbose=1):
    """
    Mine data from a .3o EQ3 output file.

    Parameters
    ----------
    this_file : str
        Filename of .3o file
    this_pressure : float
        Pressure value
    get_aq_dist : bool
        Extract aqueous distribution data
    get_mass_contribution : bool
        Extract mass contribution data
    get_mineral_sat : bool
        Extract mineral saturation data
    get_redox : bool
        Extract redox data
    get_charge_balance : bool
        Extract charge balance data
    get_ion_activity_ratios : bool
        Extract ion activity ratios
    get_fugacity : bool
        Extract fugacity data
    get_basis_totals : bool
        Extract basis totals
    get_solid_solutions : bool
        Extract solid solution data
    mass_contribution_other : bool
        Include "Other" category in mass contribution
    verbose : int
        Verbosity level

    Returns
    -------
    dict
        Dictionary containing mined data
    """

    # Set directory to rxn_3o folder where .3o files are kept
    os.chdir("rxn_3o")

    # Read .3o file as a string
    with open(this_file, 'r') as f:
        extractme = f.read()

    # Get sample name
    this_name = trimspace(isolate_block(extractme, begin_str=r'^.*\|Sample:\s+', end_str=r'\|\n\|.*$'))

    if verbose > 1:
        print(f"Processing EQ3 output for {this_name}")

    # Check if file experienced errors. If so, skip processing the file:
    if "Normal exit" not in extractme or "* Error" in extractme:
        os.chdir("../")
        return {}

    sample_3o = {}

    sample_3o["filename"] = this_file
    sample_3o["name"] = this_name

    ### Begin mining temperature, pressure, water properties

    # Mine params
    sample_3o["temperature"] = isolate_block(extractme, begin_str=r'^.*Temperature=\s+', end_str=r'\s+.*$')
    sample_3o["pressure"] = this_pressure
    sample_3o["logact_H2O"] = isolate_block(extractme, begin_str=r'^.*Log activity of water=\s+', end_str=r'\s+.*$')
    sample_3o["H2O_density"] = isolate_block(extractme, begin_str=r'^.*Solution density =\s+', end_str=r'\s+.*$')
    sample_3o["H2O_molality"] = 55.348 / float(sample_3o["H2O_density"])
    sample_3o["H2O_log_molality"] = np.log10(sample_3o["H2O_molality"])
    sample_3o["ionic_strength"] = isolate_block(extractme, begin_str=r'^.*Ionic strength \(I\)=\s+', end_str=r'\s+.*$')

    ### Begin extracting 'Distribution of Aqueous Solute Species'
    if get_aq_dist:
        # String to isolate the aqueous species distribution section:
        front_trim = r"^.*\n\n\n\n                --- Distribution of Aqueous Solute Species ---\n\n    Species                  Molality    Log Molality   Log Gamma  Log Activity\n\n\s+"

        # Isolate species distribution block
        species_block = isolate_block(extractme, begin_str=front_trim, end_str=r"\n\n.*$")

        # Split into substrings, each representing a separate row in the table
        species_block = species_block.split("\n")

        # Create an empty data frame to store results
        df_data = []

        # Convert into dataframe
        for this_row in species_block:
            # Mine row data
            this_row = trimspace(this_row)
            this_row_data = this_row.split()

            if len(this_row_data) >= 5:
                df_data.append({
                    'species': this_row_data[0],
                    'molality': this_row_data[1],
                    'log_molality': this_row_data[2],
                    'log_gamma': this_row_data[3],
                    'log_activity': this_row_data[4]
                })

        if len(df_data) > 0:
            df = pd.DataFrame(df_data)

            if "H2O" not in df['species'].values:
                # Add a row for water
                df = pd.concat([df, pd.DataFrame([{
                    'species': 'H2O',
                    'molality': sample_3o["H2O_molality"],
                    'log_molality': sample_3o["H2O_log_molality"],
                    'log_gamma': 1,
                    'log_activity': sample_3o["logact_H2O"]
                }])], ignore_index=True)

            # Set index as species names
            df = df.set_index('species')
        else:
            # Create DataFrame with just water if no species found
            df = pd.DataFrame([{
                'species': 'H2O',
                'molality': sample_3o["H2O_molality"],
                'log_molality': sample_3o["H2O_log_molality"],
                'log_gamma': 1,
                'log_activity': sample_3o["logact_H2O"]
            }])
            df = df.set_index('species')

        # Add aqueous block to this sample data
        sample_3o["aq_distribution"] = df

    # End of 'aqueous distribution' extraction

    if get_mass_contribution:
        ### Begin extracting 'Major Species by Contribution to Aqueous Mass Balances'

        # String to isolate the species saturation section:
        front_trim = r"^.*\n\n\n      --- Major Species by Contribution to Aqueous Mass Balances ---\n\n\n"

        # Isolate contribution block
        contrib_block = isolate_block(extractme, begin_str=front_trim, end_str=r"\n\n\n\n.*$")

        # Split into substrings, each representing a separate row in the table
        contrib_block = contrib_block.split("\n")
        # Remove blank lines
        contrib_block = [line for line in contrib_block if line != ""]

        # Loop through rows in this block and mine contributions
        mine_vals = False
        mass_contribution = {}
        for this_row in contrib_block:
            if "Accounting for" in this_row:
                # Get basis species for this block
                this_basis = this_row.replace(" Species Accounting for 99% or More of Aqueous ", "")
            elif "Per Cent" in this_row:
                # Get ready to mine data for this basis species
                mine_vals = True
                df_basis_data = []
            elif mine_vals and " - - - - - - - - -" not in this_row:
                # Mine data from this row
                row_data = trimspace(this_row)
                row_data = row_data.split()
                if len(row_data) >= 4:
                    df_basis_data.append({
                        'species': row_data[0],
                        'factor': row_data[1],
                        'molality': row_data[2],
                        'percent': row_data[3]
                    })
            elif " - - - - - - - - -" in this_row:
                # Stop mining for this basis species
                mine_vals = False
                df_basis = pd.DataFrame(df_basis_data)
                # Specify index for this contribution block
                df_basis = df_basis.set_index('species')
                # Add contribution data to list of sample data
                mass_contribution[this_basis] = df_basis

        sample_3o["mass_contribution"] = mass_contribution

    # End 'aqueous contribution' extraction

    ### Begin mining mineral saturation section
    if get_mineral_sat:

        # String to isolate the mineral saturation section:
        front_trim = r"^.*\n\n\n\n           --- Saturation States of Pure Solids ---\n\n       Phase                      Log Q/K    Affinity, kcal\n\n\s+"

        # Isolate mineral block
        mineral_block = isolate_block(extractme, begin_str=front_trim, end_str=r"\n\n.*$")

        # Split into substrings, each representing a separate row in the table
        mineral_block = mineral_block.split("\n")

        # Create an empty data frame to store results
        df_data = []

        # Convert into dataframe
        for this_row in mineral_block:
            # Get row data
            this_row_data = trimspace(this_row).split()

            if len(this_row_data) >= 3:
                df_data.append({
                    'mineral': this_row_data[0],
                    'logQoverK': this_row_data[1],
                    'affinity': this_row_data[2]
                })

        if len(df_data) > 0:
            df = pd.DataFrame(df_data)
            df = df.set_index('mineral')
        else:
            # Create empty DataFrame with correct structure
            df = pd.DataFrame(columns=['logQoverK', 'affinity'])
            df.index.name = 'mineral'

        # Add mineral saturation block to this sample data
        sample_3o["mineral_sat"] = df

        if get_solid_solutions:
            if "--- Saturation States of Hypothetical Solid Solutions ---" in extractme:
                # String to isolate the solid solution saturation section:
                front_trim = r"^.*\n\n\n                --- Saturation States of Hypothetical Solid Solutions ---\n\n"

                # Isolate solid solution block
                ss_block = isolate_block(extractme, begin_str=front_trim, end_str=r"\n\n                     --- Fugacities ---.*$")

                if ss_block != " None":

                    # Split into substrings, each representing a separate solid solution
                    ss_block = ss_block.split("\n\n\n                --- ")

                    ss_entries = {}
                    for ss_entry in ss_block:
                        ss_entry_split = ss_entry.split(" ---\n\n   ")
                        ss_name = ss_entry_split[0]
                        ss_name = ss_name.replace("\n                --- ", "")  # Clean up first entry name

                        if len(ss_entry_split) > 1:
                            ss_data = ss_entry_split[1]
                            ss_data = ss_data.replace("Ideal solution\n\n    Component                    x           Log x   Log lambda  Log activity\n\n", "")

                            ss_split = ss_data.split("\n\n\n    Mineral                       Log Q/K         Aff, kcal    State\n\n")

                            ss_dict = {}

                            # Process ideal solution data
                            if len(ss_split) > 0:
                                ideal_lines = ss_split[0].strip().split("\n")
                                ideal_data = []
                                for line in ideal_lines:
                                    line = line.lstrip()
                                    parts = re.split(r'\s{2,}', line)
                                    if len(parts) >= 5:
                                        # Skip header lines by checking if we can convert to float
                                        try:
                                            x_val = parts[1]
                                            # Try to convert x_val to float to check if it's a data line
                                            try:
                                                x_val = float(x_val)
                                            except:
                                                x_val = 0
                                            # Verify parts[2] is also numeric (not "Log x" header)
                                            log_x = float(parts[2])
                                            ideal_data.append({
                                                'component': parts[0],
                                                'x': x_val,
                                                'Log x': log_x,
                                                'Log lambda': float(parts[3]),
                                                'Log activity': float(parts[4])
                                            })
                                        except (ValueError, IndexError):
                                            # Skip header or invalid lines
                                            continue
                                ss_dict["ideal solution"] = pd.DataFrame(ideal_data)

                            # Process mineral data
                            if len(ss_split) > 1:
                                mineral_lines = ss_split[1].strip().split("\n")
                                mineral_data = []
                                for line in mineral_lines:
                                    line = line.lstrip()
                                    parts = re.split(r'\s{2,}', line)
                                    if len(parts) >= 3:
                                        try:
                                            # Try to convert to float to verify it's a data line, not header
                                            log_qk = float(parts[1])
                                            aff_kcal = float(parts[2])
                                            state = parts[3] if len(parts) >= 4 else ""
                                            mineral_data.append({
                                                'mineral': parts[0],
                                                'Log Q/K': log_qk,
                                                'Aff, kcal': aff_kcal,
                                                'State': state
                                            })
                                        except (ValueError, IndexError):
                                            # Skip header or invalid lines
                                            continue
                                ss_dict["mineral"] = pd.DataFrame(mineral_data)

                            ss_entries[ss_name] = ss_dict

                    sample_3o["solid_solutions"] = ss_entries
                else:
                    sample_3o["solid_solutions"] = None
            else:
                sample_3o["solid_solutions"] = None

    # End 'mineral saturation affinity' extraction

    ### Begin mining redox data
    if get_redox:
        # String to isolate the redox section:
        front_trim = r"^.*\n\n\n\n                --- Aqueous Redox Reactions ---\n\n   Couple                           Eh, volts      pe-      log fO2   Ah, kcal\n\n\s+"

        # Isolate redox block
        redox_block = isolate_block(extractme, begin_str=front_trim, end_str=r"\n\n.*$")

        # Split into substrings, each representing a separate row in the table
        redox_block = redox_block.split("\n")

        # Create an empty data frame to store results
        df_data = []

        # Convert into dataframe
        for this_row in redox_block:

            # Get row data
            this_row_data = trimspace(this_row).split()

            if len(this_row_data) >= 5:
                df_data.append({
                    'couple': this_row_data[0],
                    'Eh': this_row_data[1],
                    'pe': this_row_data[2],
                    'logfO2': this_row_data[3],
                    'Ah': this_row_data[4]
                })

        if len(df_data) > 0:
            df = pd.DataFrame(df_data)
            df = df.set_index('couple')
        else:
            # Create empty DataFrame with correct structure
            df = pd.DataFrame(columns=['Eh', 'pe', 'logfO2', 'Ah'])
            df.index.name = 'couple'

        # Add redox block to this sample data
        sample_3o["redox"] = df

    # End redox extraction

    ### Begin mining charge balance data
    if get_charge_balance:
        # String to isolate ionic strength:
        front_trim = r"^.*Ionic strength \(I\)=\s+"

        # Isolate ionic strength
        IS = isolate_block(extractme, begin_str=front_trim, end_str=r"\s+.*$")

        # String to isolate stoichiometric ionic strength:
        front_trim = r"^.*Stoichiometric ionic strength=\s+"

        IS_stoich = isolate_block(extractme, begin_str=front_trim, end_str=r"\s+.*$")

        # String to isolate the electrical balance section:
        front_trim = r"^.*Sigma\(mz\) cations=\s+"

        elec_block = isolate_block(extractme, begin_str=front_trim, end_str=r"\n\n.*$")

        # Split electrical block into strings and numerics
        elec_parts = re.split(r'=\s+|\n\s+', elec_block)

        elec_dict = {
            "sigma(mz) cations": elec_parts[0],
            "sigma(mz) anions": elec_parts[2],
            "total charge": elec_parts[4],
            "mean charge": elec_parts[6],
            "charge imbalance": elec_parts[8]
        }

        # String to isolate charge balance:
        front_trim = r"^.*The electrical imbalance is:\n\n\s+"

        cbal_bal = isolate_block(extractme, begin_str=front_trim, end_str=r"\n\n.*$")

        # Split electrical block into strings and numerics
        cbal_parts = re.split(r' per cent|\n\s+', cbal_bal)

        cbal_dict = {
            "charge imbalance % of total charge": cbal_parts[0],
            "charge imbalance % of mean charge": cbal_parts[2]
        }

        charge_balance_dict = {
            "ionic strength": IS,
            "stoichiometric ionic strength": IS_stoich,
            **elec_dict,
            **cbal_dict
        }

        sample_3o["charge_balance"] = charge_balance_dict

    # End charge balance extraction

    if get_ion_activity_ratios:
        ion_ratio_block = isolate_block(extractme, r"^.*--- Ion-H\+ Activity Ratios ---\n\n", r"\n\n.*$")
        ion_ratio_block_split = ion_ratio_block.split("\n")
        ion_ratio_block_split = [line.split("=") for line in ion_ratio_block_split]

        if len(ion_ratio_block_split) > 0 and len(ion_ratio_block_split[0]) > 1:

            ion_ratio_logs = [trimspace(item[0]) for item in ion_ratio_block_split if len(item) > 1]
            ion_ratio_values = []
            for item in ion_ratio_block_split:
                if len(item) > 1:
                    try:
                        ion_ratio_values.append(float(item[1]))
                    except:
                        ion_ratio_values.append(np.nan)

            which_to_divide = ["/" in log for log in ion_ratio_logs]

            # Extract hydrogen exponents and ions
            hydrogen_exponents = []
            ions = []
            for log in ion_ratio_logs:
                cleaned = log.replace("Log ( a(", "").replace(" )", "")
                parts = re.split(r'\)xx ', cleaned)
                if len(parts) > 1:
                    hydrogen_exponents.append(float(parts[1]))
                else:
                    hydrogen_exponents.append(np.nan)

                ion_parts = re.split(r'\) [x|/] a\(', cleaned)
                if len(ion_parts) > 0:
                    ions.append(ion_parts[0])
                else:
                    ions.append("")

            df_data = []
            for i in range(len(ion_ratio_values)):
                df_data.append({
                    'values': ion_ratio_values[i],
                    'H_exponent': hydrogen_exponents[i],
                    'divide': which_to_divide[i],
                    'ion': ions[i]
                })

            df = pd.DataFrame(df_data)
            sample_3o["ion_activity_ratios"] = df

    ### Begin fugacity mining
    if get_fugacity:
        fugacity_block = isolate_block(extractme, r"^.*--- Fugacities ---\n\n", r"\n\n\n.*$")
        lines = fugacity_block.split("\n")
        lines = [re.sub(r'\s+', ' ', line.strip()) for line in lines]
        lines = lines[2:]  # Skip header lines

        df_data = []
        for line in lines:
            parts = line.split()
            if len(parts) >= 2:
                df_data.append({
                    'gas': parts[0],
                    'log_fugacity': float(parts[1])
                })

        if len(df_data) > 0:
            df = pd.DataFrame(df_data)
            df = df.set_index('gas')
        else:
            # Create empty DataFrame with correct structure
            df = pd.DataFrame(columns=['log_fugacity'])
            df.index.name = 'gas'

        sample_3o["fugacity"] = df

    ### Begin sensible composition mining ("basis totals")
    if get_basis_totals:
        sc_block = isolate_block(extractme, r"^.*--- Sensible Composition of the Aqueous Solution ---\n\n", r"\n\n   The above data have.*$")
        lines = sc_block.split("\n")
        lines = [re.sub(r'\s+', ' ', line.strip()) for line in lines]
        lines = lines[2:]  # Skip header lines

        df_data = []
        for line in lines:
            parts = line.split()
            if len(parts) >= 5:
                df_data.append({
                    'species': parts[0] + "_total",
                    'mg/L': float(parts[1]),
                    'mg/kg.sol': float(parts[2]),
                    'molarity': float(parts[3]),
                    'molality': float(parts[4])
                })

        if len(df_data) > 0:
            df = pd.DataFrame(df_data)
            df = df.set_index('species')
        else:
            # Create empty DataFrame with correct structure
            df = pd.DataFrame(columns=['mg/L', 'mg/kg.sol', 'molarity', 'molality'])
            df.index.name = 'species'

        sample_3o["basis_totals"] = df

    os.chdir("../")

    return sample_3o


def melt_mass_contribution(batch_3o, other=False, verbose=1):
    """
    Melt aqueous contribution data from multiple samples into a single dataframe.

    Parameters
    ----------
    batch_3o : dict
        Batch data from multiple .3o files
    other : bool
        Include "Other" category
    verbose : int
        Verbosity level

    Returns
    -------
    DataFrame
        Melted mass contribution data
    """

    # Initialize empty list for data
    df_aq_cont_data = []

    # Get all aqueous contribution data
    mass_contributions = {sample: data.get('mass_contribution', {})
                         for sample, data in batch_3o["sample_data"].items()}

    # Loop through each sample and basis species
    for sample in mass_contributions:
        if verbose > 1:
            print(f"Processing mass contribution of basis species in {sample}...")
        for basis in mass_contributions[sample]:
            df = mass_contributions[sample][basis].copy()
            df['basis'] = basis
            df['sample'] = sample
            df['species'] = df.index
            df = df.reset_index(drop=True)

            if other:
                percent = round(100 - sum(pd.to_numeric(df['percent'], errors='coerce')), 2)
                other_row = pd.DataFrame([{
                    'sample': sample,
                    'basis': basis,
                    'species': 'Other',
                    'factor': np.nan,
                    'molality': np.nan,
                    'percent': str(percent)
                }])
                df = pd.concat([df, other_row], ignore_index=True)

            df_aq_cont_data.append(df)

    if len(df_aq_cont_data) > 0:
        df_aq_cont = pd.concat(df_aq_cont_data, ignore_index=True)
        df_aq_cont = df_aq_cont[['sample', 'basis', 'species', 'factor', 'molality', 'percent']]
    else:
        df_aq_cont = pd.DataFrame(columns=['sample', 'basis', 'species', 'factor', 'molality', 'percent'])

    return df_aq_cont


def create_report_df(data, category, out_type):
    """
    Create report versions of data categories.

    Parameters
    ----------
    data : dict
        Sample data dictionary
    category : str
        Category name
    out_type : int or str
        Output type (column index or name)

    Returns
    -------
    DataFrame
        Report dataframe
    """

    df_cat = {sample: data[sample].get(category) for sample in data if category in data[sample]}

    # Get all unique species/columns
    all_species = set()
    for sample_df in df_cat.values():
        if sample_df is not None:
            if isinstance(sample_df, pd.DataFrame):
                # Check if this is ion_activity_ratios with an 'ion' column
                if 'ion' in sample_df.columns:
                    all_species.update(sample_df['ion'].tolist())
                else:
                    all_species.update(sample_df.index.tolist())
            elif isinstance(sample_df, dict):
                all_species.update(sample_df.keys())

    all_species = sorted(list(all_species))

    # Create result dataframe
    result_data = []
    for sample in df_cat:
        row_data = {'sample': sample}
        sample_df = df_cat[sample]

        if sample_df is not None:
            if isinstance(sample_df, pd.DataFrame):
                # Check if this DataFrame has an 'ion' column (like ion_activity_ratios)
                if 'ion' in sample_df.columns:
                    # Use the 'ion' column to identify rows
                    for species in all_species:
                        matching_rows = sample_df[sample_df['ion'] == species]
                        if len(matching_rows) > 0:
                            if isinstance(out_type, int):
                                row_data[species] = matching_rows.iloc[0].iloc[out_type]
                            else:
                                row_data[species] = matching_rows.iloc[0][out_type] if out_type in sample_df.columns else np.nan
                        else:
                            row_data[species] = np.nan
                else:
                    # Use the index like before
                    for species in all_species:
                        if species in sample_df.index:
                            if isinstance(out_type, int):
                                row_data[species] = sample_df.loc[species].iloc[out_type] if len(sample_df.loc[species].shape) > 0 else sample_df.loc[species]
                            else:
                                row_data[species] = sample_df.loc[species, out_type] if out_type in sample_df.columns else np.nan
                        else:
                            row_data[species] = np.nan
            elif isinstance(sample_df, dict):
                for species in all_species:
                    row_data[species] = sample_df.get(species, np.nan)

        result_data.append(row_data)

    df_result = pd.DataFrame(result_data)

    # Sort columns alphabetically (except 'sample')
    cols = ['sample'] + sorted([col for col in df_result.columns if col != 'sample'])
    df_result = df_result[cols]

    return df_result


def compile_report(data, csv_filename, aq_dist_type, mineral_sat_type,
                   redox_type, get_aq_dist, get_mineral_sat, get_redox,
                   get_charge_balance, get_ion_activity_ratios, get_fugacity,
                   get_basis_totals, input_processed_df,
                   df_input_processed_names):
    """
    Compile a report from mined .3o data.

    Returns
    -------
    dict
        Report dictionary with 'report' and 'divs' keys
    """

    report_list = {}
    report_list["divs"] = {}

    # Initialize report with processed input file
    report = input_processed_df.copy()
    report.columns = df_input_processed_names

    report_list["divs"]["input"] = list(report.columns[1:])  # Exclude "Sample" column

    # Create report versions of EQ3 output blocks
    # Use suffixes=('','') to prevent pandas from adding _x, _y suffixes
    # Instead, we'll handle duplicates by keeping only the first occurrence

    if get_aq_dist:
        aq_distribution = create_report_df(data=data, category='aq_distribution', out_type=aq_dist_type)
        report = report.merge(aq_distribution, left_on='Sample', right_on='sample', how='inner', suffixes=('', ''))
        report = report.drop('sample', axis=1)
        report_list["divs"]["aq_distribution"] = list(aq_distribution.columns[1:])

    if get_mineral_sat:
        mineral_sat = create_report_df(data=data, category='mineral_sat', out_type=mineral_sat_type)
        report = report.merge(mineral_sat, left_on='Sample', right_on='sample', how='inner', suffixes=('', ''))
        report = report.drop('sample', axis=1)
        report_list["divs"]["mineral_sat"] = list(mineral_sat.columns[1:])

    if get_redox:
        redox = create_report_df(data=data, category='redox', out_type=redox_type)
        report = report.merge(redox, left_on='Sample', right_on='sample', how='inner', suffixes=('', ''))
        report = report.drop('sample', axis=1)
        report_list["divs"]["redox"] = list(redox.columns[1:])

    if get_charge_balance:
        charge_balance = create_report_df(data=data, category='charge_balance', out_type=0)
        report = report.merge(charge_balance, left_on='Sample', right_on='sample', how='inner', suffixes=('', ''))
        report = report.drop('sample', axis=1)
        report_list["divs"]["charge_balance"] = list(charge_balance.columns[1:])

    if get_ion_activity_ratios:
        if any('ion_activity_ratios' in data[sample] for sample in data):
            ion_activity_ratios = create_report_df(data=data, category='ion_activity_ratios', out_type=0)
            # Rename columns to add the suffix before merging to avoid conflicts
            rename_dict = {col: col + '_Log_ion-H+_activity_ratio'
                          for col in ion_activity_ratios.columns if col != 'sample'}
            ion_activity_ratios = ion_activity_ratios.rename(columns=rename_dict)
            report = report.merge(ion_activity_ratios, left_on='Sample', right_on='sample', how='inner', suffixes=('', ''))
            report = report.drop('sample', axis=1)
            report_list["divs"]["ion_activity_ratios"] = [col + '_Log_ion-H+_activity_ratio'
                                                           for col in list(create_report_df(data=data, category='ion_activity_ratios', out_type=0).columns[1:])]

    if get_fugacity:
        fugacities = create_report_df(data=data, category='fugacity', out_type=0)
        report = report.merge(fugacities, left_on='Sample', right_on='sample', how='inner', suffixes=('', ''))
        report = report.drop('sample', axis=1)
        report_list["divs"]["fugacity"] = list(fugacities.columns[1:])

    if get_basis_totals:
        sc = create_report_df(data=data, category='basis_totals', out_type=3)  # 3 is the molality column
        report = report.merge(sc, left_on='Sample', right_on='sample', how='inner', suffixes=('', ''))
        report = report.drop('sample', axis=1)
        report_list["divs"]["basis_totals"] = list(sc.columns[1:])

    report = report.set_index('Sample')

    report_list["report"] = report

    return report_list


### Main
def main_3o_mine(files_3o,
                 get_aq_dist,
                 get_mass_contribution,
                 get_mineral_sat,
                 get_redox,
                 get_charge_balance,
                 get_ion_activity_ratios,
                 get_fugacity,
                 get_basis_totals,
                 get_solid_solutions,
                 mass_contribution_other,
                 csv_filename,
                 aq_dist_type,
                 mineral_sat_type,
                 redox_type,
                 input_filename,
                 input_pressures,
                 batch_3o_filename,
                 df_input_processed,
                 df_input_processed_names,
                 verbose):
    """
    Main function to mine multiple .3o files.

    Returns
    -------
    dict
        Batch data from all .3o files
    """

    start_time = time.time()

    # Instantiate an empty object to store data from all 3o files
    batch_3o = {}
    batch_3o["sample_data"] = {}

    if verbose > 1:
        print("Now processing EQ3 output files...")

    # Create dict mapping files to pressures
    pressure_dict = dict(zip(files_3o, input_pressures))

    # Process each .3o file
    for file in files_3o:

        # Add this sample's aqueous data to list of all sample data
        sample_3o = mine_3o(file,
                           this_pressure=pressure_dict[file],
                           get_aq_dist=get_aq_dist,
                           get_mass_contribution=get_mass_contribution,
                           get_mineral_sat=get_mineral_sat,
                           get_redox=get_redox,
                           get_charge_balance=get_charge_balance,
                           get_ion_activity_ratios=get_ion_activity_ratios,
                           get_fugacity=get_fugacity,
                           get_basis_totals=get_basis_totals,
                           get_solid_solutions=get_solid_solutions,
                           mass_contribution_other=mass_contribution_other,
                           verbose=verbose)

        # If this file could be processed, add its data to the batch_3o object
        if len(sample_3o) > 1:
            batch_3o["sample_data"][sample_3o["name"]] = sample_3o

    if verbose > 1:
        print("Finished processing EQ3 output files...")

    # Compile aqueous contribution data into a single melted dataframe and
    # append it to the batch_3o object.
    if get_mass_contribution and len(batch_3o["sample_data"]) > 0:
        if verbose > 1:
            print("Now processing mass contribution data...")
        batch_3o["mass_contribution"] = melt_mass_contribution(batch_3o=batch_3o,
                                                               other=mass_contribution_other,
                                                               verbose=verbose)
        if verbose > 1:
            print("Finished processing mass contribution data...")

    if len(batch_3o["sample_data"]) > 0:
        # Create a report summarizing 3o data from all samples
        report_list = compile_report(data=batch_3o["sample_data"],
                                     csv_filename=csv_filename,
                                     aq_dist_type=aq_dist_type,
                                     mineral_sat_type=mineral_sat_type,
                                     redox_type=redox_type,
                                     get_aq_dist=get_aq_dist,
                                     get_mineral_sat=get_mineral_sat,
                                     get_redox=get_redox,
                                     get_charge_balance=get_charge_balance,
                                     get_ion_activity_ratios=get_ion_activity_ratios,
                                     get_fugacity=get_fugacity,
                                     get_basis_totals=get_basis_totals,
                                     input_processed_df=df_input_processed,
                                     df_input_processed_names=df_input_processed_names)

        # Add the report to the batch_3o object
        report = report_list["report"]
        batch_3o["report"] = report
        batch_3o["report_divs"] = report_list["divs"]
    else:
        return {}

    # Store user input file data
    batch_3o["input"] = pd.read_csv(input_filename)

    # Save the batch_3o object (would use pickle in Python)
    if batch_3o_filename is not None:
        import pickle
        with open(batch_3o_filename, 'wb') as f:
            pickle.dump(batch_3o, f)

    time_elapsed = time.time() - start_time
    if verbose > 1:
        print(f"Finished mining .3o files. Time elapsed: {round(time_elapsed, 2)} seconds")

    return batch_3o
