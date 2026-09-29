"""
Speciation groups and limiting reactant selection shared by energy supply
calculations.

Before an energy supply can be reported, the code has to decide how much of
each reactant is actually available. That decision involves speciation
groups: a microbe oxidizing sulfide can draw on every dissolved form of
sulfide (H2S, HS-, and the sulfide complexes), not just the one form that
happens to appear in the reaction the user wrote.

Both `Speciation` (AqSpeciation.py) and `Mass_Transfer` (MassTransfer.py)
need this, so the logic lives here as a mixin. Keeping it in one place is
what prevents the two classes from disagreeing about the energy supply of
the same reaction written two different ways.
"""

import collections
import copy

import pandas as pd
import pychnosz

from wormutils import chemlabel, format_equation, import_package_file, isnotebook


class SpeciationGroups(object):

    """
    Mixin providing speciation group lookups and limiting reactant selection.

    Classes mixing this in are expected to provide:
      - `self.thermo.csv_db` : WORM-style thermodynamic database dataframe
      - `self.err_handler`   : a wormutils Error_Handler
      - `self.verbose`       : verbosity level
    """

    # Class-level defaults so classes that never set these in their __init__
    # (e.g. Mass_Transfer) still work. Instance attributes shadow these.
    custom_grouping_filepath = None
    reactant_dict_scalar = None
    speciation_group_dict_unpacked = None

    # The grouping file to fall back on, set from the `custom_grouping_filepath`
    # parameter of AqEquil. None means the built-in speciation groups.
    default_grouping_filepath = None

    # {group key: (category, moiety count)}, filled in alongside the dictionaries
    _group_parts = {}

    # the grouping file the current dictionaries were built from
    _speciation_group_source = None

    def _set_speciation_groups(self, custom_grouping_filepath=None, reset=False):

        """
        Make sure the speciation group dictionaries exist and are up to date.

        They are rebuilt if they have never been built, or if the grouping file
        they were built from is not the one being asked for now. Without the
        latter check, a `custom_grouping_filepath` passed to a second energy
        calculation would be silently ignored.

        Parameters
        ----------
        custom_grouping_filepath : str, optional
            Filepath of a TXT file of customized speciation groups.

        reset : bool, default False
            If True, `custom_grouping_filepath=None` means "go back to the
            grouping file this object defaults to" instead of "keep whatever
            grouping file was set before". That default is the
            `custom_grouping_filepath` given to AqEquil, or the built-in
            speciation groups if none was given.
        """

        if isinstance(custom_grouping_filepath, str):
            self.custom_grouping_filepath = custom_grouping_filepath
        elif reset:
            self.custom_grouping_filepath = self.default_grouping_filepath

        if (not isinstance(self.reactant_dict_scalar, dict)
                or self._speciation_group_source != self.custom_grouping_filepath):
            self._make_speciation_group_dict()
            self._speciation_group_source = self.custom_grouping_filepath

    def _make_speciation_group_dict(self):

        if isinstance(self.custom_grouping_filepath, str):
            with open(self.custom_grouping_filepath) as file:
                lines = [line.rstrip() for line in file]
        else:
            content = import_package_file('aqequil.databases', 'speciation_groups_WORM.txt')
            content = content.split("\n")
            lines = [line.rstrip() for line in content]

        # Blank lines, and comments: a grouping file generated from a thermodynamic
        # database (e.g. by the WORM correlator's Complicator export) records at the
        # top which database it goes with, and that is worth being able to read.
        lines = [l for l in lines
                 if l.strip() != "" and not l.strip().startswith("#")]

        bad_line_indices = [i+1 for i,l in enumerate(lines) if ":" not in l]
        if len(bad_line_indices) > 0:
            bad_lines = [l for l in lines if ":" not in l]
            self.err_handler.raise_exception("Line(s) in the speciation group file do not contain a colon ':'. Fix the formatting on line(s):\n"+str(bad_line_indices)+"\nLines to blame:\n"+str(bad_lines))

        # Creates a dictionary with this format:
        #
        # {...
        #  "sulfide 1": ["H2S", "HS-"],
        #  "sulfide 2": ["Pb(HS)2", "Ag(HS)2-", "Au(HS)2-"],
        #  "sulfide 3": ["Pb(HS)3-"],
        #  "iron(II) 1": [...],
        #  ...
        # }
        #
        # Used by calculate_energy() to calculate concentrations of limiting reactants.
        #
        # A group key repeated on two lines is MERGED. It used to overwrite, and a
        # dropped line is invisible: the shipped file repeats eight keys, so every
        # metal succinate, malonate, glutarate, adipate, pimelate, glycolate and
        # lactate complex was being pooled with nothing at all. Species repeated
        # within a group are dropped, because _grouped_molality sums the group's
        # columns and would otherwise count one of them twice.
        speciation_group_dict = collections.OrderedDict()
        self._group_parts = {}
        for line in lines:
            group_name, _, rest = line.strip().partition(":")
            group_name = group_name.strip()
            group_species = speciation_group_dict.setdefault(group_name, [])
            for sp in rest.split():
                if sp not in group_species:
                    group_species.append(sp)
            self._group_parts[group_name] = self._split_group_name(group_name)

        self.reactant_dict_scalar = speciation_group_dict

        # Creates a dictionary with this format:
        # {...
        #  "H2S": ["H2S", "HS-"],
        #  "HS-": ["H2S", "HS-"],
        #  "Pb(HS)2": ["Pb(HS)2", "Ag(HS)2-", "Au(HS)2-"],
        #  "Ag(HS)2-": ["Pb(HS)2", "Ag(HS)2-", "Au(HS)2-"],
        #  "Au(HS)2-": ["Pb(HS)2", "Ag(HS)2-", "Au(HS)2-"],
        #  ...
        # }
        #
        # Used to switch a user-specified limiting reactant to one within the same scalar group.
        speciation_group_dict_unpacked = {}
        for group_species in speciation_group_dict.values():
            for sp in group_species:
                speciation_group_dict_unpacked[sp] = group_species

        self.speciation_group_dict_unpacked = speciation_group_dict_unpacked

        # Creates a dictionary with this format:
        # {...
        #  "sulfate": ['SO4-2','HSO4-','PdSO4','RhSO4', ..., 'Ru(SO4)3-4'],
        #  "formate": ['formic-acid', 'formate', 'Am(For)+2', ..., 'Zn(For)2'],
        #  ...
        # }
        # Used when the user wants to exclude all organics except for acetate (e.g.)
        # and all species in the acetate speciation group
        sp_dict_groups = {}
        for group_name, group_species in speciation_group_dict.items():
            category = self._group_parts[group_name][0]
            pooled = sp_dict_groups.setdefault(category, [])
            for sp in group_species:
                if sp not in pooled:
                    pooled.append(sp)
        self.sp_dict_groups = copy.deepcopy(sp_dict_groups)


    @staticmethod
    def _split_group_name(group_name):
        """
        A group key split into its category and the moiety count, e.g.
        "sulfide 2" -> ("sulfide", 2.0) and "PO2 anion 1" -> ("PO2 anion", 1.0).

        The count is the last space-separated token, so a category name may contain
        spaces but must not end in a number. A key whose last token is not a number
        gets a count of None and is left out of every group total rather than
        bringing the calculation down.
        """
        category, _, scalar = group_name.rpartition(" ")
        try:
            scalar = float(scalar)
        except ValueError:
            return group_name, None
        return category, scalar


    def _match_grouped_species(self, s):
        """
        Match whether a species is in a speciation group, and get a list of relevant groups and their scalars.

        e.g.,

        self.reactant_dict_scalar = {...
                                     "sulfides 1": ["H2S", "HS-"],
                                     "sulfides 2": ["Pb(HS)2", "Ag(HS)2-", "Au(HS)2-"],
                                     "sulfides 3": ["Pb(HS)3-"],
                                     "ferrous iron 1": [...],
                                     ...
                                     }
        If `s` = "Ag(HS)2-"
        then `scalars` = [1, 2, 3]
        and `groups` = [["H2S", "HS-"], ["Pb(HS)2", ...], ["Pb(HS)3-"]]

        A species in two categories (830 of the 1487 species in the shipped file
        are in both a ligand group and a metal group) matches the FIRST of them in
        file order, so whether Fe(Ac)+ totals all acetate or all ferrous iron is
        decided by the order of the grouping file.
        """
        scalars = []
        groups = []

        for group_name, grp_list in self.reactant_dict_scalar.items():
            if s not in grp_list:
                continue
            category = self._group_parts[group_name][0]
            # Every count this category has, matched on the category EXACTLY. It
            # used to be a substring test over the group keys, which read the "2"
            # of "H-citrate 2" as a count of citrate and then raised KeyError
            # looking up the "citrate 2" that does not exist -- so one category
            # name containing another was enough to break a calculation.
            for other_name, other_list in self.reactant_dict_scalar.items():
                other_category, other_scalar = self._group_parts[other_name]
                if other_category != category or other_scalar is None:
                    continue
                scalars.append(other_scalar)
                groups.append(other_list)
            break

        return scalars, groups


    def _switch_limiting(self, limiting, stoich, species, lenient=False):
        if limiting != None:
            reactant_idx = [1 if i<0 else 0 for i in stoich]
            reactants = [species[i] for i,idx in enumerate(reactant_idx) if idx == 1]

            if limiting not in reactants and limiting in list(self.speciation_group_dict_unpacked.keys()):
                # if the user specifies a limiting reactant like "HCO3-" but the
                # reaction has "CO2" as a reactant, check the dict that
                # contains speciated groups and switch the limiting reactant to
                # the relevant limiting reactant to appear in reports,
                # e.g., "HCO3-" -> "CO2"
                lim_species_group_list = self.speciation_group_dict_unpacked[limiting]
                for s in lim_species_group_list:
                    if s in reactants:
                        if self.verbose > 0:
                            print("The specified limiting reactant", str(limiting),
                                  "has been switched to", str(s), "because the latter",
                                  "appears as a reactant in the reaction:")
                            if isnotebook() and hasattr(self, "format_reaction"):
                                _ = self.format_reaction(coeffs=stoich,
                                   names=species,
                                   formatted=True,
                                   charge_sign_at_end=True,
                                   show=True)
                            else:
                                equation = format_equation(species, stoich,
                                                      charge_sign_at_end=True)
                                try:
                                    print(equation)
                                except UnicodeEncodeError:
                                    # e.g., a Windows console that cannot encode
                                    # the rightwards arrow
                                    print(equation.replace("→", "->"))
                        limiting = s
                        break

            if limiting not in reactants and lenient:
                # if the user specifies "CO2" as the limiting reactant during a
                # batch calculation of many different reactions, some reactions
                # won't actually have "CO2" as a reactant. In this case, set
                # limiting to None so that CO2 will be limiting when applicable.
                limiting = None

        return limiting


    def _species_state(self, s):
        """
        Physical state ('aq', 'cr', 'liq', 'gas') of a species according to the
        thermodynamic database CSV.
        """
        csv_db = self.thermo.csv_db
        return list(csv_db[csv_db["name"] == s]["state"])[0]


    @staticmethod
    def _mineral_amounts_known(grams_minerals):
        """
        Whether `grams_minerals` says enough about how much mineral or liquid
        reactant is present for those reactants to be treated as potentially
        limiting.
        """
        return isinstance(grams_minerals, (pd.DataFrame, dict, float, int))


    def _mineral_molality(self, s, n_rows, grams_minerals):

        """
        Moles of a mineral or liquid reactant available at each sample or step
        of reaction progress.

        Unlike an aqueous species, a mineral has no molality in the EQ3/6
        output; how much of it a microbe has to work with is something the user
        has to state. `grams_minerals` is that statement, and it is converted to
        moles here using the formula mass from the thermodynamic database.

        Parameters
        ----------
        s : str
            Name of the mineral or liquid species.

        n_rows : int
            Number of samples (`Speciation`) or steps of reaction progress
            (`Mass_Transfer`) being evaluated.

        grams_minerals : float, int, dict, pandas.DataFrame, or None
            Grams of mineral present. A float or int applies to every mineral
            reactant; a dict gives a value per mineral; a dataframe with a
            column per mineral gives a value for each row. Anything else (e.g.
            None) means the amount is unknown, in which case NaN is returned and
            the mineral cannot be limiting. Note that this is not the same as a
            mass of 0 grams, which makes the mineral limiting straight away.

        Returns
        ----------
        list of float
            Moles of `s` available, one value per row.
        """

        if not self._mineral_amounts_known(grams_minerals):
            return [float("NaN")]*n_rows

        if pychnosz.thermo().element is None:
            # touching the thermodynamic system initializes it, and that is what
            # loads the element masses that pychnosz.mass() needs. Nothing can
            # be lost by doing this here because an uninitialized system has no
            # thermodynamic data in it yet.
            _ = pychnosz.thermo().OBIGT

        csv_db = self.thermo.csv_db
        sp_formula = list(csv_db[csv_db["name"] == s]["formula"])[0]
        sp_mass = pychnosz.mass(sp_formula)

        if isinstance(grams_minerals, pd.DataFrame):
            if s not in grams_minerals.columns:
                return [0]*n_rows
            sp_grams_list = list(grams_minerals[s])
            if len(sp_grams_list) != n_rows:
                self.err_handler.raise_exception("The '"+str(s)+"' column of "
                        "grams_minerals has "+str(len(sp_grams_list))+" values "
                        "but "+str(n_rows)+" are needed, one for each row of "
                        "results.")
            return [sp_grams/sp_mass for sp_grams in sp_grams_list]

        if isinstance(grams_minerals, dict):
            sp_grams = grams_minerals.get(s, 0)
        else:
            sp_grams = grams_minerals

        return [sp_grams/sp_mass]*n_rows


    def _grouped_molality(self, s, molal_df, as_written=False):

        """
        Total molality available to a reactant, summed over its speciation
        group and weighted by how many atoms of the grouped element each
        species contributes.

        For example, the sulfide group is
        1*(H2S, HS-) + 2*(Pb(HS)2, Ag(HS)2-, Au(HS)2-) + 3*(Pb(HS)3-),
        so a reaction written with either "H2S" or "HS-" draws on the same
        total inventory of dissolved sulfide.

        Parameters
        ----------
        s : str
            Name of the species as written in the reaction.

        molal_df : pandas.DataFrame
            Table of molalities with species names as columns. One row per
            sample (`Speciation`) or per step of reaction progress
            (`Mass_Transfer`).

        as_written : bool, default False
            If True, do not sum over the speciation group; report only the
            molality of `s` itself, scaled by its group's scalar. Use this
            when the reaction is meant to be limited by one particular form
            of a species rather than by the total.

        Returns
        ----------
        list of float
            Molality available to `s`, one value per row of `molal_df`.
        """

        self._set_speciation_groups()

        # check whether the species matches any of the groups in
        # reactant_dict_scalar and retrieve scalars and groups
        scalars, groups = self._match_grouped_species(s)

        if as_written:
            # if energy supplies are to be calculated as written with no
            # grouping, then use the current species and its scalar
            for i, sc in enumerate(scalars):
                if s in groups[i]:
                    scalars = [sc]
                    groups = [[s]]
                    break

        if len(scalars) == 0:
            # `s` does not belong to any speciation group (e.g. O2, H2, CH4)
            return list(molal_df[s])

        total_summed_scaled = [0]*molal_df.shape[0]

        for i,scalar in enumerate(scalars):
            col_subset = [col for col in groups[i] if col in molal_df.columns]
            scaled_df = molal_df[col_subset].apply(lambda x: x*scalar)
            summed_scaled = list(scaled_df.sum(axis=1, numeric_only=True))
            total_summed_scaled = [ii+iii for ii,iii in zip(total_summed_scaled, summed_scaled)]

        return total_summed_scaled


    def _select_limiting_reactant(self, species, stoich, s_molal_dict, step,
                                        limiting=None, allow_minerals=False,
                                        charge_sign_at_end=False):

        """
        Identify the limiting reactant of a reaction at one sample or one step
        of reaction progress.

        Parameters
        ----------
        species, stoich : lists
            The reaction, as passed to the energy calculation.

        s_molal_dict : dict
            Maps each species name to a list of molalities (already summed
            over speciation groups by `_grouped_molality`).

        step : int
            Index into the lists in `s_molal_dict`; the sample or step of
            reaction progress being evaluated.

        limiting : str, optional
            Force this species to be the limiting reactant instead of
            choosing one by concentration and stoichiometry.

        allow_minerals : bool, default False
            Permit 'cr' and 'liq' reactants to be limiting. Only meaningful
            when mineral abundances are known (see the `grams_minerals`
            parameter of `Speciation.calculate_energy`).

        charge_sign_at_end : bool, default False
            Display charge with sign after the number (e.g. SO4 2-) in the
            reported limiting reactant name.

        Returns
        ----------
        tuple
            (lr_name, lr_reported, lr_concentration, lr_stoich). All four are
            None/NaN when the reaction has no valid limiting reactant, or when
            a limiting reactant concentration is unknown (NaN).
        """

        none_result = (None, float("nan"), float("nan"), float("nan"))

        if isinstance(limiting, str):
            lrc_dict = {}
            lr_list = [limiting]
        else:
            lrc_dict = {}
            for i_s,s in enumerate(species):
                # identify valid limiting reactants and record concentrations
                # 1. negative coefficient (reactant)
                # 2. can't be OH-, H+, H2O
                # 3. can't be cr or liq unless mineral abundances are known
                if stoich[i_s] >= 0:
                    continue
                if s in ["H2O", "H+", "OH-"]:
                    continue
                if not allow_minerals and self._species_state(s) in ["cr", "liq"]:
                    continue
                lrc_dict[s] = s_molal_dict[s][step]/abs(stoich[i_s])

            if len(lrc_dict) == 0:
                return none_result

            lr_val = min(lrc_dict.values())

            # handle situations where there might be multiple limiting reactants
            lr_list = [str(k) for k,v in lrc_dict.items() if v == lr_val]

        # a reactant of unknown concentration makes the energy supply undefined
        if any([pd.isna(v) for v in lrc_dict.values()]):
            return none_result

        lr_list_formatted = [chemlabel(lr_name, charge_sign_at_end=charge_sign_at_end) for lr_name in lr_list]
        if len(lr_list_formatted) > 1:
            lr_reported = ", ".join(lr_list_formatted)
        else:
            lr_reported = lr_list[0]

        lr_name = lr_list[0] # doesn't matter which lr is used to calculate
        lr_concentration = s_molal_dict[lr_name][step]
        lr_stoich = -stoich[species.index(lr_name)]

        return lr_name, lr_reported, lr_concentration, lr_stoich
