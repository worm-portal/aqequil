from typing import List, Dict, Tuple, Optional
import pandas as pd
import os
from dataclasses import dataclass
from collections import Counter


@dataclass
class ValidationResult:
    """Container for validation results"""
    is_valid: bool
    errors: List[str]
    warnings: List[str]
    sample_temps: List[float]
    sample_pressures: List[str]
    redox_flag: str


class ValidationError(Exception):
    """Custom exception for validation errors"""
    def __init__(self, errors: List[str]):
        self.errors = errors
        super().__init__(f"Validation failed with {len(errors)} errors")


class InputValidator:
    """
    Validates sample input files for aqueous speciation calculations.
    
    Separates validation logic from the main AqEquil class for better 
    maintainability and testability.
    """
    
    # Class constants
    VALID_SUBHEADERS = [
        "degC", "ppm", "ppb", "Suppressed", "Molality", "Molarity", 
        "mg/L", "mg/kg.sol", "Alk., eq/kg.H2O", "Alk., eq/L", 
        "Alk., eq/kg.sol", "Alk., mg/L CaCO3", "Alk., mg/L HCO3-", 
        "Log activity", "Log act combo", "Log mean act", "pX", "pH", 
        "pe", "pHCl", "pmH", "pmX", "Hetero. equil.", "Homo. equil.", 
        "Make non-basis", "logfO2", "Mineral", "bar", "volts"
    ]
    
    FIXED_SPECIES = ["H2O", "H+", "O2(g)", "water", "e-", "OH-", "O2", "H2O(g)"]
    SYSTEM_COLUMNS = ['Temperature', 'logfO2', 'pH', 'Pressure', 'Eh', 'pe']
    REQUIRED_COLUMNS = ["Temperature"]
    
    REDOX_COLUMNS = {
        "logfO2": "logfO2",
        "Eh": "volts", 
        "pe": "pe",
        "O2(g)": "Hetero. equil."
    }

    def __init__(self, verbose: int = 1):
        self.verbose = verbose
        self.errors: List[str] = []
        self.warnings: List[str] = []

    def validate_input_file(self, 
                          filename: str,
                          exclude: List[str],
                          db_species: List[str],
                          charge_balance_on: str = "none",
                          redox_flag: str = "auto") -> ValidationResult:
        """
        Main validation method that orchestrates all validation checks.
        
        Returns:
            ValidationResult: Contains validation status and extracted data
        """
        self._reset_validation_state()
        
        # Basic file validation
        df = self._validate_file_structure(filename)
        
        # Header and row validation
        self._validate_headers_and_rows(df)
        
        # Column validation after exclusions
        df_clean = self._apply_exclusions(df, exclude)
        
        # Species validation
        self._validate_species_in_database(df_clean, db_species)
        
        # Subheader validation
        self._validate_subheaders(df_clean)
        
        # Required columns validation
        self._validate_required_columns(df_clean, exclude)
        
        # Charge balance validation
        self._validate_charge_balance(df_clean, charge_balance_on)
        
        # Extract sample data
        sample_temps, sample_pressures = self._extract_sample_conditions(df)
        
        # Redox flag validation and auto-detection
        redox_flag = self._validate_and_detect_redox_flag(df_clean, redox_flag)
        
        # Compile results
        is_valid = len(self.errors) == 0
        
        if not is_valid and self.verbose > 0:
            self._print_validation_summary(filename)
            
        return ValidationResult(
            is_valid=is_valid,
            errors=self.errors.copy(),
            warnings=self.warnings.copy(),
            sample_temps=sample_temps,
            sample_pressures=sample_pressures,
            redox_flag=redox_flag
        )

    def _reset_validation_state(self):
        """Reset validation state for new validation run"""
        self.errors.clear()
        self.warnings.clear()

    def _validate_file_structure(self, filename: str) -> pd.DataFrame:
        """Validate basic file structure and load data"""
        if not self._file_exists(filename):
            self.errors.append(f"Cannot locate input file {filename}")
            return pd.DataFrame()
            
        try:
            df = pd.read_csv(filename, header=None)
        except Exception as e:
            self.errors.append(f"Cannot read CSV file {filename}: {str(e)}")
            return pd.DataFrame()
            
        if df.shape[0] <= 2:
            self.errors.append(
                f"File {filename} must contain at least three rows: "
                "column names, subheaders, and sample data"
            )
            
        return df

    def _validate_headers_and_rows(self, df: pd.DataFrame):
        """Validate column headers and sample names"""
        if df.empty:
            return
            
        # Check for blank headers (except first column)
        col_list = list(df.iloc[0, 1:])
        if any(pd.isna(x) for x in col_list):
            self.errors.append(
                "One or more columns have blank headers. "
                "Remove empty columns or provide names for all headers."
            )
        
        # Check for duplicate headers
        duplicate_headers = self._find_duplicates(col_list)
        if duplicate_headers:
            self.errors.append(
                f"Duplicate column names found: {duplicate_headers}"
            )
        
        # Set proper column names and validate rows
        df.columns = df.iloc[0]
        sample_names = list(df.iloc[1:, 0])
        
        # Check for blank sample names
        if any(pd.isna(x) for x in sample_names):
            self.errors.append(
                "One or more samples have blank names. "
                "Provide names for all samples in the first column."
            )
        
        # Check for duplicate sample names
        duplicate_samples = self._find_duplicates(sample_names)
        if duplicate_samples:
            self.errors.append(
                f"Duplicate sample names found: {duplicate_samples}"
            )
            
        # Check for leading/trailing spaces in sample names
        invalid_names = [name for name in sample_names 
                        if isinstance(name, str) and (name.startswith(' ') or name.endswith(' '))]
        if invalid_names:
            self.errors.append(
                f"Sample names have leading/trailing spaces: {invalid_names}"
            )

    def _apply_exclusions(self, df: pd.DataFrame, exclude: List[str]) -> pd.DataFrame:
        """Apply column exclusions and return cleaned dataframe"""
        if df.empty:
            return df
            
        df_clean = df.iloc[:, 1:].copy()  # Remove first column (sample names)
        
        for col in exclude:
            if col == df.columns[0]:  # Skip sample column
                continue
            if col in df_clean.columns:
                df_clean = df_clean.drop(col, axis=1)
            else:
                self.errors.append(f"Cannot exclude column '{col}' - not found in file")
                
        return df_clean

    def _validate_species_in_database(self, df_clean: pd.DataFrame, db_species: List[str]):
        """Validate that column species exist in the database"""
        if df_clean.empty or not db_species:
            return
            
        all_valid_species = set(db_species + self.SYSTEM_COLUMNS + self.FIXED_SPECIES)
        
        for species in df_clean.columns:
            if species not in all_valid_species:
                if species == 'pH':
                    self.errors.append(
                        "Please rename 'pH' column to 'H+' with subheader unit 'pH'"
                    )
                else:
                    self.errors.append(
                        f"Column '{species}' not found in thermodynamic database. "
                        f"Add to 'exclude' parameter if this is metadata."
                    )

    def _validate_subheaders(self, df_clean: pd.DataFrame):
        """Validate subheader units"""
        if df_clean.empty:
            return
            
        subheaders = df_clean.iloc[0]
        for i, subheader in enumerate(subheaders):
            if subheader not in self.VALID_SUBHEADERS:
                column_name = df_clean.columns[i]
                self.errors.append(
                    f"Invalid subheader '{subheader}' for column '{column_name}'. "
                    f"Valid subheaders: {self.VALID_SUBHEADERS}"
                )

    def _validate_required_columns(self, df_clean: pd.DataFrame, exclude: List[str]):
        """Validate that required columns are present"""
        for required_col in self.REQUIRED_COLUMNS:
            if required_col not in df_clean.columns and required_col not in exclude:
                self.errors.append(
                    f"Required column '{required_col}' not found. "
                    f"Include column with '{required_col}' header and 'degC' subheader."
                )

    def _validate_charge_balance(self, df_clean: pd.DataFrame, charge_balance_on: str):
        """Validate charge balance settings"""
        invalid_charge_balance = ['pH', 'Temperature', 'logfO2']
        
        if charge_balance_on == 'pH':
            self.errors.append("To balance charge on pH, use charge_balance_on='H+'")
        elif charge_balance_on in invalid_charge_balance:
            self.errors.append(f"Cannot balance charge on {charge_balance_on}")
        elif (charge_balance_on != "none" and 
              charge_balance_on not in df_clean.columns):
            self.errors.append(
                f"Charge balance species '{charge_balance_on}' not found in input file"
            )

    def _validate_and_detect_redox_flag(self, df_clean: pd.DataFrame, redox_flag: str) -> str:
        """Validate redox flag settings and auto-detect if needed"""
        if redox_flag == "auto":
            redox_flag = self._auto_detect_redox_flag(df_clean)
        elif redox_flag == "O2(g)":
            self._validate_o2g_redox_flag(df_clean)
            
        return redox_flag

    def _auto_detect_redox_flag(self, df_clean: pd.DataFrame) -> str:
        """Auto-detect appropriate redox flag from available columns"""
        for redox_col, expected_subheader in self.REDOX_COLUMNS.items():
            if redox_col in df_clean.columns:
                if redox_col == "O2":
                    if self.verbose > 0:
                        print("Using O2 column for redox state")
                    return "logfO2"
                elif (redox_col not in ["logfO2", "O2(g)"] or 
                      df_clean[redox_col].iloc[0] == expected_subheader):
                    if self.verbose > 0:
                        print(f"Using '{redox_col}' column for redox state")
                    return redox_col
        
        if self.verbose > 0:
            self.warnings.append(
                "No redox column detected. Available options: " +
                ", ".join(f"'{k}' (subheader: '{v}')" for k, v in self.REDOX_COLUMNS.items())
            )
        return "logfO2"

    def _validate_o2g_redox_flag(self, df_clean: pd.DataFrame):
        """Validate O2(g) redox flag setup"""
        if "O2(g)" not in df_clean.columns:
            self.errors.append("Redox flag 'O2(g)' selected but no O2(g) column found")
        elif df_clean["O2(g)"].iloc[0] != "Hetero. equil.":
            current_subheader = df_clean["O2(g)"].iloc[0]
            self.errors.append(
                f"O2(g) column subheader is '{current_subheader}' but should be "
                f"'Hetero. equil.' for redox flag 'O2(g)'"
            )

    def _extract_sample_conditions(self, df: pd.DataFrame) -> Tuple[List[float], List[str]]:
        """Extract temperature and pressure data from samples"""
        if df.empty or "Temperature" not in df.columns:
            return [], []
            
        try:
            temps = [float(t) for t in df["Temperature"].iloc[1:]]
        except (ValueError, TypeError):
            self.errors.append("Invalid temperature values found")
            temps = []
            
        if "Pressure" in df.columns:
            try:
                pressures = [float(p) if str(p).lower() != 'psat' else 'psat' 
                           for p in df["Pressure"].iloc[1:]]
            except (ValueError, TypeError):
                self.errors.append("Invalid pressure values found")
                pressures = ['psat'] * len(temps)
        else:
            pressures = ['psat'] * len(temps)
            
        return temps, pressures

    def _find_duplicates(self, items: List) -> List:
        """Find duplicate items in a list"""
        counter = Counter(items)
        return [item for item, count in counter.items() if count > 1]

    def _file_exists(self, filename: str) -> bool:
        """Check if file exists and is readable"""
        return os.path.exists(filename) and os.path.isfile(filename)

    def _print_validation_summary(self, filename: str):
        """Print validation error summary"""
        print(f"\nValidation errors found in {filename}:")
        for i, error in enumerate(self.errors, 1):
            print(f"{i}. {error}")
        
        if self.warnings:
            print(f"\nWarnings:")
            for warning in self.warnings:
                print(f"- {warning}")


# Usage example:
def validate_sample_input(filename: str, 
                         exclude: List[str],
                         db_species: List[str],
                         charge_balance_on: str = "none",
                         redox_flag: str = "auto",
                         verbose: int = 1) -> ValidationResult:
    """
    Convenience function to validate input file
    
    Returns:
        ValidationResult with validation status and extracted data
    """
    validator = InputValidator(verbose=verbose)
    result = validator.validate_input_file(
        filename=filename,
        exclude=exclude,
        db_species=db_species,
        charge_balance_on=charge_balance_on,
        redox_flag=redox_flag
    )
    
    if not result.is_valid:
        raise ValidationError(result.errors)
        
    return result