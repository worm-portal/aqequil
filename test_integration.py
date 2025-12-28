#!/usr/bin/env python
"""
Integration test for aqequil package.

This test verifies that:
1. EQ3/6 executables are properly bundled and accessible
2. The package can perform basic speciation calculations
3. Results can be accessed and plotted

This is run during the wheel build process to ensure the package works correctly.
"""

import sys
import os


def test_bundled_executables():
    """Test that bundled EQ3/6 executables are found."""
    print("=" * 60)
    print("Test 1: Checking for bundled EQ3/6 executables")
    print("=" * 60)

    import aqequil
    from aqequil.AqSpeciation import _get_bundled_exe_path

    bin_path = _get_bundled_exe_path()
    if bin_path is None:
        print("[FAIL] No bundled executables found!")
        return False

    print(f"[OK] Found bundled executables at: {bin_path}")

    # Check for required executables
    if sys.platform == "win32":
        exe_names = ['eq3nr.exe', 'eq6.exe', 'eqpt.exe']
    else:
        exe_names = ['eq3nr', 'eq6', 'eqpt']

    for exe in exe_names:
        exe_path = os.path.join(bin_path, exe)
        if os.path.isfile(exe_path):
            size_mb = os.path.getsize(exe_path) / (1024 * 1024)
            print(f"  [OK] {exe} ({size_mb:.2f} MB)")
        else:
            print(f"  [FAIL] {exe} not found!")
            return False

    return True


def test_import_and_basic_usage():
    """Test basic aqequil import and usage."""
    print("\n" + "=" * 60)
    print("Test 2: Testing aqequil import and basic usage")
    print("=" * 60)

    try:
        import aqequil
        print("[OK] Successfully imported aqequil")

        # Get test data path
        from aqequil.test_data import get_test_data_path
        test_csv = get_test_data_path("input_example_wrm.csv")

        if not os.path.isfile(test_csv):
            print(f"[FAIL] Test data not found at: {test_csv}")
            return False
        print(f"[OK] Found test data at: {test_csv}")

        return True

    except Exception as e:
        print(f"[FAIL] Error during import/setup: {e}")
        import traceback
        traceback.print_exc()
        return False


def test_speciation():
    """Test running a speciation calculation."""
    print("\n" + "=" * 60)
    print("Test 3: Running speciation calculation")
    print("=" * 60)

    try:
        import aqequil
        from aqequil.test_data import get_test_data_path
        import pandas as pd

        print("Creating AqEquil instance...")
        ae = aqequil.AqEquil(db="WORM", exclude_organics=True, verbose=0)
        print("[OK] AqEquil instance created")

        test_csv = get_test_data_path("input_example_wrm.csv")
        print(f"Running speciation on {test_csv}...")

        # Run speciation
        speciation = ae.speciate(
            input_filename=test_csv,
            exclude=["Year", "Area"]
        )
        print("[OK] Speciation completed")

        # Check that we got results
        if not hasattr(speciation, 'sample_data'):
            print("[FAIL] Speciation object missing 'sample_data' attribute")
            return False

        if "Bison Pool" not in speciation.sample_data:
            print("[FAIL] 'Bison Pool' sample not found in results")
            return False

        print("[OK] Found 'Bison Pool' sample in results")

        # Check that aq_distribution is accessible and is a DataFrame
        aq_dist = speciation.sample_data["Bison Pool"]["aq_distribution"]
        if not isinstance(aq_dist, pd.DataFrame):
            print(f"[FAIL] aq_distribution is not a DataFrame, got {type(aq_dist)}")
            return False

        print(f"[OK] aq_distribution is a DataFrame with {len(aq_dist)} rows")

        # Test plotting (just check it doesn't crash, don't display)
        print("Testing plot_mass_contribution...")
        fig = speciation.plot_mass_contribution("HCO3-")
        if fig is None:
            print("[FAIL] plot_mass_contribution returned None")
            return False

        print("[OK] plot_mass_contribution executed successfully")

        return True

    except Exception as e:
        print(f"[FAIL] Error during speciation: {e}")
        import traceback
        traceback.print_exc()
        return False


def main():
    """Run all integration tests."""
    print("\n" + "=" * 60)
    print("AQEQUIL INTEGRATION TEST SUITE")
    print("=" * 60)
    print(f"Python: {sys.version.split()[0]}")
    print(f"Platform: {sys.platform}")
    if hasattr(os, 'uname'):
        print(f"Architecture: {os.uname().machine}")
    elif sys.platform == 'win32':
        import platform
        print(f"Architecture: {platform.machine()}")
    print("=" * 60)

    tests = [
        ("Bundled Executables", test_bundled_executables),
        ("Import and Basic Usage", test_import_and_basic_usage),
        ("Speciation Calculation", test_speciation),
    ]

    results = []
    for test_name, test_func in tests:
        try:
            result = test_func()
            results.append((test_name, result))
        except Exception as e:
            print(f"\n[FAIL] Test '{test_name}' crashed with exception: {e}")
            import traceback
            traceback.print_exc()
            results.append((test_name, False))

    # Print summary
    print("\n" + "=" * 60)
    print("TEST SUMMARY")
    print("=" * 60)
    for test_name, result in results:
        status = "[PASS]" if result else "[FAIL]"
        print(f"{status} {test_name}")

    print("=" * 60)

    # Return exit code
    all_passed = all(result for _, result in results)
    if all_passed:
        print("\n✓ All tests passed!")
        return 0
    else:
        print("\n✗ Some tests failed!")
        return 1


if __name__ == "__main__":
    sys.exit(main())
