#!/usr/bin/env python
"""
Compile EQ3/6 executables for the current platform.

This script is called during the wheel build process to compile the EQ3/6
Fortran executables into platform-specific binaries.
"""

import os
import sys
import subprocess
import shutil
import platform
from pathlib import Path


def find_gfortran():
    """Find gfortran compiler on the system."""
    # Check if explicit gfortran path is provided via environment variable
    explicit_path = os.environ.get('GFORTRAN_PATH')
    if explicit_path and os.path.isfile(explicit_path):
        print(f"Using explicit gfortran path from GFORTRAN_PATH: {explicit_path}")
        return explicit_path

    # Try common gfortran locations
    candidates = ['gfortran', 'gfortran-11', 'gfortran-12', 'gfortran-13', 'gfortran-14']

    for compiler in candidates:
        compiler_path = shutil.which(compiler)
        if compiler_path:
            return compiler_path

    raise RuntimeError(
        "gfortran compiler not found. Please install gfortran:\n"
        "  Ubuntu/Debian: sudo apt-get install gfortran\n"
        "  macOS: brew install gcc\n"
        "  Windows: Install MinGW-w64 with gfortran support"
    )


def compile_eq3_6():
    """Compile EQ3/6 source code to executables."""
    # Determine paths
    script_dir = Path(__file__).parent
    source_dir = script_dir / "aqequil" / "eq3_6_src"
    output_dir = script_dir / "aqequil" / "bin"

    # Verify source directory exists
    if not source_dir.exists():
        raise FileNotFoundError(f"EQ3/6 source not found: {source_dir}")

    # Create output directory if it doesn't exist
    output_dir.mkdir(parents=True, exist_ok=True)

    # Find gfortran
    gfortran = find_gfortran()
    print(f"Using compiler: {gfortran}")

    # Get the path to gfortran directory for make
    gfortran_dir = Path(gfortran).parent
    print(f"Gfortran directory: {gfortran_dir}")

    # Build EQ3/6 using Make
    print(f"Building EQ3/6 from {source_dir}")

    # Set environment variables for compilation
    env = os.environ.copy()

    # Ensure gfortran is in PATH
    if str(gfortran_dir) not in env.get('PATH', ''):
        env['PATH'] = f"{gfortran_dir}:{env.get('PATH', '')}"

    # Set FC to explicitly use gfortran
    env['FC'] = gfortran

    # Platform-specific FFLAGS
    if sys.platform == "darwin":
        # macOS: detect architecture
        archflags = os.environ.get('ARCHFLAGS', '')
        print(f"ARCHFLAGS environment variable: '{archflags}'")

        if 'x86_64' in archflags:
            env['FFLAGS'] = '-O3 -arch x86_64'
            print("Building for macOS x86_64")
        elif 'arm64' in archflags:
            env['FFLAGS'] = '-O3 -arch arm64'
            print("Building for macOS ARM64")
        else:
            # Fallback to machine architecture
            machine = platform.machine()
            arch = 'arm64' if machine == 'arm64' else 'x86_64'
            env['FFLAGS'] = f'-O3 -arch {arch}'
            print(f"Building for macOS {arch} (detected)")
    elif sys.platform == "win32":
        # Windows: static linking to avoid runtime dependencies
        env['FFLAGS'] = '-O3 -static'
        print("Building for Windows with static linking")
    else:
        # Linux
        env['FFLAGS'] = '-O3'
        print("Building for Linux")

    try:
        # Clean any previous builds
        print("Cleaning previous builds...")
        subprocess.run(
            ['make', 'clean'],
            cwd=source_dir,
            env=env,
            check=False,  # Don't fail if clean fails
            capture_output=True
        )

        # Build all executables
        print("Compiling EQ3/6 executables...")
        result = subprocess.run(
            ['make', 'all'],
            cwd=source_dir,
            env=env,
            check=True,
            capture_output=True,
            text=True
        )

        if result.stdout:
            print(result.stdout)
        if result.stderr:
            print(result.stderr, file=sys.stderr)

        # Copy executables to output directory
        source_bin_dir = source_dir / "bin"
        if not source_bin_dir.exists():
            raise FileNotFoundError(f"Build failed: {source_bin_dir} not found")

        print(f"Copying executables to {output_dir}")

        # Define expected executables
        if sys.platform == "win32":
            executables = ['eq3nr.exe', 'eq6.exe', 'eqpt.exe', 'xcon3.exe', 'xcon6.exe']
        else:
            executables = ['eq3nr', 'eq6', 'eqpt', 'xcon3', 'xcon6']

        # Copy each executable
        for exe in executables:
            src = source_bin_dir / exe
            dst = output_dir / exe
            if src.exists():
                shutil.copy2(src, dst)
                # Make executable on Unix-like systems
                if sys.platform != "win32":
                    dst.chmod(0o755)
                size_kb = dst.stat().st_size / 1024
                print(f"  [OK] {exe} ({size_kb:.1f} KB)")
            else:
                print(f"  [WARNING] {exe} not found, skipping")

        # Verify at least the main executables exist
        eq3nr = output_dir / ('eq3nr.exe' if sys.platform == "win32" else 'eq3nr')
        eq6 = output_dir / ('eq6.exe' if sys.platform == "win32" else 'eq6')

        if not eq3nr.exists() or not eq6.exists():
            raise RuntimeError("Build failed: eq3nr or eq6 not found")

        print(f"[OK] Successfully compiled EQ3/6 executables")
        return True

    except subprocess.CalledProcessError as e:
        print(f"[ERROR] Compilation failed!", file=sys.stderr)
        print(f"Exit code: {e.returncode}", file=sys.stderr)
        if e.stdout:
            print(f"stdout: {e.stdout}", file=sys.stderr)
        if e.stderr:
            print(f"stderr: {e.stderr}", file=sys.stderr)
        raise RuntimeError(f"Failed to compile EQ3/6: {e}")


def main():
    """Main entry point."""
    print("=" * 60)
    print("Compiling EQ3/6 executables for aqequil")
    print("=" * 60)
    print(f"Platform: {sys.platform}")
    print(f"Python: {sys.version.split()[0]}")
    print(f"Machine: {platform.machine()}")

    try:
        compile_eq3_6()
        print("=" * 60)
        print("Compilation successful!")
        print("=" * 60)
        return 0
    except Exception as e:
        print("=" * 60, file=sys.stderr)
        print(f"ERROR: {e}", file=sys.stderr)
        print("=" * 60, file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
