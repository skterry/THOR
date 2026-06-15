"""
THOR — build/install hooks.

Running ``pip install -e .`` from the top-level directory triggers two custom
setup steps, in addition to the normal (no-op) Python install:

  1. Compile the Fortran sources in ``src/thor/`` (hst1pass.F, thor_go.F) into
     bare executables, using whatever Fortran compiler the user has installed
     (gfortran, ifort/ifx, flang, ...). This is delegated to the Makefile in
     that directory, which auto-detects the compiler.

  2. Download the standard geometric-distortion correction (GDC) reference
     files from STScI directly into ``src/thor/GDCs/``.

  3. Download the two large (~1 GB each) data archives from Google Drive into
     ``data/``, showing a progress bar for each.

All steps are idempotent: already-built executables and already-downloaded
files are left alone, so re-running the install is cheap.

These steps run from the ``build_py`` command (used by both regular and
editable installs under modern pip/setuptools) and from the ``develop`` command
(used by the legacy ``setup.py develop`` editable path).
"""

import os
import shutil
import subprocess
import urllib.request
import zipfile
from pathlib import Path

from setuptools import setup
from setuptools.command.build_py import build_py
from setuptools.command.develop import develop

# --------------------------------------------------------------------------- #
# Configuration
# --------------------------------------------------------------------------- #

PROJECT_ROOT = Path(__file__).resolve().parent
THOR_DIR = PROJECT_ROOT / "src" / "thor"
HAMRR_DIR = PROJECT_ROOT / "src" / "hamrr"
GDC_DIR = THOR_DIR / "GDCs"
DATA_DIR = PROJECT_ROOT / "data"
EXAMPLE_THOR_DIR = PROJECT_ROOT / "example" / "thor_HD138"
EXAMPLE_ACS_DIR = EXAMPLE_THOR_DIR / "ACS.XYM"
EXAMPLE_WFC3_DIR = EXAMPLE_THOR_DIR / "WFC3.XYM"
EXAMPLE_HAMRR_DIR = PROJECT_ROOT / "example" / "hamrr"

# Fortran sources -> output executable name (no extension).
FORTRAN_TARGETS = [
    ("hst1pass.F", "hst1pass"),
    ("thor_go.F", "thor_go"),
]

# Files to copy from src/thor/ after compilation, as (filename, destination dir)
# pairs. The compiled executables land in example/thor_HD138/, while the reduce
# scripts go into per-detector sub-directories; collate_thor.src is needed in
# both, so it is listed twice.
THOR_EXAMPLE_FILES = [
    ("hst1pass", EXAMPLE_THOR_DIR),
    ("thor_go", EXAMPLE_THOR_DIR),
    ("reduce_acs.src", EXAMPLE_ACS_DIR),
    ("reduce_wfc3.src", EXAMPLE_WFC3_DIR),
    ("collate_thor.src", EXAMPLE_ACS_DIR),
    ("collate_thor.src", EXAMPLE_WFC3_DIR),
]

# Files to copy from src/hamrr/ -> example/hamrr/.
HAMRR_EXAMPLE_FILES = [
    "hamrr.py",
    "image_cut.py",
    "params.in",
]

# Fortran compilers to look for, in order of preference, when falling back to
# direct compilation (i.e. when `make` is unavailable). The Makefile uses the
# same list.
FORTRAN_COMPILERS = ["gfortran", "ifx", "ifort", "flang-new", "flang", "g77", "f77"]

# Standard geometric-distortion correction (GDC) reference files. These are
# direct download links (plain HTTP, not Google Drive), fetched into
# src/thor/GDCs/. The destination filename is taken from the URL's basename.
GDC_BASE_URL = "https://www.stsci.edu/~jayander/HST1PASS/LIB/GDCs/STDGDCs"
GDC_URLS = [
    f"{GDC_BASE_URL}/WFC3UV/STDGDC_WFC3UV_F606W.fits",
    f"{GDC_BASE_URL}/WFC3UV/STDGDC_WFC3UV_F814W.fits",
    f"{GDC_BASE_URL}/ACSWFC/STDGDC_OFFICIAL_JFRAME_ACSWFC_F606W.fits",
    f"{GDC_BASE_URL}/ACSWFC/STDGDC_OFFICIAL_JFRAME_ACSWFC_F814W.fits",
]

# (Google Drive file id, destination filename under data/, extracted path)
# When the third element is non-None, the downloaded archive is extracted into
# data/ and the .zip is deleted afterwards. The extracted path is the entry the
# archive expands to (relative to data/); its presence is what makes the
# download+extract idempotent, since the .zip itself is removed.
DATA_FILES = [
    ("17peXwwP6HzZrwOqFuYp-16odfVQkTJof", "field_HD138.zip", "field_HD138"),
    ("1VbbdWCv8LdSW3Ph50Rt0USqdoIZsrAhF", "thor_hst_wfc3_acs_bulge_v0.2_cat.fits.zip", None),
]

# Run the custom steps at most once per interpreter process, even if more than
# one command (build_py + develop) is invoked.
_STEPS_DONE = False


# --------------------------------------------------------------------------- #
# Helpers
# --------------------------------------------------------------------------- #

def _banner(msg):
    line = "=" * 70
    print("\n" + line + "\n" + msg + "\n" + line, flush=True)


def _find_fortran_compiler():
    for compiler in FORTRAN_COMPILERS:
        path = shutil.which(compiler)
        if path:
            return compiler
    return None


def compile_fortran():
    """Compile the Fortran sources into executables (idempotent)."""
    _banner("THOR setup [1/4]: compiling Fortran sources in src/thor/")

    if not THOR_DIR.is_dir():
        print(f"  WARNING: {THOR_DIR} not found; skipping Fortran build.")
        return

    # Prefer the Makefile (it auto-detects the compiler and only rebuilds when a
    # source is newer than its executable). Fall back to compiling directly if
    # `make` is not on PATH.
    if shutil.which("make"):
        try:
            subprocess.check_call(["make", "-C", str(THOR_DIR)])
            return
        except subprocess.CalledProcessError as exc:
            print(f"  WARNING: `make` failed (exit {exc.returncode}); "
                  f"trying direct compilation.")

    compiler = _find_fortran_compiler()
    if not compiler:
        print("  WARNING: no Fortran compiler found on PATH "
              f"(looked for: {', '.join(FORTRAN_COMPILERS)}).")
        print("           Install one (e.g. gfortran) and re-run "
              "`pip install -e .`, or run `make` in src/thor/ yourself.")
        return

    print(f"  Using Fortran compiler: {compiler}")
    for source, exe in FORTRAN_TARGETS:
        src_path = THOR_DIR / source
        exe_path = THOR_DIR / exe
        if not src_path.is_file():
            print(f"  WARNING: source {src_path} missing; skipping.")
            continue
        # Skip if up to date.
        if exe_path.exists() and exe_path.stat().st_mtime >= src_path.stat().st_mtime:
            print(f"  {exe} is up to date.")
            continue
        print(f"  Compiling {source} -> {exe}")
        try:
            subprocess.check_call([compiler, source, "-o", exe], cwd=str(THOR_DIR))
        except subprocess.CalledProcessError as exc:
            print(f"  WARNING: failed to compile {source} (exit {exc.returncode}).")


def copy_example_files():
    """Copy runtime files into the example directories (idempotent)."""
    _banner("THOR setup [2/4]: copying files into example/")

    # Flatten into a single list of (source dir, filename, destination dir)
    # copies. THOR files carry their own per-file destination; HAMRR files all
    # share one destination directory.
    copies = [(THOR_DIR, name, dst_dir) for name, dst_dir in THOR_EXAMPLE_FILES]
    copies += [(HAMRR_DIR, name, EXAMPLE_HAMRR_DIR) for name in HAMRR_EXAMPLE_FILES]

    for src_dir, name, dst_dir in copies:
        if not src_dir.is_dir():
            print(f"  WARNING: {src_dir} not found; skipping.")
            continue
        src = src_dir / name
        if not src.exists():
            print(f"  WARNING: {src} not found; skipping.")
            continue
        dst_dir.mkdir(parents=True, exist_ok=True)
        dst = dst_dir / name
        shutil.copy2(src, dst)
        print(f"  Copied {src.relative_to(PROJECT_ROOT)} -> "
              f"{dst.relative_to(PROJECT_ROOT)}")


def _extract_and_cleanup(archive):
    """Unzip ``archive`` into data/, then delete the archive.

    A failed extraction leaves the archive in place so a re-run can retry.
    """
    print(f"  Extracting {archive.name} into {DATA_DIR}/ ...")
    try:
        with zipfile.ZipFile(archive) as zf:
            # Skip __MACOSX resource-fork entries added by macOS zip tools.
            members = [m for m in zf.namelist() if not m.startswith("__MACOSX/")]
            zf.extractall(DATA_DIR, members=members)
    except (zipfile.BadZipFile, OSError) as exc:
        print(f"  WARNING: failed to extract {archive.name}: {exc}")
        return

    try:
        archive.unlink()
        print(f"  Extracted and removed {archive.name}.")
    except OSError as exc:
        print(f"  WARNING: extracted {archive.name} but could not remove it: {exc}")


def _progress_hook(filename):
    """Return a urlretrieve reporthook that prints download progress."""
    def hook(block_num, block_size, total_size):
        if total_size <= 0:
            return
        got = min(block_num * block_size, total_size)
        pct = 100.0 * got / total_size
        # \r keeps it on one line for ttys; harmless (extra newlines) under pip.
        print(f"\r  {filename}: {pct:5.1f}% "
              f"({got / 1e6:.0f}/{total_size / 1e6:.0f} MB)", end="", flush=True)
        if got >= total_size:
            print()
    return hook


def download_gdcs():
    """Download the GDC reference files from STScI into src/thor/GDCs/ (idempotent)."""
    _banner("THOR setup [3/4]: downloading GDC reference files to src/thor/GDCs/ "
            "(~300 MB each)")

    GDC_DIR.mkdir(parents=True, exist_ok=True)

    for url in GDC_URLS:
        filename = url.rsplit("/", 1)[-1]
        dest = GDC_DIR / filename

        if dest.exists() and dest.stat().st_size > 0:
            print(f"  {filename} already present "
                  f"({dest.stat().st_size / 1e6:.1f} MB); skipping.")
            continue

        print(f"  Downloading {filename} ...", flush=True)
        tmp = dest.with_suffix(dest.suffix + ".part")
        try:
            urllib.request.urlretrieve(url, tmp, reporthook=_progress_hook(filename))
        except Exception as exc:  # noqa: BLE001 - keep install resilient
            print(f"  WARNING: failed to download {filename}: {exc}")
            tmp.unlink(missing_ok=True)
            continue
        # Only publish the final file once the download finishes cleanly, so an
        # interrupted run never leaves a truncated .fits in place.
        tmp.replace(dest)
        print(f"  Downloaded {filename} ({dest.stat().st_size / 1e6:.1f} MB).",
              flush=True)


def download_data():
    """Download the large data archives from Google Drive (idempotent)."""
    _banner("THOR setup [4/4]: downloading data archives to data/ "
            "(~1 GB each — this may take a while)")

    DATA_DIR.mkdir(parents=True, exist_ok=True)

    try:
        import gdown
    except ImportError:
        print("  WARNING: `gdown` is not installed, so the data files cannot be "
              "downloaded automatically.")
        print("           Install it with `pip install gdown` and re-run "
              "`pip install -e .`.")
        return

    for file_id, filename, extracted in DATA_FILES:
        dest = DATA_DIR / filename
        extracted_path = (DATA_DIR / extracted) if extracted else None

        # Already-done checks. For extracted archives the .zip is gone, so we key
        # off the extracted path instead of the archive itself.
        if extracted_path is not None and extracted_path.exists():
            print(f"  {filename} already downloaded and extracted ({extracted}/); "
                  "skipping.")
            continue
        if extracted_path is None and dest.exists() and dest.stat().st_size > 0:
            print(f"  {filename} already present ({dest.stat().st_size / 1e6:.0f} MB); "
                  "skipping.")
            continue

        # Download the archive unless it's already sitting on disk (e.g. a
        # previous run downloaded it but didn't get to extract it).
        if not (dest.exists() and dest.stat().st_size > 0):
            print(f"  Downloading {filename} ... (this may take several minutes)",
                  flush=True)
            try:
                # quiet=False requests a tqdm progress bar; note that tqdm
                # suppresses the bar when stdout is not a tty (e.g. under pip),
                # so we print file size afterwards as a fallback confirmation.
                gdown.download(id=file_id, output=str(dest), quiet=False)
            except Exception as exc:  # noqa: BLE001 - keep install resilient
                print(f"  WARNING: failed to download {filename}: {exc}")
                # Remove any partial file so a re-run starts clean.
                if dest.exists() and dest.stat().st_size == 0:
                    dest.unlink(missing_ok=True)
                continue
            if dest.exists():
                print(f"  Downloaded {filename} ({dest.stat().st_size / 1e6:.0f} MB).",
                      flush=True)

        # Extract into data/ and drop the archive, if requested.
        if extracted_path is not None:
            _extract_and_cleanup(dest)


def run_thor_setup():
    """Run both custom setup steps once."""
    global _STEPS_DONE
    if _STEPS_DONE:
        return
    _STEPS_DONE = True

    # Allow opting out (e.g. for CI or quick reinstalls).
    if os.environ.get("THOR_SKIP_SETUP"):
        print("THOR_SKIP_SETUP set; skipping Fortran build and data download.")
        return

    compile_fortran()
    copy_example_files()
    download_gdcs()
    download_data()
    _banner("THOR setup complete.")


# --------------------------------------------------------------------------- #
# Custom commands
# --------------------------------------------------------------------------- #

class BuildPyCommand(build_py):
    """Standard + modern editable installs route through build_py."""

    def run(self):
        run_thor_setup()
        super().run()


class DevelopCommand(develop):
    """Legacy `setup.py develop` / older editable installs."""

    def run(self):
        run_thor_setup()
        super().run()


# The PEP 517 build backend (see thor_build.py / pyproject.toml) calls
# run_thor_setup() directly, which is what makes `pip install -e .` trigger the
# Fortran build and data download under modern pip. The cmdclass overrides below
# additionally cover the legacy `python setup.py develop` path; thanks to the
# idempotent steps, running via both paths is harmless.
#
# Guarded so that `import setup` (done by the build backend to reach
# run_thor_setup) does not kick off an actual build.
if __name__ == "__main__":
    setup(
        cmdclass={
            "build_py": BuildPyCommand,
            "develop": DevelopCommand,
        },
    )
