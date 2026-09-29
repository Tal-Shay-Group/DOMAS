"""Fix the timestamp matplotlib stamps into every PDF it writes.

Without this each regeneration rewrites all 130 golden PDFs with the current
clock - seven bytes in a 76KB file - so git reports every one of them as
changed and the diff that matters, the CSVs and the summaries, is buried in
the noise. The tests never noticed: they compare the extracted-text manifest,
not the bytes. SOURCE_DATE_EPOCH is the reproducible-builds convention for
exactly this, and matplotlib honours it.

Test-side only. A real DOMAS run still records the time it actually ran.
"""
import os

# 2024-01-01T00:00:00Z. Any fixed instant would do; a plausible date reads
# better than the epoch when the metadata is inspected.
os.environ.setdefault("SOURCE_DATE_EPOCH", "1704067200")

import pytest

# The database is not in this repository and its location differs per install, so
# there is no default worth baking in: DOMAS_TEST_DB names it, or --dochap does.
# Same variable generate_reference_outputs.py reads, so one export serves both.
DEFAULT_DB_PATH = os.environ.get('DOMAS_TEST_DB')


def pytest_addoption(parser):
    parser.addoption(
        "--keep-test-output",
        action="store_true",
        default=False,
        help="Keep the generated CSV/PDF output directories created by the "
             "reference-comparison tests instead of deleting them after a "
             "successful comparison.",
    )
    # pytest options must be double-dashed (its parser rejects single-dash long
    # names), so the DOMAS CLI's -dochap is spelled --dochap here. --db-path is
    # the older name for the same thing and still works.
    parser.addoption(
        "--dochap", "--db-path",
        action="store",
        dest="dochap",
        default=DEFAULT_DB_PATH,
        help="Path to the DoChaP merged sqlite database to run the tests "
             f"against. Defaults to $DOMAS_TEST_DB{' (currently unset)' if not DEFAULT_DB_PATH else ''}.",
    )
    parser.addoption(
        "--ioe-input-dir",
        action="store",
        default=None,
        help="Directory of SUPPA .ioe files to run the full-scale ioe comparison "
             "against (e.g. domas_extra/external_data/H_sapiens). The test is "
             "skipped unless this and --ioe-output-file are both given.",
    )
    parser.addoption(
        "--ioe-output-file",
        action="store",
        default=None,
        help="CSV the full-scale ioe run is compared against. Created from the "
             "run when it does not exist yet, so the first invocation "
             "bootstraps a baseline and later ones compare to it.",
    )
    parser.addoption(
        "--ioe-specie",
        action="store",
        default=None,
        help="Species the --ioe-input-dir data came from (human/mouse/rat). "
             "Required whenever --ioe-input-dir is given: SUPPA files carry no "
             "species field, so a default would just be a guess that DOMAS aborts "
             "on minutes into the run.",
    )
    parser.addoption(
        "--ioe-max-clusters",
        action="store",
        type=int,
        default=0,
        help="Cap the full-scale ioe run at the first N clusters (0 = no cap). "
             "A whole H_sapiens directory is a multi-hour run; this makes the "
             "comparison something you can do in a coffee break.",
    )
    parser.addoption(
        "--run-slow",
        action="store_true",
        default=False,
        help="Run tests marked @pytest.mark.slow (e.g. the full-dataset "
             "leafcutter golden comparison, ~17k clusters). Skipped by default.",
    )


def pytest_configure(config):
    config.addinivalue_line(
        "markers", "slow: slow test, opt in with --run-slow")

    # Fail at argument-parse time rather than minutes into a run, where the
    # species guard would otherwise abort it.
    if config.getoption("--ioe-input-dir") and not config.getoption("--ioe-specie"):
        raise pytest.UsageError(
            "--ioe-specie is required with --ioe-input-dir (human/mouse/rat): "
            "SUPPA .ioe files carry no species field, and DOMAS aborts when the "
            "stated species contradicts the gene ids.")


def pytest_collection_modifyitems(config, items):
    if config.getoption("--run-slow"):
        return
    skip_slow = pytest.mark.skip(reason="need --run-slow option to run")
    for item in items:
        if "slow" in item.keywords:
            item.add_marker(skip_slow)


@pytest.fixture
def keep_test_output(request):
    return request.config.getoption("--keep-test-output")


@pytest.fixture(scope='module')
def db_path(request):
    return request.config.getoption("dochap")


@pytest.fixture(scope='module')
def ioe_input_dir(request):
    return request.config.getoption("--ioe-input-dir")


@pytest.fixture(scope='module')
def ioe_output_file(request):
    return request.config.getoption("--ioe-output-file")


@pytest.fixture(scope='module')
def ioe_max_clusters(request):
    return request.config.getoption("--ioe-max-clusters")


@pytest.fixture(scope='module')
def ioe_specie(request):
    return request.config.getoption("--ioe-specie")
