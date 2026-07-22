import pytest

DEFAULT_DB_PATH = '/Users/arielmelchior/Documents/projects/DoChaP/DoChaP-web/DB_merged.sqlite'


def pytest_addoption(parser):
    parser.addoption(
        "--keep-test-output",
        action="store_true",
        default=False,
        help="Keep the generated CSV/PDF output directories created by the "
             "reference-comparison tests instead of deleting them after a "
             "successful comparison.",
    )
    parser.addoption(
        "--db-path",
        action="store",
        default=DEFAULT_DB_PATH,
        help="Path to the DoChaP merged sqlite database to run the tests "
             f"against. Defaults to {DEFAULT_DB_PATH}.",
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
    return request.config.getoption("--db-path")
