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


@pytest.fixture
def keep_test_output(request):
    return request.config.getoption("--keep-test-output")


@pytest.fixture(scope='module')
def db_path(request):
    return request.config.getoption("--db-path")
