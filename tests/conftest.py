import pytest


def pytest_addoption(parser):
    parser.addoption(
        "--keep-test-output",
        action="store_true",
        default=False,
        help="Keep the generated CSV/PDF output directories created by the "
             "reference-comparison tests instead of deleting them after a "
             "successful comparison.",
    )


@pytest.fixture
def keep_test_output(request):
    return request.config.getoption("--keep-test-output")
