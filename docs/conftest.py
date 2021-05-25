"""pytest fixtures for documentation tests using `--doctest-modules`."""

import warnings

import matplotlib
import pytest


@pytest.fixture(scope='session', autouse=True)
def set_mpl():
    """Avoid matplotlib windows popping up."""
    matplotlib.use('agg', force=True)
    warnings.filterwarnings(
        'ignore',
        category=UserWarning,
        message=(
            'Matplotlib is currently using agg, which is a non-GUI backend, so cannot show the '
            'figure.'
        ),
    )
