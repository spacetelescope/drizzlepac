import crds
import os
import pytest
import re

__all__ = ['runslow', 'require_bigdata', 'slow',
           'not_under_travis', 'require_crds_context']

# pytest marker to mark resource-intensive tests that should not be
# executed with every commit.
slow = pytest.mark.slow

# Decorator to indicate slow tests
runslow = pytest.mark.skipif(
    not pytest.config.getoption("--runslow"),
    reason="need --runslow option to run"
)

# Decorator to indicate BIGDATA required
require_bigdata = pytest.mark.skipif(
    not pytest.config.getoption('--bigdata'),
    reason='need --bigdata option to run'
)

# Decorator to skip test if running under a TravisCI
not_under_travis = pytest.mark.skipif(
    "TRAVIS" in os.environ and os.environ["TRAVIS"] == "true",
    reason='Temporarily disable due to performance issues'
)


# Decorator to skip if CRDS_CONTEXT is not at lest a certain level.
def require_crds_context(required_context):
    """Ensure CRDS context is a certain level

    Parameters
    ----------
    level: int
        The minimal level required

    Returns
    -------
    pytest.mark.skipif decorator
    """
    current_context_string = crds.get_context_name('jwst')
    match = re.match('jwst_(\d\d\d\d)\.pmap', current_context_string)
    current_context = int(match.group(1))
    return pytest.mark.skipif(
        current_context < required_context,
        reason='CRDS context {} less than required context {}'.format(
            current_context_string, required_context
        )
    )
