import crds
import os
import pytest
import re

__all__ = ['runslow', 'slow', 'require_bigdata',
           'require_crds_context']


# Decorator to indicate slow tests
runslow = pytest.mark.skipif(
    not pytest.config.getoption("--runslow"),
    reason="need --runslow option to run"
)

# pytest marker to mark resource-intensive tests that should not be
# executed with every commit.
slow = runslow

# Decorator to indicate BIGDATA required
require_bigdata = pytest.mark.skipif(
    not pytest.config.getoption('--bigdata'),
    reason='need --bigdata option to run'
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
