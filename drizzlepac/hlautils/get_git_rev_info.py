#!/bin/env python
"""
Displays revision information on the specified git-controlled repository

Path
----
hlapipeline/hlapipeline/utils/get_git_rev_info.py

Dependencies
------------
None.

Inputs
------
* Required input
    1: *code_path*
        * Path to query for git revision information.

Example
-------
Get git revision info on git repository "foo"::

    hlapipeline/hlapipeline/utils/get_git_rev_info.py foo
"""
import os
import sys
from stsci.tools import logutil
import traceback

__taskname__ = 'get_git_rev_info'

log = logutil.create_logger(__name__, level=logutil.logging.NOTSET, stream=sys.stdout)

# -------------------------------------------------------------------------------------------------
def print_rev_id(local_repo_path):
    """prints information about the specified local repository to STDOUT. Expected method of execution: command-line or
    shell script call

    Parameters
    ----------
    local_repo_path: string
        Local repository path.

    Returns
    =======
    Nothing as such. subroutine will exit with a state of 0 if everything ran OK, and a value of '111' if
    something went wrong.
    """
    start_path = os.getcwd()
    try:
        log.info("Local repository path: {}".format(local_repo_path))
        os.chdir(local_repo_path)

        log.info("\n===== Remote URL INFO =====")
        instream = os.popen('git remote -v')
        for streamline in instream:
            log.info("{}".format(streamline.strip()))

        log.info("\n===== Local Branches =====")
        instream = os.popen('git branch')
        for streamline in instream:
            log.info("{}".format(streamline.strip()))

        log.info("\n===== Most Recent Commit =====")
        instream = os.popen('git log |head -1')
        for streamline in instream:
            log.info("{}\n".format(streamline.strip()))

        rv = 0

    except Exception:
        exc_type, exc_value, exc_tb = sys.exc_info()
        traceback.print_exception(exc_type, exc_value, exc_tb, file=sys.stdout)
        rv = 111
        log.info("WARNING! get_git_rev_info.print_rev_id() encountered a problem and cannot continue.")
    finally:
        os.chdir(start_path)
        if rv != 0:
            sys.exit(rv)
# -------------------------------------------------------------------------------------------------
def get_rev_id(local_repo_path):
    """returns the current full git revision id of the specified local repository. Expected method of execution: python
    subroutine call

    Parameters
    ----------
    local_repo_path: string
        Local repository path.

    Returns
    =======
    full git revision ID of the specified repository if everything ran OK, and "FAILURE" if something went
    wrong.
    """
    start_path = os.getcwd()
    try:
        os.chdir(local_repo_path)

        instream = os.popen("git --no-pager log --max-count=1 | head -1")
        for streamline in instream.readlines():
            streamline = streamline.strip()
            if streamline.startswith("commit "):
                rv = streamline.replace("commit ", "")
            else:
                raise
    except Exception:
        rv = "FAILURE: git revision info not found"
    finally:
        os.chdir(start_path)

    return(rv)
# -------------------------------------------------------------------------------------------------
if(__name__ == '__main__'):
    local_repo_path = sys.argv[1]
    print_rev_id(local_repo_path)
