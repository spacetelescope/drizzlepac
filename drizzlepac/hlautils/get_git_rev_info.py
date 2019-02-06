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
import os, sys
#-----------------------------------------------------------------------------------------------------------------------
def print_rev_id(localRepoPath):
	"""prints information about the specified local repository to STDOUT. Expected method of execution: command-line or
	shell script call

	Parameters
    ----------
	localRepoPath: string
		Local repository path.

	Returns
    =======
    Nothing as such. subroutine will exit with a state of 0 if everything ran OK, and a value of '111' if
	something went wrong.
	"""
	start_path = os.getcwd()
	try:
		print("Local repository path: {}".format(localRepoPath))
		os.chdir(localRepoPath)
		print("\n== Remote URL")
		os.system('git remote -v')
                
		# print("\n== Remote Branches")
		# os.system("git branch -r")

		print("\n== Local Branches")
		os.system("git branch")

		print("\n== Most Recent Commit")
		os.system("git log |head -1")
		rv = 0
	except:
		rv = 111
		print("WARNING! get_git_rev_info.print_rev_id() encountered a problem and cannot continue.")
	finally:
		os.chdir(start_path)
		if rv != 0:
			sys.exit(rv)
#-----------------------------------------------------------------------------------------------------------------------
def get_rev_id(localRepoPath):
	"""returns the current full git revision id of the specified local repository. Expected method of execution: python
    subroutine call

	Parameters
    ----------
	localRepoPath: string
		Local repository path.

	Returns
    =======
    full git revision ID of the specified repository if everything ran OK, and "FAILURE" if something went
    wrong.
    """
	start_path = os.getcwd()
	try:
		os.chdir(localRepoPath)

		instream = os.popen("git --no-pager log --max-count=1 | head -1")
		for streamline in instream.readlines():
			streamline = streamline.strip()
			if streamline.startswith("commit "):
				rv = streamline.replace("commit ","")
			else:
				raise
	except:
		rv = "FAILURE: git revision info not found"
	finally:
		os.chdir(start_path)

	return(rv)
#-----------------------------------------------------------------------------------------------------------------------
if(__name__ == '__main__'):
	localRepoPath = sys.argv[1]
	print_rev_id(localRepoPath)