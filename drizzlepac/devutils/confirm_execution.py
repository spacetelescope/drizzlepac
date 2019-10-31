#!/usr/bin/env python
import random
import string
import sys
def confirm_execution():
    """
    This subroutine prevents accidental execution by requiring the user to type in a randomly generated 4-character
    confirmation string. If the string is typed in incorrectly, the script will simply exit to the command line.

    :return: nothing
    """
    confirm_string=''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(4))
    foo = input("Confirm execution by entering the following randomized text string: {} \n".format(confirm_string))
    if foo != confirm_string:
        sys.exit("Execution aborted.")
    if foo == confirm_string:
        print("Execution confirmed.")
