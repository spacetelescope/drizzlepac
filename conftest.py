"""Project default for pytest"""


# Utilities to support command line passing of arguments to test_randomlist.py
def pytest_addoption(parser):
    """Support options for the pytest test_randomlist.py."""
    parser.addoption("--start_row", action="store", default=0)
    parser.addoption("--num_rows", action="store", default=1)
    parser.addoption("--master_list", action="store", default="ACSWFC3ListDefault50.csv")
    parser.addoption("--svm_list", action="store", default="svm_input.lst")
