"""Function to extract random entries (lines) from a CSV file 

The function, randomSelectFromCSV, allows the user to specify the desired
number of entries to obtain from a comma-separated values (CSV) file
which contains data extracted from the STScI archive for ACS, WFC3, and 
WFPC2 instruments. The data are comprised of the following information
as identified by the database name: observationID, trgName, trgposRA, 
trgPosDec, detector, aperture, filters, enrBandpassName, timMin, timMax, 
timExposure,and asnID.  The entries from the input table are chosen
at random with no duplication, and returned as an Astropy table.
"""

from random import seed
from random import sample

from astropy.table import Table

__all__ =['randomSelectFromCSV']

def randomSelectFromCSV(tableName, numEntries, seedValue):
    """Function to extract random entries (lines) from a CSV file 
    
    Parameters
    ==========
    filename : str
        Filename of the input master CSV file containing individual 
        images or association names, as well as observational 
        information regarding the images

    numEntries : int
        Number of entries/rows to extract from the master input CSV file

    seedValue : int
        Value used to initialize the random number generator for the
        selection of random entries

    Returns
    =======
    outputTable : object
        Astropy Table object
    """

    # Initialize the random number generator
    seed(seedValue)

    # Get the contents of the table
    dataTable = Table.read(tableName, format='ascii.csv')
    numRows   = len(dataTable)

    # Generate a sequence of integers the size of the table, and then 
    # obtain a random subset of the sequence with no duplicate selections
    sequence = list(range(numRows))
    subset   = sample(sequence, numEntries)

    # Extract the subset rows...
    outputTable = dataTable[subset]
    #outputTable = dataTable[0:numEntries]

    # Returns the outputTable which is an Astropy Table object
    return(outputTable)
