#!/usr/bin/env python
import sys
import pdb
from astropy.table import Table
import numpy

# ======================================================================================================================


def deconstruct_flag(flagval):
    """Breaks down an integer flag value into individual component bit values.

    Parameters
    ----------
    flagval : int
        Flag value to deconstruct

    Returns
    -------
    out_idx_list : list
        a 9-element numpy array of 0s and 1s. Each element of the array represents the presence of a particular
        bit value (element 0 = bit 0, element 1 = bit 1, ..., element 3 = bit 4 and so on...)
    """
    bitlist = [1, 2, 4, 8, 16, 32, 64, 128]
    flagval = int(flagval)
    # out_bit_list = []
    out_idx_list = numpy.zeros(9, dtype=int)
    if flagval == 0:
        # out_bit_list = [0]
        out_idx_list[0] = 1
    if flagval > 0:
        idx = 1
        for bit in bitlist:
            if flagval & bit > 0:
                # out_bit_list.append(bit)
                out_idx_list[idx] = 1
            if bit > flagval:
                break
            idx += 1
    return out_idx_list


# ======================================================================================================================

def display_catalog_bit_populations(flag_data):
    """Breaks all input flag values down into their constituent bit values and displays a bit-by-bit population summary

    Parameters
    ----------
    flag_data : astropy.table.column.Column object
        'Flags' column of a given sourcelist to analyze

    Returns
    -------
    Nothing.
    """
    bit_list = [0, 1, 2, 4, 8, 16, 32, 64, 128]
    flag_meanings = ['Point Source', 'Extended Source', 'Single-Pixel Saturation', 'Multi-Pixel Saturation',
                     'Faint Magnitude Limit', 'Hot Pixel', 'Swarm Detection', 'Edge and Chip Gap',
                     'Bleeding and Cosmic Rays']
    flag_counts = numpy.zeros(9, dtype=int)
    n_sources = len(flag_data)
    for flagval in flag_data:
        flag_counts += deconstruct_flag(flagval)
    max_length = 5
    for bitval in flag_counts:
        max_length = max([max_length, len(str(bitval))])
    print("{}".format("-"*60))
    print("{}FLAG BIT VALUE POPULATION SUMMARY".format(" "*13))
    print("Bit   Meaning{}Count Percentage".format(" "*20))
    fill_char = " "
    for ctr in range(0, len(bit_list)):
        bit_val = bit_list[ctr]
        pct_val = 100.0*(float(flag_counts[ctr])/float(n_sources))
        padding1 = 6 - len(str(bit_val))
        padding2 = 27 - len(flag_meanings[ctr])
        padding3 = max_length-len(str(flag_counts[ctr]))
        if pct_val == 100.:
            padding4 = 3
        elif pct_val >= 10.:
            padding4 = 4
        else:
            padding4 = 5
        print("{}{}{}{}{}{}{}{:.3f}%".format(bit_val, fill_char*padding1, flag_meanings[ctr], padding2*fill_char,
                                             fill_char*padding3, flag_counts[ctr], fill_char*padding4, pct_val))
    print("{}".format(" -- " * 15))
    print("NOTE: As the flag value for a given source can be composed ")
    print("of multiple bits, the above percentage values need not add")
    print("up to 100%.")
    print("{}\n".format("-" * 60))

# ======================================================================================================================
if __name__ == "__main__":
    table_data = Table.read(sys.argv[1], format='ascii.ecsv')
    display_catalog_bit_populations(table_data['Flags'])