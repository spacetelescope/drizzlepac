#!/usr/bin/env python
"""
This script simply breaks down the flagging by bit value for a user-specified sourcelist.
"""
import sys

from astropy.table import Table
import numpy as np
bit_list = [0, 1, 2, 4, 8, 16, 32, 64, 128]
flag_meanings = ['Point Source', 'Extended Source', 'Single-Pixel Saturation', 'Multi-Pixel Saturation',
                 'Faint Magnitude Limit', 'Hot Pixel', 'Swarm Detection', 'Edge and Chip Gap',
                 'Bleeding and Cosmic Rays']

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
    out_idx_list = np.zeros(9, dtype=int)
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

if __name__ == "__main__":
    sl_name = sys.argv[1]
    table_data = Table.read(sl_name, format='ascii.ecsv')
    flag_counts = np.zeros(9, dtype=int)
    for table_line in table_data:
        flagval = table_line["Flags"]
        flag_bits = deconstruct_flag((flagval))
        flag_counts += flag_bits
    print("Total Sourcelist Length{}{}".format("."*20, len(table_data["Flags"])))
    for ctr in range(0, len(bit_list)):
        bit_val = bit_list[ctr]
        padding1 = 6 - len(str(bit_val))
        padding2 = 27 - len(flag_meanings[ctr])
        text_stuff = "Bit value {}{}{}{}{}".format(bit_val, "."*padding1, flag_meanings[ctr], padding2*".", flag_counts[ctr])
        print(text_stuff)