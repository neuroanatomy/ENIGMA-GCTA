#!/usr/bin/env python3

"""Extract related pairs of individuals from GRM
"""

import sys
import math
from pathlib import Path
from functools import partial
import numpy as np


def file_float_iterator(path):
    """given a path, return an iterator over the file
    that lazily loads the file
    """
    path = Path(path)
    with path.open('rb') as file:
        reader = partial(file.read1, 2**20)
        file_iterator = iter(reader, bytes())
        for chunk in file_iterator:
            yield np.frombuffer(chunk, np.float32)


def coord_from_index(index):
    """Get coordinates in the GRM from index"""

    row = math.floor(math.sqrt(index*2 - 1) - 0.5) + 1
    triangular_number = row * (row - 1) // 2
    column = index - triangular_number
    return (row, column)


def main(grm_file, cut_off, out_file):
    """Entry point if called as an executable"""
    
    cut_off = float(cut_off)
    with open(out_file, 'w') as file_obj:
        counter = 1
        for rel_array in file_float_iterator(grm_file + '.grm.bin'):
            indexes = np.where(rel_array >= cut_off)[0] + counter
            for index in indexes:
                ind1, ind2 = coord_from_index(index)
                if ind1 != ind2:
                    file_obj.write("\t".join(map(str, [ind1, ind2, rel_array[index-counter]])) +  '\n')
            counter += rel_array.size


if __name__ == '__main__':
    main(*sys.argv[1:4])


