#!/usr/bin/env python3

"""Try to remove the smallest possible number of individuals from list of related individuals
    to keep only unrelated individuals
"""

import sys
import datetime
import numpy as np


def print_with_date(msg):
    """Print message with date"""

    now_text = datetime.datetime.now().isoformat(sep=' ', timespec='milliseconds')
    print("{} {}".format(now_text, msg))


def main(grm_file, rel_file, remove_file):
    """Entry point if called as an executable"""

    print_with_date("loading relatives table...")
    rel_table = np.loadtxt(rel_file, dtype=np.int32, usecols=(0, 1))
    sub_array, num_rel = np.unique(rel_table, return_counts=True)
    print_with_date("done")

    sub_array2 = np.arange(max(sub_array) + 1, dtype=np.int32)
    num_rel2 = np.zeros(sub_array2.size, dtype=np.int32)
    num_rel2[sub_array] = num_rel

    print_with_date("listing neighbors per subject...")
    sub_array_col1, num_rel_col1 = np.unique(rel_table[:, 0], return_counts=True)
    sub_array_col2, num_rel_col2 = np.unique(rel_table[:, 1], return_counts=True)
    num_rel2_col = np.zeros((sub_array2.size, 2), dtype=np.int32)
    num_rel2_col[sub_array_col1, 0] = num_rel_col1
    num_rel2_col[sub_array_col2, 1] = num_rel_col2
    num_cum_col = np.cumsum(num_rel2_col, axis=0)
    rel_table_sorted = np.column_stack((rel_table[np.argsort(rel_table[:, 0]), 1],
                                        rel_table[np.argsort(rel_table[:, 1]), 0]))
    neighbors = {}
    for sub in sub_array2:
        neighbors[sub] = np.append(rel_table_sorted[num_cum_col[sub-1, 0]:num_cum_col[sub, 0], 0],
                                   rel_table_sorted[num_cum_col[sub-1, 1]:num_cum_col[sub, 1], 1])
    print_with_date("done")

    def exclude_subject(subject):
        num_rel2[neighbors[subject]] -= 1
        num_rel2[subject] = 0
        excluded_subjects.append(subject)

    excluded_subjects = []
    subject_codes = np.loadtxt(grm_file + ".grm.id", dtype=str)
    # remove subjects with negative value as id
    negative_codes = np.where(np.char.startswith(subject_codes[:, 0], '-'))[0] + 1
    for sub in negative_codes:
        exclude_subject(sub)

    print_with_date("excluding related subjects...")
    counter = 0
    while True:
        counter += 1
        if counter % 100000 == 0:
            print(counter)
        sing = np.where(num_rel2 == 1)[0]
        if sing.size > 0:
            ex_sub = neighbors[sing[0]][np.where(num_rel2[neighbors[sing[0]]] > 0)[0][0]]
        else:
            ex_sub = sub_array2[np.argmax(num_rel2)]
            if num_rel2[ex_sub] == 0:
                break
        exclude_subject(ex_sub)
    print_with_date("done")

    excluded_subjects = np.array(excluded_subjects)
    excluded_codes = subject_codes[excluded_subjects - 1, :]
    np.savetxt(remove_file, excluded_codes, fmt="%s")


if __name__ == '__main__':
    main(*sys.argv[1:4])
