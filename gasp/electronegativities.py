# coding: utf-8
# Copyright (c) Henniggroup.
# Distributed under the terms of the MIT License.

"""
Electronegativities module:

This module contains the electronegativities of elements, as a dictionary based
on their atomic number and co-ordination number.
"""

def get_EN(Z, CN):
    # Look up electronegativities
    #df_EN = pd.read_excel("EN_lookup.xlsx", header = None, names = range(1,13))
    #df_EN.index = np.arange(1, len(df_EN)+1)

    EN_dict = {
    1: {1: float('nan'),  2: float('nan'),  3: float('nan'),  4: float('nan'),  5: float('nan'),  6: float('nan'),  7: float('nan'),  8: float('nan'),  9: float('nan'),  10: float('nan'),  11: float('nan'),  12: float('nan')},
    2: {1: float('nan'),  2: float('nan'),  3: float('nan'),  4: float('nan'),  5: float('nan'),  6: float('nan'),  7: float('nan'),  8: float('nan'),  9: float('nan'),  10: float('nan'),  11: float('nan'),  12: float('nan')},
    3: {1: float('nan'),  2: float('nan'),  3: float('nan'),  4: float('nan'),  5: float('nan'),  6: float('nan'),  7: float('nan'),  8: float('nan'),  9: float('nan'),  10: float('nan'),  11: float('nan'),  12: float('nan')},
    4: {1: float('nan'),  2: float('nan'),  3: float('nan'),  4: float('nan'),  5: float('nan'),  6: float('nan'),  7: float('nan'),  8: float('nan'),  9: float('nan'),  10: float('nan'),  11: float('nan'),  12: float('nan')},
    5: {1: 1.641),  2: 1.641,  3: 1.641,  4: 1.641,  5: float('nan'),  6: 1.415,  7: 1.415,  8: float('nan'),  9: 1.415,  10: float('nan'),  11: float('nan'),  12: float('nan')},
    6: {1: 2.619,  2: 2.619,  3: 2.619,  4: 2.5,  5: float('nan'),  6: 2.347,  7: float('nan'),  8: float('nan'),  9: float('nan'),  10: float('nan'),  11: float('nan'),  12: float('nan')},
    7: {1: 3.565,  2: float('nan'),  3: 3.565,  4: 3.437,  5: float('nan'),  6: 3.085,  7: float('nan'),  8: float('nan'),  9: float('nan'),  10: float('nan'),  11: float('nan'),  12: float('nan')},
    8: {1: 4.561,  2: float('nan'),  3: 4.561,  4: 4.375,  5: float('nan'),  6: 3.75,  7: float('nan'),  8: float('nan'),  9: float('nan'),  10: float('nan'),  11: float('nan'),  12: float('nan')},
    9: {1: float('nan'),  2: float('nan'),  3: float('nan'),  4: float('nan'),  5: float('nan'),  6: float('nan'),  7: float('nan'),  8: float('nan'),  9: float('nan'),  10: float('nan'),  11: float('nan'),  12: float('nan')},
    10: {1: float('nan'),  2: float('nan'),  3: float('nan'),  4: float('nan'),  5: float('nan'),  6: float('nan'),  7: float('nan'),  8: float('nan'),  9: float('nan'),  10: float('nan'),  11: float('nan'),  12: float('nan')},
    11: {1: float('nan'),  2: float('nan'),  3: float('nan'),  4: float('nan'),  5: float('nan'),  6: float('nan'),  7: float('nan'),  8: float('nan'),  9: float('nan'),  10: float('nan'),  11: float('nan'),  12: float('nan')},
    12: {1: float('nan'),  2: float('nan'),  3: float('nan'),  4: float('nan'),  5: float('nan'),  6: float('nan'),  7: float('nan'),  8: float('nan'),  9: float('nan'),  10: float('nan'),  11: float('nan'),  12: float('nan')},
    13: {1: float('nan'),  2: float('nan'),  3: float('nan'),  4: float('nan'),  5: float('nan'),  6: float('nan'),  7: float('nan'),  8: float('nan'),  9: float('nan'),  10: float('nan'),  11: float('nan'),  12: float('nan')},
    14: {1: 1.645,  2: 1.645,  3: 1.645,  4: 1.645,  5: float('nan'),  6: 1.565,  7: float('nan'),  8: float('nan'),  9: float('nan'),  10: float('nan'),  11: float('nan'),  12: float('nan')},
    15: {1: 2.187,  2: float('nan'),  3: 2.187,  4: 2.187,  5: float('nan'),  6: 2.074,  7: float('nan'),  8: float('nan'),  9: float('nan'),  10: float('nan'),  11: float('nan'),  12: float('nan')},
    16: {1: float('nan'),  2: float('nan'),  3: float('nan'),  4: float('nan'),  5: float('nan'),  6: float('nan'),  7: float('nan'),  8: float('nan'),  9: float('nan'),  10: float('nan'),  11: float('nan'),  12: float('nan')},
    17: {1: float('nan'),  2: float('nan'),  3: float('nan'),  4: float('nan'),  5: float('nan'),  6: float('nan'),  7: float('nan'),  8: float('nan'),  9: float('nan'),  10: float('nan'),  11: float('nan'),  12: float('nan')},
    18: {1: float('nan'),  2: float('nan'),  3: float('nan'),  4: float('nan'),  5: float('nan'),  6: float('nan'),  7: float('nan'),  8: float('nan'),  9: float('nan'),  10: float('nan'),  11: float('nan'),  12: float('nan')},
    19: {1: float('nan'),  2: float('nan'),  3: float('nan'),  4: float('nan'),  5: float('nan'),  6: float('nan'),  7: float('nan'),  8: float('nan'),  9: float('nan'),  10: float('nan'),  11: float('nan'),  12: float('nan')},
    20: {1: float('nan'),  2: float('nan'),  3: float('nan'),  4: float('nan'),  5: float('nan'),  6: float('nan'),  7: float('nan'),  8: float('nan'),  9: float('nan'),  10: float('nan'),  11: float('nan'),  12: float('nan')},
    21: {1: float('nan'),  2: float('nan'),  3: float('nan'),  4: float('nan'),  5: float('nan'),  6: float('nan'),  7: float('nan'),  8: float('nan'),  9: float('nan'),  10: float('nan'),  11: float('nan'),  12: float('nan')},
    22: {1: float('nan'),  2: float('nan'),  3: float('nan'),  4: float('nan'),  5: float('nan'),  6: float('nan'),  7: float('nan'),  8: float('nan'),  9: float('nan'),  10: float('nan'),  11: float('nan'),  12: float('nan')},
    23: {1: float('nan'),  2: float('nan'),  3: float('nan'),  4: float('nan'),  5: float('nan'),  6: float('nan'),  7: float('nan'),  8: float('nan'),  9: float('nan'),  10: float('nan'),  11: float('nan'),  12: float('nan')},
    24: {1: 0.659,  2: float('nan'),  3: float('nan'),  4: float('nan'),  5: float('nan'),  6: float('nan'),  7: 0.659,  8: float('nan'),  9: 0.659,  10: float('nan'),  11: float('nan'),  12: 0.659},
    25: {1: float('nan'),  2: float('nan'),  3: float('nan'),  4: float('nan'),  5: float('nan'),  6: float('nan'),  7: float('nan'),  8: float('nan'),  9: float('nan'),  10: float('nan'),  11: float('nan'),  12: float('nan')},
    26: {1: float('nan'),  2: float('nan'),  3: float('nan'),  4: float('nan'),  5: float('nan'),  6: float('nan'),  7: float('nan'),  8: float('nan'),  9: float('nan'),  10: float('nan'),  11: float('nan'),  12: float('nan')},
    27: {1: float('nan'),  2: float('nan'),  3: float('nan'),  4: float('nan'),  5: float('nan'),  6: float('nan'),  7: float('nan'),  8: float('nan'),  9: float('nan'),  10: float('nan'),  11: float('nan'),  12: float('nan')},
    28: {1: float('nan'),  2: float('nan'),  3: float('nan'),  4: float('nan'),  5: float('nan'),  6: float('nan'),  7: float('nan'),  8: float('nan'),  9: float('nan'),  10: float('nan'),  11: float('nan'),  12: float('nan')},
    29: {1: float('nan'),  2: float('nan'),  3: float('nan'),  4: float('nan'),  5: float('nan'),  6: float('nan'),  7: float('nan'),  8: float('nan'),  9: float('nan'),  10: float('nan'),  11: float('nan'),  12: float('nan')},
    30: {1: float('nan'),  2: float('nan'),  3: float('nan'),  4: float('nan'),  5: float('nan'),  6: float('nan'),  7: float('nan'),  8: float('nan'),  9: float('nan'),  10: float('nan'),  11: float('nan'),  12: float('nan')},
    31: {1: float('nan'),  2: float('nan'),  3: float('nan'),  4: float('nan'),  5: float('nan'),  6: float('nan'),  7: float('nan'),  8: float('nan'),  9: float('nan'),  10: float('nan'),  11: float('nan'),  12: float('nan')},
    32: {1: float('nan'),  2: float('nan'),  3: float('nan'),  4: float('nan'),  5: float('nan'),  6: float('nan'),  7: float('nan'),  8: float('nan'),  9: float('nan'),  10: float('nan'),  11: float('nan'),  12: float('nan')},
    33: {1: 2.039,  2: float('nan'),  3: 2.039,  4: 2.039,  5: float('nan'),  6: 1.925,  7: float('nan'),  8: float('nan'),  9: float('nan'),  10: float('nan'),  11: float('nan'),  12: float('nan')},
    34: {1: float('nan'),  2: float('nan'),  3: float('nan'),  4: float('nan'),  5: float('nan'),  6: float('nan'),  7: float('nan'),  8: float('nan'),  9: float('nan'),  10: float('nan'),  11: float('nan'),  12: float('nan')},
    35: {1: float('nan'),  2: float('nan'),  3: float('nan'),  4: float('nan'),  5: float('nan'),  6: float('nan'),  7: float('nan'),  8: float('nan'),  9: float('nan'),  10: float('nan'),  11: float('nan'),  12: float('nan')},
    36: {1: float('nan'),  2: float('nan'),  3: float('nan'),  4: float('nan'),  5: float('nan'),  6: float('nan'),  7: float('nan'),  8: float('nan'),  9: float('nan'),  10: float('nan'),  11: float('nan'),  12: float('nan')},
    37: {1: float('nan'),  2: float('nan'),  3: float('nan'),  4: float('nan'),  5: float('nan'),  6: float('nan'),  7: float('nan'),  8: float('nan'),  9: float('nan'),  10: float('nan'),  11: float('nan'),  12: float('nan')},
    38: {1: float('nan'),  2: float('nan'),  3: float('nan'),  4: float('nan'),  5: float('nan'),  6: float('nan'),  7: float('nan'),  8: float('nan'),  9: float('nan'),  10: float('nan'),  11: float('nan'),  12: float('nan')},
    39: {1: float('nan'),  2: float('nan'),  3: float('nan'),  4: float('nan'),  5: float('nan'),  6: float('nan'),  7: float('nan'),  8: float('nan'),  9: float('nan'),  10: float('nan'),  11: float('nan'),  12: float('nan')},
    40: {1: float('nan'),  2: float('nan'),  3: float('nan'),  4: float('nan'),  5: float('nan'),  6: float('nan'),  7: float('nan'),  8: float('nan'),  9: float('nan'),  10: float('nan'),  11: float('nan'),  12: float('nan')},
    41: {1: float('nan'),  2: float('nan'),  3: float('nan'),  4: float('nan'),  5: float('nan'),  6: float('nan'),  7: float('nan'),  8: float('nan'),  9: float('nan'),  10: float('nan'),  11: float('nan'),  12: float('nan')},
    42: {1: 2.171,  2: float('nan'),  3: float('nan'),  4: float('nan'),  5: float('nan'),  6: float('nan'),  7: float('nan'),  8: float('nan'),  9: float('nan'),  10: 2.171,  11: float('nan'),  12: float('nan')},
    43: {1: float('nan'),  2: float('nan'),  3: float('nan'),  4: float('nan'),  5: float('nan'),  6: float('nan'),  7: float('nan'),  8: float('nan'),  9: float('nan'),  10: float('nan'),  11: float('nan'),  12: float('nan')},
    44: {1: float('nan'),  2: float('nan'),  3: float('nan'),  4: float('nan'),  5: float('nan'),  6: float('nan'),  7: float('nan'),  8: float('nan'),  9: float('nan'),  10: float('nan'),  11: float('nan'),  12: float('nan')},
    45: {1: float('nan'),  2: float('nan'),  3: float('nan'),  4: float('nan'),  5: float('nan'),  6: float('nan'),  7: float('nan'),  8: float('nan'),  9: float('nan'),  10: float('nan'),  11: float('nan'),  12: float('nan')},
    46: {1: float('nan'),  2: float('nan'),  3: float('nan'),  4: float('nan'),  5: float('nan'),  6: float('nan'),  7: float('nan'),  8: float('nan'),  9: float('nan'),  10: float('nan'),  11: float('nan'),  12: float('nan')},
    47: {1: float('nan'),  2: float('nan'),  3: float('nan'),  4: float('nan'),  5: float('nan'),  6: float('nan'),  7: float('nan'),  8: float('nan'),  9: float('nan'),  10: float('nan'),  11: float('nan'),  12: float('nan')},
    48: {1: float('nan'),  2: float('nan'),  3: float('nan'),  4: float('nan'),  5: float('nan'),  6: float('nan'),  7: float('nan'),  8: float('nan'),  9: float('nan'),  10: float('nan'),  11: float('nan'),  12: float('nan')},
    49: {1: float('nan'),  2: float('nan'),  3: float('nan'),  4: float('nan'),  5: float('nan'),  6: float('nan'),  7: float('nan'),  8: float('nan'),  9: float('nan'),  10: float('nan'),  11: float('nan'),  12: float('nan')},
    50: {1: float('nan'),  2: float('nan'),  3: float('nan'),  4: float('nan'),  5: float('nan'),  6: float('nan'),  7: float('nan'),  8: float('nan'),  9: float('nan'),  10: float('nan'),  11: float('nan'),  12: float('nan')},
    51: {1: float('nan'),  2: float('nan'),  3: float('nan'),  4: float('nan'),  5: float('nan'),  6: float('nan'),  7: float('nan'),  8: float('nan'),  9: float('nan'),  10: float('nan'),  11: float('nan'),  12: float('nan')},
    52: {1: float('nan'),  2: float('nan'),  3: float('nan'),  4: float('nan'),  5: float('nan'),  6: float('nan'),  7: float('nan'),  8: float('nan'),  9: float('nan'),  10: float('nan'),  11: float('nan'),  12: float('nan')},
    53: {1: float('nan'),  2: float('nan'),  3: float('nan'),  4: float('nan'),  5: float('nan'),  6: float('nan'),  7: float('nan'),  8: float('nan'),  9: float('nan'),  10: float('nan'),  11: float('nan'),  12: float('nan')},
    54: {1: float('nan'),  2: float('nan'),  3: float('nan'),  4: float('nan'),  5: float('nan'),  6: float('nan'),  7: float('nan'),  8: float('nan'),  9: float('nan'),  10: float('nan'),  11: float('nan'),  12: float('nan')},
    55: {1: float('nan'),  2: float('nan'),  3: float('nan'),  4: float('nan'),  5: float('nan'),  6: float('nan'),  7: float('nan'),  8: float('nan'),  9: float('nan'),  10: float('nan'),  11: float('nan'),  12: float('nan')},
    56: {1: float('nan'),  2: float('nan'),  3: float('nan'),  4: float('nan'),  5: float('nan'),  6: float('nan'),  7: float('nan'),  8: float('nan'),  9: float('nan'),  10: float('nan'),  11: float('nan'),  12: float('nan')},
    57: {1: float('nan'),  2: float('nan'),  3: float('nan'),  4: float('nan'),  5: float('nan'),  6: float('nan'),  7: float('nan'),  8: float('nan'),  9: float('nan'),  10: float('nan'),  11: float('nan'),  12: float('nan')},
    58: {1: float('nan'),  2: float('nan'),  3: float('nan'),  4: float('nan'),  5: float('nan'),  6: float('nan'),  7: float('nan'),  8: float('nan'),  9: float('nan'),  10: float('nan'),  11: float('nan'),  12: float('nan')},
    59: {1: float('nan'),  2: float('nan'),  3: float('nan'),  4: float('nan'),  5: float('nan'),  6: float('nan'),  7: float('nan'),  8: float('nan'),  9: float('nan'),  10: float('nan'),  11: float('nan'),  12: float('nan')},
    60: {1: float('nan'),  2: float('nan'),  3: float('nan'),  4: float('nan'),  5: float('nan'),  6: float('nan'),  7: float('nan'),  8: float('nan'),  9: float('nan'),  10: float('nan'),  11: float('nan'),  12: float('nan')},
    61: {1: float('nan'),  2: float('nan'),  3: float('nan'),  4: float('nan'),  5: float('nan'),  6: float('nan'),  7: float('nan'),  8: float('nan'),  9: float('nan'),  10: float('nan'),  11: float('nan'),  12: float('nan')},
    62: {1: float('nan'),  2: float('nan'),  3: float('nan'),  4: float('nan'),  5: float('nan'),  6: float('nan'),  7: float('nan'),  8: float('nan'),  9: float('nan'),  10: float('nan'),  11: float('nan'),  12: float('nan')},
    63: {1: float('nan'),  2: float('nan'),  3: float('nan'),  4: float('nan'),  5: float('nan'),  6: float('nan'),  7: float('nan'),  8: float('nan'),  9: float('nan'),  10: float('nan'),  11: float('nan'),  12: float('nan')},
    64: {1: float('nan'),  2: float('nan'),  3: float('nan'),  4: float('nan'),  5: float('nan'),  6: float('nan'),  7: float('nan'),  8: float('nan'),  9: float('nan'),  10: float('nan'),  11: float('nan'),  12: float('nan')},
    65: {1: float('nan'),  2: float('nan'),  3: float('nan'),  4: float('nan'),  5: float('nan'),  6: float('nan'),  7: float('nan'),  8: float('nan'),  9: float('nan'),  10: float('nan'),  11: float('nan'),  12: float('nan')},
    66: {1: float('nan'),  2: float('nan'),  3: float('nan'),  4: float('nan'),  5: float('nan'),  6: float('nan'),  7: float('nan'),  8: float('nan'),  9: float('nan'),  10: float('nan'),  11: float('nan'),  12: float('nan')},
    67: {1: float('nan'),  2: float('nan'),  3: float('nan'),  4: float('nan'),  5: float('nan'),  6: float('nan'),  7: float('nan'),  8: float('nan'),  9: float('nan'),  10: float('nan'),  11: float('nan'),  12: float('nan')},
    68: {1: float('nan'),  2: float('nan'),  3: float('nan'),  4: float('nan'),  5: float('nan'),  6: float('nan'),  7: float('nan'),  8: float('nan'),  9: float('nan'),  10: float('nan'),  11: float('nan'),  12: float('nan')},
    69: {1: float('nan'),  2: float('nan'),  3: float('nan'),  4: float('nan'),  5: float('nan'),  6: float('nan'),  7: float('nan'),  8: float('nan'),  9: float('nan'),  10: float('nan'),  11: float('nan'),  12: float('nan')},
    70: {1: float('nan'),  2: float('nan'),  3: float('nan'),  4: float('nan'),  5: float('nan'),  6: float('nan'),  7: float('nan'),  8: float('nan'),  9: float('nan'),  10: float('nan'),  11: float('nan'),  12: float('nan')},
    71: {1: float('nan'),  2: float('nan'),  3: float('nan'),  4: float('nan'),  5: float('nan'),  6: float('nan'),  7: float('nan'),  8: float('nan'),  9: float('nan'),  10: float('nan'),  11: float('nan'),  12: float('nan')},
    72: {1: float('nan'),  2: float('nan'),  3: float('nan'),  4: float('nan'),  5: float('nan'),  6: float('nan'),  7: float('nan'),  8: float('nan'),  9: float('nan'),  10: float('nan'),  11: float('nan'),  12: float('nan')},
    73: {1: float('nan'),  2: float('nan'),  3: float('nan'),  4: float('nan'),  5: float('nan'),  6: float('nan'),  7: float('nan'),  8: float('nan'),  9: float('nan'),  10: float('nan'),  11: float('nan'),  12: float('nan')},
    74: {1: 2.171,  2: float('nan'),  3: float('nan'),  4: float('nan'),  5: float('nan'),  6: float('nan'),  7: float('nan'),  8: float('nan'),  9: float('nan'),  10: 2.171,  11: float('nan'),  12: float('nan')},
    75: {1: 0.71,  2: float('nan'),  3: float('nan'),  4: float('nan'),  5: float('nan'),  6: float('nan'),  7: float('nan'),  8: 0.71,  9: float('nan'),  10: float('nan'),  11: float('nan'),  12: float('nan')},
    76: {1: 1.504,  2: float('nan'),  3: float('nan'),  4: float('nan'),  5: float('nan'),  6: float('nan'),  7: float('nan'),  8: 1.504,  9: float('nan'),  10: float('nan'),  11: float('nan'),  12: float('nan')}
    }

    return EN_dict[Z][CN]
