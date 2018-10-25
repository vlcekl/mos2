#!/usr/bin/python
#
# File name:   cfg2xyz.py
# Date:        2016/05/24 08:46
# Author:      Lukas Vlcek
#
# Description: 
#

import sys
import string as s
import re
import numpy as np
from scipy.constants import physical_constants as phc

if __name__ == "__main__":

    f = open(sys.argv[1], 'r')
#    data = np.array([map(float, re.findall('\S+', ln)[:]) for ln in f])

    line = f.readline()
    sys.stdout.write(line)
    line = f.readline()
    sys.stdout.write(line)

    while 1:
        line = f.readline()
        if not line:
            break
        sarr = re.findall('\S+', line)
        i = int(sarr[0])
        if i == 1:
            a = 'Ca'
        elif i == 2:
            a = 'N'
        elif i == 3:
            a = 'O'
        elif i == 6:
            a = 'Om'
        elif i == 4:
            a = 'Mn'
        elif i == 5:
            a = 'Mn'

        print a, sarr[1], sarr[2], sarr[3]


    f.close()

# end of cfg2xyz.py 
