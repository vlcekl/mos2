#!/usr/bin/python
#
# File name:   img_ana.py
# Date:        2016/04/19 23:32
# Author:      Lukas Vlcek
#
# Description: 
#

import sys
import string as s
import re
import numpy as np

if __name__ == "__main__":

    f = open(sys.argv[1], 'r')

    for line in f:
        sarr = re.findall('\S+', line)

    f.close()

# end of img_ana.py 
