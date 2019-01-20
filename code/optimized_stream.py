# -*- coding: utf-8 -*-
""" Script for generating an optimized scan pattern, which is 1 um by 1 um.
"""

import sys


n_pixels = 10
XMAX = 1000
XMIN = 0
dwell = 0.1

pixel = list()
time = list()

spacing = (XMAX - XMIN) / n_pixels

for n in range(n_pixels):
    pixel.append(XMIN + n*spacing)
    time.append(dwell)

text_file = open("stream_scan.txt", "w")
for n in range(n_pixels):
    text_file.write("%f    %f\n" %(pixel[n], time[n]))

text_file.close()
print("done")