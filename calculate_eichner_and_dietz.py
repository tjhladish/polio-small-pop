#!/usr/bin/python
from sys import argv

circ_ivl_filename = argv[1] # pass filename with circulation intervals as command line argument
circ_intervals = []
intercase_intervals = []
num_replicates = 0
required_pcases = 0
ed_stat_output_filename = 'e_and_d_values-2pcase_filter-alt.out'

for replicate in file(circ_ivl_filename):
    times = map(float, replicate.strip().split(','))
    intervals = [times[i] - times[i-1] for i in range(1,len(times))]
    intervals = intervals[required_pcases:]
    circ_intervals = circ_intervals + intervals
    intercase_intervals = intercase_intervals + intervals[:-1]
    if len(intervals) > 0:
        num_replicates += 1

circ_intervals.sort()
intercase_intervals.sort()

num_ci = len(circ_intervals)
num_ii = len(intercase_intervals)

fo = open(ed_stat_output_filename, 'w')
last_icase_idx = 0
icase_idx = 0

for i in range(num_ci):
    numerator = float(num_ci - i)
    interval_duration = circ_intervals[i]
    for icase_idx in range(last_icase_idx, len(intercase_intervals)):
        if intercase_intervals[icase_idx] >= interval_duration:
            break
    last_icase_idx = icase_idx
    denominator = num_replicates + num_ii - icase_idx
    fo.write(str(interval_duration) + ' ' + str(numerator/denominator) + '\n')

fo.close()
