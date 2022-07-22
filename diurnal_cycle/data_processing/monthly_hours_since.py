"""Prints number of hours of month-start since a fixed reference time.

Usage:
    $ python monthly_hours_since.py ref_time first_month last_month

Three command line arguments needed, all formatted as strings
that can be passed to numpy datetime64 calls:
    ref_time -- reference time from which other times are calculated as
                hours since. E.g., 1980-01-02T00:00
    first_month -- the first month for which "hours since" is needed.
                   E.g., 2002-01
    last_month -- the last month for which "hours since" is needed.
                  The returned list of "hours since" includes one more
                  data point after this month, representing the first 
                  hour of the following month. E.g., 2005-12

Returns a list of hours since ref_time for the first hour in each month 
from first_month to (last_month + 1 month), inclusive.
"""

import sys
import numpy as np

ref_time = np.datetime64(sys.argv[1])

st_time = np.datetime64(sys.argv[2])
end_time = np.datetime64(sys.argv[3]) + np.timedelta64(2, 'M')

time_array = np.arange(st_time, end_time, dtype='datetime64[M]')

for t in time_array.astype('datetime64[s]'):
    delta_t_hours = ((t - ref_time)/np.timedelta64(1,'h')).astype(np.int64)
    print(delta_t_hours)
