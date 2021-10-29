#!/usr/bin/python

import argparse

# Parse of command line options and provides help to user
parser = argparse.ArgumentParser(description='This code is intended to extract data from Illumina NGS runs and to export it to a .csv file.', prefix_chars='-+')
parser.add_argument("-p", help="the path to the Illumina run folder (no '/' at the end of the path).")
parser.add_argument("-o", help="the path to the output file.")

# Actually parse command line options
args = parser.parse_args()
MyInPath = args.p
MyOutFile = args.o

run_folder = MyInPath

print run_folder

from interop import py_interop_run_metrics, py_interop_run, py_interop_table
import numpy
import pandas as pd

run_metrics = py_interop_run_metrics.run_metrics()

valid_to_load = py_interop_run.uchar_vector(py_interop_run.MetricCount, 0)
py_interop_table.list_imaging_table_metrics_to_load(valid_to_load)

run_metrics.read(run_folder, valid_to_load)

columns = py_interop_table.imaging_column_vector()
py_interop_table.create_imaging_table_columns(run_metrics, columns)

headers = []
for i in range(columns.size()):
    column = columns[i]
    if column.has_children():
        headers.extend([column.name()+"("+subname+")" for subname in column.subcolumns()])
    else:
        headers.append(column.name())


column_count=7

print headers[:column_count]

column_count = py_interop_table.count_table_columns(columns)
row_offsets = py_interop_table.map_id_offset()
py_interop_table.count_table_rows(run_metrics, row_offsets)

row_count=len(row_offsets)

data = numpy.zeros((row_offsets.size(), column_count), dtype=numpy.float32)
py_interop_table.populate_imaging_table_data(run_metrics, columns, row_offsets, data.ravel())

d = []
for col, label in enumerate(headers[:column_count]):
    d.append( (label, pd.Series([val for val in data[:row_count, col]], index=[tuple(r) for r in  data[:row_count, :3]])))

df = pd.DataFrame.from_dict(dict(d))
#print(df.to_string(index=False)) # print to string format (better for screen or text editor visualization)
print(df.to_csv(index=False))    # print to .csv format ( better for spreadsheet based visualization)
# Write results to the output file
myoutfile=open(MyOutFile, "wt")
myoutfile.write(df.to_csv(index=False)) 

