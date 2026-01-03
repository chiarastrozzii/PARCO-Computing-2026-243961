import re
import csv

input_file = "report_CSR.txt"
output_file = "report_CSR.csv"

pattern_seq = re.compile(
    r"matrix considered: (.+?)\n"
    r"the mode used was: (SEQ)\n"
    r"average time took: ([0-9.]+) ms"
)

pattern_par = re.compile(
    r"matrix considered: (.+?)\n"
    r"the mode used was: (PAR)\n"
    r"number of threads used: (\d+) \| type of schedule: (\w+) \| chunck size considered: (\d+)\n"
    r"average time took: ([0-9.]+) ms"
)

pattern_simd = re.compile(
    r"matrix considered: (.+?)\n"
    r"the mode used was: (SIMD)\n"
    r"average time took: ([0-9.]+) ms"
)

rows = []

with open(input_file, "r") as f:
    data = f.read()

for m in pattern_seq.finditer(data):
    rows.append([m.group(1), "SEQ", 1, "none", 0, m.group(3)])

for m in pattern_par.finditer(data):
    rows.append([
        m.group(1), "PAR", m.group(3), m.group(4), m.group(5), m.group(6)
    ])

for m in pattern_simd.finditer(data):
    rows.append([m.group(1), "SIMD", 1, "none", 0, m.group(3)])

with open(output_file, "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["matrix", "mode", "threads", "schedule", "chunk", "time_ms"])
    writer.writerows(rows)

print("CSR CSV generated:", output_file)
