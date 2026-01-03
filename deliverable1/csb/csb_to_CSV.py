import re
import csv

INPUT = "report_CSB.txt"
OUTPUT = "report_CSB.csv"

pattern = re.compile(
    r"matrix considered: (.+)\n"
    r"the mode used was: (.+)\n"
    r"(?:number of threads used: (\d+) \| schedule: ([a-zA-Z]+) \| chunk: (\d+)\n)?"
    r"average time took: ([\d\.]+) ms",
    re.MULTILINE
)

rows = []
with open(INPUT, "r") as f:
    text = f.read()

for m in pattern.finditer(text):
    matrix = m.group(1)
    mode = m.group(2)

    threads = m.group(3) if m.group(3) else ""
    schedule = m.group(4) if m.group(4) else ""
    chunk = m.group(5) if m.group(5) else ""
    time = m.group(6)

    rows.append([matrix, mode, threads, schedule, chunk, time])

with open(OUTPUT, "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["matrix", "mode", "threads", "schedule", "chunk", "time_ms"])
    writer.writerows(rows)

print("Generated:", OUTPUT)
