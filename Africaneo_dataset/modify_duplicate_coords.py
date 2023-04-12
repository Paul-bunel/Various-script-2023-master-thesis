import csv
import locale

locale.setlocale(locale.LC_ALL, '')

data = []
with open("Africaneo_dataset_simplified.csv", newline="") as csvfile:
    r = csv.reader(csvfile, delimiter=";")
    for row in r:
        data.append(row)

seen = []
for d in data[1:]:
    d[3] = locale.atof(d[3])
    d[4] = locale.atof(d[4])
    if d[3:] in seen:
        if [d[3] + 1, d[4]] in seen:
            if [d[3], d[4] + 1] in seen:
                if [d[3] + 1, d[4] + 1] in seen:
                    d[3:] = [d[3] - 1, d[4]]
                else:
                    d[3:] = [d[3] + 1, d[4] + 1]
            else:
                d[3:] = [d[3], d[4] + 1]
        else:
            d[3:] = [d[3] + 1, d[4]]
    seen.append(d[3:])

with open("Africaneo_dataset_simplified_c.csv", "w", newline="") as csvfile:
    w = csv.writer(csvfile)
    for d in data:
        w.writerow(d)
