import csv
import locale

locale.setlocale(locale.LC_ALL, '')

data = []
with open("Africaneo_dataset_geographical_coordinates.csv", newline="") as csvfile:
    r = csv.reader(csvfile, delimiter=";")
    for row in r:
        data.append(row)

seen = []
for d in data[1:]:
    d[2] = locale.atof(d[2])
    d[3] = locale.atof(d[3])
    if d[2:] in seen:
        if [d[2] + 1, d[3]] in seen:
            if [d[2], d[3] + 1] in seen:
                if [d[2] + 1, d[3] + 1] in seen:
                    d[2:] = [d[2] - 1, d[3]]
                else:
                    d[2:] = [d[2] + 1, d[3] + 1]
            else:
                d[2:] = [d[2], d[3] + 1]
        else:
            d[2:] = [d[2] + 1, d[3]]
    seen.append(d[2:])

with open("Africaneo_dataset_geographical_coordinates_c.csv", "w", newline="") as csvfile:
    w = csv.writer(csvfile)
    for d in data:
        w.writerow(d)
