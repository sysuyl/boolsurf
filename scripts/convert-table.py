import sys


def main():
    csv_filename = sys.argv[1]
    latex_table = []

    with open(csv_filename) as stats:
        lines = stats.readlines()

        for i in range(1, len(lines)):
            record_list = lines[i].strip().split(", ")
            record = [record_list[0].replace(
                "/", "\\").split("\\")[-1].split(".")[0]]

            timings = [round(float(time), 3) for time in record_list[12:]]
            perc = [round((time/timings[-1]) * 100, 0) for time in timings]

            timings_str = [
                f'{timings[i]} ({perc[i]}\%)' for i in range(len(timings))]

            record += [record_list[1]] + [record_list[3]] + record_list[5:10] + \
                [str(timings[-1])] + timings_str[:2]

            record = "& ".join(record)
            latex_table.append(record)

    for line in latex_table:
        print(f'{line}\\\\')


if __name__ == '__main__':
    main()
