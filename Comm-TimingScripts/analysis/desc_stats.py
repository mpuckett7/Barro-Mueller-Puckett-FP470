import statistics

for nodes in [1, 2, 4, 8]:
    for tpn in [1, 2, 4]:
        print(f"-------------------\nOpening raw-{nodes}-{tpn}.txt")
        f = open(f"raw-{nodes}-{tpn}.txt", 'r')
        data = []
        for line in f:
            data.append(float(line))
        print(f"Average: {statistics.mean(data)}")
        print(f"Std. Deviation: {statistics.stdev(data)}")
        print(f"Total Time: {sum(data)}")
        f.close()
        