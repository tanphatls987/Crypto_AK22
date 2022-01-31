import matplotlib.pyplot as plt
import csv
import numpy as np
import sys 

def gaussian(x, mu, sig2):
    return 1./(2.*np.sqrt(2.*np.pi)*np.sqrt(sig2))*np.exp(-np.power(x - mu, 2.) / (2. * sig2))

def make_plot(rows):
    n = int(rows[0][0])
    k = int(rows[0][1])
    h = int(rows[0][2])
    l = int(rows[0][3])
    times = int(rows[0][4])
    seed = int(rows[0][5])
    #print(n, k, h, l, times, seed)

    rho = n // k

    x_values = np.linspace(0, rho, times)
    mu0 = float(rows[2][0])
    sig20 = float(rows[2][1])
    mu1 = float(rows[4][0])
    sig21 = float(rows[4][1])

    freq0 = [float(x) for x in rows[1]]
    freq1 = [float(x) for x in rows[3]]
    plt.figure()
    plt.plot(freq0, label='0 bit', color='green')
    plt.plot(freq1, label='1 bit', color='red')

    #plt.plot(x_values, gaussian(x_values, mu0, sig20), label='0 bit Gauss')
    #plt.plot(x_values, gaussian(x_values, mu1, sig21), label='1 bit Gauss')

    plt.axvline(x = rho // 2, color='r', linestyle='dotted')
    #plt.legend()
    plt.xlabel("#1s")
    plt.ylabel("Frequency")
    plt.title("N = " + str(n) + ", K = " + str(k) + ", H = " + str(h) + ", L = " + str(l) + ", Times = " + str(times) + ", Seed = " + str(seed))
    #plt.show()
    figname = "figure_" + str(n) + "_" + str(k) + "_" + str(h) + "_" + str(l) + ".png"
    plt.savefig(figname)


if (len(sys.argv) != 2):
    sys.exit('Required 1 argument as name of csv file')

file_name = sys.argv[1]
with open(file_name, newline='\n') as file:
    rows = list(csv.reader(file, delimiter=','))
    total = len(rows)
    for x in range(0, total, 5):
        make_plot(rows[x:x + 5])
