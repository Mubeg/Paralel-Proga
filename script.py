from tqdm import tqdm as tqdm
import subprocess

filename = "output.txt"
filename2 = "better_output.txt"

#res = ["n", "m", "k", "h", "tau", "Error", "Timing"]

K = [2**(i) for i in range(8, 15)]
M = [2**(i) for i in range(8, 15)]

arr = []

for n in range(1, 5):
        for m in M:
            for k in K:
                arr.append([n, m, k]) 

with open(filename2, "w") as out:
    out.write("n, m, k, h, tau, Error, Timing\n")
    for n, m, k in tqdm(arr):
                h = -1
                tau = -1
                err = -1
                timing = -1

                subprocess.run(["mpirun", "-np", "{}".format(n), "lab1_2", "{}".format(m), "{}".format(k)])
                with open(filename, "r") as f:
                    last_line = f.readlines()[-1].split(',')
                    h = float(last_line[0])
                    tau = float(last_line[1])
                    err = float(last_line[2])
                timing = 0
                for _ in range(0, 10):
                    subprocess.run(["mpirun", "-np", "{}".format(n), "lab1_2_timing", "{}".format(m), "{}".format(k)])
                    with open(filename, "r") as f:
                        last_line = f.readlines()[-1].split(',')
                        timing += int(last_line[3])
                subprocess.run(["rm", "-f", filename])
                timing /= 10
                out.write("{}, {}, {}, {}, {}, {}, {}\n".format(n, m, k, h, tau, err, timing))
                #res.append([n, m, k, h, tau, err, timing])
        

