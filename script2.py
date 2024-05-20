from tqdm import tqdm as tqdm
import subprocess

arr = [i + 1 for i in range(8)]

print(arr)

subprocess.run(["rm", "-f", "output.txt"])
for item in tqdm(arr):
    for _ in range(0, 10):
        subprocess.run(["./integral_script", "0.2", "{}".format(item)])

        

