'''Test-Script to see the difference of sputtering with an without redeposition

reads the yamamura.srf (without redeposition) and the yamamura_redep.srf (with redeposition)
and plots those files with solid and dashed lines, respectivly.
This file is not needed for the miniTopSim itself!
'''
import numpy as np
import matplotlib.pyplot as plt

data = []
data_new = []
time_step_data = []
i = -1

filenames = ["yamamura.srf", "yamamura_redep.srf"]
with open(filenames[0], "r") as file:
    for line in file:
        if "surface" in line:
            data.append(time_step_data)
            time_step_data = np.ndarray(shape=(0,2))
            continue
        tmp = [float(i) for i in line.split(" ")]
        print(tmp)
        time_step_data = np.vstack((time_step_data, [tmp]))
with open(filenames[1], "r") as file:
    for line in file:
        if "surface" in line:
            data_new.append(time_step_data)
            time_step_data = np.ndarray(shape=(0,2))
            continue
        tmp = [float(i) for i in line.split(" ")]
        print(tmp)
        time_step_data = np.vstack((time_step_data, [tmp]))

data = data[1:]
data_new = data_new[1:]
data1 = np.array(data)
data2 = np.array(data_new)

colors = ["blue", "orange", "red", "brown", "green", "purple", "gray", "black", "magenta", "yellow", "red"]
i = 0
for d in data1:
    plt.plot(d[:,0], d[:,1])
    i += 1

i = 0
for d in data2:
    plt.plot(d[:,0], d[:,1], "--")
    i += 1

plt.show()
