import sys
import numpy as np
import matplotlib.pyplot as plt

routine = sys.argv[1]
n = len(sys.argv) - 2

if n <= 0:
  print('Args: <routine> <indices> ...')
  sys.exit(0)
 
data = dict()
length = int(1e8)

for i in range(n):
  index = sys.argv[i + 2]
  fn = f'stat/log_{routine}_{index}.txt'
  with open(fn, 'r') as log_file:
    data[index] = dict()
    tmp = dict()
    for line in log_file.readlines():
      pid = line.split()[0]
      tim = line.split()[1]
      if pid not in tmp:
        tmp[pid] = list()
      tmp[pid].append(float(tim))
    data[index][pid] = np.array(tmp[pid])
  for pid in data[index]:
    length = min(length, data[index][pid].size)
    print(f'\033[32mRead {length} data from log file {fn} - pid={pid}\033[0m')

for index in data:
  for pid in data[index]:
    data[index][pid] = data[index][pid][data[index][pid].size - length::]
print(f'\033[31mTruncate to length of {length}\033[0m')

indices = np.arange(0, length)

plt.figure(figsize=(12, 8))
plt.title('Time Costs vs Sample')
plt.ylabel('Cost (ns)')
plt.xlabel('Sample')

colors = ['blue', 'orange', 'green', 'red', 'purple', 'brown']
ci = 0

for index in data:
  for pid in data[index]:
    plt.plot(indices, data[index][pid], label=f'Core {index} - {pid}', \
      color=colors[ci], marker='.', linestyle='')
    ci += 1

plt.legend()
plt.show()


