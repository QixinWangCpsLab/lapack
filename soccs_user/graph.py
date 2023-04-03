import sys
import numpy as np
import matplotlib.pyplot as plt

routine = sys.argv[1]
n = len(sys.argv) - 3
lt = int(sys.argv[2])

if n <= 0:
  print('Args: <routine> <indices> ...')
  sys.exit(0)
 
data = dict()
length = int(1e8)

for i in range(n):
  index = sys.argv[i + 3]
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
    for pid in tmp:
      data[index][pid] = np.array(tmp[pid])
  for pid in data[index]:
    dip = data[index][pid]
    length = min(length, dip.size)
    print(f'\n{dip.size} from log file {fn} - pid={pid}')
    print(f'\033[32m\tAvg: {np.mean(dip):<9.1f}\tStd: {np.std(dip):<9.1f}')
    print(f'\tMin: {np.min(dip):<9.1f}\tMax: {np.max(dip):<9.1f}\033[0m')

lt = min(lt, length)
for index in data:
  for pid in data[index]:
    data[index][pid] = data[index][pid][data[index][pid].size - lt::]
print(f'\033[31mTruncate to length of {lt}\033[0m')

print('\n\033[31mSample Statistics -------------\033[0m')
for i in range(n):
  index = sys.argv[i + 3]
  for pid in data[index]:
    dip = data[index][pid]
    print(f'\n{dip.size} from - pid={pid}')
    print(f'\033[32m\tAvg: {np.mean(dip):<9.1f}\tStd: {np.std(dip):<9.1f}')
    print(f'\tMin: {np.min(dip):<9.1f}\tMax: {np.max(dip):<9.1f}\033[0m')

indices = np.arange(0, lt)

plt.figure(figsize=(12, 8))
plt.title('Time Costs vs Sample')
plt.ylabel('Cost (ns)')
plt.xlabel('Sample')

colors = ['blue', 'orange', 'green', 'red', 'purple', 'brown', 'cyan', 'violet', 'black', 'magenta']
ci = 0

for index in data:
  for pid in data[index]:
    if index == '-1':
      mk, ms, ls, lw = '.', 1, '', 0
    else:
      mk, ms, ls, lw = 'D', 0.5, '-', 0.25
    plt.plot(indices, data[index][pid], label=f'Core {index} - {pid}', \
      color=colors[ci], marker=mk, markersize=ms, linestyle=ls, linewidth=lw)
    ci += 1
plt.semilogy()
plt.legend()
plt.show()


