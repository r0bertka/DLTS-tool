import pandas as pd
import os
import seaborn as sns
import matplotlib.pyplot as plt


def find_start_end_of_monotonous_section(data):
    current_min = 1000
    min_index = None
    current_max = 0
    max_index = None

    for index, temp in data['T'].items():
        if temp < current_min:
            current_min = temp
            min_index = index
        elif temp > current_max:
            current_max = temp
            max_index = index

    return min_index, max_index


col_index = range(17)
names = ['T', 'Np', 'Cp', 'Cb','G',
         'GS2_W1', 'GS2_W2', 'GS2_W3', 'GS2_W4', 'GS2_W5', 'GS2_W6',
         'GS4_W1', 'GS4_W2', 'GS4_W3', 'GS4_W4', 'GS4_W5', 'GS4_W6']

path = os.path.join('..', 'data', 'He6_#2_DLTS_after300C')
data = pd.read_csv(path, delimiter='    ', skiprows=2, header=None, names=names, usecols=col_index)

min_index, max_index = find_start_end_of_monotonous_section(data)

data = data[min_index:max_index]



# sns.lineplot(x='T', y='G', data=data, ci=None)
data.plot(x='T', y=names[5:11])
plt.show()

