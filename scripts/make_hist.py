from argparse import ArgumentParser
import os
import sys
from matplotlib import pyplot as plt
import pandas as pd
import numpy as np

def main(input, out_folder, value='value'):
    if not os.path.exists(out_folder):
        os.makedirs(out_folder)
            
    for i, line in enumerate(input):
        line = line.strip()
        if i == 0:
            samples = [s for s in line.split('\t')]
            sample_hists = np.zeros((1000, len(samples)))
            continue
        
        for i, value in enumerate(line.split('\t')):
            if value == '.':
                value = 0
            if value > sample_hists.shape[0]:
                sample_hists = np.vstack((sample_hists, np.zeros((100, sample_hists.shape[1]))))
            
            sample_hists[value, i] += 1
    
    i = sample_hists.shape[0]
    while np.sum(sample_hists[i, :]) == 0:
        i -= 1
    sample_hists = sample_hists[:i+1, :]
    max_value = i
                        
    for i, sample in enumerate(samples):
        plt.bar(np.arange(max_value), sample_hists[i])
        plt.title(sample)
        plt.xlabel(value)
        plt.ylabel('Count')
        plt.savefig(os.path.join(out_folder, sample + '.png'))
        plt.close()

    pd.DataFrame(sample_hists, columns=samples).to_csv(os.path.join(out_folder, 'raw_data.csv'), index=True)

if __name__ == '__main__':
    parser = ArgumentParser()
    
    parser.add_argument('out_folder')
    parser.add_argument('--value', default='value')
    
    args = parser.parse_args()
    
    main(sys.stdin, args.out_folder, args.value)