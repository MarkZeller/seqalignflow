#!/usr/bin/env python3

import argparse
import seaborn as sns
import numpy as np
from matplotlib import (pyplot as plt,lines)

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", type=str, required=True, help='samtools depth file')
parser.add_argument("-o", "--output", type=str, required=True, help='output file')
parser.add_argument("-s", "--sample_id", type=str, required=True, help='graph title')
parser.add_argument("-d", "--depth_cut_off", type=str, required=True, help='depth threshold')
parser.add_argument("-n", "--normalize",  action='store_true', help='will normalize depth')
args = parser.parse_args()


def parse_depth(depth_input):
    lines_in_file = open(depth_input, 'r').readlines()
    genome_size = len(lines_in_file) + 1
    depth = [0] * genome_size
    references = set()
        
    with open(depth_input) as depth_object:
        for row in depth_object:
            genome_id, position, depth_count = row.split()
 
            references.add(genome_id)
 
            if len(references) > 1:
                raise Exception(' This script only handles one genome - contig.')
 
            depth[int(position)] = int(depth_count)
        return depth

def plot_depth(depth_report, output_name, plot_title, depth_cut_off, normalize):
    data = parse_depth(depth_report)
 
    y_label = "Normalized Depth" if normalize else "Depth"
    data = [xx / max(data) for xx in data] if normalize else data
    genome_size = (len(data)) + 1
    depth_cut_off = int(depth_cut_off)

 
    sns.set(color_codes=True)
    plt.title(plot_title)
    ax = plt.subplot(111)
 
    sns_plot = sns.lineplot(x=range(len(data)), y=data)
    sns_plot.set(xlabel='Genome Position (bp)', ylabel=y_label)
 
    if not normalize:
        ax.add_line(lines.Line2D([0, genome_size + 1], [depth_cut_off], color="r"))
 
    plt.savefig(output_name, bbox_inches='tight', dpi=600)
    plt.close()


if __name__ == "__main__":
    plot_depth(args.input, args.output, args.sample_id, args.depth_cut_off, args.normalize)
