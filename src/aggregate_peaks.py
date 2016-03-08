#!/usr/bin/env python

import os
import sys
import argparse
import xml.etree.ElementTree as ET

# recommended use: aggregate_peaks.py -p data/nmr -s '_' > data/processed/peaks.tsv

def aggregate_peaks(path, split_on, dest):
    subdirs = next(os.walk(path))[1] # get list of only immediate subdirectories
    for subdir in subdirs:
        if split_on is not None:
            filename_data = subdir.split(args.split_on)
        else:
            filename_data = []
        path_to_peaks = os.path.join(args.path, subdir, "1", "pdata", "1", "peaklist.xml")
        if os.path.exists(path_to_peaks):
            tree = ET.parse(path_to_peaks)
            peaks = tree.findall(".//Peak1D")
            for peak in peaks:
                chemical_shift = peak.attrib.get("F1")
                intensity = peak.attrib.get("intensity")
                dest.write('\t'.join([subdir] + filename_data + [chemical_shift] + [intensity])+'\n')
        else:
            sys.stderr.write('Missing file:' + path_to_peaks + '\n')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Extract a text file of NMR peaks from a directory of spectra')
    parser.add_argument('-p','--path', dest='path',
                       type=str, help='Path to folder full of folders, one for each NMR tube')
    parser.add_argument('-s', '--split_on', dest='split_on',
                       type=str, help='Character(s) to split filename on')
    parser.add_argument('-o', '--out', nargs='?', type=argparse.FileType('w'),
                       default=sys.stdout)
    args = parser.parse_args()
    aggregate_peaks(path=args.path,split_on=args.split_on,dest=args.out)

