import gzip
import pickle
import json
import os
import uproot
import matplotlib.pyplot as plt
import numpy as np
#import coffea
#from coffea import hist, processor
#from coffea.hist import plot
import os, sys
import time

import datetime

#from topcoffea.plotter.make_html import make_html
import topcoffea
#import topcoffea.modules.HistEFT

from coffea import hist
from topcoffea.modules.HistEFT import HistEFT


def get_hist_from_pkl(path_to_pkl):
    h = pickle.load( gzip.open(path_to_pkl) )
    return h

def main():

    #axis_names = {"jets": ["a","b"]}

    fpath_default  = "histos/plotsTopEFT.pkl.gz"
    hin1 = get_hist_from_pkl(fpath_default)
    save_dir = "../plots"
    for x in hin1["jets"].axes():
        print(x)
    h = hin1["jets"]
   

    fig, ax = plt.subplots(1, 1, figsize=(7,7))
    hist.plot1d(h)
    ax.autoscale(axis='y')
    fig.savefig(os.path.join(save_dir,'xjet_pt_eta_histogram'))
    ax.clear()
    #make_html(save_dir_full)

main()
