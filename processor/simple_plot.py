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


    fpath_default  = "histos/plotsTopEFT.pkl.gz"
    hin1 = get_hist_from_pkl(fpath_default)
    
    save_dir = "../plots"
    for x in hin1["jet_pt"].axes():
        print(x)
    h1 = hin1["jet_pt"]
    
    for x in hin1["electron_pt"].axes():
        print(x)
    h2 = hin1["electron_pt"]
    
    for x in hin1["muon_pt"].axes():
        print(x)
    h3 = hin1["muon_pt"]
    
    for x in hin1["jet_eta"].axes():
        print(x)
    h4 = hin1["jet_eta"]
    
    for x in hin1["jet_phi"].axes():
        print(x)
    h5 = hin1["jet_phi"]
    
    for x in hin1["jet_high_pt"].axes():
        print(x)
    h6 = hin1["jet_high_pt"]

    fig1, ax1 = plt.subplots(1, 1, figsize=(7,7))
    hist.plot1d(h1)
    ax1.autoscale(axis='y')
    fig1.savefig(os.path.join(save_dir,'xjet_pt_histogram'))
    ax1.clear()
    
    fig2, ax2 = plt.subplots(1, 1, figsize=(7,7))
    hist.plot1d(h2)
    ax2.autoscale(axis='y')
    fig2.savefig(os.path.join(save_dir,'electron_pt_histogram'))
    ax2.clear()
    
    fig3, ax3 = plt.subplots(1, 1, figsize=(7,7))
    hist.plot1d(h3)
    ax3.autoscale(axis='y')
    fig3.savefig(os.path.join(save_dir,'muon_pt_histogram'))
    ax3.clear()
    
    fig4, ax4 = plt.subplots(1, 1, figsize=(7,7))
    hist.plot1d(h4)
    ax4.autoscale(axis='y')
    fig4.savefig(os.path.join(save_dir,'xjet_eta_histogram'))
    ax4.clear()
    
    fig5, ax5 = plt.subplots(1, 1, figsize=(7,7))
    hist.plot1d(h5)
    ax5.autoscale(axis='y')
    fig5.savefig(os.path.join(save_dir,'xjet_phi_histogram'))
    ax5.clear()
    
    fig6, ax6 = plt.subplots(1, 1, figsize=(7,7))
    hist.plot1d(h6)
    ax6.autoscale(axis='y')
    fig6.savefig(os.path.join(save_dir,'xjet_high_pt_histogram'))
    ax6.clear()

main()
