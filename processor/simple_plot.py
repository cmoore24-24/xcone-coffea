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

    for x in hin1['electron_eta'].axes():
        print(x)
    h7 = hin1['electron_eta']

    for x in hin1['electron_phi'].axes():
        print(x)
    h8 = hin1['electron_phi']

    for x in hin1['muon_eta'].axes():
        print(x)
    h9 = hin1['muon_eta']

    for x in hin1['muon_phi'].axes():
        print(x)
    h10 = hin1['muon_phi']

    for x in hin1['jet_counts'].axes():
        print(x)
    h11 = hin1['jet_counts']

    for x in hin1['missing_pt'].axes():
        print(x)
    h12 = hin1['missing_pt']

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

    fig7, ax7 = plt.subplots(1, 1, figsize=(7,7))
    hist.plot1d(h7)
    ax7.autoscale(axis='y')
    fig7.savefig(os.path.join(save_dir,'electron_eta_histogram'))
    ax7.clear()

    fig8, ax8 = plt.subplots(1, 1, figsize=(7,7))
    hist.plot1d(h8)
    ax8.autoscale(axis='y')
    fig8.savefig(os.path.join(save_dir,'electron_phi_histogram'))
    ax8.clear()

    fig9, ax9 = plt.subplots(1, 1, figsize=(7,7))
    hist.plot1d(h9)
    ax9.autoscale(axis='y')
    fig9.savefig(os.path.join(save_dir,'muon_eta_histogram'))
    ax9.clear()

    fig10, ax10 = plt.subplots(1, 1, figsize=(7,7))
    hist.plot1d(h10)
    ax10.autoscale(axis='y')
    fig10.savefig(os.path.join(save_dir,'muon_phi_histogram'))
    ax10.clear()

    fig11, ax11 = plt.subplots(1, 1, figsize=(7,7))
    hist.plot1d(h11)
    ax11.autoscale(axis='y')
    fig11.savefig(os.path.join(save_dir,'jet_counts_histogram'))
    ax11.clear()

    fig12, ax12 = plt.subplots(1, 1, figsize=(7,7))
    hist.plot1d(h12)
    ax12.autoscale(axis='y')
    fig12.savefig(os.path.join(save_dir,'missing_pt_histogram'))
    ax12.clear()

main()
