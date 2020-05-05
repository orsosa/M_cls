#!/bin/bash

./get_Multi.py -i DH_M_MxCut/RGA/pippim/zxQ2_bins_200504/M_results.root -o RGAresult.root -b binning_info.txt -s DH_M_MxCut/SIM/pippim/zxQ2_bins_200504/M_results.root -e DH_M_elec/RGA/zxQ2_bins_elec_200417/M_elec.root -d DH_M_elec/SIM/zxQ2_bins_elec_200417/M_elec.root -a

./get_Multi.py -i DH_M_MxCut/RGB/pippim/zxQ2_bins_200504/M_results.root -o RGBresult.root -b binning_info.txt -s DH_M_MxCut/SIM/pippim/zxQ2_bins_200504/M_results.root -e DH_M_elec/RGB/zxQ2_bins_elec_200417/M_elec.root -d DH_M_elec/SIM/zxQ2_bins_elec_200417/M_elec.root -a
