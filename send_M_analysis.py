#!/usr/bin/env python
import glob
import os
from sys import argv
import subprocess
############## global setting ##########
OPTINIT=""
WF="DH_M"
ACCOPT = ["rhoCut","noMxCut","excl","DZCut","noPi0Cut"]
#######################################
for opt in argv[1:]:
    if opt not in ACCOPT:
        print "options accepted: {}".format(ACCOPT)
        print "by default, only MxCut is used."
        exit (1)
    OPTINIT += opt + " "
if "noMxCut" not in OPTINIT and "excl" not in OPTINIT:
    WF += "_MxCut"
if OPTINIT != "":
    WF += "_" + reduce(lambda x,y:x+"_"+y,OPTINIT.split())

#OPTINIT=""
#OPTINIT="rhoCut"
#OPTINIT="noMxCut"
#OPTINIT="noMxCut rhoCut"
#OPTINIT="excl"
##################################
outdir="/work/clas12/users/osoto/" + WF
app="/home/osoto/dihadron/M_cls/run_M"
script="/home/osoto/dihadron/M_cls/run_Mapp.sh"
BINFO={"pippim":"/home/osoto/dihadron/M_cls/binning_info.txt","pipi0":"/home/osoto/dihadron/M_cls/binning_info_pi0.txt"}

REGEXP = ["/volatile/clas12/users/osoto/data/mix/RGB/pass0/v16/102/*_pippim_all.root", 
"/volatile/clas12/users/osoto/data/mix/RGB/pass0/v16/102/*_pipaa_all.root",
"/volatile/clas12/users/osoto/data/mix/RGB/pass0/v16/102/*_pimaa_all.root",
"/volatile/clas12/users/osoto/data/mix/RGB/pass0/v16/102/*_pippip_all.root",
"/volatile/clas12/users/osoto/data/mix/RGB/pass0/v16/102/*_pimpim_*.root",
"/volatile/clas12/users/osoto/data/mix/RGB/pass0/v16/106/*_pippim_all.root",
"/volatile/clas12/users/osoto/data/mix/RGB/pass0/v16/106/*_pipaa_all.root",
"/volatile/clas12/users/osoto/data/mix/RGB/pass0/v16/106/*_pimaa_all.root",
"/volatile/clas12/users/osoto/data/mix/RGB/pass0/v16/106/*_pippip_all.root",
"/volatile/clas12/users/osoto/data/mix/RGB/pass0/v16/106/*_pimpim_*.root",
"/volatile/clas12/users/osoto/data/mix/RGA/trains/v16_v2/*_pippim_all.root",
"/volatile/clas12/users/osoto/data/mix/RGA/trains/v16_v2/*_pipaa_all.root",
"/volatile/clas12/users/osoto/data/mix/RGA/trains/v16_v2/*_pimaa_all.root",
"/volatile/clas12/users/osoto/data/mix/RGA/trains/v16_v2/*_pippip_all.root",
"/volatile/clas12/users/osoto/data/mix/RGA/trains/v16_v2/*_pimpim_*.root",
"/volatile/clas12/users/osoto/data/mix/RGA/trains/pass1/v1_4/*_pippim_all.root",
"/volatile/clas12/users/osoto/data/mix/RGA/trains/pass1/v1_4/*_pipaa_all.root",
"/volatile/clas12/users/osoto/data/mix/RGA/trains/pass1/v1_4/*_pimaa_all.root",
"/volatile/clas12/users/osoto/data/mix/RGA/trains/pass1/v1_4/*_pippip_all.root",
"/volatile/clas12/users/osoto/data/mix/RGA/trains/pass1/v1_4/*_pimpim_*.root",
"/volatile/clas12/users/osoto/data/mix/RGB/pass0/v16/*_pippim_all.root",
"/volatile/clas12/users/osoto/data/mix/RGB/pass0/v16/*_pipaa_all.root",
"/volatile/clas12/users/osoto/data/mix/RGB/pass0/v16/*_pimaa_all.root",
"/volatile/clas12/users/osoto/data/mix/RGB/pass0/v16/*_pippip_all.root",
"/volatile/clas12/users/osoto/data/mix/RGB/pass0/v16/*_pimpim_*.root",
"/volatile/clas12/users/osoto/sim/mix/*_pippim_*.root",
"/volatile/clas12/users/osoto/sim/mix/*_pipaa_*.root",
"/volatile/clas12/users/osoto/sim/mix/*_pimaa_*.root",
"/volatile/clas12/users/osoto/sim/mix/*_pippip_*.root",
"/volatile/clas12/users/osoto/sim/mix/*_pimpim_*.root"]

DEBUG=False
TRACK="analysis"
####### check dir ####
def checkdir(path):
    if not os.path.isdir(path):
        print "The directory " + path + " does not exist"
        print "creating it"
        cmd = "mkdir -p " + path
        subprocess.call(cmd,shell=True)
        return 0

####### get file list in path given the regexp ####
def get_file_list(path,regexp="*.root"):
    if DEBUG : print path 
    flist=glob.glob(path + "/" + regexp)
    return flist

###### import file list from file ############
def import_flist(fname):
    infile=open(fname)
    flist=[]
    for line in infile:
        flist.append(line.strip())
    return flist

##### add job to workflow
def add_job(wf,c=0):
    global outdir, app, script, REGEXP, BINFO
    size=0
    indir=reduce(lambda x,y: x +"/"+ y,REGEXP[c].split('/')[:-1])
    fnl = get_file_list(indir,REGEXP[c].split('/')[-1])
    if DEBUG : print fnl
    for fname in fnl:
        size += os.path.getsize(fname)/1024./1024.*1.15
    jname = wf + "_" + str(c)
    cmd  = "swif add-job -workflow " + wf + " -time 15h -cores 1 -ram 1500mb -project clas12 -track " + TRACK +" -disk " + "{0:.0f}".format(size) + "mb "
    cmd += " -name " + jname
    cmd += " -input " + app.split("/")[-1] + " file:" + app
    cmd += " -input " + script.split("/")[-1] + " file:" + script
    for fname in fnl:
        if ("noPi0Cut" not in OPTINIT):
            cmd += " -input " + fname.split("/")[-1] + " file:" + fname
        elif ("aa_" in fname) :
            cmd += " -input " + fname.split("/")[-1] + " file:" + fname

    outdir_local=outdir

    opt=OPTINIT
    if "RGB" in REGEXP[c]:
        opt += " rgb"
        if "102" in REGEXP[c]:
            outdir_local=outdir_local+"/RGB102"
        elif "106" in REGEXP[c]:
            outdir_local=outdir_local+"/RGB106"
        else :
            outdir_local=outdir_local+"/RGB"

    elif "pass1/v1_4" in REGEXP[c]:
        outdir_local=outdir_local+"/RGAOB"
        opt += " rgaob"
    elif "sim" in REGEXP[c]:
        opt = "sim"
        outdir_local=outdir_local+"/SIM"
    else:
        outdir_local=outdir_local+"/RGA"

    if "aa" in REGEXP[c]:
        cmd += " -input " + BINFO['pipi0'].split("/")[-1] + " file:" + BINFO['pipi0']
        if "pip" in REGEXP[c]:
            outdir_local=outdir_local+"/pippi0"
        elif "pim" in REGEXP[c]:
            outdir_local=outdir_local+"/pimpi0"
            opt += " rho-"
        opt+=" binfo:"+BINFO['pipi0'].split("/")[-1]
    else:
        outdir_local=outdir_local+"/pippim"
        cmd += " -input " + BINFO['pippim'].split("/")[-1] + " file:" + BINFO['pippim']
        opt+=" binfo:"+BINFO['pippim'].split("/")[-1]

    checkdir(outdir_local)
    if DEBUG : print "opt in python: "+ opt
    cmd += " ./" + script.split("/")[-1] + " " + outdir_local + " " + opt

    if DEBUG : print (cmd)
    subprocess.call(cmd,shell=True)

######### create workflow #####
def create_wf():
    global WF
    cmd="swif create " + WF
    if DEBUG : print (cmd)
    subprocess.call(cmd,shell=True)

#### run workflow ############
def run_wf():
    global WF
    cmd="swif run -workflow " + WF
    subprocess.check_call(cmd.split(" "))

######## main routine ####

def main():
    global WF
    create_wf()
    c=0
    for c in range(len(REGEXP)):
        add_job(WF,c)
        c = c+1
if __name__=="__main__":
    main()
    run_wf()

