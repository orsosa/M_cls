#!/usr/bin/env python
from optparse import OptionParser

parser = OptionParser()

parser.add_option("-o","--output-file", action="store", type="string",dest="outfile", help="output file", metavar="<FILE-OUT>")
parser.add_option("-i","--input-file", action="store", type="string",dest="infile", help="input file with M_survey results", metavar="<FILE-IN>")
parser.add_option("-s","--sim-file", action="store", type="string",dest="simfile", help="input file with M_survey results for simulations", metavar="<SIM-FILE>")
parser.add_option("-b","--binning-file", action="store", type="string",dest="binfo", help="input file with binning info", metavar="<BINFO>")
parser.add_option("-e","--electron-file", action="store", type="string",dest="efile", help="input file with M_elec result", metavar="<EFILE>")
parser.add_option("-d","--sim-electron-file", action="store", type="string",dest="esimfile", help="input file with M_elec result for simulations", metavar="<ESIMFILE>")

parser.add_option("-a","--not-use-acceptance-in-file", action="store_false", dest="acceptance", help="don't use the acceptance in sim file, make it and save it in the simfile", default=True)


(opt, args) = parser.parse_args() 

from ROOT import TFile, TObject, TCanvas, gROOT, gStyle, TH1D, TH2D, TProfile, TLine, kDashed, kRed, kGray, TF1, TGaxis, kBlue, kRed, kBlack, kMagenta, kCyan, AddressOf, TPad, TPaveText, TLatex, TLegend, TLegendEntry
from sys import argv
import os
import subprocess
from array import array 
import copy

gROOT.Macro(os.path.expanduser( '~/rootlogon.C' ) )

gStyle.SetOptStat(0)
gStyle.SetOptFit(0)

gStyle.Print()

def CanvasPartition(C, Nx, Ny, lMargin=0.15, rMargin=0.05, bMargin=0.15, tMargin=0.05):
    if (not C):
        return(0)
    #Setup Pad layout:
    vSpacing = 0.0;
    vStep  = (1.- bMargin - tMargin - (Ny-1) * vSpacing) / Ny
    hSpacing = 0.0
    hStep = (1.- lMargin - rMargin - (Nx-1) * hSpacing) / Nx
    vposd = 0
    vposu = 0
    vmard = 0
    vmaru = 0
    vfactor = 0
    hposl = 0
    hposr = 0
    hmarl = 0
    hmarr = 0
    hfactor = 0
    for i in range (Nx) :
        if (i==0): 
            hposl = 0.0
            hposr = lMargin + hStep
            hfactor = hposr-hposl
            hmarl = lMargin / hfactor
            hmarr = 0.0
        elif (i == Nx-1) :
            hposl = hposr + hSpacing
            hposr = hposl + hStep + rMargin
            hfactor = hposr-hposl
            hmarl = 0.0
            hmarr = rMargin / (hposr-hposl)
        else: 
            hposl = hposr + hSpacing
            hposr = hposl + hStep
            hfactor = hposr-hposl
            hmarl = 0.0
            hmarr = 0.0
         
            for j in range(Ny):
                if (j==0): 
                    vposd = 0.0
                    vposu = bMargin + vStep
                    vfactor = vposu-vposd
                    vmard = bMargin / vfactor
                    vmaru = 0.0
                elif (j == Ny-1):
                    vposd = vposu + vSpacing
                    vposu = vposd + vStep + tMargin
                    vfactor = vposu-vposd
                    vmard = 0.0
                    vmaru = tMargin / (vposu-vposd)
                else:
                    vposd = vposu + vSpacing
                    vposu = vposd + vStep
                    vfactor = vposu-vposd
                    vmard = 0.0
                    vmaru = 0.0

                C.cd(0)
                name = "pad_"+ str(i) + "_" + str(j)
                #pad = gROOT.FindObject(name)
                pad = TPad(name,"",hposl,vposd,hposr,vposu)
                pad.SetLeftMargin(hmarl)
                pad.SetRightMargin(hmarr)
                pad.SetBottomMargin(vmard)
                pad.SetTopMargin(vmaru)
                pad.SetFrameBorderMode(0)
                pad.SetBorderMode(0)
                pad.SetBorderSize(0)
                pad.Draw()

def store_sim(simfile):
    for obj in simfile.Get("M").GetListOfKeys():
        hrs = obj.ReadObj()
        bcode = obj.GetName().split("_")[-2] + "_" + obj.GetName().split("_")[-1]
        ofile.Get("rs/M").cd()
        hrs.Write("",TObject.kOverwrite)
        hmc = simfile.Get("M_mc/hM_mc_" + bcode)
        ofile.Get("mc/M").cd()
        hmc.Write("",TObject.kOverwrite)


def get_acceptance(simfile):
    d = ""
    if not simfile.Get("M_ac"):
        d = simfile.mkdir("M_ac")
    else:
        d = simfile.Get("M_ac") 
    d.cd()    
    for obj in simfile.Get("M").GetListOfKeys():
        hrs = obj.ReadObj()
        bcode = obj.GetName().split("_")[-2] + "_" + obj.GetName().split("_")[-1]
        hmc = simfile.Get("M_mc/hM_mc_" + bcode)
        hacc = hmc.Clone("hM_ac_" + bcode)
        hacc.Divide(hrs,hmc,1,1,"b")
        hacc.Write("",TObject.kOverwrite)
        
def plot_dir(dlist,binfo,outfile,prefix,infix,color=[kBlack],legend=[],norm=False):
    c = TCanvas()
    c.SetGrid(0)
    bnarr = binfo.keys()
    bnarr.sort()
    bncode = reduce(lambda x,y: x + y, bnarr)
    namex = bnarr[-2]
    namey = bnarr[-1]
    xarr = binfo[namex]
    yarr = binfo[namey]

    c.GetPad(0).SetMargin(0.2,0.15,0.2,0.15)
    c.Divide(len(xarr)-1,len(yarr)-1,0.0,0.0)
    #CanvasPartition(c,len(xarr)-1,len(yarr)-1)
    Ntot = 1
    for k in range(len(bnarr)-2):
        Ntot *= len(binfo[bnarr[k]]) - 1
        
    for k in range(Ntot):
        indarr = []
        fact = 1
        for n in range(len(bnarr)-2):
            edg = binfo[bnarr[n]]
            indarr.append( (k/fact) % (len(edg)-1) )
            fact *= len(binfo[bnarr[n]]) - 1

            
        for j in range(len(yarr)-1,0,-1):
            YMAX = -100000000.
            YMIN = 100000000.
            row = j
            for i in range(len(xarr)-1):
                col = i+1
                bcode = reduce(lambda x,y: str(x) + str(y), indarr)
                bcode = str(bcode) + str(col-1) + str(row-1)
                cnt = 0
                for d in dlist:
                    h = d.Get(prefix[cnt] + bncode + "_" + bcode)
                    maxb = h.GetMaximumBin()
                    minb = h.GetMinimumBin()
                    ymaxtest = h.GetBinContent(maxb)
                    ymintest = h.GetBinContent(minb)
#                    if True:
#                        print "before::" + d.GetPath() + "/" + h.GetName() + "::row::" + str(row) + ":: ymax / integral / bin / content(bin): " + str(ymaxtest) + " / " + str(h.Integral()) + " / " + str( maxb) +  " / " + str( h.GetBinContent(maxb))
                    if norm and h.Integral()!=0:
                        ymaxtest = ymaxtest/h.Integral()
                        ymintest = ymintest/h.Integral()
                    if norm and h.Integral()==0:
                        ymaxtest = 0
                        ymintest = 0
                    if ymaxtest > YMAX:
                        YMAX = ymaxtest
                    if ymintest < YMIN:
                        YMIN = ymintest
                    cnt +=1
 #               if norm:
 #                   print bncode + "::row::" + str(row) + "::ymax: " + str(YMAX)

            for i in range(len(xarr)-1):
                col = i+1
                
                pad = c.cd(col + (len(yarr)-1 - row)*(len(xarr)-1))
                pad.SetGridy(0)
                #pad.SetGridx(0)
                pad.SetTicks(0,0)
                bcode = reduce(lambda x,y: str(x) + str(y), indarr)
                bcode = str(bcode) + str(col-1) + str(row-1)
                cnt = 0
                leg = TLegend(0.77,0.64,0.95,0.89)
                leg.SetTextSize(0.15)
         
                for d in dlist:
                    h = d.Get(prefix[cnt] + bncode + "_" + bcode)
                    h.GetYaxis().SetRangeUser(YMIN,YMAX*1.1)
                    if "ac" in prefix[cnt] or "Acceptance" in d.GetPath(): h.GetYaxis().SetRangeUser(0,1)
                    h.GetYaxis().SetLabelSize(0.15)
                    h.GetXaxis().SetLabelSize(0.1)
                    h.GetYaxis().SetLabelSize(0.15)
                    xtitle = str(binfo[namex][col-1]) + '<'+namex+'<'+str(binfo[namex][col]) + "     " + prefix[cnt].split('_')[0][1:]
                    h.GetXaxis().SetTitle(xtitle)
                    h.GetXaxis().CenterTitle()
                    h.GetXaxis().SetTitleSize(0.12)
                    ytitle = str(binfo[namey][row-1])+'<'+namey+'<'+str(binfo[namey][row])
                    h.GetXaxis().SetNdivisions(505)
                    h.GetYaxis().SetNdivisions(505)
                    h.SetLineWidth(3)
                    h.SetLineColor(color[cnt])
                    h.SetMarkerColor(color[cnt])
                    if (legend != []):
                        le = leg.AddEntry(h,legend[cnt],'lp')
                        le.SetTextColor(color[cnt])
                    if cnt==0:
                        if not norm: h.DrawCopy()
                        elif (h.Integral()!=0):
                            hnorm = h.DrawNormalized("e")
                            hnorm.GetYaxis().SetRangeUser(YMIN,YMAX*1.1)
                        else: h.DrawCopy()
                        if 'hM_' in  h.GetName():
                            rhopeak = TLine (0.77,0,0.77,YMAX)
                            rhopeak.SetLineStyle(kDashed)
                            rhopeak.SetLineWidth(1)
                            rhopeak.SetLineColor(kRed)
                            rhopeak.DrawLine(0.77,0,0.77,YMAX)
                        pad.Update()
                        pt = pad.FindObject("title")
                        pt.SetX1NDC(0.15)
                        pt.SetTextSizePixels(12)
                        if col == len(xarr)-1:
                            yaxis = TGaxis(pad.GetUxmax(),pad.GetUymin(),pad.GetUxmax(),pad.GetUymax(),0,1,1,"U+-")
                            yaxis.SetTitle(ytitle+"     ")
                            yaxis.SetTitleSize(0.12)
                            yaxis.SetTitleOffset(-0.25)
                            yaxis.SetName(h.GetName()+"_ax_y")
                            yaxis.Draw("")
                            pad.GetListOfPrimitives().Add(yaxis)
                            pad.Update()
                    else:
                        if not norm: h.DrawCopy("same")
                        elif (h.Integral()!=0):
                            hnorm = h.DrawNormalized("samee")
                            hnorm.GetYaxis().SetRangeUser(YMIN,YMAX*1.1)
                        else: h.DrawCopy("samee")
                    cnt += 1
                if (legend != []):
                    leg.Draw()
                    pad.GetListOfPrimitives().Add(leg)
        outfile.cd()
        bcode = reduce(lambda x,y: str(x) + str(y), indarr)
        cncode = reduce(lambda x,y: x + y, bnarr[-3])
        c.Write("can_" + infix + "_" + str(cncode) + "_" + str(bcode),TObject.kOverwrite)

def get_depend(ofile,d,binfo,bn,hnpref="hM_",odir="Multiplicities"):
    ofile.mkdir(odir+'/'+bn)
    edg = binfo[bn]
    bnarr = binfo.keys()
    bnarr.sort()
    bncode = reduce(lambda x,y: x + y, bnarr)
    bnd = dict( zip( bnarr, range(len(bnarr)) ) )
    medg = []
    hm = d.Get(hnpref+bncode+"_000")
    for i in range (hm.GetNbinsX()):
        medg.append(hm.GetBinLowEdge(i+1))

    medg.append(hm.GetXaxis().GetBinUpEdge(hm.GetNbinsX()))

    binfoM = copy.deepcopy(binfo)
    binfoM.pop(bn)
    binfoM["m"] = medg
    hd = {}
    Nh = 1
    for v in binfoM.values():
        Nh *= len(v)-1;

    bnarrM = binfoM.keys()
    bnarrM.sort()
    bncodeM = reduce(lambda x,y: x + y, bnarrM)
    
    for k in range(Nh):
        fact   = 1
        bcodeM = ""
        bcut   = ''
        bid    = {}
        for n in bnarrM:
            ind = (k/fact) % (len(binfoM[n])-1)
            bcodeM += str(ind)
            bcut += '({0:.2f}<'.format(binfoM[n][ind])+n+'<{0:.2f})'.format(binfoM[n][ind+1])
            fact *= len(binfoM[n])-1
            if n != "m":
                bid[n] = ind
            else:
                bim = ind
        
        ofile.cd(odir + '/' + bn)
        h = TH1D("h"+bn+'_'+bncodeM+"_"+bcodeM,bn+": "+bcut,len(edg)-1,array('d',edg))
        
        for b in range(h.GetNbinsX()):
            bcode = [-1 for k in range(len(bnarrM))]
            for n,v in bid.items():
                bcode[bnd[n]]= bid[n]
            bcode[bnd[bn]] = b
            bcode = reduce(lambda x,y: str(x) + str(y),bcode)
            hname = hnpref + bncode + "_" + bcode
            hm = d.Get(hname)
            h.SetBinContent(b+1,hm.GetBinContent(bim+1))
            h.SetBinError(b+1,hm.GetBinError(bim+1))
        h.Write("",TObject.kOverwrite)
    return binfoM
        
    
        
def get_Q2x_bin(h,bnarr):
    bcode = h.GetName().split("_")[-1]
    bnd = dict( zip( bnarr, range(len(bnarr)-1,-1,-1) ) )
    q2bin = (int(bcode) / pow(10,(bnd["Q2"]))) % 10
    xbin = (int(bcode) / pow(10,(bnd["x"]))) % 10
    return (q2bin+1,xbin+1)

def get_bin_factor(h,binfo):
    bcode = int(h.GetName().split("_")[-1])
    bnarr = binfo.keys()
    bnarr.sort()
    fact = 1
    for k in range(len(bnarr)):
        i =  (bcode / pow(10,len(bnarr) -1 -k) ) % 10
        fact *= binfo[bnarr[k]][i+1]- binfo[bnarr[k]][i]
    return fact

def scaleHist(h):
    for k in range(h.GetNbinsX()):
        h.SetBinContent(k+1,h.GetBinContent(k+1)/h.GetBinWidth(k+1))
        h.SetBinError(k+1,h.GetBinError(k+1)/h.GetBinWidth(k+1))
        
infile = TFile(opt.infile)
esimfile = TFile(opt.esimfile)
simfile = TFile(opt.simfile,"update")
bfile = open(opt.binfo)
efile = TFile(opt.efile)

gROOT.SetBatch(True)
gStyle.SetOptStat(0)

binfo = {}
for line in bfile:
    if line=="" or len(line) <2:
        continue
    if line[0]=='#' or "sim:" in line.split()[1]:
        continue
    bn = line.split()[0]
    binfo[bn] = map(lambda x: float(x),line.split()[2:])
    
#print binfo
if (not opt.acceptance):
    get_acceptance(simfile)
ofile = TFile(opt.outfile,"recreate")
ofile.mkdir("Multiplicities/M")
plot_dir([simfile.Get("M_ac")],binfo,ofile,["hM_ac_"],"ac_m")

bnarr = binfo.keys()
bnarr.sort()
bncode = reduce(lambda x,y: x + y,bnarr)
print "bncode: " + bncode

hname = "M/hM_" + bncode
he_rs = esimfile.Get("kinematics/hQ2x")
he_mc = esimfile.Get("kinematics/hQ2x_mc")
he = efile.Get("kinematics/hQ2x")
he_rs.Sumw2()
he_mc.Sumw2()
he.Sumw2()

he_ac = he_rs.Clone("he_ac")
he_ac.Divide(he_rs,he_mc,1,1,"b")

ofile.cd()
he_ac.Write("he_Q2x_ac",TObject.kOverwrite)
he_rs.Write("he_Q2x_rs",TObject.kOverwrite)
he_mc.Write("he_Q2x_mc",TObject.kOverwrite)
he.Write("he_Q2x_data",TObject.kOverwrite)
he.Divide(he_ac)
he.Write("he_Q2x_data_ac",TObject.kOverwrite)
plot_dir([infile.Get("M")],binfo,ofile,["hM_"],"data")

ofile.mkdir("Acceptance/M")
ofile.mkdir("data/M")
ofile.mkdir("dataAC/M")
ofile.mkdir("rs/M")
ofile.mkdir("mc/M")


for obj in infile.Get("M").GetListOfKeys():
    obj = obj.ReadObj()
    obj.Sumw2()
    ofile.cd("data/M")
    obj.Write("",TObject.kOverwrite)
    (q2bin,xbin) = get_Q2x_bin(obj,bnarr)
    bfact = get_bin_factor(obj,binfo) 
    NeDIS_fact = float(he.GetBinContent(xbin,q2bin))/he.GetXaxis().GetBinWidth(xbin)/he.GetYaxis().GetBinWidth(q2bin)
    hac = simfile.Get("M_ac/hM_ac_" + bncode + "_" + obj.GetName().split("_")[-1])
    ofile.cd('Acceptance/M')
    hac.Write('',TObject.kOverwrite)

    obj.Divide(hac) # acceptance correction
    ofile.cd("dataAC/M")
    obj.Write("",TObject.kOverwrite)

    if NeDIS_fact != 0:
        obj.Scale(1./NeDIS_fact/bfact)
    else:
        obj.Scale(0)

    scaleHist(obj)
    
    ofile.Get("Multiplicities/M").cd()
    obj.Write("",TObject.kOverwrite)
    
print 'creating canvas Multiplicities (m)...'
plot_dir([ofile.Get("Multiplicities/M")],binfo,ofile,["hM_"],"Mult_m")

print 'creating canvas with data (m)...'
plot_dir([ofile.Get("data/M")],binfo,ofile,["hM_"],"data_m")

print 'creating canvas with data + AC (m)...'
plot_dir([ofile.Get("dataAC/M")],binfo,ofile,["hM_"],"data_ac_m")


# Multiplicities
print 'generating Multiplicities (z) ...'
binfoM = get_depend(ofile,ofile.Get("Multiplicities/M"),binfo,"z")
plot_dir([ofile.Get("Multiplicities/z")],binfoM,ofile,["hz_"],"Mult_z")

print 'generating Multiplicities (x) ...'
binfoM = get_depend(ofile,ofile.Get("Multiplicities/M"),binfo,"x")
plot_dir([ofile.Get("Multiplicities/x")],binfoM,ofile,["hx_"],"Mult_x")

# acceptance
print 'generating acceptance (z) ...'
binfoM = get_depend(ofile,simfile.Get("M_ac"),binfo,"z","hM_ac_","Acceptance")
plot_dir([ofile.Get("Acceptance/z")],binfoM,ofile,["hz_"],"ac_z")

print 'generating acceptance (x) ...'
binfoM = get_depend(ofile,simfile.Get("M_ac"),binfo,"x","hM_ac_","Acceptance")
plot_dir([ofile.Get("Acceptance/x")],binfoM,ofile,["hx_"],"ac_x")

# data
print 'generating data (z) ...'
binfoM = get_depend(ofile,ofile.Get("data/M"),binfo,"z","hM_","data")
plot_dir([ofile.Get("data/z")],binfoM,ofile,["hz_"],"data_z")

print 'generating data (x) ...'
binfoM = get_depend(ofile,ofile.Get("data/M"),binfo,"x","hM_","data")
plot_dir([ofile.Get("data/x")],binfoM,ofile,["hx_"],"data_x")


# data + AC
print 'generating data + ac (z) ...'
binfoM = get_depend(ofile,ofile.Get("dataAC/M"),binfo,"z","hM_","dataAC")
plot_dir([ofile.Get("dataAC/z")],binfoM,ofile,["hz_"],"data_ac_z")

print 'generating data + ac (x) ...'
binfoM = get_depend(ofile,ofile.Get("dataAC/M"),binfo,"x","hM_","dataAC")
plot_dir([ofile.Get("dataAC/x")],binfoM,ofile,["hx_"],"data_ac_x")

#store simfile
store_sim(simfile)

# data/rs comparison
print 'generating data/rs comp (m) ...'
plot_dir([ofile.Get("data/M"),ofile.Get("rs/M")],binfo,ofile,["hM_","hM_"],"data_rs_cmp_m",[kBlack,kRed],["data","RS"],norm=True)

# simulation proyections.
print 'generating data/rs comp (z) ...'
binfoM = get_depend(ofile,simfile.Get("M"),binfo,"z","hM_","rs")
plot_dir([ofile.Get("data/z"),ofile.Get("rs/z")],binfoM,ofile,["hz_","hz_"],"data_rs_cmp_z",[kBlack,kRed],["data","RS"],norm=True)

print 'generating data/rs comp (x) ...'
binfoM = get_depend(ofile,simfile.Get("M"),binfo,"x","hM_","rs")
plot_dir([ofile.Get("data/x"),ofile.Get("rs/x")],binfoM,ofile,["hx_","hx_"],"data_rs_cmp_x",[kBlack,kRed],["data","RS"],norm=True)

ofile.Close()
    
print opt.outfile + " Done!"
