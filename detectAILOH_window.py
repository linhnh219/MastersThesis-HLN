import sys
import os
import getopt
import statistics
import pandas
import matplotlib.pyplot as plt
import random
import numpy as np
from functools import reduce
import copy
import datetime
import subprocess
from plotnine import ggplot, geom_point, aes
import multiprocessing
from joblib import Parallel, delayed
import plotBAF

class SNPRecord:
    def __init__(self,chrom,pos):
        self.chrom = chrom
        self.pos = pos

    def newSNP(self,snprow):
        self.ref = snprow[2]
        self.alt = snprow[3]
        self.totalDp = int(snprow[4])
        self.gqual = int(snprow[5])
        self.baf = float(snprow[-1])
        self.mbaf = abs(self.baf - 0.5)

    def newSNPfromPop(self,snprow):
        self.bafs = []
        for i in range(1,len(snprow)):
            self.bafs.append(float(snprow[i]))

class SlidingWindow:
    def __init__(self,chrom,window_size,window_step):
        self.size = window_size
        self.step = window_step
        self.refSet = []
        self.chrom = chrom
        self.currentStart = None
        self.currentEnd = None
        self.currentCount = 0
        self.blocks = []
        self.regionStart = None
        self.regionEnd = None
        self.regionDetect = None

    def newWindow(self,start,end):
        self.currentStart = start
        self.currentEnd = end

    def newBlock(self,blockStart,blockEnd):
        self.blocks.append((blockStart,blockEnd,[]))

    def extend(self,refSNP):
        self.refSet.append(refSNP)

    def lohDetect(self,testSet):
        snp = [var.mbaf for var in testSet if var.baf > 0.1]
        het = [var.mbaf for var in testSet if var.mbaf <= 0.4]

        snp_mean = 0
        het_mean = 0
        bal_mean = 0

        if not snp == []:
            snp_mean = statistics.mean(snp)
        if not het == []:
            het_mean = statistics.mean(het)

        if snp == [] or snp_mean > 0.4:
            snpCount = len(snp)
            hetCount = len(het)
            validate = self.validateLOH(snpCount,hetCount)
            if validate:
                return ('LOH')
            else:
                return ('ROH')
        elif het_mean > 0.15:
            return ('AI')
        else:
            return ('N')

    def validateLOH(self,snpCount,hetCount):
        hetRatio = 0
        if not snpCount == 0:
            hetRatio = hetCount/snpCount

        def checkOutlier(list,x):
            list = list + [x]
            a = np.array(list)
            upper_quartile = np.percentile(a,75)
            lower_quartile = np.percentile(a,25)
            IQR = (upper_quartile - lower_quartile)*1.5

            if x <= (lower_quartile-IQR):
                return (True)
            else:
                return(False)
            
            return(resultList)

        def confin(mean,sd,sampleCount):
            return((mean-3*sd,mean+3*sd))

        if len(self.refSet) > 0:
            sampleCount = len(self.refSet[0].bafs)
            ref_hetCounts = [0]*sampleCount
            ref_snpCounts = [0]*sampleCount
            
            for snp in self.refSet:
                for i in range(sampleCount):
                    if snp.bafs[i] > 0.1:
                        ref_snpCounts[i] += 1
                        if snp.bafs[i] < 0.9:
                            ref_hetCounts[i] += 1

            ref_hetRatios = []
            for i in range(0,len(ref_snpCounts)):
                if not ref_snpCounts[i] == 0:
                    ref_hetRatios.append(ref_hetCounts[i]/ref_snpCounts[i])
                else:
                    ref_hetRatios.append(0)
            
        else:
            ref_snpCounts = []
            ref_hetRatios = []

        if checkOutlier(ref_snpCounts,snpCount) or checkOutlier(ref_hetRatios,hetRatio):
            return (True)
        else:
            return (False)

    def slide(self,nextStart,nextEnd,snpSet,out,outR):
        if self.blocks == []:
            self.newBlock(self.currentStart,nextStart)

        reportSet = [snp for snp in snpSet if snp.pos >= self.currentStart and snp.pos < self.currentEnd]

        for block in self.blocks:
            block[2].append(self.lohDetect(reportSet))

        popblock = self.blocks[0]

        outSet = [snp for snp in snpSet if snp.pos >= popblock[0] and snp.pos < popblock[1]]
        regionDetect = max(set(popblock[2]))
        self.writeRegion(popblock[0],popblock[1],regionDetect,snpSet,outR)
        self.blocks.remove(popblock)

        self.newBlock(self.currentEnd,nextEnd)

        while not self.refSet == [] and self.refSet[0].pos  < nextStart:
            self.refSet.remove(self.refSet[0])

        for snp in outSet:         
            out.write('chr%s\t%i\t%s\t%s\t%i\t%i\t%.2f\t%s\n' %(snp.chrom,snp.pos,snp.ref,snp.alt,snp.totalDp,snp.gqual,snp.baf,regionDetect))

        self.newWindow(nextStart,nextEnd)

    def writeRegion(self,start,end,regionDetect,snpSet,outR):
        if not self.regionStart and not self.regionEnd:
            self.regionStart = start
            self.regionEnd = end

        if not self.regionDetect:
            self.regionDetect = regionDetect

        else:
            if not regionDetect == self.regionDetect :
                outSet = [snp for snp in snpSet if snp.pos >= self.regionStart and snp.pos < self.regionEnd]

                if not self.regionDetect == 'N':
                    outR.write('chr%s\t%i\t%i\t%i\t%s\n' %(self.chrom,self.regionStart,self.regionEnd,len(outSet),self.regionDetect))
                self.regionStart = start
                self.regionEnd = end
                self.regionDetect = regionDetect
            else:
                self.regionEnd = end
        
    def empty(self,snpSet,out,outR):
        for block in self.blocks:
            block[2].append(self.lohDetect(snpSet))

        reportSet = [snp for snp in snpSet if snp.pos >= block[0] and snp.pos < block[1]]

        for block in self.blocks:
            regionDetect = max(set(block[2]))
            self.writeRegion(block[0],block[1],regionDetect,snpSet,outR)

        for snp in reportSet:
            out.write("chr%s\t%i\t%s\t%s\t%i\t%i\t%.2f\t%s\n" %(snp.chrom,snp.pos,snp.ref,snp.alt,snp.totalDp,snp.gqual,snp.baf,regionDetect))
            
        self.window = []
        self.currentStart = None
        self.currentEnd = None
        #print('\t%s - Calculation completed for %s' %(datetime.datetime.now().time(),currentChr))

def detectAI(name,chrom,chrStart,chrEnd,bafdir,refdir,outdir,window_size,window_step,minDp=0,mingQual=0):
    infile = '%s/%s_chr%s.baf' %(bafdir,name,chrom)
    if refdir:
        refSNPFile = '%s/chr%s.SNPset' %(refdir,chrom)

    out = open('%s/%s_chr%s.bafSeg' %(outdir,name,chrom),'w')
    outR = open('%s/%s_chr%s.segments' %(outdir,name,chrom),'w')

    print ('%s - Calculating BAFs and detecting ROH for chr%s' %(datetime.datetime.now().time(),str(chrom)))

    window = SlidingWindow(chrom,window_size,window_step)
    hetRatio = 0.6
    
    snpSet = []
    for line in open(infile,'r'):
        if not line.startswith('#'):
            line = line.strip('\n').split('\t')
            pos = int(line[1])
            snp = SNPRecord(chrom,pos)
            snp.newSNP(line)
            if snp.totalDp >= minDp and snp.gqual >= mingQual:
                snpSet.append(snp)

    if not window.currentStart:
        window.newWindow(chrStart,chrStart+window.size)

    if refdir:
        refSNPFile = '%s/chr%s.SNPset' %(refdir,chrom)
        sampleCount = None
        for line in open(refSNPFile,'r'):
            line = line.strip('\n').split('\t')
            pos = int(line[0])
            refSNP = SNPRecord(chrom,pos)
            refSNP.newSNPfromPop(line)

            while pos > window.currentEnd:
                nextStart = window.currentStart + window.step
                nextEnd = nextStart + window.size
                window.slide(nextStart,nextEnd,snpSet,out,outR)
            window.extend(refSNP)

    while window.currentEnd < chrEnd:
        nextStart = window.currentStart + window.step
        if (window.currentEnd + window.step) <= chrEnd:
            nextEnd = nextStart + window.size
        else:
            nextEnd = chrEnd
        window.slide(nextStart,nextEnd,snpSet,out,outR)

    window.empty(snpSet,out,outR)
    print ('\t%s - Calculation completed for chr%s' %(datetime.datetime.now().time(),str(chrom)))

def snpPlot(name,chrom,bafdir,outdir,figdir):
    chrom,chrStart,chrEnd = chrom
    bafFile ='%s/%s_chr%s.bafSeg' %(outdir,name,str(chrom))
    print('%s - Generating BAF scatter plots for chr%s' %(datetime.datetime.now().time(),str(chrom)))
    plotBAF.plotChromosome(name,chrom,bafFile,figdir)
    plotBAF.plotSegments(name,chrom,bafFile,figdir,chrStart,chrEnd)

def main(argv):
    usage = """usage: python detectAI.py -n <sample_name> -b <BAF_directory> 
                                        -o <output_directory]
                                        --windowSize <window_size>  
                                        --windowStep <window_step>"""
          
    try:
        opts, args = getopt.getopt(argv, 'hn:g:b:o:v', ["ref=","windowSize=", "windowStep=", "num-cores="])
    except getopt.GetoptError:
        print (usage)
        sys.exit()
         
    window_size = 50000
    visualize = False
    ncores = 1
    window_step = window_size
    refdir = None
    vcffile = None
    bafdir = None
    name = 'sample'

    for opt, arg in opts:
        if opt == '-h':
            print (usage)
            sys.exit()
        elif opt == '-g':
            chrCoord = arg
        elif opt == '-n':
            name = arg
        elif opt == '-b':
            bafdir = arg
        elif opt == '--ref':
            refdir = arg
        elif opt == '-o':
            outdir = arg
        elif opt == '--num-cores':
            ncores = int(arg)
        elif opt == '-v':
            visualize = True            
        elif opt == '--windowSize':
            window_size = int(arg)
        elif opt == '--windowStep':
            window_step = int(arg)

    if not window_step:
        window_step = window_size

    vcfdir = None
    if vcffile:
        indir = os.path.dirname(os.path.abspath(vcffile))
        vcfdir = '%s/vcf_split_by_chr' %indir

    figdir = '%s/baf_visualization' %outdir
    if not os.path.exists(outdir):
        subprocess.call('mkdir %s' %outdir,shell=True)
    if vcfdir and not os.path.exists(vcfdir):
        subprocess.call('mkdir %s' %vcfdir,shell=True)

    if not os.path.exists(figdir):
        subprocess.call('mkdir %s' %figdir,shell=True)

    if vcffile:
        subprocess.call('bcftools index %s' %vcffile,shell=True)

    chroms = []
    for line in open(chrCoord,'r'):
        line = line.strip('\n').split('\t')
        chrom = line[0].strip('chr')
        chrStart = int(line[1])
        chrEnd = int(line[2])
        chroms.append((chrom,chrStart,chrEnd))

        if vcffile and not os.path.exists('%s/%s_chr%s.vcf' %(vcfdir,name,str(chrom))):
            subprocess.call('bcftools view -r chr%s %s -o %s/%s_chr%s.vcf' %(str(chrom),vcffile,vcfdir,name,str(chrom)),shell=True)

    Parallel(n_jobs=ncores)(delayed(detectAI)(name,chrom[0],chrom[1],chrom[2],bafdir,refdir,outdir,window_size,window_step) for chrom in chroms)
    
    if visualize:
        Parallel(n_jobs=ncores)(delayed(snpPlot)(name,chrom,bafdir,outdir,figdir) for chrom in chroms)

if __name__ == "__main__":
    main(sys.argv[1:])
