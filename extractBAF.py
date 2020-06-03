import sys
import os
import getopt
import argparse
import statistics
import pandas
import matplotlib.pyplot as plt
import random
import numpy as np
import copy
import datetime
import subprocess
from plotnine import ggplot, geom_point, aes
import multiprocessing
from joblib import Parallel, delayed
import plotBAF

class VCFRecord:
    def __init__(self,chrom,pos,posEnd):
        self.chrom = chrom
        self.pos = pos
        self.posEnd = posEnd
        self.ref = '.'
        self.alt = '.'
        self.varType = '.'
        self.totalDp = 0
        self.baf = 0
        self.mbaf = 0
        self.qual = 0
        self.gqual = 0
        self.hom = 'N'

    def newSNP(self ,vcfrow):
        self.ref = vcfrow[3]
        self.alt = '.'
        self.filter = vcfrow[7]
        format_ = vcfrow[8].split(':')
        info = vcfrow[9].split(':')
        gtype = info[format_.index('GT')]

        if 'DP' in format_:
            self.totalDp = int(info[format_.index('DP')])

        self.refDp = self.totalDp
        self.altDp = 0
        self.baf = 0

        if 'GQ' in format_:
            self.gqual = int(info[format_.index('GQ')])

        if not vcfrow[5] == '.':
            self.qual = float(vcfrow[5])

        if not gtype == '0/0':
            self.alt = vcfrow[4].split(',')[0]
            if len(self.alt) == len(self.ref):
                self.varType = 'SNV'
            else:
                self.varType = 'INDEL'
            if 'AD' in format_:
                depth = info[format_.index('AD')]
                self.refDp = int(depth.split(',')[0]) # A-allele depth
                self.altDp = int(depth.split(',')[1]) # B-allele depth

                if self.totalDp > 0:
                    self.baf = float(self.altDp)/self.totalDp
                    self.mbaf = abs(self.baf - 0.5)

    def qualCheck(self,minDp=10,mingQual=0,minQual=0):
        if self.totalDp >= minDp and self.qual > minQual and self.gqual > mingQual: 
            self.filterPass = True
        else:
            self.filterPass = False
			
def bafExtraction(name,chrom,vcfdir,bafdir,stats,minDp,mingQual):
    vcffile = '%s/%s_chr%s.vcf' %(vcfdir,name,chrom)
    out = open('%s/%s_chr%s.baf.raw' %(bafdir,name,chrom),'w')
    outFiltered = open('%s/%s_chr%s.baf' %(bafdir,name,chrom),'w')
    outFlag = open('%s/%s_chr%s.poorCovFlag' %(bafdir,name,chrom),'w')

    print ('%s - Calculating BAFs for chr%s' %(datetime.datetime.now().time(),str(chrom)))

    gqual = []
    baf = []

    flagstart = None
    flagend = None
    prevVar = None
    count = 0
    for vcfrow in open(vcffile,'r'):
        if not vcfrow.startswith('#'):
            vcfrow = vcfrow.split('\t')
            pos = int(vcfrow[1])
            if 'END' in vcfrow[7]:
                posEnd = int(vcfrow[7].split('=')[1])
            else:
                posEnd = pos

            var = VCFRecord(chrom,pos,posEnd)
            var.newSNP(vcfrow)

            if var.totalDp < 10:
                if not flagstart:
                    flagstart = pos
                flagend = posEnd
            else:
                if flagstart and flagend and flagend - flagstart >= 100:
                    outFlag.write('chr%s\t%i\t%i\tPoorCov\n' %(chrom,flagstart,flagend))
                flagstart = None
                flagend = None

            if prevVar and var.pos - prevVar.posEnd >= 100:
                outFlag.write('chr%s\t%i\t%i\tNoCov\n' %(chrom,prevVar.posEnd,var.pos-1))

            if not var.baf == 0 and var.varType == 'SNV':
                out.write('chr%s\t%i\t%s\t%s\t%i\t%i\t%.2f\n' %(var.chrom,var.pos,var.ref,var.alt,var.totalDp,var.gqual,var.baf))
                if var.totalDp >= minDp and var.gqual >= mingQual:
                    outFiltered.write('chr%s\t%i\t%s\t%s\t%i\t%i\t%.2f\n' %(var.chrom,var.pos,var.ref,var.alt,var.totalDp,var.gqual,var.baf))
                    count += 1

            prevVar = var

    stats.write('chr%s\t%i\n' %(chrom,count))
    print ('\t%s - Calculation completed for chr%s' %(datetime.datetime.now().time(),str(chrom)))

def snpPlot(name,chrom,bafdir,figdir):
    baffile = '%s/%s_chr%s.baf' %(bafdir,name,str(chrom))
    print('%s - Generating BAF scatter plots' %datetime.datetime.now().time())
    plotBAF.plotChromosome(name,chrom,baffile,figdir)
    plotBAF.plotQual(name,chrom,baffile,figdir)

def snpVisualize(name,chrom,bafdir,figdir):
    infile = '%s/%s_chr%s.baf' %(bafdir,name,str(chrom))
    print ('%s - Generating BAF scatter plots' %datetime.datetime.now().time())
    posList = []
    bafList = []
    regList = []
    qualList = []

    lines = open(infile,'r').readlines()[1:]
    for line in lines:
        if not line.startswith('#'):
            line = line.strip('\n').split('\t')
        else:
            continue

        baf = float(line[-1])
        if baf == 0 or baf == 1:
            next

        pos = int(line[1])
        bafList.append(baf)
        posList.append(pos)

    df = pandas.DataFrame({'BAF': bafList, 'Pos': posList})

    fig = plt.figure(figsize=(25,8))
    plt.plot(posList,bafList,'ko',markersize=0.1)
    plt.show()
    fig.savefig('%s/%s_%s_baf.png' %(figdir,name,chrom), dpi=fig.dpi)
    plt.close()
    print ('\t%s - Plot generated for %s' %(datetime.datetime.now().time(),chrom))
        
def main(argv):
    usage = """usage: python extractBAF.py  -n <sample_name> -i <VCF_File> 
                                            -o <output_directory> --minDp <min_depth>
                                            --mingQual <min_genotyping_quality> -v"""

    try:
        opts, args = getopt.getopt(argv, 'hn:i:o:v', ["minDp=", "mingQual=", "num-cores="])
    except getopt.GetoptError:
        print(usage)
        sys.exit()
         
    visualize = False
    ncores = 1

    minDp = 0
    mingQual = 0

    for opt, arg in opts:
        if opt == '-h':
            print (usage)
            sys.exit()
        elif opt == '-i':
            infile = arg
        elif opt == '-o':
            outdir = arg
        elif opt == '-n':
            name = arg
        elif opt == '-v':
            visualize = True			
        elif opt == '--minDp':
            minDp = int(arg)
        elif opt == '--mingQual':
            mingQual = int(arg)

    indir = os.path.dirname(os.path.abspath(infile))
    chromosomes = list(range(1,23)) + ['X','Y']
    vcfdir = '%s/vcf_split_by_chr' %indir
    bafdir = '%s/baf' %outdir
    figdir = '%s/baf_visualization' %outdir
    print (indir)
	
    if not os.path.exists(outdir):
        subprocess.call('mkdir %s' %outdir,shell=True)
    if not os.path.exists(vcfdir):
        subprocess.call('mkdir %s' %vcfdir,shell=True)
    if not os.path.exists(bafdir):
        subprocess.call('mkdir %s' %bafdir,shell=True)
    if not os.path.exists(figdir):
        subprocess.call('mkdir %s' %figdir,shell=True)
    subprocess.call('bcftools index %s' %infile,shell=True)

    for chrom in chromosomes:
        if not os.path.exists('%s/%s_chr%s.vcf' %(vcfdir,name,str(chrom))):
            subprocess.call('bcftools view -r chr%s %s -o %s/%s_chr%s.vcf' %(str(chrom),infile,vcfdir,name,str(chrom)),shell=True)

    stats = open('%s/%s.snpCounts' %(bafdir,name),'w')
    stats.write('Min depth = %i; Min genotyping quality = %i\n' %(minDp,mingQual))
    Parallel(n_jobs=ncores)(delayed(bafExtraction)(name,chrom,vcfdir,bafdir,stats,minDp,mingQual) for chrom in chromosomes)

    if visualize:
        Parallel(n_jobs=ncores)(delayed(snpPlot)(name,chrom,bafdir,figdir) for chrom in chromosomes)

if __name__ == "__main__":
    main(sys.argv[1:])
