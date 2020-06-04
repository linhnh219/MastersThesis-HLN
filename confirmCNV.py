import sys
import os
import getopt
import subprocess
from statistics import median
import pandas
from plotnine import ggplot, geom_point, aes
import matplotlib.pyplot as plt

def checkResults(intersectFile):
    outfile = open('%s.correctCalls' %intersectFile,'w')
    curCall = None
    segList = []
    segLen = 0
	
    for line in open(intersectFile,'r'):
        line = line.strip('\n').split('\t')
        chrom = line[0]
        start = int(line[1])
        end = int(line[2])
        snpCount = int(line[3])
        callType = line[4]
        segStart = int(line[7])
        segEnd = int(line[8])
        stype = line[9]
        seg = ((segStart,segEnd,stype))
        call = ((start,end,snpCount,callType))

        if curCall and not call == curCall:
            callLen = curCall[1] - curCall[0]
            if not segList == [] and segLen >= (0.6*callLen):
                for seg_ in segList:
                    outfile.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' %(chrom,curCall[0],curCall[1],curCall[2],curCall[3],seg_[0],seg_[1],seg_[2]))
            curCall = None
            segList = []
            segLen = 0

        if not curCall:
            curCall = call

        loh = False
        ai = False
        if stype == 'Loss' and curCall[3] in ['LOH','ROH']:
            loh = True
        
        if stype == 'Gain' and curCall[3] == 'AI':
            ai = True

        if loh or ai:
            segList.append(seg)
            x = sorted([curCall[0],curCall[1],seg[0],seg[1]])
            segLen += x[2]-x[1]
    
def checksegFile(segFile):
    newsegFile = open('%s_' %segFile,'w')
    for line in open(segFile,'r'):
        line_ = line.split('\t')
        if not int(line_[2]) < int(line_[1]):
            newsegFile.write(line)
    subprocess.call('rm %s' %segFile,shell=True)
    subprocess.call('mv %s_ %s' %(segFile,segFile),shell=True)
    
def main(argv):
    try:
        opts, args = getopt.getopt(argv, 'hn:c:o:s:')
    except getopt.GetoptError:
        sys.exit()

    for opt, arg in opts:
        if opt == '-n':
            name = arg
        elif opt == '-c':
            callFile = arg
        elif opt == '-s':
            segFile = arg
        elif opt == '-o':
            outdir = arg

    chroms = list(range(1,23)) + ['X','Y']
    for chrom in chroms:
        chrom = 'chr%s' %str(chrom)
        callFile = '%s_%s.cbs.segments' %(name,chrom)
        checksegFile(segFile)
        intersectFile = '%s_%s.intersect' %(callFile,chrom)
        subprocess.call('bedtools intersect -a %s -b %s -wa -wb > %s/%s' %(callFile,segFile,outdir,intersectFile),shell=True)
        subprocess.call('sort %s/%s > %s/%s.sorted' %(outdir,intersectFile,outdir,intersectFile),shell=True)
        intersectFile = '%s_%s.intersect.sorted' %(callFile,chrom)
        checkResults(intersectFile)

if __name__ == "__main__":
    main(sys.argv[1:])



                    


