import sys
import os
import getopt
import subprocess

def getInput(sampleName,baffile,bafdir,minDp=0,mingQual=0,hom_cutoff=0.1):
    inputFile = '%s.forCBS' %(baffile)
    input = open('%s.forCBS' %(baffile),'w')
    input.write('Name\tChr\tPosition\t%s.B Allele Freq\t%s.Log R Ratio\n' %(sampleName,sampleName))
    curPos = None
	
    for line in open(baffile,'r'):
        line = line.strip('\n').split('\t')
        chrom = line[0].strip('chr')
        pos = int(line[1])
        baf = float(line[6])
        totalDp = int(line[4])
        gQual = int(line[5])
        if curPos and curPos - pos > 10000:
            i = pos+1000
            while i < curPos:
                input.write('.\t%s\t%i\t0.00\t0\n' %(chrom,i))
                i += 1000
        if baf >= hom_cutoff and totalDp >= minDp and gQual >= mingQual:
            input.write('.\t%s\t%s\t%.2f\t0\n' %(chrom,pos,baf))
        curPos = pos
    return(inputFile)

def runcbs(inputFile,outdir,hom_cutoff=0.9,triplet_sum=1.5,ai_cutoff=0.65,snp_count=10):
    subprocess.call('perl split_samples.pl --data_file=%s' %inputFile,shell=True)
    subprocess.call('perl BAF_segment_samples.pl --non_informative=%.2f --ai_threshold=%.2f --triplet=%.2f --ai_size=%i' %(hom_cutoff,ai_cutoff,triplet_sum,snp_count), shell=True)

def getOutput(sampleName,outputFile,event,ai_cutoff=0.65,loh_cutoff=0.9):
    cbsfile = 'segmented/AI_regions.txt'

    curStart = None
    curEnd = None
    for line in open(cbsfile,'r'):
        if not line.startswith('Assay'):
            line = line.strip('\n').split('\t')
            chrom = line[1]
            start = int(line[2])
            end = int(line[3])
            snpCount = int(line[10])
            mbaf = float(line[6])

            if mbaf > loh_cutoff:
                segtype = 'LOH'
            elif mbaf > ai_cutoff:
                segtype = 'AI'

            if segtype == event:
                outputFile.write('chr%s\t%i\t%i\t%i\t%s\t%.2f\n' %(chrom,start,end,snpCount,segtype,mbaf))

            if not curStart:
                curStart = start
            if not curEnd:
                curEnd = end

    return (outputFile)
    
def main(argv):
    try:
        opts, args = getopt.getopt(argv, 'hn:g:b:o:v', ["cbs=", "ai=", "loh=", "hom="])
    except getopt.GetoptError:
        sys.exit()

    minDp = 0
    mingQual = 0
    ai_cutoff = 0.65
    loh_cutoff = 0.9

    for opt,arg in opts:
        if opt == '-h':
            sys.exit()
        elif opt == '-g':
            chrCoord = arg
        elif opt == '-n':
            name = arg
        elif opt == '-b':
            bafdir = arg
        elif opt == '-o':
            outdir = arg
        elif opt == '--cbs':
            cbsdir = arg
        elif opt == '-v':
            visualize = True
        elif opt == '--hom':
            hom_cutoff = float(arg)
        elif opt == '--ai':
            ai_cutoff = float(arg)
        elif opt == '--loh':
            loh_cutoff = float(arg)

    if not os.path.exists(outdir):
        subprocess.call('mkdir %s' %outdir,shell=True)

    chroms = list(range(1,23)) + ['X','Y']
    for chrom in chroms:
        sampleName = '%s_chr%s' %(name,chrom)
        baffile = '%s/%s_chr%s.baf' %(bafdir,name,chrom)
        inputFile = getInput(sampleName,baffile,bafdir,minDp,mingQual)
        outputfile = open('%s/%s.cbs.segments' %(outdir,sampleName),'w')
        runcbs(inputFile,outdir,hom_cutoff)
        getOutput(sampleName,outputfile,'AI',ai_cutoff,loh_cutoff)
        runcbs(inputFile,outdir)
        getOutput(sampleName,outputfile,'LOH',ai_cutoff,loh_cutoff)

if __name__ == "__main__":
    main(sys.argv[1:])

