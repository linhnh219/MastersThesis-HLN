import sys
import os
import getopt
import subprocess

def getInput(sampleName,baffile,bafdir,minDp=0,mingQual=0,hom_cutoff=0.1):
    inputFile = '%s.forCBS' %(baffile)
    input = open('%s.forCBS' %(baffile),'w')
    curPos = None
	
    for line in open(baffile,'r'):
        line = line.strip('\n').split('\t')
        chrom = line[0].strip('chr')
        pos = int(line[1])
        baf = float(line[6])
        ref = line[2]
        alt = line[3]
        totalDp = int(line[4])
        gQual = int(line[5])

        if curPos and pos-curPos > 10000:
            i = curPos + 1000
            while i < pos:
                input.write('%s\t%i\t.\t.\t.\t.\t0.00\n' %(chrom,i))
                i += 1000
        if baf >= hom_cutoff and totalDp >= minDp and gQual >= mingQual:
            input.write('%s\t%i\t%s\t%s\t%i\t%i\t%.2f\n' %(chrom,pos,ref,alt,totalDp,gQual,baf))
        curPos = pos
		
    return(inputFile)

def runcbs(inputFile,outdir,hom_cutoff):
    rscript = 'cbsBAF.R'
    FNULL = open(os.devnull,'w')
    subprocess.call('Rscript --vanilla %s %s %.2f' %(rscript,inputFile,hom_cutoff),stdout=FNULL,shell=True)

def getOutput(sampleName,outdir,outputFile,event,ai_cutoff,loh_cutoff):
    cbsfile = 'dnacopy_results.txt'
    curStart = None
    curEnd = None
	
    for line in open(cbsfile,'r'):
        if not line.startswith('Assay'):
            line = line.strip('\n').split('\t')
            chrom = line[0]
            start = int(line[1])
            end = int(line[2])
            snpCount = int(line[3])
            mbaf = float(line[4])

            if mbaf > loh_cutoff:
                segtype = 'LOH'
            elif mbaf > ai_cutoff:
                segtype = 'AI'
            else:
                segtype = 'N'

            if segtype == event:
                outputFile.write('chr%s\t%i\t%i\t%i\t%s\t%.2f\n' %(chrom,start,end,snpCount,segtype,mbaf))

            if not curStart:
                curStart = start
            if not curEnd:
                curEnd = end

    return (outputFile)
    
def main(argv):
    try:
        opts, args = getopt.getopt(argv, 'hn:g:b:o:v', ["cbs=", "minDp=", "mingQual=", "ai=", "loh=", "hom="])
    except getopt.GetoptError:
        sys.exit()

    minDp = 0
    mingQual = 0
    ai_cutoff = 0.15
    loh_cutoff = 0.4

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
        elif opt == '--minDp':
            minDp = int(arg)
        elif opt == '--mingQual':
            mingQual = int(arg)
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
        print ('Segmenting chr%s' %chrom)
        sampleName = '%s_chr%s' %(name,chrom)
        baffile = '%s/%s_chr%s.baf' %(bafdir,name,chrom)
        subprocess.call('sort -k2 -n %s > %s_' %(baffile,baffile),shell=True)
        subprocess.call('rm %s' %baffile,shell=True)
        subprocess.call('mv %s_ %s' %(baffile,baffile),shell=True)
		
        inputFile = getInput(sampleName,baffile,bafdir,minDp,mingQual)
        outputfile = open('%s/%s.cbs.segments' %(outdir,sampleName),'w')
        runcbs(inputFile,outdir,hom_cutoff)
        getOutput(sampleName,outdir,outputfile,'AI',ai_cutoff,loh_cutoff)
        runcbs(inputFile,outdir,hom_cutoff=1)
        getOutput(sampleName,outdir,outputfile,'LOH',ai_cutoff,loh_cutoff)

if __name__ == "__main__":
    main(sys.argv[1:])

