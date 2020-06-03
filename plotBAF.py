import pandas
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from plotnine import *
from plotnine import ggplot, geom_point, aes
import warnings

def plotCNVCall(name,chrom,start,end,baffile,outdir):
    len = end - start
    posList = []
    bafList = []

    for line in open(baffile,'r'):
        line = line.strip('\n').split('\t')
        baf = float(line[6])
        pos = int(line[1])
        mbaf = abs(0.5-baf)

        if pos < start-0.5*len:
            next            
        if pos > end+0.5*len:
            outfile = '%s/%s_%s_%i-%i.png' %(outdir,name,chrom,start,end)

            fig = plt.figure(figsize=(30,10))
            plt.ylim(-0.02,1.02)
            plt.xlim(start-0.5*len,end+0.5*len)
            plt.plot(posList,bafList,'ko',markersize=8.0)
            plt.axvline(x=start)
            plt.axvline(x=end)
            plt.xlabel('POS',fontsize=25)
            plt.ylabel('BAF',fontsize=25)
            plt.show()
            fig.savefig(outfile, dpi=fig.dpi)
            plt.close()

            posList = []
            bafList = []
            break
        
        if pos in range(int(start-0.5*len),int(end+0.5*len)):
            posList.append(pos)
            bafList.append(baf)

def plotSegments(name,chrom,baffile,outdir,chrStart,chrEnd):
    segments = []
    chrom = 'chr%s' %chrom

    i = chrStart
    posList = []
    bafList = []
    typeList = []
    df = pandas.DataFrame(columns=['Pos','BAF','Type'])

    for line in open(baffile,'r'):
        line = line.strip('\n').split('\t')
        baf = float(line[-2])
        pos = int(line[1])
        segtype = line[-1]

        if pos > i+2000000:
            df = pandas.DataFrame({'BAF': bafList, 'Pos': posList, 'Type': typeList})
            if (i+2000000) < chrEnd:
                e = i + 2000000

            else:
                e = chrEnd
            outfile = '%s/%s_%s_%i-%i.png' %(outdir,name,chrom,i,e)
            plotPlt(df,outfile)
            i += 2000000
            bafList = []
            posList = []
            typeList = []

        if pos in range(i,i+2000000):
            bafList.append(baf)
            posList.append(pos)
            typeList.append(segtype)

def plotPlt(df,outfile):
    subfig, subax = plt.subplots(figsize=(25,10))
    plt.ylim(-0.02,1.02)
    colors = {'LOH': '#CC0000', 'AI': '#009900', 'N': 'black', 'ROH': 'grey'}
    grouped = df.groupby('Type')
    for key, group in grouped:
        group.plot(ax=subax, kind='scatter', x='Pos', y='BAF', label=key, color=colors[key], s=20)
    subfig.savefig(outfile, dpi=subfig.dpi)
    plt.close()

def gg_plot(outfile,df,s,e):
    base_plot = [aes(x='Pos',y='BAF'), aes(ymin=0), geom_point(aes(color='Type'),size=1),]
    plot = ggplot(df) + base_plot
    plot.save(outfile,height=8,width=25)
    warnings.simplefilter("ignore", category=matplotlib.cbook.mplDeprecation)
    warnings.simplefilter("ignore", category=UserWarning)

def plotChromosome(name,chrom,baffile,outdir):
    bafList = []
    posList = []
    for line in open(baffile,'r'):
        line = line.strip('\n').split('\t')
        baf = float(line[6])
        pos = int(line[1])
        totalDp = int(line[4])
        gQual = int(line[5])

        bafList.append(baf)
        posList.append(pos)

    outfile = '%s/%s_chr%s.png' %(outdir,name,chrom)
    fig = plt.figure(figsize=(30,10))
    plt.ylim(-0.02,1.02)
    plt.plot(posList,bafList,'ko',markersize=1.0)
    plt.show()
    fig.savefig(outfile, dpi=fig.dpi)
    plt.close()

def plotQual(name,chrom,baffile,outdir):
    bafList = []
    qualList = []
    dpList = []
    for line in open(baffile,'r'):
        line = line.strip('\n').split('\t')
        baf = float(line[6])
        pos = int(line[1])
        gQual = int(line[5])
        totalDp = int(line[4])

        bafList.append(baf)
        qualList.append(gQual)
        dpList.append(totalDp)

    fig1 = plt.figure()
    plt.plot(bafList,qualList,'ko',markersize=0.3)
    plt.show()
    fig1.savefig('%s/%s_chr%s.bafvsqual.png' %(outdir,name,chrom),dpi=fig1.dpi)
    plt.close()

    fig2 = plt.figure()
    plt.plot(bafList,dpList,'ko',markersize=0.3)
    plt.ylim(0,200)
    plt.show()
    fig2.savefig('%s/%s_chr%s.bafvsdepth.png' %(outdir,name,chrom),dpi=fig2.dpi)
    plt.close()

    fig3 = plt.figure()
    plt.hist(qualList,bins=[0,10,20,30,40,50,60,70,80,90,100])
    fig3.savefig('%s/%s_chr%s.qualHist.png' %(outdir,name,chrom),dpi=fig3.dpi)
    plt.close()

    fig4 = plt.figure()
    plt.hist(dpList,bins=[0,10,20,30,40,50,60,70,80,90,100,200])
    fig4.savefig('%s/%s_chr%s.depthHist.png' %(outdir,name,chrom),dpi=fig4.dpi)
    plt.close()



