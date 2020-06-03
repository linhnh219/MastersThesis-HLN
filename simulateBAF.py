import random
import numpy as np
import statistics

def generateSegments(chrStart,chrEnd,refSNPFile,n):
    S = []
    s = chrStart
    T = np.array([0,1,2,3])
    C = []
    E = list(np.random.randint(chrStart,chrEnd,n-1))
    E = sorted(E)
    E.append(chrEnd)

    i = 0
    snpCounts = None
    ai_rate = 0.01
    loh_rate = 0.01
    for line in open(refSNPFile,'r'):
        line = line.strip('\n').split('\t')
        pos = int(line[0])
        if not snpCounts:
            snpCounts = [0]*(len(line)-1)

        if i < len(E):
            e = E[i]
            t = int(np.random.choice([1,2,3],1,p=[loh_rate,1-loh_rate-ai_rate,ai_rate]))
            l = e-s+1
  
            if pos > e:
                a = np.array(snpCounts)
                x = []
                upper_quartile = np.percentile(a,75)
                lower_quartile = np.percentile(a,25)
                IQR = (upper_quartile - lower_quartile)*1.5
                for snpCount in snpCounts:
                    if snpCount > (lower_quartile-IQR) or snpCount < (upper_quartile+IQR):
                        x.append(snpCount)
                if not x == []:
                    c = random.randrange(min(x),max(x)+1)
                else:
                    c = abs(int(l*np.random.normal(0,0.00001)))

                S.append((s,e,t,c))
                s = e+1
                i += 1
                snpCounts = None
            
            elif pos in range(s,e+1):
                for j in range(1,len(line)):
                    if float(line[j]) > 0.1:
                        snpCounts[j-1] += 1

    if i < len(E):
        for j in range(i,len(E)):
            e = E[j]
            l = e-s+1
            c = abs(int(l*np.random.normal(0.0,0.00001)))
            t = 0
            S.append((s,e,t,c))
            s = e+1

    return (S)

def generateBaf(start,end,rtype,snpCount):
    B = []
    hetCount = -1
    aiCount = -1

    if rtype == 2 or rtype == 3:
        while hetCount < 0:
            hetCount = int((0.6 + np.random.normal(0,0.1))*snpCount)

        if rtype == 3:
            while aiCount < 0:
                aiCount = int((1 - abs(np.random.normal(0,0.05)))*hetCount)
                hetCount = hetCount - aiCount

    elif rtype == 1 or rtype == 0: # Partial loss
        while hetCount < 0:
            hetCount = abs(int(np.random.normal(0,0.05)*snpCount))

    refhomCount = np.random.normal(0,0.01)*(snpCount - hetCount)
    homCount = snpCount - hetCount - refhomCount

    for i in range(int(hetCount)):
        B.append(0.5 + np.random.normal(0,0.075))

    for i in range(int(homCount)):
        B.append(1 - abs(np.random.normal(0,0.01)))

    for i in range(int(refhomCount)):
        B.append(0 + abs(np.random.normal(0,0.01)))

    if aiCount > 0:
        cn = random.randrange(3,7)
        ai_t = float(1)/cn
        for i in range(int(aiCount)):
            B.append(random.choice([1-ai_t,ai_t]) + np.random.normal(0,0.075))

    homRefCount = end-start+1-snpCount

    for i in range(homRefCount):
        B.extend([0])

    B = np.array(B)
    B = np.random.permutation(B)
    return (B)
            
def main():
    #chrCoord = "/home/halinh_nguyen/data/chrCoord_hg19.txt"

    for line in open(chrCoord,'r'):
        line = line.strip('\n').split('\t')
        chrom = line[0]
        chrStart = int(line[1])
        chrEnd = int(line[2])

        n = 10000
        refSNPFile = "/home/halinh_nguyen/data/1000KG/%s.SNPset" %chrom
        S = generateSegments(chrStart,chrEnd,refSNPFile,n)
        segout = open("simulated_%s.segments" %chrom,'w')
        bafout = open("simulated_%s.baf" %chrom,'w')
        print(chrom)
        for segment in S:
            s,e,t,c = segment
            segout.write('%s\t%i\t%i\t%s\n' %(chrom,s,e,t))
            B = generateBaf(s,e,t,c)
            for i in range(len(B)):
                if not B[i] == 0:
                    bafout.write('%s\t%i\t.\t.\t40\t99\t%.2f\n' %(chrom,s+i,B[i]))

if __name__ == '__main__':
    main()
