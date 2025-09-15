import importlib
import math
from matplotlib import pyplot
import matplotlib
import pandas as pd
import random
import numpy
import numpy as np
import numpy.matlib
from datetime import datetime
import seaborn
# import bokeh
import pylab
import matplotlib.ticker as mtick
from collections import Counter
from fractions import Fraction  
import decimal
from matplotlib.ticker import PercentFormatter
import matplotlib.pylab as plt
from matplotlib import pylab as plt
import matplotlib.ticker as tkr     # has classes for tick-locating and -formatting

 

def SimulationFunc(gennum,cell_num,ProbBind,ProbMut, gRNAFile, GFPseqF_original, Uniformity,Shift):
        #import gRNAs from CSV file. DEFINE gRNAs
    gfplist=[]
    gRNAsDF = pd.read_csv(gRNAFile, sep = ',')
    # count number of possible mutation on gRNAs
    mutation_limit=0
    grna=mutation_limit
    GFPforCounting=GFPseqF_original
    for index, row in gRNAsDF.iterrows():   
        Start = row['Start']-Shift
        End = row['End']-Shift
        Orientation = row['Orientation']
        gRNAseq = GFPforCounting[Start-1:End]
        gRNAseqList = list(gRNAseq)

        x=0 
        while x < len(gRNAseqList): 
            if Orientation == '+'and gRNAseqList[x]=="c":
                    mutation_limit+=1
                    gRNAseqList[x]='*'
            if Orientation == '-' and gRNAseqList[x]=="g":
                    mutation_limit+=1
                    gRNAseqList[x]='*'
            x=x+1
        GFPforCounting = GFPforCounting.replace(gRNAseq,''.join(gRNAseqList))
    #print ("mutation limit number of C = " + str(mutation_limit))
    # = mutation_limit
    Mutation_matrix = numpy.zeros([gennum,cell_num])

    # accessing gRNAs in GFP
    cell=0
    while cell<cell_num:
        GFPseqF = GFPseqF_original
        
        loop=0
        listforgraph=[]
        while loop < gennum:
            GFPseqR = GFPseqF[:]
        # gRNA locating function
            for index, row in gRNAsDF.iterrows():   
                Start = row['Start']-Shift
                End = row['End']-Shift
                Orientation = row['Orientation']
                if Orientation=='+': 
                    gRNAseq = GFPseqF[Start-1:End] #Forward
                    UniformityNew = Uniformity
                if Orientation=='-': 
                    gRNAseq = GFPseqR[Start-1:End] #Reversed
                    UniformityNew = Uniformity[::-1]
                gRNAseqList = list(gRNAseq)
                mutation_counter=[]
                mutation_counter.extend([i for i,x in enumerate(gRNAseqList) if x == '*']) #adds (in the nonmath sense) function output to mutation counter

            # PROBABILITY OF Cas9 BINDING
                x=0
                if ProbBind > random.random() and len(mutation_counter) == 0:
#                     if len(gRNAseqList != len(Uniformity):
                           
                    while x < len(gRNAseqList): 
                        # PROBABILITY OF AID MUTATING C
#                         if Orientation == '-':
#                             AdditionalProb = Uniformity[max(-x,-len(Uniformity)+1)]
#                         else:
#                             AdditionalProb = Uniformity[min(x,len(Uniformity)-1)]
                        AdditionalProb = UniformityNew[x]
                        if Orientation == '+' and gRNAseqList[x]=="c" and random.random()<ProbMut*AdditionalProb:
                            gRNAseqList[x]='*'
                        if Orientation == '-' and gRNAseqList[x]=="g" and random.random()<ProbMut*AdditionalProb:
                            gRNAseqList[x]='*'
                        x=x+1
                    if Orientation=='+': GFPseqF = GFPseqF.replace(gRNAseq,''.join(gRNAseqList))
                    if Orientation=='-': GFPseqR = GFPseqR.replace(gRNAseq,''.join(gRNAseqList)) 
                
                # STORE MUTATIONS. replace C>T original GFP (gRNA check function). easier alternative replace C > *. next loop counts how many * 
            #Counting total number of mutations after gRNA mutations:
            if random.random()<0.5:
                GFPseqF=GFPseqF[:]
                                            #Equivelent to cell division (trying to keep forward and reverse)
            else:
                GFPseqF=GFPseqR[:]
            mutation_index=[]
            #all * positions:
            mutation_index.extend([i for i,x in enumerate(GFPseqF) if x == '*'])
            #print(mutation_index)
            mutation_number=len(mutation_index)
            #print (mutation_number)
            listforgraph.append(mutation_number/mutation_limit)
            Mutation_matrix[loop,cell] = mutation_number       

            loop=loop + 1

        gfplist.extend([GFPseqF])
        cell=cell+1
    #printing mutated GFP sequence after one cell at the end of all the generation number.
    #print("simulation done")
    return(Mutation_matrix,mutation_limit,gRNAsDF,gfplist)


def PlotFunct(cell_num,fullpath,Mutation_matrix,mutation_limit,ProbMut,ProbBind,gRNAsDF,gennum,gfplist,GFPseqF_original,cfig,ax,cfig2,ax2,figurehistogram,axhist,figureplot,axplot,MutIndex,BindIndex,gRNAFile,dfAvMut):
    PlotIndex = MutIndex + BindIndex
    now = datetime.now() #initial variable for time sys
    current_time = now.strftime("%H--%M--%S.%f"+"_" + "%m-%d-%Y") #time format: Hour--Minute--Second and femtosecond_month-day-year
    cfig.tight_layout()
    avgmatrix=Mutation_matrix/mutation_limit
    df = pd.DataFrame(avgmatrix)

    Legend = "Pmut=" + str(ProbMut)

    dataframevar = numpy.mean(avgmatrix, 1)
    df2 = pd.DataFrame (numpy.transpose(dataframevar), columns=[Legend]) 
    matplotlib.pyplot.figure(cfig, set_visible = False)
    df2.plot(color='r', linewidth=2.0, marker="*", ax=ax[PlotIndex])
    from matplotlib import rcParams
    rcParams['axes.titlepad'] = 20 
    #matplotlib.pyplot.subplots_adjust(hspace=0)
    df.plot(ylabel='Mutation', xlabel='Generation Number',title= 'MutP=' + str(ProbMut) + ' | ' + 'BindP='+ str(ProbBind) + ' | ' + 'gRNA=' + str(len(gRNAsDF)), alpha=0.2,color = 'k', ax=ax[PlotIndex])
    ax[PlotIndex].get_legend().remove()
    matplotlib.pyplot.ylim(0, 1)

    matplotlib.pyplot.figure(cfig2, set_visible = False)
    cfig2.tight_layout()
    from matplotlib import rcParams
    rcParams['axes.titlepad'] = 20 
    #matplotlib.pyplot.subplots_adjust(hspace=0)
    #print(df2)
    df2.plot(linewidth=2.0, ylabel='Mutation', xlabel='Generation Number',title=  'Pbind='+ str(ProbBind) + ' | ' + 'gRNA=' + str(len(gRNAsDF)), ax=ax2[BindIndex])
    matplotlib.pyplot.ylim(0, 1)

    dfAvMut.loc[PlotIndex] = [ProbMut,dataframevar[-1]]


    #########################################
 
    import matplotlib.pyplot as plt
    #----
    #Histogram: 
    #For each simulated cell we want to count the number of mutations it has
    #and plot the distribution of mutations in a histogram 
    Mutation_matrix #each row=generation, each column is a cell
    matplotlib.pyplot.figure(figurehistogram, set_visible = False)
    matplotlib.pyplot.axes(axhist[PlotIndex])
    #Select last row for mutation number in cells
    data=Mutation_matrix[-1,:]
    #axhist[PlotIndex]=seaborn.displot(data=data, stat='percent', height=4, kde=True, kind='hist',binwidth=1)
    figurehistogram.tight_layout()
    try:
        seabornplot=seaborn.histplot(data=data, stat='percent',binwidth=1,kde=False,legend=True,edgecolor='white',linewidth = 1,color="blue",alpha=0.5)
    except:
        seabornplot=seaborn.histplot(data=data, stat='percent',bins=1,kde=False,legend=True,edgecolor='white',linewidth = 1,color="blue",alpha=0.5)
    axhist[PlotIndex].set(xlabel='Mutation Number',
        ylabel='Percentage',
        title='MutP=' + str(ProbMut) + ' | ' + 'BindP='+ str(ProbBind) + ' | ' + 'gRNA=' + str(len(gRNAsDF))+'\n' ,xlim=(0,mutation_limit),ylim=(0,30)) 
    

    
    ######################################

    starlisterC=[]
    starlisterG=[]
    #For eACH GFP position what proportion of the cells have a mutation in that position
    #gfplist
    #* are mutations
    #print(gfplist)
    counter=0
    while counter < len(gfplist):
        GFPseqF=gfplist[counter]
#         GFPseqF.count('*') #the number of mutations
#         mutation_GFPseqF_list_ofwhatsleft = set(GFPseqF).symmetric_difference(GFPseqF_original)
        #print(list(mutation_GFPseqF_list_ofwhatsleft))
        for position,char in enumerate(GFPseqF):
            OriginalBase = GFPseqF_original[position]
            if char=='*' and OriginalBase=='c':
                #print(position)
                starlisterC.append(position) #getting position of *s
            if char=='*' and OriginalBase=='g':
                #print(position)
                starlisterG.append(position) #getting position of *s
        counter=counter+1
    starlisterC.sort()
    starlisterG.sort()
    Cs = Counter(starlisterC)
    Gs = Counter(starlisterG)
    
    starlisterdictionaryC = {}
    starlisterdictionaryG = {}

    newcount=0
    while len(starlisterC)>newcount:
        if newcount==0:
             starlisterdictionaryC = { starlisterC[newcount]: Cs[ starlisterC[newcount]]/cell_num*100}
        else:
             starlisterdictionaryC.update({ starlisterC[newcount]: Cs[ starlisterC[newcount]]/cell_num*100})
        newcount=newcount+1
        

    newcount=0
    while len(starlisterG)>newcount:
        if newcount==0:
             starlisterdictionaryG = { starlisterG[newcount]: Gs[ starlisterG[newcount]]/cell_num*100}
        else:
             starlisterdictionaryG.update({ starlisterG[newcount]: Gs[ starlisterG[newcount]]/cell_num*100})
        newcount=newcount+1
    
    namesC = list(starlisterdictionaryC.keys())
    valuesC = list(starlisterdictionaryC.values())
    namesG = list(starlisterdictionaryG.keys())
    valuesG = list(starlisterdictionaryG.values())
    #
    matplotlib.pyplot.figure(figureplot, set_visible = False)
    matplotlib.pyplot.axes(axplot[PlotIndex])
    matplotlib.pyplot.xlabel("Position") 
    matplotlib.pyplot.ylabel("Percent Mutation")
    #figureplot, axplot[PlotIndex] = matplotlib.pylab.subplots(1, tight_layout=True)
    axplot[PlotIndex].bar(namesC, valuesC, edgecolor='none', color = '#FED884')
    axplot[PlotIndex].bar(namesG, valuesG, edgecolor='none', color = '#7CCAA5')

    #matplotlib.pyplot.bar(names, values)
    #matplotlib.pylab.xticks(rotation=90)
    matplotlib.pyplot.xlim(0,len(GFPseqF))
    matplotlib.pyplot.ylim(0,100)
    matplotlib.pyplot.title('MutP=' + str(ProbMut) + ' | ' + 'BindP='+ str(ProbBind) + ' | ' + 'gRNA=' + str(len(gRNAsDF)) + '\n AvMut = '+str(round(dataframevar[-1]*100,1))+'%')
    #plt.label(xlabel='Mutation Number', ylabel='Percentage')
    figureplot.tight_layout()

    #print("Plotting Done")

    return(cfig,ax,cfig2,ax2,figurehistogram,axhist,figureplot,axplot,dfAvMut)
