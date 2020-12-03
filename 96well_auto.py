"""
Calculation and visualization for 96-well MTP screening
Author: Paul Zurek, pjz26@cam.ac.uk
From v2.1 - less convoluted/functional as v1, but automatic linear range!
22/10/2018 added argparse
23/10/2018 added xrange select
24/10/2018 added jump detection
09/04/2010 changed to pandas to read files, so that incomplete plates can be analysed
v2.2a
Python 3
"""
import numpy as np
import pandas as pd
import scipy.stats as sp_stats
import matplotlib.pyplot as plt
import seaborn as sns
import argparse


parser = argparse.ArgumentParser(description="""96-well plate screening kinetic analysis\n
                                 Author: Paul Zurek (pjz26@cam.ac.uk)
                                 v2.2 24/10/2019""")

parser.add_argument('input_file', help='96 well data\nColumns: Time (hh:mm:ss) followed by well data')
parser.add_argument('-s', '--single', help='Plot single wells (e.g. --single A4 G7)',
                    nargs='*')
parser.add_argument('--min_r2', help='Set min r2 for auto linear range', 
                    type=float, default=0.99)
parser.add_argument('--min_n', help='Set min number of points for auto linear range', 
                    type=int, default=5)
parser.add_argument('--break_thresh', type=float, default=7,
                    help='Set threshold of break detection')
parser.add_argument('--crop', help='Crop by time (in minutes, e.g. --crop 0 5 to select first 5 minutes)', 
                    type=float, nargs=2)
parser.add_argument('--full', action='store_true',
                    help='No slope adjustment, full data fitting')
parser.add_argument('--negative', action='store_true',
                    help='Negative sloped assay')
parser.add_argument('-o', '--output', help='flag to save slopes to file', 
                    action='store_true')
parser.add_argument('-v', '--version', action='version', version='2.2a')

args = parser.parse_args()

fileinput = args.input_file
name = fileinput[:-4]

singles = args.single
singles_flag = False
if singles is not None: singles_flag = True

min_measure = args.min_n
min_r2 = args.min_r2
scale_factor = 60   #Scales to slope h-1

full_flag = args.full
out_flag = args.output

neg_flag = args.negative
neg_point = 2

crop = args.crop

break_thresh = args.break_thresh

#Convert well name into well number
def convert_well(cs):
    l = "ABCDEFGH"
    conv = []
    for e in cs:
        conv.append(l.index(e[0]) * 12 + int(e[1:])-1)
    return conv


#Get Data
data = pd.read_csv(fileinput, sep="\t", index_col=False).fillna(0)

#Just the columns that have the wells
wells = data.columns.values[1:]  
wells = list(wells)
#Get data values and transpose
values = np.transpose(data[wells].values)
values = [list(l) for l in values]

#Convert if negative assay
if neg_flag:
    for i in range(len(values)):
        for j in range(len(values[i])):
            values[i][j] = neg_point - values[i][j]
            

#Extract times: Here format is hh:mm:ss (or only mm:ss)
#Adjust if necessary
timestr = data.iloc[:,0].values  
times = []
for time in timestr:
    t = time.split(":")
    try:
        if len(t) == 2:
            times.append(int(t[0]) + (int(t[1])/60.0))   #-> convert time to minutes
        elif len(t) == 3:
            times.append((int(t[0])*60.0) + int(t[1]) + (int(t[2])/60.0))
        else:
            raise Exception('Could not extract times')
    except ValueError as error:
        print('Could not extract times')
        print(error)


#Crop xrange to selected time
if crop is not None:
    print('Selecting data from %.1f to %.1f min' % (crop[0], crop[1]))
    #Get indices
    for i in range(len(times)):
        if times[i] >= crop[0]:
            l = i
            break
    for i in range(len(times)):
        if times[i] > crop[1]:
            r = i
            break
    #Crop
    times = times[l:r]
    for i in range(len(values)):
        values[i] = values[i][l:r]


#Need loads of times vectors
full_values = list(values)
full_time = list(times)
times = [list(times) for i in range(len(values))]

#Flag for coloring the border of the well if it was adjusted
adjusted = [0 for i in range(len(values))]
#0: no adjust. 1: Jump. 2: Slope. 3: both.

        
#Calculate linear regressions!
#Full length linear range
if full_flag:
    regress = []
    for i in range(len(values)):
        reg = sp_stats.linregress(times[i], values[i])
        regress.append(reg)
    print('Full length regressions calculated')
else:
    #Detect jumps in data
    c = 0
    for i in range(len(values)):
        #Get average point to point _variation_
        #(Just average increase doesn't work because data might be unsteady but even)
        av = sum([abs(values[i][j+1] - values[i][j]) for j in range(len(values[i])-1)])
        av = av / (len(values[i]) - 1)
        
        #Check if any datapoint is much different from av
        breakp = []
        for j in range(1,len(values[i])):
            d = abs(values[i][j] - values[i][j-1])  #Difference with next datapoint
            if d > break_thresh * av:
                breakp.append(j)
                
        if len(breakp) > 0:
            adjusted[i] = 1
            c += 1
            
            #Split array by breakp (from np arrray_split)
            Ntotal = len(values[i])
            Nsections = len(breakp) + 1
            div_points = [0] + breakp + [Ntotal]
            
            split = []
            for n in range(Nsections):
                st = div_points[n]
                end = div_points[n + 1]
                split.append([st, end])
            #Get longest section
            ls = [s[1] - s[0] for s in split]
            split = split[np.argmax(ls)]
            v = values[i][split[0]:split[1]]
            t = times[i][split[0]:split[1]]
            times[i] = t
            values[i] = v
            #Better take first if > min_n, otherwise take second?
    
    print("%d out of %d sets optimized for breaks" % (c, len(values)))    
    
    
    #Calculate linear range automatically!
    #if r2 is not improved above threshold then probably dead variant so full length taken
    c = 0
    regress = []
    for j in range(len(values)):
        for i in reversed(range(min_measure, len(values[j]))):
            t_reg = sp_stats.linregress(times[j][0:i],values[j][0:i])
            r2 = t_reg.rvalue ** 2
            if r2 > min_r2:
                regress.append(t_reg)
                if i < (len(values[j])-1):
                    adjusted[j] += 2
                    times[j] = times[j][:i]
                    values[j] = values[j][:i]
                    c += 1
                break
            if i == min_measure:
                regress.append(sp_stats.linregress(times[j],values[j]))
                
    print("%d out of %d sets optimized for linear range\n" % (c, len(values)))

#Scale for more interpretable numbers
slopes = [l.slope * scale_factor for l in regress]


#VISUALIZATION
sns.set_style("white")
sns.set_context("talk")
#Get colors
palette = sns.color_palette("Reds", 100)  #100 Colors, from light to dark red (Try YlOrRd??
palette[0] = "white"
scaled = []
colors = []
for i in range(len(slopes)):
    scaled.append(int(round(slopes[i] / max(slopes) * 99))) #Scale between 0 and 99 for coloration!
    if scaled[i] >= 0:
        colors.append(palette[scaled[i]])    #Add the corresponding color from the palette
    else:
        colors.append("white")    #White for negative values
        
#Plot
f, axarr = plt.subplots(8, 12, sharex=True, sharey=True, figsize=(12,8))
f.subplots_adjust(hspace = 0.5, wspace = 0.2)
spine_color = ["#636363", "#66c2a5", "#8da0cb", "#fc8d62"]  #gray, green, blue, yellow
               
#Label row and column
MTP_row = ["A","B","C","D","E","F","G","H"]
for i in range(8):
    f.text(0.09 ,(0.845 - 0.1 * i), MTP_row[i], transform=f.transFigure, va='center', ha='center')
for i in range(12):
    f.text((0.15 + 0.0655 * i), 0.93, str(i+1), transform=f.transFigure, va='center', ha='center')
    
#Plot 96 MTP
i = 0
for a in range(8):
    for b in range(12):
        axarr[a, b].plot(full_time, full_values[i], marker="o", 
             lw=0, ms=1, mfc="gray", mec="gray")
        axarr[a, b].plot(times[i], values[i], marker="o",
             lw=0, ms=3, mfc="g", mec="k")
        xreg = [times[i][0], times[i][-1]]
        yreg = [regress[i].slope * x + regress[i].intercept for x in xreg]
        axarr[a, b].plot(xreg, yreg, color="k", lw=1)
        axarr[a, b].set_title("%.2f" % (slopes[i]), y=0.9, fontsize=12)
        axarr[a, b].set_facecolor(colors[i])
        axarr[a, b].set_yticklabels([])
        axarr[a, b].set_xticklabels([])
        for spine in axarr[a, b].spines.values():
            spine.set_edgecolor(spine_color[adjusted[i]])  
            spine.set_linewidth(1.5)
        i = i + 1


#File output
if out_flag:
    plt.savefig('out_'+name+'.pdf', bbox_inches='tight')
    out = open('out_'+name+'.txt', "w")
    out.write("well \t slope \t r_squared \n")
    for i in range(len(values)):
        out.write(wells[i] + "\t" + str(slopes[i]) +"\t" + 
                  str(regress[i].rvalue**2) + "\n")
    out.close()

#Plot single graph for decision on range for regression   
if singles_flag:
    plt.figure(figsize=(10,6))
    pl = convert_well(singles)
    i = 0
    for e in pl:
        xreg = [times[e][0], times[e][-1]]
        yreg = [regress[e].slope * x + regress[e].intercept for x in xreg] 
        plt.plot(full_time, full_values[e], marker="o", ms=4, 
                 lw=0, mfc='gray', mec='gray')
        plt.plot(times[e], values[e], marker="o", ms=7, lw=0, 
                 label=singles[pl.index(e)])
        plt.plot(xreg, yreg, color="k",
                 label=str(round(slopes[e],1))+" / "
                 +str(round(regress[e].rvalue**2,3)))
        print("R2 of well {}: {}".format(singles[pl.index(e)],
                                         regress[e].rvalue**2))
        i += 1
    plt.legend()
    plt.ylabel("Absorbance")
    plt.xlabel("Time [min]")


plt.show()

