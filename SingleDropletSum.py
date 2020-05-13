#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Reads CSV of single droplet data and integrates droplet signals
Expects at least 1 droplets signal, and a 'background' ion signal

@author: meganwillis

"""

################################
import numpy as np  
import pandas as pd                      
import matplotlib.pyplot as plt  
################################

### -- dictionary of inputs
inputs_dict = {'path': '/Users/meganwillis/GoogleDriveLBL/Data/20190917/',
               'file': '20190917_GD_30um_soln2_temprange.csv',
               't_start':175,
               't_stop':180,
               'dropletspecies': ['maleic','malonic'],
               'backgroundspecies': ['lactic'],
               'droplet_mz': 'malonic',
               'todrop': [],
               'toratio': ['maleic','malonic’],
               'scanrate': 8.5,
               'lag': 60,
               'threshold':8,
               'influence':0.0009,
               'npoint_low': 10,
               'npoint_high': 150,
               'saveevents': False,
               'events_outfname':'20190917_GD_temprange_events.csv',
               'saveeventstats': False,
               'eventstats_outfname':'20190917_GD_temprange_eventstats.csv'
        }

#########################
#####--INPUTS:--#########
#########################
### -- full path to file
p_file = inputs_dict['path']+inputs_dict['file']
### --- time range in the file to look at
t_start = inputs_dict['t_start']
t_stop = inputs_dict['t_stop']
### --- select the species to look at in droplet events
dropletspecies = inputs_dict['dropletspecies']
background = inputs_dict['backgroundspecies']
### --- choose one species to use as an indictor of droplet events (least volatile is best)
droplet_mz = inputs_dict['droplet_mz']
### --- row indices of droplet events to ignore, zero indexed (leave empty list if none)
todrop = inputs_dict['todrop']
### --- select two species to ratio as a diagnostic (i.e., toratio[0]/toratio[1])
toratio = inputs_dict['toratio']
### --- MS scan rate to calculate droplet event duration, units of 1/s
scanrate = inputs_dict['scanrate']
### ---thresholding algorithm parameters
lag = inputs_dict['lag']
threshold = inputs_dict['threshold']
influence = inputs_dict['influence']
### --- limit the number of points in a 'good' droplet event
npoint_low = inputs_dict['npoint_low'] 
npoint_high = inputs_dict['npoint_high'] 
### -- save someoutputs
saveevents = inputs_dict['saveevents']  #save information on all individaul events in the file
outfile0 = inputs_dict['path'] + inputs_dict['events_outfname'] 
saveeventstats = inputs_dict['saveeventstats'] #save statistics across the events
outfile1 = inputs_dict['path'] + inputs_dict['eventstats_outfname']

#########################
###---FUNCTIONs--:#######
#########################

def thresholding_algo(y, lag, threshold, influence):
    # Implementation of algorithm from https://stackoverflow.com/a/22640362/6029703
    # https://stackoverflow.com/questions/22583391/peak-signal-detection-in-realtime-timeseries-data/43512887#43512887
    '''Uses a running mean of the baseline to detect peaks that are a specified number of standard deviations above the baseline'''
    signals = np.zeros(len(y))
    filteredY = np.array(y)
    avgFilter = [0]*len(y)
    stdFilter = [0]*len(y)
    avgFilter[lag - 1] = np.mean(y[0:lag])
    stdFilter[lag - 1] = np.std(y[0:lag])
    for i in range(lag, len(y)):
        if abs(y[i] - avgFilter[i-1]) > threshold * stdFilter [i-1]:
            if y[i] > avgFilter[i-1]:
                signals[i] = 1
            else:
                signals[i] = -1

            filteredY[i] = influence * y[i] + (1 - influence) * filteredY[i-1]
            avgFilter[i] = np.mean(filteredY[(i-lag+1):i+1])
            stdFilter[i] = np.std(filteredY[(i-lag+1):i+1])
        else:
            signals[i] = 0
            filteredY[i] = y[i]
            avgFilter[i] = np.mean(filteredY[(i-lag+1):i+1])
            stdFilter[i] = np.std(filteredY[(i-lag+1):i+1])

    return dict(signals = np.asarray(signals),
                avgFilter = np.asarray(avgFilter),
                stdFilter = np.asarray(stdFilter))

def plotoutput():
    fig, ax1 = plt.subplots(figsize=(16,6))
    for item in dropletspecies:
        ax1.plot(data.Time_min, data[item])
    ax1.plot(data.Time_min, result.avgFilter, color="darkgrey")
    ax1.fill_between(data.Time_min,result.avgFilter-result.stdFilter, result.avgFilter+result.stdFilter, alpha=0.5, facecolor='darkgrey', edgecolor='None')
    plt.ylim([-1e5,data[item].max()/3])
    ax1.set_xlabel('Time (minutes)', fontsize = 14)
    ax1.set_ylabel('Signal Intensity (Hz)', fontsize = 14)
    ax2=ax1.twinx()
    ax2.plot(data.Time_min, data.signals, marker=".", color="black", markersize=1)
    ax2.plot(eventsumdf.Time_min, eventsumdf.ratio, marker=".", color="blue", linestyle="None", markersize=10)
    #plt.xlim([1,20])
    ax2.set_ylabel('Ratio of peak areas', fontsize = 14)
    
    plt.show()

#########################
####----MAIN:----########
#########################

##### --- load data
df = pd.read_csv(p_file)
df_sub = df.loc[(df.Time_min >= t_start) & (df.Time_min <= t_stop)] # slice out a time range
data = df_sub.copy() # create a copy so we don't modify original df

#### --- choose species and sum up fragment signals
species = dropletspecies + background
if 'lactic' in species:
    data['lactic'] = data.mz89
if 'bisulfate' in species:
    data['bisulfate'] = data.mz97
if 'maleic' in species:
    data['maleic'] = data.mz115 #+ data.mz71
if 'maleic_frag' in species:
    data['maleic_frag'] = data.mz71
if 'succinic' in species:
    data['succinic'] = data.mz117 #+ data.mz73
if 'malonic' in species:
    data['malonic'] = data.mz103 #+ data.mz59
if 'malonic_frag' in species:
    data['malonic_frag'] = data.mz59
if 'tartaric' in species:
    data['tartaric'] = data.mz149
if 'citric' in species:
    data['citric'] = data.mz111 + data.mz173 + data.mz191


#### -- run threshold algorithm  -- inputs= y, lag, threshold, influence
result = pd.DataFrame.from_dict(thresholding_algo(np.array(data[droplet_mz]), lag, threshold, influence)) 
## -- notes on input values:
#small changes in threshold std can give quite different results in some cases, make sure you check the impact

### -- add thresholding algo results to data array
signals = np.array(result.signals)
data['signals'] = signals

### --- integrate droplet events
outcols = ['Time_min','npoints'] + species # declare column names for output data
datacolumns = species # specify only the data containing columns to loop over later
tempdf = pd.DataFrame(columns=outcols) # initialize empty temporary df
eventsumdf = pd.DataFrame(columns=outcols, dtype='float64') # initialize df to store final data (peak areas, event start time and number of mass spectra)

for idx, row in data.iterrows(): #loop over all the rows in the data array

    if row.signals == 1: # select the rows where we have detected a droplet event
        tempdf = tempdf.append(data.loc[idx].copy(deep=True)) # add a row to our temporary df

    else :
            if (len(tempdf) > npoint_low) & (len(tempdf) < npoint_high): # limit an 'event' to have a certain number of data points, helps with noisy data
                print(len(tempdf)) # print the number of data points in the droplet event
                temp = np.empty(len(outcols)) # create an empty row the sie of outcols
                temp[:]=np.nan # fill it with np.nan
                emptyrow = pd.Series(temp, index=eventsumdf.columns) # convert to a pandas Series so we can append it easily to a pandas df
                eventsumdf = eventsumdf.append(emptyrow, ignore_index=True) # append to eventsumdf
                
                for item in datacolumns:
                    eventsumdf.iloc[-1, eventsumdf.columns.get_loc(item)] = tempdf[item].sum() # loop through datacolumns and sum signals
                
                eventsumdf.iloc[-1, eventsumdf.columns.get_loc('Time_min')] = tempdf['Time_min'].iloc[0] # set the start time of the droplet event
                eventsumdf.iloc[-1, eventsumdf.columns.get_loc('npoints')] = len(tempdf) # record the number of points in the droplet evet
                
                tempdf = pd.DataFrame(columns=outcols) # empty tempdf
            
            else: #if the number of point is smaller or larger than a specific value, ignore it
                tempdf = pd.DataFrame(columns=outcols) # empty tempdf
# this is assuming that time ranges always end with some baseline and not with a droplet event

### -- calculate the ratio of droplet peak areas
eventsumdf['ratio'] = eventsumdf[toratio[0]]/eventsumdf[toratio[1]]

### -- remove specific droplet events (e.g., those that were multiple droplets)
if len(todrop) > 0:
    eventsumdf.drop(todrop, inplace=True)

### -- calculate stats across all good droplet events
event_stats = eventsumdf.describe() 

### -- print some genral outputs
print('Number of droplet events =', eventsumdf.shape[0])
print('Mean points/event = ', "{:.1f}".format(event_stats.npoints['mean']), '+/-', "{:.1f}".format(event_stats.npoints['std']))
print('Mean droplet event time = ',"{:.1f}".format((1/scanrate)*event_stats.npoints['mean']), '+/-', "{:.1f}".format((1/scanrate)*event_stats.npoints['std']),'seconds')

### -- print stats on each item in species list
for item in species:
    #print('Mean Area ', item, ' = ', event_stats[item]['mean'], '+/-', event_stats[item]['std'])
    print('Peak Area RSD for', item, ' = ', "{:.1f}".format((event_stats[item]['std']/event_stats[item]['mean'])*100), '%')

### -- print info on the calculated peak area ratio
print('Peak Area Ratio RSD', toratio[0],'/',toratio[1],'=',"{:.1f}".format((event_stats['ratio']['std']/event_stats['ratio']['mean'])*100), '%')
print('Peak Area Ratio Mean', toratio[0],'/',toratio[1],'=',"{:.2f}".format(event_stats['ratio']['mean']))

### -- plot the output
plotoutput()

#### -- save event_stats as a CSV file
if saveevents == True:
    eventsumdf.to_csv(outfile0)
if saveeventstats == True:
    event_stats.to_csv(outfile1)
    

### -- save figure
#outpath = '/Users/meganwillis/GoogleDriveLBL/Data/20190528/'
#outfile = outpath + 'DropletSignalDetection.jpg'
#plt.savefig(outfile, format='jpg', dpi=600)

