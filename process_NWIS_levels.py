# program to QC groundwater levels from NWIS database

# Required Inputs:
# NWIS Site File (with station info; tab delimited)
# NWIS groundwater levels file (with levels and dates; tab delimited)
# required NWIS fields: agency_cd, site_no, alt_va, alt_acy_va, lev_dt, lev_va, lev_status_cd, well_depth_va, qw_count_nu, reliability_cd
#
# Coordinates File (with coordinates in model coordinate system for each well)
# This can be created by exporting csv of stations and WGS84 coordinates into ArcMap
# A character needs to be added to ends of station numbers (e.g. '434238088592501n') so that Arc treats them as strings
# Coordinates exported from Arc should have columns site_no2,POINT_X,POINT_Y
#
# The sorting algorithm below might be cleaner and more transparent if reformulated in terms of cumulative error
# e.g. 5 ft. + 3.5 ft. for alt_accuracy of <=5 ft., and std deviation in measurements of 3.5
# would still have to come up with quantities to represent the lower perceived error for wells with Water Quality measurements, artesian, etc.

import numpy as np
from collections import defaultdict
import datetime as dt
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.backends.backend_pdf import PdfPages

# Input files
infofile='Columbia_header.txt'
levelsfile='Columbia_levels.txt'
coordsfile= 'Columbia_NWIS_headsWTM.csv' 

# Outfiles
pdffile='extended_records.pdf'

# Settings
mode='GFLOW' # GFLOW or MODFLOW; writes either a tp file, or .hob file for MF2k observation process
discard_lev_status_cd=['P','R','S','T','Z'] # list of level_status_cds to discard (e.g. if well was being pumped)
screen_length=50 # assumed open interval length- NWIS only has well bottoms, so this will be added to bottom to get an open interval for the target. This is also important for avoiding wells being placed in thin (e.g. 1 ft.) "dummy" layers.


print "getting well info, water levels, and coordinates..."

# get header info

def getheader(filename,indicator,delimiter):
    headervar=0
    infile=open(filename,'r').readlines()
    for line in infile:
        cline=line.strip().split(delimiter)
        if cline[0]=='agency_cd':
            break
        else:
            headervar+=1
    return(headervar)

info_header=getheader(infofile,'agency_cd','\t')
levels_header=getheader(levelsfile,'agency_cd','\t')
        
info=np.genfromtxt(infofile,delimiter='\t',skiprows=info_header,names=True,dtype=None)[1:]
levelsdata=np.genfromtxt(levelsfile,delimiter='\t',skiprows=levels_header,names=True,dtype=None)[1:]

wells=np.unique(levelsdata['site_no'])

# build dictionaries of levels,dates and codes by USGS well no.

levels=defaultdict(list)
dates=defaultdict(list)
codes=defaultdict(list)
WellDepth_elev=defaultdict()
pumping_count=defaultdict(list)

discarded=open('discarded_wells.txt','w')
# Go through each level in levels file
for line in levelsdata:
    wellnum=line['site_no']
    status=line['lev_status_cd']
    
    # toss out water levels with bad status code (see above)
    discard=False
    for code in discard_lev_status_cd:
        if status==code:
            discard=True
            pumping_count[wellnum].append(1) # keep track of number of wells and levels influenced by pumping
            break
        
    if discard:
        continue # skip level
        
    info_ind=np.where(info['site_no']==wellnum)[0][0] 
    elevation=float(info['alt_va'][info_ind].strip())
    try:
        date=dt.datetime.strptime(line['lev_dt'],'%Y-%m-%d')
    except ValueError:
        continue
    level=line['lev_va']
    try:
        level=elevation-float(level)
    except ValueError:
        level=None
    try:
        welldepth_elev=elevation-float(info['well_depth_va'][info_ind].strip())   
    except ValueError:
        welldepth_elev=None
    dates[wellnum].append(date)
    levels[wellnum].append(level)
    codes[wellnum].append(status)
    WellDepth_elev[wellnum]=welldepth_elev

pumping_discarded_levels=len([item for sublist in pumping_count.values() for item in sublist])
print "Discarded %s pumping-influenced levels from %s wells..." %(pumping_discarded_levels,len(pumping_count.keys()))
discarded.write('well,n_discarded\n')
for key in pumping_count.iterkeys():
    discarded.write('%s,%s\n' %(key,len(pumping_count[key])))

# get coordinates from coordsfile
coordsdata=np.genfromtxt(coordsfile,delimiter=',',names=True,dtype=None)
coords=defaultdict(list)
for well in coordsdata:
    site_no=well['site_no2'][1:-2] # had to modify original site_nos by adding a "n", so that Arc would treat as string
    coords[site_no]=[well['POINT_X'],well['POINT_Y']]

if len(coords)<>len(wells):
    raise Warning("Number of coordinates does not match number of wells!")
    
print "sorting wells based on QC criteria..."
names=defaultdict(list)
rejects=defaultdict(list)
wells2plot=[]

# writeout information on "poor" wells that didn't meet any of the quality criteria
discarded.write('\n\nwell,num_measurements,reliability_code,alt_accuracy\n')
num_artesian=0
for well in wells:    
    print well
    
    # reset variables
    name=None
    dateslist=[]
    maxmin=None
    n=0
    artesian=False
    
    # get list of measurement dates
    # toss wells with no dates
    dateslist=dates[well]
    if len(dateslist)==0:
        continue    
    
    # get info from header file
    info_ind=np.where(info['site_no']==well)[0][0]
    alt_acc=float(info['alt_acy_va'][info_ind].strip())
    numWQ=float(info['qw_count_nu'][info_ind].strip())
    Drely=info['reliability_cd'][info_ind]
    elevation=float(info['alt_va'][info_ind].strip())
    
    # before sorting, identify wells with more than one measurement
    # extract single value from list for wells with one measurement
    # for artesian wells, set GW elevation to wellhead elevation if no value
    
    if levels[well][0]==None:
        if 'F' in codes[well]:
            artesian=True
            levels[well]=elevation
            num_artesian+=1
        elif 'E' in codes[well]:
            artesian=True
            levels[well]=elevation
            num_artesian+=1        
    elif len(levels[well])>1:
        wells2plot.append(well)
        maxmin=np.max(levels[well])-np.min(levels[well])
    elif len(levels[well])==1:
        try:
            levels[well]=levels[well][0]
        except TypeError: # if no level, value obtained from file might be ''
            if 'F' in codes[well]:
                artesian=True
                levels[well]=elevation
                num_artesian+=1
            elif 'E' in codes[well]:
                artesian=True
                levels[well]=elevation
                num_artesian+=1
                
    # keep only wells that have been checked by a reporting agency
    if Drely=='C':
        
        # get water levels and std
        WLs=levels[well]
        try:
            std=np.std(WLs)
        except TypeError: # means there is only one measurement
            std=None
            
        # get water levels after 1970
        after1970_inds=list(np.where(np.array(dateslist)>dt.datetime(1970,1,1,0,0))[0])
        n=len(after1970_inds)
        
        # function to assign weight category based on standard dev
        def assignby_std(std,name,tier):
            # tier =1,2 (best,good); sets a ceiling on the highest name to assign
            if tier==1 and std<=5:
                name=name+'_best'
            elif tier>=1 and std<=10:
                name=name+'_good'
            else: # std is >10
                name=name+'_fair'
            return name
        
        # if there are measurements post1970
        if n>0:
            
            # cull levels to post1970 (just for weight category assignment)
            try:
                WLs=[i for j, i in enumerate(levels[well]) if j in after1970_inds]
                std=np.std(WLs)
            except TypeError: # means there is only one measurement
                std=None  
            
            # more than two measurements after 1970
            if n>2 and alt_acc<=10:
                
                if alt_acc<=5 and n>2:
                    name=assignby_std(std,well[5:],1)
                elif alt_acc<=10 and n>30:
                    name=assignby_std(std,well[5:],1)
                elif alt_acc<=10 and n>2:
                    name=assignby_std(std,well[5:],2)
                else:
                    name=well[5:]+'_fair'
            
            # Just two
            elif n==2 and alt_acc<=5:
                name=assignby_std(std,well[5:],2)
                
            # Only one measurement post 1970
            elif n==1 and alt_acc<5:
                name=well[5:]+'_good'
            elif n==1 and numWQ>0 and alt_acc<=5:
                name=well[5:]+'_good'
            elif n==1 and artesian and alt_acc<=5:
                name=well[5:]+'_good'
            elif n==1 and alt_acc<=10: 
                name=well[5:]+'_fair'
            else:
                name=well[5:]+'_poor'
                rejects[well].append([n,Drely,alt_acc,name])
                discarded.write('%s,%s,%s,%s,%s\n' %(well,n,Drely,alt_acc,name))            
        
        # No measurements after 1970
        else:
            n=len(dateslist)
            
            # if >2 measurements and good alt accuracy
            if n>2 and alt_acc<=5:
                name=assignby_std(std,well[5:],2)
            elif alt_acc<=10:
                name=well[5:]+'_fair'
                
            # only one pre1970 measurement or poor alt accuracy
            else:
                name=well[5:]+'_poor'
                rejects[well].append([n,Drely,alt_acc,name])
                discarded.write('%s,%s,%s,%s,%s\n' %(well,n,Drely,alt_acc,name))
                
    # measurement has not been checked by reporting agency
    # keep only post 1970, and then only wells with 2+ measurements, std<=10, and alt accuracy <=10
    
    elif Drely=='U':
        
        after1970_inds=list(np.where(np.array(dateslist)>dt.datetime(1970,1,1,0,0))[0])
        n=len(after1970_inds)        
        try:
            WLs=[i for j, i in enumerate(levels[well]) if j in after1970_inds]
            std=np.std(WLs)
        except TypeError: # means there is only one measurement
            std=None
            
        if n>2 and alt_acc<=5:
            if std<=5:
                name=well[5:]+'_good'
        elif n>30 and alt_acc<=10:
            if std<=5:
                name=well[5:]+'_good'
        elif n>2 and alt_acc<=10:
            if std<=10:
                name=well[5:]+'_fair'
        else:
            name=well[5:]+'_poor'
            rejects[well].append([n,Drely,alt_acc,name])
            discarded.write('%s,%s,%s,%s,%s\n' %(well,n,Drely,alt_acc,name))            
            
    # discarded because not reliable
    else:
        name=well[5:]+'_poor'
        rejects[well].append([n,Drely,alt_acc,name])
        discarded.write('%s,%s,%s,%s,%s\n' %(well,n,Drely,alt_acc,name))
    names[well]=name
discarded.close()

# all names are assigned; build list and check for duplicates
# if duplicates, use different 10-digit string of site_no
# continue, changing 10 digit string, until no duplicates
# reduce digits taken from site number if necessary
nameslist=[]

dupscount=0
print "modifying any duplicate names by choosing new 10-digit strings from site numbers..."
print "(duplicate names occur because only a 10-digit subset of the 15-digit ID is used for names;\nsome programs like GFLOW or PEST have observation name length limits)"

for name in names.itervalues():
    nameslist.append(name)
num_dups=len(nameslist)-len(np.unique(np.array(nameslist)))
i=0
while num_dups>0:
    dupscount=0
    i+=1
    wells2delete=[]
    for well in names.iterkeys():
        if len(names[well])==0:
            wells2delete.append(well)
            continue
        count=nameslist.count(names[well])
        if count>1:
            dupscount+=1
            oldname=names[well]
            if i<=4:
                names[well]=well[4-i:(-1-i)]+oldname[-5:]
                print names[well]
            else:
                names[well]=well[i:]+oldname[-5:]
    print "fixed %s duplicate names!" %(dupscount)
    for well in wells2delete:
        print 'deleting well: %s\n' %(well)
        del names[well]
    nameslist=[]
    for name in names.itervalues():
        nameslist.append(name)    
    num_dups=len(nameslist)-len(np.unique(np.array(nameslist)))
    if num_dups>0:
        print "still duplicates; trying new 10-digit strings from site numbers.."
    else:
        break
print "done fixing duplicates!"

# For wells with multiple measurements, calculate average values for post and pre-1970
# Plot out measurements and average values for comparison
# replace multiple levels with new average level
pdf=PdfPages(pdffile)

print 'calculating average values; plotting levels for wells with multiple levels:'
for well in wells2plot:
    
    print well
    info_ind=np.where(info['site_no']==well)[0][0]
    alt_acc=float(info['alt_acy_va'][info_ind].strip())
    numWQ=float(info['qw_count_nu'][info_ind].strip())
    codez=[c for c in codes[well] if c<>'']
    dates2plot=list(mdates.date2num(dates[well]))
    WLs=levels[well]
    use_post1970=True
    
    if len(dates2plot)<>len(WLs):
        raise ValueError("Dates and WLs for %s have different lengths!" %(well))
    
    cutoff=mdates.date2num(dt.datetime(1970,1,1,0,0))
    
    # post 1970 average
    post1970=[d for d in dates2plot if d>cutoff]
    num2skip=len(WLs)-len(post1970)
    post1970WLs=WLs[num2skip:]
    avg_post1970=[np.mean(post1970WLs)]*len(post1970)    
    std_post1970=np.std(post1970WLs)
    
    # pre 1970 average
    pre1970=[d for d in dates2plot if d<cutoff]
    pre1970WLs=WLs[:len(pre1970)]
    avg_pre1970=[np.mean(pre1970WLs)]*len(pre1970)
    std_pre1970=np.std(pre1970WLs)
    
    # decide which set of water levels to use
    # replace multiple values in levels file with average
    if len(post1970)==0:
        use_post1970=False
        plot_title='WLs in well ' + names[well] + '; using pre-1970 average'
        levels[well]=np.mean(pre1970WLs)
    if use_post1970:
        plot_title='WLs in well ' + names[well] + '; using post-1970 average'
        levels[well]=np.mean(post1970WLs)
        
    fig, ax1 = plt.subplots(1,1,sharex=True,sharey=False)
    ax1.grid(True)
    p1=ax1.plot_date(dates2plot,WLs,'bo',label='WLs')
    plt.xticks(rotation=45,fontsize=10)
    plt.yticks(fontsize=10)
    # get axis limits
    # set plotting range for average values
    x1,x2,y1,y2=ax1.axis()
    if x1>cutoff:
        start=x1
    else:
        start=cutoff
    if x2>cutoff:
        end=cutoff
    else:
        end=x2
    
    # plot average values
    if len(post1970)>0:
        #ax1.axhline(y=avg_post1970, xmin=start, color='r',label='post-1970 avg')
        p2=ax1.plot_date(post1970,avg_post1970,'r',label='post-1970 avg')
        ax1.text(0.95,0.05,'stdev: %s\nalt accuracy: %s\nnum WQ: %s\nWell codes: %s' %(round(std_post1970,2),alt_acc,numWQ,codez),verticalalignment='bottom',horizontalalignment='right',transform=ax1.transAxes)    
    if len(pre1970)>0:
        #ax1.axhline(y=avg_pre1970, xmax=end, color='g',label='pre-1970 avg')
        p3=ax1.plot_date(pre1970,avg_pre1970,'g',label='pre-1970 avg')
        ax1.text(0.05,0.05,'stdev: %s\nalt accuracy: %s\nnum WQ: %s\nWell codes: %s' %(round(std_pre1970,2),alt_acc,numWQ,codez),verticalalignment='bottom',horizontalalignment='left',transform=ax1.transAxes)
    handles, labels = ax1.get_legend_handles_labels()
    ax1.legend(handles,labels)
    ax1.set_title(plot_title)
    pdf.savefig()
pdf.close()
plt.close('all')
print "Done plotting, see %s for results" %(pdffile)

print "writing testpoint files..."

if mode=='GFLOW':
    ofp2=open('heads_ALL.tp','w')
    for category in ['best','good','fair','poor']:
        # create tp file for each
        fname=category+'_heads.tp'
        print fname
        ofp=open(fname,'w')
        for well in wells:
            if category in names[well]: # note: for some reason, there were a handful of well for Columbia that remained unnamed; unclear why. Looking for the category weeds them out, but doesn't explain why they weren't capture by the last "else" statement in the sorting algorithm.
                ofp.write('%s,%s,%s,0,Piezometer,%s\n' %(coords[well][0],coords[well][1],levels[well],names[well]))
                ofp2.write('%s,%s,%s,0,Piezometer,%s\n' %(coords[well][0],coords[well][1],levels[well],names[well]))
        ofp.close()
    ofp2.close()    
    print 'heads_ALL.tp'    

if mode=='MODFLOW':
    ofp=open('NWIS_MFhobs_export.csv','w')
    print 'NWIS_MFhobs_export.csv'
    ofp.write('Name,POINT_X,POINT_Y,WL,ScreenTop,ScreenBot\n')
    for category in ['best','good','fair','poor']:
        for well in wells:
            if category in names[well]:
                if WellDepth_elev[well]==None:
                    continue
                else:
                    ofp.write('%s,%s,%s,%s,%s,%s\n' %(names[well],coords[well][0],coords[well][1],levels[well],WellDepth_elev[well]+screen_length,WellDepth_elev[well]))
    ofp.close()
print "finished OK"