#!/usr/bin/python

import os, sys
import glob
import subprocess
import math
import re
import numpy
from pandas import *
from scipy import stats
from scipy.stats import ttest_ind

import rpy2.robjects as robjects
import pandas.rpy.common as com
from rpy2.robjects.packages import importr
from rpy2.robjects.lib import grid
from rpy2.robjects.lib import ggplot2

grdevices = importr('grDevices')

mytheme = {
            'panel.background':ggplot2.element_rect(fill='white',colour='white'),
            'axis.text':ggplot2.element_text(colour="black",size=15),
            'axis.title':ggplot2.element_text(colour="black",size=15),
            'plot.title':ggplot2.element_text(face="bold", size=20,colour="black"),
            'strip.text.y':ggplot2.element_text(colour="black",face="bold",size=15),
            'strip.text.x':ggplot2.element_text(colour="black",face="bold",size=15)
            }
            #'panel.grid.major':ggplot2.theme_line(colour = "grey90"),
            #'panel.grid.minor':ggplot2.theme_line(colour = "NA"),
            #'axis.line':ggplot2.theme_line(size = 1.2, colour="black"),
            #'axis.ticks':ggplot2.theme_line(colour="black"),

################################################################################
# filesInDir
################################################################################
def filesInDir( path, filetype, prefix=None ):
    filelist = []
    if prefix is not None :
        filelist = glob.glob( os.path.join( path, "%s*%s" %(prefix, filetype)))
    else :
        filelist = glob.glob( os.path.join( path, "*."+filetype ) )
    if filetype[-5:] == "fastq" :
        prefixdict = {}
        for file in sorted(filelist) :
            prefix = getBasename( file )
            if not prefixdict.has_key(prefix[:-3]) :
                prefixdict[prefix[:-3]] = []
            prefixdict[prefix[:-3]].append( file )
        return prefixdict.values()
    return filelist
# End filesInDi

def filesInDirRecursive( rootdir, filetype ):
    fileList = []
    for root, subFolders, files in os.walk(rootdir):
        for cfile in files:
            fname, fext = os.path.splitext(cfile)
            if filetype == fext :
                fileList.append(os.path.join(root,cfile))

    filedata = []
    for cfile in fileList :
        filepath, basename, ext = getBasename( cfile )
        cleanbase = re.sub("\W+","_",basename)
        #print basename, cleanbase
        drugmatch = re.match( "(\w+)-(\d)[ _]+(.*)", basename )
        regmatch = re.match( "(\w+)-(\d)[ _]b[ _](.*)", basename )
        if regmatch is not None :
            fishtype,rep,dosage = regmatch.groups()
            drug = "none"
        elif drugmatch is not None :
            fishtype,rep,dosage = drugmatch.groups()
            drug = "drug"
        else : print "Error: failed",basename; continue
        filedata.append( Series(data=[filepath, cfile, fishtype, rep, dosage, drug],
                                index=["filepath","filename","fishtype",
                                       "sample","dosage","drug"]))
    filedata = DataFrame(filedata)
    filedata.loc[filedata.fishtype == "ctl","fishtype"] = "ctrl"
    filedata.loc[filedata.fishtype == "mol","fishtype"] = "mo"
    #filedata.loc[filedata.dosage == "whole brain","dosage"] = "whole"
    #filedata["dosage"] = [x.replace(" ","") for x in filedata.dosage]
    return filedata
# END filesInDirRecursive

################################################################################
# getBasename
################################################################################
def getBasename( originalfilename ) :
    filepath, filename = os.path.split( originalfilename )
    #print "path",filepath,"name", filename
    if originalfilename[-3:] == ".gz" :
        basename, suffix = os.path.splitext( filename[:-3] )
        #print "Base:", basename, "Suff:",suffix
    else :
        basename, suffix = os.path.splitext( filename )
    return filepath, basename, suffix
# End getBasename

def fixRLevels( r_dataframe, column, tlevels ):
    replace = robjects.FactorVector( r_dataframe.rx2(column), 
                                    levels=robjects.StrVector(tlevels) )
    allcolumns = list(r_dataframe.colnames)
    allcolumns = [ x+"_old" if x == column else x for x in allcolumns]
    new_r_df = r_dataframe.cbind(replace)
    new_r_df.colnames = robjects.StrVector(allcolumns+[column])
    return new_r_df
# END fixRLevels

################################################################################
# smooth
#,window='hanning'
#if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
#raise ValueError, "Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"
#if window == 'flat': #moving average
    #w=numpy.ones(window_len,'d')
################################################################################
def smooth(tlist,window_len=11):
    # This is a hack.... Sorry about that
    window_len = 100 if tlist.size < 800 else window_len
    #print "Using window_len:", window_len
    if tlist.ndim != 1:
        raise ValueError, "smooth only accepts 1 dimension arrays."
    if tlist.size < window_len:
        print "list size:",tlist.size, "Window:",window_len
        raise ValueError, "Input vector needs to be bigger than window size."
    if window_len<3:
        return tlist
    s=numpy.r_[2*tlist[0]-tlist[window_len-1::-1],tlist,2*tlist[-1]-tlist[-1:-window_len:-1]]
    #w=numpy.hanning(window_len)
    w=numpy.ones(window_len,'d')
    #y=numpy.convolve(w/w.sum(),s,mode='same')
    #return y[window_len:-window_len+1]
    windows = [[min(len(tlist)-window_len/2,max(0,math.floor(i-window_len/2))), 
                min(len(tlist),min(len(tlist),math.floor(i+window_len/2)))] 
               for i in range(0,len(tlist)) ]
    #print windows[:10]
    #print windows[-10:]
    #print [windows[i] for i in range(0,len(tlist))]
    #for lt in [tlist[windows[i][0]:windows[i][1]] for i in range(0,len(tlist))] :
        #print lt
    return [numpy.median(tlist[windows[i][0]:windows[i][1]]) for i in range(0,len(tlist))]
# END smooth

#------------------------------------------------------------------------------#
#                                Calculations
#------------------------------------------------------------------------------#
def windowZscore( x, window_len=50 ) :
    win_zscores = []
    for i in range(len(x)) :
        if i < window_len/2 : win = (0,i,window_len)
        elif (len(x)-i ) < window_len/2 : 
            win = (len(x)-window_len,i-(len(x)-window_len),len(x))
        else : win = (i-window_len/2,window_len/2,i+window_len/2)
        czscore = stats.zscore(x[win[0]:win[2]])
        win_zscores.append(czscore[win[1]])
    return win_zscores
# ENd windowZscore

#------------------------------------------------------------------------------#
#                                   Plots
#------------------------------------------------------------------------------#
def plotTimecourse( timecourse, figdir, fighead ) :
    r_dataframe = com.convert_to_r_dataframe(timecourse)
    p = (ggplot2.ggplot(r_dataframe) +
                ggplot2.aes_string(x="time",y="zscore" ) +
                ggplot2.geom_point(ggplot2.aes_string(colour="factor(peak)")) +
                ggplot2.geom_line() +
                ggplot2.geom_abline(intercept=2, slope=0, linetype="dashed",colour="black") +
                ggplot2.scale_y_continuous("Z-score") +
                ggplot2.scale_x_continuous("Time (seconds)") +
                ggplot2.ggtitle("Timecourse After Normalization\n"+fighead) +
                ggplot2.theme(**mytheme))
                #ggplot2.theme(**mytheme) + \, colour="factor(peak)"
                #ggplot2.theme(**{'axis.text.x': ggplot2.element_text(angle = 45)})  + \
                #ggplot2.geom_boxplot(ggplot2.aes_string(fill="Continent")), notch=True ) + \
    imagefile = "%s/samp%s_timecourse.png" % (figdir,fighead)
    print "Writing file:",imagefile
    grdevices.png(imagefile)
    p.plot()
    grdevices.dev_off()
# END plotTimecourse

def plotRawTimecourse( timecourse, outprefix, fighead ) :
    print "Running plotRawTimecourse"
    print timecourse.head()
    baselinedf = timecourse[["time","smooth"]]
    mean = timecourse.intensity.mean()
    r_dataframe = com.convert_to_r_dataframe(timecourse)
    r_base = com.convert_to_r_dataframe(baselinedf)
    p = (ggplot2.ggplot(r_dataframe) +
                ggplot2.aes_string(x="time",y="intensity", 
                                   colour='factor(dosage)') +
                ggplot2.geom_point() +
                ggplot2.geom_line(colour="black") +
                ggplot2.geom_line(ggplot2.aes_string(x="time",y="smooth"),data=r_base,colour="red") +
                ggplot2.geom_abline(intercept=mean, slope=0, linetype="dashed",colour="black") +
                ggplot2.scale_y_continuous("Raw Intensities") +
                ggplot2.scale_colour_discrete( name="Dosage",
                                      breaks=robjects.StrVector(["baseline", "5mmptz", "62mmptz"]),
                                      labels=robjects.StrVector(["Baseline", "5mM PTZ", "62mM PTZ"])) +
                ggplot2.scale_x_continuous("Time (seconds)") +
                ggplot2.ggtitle("Timecourse Raw data\n"+fighead) +
                ggplot2.theme(**mytheme))
                #,labels=c("Baseline", "5mM PTZ", "62mM PTZ")
                #ggplot2.theme(**{'legend.title':ggplot2.element_blank()}) +
    figname = "%s/%s_rawtimecourse.png" % (outprefix,fighead)
    print "Writing file:",figname
    grdevices.png(figname)
    p.plot()
    grdevices.dev_off()
# END plotRawTimecourse

def plotAllRawTimecourse( timecourse, outprefix, datagrp="all" ) :
    print "Running plotAllRawTimecourse"
    timecourse["fishtype_clean"] = ["Control" if x == "ctrl" else "Morpholino" 
                                   for x in timecourse.fishtype]
    print timecourse.head()
    #baselinedf = timecourse[["time","smooth"]]
    tlevels = [x for x in ["baseline","5mmptz","62mmptz"]
                                  if x in timecourse.dosage.unique() ]
    mean = timecourse.intensity.mean()
    r_dataframe = com.convert_to_r_dataframe(timecourse)
    r_dataframe = fixRLevels(r_dataframe, "dosage", tlevels )
    #r_base = com.convert_to_r_dataframe(baselinedf)
    p = (ggplot2.ggplot(r_dataframe) +
                ggplot2.aes_string(x="time",y="intensity", 
                                   colour='factor(dosage)') +
                ggplot2.geom_point() +
                ggplot2.geom_line(colour="black") +
                ggplot2.scale_y_continuous("Raw Average Intensities (Whole Brain)") +
                ggplot2.scale_x_continuous("Time (seconds)") +
                ggplot2.ggtitle("Intensity Timeseries for all Samples\n("+datagrp+")") +
                ggplot2.scale_colour_discrete( name="Dosage",
                                      breaks=robjects.StrVector(["baseline", "5mmptz", "62mmptz"]),
                                      labels=robjects.StrVector(["Baseline", "5mM PTZ", "62mM PTZ"])) +
                ggplot2.facet_grid(robjects.Formula('sample ~ fishtype_clean')) +
                ggplot2.theme(**mytheme))
                #,labels=c("Baseline", "5mM PTZ", "62mM PTZ")
                #ggplot2.geom_abline(intercept=mean, slope=0, linetype="dashed",colour="black") +
                #ggplot2.theme(**{'legend.title':ggplot2.element_blank()}) +
                #ggplot2.guides(fill=ggplot2.guide_legend(title="Dosage")) +
                #ggplot2.geom_line(ggplot2.aes_string(x="time",y="smooth"),data=r_base,colour="red") +
    figname = "%s/%s_rawtimecourse.pdf" % (outprefix, datagrp)
    print "Writing file:",figname
    grdevices.pdf(figname, width=8, height=10)
    p.plot()
    grdevices.dev_off()
# END plotAllRawTimecourse

def plotAllTimecourse( timecourse, figdir, zscorecut=1.7, fprefix="all" ) :
    timecourse["fishtype_clean"] = ["Control" if x == "ctrl" else "Morpholino" 
                                   for x in timecourse.fishtype]
    timecourse["peak_clean"] = ["Norm" if x == "Norm" else "Peak" 
                                   for x in timecourse.peak]
    timecourse["fishtime"] = timecourse["fishframe"]/3.
    print timecourse.head(10)
    lastframes = timecourse.groupby( ["sample","fishtype_clean","dosage"]
                                   )["fishtime"].min().reset_index()
    lastframes = lastframes[lastframes.dosage != "baseline"]
    print lastframes.head()
    
    r_dataframe = com.convert_to_r_dataframe(timecourse)
    r_lastframes = com.convert_to_r_dataframe(lastframes)
    r_dataframe = fixRLevels( r_dataframe, "peak_clean", ["Max","Peak","Norm"] )
    #r_dataframe = fixRLevels( r_dataframe, "peak", ["Max","Peak","Norm"] )
    #r_lastframes = fixRLevels( r_lastframes, "peak", ["Max","Peak","Norm"] )
    p = (ggplot2.ggplot(r_dataframe) +
                ggplot2.aes_string(x="fishtime",y="zscore" ) +
                ggplot2.geom_point(ggplot2.aes_string(colour="factor(peak_clean)")) +
                ggplot2.geom_line(alpha=.7) +
                #ggplot2.geom_ribbon(ggplot2.aes_string(ymax="zscore", ymin=str(zscorecut))) +
                ggplot2.geom_abline( colour="black", intercept=zscorecut, slope=0, linetype="dashed") + 
                ggplot2.geom_segment(ggplot2.aes_string(
                                    x="fishtime", y="-10", xend="fishtime",yend="50", colour="dosage"), 
                                    linetype="solid", #arrow=grid.arrow(length = unit(0.5, "cm")), 
                                    data=r_lastframes) +
                ggplot2.scale_y_continuous("Normalized Intensity (Z-scores)") +
                ggplot2.scale_x_continuous("Time (seconds)") +
                ggplot2.scale_colour_manual( values=robjects.StrVector(
                    ["purple","#999999", "#56B4E9"]),#"red", 
                    breaks=robjects.StrVector(["5mmptz","Peak", "Norm" ]), #"62mmptz",
                    labels=robjects.StrVector(["start 5mM PTZ","Above threshold", 
                                               "Background"])) + #,"start 62mM PTZ"
                ggplot2.ggtitle("Peak Identification") +
                ggplot2.theme(**{'legend.title':ggplot2.element_blank()}) +
                ggplot2.facet_grid(robjects.Formula('sample ~ fishtype_clean')) +
                ggplot2.theme(**mytheme))
                #ggplot2.scale_colour_manual( values=robjects.StrVector(
                    #["purple","red","#E69F00", "#999999", "#56B4E9"]),
                    #breaks=robjects.StrVector(["5mmptz","62mmptz","Max", "Norm", "Peak"]),
                    #labels=robjects.StrVector(["start 5mM PTZ","start 62mM PTZ","Peak Max", 
                                               #"Background","Above threshold"])) +
                #ggplot2.scale_y_continuous("Normalized Intensity (Z-scores)") +
                #values=c("#999999", "#E69F00", "#56B4E9"), 
                #ggplot2.theme(**mytheme) + \, colour="factor(peak)"
                #ggplot2.theme(**{'axis.text.x': ggplot2.element_text(angle = 45)})  + \
                #ggplot2.geom_boxplot(ggplot2.aes_string(fill="Continent")), notch=True ) + \
    imagefile = "%s/%s_timecourse.pdf" % (figdir, fprefix )
    print "Writing file:",imagefile
    grdevices.pdf(imagefile, width=8, height=10)
    #grdevices.svg(imagefile, width=4, height=5)
    p.plot()
    grdevices.dev_off()
# END plotAllTimecourse

def plotGrpTimecourse( timecourse, figdir, zscorecut=1.7, datagrp="all" ) :
    timecourse["fishtype_clean"] = ["Control" if x == "ctrl" else "Morpholino" 
                                   for x in timecourse.fishtype]
    timecourse["peak_clean"] = ["Norm" if x == "Norm" else "Peak" 
                                   for x in timecourse.peak]
    print timecourse.head(10)
    #lastframes = timecourse.groupby( ["sample","fishtype_clean","dosage"]
                                   #)["time"].min().reset_index()
    #lastframes = lastframes[lastframes.dosage != "baseline"]
    #print lastframes.head()
    
    r_dataframe = com.convert_to_r_dataframe(timecourse)
    #r_lastframes = com.convert_to_r_dataframe(lastframes)
    r_dataframe = fixRLevels( r_dataframe, "peak_clean", ["Max","Peak","Norm"] )
    p = (ggplot2.ggplot(r_dataframe) +
                ggplot2.aes_string(x="time",y="zscore" ) +
                ggplot2.geom_point(ggplot2.aes_string(colour="factor(peak_clean)")) +
                ggplot2.geom_line(alpha=.7) +
                #ggplot2.geom_ribbon(ggplot2.aes_string(ymax="zscore", ymin=str(zscorecut))) +
                ggplot2.geom_abline(intercept=zscorecut, slope=0, linetype="dashed",colour="black") +
                #ggplot2.geom_segment(ggplot2.aes_string(
                                    #x="time", y="-3", xend="time",yend="11", colour="dosage"), 
                                    #linetype="dashed", data=r_lastframes) +
                ggplot2.scale_y_continuous("Normalized Intensity (Z-scores)") +
                ggplot2.scale_x_continuous("Time (seconds)") +
                ggplot2.scale_colour_manual( values=robjects.StrVector(["#56B4E9", "#999999"]),
                                  breaks=robjects.StrVector(["Peak", "Norm"]),
                                  labels=robjects.StrVector(["Above threshold","Background"])) +
                ggplot2.ggtitle("Peak Identification\n("+datagrp+")") +
                ggplot2.theme(**{'legend.title':ggplot2.element_blank()}) +
                ggplot2.facet_grid(robjects.Formula('sample ~ fishtype_clean')) +
                ggplot2.theme(**mytheme))
                #ggplot2.scale_colour_manual( values=robjects.StrVector(["#E69F00", "#56B4E9", "#999999"]),
                                  #breaks=robjects.StrVector(["Max", "Norm", "Peak"]),
                                  #labels=robjects.StrVector(["Peak Max","Above threshold","Background"])) +
                #ggplot2.scale_y_continuous("Normalized Intensity (Z-scores)") +
                #values=c("#999999", "#E69F00", "#56B4E9"), 
                #ggplot2.theme(**mytheme) + \, colour="factor(peak)"
                #ggplot2.theme(**{'axis.text.x': ggplot2.element_text(angle = 45)})  + \
                #ggplot2.geom_boxplot(ggplot2.aes_string(fill="Continent")), notch=True ) + \
    imagefile = "%s/%s_timecourse.pdf" % (figdir, datagrp)
    print "Writing file:",imagefile
    grdevices.pdf(imagefile, width=8, height=10)
    #grdevices.svg(imagefile, width=4, height=5)
    p.plot()
    grdevices.dev_off()
# END plotGrpTimecourse

def plotTimes( prefix, bdosage, ptimedf, factor="Drug", figout="figures" ) :
    r_dataframe = com.convert_to_r_dataframe(ptimedf)
    p = ggplot2.ggplot(r_dataframe) + \
                ggplot2.aes_string(x = "factor(sample)",y="time" ) + \
                ggplot2.geom_boxplot() + \
                ggplot2.scale_y_continuous("Time between Peaks (seconds)") + \
                ggplot2.scale_x_discrete("sample") + \
                ggplot2.ggtitle("Peak Times") + \
                ggplot2.theme(**mytheme) + \
                ggplot2.facet_grid( robjects.Formula(factor+' ~ .') )
                #ggplot2.theme(**{'axis.text.x': ggplot2.element_text(angle = 45)})  + \
                #ggplot2.geom_boxplot(ggplot2.aes_string(fill="Continent")), notch=True ) + \
    figname = "%s/%s_%s_times.png" % (figout,prefix,bdosage.replace(",",""))
    print "Writing file:",figname
    grdevices.png(figname)
    p.plot()
    grdevices.dev_off()
# END plotTimes

def plotIntensities( prefix, bdosage, pintensitydf, factor="Drug",figout="figures" ) :
    r_dataframe = com.convert_to_r_dataframe(pintensitydf)
    p = ggplot2.ggplot(r_dataframe) + \
                ggplot2.aes_string(x = "factor(Rep)",y="pintensity" ) + \
                ggplot2.geom_boxplot() + \
                ggplot2.scale_y_continuous("Boxplot of Normalized Intensity\n(zscores)") + \
                ggplot2.scale_x_discrete("Rep") + \
                ggplot2.ggtitle("Peak Intensities") + \
                ggplot2.theme(**mytheme) + \
                ggplot2.facet_grid( robjects.Formula(factor+' ~ .') )
                #, colour="factor(peak)"
                #ggplot2.theme(**{'axis.text.x': ggplot2.element_text(angle = 45)})  + \
                #ggplot2.geom_boxplot(ggplot2.aes_string(fill="Continent")), notch=True ) + \
    figname = "%s/%s_%s_intensity.png" % (figout,prefix,bdosage.replace(",",""))
    print "Writing file:",figname
    grdevices.png(figname)
    p.plot()
    grdevices.dev_off()
# END plotIntensities

def plotIntensities2( prefix, bdosage, pintensitydf, factor="Drug",figout="figures" ) :
    r_dataframe = com.convert_to_r_dataframe(pintensitydf)
    p = (ggplot2.ggplot(r_dataframe) +
                ggplot2.aes_string(x = "factor("+factor+")",y="pintensity" ) +
                ggplot2.geom_boxplot() +
                ggplot2.scale_y_continuous("Boxplot of Normalized Intensity\n(zscores)") +
                ggplot2.scale_x_discrete("Rep") +
                ggplot2.ggtitle("Peak Intensities") +
                ggplot2.theme(**mytheme))
                #ggplot2.facet_grid( robjects.Formula(factor+' ~ .') )
                #, colour="factor(peak)"
                #ggplot2.theme(**{'axis.text.x': ggplot2.element_text(angle = 45)})  + \
    figname = "%s/%s_%s_intensity2.png" % (figout,prefix,bdosage.replace(",",""))
    print "Writing file:",figname
    grdevices.png(figname)
    p.plot()
    grdevices.dev_off()
# END plotIntensities2

def plotPeaks( peaksdf, figout="figures" ) :
    #prefix, bdosage, 
    r_dataframe = com.convert_to_r_dataframe(peaksdf)
    p = (ggplot2.ggplot(r_dataframe) +
                ggplot2.aes_string(x="minute", y="cumsum", group="factor(label)", colour="factor(dosage)" ) +
                ggplot2.geom_point() + ggplot2.geom_line() +
                ggplot2.scale_y_continuous("Cumulative sum of peaks") +
                ggplot2.scale_x_discrete("Time (minute)") +
                ggplot2.ggtitle("Peaks per minute\n(cumsum)") +
                ggplot2.facet_grid( robjects.Formula('. ~ fishtype') ) +
                ggplot2.theme(**mytheme) )
                #, colour="factor(peak)"
                #ggplot2.theme(**{'axis.text.x': ggplot2.element_text(angle = 45)})  + \
    figname = "%s/sumpeaks.png" % (figout)
    print "Writing file:",figname
    grdevices.png(figname)
    p.plot()
    grdevices.dev_off()
# END plotPeaks

def plotSigIntensities( sigintensities, figdir )  :
    print sigintensities.head(10)

    tlevels = robjects.StrVector([x for x in ["baseline","5mmptz","62mmptz"]
                                  if x in sigintensities.dosage.unique() ])
    print tlevels

    sigintensities["logged"] = [math.log10(x) for x in sigintensities["zscore"]]
    r_dataframe = com.convert_to_r_dataframe(sigintensities)
    s = robjects.FactorVector( r_dataframe.rx2("dosage"), levels=tlevels )
    t = robjects.FloatVector( r_dataframe.rx2("logged") )
    new_r_df = r_dataframe.cbind(s,t)
    new_r_df.colnames = robjects.StrVector( sigintensities.columns.tolist() +
                                             ['dosage_new',"logged_new"] )
    p = (ggplot2.ggplot(new_r_df) +
                ggplot2.aes_string(x="factor(label)", y="time", fill="factor(fishtype)" ) +
                ggplot2.geom_point() + #ggplot2.aes_string( size="logged_new" )) +
                ggplot2.ggtitle("Peak Intensities") +
                ggplot2.theme(**{'axis.text.x': ggplot2.element_text(angle = 45)}) +
                ggplot2.facet_grid(robjects.Formula('. ~ dosage_new')) +
                ggplot2.theme(**mytheme))
                #ggplot2.scale_y_continuous("Normalized Peak Intensity\n(zscores)") +
                #ggplot2.scale_x_discrete("Sample") +
                #, size="zscore" 
    figname = os.path.join(figdir,"all_sig.png")
    print "Writing file:",figname
    grdevices.png(figname)
    p.plot()
    grdevices.dev_off()
# END plotSigIntensities

def calcPairwisePvalues( tdf, groupcol, factcol, onetailed=False, localmax=False ) :
    maxtime = tdf[factcol].max()
    pvals = []
    for grp in tdf[groupcol].unique(): 
        ctrl = tdf[(tdf.fishtype=="ctrl") & (tdf[groupcol] ==grp)][factcol].tolist()
        mo = tdf[(tdf.fishtype=="mo") & (tdf[groupcol] ==grp)][factcol].tolist()
        if localmax : maxtime = max(ctrl+mo)
        ts,pv = ttest_ind(ctrl, mo)
        #p/2 < alpha and t < 0
        if onetailed : print "Todo"
        pvals.append( Series(data=[grp,pv,ts,"ctrl",maxtime], 
                             index=[groupcol,"pvalue","T-statistic","x","y"]))
    pvals = DataFrame(pvals)
    pvals["Pvalue"] = ["pval: %.3g" %x for x in pvals["pvalue"].tolist()]
    print pvals.head()
    return pvals
# END calcPairwisePvalues

def plotAllIntensities( allpeakseries, figdir ) :
    #pintensitydf = []
    #for index, row in allpeakseries.iterrows() :
        #pintensitydf.append( 
            #DataFrame({ "Rep":row["Rep"],
                       #"fishtype":row["fishtype"],
                       #"dosage":row["dosage"],
                       #"pintensity":row["pintensity"]}))
    #pintensitydf = concat(pintensitydf).reset_index(drop=True)
    pintensitydf = allpeakseries[allpeakseries.peak == "Max"]
    pintensitydf.to_csv("intensities.txt",sep="\t")
    pvals = calcPairwisePvalues( pintensitydf, "dosage", "zscore" )

    print pintensitydf[pintensitydf.zscore > 15].head(10)
    #sys.exit(1)
    tlevels = [x for x in ["baseline","5mmptz","62mmptz"]
                                  if x in pintensitydf.dosage.unique() ]
    print pvals
    r_dataframe = com.convert_to_r_dataframe(pintensitydf)
    r_dataframe = fixRLevels( r_dataframe, "dosage", tlevels )
    r_pvals = com.convert_to_r_dataframe(pvals)
    r_pvals = fixRLevels( r_pvals, "dosage", tlevels )

    p = (ggplot2.ggplot(r_dataframe) +
                ggplot2.aes_string(x = "factor(fishtype)",y="zscore" ) +
                ggplot2.geom_boxplot(ggplot2.aes_string(fill="factor(fishtype)")) + #,notch=True
                #ggplot2.geom_jitter(colour="black") + #,notch=True
                ggplot2.geom_text(ggplot2.aes_string(label="Pvalue", x="x", y="y"),
                                  hjust=0, data = r_pvals ) +
                ggplot2.scale_y_continuous("Normalized Peak Intensity\n(zscores)") +
                ggplot2.scale_x_discrete("Condition") +
                ggplot2.ggtitle("Peak Intensities") +
                ggplot2.theme(**{'axis.text.x': ggplot2.element_text(angle = 45),
                                 'legend.position':"NA"}) +
                ggplot2.facet_grid(robjects.Formula('. ~ dosage')) +
                ggplot2.theme(**mytheme))
            #, fill="factor(fishtype)"
    figname = os.path.join(figdir,"all_intensity.png")
    print "Writing file:",figname
    grdevices.png(figname)
    p.plot()
    grdevices.dev_off()
# END plotAllIntensities

def plotAllIntensities2( allpeakseries, figdir, zscorecut=2 ) :
    #pintensitydf = allpeakseries[allpeakseries.peak == "Max"]
    pintensitydf = allpeakseries#[allpeakseries.zscore > zscorecut]
    #pintensitydf = allpeakseries[allpeakseries.dosage == "62mmptz"]
    print "$"*40
    print allpeakseries.head(10)
    print "$"*40
    #pintensitydf.to_csv("intensities.txt",sep="\t")
    pvals = calcPairwisePvalues( pintensitydf, "dosage", "zscore", localmax=True )

    tlevels = [x for x in ["baseline","5mmptz","62mmptz"]
                                  if x in pintensitydf.dosage.unique() ]
    print pvals
    r_dataframe = com.convert_to_r_dataframe(pintensitydf)
    r_dataframe = fixRLevels( r_dataframe, "dosage", tlevels )
    r_pvals = com.convert_to_r_dataframe(pvals)
    r_pvals = fixRLevels( r_pvals, "dosage", tlevels )

                #ggplot2.aes_string(x = "factor(fishtype)",y="zscore" ) +
    p = (ggplot2.ggplot(r_dataframe) +
                ggplot2.aes_string( x="zscore" ) +
                ggplot2.stat_density(ggplot2.aes_string(ymax="..density..", ymin="-..density..",
                                                        fill="factor(fishtype)"),
                                     geom="ribbon", position="identity") + #,notch=True 
                ggplot2.geom_text(ggplot2.aes_string(label="Pvalue", y="0.0", x="y"),
                                 hjust=0, data = r_pvals ) +
                #ggplot2.geom_jitter(colour="black") + #,notch=True
                #ggplot2.scale_x_continuous("Normalized Peak Intensity\n(zscores)") +
                #ggplot2.scale_y_continuous("Density") +
                ggplot2.ggtitle("Peak Intensities") +
                #ggplot2.theme(**{'axis.text.x': ggplot2.element_text(angle = 45),
                                 #'legend.position':"NA"}) +
                ggplot2.facet_grid(robjects.Formula('dosage ~ fishtype'), scales="free") +
                ggplot2.coord_flip() +
                ggplot2.theme(**mytheme))
                #ggplot2.scale_x_discrete("Condition") +
                #ggplot2.geom_density(ggplot2.aes_string(fill="factor(fishtype)")) + #,notch=True
            #, fill="factor(fishtype)"
    figname = os.path.join(figdir,"all_intensityvals.png")
    print "Writing file:",figname
    grdevices.png(figname)
    p.plot()
    grdevices.dev_off()

                #ggplot2.aes_string(x = "factor(fishtype)",y="zscore" ) +
    p = (ggplot2.ggplot(r_dataframe) +
                ggplot2.aes_string( x="zscore" ) +
                ggplot2.geom_density(ggplot2.aes_string(colour="factor(fishtype)")) + #,notch=True
                ggplot2.geom_text(ggplot2.aes_string(label="Pvalue", y="y", x="0.0"), 
                                  hjust=0, data = r_pvals ) +
                #ggplot2.geom_jitter(colour="black") + #,notch=True
                #ggplot2.scale_x_continuous("Normalized Peak Intensity\n(zscores)") +
                #ggplot2.scale_y_continuous("Density") +
                ggplot2.ggtitle("Peak Intensities") +
                #ggplot2.theme(**{'axis.text.x': ggplot2.element_text(angle = 45),
                                 #'legend.position':"NA"}) +
                ggplot2.facet_grid(robjects.Formula('dosage ~ .'), scales="free") +
                #ggplot2.coord_flip() +
                ggplot2.theme(**mytheme))

    figname = os.path.join(figdir,"all_intensityvals_density.png")
    print "Writing file:",figname
    grdevices.png(figname)
    p.plot()
    grdevices.dev_off()

# END plotAllIntensities2

def plotAllPeaktimes( allpeakseries, figdir ) :
    tlevels = [x for x in ["baseline","5mmptz","62mmptz"]
                                  if x in allpeakseries.dosage.unique() ]

    print allpeakseries.head()
    subdata = allpeakseries[allpeakseries.framecount >150][["fishtype","dosage","peakcount"]]
    print subdata.head()
    pvals = calcPairwisePvalues( subdata, "dosage", "peakcount" )
    r_dataframe = com.convert_to_r_dataframe(subdata)
    r_dataframe = fixRLevels( r_dataframe, "dosage", tlevels )

    r_pvals = com.convert_to_r_dataframe(pvals)
    r_pvals = fixRLevels( r_pvals, "dosage", tlevels )

    p = (ggplot2.ggplot(r_dataframe) +
                ggplot2.aes_string(x = "factor(fishtype)",y="peakcount" ) +
                ggplot2.geom_boxplot(ggplot2.aes_string(fill="factor(fishtype)")) +
                #ggplot2.geom_jitter(colour="black") +
                ggplot2.geom_text(ggplot2.aes_string(label="Pvalue", x="x", y="y"), 
                                  hjust=0, data = r_pvals ) +
                ggplot2.scale_x_discrete("Sample") +
                ggplot2.scale_y_continuous("Number of Peaks per Second") +
                ggplot2.ggtitle("Peaks per Second") +
                ggplot2.theme(**{'axis.text.x': ggplot2.element_text(angle = 45),
                                 'legend.position':"NA"}) +
                ggplot2.facet_grid(robjects.Formula('. ~ dosage')) +
                ggplot2.theme(**mytheme))
                #ggplot2.facet_grid( robjects.Formula(factor+' ~ .') )
                #ggplot2.theme(**{'axis.text.x': ggplot2.element_text(angle = 45)})  + \

    figname = os.path.join(figdir,"all_peaktimes.png")
    print "Writing file:",figname
    grdevices.png(figname)
    p.plot()
    grdevices.dev_off()
# END plotAllPeaktimes

def plotPeaksPerSecond( allpeakseries, figdir ) :
    tlevels = [x for x in ["baseline","5mmptz","62mmptz"]
                                  if x in allpeakseries.dosage.unique() ]
    groupingcols = ["sid","fishtype","sample","dosage"]

    ppersec = []
    for metainfo, subseries in allpeakseries.groupby(groupingcols):
        metainfo = Series( metainfo, index=groupingcols )
        framecount = DataFrame({'framecount':subseries.groupby(['peakgrp']).size()}).reset_index()
        pinacles = subseries[subseries.peak == "Max"].sort("fishframe")
        pinacles = merge(pinacles, framecount, on="peakgrp" )
        peakave = (len(pinacles) / ((subseries.fishframe.max())/3.))
        #print metainfo
        #print "Dosage:",metainfo["dosage"],"Peaks per second avg", peakave
        ppersec.append( metainfo.append( Series(data=[peakave], index=["peakave"]) ) )
    ppersec = DataFrame(ppersec)

    pvals = calcPairwisePvalues( ppersec, "dosage", "peakave" )
    r_dataframe = com.convert_to_r_dataframe(ppersec)
    r_dataframe = fixRLevels( r_dataframe, "dosage", tlevels )

    r_pvals = com.convert_to_r_dataframe(pvals)
    r_pvals = fixRLevels( r_pvals, "dosage", tlevels )

    p = (ggplot2.ggplot(r_dataframe) +
                ggplot2.aes_string(x = "factor(fishtype)",y="peakave" ) +
                ggplot2.geom_boxplot(ggplot2.aes_string(fill="factor(fishtype)")) +
                #ggplot2.geom_jitter(colour="black") +
                ggplot2.geom_text(ggplot2.aes_string(label="Pvalue", x="x", y="y"), 
                                  hjust=0, data = r_pvals ) +
                ggplot2.scale_x_discrete("Sample") +
                ggplot2.scale_y_continuous("Average Peaks per Second") +
                ggplot2.ggtitle("Peaks per Second") +
                ggplot2.theme(**{'axis.text.x': ggplot2.element_text(angle = 45),
                                 'legend.position':"NA"}) +
                ggplot2.facet_grid(robjects.Formula('. ~ dosage')) +
                ggplot2.theme(**mytheme))
                #ggplot2.facet_grid( robjects.Formula(factor+' ~ .') )
                #ggplot2.theme(**{'axis.text.x': ggplot2.element_text(angle = 45)})  + \

    figname = os.path.join(figdir,"all_peakpermin.png")
    print "Writing file:",figname
    grdevices.png(figname)
    p.plot()
    grdevices.dev_off()
# END plotPeaksPerSecond

def plotDuration( allpeakseries, figdir ) :
    tlevels = [x for x in ["baseline","5mmptz","62mmptz"]
                                  if x in allpeakseries.dosage.unique() ]
    print tlevels

    framesdf = allpeakseries[ (allpeakseries.peak=="Max") & (allpeakseries.peakgrp != 0) ]

    framesdf["peakdecay"] = (framesdf["peakend"] - framesdf["fishframe"]) / 3.

    pvals = calcPairwisePvalues( framesdf, "dosage", "peakduration" )
    r_dataframe = com.convert_to_r_dataframe(framesdf)
    r_dataframe = fixRLevels( r_dataframe, "dosage", tlevels )
    r_pvals = com.convert_to_r_dataframe(pvals)
    r_pvals = fixRLevels( r_pvals, "dosage", tlevels )
    p = (ggplot2.ggplot(r_dataframe) +
                ggplot2.aes_string(x = "factor(fishtype)",y="peakduration" ) +
                ggplot2.geom_boxplot(ggplot2.aes_string(fill="factor(fishtype)")) + #, notch=True
                ggplot2.geom_text(ggplot2.aes_string(label="Pvalue", x="x", y="y"), 
                                  hjust=0, data = r_pvals ) +
                ggplot2.scale_x_discrete("Condition") +
                ggplot2.scale_y_continuous("Time (seconds)") +
                ggplot2.ggtitle("Duration of Events") +
                ggplot2.theme(**{'axis.text.x': ggplot2.element_text(angle = 45),
                                 'legend.position':"NA"}) +
                ggplot2.facet_grid(robjects.Formula('. ~ dosage')) +
                ggplot2.theme(**mytheme))
                #ggplot2.guides(fill=guide_legend(title="Condition") +
                #ggplot2.facet_grid( robjects.Formula(factor+' ~ .') )
                #ggplot2.theme(**{'axis.text.x': ggplot2.element_text(angle = 45)})  + \

    figname = os.path.join(figdir,"all_frames.png")
    print "Writing file:",figname
    grdevices.png(figname)
    p.plot()
    grdevices.dev_off()

    pvals = calcPairwisePvalues( framesdf, "dosage", "peakdecay" )
    r_dataframe = com.convert_to_r_dataframe(framesdf)
    r_dataframe = fixRLevels( r_dataframe, "dosage", tlevels )
    r_pvals = com.convert_to_r_dataframe(pvals)
    r_pvals = fixRLevels( r_pvals, "dosage", tlevels )
    p = (ggplot2.ggplot(r_dataframe) +
                ggplot2.aes_string(x = "factor(fishtype)",y="peakdecay" ) +
                ggplot2.geom_boxplot(ggplot2.aes_string(fill="factor(fishtype)")) + #, notch=True
                ggplot2.geom_text(ggplot2.aes_string(label="Pvalue", x="x", y="y"), 
                                  hjust=0, data = r_pvals ) +
                ggplot2.scale_x_discrete("Condition") +
                ggplot2.scale_y_continuous("Time (seconds)") +
                ggplot2.ggtitle("Time to Decay for Firing Events") +
                ggplot2.theme(**{'axis.text.x': ggplot2.element_text(angle = 45),
                                 'legend.position':"NA"}) +
                ggplot2.facet_grid(robjects.Formula('. ~ dosage')) +
                ggplot2.theme(**mytheme))
                #ggplot2.guides(fill=guide_legend(title="Condition") +
                #ggplot2.facet_grid( robjects.Formula(factor+' ~ .') )
                #ggplot2.theme(**{'axis.text.x': ggplot2.element_text(angle = 45)})  + \

    figname = os.path.join(figdir,"firingdecay.png")
    print "Writing file:",figname
    grdevices.png(figname)
    p.plot()
    grdevices.dev_off()

    framesdf["peakonset"] = (framesdf["fishframe"] - framesdf["peakstart"]) / 3.
    print "Max peakonset"
    print framesdf.ix[framesdf["peakonset"].idxmax()]
    pvals = calcPairwisePvalues( framesdf, "dosage", "peakonset" )
    r_dataframe = com.convert_to_r_dataframe(framesdf)
    r_dataframe = fixRLevels( r_dataframe, "dosage", tlevels )
    r_pvals = com.convert_to_r_dataframe(pvals)
    r_pvals = fixRLevels( r_pvals, "dosage", tlevels )
    p = (ggplot2.ggplot(r_dataframe) +
                ggplot2.aes_string(x = "factor(fishtype)",y="peakonset" ) +
                ggplot2.geom_boxplot(ggplot2.aes_string(fill="factor(fishtype)")) + #, notch=True
                ggplot2.geom_text(ggplot2.aes_string(label="Pvalue", x="x", y="y"), 
                                  hjust=0, data = r_pvals ) +
                ggplot2.scale_x_discrete("Condition") +
                ggplot2.scale_y_continuous("Time (seconds)") +
                ggplot2.ggtitle("Time to Intensity Maximum for Firing Events") +
                ggplot2.theme(**{'axis.text.x': ggplot2.element_text(angle = 45),
                                 'legend.position':"NA"}) +
                ggplot2.facet_grid(robjects.Formula('. ~ dosage')) +
                ggplot2.theme(**mytheme))
                #ggplot2.guides(fill=guide_legend(title="Condition") +
                #ggplot2.facet_grid( robjects.Formula(factor+' ~ .') )
                #ggplot2.theme(**{'axis.text.x': ggplot2.element_text(angle = 45)})  + \

    figname = os.path.join(figdir,"firingonset.png")
    print "Writing file:",figname
    grdevices.png(figname)
    p.plot()
    grdevices.dev_off()
# END PlotDuration
    #framesdf = []
#for index, row in allpeakseries.iterrows() :
    #framesdf.append( 
        #DataFrame({ "Rep":row["Rep"],
                   #"fishtype":row["fishtype"],
                   #"dosage":row["dosage"],
                   #"frames":row["frames"]}))
#framesdf = concat(framesdf).reset_index(drop=True)
#framesdf["time"] = framesdf.frames.map(float) / 3
    #maxtime = framesdf.time.max()
    #pvals = []
    #for dosage in tlevels: 
        #ctrl = framesdf[(framesdf.fishtype=="ctrl") & (framesdf.dosage ==dosage)].time.tolist()
        #mo = framesdf[(framesdf.fishtype=="mo") & (framesdf.dosage ==dosage)].time.tolist()
        #ts,pv = ttest_ind(ctrl, mo)
        #pvals.append( Series(data=[dosage,pv,ts,"ctrl",maxtime], 
                             #index=["dosage_new","pvalue","T-statistic","x","y"]))
    #pvals = DataFrame(pvals)
    #pvals["Pvalue"] = ["pval: %.3f" %x for x in pvals["pvalue"].tolist()]
    #print pvals.head()


#------------------------------------------------------------------------------#
#                                 Process Data
#------------------------------------------------------------------------------#
def compareAttributes( allpeakseries, figdir ) :
    # Compare between Drugs
    for typedosage, data in allpeakseries.groupby(["fishtype","dosage"]) :
        fishtype, bdosage = typedosage
        ptimedf = []
        pintensitydf = []
        
        for index, row in data.iterrows() :
            ptimedf.append( 
                DataFrame({"Drug":row["Drug"],"Rep":row["Rep"],"ptime":row["ptime"]})
            )
            pintensitydf.append( 
                DataFrame({"Drug":row["Drug"],"Rep":row["Rep"],"pintensity":row["pintensity"]})
            )
        ptimedf = concat(ptimedf).reset_index(drop=True)
        pintensitydf = concat(pintensitydf).reset_index(drop=True)
        plotTimes( fishtype, bdosage, ptimedf, figout=figdir )
        plotIntensities( fishtype, bdosage, pintensitydf, figout=figdir )

    # Compare between fishtype
    for typedosage, data in allpeakseries.groupby(["Drug","dosage"]) :
        drug, bdosage = typedosage
        ptimedf = []
        pintensitydf = []
        
        for index, row in data.iterrows() :
            ptimedf.append( 
                DataFrame({"fishtype":row["fishtype"],"Rep":row["Rep"],"ptime":row["ptime"]})
            )
            pintensitydf.append( 
                DataFrame({"fishtype":row["fishtype"],"Rep":row["Rep"],"pintensity":row["pintensity"]})
            )
        ptimedf = concat(ptimedf).reset_index(drop=True)
        pintensitydf = concat(pintensitydf).reset_index(drop=True)
        plotTimes( drug, bdosage, ptimedf, factor="fishtype", figout=figdir )
        plotIntensities( drug, bdosage, pintensitydf, factor="fishtype", figout=figdir )
        plotIntensities2( drug, bdosage, pintensitydf, factor="fishtype", figout=figdir )
# END compareAttributes

def normalizeTimeSeries( timecourse, normalization ):
    #print timecourse.head(10)
    timecourse["time"] = timecourse.frame.map(int) / 3.
    timecourse["minutes"] = (timecourse.time / 60.).map(int)
    timecourse["RawZscore"] = stats.zscore(timecourse.intensity)
    if normalization == "smooth50" :
        timecourse["smooth"] = smooth(numpy.array(timecourse.intensity),
                                      window_len=50)
        timecourse["norm"] = timecourse.intensity - timecourse.smooth
    elif normalization == "smooth100" :
        timecourse["smooth"] = smooth(numpy.array(timecourse.intensity),
                                      window_len=100)
        timecourse["norm"] = timecourse.intensity - timecourse.smooth
    elif normalization == "smooth200" :
        timecourse["smooth"] = smooth(numpy.array(timecourse.intensity),
                                      window_len=200)
        timecourse["norm"] = timecourse.intensity - timecourse.smooth
    elif normalization == "smooth300" :
        timecourse["smooth"] = smooth(numpy.array(timecourse.intensity),
                                      window_len=300)
        timecourse["norm"] = timecourse.intensity - timecourse.smooth
    elif normalization == "rawpeaks" :
        timecourse["norm"] = timecourse.intensity
    else : 
        print "Error: normalization method unknown:", normalization
        sys.exit(1)

    timecourse["zscore"] = stats.zscore(timecourse.norm)
    return timecourse
# END normalizeTimeSeries

# Can also use cut of 2 for 97.5%
def identifyPeaks( timeseries, zscorecut=1.7 ) : 
    #timecourse["peak"] = ["Peak" if ispeak else "Norm" 
                          #for ispeak in (timecourse.zscore > zscorecut).tolist() ]
    timeseries["peak"] = ["Peak" if ispeak else "Norm" 
                          for ispeak in (timeseries.zscore > zscorecut).tolist() ]
    timeseries["peakgrp"] = 0
    grpnum = 1
    prev = False
    for index, row in timeseries.iterrows() :
        if row["peak"] == "Peak" :  
            timeseries.loc[index,"peakgrp"] = grpnum
            prev=True
        elif prev : grpnum += 1; prev=False 

    tokeep = (timeseries[timeseries.peakgrp > 0]
              .groupby(['peakgrp']).apply(lambda x: x['zscore'].idxmax()))
    timeseries.loc[tokeep,"peak"] = "Max"
    framecount = DataFrame({'framecount':timeseries.groupby(['peakgrp']).size(),
                            'peakstart':timeseries.groupby(['peakgrp'])["fishframe"].min(),
                            'peakend':timeseries.groupby(['peakgrp'])["fishframe"].max()}
                            #'maxintensity':timeseries.groupby(['peakgrp'])["intensity"].max(),
                            #'maxtime':timeseries.groupby(['peakgrp'])["intensity"].max(),
                          ).reset_index()
    print framecount.head()
    timeseries = merge(timeseries, framecount, on="peakgrp" )
    timeseries["peakduration"] = timeseries["framecount"] / 3.
    print timeseries.head()
    return timeseries
# END identifyPeaks

def calculatePeaktimes( timecourse ):
    grpcols = ["sid","fishtype","sample","dosage","minutes"]
    tmp = timecourse.groupby(grpcols)
    ppersec = tmp.size().reset_index()
    ppersec.columns = grpcols+["framecount"]
    ppersec["peakcount"] = [sum(data.peak == "Max") for grp,data in tmp]
    return ppersec
# ENd calculatePeaktimes

def readAllData( allfiles, normmethod, figdir, zscorecut=2) :
    dosageorder = {"baseline":1,"5mmptz":2,"62mmptz":3}
    allfiles["dosageorder"] = [dosageorder[x] for x in allfiles.dosage]
    allfiles.sort(["dosageorder","sample"],inplace=True)
    alltimeseries = []
    for idx, fileinfo in allfiles.groupby(["sample","fishtype"]) :
        samp, fishtype = idx
        print "Sample",samp
        timecourse = []
        for filename,dosage in fileinfo[["filename","dosage"]].values :
            print "Reading file:",filename
            tseries = read_csv(filename,header=None,names=["intensity","empty"],sep="\t") 
            tseries["dosage"] = dosage
            tseries["frame"] = tseries.index
            if dosage == "baseline" :
                normtseries = normalizeTimeSeries( tseries, "smooth100" )
            else :
                normtseries = normalizeTimeSeries( tseries, normmethod )
            timecourse.append(normtseries)
        timecourse = concat( timecourse ).reset_index(drop=True)
        timecourse["fishframe"] = timecourse.index
        timecourse["sample"] = samp
        timecourse["fishtype"] = fishtype
        sid = "%s-%s" % (fishtype,samp)
        timecourse["sid"] = sid
        alltimeseries.append( timecourse )
        #normtimeseries = normalizeTimeSeries( timecourse, normmethod )
        #plotRawTimecourse( normtimeseries, figdir, sid ) 
        #alltimeseries.append( normtimeseries )
    alltimeseries = concat( alltimeseries ).reset_index(drop=True)
    #alltimeseries.head(10)

    # Create two distributions
    #allzscores = []
    #for grp, data in alltimeseries.groupby( "fishtype" ) :
        #tmp = stats.zscore(data.norm)
        #allzscores.append(DataFrame( {'zscore':tmp}, index=data.index ))
    #allzscores = concat(allzscores)
    #alltimeseries.update( allzscores )

    # Create dosage distributions
    #allzscores = []
    #for grp, data in alltimeseries.groupby( "dosage" ) :
        #tmp = stats.zscore(data.norm)
        #allzscores.append(DataFrame( {'zscore':tmp}, index=data.index ))
    #allzscores = concat(allzscores)
    #alltimeseries.update( allzscores )

    # Control baseline normalization
    ctrlmean = numpy.mean(alltimeseries[(alltimeseries.dosage == "baseline") & 
                       (alltimeseries.fishtype =="ctrl")].norm)
    ctrlstd = numpy.std(alltimeseries[(alltimeseries.dosage == "baseline") & 
                       (alltimeseries.fishtype =="ctrl")].norm)

    print "Control values:", ctrlmean, ctrlstd

    alltimeseries["zscore"] = [ (x - ctrlmean)/ ctrlstd for x in alltimeseries.norm ]

    # Use a single distribution
    #alltimeseries["zscore"] = stats.zscore(alltimeseries.norm)

    allpeakseries = []
    peakspersec = []
    for samp, timecourse in alltimeseries.groupby("sid") :
        #print "Found",fishtype,rep,dosage
        timecoursewpeaks = identifyPeaks( timecourse, zscorecut )
        allpeakseries.append(timecoursewpeaks)
        #plotTimecourse( timecourse, figdir, samp ) 
        ppersec = calculatePeaktimes( timecourse )
        peakspersec.append(ppersec)

    allpeakseries = concat( allpeakseries ).reset_index()
    peakspersec = concat( peakspersec ).reset_index()
    return allpeakseries, peakspersec
# END readAllData

################################################################################
# MAIN
################################################################################
if __name__ == "__main__" :
    #movfiles = filesInDir( '.', ".mov" )
    #ctrlfiles = filesInDirRecursive( 'rawdata/4.16.14', ".txt" )
    #ctrlfiles = filesInDirRecursive( 'rawdata/5.20.14', ".txt" )
    #ctrlfiles = filesInDirRecursive( 'rawdata/merge', ".txt" )
    allfiles = filesInDirRecursive( 'rawdata/merge1', ".txt" )

    allfiles = allfiles[allfiles.dosage != "62mmptz"]
    print allfiles.head()

    #normmethod = "smooth100"
    normmethod = "smooth200"
    #normmethod = "smooth300"
    zscorecut = 4
    #normmethod = "smooth50"
    #normmethod = "rawpeaks"
    figdir = os.path.abspath(normmethod)
    print figdir
    assert os.path.exists(figdir)

    alltimeseries, peakspersec = readAllData( allfiles, normmethod, figdir, zscorecut=zscorecut )

    # plot rawvalues
    #plotAllRawTimecourse( alltimeseries, figdir ) 
    #for dose, timeseries in alltimeseries.groupby("dosage") :
        #plotAllRawTimecourse( timeseries, figdir, dose ) 

    #Ok, It was 3, 7 and 8 (ctrl) and 6, 7, 8 (MO) for the baseline
    #and 6, 7, 8 for both control and mo for the 5mM

    plotAllTimecourse( alltimeseries, figdir, zscorecut ) 

    subtimeseries = alltimeseries[ alltimeseries.sample.isin(["2","7","8"]) ]

    print subtimeseries.head(10)
    plotAllTimecourse( subtimeseries, figdir, zscorecut, fprefix="subset" ) 

    for dose, peakseries in alltimeseries.groupby("dosage") :
        plotGrpTimecourse( peakseries, figdir, zscorecut, dose ) 
        subpeakseries =  peakseries[ peakseries.sample.isin(["2","7","8"]) ]
        plotGrpTimecourse( subpeakseries, figdir, zscorecut, "subset"+dose ) 
   
    print peakspersec.head(10)
    print alltimeseries.head(10)

    #compareAttributes( alltimeseries, figdir )

    plotAllIntensities( alltimeseries, figdir )
    plotAllIntensities2( alltimeseries, figdir, zscorecut=zscorecut )
    #plotAllPeaktimes( peakspersec, figdir )
    plotPeaksPerSecond( alltimeseries, figdir )
    plotDuration( alltimeseries, figdir )

    #plotPeaks( peakspersec, figout=figdir )
    #plotSigIntensities( sigintensities, figdir )
    #filedata["newid"] = filedata["fishtype"]+" - "+filedata["dosage"]
    #firsttype = "baseline" 
    #secondtype = "62mmptz"
    #secondtype = "5mmptz"
    #calcPvals( subdata, pintensitydf, firsttype, secondtype )


# END Main
################################################################################



