# Segment nuclei in C1, measure intensity statistics in all other channels.
#	- by Richard Butler, Gurdon Institute Imaging Facility

import math as maths

from ij import IJ, ImagePlus, ImageStack
from ij.gui import Roi, TextRoi, ShapeRoi, Overlay
from ij.process import StackStatistics, Blitter, ImageProcessor, ByteProcessor, AutoThresholder, FloodFiller
from ij.plugin import Duplicator
from ij.plugin.filter import ThresholdToSelection, MaximumFinder, EDM
from ij.measure import ResultsTable

from java.awt import Color, Font



minA = 12.0 # µm²
maxA = 150.0

sigma = 15.0	# µm		- 63X images
k = 5.0

tolerance = 0.7

channelNames = ("", "DAPI", "5mc", "Tomato", "OCT4", "", "")


def DoG(ip, sigma, k):
	proc = ip.duplicate()
	sub = ip.duplicate()
	proc.blurGaussian(sigma)
	sub.blurGaussian(k*sigma)
	proc.copyBits(sub, 0,0, Blitter.SUBTRACT)
	return proc


def fillHoles(mask):
	width = mask.getWidth()
	height = mask.getHeight()
	ff = FloodFiller(mask)
	mask.setColor(127)
	foreground = 127
	background = 0
	for y in range(height):
	    if mask.getPixel(0,y)==background:
	    	ff.fill(0, y)
	    if mask.getPixel(width-1,y)==background:
	    	ff.fill(width-1, y)
	for x in range(width):
	    if mask.getPixel(x,0)==background:
	    	ff.fill(x, 0)
	    if mask.getPixel(x,height-1)==background:
	    	ff.fill(x, height-1)
	n = width*height
	for i in range(n):
		if mask.get(i)==127:
		    mask.set(i, 0)
		else:
		    mask.set(i, 255)


def watershed(ip, tol):
	floatEdm = EDM().makeFloatEDM(ip, 0, False)
	maxIp = MaximumFinder().findMaxima(floatEdm, tol, ImageProcessor.NO_THRESHOLD, MaximumFinder.SEGMENTED, False, True)
	if (maxIp != None):
		ip.copyBits(maxIp, 0, 0, Blitter.AND)


def getMask(ip, method):
	stats = ip.getStatistics()
	hist = ip.getHistogram(256)
	thresh = AutoThresholder().getThreshold( AutoThresholder.Method.Huang, hist )
	thresh = (thresh/float(255)) * (stats.max-stats.min) + stats.min
	ip.threshold(int(thresh))
	mask = ip.convertToByte(False)
	return mask


def getRois(mask):
	rois = []
	mask.setThreshold(255, 255, ImageProcessor.NO_LUT_UPDATE)
	composite = ThresholdToSelection().convert(mask)
	rois = ShapeRoi(composite).getRois()
	return rois


imp = IJ.getImage()
W = imp.getWidth()
H = imp.getHeight()
Z = imp.getNSlices()
title = imp.getTitle()
cal = imp.getCalibration()
ol = Overlay()

proc = DoG(imp.getStack().getProcessor(1), sigma, k)
mask = getMask(proc, AutoThresholder.Method.Otsu)
fillHoles(mask)
watershed(mask, tolerance)
rois = getRois(mask)
rt = ResultsTable()
for roi in rois:
	roiA = roi.getStatistics().area * cal.pixelWidth * cal.pixelHeight
	if roiA >= minA and roiA <= maxA:
		row = rt.getCounter()
		bounds = roi.getBounds()
		rt.setValue("X", row, (bounds.x+(bounds.width/2.0))*cal.pixelWidth )
		rt.setValue("Y", row, (bounds.y+(bounds.height/2.0))*cal.pixelHeight )
		rt.setValue("Area", row, roiA )
		for c in range(2,imp.getNChannels()+1):
			measure = imp.getStack().getProcessor(c)
			measure.setRoi(roi)
			roiStats = measure.getStatistics()
			rt.setValue(channelNames[c]+" (C"+str(c)+") Mean", row, roiStats.mean )
			rt.setValue(channelNames[c]+" (C"+str(c)+") StdDev", row, roiStats.stdDev )
			rt.setValue(channelNames[c]+" (C"+str(c)+") Min", row, roiStats.min )
			rt.setValue(channelNames[c]+" (C"+str(c)+") Max", row, roiStats.max )
		roi.setStrokeColor(Color.MAGENTA)
		ol.add(roi)
imp.setOverlay(ol)
rt.show(title+"_Nuclei")
