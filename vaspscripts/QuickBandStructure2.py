#!/usr/bin/env python

# region: Imports

import csv;
import math;
import os;

import numpy as np;
import matplotlib as mpl;
import matplotlib.pyplot as plt;
import matplotlib.ticker as tck;

from scipy.interpolate import splrep, splev

#EndRegion


#Region:Parameters

InputFiles = [
	("./", "OUTCAR"),
	]; # This code allows multiple files to be read in. 
           # Each file read in will have a different coloour band structure.

OutputFileName = "bandstructure.png";

BandPaths = [
	[
                ## FCC cubic lattice
	#	((0.0, 0.0, 0.0), r"$\Gamma$"),
	#	((0.5, 0.0, 0.5), "X"),
	#	((0.5, 0.25, 0.75), "W"),
	#	((0.375, 0.375, 0.75), "K"),
	#	((0.0, 0.0, 0.0), r"$\Gamma$"),
	#	((0.5, 0.5, 0.5), "L"),
	#	((0.625, 0.25, 0.625), "U"),
	#	((0.5, 0.25, 0.75), "W"),
	#	((0.5, 0.5, 0.5), "L"),
	#	((0.625, 0.625, 0.75), "K")

                ## Tetragonal lattice
	#	((0.0, 0.0, 0.0), r"$\Gamma$"),
	#	((0.0, 0.5, 0.00), "X"),
	#	((0.5, 0.5, 0.0), "M"),
             #   ((0.0, 0.0, 0.0), r"$\Gamma$"),
            #    ((0.0,0.0,0.5), "Z"),
           #     ((0.0,0.5,0.5), "R"),
          #      ((0.5,0.5,0.5), "A"),
         #       ((0.0,0.0,0.5),"Z")
                
                ## Simple cubic
                ((0.0, 0.0, 0.0), r"$\Gamma$"),
                ((0.0, 0.5, 0.0), "X"),
                ((0.5, 0.5, 0.0), "M"),
                ((0.0, 0.0, 0.0),r"$\Gamma$"),
                ((0.5, 0.5, 0.5), "R"),
                ((0.0, 0.5, 0.0), "X")
		],
#	 [
#		((0.625, 0.250, 0.625), "U"),
#		((0.500, 0.000, 0.500), "X")
#		]
	];

ReciprocalLatticeVectors = [
        (0.99889751, -0.00001237,  0.00629258),
	(-0.,          1.00148158,  0.00002446),
	(0.,          0.,          0.99781961),
	];

UseInterpolation = False; #This uses a scipy routine to interpolate between k-points

Colors = [(0, 0, 255), (128, 0, 255), (255, 0, 255), (255, 0, 128), (255, 0, 0), (255, 128, 0)] 
# If only one file being plotted, amend first entry to change colour of plot

#EndRegion


#Region: Functions

def _ReadBandEnergies(filePath):
	bandEnergies = [];
	
	inputReader = open(filePath, 'rU');
	inputReaderCSV = csv.reader(inputReader);
	
	next(inputReaderCSV);
	next(inputReaderCSV);
	
	for row in inputReaderCSV:
		bandEnergies.append([float(row[i]) for i in range(1, len(row), 2)]);
	
	inputReader.close();
	
	return bandEnergies;

def _ReadKPoints(filePath):
	kPointCoordinates = [];
	
	inputReader = open(filePath, 'rU');
	inputReaderCSV = csv.reader(inputReader);
	
	next(inputReaderCSV);
	
	for row in inputReaderCSV:
		kPointCoordinates.append((float(row[1]), float(row[2]), float(row[3])));
	
	inputReader.close();
	
	return kPointCoordinates;

#EndRegion


#Region: Main

#COMMENT: Unpack the reciprocal lattice vectors - used in the next two subsections.

(r1X, r1Y, r1Z), (r2X, r2Y, r2Z), (r3X, r3Y, r3Z) = ReciprocalLatticeVectors;

#COMMENT: Prepare the x-axis ticks positions and labels.

specialPointDistances = [0.0];
specialPointLabels = [];

for bandPath in BandPaths:
	for i in range(1, len(bandPath)):
		(kx1, ky1, kz1), label1 = bandPath[i - 1];
		(kx2, ky2, kz2), label2 = bandPath[i];
		
		kDistanceX = (kx2 * r1X + ky2 * r2X + kz2 * r3X) - (kx1 * r1X + ky1 * r2X + kz1 * r3X);
		kDistanceY = (kx2 * r1Y + ky2 * r2Y + kz2 * r3Y) - (kx1 * r1Y + ky1 * r2Y + kz1 * r3Y);
		kDistanceZ = (kx2 * r1Z + ky2 * r2Z + kz2 * r3Z) - (kx1 * r1Z + ky1 * r2Z + kz1 * r3Z);
		
		pathLength = math.sqrt(kDistanceX ** 2 + kDistanceY ** 2 + kDistanceZ ** 2);
		
		specialPointDistances.append(specialPointDistances[-1] + pathLength);
		
		if i == 1:
			if len(specialPointLabels) == 0:
				specialPointLabels.append(label1);
				specialPointLabels.append(label2);
			else:
				specialPointLabels[-1] = specialPointLabels[-1] + "|" + label1;
				specialPointLabels.append(label2);
		else:
			specialPointLabels.append(label2);

#COMMENT: Prepare sets of data for plotting.

plotSets = [];

for i in range(0, len(InputFiles)):
	basePath, prefix = InputFiles[i];
	
	kPoints = _ReadKPoints(os.path.join(basePath, "{0} - K-Points.csv".format(prefix)));
	bandEnergies = _ReadBandEnergies(os.path.join(basePath, "{0} - Band Energies.csv".format(prefix)));
	
	kPointDistances = [0.0];
	
	for j in range(1, len(kPoints)):
		kx1, ky1, kz1 = kPoints[j - 1];
		kx2, ky2, kz2 = kPoints[j];
		
		kDistanceX = (kx2 * r1X + ky2 * r2X + kz2 * r3X) - (kx1 * r1X + ky1 * r2X + kz1 * r3X);
		kDistanceY = (kx2 * r1Y + ky2 * r2Y + kz2 * r3Y) - (kx1 * r1Y + ky1 * r2Y + kz1 * r3Y);
		kDistanceZ = (kx2 * r1Z + ky2 * r2Z + kz2 * r3Z) - (kx1 * r1Z + ky1 * r2Z + kz1 * r3Z);
		
		pathLength = math.sqrt(kDistanceX ** 2 + kDistanceY ** 2 + kDistanceZ ** 2);
		
		kPointDistances.append(kPointDistances[-1] + pathLength);
	
	interpolatedKPointDistances = None;
	interpolatedBandEnergies = None;
	
	if UseInterpolation:
		interpolatedKPointDistances = np.linspace(specialPointDistances[0], specialPointDistances[-1], len(kPointDistances) * 10);
		interpolatedBandEnergies = [];
		
		for band in bandEnergies:
			tck = splrep(kPointDistances, band);
			interpolatedBandEnergies.append(splev(interpolatedKPointDistances, tck));
	
	plotSets.append((kPointDistances, bandEnergies, interpolatedKPointDistances, interpolatedBandEnergies));

#COMMENT: Plot and save the data.

mpl.rc('font', **{ 'family' : 'serif', 'size' : 14, 'serif' : 'Arial' });

plt.figure(figsize =  (12/ 2.54, 12 / 2.54)); #2.54 cm in an inch - matplotlib uses inches

for i in range(0, len(plotSets)):
                print len(plotSets)
		kPointDistances, bandEnergies, interpolatedKPointDistances, interpolatedBandEnergies = plotSets[i];
		r, g, b = Colors[i];
                if interpolatedKPointDistances != None:
			for j in range(0, len(bandEnergies)):
		#		plt.scatter(kPointDistances, bandEnergies[j], 3, facecolor = 'none', edgecolor = (r / 255.0, g / 255.0, b / 255.0), linewidth = 0.5);
				plt.plot(interpolatedKPointDistances, interpolatedBandEnergies[j], color = (r / 255.0, g / 255.0, b / 255.0), linewidth = 1.0);
		else:
			for band in bandEnergies:
				plt.plot(kPointDistances, band, color = (r / 255.0, g / 255.0, b / 255.0), linewidth = 1.0); #matplotlib likes rgb colors as a decimal

plt.xlim((specialPointDistances[0], specialPointDistances[-1]));
plt.ylim((-5, 5)); # change to plot over larger energy interval

plt.xticks(specialPointDistances, specialPointLabels);
plt.yticks([-5.0, -2.5, 0, 2.5, 5]);

plt.xlabel("Wavevector") # 
plt.ylabel ("Energy (eV)") # you can use latex (r"$\mathit{E - E_F}$ / eV")
# plt.title ("CZTS hybrid", fontweight = 'bold')
axes = plt.gca();

axes.get_xaxis().set_tick_params(**{ 'direction' : 'outward', 'top' : 'off'});
axes.get_yaxis().set_tick_params(**{ 'direction' : 'outward', 'right' : 'off'});

axes.get_xaxis().set_tick_params(width = 0.5);
axes.get_yaxis().set_tick_params(width = 0.5);

axes.get_xaxis().grid(color = (211 / 255.0, 211 / 255.0, 211 / 255.0), linestyle = '-', linewidth = 0.5);
# axes.get_yaxis().set_major_formatter(tck.FuncFormatter(lambda value, pos : "{0:.1f}".format(value)));

axes.set_axisbelow(True);

for spine in axes.spines.values():
	spine.set_linewidth(0.5);

plt.tight_layout();

plt.savefig(OutputFileName, format = 'png', dpi = 300);

plt.close();

#EndRegion
