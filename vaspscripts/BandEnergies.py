#!/usr/bin/env python3

#BandEnergies.py


#Region: Parameters

#COMMENT: Path to an OUTCAR file to read; defaults to "OUTCAR" if set to None.

InputFile = "OUTCAR";

#COMMENT: Unless set to an empty string, all output files will have Prefix prepended to them, e.g. "{Prefix} - Bands.csv", "{Prefix} - KPoints.csv"; defaults to the name of the input file (minus any file extension) if set to None.

Prefix = None;

#COMMENT: This tolerance sets the (non-zero) occupancy below which a band should be regarded as empty; this is intended to provide a means to avoid gaps being calculated incorrectly due to smearing leading to very small occupation of the lowest unoccupied band.

PartialOccupancyThreshold = 0.1;

#EndRegion


#Region: Imports

import csv;
import math;
import os;
import re;
from IPython import embed
#EndRegions


#Region: ClassDefinitions

class BandAnalyser:
	#COMMENT: Constructor
	
	def __init__(self):
		self._isSpinPolarisedCalculation = None;
		self._isSpinOrbitCouplingCalculation = None;
		
		self._numberOfElectrons = None;
		
		self._fermiEnergy = None;
		self._kPointCoordinates = None;
		self._kPointWeights = None;
		self._kPointBandData = None;
		
		self._directBandgap = None;
		self._indirectBandgap = None;
		self._kPointValenceConductionBands = None;
	
	#COMMENT: Properties
	
	@property
	def SpinPolarisedCalculation(self):
		return self._isSpinPolarisedCalculation;
	
	@property
	def SpinOrbitCouplingCalculation(self):
		return self._isSpinOrbitCouplingCalculation;
	
	@property
	def NumberOfElectrons(self):
		return self._numberOfElectrons;
	
	@property
	def FermiEnergy(self):
		return self._fermiEnergy;
	
	@property
	def KPointCoordinates(self):
		return self._kPointCoordinates;
	
	@property
	def KPointWeights(self):
		return self._kPointWeights;
	
	@property
	def KPointBandData(self):
		return self._kPointBandData;
	
	@property
	def DirectBandgap(self):
		return self._directBandgap;
	
	@property
	def IndirectBandgap(self):
		return self._indirectBandgap;
	
	@property
	def KPointValenceConductionBands(self):
		return self._kPointValenceConductionBands;
	
	#COMMENT: PrivateMethods
	
	def _GetZeroWeightKPointIndices(self):
		kPointWeights = self._kPointWeights;
		
		zeroWeightedKPointIndices = [];
		
		for i in range(0, len(kPointWeights)):
			if kPointWeights[i] == 0.0:
				zeroWeightedKPointIndices.append(i);
		
		return zeroWeightedKPointIndices;
	
	#COMMENT: PublicMethods
	
	def ReadOUTCARFile(self, filePath = "OUTCAR", performBandgapAnalysis = True, bandgapAnalysisOnlyAnalyseZeroWeightKPoints = True, bandgapAnalysisPartialOccupancyThreshold = 0.01, quiet = False):
		if not quiet:
			print("Reading \"{0}\"...".format(filePath));
		
		inputFile = open(filePath, 'r');
		
		ionicStepCounter = 0;
		
		#COMMENT: Run through the input file line by line, searching for, in the order in which they should appear, the values of ISPIN, LSORBIT and NELECT, the list of k-point coordinates (in reciprocal lattice vectors) and weights, and the Fermi energy.
		#COMMENT: Finding the last of these identifies the start of a block of k-point/band data, and the data is then read in.
		#COMMENT: The use of continue statements rather than a big nested if/else block is probably a bit odd, but I prefer the layout to a huge number of levels of nesting.
		
		kPointCoordinates = [];
		kPointWeights = [];
		
		for line in inputFile:
			match = BandAnalyser.ISPINRegex.search(line);
			
			if match:
				self._isSpinPolarisedCalculation = int(match.group("ispin_value")) == 2;
				continue;
			
			match = BandAnalyser.LSORBITRegex.search(line);
			
			if match:
				self._isSpinOrbitCouplingCalculation = match.group("lsorbit_value") == 'T';
				continue;
			
			match = BandAnalyser.NELECTRegex.search(line);
			
			if match:
				self._numberOfElectrons = int(match.group("nelect_value"));
				continue;
			
			if "k-points in reciprocal lattice and weights" in line:
				line = next(inputFile);
				
				while line.strip() != "":
					elements = SplitLine(line);
					
					if len(elements) == 4:
						kPointCoordinates.append((float(elements[0]), float(elements[1]), float(elements[2])));
						kPointWeights.append(float(elements[3]));
					else:
						raise Exception("Error: Input file does not follow the expected format.");
					
					line = next(inputFile);
				
				self._kPointCoordinates = kPointCoordinates;
				self._kPointWeights = kPointWeights;
				
				continue;
			
			match = BandAnalyser.EfermiRegex.search(line);
			
			if match:
				ionicStepCounter = ionicStepCounter + 1;
				
				if ionicStepCounter > 1:
					print("  -> WARNING: The OUTCAR file contains band data from multiple ionic steps; the BandAnalyser instance will be populated with data from step {0}.".format(ionicStepCounter));
				
				self._fermiEnergy = float(match.group("efermi_value"));
				
				if not quiet:
					if ionicStepCounter == 1:
						print("  -> INFO: Spin Polarised? {0}".format("Yes" if self._isSpinPolarisedCalculation else "No"));
						print("  -> INFO: Spin-Orbit Coupling? {0}".format("Yes" if self._isSpinOrbitCouplingCalculation else "No"));
						print("  -> INFO: Number Of Electrons = {0}".format(self._numberOfElectrons));
					
					print("  -> INFO: Fermi Energy = {0:.3f}".format(self._fermiEnergy));
				
				kPointBandDataUp, kPointBandDataDown = [], [] if self._isSpinPolarisedCalculation else None;
				
				firstSpinComponent = None;
				
				#COMMENT: The end of a block of data is marked by a line of "-" characters - the script will keep processing lines until it encounters one of these.
				
				while not BandAnalyser.SectionEndRegex.match(line):
					if "spin component" in line:
						firstSpinComponent = True if firstSpinComponent == None else False;
					
					if "k-point" in line:
						bandData = [];
						
						#COMMENT: Skip the band data header row, and read in the first line of band data following it.
						
						for i in range(0, 2):
							line = next(inputFile);
						
						#COMMENT: Read the band data lines until a blank line, which marks the end of the block, is encountered.
						
						while line.strip() != "":
							elements = SplitLine(line);
							
							if len(elements) == 3:
								bandData.append((int(elements[0]), float(elements[1]), float(elements[2])));
							else:
								raise Exception("Error: Input file does not follow the expected format.");
							
							line = next(inputFile);
						
						#COMMENT: Append the block of band data to the k-point list for the current spin component.
						
						if firstSpinComponent == None or firstSpinComponent:
							kPointBandDataUp.append(bandData);
						else:
							kPointBandDataDown.append(bandData);
					
					line = next(inputFile);
				
				if not quiet and ionicStepCounter == 1:
					print("  -> INFO: Number Of k-Points = {0}".format(len(kPointCoordinates)));
					print("  -> INFO: Number Of Bands = {0}".format(len(kPointBandDataUp[0])));
				
				
				self._kPointBandData = (kPointBandDataUp, kPointBandDataDown);
		
		inputFile.close();
		
		#COMMENT: Reset the bandgap analysis properties, so that an analysis has been done on a previously-loaded file, it gets cleared.
		
		self._directBandgap = None;
		self._indirectBandgap = None;
		self._kPointValenceConductionBands = None;
		
		if performBandgapAnalysis:
			if not quiet:
				print();
			
			self.PerformBandgapAnalysis(bandgapAnalysisOnlyAnalyseZeroWeightKPoints, bandgapAnalysisPartialOccupancyThreshold, quiet);
	
	def PerformBandgapAnalysis(self, onlyAnalyseZeroWeightKPoints = True, partialOccupancyThreshold = 0.01, quiet = False):
		#COMMENT: Check before starting the analysis that an OUTCAR file has been successfully loaded.
		#COMMENT: This could have been a one-liner, but this more verbose format makes it easier to see what the code block is doing.
		
		notLoaded = False;
		
		notLoaded = notLoaded or self._isSpinPolarisedCalculation == None;
		notLoaded = notLoaded or self._isSpinOrbitCouplingCalculation == None;
		
		notLoaded = notLoaded or self._numberOfElectrons == None;
		
		notLoaded = notLoaded or self._fermiEnergy == None;
		notLoaded = notLoaded or self._kPointCoordinates == None;
		notLoaded = notLoaded or self._kPointWeights == None;
		notLoaded = notLoaded or self._kPointBandData == None;
		
		if notLoaded:
			raise Exception("Error: A bandgap analysis cannot be performed before an OUTCAR file has been successfully loaded.");
		
		if not quiet:
			print("Performing Bandgap Analysis...");
		
		kPointWeights = self._kPointWeights;
		kPointBandDataUp, kPointBandDataDown = self._kPointBandData;
		
		if onlyAnalyseZeroWeightKPoints:
			zeroWeightedKPointIndices = self._GetZeroWeightKPointIndices();
			
			if len(zeroWeightedKPointIndices) > 0:
				if not quiet:
					print("  -> INFO: The k-point list contains some entries with zero weight - bandgap analysis will be performed only with this subset, and k-point indices will be reported with respect to the number of points it contains.")
					print("  -> INFO: Zero-weight k-points are always assigned zero occupancies - at present, the bandgap analysis determines this based on the number of electrons in the system.");
				
				kPointBandDataUp, kPointBandDataDown = [kPointBandDataUp[index] for index in zeroWeightedKPointIndices], [kPointBandDataDown[index] for index in zeroWeightedKPointIndices] if self._isSpinPolarisedCalculation else None;
			else:
				if not quiet:
					print("  -> INFO: The k-point list does not contain any entries with zero weight - all points will be considered during the bandgap analsis.")
		
		valenceConductionBandsUp, valenceConductionBandsDown = [], [] if self._isSpinPolarisedCalculation else None;
		
		for input, output in [(kPointBandDataUp, valenceConductionBandsUp), (kPointBandDataDown, valenceConductionBandsDown)] if self._isSpinPolarisedCalculation else [(kPointBandDataUp, valenceConductionBandsUp)]:
			for i in range(0, len(input)):
				bandData = input[i];
				
				j = len(bandData) - 1;
				if onlyAnalyseZeroWeightKPoints or kPointWeights[i] == 0.0:
					#COMMENT: For k-points with zero weight, VASP does not print out occupations, so assume occupancies based on the values of NELECT and ISPIN/LSORBIT read from the OUTCAR file.
					
					#TODO: This could be handled better - e.g. copying occupations from the nearest non-zero weighted k-point.
					
					cbIndex = int(math.ceil(self._numberOfElectrons / 2.0)) 
                    #if not (self._isSpinPolarisedCalculation or self._isSpinOrbitCouplingCalculation) else self._numberOfElectrons;
					
					vbIndex, vbEnergy, vbOccupation = bandData[cbIndex - 1];
					cbIndex, cbEnergy, cbOccupation = bandData[cbIndex];
					
					vbOccupation = None;
					
					if not (self._isSpinPolarisedCalculation or self._isSpinOrbitCouplingCalculation):
						vbOccupation = 2.0 if self._numberOfElectrons % 2 == 0 else 1.0;							
					else:
						vbOccupation = 1.0;
					
					output.append((vbIndex, vbEnergy, vbOccupation, cbIndex, cbEnergy, cbOccupation, cbEnergy - vbEnergy));
				else:
					while True:
						index, energy, occupation = bandData[j];
						
						if occupation > partialOccupancyThreshold:
							index2, energy2, occupation2 = bandData[j + 1];
							
							output.append((index, energy, occupation, index2, energy2, occupation2, energy2 - energy));
							
							#TODO: Possibly validate the band indices against self._numberOfElectrons to warn about metallic systems/band occupancies being set manually?
							
							break;
						
						j = j - 1;
		
		directBandgapUp, directBandgapDown = None, None;
		indirectBandgapUp, indirectBandgapDown = None, None;
		
		for firstSpinChannel in [True, False] if self._isSpinPolarisedCalculation else [True]:
			valenceConductionBandSet = valenceConductionBandsUp if firstSpinChannel else valenceConductionBandsDown;
			
			vbIndex, vbEnergy, vbOccupation, cbIndex, cbEnergy, cbOccupation, directBandgap = valenceConductionBandSet[0];
			
			minDirectBandgap = directBandgap;
			minDirectBandgapVBIndex = vbIndex;
			minDirectBandgapCBIndex = cbIndex;
			minDirectBandgapKPoint = 0;
			
			for i in range(1, len(valenceConductionBandSet)):
				vbIndex, vbEnergy, vbOccupation, cbIndex, cbEnergy, cbOccupation, directBandgap = valenceConductionBandSet[i];
				
				if directBandgap < minDirectBandgap:
					minDirectBandgap = directBandgap;
					minDirectBandgapVBIndex = vbIndex;
					minDirectBandgapCBIndex = cbIndex;
					minDirectBandgapKPoint = i;
			
			if firstSpinChannel:
				directBandgapUp = (minDirectBandgap, minDirectBandgapVBIndex, minDirectBandgapCBIndex, minDirectBandgapKPoint);
			else:
				directBandgapDown = (minDirectBandgap, minDirectBandgapVBIndex, minDirectBandgapCBIndex, minDirectBandgapKPoint);
			
			if not quiet:
				if self._isSpinPolarisedCalculation:
					print("  -> INFO: Direct bandgap (spin channel {0}) is {1:.3f} eV ({2:.0f} nm) between bands {3} and {4} at k-point {5}".format(1 if firstSpinChannel else 2, minDirectBandgap, 1240.0 / minDirectBandgap if minDirectBandgap != 0.0 else 0.0, minDirectBandgapVBIndex, minDirectBandgapCBIndex, minDirectBandgapKPoint + 1));
				else:
					print("  -> INFO: Direct bandgap is {0:.3f} eV ({1:.0f} nm) between bands {2} and {3} at k-point {4}".format(minDirectBandgap, 1240.0 / minDirectBandgap if minDirectBandgap != 0.0 else 0.0, minDirectBandgapVBIndex, minDirectBandgapCBIndex, minDirectBandgapKPoint + 1));
			
			if len(valenceConductionBandSet) > 1:
				vbIndex, vbEnergy, vbOccupation, cbIndex, cbEnergy, cbOccupation, directBandgap = valenceConductionBandSet[0];
				
				maxVBEnergy = vbEnergy;
				maxVBIndex = vbIndex;
				maxVBKPoint = 0;
				
				minCBEnergy = cbEnergy;
				minCBIndex = cbIndex;
				minCBKPoint = 0;
				
				vbIndex, vbEnergy, vbOccupation, cbIndex, cbEnergy, cbOccupation, directBandgap = valenceConductionBandSet[1];
				
				maxVBEnergy2 = vbEnergy;
				maxVBIndex2 = vbIndex;
				maxVBKPoint2 = 1;
				
				minCBEnergy2 = cbEnergy;
				minCBIndex2 = cbIndex;
				minCBKPoint2 = 1;
				
				for i in range(1, len(valenceConductionBandSet)):
					vbIndex, vbEnergy, vbOccupation, cbIndex, cbEnergy, cbOccupation, directBandgap = valenceConductionBandSet[i];
					
					if vbEnergy > maxVBEnergy:
						maxVBEnergy2 = maxVBEnergy;
						maxVBIndex2 = maxVBIndex;
						maxVBKPoint2 = maxVBKPoint;
						
						maxVBEnergy = vbEnergy;
						maxVBIndex = vbIndex;
						maxVBKPoint = i;
					elif vbEnergy > maxVBEnergy2:
						maxVBEnergy2 = vbEnergy;
						maxVBIndex2 = vbIndex;
						maxVBKPoint2 = i;
					
					if cbEnergy < minCBEnergy:
						minCBEnergy2 = minCBEnergy;
						minCBIndex2 = minCBIndex;
						minCBKPoint2 = minCBKPoint;
						
						minCBEnergy = cbEnergy;
						minCBIndex = cbIndex;
						minCBKPoint = i;
					elif cbEnergy < minCBEnergy2:
						minCBEnergy2 = cbEnergy;
						minCBIndex2 = cbIndex;
						minCBKPoint2 = i;
				
				indirectBandgap = None;
				
				if maxVBKPoint != minCBKPoint:
					indirectBandgap = (minCBEnergy - maxVBEnergy, maxVBIndex, maxVBKPoint, minCBIndex, minCBKPoint);
				else:
					if minCBEnergy - maxVBEnergy2 <= minCBEnergy2 - maxVBEnergy:
						indirectBandgap = (minCBEnergy - maxVBEnergy2, maxVBIndex2, maxVBKPoint2, minCBIndex, minCBKPoint);
					else:
						indirectBandgap = (minCBEnergy2 - maxVBEnergy, maxVBIndex, maxVBKPoint, minCBIndex2, minCBKPoint2);	
			
				if firstSpinChannel:
					indirectBandgapUp = indirectBandgap;
				else:
					indirectBandgapDown = indirectBandgap;
				
				if not quiet:
					bandgap, maxVBIndex, maxVBKPoint, minCBIndex, minCBKPoint = indirectBandgap;
					
					if self._isSpinPolarisedCalculation:
						print("  -> INFO: Indirect bandgap (spin channel {0}) is {1:.3f} eV ({2:.0f} nm) between band {3} at k-point {4} and band {5} at k-point {6}".format(1 if firstSpinChannel else 2, bandgap, 1240 / bandgap if bandgap != 0.0 else 0.0, maxVBIndex, maxVBKPoint + 1, minCBIndex, minCBKPoint + 1));
					else:
						print("  -> INFO: Indirect bandgap is {0:.3f} eV ({1:.0f} nm) between band {2} at k-point {3} and band {4} at k-point {5}".format(bandgap, 1240 / bandgap if bandgap != 0.0 else 0.0, maxVBIndex, maxVBKPoint + 1, minCBIndex, minCBKPoint + 1));
			else:
				if not quiet and firstSpinChannel:
					print("  -> INFO: Only one k-point - cannot calculate an indirect bandgap");
		
		self._directBandgap = (directBandgapUp, directBandgapDown);
		self._indirectBandgap = (indirectBandgapUp, indirectBandgapDown);
		self._kPointValenceConductionBands = (valenceConductionBandsUp, valenceConductionBandsDown);
	
	def SaveKPointData(self, filePath, quiet = False):
		fermiEnergy = self._fermiEnergy;
		kPointCoordinates = self._kPointCoordinates;
		
		if fermiEnergy == None or kPointCoordinates == None:
			raise Exception("Error: k-point data cannot be saved before an OUTCAR file has been successfully loaded.");
		
		kPointValenceConductionBandsUp, kPointValenceConductionBandsDown = self._kPointValenceConductionBands if self._kPointValenceConductionBands != None else (None, None);
		
		if not quiet:
			print("Writing k-point data to \"{0}\"...".format(os.path.split(filePath)[1]));
			
			if kPointValenceConductionBandsUp != None:
				print("  -> INFO: {0}".format("A bandgap analysis has been performed, so information about the direct bandgap at each k-point will be saved for both spin channels." if kPointValenceConductionBandsDown != None else "A bandgap analysis has been performed, so information about the direct bandgap at each k-point will be saved."));
				
				if len(kPointValenceConductionBandsUp) != len(kPointCoordinates):
					print("  -> INFO: The bandgap analysis was performed with only zero-weighted k-points - only these points will be output.");
					
					kPointCoordinates = [kPointCoordinates[index] for index in self._GetZeroWeightKPointIndices()];
		
		outputFile = open(filePath, 'w');
		outputFileCSV = csv.writer(outputFile, delimiter = ',', quotechar = '\"', quoting = csv.QUOTE_ALL);
		
		headerRow = ["Index", "kx", "ky", "kz"];
		
		if kPointValenceConductionBandsUp != None:
			headerRow = headerRow + (["VB (1)", "EV (1) / eV", "VB Occ (1)", "CB (1)", "EC (1) / eV", "CB Occ (1)", "Eg (Direct, 1) / eV", "VB (2)", "EV (2) / eV", "VB Occ (2)", "CB (2)", "EC (2) / eV", "CB Occ (2)", "Eg (Direct, 2) / eV"] if kPointValenceConductionBandsDown != None else ["VB", "EV / eV", "VB Occ", "CB", "EC / eV", "CB Occ", "Eg (Direct) / eV"]);
		
		outputFileCSV.writerow(headerRow);
		
		for i in range(0, len(kPointCoordinates)):
			kx, ky, kz = kPointCoordinates[i];
			
			dataRow = [i + 1, kx, ky, kz];
			
			if kPointValenceConductionBandsUp != None:
				vbIndex, vbEnergy, vbOccupation, cbIndex, cbEnergy, cbOccupation, directBandgap = kPointValenceConductionBandsUp[i];
				dataRow = dataRow + [vbIndex, vbEnergy - fermiEnergy, vbOccupation, cbIndex, cbEnergy - fermiEnergy, cbOccupation, directBandgap];
				
				if kPointValenceConductionBandsDown != None:
					vbIndex, vbEnergy, vbOccupation, cbIndex, cbEnergy, cbOccupation, directBandgap = kPointValenceConductionBandsDown[i];
					dataRow = dataRow + [vbIndex, vbEnergy - fermiEnergy, vbOccupation, cbIndex, cbEnergy - fermiEnergy, cbOccupation, directBandgap];
			
			outputFileCSV.writerow(dataRow);
		
		outputFile.close();
	
	def SaveBandData(self, filePath, onlyOutputZeroWeightKPoints = True, quiet = False):
		fermiEnergy = self._fermiEnergy;
		kPointWeights = self._kPointWeights;
		kPointBandData = self._kPointBandData;
		
		if fermiEnergy == None or kPointWeights == None or kPointBandData == None:
			raise Exception("Error: Band data cannot be saved before an OUTCAR file has been successfully loaded.");
		
		if not quiet:
			print("Writing band data to \"{0}\"...".format(os.path.split(filePath)[1]));
		
		kPointBandDataUp, kPointBandDataDown = kPointBandData;
		
		if onlyOutputZeroWeightKPoints:
			zeroWeightKPoints = self._GetZeroWeightKPointIndices();
			
			if len(zeroWeightKPoints) > 0:
				if not quiet:
					print("  -> INFO: Only band data for k-points with zero weight will be output.");
					print("  -> INFO: Note that VASP does not calculate the occupation of these points, so all occupations will appear as zeros in the output file.");
				
				kPointBandDataUp, kPointBandDataDown = [kPointBandDataUp[index] for index in zeroWeightKPoints], [kPointBandDataDown[index] for index in zeroWeightKPoints] if kPointBandDataDown != None else None;
			else:
				if not quiet:
					print("  -> INFO: The k-point list does not contain any zero-weight points - data for all k-points will be output.");
		
		outputFile = open(filePath, 'w');
		outputFileCSV = csv.writer(outputFile, delimiter = ',', quotechar = '\"', quoting = csv.QUOTE_ALL);
		
		headerRow1 = [""];
		headerRow2 = ["Band"];
		
		if kPointBandDataDown != None:
			for i in range(0, len(kPointBandDataUp)):
				headerRow1 = headerRow1 + ["KPoint {0}".format(i + 1), "", "", ""];
				headerRow2 = headerRow2 + ["E (1) - EF / eV", "Occ (1)", "E (2) - EF / eV", "Occ (2)"];
		else:
			for i in range(0, len(kPointBandDataUp)):
				headerRow1 = headerRow1 + ["KPoint {0}".format(i + 1), ""];
				headerRow2 = headerRow2 + ["E - EF / eV", "Occ"];
		
		outputFileCSV.writerow(headerRow1);
		outputFileCSV.writerow(headerRow2);
		
		for i in range(0, len(kPointBandDataUp[0])):
			dataRow = [i + 1];
			
			for j in range(0, len(kPointBandDataUp)):
				bandIndex, bandEnergy, occupation = kPointBandDataUp[j][i];
				dataRow = dataRow + [bandEnergy - fermiEnergy, occupation];
				
				if kPointBandDataDown != None:
					bandIndex, bandEnergy, occupation = kPointBandDataDown[j][i];
					dataRow = dataRow + [bandEnergy - fermiEnergy, occupation];
			
			outputFileCSV.writerow(dataRow);
		
		outputFile.close();
	
	#COMMENT: StaticFields
	
	ISPINRegex = re.compile("ISPIN\s*=\s*(?P<ispin_value>\d)");
	LSORBITRegex = re.compile("LSORBIT\s*=\s*(?P<lsorbit_value>[TF])");
	NELECTRegex = re.compile("NELECT\s*=\s*(?P<nelect_value>\d+)");
	
	EfermiRegex = re.compile("E-fermi\s*\:\s*(?P<efermi_value>[+-]?(\d+)?\.\d+)");
	DecimalRegex = re.compile("(?P<value>[+-]?\d+\.\d+)")
	
	SectionEndRegex = re.compile("^-+$");
	
	#COMMENT: StaticMethods
	
def SplitLine(line):
	elements = [];
		
	for item in line.strip().split(' '):
		strippedItem = item.strip();
			
		if strippedItem != "":
			elements.append(strippedItem);
		
	return elements;

#EndRegion


#Region: Main

if __name__ == "__main__":
	print("BandEnergies.py");
	print();
	
	print("WARNING: This code is supplied as is, and in particular is not guaranteed to be bug free - please check the results carefully, especially if they are unexpected.");
	print();
	
	bandAnalyser = BandAnalyser();
	
	bandAnalyser.ReadOUTCARFile(InputFile, bandgapAnalysisPartialOccupancyThreshold = PartialOccupancyThreshold);
	print();
	
	prefix = Prefix;
	
	if prefix == None or prefix == "":
		head, tail = os.path.split(InputFile);
		root, ext = os.path.splitext(tail);
		
		prefix = root;
	
	bandAnalyser.SaveKPointData("{0} - K-Points.csv".format(prefix));
	print();
	
	bandAnalyser.SaveBandData("{0} - Band Energies.csv".format(prefix));
	print();
	
	print("Done!");
	print();

#EndRegion
