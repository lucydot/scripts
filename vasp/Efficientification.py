import os;
import shutil;
import sys;


RelaxationBase = "Relaxation";

MaxSteps = 100;


breakNext = False;

while os.path.exists("INCAR"):
	os.system("aprun -n 240 /homeappl/home/pr1ujms0/appl_sisu/vasp.5.4.1.05Feb16/bin/vasp_gam > vasp.out");
	
	fileList = [];
	lastRelaxationNumber = 0;
	
	for entry in os.listdir("./"):
		if os.path.isfile(entry):
			root, ext = os.path.splitext(entry);
			
			if ext == "" or ext == ".xml":
				fileList.append(entry);
		elif os.path.isdir(entry):
			if RelaxationBase in entry:
				relaxationNumber = int(entry.replace("{0} - Run ".format(RelaxationBase), ""));
				
				if relaxationNumber > lastRelaxationNumber:
					lastRelaxationNumber = relaxationNumber;
	
	if os.path.exists("OUTCAR"):
		jobCompleted = False;
		
		outcarFile = open("OUTCAR", 'r');
		
		for line in outcarFile:
			if "Voluntary context switches" in line:
				jobCompleted = True;
				break;
		
		outcarFile.close();
		
		if not jobCompleted:
			print("WARNING: Job did not complete!");
			break;
	else:
		raise Exception("ERROR: OUTCAR not created - something went horribly wrong!");
	
	stepCount = 0;
	
	if os.path.exists("OSZICAR"):
		oszicarFile = open("OSZICAR", 'r');
		
		for line in oszicarFile:
			if "E0" in line:
				stepCount = stepCount + 1;
		
		oszicarFile.close();
	else:
		raise Exception("ERROR: OSZICAR not created - something probably went horribly wrong!");
	
	relaxationDir = "{0} - Run {1}".format(RelaxationBase, lastRelaxationNumber + 1);
	
	os.mkdir(relaxationDir);
	
	for item in fileList:
		shutil.move(item, os.path.join(relaxationDir, item));
	
	if os.path.exists("vasp.out"):
		os.rename("vasp.out", "{0}.out".format(relaxationDir));
	else:
		raise Exception("ERROR: vasp.out not created - something went horribly wrong!");
	
	print("Run: \"{0}\"".format(relaxationDir));
	
	shutil.copy(os.path.join(relaxationDir, "CONTCAR"), "POSCAR");
	
	for inputFile in ["INCAR", "KPOINTS", "POTCAR"]:
		shutil.copy(os.path.join(relaxationDir, inputFile), inputFile);
	
	#sys.exit(0);
	
	if breakNext and stepCount < MaxSteps:
		break;
	else:
		if stepCount < MaxSteps:
			breakNext = True;
		else:
			breakNext = False;
