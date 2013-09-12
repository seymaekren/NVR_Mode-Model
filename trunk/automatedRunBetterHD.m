function automatedRunBetterHD

%normalModeHD: This script runs HD on 1EF1 homology models. 

dbstop if error
dbstop if warning


addpath('/Users/student/Desktop/NVR')
addpath('/Users/student/Desktop/NVR/trunk')
addpath('/Users/student/Desktop/NVR/InputFiles/1EF1')
addpath('/Users/student/Desktop/NVR/OutputFiles/1EF1')

modeIndex        = inputModeIndex;

modelIndex       = inputModelIndex;

outputFilename   = sprintf('Mode%d_Model%d_HD_vs_NMA.txt', ...
			   modeIndex, modelIndex);			
  
if (checkIfFileExists(outputFilename))
  error ('please erase the %s\n',outputFilename);
end

fout                    = fopen(outputFilename, 'w');

fprintf(fout, 'useOrigData = %d\n', useOrigData);

SHIFTX_Filename   = sprintf('SHIFTX/MySHIFTX.%d.model%d',modeIndex, modelIndex);

SHIFTS_Filename   = sprintf('SHIFTS/MySHIFTS.%d.model%d',modeIndex, modelIndex);

%bidirectionals.
if (modeIndex == 78)
    modelDataFilename       = sprintf('ModelDataGeneration/ModelDataFiles/Mode.7.8.coords%d',  modelIndex);
else
    modelDataFilename       = sprintf('ModelDataGeneration/ModelDataFiles/Mode%d.coords%d',  modeIndex, modelIndex);
end


[VECTORS,TYPES,RESNUMS,SSTRUCT, HBOND, ALLDISTS,IALLDISTS] = loadmodeldata(modelDataFilename);

[HD_SCORE, assignmentAccuracy, assignments] = ...%HD
    NVR(HSQCDATA(:,2:3),HSQCDATA(:,4:5),HSQCDATA(:,6), ...
	HSQCDATA(:,1),NOES,VECTORS,TYPES,...
	RESNUMS,SSTRUCT, HBOND, ALLDISTS, IALLDISTS, ...
	SHIFTS_Filename, SHIFTX_Filename);




fclose(fout);



