addpath('/net/mist/apaydin/ProteinDynamics/Matlab')

outputFilename = 'templateHDResult.txt';
SHIFTX_Filename = 'templateSHIFTX';
SHIFTS_Filename = 'templateSHIFTS';
MODEL_Filename = 'templateModelData';
runHD(outputFilename,  MODEL_Filename, SHIFTS_Filename, SHIFTX_Filename);
