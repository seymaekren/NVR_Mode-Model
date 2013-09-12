addpath('/home/home4/apaydin/Mist/Matlab/FileProcessing');
templateShiftxFile = '../../templateSHIFTX.m'; 
cd SHIFT__FileGeneration/SHIFTX/
for modeIndex = 7:11 
  for modelIndex = 12:21 
    shiftxFilename = sprintf('MySHIFTX.%d.model%d.parsed',modeIndex,modelIndex); 
    adjustAndSSE_InfoToSHIFTX_File(templateShiftxFile, shiftxFilename); 
  end 
end 
