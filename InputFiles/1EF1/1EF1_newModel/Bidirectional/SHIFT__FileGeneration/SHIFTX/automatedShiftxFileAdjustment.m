addpath('/home/home4/apaydin/Mist/Matlab/FileProcessing');
templateShiftxFile = '../../templateSHIFTX.m'; 
cd SHIFT__FileGeneration/SHIFTX/
for modelIndex = 1:121 
  shiftxFilename = sprintf('MySHIFTX.7.8.bidirectional.model%d.parsed',modelIndex); 
  adjustAndSSE_InfoToSHIFTX_File(templateShiftxFile, shiftxFilename); 
end 
