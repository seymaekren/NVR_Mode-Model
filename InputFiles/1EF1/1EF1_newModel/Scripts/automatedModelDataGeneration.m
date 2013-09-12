addpath('/home/home4/apaydin/Mist/Matlab/FileProcessing')
modelDataFilename = '../templateModelData.m';
cd Pdb
for modeIndex = 7:11 
  parsedPdbFilename        = sprintf('1EF1-NMA.%d.parsedPDB', modeIndex); 
  matlabFilename           = sprintf('Mode%d.all',modeIndex); 
  generateOneModelData (modelDataFilename, parsedPdbFilename, matlabFilename); 
end 
