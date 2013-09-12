addpath('/home/home4/apaydin/Mist/Matlab/FileProcessing')
modelDataFilename = '../templateModelData.m';
cd Pdb
parsedPdbFilename        = '1EF1.7.8.bidirectional.parsedPDB'; 
matlabFilename           = 'Mode7.8.all'; 
generateOneModelData (modelDataFilename, parsedPdbFilename, matlabFilename); 
