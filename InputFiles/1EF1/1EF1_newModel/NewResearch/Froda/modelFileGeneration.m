addpath('/home/home4/apaydin/Mist/Matlab/FileProcessing')

dbstop if error;

modelDataFilename = '../../templateModelData.m';

for (modelIndex = 100:100:10000)

  modelIndex
  
  if (modelIndex < 1000)
    
    fileIndex                = sprintf('1EF1-NMA.7.model6_froda_00000%d.pdb',modelIndex);
    parsedPdbFilename        = sprintf('1EF1-NMA.7.model6_froda_00000%d.parsedPDB', ...
				       modelIndex); 
    
  elseif (modelIndex < 10000)

      fileIndex              = sprintf('1EF1-NMA.7.model6_froda_0000%d.pdb',modelIndex);
      parsedPdbFilename      = sprintf('1EF1-NMA.7.model6_froda_0000%d.parsedPDB', ...
				       modelIndex); 
      
  else
    
    fileIndex                = sprintf('1EF1-NMA.7.model6_froda_000%d.pdb',modelIndex);
    parsedPdbFilename        = sprintf('1EF1-NMA.7.model6_froda_000%d.parsedPDB', ...
				       modelIndex); 
    
  end
  
  matlabFilename             = sprintf('Model%d.coords',modelIndex); 
  
  generateOneModelData (modelDataFilename, parsedPdbFilename, matlabFilename); 
  
end
