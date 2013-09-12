function automatedRunModeModel

%normalModeHD: This script runs HD on 1EF1 homology models. 

dbstop if error
dbstop if warning


addpath('/Users/student/Desktop/NVR')
addpath('/Users/student/Desktop/NVR/trunk')
addpath('/Users/student/Desktop/NVR/InputFiles/1EF1')
addpath('/Users/student/Desktop/NVR/OutputFiles/1EF1')

modeIndex                   = [7, 8, 9, 10, 11, 78];
modelIndex                  = (1 :11 );
modelIndexBidirectional     = (1 :121 );
		
  
for i = modeIndex;
    if(i == 78)
        for j = modelIndexBidirectional;
            runBetter(i,j);
        end
    else
        for j = modelIndex;
            runBetter(i,j);
        end
    end
end




