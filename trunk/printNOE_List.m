function printNOE_List(NOES, inputModeIndex, inputModelIndex)

fid = fopen('NOE_List.txt','w');
fprintf(1, 'printing to %s\n','NOE_List.txt');

fid2 = fopen('rawNOE_List_peakIndices.txt','w');
fprintf(1, 'printing to %s\n','rawNOE_List_peakIndices.txt');


for peak1Index = 1:size(NOES,1)
    listOfNOEs = find(NOES(peak1Index,:),1);
    if (~isempty(listOfNOEs))
        printThisPeak             = 0;
        listOfPeaksToBePrinted    = [];
        for peak2Index = 1:size(NOES,2)
            if (peak2Index <= peak1Index)
                continue;
            end
            if (NOES(peak1Index,peak2Index) == 1)
                if (printThisPeak == 0)
                    listOfPeaksToBePrinted = [peak1Index];
                    printThisPeak          = 1;
                end
                listOfPeaksToBePrinted = [listOfPeaksToBePrinted peak2Index];
                fprintf(fid2, '%d %d\n', peak1Index, peak2Index);
            end
        end
        
        if (printThisPeak)
            assert (length(listOfPeaksToBePrinted)>1);
            fprintf(fid, '%d %d ', length(listOfPeaksToBePrinted)-1, listOfPeaksToBePrinted(1));
            for peaksToBePrintedIndex = 2:length(listOfPeaksToBePrinted)
                fprintf(fid, '%d ', listOfPeaksToBePrinted(peaksToBePrintedIndex));
            end
            fprintf(fid, '\n');
        end
    end
end

fclose(fid);
fclose(fid2);
