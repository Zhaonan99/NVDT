clc
clear
%%
filename_protein_sequence='D:\matlab\Sequence.txt'
[Header, Se]=fastaread(filename_protein_sequence);
R=length(Se);
features_nucleotide(filename_protein_sequence);
features_dinucleotide(filename_protein_sequence);
features_triplet_nucleotide(filename_protein_sequence);
A1=csvread ('nature_vector_nucleotide.csv',0,0,[0 0 R-1 11]); 
A2=csvread ('nature_vector_2bases.csv',0,0,[0 0 R-1 15]);
A3=csvread ('nature_vector_3bases.csv',0,0,[0 0 R-1 63]);
A0=num2cell([A1,A2,A3]);
ID=Header';
G=[ID,A0];
T=cell2table(G,'VariableNames',{'Pro_ID','nA','nC','nG','nT','uA','uC','uG','uT','D2A','D2C','D2G','D2T','nAA','nAC','nAG','nAT','nCA','nCC','nCG','nCT','nGA','nGC','GG','nGT','nTA','nTC','nTG','nTT','nAAA','nAAC','nAAG','nAAT','nACA','nACC','nACG','nACT','nAGA','nAGC','nAGG','nAGT','nATA','nATC','nATG','nATT','nCAA','nCAC','nCAG','nCAT','nCCA','nCCC','nCCG','nCCT','nCGA','nCGC','nCGG','nCGT','nCTA','nCTC','nCTG','nCTT','nGAA','nGAC','nGAG','nGAT','nGCA','nGCC','nGCG','nGCT','nGGA','nGGC','nGGG','nGGT','nGTA','nGTC','nGTG','nGTT','nTAA','nTAC','nTAG','nTAT','nTCA','nTCC','nTCG','nTCT','nTGA','nTGC','nTGG','nTGT','nTTA','nTTC','nTTG','nTTT'});
writetable(T,'pro_data.xlsx')

%%
filename_protein_pairs='D:\matlab\protein_pairs_ID.xlsx';

[data,str] = xlsread(filename_protein_pairs);
pairs=size(str,1)
F0=[];
for i=1:pairs
    F1=[];
    F2=[];
    for k=0:R-1
        if strcmp(G(k+1),str(i,1))
            G1=G(k+1,2:93);          
        end
    end       
    for  h=0:R-1
        if strcmp(G(h+1),str(i,2))                
            G2=G(h+1,2:93);
        end
    end
F1={G1,G2};
F2={G2,G1};
f1=[F1{:}];
f2=[F2{:}];
f=[f1;f2];
F0=[F0;f];
end
xlswrite('protein_pairs_features.xls',F0);

