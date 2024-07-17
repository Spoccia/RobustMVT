clear;
clc;
FeatureName='feature_';%'feature';%

DestFolder='D:\Mocap _ RMT2\';
sigmad = 4:6;
for k=1:1%1:3
DestName= ['Features 3 octave  SD 0_',num2str(sigmad(2)), ' ST 2_8'];%['Features Octave3 sd 0_',num2str(sigmad(2)), ' st 2_8'];%
FeaturePath=['D:\Mocap _ RMT2\',DestName,'\'];
datasetsize=184;
N_OD = 3;
N_OT = 3;
AllDatastFeatures=[];
for i=1:datasetsize
    load([FeaturePath,FeatureName,num2str(i),'.mat']);
    AllDatastFeatures=[AllDatastFeatures,frame1];
end


HowManyFeaturesForOctave=[];
for OT =1 : N_OT
 %   OD=OT;
    for OD =1:N_OD
        index_Of_Octave= (AllDatastFeatures(5,:)==OD & AllDatastFeatures(6,:)==OT );
        NumFOctave = sum(index_Of_Octave);
        HowManyFeaturesForOctave= [HowManyFeaturesForOctave,[OT;OD;NumFOctave;NumFOctave/datasetsize]];
    end
end

xlswrite([DestFolder,'Count_',DestName,'.xls'],HowManyFeaturesForOctave);
end