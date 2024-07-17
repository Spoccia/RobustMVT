close all;
clc;
clear;

PATH_dataset ='/Users/sicongliu/Desktop/data/mocap/';
PATH_coordinates='/Users/sicongliu/Desktop/LocationMatrix/';
saveFeaturesPath='/Users/sicongliu/Desktop/features/';
%% sift parameters
% x - variate
% y - time
% oframes - octaves
% sigmad - sigma dependency (variate)
% sigmat - sigma time (time)
% pricur - principle curvature
DeOctTime = 2;
DeOctDepd = 2;
DeLevelTime = 4;%6;%
DeLevelDepd = 4;%6;%
DeSigmaDepd = 0.4;%0.4;%0.6;%0.4;%0.5;%
DeSigmaTime = 4*sqrt(2);%1.6*2^(1/(DeLevelTime));%(1.6*2^(1/DeLevelTime))/2;%
%4*sqrt(2);%2.5*2^(1/DeLevelTime);%1.6*2^(1/DeLevelTime);%4*sqrt(2);%2*1.6*2^(1/DeLevelTime);%  8;%4*sqrt(2);%1.2*2^(1/DeLevelTime);%
thresh = 0.04 / (DeLevelTime) / 2 ;%0.04;%
DeGaussianThres = 0.1;%0.001;%0.7;%0.3;%1;%0.6;%2;%6; % TRESHOLD with the normalization of hte distance matrix should be  between 0 and 1
DeSpatialBins = 4; %NUMBER OF BINs
r= 10; %5 threshould variates

featureExtractionGaussian = zeros(1, 184);
for TEST =1:184
    
    data = csvread ([PATH_dataset,num2str(TEST),'.csv']);%
    coordinates=csvread(strcat(PATH_coordinates,'LocationMatrixMocap.csv'))';
    p = tic;
    sBoundary=1;
    eBoundary=size(data',1);
    [frames1,descr1,gss1,dogss1,depd1,idm1, time, timee, timeDescr] = sift_gaussianSmooth_Silv_Xiaolan(data,coordinates', DeOctTime, DeOctDepd,...
        DeLevelTime, DeLevelDepd, DeSigmaTime ,DeSigmaDepd,...
        DeSpatialBins, DeGaussianThres, r, sBoundary, eBoundary);
    
    while(size(frames1,2)==0)
        frames1 = zeros(4,1);
        descr2 = zeros(128,1);
    end
    frame1 = [frames1;descr1];
    
    feature = frame1;
    
    featureExtractionGaussian(1, TEST) = featureExtractionGaussian(1, TEST) + toc(p);
    
    savepath1 = [saveFeaturesPath,'feature_',TS_name,'.mat'];
    savepath2 = [saveFeaturesPath,'idm_',TS_name,'.mat'];
    savepath3 = [saveFeaturesPath,'MetaData_',TS_name,'.mat'];
    savepath5 = [saveFeaturesPath,'GaussianSmoothing/DepdMatrix_',TS_name,'.mat'];
    
    savepath6 = [saveFeaturesPath,'/ComparisonTime_',TS_name,'.csv'];
    savepath7 = [saveFeaturesPath,'/ScaleTime_',TS_name,'.csv'];
    savepath8 = [saveFeaturesPath,'/DescrTime_',TS_name,'.csv'];
    
    save(savepath1,'data', 'gss1', 'frame1','depd1');
    save(savepath2,'idm1');
    save(savepath3,'DeOctTime', 'DeOctDepd', 'DeSigmaTime','DeSigmaDepd', 'DeLevelTime','DeLevelDepd', 'DeGaussianThres', 'DeSpatialBins', 'r', 'descr1' );
    save(savepath5, 'depd1');
    
    csvwrite(savepath5, time);
    csvwrite(savepath6, timee);
    csvwrite(savepath7, timeDescr);
    
    
    clear frames2 descr2 gss2 dogss2 depd2 idm2 I2;
end