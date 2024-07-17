close all;
clc;
clear;

PATH_dataset ='E:\RMT\Dataset\BirdSong\data\';
%PATH_coordinates='D:\Mocap _ RMT2\Location matrix\';
saveFeaturesPath='E:\RMT\Dataset\BirdSong\Features\';
%% sift parameters
% x - variate
% y - time
% oframes - octaves
% sigmad - sigma dependency (variate)
% sigmat - sigma time (time)
% pricur - principle curvature
DeOctTime = 3;
DeOctDepd = 3;
DeLevelTime = 4;%6;%
DeLevelDepd = 4;%6;%
DeSigmaDepd = 0.5;%0.4;%0.6;
DeSigmaTime = 1.6149;%4*sqrt(2)/2;%3.2298;%4*sqrt(2);%1.6*2^(1/(DeLevelTime));%(1.6*2^(1/DeLevelTime))/2;%1.6*2^(1/(DeLevelTime));%4*sqrt(2)/2;%
%4*sqrt(2);%2.5*2^(1/DeLevelTime);%1.6*2^(1/DeLevelTime);%4*sqrt(2);%2*1.6*2^(1/DeLevelTime);%  8;%4*sqrt(2);%1.2*2^(1/DeLevelTime);%
thresh = 0.04 / (DeLevelTime) / 2 ;%0.04;%
DeGaussianThres = 6;%0.1;%0.001;%0.7;%0.3;%1;%0.6;%2;%6; % TRESHOLD with the normalization of hte distance matrix should be  between 0 and 1
DeSpatialBins = 4; %NUMBER OF BINs
r= 10; %5 threshould variates

% manually create location matrix
%% set up location matrix
IDM1 = [1:13];
IDM2 = [2,2,3,3,4,4,1,5,5,6,6,7,7];%[1,1,2,2,3,3,4,5,5,6,6,7,7];%[1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7];
IDM3 = [2,2,3,1,3,4,4];%[1,1,2,3,2,4,4];%[1, 2, 3, 3, 4, 4];
idm2{1} = IDM1;
idm2{2} = IDM2;
idm2{3} = IDM3;
%%% first octave location matrix
LocM1 = zeros(13, 13);
for i = 1 : 13
    LocM1(i, i) = 1;
    if(i == 1)
        LocM1(i, i+1) = 1;
        LocM1(i+1, i) = 1;
%         LocM1(i+2, i) = 1;%silv
%         LocM1(i, i+2) = 1;%silv
    elseif(i == 13)
        LocM1(i, i-1) = 1;
        LocM1(i-1, i) = 1;
%         LocM1(i, i-2) = 1;%silv
%         LocM1(i-2, i) = 1;%silv
    else
        LocM1(i-1, i) = 1;
        LocM1(i+1, i) = 1;
        LocM1(i, i-1) = 1;
        LocM1(i, i+1) = 1;
%         if(i>2 & i<12)%silv
%             LocM1(i-2, i) = 1;
%             LocM1(i+2, i) = 1;
%             LocM1(i, i-2) = 1;
%             LocM1(i, i+2) = 1;
%         end
    end
end
LocM1_1 =  [1	1	0	0	0	0	0	0	0	0	0	0	0
            1	1	0	0	0	0	0	0	0	0	0	0	0
            0	0	1	1	0	0	0	0	0	0	0	0	0
            0	0	1	1	0	0	0	0	0	0	0	0	0
            0	0	0	0	1	1	0	0	0	0	0	0	0
            0	0	0	0	1	1	0	0	0	0	0	0	0
            0	0	0	0	0	0	1	0	0	0	0	0	0
            0	0	0	0	0	0	0	1	1	0	0	0	0
            0	0	0	0	0	0	0	1	1	0	0	0	0
            0	0	0	0	0	0	0	0	0	1	1	0	0
            0	0	0	0	0	0	0	0	0	1	1	0	0
            0	0	0	0	0	0	0	0	0	0	0	1	1
            0	0	0	0	0	0	0	0	0	0	0	1	1
            ];
LocM1_2 =  [1	1	1	0	0	0	0	0	0	0	0	0	0
            1	1	1	0	0	0	0	0	0	0	0	0	0
            1	1	1	0	0	0	0	0	0	0	0	0	0
            0	0	0	1	1	1	0	0	0	0	0	0	0
            0	0	0	1	1	1	0	0	0	0	0	0	0
            0	0	0	1	1	1	0	0	0	0	0	0	0
            0	0	0	0	0	0	1	0	0	0	0	0	0
            0	0	0	0	0	0	0	1	1	1	0	0	0
            0	0	0	0	0	0	0	1	1	1	0	0	0
            0	0	0	0	0	0	0	1	1	1	0	0	0
            0	0	0	0	0	0	0	0	0	1	1	1	1
            0	0	0	0	0	0	0	0	0	0	1	1	1
            0	0	0	0	0	0	0	0	0	0	1	1	1
            ];
LocM1_3 =  [1	1	1	1	0	0	0	0	0	0	0	0	0
            1	1	1	1	0	0	0	0	0	0	0	0	0
            1	1	1	1	1	1	0	0	0	0	0	0	0
            1	1	1	1	1	1	0	0	0	0	0	0	0
            0	0	1	1	1	1	0	0	0	0	0	0	0
            0	0	1	1	1	1	0	0	0	0	0	0	0
            0	0	0	0	0	0	1	0	0	0	0	0	0
            0	0	0	0	0	0	0	1	1	1	1	0	0
            0	0	0	0	0	0	0	1	1	1	1	0	0
            0	0	0	0	0	0	0	1	1	1	1	1	1
            0	0	0	0	0	0	0	1	1	1	1	1	1
            0	0	0	0	0	0	0	0	0	1	1	1	1
            0	0	0	0	0	0	0	0	0	1	1	1	1
            ];        
 LocM1 = LocM1_1 - eye([13 13]);

%%% second octave location matrix
LocM2 = zeros(7, 7);
for i = 1 : 7
    LocM2(i, i) = 1;
    if(i == 1)
        LocM2(i, i+1) = 1;
        LocM2(i+1, i) = 1;
%         LocM2(i+2, i) = 1;%silv
%         LocM2(i, i+2) = 1;%silv
    elseif(i == 7)
        LocM2(i, i-1) = 1;
        LocM2(i-1, i) = 1;
%         LocM2(i, i-2) = 1;%silv
%         LocM2(i-2, i) = 1;%silv
    else
        LocM2(i-1, i) = 1;
        LocM2(i+1, i) = 1;
        LocM2(i, i-1) = 1;
        LocM2(i, i+1) = 1;
    end
end
LocM2_1 =  [1	1	0	0	0	0	0
            1	1	0	0	0	0	0
            0	0	1	0	1	0	0
            0	0	0	1	0	0	0
            0	0	1	0	1	0	0
            0	0	0	0	0	1	1
            0	0	0	0	0	1	1
            ];
LocM2_2 =  [1	1	1	0	1	0	0
            1	1	1	0	1	0	0
            1	1	1	0	1	1	1
            0	0	0	1	0	0	0
            1	1	1	0	1	1	1
            0	0	1	0	1	1	1
            0	0	1	0	1	1	1
            ];
% LocM2 = LocM2 - eye([7 7]);
LocM2 = LocM2_1 - eye([7 7]);

LocM3=zeros(4,4);
for i=1:4
    LocM3(i,i)=1;
    if(i==1)
        LocM3(i,i+1)=1;
        LocM3(i+1,i)=1;      
%         LocM3(i+2, i) = 1;%silv
%         LocM3(i, i+2) = 1;%silv
    elseif(i==4)
        LocM3(i,i-1)=1;
        LocM3(i-1,i)=1;
%         LocM3(i,i-2)=1;%silv
%         LocM3(i-2,i)=1;%silv
    else
        LocM3(i,i+1)=1;
        LocM3(i+1,i)=1;
        LocM3(i,i-1)=1;
        LocM3(i-1,i)=1;
%         if(i>2 & i<12)%silv
%         LocM3(i,i+2)=1;
%         LocM3(i+2,i)=1;
%         LocM3(i,i-2)=1;
%         LocM3(i-2,i)=1;
%         end
    end
end
LocM3_1 = [1	1	0	0
           1	1	1	1 
           0	1	1	1
           0	1	1	1
                          ];
LocM3_2 =  [1	1	0	0
            1	1	0	1
            0	0	1	0
            0	1	0	1
            ];         
LocM3_3 =  [1	1	0	1
            1	1	0	1
            0	0	1	0
            1	1	0	1
            ];             
% LocM3 = LocM3 - eye([4 4]);
 LocM3 = LocM3_3-eye([4 4]);
featureExtractionGaussian = zeros(1, 154);
for TEST =1:154
    TS_name=num2str(TEST)
    %load([PATH_dataset,num2str(TEST),'.mat']); 
    data = csvread ([PATH_dataset,num2str(TEST),'.csv']);%
    data(isnan(data))=0;
    %data = data'
%     coordinates=csvread(strcat(PATH_coordinates,'LocationMatrixMocap.csv'))';
    p = tic;
    sBoundary=1;
    eBoundary=size(data,1);
     [frames1,descr1,gss1,dogss1,depd1,idm1, time, timee, timeDescr] = sift_gaussianSmooth_BirdSong(data,...
         LocM1 ,LocM2,LocM3,IDM1, IDM2, IDM3, DeOctTime, DeOctDepd,...
        DeLevelTime, DeLevelDepd, DeSigmaTime ,DeSigmaDepd,...
        DeSpatialBins, DeGaussianThres, r, sBoundary, eBoundary);
%     [frames1,descr1,gss1,dogss1,depd1,idm1, time, timee, timeDescr] = sift_gaussianSmooth_Silv(data,coordinates', DeOctTime, DeOctDepd,...
%         DeLevelTime, DeLevelDepd, DeSigmaTime ,DeSigmaDepd,...
%         DeSpatialBins, DeGaussianThres, r, sBoundary, eBoundary);
     
    while(size(frames1,2)==0)
        frames1 = zeros(4,1);
        descr2 = zeros(128,1);
    end
    frame1 = [frames1;descr1];
    if( isnan(sum(descr1(:))))
        TS_name
       nanIDX=  isnan(sum(descr1));
       frame1(:,nanIDX)= [];
       descr1(:,nanIDX)= [];
       frames1(:,nanIDX)= [];
    end
    feature = frame1;
    
    featureExtractionGaussian(1, TEST) = featureExtractionGaussian(1, TEST) + toc(p);
    
    savepath1 = [saveFeaturesPath,'feature_',TS_name,'.mat'];
    savepath2 = [saveFeaturesPath,'idm_',TS_name,'.mat'];
    savepath3 = [saveFeaturesPath,'MetaData_',TS_name,'.mat'];
    savepath5 = [saveFeaturesPath,'GaussianSmoothing\DepdMatrix_',TS_name,'.mat'];
    
    savepath6 = [saveFeaturesPath,'\ComparisonTime_',TS_name,'.csv'];
    savepath7 = [saveFeaturesPath,'\ScaleTime_',TS_name,'.csv'];
    savepath8 = [saveFeaturesPath,'\DescrTime_',TS_name,'.csv'];
    
    save(savepath1, 'gss1', 'frame1','depd1');%'data',
    save(savepath2,'idm1');
    %save(savepath3,'DeOctTime', 'DeOctDepd', 'DeSigmaTime','DeSigmaDepd', 'DeLevelTime','DeLevelDepd', 'DeGaussianThres', 'DeSpatialBins', 'r', 'descr1' );
    %save(savepath5, 'depd1');
    
    csvwrite(savepath5, time);
    csvwrite(savepath6, timee);
    csvwrite(savepath7, timeDescr);
    
    
    clear frames1 descr1 gss1 dogss1 depd1 idm1 data;
end