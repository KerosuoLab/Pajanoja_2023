 
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%   Maintenance of pluripotency in entire post-gastrulation ectoderm enables neural crest formation
%
%   Ceren Pajanoja, 2022 

%   The following script includes the code used to execute STEP-1 of scMST pipeline

%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%

%% 1) read in Ilastik file (for membrane staining, h5 format)
addpath(genpath('\path\to\data\folder\'));
ilastik_filename = 'SampleName_Probabilities.h5';% type membrane ilastik file name here
ilastik_file = h5read(ilastik_filename,'/exported_data/');
pred = squeeze(ilastik_file(2,:,:,:));                                       
pred = permute(pred,[2,1,3]);          

figure
imshow(sum(pred,3),[])
title('Ilastik Prediction Map View, z-Projection')

%% 2) read in the original membrane image (tiff format)
imagename_original = 'RGB.tif';
original_img = 0*pred;                                    
for z = 1 : size(pred,3)                                         
    temp = imread(imagename_original,z);        
    original_img(:,:,z) = temp(:,:,1);  
end                                                                                       
figure
imshow(sum(original_img,3),[])
title('Membrane Image, z-Projection')

%% 3) Blurred version of the membrane img
im_blurr = imdilate(original_img,strel3D('sphere',3));
figure
imshow(sum(im_blurr,3),[])
title('Blurred Membrane Image, z-Projection')

%% 4) Threshold the Ilastik map pixel value
seg = pred>0.99; 
seg = bwareaopen(seg,50);         % min size 
seg = seg - bwareaopen(seg,3000000); % max size 
figure
imshow(sum(single(seg),3),[])
title('Ilastik map pixel value thresholding, z-Projection')

%% 5) Watershed segmentation 
seed = imimposemin(im_blurr,seg);
Label = watershed(seed);

%% 6) Size hreshold Label matrix and check in Volume viewer
vol_min = 50; 
vol_max = 200000; 
Label2 = bwareaopen(Label,vol_min);             % min size
Label2 = Label2 - bwareaopen(Label2,vol_max);   % max size
Label2 = bwlabeln(Label2);
stats = regionprops3(Label2);
volume = cat(1,stats.Volume);

%% 7) Save Label matrices to disk; 
for z = 1 : size(Label2,3)
     imwrite(Label2(:,:,z),'Label2_WatershedMatrix.tif','compression','none','WriteMode','append');
end
save('Label2_WatershedMatrix.mat','Label2');

%% Visually check the result of watershed segmentation slice by slice
for z = 1 : size(original_img,3)
    
    temp = zeros(size(original_img,1),size(original_img,2),3,'uint8');
    per = Label2(:,:,z) == 0; %
    temp(:,:,1) = original_img(:,:,z);  % original img on the back
    temp(:,:,2) = uint8(per)*10;       % label (segmentation lines)
    imwrite(temp,'seg_RGB.tif','compression','none','WriteMode','append');
end

%% HYB ROUNDS ANALYSIS
%-------------------------------------------------------------------------%

%   Below start analyzing dot (transcript) ilastik files 
%   and signal processing 

%-------------------------------------------------------------------------%

%% 1) Scale image dimensions depending on the scope and objective used

AimedMag  = 0.1135; % pixel size in microns used to analyse the data;

nChan    = 5; % number of channels used
hybnum   = 2; %remember to change here !!
somite   = 7; %number of somites for the embryo
sample   = 'e5-RightTop-3';
scope    = 'dragonfly2';%'olympus';
imageName_allDots = '280621_7som_hyb2_e5-RightTop-3.ome.tif';
scaleFactorX      = 1; % or AimedMag/dx;
scaleFactorZ      = 1;

% settings of the two different types of images. These values might need to
% be changed, depending on the acquisition used.
if strcmp(scope,'dragonfly')
    mag = 63;
    pixelSize = 0.184;     % in microns;
    dx  = pixelSize/mag;
    dz  = 0.356;           % in microns
elseif strcmp(scope,'dragonfly2')
    mag = 63;
    pixelSize = 0.155;   % in microns
    dx  = pixelSize/mag;
    dz  = 0.356;           % in microns 
end

stackSize = (length(imfinfo(imageName_allDots)))/nChan;
im_dots   = cell(1,nChan); % pre-allocating
for c = 1:nChan
    for z = 1:stackSize       
        if c < 6
            im_dots{c}(:,:,z) = (double(imread(imageName_allDots,c+(z-1)*nChan)));
        else
            im_dots{c}(:,:,z) = double(imread(imageName_allDots,c+(z-1)*nChan));
        end
    end
end

% scale the image to the Aimed PixelSize;
im_dots_scaled = cell(size(im_dots));
[X,Y] = meshgrid(1:scaleFactorX:size(im_dots{1},1),1:scaleFactorX:size(im_dots{1},2));
newNSlices = round(size(im_dots{1},3)/scaleFactorZ);
for c = 1 : nChan
    im_dots_scaled{c} = imresize(im_dots{c},1/scaleFactorX);
    scaled = imresize(permute(im_dots_scaled{c},[3,1,2]),[newNSlices,size(permute(im_dots_scaled{c},[3,1,2]),2)]);
    im_dots_scaled{c} = ipermute(scaled,[3,1,2]);
end

%% 1.extra) Write scaled images to disc
% separates channels and saves individually
%-----------------------------------------------------------
% IMPORTANT: Run this only ONE TIME for image with multiple channels 

for c = 1 : nChan
    for z = 1 : size(im_dots_scaled{c},3)
        imwrite(uint16(im_dots_scaled{c}(:,:,z)),sprintf('%dsom_hyb%d_%s_chan_%d.tif',somite,hybnum,sample,c),'tiff','compression','none','WriteMode','append');
    end
end

%%  2) Read in membrane image (unbinned) this is RGB.
imagename_RGB = 'RGB.tif';            % write there the RGB image name
RGB_StackSize = size(original_img,3); 
im_RGB = cell(1,2);                   % pre-allocating
for z = 1 : RGB_StackSize
    temp = imread(imagename_RGB,z);
    im_RGB{1}(:,:,z) = temp(:,:,1);   
end

%% 3.1) Manual rotation from im_dots_scaled into im_RGB frame of reference

chan  = 5; %choose a channel for alignment 
angle = 0; 
im_dots_rot = imrotate(im_dots_scaled{chan},angle,'crop');

[shift,image1,image2] = xcorr3fft(2*im_RGB{1}(:,:,:),uint8(255*mat2gray(im_dots_rot)));

%% 3.2)takes image2, and shifts it according to shift. 
close all; 

image2_shift = 0*image2; % pre-allocating
if shift(1)<0 && shift(2)>0
    image2_shift(-shift(1)+(1:(size(image2,1)+shift(1))),shift(2)+(1:size(image2,2)-shift(2)),:) = ...
        image2(1:(size(image2,1)+shift(1)),1:size(image2,2)-shift(2),:);
elseif shift(1)>0 && shift(2)<0
    image2_shift((1:size(image2,1)-shift(1)),1:size(image2,2)+shift(2),:) = ...
        image2((shift(1)+1):size(image2,1),-shift(2)+(1:size(image2,2)+shift(2)),:);
elseif shift(1)<=0 && shift(2)<=0
    image2_shift(-shift(1)+(1:size(image2,1)+shift(1)),1:size(image2,2)+shift(2),:) = ...
        image2((1:size(image2,1)+shift(1)),-shift(2)+(1:size(image2,2)+shift(2)),:);
elseif shift(1)>=0 && shift(2)>=0
    image2_shift((1:size(image2,1)-shift(1)),shift(2)+(1:size(image2,2)-shift(2)),:) = ...
         image2((shift(1)+1):size(image2,1),1:size(image2,2)-shift(2),:);
elseif shift(1) == 0 && shift(2) == 0
    image2_shift = image2;
end

% visualize the result of the shift

tmp        = zeros(size(image1,1),size(image1,2),3,'double');
tmp(:,:,1) = double(sum(image1,3))/size(image1,3);             % (membrane) mean intensity projection
tmp(:,:,2) = (sum(double(image2_shift)/size(image2,3),3))*5;  % green (shifted version of image2)
tmp        = uint8(255*mat2gray(tmp));
imshow(tmp);

%type in command window -> shift = [x,y,z] in order to adjust the shift and run this section again
% x (+up, -down) & y (+right, -left)

%%
%---------------------------------------
% If alignment is good, move to step 4
%---------------------------------------
%
%% 4) Summarize the dots in each cell that were found using waterhed label

% First segment the dots using ilastk prediction map.
% Points: summarizes information about cells, such as Volume, centre of mass
% position (centroid), Mean Intensity in the segmentation of the dots,
% total intensity; To get the total number of spots multiply by volume;
points = struct('volume',[],'centroid',[],'MeanIntensity',[],'Intensity',[],...
    'NSpots',[],'NSpotsinVol',[],'NSpots_Cleared',[],'NSpotsinVol_Cleared',[]);

%% 5) Run this for each channel

for chan = 1
    ilastik_file_dots = h5read(sprintf('%dsom_hyb%d_%s_chan_%d_Probabilities.h5',somite,hybnum,sample,chan),'/exported_data/');
    pred = squeeze(ilastik_file_dots(2,:,:,:));  
    pred = permute(pred,[2,1,3]);
    image2 = pred>0.5;                      % threshold value(can be low typically .5 or .3)
    image2 = imrotate(image2,angle,'crop'); % rotates the segmentation 
    
    si1 = size(image1);
    si2 = size(image1);
    desiredSize = max(si1,si2);
    padded_image = zeros(desiredSize,'uint8');
    padded_image(1:size(image2,1),1:size(image2,2),1:size(image2,3)) = image2;
    image2 = padded_image;
    
    %shift the segmentation;
    image2_shift = 0*image2;
    if shift(1)<0 && shift(2)>0
        image2_shift(-shift(1)+(1:(size(image2,1)+shift(1))),shift(2)+(1:size(image2,2)-shift(2)),:) = ...
            image2(1:(size(image2,1)+shift(1)),1:size(image2,2)-shift(2),:);
    elseif shift(1)>0 && shift(2)<0
        image2_shift((1:size(image2,1)-shift(1)),1:size(image2,2)+shift(2),:) = ...
            image2((shift(1)+1):size(image2,1),-shift(2)+(1:size(image2,2)+shift(2)),:);
    elseif shift(1)<=0 && shift(2)<=0
        image2_shift(-shift(1)+(1:size(image2,1)+shift(1)),1:size(image2,2)+shift(2),:) = ...
            image2((1:size(image2,1)+shift(1)),-shift(2)+(1:size(image2,2)+shift(2)),:);
    elseif shift(1)>=0 && shift(2)>=0
        image2_shift((1:size(image2,1)-shift(1)),shift(2)+(1:size(image2,2)-shift(2)),:) = ...
             image2((shift(1)+1):size(image2,1),1:size(image2,2)-shift(2),:);
    elseif shift(1) == 0 && shift(2) == 0
        image2_shift = image2;
    end
    
    tmp        = zeros(size(image1,1),size(image1,2),3,'double');
    tmp(:,:,1) = double(sum(image1,3))/size(image1,3);
    tmp(:,:,2) = (sum(double(mat2gray(image2_shift)*255)/size(image2,3),3))*.3;
    tmp        = uint8(tmp);
    figure     % visualize the shift
    imshow(tmp);
    
    stats1 = regionprops3(Label2(:,:,:),image2_shift,'centroid','volume','MeanIntensity','VoxelValues');

    points(chan).volume        = cat(1,stats1.Volume);
    points(chan).centroid      = cat(1,stats1.Centroid);
    points(chan).MeanIntensity = cat(1,stats1.MeanIntensity);
    points(chan).Intensity     = points(chan).MeanIntensity.*points(chan).volume;

    stats2 = regionprops3(Label2(:,:,:),im_dots_scaled{chan},'MeanIntensity');
    points(chan).MeanDotIntensity = cat(1,stats2.MeanIntensity);

    Label_dots = bwlabeln(1-image2_shift);  % in the thresholded ilastik prediction map, which was rotated and shifted,
    % we now look for connected groups of pixels. Each connected group gets a unique label, which we interpret as a hotspot. 
    
    stats1 = regionprops3(Label2(:,:,:),Label_dots,'centroid','volume','MeanIntensity','VoxelValues'); 
    
    NSpots = zeros(size(stats1(:,1)));     % pre-allocating
    for z = 1 : height(stats1)
        temp = stats1.VoxelValues{z,:}; % what labels are around in the cell k?
        temp = unique(temp);          % extract the unique label values
        NSpots(z) = length(temp);     % remove contribution from 0. 
    end
    points(chan).NSpots      = NSpots;
    points(chan).NSpotsinVol = NSpots./cat(1,stats1.Volume);

    tmp   = imrotate(im_dots_scaled{chan},angle,'crop'); % rotate the original image 
    size1 = size(image1);
    size2 = size(image1);
    desiredSize  = max(size1,size2);
    padded_image = zeros(desiredSize,'uint16');
    padded_image(1:size(image2,1),1:size(image2,2),1:size(image2,3)) = tmp;
    tmp = padded_image;
    
    tmp_shift = 0*tmp; % pre-allocating
    % shifting the two images:
    if shift(1)<0 && shift(2)>0
        tmp_shift(-shift(1)+(1:(size(image2,1)+shift(1))),shift(2)+(1:size(image2,2)-shift(2)),:) = ...
            tmp(1:(size(image2,1)+shift(1)),1:size(image2,2)-shift(2),:);
    elseif shift(1)>0 && shift(2)<0
        tmp_shift((1:size(image2,1)-shift(1)),1:size(image2,2)+shift(2),:) = ...
            tmp((shift(1)+1):size(image2,1),-shift(2)+(1:size(image2,2)+shift(2)),:);
    elseif shift(1)<=0 && shift(2)<=0
        tmp_shift(-shift(1)+(1:size(image2,1)+shift(1)),1:size(image2,2)+shift(2),:) = ...
            tmp((1:size(image2,1)+shift(1)),-shift(2)+(1:size(image2,2)+shift(2)),:);
    elseif shift(1)>=0 && shift(2)>=0
        tmp_shift((1:size(image2,1)-shift(1)),shift(2)+(1:size(image2,2)-shift(2)),:) = ...
            tmp((shift(1)+1):size(image2,1),1:size(image2,2)-shift(2),:);
    elseif shift(1) == 0 && shift(2) == 0
        tmp_shift = image2;
    end
      
    Label_dots       = bwlabeln(1-image2_shift,6);  % update the dot label matrix
    stats2           = regionprops3(Label_dots,tmp_shift,'volume','MeanIntensity','centroid','VoxelValues');
    DotVolumes       = cat(1,stats2.Volume);
    DotCentre        = cat(1,stats2.Centroid);
    MeanDotIntensity = cat(1,stats2.MeanIntensity);
    
    for h = 1 : height(stats2)
    stats2.StdIntensity{h,:} = std(double(stats2.VoxelValues{h,:}));  
    end
    StdDotIntensity  = cell2mat(cat(1,stats2.StdIntensity));

    %Kmeans cluster ::::::::::::::::::::::::::::::::::::::::::::::::::::::
    % Type number of total clusters below
    rng('default')
    [assignedClass, clusterCenters]= kmeans(MeanDotIntensity,8);
   
    [sortedDistances, sortOrder] = sort(clusterCenters, 'ascend'); % Sort the clusters according to how far each cluster center is from the origin.
    newClassNumbers = zeros(size(stats1(:,1)));
    for k = 1 : size(clusterCenters, 1)
    currentClassLocations = assignedClass == k; % First find out what points have this current class & where they are
    % Now assign all of those locations to their new class.
    newClassNumber = find(k == sortOrder);	% Find index in sortOrder where this class number appears.
    newClassNumbers(currentClassLocations) = newClassNumber; % Do the relabeling right here:
    end
    group = newClassNumbers;
    
    figure();
    datacursormode on
    gscatter(MeanDotIntensity,StdDotIntensity,group); %comment this line alone out in case of low signal
    xlabel('Mean Dot Intensity') 
    ylabel('Std Intensity')
    points(chan).DotVolumes = DotVolumes;
    points(chan).MeanIntensityOfDots = MeanDotIntensity;
    points(chan).MeanIntensityOfDotsNormalized = group;
    % Change this (below) according to chosen clusters ::::::::::::::::::::
    ind = intersect(find(group >6),find(group <9));
    % :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    
    points(chan).index_filter = ind;
    figure
    imshow(max(tmp_shift,[],3),[])               
    hold on 
    plot(DotCentre(:,1),DotCentre(:,2),'b.')     % All dots shown in blue
    plot(DotCentre(ind,1),DotCentre(ind,2),'r.') % Selected dots shown in red. 

    % now use the ind to determine which of the dots are in which cells; 
    % to achieve this, we will remove all the labels of the label_dots matrix, that are
    % not member of ind: Make it a label matrix again, call Label_dots2 and
    % gets its statistics:
    
    Label_dots2 = bwlabeln(ismember(Label_dots, ind),6);
    % update the stats1 to correspond the updated Label_dots2:
    stats1 = regionprops3(Label2(:,:,:),Label_dots2,'centroid','volume','MeanIntensity','VoxelValues');
    
    NSpots = zeros(size(stats1(:,1)));     % pre-allocating
    for z = 1 : height(stats1)
        temp = stats1.VoxelValues{z,:}; % labels that are in the cell z
        temp = unique(temp);          % extract the unique label values
        NSpots(z) = length(temp);  
    end
    points(chan).NSpots_Cleared       = NSpots;
    points(chan).NSpotsinVol_Cleared  = NSpots./cat(1,stats1.Volume);
end
%% 6) Run this only after you complete running all 5 channels for that particular hyb

pointssave = sprintf('hyb%d_points.mat',hybnum);  %run this after you are done with all channels for this hyb
save(pointssave,'points');

%% the end 
%
%
%
%
%% Functions

function se = strel3D(shape, size)
    % 3D version of matlabs 2d strel function
    % Implements 'sphere' and 'cube'
    % strel3d(shap,size)
    % Copyright 2015 Idse Heemskerk and Sebastian Streichan    
    N = size;
    if strcmp(shape, 'sphere')
        se = false([2*N+1 2*N+1 2*N+1]);
        [X,Y,Z] = meshgrid(-N:N, -N:N, -N:N);
        se(X.^2 + Y.^2 + Z.^2 <= N^2) = 1;
    elseif strcmp(shape, 'cube')
        se = true([2*N+1 2*N+1 2*N+1]);
    else 
        error('strel type not recognized');
    end
end
function [shift,image1,image2] = xcorr3fft(image1,image2)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    %   xcorr2fft computes offsets between images image1 and image2 based
    %   on Phase Correlation method. image1 & image2 are assumed to have
    %   the same dimensions.   
    %   Written by: Sebastian J Streichan, EMBL, February 29, 2012
    %   Extended to 3D and bug corrected by: Stefan Gunther, EMBL, March, 20, 2012
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    F     = fftn(image1);
    Fc    = conj(fftn(image2)); 
    R     = F.*Fc;    % Phase correlation, Compute fft of image1 and conjugate fft of image2, elementwise multiply and normalise. 
    c     = ifftn(R); % Inverse fft of Phase correlation gives cross correlation map. 
    [~,i] = max(c(:));
    [I,J,~] = ind2sub(size(c),i);
    if abs(I-1)<abs(size(image1,1)-I+1)
       shiftx = -I+1;
    else
       shiftx =  size(image1,1)-I+1;
    end
    if abs(J-1)<abs(size(image1,2)-J+1)
        shifty = -J+1;
    else
        shifty = size(image1,2)-J+1;
    end
    shift=[shiftx,shifty,0]; 
end