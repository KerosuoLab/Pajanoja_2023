 
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%   Maintenance of pluripotency in entire post-gastrulation ectoderm enables neural crest formation
%
%   Ceren Pajanoja, 2022 

%   The following script includes the code used to execute SECTION-3 of scMST pipeline

%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%

%% IMPORTANT = in order to run "centroids" part of this code, you need to change the root in matlab

% choose new sub groups from the big heatmap called "heatm"
figure, 
indices = [5397 5408 5389 5336];


cmp = jet(length(indices)); % why does this need to be in the same length as indices?
% to assign colors of choice: 
cmp(1,:)  = [0.9,0,0];   %red
cmp(2,:)  = [0,0.1,1];   %blue 
cmp(3,:)  = [0,0.7,0.1]; %green 
cmp(4,:)  = [0.8,1,0.1]; %yellow
cmp(5,:)  = [1,0.7,1];   %pink 
cmp(6,:)  = [0.9,0.6,0]; %orange  
cmp(7,:)  = [0.4,1,1];   %Light blue
cmp(8,:)  = [0.5,0,0.6]; %purple

%--------------------------------------------------------
% experiment that you'd like to inspect
% in rest of the script ONLY change "xp_ref_num"
%--------------------------------------------------------


xp_ref_num =15 ;        %type here number of the order of your sample is in (example: sample1 is 1, sample2 is 2 etc)
nucleusImage = fullfile(root_path,folders{xp_ref_num},'RGB.tif'); % this is membrane and dapi staining
NucleusStackSize = length(imfinfo(nucleusImage));
imNuc = cell(1,2);
for k = 1 : NucleusStackSize
    temp = imread(nucleusImage,k);
    imNuc{1}(:,:,k) = temp(:,:,1);
    imNuc{2}(:,:,k) = temp(:,:,2);
end
image1 = imNuc{1};
imshow(image1(:,:,40),[])

counter = 1;   %this is for the color 
sub_cluster_BIG = struct('cells_of_interest',[],'cells_of_interest_small',[]);
CellIdentities = cat(1,xp.Big_ID); 
bds = [0;cumsum(cat(1,xp.M_length))];

for num = indices % can range between 1, and nCells-1;
bla = clusterGroup(heatm, num, 'col'); % this is where you write which heatmap you take the indice from
group_of_interest= bla;

Col_Labels = group_of_interest.ColumnLabels; % here we have extracted some column labels from the clustergram; 
Double_Labels = cell(0);
for i = 1 : length(Col_Labels)
    Double_Labels{i} = str2num(Col_Labels{i});
end

goi = cell2mat(Double_Labels);
% cell identities: fill the new cell ids into struct
cells_of_Interest_Big = CellIdentities(goi);

sub_cluster_BIG(counter).cells_of_interest = cells_of_Interest_Big;
% sort through cells_of_interest, identify corresponding xp_num, and
% identify in that. 
for xp_num = 1 : length(xp) % Dont change this value!!
    % identify the cells that come from expertiment xp_num
    in = intersect(find(cells_of_Interest_Big>bds(xp_num)),find(cells_of_Interest_Big<bds(xp_num+1)));
    % find the corresponding identities in that experiment 
    template = xp(xp_num).Big_ID; 
    pattern = cells_of_Interest_Big(in);
    D = pdist2(template,pattern);
    in_xp_num = find(min(D,[],2)==0);
    small_id = xp(xp_num).CellIdentities(in_xp_num);
    sub_cluster_BIG(counter).cells_of_interest_small{xp_num} = small_id; 
 end

hold on

plot(xp(xp_ref_num).points_hyb1.points(1).centroid(sub_cluster_BIG(counter).cells_of_interest_small{xp_ref_num},1),...
    xp(xp_ref_num).points_hyb1.points(1).centroid(sub_cluster_BIG(counter).cells_of_interest_small{xp_ref_num},2),'*','color',cmp(counter,:),'LineWidth',10);
counter = counter+1;
end
% here the end result is image of membrain staining with centroids for
% selected subgroups

%% Getting the pattern (no need to modify below, just run through it)
sub_cluster_BIG_pruned = sub_cluster_BIG;
% prune the sub_cluster
for i = 1 : length(sub_cluster_BIG)
    % check in all others, if fully contained
    template = sub_cluster_BIG(i).cells_of_interest;% template
    for j = 1 : length(sub_cluster_BIG)
        pattern = sub_cluster_BIG(j).cells_of_interest;% find pattern in template
        if isempty(setdiff(pattern,template)) &&(i~=j)
            sub_cluster_BIG_pruned(i).cells_of_interest = setdiff(template,pattern); % remove the pattern from template;
        end
    end

    cells_of_Interest_Big = sub_cluster_BIG_pruned(i).cells_of_interest; 
    for xp_num = 1 : length(xp)
        % identify the cells that come from expertiment xp_num
        in = intersect(find(cells_of_Interest_Big>bds(xp_num)),find(cells_of_Interest_Big<bds(xp_num+1)));
        % find the corresponding identities in that experiment 
        template = xp(xp_num).Big_ID; 
        pattern = cells_of_Interest_Big(in);
        D = pdist2(template,pattern);
        in_xp_num = find(min(D,[],2)==0);
        small_id = xp(xp_num).CellIdentities(in_xp_num);
        sub_cluster_BIG_pruned(i).cells_of_interest_small{xp_num} = small_id; 
    end
end

%% Crop out only the cells of interest from the whole embryo below: 

Label_temp = load(fullfile(root_path,folders{xp_ref_num},'Label2_WatershedMatrix.mat')); 
Label2 = Label_temp.Label2;
% Label_sub is a subset of the label matrix, subset is chosen according to cells_of_interest

Label_sub = 0*Label2;
for c = 1 : length(sub_cluster_BIG_pruned)
    cells_of_Interest = sub_cluster_BIG_pruned(c).cells_of_interest_small{xp_ref_num};

    for i = 1 : length(cells_of_Interest) 
        Label_sub(Label2==cells_of_Interest(i)) = c;
    end
end

%% Here we save it as tiff file to be used later in Fiji/ImageJ

for k = 1 : size(Label_sub,3)
    temp  = zeros(size(Label_sub,1),size(Label_sub,2),3,'uint8');
    for c = 1 : length(sub_cluster_BIG_pruned)
        tmp = Label_sub(:,:,k);
        tmp(tmp~=c) = 0;
        tmp = tmp>0;
        temp(:,:,1) = uint8(cmp(c,1).*255*double(tmp))+temp(:,:,1);
        temp(:,:,2) = uint8(cmp(c,2).*255*double(tmp))+temp(:,:,2);
        temp(:,:,3) = uint8(cmp(c,3).*255*double(tmp))+temp(:,:,3);
        end
    imwrite(temp,['WriteFileNameHere_',num2str(xp_ref_num),'.tif'],'tiff','Compression','none','WriteMode','append');
end

%--------------------------------------------------------------------------
% End of Section 3
% Now move to Fiji/ImageJ to create pseudocolored images
%--------------------------------------------------------------------------

