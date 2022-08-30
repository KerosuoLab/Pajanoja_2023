 
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%   Maintenance of pluripotency in entire post-gastrulation ectoderm enables neural crest formation
%
%   Ceren Pajanoja, 2022 

%   The following script includes the code used to execute SECTION-2 of scMST pipeline
%   After the initial analysis is done, here we pool everything 
%   and apply Z-scoring 
%

%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%

%% Load All Sample Files
root_path = '\path\to\data\folder\';
folders   = {'e1-sample-1','e1-sample-2','e1-sample-3',...%1-2-3
             'e2-sample-1','e2-sample-1','e2-sample-1',...%4-5-6
             'e3-sample-1','e3-sample-2','e3-sample-3',...%7-8-9
             'e4-sample-1','e4-sample-2','e4-sample-3',...%10-11-12
             'e5-sample-1','e5-sample-2','e5-sample-3',...%13-14-15
             'e6-sample-1','e6-sample-2','e6-sample-3'};  %16-17-18
                  
xp = struct('points_hyb1',[],'points_hyb2',[],'points_hyb3',[],...
            'points_hyb4',[],'points_hyb5',[],'points_hyb6',[],...
            'points_hyb7',[],'M',[],'VolFilt',[],'CellIdentities',[],...
            'Centroid',[],'M_length',[],'Big_ID',[]); 
%% Pool All Samples 
counter = 0;
for xp_num = 1:length(folders)
    xp(xp_num).points_hyb1 = load(fullfile(root_path,folders{xp_num},'hyb1_points.mat'));
    xp(xp_num).points_hyb2 = load(fullfile(root_path,folders{xp_num},'hyb2_points.mat'));
    xp(xp_num).points_hyb3 = load(fullfile(root_path,folders{xp_num},'hyb3_points.mat'));
    xp(xp_num).points_hyb4 = load(fullfile(root_path,folders{xp_num},'hyb4_points.mat'));
    xp(xp_num).points_hyb5 = load(fullfile(root_path,folders{xp_num},'hyb5_points.mat'));
    xp(xp_num).points_hyb6 = load(fullfile(root_path,folders{xp_num},'hyb6_points.mat'));
    xp(xp_num).points_hyb7 = load(fullfile(root_path,folders{xp_num},'hyb7_points.mat'));
    %First M is filled with vol,NspotscCleared for hyb1-7
    xp(xp_num).M = cat(2,xp(xp_num).points_hyb1.points(1).volume,xp(xp_num).points_hyb1.points.NSpotsinVol_Cleared,...
           xp(xp_num).points_hyb2.points.NSpotsinVol_Cleared,xp(xp_num).points_hyb3.points.NSpotsinVol_Cleared,...
           xp(xp_num).points_hyb4.points.NSpotsinVol_Cleared,xp(xp_num).points_hyb5.points.NSpotsinVol_Cleared,...
           xp(xp_num).points_hyb6.points.NSpotsinVol_Cleared,xp(xp_num).points_hyb7.points.NSpotsinVol_Cleared);
    xp(xp_num).Centroid = xp(xp_num).points_hyb1.points(1).centroid;
    %volume filt is applied to first column of M matrix= which is the volume= also created new field called VolFilt     
    xp(xp_num).VolFilt = ~(xp(xp_num).M(:,1)<10000 | xp(xp_num).M(:,1)>200000); %volumetric filter between these values
    xp(xp_num).M(~xp(xp_num).VolFilt,:) = []; 
    xp(xp_num).M(isnan(xp(xp_num).M))=0;
    %xp(xp_num).M = zscore(xp(xp_num).M);
    %in = ~isnan(M(:,2));
    %M(~in,:)=[];
    xp(xp_num).CellIdentities = find(xp(xp_num).VolFilt);
    xp(xp_num).Mt = xp(xp_num).M.';
    % xlswrite('points_Volfiltered.xlsx',xp(xp_num).M);
    xp(xp_num).M_length = size(xp(xp_num).M,1);
    xp(xp_num).Big_ID = counter+(1:xp(xp_num).M_length)';
    counter = counter+xp(xp_num).M_length;
end

%% Zscore embryos (3 field of views per embryo)
% Example: 18 images = 6 embryos (3 fieldofviews each)

j=1;
for k = 1:3:length(xp)
    embryos(j).embs = zscore(vertcat(xp(k).M,xp(k+1).M,xp(k+2).M));
    j = j+1;
end
Big_M = zscore(cat(1,embryos(1:6).embs));
M_cols = Big_M;

%% Put genes in order in genesPerm variable and create Heatmap

genes = {'Foxd3','Bcl2','Sox9','Sip1','Dlx5','MycC','Msi1','Msx2','Klf4',...
    'Tfap2A','Sox10','Nanog','Pax7','Sox2','Snai2','MycN','Pax6','Ets1','PouV',...
    'Msx1','Axud1','Sox8','Tfap2B','Six1','Ccnd1','Runx2','Lin28','Mitf','Col2a1',...
    'Eya2','Krt14','Fabp7','Nestin','Krt19','Gata3'};
PermGenes = [18 21 1 6 11 15 3 23 8 13 20 10 30 24 35 17 22 5 4 16 33 7 14 19 9 12 27 25 2 28 32 29 34 31 26];
            
genesPerm = genes(:,PermGenes);
MPerm = M_cols(:,PermGenes);
MtPerm = MPerm.';    
genenum = size(genesPerm(:));
heatm = clustergram(MtPerm(1:genenum,:),'RowLabels',genesPerm,...
    'RowPDist','cosine','ColumnPDist','cosine','linkage','average',...
    'DisplayRange',3,'Colormap',redbluecmap,'Cluster',3);

%%

%--------------------------------------------------------------------------
% End of Section 2
% Now move to Fiji/ImageJ to create pseudocolored images
%--------------------------------------------------------------------------

