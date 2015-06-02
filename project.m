%% Download the data if necessary.
if ~exist('GSE44772.txt','file')
    % Only download data from net if file not already in directory.
    disp('Retrieving GSE data...')
    gsedata=getgeodata('GSE44772', 'ToFile','GSE44772.txt');
    disp('done!')
elseif ~exist('gsedata', 'var')
    % Only load from file if not already in workspace.
    disp('Loading GSE data...')
    gsedata=geoseriesread('GSE44772.txt');
    disp('done!')
else
    disp('Skipping download, reload of gene expression data')
end

if ~exist('genes.txt','file')
    % Only download data from net if file not already in directory.
    disp('Retrieving platform data...')
    genes=getgeodata(gsedata.Header.Series.platform_id,'ToFile','genes.txt');
    disp('done!')
elseif ~exist('genes', 'var')
    % Only load from file if not already in workspace.
    disp('Loading platform data...')
    genes=geosoftread('genes.txt');
    disp('done!')
else
    disp('Skipping download, reload of platform data')
end

%% filter data
expvalues = gsedata.Data;
gene = genes.Data;
[mask,Fdata] = genelowvalfilter(expvalues,'absval',log2(2));
gene = gene(mask, :);
[Fmask,fildata] = geneentropyfilter(Fdata,'Percentile',30);
gene = gene(Fmask, :);
[h,p] = ttest(fildata');
i_sigp = p<0.001;
Filtdata = fildata(i_sigp,:);
gene = gene(i_sigp, :);

%% Extract metadata

% Extract subject ids and tissue type.
reg_title = '(\d*)_(.*)';
subj_id = extract_meta(gsedata.Header.Samples.title, reg_title, 1);
tissue_type = extract_meta(gsedata.Header.Samples.title, reg_title, 2);

% Extract disease state.
reg_disease = 'disease status: (.*)';
disease_state = extract_meta(...
    gsedata.Header.Samples.characteristics_ch2(4,:), reg_disease, 1);
% Extract age.
reg_age = 'age: (\d*)';
age = str2double(extract_meta(...
    gsedata.Header.Samples.characteristics_ch2(1,:), ...
    reg_age, 1));
% Extract gender.
reg_gender = 'gender: ([MF])';
gender = extract_meta( gsedata.Header.Samples.characteristics_ch2(3,:), ...
    reg_gender, 1);

%% Find unique subject ids.
subj_id_unique = unique(subj_id);
[n_gene, ~] = size(Filtdata);

%% Stack data from 3 brain regions.
stacked_data = zeros(3*n_gene, numel(subj_id_unique));
for i = 1:numel(subj_id_unique)
    datarows = Filtdata(:, strcmp(subj_id, subj_id_unique(i)));
    for j = 1:3
        stacked_data(((j-1)*n_gene+1):j*n_gene, i) = ...
            datarows(:, j);
    end
end

%% k-means clustering
opts = statset('Display','final');
[idx,ctrs] = kmeans(stacked_data,2,'Options',opts);
silhouette(stacked_data,idx)
figure
plot(stacked_data(idx==1,1),stacked_data(idx==1,2),'r.','MarkerSize',12)
hold on
plot(stacked_data(idx==2,1),stacked_data(idx==2,2),'b.','MarkerSize',12)
plot(ctrs(:,1),ctrs(:,2),'kx',...
     'MarkerSize',12,'LineWidth',2)
plot(ctrs(:,1),ctrs(:,2),'ko',...
     'MarkerSize',12,'LineWidth',2)
legend('Cluster 1','Cluster 2','Centroids',...
       'Location','NW')
hold off
title('K-means Clustering of Samples');

%% Seperate data into diseased and control groups
for i = 1:length(subj_id_unique)
   unique_id(i) = str2num(cell2mat(subj_id_unique(i)));
end
grps=disease_state(unique_id);
disease_control_tree = fitctree(stacked_data',grps');
resuberror = resubLoss(disease_control_tree)
figure
view(disease_control_tree,'Mode','graph')
