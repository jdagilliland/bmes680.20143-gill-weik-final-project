%% Close existing windows
close all
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

%% Stack data
gene_table = genes.Data;
% [n_gene, n_genetypes] = size(gene);
expvalues = gsedata.Data;
gene_id = rownames(expvalues);
n_gene = numel(gene_id);

%% Stack data from 3 brain regions.
stacked_data = zeros(3*n_gene, numel(subj_id_unique));
stacked_genes = cell(3*n_gene, 1);
stacked_tissue = cell(3*n_gene, 1);
tissue_type_unique = {'CR', 'PFC', 'VC'};
for j = 1:3
     stacked_tissue(((j-1)*n_gene+1):j*n_gene) = tissue_type_unique(j);
     stacked_genes(((j-1)*n_gene+1):j*n_gene) = gene_id;
end
for i = 1:numel(subj_id_unique)
    datarows = expvalues(:, strcmp(subj_id, subj_id_unique(i)));
    for j = 1:3
        stacked_data(((j-1)*n_gene+1):j*n_gene, i) = ...
            datarows(:, j);
    end
end

%% filter data
[mask,stacked_data] = genevarfilter(stacked_data,'Percentile',30);
stacked_genes = stacked_genes(mask, :);
stacked_tissue = stacked_tissue(mask);
%[mask,stacked_data] = geneentropyfilter( ...
%   stacked_data,'Percentile',30);
%stacked_genes = stacked_genes(mask, :);
%stacked_tissue = stacked_tissue(mask);
for i = 1:length(subj_id_unique)
    unique_id(i) = str2num(cell2mat(subj_id_unique(i)));
end
grps=disease_state(unique_id);
slct_normal=strcmp(grps,'normal');
slct_alz=strcmp(grps,'Alzheimer''s disease');
[h,p] = ttest2(stacked_data(:, slct_normal)', stacked_data(:, slct_alz)');
i_sigp = p<1e-6;
stacked_data = stacked_data(i_sigp,:);
stacked_genes = stacked_genes(i_sigp, :);
stacked_tissue = stacked_tissue(i_sigp);
% i_nonan = sum(isnan(Filtdata), 1) == 0;
% stacked_data = stacked_data(i_nonan,:);
% stacked_genes = stacked_genes(i_nonan, :);

%% Check distribution of significantly regulated genes across tissues.
n_sig_CR = nnz(strcmp(stacked_tissue, 'CR'));
n_sig_PFC = nnz(strcmp(stacked_tissue, 'PFC'));
n_sig_VC = nnz(strcmp(stacked_tissue, 'VC'));

%% k-means clustering
opts = statset('Display','final');
[coeff, score, latent,~,explained,~]=pca(stacked_data');
[idx,ctrs] = kmeans(stacked_data',2,'Options',opts);
figure
plot(score(idx==1,1),score(idx==1,2),'r.','MarkerSize',12)
hold on
plot(score(idx==2,1),score(idx==2,2),'b.','MarkerSize',12)
legend('Cluster 1','Cluster 2',...
       'Location','NW')
hold off
title('K-means Clustering of Samples');
%% Cluster by disease and age
disease_state=strcmp(grps,'normal');
%% Assess quality of K-means clustering.
[R_alz_kmean, P_alz_kmean] = corrcoef( ...
   disease_state(~isnan(idx)), idx(~isnan(idx)))
confmat = confusionmat(double(~disease_state(~isnan(idx))), ...
   idx(~isnan(idx))-1)

%% Plot PCA colored by Alzheimer's.
figure
plot(score(disease_state==1,1),score(disease_state==1,2),'r.','MarkerSize',12)
hold on
plot(score(disease_state==0,1),score(disease_state==0,2),'b.','MarkerSize',12)
legend('Normal','Alzheimer''s Disease',...
       'Location','NW')
   hold off
title('Cluster by Disease State and PCA score')
age_unique=age(unique_id);
[decsage,decageindx]=sort(age_unique);
diseaseage=disease_state(decageindx);
colorVec=winter(5);
group=repmat([1 2 3 4 5], 46, 1);
scoreage=score(decageindx,:);
figure
for i=1:230
    if diseaseage(1,i)==1
        norm=plot(scoreage(i,1),scoreage(i,2),'*','Color', ...
            colorVec(group(i),:),'MarkerSize',12);
    else
        alz=plot(scoreage(i,1),scoreage(i,2),'o','Color', ...
            colorVec(group(i),:),'MarkerSize',12);
    end
    hold on
end
legend([norm,alz],'Normal','Alzheimer''s disease',...
       'Location','NW')
title('Correlation between Age and Disease')
hold off;
%% Seperate data into diseased and control groups
disease_control_tree = fitctree(stacked_data',grps','MaxCat',2);
resuberror = resubLoss(disease_control_tree)
view(disease_control_tree,'Mode','graph')
% Extract what the predictors are.
imp = predictorImportance(disease_control_tree);
predictors = imp > 0;
[maximps, imaximps] = sort(imp, 'descend');
predictorgenes = extract_gene_info(gene_table, stacked_genes(imaximps(1:6)));

%% Correlate age and genes
% Assumes that missing data are normal.
Alz_score = strcmp(grps, 'Alzheimer''s disease');
Norm_score = strcmp(grps, 'normal');
Norm_age=age_unique(Norm_score);
Alz_age=age_unique(Alz_score);
[n_feature,~]=size(stacked_genes);
R_alz = zeros(n_feature, 1);
R_norm = zeros(n_feature, 1);
P_alz = zeros(n_feature, 1);
P_norm = zeros(n_feature, 1);
for i = 1:n_feature
    slct_samps = ~isnan(age_unique) & ~isnan(stacked_data(i, :));
    slct_alz = slct_samps & Alz_score;
    slct_norm = slct_samps & Norm_score;
    [R_alz_i, P_alz_i] = corrcoef(age_unique(slct_alz), ...
       stacked_data(i, slct_alz));
    [R_norm_i, P_norm_i] = corrcoef(age_unique(slct_norm), ...
       stacked_data(i, slct_norm));
    R_alz(i) = R_alz_i(2);
    R_norm(i) = R_norm_i(2);
    P_alz(i) = P_alz_i(2);
    P_norm(i) = P_norm_i(2);
end
Alz_sigp = P_alz<1e-6;
Norm_sigp = P_norm<1e-6;
% Normrelated=find(Age_sigp & Norm_sigp);
% Alzrelated=find(Age_sigp & Alz_sigp);
genediff1=Alz_sigp & ~Norm_sigp;
genediff2=Norm_sigp & ~Alz_sigp;
Age_Alz_genes_table = extract_gene_info(gene_table, stacked_genes(genediff1));
Age_Norm_genes_table = extract_gene_info(gene_table, stacked_genes(genediff2));
Age_Alz_genes=Age_Alz_genes_table(:, 4);
Age_Norm_genes=Age_Norm_genes_table(:, 4);
Age_Alz_tissues = stacked_tissue(genediff1);
Age_Norm_tissues = stacked_tissue(genediff2);
% Count the genes that belong to this set of ALZ but not age related by
% brain region.
n_age_alz_pfc = nnz(strcmp(Age_Alz_tissues, 'PFC'))
n_age_alz_cr = nnz(strcmp(Age_Alz_tissues, 'CR'))
n_age_alz_vc = nnz(strcmp(Age_Alz_tissues, 'VC'))
n_age_norm_pfc = nnz(strcmp(Age_Norm_tissues, 'PFC'))
n_age_norm_cr = nnz(strcmp(Age_Norm_tissues, 'CR'))
n_age_norm_vc = nnz(strcmp(Age_Norm_tissues, 'VC'))
%% Clustergram
groups=double(disease_state);
cm=struct('ColumnLabels',{1,0},'Color',{'b','m'});
cgo_all = clustergram(stacked_data(1:30,1:30),'ColumnLabels',groups(1,1:30))
set(cgo_all,'ColumnLabelsColor',cm)
