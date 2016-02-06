%==========================================================================
% Script of feature selection with Cox regression & permutation logrank test
%--------------------------------------------------------------------------
% [Version] V1.0; by Percy Zhao. [Date]: 2015-5-16
%--------------------------------------------------------------------------
% [Instruction]: 
%  ...This method applies Cox regression by enumerating features or feature
%  ...groups one by one. Then, the median of risk score is presented to 
%  ...establish high-risk and low-risk group. PERMUTATION Logrank test is
%  ...used to test whether these two groups are significantly different.
% [Input]:
%  ...expression profiles.
% [Output]:
%  ...significant features for apparting the low-risk and high-risk group.
%--------------------------------------------------------------------------
% Function 1 <- S3_parallel_kernel_permutation.m
%==========================================================================
%% 
%==================%
% I.Initialization %
%==================%
clear;
clc;
%%%---I.0-paramters---%%%
%---values---%
core = 12;
iter_round = 1e2;   % for permutation
%---switches---%
debug_open = 0;
disp_round = 1;
appoint_or_not = 0;
pathway_or_not = 'no';   % yes/no
%%%---I.1-Select treated object---%%%
%---cancer type---%
cancer = 'GBM';
%---feature type---%
string_source = 'miRNA';
%---pathway---%
if strcmp(pathway_or_not,'yes')==1
    pathway = '\hedgehog';
end
%---label type---%
label_type = 'survival';
%---offset---%
row_head = 1;       % Head of table for genes
col_head = 1;       % Head of table for samples
offset = 2;
%%%---I.2-I/O input set---%%%
%---working directory---%
folder_source = cd;
%---input data set---%
if strcmp(pathway_or_not,'yes')==1
    %filename = ['input' pathway '\S2_pathway_selection.xlsx'];
    filename = ['input' pathway '\S2_pathway_selection_without_30days.xlsx'];
else
    %filename = ['input' '\S2_combination.xlsx'];
    filename = '\S2_combination.xlsx';
end
[Data,Head,~]=xlsread(filename);
N = size(Head,2);
Label = zeros(1,N-col_head);
for col = col_head+1:N
    if strcmp(Head{row_head+1,col},'Dead')
        Label(1,col-col_head) = 1; % 1-dead; 0-censored
    end
end
Head = Head(:,1);
Data = [Label;Data];
Time = Data(2,:);
[n,N] = size(Data);   % n-gene size; N-sample size
%---iter_num---%
iter_num = input('Please input the size of features for enumeration(>=1): ');
while iter_num < 1 && isnumeric(iter_num) == 1 && ceil(iter_num)~=iter_num
    iter_num = input('Wrong...please input the size of features for enumeration(>=1): ');
end
%---loading broad screening selected gene groups---%
if appoint_or_not == 0
    features = nchoosek(1:n-offset,iter_num);
else
   file_features = input('Please input the table containing the feature you select: ','s');
   [features,~,~]=xlsread(file_features);
end
first_num = size(features,1);
%%%---I.3-I/O output set---%%%
%---output path---%
if appoint_or_not == 0
    if strcmp(pathway_or_not,'yes')==1
        o_folder = [folder_source '\output\'  cancer '\' string_source pathway '\' label_type '_with_' string_source '_iter_' num2str(iter_num) '\permutation_all'];
    else
        o_folder = [folder_source '\output\'  cancer '\' string_source '\' label_type '_with_' string_source '_iter_' num2str(iter_num) '\permutation_all'];
    end
else
    if strcmp(pathway_or_not,'yes')==1
        o_folder = [folder_source '\output\'  cancer '\' string_source pathway '\' label_type '_with_' string_source '_iter_' num2str(iter_num) '\permutation_select'];
    else
        o_folder = [folder_source '\output\'  cancer '\' string_source '\' label_type '_with_' string_source '_iter_' num2str(iter_num) '\permutation_select'];
    end
end    
mkdir(o_folder);
%======I.end=======%
save('Step_1.mat');
disp('Step 1-initialization is over.');
%%
%===========================%
% II.Parametric enumeration %
%===========================%
%---X.open matlabpool---%
if 0 == matlabpool('size')
    matlabpool local
end
%for k = 1 : core
parfor k = 1 : core
    S3_parallel_kernel_permutation(k);   %---<-Function 1---%
end
%---X.close matlabpool---%
matlabpool close
%=========II.end============%
disp('Step 2-Parametric enumeration is over.');
%%
%===========%
% III.merge %
%===========%
Acc_P_save = zeros(first_num,2+iter_num*4); % features,z(logrank),p(logrank),b(cox),z(cox),p(cox)
num = zeros(1,core);
index_sum = 1;
for i = 1 : core
    filename = [o_folder '\P_core-' num2str(i) '.mat'];
    load(filename);
    num(i) = min(i*ceil(first_num/core),first_num)-(i-1)*ceil(first_num/core);
    Acc_P_save(index_sum:index_sum+num(i)-1,:) = P_save;
    index_sum = index_sum + num(i);
end
%==III.end==%
disp('Step 3-Merge is over.');
%%
%==========%
% V.Output %
%==========%
%---V.1-output storage pool---%
Table_head = cell(1,5*iter_num+2);
for i = 1 : iter_num
    Table_head(1,i) = {[string_source ' probe']};
    Table_head(1,i+iter_num) ={'Row number'};
    Table_head(1,i+2*iter_num+2) = {'Coef(Cox)'};
    Table_head(1,i+3*iter_num+2) = {'Z(Cox)'};
    Table_head(1,i+4*iter_num+2) = {'P(Cox)'};
end
Table_head(1,2*iter_num+1) = {'Z(log-rank)'};
Table_head(1,2*iter_num+2) = {'P(log-rank)'};
%---V.2-import gene symbol---%
Table = cell(first_num+1,5*iter_num+2);
Table(1,:) = Table_head;
for i = 1 : first_num
    for j = 1 : iter_num
        Table(i+row_head,j)=Head(row_head+offset+Acc_P_save(i,j));
    end
end
%---V.3-import data---%
Table(2:end,iter_num+1:end) = num2cell(Acc_P_save);
%---V.4-output---%
%filename = [o_folder '\' label_type '_with_' string_source '_iter_' num2str(iter_num) '_permutation_' num2str(iter_round) '.xlsx'];
filename = [o_folder '\' label_type '_with_' string_source '_iter_' num2str(iter_num) '_' num2str(iter_round) '.xlsx'];
xlswrite(filename,Table);
%==V.end===%
disp('Step 4-Output is over.');