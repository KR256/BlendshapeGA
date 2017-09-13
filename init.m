
clear;
close all;
clc;

TEST_IDs = [6,8,11,16];
EXPRESSION_TYPE = 'smile';
DATASET_NAME = strcat('Resources/',EXPRESSION_TYPE,'Dataset.mat');
PPCA_NumVars = 25;
numBS = 52;

%% Load all OBJ meshes and save vertices so they can be easily accessed
% loadOBJfiles('OBJs_neutrals','neutralDataset');
% loadOBJfiles('OBJs_smiles','smileDataset');

expressions = cell2mat(struct2cell(load(DATASET_NAME)));
neutrals = cell2mat(struct2cell(load('Resources/neutralDataset.mat')));
NUM_FACES = size(neutrals,1);

%% Save Expression deltas

exprDeltas = expressions - neutrals;
save(strcat('Resources/',EXPRESSION_TYPE,'DeltasFull.mat'),...
    strcat('exprDeltas'));
exprDeltas(TEST_IDs,:) = [];
save(strcat('Resources/',EXPRESSION_TYPE,'DeltasTraining.mat'),...
    strcat('exprDeltas'));

exprDeltas = cell2mat(struct2cell(load(strcat('Resources/',EXPRESSION_TYPE,'DeltasTraining.mat'))));


%% Create Blendshape Targets for Expression Deltas for Test subjects


for t=1:length(TEST_IDs)
    loadDeltasAsBs(TEST_IDs(t),neutrals,exprDeltas,EXPRESSION_TYPE)
end


%% Perform PPCA for identity
% [coeff,score,pcvar,mu,v,S] = ppca(neutrals,PPCA_NumVars);
% W = S.W;
% save('Resources/PPCA_neutral.mat','coeff','score','pcvar','mu','v','S','W');
% 

%% Create average identity obj
writeMesh(mu,'Resources/averageIdentity.obj');
neutralBSDeltas = bsxfun(@minus,neutrals,mu);
save('Blendshapes/neutral/neutralBSdeltas.mat',neutralBSDeltas);

%% Perform PPCA for Expression
% [coeff,score,pcvar,mu,v,S] = ppca(exprDeltas,PPCA_NumVars);
% W = S.W;
% save(strcat('Resources/PPCA_',EXPRESSION_TYPE,'.mat'),'coeff','score','pcvar','mu','v','S','W');

%% Create Blendshape Targets for PPCA Deltas for Test subjects

loadPPCADeltasAsBsNeutrals()

for t=1:length(TEST_IDs)
    loadPPCADeltasAsBsExpressions(TEST_IDs(t),neutrals,EXPRESSION_TYPE)
end




%% ------ SOLVED BLENDSHAPE WEIGHT METHODS ----- %%%

%% Deformation Transfer of all Training set to Ted rig
% % TAKES A LONG TIME !!!!

% Blendshapes to be ignored (as not related to smiles)
% ignoreList = [1:4 6:8 10 13:23 25 27:40 44:46 50:51 56:58 63 67:71 ... 
%     76:77 83 85:87 90:99 101:102 104:117 121:123 127:128 133:134 139];

% bsNeutrals = zeros(numBS,size(neutrals,2));
% for f=1:NUM_FACES
%     if(any(ismember(TEST_IDs,f)))      
%         neut = DeformationTransferBlendshapes(f,'test',ignoreList,EXPRESSION_TYPE);
%     else
%         neut = DeformationTransferBlendshapes(f,'test',ignoreList,EXPRESSION_TYPE);
%     end
%     bsNeutrals(f,:) = neut;
% end

% save('Blendshapes/neutral/bsNeutrals.mat','bsNeutrals');



%% Solve Each Training Smile in terms of its Blendshape Rig

% blendshapeSolve(expressions,neutrals,TEST_IDs,numBS,EXPRESSION_TYPE);

%% Perform PPCA on blendshape solves
% load(strcat('Blendshapes/',EXPRESSION_TYPE,'/weightFMatrix.mat'));
% [coeff,score,pcvar,mu,v,S] = ppca(weightFMatrix',PPCA_NumVars);
% W = S.W;
% save(strcat('Blendshapes/',EXPRESSION_TYPE,'weightMatrixPPCA.mat'),'coeff','score','pcvar','mu','v','S','W');
