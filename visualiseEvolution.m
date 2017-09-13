
clear;
close all;

TEST_ID = 6;

% load('r3ds127_neutrals.mat')
% load('r3ds127_smiles.mat')
% neutrals = neutrals([1:5,7,9:10,12:15,17:42,83:96],:);
% smiles = smiles([1:5,7,9:10,12:15,17:42,83:96],:);
% smileDeltas56 = smiles - neutrals;
% save('smileDeltas56.mat','smileDeltas56');
load('smileDeltas56.mat','smileDeltas56');
load('smileDeltas.mat','smileDeltas');

load('SmileLog_FACE6_ppcaSmiles_Kyle.mat');
load('3DMMF_smiles.mat');
Mycell = num2cell(SmileLog,[3 2]);

SmileLog = cell(Mycell);

for i=1:size(Mycell,1)
    SmileLog{i} = reshape(Mycell{i},size(Mycell{i},1)*size(Mycell{i},2),size(Mycell{i},3))';
end

% load('SmileChosenLog_FACE6_ppcaSmiles_Kyle.mat');

testSmileDelta = smileDeltas56(TEST_ID,:);
testSmileScore = (testSmileDelta-mu) * coeff(:,1:3);


tSmileScore = zeros(size(smileDeltas,1),3);
for i=1:size(smileDeltas,1)
    tSmileScore(i,:) = (smileDeltas(i,:)-mu) * coeff(:,1:3);
end

firstGen = SmileLog{1};
sampledSmileScore = zeros(size(firstGen,1),3);
for i=1:size(firstGen,1)
    sampledSmileScore(i,:) = firstGen(i,1:3);
end

subplot(1,3,1);
scatter(tSmileScore(:,1),tSmileScore(:,2),'filled');hold on;
% scatter(score(:,1),score(:,2),'filled');hold on;
scatter(testSmileScore(:,1),testSmileScore(:,2),'filled');
scatter(sampledSmileScore(:,1),sampledSmileScore(:,2),'filled');
xlabel('PC1');
ylabel('PC2');
% legend('Training Smiles','Target Smiles','Sampled Smiles');
title('PC1 v PC2 for Smiles');
subplot(1,3,2);
scatter(tSmileScore(:,1),tSmileScore(:,3),'filled');hold on;
% scatter(score(:,1),score(:,3),'filled');hold on;
scatter(testSmileScore(:,1),testSmileScore(:,3),'filled');
scatter(sampledSmileScore(:,1),sampledSmileScore(:,3),'filled');
xlabel('PC1');
ylabel('PC3');
% legend('Training Smiles','Target Smiles','Sampled Smiles');
title('PC1 v PC3 for Smiles');
subplot(1,3,3);
scatter(tSmileScore(:,2),tSmileScore(:,3),'filled');hold on;
% scatter(score(:,2),score(:,3),'filled');hold on;
scatter(testSmileScore(:,2),testSmileScore(:,3),'filled');
scatter(sampledSmileScore(:,2),sampledSmileScore(:,3),'filled');
xlabel('PC2');
ylabel('PC3');
legend('Training Smiles','Target Smiles','Sampled Smiles');
title('PC2 v PC3 for Smiles');

for gen=2:length(SmileLog)
    
    genScores = SmileLog{gen};
    sampledSmileScore = zeros(size(genScores,1),3);
    for i=1:size(genScores,1)
        sampledSmileScore(i,:) = genScores(i,1:3);
    end
    
    subplot(1,3,1);
    scatter(sampledSmileScore(:,1),sampledSmileScore(:,2));
    subplot(1,3,2);
    scatter(sampledSmileScore(:,1),sampledSmileScore(:,3));
    subplot(1,3,3);
    scatter(sampledSmileScore(:,2),sampledSmileScore(:,3));
    legend('Training Smiles','Target Smiles','Sampled Smiles');
    
end
