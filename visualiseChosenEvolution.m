
clear;
close all;

TEST_ID = 16;
TEST_NAME = 'ppcaIdentityBSv_Kyle.mat';

load('Resources/neutralDataset.mat','neutrals');
smileDeltasFull = neutrals;
smileDeltas = neutrals;
smileDeltas([6,8,11,16],:) = [];
load(strcat('resultingLogs/SmileLog_FACE16_',TEST_NAME));
load('Resources/PPCA_neutral.mat');

SmileLog2 = cell(size(SmileLog,1),1);

for i=1:size(SmileLog,1)
    SmileLog2{i} = vertcat(SmileLog{i,:});
end

load(strcat('resultingLogs/SmileChosenLog_FACE16_',TEST_NAME));

testSmileDelta = smileDeltasFull(TEST_ID,:);
testSmileScore = (testSmileDelta-mu) * coeff(:,1:3);



tSmileScore = zeros(size(smileDeltas,1),3);
for i=1:size(smileDeltas,1)
    tSmileScore(i,:) = (smileDeltas(i,:)-mu) * coeff(:,1:3);
end

firstGen = SmileLog2{1};
sampledSmileScore = zeros(size(firstGen,1),3);
for i=1:size(firstGen,1)
    sampledSmileScore(i,:) = firstGen(i,1:3);
end

visualiseGen(tSmileScore,testSmileScore,sampledSmileScore,1:3);

pcLims = [1,5,10,20,25];
idealFaces = [];

for gen=1:length(SmileChosenLog)
    
    if(pcLims(gen+1)==25)
        pcLims = [pcLims 25];
    end
    
    genScores = SmileChosenLog{gen};
    sampledSmileScore = genScores(2:end,1:pcLims(gen+1));
    
    testSmileDelta = smileDeltasFull(TEST_ID,:);
    testSmileScore = (testSmileDelta-mu) * coeff(:,1:pcLims(gen+1));


    tSmileScore = zeros(size(smileDeltas,1),pcLims(gen+1));
    for i=1:size(smileDeltas,1)
        tSmileScore(i,:) = (smileDeltas(i,:)-mu) * coeff(:,1:pcLims(gen+1));
    end

%     visualiseGen(tSmileScore,testSmileScore,sampledSmileScore,pcLims(gen):pcLims(gen)+2);
    visualiseGen(tSmileScore,testSmileScore,sampledSmileScore,1:3);
    
%     subplot(1,3,1);
%     scatter(genScores(1,pcLims(gen)),genScores(1,pcLims(gen)+1),200,'filled','m','p');
%     subplot(1,3,2);
%     scatter(genScores(1,pcLims(gen)),genScores(1,pcLims(gen)+2),200,'filled','m','p');
%     subplot(1,3,3);
%     scatter(genScores(1,pcLims(gen)+1),genScores(1,pcLims(gen)+2),200,'filled','m','p');
%     legend('Training Smiles','Target Smiles','Sampled Smiles','Best Chosen Smile');
    
    subplot(1,3,1);
    scatter(genScores(1,1),genScores(1,2),200,'filled','m','p');
    subplot(1,3,2);
    scatter(genScores(1,1),genScores(1,3),200,'filled','m','p');
    subplot(1,3,3);
    scatter(genScores(1,2),genScores(1,3),200,'filled','m','p');
    legend('Training Smiles','Target Smiles','Sampled Smiles','Best Chosen Smile');
    
    idealFaces = [idealFaces; genScores(1,1:3)];
    if(gen > 1)
       for t=1:gen-1
           p1 = idealFaces(t,1:2);
           p2 = idealFaces(t,[1,3]);
           p3 = idealFaces(t,2:3);
           pt1 = idealFaces(t+1,1:2);
           pt2 = idealFaces(t+1,[1,3]);
           pt3 = idealFaces(t+1,2:3);
           
           subplot(1,3,1);
           quiver(p1(1),p1(2),pt1(1)-p1(1),pt1(2)-p1(2),0,'LineWidth',2)
           subplot(1,3,2);
           quiver(p2(1),p2(2),pt2(1)-p2(1),pt2(2)-p2(2),0,'LineWidth',2)
           subplot(1,3,3);
           quiver(p3(1),p3(2),pt3(1)-p3(1),pt3(2)-p3(2),0,'LineWidth',2)
       end
       legend('Training Smiles','Target Smiles','Sampled Smiles','Best Chosen Smile','Best Smile Path');
    end
end
