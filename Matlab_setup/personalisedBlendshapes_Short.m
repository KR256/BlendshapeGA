clear

TEST_IDS = [1:5 7:11];
EXPR_NAME = 'smile';
DATASET_NAME = ['Resources/' EXPR_NAME 'Dataset.mat' ];
NEUTRAL_TOLERANCE = 0.95;
EXPR_TOLERANCE = 0.95;

expressions = cell2mat(struct2cell(load(DATASET_NAME)));
neutrals = cell2mat(struct2cell(load('Resources/neutralDataset.mat')));
exprDeltas = expressions - neutrals;

for test=1:length(TEST_IDS)
    
    testID = TEST_IDS(test);
    test_neutral = neutrals(testID,:);
    test_expr = exprDeltas(testID,:);
    tempTraining_neutral = neutrals;
    tempTraining_expressions = exprDeltas;
    tempTraining_neutral(testID,:) = [];
    tempTraining_expressions(testID,:) = [];
    mean_smile = mean(tempTraining_expressions);
    tempTraining_expressionsC = bsxfun(@minus,tempTraining_expressions,mean_smile);
    
    [coeff_neutral,score_neutral,latent_neutral] = pca(tempTraining_neutral);
    [coeff_expr,score_expr,latent_expr] = pca(tempTraining_expressionsC);
    
    neutral_PCLIM = find((cumsum(latent_neutral)./sum(latent_neutral)) ...
        >= NEUTRAL_TOLERANCE,1);
    
    expr_PCLIM = find((cumsum(latent_expr)./sum(latent_expr)) ...
        >= EXPR_TOLERANCE,1);
    
    Regression = linsolve(score_neutral(:,1:neutral_PCLIM),...
        score_expr(:,1:expr_PCLIM));
    
    Regression2 = linsolve(score_neutral,score_expr);
    Regression3 = linsolve(score_neutral(:,1:neutral_PCLIM),score_expr);
    Regression4 = linsolve(score_neutral,score_expr(:,1:expr_PCLIM) );
    
    test_neutralScore = (test_neutral - mean(tempTraining_neutral)) * coeff_neutral;
    test_neutralScore2 = (test_neutral - mean(tempTraining_neutral)) * coeff_neutral(:,1:neutral_PCLIM);
    
    exprScore1 = test_neutralScore(:,1:neutral_PCLIM) * Regression;
    exprScore2 = test_neutralScore * Regression2(:,1:expr_PCLIM);
    exprScore3 = test_neutralScore(:,1:neutral_PCLIM) * Regression3;
    exprScore5 = test_neutralScore * Regression2;

    exprVector1_1 = exprScore1 * coeff_expr(:,1:expr_PCLIM)'; 
    exprVector2_1 = exprScore2 * coeff_expr(:,1:expr_PCLIM)';
    exprVector3_1 = exprScore3 * coeff_expr';
    exprVector5_1 = exprScore5 * coeff_expr';
    
    predFace1 = test_neutral + exprVector1_1 + mean_smile;
    predFace2 = test_neutral + exprVector2_1 + mean_smile;
    predFace3 = test_neutral + exprVector3_1 + mean_smile;
    predFace5 = test_neutral + exprVector5_1 + mean_smile;
    
    avgFace = test_neutral + mean_smile;
    
%     writeMesh(avgFace,strcat('avgSmile',num2str(testID),'.obj'));
    writeMesh(predFace1,strcat('TEST1FACE',num2str(testID),'.obj'));
    writeMesh(predFace2,strcat('TEST2FACE',num2str(testID),'.obj'));
    writeMesh(predFace3,strcat('TEST3FACE',num2str(testID),'.obj'));
    writeMesh(predFace5,strcat('TEST4FACE',num2str(testID),'.obj'));
end