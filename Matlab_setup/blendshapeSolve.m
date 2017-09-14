% clear;

function blendshapeSolve(expressions,neutrals,testIDs,numBs,EXPRESSION_TYPE)

weightFMatrix = zeros(numBs,size(neutrals,1));
load('Blendshapes/neutral/bsNeutrals.mat');

for faceID=1:size(expressions,1)
    
    expr_vts = expressions(faceID,:);    
    neutral_vts = neutrals(faceID,:);
    num_vertices = size(neutrals,2);
    
    if(any(ismember(testIDs,faceID))) 
        load(strcat('Blendshapes/',EXPRESSION_TYPE,'/deformedTestTargets/FACE_',num2str(faceID),...
            '/bsVerts52_',num2str(faceID),'.mat'));
    else
        load(strcat('Blendshapes/',EXPRESSION_TYPE,'/deformedTrainingTargets/bsVerts52_',...
            num2str(faceID),'.mat'));
    end
    bSmatrix = bsVertsMat';
    
    neutralBs_vts = bsNeutrals(faceID,:);
    
    %%F(w) - b0 = B_hat * w .... ( Target - neutral) = Sum weighted Blendshapes
    B_hat = bsxfun(@minus, bSmatrix, neutralBs_vts);
    FwMinb0 = expr_vts - neutral_vts;
    
    lb = 0*ones(3*num_vertices,1);
    ub = 1*ones(3*num_vertices,1);
    [w3,resnorm2] = lsqlin(B_hat, FwMinb0,[], [], [],[],lb,ub);
    w3(w3<0.05) = 0;
    bb3 = bsxfun(@plus, B_hat*w3, neutral_vts');
    writeMesh(bb3',strcat('outSolves/',EXPRESSION_TYPE,'/weightFReconstructed_',num2str(faceID),'.obj'));
    
    weightFMatrix(:,faceID) = w3;
    
end

save(strcat('Blendshapes/',EXPRESSION_TYPE,'/weightFMatrix.mat'),'weightFMatrix');

