
function loadDeltasAsBs(TEST_FACE,neutrals,deltas,exprName)

neutralFace = neutrals(TEST_FACE,:);

for pc=1:size(deltas,1)
       
   exprOut = deltas(pc,:) + neutralFace;
   fileName = strcat('Blendshapes/',exprName,'/',exprName,'Pop/TEST',num2str(TEST_FACE),'/FACE',num2str(pc,'%02d'),'.obj');
   writeMesh( exprOut ,fileName );
      
   fprintf('PCs %i written\n',pc);
end

end