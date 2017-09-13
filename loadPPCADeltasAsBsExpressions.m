
function loadPPCADeltasAsBsExpressions(TEST_ID,neutrals,exprName)

load(strcat('Resources/PPCA_',exprName,'.mat'),'coeff','mu','pcvar');
neutralFace = neutrals(TEST_ID,:);
ppcaExprBSdeltas = zeros(size(coeff'));
for pc=1:size(coeff,2)
   pcCoeff = coeff(:,pc);
   sdPlus = 4  * pcCoeff * sqrt(pcvar(pc));

   ppcaExprBSdeltas(pc,:) = sdPlus';

   sdPlus = sdPlus' + mu + neutralFace;
   fileName = strcat('Blendshapes/',exprName,'/',exprName,'PPCA/TEST',num2str(TEST_ID),'/PC',num2str(pc,'%02d'),'.obj');
   writeMesh( sdPlus ,fileName );

   fprintf('PCs %i written\n',pc);
end
    
save(strcat('Blendshapes/',exprName,'/',exprName,'PPCA/ppca',exprName,'BSdeltas.mat'),'ppcaExprBSdeltas');


end

