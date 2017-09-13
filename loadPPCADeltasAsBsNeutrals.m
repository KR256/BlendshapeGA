
function loadPPCADeltasAsBsNeutrals()

load('Resources/PPCA_neutral.mat','coeff','mu','pcvar');
ppcaIdentityBSdeltas = zeros(size(coeff'));
for pc=1:size(coeff,2)
   pcCoeff = coeff(:,pc);
   sdPlus = 4  * pcCoeff * sqrt(pcvar(pc));

   ppcaIdentityBSdeltas(pc,:) = sdPlus';

   sdPlus = sdPlus' + mu;
   fileName = strcat('Blendshapes/neutral/identityPPCA/PC',num2str(pc,'%02d'),'.obj');
   writeMesh( sdPlus ,fileName );

   fprintf('PCs %i written\n',pc);
end
    
save('Blendshapes/neutral/identityPPCA/ppcaIdentityBSdeltas.mat','ppcaIdentityBSdeltas');


end