function visualiseGen(tSmileScore,testSmileScore,genArray,pcsUse)


    figure;hold on

    string1 = strcat('PC',num2str(pcsUse(1)));
    string2 = strcat('PC',num2str(pcsUse(2)));
    string3 = strcat('PC',num2str(pcsUse(3)));
    subplot(1,3,1);
    scatter(tSmileScore(:,pcsUse(1)),tSmileScore(:,pcsUse(2)),'filled');hold on;
    % scatter(score(:,1),score(:,2),'filled');hold on;
    scatter(testSmileScore(:,pcsUse(1)),testSmileScore(:,pcsUse(2)),250,'filled','p');
    scatter(genArray(:,pcsUse(1)),genArray(:,pcsUse(2)),'m','filled');
    xlabel(string1);
    ylabel(string2);
    % legend('Training Smiles','Target Smiles','Sampled Smiles');
    subplot(1,3,2);
    scatter(tSmileScore(:,pcsUse(1)),tSmileScore(:,pcsUse(3)),'filled');hold on;
    % scatter(score(:,1),score(:,3),'filled');hold on;
    scatter(testSmileScore(:,pcsUse(1)),testSmileScore(:,pcsUse(3)),250,'filled','p');
    scatter(genArray(:,pcsUse(1)),genArray(:,pcsUse(3)),'m','filled');
    xlabel(string1);
    ylabel(string3);
    % legend('Training Smiles','Target Smiles','Sampled Smiles');
    subplot(1,3,3);
    scatter(tSmileScore(:,pcsUse(2)),tSmileScore(:,pcsUse(3)),'filled');hold on;
    % scatter(score(:,2),score(:,3),'filled');hold on;
    scatter(testSmileScore(:,pcsUse(2)),testSmileScore(:,pcsUse(3)),250,'filled','p');
    scatter(genArray(:,pcsUse(2)),genArray(:,pcsUse(3)),'m','filled');
    xlabel(string2);
    ylabel(string3);
    legend('Training Smiles','Target Smiles','Sampled Smiles');
    

end