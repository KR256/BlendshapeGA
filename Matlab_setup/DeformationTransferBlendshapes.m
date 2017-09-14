
function neutralVec = DeformationTransferBlendshapes(INPUT_FACE_ID,trainOrTest,...
    ignoreList,EXPRESSION_TYPE)
% Deformation Transfer

if(strcmp(trainOrTest,'test')==1)
    OUTPUT_PATH = strcat('Blendshapes/',EXPRESSION_TYPE,'/deformedTestTargets/FACE_',num2str(INPUT_FACE_ID),'/');
else
    OUTPUT_PATH = strcat('Blendshapes/',EXPRESSION_TYPE,'/deformedTrainingTargets/');
end

fprintf('Performing Deformation Transfer\n');

tedBlendshapeObjs = dir('TedTargets/*.obj');

load('TedTargets/neutralTed.mat','Source1');

identityObjs = dir('Dataset/OBJs_neutrals/*.obj');
SOURCE_FILENAME = identityObjs(INPUT_FACE_ID).name;
fprintf('Reading %s\n', SOURCE_FILENAME); 
Source2 = read_wobj(strcat('Dataset/OBJs_neutrals/',SOURCE_FILENAME));

NumberOfVertices = size(Source1.vertices,1);

NumberOfTriangles = size(Source1.objects(4).data.vertices,1); 

UNMOVED_VERTEX = 1;


CoefficientMatrix = zeros(NumberOfTriangles * 3, NumberOfVertices);


    %Create the Coefficient Matrix for solving the linear equation 
    %to get the vertices of the blendshapes
    for faceNumber = 1:NumberOfTriangles

        % calculate Mj
        Index1 = Source1.objects(4).data.vertices(faceNumber,1);
        Index2 = Source1.objects(4).data.vertices(faceNumber,2);
        Index3 = Source1.objects(4).data.vertices(faceNumber,3);
        
        rowOffset = (faceNumber - 1) *3;

        CoefficientMatrix(rowOffset + 1,Index2) = 1;
        CoefficientMatrix(rowOffset + 1,Index1) = -1;

        CoefficientMatrix(rowOffset + 2,Index3) = 1;
        CoefficientMatrix(rowOffset + 2,Index1) = -1;        
        

    end
    
    
    CoefficientMatrixTranspose = single(transpose(CoefficientMatrix));
    CoefficientMatrix_new = CoefficientMatrixTranspose * single(CoefficientMatrix);    
    
    clear CoefficientMatrix;
    

%     SOURCE_FILENAME = IDENTITY_FILENAME;
%     fprintf('Reading %s\n', SOURCE_FILENAME); 
%     Source2 = read_wobj(strcat('processedNeutrals/',SOURCE_FILENAME));

bsVertsMat = zeros(length(tedBlendshapeObjs) - length(ignoreList)-1,...
    NumberOfVertices*3);
% Source1.vertices(:,3) = Source1.vertices(:,3) + 200;
% Source2.vertices(:,3) = Source2.vertices(:,3) + 200;
% ignoreList = [10 13:23 27:38 57:58 63 67:70 92:99 104:115 134 139];

for fileNumber = 1:length(tedBlendshapeObjs)
    if(ismember(fileNumber,ignoreList))
        continue;
    end
    
    BLENDSHAPE_FILENAME = tedBlendshapeObjs(fileNumber).name;
    

%     cd('TedTargets');
    
    
    Target1 = read_wobj(strcat('TedTargets/',BLENDSHAPE_FILENAME));
    Target1.vertices = Target1.vertices * 10;
    Target1.vertices(:,2) = Target1.vertices(:,2) - 1749.572;

    Output = Target1;
    
    
    DifferenceInVertices = Source1.vertices - Source2.vertices;
    DifferenceInVertices = sum(DifferenceInVertices,2);
%     IndexOfUnmovedVertex = find(abs(DifferenceInVertices) < 0.5);
%     UNMOVED_VERTEX = IndexOfUnmovedVertex(1);
    
    NumberOfTriangles = size(Source1.objects(4).data.vertices,1);

    Target_Blendshape_Gradient_Matrix = zeros(NumberOfTriangles * 3, 3);


    disp('Calculating B Matrix');

    % go through all triangles and create the B matrix
    for vertexNumber = 1:NumberOfTriangles                

        v1= zeros(1,3);
        v2= zeros(1,3);
        v3= zeros(1,3);    

        index = Source1.objects(4).data.vertices(vertexNumber,1);
        v1(1,1) = Source1.vertices(index,1);
        v1(1,2) = Source1.vertices(index,2);
        v1(1,3) = Source1.vertices(index,3);

        index = Source1.objects(4).data.vertices(vertexNumber,2);
        v2(1,1) = Source1.vertices(index,1);
        v2(1,2) = Source1.vertices(index,2);
        v2(1,3) = Source1.vertices(index,3);

        index = Source1.objects(4).data.vertices(vertexNumber,3);
        v3(1,1) = Source1.vertices(index,1);
        v3(1,2) = Source1.vertices(index,2);
        v3(1,3) = Source1.vertices(index,3);

        v4 = v1 + cross((v2 - v1),(v3 - v1))/sqrt(norm(cross((v2-v1),(v3-v1))));



        VMatrix = zeros(3,3);
        firstRow = v2-v1;
        secondRow = v3-v1;
        thirdRow = v4-v1;

        VMatrix(1,1) = firstRow(1,1);
        VMatrix(2,1) = firstRow(1,2);
        VMatrix(3,1) = firstRow(1,3);

        VMatrix(1,2) = secondRow(1,1);
        VMatrix(2,2) = secondRow(1,2);
        VMatrix(3,2) = secondRow(1,3);

        VMatrix(1,3) = thirdRow(1,1);
        VMatrix(2,3) = thirdRow(1,2);
        VMatrix(3,3) = thirdRow(1,3);  

        %VInverse = inv(VMatrix);              


        v1Bar= zeros(1,3);
        v2Bar= zeros(1,3);
        v3Bar= zeros(1,3);    

        index = Source2.objects.data.vertices(vertexNumber,1);
        v1Bar(1,1) = Source2.vertices(index,1);
        v1Bar(1,2) = Source2.vertices(index,2);
        v1Bar(1,3) = Source2.vertices(index,3);

        index = Source2.objects.data.vertices(vertexNumber,2);
        v2Bar(1,1) = Source2.vertices(index,1);
        v2Bar(1,2) = Source2.vertices(index,2);
        v2Bar(1,3) = Source2.vertices(index,3);

        index = Source2.objects.data.vertices(vertexNumber,3);
        v3Bar(1,1) = Source2.vertices(index,1);
        v3Bar(1,2) = Source2.vertices(index,2);
        v3Bar(1,3) = Source2.vertices(index,3);

        v4Bar = v1Bar + cross((v2Bar - v1Bar),(v3Bar - v1Bar))/sqrt(norm(cross((v2Bar-v1Bar),(v3Bar-v1Bar))));


        VBarMatrix = zeros(3,3);
        firstRow = v2Bar-v1Bar;
        secondRow = v3Bar-v1Bar;
        thirdRow = v4Bar-v1Bar;

        VBarMatrix(1,1) = firstRow(1,1);
        VBarMatrix(2,1) = firstRow(1,2);
        VBarMatrix(3,1) = firstRow(1,3);

        VBarMatrix(1,2) = secondRow(1,1);
        VBarMatrix(2,2) = secondRow(1,2);
        VBarMatrix(3,2) = secondRow(1,3);

        VBarMatrix(1,3) = thirdRow(1,1);
        VBarMatrix(2,3) = thirdRow(1,2);
        VBarMatrix(3,3) = thirdRow(1,3);  

        QMatrix = VBarMatrix / VMatrix;

        %QTranspose = transpose(QMatrix);
        
        
        
        Target_v1= zeros(1,3);
        Target_v2= zeros(1,3);
        Target_v3= zeros(1,3);    

        index1 = Target1.objects.data.vertices(vertexNumber,1);
        Target_v1(1,1) = Target1.vertices(index1,1);
        Target_v1(1,2) = Target1.vertices(index1,2);
        Target_v1(1,3) = Target1.vertices(index1,3);

        index2 = Target1.objects.data.vertices(vertexNumber,2);
        Target_v2(1,1) = Target1.vertices(index2,1);
        Target_v2(1,2) = Target1.vertices(index2,2);
        Target_v2(1,3) = Target1.vertices(index2,3);

        index3 = Target1.objects.data.vertices(vertexNumber,3);
        Target_v3(1,1) = Target1.vertices(index3,1);
        Target_v3(1,2) = Target1.vertices(index3,2);
        Target_v3(1,3) = Target1.vertices(index3,3);

        Target_v4 = Target_v1 + cross((Target_v2 - Target_v1),(Target_v3 - Target_v1))/sqrt(norm(cross((Target_v2-Target_v1),(Target_v3-Target_v1))));


        VMatrix = zeros(3,3);
        firstRow = Target_v2-Target_v1;
        secondRow = Target_v3-Target_v1;
        thirdRow = Target_v4-Target_v1;

        VMatrix(1,1) = firstRow(1,1);
        VMatrix(2,1) = firstRow(1,2);
        VMatrix(3,1) = firstRow(1,3);

        VMatrix(1,2) = secondRow(1,1);
        VMatrix(2,2) = secondRow(1,2);
        VMatrix(3,2) = secondRow(1,3);

        VMatrix(1,3) = thirdRow(1,1);
        VMatrix(2,3) = thirdRow(1,2);
        VMatrix(3,3) = thirdRow(1,3);  

        
        
        Target_Blendshape_Gradient_Triangle =  QMatrix * VMatrix;
        Target_Blendshape_Gradient_Triangle = Target_Blendshape_Gradient_Triangle';
                       

        offset = (vertexNumber-1) * 3 ;

        for VInvColNumber = 1:3
            for col = 1:3
                Target_Blendshape_Gradient_Matrix(VInvColNumber + offset ,col) = Target_Blendshape_Gradient_Triangle(VInvColNumber,col);
            end
        end        


    end




    
    
    for count = 1:NumberOfTriangles * 3
        if(mod(count,3) == 0)        
            for col = 1:3
                Target_Blendshape_Gradient_Matrix(count,col) = 0;
            end        
        end
    end
    
    
     
    
    
    Blend_new = CoefficientMatrixTranspose * Target_Blendshape_Gradient_Matrix;
    
    Vertices = linsolve(CoefficientMatrix_new, Blend_new);   
    

    
    
    Output.vertices = Vertices;
    
    displacement = Target1.vertices(UNMOVED_VERTEX,:) - Output.vertices(UNMOVED_VERTEX,:);
    
    %[error,Reallignedsource,transform]=rigidICP(Target1.vertices,Output.vertices,1);
    
    %Output.vertices = Reallignedsource;
    
    Output.vertices(:,1) = Output.vertices(:,1) + displacement(1,1);
    Output.vertices(:,2) = Output.vertices(:,2) + displacement(1,2);
    Output.vertices(:,3) = Output.vertices(:,3) + displacement(1,3);
   
    if(strcmp(BLENDSHAPE_FILENAME,'Neutral.obj')==1)
        write_wobj(Output,strcat('Blendshapes/neutral/bsNeutrals/FACE',num2str(INPUT_FACE_ID),'_',tedBlendshapeObjs(fileNumber).name));
        neutralVec = vec2mat(Output.vertices,1,3);
    else
        bsVertsMat(fileNumber,:) = reshape(Output.vertices',1,[]);
    end
    
    if(strcmp(trainOrTest,'test')==1)
        write_wobj(Output,strcat(OUTPUT_PATH,'FACE',num2str(INPUT_FACE_ID),'_',tedBlendshapeObjs(fileNumber).name));
    end
    
    
end

% cd(OUTPUT_PATH);
% bsVertsMat( ~any(bsVertsMat,2),: ) = [];
save(strcat(OUTPUT_PATH,'bsVerts52_',num2str(INPUT_FACE_ID),'.mat'),'bsVertsMat');