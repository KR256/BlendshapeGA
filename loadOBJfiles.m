function loadOBJfiles(OBJ_folderName,outFilename)

objFilenames = dir(strcat('Dataset/',OBJ_folderName,'/*.obj'));
numFiles = length(objFilenames);

dataset = [];
for f=1:numFiles
    filename = objFilenames(f).name;
    fprintf('Reading %s\n', filename); 
    obj = read_wobj(filename);
    dataset = [dataset ; reshape(obj.vertices',1,[])];
    
end

save(strcat('Resources/',outFilename,'.mat'),'dataset');

end