
function writeMesh(vertices,fileName)

load('Tools/OutputShape.mat','OutputShape');

outFace = vertices;
outputMat = vec2mat(outFace,3);
OutputShape.vertices = outputMat;
OutputShape.vertices_normal = outputMat;

write_wobj(OutputShape,fileName);

end