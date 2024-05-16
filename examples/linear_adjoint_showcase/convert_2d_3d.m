clear all
clc

addpath('/scratch/nobis/readwritenek')
[data_2D,lr1,elmap,time,istep,fields,emode,wdsz,etag,header,status] = readnek('baseflow0.f00000');
[data_3D,lr1,elmap,time,istep,fields,emode,wdsz,etag,header,status] = readnek(   'field0.f00000');

for i=1:size(data_3D,1)
    data_3D(i,:,4) = repmat(data_2D(i,:,3),1,6); %u
    data_3D(i,:,5) = repmat(data_2D(i,:,4),1,6); %v
end


writenek('harry0.f00001',data_3D,lr1,elmap,time,istep,fields,emode,wdsz,etag)