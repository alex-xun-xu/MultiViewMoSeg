%% function to L2 normalize feature

function Feat = func_L2Normalize(Feat,dim)

if dim == 1
    
    Feat = Feat./repmat(sqrt(sum(Feat.^2,dim))+eps,size(Feat,1),1);
    
elseif dim == 2
    
    Feat = Feat./repmat(sqrt(sum(Feat.^2,dim))+eps,1,size(Feat,2));
end