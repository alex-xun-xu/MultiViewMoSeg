function plotCoef(CKSym, grp,h2axes)

if ~exist('h2axes','var') || isempty(h2axes)
    figure;
    h2axes = axes;
end

CKSym_per = CKSym;
grp_num = max(grp);
start_pos = 1;
for i = 1:grp_num
    ind = find(grp==i);
    for j = 1:length(ind)
        tmp = CKSym_per(start_pos,:);
        CKSym_per(start_pos,:) = CKSym_per(ind(j),:);
        CKSym_per(ind(j),:) = tmp;
        
        tmp = CKSym_per(:,start_pos);
        CKSym_per(:,start_pos) = CKSym_per(:,ind(j));
        CKSym_per(:,ind(j)) = tmp;
        
        tmp = grp(start_pos);
        grp(start_pos) = i;
        grp(ind(j)) = tmp;
        
        start_pos = start_pos + 1;
    end
end

axes(h2axes);
imshow(-CKSym_per,[]);

end