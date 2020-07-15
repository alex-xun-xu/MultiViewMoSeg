function [ r P ] = fundamental_degen(X)
        
Fvgg = fundamental_fit(X);

if isempty(Fvgg)
    % Empty output implies degeneracy.
    r = 1;
    P = [];
else    
    % Store the (potentially) 3 solutions in a cell array
    r = 0;
    [~,~,Nsolutions] = size(Fvgg);
    P = cell(Nsolutions,1);
    for n = 1:Nsolutions
        P{n} = Fvgg(:,:,n);
    end
end

end
