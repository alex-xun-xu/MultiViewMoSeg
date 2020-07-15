%% function to adaptively e-nn sparsify affinity matrix
%
%   The threshold is chosen as the row-wise entropy. Details refer to [1]
%   
%   [1] T. Lai, et al. Motion Segmentation Via a Sparsity Constraint. IEEE
%   Transactions on Intelligent Transportation Systems, 2017
%
%   By: Xun Xu, Sep 2018
%   

function [A_sparse,SparseMask] = func_Adapt_eNN(A,alpha)

if ~exist('alpha','var')
    alpha = 5;
end

%% Matrix Computation
sigma = repmat(max(A.^alpha,[],2),1,size(A,2)) - A.^alpha;
Prob = sigma./repmat(sum(sigma,2),1,size(A,2));

E = sum(Prob.*log(Prob+eps),2);

SparseMask = log(Prob+eps) < repmat(E,1,size(A,2));

SparseMask = (SparseMask+SparseMask')/2;

A_sparse = SparseMask.*A;