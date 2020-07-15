%% function to do co-regularization
%
%   Input: K - cell arrays of kernels
%          nMotion - number of motions
%          lambda - co-regularization parameter
%
function [U_CoReg, itr ,loss, exitinfo] = func_CoRegularize_eig(K,nMotion,lambda,epsilon,MaxItr)

%% Convergence Parameters
if ~exist('MaxItr','var')
    MaxItr = 30;
end

if ~exist('epsilon','var')
    epsilon = 1e-8;
end

nKernel = size(K,3);

%% Intialize Each Spectral Embedding
U = [];
for k_i = 1:nKernel
    
    D = diag(sum(K(:,:,k_i),1));
    L(:,:,k_i) = eye(size(D)) - D^-0.5 * K(:,:,k_i) *D^-0.5;
    
    [U_tmp,S,V] = svd(L(:,:,k_i));
    
    U(:,:,k_i) = U_tmp(:,end-nMotion+1:end);
    
end

%% Do Co-Regularization
exitinfo.reason = '';

U_CoReg = U;

UUt = [];
for k_i = 1:nKernel
    UUt(:,:,k_i) = U(:,:,k_i)*U(:,:,k_i)';
end

itr = 1;

while true
    
    for k_i = 1:nKernel
        
        %%% Compute CoReg Kernels
        
        UUt_k = sum(UUt,3) - U_CoReg(:,:,k_i)*U_CoReg(:,:,k_i)';
%         L_CoReg(:,:,k_i) = L(:,:,k_i) - lambda/nMotion * UUt_k;
        L_CoReg(:,:,k_i) = L(:,:,k_i) - lambda * UUt_k;

%         [U_tmp,S,V] = svd(L_CoReg(:,:,k_i));
        
        [Vec,D] = eig((L_CoReg(:,:,k_i)'+L_CoReg(:,:,k_i))/2);
        
        [~,idx] = sort(diag(D),'ascend');
        
        U_CoReg(:,:,k_i) = Vec(:,idx(1:nMotion));
        
%         U_CoReg(:,:,k_i) = U_tmp(:,end-nMotion+1:end);
        
        UUt(:,:,k_i) = U_CoReg(:,:,k_i)*U_CoReg(:,:,k_i)';
        
    end
    
    %% Evaluate Loss
    loss(itr) = Loss(U_CoReg,L_CoReg);
    
    %% Check Convergence
    %%% Check Loss Change
    if itr > 1 && abs(loss(end-1)-loss(end)) < epsilon
        exitinfo.reason = 'converge';
        break;
    end
    
    %%% Exceed Max Iteration
    if itr > MaxItr
        exitinfo.reason = 'timeout';
        break;
    end
    
    %% Display Results
    if itr == 1
        string = sprintf('loss = %.5f; change of loss = %.5f\n',loss(itr),0);
        fprintf('%s',string);
    else
        %         for i = 1:length(string)
        %             fprintf('\b');
        %         end
        string = sprintf('loss = %.5f; change of loss = %.5f\n',loss(itr),loss(itr)-loss(itr-1));
        fprintf('%s',string);
    end
    itr = itr+1;
    
end


function L = Loss(U_CoReg,L_CoReg)

L = 0;

for k_i = 1:size(U_CoReg,3)
    
    L = L + trace(U_CoReg(:,:,k_i)'*L_CoReg(:,:,k_i)*U_CoReg(:,:,k_i));
    
end





