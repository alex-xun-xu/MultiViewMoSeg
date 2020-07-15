%% function to do subset constrained clustering
%
%   Input: K - cell arrays of kernels
%             nMotion - number of motions
%             gamma - subset-regularization parameter
%
%
function [U_Subset,itr,loss,exitinfo] = func_Subset_eig(K,nMotion,gamma,epsilon,MaxItr)

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
L = [];

for k_i = 1:nKernel
    
    D = diag(sum(K(:,:,k_i),1));
    L(:,:,k_i) = eye(size(D)) - D^-0.5 * K(:,:,k_i) *D^-0.5;

    [Vec,D] = eig((L(:,:,k_i)+L(:,:,k_i)')/2);
    
    [~,idx] = sort(diag(D),'ascend');
    
    U(:,:,k_i) = Vec(:,idx(1:nMotion));
end

%% Do Subset Constrained Clustering
U_Subset = U;

K_hat = [];
for k_i = 1:nKernel
    K_hat(:,:,k_i) = U(:,:,k_i)*U(:,:,k_i)';
end

itr = 1;

while true
    
    for k_i = 1:nKernel
        
        %%% Compute CoReg Kernels
        
        %         UUt_k = sum(UUt,3) - U_CoReg(:,:,k_i)*U_CoReg(:,:,k_i)';
        if k_i == 1
            tmp_K_neg = K_hat(:,:,k_i+1);
            Q_k = tmp_K_neg.*(tmp_K_neg<0);
        elseif k_i == nKernel
            tmp_K_pos = K_hat(:,:,k_i-1);
            Q_k = tmp_K_pos.*(tmp_K_pos>0);
        else
            tmp_K_neg = K_hat(:,:,k_i+1);
            tmp_K_pos = K_hat(:,:,k_i-1);
            Q_k = tmp_K_pos.*(tmp_K_pos>0)+tmp_K_neg.*(tmp_K_neg<0);
        end
        
%         L_tilde(:,:,k_i) = L(:,:,k_i) - gamma/(sqrt(nMotion)*norm(Q_k,'fro')+eps) * Q_k;
        L_tilde(:,:,k_i) = L(:,:,k_i) - gamma * Q_k;

        L_tilde(:,:,k_i) = (L_tilde(:,:,k_i)+L_tilde(:,:,k_i)')/2;
        
        [Vec,D] = eig(L_tilde(:,:,k_i));
        
        [~,idx] = sort(diag(D),'ascend');
        
        U_Subset(:,:,k_i) = Vec(:,idx(1:nMotion));
        
        
        %         [U_tmp,S,V] = svd(L_tilde(:,:,k_i));
        %
        %         U_Subset(:,:,k_i) = U_tmp(:,end-nMotion+1:end);
        
        K_hat(:,:,k_i) = U_Subset(:,:,k_i)*U_Subset(:,:,k_i)';
        
    end
    
    %% Evaluate Loss
    loss(itr) = Loss(U_Subset,L,gamma);
    
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


function Res = Loss(U_Subset,L,gamma)

Res = 0;

nKernel = size(U_Subset,3);
nMotion = size(U_Subset,2);

K_hat = [];
for k_i = 1:nKernel
    K_hat(:,:,k_i) = U_Subset(:,:,k_i)*U_Subset(:,:,k_i)';
end

for k_i = 1:size(U_Subset,3)
    
    if k_i == 1
        tmp_K_neg = K_hat(:,:,k_i+1);
        Q_k = tmp_K_neg.*(tmp_K_neg<0);
    elseif k_i == nKernel
        tmp_K_pos = K_hat(:,:,k_i-1);
        Q_k = tmp_K_pos.*(tmp_K_pos>0);
    else
        tmp_K_neg = K_hat(:,:,k_i+1);
        tmp_K_pos = K_hat(:,:,k_i-1);
        Q_k = tmp_K_pos.*(tmp_K_pos>0)+tmp_K_neg.*(tmp_K_neg<0);
    end
    
    L_tilde(:,:,k_i) = L(:,:,k_i) - gamma * Q_k;
    
    
    Res = Res + trace(U_Subset(:,:,k_i)'*L_tilde(:,:,k_i)*U_Subset(:,:,k_i));
    
end





