function [cost, grad] = softmaxCost(theta, numClasses, inputSize, lambda, data, labels)

% numClasses - the number of classes 
% inputSize - the size N of the input vector
% lambda - weight decay parameter
% data - the N x M input matrix, where each column data(:, i) corresponds to
%        a single test set
% labels - an M x 1 matrix containing the labels corresponding for the input data
%

% Unroll the parameters from theta
theta = reshape(theta, numClasses, inputSize);

numCases = size(data, 2);

groundTruth = full(sparse(labels, 1:numCases, 1,numClasses,numCases));
cost = 0;

thetagrad = zeros(numClasses, inputSize);

%% ---------- YOUR CODE HERE --------------------------------------
%  Instructions: Compute the cost and gradient for softmax regression.
%                You need to compute thetagrad and cost.
%                The groundTruth matrix might come in handy.

M1 = theta * data;
M1 = bsxfun(@minus, M1, max(M1, [], 1)); %%Preventing overflows 
M2 = exp(M1);
M2 = bsxfun(@rdivide, M2, sum(M2));

cost = cost -  mean(sum(groundTruth.*log( M2 )));
cost = cost + lambda/2 * sum(sum(theta.^2)); %权重衰减

%%%方法一
% for j=1:numClasses  %%% 当numClasses较大的时候，可以考虑将for改为parfor进行并行计算加速；当numClasses较小的时候（如10）就没有必要并行，并行会慢些
%     thetagrad(j,:) = -mean( bsxfun(@times,groundTruth(j,:)-M2(j,:),data),2 )'+lambda*theta(j,:);  
% end


%%%方法二
% temp1=repmat(groundTruth-M2,[1,1,inputSize]) .* permute( repmat(data,[1,1,numClasses]),[3 2 1] );
% temp2=-mean(temp1,2);
% temp3=squeeze(temp2);
% thetagrad1 = temp3 + lambda*theta;


%%%方法三
%%%本来以为方法二相对方法一会更快些，经过测试，当size(data)为784*60000，类别数为10时，方法一1s左右，方法二15s左右，方法三0.07s左右。三种方法速度不同，但是结果是完全相同的
%%%因为方法三最快，因此用方法三
thetagrad = -1/numCases * (groundTruth - M2) * data' + lambda * theta;


%length(find(thetagrad-thetagrad))


% ------------------------------------------------------------------
% Unroll the gradient matrices into a vector for minFunc
grad = [thetagrad(:)];
end

