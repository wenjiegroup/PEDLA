function [data] = unsupervisedstackedAEPredict(theta, netconfig, data)
                                         
% stackedAEPredict: Takes a trained theta and a test data set,
% and returns the predicted labels for each example.
                                         
% theta: trained weights from the autoencoder
% hiddenSize:  the number of hidden units *at the 2nd layer*
% numClasses:  the number of categories
% data: Our matrix containing the training data as columns.  So, data(:,i) is the i-th training example. 

% Your code should produce the prediction matrix 
% pred, where pred(i) is argmax_c P(y(c) | x(i)).
 
%% Unroll theta parameter

% Extract out the "stack"
stack = params2stack(theta, netconfig);

%% ---------- YOUR CODE HERE --------------------------------------
%  Instructions: Compute pred using theta assuming that the labels start 
%                from 1.

depth = numel(stack);
if rem(depth,2) ~= 0
    error('the depth is not correct!');
end
for d=1:depth/2
    data = sigmoid( bsxfun(@plus,stack{d}.w*data,stack{d}.b) );
end

% pred = softmaxTheta * data;
% [~,pred] = max(pred);
% 
% M1 = softmaxTheta * data;
% M1 = bsxfun(@minus, M1, max(M1, [], 1)); %%Preventing overflows 
% M2 = exp(M1);
% M2 = bsxfun(@rdivide, M2, sum(M2));
% [~,pred] = max(M2);






% -----------------------------------------------------------

end


% You might find this useful
function sigm = sigmoid(x)
    sigm = 1 ./ (1 + exp(-x));
end
