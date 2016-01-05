function [ path, logP ] = lf_HMMViterbi( PI, A, B )
%LF_HMMVITERBI Find the most-probable (Viterbi) path through the HMM state trellis
%   Detailed explanation goes here
% Inputs:
% prior(i) = Pr(Q(1) = i)
% transmat(i,j) = Pr(Q(t+1)=j | Q(t)=i)
% obslik(i,t) = Pr(y(t) | Q(t)=i)
%
% Outputs:
% path(t) = q(t), where q1 ... qT is the argmax of the above expression.
% logP:log of p
%
% Intermediate variable:
% delta(j,t) = prob. of the best sequence of length t-1 and then going to state j, and O(1:t)
% phi(j,t) = log of delta(j,t) for scaling
% psi(j,t) = the best predecessor state, given that we ended up in state j at t
% N:number of states;
% T:length of squence


% A must be square
N = size(A,1);
checkA = size(A,2);
if checkA ~= N
    error(message('stats:lf_HMMViterbi:BadTransitions'));
end

% check that PI is a row vector, and that number of rows of PI must be same as number of states
if ~isvector(PI)
    error(message('stats:lf_HMMViterbi:PI must be a vector'));
end
PI=PI(:);  % make sure PI is a row vector
checkPI = size(PI,1);
if checkPI ~= N
    error(message('stats:lf_HMMViterbi:PI InputSizeMismatch'));
end

% number of rows of B must be same as number of states
checkB = size(B,1);
if checkB ~= N
    error(message('stats:lf_HMMViterbi:B InputSizeMismatch'));
end



T = size(B, 2);
phi = zeros(N,T);   % log of delta
psi = zeros(N,T);
path = zeros(1,T);

t=1;
phi(:,t) = log(PI) + log(B(:,t));   %  use log for scaling
psi(:,t) = 0; % arbitrary value, since there is no predecessor to t=1
for t=2:T
    [phi(:,t), psi(:,t)] = max( bsxfun(@plus,phi(:,t-1),log(A)) );
    phi(:,t) = phi(:,t) + log(B(:,t));    
end

[logP, path(T)] = max(phi(:,T));
for t=T-1:-1:1
  path(t) = psi(path(t+1),t+1);
end





end

