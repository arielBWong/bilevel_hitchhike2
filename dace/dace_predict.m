function [y, sigma] = dace_predict(X, dmodel)
% DACEPREDICT - To predict using kriging model
%
% Call
%    y = dace_predict(X, rmodel)
%
% Input
% X      : Data Points
% rmodel : Kriging model obtained using dace_model
%
% Output:
% y   : Predicted response
%

% Check arguments
if nargin ~= 2
	error('dace_predict requires 2 input arguments')
end
samples = size(X,1);
if samples > 1
    [f, ssqr] = predictor(X, dmodel);
else 
    [f, ~, ssqr] = predictor(X, dmodel);
end
if size(f,1) == size(X,1)
	y = f;
else
	y = f';
end
sigma = sqrt(abs(ssqr));
return
