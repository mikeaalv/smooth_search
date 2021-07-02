function [obj]=obj_sigma(y_raw,A_hat,sigma)
% Calculting the objective function for smoothing spline when sigma is known.
% when sigma is known, refer to GCV method
% Reference:
% Craven, P. & Wahba, G. Smoothing Noisy Data with Spline Functions - Estimating the Correct Degree of Smoothing by the Method of Generalized Cross-Validation. Numer Math 31, 377-403 (1979).
%
% Argument:
%     y_raw: 1d numeric column vector. The raw y value. Must be provided
%     A_hat: numeric matirx, [n n] n=length(y_raw). The hat matrix that transform y to y_estmated (y_smooth). Must be provided
%     sigma: numeric value. The sigma of noise distribution in y. Must be provided
%
% Return:
%     obj: a numeric value measuring the deviation of smoothed curve and ground truth function. the lower, the better.
% Test:
% NEED TO BE ADDED
% Yue Wu 06112021

if ~exist('y_raw','var')
  error('y_raw is needed');
end
if ~exist('A_hat','var')
  error('A_hat is needed');
end
if ~exist('sigma','var')
  error('sigma is needed');
end

%
n=length(y_raw);
I=eye(n);
part1=sum(((I-A_hat)*y_raw).^2)/n;
part2=sigma^2*(trace(I-A_hat).^2)/n;
part3=sigma^2*(trace(A_hat^2))/n;
obj=part1-part2+part3;
end
