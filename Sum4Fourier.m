%==========================================================================
%                         FUNCTION Sum4Fourier
%==========================================================================
%
% PURPOSE:
%   Calculates a coordinate value (X, Y, or Z) at a specified distance
%   along a curve by summing the terms of its Fourier series. This function
%   effectively reconstructs the shape from its Fourier coefficients.
%
% INPUTS:
%   TDist (double):
%       The total arc length of the curve.
%   Coeffs (matrix):
%       An M-by-2 matrix of Fourier coefficients for a single dimension
%       (e.g., CoefsX). Column 1 is sine, column 2 is cosine.
%   CuDist (double):
%       The specific cumulative distance along the curve at which to
%       calculate the coordinate value.
%   Wn (vector):
%       A weighting vector (smoothing window) of length M. Each weight Wn(i)
%       is applied to the i-th Fourier coefficient, allowing for smoothing.
%
% OUTPUTS:
%   Val (double):
%       The calculated coordinate value at the specified distance.
%
% AUTHORS:
%   Gadi Herzlinger and Uzy Smilansky.
%
%==========================================================================
function [Val] = Sum4Fourier(TDist, Coeffs, CuDist, Wn)
%Sum4Fourier sums the terms in the Fourier Transform function
%   It receives the total normalized distance (=1), the coefficients
%   (either X or Y) and the cumulative distance along the curve. It also
%   recives the smoothing factor Wn. It returns the sum of all terms in the equation.
%   It is operated by FourierCoords.m

% Loop through each coefficient (frequency).
for i = 1:size(Coeffs,1)
    % Calculate the value of the i-th term in the Fourier series,
    % applying the corresponding smoothing weight from Wn.
    % The formula is: Wn(i) * (A_i*sin(...) + B_i*cos(...))
    Eq(i,1) = Wn(i)*(Coeffs(i,1)*sin(2*pi/TDist*(i-1)*CuDist)+Coeffs(i,2)*cos(2*pi/TDist*(i-1)*CuDist));
end

% Sum all terms to get the final coordinate value.
Val = sum(Eq);

end