%==========================================================================
%                          FUNCTION CurvCalc
%==========================================================================
%
% PURPOSE:
%   Calculates the 2D curvature at one or more points along a curve that is
%   represented by a Fourier series. The calculation is performed
%   analytically using the derivatives of the Fourier series.
%
% METHODOLOGY:
%   The function uses the standard formula for 2D curvature (K):
%   K = (x'*y'' - y'*x'') / (x'^2 + y'^2)^(3/2)
%   where x' and x'' are the first and second derivatives of the x-coordinate
%   with respect to the arc length parameter 's'.
%
%   Because the curve is parameterized by arc length, the denominator term
%   sqrt(x'^2 + y'^2), which is the magnitude of the tangent vector, is
%   equal to 1. The formula thus simplifies to:
%   K = x'*y'' - y'*x''
%   (The code uses the negative of this, K = -(y'*x'' - x'*y''), due to a
%   sign convention for concavity/convexity).
%
% INPUTS:
%   TDist (double):
%       The total arc length of the curve.
%   Dists (vector):
%       A vector of cumulative distances along the curve at which to
%       calculate the curvature.
%   Wn (vector):
%       A weighting vector (smoothing window) applied to the Fourier
%       coefficients, allowing for smoothed curvature calculation.
%   CoeffsX, CoeffsY (matrix):
%       M-by-2 matrices of Fourier coefficients for the X and Y dimensions.
%
% OUTPUTS:
%   K (vector):
%       A vector of curvature values corresponding to each point in 'Dists'.
%
% AUTHORS:
%   Gadi Herzlinger and Uzy Smilansky.
%
%==========================================================================
function [K] = CurvCalc(TDist,Dists,Wn,CoeffsX, CoeffsY)

% Loop through each distance specified in the input vector 'Dists'.
for j = 1:size(Dists,1)
    % Initialize the summed derivatives for the current point to zero.
    Sdxds=0;    % Sum of first derivative of X
    Sdyds=0;    % Sum of first derivative of Y
    Sd2xds2=0;  % Sum of second derivative of X
    Sd2yds2=0;  % Sum of second derivative of Y
    
    % Loop through each Fourier coefficient (frequency) to calculate its contribution.
    for i = 1:size(CoeffsX,1)
        % Analytically calculate the 1st derivatives for the i-th term.
        dxds = -2*pi*(i-1)*Wn(i)*(CoeffsX(i,2)*sin(2*pi*(i-1)*(Dists(j)./TDist)) - CoeffsX(i,1)*cos(2*pi*(i-1)*(Dists(j)./TDist)))./TDist;
        dyds = -2*pi*(i-1)*Wn(i)*(CoeffsY(i,2)*sin(2*pi*(i-1)*(Dists(j)./TDist)) - CoeffsY(i,1)*cos(2*pi*(i-1)*(Dists(j)./TDist)))./TDist;
        
        % Analytically calculate the 2nd derivatives for the i-th term.
        d2xds2 = -4*pi.^2*(i-1).^2*Wn(i)*(CoeffsX(i,1)*sin(2*pi*(i-1)*(Dists(j)./TDist)) + CoeffsX(i,2)*cos(2*pi*(i-1)*(Dists(j)./TDist)))./TDist.^2;
        d2yds2 = -4*pi.^2*(i-1).^2*Wn(i)*(CoeffsY(i,1)*sin(2*pi*(i-1)*(Dists(j)./TDist)) + CoeffsY(i,2)*cos(2*pi*(i-1)*(Dists(j)./TDist)))./TDist.^2;
        
        % Add the contribution of the i-th term to the total sum.
        Sdxds = Sdxds+dxds;
        Sdyds = Sdyds+dyds;
        Sd2xds2 = Sd2xds2+d2xds2;
        Sd2yds2 = Sd2yds2+d2yds2;
    end
    
    % Calculate the final curvature at point j using the summed derivatives.
    % This implements the simplified formula K = -(y'*x'' - x'*y'').
    K(j,1) = -(Sdyds*Sd2xds2-Sdxds*Sd2yds2)/(Sdxds.^2+Sdyds.^2);
end

end