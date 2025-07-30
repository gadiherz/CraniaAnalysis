%==========================================================================
%                           FUNCTION TanCalc
%==========================================================================
%
% PURPOSE:
%   Calculates the tangent angle (theta, in radians) at one or more points
%   along a 2D curve that is represented by a Fourier series. The calculation
%   is performed analytically using the first derivatives of the series.
%
% METHODOLOGY:
%   The function uses the relationship Theta = atan(dy/dx), where dy/dx is
%   the slope of the tangent line. The derivatives dy/ds and dx/ds (where 's'
%   is the arc length parameter) are calculated by summing the analytical
%   derivatives of each term in the Fourier series.
%
% INPUTS:
%   TDist (double):
%       The total arc length of the curve.
%   Dists (vector):
%       A vector of cumulative distances along the curve at which to
%       calculate the tangent angle.
%   Wn (vector):
%       A weighting vector (smoothing window) applied to the Fourier
%       coefficients, allowing for smoothed tangent calculation.
%   CoeffsX, CoeffsY (matrix):
%       M-by-2 matrices of Fourier coefficients for the X and Y dimensions.
%
% OUTPUTS:
%   Theta (vector):
%       A vector of tangent angles (in radians) corresponding to each
%       point in 'Dists'.
%
% AUTHORS:
%   Gadi Herzlinger and Uzy Smilansky.
%
%==========================================================================
function [Theta] = TanCalc(TDist,Dists,Wn,CoeffsX, CoeffsY)

%% --- DERIVATIVE CALCULATION LOOP --- %%

% Loop through each distance specified in the input vector 'Dists'.
for j = 1:size(Dists,1)
    % Initialize the summed derivatives for the current point to zero.
    Sdxds=0; % Sum of first derivative of X (dx/ds)
    Sdyds=0; % Sum of first derivative of Y (dy/ds)

    % Loop through each Fourier coefficient (frequency) to calculate its contribution.
    for i = 1:size(CoeffsX,1)
        % Analytically calculate the 1st derivatives for the i-th term.
        dxds = -2*pi*(i-1)*Wn(i)*(CoeffsX(i,2)*sin(2*pi*(i-1)*(Dists(j)./TDist)) - CoeffsX(i,1)*cos(2*pi*(i-1)*(Dists(j)./TDist)))./TDist;
        dyds = -2*pi*(i-1)*Wn(i)*(CoeffsY(i,2)*sin(2*pi*(i-1)*(Dists(j)./TDist)) - CoeffsY(i,1)*cos(2*pi*(i-1)*(Dists(j)./TDist)))./TDist;

        % Add the contribution of the i-th term to the total sum.
        Sdxds = Sdxds+dxds;
        Sdyds = Sdyds+dyds;
    end
    
    %% --- ANGLE CALCULATION --- %%
    
    % The tangent angle (theta) is the arctangent of the quotient of the
    % y-derivative and the x-derivative (atan(dy/dx)).
    Theta(j,1) = atan(Sdyds/Sdxds);
end
end