%==========================================================================
%                         FUNCTION FourierCoefs
%==========================================================================
%
% PURPOSE:
%   Calculates the Fourier series coefficients for a 2D or 3D curve defined
%   by an ordered list of coordinates. The function parameterizes the curve
%   by its arc length and can handle both open and closed curves.
%
% METHODOLOGY:
%   1. Resamples the input coordinates to create points that are equidistant
%      along the curve's arc length.
%   2. For closed 2D curves, it standardizes the start point by rotating
%      the curve to a consistent, data-driven position.
%   3. For open curves, it subtracts a 5th-order polynomial that matches the
%      curve's endpoints (position, velocity, and acceleration). This transforms
%      the non-periodic open curve into a periodically smooth function suitable
%      for Fourier analysis, minimizing artifacts like the Gibbs phenomenon.
%   4. It then computes the Fourier coefficients by solving a system of
%      linear equations.
%
% INPUTS:
%   Bcoords (matrix):
%       An N-by-2 or N-by-3 matrix of ordered coordinates along a curve.
%   OpenCurv (logical):
%       A flag indicating if the curve is open (true) or closed (false).
%       Defaults to false.
%
% OUTPUTS:
%   TotLen (double):
%       The total arc length of the final, processed curve.
%   CoefsX, CoefsY, CoefsZ (matrix):
%       M-by-2 matrices of Fourier coefficients for each dimension. Column 1
%       contains sine coefficients; column 2 contains cosine coefficients.
%   GCoefs (matrix):
%       For open curves, a 6-by-3 matrix of the polynomial coefficients
%       used for endpoint matching. Returns NaN for closed curves.
%
% AUTHORS:
%   Gadi Herzlinger and Uzy Smilansky.
%
%==========================================================================
function [TotLen, CoefsX, CoefsY, CoefsZ, GCoefs] = FourierCoefs(Bcoords, OpenCurv)
%% --- INITIALIZATION & RESAMPLING --- %%

% If the 'OpenCurv' flag is not provided, assume the curve is closed by default.
if ~exist('OpenCurv','var')
    OpenCurv = false;
end

% Resample the input coordinates to generate equidistantly spaced points,
% which is a prerequisite for this Fourier analysis method.
% B2: The new set of equidistant points.
% CuDist: The cumulative arc length at each point.
[B2, CuDist] = EquiDistPT(Bcoords, OpenCurv);


%% --- START-POINT STANDARDIZATION (CLOSED 2D CURVES) --- %%

% This loop is only applied to closed 2D curves. It rotates the curve to
% a consistent starting position to ensure rotational invariance.
if size(Bcoords,2) == 2 && ~OpenCurv
    % Center the curve at the origin.
    B3 = B2 - mean(B2);
    % Find all points where the curve crosses the Y-axis (where X changes sign).
    PotPts = [];
    if (B3(end,1) < 0 && B3(1,1) > 0) || (B3(end,1) > 0 && B3(1,1)<0)
        PotPts = [size(B3,1),1];
    end
    iPotPts = [find(any([all([B3(2:end,1) < 0, B3(1:end-1,1) > 0],2),all([B3(2:end,1) > 0, B3(1:end-1,1) < 0],2)],2)),...
        find(any([all([B3(2:end,1) < 0, B3(1:end-1,1) > 0],2),all([B3(2:end,1) > 0, B3(1:end-1,1) < 0],2)],2))+1];
    PotPts = [PotPts; iPotPts];
    % Select the crossing point with the minimum Y value as the reference.
    RC = find(mean(reshape(B3(PotPts,2),size(PotPts)),2)==min(mean(reshape(B3(PotPts,2),size(PotPts)),2)));
    % Shift and flip the curve so the reference point becomes the new start point.
    if B3(PotPts(RC,1),1) < 0 && B3(PotPts(RC,2),1) > 0
        B3 = circshift(B3,-PotPts(RC,2));
    else
        B3 = flip(circshift(B3,-PotPts(RC,2)));
    end
    % Translate the curve back to its original position.
    B2 = B3 + mean(B2);
    clear  PotPts iPotPts RC B3
end

%% --- CURVE END-POINT HANDLING --- %%

if OpenCurv
    % --- For OPEN curves: Use polynomial subtraction to enforce periodic smoothness ---
    % This process removes discontinuities at the endpoints, which is
    % necessary for accurate calculation of derivatives from the Fourier series.
    
    % Ensure data is 3D, adding a NaN Z-column if needed.
    if size(B2,2)==2
        B2 = [B2,nan(size(B2,1))];
    end
    % Trim curve to calculate numerical derivatives at the endpoints.
    CuDist = CuDist(1:end-2);
    
    % Define endpoint conditions for a 5th-order polynomial: position (F),
    % first derivative (Ft), and second derivative (Ftt) at start (0) and end (L).
    F0 = [B2(2,1), B2(2,2), B2(2,3)];
    FL = [B2(end-1,1), B2(end-1,2), B2(end-1,3)];
    Ft0 = [(B2(3,1)-B2(1,1))./(2.*CuDist(2)),(B2(3,2)-B2(1,2))./(2.*CuDist(2)),(B2(3,3)-B2(1,3))./(2.*CuDist(2))];
    FtL = [(B2(end,1)-B2(end-2,1))./(2.*CuDist(2)),(B2(end,2)-B2(end-2,2))./(2.*CuDist(2)),...
        (B2(end,3)-B2(end-2,3))./(2.*CuDist(2))];
    Ftt0 = [(B2(3,1)-2.*B2(2,1)+B2(1,1))./CuDist(2).^2,(B2(3,2)-2.*B2(2,2)+B2(1,2))./CuDist(2).^2,(B2(3,3)-2.*B2(2,3)+B2(1,3))./CuDist(2).^2];
    FttL = [(B2(end-2,1)-2.*B2(end-1,1)+B2(end,1))./CuDist(2).^2,...
        (B2(end-2,2)-2.*B2(end-1,2)+B2(end,2))./CuDist(2).^2,...
        (B2(end-2,3)-2.*B2(end-1,3)+B2(end,3))./CuDist(2).^2,];
    
    % Solve the system of linear equations to find the polynomial coefficients.
    AMat = [1,CuDist(1),CuDist(1).^2,CuDist(1).^3,CuDist(1).^4,CuDist(1).^5;...
        0,1,2*CuDist(1),3*CuDist(1).^2,4*CuDist(1).^3,5.*CuDist(1).^4;...
        0,0,2,6.*CuDist(1),12.^CuDist(1).^2,20.*CuDist(1).^3;...
        1,CuDist(end),CuDist(end).^2,CuDist(end).^3,CuDist(end).^4,CuDist(end).^5;...
        0,1,2*CuDist(end),3*CuDist(end).^2,4*CuDist(end).^3,5.*CuDist(end).^4;...
        0,0,2,6.*CuDist(end),12.*CuDist(end).^2,20.*CuDist(end).^3];
    RVec = [F0;Ft0;Ftt0;FL;FtL;FttL];
    GCoX = AMat\RVec(:,1);
    GCoY = AMat\RVec(:,2);
    GCoZ = AMat\RVec(:,3);
    
    % Store the polynomial coefficients for output.
    GCoefs = [GCoX,GCoY,GCoZ];
    
    % Calculate the polynomial's values at each point along the curve.
    GXs = GCoX(1) + GCoX(2)*CuDist + GCoX(3)*CuDist.^2 + GCoX(4)*CuDist.^3 + GCoX(5)*CuDist.^4 + GCoX(6)*CuDist.^5;
    GYs = GCoY(1) + GCoY(2)*CuDist + GCoY(3)*CuDist.^2 + GCoY(4)*CuDist.^3 + GCoY(5)*CuDist.^4 + GCoY(6)*CuDist.^5;
    GZs = GCoZ(1) + GCoZ(2)*CuDist + GCoZ(3)*CuDist.^2 + GCoZ(4)*CuDist.^3 + GCoZ(5)*CuDist.^4 + GCoZ(6)*CuDist.^5;
    clear F0 FL Ft0 FtL Ftt0 FttL AMat RVec GCoX GCoY GCoZ
    
    % Define the final curves for Fourier analysis by subtracting the polynomial.
    Xs = B2(2:end-2,1)-GXs(1:end-1);
    Ys = B2(2:end-2,2)-GYs(1:end-1);
    if size(Bcoords,2) == 3
        Zs = B2(2:end-2,3)-GZs(1:end-1);
    end
else
    % --- For CLOSED curves: Simply remove the duplicate endpoint ---
    Xs = B2(1:end-1,1);
    Ys = B2(1:end-1,2);
    if size(Bcoords,2) == 3
        Zs = B2(1:end-1,3);
    end
    % No polynomial is needed for closed curves.
    GCoefs = NaN(6,3);
end

%% --- FOURIER COEFFICIENT CALCULATION --- %%

% Store the total length of the processed curve for output.
TotLen = CuDist(end);
% Remove the last point from CuDist to match the length of Xs, Ys, etc.
CuDist = CuDist(1:end-1);


% The Fourier coefficients are found by fitting the basis functions (sines
% and cosines) to the curve data via a system of linear equations.

% Construct the design matrix 'Mat' with basis functions as columns.
Mat = zeros(length(Xs),size(Xs,1));
for i=1:(size(Xs,1)-1)/2
    Mat(:,2*i) = sin(2*pi*i*CuDist/TotLen);
    Mat(:,2*i+1) = cos(2*pi*i*CuDist/TotLen);
end
Mat(:,1) = ones(length(Xs),1); % The n=0 (constant) term.

% Solve for the coefficients for each dimension using the matrix inverse.
CoefsX = inv(Mat)*Xs;
CoefsY = inv(Mat)*Ys;
if size(Bcoords,2) == 3
    CoefsZ = inv(Mat)*Zs;
    % Reshape the Z coefficients into the [sine, cosine] format.
    CoefsZ = [[0;CoefsZ(2:2:end)],CoefsZ(1:2:end)];
end

% Reshape the X and Y coefficients into the final [sine, cosine] format.
CoefsX = [[0;CoefsX(2:2:end)],CoefsX(1:2:end)];
CoefsY = [[0;CoefsY(2:2:end)],CoefsY(1:2:end)];

% If the original curve was 2D, return NaN for Z coefficients.
if size(Bcoords,2) ~= 3
    Zs = NaN;
    CoefsZ = ones(size(CoefsX)) * NaN;
end

end