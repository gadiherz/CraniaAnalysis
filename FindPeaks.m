%==========================================================================
%                           FUNCTION FindPeaks
%==========================================================================
%
% PURPOSE:
%   Identifies positive (maxima) and negative (minima) peaks in a data
%   vector. A point is considered a peak if its value is strictly greater
%   (or less) than all other points within a defined neighborhood. This
%   function is used to find key anatomical landmarks based on curvature.
%
% METHODOLOGY:
%   The function vectorizes the peak-finding process by creating matrices
%   of indices for the neighborhood of each point. It then performs a
%   single logical comparison to find all points that satisfy the peak
%   condition. This is more efficient than looping through each point.
%
% INPUTS:
%   Y (vector):
%       The data values in which to find peaks (e.g., curvature).
%   X (vector):
%       The corresponding position values for each data point (e.g., arc length).
%   Closed (logical):
%       Flag indicating if the data is from a closed, periodic curve. If
%       true, the data is "wrapped" to correctly find peaks near the endpoints.
%   Abs (logical):
%       Flag to enforce absolute peak definition. If true, positive peaks
%       must have a value > 0, and negative peaks must have a value < 0.
%   Neigh (integer):
%       The size of the neighborhood (number of points on each side) to
%       check against when defining a peak.
%
% OUTPUTS:
%   PIndPos, PIndNeg (matrix):
%       N-by-3 matrices containing peak data. Column 1 is the point's index,
%       column 2 is its position (from X), and column 3 is its value (from Y).
%
% AUTHORS:
%   Gadi Herzlinger and Uzy Smilansky.
%
%==========================================================================
function [PIndPos,PIndNeg] = FindPeaks(Y, X,Closed,Abs,Neigh)

%% --- HANDLE CLOSED CURVE WRAP-AROUND --- %%

% If the curve is closed, pad the Y vector by wrapping points from the
% opposite end. This ensures that points near the start/end have a full,
% periodic neighborhood to be compared against.
if Closed
    if abs(Y(1) - Y(end)) > abs(Y(end).*10e4)
        Y = [Y;Y(1)];
    end
    Y = [Y(end-Neigh:end-1);Y;Y(2:Neigh+1)];
end

%% --- NEIGHBORHOOD COMPARISON --- %%

% Create index matrices to get the values of all neighbors for all points
% in a single, vectorized operation.
% 'indB' holds the indices for the neighborhood BEFORE each point.
indB = repmat(1:Neigh,size(Y,1)-2*Neigh,1);
indB = indB + repmat((1:size(indB,1))'-1,1,Neigh);
Bef = Y(indB);

% 'indA' holds the indices for the neighborhood AFTER each point.
indA = repmat(Neigh+2:2*Neigh+1,size(Y,1)-2*Neigh,1);
indA = indA + repmat((1:size(indA,1))'-1,1,Neigh);
Aft = Y(indA);

% Perform logical comparisons to find the peaks.
% A positive peak is a point greater than ALL points in its 'before' and 'after' neighborhoods.
PeaksP = all(Y(Neigh+1:end-Neigh)>Bef,2) & all(Y(Neigh+1:end-Neigh)>Aft,2);
% A negative peak is a point less than ALL points in its 'before' and 'after' neighborhoods.
PeaksN = all(Y(Neigh+1:end-Neigh)<Bef,2) & all(Y(Neigh+1:end-Neigh)<Aft,2);

% For an open curve, the first and last 'Neigh' points cannot be peaks by
% definition, as they lack a full neighborhood. This pads the boolean
% vectors to match the original 'Y' vector size.
if ~Closed
    PeaksP = [false(Neigh,1);PeaksP;false(Neigh,1)];
    PeaksN = [false(Neigh,1);PeaksN;false(Neigh,1)];
else
    % If closed, trim the padded Y vector back to its original size.
    Y = Y(Neigh+1:end-Neigh);
end

% Get the indices of the points that are peaks.
PeaksP = find(PeaksP);
PeaksN = find(PeaksN);


%% --- FORMAT OUTPUT --- %%

% If requested, filter out non-absolute peaks (e.g., a "positive" peak
% that has a negative value).
if Abs
    PeaksP = PeaksP(Y(PeaksP)>0);
    PeaksN = PeaksN(Y(PeaksN)<0);
end

% Combine the peak indices with their corresponding position (X) and
% value (Y) into the final output matrices.
PIndPos = [PeaksP,X(PeaksP),Y(PeaksP)];
PIndNeg = [PeaksN,X(PeaksN),Y(PeaksN)];

end