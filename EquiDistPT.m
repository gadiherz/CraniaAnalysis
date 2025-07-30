%==========================================================================
%                         FUNCTION EquiDistPT
%==========================================================================
%
% PURPOSE:
%   Resamples an ordered set of 2D or 3D coordinates to produce a new set
%   of points that are spaced at equal distances along the original curve's
%   arc length. This function is a crucial dependency for 'FourierCoefs.m'.
%
% METHODOLOGY:
%   1. Calculates the total arc length of the input curve.
%   2. Determines an appropriate sampling distance and ensures the resulting
%      number of points is suitable for the subsequent Fourier analysis.
%   3. "Walks" along the segments of the original curve, placing new points
%      at each required distance interval through interpolation.
%   4. Outputs the new, equidistantly sampled points and their corresponding
%      cumulative distances along the curve.
%
% INPUTS:
%   Bcoords (matrix):
%       An N-by-2 or N-by-3 matrix of ordered coordinates along a curve.
%   OpenCurv (logical):
%       A flag indicating if the curve is open (true) or closed (false).
%
% OUTPUTS:
%   B2 (matrix):
%       A new M-by-2 or M-by-3 matrix of equidistantly spaced points.
%   CuDist (vector):
%       A vector of cumulative arc length values for each point in B2.
%
% AUTHORS:
%   Gadi Herzlinger and Uzy Smilansky.
%
%==========================================================================
function [B2, CuDist] = EquiDistPT(Bcoords, OpenCurv)

%% --- INITIAL CALCULATIONS --- %%

% Calculate the distance between each consecutive point in the original curve.
PDist = vecnorm(Bcoords(2:end,:)-Bcoords(1:end-1,:),2,2);
% Sum these distances to get the total arc length.
TDist = sum(PDist);

%% --- DETERMINE SAMPLING INTERVAL --- %%

% Calculate the initial required distance based on the average point-to-point distance.
RqDist = TDist/(size(Bcoords,1)-1);

% Enforce a minimum threshold for the sampling distance to prevent oversampling.
if RqDist < 7.5e-4
    RqDist = 7.5e-4;
end

% Determine the number of points for the new curve.
RqPts = fix(TDist./RqDist);
% Ensure the number of points is ODD. This results in an even number of
% intervals when the curve is closed (since the last point duplicates the
% first), which is a requirement for the 'FourierCoefs' function.
if rem(RqPts,2) == 0
    RqPts = RqPts-1;
end
% Recalculate the final, precise required distance based on the adjusted number of points.
RqDist = TDist/RqPts;

%% --- PLACE EQUIDISTANT POINTS --- %%

% This loop iterates along the original curve (B1) and places new,
% equidistantly spaced points into the output curve (B2).
B1=Bcoords;
B2 = B1(1,:); % Start the new curve with the first point of the old one.
n=2; % Start tracking from the second point of the original curve.

while size(B2,1) < RqPts
    % Calculate the straight-line distance from the last placed point in B2
    % to the current point being checked in B1.
    DistTrk = vecnorm(B1(n,:) - B2(end,:),2,2);
    
    % Check if the required distance (RqDist) falls beyond the current point.
    if sum(DistTrk) < RqDist
        % If so, keep advancing along the original curve until we "overshoot"
        % the target distance.
        while sum(DistTrk) < RqDist && n < size(B1,1)
            n=n+1;
            DistTrk = [DistTrk; vecnorm(B1(n-1,:) - B1(n,:),2,2)];
        end
        % Once overshoot occurs, interpolate backwards along the last segment
        % to place the new point at the exact required distance.
        DirVec = (B1(n,:) - B1(n-1,:))./vecnorm(B1(n,:) - B1(n-1,:),2,2);
        DistVec = RqDist-sum(DistTrk(1:end-1));
        B2(end+1,:) = B1(n-1,:) + DirVec*DistVec;
    else
        % If the required distance falls within the current segment,
        % simply place the new point along the vector towards the current B1 point.
        DirVec = (B1(n,:) - B2(end,:))./vecnorm(B1(n,:) - B2(end,:),2,2);
        B2(end+1,:) = B2(end,:) + DirVec*RqDist;
    end
end

%% --- FINALIZE OUTPUTS --- %%

% Add the final point to the curve.
if OpenCurv
    % For an open curve, the last point is the last point of the original curve.
    B2 = [B2;B1(end,:)];
else
    % For a closed curve, the last point is a duplicate of the first point.
    B2 = [B2;B2(1,:)];
end

% Create the vector of cumulative distances for the new set of points.
% Since points are equidistant, this is a simple cumulative sum of RqDist.
CuDist = [0; cumsum(repmat(RqDist,size(B2,1)-1,1))];


% The following code is commented out but retained for debugging purposes.
% It allows for visualization of the original and resampled curves.
% Fig1=figure();
% Ax1 = axes(Fig1);
% hold(Ax1,'on');
% BcSc = scatter3(Bcoords(:,1),Bcoords(:,2),Bcoords(:,3),'ob');
% BcPl = plot3(Bcoords(:,1),Bcoords(:,2),Bcoords(:,3),'b');
% B2Sc = scatter3(B2(:,1),B2(:,2),B2(:,3),'ok','filled');
% B2Pl = plot3(B2(:,1),B2(:,2),B2(:,3),'k');
end