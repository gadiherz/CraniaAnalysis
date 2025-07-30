%==========================================================================
%                           FUNCTION SkulSupView
%==========================================================================
%
% PURPOSE:
%   Analyzes the 2D geometry of a hominin skull from a superior view image
%   to quantify the globularity of the frontal bone. It automatically
%   identifies the glabella region (the most anterior point on the midline),
%   calculates curvature along this segment, and derives metrics of total
%   curvature and curvature variance.
%
% INPUTS:
%   PathName (string):
%       The file path to the directory containing the image file.
%   FilesName (string):
%       The name of the image file to be analyzed.
%   Name (string):
%       An identifier or name for the specimen being analyzed.
%
% OUTPUTS:
%   GlbPtsData (matrix):
%       An N-by-3 matrix containing the X,Y coordinates and curvature (K)
%       for each of the 20 points analyzed in the glabella region.
%   Kzr (double):
%       The specific curvature value at the glabella zero-crossing point.
%   GlbCrvSDN (double):
%       The standard deviation of curvature, normalized by arc length.
%   GlbCrvSD (double):
%       The standard deviation of curvature values in the glabella region.
%       This is a key metric for globularity (lower SD = more globular).
%   GlbCrvN (double):
%       The total curvature, normalized by arc length.
%   GlbCrv (double):
%       The sum of curvature values across the glabella region.
%   Fig1 (handle):
%       The handle to the output figure visualizing the analysis.
%
% DEPENDENCIES:
%   - FourierCoefs.m: Computes Fourier coefficients from 2D coordinates.
%   - Sum4Fourier.m: Reconstructs coordinates from Fourier coefficients.
%   - CurvCalc.m: Calculates curvature from Fourier coefficients.
%
% Author: Gadi Herzlinger and Uzy Smilansky (2023). 
%
% Associated with: 
% Mishol N., Herzlinger G., Rak Y., Smilansky U., Carmel L., Gokhman D. 
% (2025). Candidate Denisovan fossils identified through gene regulatory phenotyping.
%==========================================================================
function [GlbPtsData,Kzr,GlbCrvSDN,GlbCrvSD,GlbCrvN,GlbCrv,Fig1] = SkulSupView(PathName, FilesName, Name)

%% --- IMAGE PRE-PROCESSING AND POINT EXTRACTION --- %%

% Read the specified image file.
RGBread = imread([PathName,FilesName]);
% Isolate non-white pixels, assuming a dark outline on a white background.
PtMat = sum(double(RGBread) - repmat(255,size(RGBread)),3);
[r,c] = find(PtMat < -20);
% Create a matrix containing only the RGB values of the detected pixels2].
Mat1=nan(size(RGBread));
for i=1:(length(r))
    Mat1(r(i),c(i),:) = RGBread(r(i),c(i),:);
end
AMap = ones(size(Mat1(:,:,1)));
AMap(~isnan(Mat1(:,:,1))) = 0;

% Convert pixel indices to a standard Cartesian coordinate system.
Y = -r+abs(min(-r));
X = c;
clear r c;
% Center the coordinates by subtracting the mean (centroid).
Orig = mean([X,Y]);
X = X - Orig(1);
Y = Y - Orig(2);

% Extract the boundary points and compute their Fourier coefficients.
SklPt = boundary([X,Y],0.9);
SklBn = [X(SklPt),Y(SklPt)];
[TotLen, CoefsX, CoefsY, ~, ~] = FourierCoefs(SklBn, false);

%% --- AUTOMATIC IDENTIFICATION OF THE GLABELLA REGION --- %%

% Divide the entire outline into 49 small, overlapping segments.
StrtPts = linspace(0,TotLen,50);
StrtPts = [StrtPts(1:end-1)',StrtPts(2:end)'];
% Identify all segments that cross the horizontal midline (Y=0).
PotPtInds = nan(49,2);
for i=1:49
    PotPtInds(i,:) = [Sum4Fourier(TotLen,CoefsY,StrtPts(i,1),ones(size(CoefsY,1))),Sum4Fourier(TotLen,CoefsY,StrtPts(i,2),ones(size(CoefsY,1)))];
end
PotIntvl = StrtPts(find(sign(PotPtInds(:,1)).*sign(PotPtInds(:,2))<0),:);

% Of the midline-crossing segments, find the one that is most anterior (max X value).
% This segment contains the glabella.
for i=1:size(PotIntvl,1)
    Xval(i) = Sum4Fourier(TotLen,CoefsX,PotIntvl(i,1),ones(size(CoefsX,1)));
end
Invl = PotIntvl(Xval == max(Xval),:);

% Find the precise arc-length point ('ZrPt') where Y=0 within this segment.
ZrPt = fzero(@(x) Sum4Fourier(TotLen,CoefsY,x,ones(size(CoefsY,1))),Invl);
% Define a new, smaller interval centered on this glabella point for analysis.
Invl = [ZrPt-(TotLen/25),ZrPt+(TotLen/25)];

%% --- CURVATURE ANALYSIS OF THE GLABELLA REGION --- %%

% Sample 20 points along the glabella interval.
Invl = linspace(Invl(1),Invl(2),20)';
% Define a smoothing window (low-pass filter) for the Fourier coefficients.
WnK = 1./(1+exp((1:1:size(CoefsX))-20)/3)';

% Reconstruct the smoothed XY coordinates for the 20 sample points.
for i =1:20
    InvlCrd(i,:) = [Sum4Fourier(TotLen,CoefsX,Invl(i),WnK),...
        Sum4Fourier(TotLen,CoefsY,Invl(i),WnK)];
end
% Reconstruct the coordinates for the precise glabella zero-crossing point.
ZrCrds = [Sum4Fourier(TotLen,CoefsX,ZrPt,WnK),...
        Sum4Fourier(TotLen,CoefsY,ZrPt,WnK)];
        
% Calculate the curvature (K) at each of the 20 points.
K = CurvCalc(TotLen,Invl,WnK,CoefsX,CoefsY);
% Store the coordinates and their corresponding curvature values.
GlbPtsData = [InvlCrd,K];
% Calculate the specific curvature at the glabella point.
Kzr = CurvCalc(TotLen,ZrPt,WnK,CoefsX,CoefsY);

% Calculate globularity metrics.
GlbCrv = sum(K); % Total curvature in the region
GlbCrvN = GlbCrv./(Invl(end)-Invl(1)); % Total curvature normalized by arc length
GlbCrvSD = std(K); % Standard deviation of curvature (a measure of globularity)
GlbCrvSDN = GlbCrvSD/(Invl(end)-Invl(1)); % Normalized standard deviation of curvature

%% --- VISUALIZATION --- %%

% Generate a high-resolution version of the full skull outline for plotting
CuDist = (0:TotLen/299:TotLen)';
WnGrph = 1./(1+exp((1:1:size(CoefsX))-20)/3)';
for i=1:size(CuDist)
    XpredG(i,1) = Sum4Fourier(TotLen,CoefsX,CuDist(i),WnGrph);
    YpredG(i,1) = Sum4Fourier(TotLen,CoefsY,CuDist(i),WnGrph);
end

% Create a new figure to display the results.
Fig1=figure();
Fig1.NumberTitle = 'off';
Fig1.Name = Name;
Ax = axes(Fig1);
hold(Ax,"on");
axis(Ax,'equal');
% Plot the full, high-resolution skull outline
HResPl = plot(Ax,XpredG,YpredG,'k');
% Highlight the analyzed glabella interval with a thick line
InvlPl = plot(Ax,InvlCrd(:,1),InvlCrd(:,2),'k','LineWidth',3);
% Use a scatter plot to color the points on the interval by their curvature value
InvlSC = scatter(Ax,InvlCrd(:,1),InvlCrd(:,2),25,K,'o','filled');
% Mark the precise glabella point with a large, color-coded pentagram
KzrPl = scatter(Ax,ZrCrds(1),ZrCrds(2),100,Kzr,'p','filled','MarkerEdgeColor','k');
% Set the plot title and add a color bar for the curvature scale
title(Ax,Name,'Interpreter','none');
clim(Ax,[-5e-3,5e-3]);
colorbar(Ax);
% Reorder the plotted elements to ensure correct layering
Child = Ax.Children;
Ax.Children = [Child(4),Child(1),Child(2),Child(3)];

end