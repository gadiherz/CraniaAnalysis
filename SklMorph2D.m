%==========================================================================
%                           FUNCTION SklMorph2D
%==========================================================================
%
% PURPOSE:
%   Analyzes a 2D right-lateral image of a hominin skull to compute metrics
%   for skull top flatness and forehead projection. The workflow involves
%   Fourier analysis of the skull outline, user-guided selection of a key
%   anatomical landmark, and a parametric curve fitting of the skull vault
%   to derive robust measurements.
%
% METHODOLOGY:
%   1. Extracts the skull outline from an image and computes its Fourier representation.
%   2. Launches the 'PtSelector' app for the user to identify the deepest
%      point of the concavity superior to the supraorbital ridge.
%   3. Identifies the forehead and skull vertex points based on the tangent
%      of the curve, providing a robust, data-driven landmarking approach.
%   4. Fits a parametric function to the skull vault for a smooth representation.
%   5. Calculates "skull flatness" by integrating the curvature of this fitted
%      function and measures forehead height/angle relative to the defined baseline.
%
% INPUTS:
%   PathName (string):
%       The file path to the directory containing the image file.
%   FilesName (string):
%       The name of the image file to be analyzed (e.g., 'specimen.jpg').
%   Name (string):
%       An identifier or name for the specimen, used for plot titles.
%
% OUTPUTS:
%   CrvLnTp (double): Length of the curve along the skull top.
%   SkFltns (double): Skull top flatness (integral of absolute curvature).
%   FitErr (double): Root Mean Square Error of the fitted parametric model.
%   FhCY... (double): Curvature-based forehead metrics (largely deprecated).
%   FhTYDiff (double): Tangent-based forehead height above the baseline.
%   FhTYRt (double): Ratio of tangent-based forehead height to skull vertex height.
%   FhTYTrRt (double): Tangent of the angle of the vector from glabella to forehead.
%   Fig1 (handle): The handle to the output figure visualizing the analysis.
%
% DEPENDENCIES:
%   - PtSelector.m: A MATLAB App for interactive point selection.
%   - FourierCoefs.m, Sum4Fourier.m, TanCalc.m: Helper functions for
%     Fourier analysis and geometry calculations.
%
% Author: Gadi Herzlinger and Uzy Smilansky (2023). 
%
% Associated with: 
% Mishol N., Herzlinger G., Rak Y., Smilansky U., Carmel L., Gokhman D. 
% (2025). Candidate Denisovan fossils identified through gene regulatory phenotyping.
%==========================================================================
function [CrvLnTp, SkFltns, FitErr, FhCYDiff, FhCYRt, FhCYTrRt, FhTYDiff, FhTYRt, FhTYTrRt,Fig1] = SklMorph2D(PathName, FilesName, Name)

%% --- IMAGE PRE-PROCESSING AND FOURIER ANALYSIS --- %%

% Read the specified image file.
RGBread = imread([PathName,FilesName]);
% Isolate non-white pixels. This assumes a dark outline on a white background.
PtMat = sum(double(RGBread) - repmat(255,size(RGBread)),3);
[r,c] = find(PtMat < -20);
Mat1=nan(size(RGBread));
for i=1:(length(r))
    Mat1(r(i),c(i),:) = RGBread(r(i),c(i),:);
end

% Convert pixel indices to a standard Cartesian coordinate system.
Y = -r+abs(min(-r));
X = c;
clear r c;

% Center the coordinates by subtracting the mean (centroid).
Orig = mean([X,Y]);
X = X - Orig(1);
Y = Y - Orig(2);
% Extract the boundary points of the skull shape.
SklPt = boundary([X,Y],0.9);
SklBn = [X(SklPt),Y(SklPt)];
% Compute the Fourier coefficients of the boundary shape.
[TotLen, CoefsX, CoefsY, ~, ~] = FourierCoefs(SklBn, false);

%% --- INTERACTIVE POINT SELECTION --- %%

% Launch the interactive point selection app (PtSelector).
% This app allows the user to select the deepest point of the concavity
% above the supraorbital ridge, which defines the anterior landmark.
Ap = PtSelector(X,Y,CoefsX,CoefsY,TotLen,Name);

% Pause execution and wait for the user to finalize selection in the app.
while ~Ap.Done
    pause(0.05);
end

%% --- DATA RETRIEVAL FROM APP --- %%

% Retrieve the processed data from the app. Data is scaled by Centroid Size.
% The alternative scaling by Glabella-Opisthocranion Length (GOL) is
% calculated in the app but not used in this analysis.
PtInd = Ap.Cnt.PtsInd(:,1);     % Anterior/posterior landmark positions (arc length)
TotLen = Ap.Cnt.TotLen;         % Total length of the scaled curve
CoefsX = Ap.Cnt.CoefsX;         % Scaled Fourier coefficients for X
CoefsY = Ap.Cnt.CoefsY;         % Scaled Fourier coefficients for Y
FhPt = Ap.Cnt.FhPt;             % Curvature-based forehead point (not used)
OpiT = Ap.Cnt.OpiT;             % Opisthocranion point (not used)
CntSz = Ap.Cnt.Size;            % The centroid size scaling factor

% Retrieve GOL-scaled data (recorded but unused).
PtIndG = Ap.GOL.PtsInd(:,1);
TotLenG = Ap.GOL.TotLen;
CoefsXG = Ap.GOL.CoefsX;
CoefsYG = Ap.GOL.CoefsY;
FhPtG = Ap.GOL.FhPt;
OpiTG = Ap.GOL.OpiT;
GOLSz = Ap.GOL.Size;

% Retrieve curvature-based forehead metrics (deprecated method, but recorded).
FhCYDiffNS = Ap.FHYDiff;        % Unscaled height difference
FhCYRt = Ap.FHYRt;              % Height ratio
FhCYTrRt = Ap.FHTrRt;           % Tangent ratio

% Close and delete the app object to free up resources.
Ap.CloseWindow;
delete(Ap);

%% --- TANGENT-BASED FOREHEAD & VERTEX IDENTIFICATION --- %%

% Define the curve segment for the top of the skull (from anterior to posterior landmarks).
CrvLnTp = (PtInd(2)-PtInd(1));
TpSkCrv = (PtInd(1):(PtInd(2)-PtInd(1))/299:PtInd(2))';
% Reconstruct the coordinates of this top curve segment.
for j=1:size(TpSkCrv)
    XpredT(j,1) = Sum4Fourier(TotLen,CoefsX,TpSkCrv(j),ones(size(CoefsX,1)));
    YpredT(j,1) = Sum4Fourier(TotLen,CoefsY,TpSkCrv(j),ones(size(CoefsY,1)));
end

% Define a smoothing filter and calculate the tangent along the skull top.
WnF = 1./(1+exp(1+(((1:1:size(CoefsX))-14)/1)))';
FunTanF = TanCalc(TotLen,TpSkCrv,WnF,CoefsX,CoefsY);

% Identify the forehead point by finding where the tangent angle approaches 45 degrees (pi/4).
% This logic robustly finds the most likely point in the anterior portion of the curve.
PotFh = find(islocalmin(abs(abs(FunTanF(1:150))-0.25*pi)));
if numel(PotFh)>1
    % If multiple candidates exist, a cost function selects the best one.
    PotFhCstFun = abs(abs(FunTanF(PotFh))-0.25*pi).*abs(75-PotFh);
    FhT2 = PotFh(PotFhCstFun==min(PotFhCstFun));
elseif numel(PotFh)==1
    FhT2 = PotFh;
else
    % If no local minimum is found, take the point of closest approach.
    FhT2 = find(abs(abs(FunTanF(1:150))-0.25*pi)==min(abs(abs(FunTanF(1:150))-0.25*pi)),1);
end

% Identify the skull vertex by finding where the tangent angle is closest to 0.
SkVr = find(abs(FunTanF(100:200))==min(abs(FunTanF(100:200))),1);
% Get the precise tangent value and coordinates for the identified forehead and vertex points.
FhTVal = FunTanF(FhT2);
FhPt = [Sum4Fourier(TotLen,CoefsX,TpSkCrv(FhT2),WnF),Sum4Fourier(TotLen,CoefsY,TpSkCrv(FhT2),WnF)];
SkVPt = [Sum4Fourier(TotLen,CoefsX,TpSkCrv(SkVr),WnF),Sum4Fourier(TotLen,CoefsY,TpSkCrv(SkVr),WnF)];

%% --- PARAMETRIC FITTING AND FLATNESS CALCULATION --- %%

% Re-parameterize the top-of-skull curve by its normalized arc length.
PDists = vecnorm([XpredT(2:end),YpredT(2:end)]-[XpredT(1:end-1),YpredT(1:end-1)],2,2);
ArcLen = [0;cumsum(PDists)];
ArcTotLen = ArcLen(end);
ArcLen = ArcLen./ArcTotLen;

% Fit the curve to a parametric function using least squares.
% This models the shape with a small number of trigonometric terms.
AMatC = [cos(ArcLen), sin(ArcLen), cos(ArcLen).^3, sin(ArcLen).^3];
[XcofsC, ~] = lsqr(AMatC,XpredT);
[YcofsC, ~] = lsqr(AMatC,YpredT);


% Define anonymous functions for the fitted curve and its derivatives.
Xfun = @(t) (XcofsC(1).*cos(t) + XcofsC(2).*sin(t) + XcofsC(3).*cos(t).^3 + XcofsC(4).*sin(t).^3);
Yfun = @(t) (YcofsC(1).*cos(t) + YcofsC(2).*sin(t) + YcofsC(3).*cos(t).^3 + YcofsC(4).*sin(t).^3);
% 1st derivatives (for tangent).
dxdt = @(t) 3.*XcofsC(4).*cos(t).*sin(t).^2 + (-3.*XcofsC(3).*cos(t).^2 - XcofsC(1)).*sin(t) + XcofsC(2).*cos(t);
dydt = @(t) 3.*YcofsC(4).*cos(t).*sin(t).^2 + (-3.*YcofsC(3).*cos(t).^2 - YcofsC(1)).*sin(t) + YcofsC(2).*cos(t);
% 2nd derivatives (for curvature).
d2xdt2 = @(t) -3.*XcofsC(4).*sin(t).^3 + 6.*XcofsC(3).*cos(t).*sin(t).^2 + (6.*XcofsC(4).*cos(t).^2 - XcofsC(2)).*sin(t) - 3.*XcofsC(3).*cos(t).^3 - XcofsC(1).*cos(t);
d2ydt2 = @(t) -3.*YcofsC(4).*sin(t).^3 + 6.*YcofsC(3).*cos(t).*sin(t).^2 + (6.*YcofsC(4).*cos(t).^2 - YcofsC(2)).*sin(t) - 3.*YcofsC(3).*cos(t).^3 - YcofsC(1).*cos(t);

% Define tangent and curvature functions based on the parametric model.
FunTan = @(t) acos(dxdt(t)./sqrt(dxdt(t).^2+dydt(t).^2))-pi;
FunCurv = @(t) -(dydt(t).*d2xdt2(t) - dxdt(t).*d2ydt2(t))./(dxdt(t).^2 + dydt(t).^2);

% Calculate skull flatness by integrating the absolute curvature over the fitted vault.
SkFltns = integral(@(t) abs(FunCurv(t)), ArcLen(1), ArcLen(end));

%% --- FINAL METRIC CALCULATION --- %%

% Calculate final forehead height metrics from the tangent-based method.
FhCYDiff = FhCYDiffNS./CntSz; % Scale the (deprecated) curvature-based height.
FhTYDiff = FhPt(2) - Sum4Fourier(TotLen,CoefsY,TpSkCrv(1),ones(size(CoefsX,1))); % Tangent-based height difference.
FhTYRt = FhTYDiff./ (SkVPt(2) - Sum4Fourier(TotLen,CoefsY,TpSkCrv(1),ones(size(CoefsX,1)))); % Ratio of forehead height to vertex height.
% Vector from anterior landmark to the vertex, and its angle (as a tangent).
FhVec = [SkVPt(1),SkVPt(2)] - [Sum4Fourier(TotLen,CoefsX,TpSkCrv(1),ones(size(CoefsX,1))),Sum4Fourier(TotLen,CoefsY,TpSkCrv(1),ones(size(CoefsX,1)))];
FhVec = FhVec./vecnorm(FhVec,2,2);
FhTYTrRt = tan(acos(dot(FhVec,[-1,0],2)));

% Calculate the fit error of the parametric model as a Root Mean Square Error.
FitErr = sqrt(mean(vecnorm([XpredT,YpredT] - [Xfun(ArcLen),Yfun(ArcLen)],2,2).^2));

%% --- VISUALIZATION --- %%

% Prepare data for plotting.
CuDist = (0:TotLen/299:TotLen)';
CurvL = FunCurv(ArcLen);
% Reconstruct the full skull outline with a smooth filter for graphics.
WnGrph = 1./(1+exp((1:1:size(CoefsX))-fix(0.8*size(CoefsX,1)))/fix(0.2*size(CoefsX,1)))';
for i=1:size(CuDist)
    XpredG(i,1) = Sum4Fourier(TotLen,CoefsX,CuDist(i),WnGrph);
    YpredG(i,1) = Sum4Fourier(TotLen,CoefsY,CuDist(i),WnGrph);
end
% Get coordinates for the anterior and posterior landmarks.
PtC = [[Sum4Fourier(TotLen,CoefsX,PtInd(1),WnGrph),Sum4Fourier(TotLen,CoefsY,PtInd(1),WnGrph)];...
    [Sum4Fourier(TotLen,CoefsX,PtInd(2),WnGrph),Sum4Fourier(TotLen,CoefsY,PtInd(2),WnGrph)]];
% Get coordinates for plotting the tangent line at the forehead point.
TanPl2 = [Sum4Fourier(TotLen,CoefsX,TpSkCrv(FhT2),WnF) + cos(FhTVal),...
    Sum4Fourier(TotLen,CoefsY,TpSkCrv(FhT2),WnF) + sin(FhTVal);...
    Sum4Fourier(TotLen,CoefsX,TpSkCrv(FhT2),WnF) - cos(FhTVal),...
    Sum4Fourier(TotLen,CoefsY,TpSkCrv(FhT2),WnF) - sin(FhTVal)];

% Create figure and axes.
Fig1=figure();
Fig1.NumberTitle = 'off';
Fig1.Name = Name;
Ax = axes(Fig1);
hold(Ax,"on");
axis(Ax,'equal');

% Plot the various components of the analysis.
HResPl = plot(Ax,XpredG,YpredG,'k'); % Full skull outline
PtsSc = scatter(Ax,PtC(:,1),PtC(:,2),120,'pr','filled'); % Base landmarks
CSLPl = plot(Ax,PtC(:,1),PtC(:,2),'r','LineWidth',2); % Baseline
TpCrv = scatter(Ax,Xfun(ArcLen),Yfun(ArcLen),45,CurvL,'o','filled'); % Vault curve colored by curvature
FhTPtSc = scatter(Ax,FhPt(1),FhPt(2),150,[19,228,228]./256,'d','filled'); % Forehead landmark
FhTPl = plot(TanPl2(:,1),TanPl2(:,2),'k','LineWidth',2); % Tangent line at forehead
FHTAngPl = plot([PtC(1,1);FhPt(1)],[PtC(1,2);FhPt(2)],'Color',[19,228,228]./256,'LineWidth',2); % Forehead angle line
FHYTDiffPl = plot([FhPt(1);FhPt(1)],[FhPt(2);mean([PtC(1,2);PtC(2,2)])],'Color',[19,228,228]./256,'LineWidth',1.5,'LineStyle','--'); % Forehead height line

% Finalize plot details.
clim(Ax,[0,4]);
colorbar(Ax);
title(Ax,Name,'Interpreter','none');
legend([PtsSc,FhTPtSc,TpCrv,FhTPl,FHTAngPl],{'Supraorbital ridge concavity - posterior reflection','Forehead Landmark',...
      'Skull top flatness','Forehead tangent angle','Forehead angle and height'},"Location","southwest");
% Reorder the plotted elements to ensure correct layering.
Child = Ax.Children;
Ax.Children = [Child(3),Child(2),Child(7),Child(4),Child(8),Child(6),Child(5),Child(1)];

end