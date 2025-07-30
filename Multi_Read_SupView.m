%==========================================================================
%                       FUNCTION Multi_Read_SupView
%==========================================================================
%
% PURPOSE:
%   Provides a user-friendly batch processing workflow for analyzing multiple
%   skull images for globularity from a superior view. It opens a file
%   dialog and displays a progress bar during analysis. For each skull, it
%   calls the 'SkulSupView' function and saves the resulting plot. It then
%   exports two sets of data: a detailed spreadsheet of point-by-point
%   curvature data for EACH specimen, and a single summary spreadsheet
%   containing the main globularity metrics for ALL specimens.
%
% OUTPUTS:
%   SkullResSup (struct):
%       A structure array containing the globularity analysis results for
%       each processed skull.
%   Files Generated:
%       - A summary Excel file: 'Results.xlsx'.
%       - A detailed data file for each specimen: '[Name]_PtData.xlsx'.
%       - An SVG image file for each specimen: '[Name]_result.svg'.
%       - An SVG image file without a title for each specimen: '[Name]_result_NT.svg'.
%
% DEPENDENCIES:
%   - SkulSupView.m: The main analysis function for a single skull.
%
% Author: Gadi Herzlinger and Uzy Smilansky (2023). All rights reserved.
%
% Associated with: 
% Mishol N., Herzlinger G., Rak Y., Smilansky U., Carmel L., Gokhman D. 
% (2025). Candidate Denisovan fossils identified through gene regulatory phenotyping.
%
%==========================================================================
function [SkullResSup] = Multi_Read_SupView()
%% --- FILE SELECTION --- %%

% Open a UI to let the user select one or more image files for analysis.
[FilesName, PathName] = uigetfile('C:\Users\Gadi Herzlinger\OneDrive - Haifa University\My Drive\Ongoing Projects\SkullGlobularity\*.*;','select file', 'MultiSelect', 'on'); 

% Handle cases for file selection.
% If the user cancels the dialog ('Done' is clicked with no selection), exit the function.
if ~ischar(FilesName) && ~iscell(FilesName) 
    clear
    return
% If only one file is selected, its name is returned as a char array.
% Convert it to a cell array to ensure it can be processed by the loop.
elseif ~iscell(FilesName) 
    FilesName = {FilesName}; 
end

%% --- BATCH PROCESSING LOOP --- %%

% Initialize a structure to hold the summary results for all skulls.
SupView = struct(); 
% Create and display a waitbar to provide feedback on the processing progress.
WaitB = waitbar(1/size(FilesName,2),'Processing...'); 

% Loop through each selected file.
for i =1:size(FilesName,2)
    % Extract the base name of the file (without extension) to use as an identifier.
    Name = string(FilesName{i}(1:end-4)); 
    % Store the name in the results structure.
    SupView(i).Name = Name; 
    
    % --- ANALYSIS & DATA STORAGE --- %
    % Call the main analysis function for the current file.
    [GlbPtsData,Kzr,GlbCrvSDN,GlbCrvSD, GlbCrvN,GlbCrv, Fig] = SkulSupView(PathName, FilesName{i}, Name); 
    
    % Assign the returned summary metrics to the results structure.
    SupView(i).GlbCrv = GlbCrv; 
    SupView(i).GlbCrvN = GlbCrvN; 
    SupView(i).GlbCrvSD = GlbCrvSD; 
    SupView(i).GlbCrvSDN = GlbCrvSDN; 
    SupView(i).Kzr = Kzr; 
    
    % --- SAVE FIGURES AND DATA --- %
    % Maximize the figure window and save it as a Scalable Vector Graphics (SVG) file.
    Fig.WindowState = 'maximized'; 
    saveas(Fig,[PathName,char(Name),'_result'],'svg'); 
    % Remove the title from the figure and save a second version.
    Fig.Children(2).Title.String = ''; 
    saveas(Fig,[PathName,char(Name),'_result_NT'],'svg'); 
    % Close the figure to avoid screen clutter during batch processing.
    delete(Fig); 
    
    % Update the waitbar to show the current progress.
    waitbar(i/size(FilesName,2),WaitB,['Processing... ',num2str(i),' / ',num2str(size(FilesName,2))]); 
    
    % For each specimen, save its detailed point-by-point data (Coords + Curvature) to a separate Excel file.
    writematrix(GlbPtsData,[PathName,char(Name),'_PtData.xlsx']); 
end

%% --- EXPORT FINAL SUMMARY --- %%

% Convert the final summary structure into a MATLAB table for easy export.
SkullResT = struct2table(SupView); 
% Write the summary table to a single Excel file in the same directory as the images.
writetable(SkullResT,[PathName,'Results.xlsx']); 
% Close and delete the waitbar window now that processing is complete.
delete(WaitB) 


end