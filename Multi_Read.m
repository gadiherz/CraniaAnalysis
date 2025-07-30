%==========================================================================
%                           FUNCTION Multi_Read
%==========================================================================
%
% PURPOSE:
%   Provides a batch processing workflow for the 'SklMorph2D' analysis.
%   It opens a file dialog for the user to select one or more image files.
%   It then iterates through each file, calling the main analysis function,
%   and aggregates all results. Finally, it saves the output figures and
%   compiles the numerical data into a single Excel spreadsheet.
%
% OUTPUTS:
%   SkullRes (struct):
%       A structure array containing the full set of analysis results for
%       each processed skull. This is also written to 'Results.xlsx'.
%   Files Generated:
%       - A summary Excel file: 'Results.xlsx'.
%       - An SVG image file for each specimen: '[Name]_Result.svg'.
%       - An SVG image file without a title/legend for each specimen:
%         '[Name]_Result_NT.svg'.
%
% DEPENDENCIES:
%   - SklMorph2D.m: The main analysis function for a single skull.
%
% Author: Gadi Herzlinger and Uzy Smilansky (2023). 
%
% Associated with: 
% Mishol N., Herzlinger G., Rak Y., Smilansky U., Carmel L., Gokhman D. 
% (2025). Candidate Denisovan fossils identified through gene regulatory phenotyping.
%==========================================================================
function [SkullRes] = Multi_Read()

%% --- FILE SELECTION --- %%

% Open a UI to let the user select multiple image files.
[FilesName, PathName] = uigetfile('G:\My Drive\Ongoing Projects\SkullGlobularity\*.*;','select file', 'MultiSelect', 'on');

% Handle cases for file selection.
% If the user cancels the dialog, exit the function.
if ~ischar(FilesName) && ~iscell(FilesName)
    clear
    return
% If only one file is selected, convert its name to a cell array for consistent looping.
elseif ~iscell(FilesName)
    FilesName = {FilesName};
end

%% --- BATCH PROCESSING LOOP --- %%

% Initialize a structure to hold the results for all skulls.
SkullRes = struct();

% Loop through each selected file.
for i =1:size(FilesName,2)
    
    % Extract the base name of the file (without extension) to use as an identifier.
    Name = string(FilesName{i}(1:end-4));
    
    % Store the name in the results structure.
    SkullRes(i).Name = Name;
    
    % Call the main analysis function for the current file.
    [CrvLnTp, SkFltns, FitErr, FhCYDiff, FhCYRt, FhCYTrRt, FhTYDiff, FhTYRt, FhTYTrRt, Fig] = SklMorph2D(PathName, FilesName{i}, Name);
    
    % Assign all returned metrics from the analysis to the results structure.
    SkullRes(i).CrvLnTpCS = CrvLnTp;
    SkullRes(i).SkFltnsCS = SkFltns;
    SkullRes(i).FitErr = FitErr;
    SkullRes(i).FhCYDiff = FhCYDiff;
    SkullRes(i).FhCYRt = FhCYRt;
    SkullRes(i).FhCYTrRt = FhCYTrRt;
    SkullRes(i).FhTYDiff = FhTYDiff;
    SkullRes(i).FhTYRt = FhTYRt;
    SkullRes(i).FhTYTrRt = FhTYTrRt;
    
    
    % --- SAVE FIGURE --- %
    % Maximize the figure window for better visibility.
    Fig.WindowState = 'maximized';
    % Save the figure as a Scalable Vector Graphics (SVG) file.
    saveas(Fig,[PathName,char(Name),'_Result'],'svg');
    % Create a second version of the figure with the legend and title removed for cleaner presentation.
    Fig.Children(1).Visible = "off";
    Fig.Children(3).Title.Visible = "off";
    saveas(Fig,[PathName,char(Name),'_Result_NT'],'svg');
    % Close the figure to avoid screen clutter during batch processing.
    % Note: delete is called twice, which is redundant but harmless.
    delete(Fig);
    delete(Fig);
end

%% --- EXPORT FINAL RESULTS --- %%

% Convert the results structure into a MATLAB table for easy export.
SkullResT = struct2table(SkullRes);
% Write the table to an Excel file in the same directory as the images.
writetable(SkullResT,[PathName,'Results.xlsx']);
end