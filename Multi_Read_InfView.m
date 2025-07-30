%==========================================================================
%                           FUNCTION Multi_Read_InfView
%==========================================================================
%
% PURPOSE:
%   Provides a batch processing workflow for analyzing multiple skull images.
%   It opens a file dialog for the user to select one or more image files.
%   It then iterates through each file, calling the main 'SkulInfView'
%   function to perform the geometric analysis. Finally, it aggregates all
%   results, saves the output figures, and compiles the numerical data into
%   an Excel spreadsheet.
%
% OUTPUTS:
%   SkullResInf (struct):
%       A structure array containing the full set of analysis results for
%       each processed skull. This is also written to an Excel file.
%
% PROCESS:
%   1. Prompts user to select image files via a 'uigetfile' dialog.
%   2. Initializes a structure to store results.
%   3. Loops through each selected file:
%      a. Calls 'SkulInfView' to analyze the image.
%      b. Stores the returned measurements in the results structure.
%      c. Saves the generated result figure in SVG format.
%   4. Converts the final results structure to a MATLAB table.
%   5. Writes the table to an 'Results.xlsx' file in the image directory.
%
% DEPENDENCIES:
%   - SkulInfView.m: The main analysis function for a single skull.
%
% Author: Gadi Herzlinger and Uzy Smilansky (2023). All rights reserved.
%
% Associated with: 
% Mishol N., Herzlinger G., Rak Y., Smilansky U., Carmel L., Gokhman D. 
% (2025). Candidate Denisovan fossils identified through gene regulatory phenotyping.
%==========================================================================
function [SkullResInf] = Multi_Read_InfView() 
%% --- FILE SELECTION --- %%

% Open a UI to let the user select multiple image files.
[FilesName, PathName] = uigetfile('C:\Users\Gadi Herzlinger\OneDrive - Haifa University\My Drive\Ongoing Projects\SkullGlobularity\*.*;','select file', 'MultiSelect', 'on'); 

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
    [CntSz,CrvnsL,CrvnsLNrm,CurvLPPos,CurvLPAngDir,CurvLPAngFit,CurvLPAngFit1,CurvLPAngFit2,...
        FitErrL,CrvnsR,CrvnsRNrm,CurvRPPos,CurvRPAngDir,CurvRPAngFit,CurvRPAngFit1,CurvRPAngFit2,...
        FitErrR,Fig] = SkulInfView(PathName, FilesName{i}, Name); 
        
    % --- STORE RESULTS --- %
    % Assign all returned metrics from the analysis to the results structure.
    SkullRes(i).CntSz = CntSz; 
    SkullRes(i).CrvnsL = CrvnsL; 
    SkullRes(i).CrvnsLNrm = CrvnsLNrm; 
    SkullRes(i).CurvLPPos = CurvLPPos; 
    SkullRes(i).CurvLPAngDir = CurvLPAngDir; 
    SkullRes(i).CurvLPAngFit = CurvLPAngFit; 
    SkullRes(i).CurvLPAngFit1 = CurvLPAngFit1; 
    SkullRes(i).CurvLPAngFit2 = CurvLPAngFit2; 
    SkullRes(i).FitErrL = FitErrL; 
    SkullRes(i).CrvnsR = CrvnsR; 
    SkullRes(i).CrvnsRNrm = CrvnsRNrm; 
    SkullRes(i).CurvRPPos = CurvRPPos; 
    SkullRes(i).CurvRPAngDir = CurvRPAngDir; 
    SkullRes(i).CurvRPAngFit = CurvRPAngFit; 
    SkullRes(i).CurvRPAngFit1 = CurvRPAngFit1; 
    SkullRes(i).CurvRPAngFit2 = CurvRPAngFit2; 
    SkullRes(i).FitErrR = FitErrR; 
    
    % --- SAVE FIGURE --- %
    % Maximize the figure window for better visibility.
    Fig.WindowState = 'maximized'; 
    % Save the figure as a Scalable Vector Graphics (SVG) file.
    saveas(Fig,[PathName,char(Name),'_Result'],'svg'); 
    % Create a second version of the figure with the legend and title removed for cleaner presentation.
    Fig.Children(1).Visible = "off"; 
    Fig.Children(2).Title.String = ""; 
    saveas(Fig,[PathName,char(Name),'_Result_NT'],'svg'); 
    % Close the figure to avoid screen clutter during batch processing.
    delete(Fig); 
end 

%% --- EXPORT FINAL RESULTS --- %%

% Convert the results structure into a MATLAB table for easy export.
SkullResT = struct2table(SkullRes); 
% Write the table to an Excel file in the same directory as the images.
writetable(SkullResT,[PathName,'Results.xlsx']); 
end