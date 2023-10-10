%% Initialization
clear;
close;
clc;

% Defining Paths  **BEFRORE THIS CODE PLEASE RUN
% THIRTYSECONDSEGMENTATION.MAT

segmentPath = 'E:\Ph.D\Microstates_Classification\MS_Analysis\raw data';
savePath     = 'E:\Ph.D\Microstates_Classification\MS_Analysis\Results';

% Defining clustering parameters
ClustPars.MinClasses = 7;         % minimum number of microstates
ClustPars.MaxClasses = 8;         % maximum number of microstates
ClustPars.GFPPeaks = true;        % use GFP peaks as inputs of clustering algorithm
ClustPars.IgnorePolarity = true;
ClustPars.MaxMaps = inf;          % number of time sample data feeded to clusterig algorithm 
ClustPars.Restarts = 20;          % number of repeatition for modified k-means 
ClustPars.UseAAHC = false;        % if true, TAAGC will be applied, else modified k-menas.

% Defining fitting parameters
FitPars.b = 30;                   % the window size, if labels are smoothed
FitPars.lambda = 1;               % the non-smoothness penalty factor
FitPars.PeakFit = true;           % true if fitting is only done on GFP peaks and the assignment is interpolated in between, false otherwise
FitPars.BControl = true;

eeglab

%% Reading the data
segmentIndex = [];
for i = 1:2009
    % Basic file read
    EEG = pop_loadset([segmentPath '\Seg' num2str(i) '.set']);  
    % And make this a new set
    [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, 0, 'gui','off');
    % Store the thing
    ALLEEG = eeg_store(ALLEEG, EEG, CURRENTSET); 
    % And keep track of the segments
    segmentIndex = [segmentIndex CURRENTSET];                                           
end

eeglab redraw

%% Microstates for each segment
% Loop across all segments to identify the segment-wise clusters
for i = 1:numel(segmentIndex) 
    % the EEG we want to work with
    tmpEEG = eeg_retrieve(ALLEEG,segmentIndex(i)); 
    % Z-transform normalization
    tempEEG = computeZtransform(tmpEEG);
    % Some info for the impatient user
    fprintf(1,'Clustering dataset %s (%i/%i)\n',tempEEG.setname,i,numel(segmentIndex));
    % This is the actual clustering within subjects
    tempEEG = pop_FindMSTemplates(tempEEG, ClustPars); 
    % Done, we just need to store this
    ALLEEG = eeg_store(ALLEEG, tempEEG, segmentIndex(i)); 
end

eeglab redraw

%% Combining the microstate maps across segments
% The mean accross all segments
EEG = pop_CombMSTemplates(ALLEEG, segmentIndex, 0, 0, 'GrandMean');
[ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, numel(ALLEEG)+1, 'gui', 'off'); 
[ALLEEG,EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET); 
Grand_Mean_Index = CURRENTSET;

%% Sorting things out over means and subjects
% Then, we sort the individuals based on the group mean
ALLEEG = pop_SortMSTemplates(ALLEEG, segmentIndex, 0, Grand_Mean_Index);

eeglab redraw

%% Savings
for i = 1:numel(ALLEEG)
    EEG = eeg_retrieve(ALLEEG,i);
    pop_saveset(EEG, 'filename', EEG.setname, 'filepath', savePath, 'savemode', 'onefile');
    i
end

%% Computing Parameters
for i = ClustPars.MinClasses:ClustPars.MaxClasses 
    FitPars.nClasses = i;
    FileNameGGM = ['Grand Mean Template Maps' num2str(i) ' Lambda' num2str(FitPars.lambda) ' b' num2str(FitPars.b) '.csv'];
    pop_QuantMSTemplates(ALLEEG,  segmentIndex, 1, FitPars, numel(segmentIndex)+1, fullfile(savePath,FileNameGGM));
    i
end

%% Z-transform
function EEG = computeZtransform(tempEEG)
    EEG = tempEEG;
    EEG.data = (tempEEG.data - mean(tempEEG.data, 2)) ./ std(tempEEG.data, 0, 2);
end

%% Labeling 
% for i = 1:2009
%     ALLEEG(i).msinfo.FitParams.b = FitPars.b;                
%     ALLEEG(i).msinfo.FitParams.lambda = FitPars.lambda;      
%     ALLEEG(i).msinfo.FitParams.PeakFit = FitPars.PeakFit;    
%     ALLEEG(i).msinfo.FitParams.BControl = FitPars.BControl;
%     for k = ClustPars.MinClasses:ClustPars.MaxClasses
%         [ALLEEG(i).msinfo.L_all{k,1},ALLEEG(i).msinfo.gfp(1,:), ALLEEG(i).msinfo.fit(k,:)] = AssignMStates(ALLEEG(i), ALLEEG(2010).msinfo.MSMaps(k).Maps, ALLEEG(i).msinfo.FitParams, true);
%     end
%     i
% end