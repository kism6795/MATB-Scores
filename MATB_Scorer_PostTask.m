%% MATB SA Quiz Scorer
%{ 
    Written by Kieran J Smith, kieran.smith@colorado.edu on August 23, 2022 
    as part of a project in conjunction with Draper Labs to score
    performance on a MATB task
%}
clear; clc;

%% set important folders & parameters
data_path = 'C:\Users\kiera\Documents\Kieran\CU\Research\SA\Subject Data';

subjString = {'001-0922','002-0923','003-0926','004-1025','005-1110',...
              '006-1111','007-1113','008-1215','009-1219','010-1221',...
              '011-0123','012-0125','013-0126','014-0130','015-0202',...
              '016-0208','017-0222','018-0226','019-0227','020-0229',...
              '021-0309','022-0404','023-0412','024-0416','025-0426',...
              '026-0717','027-0724','028-0729','029-0730','030-0802',...
              '031-0806','032-0807','033-0810','034-0823','035-0916'};

nSubs = length(subjString);
nTrials = 12;

adjMATBscores =  nan(nSubs,nTrials);

% Read in subject notes
response_path = "C:\Users\kiera\Documents\Kieran\CU\Research\SA\" + ...
    "Subject Data\Questionnaires\";
notes_file = "All_Subject_Data_Notes.xlsx";

% Import notes without matlab warning
opts = detectImportOptions(fullfile(response_path,notes_file));
opts.SelectedVariableNames = {'ID', 'Trial', 'SAAssessment', 'MATBData', ...
    'fNIRS', 'EEG', 'EyeTracking', 'EDA', 'ECG', 'RSP', 'SAScores', ...
    'MWP', 'SAM', 'Span', 'TLX', 'PVT', 'TaskQuiz', 'NbackTask'};
data_notes = readtable(fullfile(response_path,notes_file), opts);


adj_sub_scores = nan(nTrials,4,nSubs);
sub_penalties = nan(nTrials,4,nSubs);
max_sub_penalties = nan(nTrials,4,nSubs);

flow_rates = [1000, 500, 800, 400, 500, 500, 700, 800, 800];

for s = 5:nSubs  
    fprintf('Calculating MATB Scores for Subject %d / %d\n', s, nSubs);
    subject_folder = sprintf('Subject-%s',subjString{s});

    % Find subject notes
    subj_notes = data_notes(data_notes.ID == s,:);

    % score each trial separately
    penalty_total = zeros(1,nTrials);
    penalty_max = zeros(1,nTrials);

    % Get MATB Performance by Trial
    for i = 1:nTrials

        % Pull data from All_Subject_Data_Notes.xlsx
        trial_notes = subj_notes(subj_notes.Trial == i,:);
        
        % MATB Data is either:
        % 1: Properly collected
        if trial_notes.MATBData == 1

            % Pull info for trials collected correctly
            mb_notes = subj_notes(subj_notes.MATBData == 1,:);

            % Trial number in MATB Data (Nth rate data) in the correct file
            mb_trial_num = find(mb_notes.Trial == i);

            trial_folder = subject_folder;

        % 0: Missing entirely
        elseif trial_notes.MATBData == 0
            fprintf("\tSkipping trial %d.\n", i);
            continue;
            
        % NaN: In a 'missing-trials' sub folder
        elseif isnan(trial_notes.MATBData)
            % Pull info for trials collected on a separate file
            mb_notes = subj_notes(isnan(subj_notes.MATBData),:);

            % Locate missing data folder
            trial_folder = fullfile(subject_folder, 'missing-trials');
            fprintf("\tPulling trial %d MATB Data from missing-trials " + ...
                "folder.\n", i);

        % Some unexpected Other 
        else
            error("Unexpected value in %s.\n Check %s", notes_file, ...
                response_path);
        end

        % Trial number in MATB Data (Nth rate data) in the specified file
        mb_trial_num = find(mb_notes.Trial == i);

        % Set last datum if on first trial in a matb data set
        if mb_trial_num == 1
            last_datum = ones(1,4);
            current_datum = ones(1,4);
        end

        % Import generated XML file
        load(fullfile(data_path,subject_folder,'events.mat'));
        events = events2; clear events2
    
        % Import MATB Data
        [rate,sysmon,track,comm,resman,matb] = getMATBdata( ...
            fullfile(data_path, subject_folder) ...
            );

        % determine current time
        current_time = rate.times(mb_trial_num,1)*60 ...
            + rate.times(mb_trial_num,2);
        current_datum(1) = getLastDatum(sysmon.times, current_time);
        current_datum(2) = getLastDatum(track.times, current_time);
        current_datum(3) = getLastDatum(comm.times, current_time);
        current_datum(4) = getLastDatum(resman.times, current_time);
        
        % calculate trial penalties and trial max possible penalties
        % subtask penalties are ordered: sysm track comm resm 
        [penalty_total(i), sub_penalties(i,:,s)] = scoreTrial( ...
            rate, sysmon, track, comm, resman, ...
            matb, events, flow_rates, last_datum, current_datum ...
            );

        [penalty_max(i), max_sub_penalties(i,:,s)] = worstTrialPerformance( ...
            rate, sysmon, track, comm, resman, ...
            matb, events, flow_rates, last_datum, current_datum ...
            );

        % calculate trial penalties out of max possible penalties
        adj_sub_scores(i,:,s) = sub_penalties(i,:,s)...
                                ./max_sub_penalties(i,:,s);
        adjMATBscores(s,i) =  penalty_total(i)/penalty_max(i);

        % save last datapoint for score calculation
        last_datum = current_datum;
    end
    % 
    % % shift subject 25 scores to be the LAST 8 trials
    % if s == 25
    %     missing_trials = sub_penalties(nT(s)+1:end,:,s);
    %     last_trials = sub_penalties(1:nT(s),:,s);
    %     sub_penalties(:,:,s) = [missing_trials; last_trials];
    % 
    %     missing_trials = adj_sub_scores(nT(s)+1:end,:,s);
    %     last_trials = adj_sub_scores(1:nT(s),:,s);
    %     adj_sub_scores(:,:,s) = [missing_trials; last_trials];
    % 
    %     missing_trials = adjMATBscores(s,nT(s)+1:end);
    %     last_trials = adjMATBscores(s,1:nT(s));
    %     adjMATBscores(s,:) = [missing_trials, last_trials];
    % end
    
end

% Calculate Raw MATB Scores
raw_matb_scores = sub_penalties(:,1,:) ...
                  + sub_penalties(:,2,:) ...
                  + sub_penalties(:,3,:) ...
                  + sub_penalties(:,4,:);

raw_matb_scores = reshape(raw_matb_scores,[],1);

% Calculate Total Scores
adj_sub_scores(:,5,:) = adj_sub_scores(:,1,:) ...
                        + adj_sub_scores(:,2,:) ...
                        + adj_sub_scores(:,3,:) ...
                        + adj_sub_scores(:,4,:);
% save(fullfile(data_path,'adjMATBscores.mat'),'adjMATBscores');

save('adjMATBscores.mat','adjMATBscores');
save('adj_MATB_sub_scores.mat','adj_sub_scores');

%% Z-Score subtask performances
% order of penalties is: [sysmon , track, comms, resman]
order = {'sysmon', 'track', 'comms', 'resman'};
z_adj_sub_scores = nan*adj_sub_scores;
for i = 1:4
    subtask_scores = adj_sub_scores(:,i,:);
    sz = size(subtask_scores);
    subtask_zscores = reshape(subtask_scores,sz(1)*sz(3),sz(2));
    subtask_zscores(~isnan(subtask_zscores)) = zscore(subtask_zscores(~isnan(subtask_zscores)));
    z_adj_sub_scores(:,i,:) = reshape(subtask_zscores,sz);
end

z_adj_sub_scores(:,5,:) = z_adj_sub_scores(:,1,:) ...
                        + z_adj_sub_scores(:,2,:) ...
                        + z_adj_sub_scores(:,3,:) ...
                        + z_adj_sub_scores(:,4,:);

z_adj_matb_scores = reshape(z_adj_sub_scores(:,5,:),[],1);

save('z_adj_MATB_sub_scores.mat','z_adj_sub_scores');
save(fullfile(data_path, sprintf('z_adj_MATB_sub_scores_5-%d.mat',nSubs)),...
     'z_adj_sub_scores');
save(fullfile(data_path, sprintf('raw_matb_scores_5-%d.mat',nSubs)),...
     'raw_matb_scores');

%% Correlate z-scored with raw MATB scores
figure;
plot(reshape(adj_sub_scores(:,5,:),[],1), z_adj_matb_scores, 'ok');
xlabel('Adjusted MATB Penalties');
ylabel('Z-Score Adjusted Penalties Combined')


figure;
plot(raw_matb_scores, z_adj_matb_scores, 'ok');
xlabel('Raw MATB Penalties');
ylabel('Z-Score Adjusted Penalties Combined')

[r,p] = corr(raw_matb_scores(~isnan(z_adj_matb_scores)), ...
     z_adj_matb_scores(~isnan(z_adj_matb_scores)));

%% Plot Histograms of subtask performance
% order of penalties is: [sysmon , track, comms, resman]
order = {'sysmon', 'track', 'comms', 'resman'};

figure;
for i = 1:4
    subplot(2,2,i);
    histogram(adj_sub_scores(:,i,:))
    if kstest(adj_sub_scores(:,i,:))
        title('NOT normally distributed');
    else
        title('normally distributed');
    end
    ylabel('count');
    xlabel(sprintf('%s penalty',order{i}));
end
sgtitle('adj matb scores')

%% Plot Histograms of subtask performance
% order of penalties is: [sysmon , track, comms, resman]
order = {'sysmon', 'track', 'comms', 'resman'};

figure;
for i = 1:4
    subplot(2,2,i);
    histogram(z_adj_sub_scores(:,i,:))
    if kstest(z_adj_sub_scores(:,i,:))
        title('NOT normally distributed');
    else
        title('normally distributed');
    end
    ylabel('count');
    xlabel(sprintf('%s penalty',order{i}));
end
sgtitle('z adj matb scores')

%%
figure;
for i = 1:4
    subplot(2,2,i);
    histogram(log(adj_sub_scores(:,i,:)))
    if kstest(log(adj_sub_scores(:,i,:)))
        title('NOT normally distributed');
    else
        title('normally distributed');
    end
    ylabel('count');
    xlabel(sprintf('%s penalty',order{i}));
end
sgtitle('log matb scores')