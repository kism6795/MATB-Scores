%% Calculate Worst Possible Trial Performance given event data
function [penalty_max, penalties] = worstTrialPerformance( ...
    rate_data, sysmon_data, track_data, comm_data, resman_data, ...
    matb_data, events, flow_rates, last_datum, current_datum ...
    )
penalty_persec = 0.1;
penalty_max = 0;

%% scoreSYSMON
% add 15-second penalty for RT Opportunity
sysm_penalty = penalty_persec * 15 * ...
    height(sysmon_data.RTs(last_datum(1):current_datum(1)));
penalty_max = sysm_penalty;

%% scoreTRACK
% putting in 54.7284 --> identifies outliers in first three participants' data
% i.e. = 3 scaled median absolute deviations from the median
% https://www.mathworks.com/help/matlab/ref/isoutlier.html
track_penalty = penalty_persec*54.7284;
penalty_max = penalty_max + track_penalty; 

%% scoreCOMM
% add 15s penalty for each comm
temp_comm_data = comm_data(last_datum(3):current_datum(3),:);
temp_comm_ships = comm_data.ship_exp(last_datum(3):current_datum(3));
own_comms = temp_comm_data(cellfun(@(x) isequal(x,'OWN'), ...
    temp_comm_ships),:);  
n_own_comms = height(own_comms);
comm_penalty = penalty_persec * 15 * n_own_comms;
penalty_max = penalty_max + comm_penalty;

%% scoreRESMAN
% initialize fuel levels and pump status for tanks A/B
count = 1;
fuel_levels = zeros(length(events.times),2);
fuel_levels(1,:) = [2500, 2500];
pump_status = zeros(length(resman_data.times),8);

% set initial flow rates for tanks A & B
pumps_on = pump_status(1,:) == 1;
flow_A(1) = -800;
flow_B(1) = -800;

% start time
last_current_time = 2;

% determine current time
temp_time = resman_data.times(current_datum(4),:);
current_time = floor(temp_time(1)*60+temp_time(2));

% can't go to event file length bc scores each trial...    
for j = last_current_time:current_time  

    if count <= length(resman_data.times) && j>=(resman_data.times(count,1)*60+resman_data.times(count,2))  
        % if we hit an event, update flow-rates and fuel levels
        fuel_levels(j,:) = resman_data.fuel_levels(count,1:2);
        % update pump status according to event details
        switch resman_data.actions{count}
            case 'On'
                % Don't turn any pumps on  (worst case)
                % pump_status(count,resman_data.pumps(count)) = 1;
            case 'Off'
                % don't turn off pumps (no input mode)
                % I believe this isn't counting auto-shut-offs
                % pump_status(count,resman_data.pumps(count)) = 0;
            case 'Fail'
                pump_status(count,resman_data.pumps(count)) = 2;
            case 'Fix'
                pump_status(count,resman_data.pumps(count)) = 0;
        end
        % determine total flow rates for tanks A & B
        pumps_on = pump_status(count,:) == 1;
        flow_A(j) = pumps_on(1)*flow_rates(1)...
                    + pumps_on(2)*flow_rates(2)... 
                    + pumps_on(8)*flow_rates(8)...
                    - pumps_on(7)*flow_rates(7)...
                    - 800;
        flow_B(j) = pumps_on(3)*flow_rates(3)...
                    + pumps_on(4)*flow_rates(4)...
                    + pumps_on(7)*flow_rates(7)...
                    - pumps_on(8)*flow_rates(8)...
                    - 800;

        % increment count
        count = count+1;
    else
        % if no event, flow rate stays unchanged
        flow_A(j) = flow_A(j-1);
        flow_B(j) = flow_B(j-1);
        fuel_levels(j,:) = fuel_levels(j-1,:) + [flow_A(j), flow_B(j)].*1/60;
    end
end

temp_data1 = fuel_levels(last_current_time:current_time,1);
temp_data2 = fuel_levels(last_current_time:current_time,2);

nViolations = 0;
nViolations = nViolations + length(temp_data1(temp_data1>3000));
nViolations = nViolations + length(temp_data2(temp_data2>3000));
nViolations = nViolations + length(temp_data1(temp_data1<2000));
nViolations = nViolations + length(temp_data2(temp_data2<2000));
last_current_time = current_time;
clear temp_data1 temp_data2

resm_penalty = penalty_persec*sum(nViolations);

% ensure that a divide by zero doesn't occur
% through subj 21, no one has gotten any penalty on a round where 'no
% input' would have resulted in no penalty.
if resm_penalty == 0
    resm_penalty = 0.000001;
end
penalty_max = penalty_max + resm_penalty;


% compile penalties for output
penalties = [sysm_penalty, track_penalty, comm_penalty, resm_penalty];

end