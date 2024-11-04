%% Calculate Trial Performance
function [penalty_total, penalties] = scoreTrial( ...
    rate_data, sysmon_data, track_data, comm_data, resman_data, ...
    matb_data, events, flow_rates, last_datum, current_datum ...
    )
%% scoreSYSMON
penalty_persec = 0.1;
penalty_total = 0;

% add penalty for every second the RT took
temp_data = sysmon_data.RTs(last_datum(1):current_datum(1));
sysm_penalty = penalty_persec*sum(abs(temp_data(~isnan(temp_data))));
penalty_total = penalty_total + sysm_penalty;
clear temp_data;
%penalty_total = penalty_total + penalty_persec*sum(abs(sysmon_data.RTs(~isnan(sysmon_data.RTs(last_datum(1):current_datum(1))))));

%% scoreTRACK
% add penalty for RMSD during each tracking session
%penalty_total = penalty_total + penalty_persec*sum(abs(track_data.RMSD(~isnan(track_data.RMSD(last_datum(2):current_datum(2))))));

temp_data = track_data.RMSD(last_datum(2):current_datum(2));
track_penalty = penalty_persec*sum(abs(temp_data(~isnan(temp_data))));
penalty_total = penalty_total + track_penalty;
clear temp_data;

%% scoreCOMM
% pull trial comm_data;
temp_comm_data = comm_data(last_datum(3):current_datum(3),:);

% find all correct comm entries 
correct_comms = comm_data.scores(last_datum(3):current_datum(3),:); 
correct_comms = (correct_comms(:,1).*correct_comms(:,2)); 

% find all ships (OWN or OTHER)
comm_ships = comm_data.ship_exp(last_datum(3):current_datum(3));

% find all correct responses to OWN calls & identify event numbers
correct_responses = temp_comm_data(cellfun(@(x) isequal(x,'OWN'), ...
    comm_ships) & correct_comms,:);  
correct_event_nums = correct_responses.event_num;
matb_events = matb_data(ismember(matb_data.event_num, ...
    correct_event_nums),:); 

% penalize every second until the correct responses to OWN calls
time_delays = correct_responses.times-matb_events.times;            % subtract event start times from response times
time_delays = time_delays(:,1)*60+time_delays(:,2);                 % convert min:sec to sec
comm_penalty = penalty_persec*sum(time_delays);                     % sum up delay

% penalize 15 seconds for every incorrect change or lack of change
comm_penalty = comm_penalty + penalty_persec*15*sum(~correct_comms); 
penalty_total = penalty_total + comm_penalty;
clear correct_comms comm_ships time_delays matb_events correct_event_nums...
    correct_responses;

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

% easier to just start from the beginning here than remember last trial's
last_current_time = 2;

% determine current time
temp_time = resman_data.times(current_datum(4),:);
current_time = floor(temp_time(1)*60+temp_time(2));

% iterate through times from beginning to 'current moment'
for j = last_current_time:current_time  

    if count <= length(resman_data.times) && j>=(resman_data.times(count,1)*60+resman_data.times(count,2))  
        % if we hit an event, update flow-rates and fuel levels
        fuel_levels(j,:) = resman_data.fuel_levels(count,1:2);
        % update pump status according to event details
        switch resman_data.actions{count}
            case 'On'
                pump_status(count,resman_data.pumps(count)) = 1;
            case 'Off'
                pump_status(count,resman_data.pumps(count)) = 0;
            case 'Fail'
                pump_status(count,resman_data.pumps(count)) = 2;
            case 'Fix'
                pump_status(count,resman_data.pumps(count)) = 0;
        end

        % determine total flow rates for tanks A & B
        pumps_on = pump_status(count,:) == 1;
        flow_A(j) = pumps_on(1)*flow_rates(1) ...
                    + pumps_on(2)*flow_rates(2) ...
                    + pumps_on(8)*flow_rates(8) ...
                    - pumps_on(7)*flow_rates(7) ...
                    - 800;
        flow_B(j) = pumps_on(3)*flow_rates(3) ...
                    + pumps_on(4)*flow_rates(4) ...
                    + pumps_on(7)*flow_rates(7) ...
                    - pumps_on(8)*flow_rates(8) ...
                    - 800;

        % increment count
        count = count+1; 
    else
        % if no event, flow rate stays unchanged
        flow_A(j) = flow_A(j-1);
        flow_B(j) = flow_B(j-1);
        fuel_levels(j,:) = fuel_levels(j-1,:) ...
                           + [flow_A(j), flow_B(j)].*1/60;
    end
end

tankA_fuel = fuel_levels(last_current_time:current_time,1);
tankB_fuel = fuel_levels(last_current_time:current_time,2);

% violations for any second that either tank is out of bounds
nViolations = 0;
nViolations = nViolations + length(tankA_fuel(tankA_fuel>3000));
nViolations = nViolations + length(tankB_fuel(tankB_fuel>3000));
nViolations = nViolations + length(tankA_fuel(tankA_fuel<2000));
nViolations = nViolations + length(tankB_fuel(tankB_fuel<2000));
clear tankA_fuel tankB_fuel

last_current_time = current_time;

% add resman penalty
resm_penalty = penalty_persec*sum(nViolations);
penalty_total = penalty_total + resm_penalty;

% compile penalties for output
penalties = [sysm_penalty, track_penalty, comm_penalty, resm_penalty];

end