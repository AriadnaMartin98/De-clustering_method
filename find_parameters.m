
function [SED,soft_margin] = find_parameters(dataframe, ~, ~)
    
    % This functions provides the parameters needed to de-cluster the time series in hours, one can change to time steps if needed. 
    
    % INPUTS:
    % dataframe: Must be an array [time,data] (time in ISO calendar)

    %---Optional inputs: 

    % th: Threshold (99th percentile by default).It can be change depending on the preference of the type of event to be studied (more or less extreme)
    % W: Sampling window (6-days by default). It has to be widht enough to capture an entire event. Must be in days. 

    % OUTPUTS:
    % SED: 1st parameter to decluster (standard event duration). (in time steps)
    % soft_margin: 2st parameter to decluster. (in time steps)


    % For more information about the methodology: PUT PAPER

    %----------------------------------------------------------------------
    % EXAMPLE
    % [SED,soft_margin] = find_parameters(dataframe)
    % If SED = 5 it means 2.5 time steps before and after the peak of the event. 
    % If soft_margin = 2 means, to keep at least 2 time steps betweenn the end of an event and the begginins of the next one. 
    
    %----------------------------------------------------------------------
    %---Variables

    time = dataframe(:,1);
    data = dataframe(:,2);

    th =  quantile(data, .99);
    W = 6; %--day

    %---Events over threshold

    EOT = data;
    EOT(EOT< th) = nan;
    EOT = [time EOT];
    EOT(isnan(EOT(:,2)),:) = [];

    %---Keep max peak in the window

    xnan = sum(isnan(EOT(:,2)));
    POT = NaN(length(EOT),2);

    while xnan< size(EOT,1)

        [~,fmax] = max(EOT(:,2));

        POT(fmax,:) = EOT(fmax,:);

        dec_wind=  find(EOT(:,1) >= EOT(fmax,1)-W/2 & EOT(:,1) <= EOT(fmax,1)+W/2);

        EOT(dec_wind,2) = nan;

        xnan = sum(isnan(EOT(:,2)));

    end

    POT(isnan(POT(:,2)),:) = [];

    %---Time series for each peak

    sl = [];

    for wd = 1:length(POT)

        cent = POT(wd,1);
    
        window = find(time >= cent-W/2 & time <= cent+ W/2);
    
        sl(1:length(data(window)),wd) = data(window);

    end
    
    %---Correlation function
    
    cor = corrcoef(sl.','rows','complete');
    M = nanmean(cor,2);
    S = nanstd(cor,[],2);

    C = length(M)/2;
    
    mx = find(max(M) == M);
    mx_s = find(max(M+S) == (M+S));

    SED = abs(mx - C)*2;
    soft_margin = abs(SED/2 - abs(mx_s - C));

end