
function [Ind_events] = de_cluster(dataframe,SED,soft_margin,~)

    % This functions de-cluster the time series (providing the independent events)

    % INPUTS:
    % dataframe: Must be an array [time,data] (time in ISO calendar)
    % SED,soft_margin: parameters found using "find_parameter" function.

    %---Optional inputs: 
    % th: Threshold (99th percentile by default).It can be change depending on the preference of the type of event to be studied (more or less extreme)
    % If you want to compute return levels and used them as thresholds. The "find_parameters" function can be used first to find the window to use in the EVA. 
    % Then those return levels can be include in this function as thresholds. 

    % OUTPUTS:
    % Ind_events: Array [time, event] of the de-clustered time series. 

    % For more information about the methodology: PUT PAPER

    %----------------------------------------------------------------------
    %---Variables

    time = dataframe(:,1);
    data = dataframe(:,2);

    th =  quantile(data, .99);
    W = (SED/24)/2; %--- centered (dys)

    %---1) De-cluster using SED 

    EOT = data;
    EOT(EOT< th) = nan;
    EOT = [time EOT];
    EOT(isnan(EOT(:,2)),:) = [];
    
    xnan = sum(isnan(EOT(:,2)));
    POT = NaN(length(EOT),2);

    while xnan< size(EOT,1)

        [~,fmax] = max(EOT(:,2));

        POT(fmax,:) = EOT(fmax,:);

        dec_wind=  find(EOT(:,1) >= EOT(fmax,1)-W & EOT(:,1) <= EOT(fmax,1)+W);

        EOT(dec_wind,2) = nan;

        xnan = sum(isnan(EOT(:,2)));

    end

    POT(isnan(POT(:,2)),:) = [];

    clear EOT 

    %---2) Applying soft margin
    
    %---Time between events
    s_event= POT(:,1) - W;
    e_event= POT(:,1) + W;

    dif_events = (s_event(2:end) - e_event(1:end-1)) * 24; %---hours CHANGE THAT IF THE TIME EXIT FROM FUNCTION IS DIFFERENT THAN HRS

    t=find(round(dif_events) <= soft_margin); %---those within the soft_margin

    %---Merge events in t
    xnan=length(t);
    
    if xnan~=0

        EOT = NaN(length(e_event),2);
        EOT(:,2) = POT(:,2);
        EOT(1:end-1,1)= dif_events;

        N_POT = POT;
               
        for k=1:xnan
             P=t(k);
             
             fmin=min(EOT(P:P+1,2));
               
             h=find(EOT(:,2)==fmin);

             fmax(k,1)=find(EOT(:,2)==max(EOT(P:P+1,2)));

             fmax(k,2)=h(end);
        
             EOT(h,2)=0;
        
        end

        N_POT(:,2)=EOT(:,2); N_POT(N_POT(:,2)==0,:)=[];
               
        Ind_events(:,1) = N_POT(:,1);
        Ind_events(:,2) = N_POT(:,2);
    
    else
        
        Ind_events(:,1) = POT(:,1);
        Ind_events(:,2) = POT(:,2);

    end

end









