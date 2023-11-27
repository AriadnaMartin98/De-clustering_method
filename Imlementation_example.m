clc; clear all; 
keyboard

%% Implementing the de-cluster method
%--Example for a tide gauge time series of storm surge. 

%---1) Import data time series (time in ISO calendar). 
load("Example_abashiri.mat"); 

%---2) To find the de-cluster parameters
[SED,soft_margin] = find_parameters(t_series); % To specify a different threshold or window (dys): find_parameters(t_series,th,W)

%--3) To de-cluster the time series using those parameters
[Ind_events] = de_cluster(t_series,SED,soft_margin); % To specify a different threshold:  de_cluster(dataframe,SED,soft_margin,th)



















