%% This function serves to calculate the work hardening rate, or instantaneous slope of the true stress strain curve
%% Inputs 
% TSS: 2-by-x double array of true stress-strain data with strain [mm/mm] in the first 
% column and stress [MPa] in the second column

%% Outputs
% wh: 1-by-(x-2) double array of work hardening rate [MPa]
% TSS_wh: 2-by-(x-2) double array of true strain [mm/mm] and true stress [MPa]
% with the values matching the index of wh. True strain is in column 1 and
% true stress is in column 2.

function[wh,TSS]=wrkhard2(TSS)
%% Check Inputs
if class(TSS)~="double"
    error('stress-strain (SS) input is of the wrong data type.')
elseif size(TSS,2)~=2
    error('stress-strain (SS) input is of the wrong dimensions.')
end
%% Calculate instantaneous slope of each datapoint from points adjacent datapoints
wh=diff(TSS(:,2))./diff(TSS(:,1));
TSS=TSS(1:(end-1),:);
end