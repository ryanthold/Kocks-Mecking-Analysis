%% This function serves to remove the points at the end of the stress
% strain curve, which are a result of post-failure.  These points tend
% to have a behavior characterized by a sharp decrease in stress and
% a stress jump of 10 times the median is used to determine the cutoff.

%% Inputs 
% SS: 2-by-x double array of stress-strain data with strain in the first 
% column and stress in the second column

%% Outputs
% SS_cut: 2-by-x double array of cut stress-strain data with strain in the
% first column and stress in the second column

function[SS_c]=cutend(SS)
%% Check Inputs
if class(SS)~="double"
    error('stress-strain (SS) input is of the wrong data type.')
elseif size(SS,2)~=2
    error('stress-strain (SS) input is of the wrong dimensions.')
end
%% Delete end data after fracture
% Determine stress difference between datapoints in last half of data and median value
d=abs(diff(SS(round(length(SS)/2)+1:end,2)));
m=median(d);

% Use a cutoff of 20 times the median difference
Idx_d=find((d<20*m)==0);
Idx_SS=Idx_d+round(length(SS)/2);
if ~isempty(Idx_d)
    SS_c=SS(1:Idx_SS(1),:);
else
    SS_c=SS;
end
% Uncomment below to plot new dataset on top of raw dataset
% figure
% hold on
% plot(SS(:,1),SS(:,2))
% plot(SS_c(:,1),SS_c(:,2))
% xlabel('Strain [mm/mm]')
% ylabel('Stress [MPa]')
% legend('SS','SS_cut')
end