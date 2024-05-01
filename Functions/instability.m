%% This function is used to determine the Instability Limit through Considere's Criterion
%% Inputs
% ESS: 2-by-x double array of engineering stress[MPa]-strain[mm/mm] data with strain in the first column and stress in the second column
% n: an odd double greater than 1 that specifies the number of points used
% for calculation of the instantaneous slopes.
%% Outputs
% I_e: 2-by-1 double array representing the instability limit with engineering strain [mm/mm] in the first column and engineering stress[MPa] in the second column
% I_t: 2-by-1 double array representing the instability limit with true strain [mm/mm] in the first column and true stress[MPa] in the second column

function[I_e,I_t]=instability(ESS,n)
%% Check Inputs
if class(ESS)~="double"
    error('engineering stress-strain (ESS) input is of the wrong data type.')
elseif class(n)~="double"
    error('n input value is of the wrong data type.')
elseif size(ESS,2)~=2
    error('engineering stress-strain (ESS) input is of the wrong dimensions.')
elseif size(n)~=1
    error('n input value is of the wrong dimensions.')
elseif n<2 || round(n/2)==(n/2)
    error('n input value must be odd and greater than 1')
end

%% Create engineering strain vs true stress matrix for construction and find maximum. Eliminate data after maximum value.
TSS=ESS_to_TSS(ESS);
data=[ESS(:,1) TSS(:,2)];
sts_max_lim=max(data(:,2));
[Idx_max,~]=find(data==sts_max_lim);
data_c=data(1:Idx_max,:);
slope_list=ones(length(ESS),1).*100;
sub_list=ones(length(ESS),1).*100;
%% Use n points to calculate slope of each point and find where it is equal to the subtangent to -1.
for i=ceil(n/2):(length(data_c)-floor(n/2))
    strain_n=data_c(i,1);
    stress_n=data_c(i,2);
    xn=data_c(i-floor(n/2):i+floor(n/2),1);
    yn=data_c(i-floor(n/2):i+floor(n/2),2);
    tan=polyfit(xn,yn,1);
    slope_n=tan(1);
    %The slope is then subtracted by the right side of Considere's equation for each point
    sub=abs(slope_n-(stress_n/(1+strain_n)));
    slope_list(i)=slope_n;
    sub_list(i)=sub;
end
%The minimum value of the slope subtracted by the right side of Considere's
%equation is determined to be the instability point.
    
%Goes through this loop to check if the points surrounding the minimum are close to
%being subtangent as well.  Does this for 10 closest points and finds
%overall minimum.
Idx_list=zeros(10,1);
sub_ard=zeros(10,1);
for i=1:10
    Min=min(sub_list);
    Idx=find(sub_list==Min);
    if length(Idx)>1
        for j=1:length(Idx)
            if Idx(j)<(length(ESS)/2)
                Idx(j)=0;
            end
        end
        Idx=Idx(Idx~=0);
        Idx=Idx(1);
    end
    if Idx>3 && Idx<(length(sub_list)-3)
        V5=sub_list(Idx-3)+sub_list(Idx-2)+sub_list(Idx-1)+sub_list(Idx)+sub_list(Idx+1)+sub_list(Idx+2)+sub_list(Idx+3);
        Idx_list(i)=Idx;
        sub_ard(i)=V5;
    end
    sub_list(Idx)=100;
end
V5_min_idx=sub_ard==min(sub_ard);
Idx=Idx_list(V5_min_idx);
min_slp=slope_list(Idx);
    
%Instability point and UTS are calculated and plotted
strain_e=data_c(Idx,1);
strain_t=log(strain_e+1);
stress_t=data_c(Idx,2);
stress_e=stress_t./(1+strain_e);

I_e=[strain_e,stress_e];
I_t=[strain_t,stress_t];
    
% %Uncomment below for plot of tangent construction and calculated
% %instability point
% f=@(x) (min_slp)*(x-strain_e)+stress_t;
% figure
% hold on
% grid on
% plot(data(:,1),data(:,2))
% fplot(f,[-1 1],'r')
% scatter(strain_e,stress_t,'ok','filled')    
% xlabel('Engineering Strain [mm/mm]')
% ylabel('True Stress [MPa]')
% legend('True Stress Eng Strain','Considere Construction','Instability Point')
end