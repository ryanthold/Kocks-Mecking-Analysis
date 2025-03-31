%% This function is used to determine the Proportional Limit, 0.2% Offset, and Elastic Modulus.
%% Inputs
% ESS: 2-by-x double array of engineering stress[MPa]-strain[mm/mm] data with strain in the first column and stress in the second column
% R: double between 0 and 1 representing allowable R^2 value for Elastic
% Region
% S: double representing the stress value [MPa] below which the data is deleted.
% This is done to eliminate fluctuations at the start of loading.  If not
% entered, will be set to 0 [MPa]. 

%% Outputs
% EM: double representing the Elastic Modulus [MPa]
% PL: 2-by-1 double array representing the proportional limit with strain [mm/mm] in the first column and stress[MPa] in the second column
% off_02: 2-by-1 double array representing the 0.2% offset yield point with strain[mm/mm] in the first column and stress[MPa] in the second column

function[EM,PL,off_02]=PL02_4(ESS,R,R_PL,S)
%% Check Inputs
if exist('S','var')==0
    S=0;
elseif size(S)~=1
    error('minimum stress value (S) for calculation is of the wrong dimensions.')
elseif S<0
    error('minimum stress value (S) for calculation must be greater or equal than 0')
end
if class(ESS)~="double"
    error('engineering stress-strain (ESS) input is of the wrong data type.')
elseif class(R)~="double"
    error('R^2 input value for modulus calculation (R) is of the wrong data type.')
elseif class(S)~="double"
    error('minimum stress value (S) for calculation is of the wrong data type.')
elseif size(ESS,2)~=2
    error('engineering stress-strain (ESS) input is of the wrong dimensions.')
elseif size(R)~=1
    error('R^2 input value for modulus calculation (R) is of the wrong dimensions.')
elseif R>1 || R<0
    error('R^2 input value for modulus calculation (R) must be a value between 0 and 1.')
end

R=round(R,3);
%% Determine Elastic Modulus of ESS curve
% Delete initial data less than S[MPa]
ESS_c=ESS(ESS(:,2)>S,:);

% Fit point by point before going below R2 value starting with 10 points
Rsq=1;
p0=10;
while Rsq>=R
    p0=p0+1;
    data_check=ESS_c(1:p0,:);
    [~,D]=polyfit(data_check(:,1),data_check(:,2),1);
    Rsq=1-(D.normr/norm(data_check(:,2)-mean(data_check(:,2))))^2;
    %test next data point
    if Rsq<R
        p0=p0+1;
        data_check=ESS_c(1:p0,:);
        [~,D]=polyfit(data_check(:,1),data_check(:,2),1);
        Rsq=1-(D.normr/norm(data_check(:,2)-mean(data_check(:,2))))^2;
    end
end
ES=polyfit(data_check(1:end-2,1),data_check(1:end-2,2),1);
EM=ES(1);
data_check_e=data_check;
yfit_e=polyval(ES,data_check_e(:,1));

%% Determine 0.2% Offset Yield
f0_2=@(x) EM*(x-0.002);
ESS_c(:,3)=abs(ESS_c(:,2)-f0_2(ESS_c(:,1)));
Idx_YS_off=find(ESS_c(:,3)==min(ESS_c(:,3)));
YS_off=ESS_c(Idx_YS_off(1),2);
Strain_off=ESS_c(Idx_YS_off(1),1);
off_02=[Strain_off YS_off];

%% Determine Proportional Limit from True Stress Strain Curve
TSS_c=ESS_to_TSS(ESS_c);
% Determine Elastic Modulus from TSS Curve
Rsq=1;
p0=10;
while Rsq>=R
    p0=p0+1;
    data_check=TSS_c(1:p0,:);
    [~,D]=polyfit(data_check(:,1),data_check(:,2),1);
    Rsq=1-(D.normr/norm(data_check(:,2)-mean(data_check(:,2))))^2;
    %test next data point
    if Rsq<R
        p0=p0+1;
        data_check=TSS_c(1:p0,:);
        [~,D]=polyfit(data_check(:,1),data_check(:,2),1);
        Rsq=1-(D.normr/norm(data_check(:,2)-mean(data_check(:,2))))^2;
    end
end
ES_t=polyfit(data_check(1:end-1,1),data_check(1:end-1,2),1);

% Determine proportional limit as first point with an R^2 below R_PL when
% compared to Curve of Elastic Modulus
Rsq=1;
while Rsq>R_PL
    p0=p0+1;
    data_check=TSS_c(1:p0,:);
    yfit=polyval(ES_t,data_check(:,1));
    SStot_t=sum((data_check(:,2)-mean(data_check(:,2))).^2);
    SSres_t=sum((data_check(:,2)-yfit).^2);
    Rsq=1-SSres_t/SStot_t;
    if Rsq<R_PL
        p0=p0+1;
        data_check=TSS_c(1:p0,:);
        yfit=polyval(ES_t,data_check(:,1));
        SStot_t=sum((data_check(:,2)-mean(data_check(:,2))).^2);
        SSres_t=sum((data_check(:,2)-yfit).^2);
        Rsq=1-SSres_t/SStot_t;
    end
end
PL= data_check(end-1,:);
data_check_t=data_check;
yfit_t=polyval(ES_t,data_check_t(:,1));


% %% Plots
% % Uncomment below for plot of engineering stress-strain curve with elastic
% % modulus line, engineering proportional limit, and 0.2% offset 
% figure
% hold on
% grid on
% plot(ESS(:,1),ESS(:,2))
% plot(data_check_e(:,1),data_check_e(:,2),'k')
% plot(data_check_e(:,1),yfit_e,'--r')
% fplot(f0_2,[0.002 Strain_off],'--g')
% scatter(off_02(:,1),off_02(:,2),'or','filled')
% legend('Engineering Stress-Strain curve','Elastic Modulus Line','Offset Elastic Modulus','0.2% Offset')
% xlabel('Engineering Strain [mm/mm]')
% ylabel('Engineering Stress [MPa]')
% 
% % Uncomment below for plot of true stress-strain curve with elastic modulus
% % line and proportional limit
% figure
% hold on
% grid on
% plot(TSS(:,1),TSS(:,2))
% plot(data_check_t(:,1),yfit_t,'--r')
% scatter(PL(:,1),PL(:,2),'ok','filled')
% legend('True Stress-Strain curve','Elastic Modulus Line','Proportional Limit')
% xlabel('True Strain [mm/mm]')
% ylabel('True Stress [MPa]')
end
