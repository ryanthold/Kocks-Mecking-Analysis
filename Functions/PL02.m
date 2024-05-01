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

function[EM,PL,off_02]=PL02(ESS,R,S)
%% Check Inputs
if exist('S','var')==0
    S=0;
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
elseif size(S)~=1
    error('minimum stress value (S) for calculation is of the wrong dimensions.')
elseif R>1 || R<0
    error('R^2 input value for modulus calculation (R) must be a value between 0 and 1.')
elseif S<0
    error('minimum stress value (S) for calculation must be ')
end
R=round(R,3);
%% Determine Elastic Modulus of ESS curve
% Delete initial data less than S[MPa]
ESS_c=ESS(ESS(:,2)>S,:);

%Fit first section of data starting at 50% of pts to find the elastic slope
i=0.5;
p0_1=round(length(ESS)*i);
data_begin=ESS_c(1:p0_1,:);
data_begin2=data_begin;
lin_fit_begin=polyfit(data_begin(:,1),data_begin(:,2),1);
yfit=polyval(lin_fit_begin,data_begin(:,1));
SStot=sum((data_begin(:,2)-mean(data_begin(:,2))).^2);
SSres=sum((data_begin(:,2)-yfit).^2);
Rsq=1-SSres/SStot;

%Decreasing point by point until an R^2 value of L is reached and the
%elastic modulus is defined as the slope of this region.
while Rsq<R && p0_1>1
    p0_1=p0_1-1;
    data_begin=ESS_c(1:p0_1,:);
    lin_fit_begin=polyfit(data_begin(:,1),data_begin(:,2),1);
    yfit=polyval(lin_fit_begin,data_begin(:,1));
    SStot=sum((data_begin(:,2)-mean(data_begin(:,2))).^2);
    SSres=sum((data_begin(:,2)-yfit).^2);
    Rsq=1-SSres/SStot;
    if p0_1==2
        error('Data does not reach the designated R^2 value.  Lower L or further smooth data.')
    end
end
EM=lin_fit_begin(1);
ePL= data_begin(end,:);

%% Determine 0.2% Offset Yield
f0_2=@(x) EM*(x-0.002);
data_begin2(:,3)=abs(data_begin2(:,2)-f0_2(data_begin2(:,1)));
Idx_YS_off=find(data_begin2(:,3)==min(data_begin2(:,3)));
YS_off=data_begin2(Idx_YS_off(1),2);
Strain_off=data_begin2(Idx_YS_off(1),1);
off_02=[Strain_off YS_off];

%% Determine Proportional Limit from True Stress Strain Curve
TSS=ESS_to_TSS(ESS);
TSS_c=ESS_to_TSS(ESS_c);
% Determine Elastic Modulus from True Stress Strain Curve
%Fit first section of data starting at 50% of pts to find the elastic slope
j=0.5;
p0_1_t=round(length(TSS)*j);
data_begin_t=TSS_c(1:p0_1_t,:);
lin_fit_begin_t=polyfit(data_begin_t(:,1),data_begin_t(:,2),1);
yfit_t=polyval(lin_fit_begin_t,data_begin_t(:,1));
SStot=sum((data_begin_t(:,2)-mean(data_begin_t(:,2))).^2);
SSres=sum((data_begin_t(:,2)-yfit_t).^2);
Rsq=1-SSres/SStot;

%Decreasing point by point until an R^2 value of L is reached and the
%elastic modulus is defined as the slope of this region.
while Rsq<R && p0_1_t>1
    p0_1_t=p0_1_t-1;
    data_begin_t=TSS_c(1:p0_1_t,:);
    lin_fit_begin_t=polyfit(data_begin_t(:,1),data_begin_t(:,2),1);
    yfit_t=polyval(lin_fit_begin_t,data_begin_t(:,1));
    SStot=sum((data_begin_t(:,2)-mean(data_begin_t(:,2))).^2);
    SSres=sum((data_begin_t(:,2)-yfit_t).^2);
    Rsq=1-SSres/SStot;
    
     %check next point if passes
    if Rsq>=R
        p0_1c=p0_1_t-1;
        data_beginc=TSS_c(1:p0_1c,:);
        lin_fit_beginc=polyfit(data_beginc(:,1),data_beginc(:,2),1);
        yfitc=polyval(lin_fit_beginc,data_beginc(:,1));
        SStot=sum((data_beginc(:,2)-mean(data_beginc(:,2))).^2);
        SSres=sum((data_beginc(:,2)-yfitc).^2);
        Rsq=1-(SSres/SStot);
        p0_1_t=p0_1c;
        Rsq=round(Rsq,3);
    end
    if p0_1_t==2
        error('Data does not reach the designated R^2 value.  Lower L or further smooth data.')
    end

end
PL= data_begin_t(end,:);

%% Plots
% % Uncomment below for plot of engineering stress-strain curve with elastic
% % modulus line, engineering proportional limit, and 0.2% offset 
% figure
% hold on
% grid on
% plot(ESS(:,1),ESS(:,2))
% plot(data_begin(:,1),yfit,'--r')
% fplot(f0_2,[0.002 Strain_off],'--g')
% scatter(ePL(:,1),ePL(:,2),'ok','filled')
% scatter(off_02(:,1),off_02(:,2),'or','filled')
% legend('Engineering Stress-Strain curve','Elastic Modulus Line','Offset Elastic Modulus','Engineering Proportional Limit','0.2% Offset')
% xlabel('Engineering Strain [mm/mm]')
% ylabel('Engineering Stress [MPa]')

% Uncomment below for plot of true stress-strain curve with elastic modulus
% line and proportional limit
figure
hold on
grid on
plot(TSS(:,1),TSS(:,2))
plot(data_begin_t(:,1),yfit_t,'--r')
scatter(PL(:,1),PL(:,2),'ok','filled')
legend('True Stress-Strain curve','Elastic Modulus Line','Proportional Limit')
xlabel('True Strain [mm/mm]')
ylabel('True Stress [MPa]')
end