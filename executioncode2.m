%% This function is used to execute all functions to acquire mechanical properties and cbi values starting from engineering stress[MPa]-strain[mm/mm] data 
%% Inputs
% ESS: a 1-by-x cell array with each component containing x-by-2 double arrays with engineering strain[mm/mm] in the first column and engineering stress [MPa] in the second column.
% p: 1-by-x array with each value between 0 and 1 representing the smoothing parameter
% div: a 1-by-x array whose size is 1 smaller than p with each value representing strain values to distinguish the regions each spline fit occurs on
% R: double between 0 and 1 representing allowable R^2 value for Elastic Region
% S: double representing the stress value [MPa] below which the data is deleted. This is done to eliminate fluctuations at the start of loading.  If not entered, will be set to 0 [MPa]. 
% n: an odd double greater than 1 that specifies the number of points used for calculation of the instantaneous slopes.
% dev3: a double between 0 and 1 representing the allowable deviation in slope between sections to be included in the overall slope measurement of region 3
% Ss3: a double between 0 and 1 representing the percentage of total plastic stress to be used for section size for region 3 calculation
% p3: a double between 0 and 1 representing the spline smoothing parameter for the curve of the derivative of the KM curve
% dev4: a double between 0 and 1 representing the allowable deviation in slope between regions to be included in the overall slope measurement of region 4
% Ss4: a double between 0 and 1 representing the percentage of total stress to be used for section size for region 4 calculation
% rp4: a double between 0 and 1 representing the percentage between the PL and I_t to be used as the reference point for region 4

%% Outputs
% Mech_Props: Table of Mechanical properties with everything in MPa and
%  unitless strain calculated from the TSS curve except for the 0.2% offset,
%  Elongation and UTS, which is from the ESS curve.
% ESS_cs: a 1-by-x cell array where each component is the smoothed and cut
%  engineering stress[MPa]-strain[mm/mm] curves produced. Strain is in the
%  first column and stress is in the second column.
% TSS_cs: a 1-by-x cell array where each component is the smoothed and cut
%  true stress[MPa]-strain[mm/mm] curves produced. Strain is in the first column and
%  stress is in the second column.
% KM: a 1-by-x cell array where each component is the KM curve produced
%  from the smoothed and cut stress-strain curve.  True stress [MPa] is in
%  the first column and work hardening rate [MPa] is in the second column.



function[Mech_Props,ESS_cs,TSS_cs,KM]=executioncode2(ESS,p,div,R,S,n,dev3,Ss3,p3,dev4,Ss4,rp4)
%% Check inputs
if class(ESS)~="cell"
    error('engineering stress-strain (ESS) input must be cell data type.')
elseif size(ESS,1)~=1
    error('engineering stress-strain (ESS) input is of the wrong dimensions.')
end
if class(p)~="double"
    error('p input is of the wrong data type')
end
if class(div)~="double"
    error('div input is of the wrong data type.')
end
if length(p)-length(div)~=1
    error('length of p is not 1 greater than div.')
end

%% Run through loop to calculate mechanical properties
E=ones(length(ESS),1);
PL=ones(length(ESS),1);
off_02=ones(length(ESS),1);
I=ones(length(ESS),1);
ep=ones(length(ESS),1);
Elong=ones(length(ESS),1);
UTS=ones(length(ESS),1);
cb3=ones(length(ESS),1);
sigma_03=ones(length(ESS),1);
cb4=ones(length(ESS),1);
sigma_04=ones(length(ESS),1);
ESS_cs=cell(length(ESS),1);
TSS_cs=cell(length(ESS),1);
KM=cell(length(ESS),1);

for i=1:length(ESS)
    % Create smoothed ESS Curve
    ESSi=ESS{i};
    if sum(ESSi(:,1)>1)>0
        error('make sure strain is unitless [mm/mm].')
    end
    ESS_ci=cutend(ESSi);
    if isempty(div)
        ESS_cip=ESS_ci(:,1).*100;
        [~,gof,ESS_csi]=spline_fit([ESS_cip(:,1),ESS_ci(:,2)],p);
        ESS_csi(:,1)=ESS_csi(:,1)./100;
        rsquare=gof.rsquare;
        % if rsquare<0.99
        %     error("spline fit below rsquare of 0.99 for ESS curve #"+num2str(i))
        % end
        ESS_cs{i}=ESS_csi;        
    else
        id=0;
        ESS_csi=[];
        for j=1:length(div)
            idb=id+1;
            id=find(abs(ESS_ci(:,1)-div(j))==min(abs(ESS_ci(:,1)-div(j))),1);
            ESS_cij{j}=ESS_ci(idb:id,:);
            if j==length(div)
                ESS_cij{j+1}=ESS_ci(id+1:end,:);
            end
        end
        for k=1:length(ESS_cij)
            ESS_cijp=ESS_cij{k}(:,1).*100;
            [~,gof,ESS_k]=spline_fit([ESS_cijp(:,1),ESS_cij{k}(:,2)],p(k));
            rsquare=gof.rsquare;
            % if rsquare<0.99
            %     error("spline fit below rsquare of 0.99 for ESS curve #"+num2str(i))
            % end
            ESS_k(:,1)=ESS_k(:,1)./100;
            ESS_csi=[ESS_csi;ESS_k];
        end
        ESS_cs{i}=ESS_csi;        
    end

    Elong(i)=max(ESS_csi(:,1));
    UTS(i)=max(ESS_csi(:,2));


    % Calculate proportional limit, 0.2% offset, and instability values
    R_PL=0.99;
    [E(i),PLi,off_02i]=PL02_4(ESS_csi,R,R_PL,S);
    [~,I_ti]=instability(ESS_csi,n);

    % Calculate work hardening rate and cbi values
    TSS_cs{i}=ESS_to_TSS(ESS_csi);
    [whi,TSS_i]=wrkhard(TSS_cs{i});
    KM{i}=[TSS_i(:,2),whi];
    [cb3(i),KM{i},reg_num,~,sigma_03(i)]=cb3_calc(whi,TSS_i,PLi,off_02i,I_ti,dev3,Ss3,p3);
    disp(reg_num)
    [cb4(i),~,~,~,sigma_04(i)]=cb4_calc(whi,TSS_i,PLi,I_ti,dev4,Ss4,rp4);

    % Assign values
    PL(i)=PLi(2);
    PLi
    off_02(i)=off_02i(2);
    I(i)=I_ti(2);
    I_ti
    ep(i)=I_ti(1)-PLi(1);

end
%% Create Table 
Mech_Props=table(E,PL,off_02,I,ep,UTS,Elong,cb3,sigma_03,cb4,sigma_04);
end