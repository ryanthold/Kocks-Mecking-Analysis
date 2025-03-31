%% This function is used to determine the values of cb4 through slope construction on the KM curve
%% Inputs
% wh: 1-by-(x-2) double array of work hardening rate [MPa]
% TSS: 2-by-x double array of true stress[MPa]-strain[mm/mm] data with strain in the first column and stress in the second column
% PL: 1-by-2 double array representing the proportional limit with strain [mm/mm] in the first column and stress[MPa] in the second column
% I_t: 1-by-2 double array representing the instability with true strain [mm/mm] in the first column and true stress[MPa] in the second column
% dev4: a double between 0 and 1 representing the allowable deviation in slope between regions to be included in the overall slope measurement of region 4
% Ss: a double between 0 and 1 representing the percentage of total stress to be used for section size
% rp4: a double between 0 and 1 representing the percentage between the PL and I_t to be used as the reference point for region 4
%% Outputs
% cb4: a double with the value of the slope of region 4
% KM: a 2-by-x double with stress[MPa] in column 1 and work hardening[GPa] in column 2
%regions_t4: outputs the section numbers of region 4 used for slope calculation
%StdSt_reg4: outputs points of region 4 used for slope calculation
%sigma_04: the y-intercept of the cb4 line construct

function[cb4,KM,regions_t4,StdSt_reg4,sigma_04]=cb4_calc(wh,TSS,PL,I_t,dev4,Ss,rp4)
%% Check Inputs
if class(wh)~="double"
    error('work hardening rate (wh) input is of the wrong data type.')
elseif class(TSS)~="double"
    error('true stress-strain (TSS) input is of the wrong data type.')
elseif class(PL)~="double"
    error('proportional limit (PL) input is of the wrong data type.')
elseif class(I_t)~="double"
    error('Instability (I_t) input is of the wrong data type.')
elseif class(dev4)~="double"
    error('deviation value for region 4 (dev4) input is of the wrong data type.')
elseif class(Ss)~="double"
    error('section size value (Ss) input is of the wrong data type.')
elseif class(rp4)~="double"
    error('rp4 input is of the wrong data type.')
elseif size(wh,2)~=1
    error('work hardening rate (wh) input is of the wrong dimensions.')
elseif size(TSS,2)~=2
    error('true stress-strain (TSS) input is of the wrong dimensions.')
elseif size(PL,2)~=2
    error('proportional limit (PL) input is of the wrong dimensions.')
elseif size(I_t,2)~=2
    error('instability (I_t) input is of the wrong dimensions.')
elseif size(dev4)~=1
    error('dev4 input value is of the wrong dimensions.')
elseif size(Ss)~=1
    error('section size (Ss) input value is of the wrong dimensions.')    
elseif size(rp4)~=1
    error('rp4 input value is of the wrong dimensions.')    
elseif dev4<0 || dev4>1
    error('dev4 input value must be between 0 and 1')
elseif Ss<0 || Ss>1
    error('section size (Ss) input value must be between 0 and 1')
elseif rp4<0 || rp4>1
    error('rp4 input value must be between 0 and 1')
end

%% Determine reference point for region 4
KMs=[TSS(:,2) wh TSS(:,1)];
KMs(any(isnan(KMs),2),:)=[];
KM=KMs;
KMs(KMs(:,2)<0,:)=[];
R4=(I_t(:,2)+PL(:,2))*(rp4);
Idx_R4=abs(KMs(:,1)-R4)==min(abs(KMs(:,1)-R4));
if sum(Idx_R4)>1
    r=abs(KMs(Idx_R4,1)-I_t)==min(abs(KMs(Idx_R4,1)-I_t));
    Idx_R4=Idx_R4(r);
end

%% Split KM plot from PL to I_t into sections
Idx_i=find(abs(KMs(:,1)-PL(2))==min(abs(KMs(:,1)-PL(2))));
Idx_f=length(KMs);
KMcut=KMs(Idx_i:Idx_f,:);
Range=max(KMcut(:,1))-min(KMcut(:,1));
Step=Ss*Range;
Strs=ones(floor(Range/Step)+2,1).*KMcut(1,1);
Idxs=ones(floor(Range/Step)+2,1);
regions=cell(floor(Range/Step)+1,1);

for i=1:(floor(Range/Step)+1)
    if i==1
        Idx1=1;
    else
        Idx1=Idx2;
    end
    Idx2=find(abs(KMcut(:,1)-(KMcut(1,1)+Step*(i)))==min(abs(KMcut(:,1)-(KMcut(1,1)+Step*(i)))));
    regions{i}=KMcut(Idx1:Idx2,[1 2]);
    Idxs(i+1)=Idx2(1);
    Strs(i+1)=KMcut(Idx2(1),1);
end

%% Calculate region where reference point resides
tstress2=KMs(:,1);
Strs2=sort([Strs;tstress2(Idx_R4)]);
region_r=find(Strs2==tstress2(Idx_R4))-1;
if length(region_r)>1
    region_r=region_r(2);
end

%% Calculate deviation from this slope and adjacent sections
% Going left from reference section
l=0;
r=0;
fit_r=polyfit(KMcut(Idxs(region_r):Idxs(region_r+1),1),KMcut(Idxs(region_r):Idxs(region_r+1),2),1);
for i=1:region_r-1
    rg=region_r-i;
    fitl=polyfit(KMcut(Idxs(rg):Idxs(rg+1),1),KMcut(Idxs(rg):Idxs(rg+1),2),1);
    if abs(fitl(1)-fit_r(1))<=abs(fit_r(1).*dev4)
        l=l+1;
    else
        break
    end    
end


%Calculate deviation from this slope in regions going right
for i=region_r+1:length(regions)
    rg=i;
    fitr=polyfit(KMcut(Idxs(rg):Idxs(rg+1),1),KMcut(Idxs(rg):Idxs(rg+1),2),1);
    if abs(fitr(1)-fit_r(1))<=abs(fit_r(1).*dev4)
        r=r+1;
    else
        break
    end    
end
regions_t4=(region_r-l):(region_r+r);

%% From Regions Found, Calculate Steady State Slope
l_regs=ones(length(regions_t4),1);
for i=1:length(regions_t4)
    l_regs(i)=length(regions{regions_t4(i)});
end
StdSt_reg4=ones(sum(l_regs)-1,2).*8;    
for i=1:length(regions_t4)
    if i==1
        st=1;
    else
        st=sum(l_regs(1:i-1));
    end
    StdSt_reg4(st:st+l_regs(i)-1,:)=regions{regions_t4(i)};
end
fit_ss=polyfit(StdSt_reg4(:,1),StdSt_reg4(:,2),1);
ss_slp=fit_ss(1);
sigma_04=fit_ss(2);
cb4=ss_slp.*-1;
KM=KM(:,1:2);

% %Uncomment to plot regions, bounds of region for slope, and KM construction
% figure
% hold on
% leg=cell(1,length(regions)+3);
% leg{end}="KM construction";
% leg{end-1}="bounds for slope calculation";
% leg{end-2}="reference point";
% for i=1:length(regions)
%     plot(regions{i}(:,1),regions{i}(:,2))
%     leg{i}="segment "+num2str(i);
% end
% scatter(KMs(Idx_R4,1),KMs(Idx_R4,2),'dm','filled')
% scatter([StdSt_reg4(1,1) StdSt_reg4(end,1)],[StdSt_reg4(1,2) StdSt_reg4(end,2)],'ok','filled')
% f=@(x) cb4*-1*x+sigma_04;
% fplot(f,[0 max(KMs(:,1))],'--r')
% legend(string(leg));
end