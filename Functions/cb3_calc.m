%% This function is used to determine the values of cb3 through slope construction on the KM curve
%% Inputs
% wh: 1-by-(x-2) double array of work hardening rate [MPa]
% TSS: 2-by-x double array of true stress[MPa]-strain[mm/mm] data with strain in the first column and stress in the second column
% PL: 1-by-2 double array representing the proportional limit with strain [mm/mm] in the first column and stress[MPa] in the second column
% off_02: 1-by-2 double array representing the 0.2% offset yield point with strain[mm/mm] in the first column and stress[MPa] in the second column
% I_t: 1-by-2 double array representing the instability with true strain [mm/mm] in the first column and true stress[MPa] in the second column
% dev3: a double between 0 and 1 representing the allowable deviation in slope between sections to be included in the overall slope measurement of region 3
% Ss: a double between 0 and 1 representing the percentage of total plastic stress to be used for section size
% p3: a double between 0 and 1 representing the spline smoothing parameter for the curve of the derivative of the KM curve
%% Outputs
% cb3: a double with the value of the slope of region 3
% KM: a 2-by-x double with stress[MPa] in column 1 and work hardening[MPa] in column 2
%regnum: number of sections used
%regions_t: outputs the section numbers of region 3 used for slope calculation
%StdSt_reg3: outputs points of region 3 used for slope calculation
%sigma_03: the y-intercept of the cb3 line construct

function[cb3,KM,regions_t3,StdSt_reg3,sigma_03]=cb3_calc(wh,TSS,PL,off_02,I_t,dev3,Ss,p3)
%% Check Inputs
if class(wh)~="double"
    error('work hardening rate (wh) input is of the wrong data type.')
elseif class(TSS)~="double"
    error('true stress-strain (TSS) input is of the wrong data type.')
elseif class(PL)~="double"
    error('proportional limit (PL) input is of the wrong data type.')
elseif class(off_02)~="double"
    error('0.2% offset (off_02) input is of the wrong data type.')
elseif class(I_t)~="double"
    error('Instability (I_t) input is of the wrong data type.')
elseif class(dev3)~="double"
    error('deviation value for region 3 (dev3) input is of the wrong data type.')
elseif class(Ss)~="double"
    error('section size value (Ss) input is of the wrong data type.')
elseif class(p3)~="double"
    error('spline smoothing parameter (p3) input is of the wrong data type.')
elseif size(wh,2)~=1
    error('work hardening rate (wh) input is of the wrong dimensions.')
elseif size(TSS,2)~=2
    error('true stress-strain (TSS) input is of the wrong dimensions.')
elseif size(PL,2)~=2
    error('proportional limit (PL) input is of the wrong dimensions.')
elseif size(off_02,2)~=2
    error('0.2% offset (off_02) input is of the wrong dimensions.')
elseif size(I_t,2)~=2
    error('instability (I_t) input is of the wrong dimensions.')
elseif size(dev3)~=1
    error('dev3 input value is of the wrong dimensions.')
elseif size(Ss)~=1
    error('section size (Ss) input value is of the wrong dimensions.')
elseif size(p3)~=1
    error('spline smoothing parameter (p3) input value is of the wrong dimensions.') 
elseif Ss<0 || Ss>1
    error('section size (Ss) input value must be between 0 and 1')
elseif p3<0 || p3>1
    error('spline smoothing parameter (p3) input value must be between 0 and 1')
end
%% Determine reference point for region 3
off_02_t=ESS_to_TSS(off_02);
% take derivative of KM plot
KM=[TSS(:,2) wh];
KM(any(isnan(KM),2),:)=[];
[dwh,KM1]=wrkhard(KM);
[~,~,dKM]=spline_fit([KM1(:,1)./100,dwh],p3);
dKM(:,1)=dKM(:,1).*100;
dKM=sortrows(dKM);
% find local minima of derivative plot and find closest to midpoint between local PL and 0.2% offset
[~,locs]=findpeaks(-1.*dKM(:,2));
R3_e=(PL(2)+off_02_t(2))/2;
St_R3=dKM(locs(abs(dKM(locs,1)-R3_e)==min(abs(dKM(locs,1)-R3_e))),:);
Idx_R3=find(abs(KM(:,1)-St_R3(1))==min(abs(KM(:,1)-St_R3(1))));

% % Uncomment to plot KM derivative and point of local min
% figure
% hold on
% plot(dKM(:,1),dKM(:,2))
% scatter(dKM(locs,1),-1.*peaks,'or','filled')
% scatter(St_R3(1),St_R3(2),'dg','filled')
% legend('deverivative of KM vs True Stress','Local Minima','Minima corresponding to inflection point in region 3')
% figure
% hold on
% plot(KM1(:,1),dwh)
% scatter(St_R3(1),St_R3(2),'ok','filled')
% xlabel('Stress [MPa]','FontSize',15)
% ylabel('d^2\sigma / d\epsilon^2 [MPa]')

%% Split KM plot from half PL to midpoint between PL and I into sections
Idx_PL=find(abs(KM(:,1)-PL(2))==min(abs(KM(:,1)-PL(2))));
Idx_t=length(KM);
hPL=PL(2)/2;
Idx_i=find(abs(KM(:,1)-hPL)==min(abs(KM(:,1)-hPL)));
hPLI=(PL(2)+I_t(2))/2;
Idx_f=find(abs(KM(:,1)-hPLI)==min(abs(KM(:,1)-hPLI)));
KM4=KM(Idx_PL:Idx_t,:);
KMcut=KM(Idx_i:Idx_f,:);
Range=max(KM4(:,1))-min(KM4(:,1));
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
tstress2=KM(:,1);
Strs2=sort([Strs;tstress2(Idx_R3)]);
region_r=find(Strs2==tstress2(Idx_R3))-1;
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
    if abs(fitl(1)-fit_r(1))<=abs(fit_r(1).*dev3)
        l=l+1;
    else
        break
    end    
end


%Calculate deviation from this slope in regions going right
for i=region_r+1:length(regions)
    rg=i;
    fitr=polyfit(KMcut(Idxs(rg):Idxs(rg+1),1),KMcut(Idxs(rg):Idxs(rg+1),2),1);
    if abs(fitr(1)-fit_r(1))<=abs(fit_r(1).*dev3)
        r=r+1;
    else
        break
    end    
end
regions_t3=(region_r-l):(region_r+r);

%% From Regions Found, Calculate Steady State Slope
l_regs=ones(length(regions_t3),1);
for i=1:length(regions_t3)
    l_regs(i)=length(regions{regions_t3(i)});
end
StdSt_reg3=ones(sum(l_regs)-1,2).*8;   
for i=1:length(regions_t3)
    if i==1
        st=1;
    else
        st=sum(l_regs(1:i-1));
    end
    StdSt_reg3(st:st+l_regs(i)-1,:)=regions{regions_t3(i)};
end
fit_ss=polyfit(StdSt_reg3(:,1),StdSt_reg3(:,2),1);
ss_slp=fit_ss(1);
sigma_03=fit_ss(2);
cb3=ss_slp.*-1;

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
% scatter(KM(Idx_R3,1),KM(Idx_R3,2),'dm','filled')
% scatter([StdSt_reg3(1,1) StdSt_reg3(end,1)],[StdSt_reg3(1,2) StdSt_reg3(end,2)],'ok','filled')
% f=@(x) cb3*-1*x+sigma_03;
% fplot(f,[PL(2)/2 max(KM(:,1))/2],'--r')
% legend(string(leg));

end