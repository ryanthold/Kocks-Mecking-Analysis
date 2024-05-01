%%This function serves to create a smoothing spline fit of stress-strain data
%% Inputs 
% SS: 2-by-x double array of stress[MPa]-strain[%] data with strain in the first column and stress in the second column
% p: double with a value between 0 and 1 representing the smoothing parameter

% Note: Input strain units will require adjusted smoothing parameter.  It
% is suggested to use % strain.
%% Outputs:
% fitresult: cfit smoothing spline with p coefficient structure
% gof: struct array with fields sse, rsquare, dfe, adjrsquare, and rmse
% showing these various goodness of fit values
% SS_sp: 2-by-x double array with spline fit curve using strain values of
% input SS.

function [fitresult,gof,SS_sp] = spline_fit(SS,p)
%% Initialization
% Check Inputs
if class(SS)~="double"
    error('stress-strain (SS) input is of the wrong data type.')
elseif class(p)~="double"
    error('smoothing parameter (p) input is of the wrong data type.')
elseif size(SS,2)~=2
    error('stress-strain (SS) input is of the wrong dimensions.')
elseif size(p)~=1
    error('smoothing parameter (p) input is of the wrong dimensions.')
elseif p>1 || p<0
    error('smoothing parameter (p) must be a value between 0 and 1.')
end

% Setting stress and strain values to individual variables.
SS(any(isinf(SS),1),:)=[];
SS(any(isnan(SS),1),:)=[];
SS(any(isinf(SS),2),:)=[];
SS(any(isnan(SS),2),:)=[];
strain=SS(:,1);
stress=SS(:,2);

% Initialize arrays to store fits and goodness-of-fit.
gof = struct( 'sse', cell( 1, 1 ), ...
    'rsquare', [], 'dfe', [], 'adjrsquare', [], 'rmse', [] );

%% Fit
[xData, yData] = prepareCurveData( strain, stress );

% Set up fittype and options.
ft = fittype( 'smoothingspline' );
opts = fitoptions( 'Method', 'SmoothingSpline' );
opts.SmoothingParam = p;

% Fit model to strain values
[fitresult, gof(1)] = fit( xData, yData, ft, opts );
stress=fitresult(strain);
SS_sp=[strain stress];

% % Uncomment below to fit spline fit data on top of input data
% figure
% hold on
% plot(SS(:,1),SS(:,2))
% plot(SS_sp(:,1),SS_sp(:,2))
% xlabel('Strain [%]')
% ylabel('Stress [MPa]')
% legend('SS','SS_s_p')
end