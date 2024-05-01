classdef SS_KMsmoothapp < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                      matlab.ui.Figure
        GridLayout                    matlab.ui.container.GridLayout
        LeftPanel                     matlab.ui.container.Panel
        FramelengthLabel              matlab.ui.control.Label
        OrderLabel                    matlab.ui.control.Label
        SavitzkyGolaySmoothSpinner_2  matlab.ui.control.Spinner
        Select1ListBox                matlab.ui.control.ListBox
        Select1Label                  matlab.ui.control.Label
        UITable                       matlab.ui.control.Table
        ExportSmoothedDataButton      matlab.ui.control.Button
        PrintSmoothButton             matlab.ui.control.Button
        SavitzkyGolaySmoothSpinner    matlab.ui.control.Spinner
        SavitzkyGolaySmoothSpinnerLabel  matlab.ui.control.Label
        LoessSmoothSpinner            matlab.ui.control.Spinner
        LoessSmoothSpinnerLabel       matlab.ui.control.Label
        MovingAvgStrainRangeSpinner   matlab.ui.control.Spinner
        MovingAvgStrainRangeLabel     matlab.ui.control.Label
        SplineSmoothSpinner           matlab.ui.control.Spinner
        SplineSmoothSpinnerLabel      matlab.ui.control.Label
        DataReductionSpinner          matlab.ui.control.Spinner
        DataReductionLabel            matlab.ui.control.Label
        DatasetNameDropDown           matlab.ui.control.DropDown
        DatasetNameDropDownLabel      matlab.ui.control.Label
        RightPanel                    matlab.ui.container.Panel
        UIAxes_2                      matlab.ui.control.UIAxes
        UIAxes                        matlab.ui.control.UIAxes
    end

    % Properties that correspond to apps with auto-reflow
    properties (Access = private)
        onePanelWidth = 576;
    end

    
    properties (Access = private)
        ESS_out % ESS smoothed data to be output
        KM_out % KM smoothed data to be output
    end
    
    methods (Access = private)
        % The following function reduces data by closest amount to that
        % specified by DataReductionSpinner
        % Inputs:
        % ESS - engineering stress strain double array with strain[mm/mm]
        % in the first column and stress[MPa] in the second column.
        % rd - double value of how much data is set to be excluded.
        % Outputs:
        % ESS_dr - engineering stress strain double array produced after
        % data reduction with strain[mm/mm] in the first column and
        % stress[MPa] in the second column.
        % cv - double representing the actual value of how much data was
        % excluded.
        function [ESS_rd,cv] = datared(~,ESS,rd)
            if rd>=50
                dr=length(ESS)*rd/100;
                ds=length(ESS)-dr;
                ESS_rd=ESS(1:round(length(ESS)/ds):end,:);
                cv=100-length(ESS_rd)/length(ESS)*100;
            else
                revs=1-rd/100;
                dr=length(ESS)*revs;
                ds=length(ESS)-dr;
                del=round(length(ESS)/ds);
                ESS_rd=ESS;
                ESS_rd(del:del:end,:)=[];
                cv=100-length(ESS_rd)/length(ESS)*100;
            end            
        end
        
        % This function is used to calculate the goodness of fit.
        % Inputs:
        % Raw - raw data used for fit
        % Fit - fit data
        % Outputs:
        % gof - a struct containing fields for sse, rsquare, adjrsquare, and rmse.
        function gof = goodoffit(~,raw,fit)
            % Initializing arrays to store goodness of fit
            gof=struct('sse',cell(1,1),'rsquare',[], 'adjrsquare',[],'rmse', [] );
            if length(raw)~=length(fit)
                [fitresult_r,~,~]=spline_fit(raw,1);
                [fitresult_f,~,~]=spline_fit(fit,1);
                raw=fitresult_r(fit(:,1));
                fit=fitresult_f(fit(:,1));
            end
            % Calculating sse
                net=fitnet(25);
                gof.sse=sse(net,raw,fit);

            % Calculating rsquare
            SSres=sum((raw-fit).^2);
            SStot=sum((raw-mean(raw)).^2);
            gof.rsquare=1-SSres/SStot;

            %Calculating adjrsquare
            gof.adjrsquare=1-SSres/SStot*(length(raw)-1)/(length(raw)-2);

            % Calculating rmse
            RMSE=rmse(raw,fit);
            RMSE(RMSE==0)=[];
            gof.rmse=RMSE;          
        end

        % This function conducts a smoothing average fit of ESS data
        % Inputs:
        % ESS - engineering stress strain double array with strain[mm/mm]
        % in the first column and stress[MPa] in the second column.
        % p - double representing the smoothing parameter
        % Outputs:
        % ESS_mva - engineering stress strain double array produced after
        % smoothing operation with strain[%] in the first column and
        % stress[MPa] in the second column.        
        function [ESS_mva] = mvgavg(~,ESS,p)
           % Producing fit
            ESS_mva=zeros(size(ESS,1),size(ESS,2));
            for i=1:length(ESS)
                if ESS(i,1)>p*100 && ESS(i,1)<max(ESS(:,1))-p*100
                    Idx_s=find(abs(ESS(:,1)-(ESS(i,1)-p*100))==min(abs(ESS(:,1)-(ESS(i,1)-p*100))));
                    Idx_f=find(abs(ESS(:,1)-(ESS(i,1)+p*100))==min(abs(ESS(:,1)-(ESS(i,1)+p*100))));
                    ESS_mva(i,:)=mean(ESS(Idx_s:Idx_f,:));
                end
            end
            ESS_mva=ESS_mva(any(ESS_mva,2),:);
        end
        
    end
    

    % Callbacks that handle component events
    methods (Access = private)

        % Drop down opening function: DatasetNameDropDown
        function DatasetNameDropDownOpening(app, ~)
            vars = evalin('base', 'whos');
            cell_array = cell(size(vars));
            for i=1:length(vars)
                cell_array{i} = vars(i).name;
            end
            app.DatasetNameDropDown.Items = cell_array;
        end

        % Button pushed function: PrintSmoothButton
        function PrintSmoothButtonPushed(app, ~)
            % Assign Input Data
            cla(app.UIAxes)
            cla(app.UIAxes_2)
            ESS_name=app.DatasetNameDropDown.Value;
            ESS=evalin('base', ESS_name);
            if size(ESS,2)~=2
                error('Input is not the proper dimensions.')
            end
            TSS=ESS_to_TSS(ESS);
            [wh_i,TSS_i]=wrkhard(TSS);
            KM=[TSS_i(:,2),wh_i];

            % Conduct Data Reduction
            if app.DataReductionSpinner.Value==0
                ESS_rd=cutend(ESS);
            else
                ESS=cutend(ESS);
                [ESS_rd,cv]=datared(app,ESS,app.DataReductionSpinner.Value);
                app.DataReductionSpinner.Value=cv;            
            end

            % Conduct Smooth Selected
            if app.Select1ListBox.ValueIndex==1
                p=app.SplineSmoothSpinner.Value;
                ESS_rd(:,1)=ESS_rd(:,1).*100;
                [~,gof,ESS_f] = spline_fit(ESS_rd,p);
                gof=rmfield(gof,'dfe');
                ESS_rd(:,1)=ESS_rd(:,1)./100;
                ESS_f(:,1)=ESS_f(:,1)./100;
                TSS_f=ESS_to_TSS(ESS_f);
                [wh,TSS_whf]=wrkhard(TSS_f);
                KM_f=[TSS_whf(:,2),wh];

            elseif app.Select1ListBox.ValueIndex==2
                p=app.MovingAvgStrainRangeSpinner.Value.*100;
                ESS_rd(:,1)=ESS_rd(:,1).*100;
                [ESS_f] = mvgavg(app,ESS_rd,p);
                gof=goodoffit(app,ESS_rd,ESS_f);
                ESS_rd(:,1)=ESS_rd(:,1)./100;
                ESS_f(:,1)=ESS_f(:,1)./100;
                TSS_f=ESS_to_TSS(ESS_f);
                [wh,TSS_whf]=wrkhard(TSS_f);
                KM_f=[TSS_whf(:,2),wh];

            elseif app.Select1ListBox.ValueIndex==3
                p=app.LoessSmoothSpinner.Value.*100;
                ESS_rd(:,1)=ESS_rd(:,1).*100;
                [ESS_f] = [ESS_rd(:,1) smooth(ESS_rd(:,1),ESS_rd(:,2),p,'rloess')];
                gof=goodoffit(app,ESS_rd,ESS_f);
                ESS_rd(:,1)=ESS_rd(:,1)./100;
                ESS_f(:,1)=ESS_f(:,1)./100;
                TSS_f=ESS_to_TSS(ESS_f);
                [wh,TSS_whf]=wrkhard(TSS_f);
                KM_f=[TSS_whf(:,2),wh];

            elseif app.Select1ListBox.ValueIndex==4
                p1=app.SavitzkyGolaySmoothSpinner.Value;
                p2=app.SavitzkyGolaySmoothSpinner_2.Value;
                ESS_rd(:,1)=ESS_rd(:,1).*100;
                x = ESS_rd(:,1);
                y=sgolayfilt(ESS_rd(:,2),p1,p2);
                ESS_f=[x y];
                gof=goodoffit(app,ESS_rd,ESS_f);
                ESS_rd(:,1)=ESS_rd(:,1)./100;
                ESS_f(:,1)=ESS_f(:,1)./100;
                TSS_f=ESS_to_TSS(ESS_f);
                [wh,TSS_whf]=wrkhard(TSS_f);
                KM_f=[TSS_whf(:,2),wh];

            else
                error('Select Smoothing Operator.')
            end

            % Plot original ESS & KM and fit curves
            hold(app.UIAxes,"on")
            scatter(app.UIAxes,ESS_rd(:,1),ESS_rd(:,2),2,'k','filled')
            plot(app.UIAxes,ESS_f(:,1),ESS_f(:,2),"LineWidth",1,"Color",'r')
            legend(app.UIAxes,'ESS Reduce','ESS fit')
            
            hold(app.UIAxes_2,"on")
            scatter(app.UIAxes_2,KM(:,1),KM(:,2),2,'k','filled')
            plot(app.UIAxes_2,KM_f(:,1),KM_f(:,2),"LineWidth",1,"Color",'r')
            legend(app.UIAxes_2,'KM Reduced','KM fit')

            % Fill in Table
            app.UITable.Data=struct2cell(gof(1,:))';

            % Set outputs
            app.ESS_out=ESS_f;
            app.KM_out=KM_f;




        end

        % Button pushed function: ExportSmoothedDataButton
        function ExportSmoothedDataButtonPushed(app, ~)
            if exist(app.ESS_out,'var')==1
                assignin('base','ESS_out',app.ESS_out)
                assignin('base','KM_out',app.KM_out)
            end
        end

        % Changes arrangement of the app based on UIFigure width
        function updateAppLayout(app, ~)
            currentFigureWidth = app.UIFigure.Position(3);
            if(currentFigureWidth <= app.onePanelWidth)
                % Change to a 2x1 grid
                app.GridLayout.RowHeight = {480, 480};
                app.GridLayout.ColumnWidth = {'1x'};
                app.RightPanel.Layout.Row = 2;
                app.RightPanel.Layout.Column = 1;
            else
                % Change to a 1x2 grid
                app.GridLayout.RowHeight = {'1x'};
                app.GridLayout.ColumnWidth = {223, '1x'};
                app.RightPanel.Layout.Row = 1;
                app.RightPanel.Layout.Column = 2;
            end
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.AutoResizeChildren = 'off';
            app.UIFigure.Position = [100 100 640 480];
            app.UIFigure.Name = 'MATLAB App';
            app.UIFigure.SizeChangedFcn = createCallbackFcn(app, @updateAppLayout, true);

            % Create GridLayout
            app.GridLayout = uigridlayout(app.UIFigure);
            app.GridLayout.ColumnWidth = {223, '1x'};
            app.GridLayout.RowHeight = {'1x'};
            app.GridLayout.ColumnSpacing = 0;
            app.GridLayout.RowSpacing = 0;
            app.GridLayout.Padding = [0 0 0 0];
            app.GridLayout.Scrollable = 'on';

            % Create LeftPanel
            app.LeftPanel = uipanel(app.GridLayout);
            app.LeftPanel.Layout.Row = 1;
            app.LeftPanel.Layout.Column = 1;

            % Create DatasetNameDropDownLabel
            app.DatasetNameDropDownLabel = uilabel(app.LeftPanel);
            app.DatasetNameDropDownLabel.HorizontalAlignment = 'right';
            app.DatasetNameDropDownLabel.Position = [18 445 82 22];
            app.DatasetNameDropDownLabel.Text = 'Dataset Name';

            % Create DatasetNameDropDown
            app.DatasetNameDropDown = uidropdown(app.LeftPanel);
            app.DatasetNameDropDown.DropDownOpeningFcn = createCallbackFcn(app, @DatasetNameDropDownOpening, true);
            app.DatasetNameDropDown.Position = [115 445 70 26];

            % Create DataReductionLabel
            app.DataReductionLabel = uilabel(app.LeftPanel);
            app.DataReductionLabel.HorizontalAlignment = 'center';
            app.DataReductionLabel.Position = [5 406 80 30];
            app.DataReductionLabel.Text = {'Data '; 'Reduction [%]'};

            % Create DataReductionSpinner
            app.DataReductionSpinner = uispinner(app.LeftPanel);
            app.DataReductionSpinner.Step = 5;
            app.DataReductionSpinner.Limits = [0 95];
            app.DataReductionSpinner.Position = [100 410 100 22];

            % Create SplineSmoothSpinnerLabel
            app.SplineSmoothSpinnerLabel = uilabel(app.LeftPanel);
            app.SplineSmoothSpinnerLabel.HorizontalAlignment = 'right';
            app.SplineSmoothSpinnerLabel.Position = [6 277 83 22];
            app.SplineSmoothSpinnerLabel.Text = 'Spline Smooth';

            % Create SplineSmoothSpinner
            app.SplineSmoothSpinner = uispinner(app.LeftPanel);
            app.SplineSmoothSpinner.Step = 0.0001;
            app.SplineSmoothSpinner.Limits = [0 1];
            app.SplineSmoothSpinner.Position = [117 277 100 22];
            app.SplineSmoothSpinner.Value = 1;

            % Create MovingAvgStrainRangeLabel
            app.MovingAvgStrainRangeLabel = uilabel(app.LeftPanel);
            app.MovingAvgStrainRangeLabel.HorizontalAlignment = 'right';
            app.MovingAvgStrainRangeLabel.Position = [14 237 75 30];
            app.MovingAvgStrainRangeLabel.Text = {'Moving Avg. '; 'Strain Range'};

            % Create MovingAvgStrainRangeSpinner
            app.MovingAvgStrainRangeSpinner = uispinner(app.LeftPanel);
            app.MovingAvgStrainRangeSpinner.Step = 0.0001;
            app.MovingAvgStrainRangeSpinner.Limits = [0 1];
            app.MovingAvgStrainRangeSpinner.Position = [117 241 100 22];
            app.MovingAvgStrainRangeSpinner.Value = 0.0002;

            % Create LoessSmoothSpinnerLabel
            app.LoessSmoothSpinnerLabel = uilabel(app.LeftPanel);
            app.LoessSmoothSpinnerLabel.HorizontalAlignment = 'right';
            app.LoessSmoothSpinnerLabel.Position = [7 207 82 22];
            app.LoessSmoothSpinnerLabel.Text = 'Loess Smooth';

            % Create LoessSmoothSpinner
            app.LoessSmoothSpinner = uispinner(app.LeftPanel);
            app.LoessSmoothSpinner.Step = 0.0001;
            app.LoessSmoothSpinner.Limits = [0 1];
            app.LoessSmoothSpinner.Position = [117 207 100 22];
            app.LoessSmoothSpinner.Value = 0.0001;

            % Create SavitzkyGolaySmoothSpinnerLabel
            app.SavitzkyGolaySmoothSpinnerLabel = uilabel(app.LeftPanel);
            app.SavitzkyGolaySmoothSpinnerLabel.HorizontalAlignment = 'center';
            app.SavitzkyGolaySmoothSpinnerLabel.Position = [7 154 88 30];
            app.SavitzkyGolaySmoothSpinnerLabel.Text = {'Savitzky-Golay '; 'Smooth'};

            % Create SavitzkyGolaySmoothSpinner
            app.SavitzkyGolaySmoothSpinner = uispinner(app.LeftPanel);
            app.SavitzkyGolaySmoothSpinner.Limits = [1 6];
            app.SavitzkyGolaySmoothSpinner.RoundFractionalValues = 'on';
            app.SavitzkyGolaySmoothSpinner.ValueDisplayFormat = '%11g';
            app.SavitzkyGolaySmoothSpinner.Position = [169 169 46 22];
            app.SavitzkyGolaySmoothSpinner.Value = 1;

            % Create PrintSmoothButton
            app.PrintSmoothButton = uibutton(app.LeftPanel, 'push');
            app.PrintSmoothButton.ButtonPushedFcn = createCallbackFcn(app, @PrintSmoothButtonPushed, true);
            app.PrintSmoothButton.FontSize = 14;
            app.PrintSmoothButton.Position = [6 94 78 42];
            app.PrintSmoothButton.Text = {'Print '; 'Smooth'};

            % Create ExportSmoothedDataButton
            app.ExportSmoothedDataButton = uibutton(app.LeftPanel, 'push');
            app.ExportSmoothedDataButton.ButtonPushedFcn = createCallbackFcn(app, @ExportSmoothedDataButtonPushed, true);
            app.ExportSmoothedDataButton.FontSize = 14;
            app.ExportSmoothedDataButton.Position = [90 94 113 42];
            app.ExportSmoothedDataButton.Text = {'Export '; 'Smoothed Data'};

            % Create UITable
            app.UITable = uitable(app.LeftPanel);
            app.UITable.ColumnName = {'SSE'; 'R2'; 'ADR2'; 'RMSE'};
            app.UITable.RowName = {};
            app.UITable.Position = [0 6 206 82];

            % Create Select1Label
            app.Select1Label = uilabel(app.LeftPanel);
            app.Select1Label.HorizontalAlignment = 'right';
            app.Select1Label.Position = [14 369 52 22];
            app.Select1Label.Text = 'Select 1:';

            % Create Select1ListBox
            app.Select1ListBox = uilistbox(app.LeftPanel);
            app.Select1ListBox.Items = {'Spline Fit', 'Moving Average', 'Loess Fit', 'Savitzky-Golay'};
            app.Select1ListBox.Position = [81 310 119 83];
            app.Select1ListBox.Value = 'Spline Fit';

            % Create SavitzkyGolaySmoothSpinner_2
            app.SavitzkyGolaySmoothSpinner_2 = uispinner(app.LeftPanel);
            app.SavitzkyGolaySmoothSpinner_2.Step = 2;
            app.SavitzkyGolaySmoothSpinner_2.Limits = [3 Inf];
            app.SavitzkyGolaySmoothSpinner_2.RoundFractionalValues = 'on';
            app.SavitzkyGolaySmoothSpinner_2.ValueDisplayFormat = '%11g';
            app.SavitzkyGolaySmoothSpinner_2.Position = [170 142 46 22];
            app.SavitzkyGolaySmoothSpinner_2.Value = 3;

            % Create OrderLabel
            app.OrderLabel = uilabel(app.LeftPanel);
            app.OrderLabel.Position = [115 168 39 22];
            app.OrderLabel.Text = 'Order:';

            % Create FramelengthLabel
            app.FramelengthLabel = uilabel(app.LeftPanel);
            app.FramelengthLabel.Position = [96 142 76 22];
            app.FramelengthLabel.Text = 'Framelength:';

            % Create RightPanel
            app.RightPanel = uipanel(app.GridLayout);
            app.RightPanel.Layout.Row = 1;
            app.RightPanel.Layout.Column = 2;

            % Create UIAxes
            app.UIAxes = uiaxes(app.RightPanel);
            title(app.UIAxes, 'ESS Curve')
            xlabel(app.UIAxes, 'Engineering Strain [mm/mm]')
            ylabel(app.UIAxes, 'Engineering Stress [MPa]')
            zlabel(app.UIAxes, 'Z')
            app.UIAxes.Position = [49 262 323 199];

            % Create UIAxes_2
            app.UIAxes_2 = uiaxes(app.RightPanel);
            title(app.UIAxes_2, 'KM Curve')
            xlabel(app.UIAxes_2, 'True Stress [MPa]')
            ylabel(app.UIAxes_2, 'Work Hardening Rate [MPa]')
            zlabel(app.UIAxes_2, 'Z')
            app.UIAxes_2.Position = [43 32 334 206];

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = SS_KMsmoothapp

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.UIFigure)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.UIFigure)
        end
    end
end