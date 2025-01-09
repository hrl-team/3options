function [Nbar,Nsub] = bars_datamodel(DataCell,Model_DataCell,Colors,Yinf,Ysup)

% Sophie Bavard - December 2018
% Creates a violin plot with mean, error bars, confidence interval, kernel density.
% Warning: the function can accept any number of arguments > 9.
% After the Title, LabelX, LabelY : varargin for bar names under X-axis

% transforms the Data matrix into cell format if needed
if iscell(DataCell)==0
    DataCell = num2cell(DataCell,2);
end
if iscell(Model_DataCell)==0
    Model_DataCell = num2cell(Model_DataCell', 1)';
end

% number of factors/groups/conditions
Nbar = size(DataCell,1);
% number of models
Nmod = size(Model_DataCell,3);
% bar size
Wbar = 0.75;
% model spacing
if Nmod>1
    space = zscore(1:Nmod)*Wbar/Nmod;
    shape = ['o','o','o','o'];
else
    space = 0;
    shape = 'o';
end

% confidence interval
ConfInter = 0.95;

% color of the box + error bar
trace = [0.5 0.5 0.5];

for n = 1:Nbar
    
    clear DataMatrix
    clear jitter jitterstrength
    DataMatrix = DataCell{n,:}';
    
    % number of subjects
    Nsub = length(DataMatrix(~isnan(DataMatrix)));
    
    curve = nanmean(DataMatrix);
    sem   = nanstd(DataMatrix')'/sqrt(Nsub);
    std   = nanstd(DataMatrix')';
    conf  = tinv(1 - 0.5*(1-ConfInter),Nsub);

    % COLORED BARS
    bar(n,curve,...
        'FaceColor',Colors(n,:),...
        'EdgeColor',Colors(n,:),...
        'BarWidth',0.75,...
        'LineWidth',1, ...
        'FaceAlpha',0.7);
    hold on
        
    % S.T.D. RECTANGLE
%     rectangle('Position',[n-Wbar/2, curve - std, Wbar, std*2],...
%         'EdgeColor','none',...
%         'FaceColor',[Colors(n,:) 0.2],...
%         'LineWidth',1);
%     hold on
    
    % CONFIDENCE INTERVAL RECTANGLE
%     if Nsub>1
%         rectangle('Position',[n-Wbar/2, curve - sem*conf, Wbar, sem*conf*2],...
%             'EdgeColor','none',...
%             'FaceColor',[Colors(n,:) 0.5],...
%             'LineWidth',1);
%         hold on
%     end
    
    % S.E.M. RECTANGLE
%     rectangle('Position',[n-Wbar/2, curve - sem, Wbar, sem*2],...
%         'EdgeColor','none',...
%         'FaceColor',Colors(n,:),...
%         'LineWidth',1);
%     hold on
    
    % MEAN HORIZONTAL BAR
%     xMean = [n - Wbar/2 ; n + Wbar/2];
%     yMean = [curve; curve];
%     plot(xMean,yMean,'-','LineWidth',1,'Color','k');
%     hold on 

    % ERROR BARS
%     errorbar(n,curve,sem,...
%         'Color',[0 0 0],...
%         'LineStyle','none',...
%         'LineWidth',1);
%     hold on
    
    for m=1:Nmod
        
        if size(Model_DataCell,1)>1
            Model_Matrix = Model_DataCell(:,:,m);
        else
            Model_Matrix = Model_DataCell{:,:,m};
        end
        if iscell(Model_Matrix)==0
            Model_Matrix = num2cell(Model_Matrix', 1)';
        end
        
        % PLOT THE MODEL SIMULATIONS
        if Nmod > 1
            errorbar(n+space(m),nanmean(Model_Matrix{n,:}'),conf*nanstd(Model_Matrix{n,:}')'/sqrt(Nsub),...
                'Color',[0 0 0],...
                'LineWidth',1,...
                'LineStyle','none',...
                'Marker',shape(m),...
                'MarkerFaceColor',1-[(m-1)/(Nmod-1) (m-1)/(Nmod-1) (m-1)/(Nmod-1)],...
                'MarkerSize',15/Nmod);
            hold on
        elseif Nmod == 1
            errorbar(n+space(m),nanmean(Model_Matrix{n,:}'),conf*nanstd(Model_Matrix{n,:}')'/sqrt(Nsub),...
                'Color',[0 0 0],...
                'LineWidth',1,...
                'LineStyle','none',...
                'Marker',shape(m),...
                'MarkerFaceColor','k',...
                'MarkerSize',7.5);
            hold on
        end
        
    end
end


ylim([Yinf Ysup]);
xlim([0 Nbar+1]);











