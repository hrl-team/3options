clear

%% Load data

whichexpe = listdlg('SelectionMode','single','ListString',{'expe1','expe2','expe3'});

if whichexpe == 1
    simus = load('simulation_model4_expe1_learning.mat');
elseif whichexpe == 2
    simus = load('simulation_model4_expe2_learning.mat');
elseif whichexpe == 3
    simus = load('simulation_model5_expe3_learning.mat');
end

uti = simus.parameters(:,4);

if whichexpe == 3
    utiU = simus.parameters(:,5);
end

x = 0:0.001:1;

Colors = [179 0 242]/255;


%% Plot curves

for sub = simus.subjects
    curve(:,sub) = x.^uti(sub)';
end

mcb = mean(curve');

figure;

for sub = simus.subjects
    plot(x.*size(curve,1), x.^uti(sub),'Color',[0 0 0 .1]);
    hold on
end

SurfaceCurvePlot(curve,NaN,Colors,1,.3,0,1,12,'','','');  % [0 182 235]/255

hold on

line(x.*size(curve,1),x,'Color','k','LineStyle','--');

line([size(curve,1)/2 size(curve,1)/2],[0 mcb(round(size(curve,1)/2))],'Color','k','LineStyle',':');
line([0 size(curve,1)/2],[mcb(round(size(curve,1)/2)) mcb(round(size(curve,1)/2))],'Color','k','LineStyle',':');

xticks([0 size(curve,1)/2 size(curve,1)]);
xticklabels({'0','0.5','1'});

yticks([0 mcb(round(size(curve,1)/2)) 0.5 1]);
yticklabels({'0',num2str(round(mcb(round(size(curve,1)/2)),2)),'0.5','1'});

pbaspect([1 1 1]);

%% Unchosen

if whichexpe == 3

    for sub = simus.subjects
        curveU(:,sub) = x.^utiU(sub)';
    end

    mcs = mean(curveU');

    figure;

    for sub = data.subjects
        plot(x.*size(curveU,1), x.^utiU(sub),'Color',[0 0 0 .1]);
        hold on
    end

    SurfaceCurvePlot(curveU,NaN,Colors,1,.3,0,1,12,'','','');  % [0 182 235]/255

    hold on

    line(x.*size(curveU,1),x,'Color','k','LineStyle','--');

    line([size(curveU,1)/2 size(curveU,1)/2],[0 mcs(round(size(curveU,1)/2))],'Color','k','LineStyle',':');
    line([0 size(curveU,1)/2],[mcs(round(size(curveU,1)/2)) mcs(round(size(curveU,1)/2))],'Color','k','LineStyle',':');

    xticks([0 size(curveU,1)/2 size(curveU,1)]);
    xticklabels({'0','0.5','1'});

    yticks(sort([0 mcs(round(size(curveU,1)/2)) 0.5 1]));
    yticklabels({'0',round(mcs(round(size(curveU,1)/2)),2),'0.5','1'});

    pbaspect([1 1 1]);

end
