
clear


%% Load data

data =  1;
model1= 1;
model2= 3;
model3= 5;
model4= 3;

% Expe 1-2
% 1 Q-LEARNING
% 2 DIVISIVE
% 3 RANGE
% 4 RANGE uti
% 5 DIVISIVE uti
% 6 WEBB 2020
% 7 HYBRID D-R
% 8 HYBRID Q-R
% 9 HABIT

% Expe 3
% 1 Q-LEARNING
% 2 RANGE
% 3 RANGE 1 omega 2 alphas
% 4 RANGE 1 omega 1 alpha
% 5 RANGE 2 omega 2 alpha
% 6 RANGE 2 omega 1 alpha

whichexpe = listdlg('ListString',{'expe1','expe2','expe3'});

sim_data = load(strcat('simulation_model',num2str(data),'_expe',num2str(whichexpe),'_learning'));
sim_mod1 = load(strcat('simulation_model',num2str(model1),'_expe',num2str(whichexpe),'_learning'));
sim_mod2 = load(strcat('simulation_model',num2str(model2),'_expe',num2str(whichexpe),'_learning'));
sim_mod3 = load(strcat('simulation_model',num2str(model3),'_expe',num2str(whichexpe),'_learning'));
sim_mod4 = load(strcat('simulation_model',num2str(model4),'_expe',num2str(whichexpe),'_learning'));


%% Colors for figure

if whichexpe < 3

    order = 10:-1:1;

    Colors(1,:) = [135 62 35] /255;
    Colors(2,:) = [21 76 121] /255;
    Colors(3,:) = [226 135 67]/255;
    Colors(4,:) = [30 129 176] /255;

    Colors2(1,:) = Colors(1,:);
    Colors2(2,:) = Colors(1,:);
    Colors2(3,:) = Colors(1,:);
    Colors2(4,:) = Colors(2,:);
    Colors2(5,:) = Colors(2,:);
    Colors2(6,:) = Colors(3,:);
    Colors2(7,:) = Colors(3,:);
    Colors2(8,:) = Colors(3,:);
    Colors2(9,:) = Colors(4,:);
    Colors2(10,:) = Colors(4,:);

else

    order = 12:-1:1;

    Colors(1,:) = [29 81 90]  /255;
    Colors(2,:) = [97 179 159] /255;
    Colors(3,:) = [159 198 148] /255;
    Colors(4,:) = [224 133 89] /255;

    Colors2(1,:) = Colors(1,:);
    Colors2(2,:) = Colors(1,:);
    Colors2(3,:) = Colors(1,:);
    Colors2(4,:) = Colors(2,:);
    Colors2(5,:) = Colors(2,:);
    Colors2(6,:) = Colors(2,:);
    Colors2(7,:) = Colors(3,:);
    Colors2(8,:) = Colors(3,:);
    Colors2(9,:) = Colors(3,:);
    Colors2(10,:) = Colors(4,:);
    Colors2(11,:) = Colors(4,:);
    Colors2(12,:) = Colors(4,:);

end

%% Reorganize simulation variables

% Learning accuracy
perfLT_m{:,:,1} = sim_mod1.perf_m;
perfLT_m{:,:,2} = sim_mod2.perf_m;
perfLT_m{:,:,3} = sim_mod3.perf_m;
perfLT_m{:,:,4} = sim_mod4.perf_m;

% Transfer accuracy
perfPT_m{:,:,1} = sim_mod1.perfPT_m;
perfPT_m{:,:,2} = sim_mod2.perfPT_m;
perfPT_m{:,:,3} = sim_mod3.perfPT_m;
perfPT_m{:,:,4} = sim_mod4.perfPT_m;

% Transfer ratings
rating_m{:,:,1} = sim_mod1.rating_m;
rating_m{:,:,2} = sim_mod2.rating_m;
rating_m{:,:,3} = sim_mod3.rating_m;
rating_m{:,:,4} = sim_mod4.rating_m;

%%

figure;
bars_datamodel(rating_m{:,:,1}(order,:),sim_data.rating(order,:),Colors2(order,:),0,1);
yline(0.5,'k:');
yline(1/3,'k:');
set(gca, 'Box', 'off') ;

figure;
bars_datamodel(rating_m{:,:,2}(order,:),sim_data.rating(order,:),Colors2(order,:),0,1);
yline(0.5,'k:');
yline(1/3,'k:');
set(gca, 'Box', 'off') ;

figure;
bars_datamodel(rating_m{:,:,3}(order,:),sim_data.rating(order,:),Colors2(order,:),0,1);
yline(0.5,'k:');
yline(1/3,'k:');
set(gca, 'Box', 'off') ;

figure;
bars_datamodel(rating_m{:,:,4}(order,:),sim_data.rating(order,:),Colors2(order,:),0,1);
yline(0.5,'k:');
yline(1/3,'k:');
set(gca, 'Box', 'off') ;
