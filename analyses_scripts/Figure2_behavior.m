clear


%% Load data

whichexpe = listdlg('SelectionMode','single','ListString',{'expe1','expe2'});

if whichexpe == 1
    load data_fig_expe1
elseif whichexpe == 2
    load data_fig_expe2
end

expe_a = 1:50;
expe_b = 51:100;
expe_c = 101:150;


%% Choose colors for figures

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


%% EXPE A (no forced trials)

figure;

set(gcf,'Position',[100 200 800 400])

subplot(1,3,1)
violinplotSB([mean(cond4(:,expe_a));mean(cond3(:,expe_a));mean(cond2(:,expe_a));mean(cond1(:,expe_a))],Colors([4 3 2 1],:),0,1);
yline(1/2,'k:'); yline(1/3,'k:');
set(gca, 'Box', 'off') ;

subplot(1,3,2:3)
violinplotSB(rating([10 9 8 7 6 5 4 3 2 1],expe_a),Colors2([10 9 8 7 6 5 4 3 2 1],:),0,1);
yline(1/2,'k:'); yline(1/3,'k:');
set(gca, 'Box', 'off') ;


%% EXPE B (forced trials / complete feedback)

figure;

set(gcf,'Position',[100 200 800 400])

subplot(1,3,1)
violinplotSB([mean(cond4(:,expe_b));mean(cond3(:,expe_b));mean(cond2(:,expe_b));mean(cond1(:,expe_b))],Colors([4 3 2 1],:),0,1);
yline(1/2,'k:'); yline(1/3,'k:');
set(gca, 'Box', 'off') ;

subplot(1,3,2:3)
violinplotSB(rating([10 9 8 7 6 5 4 3 2 1],expe_b),Colors2([10 9 8 7 6 5 4 3 2 1],:),0,1);
yline(1/2,'k:'); yline(1/3,'k:');
set(gca, 'Box', 'off') ;


%% EXPE C (forced trials / partial feedback)

figure;

set(gcf,'Position',[100 200 800 400])

subplot(1,3,1)
violinplotSB([mean(cond4(:,expe_c));mean(cond3(:,expe_c));mean(cond2(:,expe_c));mean(cond1(:,expe_c))],Colors([4 3 2 1],:),0,1);
yline(1/2,'k:'); yline(1/3,'k:');
set(gca, 'Box', 'off') ;

subplot(1,3,2:3)
violinplotSB(rating([10 9 8 7 6 5 4 3 2 1],expe_c),Colors2([10 9 8 7 6 5 4 3 2 1],:),0,1);
yline(1/2,'k:'); yline(1/3,'k:');
set(gca, 'Box', 'off') ;