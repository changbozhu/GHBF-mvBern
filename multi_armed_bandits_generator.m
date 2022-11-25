close all;
clear all;
clc;

A=[0.5 0.4, 0.6, 0.2, 0.5, 0.9, 0.1, 0.9, 0.1, 0.5, 0.9, 0.1, 0.5, 0.1, 0.7, 0.5, 0.3];
B=[0.5 0.6, 0.4, 0.8, 0.5, 0.1, 0.9, 0.1, 0.9, 0.5, 0.9, 0.1, 0.5, 0.1, 0.7, 0.5, 0.3];
NUM_TRIALS_PER_STAGE=15;
NUM_STAGES=length(A) ;
muStepA = [A; NUM_TRIALS_PER_STAGE*ones( 1,NUM_STAGES) ];
muStepB = [B; NUM_TRIALS_PER_STAGE*ones( 1,NUM_STAGES ) ];


muA=[];
muB=[];
for idx = 1: size(muStepA,2)
    muA=[muA;  muStepA(1,idx)*ones(muStepA(2,idx),1) ];
end
for idx = 1: size(muStepB,2)
    muB=[muB;  muStepB(1,idx)*ones(muStepB(2,idx),1) ];
end
fz  =14;
tfz =14;
afz =14;
outerpos = [1,1,1280, 640];
set(0,'defaultfigurecolor','w')
fig=figure(...
    'OuterPosition', outerpos, ...
    'Visible','on' );
subplot(211)
plot(1:length(muA), muA, 'black-', 'LineWidth', 2);
hold on;
set(gca,  'fontsize', afz)
xlim([1 length(muA)])
ylim([-0.05, 1.05])
xlabel('Trial.',  'fontsize', fz)
ylabel('$P(x_0^{(1)}=1)$', 'interpreter','latex' ,  'fontsize', fz)
grid on;
%box off;



subplot(212)
plot(1:length(muB), muB, 'black-', 'LineWidth', 2);
hold on;
xlim([1 length(muB)])
ylim([-0.05, 1.05])
set(gca,  'fontsize', afz)
xlabel('Trial.',  'fontsize', fz)
ylabel('$P(x_0^{(2)}=1)$',  'interpreter','latex' , 'fontsize', fz)
ylim([-0.05 1.05 ])
grid on;
%box off;
print('-depsc', '-painters', '-r300', ['.' filesep 'MAB_MU.eps'])
print('-dpng', '-painters', '-r300',  ['.' filesep 'MAB_MU.png'])


for i=1:100
    savepath=['.' filesep 'env' num2str(i)];
    mkdir(savepath);
    multiArmedBanditsGenerator(savepath, muStepA, muStepB);
end

function [uA, uB, muA, muB]=multiArmedBanditsGenerator(save_path, muStepA, muStepB)
uA=[];
muA=[];
for idx = 1: size(muStepA,2)
    stage_uA=zeros(muStepA(2,idx),1);
    idx_ones=randperm(muStepA(2,idx));
    idx_ones=idx_ones(1:muStepA(1,idx)*muStepA(2,idx));
    stage_uA(idx_ones) = 1;
    uA = [uA; stage_uA];
    muA=[muA;  muStepA(1,idx)*ones(muStepA(2,idx),1) ];
end
% line_color ={};
% %decision 1
% uB=uA(1:20);
% line_color = [line_color  'green'];
% muB = muA(1:20);
% % decision 0
% uB =[uB; 1 - uA(21:40)];
% muB =[muB; 1-muA(21:40)];
% line_color = [ line_color 'yellow'];
% 
% % decision 1
% uB =[uB;  uA(41:60)];
% muB =[muB; muA(41:60)];
% line_color = [line_color 'green'];
% 
% % decision 1
% uB  =[uB;   uA(61:80)];
% muB =[muB;  muA(61:80)];
% line_color = [ line_color 'green'];
% 
% % random decision 1
% uB  =[uB;  binornd(1, muStepA(1,5), [20,1])];
% muB =[muB; muStepA(1,5)*ones(20,1)];
% line_color = [ line_color  'black'];
% 
% % radom decision 0.5
% uB  =[uB;  binornd(1, 0.5, [20,1])];
% muB =[muB;  0.5*ones(20,1)];
% line_color = [ line_color  'black'];
% 
% uB= 1 - uA;
muB=[];
uB=[];
for idx = 1: size(muStepB,2)
    stage_uB=zeros(muStepB(2,idx),1);
    idx_ones=randperm(muStepB(2,idx));
    idx_ones=idx_ones(1:muStepB(1,idx)*muStepB(2,idx));
    stage_uB(idx_ones) = 1;
    uB = [uB; stage_uB];
    muB=[muB;  muStepB(1,idx)*ones(muStepB(2,idx),1) ];
end


Umab = struct;
Umab.muStepA = muStepA;
Umab.muStepB = muStepB;
Umab.uB=uB;
Umab.muB=muB;
Umab.uA=uA;
Umab.muA=muA;


fz  =14;
tfz =14;
afz =14;
outerpos = [1,1,1280, 640];
set(0,'defaultfigurecolor','w')
fig=figure(...
    'OuterPosition', outerpos, ...
    'Visible','off' );
subplot(211)
plot(1:length(Umab.muA), Umab.muA, 'black-', 'LineWidth', 2);
hold on;
%plot(1:length(Umab.muA), 0.5*ones(1,length(Umab.muA)), 'black--', 'LineWidth', 2);
%hold on;
scatter(1:length(Umab.uA), Umab.uA, 15, 'green', 'filled', 'Marker', 'o')
set(gca,  'fontsize', afz)
xlim([1 length(Umab.muA)])
ylim([-0.05, 1.05])
xlabel('Trial.',  'fontsize', fz)
ylabel('$P(x_0^{(1)}=1)$', 'interpreter','latex' ,  'fontsize', fz)
%title('Armed Bandit A','FontWeight', 'Bold', 'fontsize', tfz)
xlim([1 length(Umab.uA)])
grid on;
%box off;



subplot(212)
plot(1:length(Umab.muB), Umab.muB, 'black-', 'LineWidth', 2);
hold on;
scatter(1:length(Umab.uB), Umab.uB, 15, 'green', 'filled', 'Marker', 'o')
hold on;
%plot(1:length(Umab.muB), 0.5*ones(1,length(Umab.muB)), 'black--', 'LineWidth', 2);
%hold on;

xlim([1 length(Umab.muB)])
ylim([-0.05, 1.05])
set(gca,  'fontsize', afz)
xlabel('Trial.',  'fontsize', fz)
ylabel('$P(x_0^{(2)}=1)$',  'interpreter','latex' , 'fontsize', fz)
%title('Armed Bandit B',  'FontWeight', 'Bold', 'fontsize', tfz)
xlim([1 length(Umab.uB)])
ylim([-0.05 1.05 ])
grid on;
%box off;
print('-depsc', '-painters', '-r300', [save_path filesep 'Umab.eps'])
print('-dpng', '-painters', '-r300',  [save_path filesep 'Umab.png'])




% 
save([save_path filesep 'Umab.mat'], 'Umab')
% 
rewards= randi(4, 2, 15 * 17); %ones(2,120);
save([save_path filesep 'rewards.mat'], 'rewards')

end
