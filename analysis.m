close all;
clear all;
clc;
afz=12;
fz=12;
tfz=12;

data_root='./sim2';
IDX=1:1:100;
ghbf_R = NaN(length(IDX),3);
BF=NaN(length(IDX), 1); % Bayesian Factor (BF)
BF_BIC=NaN(length(IDX), 1); % Bayesian Factor (BF)
pred_state_error=NaN(length(IDX), 2);

T=255;
ghbf_mu = zeros(2,T,100);
RWM_mu=zeros(2,T,100);
d_xi_ghbf=14;
d_xi_rw=2;
for i=1:length(IDX)
    idx=IDX(i);
    miniGHBF_mvBernBrownian_fit_r = load([data_root filesep 'env' num2str(idx) filesep 'miniGHBF_mvBernBrownian_fit_r.mat'], 'fit_r');
    miniGHBF_RW_mvMAB_HP_fit_r = load([data_root filesep 'env' num2str(idx) filesep 'miniGHBF_RW_mvMAB-HP_fit_r.mat'], 'fit_r');
    
    ghbf_mu(1,:,idx)=sgm(miniGHBF_mvBernBrownian_fit_r.fit_r.prc_states.mu1(1,:), 1);
    ghbf_mu(2,:,idx)=sgm(miniGHBF_mvBernBrownian_fit_r.fit_r.prc_states.mu1(2,:), 1);
    RWM_HP_mu(:,:,idx)=miniGHBF_RW_mvMAB_HP_fit_r.fit_r.prc_states.mu0;

    % Store response parameters
    % Please reference to <resp_model>_namespace.m for more details
    sim_r = struct('rewards', miniGHBF_mvBernBrownian_fit_r.fit_r.rewards);
    resp_simconfig = str2func([miniGHBF_mvBernBrownian_fit_r.fit_r.resp_conf.resp_model, '_simconfig']);
    resp_para = resp_simconfig(sim_r);
    
    % Get the handle to response model
    resp_fun = str2func([miniGHBF_mvBernBrownian_fit_r.fit_r.resp_conf.resp_model '_sim']);
        
    % Simulation for decision-making
    prc_states=struct('hatmu0', miniGHBF_mvBernBrownian_fit_r.fit_r.mu);
    optimal_mab_resp_r = resp_fun(sim_r, prc_states, resp_para);
    model_R=ghbf_R(:, [2,3]);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%
    %%%     optimal agent
    %%%
    if idx == 65
        v='on';
    else
        v='off';
    end
    outerpos = [100, 100, 640, 500];
    set(0,'defaultfigurecolor','w')
    fig=figure('OuterPosition', outerpos,  'Visible', v,  'Name', 'Regret' );
    set(fig,'Name','Regret')
    subplot(2,1,1)
    plot(1:T,optimal_mab_resp_r.prob_a1, 'b--', 'LineWidth', 1)
    hold on;
    plot(1:T,miniGHBF_mvBernBrownian_fit_r.fit_r.resp_r.prob_a1, 'r-', 'LineWidth', 1)
    set(gca, 'fontname', 'Arial', 'fontsize', afz )
    xlim([1,T])
    ylim([0, 1])
    ylabel('$P(a=1)$', 'Interpreter','latex','fontname', 'Arial', 'FontSize', afz)
    xlabel('Trial Number', 'Interpreter','latex','fontname', 'Arial', 'FontSize', afz)
    box on;
    grid on;
    hold on;

    subplot(2,1,2)
    opt_a = 1-abs(miniGHBF_mvBernBrownian_fit_r.fit_r.u(3,:) - miniGHBF_mvBernBrownian_fit_r.fit_r.u(2,:));
    ghbf_a=miniGHBF_mvBernBrownian_fit_r.fit_r.resp_r.a;
    random_agent_a=optimal_mab_resp_r.a;
    rp=[];
    rs=[];
    opt_rp =[];
    opt_rs = [];
    for k=1:T
        opt_rp = [opt_rp, miniGHBF_mvBernBrownian_fit_r.fit_r.resp_r.rewards(opt_a(k)+1, k)];
        if ghbf_a(k)==opt_a(k)
            rp=[rp, (1-ghbf_a(k))*miniGHBF_mvBernBrownian_fit_r.fit_r.resp_r.rewards(1,k) + ghbf_a(k)*miniGHBF_mvBernBrownian_fit_r.fit_r.resp_r.rewards(2,k)];
        else
            rp=[rp, 0];
        end
        rs=[rs, sum(rp(1:k))];
        opt_rs = [opt_rs, sum(opt_rp(1:k))];
    end
    plot(1:T, rs, 'r-', 'LineWidth', 1)
    hold on;
    plot(1:T, opt_rs, 'g--', 'LineWidth', 1)
    hold all
    rp=[];
    rs=[];
    opt_rp =[];
    opt_rs = [];
    for k=1:T
        opt_rp = [opt_rp, optimal_mab_resp_r.rewards(opt_a(k)+1, k)];
        if random_agent_a(k)==opt_a(k)
            rp=[rp, (1-random_agent_a(k))*optimal_mab_resp_r.rewards(1,k) + random_agent_a(k)*optimal_mab_resp_r.rewards(2,k)];
        else
            rp=[rp, 0];
        end
        rs=[rs, sum(rp(1:k))];
        opt_rs = [opt_rs, sum(opt_rp(1:k))];
    end
    plot(1:T, rs, 'b--', 'LineWidth', 1)
    set(gca, 'fontsize', afz )
    %title({'\bf{Score (red) and ideal score (green)}'} , 'interpreter', 'latex', 'fontsize', fz);
    ylabel('Score', 'interpreter', 'latex', 'fontsize', afz );
    xlabel('Trial number', 'fontsize', afz);
    ylim([0, max(opt_rs)+1])
    xlim([0 T+1])



    print('-depsc', '-painters', '-r300', [data_root filesep 'env' num2str(idx) filesep 'regret.eps'])
    print('-dpng', '-painters', '-r300',  [data_root filesep 'env' num2str(idx) filesep 'regret.png'])

    ghbf_R(idx, 1) = sum(abs(optimal_mab_resp_r.a - optimal_mab_resp_r.prob_a1));
    ghbf_R(idx, 2) = sum(abs(optimal_mab_resp_r.prob_a1 - miniGHBF_mvBernBrownian_fit_r.fit_r.resp_r.prob_a1));
    ghbf_R(idx, 3) = sum(abs(optimal_mab_resp_r.prob_a1 - miniGHBF_RW_mvMAB_HP_fit_r.fit_r.resp_r.prob_a1));
    
    BF(idx) = exp(miniGHBF_mvBernBrownian_fit_r.fit_r.F - miniGHBF_RW_mvMAB_HP_fit_r.fit_r.F);
    BF_BIC(idx)  = exp(miniGHBF_mvBernBrownian_fit_r.fit_r.F - miniGHBF_RW_mvMAB_HP_fit_r.fit_r.F - (d_xi_ghbf - d_xi_rw)/2*log(T) );
    pred_state_error(idx,1) = sum(sqrt(sum( (miniGHBF_mvBernBrownian_fit_r.fit_r.prc_states.PRED_mu0 - miniGHBF_mvBernBrownian_fit_r.fit_r.mu).^2, 1 )));
    pred_state_error(idx,2) = sum(sqrt(sum( (miniGHBF_RW_mvMAB_HP_fit_r.fit_r.prc_states.hatmu0 - miniGHBF_RW_mvMAB_HP_fit_r.fit_r.mu).^2, 1 )));
end
mu = miniGHBF_mvBernBrownian_fit_r.fit_r.mu;
ghbf_mean_Tvariant_mu=mean(ghbf_mu, 3);
ghbf_std_Tvariant_mu=std(ghbf_mu,0, 3);

RWM_HP_mean_Tvariant_mu = mean(RWM_HP_mu,3);
RWM_HP_std_Tvariant_mu = std(RWM_HP_mu, 0, 3);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%     mean states ghbf vs. rw-HP
%%%
outerpos = [100,100,740, 680];
set(0,'defaultfigurecolor','w')
fig=figure('OuterPosition', outerpos,  'Visible','on');
set(fig,'Name','mean states (ghbf vs. rw)');
subplot(2, 1, 1);
plot(1:T, ghbf_mean_Tvariant_mu(1,:), '-r', 'LineWidth', 1);
hold all;
fill([(1:T)'; flipud((1:T)')], [ghbf_mean_Tvariant_mu(1,:)'+ghbf_std_Tvariant_mu(1,:)';  flipud(ghbf_mean_Tvariant_mu(1,:)'-ghbf_std_Tvariant_mu(1,:)')], ...
         'r-', 'EdgeAlpha', 0, 'FaceAlpha', 0.25);   
hold all;
plot(1:T, RWM_HP_mean_Tvariant_mu(1,:), '--b', 'LineWidth', 1);
hold all;
fill([(1:T)'; flipud((1:T)')], [RWM_HP_mean_Tvariant_mu(1,:)'+RWM_HP_std_Tvariant_mu(1,:)';  flipud(RWM_HP_mean_Tvariant_mu(1,:)'-RWM_HP_std_Tvariant_mu(1,:)')], ...
         'b-', 'EdgeAlpha', 0, 'FaceAlpha', 0.4);   
hold all;
plot(1:T, mu(1,:), '-black', 'LineWidth', 1);
hold all;
set(gca, 'fontsize', afz )
% title({'\bf{States of armed bandit A (green)}'; '\bf{posterior expectation of armed bandit A (red)}'}, ...
%      'interpreter', 'latex', 'fontsize', fz);

%plot(1:T, 0.5*ones(1,T), 'k--', 'LineWidth', 0.5);
ylabel('$\mu_0^{(1)}$', 'interpreter', 'latex', 'FontName', 'Arial', 'fontsize', afz)
xlabel('Trial number', 'FontName', 'Arial', 'fontsize', afz);
ylim([0, 1])
xlim([1 T])
grid on;
hold off;

subplot(2,1,2);
plot(1:T, ghbf_mean_Tvariant_mu(2,:), '-r', 'LineWidth', 1);
hold all;
fill([(1:T)'; flipud((1:T)')], [ghbf_mean_Tvariant_mu(2,:)'+ghbf_std_Tvariant_mu(2,:)';  flipud(ghbf_mean_Tvariant_mu(2,:)'-ghbf_std_Tvariant_mu(2,:)')], ...
         'r-', 'EdgeAlpha', 0, 'FaceAlpha', 0.25);
hold all;
plot(1:T, RWM_HP_mean_Tvariant_mu(2,:), '--b', 'LineWidth', 1);
hold all;
fill([(1:T)'; flipud((1:T)')], [RWM_HP_mean_Tvariant_mu(2,:)'+RWM_HP_std_Tvariant_mu(2,:)';  flipud(RWM_HP_mean_Tvariant_mu(2,:)'-RWM_HP_std_Tvariant_mu(2,:)')], ...
         'b-', 'EdgeAlpha', 0, 'FaceAlpha', 0.4);  
hold all;
plot(1:T, mu(2,:), '-black', 'LineWidth', 1);
hold all;
set(gca, 'fontsize', afz )
% title({'\bf{States of armed bandit B (green)}'; '\bf{posterior expectation of armed bandit B (red)}'}, ...
%      'interpreter', 'latex', 'fontsize', fz);
%plot(1:T, 0.5*ones(1,T), 'k--', 'LineWidth', 0.5);
ylabel('$ \mu_0^{(2)} $', 'interpreter', 'latex', 'FontName', 'Arial', 'fontsize', afz)
xlabel('Trial number', 'FontName', 'Arial', 'fontsize', afz);
ylim([0, 1])
xlim([1 T])
grid on;
hold off;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%     Regret
%%%
model_R=ghbf_R(:, [2,3]);
outerpos = [100,100,740, 740];
set(0,'defaultfigurecolor','w')
fig=figure('OuterPosition', outerpos,  'Visible','on',  'Name', 'Regret' );
set(fig,'Name','Regret')
X = categorical({'GHBF', 'RW'});
X = reordercats(X,{'GHBF','RW'});
b=bar(X, mean(model_R, 1),'FaceColor','black');
hold all;
errorbar(mean(model_R, 1), std(model_R,1), 'blue', 'Linestyle', 'None', 'LineWidth', 2)
hold on;
xtips1 = b(1).XEndPoints;
ytips1 = b(1).YEndPoints+std(model_R,1);
labels1 = ["$"  "$"] + string(b(1).YData) + ["\pm"  "\pm"] + string( std(model_R,1) ) + ["$" "$"]  ;
text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom', 'Interpreter', 'latex' ,'FontSize', fz)

set(gca, 'fontname', 'Arial', 'fontsize', afz )
ylim([0, max(mean(model_R, 1)+std(model_R,1)) *1.1])
ylabel('Regret', 'Interpreter','latex','fontname', 'Arial', 'FontSize', afz)
box on;
grid off;
hold off;
print('-depsc', '-painters', '-r300', [data_root filesep 'regret.eps'])
print('-dpng', '-painters', '-r300',  [data_root filesep 'regret.png'])

outerpos = [100,100,860, 500];
set(0,'defaultfigurecolor','w')
fig=figure('OuterPosition', outerpos,  'Visible','on',  'Name', 'Bayesian Factor1' );
set(fig,'Name','ghbf_vs_rw_HP-Bayesian Factor');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%                 Bayesian Factor
%%%

%edges=[10-7, 0.01, 0.1, 0.333, 1, 3, 10, 100  max(BF)];
subplot(1,2,1)
numBF0=sum(BF<0.01);
numBF1=sum(BF<0.1)-sum(BF<0.01);
numBF2=sum(BF<0.333)-sum(BF<0.1);
numBF3=sum(BF<1)-sum(BF<0.333);
numBF4=sum(BF<3)-sum(BF<1);
numBF5=sum(BF<10)-sum(BF<3);
numBF6=sum(BF<100)-sum(BF<10);
numBF7=sum(BF>100);
probBFcat=[numBF0, numBF1, numBF2, numBF3, numBF4, numBF5, numBF6, numBF7]./length(IDX);

X = categorical({'BF<0.01', '0.01<BF<0.1', '0.1<BF<0.333', '0.333<BF<1', '1<BF<3', '3<BF<10', '10<BF<100', 'BF>100'});
X = reordercats(X,{'BF<0.01', '0.01<BF<0.1', '0.1<BF<0.333', '0.333<BF<1', '1<BF<3', '3<BF<10', '10<BF<100', 'BF>100'});
b=bar(X, probBFcat,'FaceColor','black');


set(gca, 'fontsize', afz )
ylim([0, 1])
ylabel('Prob.', 'Interpreter','latex', 'FontSize', afz)
box on;
grid on;
hold off;
% print('-depsc', '-painters', '-r300', [data_root filesep 'BF_GHBF_RW_HP.eps'])
% print('-dpng', '-painters', '-r300',  [data_root filesep 'BF_GHBF_RW_HP.png'])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%                 Bayesian Factor_BIC
%%%
% outerpos = [100,100,740, 740];
% set(0,'defaultfigurecolor','w')
% fig=figure('OuterPosition', outerpos,  'Visible','on',  'Name', 'Bayesian Factor with BIC' );
% set(fig,'Name','ghbf_vs_rw_HP-Bayesian Factor');
%edges=[10-7, 0.01, 0.1, 0.333, 1, 3, 10, 100  max(BF)];
subplot(1,2,2)
numBF0=sum(BF_BIC<0.01);
numBF1=sum(BF_BIC<0.1)-sum(BF_BIC<0.01);
numBF2=sum(BF_BIC<0.333)-sum(BF_BIC<0.1);
numBF3=sum(BF_BIC<1)-sum(BF_BIC<0.333);
numBF4=sum(BF_BIC<3)-sum(BF_BIC<1);
numBF5=sum(BF_BIC<10)-sum(BF_BIC<3);
numBF6=sum(BF_BIC<100)-sum(BF_BIC<10);
numBF7=sum(BF_BIC>100);
probBFcat=[numBF0, numBF1, numBF2, numBF3, numBF4, numBF5, numBF6, numBF7]./length(IDX);

X = categorical({'BF_{BIC}<0.01', '0.01<BF_{BIC}<0.1', '0.1<BF_{BIC}<0.333', '0.333<BF_{BIC}<1', '1<BF_{BIC}<3', '3<BF_{BIC}<10', '10<BF_{BIC}<100', 'BF_{BIC}>100'});
X = reordercats(X,{'BF_{BIC}<0.01', '0.01<BF_{BIC}<0.1', '0.1<BF_{BIC}<0.333', '0.333<BF_{BIC}<1', '1<BF_{BIC}<3', '3<BF_{BIC}<10', '10<BF_{BIC}<100', 'BF_{BIC}>100'});
b=bar(X, probBFcat,'FaceColor','black');


set(gca, 'fontsize', afz )
ylim([0, 1])
ylabel('Prob.', 'Interpreter','latex', 'FontSize', afz)
box on;
grid on;
hold off;
print('-depsc', '-painters', '-r300', [data_root filesep 'BF_BIC_GHBF_RW_HP.eps'])
print('-dpng', '-painters', '-r300',  [data_root filesep 'BF_BIC_GHBF_RW_HP.png'])