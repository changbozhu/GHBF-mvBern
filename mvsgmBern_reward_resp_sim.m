function resp_r = mvsgmBern_reward_resp_sim(sim_r, prc_states, resp_para)
% Get states of percepture model for action response
if isfield(prc_states,'PRED_mu0')
    PRED_mu0 = prc_states.PRED_mu0;
    pred_mu0_1 = PRED_mu0(1,:);
    pred_mu0_2 = PRED_mu0(2,:);
elseif isfield(prc_states,'hatmu0')
    hatmu0 = prc_states.hatmu0;
    pred_mu0_1 = hatmu0(1,:);
    pred_mu0_2 = hatmu0(2,:);
end
%load('./Experiments/SimMAB/rewards.mat', 'rewards');
rewards=sim_r.rewards;
r0 = rewards(1,:);
r1 = rewards(2,:);
zeta_a = resp_para.zeta_a;
Qa1 = r1.*( pred_mu0_1.*pred_mu0_2 + (1-pred_mu0_1).*(1-pred_mu0_2) );
Qa0 = r0.*( pred_mu0_1.*(1-pred_mu0_2) + (1-pred_mu0_1).*pred_mu0_2 );
Pa=sgm(Qa1-Qa0);
Pa=Pa.^zeta_a ./( Pa.^zeta_a  + (1- Pa).^zeta_a );
a=binornd(1, Pa);
resp_r = struct('prob_a1',Pa, 'a', a, 'rewards', rewards);
return;
