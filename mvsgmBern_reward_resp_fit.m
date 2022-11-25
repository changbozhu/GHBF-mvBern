function logpdf = mvsgmBern_reward_resp_fit(fit_r, prc_states, logp)

% Get states of percepture model for action response
if isfield(prc_states,'PRED_mu0')
    PRED_mu0 = prc_states.PRED_mu0;
elseif isfield(prc_states,'hatmu0')
    PRED_mu0 = prc_states.hatmu0;
end
%load('./Experiments/SimMAB/rewards4.mat', 'rewards');
rewards = fit_r.rewards;
r0 = rewards(1,:);
r1 = rewards(2,:);
[d, T] = size(PRED_mu0);
a= fit_r.a;
% Transform Uncertainty or Covariance C to its native space
zeta_a = exp(logp);
% Initialize the result log-pdf as NaNs so that NaN is
% returned for all irregualar or ignored trials

pred_mu0_1 = PRED_mu0(1,:);
pred_mu0_2 = PRED_mu0(2,:);
Qa1 = r1.*( pred_mu0_1.*pred_mu0_2 + (1-pred_mu0_1).*(1-pred_mu0_2) );
Qa0 = r0.*( pred_mu0_1.*(1-pred_mu0_2) + (1-pred_mu0_1).*pred_mu0_2 );
Pa=sgm(Qa1-Qa0);
Pa=Pa.^zeta_a ./( Pa.^zeta_a  + (1- Pa).^zeta_a );

% Weed irregular trials out from inferred states and responses
%logpdf(not(ismember(1:length(logpdf),fit_r.ignore)))
% Calculate log-probabilities for non-irregular trials
a=reshape(a, length(a),1);
logpdf= a.*log(Pa) + (1-a).*log(1-Pa);

return;
