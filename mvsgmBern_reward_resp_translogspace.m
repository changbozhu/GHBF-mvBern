function [pmus, resp_para] = mvsgmBern_reward_resp_translogspace(fit_r, logp)
pnas = fit_r.resp_conf.pnas;
pdims = fit_r.resp_conf.pdims;

logpmus  = logp(1:sum(pdims));
pmus     = NaN(length(logpmus),1);

resp_para = struct;
%1.zeta_a
i =1;
j =sum(pdims(1:1));
pmus(i:j)  =exp(logpmus( i:j ));
zeta_a = pmus( i:j );
resp_para.zeta_a = zeta_a;

return;
