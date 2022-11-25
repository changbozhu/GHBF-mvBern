function resp_conf = mvsgmBern_reward_resp_fitconfig(da)

% Create config structure and Store name of percepture model
resp_conf = struct;
resp_conf.resp_model = 'mvsgmBern_reward_resp';

% Store dimensions of each level
resp_conf.xdims = da;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Configure prior of all structural parameters.
%
dze=1;
logzeta_mu0=log(2);
logzeta_aC0=1e-10;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Gather prior settings in four vectors and one cell
resp_conf.pdims = [...
                     dze...
                      ];
resp_conf.pnas = {...
                    'zeta_a'...
                      };
resp_conf.opt_pbos = [...
                     0....
                        ];
resp_conf.priormus  = [...
                     logzeta_mu0...
                      ];
resp_conf.priorCs  =  [...
                     logzeta_aC0...
                      ];
num_paras = sum(resp_conf.pdims)+sum(resp_conf.pdims.^2);
resp_conf.num_paras = num_paras;

%--------------------------------------------------------------------------
% Model function handle
resp_conf.model = str2func([resp_conf.resp_model  '_fit']);

% Handle to function that transforms perceptual parameters to their native space
% from the space they are estimated in
resp_conf.translogspace = str2func([resp_conf.resp_model '_translogspace']);
disp('RESPONSE-MODEL-INFORMATION:  ')
disp(['--resp_model:  ' resp_conf.resp_model ])
disp(['--The number of optimized parameters: ' num2str(sum(resp_conf.opt_pbos))])
disp(['--Total Dimensions of all parameters: ' num2str(sum(resp_conf.pdims))])
disp(['--Total Dimensions of optimized parameters: ' num2str(sum(resp_conf.opt_pbos'*resp_conf.pdims))])
return;
