function [assignments clusterCount] = vdpgmEstimate(colorData)

cd /vega/stats/users/us2157/bb/vdpgm
vdpgmopts                               = mkopts_avdp;
vdpgmopts.get_q_of_z                    = 1;
vdpgmopts.max_target_ratio              = 0.999;
%vdpgmopts.initial_depth                 = 6;
vdpgmopts.recursive_expanding_threshold = 1.0e-2;
[~, result]                             = evalc('vdpgm(transpose(colorData), vdpgmopts)');
clusterCount                            = result.K;
[~, assignments]                        = max(result.q_of_z, [], 2);
cd /vega/stats/users/us2157/bb/bbSimulation

