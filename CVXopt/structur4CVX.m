% structure for iterations
% N = 1000;
CVXstruct.iter = 0;
CVXstruct.Rp = zeros(K, 1);
CVXstruct.Rc = zeros(K, 1);
CVXstruct.Cc = zeros(K, 1);
CVXstruct.Rov = zeros(1);
CVXstruct.P = zeros(2, 1); % (||Pc||, ||Pk||)

%%

CVXstruct_ps = CVXstruct;
CVXstruct_ps.Cc = nan;
CVXstruct_ps.P = nan;
CVXstruct_ps.s = zeros(N_R, 1);
CVXstruct_ps.sgt1 = 0; % s greater than 1;

