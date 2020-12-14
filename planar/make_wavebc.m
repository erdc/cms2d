function g = make_wavebc(g)

g.bc.Hrms = 1;
g.bc.T = 8;
g.bc.Dir_STWAVE_convention = 10;
g.bc.Dir_STWAVE_convention = 0;

in.name = g.name;
in.num_f_bins = 30; % number of freqency bins
in.f_min=1/(3*g.bc.T); % lowest freq 1/s
in.f_max=1/(g.bc.T/3); % highest freq 1/s
in.s = 10^10; % directional spread param, larger s = smaller directional spread
in.gamma = 3.3; % peak param
                %in.gamma = 10^10; % peak param
in.dF = (in.f_max-in.f_min)/(2*in.num_f_bins);
in.F = in.f_min:in.dF:in.f_max;
in.h = -mean(g.grid.zb_c_init(:,1)); %depth
in.Hrms = g.bc.Hrms;
in.T = g.bc.T;
in.alpha_p = g.bc.Dir_STWAVE_convention;
in.fp=1./in.T;
in = tma(in);
in.na = 35; % number of angle bins

g.bc.in=in;


%plotbc