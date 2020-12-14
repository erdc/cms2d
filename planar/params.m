
% general params
g.params.iverbose = 0;   % noisy output

% Model simulation timing params
g.time.begin_date = datenum(2016,1,1,0,0,0);
g.time.begin_date = datenum(2000,1,1,0,0,0);
g.time.duration_run = 3; %hours
g.time.dt_bc = 1/24; % hours

% grid stuff
g.grid.dx=1;
g.grid.dx=2;
%g.grid.dx = .5
g.grid.dx=2;


%current stuff
g.fw = 0.05; 
g.n = 0.03; 


% sed stuff
g.transport_formula = 'CSHORE';
g.transport_formula = 'LUND-CIRP';
g.cseffb = .00;
%g.cseffb = .005
g.csblp  = .00;
%g.csblp  = .00
g.csslp  = .5;