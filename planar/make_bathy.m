function g = make_bathy(g)

  
% make the domain 
  xmin = 0;
  xmax = 100;
  ymin = 0;
  ymax = 200;
  zboff = -3;
  zbon = 1;
  
  
  x1d = [xmin:g.grid.dx:xmax];
  y1d = [ymin:g.grid.dx:ymax]';
  x2d = repmat(x1d,length(y1d),1);
  y2d = repmat(y1d,1,length(x1d));

  g.grid.x_walls  = x2d;
  g.grid.y_walls  = y2d;
  
  g.grid.x_c=.5*(g.grid.x_walls(2:end,1:end-1)+g.grid.x_walls(2:end,2:end));
  g.grid.y_c=.5*(g.grid.y_walls(2:end,2:end)+g.grid.y_walls(1:end-1,2:end));
  
  g.grid.zb_walls_init  = zboff+(zbon-zboff)*(x2d-xmin)./(xmax-xmin);
  g.grid.zb_c_init = interp2(g.grid.x_walls ,g.grid.y_walls,  g.grid.zb_walls_init,g.grid.x_c ,g.grid.y_c)  ;
  g.grid.activity = g.grid.zb_c_init<Inf;
 
  
  
  if g.plot.initialbathy
    figure
    surf(g.grid.x_c,g.grid.y_c,g.grid.zb_c_init)
    title('Initial Bathy','interpreter','latex','fontsize',24)
    xlabel('$x [m]$','interpreter','latex','fontsize',24)
    ylabel('$y [m]$','interpreter','latex','fontsize',24)
    zlabel('$z_b [m]$','interpreter','latex','fontsize',24)
    shading flat 
    shading interp
    if g.plot.iprint
      print('-dpng','-r300',[g.name,'/zb_init.png'])
    end

  end