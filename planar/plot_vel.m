if g.plot.vel
  fs = 24;
  hh=figure;
  quiver(g.grid.x_c,g.grid.y_c,out.vel.data2d.x(:,:,end),out.vel.data2d.y(:,:,end))
  xlabel('$x [m]$','interpreter','latex','fontsize',fs)
  ylabel('$y [m]$','interpreter','latex','fontsize',fs)
  axis equal
  
end