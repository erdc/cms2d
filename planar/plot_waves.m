if g.plot.waves
  fs = 30;
  ind_center = round(size(g.grid.x_c,1)/2);
  figure
  surf(g.grid.x_c,g.grid.y_c,out.Hs.data2d(:,:,end))
  title('Wave Height','interpreter','latex','fontsize',fs)
  xlabel('$x$','interpreter','latex','fontsize',fs)
  ylabel('$y$','interpreter','latex','fontsize',fs)
  zlabel('$H_s$','interpreter','latex','fontsize',fs)
  shading flat 
  shading interp
  if g.plot.iprint
    print('-dpng','-r300',[g.name,'/Hs_2d.png'])
  end

end