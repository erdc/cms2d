if g.plot.eta
  fs = 30;
  
  for i = 1:length(out.eta.time)
    figure
    eta = out.eta.data2d(:,:,i);
    eta(eta<-900) = NaN;
    surf(g.grid.x_c,g.grid.y_c,eta )
    title(['Free surface position, time = ',num2str(out.eta.time(i)),' ',num2str(mean(eta(:,1)))])
    xlabel('$x [m]$','interpreter','latex','fontsize',fs)
    ylabel('$y [m]$','interpreter','latex','fontsize',fs)
    zlabel('$z [m]$','interpreter','latex','fontsize',fs)
    shading flat 
    shading interp

    if g.plot.iprint
      print('-dpng','-r300',[g.name,'/eta_2d.png'])
    end
  end
  
end