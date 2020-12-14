if g.plot.morpho
  extratext = 'yeswavesedtrans';
  extratext = '';
  fs = 24;
  x = out.xy.x;
  y = out.xy.y;
  dzbtot = sum(out.dzb.data2d,3);

  % bed conservation check
  x0 = 20;
  x1 = 40;
  y0 = 100;
  y1 = 140;
  dx = .1;
  
  
  figure
  quiver(x,y,out.Qt.data2d.x(:,:,end),out.Qt.data2d.y(:,:,end),'linewidth',3);hold on
  title(['Total Sediment Transport ',extratext])
  xlabel('$x [m]$','interpreter','latex','fontsize',fs)
  ylabel('$y [m]$','interpreter','latex','fontsize',fs)
  plot([x0 x1 x1 x0 x0],[y0 y0 y1 y1 y0],'r')
  if g.plot.iprint
    print('-dpng','-r300',[g.name,'/transport_vec_',extratext,'.png'])
  end
  
  xi = x0:dx:x1;  
  yi = [y0:dx:y1]';
  xi2d = repmat(xi,length(yi),1);
  yi2d = repmat(yi,1,length(yi));
  dt = out.Qt.time(2)-out.Qt.time(1);
  for i = 1:length(out.dzb.time)
    Qx0(i) = dx*sum(interp2(x,y,out.Qt.data2d.x(:,:,i),x0,yi)/2650);
    Qx1(i) = dx*sum(interp2(x,y,out.Qt.data2d.x(:,:,i),x1,yi)/2650);
    Qy0(i) = dx*sum(interp2(x,y,out.Qt.data2d.y(:,:,i),xi,y0)/2650);
    Qy1(i) = dx*sum(interp2(x,y,out.Qt.data2d.y(:,:,i),xi,y1)/2650);
    zbi = interp2(x,y,-out.dep.data2d(:,:,i),xi,yi);
    Vol_nopores(i) = .6*dx^2*sum(sum(zbi));
  end
  
  figure
  surf(x,y,dzbtot)
  title('Profile Change')
  xlabel('$x [m]$','interpreter','latex','fontsize',fs)
  ylabel('$y [m]$','interpreter','latex','fontsize',fs)
  if g.plot.iprint
    print('-dpng','-r300',[g.name,'/bed_change.png'])
  end

  figure
  clear hh;
  hh(1) = plot(out.dep.time,(Vol_nopores-Vol_nopores(1)),'-rs','markerfacecolor','k');hold all
  hh(2) = plot(out.dep.time,dt*3600*cumsum(Qx0),'-bs','markerfacecolor','k');
  hh(3) = plot(out.dep.time,-dt*3600*cumsum(Qx1),'--bs','markerfacecolor','k');
  hh(4) = plot(out.dep.time,dt*3600*cumsum(Qy0),'-ks','markerfacecolor','k');
  hh(5) = plot(out.dep.time,-dt*3600*cumsum(Qy1),'--ks','markerfacecolor','k');
  hh(6) = plot(out.dep.time,dt*3600*cumsum(Qx0-Qx1+Qy0-Qy1),'--r','linewidth',3);

  title(['Control Volume ',extratext])
  xlabel('$t [hr]$','interpreter','latex','fontsize',fs)
  ylabel('$Vol [m^3]$','interpreter','latex','fontsize',fs)
  grid on
  legend(hh,'Sand in Bed','V from Qx0','V from -Qxl','V from Qy0','V from -Qyl','V from all bnd','location','southwest')
  if g.plot.iprint
    print('-dpng','-r300',[g.name,'/control_vol_',extratext,'.png'])
  end
  
  
  
  
  
  

end