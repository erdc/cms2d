if g.plot.sed
  
  fs = 24;
  figpos = [600 500 1000 800];
  x = out.xy.x;
  y = out.xy.y;
  ind_center = round(size(x,1)/2);
  
  hh=figure;
  set(hh, 'Position', figpos)
  h         = out.dep.data2d(ind_center,:,end);
  ctstar = out.Ctstar.data2d(ind_center,:,end);
  Vs = h.*ctstar/(1000*2.65);
  x_center = x(ind_center,:);
  plot(x_center,Vs,'linewidth',2);hold all
  %  plot(eta.x,sqrt(2)*eta.Hrms,'rs','markerfacecolor','k')
  
  title('Integrated Sediment Volume')
  xlabel('$x [m]$','interpreter','latex','fontsize',fs)
  ylabel('$V_s [m]$','interpreter','latex','fontsize',fs)
  if g.plot.iprint
    set(gcf,'PaperPositionMode','auto')
    print('-dpng','-r300',[g.name,'/Vs.png'])
  end
  
  
  figure
  quiver(x,y,out.Qt.data2d.x(:,:,end),out.Qt.data2d.y(:,:,end),'linewidth',3)
  title('Total Sediment Transport')
  xlabel('$x [m]$','interpreter','latex','fontsize',fs)
  ylabel('$y [m]$','interpreter','latex','fontsize',fs)
  if g.plot.iprint
    set(gcf,'PaperPositionMode','auto')
    print('-dpng','-r300',[g.name,'/Q.png'])
  end
  
  if exist('fort.1234')==2
    load fort.1234
    numtimes = length(find(fort(:,1)==1));
    N = size(fort,1)/numtimes;
    dum = fort(end-N+1:end,:);
    qbx = reshape(dum(:,2),size(x,2),size(x,1))';
    qsx = reshape(dum(:,3),size(x,2),size(x,1))';
    qtx = reshape(dum(:,4),size(x,2),size(x,1))';
    load fort.1235
    numtimes = length(find(fort(:,1)==1));
    N = size(fort,1)/numtimes;
    dum = fort(end-N+1:end,:);
    qby = reshape(dum(:,2),size(x,2),size(x,1))';
    qsy = reshape(dum(:,3),size(x,2),size(x,1))';
    qty = reshape(dum(:,4),size(x,2),size(x,1))';

    hh=figure;
    %set(hh, 'Position', figpos)
    set(hh, 'Position', [680   558   800   300])
    dec = 3;
    zb = g.grid.zb_c_init(ind_center,:);
    dzbdx=[(zb(2:end)-zb(1:end-1))/g.grid.dx' 0];
    dzbdx2=dzbdx(1:dec:end);
    x2 = x_center(1:dec:end);
    zb2 = zb(1:dec:end);
    qbx2=qbx(ind_center,1:dec:end);
    qsx2=qsx(ind_center,1:dec:end);
    fill([x_center x_center(end) x_center(1)],[zb -.9 -.9],[.8 .8 .5]);hold on
    amp = 300;
    quiver(x2,zb2+.05,amp*qbx2,amp*qbx2.*dzbdx2,0,'r','linewidth',2,'maxheadsize',.01)
    quiver(x2,.5*zb2,amp*qsx2,amp*qbx2*0,0,'k','linewidth',2,'maxheadsize',.01)
    plot([0 21.5],[0 0],'b')
    plot([x2;x2],[0*zb2;zb2],'k')
    set(gca,'ylim',[-.9 .2])
    set(gca,'xlim',[0 24])
    title('Bedload and Suspended load')
    xlabel('$x [m]$','interpreter','latex','fontsize',fs)
    ylabel('$z [m]$','interpreter','latex','fontsize',fs)
    if g.plot.iprint
      %set(gcf, 'Position', [680   558   800   300])
    set(gcf,'PaperPositionMode','auto')
    print('-dpng','-r300',[g.name,'/transport_2d.png'])
    end

    hh=figure;
    set(hh, 'Position', figpos)
  
    plot(x_center,qbx(ind_center,:)/2650,'r-','linewidth',2);hold on
    plot(x_center,qsx(ind_center,:)/2650,'k-','linewidth',2)
    plot(x_center,qtx(ind_center,:)/2650,'k--','linewidth',2)
    plot(x_center,out.Qt.data2d.x(ind_center,:,end)/2650,'k-x','linewidth',2)
    plot(x,0*x,'b','linewidth',2)
    text(1,.004/2650,'Bedload','fontsize',16,'color','r')
    text(1,-.005/2650,'Suspended load','fontsize',16,'color','k')
    title('Bedload and Suspended load')
    xlabel('$x [m]$','interpreter','latex','fontsize',fs)
    ylabel('$Q [\frac{m^2}{s}]$','interpreter','latex','fontsize',fs)
    if g.plot.iprint
    set(gcf,'PaperPositionMode','auto')
    print('-dpng','-r300',[g.name,'/transport.png'])
    end
    
      
  end
end