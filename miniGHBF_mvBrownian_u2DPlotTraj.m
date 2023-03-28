function miniGHBF_mvBrownian_u2DPlotTraj(sim_r)
% Plot dynamic trajectories of GHBF which accepts a series
% of 2D continous observations u. 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % AUTHER:  Changbo Zhu
 % E-mail:  changbozhu@outlook.com
 % DATE:    Nov. 18th 2020
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %
 prc_states=sim_r.prc_states;
 %
 % Get total times 
 u = sim_r.u;
 % Adding dummy "zeroth" trial
 u = [ zeros(size(u, 1),1), u ];
 % Number of trials (including prior)
 T = size(u,2)-1
 PriorC1 = prc_states.PRED_C1;
 priorcorr1 = squeeze(PriorC1(2,1,:)./ sqrt(PriorC1(1,1,:).*PriorC1(2,2,:)));
 postcorr1 = squeeze(prc_states.C1(2,1,:) ./ sqrt(prc_states.C1(1,1,:).*prc_states.C1(2,2,:)));
 eigval = [0,0];
 for k = 2:T
    C1 = squeeze(prc_states.C1(:,:,k));
    eigval = [eigval; eig( C1 )'] ;
 end
 eigval_max_min = min([max(eigval(:,2)), max(eigval(:,1))])*5;
 % sigmxy = [sqrt(max(C1(1,1,:))), sqrt(max(C1(2,2,:)))]
  
 %configure parameters in the 1st level
 NumSamples = 400;
 X1Min=min(min(prc_states.mu1))  - 0.1*abs(max(max(prc_states.mu1)) - min(min(prc_states.mu1)));
 X1Max=max(max(prc_states.mu1))  + 0.1*abs(max(max(prc_states.mu1)) - min(min(prc_states.mu1)));
 uxMin = min(min(prc_states.u_x(:,2:end) ))  - 0.1*abs(max(max(prc_states.u_x(:,2:end) )) - min(min(prc_states.u_x(:,2:end) ))); 
 uxMax=max(max(prc_states.u_x(:,2:end)))  + 0.1*abs(max(max(prc_states.u_x(:,2:end))) - min(min(prc_states.u_x(:,2:end))));
 X1rang = [min([X1Min  uxMin])  max([X1Max uxMax])];
 [X, Y] = meshgrid(linspace(X1rang(1), X1rang(2), NumSamples), linspace(X1rang(1), X1rang(2), NumSamples));
 XY = [X(:) Y(:)];
 
 PE1Min=min(min(prc_states.PE1))  - 0.1*abs(max(max(prc_states.PE1)) - min(min(prc_states.PE1)));
 PE1Max=max(max(prc_states.PE1))  + 0.1*abs(max(max(prc_states.PE1)) - min(min(prc_states.PE1)));
 %
 uMin =min(min(u(2:3,2:end)))  - 0.1*abs(max(max(u(2:3,2:end))) - min(min(u(2:3,2:end))));
 uMax=max(max(u(2:3,2:end)))  + 0.1*abs(max(max(u(2:3,2:end))) - min(min(u(2:3,2:end))));
 %
 PEuMin=min(min(prc_states.PEu))  - 0.1*abs(max(max(prc_states.PEu)) - min(min(prc_states.PEu)));
 PEuMax=max(max(prc_states.PEu))  + 0.1*abs(max(max(prc_states.PEu)) - min(min(prc_states.PEu)));
 %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 % plot states in English
 % Set up display 
 scrsz = get(0,'screenSize');
 imgsz = min([scrsz(3), scrsz(4)]);
 %outerpos = [320,60,1600,1020];
 outerpos = [1,1,imgsz*1.151,imgsz];
 set(0,'defaultfigurecolor','w')
 figure(...
    'OuterPosition', outerpos,...
    'Name', 'General Hierarchical Brownian Filter (GHBF)');
 for k =1:1:T
    clf;
    subplot(4,4,[1,2,5,6])
    titleText = 'Phase Map: Expectation $\mu_1$ of r.v. $x_1$';
    % plot bump in the 1st level
    pz = mvnpdf( XY, prc_states.mu1(:,k)', prc_states.C1(:,:,k)./eigval_max_min);
    %colorbar('YTickLabel',{'0',      '1'});
    colormap hot;
    %contourf(X, Y, reshape(pz,size(X)),0,':','LineStyle','none')
    imagesc(linspace(X1rang(1), X1rang(2), NumSamples),...
        linspace(X1rang(1), X1rang(2), NumSamples),...
        reshape(pz,size(X)),[min(pz(:)) max(pz(:))+1e-3])
    axis equal;
    axis xy ;
    colorbar('location','eastoutside');  
    hold on;
    if k > 1
        plot(prc_states.u_x(1,2:k),prc_states.u_x(2,2:k),'b+','LineWidth',2,'markersize',3)
        hold on;
    end
    axis equal;
    hold on;
    plot(prc_states.mu1(1,1:k),prc_states.mu1(2,1:k),'g-','LineWidth',1,'markersize',3);
    hold on;
    xlim(X1rang);
    ylim(X1rang);
    xlabel('$x_1^{(1)}$','Interpreter','latex');
    ylabel('$x_1^{(2)}$','Interpreter','latex');
    axis equal
    hti=title(titleText);
    set(hti,'Interpreter','latex')
    %axis off;
    hold off;
    %
    % Plot Sensory input u and prediction error PEu  
    subplot(443)
    hold on;
    if k >2
        %plot(prc_states.PEu(1,2:k),prc_states.PEu(2,2:k),'y*','LineWidth',1,'markersize',2);
        %hold on;
        plot(u(2,2:k),u(3,2:k),'b+','LineWidth',1,'markersize',3);
        hold on;
        plot(prc_states.PRED_sens(1,2:k),prc_states.PRED_sens(2,2:k),'r*','LineWidth',1,'markersize',3);
        hold on;
    end
    xlim([uMin uMax]);
    ylim([uMin uMax]);
    hti=title({'Sensory Input $u$ (blue points )';...
        ['Prediction ' '$\hat{u}$' '(red points ) ']});
    set(hti,'Interpreter','latex')
    xlabel('$u^{(1)}$','Interpreter','latex');
    ylabel('$u^{(2)}$','Interpreter','latex');
    axis equal;    
    hold off;
    
    subplot(447)
    hold on;
    if k >2
        plot(prc_states.PEu(1,2:k),prc_states.PEu(2,2:k),'k*','LineWidth',1,'markersize',3);
        hold on;
    end
    xlim([PEuMin PEuMax]);
    ylim([PEuMin PEuMax]);
    title('Prediction Error $PE_u$','Interpreter','latex');
    xlabel('$u^{(1)}$','Interpreter','latex');
    ylabel('$u^{(2)}$','Interpreter','latex');
    axis equal;
    hold off;
    
    % Plot prediction error PE1    
    subplot(444)
    hold on;
    plot(prc_states.PE1(1,1:k),prc_states.PE1(2,1:k),'g*','LineWidth',1,'markersize',3)
    xlim([PE1Min PE1Max]);
    ylim([PE1Min PE1Max]);
    title('Prediction Error $PE_1$','Interpreter','latex');
    xlabel('$x_1^{(1)}$','Interpreter','latex');
    ylabel('$x_1^{(2)}$','Interpreter','latex');
    axis equal;
    hold off;
    
    %plot uncertainties or conffidence epllisoid volume
    subplot(4,4,[9,10])
    plot( 0:k-1,prc_states.EllipsoidVol(1,1:k),'g-' );
    %hold on;
    %plot( 0:k-1,prc_states.EllipsoidVol(2,1:k),'r-' );
    hold on;
    %plot( 1:k,prc_states.EllipsoidVol(3,1:k),'b-' ); 
    %hold on;
    xlim([0 T-1]);
    ylim([min(min(prc_states.EllipsoidVol(1,2:end))),...
        max(max(prc_states.EllipsoidVol(1,2:end)))]);
    %h=legend('$x_1-space$','$x_2-space$','FontSize',10);
    %set(h,'Interpreter','latex','Location','Northeast') %设置legend为latex解释器显示分式
    title(['Uncertainties and Confidence Epllisoid Volume in ' '$x_1$' '-space'],'Interpreter','latex' )
    hold off;
    
    %plot uncertainties or conffidence epllisoid volume
    subplot(4,4,[11,12])
    %plot( 0:k-1,prc_states.EllipsoidVol(1,1:k),'g-' );
    %hold on;
    plot( 0:k-1,prc_states.EllipsoidVol(2,1:k),'r-' );
    hold on;
    %plot( 1:k,prc_states.EllipsoidVol(3,1:k),'b-' ); 
    %hold on;
    xlim([0 T-1]);
    ylim([min(min(prc_states.EllipsoidVol(2,:)))-0.01...
        max(max(prc_states.EllipsoidVol(2,:))) + ...
       0.3*(max(max(prc_states.EllipsoidVol(2,:))) -   min(min(prc_states.EllipsoidVol(2,:))) )]);
    title(['Uncertainties and Confidence Epllisoid Volume in ' '$x_2$' '-space'],'Interpreter','latex' )
    hold off;
    
    % Plot Prediction Prior correlation 
    subplot(4,4,[13,14])
    %plot( 0:k-1,  postcorr1(1:k) ,'g-' );
    %hold on;
    if k > 1
        plot(1:k-1, priorcorr1(2:k), 'r-');
        hold on;
    end
    %plot(0,postcorr1(1) ,'go' );
    xlim([0 T-1]);
    %ylim([-1 2.2]);
    ylim( [min(priorcorr1)   max(priorcorr1)+1e-6 ])
    title('Prediction Prior Correlation between $\mu^{(1)}_1$ and $\mu^{(2)}_1$ ','Interpreter','latex')
    %h=legend('Posterior-$\rho_1$','Prior-$\rho_1$','FontSize',10);
    %set(h,'Interpreter','latex','Location','Northeast') %设置legend为latex解释器显示分式
    
    % Plot Posterior correlation
    subplot(4,4,[15,16])
    plot( 0:k-1,  postcorr1(1:k) ,'g-' );
    hold on;
    plot(0,postcorr1(1) ,'go' );
    xlim([0 T-1]);
    ylim( [min(postcorr1)   max(postcorr1)+1e-6 ])
    title('Posterior Correlation between $\mu^{(1)}_1$ and $\mu^{(2)}_1$','Interpreter','latex')
    
    frame=getframe(gcf);
    imind=frame2im(frame);
    [imind,cm] = rgb2ind(imind,256);
    if k==1
         imwrite(imind,cm,'./sim1/ghbf_mvBern_sim1.gif','gif', 'Loopcount',inf,'DelayTime',0.2);
    else
         imwrite(imind,cm,'./sim1/ghbf_mvBern_sim1.gif','gif','WriteMode','append','DelayTime',0.2);
    end
    pause(0.02)
 %save as a picture
 %print('-depsc','-painters','-r300','./miniGHBF_multiGauss_u2DPlotTraj.eps')
 %print( '-dpng','-painters','-r300','./miniGHBF_multiGauss_u2DPlotTraj.png')
 %saveas(gcf,'./miniGHBF_multiGauss_u2DPlotTraj.pdf')
 end 
 end


 
 
