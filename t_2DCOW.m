clear
load gndfDataRed
Y = cat(3,X.data{:});
figure
ref = 7; sam = 1;
showCh      = 1; % which "channel" is shown 
band        = [5 10];
optimPars   = [20 10 2 1];
[xAOr,wPOr] = TwoDCOWorig(Y(1:221,:,sam(showCh)),Y(1:221,:,ref(showCh)),optimPars(1:2),optimPars(3:4));
[xAKU,wPKU] = TwoDCOW(Y(1:221,:,sam(showCh)),Y(1:221,:,ref(showCh)),optimPars(1:2),optimPars(3:4),band);
ah          = gobjects(3,1);
ph          = gobjects(3,2);
v           = 2.^linspace(18,24,7);
ah(1)       = subplot(2,2,1);
ah(2)       = subplot(2,2,2);
ah(3)       = subplot(2,2,3);
set(ah,'NextPlot','add');
[~, ph(1,2)] = contour(ah(1),Y(:,:,ref(showCh)),v,'LineColor','b');
[~, ph(2,1)] = contour(ah(2),Y(:,:,ref(showCh)),v,'LineColor','b');
[~, ph(3,1)] = contour(ah(3),Y(:,:,ref(showCh)),v,'LineColor','b');
[~, ph(2,2)] = contour(ah(1),Y(:,:,sam(showCh)),v,'LineColor','r');
[~, ph(2,2)] = contour(ah(2),xAOr,v,'LineColor','r');
[~, ph(3,2)] = contour(ah(3),xAKU(:,:,showCh),v,'LineColor','r');
ah(1).Title.String = 'Raw';
ah(2).Title.String = 'Original';
ah(3).Title.String = 'KU'; 

ah = gobjects(4,1);
figure
ah(1) = subplot(221);
imagesc(bsxfun(@minus,wPKU.xWarp,(1:size(xAKU,2)))); colorbar
ah(2) = subplot(222);
imagesc(bsxfun(@minus,wPKU.yWarp,(1:size(xAKU,1))'));colorbar
cL    = [min(ah(1).CLim(1),ah(2).CLim(1)),max(ah(1).CLim(2),ah(2).CLim(2))];
linkprop(ah(1:2),'CLim');
ah(1).CLim = cL;
set(ah(1:2),'YDir','normal');
colormap gray
ah(3) = subplot(223);
ah(3).NextPlot = 'add';
qh = quiver(wPKU.x,wPKU.y,wPKU.xR - wPKU.x,wPKU.yR - wPKU.y,'Marker','o','LineWidth',1);
axis([-1 272 -1 223])
ah(3).XTick = wPKU.xR(1,1:fix(size(wPKU.xR,2)/6):size(wPKU.xR,2));
ah(3).YTick = wPKU.yR(1:fix(size(wPKU.yR,1)/6):size(wPKU.yR,1),1);
grid on

% ph2 = plot(wPKU.x',wPKU.y','-k');
% ph3 = plot(wPKU.x,wPKU.y,'Color','k','Marker','o','MarkerFaceColor','w','MarkerEdgeColor','r');
% axis(ah(3),'tight')
ah(4) = subplot(224);
ah(4).NextPlot = 'add';
[~, ph4(1)] = contour(ah(4),Y(1:221,:,ref(showCh)),v,'LineColor','b');
[~, ph4(2)] = contour(ah(4),Y(1:221,:,sam(showCh)),v,'LineColor','r');
[~, ph4(3)] = contour(ah(4),xAKU(:,:,showCh),v,'LineColor',[0 0.8 0]);
legend(ph4,'Target','Raw','Warped')
linkaxes(ah,'xy')


Yp = permute(Y,[1 2 4 3]);