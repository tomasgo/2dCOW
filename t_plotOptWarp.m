hFig = figure;
ah = axes('NextPlot','add');
p = 1;
q = diff(Diagnos.rangeP(:,:,1),1,2) + 1;
s = [0;cumsum(q)];
Pointer = Diagnos.table(1,end);
hFig.Name = sprintf('Seg100-Sl%i',6 - p);
h = gobjects(size(s,1) - 2,1);
c = 1;
bkg = NaN(max(q),length(q),2);
Table = Diagnos.table;
bkg(1,1,:) = [0,Table(2 * p,1)];
Table(p * 2:p * 2 + 1,1 + find(Table(p * 2 + 1,2:end) == 0)) = NaN;
for (i = length(s) - 1:-1:2)
    t = s(i) + 1:s(i + 1); 
    for (j = 1:length(t))
        k = t(j); 
        u = Diagnos.indexP([i,i - 1])';
        if (isnan(Table(1 + 2 * p,k))), v = [Table(1,k);NaN]; else v = [Table(1,k);Table(1,Table(1 + 2 * p,k))]; end
        w = v - u;
        bkg(j,i,1) = w(1);
        bkg(j,i,2) = Table(2 * p,k);
        ph = plot(u,w,'Color','k','LineStyle','-','Marker','o','MarkerFaceColor',[0.4 0.4 0.4]);
        if (v(1) == Pointer)
            set(ph,'Color',[.8 0 0],'MarkerFaceColor','r','LineWidth',2);
            Pointer = v(2);
            h(c) = ph;
            c = c + 1;
        end
    end
end
box on
ah.Title.String = sprintf('Segment: 100; Slack: %i',6 - p);
xlabel('Points')
ylabel('Correction [\Deltapts]')
PrintFigure(hFig,400,'',[32 20],'png')

% [~,m] = max(sum(~isnan(bkg(:,:,1))));
% [a,~] = find(bkg(:,:,1) == 0);
% d     = median(a) - a;
% for (i = 1:length(q)), bkg(:,i,:) = circshift(squeeze(bkg(:,i,:)),d(i)); end
% hPC = pcolor(Diagnos.indexP(ones(size(bkg,1),1),:),bkg(:,:,1),bkg(:,:,2));
% ah.Children = [h;setdiff(ah.Children,[hPC;h]);hPC];
% 
    