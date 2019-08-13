function [h] = plotMGSM(mgsm_model, indm, indc, inds, nbrs, h, is_complex, varSZ, covSZ)

if ~exist('h', 'var')
    h = figure;
else
    set( h)
end
if ~exist('varSZ' , 'var')
    varSZ = 10;
end
if ~exist('covSZ' , 'var')
    covSZ = 7;
end


COVc1 = mgsm_model.COVc/max(diag(mgsm_model.COVc));
COVs1 = mgsm_model.COVs/max(diag(mgsm_model.COVs));
COVcs1 = mgsm_model.COVcs/max(diag(mgsm_model.COVcs));

xlimval = max(abs(nbrs(:,4)))*[-1.5 1.5];
ylimval = max(abs(nbrs(:,3)))*[-1.5 1.5];

subplot(1,3,1); hold on;
neigh_colors = 'rgbkymc';
main_marker = 'o';
neigh_markers = 'sdt*+'; 
for i = 1:numel(indc)
    ind = indc(i);
    x = nbrs(ind,4);
    y = nbrs(ind,3);
    plot(x, y,'ko', 'MarkerSize',varSZ*COVc1(i,i));
    set(gca,'YDir', 'reverse')
end
axis square; xlim(xlimval);ylim(ylimval);
subplot(1,3,2); hold on

for i = 1 : (numel(inds) - 1)
    for j = i+1:numel(inds)
    ind1 = inds(i);
    ind2 = inds(j);
    x1 = nbrs(ind1,4);
    y1 = nbrs(ind1,3);
    x2 = nbrs(ind2,4);
    y2 = nbrs(ind2,3);
    if abs(COVs1(i,j)) > eps
      plot([x1, x2], [y1, y2],'b-', ...
	   'LineWidth', covSZ*abs(COVs1(i,j)/sqrt(COVs1(i,i)*COVs1(j,j))));
    end
%     set(gca,'YDir','reverse')
    end
end
for i = 1:numel(inds)
    ind = inds(i);
    x = nbrs(ind,4);
    y = nbrs(ind,3);
    plot(x, y,'ko', 'MarkerSize',varSZ*COVs1(i,i), 'LineWidth', 3);
    set(gca,'YDir','reverse')
end
axis square; xlim(xlimval);ylim(ylimval);
subplot(1,3,3); hold on
for i = 1 : size(COVcs1,1)
    for j = i+1:size(COVcs1,1)
    ind1 = i;
    ind2 = j;
    x1 = nbrs(ind1,4);
    y1 = nbrs(ind1,3);
    x2 = nbrs(ind2,4);
    y2 = nbrs(ind2,3);
    if  abs(COVcs1(ind1,ind2)) > eps
      plot([x1, x2], [y1, y2],'b-', ...
	   'LineWidth', covSZ*abs(COVcs1(ind1,ind2)/sqrt(COVcs1(i,i)*COVcs1(j,j))));
    end
%     set(gca,'YDir','reverse')
    end
end
for i=1:size(COVcs1,1)
    ind = i;
    x = nbrs(ind,4);
    y = nbrs(ind,3);
    plot(x, y,'ko', 'MarkerSize',varSZ*COVcs1(i,i), 'LineWidth', 3);
    set(gca,'YDir','reverse')
    % OS just reversed? Changed to cos first
end
axis square; xlim(xlimval);ylim(ylimval);
