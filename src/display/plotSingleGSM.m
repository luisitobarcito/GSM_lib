function plotSingleGSM(gsm_COV, nbrs, c_idx, varSZ, covSZ, only_center)%, h, is_complex, varSZ, covSZ)


if ~exist('varSZ' , 'var')
    varSZ = 10;
end
if ~exist('covSZ' , 'var')
    covSZ = 10;
end

if ~exist('only_center', 'var')
    only_center = false;
end

if ~exist('c_idx', 'var')
    COV = gsm_COV/max(diag(gsm_COV));
else
    COV = gsm_COV/gsm_COV(c_idx, c_idx);
end 

xlimval = max(abs(nbrs(:,4)))*[-1.5 1.5];
ylimval = max(abs(nbrs(:,3)))*[-1.5 1.5];

hold on;
neigh_colors = 'rgbkymc';
main_marker = 'o';
neigh_markers = 'sdt*+'; 
for i = 1 : size(COV,1)
    for j = i+1:size(COV,1)
        ind1 = i;
        ind2 = j;
        x1 = nbrs(ind1,4);
        y1 = nbrs(ind1,3);
        x2 = nbrs(ind2,4);
        y2 = nbrs(ind2,3);
        if  abs(COV(ind1,ind2)) > eps
            if COV(ind1,ind2) > 0
              plot([x1, x2], [y1, y2],'-', ...
		   'LineWidth',covSZ*abs(COV(ind1,ind2)/sqrt(COV(i,i)*COV(j,j))), ...
		   'Color', [0.5 0 0.5] + [0.5 0 -0.5]*abs(COV(ind1,ind2)/sqrt(COV(i,i)*COV(j,j))));
            else
              plot([x1, x2], [y1, y2],'-', ...
		   'LineWidth',covSZ*abs(COV(ind1,ind2)/sqrt(COV(i,i)*COV(j,j))), ...
		   'Color', [0.5 0 0.5] + [-0.5 0 0.5]*abs(COV(ind1,ind2)/sqrt(COV(i,i)*COV(j,j))));              
            end
            %     set(gca,'YDir','reverse')
        end
    end
end
for i=1:size(COV,1)
    ind = i;
    x = nbrs(ind,4);
    y = nbrs(ind,3);
    if only_center == true && (x == 0  && y == 0) && (i ~= c_idx) 
        continue
    else
        plot(x, y,'ko', 'MarkerSize',varSZ*COV(i,i), 'LineWidth', 3);
    end
    set(gca,'YDir','reverse')
    % OS just reversed? Changed to cos first
end
set(gca, 'Visible', 'off')
hold off
axis square; xlim(xlimval);ylim(ylimval);
