function hh = quiverc(x,y,u,v,CC,lw,color_bar_enable,color_bar_pos,vecscale,caluv)
alpha = 0.33; % Size of arrow head relative to the length of the vector
beta = 0.23;  % Width of the base of the arrow head relative to the length
plotarrows = 1; % Plot arrows



sqrt(u.^2+v.^2);
%----------------------------------------------
% Define colormap 
vr=sqrt(u.^2+v.^2);
vrn=round(vr/max(vr(:))*64);


x = x(:).';y = y(:).';
u = u(:).';v = v(:).';
vrn=vrn(:).';
uu = [x;x+u];
vv = [y;y+v];
vrn1 = vrn;
hold on
vrn1 =int8(round(vrn1));
vrn1(vrn1==0) = 1;
c = CC(vrn1,:);
vectors = plot(uu,vv,'linewidth',lw);
c = mat2cell(c,ones(1,size(c,1)),3);
set(vectors, {'color'}, c);
%----------------------------------------------
% Make arrow heads and plot them
if plotarrows
 
  hu = [x+u-alpha*(u+beta*(v+eps)) x+u-alpha*(u-beta*(v+eps));x+u x+u];
  hv = [y+v-alpha*(v-beta*(u+eps)) y+v-alpha*(v+beta*(u+eps));y+v y+v];
  vrn2= [vrn vrn];
  vrn2 = vrn2(:);

vrn2=int8(round(vrn2));
vrn2(vrn2==0) = 1;
c = CC(vrn2,:);
arrows = plot(hu,hv,'linewidth',lw);
c = mat2cell(c,ones(1,size(c,1)),3);
set(arrows, {'color'}, c);

end
if color_bar_enable == 1
    if caluv==1
        name='Velocity magnitude [pix/frame]';
    else
        name='Velocity magnitude [m/s]';
    end
    coloobj = colorbar(color_bar_pos,'FontWeight','bold','Fontsize',12,'color','white');
    caxis([0 max(vr(:))*caluv/vecscale]);
    if strcmp(color_bar_pos,'East')==1 || strcmp(color_bar_pos,'West')==1
%         set(coloobj,'YTickLabel',num2str(get(coloobj,'YTick')','%5.5g'))
        ylabel(coloobj,name,'fontsize',9,'fontweight','normal');
    end
    if strcmp(color_bar_pos,'North')==1 || strcmp(color_bar_pos,'South')==1
%         set(coloobj,'XTickLabel',num2str(get(coloobj,'XTick')','%5.5g'))
        xlabel(coloobj,name,'fontsize',9,'fontweight','normal');
    end
end
hh = [vectors ; arrows];
