function [] = plot2curve(h,hs)
global X Y length h0

vvv=[1 .8 .7];

ax1=axes;
surfl(ax1,X,Y,hs);
shading interp;
alpha 1
axis([0 length 0 length -.25*h0 2*h0]);
view(vvv)
hold on

ax2=axes;
surfl(ax2,X,Y,h);
shading interp;
alpha 0.7
axis([0 length 0 length -.25*h0 2*h0]);
view(vvv)
hold on

linkaxes([ax1,ax2])

%% Hide the top axes
ax2.Visible = 'off';
ax2.XTick = [];
ax2.YTick = [];
%% Give each one its own colormap
colormap(ax1,'copper')
colormap(ax2,'default')

hold off
end

