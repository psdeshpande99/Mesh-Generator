%% Domain Size %%
prompt1 = {'Enter domain radius required from centre of breakwater'};
dfin  = {'1000'};
r = inputdlg(prompt1,'Input',[1,60],dfin);
r = str2double(r);
prompt1 = {'No. of parts circumferentially (enter even no)'};
dfin  = {'48'};
a = inputdlg(prompt1,'Input',[1,60],dfin);
a = str2double(a);
prompt1 = {'No. of parts radially'};
dfin  = {'26'};
b = inputdlg(prompt1,'Input',[1,60],dfin);
b = str2double(b);
%% Breakwater size and coordinates %%
prompt1 = {'Enter number of corner points for the breakwater'};
dfin  = {'4'};
m = inputdlg(prompt1,'Input',[1,60],dfin);
m = str2double(m);
breakwater_coord = zeros(m,2);
for i=1:m
    c = 1:2;
    c = arrayfun(@num2str, c, 'uni',0);
    prompt2 = c;
    inbuiltip = zeros(1,2);
    inbuiltip = arrayfun(@num2str, inbuiltip, 'uni',0);
    dfin2  = inbuiltip;
    dlgtitle = (['Point Number',num2str(i)]);
    mat = inputdlg(prompt2,dlgtitle,[1,60],dfin2);
    breakwater_coord(i,1:2) = str2double(mat');
    breakwater_coord(i,3) = sqrt((breakwater_coord(i,1))^2+(breakwater_coord(i,2))^2);
end

%% Node value Retaintion %%
min_radius = min(abs(breakwater_coord(:,3)));
theta = linspace( pi/2, -pi/2, a+1);
radial_radius = linspace(r,min_radius,b+1);
x = radial_radius'*cos(theta);
y = radial_radius'*sin(theta);

%% Removal from Breakwater %%
xnu = [];
for i = 1:b+1
    xnu = [xnu x(i,:)];
end
ynu = [];
for i = 1:b+1
    ynu = [ynu y(i,:)];
end
[in,on] = inpolygon(xnu,ynu,breakwater_coord(:,1)',breakwater_coord(:,2)');
x_in = [];
y_in = [];
index_in = [];
index_circum = [];
index_radial = [];
for i=1:length(in)
    if(in(i)==1)
        x_in = [x_in xnu(i)];
        y_in = [y_in ynu(i)];
        index_in = [index_in i];
        if(rem(i,(a+1))==0)
            index_circum = [index_circum a+1];
            index_radial = [index_radial (floor(i/(a+1)))];
        else
            index_circum = [index_circum rem(i,a+1)];
            index_radial = [index_radial (floor(i/(a+1))+1)];
        end
    end 
end
index_circum_radial(1,:) = index_circum(1,:);
index_circum_radial(2,:) = index_radial(1,:);
% poly1 = polyshape(breakwater_coord(:,1)',breakwater_coord(:,2)');
% out1_new = [];
% for i=1:length(x_in)
%     lineseg = [0 0; x_in(i) y_in(i)];
%     [in1,out1] = intersect(poly1,lineseg);
%     out1_new = [out1_new out1];
% end

%% Point intersection with breakwater line %%
upper_angle = atand(breakwater_coord(2,2)/breakwater_coord(2,1));
lower_angle = atand(breakwater_coord(3,2)/breakwater_coord(3,1));
angles = atand(y_in./x_in);
y_extend =[];
x_extend =[];
for i = 1:length(x_in)
    if (angles(i)>upper_angle)
        y_extend = [y_extend breakwater_coord(1,2)];
        x_extend = [x_extend breakwater_coord(1,2)/tand(angles(1,i))];
    elseif(angles(i)<=upper_angle && angles(i)>=lower_angle)
        y_extend = [y_extend breakwater_coord(2,1)*tand(angles(1,i))];
        x_extend = [x_extend breakwater_coord(2,1)];
    else
        y_extend = [y_extend breakwater_coord(4,2)];
        x_extend = [x_extend breakwater_coord(4,2)/tand(angles(1,i))];
    end
end

%% Rectification %%
for i=1:length(x_extend)
    x(index_circum_radial(2,i),index_circum_radial(1,i)) = x_extend(1,i);
    y(index_circum_radial(2,i),index_circum_radial(1,i)) = y_extend(1,i);
end

%% Triangular element formation%%
nodal = zeros(3*b,4*a);
for j=1:b
    for i=1:a
        if(rem(i,2)==0)
            nodal(3*j-2,4*i-3) = x(j,i);
            nodal(3*j-2,4*i-2) = y(j,i);
            nodal(3*j-1,4*i-3) = x(j,i+1);
            nodal(3*j-1,4*i-2) = y(j,i+1);
            nodal(3*j,4*i-3)   = x(j+1,i);
            nodal(3*j,4*i-2)   = y(j+1,i);
            nodal(3*j-2,4*i-1) = x(j,i+1);
            nodal(3*j-2,4*i)   = y(j,i+1);
            nodal(3*j-1,4*i-1) = x(j+1,i+1);
            nodal(3*j-1,4*i)   = y(j+1,i+1);
            nodal(3*j,4*i-1)   = x(j+1,i);
            nodal(3*j,4*i)     = y(j+1,i);
        else
            nodal(3*j-2,4*i-3) = x(j,i);
            nodal(3*j-2,4*i-2) = y(j,i);
            nodal(3*j-1,4*i-3) = x(j+1,i+1);
            nodal(3*j-1,4*i-2) = y(j+1,i+1);
            nodal(3*j,4*i-3)   = x(j+1,i);
            nodal(3*j,4*i-2)   = y(j+1,i);
            nodal(3*j-2,4*i-1) = x(j,i);
            nodal(3*j-2,4*i)   = y(j,i);
            nodal(3*j-1,4*i-1) = x(j,i+1);
            nodal(3*j-1,4*i)   = y(j,i+1);
            nodal(3*j,4*i-1)   = x(j+1,i+1);
            nodal(3*j,4*i)     = y(j+1,i+1);
        end
    end
end

%% Plot the domain along with Breakwater %%
f = figure('Name','FEM_Breakwater','NumberTitle','off','Units','normalized','Position',[0 0 0.4 0.4]);
line(breakwater_coord(:,1),breakwater_coord(:,2),'color','k','linewidth',2);
hold on;
line([0,0],[-r,breakwater_coord(end,2)],'color','r','linewidth',2);
hold on;
line([0,0],[breakwater_coord(1,2),r],'color','r','linewidth',2);
hold on;
plot(x,y,'*');
hold on;
for i=1:b+1
    plot(x(i,:),y(i,:),'b','linewidth',0.01); 
    hold on;
end
hold on;
for i=1:b
    for j=1:a
        line(nodal(3*i-2:3*i,4*j-3),nodal(3*i-2:3*i,4*j-2),'color','b','linewidth',0.01);
        line(nodal(3*i-2:3*i,4*j-1),nodal(3*i-2:3*i,4*j),'color','b','linewidth',0.01);
    end
end
xlim([0-r*.1 r*1.1]);
ylim([-r*1.1 r*1.1]);
axis square;
