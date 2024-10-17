% x1= [0, 0.866];
% x2=[1,-1];
% x3=[-1,-1];
% e=0.1;
% d1=[e/(2*x1(1)^2+x1(2)), e/(6*x1(2)^2-1.5)];
% d2=[e/(2*x2(1)^2+x1(2)), e/(6*x2(2)^2-1.5)];
% d3=[e/(2*x3(1)^2+x1(2)), e/(6*x3(2)^2-1.5)];
% 
% [X,Y]=meshgrid(-1.5:0.1:1.5, -1.5:0.1:1.5);
% 
% Z =  -log(1+0.5679*(1/sqrt(2*pi*d1(1))*exp(-(X-x1(1)).^2/(2*d1(1)))).*(1/sqrt(2*pi*d1(2))*exp(-(Y-x1(2)).^2/(2*d1(2)))) + 0.2222*(1/sqrt(2*pi*d2(1))*exp(-(X-x2(1)).^2/(2*d2(1)))).*(1/sqrt(2*pi*d2(2))*exp(-(Y-x2(2)).^2/(2*d2(2)))) + 0.2099*(1/sqrt(2*pi*d3(1))*exp(-(X-x3(1)).^2/(2*d3(1)))).*(1/sqrt(2*pi*d3(2))*exp(-(Y-x3(2)).^2/(2*d3(2)))));
% %% 绘制潜势图谱
% figure;
% surfc(X, Y, Z);
% xlabel('CHEK1');
% ylabel('P53');

[X,Y] = meshgrid(0:0.01:3, 3:0.01:8);
vertices = [2.09610945745910,2,1.90000000000000,1.85424639310745,1.80000000000000,1.70000000000000,1.60776649648111,1.60000000000000,1.50000000000000,1.49516842924028,1.41191186023602,1.40000000000000,1.35649380390473,1.32490689600531,1.31599275948493,1.31835463945444,1.33461276017114,1.38217669851717,1.40000000000000,1.40194542504546,1.40000000000000,1.30000000000000,1.20000000000000,1.15113109748833,1.10000000000000,1.01925710851787,1,0.934801155541044,0.901597200455951,0.900000000000000,0.859314148746271,0.845751434629639,0.827562922838718,0.800000000000000,0.700000000000000,0.668011318261220,0.600000000000000,0.527664018396188,0.500000000000000,0.436976395686492,0.400000000000000,0.397993067721465,0.338023058030910,0.319730751068124,0.316203233656117,0.323197781674981,0.349462353254393,0.400000000000000,0.406193963975808,0.464892040487186,0.500000000000000,0.564577706940749,0.600000000000000,0.700000000000000,0.732909456462341,0.800000000000000,0.900000000000000,1,1.10000000000000,1.20000000000000,1.27974173594017,1.30000000000000,1.40000000000000,1.45227468891528,1.50000000000000,1.55337199636759,1.60000000000000,1.60000000000003,1.66646277246979,1.68629009979088,1.69160848260707,1.68945086559864,1.68304701552352,1.70000000000000,1.73222018420470,1.80000000000000,1.85980752076285,1.90000000000000,1.94001162203375,1.98683326545417,2,2.00513120402143,2.02274617741634,2.00716990808662,2,1.98942630602442,1.98097915561640,2,2.10000000000000,2.20000000000000,2.26485945361721,2.30000000000000,2.40000000000000,2.40715964737049,2.49666495395582,2.50000000000000,2.57194166576335,2.59561249260112,2.60000000000000,2.61851531559254,2.62527146537037,2.60000000000000,2.59997788155227,2.58611594152067,2.53900837507808,2.50000000000000,2.46546371730647,2.40000000000000,2.34620149241047,2.30000000000000,2.20000000000000,2.10000000000000,2.09610945745910;3.90000000000000,3.87577662915638,3.88280131479451,3.90000000000000,3.90472418394430,3.92881083976050,4,4.00385273486922,4.09499335350114,4.10000000000000,4.20000000000000,4.23689874446545,4.30000000000000,4.40000000000000,4.50000000000000,4.60000000000000,4.70000000000000,4.80000000000000,4.88369332275121,4.90000000000000,4.90123972678489,4.92323379428242,4.96377408394794,5,5.02422184436094,5.10000000000000,5.12351142170099,5.20000000000000,5.30000000000000,5.31027202429083,5.40000000000000,5.50000000000000,5.60000000000000,5.61669585338562,5.67415747495324,5.70000000000000,5.73195129832038,5.80000000000000,5.83048183869495,5.90000000000000,5.99734818805794,6,6.10000000000000,6.20000000000000,6.30000000000000,6.40000000000000,6.50000000000000,6.57650157435951,6.60000000000000,6.70000000000000,6.73760337696478,6.80000000000000,6.83266676842901,6.89197961475306,6.90000000000000,6.94332319320734,6.97139569735412,6.97847408036956,6.97330648747747,6.94962266016063,6.90000000000000,6.89546259628882,6.84754644431066,6.80000000000000,6.75664608608514,6.70000000000000,6.60000000000004,6.60000000000000,6.50000000000000,6.40000000000000,6.30000000000000,6.20000000000000,6.10000000000000,6.03710136912035,6,5.95777110490617,5.90000000000000,5.84734005005073,5.80000000000000,5.70000000000000,5.61120118642450,5.60000000000000,5.50000000000000,5.40000000000000,5.38220825823767,5.30000000000000,5.20000000000000,5.19249871182765,5.17687454194564,5.14479334655988,5.10000000000000,5.08865358552378,5.00690870121911,5,4.90000000000000,4.89327373426607,4.80000000000000,4.70000000000000,4.64944009475821,4.60000000000000,4.50000000000000,4.40043664699412,4.40000000000000,4.30000000000000,4.20000000000000,4.15346205928389,4.10000000000000,4.03456984853665,4,3.96131767552236,3.91410012566502,3.90019488108157,3.90000000000000];
vertices = vertices';
[in,on] = inpolygon(X,Y,vertices(:,1),vertices(:,2));
Y1=Y(in==1);
X1=X(in==1);
x1 = [3.08400000000000,1.09900000000000,1.09940000000000,1.09530000000000,6.28470000000000,1.14220000000000,1.09940000000000,1.00720000000000,5.12540000000000,3.12370000000000,3.03360000000000,0.164300000000000,1.13710000000000,0];
x2 = [3.15950000000000,1.09900000000000,1.09920000000000,1.09570000000000,5.49760000000000,1.14280000000000,1.09920000000000,1.44000000000000,5.12590000000000,3.12410000000000,3.03380000000000,0.164000000000000,1.13660000000000,0.432900000000000];
x3 = [3.18550000000000,1.09900000000000,1.09850000000000,1.09580000000000,4.52740000000000,1.14440000000000,1.09850000000000,1.97320000000000,5.12620000000000,3.12460000000000,3.03410000000000,0.164400000000000,1.13630000000000,0.966100000000000];
e=0.03;
a=1.1;
b=2;
k=1;
n=3;
s=0.5;
d1 = [e,e,e,e,e,e,e,e,e,e/(k-b*s^n*(n*x1(10)^(n-1)/(s^n+x1(10)^n)^2)),e/(k-b*s^n*(n*x1(11)^(n-1)/(s^n+x1(11)^n)^2)),e, e, e/(k-a*(n*x1(14)^(n-1)*(s^n+x1(14)^n)-n*x1(14)^n*x1(14)^(n-1))/(s^n+x1(14)^n)^2)];
d2 = [e,e,e,e,e,e,e,e,e,e/(k-b*s^n*(n*x2(10)^(n-1)/(s^n+x2(10)^n)^2)),e/(k-b*s^n*(n*x2(11)^(n-1)/(s^n+x2(11)^n)^2)),e, e, e/(k-a*(n*x2(14)^(n-1)*(s^n+x2(14)^n)-n*x2(14)^n*x2(14)^(n-1))/(s^n+x2(14)^n)^2)];
d3 = [e,e,e,e,e,e,e,e,e,e/(k-b*s^n*(n*x3(10)^(n-1)/(s^n+x3(10)^n)^2)),e/(k-b*s^n*(n*x3(11)^(n-1)/(s^n+x3(11)^n)^2)),e, e, e/(k-a*(n*x3(14)^(n-1)*(s^n+x3(14)^n)-n*x3(14)^n*x3(14)^(n-1))/(s^n+x3(14)^n)^2)];
x4 = [3.08400000000000,1.09900000000000,1.09940000000000,1.09530000000000,6.28470000000000,1.14220000000000,1.09940000000000,1.00720000000000,5.12540000000000,3.12370000000000,3.03360000000000,0.164300000000000,1.13710000000000,0];
x5 = [3.15970000000000,1.09900000000000,1.09920000000000,1.09570000000000,5.49430000000000,1.14280000000000,1.09920000000000,1.44180000000000,5.12590000000000,3.12410000000000,3.03380000000000,0.164000000000000,1.13660000000000,0.433900000000000];
x6 = [3.18550000000000,1.09900000000000,1.09850000000000,1.09580000000000,4.52530000000000,1.14440000000000,1.09850000000000,1.97440000000000,5.12620000000000,3.12460000000000,3.03410000000000,0.164400000000000,1.13630000000000,0.969300000000000];
d4 = [0.0300000000000000,0.0300000000000000,0.0300000000000000,0.0300000000000000,0.0301000000000000,0.0300000000000000,0.0300000000000000,0.0300000000000000,0.0300000000000000,0.0302000000000000,0.0303000000000000,0.0300000000000000,0.0300000000000000,0.0300000000000000];
d5 = [0.0300000000000000,0.0300000000000000,0.0300000000000000,0.0300000000000000,0.0300000000000000,0.0300000000000000,0.00330000000000000,0.0300000000000000,0.0300000000000000,0.0302000000000000,0.0303000000000000,0.0300000000000000,0.0300000000000000,0.0300000000000000];
d6 = [0.0300000000000000,0.0300000000000000,0.0300000000000000,0.0300000000000000,0.0300000000000000,0.0300000000000000,0.0349000000000000,0.0300000000000000,0.0300000000000000,0.0302000000000000,0.0303000000000000,0.0300000000000000,0.0300000000000000,0.0300000000000000];
Z1 = -log(0.6153*(1/sqrt(2*pi*d1(8))*exp(-(X1-x1(8)).^2/(2*d1(8)))).*(1/sqrt(2*pi*d1(5))*exp(-(Y1-x1(5)).^2/(2*d1(5)))) + 0.0779*(1/sqrt(2*pi*d2(8))*exp(-(X1-x2(8)).^2/(2*d2(8)))).*(1/sqrt(2*pi*d2(5))*exp(-(Y1-x2(5)).^2/(2*d2(5)))) + 0.3066*(1/sqrt(2*pi*d3(8))*exp(-(X1-x3(8)).^2/(2*d3(8)))).*(1/sqrt(2*pi*d3(5))*exp(-(Y1-x3(5)).^2/(2*d3(5)))));
Z2 = -log(0.4015*(1/sqrt(2*pi*d1(8))*exp(-(X1-x1(8)).^2/(2*d1(8)))).*(1/sqrt(2*pi*d1(5))*exp(-(Y1-x1(5)).^2/(2*d1(5)))) + 0.0774*(1/sqrt(2*pi*d2(8))*exp(-(X1-x2(8)).^2/(2*d2(8)))).*(1/sqrt(2*pi*d2(5))*exp(-(Y1-x2(5)).^2/(2*d2(5)))) + 0.521*(1/sqrt(2*pi*d3(8))*exp(-(X1-x3(8)).^2/(2*d3(8)))).*(1/sqrt(2*pi*d3(5))*exp(-(Y1-x3(5)).^2/(2*d3(5)))));
Z = Z1 - Z2;
T= delaunay(X1, Y1);
% 提取三角形的顶点索引
points = [X1,Y1];

% 提取三角形顶点
triPoints = points(T, :);
% 计算三角形边长
edges = [T(:,[1 2]); T(:,[2 3]); T(:,[3 1])];
edgeLengths = sqrt(sum((triPoints(edges(:,1),:) - triPoints(edges(:,2),:)).^2, 2));

% 设置边长阈值
edgeLengthThreshold = 0.3; % 这个值可以根据你的需求进行调整

validTri = [];

% 遍历每个三角形
for i = 1:size(T, 1)
    % 提取三角形的边长
    triangleEdges = edgeLengths(T(i,:));
    point1 = points(T(i,1),:) ;
    point2 = points(T(i,2),:) ;
    point3 = points(T(i,3),:) ;
    edge12 = sqrt((point1(1) - point2(1))^2 + (point1(2) - point2(2))^2);
    edge13 = sqrt((point1(1) - point3(1))^2 + (point1(2) - point3(2))^2);
    edge23 = sqrt((point2(1) - point3(1))^2 + (point2(2) - point3(2))^2);
    dis = edge12 + edge13 + edge23;
    % 如果所有边长都小于阈值，则保留该三角形
    if dis < edgeLengthThreshold
        validTri = [validTri; T(i,:)];
    end
end

% 可选：重新三角化剩余的点
% 如果删除了很多三角形，你可能希望重新三角化剩余的点集以得到一个更合理的剖分
% 注意：这一步可能不是必须的，取决于你的应用需求
% remainingPoints = points(unique(validTri), :);
% reTri = delaunay(remainingPoints(:,1), remainingPoints(:,2));

figure;
trisurf(validTri, X1,Y1, Z);
% tricontour(validTri, X1,Y1, Z,5);
hold on;
scatter3(x1(8), x1(5),-0.426902420359625, 200, 'MarkerEdgeColor','k', 'MarkerFaceColor','r');
hold on;
scatter3(x2(8), x2(5),-0.006439746977691, 200, 'MarkerEdgeColor','k', 'MarkerFaceColor','g');
hold on;
scatter3(x3(8), x3(5),0.530206075315647, 200, 'MarkerEdgeColor','k', 'MarkerFaceColor','y');
xlabel('CHEK1');
ylabel('P53');
% trimesh(T, X1,Y1, Z);
% tricontour(T, X1,Y1, Z);

function [c,h] = tricontour(tri,x,y,z,nv)
%TRICONTOUR Triangular Contour Plot.
% TRICONTOUR(TRI,X,Y,Z,N) draws scalar N contour lines treating the values
% in Z as heights above a plane. TRI,X,Y,and Z define a triangulation where
% the triangles are defined by the M-by-3 face matrix TRI, such as that
% returned by DELAUNAY. Each row of TRI contains indices into the X,Y, and
% Z vertex vectors to define a single triangular face. Contours are
% computed directly from the triangulation rather than interpolating back
% to a cartesian grid using GRIDDATA.
% TRICONTOUR(TRI,X,Y,Z,V) draws length(V) contour lines at the values
% specified in vector V.
% TRICONTOUR(TRI,X,Y,Z,[v v]) draws a single contour line at the level v.
%
% [C,H] = TRICONTOUR(...) returns contour matrix C as described in CONTOURC
% and a vector of handles H to the created patch objects.
% H can be used to set patch properties.
% CLABEL(C) or CLABEL(C,H) labels the contour levels.
%
% Example:
%           x=linspace(-3,3,39);
%           y=linspace(-2.5,2.5,49);
%           [xx,yy]=meshgrid(x,y);
%           zz=peaks(xx,yy);
%           v=-3:1:5; % contour levels
%           subplot(1,2,1)
%           [C,h]=contour(xx,yy,zz,v);   % standard contour for comparison
%           clabel(C)
%           title Contour
% 
%           idx=randperm(numel(zz));     % grab some scattered indices
%           n=idx(1:ceil(numel(zz)/2))'; % one half of them
%           x=xx(n);                     % get scattered data
%           y=yy(n);
%           z=zz(n);
%           tri=delaunay(x,y);           % triangulate scattered data
%           subplot(1,2,2)
%           [C,h]=tricontour(tri,x,y,z,v);
%           clabel(C,h)
%           title TriContour
%
% view(3) displays the contour in 3-D.
%
% See also DELAUNAY, CONTOUR, TRIMESH, TRISURF, TRIPLOT, PATCH.
% D.C. Hanselman, University of Maine, Orono, ME 04469
% MasteringMatlab@yahoo.com
% Mastering MATLAB 7
% 2006-05-07, 2006-05-16, 2006-07-25
if nargin<5
	error('Not Enough Input Arguments.')
end
x=x(:);	% convert input data into column vectors
y=y(:);
z=z(:);
xlen=length(x);
if ~isequal(xlen,length(y),length(z))
   error('X, Y, and Z Must Have the Same Number of Elements.')
end
if size(tri,2)~=3 || any(tri(:)<0) || any(tri(:)>xlen)
   error('TRI Must Be a Valid Triangulation of the Data in X, Y, Z.')
end
zs=z(tri);
zmax=max(max(zs));              % find max and min in z data that is in tri
zmin=min(min(zs));
if length(nv)==1                                 % nv is number of contours
   zlev=linspace(zmax,zmin,nv+2);
elseif length(nv)==2 && nv(1)==nv(2)              % nv is one contour level
   zlev=nv(1);
else                                       % nv is vector of contour levels
   zlev=sort(nv,'descend');
end
zlev(zlev>=zmax | zlev<=zmin)=[];  % eliminate contours outside data limits
nlev=length(zlev);
if nlev==0
   error('No Contours to Plot. Chosen Contours Outside Limits of Data.')
end
% precondition the input data
[zs,zidx]=sort(zs,2);         % sort vertices by z value ascending
for k=1:size(zs,1)            % shuffle triangles to match
   tri(k,:)=tri(k,zidx(k,:));
end
hax=newplot;                  % create new axis if needed
h=[];                         % patch handle storage
C=zeros(2,0);                 % Clabel data storage
cs=[2 1];                     % column swap vector cs(1)=2, cs(2)=1;
% Main Loop ---------------------------------------------------------------
for v=1:nlev                  % one contour level at a time
   zc=zlev(v);                % chosen level
   above=zs>=zc;              % true for vertices above given contour
   numabove=sum(above,2);     % number of triangle vertices above contour
   tri1=tri(numabove==1,:);   % triangles with one vertex above contour
   tri2=tri(numabove==2,:);   % triangles with two vertices above contour
   n1=size(tri1,1);           % number with one vertex above
   n2=size(tri2,1);           % number with two vertices above
   edge=[tri1(:,[1 3])        % first column is indices below contour level
         tri1(:,[2 3])        % second column is indices above contour level
         tri2(:,[1 2])
         tri2(:,[1 3])];
   if n1==0                   % assign edges to triangle number
      n=[1:n2 1:n2]';
   elseif n2==0
      n=[1:n1 1:n1]';
   else
      n=[1:n1 1:n1 n1+(1:n2) n1+(1:n2)]';
   end
   [edge,idx]=sortrows(edge);    % put shared edges next to each other
   n=n(idx);                     % shuffle triangle numbers to match
   idx=all(diff(edge)==0,2);     % find shared edges
   idx=[idx;false]|[false;idx];  % True for all shared edges
   
   % eliminate redundant edges, two triangles per interior edge
   edgeh=edge(~idx,:);           % hull edges
   nh=n(~idx);                   % hull triangle numbers
   if ~isempty(nh)
      nh(end,2)=0;               % zero second column for hull edges
   end
   edges=edge(idx,:);            % shared edges
   edges=edges(1:2:end-1,:);     % take only unique edges
   ns=n(idx);                    % interior triangle numbers
   ns=[ns(1:2:end) ns(2:2:end)]; % second column is second triangle
   edge=[edgeh;edges];           % unique edges
   nn=[nh;ns];                   % two columns of triangle numbers
   ne=size(edge,1);              % number of edges
   
   flag=true(ne,2);              % true for each unused edge per triangle
   tmp=zeros(ne+1,1);            % contour data temporary storage
   
   xe=x(edge);                   % x values at vertices of edges
   ye=y(edge);                   % y values at  vertices of edges
   ze=z(edge);                   % z data at  vertices of edges
   alpha=(zc-ze(:,1))./(ze(:,2)-ze(:,1)); % interpolate all edges
   xc=alpha.*(xe(:,2)-xe(:,1)) + xe(:,1); % x values on this contour
   yc=alpha.*(ye(:,2)-ye(:,1)) + ye(:,1); % y values on this contour
   while any(flag)	% while there are still unused edges -----------------
      
      xtmp=tmp;
      ytmp=tmp;
      [ir,ic]=find(flag,1);            % find next unused edge
      flag(ir,ic)=false;               % mark this edge used
      
      k=1;                             % first data point in subcontour
      xtmp(k)=xc(ir);                  % store data from this edge
      ytmp(k)=yc(ir);
      
      while true     % complete this subcontour ---------------------------
         
         [ir,ic]=find(flag&nn(ir,ic)==nn,1);% find other edge of triangle
         flag(ir,ic)=false;            % mark this edge used
         k=k+1;
         xtmp(k)=xc(ir);               % store data from this edge
         ytmp(k)=yc(ir);
         
         ic=cs(ic);                    % other triangle that shares edge
         if nn(ir,ic)==0               % reached hull, subcontour complete
            k=k+1;
            xtmp(k)=nan;               % don't let subcontour close
            ytmp(k)=nan;
            break
         elseif ~flag(ir,ic)           % complete closed subcontour
            break
         else                          % more points remain on subcontour
            flag(ir,ic)=false;         % mark this edge used
         end
      end % while true ----------------------------------------------------
      xtmp(k+1:end)=[];                % throw away unused storage
      ytmp(k+1:end)=[];                % xtmp,ytmp contain subcontour
      
      if nargout<2                     % plot the subcontour
         patch('XData',xtmp,'YData',ytmp,'CData',repmat(zc,k,1),...
               'Parent',hax,'FaceColor','none','EdgeColor','flat',...
               'UserData',zc)
         C=horzcat(C,[zc xtmp';k ytmp']); % contour label data
      else                             % plot subcontour and create output
         h=[h;patch('XData',xtmp,'YData',ytmp,'CData',repmat(zc,k,1),...
         'Parent',hax,'FaceColor','none','EdgeColor','flat',...
         'UserData',zc)]; %#ok
         C=horzcat(C,[zc xtmp';k ytmp']); % contour label data
      end
   end % while any(flag) --------------------------------------------------
end % for v=1:nlev
if nargout
   c=C;
end
end

