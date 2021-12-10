function saveTightFigure(h,outfilename,varargin)
% SAVETIGHTFIGURE(H,OUTFILENAME) Saves figure H in file OUTFILENAME without
%   the white space around it. 
%
% by ``a grad student"
% http://tipstrickshowtos.blogspot.com/2010/08/how-to-get-rid-of-white-margin-in.html

p = inputParser();
addParameter(p,'EqualX',false,@(x) x==0 || x==1);
addParameter(p,'XShift',0,@(x) isscalar(x) && isreal(x));
addParameter(p,'SubplotYspace',0,@(x) isscalar(x) && isreal(x)); 
addParameter(p,'SubplotXspace',0,@(x) isscalar(x) && isreal(x));
addParameter(p,'FigPosOffset',[0 0 0 0],@(x) isvector(x) && length(x)==4);
addParameter(p,'AxPosOffset',[0 0 0 0],@(x) isvector(x) && length(x)==4);
parse(p,varargin{:});


% get all axes
a = findall(h,'type','axes');

% remove axis corresponding to legends
legendIndex = zeros(length(a),1);
for i = 1:length(a)
if(strcmp(get(a(i),'Tag'),'legend'))
legendIndex(i) = 1;
end
end
a(legendIndex==1) = [];

% expand plot view
for i=1:length(a)
    ti = get(a(i),'Tightinset');
    op = get(a(i),'OuterPosition');
    newPos(i,:) = [op(1)+ti(1) op(2)+ti(2) op(3)-ti(3)-ti(1) op(4)-ti(4)-ti(2)];
end
newPos = bsxfun(@plus,newPos,p.Results.AxPosOffset(:)'); 

if p.Results.EqualX
    dx = max(newPos(:,1)) - min(newPos(:,1));
    newPos(:,1) = max(newPos(:,1)); 
    newPos(:,3) = min(newPos(:,3))-dx; 
end
newPos(:,1) = newPos(:,1)+p.Results.XShift;
if p.Results.SubplotYspace
    newPos(:,2) = newPos(:,2)-arrayfun(@(nsp) nsp*p.Results.SubplotYspace,1:length(a))'; 
end
if p.Results.SubplotXspace
    newPos(:,1) = newPos(:,1)-arrayfun(@(nsp) nsp*p.Results.SubplotXspace,length(a):-1:1)'; 
end


for i = 1:length(a)
    set(a(i),'Position',newPos(i,:));
end

% calculate papersize
set(a,'units','centimeters');
xpapermax =-inf; xpapermin=+inf;
ypapermax =-inf; ypapermin=+inf;
for i=1:length(a)
    pos = get(a(i),'Position');
    ti = get(a(i),'Tightinset');
    if( pos(1)+pos(3)+ti(1)+ti(3) > xpapermax); xpapermax = pos(1)+pos(3)+ti(1)+ti(3);end
    if( pos(1) < xpapermin); xpapermin = pos(1);end
    if( pos(2)+pos(4)+ti(2)+ti(4) > ypapermax); ypapermax = pos(2)+pos(4)+ti(2)+ti(4);end
    if( pos(2) < ypapermin); ypapermin = pos(2);end
end
paperwidth = xpapermax - xpapermin;
paperhight = ypapermax - ypapermin;

% adjust the papersize
set(h, 'PaperUnits','centimeters');
set(h, 'PaperSize', [paperwidth paperhight]);
set(h, 'PaperPositionMode', 'manual');
set(h, 'PaperPosition',[0 0 paperwidth paperhight]);
set(h,'Position',get(h,'Position')+p.Results.FigPosOffset(:)');


saveas(h,outfilename);