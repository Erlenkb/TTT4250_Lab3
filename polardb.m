%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% polardb.m
% modified from Matlab's polar.m by K. Bell
% Last updated 8/30/00 by K. Bell
% downloaed from http://gunston.gmu.edu/demt/oap/contents.htm
% changed 2/17/2012 by Bo Peng 
% updated 2/25/2012 by bo peng
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function hpol=polardb(theta,rho,lim,NN,line_style)
% polardb(theta,rho,lim,linestyle)
% Polardb is a modified version of the regular 'polar' function.
% Plotting range is -180 to 180 deg with zero at top.
% Theta increases in clockwise manner.
% Inputs:
%   theta - angles in degrees ( axes labeled in degrees)
%   rho   - plot value in dB
%   lim   - lower limit for plot in dB (e.g. -10)
%   NN    - resolution in magnitude (unit dB) 
%   line_style - string indicating line style, e.g. '-g' (optional)
%
%   Example polardb(theta,beampatt,-10,3,'*r')
%
%POLAR	Polar coordinate plot.
%	POLAR(THETA, RHO) makes a plot using polar coordinates of
%	the angle THETA, in radians, versus the radius RHO.
%	POLAR(THETA,RHO,S) uses the linestyle specified in string S.
%	See PLOT for a description of legal linestyles.
%	See also PLOT, LOGLOG, SEMILOGX, SEMILOGY.
%	Copyright (c) 1984-94 by The MathWorks, Inc.

const=1;
theta=theta/180*pi;
if nargin < 4
	error('Requires 4 or 5 input arguments.')
elseif nargin == 4 
        [n,n1]=size(theta);
	if isstr(rho)
		line_style = rho;
		rho = theta;
		[mr,nr] = size(rho);
		if mr == 1
			theta = 1:nr;
		else
			th = (1:mr)';
			theta = th(:,ones(1,nr));
		end
	else
		line_style = 'auto';
	end
elseif nargin == 3
	line_style = 'auto';
	rho = theta;
	[mr,nr] = size(rho);
	if mr == 1
		theta = 1:nr;
	else
		th = (1:mr)';
		theta = th(:,ones(1,nr));
	end
end
if isstr(theta) | isstr(rho)
	error('Input arguments must be numeric.');
end
if any(size(theta) ~= size(rho))
	error('THETA and RHO must be the same size.');
end
% NN=1;
nr = size(rho,2);
tck = -floor(lim/NN);
lim = -tck*NN;
I=find(rho<lim);
ni=size(I,2);
rho(I)=lim*ones(1,ni);
rho = rho/NN+tck*ones(1,nr);

% get hold state
cax = newplot;
next = lower(get(cax,'NextPlot'));
hold_state = ishold;

% get x-axis text color so grid is in same color
tc = get(cax,'xcolor');

% Hold on to current Text defaults, reset them to the
% Axes' font attributes so tick marks use them.
fAngle  = get(cax, 'DefaultTextFontAngle');
fName   = get(cax, 'DefaultTextFontName');
fSize   = get(cax, 'DefaultTextFontSize');
fWeight = get(cax, 'DefaultTextFontWeight');
set(cax, 'DefaultTextFontAngle',  get(cax, 'FontAngle'), ...
	'DefaultTextFontName',   get(cax, 'FontName'), ...
	'DefaultTextFontSize',   get(cax, 'FontSize'), ...
	'DefaultTextFontWeight', get(cax, 'FontWeight') )

% only do grids if hold is off
if ~hold_state

% make a radial grid
	hold on;
	hhh=plot([0 max(theta(:))],[0 max(abs(rho(:)))]);
	v = [get(cax,'xlim') get(cax,'ylim')];
	ticks = length(get(cax,'ytick'));
	delete(hhh);
% check radial limits and ticks
	rmin = 0; rmax = v(4); rticks = ticks-1;
	if rticks > 5	% see if we can reduce the number
		if rem(rticks,2) == 0
			rticks = rticks/2;
		elseif rem(rticks,3) == 0
			rticks = rticks/3;
		end
	end

% define a circle
	th = 0:pi/50:2*pi;
	xunit = cos(th);
	yunit = sin(th);
% now really force points on x/y axes to lie on them exactly
    inds = [1:(length(th)-1)/4:length(th)];
    xunits(inds(2:2:4)) = zeros(2,1);
    yunits(inds(1:2:5)) = zeros(3,1);
    if ~isstr(get(cax,'color')),
       patch('xdata',xunit*rmax,'ydata',yunit*rmax, ...
             'edgecolor',tc,'facecolor',get(gca,'color'),...
             'handlevisibility','off');
    end
%	rinc = (rmax-rmin)/rticks;
        rinc = const;
%	for i=(rmin+rinc):rinc:rmax
	for i=[1:1:tck-1]
%	for i=rmax:rmax
		plot(const*xunit*i,const*yunit*i,'--','color',tc,'linewidth',0.5,'handlevisibility','off');
%		text(0,i+rinc/20,['  ' num2str(10*(i-tck))],'verticalalignment','bottom' );
		text(const*(-i+rinc/(1000)),0,[' ' num2str(NN*(i-tck))],'verticalalignment','bottom','fontsize',9 );
	end
plot(const*xunit*tck,const*yunit*tck,'-','color',tc,'linewidth',0.5,'handlevisibility','off');
text(const*(-tck+rinc/(1000)),0,[' ' num2str(NN*(tck-tck))],'verticalalignment','bottom' ,'fontsize',9);
% plot spokes
	th = (1:6)*2*pi/12;
%	th = (1:2)*2*pi/4;
	cst = cos(th); snt = sin(th);
	cs = [-cst; cst];
	sn = [-snt; snt];
	plot(const*rmax*cs,const*rmax*sn,'--','color',tc,'linewidth',0.5);
       
% annotate spokes in degrees
%	rt = 1.1*rmax;
	rt = 1.15*rmax;
	for i = 1:max(size(th))
	        text(const*rt*snt(i),const*rt*cst(i),int2str(i*30),'horizontalalignment','center' );
	        
  		loc = int2str(i*30-180);
%		if i == max(size(th))
%			loc = int2str(0);
% 		end
		text(-const*rt*snt(i),-const*rt*cst(i),loc,'horizontalalignment','center' );
	end

% set viewto 2-D
	view(0,90);
% set axis limits
	axis(rmax*[-1 1 -1.1 1.1]);
end

% Reset defaults.
set(cax, 'DefaultTextFontAngle', fAngle , ...
	'DefaultTextFontName',   fName , ...
	'DefaultTextFontSize',   fSize, ...
	'DefaultTextFontWeight', fWeight );

% transform data to Cartesian coordinates.
yy = const*rho.*cos(theta);
xx = const*rho.*sin(theta);

% plot data on top of grid
if strcmp(line_style,'auto')
	q = plot(xx,yy);
else
	q = plot(xx,yy,line_style);
end
set(q,'LineWidth',1.0);
if nargout > 0
	hpol = q;
end
if ~hold_state
	axis('equal');axis('off');
end

% reset hold state
if ~hold_state, set(cax,'NextPlot',next); end





