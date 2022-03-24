clc;
clear all;
close all;
%% Task 1 Frequency

% Data Task 1
angle_task1 = [0 45 90];            % in degrees
freq = [8 9 10 11 12 13 14 15 16];  % in kHz
peaks_0 = [2.09 2.12 2.77 3.12 3.19 ...
    3.05 2.79 2.295 1.455];         % peak to peak average at 0 degrees
peaks_45 = [0.855 1.085 1.33 1.75 2.165 ...
    2.385 2.35 1.985 1.42];         % peak to peak average at 45 degrees
peaks_90 = [1.455 1.76 2.17 2.485 2.59 ...
    2.61 2.35 2.035 1.625];         % peak to peak average at 90 degrees

% Data Task 2
angle_task2 = [90 75 60 45 30 15 0 -15 ...
    -30 -45 -60 -75 -90];               % in degrees
peaks_8 = [1.655 1.59 1.41 1.37 1.385 1.45 1.52 1.47 ...
    1.415 1.415 1.47 1.57 1.555];       % peak to peak average at 8 kHz
peaks_12 = [2.66 2.635 2.385 2.07 2.02 2.89 3.115 2.69 ...
    2.045 2.34 2.585 2.52 2.27];        % peak to peak average at 12 kHz
peaks_16 = [1.46 1.385 1.285 0.83 0.775 1.89 2.15 1.33 ...
    0.655 1.06 1.385 1.535 1.495];      % peak to peak average at 16 kHz

%% Normalizing the data
max_task1 = max([peaks_0 peaks_45 peaks_90]);
max_task2 = max([peaks_8 peaks_12 peaks_16]);

% Normalized Data Task 1
peaks_0_norm = db(peaks_0./max_task1);
peaks_45_norm = db(peaks_45./max_task1);
peaks_90_norm = db(peaks_90./max_task1);

% Normalized Data Task 2
peaks_8_norm = db(peaks_8./max_task2);
peaks_12_norm = db(peaks_12./max_task2);
peaks_16_norm = db(peaks_16./max_task2);

%% Frequency response for 3 angles
figure()
plot(freq, peaks_0_norm)
hold on;
plot(freq, peaks_45_norm)
plot(freq, peaks_90_norm)
hold off;
grid on;
ylabel('Amplitude in dB')
xlabel('Frequency in kHz')
legend('Angle 0°','Angle 45°','Angle 90°',Location='southeast')
title('Frequency Response at different angles')

% %% Beam Pattern
% figure()
% plot(angle_task2, peaks_8_norm)
% hold on;
% plot(angle_task2, peaks_12_norm)
% plot(angle_task2, peaks_16_norm)
% hold off;
% grid on;
% ylabel('Amplitude in dB')
% xlabel('Angle in degrees')
% xticks(-90:15:90)
% legend('Frequency 8kHz','Frequency 12kHz','Frequency 16kHz', ...
%     Location='southoutside')
% title('Beam Pattern at different frequencies')
% 
% 
% figure()
% polarplot(angle_task2./(180*pi), peaks_8_norm)
% hold on;
% polarplot(angle_task2./(180*pi), peaks_12_norm)
% polarplot(angle_task2./(180*pi), peaks_16_norm)
% legend('Frequency 8kHz','Frequency 12kHz','Frequency 16kHz', ...
%     Location='southoutside')
% title('Beam Pattern at different frequencies2')


%% Directivity
Theta = angle_task2.*(pi/180);      % angle in rad
Theta_right = Theta(1:7);          % angles at the right side
Theta_left = Theta(7:13);          % angles at the left side
d = 15*(pi/180);                    % d(Theta) in rad

% Directivity at 8kHz
% peak-to-peak amplitude at angles on each side
B_8_right = peaks_8_norm(1:7);
B_8_left = peaks_8_norm(7:13);
% denominator to integrate at each side
integral_8_right = sum(B_8_right.^2 .* sin(Theta_right).*d);
integral_8_left = sum(B_8_left.^2 .* sin(Theta_left).*d);
D_8_right = 2./ integral_8_right;   % Directivity D at right quadrant
D_8_left = 2./ integral_8_left;     % Directivity D at left quadrant
DI_8_right = 10*log10(D_8_right);   % DI at right quadrant
DI_8_left = 10*log10(D_8_left);     % DI at left quadrant

% Directivity at 12kHz
% peak-to-peak amplitude at angles on each side
B_12_right = peaks_12_norm(1:7);
B_12_left = peaks_12_norm(7:13);
% denominator to integrate at each side
integral_12_right = sum(B_12_right.^2 .* sin(Theta_right).*d);
integral_12_left = sum(B_12_left.^2 .* sin(Theta_left).*d);
D_12_right = 2./ integral_12_right;   % Directivity D at right quadrant
D_12_left = 2./ integral_12_left;     % Directivity D at left quadrant
DI_12_right = 10*log10(D_12_right);   % DI at right quadrant
DI_12_left = 10*log10(D_12_left);     % DI at left quadrant

% Directivity at 16kHz
% peak-to-peak amplitude at angles on each side
B_16_right = peaks_16_norm(1:7);
B_16_left = peaks_16_norm(7:13);
% denominator to integrate at each side
integral_16_right = sum(B_16_right.^2 .* sin(Theta_right).*d);
integral_16_left = sum(B_16_left.^2 .* sin(Theta_left).*d);
D_16_right = 2./ integral_16_right;   % Directivity D at right quadrant
D_16_left = 2./ integral_16_left;     % Directivity D at left quadrant
DI_16_right = 10*log10(D_16_right);   % DI at right quadrant
DI_16_left = 10*log10(D_16_left);     % DI at left quadrant

% Print matrix of Ds and DIs per frequency
D_left = [D_8_left D_12_left D_16_left]
D_right = [D_8_right D_12_right D_16_right]
DI_left = [DI_8_left DI_12_left DI_16_left]
DI_right = [DI_8_right DI_12_right DI_16_right]

%% polardb
theta = angle_task2;                % angles in degrees
rho8 = peaks_8_norm;                % plot value in dB at 8kHz
rho12 = peaks_12_norm;              % plot value in dB at 12kHz
rho16 = peaks_16_norm;              % plot value in dB at 16kHz
lim = -10;                          % lower limit for plot in dB
NN = 5;                             % resolution in magnitude in dB
line_style = '-g';                  % string indicating line style

figure()
subplot(1,2,1)
first = polardb(theta,rho16,lim,NN,line_style)
hold on;
sec = polardb(theta,rho12,lim,NN,line_style)
thi = polardb(theta,rho8,lim,NN,line_style)
hold  off;
legend('','','','',Location='best')
set(gca,'FontSize',12,'FontWeight','bold')
title('Beam Pattern')
% hpol=polardb(Theta,rho,lim,NN,line_style);

% function hpol=polardb(theta,rho,lim,NN,line_style)
% 
% const=1;
% theta=theta/180*pi;
% if nargin < 4
% 	error('Requires 4 or 5 input arguments.')
% elseif nargin == 4 
%         [n,n1]=size(theta);
% 	if isstr(rho)
% 		line_style = rho;
% 		rho = theta;
% 		[mr,nr] = size(rho);
% 		if mr == 1
% 			theta = 1:nr;
% 		else
% 			th = (1:mr)';
% 			theta = th(:,ones(1,nr));
% 		end
% 	else
% 		line_style = 'auto';
% 	end
% elseif nargin == 3
% 	line_style = 'auto';
% 	rho = theta;
% 	[mr,nr] = size(rho);
% 	if mr == 1
% 		theta = 1:nr;
% 	else
% 		th = (1:mr)';
% 		theta = th(:,ones(1,nr));
% 	end
% end
% if isstr(theta) | isstr(rho)
% 	error('Input arguments must be numeric.');
% end
% if any(size(theta) ~= size(rho))
% 	error('THETA and RHO must be the same size.');
% end
% % NN=1;
% nr = size(rho,2);
% tck = -floor(lim/NN);
% lim = -tck*NN;
% I=find(rho<lim);
% ni=size(I,2);
% rho(I)=lim*ones(1,ni);
% rho = rho/NN+tck*ones(1,nr);
% 
% % get hold state
% cax = newplot;
% next = lower(get(cax,'NextPlot'));
% hold_state = ishold;
% 
% % get x-axis text color so grid is in same color
% tc = get(cax,'xcolor');
% 
% % Hold on to current Text defaults, reset them to the
% % Axes' font attributes so tick marks use them.
% fAngle  = get(cax, 'DefaultTextFontAngle');
% fName   = get(cax, 'DefaultTextFontName');
% fSize   = get(cax, 'DefaultTextFontSize');
% fWeight = get(cax, 'DefaultTextFontWeight');
% set(cax, 'DefaultTextFontAngle',  get(cax, 'FontAngle'), ...
% 	'DefaultTextFontName',   get(cax, 'FontName'), ...
% 	'DefaultTextFontSize',   get(cax, 'FontSize'), ...
% 	'DefaultTextFontWeight', get(cax, 'FontWeight') )
% 
% % only do grids if hold is off
% if ~hold_state
% 
% % make a radial grid
% 	hold on;
% 	hhh=plot([0 max(theta(:))],[0 max(abs(rho(:)))]);
% 	v = [get(cax,'xlim') get(cax,'ylim')];
% 	ticks = length(get(cax,'ytick'));
% 	delete(hhh);
% % check radial limits and ticks
% 	rmin = 0; rmax = v(4); rticks = ticks-1;
% 	if rticks > 5	% see if we can reduce the number
% 		if rem(rticks,2) == 0
% 			rticks = rticks/2;
% 		elseif rem(rticks,3) == 0
% 			rticks = rticks/3;
% 		end
% 	end
% 
% % define a circle
% 	th = 0:pi/50:2*pi;
% 	xunit = cos(th);
% 	yunit = sin(th);
% % now really force points on x/y axes to lie on them exactly
%     inds = [1:(length(th)-1)/4:length(th)];
%     xunits(inds(2:2:4)) = zeros(2,1);
%     yunits(inds(1:2:5)) = zeros(3,1);
%     if ~isstr(get(cax,'color')),
%        patch('xdata',xunit*rmax,'ydata',yunit*rmax, ...
%              'edgecolor',tc,'facecolor',get(gca,'color'),...
%              'handlevisibility','off');
%     end
% %	rinc = (rmax-rmin)/rticks;
%         rinc = const;
% %	for i=(rmin+rinc):rinc:rmax
% 	for i=[1:1:tck-1]
% %	for i=rmax:rmax
% 		plot(const*xunit*i,const*yunit*i,'--','color',tc, ...
%             'linewidth',0.5,'handlevisibility','off');
% %		text(0,i+rinc/20,['  ' num2str(10*(i-tck))],...
% % 'verticalalignment','bottom' );
% 		text(const*(-i+rinc/(1000)),0,...
%             [' ' num2str(NN*(i-tck))],'verticalalignment', ...
%             'bottom','fontsize',9 );
% 	end
% plot(const*xunit*tck,const*yunit*tck,'-','color',tc, ...
%     'linewidth',0.5,'handlevisibility','off');
% text(const*(-tck+rinc/(1000)),0,[' ' num2str(NN*(tck-tck))], ...
%     'verticalalignment','bottom' ,'fontsize',9);
% % plot spokes
% 	th = (1:6)*2*pi/12;
% %	th = (1:2)*2*pi/4;
% 	cst = cos(th); snt = sin(th);
% 	cs = [-cst; cst];
% 	sn = [-snt; snt];
% 	plot(const*rmax*cs,const*rmax*sn,'--','color',tc,'linewidth',0.5);
%        
% % annotate spokes in degrees
% %	rt = 1.1*rmax;
% 	rt = 1.15*rmax;
% 	for i = 1:max(size(th))
% 	        text(const*rt*snt(i),const*rt*cst(i),int2str(i*30), ...
%                 'horizontalalignment','center' );
% 	        
%   		loc = int2str(i*30-180);
% %		if i == max(size(th))
% %			loc = int2str(0);
% % 		end
% 		text(-const*rt*snt(i),-const*rt*cst(i),loc, ...
%             'horizontalalignment','center' );
% 	end
% 
% % set viewto 2-D
% 	view(0,90);
% % set axis limits
% 	axis(rmax*[-1 1 -1.1 1.1]);
% end
% 
% % Reset defaults.
% set(cax, 'DefaultTextFontAngle', fAngle , ...
% 	'DefaultTextFontName',   fName , ...
% 	'DefaultTextFontSize',   fSize, ...
% 	'DefaultTextFontWeight', fWeight );
% 
% % transform data to Cartesian coordinates.
% yy = const*rho.*cos(theta);
% xx = const*rho.*sin(theta);
% 
% % plot data on top of grid
% if strcmp(line_style,'auto')
% 	q = plot(xx,yy);
% else
% 	q = plot(xx,yy,line_style);
% end
% set(q,'LineWidth',1.0);
% if nargout > 0
% 	hpol = q;
% end
% if ~hold_state
% 	axis('equal');axis('off');
% end
% 
% % reset hold state
% if ~hold_state, set(cax,'NextPlot',next); 
% end
% 
% end