%{
Author: S. Watanabe
Last updated; 7/22/2021

This program is for the measurement of the size and volume of dendritic spines.


Usage

Segment 1
Load a bright field image (data_in; either color or B/W image). 
Set the pixel size. 
The pixel values will be inverted for ease of calculation.

Segment 2
Set a background ROI near the target spine. 
This uses the ginput function: click at the edges of the ROI to be drawn near the spine, and hit the Enter key to complete.
The background ROI will be displayed in the image window. 
The background intensity will be displayed in the command window.

Segment 3
Set a spine ROI surrounding the target spine. 
This uses the ginput function: click at the edges of the ROI to be drawn that surrounds the spine, and hit the Enter key to complete.
The spine ROI will be displayed in the image window. 
The intensity (sum of the ROI value minus background) will be displayed in the command window. 

Segment 4
Set a line crossing the spine. 
This uses the ginput function: click at the ends of the line to be drawn, and hit the Enter key to complete.
The line will be displayed in the image window. 
The profile and Gaussian fitting will be displayed in a graph window.
The spine width (2 times sigma of the Gaussian fit) will be displayed in the command window.

%}

%% Segment 1 - import data

pixel_size = 0.036;

%for 8 bit color brightfield image
data = 255 - double(data_in(:, :, 2));
data_disp = data_in / 255;

%for 12 bit B/W brightfield image
% data = 4096 - double(data_in);
% data_disp = zeros(size(data_in, 1), size(data_in, 2), 3);
% data_disp(:, :, 1) = data_in / 4096;
% data_disp(:, :, 2) = data_disp(:, :, 1);
% data_disp(:, :, 3) = data_disp(:, :, 1);


figure(11); clf; hold on;
image(data_disp);


%% Segment 2 - background level

[xb,yb] = ginput;

for i=1:length(xb)
	plot([xb; xb(1)], [yb; yb(1)], 'r');
end

[x,y] = meshgrid(1:size(data,2), 1:size(data,1));
pat_b = inpolygon(x,y, [xb; xb(1)], [yb; yb(1)])';
val_b = sum(sum(data .* pat_b')) / sum(sum(pat_b));

fprintf('background = %.2f\n', val_b);

%% Segment 3 - spine intensity

[xa,ya] = ginput;

for i=1:length(xa)
	plot([xa; xa(1)], [ya; ya(1)], 'r');
end

[x,y] = meshgrid(1:size(data,2), 1:size(data,1));
pat_a = inpolygon(x,y, [xa; xa(1)], [ya; ya(1)])';
val_a = sum(sum(data .* pat_a'));
area_a = sum(sum(pat_a));

sp_intensity = val_a - area_a * val_b;

fprintf('spine intensity = %.0f\n', sp_intensity);

%% Segment 4 - spine profile
[xs,ys] = ginput;
plot(xs, ys, 'r');

linelength = sqrt(diff(xs)^2 + diff(ys)^2);
nlinept = floor(linelength) +1;
lineptwt = zeros(2,2, nlinept);

if abs(diff(xs)) >0
	sl = diff(ys) / diff(xs);
	dx = 1 / sqrt(1 + sl^2);
	dy = sl / sqrt(1 + sl^2);
else
	sl = Inf;
	dx = 0;
	dy = 1;
end

if abs(diff(xs)) >0
	if xs(1) < xs(2)
		xp = xs(1) + (0:nlinept-1)' * dx;
		yp = ys(1) + (0:nlinept-1)' * dy;
	else
		xp = xs(2) + (0:nlinept-1)' * dx;
		yp = ys(2) + (0:nlinept-1)' * dy;
	end
else
	if y(1) < y(2)
		xp = xs(1) + (0:nlinept-1)' * dx;
		yp = ys(1) + (0:nlinept-1)' * dy;
	else
		xp = xs(2) + (0:nlinept-1)' * dx;
		yp = ys(2) + (0:nlinept-1)' * dy;
	end
end

linept = round([xp, yp]);
line_prof = zeros(nlinept, 1);
for i=1:length(linept)
	line_prof(i) = data(linept(i,2), linept(i,1));
end


figure(12); clf; hold on;
plot(line_prof)

[params, sse] = gaussian_fit_fun(line_prof, 'offset', val_b, 'amp', max(line_prof));

offset = params(1);
amp = params(2);
mu = params(3);
sigma = params(4);

fplot(@(x) offset + amp * exp(-(x-mu)^2 / (2 * sigma^2)), [1, length(line_prof)], 'r');

sigma_times2_actual = 2 * sigma * pixel_size;

fprintf('spine 2*sigma = %.4f\n', sigma_times2_actual);

