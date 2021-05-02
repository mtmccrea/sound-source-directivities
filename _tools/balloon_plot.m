function [coefficients, D] = balloon_plot(coefficients, order, fs, sph_definition, f_to_plot)
% Creates a ballon plot from spherical harmonic coefficients
%  f_to_plot: vector frequencies to plot in Hz (default: [500 1000 2000
%                                                                  4000])

Nfft = size(coefficients, 1);

if (nargin < 5)
    % plot 4 different frequencies
    f_to_plot    = [500 1000 2000 4000]; % in Hz
end

f_to_plot    = f_to_plot(f_to_plot < fs/2); % sort out
bins_to_plot = round(f_to_plot/(fs/2) * Nfft);

% update values
f_to_plot = bins_to_plot/Nfft * fs/2;

% remove all unused data
coefficients = coefficients(bins_to_plot, :); % [NplotBins x Nsh]

[Nfreq2plot, Nsh] = size(coefficients);

% set up a spatial grid 
resolution        = 4*order; 
alpha_v           = linspace(0, 2*pi, 2*resolution);    % azimuth
beta_v            = linspace(0,   pi,   resolution+1);  % colatitude
[alpha_m, beta_m] = meshgrid(alpha_v, beta_v);

colatitude = beta_m(:).';
azimuth    = alpha_m(:).';
Npnt2plot = numel(azimuth); % number of points to plot

% compute the directivity D on the spatial grid
D = zeros(Nfreq2plot, Npnt2plot); % [NplotBins , NplotPnts]

% % original
% for l = 0 : order
%     for m = -l : l
%          D = D + repmat(coefficients(:, l^2+l+m+1), [1 size(azimuth, 2)]) .* repmat(sphharm(l, m, colatitude, azimuth, sph_definition), [size(coefficients, 1) 1]); 
%     end
% end

% % annotated version - mtm
% for l = 0 : order
%     for m = -l : l
%         idx_coeff = l^2+l+m+1;
%         coeff_lm = coefficients(:, idx_coeff);
%         % repeat each modal coeff for each plot direction
%         coeff_x_pltpnt = repmat(coeff_lm, [1 Npnt2plot]); % [NplotBins NplotPnts]
%         % SH probe of this SH mode in every plot direction
%         sphHarm_pltdir = sphharm(l, m, colatitude, azimuth, sph_definition); % [1 x NplotPnts]
%         % expand to number of plots
%         sphHarm_pltdir = repmat(sphHarm_pltdir, [Nfreq2plot 1]);  % [NplotBins x NplotPnts]
%         % scale this sampled harmonic by the modeal coefficient weight
%         mode_scaled = coeff_x_pltpnt .* sphHarm_pltdir; % [NplotBins x NplotPnts] .* [NplotBins x NplotPnts]
%         % sum this mode with the other modes for overall directivity
%         D = D + mode_scaled;
%     end
% end

% % plot-by-plot version - mtm
for bini = 1:Nfreq2plot
    for l = 0 : order
        for m = -l : l
            % select the complex weight for this mode
            idx_coeff = l^2+l+m+1;
            weight_lm = coefficients(bini, idx_coeff);
            
            % SH probe of this SH mode in every plot direction
            lm_pltdir = sphharm(l, m, colatitude, azimuth, sph_definition); % [1 x NplotPnts]
            
            % scale this sampled mode by the coefficient weight
            mode_dirs_weighted = lm_pltdir * weight_lm; % [1 x NplotPnts] .* [1]
            
            % sum this mode with the other modes for overall directivity
            D(bini, :) = D(bini, :) + mode_dirs_weighted;
        end
    end
end

% prepare data for plotting
alpha_m = reshape(azimuth,    resolution+1, []);
beta_m  = reshape(colatitude, resolution+1, []);

%  --------------------- finally, plot data -------------------------------
figure('Position', [100 100 500 500]);
set(gcf, 'Color', [1 1 1]);
Nrowcol = ceil(sqrt(Nfreq2plot));

tiledlayout(Nrowcol, Nrowcol, 'tilespacing', 'compact');
for n = 1 : Nfreq2plot
    nexttile
    plot_it(abs(reshape(D(n, :), resolution+1, [])),  alpha_m, beta_m, f_to_plot(n));
end

% for n = 1 : Nfreq2plot
%     subplot(Nrowcol, Nrowcol, n);
%     plot_it(abs(reshape(D(n, :), resolution+1, [])),  alpha_m, beta_m, f_to_plot(n));
% end

end

function plot_it(data_to_plot, alpha_m, beta_m, f)

view_angle = [40 20];

color1 = [0.9769    0.9839    0.0805];
color2 = [0.2440    0.4358    0.9988];
color3 = [0.2422    0.1504    0.6603];

% avoid a hole in the hull
alpha_m      = [alpha_m, alpha_m(:, 1)];
beta_m       = [beta_m, beta_m(:, 1)];
data_to_plot = [data_to_plot, data_to_plot(:, 1)];

[Xm, Ym, Zm] = sph2cart(alpha_m, pi/2-beta_m, data_to_plot);

surf(Xm, Ym, Zm); 

plot_max = max(abs([Xm(:); Ym(:); Zm(:)]));

% plot coordinate axes
hold on;
line([-1 1] * plot_max * .9, [.0 .0], [ 0 0], 'Marker', '.', 'LineStyle', '-', 'Color', [.5 .5 .5], 'LineWidth', 1 );
line([ 0 0], [-1  1] * plot_max * .9, [ 0 0], 'Marker', '.', 'LineStyle', '-', 'Color', [.5 .5 .5], 'LineWidth', 1 );
line([ 0 0], [ 0  0], [-1 1] * plot_max * .9, 'Marker', '.', 'LineStyle', '-', 'Color', [.5 .5 .5], 'LineWidth', 1 );
hold off;

view(view_angle(1), view_angle(2));

box on;

xlabel( '$x$', 'Interpreter', 'latex' );
ylabel( '$y$', 'Interpreter', 'latex' );
zlabel( '$z$', 'Interpreter', 'latex' );

axis equal
axis([-1 1 -1 1 -1 1] * plot_max);

lighting phong
shading interp

light('Position', [1 0 0], 'Color', color1);
light('Position', [0 -1 0], 'Color', color2);
light('Position', [0 0 1], 'Color', color3);
set(gca, 'Projection', 'Perspective', 'FontSize', 10);

title(sprintf('%d Hz', round(f)), 'fontsize', 16);

end
