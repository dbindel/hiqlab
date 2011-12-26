% function plot_bode(freq, H, opts)
%
% Plot a Bode plot for the given frequency range and H.
% Options are:
%   usehz     ( 0 ) - Assume input frequency is in Hz (otherwise radians/s)
%   logf      ( 0 ) - Use log frequency for x axis
%   magnitude ( 0 ) - Plot magnitude only
%   visualQ   ( 0 ) - Visually determine Q value
%   lstyle    ('b') - Line style

function plot_bode(freq, H, opt)

% HiQLab
% Copyright (c): Regents of the University of California

if nargin < 2,            error('Need both freq and H');  end
if nargin < 3,            opt = [];                       end

usehz     = qoptdefault(opt, 'usehz',     0 );
logf      = qoptdefault(opt, 'logf',      0 );
magnitude = qoptdefault(opt, 'magnitude', 0 );
visualQ   = qoptdefault(opt, 'visualQ',   0 );
lstyle    = qoptdefault(opt, 'lstyle',   'b');
holdit    = ishold;

if ~usehz,      freq = freq / 2 / pi;  end
if ~magnitude,  subplot(2,1,1);        end
if holdit,      hold on;               end

dB = 20*log10(abs(H));
if logf
  semilogx( freq, dB, lstyle );
else
  plot( freq, dB, lstyle );
end
xlabel('Frequency (Hz)');
ylabel('Transfer (dB)');

if ~magnitude
  subplot(2,1,2);
  if holdit, hold on; end
  if logf
    semilogx( freq, angle(H)*180/pi, lstyle );
  else
    plot( freq, angle(H)*180/pi, lstyle );
  end
  xlabel('Frequency (Hz)');
  ylabel('Phase (degrees)');
end

if visualQ

  if ~magnitude, subplot(2,1,1); end
  hold on;

  lh1 = line([freq(1), freq(end)], [max(dB), max(dB)]);
  lh2 = line([freq(1), freq(end)], [max(dB), max(dB)]-3);
  set(lh1, 'LineStyle', '--', 'Color', 'r');
  set(lh2, 'LineStyle', '--', 'Color', 'r');
  xlabel('Frequency (Hz)');
  ylabel('Transfer (dB)');
  [dBmax, Imax] = max(dB);
  I3dB = find(dB >= dBmax-3);
  title(sprintf('f0 = %g, Q = %g\n', freq(Imax), ...
                freq(Imax)/( freq(I3dB(end)) - freq(I3dB(1)) ) ));

  hold off;

end

if holdit, hold on; end
