% plotmesh_bode(mesh, f, H, fcurrent, opt)
%
% Simultaneously show a deformed mesh and a Bode plot.
% f is the frequency list, H is the transfer function,
% and fcurrent is the marked current frequency.

function plotmesh_bode(mesh, f, H, fcurrent, opt)

if nargin < 5, opt = []; end
do_clf = qoptdefault(opt, 'clf', 1);

if do_clf, clf; end
f = f / 2/pi;
fcurrent = fcurrent/2/pi;

axes('position', [0, 0.2, 0.5, 0.8]);
opt.forces = [];
opt.deform = 1;
opt.clf = 0;
plotmesh(mesh,opt);
axis equal

subplot(2,2,2);
plot( f, 20*log10(abs(H)) );
xlabel('Frequency (Hz)');
ylabel('Amplitude (dB)');
v = axis;
l = line([fcurrent fcurrent], [v(3) v(4)]);
set(l, 'Color', 'r');

subplot(2,2,4);
plot( f, angle(H) );
xlabel('Frequency (Hz)');
ylabel('Phase (rad)');
v = axis;
l = line([fcurrent fcurrent], [v(3) v(4)]);
set(l, 'Color', 'r');
