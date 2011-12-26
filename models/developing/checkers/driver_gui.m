function driver_gui(selector, mesh, Mk, Kk, Pk, Lk, f, H)

if nargin == 0, selector = 0; end
fscale = 2*pi*1e6;

switch selector
case 0 % Initialize GUI application

  fig = figure( 'CloseRequestFcn','driver_gui(4)', 'Resize','off');

  info.mhz   = f(1) / fscale;
  info.mesh  = mesh;
  info.f     = f;
  info.H     = H;
  info.Mk    = Mk;
  info.Kk    = Kk;
  info.Pk    = Pk;
  info.Lk    = Lk;

  info.u     = Lk * ((Kk-(info.mhz*fscale)^2*Mk)\Pk);
  info.u     = info.u * 4e-6/max(abs(info.u));
  Mesh_set_u(info.mesh, info.u);
  plotmesh_bode(info.mesh, info.f, info.H, info.mhz*fscale);

  info.editbox = uicontrol(...
    'Style','edit',...
    'Position',[10 10 60 20],...
    'String',num2str(info.mhz),...
    'CallBack','driver_gui(3)');

  info.slider1=uicontrol(gcf,...
    'Style','slider',...
    'Min' ,f(1)/fscale,'Max',f(end)/fscale, ...
    'Position',[70,10,200,20], ...
    'Value', info.mhz,...
    'SliderStep',[0.01 0.1], ...
    'BackgroundColor',[0.8,0.8,0.8],...
    'CallBack', 'driver_gui(5)');

  info.t = 0;

  set(fig,'UserData',info, 'HandleVisibility','callback')
	
case 3 % Enter value in edit box

  fig = gcf;
  info = get(fig,'UserData');
  newcount = sscanf(get(info.editbox,'String'),'%d');
  if ~isempty(newcount)
    info.mhz = newcount;
  end
  set(info.editbox,'String',num2str(info.mhz))
  set(info.slider1,'Value',info.mhz)
  set(fig,'UserData',info)
  updategraph(info, fscale);

case 5 % Slider control

  fig = gcf;
  info = get(fig,'UserData');
  newcount = get(info.slider1,'Value');
  if ~isempty(newcount)
    info.mhz= newcount;
  end
  set(info.editbox,'String',num2str(info.mhz))
  set(fig,'UserData',info)
  updategraph(info, fscale);

case 4 % Close requested

  fig = gcf;
  info = get(fig,'UserData');
  delete(gcf)
	
end


function updategraph(info, fscale)

Lk  = info.Lk;
Mk  = info.Mk;
Kk  = info.Kk;
Pk  = info.Pk;
info.u   = Lk * ((Kk-(info.mhz*fscale)^2*Mk)\Pk);
info.u   = info.u * 4e-6/max(abs(info.u));
children = get(gcf, 'Children');
children = setdiff(children, [info.editbox, info.slider1]);
delete(children);
Mesh_set_u(info.mesh, info.u);
opt.clf = 0;
plotmesh_bode(info.mesh, info.f, info.H, info.mhz*fscale, opt);

