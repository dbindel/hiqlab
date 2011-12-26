addpath([pwd '/src/mfiles']);
addpath([pwd '/src/mfiles/mgraph']);
addpath([pwd '/src/mex']);
addpath([pwd '/models/util']);

hiqaddpath([pwd '/src/lua/?']);
hiqaddpath([pwd '/src/lua/?.lua']);
hiqaddpath([pwd '/models/?']);
hiqaddpath([pwd '/models/?.lua']);
hiqaddpath([pwd '/models/util/?']);
hiqaddpath([pwd '/models/util/?.lua']);
hiqaddpath('?');
hiqaddpath('?.lua');

try
  print_hiqlab_banner;
catch
  disp(' ');
  disp('Could not access the HiQLab MEX file.  You may not have the file');
  disp('compiled for your platform.');
  disp(' ');
  fprintf('Error message: %s\n', lasterr);
end
