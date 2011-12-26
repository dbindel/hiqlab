HOME='/Users/dbindel/work/hiqlab-new'
function addhpath(p)  addpath(HOME .. p)  end

addhpath('/src/lua/?');
addhpath('/src/lua/?.lua');
addhpath('/models/?');
addhpath('/models/?.lua');
addhpath('/models/util/?');
addhpath('/models/util/?.lua');
addpath('?');
addpath('?.lua');
