% function [H,freq] = second_order_bode(mesh,wc,force_pat,sense_pat,opt)
%
% Plots forced response
% Output:
%   freq          - Frequency array [rad/s]
%   H             - Output
%
% Input:
%   mesh          - the mesh input
%   wc            - center frequency[rad/s]
%   force_pat     - Forcing pattern vector
%   sense_pat     - Sensing pattern vector
%   opt
%     freq          - Predefined array of freq
%     mkc(0)        - Include damping or not
%     wr_min(0.90)  - left  value mag of bode plot
%     wr_max(1.10)  - right value mag of bode plot
%     w_ndiv(50)    - number divisions in bode plot
%     w_type('lin') - division type (linspae or logspace)
%     kmax(0)       - number of arnoldi iterations
%     realbasis(0)  - use real basis??
%     skewproj(0)   - use skew projection on complex??
%     structurep(0) - use structure preserving basis??
%     use_umfpack   - use UMFPACK?? (Default: use if exist)
%     viewmatrices(0)- view reduced matrices ?
% If kmax and w0 are given, use model reduction via an Arnoldi
% expansion about the shift w0.  If the shift frequency w0 is omitted,
% use w0 = mean(freq).

function [H,freq] = second_order_bode(mesh,wc,force_pat,sense_pat,opt)

mkc       = qoptdefault(opt, 'mkc',          0);
wr_max    = qoptdefault(opt, 'wr_max'   , 1.10);
wr_min    = qoptdefault(opt, 'wr_min'   , 0.90);
w_ndiv    = qoptdefault(opt, 'w_ndiv'   ,   50);
w_type    = qoptdefault(opt, 'w_type'   ,'lin');
kmax      = qoptdefault(opt, 'kmax'     ,    0);
opt.realbasis   = qoptdefault(opt, 'realbasis',    0);
opt.skewproj    = qoptdefault(opt, 'skewproj',     0);
opt.structurep  = qoptdefault(opt, 'structurep',   0);
opt.use_umfpack = qoptdefault(opt, 'use_umfpack',  exist('umfpack')==3);
opt.viewmatrices= qoptdefault(opt, 'viewmatrices',   0);

if strcmp(w_type,'lin')
  freq = linspace(wr_min,wr_max,w_ndiv)*wc;
elseif strcmp(w_type,'log')
  freq = logspace(wr_min,wr_max,w_ndiv)*wc;
end
if isfield('opt','freq')
  freq = opt.freq;
end
H = zeros(size(sense_pat,2),length(freq));

F = force_pat;
L = sense_pat;
if mkc
  [M,K,C,su,sf] = Mesh_assemble_mkc_nd(mesh);
  F = spdiag(1./sf) * F;
  L = spdiag(su)    * L;
  if kmax > 0
    opt.mechdof = ted_block_mesh(mesh);
    [M,C,K,L,F] = rom_soar(M,C,K,L,F, kmax, wc, opt);
  end
else
  [M,K,su,sf] = Mesh_assemble_mk_nd(mesh);
  F = spdiag(1./sf) * F;
  L = spdiag(su)    * L;
  numid = Mesh_get_numid(mesh);
  C     = spalloc(numid,numid,0);
  if kmax > 0
    [M,K,L,F] = rom_arnoldi(M,K,L,F, kmax, wc, opt);
    C         = zeros(size(M,1),size(M,2));
  end
end

if (kmax > 0)
  if mkc, rom_name = 'SOAR'; else rom_name = 'ARNOLDI'; end;
  fprintf('ROM type             ;%s\n',rom_name);
  fprintf('No Arnoldi iterations:%d\n',kmax);
  fprintf('Size of ROM          :%d\n', size(M,1));
  fprintf('        REALBASIS    :%d\n', opt.realbasis);
  fprintf('        SKEWPROJ     :%d\n', opt.skewproj);
  fprintf('        STRUCTUREP   :%d\n', opt.structurep);
  if(opt.viewmatrices==1)
    disp('Structure of reduced matrices');
    subplot(1,3,1); spy(M); title('M');
    subplot(1,3,2); spy(C); title('C');
    subplot(1,3,3); spy(K); title('K');
    pause(1);
  end
end

for idx = 1:length(freq)

  w = freq(idx);
  if (kmax == 0)
    fprintf('%d: %d Hz\n', idx, freq(idx)/2/pi);
  end

  if opt.use_umfpack & kmax==0
    H(:,idx) = L' * umfpack(K + complex(0,w)*C - w^2*M, '\', F);
  else
    H(:,idx) = L' *       ((K + complex(0,w)*C - w^2*M)  \   F);
  end

end

