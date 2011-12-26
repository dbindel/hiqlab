function p = hiqpath(varargin)

global HIQ_PATH;
if nargin > 0
  HIQ_PATH = [];
  for k = 1:nargin
    if ~isempty(varargin{k})
      if ~ischar(varargin{k}), error('Path component must be a string'); end
      if isempty(HIQ_PATH)
        HIQ_PATH = varargin{k};
      else
        HIQ_PATH = [HIQ_PATH ';' varargin{k}];
      end
    end
  end
end
p = HIQ_PATH;
