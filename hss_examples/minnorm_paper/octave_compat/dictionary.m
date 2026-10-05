classdef dictionary < handle
  % Minimal Octave stand-in for MATLAB's dictionary, covering only the API
  % used by the hss class: d({key}) = {val}; v = d({key}); isKey; isConfigured.
  % Octave only; MATLAB has dictionary built in.
  properties
    map
  end
  methods
    function obj = dictionary()
      obj.map = containers.Map('KeyType','char','ValueType','any');
    end
    function varargout = subsref(obj, s)
      switch s(1).type
        case '()'
          key = dictionary.k2s(s(1).subs{1});
          if ~isKey(obj.map, key)
            error('dictionary: key %s not found', key);
          end
          out = {obj.map(key)};
          if numel(s) > 1
            out = subsref(out, s(2:end));
          end
          varargout{1} = out;
        otherwise
          varargout{1} = builtin('subsref', obj, s);
      end
    end
    function obj = subsasgn(obj, s, val)
      switch s(1).type
        case '()'
          key = dictionary.k2s(s(1).subs{1});
          if iscell(val), val = val{1}; end
          obj.map(key) = val;
        otherwise
          obj = builtin('subsasgn', obj, s, val);
      end
    end
    function tf = isKey(obj, key)
      tf = isKey(obj.map, dictionary.k2s(key));
    end
    function tf = isConfigured(obj)
      tf = obj.map.Count > 0;
    end
  end
  methods (Static)
    function s = k2s(k)
      if iscell(k), k = k{1}; end
      s = sprintf('%d_', k);
    end
  end
end
