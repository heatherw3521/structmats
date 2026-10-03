"""Create an Octave-runnable COPY of the repo's @hss and +hssutil folders in OUTDIR.

Only needed to run the minnorm_paper suite in GNU Octave (>= 8). MATLAB users do not
need this.  The repo itself is never modified.  Changes made in the copy:
  * hss_constructor.m, inter_decompv3.m: the `arguments` blocks are replaced by
    explicit name/value parsing (Octave cannot parse validation functions), and
    calls of the form f(..., name = value) become f(..., 'name', value);
  * inter_decompv3.m: qr(A,'econ','vector') -> qr(A,0) (Octave syntax, same result);
  * +hssutil/sparsesign.m: plain-MATLAB stand-in for the MEX file (same distribution).
Add  octave_compat/  (for dictionary.m) and OUTDIR to the Octave path.
Usage:  python3 make_octave_copy.py /path/to/structmats /path/to/OUTDIR
"""
import re, shutil, sys, pathlib

repo, outdir = pathlib.Path(sys.argv[1]), pathlib.Path(sys.argv[2])
outdir.mkdir(parents=True, exist_ok=True)
for d in ['@hss', '+hssutil']:
    if (outdir / d).exists():
        shutil.rmtree(outdir / d)
    shutil.copytree(repo / d, outdir / d)
here = pathlib.Path(__file__).parent
shutil.copy(here / 'sparsesign_octave.m', outdir / '+hssutil' / 'sparsesign.m')
kw = ['ctype', 'cval', 'orientation', 'oversample', 'powerits', 'escalatemargin', 'escalatepowerits', 'sketcher']

def namevalue(s):
    for w in kw:
        s = re.sub(r"(?<![\.\w])%s\s*=\s*" % w, "'%s', " % w, s)
    return s

p = outdir / '@hss/private/hss_constructor.m'
s = p.read_text()
s = s.replace('function H = hss_constructor(H,A,options)', 'function H = hss_constructor(H,A,varargin)')
i0 = s.index('    arguments'); i1 = s.index('    end', i0) + len('    end')
s = s[:i0] + """
    options = struct('blocksize',200,'sizeA',size(A),'cutrule',@(k) ceil(k/2),'tol',1e-12, ...
        'k',0,'decomp','ID','powerits',0,'escalatemargin',10,'escalatepowerits',1);
    for ii = 1:2:numel(varargin)
        options.(varargin{ii}) = varargin{ii+1};
    end
""" + s[i1:]
p.write_text(namevalue(s))

p = outdir / '+hssutil/inter_decompv3.m'
s = p.read_text()
s = s.replace('function [Z, rows, varargout] = inter_decompv3(A, options)', 'function [Z, rows, varargout] = inter_decompv3(A, varargin)')
i0 = s.index('\narguments'); i1 = s.index('\nend', i0) + len('\nend')
s = s[:i0] + """
options = struct('ctype','threshold','cval',1e-12,'orientation',0,'oversample',10, ...
    'powerits',0,'sketcher','sparsestack','escalatemargin',10,'escalatepowerits',1);
for ii = 1:2:numel(varargin)
    options.(varargin{ii}) = varargin{ii+1};
end
""" + s[i1:]
s = re.sub(r"qr\(([^,]*),'econ','vector'\)", r"qr(\1,0)", s)
s = re.sub(r"qr\(([^,]*),'econ'\)", r"qr(\1,0)", s)
p.write_text(namevalue(s))
print('Octave copy written to', outdir)
