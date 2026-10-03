function mp_write_table(file, S, cols, note)
%MP_WRITE_TABLE  Write a struct array as a plain-text CSV table.
%   mp_write_table(file, S, cols, note) writes two '#' comment lines (note,
%   and the MATLAB/Octave version, platform, thread count and date), a header
%   line with the field names cols, and one line per element of S. Numbers are
%   written with 6 significant digits, strings as they are (no commas).
d = fileparts(file);
if ~isempty(d) && ~exist(d, 'dir'), mkdir(d); end
fid = fopen(file, 'w');
if fid < 0, error('mp_write_table: cannot open %s', file); end
try, thr = sprintf('%d', maxNumCompThreads); catch, thr = 'n/a'; end
fprintf(fid, '# %s\n', note);
fprintf(fid, '# %s | %s | threads %s | %s\n', version, computer, thr, datestr(now));
fprintf(fid, '%s\n', strjoin(cols, ','));
for i = 1:numel(S)
    v = cell(1, numel(cols));
    for j = 1:numel(cols)
        x = S(i).(cols{j});
        if ischar(x)
            v{j} = strrep(x, ',', ';');
        elseif isempty(x)
            v{j} = 'NaN';
        else
            v{j} = sprintf('%.6g', x);
        end
    end
    fprintf(fid, '%s\n', strjoin(v, ','));
end
fclose(fid);
fprintf('wrote %s (%d rows)\n', file, numel(S));
end
