function mp_save(resdir, name, S)
%MP_SAVE  save a results struct with machine information.
if ~exist(resdir, 'dir'), mkdir(resdir); end
S.machine.version = version;
S.machine.computer = computer;
try, S.machine.threads = maxNumCompThreads; catch, S.machine.threads = NaN; end
S.machine.date = datestr(now);
save(fullfile(resdir, [name '.mat']), '-struct', 'S');
end
