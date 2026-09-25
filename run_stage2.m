% run_stage2.m
%
% PURPOSE: Runs connectome_harmonics.m unattended and logs its output and
%          figures.
% ROLE:    Optional helper for pipeline step 1.
%
% Stage 2 is slow and memory-hungry, so this runs it without supervision:
% console output goes to stage2_log.txt, figures are saved as PNGs, and
% the log records whether it finished or where it failed.
diary('stage2_log.txt');
try
    connectome_harmonics
    figs = findobj('Type','figure');
    for i = 1:numel(figs)
        saveas(figs(i), sprintf('stage2_fig%d.png', figs(i).Number));
    end
    disp('SUCCESS: full Stage 2 run complete');
catch ME
    fprintf('ERROR: %s\n', ME.message);
    for k = 1:numel(ME.stack)
        fprintf('  %s at line %d\n', ME.stack(k).name, ME.stack(k).line);
    end
end
diary off
