%% Parameters
variables = struct(...
    'NBASIS_STAT', {{nan}}, ...
    'NBASIS_EIGS', {{nan}}, ...
    'NBASIS_ADJT', {{nan}}, ...
    'CA_START', {{inf}}, ...
    'CHOL_UPDATE_PACE', {{nan}}, ...
    'domain', {{'column', 'spire', 'twobar'}}, ...
    'helem', {{0.5}}, ...
    'f', {{2.0}}, ...
    'rmin', {{3.00}}, ...
    'filename', {{'results.mat'}});
jobdir = 'src/simulate';
jobname = 'minCstVBLF_study.m';
jobfarmParentdir = 'processed_data';
jobfarmdirNamebase = 'ref';
generate_jobs(jobdir, jobname, jobfarmParentdir, jobfarmdirNamebase, variables);