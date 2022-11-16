%% Parameters
variables = struct(...
    'domain', {{'column', 'spire', 'twobar'}}, ...
    'nevals', {{6}}, ...
    'helem', {{0.50}}, ...
    'sf', {{2.50}}, ...
    'rmin', {{3.00}}, ...
    'NBASIS_STAT', {{8}}, ...
    'NBASIS_EIGS', {{3}}, ...
    'NBASIS_ADJT', {{8}}, ...
    'CA_START', {{10}}, ...
    'CHOL_UPDATE_PACE', {{10}}, ...
    'CA_ANG_TOL', {{1-5e-4, 1-10e-4, 1-20e-4}}, ...
    'SHOW_DESIGN', {{1}});

jobdir = 'src/simulate';
jobname = 'minCstVBLF.m';
jobfarmParentdir = 'processed_data';
jobfarmdirNamebase = 'rom';
generate_jobs(jobdir, jobname, jobfarmParentdir, jobfarmdirNamebase, variables);