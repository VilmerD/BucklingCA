%% Parameters
variables = struct(...
    'domain',           {{'column', 'spire', 'twobar'}}, ...
    'helem',            {{0.50}}, ...
    'nevals',           {{6.00}}, ...
    'sf',               {{2.50}}, ...
    'minloop',          {{320}}, ...
    'maxloop',          {{400}}, ...
    'pphysstart',       {{2}}, ...
    'pphysend',         {{6}}, ...
    'dpphys',           {{0.20}}, ...
    'pNstart',          {{8}}, ...
    'pNend',            {{24}}, ...
    'betastart',        {{1}}, ...
    'betaend',          {{1}}, ...
    'contstart',        {{20}}, ...
    'NBASIS_STAT',      {{6.00}}, ...
    'NBASIS_EIGS',      {{3.00}}, ...
    'NBASIS_ADJT',      {{6.00}}, ...
    'CA_START',         {{10.00}}, ...
    'CHOL_UPDATE_PACE', {{10.0}}, ...
    'CA_ANG_TOL',       {{1-5e-4}}, ...
    'SHOW_DESIGN',      {{1.00}}, ...
    'mpara',            {{[2e-1, 2e5, 0.30]}}, ...
    'rmin',             {{1.50, 2.00, 2.50}}, ...
    'xi',               {{1}}, ...
    'x0strat',          {{'hom'}});

jobdir = 'src/simulate';
jobname = 'minCstVBLF.m';
jobfarmParentdir = 'processed_data';
jobfarmdirNamebase = 'rom_rmin';
generate_jobs(jobdir, jobname, jobfarmParentdir, jobfarmdirNamebase, variables);