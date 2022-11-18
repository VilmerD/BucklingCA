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
    'betastart',        {{6}}, ...
    'betaend',          {{6}}, ...
    'contstart',        {{20}}, ...
    'NBASIS_STAT',      {{nan}}, ...
    'NBASIS_EIGS',      {{nan}}, ...
    'NBASIS_ADJT',      {{nan}}, ...
    'CA_START',         {{inf}}, ...
    'CHOL_UPDATE_PACE', {{nan}}, ...
    'CA_ANG_TOL',       {{1.00}}, ...
    'SHOW_DESIGN',      {{1.00}}, ...
    'mpara',            {{[2e-1, 2e5, 0.30]}}, ...
    'rmin',             {{3.00}}, ...
    'xi',               {{1}}, ...
    'x0strat',          {{'hom'}});

jobdir = 'src/simulate';
jobname = 'minCstVBLF.m';
jobfarmParentdir = 'processed_data';
jobfarmdirNamebase = 'ref';
generate_jobs(jobdir, jobname, jobfarmParentdir, jobfarmdirNamebase, variables);
