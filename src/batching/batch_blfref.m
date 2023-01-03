%% Parameters
variables = struct(...
    'domain',           {{'column', 'spire', 'twobar'}}, ...
    'helem',            {{0.50}}, ...
    'nevals',           {{6.00}}, ...
    'sf',               {{3.00}}, ...
    'minloop',          {{350}}, ...
    'maxloop',          {{400}}, ...
    'pphysstart',       {{2}}, ...
    'pphysend',         {{6}}, ...
    'dpphys',           {{0.25}}, ...
    'contpace',         {{20}}, ...
    'contstart',        {{30}}, ...
    'pNstart',          {{8}}, ...
    'pNend',            {{16}}, ...
    'dpN',              {{0.50}}, ...
    'betastart',        {{6}}, ...
    'betaend',          {{6}}, ...
    'dbeta',            {{0}}, ...
    'betapace',         {{1}}, ...
    'NBASIS_STAT',      {{nan}}, ...
    'NBASIS_EIGS',      {{nan}}, ...
    'NBASIS_ADJT',      {{nan}}, ...
    'CA_START',         {{inf}}, ...
    'CHOL_UPDATE_PACE', {{nan}}, ...
    'CA_ANG_TOL',       {{1.00}}, ...
    'CAorthotype',      {{'none'}}, ...
    'SHOW_DESIGN',      {{1}}, ...
    'mpara',            {{[2e-1, 2e5, 0.30]}}, ...
    'rmin',             {{3.00}}, ...
    'xi',               {{1}}, ...
    'x0strat',          {{'comp'}}, ...
    'mmalen',           {{0.200}}, ...
    'PROFILE',          {{true}}, ...
    'SAVEX',            {{true}});

jobdir = 'src/simulate';
jobname = 'minCstVBLF.m';
jobfarmParentdir = 'processed_data';
jobfarmdirNamebase = 'ref';
generate_jobs(jobdir, jobname, jobfarmParentdir, jobfarmdirNamebase, variables);
