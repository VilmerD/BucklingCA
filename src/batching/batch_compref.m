%% Parameters
variables = struct(...
    'domain',           {{'twobar'}}, ...
    'helem',            {{0.5}}, ...
    'nevals',           {{8}}, ...
    'minloop',          {{50}}, ...
    'maxloop',          {{100}}, ...
    'pphysstart',       {{3}}, ...
    'pphysend',         {{6}}, ...
    'dpphys',           {{0.5}}, ...
    'pNstart',          {{8}}, ...
    'pNend',            {{8}}, ...
    'betastart',        {{6}}, ...
    'betaend',          {{6}}, ...
    'contstart',        {{20}}, ...
    'mpara',            {{[2e-1, 2e5, 0.30]}}, ...
    'rmin',             {{3.00}}, ...
    'xi',               {{1}});

jobdir = 'src/simulate';
jobname = 'minCstV.m';
jobfarmParentdir = 'processed_data';
jobfarmdirNamebase = 'compref';
generate_jobs(jobdir, jobname, jobfarmParentdir, jobfarmdirNamebase, variables);