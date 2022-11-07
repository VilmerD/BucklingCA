% Solves all the reference problems
domains = {'column', 'spire', 'twobar'};
for k = 1:3
    clearvars -except domains k
    close all
    domain = domains{k};
    minCstV_study;
end