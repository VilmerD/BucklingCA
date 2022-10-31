%% Plotting designs of reference problem using exact analysis
jobs = dir('processed_data/buckling_study_reference1/job*');

for k = 1:3
    results_file = fullfile(jobs(k).folder, jobs(k).name, 'results.mat');
    load(results_file, 'xPhys', 'i_img', 'j_img');
    outfile = sprintf('design%i.png', k);
    plotDesignAndSave(xPhys, i_img, j_img, outfile);
end

function plotDesignAndSave(x, i_img, j_img, outfile)
% Plots design x on i_img and j_img and saves figure
fig = figure();
ax = axes(fig, 'Position', [0, 0, 1, 1]);
axis(ax, 'tight');
axis(ax, 'off');
hold(ax, 'on');

h = max(i_img);
w = max(j_img);
fig.Position(3:4) = [w, h];

spimg = sparse(i_img, j_img, x);
imagesc(spimg);

saveas(fig, outfile);
close(fig);
end