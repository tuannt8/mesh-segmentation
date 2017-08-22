%%
d = 'cement_1';

ms = log(importdata(['MS_' d '.txt'], '\n'));
variant = log(importdata(['variant_' d '.txt'], '\n'));

%%
h = figure;
hist([ms variant], 50);
legend('Mumford-Shah', 'Variance');
% title(d);
set(h, 'Position', [300,300,400,150]);