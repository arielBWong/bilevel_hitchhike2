function processplot_ea(fighn, pop, funh_obj, trainy, krg)
clf(fighn);
%---
% plot EI over range

testdata = linspace(-5, 5, 100);
testdata = testdata';

[ynorm, ~, ~] = zscore(trainy);
ynorm_min = min(ynorm);
fit = EIM_eval(testdata, ynorm_min,  krg, []);
fit = -fit;
plot(testdata, fit, 'k--');


f = funh_obj(pop);
scatter(pop, f, 'ro', 'filled');

pause(1);


end