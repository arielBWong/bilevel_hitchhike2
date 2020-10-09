figure(1);
x = 1:10
y1 = log(x)
y2 = log10(x);
plot(x, y1); hold on;
plot(x, y2); hold on;


x = [x, fliplr(x)];
y = [y1, fliplr(y2)];

fill(x, y, 'g');