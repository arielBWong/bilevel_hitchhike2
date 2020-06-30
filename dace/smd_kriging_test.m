%%
clc;
clear all;
train_x = csvread('train_x.csv');
[train_x,id]=unique(round(train_x,3),'rows');
train_y = csvread('train_y.csv');
train_y = train_y(id);

test_x = csvread('test_x.csv');
model = dace_train(train_x, train_y);
test_y = dace_predict(test_x, model);
colormap(jet)
scatter3(test_x(:,1), test_x(:,2),test_y(:), 100, test_y(:),'filled');
figure(2);
scatter3(train_x(:,1), train_x(:,2),  train_y(:), 40, train_y(:),'filled');
