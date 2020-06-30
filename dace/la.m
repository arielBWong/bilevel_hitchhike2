%%

close all;clc;
horizen = {'hv min','hv max','hv ave','hv std',... 
                             'eim min','eim max','eim ave','eim std',...
                             'hvr min','hvr max','hvr ave','hvr std'};
problem = {'ZDT1', 'ZDT2', 'ZDT3' , 'DTLZ2','DTLZ4', 'DTLZ1'};


x = csvread('_cheat_100000hv_eim_hvr.csv');
size(x);
x(:,end) = [];

latextable(x,'Horiz',horizen, 'Vert',problem,'format','%.4f')