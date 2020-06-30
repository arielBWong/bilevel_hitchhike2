%%

f = [0.38061, 0.05548];
y = [0.34402, 0.60999];
sig = [0.49642, 1.19668];
img = (f-y)./sig;
gausscdf(img)
eim = (f-y).*gausscdf(img) + sig.*gausspdf(img)
r = [1.1, 1.1]
a=r(ones(6, 1), :)
prod(r-f+eim) 
prod(r-f)




function y = gausscdf(x)
y = 0.5 * (1 +  erf(x/sqrt(2)));
end

function y = gausspdf(x)
y = 1/sqrt(2*pi)*exp(-x.^2/2);
end