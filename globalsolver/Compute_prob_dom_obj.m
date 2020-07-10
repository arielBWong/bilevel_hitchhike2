function [p]= Compute_prob_dom_obj(mu1,mu2,sigma1,sigma2,prob)
for k = 1:prob.nf
    if((sigma1(1,k)==0) && (sigma2(1,k)==0))
        if(mu1(1,k)<mu2(1,k))
            p(k)=3;
        end
        if(mu1(1,k)==mu2(1,k))
            p(k)=2;
        end
        if(mu1(1,k)>mu2(1,k))
            p(k)=1;
        end
    else
        p(k)=1+2*(0.5+0.5*erf((mu2(1,k)-mu1(1,k))/sqrt(2*(sigma1(1,k)^2+sigma2(1,k)^2))));
    end
end
return