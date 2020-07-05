
%%
% this snappet shows how dace internally deal with 
% repeated x value
num_xu = size(xu, 1);
      dismatrix = zeros(num_xu);
      for ii = 1:num_xu
          for jj = 1:num_xu
              if ii==jj
                  dismatrix(ii, jj) = 10;
              else
                  dismatrix(ii, jj) = sum(abs(xu(ii, :)-xu(jj, :)), 2);
              end
          end
      end
      min_d = min(dismatrix(:));
      disp(min_d);