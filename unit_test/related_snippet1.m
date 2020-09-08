
%%
num_point = 100;
 sep = pi/(2 *(num_point-1));
            pf = [1, 0];
            
            deg = 0;
            for i = 1:num_point-1
                one = [cos(deg + i * sep), sin(deg + i * sep)];
                pf = [pf; one];
            end
            

scatter(pf(:, 1), pf(:,2)); drawnow;
