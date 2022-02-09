function [turn,nei_all,find_time,count,nei_dist] = judge_nei(pair,dist,count,dc,dist_ther)

    turn = 0;
    find_time = 0;
    nei_all = [pair(1)];
    nei_new = nei_all;
    nei_dist = 0;
    while turn == 0

        find_time = find_time+1;
        nei_temp = [];
         if  dist(pair(1),pair(2)) < dc*2
            turn = 1; 
            count(1) = count(1)+1;
         end
        

        for k = 1: size(nei_new,2)
            nei_add = find(dist(nei_new(k),:) < dist_ther);
            nei_dist_max = mean(dist(nei_new(k),nei_add));
            if nei_dist_max > nei_dist
                nei_dist = nei_dist_max;
            end
            nei_temp = [nei_temp nei_add];
        end

        nei_new = setdiff(nei_temp, nei_all);
       
        if ismember(pair(2), nei_all) %% 判断为连通
            turn = 2; 
            count(2) = count(2)+1;   
            
        elseif length(nei_new)<1  %% 判断不连通
            turn = 3;
            count(3) = count(3)+1;
        end
        nei_all = [nei_all,nei_new];
        
    end  