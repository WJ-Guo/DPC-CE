				clear all;  
                clc;
                close all;
				load fourlines.txt;
				tic;
				data = fourlines(:,1:(end-1));

				%% centralize and scale the data
				data = data - repmat(mean(data),size(data,1),1);
				data = data/max(max(abs(data)));
			
				true_label = fourlines(:,end);
                
				halo = true_label;
			    NClusters = max(halo);
				% twofourlines = [[two_fourlines(xx1,[1,2]),yy1];[two_fourlines(xx2,[1,2]),-yy2]];
				    figure
					hold on
				    cmap=colormap;
					for i=1:NClusters
					   ic=int8((i*64.)/(NClusters*1.));
					   
					   if NClusters<=12
						   switch i
							case 1
								  plot(data(halo==i,1),data(halo==i,2),'o','MarkerSize',5,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
							case 2
								  plot(data(halo==i,1),data(halo==i,2),'*','MarkerSize',5,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
							case 3
								  plot(data(halo==i,1),data(halo==i,2),'x','MarkerSize',5,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
							case 4
								  plot(data(halo==i,1),data(halo==i,2),'+','MarkerSize',5,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
							case 5
								  plot(data(halo==i,1),data(halo==i,2),'s','MarkerSize',5,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
							case 6
								  plot(data(halo==i,1),data(halo==i,2),'d','MarkerSize',5,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
							case 7
								  plot(data(halo==i,1),data(halo==i,2),'v','MarkerSize',5,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
							case 8
								  plot(data(halo==i,1),data(halo==i,2),'^','MarkerSize',5,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
							case 9
								  plot(data(halo==i,1),data(halo==i,2),'<','MarkerSize',5,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
							case 10
								  plot(data(halo==i,1),data(halo==i,2),'>','MarkerSize',5,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
							case 11
								  plot(data(halo==i,1),data(halo==i,2),'p','MarkerSize',5,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
							case 12
								  plot(data(halo==i,1),data(halo==i,2),'h','MarkerSize',5,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
						   end
						   	 
					    else
						  plot(data(halo==i,1),data(halo==i,2),'.','MarkerSize',5,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));	
                        end % end of if 					  
					    
				  end %end of for
				  title('Ground-truth','FontSize',12.0)				  					
						
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            	dim = size(data,2);
				num = size(data,1);   
			%% compute and save diatance matrix if dont exists							
				mdist = [];
				
                for i=1:num  
				  mdist_i = [];
				  for j=i+1:num 
				        mdist_i = [mdist_i;[i,j, sqrt(sum((data(i,:) - data(j,:)) .^ 2)) ]];
                  end
                  mdist = [mdist; mdist_i];
                end
            save fourlinesDPCmdist.dist   mdist -ascii  
				
           %% load diatance matrix if exists				
			load fourlinesDPCmdist.dist;	
            xx = fourlinesDPCmdist;

            ND=max(xx(:,2));
			NL=max(xx(:,1));
			if (NL>ND)
			  ND=NL;  %% num of points
			end
			 
			N=size(xx,1); %% ND*(ND-1)/2,the number of pair distances
			 	
			
           %% matrix of distance?ND*ND
			for i=1:ND
			  for j=1:ND
				dist(i,j)=0;
			  end
			end
			%% matrix of distance?ND*ND, symmetry
			for i=1:N
			  ii=xx(i,1);
			  jj=xx(i,2);
			  dist(ii,jj)=xx(i,3);
			  dist(jj,ii)=xx(i,3);
			end

		%% compute dc: 2% dist
			 
			percent=2.0;
			fprintf('average percentage of neighbours (hard coded): %5.6f\n', percent);
			 
			position=round(N*percent/100); % the position of dc distance
			[sda,id]=sort(xx(:,3)); %% sort the distance
			dc=sda(position);
            
			fprintf('Computing Rho with gaussian kernel of radius: %12.6f\n', dc);
            		    

       %% initialize?rho: density of points 
			for i=1:ND
			  rho(i)=0.;
			end
			 
		% Gaussian kernel within dc and cut_off
			for i=1:ND-1
			  for j=i+1:ND
			     if (dist(i,j)<dc)
				   rho(i)=rho(i)+exp(-(dist(i,j)/dc)*(dist(i,j)/dc));
				   rho(j)=rho(j)+exp(-(dist(i,j)/dc)*(dist(i,j)/dc));
				 end  
			  end
			end
			
        %% find the max distance?
			maxd=max(max(dist));
			 
		%% rank rho by descend?
			[rho_sorted,ordrho]=sort(rho,'descend');
            rho_or = rho_sorted;
            ordrho_or = ordrho;
          
		%% deal with point with max rho?
			delta(ordrho(1))=-1.;
			nneigh(ordrho(1))=0; 

		%% compute the delta(relative distance), find nneigh for points
			for ii=2:ND
			   delta(ordrho(ii))=maxd;
			   for jj=1:ii-1
				 if(dist(ordrho(ii),ordrho(jj))<delta(ordrho(ii)))
					delta(ordrho(ii))=dist(ordrho(ii),ordrho(jj));
					nneigh(ordrho(ii))=ordrho(jj);
				   % nneigh record the point with higher rho and most close to ordrho(ii)
				 end
			   end
			end
			 
		%% give max rho point max delta 
			delta(ordrho(1))=max(delta(:));
            
        %% CES  £¨¾àÀë²¹³¥£©
            count= [0 0 0]; % record punishment
            choose_num = 20;
            dist_ther1_ratio = 0.25;
            punish_ratio = 0.3;
            
            [delta_sorted,orddelta]=sort(delta,'descend');
            choose_point = orddelta(1:choose_num);
            rho_choose_point = rho(choose_point);
            [~,ordrho_choose]=sort(rho_choose_point,'descend');
            for i = 1:size(choose_point,2)
                pair_all(i) = choose_point(ordrho_choose(i));
            end
            condition = [];
            for i = 2:size(choose_point,2)   
                pair_use = [];
                for j = 1 : i
                    if j == 1
                        pair = [pair_all(i), nneigh(pair_all(i))];
                        pair_use = [pair_use,nneigh(pair_all(i))];
                    else
                        if   ismember(pair_all(j-1),pair_use)
                               continue
                        else
                            pair = [pair_all(i), pair_all(j-1)];
                            pair_use = [pair_use, pair_all(j-1)];
                            
                        end
                    end
                    
                    dist_ther1 = dist(pair(1),pair(2));
                    dist_ther = dist_ther1 * dist_ther1_ratio;
                    [turn,nei_all,find_time,count,nei_dist] = judge_nei(pair,dist,count,dc, dist_ther);        
                           
                    if turn == 2
                        punish_time = find_time-5;
                        dist_new = nei_dist  + dist_ther * (1+ punish_time * punish_ratio);
                        dist(pair(1),pair(2)) = dist_new;
                        dist(pair(2),pair(1)) = dist(pair(1),pair(2));
                        if dist_new < delta(pair_all(i))
                            delta(pair_all(i)) = dist_new;
                            nneigh(pair_all(i)) = pair(2);
                        end
                    end

                     if turn == 3 && j == 1
                         % max value between 2 and 1's all neighbor
                            max_dist_nei = max(dist(pair(2),nei_all));
                            dist(pair(1),pair(2)) =  max_dist_nei * 1.1;
                            dist(pair(2),pair(1)) = dist(pair(1),pair(2));
                            delta(pair_all(i)) = dist(pair(1),pair(2));
                    end                    
                   
%                     if pair == [66, 47]
%                         turn
%                         dist_ther
%                         nei_all
%                     end       

                    condition1 = [pair,turn,find_time,nei_dist,dist_ther1,dist(pair(1),pair(2))];
                    condition = [condition; condition1]; 
                end
                
             end                          

maxd=max(max(dist));            
delta(ordrho(1))=max(delta(:));            

                       
			toc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%			
		%% decision graph 
			 
			disp('Generated file:DECISION GRAPH')
			disp('column 1:Density')
			disp('column 2:Delta')
			 
			fid = fopen('DECISION_GRAPH', 'w');
			for i=1:ND
			   fprintf(fid, '%6.2f %6.2f\n', rho(i),delta(i));
			end
			 
		% Select a rectangle enclosing cluster centers
			disp('Select a rectangle enclosing cluster centers')

			scrsz = get(0,'ScreenSize');			 
			 
			fig = figure;
			tt=plot(rho(:),delta(:),'o','MarkerSize',5,'MarkerFaceColor','k','MarkerEdgeColor','k');
			title ('Decision Graph of DPC-CE','FontSize',12.0)
			box off;
			xlabel ('\rho','FontSize',14.0)
			ylabel ('\delta','FontSize',14.0)
			 
		
			rect = getrect(fig);
			rhomin=rect(1);
			deltamin=rect(2); 
			 
		% initialize number of cluster
			NCLUST=0;

%%%%%%%%%%%%%%%%%%%%%%% cluster center%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
			
		%% cl(i)=j: point(i) belong to cluster j, initialize cl = -1;
			for i=1:ND
			  cl(i)=-1;
            end
			 
		   % find cluster center?
			for i=1:ND
			  if ( (rho(i)>rhomin) && (delta(i)>deltamin))
				 NCLUST=NCLUST+1;
				 cl(i)=NCLUST; %% cl(i) is center
				 icl(NCLUST)=i;%% icl record index of cluster center
			  end
			end
			 
			fprintf('NUMBER OF CLUSTERS: %i \n', NCLUST);
			a = axis;
            
			% limits of current figure
			xmin=a(1);
			xmax=a(2);
			ymin=a(3);
			ymax=a(4);
			  
			% makes grid 
			options.gridx = 50;
			options.gridy = 50;
			[X,Y] = meshgrid(xmin:(xmax-xmin)/options.gridx:xmax,...
							 ymin:(ymax-ymin)/options.gridy:ymax);

			% make testing patterns covering whole grid
			tst_data=[reshape(X',1,prod(size(X)));reshape(Y',1,prod(size(Y)))];
			dec_fun= tst_data(1,:).*tst_data(2,:);
			% reshape dec_fun
            Z = reshape(dec_fun,size(X,1),size(X,2))';
			% smooth shading
			hold on
            %contour(X,Y,Z,1,'k');
			
			% draw cluster centers with different color
            cmap=colormap;
			for i=1:NCLUST
			   ic=int8((i*64.)/(NCLUST*1.));
			   hold on
			   plot(rho(icl(i)),delta(icl(i)),'o','MarkerSize',8,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
			   hold on
			   cmap(ic,:)
               contour(X,Y,Z,[rho(icl(i))*delta(icl(i)) rho(icl(i))*delta(icl(i)) rho(icl(i))*delta(icl(i))],'linecolor',cmap(ic,:));
			end
			
		
	
			  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%		
		tic;
		    disp('Performing assignation')
					 
					%% clusting non_center points, traverse by rho?
					for i=1:ND
					  if (cl(ordrho(i))==-1)
						cl(ordrho(i))=cl(nneigh(ordrho(i)));
					  end
                    end		
		
		 %% deal with halo
					for i=1:ND
					  halo(i)=cl(i);
					end
															
					if (NCLUST>1)
					 
							  % initialize bord_rho for every cluster = 0
							  for i=1:NCLUST
								bord_rho(i)=0.;
							  end
							 					 
							  % calculate average bird_rho for each cluster
							  for i=1:ND-1
								for j=i+1:ND
								  
							% two close points but in different cluster
								  if ((cl(i)~=cl(j)) && (dist(i,j)<=dc))
									rho_aver=(rho(i)+rho(j))/2.; %% the average rho is theshold  ?
									
									if (rho_aver>bord_rho(cl(i)))
									  bord_rho(cl(i))=rho_aver;
									end
									
									if (rho_aver>bord_rho(cl(j)))
									  bord_rho(cl(j))=rho_aver;
									end
								  end
								end
							  end

						%% for outline(noise points)
% 							  for i=1:ND
% 								if (rho(i)<bord_rho(cl(i)))
% 								  halo(i)=0;
% 								end
% 							  end
					 
					end
					
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%halo(i)=0,outlier, halo(i)=1,normal%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
					
					 %% for every cluster
							for i=1:NCLUST
								  nc=0; %% the num of all points in cluster
								  nh=0; %% the num of core points in cluster 
								  for j=1:ND
									if (cl(j)==i)
									  nc=nc+1;
									end
									if (halo(j)==i) % non_outlier
									  nh=nh+1;
									end
								  end
								 
								  fprintf('CLUSTER: %i CENTER: %i ELEMENTS: %i CORE: %i HALO: %i \n', i,icl(i),nc,nh,nc-nh);
								 
							end
												
					 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%end of DPC %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	

			    if (dim==2)
				    figure
					hold on
				    cmap=colormap;
					for i=1:NCLUST
					   ic=int8((i*64.)/(NCLUST*1.));
					   if NCLUST<=12
						   switch i
							case 1
								  plot(data(icl(i),1),data(icl(i),2),'h','MarkerSize',8,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
								  plot(data(halo==i,1),data(halo==i,2),'o','MarkerSize',5,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
							case 2
								  plot(data(icl(i),1),data(icl(i),2),'h','MarkerSize',8,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
								  plot(data(halo==i,1),data(halo==i,2),'*','MarkerSize',5,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
							case 3
								  plot(data(icl(i),1),data(icl(i),2),'o','MarkerSize',8,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
								  plot(data(halo==i,1),data(halo==i,2),'x','MarkerSize',5,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
							case 4
								  plot(data(icl(i),1),data(icl(i),2),'o','MarkerSize',8,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
								  plot(data(halo==i,1),data(halo==i,2),'+','MarkerSize',5,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
							case 5
								  plot(data(icl(i),1),data(icl(i),2),'o','MarkerSize',8,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
								  plot(data(halo==i,1),data(halo==i,2),'s','MarkerSize',5,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
							case 6
								  plot(data(icl(i),1),data(icl(i),2),'o','MarkerSize',8,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
								  plot(data(halo==i,1),data(halo==i,2),'d','MarkerSize',5,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
							case 7
								  plot(data(icl(i),1),data(icl(i),2),'o','MarkerSize',8,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
								  plot(data(halo==i,1),data(halo==i,2),'v','MarkerSize',5,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
							case 8
								  plot(data(icl(i),1),data(icl(i),2),'o','MarkerSize',8,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
								  plot(data(halo==i,1),data(halo==i,2),'^','MarkerSize',5,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
							case 9
								  plot(data(icl(i),1),data(icl(i),2),'o','MarkerSize',8,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
								  plot(data(halo==i,1),data(halo==i,2),'<','MarkerSize',5,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
							case 10
								  plot(data(icl(i),1),data(icl(i),2),'o','MarkerSize',8,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
								  plot(data(halo==i,1),data(halo==i,2),'>','MarkerSize',5,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
							case 11
								  plot(data(icl(i),1),data(icl(i),2),'o','MarkerSize',8,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
								  plot(data(halo==i,1),data(halo==i,2),'p','MarkerSize',5,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
							case 12
								  plot(data(icl(i),1),data(icl(i),2),'o','MarkerSize',8,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
								  plot(data(halo==i,1),data(halo==i,2),'h','MarkerSize',5,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
							otherwise
								  plot(data(icl(i),1),data(icl(i),2),'o','MarkerSize',8,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
								  plot(data(halo==i,1),data(halo==i,2),'h','MarkerSize',5,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
						   end
						  
					  else
						   plot(data(icl(i),1),data(icl(i),2),'o','MarkerSize',8,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
						   plot(data(halo==i,1),data(halo==i,2),'.','MarkerSize',5,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));	
                      end 					  
					    
				  end
				  plot(data(halo==0,1),data(halo==0,2),'.','MarkerSize',8,'MarkerFaceColor','k','MarkerEdgeColor','k');	
				 
			end
			title('DPC-CE','FontSize',12.0)
toc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  FMI,ARI,NMI  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
				
                DPCFMI = FMI(true_label(:),cl(:));
				fprintf('FMI value of DPC on fourlines dataset %i \n', DPCFMI);
				[cluster_acc,randindex,DPCARI] = ARI(true_label(:),cl(:));
				fprintf('ARI value of DPC on fourlines dataset %i \n', DPCARI);
				DPCNMI = NMI(true_label(:),cl(:));
				fprintf('NMI value of DPC on fourlines dataset %i \n', DPCNMI);
								 
				