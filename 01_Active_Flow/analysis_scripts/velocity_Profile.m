% Bulid histogram of filament minus-end velocities. Input has filament
% labels that are used to construct the sequential time path for each
% filament - much like tracking. The global array holds all filament
% data together for all time, which is then sectioned into the global_cell
% structure , where each cell element is a single filament.

tic

N_dim = 3;

min_steps = 5;  % only consider filaments that live longer then this 
delta_t   = 4;  % time increments between hist accumulations (defines distance travelled)      

% Manually construct a histogram
N_bins       = 30;
Hist_start   = 0e-6;
Hist_end     = 15.0e-6;
Hist_inc     = (Hist_end - Hist_start)/N_bins;
hist_x_bounds = Hist_start : Hist_inc : Hist_end; 
hist_x_vec    = Hist_start+Hist_inc : Hist_inc : Hist_end;
hist_x_plot   = hist_x_vec - ones(1,length(hist_x_vec))*Hist_inc/2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% use a cell to accumuate the velocity steps per bin
hist_cell = cell(length(hist_x_vec),1);

traj_array = dlmread('../trajectory_id.txt','\t'); size(traj_array)

t      = traj_array(1:end-3,1); 
N_vec  = traj_array(1:end-3,2); 
counts = traj_array(1:end-3,3); % total number of filaments nucleated

dt_sim = t(2)-t(1); % time step from simulation data output

% remove unneeded columns from input array
traj_reduced = traj_array(:,4:end); 

clear traj_array;

count_final = counts(end);

global_array = zeros(count_final,length(t),N_dim + 1); % + 1 for the time array 
traj_count   = zeros(count_final,1);

for i = 1:length(t)
    for j = 1:N_vec(i)

       count_location = traj_reduced(i,(j-1)*(N_dim+1) + 1) + 1;

       global_array(count_location,i,1) = t(i);
       global_array(count_location,i,2) = traj_reduced(i,(j-1)*(N_dim+1) + 2);
       global_array(count_location,i,3) = traj_reduced(i,(j-1)*(N_dim+1) + 3);
       global_array(count_location,i,4) = traj_reduced(i,(j-1)*(N_dim+1) + 4); 
    end
end

clear traj_reduced;

fil_count = 0;
cell_switch = 0;

% Loop through each filament
for i = 1:count_final

    t_count = 0;
    T = 0; X = 0; Y = 0; Z = 0;
    
    for j = 1:length(t)        
        if(global_array(i,j,1) ~= 0)
            if(t_count == 0)
                t_count = t_count + 1;
                T = global_array(i,j,1);
                X = global_array(i,j,2);
                Y = global_array(i,j,3);
                Z = global_array(i,j,4);
            elseif(t_count > 0)
                t_count = t_count + 1;
                T = [T,global_array(i,j,1)];
                X = [X,global_array(i,j,2)];
                Y = [Y,global_array(i,j,3)];
                Z = [Z,global_array(i,j,4)];
            end
        end
    end

    if(t_count > 0)

        loc_vec = 1;
        for dt_check = 1:length(T)-1
            t_inc = round((T(dt_check + 1) - T(dt_check))*10^4)/10^4;
            if(t_inc ~= dt_sim)
                loc_vec = [loc_vec,dt_check];
            end
        end            

        fil_count = fil_count + 1;
        if(cell_switch == 0)
           fil_cell = {t_count,T,X,Y,Z};
           global_cell = {fil_cell};
           cell_switch = 1;

        elseif(cell_switch == 1)
            fil_cell = {t_count,T,X,Y,Z};
            global_cell{end + 1} = fil_cell;
        end   

    end
end

%% 

display(fil_count)

hist_accum     = zeros(1,length(hist_x_vec));
hist_accum_std = zeros(1,length(hist_x_vec));
hist_accum_ste = zeros(1,length(hist_x_vec));
n_bin_vec      = zeros(1,length(hist_x_vec));

% Accumulate values into a histogram
for i = 1:fil_count
    if(global_cell{i}{1} > delta_t) % make sure the time domain is greater then delta t

        loop_lim = floor(global_cell{i}{1}/delta_t); % round down the time steps

        for j = 1:loop_lim

            fil_vel = (global_cell{i}{4}((j-1)*delta_t + delta_t) - global_cell{i}{4}(1 + (j-1)*delta_t))/(global_cell{i}{2}(delta_t) - global_cell{i}{2}(1));
            for bin = 1:N_bins
                mid_point = global_cell{i}{4}(1 + (j-1)*delta_t) + (global_cell{i}{4}((j-1)*delta_t + delta_t) - global_cell{i}{4}(1 + (j-1)*delta_t))/2.0;
                if(mid_point > hist_x_bounds(bin) && mid_point < hist_x_bounds(bin+1))
                    hist_cell{bin} = [hist_cell{bin},fil_vel];
                    n_bin_vec(bin)  = n_bin_vec(bin) + 1;
                    break;
                end
            end
        end
    end
end

for bin = 1:N_bins
    if(length(hist_cell{bin}) >= 5)
        hist_accum(bin) = mean(nonzeros(hist_cell{bin}));
    end
end


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plotter = figure;

width   = 12.0;
hight = 9;
set(gcf,'units',' centimeters', 'PaperUnits', 'centimeters','Position',[0,0,width,hight],...
                                       'paperPosition',[0,0,width,hight],'papersize',[width,hight])       
%%%%%%%%% Velocity Spectrum %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot_width  = 8;
plot_height = plot_width/1.3;
horz        = 2.8;
vert        = 1.8; 
axes('Parent',plotter,'units',' centimeters','Position', ...   
        [horz,vert,plot_width,plot_height],'fontweight','normal',...
        'fontname', 'times new roman'); color_map = color_map_array_red_blue;    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hold on; box on

% set(gca,'fontsize',16,'fontname', 'times new roman')
set(gca,'fontsize',16)

plot(hist_x_vec*1e6,60*1e6*hist_accum,'.-','linewidth',1.0,'markersize',7)

axis([0 15 0 20])

xlabel('Channel position (\mum)','FontSize',16)
ylabel('Velocity (\mum/min)','FontSize',16)

print('-djpeg','-r1000','velocity_Profile')

clear hist_cell global_cell global_array fil_cell;

toc
