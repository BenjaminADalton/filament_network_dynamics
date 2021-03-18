% 14-03-2021 - Benjamin A. Dalton
% Script to read in pinning field and trapped filament trajectory and plot
% Trapped trajectory is shifted by its value at t=0 to show the driving 
% field and signal together. 
% If the trajectory of the trapped filament does not follow closely to the
% pinning field, then either the cross-linker concentration (N_CLS) or 
% the filament density is too low (increase k_n), i.e. the system is inelastic 

trajectory_read = dlmread('../trajectory_trap.txt');

time            = trajectory_read(1:end-1,1);
pinning_field_y = trajectory_read(1:end-1,3);
trap_y_abs      = trajectory_read(1:end-1,6);
trap_y_delta    = trap_y_abs - trap_y_abs(1);

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


plot(time,1e6*pinning_field_y,'r','linewidth',1.0)
plot(time,1e6*trap_y_delta,'k','linewidth',1.0)

axis([0 time(end) -1.5 1.5])

leg = legend('pinning field','trapped filament','location','southwest','fontsize',10,'box','off');
leg.ItemTokenSize = [10,10];

set(gca,'fontsize',15,'fontname', 'Avenir Next')

xlabel('time (s)','FontSize',16)
ylabel('channel position (\mum)','FontSize',16)

print('-djpeg','-r1000','AMR_Signal') 


