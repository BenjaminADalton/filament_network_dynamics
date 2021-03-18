% 14-03-2021 - Benjamin A. Dalton
% Show the number of filaments and cross-linkers/motors as as
% function of time.

N_PIN = 16;
dt    = 2.5e-6;

filament_read    = dlmread('../trajectory.txt');
crosslinker_read = dlmread('../trajectory_cl.txt');



time_fil = dt*filament_read(1:end-1,1);
N_FILS   = filament_read(1:end-1,2);

time_cl = dt*crosslinker_read(1:end-1,1);
N_CL_s1 = crosslinker_read(1:end-1,2);
N_CL_s2 = crosslinker_read(1:end-1,3);


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


plot(time_fil,N_FILS,'k','linewidth',1.0)
plot(time_cl,N_CL_s1,'r','linewidth',1.0)
plot(time_cl,N_CL_s2,'b','linewidth',1.0)
%plot(time,1e6*trap_y_delta,'k','linewidth',1.0)

% axis([0 time(end) -1.5 1.5])

leg = legend('MTs','state-1 CLs','state-2 CLs','location','northwest','fontsize',8,'box','off');
leg.ItemTokenSize = [10,10];

set(gca,'fontsize',15,'fontname', 'Avenir Next')

xlabel('time (s)','FontSize',16)
ylabel('Number','FontSize',16)

print('-djpeg','-r1000','numbers_vs_time') 


