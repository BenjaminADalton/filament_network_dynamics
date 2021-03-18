% plot time and length data produced using the filament program

tic

N_dim = 3;

dir_plot = 1;

file_string = 'AMR_speckle';

% Set box dimensions from simulation data
L_x = 1.0e-6;
L_y = 30.0e-6;
L_z = 1e-6;

N_PIN       = 16;
N_AMR_Probe = 1;

% Set figure aspect ratio
Y_low  = 0;
Y_high = L_y;

x_asp = L_x*1e6;
y_asp = (Y_high - Y_low)*1e6;
z_asp = L_x*1e6;

% read in simulation trajectory (see ReadMe for formatting description)
traj_array = dlmread('../trajectory.txt','\t');

t = traj_array(:,1);

display(['Size MT trajectory - ',num2str(size(traj_array))])

N_vec  = traj_array(:,2); % number of filaments

set(0,'DefaultFigureVisible','on')
    
figure, set(gcf, 'Color','white','unit','centimeters','position',[0,0,48,8],'papersize',[48,8],'paperunits','centimeters')
axes('Parent',gcf,'fontsize',200,'unit','centimeters','Position',[2,0.5,44,8])

% create AVI objec
vid_object = VideoWriter(file_string,'MPEG-4');
vid_object.Quality = 100;
vid_object.FrameRate = 25;
open(vid_object);  

 for i = 1:1:length(t)-1

disp(['printing: ',num2str(i)]);

%%%%%%%%%%%%%%%% read data into arrays for plotting %%%%%%%%%%%%%%%%
L_vec     = zeros(N_vec(i),1);
r0        = zeros(N_dim*N_vec(i),1);
u0        = zeros(N_dim*N_vec(i),1);

% read filament data
k=0;
for n = 1:N_vec(i)
    L_vec(n) = traj_array(i,(3 + (n-1)*(2*N_dim + 1)));
    for j = 1:N_dim
    k = k+1;
    r0(k) = traj_array(i,(3+(n-1)*(2*N_dim+1) + j));
    u0(k) = traj_array(i,(3 + N_dim +(n-1)*(2*N_dim+1) + j));
    end
end    
    
plot3(0,0,0,'k','linewidth',2)
    
hold on; box on

pbaspect([x_asp y_asp z_asp])

% plot filaments
for n = 1:N_PIN+N_AMR_Probe
    if(u0((n-1)*N_dim+2)<=0)
    plot3([r0((n-1)*N_dim+1)-0.5*L_vec(n)*u0((n-1)*N_dim+1),r0((n-1)*N_dim+1)+0.5*L_vec(n)*u0((n-1)*N_dim+1)],...
	[r0((n-1)*N_dim+2)-0.5*L_vec(n)*u0((n-1)*N_dim+2),r0((n-1)*N_dim+2)+0.5*L_vec(n)*u0((n-1)*N_dim+2)],...
    [r0((n-1)*N_dim+3)-0.5*L_vec(n)*u0((n-1)*N_dim+3),r0((n-1)*N_dim+3)+0.5*L_vec(n)*u0((n-1)*N_dim+3)],'k','linewidth',1)
    elseif(u0((n-1)*N_dim+2)>0)
    plot3([r0((n-1)*N_dim+1)-0.5*L_vec(n)*u0((n-1)*N_dim+1),r0((n-1)*N_dim+1)+0.5*L_vec(n)*u0((n-1)*N_dim+1)],...
	[r0((n-1)*N_dim+2)-0.5*L_vec(n)*u0((n-1)*N_dim+2),r0((n-1)*N_dim+2)+0.5*L_vec(n)*u0((n-1)*N_dim+2)],...
    [r0((n-1)*N_dim+3)-0.5*L_vec(n)*u0((n-1)*N_dim+3),r0((n-1)*N_dim+3)+0.5*L_vec(n)*u0((n-1)*N_dim+3)],'color',[0,0,0]+0.5,'linewidth',0.9)            
    end
end


for n = N_PIN+N_AMR_Probe+1:N_vec(i)
    plot3(r0((n-1)*N_dim+1)-0.5*L_vec(n)*u0((n-1)*N_dim+1),...
	r0((n-1)*N_dim+2)-0.5*L_vec(n)*u0((n-1)*N_dim+2),...
    r0((n-1)*N_dim+3)-0.5*L_vec(n)*u0((n-1)*N_dim+3),'.k','markersize',8)
end


scale = 0.1;

axis([-scale*L_x (1+scale)*L_x Y_low Y_high -scale*L_z (1+scale)*L_z])

view(90 ,90)

xlabel('X')
ylabel('Y')
zlabel('Z')

set(gca,'fontsize',12,'fontname', 'times new roman')
writeVideo(vid_object, getframe(gcf));    

pause(0.01)
hold off

end

close(gcf)
close(vid_object);  

toc