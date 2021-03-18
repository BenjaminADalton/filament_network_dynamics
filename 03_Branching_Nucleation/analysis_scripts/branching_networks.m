% 14-03-2021 - Benjamin A. Dalton
% Generate a movie to visualize the interaction of two actively driven,
% branching filament networks, interacting via active motors in
% confinement. The movie is a 2D projection of a 3D system. To view
% full 3D system, comment-out view(90 ,90), and consider increasing
% the figure height and paper size.

tic

N_dim = 3;

dir_plot = 1;

file_string = 'branching_networks2';

% Set box dimensions from simulation data
L_x = 1.0e-6;
L_y = 6.0e-6;
L_z = 1e-6;

% Set figure aspect ratio
Y_low  = -L_y;
Y_high = L_y;

x_asp = L_x*1e6;
y_asp = (Y_high - Y_low)*1e6;
z_asp = L_x*1e6;

% read in simulation trajectory (see ReadMe for formatting description)
traj_array = dlmread('../trajectory.txt','\t');
cl_array   = dlmread('../trajectory_cl.txt','\t');
br_array   = dlmread('../trajectory_branch.txt','\t');

t = traj_array(:,1);

display(['Size MT trajectory - ',num2str(size(traj_array))])
display(['Size CL trajectory - ',num2str(size(cl_array))])
display(['Size Br trajectory - ',num2str(size(br_array))])

N_vec  = traj_array(:,2); % number of filaments
N_cl_st1 = cl_array(:,2); % number of state-1 CLs
N_cl_st2 = cl_array(:,3); % number of state-1 CLs
N_br     = br_array(:,3); % number of branch connections

set(0,'DefaultFigureVisible','on')
    
figure, set(gcf, 'Color','white','unit','centimeters','position',[0,0,48,8],'papersize',[48,8],'paperunits','centimeters')
axes('Parent',gcf,'fontsize',200,'unit','centimeters','Position',[1.5,1.5,44,5])

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
r0_cl_st1 = zeros(N_dim*N_cl_st1(i),1);
r0_cl_st2 = zeros(2*N_dim*N_cl_st2(i),1);
cl_st2_id = zeros(N_cl_st2(i),1);
r0_br     = zeros(2*N_dim*N_br(i),1);

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

% read state-1 and state-2 CL data
k=0; l=0; id = 0; cl_st_inc = 4;

for n = 1:N_cl_st1(i) + N_cl_st2(i)
    if(cl_array(i,cl_st_inc)==1)
        for j = 1:N_dim
            k=k+1;
            r0_cl_st1(k) = cl_array(i,(cl_st_inc + 1 + j));
        end
        cl_st_inc = cl_st_inc + 1*N_dim + 2;
    end
    if(cl_array(i,cl_st_inc)==2)
    id = id + 1;    
    cl_st2_id(id) = cl_array(i,(cl_st_inc + 1));
        for j = 1:2*N_dim
            l=l+1;
            r0_cl_st2(l) = cl_array(i,(cl_st_inc + 1 + j));
        end
        cl_st_inc = cl_st_inc + 2*N_dim + 2;
    end    
end

% read branching CL data
k=0; l=0; br_st_inc = 4;

for n = 1:(N_br(i))
    if(br_array(i,br_st_inc) == 2)
        for j = 1:2*N_dim
            l = l + 1;
            r0_br(l) = br_array(i,(br_st_inc + j));
        end
        br_st_inc = br_st_inc + 2*N_dim+1;
    end
end    
    
plot3(0,0,0,'k','linewidth',2)
    
hold on; box on

pbaspect([x_asp y_asp z_asp])

% plot filaments
for n = 1:N_vec(i)
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

% plot branching connections
for n = 1:2:N_br(i)
    if(abs(r0_br(2*(n-1)*N_dim+1) - r0_br(2*(n-1)*N_dim+4)) < L_x/2 && ...
            abs(r0_br(2*(n-1)*N_dim+2) - r0_br(2*(n-1)*N_dim+5)) < L_y/2 && ...
            abs(r0_br(2*(n-1)*N_dim+3) - r0_br(2*(n-1)*N_dim+6)) < L_z/2)
        
        plot3([r0_br(2*(n-1)*N_dim+1),r0_br(2*(n-1)*N_dim+4)],... 
              [r0_br(2*(n-1)*N_dim+2),r0_br(2*(n-1)*N_dim+5)],...
              [r0_br(2*(n-1)*N_dim+3),r0_br(2*(n-1)*N_dim+6)],'b','linewidth',0.9)
    end
end

% plot state-1 CLs
% for n = 1:N_cl_st1(i)
%     plot3(r0_cl_st1((n-1)*N_dim+1),r0_cl_st1((n-1)*N_dim+2),r0_cl_st1((n-1)*N_dim+3),'.r','markersize',6); 
% end


for n = 1:N_cl_st2(i)
    if(cl_st2_id(n) == 1)
        if(abs(r0_cl_st2(2*(n-1)*N_dim+1) - r0_cl_st2(2*(n-1)*N_dim+4)) < L_x/2 && ...
                abs(r0_cl_st2(2*(n-1)*N_dim+2) - r0_cl_st2(2*(n-1)*N_dim+5)) < L_y/2 && ...
                abs(r0_cl_st2(2*(n-1)*N_dim+3) - r0_cl_st2(2*(n-1)*N_dim+6)) < L_z/2)
            plot3([r0_cl_st2(2*(n-1)*N_dim+1),r0_cl_st2(2*(n-1)*N_dim+4)],...
                  [r0_cl_st2(2*(n-1)*N_dim+2),r0_cl_st2(2*(n-1)*N_dim+5)],... 
                  [r0_cl_st2(2*(n-1)*N_dim+3),r0_cl_st2(2*(n-1)*N_dim+6)],'r','linewidth',0.9)
        end
    elseif(cl_st2_id(n) == 2)
        if(abs(r0_cl_st2(2*(n-1)*N_dim+1) - r0_cl_st2(2*(n-1)*N_dim+4)) < L_x/2 && ...
                abs(r0_cl_st2(2*(n-1)*N_dim+2) - r0_cl_st2(2*(n-1)*N_dim+5)) < L_y/2 && ...
                abs(r0_cl_st2(2*(n-1)*N_dim+3) - r0_cl_st2(2*(n-1)*N_dim+6)) < L_z/2)
            plot3([r0_cl_st2(2*(n-1)*N_dim+1),r0_cl_st2(2*(n-1)*N_dim+4)],...
                  [r0_cl_st2(2*(n-1)*N_dim+2),r0_cl_st2(2*(n-1)*N_dim+5)],... 
                  [r0_cl_st2(2*(n-1)*N_dim+3),r0_cl_st2(2*(n-1)*N_dim+6)],'m','linewidth',0.9)
        end
    end
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
