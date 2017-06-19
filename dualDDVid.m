function [] = dualDDVid(DinT_real,DinT_gen,dir_ref)

filename = [dir_ref,'/DegreeDistribution.avi'];
contact_time = 20;
num_times1 = size(DinT_real,2);
num_times2 = size(DinT_gen,2);
num_times = min(num_times1,num_times2);
frame(num_times) = struct('cdata',[],'colormap',[]);

for m=1:num_times
    current_time = (m-1)*contact_time;
    thisframe = figure();
    subplot(1,2,1)
    ecdf(DinT_real(:,m));
    str = sprintf('Time: %d', current_time);
    text(4,0.1,str);
    axis([0 5 0 1]);
    subplot(1,2,2)
    ecdf(DinT_gen(:,m));
    str = sprintf('Time: %d', current_time);
    text(4,0.1,str);
    axis([0 5 0 1]);
    frame(m) = getframe(thisframe);
    close(thisframe);
end

v = VideoWriter(filename);
open(v)
writeVideo(v,frame)
close(v)

end