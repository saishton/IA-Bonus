function [] = DDVid(DinT,dir_ref)

filename = [dir_ref,'/DegreeDistribution.avi'];
contact_time = 20;
num_times = size(DinT,2);
frame(num_times) = struct('cdata',[],'colormap',[]);

for m=1:num_times
    current_time = (m-1)*contact_time;
    thisframe = figure();
    ecdf(DinT(:,m));
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