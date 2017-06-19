function [] = create_avi(realdata,gendata,dir_ref)

mapfilename = [dir_ref,'/create_avi-map-vid.avi'];
linksfilename = [dir_ref,'/create_avi-links-vid.avi'];

contact_time = 20;
close all
scale = 20;

all_people1 = [realdata(:,2)' realdata(:,3)'];
num_times1 = size(unique(realdata(:,1)),1);
data_length1 = size(realdata(:,1),1);
num_people1 = max(all_people1);
coords1 = zeros(num_people1,2);
theta1 = 2*pi/num_people1;

all_people2 = [gendata(:,2)' gendata(:,3)'];
num_times2 = size(unique(gendata(:,1)),1);
data_length2 = size(gendata(:,1),1);
num_people2 = max(all_people2);
coords2 = zeros(num_people2,2);
theta2 = 2*pi/num_people2;

num_times = min(num_times1,num_times2);

parfor n=1:num_people1
    coords1(n,:) = scale*[sin(n*theta1) cos(n*theta1)];
end

parfor n=1:num_people2
    coords2(n,:) = scale*[sin(n*theta2) cos(n*theta2)];
end

step_nf1 = zeros(num_people1,1);
line_freq1 = zeros(num_people1,num_people1);
step_nf2 = zeros(num_people2,1);
line_freq2 = zeros(num_people2,num_people2);

links(num_times) = struct('cdata',[],'colormap',[]);
map(num_times) = struct('cdata',[],'colormap',[]);

for m=1:num_times
    thisadj1 = zeros(num_people1);
    thisadj2 = zeros(num_people2);
    current_time = (m-1)*contact_time;
    for i=1:data_length1
        test_time = realdata(i,1);
        if test_time==current_time
            person1 = realdata(i,2);
            person2 = realdata(i,3);
            thisadj1(person1,person2) = 1;
            thisadj1(person2,person1) = 1;
        end
    end
    for i=1:data_length2
        test_time = gendata(i,1);
        if test_time==current_time
            person1 = gendata(i,2);
            person2 = gendata(i,3);
            thisadj2(person1,person2) = 1;
            thisadj2(person2,person1) = 1;
        end
    end
    step_nf1 = step_nf1+sum(thisadj1,2);
    step_nf2 = step_nf2+sum(thisadj2,2);
    
    %==Create Adj Frame==%
    map_fig = figure();
    subplot(1,2,1)
    gplot(thisadj1,coords1,'-*');
    str = sprintf('Time: %d', current_time);
    text(0,-1.2*scale,str);
    axis([-1.5 1.5 -1.5 1.5]*scale);
    axis off;
    set(gcf,'color','w');
    subplot(1,2,2)
    gplot(thisadj2,coords2,'-*');
    str = sprintf('Time: %d', current_time);
    text(0,-1.2*scale,str);
    axis([-1.5 1.5 -1.5 1.5]*scale);
    axis off;
    set(gcf,'color','w');
    map(m) = getframe(map_fig);
    close(map_fig);
    %==End Frame Creation==%
    
    this_rel_node_freq1 = step_nf1/max(step_nf1);
    this_rel_node_freq2 = step_nf2/max(step_nf2);
	line_freq1 = line_freq1+thisadj1;
    line_freq2 = line_freq2+thisadj2;
    thisRLF1 = line_freq1/(max(max(line_freq1)));
    thisRLF2 = line_freq2/(max(max(line_freq2)));
    
    %==Create Activity Frame==%
	links_fig = figure();
    subplot(1,2,1)
    link_size1 = this_rel_node_freq1*50;
    tempadj1 = logical(thisRLF1);
    [row1,col1] = find(tempadj1);
    tempcoords1 = zeros(num_people1,2);
    hold on
    for i=1:length(row1)
        thisrow = row1(i);
        thiscol = col1(i);
        if thisrow >= thiscol
            line_col = (1-thisRLF1(thisrow,thiscol))*[1 1 1];
            x1 = coords1(thisrow,1);
            y1 = coords1(thisrow,2);
            x2 = coords1(thiscol,1);
            y2 = coords1(thiscol,2);
            line([x1 x2],[y1 y2],'Color',line_col)
        end
        tempcoords1(i,:) = coords1(thisrow,:);
    end
    tempcoords1 = tempcoords1(any(tempcoords1,2),:);
    tempcoords1 = unique(tempcoords1,'rows','stable');
    link_size1(link_size1==0) = [];
    scatter(tempcoords1(:,1),tempcoords1(:,2),link_size1,'filled')
    hold off
    str = sprintf('Time: %d', current_time);
    text(0,-1.2*scale,str);
    axis([-1.5 1.5 -1.5 1.5]*scale);
    axis off;
    set(gcf,'color','w');
    subplot(1,2,2)
    link_size2 = this_rel_node_freq2*50;
    tempadj2 = logical(thisRLF2);
    [row2,col2] = find(tempadj2);
    tempcoords2 = zeros(num_people2,2);
    hold on
    for i=1:length(row2)
        thisrow = row2(i);
        thiscol = col2(i);
        if thisrow >= thiscol
            line_col = (1-thisRLF2(thisrow,thiscol))*[1 1 1];
            x1 = coords2(thisrow,1);
            y1 = coords2(thisrow,2);
            x2 = coords2(thiscol,1);
            y2 = coords2(thiscol,2);
            line([x1 x2],[y1 y2],'Color',line_col)
        end
        tempcoords2(i,:) = coords2(thisrow,:);
    end
    tempcoords2 = tempcoords2(any(tempcoords2,2),:);
    tempcoords2 = unique(tempcoords2,'rows','stable');
    link_size2(link_size2==0) = [];
    scatter(tempcoords2(:,1),tempcoords2(:,2),link_size2,'filled')
    hold off
    str = sprintf('Time: %d', current_time);
    text(0,-1.2*scale,str);
    axis([-1.5 1.5 -1.5 1.5]*scale);
    axis off;
    set(gcf,'color','w');
    links(m) = getframe(links_fig);
    close(links_fig);
end

v = VideoWriter(mapfilename);
open(v)
writeVideo(v,map)
close(v)

v = VideoWriter(linksfilename);
open(v)
writeVideo(v,links)
close(v)

end