function DinT = deg_in_time(data)

number_times = (max(data(:,1))/20)+1;
number_people = max([data(:,2); data(:,3)]);
contact_time = 20;

DinT = zeros(number_people,number_times);

parfor i=1:number_times
    this_DinT = zeros(number_people,1);
    this_time = (i-1)*contact_time;
    this_idx = (data(:,1) == this_time);
    this_data = data(this_idx,:);
    this_active = [this_data(:,2); this_data(:,3)];
    if isempty(this_active)==0
        for j=1:length(this_active)
            ID = this_active(j);
            this_DinT(ID,1) = this_DinT(ID,1)+1;
        end
    end
    DinT(:,i) = this_DinT;
end
end