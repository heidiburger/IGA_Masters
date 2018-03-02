for n_elem=1:5
p=4;

x_start=0;
x_end=1;

knot_vec_int = linspace(x_start,x_end,n_elem+1);

if p == 1
    Xi_to_plot = linspace(x_start,x_end,p+n_elem);
elseif p == 2
    Xi_to_plot = [x_start,(knot_vec_int(2:end)-knot_vec_int(1:end-1))/2+knot_vec_int(1:end-1),x_end];
elseif p == 3
    if n_elem == 1
        Xi_to_plot = [x_start,x_start+(x_end-x_start)/3,x_end - (x_end-x_start)/3,x_end];
    else
        Xi_to_plot = [x_start,x_start+(knot_vec_int(3)-knot_vec_int(2))/3,knot_vec_int(2:end-1),x_end-((knot_vec_int(3)-knot_vec_int(2))/3),x_end];
    end
elseif p == 4
    if n_elem == 1
        Xi_to_plot = [x_start,x_start+(x_end-x_start)/4,x_start+(x_end-x_start)/2,x_end - (x_end-x_start)/4,x_end];
    else
        if mod(n_elem,2) == 0
            knot_spacing = (knot_vec_int(3)-knot_vec_int(2));

            Xi_to_plot = [x_start,x_start+knot_spacing/4,x_start+knot_spacing*3/4, knot_vec_int(2:floor(end/2))+knot_spacing/2,...
                knot_vec_int(floor(end/2)+1:end-2)+knot_spacing/2,x_end-knot_spacing*3/4,x_end-knot_spacing/4,x_end];
            
        else
            knot_spacing = (knot_vec_int(3)-knot_vec_int(2));
            Xi_to_plot = [x_start,x_start+knot_spacing/4,x_start+knot_spacing*3/4, knot_vec_int(3:end/2)-knot_spacing/2,...
                x_start+(x_end-x_start)/2,knot_vec_int(end/2+1:end-2)+knot_spacing/2,x_end-knot_spacing*3/4,x_end-knot_spacing/4,x_end];
        end
    end
end

Xi_to_plot
end

