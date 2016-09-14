function newRadius=getNewRadius(oldRadius)
% GETNEWRADIUS   Ask user for new dilator radius in range 0-10 and print
% used values to the command window
    
    need_new_radius=true;
    while need_new_radius
        user_radius=input(['Old radius was ' num2str(oldRadius) ...
                            '\nEnter new refinement radius: '],'s');
        user_radius = str2double(user_radius);
        if (user_radius>0 && user_radius < 10)
            need_new_radius=false;
            disp(['Using new refinement radius = ' num2str(user_radius)]);
            newRadius=user_radius;
            clear user_radius;
        else
            disp('Invalid response. Radius must be an integer, between 0 and 10')
        end
    end
    clear need_new_radius
end
