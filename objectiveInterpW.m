%% Weights of objective Function 
% Returns the weights for the interpolation of the position between the two
% relevent indicies
%%
function [wobj_1, wobj_2] = objectiveInterpW(relevant_pos, x_p)


pos_2 = relevant_pos(2) - x_p;
pos_1 = x_p - relevant_pos(1);

dist = relevant_pos(2) - relevant_pos(1);

wobj_1 = 1 - pos_1/dist;
wobj_2 = 1 - pos_2/dist;


end