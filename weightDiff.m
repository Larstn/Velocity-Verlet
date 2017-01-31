%% Differentiation of the weights
% Returns the differentiation of the weights with respect to the two
% relevant positions
%%
function [dw1dx1, dw1dx2, dw2dx1, dw2dx2] = weightDiff(relevant_pos, x_p)

    dist = relevant_pos(2) - relevant_pos(1);

    dw1dx1 = (relevant_pos(2) - x_p) / (dist^2);
    dw1dx2 = (x_p -relevant_pos(1)) / (dist^2);

    dw2dx1 = (x_p - relevant_pos(2)) / (dist^2);
    dw2dx2 = (relevant_pos(1) - x_p) / (dist^2);

end