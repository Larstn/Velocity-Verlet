%% Return Indicies
% returns the indicies of the x-v-vector corresponding to the positions and
% velocities, as well as the indicies for the acceleration in all three
% dimensions for the acceleartion vector
%%
function [id_x1, id_x2, id_x3, id_v1, id_v2, id_v3, id_a1, id_a2, id_a3]= get_Index3D(nParticle)

id_x1 = 1:nParticle;
id_x2 = id_x1 + nParticle;
id_x3 = id_x2 + nParticle;
id_v1 = id_x3 + nParticle;
id_v2 = id_v1 + nParticle;
id_v3 = id_v2 + nParticle;

id_a1 = id_x1;
id_a2 = id_x2;
id_a3 = id_x3;
end