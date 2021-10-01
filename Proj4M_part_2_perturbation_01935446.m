%take a matrix phi and then plot dphi/dx and dphi/dy at y=0
function u=Proj4_part_2_perturbation_01935446(phi,rho,h,q,s)
   layer=phi(1,:);
   layer_2=phi(2,:);
   layer_u=[];
   layer_v=[];
   for i=2:1:length(layer)-1
          layer_u=[layer_u,(layer(i+1)-layer(i-1))/(2*h)];
   end
   layer_u=[0,layer_u,0];
   for i=2:1:length(layer)-1
          layer_v=[layer_v,(layer_2(i)-layer(i))/(h)];
   end
   layer_v=[0,layer_v,0];
   plot(-q:h:s,layer_u);
end