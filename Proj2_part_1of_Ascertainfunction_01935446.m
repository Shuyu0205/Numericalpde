%take a matrix phi and then plot dphi/dx and dphi/dy at y=0

   layer=phi(1,:);
   layer_2=phi(2,:);
   layer_u=[];
   layer_v=[];
   for i=2:1:length(layer)-1
          layer_u=[layer_u,(layer(i+1)-layer(i-1))/(2*h)];
   end
   for i=2:1:length(layer)-1
          layer_v=[layer_v,(layer_2(i)-layer(i))/(h)];
   end
