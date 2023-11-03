function G = GLV(images)

[width,height,img_no]=size(images);
G = ones(width,height,img_no);



   for k=1:img_no
     for i=2:width-1
        for j=2:height-1
       
            Mean = (images(i-1,j-1,k) + images(i,j-1,k) + images(i+1,j-1,k) +...
                        images(i-1,j,k) + images(i,j,k) + images(i+1,j,k) +...
                        images(i-1,j+1,k) + images(i,j+1,k) + images(i+1,j+1,k))/9;
            G(i,j,k) = ((images(i-1,j-1,k)-Mean)^2 + (images(i,j-1,k)-Mean)^2 + (images(i+1,j-1,k)-Mean)^2 +...
                        (images(i-1,j,k)-Mean)^2 + (images(i,j,k)-Mean)^2 + (images(i+1,j,k)-Mean)^2 +...
                        (images(i-1,j+1,k)-Mean)^2 + (images(i,j+1,k)-Mean)^2 + (images(i+1,j+1,k)-Mean)^2)/(9-1);
        end
     end
   end
G=smooth3(G,'box',5);
end