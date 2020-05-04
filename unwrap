% ======================= Unwrap Phase In Matlab
function output = unwrap(imput)
 NUM = length(imput);
 Add = zeros(NUM,1);
 for index = 1:NUM-1
     if imput(index+1,1)*imput(index,1)<0
         if pi-abs(imput(index))-abs(imput(index+1))<0
             if imput(1,1)<0
                 Add(index+1:NUM,1) = -2*pi*ones(NUM-index,1);
             else
                 Add(index+1:NUM,1) = 2*pi*ones(NUM-index,1);
             end
         end
     end
 end
 output = imput + Add;
%  figure(1)
%  plot([1:NUM],output,'Linewidth',1.5);
%  figure(2)
%  plot([1:NUM],imput,'Linewidth',1.5);
