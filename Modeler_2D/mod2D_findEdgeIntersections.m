function [t_a,t_b,La,Lb] = mod2D_findEdgeIntersections(a0,a1,b0,b1)

da = permute(a1 - a0,[3 4 1 2]);
db = permute(b1 - b0),[3 4 2 1]);

La = permute(sqrt(sum(da.^2,2)),[3 4 1 2]);
Lb = permute(sqrt(sum(db.^2,2)),[3 4 2 1]);

% (x_a,y_a) = a0 + t_a*(da/La);
% (x_b,y_b) = b0 + t_b*(db/Lb);

% t_a*x_da/La - t_b*x_db/Lb = -(x_a0 - x_b0)
% t_a*y_da/La - t_b*y_db/Lb = -(y_a0 - y_b0)
  
%     | x_da/La -x_db |
% A = |               |
%     | y_da/La -y_db |

%     | -(x_a0 - x_b0) |
% b = |                |
%     | -(y_a0 - y_b0) |

% dim map:
% #1 - A matrix column index (i)
% #2 - A matrix row index (j)
% #3 - a vector index
% #4 - b vector index
A = [repmat([da(1,1,:,1)./La da(1,1,:,2)./La],[1 1 1 size(db,4)]) ; ... 
     repmat([db(1,1,1,:)./Lb db(1,1,2,:)./Lb],[1 1 size(da,3) 1])];

% Solution vector     
b = -bsxfun(@minus,a0,b0)

% The solution is given with the inversion of A
det_A = A(1,1,:,:).*A(2,2,:,:) - A(1,2,:,:).*A(2,1,:,:);
invA = bsxfun(@times,[[A(2,2,:,:) -A(1,2,:,:)] ; [-A(2,1,:,:) A(1,1,:,:)]],1./det_A);

t = [(invA(1,1,:,:).*b(1,1,:,:) + invA(1,2,:,:).*b(2,1,:,:)) ; ...
     (invA(2,1,:,:).*b(1,1,:,:) + invA(2,2,:,:).*b(2,1,:,:))];
     
% Remap dimension such that:
% #1 - a vector index
% #2 - b vector index
t_a = permute(t(1,1,:,:),[3 4 1 2]);
t_b = permute(t(2,1,:,:),[3 4 1 2]);

% Hence, transpose the Lb vector.
Lb = Lb.';

end