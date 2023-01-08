function isin = isinside2Dset(V, P)
%% isinside2Dset : function to check if a vertex is located inside or outside a given
% 2D set, boundary not included (opened set).
%
% Author & support : nicolas.douillet (at) free.fr, 2023.
%
%
% Syntax
%
% isin = isinside2Dset(V, P);
%
%
% Description
%
% isin = isinside2Dset(V, P) computes the boolean isin which is true/logical 1
% in the case the vertex P belongs to the opened convex set V. isin is false/logical
% 0 in the case vertex P belongs to the complementary set or the boundary.
%
%
% Input arguments
%
%       [ |  |  ]
% - V = [ Vx Vy ], real matrix double, the convex set, with size(V,2) = 2.
%       [ |  |  ]
%
%       [ |  | ]
% - P = [Px Py ], real row vector or matrix double, the coordinates of the vertex / vertices  to check. Size(P,2) = size(V,2).
%       [ |  | ]
%
%
% Output argument
%
%          [      |      ]
% - isin = [logical 1 / 0], logical true (1)/false (0) scalar / column vector. The boolean result. Size(isin) = [size(P,1),1].
%          [      |      ]
%
%
% Example #1 : 2D random point cloud
% N = 16;
% V = 2*(rand(N,2)-0.5);
% G = mean(V,1);
% V = V - G;
% theta = atan2(V(:,2),V(:,1));
% [~,i] = sort(theta);
% V = V(i,:);
% [A B] = meshgrid(-1:0.1:1);
% P = cat(2,A(:),B(:));
% 
% figure
% line([V(:,1); V(1,1)],[V(:,2); V(1,2)],'Color',[0 0 1],'LineWidth',2), hold on;
% set(gcf,'Color',[0 0 0]), set(gca,'Color',[0 0 0],'XColor',[1 1 1],'YColor',[1 1 1],'FontSize',16);
% xlabel('X'), ylabel('Y');
% 
% isin = cell2mat(cellfun(@(c) isinside2Dset(V,c),num2cell(P,2),'un',0));
% ColorSpec = cell2mat(cellfun(@(c) cat(2,~c,c,0),num2cell(isin,2),'un',0));
% 
% cellfun(@(r1,r2) plot(r1(1,1),r1(1,2),'+','Color',r2,'MarkerSize',4,'LineWidth',4),num2cell(P,2),num2cell(ColorSpec,2),'un',0);
% axis equal, axis tight;


%% Input parsing
assert(nargin > 1,'Not enought input arguments.');
assert(nargin < 3,'Too many input arguments.');
assert(isequal(size(V,2),size(P,2),2),'All the inputs must have the same number of colums (two dimensions here).');

%% Body
nb_test_vtx = size(P,1);
G = mean(V,1);
nb_vtx = size(V,1);

V1 = cat(2,V,zeros(nb_vtx,1));
V2 = circshift(V1,1,1);

% Hyperplane normals
ui = V1 - V2;
ni = cat(2,-ui(:,2),ui(:,1));
ni = ni ./ sqrt(sum(ni.^2,2));
ni = ni(:,1:2);

% Coherently orient outward hyperplane normals
orientation = sign(dot(ni,V-repmat(G,[nb_vtx,1]),2));

if isequal(orientation,-ones(nb_vtx,1)) || isequal(orientation,ones(nb_vtx,1))    
    ni = ni.*orientation;    
end

ni = cat(2,ni,zeros(size(ni,1),1));

% Vertex normals
mi = 0.5*(ni+circshift(ni,-1,1)); 
mi = mi ./ sqrt(sum(mi.^2,2));
                         
% (V,u) lines corresponding to [Vi; Vi+1] segments
ui = V1 - V2; 

% Hi, projections of P on each segment support line
P = cat(2,P,zeros(nb_test_vtx,1));

% TODO : use arrayfun instead
[d2Hi,Hi] = cellfun(@(r1,r2) point_to_line_distance(P,r1,r2),num2cell(ui,2),num2cell(V1,2),'un',0);    

d2Hi = cell2mat(d2Hi);
Hi = cell2mat(Hi);

% Check if H_i belongs to [Vi; Vi+1 segment]
dot_prod_sign = sign(dot(V1-Hi,V2-Hi,2));

dst_mat = sqrt(sum((P(:,1:2)-V).^2,2));
[min_dst2V,nrst_vtx_idx] = min(dst_mat);
Mi = mi(nrst_vtx_idx,:);
Vi = V(nrst_vtx_idx,:);


f = find(1-dot_prod_sign);

if ~isempty(f)
    
    Hi = Hi(f,:);
    di = d2Hi(f,:);
    Ni = ni(f,:);
    [min_dst2H,i] = min(di);
    H = Hi(i,:);
    N = Ni(i,:);
    
    if min_dst2H <= min_dst2V
        
        % 'hyperplane'
        isin = sign(1+sign(dot(H-P,N)));
        
    else % if min_dst2H > min_dst2V
        
        % 'vertex'
        isin = sign(1+sign(dot(Vi-P(:,1:2),Mi(:,1:2))));
        
    end
    
else % if isempty(f) % cas où f vide (le point ne se projette sur aucun segment);
    
    isin = 0; % out by default
    
    if sign(dot(Vi-P(:,1:2),Mi(:,1:2))) > 0 % the closest vertex only
        
        isin = 1;
        
    end
    
end


end % isinside2Dset


%% point_to_line_distance subfunction
function [d2H, H] = point_to_line_distance(P, u, I0)
%
% Author & support : nicolas.douillet (at) free.fr, 2019-2023.


dimension = size(P,2);  
nb_pts    = size(P,1);

% Body
t_H = (u(1)*(P(:,1)-repmat(I0(1),[nb_pts,1])) + ...
       u(2)*(P(:,2)-repmat(I0(2),[nb_pts,1])) + ...
       u(3)*(P(:,3)-repmat(I0(3),[nb_pts,1])) ) / ...
       sum(u.^2); 

x_H = I0(1) + t_H*u(1);
y_H = I0(2) + t_H*u(2);
z_H = I0(3) + t_H*u(3);


% Orthogonal projected point
H = zeros(size(P,1),dimension);
H(:,1) = x_H;
H(:,2) = y_H;
H(:,3) = z_H;


% Distance
d2H = sqrt(sum((P-H).^2,2));
H = H(:,1:dimension);


end % point_to_line_distance