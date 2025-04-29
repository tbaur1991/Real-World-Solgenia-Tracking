function [iou,h,min_z,max_z] = iouSH(X,X_Ref,pts_ref,nf)
    % get shape points for estimate
    [pts,min_z,max_z] = getShapePoints2dTop(X,nf);
    h = max_z - min_z;

    % translate and rotate reference points
    pos_ref = X_Ref(1:2); or_ref = X_Ref(5);
    R = [cos(or_ref) -sin(or_ref); sin(or_ref ) cos(or_ref )];
    pts_ref = R*pts_ref + pos_ref;

    % calculate intersection over union using alpha shapes and convex hull
    iou = calculateIoU(pts, pts_ref);
end

function [pts,min_z,max_z] = getShapePoints2dTop(X,nf)
    % get pose and shape parameters
    pos = X(1:2); pos_z = X(3); or = X(5); shape = X(7:end);

    % calculate shape points
    phis = 0:2*pi/50:2*pi; thetas = 0:pi/25:pi;
    rs = zeros(length(thetas),length(phis));
    for i = 1:length(thetas)
        for j = 1:length(phis)
            rs(i,j) = sph_harmonics(shape,thetas(i),phis(j),nf,length(shape));
        end
    end

    % calculate height
    zs = rs.*cos(thetas)' + pos_z;
    min_z = min(zs(:));
    max_z = max(zs(:));
    
    % calculate contour
    xs = rs.*sin(thetas)'.*cos(phis); ys = rs.*sin(thetas)'.*sin(phis);
    ns = sqrt(xs.^2 + ys.^2); nsm = max(ns);
    pts = [nsm.*cos(phis); nsm.*sin(phis)];
    
    % translate and rotate points
    R = [cos(or) -sin(or); sin(or) cos(or)];
    pts = R*pts + pos;
end

function iou = calculateIoU(pts,pts_ref)
    % polyshapes
    poly = polyshape(pts(1,:),pts(2,:));
    poly_ref = polyshape(pts_ref(1,:),pts_ref(2,:));

    % intersection area
    poly_int = intersect(poly,poly_ref);
    area_int = area(poly_int);

    % union area
    poly_union = union(poly,poly_ref);
    area_union = area(poly_union);

    % intersection over union
    iou = area_int/area_union;
    if isinf(iou) || isnan(iou)
        iou = 0;
    end
end