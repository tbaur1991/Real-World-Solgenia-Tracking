function iou = iouFCDS(X,X_Ref,pts_ref,nu,nth)
    % get shape points for estimate
    pts = getShapePoints2dTop(X,nu,nth);

    % translate and rotate reference points
    pos_ref = X_Ref(1:2); or_ref = X_Ref(5);
    R = [cos(or_ref) -sin(or_ref); sin(or_ref ) cos(or_ref )];
    pts_ref = R*pts_ref + pos_ref;

    % calculate intersection over union using alpha shapes and convex hull
    iou = calculateIoU(pts, pts_ref);
end

function pts = getShapePoints2dTop(X,nu,nth)
    % get pose and shape parameters
    pos = X(1:2); or = X(5); shape = X(8:end);

    % calculate shape points
    us = -1:0.1:1; thetas = linspace(0,2*pi,100);
    rs = zeros(length(us),length(thetas));
    for i = 1:length(us)
        for j = 1:length(thetas)
            rs(i,j) = fourier_chebychev_series(shape,thetas(j),us(i),nu,nth);
        end
    end
    rs = max(rs);
    pts = [rs.*cos(thetas); rs.*sin(thetas)];

    % translate and rotate points
    R = [cos(or) -sin(or); sin(or) cos(or)];
    pts = R*pts + pos;
end

function iou = calculateIoU(pts, pts_ref)
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