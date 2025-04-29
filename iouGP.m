function [iou,h,min_z,max_z] = iouGP(X,quat,X_Ref,pts_ref,basisAngles)
    % get shape points for estimate
    [pts,min_z,max_z] = getShapePoints2dTop(X,quat,basisAngles);
    h = max_z - min_z;

    % translate and rotate reference points
    pos_ref = X_Ref(1:2); or_ref = X_Ref(5);
    R = [cos(or_ref) -sin(or_ref); sin(or_ref ) cos(or_ref )];
    pts_ref = R*pts_ref + pos_ref;

    % calculate intersection over union using alpha shapes and convex hull
    iou = calculateIoU(pts, pts_ref);
end

function [pts,min_z,max_z] = getShapePoints2dTop(X,quat,basisAngles)
    % get pose and shape parameters
    pos = X(1:3); shape = X(13:end);

    % calculate the rotation matrix from Local to Global by the quternions
    R = rotation_matrix_from_global_to_local(quat); R = R';
    % calculate points in cartesian coordinates
    [x,y,z] = sph2cart(basisAngles(:,1),basisAngles(:,2),shape);
    pts = [x';y';z'];
    % transform the extent to the Global frame
    pts = pos + R*pts;
    
    % calculate height
    min_z = min(pts(3,:));
    max_z = max(pts(3,:));

    % 2D convex hull
    k = convhull(pts(1:2,:)');
    pts = pts(1:2,k);
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