function iou = iouGP_full(X,quat,min_z,max_z,X_Ref,pts_ref,basisAngles,nIOUS)
    % suppress warnings
    warning('off')

    % translate and rotate reference points
    pos_ref = X_Ref(1:3); or_ref = X_Ref(5);
    R = [cos(or_ref) -sin(or_ref) 0; sin(or_ref) cos(or_ref) 0; 0 0 1];
    pts_ref = R*pts_ref + pos_ref;

    % get min max of shape union
    minh_union = min(min_z,X_Ref(3)-X_Ref(7)/2);
    maxh_union = max(max_z,X_Ref(3)+X_Ref(7)/2);

    % calculate multiple IoUs
    ious = zeros(1,nIOUS);
    hs = linspace(minh_union,maxh_union,nIOUS);
    for i = 1:nIOUS
        % get shape points of estimate
        if hs(i) < min_z || hs(i) > max_z
            pts = [0;0];
        else
            pts = getShapePoints2dTop(X,quat,basisAngles,hs(i));
        end

        % get shape points of reference
        if hs(i) < X_Ref(3)-X_Ref(7)/2 || hs(i) > X_Ref(3)+X_Ref(7)/2
            pts_r = [0;0];
        else
            idx = pts_ref(3,:) < hs(i)+0.1 & pts_ref(3,:) > hs(i)-0.1;
            pts_r = pts_ref(:,idx);
            k = convhull(pts_r(1:2,:)');
            pts_r = pts_r(:,k);
        end
    
        % calculate intersection over union
        ious(i) = calculateIoU(pts, pts_r);
    end

    % mean value
    iou = mean(ious,2);
end

function pts = getShapePoints2dTop(X,quat,basisAngles,h)
    % get pose and shape parameters
    pos = X(1:3); shape = X(13:end);

    % calculate the rotation matrix from Local to Global by the quternions
    R = rotation_matrix_from_global_to_local(quat); R = R';
    % calculate points in cartesian coordinates
    [x,y,z] = sph2cart(basisAngles(:,1),basisAngles(:,2),shape);
    pts = [x';y';z'];
    % transform the extent to the Global frame
    pts = pos + R*pts;

    % calculate height points
    idx = pts(3,:) > h-0.1 & pts(3,:) < h+0.1;
    pts = pts(1:2,idx);
    try
        idx_pts = convhull(pts');
        pts = pts(:,idx_pts);
    end
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