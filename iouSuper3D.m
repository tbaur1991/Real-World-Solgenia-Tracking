function iou = iouSuper3D(X,X_Ref,pts_ref)
    % get shape points for estimate
    pts = getShapePoints2dTop(X,true);

    % translate and rotate reference points
    pos_ref = X_Ref(1:2); or_ref = X_Ref(5);
    R = [cos(or_ref) -sin(or_ref); sin(or_ref ) cos(or_ref )];
    pts_ref = R*pts_ref + pos_ref;

    % calculate intersection over union using alpha shapes and convex hull
    iou = calculateIoU(pts, pts_ref, X, X_Ref);
end

function pts = getShapePoints2dTop(X,taper)
    % get pose and shape parameters
    pos = X(1:2); or = X(5); a = c1(X(7),0,'lower'); b = c1(X(8),0,'lower');
    e = c1(X(10),1,'lower');

    % calculate shape points
    theta = 0:0.1:2*pi;
    pts(1,:) = a*sign(cos(theta)).*abs(cos(theta)).^(2/e);
    pts(2,:) = b*sign(sin(theta)).*abs(sin(theta)).^(2/e);
    if taper
        ty = 2/pi*atan(X(11));
        pts(2,:) = (ty*pts(1,:)/a + 1).*pts(2,:);
    end

    % translate and rotate points
    R = [cos(or) -sin(or); sin(or) cos(or)];
    pts = R*pts + pos;
end

function iou = calculateIoU(pts, pts_ref, X, X_Ref)
    % get state variables
    z = X(3); z_Ref = X_Ref(3);
    h = c1(X(9),0,'lower'); h_Ref = X_Ref(7);

    % polyshapes in top view
    poly = polyshape(pts(1,:),pts(2,:));
    poly_ref = polyshape(pts_ref(1,:),pts_ref(2,:));

    % intersection area
    poly_int = intersect(poly,poly_ref);
    area_int = area(poly_int);

    % intersection height and volume
    i1 = fixed.Interval(z-h/2,z+h/2);
    i2 = fixed.Interval(z_Ref-h_Ref/2,z_Ref+h_Ref/2);
    h_int = intersect(i1,i2);
    if ~isempty(h_int)
        vol_int = area_int*(h_int.RightEnd - h_int.LeftEnd);
    else
        vol_int = 0;
    end

    % union height and volume
    shape = [pts pts; (z-h/2)*ones(1,size(pts,2)) (z+h/2)*ones(1,size(pts,2))];
    [~, vol_un] = convhull([shape'; [pts_ref pts_ref; (z_Ref-h_Ref/2)*ones(1,size(pts_ref,2)) (z_Ref+h_Ref/2)*ones(1,size(pts_ref,2))]']);

    % intersection over union
    iou = double(vol_int/vol_un);
    if isinf(iou) || isnan(iou)
        iou = 0;
    end
end