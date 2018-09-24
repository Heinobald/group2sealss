function distance = geocalc(lat1, lon1, lat2, lon2)
    % calculates the distance between a pair of coordinates, 
    % used to calculate cumulative distance of the tours
    lat1 = deg2rad(double(lat1));
    lon1 = deg2rad(double(lon1));
    lat2 = deg2rad(double(lat2));
    lon2 = deg2rad(double(lon2));

    dlon = lon1 - lon2;
    EARTH_R = 6372.8;

    y = sqrt((cos(lat2)*sin(dlon))^2 + (cos(lat1) * sin(lat2) - sin(lat1) * cos(lat2) * cos(dlon))^2);
    x = sin(lat1) * sin(lat2) + cos(lat1) * cos(lat2) * cos(dlon);
    c = atan2(y, x);
    distance = EARTH_R * c;
    %return distance